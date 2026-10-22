"""Exercise the browser pager with the embedded renderer and live widget messages."""

import itertools
import re
import sys
import unicodedata
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import patch

from symbolica import E, MathMLViewer, PrintMode, S
from symbolica.community.tensor import (
    DisplaySettings,
    Representation,
    TensorExpression,
    TensorName,
)


@unittest.skipIf(sys.platform in ("emscripten", "wasi"), "native compiler mock")
class LeafCompilationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        bootstrap = TensorName("leaf_compilation::bootstrap")(
            Representation.mink(4)("mu")
        ).paged()
        try:
            bootstrap._page()
        finally:
            bootstrap.close()
        import _spenso_paging

        cls.paging = _spenso_paging

    def setUp(self):
        self.document = ""
        self.payload = {
            "tree": b"",
            "notation": "",
            "settings": "(:)",
            "power_base": False,
        }

        def compile_native(root, timeout):
            (root / "out.html").write_text(self.document)

        self.compiler = self.enterContext(
            patch.object(self.paging, "_compile_native", side_effect=compile_native)
        )

    def test_indented_leaf_matches_compact_math_with_phantoms_and_alignment(self):
        compact = (
            "<mrow><msubsup><mphantom><mi>A</mi></mphantom>"
            '<mrow><maligngroup groupalign="left"/><mi>μ</mi>'
            '<malignmark edge="left"/></mrow><mspace width="0.5em"/>'
            "</msubsup></mrow>"
        )
        expected = ET.tostring(ET.fromstring(compact), encoding="unicode")
        for indented in (False, True):
            with self.subTest(indented=indented):
                root = ET.fromstring(f"<math>{compact}</math>")
                if indented:
                    ET.indent(root)
                self.document = ET.tostring(root, encoding="unicode")
                self.assertEqual(self.paging.compile_leaf(self.payload), expected)

    def test_token_text_and_non_xml_spaces_are_preserved(self):
        for namespace in ("", ' xmlns="http://www.w3.org/1998/Math/MathML"'):
            for text in (" \t\n ", "  a b \t ", "\u00a0", "\u2003"):
                with self.subTest(namespace=namespace, text=text):
                    self.document = (
                        f"<math{namespace}>\n <mrow>\n"
                        + "".join(
                            f"<{tag}>{text}<malignmark/>{text}</{tag}>\n"
                            for tag in ("mi", "mo", "mn", "mtext", "ms")
                        )
                        + '<annotation encoding="text/plain"> \t </annotation>\n'
                        + "<mrow>\u00a0<mspace width='1em'/>\u2003</mrow>\n"
                        + "</mrow>\n</math>"
                    )
                    fragment = self.paging.compile_leaf(self.payload)
                    row = ET.fromstring(fragment)
                    self.assertIsNone(row.text)
                    self.assertIsNone(row.tail)
                    for token in row[:5]:
                        self.assertEqual(token.text, text)
                        self.assertEqual(token[0].tail, text)
                        self.assertIsNone(token.tail)
                    self.assertEqual(row[5].text, " \t ")
                    self.assertEqual(row[6].text, "\u00a0")
                    self.assertEqual(row[6][0].tail, "\u2003")

    def test_xml_space_inheritance_resets_and_parent_owned_tails(self):
        self.document = (
            "<math>\n<mrow>\n"
            '<mrow xml:space="preserve">  <mrow> \t<mi>x</mi> \n</mrow> \t'
            '<mrow xml:space="default"> \n<mi>y</mi> \t</mrow>  </mrow>\n'
            '<mrow xml:space="preserve"> \n<mi>z</mi> \t</mrow>\n'
            "</mrow>\n</math>"
        )
        row = ET.fromstring(self.paging.compile_leaf(self.payload))
        self.assertIsNone(row.text)
        self.assertIsNone(row.tail)
        preserved, sibling = row
        self.assertEqual(preserved.text, "  ")
        self.assertIsNone(preserved.tail)
        inherited, reset = preserved
        self.assertEqual(inherited.text, " \t")
        self.assertEqual(inherited[0].tail, " \n")
        self.assertEqual(inherited.tail, " \t")
        self.assertIsNone(reset.text)
        self.assertIsNone(reset[0].tail)
        self.assertEqual(reset.tail, "  ")
        self.assertEqual(sibling.text, " \n")
        self.assertEqual(sibling[0].tail, " \t")
        self.assertIsNone(sibling.tail)
        space = "{http://www.w3.org/XML/1998/namespace}space"
        self.assertEqual(preserved.get(space), "preserve")
        self.assertEqual(reset.get(space), "default")
        self.document = (
            '<math xml:space="preserve"><mrow> \n<mi>x</mi> \t</mrow> \n'
            '<mrow xml:space="default"> \t<mi>y</mi> \n</mrow> \t</math>'
        )
        root = ET.fromstring(f"<math>{self.paging.compile_leaf(self.payload)}</math>")
        self.assertEqual(root[0].text, " \n")
        self.assertEqual(root[0][0].tail, " \t")
        self.assertEqual(root[0].tail, " \n")
        self.assertIsNone(root[1].text)
        self.assertIsNone(root[1][0].tail)
        self.assertEqual(root[1].tail, " \t")

    def test_annotation_xml_content_is_preserved(self):
        content = (
            '<annotation-xml encoding="application/xhtml+xml"> \t'
            '<div xmlns="http://www.w3.org/1999/xhtml"> \n<span/> \t</div> \n'
            "</annotation-xml>"
        )
        self.document = f"<math>\n<mrow>\n<mi>x</mi>\n{content}\n</mrow>\n</math>"
        row = ET.fromstring(self.paging.compile_leaf(self.payload))
        self.assertIsNone(row.text)
        self.assertIsNone(row[0].tail)
        self.assertIsNone(row[1].tail)
        self.assertEqual(
            ET.tostring(row[1], encoding="unicode"),
            ET.tostring(ET.fromstring(content), encoding="unicode"),
        )

    def test_fragment_budget_is_applied_after_compaction_in_utf8_bytes(self):
        for size in (4095, 4096, 4097):
            with self.subTest(bytes=size):
                # μ occupies two bytes; significant token text must survive.
                text = "μ" * 1000 + "x" * (size - 2028)
                compact = f"<mrow><mtext>{text}</mtext></mrow>"
                self.assertEqual(len(compact.encode()), size)
                self.document = f"<math>\n{' ' * 5000}<mrow>\n  <mtext>{text}</mtext>\n</mrow>\n</math>"
                fragment = self.paging.compile_leaf(self.payload)
                self.assertEqual(fragment, compact if size <= 4096 else None)

    def test_compaction_keeps_path_removal_slash_styles_and_errors(self):
        self.document = (
            '<math data-symbolica-path="0">\n<mrow data-symbolica-path="1">\n'
            '<mi data-symbolica-path="2" data-spenso-cancel="updiagonal" '
            'style="color:red;">p</mi>\n</mrow>\n</math>'
        )
        row = ET.fromstring(self.paging.compile_leaf(self.payload))
        self.assertTrue(all("data-symbolica-path" not in n.attrib for n in row.iter()))
        self.assertEqual(row[0].get("style"), f"color:red;{self.paging.CANCEL_STYLE}")
        self.assertEqual(row[0].get("data-spenso-cancel"), "updiagonal")
        self.document = "<html/>"
        with self.assertRaisesRegex(RuntimeError, "Renderer returned no MathML"):
            self.paging.compile_leaf(self.payload)
        self.document = "<math><mi></math>"
        with self.assertRaises(ET.ParseError):
            self.paging.compile_leaf(self.payload)
        self.compiler.side_effect = RuntimeError("bad notation")
        with self.assertRaisesRegex(RuntimeError, "bad notation"):
            self.paging.compile_leaf(self.payload)


class BrowserPagingTests(unittest.TestCase):
    def setUp(self):
        self.value = TensorExpression(
            E("+".join(f"{i + 1}*wasm_paging::x^{i}" for i in range(301)))
        )
        self.original = self.value.to_expression()
        self.pager = self.value.paged()
        self.addCleanup(self.pager.close)
        self.assertIs(type(self.pager), MathMLViewer)
        # Also exercise this path in native CI. Real Pyodide needs no emulation.
        self.enterContext(patch.object(sys, "platform", "emscripten"))
        self.enterContext(patch.dict(sys.modules, {"typst": None}))
        self.enterContext(
            patch("subprocess.run", side_effect=AssertionError("WASM has no processes"))
        )

    def assertPage(self, page):
        self.assertNotIn("Page rendering failed", page["html"])
        self.assertLessEqual(len(page["html"].encode()), 256 * 1024)
        math = re.findall(r"<math\b.*?</math>", page["html"], re.DOTALL)
        self.assertTrue(math)
        roots = [ET.fromstring(expression) for expression in math]
        self.assertLessEqual(sum(sum(1 for _ in root.iter()) for root in roots), 10_000)
        marker_ids = {
            int(node.attrib["data-symbolica-marker"])
            for root in roots
            for node in root.iter()
            if "data-symbolica-marker" in node.attrib
        }
        self.assertEqual(marker_ids, set(page["marker_ids"]))
        for root in roots:
            root.attrib.pop("data-spenso-math", None)
        return roots

    def test_embedded_pages_are_contiguous_and_reversible(self):
        pages = []
        while True:
            page = self.pager._page()
            self.assertPage(page)
            pages.append((page["start"], page["end"], page["html"]))
            if page["next"] is None:
                break
            self.pager._action({"action": "next"})
        self.assertEqual(pages[0][0], 0)
        self.assertEqual(pages[-1][1], 301)
        self.assertTrue(all(a[1] == b[0] for a, b in itertools.pairwise(pages)))
        for start, end, html in reversed(pages[:-1]):
            page = self.pager._action({"action": "previous"})
            self.assertEqual(
                (page["start"], page["end"], page["html"]), (start, end, html)
            )
        self.assertEqual(self.value.to_expression(), self.original)

    def test_live_widget_navigation_and_page_size(self):
        widget = self.pager._get_widget()
        self.pager._on_message(widget, {"action": "attach", "view": "browser"}, [])
        self.assertTrue(widget.page["supports_expansion"])
        self.assertPage(widget.page)
        for action, expected in (("next", 25), ("previous", 0)):
            self.pager._on_message(widget, {"action": action, "request": action}, [])
            self.assertEqual(widget.page["start"], expected)
            self.assertPage(widget.page)
        self.pager._on_message(widget, {"action": "size", "value": 100}, [])
        self.assertEqual(widget.page["page_size"], 100)
        self.assertEqual(widget.page["end"], 100)
        self.assertPage(widget.page)
        self.pager._on_message(widget, {"action": "dispose", "view": "browser"}, [])
        self.assertIsNone(self.pager._widget)

    def test_nested_sums_keep_subexpression_navigation(self):
        value = TensorExpression(
            (E("x+1") ** 100).expand() * (E("y+1") ** 100).expand()
        )
        original = value.to_expression()
        pager = value.paged()
        self.addCleanup(pager.close)
        parent = pager._page()
        self.assertPage(parent)
        compact = parent["layouts"][-1]
        self.assertEqual(parent["kind"], "Product")
        self.assertEqual(parent["total"], 2)
        target = next(path for path, label in compact["markers"] if "sum" in label)
        self.assertEqual(len(target), 1)
        child = pager._action(
            {"action": "open", "value": target, "layout": compact["layout"]}
        )
        self.assertPage(child)
        self.assertTrue(child["back"])
        self.assertPage(pager._action({"action": "next"}))
        self.assertEqual(pager._action({"action": "back"})["html"], parent["html"])
        self.assertEqual(value.to_expression(), original)

    def test_markers_share_hover_inspection_and_navigation(self):
        page = self.pager._page()
        target = page["markers"][-1][0]
        field, preview = self.pager._message({"action": "preview", "value": target})
        self.assertEqual(field, "preview")
        self.assertPage(preview)
        self.assertEqual(preview["start"], page["end"])
        self.assertTrue(preview["can_expand"])
        field, inspected = self.pager._message(
            {"action": "preview", "value": target, "popup": True, "size": 10}
        )
        self.assertEqual(field, "preview")
        self.assertPage(inspected)
        self.assertEqual((inspected["start"], inspected["end"]), (0, 10))
        second = self.pager._action({"action": "next"})
        self.assertPage(second)
        self.assertEqual(
            len({tuple(p) for p, _ in second["markers"]}), len(second["markers"])
        )
        for target, start in zip(second["markers"], second["marker_starts"]):
            _, preview = self.pager._message({"action": "preview", "value": target[0]})
            self.assertEqual(preview["start"], start)

    def test_small_tensor_pages_preserve_permuted_ports_and_custom_printers(self):
        rep = Representation.euc(3)
        value = TensorName("paging_ports::A", print={"typst": "macron(A)"})(
            S("paging_ports::parameter"), 7, rep("i"), Representation.mink(4)("mu")
        ).permute_axes([1, 0])
        original, axes = value.to_expression(), value.axes
        source = bytes(value.__getstate__())
        structure, name, arguments = (
            value.structure,
            value.name.to_expression(),
            value.arguments,
        )
        pager = value.paged(page_size=7)
        self.addCleanup(pager.close)
        page = pager._page()
        self.assertPage(page)
        self.assertIn("mover", page["html"])
        for html in (value.to_html(), page["html"]):
            root = ET.fromstring(
                re.search(r"<math\b.*?</math>", html, re.DOTALL).group()
            )
            hidden = {
                child
                for phantom in root.iter()
                if phantom.tag.rsplit("}", 1)[-1] == "mphantom"
                for child in phantom.iter()
            }
            indices = [
                unicodedata.normalize("NFKC", node.text or "")
                for node in root.iter()
                if node not in hidden
                and node.tag.rsplit("}", 1)[-1] in ("mi", "mtext")
                and unicodedata.normalize("NFKC", node.text or "") in ("mu", "i")
            ]
            self.assertEqual(indices, ["mu", "i"])
        self.assertEqual(value.axes, axes)
        self.assertEqual(value.to_expression(), original)
        self.assertEqual(bytes(value.__getstate__()), source)
        self.assertEqual(value.structure, structure)
        self.assertEqual(value.name.to_expression(), name)
        self.assertEqual(value.arguments, arguments)

    def test_paged_leaves_retain_rich_settings(self):
        value = TensorName("paging_settings::A")(
            Representation.euc(3)("i"), Representation.mink(4)("mu")
        )
        source, axes, structure = (
            bytes(value.__getstate__()),
            value.axes,
            value.structure,
        )
        ports = value.paged(settings=DisplaySettings(tensor_layout="ports"))
        call = value.paged(
            settings=DisplaySettings(tensor_layout="call", show_dimensions=True)
        )
        self.addCleanup(ports.close)
        self.addCleanup(call.close)
        first = ports._page()
        explicit = call._page()
        first_text = "".join(self.assertPage(first)[0].itertext())
        explicit_text = "".join(self.assertPage(explicit)[0].itertext())
        self.assertNotEqual(first_text, explicit_text)
        self.assertIn("3", explicit_text)
        self.assertIn("4", explicit_text)
        self.assertEqual(ports._page()["html"], first["html"])
        self.assertEqual(bytes(value.__getstate__()), source)
        self.assertEqual(value.axes, axes)
        self.assertEqual(value.structure, structure)

    def test_leaf_cache_distinguishes_product_and_power_contexts(self):
        def printer(expression, *, mode, in_product, in_exp_base, **options):
            if mode != PrintMode.Typst:
                return None
            label = "base" if in_exp_base else "factor" if in_product else "single"
            return f'upright("{label}")'

        head = S(
            "paging_context::F",
            is_scalar=True,
            tags=["spenso::print-callback"],
            print=printer,
        )
        leaf = head(1)
        value = TensorExpression(leaf + S("paging_context::x") * leaf + leaf**2)
        source = bytes(value.__getstate__())
        viewer = value.paged()
        self.addCleanup(viewer.close)
        for layout in (0, 1):
            text = "".join(self.assertPage(viewer._render_page(layout))[0].itertext())
            for label in ("single", "factor", "base"):
                self.assertIn(label, text)
        self.assertEqual(bytes(value.__getstate__()), source)

    def test_paged_slash_carries_its_own_style(self):
        gamma = TensorExpression.dirac_gamma(4)
        p = TensorName.vector("paging_slash::p")(Representation.mink(4))
        value = (gamma("a", "b", "mu") * p("mu")).contract()
        source, axes = bytes(value.__getstate__()), value.axes
        viewer = value.paged()
        self.addCleanup(viewer.close)
        for layout in (0, 1):
            root = self.assertPage(viewer._render_page(layout))[0]
            strokes = [
                node
                for node in root.iter()
                if node.get("data-spenso-cancel") == "updiagonal"
            ]
            self.assertEqual(len(strokes), 1)
            self.assertIn("linear-gradient", strokes[0].get("style", ""))
            self.assertIn("currentColor", strokes[0].get("style", ""))
        self.assertEqual(bytes(value.__getstate__()), source)
        self.assertEqual(value.axes, axes)

    def test_inline_expansion_keeps_tensor_notation_and_source(self):
        rep = Representation.mink(4)
        value = (
            TensorName("wasm_expansion::T")(rep("mu"), rep("nu"))
            * (E("wasm_expansion::x+1") ** 60).expand()
        )
        original, axes = value.to_expression(), value.axes
        viewer = value.paged()
        self.addCleanup(viewer.close)
        parent = viewer._page()
        compact = parent["layouts"][-1]
        target = next(path for path, label in compact["markers"] if "sum" in label)
        self.assertTrue(parent["supports_expansion"])
        self.assertEqual(len(target), 1)
        expanded = viewer._action(
            {
                "action": "expand",
                "value": target,
                "layout": compact["layout"],
                "size": 25,
            }
        )
        self.assertPage(expanded)
        self.assertFalse(expanded["back"])
        self.assertIn("data-symbolica-path", expanded["html"])
        self.assertIn("mu", expanded["html"])
        self.assertIn("nu", expanded["html"])
        expanded_compact = viewer._render_page(compact["layout"])
        self.assertPage(expanded_compact)
        self.assertGreater(
            expanded_compact["html"].count("<"), compact["html"].count("<")
        )
        self.assertEqual(value.to_expression(), original)
        self.assertEqual(value.axes, axes)

    def test_oversized_tensor_output_defers_to_native_argument_layout(self):
        rep = Representation.mink(4)
        bootstrap = TensorName("wasm_paging::bootstrap")(rep("mu")).paged()
        bootstrap._page()
        bootstrap.close()
        import _spenso_paging

        value = TensorName("wasm_paging::oversized")(rep("mu"), rep("nu"))
        original = value.to_expression()
        pager = value.paged()
        self.addCleanup(pager.close)
        with patch.object(
            _spenso_paging, "compile_typst", return_value=b"x" * (256 * 1024 + 1)
        ):
            page = pager._page()
        self.assertPage(page)
        native = original.mathml_viewer().page()
        self.assertEqual(
            [ET.tostring(root) for root in self.assertPage(page)],
            [ET.tostring(root) for root in self.assertPage(native)],
        )
        self.assertEqual((page["kind"], page["unit"]), ("Function", "arguments"))
        self.assertEqual(value.to_expression(), original)

    def test_renderer_errors_propagate_from_the_bounded_preview(self):
        rep = Representation.mink(4)
        bootstrap = TensorName("wasm_errors::bootstrap")(rep("mu")).paged()
        bootstrap._page()
        bootstrap.close()
        import _spenso_paging

        value = TensorName("wasm_errors::T")(rep("mu"))
        pager = value.paged()
        self.addCleanup(pager.close)
        with patch.object(
            _spenso_paging, "compile_typst", side_effect=RuntimeError("bad notation")
        ):
            with self.assertRaisesRegex(RuntimeError, "bad notation"):
                pager._page()
        self.assertEqual(self.value.to_expression(), self.original)

    def test_scalar_arithmetic_uses_identical_native_responsive_layouts(self):
        summands = "+".join(f"layout::a{i}" for i in range(10))
        nested = {
            "product": f"z*({summands})",
            "fraction": f"({summands})/(z+w)",
            "power": f"({summands})^2",
            "function": f"f({summands})",
        }
        for label, expression in nested.items():
            with self.subTest(label=label):
                value = TensorExpression(E(expression))
                original = value.to_expression()
                pager, native = value.paged(), original.mathml_viewer()
                self.addCleanup(pager.close)
                self.addCleanup(native.close)
                for layout in range(8):
                    page = pager._render_page(layout)
                    self.assertPage(page)
                    self.assertEqual(
                        [ET.tostring(root) for root in self.assertPage(page)],
                        [
                            ET.tostring(root)
                            for root in self.assertPage(native._render_page(layout))
                        ],
                    )
                self.assertEqual(value.to_expression(), original)

    def test_embedded_tensor_leaf_keeps_semantic_annotations(self):
        rep = Representation.mink(4)
        bootstrap = TensorName("wasm_annotation::bootstrap")(rep("mu")).paged()
        bootstrap._page()
        bootstrap.close()
        import _spenso_paging

        document = (
            "<html><body><math><semantics><mstyle><mrow>"
            "<mi>T</mi><msub><mi>v</mi><mi>mu</mi></msub>"
            "</mrow></mstyle><annotation encoding='text/plain'>tensor</annotation>"
            + "</semantics></math></body></html>"
        )
        pager = TensorName("wasm_annotation::T")(rep("mu")).paged()
        self.addCleanup(pager.close)
        with patch.object(
            _spenso_paging, "compile_typst", return_value=document.encode()
        ):
            page = pager._page()
        self.assertPage(page)
        root = ET.fromstring(
            re.search(r"<math\b.*?</math>", page["html"], re.DOTALL).group()
        )
        annotations = [
            n for n in root.iter() if n.tag.rsplit("}", 1)[-1] == "annotation"
        ]
        self.assertEqual([node.text for node in annotations], ["tensor"])
        self.assertTrue(any(n.tag.rsplit("}", 1)[-1] == "msub" for n in root.iter()))

    def test_generated_aliases_reserve_free_labels_and_keep_leaf_payloads_local(self):
        rep = Representation.mink(4)
        for name in ("F", "v", "w"):
            TensorName(f"alias_history::{name}")
        free, left, right, mink = S(
            "alias_history::F", "alias_history::v", "alias_history::w", "spenso::mink"
        )
        free_index, other_free_index = S("μ", "ρ")
        expression = sum(
            E(f"alias_history::z^{i}")
            * free(mink(4, free_index), mink(4, other_free_index))
            * left(mink(4, S(f"alias_history::i{i}", tags=["spenso::dummy-index"])))
            * right(mink(4, S(f"alias_history::i{i}", tags=["spenso::dummy-index"])))
            for i in range(100)
        )
        interface = TensorName("alias_history::F")(
            rep(free_index), rep(other_free_index)
        ).structure
        value = TensorExpression(expression, structure=interface)
        original = value.to_expression()
        viewer = value.paged(settings=DisplaySettings(index_style="alphabet"))
        self.addCleanup(viewer.close)
        bootstrap = TensorName("alias_history::bootstrap")(rep("a")).paged()
        bootstrap._page()
        bootstrap.close()
        import _spenso_paging

        with patch.object(
            _spenso_paging, "compile_leaf", wraps=_spenso_paging.compile_leaf
        ) as render:
            page = viewer._page()
            first = page["html"]
            while True:
                roots = self.assertPage(page)
                hidden = {
                    child
                    for root in roots
                    for phantom in root.iter()
                    if phantom.tag.rsplit("}", 1)[-1] == "mphantom"
                    for child in phantom.iter()
                }
                numbered_aliases = {
                    child
                    for root in roots
                    for subscript in root.iter()
                    if subscript.tag.rsplit("}", 1)[-1] == "msub"
                    and len(subscript) == 2
                    and unicodedata.normalize("NFKC", subscript[0].text or "")
                    in ("μ", "ν", "ρ", "σ")
                    and subscript[1].tag.rsplit("}", 1)[-1] == "mn"
                    for child in subscript.iter()
                }
                labels = [
                    unicodedata.normalize("NFKC", node.text or "")
                    for root in roots
                    for node in root.iter()
                    if node not in hidden
                    and node not in numbered_aliases
                    and node.tag.rsplit("}", 1)[-1] in ("mi", "mtext")
                ]
                # Every shown term owns one explicit μ and ρ. Generated labels
                # skip both, including when advancing past an earlier alias ν.
                self.assertEqual(labels.count("μ"), page["end"] - page["start"])
                self.assertEqual(labels.count("ρ"), page["end"] - page["start"])
                if page["start"] == 0:
                    self.assertIn("ν", labels)
                    self.assertIn("σ", labels)
                if page["end"] >= 75 or page["next"] is None:
                    break
                page = viewer._action({"action": "next"})
            self.assertGreaterEqual(page["end"], 75)
            for call in render.call_args_list:
                aliases = [
                    line
                    for line in call.args[0]["settings"].splitlines()
                    if "index-aliases:" in line
                ][-1]
                self.assertLessEqual(aliases.count('"spenso::mink"'), 2)
                self.assertLess(len(aliases.encode()), 1024)
            while page["previous"]:
                page = viewer._action({"action": "previous"})
            self.assertEqual(page["html"], first)
        self.assertEqual(value.to_expression(), original)
        self.assertEqual(value.structure, interface)


if __name__ == "__main__":
    unittest.main()
