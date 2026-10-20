"""Exercise the browser pager with the embedded renderer and live widget messages."""

import itertools
import re
import sys
import unicodedata
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import patch

from symbolica import E, MathMLViewer, S
from symbolica.community.tensor import (
    DisplaySettings,
    Representation,
    TensorExpression,
    TensorName,
)


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
            rep("i"), Representation.mink(4)("mu")
        ).permute_axes([1, 0])
        original, axes = value.to_expression(), value.axes
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

    def test_renderer_errors_stay_in_the_bounded_preview(self):
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
            page = pager._page()
        self.assertIn("Page rendering failed: RuntimeError: bad notation", page["html"])
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
