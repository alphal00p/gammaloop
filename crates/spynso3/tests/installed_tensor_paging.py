"""Exercise direct tensor MathML in the shared responsive expression viewer."""

import builtins
import itertools
import re
import sys
import unicodedata
import unittest
import xml.etree.ElementTree as ET
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

from symbolica import E, MathMLViewer, PrintMode, S
from symbolica.community.tensor import (
    AUTO,
    DisplaySettings,
    Nc,
    Representation,
    TensorExpression,
    TensorName,
    dot,
    trace,
)


class PagingAssertions(unittest.TestCase):
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

    @staticmethod
    def visibleNodes(node):
        if node.tag.rsplit("}", 1)[-1] not in (
            "mphantom",
            "annotation",
            "annotation-xml",
        ):
            yield node
            for child in node:
                yield from PagingAssertions.visibleNodes(child)

    @staticmethod
    def visible(node):
        text = "".join(n.text or "" for n in PagingAssertions.visibleNodes(node))
        normalized = unicodedata.normalize("NFKC", text).translate(
            dict.fromkeys(map(ord, "\u2061\u2062\u2063\u2064"))
        )
        return "".join(normalized.split())

    @staticmethod
    def tensorState(value):
        return (
            bytes(value.__getstate__()),
            value.to_expression(),
            value.axes,
            value.structure,
            value.name.to_expression() if value.name else None,
            value.arguments,
            value.reduction_status,
        )

    def notationRoots(self, value, settings=None):
        original = self.tensorState(value)
        viewer = value.paged(settings=settings)
        self.addCleanup(viewer.close)
        roots = [self.assertPage(viewer._render_page(layout))[0] for layout in (0, 1)]
        self.assertEqual(self.tensorState(value), original)
        return roots


class DirectLeafPathTests(PagingAssertions):
    def setUp(self):
        original_import = builtins.__import__

        def import_direct(name, *args, **kwargs):
            if name.split(".", 1)[0] in ("typst", "_spenso_paging"):
                raise AssertionError(f"direct paging imported {name}")
            return original_import(name, *args, **kwargs)

        self.enterContext(patch("builtins.__import__", side_effect=import_direct))
        for entry in ("subprocess.run", "subprocess.Popen"):
            self.enterContext(
                patch(entry, side_effect=AssertionError("paging spawned a process"))
            )
        self.compiler = self.enterContext(
            patch(
                "symbolica.community.tensor_native._compile_typst",
                side_effect=AssertionError("paging invoked the Typst compiler"),
            )
        )

    def test_native_and_browser_pages_need_no_typst_module_worker_or_compiler(self):
        gamma = TensorExpression.dirac_gamma(4)
        for platform in (sys.platform, "emscripten"):
            with (
                self.subTest(platform=platform),
                patch.object(sys, "platform", platform),
            ):
                value = gamma("a", "b", AUTO)
                viewer = value.paged()
                self.addCleanup(viewer.close)
                self.assertIs(type(viewer), MathMLViewer)
                root = self.assertPage(viewer._page())[0]
                self.assertEqual(self.visible(root), "γab□")
                self.assertPage(viewer._action({"action": "size", "value": 7}))
        self.compiler.assert_not_called()
        self.assertNotIn("_spenso_paging", sys.modules)

    def test_callable_notation_uses_mathml_without_requesting_typst(self):
        modes = []

        def printer(expression, *, mode, **options):
            modes.append(mode)
            if mode != PrintMode.MathML:
                raise AssertionError("pager requested a non-MathML custom printer")
            return "<mover><mi>J</mi><mo>¯</mo></mover>"

        head = S(
            "paging_direct::OnlyMathML",
            is_scalar=True,
            tags=["spenso::print-callback"],
            print=printer,
        )
        value = TensorExpression(head(7))
        for root in self.notationRoots(value):
            self.assertIn("J", self.visible(root))
            self.assertTrue(
                any(n.tag.rsplit("}", 1)[-1] == "mover" for n in root.iter())
            )
        self.assertTrue(modes)
        self.assertEqual(set(modes), {PrintMode.MathML})

    def test_paging_rejects_every_non_none_custom_typst_notation(self):
        for value in (TensorExpression(E("x+1")), TensorExpression.dirac_gamma(4)):
            original = self.tensorState(value)
            for notation in ("", '#panic("custom-notation-must-stay-explicit")'):
                with self.subTest(notation=notation, rank=value.rank):
                    with self.assertRaises(ValueError) as error:
                        value.paged(notation_source=notation)
                    message = str(error.exception)
                    self.assertIn("notation_source", message)
                    self.assertRegex(message, r"to_html|to_svg|export")
            self.assertEqual(self.tensorState(value), original)

    def test_oversized_callable_keeps_the_tensor_and_all_scalar_arguments(self):
        calls = []

        def printer(expression, *, mode, **options):
            if mode == PrintMode.MathML:
                calls.append(expression)
                return "<mtext>" + "x" * 5000 + "</mtext>"
            return None

        value = TensorName("paging_direct::Fallback", print=printer)(1, 2)
        for root in self.notationRoots(value):
            text = self.visible(root)
            self.assertIn("Fallback", text)
            self.assertIn("1", text)
            self.assertIn("2", text)
            self.assertNotIn("⋯", text)
            self.assertFalse(
                any("data-symbolica-marker" in n.attrib for n in root.iter())
            )
        self.assertTrue(calls)

    def test_leaf_cache_distinguishes_exponents_from_function_arguments(self):
        def printer(expression, *, mode, in_exp, **options):
            if mode == PrintMode.MathML:
                return f"<mtext>{'exponent' if in_exp else 'ordinary'}</mtext>"
            return None

        leaf = TensorName("paging_exponent::F", print=printer)(1).to_expression()
        function = S("paging_exponent::f", is_scalar=True)
        value = TensorExpression(function(leaf) + S("paging_exponent::x") ** leaf)
        for root in self.notationRoots(value):
            text = self.visible(root)
            self.assertIn("ordinary", text)
            self.assertIn("exponent", text)

    def test_cooked_function_indices_keep_their_public_source_labels(self):
        index = S("paging_cooked::i")
        vector = TensorName.vector("paging_cooked::v")(Representation.euc(3))
        source = (
            vector(index).to_expression().replace(index, S("paging_cooked::edge")(7, 2))
        )
        value = TensorExpression(source, intern="indices")
        for root in self.notationRoots(value, DisplaySettings(index_style="raw")):
            text = self.visible(root)
            self.assertIn("edge", text)
            self.assertIn("7", text)
            self.assertIn("2", text)
            self.assertNotIn("cooked_", text)

    def test_native_custom_mathml_preserves_annotations_and_significant_text(self):
        text = "  a b \t \u00a0\u2003"
        fragment = (
            '<semantics><mrow><mtext xml:space="preserve">' + text + "</mtext>"
            "<msub><mi>v</mi><mi>mu</mi></msub></mrow>"
            '<annotation encoding="text/plain">  tensor \t </annotation>'
            '<annotation-xml encoding="application/xhtml+xml">'
            '<span xmlns="http://www.w3.org/1999/xhtml">  label </span>'
            "</annotation-xml></semantics>"
        )

        def printer(expression, *, mode, **options):
            return fragment if mode == PrintMode.MathML else None

        head = S(
            "paging_direct::Annotated",
            is_scalar=True,
            tags=["spenso::print-callback"],
            print=printer,
        )
        for root in self.notationRoots(TensorExpression(head(1))):
            by_tag = {n.tag.rsplit("}", 1)[-1]: n for n in root.iter()}
            self.assertEqual(by_tag["mtext"].text, text)
            self.assertEqual(by_tag["annotation"].text, "  tensor \t ")
            self.assertEqual(by_tag["span"].text, "  label ")
            self.assertIn("msub", by_tag)
            self.assertNotIn("tensor", self.visible(root))


class BrowserPagingTests(PagingAssertions):
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
        value = TensorName(
            "paging_ports::A",
            print={
                "typst": "macron(A)",
                "mathml": "<mover><mi>A</mi><mo>¯</mo></mover>",
            },
        )(
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
        first_text = self.visible(self.assertPage(first)[0])
        explicit_text = self.visible(self.assertPage(explicit)[0])
        self.assertNotEqual(first_text, explicit_text)
        self.assertIn("3", explicit_text)
        self.assertIn("4", explicit_text)
        self.assertEqual(ports._page()["html"], first["html"])
        self.assertEqual(bytes(value.__getstate__()), source)
        self.assertEqual(value.axes, axes)
        self.assertEqual(value.structure, structure)

    def test_leaf_cache_distinguishes_product_and_power_contexts(self):
        def printer(expression, *, mode, in_product, in_exp_base, **options):
            if mode != PrintMode.MathML:
                return None
            label = "base" if in_exp_base else "factor" if in_product else "single"
            return f"<mtext>{label}</mtext>"

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
            text = self.visible(self.assertPage(viewer._render_page(layout))[0])
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
        value = TensorName(
            "paging_direct::Oversized",
            print={"mathml": "<mtext>" + "μ" * 3000 + "</mtext>"},
        )(Representation.mink(4)("mu"), Representation.euc(3)("i"))
        original = self.tensorState(value)
        viewer = value.paged()
        self.addCleanup(viewer.close)
        page = viewer._page()
        self.assertPage(page)
        self.assertEqual((page["kind"], page["unit"]), ("Function", "arguments"))
        self.assertNotIn("μ" * 3000, page["html"])
        self.assertEqual(self.tensorState(value), original)

    def test_callable_errors_follow_the_native_symbolica_fallback(self):
        def printer(expression, *, mode, **options):
            if mode == PrintMode.MathML:
                raise RuntimeError("intentional-native-notation-failure")

        head = S(
            "paging_direct::Broken",
            is_scalar=True,
            tags=["spenso::print-callback"],
            print=printer,
        )
        value = TensorExpression(head(1))
        original = self.tensorState(value)
        viewer = value.paged()
        self.addCleanup(viewer.close)
        # Symbolica's registered Python callback contract logs exceptions and
        # declines the override. The scoped pager preserves that native fallback.
        native = value.to_expression().mathml_viewer()
        self.addCleanup(native.close)
        for layout in (0, 1):
            page = viewer._render_page(layout)
            expected = native._render_page(layout)
            self.assertEqual(page["html"], expected["html"])
            self.assertIn("Broken", self.visible(self.assertPage(page)[0]))
        self.assertEqual(self.tensorState(value), original)

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
        fragment = (
            "<semantics><mstyle><mrow><mi>T</mi>"
            "<msub><mi>v</mi><mi>mu</mi></msub></mrow></mstyle>"
            "<annotation encoding='text/plain'>tensor</annotation></semantics>"
        )

        def printer(expression, *, mode, **options):
            return fragment if mode == PrintMode.MathML else None

        head = S(
            "paging_annotations::T",
            is_scalar=True,
            tags=["spenso::print-callback"],
            print=printer,
        )
        for root in self.notationRoots(TensorExpression(head(1))):
            annotations = [
                n for n in root.iter() if n.tag.rsplit("}", 1)[-1] == "annotation"
            ]
            self.assertEqual([node.text for node in annotations], ["tensor"])
            self.assertTrue(
                any(n.tag.rsplit("}", 1)[-1] == "msub" for n in root.iter())
            )
            self.assertEqual(self.visible(root), "Tvmu")

    def test_generated_aliases_reserve_free_labels_and_stay_viewer_local(self):
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
        while page["previous"]:
            page = viewer._action({"action": "previous"})
        self.assertEqual(page["html"], first)
        self.assertEqual(value.to_expression(), original)
        self.assertEqual(value.structure, interface)
        fresh = value.paged(settings=DisplaySettings(index_style="alphabet"))
        self.addCleanup(fresh.close)
        self.assertEqual(fresh._page()["html"], first)


class TensorNotationTests(PagingAssertions):
    @staticmethod
    def notationCases():
        mink, spin = Representation.mink(4), Representation.bis(4)
        gamma = TensorExpression.dirac_gamma(4)
        p, q = [
            TensorName.vector(f"paging_matrix::{name}")(mink) for name in ("p", "q")
        ]
        left = (gamma(AUTO, AUTO, "mu") * p("mu")).contract()
        right = (gamma(AUTO, AUTO, "nu") * q("nu")).contract()
        dual = Representation("paging_matrix::V", 3, is_self_dual=False)
        bra = TensorName.vector("paging_matrix::a")(dual.dual())
        ket = TensorName.vector("paging_matrix::b")(dual)
        matrix = TensorName("paging_matrix::M")
        tensor = TensorName("paging_matrix::T")(mink, mink)
        cases = [
            ("gamma-indexed", gamma("a", "b", "mu"), None, "γabmu", 0),
            ("gamma-open", gamma("a", "b", AUTO), None, "γab□", 0),
            ("gamma-zero", TensorExpression.gamma0(4)("a", "b"), None, "γ0ab", 0),
            ("gamma-five", TensorExpression.gamma5(4)("a", "b"), None, "γ5ab", 0),
            ("projector-plus", TensorExpression.projp(4)("a", "b"), None, "P+ab", 0),
            ("projector-minus", TensorExpression.projm(4)("a", "b"), None, "P−ab", 0),
            ("slash-indexed", left("a", "b"), None, "pab", 1),
            ("slash-open", left, None, "p□0□1", 1),
            ("slash-chain", left * right, None, "[pq]□0□1", 2),
            ("slash-trace", trace(spin, left, right), None, "Tr(pq)", 2),
            ("gamma-chain", gamma * gamma, None, "[γ□γ□]□0□1", 0),
            (
                "bra-ket",
                matrix(bra.to_expression(), ket.to_expression(), dual),
                None,
                "⟨a|M▶◀□|b⟩",
                0,
            ),
            (
                "supplied-port",
                (tensor("mu", AUTO) * p("mu")).contract(),
                None,
                "⟨p|T■□",
                0,
            ),
            (
                "schoonschip",
                (tensor("mu", AUTO) * p("mu")).contract(),
                DisplaySettings.schoonschip(),
                "Tp□",
                0,
            ),
        ]
        fundamental, adjoint = Representation.cof(3), Representation.coad(8)
        for name, expression, settings in (
            ("casimir-fundamental", fundamental.casimir(), None),
            ("casimir-adjoint", adjoint.casimir(), None),
            ("dynkin-index", fundamental.dynkin_index(), None),
            ("color-gram", fundamental.gram(4, adjoint), None),
            (
                "color-dimensions",
                fundamental.casimir(),
                DisplaySettings(show_dimensions=True),
            ),
            (
                "color-explicit",
                fundamental.casimir(),
                DisplaySettings(invariant_style="explicit"),
            ),
            ("color-unscripted", Nc(), DisplaySettings(symbol_scripts=False)),
        ):
            cases.append((name, TensorExpression(expression), settings, None, 0))
        euc = Representation.euc(3)
        indexed = TensorName("paging_matrix::A")(euc("i"), mink("mu"))
        component = TensorExpression(
            TensorName("paging_matrix::C")(euc, euc).components()[1]
        )
        custom = TensorName(
            "paging_matrix::Custom",
            print={
                "mathml": "<mover><mi>A</mi><mo>¯</mo></mover>",
                "typst": "macron(A)",
            },
        )(euc("i"))
        graph_index = S(
            "paging_indices::hedge",
            tags=["spenso::index", "spenso::index-label:h"],
        )
        TensorName("paging_matrix::Scoped")
        scoped_source = S("paging_matrix::Scoped")(
            S("spenso::mink")(4, graph_index(2, 1))
        )
        scoped = TensorExpression(scoped_source).wrap_indices(S("paging_matrix::scope"))
        cases.extend(
            [
                ("metric", TensorExpression.g(mink("i"), mink("j")), None, None, 0),
                (
                    "representation",
                    TensorExpression(mink.to_expression()),
                    None,
                    None,
                    0,
                ),
                ("layout-ports", indexed, DisplaySettings.ports(), None, 0),
                ("layout-call", indexed, DisplaySettings.call(), None, 0),
                (
                    "layout-dimensions",
                    indexed,
                    DisplaySettings(show_dimensions=True),
                    None,
                    0,
                ),
                (
                    "index-spacing",
                    indexed,
                    DisplaySettings(index_gap="0.37em"),
                    None,
                    0,
                ),
                (
                    "factor-spacing",
                    left * right,
                    DisplaySettings(factor_gap="0.51em"),
                    None,
                    2,
                ),
                ("component-superscript", component, None, None, 0),
                (
                    "component-array",
                    component,
                    DisplaySettings(component_style="array"),
                    None,
                    0,
                ),
                ("custom-mathml-head", custom, None, None, 0),
                ("power-grouping", dot(p, q) ** 2, None, None, 0),
                ("typed-zero", 0 * indexed, None, "0", 0),
            ]
        )
        for style in ("alphabet", "graph", "raw"):
            cases.append(
                (f"index-{style}", scoped, DisplaySettings(index_style=style), None, 0)
            )
        return cases

    def test_full_notation_matrix_in_responsive_layouts(self):
        for name, value, settings, expected, strokes in self.notationCases():
            with self.subTest(name=name):
                for root in self.notationRoots(value, settings):
                    if expected is not None:
                        self.assertEqual(self.visible(root), expected)
                    self.assertTrue(self.visible(root))
                    cancelled = [
                        n
                        for n in self.visibleNodes(root)
                        if n.get("data-spenso-cancel") == "updiagonal"
                    ]
                    self.assertEqual(len(cancelled), strokes)
                    for node in cancelled:
                        self.assertIn("linear-gradient", node.get("style", ""))
                        self.assertIn("currentColor", node.get("style", ""))

    def test_metrics_representations_and_permuted_logical_ports(self):
        mink, euc = Representation.mink(4), Representation.euc(3)
        dual = Representation("paging_metric::V", 5, is_self_dual=False)
        for rep in (mink, euc, dual):
            paired = rep.dual() if rep == dual else rep
            for first, second in ((rep("i"), paired("j")), (rep, paired)):
                metric = TensorExpression.g(first, second)
                for root in self.notationRoots(metric):
                    text = self.visible(root)
                    self.assertTrue("g" in text or "δ" in text, text)
                    if first is rep:
                        ports = [
                            n
                            for n in self.visibleNodes(root)
                            if "data-spenso-open-axis" in n.attrib
                        ]
                        self.assertEqual(len(ports), 2)
                    else:
                        self.assertIn("i", text)
                        self.assertIn("j", text)
        # The factory retains AUTO ports; permutation materializes their
        # assignments. The renderer follows each carrier's established metadata.
        value = TensorExpression.dirac_gamma(4)
        for root in self.notationRoots(value):
            ports = [
                n
                for n in self.visibleNodes(root)
                if "data-spenso-open-axis" in n.attrib
            ]
            self.assertCountEqual(
                [n.get("data-spenso-open-axis") for n in ports], ["0", "1", "2"]
            )
            self.assertTrue(all("spenso::" in n.get("title", "") for n in ports))
        permuted = value.permute_axes([2, 0, 1])
        for root in self.notationRoots(permuted):
            self.assertEqual(self.visible(root), "γμab")
            self.assertFalse(
                any("data-spenso-open-axis" in n.attrib for n in root.iter())
            )
        labelled = TensorName("paging_metric::Labelled")(
            S("paging_metric::x"), 7, mink("mu"), euc("i")
        )
        for root in self.notationRoots(labelled.permute_axes([1, 0])):
            self.assertIn("7", self.visible(root))
        for rep, family, dimension in ((mink, "mink", "4"), (euc, "euc", "3")):
            value = TensorExpression(rep.to_expression())
            for root in self.notationRoots(value):
                self.assertIn(family, self.visible(root))
                self.assertIn(dimension, self.visible(root))

    def test_color_invariants_keep_degree_representations_and_settings(self):
        fundamental, adjoint = Representation.cof(3), Representation.coad(8)
        cases = [
            (fundamental.casimir(), None, "CF"),
            (adjoint.casimir(), None, "CA"),
            (fundamental.dynkin_index(), None, "TR"),
            (Nc(), None, "Nc"),
            (
                fundamental.casimir(),
                DisplaySettings(invariant_style="explicit"),
                "C2(F)",
            ),
            (fundamental.casimir(4), None, "C4(F)"),
            (fundamental.dynkin_index(4), None, "I4(F)"),
            (fundamental.gram(4, adjoint), None, "G4(F,A)"),
            (fundamental.casimir(), DisplaySettings(show_dimensions=True), "C2(F3)"),
            (fundamental.casimir(), DisplaySettings(symbol_scripts=False), "CF"),
        ]
        for expression, settings, expected in cases:
            with self.subTest(expected=expected, settings=settings):
                for root in self.notationRoots(TensorExpression(expression), settings):
                    self.assertEqual(self.visible(root), expected)

    def test_layout_spacing_commas_dimensions_and_component_styles(self):
        rep = Representation.euc(3)
        value = TensorName("paging_layout::A")(rep("i"), rep("j"))
        settings = DisplaySettings(
            tensor_layout="call", commas=True, show_dimensions=True
        )
        for root in self.notationRoots(value, settings):
            text = self.visible(root)
            self.assertIn("(", text)
            self.assertIn(",", text)
            self.assertIn("3", text)
        without = DisplaySettings(tensor_layout="call", commas=False)
        for root in self.notationRoots(value, without):
            self.assertNotIn(",", self.visible(root))
        spaced = DisplaySettings(index_gap="0.37em", factor_gap="0.51em")
        indexed = self.notationRoots(value, spaced)[0]
        widths = [
            n.get("width")
            for n in indexed.iter()
            if n.tag.rsplit("}", 1)[-1] == "mspace"
        ]
        self.assertIn("0.37em", widths)
        p, q = [TensorName.vector(f"paging_layout::{name}")(rep) for name in ("p", "q")]
        contraction = dot(p, q)
        for root in self.notationRoots(contraction**2, spaced):
            powers = [
                n for n in self.visibleNodes(root) if n.tag.rsplit("}", 1)[-1] == "msup"
            ]
            self.assertEqual(len(powers), 1)
            self.assertIn(self.visible(powers[0][0]), ("(p·q)", "(p⋅q)"))
            self.assertEqual(self.visible(powers[0][1]), "2")
        gamma = TensorExpression.dirac_gamma(4)
        # Factor gaps apply inside semantic chains. Arithmetic between distinct
        # scalar factors remains owned by Symbolica's native renderer.
        product = gamma * gamma
        for root in self.notationRoots(product, spaced):
            widths = [
                n.get("width")
                for n in root.iter()
                if n.tag.rsplit("}", 1)[-1] == "mspace"
            ]
            self.assertIn("0.51em", widths)
        component = TensorExpression(value.components()[1])
        superscript = self.notationRoots(component)[0]
        array = self.notationRoots(component, DisplaySettings(component_style="array"))[
            0
        ]
        self.assertTrue(
            any(n.tag.rsplit("}", 1)[-1] == "msup" for n in superscript.iter())
        )
        self.assertIn("[", self.visible(array))
        self.assertIn("0", self.visible(array))
        self.assertIn("1", self.visible(array))

    def test_index_styles_scopes_and_typed_zero_preserve_all_source_metadata(self):
        TensorName("paging_scopes::T")
        graph_index = S(
            "paging_indices::hedge",
            tags=["spenso::index", "spenso::index-label:h"],
        )
        value = TensorExpression(
            S("paging_scopes::T")(
                S("spenso::mink")(4, graph_index(2, 1)),
                S("spenso::mink")(4, S("free")),
            )
        )
        wrapped = value.wrap_indices(S("paging_scopes::scope"))
        texts = {}
        for style in ("alphabet", "graph", "raw"):
            roots = self.notationRoots(wrapped, DisplaySettings(index_style=style))
            texts[style] = self.visible(roots[0])
            self.assertIn("free", texts[style])
        self.assertIn("μ", texts["alphabet"])
        self.assertIn("′", texts["alphabet"])
        self.assertNotIn("hedge", texts["alphabet"])
        self.assertIn("2", texts["graph"])
        self.assertIn("1", texts["graph"])
        self.assertNotEqual(texts["alphabet"], texts["raw"])
        zero = (0 * value).with_name("paging_scopes::NamedZero")
        for root in self.notationRoots(zero):
            self.assertEqual(self.visible(root), "0")
        self.assertEqual(zero.rank, 2)
        self.assertEqual(zero.axes, value.axes)

    def test_cached_fragments_are_path_neutral_and_reused_across_layouts(self):
        fragment = '<mrow data-symbolica-path="97"><mi>A</mi><mo>′</mo></mrow>'
        tensor = TensorName("paging_paths::T", print={"mathml": fragment})(
            Representation.mink(4)("mu")
        )
        f, g = S("paging_paths::F", "paging_paths::G", is_scalar=True)
        value = TensorExpression(f(tensor.to_expression()) + g(tensor.to_expression()))
        viewer = value.paged()
        self.addCleanup(viewer.close)
        for layout in (0, 1, 0):
            page = viewer._render_page(layout)
            root = self.assertPage(page)[0]
            self.assertEqual(self.visible(root).count("A"), 2)
            self.assertNotIn("97", [n.get("data-symbolica-path") for n in root.iter()])

    def test_multi_root_custom_head_stays_inside_component_power_base(self):
        head = TensorName(
            "paging_grouping::A", print={"mathml": "<mi>A</mi><mo>′</mo>"}
        )
        component = TensorExpression(head(Representation.euc(3)).components()[1])
        for root in self.notationRoots(component**2):
            powers = [
                node
                for node in self.visibleNodes(root)
                if node.tag.rsplit("}", 1)[-1] == "msup"
                and len(node) == 2
                and self.visible(node[1]) == "2"
            ]
            self.assertEqual(len(powers), 1)
            base = self.visible(powers[0][0])
            self.assertIn("A′", base)
            self.assertIn("1", base)
            self.assertEqual(self.visible(root).count("′"), 1)

    def test_native_leaf_cache_reuses_recent_entries_and_evicts_old_entries(self):
        calls = {}

        def printer(expression, *, mode, **options):
            if mode != PrintMode.MathML:
                return None
            key = expression.format_plain()
            calls[key] = calls.get(key, 0) + 1
            return "<mi>leaf</mi>"

        head = TensorName("paging_cache::F", print=printer)
        value = TensorExpression(sum(head(i) for i in range(140)))
        original = self.tensorState(value)
        viewer = value.paged(page_size=1)
        self.addCleanup(viewer.close)
        first = viewer._render_page(0)
        self.assertPage(first)
        recent = dict(calls)
        self.assertTrue(recent)
        self.assertPage(viewer._render_page(0))
        self.assertEqual(calls, recent)
        page = viewer._page()
        while page["next"] is not None:
            page = viewer._action({"action": "next"})
            self.assertPage(page)
        self.assertGreaterEqual(len(calls), 140)
        while page["previous"]:
            page = viewer._action({"action": "previous"})
        self.assertEqual(page["html"], first["html"])
        self.assertTrue(any(calls[key] > count for key, count in recent.items()))
        self.assertEqual(self.tensorState(value), original)


class ExplicitNotationExportTests(PagingAssertions):
    def test_custom_typst_notation_remains_complete_in_explicit_exports(self):
        assets = Path(__file__).resolve().parents[1] / "typst"
        notation = (assets / "notation.typ").read_text() + (
            "\n#let default-notation(settings: (:)) = notation(settings: settings, "
            'heads: ("paging_exports::A": $macron(Z)$))\n'
        )
        head = TensorName(
            "paging_exports::A",
            print={"typst": "macron(A)", "mathml": "<mi>native</mi>"},
        )
        value = TensorExpression(sum(head(i) for i in range(35)))
        original = self.tensorState(value)
        self.assertIn("macron(A)", value.to_typst())
        html = value.to_html(notation_source=notation)
        root = ET.fromstring(re.search(r"<math\b.*?</math>", html, re.DOTALL)[0])
        self.assertEqual(self.visible(root).count("Z"), 35)
        self.assertTrue(any(n.tag.rsplit("}", 1)[-1] == "mover" for n in root.iter()))
        self.assertIn("<svg", value.to_svg(notation_source=notation))
        formatted = value.formatted(notation_source=notation)
        self.assertIsNone(formatted.mathml_viewer())
        formatted_root = ET.fromstring(
            re.search(r"<math\b.*?</math>", formatted._repr_html_(), re.DOTALL)[0]
        )
        self.assertEqual(self.visible(formatted_root).count("Z"), 35)
        for render in (value.to_html, value.to_svg):
            with self.assertRaisesRegex(RuntimeError, "explicit-notation-failure"):
                render(notation_source='#panic("explicit-notation-failure")')
        self.assertEqual(self.tensorState(value), original)


def browser_fixtures():
    """Capture fresh native pages and hover responses only when called explicitly."""
    cases = []
    with ExitStack() as cleanup:
        for (
            name,
            value,
            settings,
            expected,
            strokes,
        ) in TensorNotationTests.notationCases():
            original = PagingAssertions.tensorState(value)
            viewer = value.paged(settings=settings)
            cleanup.callback(viewer.close)
            cases.append(
                {
                    "name": name,
                    "page": viewer._page(),
                    "expected": expected,
                    "slash_count": strokes,
                }
            )
            assert PagingAssertions.tensorState(value) == original

        rep = Representation.mink(4)
        value = (
            TensorName("paging_browser::T")(rep("mu"), rep("nu"))
            * (E("paging_browser::x+1") ** 300).expand()
            * (E("paging_browser::y+1") ** 250).expand()
        )
        original = PagingAssertions.tensorState(value)
        viewer = value.paged()
        cleanup.callback(viewer.close)
        page = viewer._page()
        esm = viewer._get_widget()._esm
        previews = {}
        for layout in [page, *page["layouts"]]:
            for (path, label), start in zip(layout["markers"], layout["marker_starts"]):
                key = f"{layout['layout']}:{','.join(map(str, path))}:{start}"
                field, hover = viewer._message(
                    {"action": "preview", "value": path, "layout": layout["layout"]}
                )
                inspect_field, inspected = viewer._message(
                    {
                        "action": "preview",
                        "value": path,
                        "layout": layout["layout"],
                        "popup": True,
                        "size": page["page_size"],
                    }
                )
                previews[key] = {
                    "path": path,
                    "start": start,
                    "label": label,
                    "layout": layout["layout"],
                    "field": field,
                    "page": hover,
                    "inspect_field": inspect_field,
                    "inspect_page": inspected,
                }
        compact = page["layouts"][-1]
        target = next(path for path, label in compact["markers"] if "sum" in label)
        opened = viewer._action(
            {"action": "open", "value": target, "layout": compact["layout"]}
        )
        advanced = viewer._action({"action": "next"})
        restored = viewer._action({"action": "back"})
        assert restored["html"] == page["html"]
        assert PagingAssertions.tensorState(value) == original
        cases.append(
            {
                "name": "nested-pagination",
                "page": page,
                "previews": previews,
                "opened": opened,
                "next": advanced,
                "back": restored,
            }
        )
    return {"esm": esm, "cases": cases}


if __name__ == "__main__":
    unittest.main()
