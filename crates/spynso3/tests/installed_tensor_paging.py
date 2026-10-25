"""Exercise direct tensor MathML in the shared responsive expression viewer."""

import builtins
import importlib.util
import itertools
import json
import re
import sys
import tempfile
import unicodedata
import unittest
import xml.etree.ElementTree as ET
from contextlib import ExitStack
from html.parser import HTMLParser
from pathlib import Path
from types import ModuleType
from unittest.mock import Mock, patch

from symbolica import AtomType, E, Expression, MathMLViewer, PrintMode, S
from symbolica.community.tensor import (
    AUTO,
    DisplaySettings,
    Nc,
    Representation,
    Tensor,
    TensorExpression,
    TensorName,
    dot,
    formatted,
    to_html,
    to_svg,
    to_typst,
    to_typst_project,
    trace,
)
from symbolica.community.tensor_native import _compile_typst

try:
    from anywidget import AnyWidget
except ImportError:
    AnyWidget = None


def graph_index_print(expression, *, mode, **options):
    if (
        mode != PrintMode.MathML
        or expression.get_type() != AtomType.Fn
        or "spenso" not in options.get("custom_print_mode", {})
    ):
        return None
    arguments = [str(argument) for argument in expression]
    assert all(argument.isdecimal() for argument in arguments)
    body = "<mo>.</mo>".join(f"<mn>{argument}</mn>" for argument in arguments)
    return f"<msub><mi>h</mi><mrow>{body}</mrow></msub>"


def graph_index_symbol():
    return TensorName(
        "paging_indices::hedge", tags=["spenso::index"], print=graph_index_print
    ).to_expression()


class TableCells(HTMLParser):
    def __init__(self, document):
        super().__init__()
        self.rows = []
        self.tables = 0
        self.depth = 0
        self.hidden = 0
        self.row = None
        self.cell = None
        self.feed(document)
        self.close()

    def handle_starttag(self, tag, attrs):
        if tag in ("style", "script", "mphantom", "annotation", "annotation-xml"):
            self.hidden += 1
        if tag == "table":
            self.depth += 1
            self.tables += self.depth == 1
        if self.depth == 1:
            if tag == "tr":
                self.row = []
            elif tag == "td":
                self.cell = []

    def handle_endtag(self, tag):
        if tag in ("style", "script", "mphantom", "annotation", "annotation-xml"):
            self.hidden -= 1
        if self.depth == 1:
            if tag == "td":
                self.row.append("".join(self.cell).strip())
                self.cell = None
            elif tag == "tr":
                self.rows.append(self.row)
                self.row = None
        if tag == "table":
            self.depth -= 1

    def handle_data(self, data):
        if self.cell is not None and not self.hidden:
            self.cell.append(data)


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

    def test_loaded_legacy_mathml_is_literal_and_forbidden_elements_are_rejected(self):
        payload = '<img src="x" onerror="legacy_loaded_xss=1"/>'
        head = S(
            "paging_xss::" + json.dumps(payload),
            tags=["spenso::tensor", "spenso::print-head:mathml:" + payload],
        )
        with tempfile.TemporaryDirectory() as directory:
            filename = str(Path(directory) / "expression.bin")
            head(1).save(filename)
            loaded = Expression.load(filename)
        value = TensorExpression(loaded)
        html = value.to_html()
        self.assertNotIn("<img", html)
        self.assertIn("&lt;img", html)
        root = ET.fromstring(re.search(r"<math\b.*?</math>", html, re.DOTALL).group())
        self.assertIn(json.dumps(payload), "".join(root.itertext()))
        self.assertEqual(value.to_expression(), loaded)
        for forbidden in (payload, "<script>alert(1)</script>", "<iframe/>"):
            with self.subTest(forbidden=forbidden), self.assertRaises(ValueError):
                TensorName("paging_xss::Rejected", print={"mathml": forbidden})

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

    def test_complete_html_contains_every_term_beyond_the_native_page_size(self):
        names = [f"leaf{i:03d}" for i in range(301)]
        value = TensorExpression(sum(S(f"complete_html::{name}") for name in names))
        original = self.tensorState(value)
        viewer = value.paged()
        self.addCleanup(viewer.close)
        page = viewer._page()
        self.assertEqual((page["start"], page["end"], page["total"]), (0, 25, 301))
        self.assertPage(page)
        for document in (value.to_html(), to_html(value.to_expression())):
            root = ET.fromstring(
                re.search(r"<math\b.*?</math>", document, re.DOTALL)[0]
            )
            text = self.visible(root)
            self.assertCountEqual(re.findall(r"leaf\d{3}", text), names)
            self.assertNotRegex(text, "[⋯…]")
            self.assertFalse(
                any("data-symbolica-marker" in n.attrib for n in root.iter())
            )
        self.assertEqual(self.tensorState(value), original)

    def test_complete_html_walks_deep_and_wide_expression_children(self):
        names = [f"leaf{i:03d}" for i in range(67)]
        leaves = [S(f"complete_html_tree::{name}") for name in names]
        call = S("complete_html_tree::wrap", is_scalar=True)
        deep = sum(leaves)
        for _ in range(48):
            deep = call(deep)
        cases = {
            "product": E("1"),
            "power": sum(leaves) ** 2,
            "fraction": sum(leaves) / S("complete_html_tree::denominator"),
            "wide-call": call(*leaves),
            "deep-call": deep,
        }
        for leaf in leaves:
            cases["product"] *= leaf
        for name, expression in cases.items():
            with self.subTest(name=name):
                value = TensorExpression(expression)
                original = self.tensorState(value)
                document = value.to_html()
                root = ET.fromstring(
                    re.search(r"<math\b.*?</math>", document, re.DOTALL)[0]
                )
                text = self.visible(root)
                self.assertCountEqual(re.findall(r"leaf\d{3}", text), names)
                self.assertNotRegex(text, "[⋯…]")
                if name == "deep-call":
                    self.assertEqual(text.count("wrap"), 48)
                if name == "wide-call":
                    self.assertEqual(text.count(","), len(leaves) - 1)
                self.assertEqual(self.tensorState(value), original)

    def test_complete_html_preserves_custom_mathml_text_and_annotations(self):
        modes = []
        text = "  a b \t \u00a0\u2003"
        fragment = (
            '<semantics><mrow><mtext xml:space="preserve">' + text + "</mtext>"
            "<mover><mi>J</mi><mo>¯</mo></mover></mrow>"
            '<annotation encoding="text/plain">  tensor \t </annotation>'
            '<annotation-xml encoding="application/xhtml+xml">'
            '<span xmlns="http://www.w3.org/1999/xhtml">  label </span>'
            "</annotation-xml></semantics>"
        )

        def printer(expression, *, mode, **options):
            modes.append(mode)
            if mode != PrintMode.MathML:
                raise AssertionError("HTML requested a non-MathML printer")
            return fragment

        value = TensorName("complete_html::Annotated", print=printer)(7)
        original = self.tensorState(value)
        for document in (value.to_html(), value.formatted()._repr_html_()):
            root = ET.fromstring(
                re.search(r"<math\b.*?</math>", document, re.DOTALL)[0]
            )
            nodes = {n.tag.rsplit("}", 1)[-1]: n for n in root.iter()}
            self.assertEqual(nodes["mtext"].text, text)
            self.assertEqual(nodes["annotation"].text, "  tensor \t ")
            self.assertEqual(nodes["span"].text, "  label ")
            self.assertIn("mover", nodes)
            self.assertNotIn("tensor", self.visible(root))
        self.assertTrue(modes)
        self.assertEqual(set(modes), {PrintMode.MathML})
        self.assertEqual(self.tensorState(value), original)

    def test_concrete_html_cells_follow_logical_axes_after_storage_permutation(self):
        descriptor = TensorName("complete_html::Matrix")(
            Representation.mink(2)("mu"), Representation.euc(3)("i")
        )
        tensor = Tensor.dense(descriptor, [1.0, 20.0, 300.0, 4000.0, 50.0, 6.0])
        for value, expected in (
            (tensor, [["1", "20", "300"], ["4000", "50", "6"]]),
            (tensor.permute_axes([1, 0]), [["1", "4000"], ["20", "50"], ["300", "6"]]),
        ):
            with self.subTest(shape=value.shape):
                original = (
                    value[:],
                    value.structure,
                    value.expression().to_expression(),
                )
                document = value.to_html(settings=DisplaySettings(tensor_view="matrix"))
                table = TableCells(document)
                self.assertEqual(table.tables, 1)
                self.assertEqual(table.rows, expected)
                self.assertNotIn("<iframe", document)
                self.assertEqual(
                    (value[:], value.structure, value.expression().to_expression()),
                    original,
                )

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

    def test_cooked_raw_index_metadata_renders_the_source_without_proxy_symbols(self):
        index = S("review_source::i")
        vector = TensorName.vector("review_source::v")(Representation.euc(3))
        raw = (
            vector(index).to_expression().replace(index, S("review_source::edge")(7, 2))
        )
        value = TensorExpression(raw, intern="indices")
        original = self.tensorState(value)
        settings = DisplaySettings(index_style="raw")
        document = value.structure.to_html(settings=settings)
        roots = [
            ET.fromstring(source)
            for source in re.findall(r"<math\b.*?</math>", document, re.DOTALL)
        ]
        self.assertTrue(roots)
        self.assertEqual([self.visible(root) for root in roots], ["edge(7,2)"])
        for root in roots:
            self.assertNotRegex(
                ET.tostring(root, encoding="unicode"), r"display_source|cooked_"
            )
        self.assertEqual(self.tensorState(value), original)

    def test_cooked_raw_index_source_exceeding_export_depth_raises_limit_error(self):
        index = S("review_deep_source::i")
        vector = TensorName.vector("review_deep_source::v")(Representation.euc(3))
        edge = S("review_deep_source::edge")
        label = S("review_deep_source::leaf")
        for _ in range(257):
            label = edge(label)
        raw = vector(index).to_expression().replace(index, label)
        value = TensorExpression(raw, intern="indices")
        original = self.tensorState(value)
        with self.assertRaisesRegex(ValueError, r"(?i)(export|resource).*limit"):
            value.to_html(settings=DisplaySettings(index_style="raw"))
        self.assertEqual(self.tensorState(value), original)
        self.assertIn("<math", vector(index).to_html())

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

    @unittest.skipUnless(
        AnyWidget is not None,
        "live widget navigation requires the optional anywidget package",
    )
    def test_live_widget_navigation_and_page_size(self):
        widget = self.pager._get_widget()
        self.pager._on_message(widget, {"action": "attach", "view": "browser"}, [])
        self.assertIsInstance(widget, AnyWidget)
        self.assertIs(widget.page["supports_expansion"], True)
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

    def test_nested_permuted_leaf_matches_complete_export(self):
        rep = Representation.mink(4)
        leaf = TensorName("paging_nested_ports::P")(
            Representation.euc(3)("i"), rep("mu")
        ).permute_axes([1, 0])
        value = S("paging_nested_ports::y") * leaf
        original = self.tensorState(value)
        for root in self.notationRoots(value):
            labels = [
                unicodedata.normalize("NFKC", node.text or "")
                for node in self.visibleNodes(root)
                if node.tag.rsplit("}", 1)[-1] in ("mi", "mtext")
                and unicodedata.normalize("NFKC", node.text or "") in ("mu", "i")
            ]
            self.assertEqual(labels, ["mu", "i"])
        static = ET.fromstring(
            re.search(r"<math\b.*?</math>", value.to_html(), re.DOTALL).group()
        )
        self.assertLess(
            self.visible(static).index("mu"), self.visible(static).index("i")
        )
        self.assertEqual(self.tensorState(value), original)

    def test_generated_aliases_reserve_contracted_explicit_labels_viewer_wide(self):
        Representation.mink(4)
        names = [
            TensorName(f"paging_internal_aliases::{name}")
            for name in ("A", "B", "v", "w")
        ]
        a, b, v, w = [name.to_expression() for name in names]
        mink = S("spenso::mink")
        mu = S("mu")
        dummy = S("paging_internal_aliases::d_0", tags=["spenso::dummy-index"])
        value = TensorExpression(
            a(mink(4, mu)) * b(mink(4, mu)) * v(mink(4, dummy)) * w(mink(4, dummy))
        )
        original = self.tensorState(value)
        settings = DisplaySettings(index_style="alphabet")
        static = ET.fromstring(
            re.search(
                r"<math\b.*?</math>", value.to_html(settings=settings), re.DOTALL
            ).group()
        )
        self.assertEqual(self.visible(static).count("ν"), 2)
        for root in self.notationRoots(value, settings):
            self.assertEqual(self.visible(root).count("ν"), 2)
        self.assertEqual(self.tensorState(value), original)

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

    def test_paged_slash_has_semantic_mathml(self):
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
            self.assertEqual(strokes[0].tag.rsplit("}", 1)[-1], "menclose")
            self.assertEqual(strokes[0].get("notation"), "updiagonalstrike")
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
            # A bounded source can still exceed the leaf output bound after escaping.
            print={"mathml": "<mtext>" + "'" * 1000 + "</mtext>"},
        )(Representation.mink(4)("mu"), Representation.euc(3)("i"))
        original = self.tensorState(value)
        viewer = value.paged()
        self.addCleanup(viewer.close)
        page = viewer._page()
        self.assertPage(page)
        self.assertEqual((page["kind"], page["unit"]), ("Function", "arguments"))
        self.assertNotIn("&#39;" * 1000, page["html"])
        self.assertEqual(self.tensorState(value), original)

    def test_callable_errors_follow_the_native_symbolica_fallback(self):
        def printer(expression, *, mode, **options):
            if mode == PrintMode.MathML:
                raise RuntimeError("intentional-native-notation-failure")

        head = S(
            "paging_direct::Broken",
            is_scalar=True,
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
        graph_index = graph_index_symbol()
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
                        self.assertEqual(node.tag.rsplit("}", 1)[-1], "menclose")
                        self.assertEqual(node.get("notation"), "updiagonalstrike")
                        self.assertIn("linear-gradient", node.get("style", ""))
                        self.assertIn("currentColor", node.get("style", ""))

    def test_complete_html_keeps_the_full_notation_matrix_and_source_bytes(self):
        for name, value, settings, expected, strokes in self.notationCases():
            with self.subTest(name=name):
                original = self.tensorState(value)
                document = value.to_html(settings=settings)
                root = ET.fromstring(
                    re.search(r"<math\b.*?</math>", document, re.DOTALL)[0]
                )
                reference = self.notationRoots(value, settings)[0]
                if name == "index-raw":
                    # Raw disables generated aliases. The scoped leaf declines
                    # bounded rendering, so paging retains native call syntax;
                    # complete HTML can keep the tagged index's own notation.
                    self.assertEqual(
                        self.visible(reference),
                        "Scoped(mink(4,index_scope(scope,h2.1)))",
                    )
                    self.assertEqual(
                        self.visible(root), "Scopedindex_scope(scope,h2.1)"
                    )
                else:
                    self.assertEqual(self.visible(root), self.visible(reference))
                if expected is not None:
                    self.assertEqual(self.visible(root), expected)
                cancelled = [
                    n
                    for n in self.visibleNodes(root)
                    if n.get("data-spenso-cancel") == "updiagonal"
                ]
                self.assertEqual(len(cancelled), strokes)
                for node in cancelled:
                    self.assertEqual(node.tag.rsplit("}", 1)[-1], "menclose")
                    self.assertEqual(node.get("notation"), "updiagonalstrike")
                    self.assertIn("linear-gradient", node.get("style", ""))
                    self.assertIn("currentColor", node.get("style", ""))
                self.assertEqual(self.tensorState(value), original)

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
        graph_index = graph_index_symbol()
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
        tensor = TensorName(
            "paging_paths::T",
            print=lambda *args, **kwargs: (
                fragment if kwargs.get("mode") == PrintMode.MathML else None
            ),
        )(Representation.mink(4)("mu"))
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
    def setUp(self):
        # This fake checks the Python compiler boundary, not Typst rendering.
        self.typst = ModuleType("typst")
        self.output = Mock(return_value=b'<svg xmlns="http://www.w3.org/2000/svg"/>')
        self.staged = []

        def compile_project(main, *, root, format):
            self.assertIs(type(main), str)
            self.assertIs(type(root), str)
            main, root = Path(main), Path(root)
            self.assertTrue(main.is_absolute())
            self.assertTrue(root.is_absolute())
            self.assertEqual(main, root / "main.typ")
            self.assertTrue(main.is_file())
            self.assertEqual(format, "svg")
            files = {
                path.relative_to(root).as_posix(): path.read_bytes()
                for path in root.rglob("*")
                if path.is_file()
            }
            self.staged.append((root, files))
            return self.output()

        self.typst.compile = Mock(side_effect=compile_project)
        self.addCleanup(self.assertStagingRemoved)
        self.enterContext(patch.dict(sys.modules, {"typst": self.typst}))

    def assertStagingRemoved(self):
        for root, _ in self.staged:
            self.assertFalse(
                root.exists(), f"compiler staging directory remains: {root}"
            )

    def assertCompilerFiles(self, project):
        self.typst.compile.assert_called_once()
        self.assertEqual(len(self.staged), 1)
        self.assertEqual(self.staged[0][1], project)
        self.assertStagingRemoved()

    def test_custom_notation_is_raw_source_and_an_exact_svg_compiler_input(self):
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
        bare = value.to_typst()
        self.assertEqual(bare.count("macron(A)"), 35)
        project = value.to_typst_project(notation_source=notation)
        self.assertIs(type(project), dict)
        self.assertEqual(
            set(project), {"main.typ", "render.typ", "notation.typ", "tree.cbor"}
        )
        self.assertTrue(
            all(
                type(name) is str and type(data) is bytes
                for name, data in project.items()
            )
        )
        self.assertEqual(project["notation.typ"], notation.encode())
        self.assertEqual(project["render.typ"], (assets / "render.typ").read_bytes())
        self.assertIn(b'cbor(read("tree.cbor", encoding: none))', project["main.typ"])
        self.assertTrue(project["tree.cbor"])
        default_project = value.to_typst_project()
        self.assertEqual(
            default_project["notation.typ"], (assets / "notation.typ").read_bytes()
        )
        for name in ("main.typ", "render.typ", "tree.cbor"):
            self.assertEqual(default_project[name], project[name])
        self.typst.compile.assert_not_called()
        svg = value.to_svg(notation_source=notation)
        self.assertEqual(svg, self.output.return_value.decode())
        self.assertCompilerFiles(project)
        self.assertEqual(value.to_typst(), bare)
        self.assertEqual(self.tensorState(value), original)

    def test_raw_expression_and_concrete_exports_share_the_same_file_contract(self):
        expression = E("portable_exports::x+1")
        descriptor = TensorName("portable_exports::Matrix")(
            Representation.mink(2)("mu"), Representation.euc(3)("i")
        )
        tensor = Tensor.dense(descriptor, [1.0, 20.0, 300.0, 4000.0, 50.0, 6.0])
        notation = "// exact UTF-8: μ, λ\n#let default-notation(settings: (:)) = none\n"
        for project_export, svg_export in (
            (
                lambda **kwargs: to_typst_project(expression, **kwargs),
                lambda **kwargs: to_svg(expression, **kwargs),
            ),
            (tensor.to_typst_project, tensor.to_svg),
        ):
            with self.subTest(export=svg_export):
                self.typst.compile.reset_mock()
                self.staged.clear()
                project = project_export(notation_source=notation)
                self.assertEqual(
                    set(project),
                    {"main.typ", "render.typ", "notation.typ", "tree.cbor"},
                )
                self.assertTrue(all(type(data) is bytes for data in project.values()))
                self.assertEqual(project["notation.typ"], notation.encode())
                self.typst.compile.assert_not_called()
                self.assertEqual(
                    svg_export(notation_source=notation),
                    self.output.return_value.decode(),
                )
                self.assertCompilerFiles(project)

    def test_html_formatted_and_paged_reject_every_supplied_notation_source(self):
        expression = E("explicit_rejection::x+1")
        value = TensorExpression(expression)
        tensor = Tensor.dense(
            TensorName("explicit_rejection::A")(Representation.euc(2)), [1.0, 2.0]
        )
        original = self.tensorState(value)
        for render in (
            value.to_html,
            value.formatted,
            value.paged,
            tensor.to_html,
            tensor.formatted,
            lambda **kwargs: to_html(expression, **kwargs),
            lambda **kwargs: formatted(expression, **kwargs),
        ):
            for notation in ("", '#panic("must-not-compile")'):
                with (
                    self.subTest(render=render, notation=notation),
                    self.assertRaisesRegex(ValueError, "notation_source"),
                ):
                    render(notation_source=notation)
        self.typst.compile.assert_not_called()
        self.assertEqual(self.tensorState(value), original)

    def test_svg_accepts_bytes_or_one_page_and_rejects_other_page_counts(self):
        value = TensorExpression(E("compiler_pages::x+1"))
        for output in (b"<svg/>", [b"<svg/>"]):
            with self.subTest(output=output):
                self.output.return_value = output
                self.assertEqual(value.to_svg(), "<svg/>")
                self.assertStagingRemoved()
        for pages in ([], [b"<svg/>", b"<svg/>"]):
            with self.subTest(pages=len(pages)):
                self.output.return_value = pages
                with self.assertRaisesRegex(RuntimeError, "exactly one page"):
                    value.to_svg()
                self.assertStagingRemoved()

    def test_svg_reports_compiler_diagnostics_and_invalid_output(self):
        value = TensorExpression(E("compiler_errors::x+1"))
        original = self.tensorState(value)
        self.output.side_effect = ValueError("explicit-notation-failure")
        with self.assertRaisesRegex(RuntimeError, "explicit-notation-failure"):
            value.to_svg()
        self.assertStagingRemoved()
        self.output.side_effect = None
        for output in ("<svg/>", ["<svg/>"]):
            with self.subTest(output=output):
                self.output.return_value = output
                with self.assertRaises(TypeError):
                    value.to_svg()
                self.assertStagingRemoved()
        self.output.return_value = b"\xff"
        with self.assertRaisesRegex(RuntimeError, "invalid UTF-8"):
            value.to_svg()
        self.assertStagingRemoved()
        self.output.return_value = b"<svg/>"
        self.assertEqual(value.to_svg(), "<svg/>")
        self.assertStagingRemoved()
        self.assertEqual(self.tensorState(value), original)

    def test_private_compiler_hook_is_svg_only_and_preserves_raw_source(self):
        source, notation, tree = "$ x $", "// custom μ\n", b"\xa0"
        for format in ("html", "pdf", "png", "unknown"):
            with (
                self.subTest(format=format),
                self.assertRaisesRegex(ValueError, "SVG"),
            ):
                _compile_typst(source, format, notation, tree)
        self.typst.compile.assert_not_called()
        output = _compile_typst(source, "svg", notation, tree)
        self.assertEqual(bytes(output), self.output.return_value)
        self.typst.compile.assert_called_once()
        files = self.staged[0][1]
        self.assertEqual(
            set(files), {"main.typ", "render.typ", "notation.typ", "tree.cbor"}
        )
        self.assertEqual(files["main.typ"], source.encode())
        self.assertEqual(files["notation.typ"], notation.encode())
        self.assertEqual(files["tree.cbor"], tree)
        self.assertStagingRemoved()

    def test_missing_optional_compiler_affects_only_explicit_svg(self):
        expression = E("missing_compiler::x+1")
        value = TensorExpression(expression)
        tensor = Tensor.dense(
            TensorName("missing_compiler::A")(Representation.euc(2)), [1.0, 2.0]
        )
        original = self.tensorState(value)
        original_import = builtins.__import__
        attempts = []

        def without_typst(name, *args, **kwargs):
            if name.split(".", 1)[0] == "typst":
                attempts.append(name)
                raise ModuleNotFoundError("No module named 'typst'", name="typst")
            return original_import(name, *args, **kwargs)

        with (
            patch.dict(sys.modules, {"typst": None}),
            patch("builtins.__import__", side_effect=without_typst),
        ):
            self.assertIn("<math", value.to_html())
            self.assertIn("<math", value.formatted()._repr_html_())
            self.assertIn("<math", to_html(expression))
            self.assertIn("<math", formatted(expression)._repr_html_())
            self.assertIn(
                "<table", tensor.to_html(settings=DisplaySettings(tensor_view="matrix"))
            )
            self.assertEqual(to_typst(expression), value.to_typst())
            for project in (
                value.to_typst_project(),
                to_typst_project(expression),
                tensor.to_typst_project(),
            ):
                self.assertEqual(
                    set(project),
                    {"main.typ", "render.typ", "notation.typ", "tree.cbor"},
                )
            viewer = value.paged()
            self.addCleanup(viewer.close)
            self.assertPage(viewer._page())
            self.assertEqual(attempts, [])
            for render in (value.to_svg, lambda: to_svg(expression), tensor.to_svg):
                with self.assertRaisesRegex(
                    ImportError, r"(?s)typst.*pip install typst"
                ):
                    render()
            self.assertIn("<math", value.to_html())
            self.assertTrue(value.to_typst_project()["tree.cbor"])
        self.typst.compile.assert_not_called()
        self.assertEqual(self.tensorState(value), original)


@unittest.skipUnless(
    importlib.util.find_spec("typst") is not None,
    "optional real SVG integration requires the typst package; compiler contracts run unconditionally",
)
class RealCompilerExportTests(PagingAssertions):
    def test_portable_custom_cancellation_retains_its_html_stroke(self):
        import typst

        assets = Path(__file__).resolve().parents[1] / "typst"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "render.typ").write_bytes((assets / "render.typ").read_bytes())
            # Exercise the portable custom-display AST independently of gamma's
            # specialized renderer, including Typst's lexical helper scope.
            (root / "notation.typ").write_text(
                (assets / "notation.typ").read_text()
                + '\n#let custom-cancel() = _display-node(("math", "cancel", (("symbol", "J"),)))\n'
            )
            source = root / "main.typ"
            source.write_text(
                '#import "notation.typ": custom-cancel\n$ #custom-cancel() $'
            )
            document = typst.compile(
                str(source), format="html", root=directory
            ).decode()
            svg = typst.compile(str(source), format="svg", root=directory).decode()
        root = ET.fromstring(re.search(r"<math\b.*?</math>", document, re.DOTALL)[0])
        strokes = [node for node in root.iter() if "data-spenso-cancel" in node.attrib]
        self.assertEqual(len(strokes), 1)
        self.assertEqual(strokes[0].attrib["data-spenso-cancel"], "updiagonal")
        self.assertEqual(
            unicodedata.normalize("NFKC", "".join(strokes[0].itertext())).strip(), "J"
        )
        self.assertIn("<svg", svg)
        self.assertNotIn("data-spenso-cancel", svg)

    def test_custom_typst_notation_compiles_complete_real_svg_projects(self):
        assets = Path(__file__).resolve().parents[1] / "typst"
        notation = (assets / "notation.typ").read_text() + (
            "\n#let default-notation(settings: (:)) = notation(settings: settings, "
            'heads: ("real_exports::A": $macron(Z)$))\n'
        )
        head = TensorName(
            "real_exports::A", print={"typst": "macron(A)", "mathml": "<mi>native</mi>"}
        )
        value = TensorExpression(sum(head(i) for i in range(35)))
        original = self.tensorState(value)
        self.assertEqual(value.to_typst().count("macron(A)"), 35)
        svg = ET.fromstring(value.to_svg(notation_source=notation))
        self.assertEqual(svg.tag.rsplit("}", 1)[-1], "svg")
        with self.assertRaisesRegex(RuntimeError, "explicit-notation-failure"):
            value.to_svg(notation_source='#panic("explicit-notation-failure")')
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
    from symbolica.community import tensor as spenso

    expressions = [
        "(a/(b+c)+d/(e+f))*(g+h)+z",
        "+".join(f"a{i}/(b{i}+c{i})" for i in range(12)),
        "*".join(f"x{i}" for i in range(30)) + "+z",
        "*".join(f"x{i}" for i in range(30)),
        "(" + "*".join(f"x{i}" for i in range(30)) + "+z)/(D+1)",
    ]
    html = [TensorExpression(E(source)).to_html() for source in expressions]
    library = spenso.TensorLibrary.hep_lib_atom()
    explorer = (
        '<section id="standalone">'
        + library[TensorExpression.dirac_gamma(4)].to_html()
        + '</section><section id="library">'
        + library.to_html()
        + "</section>"
    )
    return {
        "esm": esm,
        "cases": cases,
        "font": [html[0], spenso.load_math_font(), *html[1:]],
        "explorer": explorer,
    }


if __name__ == "__main__":
    unittest.main()
