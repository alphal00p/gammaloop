"""Exercise the browser pager with the embedded renderer and live widget messages."""

import itertools
import re
import sys
import unicodedata
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import patch

from symbolica import E, MathMLViewer
from symbolica.community.tensor import Representation, TensorExpression, TensorName


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
        self.assertLessEqual(
            sum(sum(1 for _ in ET.fromstring(m).iter()) for m in math), 10_000
        )
        marker_ids = {
            int(node.attrib["data-symbolica-marker"])
            for expression in math
            for node in ET.fromstring(expression).iter()
            if "data-symbolica-marker" in node.attrib
        }
        self.assertEqual(marker_ids, set(page["marker_ids"]))

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
        self.assertFalse(widget.page["supports_expansion"])
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
        target = parent["markers"][0][0]
        child = pager._action({"action": "open", "value": target})
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
        self.assertFalse(preview["can_expand"])
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
        self.assertIn("<mover", page["html"])
        for html in (value.to_html(), page["html"]):
            root = ET.fromstring(
                re.search(r"<math\b.*?</math>", html, re.DOTALL).group()
            )
            hidden = {
                child for phantom in root.iter("mphantom") for child in phantom.iter()
            }
            indices = [
                unicodedata.normalize("NFKC", node.text or "")
                for node in root.iter()
                if node not in hidden
                and node.tag in ("mi", "mtext")
                and unicodedata.normalize("NFKC", node.text or "") in ("mu", "i")
            ]
            self.assertEqual(indices, ["mu", "i"])
        self.assertEqual(value.axes, axes)
        self.assertEqual(value.to_expression(), original)

    def test_oversized_output_becomes_an_inspectable_placeholder(self):
        bootstrap = TensorExpression(0).paged()
        bootstrap._page()
        bootstrap.close()
        import _spenso_paging

        value = TensorExpression(E("wasm_paging::f(a+b+c)"))
        pager = value.paged()
        self.addCleanup(pager.close)
        with patch.object(
            _spenso_paging, "compile_typst", return_value=b"x" * (256 * 1024 + 1)
        ):
            page = pager._page()
        self.assertIn("exceeds the display budget", page["html"])
        self.assertTrue(page["markers"])
        self.assertEqual(page["end"], 1)
        self.assertEqual(value.to_expression(), E("wasm_paging::f(a+b+c)"))

    def test_renderer_errors_stay_in_the_bounded_preview(self):
        bootstrap = TensorExpression(0).paged()
        bootstrap._page()
        bootstrap.close()
        import _spenso_paging

        with patch.object(
            _spenso_paging, "compile_typst", side_effect=RuntimeError("bad notation")
        ):
            page = self.pager._page()
        self.assertIn("Page rendering failed: RuntimeError: bad notation", page["html"])
        self.assertEqual(self.value.to_expression(), self.original)

    def test_wrapped_terms_never_shrink_or_wrap_inside_math_children(self):
        summands = "+".join(f"layout::a{i}" for i in range(10))
        nested = {
            "product": f"z*({summands})",
            "fraction": f"({summands})/(z+w)",
            "power": f"({summands})^2",
            "function": f"f({summands})",
        }
        for label, expression in nested.items():
            with self.subTest(label=label):
                pager = TensorExpression(E(expression)).paged()
                self.addCleanup(pager.close)
                page = pager._page()
                self.assertPage(page)
                roots = [
                    ET.fromstring(m)
                    for m in re.findall(r"<math\b.*?</math>", page["html"], re.DOTALL)
                ]
                self.assertFalse(
                    any(
                        "data-spenso-sum" in node.attrib
                        for root in roots
                        for node in root.iter()
                    )
                )
        pager = TensorExpression(E(summands)).paged()
        self.addCleanup(pager.close)
        page = pager._page()
        self.assertPage(page)
        root = ET.fromstring(
            re.search(r"<math\b.*?</math>", page["html"], re.DOTALL).group()
        )
        self.assertIn("data-spenso-sum", root.attrib)
        self.assertTrue(all(group.get("style") == "flex:0 0 auto" for group in root))

    def test_outer_sum_wrapper_chain_preserves_semantic_annotations(self):
        bootstrap = TensorExpression(0).paged()
        bootstrap._page()
        bootstrap.close()
        import _spenso_paging

        terms = "<mo>+</mo>".join(f"<mi>a{i}</mi>" for i in range(10))
        document = (
            "<html><body><math><semantics><mstyle><mrow>"
            + terms
            + "</mrow></mstyle><annotation encoding='text/plain'>sum</annotation>"
            + "</semantics></math></body></html>"
        )
        pager = TensorExpression(E("x")).paged()
        self.addCleanup(pager.close)
        with patch.object(
            _spenso_paging, "compile_typst", return_value=document.encode()
        ):
            page = pager._page()
        self.assertPage(page)
        root = ET.fromstring(
            re.search(r"<math\b.*?</math>", page["html"], re.DOTALL).group()
        )
        row = root.find("semantics/mstyle/mrow")
        self.assertIn("data-spenso-sum", row.attrib)
        self.assertEqual(root.find("semantics/annotation").text, "sum")
        self.assertTrue(all(group.get("style") == "flex:0 0 auto" for group in row))


if __name__ == "__main__":
    unittest.main()
