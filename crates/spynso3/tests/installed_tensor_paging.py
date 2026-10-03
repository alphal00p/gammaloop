"""Exercise the browser pager with the embedded renderer and live widget messages."""

import itertools
import re
import sys
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import patch

from symbolica import E
from symbolica.community.tensor import TensorExpression


class BrowserPagingTests(unittest.TestCase):
    def setUp(self):
        self.value = TensorExpression(
            E("+".join(f"{i + 1}*wasm_paging::x^{i}" for i in range(301)))
        )
        self.original = self.value.to_expression()
        self.pager = self.value.paged()
        self.addCleanup(self.pager.close)
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
        self.assertLessEqual(len(self.pager._cache), 3)
        for start, end, html in reversed(pages[:-1]):
            page = self.pager._action({"action": "previous"})
            self.assertEqual(
                (page["start"], page["end"], page["html"]), (start, end, html)
            )
        self.assertEqual(self.value.to_expression(), self.original)

    def test_live_widget_navigation_and_page_size(self):
        widget = self.pager._get_widget()
        self.pager._on_message(widget, {"action": "attach", "view": "browser"}, [])
        self.assertEqual(widget.page["connected"], "browser")
        self.assertPage(widget.page)
        for action, expected in (("next", 25), ("previous", 0)):
            self.pager._on_message(widget, {"action": action, "request": action}, [])
            self.assertEqual(widget.page["start"], expected)
            self.assertEqual(widget.page["request"], action)
            self.assertPage(widget.page)
        self.pager._on_message(widget, {"action": "size", "value": 100}, [])
        self.assertEqual(widget.page["page_size"], 100)
        self.assertEqual(widget.page["end"], 100)
        self.assertPage(widget.page)
        self.pager._on_message(widget, {"action": "dispose", "view": "browser"}, [])
        self.assertIsNone(self.pager._widget)
        self.assertFalse(self.pager._cache)

    def test_nested_sums_keep_subexpression_navigation(self):
        value = TensorExpression(
            (E("x+1") ** 100).expand() * (E("y+1") ** 100).expand()
        )
        original = value.to_expression()
        pager = value.paged()
        self.addCleanup(pager.close)
        parent = pager._page()
        self.assertPage(parent)
        target = next(target for target, _ in parent["holes"] if target != 0)
        child = pager._action({"action": "open", "value": target})
        self.assertPage(child)
        self.assertTrue(child["breadcrumbs"])
        self.assertPage(pager._action({"action": "next"}))
        self.assertEqual(pager._action({"action": "back"})["html"], parent["html"])
        self.assertEqual(value.to_expression(), original)

    def test_oversized_output_becomes_an_inspectable_placeholder(self):
        import _spenso_paging

        value = TensorExpression(E("wasm_paging::f(a+b+c)"))
        pager = value.paged()
        self.addCleanup(pager.close)
        with patch.object(
            _spenso_paging, "compile_typst", return_value=b"x" * (256 * 1024 + 1)
        ):
            page = pager._page()
        self.assertIn("exceeds the display budget", page["html"])
        self.assertTrue(page["holes"])
        self.assertEqual(page["end"], 1)
        self.assertEqual(value.to_expression(), E("wasm_paging::f(a+b+c)"))

    def test_renderer_errors_stay_in_the_bounded_preview(self):
        import _spenso_paging

        with patch.object(
            _spenso_paging, "compile_typst", side_effect=RuntimeError("bad notation")
        ):
            page = self.pager._page()
        self.assertIn("Page rendering failed: bad notation", page["html"])
        self.assertEqual(self.value.to_expression(), self.original)


if __name__ == "__main__":
    unittest.main()
