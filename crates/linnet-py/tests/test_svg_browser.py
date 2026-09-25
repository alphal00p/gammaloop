"""Exercise the shared graph camera and inspector in Chromium and WebKit."""

import importlib.util
import unittest

import linnet as lp


@unittest.skipUnless(
    importlib.util.find_spec("playwright"), "Playwright is not installed"
)
class SvgBrowserTests(unittest.TestCase):
    def test_navigation_inspection_and_independent_views(self):
        from playwright.sync_api import sync_playwright

        left, right = lp.node("left"), lp.node("right")
        graph = lp.build(left, right, lp.edge(lp.source(left), "p", lp.sink(right)))
        drawing = graph.to_svg()
        document = (
            '<!doctype html><meta charset="utf-8">'
            f'<div id="first">{drawing}</div><div id="second">{drawing}</div>'
        )
        with sync_playwright() as playwright:
            for engine in (playwright.chromium, playwright.webkit):
                with self.subTest(browser=engine.name):
                    browser = engine.launch()
                    try:
                        page = browser.new_page(
                            viewport={"width": 1000, "height": 1000}
                        )
                        errors = []
                        page.on(
                            "pageerror",
                            lambda error, errors=errors: errors.append(str(error)),
                        )
                        page.set_content(document)
                        first = page.locator("#first")
                        second = page.locator("#second")
                        camera = first.locator(".linnet-viewport")
                        node = first.locator(
                            '[data-linnet-kind="node"][data-linnet-id="0"]'
                        )
                        other = first.locator(
                            '[data-linnet-kind="node"][data-linnet-id="1"]'
                        )
                        panel = first.locator(".linnet-inspector")
                        node.hover()
                        self.assertIn("Node 0", panel.inner_text())
                        self.assertGreater(
                            panel.bounding_box()["x"],
                            camera.bounding_box()["x"] + camera.bounding_box()["width"],
                        )
                        node.click()
                        other.hover()
                        self.assertIn("Node 0", panel.inner_text())
                        self.assertEqual(
                            first.locator(".linnet-toolbar button").count(), 0
                        )
                        initial = camera.get_attribute("viewBox")
                        camera.focus()
                        page.keyboard.press("+")
                        self.assertNotEqual(camera.get_attribute("viewBox"), initial)
                        self.assertEqual(first.locator("output").inner_text(), "105%")
                        self.assertEqual(second.locator("output").inner_text(), "100%")
                        self.assertIn("Node 0", panel.inner_text())
                        page.keyboard.press("0")
                        self.assertEqual(camera.get_attribute("viewBox"), initial)

                        # A drag beginning on a node pans without changing the selection.
                        bounds = node.bounding_box()
                        x = bounds["x"] + bounds["width"] / 2
                        y = bounds["y"] + bounds["height"] / 2
                        page.mouse.move(x, y)
                        page.keyboard.down("Shift")
                        page.mouse.down()
                        page.mouse.move(x + 50, y + 30, steps=6)
                        page.mouse.up()
                        page.keyboard.up("Shift")
                        self.assertNotEqual(camera.get_attribute("viewBox"), initial)
                        self.assertEqual(
                            first.locator("svg[data-linnet-interactive]").evaluate(
                                "svg => svg.linnetSelection.nodes"
                            ),
                            [],
                        )
                        camera.focus()
                        page.keyboard.press("0")
                        self.assertEqual(camera.get_attribute("viewBox"), initial)
                        page.keyboard.press("+")
                        self.assertNotEqual(camera.get_attribute("viewBox"), initial)
                        page.keyboard.press("0")
                        first.get_by_role("button", name="Close graph details").click()
                        other.hover()
                        self.assertIn("Node 1", panel.inner_text())
                        other.click(modifiers=["Shift"])
                        self.assertTrue(
                            first.locator(".linnet-inspector-construction").is_visible()
                        )
                        self.assertTrue(
                            second.locator(".linnet-inspector-empty").is_visible()
                        )

                        # Wheel zoom is opt-in so ordinary notebook scrolling still works.
                        node.hover()
                        page.mouse.wheel(0, -40)
                        self.assertEqual(camera.get_attribute("viewBox"), initial)
                        node.hover()
                        page.keyboard.down("Control")
                        page.mouse.wheel(0, -40)
                        page.keyboard.up("Control")
                        page.wait_for_function(
                            "initial => document.querySelector('.linnet-viewport').getAttribute('viewBox') !== initial",
                            arg=initial,
                        )
                        self.assertIn("Node 1", panel.inner_text())
                        for theme in ("light", "dark"):
                            page.emulate_media(color_scheme=theme)
                            page.set_viewport_size({"width": 320, "height": 1000})
                            page.wait_for_function(
                                """() => {
                                    const graph = document.querySelector('.linnet-viewport').getBoundingClientRect();
                                    const panel = document.querySelector('.linnet-inspector').getBoundingClientRect();
                                    return panel.top >= graph.bottom;
                                }"""
                            )
                            self.assertTrue(
                                page.evaluate(
                                    "document.documentElement.scrollWidth <= innerWidth"
                                )
                            )
                        self.assertEqual(errors, [])
                    finally:
                        browser.close()


if __name__ == "__main__":
    unittest.main()
