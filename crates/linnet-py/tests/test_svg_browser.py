"""Exercise the shared graph camera and inspector in Chromium and WebKit."""

import importlib.util
import unittest

from symbolica.community import graph as lp


@unittest.skipUnless(
    importlib.util.find_spec("playwright"), "Playwright is not installed"
)
class SvgBrowserTests(unittest.TestCase):
    def test_physical_scale_is_shared_and_fitted_drawings_are_not_cropped(self):
        from playwright.sync_api import sync_playwright

        sizes = [(120, 80), (360, 80), (120, 600), (900, 80)]
        drawings = []
        for width, height in sizes:
            source = (
                f"#set page(width: {width}pt, height: {height}pt, margin: 0pt)\n"
                '#link("#linnet-node-0")[#text(size: 12pt)[Same label]]'
            )
            drawings.append(
                lp.DiagramRender.from_sources({"main.typ": source.encode()}).to_html()
            )
        document = "<style>figure { margin: 0; }</style>" + "".join(
            f'<div id="graph-{index}" style="width:520px">{drawing}</div>'
            for index, drawing in enumerate(drawings)
        )
        with sync_playwright() as playwright:
            for engine in (playwright.chromium, playwright.webkit):
                with self.subTest(browser=engine.name):
                    browser = engine.launch()
                    try:
                        page = browser.new_page(
                            viewport={"width": 1100, "height": 1000}
                        )
                        page.set_content(document)
                        labels = []
                        for index, (width, height) in enumerate(sizes):
                            owner = page.locator(f"#graph-{index}")
                            camera = owner.locator(".linnet-viewport")
                            scale = camera.evaluate("s => s.getScreenCTM().a")
                            expected = min(4 / 3, 520 / width)
                            self.assertAlmostEqual(scale, expected, places=5)
                            self.assertEqual(
                                owner.locator("output").inner_text(),
                                f"{round(100 * expected / (4 / 3))}%",
                            )
                            x, y, w, h = map(
                                float, camera.get_attribute("viewBox").split()
                            )
                            self.assertLessEqual(x, 1e-6)
                            self.assertLessEqual(y, 1e-6)
                            self.assertGreaterEqual(x + w, width - 1e-6)
                            self.assertGreaterEqual(y + h, height - 1e-6)
                            labels.append(
                                owner.locator(
                                    '[data-linnet-kind="node"]'
                                ).bounding_box()
                            )
                        # Identical authored text has identical displayed dimensions,
                        # despite different drawing widths and a tall viewport.
                        for label in labels[1:3]:
                            self.assertAlmostEqual(
                                label["height"], labels[0]["height"], places=3
                            )
                            self.assertAlmostEqual(
                                label["width"], labels[0]["width"], places=3
                            )
                        self.assertGreater(
                            page.locator("#graph-2 > figure > svg").bounding_box()[
                                "height"
                            ],
                            800,
                        )

                        first = page.locator("#graph-0")
                        camera = first.locator(".linnet-viewport")
                        camera.focus()
                        page.keyboard.press("+")
                        first.evaluate("e => e.style.width = '300px'")
                        page.wait_for_function(
                            "document.querySelector('#graph-0 svg').viewBox.baseVal.width === 300"
                        )
                        self.assertAlmostEqual(
                            camera.evaluate("s => s.getScreenCTM().a"), 1.4, places=5
                        )
                        self.assertEqual(first.locator("output").inner_text(), "105%")
                        page.keyboard.press("0")
                        self.assertEqual(first.locator("output").inner_text(), "100%")
                    finally:
                        browser.close()

    def test_neutral_palette_follows_media_and_explicit_theme(self):
        from playwright.sync_api import sync_playwright

        left, right = lp.node("left"), lp.node("right")
        graph = lp.build(left, right, lp.edge(lp.source(left), "p", lp.sink(right)))
        with sync_playwright() as playwright:
            for engine in (playwright.chromium, playwright.webkit):
                with self.subTest(browser=engine.name):
                    browser = engine.launch()
                    try:
                        page = browser.new_page()
                        page.set_content(graph.render().to_html())
                        paper = page.locator(
                            'svg[data-linnet-interactive] [fill="#ffffff"]'
                        ).first
                        ink = page.locator(
                            'svg[data-linnet-interactive] [fill="#000000"]'
                        ).first
                        self.assertGreater(paper.count(), 0)
                        self.assertGreater(ink.count(), 0)
                        for theme, background, foreground in (
                            ("light", "rgb(255, 255, 255)", "rgb(0, 0, 0)"),
                            ("dark", "rgb(28, 32, 37)", "rgb(230, 235, 241)"),
                        ):
                            page.emulate_media(color_scheme=theme)
                            self.assertEqual(
                                paper.evaluate("e => getComputedStyle(e).fill"),
                                background,
                            )
                            self.assertEqual(
                                ink.evaluate("e => getComputedStyle(e).fill"),
                                foreground,
                            )
                        page.locator("html").evaluate("e => e.dataset.theme = 'light'")
                        self.assertEqual(
                            paper.evaluate("e => getComputedStyle(e).fill"),
                            "rgb(255, 255, 255)",
                        )
                        self.assertEqual(
                            ink.evaluate("e => getComputedStyle(e).fill"),
                            "rgb(0, 0, 0)",
                        )
                    finally:
                        browser.close()

    def test_collection_viewport_keeps_hover_and_pinned_details_inside_graph(self):
        from playwright.sync_api import sync_playwright

        left, right = lp.node("left"), lp.node("right")
        graph = lp.build(left, right, lp.edge(lp.source(left), "p", lp.sink(right)))
        drawing = graph.to_svg().replace(
            'data-linnet-interactive="true"',
            'data-linnet-interactive="true" data-linnet-viewport-height="360"',
            1,
        )
        with sync_playwright() as playwright:
            for engine in (playwright.chromium, playwright.webkit):
                with self.subTest(browser=engine.name):
                    browser = engine.launch()
                    try:
                        page = browser.new_page(viewport={"width": 420, "height": 800})
                        page.set_content('<meta charset="utf-8">' + drawing)
                        svg = page.locator("svg[data-linnet-interactive]")
                        target = page.locator('[data-linnet-kind="node"]').first
                        self.assertEqual(svg.bounding_box()["height"], 400)
                        target.dispatch_event("pointerover")
                        self.assertEqual(svg.bounding_box()["height"], 400)
                        panel = page.locator(".linnet-inspector")
                        self.assertTrue(panel.is_visible())
                        target.dispatch_event("click")
                        target.dispatch_event("pointerout")
                        self.assertIn("Pinned", panel.inner_text())
                        self.assertEqual(svg.bounding_box()["height"], 400)
                        page.get_by_role("button", name="Close graph details").click()
                        self.assertFalse(panel.is_visible())
                    finally:
                        browser.close()

    def test_preview_only_while_hovered_or_focused_unless_pinned(self):
        from playwright.sync_api import sync_playwright

        left, right = lp.node("left"), lp.node("right")
        graph = lp.build(left, right, lp.edge(lp.source(left), "p", lp.sink(right)))
        with sync_playwright() as playwright:
            for engine in (playwright.chromium, playwright.webkit):
                with self.subTest(browser=engine.name):
                    browser = engine.launch()
                    try:
                        page = browser.new_page(viewport={"width": 1000, "height": 800})
                        page.set_content(
                            '<!doctype html><meta charset="utf-8">' + graph.to_svg()
                        )
                        panel = page.locator(".linnet-inspector")
                        self.assertFalse(panel.is_visible())
                        for kind in ("node", "edge", "halfedge"):
                            target = page.locator(f'[data-linnet-kind="{kind}"]').first
                            target.dispatch_event("pointerover")
                            self.assertTrue(panel.is_visible())
                            target.dispatch_event("pointerout")
                            self.assertFalse(panel.is_visible())
                        node = page.locator('[data-linnet-kind="node"]').first
                        node.focus()
                        self.assertTrue(panel.is_visible())
                        page.locator(".linnet-viewport").focus()
                        self.assertFalse(panel.is_visible())
                        node.click()
                        page.mouse.move(990, 790)
                        self.assertTrue(panel.is_visible())
                        page.get_by_role("button", name="Close graph details").click()
                        self.assertFalse(panel.is_visible())
                        node.hover()
                        self.assertTrue(panel.is_visible())
                        page.mouse.move(990, 790)
                        self.assertFalse(panel.is_visible())
                    finally:
                        browser.close()

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
                        self.assertFalse(
                            second.locator(".linnet-inspector").is_visible()
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
