"""Exercise the shared graph camera and inspector in Chromium and WebKit."""

import importlib.util
import math
import os
import struct
import unittest
import xml.etree.ElementTree as ET

from rendering_metrics import external_metrics_snapshot
from symbolica.community import graph as lp


def pinned_drawing(x, y):
    """Render a labelled origin joined to a node pinned at (x, y)."""

    first = lp.node("first", label=lp.TextLabel("Same label"))
    second = lp.node("second")
    first.drawing.placement = {"x": 0.0, "y": 0.0, "mode": lp.Placement.Pin}
    second.drawing.placement = {"x": x, "y": y, "mode": lp.Placement.Pin}
    graph = lp.build(first, second, lp.edge(lp.source(first), "line", lp.sink(second)))
    return external_metrics_snapshot(graph.render())


@unittest.skipUnless(
    importlib.util.find_spec("playwright"), "Playwright is not installed"
)
class SvgBrowserTests(unittest.TestCase):
    def test_readable_scale_and_fitted_drawings_are_not_cropped(self):
        from playwright.sync_api import sync_playwright

        renders = [pinned_drawing(x, y) for x, y in ((0, 2), (22, 2), (0, 42), (64, 2))]
        sizes = []
        for render in renders:
            root = ET.fromstring(render.to_svg())
            sizes.append(
                tuple(
                    float(root.get(side).removesuffix("pt"))
                    for side in ("width", "height")
                )
            )
        self.assertLess(sizes[1][0], 390)
        self.assertGreater(sizes[2][1], 600)
        self.assertGreater(sizes[3][0], 390)
        # Standalone SVGs use pt (4/3 CSS px); the readable preview is bounded
        # by 450px width, 300px height, and a total enlargement of 3.
        natural_scales = [
            max(4 / 3, min(3, 450 / width, 300 / height))
            for width, height in sizes
        ]
        # The label comparisons use unfitted previews, and reset must still
        # reach 100% after the first preview's column narrows to 300px.
        for (width, _), natural in zip(sizes[:3], natural_scales[:3]):
            self.assertLessEqual(width * natural, 520)
        self.assertLessEqual(sizes[0][0] * natural_scales[0], 300)
        drawings = [render.to_html() for render in renders]
        document = "<style>figure { margin: 0; }</style>" + "".join(
            f'<div id="graph-{index}" style="width:520px">{drawing}</div>'
            for index, drawing in enumerate(drawings)
        )

        def outer_scale(owner, column_width, expected_height):
            geometry = owner.locator(":scope > figure > svg").evaluate(
                """svg => {
                    const m = svg.getScreenCTM();
                    const box = svg.getBoundingClientRect();
                    const view = svg.viewBox.baseVal;
                    return {a:m.a, b:m.b, c:m.c, d:m.d,
                        width:box.width, height:box.height,
                        view_width:view.width, view_height:view.height};
                }"""
            )
            self.assertGreater(geometry["a"], 0)
            self.assertEqual(geometry["a"], geometry["d"])
            self.assertEqual(geometry["b"], 0)
            self.assertEqual(geometry["c"], 0)
            self.assertGreater(geometry["view_height"], 0)
            self.assertEqual(geometry["view_width"], column_width)
            # Five binary32 roundings bound authored dimensions, outer
            # viewBox and screen matrix; one layout unit is 1/64 of a CSS
            # pixel.
            unit_roundoff = 2**-24
            representation = expected_height * (
                5 * unit_roundoff / (1 - 5 * unit_roundoff)
            )
            self.assertAlmostEqual(geometry["width"], column_width, delta=1 / 64)
            self.assertAlmostEqual(
                geometry["height"],
                expected_height,
                delta=1 / 64 + representation,
            )
            self.assertAlmostEqual(
                geometry["view_height"], expected_height, delta=representation
            )
            # Default meet uses the smaller viewport/viewBox ratio. Check it
            # independently so normalization cannot conceal CSS scaling.
            meet = min(
                geometry["width"] / geometry["view_width"],
                geometry["height"] / geometry["view_height"],
            )
            self.assertAlmostEqual(geometry["a"], meet, delta=2**-23)
            self.assertLessEqual(
                abs(geometry["a"] - 1) * geometry["view_height"],
                1 / 64 + representation,
            )
            return geometry["a"]

        with sync_playwright() as playwright:
            selected = os.environ.get("LINNET_TEST_BROWSER_ENGINES", "chromium,webkit")
            for name in selected.split(","):
                engine = getattr(playwright, name.strip())
                with self.subTest(browser=engine.name):
                    executable = (
                        os.environ.get("LINNET_TEST_CHROMIUM")
                        if name == "chromium"
                        else None
                    )
                    browser = engine.launch(
                        executable_path=executable,
                        args=["--no-sandbox", "--disable-dev-shm-usage"]
                        if name == "chromium"
                        else [],
                    )
                    try:
                        page = browser.new_page(
                            viewport={"width": 1100, "height": 1000}
                        )
                        page.set_content(document)
                        labels = []
                        for index, (width, height) in enumerate(sizes):
                            owner = page.locator(f"#graph-{index}")
                            camera = owner.locator(".linnet-viewport")
                            natural = natural_scales[index]
                            expected = min(natural, 520 / width)
                            outer = outer_scale(owner, 520, expected * height)
                            scale = camera.evaluate("s => s.getScreenCTM().a") / outer
                            self.assertAlmostEqual(scale, expected, delta=1e-4)
                            self.assertEqual(
                                owner.locator("output").inner_text(),
                                f"{math.floor(100 * expected / natural + 0.5)}%",
                            )
                            x, y, w, h = map(
                                float, camera.get_attribute("viewBox").split()
                            )
                            # The camera's box comes from float32 SVGRect values.
                            width32, height32 = struct.unpack(
                                "!2f", struct.pack("!2f", width, height)
                            )
                            self.assertLessEqual(x, 1e-6)
                            self.assertLessEqual(y, 1e-6)
                            self.assertGreaterEqual(x + w, width32 - 1e-6)
                            self.assertGreaterEqual(y + h, height32 - 1e-6)
                            label = owner.locator(
                                '[data-linnet-kind="node"][data-linnet-id="0"]'
                            ).bounding_box()
                            labels.append(
                                {side: label[side] / outer for side in ("width", "height")}
                            )
                        # The same node-0 label target has fixed authored dimensions;
                        # its screen dimensions grow with the natural scale.
                        for index, label in enumerate(labels[1:3], start=1):
                            ratio = natural_scales[0] / natural_scales[index]
                            self.assertAlmostEqual(
                                label["height"] * ratio, labels[0]["height"], places=2
                            )
                            self.assertAlmostEqual(
                                label["width"] * ratio, labels[0]["width"], places=2
                            )
                        # Natural height has no toolbar slack; allow one CSS
                        # layout unit (1/64px) for the outer SVG's snapped height.
                        self.assertAlmostEqual(
                            page.locator("#graph-2 > figure > svg").bounding_box()[
                                "height"
                            ],
                            4 / 3 * sizes[2][1],
                            delta=1 / 64,
                        )

                        first = page.locator("#graph-0")
                        camera = first.locator(".linnet-viewport")
                        camera.focus()
                        page.keyboard.press("+")
                        first.evaluate("e => e.style.width = '300px'")
                        page.wait_for_function(
                            "document.querySelector('#graph-0 svg').viewBox.baseVal.width === 300"
                        )
                        # Resize can change the outer meet transform; read it again.
                        outer = outer_scale(first, 300, natural_scales[0] * sizes[0][1])
                        self.assertAlmostEqual(
                            camera.evaluate("s => s.getScreenCTM().a") / outer,
                            1.05 * natural_scales[0],
                            delta=1e-4,
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
                        page.set_content(
                            external_metrics_snapshot(graph.render()).to_html()
                        )
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
        drawing = (
            external_metrics_snapshot(graph.render())
            .to_svg()
            .replace(
                'data-linnet-interactive="true"',
                'data-linnet-interactive="true" data-linnet-viewport-height="360"',
                1,
            )
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
                        self.assertEqual(svg.bounding_box()["height"], 360)
                        target.dispatch_event("pointerover")
                        self.assertEqual(svg.bounding_box()["height"], 360)
                        panel = page.locator(".linnet-inspector")
                        self.assertTrue(panel.is_visible())
                        target.dispatch_event("click")
                        target.dispatch_event("pointerout")
                        self.assertIn("Pinned", panel.inner_text())
                        self.assertEqual(svg.bounding_box()["height"], 360)
                        bounds = panel.bounding_box()
                        self.assertGreaterEqual(bounds["x"], 8)
                        self.assertGreaterEqual(bounds["y"], 8)
                        self.assertLessEqual(bounds["x"] + bounds["width"], 412)
                        self.assertLessEqual(bounds["y"] + bounds["height"], 792)
                        help = svg.locator(".linnet-navigation-help")
                        help.hover()
                        self.assertTrue(
                            svg.locator(".linnet-navigation-hint").is_visible()
                        )
                        page.mouse.move(0, 0)
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
                            '<!doctype html><meta charset="utf-8">'
                            + external_metrics_snapshot(graph.render()).to_svg()
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
        drawing = external_metrics_snapshot(graph.render()).to_svg()
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
                        panel = page.locator(".linnet-inspector")
                        node.hover()
                        self.assertIn("Node 0", panel.inner_text())
                        self.assertGreater(
                            panel.bounding_box()["y"],
                            node.bounding_box()["y"],
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
                        page.get_by_role("button", name="Close graph details").click()
                        other.hover()
                        self.assertIn("Node 1", panel.inner_text())
                        other.click(modifiers=["Shift"])
                        self.assertTrue(
                            page.locator(".linnet-inspector-construction").is_visible()
                        )
                        self.assertEqual(page.locator(".linnet-inspector").count(), 1)
                        self.assertEqual(
                            second.locator("svg[data-linnet-interactive]").evaluate(
                                "svg => svg.linnetSelection.nodes"
                            ),
                            [],
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
                                    const host = document.querySelector('.linnet-inspector-popup');
                                    const panel = host.getBoundingClientRect();
                                    return panel.left >= 8 && panel.right <= innerWidth - 8
                                      && panel.top >= 8 && panel.bottom <= innerHeight - 8;
                                }"""
                            )
                            self.assertTrue(
                                page.evaluate(
                                    "document.documentElement.scrollWidth <= innerWidth"
                                )
                            )
                            self.assertEqual(
                                panel.evaluate(
                                    "e => getComputedStyle(e).backgroundColor"
                                ),
                                "rgb(37, 33, 43)"
                                if theme == "dark"
                                else "rgb(255, 253, 249)",
                            )
                        self.assertEqual(errors, [])
                    finally:
                        browser.close()


if __name__ == "__main__":
    unittest.main()
