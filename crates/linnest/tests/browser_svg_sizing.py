"""Check notebook iframe sizing with Playwright; --chromium selects the browser."""

import argparse
import json
from pathlib import Path

from playwright.sync_api import sync_playwright


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromium")
    args = parser.parse_args()
    source = Path(__file__).parents[1] / "src/svg"
    drawing = f"""
        <style>
          body {{ margin: 0; }}
          figure {{ margin: 8px 0; }}
          figcaption {{ font: 14px/20px sans-serif; }}
        </style>
        <figure>
          <div style="max-width:100%;overflow-x:auto">
            <svg data-linnet-interactive viewBox="0 0 1000 700">
              <rect width="1000" height="700" fill="lightblue"/>
              <style>{(source / "interactive.css").read_text()}</style>
            </svg>
          </div>
          <figcaption>Graph sizing regression</figcaption>
        </figure>
        <script>{(source / "interactive.js").read_text()}</script>
    """
    with sync_playwright() as driver:
        # Headless Playwright hides scrollbars by default, masking the width
        # change that triggers this regression with desktop scrollbars.
        browser = driver.chromium.launch(
            executable_path=args.chromium, ignore_default_args=["--hide-scrollbars"]
        )
        try:
            page = browser.new_page(viewport={"width": 1200, "height": 800})
            errors = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.set_content("""
                <style>
                  body { margin: 0; }
                  #output { width: 1000px; max-height: 580px; overflow: auto; }
                  iframe { width: 100%; border: 0; }
                </style>
                <div id="output"><iframe></iframe></div>
            """)
            page.locator("iframe").evaluate(
                "(frame, content) => { frame.srcdoc = content; }", drawing
            )
            page.frame_locator("iframe").locator(".linnet-viewport").wait_for()
            # Cross the scrollbar threshold in both directions at wide and
            # narrow notebook column widths.
            for width, limit in (
                (1000, 580),
                (1000, 550),
                (1000, 580),
                (1000, 570),
                (1000, 560),
                (610, 300),
                (590, 450),
                (1000, 580),
            ):
                samples = page.evaluate(
                    """async ([width, limit]) => {
                        const output = document.querySelector('#output');
                        const frame = document.querySelector('iframe');
                        output.style.width = width + 'px';
                        output.style.maxHeight = limit + 'px';
                        const samples = [];
                        for (let i = 0; i < 30; i++) {
                            await new Promise(requestAnimationFrame);
                            samples.push([output.clientWidth, output.scrollHeight,
                                frame.getBoundingClientRect().height]);
                        }
                        return samples;
                    }""",
                    [width, limit],
                )
                assert len({tuple(sample) for sample in samples[-10:]}) == 1, (
                    width,
                    limit,
                    samples,
                )
                print(f"width={width}, limit={limit}: stable at {samples[-1]}")

            # Leave exactly enough room for the iframe. Inline baselines must
            # not add a few pixels of overflow outside it or below the SVG.
            measurements = page.evaluate("""() => {
                const output = document.querySelector('#output');
                const frame = document.querySelector('iframe');
                const doc = frame.contentDocument;
                const svg = doc.querySelector('svg');
                const caption = doc.querySelector('figcaption');
                output.style.maxHeight = frame.style.height;
                return {
                    overflow: output.scrollHeight - output.clientHeight,
                    baseline: svg.parentElement.getBoundingClientRect().height
                        - svg.getBoundingClientRect().height,
                    slack: frame.getBoundingClientRect().height
                        - caption.getBoundingClientRect().bottom,
                };
            }""")
            assert measurements["overflow"] == 0, measurements
            assert measurements["baseline"] == 0, measurements
            assert 0 <= measurements["slack"] < 1, measurements
            assert not errors, errors
            print("Exact fit:", json.dumps(measurements))

            # Small standalone drawings get a readable baseline. Large and
            # fixed-height drawings retain their scale and are never cropped.
            # Include a tall graph exceeding the old 520px height cap.
            page.set_content(
                f"<style>{(source / 'interactive.css').read_text()}</style>"
                + "".join(
                    f'<div style="width:550px"><svg id="{name}" '
                    f'data-linnet-interactive viewBox="0 0 {width} {height}" '
                    f'width="{width}pt" height="{height}pt" '
                    f'data-linnet-viewport-height="{fixed_height}">'
                    f'<rect width="{width}" height="{height}" fill="lightblue"/>'
                    '<circle cx="20" cy="20" r="5"/></svg></div>'
                    for name, width, height, fixed_height in (
                        ("tall", 100, 600, 0),
                        ("wide", 400, 100, 0),
                        ("small", 50, 50, 0),
                        ("fixed", 100, 600, 200),
                    )
                )
            )
            page.add_script_tag(content=(source / "interactive.js").read_text())
            measure = """svg => {
                const viewport = svg.querySelector('.linnet-viewport');
                const drawing = viewport.querySelector('rect');
                const bounds = drawing.getBoundingClientRect();
                const outer = svg.getBoundingClientRect();
                return {
                    scale: drawing.getScreenCTM().a,
                    label: svg.querySelector('output').textContent,
                    diameter: viewport.querySelector('circle').getBoundingClientRect().width,
                    contained: bounds.left >= outer.left - .01
                        && bounds.right <= outer.right + .01
                        && bounds.top >= outer.top - .01
                        && bounds.bottom <= outer.bottom + .01,
                };
            }"""
            # Allow browser subpixel rounding of the outer SVG dimensions.
            for name, expected in (("tall", 4 / 3), ("wide", 4 / 3), ("small", 3)):
                result = page.locator(f"#{name}").evaluate(measure)
                assert abs(result["scale"] - expected) < 1e-3, result
                assert abs(result["diameter"] - 10 * expected) < 0.02, result
                assert result["label"] == "100%", result
                assert result["contained"], result
            fixed = page.locator("#fixed").evaluate(measure)
            assert abs(fixed["scale"] - 1 / 3) < 1e-3, fixed
            assert fixed["label"] == "25%" and fixed["contained"], fixed

            # Shrink to notebook width, recover native scale on expansion,
            # and keep the displayed percentage truthful after zoom/reset.
            wide = page.locator("#wide")
            for width, label in ((120, "23%"), (550, "100%")):
                wide.evaluate(
                    "(svg, width) => svg.parentElement.style.width = width + 'px'",
                    width,
                )
                page.wait_for_function(
                    "label => document.querySelector('#wide output').textContent === label",
                    arg=label,
                )
                result = wide.evaluate(measure)
                assert result["contained"], result
                assert abs(result["scale"] - min(width / 400, 4 / 3)) < 1e-3, result
            viewport = wide.locator(".linnet-viewport")
            viewport.press("+")
            assert wide.locator("output").inner_text() == "105%"
            assert viewport.evaluate("e=>getComputedStyle(e).outlineStyle") == "none"
            assert (
                wide.locator(".linnet-pan-surface").evaluate(
                    "e=>getComputedStyle(e).strokeOpacity"
                )
                == "0.5"
            )
            viewport.press("0")
            assert wide.locator("output").inner_text() == "100%"
            assert wide.evaluate(measure)["contained"]
            assert not errors, errors
            print(
                "Readable baseline, width fitting, fixed-height fitting, and zoom/reset passed"
            )
            # A small notebook diagram should be readable and centered above
            # its caption, without blank space reserved for hidden details.
            page.set_content(
                f"<style>{(source / 'interactive.css').read_text()}</style>"
                + """
                <figure id="compact" style="width:800px;margin:0;text-align:center">
                  <svg data-linnet-interactive viewBox="0 0 150 75"
                       width="150pt" height="75pt">
                    <rect width="150" height="75" fill="lightblue"/>
                    <a data-linnet-kind="node" data-linnet-id="0"
                       data-linnet-detail='{"edges":[0]}'>
                      <text x="75" y="40" font-size="9" text-anchor="middle">g</text>
                    </a>
                  </svg>
                  <figcaption>One-loop diagram</figcaption>
                </figure>
                """
            )
            page.add_script_tag(content=(source / "interactive.js").read_text())
            compact = page.locator("#compact")
            svg = compact.locator("svg[data-linnet-interactive]")
            camera = compact.locator(".linnet-viewport")
            label = camera.locator("text")
            assert svg.locator("output").inner_text() == "100%"
            assert 220 < svg.bounding_box()["height"] < 230
            assert abs(label.evaluate("e => e.getScreenCTM().a") * 9 - 27) < 0.01
            for width, expected in ((320, "71%"), (220, "49%"), (800, "100%")):
                compact.evaluate("(e, width) => e.style.width = width + 'px'", width)
                page.wait_for_function(
                    "expected => document.querySelector('#compact output').textContent === expected",
                    arg=expected,
                )
                result = svg.evaluate("""svg => {
                    const art = svg.querySelector('.linnet-viewport > rect').getBoundingClientRect();
                    const caption = svg.parentElement.querySelector('figcaption').getBoundingClientRect();
                    const outer = svg.getBoundingClientRect();
                    const viewport = svg.querySelector('.linnet-viewport');
                    const controls = svg.querySelector('.linnet-toolbar').getBoundingClientRect();
                    return {
                        centered: Math.abs((art.left + art.right - caption.left - caption.right) / 2) < .1,
                        contained: art.left >= outer.left - .1 && art.right <= outer.right + .1
                            && art.top >= outer.top - .1 && art.bottom <= outer.bottom + .1,
                        noToolbarRow: viewport.y.baseVal.value === 0
                            && Math.abs(viewport.height.baseVal.value - outer.height) < .1,
                        compactControls: controls.width < 84 && controls.height <= 28
                            && controls.right <= outer.right && controls.top >= outer.top
                            && controls.left > outer.right - 88,
                        captionGap: caption.top - outer.bottom,
                    };
                }""")
                assert result["centered"] and result["contained"], result
                assert result["noToolbarRow"] and result["compactControls"], result
                assert 0 <= result["captionGap"] <= 4, result
                # The hint stays legible at narrow widths without creating a
                # toolbar row or adding output height, and remains accessible
                # by keyboard as well as hover.
                help = compact.locator(".linnet-navigation-help")
                hint = compact.locator(".linnet-navigation-hint")
                height = svg.bounding_box()["height"]
                help.focus()
                assert hint.is_visible()
                assert "0 to fit" in help.get_attribute("aria-label")
                hint_box, svg_box = hint.bounding_box(), svg.bounding_box()
                assert hint_box["x"] >= svg_box["x"] - 0.1
                assert (
                    hint_box["x"] + hint_box["width"]
                    <= svg_box["x"] + svg_box["width"] + 0.1
                )
                assert svg.bounding_box()["height"] == height
                initial_camera = camera.get_attribute("viewBox")
                help.dispatch_event("wheel", {"deltaY": -100, "ctrlKey": True})
                assert camera.get_attribute("viewBox") == initial_camera
                camera.focus()
                assert not hint.is_visible()
            before = label.bounding_box()
            height = svg.bounding_box()["height"]
            target = compact.locator('[data-linnet-kind="node"]')
            target.hover()
            panel = page.locator(".linnet-inspector")
            assert panel.is_visible()
            target_box, popup_box = target.bounding_box(), panel.bounding_box()
            assert popup_box["x"] >= target_box["x"] + target_box["width"] / 2
            assert popup_box["y"] >= target_box["y"] + target_box["height"] / 2
            assert svg.bounding_box()["height"] == height
            assert label.bounding_box() == before
            target.click()
            page.mouse.move(1000, 700)
            assert panel.is_visible() and label.bounding_box() == before
            page.get_by_role("button", name="Close graph details").click()
            assert not panel.is_visible() and svg.bounding_box()["height"] == height
            camera.press("+")
            assert svg.locator("output").inner_text() == "105%"
            assert abs(label.evaluate("e => e.getScreenCTM().a") - 3.15) < 1e-3
            camera.press("0")
            assert svg.locator("output").inner_text() == "100%"
            assert abs(label.evaluate("e => e.getScreenCTM().a") - 3) < 1e-3
            assert label.bounding_box() == before
            assert svg.bounding_box()["height"] == height
            # At the reported 223% zoom, clipping is confined to the canvas;
            # the output frame and caption still occupy their original bounds.
            bounds = svg.bounding_box()
            for _ in range(8):
                camera.dispatch_event(
                    "wheel",
                    {
                        "deltaY": -100,
                        "ctrlKey": True,
                        "clientX": bounds["x"] + bounds["width"] / 2,
                        "clientY": bounds["y"] + bounds["height"] / 2,
                    },
                )
            assert svg.locator("output").inner_text() == "223%"
            assert svg.bounding_box()["height"] == height
            assert camera.evaluate("e => getComputedStyle(e).overflow") == "hidden"
            assert camera.get_attribute("y") == "0"
            camera.press("0")
            assert label.bounding_box() == before
            assert not errors, errors
            print(
                "Readable compact preview, centered resizing, stable inspection, and baseline reset passed"
            )
            # Graph-only notebook frames shed the browser's default body inset,
            # while explicit page margins and the notebook's padding remain.
            for body_style, inset in (("", 0), ("body { margin: 12px; }", 12)):
                page.set_content(
                    '<div style="width:600px;padding:16px"><iframe style="width:100%;border:0"></iframe></div>'
                )
                graph_only = drawing.replace("body { margin: 0; }", body_style).replace(
                    "margin: 8px 0", "margin: 0"
                )
                page.locator("iframe").evaluate(
                    "(frame, content) => frame.srcdoc = content", graph_only
                )
                page.frame_locator("iframe").locator(".linnet-viewport").wait_for()
                page.wait_for_function("document.querySelector('iframe').style.height")
                state = page.evaluate("""() => {
                    const frame = document.querySelector('iframe');
                    const doc = frame.contentDocument;
                    const svg = doc.querySelector('svg');
                    return {
                        margin: frame.contentWindow.getComputedStyle(doc.body).margin,
                        canvas: svg.querySelector('.linnet-viewport').y.baseVal.value,
                        svg: svg.getBoundingClientRect().x,
                        hostPadding: getComputedStyle(frame.parentElement).padding,
                        slack: frame.getBoundingClientRect().height
                            - doc.querySelector('figcaption').getBoundingClientRect().bottom,
                    };
                }""")
                assert state["margin"] == f"{inset}px" and state["svg"] == inset, state
                assert state["canvas"] == 0 and state["hostPadding"] == "16px", state
                assert inset <= state["slack"] < inset + 1, state
            assert not errors, errors
            print(
                "Compact iframe margins and authored/notebook padding preservation passed"
            )
            # Pinned overlay details leave the floating controls and their
            # tooltip accessible without taking height from the graph canvas.
            for embedded in (False, True):
                detail = json.dumps({"properties": [["Details", "long text " * 100]]})
                graph = f"""
                    <svg id="overlay" data-linnet-interactive
                         viewBox="0 0 500 180" width="500" height="180"
                         data-linnet-viewport-height="{0 if embedded else 180}">
                      <rect width="500" height="180" fill="lightblue"/>
                      <a data-linnet-kind="node" data-linnet-id="0"
                         data-linnet-detail='{detail}'>
                        <text x="20" y="90">node</text>
                      </a>
                    </svg>
                """
                if embedded:
                    graph = f'<svg width="500" height="180">{graph}</svg>'
                page.set_content(
                    f"<style>{(source / 'interactive.css').read_text()}</style>"
                    f'<div style="width:500px">{graph}</div>'
                )
                page.add_script_tag(content=(source / "interactive.js").read_text())
                svg = page.locator("#overlay")
                canvas = svg.locator(".linnet-pan-surface")
                canvas_geometry = """e => [e.x.baseVal.value, e.y.baseVal.value,
                    e.width.baseVal.value, e.height.baseVal.value,
                    e.getScreenCTM().toString()]"""
                original_canvas = canvas.evaluate(canvas_geometry)
                original = svg.bounding_box()
                svg.locator('[data-linnet-kind="node"]').dispatch_event("click")
                panel = page.locator(".linnet-inspector")
                assert "Pinned" in panel.inner_text()
                bounds = panel.bounding_box()
                assert bounds["x"] >= 8 and bounds["y"] >= 8
                assert bounds["x"] + bounds["width"] <= 1200 - 8
                assert bounds["y"] + bounds["height"] <= 800 - 8
                page.set_viewport_size({"width": 1200, "height": 200})
                assert panel.evaluate("e => e.clientHeight < e.scrollHeight"), (
                    "Long details remain scrollable"
                )
                page.set_viewport_size({"width": 1200, "height": 800})
                assert svg.bounding_box() == original
                camera = svg.locator(".linnet-viewport")
                assert camera.get_attribute("y") == "0"
                assert camera.get_attribute("height") == "180"
                help = svg.locator(".linnet-navigation-help")
                help.hover()
                hint = svg.locator(".linnet-navigation-hint")
                assert hint.is_visible()
                assert hint.evaluate("""element => {
                    const bounds = element.getBoundingClientRect();
                    return element.contains(document.elementFromPoint(
                        bounds.x + bounds.width / 2, bounds.y + bounds.height / 2));
                }""")
                help.focus()
                assert hint.is_visible()
                camera.press("+")
                assert svg.locator("output").inner_text() == "105%"
                # Ignore focus strokes and zoomed ink in SVG bounding boxes.
                assert canvas.evaluate(canvas_geometry) == original_canvas
                page.mouse.move(0, 0)
                assert not hint.is_visible()
                page.get_by_role("button", name="Close graph details").click()
                assert not panel.is_visible()
            assert not errors, errors
            print("Pinned bounded and embedded inspector controls remain accessible")
            # A domain overlay uses the same precise native half-edge geometry
            # as hover, without changing selection or depending on camera zoom.
            page.set_content("""
                <svg id="native" data-linnet-interactive data-linnet-renderer="native"
                     viewBox="0 0 140 80" width="140pt" height="80pt"
                     data-linnet-viewport-height="180">
                  <circle data-linnet-node="0" cx="10" cy="40" r="4"/>
                  <circle data-linnet-node="1" cx="130" cy="40" r="4"/>
                  <a data-linnet-carrier data-linnet-kind="halfedge" data-linnet-id="0"
                     data-linnet-detail='{"edge":0,"half-edge":0,"node":0,"flow":"source"}'
                     style="display:none"><path d="M14 40C30 10 50 40 70 40"/></a>
                  <a data-linnet-carrier data-linnet-kind="halfedge" data-linnet-id="1"
                     data-linnet-detail='{"edge":0,"half-edge":1,"node":1,"flow":"sink"}'
                     style="display:none"><path d="M70 40C90 40 110 70 126 40"/></a>
                </svg>
                <svg id="authored" data-linnet-interactive viewBox="0 0 140 80" width="140pt"
                     data-linnet-viewport-height="180">
                  <g transform="translate(10 20)">
                    <a data-linnet-kind="halfedge" data-linnet-id="3"
                       data-linnet-detail='{"edge":2,"half-edge":3,"node":1,"flow":"source"}'>
                      <rect x="20" y="10" width="8" height="8"/>
                    </a>
                  </g>
                </svg>
            """)
            page.add_script_tag(content=(source / "interactive.js").read_text())
            shade = """svg => {
                const layer = svg.linnetShade({nodes: [0], half_edges: [0]});
                svg.querySelector('.linnet-viewport').append(layer);
                const path = layer.querySelector('path');
                const matrix = path.transform.baseVal.consolidate().matrix;
                const result = {
                    curve: path.getAttribute('d'),
                    cap: path.getAttribute('stroke-linecap'),
                    start: [path.getPointAtLength(0).x, path.getPointAtLength(0).y],
                    node: layer.querySelector('circle').dataset.linnetShadeNode,
                    halves: [...layer.querySelectorAll('path')].map(p=>p.dataset.linnetShadeHalfedge),
                    transform: [matrix.a, matrix.b, matrix.c, matrix.d, matrix.e, matrix.f],
                    selection: svg.linnetSelection,
                    edgeHalves: svg.linnetShade({edges: [0]}).querySelectorAll('path').length,
                };
                layer.remove();
                return result;
            }"""
            native = page.locator("#native")
            for key in ("0", "+", "ArrowRight"):
                native.locator(".linnet-viewport").press(key)
                result = native.evaluate(shade)
                assert result["curve"].endswith(" L14 40C30 10 50 40 70 40"), result
                assert result["cap"] == "butt", (
                    "shade extended across its half-edge boundary"
                )
                assert all(
                    abs(x - y) < 1e-6 for x, y in zip(result["start"], [10, 40])
                ), result
                assert result["halves"] == ["0"] and result["node"] == "0", result
                assert result["edgeHalves"] == 2, result
                assert all(
                    abs(x - y) < 1e-8
                    for x, y in zip(result["transform"], [1, 0, 0, 1, 0, 0])
                ), result
                assert result["selection"] == {
                    "nodes": [],
                    "edges": [],
                    "half_edges": [],
                }
            # Circlings outline the union once, with an independently tinted
            # interior. Layers are independent and do not modify selection.
            regions = native.evaluate("""svg => {
                const viewport=svg.querySelector('.linnet-viewport');
                const regions=[4, 8].map(padding=>svg.linnetShade(
                  {nodes:[0], half_edges:[0]},
                  {width:20, padding, outline:1.6, fillOpacity:.09}));
                viewport.append(...regions);
                const result=regions.map(g=>({
                  filter:g.querySelector('filter').id,
                  outline:+g.querySelector('feMorphology').getAttribute('radius'),
                  opacity:+g.style.opacity,
                  tint:+g.querySelector('feFuncA').getAttribute('slope'),
                  radius:+g.querySelector('circle').getAttribute('r'),
                  width:+g.querySelector('path').getAttribute('stroke-width'),
                  namespace:g.querySelector('feMorphology').namespaceURI,
                }));
                regions.forEach(g=>g.remove());
                return result;
            }""")
            assert regions[0]["filter"] != regions[1]["filter"]
            assert regions[1]["radius"] - regions[0]["radius"] == 4
            for region in regions:
                assert region["outline"] == 1.6 and region["opacity"] == 1
                assert region["tint"] == 0.09 and region["width"] == 20
                assert region["namespace"] == "http://www.w3.org/2000/svg"
            native.locator('[data-linnet-kind="halfedge"]').first.dispatch_event(
                "pointerover"
            )
            assert (
                native.locator(
                    '.linnet-edge-highlight [data-linnet-shade-halfedge="0"]'
                ).count()
                == 1
            )
            assert (
                native.locator(
                    '.linnet-edge-highlight [data-linnet-shade-halfedge="1"]'
                ).count()
                == 0
            )
            result = page.locator("#authored").evaluate("""svg => {
                const layer = svg.linnetShade({half_edges: [3]});
                svg.querySelector('.linnet-viewport').append(layer);
                const rect = layer.querySelector('rect'), matrix = rect.transform.baseVal.consolidate().matrix;
                return {count: layer.children.length, x: matrix.e, y: matrix.f};
            }""")
            assert result["count"] == 1, result
            assert abs(result["x"] - 10) < 1e-8 and abs(result["y"] - 20) < 1e-8, result
            assert not errors, errors
            print(
                "Native and authored half-edge shading, hover, and camera independence passed"
            )
        finally:
            browser.close()


if __name__ == "__main__":
    main()
