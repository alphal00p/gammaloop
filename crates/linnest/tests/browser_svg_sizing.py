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
            # Cross the scrollbar threshold in both directions, and cross
            # the 600px breakpoint between stacked and side-by-side details.
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

            # Different aspect ratios must preserve the same physical drawing
            # scale. Include a tall graph exceeding the old 520px height cap.
            page.set_content(
                f"<style>{(source / 'interactive.css').read_text()}</style>"
                + "".join(
                    f'<div style="width:550px"><svg id="{name}" '
                    f'data-linnet-interactive viewBox="0 0 {width} {height}" '
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
                        && bounds.top >= outer.top + 40 - .01
                        && bounds.bottom <= outer.bottom + .01,
                };
            }"""
            # Allow browser subpixel rounding of the outer SVG dimensions.
            for name in ("tall", "wide", "small"):
                result = page.locator(f"#{name}").evaluate(measure)
                assert abs(result["scale"] - 4 / 3) < 1e-3, result
                assert abs(result["diameter"] - 40 / 3) < 0.02, result
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
            viewport.press("0")
            assert wide.locator("output").inner_text() == "100%"
            assert wide.evaluate(measure)["contained"]
            assert not errors, errors
            print(
                "Consistent point scale, width fitting, fixed-height fitting, and zoom/reset passed"
            )
        finally:
            browser.close()


if __name__ == "__main__":
    main()
