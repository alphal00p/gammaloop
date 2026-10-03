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
        finally:
            browser.close()


if __name__ == "__main__":
    main()
