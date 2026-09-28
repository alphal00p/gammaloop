"""Check revealed library grids with Playwright against an installed Spenso.

Run with a Playwright-enabled Python; --python selects the Spenso interpreter.
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

from playwright.sync_api import sync_playwright


def check_grid(frame):
    frame.locator(".tile").first.wait_for(state="visible")
    measurements = frame.evaluate("""() => {
        const rect = coordinate => document.querySelector(
            `[data-coordinate="${coordinate}"]`).getBoundingClientRect();
        return Array.from({length: 4}, (_, slice) => {
            const rows = Array.from({length: 4}, (_, row) => rect([slice,row,0]));
            return rows.slice(1).map((row,i) => row.top - rows[i].bottom);
        }).flat();
    }""")
    assert frame.locator(".tile").count() == 64
    assert all(abs(gap - 4) < 0.1 for gap in measurements), measurements


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--screenshot", type=Path)
    args = parser.parse_args()
    rendered = subprocess.run(
        [
            args.python,
            "-c",
            """from symbolica.community.spenso import TensorExpression, TensorLibrary
library = TensorLibrary.hep_lib_atom()
print('<section id="standalone">' + library[TensorExpression.gamma(4)].to_html()
      + '</section><section id="library">' + library.to_html() + '</section>')
""",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    with tempfile.TemporaryDirectory() as directory, sync_playwright() as driver:
        page_path = Path(directory) / "explorers.html"
        page_path.write_text(
            '<!doctype html><meta charset="utf-8">'
            '<body style="margin:24px;color-scheme:light">' + rendered
        )
        for engine in (driver.webkit, driver.chromium):
            browser = engine.launch()
            page = browser.new_page(viewport={"width": 1000, "height": 1800})
            page.goto(page_path.as_uri())
            chip = page.locator('.sl-chip[title="spenso::gamma"]')
            assert page.locator(".sl-panel:visible").count() == 0
            chip.click()
            standalone = (
                page.locator("#standalone iframe").element_handle().content_frame()
            )
            library = (
                page.locator(".sl-panel:visible iframe")
                .element_handle()
                .content_frame()
            )
            for width in (1000, 320, 736, 1000):
                page.set_viewport_size({"width": width, "height": 1800})
                check_grid(standalone)
                check_grid(library)
            library.locator('[data-display="matrix"]').click()
            assert library.locator(".matrix td").count() == 64
            library.locator('[data-display="grid"]').click()
            check_grid(library)
            page.locator(".sl-panel:visible .sl-close").click()
            chip.click()
            check_grid(library)
            if args.screenshot and engine.name == "webkit":
                page.screenshot(path=str(args.screenshot), full_page=True)
            browser.close()
            print(f"{engine.name}: revealed grids retain 4 px row gaps")


if __name__ == "__main__":
    main()
