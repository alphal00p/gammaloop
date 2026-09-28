"""Check stretchy MathML brackets, shared font loading, and offline setup.

Run with a Playwright-enabled Python; --python selects the installed Spenso.
The CDN request is fulfilled locally so this regression test works offline.
Use --font to check a separately downloaded copy of the pinned CDN font.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--chromium")
    parser.add_argument(
        "--font",
        type=Path,
        default=Path(__file__).parents[1] / "typst/STIXTwoMath-Regular.woff2",
    )
    args = parser.parse_args()
    result = subprocess.run(
        [
            args.python,
            "-c",
            """import json
from symbolica import E
from symbolica.community import spenso
expression = spenso.TensorExpression(E('(a/(b+c)+d/(e+f))*(g+h)'))
print(json.dumps([expression.to_html(), spenso.load_math_font()]))
""",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    html, setup = json.loads(result.stdout)
    assert "data:font" not in html
    assert len(html.encode()) < 10_000
    assert setup.count("data:font/woff2;base64,") == 1
    assert "SIL OPEN FONT LICENSE" in setup

    with sync_playwright() as driver:
        browser = driver.chromium.launch(executable_path=args.chromium)
        for offline in (False, True):
            context = browser.new_context()
            page = context.new_page()
            requests = []

            def font_request(route, request, requests=requests):
                requests.append(route.request.url)
                assert "/stixfonts@2.13b171/" in route.request.url
                route.fulfill(
                    body=args.font.read_bytes(),
                    content_type="font/woff2",
                    headers={"Access-Control-Allow-Origin": "*"},
                )

            context.route("https://cdn.jsdelivr.net/**", font_request)
            # Force the webfont path even on machines with STIX installed.
            html_without_local = html.replace('local("STIX Two Math"),', "")
            page.set_content((setup if offline else "") + html_without_local * 10)
            page.evaluate("() => document.fonts.ready")
            measurements = page.evaluate("""() => {
                const math = document.querySelector('math');
                const outer = math.querySelector('mfrac').parentElement.parentElement;
                const brackets = [outer.firstElementChild, outer.lastElementChild];
                const content = outer.children[1];
                return {
                    brackets: brackets.map(x => x.getBoundingClientRect().height),
                    content: content.getBoundingClientRect().height,
                    font: [...document.fonts].filter(x => x.status === 'loaded')
                                             .map(x => x.family)
                };
            }""")
            assert all(
                height >= measurements["content"] * 0.95
                for height in measurements["brackets"]
            ), measurements
            expected_font = "Spenso Offline Math" if offline else "STIX Two Math"
            assert expected_font in measurements["font"], measurements
            assert len(requests) == (0 if offline else 1), requests
            print(f"offline={offline}: {len(requests)} requests; {measurements}")
            context.close()
        browser.close()


if __name__ == "__main__":
    main()
