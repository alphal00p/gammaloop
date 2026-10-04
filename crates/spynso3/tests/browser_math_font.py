"""Check stretchy brackets, shared fonts, offline setup, and responsive MathML.

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
from symbolica.community import tensor as spenso
expression = spenso.TensorExpression(E('(a/(b+c)+d/(e+f))*(g+h)+z'))
long_sum = spenso.TensorExpression(E('+'.join(f'a{i}/(b{i}+c{i})' for i in range(12))))
wide_term = spenso.TensorExpression(E('*'.join(f'x{i}' for i in range(30)) + '+z'))
wide_product = spenso.TensorExpression(E('*'.join(f'x{i}' for i in range(30))))
wide_fraction = spenso.TensorExpression(E('(' + '*'.join(f'x{i}' for i in range(30)) + '+z)/(D+1)'))
print(json.dumps([expression.to_html(), spenso.load_math_font(),
                  long_sum.to_html(), wide_term.to_html(),
                  wide_product.to_html(), wide_fraction.to_html()]))
""",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    html, setup, long_sum, wide_term, wide_product, wide_fraction = json.loads(
        result.stdout
    )
    assert "data:font" not in html
    assert len(html.encode()) < 10_000
    assert setup.count("data:font/woff2;base64,") == 1
    assert "SIL OPEN FONT LICENSE" in setup
    assert "data-spenso-wrap" in long_sum
    assert "<script" not in long_sum

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
            output = (
                html * 10
                + f'<section id="sum" style="width:1200px">{long_sum}</section>'
                + f'<section id="wide" style="width:240px">{wide_term}</section>'
                + f'<section id="product" style="width:240px">{wide_product}</section>'
                + f'<section id="fraction" style="width:240px">{wide_fraction}</section>'
            ).replace('local("STIX Two Math"),', "")
            page.set_content((setup if offline else "") + output)
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

            original_math = page.locator("#sum math").inner_html()
            rows = []
            for width in (1200, 280, 1200):
                layout = page.evaluate(
                    """width => {
                        const section = document.querySelector('#sum');
                        section.style.width = width + 'px';
                        const box = section.querySelector('[data-spenso-math]');
                        const terms = [...section.querySelectorAll('[data-spenso-term]')];
                        return {
                            rows: new Set(terms.map(x => Math.round(x.getBoundingClientRect().y))).size,
                            width: box.clientWidth,
                            scroll: box.scrollWidth,
                            fractions: section.querySelectorAll('mfrac').length,
                            signs: terms.slice(1).map(x => x.firstElementChild.textContent),
                        };
                    }""",
                    width,
                )
                assert layout["scroll"] <= layout["width"] + 1, layout
                assert layout["fractions"] == 12, layout
                assert layout["signs"] == ["+"] * 11, layout
                rows.append(layout["rows"])
            assert rows[0] == rows[2] == 1 and rows[1] > 1, rows
            assert page.locator("#sum math").inner_html() == original_math
            # A scrollbar alone is insufficient: centered overflowing math can
            # put its leading factors outside the reachable scroll area.
            for selector in ("#wide", "#product", "#fraction"):
                bounds = page.locator(f"{selector} [data-spenso-math]").evaluate(
                    """box => {
                        const children = [...box.querySelector('math').children];
                        box.scrollLeft = 0;
                        const left = Math.min(...children.map(
                            x => x.getBoundingClientRect().left));
                        box.scrollLeft = box.scrollWidth;
                        const right = Math.max(...children.map(
                            x => x.getBoundingClientRect().right));
                        box.scrollLeft = 0;
                        const viewport = box.getBoundingClientRect();
                        return {left, right, start: viewport.left, end: viewport.right,
                                width: box.clientWidth, scroll: box.scrollWidth};
                    }"""
                )
                assert bounds["scroll"] > bounds["width"], (selector, bounds)
                assert bounds["left"] >= bounds["start"] - 1, (selector, bounds)
                assert bounds["right"] <= bounds["end"] + 1, (selector, bounds)
            print(
                f"offline={offline}: sum rows at 1200/280/1200px = {rows}; wide term scrolls"
            )
            context.close()
        browser.close()


if __name__ == "__main__":
    main()
