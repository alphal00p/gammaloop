"""Check generalized numerator math, powers and higher-pole graph annotations."""

import sys
from pathlib import Path

from playwright.sync_api import sync_playwright


def assert_baselines(page):
    groups = page.locator(".hs-factored-equation").evaluate("""root => {
        const baseline = element => {
            const math = element.namespaceURI === 'http://www.w3.org/1998/Math/MathML';
            const probe = math
                ? document.createElementNS(element.namespaceURI, 'mspace')
                : document.createElement('span');
            if (math)
                for (const key of ['width', 'height', 'depth']) probe.setAttribute(key, '0px');
            else probe.style = 'display:inline-block;width:0;height:0';
            element.prepend(probe);
            const y = probe.getBoundingClientRect().y;
            probe.remove();
            return y;
        };
        return [
            [...root.querySelectorAll('.hs-factored-lhs,.hs-factored-fraction')],
            ...[...root.querySelectorAll('.hs-factored-denominator')].map(denominator =>
                [...denominator.querySelectorAll('.hs-factored-factor,.hs-energy math')])
        ].filter(group => group.length > 1).map(group => group.map(baseline));
    }""")
    for group in groups:
        assert max(group) - min(group) < 1, group


def assert_math_fonts(page, engine):
    # Repeated selections must not reinstall a stylesheet with every factor.
    assert not page.locator(".hs-energy style, .lp-math style").count()
    selectors = [
        ".hs-sum-math",
        ".hs-factored-lhs",
        ".hs-factored-factor",
        ".hs-energy mi",
    ]
    families = [
        page.locator(selector).first.evaluate("e => getComputedStyle(e).fontFamily")
        for selector in selectors
    ]
    assert len(set(families)) == 1, families
    if engine == "chromium":
        # Computed styles alone cannot detect silent fallback to the UI font.
        session = page.context.new_cdp_session(page)
        session.send("DOM.enable")
        session.send("CSS.enable")
        root = session.send("DOM.getDocument")["root"]["nodeId"]
        for selector in selectors:
            node = session.send(
                "DOM.querySelector", {"nodeId": root, "selector": selector}
            )
            fonts = session.send("CSS.getPlatformFontsForNode", node)["fonts"]
            assert fonts and all("STIX" in font["familyName"] for font in fonts), fonts
        session.detach()


def main(directory, ltd_directory, executable, engine="chromium"):
    with sync_playwright() as p:
        browser = getattr(p, engine).launch(
            executable_path=executable,
            args=["--no-sandbox"] if engine == "chromium" else [],
        )
        page = browser.new_page(viewport={"width": 736, "height": 1100})
        page.route(
            "**/*STIXTwoMath-Regular.woff2",
            lambda route: route.fulfill(
                path=Path(__file__).parents[2]
                / "spynso3/typst/STIXTwoMath-Regular.woff2",
                content_type="font/woff2",
                headers={"Access-Control-Allow-Origin": "*"},
            ),
        )
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        for source in sorted(directory.glob("cff-degree-*.html")):
            print(f"{engine}: {source.name}", flush=True)
            # Measure one unwrapped math row; small-screen wrapping is checked below.
            page.set_viewport_size({"width": 3200, "height": 1100})
            page.set_content(source.read_text())
            page.wait_for_selector(".feynkit-cff-result[data-ready]")
            page.evaluate("document.fonts.load('19px \"STIX Two Math\"'); true")
            page.wait_for_function("document.fonts.check('19px \"STIX Two Math\"')")
            assert_math_fonts(page, engine)
            assert_baselines(page)
            equation = page.locator(".hs-factored-equation")
            equation.evaluate("e => e.style.fontSize = '26px'")
            assert_baselines(page)
            equation.evaluate("e => e.style.removeProperty('font-size')")
            page.evaluate("""() => {
                const root = document.querySelector('.feynkit-cff-result');
                const data = JSON.parse(root.querySelector('[data-cff]').textContent);
                for (const orientation of data.orientations) {
                    const input = root.querySelector('[data-sum-input]');
                    input.value = orientation.id;
                    input.dispatchEvent(new KeyboardEvent('keydown', {key:'Enter', bubbles:true}));
                    for (const leaf of root.querySelectorAll('[data-family-leaf]')) {
                        const i = +leaf.dataset.familyLeaf, contribution = orientation.contributions[i];
                        const collect = target => [...target.querySelectorAll('[data-factor-id]')].flatMap(n=>
                            Array(Number(n.dataset.factorPower || 1)).fill(n.dataset.factorId));
                        const factors = collect(leaf);
                        for (let n=leaf.parentElement; n && !n.matches('[data-orientation-sum]'); n=n.parentElement)
                            if (n.matches('.hs-factor-node')) {
                                const own=n.querySelector(':scope > .hs-factored-fraction');
                                if(own)factors.push(...collect(own));
                            }
                        const expected = [...orientation.terms[i].map(String), ...contribution.energies.map(e=>`energy-${e}`)];
                        if(JSON.stringify(factors.sort()) !== JSON.stringify(expected.sort()))throw Error('lost factor power');
                        if(leaf.dataset.coefficient !== contribution.coefficient)throw Error('lost coefficient');
                        if(contribution.numerator !== '1' && !leaf.querySelector('math'))throw Error('missing numerator math');
                        if(!leaf.title.includes('q'))throw Error('missing numerator map explanation');
                    }
                }
            }""")
            assert not page.locator("merror").count()
            assert all(
                "OSE" not in text for text in page.locator("math").all_text_contents()
            )
            for energy in page.locator('[data-factor-id^="energy-"]').all():
                assert energy.locator("math").count() == 1
                assert "os" in energy.inner_text()
            page.locator("[data-explorer] > summary").click()
            page.locator("[data-graph] > svg").wait_for(state="visible")
            for width in (736, 360):
                page.set_viewport_size({"width": width, "height": 1100})
                assert page.evaluate(
                    "document.documentElement.scrollWidth <= innerWidth + 1"
                )
            assert not errors, errors
        page.set_viewport_size({"width": 3200, "height": 1100})
        page.set_content((ltd_directory / "ltd-raised.html").read_text())
        page.evaluate("document.fonts.load('19px \"STIX Two Math\"'); true")
        page.wait_for_function("document.fonts.check('19px \"STIX Two Math\"')")
        assert_math_fonts(page, engine)
        assert_baselines(page)
        assert page.locator("[data-power-label]").count() > 0
        page.locator("[data-residue-graph] > summary").click()
        page.locator(".lp-graph").wait_for(state="visible")
        assert page.locator('[data-pole-order="2"]').count() > 0
        for energy in page.locator('[data-factor-id^="energy-"]').all():
            assert energy.locator("math").count() == 1
            assert "os" in energy.inner_text()
        assert not errors, errors
        print(
            f"{engine}: math baselines, multiplicities, maps and raised-pole markers passed"
        )
        browser.close()


if __name__ == "__main__":
    main(
        Path(sys.argv[1]),
        Path(sys.argv[2]),
        sys.argv[3],
        sys.argv[4] if len(sys.argv) > 4 else "chromium",
    )
