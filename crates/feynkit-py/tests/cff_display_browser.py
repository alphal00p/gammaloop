"""Run after CFF_DISPLAY_TEST_OUTPUT=DIR cargo test -p feynkit-py --test cff_display.

Usage: python cff_display_browser.py DIR [browser-executable] [chromium|webkit]
"""

import html
import json
import sys
from pathlib import Path

from generalized_energy_display_browser import assert_baselines
from playwright.sync_api import sync_playwright


def payload(page):
    return json.loads(page.locator("[data-cff]").first.text_content())


def choose_orientation(page, index):
    root = page.locator(".feynkit-cff-result").first
    field = root.locator("[data-sum-input]")
    field.fill(str(index))
    field.press("Enter")


def assert_factors(page, orientation):
    families = page.locator(".feynkit-cff-result").first.evaluate("""root => {
      return [...root.querySelectorAll('[data-family-leaf]')].map(leaf => {
        const factors = [...leaf.querySelectorAll('[data-factor-id]')].flatMap(n => Array(Number(n.dataset.factorPower || 1)).fill(+n.dataset.factorId));
        for (let node=leaf.parentElement; node && !node.matches('[data-orientation-sum]'); node=node.parentElement)
          if (node.classList.contains('hs-factor-node'))
            factors.push(...[...node.querySelectorAll(':scope > .hs-factored-fraction [data-factor-id]')].flatMap(n=>Array(Number(n.dataset.factorPower || 1)).fill(+n.dataset.factorId)));
        return [+leaf.dataset.familyLeaf, factors.sort((a,b)=>a-b)];
      });
    }""")
    assert len(families) == len(orientation["terms"])
    for family, factors in families:
        assert factors == sorted(orientation["terms"][family]), (family, factors)


def main(directory, executable=None, engine="chromium"):
    sources = sorted(Path(directory).glob("cff-*.html"))
    assert sources, "Run the native CFF test first to export display fixtures"
    with sync_playwright() as playwright:
        browser = getattr(playwright, engine).launch(
            executable_path=executable,
            args=["--no-sandbox"] if engine == "chromium" else [],
        )
        page = browser.new_page(viewport={"width": 800, "height": 1100})
        errors = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        for source in sources:
            page.set_content(source.read_text())
            assert page.locator(".feynkit-cff-result[data-ready=true]").count() == 1
            data = payload(page)
            if data["orientations"]:
                assert_factors(page, data["orientations"][0])
                assert page.locator("[data-sum-input]").count() == 1
                if len(data["orientations"]) == 1:
                    assert page.locator("[data-sum-id]").count() == 1
                    assert page.locator("[data-sum-index]").count() == 0
            else:
                assert page.locator("[data-sum-input]").count() == 0
            assert not errors, errors

        # Nested sums used to wrap at every level in WebKit. Check at a normal
        # notebook cell width, not just a desktop-wide standalone document.
        # Preserve real line breaks on small screens, without overflow.
        page.set_content((Path(directory) / "cff-2.html").read_text())
        choose_orientation(page, 97)
        equation = page.locator(".hs-factored-equation")
        page.set_viewport_size({"width": 736, "height": 1100})
        tops = equation.locator("[data-family-leaf]").evaluate_all(
            "nodes=>nodes.map(n=>n.getBoundingClientRect().top)"
        )
        assert len(tops) == 5
        assert max(tops) - min(tops) < 2, tops
        height = equation.bounding_box()["height"]
        prefactor = equation.locator(".hs-factored-fraction[title]")
        assert prefactor.locator("mfrac > mtext").first.text_content() == "−1"
        # The prefactor and every family fraction share one mathematical axis,
        # even though the family captions have their own height below it.
        page.set_viewport_size({"width": 1200, "height": 1100})
        assert_baselines(page)
        page.set_viewport_size({"width": 320, "height": 1100})
        assert equation.bounding_box()["height"] > height
        assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
        page.set_viewport_size({"width": 800, "height": 1100})

        page.set_content((Path(directory) / "collection.html").read_text())
        thumbnail = page.locator(".fk-thumbnail").first
        page.wait_for_function(
            "document.querySelector('.fk-thumbnail').naturalWidth > 0"
        )
        light = thumbnail.get_attribute("src")
        page.locator("body").evaluate("body=>body.dataset.theme='dark'")
        page.wait_for_function(
            "document.querySelector('.feynkit-collection').dataset.theme==='dark'"
        )
        assert thumbnail.get_attribute("src") != light
        assert page.locator(".fk-stage > svg").count() == 1

        source = (Path(directory) / "cff-3.html").read_text()
        page.set_content(
            '<body data-theme="dark"><iframe style="width:736px;border:0" srcdoc="'
            + html.escape(source, quote=True)
            + '"></iframe></body>'
        )
        frame = page.frames[1]
        frame.wait_for_selector('.feynkit-cff-result[data-theme="dark"]')
        explorer = frame.locator("[data-explorer]")
        assert not explorer.evaluate("e=>e.open")
        assert frame.locator("[data-graph] svg").count() == 0
        collapsed_height = frame.locator(".feynkit-cff-result").bounding_box()["height"]
        # Navigation and definitions work before the graph is ever mounted.
        choose_orientation(frame, 1)
        frame.locator(".hs-options > summary").click()
        frame.locator("[data-multiple]").check()
        for button in frame.locator("[data-orientation-sum] button").all():
            if button.get_attribute("aria-pressed") == "false":
                button.click()
        page.wait_for_function(
            "Math.abs(parseFloat(document.querySelector('iframe').style.height) - (document.querySelector('iframe').contentDocument.querySelector('.feynkit-cff-result').getBoundingClientRect().bottom + 14)) < 2"
        )
        explorer.locator("summary").click()
        frame.locator("[data-graph] .linnet-viewport").wait_for(state="visible")
        page.wait_for_function(
            "parseFloat(document.querySelector('iframe').style.height) > "
            + str(collapsed_height)
        )
        assert (
            frame.locator("[data-graph] > svg")
            .get_attribute("aria-label")
            .startswith("Energy-flow orientation 1;")
        )
        explorer.locator("summary").click()
        frame.locator(".hs-options > summary").click()
        page.wait_for_function(
            "Math.abs(parseFloat(document.querySelector('iframe').style.height) - (document.querySelector('iframe').contentDocument.querySelector('.feynkit-cff-result').getBoundingClientRect().bottom + 14)) < 2"
        )
        page.set_content(source + source)
        data = payload(page)
        root = page.locator(".feynkit-cff-result").first
        other = page.locator(".feynkit-cff-result").nth(1)
        for display in (root, other):
            display.locator("[data-explorer] > summary").click()
            display.locator("[data-graph] .linnet-viewport").wait_for(state="visible")
        ids = page.locator("[id]").evaluate_all("nodes => nodes.map(n=>n.id)")
        assert len(ids) == len(set(ids)), "SVG glyph IDs collide across displays"
        # The sum is also navigation: endpoints, adjacent terms, and one ID editor.
        choose_orientation(page, 10)
        assert root.locator("[data-sum-id]").evaluate_all(
            "nodes=>nodes.map(n=>+n.dataset.sumId)"
        ) == [0, 9, 10, 11, 685]
        assert root.locator(".hs-sum-math").inner_text().count("⋯") == 2
        root.get_by_role("button", name="Show orientation 9", exact=True).click()
        assert root.locator("[data-sum-input]").input_value() == "9"
        root.get_by_role("button", name="Show orientation 10", exact=True).press(
            "Enter"
        )
        assert root.locator("[data-sum-input]").input_value() == "10"
        root.get_by_role("button", name="Show orientation 11", exact=True).click()
        assert root.locator("[data-sum-input]").input_value() == "11"
        root.get_by_role("button", name="Show orientation 685", exact=True).click()
        assert root.locator("[data-sum-id]").evaluate_all(
            "nodes=>nodes.map(n=>+n.dataset.sumId)"
        ) == [0, 684, 685]
        root.get_by_role("button", name="Show orientation 0", exact=True).click()
        assert root.locator("[data-sum-id]").evaluate_all(
            "nodes=>nodes.map(n=>+n.dataset.sumId)"
        ) == [0, 1, 685]
        field = root.locator("[data-sum-input]")
        field.fill("10")
        field.press("Escape")
        assert field.input_value() == "0"
        field.fill("10")
        root.get_by_role("button", name="Show orientation 1", exact=True).click()
        assert field.input_value() == "1", "Leaving an edit must not swallow navigation"
        for invalid in ("", "-1", "1.5", "686"):
            field.fill(invalid)
            field.press("Enter")
            assert field.input_value() == "1"
            assert root.locator('[data-status][data-error="true"]').count() == 1
        assert other.locator("[data-sum-input]").input_value() == "0"
        orientation = next(o for o in data["orientations"] if len(o["terms"]) == 5)
        choose_orientation(page, orientation["id"])
        root.locator('[data-choose-family="1"]').click()
        assert (
            root.locator('[data-choose-family="1"]').get_attribute("aria-pressed")
            == "true"
        )
        assert other.locator("[data-sum-input]").input_value() == "0"
        assert_factors(page, orientation)
        factors = orientation["terms"][1]
        root.locator(f'[data-orientation-sum] [data-surface="{factors[0]}"]').click()
        root.locator(f'[data-orientation-sum] [data-surface="{factors[1]}"]').click(
            modifiers=["Shift"]
        )
        assert root.locator("[data-definition]").count() == 2
        bands = root.locator("[data-graph] [data-boundary-half]").evaluate_all(
            "nodes=>nodes.map(n=>({surface:+n.dataset.surface,owner:+n.dataset.owner,edge:+n.dataset.boundaryHalf}))"
        )
        assert {b["surface"] for b in bands} == set(factors[:2])
        for band in bands:
            surface = data["surfaces"][band["surface"]]
            assert band["owner"] in surface["v"]
            assert band["edge"] in surface["e"] + surface["negative"]
        root.locator(".hs-options > summary").click()
        root.locator("[data-contours]").check()
        root.locator(".hs-options > summary").click()
        assert root.locator("[data-graph] filter").count() == 2
        # Nested regions need separate boundaries. Their geometry must stay put
        # when another surface is selected, and must follow only internal edges.
        outline_geometry = (
            "e=>[...e.querySelectorAll('circle,path')].map(n=>n.outerHTML)"
        )
        first_outline = root.locator(f'[data-outline-surface="{factors[0]}"]').evaluate(
            outline_geometry
        )
        for surface in factors[2:]:
            root.locator(f'[data-orientation-sum] [data-surface="{surface}"]').click(
                modifiers=["Shift"]
            )
        outlines = root.locator("[data-outline-surface]").evaluate_all("""nodes =>
          nodes.map(n => ({id:+n.dataset.outlineSurface,
            vertices:[...n.querySelectorAll('circle')].map(c=>+c.dataset.linnetShadeNode),
            edges:[...n.querySelectorAll('path')].map(p=>JSON.parse(n.closest('svg[data-linnet-interactive]').querySelector(`[data-linnet-carrier][data-linnet-id="${p.dataset.linnetShadeHalfedge}"]`).dataset.linnetDetail).edge),
            outline:+n.querySelector('feMorphology').getAttribute('radius'),
            opacity:+getComputedStyle(n).opacity,
            radius:+n.querySelector('circle').getAttribute('r')}))
        """)
        assert {o["id"] for o in outlines} == set(factors)
        for outline in outlines:
            assert outline["outline"] >= 1.5 and outline["opacity"] == 1
            members = set(data["surfaces"][outline["id"]]["v"])
            assert set(outline["vertices"]) == members
            assert set(outline["edges"]) == {
                edge["id"]
                for edge in data["edges"]
                if edge["source"] in members and edge["target"] in members
            }
            for inner in outlines:
                if set(inner["vertices"]) < members:
                    assert outline["radius"] > inner["radius"] + 1
        for surface in factors[1:]:
            root.locator(f'[data-orientation-sum] [data-surface="{surface}"]').click(
                modifiers=["Shift"]
            )
        assert root.locator("[data-outline-surface]").count() == 1
        assert (
            root.locator(f'[data-outline-surface="{factors[0]}"]').evaluate(
                outline_geometry
            )
            == first_outline
        )
        root.locator(f'[data-orientation-sum] [data-surface="{factors[1]}"]').click(
            modifiers=["Shift"]
        )
        graph = root.locator("[data-graph]")
        camera = graph.locator(".linnet-viewport")
        graph.scroll_into_view_if_needed()
        # The native camera centers the drawing without enlarging it to fill.
        geometry = root.evaluate("""root => {
          const original = root.querySelector('[data-drawing]').content.querySelector('svg').viewBox.baseVal;
          const camera = root.querySelector('.linnet-viewport');
          const box = camera.viewBox.baseVal;
          return {scale:camera.getScreenCTM().a,
            dx:(box.x+box.width/2)-(original.x+original.width/2),
            dy:(box.y+box.height/2)-(original.y+original.height/2)};
        }""")
        assert 0 < geometry["scale"] <= 4 / 3 + 1e-6, geometry
        assert abs(geometry["dx"]) < 1e-4 and abs(geometry["dy"]) < 1e-4, geometry
        # Edge bodies do not change orientation; only their arrowheads do.
        arrow = graph.locator('[data-flip-edge][aria-disabled="false"]').first
        edge = arrow.get_attribute("data-flip-edge")
        before = root.locator("[data-sum-input]").input_value()
        point = graph.locator(
            f'path[data-edge="{edge}"][data-flow="source"]'
        ).first.evaluate("""path => {
          const p=path.getPointAtLength(path.getTotalLength()*.1);
          const screen=new DOMPoint(p.x,p.y).matrixTransform(path.getScreenCTM());
          return {x:screen.x,y:screen.y};
        }""")
        page.mouse.click(**point)
        assert root.locator("[data-sum-input]").input_value() == before
        camera.focus()
        camera.press("+")
        zoomed = camera.get_attribute("viewBox")
        assert camera.evaluate("s=>s.getScreenCTM().a") > geometry["scale"]
        arrow.click()
        assert root.locator("[data-sum-input]").input_value() != before
        assert camera.get_attribute("viewBox") == zoomed
        flipped = graph.locator(f'[data-flip-edge="{edge}"]')
        assert flipped.evaluate("node=>node===document.activeElement")
        flipped.press("Enter")
        assert root.locator("[data-sum-input]").input_value() == before
        flipped.press("Space")
        assert root.locator("[data-sum-input]").input_value() != before
        assert other.locator("[data-sum-input]").input_value() == "0"
        before = root.locator("[data-sum-input]").input_value()
        bounds = flipped.bounding_box()
        page.mouse.move(
            bounds["x"] + bounds["width"] / 2, bounds["y"] + bounds["height"] / 2
        )
        page.mouse.down()
        page.mouse.move(bounds["x"] + 50, bounds["y"] + 30, steps=5)
        page.mouse.up()
        assert root.locator("[data-sum-input]").input_value() == before
        panned = camera.get_attribute("viewBox")
        assert panned != zoomed
        root.locator("[data-explorer] > summary").click()
        # Updating a collapsed expression must leave its mounted camera intact.
        root.locator("[data-orientation-sum] button").first.click()
        root.locator("[data-explorer] > summary").click()
        camera.wait_for(state="visible")
        assert camera.get_attribute("viewBox") == panned
        root.locator("[data-orientation-sum] button").first.click()
        assert camera.get_attribute("viewBox") == panned
        choose_orientation(page, orientation["id"])
        root.locator('[data-choose-family="1"]').click()
        assert camera.get_attribute("viewBox") == panned
        graph.scroll_into_view_if_needed()
        point = camera.bounding_box()
        page.mouse.move(
            point["x"] + point["width"] / 2, point["y"] + point["height"] / 2
        )
        scale = camera.evaluate("s=>s.getScreenCTM().a")
        page.keyboard.down("Control")
        page.mouse.wheel(0, -100)
        page.keyboard.up("Control")
        page.wait_for_function(
            "s=>s.getScreenCTM().a > " + str(scale), arg=camera.element_handle()
        )
        camera.focus()
        camera.press("0")
        assert abs(camera.evaluate("s=>s.getScreenCTM().a") - geometry["scale"]) < 1e-6
        root.locator("[data-sum-input]").fill("-1")
        root.locator("[data-sum-input]").press("Enter")
        assert root.locator('[data-status][data-error="true"]').count() == 1
        root.locator("[data-explorer] > summary").click()
        # Verify the displayed factorization against every native family, including
        # multiplicities, not just the selected branch or numerical samples.
        for orientation in data["orientations"]:
            choose_orientation(page, orientation["id"])
            assert_factors(page, orientation)
        choose_orientation(
            page, next(o["id"] for o in data["orientations"] if len(o["terms"]) == 5)
        )
        for width in (320, 360, 736):
            page.set_viewport_size({"width": width, "height": 1100})
            for expanded in (False, True):
                if root.locator("[data-explorer]").evaluate("e=>e.open") != expanded:
                    root.locator("[data-explorer] > summary").click()
                if expanded:
                    camera.wait_for(state="visible")
                assert page.evaluate(
                    "document.documentElement.scrollWidth <= innerWidth"
                ), (width, expanded)
        root.locator("[data-explorer] > summary").click()
        page.locator("body").evaluate("body=>body.dataset.theme='dark'")
        page.wait_for_function(
            "document.querySelector('.feynkit-cff-result').dataset.theme==='dark'"
        )
        assert root.evaluate("e=>getComputedStyle(e).colorScheme") == "dark"
        page.screenshot(
            path=str(Path(directory) / "cff-mobile-dark.png"), full_page=True
        )
        page.set_viewport_size({"width": 800, "height": 1100})
        page.locator("body").evaluate("body=>body.dataset.theme='light'")
        root.screenshot(path=str(Path(directory) / "cff-desktop.png"))
        assert not errors, errors
        browser.close()
        print(
            "CFF browser: all 3432 families, expression wrapping and prefactor, sum navigation and ID editing, arrowhead clicks, centered camera, zoom/pan persistence, multiselection, nested outlines, instance isolation, themes and mobile widths passed"
        )


if __name__ == "__main__":
    main(*sys.argv[1:])
