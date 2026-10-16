"""Check that child displays cannot navigate outside the represented object."""

import json
import sys
from pathlib import Path

from playwright.sync_api import sync_playwright


def check(directory, executable, engine="chromium"):
    with sync_playwright() as playwright:
        browser = getattr(playwright, engine).launch(
            executable_path=executable,
            args=["--no-sandbox"] if engine == "chromium" else [],
        )
        page = browser.new_page(viewport={"width": 800, "height": 1000})
        errors = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        for scope in ("orientation", "family"):
            page.set_content((Path(directory) / f"{scope}.html").read_text())
            data = json.loads(page.locator("[data-cff]").text_content())
            assert data["scope"] == scope
            assert len(data["orientations"]) == 1
            orientation = data["orientations"][0]
            assert page.locator("[data-total-sum]").is_hidden()
            assert page.locator("[data-sum-input]").count() == 0
            assert not page.locator("[data-explorer]").evaluate("node => node.open")
            equation = page.locator(".hs-factored-equation")
            assert (
                str(orientation["id"])
                in equation.locator(".hs-factored-lhs").inner_text()
            )
            page.locator("[data-explorer] > summary").click()
            page.wait_for_timeout(100)
            assert page.locator("[data-flip-edge]").count() == 0
            assert page.locator("[data-graph] > svg").count() == 1, errors
            if scope == "family":
                assert len(orientation["terms"]) == 1
                assert equation.locator("[data-family-leaf]").count() == 1
                assert page.locator("[data-family-previews]").is_hidden()
                assert page.locator(".hs-family-bar").is_hidden()
                actual = set(
                    page.locator("[data-definition]").evaluate_all(
                        "nodes => nodes.map(n => Number(n.dataset.definition))"
                    )
                )
                assert actual == set(orientation["terms"][0])
            else:
                lhs = equation.locator(".hs-factored-lhs").inner_text()
                choices = page.locator("[data-choose-family]")
                assert choices.count() > 1
                choices.nth(1).click()
                assert equation.locator(".hs-factored-lhs").inner_text() == lhs
                assert (
                    page.locator('[data-choose-family="1"]').get_attribute(
                        "aria-pressed"
                    )
                    == "true"
                )
            assert not errors, errors

        page.set_content((Path(directory) / "residue.html").read_text())
        data = json.loads(page.locator("[data-ltd]").text_content())
        assert data["scope"] == "residue" and len(data["residues"]) == 1
        residue = data["residues"][0]
        assert page.locator("[data-sum-input]").count() == 0
        assert not page.locator("[data-residue-graph]").evaluate("n => n.open")
        page.locator("[data-residue-graph] > summary").click()
        page.wait_for_timeout(100)
        lhs = page.locator(".lp-equation .hs-factored-lhs").inner_text()
        assert str(residue["id"]) in lhs
        cut = data["edges"][residue["cuts"][0]]
        cut_marker = page.locator(f'[data-graph-host] [data-pole-edge="{cut}"]')
        cut_marker.hover()
        cut_marker.click()
        assert page.locator(".lp-equation .hs-factored-lhs").inner_text() == lhs
        edge = data["edges"][int(next(iter(residue["propagator_factors"])))]
        marker = page.locator(f'[data-graph-host] [data-pole-edge="{edge}"]')
        before = page.locator(".lp-definition").inner_text()
        marker.click()
        assert page.locator(".lp-definition").inner_text() != before
        assert page.locator(".lp-equation .hs-factored-lhs").inner_text() == lhs
        assert not errors, errors
        for name in ("surface", "factor", "pair", "energy", "group"):
            page.set_content((Path(directory) / f"{name}.html").read_text())
            assert page.locator("[data-spenso-math]").count() > 0
            assert page.locator("math").count() > 0
            assert page.locator("math merror").count() == 0
            if name == "energy":
                page.locator("details > summary").click()
                assert page.locator("details [data-spenso-math]").is_visible()
            assert not errors, errors
        browser.close()
    print(
        "Scoped displays: fixed orientations, individual families, fixed residues and pair switching passed"
    )


if __name__ == "__main__":
    check(*sys.argv[1:])
