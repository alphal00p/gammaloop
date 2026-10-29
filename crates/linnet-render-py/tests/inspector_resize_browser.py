"""Check graph inspector hover with real scrollbars in inline and iframe outputs.

Use --html with a retained native DiagramView.to_html() capture. --source ROOT
refreshes only the production interactive JS/CSS, as iframe_display_browser.py
does, allowing focused browser checks without rebuilding the native extension.
"""

import argparse
import hashlib
import importlib.util
import json
import re
from pathlib import Path

import marimo
from mathml_view_browser import launch
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[3]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--html", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source = args.source or ROOT
    html = args.html.read_text()
    if args.source:
        for pattern, replacement in (
            (
                r"<script><!\[CDATA\[.*?\]\]></script>",
                "<script><![CDATA["
                + (source / "crates/linnest/src/svg/interactive.js").read_text()
                + "]]></script>",
            ),
            (
                r"<style>svg\[data-linnet-interactive\].*?</style>",
                "<style>"
                + (source / "crates/linnest/src/svg/interactive.css").read_text()
                + "</style>",
            ),
        ):
            html, count = re.subn(
                pattern,
                lambda _, replacement=replacement: replacement,
                html,
                flags=re.DOTALL,
            )
            assert count == 1, "Expected one native interactive graph"
    spec = importlib.util.spec_from_file_location(
        "iframe_adapter", source / "crates/linnet-render-py/src/view.py"
    )
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    bootstrap = (source / "crates/linnet-render-py/src/iframe.js").read_text()
    wrapper = adapter.iframe(html, bootstrap).text
    marimo_index = Path(marimo.__file__).parent / "_static/index.html"
    resize = re.search(
        r'<script data-marimo="true">(.*?)</script>',
        marimo_index.read_text(),
        re.DOTALL,
    ).group(1)
    report = {
        "passed": False,
        "html_sha256": hashlib.sha256(html.encode()).hexdigest(),
        "marimo_version": marimo.__version__,
        "cases": [],
    }
    try:
        with sync_playwright() as playwright:
            # Headless Chromium normally hides scrollbars. Keep them visible:
            # a disappearing 15px gutter is what moves the hover target.
            browser = launch(
                playwright,
                "chromium",
                ignore_default_args=["--hide-scrollbars"],
            )
            try:
                for mode in ("inline", "iframe", "shadow-iframe", "late-scroll"):
                    for theme in ("light", "dark"):
                        for width in (590, 350):
                            page = browser.new_page(
                                viewport={"width": 1100, "height": 900},
                                color_scheme=theme,
                            )
                            errors = []
                            page.on(
                                "pageerror",
                                lambda e, errors=errors: errors.append(str(e)),
                            )
                            background, color = (
                                ("#fff", "#000")
                                if theme == "light"
                                else ("#171e1a", "#eee")
                            )
                            page.set_content(
                                f'<body data-theme="{theme}"><style>'
                                f"body{{background:{background};color:{color};margin:20px}}"
                                f"#host{{width:{width}px;overflow:{'visible' if mode == 'late-scroll' else 'auto'}}}"
                                "</style><script>"
                                + resize
                                + '</script><div id="host">'
                                + (
                                    html
                                    if mode == "inline"
                                    else wrapper
                                    if mode in ("iframe", "late-scroll")
                                    else ""
                                )
                                + "</div>"
                            )
                            if mode == "shadow-iframe":
                                page.evaluate(
                                    "wrapper => host.attachShadow({mode:'open'}).innerHTML = wrapper",
                                    wrapper,
                                )
                            owner = page
                            if mode != "inline":
                                page.wait_for_function(
                                    "(host.shadowRoot || host).querySelector('iframe[data-linnet-ready]')"
                                )
                                owner = (
                                    page.locator("#host iframe")
                                    .element_handle()
                                    .content_frame()
                                )
                            owner.evaluate("document.fonts.ready")
                            owner.evaluate(
                                "() => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)))"
                            )
                            assert (
                                owner.locator(
                                    "svg[data-linnet-interactive]"
                                ).get_attribute("data-theme")
                                == theme
                            )
                            host = page.locator("#host")
                            host.evaluate(
                                "e => e.style.maxHeight = `${e.scrollHeight + 4}px`"
                            )
                            if mode == "late-scroll":
                                host.evaluate("e => e.style.overflow = 'auto'")
                                owner.evaluate(
                                    "() => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)))"
                                )
                            target = owner.locator('[data-linnet-kind="node"]').first
                            panel = owner.locator(".linnet-inspector")
                            initial_width = host.evaluate("e => e.clientWidth")
                            initial_target = target.bounding_box()
                            # Record before the pointer enters. Waiting for the
                            # inspector's final height misses single-frame jumps.
                            page.evaluate("""() => {
                              window.hoverFrames = [];
                              const iframe=(host.shadowRoot||host).querySelector('iframe');
                              const document=iframe?.contentDocument||window.document;
                              const target=document.querySelector('[data-linnet-kind="node"]');
                              const sample=()=>{
                                const box=target.getBoundingClientRect();
                                const offset=iframe?.getBoundingClientRect()||{x:0,y:0};
                                hoverFrames.push({width:host.clientWidth,target:{
                                  x:box.x+offset.x,y:box.y+offset.y,
                                  width:box.width,height:box.height}});
                                if(hoverFrames.length<24) requestAnimationFrame(sample);
                              };
                              sample();
                            }""")
                            page.mouse.move(
                                initial_target["x"] + initial_target["width"] / 2,
                                initial_target["y"] + initial_target["height"] / 2,
                            )
                            page.wait_for_function(
                                "host.scrollHeight > host.clientHeight"
                            )
                            page.wait_for_function("hoverFrames.length === 24")
                            frames = page.evaluate("hoverFrames")
                            widths = [frame["width"] for frame in frames]
                            assert widths == [initial_width] * len(widths), widths
                            assert all(
                                abs(frame["target"][key] - value) < 0.05
                                for frame in frames
                                for key, value in initial_target.items()
                            ), (mode, theme, width, initial_target, frames)
                            assert initial_width < width, (
                                "Real scrollbar gutter required"
                            )
                            assert all(
                                abs(value - initial_target[key]) < 0.05
                                for key, value in target.bounding_box().items()
                            ), (
                                mode,
                                theme,
                                width,
                                initial_target,
                                target.bounding_box(),
                            )
                            assert panel.is_visible()
                            target.dispatch_event("click")
                            page.mouse.move(1090, 890)
                            assert "Pinned" in panel.inner_text()
                            owner.get_by_role(
                                "button", name="Close graph details"
                            ).dispatch_event("click")
                            assert not panel.is_visible()
                            owner.evaluate(
                                "() => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)))"
                            )
                            assert host.evaluate("e => e.clientWidth") == initial_width
                            assert all(
                                abs(value - initial_target[key]) < 0.05
                                for key, value in target.bounding_box().items()
                            ), (
                                mode,
                                theme,
                                width,
                                initial_target,
                                target.bounding_box(),
                            )
                            assert not errors, errors
                            name = f"{mode}-{theme}-{width}"
                            page.locator("#host").screenshot(
                                path=str(args.output / f"{name}.png")
                            )
                            report["cases"].append(
                                {
                                    "name": name,
                                    "gutter": width - initial_width,
                                    "sampled_widths": widths,
                                    "hover_frames": frames,
                                    "target": initial_target,
                                    "hover_pin_close": True,
                                }
                            )
                            page.close()
            finally:
                browser.close()
        report["passed"] = True
    finally:
        (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
