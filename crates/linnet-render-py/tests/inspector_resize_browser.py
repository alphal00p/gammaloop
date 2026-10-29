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
    parser.add_argument("--engine", choices=("chromium", "firefox"), default="chromium")
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
                args.engine,
                **(
                    {"ignore_default_args": ["--hide-scrollbars"]}
                    if args.engine == "chromium"
                    else {}
                ),
            )
            try:
                for mode in (
                    "inline",
                    "iframe",
                    "shadow-iframe",
                    "late-scroll",
                    "opaque-iframe",
                ):
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
                            if mode == "opaque-iframe":
                                # A custom notebook can sandbox its frame without
                                # allow-same-origin. Only the host can size it;
                                # graph code cannot reserve the outer gutter.
                                page.evaluate(
                                    r"""({html,color}) => {
                                  const frame = document.createElement('iframe');
                                  frame.sandbox = 'allow-scripts';
                                  frame.style.cssText = 'display:block;width:100%;border:0';
                                  addEventListener('message', event => {
                                    if (event.source === frame.contentWindow && event.data.height)
                                      frame.style.height = `${event.data.height}px`;
                                  });
                                  frame.srcdoc = `<style>html{overflow:hidden}body{margin:0;color:${color}}</style>`
                                    + html + `<script>new ResizeObserver(() => {
                                      const range=document.createRange();range.selectNodeContents(document.body);
                                      parent.postMessage({height:Math.ceil(range.getBoundingClientRect().bottom+scrollY)},'*');
                                    }).observe(document.body);<\/script>`;
                                  host.append(frame);
                                }""",
                                    {"html": html, "color": color},
                                )
                            owner = page
                            if mode != "inline":
                                if mode != "opaque-iframe":
                                    page.wait_for_function(
                                        "(host.shadowRoot || host).querySelector('iframe[data-linnet-ready]')"
                                    )
                                owner = (
                                    page.locator("#host iframe")
                                    .element_handle()
                                    .content_frame()
                                )
                                owner.wait_for_selector(".linnet-toolbar")
                            owner.evaluate("document.fonts.ready")
                            owner.evaluate(
                                "() => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)))"
                            )
                            if mode != "opaque-iframe":
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
                            popup_owner = owner if mode == "opaque-iframe" else page
                            panel = popup_owner.locator(".linnet-inspector")
                            initial_width = host.evaluate("e => e.clientWidth")
                            initial_height = host.evaluate("e => e.scrollHeight")
                            initial_target = target.bounding_box()
                            local_target = target.evaluate(
                                "e => {const {x,y,width,height}=e.getBoundingClientRect();return {x,y,width,height}}"
                            )
                            # Record before the pointer enters. Waiting for the
                            # inspector's final height misses single-frame jumps.
                            owner.evaluate("""() => {
                              window.hoverFrames = [];
                              const target=document.querySelector('[data-linnet-kind="node"]');
                              const sample=()=>{
                                const box=target.getBoundingClientRect();
                                hoverFrames.push({width:document.documentElement.clientWidth,target:{
                                  x:box.x,y:box.y,
                                  width:box.width,height:box.height}});
                                if(hoverFrames.length<24) requestAnimationFrame(sample);
                              };
                              sample();
                            }""")
                            page.mouse.move(
                                initial_target["x"] + initial_target["width"] / 2,
                                initial_target["y"] + initial_target["height"] / 2,
                            )
                            owner.wait_for_function("hoverFrames.length === 24")
                            frames = owner.evaluate("hoverFrames")
                            widths = [frame["width"] for frame in frames]
                            assert widths == [widths[0]] * len(widths), widths
                            assert host.evaluate("e => e.clientWidth") == initial_width
                            assert (
                                host.evaluate("e => e.scrollHeight") == initial_height
                            )
                            assert all(
                                abs(frame["target"][key] - value) < 0.05
                                for frame in frames
                                for key, value in local_target.items()
                            ), (mode, theme, width, initial_target, frames)
                            if mode == "opaque-iframe":
                                assert initial_width == width, (
                                    "The frame cannot reach the host"
                                )
                            elif args.engine == "chromium":
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
                            page.screenshot(
                                path=str(
                                    args.output / f"{mode}-{theme}-{width}-hover.png"
                                )
                            )
                            target.dispatch_event("click")
                            page.mouse.move(1090, 890)
                            assert "Pinned" in panel.inner_text()
                            assert (
                                host.evaluate("e => e.scrollHeight") == initial_height
                            )
                            popup_box = panel.bounding_box()
                            assert popup_box["width"] > 250, popup_box
                            assert popup_box["x"] >= 0 and popup_box["y"] >= 0
                            assert popup_box["x"] + popup_box["width"] <= 1100
                            assert popup_box["y"] + popup_box["height"] <= 900
                            popup_owner.get_by_role(
                                "button", name="Close graph details"
                            ).dispatch_event("click")
                            panel.wait_for(state="hidden")
                            owner.evaluate(
                                "() => new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)))"
                            )
                            # Closing a pinned inspector can remove an outer
                            # scrollbar across a postMessage/ResizeObserver turn.
                            owner.wait_for_function(
                                """expected => {
                              const box=document.querySelector('[data-linnet-kind="node"]').getBoundingClientRect();
                              return Object.entries(expected).every(([key,value])=>Math.abs(box[key]-value)<.05);
                            }""",
                                arg=local_target,
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
                            # Hover a real sampled half-edge region, not a synthetic
                            # event: the visible highlight must stay under the cursor.
                            edge_point = owner.evaluate("""() => {
                              for (const target of document.querySelectorAll(
                                '[data-linnet-kind="halfedge"]:not([data-linnet-carrier])')) {
                                const r=target.getBoundingClientRect();
                                const x=r.x+r.width/2, y=r.y+r.height/2;
                                const hit=document.elementFromPoint(x,y)?.closest('[data-linnet-kind]');
                                if (hit === target) return {x,y,id:target.dataset.linnetId};
                              }
                            }""")
                            assert edge_point, "A real half-edge hit region is required"
                            edge = owner.locator(
                                f'[data-linnet-kind="halfedge"][data-linnet-id="{edge_point["id"]}"]:not([data-linnet-carrier])'
                            )
                            # Convert frame coordinates to page coordinates for the
                            # real pointer, including shadow-root frame wrappers.
                            offset = {"x": 0, "y": 0}
                            if mode != "inline":
                                offset = page.locator("#host iframe").bounding_box()
                            point = {
                                key: edge_point[key] + offset[key] for key in ("x", "y")
                            }
                            page.mouse.move(**point)
                            panel.wait_for(state="visible")
                            assert "Half-edge" in panel.inner_text()
                            shade = owner.locator(
                                f'.linnet-edge-highlight [data-linnet-shade-halfedge="{edge_point["id"]}"]'
                            )
                            assert shade.is_visible()
                            shade_box = shade.bounding_box()
                            assert all(
                                shade_box[key] - 12
                                <= point[key]
                                <= shade_box[key] + shade_box[size] + 12
                                for key, size in (("x", "width"), ("y", "height"))
                            ), (
                                mode,
                                point,
                                shade_box,
                                shade.evaluate("e=>e.outerHTML"),
                            )
                            popup_box = panel.bounding_box()
                            assert not (
                                popup_box["x"]
                                <= point["x"]
                                <= popup_box["x"] + popup_box["width"]
                                and popup_box["y"]
                                <= point["y"]
                                <= popup_box["y"] + popup_box["height"]
                            )
                            page.mouse.click(**point)
                            assert "Pinned" in panel.inner_text()
                            assert shade.is_visible(), (
                                "Pinning must preserve the hover highlight"
                            )
                            assert (
                                host.evaluate("e => e.scrollHeight") == initial_height
                            )
                            page.screenshot(
                                path=str(
                                    args.output / f"{mode}-{theme}-{width}-edge.png"
                                )
                            )
                            # Same-origin notebook frames share the outer top layer;
                            # a clipped host must not clip the inspector.
                            if mode != "opaque-iframe":
                                host.evaluate("e => e.style.maxHeight='50px'")
                                assert panel.is_visible()
                                assert (
                                    panel.bounding_box()["y"]
                                    + panel.bounding_box()["height"]
                                    > host.bounding_box()["y"] + 50
                                )
                            popup_owner.get_by_role(
                                "button", name="Close graph details"
                            ).click()
                            panel.wait_for(state="hidden")
                            # Keyboard inspection and dismissal use the same popup.
                            page.mouse.move(1090, 890)
                            edge.first.focus()
                            edge.first.press("Enter")
                            assert "Pinned" in panel.inner_text()
                            popup_owner.get_by_role(
                                "button", name="Close graph details"
                            ).focus()
                            page.keyboard.press("Escape")
                            panel.wait_for(state="hidden")
                            # Removing/replacing an output must remove its portal,
                            # including when the graph lived in a child document.
                            if mode != "opaque-iframe":
                                edge.first.press("Enter")
                                assert "Pinned" in panel.inner_text()
                                host.evaluate("e => e.remove()")
                                panel.wait_for(state="hidden")
                            assert not errors, errors
                            name = f"{mode}-{theme}-{width}"
                            page.screenshot(path=str(args.output / f"{name}.png"))
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
                # Two outputs share one floating inspector, but keep separate
                # selections. The portal also follows a live system-theme change.
                page = browser.new_page(
                    viewport={"width": 1100, "height": 1000}, color_scheme="light"
                )
                page.set_content(
                    '<body><div id="first" style="width:500px">'
                    + html
                    + '</div><div id="second" style="width:500px">'
                    + html
                    + "</div></body>"
                )
                first = page.locator('#first [data-linnet-kind="node"]').first
                second = page.locator('#second [data-linnet-kind="node"]').first
                panel = page.locator(".linnet-inspector")
                first.click(modifiers=["Shift"])
                assert page.locator(".linnet-inspector-construction").is_visible()
                second.hover()
                assert panel.count() == 1
                assert not page.locator(".linnet-inspector-construction").count()
                assert (
                    page.locator("#first svg[data-linnet-interactive]").evaluate(
                        "s => s.linnetSelection.nodes.length"
                    )
                    == 1
                )
                assert (
                    page.locator("#second svg[data-linnet-interactive]").evaluate(
                        "s => s.linnetSelection.nodes.length"
                    )
                    == 0
                )
                for theme, color in (
                    ("dark", "rgb(37, 33, 43)"),
                    ("light", "rgb(255, 253, 249)"),
                ):
                    page.emulate_media(color_scheme=theme)
                    page.wait_for_function(
                        """color => getComputedStyle(document.querySelector('.linnet-inspector-popup')
                      .shadowRoot.querySelector('.linnet-inspector')).backgroundColor === color""",
                        arg=color,
                    )
                page.locator("#second").evaluate("e => e.remove()")
                panel.wait_for(state="hidden")
                report["shared_popup_and_live_theme"] = True
                page.close()
            finally:
                browser.close()
        report["passed"] = True
    finally:
        (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
