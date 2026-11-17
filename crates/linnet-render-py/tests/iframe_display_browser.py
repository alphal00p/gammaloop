"""Exercise the actual Marimo iframe adapter with retained native graph captures.

Pass --html repeatedly to wrap raw captures with src/view.py, or --wrapper-html
to use installed DiagramView._display_().text captures without changing their
bootstrap. --source ROOT optionally refreshes raw captures' interactive JS/CSS.
This client never builds or imports Symbolica. MATHML_CHROMIUM_EXECUTABLE can
select Chromium; screenshots use CDP so pending fonts cannot delay the capture.
"""

import argparse
import base64
import hashlib
import importlib.util
import json
import re
from pathlib import Path

from mathml_view_browser import ROOT, launch
from playwright.sync_api import sync_playwright

INIT = r"""mode => {
  if(window===top) return;
  const state=window.lifecycle={token:Math.random(),settled:false,first:null,frames:[],timers:[]};
  new MutationObserver(()=>{
    if(state.first===null&&document.querySelector('svg'))
      state.first=getComputedStyle(frameElement).visibility;
  }).observe(document,{childList:true,subtree:true});
  if(['pending','hung','rejected'].includes(mode)) {
    const promise=mode==='rejected'?Promise.reject(new Error('expected font failure')):
      new Promise(resolve=>{if(mode==='pending')state.resolve=resolve;});
    promise.then(()=>state.settled=true,()=>state.settled=true);
    Object.defineProperty(document.fonts,'ready',{get:()=>promise});
  }
  if(mode==='raf'||mode==='stale') {
    window.requestAnimationFrame=callback=>state.frames.push(callback);
    if(mode==='stale') {
      const timer=window.setTimeout;
      window.setTimeout=(callback,delay,...args)=>delay===100&&callback.name==='reveal'
        ?state.timers.push(callback):timer(callback,delay,...args);
    }
  }
  if(mode==='network') document.addEventListener('DOMContentLoaded',()=>{
    document.fonts.load('21px "Linnet Offline Math"');
    document.fonts.ready.then(()=>state.settled=true);
  },{once:true});
}"""
STATE = """() => {
  const svg=document.querySelector('svg[data-linnet-interactive]');
  return {visible:getComputedStyle(frameElement).visibility,
    parent:getComputedStyle(frameElement.parentElement).visibility,
    ready:frameElement.hasAttribute('data-linnet-ready'),first:lifecycle.first,
    settled:lifecycle.settled,theme:svg.dataset.theme,
    controls:!!svg.querySelector('.linnet-toolbar'),
    zoom:svg.querySelector('.linnet-toolbar output')?.textContent,
    camera:svg.querySelector('.linnet-viewport')?.getAttribute('viewBox'),
    selection:svg.linnetSelection,height:frameElement.getBoundingClientRect().height,
    caption:(document.querySelector('figcaption')||svg).getBoundingClientRect().bottom};
}"""


def wrappers(args):
    if args.wrapper_html:
        assert not args.source, "--source is only valid for raw --html captures"
        return [(path, path.read_text()) for path in args.wrapper_html]
    spec = importlib.util.spec_from_file_location(
        "iframe_adapter", ROOT / "crates/linnet-render-py/src/view.py"
    )
    adapter = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adapter)
    bootstrap = (ROOT / "crates/linnet-render-py/src/iframe.js").read_text()
    result = []
    for path in args.html:
        html = path.read_text()
        if args.source:
            sources = args.source / "crates/linnest/src/svg"
            html, scripts = re.subn(
                r"<script><!\[CDATA\[.*?\]\]></script>",
                lambda _, sources=sources: (
                    "<script><![CDATA["
                    + (sources / "interactive.js").read_text()
                    + "]]></script>"
                ),
                html,
                flags=re.DOTALL,
            )
            html, styles = re.subn(
                r"<style>svg\[data-linnet-interactive\].*?</style>",
                lambda _, sources=sources: (
                    "<style>" + (sources / "interactive.css").read_text() + "</style>"
                ),
                html,
                flags=re.DOTALL,
            )
            assert scripts == styles == 1
        result.append((path, adapter.iframe(html, bootstrap).text))
    return result


def screenshot(page, locator, path):
    """Capture real paint without Playwright waiting for document.fonts.ready."""
    client = page.context.new_cdp_session(page)
    try:
        image = client.send(
            "Page.captureScreenshot",
            {"format": "png", "clip": {**locator.bounding_box(), "scale": 1}},
        )["data"]
    finally:
        client.detach()
    path.write_bytes(base64.b64decode(image))
    return page.evaluate(
        """async data=>{
          const image=await createImageBitmap(await (await fetch('data:image/png;base64,'+data)).blob());
          const canvas=document.createElement('canvas');canvas.width=image.width;canvas.height=image.height;
          const context=canvas.getContext('2d');context.drawImage(image,0,0);
          const pixels=context.getImageData(0,0,canvas.width,canvas.height).data;
          let ink=0;
          for(let i=4;i<pixels.length;i+=4)
            if(Math.max(...[0,1,2].map(c=>Math.abs(pixels[i+c]-pixels[c])))>60) ink++;
          return ink;
        }""",
        image,
    )


def mount(page, markup, *, error=False):
    return page.evaluate(
        r"""({markup,error})=>{
          const template=document.createElement('template');template.innerHTML=markup;
          const iframe=template.content.querySelector('iframe');
          if(error) iframe.srcdoc=iframe.srcdoc.replace(/<script>[\s\S]*?<\/script>/g,
            source=>source.includes('Scope state to the SVG')
              ?'<script>throw new Error("expected graph initialization failure")<\/script>':source);
          host.replaceChildren(template.content);
          return getComputedStyle(host.querySelector('iframe')).visibility;
        }""",
        {"markup": markup, "error": error},
    )


def ready(page):
    page.wait_for_function("host.querySelector('iframe[data-linnet-ready]')")
    return page.locator("#host iframe").element_handle().content_frame()


def settle(frame):
    frame.evaluate(
        "() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))"
    )


def assert_fits(frame, selector=None):
    state = frame.evaluate(STATE)
    assert state["height"] >= state["caption"] - 1, state
    if selector:
        bottom = frame.locator(selector).evaluate(
            "node=>node.getBoundingClientRect().bottom"
        )
        assert state["height"] >= bottom - 1, (state, selector, bottom)
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--html", type=Path, action="append")
    inputs.add_argument("--wrapper-html", type=Path, action="append")
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    captures = wrappers(args)
    report = {"ok": False, "installed_wrappers": bool(args.wrapper_html), "cases": []}
    report["inputs"] = {
        str(path): hashlib.sha256(markup.encode()).hexdigest()
        for path, markup in captures
    }
    with sync_playwright() as playwright:
        browser = launch(playwright, "chromium")
        try:
            for path, markup in captures:
                font_url = re.search(
                    r"data:font/woff2;base64,([A-Za-z0-9+/=]+)", markup
                )
                assert font_url, "Capture must retain its native embedded font"
                font = base64.b64decode(font_url.group(1))
                for theme in ("light", "dark"):
                    for mode in (
                        "pending",
                        "hung",
                        "rejected",
                        "network",
                        "raf",
                        "error",
                        "stale",
                    ):
                        context = browser.new_context(
                            viewport={"width": 1200, "height": 1000},
                            color_scheme="dark" if theme == "light" else "light",
                        )
                        held, errors = [], []
                        context.add_init_script(
                            "(" + INIT + ")(" + json.dumps(mode) + ");"
                        )
                        page = context.new_page()
                        page.set_default_timeout(5000)
                        page.on(
                            "pageerror",
                            lambda error, errors=errors: errors.append(str(error)),
                        )
                        page.route(
                            "**/font.woff2",
                            lambda route, _request, *, held=held: held.append(route),
                        )
                        background, color = (
                            ("#fff", "#000")
                            if theme == "light"
                            else ("#1c2025", "#e6ebf1")
                        )
                        page.set_content(
                            f'<body data-theme="{theme}" style="background:{background};color:{color}">'
                            '<div id="host" style="width:1060px"></div>'
                            "<script>window.__resizeIframe=()=>{};</script>"
                        )
                        item = {"capture": path.name, "theme": theme, "mode": mode}
                        report["cases"].append(item)
                        try:
                            source = (
                                markup.replace(
                                    font_url.group(0), "https://linnet.test/font.woff2"
                                )
                                if mode == "network"
                                else markup
                            )
                            for _ in range(
                                2
                                if mode in ("pending", "hung", "rejected", "raf")
                                else 1
                            ):
                                assert (
                                    mount(page, source, error=mode == "error")
                                    == "hidden"
                                )
                                if mode == "stale":
                                    page.wait_for_function(
                                        "host.querySelector('iframe').contentWindow.lifecycle?.timers.length"
                                    )
                                    page.evaluate(
                                        """()=>{const frame=host.querySelector('iframe');window.old=frame.contentWindow.lifecycle;frame.srcdoc=frame.srcdoc;}"""
                                    )
                                    page.wait_for_function(
                                        "host.querySelector('iframe').contentWindow.lifecycle?.token!==old.token&&host.querySelector('iframe').contentWindow.lifecycle?.timers.length"
                                    )
                                    assert page.evaluate(
                                        """()=>{old.frames.forEach(f=>f());old.timers.forEach(f=>f());const frame=host.querySelector('iframe');return getComputedStyle(frame).visibility==='hidden'&&!frame.hasAttribute('data-linnet-ready');}"""
                                    )
                                    page.evaluate(
                                        """()=>{const current=host.querySelector('iframe').contentWindow.lifecycle;current.frames.forEach(f=>f());current.timers.forEach(f=>f());}"""
                                    )
                                frame = ready(page)
                                state = frame.evaluate(STATE)
                                assert (
                                    state["visible"] == "visible"
                                    and state["parent"] == "hidden"
                                ), state
                                assert state["first"] == "hidden" and state["ready"], (
                                    state
                                )
                                assert state["height"] >= state["caption"] - 1, state
                                if mode != "error":
                                    assert (
                                        state["controls"] and state["theme"] == theme
                                    ), state
                                else:
                                    assert not state["controls"], state
                                if mode == "raf":
                                    assert frame.evaluate("lifecycle.frames.length > 0")
                                if mode in ("pending", "hung", "network"):
                                    assert not state["settled"], state
                                item["ready"] = state
                            if mode == "stale":
                                # Navigate the already visible element: the new bootstrap
                                # must hide it before any of the next document's SVG paints.
                                page.evaluate(
                                    """()=>{const frame=host.querySelector('iframe');window.visibleDocument=frame.contentWindow.lifecycle;frame.srcdoc=frame.srcdoc;}"""
                                )
                                page.wait_for_function(
                                    "host.querySelector('iframe').contentWindow.lifecycle?.token!==visibleDocument.token&&host.querySelector('iframe').contentWindow.lifecycle?.timers.length"
                                )
                                assert page.evaluate(
                                    """()=>{const frame=host.querySelector('iframe');const current=frame.contentWindow.lifecycle;return current.first==='hidden'&&getComputedStyle(frame).visibility==='hidden'&&!frame.hasAttribute('data-linnet-ready');}"""
                                )
                                page.evaluate(
                                    """()=>{const current=host.querySelector('iframe').contentWindow.lifecycle;current.frames.forEach(f=>f());current.timers.forEach(f=>f());}"""
                                )
                                frame = ready(page)
                                item["visible_navigation"] = assert_fits(frame)
                                assert (
                                    item["visible_navigation"]["visible"] == "visible"
                                )
                            if mode == "pending":
                                frame.evaluate("lifecycle.resolve()")
                                frame.wait_for_function("lifecycle.settled")
                            if mode == "network":
                                assert (
                                    held
                                    and frame.evaluate("document.fonts.status")
                                    == "loading"
                                )
                                prefix = args.output / f"{path.stem}-{theme}"
                                page.wait_for_timeout(200)
                                ink = screenshot(
                                    page,
                                    frame.locator("math").first,
                                    prefix.with_name(
                                        prefix.name + "-pending-label.png"
                                    ),
                                )
                                assert ink > 5, (
                                    "Ready graph labels remain invisible while the font loads"
                                )
                                item["pending_label_ink"] = ink
                                screenshot(
                                    page,
                                    page.locator("#host"),
                                    prefix.with_name(prefix.name + "-100.png"),
                                )
                                frame.locator(".linnet-navigation-help").focus()
                                assert frame.locator(
                                    ".linnet-navigation-hint"
                                ).is_visible()
                                settle(frame)
                                assert_fits(frame, ".linnet-navigation-hint")
                                screenshot(
                                    page,
                                    page.locator("#host"),
                                    prefix.with_name(prefix.name + "-help.png"),
                                )
                                frame.locator(
                                    '[data-linnet-kind="node"][tabindex="0"]'
                                ).first.press("Shift+Enter")
                                settle(frame)
                                item["inspector"] = assert_fits(frame)
                                popup = page.locator(".linnet-inspector").bounding_box()
                                assert popup["x"] >= 0 and popup["y"] >= 0, popup
                                assert (
                                    popup["x"] + popup["width"]
                                    <= page.viewport_size["width"]
                                ), popup
                                assert (
                                    popup["y"] + popup["height"]
                                    <= page.viewport_size["height"]
                                ), popup
                                screenshot(
                                    page,
                                    page.locator("#host"),
                                    prefix.with_name(prefix.name + "-inspector.png"),
                                )
                                unzoomed = frame.evaluate(STATE)
                                frame.locator(".linnet-viewport").press("+")
                                before = frame.evaluate(STATE)
                                assert (
                                    before["selection"]["nodes"]
                                    and before["camera"] != unzoomed["camera"]
                                ), before
                                for route in held:
                                    route.fulfill(body=font, content_type="font/woff2")
                                frame.wait_for_function("lifecycle.settled")
                                after = frame.evaluate(STATE)
                                assert all(
                                    before[key] == after[key]
                                    for key in ("zoom", "camera", "selection")
                                ), (before, after)
                                frame.locator(".linnet-viewport").press("0")
                                frame.locator(".linnet-viewport").press("Escape")
                                page.locator("#host").evaluate(
                                    "host=>host.style.width='320px'"
                                )
                                settle(frame)
                                narrow = assert_fits(frame)
                                screenshot(
                                    page,
                                    page.locator("#host"),
                                    prefix.with_name(prefix.name + "-narrow.png"),
                                )
                                item["narrow"] = narrow
                            assert errors == (
                                ["expected graph initialization failure"]
                                if mode == "error"
                                else []
                            ), errors
                            page.evaluate(
                                """()=>{const frame=host.querySelector('iframe');const callbacks=frame.contentWindow.lifecycle;host.replaceChildren();callbacks.frames.forEach(f=>f());callbacks.timers.forEach(f=>f());}"""
                            )
                            assert page.locator("#host iframe").count() == 0
                            item["ok"] = True
                        finally:
                            context.close()
                            (args.output / "report.json").write_text(
                                json.dumps(report, indent=2)
                            )
            report["ok"] = True
        finally:
            browser.close()
            (args.output / "report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({"ok": report["ok"], "cases": len(report["cases"])}))


if __name__ == "__main__":
    main()
