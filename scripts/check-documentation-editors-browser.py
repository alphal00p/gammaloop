"""Exercise the actual standalone, embedded and full Marimo source editors."""

import argparse
import functools
import http.server
import json
import threading
from pathlib import Path

from playwright.sync_api import expect, sync_playwright

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--site", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument(
    "--skip-full-app",
    action="store_true",
    help="Omit the separate full app; still check standalone and embedded editors",
)
browsers = parser.add_mutually_exclusive_group(required=True)
browsers.add_argument("--chromium", type=Path)
browsers.add_argument("--webkit", type=Path)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass


server = http.server.ThreadingHTTPServer(
    ("127.0.0.1", 0), functools.partial(Quiet, directory=str(args.site.resolve()))
)
threading.Thread(target=server.serve_forever, daemon=True).start()
report = {"editors": [], "layout": [], "console": []}
try:
    with sync_playwright() as playwright:
        browser = (
            playwright.webkit.launch(executable_path=str(args.webkit), headless=True)
            if args.webkit
            else playwright.chromium.launch(
                executable_path=str(args.chromium), headless=True, args=["--no-sandbox"]
            )
        )
        # The lazy-import routing probe must reach Playwright rather than an
        # already controlling cache worker. Cache upgrades have a separate suite.
        page = browser.new_page(
            viewport={"width": 1440, "height": 1000}, service_workers="block"
        )
        page.on(
            "console",
            lambda message: (
                report["console"].append(message.text)
                if message.type == "error"
                else None
            ),
        )
        base = f"http://127.0.0.1:{server.server_port}/"

        def replace_code(content, source):
            if not args.webkit:
                content.fill(source)
                return
            # WebKit's contenteditable fill does not update CodeMirror reliably.
            # Use the owning view; completion picks and hover still use real UI.
            content.evaluate(
                """(element, source) => {
                const view = window.checkedCM.EditorView.findFromDOM(element.closest('.cm-editor'));
                view.dispatch({changes: {from: 0, to: view.state.doc.length, insert: source}, selection: {anchor: source.length}});
                view.focus();
            }""",
                source,
            )

        def type_interactions(owner, content, name):
            print("Editing", name, flush=True)
            original = content.inner_text()
            identity = owner.evaluate(
                """async (element, full) => {
                  const configURL = new URL('/' + location.pathname.split('/').slice(1,4).join('/') + '/assets/documentation-cells-config.json', location.href);
                  const config = await fetch(configURL).then(response=>response.json());
                  const module = await import(new URL('documentation-editor.js',new URL(config.entry,configURL)));
                  const host = element.getRootNode().host;
                  const runtime = full ? new URL(document.documentElement.dataset.documentationNotebookRuntime) :
                    new URL(host?.closest('.live-notebook')?.dataset.notebookRuntime || config.runtime,configURL);
                  window.checkedCM = await module.load(runtime,full?'full':'islands');
                  const dom = element.matches('.cm-editor')?element:element.querySelector('.cm-editor');
                  window.checkedView = window.checkedCM.EditorView.findFromDOM(dom);
                  return !!window.checkedView;
                }""",
                name == "full Marimo",
            )
            assert identity
            replace_code(
                content, "from symbolica.community.graph import Graph\nGraph.from_"
            )
            page.keyboard.press("Control+Space")
            tooltip = owner.locator(".cm-tooltip-autocomplete")
            expect(tooltip).to_be_visible(timeout=60000)
            # Both native completion surfaces must stay themed while open,
            # including at phone width with a long real API signature.
            replace_code(
                content,
                "from symbolica.community.graph import LayoutStream\nLayoutStream.from_",
            )
            for width in (1440, 390):
                page.set_viewport_size({"width": width, "height": 800})
                content.scroll_into_view_if_needed()
                for theme in ("light", "dark"):
                    print("Completion", name, width, theme, flush=True)
                    content.focus()
                    page.keyboard.press("Escape")
                    page.keyboard.press("Control+Space")
                    expect(tooltip.locator(".cm-completionInfo")).to_be_visible(
                        timeout=30000
                    )
                    page.evaluate(
                        "theme=>document.documentElement.dataset.theme=theme", theme
                    )
                    expect(tooltip).to_be_visible()
                    info = tooltip.locator(".cm-completionInfo")
                    expect(info).to_be_visible(timeout=30000)
                    # A theme update can reconfigure Marimo's view. Wait for
                    # its open completion to settle before reading both panels.
                    expect(tooltip).to_be_visible()
                    page.screenshot(
                        path=str(
                            args.output
                            / f"{name.replace(' ', '-')}-{width}-{theme}-completion.png"
                        )
                    )
                    for popup in (tooltip, info):
                        expect(popup).to_have_css("border-radius", "0px")
                        metrics = popup.evaluate(
                            """e=>{const r=e.getBoundingClientRect(),s=getComputedStyle(e);return {x:r.x,y:r.y,right:r.right,bottom:r.bottom,bg:s.backgroundColor,radius:s.borderRadius}}"""
                        )
                        assert metrics["x"] >= 0 and metrics["y"] >= 0, metrics
                        assert metrics["right"] <= width + 1, metrics
                        footer = page.locator(".notebook-toolbar")
                        bottom = footer.bounding_box()["y"] if footer.count() else 800
                        assert metrics["bottom"] <= bottom + 1, metrics
                        assert metrics["radius"] == "0px", metrics
                        assert metrics["bg"] == (
                            "rgb(45, 36, 48)"
                            if theme == "dark"
                            else "rgb(249, 246, 240)"
                        ), metrics
                    selected = tooltip.locator('[aria-selected="true"]')
                    assert selected.evaluate(
                        "e=>getComputedStyle(e).backgroundColor"
                    ) == (
                        "rgb(74, 55, 80)" if theme == "dark" else "rgb(234, 223, 240)"
                    )
                    assert page.evaluate(
                        "document.documentElement.scrollWidth <= innerWidth"
                    ), page.evaluate("""() => {
                      const roots=[document,...[...document.querySelectorAll('*')].map(e=>e.shadowRoot).filter(Boolean)];
                      return {width:innerWidth,scroll:document.documentElement.scrollWidth,
                        overflow:roots.flatMap(root=>[...root.querySelectorAll('*')].filter(e=>e.scrollWidth>e.clientWidth+1).map(e=>({tag:e.tagName,cls:e.className,width:e.clientWidth,scroll:e.scrollWidth,right:e.getBoundingClientRect().right,overflow:getComputedStyle(e).overflow,position:getComputedStyle(e).position}))).slice(-30)};
                    }""")
                    assert (
                        info.locator(
                            ".code-frame, .code-line-numbers, [data-copy-target]"
                        ).count()
                        == 0
                    )
                    page.screenshot(
                        path=str(
                            args.output
                            / f"{name.replace(' ', '-')}-{width}-{theme}-completion.png"
                        )
                    )
            page.set_viewport_size({"width": 1440, "height": 1000})
            content.focus()
            page.keyboard.press("Control+Space")
            expect(tooltip).to_be_visible()
            # The pinned completion keymap ignores acceptance for 75ms after
            # opening, to protect against an Enter intended for ordinary typing.
            page.wait_for_timeout(100)
            page.keyboard.press("ArrowDown")
            page.keyboard.press("ArrowUp")
            page.keyboard.press("Enter")
            expect(content).to_contain_text("LayoutStream.from_dot")
            replace_code(
                content, "from symbolica.community.graph import Graph\nGraph.from_"
            )
            page.keyboard.press("Control+Space")
            expect(tooltip).to_be_visible(timeout=30000)
            page.screenshot(
                path=str(args.output / (name.replace(" ", "-") + "-completion.png"))
            )
            tooltip.get_by_role("option").filter(has_text="from_dot").first.click()
            expect(content).to_contain_text("Graph.from_dot")
            page.keyboard.type("(")
            signature = page.locator(".documentation-signature")
            expect(signature).to_be_visible(timeout=30000)
            assert "source: str" in signature.inner_text()
            signature_bounds = signature.bounding_box()
            assert (
                signature_bounds
                and signature_bounds["x"] >= 0
                and signature_bounds["y"] >= 0
            )
            assert (
                signature_bounds["x"] + signature_bounds["width"]
                <= page.viewport_size["width"]
            )
            assert (
                signature_bounds["y"] + signature_bounds["height"]
                <= page.viewport_size["height"]
            )
            page.keyboard.press("Escape")
            expect(signature).to_have_count(0)
            for theme in ("light", "dark"):
                page.evaluate(
                    "theme=>document.documentElement.dataset.theme=theme", theme
                )
                replace_code(
                    content, "from symbolica.community.graph import Graph\nGraph"
                )
                content.locator(".cm-line").last.locator(".tok-variableName").hover(
                    position={"x": 15, "y": 8}
                )
                hover = owner.locator(".cm-tooltip-hover")
                expect(hover).to_be_visible(timeout=30000)
                assert "Graph" in hover.inner_text()
                info = hover.locator(".documentation-type-info")
                assert "```" not in info.inner_text()
                assert info.locator("pre > code").count() > 0
                assert (
                    info.locator(
                        ".code-frame, .code-line-numbers, [data-copy-target]"
                    ).count()
                    == 0
                )
                assert "Atkinson Hyperlegible" in info.evaluate(
                    "el=>getComputedStyle(el).fontFamily"
                )
                assert "IBM Plex Mono" in info.locator("code").first.evaluate(
                    "el=>getComputedStyle(el).fontFamily"
                )
                assert owner.evaluate(r"""() => {
                  const html = window.checkedCM.renderMarkdown('**Safe**\n\n```python\nx = 1\n```\n\n---\n<script>window.unsafeHover = true</script><img src=x onerror=alert(1)>[bad](javascript:alert(1))');
                  const root = document.createElement('div'); root.innerHTML = html;
                  return !!root.querySelector('strong') && !!root.querySelector('pre code') && !!root.querySelector('hr')
                    && !root.querySelector('script, [onerror], [href^="javascript:"]');
                }""")
                page.evaluate(
                    "() => new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))"
                )
                bounds = info.bounding_box()
                viewport = page.viewport_size
                assert bounds and bounds["width"] <= viewport["width"] - 16
                assert bounds["height"] <= viewport["height"] - 16
                assert bounds["x"] >= 0 and bounds["y"] >= 0
                assert bounds["x"] + bounds["width"] <= viewport["width"]
                assert bounds["y"] + bounds["height"] <= viewport["height"] + 1, (
                    name,
                    theme,
                    bounds,
                    viewport,
                )
                page.screenshot(
                    path=str(
                        args.output
                        / (name.replace(" ", "-") + "-" + theme + "-hover.png")
                    )
                )
            assert owner.evaluate(
                "e=>window.checkedCM.EditorView.findFromDOM(e.matches('.cm-editor')?e:e.querySelector('.cm-editor'))===window.checkedView"
            )
            replace_code(content, original)
            print("Passed", name, flush=True)
            report["editors"].append(
                {
                    "name": name,
                    "completionPick": True,
                    "signature": True,
                    "hover": True,
                    "sameView": True,
                }
            )

        page.goto(base + "products/linnet/latest/quickstart/python/")
        cell = page.locator(".documentation-cell").first
        expect(cell).to_be_visible(timeout=30000)
        started = []
        wheel_requests = []
        page.on(
            "request",
            lambda request: (
                wheel_requests.append(request.url) if ".whl" in request.url else None
            ),
        )
        page.on("worker", lambda worker: started.append(worker.url))
        # Hold the first lazy import so repeated clicks exercise initialization.
        imports = []
        toggle = cell.locator(".notebook-source-toggle")

        def hold_import(route):
            imports.append(route)
            toggle.evaluate("button=>button.dataset.importHeld='true'")

        page.route("**/dist-Kxgjxakw.js", hold_import)
        with page.expect_request("**/dist-Kxgjxakw.js"):
            toggle.click()
        expect(toggle).to_have_attribute("data-import-held", "true")
        expect(toggle).to_be_disabled()
        toggle.evaluate("button => { button.click(); button.click(); }")
        # WebKit can fetch one shared module twice while loading its import graph.
        # Keep the single-editor assertion below as the initialization invariant.
        assert len({route.request.url for route in imports}) == 1
        if not args.webkit:
            assert len(imports) == 1
        report["lazyImportRequests"] = len(imports)
        for route in imports:
            route.continue_()
        page.unroute("**/dist-Kxgjxakw.js")
        content = cell.locator(".cm-content")
        expect(content).to_be_visible(timeout=30000)
        expect(toggle).to_be_enabled()
        expect(cell.locator(".cm-editor")).to_have_count(1)
        type_interactions(cell, content, "standalone")
        assert not any("documentation-python-worker" in url for url in started)
        replace_code(content, "print(6 * 7)")
        cell.get_by_role("button", name="Run example", exact=True).click()
        reset = cell.get_by_role("button", name="Reset", exact=True)
        expect(reset).to_be_disabled()
        reset.evaluate("button => button.click()")
        expect(content).to_have_text("print(6 * 7)")
        expect(cell).to_have_attribute("data-cell-state", "success", timeout=120000)
        expect(reset).to_be_enabled()
        expect(cell.locator(".documentation-cell-output")).to_contain_text("42")
        installed = len(wheel_requests)
        assert installed > 0
        content.click()
        page.keyboard.press("Control+z")
        assert content.inner_text() != "print(6 * 7)"
        cell.locator(".notebook-source-toggle").click()
        expect(content).to_be_hidden()
        cell.locator(".notebook-source-toggle").click()
        expect(content).to_be_visible()
        for width in [1440, 832, 390, 360]:
            page.set_viewport_size({"width": width, "height": 1000})
            for theme in ["light", "dark"]:
                page.evaluate(
                    "theme=>document.documentElement.dataset.theme=theme", theme
                )
                assert page.evaluate(
                    "document.documentElement.scrollWidth <= innerWidth"
                )
                metrics = cell.locator(".cm-editor").evaluate("""editor=>{
                  const style=getComputedStyle(editor.querySelector('.cm-scroller'));
                  const gutter=editor.querySelector('.cm-gutters').getBoundingClientRect();
                  const line=editor.querySelector('.cm-line').getBoundingClientRect();
                  const frame=editor.closest('.reference-copy-row').getBoundingClientRect();
                  return {font:style.fontFamily,size:style.fontSize,leading:style.lineHeight,
                    outset:gutter.right<=frame.left+1,lineInset:line.left-frame.left};
                }""")
                assert "IBM Plex Mono" in metrics["font"]
                assert metrics["size"] == "14.4px" and metrics["outset"], metrics
                report["layout"].append({"width": width, "theme": theme, **metrics})
                page.screenshot(
                    path=str(args.output / f"standalone-{width}-{theme}.png")
                )
        page.set_viewport_size({"width": 1440, "height": 1000})
        cell.get_by_role("button", name="Reset", exact=True).click()
        graph = page.locator(".documentation-cell").nth(1)
        graph.get_by_role("button", name="Run example", exact=True).click()
        expect(graph).to_have_attribute("data-cell-state", "success", timeout=120000)
        assert len(wheel_requests) == installed
        report["warmReuse"] = {"wheelRequestsOnSecondRun": 0}
        drawing = graph.frame_locator("iframe").locator("svg[data-linnet-interactive]")
        expect(drawing).to_be_visible(timeout=30000)
        target = drawing.get_by_role("button", name="Edge 0", exact=True).first
        target.focus()
        expect(drawing.locator("[data-linnet-active]").first).to_be_visible()
        target.press("Shift+Space")
        expect(target).to_have_attribute("aria-pressed", "true")
        for theme in ["light", "dark"]:
            page.evaluate("theme=>document.documentElement.dataset.theme=theme", theme)
            expect(drawing).to_have_attribute("data-theme", theme)
        page.screenshot(path=str(args.output / "direct-graph.png"))
        # Exercise the LaTeX-only MIME path through the real MiTeX compiler.
        replace_code(
            content,
            "from symbolica.community import hepkit as hep\nhep.ThreeMomentum(1.0, 0.0, 0.0).cross(hep.ThreeMomentum(0.0, 1.0, 0.0))",
        )
        cell.get_by_role("button", name="Run example", exact=True).click()
        expect(cell).to_have_attribute("data-cell-state", "success", timeout=120000)
        expect(cell.locator(".documentation-cell-output img")).to_be_visible()
        cell.locator(".documentation-cell-output img").evaluate(
            "image => image.decode()"
        )
        page.screenshot(path=str(args.output / "direct-math.png"))
        page.goto(base + "products/linnet/latest/quickstart/typst/")
        typst = page.locator(".documentation-cell").first
        expect(typst).to_be_visible(timeout=30000)
        static = page.locator(
            'pre:has(code[data-lang="sh"]),pre:has(code[data-lang="shell"])'
        ).first
        assert (
            static.locator(
                'xpath=ancestor::*[contains(@class,"documentation-cell")]'
            ).count()
            == 0
        )
        typst.get_by_role("button", name="Compile example", exact=True).click()
        expect(typst).to_have_attribute("data-cell-state", "success", timeout=120000)
        expect(typst.locator(".documentation-cell-output img")).to_be_visible()
        typst.locator(".documentation-cell-output img").evaluate(
            "image => image.decode()"
        )
        typst.locator(".notebook-source-toggle").click()
        source = typst.locator(".cm-content")
        expect(source.locator(".tok-keyword").first).to_be_visible()
        source.fill("#let =")
        typst.get_by_role("button", name="Compile example", exact=True).click()
        expect(typst).to_have_attribute("data-cell-state", "error", timeout=30000)
        expect(typst.get_by_role("button", name="Reset", exact=True)).to_be_enabled()
        typst.get_by_role("button", name="Reset", exact=True).click()
        typst.get_by_role("button", name="Compile example", exact=True).click()
        expect(typst).to_have_attribute("data-cell-state", "success", timeout=30000)
        typst.locator(".documentation-cell-output img").evaluate(
            "image => image.decode()"
        )
        page.screenshot(path=str(args.output / "direct-typst.png"))
        print("Passed direct graph, math, Typst and recovery", flush=True)
        report["richOutput"] = {
            "graphInteraction": True,
            "graphTheme": True,
            "latex": True,
            "typstCompile": True,
            "typstErrorRecovery": True,
            "staticHasNoRun": True,
        }
        hydration = []
        page.route(
            "**/notebook-runtime/**/main.js", lambda route: hydration.append(route)
        )
        with page.expect_request("**/notebook-runtime/**/main.js", timeout=30000):
            page.goto(base + "products/feynkit/latest/guides/showcases/first-diagram/")
            notebook = page.locator(".live-notebook")
            notebook.scroll_into_view_if_needed()
        expect(
            notebook.locator(".live-notebook-code-preview:visible").first
        ).to_be_visible(timeout=30000)
        page.wait_for_function(
            "document.querySelector('.live-notebook')?.dataset.notebookState === 'loading'"
        )
        # The held import proves this is the fallback during runtime hydration.
        for route in hydration:
            route.continue_()
        page.unroute("**/notebook-runtime/**/main.js")
        expect(notebook).to_have_attribute(
            "data-notebook-state", "active", timeout=120000
        )
        wrapper = page.locator(
            ".notebook-cell:not([hidden]):has(marimo-code-editor[data-notebook-source])"
        ).first
        toggle = wrapper.locator(".notebook-source-toggle")
        if toggle.get_attribute("aria-expanded") != "true":
            toggle.click()
        native = wrapper.locator("marimo-code-editor[data-notebook-source]")
        content = native.locator(".cm-content")
        expect(content).to_be_visible(timeout=30000)
        # Playwright locators cross this editor's shadow root. Evaluate its own
        # root explicitly when checking the native view identity.
        owner = native.locator(".cm-editor")
        type_interactions(owner, content, "embedded Marimo")
        styles = page.locator('link[rel="stylesheet"][title]').evaluate_all(
            "links=>links.map(link=>({title:link.title,loaded:!!link.sheet,disabled:link.sheet?.disabled}))"
        )
        assert {sheet["title"] for sheet in styles} == {"marimo"}, styles
        assert all(sheet["loaded"] and not sheet["disabled"] for sheet in styles), (
            styles
        )
        report["stylesheetSet"] = styles
        content.click()
        page.keyboard.press("Control+k")
        expect(page.locator("[data-search-dialog]")).to_be_visible()
        page.keyboard.press("Escape")
        expect(content).to_be_focused()
        page.screenshot(path=str(args.output / "embedded.png"))
        # ty inserts a notice between the native editor and its loading source.
        # Visibility must follow the hydrated editor even after that insertion.
        expect(
            native.locator(
                "xpath=following-sibling::span[contains(@class,'documentation-type-status')]"
            )
        ).to_have_count(1)
        preview = wrapper.locator(".live-notebook-code-preview")
        expect(preview).to_be_hidden()
        original = content.inner_text()
        replace_code(content, 'print("retained source")')
        selection = page.evaluate(
            "JSON.stringify(window.checkedView.state.selection.toJSON())"
        )
        toggle.click()
        expect(content).to_be_hidden()
        expect(preview).to_be_hidden()
        toggle.click()
        expect(content).to_have_text('print("retained source")')
        assert owner.evaluate(
            "e=>window.checkedCM.EditorView.findFromDOM(e)===window.checkedView"
        )
        assert (
            page.evaluate("JSON.stringify(window.checkedView.state.selection.toJSON())")
            == selection
        )
        content.click()
        page.keyboard.press("Control+z")
        assert content.inner_text() == original
        setup = notebook.locator("[data-notebook-setup]").first
        expect(setup).to_be_hidden()
        page.locator(".notebook-toolbar summary").click()
        page.get_by_role("button", name="Show setup cells", exact=True).click()
        expect(setup).to_be_visible()
        expect(setup.locator(".live-notebook-code-preview")).to_be_hidden()
        page.get_by_role("button", name="Hide setup cells", exact=True).click()
        expect(setup).to_be_hidden()
        report["sourceVisibility"] = {
            "afterTyNotice": True,
            "togglePreservesUndo": True,
            "setupReveal": True,
        }

        def markdown_copy():
            frame = (
                page.locator(".codehilite-wrapper:has(.codehilite.code-frame)")
                .filter(has_text="ufo-model-loader==1.0.1")
                .first
            )
            expect(frame).to_be_visible(timeout=180000)
            button = frame.get_by_role("button", name="Copy to clipboard", exact=True)
            expect(button).to_have_count(1)
            expected = frame.locator("pre code").inner_text()
            for width in [1440, 390]:
                page.set_viewport_size({"width": width, "height": 1000})
                for theme in ["light", "dark"]:
                    if "/full/" in page.url:
                        page.locator('button[aria-label="Config"]:visible').last.click()
                        page.get_by_text("User settings", exact=True).click()
                        page.get_by_text("Display", exact=True).last.click()
                        page.get_by_test_id("theme-select").select_option(theme)
                        page.get_by_role(
                            "button", name="Close", exact=True
                        ).last.click()
                        page.keyboard.press("Escape")
                        expect(page.locator("html")).to_have_attribute(
                            "data-code-theme", theme
                        )
                    else:
                        page.evaluate(
                            "theme=>document.documentElement.dataset.theme=theme", theme
                        )
                    frame.scroll_into_view_if_needed()
                    metrics = frame.evaluate("""e=>{
                      const f=e.querySelector('.code-frame').getBoundingClientRect();
                      const b=e.querySelector('button').getBoundingClientRect();
                      const p=e.querySelector('pre');
                      const g=e.querySelector('.code-line-numbers');
                      return {bounded:b.left>=f.left && b.right<=f.right && b.top>=f.top && b.bottom<=f.bottom,
                        aligned:Math.abs(p.getBoundingClientRect().top-g.getBoundingClientRect().top)<1,
                        padding:getComputedStyle(p).paddingTop,gutterPadding:getComputedStyle(g).paddingTop};
                    }""")
                    assert (
                        metrics["bounded"]
                        and metrics["aligned"]
                        and metrics["padding"] == metrics["gutterPadding"] == "12px"
                    ), metrics
                    report["layout"].append(
                        {
                            "markdown": page.url,
                            "width": width,
                            "theme": theme,
                            **metrics,
                        }
                    )
                    page.screenshot(
                        path=str(
                            args.output
                            / f"markdown-{'full' if '/full/' in page.url else 'embedded'}-{width}-{theme}.png"
                        )
                    )
            button.click()
            assert page.evaluate("navigator.clipboard.readText()") == expected
            expect(button.locator("svg.lucide-check")).to_be_visible()
            page.set_viewport_size({"width": 1440, "height": 1000})

        if args.webkit:
            # Playwright WebKit cannot grant these permissions. Check the copy
            # API contract without claiming coverage of the system clipboard.
            page.add_init_script("""Object.defineProperty(navigator, 'clipboard', {value: {
              writeText: async text => { window.copiedText = text; },
              readText: async () => window.copiedText,
            }});""")
        else:
            page.context.grant_permissions(["clipboard-read", "clipboard-write"])
        page.goto(base + "products/feynkit/latest/guides/showcases/ufo/")
        notebook = page.locator(".live-notebook")
        notebook.scroll_into_view_if_needed()
        expect(notebook).to_have_attribute(
            "data-notebook-state", "active", timeout=120000
        )
        markdown_copy()
        assert not page.evaluate("""() => performance.getEntriesByType('resource')
          .some(entry => new URL(entry.name).pathname.endsWith('/cells.json'))""")
        for width in [1440, 832, 390, 360]:
            page.set_viewport_size({"width": width, "height": 1000})
            for theme in ["dark", "light", "dark"]:
                page.evaluate(
                    "theme=>document.documentElement.dataset.theme=theme", theme
                )
                metrics = page.evaluate("""() => {
                  const bar=document.querySelector('.notebook-toolbar').getBoundingClientRect();
                  const codes=[...document.querySelectorAll('.live-notebook .prose code')]
                    .filter(code=>!code.closest('pre'));
                  const probe=document.createElement('span');
                  probe.style.color='var(--ink)';document.body.append(probe);
                  const ink=getComputedStyle(probe).color;probe.remove();
                  return {left:bar.left,right:bar.right,bottom:bar.bottom,
                    width:innerWidth,height:innerHeight,
                    codes:codes.length,readable:codes.every(code=>getComputedStyle(code).color===ink)};
                }""")
                assert (
                    abs(metrics["left"]) < 1
                    and abs(metrics["right"] - width) < 1
                    and abs(metrics["bottom"] - 1000) < 1
                    and metrics["codes"]
                    and metrics["readable"]
                ), metrics
                report["layout"].append({"footer": True, "theme": theme, **metrics})
        page.set_viewport_size({"width": 1440, "height": 1000})
        page.goto(
            base + "products/feynkit/latest/guides/showcases/four-loop-numerator/"
        )
        page.locator(".live-notebook").scroll_into_view_if_needed()
        dropdown = page.locator("marimo-dropdown").first
        trigger = dropdown.get_by_role("button")
        expect(trigger).to_be_visible(timeout=180000)
        report["dropdowns"] = []
        for width in (1440, 390):
            page.set_viewport_size({"width": width, "height": 800})
            trigger.click()
            menu = page.locator('[role="dialog"]:has([cmdk-root])')
            search = menu.get_by_role("combobox")
            expect(search).to_be_focused()
            for theme in ("light", "dark", "light"):
                page.evaluate(
                    "theme=>document.documentElement.dataset.theme=theme", theme
                )
                expect(dropdown.locator(".marimo > .contents")).to_have_attribute(
                    "data-theme", theme
                )
                expect(menu).to_have_css(
                    "background-color",
                    "rgb(45, 36, 48)" if theme == "dark" else "rgb(249, 246, 240)",
                )
                expect(menu).to_have_css("border-radius", "0px")
                expect(search).to_be_focused()
                bounds = menu.bounding_box()
                assert (
                    bounds["x"] >= 0 and bounds["x"] + bounds["width"] <= width + 1
                ), bounds
                assert (
                    bounds["y"] >= 0
                    and bounds["y"] + bounds["height"]
                    <= page.locator(".notebook-toolbar").bounding_box()["y"] + 1
                ), bounds
                page.screenshot(path=str(args.output / f"dropdown-{width}-{theme}.png"))
                report["dropdowns"].append(
                    {"width": width, "theme": theme, "bounds": bounds}
                )
            search.fill("no-such-topology")
            expect(menu.locator("[cmdk-empty]")).to_be_visible()
            search.fill("D0001")
            expect(menu.get_by_role("option")).to_have_count(1)
            page.keyboard.press("ArrowDown")
            page.keyboard.press("Enter")
            expect(menu).to_have_count(0)
            expect(trigger).to_be_focused()
            expect(trigger).to_contain_text("D0001")
            expect(
                page.frame_locator('iframe[title="Graph diagram"]').last.locator(
                    "figcaption"
                )
            ).to_contain_text("D0001", timeout=120000)
            trigger.click()
            page.keyboard.press("Escape")
            expect(trigger).to_be_focused()
            assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
        # React replaces the dropdown output on rerun. The new shadow boundary
        # and newly created portal must pick up the current theme again.
        page.evaluate("document.documentElement.dataset.theme='dark'")
        dropdown_cell = page.locator(".notebook-cell").filter(has=dropdown)
        dropdown_cell.get_by_role("button", name="Run cell", exact=True).first.click()
        expect(dropdown.locator(".marimo > .contents")).to_have_attribute(
            "data-theme", "dark"
        )
        trigger.click()
        expect(page.locator('[role="dialog"]:has([cmdk-root])')).to_have_css(
            "background-color", "rgb(45, 36, 48)"
        )
        page.keyboard.press("Escape")
        page.goto(
            base + "products/feynkit/latest/guides/showcases/numerical-integration/"
        )
        page.locator(".live-notebook").scroll_into_view_if_needed()
        native_select = page.locator("marimo-dropdown select").first
        expect(native_select).to_be_visible(timeout=180000)
        for theme in ("light", "dark"):
            page.evaluate("theme=>document.documentElement.dataset.theme=theme", theme)
            expect(native_select).to_have_css("color-scheme", theme)
            expect(native_select).to_have_css("border-radius", "0px")
        orientation = page.locator('iframe[title="Graph diagram"]').nth(4)
        old_drawing = orientation.get_attribute("srcdoc")
        native_select.select_option(label="Crossed double box · 2 loops")
        expect(native_select).to_have_value("Crossed double box · 2 loops")
        expect(orientation).not_to_have_attribute("srcdoc", old_drawing, timeout=120000)
        print("Passed searchable and native dropdowns", flush=True)
        page.set_viewport_size({"width": 1440, "height": 1000})
        if not args.skip_full_app:
            page.goto(
                base
                + "products/feynkit/latest/assets/notebooks/full/04_ufo_loading_marimo.html"
            )
            expect(
                page.locator("[data-documentation-type-ready]:visible").first
            ).to_be_visible(timeout=90000)
            page.locator("[data-testid=run-button]").last.click()
            markdown_copy()
        report["nativeMarkdownCopy"] = {
            "embedded": True,
            "full": not args.skip_full_app,
            "clipboard": "browser API stub" if args.webkit else True,
        }
        if not args.skip_full_app:
            page.goto(
                base
                + "products/feynkit/latest/assets/notebooks/full/00_quickstart_marimo.html",
                wait_until="domcontentloaded",
            )
            full = page.locator('[data-documentation-type-ready="true"]:visible').first
            expect(full).to_be_visible(timeout=90000)
            content = full.locator(".cm-content")
            type_interactions(full, content, "full Marimo")
            page.screenshot(path=str(args.output / "full.png"))
        assert not any(
            "Field is not present" in error for error in report["console"]
        ), report["console"]
        browser.close()
finally:
    server.shutdown()
    server.server_close()
    (args.output / "editors.json").write_text(json.dumps(report, indent=2))
print(
    "Editor interactions passed:",
    ", ".join(editor["name"] for editor in report["editors"]),
)
