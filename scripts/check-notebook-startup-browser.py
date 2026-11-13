"""Exercise real notebook startup, an existing cache upgrade, and worker failures."""

import argparse
import functools
import http.server
import json
import threading
from pathlib import Path

from playwright.sync_api import expect, sync_playwright

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--site", required=True, type=Path)
parser.add_argument("--output", required=True, type=Path)
parser.add_argument("--previous-cache-worker", type=Path)
engines = parser.add_mutually_exclusive_group(required=True)
engines.add_argument("--webkit", type=Path)
engines.add_argument("--chromium", type=Path)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
mode = "cold"
route = "/products/feynkit/latest/guides/showcases/four-loop-numerator/"


class Handler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *_args):
        pass

    def do_GET(self):
        path = self.path.split("?")[0]
        if path.endswith("/assets/worker-D7hWi02M.js") and mode == "module-load":
            self.send_error(503)
            return
        if path.endswith("/assets/worker-D7hWi02M.js") and mode in {
            "worker-rejection",
            "dependency-rejection",
        }:
            data = Path(self.translate_path(path)).read_bytes()
            if mode == "worker-rejection":
                data = (
                    b"Promise.reject(new Error('injected startup rejection'));\n" + data
                )
            else:
                boundary = b"async loadNotebookDeps(e,t){"
                assert data.count(boundary) == 1
                data = data.replace(
                    boundary, boundary + b"throw Error('injected dependency failure');"
                )
            self.send_response(200)
            self.send_header("Content-Type", "text/javascript")
        elif path == "/notebook-startup-previous-worker.js":
            data = args.previous_cache_worker.read_bytes()
            self.send_response(200)
            self.send_header("Content-Type", "text/javascript")
        elif mode == "previous-cache" and path.endswith(
            "05_four_loop_numerator_marimo.json"
        ):
            payload = json.loads(Path(self.translate_path(path)).read_text())
            payload["cache_worker"] = "/notebook-startup-previous-worker.js"
            data = json.dumps(payload).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
        else:
            super().do_GET()
            return
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)


server = http.server.ThreadingHTTPServer(
    ("127.0.0.1", 0), functools.partial(Handler, directory=str(args.site.resolve()))
)
thread = threading.Thread(target=server.serve_forever, daemon=True)
thread.start()
results = []
try:
    with sync_playwright() as driver:
        browser = (
            driver.webkit.launch(executable_path=str(args.webkit))
            if args.webkit
            else driver.chromium.launch(
                executable_path=str(args.chromium), args=["--no-sandbox"]
            )
        )
        try:
            context = browser.new_context()
            trials = (
                ["previous-cache", "cache-upgrade", "reload"]
                if args.previous_cache_worker
                else ["cold", "reload"]
            )
            for mode in [
                *trials,
                "module-load",
                "worker-rejection",
                "dependency-rejection",
            ]:
                failing = mode in {
                    "module-load",
                    "worker-rejection",
                    "dependency-rejection",
                }
                current = browser.new_context() if failing else context
                page = current.new_page()
                errors = []
                page.on(
                    "pageerror", lambda error, errors=errors: errors.append(str(error))
                )
                page.add_init_script("""window.startupStates=[];
                  new MutationObserver(()=>{
                    const s=document.querySelector('.notebook-toolbar [role=status]');
                    if(s && window.startupStates.at(-1)?.text!==s.textContent)
                      window.startupStates.push({ms:performance.now(),text:s.textContent});
                  }).observe(document,{subtree:true,childList:true,characterData:true});""")
                page.goto(f"http://127.0.0.1:{server.server_port}{route}")
                page.locator("[data-notebook]").scroll_into_view_if_needed()
                page.wait_for_function(
                    "state=>document.querySelector('.notebook-toolbar [role=status]')?.dataset.state===state",
                    arg="error" if failing else "ready",
                    timeout=30000 if failing else 180000,
                )
                if not failing:
                    expect(page.locator("[data-notebook-ready]")).to_have_count(1)
                    expect(page.locator("[data-notebook-ready]")).to_contain_text(
                        "Generated 4,970 graphs"
                    )
                    assert not errors, errors
                    page.wait_for_function("navigator.serviceWorker.controller!==null")
                if mode == "reload":
                    for rerun in (False, True):
                        if rerun:
                            page.get_by_role(
                                "button", name="Run all", exact=True
                            ).click()
                            expect(
                                page.locator(".notebook-toolbar [role=status]")
                            ).to_have_attribute("data-state", "ready", timeout=120000)
                        drawing = page.frame_locator(
                            'iframe[title="Graph diagram"]'
                        ).last.locator("svg[data-linnet-interactive]")
                        expect(drawing).to_be_visible(timeout=60000)
                        assert drawing.evaluate("""svg => {
                            const doc = svg.ownerDocument, win = doc.defaultView;
                            const button = doc.createElement('button');
                            let trusted = 0;
                            button.addEventListener('click', () => trusted++);
                            button.setAttribute('onclick', 'window.untrustedGraphCallback = true');
                            const script = doc.createElement('script');
                            script.textContent = 'window.untrustedGraphCallback = true';
                            doc.body.append(button, script); button.click();
                            button.remove(); script.remove();
                            return trusted === 1 && !win.untrustedGraphCallback;
                        }""")
                        for theme in ("light", "dark"):
                            page.evaluate(
                                "theme=>document.documentElement.dataset.theme=theme",
                                theme,
                            )
                            expect(drawing).to_have_attribute("data-theme", theme)
                            for kind in ("edge", "node"):
                                target = drawing.locator(
                                    f'[data-linnet-kind="{kind}"]:visible'
                                ).first
                                target.hover()
                                popup = page.locator(".linnet-inspector-popup")
                                expect(popup).to_be_visible()
                                expect(
                                    popup.locator(".linnet-inspector-header")
                                ).to_contain_text("Edge" if kind == "edge" else "Node")
                                target.click()
                                expect(popup).to_have_attribute("data-pinned", "")
                                page.mouse.move(0, 0)
                                expect(popup).to_be_visible()
                                popup.get_by_role(
                                    "button", name="Close graph details"
                                ).click()
                                target.focus()
                                expect(popup).to_be_visible()
                                target.press("Enter")
                                expect(popup).to_have_attribute("data-pinned", "")
                                popup.get_by_role(
                                    "button", name="Close graph details"
                                ).click()
                    assert not errors, errors
                status = page.locator(".notebook-toolbar [role=status]").inner_text()
                if failing:
                    assert {
                        "module-load": "failed to load",
                        "worker-rejection": "injected startup rejection",
                        "dependency-rejection": "injected dependency failure",
                    }[mode] in status
                results.append(
                    {
                        "mode": mode,
                        "status": status,
                        "states": page.evaluate("window.startupStates"),
                        "controller": page.evaluate(
                            "navigator.serviceWorker.controller?.scriptURL"
                        ),
                        "errors": errors,
                    }
                )
                print("Passed", mode, flush=True)
                if failing:
                    fault = mode
                    mode = "recovery"
                    page.reload()
                    page.locator("[data-notebook]").scroll_into_view_if_needed()
                    expect(
                        page.locator(".notebook-toolbar [role=status]")
                    ).to_have_attribute("data-state", "ready", timeout=180000)
                    expect(page.locator("[data-notebook-ready]")).to_contain_text(
                        "Generated 4,970 graphs"
                    )
                    results.append(
                        {
                            "mode": fault + "-recovery",
                            "result": page.locator(
                                "[data-notebook-ready]"
                            ).inner_text(),
                        }
                    )
                    print("Passed", fault + "-recovery", flush=True)
                page.close()
                if failing:
                    current.close()
            context.close()
        finally:
            browser.close()
finally:
    server.shutdown()
    server.server_close()
    thread.join()
    (args.output / "startup.json").write_text(json.dumps(results, indent=2))
