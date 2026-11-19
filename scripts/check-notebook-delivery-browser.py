"""Check negotiated delivery, real startup, and verified old-cache migration."""

import argparse
import http.client
import importlib.util
import json
import os
import socket
import subprocess
import tempfile
import time
from pathlib import Path

from playwright.sync_api import expect, sync_playwright

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--site", required=True, type=Path)
parser.add_argument("--output", required=True, type=Path)
parser.add_argument("--caddy", required=True, type=Path)
parser.add_argument("--previous-runtime", required=True)
parser.add_argument("--previous-cache-worker", required=True)
engines = parser.add_mutually_exclusive_group(required=True)
engines.add_argument("--chromium", type=Path)
engines.add_argument("--webkit", type=Path)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location(
    "profile_wasm", Path("crates/linnet-py/examples/profile_wasm.py")
)
profile = importlib.util.module_from_spec(spec)
spec.loader.exec_module(profile)
route = "/products/feynkit/latest/guides/showcases/four-loop-numerator/"
payload = json.loads(
    (
        args.site
        / "products/feynkit/latest/assets/notebooks/05_four_loop_numerator_marimo.json"
    ).read_text()
)
runtime = next(
    path.parent
    for path in (args.site / "assets/notebook-runtime").glob("*/runtime-manifest.json")
    if json.loads(path.read_text()).get("shared_files")
)
manifest = json.loads((runtime / "runtime-manifest.json").read_text())
shared = manifest["shared_files"]
wasm = f"/assets/notebook-runtime-files/{shared['symbolica-core.wasm']}/symbolica-core.wasm"
archives = [
    {
        "old": f"/assets/notebook-runtime/{args.previous_runtime}/{name}",
        "new": f"/assets/notebook-runtime-files/{digest}/{Path(name).name}",
        "sha256": digest,
    }
    for name, digest in shared.items()
    if name.endswith((".whl", ".zip"))
]
with socket.socket() as listener:
    listener.bind(("127.0.0.1", 0))
    port = listener.getsockname()[1]
base = f"http://127.0.0.1:{port}"
results = {"runtime": manifest["identity"], "trials": []}
access = args.output.resolve() / "access.jsonl"
with tempfile.TemporaryDirectory() as temporary:
    config = Path(temporary) / "Caddyfile"
    config.write_text(
        Path("scripts/website-preview.Caddyfile")
        .read_text()
        .replace("\troot *", f"\tlog {{\n\t\toutput file {access}\n\t}}\n\troot *")
    )
    environment = dict(
        os.environ,
        WEBSITE_PREVIEW_PORT=str(port),
        WEBSITE_PREVIEW_ROOT=str(args.site.resolve()),
        XDG_CONFIG_HOME=str(Path(temporary) / "config"),
        XDG_DATA_HOME=str(Path(temporary) / "data"),
    )
    with (args.output / "server.log").open("w") as log:
        server = subprocess.Popen(
            [str(args.caddy), "run", "--config", str(config), "--adapter", "caddyfile"],
            env=environment,
            stdout=log,
            stderr=log,
        )
        try:
            for _ in range(100):
                try:
                    with socket.create_connection(("127.0.0.1", port), timeout=0.1):
                        break
                except OSError:
                    if server.poll() is not None:
                        raise RuntimeError("Preview server exited")
                    time.sleep(0.05)
            connection = http.client.HTTPConnection("127.0.0.1", port)
            checks = []
            for path, headers in [
                (wasm, {"Accept-Encoding": "gzip, br"}),
                (wasm, {"Accept-Encoding": "identity"}),
                (wasm + ".gz", {"Accept-Encoding": "gzip"}),
                ("/assets/site.js", {}),
            ]:
                connection.request("HEAD", path, headers=headers)
                response = connection.getresponse()
                checks.append(
                    {
                        "path": path,
                        "status": response.status,
                        "headers": dict(response.getheaders()),
                    }
                )
                response.read()
            encoded, raw, explicit, mutable = [item["headers"] for item in checks]
            assert (
                encoded["Content-Encoding"] == "gzip"
                and encoded["Content-Type"] == "application/wasm"
            )
            assert int(encoded["Content-Length"]) < int(raw["Content-Length"])
            assert (
                encoded["Vary"] == "Accept-Encoding"
                and "immutable" in encoded["Cache-Control"]
            )
            assert (
                "Content-Encoding" not in explicit
                and explicit["Content-Length"] == encoded["Content-Length"]
            )
            assert mutable["Cache-Control"] == "no-cache"
            connection.request("HEAD", wasm + ".missing")
            response = connection.getresponse()
            assert response.status == 404 and "immutable" not in response.getheader(
                "Cache-Control", ""
            )
            response.read()
            connection.request(
                "GET",
                wasm,
                headers={"Accept-Encoding": "gzip", "If-None-Match": encoded["Etag"]},
            )
            response = connection.getresponse()
            assert response.status == 304 and not response.read()
            connection.request(
                "GET",
                wasm,
                headers={"Accept-Encoding": "identity", "Range": "bytes=0-7"},
            )
            response = connection.getresponse()
            assert response.status == 206 and response.read() == b"\0asm\x01\0\0\0"
            connection.close()
            results["http"] = checks
            results["conditional304"] = results["range206"] = True
            with sync_playwright() as driver:
                browser = (
                    driver.webkit.launch(executable_path=str(args.webkit))
                    if args.webkit
                    else driver.chromium.launch(
                        executable_path=str(args.chromium), args=["--no-sandbox"]
                    )
                )
                try:
                    # Block all service workers: native compression must still
                    # work through the ordinary HTTP representation negotiation.
                    context = browser.new_context(service_workers="block")
                    page = context.new_page()
                    page.add_init_script(profile.INSTRUMENT)
                    page.add_init_script("""window.states=[];new MutationObserver(()=>{
                      const s=document.querySelector('.notebook-toolbar [role=status]');
                      if(s&&window.states.at(-1)!==s.textContent)window.states.push(s.textContent);
                    }).observe(document,{subtree:true,childList:true,characterData:true});""")

                    def delayed(request):
                        # A controlled five-second unresolved fetch distinguishes
                        # response-header waiting from measured body progress.
                        time.sleep(5)
                        request.continue_()

                    page.route("**/python_stdlib.zip", delayed)
                    page.goto(base + route)
                    page.locator("[data-notebook]").scroll_into_view_if_needed()
                    expect(
                        page.locator(".notebook-toolbar [role=status]")
                    ).to_have_attribute("data-state", "ready", timeout=180000)
                    expect(page.locator("[data-notebook-ready]")).to_contain_text(
                        "Generated 4,970 graphs"
                    )
                    states = page.evaluate("window.states")
                    assert any("awaiting response headers" in state for state in states)
                    assert not any("receiving bytes" in state for state in states)
                    resources = page.evaluate("window.notebookProfile.workerResources")
                    native = [
                        entry
                        for entry in resources
                        if entry.get("name", "").endswith("/symbolica-core.wasm")
                    ]
                    assert native, "Native resource completion must be observed"
                    assert any(
                        entry["encodedBodySize"] == int(encoded["Content-Length"])
                        and entry["decodedBodySize"] == int(raw["Content-Length"])
                        for entry in native
                    ), native
                    results["trials"].append(
                        {
                            "mode": "service-workers-blocked",
                            "native": native,
                            "states": states,
                            "scientific_default": page.locator(
                                "[data-notebook-ready]"
                            ).inner_text(),
                        }
                    )
                    context.close()
                    context = browser.new_context()
                    page = context.new_page()
                    page.goto(base + "/")
                    page.evaluate(
                        """async worker => {
                      await navigator.serviceWorker.register(worker);await navigator.serviceWorker.ready;
                      if(!navigator.serviceWorker.controller)await new Promise(resolve=>navigator.serviceWorker.addEventListener('controllerchange',resolve,{once:true}));
                    }""",
                        "/" + args.previous_cache_worker,
                    )
                    page.evaluate(
                        """async archives=>{
                      const cache=await caches.open('alphal00p-notebook-runtime-v1');
                      for(const archive of archives) {
                        const response=await fetch(archive.old);
                        if(!response.ok)throw Error(archive.old);
                        await cache.put(archive.old,response);
                      }
                    }""",
                        archives,
                    )
                    old_controller = page.evaluate(
                        "navigator.serviceWorker.controller.scriptURL"
                    )
                    page.evaluate(
                        """async worker=>{
                      const changed=new Promise(resolve=>navigator.serviceWorker.addEventListener('controllerchange',resolve,{once:true}));
                      await navigator.serviceWorker.register(worker);await changed;
                    }""",
                        "/" + Path(payload["cache_worker"]).name,
                    )
                    migrated = page.evaluate(
                        """async archives=>{
                      const result=[];
                      for(const archive of archives) {
                        const response=await fetch(archive.new);
                        const hash=[...new Uint8Array(await crypto.subtle.digest('SHA-256',await response.arrayBuffer()))].map(n=>n.toString(16).padStart(2,'0')).join('');
                        if(hash!==archive.sha256)throw Error('Migration hash mismatch');
                        result.push(archive.new);
                      }
                      return result;
                    }""",
                        archives,
                    )
                    new_paths = set(migrated)
                    logged = [
                        json.loads(line) for line in access.read_text().splitlines()
                    ]
                    # Earlier blocked startup downloaded packages. Only examine
                    # the log after the final old URL used to seed the cache.
                    last_old = max(
                        i
                        for i, entry in enumerate(logged)
                        if entry.get("request", {}).get("uri") == archives[-1]["old"]
                    )
                    assert not [
                        entry
                        for entry in logged[last_old + 1 :]
                        if entry.get("request", {}).get("uri") in new_paths
                    ]
                    results["migration"] = {
                        "old_controller": old_controller,
                        "new_controller": page.evaluate(
                            "navigator.serviceWorker.controller.scriptURL"
                        ),
                        "verified_archives": len(migrated),
                        "network_requests_for_migrated_archives": 0,
                        "urls": migrated,
                    }
                    for mode in ("active-upgraded-cache", "warm"):
                        page.goto(base + route)
                        page.locator("[data-notebook]").scroll_into_view_if_needed()
                        expect(
                            page.locator(".notebook-toolbar [role=status]")
                        ).to_have_attribute("data-state", "ready", timeout=180000)
                        expect(page.locator("[data-notebook-ready]")).to_contain_text(
                            "Generated 4,970 graphs"
                        )
                        results["trials"].append(
                            {
                                "mode": mode,
                                "scientific_default": page.locator(
                                    "[data-notebook-ready]"
                                ).inner_text(),
                            }
                        )
                    context.close()
                finally:
                    browser.close()
        finally:
            server.terminate()
            server.wait(timeout=10)
(args.output / "delivery.json").write_text(json.dumps(results, indent=2))
print(
    "Verified compressed HTTP delivery, conditional/range requests, blocked/active/warm scientific startup and verified old-cache migration",
    flush=True,
)
