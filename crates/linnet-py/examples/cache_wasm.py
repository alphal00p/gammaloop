"""Exercise cached reloads, cross-product navigation and cache invalidation."""

import hashlib
import json
import shutil
import time
import urllib.parse
import uuid
from pathlib import Path

from profile_wasm import INSTRUMENT


def cache_navigation_smoke(
    browser, base_url, site, report, timeout, tracking, wait_for_cells
):
    context = browser.new_context(service_workers="allow")
    context.add_init_script(tracking + ";\n" + INSTRUMENT)
    page = context.new_page()
    errors = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    observations = []
    fixtures = []
    # A pending CacheStorage write must not hold the response needed by
    # compileStreaming. Use a controlled storage gate rather than timing a
    # particular machine's disk or network.
    page.evaluate(
        """async source => {
        let handler, releaseWrite, wrote;
        const gate = new Promise(resolve => releaseWrite = resolve);
        const writing = new Promise(resolve => wrote = resolve);
        const bytes = new Uint8Array([0,97,115,109,1,0,0,0]);
        const compressed = await new Response(new Blob([bytes]).stream()
            .pipeThrough(new CompressionStream('gzip'))).arrayBuffer();
        let stored;
        const cache = {
            match: async () => stored,
            put: async (_request, response) => {
                await response.arrayBuffer(); wrote(); await gate;
            },
        };
        const scope = {
            location: new URL('https://fixture.test/notebook-cache-worker.js'),
            addEventListener: (name, callback) => {if (name === 'fetch') handler = callback;},
        };
        new Function('self', 'caches', 'fetch', source)(scope,
            {open: async () => cache}, async () => new Response(compressed));
        const request = new Request('https://fixture.test/assets/notebook-runtime/' + 'a'.repeat(64) + '/symbolica-core.wasm');
        const pending = [];
        const dispatch = () => new Promise(resolve => handler({request,
            respondWith: promise => resolve(promise), waitUntil: promise => pending.push(promise)}));
        const response = await Promise.race([dispatch(), new Promise((_, reject) =>
            setTimeout(() => reject(new Error('Cache write blocked streaming response')), 2000))]);
        if (response.headers.get('content-type') !== 'application/wasm') throw new Error('Wrong WASM MIME');
        await WebAssembly.compileStreaming(Promise.resolve(response));
        await writing; releaseWrite(); await Promise.all(pending);
        stored = new Response(bytes, {headers: {'content-type': 'application/wasm'}});
        if (await dispatch() !== stored) throw new Error('Cached response metadata was rewrapped');
    }""",
        Path(__file__).with_name("notebook-cache-worker.js").read_text(),
    )

    def navigate(url, kind, *, reload=False):
        start = time.perf_counter()
        if reload:
            page.reload(wait_until="domcontentloaded")
        else:
            page.goto(url, wait_until="domcontentloaded")
        page.locator("[data-notebook]").scroll_into_view_if_needed()
        wait_for_cells(page, kind)
        page.wait_for_function(
            "document.querySelector('[data-notebook]>[role=status]')?.dataset.state === 'ready'",
            timeout=timeout,
        )
        profile = page.evaluate("window.notebookProfile")
        measures = [
            entry
            for entry in profile["workerResources"]
            if entry.get("entryType") == "measure"
        ]
        observations.append(
            {
                "kind": kind,
                "url": page.url,
                "seconds": time.perf_counter() - start,
                "measures": measures,
                "worker_resources": profile["workerResources"],
            }
        )
        assert not page.locator(".text-error").count()

    def execute_probe(code):
        cell = page.locator(
            'marimo-island:has(marimo-code-editor):has(button[aria-label="Run cell"])'
        ).first
        editor = cell.locator("marimo-code-editor")
        original = json.loads(editor.get_attribute("data-initial-value"))
        previous = page.evaluate("window.notebookCompletedRuns")
        editor.locator(".cm-content").fill(original + "\n" + code)
        cell.locator("..").get_by_role("button", name="Run cell", exact=True).click()
        page.wait_for_function(
            "previous => window.notebookCompletedRuns > previous",
            arg=previous,
            timeout=timeout,
        )
        wait_for_cells(page, "Python state isolation probe")

    try:
        first = urllib.parse.urljoin(base_url, "guides/showcases/tensor-reduction/")
        second = urllib.parse.urljoin(
            base_url, "../../spenso/latest/guides/tensor-explorer/"
        )
        navigate(first, "fresh-context")
        assert page.evaluate("Boolean(navigator.serviceWorker.controller)")
        execute_probe(
            "import builtins as _probe\n_probe._alphal00p_navigation_probe = True"
        )
        navigate(first, "reload-new-worker", reload=True)
        execute_probe(
            "import builtins as _probe\nassert not hasattr(_probe, '_alphal00p_navigation_probe')\n_probe._alphal00p_navigation_probe = True"
        )
        navigate(second, "switch-product-new-worker")
        execute_probe(
            "import builtins as _probe\nassert not hasattr(_probe, '_alphal00p_navigation_probe')"
        )
        assert not errors, errors
        # Both products must request one shared immutable support wheel.
        first_payload = json.loads(
            (
                site
                / "products/feynkit/latest/assets/notebooks/07_tensor_reduction_marimo.json"
            ).read_text()
        )
        second_payload = json.loads(
            (
                site
                / "products/spenso/latest/assets/notebooks/spenso_tensor_explorer.json"
            ).read_text()
        )
        # Resolve against the actual product asset directories, not the guide URLs.
        origin = urllib.parse.urljoin(base_url, "../../..").rstrip("/")
        first_url = urllib.parse.urljoin(
            origin + "/products/feynkit/latest/assets/notebooks/x.json",
            first_payload["runtime_wheel"],
        )
        second_url = urllib.parse.urljoin(
            origin + "/products/spenso/latest/assets/notebooks/x.json",
            second_payload["runtime_wheel"],
        )
        assert first_url == second_url, (first_url, second_url)
        cache_urls = page.evaluate("""async () => {
            const cache = await caches.open('alphal00p-notebook-runtime-v1');
            return (await cache.keys()).map(request => request.url);
        }""")
        assert first_url in cache_urls
        assert all("/assets/notebook-runtime/" in url for url in cache_urls)
        # CacheStorage remains effective independently of server cache headers.
        context.set_offline(True)
        offline = page.evaluate(
            """async url => {
            const response = await fetch(url);
            return {status: response.status, bytes: (await response.arrayBuffer()).byteLength,
                    cacheControl: response.headers.get('cache-control')};
        }""",
            first_url,
        )
        context.set_offline(False)
        assert offline["status"] == 200 and offline["bytes"] > 0
        assert "immutable" in offline["cacheControl"]
        # Verify immutable versions cannot leak bodies into one another, and a
        # transient 404 does not poison a later successful request.
        urls = []
        for contents in ("first version", "second version", None):
            identity = hashlib.sha256(uuid.uuid4().bytes).hexdigest()
            directory = site / "assets/notebook-runtime" / identity
            directory.mkdir()
            fixtures.append(directory)
            if contents is not None:
                (directory / "cache-probe.whl").write_text(contents)
            urls.append(
                urllib.parse.urljoin(
                    origin + "/",
                    "assets/notebook-runtime/" + identity + "/cache-probe.whl",
                )
            )
        fetch = """async url => {const response=await fetch(url,{cache:'reload'});return {status:response.status,text:await response.text()};}"""
        assert page.evaluate(fetch, urls[0])["text"] == "first version"
        page.wait_for_function(
            "async url => Boolean(await (await caches.open('alphal00p-notebook-runtime-v1')).match(url))",
            arg=urls[0],
            timeout=timeout,
        )
        (fixtures[0] / "cache-probe.whl").write_text("changed at an immutable URL")
        assert page.evaluate(fetch, urls[0])["text"] == "first version"
        assert page.evaluate(fetch, urls[1])["text"] == "second version"
        assert page.evaluate(fetch, urls[2])["status"] == 404
        (fixtures[2] / "cache-probe.whl").write_text("recovered")
        assert page.evaluate(fetch, urls[2]) == {"status": 200, "text": "recovered"}
        report = Path(report)
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(
            json.dumps(
                {
                    "observations": observations,
                    "shared_wheel_url": first_url,
                    "cached_runtime_urls": cache_urls,
                    "fresh_python_state_on_reload_and_product_switch": True,
                    "immutable_version_isolation": True,
                    "failed_responses_not_cached": True,
                    "errors": errors,
                },
                indent=2,
            )
        )
        print(
            "Notebook reload, cross-product cache, state isolation and error/version invalidation passed",
            flush=True,
        )
    finally:
        context.close()
        for directory in fixtures:
            shutil.rmtree(directory)
