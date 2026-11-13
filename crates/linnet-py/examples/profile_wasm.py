"""Record browser startup phases and network requests without changing notebooks.

Run with the Nix notebook exporter's Python and Chromium. Cold means a fresh
browser context; warm means a second navigation in that same context. Both sites
use the same browser process, host network, viewport and instrumentation.
"""

import argparse
import json
import os
import select
import socket
import socketserver
import threading
import time
from pathlib import Path

from playwright.sync_api import Error, sync_playwright

INSTRUMENT = r"""(() => {
  window.notebookProfile = {events: [], fetches: [], workerResources: []};
  document.addEventListener('DOMContentLoaded', () => {
    window.notebookProfile.presentation = {};
    const inspect = setInterval(() => {
      const milestones = window.notebookProfile.presentation;
      const roots = [document];
      for (let i = 0; i < roots.length; i++) {
        for (const element of roots[i].querySelectorAll('*')) {
          if (element.shadowRoot) roots.push(element.shadowRoot);
        }
      }
      for (const [name, selector] of Object.entries({
        firstOutput: '.output:not(:empty)',
        firstWidgetHost: 'iframe',
        firstPlot: '.output img, marimo-matplotlib',
        integrationControl: 'button[type=submit]',
      })) {
        if (!milestones[name] && roots.some(root => [...root.querySelectorAll(selector)].some(e => e.getBoundingClientRect().height)))
          milestones[name] = performance.now();
      }
    }, 100);
    addEventListener('pagehide', () => clearInterval(inspect), {once:true});
  });
  const record = (event) => {
    const rpc = event.data?.data || event.data;
    if (rpc?.profileResource) {
      window.notebookProfile.workerResources.push(rpc.profileResource);
      return;
    }
    if (rpc?.notebookProfile) {
      window.notebookProfile.fetches.push(rpc.notebookProfile);
      return;
    }
    const entry = {time: performance.now(), id: rpc?.id, method: rpc?.method,
                   type: rpc?.type, success: rpc?.success};
    if (rpc?.id === 'initializingMessage') entry.phase = rpc.payload?.message || rpc.payload;
    if (rpc?.id === 'kernelMessage') {
      try {
        const message = JSON.parse(rpc.payload.message);
        entry.op = message.op;
        entry.status = message.data?.status;
        entry.cell = message.data?.cell_id;
        entry.outputType = message.data?.output?.mimetype;
        if (entry.outputType === 'application/vnd.marimo+error') entry.error = message.data.output.data;
      } catch {}
    }
    window.notebookProfile.events.push(entry);
  };
  // Upstream revokes its worker Blob immediately after construction. The wrapper
  // imports that Blob asynchronously, so retain it until import has completed.
  const revoke = URL.revokeObjectURL;
  URL.revokeObjectURL = url => setTimeout(() => revoke(url), 30000);
  const Worker = window.Worker;
  window.Worker = class extends Worker {
    constructor(url, options) {
      const instrumentation = `
        new PerformanceObserver(list => {
          for (const entry of list.getEntries()) {
            const resource = entry.toJSON();
            if (entry.entryType === "measure") resource.detail = entry.detail;
            self.postMessage({profileResource: JSON.parse(JSON.stringify(resource))});
          }
        }).observe({entryTypes: ["resource", "measure"]});
        const originalFetch = self.fetch;
        self.fetch = async (...args) => {
          const start = performance.now();
          const response = await originalFetch(...args);
          const address = args[0]?.url || String(args[0]);
          self.postMessage({notebookProfile: {url: address, start,
            headers: performance.now(), status: response.status,
            contentLength: response.headers.get('content-length'), timeOrigin: performance.timeOrigin}});
          return response;
        };
        const queued = [];
        const hold = event => { event.stopImmediatePropagation(); queued.push(event.data); };
        self.addEventListener('message', hold);
        await import(${JSON.stringify(String(url))});
        self.removeEventListener('message', hold);
        for (const data of queued) self.dispatchEvent(new MessageEvent('message', {data}));`;
      const wrapper = URL.createObjectURL(new Blob([instrumentation], {type:'application/javascript'}));
      super(wrapper, {...options, type:'module'});
      URL.revokeObjectURL(wrapper);
      this.addEventListener('message', record);
      const post = this.postMessage;
      this.postMessage = (...args) => {
        record({data: args[0]});
        return post.apply(this, args);
      };
    }
  };
})();"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("site")
    parser.add_argument("chromium")
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--only-alpha", action="store_true")
    parser.add_argument(
        "--route",
        default="tensor-reduction",
        choices=["tensor-reduction", "numerical-integration"],
    )
    parser.add_argument("--latency-ms", type=float, default=0)
    parser.add_argument("--download-mbps", type=float, default=0)
    parser.add_argument("--read-delay", type=float, default=0)
    parser.add_argument("--viewport-width", type=int, default=1280)
    parser.add_argument("--viewport-height", type=int, default=800)
    parser.add_argument(
        "--stdlib-url",
        help="Override only the prepared worker's stdlib URL for a same-origin source/bytecode comparison",
    )
    parser.add_argument(
        "--upstream-runtime",
        action="store_true",
        help="Compare the current Python exporter's unmodified upstream islands runtime on the minimal workload",
    )
    parser.add_argument(
        "--no-stream",
        action="store_true",
        help="Compare full-wheel loading with the streaming loader on the minimal workload",
    )
    parser.add_argument(
        "--minimal",
        action="store_true",
        help="Use one checked rank-two tensor computation on alphal00p; public workload stays unchanged",
    )
    options = parser.parse_args()
    if options.route != "tensor-reduction" and (
        not options.only_alpha or options.minimal
    ):
        parser.error(
            "The numerical route requires --only-alpha and its complete published workload"
        )
    if min(options.latency_ms, options.download_mbps, options.read_delay) < 0:
        parser.error("Network simulation and reading delay must be nonnegative")
    if options.stdlib_url and (not options.minimal or options.upstream_runtime):
        parser.error("--stdlib-url requires --minimal and the prepared runtime")
    if options.upstream_runtime and not options.minimal:
        parser.error("--upstream-runtime requires --minimal")
    if options.no_stream and not options.minimal:
        parser.error("--no-stream requires --minimal")
    if options.repetitions < 1:
        parser.error("--repetitions must be positive")
    options.report.mkdir(parents=True, exist_ok=True)
    results = []

    class NetworkProxy(socketserver.StreamRequestHandler):
        def handle(self):
            method, address, _ = self.rfile.readline().decode().split()
            while self.rfile.readline().strip():
                pass
            if method != "CONNECT":
                self.wfile.write(b"HTTP/1.1 405 Method Not Allowed\r\n\r\n")
                return
            host, port = address.rsplit(":", 1)
            with socket.create_connection((host, int(port))) as upstream:
                time.sleep(options.latency_ms / 1000)
                self.wfile.write(b"HTTP/1.1 200 Connection Established\r\n\r\n")
                while True:
                    readable, _, _ = select.select(
                        [self.connection, upstream], [], [], 60
                    )
                    if not readable:
                        return
                    for source in readable:
                        data = source.recv(16384)
                        if not data:
                            return
                        if source is upstream and options.download_mbps:
                            time.sleep(len(data) / (options.download_mbps * 125000))
                        if source is upstream:
                            proxy.download_bytes += len(data)
                        (self.connection if source is upstream else upstream).sendall(
                            data
                        )

    proxy = None
    if options.latency_ms or options.download_mbps:
        proxy = socketserver.ThreadingTCPServer(("127.0.0.1", 0), NetworkProxy)
        proxy.daemon_threads = True
        proxy.download_bytes = 0
        threading.Thread(target=proxy.serve_forever, daemon=True).start()
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=options.chromium,
            proxy={"server": f"http://127.0.0.1:{proxy.server_address[1]}"}
            if proxy
            else None,
        )
        sites = [
            (
                "alphal00p",
                options.site.rstrip("/")
                + f"/products/feynkit/latest/guides/showcases/{options.route}/",
            ),
            (
                "symbolica",
                "https://symbolica.io/gallery/notebooks/tensor_reduction.html",
            ),
        ]
        if options.only_alpha:
            sites = sites[:1]
        for name, url, repetition in (
            (name, url, repetition)
            for repetition in range(options.repetitions)
            for name, url in sites
        ):
            context = browser.new_context(
                viewport={
                    "width": options.viewport_width,
                    "height": options.viewport_height,
                },
                service_workers="block" if options.minimal else "allow",
            )
            if options.minimal and name == "alphal00p":
                import marimo as mo
                from export_wasm import browser_install_code

                generator = mo.MarimoIslandGenerator()
                generator.add_code(browser_install_code(["__NOTEBOOK_WHEEL_URL__"]))
                generator.add_code(
                    """from symbolica.community.tensor import Representation, TensorName
rep = Representation.euc(2)
a = TensorName("startup_probe::A")(rep)
assert a.rank == 1
result = a("i") * a("i")
assert result.rank == 0
result
""",
                    display_code=True,
                )
                body = generator.render_body(
                    include_init_island=False, include_payload=True
                )

                head = generator.render_head()

                def minimal_payload(route, _request, body=body, head=head):
                    response = route.fetch()
                    payload = response.json()
                    payload["body"] = body
                    if options.no_stream or options.upstream_runtime:
                        payload["runtime_wheel"] = None
                    if options.upstream_runtime:
                        payload["runtime"] = None
                        payload["head"] = head
                    route.fulfill(response=response, json=payload)

                context.route("**/07_tensor_reduction_marimo.json", minimal_payload)
                if options.no_stream or options.stdlib_url:

                    def runtime_worker(route):
                        response = route.fetch()
                        source = response.text()
                        if options.no_stream:
                            source = source.replace(
                                'import {mountNative} from "../native-loader.js";\n', ""
                            )
                            source = source.replace("await mountNative(r);", "")
                            source = source.replace(
                                "../wheels/symbolica-", "../../public/wheels/symbolica-"
                            )
                        if options.stdlib_url:
                            expected = "stdLibURL:new URL('../python_stdlib.zip',import.meta.url).href"
                            if source.count(expected) != 1:
                                raise ValueError(
                                    "Review the pinned stdlib override boundary"
                                )
                            source = source.replace(
                                expected, "stdLibURL:" + json.dumps(options.stdlib_url)
                            )
                        route.fulfill(response=response, body=source)

                    context.route("**/worker-D7hWi02M.js", runtime_worker)
            for cache in ("cold", "warm"):
                page = context.new_page()
                page.add_init_script(INSTRUMENT)
                cdp = context.new_cdp_session(page)
                cdp.send("Network.enable")
                cache_events = {}
                served_from_cache = set()

                def cached_request(
                    event, entries=cache_events, cached=served_from_cache
                ):
                    cached.add(event["requestId"])
                    if event["requestId"] in entries:
                        entries[event["requestId"]]["servedFromCache"] = True

                def response_received(
                    event, entries=cache_events, cached=served_from_cache
                ):
                    response = event["response"]
                    keep = {
                        "cache-control",
                        "etag",
                        "last-modified",
                        "age",
                        "expires",
                        "content-length",
                        "content-encoding",
                        "content-type",
                        "vary",
                    }
                    entries[event["requestId"]] = {
                        "url": response["url"],
                        "status": response["status"],
                        "fromDiskCache": response.get("fromDiskCache", False),
                        "fromServiceWorker": response.get("fromServiceWorker", False),
                        "fromPrefetchCache": response.get("fromPrefetchCache", False),
                        "servedFromCache": event["requestId"] in cached,
                        "headers": {
                            key.lower(): value
                            for key, value in response["headers"].items()
                            if key.lower() in keep
                        },
                    }

                def loading_finished(event, entries=cache_events):
                    if event["requestId"] in entries:
                        entries[event["requestId"]]["encodedDataLength"] = event[
                            "encodedDataLength"
                        ]

                cdp.on("Network.requestServedFromCache", cached_request)
                cdp.on("Network.responseReceived", response_received)
                cdp.on("Network.loadingFinished", loading_finished)
                network = []
                errors = []
                console_errors = []
                page.on(
                    "console",
                    lambda message, console_errors=console_errors: (
                        console_errors.append(message.text)
                        if message.type == "error"
                        else None
                    ),
                )
                page.on(
                    "pageerror", lambda error, errors=errors: errors.append(str(error))
                )

                def request_finished(request, entries=network):
                    item = {"url": request.url, "timing": request.timing}
                    response = request.response()
                    item["fromServiceWorker"] = bool(
                        response and response.from_service_worker
                    )
                    item["serviceWorkerRequest"] = request.service_worker is not None
                    try:
                        item["sizes"] = request.sizes()
                    except Error as error:
                        item["sizes"] = {"unavailable": str(error)}
                    entries.append(item)

                context.on("requestfinished", request_finished)
                started = time.perf_counter()
                host_cpu_start = list(
                    map(
                        int, Path("/proc/stat").read_text().splitlines()[0].split()[1:9]
                    )
                )
                proxy_start = proxy.download_bytes if proxy else None
                result = {
                    "site": name,
                    "url": url,
                    "cache": cache,
                    "repetition": repetition,
                    "workload": "minimal-rank-two"
                    if options.minimal and name == "alphal00p"
                    else "published",
                }
                try:
                    page.goto(url, wait_until="domcontentloaded", timeout=120000)
                    if name == "alphal00p":
                        result["notebookTopBeforeScroll"] = page.locator(
                            "[data-notebook]"
                        ).evaluate("e => e.getBoundingClientRect().top")
                        page.wait_for_timeout(options.read_delay * 1000)
                        result["scrollTime"] = page.evaluate("performance.now()")
                        page.locator("[data-notebook]").scroll_into_view_if_needed()
                    else:
                        # Kernel initialization emits completed-run even when code
                        # remains unexecuted. Wait for the kernel, then request a run.
                        page.wait_for_function(
                            """() => window.notebookProfile.events.some(
                            e => e.op === 'kernel-ready')""",
                            timeout=120000,
                        )
                        # The published full editor starts with unexecuted cells.
                        page.locator('[data-testid="run-button"]').last.click(
                            timeout=120000
                        )
                    page.wait_for_function(
                        """() => {
                        const events = window.notebookProfile.events;
                        const running = events.filter(e => e.status === 'running').at(-1)?.time;
                        const completed = events.filter(e => e.op === 'completed-run').at(-1)?.time;
                        return running !== undefined && completed > running;
                    }""",
                        timeout=300000,
                    )
                    if name == "alphal00p":
                        page.wait_for_function(
                            """() => document.querySelector(
                            '[data-notebook] > [role=status]')?.dataset.state === 'ready'""",
                            timeout=300000,
                        )
                        result["initialCellsSeconds"] = time.perf_counter() - started
                        if options.route == "numerical-integration":
                            # Kernel completion precedes async widget mounting.
                            # This workload has six diagram/expression views
                            # (SVG or MathML) and two initial Matplotlib plots.
                            for index in range(6):
                                page.frame_locator("marimo-anywidget iframe").nth(
                                    index
                                ).locator("svg, math").first.wait_for(
                                    state="visible", timeout=120000
                                )
                                if index == 0:
                                    result["renderedDiagramSeconds"] = (
                                        time.perf_counter() - started
                                    )
                            page.locator("marimo-island img").nth(1).wait_for(
                                state="visible", timeout=120000
                            )
                            result["renderedOutputsSeconds"] = (
                                time.perf_counter() - started
                            )
                    result["completed"] = True
                except Error as error:
                    result["completed"] = False
                    result["failure"] = str(error)
                result["wallSeconds"] = time.perf_counter() - started
                result["proxyDownloadBytesIncludingTLS"] = (
                    proxy.download_bytes - proxy_start if proxy else None
                )
                result["hostLoadAverage"] = os.getloadavg()
                host_cpu_end = list(
                    map(
                        int, Path("/proc/stat").read_text().splitlines()[0].split()[1:9]
                    )
                )
                cpu_delta = [
                    end - start for start, end in zip(host_cpu_start, host_cpu_end)
                ]
                result["hostCpu"] = {
                    "logicalCPUs": os.cpu_count(),
                    "busyFraction": 1 - (cpu_delta[3] + cpu_delta[4]) / sum(cpu_delta),
                }
                print(
                    f"{name} {cache} startup checkpoint: {result['wallSeconds']:.3f}s",
                    flush=True,
                )
                result["raw"] = page.evaluate("window.notebookProfile")
                result["cells"] = (
                    page.evaluate("""() => [...document.querySelectorAll('marimo-island')].map(cell => ({
                    id: cell.cellId, phase: cell.dataset.notebookPhase,
                    code: JSON.parse(cell.querySelector('marimo-code-editor')?.dataset.initialValue || '""').slice(0,120)
                }))""")
                )
                result["resources"] = page.evaluate(
                    "performance.getEntriesByType('resource').map(e => e.toJSON())"
                )
                result["network"] = list(network)
                result["httpCache"] = list(cache_events.values())
                result["cacheConfiguration"] = {
                    "contextRouteInterception": bool(
                        options.minimal and name == "alphal00p"
                    ),
                    "cdpCacheOverride": None,
                    "workerAndPythonRecreated": True,
                    "liveKernelRerunMeasuredSeparately": name == "alphal00p",
                }
                result["errors"] = list(errors)
                result["consoleErrors"] = list(console_errors)
                result["globalStatus"] = page.locator(
                    "[data-notebook] > [role=status]"
                ).all_text_contents()
                result["pythonErrors"] = page.locator(".text-error").all_text_contents()
                if result["pythonErrors"] or errors:
                    result["completed"] = False
                if name == "alphal00p" and result["completed"]:
                    # A fresh notebook session retains installed Python packages.
                    # Measure this separately from a navigation's fresh interpreter.
                    marker = page.evaluate("performance.now()")
                    network_start = len(network)
                    rerun_start = time.perf_counter()
                    page.get_by_role("button", name="Run all", exact=True).click()
                    page.wait_for_function(
                        """marker => window.notebookProfile.events.some(
                            e => e.op === 'completed-run' && e.time > marker)""",
                        arg=marker,
                        timeout=300000,
                    )
                    page.wait_for_function(
                        """() => document.querySelector(
                            '[data-notebook] > [role=status]')?.dataset.state === 'ready'""",
                        timeout=300000,
                    )
                    result["rerun"] = {
                        "errors": errors[len(result["errors"]) :],
                        "pythonErrors": page.locator(".text-error").all_text_contents(),
                        "wallSeconds": time.perf_counter() - rerun_start,
                        "raw": page.evaluate("window.notebookProfile"),
                        "network": network[network_start:],
                        "marker": marker,
                    }
                    if result["rerun"]["errors"] or result["rerun"]["pythonErrors"]:
                        result["completed"] = False
                    if options.route == "numerical-integration":
                        integration_start = time.perf_counter()
                        marker = page.evaluate("performance.now()")
                        page.get_by_role(
                            "button", name="Run integration comparison", exact=True
                        ).click()
                        page.wait_for_function(
                            """marker => window.notebookProfile.events.some(
                            e => e.op === 'completed-run' && e.time > marker)""",
                            arg=marker,
                            timeout=300000,
                        )
                        page.wait_for_function(
                            """() => document.querySelector(
                            '[data-notebook] > [role=status]')?.dataset.state === 'ready'""",
                            timeout=300000,
                        )
                        result["integration"] = {
                            "wallSeconds": time.perf_counter() - integration_start,
                            "pythonErrors": page.locator(
                                ".text-error"
                            ).all_text_contents(),
                        }
                        page.locator("marimo-island img").nth(2).wait_for(
                            state="visible", timeout=120000
                        )
                        result["integration"]["renderedWallSeconds"] = (
                            time.perf_counter() - integration_start
                        )
                        result["integration"]["presentation"] = page.evaluate(
                            "window.notebookProfile.presentation"
                        )
                        if result["integration"]["pythonErrors"]:
                            result["completed"] = False
                results.append(result)
                (options.report / "startup.json").write_text(
                    json.dumps(
                        {
                            "browser": browser.version,
                            "minimal_interception_disables_http_cache": options.minimal,
                            "streaming_disabled_for_comparison": options.no_stream,
                            "upstream_runtime": options.upstream_runtime,
                            "stdlib_override": options.stdlib_url,
                            "exporter_marimo_version": __import__("marimo").__version__,
                            "method": "Fresh context then same-context second navigation; server-local Chromium, not the user browser. Optional CONNECT proxy limits downloads per TLS connection (including workers), with a connection delay; not a complete last-mile network simulation.",
                            "simulation": {
                                "latencyMs": options.latency_ms,
                                "downloadMbps": options.download_mbps,
                                "readDelaySeconds": options.read_delay,
                            },
                            "runs": results,
                        },
                        indent=2,
                    )
                )
                print(
                    f"{name} {cache}: {result['wallSeconds']:.3f}s, completed={result['completed']}",
                    flush=True,
                )
                page.close()
                context.remove_listener("requestfinished", request_finished)
            context.close()
        browser.close()
    if proxy:
        proxy.shutdown()
        proxy.server_close()
    if not all(result["completed"] for result in results):
        raise RuntimeError("Startup profiling encountered errors; inspect startup.json")


if __name__ == "__main__":
    main()
