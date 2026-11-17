"""Opt-in browser regressions; run manually with an installed Playwright browser.

Use --ready with the existing mathml_view_browser_producer. Its eight native
views supply requests, bundled STIX CSS and immutable captures. Native jobs
finish first; cache stress then runs independently of the finite producer.
--fixture reuses the retained fixture for focused reruns without native work.
This client never builds, imports Symbolica, or invents native SVG geometry.
"""

import argparse
import hashlib
import importlib.util
import json
import signal
import sys
import time
import traceback
import urllib.parse
import urllib.request
from pathlib import Path

from mathml_view_browser import launch, post
from playwright.sync_api import sync_playwright

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CASES = (
    "reuse",
    "concurrent",
    "isolation",
    "font-protocol",
    "failed-font",
    "native-error-lease",
    "metric-limit",
    "metric-bytes",
    "font-lifetime",
)


def wrapper_source(javascript):
    """Use the actual Python adapter's wrapper, rather than a parallel loader."""
    source = ROOT / "crates/linnet-render-py/src/view.py"
    spec = importlib.util.spec_from_file_location("graph_reuse_adapter", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return str(module._widget_class(javascript)._esm)


def fetch_json(url):
    with urllib.request.urlopen(url, timeout=10) as response:
        return json.load(response)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--ready", type=Path)
    inputs.add_argument("--fixture", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--case", action="append", choices=CASES, dest="cases")
    parser.add_argument("--engine", action="append", choices=("chromium", "webkit"))
    args = parser.parse_args()
    assert not args.ready or not args.engine, (
        "--ready consumes all eight jobs in both engines"
    )
    output = args.output.resolve()
    allowed = (
        ROOT / ".local/graph-mathml-labels/browser",
        ROOT / ".local/graph-mathml-reuse",
    )
    assert any(output.is_relative_to(path.resolve()) for path in allowed), output
    output.mkdir(parents=True, exist_ok=True)
    javascript = (ROOT / "crates/linnet-render-py/src/view.js").read_text(
        encoding="utf-8"
    )
    wrapper = wrapper_source(javascript)
    source_hash = hashlib.sha256(javascript.encode()).hexdigest()
    report = {"schema": 1, "source_sha256": source_hash, "matrix": [], "ok": False}
    producer = None
    producer_finished = False
    started = time.monotonic()

    def deadline(_signum, _frame):
        raise TimeoutError("graph reuse browser matrix exceeded 360 seconds")

    old_alarm = signal.signal(signal.SIGALRM, deadline)
    signal.alarm(360)
    try:
        if args.ready:
            producer = json.loads(args.ready.read_text(encoding="utf-8"))
            origin = producer["url"]
            parsed = urllib.parse.urlparse(origin)
            assert (
                parsed.scheme == "http"
                and parsed.hostname == "127.0.0.1"
                and parsed.port
            )
            manifest = fetch_json(origin + "/manifest")
            assert manifest["native"] and len(manifest["requests"]) == 8
            assert manifest["source_sha256"] == source_hash, (
                "producer compiled stale production JavaScript"
            )
            fixtures = {}
        else:
            retained = json.loads(args.fixture.read_text(encoding="utf-8"))
            assert retained["source_sha256"] == source_hash, (
                "retained fixture belongs to another source freeze"
            )
            fixtures = retained["fixtures"]
            origin = "http://127.0.0.1:19731"

        with sync_playwright() as driver:
            browsers = {}
            try:
                for engine in args.engine or ("chromium", "webkit"):
                    browsers[engine] = launch(driver, engine)

                def run_case(engine, name, request, width=350):
                    browser = browsers[engine]
                    context = browser.new_context(
                        viewport={"width": width, "height": 1000}
                    )
                    denied, errors = [], []

                    def route_local(route):
                        address = route.request.url
                        path = urllib.parse.urlparse(address).path
                        if address == origin + "/":
                            route.fulfill(
                                content_type="text/html",
                                body="<!doctype html><body></body>",
                            )
                        elif address == origin + "/graph-reuse.js":
                            route.fulfill(
                                content_type="text/javascript",
                                body=(HERE / "graph_reuse_browser.js").read_text(
                                    encoding="utf-8"
                                ),
                            )
                        elif address == origin + "/harness.js":
                            route.fulfill(
                                content_type="text/javascript",
                                body=(HERE / "mathml_view_browser.js").read_text(
                                    encoding="utf-8"
                                ),
                            )
                        elif address == origin + "/view.js":
                            route.fulfill(
                                content_type="text/javascript", body=javascript
                            )
                        elif (
                            producer
                            and not producer_finished
                            and address.startswith(origin + "/")
                            and (
                                path == "/manifest"
                                or path.startswith(("/receipt/", "/saved/"))
                            )
                        ) or address.startswith(("blob:", "data:")):
                            route.continue_()
                        else:
                            denied.append(address)
                            route.abort()

                    context.route("**/*", route_local)
                    page = context.new_page()
                    page.on("pageerror", lambda error: errors.append(str(error)))
                    try:
                        page.goto(origin + "/", wait_until="load")
                        result = page.evaluate(
                            """async input => {
                                const suite = await import('/graph-reuse.js');
                                return suite.run(input);
                            }""",
                            {
                                "name": name,
                                "request": request,
                                "wrapper": wrapper,
                                "javascript": javascript,
                                "engine": engine,
                                "width": width,
                            },
                        )
                        assert not denied, f"unexpected network: {denied}"
                        assert not errors, f"uncaught browser failures: {errors}"
                        report["matrix"].append(
                            {
                                "engine": engine,
                                "version": browser.version,
                                "width": width,
                                **result,
                            }
                        )
                        print(f"PASS {engine} {width} {name}", flush=True)
                    except Exception:
                        page.screenshot(
                            path=str(output / f"{engine}-{width}-{name}-failure.png")
                        )
                        report["failure_context"] = {
                            "engine": engine,
                            "width": width,
                            "case": name,
                            "denied": denied,
                            "errors": errors,
                        }
                        raise
                    finally:
                        context.close()

                if producer:
                    for engine in browsers:
                        for width in (350, 1400):
                            request = manifest["requests"][f"{engine}-{width}-a"]
                            run_case(engine, "native", request, width)
                        fixtures[engine] = {
                            **manifest["requests"][f"{engine}-350-a"],
                            "completed_html": (
                                args.ready.parent / f"{engine}-350-a.html"
                            ).read_text(encoding="utf-8"),
                        }
                    post(
                        origin + "/finish",
                        {"ok": True, "native_capture_only": True, "matrix_jobs": 4},
                    )
                    producer_finished = True
                    (output / "fixture.json").write_text(
                        json.dumps(
                            {
                                "source_sha256": source_hash,
                                "fixtures": fixtures,
                            }
                        ),
                        encoding="utf-8",
                    )
                for engine in browsers:
                    for name in args.cases or CASES:
                        run_case(engine, name, fixtures[engine])
            finally:
                for browser in browsers.values():
                    browser.close()
        assert "symbolica" not in sys.modules, "browser-only client imported Symbolica"
        report["ok"] = True
    except Exception as error:  # noqa: BLE001 -- retain finite producer/browser failure evidence.
        report.update(error=str(error), traceback=traceback.format_exc())
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_alarm)
        if producer and not producer_finished:
            try:
                post(
                    producer["url"] + "/finish",
                    {"ok": False, "error": report.get("error")},
                )
            except (OSError, ValueError) as error:
                report["finish_error"] = str(error)
        report["elapsed_seconds"] = time.monotonic() - started
        (output / "graph-reuse-report.json").write_text(
            json.dumps(report, indent=2), encoding="utf-8"
        )
        print(
            json.dumps(
                {
                    "ok": report["ok"],
                    "report": str(output / "graph-reuse-report.json"),
                    "error": report.get("error"),
                }
            )
        )
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
