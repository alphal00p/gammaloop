"""Explicit browser client for the ignored native integration-test producer.

This program never builds, imports Symbolica, creates fixtures or invents SVGs.
Launch the ignored producer separately, then pass its ready.json path.

Playwright discovers its installed browsers by default. Override executables with
MATHML_CHROMIUM_EXECUTABLE or MATHML_WEBKIT_EXECUTABLE. Optional extra library
paths in MATHML_NIX_LD_LIBRARY_PATH are prepended to NIX_LD_LIBRARY_PATH.
Artifacts must stay under the repository's ignored .local/graph-mathml-labels/browser.
"""

import argparse
import base64
import datetime
import hashlib
import json
import os
import signal
import sys
import time
import traceback
import urllib.parse
import urllib.request
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[3]
ARTIFACTS = ROOT / ".local/graph-mathml-labels/browser"


def runtime_environment():
    environment = dict(os.environ)
    paths = environment.get("MATHML_NIX_LD_LIBRARY_PATH", "").split(":")
    paths += environment.get("NIX_LD_LIBRARY_PATH", "").split(":")
    environment["NIX_LD_LIBRARY_PATH"] = ":".join(dict.fromkeys(filter(None, paths)))
    environment["TMPDIR"] = "/tmp"
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return environment


def launch(playwright, engine):
    options = {
        "headless": True,
        "timeout": 20000,
        "env": runtime_environment(),
    }
    if executable := os.environ.get(f"MATHML_{engine.upper()}_EXECUTABLE"):
        options["executable_path"] = executable
    if engine == "chromium":
        options["args"] = ["--no-sandbox", "--disable-dev-shm-usage"]
    return getattr(playwright, engine).launch(**options)


def post(url, content):
    request = urllib.request.Request(
        url,
        data=json.dumps(content).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=10) as response:
        return json.load(response)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ready", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    assert output.is_relative_to(ARTIFACTS.resolve()), output
    output.mkdir(parents=True, exist_ok=True)
    assert args.preflight or args.ready, "--ready is required outside blank preflight"
    report = {
        "schema": 1,
        "mode": "blank-preflight" if args.preflight else "retained-native-roundtrip",
        "started_at": datetime.datetime.now(datetime.UTC).isoformat(),
        "symbolica_imported": False,
        "client_python_version_info": list(sys.version_info),
        "engines": [],
        "matrix": [],
        "source_sha256": {
            name: hashlib.sha256(
                (Path(__file__).parent / name).read_bytes()
            ).hexdigest()
            for name in (
                "mathml_view.rs",
                "mathml_view_support.py",
                "mathml_view_browser.py",
                "mathml_view_browser.js",
            )
        },
    }
    report["production_widget_sha256"] = hashlib.sha256(
        (ROOT / "crates/linnet-render-py/src/view.js").read_bytes()
    ).hexdigest()
    producer = None
    deadline = time.monotonic() + 210

    def time_limit(_signum, _frame):
        raise TimeoutError("browser client exceeded its 225 second limit")

    previous_alarm = signal.signal(signal.SIGALRM, time_limit)
    signal.alarm(225)
    try:
        if not args.preflight:
            producer = json.loads(args.ready.read_text())
            url = producer["url"]
            parsed = urllib.parse.urlparse(url)
            assert (
                parsed.scheme == "http"
                and parsed.hostname == "127.0.0.1"
                and parsed.port
            )
            assert len(producer["jobs"]) == 8
            report["producer"] = producer
        with sync_playwright() as playwright:
            for engine in ("chromium", "webkit"):
                browser = launch(playwright, engine)
                try:
                    expected = "151." if engine == "chromium" else "27.2"
                    assert browser.version.startswith(expected), browser.version
                    report["engines"].append(
                        {"engine": engine, "version": browser.version}
                    )
                    if args.preflight:
                        page = browser.new_page()
                        assert page.url == "about:blank"
                        page.close()
                        continue
                    for width in (350, 1400):
                        assert time.monotonic() < deadline, (
                            "browser matrix exceeded deadline"
                        )
                        context = browser.new_context(
                            viewport={"width": width, "height": 1000}
                        )
                        network = []
                        errors = []

                        def local_only(route, *, denied=network, origin=url):
                            request_url = route.request.url
                            if request_url.startswith(
                                (origin + "/", "blob:" + origin + "/")
                            ):
                                route.continue_()
                            else:
                                denied.append(request_url)
                                route.abort()

                        context.route("**/*", local_only)
                        page = context.new_page()
                        page.set_default_timeout(15000)
                        page.on(
                            "pageerror",
                            lambda error, captured=errors: captured.append(str(error)),
                        )
                        try:
                            page.goto(url + "/", wait_until="load")
                            result = page.evaluate(
                                """async ({engine,width}) => {
                                    window.harness=await import('/harness.js');
                                    return await harness.mountPair(engine,width);
                                }""",
                                {"engine": engine, "width": width},
                            )
                            result.update({"engine": engine, "width": width})
                            report["matrix"].append(result)
                            assert (
                                result["source_sha256"]
                                == report["production_widget_sha256"]
                            ), "native producer embeds a stale widget module"
                            assert result["completion_order"] == [
                                f"{engine}-{width}-b",
                                f"{engine}-{width}-a",
                            ]
                            page.screenshot(
                                path=str(output / f"{engine}-{width}.png"),
                                full_page=True,
                            )
                            result["lifecycle"] = page.evaluate(
                                "harness.disposeAndRemount()"
                            )
                            result["standalone"] = []
                            for slot in ("a", "b"):
                                saved = context.new_page()
                                try:
                                    name = f"{engine}-{width}-{slot}"
                                    response = saved.goto(
                                        url + "/saved/" + name, wait_until="load"
                                    )
                                    assert response and response.ok
                                    assert (
                                        response.headers["content-type"]
                                        == "image/svg+xml"
                                    )
                                    original = next(
                                        item
                                        for item in result["captures"]
                                        if item["name"] == name
                                    )
                                    # Parse and execute the actual saved XML SVG document.
                                    capture = saved.evaluate(
                                        """async request => {
                                            const harness=await import('/harness.js');
                                            return await harness.inspect(document,request);
                                        }""",
                                        {"labels": original["source_labels"]},
                                    )
                                    assert (
                                        capture["font_embedded"]
                                        and capture["font_loaded"]
                                    )
                                    assert [
                                        item["authored"] for item in capture["labels"]
                                    ] == [
                                        item["authored"]
                                        for item in original["geometry"]["labels"]
                                    ]
                                    result["standalone"].append(
                                        {"name": name, "geometry": capture}
                                    )
                                    saved.screenshot(
                                        path=str(output / f"{name}-standalone.png"),
                                        # Raw SVG XML has no HTML body. Playwright's
                                        # full-page helper waits for one indefinitely.
                                        full_page=False,
                                    )
                                finally:
                                    saved.close()
                                if width == 350 and slot == "a":
                                    image_page = context.new_page()
                                    image_page.set_default_timeout(15000)
                                    try:
                                        image_page.goto(url + "/", wait_until="load")
                                        image_capture = image_page.evaluate(
                                            """async name => {
                                                window.harness=await import('/harness.js');
                                                return await harness.mountSavedImage(name);
                                            }""",
                                            name,
                                        )
                                        result["image_embedding"] = image_capture
                                        reference_png = (
                                            image_page.frame_locator("#image-reference")
                                            .locator("svg")
                                            .screenshot(
                                                path=str(
                                                    output
                                                    / f"{name}-image-reference.png"
                                                )
                                            )
                                        )
                                        image_page.evaluate("harness.showSavedImage()")
                                        image_png = image_page.locator(
                                            "#saved-svg-image"
                                        ).screenshot(
                                            path=str(output / f"{name}-image.png")
                                        )
                                        image_page.evaluate(
                                            "harness.showSavedImage(true)"
                                        )
                                        negative_png = image_page.locator(
                                            "#saved-svg-negative-control"
                                        ).screenshot(
                                            path=str(
                                                output
                                                / f"{name}-image-negative-control.png"
                                            )
                                        )
                                        comparison = image_page.evaluate(
                                            """async ({reference,image,negative}) =>
                                                await harness.compareSavedImagePixels(reference,image,negative)""",
                                            {
                                                "reference": base64.b64encode(
                                                    reference_png
                                                ).decode("ascii"),
                                                "image": base64.b64encode(
                                                    image_png
                                                ).decode("ascii"),
                                                "negative": base64.b64encode(
                                                    negative_png
                                                ).decode("ascii"),
                                            },
                                        )
                                        image_capture["pixels"] = comparison
                                        assert comparison["negative_control"][
                                            "passed"
                                        ], (
                                            "saved SVG <img> MathML ink/unchanged-vector "
                                            "negative control failed: "
                                            + json.dumps(comparison)
                                        )
                                    finally:
                                        image_page.close()
                            assert not network, (
                                f"non-producer network requests: {network}"
                            )
                            assert not errors, errors
                            result["ok"] = True
                        except Exception:
                            page.screenshot(
                                path=str(output / f"{engine}-{width}-failure.png"),
                                full_page=True,
                            )
                            (output / f"{engine}-{width}-failure.html").write_text(
                                page.content()
                            )
                            report["failure_context"] = {
                                "engine": engine,
                                "width": width,
                                "network": network,
                                "page_errors": errors,
                                "geometry": page.evaluate(
                                    """() => [...document.querySelectorAll('iframe')].map(frame=>({
                                        title:frame.title,
                                        geometry:frame.contentWindow?.linnetGeometryDiagnostic || null,
                                    }))"""
                                ),
                            }
                            raise
                        finally:
                            context.close()
                finally:
                    browser.close()
        assert "symbolica" not in sys.modules
        report["ok"] = True
    except Exception as error:  # noqa: BLE001 -- preserve every failure and notify the finite producer.
        report["ok"] = False
        report["error"] = str(error)
        report["traceback"] = traceback.format_exc()
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous_alarm)
        if producer:
            try:
                post(
                    producer["url"] + "/finish",
                    {
                        "ok": report.get("ok", False),
                        "error": report.get("error"),
                        "matrix_jobs": len(report["matrix"]),
                    },
                )
            except (OSError, ValueError) as error:
                report["finish_error"] = str(error)
                report["ok"] = False
        report["finished_at"] = datetime.datetime.now(datetime.UTC).isoformat()
        (output / "browser-report.json").write_text(json.dumps(report, indent=2))
        print(
            json.dumps(
                {
                    "ok": report["ok"],
                    "mode": report["mode"],
                    "report": str(output / "browser-report.json"),
                    "error": report.get("error"),
                }
            )
        )
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
