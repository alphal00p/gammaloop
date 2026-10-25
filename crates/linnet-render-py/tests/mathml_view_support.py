"""Executed inside the Rust integration binary with real native fixture owners."""
# ruff: noqa: F821 -- native, adapter_source, javascript, test_directory are injected by Rust.

import builtins
import contextlib
import copy
import hashlib
import json
import math
import os
import sys
import tempfile
import time
import types
import xml.etree.ElementTree as ET
from collections import Counter
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

SVG = "{http://www.w3.org/2000/svg}"
MATH = "{http://www.w3.org/1998/Math/MathML}"


@contextlib.contextmanager
def no_optional_imports():
    original = builtins.__import__
    attempted = []

    def guarded(name, *args, **kwargs):
        if name.split(".")[0] in {
            "typst",
            "playwright",
            "anywidget",
            "marimo",
            "symbolica",
        }:
            attempted.append(name)
            raise AssertionError(f"unexpected optional import: {name}")
        return original(name, *args, **kwargs)

    builtins.__import__ = guarded
    try:
        yield attempted
    finally:
        builtins.__import__ = original
    assert not attempted, attempted


def raises(exception, fragment, action):
    try:
        action()
    except exception as error:
        assert fragment.lower() in str(error).lower(), str(error)
    else:
        raise AssertionError(f"expected {exception}: {fragment}")


def protocol_receipt(request):
    # Deliberate protocol inputs, never used as browser geometry evidence.
    return {
        "id": request["id"],
        "metrics": [
            {"id": label["id"], "width": 35.25, "height": 19.5, "baseline": 12.25}
            for label in request["labels"]
        ],
    }


def geometry(element):
    """Compare all authored geometry, ignoring only the interactive decorator."""
    return (
        element.tag,
        tuple(
            sorted(
                (key, value)
                for key, value in element.attrib.items()
                if not key.startswith(("data-linnet-", "aria-"))
                and key not in {"href", "role", "tabindex"}
            )
        ),
        element.text if element.tag.startswith(MATH) else None,
        tuple(
            geometry(child)
            for child in element
            if child.tag not in {SVG + "style", SVG + "script", SVG + "title"}
        ),
    )


def verify_capture(view, marker, receipt):
    request = view.measurement_request()
    result = view.measured(measurements=receipt)
    svg = result.to_svg()
    assert svg == view.to_svg() == view.measured().to_svg()
    assert result.to_html() == view.to_html()
    assert result.to_svg_pages() == [svg]
    assert "<!--linnet-svg-" not in result.to_html()
    assert geometry(ET.fromstring(svg)) == geometry(
        ET.fromstring(native.native_oracle(marker, json.dumps(receipt)))
    ), "retained scene changed native geometry or source"

    root = ET.fromstring(svg)
    objects = root.findall(".//" + SVG + "foreignObject")
    assert len(objects) == 9  # title, three nodes, five edge labels
    counts = Counter()
    for obj in objects:
        mathml = obj.find(".//" + MATH + "math")
        assert mathml is not None
        matches = [
            label["id"]
            for label in request["labels"]
            if geometry(mathml)
            == geometry(
                ET.fromstring(
                    '<math xmlns="http://www.w3.org/1998/Math/MathML">'
                    + label["mathml"]
                    + "</math>"
                )
            )
        ]
        assert len(matches) == 1
        label_id = matches[0]
        counts[label_id] += 1
        metric = receipt["metrics"][label_id]
        assert float(obj.attrib["width"]) == metric["width"]
        assert float(obj.attrib["height"]) == metric["height"]
        assert "overflow:visible" not in obj.attrib.get("style", "").replace(" ", "")
        assert request["labels"][label_id]["mathml"] in svg
    assert counts == {0: 1, 1: 1, 2: 1, 3: 1, 4: 4, 5: 1}, counts

    targets = root.findall(".//*[@data-linnet-kind]")
    identities = {
        kind: {
            int(target.attrib["data-linnet-id"])
            for target in targets
            if target.attrib["data-linnet-kind"] == kind
        }
        for kind in ("node", "edge", "halfedge")
    }
    assert identities == {
        "node": {0, 1, 2},
        "edge": set(range(5)),
        "halfedge": set(range(8)),
    }
    endpoints = {0: (0, 0), 1: (0, 1), 2: (0, 1), 3: (None, 0), 4: (1, None)}
    for target in targets:
        detail = json.loads(target.attrib["data-linnet-detail"])
        assert detail["fixture"] == marker
        kind = target.attrib["data-linnet-kind"]
        identity = int(target.attrib["data-linnet-id"])
        if kind == "edge":
            assert (detail["source"], detail["sink"]) == endpoints[identity]
            assert detail["particle"] == f"p{identity}"
        elif kind == "node":
            assert set(detail["edges"]) == (
                {0, 1, 2, 3} if identity == 0 else {1, 2, 4} if identity == 1 else set()
            )
    assert "#246813" in svg and "#edf3e7" in svg
    assert any(
        [
            float(value)
            for value in node.attrib.get("stroke-dasharray", "")
            .replace(",", " ")
            .split()
        ]
        == [3.0, 2.0]
        for node in root.iter()
    )
    assert "data:font/woff2;base64," in svg and "STIX" in svg
    return result


def verify_native():
    with no_optional_imports():
        view = native.fresh_view("protocol")
        request = view.measurement_request()
        assert request == view.measurement_request()
        assert "points_per_css_pixel" not in request
        assert request["coordinate_system"] == "svg-user-units"
        assert request["font_style"]["fontSize"] == "9px"
        assert "data:font/woff2;base64," in request["font_css"]
        assert [label["id"] for label in request["labels"]] == list(range(6))
        assert "Native protocol" in request["labels"][0]["mathml"]
        assert all("<math" not in label["mathml"] for label in request["labels"])
        raises(RuntimeError, "measurement is pending", view.measured)
        assert "<svg" in view.to_svg()
        assert "<svg" in view.to_html()
        assert len(view.to_svg_pages()) == 1
        raises(RuntimeError, "measurement is pending", view.measured)
        raises(ValueError, "chromium", lambda: view.measured(browser="webkit"))
        valid = protocol_receipt(request)
        raises(
            ValueError,
            "not both",
            lambda: view.measured(measurements=valid, browser="chromium"),
        )
        invalid = []
        stale = copy.deepcopy(valid)
        stale["id"] = "stale-request"
        invalid.append((stale, "stale"))
        missing = copy.deepcopy(valid)
        missing["metrics"].pop()
        invalid.append((missing, "every requested label"))
        extra = copy.deepcopy(valid)
        extra["metrics"].append(copy.deepcopy(extra["metrics"][0]))
        invalid.append((extra, "every requested label"))
        for field, value in (
            ("id", 1),
            ("width", -1),
            ("height", -1),
            ("baseline", 20),
            ("baseline", -1),
            ("width", math.nan),
            ("height", math.inf),
            ("baseline", -math.inf),
        ):
            bad = copy.deepcopy(valid)
            bad["metrics"][0][field] = value
            invalid.append((bad, "invalid" if not math.isfinite(value) else "label"))
        unknown = copy.deepcopy(valid)
        unknown["metrics"][0]["guessed"] = True
        invalid.append((unknown, "unknown field"))
        for bad, fragment in invalid:
            raises(
                ValueError, fragment, lambda bad=bad: view.measured(measurements=bad)
            )
            raises(RuntimeError, "pending", view.measured)
        result = verify_capture(view, "protocol", valid)
        changed = copy.deepcopy(valid)
        changed["metrics"][0]["width"] += 1
        raises(
            ValueError,
            "different label measurements",
            lambda: view.measured(measurements=changed),
        )
        assert result.to_svg() == view.to_svg()
        with tempfile.TemporaryDirectory(prefix="linnet-mathml-") as directory:
            directory = Path(directory)
            for suffix, text in (("svg", result.to_svg()), ("html", result.to_html())):
                path = directory / ("native." + suffix)
                result.save(path)
                assert path.read_text() == text
            for suffix in ("pdf", "png", "typ"):
                path = directory / ("guarded." + suffix)
                raises(
                    (RuntimeError, ValueError),
                    "outlined",
                    lambda path=path: result.save(path),
                )
                assert not path.exists()
        captured = native.DiagramRender  # registration includes the frozen render API
        assert isinstance(result, captured)
        raises(AttributeError, "", lambda: setattr(result, "geometry", "mutable"))
    # Simulate missing optional packages through import failures. No browser or
    # compiler is started, and no measurement or final SVG is mocked.
    original = builtins.__import__
    attempted = []

    def missing_optional(name, *args, **kwargs):
        if name.split(".")[0] in {"anywidget", "playwright", "typst"}:
            attempted.append(name)
            raise ImportError(f"test intentionally omits {name}")
        return original(name, *args, **kwargs)

    adapter = types.ModuleType("mathml_test_adapter")
    exec(adapter_source, adapter.__dict__)  # noqa: S102 -- trusted source embedded by Rust.
    pending = native.fresh_view("opt-in")
    builtins.__import__ = missing_optional
    try:
        raises(RuntimeError, "pending", pending.measured)
        assert not attempted
        raises(ImportError, "anywidget", pending.widget)
        assert attempted == ["anywidget"]
        attempted.clear()
        raises(ImportError, "Playwright", lambda: pending.measured(browser="chromium"))
        assert attempted == ["playwright.sync_api"]
        attempted.clear()
        raises(ImportError, "Playwright", lambda: adapter.headless(request, javascript))
        assert attempted == ["playwright.sync_api"]
    finally:
        builtins.__import__ = original


def verify_outlined():
    # This opt-in test requires the optional Typst package. It is separate from
    # receipt/browser tests so measurement can never acquire a compiler fallback.
    view = native.fresh_view("outlined")
    snapshot = view.outlined()
    assert "linnet-math-label" not in snapshot.to_svg()
    assert "<foreignObject" not in snapshot.to_svg()
    assert view.outlined().to_svg() == snapshot.to_svg()
    raises(RuntimeError, "pending", view.measured)
    with tempfile.TemporaryDirectory(prefix="linnet-outlined-") as directory:
        for suffix, prefix in (("pdf", b"%PDF"), ("png", b"\x89PNG\r\n\x1a\n")):
            path = Path(directory) / ("native." + suffix)
            snapshot.save(path)
            assert path.read_bytes().startswith(prefix)


def serve_browser():
    output = Path(os.environ["LINNET_MATHML_BROWSER_DIR"]).resolve()
    allowed = Path(
        os.environ.get(
            "LINNET_MATHML_BROWSER_ROOT",
            Path(test_directory).parents[2] / ".local/graph-mathml-labels/browser",
        )
    )
    assert output.is_relative_to(allowed.resolve()), output
    output.mkdir(parents=True, exist_ok=True)
    python = {
        "version_info": list(sys.version_info),
        "version": sys.version,
        "executable": sys.executable,
        "prefix": sys.prefix,
        "base_prefix": sys.base_prefix,
        "pyo3_config_file": os.environ.get("PYO3_CONFIG_FILE"),
        "python_home": os.environ.get("PYTHONHOME"),
    }
    engines = os.environ.get("LINNET_MATHML_BROWSER_ENGINES", "chromium,webkit").split(
        ","
    )
    assert engines and set(engines) <= {"chromium", "webkit"}
    jobs = {
        f"{engine}-{width}-{slot}": native.fresh_view(f"{engine}-{width}-{slot}")
        for engine in engines
        for width in (350, 1400)
        for slot in ("a", "b")
    }
    requests = {name: view.measurement_request() for name, view in jobs.items()}
    widgets = {name: view.widget() for name, view in jobs.items()}
    assert all(widget.scripts for widget in widgets.values())
    assert len({request["id"] for request in requests.values()}) == len(jobs)
    completed = {}
    finish = None

    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.0"
        handled = 0

        def setup(self):
            super().setup()
            self.connection.settimeout(2)
            type(self).handled += 1

        def log_message(self, *_args):
            pass

        def reply(self, content, kind="application/json", status=200):
            payload = (
                content.encode()
                if isinstance(content, str)
                else json.dumps(content).encode()
            )
            self.send_response(status)
            self.send_header("Content-Type", kind)
            self.send_header("Content-Length", str(len(payload)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(payload)

        def do_GET(self):
            if self.path == "/":
                self.reply(
                    '<!doctype html><meta charset="utf-8"><style>body{margin:0}.host{margin:0}</style><main id="hosts"></main>',
                    "text/html",
                )
            elif self.path == "/manifest":
                self.reply(
                    {
                        "requests": requests,
                        "scripts": {
                            name: widget.scripts for name, widget in widgets.items()
                        },
                        "native": True,
                        "python": python,
                        "source_sha256": hashlib.sha256(
                            javascript.encode()
                        ).hexdigest(),
                    }
                )
            elif self.path == "/view.js":
                self.reply(javascript, "text/javascript")
            elif self.path == "/harness.js":
                self.reply(
                    (Path(test_directory) / "mathml_view_browser.js").read_text(),
                    "text/javascript",
                )
            elif (
                self.path.startswith("/saved/")
                and self.path.removeprefix("/saved/") in completed
            ):
                self.reply(
                    (output / (self.path.removeprefix("/saved/") + ".svg")).read_text(),
                    "image/svg+xml",
                )
            else:
                self.reply({"error": "unknown local fixture route"}, status=404)

        def do_POST(self):
            nonlocal finish
            try:
                length = int(self.headers.get("Content-Length", "0"))
                assert 0 < length <= 2_000_000
                payload = json.loads(self.rfile.read(length))
                if self.path == "/finish":
                    finish = payload
                    self.reply({"accepted": True})
                    return
                name = self.path.removeprefix("/receipt/")
                assert self.path.startswith("/receipt/") and name in jobs
                assert name not in completed, "duplicate native completion"
                widget = widgets[name]
                assert widget.request["pending"]
                widget._receive(
                    widget, {"action": "measure", "measurements": payload}, []
                )
                assert widget.request["pending"] is False
                with no_optional_imports():
                    result = verify_capture(jobs[name], name, payload)
                    result.save(output / (name + ".svg"))
                    result.save(output / (name + ".html"))
                (output / (name + ".receipt.json")).write_text(
                    json.dumps(payload, indent=2)
                )
                completed[name] = payload
                self.reply(
                    {
                        "id": requests[name]["id"],
                        "html": result.to_html(),
                        "svg": result.to_svg(),
                        "native": True,
                        "request": widget.request,
                    }
                )
            except (
                AssertionError,
                ValueError,
                RuntimeError,
                TypeError,
                OSError,
            ) as error:
                self.reply({"error": f"{type(error).__name__}: {error}"}, status=400)

    started = time.monotonic()
    with HTTPServer(("127.0.0.1", 0), Handler) as server:
        server.timeout = 0.5
        # Bound slow or abandoned clients as well as total requests/time.
        ready = {
            "url": f"http://127.0.0.1:{server.server_port}",
            "jobs": list(jobs),
            "pid": os.getpid(),
            "max_seconds": 240,
            "python": python,
        }
        (output / "ready.json").write_text(json.dumps(ready, indent=2))
        print("LINNET_MATHML_READY " + json.dumps(ready), flush=True)
        while (
            finish is None
            and time.monotonic() - started < 240
            and Handler.handled < 256
        ):
            server.handle_request()
    for widget in widgets.values():
        widget.close()
    summary = {
        "completed": list(completed),
        "expected": list(jobs),
        "finish": finish,
        "elapsed_seconds": time.monotonic() - started,
        "python": python,
    }
    (output / "producer-summary.json").write_text(json.dumps(summary, indent=2))
    assert finish is not None, "browser producer timed out or exhausted request budget"
    assert finish.get("ok") is True, finish
    assert completed.keys() == jobs.keys(), summary
