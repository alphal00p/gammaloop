"""Compile bounded tensor notation leaves for Symbolica's MathML renderer."""

import atexit
import json
import queue
import subprocess
import sys
import tempfile
import threading
import time
import xml.etree.ElementTree as ET
from pathlib import Path

MAX_HTML = 256 * 1024
MAX_ELEMENTS = 1000

_MAIN = """#import "notation.typ" as n
#set page(width: auto, height: auto, margin: 0pt)
#let tree = cbor(read("tree.cbor", encoding: none))
#let settings = {settings}
#let notation = n.default-notation(settings: settings)
$ #n.render(tree, notation: notation, power-base: {power_base}) $
"""
_WORKER = """import json, os, sys
from pathlib import Path
protocol = os.fdopen(os.dup(sys.stdout.fileno()), 'w', buffering=1)
os.dup2(sys.stderr.fileno(), sys.stdout.fileno())
from symbolica.community.tensor_native import _compile_typst
print('ready', file=protocol, flush=True)
for request in sys.stdin:
    root = Path(json.loads(request))
    try:
        (root/'out.html').write_bytes(bytes(_compile_typst((root/'main.typ').read_text(), 'html', (root/'notation.typ').read_text(), (root/'tree.cbor').read_bytes())))
        print('ok', file=protocol, flush=True)
    except Exception as error:  # noqa: BLE001 - display errors belong inside the viewer
        (root/'errors').write_text(str(error)[:2048])
        print('error', file=protocol, flush=True)
"""

_worker = None
_worker_errors = None
_worker_replies = None
_worker_reader = None
_worker_lock = threading.Lock()


def _stop_worker():
    """Release the module-owned compiler; callers serialize active requests."""
    global _worker, _worker_errors, _worker_replies, _worker_reader
    worker, errors, reader = _worker, _worker_errors, _worker_reader
    _worker = _worker_errors = _worker_replies = _worker_reader = None
    try:
        if worker is not None:
            if worker.poll() is None:
                worker.kill()
            worker.wait()
    finally:
        if worker is not None:
            worker.stdin.close()
            if reader is not None:
                reader.join(timeout=1)
            worker.stdout.close()
        if errors is not None:
            errors.close()


def _read_reply(deadline):
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise TimeoutError("Tensor compiler exceeded the page deadline")
    try:
        message = _worker_replies.get(timeout=remaining)
    except queue.Empty:
        raise TimeoutError("Tensor compiler exceeded the page deadline") from None
    if not message:
        _worker_errors.seek(0)
        detail = _worker_errors.read(2048).decode(errors="replace")
        raise RuntimeError(detail or "Tensor compiler closed its output")
    if len(message) <= 64 and message.endswith(b"\n"):
        reply = message.decode().strip()
        if reply in ("ready", "ok", "error"):
            return reply
    raise RuntimeError("Tensor compiler returned an invalid reply")


def _compile_native(root, timeout):
    global _worker, _worker_errors, _worker_replies, _worker_reader
    deadline = time.monotonic() + min(float(timeout), 10.0)
    remaining = deadline - time.monotonic()
    if remaining <= 0 or not _worker_lock.acquire(timeout=remaining):
        raise TimeoutError("Tensor compiler queue exceeded the page deadline")
    try:
        try:
            if _worker is not None and _worker.poll() is not None:
                _stop_worker()
            if _worker is None:
                _worker_errors = tempfile.TemporaryFile()  # noqa: SIM115 - closed with the worker
                _worker = subprocess.Popen(
                    [sys.executable, "-c", _WORKER],
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=_worker_errors,
                    bufsize=0,
                )
                _worker_replies = queue.Queue(maxsize=2)
                stdout, replies = _worker.stdout, _worker_replies

                def collect_replies():
                    try:
                        while True:
                            message = stdout.readline(65)
                            try:
                                replies.put_nowait(message)
                            except queue.Full:
                                return
                            if (
                                not message
                                or not message.endswith(b"\n")
                                or len(message) > 64
                            ):
                                return
                    except (OSError, ValueError):
                        try:
                            replies.put_nowait(b"")
                        except queue.Full:
                            pass

                _worker_reader = threading.Thread(target=collect_replies, daemon=True)
                _worker_reader.start()
                if _read_reply(deadline) != "ready":
                    raise RuntimeError("Tensor compiler did not become ready")
            _worker_errors.seek(0)
            _worker_errors.truncate()
            _worker.stdin.write((json.dumps(str(root)) + "\n").encode())
            _worker.stdin.flush()
            reply = _read_reply(deadline)
            if reply == "error":
                with (root / "errors").open("rb") as errors:
                    raise RuntimeError(errors.read(2048).decode(errors="replace"))
            if reply != "ok":
                raise RuntimeError("Tensor compiler returned an unexpected reply")
        except BaseException:
            _stop_worker()
            raise
    finally:
        _worker_lock.release()


atexit.register(_stop_worker)


def compile_leaf(payload):
    deadline = time.monotonic() + min(float(payload.get("timeout", 10)), 10.0)
    source = _MAIN.format(
        settings=payload["settings"], power_base=str(payload["power_base"]).lower()
    )
    if sys.platform in ("emscripten", "wasi"):
        document = bytes(
            compile_typst(  # noqa: F821 - supplied by the native module
                source, "html", payload["notation"], payload["tree"]
            )
        )
        if len(document) > MAX_HTML:
            return None
        html = document.decode()
    else:
        with tempfile.TemporaryDirectory(prefix="spenso-leaf-") as directory:
            root = Path(directory)
            (root / "tree.cbor").write_bytes(payload["tree"])
            (root / "notation.typ").write_text(payload["notation"])
            (root / "main.typ").write_text(source)
            _compile_native(root, deadline - time.monotonic())
            path = root / "out.html"
            if path.stat().st_size > MAX_HTML:
                return None
            html = path.read_text()
    import re

    match = re.search(r"<math\b.*?</math>", html, flags=re.DOTALL)
    if match is None:
        raise RuntimeError("Renderer returned no MathML")
    math = ET.fromstring(match.group())
    if sum(1 for _ in math.iter()) > MAX_ELEMENTS:
        return None
    # Symbolica owns outer MathML, arithmetic rows, paths and omissions.
    fragment = "".join(ET.tostring(child, encoding="unicode") for child in math)
    return fragment if len(fragment.encode()) <= 4096 else None
