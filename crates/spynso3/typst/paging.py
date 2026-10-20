"""Compile bounded tensor notation leaves for Symbolica's MathML renderer."""

import subprocess
import sys
import tempfile
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
_WORKER = """import sys, typst
from pathlib import Path
root=Path(sys.argv[1])
try:
    typst.compile(str(root/'main.typ'), output=str(root/'out.html'), format='html', root=str(root), pretty=False)
except Exception as error:  # noqa: BLE001 - display errors belong inside the viewer
    print(getattr(error, 'diagnostic', str(error)), file=sys.stderr)
    sys.exit(1)
"""


def compile_leaf(payload):
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
            (root / "render.typ").write_text(payload["render"])
            (root / "main.typ").write_text(source)
            with (root / "errors").open("wb") as errors:
                job = subprocess.run(
                    [sys.executable, "-c", _WORKER, directory],
                    stdout=errors,
                    stderr=errors,
                    timeout=payload.get("timeout", 10),
                    check=False,
                )
            if job.returncode:
                with (root / "errors").open("rb") as errors:
                    raise RuntimeError(errors.read(2048).decode(errors="replace"))
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
