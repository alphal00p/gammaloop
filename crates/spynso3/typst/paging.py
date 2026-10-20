"""Compile bounded tensor pages for Symbolica's shared MathML viewer."""

import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from html import escape
from pathlib import Path

MAX_HTML = 256 * 1024
MAX_ELEMENTS = 10_000
NOTEBOOK_STYLE = ""

_MAIN = """#import "notation.typ" as n
#set page(width: auto, height: auto, margin: 4pt)
#let tree = cbor(read("tree.cbor", encoding: none))
#let settings = {settings}
#let notation = n.merge-notation(n.default-notation(settings: settings),
  (fallback-node: ctx => if ctx.kind == "omission" {{ text(ctx.node.at("label", default: "⋯")) }} else if ctx.kind == "inspection" {{ ctx.node.items.map(ctx.render-visual).join($,quad$) }} else {{ ctx.default() }}))
$ #n.render(tree, notation: notation) $
#if tree.at("page", default: none) != none {{
  $ #n.render((..tree, root: tree.page), notation: notation) $
}}
"""
# Output stays in a file, so a bad custom notation cannot flood a PIPE buffer.
_WORKER = """import sys, typst
from pathlib import Path
root=Path(sys.argv[1])
try:
    typst.compile(str(root/'main.typ'), output=str(root/'out.html'), format='html', root=str(root), pretty=False)
except Exception as error:  # noqa: BLE001 - display errors belong inside the viewer
    print(getattr(error, 'diagnostic', str(error)), file=sys.stderr)
    sys.exit(1)
"""


def compile_page(payload):
    if sys.platform in ("emscripten", "wasi"):
        # Browser Python cannot launch a worker process. Reuse the embedded
        # export compiler with the same bounded tree and output budgets.
        document = bytes(
            compile_typst(  # noqa: F821 - supplied by the native module
                _MAIN.format(settings=payload["settings"]),
                "html",
                payload["notation"],
                payload["tree"],
            )
        )
        if len(document) > MAX_HTML:
            return None
        html = document.decode()
    else:
        with tempfile.TemporaryDirectory(prefix="spenso-page-") as directory:
            root = Path(directory)
            (root / "tree.cbor").write_bytes(payload["tree"])
            (root / "notation.typ").write_text(payload["notation"])
            (root / "render.typ").write_text(payload["render"])
            (root / "main.typ").write_text(_MAIN.format(settings=payload["settings"]))
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

    markers = dict(payload.get("markers", ()))
    seen = set()

    # Wrap only a standalone outer sum. Flex layout inside a MathML product,
    # fraction or fence breaks its native geometry and can overlap neighbours.
    # Nested sums remain available through the native subexpression controls.
    def wrap(match):
        root = ET.fromstring(match.group())
        # The Typst source contains ordinary text tokens. Only after rendering
        # replace their MathML leaves with the trusted shared-viewer markers.
        for parent in root.iter():
            for index, child in enumerate(list(parent)):
                if not len(child) and child.text in markers:
                    replacement = ET.fromstring(markers[child.text])
                    replacement.tail = child.tail
                    parent[index] = replacement
                    seen.add(child.text)
        row = root
        while True:
            content = [
                child
                for child in row
                if child.tag not in ("annotation", "annotation-xml")
            ]
            if len(content) != 1 or content[0].tag not in (
                "mrow",
                "mstyle",
                "semantics",
            ):
                break
            row = content[0]
        children = list(row)
        signs = [
            i
            for i, child in enumerate(children)
            if i and child.tag == "mo" and child.text in ("+", "−", "-")
        ]
        if len(signs) < 8 or any(
            child.tag == "mo" and child.text in ("(", ")", "[", "]", "|")
            for child in children
        ):
            return ET.tostring(root, encoding="unicode")
        groups = []
        group = ET.Element("mrow", {"style": "flex:0 0 auto"})
        for i, child in enumerate(children):
            if i in signs:
                groups.append(group)
                group = ET.Element("mrow", {"style": "flex:0 0 auto"})
                child.set("form", "infix")
            group.append(child)
        groups.append(group)
        row[:] = groups
        row.set("data-spenso-sum", "")
        row.set(
            "style",
            "display:inline-flex;flex-wrap:wrap;align-items:baseline;max-width:min(650px,80vw);row-gap:.5em",
        )
        return ET.tostring(root, encoding="unicode")

    html = re.sub(r"<math\b.*?</math>", wrap, html, flags=re.DOTALL)
    if seen != markers.keys():
        raise RuntimeError("Renderer did not preserve omission markers")
    math = re.findall(r"<math\b.*?</math>", html, flags=re.DOTALL)
    if not math:
        raise RuntimeError("Renderer returned no MathML")
    if sum(sum(1 for _ in ET.fromstring(m).iter()) for m in math) > MAX_ELEMENTS:
        return None
    # Scope Typst's equation styling to the viewer, sharing the standard font.
    body = re.search(r"<body[^>]*>(.*?)</body>", html, flags=re.DOTALL)
    fragment = body.group(1) if body else "".join(math)
    if payload.get("selection") and len(math) == 2:
        selection = escape(payload["selection"])
        fragment = (
            '<p class="spenso-page-heading">Expression structure · numbered markers identify omitted parts</p>'
            "<div data-spenso-math>" + math[0] + "</div>"
            '<p class="spenso-page-heading">' + selection + " · displayed terms</p>"
            "<div data-spenso-math>" + math[1] + "</div>"
        )
    else:
        fragment = "<div data-spenso-math>" + fragment + "</div>"
    result = "<style>" + NOTEBOOK_STYLE + "</style>" + fragment
    return result if len(result.encode()) <= MAX_HTML else None
