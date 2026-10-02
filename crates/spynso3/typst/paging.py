"""Live, bounded MathML pages shared by marimo and Jupyter.

Only the native source owns the expression. No complete expression or hidden
page is serialized into widget state or exported notebook HTML.
"""

import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from collections import OrderedDict
from html import escape
from pathlib import Path

MAX_HTML = 256 * 1024
MAX_ELEMENTS = 10_000
WIDGET_ESM = ""  # supplied by the native module
NOTEBOOK_STYLE = ""

_MAIN = """#import "notation.typ" as n
#set page(width: auto, height: auto, margin: 4pt)
#let tree = cbor(read("tree.cbor", encoding: none))
#let settings = {settings}
#let notation = n.merge-notation(n.default-notation(settings: settings),
  (fallback-node: ctx => if ctx.kind == "omission" {{ math.upright(ctx.node.at("label", default: "⋯")) }} else if ctx.kind == "inspection" {{ ctx.node.items.map(ctx.render-visual).join($,quad$) }} else {{ ctx.default() }}))
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

    # Group only binary additions at the row's own level. Internal products,
    # fractions and scripted labels retain native MathML layout.
    def wrap(match):
        root = ET.fromstring(match.group())
        for row in list(root.iter()):
            if row.tag not in ("math", "mrow"):
                continue
            children = list(row)
            signs = [
                i
                for i, child in enumerate(children)
                if i and child.tag == "mo" and child.text in ("+", "−", "-")
            ]
            if len(signs) < 8:
                continue
            # A row containing literal fences is not a bare additive sequence.
            if any(
                c.tag == "mo" and c.text in ("(", ")", "[", "]", "|") for c in children
            ):
                continue
            groups = []
            group = ET.Element("mrow")
            for i, child in enumerate(children):
                if i in signs:
                    groups.append(group)
                    group = ET.Element("mrow")
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
            '<p class="spenso-page-heading">Expression structure · numbered ellipses mark sums</p>'
            "<div data-spenso-math>" + math[0] + "</div>"
            '<p class="spenso-page-heading">' + selection + " · displayed terms</p>"
            "<div data-spenso-math>" + math[1] + "</div>"
        )
    else:
        fragment = "<div data-spenso-math>" + fragment + "</div>"
    result = "<style>" + NOTEBOOK_STYLE + "</style>" + fragment
    return result if len(result.encode()) <= MAX_HTML else None


class Pager:
    def __init__(self, source, page_size=25):
        self._source = source
        self._size = page_size
        self._cache = OrderedDict()
        self._focus = 0
        self._cursor = None
        self._previous = []
        self._breadcrumbs = []
        self._widget = None
        self._closed = False
        self._views = set()

    def _page(self):
        if self._closed:
            raise RuntimeError("This display has been closed")
        key = (self._focus, self._cursor, self._size)
        if key in self._cache:
            self._cache.move_to_end(key)
            return self._cache[key]
        limit = self._size
        deadline = time.monotonic() + 10
        while True:
            payload = self._source.page(self._focus, self._cursor, limit)
            try:
                payload["timeout"] = deadline - time.monotonic()
                if payload["timeout"] <= 0:
                    raise TimeoutError("Page rendering exceeded the 10-second deadline")
                html = compile_page(payload)
            except Exception as error:  # noqa: BLE001 - display errors belong inside the viewer
                # A timeout/error never tries the unbounded formatter.
                html = (
                    '<p role="alert">Page rendering failed: '
                    + escape(str(error)[:2048])
                    + "</p>"
                )
            if html is not None:
                break
            if limit == 1:
                payload = self._source.page(self._focus, self._cursor, 1, True)
                html = "<p>This term exceeds the display budget. Use the subexpression controls to inspect it.</p>"
                break
            limit = max(1, limit // 2)
        page = {
            k: payload[k] for k in ("start", "end", "total", "next", "holes", "unit")
        }
        page["html"] = html
        page["selection"] = payload.get("selection")
        expected = min(self._size, page["total"] - page["start"])
        page["budget"] = (
            ("rendered output" if limit < self._size else "expression complexity")
            if page["end"] - page["start"] < expected
            else None
        )
        self._cache[key] = page
        while len(self._cache) > 3:
            self._cache.popitem(last=False)
        return page

    def _state(self):
        page = dict(self._page())
        page.update(
            previous=bool(self._previous),
            page_size=self._size,
            breadcrumbs=[label for _, _, _, _, label in self._breadcrumbs],
        )
        return page

    def _action(self, message):
        page = self._page()
        action = message.get("action")
        if action == "next" and page["next"] is not None:
            self._previous.append(self._cursor)
            self._cursor = page["next"]
        elif action == "previous" and self._previous:
            self._cursor = self._previous.pop()
        elif action == "size":
            size = message.get("value")
            if size not in (25, 100, 250, 500):
                raise ValueError("Invalid page size")
            self._size = size
            self._cursor = None
            self._previous = []
        elif action == "open":
            target = message.get("value")
            links = dict(page["holes"])
            if target not in links or target == self._focus:
                raise ValueError("Unknown subexpression")
            self._breadcrumbs.append(
                (self._focus, self._cursor, self._previous, self._size, links[target])
            )
            self._focus, self._cursor, self._previous = target, None, []
        elif action == "back" and self._breadcrumbs:
            self._focus, self._cursor, self._previous, self._size, _ = (
                self._breadcrumbs.pop()
            )
        return self._state()

    def close(self):
        self._closed = True
        self._cache.clear()
        self._previous.clear()
        self._breadcrumbs.clear()
        self._source = None
        if self._widget is not None:
            widget, self._widget = self._widget, None
            widget.on_msg(self._on_message, remove=True)
            widget.close()

    def _on_message(self, widget, content, buffers):
        if content.get("action") == "attach":
            self._views.add(content.get("view"))
            widget.page = dict(self._state(), connected=content.get("view"))
            return
        if content.get("action") == "dispose":
            self._views.discard(content.get("view"))
            if not self._views:
                self._cache.clear()
                widget.on_msg(self._on_message, remove=True)
                widget.close()
                self._widget = None
            return
        request = content.get("request")
        try:
            widget.page = dict(self._action(content), request=request)
        except Exception as error:  # noqa: BLE001 - display errors belong inside the viewer
            widget.page = dict(self._state(), request=request, error=str(error)[:2048])

    def _get_widget(self):
        if self._widget is None:
            import anywidget
            import traitlets

            class MathPages(anywidget.AnyWidget):
                _esm = WIDGET_ESM
                page = traitlets.Dict().tag(sync=True)

            self._widget = MathPages(page=self._state())
            self._widget.on_msg(self._on_message)
        return self._widget

    def _display_(self):
        import marimo as mo

        try:
            return mo.ui.anywidget(self._get_widget())
        except ImportError:
            return mo.Html(self._repr_html_())

    def _repr_mimebundle_(self, include=None, exclude=None):
        def wanted(mime):
            return (include is None or mime in include) and (
                exclude is None or mime not in exclude
            )

        result = {}
        if wanted("text/plain"):
            result["text/plain"] = repr(self)
        if wanted("text/html"):
            result["text/html"] = self._repr_html_()
        if wanted("application/vnd.jupyter.widget-view+json"):
            try:
                bundle = self._get_widget()._repr_mimebundle_(
                    include=include, exclude=exclude
                )
                # Anywidget follows IPython's (data, metadata) bundle protocol;
                # ipywidgets versions also accept a bare data dictionary.
                if isinstance(bundle, tuple):
                    bundle = bundle[0]
                result.update({k: v for k, v in (bundle or {}).items() if wanted(k)})
                if wanted("text/plain"):
                    result["text/plain"] = repr(self)
            except ImportError:
                pass
        return result

    def _repr_html_(self):
        page = self._page()
        label = (
            "{} {}–{} of {}".format(
                page["unit"], page["start"] + 1, page["end"], page["total"]
            )
            if page["total"]
            else "Expression preview"
        )
        return (
            "<div><p>"
            + label
            + " · Preview; omitted portions are marked ⋯.</p>"
            + page["html"]
            + "<p>Navigation requires a live notebook with symbolica[notebook-display]. Full export is available through to_html()/to_svg().</p></div>"
        )

    def __repr__(self):
        return "Tensor MathML viewer (bounded preview; use a live notebook to navigate)"
