"""Optional notebook and headless adapters for the native prepared scene."""

from functools import lru_cache
from threading import RLock

_widget_lock = RLock()


def available(module):
    """Optional display dependencies are absent only when their import fails."""
    from importlib import import_module

    try:
        import_module(module)
    except ImportError:
        return False
    return True


def scripts(templates, interactive):
    """Only Rust-authored templates supply executable notebook controls.

    Extract before inserting any user labels or snapshot SVG. The iframe itself
    cannot execute scripts; the parent runs only these exact trusted sources.
    """
    from html.parser import HTMLParser

    class Controls(HTMLParser):
        def __init__(self):
            super().__init__()
            self.active = False
            self.parts = []
            self.sources = [interactive]

        def handle_starttag(self, tag, attrs):
            self.active = tag == "script" and not dict(attrs).get("type")

        def handle_data(self, data):
            if self.active:
                self.parts.append(data)

        def handle_endtag(self, tag):
            if tag == "script" and self.active:
                self.sources.append("".join(self.parts))
                self.parts.clear()
                self.active = False

    parser = Controls()
    for template in templates:
        parser.feed(template)
    return list(dict.fromkeys(parser.sources))


def widget(view, request, javascript, completed=None, scripts=()):
    # Imports may release the GIL. lru_cache alone permits concurrent initial
    # misses to construct separate classes for the same source.
    with _widget_lock:
        cls = _widget_class(javascript)
    return cls(view, request, completed, scripts)


@lru_cache(maxsize=1)
def _widget_class(javascript):
    try:
        import anywidget
        import traitlets
    except ImportError as error:
        raise ImportError(
            "Notebook diagram measurement requires anywidget; "
            "install it with `python -m pip install anywidget`. "
            "Scripts can explicitly use measured(browser='chromium')."
        ) from error

    import base64
    import json

    # AnyWidget loads each source through a fresh Blob URL. Importing a stable
    # inner URL gives the browser one module owner for document-scoped caches.
    module_url = "data:text/javascript;base64," + base64.b64encode(
        javascript.encode("utf-8")
    ).decode("ascii")
    module = json.dumps(module_url)

    class DiagramWidget(anywidget.AnyWidget):
        _esm = f"export {{ default, measure }} from {module};"
        request = traitlets.Dict().tag(sync=True)
        rendered = traitlets.Dict().tag(sync=True)
        font = traitlets.Dict().tag(sync=True)
        scripts = traitlets.List(traitlets.Unicode()).tag(sync=True)

        def __init__(self, owner, initial, html, scripts):
            initial = dict(initial)
            self._font_css = initial.pop("font_css")
            super().__init__(request=initial, scripts=list(scripts))
            self._view = owner
            self._views = set()
            if html is not None:
                self.rendered = {
                    "id": self.request["id"],
                    "html": self._notebook_html(html),
                }
            self.on_msg(self._receive)

        def _receive(self, widget, content, buffers):
            action = content.get("action")
            if action == "attach":
                self._views.add(content["view"])
            elif action == "dispose":
                # Detaching a frontend does not close a retained widget's comm.
                self._views.discard(content["view"])
            elif action == "font" and self._view is not None:
                if content.get("key") == self.request["font_key"]:
                    self.font = {"key": self.request["font_key"], "css": self._font_css}
            elif action == "measure" and self._view is not None:
                receipt = content.get("measurements", {})
                if (
                    not isinstance(receipt, dict)
                    or not self.request.get("pending", True)
                    or receipt.get("id") != self.request["id"]
                ):
                    return
                try:
                    snapshot = self._view.measured(measurements=receipt)
                    self.request = {**self.request, "pending": False}
                    self.rendered = {
                        "id": self.request["id"],
                        "html": self._notebook_html(snapshot.to_html()),
                    }
                except Exception as error:  # noqa: BLE001 -- Report failures in the widget.
                    # A reentrant/concurrent completion may already have won.
                    # Keep the initial static snapshot on malformed receipts.
                    if not self.rendered.get("html"):
                        self.rendered = {"id": self.request["id"], "error": str(error)}

        def _notebook_html(self, html):
            # Only the exact block owned by native capture is transported as a
            # separate resource. Immutable exports retain their embedded font.
            if self._font_css:
                html = html.replace(
                    '<style data-linnet-font="true"><![CDATA['
                    + self._font_css
                    + "]]></style>",
                    "",
                )
                html = html.replace(
                    '<style data-linnet-font="true">' + self._font_css + "</style>", ""
                )
            return html

        def close(self):
            if getattr(self, "_view", None) is not None:
                self.on_msg(self._receive, remove=True)
                self._view = None
                self._font_css = ""
            super().close()

    return DiagramWidget


def headless(request, javascript):
    """Measure explicitly; a five-second browser deadline raises TimeoutError."""
    import asyncio

    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return _headless(request, javascript)

    # Sync Playwright owns an event loop. A notebook may already own this
    # thread's loop, so finish the explicit request in a temporary worker.
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(
        max_workers=1, thread_name_prefix="linnet-measure"
    ) as worker:
        return worker.submit(_headless, request, javascript).result()


def _headless(request, javascript):
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as error:
        raise ImportError(
            "Headless diagram measurement requires Playwright and Chromium; "
            "install `playwright` and run `python -m playwright install chromium`."
        ) from error

    with sync_playwright() as driver:
        browser = driver.chromium.launch()
        try:
            page = browser.new_page()
            # Every required byte is in the request/module. Never resolve a CDN.
            # The inline font needs no widget transfer allowance; bound browser
            # layout by the same five seconds used after notebook font transfer.
            page.route("**/*", lambda route: route.abort())
            receipt = page.evaluate(
                """async ({request, javascript}) => {
                    const url = URL.createObjectURL(new Blob([javascript],
                        {type: 'text/javascript'}));
                    let timer;
                    try {
                        const module = await import(url);
                        return await Promise.race([
                            module.measure(request),
                            new Promise((_, reject) => { timer = setTimeout(() => {
                                const error = new Error('Diagram measurement exceeded 5 seconds');
                                error.name = 'LinnetMeasurementTimeout';
                                reject(error);
                            }, 5000); }),
                        ]);
                    } catch (error) {
                        if (error.name === 'LinnetMeasurementTimeout')
                            return {__linnet_timeout: error.message};
                        throw error;
                    } finally { clearTimeout(timer); URL.revokeObjectURL(url); }
                }""",
                {"request": request, "javascript": javascript},
            )
            if "__linnet_timeout" in receipt:
                raise TimeoutError(
                    "Diagram measurement exceeded 5 seconds"
                ) from RuntimeError(receipt["__linnet_timeout"])
            return receipt
        finally:
            browser.close()
