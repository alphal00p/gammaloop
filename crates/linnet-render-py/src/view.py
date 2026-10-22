"""Optional notebook and headless adapters for the native prepared scene."""


def widget(view, request, javascript, completed=None):
    try:
        import anywidget
        import traitlets
    except ImportError as error:
        raise ImportError(
            "Notebook diagram measurement requires anywidget; "
            "install it with `python -m pip install anywidget`. "
            "Scripts can explicitly use measured(browser='chromium')."
        ) from error

    class DiagramWidget(anywidget.AnyWidget):
        _esm = javascript
        request = traitlets.Dict().tag(sync=True)
        rendered = traitlets.Dict().tag(sync=True)

        def __init__(self, owner, initial, html):
            super().__init__(request=initial)
            self._view = owner
            self._views = set()
            if html is not None:
                self.rendered = {"id": self.request["id"], "html": html}
            self.on_msg(self._receive)

        def _receive(self, widget, content, buffers):
            action = content.get("action")
            if action == "attach":
                self._views.add(content["view"])
            elif action == "dispose":
                # Detaching a frontend does not close a retained widget's comm.
                self._views.discard(content["view"])
            elif action == "measure" and self._view is not None:
                try:
                    snapshot = self._view.measured(measurements=content["measurements"])
                    self.rendered = {
                        "id": self.request["id"],
                        "html": snapshot.to_html(),
                    }
                except Exception as error:  # noqa: BLE001 -- Report failures in the widget.
                    self.rendered = {"id": self.request["id"], "error": str(error)}

        def close(self):
            if getattr(self, "_view", None) is not None:
                self.on_msg(self._receive, remove=True)
                self._view = None
            super().close()

    return DiagramWidget(view, request, completed)


def headless(request, javascript):
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
            page.route("**/*", lambda route: route.abort())
            return page.evaluate(
                """async ({request, javascript}) => {
                    const url = URL.createObjectURL(new Blob([javascript],
                        {type: 'text/javascript'}));
                    try {
                        const module = await import(url);
                        return await module.measure(request);
                    } finally { URL.revokeObjectURL(url); }
                }""",
                {"request": request, "javascript": javascript},
            )
        finally:
            browser.close()
