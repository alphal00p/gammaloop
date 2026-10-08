"""Export real native displays with Chromium measurements for browser CI."""
# ruff: noqa: F821 -- test_directory is injected by the Rust producer.

import importlib.util
import os
from pathlib import Path

from playwright.sync_api import sync_playwright

root = Path(os.environ["FEYNKIT_BROWSER_FIXTURES"])
javascript = (
    Path(test_directory).parents[1] / "linnet-render-py/src/view.js"
).read_text(encoding="utf-8")
with sync_playwright() as driver:
    browser = driver.chromium.launch(executable_path=os.environ["LINNET_TEST_CHROMIUM"])
    try:
        page = browser.new_page()
        page.route("**/*", lambda route: route.abort())

        def measured(view):
            receipt = page.evaluate(
                """async ({request, javascript}) => {
                    const url = URL.createObjectURL(new Blob([javascript], {type:'text/javascript'}));
                    try {
                        const module = await import(url);
                        return await module.measure(request);
                    } finally {
                        URL.revokeObjectURL(url);
                    }
                }""",
                {"request": view.measurement_request(), "javascript": javascript},
            )
            return view.measured(measurements=receipt)

        for name, entry, directory in (
            ("installed_ltd_display", "main", "ltd"),
            ("installed_energy_scoped_displays", "check", "scoped"),
            ("installed_generalized_cff", "main", "generalized"),
        ):
            spec = importlib.util.spec_from_file_location(
                name, Path(test_directory) / f"{name}.py"
            )
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            # Browser evidence uses actual production font/layout measurements;
            # the installed non-browser tests retain their declared metric fixture.
            module.external_metrics_snapshot = measured
            getattr(module, entry)(root / directory)
    finally:
        browser.close()
