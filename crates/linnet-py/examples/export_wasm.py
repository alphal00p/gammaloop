#!/usr/bin/env python3
"""Lint, export, and smoke-test the editable Linnet Marimo notebooks."""

from __future__ import annotations

import argparse
import contextlib
import filecmp
import functools
import http.server
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading
import urllib.parse
import urllib.request
import zipfile
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
PUBLISHED_REQUIREMENT = '#     "linnet-py==0.1.0",'


@dataclass(frozen=True)
class Notebook:
    filename: str
    ready_value: str

    @property
    def source(self) -> Path:
        return EXAMPLES_DIR / self.filename

    @property
    def output_name(self) -> str:
        return f"{Path(self.filename).stem}.html"

    @property
    def ready_selector(self) -> str:
        return f'[data-linnet-render-ready="{self.ready_value}"] svg'


NOTEBOOKS = (
    Notebook("rendering_api.py", "generic"),
    Notebook("physics_render_settings.py", "physics"),
)


def validate_wasm_wheel(wheel: Path) -> Path:
    wheel = wheel.expanduser().resolve()
    filename = wheel.name.lower()
    if not wheel.is_file() or wheel.suffix != ".whl":
        raise ValueError(f"Linnet WASM wheel does not exist: {wheel}")
    if not filename.startswith("linnet_py-"):
        raise ValueError(f"Expected a linnet-py wheel, got {wheel.name}")
    if re.fullmatch(r"[a-z0-9_.+\-]+\.whl", filename) is None:
        raise ValueError(f"Wheel filename is not safe for PEP 723: {wheel.name}")

    platform = filename.removesuffix(".whl").rsplit("-", 1)[-1]
    if "emscripten" not in platform and "wasm32" not in platform:
        raise ValueError(
            "The local wheel must target Emscripten/wasm32, not the host "
            f"platform: {wheel.name}"
        )

    try:
        with zipfile.ZipFile(wheel) as archive:
            metadata = [
                name for name in archive.namelist() if name.endswith(".dist-info/WHEEL")
            ]
            if len(metadata) != 1:
                raise ValueError(
                    f"Expected one .dist-info/WHEEL in {wheel.name}, "
                    f"found {len(metadata)}"
                )
            wheel_metadata = archive.read(metadata[0]).decode("utf-8")
    except (OSError, UnicodeError, zipfile.BadZipFile) as error:
        raise ValueError(f"Invalid wheel archive: {wheel}") from error

    tags = [
        line.partition(":")[2].strip().lower()
        for line in wheel_metadata.splitlines()
        if line.lower().startswith("tag:")
    ]
    if not any("emscripten" in tag or "wasm32" in tag for tag in tags):
        raise ValueError(f"Wheel metadata has no Emscripten/wasm32 tag: {wheel}")
    return wheel


def with_local_wheel(source: str, wheel_name: str) -> str:
    """Replace the published dependency in a temporary notebook copy."""

    if source.count(PUBLISHED_REQUIREMENT) != 1:
        raise ValueError(
            "Expected exactly one pinned linnet-py dependency in notebook "
            "PEP 723 metadata"
        )
    local_requirement = f'#     "linnet-py @ ./wheels/{wheel_name}",'
    return source.replace(PUBLISHED_REQUIREMENT, local_requirement)


@contextlib.contextmanager
def staged_notebooks(
    wheel: Path | None,
) -> Iterator[tuple[tuple[Notebook, Path], ...]]:
    """Stage copies so a local wheel override never edits the notebooks."""

    with tempfile.TemporaryDirectory(prefix="linnet-marimo-wasm-") as temporary:
        stage = Path(temporary)
        if wheel is not None:
            wheels = stage / "wheels"
            wheels.mkdir()
            shutil.copy2(wheel, wheels / wheel.name)

        staged = []
        for notebook in NOTEBOOKS:
            source = notebook.source.read_text(encoding="utf-8")
            if wheel is not None:
                source = with_local_wheel(source, wheel.name)
            path = stage / notebook.filename
            path.write_text(source, encoding="utf-8")
            staged.append((notebook, path))
        yield tuple(staged)


def marimo(*arguments: str) -> None:
    command = (sys.executable, "-m", "marimo", *arguments)
    print(f"+ {shlex.join(command)}", flush=True)
    subprocess.run(command, check=True)


def lint(staged: Sequence[tuple[Notebook, Path]]) -> None:
    marimo(
        "check",
        "--strict",
        "--format",
        "full",
        "--select",
        # MW003 queries remote package indexes. The local Emscripten wheel is
        # validated separately, so keep export validation deterministic.
        "MW001,MW002",
        *(str(path) for _, path in staged),
    )


def export(
    staged: Sequence[tuple[Notebook, Path]], output: Path
) -> tuple[tuple[Notebook, Path], ...]:
    output.mkdir(parents=True, exist_ok=True)
    artifacts = []
    for notebook, source in staged:
        artifact = output / notebook.output_name
        marimo(
            "export",
            "html-wasm",
            str(source),
            "--output",
            str(artifact),
            "--mode",
            "edit",
            "--no-sandbox",
            "--no-execute",
            "--force",
        )
        artifacts.append((notebook, artifact))
    return tuple(artifacts)


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format: str, *args: object) -> None:
        pass


@contextlib.contextmanager
def serve(directory: Path) -> Iterator[str]:
    handler = functools.partial(QuietHandler, directory=str(directory))
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        host, port = server.server_address
        yield f"http://{host}:{port}/"
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def http_smoke(
    base_url: str,
    artifacts: Sequence[tuple[Notebook, Path]],
    output: Path,
    wheel: Path | None,
) -> None:
    assets = output / "assets"
    if not assets.is_dir() or not any(path.is_file() for path in assets.rglob("*")):
        raise RuntimeError(f"Marimo export did not produce assets under {assets}")

    if wheel is not None:
        hosted_wheel = output / "public" / "wheels" / wheel.name
        if not hosted_wheel.is_file():
            raise RuntimeError(f"Marimo did not host the local wheel at {hosted_wheel}")
        if not filecmp.cmp(wheel, hosted_wheel, shallow=False):
            raise RuntimeError(f"Hosted wheel differs from {wheel}")

    for notebook, artifact in artifacts:
        url = urllib.parse.urljoin(base_url, urllib.parse.quote(artifact.name))
        with urllib.request.urlopen(url, timeout=10) as response:
            body = response.read().decode("utf-8")
            if response.status != 200:
                raise RuntimeError(f"Invalid Marimo WASM response from {url}")
            required = {
                "embedded notebook source": "<marimo-code",
                "edit-mode configuration": '"mode": "edit"',
                "render readiness marker": (
                    f'data-linnet-render-ready=\\"{notebook.ready_value}\\"'
                ),
            }
            missing = [
                label for label, marker in required.items() if marker not in body
            ]
            if missing:
                raise RuntimeError(f"{artifact.name} is missing " + ", ".join(missing))
            if wheel is not None and wheel.name not in body:
                raise RuntimeError(
                    f"{artifact.name} does not reference the staged wheel"
                )
            if 'language=\\"typst\\"' in body:
                raise RuntimeError(
                    f"{artifact.name} requests Marimo's unsupported Typst highlighter"
                )


def browser_smoke(
    base_url: str,
    artifacts: Sequence[tuple[Notebook, Path]],
    timeout_seconds: float,
) -> None:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as error:
        raise RuntimeError(
            "--browser-smoke requires Playwright; install it and Chromium "
            "before running this check"
        ) from error

    timeout = timeout_seconds * 1000
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            for notebook, artifact in artifacts:
                page = browser.new_page()
                errors: list[str] = []
                page.on(
                    "pageerror",
                    lambda error, errors=errors: errors.append(str(error)),
                )
                url = urllib.parse.urljoin(base_url, urllib.parse.quote(artifact.name))
                response = page.goto(
                    url,
                    wait_until="domcontentloaded",
                    timeout=timeout,
                )
                if response is None or not response.ok:
                    raise RuntimeError(f"Browser failed to load {url}")
                page.locator(".cm-editor").first.wait_for(
                    state="visible",
                    timeout=timeout,
                )
                # Editable exports intentionally start with unexecuted cells.
                page.locator('[data-testid="run-button"]').last.click()
                page.locator(notebook.ready_selector).wait_for(
                    state="visible",
                    timeout=timeout,
                )
                if errors:
                    raise RuntimeError(
                        f"Browser errors while loading {artifact.name}: "
                        + "; ".join(errors)
                    )
                page.close()
        finally:
            browser.close()


def parse_args(arguments: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Directory for the two editable HTML-WASM notebooks",
    )
    parser.add_argument(
        "--wheel",
        type=Path,
        help=(
            "Local linnet-py Emscripten wheel. Without this option the "
            "published linnet-py==0.1.0 dependency is used."
        ),
    )
    parser.add_argument(
        "--browser-smoke",
        action="store_true",
        help="Open both exports in headless Chromium and wait for their SVGs",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=180,
        help="Per-page browser timeout in seconds (default: 180)",
    )
    return parser.parse_args(arguments)


def main(arguments: Sequence[str] | None = None) -> int:
    options = parse_args(arguments)
    if options.timeout <= 0:
        raise ValueError("--timeout must be greater than zero")
    wheel = validate_wasm_wheel(options.wheel) if options.wheel else None
    output = options.output.expanduser().resolve()

    with staged_notebooks(wheel) as staged:
        lint(staged)
        artifacts = export(staged, output)

    with serve(output) as base_url:
        http_smoke(base_url, artifacts, output, wheel)
        if options.browser_smoke:
            browser_smoke(base_url, artifacts, options.timeout)

    print(f"Editable Marimo WASM notebooks exported to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
