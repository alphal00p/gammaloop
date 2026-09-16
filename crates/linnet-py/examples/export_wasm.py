#!/usr/bin/env python3
"""Lint, export, and smoke-test the editable documentation Marimo notebooks."""

from __future__ import annotations

import argparse
import base64
import contextlib
import filecmp
import functools
import http.server
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import threading
import tomllib
import urllib.parse
import urllib.request
import uuid
import zipfile
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

EXAMPLES_DIR = Path(__file__).resolve().parent
PUBLISHED_REQUIREMENTS = {
    "linnet-py": '#     "linnet-py==0.1.0",',
    "symbolica": '#     "symbolica==3.0.0",',
    "ufo-model-loader": '#     "ufo-model-loader @ git+https://github.com/alphal00p/ufo_model_loader.git@70ddee6b416f8c8b340e0d087646d77095c5d24b",',
}


@dataclass(frozen=True)
class Notebook:
    filename: str
    ready_value: str
    docs_product: str
    docs_route: str
    package: str = "linnet-py"

    @property
    def source(self) -> Path:
        return EXAMPLES_DIR / self.filename

    @property
    def output_name(self) -> str:
        return f"{Path(self.filename).stem}.html"

    @property
    def ready_selector(self) -> str:
        if self.package == "symbolica":
            return f'[data-notebook-ready="{Path(self.filename).stem}"]'
        if self.ready_value == "quickstart":
            return '[data-notebook="python_quickstart"] svg[width$="pt"]'
        return f'[data-linnet-render-ready="{self.ready_value}"] svg'


NOTEBOOKS = (
    Notebook("rendering_api.py", "generic", "linnet", "guides/python-rendering/"),
    Notebook("physics_render_settings.py", "physics", "gammaloop", "guides/dot-input/"),
    Notebook("layout_stream.py", "stream", "linnet", "playground/"),
    Notebook(
        "../../../examples/notebooks/spenso_idenso_display.py",
        "spenso_idenso_display",
        "spenso",
        "guides/showcase/",
        "symbolica",
    ),
    *(
        Notebook(
            f"../../../examples/notebooks/feynkit/{filename}.py",
            filename,
            "feynkit",
            f"guides/showcases/{route}/",
            "symbolica",
        )
        for filename, route in (
            ("00_quickstart_marimo", "first-diagram"),
            ("01_models_and_diagrams_marimo", "models-and-diagrams"),
            ("02_cff_and_symbolica_marimo", "cff"),
            ("03_kinematics_and_jets_marimo", "kinematics"),
            ("04_ufo_loading_marimo", "ufo"),
            ("07_tensor_reduction_marimo", "tensor-reduction"),
        )
    ),
)


def validate_notebook_wheel(wheel: Path, package: str = "linnet-py") -> Path:
    wheel = wheel.expanduser().resolve()
    filename = wheel.name.lower()
    if not wheel.is_file() or wheel.suffix != ".whl":
        raise ValueError(f"WASM wheel does not exist: {wheel}")
    if not filename.startswith(package.replace("-", "_") + "-"):
        raise ValueError(f"Expected a {package} wheel, got {wheel.name}")
    if re.fullmatch(r"[a-z0-9_.+\-]+\.whl", filename) is None:
        raise ValueError(f"Wheel filename is not safe for PEP 723: {wheel.name}")

    platform = filename.removesuffix(".whl").rsplit("-", 1)[-1]
    if package == "ufo-model-loader":
        if not filename.endswith("-py3-none-any.whl"):
            raise ValueError("The UFO loader must be a pure Python wheel")
    elif "emscripten" not in platform and "wasm32" not in platform:
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
    if package == "ufo-model-loader":
        if "py3-none-any" not in tags:
            raise ValueError("The UFO loader metadata must declare py3-none-any")
    elif not any("emscripten" in tag or "wasm32" in tag for tag in tags):
        raise ValueError(f"Wheel metadata has no Emscripten/wasm32 tag: {wheel}")
    return wheel


def with_local_wheel(source: str, wheel_name: str, package: str = "linnet-py") -> str:
    """Replace the published dependency in a temporary notebook copy."""

    published_requirement = PUBLISHED_REQUIREMENTS[package]
    if source.count(published_requirement) != 1:
        raise ValueError(
            f"Expected exactly one pinned {package} dependency in notebook "
            "PEP 723 metadata"
        )
    local_requirement = f'#     "{package} @ ./wheels/{wheel_name}",'
    return source.replace(published_requirement, local_requirement)


@contextlib.contextmanager
def staged_notebooks(
    wheel: Path | None,
    notebooks: Sequence[Notebook],
    gammaloop: Path | None = None,
    dependency_wheel: Path | None = None,
) -> Iterator[tuple[tuple[Notebook, Path], ...]]:
    """Stage copies so a local wheel override never edits the notebooks."""

    with tempfile.TemporaryDirectory(prefix="docs-marimo-wasm-") as temporary:
        stage = Path(temporary)
        if wheel is not None:
            wheels = stage / "wheels"
            wheels.mkdir()
            shutil.copy2(wheel, wheels / wheel.name)
            if dependency_wheel is not None:
                shutil.copy2(dependency_wheel, wheels / dependency_wheel.name)

        staged = []
        for notebook in notebooks:
            source = notebook.source.read_text(encoding="utf-8")
            if wheel is not None:
                source = with_local_wheel(source, wheel.name, notebook.package)
            if dependency_wheel is not None:
                dependency = dependency_wheel.name.split("-")[0].replace("_", "-")
                if PUBLISHED_REQUIREMENTS[dependency] in source:
                    source = with_local_wheel(source, dependency_wheel.name, dependency)
            if notebook.ready_value == "physics":
                repository = EXAMPLES_DIR.parents[2]
                if gammaloop is None:
                    package = subprocess.check_output(
                        [
                            "nix",
                            "build",
                            "--no-link",
                            "--print-out-paths",
                            ".#gammaloop",
                        ],
                        cwd=repository,
                        text=True,
                    ).strip()
                    gammaloop = Path(package) / "bin/gammaloop"
                gammaloop = gammaloop.expanduser().resolve()
                drawing_export = stage / "gammaloop-export"
                model = repository / "assets/models/json/sm/sm.json"
                print(
                    "Generating GammaLoop drawing bundle from the Standard Model",
                    flush=True,
                )
                subprocess.run(
                    [
                        str(gammaloop),
                        "-l",
                        "warn",
                        "--state-folder",
                        str(stage / "gammaloop-state"),
                        "--no-save-state",
                        "run",
                        "-c",
                        (
                            f"import model {shlex.quote(str(model))}; "
                            f"save dot {shlex.quote(str(drawing_export))}"
                        ),
                    ],
                    cwd=repository,
                    check=True,
                )
                templates = drawing_export / "drawings/templates"
                # GammaLoop's CLI can report execution errors without a failing
                # process status. Require a complete bundle from this source tree.
                canonical = repository / "assets/embedded/drawing/templates"
                sources = [
                    (expected, templates / expected.relative_to(canonical))
                    for expected in canonical.rglob("*.typ")
                ]
                for package in ("linnest", "kurvst"):
                    relative = Path("crates") / package / "typst"
                    canonical = repository / relative
                    sources.extend(
                        (
                            expected,
                            templates / relative / expected.relative_to(canonical),
                        )
                        for expected in (
                            canonical / "typst.toml",
                            *(canonical / "src").rglob("*.typ"),
                        )
                    )
                    wasm = templates / relative / f"{package}.wasm"
                    if not wasm.is_file():
                        raise RuntimeError(
                            f"GammaLoop drawing bundle is missing {wasm}"
                        )
                for expected, actual in sources:
                    if (
                        not actual.is_file()
                        or actual.read_bytes() != expected.read_bytes()
                    ):
                        raise RuntimeError(
                            "GammaLoop exported missing or stale template "
                            f"{actual.relative_to(templates)}"
                        )
                particle_map = templates / "edge-style.typ"
                if not particle_map.is_file() or not all(
                    marker in particle_map.read_text(encoding="utf-8")
                    for marker in (
                        "#let generated-map = (",
                        '"a":',
                        '"g":',
                        '"t":',
                        '"H":',
                    )
                ):
                    raise RuntimeError(
                        "GammaLoop did not generate the Standard Model particle map"
                    )
                package_cache = os.environ.get("TYPST_PACKAGE_CACHE_PATH")
                if not package_cache:
                    raise RuntimeError(
                        "Use the flake to provide TYPST_PACKAGE_CACHE_PATH for MiTeX"
                    )
                roots = [(templates, "drawings/templates")]
                for package in ("cetz/0.5.1", "oxifmt/1.0.0", "mitex/0.2.6"):
                    package_root = (
                        (
                            Path(package_cache)
                            if package.startswith("mitex/")
                            else EXAMPLES_DIR.parent / "vendor/typst-packages"
                        )
                        / "preview"
                        / package
                    )
                    if not (package_root / "typst.toml").is_file():
                        raise RuntimeError(
                            f"Missing pinned Typst package {package_root}"
                        )
                    roots.append((package_root, "typst-packages/preview/" + package))
                assets = sorted(
                    (f"{prefix}/{asset.relative_to(root).as_posix()}", asset)
                    for root, prefix in roots
                    for asset in root.rglob("*")
                    if asset.is_file()
                )
                bundle = stage / "gammaloop-drawing.zip"
                with zipfile.ZipFile(bundle, "w") as archive:
                    for name, asset in assets:
                        entry = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                        entry.compress_type = zipfile.ZIP_DEFLATED
                        entry.create_system = 3
                        entry.external_attr = 0o644 << 16
                        archive.writestr(entry, asset.read_bytes())
                marker = "    drawing_bundle = None"
                if source.count(marker) != 1:
                    raise ValueError(
                        "Expected one drawing_bundle placeholder in the physics notebook"
                    )
                encoded = base64.b64encode(bundle.read_bytes()).decode("ascii")
                source = source.replace(marker, f"    drawing_bundle = {encoded!r}")
            if notebook.docs_product == "feynkit":
                source = source.replace(
                    "Path(__file__).resolve().parents[3]", 'Path("/feynkit-data")'
                )
            path = stage / Path(notebook.filename).name
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
    staged: Sequence[tuple[Notebook, Path]], output: Path, *, docs: str | None = None
) -> tuple[tuple[Notebook, Path], ...]:
    output.mkdir(parents=True, exist_ok=True)
    bundle = staged[0][1].parent / "gammaloop-drawing.zip"
    if bundle.is_file():
        (output / "public").mkdir(exist_ok=True)
        shutil.copy2(bundle, output / "public/gammaloop-drawing.zip")
    artifacts = []
    if docs:
        import marimo as mo
        from marimo._session.notebook.loader import load_notebook

        generators = []
        # Islands omit PEP 723 metadata. Install the staged wheel explicitly and
        # make the import cell depend on it; the browser supplies its absolute URL.
        for notebook, source in staged:
            text = source.read_text(encoding="utf-8")
            metadata = re.search(r"# /// script\n(.*?)# ///", text, re.DOTALL)
            requirements = tomllib.loads(re.sub(r"(?m)^# ?", "", metadata.group(1)))[
                "dependencies"
            ]
            requirements = [
                "__NOTEBOOK_WHEEL_URL__"
                if item.startswith(notebook.package + " @")
                else "__DEPENDENCY_WHEEL_URL__"
                if " @ ./wheels/" in item
                else item
                for item in requirements
                if not item.startswith("marimo==")
            ]
            source.write_text(
                re.sub(
                    r"(?m)^(    )(import |from )",
                    r"\1notebook_browser_ready\n\1\2",
                    text,
                ),
                encoding="utf-8",
            )
            generator = mo.MarimoIslandGenerator.from_file(
                str(source),
                display_code=notebook.ready_value == "generic"
                or notebook.package == "symbolica",
            )
            if notebook.package == "symbolica":
                # Marimo 0.24's from_file applies one display flag to every cell.
                # Preserve the notebook's prose/control cells and expose computation code.
                cells = load_notebook(str(source)).app.cell_manager.cell_data()
                for island, cell in zip(generator.stubs, cells, strict=True):
                    island._display_code = not cell.config.hide_code
            bootstrap = ""
            if notebook.docs_product == "feynkit":
                repository = EXAMPLES_DIR.parents[2]
                assets = [
                    repository
                    / "crates/feynkit-model/tests/fixtures/scalars_2p_3p.json",
                    repository / "crates/feynkit-model/tests/fixtures/sm.json",
                    *(
                        asset
                        for asset in (repository / "assets/models/ufo/scalars").rglob(
                            "*"
                        )
                        if asset.is_file() and "__pycache__" not in asset.parts
                    ),
                ]
                bundle = source.parent / "feynkit-data.zip"
                with zipfile.ZipFile(bundle, "w") as archive:
                    for asset in sorted(assets):
                        entry = zipfile.ZipInfo(
                            asset.relative_to(repository).as_posix(),
                            date_time=(1980, 1, 1, 0, 0, 0),
                        )
                        entry.compress_type = zipfile.ZIP_DEFLATED
                        archive.writestr(entry, asset.read_bytes())
                encoded = base64.b64encode(bundle.read_bytes()).decode("ascii")
                bootstrap = (
                    "import base64 as _base64, io as _io, zipfile as _zipfile\n"
                    f"with _zipfile.ZipFile(_io.BytesIO(_base64.b64decode({encoded!r}))) as _zip:\n"
                    "    _zip.extractall('/feynkit-data')\n"
                )
            generator.add_code(
                "import micropip as _micropip\n"
                f"await _micropip.install({requirements!r}, reinstall=True)\n"
                + bootstrap
                + "notebook_browser_ready = True",
            )
            generators.append((notebook, generator))

        if docs == "linnet":
            # The live quickstart uses the canonical documented example verbatim,
            # then renders the resulting graph. Keeping
            # output in the editable cell also gives it Marimo's rerun control.
            quickstart = (
                EXAMPLES_DIR.parents[2]
                / "docs/products/linnet/content/quickstart-python.typ"
            )
            match = re.search(
                r"// docs-example: compile linnet-python-quickstart\s*```python\n(.*?)\n```",
                quickstart.read_text(encoding="utf-8"),
                re.DOTALL,
            )
            if match is None:
                raise ValueError("The canonical Linnet Python quickstart was not found")
            generator = mo.MarimoIslandGenerator()
            generator.add_code(
                "import micropip as _micropip\n"
                "await _micropip.install('__NOTEBOOK_WHEEL_URL__')\n"
                "notebook_browser_ready = True",
            )
            generator.add_code(
                "notebook_browser_ready\n" + match.group(1) + "\n\n"
                "import marimo as mo\n"
                "mo.Html(graph.to_svg())",
                display_code=True,
            )
            generators.append(
                (
                    Notebook(
                        "python_quickstart.py",
                        "quickstart",
                        "linnet",
                        "quickstart/python/",
                    ),
                    generator,
                )
            )
        wheels = sorted((staged[0][1].parent / "wheels").glob("*.whl"))
        wheel = next(
            wheel
            for wheel in wheels
            if wheel.name.startswith(staged[0][0].package.replace("-", "_") + "-")
        )
        hosted_wheel = output / "public" / "wheels" / wheel.name
        hosted_wheel.parent.mkdir(parents=True, exist_ok=True)
        for local_wheel in wheels:
            shutil.copy2(local_wheel, hosted_wheel.parent / local_wheel.name)
        for notebook, generator in generators:
            body = generator.render_body(
                include_init_island=False, include_payload=True
            )
            # Static editor handles are random UUIDs in Marimo. Stable handles
            # keep repeated exports identical for immutable documentation snapshots.
            for index, identifier in enumerate(
                dict.fromkeys(re.findall(r"object-id='([^']+)'", body))
            ):
                stable = uuid.uuid5(uuid.NAMESPACE_URL, f"{notebook.filename}:{index}")
                body = body.replace(f"'{identifier}'", f"'{stable}'")
            artifact = (output / notebook.output_name).with_suffix(".json")
            artifact.write_text(
                json.dumps(
                    {
                        "head": generator.render_head(),
                        "body": body,
                        "wheel": "public/wheels/" + wheel.name,
                        "dependency_wheel": next(
                            (
                                "public/wheels/" + item.name
                                for item in wheels
                                if item != wheel
                            ),
                            None,
                        ),
                    }
                ),
                encoding="utf-8",
            )
            artifacts.append((notebook, artifact))
        return tuple(artifacts)
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
    *,
    docs: bool = False,
) -> None:
    assets = output / "assets"
    if not docs and (
        not assets.is_dir() or not any(path.is_file() for path in assets.rglob("*"))
    ):
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
            if docs:
                payload = json.loads(body)
                if (
                    "<marimo-island" not in payload["body"]
                    or "__NOTEBOOK_WHEEL_URL__" not in body
                ):
                    raise RuntimeError(
                        f"{artifact.name} is missing its executable islands"
                    )
                if payload["wheel"] != "public/wheels/" + wheel.name:
                    raise RuntimeError(f"{artifact.name} references the wrong wheel")
                continue
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
    browser_executable: Path | None,
    *,
    docs: bool = False,
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
        browser = playwright.chromium.launch(
            executable_path=str(browser_executable) if browser_executable else None
        )
        try:
            for notebook, artifact in artifacts:
                page = browser.new_page()
                errors: list[str] = []
                page.on(
                    "pageerror",
                    lambda error, errors=errors: errors.append(str(error)),
                )
                route = (
                    notebook.docs_route if docs else urllib.parse.quote(artifact.name)
                )
                url = urllib.parse.urljoin(base_url, route)
                response = page.goto(
                    url,
                    wait_until="domcontentloaded",
                    timeout=timeout,
                )
                if response is None or not response.ok:
                    raise RuntimeError(f"Browser failed to load {url}")
                if docs:
                    page.locator("[data-notebook]").scroll_into_view_if_needed()
                else:
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
                if docs:
                    page.wait_for_function(
                        """() => !document.querySelector(
                            'marimo-island[data-status="running"], marimo-island[data-status="queued"]'
                        )""",
                        timeout=timeout,
                    )
                    runtime_errors = page.locator(
                        "marimo-island .text-error"
                    ).all_text_contents()
                    if runtime_errors:
                        raise RuntimeError(
                            f"Notebook errors while loading {artifact.name}: "
                            + "; ".join(runtime_errors)
                        )
                if notebook.ready_value == "physics":
                    settings = page.get_by_role(
                        "button", name="Layout settings", exact=True
                    )
                    if settings.get_attribute("aria-expanded") != "false":
                        raise RuntimeError(
                            "Physics layout settings must start collapsed"
                        )
                    settings.click()
                    page.get_by_role("slider").first.wait_for(
                        state="visible", timeout=timeout
                    )
                    for name in (
                        "Momentum arrows",
                        "Momentum labels qₑ",
                        "Cross-section",
                    ):
                        previous = page.locator(notebook.ready_selector).evaluate(
                            "svg => svg.outerHTML"
                        )
                        if name == "Cross-section":
                            page.get_by_role(
                                "combobox", name="Example", exact=True
                            ).select_option(label=name)
                        else:
                            page.get_by_role("checkbox", name=name, exact=True).check()
                        page.wait_for_function(
                            """([selector, previous]) => {
                                const svg = document.querySelector(selector);
                                return svg && svg.outerHTML !== previous;
                            }""",
                            arg=[notebook.ready_selector, previous],
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
        help="Directory for the editable HTML-WASM notebooks",
    )
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--docs",
        choices=("linnet", "gammaloop", "spenso", "idenso", "feynkit"),
        help="Export this product's live cells into its built assets/notebooks directory (requires --wheel)",
    )
    selection.add_argument(
        "--notebook",
        choices=[Path(notebook.filename).stem for notebook in NOTEBOOKS],
        help="Export one notebook (default: all)",
    )
    parser.add_argument(
        "--wheel",
        type=Path,
        help=(
            "Local Emscripten wheel: linnet-py for Linnet/GammaLoop, or symbolica "
            "with FeynKit/Spenso/Idenso for community showcases. Without this option the "
            "published linnet-py==0.1.0 dependency is used."
        ),
    )
    parser.add_argument(
        "--dependency-wheel",
        type=Path,
        help="Linnet WASM wheel for Spenso/Idenso, or the pinned UFO loader wheel for FeynKit",
    )
    parser.add_argument(
        "--gammaloop",
        type=Path,
        help="GammaLoop executable for the physics bundle (default: build this revision with Nix)",
    )
    parser.add_argument(
        "--browser-smoke",
        action="store_true",
        help="Open the exports in headless Chromium and wait for their SVGs",
    )
    parser.add_argument(
        "--browser-executable",
        type=Path,
        help="Chromium executable for --browser-smoke (for example from Nix)",
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
    output = options.output.expanduser().resolve()
    if options.docs and options.wheel is None:
        raise ValueError("--docs requires --wheel")
    if options.docs and (output.parent.name, output.name) != ("assets", "notebooks"):
        raise ValueError(
            "Build the product docs first and use --output <product-version>/assets/notebooks"
        )
    notebooks = tuple(
        notebook
        for notebook in NOTEBOOKS
        if (
            (
                notebook.docs_product == options.docs
                or options.docs == "idenso"
                and notebook.docs_product == "spenso"
            )
            if options.docs
            else (options.notebook is None and notebook.package == "linnet-py")
            or Path(notebook.filename).stem == options.notebook
        )
    )
    packages = {notebook.package for notebook in notebooks}
    if len(packages) != 1:
        raise ValueError(
            "Select --docs or --notebook so one host wheel serves the export"
        )
    package = packages.pop()
    if package == "symbolica" and not options.docs:
        raise ValueError(
            "Community showcases require --docs and a combined Symbolica WASM wheel"
        )
    wheel = validate_notebook_wheel(options.wheel, package) if options.wheel else None
    dependency_wheel = (
        validate_notebook_wheel(
            options.dependency_wheel,
            "ufo-model-loader" if options.docs == "feynkit" else "linnet-py",
        )
        if options.dependency_wheel
        else None
    )
    if options.docs in {"spenso", "idenso", "feynkit"} and dependency_wheel is None:
        raise ValueError(
            "Community showcases require --dependency-wheel (Linnet or UFO loader)"
        )
    if dependency_wheel is not None and options.docs not in {
        "spenso",
        "idenso",
        "feynkit",
    }:
        raise ValueError("--dependency-wheel applies only to community showcases")
    with staged_notebooks(
        wheel, notebooks, options.gammaloop, dependency_wheel
    ) as staged:
        lint(staged)
        artifacts = export(staged, output, docs=options.docs)

    with serve(output.parent.parent if options.docs else output) as base_url:
        http_smoke(
            urllib.parse.urljoin(base_url, "assets/notebooks/")
            if options.docs
            else base_url,
            artifacts,
            output,
            wheel,
            docs=bool(options.docs),
        )
        if options.browser_smoke:
            browser_smoke(
                base_url,
                artifacts,
                options.timeout,
                options.browser_executable,
                docs=bool(options.docs),
            )

    print(f"Marimo WASM notebooks exported to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
