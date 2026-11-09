#!/usr/bin/env python3
"""Lint, export, and smoke-test the editable documentation Marimo notebooks."""

# ruff: noqa: EXE001  # Publication invokes this source through Python.

from __future__ import annotations

import argparse
import ast
import base64
import contextlib
import filecmp
import functools
import hashlib
import html
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
import time
import tomllib
import urllib.parse
import urllib.request
import uuid
import zipfile
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

EXECUTION_TRACKING = r"""(() => {
                        window.notebookExecutedCells = new Set();
                        window.notebookRunCompleted = false;
                        window.notebookCompletedRuns = 0;
                        window.notebookCellStarts = new Map();
                        window.notebookCellTimings = new Map();
                        const NotebookWorker = window.Worker;
                        window.Worker = class extends NotebookWorker {
                            constructor(...args) {
                                super(...args);
                                this.addEventListener("message", ({data}) => {
                                    const rpc = data.data;
                                    if (rpc?.id !== "kernelMessage") return;
                                    const message = JSON.parse(rpc.payload.message);
                                    if (message.op === "completed-run") {
                                        window.notebookRunCompleted = true;
                                        window.notebookCompletedRuns += 1;
                                    } else if (message.op === "cell-op") {
                                        const {cell_id, status, timestamp} = message.data;
                                        if (status === "queued") {
                                            window.notebookCellStarts.delete(cell_id);
                                            window.notebookCellTimings.delete(cell_id);
                                        } else if (status === "running") {
                                            window.notebookCellStarts.set(cell_id, timestamp);
                                        } else if (status === "idle" && window.notebookCellStarts.has(cell_id)) {
                                            window.notebookCellTimings.set(cell_id,
                                                (timestamp - window.notebookCellStarts.get(cell_id)) * 1000);
                                            window.notebookCellStarts.delete(cell_id);
                                        }
                                        if (message.data.status === "idle") {
                                            window.notebookExecutedCells.add(message.data.cell_id);
                                        } else if (["queued", "running"].includes(message.data.status)) {
                                            window.notebookRunCompleted = false;
                                        }
                                    }
                                });
                            }
                        };
                    })()"""

EXAMPLES_DIR = Path(__file__).resolve().parent
PUBLISHED_REQUIREMENTS = {
    "symbolica": '#     "symbolica==3.0.1",',
    "ufo-model-loader": '#     "ufo-model-loader==1.0.1",',
}


@dataclass(frozen=True)
class Notebook:
    filename: str
    ready_value: str
    docs_product: str
    docs_route: str
    package: str = "symbolica"

    @property
    def source(self) -> Path:
        return EXAMPLES_DIR / self.filename

    @property
    def output_name(self) -> str:
        return f"{Path(self.filename).stem}.html"

    @property
    def ready_selector(self) -> str:
        if self.ready_value in {"generic", "physics"}:
            return 'iframe[title="Graph diagram"]'
        if self.docs_product != "linnet":
            return f'[data-notebook-ready="{Path(self.filename).stem}"]'
        return f'[data-linnet-render-ready="{self.ready_value}"] svg'


NOTEBOOKS = (
    Notebook("rendering_api.py", "generic", "linnet", "guides/python-rendering/"),
    Notebook(
        "physics_render_settings.py",
        "physics",
        "gammaloop",
        "guides/dot-input/",
        "symbolica",
    ),
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
            f"../../../examples/notebooks/{name}.py",
            name,
            "spenso",
            f"guides/{route}/",
        )
        for name, route in (
            ("spenso_api_tour", "api-tour"),
            ("spenso_current_kernel", "current-kernel"),
            ("spenso_notation", "notation"),
            ("spenso_tensor_explorer", "tensor-explorer"),
        )
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
            ("05_four_loop_numerator_marimo", "four-loop-numerator"),
            ("07_tensor_reduction_marimo", "tensor-reduction"),
            ("08_numerical_cff_ltd_marimo", "numerical-integration"),
        )
    ),
)


def validate_notebook_wheel(wheel: Path, package: str = "symbolica") -> Path:
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


def hosted_wheel_path(wheel: Path) -> Path:
    """Keep previously published wheels immutable when their ABI name repeats."""
    identity = hashlib.sha256(wheel.read_bytes()).hexdigest()
    return Path("public/wheels") / identity / wheel.name


def with_local_wheel(source: str, wheel_name: str, package: str = "symbolica") -> str:
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
            if notebook.ready_value == "04_ufo_loading_marimo":
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


def browser_install_code(requirements: Sequence[str]) -> str:
    """Reuse installed wheel URLs and replace conflicting explicit package pins."""
    return (
        "import micropip as _micropip\n"
        "from packaging.requirements import Requirement as _Requirement\n"
        "from packaging.utils import canonicalize_name as _canonicalize_name\n"
        f"_requirements = {requirements!r}\n"
        "_installed = {_canonicalize_name(name): package for name, package in _micropip.list().items()}\n"
        "for _requirement in _requirements:\n"
        "    if _requirement.startswith(('https://', 'http://')):\n"
        "        continue\n"
        "    _constraint = _Requirement(_requirement)\n"
        "    _package = _installed.get(_canonicalize_name(_constraint.name))\n"
        "    if _package and not _constraint.specifier.contains(_package.version):\n"
        "        _micropip.uninstall(_constraint.name)\n"
        "_installed_urls = {package.source for package in _installed.values()}\n"
        "_requirements = [requirement for requirement in _requirements if requirement not in _installed_urls]\n"
        "await _micropip.install(_requirements)\n"
    )


def documentation_site(output: Path) -> Path:
    for parent in output.parents:
        if parent.name == "products":
            return parent.parent
    raise ValueError("Documentation exports require a built site's products directory")


def export(
    staged: Sequence[tuple[Notebook, Path]],
    output: Path,
    *,
    docs: str | None = None,
    runtime: Path | None = None,
) -> tuple[tuple[Notebook, Path], ...]:
    output.mkdir(parents=True, exist_ok=True)
    artifacts = []
    if docs:
        import marimo as mo
        from marimo._session.notebook.loader import load_notebook

        generators = []
        setup_cells = {}
        prose_cells = {}
        hidden_source_cells = {}
        full_sources = {}
        # Islands omit PEP 723 metadata. Install the staged wheel in the first
        # cell; the browser supplies its absolute URL.
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
            # Export the full upstream editor as an explicit alternative to islands.
            # It supplies notebook-level controls without changing the guide layout.
            generator = mo.MarimoIslandGenerator()
            bootstrap = ""
            if notebook.ready_value == "04_ufo_loading_marimo":
                repository = EXAMPLES_DIR.parents[2]
                ufo_model = repository / "assets/models/ufo/scalars"
                # Require the model assets from this source tree; missing fixtures
                # must fail export rather than leave the browser with a partial bundle.
                if not (ufo_model / "__init__.py").is_file():
                    raise FileNotFoundError(ufo_model / "__init__.py")
                assets = [
                    asset
                    for asset in ufo_model.rglob("*")
                    if asset.is_file() and "__pycache__" not in asset.parts
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
            if notebook.ready_value == "spenso_api_tour":
                stub = EXAMPLES_DIR.parents[2] / "docs/api/python/spynso3.pyi"
                shutil.copyfile(stub, output / "spynso3.pyi")
            full_source = source
            if notebook.ready_value == "spenso_api_tour":
                # The standalone editor cannot use the guide's URL substitution.
                # Embed its reviewed catalog so the full export remains portable.
                full_source = source.with_name("full_" + source.name)
                full_source.write_text(
                    text.replace(
                        'source = open_url("__SPENSO_STUB_URL__").read()',
                        "source = " + repr(stub.read_text()),
                    ),
                    encoding="utf-8",
                )
            if bootstrap:
                full_source = source.with_name("full_" + source.name)
                setup = (
                    "@app.cell(hide_code=True)\ndef _():\n"
                    + "".join("    " + line + "\n" for line in bootstrap.splitlines())
                    + "    return\n\n\n"
                )
                full_source.write_text(
                    text.replace("@app.cell", setup + "@app.cell", 1), encoding="utf-8"
                )
            # Full exports and islands have different injected setup cells.
            # Retain the actual full program's native cell identities for ty.
            full_sources[notebook.filename] = load_notebook(
                str(full_source)
            ).app.cell_manager.code_map()
            full = output / "full" / notebook.output_name
            marimo(
                "export",
                "html-wasm",
                str(full_source),
                "--output",
                str(full),
                "--mode",
                "edit",
                "--no-sandbox",
                "--no-execute",
                "--force",
            )
            # Native Marimo already supplies an override completion source.
            # Let the appended ty source own completion through the existing
            # CodeMirror facet. The static export has no native language-server
            # connection; retaining that pending source would delay every result.
            # This changes completion configuration, never editor state/history.
            patches = (
                (
                    "dist-Ccyblc8N.js",
                    "{defaultKeymap:(e,t)=>e&&t,closeOnBlur:",
                    "{override:(e,t)=>t,defaultKeymap:(e,t)=>e&&t,closeOnBlur:",
                    1,
                ),
                # Marimo's Escape handler can run without its optional AI
                # suggestion field. A non-throwing lookup preserves normal
                # Escape handling for ty help and native notebook controls.
                (
                    "RunButton-C2kFUhyO.js",
                    "let t=e.state.field(vc)?.suggestion",
                    "let t=e.state.field(vc,!1)?.suggestion",
                    3,
                ),
            )
            for filename, original, extended, expected in patches:
                asset = output / "full/assets" / filename
                content = asset.read_text()
                if (
                    content.count(original) != expected
                    and content.count(extended) != expected
                ):
                    raise ValueError(
                        f"Pinned full-notebook editor asset changed: {filename}"
                    )
                asset.write_text(content.replace(original, extended))
            # The full app retains its own editor; the website adds the same
            # type service using that app's pinned CodeMirror constructors.
            full.write_text(
                full.read_text().replace(
                    "</body>",
                    r"""<script type="module" data-documentation-editor>
const configURL = new URL('../../documentation-cells-config.json', location.href);
const config = await fetch(configURL).then(response => response.json());
const entry = new URL(config.entry, configURL);
for (const name of ['fonts.css', 'code.css', 'documentation-full.css']) {
  const link = document.createElement('link'); link.rel = 'stylesheet';
  link.href = new URL('../../' + name, location.href); document.head.append(link);
}
const module = await import(entry.href);
await module.initializeFullNotebook(configURL, new URL('./assets/', location.href),
  new URL('../' + location.pathname.split('/').pop().replace(/\.html$/, '.json'), location.href));
</script></body>""",
                ),
                encoding="utf-8",
            )
            # Keep compatible packages on reruns. Only replace preloaded packages
            # that conflict with the source's explicit requirements (for example
            # layout_stream's pinned AnyWidget); never force-reinstall the wheel.
            generator.add_code(
                browser_install_code(requirements) + bootstrap, display_code=True
            )
            # Marimo executes independent cells in notebook order. Install first,
            # then preserve source verbatim: no bootstrap markers in editable code.
            setup_cells[notebook.filename] = [0]
            prose_cells[notebook.filename] = []
            hidden_source_cells[notebook.filename] = []
            marimo_aliases = {
                alias.asname or alias.name
                for node in ast.walk(ast.parse(text))
                if isinstance(node, ast.Import)
                for alias in node.names
                if alias.name == "marimo"
            }
            for index, cell in enumerate(
                load_notebook(str(source)).app.cell_manager.cell_data(), 1
            ):
                statements = ast.parse(cell.code).body
                if (
                    cell.config.hide_code
                    and statements
                    and all(
                        isinstance(node, (ast.Import, ast.ImportFrom))
                        for node in statements
                    )
                ):
                    setup_cells[notebook.filename].append(index)
                # Standalone Markdown expressions are prose, including dynamic
                # arguments. Keep their dependency execution in the real kernel.
                prose = (
                    len(statements) == 1
                    and isinstance(statements[0], ast.Expr)
                    and isinstance(statements[0].value, ast.Call)
                    and isinstance(statements[0].value.func, ast.Attribute)
                    and statements[0].value.func.attr == "md"
                    and isinstance(statements[0].value.func.value, ast.Name)
                    and statements[0].value.func.value.id in marimo_aliases
                )
                if prose:
                    prose_cells[notebook.filename].append(index)
                if cell.config.hide_code and not prose:
                    hidden_source_cells[notebook.filename].append(index)
                generator.add_code(cell.code, display_code=not prose)
            generators.append((notebook, generator))

        wheels = sorted((staged[0][1].parent / "wheels").glob("*.whl"))
        wheel = next(
            wheel
            for wheel in wheels
            if wheel.name.startswith(staged[0][0].package.replace("-", "_") + "-")
        )
        for local_wheel in wheels:
            destination = output / hosted_wheel_path(local_wheel)
            destination.parent.mkdir(parents=True, exist_ok=True)
            if destination.exists():
                if not filecmp.cmp(local_wheel, destination, shallow=False):
                    raise ValueError(f"Immutable hosted wheel collision: {destination}")
            else:
                shutil.copyfile(local_wheel, destination)
        # Full Marimo contains two pinned editor fixes above. Address its whole
        # module tree so a browser cannot combine cached original chunks with
        # the new completion adapter. All relative module imports remain intact.
        assets = output / "full/assets"
        files = {
            path.relative_to(assets).as_posix(): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in sorted(assets.rglob("*"))
            if path.is_file()
        }
        identity = hashlib.sha256(
            json.dumps(files, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        destination = output / "full/runtime" / identity / "assets"
        if destination.exists():
            existing = {
                path.relative_to(destination).as_posix(): hashlib.sha256(
                    path.read_bytes()
                ).hexdigest()
                for path in sorted(destination.rglob("*"))
                if path.is_file()
            }
            if files != existing:
                raise ValueError(
                    f"Immutable full-notebook runtime collision: {identity}"
                )
        else:
            shutil.copytree(assets, destination)
            (destination.parent / "manifest.json").write_text(
                json.dumps({"identity": identity, "files": files})
            )
        shutil.rmtree(assets)
        for notebook, _ in generators:
            full = output / "full" / notebook.output_name
            html = full.read_text().replace(
                "./assets/", f"./runtime/{identity}/assets/"
            )
            for local_wheel in wheels:
                original = "../public/wheels/" + local_wheel.name
                # Requirements resolve beside the assets worker. Retain relative
                # links for sites hosted below a repository prefix.
                addressed = Path(
                    os.path.relpath(
                        output / hosted_wheel_path(local_wheel), destination
                    )
                ).as_posix()
                html = html.replace(original, addressed).replace(
                    urllib.parse.quote(original, safe=""),
                    urllib.parse.quote(addressed, safe=""),
                )
            full.write_text(html)
        runtime_url = None
        runtime_wheel = None
        cache_worker_url = None
        if runtime:
            manifest = json.loads((runtime / "runtime-manifest.json").read_text())
            site = documentation_site(output)
            for name, digest in manifest["shared_files"].items():
                source = runtime / name
                if hashlib.sha256(source.read_bytes()).hexdigest() != digest:
                    raise ValueError(f"Runtime file hash differs: {name}")
                shared = site / "assets/notebook-runtime-files" / digest / source.name
                shared.parent.mkdir(parents=True, exist_ok=True)
                if shared.exists() and shared.read_bytes() != source.read_bytes():
                    raise ValueError(f"Immutable runtime file collision: {shared}")
                if not shared.exists():
                    shutil.copyfile(source, shared)
                if name == "symbolica-core.wasm":
                    compressed = runtime / (name + ".gz")
                    expected = manifest["files"][name + ".gz"]
                    if hashlib.sha256(compressed.read_bytes()).hexdigest() != expected:
                        raise ValueError("Compressed native runtime hash differs")
                    encoded = shared.with_suffix(shared.suffix + ".gz")
                    if (
                        encoded.exists()
                        and encoded.read_bytes() != compressed.read_bytes()
                    ):
                        raise ValueError("Immutable compressed runtime collision")
                    if not encoded.exists():
                        shutil.copyfile(compressed, encoded)
            destination = site / "assets/notebook-runtime" / manifest["identity"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            if not destination.exists():
                with tempfile.TemporaryDirectory(
                    dir=destination.parent, prefix=".runtime-"
                ) as temporary:
                    candidate = Path(temporary) / "runtime"
                    shutil.copytree(runtime, candidate, copy_function=shutil.copyfile)
                    for directory in (candidate, *candidate.rglob("*")):
                        if directory.is_dir():
                            directory.chmod(0o755)
                    candidate.rename(destination)
            if (destination / "runtime-manifest.json").read_bytes() != (
                runtime / "runtime-manifest.json"
            ).read_bytes():
                raise ValueError(
                    "Shared immutable runtime manifest differs from the prepared input"
                )
            runtime_url = Path(os.path.relpath(destination, output)).as_posix() + "/"
            wheel_hash = manifest["shared_files"]["wheels/" + wheel.name]
            runtime_wheel = Path(
                os.path.relpath(
                    site / "assets/notebook-runtime-files" / wheel_hash / wheel.name,
                    output,
                )
            ).as_posix()
            cache_source = EXAMPLES_DIR / "notebook-cache-worker.js"
            cache_identity = hashlib.sha256(cache_source.read_bytes()).hexdigest()
            cache_worker = site / f"notebook-cache-worker-{cache_identity}.js"
            if (
                cache_worker.exists()
                and cache_worker.read_bytes() != cache_source.read_bytes()
            ):
                raise ValueError("Immutable cache worker differs from its source")
            if not cache_worker.exists():
                shutil.copyfile(cache_source, cache_worker)
            cache_worker_url = Path(os.path.relpath(cache_worker, output)).as_posix()
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
                        "full_notebook": "full/" + notebook.output_name,
                        "head": generator.render_head(),
                        "runtime": runtime_url + "main.js" if runtime_url else None,
                        "cache_worker": cache_worker_url,
                        "runtime_wheel": runtime_wheel,
                        "body": body,
                        "full_sources": full_sources[notebook.filename],
                        "setup_cells": setup_cells.get(notebook.filename, [0]),
                        "hidden_source_cells": hidden_source_cells.get(
                            notebook.filename, []
                        ),
                        "prose_cells": prose_cells.get(notebook.filename, []),
                        "wheel": hosted_wheel_path(wheel).as_posix(),
                        "dependency_wheel": next(
                            (
                                hosted_wheel_path(item).as_posix()
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
        hosted_wheel = output / (
            hosted_wheel_path(wheel) if docs else Path("public/wheels") / wheel.name
        )
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
                if payload["wheel"] != hosted_wheel_path(wheel).as_posix():
                    raise RuntimeError(f"{artifact.name} references the wrong wheel")
                continue
            required = {
                "embedded notebook source": "<marimo-code",
                "edit-mode configuration": '"mode": "edit"',
                "render readiness marker": (
                    f'data-notebook-ready=\\"{Path(notebook.filename).stem}\\"'
                    if notebook.docs_product != "linnet"
                    else f'data-linnet-render-ready=\\"{notebook.ready_value}\\"'
                ),
            }
            if notebook.ready_value in {"generic", "physics"}:
                required.pop("render readiness marker")
                required["prepared diagram view"] = ".render("
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


def notebook_loader_smoke(browser, base_url: str) -> None:
    """Exercise the site's loader against controlled startup and failure states.

    The normal browser smoke below executes the exported notebooks. These routes
    isolate network/runtime failures without depending on a CDN outage or wheel.
    """
    from playwright.sync_api import expect

    fixture_url = urllib.parse.urljoin(base_url, "notebook-loader-check/")
    runtime_url = urllib.parse.urljoin(fixture_url, "runtime.js")
    code = "print('<ready>')"
    setup_code = html.escape(json.dumps("await _micropip.install([])"), quote=True)
    payload = {
        "body": (
            '<marimo-island data-reactive="true" data-app-id="loader-check">'
            '<marimo-ui-element object-id="loader-setup"><marimo-code-editor '
            f"data-initial-value='{setup_code}'></marimo-code-editor></marimo-ui-element>"
            "</marimo-island>"
            '<marimo-island data-reactive="true" data-app-id="loader-check">'
            "<marimo-cell-output><p>Saved output</p></marimo-cell-output>"
            '<marimo-ui-element object-id="loader-check"><marimo-code-editor '
            f"data-initial-value='{html.escape(json.dumps(code), quote=True)}'>"
            "</marimo-code-editor></marimo-ui-element></marimo-island>"
            '<marimo-island data-reactive="true" data-app-id="loader-check">'
            '<marimo-cell-code>mo.md(f"Dynamic prose {value}")</marimo-cell-code>'
            '<marimo-cell-output><div class="markdown prose"><h2>Prose heading</h2>'
            "<p>Dynamic prose</p></div></marimo-cell-output></marimo-island>"
        ),
        "head": (
            f'<script type="module" src="{runtime_url}"></script>'
            '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Lora">'
        ),
        "wheel": "unused.whl",
        "prose_cells": [2],
    }
    document = f"""<!doctype html><html><head>
      <link rel="stylesheet" href="{base_url}assets/site.css">
      <script defer src="{base_url}assets/site.js"></script></head>
      <body data-docs-root="{base_url}"><button data-theme-toggle>Theme</button>
      <main class="docs-article">
      <h2 id="guide-heading">Guide heading</h2>
      <div class="live-notebook" data-notebook="loader-first">
        <p class="live-notebook-fallback">Source link</p></div>
      <div class="live-notebook" data-notebook="loader-second">
        <p class="live-notebook-fallback">Source link</p></div>
      </main></body></html>"""
    runtime = """
      export const initialized = Promise.resolve();
      export async function initialize() {}
      window.runtimeStarts = (window.runtimeStarts || 0) + 1;
      await new Promise(resolve => window.finishNotebookRuntime = resolve);
      customElements.define('marimo-code-editor', class extends HTMLElement {
        connectedCallback() {
          const shadow = this.attachShadow({mode: 'open'});
          const code = document.createElement('div');
          code.className = 'cm-content';
          code.contentEditable = 'true';
          code.textContent = JSON.parse(this.dataset.initialValue);
          shadow.append(code);
        }
      });
    """
    for failed_import in (False, True, "initialize"):
        page = browser.new_page(
            viewport={"width": 390, "height": 1200}, color_scheme="light"
        )
        errors: list[str] = []
        page.on("pageerror", lambda error, errors=errors: errors.append(str(error)))

        def route(request, *, failed_import=failed_import):
            url = request.request.url
            if url == fixture_url:
                request.fulfill(body=document, content_type="text/html")
            elif url == runtime_url:
                if failed_import is True:
                    request.abort()
                else:
                    source = (
                        runtime.replace(
                            "export async function initialize() {}",
                            "export async function initialize() {throw new Error('startup probe failure')}",
                        )
                        if failed_import == "initialize"
                        else runtime
                    )
                    request.fulfill(body=source, content_type="text/javascript")
            elif url.endswith(("/loader-first.json", "/loader-second.json")):
                request.fulfill(json=payload)
            else:
                request.continue_()

        page.route("**/*", route)
        try:
            page.goto(fixture_url, wait_until="domcontentloaded")
            previews = page.locator(
                ".notebook-cell:not([data-notebook-setup]) .live-notebook-code-preview"
            )
            expect(previews).to_have_count(2)
            expect(previews.first).to_be_visible()
            expect(previews.first).to_have_text(code)
            expect(page.get_by_text("Saved output").first).to_be_visible()
            expect(page.locator('link[href*="family=Lora"]')).to_have_count(0)
            expect(page.locator(".notebook-cell[hidden]")).to_have_count(2)
            expect(
                page.locator(".notebook-cell[hidden] .notebook-cell-state")
            ).to_have_count(2)
            toggle = page.locator(".notebook-toolbar button[aria-expanded]").first
            page.locator(".notebook-more > summary").first.click()
            toggle.click()
            expect(toggle).to_have_attribute("aria-expanded", "true")
            setup_source = page.locator(
                "[data-notebook-setup] .live-notebook-code-preview"
            ).first
            expect(setup_source).to_be_visible()
            expect(setup_source).to_have_text("await _micropip.install([])")
            toggle.click()
            expect(setup_source).to_be_hidden()
            assert page.locator(".notebook-cell[hidden]").first.evaluate(
                "e => e.getBoundingClientRect().height === 0"
            )
            if failed_import:
                if failed_import == "initialize":
                    page.wait_for_function("Boolean(window.finishNotebookRuntime)")
                    page.evaluate("window.finishNotebookRuntime()")
                expect(page.locator('[data-notebook-state="error"]')).to_have_count(2)
                page.locator(".notebook-more").first.evaluate("e => e.open = true")
                expect(
                    page.get_by_role("button", name="Restart original").first
                ).to_be_enabled()
                expect(
                    page.locator('.live-notebook-cells[aria-busy="true"]')
                ).to_have_count(0)
                expect(page.get_by_text("Source link").first).to_be_visible()
                if errors:
                    raise RuntimeError("Notebook loader errors: " + "; ".join(errors))
                continue
            expect(page.locator('[data-notebook-state="loading"]')).to_have_count(2)
            page.wait_for_function("window.runtimeStarts === 1")
            expect(page.locator('[data-state="ready"]')).to_have_count(0)
            page.evaluate("window.finishNotebookRuntime()")
            expect(page.locator('[data-notebook-state="active"]')).to_have_count(2)
            expect(previews.first).to_be_hidden()
            expect(page.locator('[data-state="ready"]')).to_have_count(0)
            prose = page.locator("[data-notebook-prose]").first
            expect(prose.locator(".notebook-cell-header")).to_be_hidden()
            expect(
                prose.locator("marimo-code-editor, .notebook-readonly-source")
            ).to_have_count(0)
            expect(prose.get_by_text("Dynamic prose", exact=True)).to_be_visible()
            assert prose.evaluate(
                """el => {
                    const text = getComputedStyle(el.querySelector('.prose'));
                    const guide = getComputedStyle(el.closest('.live-notebook'));
                    const heading = getComputedStyle(el.querySelector('h2'));
                    const guideHeading = getComputedStyle(document.querySelector('#guide-heading'));
                    return ['fontFamily', 'fontSize', 'lineHeight'].every(key => text[key] === guide[key])
                        && ['fontFamily', 'fontSize', 'fontWeight', 'lineHeight', 'color'].every(key => heading[key] === guideHeading[key]);
                }"""
            )
            island = page.locator(".notebook-cell:not([hidden]) marimo-island").first
            status = page.locator(
                '[data-notebook="loader-first"] .notebook-toolbar > [role="status"]'
            )
            island.evaluate("el => el.dataset.status = 'running'")
            expect(status).to_have_attribute("data-state", "busy")
            expect(
                page.get_by_role("button", name="Restart original").first
            ).to_be_enabled()
            island.evaluate("el => el.dataset.status = 'idle'")
            expect(status).to_have_attribute("data-state", "ready")
            expect(island.locator("..").locator(".notebook-cell-state")).to_have_text(
                "finished"
            )
            island.locator(".cm-content").fill("print('changed')")
            expect(status).to_have_attribute("data-state", "stale")
            page.locator(".notebook-more > summary").first.click()
            toggle.click()
            toggle.click()
            expect(island.locator(".cm-content")).to_have_text("print('changed')")
            expect(status).to_have_attribute("data-state", "stale")
            source_toggle = page.locator(".notebook-source-toggle").first
            source_toggle.click()
            expect(island.locator(".cm-content")).to_be_hidden()
            expect(page.get_by_text("Saved output").first).to_be_visible()
            source_toggle.click()
            expect(island.locator(".cm-content")).to_have_text("print('changed')")
            expect(status).to_have_attribute("data-state", "stale")
            island.evaluate("el => el.dataset.status = 'running'")
            expect(status).to_have_attribute("data-state", "busy")
            island.evaluate("el => el.dataset.status = 'idle'")
            # Error rendering can change a class without inserting a new node.
            island.locator("p").evaluate("el => el.classList.add('text-error')")
            expect(status).to_have_attribute("data-state", "error")
            setup = page.locator("[data-notebook-setup] marimo-island").first
            setup.evaluate(
                "el => el.innerHTML = '<div class=\"text-error\">setup probe failure</div>'"
            )
            expect(status).to_contain_text("Notebook setup failed: setup probe failure")
            page.locator(".notebook-more > summary").first.click()
            page.get_by_role("button", name="Show setup cells").first.click()
            expect(setup.locator(".text-error")).to_be_visible()
            page.locator("[data-theme-toggle]").click()
            expect(page.locator("marimo-code-editor").first).to_have_attribute(
                "data-theme", '"dark"'
            )
            if errors:
                raise RuntimeError("Notebook loader errors: " + "; ".join(errors))
        finally:
            page.close()


def browser_smoke(
    base_url: str,
    artifacts: Sequence[tuple[Notebook, Path]],
    timeout_seconds: float,
    browser_executable: Path | None,
    *,
    docs: bool = False,
    cache_report: Path | None = None,
) -> None:
    try:
        from playwright.sync_api import expect, sync_playwright
    except ImportError as error:
        raise RuntimeError(
            "--browser-smoke requires Playwright; install it and Chromium "
            "before running this check"
        ) from error

    timeout = timeout_seconds * 1000

    def wait_for_cells(page, artifact_name: str) -> None:
        page.wait_for_function(
            """() => {
                const failures = [...document.querySelectorAll('marimo-island .text-error')];
                // Dependent cells may render their ancestor error before the
                // originating cell renders its exception and diagnostic message.
                if (failures.some(error => error.textContent.trim() !== 'Ancestor raised'))
                    return true;
                if (document.querySelector('[data-notebook-state="error"]')) return true;
                if (failures.length || !window.notebookRunCompleted) return false;
                const cells = [...document.querySelectorAll('marimo-island')];
                return cells.length > 0 && cells.every(cell =>
                    window.notebookExecutedCells.has(cell.cellId)
                    && cell.dataset.status === 'idle');
            }""",
            timeout=timeout,
        )
        runtime_errors = page.locator(
            "marimo-island .text-error, "
            '[data-notebook] .notebook-toolbar > [role="status"][data-state="error"]'
        ).all_text_contents()
        if runtime_errors:
            # Marimo's error label contains only the exception type; retain the
            # cell's traceback and message so notebook regressions are actionable.
            runtime_errors = (
                page.locator(
                    'marimo-island [role="alert"]:has(.text-error)'
                ).all_inner_texts()
                or runtime_errors
            )
            raise RuntimeError(
                f"Notebook errors while executing {artifact_name}: "
                + "; ".join(runtime_errors)
            )

    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=str(browser_executable) if browser_executable else None
        )
        try:
            if docs:
                notebook_loader_smoke(browser, base_url)
            for notebook, artifact in artifacts:
                page = browser.new_page()
                started = time.perf_counter()
                if docs:
                    # React can batch a fast cell's complete run into its initial
                    # idle DOM state. Observe kernel completion messages instead.
                    page.add_init_script(EXECUTION_TRACKING)

                errors: list[str] = []
                wheel_requests: list[str] = []
                page.on(
                    "request",
                    lambda request, urls=wheel_requests: (
                        urls.append(request.url)
                        if "symbolica-" in request.url and request.url.endswith(".whl")
                        else None
                    ),
                )
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
                    wait_for_cells(page, artifact.name)
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
                if notebook.ready_value in {"generic", "physics"}:
                    page.frame_locator(notebook.ready_selector).locator(
                        "svg[data-linnet-renderer=native]"
                    ).wait_for(state="visible", timeout=timeout)
                if docs:
                    page.wait_for_function(
                        """() => document.querySelector(
                            '[data-notebook] .notebook-toolbar > [role="status"]'
                        )?.dataset.state === 'ready'""",
                        timeout=timeout,
                    )
                    # Compare displayed durations with the worker's own kernel
                    # timestamps, independently of the UI store subscription.
                    page.wait_for_function(
                        """() => [...document.querySelectorAll('.live-notebook .notebook-cell-time')].every(time => {
                            const island = time.closest('.notebook-cell').querySelector('marimo-island');
                            const elapsed = window.notebookCellTimings.get(island.cellId);
                            // mo.stop can leave dependent cells idle without executing them.
                            if (elapsed === undefined) return time.hidden && !time.hasAttribute('data-elapsed-ms');
                            return Number.isFinite(elapsed) && !time.hidden &&
                                Number(time.dataset.elapsedMs) === elapsed;
                        })""",
                        timeout=timeout,
                    )
                    assert (
                        page.locator(".live-notebook .notebook-cell-time").count() > 0
                    )
                    assert not page.locator(
                        "[data-notebook-prose] .notebook-cell-time"
                    ).count()
                    # Inspect inside the real editor shadow root. Host font
                    # variables alone do not prove the editor adopted them.
                    for theme in ("dark", "light"):
                        page.evaluate(
                            "theme => document.documentElement.dataset.theme = theme",
                            theme,
                        )
                        page.wait_for_function(
                            """() => [...document.querySelectorAll('marimo-code-editor')].every(editor => {
                                const shadow = editor.shadowRoot;
                                const scroller = shadow?.querySelector('.cm-scroller');
                                if (!scroller) return false;
                                if (shadow.querySelector('.cm-content').textContent.trim() &&
                                    !shadow.querySelector('[class*="tok-"]')) return false;
                                const style = getComputedStyle(scroller);
                                const keyword = shadow.querySelector('.tok-keyword');
                                const expected = getComputedStyle(editor).getPropertyValue('--syntax-keyword').trim();
                                const probe = document.createElement('span');
                                probe.style.color = expected;
                                document.body.append(probe);
                                const color = getComputedStyle(probe).color;
                                probe.remove();
                                return style.fontFamily.startsWith('"IBM Plex Mono"') &&
                                    style.fontSize === '14.4px' && style.lineHeight === '22.32px' &&
                                    (!keyword || getComputedStyle(keyword).color === color);
                            })""",
                            timeout=timeout,
                        )
                    setup = page.locator("[data-notebook-setup]")
                    for cell in setup.all():
                        expect(cell).to_be_hidden()
                    assert setup.evaluate_all(
                        "cells => cells.every(e => e.getBoundingClientRect().height === 0)"
                    )
                    toolbar = page.get_by_role(
                        "region", name="Notebook runtime controls"
                    )
                    expect(toolbar.get_by_role("status")).to_be_visible()
                    expect(
                        toolbar.get_by_role("button", name="Run all", exact=True)
                    ).to_be_visible()
                    expect(toolbar.locator(".notebook-help")).to_be_hidden()
                    assert toolbar.evaluate(
                        "el => getComputedStyle(el).position === 'fixed' && "
                        "Math.abs(el.getBoundingClientRect().bottom - innerHeight) < 1"
                    )
                    assert page.locator(".page-footer").evaluate(
                        "el => !['fixed', 'sticky'].includes(getComputedStyle(el).position)"
                    )
                    summary = toolbar.locator("summary")
                    summary.focus()
                    summary.press("Enter")
                    expect(
                        toolbar.get_by_role("button", name="Restart original")
                    ).to_be_enabled()
                    summary.press("Escape")
                    expect(toolbar.locator("details")).not_to_have_attribute("open", "")
                    expect(summary).to_be_focused()
                    completed = page.evaluate("window.notebookCompletedRuns")
                    toggle = page.locator(".notebook-toolbar button[aria-expanded]")
                    page.locator(".notebook-more > summary").click()
                    toggle.focus()
                    toggle.press("Control+k")
                    expect(page.locator("[data-search-dialog]")).to_be_visible()
                    page.locator("[data-search-input]").fill("notebook")
                    page.keyboard.press("Escape")
                    expect(toggle).to_be_focused()
                    expect(toggle).to_be_visible()
                    toggle.click()
                    expect(toggle).to_have_attribute("aria-expanded", "true")
                    for editor in setup.locator("marimo-code-editor").all():
                        expect(editor).to_be_visible()
                        expect(editor.locator(".cm-content")).to_have_attribute(
                            "contenteditable", "false"
                        )
                    assert "_micropip.install(" in json.loads(
                        setup.locator("marimo-code-editor").first.get_attribute(
                            "data-initial-value"
                        )
                    )
                    expect(setup.locator(".notebook-cell-state")).to_have_text(
                        ["finished"] * setup.count()
                    )
                    toggle.click()
                    for cell in setup.all():
                        expect(cell).to_be_hidden()
                    assert setup.evaluate_all(
                        "cells => cells.every(e => e.getBoundingClientRect().height === 0)"
                    )
                    assert page.evaluate("window.notebookCompletedRuns") == completed
                    source_cell = page.locator(
                        '.notebook-cell:has(.notebook-source-toggle) > marimo-island:not([data-notebook-source-hidden="true"]):has(> marimo-ui-element > marimo-code-editor):has(.output)'
                    ).first
                    if source_cell.count():
                        source_cell = page.locator(
                            f'marimo-island[data-cell-idx="{source_cell.get_attribute("data-cell-idx")}"]'
                        )
                        editor = source_cell.locator(
                            ":scope > marimo-ui-element > marimo-code-editor"
                        )
                        original_source = editor.get_attribute("data-initial-value")
                        source_toggle = source_cell.locator("..").locator(
                            ".notebook-source-toggle"
                        )
                        source_toggle.click()
                        expect(editor).to_be_hidden()
                        expect(
                            source_cell.locator(":scope > .relative > .output")
                        ).to_be_visible()
                        source_toggle.click()
                        expect(editor).to_be_visible()
                        assert (
                            editor.get_attribute("data-initial-value")
                            == original_source
                        )
                        assert (
                            page.evaluate("window.notebookCompletedRuns") == completed
                        )
                    # The numerical topology picker is authored with hide_code=True.
                    # Require its native editor before exercising the upstream Run.
                    if notebook.ready_value == "08_numerical_cff_ltd_marimo":
                        hidden = page.locator('marimo-island[data-cell-idx="9"]')
                        expect(hidden).to_have_attribute(
                            "data-notebook-source-hidden", "true"
                        )
                        editor = hidden.locator(
                            ":scope > marimo-ui-element > marimo-code-editor"
                        )
                        expect(editor).to_have_count(1)
                        expect(editor).to_be_hidden()
                        wrapper = hidden.locator("..")
                        toggle = wrapper.locator(".notebook-source-toggle")
                        toggle.click()
                        expect(editor).to_be_visible()
                        content = editor.locator(".cm-content")
                        original = json.loads(
                            editor.get_attribute("data-initial-value")
                        )
                        expect(
                            editor.get_by_role(
                                "button", name="Copy to clipboard", exact=True
                            )
                        ).to_have_count(1)
                        assert editor.locator('[class*="tok-"]').count() > 0
                        page.wait_for_function(
                            """editor => {
                            const shadow = editor.shadowRoot;
                            const number = [...shadow.querySelectorAll('.cm-lineNumbers .cm-gutterElement')].find(n => n.style.visibility !== 'hidden');
                            if (!number) return false;
                            const range = document.createRange();
                            range.selectNodeContents(number);
                            const frame = editor.closest('.notebook-cell').getBoundingClientRect();
                            return Math.abs(frame.left - range.getBoundingClientRect().right - 8) < 1;
                        }""",
                            arg=editor.element_handle(),
                        )
                        content.fill(original + "\n# hidden source regression")
                        content.press("End")
                        selection = content.evaluate(
                            "el => el.getRootNode().getSelection().anchorOffset"
                        )
                        for modifier in ("Meta", "Control"):
                            content.press(modifier + "+k")
                            expect(page.locator("[data-search-dialog]")).to_be_visible()
                            page.locator("[data-search-input]").fill("numerical")
                            page.keyboard.press("Escape")
                            expect(page.locator("[data-search-dialog]")).to_be_hidden()
                            expect(content).to_be_focused()
                            page.wait_for_function(
                                "({element, offset}) => element.getRootNode().getSelection().anchorOffset === offset",
                                arg={
                                    "element": content.element_handle(),
                                    "offset": selection,
                                },
                                timeout=5000,
                            )
                        content.press("/")
                        expect(page.locator("[data-search-dialog]")).to_be_hidden()
                        content.press("Backspace")
                        page.keyboard.down("Meta")
                        for copy in hidden.locator(
                            ":scope > .relative > .absolute button:has(svg.lucide-copy)"
                        ).all():
                            expect(copy).to_be_hidden()
                        page.keyboard.up("Meta")
                        toggle.click()
                        expect(editor).to_be_hidden()
                        toggle.click()
                        expect(content).to_contain_text("# hidden source regression")
                        completed_hidden = page.evaluate("window.notebookCompletedRuns")
                        wrapper.locator(".notebook-run").click()
                        page.wait_for_function(
                            "n => window.notebookCompletedRuns > n",
                            arg=completed_hidden,
                            timeout=timeout,
                        )
                        wait_for_cells(page, artifact.name)
                        content.fill(original)
                        wrapper.locator(".notebook-run").click()
                        page.wait_for_function(
                            "n => window.notebookCompletedRuns > n + 1",
                            arg=completed_hidden,
                            timeout=timeout,
                        )
                        wait_for_cells(page, artifact.name)
                        toggle.click()
                    page.evaluate("() => { document.activeElement?.blur(); }")
                    page.keyboard.press("Meta+k")
                    expect(page.locator("[data-search-dialog]")).to_be_visible()
                    page.keyboard.press("Escape")
                    page.evaluate("""() => {
                        const input = document.createElement('input');
                        input.dataset.notebookSearchProbe = '';
                        input.value = 'preserve form text';
                        document.body.append(input);
                        input.focus();
                    }""")
                    form_input = page.locator("[data-notebook-search-probe]")
                    page.keyboard.press("Control+k")
                    expect(page.locator("[data-search-dialog]")).to_be_visible()
                    page.keyboard.press("Escape")
                    expect(form_input).to_be_focused()
                    expect(form_input).to_have_value("preserve form text")
                    form_input.evaluate("el => el.remove()")
                    for prose in page.locator("[data-notebook-prose]").all():
                        expect(prose.locator("marimo-code-editor")).to_have_count(0)
                        expect(prose.locator(".notebook-cell-header")).to_be_hidden()
                        expect(
                            prose.locator(".notebook-readonly-source")
                        ).to_have_count(0)
                        assert prose.evaluate(
                            """el => {
                                const text = getComputedStyle(el.querySelector('.prose, .markdown') || el);
                                const guide = getComputedStyle(el.closest('.live-notebook'));
                                const output = el.querySelector(':scope > marimo-island > .relative, :scope > marimo-island > marimo-cell-output');
                                const frame = output && getComputedStyle(output);
                                return ['fontFamily', 'fontSize', 'lineHeight'].every(key => text[key] === guide[key])
                                    && (!frame || ['paddingLeft', 'paddingRight', 'paddingTop', 'paddingBottom', 'borderTopWidth'].every(key => frame[key] === '0px'));
                            }"""
                        )
                    # The native outset gutter belongs in the shared guide
                    # margin; it must never indent only notebook prose or frames.
                    assert page.locator(".live-notebook-cells").evaluate(
                        """el => {
                            const guide = el.closest('.docs-article');
                            if (!guide) return true;
                            const column = guide.getBoundingClientRect();
                            return [...el.querySelectorAll('.notebook-cell')]
                                .filter(cell => !cell.hidden)
                                .every(cell => {
                                    const frame = cell.getBoundingClientRect();
                                    return Math.abs(frame.left - column.left) < 1
                                        && Math.abs(frame.right - column.right) < 1;
                                });
                        }"""
                    )
                    # Snapshot actual code/output pairs in one DOM read. Import
                    # cells can lose their empty output while editors hydrate.
                    pairs = page.locator("marimo-island").evaluate_all("""cells =>
                        cells.flatMap(cell => {
                            const editor = cell.querySelector(
                                ':scope > marimo-ui-element:has(> marimo-code-editor)');
                            const output = cell.querySelector(':scope > .relative > .output');
                            if (!editor || !output) return [];
                            const controls = cell.parentElement.querySelector('.notebook-cell-header');
                            const controlsBox = controls.getBoundingClientRect();
                            const codeBox = editor.getBoundingClientRect();
                            const outputBox = output.getBoundingClientRect();
                            return codeBox.height && outputBox.height
                                ? [{codeY: codeBox.y, codeBottom: codeBox.bottom,
                                    controlsY: controlsBox.y, controlsBottom: controlsBox.bottom,
                                    outputY: outputBox.y}] : [];
                        })""")
                    if any(
                        pair["codeBottom"] > pair["controlsY"] + 1
                        or pair["controlsBottom"] > pair["outputY"] + 1
                        for pair in pairs
                    ):
                        raise RuntimeError(
                            f"{artifact.name} does not present source, controls, then output"
                        )
                if docs and notebook.ready_value in {
                    "07_tensor_reduction_marimo",
                }:
                    installed = len(wheel_requests)
                    completed = page.evaluate("window.notebookCompletedRuns")
                    page.get_by_role("button", name="Run all", exact=True).click()
                    page.wait_for_function(
                        "previous => window.notebookCompletedRuns > previous",
                        arg=completed,
                        timeout=timeout,
                    )
                    wait_for_cells(page, artifact.name)
                    assert len(wheel_requests) == installed, (
                        "Run all downloaded the installed wheel again"
                    )
                if docs and notebook.ready_value == "02_cff_and_symbolica_marimo":
                    # Exercise real upstream execution, not fabricated data-status.
                    run_all = page.get_by_role("button", name="Run all", exact=True)
                    expect(run_all).to_be_enabled()
                    cell = page.locator(
                        '.notebook-cell:not([data-notebook-setup]) > marimo-island:has(marimo-code-editor):has(button[aria-label="Run cell"])'
                    ).first
                    editor = cell.locator("marimo-code-editor")
                    code = json.loads(editor.get_attribute("data-initial-value"))
                    assert "notebook_browser_ready" not in code
                    content = editor.locator(".cm-content")
                    content.fill(
                        code
                        + "\nimport asyncio as _notebook_probe_asyncio"
                        + "\nawait _notebook_probe_asyncio.sleep(0.25)"
                        + "\nraise ValueError('notebook control regression probe')"
                    )
                    expect(
                        page.locator(".live-notebook .notebook-cell-state").filter(
                            has_text="edited"
                        )
                    ).to_have_count(1)
                    wrapper = cell.locator("..")
                    source_toggle = wrapper.locator(".notebook-source-toggle")
                    source_toggle.click()
                    expect(editor).to_be_hidden()
                    cell.evaluate("""el => {
                        window.cellExecutionStates = [];
                        const observer = new MutationObserver(() => {
                            window.cellExecutionStates.push({
                                state: el.dataset.status,
                                disabled: el.parentElement.querySelector('.notebook-run').disabled,
                                busy: el.parentElement.getAttribute('aria-busy'),
                                timingHidden: el.parentElement.querySelector('.notebook-cell-time').hidden
                            });
                        });
                        observer.observe(el, {attributes: true, attributeFilter: ['data-status']});
                    }""")
                    wrapper.get_by_role("button", name="Run cell", exact=True).click()
                    expect(
                        page.locator(
                            '[data-notebook] .notebook-toolbar > [role="status"]'
                        )
                    ).to_have_attribute("data-state", "error", timeout=timeout)
                    expect(wrapper.locator(".notebook-cell-state")).to_have_text(
                        "error"
                    )
                    states = page.evaluate("window.cellExecutionStates")
                    busy_states = [
                        entry
                        for entry in states
                        if entry["state"] in {"queued", "running"}
                    ]
                    assert busy_states and all(
                        entry["disabled"]
                        and entry["busy"] == "true"
                        and entry["timingHidden"]
                        for entry in busy_states
                    )
                    expect(wrapper.locator(".notebook-cell-time")).to_be_visible()
                    assert (
                        float(
                            wrapper.locator(".notebook-cell-time").get_attribute(
                                "data-elapsed-ms"
                            )
                        )
                        >= 250
                    )
                    expect(
                        wrapper.get_by_role("button", name="Run cell", exact=True)
                    ).to_be_enabled()
                    source_toggle.click()
                    expect(editor.locator(".cm-content")).to_contain_text(
                        "notebook control regression probe"
                    )
                    # Run all preserves edits; Restart original intentionally discards them.
                    expect(run_all).to_be_enabled()
                    page.evaluate(
                        "window.notebookExecutedCells.clear(); window.notebookRunCompleted = false"
                    )
                    installed_wheel_requests = len(wheel_requests)
                    run_all.click()
                    expect(page.locator("marimo-island .text-error")).to_have_count(
                        0, timeout=timeout
                    )
                    expect(
                        page.get_by_text(
                            "notebook control regression probe", exact=False
                        ).first
                    ).to_be_visible(timeout=timeout)
                    expect(
                        page.locator(
                            '[data-notebook] .notebook-toolbar > [role="status"]'
                        )
                    ).to_have_attribute("data-state", "error", timeout=timeout)
                    assert len(wheel_requests) == installed_wheel_requests, (
                        "Run all downloaded the installed wheel again"
                    )
                    expect(run_all).to_be_enabled()
                    page.evaluate(
                        "window.notebookExecutedCells.clear(); window.notebookRunCompleted = false"
                    )
                    # Pinned reset must restart even when lazy notebook content
                    # would be outside the restored viewport after a reload.
                    page.evaluate("scrollTo(0, document.documentElement.scrollHeight)")
                    page.locator(".notebook-more > summary").click()
                    page.get_by_role(
                        "button", name="Restart original", exact=True
                    ).click()
                    expect(page.locator("marimo-island .text-error")).to_have_count(
                        0, timeout=timeout
                    )
                    wait_for_cells(page, artifact.name)
                    expect(
                        page.locator(
                            '[data-notebook] .notebook-toolbar > [role="status"]'
                        )
                    ).to_have_attribute("data-state", "ready", timeout=timeout)
                    expect(
                        page.locator(".live-notebook .notebook-cell-state").first
                    ).to_have_text("finished")
                    full = page.get_by_role("link", name="Open full notebook")
                    assert full.get_attribute("href").endswith(notebook.output_name)
                    full_page = browser.new_page()
                    try:
                        full_page.goto(
                            full.get_attribute("href"), wait_until="domcontentloaded"
                        )
                        full_page.locator(".cm-editor").first.wait_for(timeout=timeout)
                        full_page.locator('[data-testid="run-button"]').last.click()
                        full_page.locator(notebook.ready_selector).wait_for(
                            timeout=timeout
                        )
                        assert not full_page.locator(".text-error").count()
                    finally:
                        full_page.close()
                if docs and notebook.ready_value == "spenso_tensor_explorer":
                    component = page.get_by_role("slider").first
                    completed = page.evaluate("window.notebookCompletedRuns")
                    component.press("End")
                    page.wait_for_function(
                        "previous => window.notebookCompletedRuns > previous",
                        arg=completed,
                        timeout=timeout,
                    )
                    wait_for_cells(page, artifact.name)
                    expect(
                        page.get_by_text(
                            "Component 255 of the 256-entry current", exact=False
                        )
                    ).to_be_visible()
                if docs and notebook.ready_value == "spenso_api_tour":
                    full_page = browser.new_page()
                    try:
                        full_page.goto(
                            page.get_by_role(
                                "link", name="Open full notebook"
                            ).get_attribute("href")
                        )
                        full_page.locator(".cm-editor").first.wait_for(timeout=timeout)
                        full_page.locator('[data-testid="run-button"]').last.click()
                        full_page.locator(notebook.ready_selector).wait_for(
                            timeout=timeout
                        )
                        expect(full_page.locator(".text-error")).to_have_count(0)
                    finally:
                        full_page.close()
                if docs and notebook.ready_value == "05_four_loop_numerator_marimo":
                    picker = page.locator("marimo-dropdown")
                    # Names and particle content identify these representatives;
                    # historical dropdown positions are not assertions.
                    for name, content in (
                        ("D2603", "11 g"),
                        ("D0713", "2 d, 3 g, 6 ghG"),
                        ("D2881", "2 d, 9 g"),
                        ("D2830", "10 g"),
                        ("D0000", "3 g, 8 ghG"),
                    ):
                        labels = json.loads(picker.get_attribute("data-options"))
                        label = next(
                            value
                            for value in labels
                            if f": {name} — {content}" in value
                        )
                        completed = page.evaluate("window.notebookCompletedRuns")
                        picker.locator("button").click()
                        page.locator("input[cmdk-input]").fill(name)
                        page.get_by_role("option", name=label, exact=True).click()
                        page.wait_for_function(
                            "previous => window.notebookCompletedRuns > previous",
                            arg=completed,
                            timeout=timeout,
                        )
                        wait_for_cells(page, artifact.name)
                if docs and notebook.ready_value == "08_numerical_cff_ltd_marimo":
                    # Every topology in the source picker must execute, including
                    # both three-loop examples. Select through the upstream UI.
                    picker = (
                        page.locator("marimo-dropdown").filter(has_text="Diagram").first
                    )
                    for topology in (
                        "Crossed double box · 2 loops",
                        "Triple box · 3 loops",
                        "Tetrahedron · 3 loops",
                        "Double box · 2 loops",
                    ):
                        completed = page.evaluate("window.notebookCompletedRuns")
                        picker.locator("select").select_option(label=topology)
                        page.wait_for_function(
                            "previous => window.notebookCompletedRuns > previous",
                            arg=completed,
                            timeout=timeout,
                        )
                        wait_for_cells(page, artifact.name)
                    # Sampling must work without replacing the independent OneLOop reference.
                    page.get_by_role(
                        "button", name="Run integration comparison", exact=True
                    ).click()
                    page.get_by_text(
                        "OneLOop reference unavailable in this browser build.",
                        exact=True,
                    ).wait_for(timeout=timeout)
                    wait_for_cells(page, artifact.name)
                    expect(
                        page.get_by_text("Uniform Cartesian", exact=True).first
                    ).to_be_visible()
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
                        "Momentum labels",
                        "Cross-section",
                    ):
                        previous = page.locator(notebook.ready_selector).get_attribute(
                            "srcdoc"
                        )
                        assert previous is not None, (
                            "Expected a rendered diagram iframe"
                        )
                        if name == "Cross-section":
                            page.get_by_role(
                                "combobox", name="Example", exact=True
                            ).select_option(label=name)
                        else:
                            page.get_by_role("checkbox", name=name, exact=True).check()
                        # AnyWidget mounts its iframe inside a shadow root;
                        # Playwright locators pierce it, document.querySelector does not.
                        expect(
                            page.locator(notebook.ready_selector)
                        ).not_to_have_attribute("srcdoc", previous, timeout=timeout)
                    preview = page.frame_locator(notebook.ready_selector)
                    svg = preview.locator("svg[data-linnet-interactive]")
                    svg.wait_for(state="visible", timeout=timeout)
                    if not svg.evaluate("svg => Boolean(svg.linnetSelection)"):
                        raise RuntimeError(
                            "FeynKit preview did not initialize SVG interaction"
                        )
                    for split in (True, False, True):
                        toggle = page.get_by_role(
                            "checkbox", name="Split initial state", exact=True
                        )
                        if toggle.is_checked() != split:
                            previous = page.locator(
                                notebook.ready_selector
                            ).get_attribute("srcdoc")
                            assert previous is not None, (
                                "Expected a rendered diagram iframe"
                            )
                            toggle.set_checked(split)
                            expect(
                                page.locator(notebook.ready_selector)
                            ).not_to_have_attribute("srcdoc", previous, timeout=timeout)
                        fragments = svg.evaluate("""svg => {
                            const edges = new Map();
                            for (const target of svg.querySelectorAll('[data-linnet-kind="edge"]')) {
                                const detail = JSON.parse(target.dataset.linnetDetail);
                                if (!edges.has(detail.edge)) edges.set(detail.edge, new Set());
                                edges.get(detail.edge).add(detail.name);
                            }
                            return [...edges.values()].map(names => names.size).sort();
                        }""")
                        expected = [1, 1, 1, 1, 2, 2] if split else [1] * 6
                        if fragments != expected:
                            raise RuntimeError(
                                f"Cross-section split mode {split}: {fragments}"
                            )
                    for theme in ("dark", "light"):
                        page.locator("html").evaluate(
                            "(root, theme) => { root.dataset.theme = theme; }", theme
                        )
                        expect(svg).to_have_attribute(
                            "data-theme", theme, timeout=timeout
                        )
                        color_scheme = page.locator(notebook.ready_selector).evaluate(
                            "frame => getComputedStyle(frame).colorScheme"
                        )
                        expect(preview.locator("html")).to_have_css(
                            "color-scheme", color_scheme, timeout=timeout
                        )
                    half_edge = svg.locator('[data-linnet-kind="halfedge"]').first
                    half_id = int(half_edge.get_attribute("data-linnet-id"))
                    half_edge.dispatch_event("click", {"shiftKey": True})
                    if svg.evaluate("svg => svg.linnetSelection.half_edges") != [
                        half_id
                    ]:
                        raise RuntimeError("FeynKit preview did not select a half-edge")
                if docs:
                    wait_for_cells(page, artifact.name)
                if errors:
                    raise RuntimeError(
                        f"Browser errors while loading {artifact.name}: "
                        + "; ".join(errors)
                    )
                print(
                    f"Browser notebook passed: {artifact.name} "
                    f"({page.locator('marimo-island').count()} cells, "
                    f"{time.perf_counter() - started:.1f}s)",
                    flush=True,
                )
                page.close()
            if cache_report:
                from cache_wasm import cache_navigation_smoke

                cache_navigation_smoke(
                    browser,
                    base_url,
                    documentation_site(artifacts[0][1].parent),
                    cache_report,
                    timeout,
                    EXECUTION_TRACKING,
                    wait_for_cells,
                )
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
            "Local Emscripten Symbolica Community wheel "
            "with FeynKit/Spenso/Idenso for GammaLoop and community showcases. Without this option the "
            "published Symbolica dependency is used."
        ),
    )
    parser.add_argument(
        "--dependency-wheel",
        type=Path,
        help="Linnet WASM wheel for GammaLoop/Spenso/Idenso, or the pinned UFO loader wheel for FeynKit",
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
        "--cache-report",
        type=Path,
        help="Test cache, reload and cross-product navigation flows",
    )
    parser.add_argument("--runtime", type=Path, help="Pinned, prepared islands runtime")
    parser.add_argument(
        "--timeout",
        type=float,
        default=180,
        help="Per-page browser timeout in seconds (default: 180)",
    )
    return parser.parse_args(arguments)


def main(arguments: Sequence[str] | None = None) -> int:
    options = parse_args(arguments)
    # Detect omitted notebooks throughout the source tree, not only FeynKit.
    repository = EXAMPLES_DIR.parents[2]
    discovered = set()
    for directory in (
        repository / "examples/notebooks",
        EXAMPLES_DIR,
        repository / "crates/feynkit-py/examples",
    ):
        for path in directory.rglob("*.py"):
            if "symbolica-host" in path.parts:
                continue
            tree = ast.parse(path.read_text(encoding="utf-8"))
            aliases = {
                alias.asname or alias.name
                for node in ast.walk(tree)
                if isinstance(node, ast.Import)
                for alias in node.names
                if alias.name == "marimo"
            }
            constructors = {
                alias.asname or alias.name
                for node in ast.walk(tree)
                if isinstance(node, ast.ImportFrom) and node.module == "marimo"
                for alias in node.names
                if alias.name == "App"
            }
            if any(
                isinstance(node, ast.Call)
                and (
                    isinstance(node.func, ast.Name)
                    and node.func.id in constructors
                    or isinstance(node.func, ast.Attribute)
                    and isinstance(node.func.value, ast.Name)
                    and node.func.value.id in aliases
                    and node.func.attr == "App"
                )
                for node in ast.walk(tree)
            ):
                discovered.add(path.resolve())
    registered = {notebook.source.resolve() for notebook in NOTEBOOKS}
    if discovered != registered:
        raise ValueError(
            f"Notebook source/registry drift: missing={discovered - registered}, obsolete={registered - discovered}"
        )
    # The upstream gallery is an independent coverage requirement: local AST
    # discovery cannot detect a canonical example never copied into this tree.
    canonical_four_loop = (
        repository / "examples/notebooks/feynkit/05_four_loop_numerator_marimo.py"
    )
    if canonical_four_loop.resolve() not in registered:
        raise ValueError("Missing canonical Community four_loop_numerator.py showcase")
    canonical_source = canonical_four_loop.read_text(encoding="utf-8")
    if "39b556c7d7dcc919f61d125f60ff83f8d49f7ac3" not in canonical_source:
        raise ValueError(
            "Four-loop showcase must identify its pinned upstream revision"
        )
    required_calculations = {
        "tensor_reduce",
        "overall_factor_expression",
        "numerator_expression",
        "with_lorentz_dimension",
        "contract",
        "to_dots",
        "simplify_algebra",
    }
    actual_calls = {
        node.func.attr
        for node in ast.walk(ast.parse(canonical_source))
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    }
    if missing := required_calculations - actual_calls:
        raise ValueError(f"Four-loop showcase omitted scientific operations: {missing}")
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
                and notebook.ready_value == "spenso_idenso_display"
            )
            if options.docs
            else (options.notebook is None and notebook.docs_product == "linnet")
            or Path(notebook.filename).stem == options.notebook
        )
    )
    packages = {notebook.package for notebook in notebooks}
    if len(packages) != 1:
        raise ValueError(
            "Select --docs or --notebook so one host wheel serves the export"
        )
    package = packages.pop()
    wheel = validate_notebook_wheel(options.wheel, package) if options.wheel else None
    dependency_wheel = (
        validate_notebook_wheel(
            options.dependency_wheel,
            "ufo-model-loader",
        )
        if options.dependency_wheel
        else None
    )
    if options.docs == "feynkit" and dependency_wheel is None:
        raise ValueError(
            "The UFO showcase requires --dependency-wheel for the UFO loader"
        )
    if dependency_wheel is not None and options.docs != "feynkit":
        raise ValueError("--dependency-wheel applies only to community showcases")
    with staged_notebooks(wheel, notebooks, dependency_wheel) as staged:
        lint(staged)
        artifacts = export(staged, output, docs=options.docs, runtime=options.runtime)

    root = documentation_site(output) if options.docs else output
    with serve(root) as base_url:
        if options.docs:
            base_url = urllib.parse.urljoin(
                base_url, output.parent.parent.relative_to(root).as_posix() + "/"
            )
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
                cache_report=options.cache_report,
            )

    print(f"Marimo WASM notebooks exported to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
