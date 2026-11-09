"""Package authored documentation cells and their explicit execution contexts."""

import argparse
import ast
import base64
import builtins
import doctest
import hashlib
import json
import os
import re
import shutil
import symtable
import tempfile
import tomllib
import zipfile
from pathlib import Path

FENCES = re.compile(r"```(python|py|typst|typ)\n([\s\S]*?)```", re.MULTILINE)


def inventory(root):
    cells = []
    registry = tomllib.loads((root / "docs/products/registry.toml").read_text())
    routes = {
        page["source"]: {"product": product["id"], "route": page["route"]}
        for product in registry["product"]
        for page in product.get("pages", [])
    }
    for path in sorted((root / "docs/products").rglob("*.typ")):
        text = path.read_text()
        for index, match in enumerate(FENCES.finditer(text)):
            cells.append(
                {
                    "id": f"{path.relative_to(root).as_posix()}:{index + 1}",
                    "file": path.relative_to(root).as_posix(),
                    "line": text[: match.start()].count("\n") + 1,
                    **routes.get(path.relative_to(root).as_posix(), {}),
                    "language": "python" if match[1] in ("py", "python") else "typst",
                    "source": match[2].rstrip(),
                    "sha256": hashlib.sha256(match[2].rstrip().encode()).hexdigest(),
                }
            )
    return cells


def python_examples(content):
    parser = doctest.DocTestParser()
    fences = re.compile(r"(?m)^[ \t]*(`{3,})python[^\n]*\n([\s\S]*?)^[ \t]*\1[ \t]*$")
    cursor = 0
    for fence in fences.finditer(content):
        for paragraph in re.split(r"\n\s*\n", content[cursor : fence.start()]):
            examples = parser.get_examples(paragraph)
            if examples:
                yield "".join(example.source for example in examples).rstrip()
        code = fence[2]
        examples = parser.get_examples(code)
        yield (
            "".join(example.source for example in examples) if examples else code
        ).rstrip()
        cursor = fence.end()
    for paragraph in re.split(r"\n\s*\n", content[cursor:]):
        examples = parser.get_examples(paragraph)
        if examples:
            yield "".join(example.source for example in examples).rstrip()


def python_reference_inventory(root):
    """Retain exported API docstrings and their explicitly named class setup."""
    entries = []
    modules = {
        "linnet-graph": "linnet",
        "feynkit-community": "feynkit",
        "spynso3": "spenso",
        "vakint-community": "vakint",
        "gammaloop-python": "gammaloop",
    }
    for stub, product in modules.items():
        path = root / f"docs/api/python/{stub}.pyi"
        tree = ast.parse(path.read_text())
        classes = {
            node.name: node for node in tree.body if isinstance(node, ast.ClassDef)
        }
        class_sources = {
            name: "\n".join(python_examples(ast.get_docstring(node) or ""))
            for name, node in classes.items()
        }
        for owner in tree.body:
            nodes = [(owner, getattr(owner, "name", "module"))]
            if isinstance(owner, ast.ClassDef):
                nodes += [
                    (node, owner.name + "." + node.name)
                    for node in owner.body
                    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                ]
            for node, qualified in nodes:
                if not isinstance(
                    node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
                ):
                    continue
                doc = ast.get_docstring(node) or ""
                content = doc
                setup_class = re.search(
                    r"Using the setup in the ``([^`]+)`` class example", content
                )
                inherited = (
                    class_sources.get(owner.name, "")
                    if isinstance(owner, ast.ClassDef) and node is not owner
                    else ""
                )
                prior = class_sources.get(setup_class[1], "") if setup_class else ""
                for index, source in enumerate(python_examples(content)):
                    if not prior and inherited:
                        symbols = symtable.symtable(source, qualified, "exec")
                        needed = {
                            item.get_name()
                            for item in symbols.get_symbols()
                            if item.is_referenced()
                            and not (item.is_assigned() or item.is_imported())
                        }
                        provided = {
                            item.get_name()
                            for item in symtable.symtable(
                                inherited, qualified, "exec"
                            ).get_symbols()
                            if item.is_assigned() or item.is_imported()
                        }
                        if (needed - set(dir(builtins))) & provided:
                            prior = inherited
                    identity = f"{path.relative_to(root).as_posix()}:{qualified}:{index}:{hashlib.sha256(source.encode()).hexdigest()[:12]}"
                    entry = {
                        "id": identity,
                        "file": identity.rsplit(":", 1)[0],
                        "owner": qualified,
                        "line": node.lineno,
                        "product": product,
                        "route": f"reference/python/{stub}/{qualified.split('.')[0]}{'-function' if isinstance(owner, (ast.FunctionDef, ast.AsyncFunctionDef)) else ''}/",
                        "language": "python",
                        "source": source,
                        "sha256": hashlib.sha256(source.encode()).hexdigest(),
                        "setup": prior,
                        "depends": [],
                    }
                    if product == "vakint":
                        entry["reason"] = (
                            "Requires the native Vakint extension and external integral evaluation programs, absent from the browser host."
                        )
                    if product == "gammaloop":
                        entry["reason"] = (
                            "Requires the native gammaloop Python extension and application state/evaluator inputs; the browser Symbolica host does not contain gammaloop."
                        )
                    native_exports = {
                        "FeynmanDiagram.to_linnest",
                        "Tensor.to_svg",
                        "TensorExpression.to_svg",
                        "TensorLibrary.to_html",
                        "TensorNetwork.to_linnest",
                        "TensorNetwork.to_svg",
                        "to_svg",
                    }
                    if product != "vakint" and qualified in native_exports:
                        entry["reason"] = (
                            "Requires the native typst-py extension for outlined vector exports, unavailable in Pyodide."
                        )
                    if 'browser="chromium"' in (prior + source):
                        entry["reason"] = (
                            "Requires a local Playwright Chromium process for explicit headless measurement; direct browser cells cannot spawn native processes."
                        )
                    if qualified == "DiagramView.widget":
                        entry["reason"] = (
                            "Requires a Jupyter or Marimo AnyWidget frontend; direct cells display the measured native graph HTML."
                        )
                    if qualified.split(".")[0] in {
                        "LoadedModel",
                        "UfoLoadDiagnostics",
                        "UfoLoader",
                    }:
                        entry["reason"] = (
                            "Requires an external UFO model directory and the native UFO loader host; this model is not included in the browser filesystem."
                        )
                    if qualified in {
                        "FeynmanDiagram.sector_decompose",
                        "IntegralFamily.sector_decompose",
                    }:
                        entry["reason"] = (
                            "Requires a host registering hepkit.sector_decompose, absent from the pinned browser assembly."
                        )
                    if qualified in {"Tensor.to_numpy", "Tensor.from_numpy"}:
                        entry["packages"] = ["numpy"]
                    entries.append(entry)
                    prior += "\n" + source
    entries = list({entry["id"]: entry for entry in entries}.values())
    groups = {}
    for entry in entries:
        groups.setdefault((entry["product"], entry["owner"]), []).append(entry)
    for (product, owner), group in groups.items():
        for index, entry in enumerate(group):
            if index:
                entry["depends"] = [group[index - 1]["id"]]
                entry["setup"] = ""
            elif "." in owner:
                parent = groups.get((product, owner.split(".")[0]), [])
                if parent:
                    parent_source = "\n".join(cell["source"] for cell in parent)
                    # Docstring examples carry one trailing newline per statement.
                    if entry["setup"].strip().startswith(parent_source.strip()):
                        entry["depends"] = [parent[-1]["id"]]
                        entry["setup"] = (
                            entry["setup"]
                            .strip()[len(parent_source.strip()) :]
                            .lstrip()
                        )
            for dependency in entry["depends"]:
                prerequisite = next(
                    cell for cell in entries if cell["id"] == dependency
                )
                if prerequisite.get("reason"):
                    entry["reason"] = prerequisite["reason"]
    return entries


def validate_contexts(entries, contexts):
    ids = {cell["id"] for cell in entries}
    if ids != set(contexts):
        raise ValueError(f"cell metadata inventory drift: {ids ^ set(contexts)}")
    visiting, complete = set(), set()

    def visit(identity):
        if identity in visiting:
            raise ValueError(f"cyclic documentation prerequisite: {identity}")
        if identity in complete:
            return
        if identity not in ids:
            raise ValueError(f"missing documentation prerequisite: {identity}")
        visiting.add(identity)
        context = contexts[identity]
        if context.get("reason") and context.get("depends"):
            raise ValueError(f"excluded cell has runnable prerequisites: {identity}")
        for dependency in context.get("depends", []):
            if contexts.get(dependency, {}).get("reason"):
                raise ValueError(f"excluded prerequisite: {dependency}")
            visit(dependency)
        visiting.remove(identity)
        complete.add(identity)

    for identity in ids:
        visit(identity)


def typst_reference_inventory(root):
    entries = []
    for package, product in (("linnest", "linnet"), ("kurvst", "gammaloop")):
        for path in sorted((root / f"crates/{package}/typst/src").glob("*.typ")):
            comments = "\n".join(
                re.sub(r"^\s*/// ?", "", line) if re.match(r"^\s*///", line) else ""
                for line in path.read_text().splitlines()
            )
            for index, match in enumerate(
                re.finditer(r"(`{3,})example[^\n]*\n([\s\S]*?)\1", comments)
            ):
                lines = match[2].rstrip().splitlines()
                displayed = "\n".join(
                    line.removeprefix("<<<")
                    for line in lines
                    if not line.startswith(">>>")
                )
                execution = "\n".join(
                    line.removeprefix(">>>")
                    for line in lines
                    if not line.startswith("<<<")
                )
                prefix = f'#import "/crates/{package}/typst/docs/{"example-scope" if package == "linnest" else "manual"}.typ": *\n'
                if package == "linnest" and path.stem == "graph":
                    prefix += "#import graph: *\n"
                entries.append(
                    {
                        "id": f"{path.relative_to(root)}:{index}",
                        "file": str(path.relative_to(root)),
                        "language": "typst",
                        "source": displayed,
                        "execution": execution,
                        "prefix": prefix,
                        "product": product,
                        "sha256": hashlib.sha256(displayed.encode()).hexdigest(),
                        "executionSha256": hashlib.sha256(
                            execution.encode()
                        ).hexdigest(),
                        "depends": [],
                    }
                )
    return entries


def content_manifest(directory):
    files = {
        path.relative_to(directory).as_posix(): hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        for path in sorted(directory.rglob("*"))
        if path.is_file() and path.relative_to(directory).as_posix() != "manifest.json"
    }
    identity = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
    return {"identity": identity, "files": files}


def export(
    root,
    output,
    compiler,
    renderer,
    pyodide,
    ty,
    wheel,
    symbolica_stub,
    typst_packages,
    font_path,
    default_fonts,
    linnest_wasm,
):
    output.mkdir(parents=True, exist_ok=True)
    licenses = output / "licenses"
    licenses.mkdir(exist_ok=True)
    for directory in ("crates/linnest/typst", "crates/kurvst/typst"):
        for notice in (root / directory).glob("LICENSE*"):
            shutil.copyfile(
                notice, licenses / (Path(directory).parent.name + "-" + notice.name)
            )
    shutil.copyfile(
        root / "crates/linnet-render-py/fonts/STIX-Two-OFL.txt",
        licenses / "STIX-Two-OFL.txt",
    )
    for name, source in (("compiler", compiler), ("renderer", renderer)):
        destination = output / name
        destination.mkdir(exist_ok=True)
        for path in (source / "pkg").iterdir():
            if path.suffix in (".mjs", ".wasm"):
                shutil.copyfile(path, destination / path.name)
    shutil.copytree(pyodide, output / "pyodide", dirs_exist_ok=True)
    shutil.copytree(ty, output / "ty", dirs_exist_ok=True)
    shutil.copyfile(
        root / "crates/linnet-render-py/src/view.js", output / "graph-measure.js"
    )
    for name in (
        "documentation-cells.js",
        "documentation-editor.js",
        "code-blocks.js",
        "documentation-python-worker.js",
        "documentation-typst-worker.js",
        "documentation-ty-worker.js",
    ):
        shutil.copyfile(root / "docs/assets" / name, output / name)
    shutil.copyfile(
        root / "crates/linnet-py/examples/notebook-cache-worker.js",
        output / "cache-worker.js",
    )
    files = {}
    for directory in (
        "crates/linnest/typst",
        "crates/kurvst/typst",
        "vendor/typst-packages",
    ):
        for path in sorted((root / directory).rglob("*")):
            if path.is_file() and (
                path.suffix in (".typ", ".toml", ".wasm")
                or path.name.startswith(("LICENSE", "NOTICE"))
            ):
                name = path.relative_to(root).as_posix()
                name = name.replace("vendor/typst-packages/", "packages/")
                files[name] = base64.b64encode(path.read_bytes()).decode()
    for path in sorted(typst_packages.rglob("*", recurse_symlinks=True)):
        if path.is_file():
            files["packages/" + path.relative_to(typst_packages).as_posix()] = (
                base64.b64encode(path.read_bytes()).decode()
            )
    for source, target in (
        ("linnest.wasm", "crates/linnest/typst/linnest.wasm"),
        ("ec-layout.wasm", "crates/linnest/typst/ec-layout.wasm"),
        ("kurvst.wasm", "crates/kurvst/typst/kurvst.wasm"),
    ):
        files[target] = base64.b64encode((linnest_wasm / source).read_bytes()).decode()
    (output / "typst-files.json").write_text(json.dumps(files, separators=(",", ":")))
    fonts = output / "fonts"
    fonts.mkdir(exist_ok=True)
    shutil.copyfile(
        root / "crates/linnet-render-py/fonts/STIXTwoMath-Regular.woff2",
        fonts / "STIXTwoMath-Regular.woff2",
    )
    for path in sorted(list(font_path.rglob("*")) + list(default_fonts.rglob("*"))):
        if path.is_file() and path.suffix in (".ttf", ".otf", ".woff2"):
            shutil.copyfile(path, fonts / path.name)
    (output / "fonts.json").write_text(
        json.dumps(["fonts/" + path.name for path in sorted(fonts.iterdir())])
    )
    shutil.copyfile(
        root / "docs/documentation-cell-checks.json", output / "checks.json"
    )
    # The wheel supplies the exact Symbolica/community module tree used to run
    # examples. Local generated stubs supersede only their matching modules.
    stubs = {
        "symbolica/__init__.pyi": "from .core import *\n",
        "symbolica/core.pyi": symbolica_stub.read_text(),
    }
    with zipfile.ZipFile(wheel) as archive:
        for name in archive.namelist():
            if name.endswith(".pyi"):
                stubs[name] = archive.read(name).decode()
    for source, target in (
        ("docs/api/python/linnet-graph.pyi", "symbolica/community/graph/__init__.pyi"),
        (
            "docs/api/python/feynkit-community.pyi",
            "symbolica/community/hepkit/__init__.pyi",
        ),
        ("docs/api/python/spynso3.pyi", "symbolica/community/tensor/__init__.pyi"),
    ):
        stubs[target] = (root / source).read_text()
    (output / "type-stubs.json").write_text(json.dumps(stubs))
    entries = inventory(root)
    python_files = {}
    for cell in entries:
        for name in re.findall(r'["\'](crates/[^"\']+\.json)["\']', cell["source"]):
            path = root / name
            if path.is_file():
                python_files[name] = path.read_text()
    (output / "python-files.json").write_text(json.dumps(python_files))
    contexts = tomllib.loads((root / "docs/documentation-cells.toml").read_text())[
        "cell"
    ]
    validate_contexts(entries, contexts)
    for cell in entries:
        cell.update(contexts[cell["id"]])
    entries.extend(python_reference_inventory(root))
    entries.extend(typst_reference_inventory(root))
    (output / "cells.json").write_text(json.dumps(entries, indent=2))
    (output / "manifest.json").write_text(json.dumps(content_manifest(output)))


def publish(site, runtime, bundle):
    manifest = json.loads((bundle / "manifest.json").read_text())
    if manifest != content_manifest(bundle):
        raise ValueError("documentation runtime content does not match its manifest")
    identity = manifest["identity"]
    runtime_manifest = json.loads((runtime / "runtime-manifest.json").read_text())
    runtime_id = runtime_manifest["identity"]
    shared_files = runtime_manifest["shared_files"]
    wheel_name = next((runtime / "wheels").glob("symbolica-*.whl")).name
    target = site / "assets/documentation-runtime" / identity
    if target.exists() and (
        content_manifest(target) != manifest
        or not (target / "manifest.json").is_file()
        or (target / "manifest.json").read_bytes()
        != (bundle / "manifest.json").read_bytes()
    ):
        raise ValueError(f"immutable documentation runtime collision: {identity}")
    runtime_target = site / "assets/notebook-runtime" / runtime_id
    assert (runtime_target / "native-loader.js").is_file(), (
        "publish notebooks before standalone cells"
    )
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=target.parent) as staging:
            staged = Path(staging) / identity
            shutil.copytree(bundle, staged)
            staged.chmod(0o755)
            staged.rename(target)
    cache_source = (bundle / "cache-worker.js").read_bytes()
    cache_worker = site / (
        "notebook-cache-worker-" + hashlib.sha256(cache_source).hexdigest() + ".js"
    )
    if not cache_worker.exists():
        cache_worker.write_bytes(cache_source)
    files = json.loads((bundle / "python-files.json").read_text())
    for assets in site.rglob("assets"):
        if not (assets / "site.js").is_file():
            continue
        config = {
            "compilerRoot": os.path.relpath(target, assets) + "/",
            "entry": os.path.relpath(target / "documentation-cells.js", assets),
            "cacheWorker": os.path.relpath(cache_worker, assets),
            "runtime": os.path.relpath(runtime_target, assets) + "/",
            "wheel": f"../../notebook-runtime-files/{shared_files['wheels/' + wheel_name]}/{wheel_name}",
            "stdlib": f"../../notebook-runtime-files/{shared_files['python_stdlib.zip']}/python_stdlib.zip",
            "fonts": json.loads((bundle / "fonts.json").read_text()),
            "files": files,
        }
        with tempfile.NamedTemporaryFile(mode="w", dir=assets, delete=False) as stream:
            stream.write(json.dumps(config))
            temporary = stream.name
        os.replace(temporary, assets / "documentation-cells-config.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--inventory", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--site", type=Path)
    parser.add_argument("--runtime", type=Path)
    parser.add_argument("--bundle", type=Path)
    for name in (
        "compiler",
        "renderer",
        "pyodide",
        "ty",
        "wheel",
        "symbolica-stub",
        "typst-packages",
        "font-path",
        "default-fonts",
        "linnest-wasm",
    ):
        parser.add_argument(f"--{name}", type=Path)
    args = parser.parse_args()
    if args.verify:
        manifest = json.loads((args.bundle / "manifest.json").read_text())
        if manifest != content_manifest(args.bundle):
            raise ValueError(
                "documentation runtime content does not match its manifest"
            )
        print(manifest["identity"])
    elif args.site:
        publish(args.site, args.runtime, args.bundle)
    elif args.inventory:
        print(json.dumps(inventory(args.root), indent=2))
    else:
        export(
            args.root,
            args.output,
            args.compiler,
            args.renderer,
            args.pyodide,
            args.ty,
            args.wheel,
            args.symbolica_stub,
            args.typst_packages,
            args.font_path,
            args.default_fonts,
            args.linnest_wasm,
        )


if __name__ == "__main__":
    main()
