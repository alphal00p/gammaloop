"""Prepare the pinned islands worker and split wheel without rebuilding Rust."""

import base64
import compileall
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import py_compile
import sys
import tempfile
import zipfile
from pathlib import Path


def cache_python_archive(source: Path, destination: Path, *, legacy: bool = False):
    """Add validated bytecode while retaining every upstream source and data file."""
    with tempfile.TemporaryDirectory() as temporary:
        stage = Path(temporary)
        with zipfile.ZipFile(source) as archive:
            archive.extractall(stage)
        compiled = compileall.compile_dir(
            stage,
            quiet=2,
            force=True,
            stripdir=str(stage),
            legacy=legacy,
            prependdir="/lib/python314.zip"
            if legacy
            else "/lib/python3.14/site-packages",
            invalidation_mode=py_compile.PycInvalidationMode.CHECKED_HASH,
        )
        if not compiled:
            raise ValueError(f"Failed to prepare runtime bytecode: {source.name}")
        record_name = next(
            (
                path.relative_to(stage).as_posix()
                for path in stage.glob("*.dist-info/RECORD")
            ),
            None,
        )
        rows = []
        with zipfile.ZipFile(
            destination, "w", zipfile.ZIP_DEFLATED, compresslevel=9
        ) as archive:
            for path in sorted(stage.rglob("*")):
                if not path.is_file():
                    continue
                name = path.relative_to(stage).as_posix()
                if name == record_name:
                    continue
                data = path.read_bytes()
                digest = (
                    base64.urlsafe_b64encode(hashlib.sha256(data).digest())
                    .decode()
                    .rstrip("=")
                )
                rows.append([name, "sha256=" + digest, str(len(data))])
                info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                archive.writestr(info, data)
            record = io.StringIO()
            if record_name:
                csv.writer(record, lineterminator="\n").writerows(
                    [*rows, [record_name, "", ""]]
                )
                info = zipfile.ZipInfo(record_name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                archive.writestr(info, record.getvalue())


def main():
    root, wheel, runtime_lock, stdlib, *runtime_wheels = map(Path, sys.argv[1:])
    wheels = root / "wheels"
    wheels.mkdir()
    with zipfile.ZipFile(wheel) as source:
        native_path = "symbolica/core.abi3.so"
        record_path = "symbolica-3.0.1.dist-info/RECORD"
        binary = source.read(native_path)
        record = source.read(record_path).decode()
        entry = next(
            row for row in record.splitlines(True) if row.startswith(native_path + ",")
        )
        with zipfile.ZipFile(
            wheels / wheel.name, "w", zipfile.ZIP_DEFLATED, compresslevel=9
        ) as support:
            for name in source.namelist():
                if name == native_path:
                    continue
                data = (
                    record.replace(entry, "").encode()
                    if name == record_path
                    else source.read(name)
                )
                info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                support.writestr(info, data)
    (root / "symbolica-core.wasm").write_bytes(binary)
    (root / "symbolica-core.wasm.gz").write_bytes(
        gzip.compress(binary, compresslevel=9, mtime=0)
    )
    descriptor = {
        "path": native_path,
        "size": len(binary),
        "record": record_path,
        "entry": entry,
        "bytecode_magic": importlib.util.MAGIC_NUMBER.hex(),
        "sha256": hashlib.sha256(binary).hexdigest(),
    }

    # URLs resolve from assets/notebook-runtime/<runtime hash>/. The exporter
    # verifies and publishes each shared file under its own SHA-256.
    def shared_url(path, *, worker=False):
        digest = hashlib.sha256((root / path).read_bytes()).hexdigest()
        return ("../../../" if worker else "../../") + (
            f"notebook-runtime-files/{digest}/{Path(path).name}"
        )

    loader = Path(__file__).with_name("native-loader.js").read_text()
    (root / "native-loader.js").write_text(
        loader.replace("__NATIVE_DESCRIPTOR__", json.dumps(descriptor)).replace(
            "__NATIVE_URL__", shared_url("symbolica-core.wasm")
        )
    )
    # Keep every source and data file; add checked CPython bytecode so fresh
    # interpreters do not recompile the same pinned runtime on each navigation.
    lock = json.loads(runtime_lock.read_text())
    for package_wheel in runtime_wheels:
        original_digest = hashlib.sha256(package_wheel.read_bytes()).hexdigest()
        packages = [
            package
            for package in lock["packages"].values()
            if package["sha256"] == original_digest
        ]
        filename = packages[0]["file_name"].rsplit("/", 1)[-1]
        cache_python_archive(package_wheel, wheels / filename)
        prepared_digest = hashlib.sha256((wheels / filename).read_bytes()).hexdigest()
        # Marimo's lock has both `marimo` and `marimo-base` aliases for one
        # archive. Every alias must use the same prepared offline bytes.
        for package in packages:
            package["file_name"] = f"__LOCAL_WHEEL__/{prepared_digest}/{filename}"
            package["sha256"] = prepared_digest
    (root / "pyodide-lock.json").write_text(json.dumps(lock))
    cache_python_archive(stdlib, root / "python_stdlib.zip", legacy=True)
    worker = root / "assets/worker-D7hWi02M.js"
    text = worker.read_text()
    replacements = {
        "indexURL:n,packageBaseUrl:n,convertNullToNone:!0": f"indexURL:n,stdLibURL:new URL('{shared_url('python_stdlib.zip', worker=True)}',import.meta.url).href,packageBaseUrl:n,convertNullToNone:!0",
        "lockFileURL:e.pyodideLockfileUrl??`https://wasm.marimo.app/pyodide-lock.json?v=${e.version}&pyodide=${e.pyodideVersion}`": "lockFileContents:await (async()=>{const lock=await (await fetch(new URL('../pyodide-lock.json',import.meta.url))).json();for(const p of Object.values(lock.packages)){if(p.file_name.startsWith('__LOCAL_WHEEL__/'))p.file_name=new URL('../../../notebook-runtime-files/'+p.file_name.slice(16),import.meta.url).href;}return lock;})()",
        "`narwhals`,`packaging`]": "`narwhals`,`packaging`,`pygments`,`pyodide-http`,`jedi`,`docutils`,`anywidget`,new URL('"
        + shared_url("wheels/" + wheel.name, worker=True)
        + "',import.meta.url).href]",
        "return this.pyodide=r,t.end(`ok`),r": "await mountNative(r);return this.pyodide=r,t.end(`ok`),r",
    }
    for old, new in replacements.items():
        if text.count(old) != 1:
            raise ValueError("Review the pinned Marimo worker boundary: " + old)
        text = text.replace(old, new)
    worker.write_text('import {mountNative} from "../native-loader.js";\n' + text)
    # The upstream 20-second RPC limit includes downloading and initializing
    # Python. On a slow connection it discards the active session before the
    # kernel responds, so Python runs while every island remains unmounted.
    # Keep a bounded five-minute startup budget, and let the host await/report
    # initialization instead of upstream's fire-and-forget console-only error.
    # Retain self-starting behavior for previously cached page JavaScript.
    main = root / "main.js"
    text = main.read_text()
    replacements = {
        "maxRequestTime: 2e4,": "maxRequestTime: 300000,",
        "RB().catch((e) => {": "const initialized = RB();\nexport { initialized };\ninitialized.catch((e) => {",
    }
    for old, new in replacements.items():
        if text.count(old) != 1:
            raise ValueError("Review the pinned Marimo initialization boundary: " + old)
        text = text.replace(old, new)
    main.write_text(text)
    manifest = {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }
    encoded = json.dumps(manifest, sort_keys=True, separators=(",", ":"))
    identity = hashlib.sha256(encoded.encode()).hexdigest()
    (root / "runtime-manifest.json").write_text(
        json.dumps(
            {
                "identity": identity,
                "files": manifest,
                "shared_files": {
                    name: digest
                    for name, digest in manifest.items()
                    if name == "python_stdlib.zip"
                    or name.startswith("wheels/")
                    or name == "symbolica-core.wasm"
                },
            }
        )
    )


if __name__ == "__main__":
    main()
