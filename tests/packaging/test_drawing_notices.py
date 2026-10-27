"""Exercise Maturin's shipped license payload without rebuilding the physics stack."""

import base64
import csv
import hashlib
import io
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import tomllib
import unittest
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[2]
HOST = ROOT / "examples/notebooks/symbolica-host"
LINNEST_LICENSE = "MIT AND Apache-2.0 AND BSD-3-Clause"
SOURCES = {
    "clarabel-LICENSE": ROOT / "vendor/clarabel/LICENSE.md",
    "networkx-LICENSE": ROOT / "crates/linnest/layout/native/LICENSE.networkx",
    "tiptoe-LICENSE": ROOT / "crates/kurvst/typst/LICENSE.tiptoe",
}
PACKAGES = {"clarabel": "linnest", "networkx": "linnest", "tiptoe": "kurvst"}
EMBED_MANIFEST = (
    '[package]\nname = "drawing-embed-fixture"\nversion = "0.1.0"\n'
    'edition = "2021"\n[workspace]\n[dependencies]\n'
    'include_dir = "0.7.4"\ncolor-eyre = "0.6"\n'
    '[build-dependencies]\nwalkdir = "2"\n'
)
WHEEL_MANIFEST = (
    '[package]\nname = "drawing-notice-fixture"\nversion = "0.1.0"\n'
    'edition = "2021"\n[workspace]\n[lib]\nname = "core"\n'
    'crate-type = ["cdylib"]\n[dependencies]\n'
    'pyo3 = { version = "0.28", features = ["abi3-py310", "extension-module"] }\n'
)


class DrawingNotices(unittest.TestCase):
    @staticmethod
    def prepare_offline_fixtures():
        """Fetch the same manifests the offline tests resolve on a clean runner."""
        scratch = ROOT / "target/notice-validation"
        scratch.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as temporary:
            for name, manifest in [
                ("embed", EMBED_MANIFEST),
                ("wheel", WHEEL_MANIFEST),
            ]:
                fixture = Path(temporary) / name
                (fixture / "src").mkdir(parents=True)
                (fixture / "src/lib.rs").write_text("")
                (fixture / "Cargo.toml").write_text(manifest)
                subprocess.run(
                    ["cargo", "fetch", "--manifest-path", str(fixture / "Cargo.toml")],
                    cwd=ROOT,
                    check=True,
                )

    def test_development_sync_upgrades_notice_payload(self):
        # Execute each recipe's actual source-copy helpers, avoiding unrelated
        # WASM rebuilds. Start with an older extracted bundle without Tiptoe.
        recipes = (ROOT / "just/drawing.just").read_text()
        scratch = ROOT / "target/notice-validation"
        scratch.mkdir(parents=True, exist_ok=True)
        for recipe in ["sync-drawing-assets", "sync-drawing-assets-cargo"]:
            with (
                self.subTest(recipe=recipe),
                tempfile.TemporaryDirectory(dir=scratch) as temporary,
            ):
                target = Path(temporary)
                for package in ["linnest", "kurvst"]:
                    destination = target / f"crates/{package}/typst"
                    shutil.copytree(ROOT / f"crates/{package}/typst", destination)
                tiptoe = target / "crates/kurvst/typst/LICENSE.tiptoe"
                tiptoe.unlink()
                self.assertFalse(tiptoe.exists())
                body = recipes.split(recipe + ":\n", 1)[1].split("\n#", 1)[0]
                helpers = re.findall(r"^    \w+\(\) \{\n.*?^    \}", body, re.M | re.S)
                self.assertEqual(len(helpers), 6)
                subprocess.run(
                    [
                        "bash",
                        "-euc",
                        "\n".join(helpers) + '\ncopy_package_sources "$1"',
                        recipe,
                        str(target),
                    ],
                    cwd=ROOT,
                    check=True,
                )
                for package in ["linnest", "kurvst"]:
                    source = ROOT / f"crates/{package}/typst"
                    for notice in source.glob("LICENSE*"):
                        self.assertEqual(
                            (
                                target / f"crates/{package}/typst" / notice.name
                            ).read_bytes(),
                            notice.read_bytes(),
                        )
                self.assertEqual(
                    (target / "crates/kurvst/typst/src/mark.typ").read_bytes(),
                    (ROOT / "crates/kurvst/typst/src/mark.typ").read_bytes(),
                )

    def test_gammaloop_sdist_payload(self):
        scratch = ROOT / "target/notice-validation"
        scratch.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as temporary:
            result = subprocess.run(
                [
                    "maturin",
                    "sdist",
                    "--manifest-path",
                    "crates/gammaloop-api/Cargo.toml",
                    "--out",
                    temporary,
                ],
                cwd=ROOT,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            archives = list(Path(temporary).glob("*.tar.gz"))
            self.assertEqual(len(archives), 1, archives)
            with tarfile.open(archives[0]) as archive:
                prefix = archive.getnames()[0].split("/")[0] + "/"
                metadata = archive.extractfile(
                    prefix + "crates/linnest/typst/typst.toml"
                ).read()
                self.assertEqual(
                    tomllib.loads(metadata.decode())["package"]["license"],
                    LINNEST_LICENSE,
                )
                for relative in [
                    "crates/linnest/typst/LICENSE.clarabel",
                    "crates/linnest/typst/LICENSE.networkx",
                    "crates/linnest/typst/LICENSE",
                    "crates/kurvst/typst/LICENSE",
                    "crates/kurvst/typst/LICENSE.tiptoe",
                ]:
                    with self.subTest(path=relative):
                        self.assertEqual(
                            archive.extractfile(prefix + relative).read(),
                            (ROOT / relative).read_bytes(),
                        )

    def test_api_embedded_payload(self):
        # Compile the actual API build script and extraction module in a small
        # Cargo fixture. This exercises their shipped bytes without the unrelated
        # physics tests and dependencies of the full API crate.
        scratch = ROOT / "target/notice-validation"
        scratch.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as temporary:
            fixture = Path(temporary)
            (fixture / "assets").symlink_to(ROOT / "assets", target_is_directory=True)
            crates = fixture / "crates"
            crates.mkdir()
            for package in ["linnest", "kurvst"]:
                (crates / package).symlink_to(
                    ROOT / "crates" / package, target_is_directory=True
                )
            api = crates / "gammaloop-api"
            (api / "src").mkdir(parents=True)
            shutil.copy(ROOT / "crates/gammaloop-api/build.rs", api / "build.rs")
            (api / "src/templates.rs").symlink_to(
                ROOT / "crates/gammaloop-api/src/templates.rs"
            )
            (api / "Cargo.toml").write_text(EMBED_MANIFEST)
            (api / "src/main.rs").write_text(
                "mod templates;\nfn main() -> color_eyre::Result<()> {\n"
                'templates::Assets::extract_templates("extracted")\n}\n'
            )
            result = subprocess.run(
                ["cargo", "run", "--offline", "--quiet"],
                cwd=api,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            extracted = api / "extracted/drawings/templates/crates"
            self.assertEqual(
                tomllib.loads((extracted / "linnest/typst/typst.toml").read_text())[
                    "package"
                ]["license"],
                LINNEST_LICENSE,
            )
            for name, source in SOURCES.items():
                suffix = name.removesuffix("-LICENSE")
                self.assertEqual(
                    (
                        extracted / f"{PACKAGES[suffix]}/typst/LICENSE.{suffix}"
                    ).read_bytes(),
                    source.read_bytes(),
                )
            for package in ["linnest", "kurvst"]:
                self.assertEqual(
                    (extracted / package / "typst/LICENSE").read_bytes(),
                    (ROOT / "crates" / package / "typst/LICENSE").read_bytes(),
                )
                self.assertEqual(
                    (extracted / package / f"typst/{package}.wasm").read_bytes(),
                    (ROOT / "crates" / package / f"typst/{package}.wasm").read_bytes(),
                )

    def test_copies_preserve_upstream_bytes(self):
        for name, source in SOURCES.items():
            with self.subTest(name=name):
                self.assertEqual(
                    (HOST / "licenses" / name).read_bytes(), source.read_bytes()
                )
                suffix = name.removesuffix("-LICENSE")
                self.assertEqual(
                    (
                        ROOT / f"crates/{PACKAGES[suffix]}/typst/LICENSE.{suffix}"
                    ).read_bytes(),
                    source.read_bytes(),
                )

    def test_host_wheel_license_payload(self):
        # Use the real host's Python tree, pyproject and licenses. Only the Rust
        # module is a small fixture: notice delivery is Maturin's responsibility
        # and does not depend on compiling Symbolica or loading a user license.
        scratch = ROOT / "target/notice-validation"
        scratch.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as temporary:
            fixture = Path(temporary)
            shutil.copy(HOST / "pyproject.toml", fixture)
            shutil.copytree(HOST / "python", fixture / "python")
            shutil.copytree(HOST / "licenses", fixture / "licenses")
            (fixture / "Cargo.toml").write_text(WHEEL_MANIFEST)
            (fixture / "src").mkdir()
            (fixture / "src/lib.rs").write_text(
                "use pyo3::{pymodule, Bound, PyResult};\n"
                "use pyo3::types::PyModule;\n"
                "#[pymodule]\nfn core(_module: &Bound<'_, PyModule>) -> PyResult<()> { Ok(()) }\n"
            )
            output = Path(
                os.environ.get("DRAWING_NOTICE_ARTIFACTS", fixture / "wheels")
            )
            output = output.resolve()
            result = subprocess.run(
                [
                    "maturin",
                    "build",
                    "--offline",
                    "--compatibility",
                    "linux",
                    "--out",
                    str(output),
                ],
                cwd=fixture,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            self.assertEqual(result.returncode, 0, result.stdout)
            wheels = list(output.glob("symbolica-3.0.1-*.whl"))
            self.assertEqual(len(wheels), 1, wheels)
            with ZipFile(wheels[0]) as wheel:
                prefix = "symbolica-3.0.1.dist-info/"
                metadata = wheel.read(prefix + "METADATA").decode()
                records = {
                    row[0]: row[1:]
                    for row in csv.reader(
                        io.StringIO(wheel.read(prefix + "RECORD").decode())
                    )
                }
                for name, source in SOURCES.items():
                    with self.subTest(name=name):
                        relative = "licenses/" + name
                        path = prefix + "licenses/" + relative
                        contents = wheel.read(path)
                        self.assertEqual(contents, source.read_bytes())
                        self.assertIn("License-File: " + relative, metadata)
                        digest = base64.urlsafe_b64encode(
                            hashlib.sha256(contents).digest()
                        )
                        self.assertEqual(
                            records[path],
                            [
                                "sha256=" + digest.decode().rstrip("="),
                                str(len(contents)),
                            ],
                        )


if __name__ == "__main__":
    if "--prepare-offline-fixtures" in sys.argv:
        sys.argv.remove("--prepare-offline-fixtures")
        DrawingNotices.prepare_offline_fixtures()
    unittest.main()
