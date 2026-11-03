"""Real cross-revision registry reuse, rustdoc, and installed Maturin wheels."""

import json
import shutil
import sys
from pathlib import Path

from test_dyndrv_integration import REPOSITORY, Harness


def qualify(harness):
    h = harness
    other_python = h.options.other_python
    h.prepare()
    for name in ("app", "crates"):
        shutil.rmtree(h.a / name)
    with h.step("registry dependency reused across different workspace sources"):
        h.write(
            h.a,
            "Cargo.toml",
            '[package]\nname="sharing-probe"\nversion="0.1.0"\nedition="2021"\n'
            '[workspace]\n[dependencies]\nitoa="=1.0.18"\n'
            '[profile.ci-optim]\ninherits="dev"\n',
        )
        h.write(
            h.a,
            "src/main.rs",
            'fn main() { println!("{}", itoa::Buffer::new().format(42)); }\n',
        )
        h.write(
            h.a, "src/lib.rs", "/// A documented value.\npub fn value() -> u32 { 42 }\n"
        )
        h.lock_and_metadata(h.a, "registry")
        _, first = h.public(h.a, "registry-first", "build", "-p", "sharing-probe")
        shutil.copytree(
            h.a,
            h.b,
            ignore=shutil.ignore_patterns("target", ".git", ".venv", "result-*"),
        )
        h.command("second-git", ["git", "init", "--quiet", h.b])
        h.write(h.b, "README.md", "Unrelated checkout content differs.\n")
        h.write(
            h.b,
            "src/main.rs",
            'fn main() { println!("{}", itoa::Buffer::new().format(43)); }\n',
        )
        _, second = h.public(h.b, "registry-second", "build", "-p", "sharing-probe")

        a = [r for r in first["receipts"] if r["crate"] == "itoa"]
        b = [r for r in second["receipts"] if r["crate"] == "itoa"]
        h.require(a and b, "Missing real registry compiler receipts")
        h.require(
            [(r["derivation"], r["output"]) for r in a]
            == [(r["derivation"], r["output"]) for r in b],
            "Unrelated workspace edits changed the registry compiler request",
        )
        for receipt in b:
            log = Path(
                receipt["path"].replace(".nix-result.json", ".nix-stderr")
            ).read_text()
            for line in log.splitlines():
                if line.startswith("@nix "):
                    event = json.loads(line[5:])
                    h.require(
                        not (
                            event.get("action") == "start" and event.get("type") == 105
                        ),
                        "Registry reuse started a compiler builder",
                    )
        for workspace, expected in ((h.a, "42"), (h.b, "43")):
            output, _, _ = h.command(
                "binary-" + expected, [workspace / "target/ci-optim/sharing-probe"]
            )
            h.require(
                output.strip() == expected, "Workspace output was shared incorrectly"
            )
    with h.step("native rustdoc exports HTML and preserves failure status"):
        h.public(h.a, "docs", "doc", "-p", "sharing-probe", "--no-deps")
        h.require(
            (h.a / "target/doc/sharing_probe/index.html").is_file(),
            "Missing exported documentation",
        )
        h.write(h.a, "src/lib.rs", "/// [missing_item]\npub fn value() -> u32 { 42 }\n")
        h.public(
            h.a,
            "docs-failure",
            "doc",
            "-p",
            "sharing-probe",
            "--no-deps",
            env=h.environment | {"RUSTDOCFLAGS": "-D rustdoc::broken_intra_doc_links"},
            expected=101,
        )
        h.write(h.a, "src/lib.rs", "/// Repaired.\npub fn value() -> u32 { 42 }\n")
        h.public(h.a, "docs-recovery", "doc", "-p", "sharing-probe", "--no-deps")
    with h.step("Python wheel uses immutable compilers and imports outside namespace"):
        h.write(
            h.a,
            "scripts/python_wheel.py",
            (REPOSITORY / "scripts/python_wheel.py").read_text(),
        )
        h.write(
            h.a,
            "Cargo.toml",
            '[workspace]\nmembers=["crates/gammaloop-api"]\nresolver="2"\n[profile.ci-optim]\ninherits="dev"\n',
        )
        h.write(
            h.a,
            "crates/gammaloop-api/Cargo.toml",
            '[package]\nname="gammaloop-api"\nversion="0.0.0"\nedition="2021"\n'
            '[lib]\nname="_gammaloop"\ncrate-type=["cdylib"]\n'
            '[features]\nstandalone=[]\npyo3-extension-module=["pyo3/extension-module"]\n'
            '[dependencies]\npyo3="=0.28.3"\nitoa="=1.0.18"\n',
        )
        h.write(
            h.a,
            "crates/gammaloop-api/src/lib.rs",
            "use pyo3::prelude::*;\n#[pyfunction]\nfn value() -> u32 { itoa::Buffer::new().format(73).parse().unwrap() }\n"
            "#[pymodule]\nfn _gammaloop(m: &Bound<'_, PyModule>) -> PyResult<()> { m.add_function(wrap_pyfunction!(value,m)?) }\n",
        )
        h.write(
            h.a,
            "pyproject.toml",
            '[build-system]\nrequires=["maturin"]\nbuild-backend="maturin"\n'
            '[project]\nname="gammaloop"\nversion="0.0.0"\n'
            '[tool.maturin]\nmodule-name="gammaloop._gammaloop"\npython-source="python"\n',
        )
        h.write(h.a, "python/gammaloop/__init__.py", "from ._gammaloop import value\n")
        h.lock_and_metadata(h.a, "python")
        signature, _, _ = h.command(
            "python-signature",
            [
                h.spec["environment"]["PYO3_PYTHON"],
                "-c",
                "import sys,struct; print(f'cpython-{sys.version_info.major}.{sys.version_info.minor}-{struct.calcsize(chr(80))*8}bit')",
            ],
        )
        # Compare equivalent compiler environments, not ordinary non-Python
        # builds: Maturin supplies these compile-affecting PyO3 variables.
        h.write(
            h.a,
            ".cargo/config.toml",
            "[env]\n"
            + "PYTHON_SYS_EXECUTABLE="
            + json.dumps(h.spec["environment"]["PYO3_PYTHON"])
            + "\n"
            + "PYO3_ENVIRONMENT_SIGNATURE="
            + json.dumps(signature.strip())
            + '\nPYO3_BUILD_EXTENSION_MODULE="1"\n',
        )
        # Stable vendor paths can leave itoa Cargo-fresh from the earlier
        # registry fixture. Require a delivery on both sides of this comparison.
        h.public(h.a, "native-extension-cold-target", "clean")
        _, native_record = h.public(
            h.a,
            "native-extension-build",
            "build",
            "-p",
            "gammaloop-api",
            "--no-default-features",
            "--features",
            "standalone,pyo3-extension-module",
        )
        h.public(h.a, "native-extension-clean", "clean")
        _, python_record = h.public(
            h.a,
            "python-wheel",
            "python",
            "gammaloop",
            "--python",
            h.spec["environment"]["PYO3_PYTHON"],
        )
        python_itoa = [r for r in python_record["receipts"] if r["crate"] == "itoa"]
        h.require(
            [(r["derivation"], r["output"]) for r in python_itoa]
            == [
                (r["derivation"], r["output"])
                for r in native_record["receipts"]
                if r["crate"] == "itoa"
            ],
            "Python changed the matching native extension dependency request",
        )
        wheels = [
            path for path in (h.a / "target/wheels").glob("*.whl") if path.is_file()
        ]
        h.require(len(wheels) == 1, "Missing exported wheel")
        root = h.a / "target/wheels" / (".nix-" + wheels[0].name)
        h.require(
            root.is_symlink(), "Wheel runtime payload lacks an independent GC root"
        )
        closure, _, _ = h.command(
            "wheel-runtime-closure",
            ["nix-store", "--query", "--requisites", root.resolve()],
        )
        required = {str(h.bootstrap)} | {
            str(Path(*Path(path).parts[:4]))
            for path in h.spec["environment"]["LD_LIBRARY_PATH"].split(":")
        }
        h.require(
            required.issubset(set(closure.splitlines())),
            "Wheel root does not retain its runtime closure",
        )
        venv = h.run_directory / "wheel-environment"
        h.command(
            "wheel-venv",
            ["uv", "venv", "--python", h.spec["environment"]["PYO3_PYTHON"], venv],
        )
        h.command(
            "wheel-install",
            [
                "uv",
                "pip",
                "install",
                "--python",
                venv / "bin/python",
                "--no-deps",
                wheels[0],
            ],
        )
        h.public(h.a, "python-clean", "clean")
        h.public(h.a, "python-private-gc", "gc")
        h.command(
            "wheel-import-after-clean",
            [
                venv / "bin/python",
                "-c",
                "import gammaloop; assert gammaloop.value() == 73",
            ],
            env={
                k: v
                for k, v in h.environment.items()
                if k not in {"LD_LIBRARY_PATH", "PYTHONPATH"}
            },
        )

    if other_python:
        with h.step("matching bootstrap for another CPython version"):
            h.command(
                "other-python-wheel",
                [
                    "bash",
                    h.a / "scripts/develop.sh",
                    "--backend=dyndrv",
                    "python",
                    "gammaloop",
                    "--python",
                    other_python,
                    "--out",
                    "target/wheels-other",
                ],
            )
            wheels = [
                path
                for path in (h.a / "target/wheels-other").glob("*.whl")
                if path.is_file()
            ]
            h.require(len(wheels) == 1, "Missing wheel for the other interpreter")
            runtime = (wheels[0].parent / (".nix-" + wheels[0].name)).resolve()
            closure, _, _ = h.command(
                "other-runtime-closure",
                ["nix-store", "--query", "--requisites", runtime],
            )
            selected = (runtime / "bootstrap").resolve()
            selected_spec = json.loads((selected / "environment.json").read_text())
            required = {str(selected)} | {
                str(Path(*Path(path).parts[:4]))
                for path in selected_spec["environment"]["LD_LIBRARY_PATH"].split(":")
            }
            h.require(
                required.issubset(set(closure.splitlines())),
                "Alternate Python bootstrap is not retained",
            )
            venv = h.run_directory / "other-wheel-environment"
            h.command(
                "other-wheel-venv", ["uv", "venv", "--python", other_python, venv]
            )
            h.command(
                "other-wheel-install",
                [
                    "uv",
                    "pip",
                    "install",
                    "--python",
                    venv / "bin/python",
                    "--no-deps",
                    wheels[0],
                ],
            )
            h.command(
                "other-wheel-import",
                [
                    venv / "bin/python",
                    "-c",
                    "import gammaloop; assert gammaloop.value() == 73",
                ],
            )


if __name__ == "__main__":
    parser = Harness.parser(__doc__, "extensions")
    parser.add_argument("--other-python", type=Path)
    sys.exit(Harness.run_cli(parser, qualify))
