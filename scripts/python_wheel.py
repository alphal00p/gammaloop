#!/usr/bin/env python3
"""Build and install native Python distribution wheels."""

import argparse
import fcntl
import json
import os
import platform
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from urllib.parse import quote

ROOT = Path(os.environ.get("GAMMALOOP_WORKSPACE_ROOT", Path.cwd())).resolve()
HOST = "alphal00p-notebook-host"


def nix_environment(env=None):
    """Keep runtime licenses and stale request arguments out of Nix processes."""
    result = dict(os.environ if env is None else env)
    for key in (
        "SYMBOLICA_LICENSE",
        "SYMBOLICA_LICENSE_SIGNED",
        "GAMMALOOP_NIX_TEST_ARGS",
        "GAMMALOOP_DEVELOP_SPEC",
    ):
        result.pop(key, None)
    return result


def supplied_source():
    """Borrow the snapshot kept alive by the shell backend."""
    supplied = os.environ.get("GAMMALOOP_FLAKE_SOURCE")
    if not supplied:
        raise ValueError(
            "Nix wheels require GAMMALOOP_FLAKE_SOURCE from the shell backend"
        )
    source = Path(supplied).resolve()
    if not source.is_dir() or not (source / "flake.nix").is_file():
        raise ValueError("GAMMALOOP_FLAKE_SOURCE must name a flake snapshot")
    return source


def run(command, *, env=None, capture=False, cwd=ROOT):
    result = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        check=False,
    )
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    return result.stdout.strip() if capture else None


def license_environment():
    env = os.environ.copy()
    if not env.get("SYMBOLICA_LICENSE") and env.get("SYMBOLICA_LICENSE_SIGNED"):
        env["SYMBOLICA_LICENSE"] = env["SYMBOLICA_LICENSE_SIGNED"]
    if not env.get("SYMBOLICA_LICENSE"):
        license_file = (
            Path(env.get("XDG_CONFIG_HOME") or Path.home() / ".config")
            / "symbolica/license"
        )
        if license_file.is_file():
            env["SYMBOLICA_LICENSE"] = license_file.read_text().strip()
    return env


def nix_build(attribute, *, spec=None, env=None, root=None):
    env = nix_environment(env)
    if spec is not None:
        env["GAMMALOOP_DEVELOP_SPEC"] = json.dumps(spec)
    command = ["nix", "build", "--impure", "--print-out-paths", "--no-write-lock-file"]
    if root is None:
        command.append("--no-link")
    else:
        root.parent.mkdir(parents=True, exist_ok=True)
        command += ["--out-link", str(root)]
    command.append("--accept-flake-config")
    source = supplied_source()
    return Path(
        run(
            command + [f"path:{quote(str(source), safe='/')}#{attribute}"],
            env=env,
            capture=True,
        ).splitlines()[-1]
    )


def retain_output(output, root, *, capture=False):
    root.parent.mkdir(parents=True, exist_ok=True)
    run(
        [
            "nix-store",
            "--realise",
            str(output),
            "--add-root",
            str(root),
            "--indirect",
        ],
        env=nix_environment(),
        capture=capture,
    )


def python_environment(requested=None, install=False):
    candidate = requested
    if not candidate:
        prefix = os.environ.get("VIRTUAL_ENV") or os.environ.get("CONDA_PREFIX")
        if prefix:
            candidate = str(Path(prefix) / "bin/python")
        elif (ROOT / ".venv/bin/python").is_file():
            candidate = str(ROOT / ".venv/bin/python")
        else:
            candidate = sys.executable
    executable = shutil.which(candidate) or str(Path(candidate).absolute())
    details = json.loads(
        run(
            [
                executable,
                "-c",
                (
                    "import json,sys,sysconfig,platform; print(json.dumps(dict(version=list(sys.version_info[:2]), "
                    "implementation=sys.implementation.name, prefix=sys.prefix, base=sys.base_prefix, "
                    "machine=platform.machine(), debug=bool(sysconfig.get_config_var('Py_DEBUG')), "
                    "free_threaded=bool(sysconfig.get_config_var('Py_GIL_DISABLED')))))"
                ),
            ],
            capture=True,
        )
    )
    conda_prefix = os.environ.get("CONDA_PREFIX")
    selected_conda = (
        conda_prefix
        and Path(details["prefix"]).resolve() == Path(conda_prefix).resolve()
    )
    if install and details["prefix"] == details["base"] and not selected_conda:
        raise ValueError(
            "Python development requires an existing virtual environment; activate one or create .venv"
        )
    return executable, details


def python_package(arguments, backend):
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("distribution", choices=("gammaloop", "community"))
    parser.add_argument("--install", action="store_true")
    parser.add_argument("--profile")
    parser.add_argument("--release", action="store_true")
    parser.add_argument("--abi3", action="store_true")
    parser.add_argument("--python")
    parser.add_argument("--out", default="target/wheels")
    parser.add_argument("--no-deps", action="store_true")
    options = parser.parse_args(arguments)
    host = options.distribution == "community"
    profile = (
        "release"
        if options.release
        else options.profile
        or os.environ.get("GAMMALOOP_BUILD_PROFILE")
        or (ROOT / ".config/build-profile").read_text().strip()
    )
    if host and profile == "dev-optim":
        raise ValueError(
            "The Community host has no dev-optim profile; use dev, ci-optim or release"
        )
    interpreter, details = python_environment(options.python, options.install)
    manifest = (
        "examples/notebooks/symbolica-host/Cargo.toml"
        if host
        else "crates/gammaloop-api/Cargo.toml"
    )
    pyproject = "examples/notebooks/symbolica-host" if host else "."
    features = (
        []
        if host
        else ["pyo3-extension-module"] + (["python_abi"] if options.abi3 else [])
    )
    output_dir = (ROOT / options.out).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    distribution = "symbolica" if host else "gammaloop"
    if backend == "nix":
        if (
            details["implementation"] != "cpython"
            or details["free_threaded"]
            or details.get("debug", False)
        ):
            raise ValueError(
                "Nix wheels currently require a standard CPython interpreter"
            )
        normalize_architecture = {"arm64": "aarch64", "AMD64": "x86_64"}
        interpreter_arch = details.get("machine", platform.machine())
        if normalize_architecture.get(
            interpreter_arch, interpreter_arch
        ) != normalize_architecture.get(platform.machine(), platform.machine()):
            raise ValueError(
                "The selected Python architecture does not match the Nix host"
            )
        version = details["version"]
        spec = {
            "package": HOST if host else "gammaloop-api",
            "manifest": manifest,
            "pyproject": pyproject,
            "profile": profile,
            "features": features,
            "default_features": host,
            "python": f"python{version[0]}{version[1]}",
            "args": [],
        }
        # Keep the output and extracted RPATH references alive for installed wheels.
        root_directory = (
            Path(details["prefix"]) / ".nix-python" if options.install else output_dir
        )
        root = root_directory / f".pending-{distribution}-{uuid.uuid4().hex}"
        try:
            output = nix_build("development-wheel", spec=spec, root=root)
        except BaseException:
            root.unlink(missing_ok=True)
            raise
        wheels = list((output / "dist").glob("*.whl"))
        if len(wheels) != 1:
            raise ValueError(f"Expected exactly one wheel in {output / 'dist'}")
        wheel = output_dir / wheels[0].name
        # Store files are read-only. Export through a writable temporary file so
        # repeat builds can atomically replace the same wheel, then retain its
        # own closure even after its development environment is removed.
        with (output_dir / ".gammaloop-wheel.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            temporary = output_dir / f".{uuid.uuid4().hex}.whl"
            try:
                shutil.copyfile(wheels[0], temporary)
                temporary.replace(wheel)
            finally:
                temporary.unlink(missing_ok=True)
            retain_output(output, output_dir / f".nix-{wheel.name}", capture=True)
        if not options.install:
            root.unlink()
        installation_wheel = wheels[0]
    else:
        command = [
            "maturin",
            "build",
            "--locked",
            "--manifest-path",
            manifest,
            "--profile",
            profile,
            "--interpreter",
            interpreter,
            "--out",
            str(output_dir),
        ]
        if not host:
            command += ["--no-default-features", "--features", ",".join(features)]
        before = {path: path.stat().st_mtime_ns for path in output_dir.glob("*.whl")}
        run(command)
        wheels = [
            path
            for path in output_dir.glob("*.whl")
            if before.get(path) != path.stat().st_mtime_ns
        ]
        if len(wheels) != 1:
            raise ValueError("Maturin did not produce exactly one new wheel")
        wheel = wheels[0]
        installation_wheel = wheel
    print(f"Wheel: {wheel}")
    if options.install:
        # uv respects the chosen environment and reinstalls only this distribution.
        command = [
            "uv",
            "pip",
            "install",
            "--python",
            interpreter,
            "--reinstall-package",
            distribution,
        ]
        if options.no_deps:
            command.append("--no-deps")
        # Serialize installation and root replacement together. Keep the old
        # runtime rooted until the new wheel is successfully installed; a failed
        # install retains its pending root for inspection or retry.
        prefix = Path(details["prefix"])
        with (prefix / ".gammaloop-install.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            if backend == "nix":
                native = sorted(
                    [
                        *(output / "lib").rglob("*.so"),
                        *(output / "lib").rglob("*.dylib"),
                    ]
                )
                if not native:
                    raise ValueError("The Nix wheel has no retained native payload")
                # Probe the target interpreter's actual loader before replacing
                # packages, including glibc/libstdc++ compatibility on non-NixOS.
                run(
                    [
                        interpreter,
                        "-c",
                        "import ctypes,sys; [ctypes.CDLL(path) for path in sys.argv[1:]]",
                        *map(str, native),
                    ],
                    env=license_environment(),
                )
            run(command + [str(installation_wheel)])
            imports = (
                "import symbolica.community.hepkit, symbolica.community.graph, symbolica.community.tensor"
                if host
                else "from gammaloop import GammaLoopAPI"
            )
            run([interpreter, "-c", imports], env=license_environment())
            if backend == "nix":
                permanent_root = prefix / ".nix-python" / distribution
                retain_output(output, permanent_root)
                root.unlink()


def main(arguments=None):
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument(
        "--backend", choices=("nix", "cargo"), action="append", required=True
    )
    options, arguments = parser.parse_known_args(arguments)
    if len(options.backend) != 1:
        parser.error("--backend may be given only once")
    python_package(arguments, options.backend[0])


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as error:
        sys.exit(error.returncode)
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(2)
