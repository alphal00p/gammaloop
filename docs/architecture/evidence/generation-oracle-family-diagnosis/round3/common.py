"""Recorded round-three diagnostic operations; no secrets in records."""

import datetime
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[4]
OUT = WORKSPACE / "target/oracle-numerator-family-round3"
SPEC = importlib.util.spec_from_file_location(
    "frozen_round2", HERE.parent / "round2/validate.py"
)
B = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(B)


class Rejection(Exception):
    def __init__(self, stage, detail):
        super().__init__(detail)
        self.stage = stage


def require(value, stage, detail):
    if not value:
        raise Rejection(stage, detail)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n")


def load(path):
    def unique(items):
        result = {}
        for k, v in items:
            require(k not in result, "json", "duplicate key " + k)
            result[k] = v
        return result

    try:
        return json.loads(
            Path(path).read_text(),
            object_pairs_hook=unique,
            parse_constant=lambda x: (_ for _ in ()).throw(
                Rejection("json", "nonfinite " + x)
            ),
        )
    except json.JSONDecodeError as error:
        raise Rejection("json", str(error)) from error


def strict_equal(a, b):
    if type(a) is not type(b):
        return False
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(strict_equal(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return len(a) == len(b) and all(strict_equal(x, y) for x, y in zip(a, b))
    return a == b


def source_hashes():
    # Executable proof/harness dependency closure. Auxiliary report, registry
    # and inventory scripts are separately hashed in the final source manifest;
    # they neither enter native pairs nor change semantic validation.
    files = [
        HERE / name for name in ("common.py", "validate.py", "suite.py", "schema.json")
    ]
    files += [HERE.parent / "round2" / name for name in ("validate.py", "controls.py")]
    return {str(p): sha(p) for p in files if p.exists()}


def files(root):
    return {
        str(p.absolute()): {"sha256": sha(p), "bytes": p.stat().st_size}
        for p in sorted(Path(root).rglob("*"))
        if p.is_file() and p.suffix != ".drv"
    }


def run(argv, destination, runtime=False, extra_env=None, bindings=None):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    for key in list(env):
        if key.startswith("FEYNKIT_ORACLE_") or key in (
            "SYMBOLICA_LICENSE",
            "SYMBOLICA_LICENSE_SIGNED",
            "CARGO_TARGET_DIR",
        ):
            env.pop(key)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    if extra_env:
        env.update(extra_env)
    if runtime:
        # A real license is injected only into this direct child, never Nix.
        license_value = os.environ.get("SYMBOLICA_LICENSE")
        if not license_value:
            config = Path(
                os.environ.get("XDG_CONFIG_HOME", str(Path.home() / ".config"))
            )
            license_value = (config / "symbolica/license").read_text().strip()
        if not license_value:
            raise RuntimeError("missing runtime license")
        env["SYMBOLICA_LICENSE"] = license_value
    argv = [str(x) for x in argv]
    started = datetime.datetime.now(datetime.UTC).isoformat()
    tick = time.monotonic()
    record = {
        "argv": argv,
        "cwd": str(WORKSPACE),
        "started_utc": started,
        "environment_policy": {
            "license": "private direct-child injection"
            if runtime
            else "both license variables absent",
            "signed_license_absent": True,
            "oracle_flags_absent": True,
            "CARGO_TARGET_DIR_absent": True,
            "PYTHONDONTWRITEBYTECODE": "1",
            "XDG_CONFIG_HOME": env.get("XDG_CONFIG_HOME", "default user config"),
        },
        "source_hashes": source_hashes(),
    }
    if bindings:
        record["bindings_before"] = bindings
    write(destination / "command.json", record)
    proc = subprocess.run(
        argv, cwd=WORKSPACE, env=env, capture_output=True, text=True, check=False
    )
    (destination / "stdout.txt").write_text(proc.stdout)
    (destination / "stderr.txt").write_text(proc.stderr)
    record.update(
        exit=proc.returncode,
        elapsed_seconds=time.monotonic() - tick,
        ended_utc=datetime.datetime.now(datetime.UTC).isoformat(),
    )
    write(destination / "command.json", record)
    return proc, record


if __name__ == "__main__":
    # Record actual top-level invocations, including their outcome and times.
    label, *argv = sys.argv[1:]
    stamp = datetime.datetime.now(datetime.UTC).strftime("%Y%m%dT%H%M%S%f")
    proc, record = run(argv, OUT / "commands" / (stamp + "-" + label))
    print(proc.stdout, end="")
    print(proc.stderr, end="", file=sys.stderr)
    sys.exit(proc.returncode)
