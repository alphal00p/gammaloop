"""Link driver.rs against explicit compatible cached Rust libraries; never rebuild them."""

import argparse
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
for name in ["idenso", "spenso", "symbolica", "serde-json"]:
    p.add_argument("--" + name, type=Path, required=True)
p.add_argument("--dependency-dir", type=Path, action="append", required=True)
p.add_argument("--output", type=Path, required=True)
p.add_argument("--rustc", default="rustc")
p.add_argument("--linker", default="cc")
p.add_argument("--cpus")
a = p.parse_args()
source = Path(__file__).with_name("driver.rs")
a.output.parent.mkdir(parents=True, exist_ok=True)
sha = lambda f: hashlib.sha256(Path(f).read_bytes()).hexdigest()
inputs = {
    str(f): sha(f) for f in [source, a.idenso, a.spenso, a.symbolica, a.serde_json]
}
command = (["taskset", "-c", a.cpus] if a.cpus else []) + [
    a.rustc,
    "--edition=2024",
    "--crate-name",
    "trusted_algebra_benchmark",
    "-C",
    "opt-level=2",
    "-C",
    "debuginfo=1",
    "-C",
    "debug-assertions=yes",
    "-C",
    "overflow-checks=yes",
    "-C",
    "linker=" + a.linker,
]
for d in a.dependency_dir:
    command += ["-L", "dependency=" + str(d)]
for name in ["idenso", "spenso", "symbolica", "serde_json"]:
    command += ["--extern", name + "=" + str(getattr(a, name))]
command += [str(source), "-o", str(a.output)]
start = time.time_ns()
run = subprocess.run(
    command,
    check=False,
    env=os.environ
    | {
        "CARGO_CRATE_NAME": "trusted_algebra_benchmark",
        "CARGO_PKG_NAME": "trusted_algebra_benchmark",
        "CARGO_PKG_VERSION": "0.0.0",
    },
)
assert all(sha(f) == h for f, h in inputs.items())
receipt = {
    "command": command,
    "exit_code": run.returncode,
    "elapsed_seconds": (time.time_ns() - start) / 1e9,
    "inputs": inputs,
}
if run.returncode == 0:
    receipt["binary_sha256"] = sha(a.output)
a.output.with_suffix(".build.json").write_text(json.dumps(receipt, indent=2) + "\n")
raise SystemExit(run.returncode)
