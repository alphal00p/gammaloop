"""Migration-only evidence operations; frozen certificate sources are read-only."""

import datetime
import importlib.util
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[3]
OUT = WORKSPACE / "target/oracle-family-migration"
APPROVED = "87fa13a93b1c0b3e5792fd56937486ad962b2b55"
CASES = ("qed_eemm_L1", "ctorder_qcd_uugg_L1", "sm_eemm_L0", "ym_g6_L0")
POOLS = (1, 8, 32)
MAIN = "crates/feynkit-generator/tests/generation_oracle/main.rs"
SPEC = importlib.util.spec_from_file_location(
    "approved_common",
    HERE.parent / "generation-oracle-family-diagnosis/round3/common.py",
)
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)
sha, load, write, require = C.sha, C.load, C.write, C.require


def source_hashes():
    return {str(p): sha(p) for p in sorted(HERE.glob("*.py"))}


def clean_env():
    return {
        k: v
        for k, v in os.environ.items()
        if k
        not in (
            "SYMBOLICA_LICENSE",
            "SYMBOLICA_LICENSE_SIGNED",
            "INSTA_UPDATE",
            "UPDATE_EXPECT",
            "EXPECTTEST_ACCEPT",
            "BLESS",
            "UPDATE_SNAPSHOTS",
            "RECORD",
            "FORCE_PASS",
        )
        and not k.startswith("FEYNKIT_ORACLE_")
    }


def run(argv, directory, *, cwd=WORKSPACE, env_extra=None, inputs=()):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    argv = [str(a) for a in argv]
    env = clean_env()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env.update(env_extra or {})
    require(
        "SYMBOLICA_LICENSE" not in env and "SYMBOLICA_LICENSE_SIGNED" not in env,
        "env",
        "license must enter only exec runtime child",
    )
    tool = Path(shutil.which(argv[0]) or argv[0])
    before = {str(p): sha(p) for p in inputs}
    for source in sorted(HERE.glob("*.py")):
        saved = OUT / "helper-versions" / sha(source) / source.name
        if not saved.exists():
            saved.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, saved)
    record = {
        "argv": argv,
        "cwd": str(cwd),
        "started_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "sources_before": source_hashes(),
        "inputs_before": before,
        "tool_sha256_before": sha(tool),
        "environment": {
            "both_license_variables": "absent; runtime_exec injects privately after sandbox entry",
            "oracle_flags": {
                k: v for k, v in env.items() if k.startswith("FEYNKIT_ORACLE_")
            },
            "CARGO_TARGET_DIR": env.get(
                "CARGO_TARGET_DIR", "unset; normal shell/cache only"
            ),
            "TMPDIR": env.get("TMPDIR"),
            "XDG_CACHE_HOME": env.get("XDG_CACHE_HOME"),
            "PYTHONDONTWRITEBYTECODE": "1",
            "update_force_flags_absent": {
                k: k not in env
                for k in (
                    "INSTA_UPDATE",
                    "UPDATE_EXPECT",
                    "EXPECTTEST_ACCEPT",
                    "BLESS",
                    "UPDATE_SNAPSHOTS",
                    "RECORD",
                    "FORCE_PASS",
                )
            },
        },
    }
    write(directory / "command.json", record)
    tick = time.monotonic()
    proc = subprocess.run(
        argv, cwd=cwd, env=env, capture_output=True, text=True, check=False
    )
    (directory / "stdout.txt").write_text(proc.stdout)
    (directory / "stderr.txt").write_text(proc.stderr)
    record.update(
        exit=proc.returncode,
        ended_utc=datetime.datetime.now(datetime.UTC).isoformat(),
        elapsed_seconds=time.monotonic() - tick,
        inputs_after={p: sha(p) for p in before},
        tool_sha256_after=sha(tool),
        sources_after=source_hashes(),
        stdout_sha256=sha(directory / "stdout.txt"),
        stderr_sha256=sha(directory / "stderr.txt"),
    )
    require(
        record["inputs_before"] == record["inputs_after"]
        and record["tool_sha256_before"] == record["tool_sha256_after"],
        "run",
        "consumed input/tool changed",
    )
    write(directory / "command.json", record)
    return proc, record


def immutable(revision, path):
    proc, record = run(
        ["jj", "file", "show", "--ignore-working-copy", "-r", revision, path],
        OUT / "sources" / revision / path,
    )
    require(proc.returncode == 0, "source", path)
    return proc.stdout, record


if __name__ == "__main__":
    label, *argv = sys.argv[1:]
    stamp = datetime.datetime.now(datetime.UTC).strftime("%Y%m%dT%H%M%S%f")
    proc, _ = run(argv, OUT / "commands" / (stamp + "-" + label))
    print(proc.stdout, end="")
    print(proc.stderr, end="", file=sys.stderr)
    sys.exit(proc.returncode)
