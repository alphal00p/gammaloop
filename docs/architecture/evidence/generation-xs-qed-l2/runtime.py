"""Actual pristine assertion and separate live-owner captures, pool1 only."""

import argparse
import re
import shutil
import sys

import common as M


def run(route):
    root = M.OUT / "runtime" / route
    root.mkdir(parents=True, exist_ok=False)
    for name in ("tmp", "cache", "state", "dumps", "captures"):
        (root / name).mkdir()
    qualification = M.OUT / "build" / route
    records = []
    runner = M.HERE.parent / "generation-oracle-family-migration/runtime_exec.py"
    for name, args in (
        (
            "generation_oracle",
            [M.CASES[0], "--exact", "--test-threads", "1", "--nocapture"],
        ),
        ("oracle_xs_qed_l2_capture", [root / "captures", M.CASES[0]]),
    ):
        elf = M.OUT / "bin" / route / name
        q = M.load(qualification / name / "qualification.json")
        M.require(M.sha(elf) == q["artifact_sha256"], "elf", name)
        argv = [
            shutil.which("bwrap"),
            "--ro-bind",
            "/",
            "/",
            "--dev",
            "/dev",
            "--proc",
            "/proc",
            "--bind",
            root,
            root,
            "--bind",
            root / "tmp",
            "/tmp",
            "--bind",
            root / "cache",
            "/home/lcnbr/.cache",
            "--unshare-net",
            "--die-with-parent",
            "--",
            sys.executable,
            "-B",
            runner,
            elf,
            *args,
        ]
        env = {
            "TMPDIR": str(root / "tmp"),
            "XDG_CACHE_HOME": str(root / "cache"),
            "GL_TEST_LOG_DIR": str(root / "state"),
        }
        if name == "generation_oracle":
            env.update(
                FEYNKIT_ORACLE_THREADS="1", FEYNKIT_ORACLE_DUMP=str(root / "dumps")
            )
        proc, record = M.run(
            argv,
            root / name,
            env_extra=env,
            inputs=[elf, runner, *[M.HERE / n for n in ("runtime.py", "common.py")]],
        )
        if name == "generation_oracle":
            found = re.findall(
                r"test result: (ok|FAILED)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out",
                proc.stdout,
            )
            M.require(len(found) == 1, "selection", "one actual libtest summary")
            selected = [int(v) for v in found[0][1:5]]
            M.require(
                selected == ([1, 0, 0, 0] if route == "historical" else [0, 1, 0, 0]),
                "assertion",
                str(found),
            )
            M.require(
                proc.returncode == (0 if route == "historical" else 101),
                "assertion-exit",
                route,
            )
            if route == "current":
                M.require(
                    f"{M.CASES[0]}: summary differs from its oracle row" in proc.stderr,
                    "assertion-boundary",
                    route,
                )
            record["actual_test_summary"] = found
        else:
            M.require(proc.returncode == 0, "capture", route)
        records.append(record)
        M.write(root / "results.json", records)
        print(route, name, proc.returncode, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("route", choices=M.REVISIONS)
    run(p.parse_args().route)
