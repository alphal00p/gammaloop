"""Direct, non-nested test/capture processes; explicit selection and pool delivery."""

import argparse
import re
import shutil
import sys

import common as M

OLD = (
    M.WORKSPACE
    / "target/oracle-numerator-family-diagnosis/pristine-current/bin/generation_oracle"
)


def sandbox(elf, args, directory, extras=None):
    directory.mkdir(parents=True, exist_ok=False)
    for name in ("tmp", "cache", "dumps", "state"):
        (directory / name).mkdir()
    bwrap = shutil.which("bwrap")
    argv = [
        bwrap,
        "--ro-bind",
        "/",
        "/",
        "--dev",
        "/dev",
        "--proc",
        "/proc",
        "--bind",
        directory,
        directory,
        "--bind",
        directory / "tmp",
        "/tmp",
        "--bind",
        directory / "cache",
        "/home/lcnbr/.cache",
        "--unshare-net",
        "--die-with-parent",
        "--",
        sys.executable,
        "-B",
        M.HERE / "runtime_exec.py",
        elf,
        *args,
    ]
    env = {
        "TMPDIR": str(directory / "tmp"),
        "XDG_CACHE_HOME": str(directory / "cache"),
        "GL_TEST_LOG_DIR": str(directory / "state"),
        **(extras or {}),
    }
    return M.run(
        argv,
        directory / "execution",
        env_extra=env,
        inputs=(
            elf,
            M.HERE / "runtime_exec.py",
            M.HERE / "runtime.py",
            M.HERE / "common.py",
        ),
    )


def tests(route):
    if route == "old":
        elf = OLD
        M.require(
            M.sha(elf)
            == "133899560973c3b824a1fde41061c0cf32995438434436567f1bfaa035d9d7ba",
            "provenance",
            "old producer",
        )
    else:
        elf = M.OUT / "bin/generation_oracle"
    rows = []
    for case in M.CASES:
        for pool in M.POOLS:
            d = M.OUT / "runtime-final" / route / case / f"pool-{pool}"
            proc, record = sandbox(
                elf,
                [case, "--exact", "--test-threads", "1", "--nocapture"],
                d,
                {
                    "FEYNKIT_ORACLE_THREADS": str(pool),
                    "FEYNKIT_ORACLE_DUMP": str(d / "dumps"),
                },
            )
            matches = re.findall(
                r"test result: (ok|FAILED)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out",
                proc.stdout,
            )
            M.require(len(matches) == 1, "selection", case + ": actual libtest summary")
            status, passed, failed, ignored, measured, filtered = matches[0]
            actual = [int(v) for v in (passed, failed, ignored, measured)]
            M.require(
                actual == ([0, 1, 0, 0] if route == "old" else [1, 0, 0, 0]),
                "test-result",
                str(actual),
            )
            M.require(
                proc.returncode == (101 if route == "old" else 0),
                "test-exit",
                str(proc.returncode),
            )
            if route == "old":
                M.require(
                    f"{case}: summary differs from its oracle row" in proc.stderr,
                    "old-failure",
                    case,
                )
            for mode in ("None", "Auto"):
                for suffix in ("json", "sources.json"):
                    M.require(
                        (
                            d / "dumps" / f"{case}.classes-{mode}-t{pool}.{suffix}"
                        ).is_file(),
                        "pool-delivery",
                        case,
                    )
            rows.append(
                {
                    "case": case,
                    "pool": pool,
                    "route": route,
                    "selected": int(passed) + int(failed),
                    "passed": int(passed),
                    "failed": int(failed),
                    "ignored": int(ignored),
                    "filtered": int(filtered),
                    "status": status,
                    "command": record,
                    "dump_files": M.C.files(d / "dumps"),
                }
            )
            M.write(M.OUT / "runtime-final" / route / "results.json", rows)
            print(route, case, pool, status, flush=True)


def captures(route):
    elf = M.OUT / "bin/oracle_family_migration_capture"
    rows = []
    for case in M.CASES:
        for pool in M.POOLS:
            d = M.OUT / "companions" / route / case / f"pool-{pool}"
            proc, record = sandbox(elf, [d / "dumps", str(pool), case], d)
            M.require(proc.returncode == 0, "companion", case)
            rows.append(
                {
                    "case": case,
                    "pool": pool,
                    "route": route,
                    "command": record,
                    "files": M.C.files(d / "dumps"),
                }
            )
            M.write(M.OUT / "companions" / route / "results.json", rows)
            print("companion", route, case, pool, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=("tests", "captures"))
    parser.add_argument("route", choices=("old", "new"))
    args = parser.parse_args()
    (tests if args.kind == "tests" else captures)(args.route)
