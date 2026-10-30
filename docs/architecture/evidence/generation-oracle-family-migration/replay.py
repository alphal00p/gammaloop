"""Unmodified approved portable entrypoints; all output belongs to migration."""

import concurrent.futures
import shutil
import sys

import common as M


def replay(kind):
    d = M.OUT / "replay" / kind
    d.mkdir(parents=True, exist_ok=False)
    (d / "cache").mkdir()
    suite = M.HERE.parent / "generation-oracle-family-diagnosis/round3/suite.py"
    proc, record = M.run(
        [
            shutil.which("bwrap"),
            "--ro-bind",
            "/",
            "/",
            "--dev",
            "/dev",
            "--proc",
            "/proc",
            "--bind",
            d,
            d,
            "--bind",
            d / "cache",
            "/home/lcnbr/.cache",
            "--unshare-net",
            "--die-with-parent",
            "--",
            sys.executable,
            "-B",
            suite,
            kind,
            "--out",
            d / "results",
        ],
        d / "execution",
        env_extra={"XDG_CACHE_HOME": str(d / "cache")},
        inputs=tuple(M.C.source_hashes()),
    )
    M.require(proc.returncode == 0, "replay", kind)
    rows = M.load(d / "results/results.json")
    M.require(rows and all(r["matched"] for r in rows), "replay", kind)
    return {
        "kind": kind,
        "count": len(rows),
        "command": record,
        "results_sha256": M.sha(d / "results/results.json"),
    }


if __name__ == "__main__":
    kinds = (
        "baseline",
        "legacy",
        "reviewer",
        "domain",
        "staged",
        "anchor",
        "infrastructure",
    )
    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        for row in executor.map(replay, kinds):
            rows.append(row)
            M.write(M.OUT / "replay/results.json", rows)
            print(row["kind"], row["count"], flush=True)
