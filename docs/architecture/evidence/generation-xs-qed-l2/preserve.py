"""Strict required prior evidence closure; no cache/secret-derived hashes."""

import argparse
import datetime

import common as M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=("before", "after"))
    phase = ap.parse_args().phase
    old = M.WORKSPACE / "target/oracle-family-migration"
    paths = set()
    for file in (
        old / "preservation/after.json",
        old / "inventory/artifact-manifest.json",
    ):
        manifest = M.load(file)
        paths.add(str(file))
        paths.update(manifest["files"])
    for file in (
        old / "freeze.json",
        M.WORKSPACE
        / ".local/review-family-migration-opus-20261030T050603Z/review.json",
    ):
        M.require(file.is_file(), "preservation", str(file))
        paths.add(str(file))
    root = M.WORKSPACE / ".local/review-family-migration-opus-20261030T050603Z"
    paths.update(
        str(p)
        for p in root.rglob("*")
        if p.is_file()
        and not any(k in p.relative_to(root).parts for k in ("cache", "tmp", "state"))
    )
    paths.update(
        str(p)
        for p in (M.HERE.parent / "generation-oracle-family-migration").iterdir()
        if p.is_file()
    )
    files = {
        p: {"sha256": M.sha(p), "bytes": M.Path(p).stat().st_size}
        for p in sorted(paths)
    }
    if phase == "after":
        M.require(
            M.load(M.OUT / "preservation/before.json")["files"] == files,
            "preservation",
            "prior evidence changed",
        )
    M.write(
        M.OUT / "preservation" / (phase + ".json"),
        {
            "utc": datetime.datetime.now(datetime.UTC).isoformat(),
            "files": files,
            "scope": "required immutable approved/migration closure plus new migration review; no private cache/tmp/state values hashed",
        },
    )
    print(phase, len(files), "prior evidence files")


if __name__ == "__main__":
    main()
