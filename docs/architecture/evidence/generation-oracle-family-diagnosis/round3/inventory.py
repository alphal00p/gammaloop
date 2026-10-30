"""Read-only closure inventories; no secret values or derivation contents."""

import argparse
import datetime

import common as C


def roots():
    w = C.WORKSPACE
    return [
        w / "target/oracle-numerator-family-diagnosis",
        w.parent
        / "generation-oracle-family-historical/target/oracle-numerator-family-diagnosis",
        w / "target/oracle-numerator-family-round2-candidate",
        w / "target/oracle-family-independent-review-r1",
        w / ".local/review-family-r2-a-20261030T010834Z",
        w / ".local/review-family-r2-b-20261030T010814",
        w / "target/oracle-numerator-families-preflight",
        w / "docs/architecture/evidence/generation-oracle-family-diagnosis/round2",
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phase", choices=("before", "after", "final"))
    args = ap.parse_args()
    if args.phase == "final":
        files = C.files(C.OUT)
        # Self and final command are excluded: they cannot hash their own result.
        files = {
            p: h
            for p, h in files.items()
            if "/inventory/closure-final.json" not in p
            and not p.startswith(str(C.OUT / "commands") + "/")
            and "/attempts/" not in p
        }
        files.update(C.files(C.HERE))
        for p in (
            C.WORKSPACE
            / "docs/architecture/generation-oracle-family-diagnosis-round3.typ",
            C.WORKSPACE / "docs/developers.toml",
        ):
            files[str(p)] = {"sha256": C.sha(p), "bytes": p.stat().st_size}
        C.write(
            C.OUT / "inventory/closure-final.json",
            {
                "files": files,
                "exclusions": "self, ongoing top-level command records and development attempts; commands are separately inventoried after completion",
            },
        )
        return
    files = {}
    for root in roots():
        C.require(
            root.is_dir(), "inventory", "missing frozen closure root " + str(root)
        )
        files.update(C.files(root))
    old = C.files(C.HERE.parent)
    files.update({p: h for p, h in old.items() if "/round3/" not in p})
    for name in (
        "generation-oracle-family-diagnosis.typ",
        "generation-oracle-family-diagnosis-round2.typ",
    ):
        p = C.WORKSPACE / "docs/architecture" / name
        files[str(p)] = {"sha256": C.sha(p), "bytes": p.stat().st_size}
    value = {
        "utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "files": files,
        "source_hashes": C.source_hashes(),
    }
    C.write(C.OUT / ("inventory/frozen-" + args.phase + ".json"), value)
    if args.phase == "after":
        before = C.load(C.OUT / "inventory/frozen-before.json")
        C.require(
            before["files"] == files,
            "frozen-changed",
            "frozen captures/reviews/prior evidence changed",
        )
        C.require(
            before["source_hashes"] == value["source_hashes"],
            "source-changed",
            "final executable sources changed",
        )
    print(args.phase, len(files), sum(f["bytes"] for f in files.values()))


if __name__ == "__main__":
    main()
