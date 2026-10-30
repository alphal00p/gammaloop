"""Immutable proof/review preservation guard (never prunes reviewer evidence)."""

import argparse
import datetime

import common as M


def roots():
    return [
        M.HERE.parent / "generation-oracle-family-diagnosis",
        M.WORKSPACE / "target/oracle-numerator-family-diagnosis",
        M.WORKSPACE.parent
        / "generation-oracle-family-historical/target/oracle-numerator-family-diagnosis",
        *[
            M.WORKSPACE / "target" / p
            for p in (
                "oracle-numerator-family-round2",
                "oracle-numerator-family-round2-final",
                "oracle-numerator-family-round2-candidate",
                "oracle-numerator-family-round3",
                "oracle-family-independent-review-r1",
                "oracle-numerator-families-preflight",
            )
        ],
        *[
            p
            for p in (M.WORKSPACE / ".local").glob("review-family-r[23]-*")
            if p.is_dir()
        ],
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("before", "after"))
    phase = parser.parse_args().phase
    files = {}
    for root in roots():
        M.require(root.is_dir(), "preserve", str(root))
        files.update(M.C.files(root))
    value = {
        "utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "files": files,
        "chronology": "global preservation snapshot; per-run guards and immutable approved manifests give strict consumed-input before/after bindings; this snapshot began after first runtime/replay starts",
    }
    M.write(M.OUT / "preservation" / (phase + ".json"), value)
    if phase == "after":
        M.require(
            M.load(M.OUT / "preservation/before.json")["files"] == files,
            "preserve",
            "original proof/review changed",
        )
    print(
        phase,
        len(files),
        "paths",
        sum(v["bytes"] for v in files.values()),
        "logical bytes including aliases",
    )


if __name__ == "__main__":
    main()
