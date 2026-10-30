"""Resolve actual Cargo dependency-unit feature/profile closures; no root shorthand."""

import argparse
from pathlib import Path

import common as M


def collect(root, unit):
    units = {}
    for path in root.glob("*/*.json"):
        stored = path.with_suffix("")
        units.setdefault(path.parent.name, []).append(
            {
                "json": str(path),
                "data": M.load(path),
                "stored_fingerprint": stored.read_text().strip()
                if stored.is_file()
                else None,
            }
        )
    index = {
        r["stored_fingerprint"]: (name, r)
        for name, rows in units.items()
        for r in rows
        if r["stored_fingerprint"]
    }
    data = M.load(unit)
    seen, todo, closure, files = (
        set(),
        [d[3] for d in data["deps"]],
        [],
        {str(unit): M.sha(unit)},
    )
    while todo:
        key = todo.pop().to_bytes(8, "little").hex()
        if key in seen:
            continue
        seen.add(key)
        M.require(key in index, "fingerprints", "missing dependency unit " + key)
        name, rec = index[key]
        d = rec["data"]
        files[rec["json"]] = M.sha(rec["json"])
        closure.append(
            {
                "unit": name,
                "kind": Path(rec["json"]).stem,
                "features": d["features"],
                "profile": d["profile"],
                "rustc": d["rustc"],
                "rustflags": d["rustflags"],
                "data": d,
            }
        )
        todo += [dep[3] for dep in d["deps"]]
    return {
        "root": data,
        "closure": sorted(closure, key=lambda r: (r["unit"], r["kind"])),
        "fingerprint_files": files,
        "closure_units": len(closure),
        "limits": "actual dependency units; path/config hashes may reflect build-root changes; compare feature/profile/compiler projections explicitly",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("artifact", type=Path)
    ap.add_argument("fingerprint", type=Path)
    args = ap.parse_args()
    before = M.sha(args.artifact)
    data = collect(args.fingerprint.parent.parent, args.fingerprint)
    data.update(artifact=str(args.artifact), artifact_sha256=before)
    M.require(
        M.sha(args.artifact) == before,
        "fingerprint",
        "artifact changed during collection",
    )
    M.write(M.OUT / "provenance" / (args.label + ".json"), data)
    print(args.label, before, data["closure_units"])


if __name__ == "__main__":
    main()
