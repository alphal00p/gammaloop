"""Derive class-order and CT occurrence records, and bind immutable sources."""

import collections
import hashlib
import importlib.util
import os
from pathlib import Path
import subprocess

SPEC = importlib.util.spec_from_file_location(
    "certificate", Path(__file__).with_name("validate.py")
)
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)
OUT = (
    Path(
        os.environ.get(
            "FAMILY_ROUND2_OUT",
            str(V.WORKSPACE / "target/oracle-numerator-family-round2"),
        )
    )
    / "qualification"
)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    rows, ct = [], {}
    for route, root in V.ROOTS.items():
        for case in V.CASES:
            np = V.load(root / case / "classes-None-t1/payload.json")
            ap = V.load(root / case / "classes-Auto-t1/payload.json")
            ka = {str(g["key"]): (i, g) for i, g in enumerate(ap["diagrams"])}
            V.require(
                len(ka) == len(ap["diagrams"]) == len(np["diagrams"]),
                "class-domain",
                case,
            )
            moved = []
            for i, g in enumerate(np["diagrams"]):
                j, other = ka[str(g["key"])]
                V.require(
                    {k: v for k, v in g.items() if k != "name"}
                    == {k: v for k, v in other.items() if k != "name"},
                    "class-payload",
                    g["id"],
                )
                ng, ag = np["groups"][i], ap["groups"][j]
                for group, graph, index in ((ng, g, i), (ag, other, j)):
                    V.require(
                        group
                        == {
                            "master": index,
                            "members": [
                                {
                                    "id": graph["id"],
                                    "name": graph["name"],
                                    "ratio": "1",
                                    "overall_factor": graph["overall_factor"],
                                    "source_index": index,
                                }
                            ],
                        },
                        "class-group-chain",
                        graph["id"],
                    )
                if i != j:
                    moved.append(
                        {
                            "none": i,
                            "auto": j,
                            "id": g["id"],
                            "none_name": g["name"],
                            "auto_name": other["name"],
                            "none_group": ng,
                            "auto_group": ag,
                        }
                    )
            nt = V.load(root / case / "classes-None-t1/typed-owners.json")
            at = V.load(root / case / "classes-Auto-t1/typed-owners.json")
            V.require(
                {g["id"]: g["owners"] for g in nt}
                == {g["id"]: g["owners"] for g in at},
                "class-typed",
                case,
            )
            rows.append(
                {"route": route, "case": case, "graphs": len(ka), "permutation": moved}
            )
            if case == "ctorder_qcd_uugg_L1":
                for mode, payload in (("None", np), ("Auto", ap)):
                    folder = root / case / f"classes-{mode}-t1"
                    ct[route + "/" + mode] = V.ct_provenance(
                        V.load(folder / "model.json"),
                        payload,
                        V.load(folder / "request.json"),
                    )
    V.write(OUT / "class-order.json", rows)
    V.write(OUT / "ct-occurrences.json", ct)
    sources = []
    for revision in (
        "2b1ef4f07a44bc061bcc97c378bbdd9c41cb0e94",
        "b9f788922341e7bde2e801a054423ae70ed7e9ae",
    ):
        for path in (
            "crates/feynkit-graph/src/symbols.rs",
            "crates/idenso/src/color/mod.rs",
            "crates/idenso/src/cook.rs",
            "Cargo.lock",
        ):
            argv = ["jj", "file", "show", "--ignore-working-copy", "-r", revision, path]
            proc = subprocess.run(argv, cwd=V.WORKSPACE, capture_output=True)
            if proc.returncode:
                raise RuntimeError(proc.stderr.decode())
            dest = OUT / "sources" / revision / path
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(proc.stdout)
            sources.append(
                {"revision": revision, "path": path, "sha256": sha(dest), "argv": argv}
            )
    V.write(OUT / "source-bindings.json", sources)
    V.write(
        OUT / "elf-bindings.json",
        [
            {
                "kernel": k,
                "path": str(p),
                "sha256": sha(p),
                "build": "unchanged round-one retained ELF; original build-log gap retained, corroborated by independent round-one review reruns",
            }
            for k, p in V.ELFS.items()
        ],
    )
    print([(k, x["rule_counts"]) for k, x in ct.items()])


if __name__ == "__main__":
    main()
