import collections
import json
from pathlib import Path

root = Path("target/frame-certificate")
heads = json.loads(
    Path(
        "../generation-integration-plan/target/wp3-drawing-transplant/frame-digest-boundary.json"
    ).read_text()
)["spellings"]


def old_chart(value):
    for current, old in heads.items():
        value = value.replace(current, old)
    return value


results = []
for case in ("comparison-qed", "representative-qed", "representative-lr"):
    old_root = root / "pair-captures-historical-graphproof" / case / "frame-pairs"
    new_root = root / "pair-captures-current-graphproof" / case / "frame-pairs"
    names = sorted(x.name for x in old_root.iterdir() if x.is_dir())
    counts = collections.Counter()
    first = {}
    for name in names:
        for stage in ("source", "comparison", "representative"):
            old = json.loads((old_root / name / f"{stage}.json").read_text())
            new = json.loads((new_root / name / f"{stage}.json").read_text())
            assert old.keys() == new.keys()
            assert len(old["vertices"]) == len(new["vertices"])
            assert len(old["edges"]) == len(new["edges"])
            assert old["loop_momentum_basis"] == new["loop_momentum_basis"]
            counts[stage + "_diagrams"] += 1
            for field in old.keys() - {
                "vertices",
                "edges",
                "numerator",
                "numerator_prefactor",
                "projector",
            }:
                if old[field] != new[field]:
                    counts["other_field_mismatch"] += 1
                    first.setdefault("other", (case, name, stage, field))
            pairs = [
                (field, old[field], new[field])
                for field in ("numerator", "numerator_prefactor", "projector")
            ]
            for i, (a, b) in enumerate(zip(old["vertices"], new["vertices"])):
                if {
                    k: v
                    for k, v in a.items()
                    if k not in ("numerator", "num_interface")
                } != {
                    k: v
                    for k, v in b.items()
                    if k not in ("numerator", "num_interface")
                }:
                    counts["vertex_physical_mismatch"] += 1
                    first.setdefault("vertex", (case, name, stage, i))
                pairs.append((f"vertex:{i}", a["numerator"], b["numerator"]))
            for i, (a, b) in enumerate(zip(old["edges"], new["edges"])):
                if a[0] != b[0] or {
                    k: v
                    for k, v in a[1].items()
                    if k not in ("numerator", "num_interface")
                } != {
                    k: v
                    for k, v in b[1].items()
                    if k not in ("numerator", "num_interface")
                }:
                    counts["edge_physical_mismatch"] += 1
                    first.setdefault("edge", (case, name, stage, i))
                pairs.append((f"edge:{i}", a[1]["numerator"], b[1]["numerator"]))
            for owner, a, b in pairs:
                counts["expression_factors"] += 1
                if a != old_chart(b):
                    counts["factor_mismatch"] += 1
                    first.setdefault("factor", (case, name, stage, owner, a, b))
    assert not first and not any("mismatch" in key for key in counts)
    results.append({"case": case, "counts": dict(counts), "first": first})
(root / "pair-expression-summary.json").write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps(results, indent=2))
