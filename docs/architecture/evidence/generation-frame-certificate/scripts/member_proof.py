import collections
import json
from pathlib import Path

root = Path("target/frame-certificate")
results = []
for side in ("historical", "current"):
    for case in ("comparison-qed", "representative-qed", "representative-lr"):
        final = json.loads((root / f"{side}-frozen" / f"{case}.json").read_text())
        pair_root = root / f"pair-captures-{side}-graphproof" / case / "frame-pairs"
        sources = {
            int(folder.name.split("-")[1]): folder
            for folder in pair_root.iterdir()
            if folder.is_dir()
        }
        members = [member for group in final["groups"] for member in group["members"]]
        assert set(sources) == {member["source_index"] for member in members}
        assert len(sources) == len(members)
        counts = collections.Counter()
        first = {}
        masters = set()
        for group in final["groups"]:
            assert group["master"] not in masters and group["master"] < len(
                final["diagrams"]
            )
            masters.add(group["master"])
            master = group["members"][0]
            master_folder = sources[master["source_index"]]
            master_cmp = json.loads((master_folder / "comparison.json").read_text())
            representative = json.loads(
                (master_folder / "representative.json").read_text()
            )
            output = final["diagrams"][group["master"]]
            assert representative["name"] == master["name"] == output["name"]
            assert representative["id"] == int(output["id"], 16)
            assert representative.keys() == output["payload"].keys() | {"id", "name"}
            for field in output["payload"]:
                if field != "overall_factor":
                    assert representative[field] == output["payload"][field], (
                        side,
                        case,
                        group["master"],
                        field,
                    )
                    counts["final_representative_fields"] += 1
            master_ports = [
                line
                for line in (master_folder / "frame.txt").read_text().splitlines()
                if line.startswith("ports\tcomparison")
            ]
            grouped_terms = []
            for member in group["members"]:
                folder = sources[member["source_index"]]
                source = json.loads((folder / "source.json").read_text())
                comparison = json.loads((folder / "comparison.json").read_text())
                ports = [
                    line
                    for line in (folder / "frame.txt").read_text().splitlines()
                    if line.startswith("ports\tcomparison")
                ]
                counts["members"] += 1
                if (
                    member["name"] != source["name"]
                    or int(member["id"], 16) != source["id"]
                ):
                    counts["source_identity_mismatch"] += 1
                    first.setdefault(
                        "identity", (group["master"], member["source_index"])
                    )
                if member["overall_factor"] != source["overall_factor"]:
                    counts["source_factor_mismatch"] += 1
                    first.setdefault(
                        "factor", (group["master"], member["source_index"])
                    )
                if member["ratio"] != "1":
                    counts["nonunit_ratio"] += 1
                    first.setdefault(
                        "ratio",
                        (group["master"], member["source_index"], member["ratio"]),
                    )
                if any(
                    comparison[field] != master_cmp[field]
                    for field in (
                        "numerator",
                        "numerator_prefactor",
                        "projector",
                        "loop_momentum_basis",
                    )
                ):
                    counts["comparison_value_mismatch"] += 1
                    first.setdefault(
                        "comparison", (group["master"], member["source_index"])
                    )
                if ports != master_ports:
                    counts["comparison_port_mismatch"] += 1
                    first.setdefault("ports", (group["master"], member["source_index"]))
                if member is not master:
                    counts["additional_members"] += 1
                grouped_terms.append(
                    f"feynkit_generator_factor::{{}}::NumeratorDependentGrouping({member['source_index']},1,{source['overall_factor']})"
                )
            if len(group["members"]) == 1:
                assert (
                    output["payload"]["overall_factor"]
                    == representative["overall_factor"]
                )
            else:
                assert output["payload"]["overall_factor"] == "+".join(grouped_terms)
        assert masters == set(range(len(final["diagrams"])))
        assert not first and not any(
            "mismatch" in key or "nonunit" in key for key in counts
        )
        results.append(
            {
                "side": side,
                "case": case,
                "groups": len(final["groups"]),
                "sources": len(sources),
                "counts": dict(counts),
                "first": first,
            }
        )
(root / "member-summary.json").write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps(results, indent=2))
