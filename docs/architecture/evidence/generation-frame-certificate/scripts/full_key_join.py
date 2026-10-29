import json
from pathlib import Path

root = Path("target/frame-certificate")
head_map = json.loads(
    Path(
        "../generation-integration-plan/target/wp3-drawing-transplant/frame-digest-boundary.json"
    ).read_text()
)["spellings"]


def old_chart(value):
    for current, old in head_map.items():
        value = value.replace(current, old)
    return value


def local_fields(payload):
    fields = [
        ("numerator", payload["numerator"]),
        ("numerator_prefactor", payload["numerator_prefactor"]),
        ("projector", payload["projector"]),
    ]
    fields.extend(
        (f"vertex:{i}", vertex["numerator"])
        for i, vertex in enumerate(payload["vertices"])
    )
    fields.extend(
        (f"edge:{i}", edge[1]["numerator"]) for i, edge in enumerate(payload["edges"])
    )
    return dict(fields)


results = []
for case in ("comparison-qed", "representative-qed", "representative-lr"):
    paired = {}
    for side in ("historical", "current"):
        final = json.loads((root / f"{side}-frozen/{case}.json").read_text())
        pair_root = root / f"pair-captures-{side}-graphproof/{case}/frame-pairs"
        sources = {
            int(path.name.split("-")[1]): path
            for path in pair_root.iterdir()
            if path.is_dir()
        }
        key = lambda diagram: json.dumps(
            diagram["key"], sort_keys=True, separators=(",", ":")
        )
        by_key = {key(diagram): diagram for diagram in final["diagrams"]}
        assert len(by_key) == len(final["diagrams"])
        masters = set()
        members = set()
        final_field_checks = 0
        for group in final["groups"]:
            index = group["master"]
            assert index not in masters and index < len(final["diagrams"])
            masters.add(index)
            first = group["members"][0]
            representative = json.loads(
                (sources[first["source_index"]] / "representative.json").read_text()
            )
            output = final["diagrams"][index]
            assert representative["name"] == output["name"] == first["name"]
            assert representative["id"] == int(output["id"], 16)
            assert representative.keys() == output["payload"].keys() | {"id", "name"}
            for field in output["payload"]:
                if field != "overall_factor":
                    assert representative[field] == output["payload"][field], (
                        side,
                        case,
                        index,
                        field,
                    )
                    final_field_checks += 1
            terms = []
            for member in group["members"]:
                source_index = member["source_index"]
                assert source_index not in members and source_index in sources
                members.add(source_index)
                folder = sources[source_index]
                source = json.loads((folder / "source.json").read_text())
                header = dict(
                    line.split("\t", 1)
                    for line in (folder / "frame.txt").read_text().splitlines()
                    if line.startswith(("source_index\t", "source_name\t"))
                )
                assert header == {
                    "source_index": str(source_index),
                    "source_name": member["name"],
                }
                assert member["name"] == source["name"]
                assert int(member["id"], 16) == source["id"]
                assert member["overall_factor"] == source["overall_factor"]
                assert member["ratio"] == "1"
                terms.append(
                    f"feynkit_generator_factor::{{}}::NumeratorDependentGrouping({source_index},1,{source['overall_factor']})"
                )
            if len(group["members"]) == 1:
                assert (
                    output["payload"]["overall_factor"]
                    == representative["overall_factor"]
                )
            else:
                assert output["payload"]["overall_factor"] == "+".join(terms)
        assert masters == set(range(len(final["diagrams"])))
        assert members == set(sources)
        paired[side] = (by_key, final, final_field_checks)

    old_by_key, old, old_checks = paired["historical"]
    current_by_key, current, current_checks = paired["current"]
    assert old_by_key.keys() == current_by_key.keys()
    assert old["groups"] == current["groups"]
    expression_checks = 0
    physical_checks = 0
    for graph_key in old_by_key:
        left = old_by_key[graph_key]
        right = current_by_key[graph_key]
        for field in ("id", "name", "aut", "basis", "overall_factor"):
            assert left[field] == right[field], (case, graph_key, field)
            physical_checks += 1
        a = left["payload"]
        b = right["payload"]
        assert a.keys() == b.keys()
        for field in a:
            if field not in (
                "edges",
                "vertices",
                "numerator",
                "numerator_prefactor",
                "projector",
            ):
                assert a[field] == b[field], (case, graph_key, field)
                physical_checks += 1
        assert len(a["vertices"]) == len(b["vertices"])
        assert len(a["edges"]) == len(b["edges"])
        for old_vertex, current_vertex in zip(a["vertices"], b["vertices"]):
            assert {
                key: value
                for key, value in old_vertex.items()
                if key not in ("numerator", "num_interface")
            } == {
                key: value
                for key, value in current_vertex.items()
                if key not in ("numerator", "num_interface")
            }
        for old_edge, current_edge in zip(a["edges"], b["edges"]):
            assert old_edge[0] == current_edge[0]
            assert {
                key: value
                for key, value in old_edge[1].items()
                if key not in ("numerator", "num_interface")
            } == {
                key: value
                for key, value in current_edge[1].items()
                if key not in ("numerator", "num_interface")
            }
        old_factors = local_fields(a)
        current_factors = local_fields(b)
        assert old_factors.keys() == current_factors.keys()
        for owner in old_factors:
            assert old_factors[owner] == old_chart(current_factors[owner]), (
                case,
                graph_key,
                owner,
            )
            expression_checks += 1
    results.append(
        {
            "case": case,
            "full_graph_key_matches": len(old_by_key),
            "unmatched_keys": 0,
            "exact_groups": len(old["groups"]),
            "source_members": sum(len(group["members"]) for group in old["groups"]),
            "historical_final_representative_field_bindings": old_checks,
            "current_final_representative_field_bindings": current_checks,
            "old_current_final_factorized_expression_bindings": expression_checks,
            "old_current_final_physical_field_bindings": physical_checks,
            "complete_final_master_and_group_weight_binding": True,
        }
    )
(root / "full-key-join-summary.json").write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps(results, indent=2))
