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


def read_frame(path):
    header = {}
    maps = collections.defaultdict(list)
    ports = collections.defaultdict(dict)
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if parts[0] in (
            "source_index",
            "source_name",
            "requested_prefactor",
            "comparison_prefactor",
            "representative_prefactor",
        ):
            header[parts[0]] = parts[1]
        elif parts[0] in ("vertex_map", "edge_map"):
            maps[(parts[0], parts[1])].append(tuple(map(int, parts[2:])))
        elif parts[0] == "final_edge_map":
            maps[("final_edge_map", "representative")].append(
                tuple(map(int, parts[1:]))
            )
        elif parts[0] in ("graph_vertex", "graph_edge"):
            continue
        elif parts[0] == "ports":
            ports[parts[1]][parts[2]] = (
                json.loads(parts[3]),
                int(parts[4]),
                int(parts[5]),
            )
        else:
            raise AssertionError(parts[0])
    return header, maps, ports


def declared(diagram, owner):
    if owner.startswith("vertex:"):
        return diagram["vertices"][int(owner.split(":")[1])].get("num_interface")
    if owner.startswith("edge:"):
        return diagram["edges"][int(owner.split(":")[1])][1].get("num_interface")
    return diagram.get(
        {
            "prefactor": "numerator_prefactor_interface",
            "projector": "projector_interface",
        }[owner]
    )


results = []
with (root / "all-pair-axis-maps.jsonl").open("w") as out:
    for case in ("comparison-qed", "representative-qed", "representative-lr"):
        old_root = root / "pair-captures-historical-graphproof" / case / "frame-pairs"
        new_root = root / "pair-captures-current-graphproof" / case / "frame-pairs"
        names = sorted(x.name for x in old_root.iterdir() if x.is_dir())
        assert names == sorted(x.name for x in new_root.iterdir() if x.is_dir())
        counts = collections.Counter()
        first = {}
        for name in names:
            old_header, old_maps, old_ports = read_frame(old_root / name / "frame.txt")
            new_header, new_maps, new_ports = read_frame(new_root / name / "frame.txt")
            assert old_maps == new_maps
            assert old_header["source_index"] == new_header["source_index"]
            assert old_header["source_name"] == new_header["source_name"]
            counts["source_pairs"] += 1
            for field in (
                "requested_prefactor",
                "comparison_prefactor",
                "representative_prefactor",
            ):
                assert old_header[field] == old_chart(new_header[field]), (
                    case,
                    name,
                    field,
                )
                counts[field + "_exact"] += 1
            for stage in ("source", "comparison", "representative"):
                old_diagram = json.loads(
                    (old_root / name / f"{stage}.json").read_text()
                )
                new_diagram = json.loads(
                    (new_root / name / f"{stage}.json").read_text()
                )
                assert old_ports[stage].keys() == new_ports[stage].keys()
                assert (
                    old_diagram["loop_momentum_basis"]
                    == new_diagram["loop_momentum_basis"]
                )
                counts[stage + "_routing_equal"] += 1
                for owner, (old_logical, old_metric, old_composite) in old_ports[
                    stage
                ].items():
                    current_logical, current_metric, current_composite = new_ports[
                        stage
                    ][owner]
                    current_in_old_chart = [old_chart(x) for x in current_logical]
                    assert (
                        old_metric == current_metric
                        and old_composite == current_composite
                    )
                    assert len(set(old_logical)) == len(old_logical)
                    assert len(set(current_in_old_chart)) == len(current_in_old_chart)
                    if collections.Counter(old_logical) != collections.Counter(
                        current_in_old_chart
                    ):
                        first.setdefault(
                            "port_multiset",
                            (
                                case,
                                name,
                                stage,
                                owner,
                                old_logical,
                                current_in_old_chart,
                            ),
                        )
                        counts["port_multiset_mismatch"] += 1
                        continue
                    permutation = [current_in_old_chart.index(x) for x in old_logical]
                    if permutation != list(range(len(permutation))):
                        counts["permuted_owners"] += 1
                    old_decl = declared(old_diagram, owner)
                    new_decl = declared(new_diagram, owner)
                    if old_decl is not None:
                        assert old_decl == old_logical
                        counts["old_explicit_interfaces"] += 1
                    if new_decl is not None:
                        assert new_decl == current_logical
                        counts["current_explicit_interfaces"] += 1
                    counts["owners"] += 1
                    counts["axes"] += len(permutation)
                    out.write(
                        json.dumps(
                            {
                                "case": case,
                                "source": name,
                                "stage": stage,
                                "owner": owner,
                                "old_logical": old_logical,
                                "current_logical_in_old_chart": current_in_old_chart,
                                "old_axis_to_current_axis": permutation,
                                "metric": old_metric,
                                "composite": old_composite,
                                "old_declared": old_decl is not None,
                                "current_declared": new_decl is not None,
                            }
                        )
                        + "\n"
                    )
        assert not first and counts["port_multiset_mismatch"] == 0
        results.append({"case": case, "counts": dict(counts), "first": first})
(root / "axis-and-frame-summary.json").write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps(results, indent=2))
