import collections
import json
from pathlib import Path

root = Path("target/frame-certificate")
case_assignments = {
    "comparison-qed": 66,
    "representative-qed": 66,
    "representative-lr": 3456,
}


def transport(source, target, frame_lines, stage):
    final_edges = {
        int(p[1]): int(p[2])
        for line in frame_lines
        if (p := line.split("\t"))[0] == "final_edge_map"
    }
    vertices = {
        int(p[2]): int(p[3])
        for line in frame_lines
        if (p := line.split("\t"))[0] == "vertex_map" and p[1] == stage
    }
    generated_edges = collections.defaultdict(dict)
    for line in frame_lines:
        p = line.split("\t")
        if p[0] == "graph_edge":
            generated_edges[p[1]][int(p[2])] = p[3:5]
    provisional = {
        int(p[2]): (int(p[3]), int(p[4]))
        for line in frame_lines
        if (p := line.split("\t"))[0] == "edge_map" and p[1] == stage
    }
    external = {}
    for old, (new, flip) in provisional.items():
        left = generated_edges["source"][old]
        right = generated_edges[stage][new]
        right = right[::-1] if flip else right
        for before, after in zip(left, right):
            if before.startswith("x:"):
                external[before] = after
    mirror = any(before != after for before, after in external.items())
    flips = {}
    for old, new in final_edges.items():
        before = source["edges"][old][0]
        after = target["edges"][new][0]
        mapped = (vertices[before["source"]], vertices[before["target"]])
        actual = (after["source"], after["target"])
        assert mapped == actual or mapped == actual[::-1]
        flips[old] = mapped != actual

    def cuts(diagram, map_source):
        result = []
        for cut in diagram["cuts"]:
            half_edges = []
            for half_edge in cut["cut"]:
                edge = half_edge["edge"]
                endpoint = half_edge["endpoint"]
                if map_source:
                    endpoint_flip = flips[edge] ^ mirror
                    edge = final_edges[edge]
                    if endpoint_flip:
                        endpoint = "target" if endpoint == "source" else "source"
                half_edges.append((edge, endpoint))
            result.append(tuple(sorted(half_edges)))
        return sorted(result)

    assert cuts(source, True) == cuts(target, False)
    return len(source["cuts"]), mirror


summary = []
full = []
for case in ("comparison-qed", "representative-qed", "representative-lr"):
    pair_roots = {
        side: root / f"pair-captures-{side}-graphproof" / case / "frame-pairs"
        for side in ("historical", "current")
    }
    source_names = {
        side: sorted(x.name for x in folder.iterdir() if x.is_dir())
        for side, folder in pair_roots.items()
    }
    assert source_names["historical"] == source_names["current"]
    counts = collections.Counter()
    for name in source_names["current"]:
        record = {"case": case, "source": name, "stages": {}}
        for stage in ("source", "comparison", "representative"):
            old = json.loads(
                (pair_roots["historical"] / name / f"{stage}.json").read_text()
            )
            current = json.loads(
                (pair_roots["current"] / name / f"{stage}.json").read_text()
            )
            assert old["cuts"] == current["cuts"]
            assert old["overall_factor"] == current["overall_factor"]
            assert old["symmetry_factor"] == current["symmetry_factor"]
            assert old["name"] == current["name"]
            counts[stage + "_cuts"] += len(current["cuts"])
            counts[stage + "_sources"] += 1
            record["stages"][stage] = {
                "cuts": current["cuts"],
                "overall_factor": current["overall_factor"],
                "symmetry_factor": current["symmetry_factor"],
                "name": current["name"],
            }
        for side in ("historical", "current"):
            frame_lines = (
                (pair_roots[side] / name / "frame.txt").read_text().splitlines()
            )
            source_diagram = json.loads(
                (pair_roots[side] / name / "source.json").read_text()
            )
            for stage in ("comparison", "representative"):
                target_diagram = json.loads(
                    (pair_roots[side] / name / f"{stage}.json").read_text()
                )
                cut_count, mirrored = transport(
                    source_diagram, target_diagram, frame_lines, stage
                )
                counts[side + "_" + stage + "_transported_cuts"] += cut_count
                counts[side + "_" + stage + "_mirrored_sources"] += mirrored
                if side == "current":
                    record["stages"][stage]["source_cut_transport_mirrored"] = mirrored
        full.append(record)
    grouped = json.loads((root / f"current-frozen/{case}.json").read_text())
    historical_grouped = json.loads(
        (root / f"historical-frozen/{case}.json").read_text()
    )
    assert grouped["groups"] == historical_grouped["groups"]
    members = [member for group in grouped["groups"] for member in group["members"]]
    assert len(members) == len(source_names["current"])
    assert {member["source_index"] for member in members} == {
        int(name.split("-")[1]) for name in source_names["current"]
    }
    summary.append(
        {
            "case": case,
            "original_interaction_assignments": case_assignments[case],
            "post_physical_selection_built_sources": len(source_names["current"]),
            "group_representatives": len(grouped["groups"]),
            "group_members": len(members),
            "removed_after_built_source_boundary": 0,
            "counts": dict(counts),
        }
    )
(root / "cut-sign-provenance-full.json").write_text(
    json.dumps(full, ensure_ascii=False, indent=2) + "\n"
)
(root / "cut-sign-provenance-summary.json").write_text(
    json.dumps(summary, indent=2) + "\n"
)
print(json.dumps(summary, indent=2))
