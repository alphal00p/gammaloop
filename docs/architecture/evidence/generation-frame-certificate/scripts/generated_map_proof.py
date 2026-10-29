import collections
import json
from pathlib import Path

root = Path("target/frame-certificate")


def read(path):
    graph = collections.defaultdict(lambda: {"vertices": {}, "edges": {}})
    maps = collections.defaultdict(dict)
    for line in path.read_text().splitlines():
        p = line.split("\t")
        if p[0] == "graph_vertex":
            graph[p[1]]["vertices"][int(p[2])] = p[3]
        elif p[0] == "graph_edge":
            graph[p[1]]["edges"][int(p[2])] = (p[3], p[4], p[5], p[6])
        elif p[0] == "vertex_map":
            maps[(p[1], "vertex")][int(p[2])] = int(p[3])
        elif p[0] == "edge_map":
            maps[(p[1], "edge")][int(p[2])] = (int(p[3]), int(p[4]))
    return graph, maps


def check(graph, maps, case):
    counts = collections.Counter()
    first = {}
    source = graph["source"]
    for stage in ("comparison", "representative"):
        target = graph[stage]
        vm = maps[(stage, "vertex")]
        em = maps[(stage, "edge")]
        assert set(vm) == set(source["vertices"]) and set(vm.values()) == set(
            target["vertices"]
        )
        assert set(em) == set(source["edges"]) and {j for j, _ in em.values()} == set(
            target["edges"]
        )
        counts[stage + "_vertex_bijections"] += 1
        counts[stage + "_edge_bijections"] += 1
        external_map = {}
        for i, (j, flipped) in em.items():
            source_ends = source["edges"][i][:2]
            target_ends = target["edges"][j][:2]
            target_ends = target_ends[::-1] if flipped else target_ends
            for before, after in zip(source_ends, target_ends):
                if before.startswith("x:"):
                    if before in external_map:
                        assert external_map[before] == after
                    external_map[before] = after
        assert len(set(external_map.values())) == len(external_map)
        labels = sorted(external_map)
        actual = tuple(external_map[label] for label in labels)
        identity = tuple(labels)
        # The left/right pass exchanges the two incoming photon labels with
        # their paired outgoing labels and their states. No other relabelling is admitted.
        if case == "representative-lr":
            assert actual in (
                identity,
                tuple(
                    f"x:{2 if n == 0 else 3 if n == 1 else 0 if n == 2 else 1}:"
                    f"{'Outgoing' if n < 2 else 'Incoming'}:ParticleId(0)"
                    for n in (int(label.split(":")[1]) for label in labels)
                ),
            )
        else:
            assert actual == identity
        if actual != identity:
            counts[stage + "_lr_swapped_sources"] += 1
        for i, j in vm.items():
            counts[stage + "_vertices"] += 1
            if source["vertices"][i] != target["vertices"][j]:
                counts[stage + "_rule_changes"] += 1
        for i, (j, flipped) in em.items():
            a = source["edges"][i]
            b = target["edges"][j]
            endpoint = lambda x, vm=vm, external_map=external_map: (
                f"v:{vm[int(x[2:])]}" if x.startswith("v:") else external_map[x]
            )
            expected = (endpoint(a[0]), endpoint(a[1]))
            if flipped:
                expected = expected[::-1]
            counts[stage + "_edges"] += 1
            counts[stage + "_flips"] += flipped
            if expected != b[:2]:
                counts[stage + "_endpoint_failures"] += 1
                first.setdefault(stage, (i, j, expected, b[:2], flipped))
            if a[2] != b[2]:
                counts[stage + "_particle_changes"] += 1
            if a[3] != b[3]:
                counts[stage + "_directed_changes"] += 1
    return counts, first


results = []
for case in ("comparison-qed", "representative-qed", "representative-lr"):
    current_root = root / "pair-captures-current-graphproof" / case / "frame-pairs"
    historical_root = (
        root / "pair-captures-historical-graphproof" / case / "frame-pairs"
    )
    current_names = sorted(x.name for x in current_root.iterdir() if x.is_dir())
    historical_names = sorted(x.name for x in historical_root.iterdir() if x.is_dir())
    assert current_names == historical_names
    counts = collections.Counter()
    first = {}
    for name in current_names:
        current_graph, current_maps = read(current_root / name / "frame.txt")
        old_graph, old_maps = read(historical_root / name / "frame.txt")
        assert current_graph == old_graph
        assert current_maps == old_maps
        checked, failures = check(current_graph, current_maps, case)
        counts.update(checked)
        if failures:
            first.setdefault("map", (name, failures))
    assert not first and not any("failures" in key for key in counts)
    results.append(
        {
            "case": case,
            "sources": len(current_names),
            "counts": dict(counts),
            "first": first,
        }
    )
(root / "generated-map-summary.json").write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps(results, indent=2))
