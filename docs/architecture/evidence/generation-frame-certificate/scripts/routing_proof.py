import collections
import json
from fractions import Fraction
from pathlib import Path

root = Path("target/frame-certificate")


def read_frame(path):
    maps = collections.defaultdict(dict)
    graphs = collections.defaultdict(dict)
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if parts[0] == "edge_map":
            maps[parts[1]][int(parts[2])] = (int(parts[3]), int(parts[4]))
        elif parts[0] == "graph_edge":
            graphs[parts[1]][int(parts[2])] = (
                parts[3],
                parts[4],
                parts[5],
                int(parts[6]),
            )
    return maps, graphs


def finalized_edges(provisional, diagram, cross_section):
    internal = [
        (index, edge)
        for index, edge in provisional.items()
        if not any(endpoint.startswith("x:") for endpoint in edge[:2])
    ]
    final_internal = [
        (index, edge)
        for index, edge in enumerate(diagram["edges"])
        if edge[1]["external"] is None
    ]
    assert len(internal) == len(final_internal)
    mapping = dict(
        zip((index for index, _ in internal), (index for index, _ in final_internal))
    )
    orientation = {}
    external_owners = collections.defaultdict(list)
    for (old, source), (new, target) in zip(internal, final_internal):
        before = tuple(int(endpoint[2:]) for endpoint in source[:2])
        after = (target[0]["source"], target[0]["target"])
        assert before == after or before == after[::-1]
        assert int(source[2].split("(")[1][:-1]) == target[1]["particle"]
        assert source[3] == int(target[1]["directed"])
    for old, edge in provisional.items():
        external = next(
            (endpoint for endpoint in edge[:2] if endpoint.startswith("x:")), None
        )
        if external is not None:
            vertex = int((edge[0] if edge[0].startswith("v:") else edge[1])[2:])
            particle = int(external.split("ParticleId(")[1][:-1])
            candidates = [
                index
                for index, target in enumerate(diagram["edges"])
                if target[1]["external"] is not None
                and vertex in target[0].values()
                and target[1]["particle"] == particle
            ]
            assert len(candidates) == 1, (old, external, candidates)
            mapping[old] = candidates[0]
            external_owners[candidates[0]].append(external)
        final = diagram["edges"][mapping[old]][0]
        if edge[0].startswith("v:") and edge[1].startswith("v:"):
            before = (int(edge[0][2:]), int(edge[1][2:]))
            after = (final["source"], final["target"])
            assert before == after or before == after[::-1]
            orientation[old] = int(before != after)
        else:
            attached = int((edge[0] if edge[0].startswith("v:") else edge[1])[2:])
            assert attached in final.values()
            orientation[old] = int(
                (final["source"] != attached)
                if edge[0].startswith("v:")
                else (final["target"] != attached)
            )
    assert set(mapping) == set(provisional)
    assert set(mapping.values()) == set(range(len(diagram["edges"])))
    for owners in external_owners.values():
        if cross_section:
            assert len(owners) == 2
            assert {owner.split(":")[2] for owner in owners} == {"Incoming", "Outgoing"}
        else:
            assert len(owners) == 1
    return mapping, orientation


def signed_final_map(
    source_final, target_final, provisional_map, source_frame, target_frame
):
    mapped = {}
    for old, (new, reversed_edge) in provisional_map.items():
        source_edge = source_final[0][old]
        target_edge = target_final[0][new]
        sign = -1 if source_final[1][old] ^ reversed_edge ^ target_final[1][new] else 1
        if source_edge in mapped:
            assert mapped[source_edge] == (target_edge, sign)
        mapped[source_edge] = (target_edge, sign)
    assert set(mapped) == set(range(len(source_frame["edges"])))
    assert {edge for edge, _ in mapped.values()} == set(
        range(len(target_frame["edges"]))
    )
    return mapped


def coordinates(basis):
    return [
        (kind, edge)
        for kind, field in (("loops", "loop_edges"), ("external", "external_edges"))
        for edge in basis[field]
    ]


def replay(source, target, signed_map):
    source_basis = source["loop_momentum_basis"]
    target_basis = target["loop_momentum_basis"]
    old_coordinates = coordinates(source_basis)
    new_coordinates = coordinates(target_basis)
    assert len(old_coordinates) == len(new_coordinates)
    transformation = {}
    for old_index, (kind, old_edge) in enumerate(old_coordinates):
        new_edge, sign = signed_map[old_edge]
        candidates = [
            index
            for index, (other_kind, edge) in enumerate(new_coordinates)
            if kind == other_kind and edge == new_edge
        ]
        assert len(candidates) == 1, (kind, old_edge, new_edge, candidates)
        transformation[candidates[0]] = (old_index, sign)
    assert len(transformation) == len(old_coordinates)
    for old_edge, (new_edge, edge_sign) in signed_map.items():
        old_signature = source_basis["edge_signatures"][str(old_edge)]
        new_signature = target_basis["edge_signatures"][str(new_edge)]
        source_vector = [
            Fraction(value)
            for value in old_signature["loops"] + old_signature["external"]
        ]
        target_vector = [
            Fraction(value)
            for value in new_signature["loops"] + new_signature["external"]
        ]
        assert len(source_vector) == len(old_coordinates)
        assert len(target_vector) == len(new_coordinates)
        mapped = [Fraction(0) for _ in old_coordinates]
        for index, coefficient in enumerate(target_vector):
            old_index, sign = transformation[index]
            mapped[old_index] += coefficient * sign
        assert mapped == [edge_sign * coefficient for coefficient in source_vector], (
            old_edge,
            new_edge,
            mapped,
            source_vector,
            transformation,
        )
    return len(signed_map), len(old_coordinates)


results = []
for case in ("comparison-qed", "representative-qed", "representative-lr"):
    roots = {
        side: root / f"pair-captures-{side}-graphproof" / case / "frame-pairs"
        for side in ("historical", "current")
    }
    names = sorted(path.name for path in roots["historical"].iterdir() if path.is_dir())
    assert names == sorted(
        path.name for path in roots["current"].iterdir() if path.is_dir()
    )
    counts = collections.Counter()
    for name in names:
        for side in ("historical", "current"):
            folder = roots[side] / name
            maps, graphs = read_frame(folder / "frame.txt")
            source = json.loads((folder / "source.json").read_text())
            receipt = json.loads(
                (
                    root / f"model-hooks-{side}" / case / "effective-bindings.json"
                ).read_text()
            )
            by_name = {record["name"]: record for record in receipt["second_build"]}
            assert source["name"] in by_name
            assert (
                source["loop_momentum_basis"]["loop_edges"]
                == by_name[source["name"]]["loop_edges"]
            )
            source_final = finalized_edges(
                graphs["source"], source, case == "representative-lr"
            )
            for stage in ("comparison", "representative"):
                target = json.loads((folder / f"{stage}.json").read_text())
                target_final = finalized_edges(
                    graphs[stage], target, case == "representative-lr"
                )
                signed_map = signed_final_map(
                    source_final, target_final, maps[stage], source, target
                )
                edges, coordinates_checked = replay(source, target, signed_map)
                counts[side + "_" + stage + "_transported_edges"] += edges
                counts[side + "_" + stage + "_basis_coordinates"] += coordinates_checked
                counts[side + "_" + stage + "_complete_signed_affine_checks"] += 1
        for stage in ("source", "comparison", "representative"):
            old = json.loads((roots["historical"] / name / f"{stage}.json").read_text())
            current = json.loads(
                (roots["current"] / name / f"{stage}.json").read_text()
            )
            assert old["loop_momentum_basis"] == current["loop_momentum_basis"]
            counts[stage + "_old_current_exact_basis_pairs"] += 1
    results.append({"case": case, "counts": dict(counts)})
(root / "routing-summary.json").write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps(results, indent=2))
