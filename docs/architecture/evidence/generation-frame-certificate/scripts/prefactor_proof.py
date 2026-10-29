import collections
import json
import re
from fractions import Fraction
from pathlib import Path

root = Path("target/frame-certificate")
heads = json.loads(
    Path(
        "../generation-integration-plan/target/wp3-drawing-transplant/frame-digest-boundary.json"
    ).read_text()
)["spellings"]
q_head = (
    "gammalooprs::{spenso::rank1,spenso::tensor,spenso::component-labels:subscript}::Q"
)
representation = "spenso::{spenso::representation,spenso::self_dual}::mink"
vertex_head = "gammalooprs::{spenso::index}::vertex"
q_call = (
    re.escape(q_head)
    + r"\((\d+),"
    + re.escape(representation)
    + r"\(4,"
    + re.escape(vertex_head)
    + r"\((\d+),9\)\)\)"
)
term = re.compile(r"(?:(-?\d+(?:/\d+)?)\*)?" + q_call + r"\*" + q_call)


def current_chart(source):
    for current, historical in heads.items():
        source = source.replace(historical, current)
    return source


def terms(source):
    values = []
    for piece in current_chart(source).split("+"):
        match = term.fullmatch(piece)
        assert match is not None, piece
        coefficient, first_edge, first_vertex, second_edge, second_vertex = (
            match.groups()
        )
        assert first_vertex == second_vertex
        values.append(
            (
                Fraction(coefficient or 1),
                tuple(sorted((int(first_edge), int(second_edge)))),
                int(first_vertex),
            )
        )
    assert len(values) == 2
    return sorted(values)


def read_frame(path):
    header = {}
    maps = collections.defaultdict(dict)
    graphs = collections.defaultdict(dict)
    ports = collections.defaultdict(dict)
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if parts[0].endswith("prefactor"):
            header[parts[0]] = parts[1]
        elif parts[0] == "vertex_map":
            maps[(parts[1], "vertex")][int(parts[2])] = int(parts[3])
        elif parts[0] == "edge_map":
            maps[(parts[1], "edge")][int(parts[2])] = (int(parts[3]), int(parts[4]))
        elif parts[0] == "graph_edge":
            graphs[parts[1]][int(parts[2])] = (
                parts[3],
                parts[4],
                parts[5],
                int(parts[6]),
            )
        elif parts[0] == "ports":
            ports[parts[1]][parts[2]] = json.loads(parts[3])
    return header, maps, graphs, ports


def signed_transport(requested, vertices, edges):
    return sorted(
        (
            coefficient * (-1 if edges[first][1] ^ edges[second][1] else 1),
            tuple(sorted((edges[first][0], edges[second][0]))),
            vertices[vertex],
        )
        for coefficient, (first, second), vertex in requested
    )


def finalization_map(provisional, diagram, cross_section):
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
    reversed_edges = {}
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
        final_endpoints = diagram["edges"][mapping[old]][0]
        if edge[0].startswith("v:") and edge[1].startswith("v:"):
            before = (int(edge[0][2:]), int(edge[1][2:]))
            after = (final_endpoints["source"], final_endpoints["target"])
            assert before == after or before == after[::-1]
            reversed_edges[old] = int(before != after)
        else:
            attached = int((edge[0] if edge[0].startswith("v:") else edge[1])[2:])
            assert attached in final_endpoints.values()
            reversed_edges[old] = int(
                (final_endpoints["source"] != attached)
                if edge[0].startswith("v:")
                else (final_endpoints["target"] != attached)
            )
    assert set(mapping) == set(provisional)
    assert set(mapping.values()) == set(range(len(diagram["edges"])))
    for external_index, owners in external_owners.items():
        if cross_section:
            assert len(owners) == 2
            assert {owner.split(":")[2] for owner in owners} == {"Incoming", "Outgoing"}
        else:
            assert len(owners) == 1
        assert diagram["edges"][external_index][1]["external"] is not None
    return mapping, reversed_edges


def finalized_terms(provisional_terms, mapping, reversed_edges):
    return sorted(
        (
            coefficient * (-1 if reversed_edges[first] ^ reversed_edges[second] else 1),
            tuple(sorted((mapping[first], mapping[second]))),
            vertex,
        )
        for coefficient, (first, second), vertex in provisional_terms
    )


results = []
for side in ("historical", "current"):
    for case in ("comparison-qed", "representative-qed", "representative-lr"):
        counts = collections.Counter()
        receipt = json.loads(
            (
                root / f"model-hooks-{side}" / case / "effective-bindings.json"
            ).read_text()
        )
        pair_root = root / f"pair-captures-{side}-graphproof" / case / "frame-pairs"
        for folder in pair_root.iterdir():
            if not folder.is_dir():
                continue
            header, maps, graphs, ports = read_frame(folder / "frame.txt")
            assert header["requested_prefactor"] == receipt["requested_prefactor"]
            requested = terms(header["requested_prefactor"])
            assert requested == [(1, (0, 4), 0), (1, (1, 6), 1)]
            for stage in ("source", "comparison", "representative"):
                provisional = (
                    requested
                    if stage == "source"
                    else signed_transport(
                        requested, maps[(stage, "vertex")], maps[(stage, "edge")]
                    )
                )
                if stage != "source":
                    assert provisional == terms(header[f"{stage}_prefactor"])
                    counts["provisional_signed_transport"] += 1
                diagram = json.loads((folder / f"{stage}.json").read_text())
                assert ports[stage]["prefactor"] == []
                assert diagram.get("numerator_prefactor_interface") is None
                mapping, reversed_edges = finalization_map(
                    graphs[stage], diagram, case == "representative-lr"
                )
                assert finalized_terms(provisional, mapping, reversed_edges) == terms(
                    diagram["numerator_prefactor"]
                ), (side, case, folder.name, stage)
                counts["complete_effective_tensor_bindings"] += 1
                counts["provisional_to_final_edges"] += len(mapping)
        results.append({"side": side, "case": case, "counts": dict(counts)})
(root / "prefactor-summary.json").write_text(json.dumps(results, indent=2) + "\n")
print(json.dumps(results, indent=2))
