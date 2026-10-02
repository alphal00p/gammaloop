"""Greedy EC-planar subgraph followed by individually optimal edge insertions.

The insertion algorithm minimizes each edge's crossings over all admissible
embeddings of the current graph. The sequence is a heuristic for the complete
graph's crossing number; it is not a global crossing-minimum algorithm.
"""

from __future__ import annotations

from ec_planarity import ECExpansion, NotECPlanar, embed_expanded_graph
from insertion import optimal_insertion, planarize_insertion


def _subgraph(graph, chosen):
    result = graph.copy()
    result.edges = {edge: ends for edge, ends in graph.edges.items() if edge in chosen}
    result.rotation = {
        node: [edge for edge in order if edge in chosen]
        for node, order in graph.rotation.items()
    }
    result.protected &= chosen
    return result


def planarize(graph):
    """Return a planarized expansion with complete original-edge provenance.

    Protected gadget/scaffold edges are never removed or crossed. All original
    expansion vertices survive the subgraph phase, including reserved ports of
    temporarily deleted edges. Existing crossings become mirror-wheel gadgets
    during subsequent insertions, retaining their alternating incidences.
    """
    try:
        embedded = embed_expanded_graph(graph)
        return {
            "embedding": embedded,
            "edge_chains": {edge: [edge] for edge in graph.edges},
            "crossing_nodes": {},
            "insertions": [],
        }
    except NotECPlanar:
        pass
    chosen = set(graph.protected)
    embedded = embed_expanded_graph(_subgraph(graph, chosen))
    removed = []
    for edge in graph.edges:
        if edge in chosen:
            continue
        try:
            candidate = embed_expanded_graph(_subgraph(graph, chosen | {edge}))
        except NotECPlanar:
            removed.append(edge)
        else:
            chosen.add(edge)
            embedded = candidate
    chains = {edge: [edge] for edge in chosen}
    crossings, insertions = {}, []
    for edge in removed:
        start, end = graph.edges[edge]
        owners = {
            segment: owner for owner, chain in chains.items() for segment in chain
        }
        constraints = {
            node: {"kind": "mirror", "children": list(embedded.rotation[node])}
            for node in crossings
        }
        preserved = ECExpansion(embedded, constraints)
        # Expand the actual embedding afresh. The DP may change its embedding,
        # but each old crossing's mirror wheel enforces the same alternation.
        constrained = preserved.embed_expanded()
        path = optimal_insertion(constrained, start, end)
        step = planarize_insertion(path, edge)
        embedded = preserved.collapse(step["embedding"])
        for node, record in step["crossing_nodes"].items():
            crossings[node] = {
                "crossed_edge": owners[record["crossed_edge"]],
                "inserted_edge": edge,
            }
        chains = {
            owner: [part for segment in chain for part in step["edge_chains"][segment]]
            for owner, chain in chains.items()
        }
        chains[edge] = step["edge_chains"][edge]
        owners = {
            segment: owner for owner, chain in chains.items() for segment in chain
        }
        for node, record in crossings.items():
            incidence = [owners[segment] for segment in embedded.rotation[node]]
            if (
                len(incidence) != 4
                or incidence[0] != incidence[2]
                or incidence[1] != incidence[3]
                or incidence[0] == incidence[1]
            ):
                raise ArithmeticError(
                    "An insertion changed an existing crossing's alternation"
                )
            if set(incidence) != set(record.values()):
                raise ArithmeticError("Crossing provenance changed during insertion")
        insertions.append(
            {
                "edge": edge,
                "crossed_edges": [
                    record["crossed_edge"]
                    for node, record in crossings.items()
                    if node in step["crossing_nodes"]
                ],
                "crossings": path["cost"],
                "states": path["states"],
            }
        )
    if set(chains) != graph.edges.keys():
        raise ArithmeticError("Planarization lost an original edge")
    return {
        "embedding": embedded.validate_embedding(),
        "edge_chains": chains,
        "crossing_nodes": crossings,
        "insertions": insertions,
    }
