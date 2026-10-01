"""Run EC embedding or optimal insertion from a topology-only JSON problem."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from ec_planarity import ECExpansion, Graph, NotECPlanar, embed_expanded_graph
from insertion import optimal_insertion, planarize_insertion
from planarization import planarize


def graph_json(graph):
    return {
        "nodes": sorted(graph.nodes),
        "edges": graph.edges,
        "rotation_cw": graph.rotation,
        "protected": sorted(graph.protected),
        "oriented_hubs": graph.hubs,
        "wheel_hubs": sorted(graph.wheel_hubs),
    }


def run(problem, mode, inserted=None):
    graph = Graph(
        set(problem["nodes"]),
        problem["edges"],
        protected=set(problem.get("protected", [])),
    )
    expansion = ECExpansion(graph, problem.get("constraints", {}))
    if mode == "embed":
        try:
            embedded = expansion.embed()
        except NotECPlanar as error:
            return {"ec_planar": False, "reason": str(error)}
        return {"ec_planar": True, "embedding": graph_json(embedded)}
    if mode == "insert":
        if inserted not in graph.edges:
            raise ValueError("--insert must name an edge of the full input graph")
        if inserted in graph.protected:
            raise ValueError("The edge selected for insertion is protected")
        # Testing zero first also handles inserting a bridge between components.
        # The positive-cost branch uses E(G+e,C) minus e, with all leaf ports kept.
        try:
            embedded = expansion.embed_expanded()
        except NotECPlanar:
            base = expansion.graph.copy()
            source, target = base.edges.pop(inserted)
            base.rotation[source].remove(inserted)
            base.rotation[target].remove(inserted)
            base = embed_expanded_graph(base)
            path = optimal_insertion(base, source, target)
            result = planarize_insertion(path, inserted)
            result["states"] = path["states"]
        else:
            result = {
                "embedding": embedded,
                "edge_chains": {edge: [edge] for edge in embedded.edges},
                "crossing_nodes": {},
                "crossings": [],
                "cost": 0,
                "states": [],
            }
    else:
        result = planarize(expansion.graph)
    result["embedding"].validate_embedding()
    return {
        **result,
        "embedding": graph_json(result["embedding"]),
        "expansion_members": {
            node: sorted(members) for node, members in expansion.members.items()
        },
        "leaf_ports": expansion.leaf_ports,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input", type=Path, help="JSON nodes/edges/constraints/protected"
    )
    parser.add_argument(
        "--mode", choices=("embed", "insert", "planarize"), default="embed"
    )
    parser.add_argument(
        "--insert", help="Original edge removed and optimally reinserted"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.insert is not None and args.mode != "insert":
        parser.error("--insert is only used with --mode insert")
    try:
        result = run(json.loads(args.input.read_text()), args.mode, args.insert)
    except (ValueError, RuntimeError) as error:
        print(str(error), file=sys.stderr)
        return 1
    encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(encoded)
    else:
        sys.stdout.write(encoded)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
