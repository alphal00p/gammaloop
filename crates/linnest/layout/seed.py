"""Physical route coordinates from an EC-constrained planarization.

The exterior constraint is a grouping tree with incoming/outgoing child groups.
It encodes consecutive exterior blocks, not metric X coordinates. The latter
are realized separately by trimming protected side-hub spokes at vertical rails.
"""

from __future__ import annotations

import math
import statistics
from collections import Counter
from itertools import pairwise

import networkx as nx
from coordinates import open_crossings
from ec_planarity import ECExpansion, Graph
from geometry import (
    draw_with_outer_edge,
    normalized_edges,
    spread_external_fans,
    validate_geometry,
)
from planarization import planarize


def _group(children):
    return (
        children[0] if len(children) == 1 else {"kind": "group", "children": children}
    )


def _side_hubs(graph, expansion, exterior, sides):
    """Keep two metric hubs even when a constraint side is only one leaf."""
    graph = graph.copy()
    ports = expansion.leaf_ports[exterior]
    hubs = {}
    for side, spokes in sides.items():
        node = ports[spokes[0]]
        if len(spokes) == 1:
            spoke = spokes[0]
            owner = graph.other(spoke, node)
            hub, reference = f"__metric_{side}_hub__", f"__metric_{side}_reference__"
            if hub in graph.nodes or reference in graph.edges:
                raise ValueError("Reserved side-hub identity is already in use")
            graph.nodes.add(hub)
            graph.edges[spoke] = (owner, hub)
            graph.edges[reference] = (hub, node)
            graph.rotation[node][graph.rotation[node].index(spoke)] = reference
            graph.rotation[hub] = [spoke, reference]
            graph.protected.add(reference)
            node = hub
        elif any(ports[spoke] != node for spoke in spokes):
            raise ArithmeticError("One side did not expand to one grouping node")
        hubs[side] = node
    centers = expansion.members[exterior] - set(hubs.values())
    if len(centers) != 1:
        raise ArithmeticError("Exterior grouping scaffold has no unique center")
    center = centers.pop()
    references = graph.rotation[center]
    if len(references) != 2 or any(edge not in graph.protected for edge in references):
        raise ArithmeticError("Exterior center is not a protected degree-two vertex")
    left, right = hubs["incoming"], hubs["outgoing"]
    bridge = "__metric_exterior_bridge__"
    for edge in references:
        other = graph.other(edge, center)
        graph.rotation[other][graph.rotation[other].index(edge)] = bridge
        del graph.edges[edge]
        graph.protected.remove(edge)
    del graph.rotation[center]
    graph.nodes.remove(center)
    graph.edges[bridge] = (left, right)
    graph.protected.add(bridge)
    graph.validate_embedding()
    return graph, hubs, bridge


def _canonical_rotation(order):
    if not order:
        return []
    start = order.index(min(order))
    return order[start:] + order[:start]


def _draw(graph, outer=None):
    """Pass the certified edge-identity rotation to the grid drawing routine."""
    # NetworkX's canonical ordering selects from a set. Integer drawing IDs
    # keep that choice independent of Python's randomized string hashes.
    name_of = dict(enumerate(sorted(graph.nodes)))
    id_of = {name: index for index, name in name_of.items()}
    simple = nx.Graph()
    simple.add_nodes_from(name_of)
    simple.add_edges_from(
        tuple(id_of[node] for node in graph.edges[edge]) for edge in sorted(graph.edges)
    )
    if simple.number_of_edges() != len(graph.edges):
        raise ValueError("Subdivide parallel drawing edges before grid drawing")
    embedding = nx.PlanarEmbedding()
    embedding.set_data(
        {
            id_of[node]: [
                id_of[graph.other(edge, node)]
                for edge in _canonical_rotation(graph.rotation[node])
            ]
            for node in sorted(graph.nodes)
        }
    )
    embedding.add_nodes_from(name_of)
    embedding.check_structure()
    if outer:
        outer = tuple(id_of[node] for node in outer)
    components = sorted(
        nx.connected_components(simple),
        key=lambda nodes: (bool(outer and outer[0] not in nodes), min(nodes)),
    )
    positions = {}
    top = 0.0
    for component in components:
        local = nx.PlanarEmbedding()
        local.set_data(
            {
                node: list(embedding.neighbors_cw_order(node))
                for node in sorted(component)
            }
        )
        local.add_nodes_from(component)
        if len(component) == 1:
            coords = {next(iter(component)): (0.0, 0.0)}
        elif outer and outer[0] in component:
            coords = draw_with_outer_edge(local, *outer)
        else:
            coords = nx.combinatorial_embedding_to_pos(local)
        if positions:
            # Components without external legs occupy a separate horizontal
            # strip within the first component's X range, above every edge.
            low, high = (
                min(p[0] for p in positions.values()),
                max(p[0] for p in positions.values()),
            )
            minimum, maximum = (
                min(p[0] for p in coords.values()),
                max(p[0] for p in coords.values()),
            )
            factor = (high - low) / (2 * max(1.0, maximum - minimum))
            center = (minimum + maximum) / 2
            coords = {
                node: ((low + high) / 2 + (x - center) * factor, top + 4 + y * factor)
                for node, (x, y) in coords.items()
            }
        positions.update(coords)
        top = max(p[1] for p in positions.values())
    return {name_of[node]: point for node, point in positions.items()}


def initialize(diagram, scale=2.4, *, external_sides=True):
    """Generate an EC seed; crossings come only from explicit edge insertions.

    The returned node/edge identities are physical identities. Route knots and
    planarization crossings are drawing geometry only. This adapter supports
    free incoming/outgoing side groups; arbitrary fixed XY coordinates remain
    the native constraint projector's responsibility.
    """
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("scale must be finite and positive")
    node_ids = {str(i): f"v:{i}" for i in range(len(diagram["vertices"]))}
    nodes = set(node_ids.values())
    edges, routes, external_ids, endpoints, route_segments = {}, {}, {}, {}, {}
    incoming, outgoing, spokes, groups = [], [], {"incoming": [], "outgoing": []}, {}
    exterior = "__ec_exterior__"
    for edge in normalized_edges(diagram):
        identity = str(edge["id"])
        source, target = edge["source"], edge["target"]
        if source is None and target is None:
            raise ValueError(f"Edge {identity} has no physical endpoint")
        for endpoint in (source, target):
            if endpoint is not None and str(endpoint) not in node_ids:
                raise ValueError(f"Unknown vertex {endpoint} on edge {identity}")
        a = f"x:{identity}" if source is None else node_ids[str(source)]
        b = f"x:{identity}" if target is None else node_ids[str(target)]
        endpoints[identity] = {"source": source, "target": target}
        if source is None or target is None:
            endpoint, owner = (a, b) if source is None else (b, a)
            flow = edge.get("state") or ("incoming" if source is None else "outgoing")
            if flow not in spokes:
                raise ValueError(f"Unknown external flow {flow}")
            (incoming if flow == "incoming" else outgoing).append(endpoint)
            external_ids[identity] = endpoint
            routes[identity] = [a, b]
            key = flow, owner
            if key not in groups:
                spoke = f"spoke:{flow}:{owner}"
                groups[key] = {"spoke": spoke, "endpoints": []}
                spokes[flow].append(spoke)
                edges[spoke] = (owner, exterior)
            groups[key]["endpoints"].append(endpoint)
        else:
            bends = 2 if a == b else 1
            route = [a, *(f"b:{identity}:{i}" for i in range(bends)), b]
            routes[identity] = route
            route_segments[identity] = []
            nodes.update(route)
            for index, pair in enumerate(pairwise(route)):
                segment = f"route:{identity}:{index}"
                edges[segment] = pair
                route_segments[identity].append(segment)
    protected = {data["spoke"] for data in groups.values()}
    constraints = {}
    side_mode = external_sides and incoming and outgoing
    if groups:
        nodes.add(exterior)
        leaves = [spoke for side in spokes.values() for spoke in side]
        constraints[exterior] = (
            _group([_group(side) for side in spokes.values() if side])
            if side_mode
            else _group(leaves)
        )
    expansion = ECExpansion(Graph(nodes, edges, protected=protected), constraints)
    result = planarize(expansion.graph)
    embedded = result["embedding"]
    if side_mode:
        embedded, hubs, _ = _side_hubs(embedded, expansion, exterior, spokes)
        positions = _draw(embedded, (hubs["incoming"], hubs["outgoing"]))
    elif groups:
        hub = next(iter(expansion.members[exterior]))
        hubs = {"radial": hub}
        positions = _draw(embedded, (hub,))
    else:
        hubs = {}
        positions = _draw(embedded)
    for identity, segments in route_segments.items():
        current = routes[identity][0]
        route = [current]
        for segment in segments:
            for part in result["edge_chains"][segment]:
                current = embedded.other(part, current)
                route.append(current)
        if current != routes[identity][-1]:
            raise ArithmeticError("Planarization changed a physical route endpoint")
        routes[identity] = route
    hub_positions = {side: positions.pop(hub) for side, hub in hubs.items()}
    if side_mode:
        left, right = hub_positions["incoming"], hub_positions["outgoing"]
        low, high = (
            min(p[0] for p in positions.values()),
            max(p[0] for p in positions.values()),
        )
        if not left[0] < low <= high < right[0]:
            raise ArithmeticError("The exterior hub edge did not bound the drawing")
        rails = {"incoming": (left[0] + low) / 2, "outgoing": (right[0] + high) / 2}
    for (flow, owner), data in groups.items():
        point = positions[owner]
        hub = hub_positions[flow] if side_mode else hub_positions["radial"]
        fraction = (rails[flow] - hub[0]) / (point[0] - hub[0]) if side_mode else 0.15
        positions[data["endpoints"][0]] = [
            hub[k] + fraction * (point[k] - hub[k]) for k in (0, 1)
        ]
    spread_external_fans(
        positions, routes, {key: data["endpoints"] for key, data in groups.items()}
    )
    positions, routes, crossings = open_crossings(
        positions, routes, result["crossing_nodes"]
    )
    center = (
        [statistics.mean(point[k] for point in positions.values()) for k in (0, 1)]
        if positions
        else [0, 0]
    )
    lengths = [
        sum(math.dist(positions[a], positions[b]) for a, b in pairwise(route))
        for route in routes.values()
    ]
    factor = scale / statistics.median(lengths) if lengths else 1.0
    positions = {
        node: [factor * (point[k] - center[k]) for k in (0, 1)]
        for node, point in positions.items()
    }
    for crossing in crossings:
        crossing["point"] = [
            factor * (crossing["point"][k] - center[k]) for k in (0, 1)
        ]
        crossing["clearance_radius"] *= factor
    geometry = validate_geometry(positions, routes)
    if any(geometry[key] for key in ("contacts", "overlaps", "degenerate")):
        raise ArithmeticError(f"EC coordinates have a geometric degeneracy: {geometry}")
    actual = Counter(
        tuple(sorted((first[0], second[0]))) for first, second in geometry["crossings"]
    )
    expected = Counter(tuple(sorted(record["edges"])) for record in crossings)
    if actual != expected:
        raise ArithmeticError(
            "Geometric crossings disagree with EC insertion provenance"
        )
    return {
        "positions": positions,
        "routes": routes,
        "node_ids": node_ids,
        "external_ids": external_ids,
        "incoming_ids": incoming,
        "outgoing_ids": outgoing,
        "edge_endpoints": endpoints,
        "crossings": crossings,
        "report": {
            "status": "ec-planarized" if crossings else "ec-planar",
            "method": "Gutwenger–Klein–Mutzel EC expansion and optimal individual edge insertion; Chrobak–Payne coordinates",
            "constraint_tree": constraints.get(exterior),
            "constraint_scope": "consecutive exterior incoming/outgoing groups; metric rails constructed separately"
            if side_mode
            else "unconstrained radial/cofacial external order",
            "external_routes_straight": True,
            "external_side_constraints": {
                "requested": bool(external_sides),
                "realized": bool(side_mode),
            },
            "insertions": result["insertions"],
            "geometry": geometry,
        },
    }
