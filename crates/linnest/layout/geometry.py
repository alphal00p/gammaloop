"""Shared identity and geometric checks for planar layout initializers."""

from __future__ import annotations

import math
from itertools import pairwise

import networkx as nx


def normalized_edges(diagram):
    """Read native FeynmanDiagram JSON or explicit {id,source,target} edges."""
    edges = []
    for index, record in enumerate(diagram["edges"]):
        if isinstance(record, list):
            endpoints, metadata = record
            external = metadata.get("external") or {}
            edge = {"id": index, **endpoints, "state": external.get("state")}
        else:
            edge = dict(record)
            edge.setdefault("id", index)
        edges.append(edge)
    if len({str(edge["id"]) for edge in edges}) != len(edges):
        raise ValueError("Physical edge IDs must be unique")
    return edges


def validate_geometry(positions, routes, tolerance=1e-8):
    """Check proper crossings, overlaps, contacts and zero-length segments.

    Shared graph endpoints are allowed; overlapping incident edges are not.
    This check is independent of NetworkX's combinatorial planarity result.
    """
    segments = [
        (edge, index, a, b)
        for edge, path in routes.items()
        for index, (a, b) in enumerate(pairwise(path))
    ]
    crossings, contacts, overlaps, degenerate = [], [], [], []

    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def on_segment(a, b, point):
        return abs(orient(a, b, point)) <= tolerance and all(
            min(a[k], b[k]) - tolerance <= point[k] <= max(a[k], b[k]) + tolerance
            for k in (0, 1)
        )

    for i, first in enumerate(segments):
        edge, index, aid, bid = first
        a, b = positions[aid], positions[bid]
        if math.dist(a, b) <= tolerance:
            degenerate.append([edge, index])
        for second in segments[i + 1 :]:
            other, j, cid, did = second
            c, d = positions[cid], positions[did]
            pair = [[edge, index], [other, j]]
            values = [
                orient(a, b, c),
                orient(a, b, d),
                orient(c, d, a),
                orient(c, d, b),
            ]
            if max(abs(value) for value in values) <= tolerance:
                axis = 0 if abs(a[0] - b[0]) >= abs(a[1] - b[1]) else 1
                extent = min(max(a[axis], b[axis]), max(c[axis], d[axis])) - max(
                    min(a[axis], b[axis]), min(c[axis], d[axis])
                )
                if extent > tolerance:
                    overlaps.append(pair)
                    continue
            signs = [
                0 if abs(value) <= tolerance else math.copysign(1, value)
                for value in values
            ]
            if signs[0] * signs[1] < 0 and signs[2] * signs[3] < 0:
                crossings.append(pair)
            elif not ({aid, bid} & {cid, did}) and (
                on_segment(a, b, c)
                or on_segment(a, b, d)
                or on_segment(c, d, a)
                or on_segment(c, d, b)
            ):
                contacts.append(pair)
    return {
        "valid": not (crossings or contacts or overlaps or degenerate),
        "crossings": crossings,
        "contacts": contacts,
        "overlaps": overlaps,
        "degenerate": degenerate,
        "segments": len(segments),
    }


def draw_with_outer_edge(embedding, first, second=None):
    """Keep the chosen edge on the outer triangle of the integer-grid drawing."""
    from types import FunctionType

    if len(embedding) == 3 and second is not None:
        # NetworkX returns a fixed triangle before consulting its outer-face
        # callback for this smallest case. Preserve the requested hub edge.
        third = next(node for node in embedding if node not in (first, second))
        return {first: (0, 0), second: (2, 0), third: (1, 1)}
    drawing = nx.combinatorial_embedding_to_pos
    triangulate = drawing.__globals__["triangulate_embedding"]

    def triangulate_with_exterior(source, fully_triangulate):
        triangulated, _ = triangulate(source, fully_triangulate=True)
        neighbor = second
        if neighbor is None:
            neighbor = next(triangulated.neighbors_cw_order(first))
        return triangulated, triangulated.traverse_face(first, neighbor)

    local_drawing = FunctionType(
        drawing.__code__,
        {**drawing.__globals__, "triangulate_embedding": triangulate_with_exterior},
        drawing.__name__,
        drawing.__defaults__,
        drawing.__closure__,
    )
    return local_drawing(embedding)


def spread_external_fans(positions, routes, groups):
    """Replace each repeated auxiliary spoke by a narrow straight planar fan.

    Moving a spoke's rail endpoint by epsilon moves its whole segment by at
    most epsilon. Nonincident segment/vertex clearance bounds protect those
    obstacles, and determinant bounds protect the angles at shared owners.
    The same bound applies simultaneously to every fan, so two moving spokes
    cannot consume each other's clearance. No physical cyclic order is fixed
    by this initializer; duplicate external legs occupy the spoke's free sector.
    """
    repeated = [endpoints for endpoints in groups.values() if len(endpoints) > 1]
    if not repeated:
        return
    # Only one representative per coincident auxiliary spoke exists yet.
    segments = [
        (a, b)
        for route in routes.values()
        if all(point in positions for point in route)
        for a, b in pairwise(route)
    ]
    radius = min(math.dist(positions[a], positions[b]) for a, b in segments) / 4

    def distance(point, a, b):
        direction = [b[k] - a[k] for k in (0, 1)]
        squared = sum(value * value for value in direction)
        fraction = max(
            0.0,
            min(1.0, sum((point[k] - a[k]) * direction[k] for k in (0, 1)) / squared),
        )
        return math.dist(point, [a[k] + fraction * direction[k] for k in (0, 1)])

    for (_, owner), endpoints in groups.items():
        if len(endpoints) == 1:
            continue
        endpoint = endpoints[0]
        a, b = positions[owner], positions[endpoint]
        u = [b[k] - a[k] for k in (0, 1)]
        for key, point in positions.items():
            if key not in (owner, endpoint):
                radius = min(radius, distance(point, a, b) / 4)
        for first, second in segments:
            if {first, second} == {owner, endpoint}:
                continue
            c, d = positions[first], positions[second]
            if owner in (first, second):
                other = d if first == owner else c
                w = [other[k] - a[k] for k in (0, 1)]
                determinant = u[0] * w[1] - u[1] * w[0]
                if abs(determinant) <= 1e-12 * math.hypot(*u) * math.hypot(*w):
                    if sum(u[k] * w[k] for k in (0, 1)) >= 0:
                        raise ArithmeticError("Auxiliary spokes overlap at their owner")
                    # Opposite rays cannot overlap while their dot product is
                    # negative. Both endpoints may move vertically by radius.
                    radius = min(radius, math.hypot(*u) / 4, math.hypot(*w) / 4)
                else:
                    # Both incident fan endpoints can move; their fixed X
                    # components make this determinant bound linear and exact.
                    radius = min(
                        radius, abs(determinant) / (4 * (abs(u[0]) + abs(w[0])))
                    )
            else:
                separation = min(
                    distance(a, c, d),
                    distance(b, c, d),
                    distance(c, a, b),
                    distance(d, a, b),
                )
                radius = min(radius, separation / 4)
    if not math.isfinite(radius) or radius <= 0:
        raise ArithmeticError("No positive planar clearance for external fans")
    for endpoints in repeated:
        x, y = positions[endpoints[0]]
        for index, endpoint in enumerate(endpoints):
            positions[endpoint] = (
                x,
                y + radius * (2 * index / (len(endpoints) - 1) - 1),
            )
