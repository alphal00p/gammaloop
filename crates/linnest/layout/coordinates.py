"""Convert planarization vertices into proper geometric route crossings.

Crossing vertices belong only to the planarized drawing graph. Physical routes
pass through them without gaining a graph vertex or sharing a movable knot.
"""

from __future__ import annotations

import math
from itertools import pairwise


def _point_segment_distance(point, first, second):
    direction = [second[k] - first[k] for k in (0, 1)]
    squared = sum(value * value for value in direction)
    if squared == 0:
        raise ValueError("The planarized drawing has a zero-length edge")
    fraction = max(
        0.0,
        min(
            1.0,
            sum((point[k] - first[k]) * direction[k] for k in (0, 1)) / squared,
        ),
    )
    return math.dist(point, [first[k] + fraction * direction[k] for k in (0, 1)])


def _proper_intersection(first, second, third, fourth):
    def cross(a, b):
        return a[0] * b[1] - a[1] * b[0]

    u = [second[k] - first[k] for k in (0, 1)]
    v = [fourth[k] - third[k] for k in (0, 1)]
    offset = [third[k] - first[k] for k in (0, 1)]
    determinant = cross(u, v)
    if determinant == 0:
        raise ValueError(
            "The two physical routes do not alternate at a crossing (parallel chords)"
        )
    a, b = cross(offset, v) / determinant, cross(offset, u) / determinant
    if not (0 < a < 1 and 0 < b < 1):
        raise ValueError("The two physical routes do not alternate at a crossing")
    return [first[k] + a * u[k] for k in (0, 1)]


def open_crossings(positions, routes, crossing_nodes):
    """Replace each crossing node by two disjoint knots on each passing route.

    The input drawing is planar after treating ``crossing_nodes`` as vertices.
    Each such vertex must occur internally exactly twice among the routes, with
    four distinct neighbors. Its disk is disjoint from every nonincident segment
    and other vertex. Connecting alternate disk-boundary points by chords then
    creates exactly the intended proper crossing and no accidental contacts.

    Returns copied positions/routes and crossing records. Input objects remain
    unchanged, and crossing keys never survive into the returned physical data.
    """
    crossing_nodes = set(crossing_nodes)
    result_positions = {key: list(point) for key, point in positions.items()}
    result_routes = {edge: list(route) for edge, route in routes.items()}
    if not crossing_nodes:
        return result_positions, result_routes, []
    if not crossing_nodes <= positions.keys():
        raise ValueError("Unknown crossing vertex")
    segments = {(a, b) for route in routes.values() for a, b in pairwise(route)}
    replacements = {}
    records = []
    for node in sorted(crossing_nodes, key=str):
        occurrences = [
            (edge, index)
            for edge, route in routes.items()
            for index, key in enumerate(route)
            if key == node
        ]
        if len(occurrences) != 2:
            raise ValueError("Every crossing vertex must belong to two route spans")
        neighbors = []
        for edge, index in occurrences:
            route = routes[edge]
            if index == 0 or index == len(route) - 1:
                raise ValueError("A crossing cannot be a physical route endpoint")
            neighbors.extend((route[index - 1], route[index + 1]))
        if len(set(neighbors)) != 4:
            raise ValueError("A crossing needs four distinct incident neighbors")
        center = positions[node]
        radius = min(
            math.dist(center, point) / 8
            for key, point in positions.items()
            if key != node
        )
        for first, second in segments:
            if node not in (first, second):
                radius = min(
                    radius,
                    _point_segment_distance(center, positions[first], positions[second])
                    / 4,
                )
        if not math.isfinite(radius) or radius <= 0:
            raise ValueError("The crossing has no positive geometric clearance")
        knots = []
        for occurrence, (edge, index) in enumerate(occurrences):
            route = routes[edge]
            pair = []
            for end, neighbor in enumerate((route[index - 1], route[index + 1])):
                direction = [positions[neighbor][k] - center[k] for k in (0, 1)]
                length = math.hypot(*direction)
                key = f"b:{edge}:cross:{node}:{occurrence}:{end}"
                if key in result_positions:
                    raise ValueError("Crossing knot ID collides with an existing point")
                result_positions[key] = [
                    center[k] + radius * direction[k] / length for k in (0, 1)
                ]
                pair.append(key)
            replacements[edge, index] = pair
            knots.append(pair)
        point = _proper_intersection(
            *(result_positions[key] for pair in knots for key in pair)
        )
        records.append(
            {
                "id": node,
                "edges": [edge for edge, _ in occurrences],
                "point": point,
                "route_points": knots,
                "clearance_radius": radius,
            }
        )
    for edge, route in routes.items():
        result_routes[edge] = [
            key
            for index, node in enumerate(route)
            for key in replacements.get((edge, index), [node])
        ]
    for node in crossing_nodes:
        del result_positions[node]
    return result_positions, result_routes, records
