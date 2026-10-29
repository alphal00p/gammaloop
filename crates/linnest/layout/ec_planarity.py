"""Embedding constraints of Gutwenger, Klein and Mutzel (JGAA 12(1), 2008).

The constraint expansion and O-hub test implement Sections 3--5.  OGDF is used
only for its planar SPQR decomposition, not as an unconstrained edge inserter.
Rotations are clockwise and contain edge identities, so parallel edges remain
distinct. Self-loops must first be subdivided by the drawing adapter.
"""

from __future__ import annotations

import copy
import json
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path


class NotECPlanar(ValueError):
    """No planar rotation satisfies all of the supplied constraint trees."""


def cyclic_equal(first, second):
    if len(first) != len(second):
        return False
    if not first:
        return True
    try:
        start = second.index(first[0])
    except ValueError:
        return False
    return all(
        value == second[(start + i) % len(second)] for i, value in enumerate(first)
    )


@dataclass
class Graph:
    nodes: set[str] = field(default_factory=set)
    edges: dict[str, tuple[str, str]] = field(default_factory=dict)
    rotation: dict[str, list[str]] = field(default_factory=dict)
    protected: set[str] = field(default_factory=set)
    hubs: dict[str, list[str]] = field(default_factory=dict)
    wheel_hubs: set[str] = field(default_factory=set)

    def __post_init__(self):
        self.nodes = set(self.nodes)
        self.edges = {edge: tuple(ends) for edge, ends in self.edges.items()}
        if any(not isinstance(key, str) for key in (*self.nodes, *self.edges)):
            raise ValueError("Vertex and edge identities must be strings")
        for edge, ends in self.edges.items():
            if len(ends) != 2 or any(node not in self.nodes for node in ends):
                raise ValueError(f"Invalid endpoints for edge {edge}")
            if ends[0] == ends[1]:
                raise ValueError("Subdivide self-loops before EC embedding")
        if not self.protected <= self.edges.keys():
            raise ValueError("Protected edges must belong to the graph")
        if not self.rotation:
            self.rotation = {node: [] for node in self.nodes}
            for edge, (a, b) in self.edges.items():
                self.rotation[a].append(edge)
                self.rotation[b].append(edge)

    def copy(self):
        return copy.deepcopy(self)

    def other(self, edge, node):
        a, b = self.edges[edge]
        if node == a:
            return b
        if node == b:
            return a
        raise ValueError(f"Vertex {node} is not incident to edge {edge}")

    def subgraph(self, edge_ids):
        chosen = set(edge_ids)
        edges = {edge: ends for edge, ends in self.edges.items() if edge in chosen}
        nodes = {node for ends in edges.values() for node in ends}
        return Graph(
            nodes,
            edges,
            {
                node: [edge for edge in self.rotation[node] if edge in chosen]
                for node in nodes
            },
            self.protected & chosen,
            {node: list(order) for node, order in self.hubs.items() if node in nodes},
            self.wheel_hubs & nodes,
        )

    def faces(self):
        """Return right-face boundaries and the face ID of every directed dart.

        A dart is ``(edge_id, origin_vertex)``.  Walking its right face takes
        the preceding edge in the clockwise rotation at its target.
        """
        previous = {}
        incidences = {node: set() for node in self.nodes}
        for edge, ends in self.edges.items():
            for node in ends:
                incidences[node].add(edge)
        for node in self.nodes:
            order = self.rotation.get(node, [])
            if len(order) != len(set(order)) or set(order) != incidences[node]:
                raise ValueError(f"Rotation at {node} is not its incidence list")
            for index, edge in enumerate(order):
                previous[(edge, node)] = order[index - 1]
        boundaries, face_of = [], {}
        for edge, ends in self.edges.items():
            for node in ends:
                start = (edge, node)
                if start in face_of:
                    continue
                boundary = []
                dart = start
                while dart not in face_of:
                    face_of[dart] = len(boundaries)
                    boundary.append(dart)
                    target = self.other(*dart)
                    dart = (previous[(dart[0], target)], target)
                if dart != start:
                    raise ValueError("Invalid face permutation")
                boundaries.append(boundary)
        return boundaries, face_of

    def validate_embedding(self):
        """Certify Euler planarity component by component and O-hub orientation."""
        boundaries, _ = self.faces()
        adjacency = {node: [] for node in self.nodes}
        for a, b in self.edges.values():
            adjacency[a].append(b)
            adjacency[b].append(a)
        seen = set()
        for root in self.nodes:
            if root in seen:
                continue
            pending, component = [root], set()
            while pending:
                node = pending.pop()
                if node in component:
                    continue
                component.add(node)
                pending.extend(adjacency[node])
            seen.update(component)
            edge_count = sum(a in component for a, _ in self.edges.values())
            face_count = sum(
                bool(face) and face[0][1] in component for face in boundaries
            )
            if edge_count and len(component) - edge_count + face_count != 2:
                raise NotECPlanar("Rotation has positive genus")
        for hub, expected in self.hubs.items():
            if not cyclic_equal(self.rotation[hub], expected):
                raise NotECPlanar(f"Incorrectly oriented O-hub {hub}")
        for boundary in boundaries:
            if any(node in self.wheel_hubs for _, node in boundary) and (
                len(boundary) != 3
                or any(edge not in self.protected for edge, _ in boundary)
            ):
                raise NotECPlanar("A wheel interior contains non-gadget edges")
        return self


def decompose(graph):
    """Obtain the embedded S/P/R skeletons of one biconnected graph."""
    binary = Path(os.environ.get("EC_SPQR_BINARY") or shutil.which("ec-spqr") or "ec-spqr")
    request = {
        "nodes": sorted(graph.nodes),
        "edges": [
            {"id": edge, "source": a, "target": b}
            for edge, (a, b) in graph.edges.items()
        ],
    }
    if not binary.is_file():
        raise RuntimeError(f"Build the pinned SPQR helper first: {binary}")
    process = subprocess.run(
        [str(binary)],
        input=json.dumps(request) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    result = json.loads(process.stdout)
    if result.get("error"):
        raise RuntimeError(result["error"])
    if not result["planar"]:
        raise NotECPlanar("The constraint expansion is nonplanar")
    if not result.get("biconnected", True):
        raise ValueError("SPQR decomposition requires a biconnected block")
    return result


def feasible_skeletons(graph, spqr):
    """Algorithm 1: orient each R skeleton or reject conflicting O-hubs."""
    from insertion import skeleton_edge_id

    skeletons = {}
    for component in spqr["components"]:
        local_names = {
            edge["id"]: edge["real_edge"]
            if edge["real_edge"] is not None
            else skeleton_edge_id(component["id"], edge["id"], graph.edges)
            for edge in component["edges"]
        }
        edges = {
            local_names[edge["id"]]: graph.edges[edge["real_edge"]]
            if edge["real_edge"] is not None
            else tuple(sorted((edge["source"], edge["target"])))
            for edge in component["edges"]
        }
        original_to_local = {
            edge["real_edge"]: local_names[edge["id"]]
            for edge in component["edges"]
            if edge["real_edge"] is not None
        }
        hubs = {
            hub: [original_to_local[edge] for edge in expected]
            for hub, expected in graph.hubs.items()
            if hub in component["nodes"]
        }
        skeleton = Graph(
            set(component["nodes"]),
            edges,
            {
                node: [local_names[edge] for edge in order]
                for node, order in component["rotation_cw"].items()
            },
            {
                original_to_local[edge]
                for edge in graph.protected
                if edge in original_to_local
            },
            hubs,
            graph.wheel_hubs & set(component["nodes"]),
        )
        orientations = set()
        for hub, expected in hubs.items():
            if component["type"] != "R":
                raise ArithmeticError("Wheel hub is not contained in an R skeleton")
            actual = skeleton.rotation[hub]
            if cyclic_equal(actual, expected):
                orientations.add(False)
            elif cyclic_equal(actual[::-1], expected):
                orientations.add(True)
            else:
                raise ArithmeticError(
                    "SPQR decomposition changed a wheel's cyclic order"
                )
        if len(orientations) > 1:
            raise NotECPlanar("An R skeleton contains conflicting O-hubs")
        if orientations == {True}:
            skeleton.rotation = {
                node: order[::-1] for node, order in skeleton.rotation.items()
            }
        skeletons[component["id"]] = skeleton
    return skeletons


def embed_expanded_graph(graph):
    """Algorithm 1 on an expansion, including independently oriented blocks."""
    from insertion import blocks, compose, glue_blocks

    embedded = []
    for edge_ids in blocks(graph):
        block = graph.subgraph(edge_ids)
        if len(block.edges) <= 2:
            embedded.append(block)
        else:
            spqr = decompose(block)
            embedded.append(compose(spqr, feasible_skeletons(block, spqr)))
    return glue_blocks(graph, embedded).validate_embedding()


class ECExpansion:
    """Section 4 wheel expansion, with explicit original incidence provenance."""

    def __init__(self, graph, constraints):
        self.original = graph.copy()
        self.constraints = copy.deepcopy(constraints)
        self.members = {node: {node} for node in graph.nodes}
        self.leaf_ports = {node: {} for node in graph.nodes}
        self.graph = Graph(set(graph.nodes) - constraints.keys(), {})
        if set(constraints) & graph.wheel_hubs:
            raise ValueError("Cannot constrain an existing wheel hub")
        self.graph.hubs = copy.deepcopy(graph.hubs)
        self.graph.wheel_hubs = set(graph.wheel_hubs)
        self._counter = 0
        self._reserved = set(graph.nodes) | set(graph.edges)
        if not constraints.keys() <= graph.nodes:
            raise ValueError("Constraint refers to an unknown vertex")
        for node, tree in constraints.items():
            incident = {edge for edge, ends in graph.edges.items() if node in ends}
            leaves = self._leaves(tree)
            if len(leaves) != len(set(leaves)) or set(leaves) != incident:
                raise ValueError(
                    f"Constraint at {node} must contain every incident edge exactly once"
                )
            self.members[node] = set()
            self._expand(node, tree)
        for edge, (a, b) in graph.edges.items():
            self._edge(
                edge, self.leaf_ports[a].get(edge, a), self.leaf_ports[b].get(edge, b)
            )
        self.graph.protected.update(graph.protected)

    @staticmethod
    def _leaves(tree):
        if isinstance(tree, str):
            return [tree]
        if not isinstance(tree, dict) or set(tree) != {"kind", "children"}:
            raise ValueError("A constraint is an edge ID or a kind/children tree")
        if tree["kind"] not in {"group", "mirror", "oriented"}:
            raise ValueError("Unknown embedding constraint kind")
        if not isinstance(tree["children"], list) or len(tree["children"]) < 2:
            raise ValueError("Constraint nodes must have at least two children")
        return [
            leaf for child in tree["children"] for leaf in ECExpansion._leaves(child)
        ]

    def _fresh(self, prefix):
        while True:
            self._counter += 1
            value = f"__ec_{prefix}_{self._counter}"
            if value not in self._reserved:
                self._reserved.add(value)
                return value

    def _node(self, owner, prefix):
        node = self._fresh(prefix)
        self.graph.nodes.add(node)
        self.graph.rotation[node] = []
        self.members[owner].add(node)
        return node

    def _edge(self, edge, a, b, protected=False):
        self.graph.edges[edge] = (a, b)
        self.graph.rotation[a].append(edge)
        self.graph.rotation[b].append(edge)
        if protected:
            self.graph.protected.add(edge)
        return edge

    def _expand(self, owner, tree, parent=None):
        if isinstance(tree, str):
            # A degree-one constraint carries no cyclic restriction.
            point = self._node(owner, "leaf")
            self.leaf_ports[owner][tree] = point
            return point
        children = tree["children"]
        degree = len(children) + (parent is not None)
        if tree["kind"] == "group":
            center = self._node(owner, "group")
            ports = [center] * degree
        else:
            hub = self._node(owner, "hub")
            self.graph.wheel_hubs.add(hub)
            rim = [self._node(owner, "rim") for _ in range(2 * degree)]
            spokes = []
            for index, point in enumerate(rim):
                spokes.append(self._edge(self._fresh("spoke"), hub, point, True))
                self._edge(
                    self._fresh("rimedge"), point, rim[(index + 1) % len(rim)], True
                )
            if tree["kind"] == "oriented":
                self.graph.hubs[hub] = spokes
            ports = rim[::2]
        if parent is not None:
            self._edge(self._fresh("treeedge"), parent, ports[-1], True)
        for index, child in enumerate(children):
            if isinstance(child, str):
                self.leaf_ports[owner][child] = ports[index]
            else:
                self._expand(owner, child, ports[index])

    def embed_expanded(self):
        return embed_expanded_graph(self.graph)

    def collapse(self, embedded):
        """Contract gadgets while retaining subsequent edge splits and insertions."""
        owner_of = {
            member: owner
            for owner, members in self.members.items()
            for member in members
        }
        nodes = {owner_of.get(node, node) for node in embedded.nodes}
        edges = {
            edge: (owner_of.get(a, a), owner_of.get(b, b))
            for edge, (a, b) in embedded.edges.items()
            if owner_of.get(a, a) != owner_of.get(b, b)
        }
        result = Graph(
            nodes,
            edges,
            protected=embedded.protected & edges.keys(),
            hubs=copy.deepcopy(self.original.hubs),
            wheel_hubs=set(self.original.wheel_hubs),
        )
        for node in embedded.nodes - owner_of.keys():
            result.rotation[node] = list(embedded.rotation[node])
        # Walking the boundary of a regular neighbourhood of each expansion
        # visits all current incident darts in their induced clockwise order.
        for owner in self.original.nodes:
            members = self.members[owner]
            boundary = []
            internal = {
                edge
                for edge, (a, b) in embedded.edges.items()
                if a in members and b in members
            }
            incident = [edge for edge, ends in edges.items() if owner in ends]
            if not incident:
                result.rotation[owner] = []
                continue
            start = incident[0]
            node = next(node for node in embedded.edges[start] if node in members)
            edge = start
            visited = set()
            while (edge, node) not in visited:
                visited.add((edge, node))
                order = embedded.rotation[node]
                edge = order[(order.index(edge) + 1) % len(order)]
                if edge in internal:
                    node = embedded.other(edge, node)
                else:
                    boundary.append(edge)
            if len(boundary) != len(incident) or set(boundary) != set(incident):
                raise ArithmeticError(
                    "Expansion boundary did not recover original incidences"
                )
            result.rotation[owner] = boundary
        result.validate_embedding()
        for node, tree in self.constraints.items():
            if set(result.rotation[node]) == set(self._leaves(tree)) and not admissible(
                tree, result.rotation[node]
            ):
                raise ArithmeticError(
                    "Collapsed embedding violates its constraint tree"
                )
        return result

    def embed(self):
        return self.collapse(self.embed_expanded())


def admissible(tree, rotation):
    """Check a cyclic order against a constraint tree without enumerating orders."""
    expected = ECExpansion._leaves(tree)
    if len(rotation) != len(expected) or set(rotation) != set(expected):
        return False

    def linear(node, order):
        if isinstance(node, str):
            return order == [node]
        children = node["children"]
        leaf_child = {
            leaf: i
            for i, child in enumerate(children)
            for leaf in ECExpansion._leaves(child)
        }
        runs = []
        for leaf in order:
            child = leaf_child[leaf]
            if not runs or runs[-1][0] != child:
                runs.append((child, []))
            runs[-1][1].append(leaf)
        indices = [index for index, _ in runs]
        if len(indices) != len(children) or len(set(indices)) != len(children):
            return False
        forward = list(range(len(children)))
        if node["kind"] == "oriented" and indices != forward:
            return False
        if node["kind"] == "mirror" and indices not in (forward, forward[::-1]):
            return False
        return all(linear(children[index], group) for index, group in runs)

    return any(linear(tree, rotation[i:] + rotation[:i]) for i in range(len(rotation)))
