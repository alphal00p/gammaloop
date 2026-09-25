"""The two-sided SPQR edge-insertion dynamic program of GKM (2008), §6.

The graph supplied here is an ec-expansion with a feasible embedding. Expansion
edges (and any additional protected edges) cannot be crossed. We keep both ways
of leaving every pertinent graph until backtracking; choosing a locally best
side loses optimality in the presence of oriented constraints.

Off-path expansions are explicit, as in the paper. For a reversible R skeleton
we rebuild its mirrored embedding while leaving child embeddings oriented. This
is equivalent to the paper's starred subpaths, and avoids mirroring O-hubs in an
attached expansion. Python graph copying means this implementation does not
claim the paper's linear running-time bound.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from itertools import pairwise
from math import inf
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ec_planarity import Graph


class NoInsertionPath(ValueError):
    """No insertion avoids every protected edge."""


def skeleton_edge_id(component, edge, reserved=()):
    candidate = f"@spqr:{component}:{edge}"
    while candidate in reserved:
        candidate = "@" + candidate
    return candidate


def blocks(graph: Graph) -> list[set[str]]:
    """Vertex-biconnected edge sets, including bridges, of a loop-free multigraph."""
    adjacency = {v: [] for v in graph.nodes}
    for edge, (u, v) in graph.edges.items():
        adjacency[u].append((edge, v))
        adjacency[v].append((edge, u))
    discovery, low, stack, result = {}, {}, [], []

    def visit(u, parent_edge=None):
        discovery[u] = low[u] = len(discovery)
        for edge, v in adjacency[u]:
            if edge == parent_edge:
                continue
            if v not in discovery:
                stack.append(edge)
                visit(v, edge)
                low[u] = min(low[u], low[v])
                if low[v] >= discovery[u]:
                    component = set()
                    while stack:
                        current = stack.pop()
                        component.add(current)
                        if current == edge:
                            break
                    result.append(component)
            elif discovery[v] < discovery[u]:
                stack.append(edge)
                low[u] = min(low[u], discovery[v])

    for v in sorted(graph.nodes):
        if v not in discovery:
            visit(v)
    return result


def _after(rotation, edge):
    position = rotation.index(edge) + 1
    return rotation[position:] + rotation[: position - 1]


def splice(parent: Graph, parent_edge: str, child: Graph, child_edge: str) -> Graph:
    """Glue oriented embeddings at matching virtual edges, removing both.

    The child rotation starts immediately after its virtual incidence. Keeping
    its clockwise order (rather than reversing the child) preserves O-hubs.
    """
    if set(parent.edges[parent_edge]) != set(child.edges[child_edge]):
        raise ValueError("virtual edges have different poles")
    if set(parent.edges).intersection(child.edges) - {parent_edge, child_edge}:
        raise ValueError("edge IDs overlap while expanding an SPQR skeleton")
    poles = set(parent.edges[parent_edge])
    if parent.nodes.intersection(child.nodes) != poles:
        raise ValueError("SPQR expansions share vertices other than their poles")
    out = parent.copy()
    out.nodes.update(child.nodes)
    out.edges.pop(parent_edge)
    out.edges.update((e, ends) for e, ends in child.edges.items() if e != child_edge)
    for vertex in child.nodes:
        if vertex in poles:
            position = out.rotation[vertex].index(parent_edge)
            out.rotation[vertex][position : position + 1] = _after(
                child.rotation[vertex], child_edge
            )
        else:
            out.rotation[vertex] = child.rotation[vertex].copy()
    out.protected.discard(parent_edge)
    out.protected.update(child.protected - {child_edge})
    out.hubs.update(child.hubs)
    out.wheel_hubs.update(child.wheel_hubs)
    return out


class _Tree:
    def __init__(self, spqr, skeletons):
        self.raw = {str(c["id"]): c for c in spqr["components"]}
        self.skeletons = {str(k): v for k, v in skeletons.items()}
        self.adjacency = {key: {} for key in self.raw}
        reserved = {
            e["real_edge"]
            for c in self.raw.values()
            for e in c["edges"]
            if e["real_edge"] is not None
        }
        for cid, component in self.raw.items():
            for edge in component["edges"]:
                twin = edge.get("twin")
                if twin is not None:
                    other = str(twin["component"])
                    self.adjacency[cid][other] = (
                        skeleton_edge_id(cid, edge["id"], reserved),
                        skeleton_edge_id(other, twin["edge"], reserved),
                    )
        self.root = str(spqr["root"])

    def expand(self, root, blocked=frozenset(), override=None):
        """Expand a connected subtree; keep boundary virtual edges intact."""
        root = str(root)

        def visit(current, parent):
            out = (
                override
                if current == root and override is not None
                else self.skeletons[current]
            ).copy()
            for neighbor, (edge, twin) in self.adjacency[current].items():
                if neighbor == parent or frozenset((current, neighbor)) in blocked:
                    continue
                out = splice(out, edge, visit(neighbor, current), twin)
            return out

        return visit(root, None)

    def allocation_path(self, start, end):
        starts = sorted(k for k, g in self.skeletons.items() if start in g.nodes)
        ends = {k for k, g in self.skeletons.items() if end in g.nodes}
        predecessors = dict.fromkeys(starts)
        queue = deque(starts)
        while queue:
            current = queue.popleft()
            if current in ends:
                path = []
                while current is not None:
                    path.append(current)
                    current = predecessors[current]
                return path[::-1]
            for neighbor in self.adjacency[current]:
                if neighbor not in predecessors:
                    predecessors[neighbor] = current
                    queue.append(neighbor)
        raise ValueError("insertion endpoint absent from SPQR decomposition")


def compose(spqr, skeletons, root=None):
    """Reconstruct the whole block without changing any skeleton orientation."""
    tree = _Tree(spqr, skeletons)
    return tree.expand(tree.root if root is None else root)


@dataclass
class _Dual:
    graph: Graph

    def __post_init__(self):
        self.faces, self.face_of = self.graph.faces()
        self.arcs = [[] for _ in self.faces]
        for edge, (u, v) in self.graph.edges.items():
            if edge not in self.graph.protected:
                left, right = self.face_of[edge, u], self.face_of[edge, v]
                if left != right:
                    self.arcs[left].append((right, edge))
                    self.arcs[right].append((left, edge))

    def at_vertex(self, vertex):
        return {self.face_of[edge, vertex] for edge in self.graph.rotation[vertex]}

    def side(self, edge, side):
        return {self.face_of[edge, self.graph.edges[edge][side]]}

    def shortest(self, sources, targets):
        predecessors = {f: None for f in sorted(sources)}
        queue = deque(sorted(sources))
        while queue:
            face = queue.popleft()
            if face in targets:
                path = []
                while predecessors[face] is not None:
                    face, edge = predecessors[face]
                    path.append(edge)
                return path[::-1]
            for neighbor, edge in self.arcs[face]:
                if neighbor not in predecessors:
                    predecessors[neighbor] = (face, edge)
                    queue.append(neighbor)
        return None


def _p_embedding(tree, cid, previous, following, side):
    graph = tree.skeletons[cid].copy()
    incoming = tree.adjacency[cid][previous][0]
    outgoing = tree.adjacency[cid][following][0]
    u, v = graph.edges[outgoing]
    middle = [e for e in graph.rotation[u] if e not in (incoming, outgoing)]
    order = (
        [outgoing, incoming, *middle] if side == 0 else [outgoing, *middle, incoming]
    )
    graph.rotation[u], graph.rotation[v] = order, order[::-1]
    return graph


def _block_insertion(graph, start, end):
    # Import lazily: ec_planarity also uses blocks/compose for Algorithm 1.
    from ec_planarity import decompose, feasible_skeletons

    if len(graph.edges) <= 2:
        return {"embedding": graph.copy(), "crossings": [], "cost": 0, "states": []}
    spqr = decompose(graph)
    skeletons = feasible_skeletons(graph, spqr)
    tree = _Tree(spqr, skeletons)
    path = tree.allocation_path(start, end)
    boundaries = {frozenset(pair) for pair in pairwise(path)}
    costs = [0, 0]
    choices, history = [], []
    for index, cid in enumerate(path):
        previous = path[index - 1] if index else None
        following = path[index + 1] if index + 1 < len(path) else None
        kind = tree.raw[cid]["type"]
        if kind == "P" and len(path) > 1:
            if previous is None or following is None:
                raise ValueError("a shortest allocation path cannot end at a P node")
            local = []
            for side in (0, 1):
                skeleton = _p_embedding(tree, cid, previous, following, side)
                graph_i = tree.expand(cid, boundaries, skeleton)
                local.append((side, [], graph_i, False))
            choices.append(local)
            history.append(
                {
                    "component": cid,
                    "type": kind,
                    "costs": [None if cost == inf else cost for cost in costs],
                }
            )
            continue
        candidates = [[], []]
        for mirrored in (False, True):
            if mirrored and (kind != "R" or tree.skeletons[cid].hubs):
                continue
            skeleton = tree.skeletons[cid].copy()
            if mirrored:
                skeleton.rotation = {
                    v: order[::-1] for v, order in skeleton.rotation.items()
                }
            graph_i = tree.expand(cid, boundaries, skeleton)
            incoming = (
                tree.adjacency[cid][previous][0] if previous is not None else None
            )
            outgoing = (
                tree.adjacency[cid][following][0] if following is not None else None
            )
            # Boundary edges are placeholders, not permissible crossings.
            graph_i.protected.update(e for e in (incoming, outgoing) if e is not None)
            dual = _Dual(graph_i)
            for entering in (0, 1):
                sources = (
                    dual.at_vertex(start)
                    if incoming is None
                    else dual.side(incoming, 1 - entering)
                )
                for leaving in (0, 1):
                    targets = (
                        dual.at_vertex(end)
                        if outgoing is None
                        else dual.side(outgoing, leaving)
                    )
                    local_path = dual.shortest(sources, targets)
                    if local_path is not None and costs[entering] != inf:
                        candidates[leaving].append(
                            (
                                costs[entering] + len(local_path),
                                entering,
                                local_path,
                                graph_i,
                                mirrored,
                            )
                        )
        best = [
            min(items, key=lambda item: item[0]) if items else (inf, 0, [], None, False)
            for items in candidates
        ]
        costs = [entry[0] for entry in best]
        choices.append([entry[1:] for entry in best])
        history.append(
            {
                "component": cid,
                "type": kind,
                "costs": [None if cost == inf else cost for cost in costs],
            }
        )
    if min(costs) == inf:
        raise NoInsertionPath("no ec-insertion path avoids protected edges")
    if costs[0] != costs[1]:
        raise AssertionError("terminal insertion costs must agree")
    selected = []
    side = 0
    for cid, local in reversed(list(zip(path, choices))):
        entering, crossed, graph_i, mirrored = local[side]
        selected.append((cid, graph_i, crossed, mirrored, entering, side))
        side = entering
    selected.reverse()
    for state, (_, _, _, mirrored, entering, leaving) in zip(history, selected):
        state.update(mirrored=mirrored, entering=entering, leaving=leaving)
    embedding = selected[0][1]
    crossings = selected[0][2].copy()
    for index in range(1, len(selected)):
        before, current = selected[index - 1][0], selected[index][0]
        parent_edge, child_edge = tree.adjacency[current][before]
        embedding = splice(selected[index][1], parent_edge, embedding, child_edge)
        crossings.extend(selected[index][2])
    embedding.protected.intersection_update(graph.protected)
    _path_witness(embedding, start, end, crossings)
    return {
        "embedding": embedding,
        "crossings": crossings,
        "cost": int(costs[0]),
        "states": history,
    }


def _path_witness(graph, start, end, crossings):
    """Find faces/corners realizing an already chosen ordered crossing list."""
    dual = _Dual(graph)
    targets = dual.at_vertex(end)
    for source in sorted(dual.at_vertex(start)):
        face, darts = source, []
        for edge in crossings:
            u, v = graph.edges[edge]
            fu, fv = dual.face_of[edge, u], dual.face_of[edge, v]
            if face == fu:
                darts.append((edge, u))
                face = fv
            elif face == fv:
                darts.append((edge, v))
                face = fu
            else:
                break
        else:
            if face in targets:
                start_corner = next(e for e, v in dual.faces[source] if v == start)
                end_corner = next(e for e, v in dual.faces[face] if v == end)
                return {
                    "darts": darts,
                    "start_corner": start_corner,
                    "end_corner": end_corner,
                }
    raise AssertionError("chosen embedding does not realize the insertion path")


def _free_corner(graph, vertex):
    faces, face_of = graph.faces()
    for edge in graph.rotation[vertex]:
        face = faces[face_of[edge, vertex]]
        if not any(v in graph.wheel_hubs for _, v in face):
            return edge
    raise ValueError("cannot attach another block inside a wheel gadget")


def _join_vertex(parent, child, vertex, parent_corner, child_corner):
    if parent.nodes.intersection(child.nodes) != {vertex}:
        raise ValueError("blocks do not meet at exactly their cut vertex")
    out = parent.copy()
    out.nodes.update(child.nodes)
    out.edges.update(child.edges)
    for node, order in child.rotation.items():
        if node != vertex:
            out.rotation[node] = order.copy()
    position = out.rotation[vertex].index(parent_corner) + 1
    child_order = child.rotation[vertex]
    after = child_order.index(child_corner) + 1
    out.rotation[vertex][position:position] = child_order[after:] + child_order[:after]
    out.protected.update(child.protected)
    out.hubs.update(child.hubs)
    out.wheel_hubs.update(child.wheel_hubs)
    return out


def _disjoint_union(parent, child):
    out = parent.copy()
    if out.nodes.intersection(child.nodes):
        raise ValueError("disjoint union has shared vertices")
    out.nodes.update(child.nodes)
    out.edges.update(child.edges)
    out.rotation.update({v: order.copy() for v, order in child.rotation.items()})
    out.protected.update(child.protected)
    out.hubs.update(child.hubs)
    out.wheel_hubs.update(child.wheel_hubs)
    return out


def glue_blocks(graph, embedded_blocks):
    """Attach feasible expanded blocks only through non-wheel face corners."""
    from ec_planarity import Graph

    pending = list(embedded_blocks)
    out = Graph()
    while pending:
        for index, current in enumerate(pending):
            shared = out.nodes.intersection(current.nodes)
            if shared:
                if len(shared) != 1:
                    raise ValueError("blocks share more than one cut vertex")
                vertex = next(iter(shared))
                out = _join_vertex(
                    out,
                    current,
                    vertex,
                    _free_corner(out, vertex),
                    _free_corner(current, vertex),
                )
                pending.pop(index)
                break
        else:
            out = _disjoint_union(out, pending.pop(0))
    out.nodes.update(graph.nodes)
    for vertex in graph.nodes:
        out.rotation.setdefault(vertex, [])
    return out


def optimal_insertion(graph: Graph, start: str, end: str):
    """Minimize crossings over all ec-embeddings, preserving protected edges.

    Implements the block-vertex-tree wrapper in Algorithm 4. Endpoints must be
    distinct and belong to the same connected component, as in the paper.
    """
    if start == end or start not in graph.nodes or end not in graph.nodes:
        raise ValueError("insertion requires two distinct graph vertices")
    edge_blocks = blocks(graph)
    incidence = {v: [] for v in graph.nodes}
    subgraphs = []
    for index, edges in enumerate(edge_blocks):
        block = graph.subgraph(edges)
        subgraphs.append(block)
        for vertex in block.nodes:
            incidence[vertex].append(index)
    source, target = ("v", start), ("v", end)
    queue, predecessors = deque([source]), {source: None}
    while queue and target not in predecessors:
        kind, current = queue.popleft()
        adjacent = (
            [("b", b) for b in incidence[current]]
            if kind == "v"
            else [("v", v) for v in sorted(subgraphs[current].nodes)]
        )
        for neighbor in adjacent:
            if neighbor not in predecessors:
                predecessors[neighbor] = (kind, current)
                queue.append(neighbor)
    if target not in predecessors:
        raise ValueError("insertion endpoints must lie in one connected component")
    sequence, cursor = [], target
    while cursor is not None:
        sequence.append(cursor)
        cursor = predecessors[cursor]
    sequence.reverse()
    embedding, crossings, used, history = None, [], set(), []
    for offset in range(1, len(sequence), 2):
        index = sequence[offset][1]
        x, y = sequence[offset - 1][1], sequence[offset + 1][1]
        result = _block_insertion(subgraphs[index], x, y)
        current, path = result["embedding"], result["crossings"]
        if embedding is None:
            embedding = current
        else:
            left = _path_witness(embedding, start, x, crossings)
            right = _path_witness(current, x, y, path)
            embedding = _join_vertex(
                embedding, current, x, left["end_corner"], right["start_corner"]
            )
        crossings.extend(path)
        used.add(index)
        history.extend(result["states"])
        _path_witness(embedding, start, y, crossings)
    # Off-path blocks can be attached to any non-wheel face at their cut vertex.
    pending = set(range(len(subgraphs))) - used
    while pending:
        progress = False
        for index in sorted(pending):
            current = subgraphs[index]
            shared = embedding.nodes.intersection(current.nodes)
            if not shared:
                continue
            if len(shared) != 1:
                raise AssertionError("invalid block decomposition")
            vertex = next(iter(shared))
            embedding = _join_vertex(
                embedding,
                current,
                vertex,
                _free_corner(embedding, vertex),
                _free_corner(current, vertex),
            )
            pending.remove(index)
            progress = True
        if not progress:
            # Start an unrelated component, then attach its remaining blocks.
            index = min(pending)
            current = subgraphs[index]
            embedding = _disjoint_union(embedding, current)
            pending.remove(index)
    embedding.nodes.update(graph.nodes)
    for vertex in graph.nodes:
        embedding.rotation.setdefault(vertex, [])
    _path_witness(embedding, start, end, crossings)
    return {
        "embedding": embedding,
        "crossings": crossings,
        "cost": len(crossings),
        "start": start,
        "end": end,
        "states": history,
    }


def planarize_insertion(result, edge_id: str):
    """Realize the selected insertion with alternating degree-four crossings."""
    graph = result["embedding"]
    start, end, crossed = result["start"], result["end"], result["crossings"]
    if edge_id in graph.edges:
        raise ValueError("inserted edge ID already exists")
    if len(crossed) != len(set(crossed)):
        raise ValueError("an optimal path must not cross one edge twice")
    witness = _path_witness(graph, start, end, crossed)
    out = graph.copy()
    chains = {edge: [edge] for edge in graph.edges}
    crossing_nodes = {}

    def fresh(base, existing):
        candidate, serial = base, 0
        while candidate in existing:
            serial += 1
            candidate = f"{base}:{serial}"
        return candidate

    nodes = [start]
    split_edges = []
    for index, (edge, entering_origin) in enumerate(witness["darts"]):
        if edge in graph.protected:
            raise AssertionError("a protected edge appears in the insertion path")
        node = fresh(f"@cross:{edge_id}:{index}", out.nodes)
        out.nodes.add(node)
        nodes.append(node)
        u, v = out.edges.pop(edge)
        first = fresh(f"{edge}:part:0", out.edges)
        out.edges[first] = (u, node)
        second = fresh(f"{edge}:part:1", out.edges)
        out.edges[second] = (node, v)
        out.rotation[u][out.rotation[u].index(edge)] = first
        out.rotation[v][out.rotation[v].index(edge)] = second
        chains[edge] = [first, second]
        back, forward = (first, second) if entering_origin == u else (second, first)
        split_edges.append((back, forward))
        crossing_nodes[node] = {"crossed_edge": edge, "inserted_edge": edge_id}
    nodes.append(end)
    insertion_edges = []
    for index, (u, v) in enumerate(pairwise(nodes)):
        segment = (
            edge_id if not crossed else fresh(f"{edge_id}:part:{index}", out.edges)
        )
        out.edges[segment] = (u, v)
        insertion_edges.append(segment)
    for index, node in enumerate(nodes[1:-1]):
        back, forward = split_edges[index]
        out.rotation[node] = [
            back,
            insertion_edges[index + 1],
            forward,
            insertion_edges[index],
        ]
    for vertex, corner, segment in (
        (start, witness["start_corner"], insertion_edges[0]),
        (end, witness["end_corner"], insertion_edges[-1]),
    ):
        # A corner edge may have been subdivided; choose its incident segment.
        if corner not in out.edges:
            corner = next(e for e in chains[corner] if vertex in out.edges[e])
        position = out.rotation[vertex].index(corner) + 1
        out.rotation[vertex].insert(position, segment)
    chains[edge_id] = insertion_edges
    return {
        "embedding": out,
        "edge_chains": chains,
        "crossing_nodes": crossing_nodes,
        "crossings": crossed.copy(),
        "cost": len(crossed),
    }
