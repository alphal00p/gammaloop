"""Independent small-rotation oracles for constrained variable edge insertion."""

import math
import random
import unittest
from collections import deque
from itertools import combinations, permutations, product

from ec_planarity import ECExpansion, Graph, decompose, embed_expanded_graph
from insertion import (
    NoInsertionPath,
    blocks,
    optimal_insertion,
    planarize_insertion,
    skeleton_edge_id,
)


def brute_insertion(
    graph, start, end, fixed=None, endpoint_orders=None, new_edge="new"
):
    """Enumerate original rotations, independently of SPQR and ec-expansion."""
    fixed = fixed or {}
    endpoint_orders = endpoint_orders or {}
    vertices = sorted(graph.nodes)
    choices = []
    for vertex in vertices:
        incident = [e for e, ends in graph.edges.items() if vertex in ends]
        if vertex in fixed:
            choices.append([fixed[vertex]])
        else:
            choices.append([[incident[0], *p] for p in permutations(incident[1:])])
    best = float("inf")
    for candidate in product(*choices):
        rotation = dict(zip(vertices, candidate))
        darts = {
            (edge, origin) for edge, ends in graph.edges.items() for origin in ends
        }
        face_of, faces = {}, []
        while darts:
            dart = next(iter(darts))
            face = []
            while dart in darts:
                darts.remove(dart)
                face_of[dart] = len(faces)
                face.append(dart)
                edge, origin = dart
                u, v = graph.edges[edge]
                target = v if origin == u else u
                order = rotation[target]
                dart = (order[(order.index(edge) - 1) % len(order)], target)
            faces.append(face)
        if len(vertices) - len(graph.edges) + len(faces) != 2:
            continue
        adjacent = [[] for _ in faces]
        for edge, (u, v) in graph.edges.items():
            if edge not in graph.protected:
                a, b = face_of[edge, u], face_of[edge, v]
                adjacent[a].append(b)
                adjacent[b].append(a)

        def allowed_faces(vertex, rotation=rotation, face_of=face_of):
            result = set()
            for index, edge in enumerate(rotation[vertex]):
                order = rotation[vertex]
                inserted = order[: index + 1] + [new_edge] + order[index + 1 :]
                expected = endpoint_orders.get(vertex)
                if expected is not None:
                    pivot = inserted.index(expected[0])
                    if inserted[pivot:] + inserted[:pivot] != expected:
                        continue
                result.add(face_of[edge, vertex])
            return result

        sources = allowed_faces(start)
        targets = allowed_faces(end)
        distances = {f: 0 for f in sources}
        queue = deque(sources)
        while queue:
            face = queue.popleft()
            if face in targets:
                best = min(best, distances[face])
                break
            for neighbor in adjacent[face]:
                if neighbor not in distances:
                    distances[neighbor] = distances[face] + 1
                    queue.append(neighbor)
    return best


def two_tetrahedra(parallel=False):
    edges = ["ac", "ad", "bc", "bd", "cd", "ae", "af", "be", "bf", "ef"]
    if parallel:
        edges += ["ab"]
    return Graph(set("abcdef"), {e: tuple(e) for e in edges})


class InsertionTests(unittest.TestCase):
    def check_result(self, result, edge="new"):
        result["embedding"].validate_embedding()
        planarization = planarize_insertion(result, edge)
        planarization["embedding"].validate_embedding()
        graph = planarization["embedding"]
        self.assertEqual(result["cost"], len(planarization["crossing_nodes"]))
        for node, crossing in planarization["crossing_nodes"].items():
            order = graph.rotation[node]
            original = set(planarization["edge_chains"][crossing["crossed_edge"]])
            inserted = set(planarization["edge_chains"][edge])
            self.assertEqual(len(order), 4)
            self.assertEqual([e in original for e in order], [True, False, True, False])
            self.assertTrue(order[1] in inserted and order[3] in inserted)
        for old, chain in planarization["edge_chains"].items():
            if old == edge:
                source, target = result["start"], result["end"]
            else:
                source, target = result["embedding"].edges[old]
            point = source
            for segment in chain:
                self.assertEqual(point, graph.edges[segment][0])
                point = graph.edges[segment][1]
            self.assertEqual(point, target)
        return planarization

    def test_rigid_one_crossing_and_protected_edges(self):
        original = Graph(
            {str(i) for i in range(5)},
            {
                f"{i}-{j}": (str(i), str(j))
                for i, j in combinations(range(5), 2)
                if (i, j) != (0, 1)
            },
        )
        embedded = embed_expanded_graph(original)
        result = optimal_insertion(embedded, "0", "1")
        self.assertEqual(result["cost"], brute_insertion(original, "0", "1"))
        self.assertEqual(result["cost"], 1)
        self.check_result(result)
        original.protected = set(original.edges) - {"2-3"}
        result = optimal_insertion(embed_expanded_graph(original), "0", "1")
        self.assertEqual(result["crossings"], ["2-3"])
        self.check_result(result)
        original.protected.add("2-3")
        with self.assertRaises(NoInsertionPath):
            optimal_insertion(embed_expanded_graph(original), "0", "1")

    def test_oriented_rr_and_rpr_match_rotation_oracle(self):
        unequal_sides = False
        saw_p = False
        costs = set()
        for parallel, reverse_d, reverse_f in product((False, True), repeat=3):
            original = two_tetrahedra(parallel)
            fixed = {"d": ["ad", "bd", "cd"], "f": ["af", "bf", "ef"]}
            if reverse_d:
                fixed["d"].reverse()
            if reverse_f:
                fixed["f"].reverse()
            constraints = {
                v: {"kind": "oriented", "children": order} for v, order in fixed.items()
            }
            expansion = ECExpansion(original, constraints)
            result = optimal_insertion(expansion.embed_expanded(), "c", "e")
            self.assertEqual(result["cost"], brute_insertion(original, "c", "e", fixed))
            costs.add(result["cost"])
            unequal_sides |= any(
                s["costs"][0] != s["costs"][1] for s in result["states"]
            )
            saw_p |= any(s["type"] == "P" for s in result["states"])
            self.check_result(result)
        self.assertTrue(
            unequal_sides, "the test must exercise the two distinct DP states"
        )
        self.assertTrue(
            saw_p, "the test must exercise the P-node permutation transition"
        )
        self.assertEqual(costs, {0, 1})

    def test_full_expansion_preserves_endpoint_insertion_slots(self):
        full = Graph(
            {str(i) for i in range(5)},
            {f"{i}-{j}": (str(i), str(j)) for i, j in combinations(range(5), 2)},
        )
        original = embed_expanded_graph(full.subgraph(set(full.edges) - {"0-1"}))
        costs = set()
        for left, right in product(range(3), repeat=2):
            orders = {}
            for vertex, slot in [("0", left), ("1", right)]:
                incident = original.rotation[vertex]
                orders[vertex] = incident[: slot + 1] + ["0-1"] + incident[slot + 1 :]
            constraints = {
                v: {"kind": "oriented", "children": order}
                for v, order in orders.items()
            }
            expansion = ECExpansion(full, constraints)
            start, end = expansion.graph.edges["0-1"]
            base = expansion.graph.subgraph(set(expansion.graph.edges) - {"0-1"})
            result = optimal_insertion(embed_expanded_graph(base), start, end)
            expected = brute_insertion(
                original, "0", "1", endpoint_orders=orders, new_edge="0-1"
            )
            self.assertEqual(result["cost"], expected)
            costs.add(expected)
            self.check_result(result, "0-1")
        self.assertGreater(
            len(costs), 1, "endpoint slots must change the insertion optimum"
        )

    def test_three_rigid_components_need_mirrored_state_transitions(self):
        edges = [
            "ac",
            "ad",
            "bc",
            "bd",
            "cd",
            "af",
            "be",
            "bf",
            "ef",
            "ag",
            "ah",
            "eg",
            "eh",
            "gh",
        ]
        original = Graph(set("abcdefgh"), {edge: tuple(edge) for edge in edges})
        saw_mirror = False
        for reverse_d, reverse_h in product((False, True), repeat=2):
            fixed = {"d": ["ad", "bd", "cd"], "h": ["ah", "eh", "gh"]}
            if reverse_d:
                fixed["d"].reverse()
            if reverse_h:
                fixed["h"].reverse()
            constraints = {
                v: {"kind": "oriented", "children": order} for v, order in fixed.items()
            }
            expanded = ECExpansion(original, constraints).embed_expanded()
            result = optimal_insertion(expanded, "c", "g")
            self.assertEqual(result["cost"], brute_insertion(original, "c", "g", fixed))
            self.assertEqual(len(result["states"]), 3)
            saw_mirror |= any(s["mirrored"] for s in result["states"])
            self.check_result(result)
        self.assertTrue(saw_mirror)

    def test_virtual_names_cannot_shadow_original_edge_ids(self):
        graph = two_tetrahedra(parallel=True)
        spqr = decompose(graph)
        virtual = next(
            (c, e)
            for c in spqr["components"]
            for e in c["edges"]
            if e["twin"] is not None
        )
        collision = skeleton_edge_id(virtual[0]["id"], virtual[1]["id"])
        graph = Graph(
            graph.nodes,
            {collision if e == "ab" else e: ends for e, ends in graph.edges.items()},
        )
        spqr = decompose(graph)
        self.assertTrue(
            any(
                skeleton_edge_id(c["id"], e["id"]) in graph.edges
                for c in spqr["components"]
                for e in c["edges"]
                if e["twin"] is not None
            )
        )
        result = optimal_insertion(embed_expanded_graph(graph), "c", "e")
        self.assertEqual(set(result["embedding"].edges), set(graph.edges))
        self.check_result(result)

    def test_reversible_parent_keeps_oriented_off_path_child(self):
        original = two_tetrahedra(parallel=True)
        constraint = {"f": {"kind": "oriented", "children": ["af", "bf", "ef"]}}
        expanded = ECExpansion(original, constraint).embed_expanded()
        for source, target in [("c", "e"), ("d", "e"), ("c", "f")]:
            # f is expanded, so only the first two pairs are original endpoints.
            if target not in expanded.nodes:
                continue
            result = optimal_insertion(expanded, source, target)
            self.assertEqual(
                result["cost"],
                brute_insertion(original, source, target, {"f": ["af", "bf", "ef"]}),
            )
            self.check_result(result)

    def test_block_path_glues_faces_and_retains_off_path_wheel(self):
        edges = {}
        for prefix, a, b in [("l", "start", "cut"), ("r", "cut", "end")]:
            nodes = [a, b, *[prefix + str(i) for i in range(3)]]
            for i, j in combinations(range(5), 2):
                if (i, j) != (0, 1):
                    edges[f"{prefix}{i}{j}"] = (nodes[i], nodes[j])
        edges.update({"ta": ("cut", "t"), "tb": ("t", "u"), "tc": ("u", "cut")})
        original = Graph({v for pair in edges.values() for v in pair}, edges)
        tree = {"t": {"kind": "oriented", "children": ["ta", "tb"]}}
        expanded = ECExpansion(original, tree).embed_expanded()
        result = optimal_insertion(expanded, "start", "end")
        self.assertEqual(result["cost"], 2)
        self.assertEqual(set(result["embedding"].edges), set(expanded.edges))
        self.check_result(result)

    def test_parallel_block_and_bridge(self):
        graph = Graph(
            set("abcd"),
            {"ab1": ("a", "b"), "ab2": ("a", "b"), "bc": ("b", "c"), "cd": ("c", "d")},
        )
        self.assertEqual(sorted(map(len, blocks(graph))), [1, 1, 2])
        result = optimal_insertion(embed_expanded_graph(graph), "a", "d")
        self.assertEqual(result["cost"], 0)
        self.check_result(result)

    def test_seeded_graphs_and_orientations_match_exact_oracle(self):
        rng = random.Random(1602008)
        tested = 0
        while tested < 40:
            n = rng.choice([5, 6])
            pairs = [(i, rng.randrange(i)) for i in range(1, n)]
            unused = [
                p
                for p in combinations(range(n), 2)
                if p not in pairs and p[::-1] not in pairs
            ]
            rng.shuffle(unused)
            pairs += unused[: rng.randint(1, 4)]
            remaining = [
                p
                for p in combinations(range(n), 2)
                if p not in pairs and p[::-1] not in pairs
            ]
            if not remaining:
                continue
            a, b = rng.choice(remaining)
            graph = Graph(
                {str(i) for i in range(n)},
                {str(i): (str(u), str(v)) for i, (u, v) in enumerate(pairs)},
            )
            if (
                math.prod(
                    math.factorial(len(order) - 1) for order in graph.rotation.values()
                )
                > 20000
            ):
                continue
            graph = embed_expanded_graph(graph)
            fixed = {
                v: order
                for v, order in graph.rotation.items()
                if v not in (str(a), str(b)) and len(order) > 2 and rng.random() < 0.5
            }
            graph.protected = {e for e in graph.edges if rng.random() < 0.12}
            expansion = ECExpansion(
                graph,
                {
                    v: {"kind": "oriented", "children": order}
                    for v, order in fixed.items()
                },
            )
            embedded = expansion.embed_expanded()
            before = embedded.copy()
            expected = brute_insertion(graph, str(a), str(b), fixed)
            if expected == float("inf"):
                with self.assertRaises(NoInsertionPath):
                    optimal_insertion(embedded, str(a), str(b))
            else:
                result = optimal_insertion(embedded, str(a), str(b))
                self.assertEqual(result["cost"], expected)
                self.check_result(result)
            self.assertEqual(
                embedded, before, "insertion must leave its input untouched"
            )
            tested += 1

    def test_small_planar_inputs_match_exhaustive_oracle(self):
        universe = list(combinations(range(5), 2))
        for absent in [(0,), (1,), (2,), (0, 1), (0, 2), (1, 5), (0, 4, 9), (1, 3, 7)]:
            edges = {
                str(i): (str(a), str(b))
                for i, (a, b) in enumerate(universe)
                if i not in absent
            }
            original = Graph({str(i) for i in range(5)}, edges)
            embedded = embed_expanded_graph(original)
            a, b = universe[absent[0]]
            expected = brute_insertion(original, str(a), str(b))
            result = optimal_insertion(embedded, str(a), str(b))
            self.assertEqual(result["cost"], expected)
            self.check_result(result)


if __name__ == "__main__":
    unittest.main()
