"""Independent checks of full-tree port reservation and insertion witnesses."""

import itertools
import json
import subprocess
import sys
import tempfile
import unittest
from collections import deque
from pathlib import Path

from ec_planarity import ECExpansion, Graph, admissible, cyclic_equal
from insertion import NoInsertionPath, optimal_insertion, planarize_insertion
from planarization import planarize


class PaperFidelity(unittest.TestCase):
    @staticmethod
    def tree_orders(tree):
        """Tiny factorial oracle, independent of the production tree checker."""
        if isinstance(tree, str):
            return [(tree,)]
        children = tree["children"]
        indices = tuple(range(len(children)))
        permutations = {
            "oriented": [indices],
            "mirror": [indices, indices[::-1]],
            "group": list(itertools.permutations(indices)),
        }[tree["kind"]]
        orders = []
        for permutation in permutations:
            for choices in itertools.product(
                *(PaperFidelity.tree_orders(children[i]) for i in permutation)
            ):
                orders.append(sum(choices, ()))
        return orders

    @staticmethod
    def fixed_rotation_cost(graph, full_orders, start, end, inserted):
        """Direct dual BFS restricted to admissible endpoint corners."""
        faces = {}
        for edge, ends in graph.edges.items():
            for vertex in ends:
                initial = (edge, vertex)
                if initial in faces:
                    continue
                face = max(faces.values(), default=-1) + 1
                dart = initial
                while dart not in faces:
                    faces[dart] = face
                    current, origin = dart
                    a, b = graph.edges[current]
                    target = b if origin == a else a
                    rotation = graph.rotation[target]
                    dart = (rotation[rotation.index(current) - 1], target)
                if dart != initial:
                    raise AssertionError("invalid face walk")
        adjacency = {face: set() for face in faces.values()}
        for edge, (a, b) in graph.edges.items():
            if edge not in graph.protected:
                first, second = faces[edge, a], faces[edge, b]
                adjacency[first].add(second)
                adjacency[second].add(first)

        def allowed(vertex):
            expected = full_orders[vertex]
            result = set()
            for index, corner in enumerate(graph.rotation[vertex]):
                order = list(graph.rotation[vertex])
                order.insert(index + 1, inserted)
                if cyclic_equal(order, expected):
                    result.add(faces[corner, vertex])
            return result

        targets = allowed(end)
        queue = deque((face, 0) for face in allowed(start))
        seen = set()
        while queue:
            face, cost = queue.popleft()
            if face in targets:
                return cost
            if face not in seen:
                seen.add(face)
                queue.extend((neighbor, cost + 1) for neighbor in adjacency[face])
        return None

    @staticmethod
    def full_k4_with_reserved_edge(first_slot, second_slot):
        pairs = list(itertools.combinations("abcd", 2))
        host = ECExpansion(
            Graph(set("abcd"), {str(i): pair for i, pair in enumerate(pairs)}), {}
        ).embed()
        graph = host.copy()
        graph.edges["pending"] = ("a", "b")
        full_orders = {node: list(order) for node, order in host.rotation.items()}
        full_orders["a"].insert(first_slot + 1, "pending")
        full_orders["b"].insert(second_slot + 1, "pending")
        graph.rotation = full_orders
        expansion = ECExpansion(
            graph,
            {
                node: {"kind": "oriented", "children": order}
                for node, order in full_orders.items()
            },
        )
        start, end = expansion.graph.edges["pending"]
        missing = expansion.graph.copy()
        del missing.edges["pending"]
        missing.rotation[start].remove("pending")
        missing.rotation[end].remove("pending")
        return host, graph, expansion, missing, start, end, full_orders

    def test_reserved_ports_match_all_nine_fixed_endpoint_corners(self):
        costs = set()
        for first, second in itertools.product(range(3), repeat=2):
            with self.subTest(first=first, second=second):
                host, full, expansion, missing, start, end, orders = (
                    self.full_k4_with_reserved_edge(first, second)
                )
                expected = self.fixed_rotation_cost(host, orders, "a", "b", "pending")
                # All wheel/rim/tree edges and both unused attachment slots survive.
                self.assertEqual(missing.protected, expansion.graph.protected)
                self.assertEqual(missing.nodes, expansion.graph.nodes)
                result = optimal_insertion(missing, start, end)
                self.assertEqual(result["cost"], expected)
                self.assertFalse(set(result["crossings"]) & missing.protected)
                self.assertLessEqual(set(result["crossings"]), set(host.edges))
                costs.add(expected)
                planarized = planarize_insertion(result, "pending")
                drawing = planarized["embedding"]
                drawing.validate_embedding()
                self.assertEqual(drawing.hubs, expansion.graph.hubs)
                self.assertEqual(drawing.protected, expansion.graph.protected)
                chains = planarized["edge_chains"]
                owners = {
                    segment: original
                    for original, chain in chains.items()
                    for segment in chain
                }
                for edge, ends in expansion.graph.edges.items():
                    cursor = ends[0]
                    for segment in chains[edge]:
                        u, v = drawing.edges[segment]
                        self.assertEqual(u, cursor)
                        cursor = v
                    self.assertEqual(cursor, ends[1])
                for crossing in planarized["crossing_nodes"]:
                    rotation = drawing.rotation[crossing]
                    self.assertEqual(len(rotation), 4)
                    physical = [owners[edge] for edge in rotation]
                    self.assertEqual(physical[0], physical[2])
                    self.assertEqual(physical[1], physical[3])
                    self.assertNotEqual(physical[0], physical[1])
                if expected == 0:
                    # With no subdivision, contraction recovers every original port.
                    collapsed = expansion.collapse(drawing)
                    self.assertEqual(collapsed.edges, full.edges)
                    for vertex, order in orders.items():
                        self.assertTrue(cyclic_equal(collapsed.rotation[vertex], order))
        self.assertEqual(costs, {0, 1})

    def test_protected_cycle_can_make_reserved_insertion_impossible(self):
        # Distinct fixed K4 faces require a crossing, regardless of R-node mirroring.
        for first, second in itertools.product(range(3), repeat=2):
            host, _, _, missing, start, end, orders = self.full_k4_with_reserved_edge(
                first, second
            )
            if self.fixed_rotation_cost(host, orders, "a", "b", "pending") == 1:
                missing.protected.update(host.edges)
                with self.assertRaises(NoInsertionPath):
                    optimal_insertion(missing, start, end)
                return
        self.fail("fixture has no separated endpoint corners")

    def test_mirror_parent_does_not_reverse_oriented_child(self):
        tree = {
            "kind": "mirror",
            "children": [
                {"kind": "oriented", "children": ["a", "b", "c"]},
                {"kind": "group", "children": ["d", "e"]},
                "f",
            ],
        }
        allowed = self.tree_orders(tree)
        for rotation in itertools.permutations("abcdef"):
            expected = any(cyclic_equal(rotation, order) for order in allowed)
            self.assertEqual(admissible(tree, list(rotation)), expected, rotation)

    def test_ec_nonplanar_k4_recovers_original_oriented_orders(self):
        physical = ECExpansion(
            Graph(
                set("abcd"),
                {
                    str(i): pair
                    for i, pair in enumerate(itertools.combinations("abcd", 2))
                },
            ),
            {},
        ).embed()
        orders = {"a": physical.rotation["a"], "b": physical.rotation["b"][::-1]}
        expansion = ECExpansion(
            physical,
            {
                node: {"kind": "oriented", "children": order}
                for node, order in orders.items()
            },
        )
        result = planarize(expansion.graph)
        self.assertGreater(len(result["crossing_nodes"]), 0)
        self.assertGreater(len(result["insertions"]), 0)
        result["embedding"].validate_embedding()
        self.assertEqual(result["embedding"].hubs, expansion.graph.hubs)
        collapsed = expansion.collapse(result["embedding"])
        owners = {
            segment: owner
            for owner, chain in result["edge_chains"].items()
            for segment in chain
        }
        for node, expected in orders.items():
            actual = [owners[segment] for segment in collapsed.rotation[node]]
            self.assertTrue(cyclic_equal(actual, expected), (actual, expected))
        self.assertTrue(set(physical.edges) <= result["edge_chains"].keys())
        self.assertLessEqual(
            set(result["embedding"].protected), expansion.graph.protected
        )

    def test_prior_crossing_mirrors_preserve_alternation_and_original_o_hub(self):
        physical = Graph(
            set("abcdef"),
            {
                str(i): pair
                for i, pair in enumerate(itertools.combinations("abcdef", 2))
            },
        )
        expected = physical.rotation["a"].copy()
        expansion = ECExpansion(
            physical, {"a": {"kind": "oriented", "children": expected}}
        )
        result = planarize(expansion.graph)
        self.assertGreaterEqual(len(result["insertions"]), 2)
        self.assertGreaterEqual(len(result["crossing_nodes"]), 3)
        drawing = expansion.collapse(result["embedding"])
        owners = {
            segment: owner
            for owner, chain in result["edge_chains"].items()
            for segment in chain
        }
        self.assertTrue(
            cyclic_equal([owners[e] for e in drawing.rotation["a"]], expected)
        )
        for node, record in result["crossing_nodes"].items():
            order = [owners[edge] for edge in drawing.rotation[node]]
            self.assertEqual(len(order), 4)
            self.assertEqual(order[0], order[2])
            self.assertEqual(order[1], order[3])
            self.assertEqual(set(order), set(record.values()))
        for edge, endpoints in physical.edges.items():
            cursor = endpoints[0]
            for segment in result["edge_chains"][edge]:
                start, end = drawing.edges[segment]
                self.assertEqual(start, cursor)
                cursor = end
            self.assertEqual(cursor, endpoints[1])
        drawing.validate_embedding()

    def test_unreachable_dp_side_and_cli_output_are_strict_json(self):
        edge_ids = ["ac", "ad", "bc", "bd", "cd", "ae", "af", "be", "bf", "ef", "ab"]
        graph = Graph(set("abcdef"), {edge: tuple(edge) for edge in edge_ids})
        graph.protected.update(graph.edges)
        expansion = ECExpansion(
            graph,
            {
                "d": {"kind": "oriented", "children": ["ad", "bd", "cd"]},
                "f": {"kind": "oriented", "children": ["ef", "bf", "af"]},
            },
        )
        result = optimal_insertion(expansion.embed_expanded(), "c", "e")
        self.assertEqual(result["cost"], 0)
        self.assertTrue(any(None in state["costs"] for state in result["states"]))
        json.dumps(result["states"], allow_nan=False)

        for first, second in itertools.product(range(3), repeat=2):
            host, graph, expansion, _, _, _, orders = self.full_k4_with_reserved_edge(
                first, second
            )
            if self.fixed_rotation_cost(host, orders, "a", "b", "pending") != 1:
                continue
            with tempfile.TemporaryDirectory() as directory:
                problem = Path(directory) / "problem.json"
                problem.write_text(
                    json.dumps(
                        {
                            "nodes": sorted(graph.nodes),
                            "edges": graph.edges,
                            "constraints": expansion.constraints,
                        }
                    )
                )
                command = [
                    sys.executable,
                    str(Path(__file__).with_name("cli.py").resolve()),
                    str(problem),
                    "--mode",
                    "insert",
                    "--insert",
                    "pending",
                ]
                process = subprocess.run(
                    command,
                    cwd=directory,
                    capture_output=True,
                    text=True,
                    check=True,
                    timeout=30,
                )
                self.assertEqual(process.stderr, "")
                response = json.loads(process.stdout)
                self.assertEqual(response["cost"], 1)
                json.dumps(response, allow_nan=False)
            break
        else:
            self.fail("fixture has no positive-cost insertion")


if __name__ == "__main__":
    unittest.main()
