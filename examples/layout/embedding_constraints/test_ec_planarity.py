"""Small exhaustive rotation oracles, independent of expansion and SPQR."""

import itertools
import random
import unittest

from ec_planarity import ECExpansion, Graph, NotECPlanar, admissible, cyclic_equal


def graph_from_pairs(pairs):
    return Graph(
        {node for ends in pairs for node in ends},
        {f"e{i}": ends for i, ends in enumerate(pairs)},
    )


def all_planar_rotations(graph, constraints=None):
    """Factorial enumeration is only a tiny-graph test oracle."""
    constraints = constraints or {}
    nodes = sorted(graph.nodes)
    choices = []
    for node in nodes:
        edges = graph.rotation[node]
        orders = (
            ([edges[0], *tail] for tail in itertools.permutations(edges[1:]))
            if edges
            else [[]]
        )
        choices.append(
            [
                order
                for order in orders
                if node not in constraints or admissible(constraints[node], order)
            ]
        )
    for orders in itertools.product(*choices):
        candidate = graph.copy()
        candidate.rotation = dict(zip(nodes, orders))
        try:
            candidate.validate_embedding()
        except NotECPlanar:
            continue
        yield candidate


def tree(kind, *children):
    return {"kind": kind, "children": list(children)}


class ConstraintSemantics(unittest.TestCase):
    def test_nested_group_mirror_and_oriented(self):
        constraint = tree(
            "oriented", tree("group", "a", "b"), tree("mirror", "c", "d", "e"), "f"
        )
        self.assertTrue(admissible(constraint, ["b", "a", "e", "d", "c", "f"]))
        self.assertTrue(admissible(constraint, ["f", "a", "b", "c", "d", "e"]))
        self.assertFalse(admissible(constraint, ["a", "b", "c", "e", "d", "f"]))
        self.assertFalse(admissible(constraint, ["a", "b", "f", "c", "d", "e"]))
        self.assertFalse(admissible(constraint, ["a", "c", "b", "d", "e", "f"]))

    def test_expansion_has_two_rim_vertices_per_incidence(self):
        graph = graph_from_pairs([("v", str(i)) for i in range(4)])
        expansion = ECExpansion(graph, {"v": tree("oriented", *graph.rotation["v"])})
        self.assertEqual(len(expansion.graph.wheel_hubs), 1)
        hub = next(iter(expansion.graph.wheel_hubs))
        self.assertEqual(len(expansion.graph.hubs[hub]), 8)
        self.assertEqual(len(expansion.graph.protected), 16)
        self.assertEqual(len(expansion.members["v"]), 9)
        self.assertEqual(set(expansion.leaf_ports["v"]), set(graph.edges))

    def test_constraint_validation(self):
        graph = graph_from_pairs([("v", "a"), ("v", "b"), ("v", "c")])
        for constraint in [
            tree("group", "e0", "e0", "e2"),
            tree("group", "e0", "e1"),
            tree("group", "e0"),
            tree("unknown", "e0", "e1", "e2"),
        ]:
            with self.assertRaises(ValueError):
                ECExpansion(graph, {"v": constraint})

    def test_parallel_darts_are_distinct(self):
        graph = graph_from_pairs([("a", "b")] * 3)
        graph.rotation["b"].reverse()
        self.assertEqual(len(graph.faces()[0]), 3)
        graph.validate_embedding()
        graph.rotation["b"].reverse()
        with self.assertRaises(NotECPlanar):
            graph.validate_embedding()


class EcPlanarity(unittest.TestCase):
    def test_oriented_conflict_inside_same_r_skeleton(self):
        graph = graph_from_pairs(list(itertools.combinations("abcd", 2)))
        oracle = next(all_planar_rotations(graph))
        constraints = {node: tree("oriented", *oracle.rotation[node]) for node in "ab"}
        result = ECExpansion(graph, constraints).embed()
        for node in constraints:
            self.assertTrue(cyclic_equal(result.rotation[node], oracle.rotation[node]))
        constraints["b"]["children"].reverse()
        self.assertEqual(list(all_planar_rotations(graph, constraints)), [])
        with self.assertRaises(NotECPlanar):
            ECExpansion(graph, constraints).embed()

    def test_mirror_can_reverse_where_oriented_cannot(self):
        graph = graph_from_pairs(list(itertools.combinations("abcd", 2)))
        oracle = next(all_planar_rotations(graph))
        constraints = {
            "a": tree("oriented", *oracle.rotation["a"]),
            "b": tree("mirror", *oracle.rotation["b"][::-1]),
        }
        result = ECExpansion(graph, constraints).embed()
        self.assertTrue(admissible(constraints["b"], result.rotation["b"]))

    def test_constraints_merge_blocks_before_planarity_test(self):
        graph = graph_from_pairs(
            [("v", "a"), ("a", "b"), ("b", "v"), ("v", "c"), ("c", "d"), ("d", "v")]
        )
        compatible = {"v": tree("oriented", "e0", "e2", "e3", "e5")}
        ECExpansion(graph, compatible).embed()
        alternating = {"v": tree("oriented", "e0", "e3", "e2", "e5")}
        self.assertFalse(any(all_planar_rotations(graph, alternating)))
        with self.assertRaises(NotECPlanar):
            ECExpansion(graph, alternating).embed()

    def test_independent_r_skeletons_can_have_different_orientations(self):
        graph = graph_from_pairs(
            list(itertools.combinations("abcd", 2))
            + [ends for ends in itertools.combinations("abef", 2) if ends != ("a", "b")]
        )
        baseline = ECExpansion(graph, {}).embed()
        constraints = {
            "c": tree("oriented", *baseline.rotation["c"]),
            "e": tree("oriented", *baseline.rotation["e"][::-1]),
        }
        result = ECExpansion(graph, constraints).embed()
        self.assertTrue(admissible(constraints["c"], result.rotation["c"]))
        self.assertTrue(admissible(constraints["e"], result.rotation["e"]))

    def test_nested_groups_and_parallel_edges(self):
        graph = graph_from_pairs([("a", "b")] * 5)
        constraints = {
            "a": tree(
                "oriented", tree("group", "e0", "e1"), tree("mirror", "e2", "e3"), "e4"
            ),
            "b": tree("oriented", "e4", "e3", "e2", "e1", "e0"),
        }
        result = ECExpansion(graph, constraints).embed()
        for node, constraint in constraints.items():
            self.assertTrue(admissible(constraint, result.rotation[node]))

    def test_empty_isolated_and_bridge_components(self):
        graph = Graph({"isolated", "a", "b", "c"}, {"ab": ("a", "b"), "bc": ("b", "c")})
        result = ECExpansion(graph, {"b": tree("oriented", "ab", "bc")}).embed()
        self.assertEqual(result.nodes, graph.nodes)
        self.assertEqual(result.edges, graph.edges)
        ECExpansion(Graph(), {}).embed().validate_embedding()

    def test_random_small_graphs_against_exhaustive_rotation_oracle(self):
        rng = random.Random(160)
        for index in range(36):
            nodes = list("abcde")
            pairs = list(itertools.combinations(nodes, 2))
            rng.shuffle(pairs)
            graph = graph_from_pairs(pairs[: rng.randrange(4, 9)])
            constraints = {}
            for node in sorted(graph.nodes):
                incident = list(graph.rotation[node])
                rng.shuffle(incident)
                if len(incident) > 1 and rng.random() < 0.7:
                    constraints[node] = tree(
                        rng.choice(["group", "mirror", "oriented"]), *incident
                    )
            expected = next(all_planar_rotations(graph, constraints), None)
            with self.subTest(index=index):
                try:
                    actual = ECExpansion(graph, constraints).embed()
                except NotECPlanar:
                    self.assertIsNone(expected)
                else:
                    self.assertIsNotNone(expected)
                    actual.validate_embedding()
                    for node, constraint in constraints.items():
                        self.assertTrue(admissible(constraint, actual.rotation[node]))


if __name__ == "__main__":
    unittest.main()
