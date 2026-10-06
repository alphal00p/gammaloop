"""Graphica algorithms preserve the public half-edge graph and its payloads."""

import copy
from collections import Counter
from itertools import permutations
import unittest
from symbolica import S
from symbolica.community import graph as g


def link(a, b, *, data=0, orientation=g.Orientation.Undirected, **kwargs):
    return g.edge(
        g.source(a), second=g.sink(b), data=data, orientation=orientation, **kwargs
    )


class AlgorithmsTests(unittest.TestCase):
    def test_graphica_generation_fixtures(self):
        gluon = g.EdgeSignature(S("g"))
        quark = g.EdgeSignature(S("q"), direction=False)
        for external, signatures, loops, expected in (
            (2, [[gluon] * 3, [quark.flip(), quark, gluon], [gluon] * 4], 3, 210),
            (4, [[gluon] * 3, [gluon] * 4], 2, 278),
        ):
            with self.subTest(external=external, loops=loops):
                graphs = g.Graph.generate(
                    [(i + 1, gluon) for i in range(external)],
                    signatures,
                    max_loops=loops,
                    max_bridges=0,
                    allow_self_loops=True,
                )
                self.assertEqual(len(graphs), expected)
                for graph, size in graphs:
                    self.assertEqual(graph.canonize()[2], size)

    def test_graphica_canonical_symmetry_fixture(self):
        # The colored nine-vertex fixture from Graphica's canonicalization suite.
        labels = [1, 0, 1, 0, 2, 0, 1, 0, 1]
        pairs = [
            (0, 1),
            (0, 3),
            (1, 2),
            (1, 3),
            (1, 4),
            (1, 5),
            (2, 5),
            (3, 4),
            (3, 6),
            (3, 7),
            (4, 5),
            (4, 7),
            (5, 7),
            (5, 8),
            (6, 7),
            (7, 8),
        ]
        graph = g.build(
            *(g.node(str(i), data=value) for i, value in enumerate(labels)),
            *(link(str(a), str(b)) for a, b in pairs),
        )
        canonical, mapping, size, orbit = graph.canonize()
        self.assertEqual(size, 8)
        self.assertEqual(sorted(Counter(orbit).values()), [1, 4, 4])
        for old, new in enumerate(mapping):
            self.assertEqual(canonical.node(new).data, labels[old])

    def test_graphica_generation_filter_fixture(self):
        gluon = g.EdgeSignature(S("g"))
        quark = g.EdgeSignature(S("q"), direction=False)

        def at_most_one_cubic_gluon_vertex(graph, completed):
            self.assertIs(type(graph), g.Graph)
            count = 0
            for node in graph.nodes()[:completed]:
                # Graphica counts each self-loop once in the incident edge list.
                edges = {port.edge.index: port.edge for port in node.incidence}
                if len(edges) == 3 and all(
                    edge.data == gluon.data for edge in edges.values()
                ):
                    count += 1
                    if count == 2:
                        return False
            return True

        graphs = g.Graph.generate(
            [(1, gluon), (2, gluon)],
            [[gluon] * 3, [quark.flip(), quark, gluon], [gluon] * 4],
            max_loops=4,
            max_bridges=0,
            allow_self_loops=True,
            filter_fn=at_most_one_cubic_gluon_vertex,
        )
        self.assertEqual(len(graphs), 845)

    def test_graphica_permutation_fixture(self):
        pairs = [(0, 0), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
        expected = None
        for order in permutations(range(5)):
            graph = g.build(
                *(g.node(str(i)) for i in range(5)),
                *(link(str(order[a]), str(order[b])) for a, b in pairs),
            ).canonize()[0]
            topology = [
                (edge.source.node.index, edge.sink.node.index) for edge in graph.edges()
            ]
            if expected is None:
                expected = topology
            self.assertEqual(topology, expected)

    def test_generation_returns_common_graphs_and_integer_group_sizes(self):
        signature = g.EdgeSignature(S("phi"))
        generated = g.Graph.generate(
            [(1, signature), (2, signature)], [[signature] * 3], max_loops=0
        )
        self.assertIsInstance(generated, list)
        for graph, size in generated:
            self.assertIs(type(graph), g.Graph)
            self.assertIsInstance(size, int)
            self.assertGreaterEqual(size, 1)
        # Four external legs in phi^3 have exactly the s, t, and u trees.
        generated = g.Graph.generate(
            [(i, signature) for i in range(1, 5)], [[signature] * 3], max_loops=0
        )
        self.assertEqual(len(generated), 3)
        self.assertEqual([size for _, size in generated], [1, 1, 1])

    def test_parallel_edges_and_self_loops_have_expected_symmetries(self):
        parallel = g.build(
            g.node("a", data=1), g.node("b", data=2), link("a", "b"), link("a", "b")
        )
        self.assertEqual(parallel.canonize()[2], 2)
        loop = g.build(g.node("a"), link("a", "a"))
        self.assertEqual(loop.canonize()[2], 2)
        directed = g.build(
            g.node("a"), link("a", "a", orientation=g.Orientation.Default)
        )
        self.assertEqual(directed.canonize()[2], 1)

    def test_dangling_half_edges_and_isolated_nodes(self):
        a = g.build(g.node("a"), g.node("isolated"), g.edge(g.source("a"), data=1))
        b = g.build(g.node("isolated"), g.node("b"), g.edge(g.source("b"), data=1))
        self.assertTrue(a.is_isomorphic(b))
        canonical, mapping, size, orbit = a.canonize()
        self.assertEqual(canonical.n_nodes, 2)
        self.assertEqual(canonical.n_half_edges, 1)
        self.assertEqual(sorted(mapping), [0, 1])
        self.assertEqual(len(orbit), 2)
        self.assertTrue(all(0 <= v < 2 for v in orbit))
        self.assertEqual(size, 1)

    def test_custom_keys_preserve_payload_identity_and_drawing(self):
        payload = {"label": 4}
        graph = g.build(
            g.node("a", data=payload, label="custom"),
            g.node("b", data={"label": 3}),
            link("a", "b"),
        )
        with self.assertRaises(TypeError):
            graph.canonize()
        result, mapping, _, _ = graph.canonize(node_key=lambda data: data["label"])
        self.assertIs(result.node(mapping[0]).data, payload)
        self.assertEqual(result.node(mapping[0]).drawing.label, "custom")
        self.assertTrue(
            graph.is_isomorphic(result, node_key=lambda data: data["label"])
        )

    def test_drawing_metadata_does_not_color_the_graph(self):
        a = g.build(g.node("a", label="first"), g.node("b"), link("a", "b"))
        b = g.build(g.node("x", label="different"), g.node("y"), link("x", "y"))
        self.assertTrue(a.is_isomorphic(b))
        self.assertEqual(a.canonize()[2], 2)

    def test_half_edge_keys_and_direction(self):
        a = g.build(
            g.node("a", data=1),
            g.node("b", data=2),
            link("a", "b", orientation=g.Orientation.Default),
        )
        b = g.build(
            g.node("a", data=1),
            g.node("b", data=2),
            link("b", "a", orientation=g.Orientation.Reversed),
        )
        self.assertTrue(a.is_isomorphic(b))
        a.half_edges()[0].data = {"port": 7}
        with self.assertRaises(TypeError):
            a.canonize()
        canonical, _, _, _ = a.canonize(half_edge_key=lambda d: d["port"] if d else 0)
        self.assertTrue(
            a.is_isomorphic(canonical, half_edge_key=lambda d: d["port"] if d else 0)
        )

    def test_copy_and_canonicalization_do_not_stale_original_views(self):
        original = g.build(g.node("a"), g.node("b"), link("a", "b"))
        view = original.node(0)
        cloned = copy.copy(original)
        original.canonize()
        self.assertIsNone(view.data)
        cloned.add_node(g.node("new"))
        self.assertEqual(original.n_nodes, 2)

    def test_callback_errors_propagate_and_generation_recovers(self):
        signature = g.EdgeSignature(0)

        def reject(graph, count):
            self.assertIs(type(graph), g.Graph)
            raise RuntimeError("callback sentinel")

        kwargs = dict(
            external_edges=[(i, signature) for i in range(4)],
            vertex_signatures=[[signature] * 3],
            max_loops=0,
        )
        with self.assertRaisesRegex(RuntimeError, "callback sentinel"):
            g.Graph.generate(**kwargs, filter_fn=reject)
        self.assertEqual(len(g.Graph.generate(**kwargs)), 3)

    def test_edge_canonicalization_ignores_fixed_vertex_payloads(self):
        payload = object()
        graph = g.build(
            g.node("a", data=payload),
            g.node("b", data=object()),
            link("a", "b", data=2),
            link("a", "b", data=1),
        )
        graph.canonize_edges()
        self.assertIs(graph.node("a").data, payload)
        self.assertEqual(sorted(edge.data for edge in graph.edges()), [1, 2])

    def test_mermaid_escapes_labels(self):
        graph = g.build(g.node("a", data='<bad "label">'))
        self.assertIn("&lt;bad &quot;label&quot;&gt;", graph.to_mermaid())


if __name__ == "__main__":
    unittest.main()
