"""Generation and canonization preserve the public half-edge graph and its payloads."""

import copy
import gc
import unittest
import weakref
from collections import Counter
from itertools import permutations

from symbolica import E, S
from symbolica.community import graph as g


def link(a, b, *, data=0, orientation=g.Orientation.Undirected, **kwargs):
    return g.Edge(source=a, sink=b, data=data, orientation=orientation, **kwargs)


class AlgorithmsTests(unittest.TestCase):
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
        graph = g.Graph(
            nodes=[g.Node(str(i), data=value) for i, value in enumerate(labels)],
            edges=[link(str(a), str(b)) for a, b in pairs],
        )
        canonical, mapping, size, orbit = graph.canonize()
        self.assertEqual(size, 8)
        self.assertEqual(sorted(Counter(orbit).values()), [1, 4, 4])
        for old, new in enumerate(mapping):
            self.assertEqual(canonical.node(new).data, labels[old])

    def test_graphica_permutation_fixture(self):
        pairs = [(0, 0), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
        expected = None
        for order in permutations(range(5)):
            graph = g.Graph(
                nodes=[g.Node(str(i)) for i in range(5)],
                edges=[link(str(order[a]), str(order[b])) for a, b in pairs],
            ).canonize()[0]
            topology = [
                (edge.source.node.index, edge.sink.node.index) for edge in graph.edges()
            ]
            if expected is None:
                expected = topology
            self.assertEqual(topology, expected)

    def test_parallel_edges_and_self_loops_have_expected_symmetries(self):
        parallel = g.Graph(
            nodes=[g.Node("a", data=1), g.Node("b", data=2)],
            edges=[link("a", "b"), link("a", "b")],
        )
        self.assertEqual(parallel.canonize()[2], 2)
        loop = g.Graph(nodes=[g.Node("a")], edges=[link("a", "a")])
        self.assertEqual(loop.canonize()[2], 2)
        directed = g.Graph(
            nodes=[g.Node("a")],
            edges=[link("a", "a", orientation=g.Orientation.Default)],
        )
        self.assertEqual(directed.canonize()[2], 1)

    def test_dangling_half_edges_and_isolated_nodes(self):
        a = g.Graph(
            nodes=[g.Node("a"), g.Node("isolated")], edges=[g.Edge(source="a", data=1)]
        )
        b = g.Graph(
            nodes=[g.Node("isolated"), g.Node("b")], edges=[g.Edge(source="b", data=1)]
        )
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
        graph = g.Graph(
            nodes=[
                g.Node("a", data=payload, label="custom"),
                g.Node("b", data={"label": 3}),
            ],
            edges=[link("a", "b")],
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
        a = g.Graph(
            nodes=[g.Node("a", label="first"), g.Node("b")], edges=[link("a", "b")]
        )
        b = g.Graph(
            nodes=[g.Node("x", label="different"), g.Node("y")], edges=[link("x", "y")]
        )
        self.assertTrue(a.is_isomorphic(b))
        self.assertEqual(a.canonize()[2], 2)

    def test_half_edge_keys_and_direction(self):
        a = g.Graph(
            nodes=[g.Node("a", data=1), g.Node("b", data=2)],
            edges=[link("a", "b", orientation=g.Orientation.Default)],
        )
        b = g.Graph(
            nodes=[g.Node("a", data=1), g.Node("b", data=2)],
            edges=[link("b", "a", orientation=g.Orientation.Reversed)],
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
        original = g.Graph(nodes=[g.Node("a"), g.Node("b")], edges=[link("a", "b")])
        view = original.node(0)
        cloned = copy.copy(original)
        original.canonize()
        self.assertIsNone(view.data)
        cloned.add_node(g.Node("new"))
        self.assertEqual(original.n_nodes, 2)

    def test_edge_canonicalization_ignores_fixed_vertex_payloads(self):
        payload = object()
        graph = g.Graph(
            nodes=[g.Node("a", data=payload), g.Node("b", data=object())],
            edges=[link("a", "b", data=2), link("a", "b", data=1)],
        )
        graph.canonize_edges()
        self.assertIs(graph.node("a").data, payload)
        self.assertEqual(sorted(edge.data for edge in graph.edges()), [1, 2])

    def test_mermaid_escapes_labels(self):
        graph = g.Graph(nodes=[g.Node("a", data='<bad "label">')])
        self.assertIn("&lt;bad &quot;label&quot;&gt;", graph.to_mermaid())


GLUON = g.EdgeSignature(S("g"))
QUARK = g.EdgeSignature(S("q"), direction=False)
QCD = [[GLUON] * 3, [QUARK.flip(), QUARK, GLUON], [GLUON] * 4]


def legs(graph):
    """`(node, label, orientation)` of every dangling half-edge."""
    return sorted(
        (
            (half.node.index, half.data, edge.orientation)
            for edge in graph.edges()
            for half in (edge.source, edge.sink)
            if half is not None and half.pair is None
        ),
        key=repr,
    )


class GenerationTests(unittest.TestCase):
    def test_counts_and_symmetry_factors_match_graphica(self):
        for external, loops, expected in ((2, range(4), 210), (4, range(3), 278)):
            with self.subTest(external=external):
                graphs = list(
                    g.Graph.generate(
                        [(i + 1, GLUON) for i in range(external)],
                        [QCD[0], QCD[2]] if external == 4 else QCD,
                        loops=loops,
                        max_bridges=0,
                        allow_self_loops=True,
                    )
                )
                self.assertEqual(len(graphs), expected)
                for graph, size in graphs:
                    self.assertEqual(graph.canonize()[2], size)

    def test_graphs_are_native_with_dangling_labelled_legs(self):
        phi = g.EdgeSignature("phi")
        generated = g.Graph.generate(
            [(i, phi) for i in range(1, 5)], [[phi] * 3], loops=0
        )
        self.assertIsInstance(generated, g.GeneratedGraphs)
        self.assertIs(iter(generated), generated)
        graphs = list(generated)
        # Four distinct legs in phi^3 have exactly the s, t and u trees.
        self.assertEqual([size for _, size in graphs], [1, 1, 1])
        for graph, size in graphs:
            self.assertIs(type(graph), g.Graph)
            self.assertIs(type(size), int)
            self.assertEqual(graph.n_nodes, 2)
            self.assertEqual([node.data for node in graph.nodes()], [0, 0])
            self.assertTrue(all(edge.data == "phi" for edge in graph.edges()))
            self.assertEqual(sorted(label for _, label, _ in legs(graph)), [1, 2, 3, 4])
            self.assertTrue(all(o == g.Orientation.Undirected for *_, o in legs(graph)))
            for edge in graph.edges():
                if edge.source is not None and edge.sink is not None:
                    self.assertIsNone(edge.source.data)
        self.assertEqual(list(generated), [])

    def test_directed_fields_keep_their_slot_directions(self):
        photon = g.EdgeSignature("a")
        electron = g.EdgeSignature("e", direction=False)
        qed = [[electron, electron.flip(), photon]]
        ((graph, size),) = g.Graph.generate(
            [(1, photon), (2, photon)], qed, loops=1, max_bridges=0
        )
        self.assertEqual(size, 1)
        fermions = [e for e in graph.edges() if e.data == "e"]
        self.assertEqual(len(fermions), 2)
        self.assertTrue(all(e.orientation == g.Orientation.Default for e in fermions))
        # e+ e- annihilation: the incoming electron slot is a leg into its vertex.
        graphs = list(
            g.Graph.generate(
                [("e-", electron), ("e+", electron.flip()), (3, photon)], qed, loops=0
            )
        )
        self.assertEqual(len(graphs), 1)
        flows = {
            half.data: half.flow
            for e in graphs[0][0].edges()
            for half in (e.source, e.sink)
            if half is not None and half.pair is None
        }
        self.assertEqual(
            flows, {"e-": g.Flow.Sink, "e+": g.Flow.Source, 3: g.Flow.Source}
        )

    def test_signatures_keep_any_hashable_field(self):
        field = (1, "x")
        signature = g.EdgeSignature(field, direction=True)
        self.assertIs(signature.data, field)
        self.assertEqual(signature.direction, True)
        self.assertIs(signature.flip().data, field)
        self.assertEqual(signature.flip().direction, False)
        self.assertEqual(repr(signature), "EdgeSignature((1, 'x'), direction=True)")
        # Fields compare like canonization keys: 1, "1" and the expression 1 never pair.
        one = g.EdgeSignature(1)

        def count(*fields):
            legs = [(i, g.EdgeSignature(f)) for i, f in enumerate(fields)]
            return len(list(g.Graph.generate(legs, [[one] * 3], loops=1)))

        self.assertEqual(count(1, 1.0), 1)
        self.assertEqual(count(1, "1"), 0)
        self.assertEqual(count(1, E("1")), 0)
        with self.assertRaisesRegex(TypeError, "use hashable fields"):
            g.Graph.generate([], [[g.EdgeSignature([1])] * 3], loops=0)

    def test_filter_sees_every_stage_and_prunes_hereditarily(self):
        stages = set()

        def generate(keep):
            return len(
                list(
                    g.Graph.generate(
                        [(1, GLUON), (2, GLUON)],
                        QCD,
                        loops=range(5),
                        max_bridges=0,
                        allow_self_loops=True,
                        filter_fn=keep,
                    )
                )
            )

        def at_most_one_cubic_gluon_vertex(view):
            # Hereditary: rules(v) is [0] only once the fields at v fix the rule.
            self.assertIs(type(view), g.GenerationView)
            stages.add(view.stage)
            return sum(view.rules(v) == [0] for v in range(view.n_vertices)) <= 1

        def on_complete_graphs(condition):
            complete = g.GenerationStage.Complete
            return lambda view: view.stage != complete or condition(view)

        def three_gluon_neighbours(view):
            # Graphica's fixture: at most one vertex with three distinct incident edges,
            # legs included and a self-loop counted once, all gluons.
            count = 0
            for v in range(view.n_vertices):
                fields = [s.data for a, b, s in view.edges() if v in (a, b)]
                fields += [s.data for _, s in view.legs(v)]
                count += len(fields) == 3 and all(f == GLUON.data for f in fields)
            return count <= 1

        pruned = generate(at_most_one_cubic_gluon_vertex)
        self.assertEqual(
            stages,
            {
                g.GenerationStage.Topology,
                g.GenerationStage.Fields,
                g.GenerationStage.Complete,
            },
        )
        self.assertEqual(
            pruned, generate(on_complete_graphs(at_most_one_cubic_gluon_vertex))
        )
        self.assertEqual(pruned, 960)
        self.assertEqual(generate(on_complete_graphs(three_gluon_neighbours)), 845)

    def test_views_describe_the_search_state(self):
        seen = []

        def check(view):
            n = view.n_vertices
            edges = view.edges()
            ends = [sum(e[:2].count(v) for e in edges) for v in range(n)]
            for v in range(n):
                attached = len(view.legs(v))
                self.assertEqual(
                    view.degree(v), attached + ends[v] + view.open_stubs(v)
                )
            placed = sum(len(view.legs(v)) for v in range(n))
            self.assertEqual(placed + view.unplaced_legs, 2)
            self.assertEqual(len(edges) + 1, n + view.loops)
            self.assertEqual(view.loop_order, 2)
            for outside, inside, side in view.bridges():
                self.assertIn(inside, side)
                self.assertNotIn(outside, side)
                self.assertLessEqual(view.loops_in(side), view.loops)
            # loops_in needs a non-empty connected set, such as every placed vertex
            with self.assertRaisesRegex(ValueError, "connected"):
                view.loops_in([])
            if n >= 1:
                self.assertEqual(view.loops_in(range(n)), view.loops)
            for u in range(n):
                for v in range(u + 1, n):
                    if view.multiplicity(u, v) == 0:
                        with self.assertRaisesRegex(ValueError, "connected"):
                            view.loops_in([u, v])
            graph = view.graph()
            self.assertEqual(graph.n_nodes, n)
            self.assertEqual(graph.n_edges, len(edges) + placed)
            complete = view.complete_graph()
            if view.stage == g.GenerationStage.Complete:
                self.assertTrue(all(signature is not None for *_, signature in edges))
                self.assertTrue(complete.is_isomorphic(graph))
                seen.append(complete)
            else:
                self.assertIsNone(complete)
            if view.stage != g.GenerationStage.Topology:
                self.assertTrue(view.is_topology_complete)
            with self.assertRaises(IndexError):
                view.degree(n)
            if n >= 1:
                with self.assertRaises(ValueError):
                    view.multiplicity(0, 0)
            return True

        generated = list(
            g.Graph.generate(
                [(1, GLUON), (2, GLUON)], QCD, loops=2, max_bridges=0, filter_fn=check
            )
        )
        self.assertEqual(len(seen), len(generated))
        for (graph, _), complete in zip(generated, seen):
            self.assertTrue(graph.is_isomorphic(complete))

    def test_views_expire_with_their_filter_call(self):
        views = []

        def keep(view):
            views.append(view)
            return True

        list(g.Graph.generate([(1, GLUON), (2, GLUON)], QCD, loops=1, filter_fn=keep))
        self.assertTrue(views)
        with self.assertRaisesRegex(ReferenceError, "filter call"):
            views[0].degree(0)

    def test_generation_streams_and_progress_can_stop_it(self):
        reports = []

        def progress(loops, topologies, graphs):
            reports.append((loops, topologies, graphs))
            return graphs < 2

        generated = g.Graph.generate(
            [(1, GLUON), (2, GLUON)], QCD, loops=range(4), progress_fn=progress
        )
        next(generated)
        # Lazy: the first step searched only up to the first graph.
        self.assertEqual(max(r[2] for r in reports), 1)
        self.assertEqual(len(list(generated)), 1)
        self.assertTrue(all(isinstance(x, int) for r in reports for x in r))

    def test_callback_errors_propagate_and_generation_recovers(self):
        signature = g.EdgeSignature(0)

        def reject(view):
            raise RuntimeError("callback sentinel")

        kwargs = {
            "external_edges": [(i, signature) for i in range(4)],
            "vertex_signatures": [[signature] * 3],
            "loops": 0,
        }
        generated = g.Graph.generate(**kwargs, filter_fn=reject)
        with self.assertRaisesRegex(RuntimeError, "callback sentinel"):
            next(generated)
        self.assertEqual(list(generated), [])
        with self.assertRaisesRegex(TypeError, "bool"):
            next(g.Graph.generate(**kwargs, progress_fn=lambda *counts: None))
        self.assertEqual(len(list(g.Graph.generate(**kwargs))), 3)

    def test_invalid_bounds_are_rejected(self):
        phi = g.EdgeSignature(0)
        with self.assertRaisesRegex(ValueError, "bound the loop order"):
            g.Graph.generate([(1, phi), (2, phi)], [[phi] * 3])
        with self.assertRaisesRegex(ValueError, "step 1"):
            g.Graph.generate([(1, phi), (2, phi)], [[phi] * 3], loops=range(0, 4, 2))
        with self.assertRaises(OverflowError):
            g.Graph.generate([(1, phi), (2, phi)], [[phi] * 3], loops=-1)
        with self.assertRaisesRegex(TypeError, "not a bool"):
            g.Graph.generate([(1, phi), (2, phi)], [[phi] * 3], loops=True)

    def test_automorphism_orders_beyond_u128_raise(self):
        # Permuting 70 interchangeable legs at their vertices gives more than 2**128 symmetries.
        phi = g.EdgeSignature(0)
        generated = g.Graph.generate(
            [(0, phi)] * 70, [[phi] * 15, [phi] * 5], loops=0, max_vertices=6
        )
        with self.assertRaisesRegex(OverflowError, "2\\*\\*128"):
            next(generated)

    def test_callbacks_referring_to_their_iterator_are_collected(self):
        phi = g.EdgeSignature(0)

        class Marker:
            pass

        def cycle(callback):
            marker = Marker()
            holder = {"marker": marker}

            def keep(*args):
                return holder is not None

            holder["generated"] = g.Graph.generate(
                [(i, phi) for i in range(4)], [[phi] * 3], loops=0, **{callback: keep}
            )
            next(holder["generated"])
            return weakref.ref(marker)

        for callback in ("filter_fn", "progress_fn"):
            with self.subTest(callback=callback):
                probe = cycle(callback)
                gc.collect()
                self.assertIsNone(probe())


if __name__ == "__main__":
    unittest.main()
