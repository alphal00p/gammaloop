"""Reusable descriptions and native insertion have distinct ownership contracts."""

import gc
import unittest
from functools import partial

from symbolica.community import graph as g


class ConstructionTests(unittest.TestCase):
    def test_specs_are_reusable_and_inferred_by_identity(self):
        for store in (g.NodeStore.Vec, g.NodeStore.Forest):
            with self.subTest(store=store):
                payload, port_payload = object(), object()
                a = g.Node("a", data=payload, label="original")
                b = g.Node(data=payload)
                port = g.HalfEdge(a, data=port_payload, compass=g.Compass.E)
                edge = g.Edge(port, b, name="ab", data=payload, label="edge")
                first = g.Graph(
                    edges=[edge, g.Edge(a, a, name="loop")], node_store=store
                )
                second = g.Graph(nodes=[a, b], edges=[edge], node_store=store)
                self.assertEqual((first.n_nodes, first.n_edges), (2, 2))
                self.assertEqual((second.n_nodes, second.n_edges), (2, 1))
                self.assertIsInstance(a, g.Node)
                self.assertIsInstance(first.node("a"), g.NodeView)
                self.assertIsInstance(first.edge("ab"), g.EdgeView)
                self.assertIsInstance(first.edge("ab").source, g.HalfEdgeView)
                self.assertIs(first.node("a").data, payload)
                self.assertIs(second.edge("ab").source.data, port_payload)
                self.assertEqual(first.edge("ab").source.drawing.compass, g.Compass.E)
                self.assertNotEqual(
                    first.edge("loop").source.index, first.edge("loop").sink.index
                )

                first.node("a").data = "first graph"
                first.node("a").drawing.label = "first graph"
                a.data = "next build"
                a.drawing.label = "next build"
                edge.drawing.label = "next edge"
                self.assertIs(second.node("a").data, payload)
                self.assertEqual(second.node("a").drawing.label, "original")
                self.assertEqual(second.edge("ab").drawing.label, "edge")
                third = g.Graph(edges=[edge], node_store=store)
                self.assertEqual(third.node("a").data, "next build")
                self.assertEqual(third.node("a").drawing.label, "next build")
                self.assertEqual(third.edge("ab").drawing.label, "next edge")

    def test_endpoint_slots_determine_flow_and_create_distinct_incidences(self):
        for store in (g.NodeStore.Vec, g.NodeStore.Forest):
            with self.subTest(store=store):
                a = g.Node("a")
                port = g.HalfEdge(a, data="shared recipe")
                graph = g.Graph(
                    edges=[
                        g.Edge(source=a, name="out"),
                        g.Edge(sink=port, name="in"),
                        g.Edge(
                            port,
                            port,
                            name="loop",
                            orientation=g.Orientation.Undirected,
                        ),
                    ],
                    node_store=store,
                )
                self.assertEqual((graph.n_nodes, graph.n_half_edges), (1, 4))
                self.assertIsNone(graph.edge("out").sink)
                self.assertEqual(graph.edge("out").source.flow, g.Flow.Source)
                self.assertIsNone(graph.edge("in").source)
                self.assertEqual(graph.edge("in").sink.flow, g.Flow.Sink)
                loop = graph.edge("loop")
                self.assertNotEqual(loop.source.index, loop.sink.index)
                self.assertEqual(loop.source.pair.index, loop.sink.index)
                self.assertEqual(loop.source.data, loop.sink.data)
                self.assertEqual(loop.orientation, g.Orientation.Undirected)
                with self.assertRaises(ValueError):
                    g.Edge()

    def test_keys_resolve_declared_nodes_and_do_not_create_them(self):
        for store in (g.NodeStore.Vec, g.NodeStore.Forest):
            with self.subTest(store=store):
                a, b = g.Node("a"), g.Node("b")
                graph = g.Graph(
                    nodes=(node for node in (a, b)),
                    edges=(edge for edge in (g.Edge("a", 1),)),
                    node_store=store,
                )
                self.assertEqual(graph.edge(0).source.node.name, "a")
                self.assertEqual(graph.edge(0).sink.node.name, "b")
                with self.assertRaises((KeyError, ValueError)):
                    g.Graph(edges=[g.Edge("missing")], node_store=store)
                with self.assertRaises((IndexError, ValueError)):
                    g.Graph(nodes=[a], edges=[g.Edge(1)], node_store=store)
                with self.assertRaisesRegex(ValueError, "duplicate node name"):
                    g.Graph(nodes=[a, g.Node("a")], node_store=store)
                with self.assertRaises(ValueError):
                    g.Graph(
                        nodes=[a], edges=[g.Edge(graph.node("a"))], node_store=store
                    )

    def test_incremental_views_and_drawing_objects_survive_append(self):
        for store in (g.NodeStore.Vec, g.NodeStore.Forest):
            with self.subTest(store=store):
                graph = g.Graph(node_store=store)
                payload = object()
                a = graph.add_node(name="a", data=payload, label="a")
                node_drawing = a.drawing
                b = graph.add_node("b", data=payload)
                edge = graph.add_edge(
                    a, g.HalfEdge(b, data="sink"), name="ab", data=payload
                )
                source, sink = edge.source, edge.sink
                edge_drawing, port_drawing = edge.drawing, source.drawing
                graph.add_node(g.Node("c"))
                graph.add_edge(a, b, name="parallel")
                graph.add_edge(source=a, name="external")
                self.assertIs(a.data, payload)
                self.assertIs(edge.data, payload)
                self.assertEqual(source.pair.index, sink.index)
                self.assertEqual(sink.data, "sink")
                self.assertEqual(len(a.incidence), 3)
                node_drawing.label = "retained node drawing"
                edge_drawing.label = "retained edge drawing"
                port_drawing.label = "retained port drawing"
                self.assertEqual(graph.node("a").drawing.label, "retained node drawing")
                self.assertEqual(
                    graph.edge("ab").drawing.label, "retained edge drawing"
                )
                self.assertEqual(
                    graph.edge("ab").source.drawing.label, "retained port drawing"
                )
                source.data = "live source"
                self.assertEqual(graph.edge("ab").source.data, "live source")

    def test_failed_insertions_leave_topology_and_all_existing_views_usable(self):
        for store in (g.NodeStore.Vec, g.NodeStore.Forest):
            with self.subTest(store=store):
                graph = g.Graph(node_store=store)
                a, b = graph.add_node("a"), graph.add_node("b")
                edge = graph.add_edge(a, b, name="ab")
                selected = graph.full_subgraph()
                revision = selected.revision
                drawing = edge.drawing
                foreign = g.Graph(nodes=[g.Node("foreign")]).node(0)
                failures = [
                    (ValueError, partial(graph.add_node, "a")),
                    (ValueError, partial(graph.add_edge, a, b, name="ab")),
                    ((KeyError, ValueError), partial(graph.add_edge, a, "missing")),
                    ((IndexError, ValueError), partial(graph.add_edge, a, 99)),
                    (ValueError, partial(graph.add_edge, a, foreign)),
                    (
                        (TypeError, ValueError),
                        partial(graph.add_node, g.Node("c"), data="conflict"),
                    ),
                    (
                        (TypeError, ValueError),
                        partial(graph.add_edge, g.Edge(a, b), data="conflict"),
                    ),
                    (
                        TypeError,
                        partial(graph.add_node, g.Node("c", data="payload"), data=None),
                    ),
                    (
                        TypeError,
                        partial(
                            graph.add_edge, g.Edge(a, b, data="payload"), data=None
                        ),
                    ),
                    (
                        (TypeError, ValueError),
                        partial(graph.add_node, "c", radius="bad radius"),
                    ),
                ]
                for expected, insert in failures:
                    with self.subTest(insert=insert), self.assertRaises(expected):
                        insert()
                    self.assertEqual(
                        (graph.n_nodes, graph.n_edges, graph.n_half_edges), (2, 1, 2)
                    )
                    self.assertEqual(selected.revision, revision)
                    self.assertEqual(selected.half_edge_indices(), [0, 1])
                    self.assertEqual((a.name, b.name, edge.name), ("a", "b", "ab"))
                    drawing.label = "still live"
                    self.assertEqual(graph.edge("ab").drawing.label, "still live")

    def test_incremental_payload_defaults_preserve_explicit_values(self):
        for store in (g.NodeStore.Vec, g.NodeStore.Forest):
            with self.subTest(store=store):
                graph = g.Graph(node_store=store)
                node = graph.add_node()
                self.assertIsNone(node.data)
                self.assertIsNone(graph.add_edge(node).data)
                for payload in (None, Ellipsis):
                    self.assertIs(graph.add_node(data=payload).data, payload)
                    self.assertIs(graph.add_edge(node, data=payload).data, payload)
                payload = object()
                self.assertIs(graph.add_node(g.Node(data=payload)).data, payload)
                self.assertIs(graph.add_edge(g.Edge(node, data=payload)).data, payload)

    def test_append_invalidates_topology_results_but_not_returned_element_views(self):
        for store in (g.NodeStore.Vec, g.NodeStore.Forest):
            for append in ("node", "edge"):
                with self.subTest(store=store, append=append):
                    a, b, c = (g.Node(name) for name in "abc")
                    graph = g.Graph(
                        edges=[g.Edge(a, b), g.Edge(b, c), g.Edge(c, a)],
                        node_store=store,
                    )
                    selected = graph.full_subgraph()
                    cycle = graph.cycle_basis()[0][0]
                    tree = graph.depth_first_traverse("a")
                    partition = graph.all_cuts(["a"], ["b"])[0]
                    cut = partition.boundary
                    node_view, edge_view = tree.nodes[0], cut.edges[0]
                    if append == "node":
                        graph.add_node("later")
                    else:
                        graph.add_edge("a", "b")
                    self.assertEqual(node_view.name, "a")
                    self.assertIsInstance(edge_view.source, g.HalfEdgeView)
                    for access in (
                        selected.half_edge_indices,
                        partial(getattr, cycle, "filter"),
                        partial(getattr, tree, "nodes"),
                        partial(getattr, partition, "source_side"),
                        partial(getattr, cut, "edges"),
                    ):
                        with self.assertRaises(ReferenceError):
                            access()

    def test_appending_in_callbacks_is_still_detected(self):
        for store in (g.NodeStore.Vec, g.NodeStore.Forest):
            for operation in ("filter", "map", "render"):
                with self.subTest(store=store, operation=operation):
                    graph = g.Graph(nodes=[g.Node("a")], node_store=store)
                    view = graph.node("a")

                    def mutate(_node, graph=graph, operation=operation):
                        graph.add_node()
                        return g.NodeDrawing() if operation == "render" else True

                    with self.assertRaises(ReferenceError):
                        if operation == "render":
                            graph.render(
                                config=g.RenderSettings(
                                    selectors=g.DrawingSelectors(node=mutate)
                                )
                            )
                        else:
                            getattr(graph, operation)(node=mutate)
                    self.assertEqual(view.name, "a")

    def test_payload_finalizers_run_after_graph_borrows_are_released(self):
        for store in (g.NodeStore.Vec, g.NodeStore.Forest):
            with self.subTest(store=store):
                graph = g.Graph(node_store=store)
                finalized = []

                class Payload:
                    def __init__(self, graph, finalized, label):
                        self.graph = graph
                        self.finalized = finalized
                        self.label = label

                    def __del__(self):
                        self.graph.add_node(name=self.label)
                        self.finalized.append(self.label)

                node = graph.add_node(data=Payload(graph, finalized, "node finalizer"))
                edge = graph.add_edge(
                    node, data=Payload(graph, finalized, "edge finalizer")
                )
                half = edge.source
                half.data = Payload(graph, finalized, "half-edge finalizer")
                node.data = None
                edge.data = None
                half.data = None
                gc.collect()
                self.assertEqual(
                    finalized,
                    ["node finalizer", "edge finalizer", "half-edge finalizer"],
                )
                self.assertEqual(graph.n_nodes, 4)

    def test_index_conversion_can_reenter_graph_without_borrowing_panics(self):
        class Index:
            def __init__(self, callback=None):
                self.callback = callback

            def __index__(self):
                if self.callback is not None:
                    self.callback()
                return 0

        for store in (g.NodeStore.Vec, g.NodeStore.Forest):
            with self.subTest(store=store):
                graph = g.Graph(node_store=store)
                node = graph.add_node("a")
                graph.add_edge(node, name="external")
                for kind, expected in (("node", "a"), ("edge", "external")):
                    with self.subTest(kind=kind):
                        added = f"{kind} lookup callback"
                        key = Index(
                            lambda added=added, graph=graph: graph.add_node(added)
                        )
                        result = getattr(graph, kind)(key)
                        self.assertEqual(result.name, expected)
                        self.assertEqual(graph.node(added).name, added)
                self.assertEqual((graph.n_nodes, graph.n_edges), (3, 1))

                # Prepare the endpoint before enabling its callback so that
                # reentry happens during insertion's endpoint resolution.
                key = Index()
                endpoint = g.HalfEdge(key)
                key.callback = lambda graph=graph: graph.add_node("endpoint callback")
                with self.assertRaises(ReferenceError):
                    graph.add_edge(endpoint, node, name="not inserted")
                self.assertEqual((graph.n_nodes, graph.n_edges), (4, 1))
                self.assertEqual(
                    graph.node("endpoint callback").name, "endpoint callback"
                )
                self.assertEqual(node.name, "a")
                self.assertEqual(graph.edge(0).name, "external")

    def test_spec_payload_finalizers_can_read_the_replacement(self):
        for specification in (g.Node("a"), g.HalfEdge("a"), g.Edge("a")):
            with self.subTest(specification=type(specification).__name__):
                observed, errors = [], []

                class Payload:
                    def __init__(self, specification, observed, errors):
                        self.specification = specification
                        self.observed = observed
                        self.errors = errors

                    def __del__(self):
                        try:
                            self.observed.append(self.specification.data)
                        except RuntimeError as error:
                            self.errors.append(error)

                specification.data = Payload(specification, observed, errors)
                specification.data = "replacement"
                gc.collect()
                self.assertEqual(errors, [])
                self.assertEqual(observed, ["replacement"])


if __name__ == "__main__":
    unittest.main()
