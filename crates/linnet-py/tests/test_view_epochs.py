"""Element identity survives append/direction changes; mutation guards do not."""

import unittest

from symbolica.community import graph as g


def graph():
    return g.build(
        g.node("a", data=7),
        g.node("b", data=9),
        g.edge(g.source("a", data=1), "ab", g.sink("b", data=2), data=4),
    )


def mutate(value, operation):
    if operation == "node":
        value.add_node(g.node("c", data=11))
    elif operation == "edge":
        value.add_edge(g.edge(g.source("a", data=3), "out", data=5))
    elif operation == "orientation":
        value.set_orientation("ab", g.Orientation.Reversed)


class ViewEpochTests(unittest.TestCase):
    def test_map_each_callback_after_noninvalidating_mutations(self):
        for operation in ("none", "node", "edge", "orientation", "reordered"):
            for kind in ("node", "edge", "source", "sink"):
                with self.subTest(operation=operation, callback=kind):
                    value = graph()
                    if operation == "reordered":
                        value.reorder_nodes([1, 0])
                    retained = (
                        value.node("a"),
                        value.edge("ab"),
                        value.edge("ab").source,
                    )
                    mutate(value, operation)
                    self.assertEqual([v.data for v in retained], [7, 4, 1])
                    seen = []

                    def callback(view):
                        seen.append(view.data)
                        return (view.data, "mapped")

                    mapped = value.map(**{kind: callback})
                    if kind == "node":
                        inputs = value.nodes()
                    elif kind == "edge":
                        inputs = value.edges()
                    else:
                        inputs = [getattr(edge, kind) for edge in value.edges()]
                        inputs = [port for port in inputs if port is not None]
                    self.assertEqual(seen, [view.data for view in inputs])
                    self.assertEqual(mapped.n_nodes, value.n_nodes)
                    self.assertEqual(mapped.n_edges, value.n_edges)
                    if kind == "node":
                        self.assertEqual(mapped.node("a").data, (7, "mapped"))
                    elif kind == "edge":
                        self.assertEqual(mapped.edge("ab").data, (4, "mapped"))
                    else:
                        port = getattr(mapped.edge("ab"), kind)
                        original = getattr(value.edge("ab"), kind)
                        self.assertEqual(port.data, (original.data, "mapped"))
                    # The derived graph starts with its own consistent identity counters.
                    self.assertEqual(
                        mapped.map(node=lambda n: n.data).n_nodes, mapped.n_nodes
                    )
                    mapped.split_edge(mapped.half_edges()[0], g.EdgeValue(data=6))

    def test_split_connect_and_traversal_accept_fresh_views_and_indices(self):
        for operation in ("none", "node", "edge", "orientation", "reordered"):
            outcomes = []
            for views in (False, True):
                with self.subTest(operation=operation, views=views):
                    value = graph()
                    if operation == "reordered":
                        value.reorder_nodes([1, 0])
                    mutate(value, operation)
                    root = value.node("a")
                    for traverse in (
                        value.breadth_first_traverse,
                        value.depth_first_traverse,
                    ):
                        tree = traverse(root if views else root.index)
                        self.assertEqual([node.data for node in tree.nodes], [7, 9])
                        self.assertTrue(
                            all(port.data in (1, 2, 3) for port in tree.half_edges)
                        )
                    before = value.edge("ab").source
                    index = before.index
                    value.split_edge(before if views else index, g.EdgeValue(data=6))
                    with self.assertRaises(ReferenceError):
                        _ = before.data
                    ports = [value.half_edges()[index], value.half_edges()[1 - index]]
                    connected = value.connect(
                        *(ports if views else [p.index for p in ports]),
                        g.EdgeValue(data=8),
                    )
                    self.assertEqual(connected.data, 8)
                    with self.assertRaises(ReferenceError):
                        _ = ports[0].data
                    outcomes.append((value.n_nodes, value.n_edges, value.n_half_edges))
            self.assertEqual(len(outcomes), 2)
            self.assertEqual(*outcomes)

    def test_render_selectors_use_element_identity_after_mutation(self):
        for operation in ("node", "edge", "orientation"):
            for kind in ("node", "edge", "source", "sink"):
                with self.subTest(operation=operation, selector=kind):
                    value = graph()
                    mutate(value, operation)
                    seen = []

                    def select(view):
                        seen.append(view.data)
                        return None

                    value.render(
                        config=g.RenderSettings(
                            selectors=g.DrawingSelectors(**{kind: select})
                        )
                    )
                    self.assertTrue(seen)
        value = graph()

        def mutate_selector(view):
            data = view.data
            if value.n_nodes == 2:
                value.add_node(g.node("c"))
            return g.NodeDrawing(label=g.TextLabel(str(data)))

        with self.assertRaisesRegex(ReferenceError, "must not mutate graph topology"):
            value.render(
                config=g.RenderSettings(
                    selectors=g.DrawingSelectors(node=mutate_selector)
                )
            )

    def test_map_mutation_guard_and_exception_propagation(self):
        for operation in ("node", "edge", "orientation", "reorder"):
            with self.subTest(operation=operation):
                value = graph()
                changed = False

                def callback(view):
                    nonlocal changed
                    data = view.data
                    if not changed:
                        changed = True
                        if operation == "reorder":
                            value.reorder_nodes([1, 0])
                        else:
                            mutate(value, operation)
                    return data

                with self.assertRaises(ReferenceError):
                    value.map(node=callback)
        sentinel = LookupError("map sentinel")

        def fail(_view):
            raise sentinel

        with self.assertRaises(LookupError) as caught:
            graph().map(node=fail)
        self.assertIs(caught.exception, sentinel)

    def test_reordering_invalidates_retained_elements(self):
        value = graph()
        retained = (value.node("a"), value.edge("ab"), value.edge("ab").source)
        value.reorder_nodes([1, 0])
        for view in retained:
            with self.assertRaises(ReferenceError):
                _ = view.data
        with self.assertRaises(ReferenceError):
            value.split_edge(retained[-1], g.EdgeValue())
        with self.assertRaises(ReferenceError):
            value.breadth_first_traverse(retained[0])
        with self.assertRaises(ValueError):
            value.depth_first_traverse(graph().node("a"))
        self.assertEqual(value.node("a").data, 7)


if __name__ == "__main__":
    unittest.main()
