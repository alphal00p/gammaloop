import copy
import gc
import json
import os
import re
import shutil
import unittest
import weakref
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import linnet_py as lp
import typst


class Payload:
    pass


@dataclass(eq=False)
class UserNodePayload:
    identifier: str


@dataclass(eq=False)
class UserEdgePayload:
    coupling: str


@dataclass(eq=False)
class UserHalfEdgePayload:
    port: str


def replace_linnest_native(dot, statement, value):
    envelope = {
        "schema": "linnet-typst-native",
        "version": 1,
        "value": value,
    }
    replacement = json.dumps(envelope, separators=(",", ":")).encode().hex()
    match = re.search(rf"{re.escape(statement)}\s*=\s*\"?([0-9a-fA-F]+)\"?", dot)
    if match is None:
        raise AssertionError(f"missing Linnest native statement {statement!r}")
    return dot[: match.start(1)] + replacement + dot[match.end(1) :]


def sample_graph(*, codec=None, render_config=None, node_store=lp.NodeStore.Vec):
    node_data = [Payload(), Payload()]
    edge_data = [Payload(), Payload()]
    half_edge_data = [Payload(), Payload(), Payload()]
    left = lp.node(
        "left",
        data=node_data[0],
        label=lp.MathSymbol("n", subscript=0),
        extensions={"role": "initial"},
    )
    right = lp.node(
        "right",
        data=node_data[1],
        label=lp.MathSymbol("n", subscript=1),
    )
    propagator = lp.edge(
        lp.source(
            left,
            data=half_edge_data[0],
            compass=lp.Compass.E,
            anchor=lp.Anchor.East,
            routing=lp.Routing.Direct,
        ),
        "propagator",
        lp.sink(right, data=half_edge_data[1], compass=lp.Compass.W),
        data=edge_data[0],
        label=lp.MathSymbol("p", subscript=0),
        routing=lp.Routing.HobbyThrough,
        extensions={"category": "primary"},
    )
    outgoing = lp.edge(
        lp.source(right, data=half_edge_data[2]),
        "outgoing",
        data=edge_data[1],
        extensions={"category": "outgoing"},
    )
    graph = lp.build(
        left,
        right,
        propagator,
        outgoing,
        name="example",
        codec=codec,
        render_config=render_config,
        node_store=node_store,
    )
    return graph, node_data, edge_data, half_edge_data


def topology_graph(*, node_store=lp.NodeStore.Vec):
    node_data = {name: UserNodePayload(name) for name in ("a", "b", "c", "d", "e", "f")}
    edge_data = {
        name: UserEdgePayload(f"coupling-{name}")
        for name in ("ab", "bc", "ca", "cd", "ef")
    }
    half_edge_data = {
        (name, role): UserHalfEdgePayload(f"{name}-{role}")
        for name in edge_data
        for role in ("source", "sink")
    }
    nodes = {name: lp.node(name, data=data) for name, data in node_data.items()}

    def internal_edge(name, source_name, sink_name):
        return lp.edge(
            lp.source(nodes[source_name], data=half_edge_data[(name, "source")]),
            name,
            lp.sink(nodes[sink_name], data=half_edge_data[(name, "sink")]),
            data=edge_data[name],
            extensions={"channel": f"channel-{name}"},
        )

    graph = lp.build(
        *nodes.values(),
        internal_edge("ab", "a", "b"),
        internal_edge("bc", "b", "c"),
        internal_edge("ca", "c", "a"),
        internal_edge("cd", "c", "d"),
        internal_edge("ef", "e", "f"),
        name="topology",
        node_store=node_store,
    )
    return graph, node_data, edge_data, half_edge_data


def dangling_graph(name, flow, *, node_store=lp.NodeStore.Vec):
    node_data = UserNodePayload(f"node-{name}")
    edge_data = UserEdgePayload(f"edge-{name}")
    half_edge_data = UserHalfEdgePayload(f"port-{name}")
    node_spec = lp.node(f"node-{name}", data=node_data)
    endpoint = (
        lp.source(node_spec, data=half_edge_data)
        if flow == lp.Flow.Source
        else lp.sink(node_spec, data=half_edge_data)
    )
    graph = lp.build(
        node_spec,
        lp.edge(endpoint, f"edge-{name}", data=edge_data),
        name=name,
        node_store=node_store,
    )
    return graph, node_data, edge_data, half_edge_data


def round_trip_codec(seen=None):
    if seen is None:
        seen = []

    def encode_node(value):
        seen.append(("encode_node", value))
        data = value.data
        return lp.DotVertexData(
            name=data["name"],
            index=data["index"],
            payload=b"node payload",
            statements={"token": data["token"]},
        )

    def decode_node(value):
        seen.append(("decode_node", value))
        return lp.NodeValue(
            data={
                "name": value.name,
                "index": value.index,
                "payload": value.payload,
                "token": value.statements["token"],
            },
            drawing=lp.NodeDrawing(label=lp.TextLabel(value.statements["token"])),
        )

    def encode_edge(value):
        seen.append(("encode_edge", value))
        data = value.data
        return lp.DotEdgeData(
            edge_id=data["index"],
            payload=b"edge payload",
            statements={
                "codec_name": data["name"],
                "token": data["token"],
            },
            local_statements={"local": "edge"},
        )

    def decode_edge(value):
        seen.append(("decode_edge", value))
        return lp.EdgeValue(
            data={
                "name": value.statements["codec_name"],
                "index": value.edge_id,
                "payload": value.payload,
                "local": value.local_statements["local"],
                "token": value.statements["token"],
            },
            drawing=lp.EdgeDrawing(
                extensions={"codec-token": value.statements["token"]}
            ),
        )

    def encode_half_edge(value):
        seen.append(("encode_half_edge", value))
        data = value.data
        return lp.DotHalfEdgeData(
            index=data["index"],
            payload=b"half-edge payload",
            statement=data["token"],
            port_label=data["port"],
            compass=data["compass"],
        )

    def decode_half_edge(value):
        seen.append(("decode_half_edge", value))
        compass = {
            "e": lp.Compass.E,
            "w": lp.Compass.W,
        }.get(value.compass)
        return lp.HalfEdgeValue(
            data={
                "index": value.index,
                "payload": value.payload,
                "token": value.statement,
                "port": value.port_label,
                "compass": value.compass,
            },
            drawing=lp.HalfEdgeDrawing(
                statement=value.statement,
                port_label=value.port_label,
                compass=compass,
            ),
        )

    return lp.DotCodec(
        encode_node=encode_node,
        decode_node=decode_node,
        encode_edge=encode_edge,
        decode_edge=decode_edge,
        encode_half_edge=encode_half_edge,
        decode_half_edge=decode_half_edge,
    )


def codec_graph(codec):
    # Deliberately put every explicit ID in reverse encounter order. Decoding
    # must apply those IDs before invoking the callbacks.
    left = lp.node("left", data={"name": "left", "index": 1, "token": "L"})
    right = lp.node("right", data={"name": "right", "index": 0, "token": "R"})
    connection = lp.edge(
        lp.source(
            left,
            data={"index": 1, "token": "source", "port": "out", "compass": "e"},
        ),
        "connection",
        lp.sink(
            right,
            data={"index": 0, "token": "sink", "port": "in", "compass": "w"},
        ),
        data={"name": "connection", "index": 0, "token": "fermion"},
    )
    return lp.build(left, right, connection, codec=codec)


class TestGraphModel(unittest.TestCase):
    def test_declarative_build_identity_views_and_named_lookup(self):
        shared = Payload()
        left = lp.node(
            "in",
            data=shared,
            label=lp.MathSymbol("n", subscript=0),
        )
        right = lp.node("out", data=shared)
        endpoint_data = Payload()
        propagator = lp.edge(
            lp.source(left, data=endpoint_data, compass=lp.Compass.E),
            "propagator",
            lp.sink(right, data=endpoint_data, compass=lp.Compass.W),
            data=shared,
            label=lp.MathSymbol("p", subscript=0),
        )

        # Specification drawing views are live before construction.
        left.drawing.placement = lp.Placement.Start
        graph = lp.build(left, right, propagator, name="amplitude")

        self.assertEqual(graph.name, "amplitude")
        self.assertEqual((len(graph), graph.n_nodes, graph.n_edges), (2, 2, 1))
        self.assertEqual(graph.n_half_edges, 2)
        self.assertIs(graph.node("in").data, shared)
        self.assertIs(graph.node("out").data, shared)
        self.assertIs(graph.edge("propagator").data, shared)
        self.assertIs(graph.edge(0).source.data, endpoint_data)
        self.assertIs(graph.edge(0).sink.data, endpoint_data)
        self.assertEqual(graph.node("in").drawing.placement, lp.Placement.Start)
        self.assertIn("MathSymbol", repr(graph.edge(0).drawing.label))

        replacement = Payload()
        graph.node("in").data = replacement
        graph.node("in").drawing.label = lp.TextLabel("incoming")
        graph.edge("propagator").drawing.extensions = {"template-field": [1, 2]}
        self.assertIs(graph.node(0).data, replacement)
        self.assertEqual(graph.node(0).drawing.label.text, "incoming")
        self.assertEqual(graph.edge(0).drawing.extensions["template-field"], (1, 2))

        edge = graph.edge(0)
        self.assertEqual(edge.orientation, lp.Orientation.Default)
        self.assertEqual(edge.source.flow, lp.Flow.Source)
        self.assertEqual(edge.sink.flow, lp.Flow.Sink)
        self.assertEqual(edge.source.node.name, "in")
        self.assertEqual(edge.sink.node.name, "out")
        self.assertEqual(edge.source.edge.index, edge.index)
        self.assertEqual(edge.source.pair.index, edge.sink.index)
        self.assertEqual(edge.sink.pair.index, edge.source.index)
        self.assertEqual([item.index for item in graph.node("in").incidence], [0])

        with self.assertRaises(AttributeError):
            graph.node(0).index = 4
        read_only_topology = (
            (edge, "index", 4),
            (edge, "orientation", lp.Orientation.Reversed),
            (edge, "source", None),
            (edge, "sink", None),
            (edge.source, "index", 4),
            (edge.source, "flow", lp.Flow.Sink),
            (edge.source, "node", graph.node("out")),
            (edge.source, "edge", graph.edge(0)),
            (edge.source, "pair", None),
        )
        for view, field, value in read_only_topology:
            with self.subTest(field=field), self.assertRaises(AttributeError):
                setattr(view, field, value)
        with self.assertRaises(KeyError):
            graph.node("missing")
        with self.assertRaises(KeyError):
            graph.edge("missing")

    def test_external_half_edge_has_no_pair(self):
        endpoint = lp.node("endpoint")
        graph = lp.build(endpoint, lp.edge(lp.sink(endpoint), "incoming"))

        half_edge = graph.edge("incoming").sink
        self.assertIsNone(graph.edge("incoming").source)
        self.assertIsNone(half_edge.pair)
        self.assertEqual(half_edge.flow, lp.Flow.Sink)

    def test_map_callbacks_and_omitted_callbacks_preserve_identity(self):
        config = lp.RenderConfig(title=lp.TextLabel("base"))
        graph, nodes, edges, half_edges = sample_graph(render_config=config)
        seen = []

        mapped = graph.map(
            node=lambda value: (
                seen.append(("node", value.name)) or {"mapped": value.data}
            ),
            edge=lambda value: (
                seen.append(("edge", value.name)) or {"mapped": value.data}
            ),
            source=lambda value: (
                seen.append(("source", value.edge.name)) or {"mapped": value.data}
            ),
            sink=lambda value: (
                seen.append(("sink", value.edge.name)) or {"mapped": value.data}
            ),
        )

        self.assertEqual(
            seen,
            [
                ("node", "left"),
                ("node", "right"),
                ("edge", "propagator"),
                ("source", "propagator"),
                ("sink", "propagator"),
                ("edge", "outgoing"),
                ("source", "outgoing"),
            ],
        )
        self.assertIs(mapped.node("left").data["mapped"], nodes[0])
        self.assertIs(mapped.node("right").data["mapped"], nodes[1])
        self.assertIs(mapped.edge("propagator").data["mapped"], edges[0])
        self.assertIs(mapped.edge("outgoing").data["mapped"], edges[1])
        self.assertIs(mapped.edge("propagator").source.data["mapped"], half_edges[0])
        self.assertIs(mapped.edge("propagator").sink.data["mapped"], half_edges[1])
        self.assertIs(mapped.edge("outgoing").source.data["mapped"], half_edges[2])

        unchanged = graph.map()
        self.assertIs(unchanged.node("left").data, nodes[0])
        self.assertIs(unchanged.node("right").data, nodes[1])
        self.assertIs(unchanged.edge("propagator").data, edges[0])
        self.assertIs(unchanged.edge("outgoing").data, edges[1])
        self.assertIs(unchanged.edge("propagator").source.data, half_edges[0])
        self.assertIs(unchanged.edge("propagator").sink.data, half_edges[1])
        self.assertIs(unchanged.edge("outgoing").source.data, half_edges[2])

        mapped.node("left").drawing.label = lp.TextLabel("mapped")
        self.assertNotEqual(
            repr(mapped.node("left").drawing.label),
            repr(graph.node("left").drawing.label),
        )
        self.assertIsNot(mapped.render_config, graph.render_config)

        replacement = Payload()

        def mutate_data(view):
            view.data = replacement
            return view.data

        mutated = graph.map(node=mutate_data)
        self.assertTrue(all(node.data is replacement for node in graph.nodes()))
        self.assertTrue(all(node.data is replacement for node in mutated.nodes()))

        topology_mutating_graph, _, _, _ = sample_graph()

        def mutate_topology(view):
            data = view.data
            topology_mutating_graph.reorder_nodes([1, 0])
            return data

        with self.assertRaises(ReferenceError):
            topology_mutating_graph.map(node=mutate_topology)

    def test_arbitrary_data_is_never_compared_hashed_or_copied(self):
        class TrapPayload:
            def __str__(self):
                raise AssertionError("data must not be stringified")

            def __repr__(self):
                raise AssertionError("data must not be represented")

            def __hash__(self):
                raise AssertionError("data must not be hashed")

            def __eq__(self, _other):
                raise AssertionError("data must not be compared")

            def __copy__(self):
                raise AssertionError("data must not be copied")

            def __deepcopy__(self, _memo):
                raise AssertionError("data must not be deep-copied")

            def __reduce_ex__(self, _protocol):
                raise AssertionError("data must not be pickled")

        payload = TrapPayload()
        left = lp.node("left", data=payload)
        right = lp.node("right", data=payload)
        graph = lp.build(
            left,
            right,
            lp.edge(
                lp.source(left, data=payload),
                "internal",
                lp.sink(right, data=payload),
                data=payload,
            ),
            lp.edge(lp.source(right, data=payload), "external", data=payload),
        )
        mapped = graph.map()
        graph.reorder_nodes([1, 0])
        graph.reorder_edges([1, 0])

        for candidate in (
            *(node.data for node in graph.nodes()),
            *(edge.data for edge in graph.edges()),
            *(half_edge.data for half_edge in graph.half_edges()),
            *(node.data for node in mapped.nodes()),
            *(edge.data for edge in mapped.edges()),
            *(half_edge.data for half_edge in mapped.half_edges()),
        ):
            self.assertIs(candidate, payload)

    def test_reorder_preserves_associations_and_invalidates_old_views(self):
        graph, nodes, edges, half_edges = sample_graph()
        old_node = graph.node("left")
        old_edge = graph.edge("propagator")
        old_half_edge = old_edge.source
        old_node_drawing = old_node.drawing
        old_edge_drawing = old_edge.drawing

        graph.reorder_nodes([1, 0])

        for operation in (
            lambda: old_node.data,
            lambda: old_edge.data,
            lambda: old_half_edge.data,
            lambda: old_node_drawing.label,
            lambda: old_edge_drawing.label,
        ):
            with self.assertRaises(ReferenceError):
                operation()
        self.assertEqual(graph.node(0).name, "right")
        self.assertIs(graph.node("left").data, nodes[0])
        self.assertIs(graph.edge("propagator").source.node.data, nodes[0])
        self.assertIs(graph.edge("propagator").sink.node.data, nodes[1])

        edge_view = graph.edge("propagator")
        source_view = edge_view.source
        source_drawing = source_view.drawing
        graph.reverse_edge("propagator")
        for operation in (
            lambda: edge_view.source,
            lambda: source_view.flow,
            lambda: source_drawing.compass,
        ):
            with self.assertRaises(ReferenceError):
                operation()
        reversed_edge = graph.edge("propagator")
        self.assertIs(reversed_edge.source.data, half_edges[1])
        self.assertIs(reversed_edge.sink.data, half_edges[0])
        self.assertEqual(reversed_edge.source.node.name, "right")
        self.assertEqual(reversed_edge.sink.node.name, "left")

        graph.edge("outgoing").source.drawing.statement = "outgoing-source"
        graph.edge("propagator").source.drawing.statement = "propagator-source"
        graph.edge("propagator").sink.drawing.statement = "propagator-sink"
        outgoing_view = graph.edge("outgoing")
        outgoing_half_edge = outgoing_view.source
        outgoing_half_edge_drawing = outgoing_half_edge.drawing
        graph.reorder_edges([1, 0])
        for operation in (
            lambda: outgoing_view.name,
            lambda: outgoing_half_edge.data,
            lambda: outgoing_half_edge_drawing.statement,
        ):
            with self.assertRaises(ReferenceError):
                operation()
        self.assertEqual(graph.edge(0).name, "outgoing")
        self.assertIs(graph.edge("outgoing").data, edges[1])
        self.assertEqual(graph.edge(0).drawing.extensions["category"], "outgoing")
        self.assertIs(graph.half_edge(0).data, half_edges[2])
        self.assertEqual(graph.half_edge(0).drawing.statement, "outgoing-source")
        self.assertIs(graph.half_edge(1).data, half_edges[1])
        self.assertEqual(graph.half_edge(1).drawing.statement, "propagator-source")
        self.assertIs(graph.half_edge(2).data, half_edges[0])
        self.assertEqual(graph.half_edge(2).drawing.statement, "propagator-sink")

    def test_graph_specs_codecs_callbacks_and_views_are_cycle_collected(self):
        def graph_cycle():
            payload = Payload()
            spec = lp.node(data=payload)
            graph = lp.build(spec)
            payload.graph = graph
            return weakref.ref(payload)

        def spec_cycle():
            payload = Payload()
            spec = lp.node(data=payload)
            payload.spec = spec
            return weakref.ref(payload)

        def view_cycle():
            payload = Payload()
            spec = lp.node(data=payload)
            graph = lp.build(spec)
            payload.view = graph.node(0)
            return weakref.ref(payload)

        def endpoint_spec_cycle():
            payload = Payload()
            graph = lp.build(lp.node("endpoint", data=payload))
            payload.spec = lp.edge(lp.source(graph.node("endpoint")), "external")
            return weakref.ref(payload)

        def cut_partition_cycle():
            payload = Payload()
            source = lp.node("source", data=payload)
            target = lp.node("target")
            graph = lp.build(
                source,
                target,
                lp.edge(lp.source(source), "edge", lp.sink(target)),
            )
            payload.partition = graph.all_cuts(["source"], ["target"])[0]
            return weakref.ref(payload)

        class Callback:
            def __call__(self, _value):
                raise AssertionError("cycle callbacks are never executed")

        def codec_cycle():
            callback = Callback()
            codec = lp.DotCodec(
                encode_node=callback,
                decode_node=callback,
                encode_edge=callback,
                decode_edge=callback,
                encode_half_edge=callback,
                decode_half_edge=callback,
            )
            callback.codec = codec
            return weakref.ref(callback)

        def selector_cycle():
            callback = Callback()
            config = lp.RenderConfig(selectors=lp.DrawingSelectors(node=callback))
            callback.config = config
            return weakref.ref(callback)

        references = [
            graph_cycle(),
            spec_cycle(),
            view_cycle(),
            endpoint_spec_cycle(),
            cut_partition_cycle(),
            codec_cycle(),
            selector_cycle(),
        ]
        gc.collect()
        self.assertTrue(all(reference() is None for reference in references))


class TestNodeStores(unittest.TestCase):
    def test_default_and_explicit_construction(self):
        self.assertEqual(lp.Graph().node_store, lp.NodeStore.Vec)
        self.assertEqual(
            lp.Graph(node_store=lp.NodeStore.Forest).node_store,
            lp.NodeStore.Forest,
        )

        for node_store in (lp.NodeStore.Vec, lp.NodeStore.Forest):
            with self.subTest(node_store=node_store):
                graph = lp.build(lp.node("only"), node_store=node_store)
                self.assertEqual(graph.node_store, node_store)
                self.assertEqual([node.name for node in graph.nodes()], ["only"])

    def test_conversion_preserves_data_identity_and_copies_drawing(self):
        graph, node_data, edge_data, half_edge_data = sample_graph()
        original_label = repr(graph.node("left").drawing.label)
        converted = graph.to_node_store(lp.NodeStore.Forest)

        self.assertEqual(graph.node_store, lp.NodeStore.Vec)
        self.assertEqual(converted.node_store, lp.NodeStore.Forest)
        self.assertEqual(
            [node.name for node in converted.nodes()],
            [node.name for node in graph.nodes()],
        )
        self.assertIs(converted.node("left").data, node_data[0])
        self.assertIs(converted.node("right").data, node_data[1])
        self.assertIs(converted.edge("propagator").data, edge_data[0])
        self.assertIs(converted.edge("outgoing").data, edge_data[1])
        self.assertIs(
            converted.edge("propagator").source.data,
            half_edge_data[0],
        )
        self.assertIs(
            converted.edge("propagator").sink.data,
            half_edge_data[1],
        )
        self.assertIs(
            converted.edge("outgoing").source.data,
            half_edge_data[2],
        )

        converted.node("left").drawing.label = "converted node"
        converted.edge("propagator").drawing.extensions = {"category": "converted"}
        converted.edge("propagator").source.drawing.statement = "converted port"
        self.assertEqual(repr(graph.node("left").drawing.label), original_label)
        self.assertEqual(
            graph.edge("propagator").drawing.extensions["category"],
            "primary",
        )
        self.assertIs(
            graph.edge("propagator").source.drawing.statement,
            lp.INHERIT,
        )

        round_trip = converted.to_node_store(lp.NodeStore.Vec)
        self.assertEqual(round_trip.node_store, lp.NodeStore.Vec)
        self.assertIs(round_trip.node("left").data, node_data[0])
        self.assertIs(round_trip.edge("propagator").data, edge_data[0])
        self.assertIs(
            round_trip.edge("propagator").source.data,
            half_edge_data[0],
        )
        self.assertEqual(round_trip.node("left").drawing.label, "converted node")

    def test_derived_graphs_preserve_the_backend(self):
        graph, node_data, edge_data, _ = topology_graph(node_store=lp.NodeStore.Forest)

        mapped = graph.map()
        fragment = graph.concretize(graph.subgraph(edges=["ef"]))
        self.assertEqual(mapped.node_store, lp.NodeStore.Forest)
        self.assertEqual(fragment.node_store, lp.NodeStore.Forest)
        self.assertIs(mapped.node("a").data, node_data["a"])
        self.assertIs(fragment.edge("ef").data, edge_data["ef"])

        extracted_from, _, extracted_edge_data, _ = topology_graph(
            node_store=lp.NodeStore.Forest
        )
        extracted = extracted_from.extract(extracted_from.subgraph(edges=["ef"]))
        self.assertEqual(extracted_from.node_store, lp.NodeStore.Forest)
        self.assertEqual(extracted.node_store, lp.NodeStore.Forest)
        self.assertIs(extracted.edge("ef").data, extracted_edge_data["ef"])

    def test_forest_reorder_and_topology_algorithms_match_vec(self):
        vec, _, _, _ = topology_graph()
        forest, node_data, edge_data, half_edge_data = topology_graph(
            node_store=lp.NodeStore.Forest
        )
        node_order = list(reversed(range(forest.n_nodes)))
        edge_order = list(reversed(range(forest.n_edges)))
        for graph in (vec, forest):
            graph.reorder_nodes(node_order)
            graph.reorder_edges(edge_order)

        self.assertIs(forest.node("a").data, node_data["a"])
        self.assertIs(forest.edge("ef").data, edge_data["ef"])
        self.assertIs(
            forest.edge("ef").source.data,
            half_edge_data[("ef", "source")],
        )

        def topology_summary(graph):
            components = {
                frozenset(edge.name for edge in graph.edges_of(component))
                for component in graph.connected_components()
            }
            cycles, _ = graph.cycle_basis()
            cycle_edges = {
                frozenset(edge.name for edge in graph.edges_of(cycle.filter))
                for cycle in cycles
            }
            return components, cycle_edges

        self.assertEqual(topology_summary(forest), topology_summary(vec))

    def test_forest_extracts_isolated_nodes(self):
        left = lp.node("left", data=UserNodePayload("left"))
        right = lp.node("right", data=UserNodePayload("right"))
        isolated_data = UserNodePayload("isolated")
        isolated = lp.node("isolated", data=isolated_data, label="isolated label")
        graph = lp.build(
            left,
            right,
            isolated,
            lp.edge(lp.source(left), "internal", lp.sink(right)),
            node_store=lp.NodeStore.Forest,
        )
        selection = graph.subgraph(nodes=["isolated"])

        copied = graph.concretize(selection)
        self.assertEqual(copied.node_store, lp.NodeStore.Forest)
        self.assertIs(copied.node("isolated").data, isolated_data)
        copied.node("isolated").drawing.label = "copied label"
        self.assertEqual(graph.node("isolated").drawing.label, "isolated label")

        extracted = graph.extract(selection)
        self.assertEqual(extracted.node_store, lp.NodeStore.Forest)
        self.assertEqual((extracted.n_nodes, extracted.n_edges), (1, 0))
        self.assertIs(extracted.node("isolated").data, isolated_data)
        self.assertEqual(graph.node_store, lp.NodeStore.Forest)
        self.assertEqual(
            [node.name for node in graph.nodes()],
            ["left", "right"],
        )
        self.assertEqual([edge.name for edge in graph.edges()], ["internal"])

    def test_dot_import_selects_the_requested_backend(self):
        codec = lp.DotCodec.topology()
        graph, _, _, _ = sample_graph(codec=codec)
        encoded = graph.to_dot()

        default = lp.Graph.from_dot(encoded, codec)
        forest = lp.Graph.from_dot(
            encoded,
            codec,
            node_store=lp.NodeStore.Forest,
        )
        self.assertEqual(default.node_store, lp.NodeStore.Vec)
        self.assertEqual(forest.node_store, lp.NodeStore.Forest)
        self.assertEqual(
            [node.name for node in forest.nodes()],
            [node.name for node in default.nodes()],
        )
        self.assertIn("digraph", forest.to_dot())

        forest_set = lp.Graph.from_dot_set(
            "digraph one { a -> b; } digraph two { c -> d; }",
            codec,
            node_store=lp.NodeStore.Forest,
        )
        self.assertTrue(
            all(item.node_store == lp.NodeStore.Forest for item in forest_set)
        )
        with TemporaryDirectory() as directory:
            path = Path(directory) / "forest graph.dot"
            path.write_text(encoded, encoding="utf-8")
            from_file = lp.Graph.from_dot_file(
                path,
                codec,
                node_store=lp.NodeStore.Forest,
            )
        self.assertEqual(from_file.node_store, lp.NodeStore.Forest)

    def test_cross_backend_append_and_join_use_the_left_backend(self):
        backend_pairs = (
            (lp.NodeStore.Forest, lp.NodeStore.Vec),
            (lp.NodeStore.Vec, lp.NodeStore.Forest),
        )
        for left_backend, right_backend in backend_pairs:
            with self.subTest(
                left_backend=left_backend,
                right_backend=right_backend,
            ):
                left, left_node, _, left_port = dangling_graph(
                    "left",
                    lp.Flow.Source,
                    node_store=left_backend,
                )
                right, right_node, right_edge, right_port = dangling_graph(
                    "right",
                    lp.Flow.Sink,
                    node_store=right_backend,
                )

                appended = left.append(right)
                self.assertEqual(appended.node_store, left_backend)
                self.assertEqual(left.node_store, left_backend)
                self.assertEqual(right.node_store, right_backend)
                self.assertIs(appended.node("node-left").data, left_node)
                self.assertIs(appended.node("node-right").data, right_node)
                self.assertIs(appended.edge("edge-right").data, right_edge)

                joined_data = UserEdgePayload("joined")

                def merge(left_half_edge, right_half_edge):
                    self.assertIs(left_half_edge.data, left_port)
                    self.assertIs(right_half_edge.data, right_port)
                    return (
                        lp.Flow.Source,
                        lp.Orientation.Default,
                        "joined",
                        lp.EdgeValue(data=joined_data),
                    )

                joined = left.join(
                    right,
                    matching=lambda *_ports: True,
                    merge=merge,
                )
                self.assertEqual(joined.node_store, left_backend)
                self.assertIs(joined.node("node-left").data, left_node)
                self.assertIs(joined.node("node-right").data, right_node)
                self.assertIs(joined.edge("joined").data, joined_data)


class TestHedgeGraphTopology(unittest.TestCase):
    def test_incremental_node_and_edge_insertion_is_transactional(self):
        for node_store in (lp.NodeStore.Vec, lp.NodeStore.Forest):
            with self.subTest(node_store=node_store):
                graph = lp.Graph(node_store=node_store)
                a_data = UserNodePayload("a")
                b_data = UserNodePayload("b")
                a_spec = lp.node("a", data=a_data, label="node a")
                b_spec = lp.node("b", data=b_data)
                first = graph.add_node(a_spec)
                self.assertIs(first.data, a_data)

                second = graph.add_node(b_spec)
                self.assertIs(second.data, b_data)
                with self.assertRaises(ReferenceError):
                    _ = first.data

                source_data = UserHalfEdgePayload("source")
                sink_data = UserHalfEdgePayload("sink")
                edge_data = UserEdgePayload("ab")
                edge_spec = lp.edge(
                    lp.source(a_spec, data=source_data, compass=lp.Compass.E),
                    "ab",
                    lp.sink(graph.node("b"), data=sink_data),
                    data=edge_data,
                    extensions={"protocol": "sync"},
                )
                inserted = graph.add_edge(edge_spec)
                self.assertIs(inserted.data, edge_data)
                self.assertIs(inserted.source.data, source_data)
                self.assertIs(inserted.sink.data, sink_data)
                self.assertEqual(inserted.drawing.extensions["protocol"], "sync")
                self.assertEqual(inserted.source.node.name, "a")
                self.assertEqual(inserted.sink.node.name, "b")

                dangling_data = UserEdgePayload("external")
                dangling = graph.add_edge(
                    lp.edge(
                        lp.source("b"),
                        "external",
                        data=dangling_data,
                        extensions={"protocol": "async"},
                    )
                )
                self.assertIs(dangling.data, dangling_data)
                self.assertIsNotNone(dangling.source)
                self.assertIsNone(dangling.sink)

                a_spec.drawing.label = "changed spec"
                edge_spec.drawing.extensions = {"protocol": "changed spec"}
                self.assertEqual(graph.node("a").drawing.label, "node a")
                self.assertEqual(
                    graph.edge("ab").drawing.extensions["protocol"],
                    "sync",
                )

                live = graph.node("a")
                with self.assertRaisesRegex(ValueError, "duplicate node name"):
                    graph.add_node(lp.node("a"))
                self.assertIs(live.data, a_data)

                foreign = lp.Graph(node_store=node_store)
                foreign.add_node(lp.node("foreign"))
                with self.assertRaisesRegex(ValueError, "different graph revision"):
                    graph.add_edge(
                        lp.edge(
                            lp.source(foreign.node("foreign")),
                            "foreign-edge",
                        )
                    )
                self.assertIs(live.data, a_data)

                stale = graph.node("a")
                graph.add_node(lp.node("c"))
                current = graph.node("a")
                with self.assertRaises(ReferenceError):
                    graph.add_edge(lp.edge(lp.source(stale), "stale-edge"))
                self.assertIs(current.data, a_data)

    def test_build_resolves_declared_endpoint_target_forms(self):
        key_graph = lp.Graph()
        key_graph.add_node(lp.node("b"))
        live_b = key_graph.node("b")
        a = lp.node("a")
        b = lp.node("b")
        graph = lp.build(
            a,
            b,
            lp.edge(lp.source("a"), "by-key", lp.sink(1)),
            lp.edge(lp.source(live_b), "by-live-node"),
        )
        self.assertEqual(graph.edge("by-key").source.node.name, "a")
        self.assertEqual(graph.edge("by-key").sink.node.name, "b")
        self.assertEqual(graph.edge("by-live-node").source.node.name, "b")

    def test_incremental_split_and_connect_preserve_endpoint_values(self):
        for node_store in (lp.NodeStore.Vec, lp.NodeStore.Forest):
            with self.subTest(node_store=node_store):
                graph = lp.Graph(node_store=node_store)
                graph.add_node(lp.node("a"))
                graph.add_node(lp.node("b"))
                source_data = UserHalfEdgePayload("source")
                sink_data = UserHalfEdgePayload("sink")
                original_data = UserEdgePayload("original")
                graph.add_edge(
                    lp.edge(
                        lp.source("a", data=source_data),
                        "ab",
                        lp.sink("b", data=sink_data),
                        data=original_data,
                    )
                )

                original = graph.edge("ab")
                selected, opposite = graph.split_edge(
                    original.source,
                    lp.EdgeValue(
                        data=UserEdgePayload("split"),
                        drawing=lp.EdgeDrawing(extensions={"operation": "split"}),
                    ),
                    name="a-external",
                )
                self.assertEqual((graph.n_edges, graph.n_half_edges), (2, 2))
                self.assertEqual(selected.name, "a-external")
                self.assertEqual(opposite.name, "ab")
                self.assertEqual(selected.drawing.extensions["operation"], "split")
                self.assertIs(selected.source.data, source_data)
                self.assertIs(opposite.sink.data, sink_data)
                self.assertIs(opposite.data, original_data)
                with self.assertRaises(ReferenceError):
                    _ = original.data

                joined_data = UserEdgePayload("joined")
                selected_view = selected
                joined = graph.connect(
                    selected.source,
                    opposite.sink,
                    lp.EdgeValue(
                        data=joined_data,
                        drawing=lp.EdgeDrawing(extensions={"operation": "joined"}),
                    ),
                    name="joined",
                )
                self.assertEqual((graph.n_edges, graph.n_half_edges), (1, 2))
                self.assertIs(joined.data, joined_data)
                self.assertEqual(joined.drawing.extensions["operation"], "joined")
                self.assertEqual(joined.source.node.name, "a")
                self.assertEqual(joined.sink.node.name, "b")
                self.assertIs(joined.source.data, source_data)
                self.assertIs(joined.sink.data, sink_data)
                with self.assertRaises(ReferenceError):
                    _ = selected_view.data

                live = graph.edge("joined")
                with self.assertRaisesRegex(ValueError, "two dangling"):
                    graph.connect(
                        live.source,
                        live.sink,
                        lp.EdgeValue(data=UserEdgePayload("invalid")),
                    )
                self.assertIs(live.data, joined_data)

    def test_subgraph_views_and_set_algebra_use_half_edge_indices(self):
        graph, _, _, _ = topology_graph()
        full = graph.full_subgraph()
        empty = graph.empty_subgraph()
        ab = graph.subgraph(edges=["ab"])
        bc = graph.subgraph(edges=["bc"])

        self.assertIsInstance(full, lp.Subgraph)
        self.assertIs(full.graph, graph)
        self.assertEqual(
            (len(full), full.size),
            (graph.n_half_edges, graph.n_half_edges),
        )
        self.assertFalse(empty)
        self.assertEqual(len(ab), 2)
        self.assertEqual(
            set(ab.half_edge_indices()),
            {graph.edge("ab").source.index, graph.edge("ab").sink.index},
        )
        self.assertEqual([edge.name for edge in graph.edges_of(ab)], ["ab"])
        self.assertEqual({node.name for node in graph.nodes_of(ab)}, {"a", "b"})
        self.assertTrue(ab.is_disjoint(bc))
        self.assertEqual(len(ab | bc), 4)
        self.assertEqual(len((ab | bc) - bc), 2)
        self.assertEqual(len(ab & bc), 0)
        self.assertEqual(len(ab ^ bc), 4)
        self.assertEqual(len(~ab), graph.n_half_edges - 2)
        self.assertTrue(ab < full)
        self.assertTrue(full >= ab)
        self.assertIn(graph.edge("ab").source.index, ab)

        other, _, _, _ = topology_graph()
        with self.assertRaisesRegex(ValueError, "different graph revision"):
            ab.union(other.full_subgraph())

    def test_construct_and_filter_union_live_views(self):
        graph, _, _, _ = topology_graph()
        at_b = graph.subgraph(nodes=["b"])
        self.assertEqual(len(at_b), 2)
        self.assertEqual(
            {edge.name for edge in graph.edges_of(at_b)},
            {"ab", "bc"},
        )
        self.assertEqual(
            {half_edge.node.name for half_edge in at_b.to_half_edges()},
            {"b"},
        )

        seen = {"node": [], "edge": [], "half_edge": []}
        selected = graph.filter(
            node=lambda node: seen["node"].append(node.name) or node.name == "a",
            edge=lambda edge: seen["edge"].append(edge.name) or edge.name == "ef",
            half_edge=lambda half_edge: (
                seen["half_edge"].append(half_edge.index) or half_edge.node.name == "d"
            ),
        )
        self.assertEqual(len(seen["node"]), graph.n_nodes)
        self.assertEqual(len(seen["edge"]), graph.n_edges)
        self.assertEqual(len(seen["half_edge"]), graph.n_half_edges)
        self.assertEqual(len(selected), 5)
        self.assertEqual(
            {edge.name for edge in graph.edges_of(selected)},
            {"ab", "ca", "cd", "ef"},
        )

        triangle = graph.subgraph(edges=["ab", "bc", "ca"])
        self.assertEqual(graph.count_connected_components(), 2)
        self.assertEqual(graph.count_connected_components(triangle), 1)
        self.assertEqual(graph.cyclomatic_number(), 1)
        self.assertEqual(graph.cyclomatic_number(triangle), 1)

        class PredicateFailure(RuntimeError):
            pass

        def fail(_edge):
            raise PredicateFailure("filter failed")

        with self.assertRaisesRegex(PredicateFailure, "filter failed"):
            graph.filter(edge=fail)

        mutated, _, _, _ = topology_graph()

        def mutate(_node):
            mutated.reverse_edge("ab")
            return True

        with self.assertRaises(ReferenceError):
            mutated.filter(node=mutate)

    def test_components_traversals_cycles_and_bridges(self):
        graph, _, _, _ = topology_graph()
        components = graph.connected_components()
        component_edges = {
            frozenset(edge.name for edge in graph.edges_of(component))
            for component in components
        }
        self.assertEqual(
            component_edges,
            {frozenset({"ab", "bc", "ca", "cd"}), frozenset({"ef"})},
        )
        self.assertFalse(graph.is_connected())
        self.assertTrue(graph.is_connected(graph.subgraph(edges=["ab", "bc", "ca"])))

        bridges = graph.bridges()
        self.assertEqual(
            {edge.name for edge in graph.edges_of(bridges)},
            {"cd", "ef"},
        )

        cycles, covered = graph.cycle_basis()
        self.assertEqual(len(cycles), 1)
        self.assertIsInstance(cycles[0], lp.Cycle)
        self.assertEqual(
            {edge.name for edge in graph.edges_of(cycles[0].filter)},
            {"ab", "bc", "ca"},
        )
        self.assertGreaterEqual(len(covered), len(cycles[0].filter))

        triangle = graph.subgraph(edges=["ab", "bc", "ca"])
        forests = graph.all_spanning_forests(triangle)
        self.assertEqual(len(forests), 3)
        self.assertTrue(all(len(forest) == 4 for forest in forests))

        cuts = graph.all_cuts(["a"], ["d"])
        self.assertTrue(cuts)
        self.assertTrue(all(isinstance(cut, lp.CutPartition) for cut in cuts))
        ordering = []
        for partition in cuts:
            cut = partition.boundary
            self.assertIsInstance(cut, lp.OrientedCut)
            self.assertEqual(cut.left, cut.side(True))
            self.assertEqual(cut.right, cut.side(False))
            self.assertTrue(partition.source_side.is_disjoint(partition.target_side))
            self.assertEqual(
                {edge.index for edge in partition.edges},
                {edge.index for edge in graph.edges_of(cut.left)},
            )
            ordering.append(
                (
                    partition.source_side.half_edge_indices(),
                    cut.left.half_edge_indices(),
                    partition.target_side.half_edge_indices(),
                )
            )
        self.assertEqual(ordering, sorted(ordering))

        with self.assertRaisesRegex(ValueError, "must be disjoint"):
            graph.all_cuts(["a"], ["a"])
        with self.assertRaisesRegex(ValueError, "boundary half-edge"):
            graph.all_cuts(["a", "b", "c", "d"], ["e"])

        isolated = lp.build(lp.node("isolated-source"), lp.node("isolated-target"))
        with self.assertRaisesRegex(ValueError, "incident half-edge"):
            isolated.all_cuts(["isolated-source"], ["isolated-target"])

        mixed, _, _, _ = topology_graph()
        mixed.add_node(lp.node("isolated"))
        with self.assertRaisesRegex(ValueError, "every source and target node"):
            mixed.all_cuts(["a", "isolated"], ["d"])

        for traverse in (graph.depth_first_traverse, graph.breadth_first_traverse):
            tree = traverse("a")
            self.assertIsInstance(tree, lp.TraversalTree)
            self.assertEqual(tree.nodes[0].name, "a")
            self.assertEqual(
                {node.name for node in tree.nodes},
                {"a", "b", "c", "d"},
            )
            self.assertEqual(len(tree.half_edges), 6)

    def test_cut_partitions_match_native_square_regression(self):
        for node_store in (lp.NodeStore.Vec, lp.NodeStore.Forest):
            with self.subTest(node_store=node_store):
                nodes = {name: lp.node(name) for name in ("a", "b", "c", "d")}

                def connection(source, sink, name):
                    return lp.edge(
                        lp.source(nodes[source]),
                        name,
                        lp.sink(nodes[sink]),
                    )

                graph = lp.build(
                    *nodes.values(),
                    connection("a", "b", "ab"),
                    connection("b", "c", "bc"),
                    connection("c", "d", "cd"),
                    connection("d", "a", "da"),
                    connection("b", "d", "bd"),
                    node_store=node_store,
                )
                partitions = graph.all_cuts(["a"], ["c"])
                self.assertEqual(len(partitions), 4)
                for partition in partitions:
                    self.assertTrue(
                        partition.source_side.includes_node(graph.node("a").index)
                    )
                    self.assertTrue(
                        partition.target_side.includes_node(graph.node("c").index)
                    )
                    self.assertEqual(
                        partition.source_side | partition.target_side,
                        graph.full_subgraph(),
                    )
                    right = set(partition.boundary.right.half_edge_indices())
                    for half_edge in partition.boundary.left.to_half_edges():
                        self.assertIn(half_edge.pair.index, right)

                partition = partitions[0]
                boundary = partition.boundary
                graph.reverse_edge("ab")
                with self.assertRaises(ReferenceError):
                    _ = partition.source_side
                with self.assertRaises(ReferenceError):
                    _ = boundary.left

    def test_bond_enumeration_uses_native_size_bounds(self):
        for node_store in (lp.NodeStore.Vec, lp.NodeStore.Forest):
            with self.subTest(node_store=node_store):
                names = ("0", "1", "2", "3", "a", "b", "c")
                nodes = {name: lp.node(name) for name in names}

                def connection(source, sink, index):
                    return lp.edge(
                        lp.source(nodes[source]),
                        f"e{index}",
                        lp.sink(nodes[sink]),
                    )

                graph = lp.build(
                    *nodes.values(),
                    connection("0", "1", 0),
                    connection("1", "2", 1),
                    connection("0", "3", 2),
                    connection("2", "a", 3),
                    connection("a", "b", 4),
                    connection("b", "c", 5),
                    connection("c", "3", 6),
                    connection("3", "0", 7),
                    connection("2", "1", 8),
                    node_store=node_store,
                )
                bonds = graph.all_bonds(min_size=2, max_size=2)
                self.assertEqual(len(bonds), 10)
                self.assertTrue(all(bond.n_half_edges == 2 for bond in bonds))
                self.assertEqual(
                    [bond.half_edge_indices() for bond in bonds],
                    sorted(bond.half_edge_indices() for bond in bonds),
                )
                with self.assertRaisesRegex(ValueError, "min_size"):
                    graph.all_bonds(min_size=3, max_size=2)

    def test_bond_defaults_restriction_and_live_revision(self):
        for node_store in (lp.NodeStore.Vec, lp.NodeStore.Forest):
            with self.subTest(node_store=node_store):
                nodes = {name: lp.node(name) for name in ("a", "b", "c", "d")}
                graph = lp.build(
                    *nodes.values(),
                    lp.edge(lp.source(nodes["a"]), "ab", lp.sink(nodes["b"])),
                    lp.edge(lp.source(nodes["c"]), "cd", lp.sink(nodes["d"])),
                    node_store=node_store,
                )
                bonds = graph.all_bonds()
                self.assertEqual(len(bonds), 2)
                self.assertTrue(all(bond.n_half_edges > 0 for bond in bonds))

                restricted = graph.all_bonds(subgraph=graph.subgraph(edges=["ab"]))
                self.assertEqual(len(restricted), 1)
                self.assertTrue(restricted[0].is_subset(graph.subgraph(edges=["ab"])))

                stale = restricted[0]
                graph.add_node(lp.node("isolated"))
                with self.assertRaises(ReferenceError):
                    _ = stale.n_half_edges

                edgeless = lp.build(lp.node("only"), node_store=node_store)
                self.assertEqual(edgeless.all_bonds(), [])
                self.assertEqual(edgeless.all_bonds(min_size=2), [])
                with self.assertRaisesRegex(ValueError, "at least 1"):
                    edgeless.all_bonds(min_size=0)

    def test_boundaries_cycles_bonds_and_tadpoles_are_backend_neutral(self):
        for node_store in (lp.NodeStore.Vec, lp.NodeStore.Forest):
            with self.subTest(node_store=node_store):
                graph, _, _, _ = topology_graph(node_store=node_store)
                triangle = graph.subgraph(edges=["ab", "bc", "ca"])
                self.assertFalse(graph.internal_boundary(triangle))
                self.assertEqual(
                    {edge.name for edge in graph.edges_of(graph.boundary(triangle))},
                    {"cd"},
                )

                split = graph.subgraph(half_edges=[graph.edge("cd").source.index])
                self.assertEqual(graph.internal_boundary(split), split)

                cycles = graph.all_cycles(max_results=1)
                self.assertEqual(len(cycles), 1)
                self.assertEqual(
                    {edge.name for edge in graph.edges_of(cycles[0].filter)},
                    {"ab", "bc", "ca"},
                )
                with self.assertRaisesRegex(
                    ValueError, "1 candidate combinations.*max_results=0"
                ):
                    graph.all_cycles(max_results=0)
                self.assertEqual(
                    graph.all_cycles(
                        max_results=0,
                        subgraph=graph.subgraph(edges=["cd"]),
                    ),
                    [],
                )
                loop_node = lp.node("loop")
                large_cycle_space = lp.build(
                    loop_node,
                    *[
                        lp.edge(
                            lp.source(loop_node),
                            f"loop-{index}",
                            lp.sink(loop_node),
                        )
                        for index in range(64)
                    ],
                    node_store=node_store,
                )
                with self.assertRaisesRegex(
                    ValueError, "powerset at rank 64 is not representable"
                ):
                    large_cycle_space.all_cycles(max_results=2**64 - 1)

                bond = graph.find_bond(min_size=1, max_size=1)
                self.assertIsInstance(bond, lp.Subgraph)
                self.assertEqual(bond.n_half_edges, 1)
                self.assertIsNone(graph.find_bond(min_size=10, max_size=10))

                terminal_names = ("terminal", "a", "b", "c", "d", "e", "f", "x")
                nodes = {name: lp.node(name) for name in terminal_names}

                def connection(source, sink, name):
                    return lp.edge(
                        lp.source(nodes[source]),
                        name,
                        lp.sink(nodes[sink]),
                    )

                tadpole_graph = lp.build(
                    *nodes.values(),
                    connection("terminal", "a", "ta"),
                    connection("a", "b", "ab"),
                    connection("b", "c", "bc"),
                    connection("c", "a", "ca"),
                    connection("d", "e", "de"),
                    connection("e", "f", "ef"),
                    connection("f", "d", "fd"),
                    node_store=node_store,
                )
                tadpoles = tadpole_graph.tadpoles([tadpole_graph.node("terminal")])
                self.assertEqual(len(tadpoles), 2)
                self.assertEqual(
                    [part.half_edge_indices() for part in tadpoles],
                    sorted(part.half_edge_indices() for part in tadpoles),
                )
                with self.assertRaisesRegex(ValueError, "at least one"):
                    tadpole_graph.tadpoles([])
                with self.assertRaisesRegex(ValueError, "more than once"):
                    tadpole_graph.tadpoles(["terminal", "terminal"])
                with self.assertRaisesRegex(ValueError, "no incident"):
                    tadpole_graph.tadpoles(["x"])

                sample, _, _, _ = sample_graph(node_store=node_store)
                external = sample.external_half_edges()
                self.assertIsInstance(external, lp.Subgraph)
                self.assertEqual(
                    [half_edge.index for half_edge in external.to_half_edges()],
                    [sample.edge("outgoing").source.index],
                )

    def test_traversal_relationships_and_cut_orientation_are_graph_bound(self):
        for node_store in (lp.NodeStore.Vec, lp.NodeStore.Forest):
            with self.subTest(node_store=node_store):
                graph, _, _, _ = topology_graph(node_store=node_store)
                graph.add_edge(lp.edge(lp.source("d"), "external"))
                tree = graph.depth_first_traverse("a")
                root = tree.nodes[0]
                self.assertEqual(root.name, "a")
                self.assertIsNone(tree.parent(root))
                self.assertEqual(tree.ancestors(root), [])

                for node in tree.nodes[1:]:
                    parent = tree.parent(node)
                    self.assertIsNotNone(parent)
                    self.assertIn(
                        node.name,
                        {child.name for child in tree.children(parent)},
                    )
                    self.assertEqual(tree.ancestors(node)[-1].name, "a")

                with self.assertRaisesRegex(ValueError, "not part"):
                    tree.parent(graph.node("e"))
                foreign, _, _, _ = topology_graph(node_store=node_store)
                with self.assertRaisesRegex(ValueError, "different graph revision"):
                    tree.children(foreign.node("a"))

                tree_indices = {half_edge.index for half_edge in tree.half_edges}
                triangle_edges = [graph.edge(name) for name in ("ab", "bc", "ca")]
                non_tree = next(
                    edge
                    for edge in triangle_edges
                    if edge.source.index not in tree_indices
                    and edge.sink.index not in tree_indices
                )
                tree_edge = next(
                    edge for edge in triangle_edges if edge.source.index in tree_indices
                )
                cycle = tree.fundamental_cycle(non_tree.source)
                self.assertIsInstance(cycle, lp.Cycle)
                self.assertIsNone(tree.fundamental_cycle(tree_edge.source))
                with self.assertRaisesRegex(ValueError, "external half-edge"):
                    tree.fundamental_cycle(graph.edge("external").source)
                with self.assertRaisesRegex(ValueError, "both edge endpoints"):
                    tree.fundamental_cycle(graph.edge("ef").source)

                partition = graph.all_cuts(["a"], ["d"])[0]
                cut = partition.boundary
                self.assertEqual(
                    [edge.index for edge in cut.edges],
                    sorted(edge.index for edge in cut.edges),
                )
                for edge in cut.edges:
                    self.assertNotEqual(
                        cut.orientation(edge), lp.Orientation.Undirected
                    )
                non_cut = next(
                    edge
                    for edge in graph.edges()
                    if edge.index not in {cut_edge.index for cut_edge in cut.edges}
                )
                with self.assertRaisesRegex(ValueError, "not part of this cut"):
                    cut.orientation(non_cut)
                self.assertIsInstance(cut.winding_number(cycle), int)
                foreign_cycle = foreign.all_cycles(max_results=1)[0]
                with self.assertRaisesRegex(ValueError, "different graph revision"):
                    cut.winding_number(foreign_cycle)

                fresh = graph.add_node(lp.node("fresh"))
                for relationship in (tree.parent, tree.children, tree.ancestors):
                    with self.assertRaises(ReferenceError):
                        relationship(fresh)
                    with self.assertRaises(ReferenceError):
                        relationship("fresh")
                with self.assertRaises(ReferenceError):
                    cut.orientation("no-longer-resolvable")
                with self.assertRaises(ReferenceError):
                    cut.winding_number(cycle)

    def test_directed_algorithms_select_underlying_or_superficial_direction(self):
        for node_store in (lp.NodeStore.Vec, lp.NodeStore.Forest):
            with self.subTest(node_store=node_store):
                nodes = {
                    name: lp.node(name, data=UserNodePayload(name)) for name in "abcx"
                }
                edge_data = {name: UserEdgePayload(name) for name in ("ab", "bc", "ac")}

                def connection(source, sink, name, orientation=lp.Orientation.Default):
                    return lp.edge(
                        lp.source(nodes[source]),
                        name,
                        lp.sink(nodes[sink]),
                        data=edge_data[name],
                        orientation=orientation,
                        label=f"label-{name}",
                    )

                graph = lp.build(
                    *nodes.values(),
                    connection("a", "b", "ab"),
                    connection("b", "c", "bc"),
                    connection("a", "c", "ac", lp.Orientation.Reversed),
                    node_store=node_store,
                )
                source_view = graph.edge("ab")
                self.assertFalse(graph.is_reachable("c", "a"))
                self.assertTrue(
                    graph.is_reachable(
                        "c",
                        "a",
                        direction=lp.DirectionBasis.Superficial,
                    )
                )
                self.assertEqual(
                    [node.name for node in graph.topological_order()],
                    ["x", "a", "b", "c"],
                )
                with self.assertRaisesRegex(ValueError, "Not a DAG"):
                    graph.topological_order(direction=lp.DirectionBasis.Superficial)

                selected = graph.subgraph(edges=["ab"])
                self.assertTrue(graph.is_reachable("a", "b", subgraph=selected))
                with self.assertRaisesRegex(ValueError, "target node"):
                    graph.is_reachable("a", "c", subgraph=selected)
                isolated = graph.subgraph(nodes=["x"])
                self.assertTrue(graph.is_reachable("x", "x", subgraph=isolated))
                self.assertEqual(
                    [node.name for node in graph.topological_order(subgraph=isolated)],
                    ["x"],
                )

                reduced = graph.transitive_reduction()
                self.assertEqual({edge.name for edge in reduced.edges()}, {"ab", "bc"})
                self.assertEqual(graph.n_edges, 3)
                self.assertIs(source_view.data, edge_data["ab"])
                self.assertIs(reduced.edge("ab").data, edge_data["ab"])
                reduced.edge("ab").drawing.label = "changed"
                self.assertEqual(graph.edge("ab").drawing.label, "label-ab")
                with self.assertRaisesRegex(ValueError, "Not a DAG"):
                    graph.transitive_reduction(direction=lp.DirectionBasis.Superficial)

                left = lp.node("left")
                right = lp.node("right")
                reversed_graph = lp.build(
                    left,
                    right,
                    lp.edge(
                        lp.source(left),
                        "edge",
                        lp.sink(right),
                        orientation=lp.Orientation.Reversed,
                    ),
                    node_store=node_store,
                )
                self.assertEqual(
                    [node.name for node in reversed_graph.topological_order()],
                    ["left", "right"],
                )
                self.assertEqual(
                    [
                        node.name
                        for node in reversed_graph.topological_order(
                            direction=lp.DirectionBasis.Superficial
                        )
                    ],
                    ["right", "left"],
                )

    def test_isolated_nodes_are_first_class_subgraph_members(self):
        a = lp.node("a", data=UserNodePayload("a"))
        b = lp.node("b", data=UserNodePayload("b"))
        x = lp.node("x", data=UserNodePayload("x"))
        y = lp.node("y", data=UserNodePayload("y"))
        graph = lp.build(
            a,
            b,
            x,
            y,
            lp.edge(lp.source(a), "ab", lp.sink(b)),
        )

        selected_x = graph.subgraph(nodes=["x"])
        selected_y = graph.filter(node=lambda node: node.name == "y")
        self.assertTrue(selected_x)
        self.assertEqual(
            (len(selected_x), selected_x.n_half_edges, selected_x.n_isolated_nodes),
            (1, 0, 1),
        )
        self.assertEqual(selected_x.isolated_node_indices(), [graph.node("x").index])
        self.assertTrue(selected_x.includes_node(graph.node("x").index))
        self.assertEqual([node.name for node in graph.nodes_of(selected_x)], ["x"])
        self.assertTrue(selected_x.is_disjoint(selected_y))
        self.assertEqual(len(selected_x | selected_y), 2)
        self.assertEqual(len(~selected_x), graph.n_half_edges + 1)

        full = graph.full_subgraph()
        self.assertEqual(
            (len(full), full.size, full.n_half_edges, full.n_isolated_nodes),
            (graph.n_half_edges + 2, graph.n_half_edges + 2, graph.n_half_edges, 2),
        )
        self.assertEqual(graph.count_connected_components(), 3)
        self.assertEqual(len(graph.connected_components()), 3)
        self.assertFalse(graph.is_connected())
        self.assertEqual(graph.count_connected_components(selected_x), 1)
        self.assertTrue(graph.is_connected(selected_x))
        self.assertEqual(graph.cyclomatic_number(selected_x), 0)

        for traverse in (graph.depth_first_traverse, graph.breadth_first_traverse):
            tree = traverse("x", subgraph=selected_x)
            self.assertEqual([node.name for node in tree.nodes], ["x"])
            self.assertEqual(tree.half_edges, [])
            self.assertEqual(tree.subgraph, selected_x)

        isolated_forests = graph.all_spanning_forests(selected_x)
        self.assertEqual(len(isolated_forests), 1)
        self.assertEqual(isolated_forests[0], selected_x)
        self.assertTrue(
            all(forest.n_isolated_nodes == 2 for forest in graph.all_spanning_forests())
        )

    def test_isolated_node_owning_transformations(self):
        a = lp.node("a", data=UserNodePayload("a"))
        b = lp.node("b", data=UserNodePayload("b"))
        x_data = UserNodePayload("x")
        y_data = UserNodePayload("y")
        x = lp.node("x", data=x_data, label="x-label")
        y = lp.node("y", data=y_data)
        graph = lp.build(
            a,
            b,
            x,
            y,
            lp.edge(lp.source(a), "ab", lp.sink(b)),
        )
        selected_x = graph.subgraph(nodes=["x"])

        copied = graph.concretize(selected_x)
        self.assertEqual((copied.n_nodes, copied.n_edges), (1, 0))
        self.assertIs(copied.node("x").data, x_data)
        copied.node("x").drawing.label = "changed"
        self.assertEqual(graph.node("x").drawing.label, "x-label")

        old_x = graph.node("x")
        extracted = graph.extract(selected_x)
        self.assertEqual((extracted.n_nodes, extracted.n_edges), (1, 0))
        self.assertIs(extracted.node("x").data, x_data)
        self.assertEqual((graph.n_nodes, graph.n_edges), (3, 1))
        self.assertIs(graph.node("y").data, y_data)
        with self.assertRaises(ReferenceError):
            _ = old_x.data

        selected_y = graph.subgraph(nodes=["y"])
        graph.delete(selected_y)
        self.assertEqual((graph.n_nodes, graph.n_edges), (2, 1))
        with self.assertRaises(KeyError):
            graph.node("y")

        left = lp.node("left")
        right = lp.node("right")
        isolated = lp.node("isolated")
        survivor = lp.node("survivor")
        mixed = lp.build(
            left,
            right,
            isolated,
            survivor,
            lp.edge(lp.source(left), "internal", lp.sink(right)),
        )
        replacement_data = UserNodePayload("replacement")
        mixed.contract(
            mixed.subgraph(edges=["internal"], nodes=["isolated"]),
            lp.NodeValue(data=replacement_data),
            name="contracted",
        )
        self.assertEqual((mixed.n_nodes, mixed.n_edges), (2, 0))
        self.assertIs(mixed.node("contracted").data, replacement_data)
        self.assertEqual(mixed.node("survivor").name, "survivor")

        isolate_only = lp.build(lp.node("u"), lp.node("v"))
        isolate_only.contract(
            isolate_only.full_subgraph(),
            lp.NodeValue(data=replacement_data),
            name="only",
        )
        self.assertEqual((isolate_only.n_nodes, isolate_only.n_edges), (1, 0))
        self.assertIs(isolate_only.node("only").data, replacement_data)

    def test_concretize_preserves_data_identity_and_copies_drawing(self):
        graph, node_data, edge_data, half_edge_data = topology_graph()
        selected = graph.subgraph(edges=["ef"])
        fragment = graph.concretize(selected)

        self.assertEqual((fragment.n_nodes, fragment.n_edges), (2, 1))
        self.assertEqual((graph.n_nodes, graph.n_edges), (6, 5))
        self.assertIs(fragment.node("e").data, node_data["e"])
        self.assertIs(fragment.node("f").data, node_data["f"])
        self.assertIs(fragment.edge("ef").data, edge_data["ef"])
        self.assertIs(
            fragment.edge("ef").source.data,
            half_edge_data[("ef", "source")],
        )
        self.assertIs(
            fragment.edge("ef").sink.data,
            half_edge_data[("ef", "sink")],
        )
        fragment.edge("ef").drawing.extensions = {"channel": "changed"}
        self.assertEqual(graph.edge("ef").drawing.extensions["channel"], "channel-ef")

    def test_extract_and_delete_are_owning_topology_mutations(self):
        graph, node_data, edge_data, _ = topology_graph()
        selected = graph.subgraph(edges=["ef"])
        old_node = graph.node("e")
        old_edge = graph.edge("ef")
        fragment = graph.extract(selected)

        self.assertEqual((graph.n_nodes, graph.n_edges), (4, 4))
        self.assertEqual((fragment.n_nodes, fragment.n_edges), (2, 1))
        self.assertIs(fragment.node("e").data, node_data["e"])
        self.assertIs(fragment.edge("ef").data, edge_data["ef"])
        with self.assertRaises(KeyError):
            graph.edge("ef")
        for stale in (
            lambda: old_node.data,
            lambda: old_edge.data,
            lambda: selected.half_edge_indices(),
        ):
            with self.assertRaises(ReferenceError):
                stale()

        delete_graph, _, _, _ = topology_graph()
        delete_selection = delete_graph.subgraph(edges=["ef"])
        delete_graph.delete(delete_selection)
        self.assertEqual((delete_graph.n_nodes, delete_graph.n_edges), (4, 4))
        with self.assertRaises(KeyError):
            delete_graph.node("e")

    def test_contract_closed_component_restores_replacement_and_isolated_nodes(self):
        for node_store in (lp.NodeStore.Vec, lp.NodeStore.Forest):
            for with_isolated in (False, True):
                with self.subTest(
                    node_store=node_store,
                    with_isolated=with_isolated,
                ):
                    self.check_contract_closed_component(
                        node_store,
                        with_isolated,
                    )

    def check_contract_closed_component(self, node_store, with_isolated):
        left_data = UserNodePayload("left")
        right_data = UserNodePayload("right")
        edge_data = UserEdgePayload("internal")
        left = lp.node("left", data=left_data)
        right = lp.node("right", data=right_data)
        items = [
            left,
            right,
            lp.edge(
                lp.source(left),
                "internal",
                lp.sink(right),
                data=edge_data,
            ),
        ]
        isolated_data = UserNodePayload("pre-existing-isolated")
        if with_isolated:
            items.append(lp.node("isolated", data=isolated_data))
        graph = lp.build(*items, node_store=node_store)
        selected = graph.subgraph(edges=["internal"])
        replacement_data = UserNodePayload("replacement")

        graph.contract(
            selected,
            lp.NodeValue(data=replacement_data),
            name="contracted",
        )

        self.assertEqual(graph.n_edges, 0)
        self.assertEqual(graph.n_nodes, 2 if with_isolated else 1)
        self.assertIs(graph.node("contracted").data, replacement_data)
        if with_isolated:
            self.assertIs(graph.node("isolated").data, isolated_data)
        with self.assertRaises(ReferenceError):
            selected.half_edge_indices()

    def test_append_and_append_mut_preserve_identity(self):
        left, left_node, left_edge, left_half_edge = dangling_graph(
            "left", lp.Flow.Source
        )
        right, right_node, right_edge, right_half_edge = dangling_graph(
            "right", lp.Flow.Sink
        )
        appended = left.append(right)

        self.assertEqual((left.n_nodes, left.n_edges), (1, 1))
        self.assertEqual((right.n_nodes, right.n_edges), (1, 1))
        self.assertEqual((appended.n_nodes, appended.n_edges), (2, 2))
        self.assertIs(appended.node("node-left").data, left_node)
        self.assertIs(appended.node("node-right").data, right_node)
        self.assertIs(appended.edge("edge-left").data, left_edge)
        self.assertIs(appended.edge("edge-right").data, right_edge)
        self.assertIs(appended.edge("edge-left").source.data, left_half_edge)
        self.assertIs(appended.edge("edge-right").sink.data, right_half_edge)

        old_view = left.node("node-left")
        left.append_mut(right)
        self.assertEqual((left.n_nodes, left.n_edges), (2, 2))
        with self.assertRaises(ReferenceError):
            _ = old_view.data

    def test_join_and_callback_errors_are_transactional(self):
        left, left_node, _, left_half_edge = dangling_graph("left", lp.Flow.Source)
        right, right_node, _, right_half_edge = dangling_graph("right", lp.Flow.Sink)
        merged_data = UserEdgePayload("merged-coupling")
        matched = []
        merged = []

        def matching(left_port, right_port):
            matched.append((left_port.data, right_port.data))
            return True

        def merge(left_port, right_port):
            merged.append((left_port.data, right_port.data))
            return (
                lp.Flow.Source,
                lp.Orientation.Default,
                "joined",
                lp.EdgeValue(data=merged_data),
            )

        joined = left.join(right, matching=matching, merge=merge)
        self.assertEqual((joined.n_nodes, joined.n_edges), (2, 1))
        self.assertEqual(matched, [(left_half_edge, right_half_edge)])
        self.assertEqual(merged, [(left_half_edge, right_half_edge)])
        self.assertIs(joined.edge("joined").data, merged_data)
        self.assertIs(joined.edge("joined").source.node.data, left_node)
        self.assertIs(joined.edge("joined").sink.node.data, right_node)
        self.assertEqual((left.n_edges, right.n_edges), (1, 1))

        class CallbackFailure(RuntimeError):
            pass

        old_left = left.edge("edge-left")

        def fail_matching(_left_port, _right_port):
            raise CallbackFailure("matcher failed")

        with self.assertRaisesRegex(CallbackFailure, "matcher failed"):
            left.join_mut(
                right,
                matching=fail_matching,
                merge=lambda *_ports: self.fail("merge must not run"),
            )

        def fail_merge(_left_port, _right_port):
            raise CallbackFailure("merge failed")

        with self.assertRaisesRegex(CallbackFailure, "merge failed"):
            left.join_mut(right, matching=lambda *_ports: True, merge=fail_merge)
        self.assertEqual(
            (left.n_nodes, left.n_edges, left.n_half_edges),
            (1, 1, 1),
        )
        self.assertIs(old_left.data, left.edge("edge-left").data)


class TestDotCodec(unittest.TestCase):
    def test_published_python_quickstart_runs(self):
        source_path = (
            Path(__file__).resolve().parents[3]
            / "docs/products/linnet/content/quickstart-python.typ"
        )
        source = source_path.read_text(encoding="utf-8")
        example = re.search(
            r"// docs-example: compile linnet-python-quickstart\s*"
            r"```python\n(.*?)\n```",
            source,
            re.DOTALL,
        )
        self.assertIsNotNone(example)
        exec(
            compile(example.group(1), str(source_path), "exec"),
            {"print": lambda *_args, **_kwargs: None},
        )

    def test_constructor_requires_six_callable_callbacks(self):
        callbacks = {
            "encode_node": lambda _value: None,
            "decode_node": lambda _value: None,
            "encode_edge": lambda _value: None,
            "decode_edge": lambda _value: None,
            "encode_half_edge": lambda _value: None,
            "decode_half_edge": lambda _value: None,
        }
        for name in callbacks:
            invalid = callbacks | {name: object()}
            with self.subTest(name=name), self.assertRaisesRegex(TypeError, name):
                lp.DotCodec(**invalid)

    def test_custom_codec_detached_values_ordering_and_stored_reuse(self):
        seen = []
        codec = round_trip_codec(seen)
        graph = codec_graph(codec)
        encoded = graph.to_dot()
        self.assertIn("digraph", encoded)

        detached_node = next(value for action, value in seen if action == "encode_node")
        detached_edge = next(value for action, value in seen if action == "encode_edge")
        detached_half_edge = next(
            value for action, value in seen if action == "encode_half_edge"
        )
        original_node_data = graph.node("left").data
        detached_node.data = Payload()
        detached_node.drawing.label = lp.TextLabel("detached node")
        detached_edge.drawing.label = "detached edge"
        detached_half_edge.drawing.statement = "detached half-edge"
        self.assertIs(graph.node("left").data, original_node_data)
        self.assertIs(graph.node("left").drawing.label, lp.INHERIT)
        self.assertIs(graph.edge("connection").drawing.label, lp.INHERIT)
        self.assertIs(graph.half_edge(0).drawing.statement, lp.INHERIT)

        decoded = lp.Graph.from_dot(encoded, codec)
        self.assertEqual([node.name for node in decoded.nodes()], ["right", "left"])
        self.assertEqual(decoded.node("right").data["token"], "R")
        self.assertEqual(decoded.node("left").data["token"], "L")
        self.assertEqual(decoded.node("left").data["payload"], b"node payload")
        self.assertEqual(decoded.edge("connection").data["token"], "fermion")
        self.assertEqual(decoded.edge("connection").data["payload"], b"edge payload")
        self.assertEqual(decoded.edge("connection").data["local"], "edge")
        self.assertEqual(decoded.edge(0).source.data["token"], "source")
        self.assertEqual(decoded.edge(0).sink.data["token"], "sink")
        self.assertEqual(decoded.edge(0).source.data["payload"], b"half-edge payload")
        self.assertEqual(decoded.edge(0).source.node.name, "left")
        self.assertEqual(decoded.edge(0).sink.node.name, "right")
        self.assertIn("TextLabel", repr(decoded.node(0).drawing.label))
        self.assertEqual(decoded.edge(0).drawing.extensions["codec-token"], "fermion")

        # Graph.from_dot stores the codec for later exports.
        self.assertIn("digraph", decoded.to_dot())
        detached = [value for action, value in seen if action.startswith("encode_")]
        self.assertTrue(detached)
        self.assertTrue(all(not hasattr(value, "index") for value in detached))

    def test_topology_codec_discards_non_topological_state(self):
        graph, _, _, _ = sample_graph()
        encoded = graph.to_dot(lp.DotCodec.topology())
        decoded = lp.Graph.from_dot(encoded, lp.DotCodec.topology())

        self.assertEqual(
            (decoded.n_nodes, decoded.n_edges, decoded.n_half_edges),
            (2, 2, 3),
        )
        self.assertTrue(all(node.data is None for node in decoded.nodes()))
        self.assertTrue(all(edge.data is None for edge in decoded.edges()))
        self.assertTrue(
            all(half_edge.data is None for half_edge in decoded.half_edges())
        )
        self.assertIn("digraph", decoded.to_dot())

    def test_linnest_codec_round_trips_canonical_drawing_fields(self):
        graph, _, _, _ = sample_graph()
        graph.node("left").drawing.style = {
            "stroke": lp.Stroke(paint=lp.Color("red"), thickness=lp.Length.pt(0.7))
        }
        graph.edge("propagator").drawing.decoration = {
            "marks": [lp.Mark.barbed()],
            "padding": lp.Insets(all=lp.Length.mm(1)),
        }
        graph.edge("propagator").source.drawing.statement = "source statement"
        graph.edge("propagator").source.drawing.port_label = "out"
        graph.edge("propagator").source.drawing.compass = lp.Compass.E
        graph.edge("propagator").sink.drawing.statement = "h_1"
        graph.edge("propagator").sink.drawing.port_label = "in"

        codec = lp.DotCodec.linnest()
        decoded = lp.Graph.from_dot(graph.to_dot(codec), codec)

        self.assertIsNone(decoded.node("left").data)
        self.assertIn("MathSymbol", repr(decoded.node("left").drawing.label))
        self.assertIn("Stroke", repr(decoded.node("left").drawing.style))
        self.assertEqual(
            decoded.edge("propagator").drawing.extensions["category"],
            "primary",
        )
        self.assertEqual(
            decoded.edge("propagator").source.drawing.statement,
            "source statement",
        )
        self.assertEqual(
            decoded.edge("propagator").source.drawing.port_label,
            "out",
        )
        self.assertEqual(
            decoded.edge("propagator").source.drawing.compass,
            lp.Compass.E,
        )
        self.assertEqual(decoded.edge("propagator").sink.drawing.statement, "h_1")
        self.assertEqual(decoded.edge("propagator").sink.drawing.port_label, "in")

    def test_linnest_codec_rejects_module_references_and_tampered_native_values(self):
        codec = lp.DotCodec.linnest()
        item = lp.node("item", extensions={"safe": "value"})
        graph = lp.build(item)
        module = lp.TypstModule.file("styles.typ")

        graph.node("item").drawing.extensions = {
            "module-value": module.value("particle_map")
        }
        with self.assertRaisesRegex(ValueError, "custom DotCodec"):
            graph.to_dot(codec)

        graph.node("item").drawing.extensions = {"safe": "value"}
        encoded = graph.to_dot(codec)
        tampered_values = {
            "unknown color": {
                "type": "dict",
                "value": {
                    "bad": {
                        "type": "color",
                        "value": {"kind": "named", "value": "not-a-color"},
                    }
                },
            },
            "invalid angle unit": {
                "type": "dict",
                "value": {"bad": {"type": "angle", "value": [1.0, "turn"]}},
            },
            "negative fraction": {
                "type": "dict",
                "value": {"bad": {"type": "fraction", "value": -1.0}},
            },
            "empty insets": {
                "type": "dict",
                "value": {"bad": {"type": "insets", "value": {}}},
            },
            "non-length dash phase": {
                "type": "dict",
                "value": {
                    "bad": {
                        "type": "dash",
                        "value": {
                            "kind": "pattern",
                            "value": {
                                "array": [
                                    {
                                        "type": "length",
                                        "value": [1.0, "pt"],
                                    }
                                ],
                                "phase": {
                                    "type": "color",
                                    "value": {
                                        "kind": "named",
                                        "value": "red",
                                    },
                                },
                            },
                        },
                    }
                },
            },
            "module reference": {
                "type": "dict",
                "value": {
                    "bad": {
                        "type": "expression",
                        "value": {
                            "kind": "symbol",
                            "value": {
                                "module": {
                                    "kind": "file",
                                    "value": "styles.typ",
                                },
                                "path": ["particle_map"],
                                "kind": "value",
                            },
                        },
                    }
                },
            },
        }
        for name, value in tampered_values.items():
            tampered = replace_linnest_native(encoded, "linnet_node_extensions", value)
            with self.subTest(name=name), self.assertRaises((TypeError, ValueError)):
                lp.Graph.from_dot(tampered, codec)

    def test_dot_requires_a_codec_and_file_and_set_parsers_store_it(self):
        graph, _, _, _ = sample_graph()
        with self.assertRaises(ValueError):
            graph.to_dot()
        with self.assertRaises(TypeError):
            lp.Graph.from_dot("digraph { a -> b; }")

        codec = lp.DotCodec.topology()
        with TemporaryDirectory() as directory:
            path = Path(directory) / "graphs with spaces.dot"
            path.write_text("digraph one { a -> b; }", encoding="utf-8")
            from_file = lp.Graph.from_dot_file(path, codec)
        self.assertEqual((from_file.n_nodes, from_file.n_edges), (2, 1))
        self.assertIn("digraph", from_file.to_dot())

        graphs = lp.Graph.from_dot_set(
            "digraph one { a -> b; } digraph two { c -> d; }", codec
        )
        self.assertEqual(len(graphs), 2)
        self.assertTrue(all("digraph" in graph.to_dot() for graph in graphs))

    def test_global_data_stays_separate_from_element_codecs(self):
        global_data = lp.GlobalData(
            "global-name",
            payload=b"global payload",
            statements={"label": "global label"},
            edge_statements={"weight": "2"},
            node_statements={"shape": "circle"},
        )
        item = lp.node('isolated "item"')
        graph = lp.build(item, name="python-graph", global_data=global_data)
        codec = lp.DotCodec.topology()
        decoded = lp.Graph.from_dot(graph.to_dot(codec), codec)

        self.assertEqual(decoded.name, "python-graph")
        self.assertEqual(decoded.n_nodes, 1)
        self.assertEqual(decoded.node(0).name, 'isolated "item"')
        self.assertEqual(decoded.global_data.name, "global-name")
        self.assertEqual(decoded.global_data.payload, b"global payload")
        self.assertEqual(decoded.global_data.statements["label"], "global label")
        self.assertEqual(decoded.global_data.edge_statements["weight"], "2")
        self.assertEqual(decoded.global_data.node_statements["shape"], "circle")

    def test_codec_callback_exceptions_propagate(self):
        graph, _, _, _ = sample_graph()

        def fail_encode(_value):
            raise RuntimeError("encode sentinel")

        pass_through = lambda _value: None
        encoding_codec = lp.DotCodec(
            encode_node=fail_encode,
            decode_node=pass_through,
            encode_edge=pass_through,
            decode_edge=pass_through,
            encode_half_edge=pass_through,
            decode_half_edge=pass_through,
        )
        with self.assertRaisesRegex(RuntimeError, "encode sentinel"):
            graph.to_dot(encoding_codec)

        def fail_decode(_value):
            raise LookupError("decode sentinel")

        decoding_codec = lp.DotCodec(
            encode_node=pass_through,
            decode_node=fail_decode,
            encode_edge=pass_through,
            decode_edge=pass_through,
            encode_half_edge=pass_through,
            decode_half_edge=pass_through,
        )
        with self.assertRaisesRegex(LookupError, "decode sentinel"):
            lp.Graph.from_dot("digraph { a -> b; }", decoding_codec)

    def test_encode_callbacks_do_not_borrow_the_live_graph(self):
        replacement = Payload()
        graph = lp.build(lp.node("only", data=Payload()))

        def encode_node(_value):
            graph.node(0).data = replacement
            return lp.DotVertexData()

        codec = lp.DotCodec(
            encode_node=encode_node,
            decode_node=lambda _value: lp.NodeValue(),
            encode_edge=lambda _value: lp.DotEdgeData(),
            decode_edge=lambda _value: lp.EdgeValue(),
            encode_half_edge=lambda _value: lp.DotHalfEdgeData(),
            decode_half_edge=lambda _value: lp.HalfEdgeValue(),
        )

        self.assertIn("digraph", graph.to_dot(codec))
        self.assertIs(graph.node(0).data, replacement)

    def test_codec_statements_escape_quotes_backslashes_and_keys(self):
        key = "hyphen-key"
        value = 'a"b\\N'
        graph = lp.build(
            lp.node('node "name', data=value),
            name='graph "name',
            global_data=lp.GlobalData(statements={key: value}),
        )
        codec = lp.DotCodec(
            encode_node=lambda item: lp.DotVertexData(statements={key: item.data}),
            decode_node=lambda item: lp.NodeValue(data=item.statements[key]),
            encode_edge=lambda _item: lp.DotEdgeData(),
            decode_edge=lambda _item: lp.EdgeValue(),
            encode_half_edge=lambda _item: lp.DotHalfEdgeData(),
            decode_half_edge=lambda _item: lp.HalfEdgeValue(),
        )

        encoded = graph.to_dot(codec)
        decoded = lp.Graph.from_dot(encoded, codec)

        self.assertEqual(decoded.name, 'graph "name')
        self.assertEqual(decoded.node(0).name, 'node "name')
        self.assertEqual(decoded.node(0).data, value)
        self.assertEqual(decoded.global_data.statements[key], value)

        unsupported = lp.build(lp.node("only", data="trailing\\"))
        with self.assertRaisesRegex(ValueError, "cannot end in a backslash"):
            unsupported.to_dot(codec)

    def test_import_rejects_duplicate_python_element_names(self):
        codec = lp.DotCodec.topology()
        duplicates = {
            "node": "digraph { a [linnet_name=78]; b [linnet_name=78]; }",
            "edge": ("digraph { a -> b [linnet_name=78]; b -> c [linnet_name=78]; }"),
        }
        for kind, source in duplicates.items():
            with (
                self.subTest(kind=kind),
                self.assertRaisesRegex(ValueError, f"duplicate imported {kind} name"),
            ):
                lp.Graph.from_dot(source, codec)

    def test_parser_errors_are_python_exceptions(self):
        codec = lp.DotCodec.topology()
        invalid_sources = [
            "digraph { a -> b [dir=sideways]; }",
            "digraph { a [id=invalid]; }",
            "digraph { a [id=0]; b [id=0]; }",
            "digraph { a -> b [id=1]; }",
            "digraph { a:0 -> b:0; }",
        ]
        for source in invalid_sources:
            with self.subTest(source=source), self.assertRaises(ValueError):
                lp.Graph.from_dot(source, codec)


class TestTypedTypstSurface(unittest.TestCase):
    def test_native_values_drawing_and_configuration_smoke(self):
        length = lp.Length.mm(2)
        ratio = lp.Ratio.from_fraction(0.25)
        relative = lp.RelativeLength(ratio, lp.Length.pt(3))
        angle = lp.Angle.degrees(30)
        fraction = lp.Fraction(1)
        color = lp.Color.rgb(10, 20, 30, 200)
        dash = lp.Dash.pattern([lp.Length.pt(2), lp.Length.pt(1)])
        stroke = lp.Stroke(
            paint=color,
            thickness=lp.Length.pt(0.5),
            cap=lp.StrokeCap.Round,
            join=lp.StrokeJoin.Bevel,
            dash=dash,
        )
        insets = lp.Insets(x=length, y=lp.Length.em(0.5))
        mark = lp.Mark(
            end=lp.MarkSymbol.Barbed,
            fill=lp.Color("red"),
            anchor=lp.Anchor.Center,
        )
        label = lp.TextLabel(
            'literal "#(not code)"',
            fill=color,
            size=lp.Length.pt(9),
            style=lp.TextStyle.Italic,
        )
        symbol = lp.MathSymbol("p", subscript="i", superscript="2")

        drawing = lp.EdgeDrawing(
            label=label,
            label_angle=angle,
            bend=angle,
            minimum_length=3,
            style={"stroke": stroke},
            decoration={"dash": dash, "insets": insets, "mark": mark},
            extensions={
                "typed-array": [
                    lp.AUTO,
                    None,
                    symbol,
                    relative,
                    fraction,
                    lp.Pattern.Wave,
                    lp.MarkPosition.Center,
                    lp.MarkOrientation.Path,
                    lp.MarkDirection.Forward,
                    lp.Dash(lp.DashPattern.Dashed),
                ]
            },
        )
        self.assertEqual(drawing.label.text, label.text)
        self.assertIn("Angle", repr(drawing.label_angle))
        self.assertIn("Angle", repr(drawing.bend))
        self.assertEqual(repr(lp.AUTO), "AUTO")
        self.assertEqual(repr(lp.INHERIT), "INHERIT")

        layouts = lp.LayoutOptions(
            algorithm=lp.LayoutAlgorithm.Force,
            direction=lp.LayoutDirection.Right,
            label_layout=lp.LabelLayout.Normal,
            rank_align=lp.RankAlignment.Center,
            roots=[0],
            rank_same=[[0, 1]],
        ).then(
            algorithm=lp.LayoutAlgorithm.Anneal,
            nodes=lp.LayoutNodes.Fixed,
        )
        self.assertEqual(layouts.pass_count, 2)
        style = lp.GraphStyleOptions(
            unit=lp.Length.cm(1),
            node_label=lp.AUTO,
            node_style={"stroke": stroke},
        )
        selectors = lp.DrawingSelectors(
            node=lambda node: lp.NodeDrawing(
                style={"fill": lp.Color("blue") if node.index == 0 else lp.Color("red")}
            ),
            edge=lambda edge: lp.EdgeDrawing(
                label=lp.MathSymbol("e", subscript=edge.index)
            ),
            source=lambda half_edge: lp.HalfEdgeDrawing(
                extensions={"selected-port": half_edge.index}
            ),
            sink=lambda _half_edge: None,
        )
        draw = lp.DrawOptions(
            unit=lp.AUTO,
            node_radius=lp.AUTO,
            edge_stroke=stroke,
            padding=insets,
            draw_node=lp.AUTO,
            edge_resolve_length=lp.EdgeLengthResolution.Min,
            debug=lp.DebugLevel.Off,
        )
        base = lp.RenderConfig(
            title=label,
            style=style,
            layouts=layouts,
            drawing=draw,
            selectors=selectors,
            template_options={
                "profile": "analysis",
                "nested": {"base": True},
                "mark": mark,
            },
        )
        overlay = lp.RenderConfig(
            title=None,
            selectors=lp.INHERIT,
            template_options={"nested": {"overlay": True}},
        )
        merged = base.overlay(overlay)
        self.assertIn("RenderConfig", repr(merged))
        self.assertEqual(
            merged.template_options["nested"],
            {"base": True, "overlay": True},
        )

        graph, _, _, _ = sample_graph(render_config=base)
        graph.render_config = merged
        self.assertIs(graph.render_config, merged)

    def test_sentinels_are_process_singletons(self):
        self.assertIs(copy.copy(lp.AUTO), lp.AUTO)
        self.assertIs(copy.deepcopy(lp.AUTO), lp.AUTO)
        self.assertIs(copy.copy(lp.INHERIT), lp.INHERIT)
        self.assertIs(copy.deepcopy(lp.INHERIT), lp.INHERIT)

        drawing = lp.NodeDrawing(extensions={"automatic": lp.AUTO})
        self.assertIs(drawing.extensions["automatic"], lp.AUTO)
        self.assertIs(lp.NodeDrawing().label, lp.INHERIT)

    def test_render_config_has_mutable_typed_properties(self):
        config = lp.RenderConfig()
        for value in (
            config.template,
            config.title,
            config.style,
            config.layouts,
            config.drawing,
            config.selectors,
            config.template_options,
        ):
            self.assertIs(value, lp.INHERIT)

        config.template = None
        self.assertIsNone(config.template)
        config.template = Path("template with spaces.typ")
        self.assertEqual(Path(config.template), Path("template with spaces.typ"))
        with self.assertRaisesRegex(TypeError, "unknown RenderConfig option"):
            lp.RenderConfig(typst_executable="typst")

        config.title = lp.TextLabel("mutable")
        self.assertEqual(config.title.text, "mutable")
        config.title = lp.AUTO
        self.assertIs(config.title, lp.AUTO)
        config.title = lp.INHERIT
        self.assertIs(config.title, lp.INHERIT)

        config.layouts = lp.LayoutOptions().then(algorithm=lp.LayoutAlgorithm.Force)
        self.assertEqual(config.layouts.pass_count, 2)
        config.drawing = lp.DrawOptions(debug=lp.DebugLevel.Off)
        self.assertIsInstance(config.drawing, lp.DrawOptions)

        class Callback:
            def __call__(self, _node):
                return lp.NodeDrawing(label="selected")

        callback = Callback()
        reference = weakref.ref(callback)
        config.selectors = lp.DrawingSelectors(node=callback)
        del callback
        gc.collect()
        self.assertIsNotNone(reference())
        self.assertIsInstance(config.selectors, lp.DrawingSelectors)
        config.selectors = None
        gc.collect()
        self.assertIsNone(reference())
        self.assertIsNone(config.selectors)

        config.style = lp.GraphStyleOptions(node_style={"fill": lp.Color("blue")})
        self.assertIsInstance(config.style, lp.GraphStyleOptions)
        config.style = None
        self.assertIsNone(config.style)
        config.template_options = {"nested": {"value": 1}}
        snapshot = config.template_options
        snapshot["nested"]["value"] = 2
        self.assertEqual(config.template_options["nested"]["value"], 1)
        config.template_options = None
        self.assertIsNone(config.template_options)

        graph = lp.build(lp.node("only"))
        graph.render_config.title = lp.TextLabel("live graph default")
        self.assertEqual(graph.render_config.title.text, "live graph default")

        for field in (
            "template",
            "layouts",
            "drawing",
            "selectors",
            "template_options",
        ):
            setattr(config, field, lp.INHERIT)
            self.assertIs(getattr(config, field), lp.INHERIT)

    def test_template_options_overlay_and_module_references_are_typed(self):
        module = lp.TypstModule.file("template options.typ")
        base = lp.RenderConfig(
            template_options={
                "nested": {"base": 1, "shared": "base"},
                "accent": module.value("accent"),
                "decorate": module.function("decorate").bind(prefix="graph"),
            }
        )
        merged = base.overlay(
            lp.RenderConfig(
                template_options={
                    "nested": {"overlay": 2, "shared": "overlay"},
                    "title": module.content("title"),
                }
            )
        )

        options = merged.template_options
        self.assertEqual(
            options["nested"],
            {"base": 1, "overlay": 2, "shared": "overlay"},
        )
        self.assertIn("TypstRef", repr(options["accent"]))
        self.assertIn("TypstBind", repr(options["decorate"]))
        self.assertIn("TypstRef", repr(options["title"]))
        preserved = base.overlay(
            lp.RenderConfig(template_options={"nested": {"base": lp.INHERIT}})
        )
        self.assertEqual(preserved.template_options["nested"]["base"], 1)
        self.assertIsNone(
            base.overlay(lp.RenderConfig(template_options=None)).template_options
        )
        with self.assertRaises(TypeError):
            lp.RenderConfig(template_options=["not", "a", "dictionary"])

    def test_drawing_fields_validate_their_declared_types(self):
        module = lp.TypstModule.file("drawing-values.typ")
        value = module.value("value")
        function = module.function("style")
        stroke = lp.Stroke(paint=lp.Color("red"))

        lp.NodeDrawing(
            label=lp.MathSymbol("n", subscript=0),
            shift={"x": 0, "y": value},
            rank=0,
            minimum_size=0.5,
            maximum_size=value,
            style=function,
            label_style={"fill": lp.Color("black")},
        )
        lp.EdgeDrawing(
            label=lp.TextLabel("edge"),
            label_position=(0, 1),
            label_offset=-0.2,
            label_angle=lp.Angle.degrees(30),
            bend=0.1,
            minimum_length=3,
            same_rank=True,
            style=[{"stroke": stroke}],
            label_style={"fill": lp.Color("black")},
            decoration=lp.Pattern.Wave,
            extensions={
                "particle": value,
                "momentum": "p0",
                "cut-id": 0,
                "momentum-style": {"stroke": stroke},
            },
        )
        lp.HalfEdgeDrawing(
            label=lp.MathSymbol("p", subscript=0),
            statement="cut-0",
            port_label="out",
            style={"stroke": stroke},
        )

        invalid = (
            lambda: lp.NodeDrawing(label=1),
            lambda: lp.NodeDrawing(shift=[0]),
            lambda: lp.NodeDrawing(rank=-1),
            lambda: lp.NodeDrawing(minimum_size=-0.1),
            lambda: lp.NodeDrawing(style=stroke),
            lambda: lp.NodeDrawing(label_style=[{}]),
            lambda: lp.EdgeDrawing(label_position=(0, 1, 2)),
            lambda: lp.EdgeDrawing(label_angle="30deg"),
            lambda: lp.EdgeDrawing(bend=lp.RelativeLength(length=lp.Length.pt(1))),
            lambda: lp.EdgeDrawing(minimum_length=lp.Fraction(1)),
            lambda: lp.EdgeDrawing(same_rank=1),
            lambda: lp.EdgeDrawing(style=stroke),
            lambda: lp.EdgeDrawing(decoration=lp.Angle.degrees(5)),
            lambda: lp.HalfEdgeDrawing(statement=lp.MathSymbol("h")),
            lambda: lp.HalfEdgeDrawing(port_label=lp.TextLabel("out")),
            lambda: lp.HalfEdgeDrawing(style=stroke),
        )
        for constructor in invalid:
            with self.subTest(constructor=constructor), self.assertRaises(TypeError):
                constructor()
        for retired_field in (
            {"particle": "electron"},
            {"momentum": "p0"},
            {"cut_id": 0},
            {"momentum_style": {"stroke": stroke}},
        ):
            with (
                self.subTest(field=retired_field),
                self.assertRaisesRegex(ValueError, "template-specific fields"),
            ):
                lp.EdgeDrawing(**retired_field)

    def test_drawing_selectors_validate_results_and_topology_stability(self):
        with TemporaryDirectory(prefix="linnet selectors ") as directory:
            root = Path(directory)
            entrypoint = root / "entrypoint.typ"
            template = root / "template.typ"
            template.write_text("#let render(config) = [ok]\n", encoding="utf-8")

            def compile_typst(input, **kwargs):
                self.assertEqual(kwargs["format"], "svg")
                entrypoint.write_text(Path(input).read_text(encoding="utf-8"))
                return b"<svg>fake</svg>"

            invalid_selectors = (
                (
                    lp.DrawingSelectors(node=lambda _node: lp.EdgeDrawing()),
                    "node drawing selector must return NodeDrawing",
                ),
                (
                    lp.DrawingSelectors(edge=lambda _edge: lp.NodeDrawing()),
                    "edge drawing selector must return EdgeDrawing",
                ),
                (
                    lp.DrawingSelectors(source=lambda _half_edge: lp.EdgeDrawing()),
                    "half-edge drawing selector must return HalfEdgeDrawing",
                ),
                (
                    lp.DrawingSelectors(sink=lambda _half_edge: lp.AUTO),
                    "half-edge drawing selector must return HalfEdgeDrawing",
                ),
            )
            for selectors, message in invalid_selectors:
                graph, _, _, _ = sample_graph()
                with (
                    self.subTest(selectors=selectors),
                    self.assertRaisesRegex(TypeError, message),
                ):
                    graph.to_svg(
                        config=lp.RenderConfig(
                            template=template,
                            selectors=selectors,
                        )
                    )

            node_data = [UserNodePayload("left"), UserNodePayload("right")]
            edge_data = UserEdgePayload("depends-on")
            source_data = UserHalfEdgePayload("producer")
            sink_data = UserHalfEdgePayload("consumer")
            left = lp.node("left", data=node_data[0], label="explicit left")
            right = lp.node("right", data=node_data[1])
            connection = lp.edge(
                lp.source(left, data=source_data, statement="explicit source"),
                "connection",
                lp.sink(right, data=sink_data),
                data=edge_data,
                label="explicit edge",
                extensions={"priority": "explicit"},
            )
            graph = lp.build(left, right, connection)
            captured = []

            def select_node(node):
                captured.append(node)
                self.assertIs(node.data, node_data[node.index])
                return lp.NodeDrawing(
                    label=f"selected {node.data.identifier}",
                    extensions={"selected-node-data": node.data.identifier},
                )

            def select_edge(edge):
                captured.append(edge)
                self.assertIs(edge.data, edge_data)
                return lp.EdgeDrawing(
                    label="selected edge",
                    extensions={
                        "priority": "selector",
                        "selected-coupling": edge.data.coupling,
                    },
                )

            def select_source(half_edge):
                captured.append(half_edge)
                self.assertIs(half_edge.data, source_data)
                return lp.HalfEdgeDrawing(
                    statement="selected source",
                    extensions={"selected-port-data": half_edge.data.port},
                )

            def select_sink(half_edge):
                captured.append(half_edge)
                self.assertIs(half_edge.data, sink_data)
                return lp.HalfEdgeDrawing(
                    statement="selected sink",
                    extensions={"selected-port-data": half_edge.data.port},
                )

            selectors = lp.DrawingSelectors(
                node=select_node,
                edge=select_edge,
                source=select_source,
                sink=select_sink,
            )
            with patch.object(typst, "compile", side_effect=compile_typst):
                self.assertEqual(
                    graph.to_svg(
                        config=lp.RenderConfig(
                            template=template,
                            selectors=selectors,
                        )
                    ),
                    "<svg>fake</svg>",
                )
            source = entrypoint.read_text(encoding="utf-8")
            self.assertIn("explicit left", source)
            self.assertNotIn("selected left", source)
            self.assertIn("selected right", source)
            self.assertIn("explicit edge", source)
            self.assertNotIn("selected edge", source)
            self.assertIn("selected-coupling", source)
            self.assertIn('("priority"): "explicit"', source)
            self.assertNotIn('("priority"): "selector"', source)
            self.assertIn("explicit source", source)
            self.assertNotIn("selected source", source)
            self.assertIn("selected sink", source)
            self.assertIn("producer", source)
            self.assertIn("consumer", source)

            graph.reorder_nodes([1, 0])
            for view in captured:
                with self.assertRaises(ReferenceError):
                    _ = view.data

            mutating, _, _, _ = sample_graph()

            def mutate_topology(_node):
                mutating.reorder_nodes([1, 0])
                return lp.NodeDrawing()

            with self.assertRaisesRegex(
                ReferenceError,
                "must not mutate graph topology",
            ):
                mutating.to_svg(
                    config=lp.RenderConfig(
                        template=template,
                        selectors=lp.DrawingSelectors(node=mutate_topology),
                    )
                )

    def test_complete_option_field_surfaces(self):
        stroke = lp.Stroke(
            paint=lp.Color("black"),
            thickness=lp.Length.pt(0.6),
            cap=lp.StrokeCap.Round,
            join=lp.StrokeJoin.Miter,
            dash=lp.Dash(lp.DashPattern.DenselyDashed),
            miter_limit=4,
        )
        style = lp.GraphStyleOptions(
            scope={"accent": lp.Color("blue")},
            unit=lp.RelativeLength(lp.Ratio(100), lp.Length.pt(1)),
            node_label=lp.TextLabel("node"),
            node_label_style={"fill": lp.Color("black")},
            node_style={"stroke": stroke},
            edge_label=lp.MathSymbol("p", subscript=0),
            edge_label_style={"fill": lp.Color("red")},
        )
        selectors = lp.DrawingSelectors(
            node=lambda _node: lp.NodeDrawing(style={"fill": lp.Color("white")}),
            edge=lambda _edge: lp.EdgeDrawing(style={"stroke": stroke}),
            source=lambda _half_edge: lp.HalfEdgeDrawing(style={"stroke": stroke}),
            sink=lambda _half_edge: lp.HalfEdgeDrawing(style={"stroke": stroke}),
        )
        layouts = lp.LayoutOptions(
            subgraph=[True, False],
            viewport_width=12,
            viewport_height=8,
            tree_dx=1.1,
            tree_dy=1.2,
            steps=20,
            seed=7,
            step=0.5,
            step_shrink=0.2,
            cool=0.8,
            accept_floor=0.1,
            early_tolerance=1e-6,
            temperature=0.3,
            delta=0.4,
            beta=40,
            spring_strength=10,
            centering_strength=0.01,
            epochs=10,
            crossing_penalty=20,
            dangling_repulsion=2,
            edge_edge_repulsion=0.2,
            directional_force=4.5,
            label_length_scale=1.1,
            label_spring=20,
            label_charge=2,
            label_steps=30,
            label_layout=lp.LabelLayout.DanglingTangent,
            label_step=0.1,
            label_early_tolerance=1e-3,
            label_max_delta_scale=0.4,
            edge_vertex_repulsion=0.05,
            epsilon=1e-4,
            incremental_energy=True,
            algorithm=lp.LayoutAlgorithm.StableLayered,
            nodes=lp.LayoutNodes.Layout,
            direction=lp.LayoutDirection.Down,
            rank_align=lp.RankAlignment.Start,
            roots=[0, 1],
            rank_same=[[0, 1]],
            route_edge_weight=0.2,
            route_exit_weight=3,
            route_label_width_scale=1.1,
            route_label_width_cap=2,
            z_spring=2,
            z_spring_growth=1.01,
            length_scale=0.4,
        )
        drawing = lp.DrawOptions(
            scope={"accent": lp.Color("blue")},
            unit=lp.Length.em(1),
            title=lp.TextLabel("all drawing fields"),
            subgraph=[[True, False], [False, True]],
            debug=lp.DebugLevel.EdgePositions,
            node_radius=[0.2, 0.3],
            node_min_radius=0.1,
            node_label_padding=0.08,
            node_fill=lp.Color("white"),
            node_stroke=stroke,
            node_outset=0.2,
            node_label_style={"fill": lp.Color("black")},
            node_style={"stroke": stroke},
            node_label=lp.TextLabel("node"),
            draw_node=lp.AUTO,
            edge_stroke=stroke,
            edge_offset=0.1,
            edge_length=4,
            edge_ratio=0.8,
            edge_resolve_length=lp.EdgeLengthResolution.Shorter,
            edge_accuracy=0.001,
            edge_optimize=True,
            source_style=[{"stroke": stroke}],
            sink_style={"stroke": stroke},
            edge_label=lp.MathSymbol("p", subscript=1),
            edge_label_style={"fill": lp.Color("red")},
            edge_omega=1,
            edge_trim_accuracy=0.001,
            padding=lp.Insets(x=lp.Length.em(0.4), y=lp.Length.em(0.2)),
            debug_edge_radius=0.08,
            debug_edge_fill=lp.Color("orange"),
            debug_edge_stroke=stroke,
            debug_edge_label_fill=lp.Color("red"),
            subgraph_edge_style={"stroke": stroke},
            subgraph_edge_underlay=True,
        )
        template_options = {
            "palette": {"accent": lp.Color("blue")},
            "orientation": {"split": True},
            "edge-decoration": {
                "offset": 0.4,
                "length": 5,
                "ratio": 0.5,
                "stroke": stroke,
                "mark": lp.Mark.barbed(),
            },
            "annotations": {
                "node-indices": True,
                "edge-indices": True,
                "half-edge-indices": True,
                "prefix": lp.MathSymbol("q"),
                "separator": lp.TextLabel(", "),
                "size": lp.Length.pt(7),
                "fill": lp.Color("red"),
            },
        }
        config = lp.RenderConfig(
            template=Path("template with spaces.typ"),
            title=lp.AUTO,
            style=style,
            layouts=layouts,
            drawing=drawing,
            selectors=selectors,
            template_options=template_options,
        )
        self.assertEqual(layouts.pass_count, 1)
        self.assertIn("RenderConfig", repr(config))

    def test_closed_native_model_rejects_arbitrary_values_and_reserved_extensions(self):
        nested = [lp.Color("red")]
        drawing = lp.NodeDrawing(style={"fills": nested})
        nested.append(Payload())
        self.assertEqual(len(drawing.style["fills"]), 1)

        with self.assertRaises(TypeError):
            lp.NodeDrawing(label=Payload())
        with self.assertRaises(TypeError):
            lp.NodeDrawing(extensions={"value": Payload()})
        with self.assertRaises((TypeError, ValueError)):
            lp.NodeDrawing(extensions={1: "not a string key"})
        with self.assertRaises(OverflowError):
            lp.NodeDrawing(extensions={"too-large": 2**63})
        with self.assertRaises(OverflowError):
            lp.NodeDrawing(extensions={"too-small": -(2**63) - 1})
        with self.assertRaises(ValueError):
            lp.NodeDrawing(extensions={"label": "reserved"})
        for name in ("node-style", "pos"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                lp.NodeDrawing(extensions={name: lp.Color("red")})
        selector_extension = lp.NodeDrawing(
            extensions={"selector-style": lp.Color("red")}
        )
        self.assertIn("Color", repr(selector_extension.extensions["selector-style"]))
        with self.assertRaises(ValueError):
            lp.NodeDrawing(extensions={"statements": {}})
        for name in ("shift", "statements"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                lp.EdgeDrawing(extensions={name: {}})
        with self.assertRaises(TypeError):
            lp.NodeDrawing(placement="start")
        lp.NodeDrawing(placement={"x": 0.0, "y": 1.0, "mode": lp.Placement.Start})
        lp.NodeDrawing(placement=lp.Placement.Start)
        with self.assertRaises(ValueError):
            lp.NodeDrawing(placement=lp.Placement.Pin)
        with self.assertRaises(TypeError):
            lp.NodeDrawing(placement={"mode": None})
        with self.assertRaises(TypeError):
            lp.NodeDrawing(placement={"x": "bad", "mode": lp.Placement.Pin})
        with self.assertRaises(TypeError):
            lp.NodeDrawing(placement={"ref": "left", "mode": lp.Placement.Pin})
        with self.assertRaises(TypeError):
            lp.LayoutOptions(algorithm="force")
        with self.assertRaises(TypeError):
            lp.LayoutOptions(seed=-1)
        with self.assertRaises(TypeError):
            lp.LayoutOptions(roots=[-1])
        with self.assertRaises(TypeError):
            lp.LayoutOptions(rank_same=[["node"]])
        with self.assertRaises(TypeError):
            lp.LayoutOptions(subgraph=[0, 1])
        lp.LayoutOptions(subgraph=[True, False], rank_same=[[0, 1]])
        module = lp.TypstModule.file("selections.typ")
        selection = module.function("selection")
        lp.LayoutOptions(
            subgraph=selection,
            rank_same=[[0, 1], selection, selection.bind(side="left")],
        )
        with self.assertRaises(TypeError):
            lp.GraphStyleOptions(unit=lp.AUTO)
        with self.assertRaises(TypeError):
            lp.DrawOptions(node_outset=lp.Length.pt(1))
        lp.DrawOptions(subgraph=[[True, False], [False, True]])
        lp.DrawOptions(
            subgraph=[
                selection,
                {
                    "subgraph": [True, False],
                    "edge-style": {"stroke": lp.Color("red")},
                },
                {
                    "subgraph": selection.bind(side="right"),
                    "edge-style": selection,
                },
            ]
        )
        with self.assertRaises(TypeError):
            lp.DrawOptions(subgraph=[{"edge-style": {}}])
        with self.assertRaises(TypeError):
            lp.DrawOptions(subgraph=[{"subgraph": [True, False], "edge-style": "red"}])
        with self.assertRaises(TypeError):
            lp.DrawOptions(subgraph=[{"subgraph": [True, False], "unknown": {}}])
        static_content_fields = {
            "draw title": lambda value: lp.DrawOptions(title=value),
            "render title": lambda value: lp.RenderConfig(title=value),
        }
        for name, static_content in static_content_fields.items():
            with self.subTest(field=name), self.assertRaises(TypeError):
                static_content(selection)
            with self.subTest(field=name), self.assertRaises(TypeError):
                static_content(selection.bind(side="left"))
            static_content(module.content("title"))
        lp.GraphStyleOptions(node_label=selection, edge_label=selection)
        lp.DrawOptions(node_label=selection, edge_label=selection.bind(side="left"))
        with self.assertRaises(TypeError):
            lp.DrawOptions(draw_node=selection.call())
        with self.assertRaises(TypeError):
            lp.MathSymbol("p", subscript=True)
        with self.assertRaises(TypeError):
            lp.RenderConfig(unknown_option=True)

    def test_module_references_calls_and_binds_are_typed_values(self):
        module = lp.TypstModule.file("styles.typ")
        function = module.function("edge_style")
        drawing = lp.EdgeDrawing(
            style=function.bind(fill=lp.Color("red")),
            decoration={
                "call": function.call(lp.Length.pt(1), fill=lp.Color("blue")),
                "content": module.content("title"),
                "value": module.value("particle_map"),
            },
        )
        self.assertIn("TypstBind", repr(drawing.style))
        self.assertIn("TypstCall", repr(drawing.decoration["call"]))
        lp.DrawOptions(edge_resolve_length=function)
        for reference in [module.value("particle_map"), module.content("title")]:
            with self.subTest(reference=reference), self.assertRaises(TypeError):
                reference.call()
            with self.subTest(reference=reference), self.assertRaises(TypeError):
                reference.bind(fill="red")

    def test_typst_package_references_validate_every_component(self):
        self.assertIn(
            "TypstModule",
            repr(lp.TypstModule.package("@preview/cetz:0.5.1")),
        )
        for specification in [
            "preview/cetz:0.5.1",
            "@preview:cetz",
            "@preview/:0.5.1",
            "@preview/cetz:",
            "@preview/cetz:0.5.1/extra",
        ]:
            with (
                self.subTest(specification=specification),
                self.assertRaises(ValueError),
            ):
                lp.TypstModule.package(specification)

    def test_legacy_names_and_render_inputs_are_gone(self):
        legacy = {
            "DotGraph",
            "DotGraphBuilder",
            "Hedge",
            "NodeIndex",
            "EdgeIndex",
            "HedgePair",
        }
        self.assertTrue(all(not hasattr(lp, name) for name in legacy))
        topology_views = {
            "Subgraph",
            "Cycle",
            "OrientedCut",
            "CutPartition",
            "TraversalTree",
        }
        self.assertTrue(all(hasattr(lp, name) for name in topology_views))
        self.assertFalse(hasattr(lp.Graph, "from_string"))
        graph, _, _, _ = sample_graph()
        self.assertFalse(hasattr(graph, "dot"))
        with self.assertRaises(TypeError):
            lp.Graph(data=Payload())
        with self.assertRaises(AttributeError):
            graph.data = Payload()
        with self.assertRaises(TypeError):
            graph.to_svg(inputs={"raw": "Typst source"})


class TestRendering(unittest.TestCase):
    def test_typst_py_version_matches_the_supported_typst_runtime(self):
        self.assertEqual(typst.__version__, "0.15.0")

    def test_render_calls_typst_py_with_format_and_environment(self):
        with TemporaryDirectory(prefix="linnet typst py ") as directory:
            root = Path(directory)
            template = root / "template.typ"
            template.write_text("#let render(config) = [ok]\n", encoding="utf-8")
            graph = lp.build(
                lp.node("only"),
                render_config=lp.RenderConfig(template=template),
            )
            calls = []

            def compile_typst(input, **kwargs):
                self.assertTrue(Path(input).is_file())
                self.assertTrue(Path(kwargs["root"]).is_dir())
                Path(kwargs["output"]).write_bytes(b"rendered")
                calls.append(kwargs)

            environment = {
                "TYPST_PACKAGE_PATH": str(root / "local packages"),
                "TYPST_PACKAGE_CACHE_PATH": str(root / "package cache"),
                "TYPST_FONT_PATHS": os.pathsep.join(
                    (str(root / "font one"), str(root / "font two"))
                ),
            }
            with (
                patch.dict(os.environ, environment),
                patch.object(typst, "compile", side_effect=compile_typst),
            ):
                for suffix in ("pdf", "svg", "png"):
                    output = root / "nested output" / f"diagram.{suffix}"
                    self.assertEqual(graph.render(output), output)
                    self.assertEqual(output.read_bytes(), b"rendered")

            self.assertEqual([call["format"] for call in calls], ["pdf", "svg", "png"])
            for call in calls:
                self.assertEqual(Path(call["package_path"]), root / "local packages")
                self.assertEqual(
                    Path(call["package_cache_path"]), root / "package cache"
                )
                self.assertEqual(
                    list(map(Path, call["font_paths"])),
                    [root / "font one", root / "font two"],
                )
            with self.assertRaisesRegex(RuntimeError, "expected a .pdf, .svg, or .png"):
                graph.render(root / "diagram.jpg")

    def test_to_svg_accepts_one_page_and_rejects_multiple_pages(self):
        with TemporaryDirectory(prefix="linnet typst pages ") as directory:
            template = Path(directory) / "template.typ"
            template.write_text("#let render(config) = [ok]\n", encoding="utf-8")
            graph = lp.build(
                lp.node("only"),
                render_config=lp.RenderConfig(template=template),
            )

            with patch.object(
                typst,
                "compile",
                return_value=[b"<svg>one page</svg>"],
            ):
                self.assertEqual(graph.to_svg(), "<svg>one page</svg>")
            with (
                patch.object(
                    typst,
                    "compile",
                    return_value=[b"<svg>one</svg>", b"<svg>two</svg>"],
                ),
                self.assertRaisesRegex(RuntimeError, "produced 2 pages"),
            ):
                graph.to_svg()

    def test_render_stages_structural_names_without_python_data(self):
        class Opaque:
            def __str__(self):
                raise AssertionError("renderer must not stringify Python data")

            def __repr__(self):
                raise AssertionError("renderer must not inspect Python data")

        with TemporaryDirectory(prefix="linnet topology ") as directory:
            root = Path(directory)
            staged_dot = root / "staged topology.dot"
            template = root / "template.typ"
            template.write_text("#let render(config) = [ok]\n", encoding="utf-8")

            def compile_typst(input, **kwargs):
                self.assertEqual(kwargs["format"], "svg")
                staged_dot.write_bytes(
                    Path(input).with_name("diagram.dot").read_bytes()
                )
                return b"<svg>fake</svg>"

            left = lp.node("left node", data=Opaque())
            right = lp.node('right "node"', data=Opaque())
            graph = lp.build(
                left,
                right,
                lp.edge(
                    lp.source(left, data=Opaque()),
                    "propagator edge",
                    lp.sink(right, data=Opaque()),
                    data=Opaque(),
                ),
                name="render graph",
                global_data=lp.GlobalData(
                    "DOT_CODEC_STATE_MUST_NOT_STAGE",
                    payload=b"DOT_CODEC_PAYLOAD_MUST_NOT_STAGE",
                    statements={"codec-global": "DOT_CODEC_VALUE_MUST_NOT_STAGE"},
                    edge_statements={"codec-edge": "DOT_CODEC_EDGE_MUST_NOT_STAGE"},
                    node_statements={"codec-node": "DOT_CODEC_NODE_MUST_NOT_STAGE"},
                ),
                render_config=lp.RenderConfig(template=template),
                node_store=lp.NodeStore.Forest,
            )

            self.assertEqual(graph.node_store, lp.NodeStore.Forest)
            with patch.object(typst, "compile", side_effect=compile_typst):
                self.assertEqual(graph.to_svg(), "<svg>fake</svg>")
            topology = staged_dot.read_text(encoding="utf-8")
            self.assertIn("left node", topology)
            self.assertIn(r"right \"node\"", topology)
            self.assertIn("linnet_render_edge_name", topology)
            self.assertIn("propagator edge", topology)
            self.assertIn("render graph", topology)
            self.assertNotIn("MUST_NOT_STAGE", topology)

    def test_default_renderer_accepts_typed_placement(self):
        left = lp.node("left", placement=lp.Placement.Start)
        right = lp.node(
            "right",
            placement={"x": 2.0, "y": 0.0, "mode": lp.Placement.Pin},
        )
        graph = lp.build(
            left,
            right,
            lp.edge(lp.source(left), "line", lp.sink(right)),
        )
        subgraph = [True] * graph.n_half_edges
        graph.render_config = lp.RenderConfig(
            layouts=lp.LayoutOptions(
                algorithm=lp.LayoutAlgorithm.StableLayered,
                roots=[0],
                rank_same=[[0, 1]],
                subgraph=subgraph,
            ).then(
                algorithm=lp.LayoutAlgorithm.Force,
                nodes=lp.LayoutNodes.Fixed,
                steps=2,
                subgraph=subgraph,
            ),
            drawing=lp.DrawOptions(subgraph=subgraph),
            selectors=lp.DrawingSelectors(
                node=lambda node: lp.NodeDrawing(
                    label=lp.MathSymbol("n", subscript=node.index)
                ),
                edge=lambda edge: lp.EdgeDrawing(
                    label=lp.MathSymbol("e", subscript=edge.index)
                ),
            ),
        )
        with patch.dict(os.environ, {"PATH": ""}):
            self.assertIn("<svg", graph.to_svg())

    def test_notebook_display_is_generic_and_specialized_options_are_custom(self):
        node_data = {
            name: UserNodePayload(name)
            for name in ("ingest", "validate", "transform", "publish", "archive")
        }
        nodes = {
            name: lp.node(name, data=data, label=lp.TextLabel(name.title()))
            for name, data in node_data.items()
        }

        def dependency(name, source_name, sink_name):
            return lp.edge(
                lp.source(
                    nodes[source_name],
                    data=UserHalfEdgePayload(f"{name}-source"),
                ),
                name,
                lp.sink(
                    nodes[sink_name],
                    data=UserHalfEdgePayload(f"{name}-sink"),
                ),
                data=UserEdgePayload(name),
                extensions={"category": "dependency"},
            )

        graph = lp.build(
            *nodes.values(),
            dependency("parse", "ingest", "validate"),
            dependency("clean", "validate", "transform"),
            dependency("retry", "transform", "validate"),
            dependency("release", "transform", "publish"),
            dependency("snapshot", "transform", "archive"),
            dependency("audit", "archive", "publish"),
            render_config=lp.RenderConfig(
                layouts=lp.LayoutOptions(
                    algorithm=lp.LayoutAlgorithm.Force,
                    seed=19,
                    steps=80,
                    label_steps=20,
                ),
            ),
        )

        neutral = graph._repr_svg_()
        self.assertEqual(neutral, graph.to_svg(config=lp.RenderConfig()))
        self.assertIn("<svg", neutral)

        with TemporaryDirectory(prefix="linnet custom options ") as directory:
            root = Path(directory)
            module_path = root / "workflow style.typ"
            module_path.write_text(
                '#let accent = rgb("#356a9a")\n'
                "#let title = [Workflow graph]\n"
                "#let decorate(body, prefix: [Graph]) = [#prefix: #body]\n",
                encoding="utf-8",
            )
            template = root / "workflow template.typ"
            template.write_text(
                "#let render(config) = {\n"
                '  let options = config.at("options")\n'
                '  let diagram = options.at("diagram")\n'
                '  let physics = options.at("physics")\n'
                '  assert.eq(diagram.at("kind"), "workflow")\n'
                '  assert.eq(config.elements.edges.at(0).at("category"), "dependency")\n'
                '  let decorate = options.at("decorate")\n'
                '  let state = if physics.at("momentum-arrows") { [enabled] } else { [disabled] }\n'
                '  set text(fill: options.at("accent"))\n'
                '  decorate([#(options.at("title")): physics #state])\n'
                "}\n",
                encoding="utf-8",
            )
            module = lp.TypstModule.file(module_path)
            graph.render_config = lp.RenderConfig(
                template=template,
                template_options={
                    "diagram": {"kind": "workflow"},
                    "physics": {"momentum-arrows": False},
                    "accent": module.value("accent"),
                    "title": module.content("title"),
                    "decorate": module.function("decorate").bind(prefix="Graph"),
                },
            )
            disabled = graph._repr_svg_()
            enabled = graph.to_svg(
                config=lp.RenderConfig(
                    template_options={"physics": {"momentum-arrows": True}}
                )
            )
            self.assertIn("<svg", disabled)
            self.assertIn("<svg", enabled)
            self.assertNotEqual(disabled, enabled)

    def test_real_typst_preserves_structural_names_and_half_edge_fields(self):
        repository = Path(__file__).resolve().parents[3]
        with TemporaryDirectory(prefix="linnet names ") as directory:
            root = Path(directory)
            package = root / "linnest"
            shutil.copytree(repository / "crates/linnest/typst/src", package / "src")
            shutil.copy2(
                repository / "crates/linnest/typst/linnest.wasm",
                package / "linnest.wasm",
            )
            graph_module = lp.TypstModule.file(package / "src/graph.typ")
            template = root / "inspect names.typ"
            template.write_text(
                "#let render(config) = {\n"
                "  let helpers = config.elements.nodes.at(0)\n"
                '  let parsed = helpers.at("inspect-parse")(read(config.at("data-path"))).first()\n'
                "  let attach = half-edge => {\n"
                "    let drawing = config.elements.hedges.at(half-edge.hedge)\n"
                "    let patch = (data: drawing)\n"
                '    for key in ("statement", "port-label", "compass") {\n'
                "      if drawing.keys().contains(key) { patch.insert(key, drawing.at(key)) }\n"
                "    }\n"
                "    patch\n"
                "  }\n"
                '  parsed = helpers.at("inspect-map")(parsed, source: attach, sink: attach)\n'
                '  let nodes = helpers.at("inspect-nodes")(parsed)\n'
                '  let edges = config.elements.edges.at(0).at("inspect-edges")(parsed)\n'
                "  let node-names = nodes.map(node => str(node.name))\n"
                "  let edge-names = edges.map(edge => str(edge.name))\n"
                '  assert(node-names == ("left node", "right \\"node\\""))\n'
                '  assert(edge-names == ("propagator edge",))\n'
                '  assert(edges.first().source.statement == "source-statement")\n'
                '  assert(edges.first().source.at("port-label") == "source-port")\n'
                '  assert(edges.first().source.compass == "e")\n'
                '  assert(edges.first().sink.statement == "sink-statement")\n'
                '  assert(edges.first().sink.at("port-label") == "sink-port")\n'
                '  assert(edges.first().sink.compass == "w")\n'
                '  assert(config.elements.nodes.at(nodes.at(0).node).token == "left-node")\n'
                '  assert(config.elements.nodes.at(nodes.at(1).node).token == "right-node")\n'
                '  assert(config.elements.edges.at(edges.first().edge).token == "edge")\n'
                '  assert(config.elements.hedges.at(edges.first().source.hedge).token == "source")\n'
                '  assert(config.elements.hedges.at(edges.first().sink.hedge).token == "sink")\n'
                "  set page(width: auto, height: auto, margin: 2mm)\n"
                "  [structural names preserved]\n"
                "}\n",
                encoding="utf-8",
            )
            left = lp.node(
                "left node",
                extensions={
                    "inspect-parse": graph_module.function("parse"),
                    "inspect-map": graph_module.function("map"),
                    "inspect-nodes": graph_module.function("nodes"),
                    "token": "left-node",
                },
            )
            right = lp.node('right "node"', extensions={"token": "right-node"})
            graph = lp.build(
                left,
                right,
                lp.edge(
                    lp.source(
                        left,
                        statement="source-statement",
                        port_label="source-port",
                        compass=lp.Compass.E,
                        extensions={"token": "source"},
                    ),
                    "propagator edge",
                    lp.sink(
                        right,
                        statement="sink-statement",
                        port_label="sink-port",
                        compass=lp.Compass.W,
                        extensions={"token": "sink"},
                    ),
                    extensions={
                        "inspect-edges": graph_module.function("edges"),
                        "token": "edge",
                    },
                ),
                render_config=lp.RenderConfig(template=template),
            )

            self.assertIn("<svg", graph.to_svg())

    def test_real_typst_v1_template_module_paths_and_high_level_methods(self):
        with TemporaryDirectory(prefix="linnet py ") as directory:
            root = Path(directory) / "project with spaces"
            root.mkdir()
            module_path = root / "drawing styles.typ"
            module_path.write_text(
                "#let title = [trusted module content]\n"
                "#let wrap(body, prefix: [prefix]) = [#prefix #body]\n"
                "#let underscore(foo_bar: [missing]) = foo_bar\n"
                '#let particle_map = (electron: "e-",)\n',
                encoding="utf-8",
            )
            module = lp.TypstModule.file(module_path)
            template = root / "custom template.typ"
            template.write_text(
                "#let render(config) = {\n"
                '  assert(config.schema == "linnet-render-config")\n'
                "  assert(config.version == 1)\n"
                "  set page(width: auto, height: auto, margin: 2mm)\n"
                '  config.elements.nodes.at(0).at("underscore-call")\n'
                "  config.elements.nodes.at(0).label\n"
                "}\n",
                encoding="utf-8",
            )
            title = module.function("wrap").bind(prefix=lp.TextLabel("typed"))
            first = lp.node(
                "first",
                label=title.call(module.content("title")),
                extensions={
                    "particle-map": module.value("particle_map"),
                    "underscore-call": module.function("underscore").call(
                        foo_bar=lp.TextLabel("underscore")
                    ),
                },
            )
            graph = lp.build(
                first,
                render_config=lp.RenderConfig(template=template),
            )

            output = root / "rendered diagram.svg"
            self.assertEqual(graph.render(output), output)
            self.assertIn("<svg", output.read_text(encoding="utf-8"))
            self.assertIn("<svg", graph.to_svg())
            self.assertIn("<svg", graph._repr_svg_())

    def test_large_config_is_staged_for_in_process_compilation(self):
        with TemporaryDirectory(prefix="linnet transport ") as directory:
            root = Path(directory)
            entrypoint = root / "entrypoint.typ"
            template = root / "template.typ"
            template.write_text("#let render(config) = [ok]\n", encoding="utf-8")
            module_path = root / "drawing styles.typ"
            module_path.write_text(
                '#let particle_map = (electron: "e-",)\n',
                encoding="utf-8",
            )
            module = lp.TypstModule.file(module_path)
            graph, _, _, _ = sample_graph(
                render_config=lp.RenderConfig(
                    template=template,
                    title=lp.TextLabel('")\n#let injected = true\n' + "x" * 100_000),
                    selectors=lp.DrawingSelectors(
                        node=lambda node: lp.NodeDrawing(
                            extensions={"selected-name": lp.TextLabel(node.name)}
                        )
                    ),
                )
            )
            graph.node("left").drawing.extensions = {
                "named-call": module.function("named").call(
                    **{"foo_bar": 1, "foo-bar": 2}
                ),
                "particle-map": module.value("particle_map"),
            }
            graph.node("right").drawing.extensions = {
                "particle-map": module.value("particle_map")
            }

            def compile_typst(input, **kwargs):
                self.assertEqual(kwargs["format"], "svg")
                entrypoint.write_text(Path(input).read_text(encoding="utf-8"))
                return b"<svg>fake</svg>"

            with patch.object(
                typst, "compile", side_effect=compile_typst
            ) as compile_mock:
                self.assertEqual(graph.to_svg(), "<svg>fake</svg>")
            compile_mock.assert_called_once()
            self.assertGreater(entrypoint.stat().st_size, 100_000)
            source = entrypoint.read_text(encoding="utf-8")
            self.assertEqual(source.count("drawing styles.typ"), 1)
            self.assertIn("foo-bar: 2", source)
            self.assertIn("foo_bar: 1", source)
            self.assertIn("selected-name", source)
            self.assertNotIn("\n#let injected = true", source)
            self.assertIn("\\n#let injected = true\\n", source)


if __name__ == "__main__":
    unittest.main()
