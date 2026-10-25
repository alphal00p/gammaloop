"""Ordinary graphs share the native SVG path used by domain drawings."""

import json
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from tempfile import TemporaryDirectory

from rendering_metrics import external_metrics_snapshot
from symbolica.community import graph


class NativeRenderingTests(unittest.TestCase):
    def setUp(self):
        self.graph = graph.build(
            graph.node("a", label=graph.TextLabel("left")),
            graph.node("b", label=graph.MathSymbol("p", subscript=1)),
            graph.node("isolated"),
            graph.edge(graph.source("a"), "ab", graph.sink("b")),
            graph.edge(graph.source("a"), "parallel", graph.sink("b")),
            graph.edge(graph.source("b"), "loop", graph.sink("b")),
            graph.edge(graph.source("a"), "open"),
        )

    def assert_native(self, result):
        svg = result.to_svg()
        self.assertIn("data-linnet-kind", svg)
        self.assertEqual(ET.fromstring(svg).get("data-linnet-renderer"), "native")
        # MathML snapshots require an explicit outlined export for Typst.
        with self.assertRaisesRegex(ValueError, "outlined labels"):
            result.to_linnest()
        root = ET.fromstring(svg)
        for target in root.iter():
            kind = target.get("data-linnet-kind")
            if kind in {"node", "edge"}:
                element = getattr(self.graph, kind)(int(target.get("data-linnet-id")))
                self.assertEqual(
                    json.loads(target.get("data-linnet-detail"))["name"], element.name
                )
        return root

    def test_default_graph_and_subgraph_use_native_geometry(self):
        for owner in (
            self.graph,
            self.graph.subgraph(edges=["ab"]),
            self.graph.subgraph(nodes=["isolated"]),
            self.graph.empty_subgraph(),
        ):
            with self.subTest(owner=type(owner).__name__):
                view = owner.render()
                self.assertIs(type(view), graph.DiagramView)
                result = external_metrics_snapshot(view)
                self.assert_native(result)
                self.assertEqual(result.to_svg(), result.to_svg_pages()[0])
                self.assertEqual(result._mime_(), ("text/html", result.to_html()))
        for owner, highlighted in (
            (self.graph.subgraph(edges=["ab"]), True),
            (self.graph.empty_subgraph(), False),
        ):
            with self.subTest(highlighted=highlighted):
                svg = external_metrics_snapshot(owner.render()).to_svg()
                painted = []
                for element in ET.fromstring(svg).iter():
                    if element.tag.rsplit("}", 1)[-1] == "style":
                        continue
                    paints = [element.get("stroke", ""), element.get("fill", "")]
                    for declaration in element.get("style", "").split(";"):
                        property_name, _, value = declaration.partition(":")
                        if property_name.strip().lower() in {"stroke", "fill"}:
                            paints.append(value)
                    if any(paint.strip().lower() == "#c58b13" for paint in paints):
                        painted.append(element)
                self.assertEqual(bool(painted), highlighted)

    def test_native_styles_selectors_and_snapshots(self):
        calls = []

        def node_style(node):
            calls.append(node.index)
            return graph.NodeDrawing(label=graph.TextLabel("selected"))

        config = graph.RenderSettings(
            layouts=graph.LayoutSettings(impred_steps=2),
            drawing=graph.DrawOptions(
                node_radius=0.2,
                node_fill=graph.Color("#123456"),
                edge_stroke=graph.Stroke(
                    paint=graph.Color("#246813"), thickness=graph.Length.pt(2)
                ),
            ),
            selectors=graph.DrawingSelectors(node=node_style),
        )
        view = self.graph.render(config=config)
        self.assertEqual(calls, [0, 1, 2])
        config.drawing = graph.DrawOptions(node_fill=graph.Color("#abcdef"))
        self.graph.node("a").drawing.label = graph.TextLabel("changed")
        self.graph.add_node(graph.node("before measurement"))
        result = external_metrics_snapshot(view)
        self.assert_native(result)
        self.assertEqual(calls, [0, 1, 2])
        svg = result.to_svg()
        self.assertIn("#123456", svg)
        self.assertIn("#246813", svg)
        self.assertIn("<mtext>selected</mtext>", svg)
        self.assertNotIn("<mtext>changed</mtext>", svg)
        self.assertNotIn("before measurement", svg)
        self.graph.node("a").drawing.label = graph.TextLabel("changed again")
        self.graph.add_node(graph.node("later"))
        self.assertEqual(result.to_svg(), svg)
        self.assertEqual(calls, [0, 1, 2])
        with TemporaryDirectory() as directory:
            for suffix in ("svg", "html"):
                path = result.save(Path(directory) / f"native.{suffix}")
                self.assertGreater(path.stat().st_size, 0)
            for suffix in ("typ", "pdf", "png"):
                with (
                    self.subTest(suffix=suffix),
                    self.assertRaisesRegex(
                        ValueError if suffix == "typ" else RuntimeError,
                        "outlined labels",
                    ),
                ):
                    result.save(Path(directory) / f"native.{suffix}")

    def test_multiple_layout_passes_use_native_geometry(self):
        advanced = graph.RenderSettings(
            layouts=graph.LayoutSettings(impred_steps=1).then(impred_steps=1),
            drawing=graph.DrawOptions(node_fill=graph.Color("#246813")),
        )
        result = external_metrics_snapshot(self.graph.render(config=advanced))
        self.assert_native(result)
        self.assertIn("#246813", result.to_svg())
        # A selection containing only one half of an internal edge must retain
        # its distinct half-edge styles instead of coloring the whole edge.
        selected = self.graph.subgraph(half_edges=[self.graph.edge("ab").source.index])
        with self.assertRaisesRegex(ValueError, "native SVG renderer"):
            selected.render()


if __name__ == "__main__":
    unittest.main()
