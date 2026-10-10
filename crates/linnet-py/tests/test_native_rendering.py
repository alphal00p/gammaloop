"""Ordinary graphs share the native SVG path used by domain drawings."""

import json
import os
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from linnet_under_test import graph


class NativeRenderingTests(unittest.TestCase):
    def setUp(self):
        self.graph = graph.Graph(
            nodes=[
                graph.Node("a", label=graph.TextLabel("left")),
                graph.Node("b", label=graph.MathSymbol("p", subscript=1)),
                graph.Node("isolated"),
            ],
            edges=[
                graph.Edge(source="a", sink="b", name="ab"),
                graph.Edge(source="a", sink="b", name="parallel"),
                graph.Edge(source="b", sink="b", name="loop"),
                graph.Edge(source="a", name="open"),
            ],
        )

    def assert_native(self, result):
        svg = result.to_svg(interactive=True)
        self.assertIn("data-linnet-kind", svg)
        self.assertEqual(ET.fromstring(svg).get("data-linnet-renderer"), "native")
        # Files and documents get the same drawing without inspection links.
        drawing = result.to_svg()
        self.assertNotIn("<a ", drawing)
        self.assertNotIn("data-linnet-kind", drawing)
        self.assertEqual(ET.fromstring(drawing).get("data-linnet-renderer"), "native")
        # Configuration remains inspectable; exports embed the rendered SVG.
        self.assertIn("_linnet_template.render", result.typst_source)
        self.assertIn("image", result.to_linnest())
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
                result = owner.render()
                self.assert_native(result)
                self.assertEqual(result.to_svg(), result.to_svg_pages()[0])
                self.assertEqual(result._mime_(), ("text/html", result.to_html()))
        self.assertIn("#c58b13", self.graph.subgraph(edges=["ab"]).render().to_svg())
        self.assertNotIn("#c58b13", self.graph.empty_subgraph().render().to_svg())

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
        result = self.graph.render(config=config)
        self.assert_native(result)
        self.assertEqual(calls, [0, 1, 2])
        svg = result.to_svg()
        self.assertIn("#123456", svg)
        self.assertIn("#246813", svg)
        config.drawing = graph.DrawOptions(node_fill=graph.Color("#abcdef"))
        self.graph.node("a").drawing.label = graph.TextLabel("changed")
        self.graph.add_node(graph.Node("later"))
        self.assertEqual(result.to_svg(), svg)
        self.assertEqual(calls, [0, 1, 2])
        with TemporaryDirectory() as directory:
            for suffix in ("svg", "html", "typ", "pdf", "png"):
                path = result.save(Path(directory) / f"native.{suffix}")
                self.assertGreater(path.stat().st_size, 0)
            self.assert_static_files(result, Path(directory))

    def assert_static_files(self, result, directory):
        # Only the HTML display is interactive; file exports carry no targets.
        self.assertEqual((directory / "native.svg").read_text(), result.to_svg())
        self.assertNotIn("linnet-edge", (directory / "native.typ").read_text())
        self.assertNotIn(b"/Subtype/Link", (directory / "native.pdf").read_bytes())
        self.assertIn(
            'data-linnet-kind="edge"', (directory / "native.html").read_text()
        )

    def test_graph_edge_style_is_the_native_edge_default(self):
        settings = graph.RenderSettings(
            style=graph.GraphStyleOptions(
                edge_style={"stroke": graph.Stroke(paint=graph.Color("#13579b"))}
            )
        )
        result = self.graph.render(config=settings)
        self.assert_native(result)
        self.assertIn("#13579b", result.to_svg())

    def test_native_renders_and_exports_stage_no_typst_project(self):
        with TemporaryDirectory() as directory:
            # Typst projects are staged in a temporary directory; without one,
            # only the Typst route can fail.
            missing = Path(directory) / "missing"
            with patch.dict(os.environ, {"TMPDIR": str(missing)}):
                for owner in (self.graph, self.graph.subgraph(edges=["ab"])):
                    result = owner.render()
                    self.assert_native(result)
                    for suffix in ("svg", "html", "typ", "pdf", "png"):
                        path = result.save(Path(directory) / f"native.{suffix}")
                        self.assertGreater(path.stat().st_size, 0)
                with self.assertRaises(RuntimeError):
                    graph.DiagramRender.from_sources({"main.typ": b"[typst]"})

    def test_advanced_requests_keep_full_typst_rendering(self):
        advanced = graph.RenderSettings(
            layouts=graph.LayoutSettings(impred_steps=1).then(impred_steps=1),
            drawing=graph.DrawOptions(node_fill=graph.Color("#246813")),
        )
        result = self.graph.render(config=advanced)
        self.assertIn("_linnet_template.render", result.typst_source)
        self.assertIn("#246813", result.to_svg())
        # The Typst route compiles inspection targets only for interactive pages.
        self.assertNotIn("#linnet-", result.to_svg())
        self.assertIn('data-linnet-kind="halfedge"', result.to_svg(interactive=True))
        with TemporaryDirectory() as directory:
            for suffix in ("svg", "html", "typ", "pdf"):
                result.save(Path(directory) / f"native.{suffix}")
            self.assert_static_files(result, Path(directory))
        # A selection containing only one half of an internal edge must retain
        # its distinct half-edge styles instead of coloring the whole edge.
        selected = self.graph.subgraph(half_edges=[self.graph.edge("ab").source.index])
        result = selected.render()
        self.assertIn("_linnet_subgraph.focus", result.typst_source)
        self.assertIn("#c58b13", result.to_svg())


if __name__ == "__main__":
    unittest.main()
