"""Authored physics documents share typed configuration and source snapshots."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from symbolica.community import graph as ln


class RenderSourcesTests(unittest.TestCase):
    def test_default_graph_render_uses_native_impred(self):
        left, right = ln.node("left"), ln.node("right")
        graph = ln.build(
            left,
            right,
            ln.edge(ln.source(left), "first", ln.sink(right)),
            ln.edge(ln.source(left), "second", ln.sink(right)),
        )
        explicit = ln.RenderSettings(
            layouts=ln.LayoutSettings(algorithm=ln.LayoutAlgorithm.Impred)
        )
        rendered = graph.to_svg()
        self.assertIn("<svg", rendered)
        self.assertEqual(rendered, graph.to_svg(config=explicit))

    def test_configuration_and_imported_modules_are_snapshotted(self):
        with TemporaryDirectory() as directory:
            module_path = Path(directory) / "labels.typ"
            module_path.write_text("#let label = [original label]\n")
            config = ln.RenderSettings(
                title=ln.TypstModule.file(module_path).content("label"),
                layouts=ln.LayoutSettings(seed=13, steps=0),
                template_options={"marker": "original"},
            )
            prepared = ln.DiagramRender.from_sources(
                {
                    "main.typ": b"""#set page(width: auto, height: auto, fill: none)
#assert.eq(_linnet_config.options.marker, "original")
#assert.eq(_linnet_config.layouts.at(0).seed, 13)
#_linnet_config.title
"""
                },
                config=config,
            )
            config.template_options = {"marker": "changed"}
            module_path.unlink()
            self.assertIn("original", prepared.typst_source)
            self.assertIn("<svg", prepared.to_svg())

    def test_half_edge_routes_validate_and_reach_native_rendering(self):
        points = [(0.5, 1.0), (1.5, 1.0)]
        drawing = ln.HalfEdgeDrawing(route_points=points)
        self.assertEqual(drawing.route_points, tuple(points))
        points.append((9, 9))
        self.assertEqual(len(drawing.route_points), 2)
        for invalid in ["bad", [1], [(1, 2, 3)], [(float("nan"), 0)]]:
            with (
                self.subTest(invalid=invalid),
                self.assertRaises((TypeError, ValueError)),
            ):
                ln.HalfEdgeDrawing(route_points=invalid)

        left, right = ln.node("left"), ln.node("right")
        graph = ln.build(
            left,
            right,
            ln.edge(ln.source(left), "edge", ln.sink(right)),
        )
        config = ln.RenderSettings(
            layouts=ln.LayoutSettings(algorithm=ln.LayoutAlgorithm.Force, steps=0),
            selectors=ln.DrawingSelectors(
                source=lambda half_edge: ln.HalfEdgeDrawing(
                    route_points=[(0.5, 1), (1.5, 1)]
                ),
                sink=lambda half_edge: ln.HalfEdgeDrawing(route_points=[(3, -1)]),
            ),
        )
        prepared = graph.render(config=config)
        self.assertIn("route-points", prepared.typst_source)
        svg = prepared.to_svg()
        self.assertIn("<svg", svg)
        self.assertNotEqual(
            svg, graph.to_svg(config=ln.RenderSettings(layouts=config.layouts))
        )

    def test_multipage_snapshot_and_exports(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            result = ln.DiagramRender.from_sources(
                {"main.typ": b"one #pagebreak() two"}
            )
            self.assertEqual(len(result.to_svg_pages()), 2)
            with self.assertRaisesRegex(ValueError, "2 pages"):
                result.to_svg()
            self.assertIsNone(result._repr_svg_())
            self.assertEqual(len(result.save_pages(root / "svg")), 2)
            self.assertTrue(
                all(
                    path.read_bytes().startswith(b"<svg")
                    for path in (root / "svg").iterdir()
                )
            )
            self.assertEqual(len(result.save_pages(root / "png", format="png")), 2)
            result.save(root / "portable.typ")
            portable = ln.DiagramRender.from_sources(
                {"main.typ": (root / "portable.typ").read_bytes()}
            )
            self.assertEqual(len(portable.to_svg_pages()), 2)
            result.save(root / "all.pdf")
            self.assertTrue((root / "all.pdf").read_bytes().startswith(b"%PDF"))

    def test_invalid_configuration_and_missing_entrypoint_are_rejected(self):
        with self.assertRaises(TypeError):
            ln.DiagramRender.from_sources({"main.typ": b"hello"}, config={})
        with self.assertRaisesRegex(ValueError, "main.typ"):
            ln.DiagramRender.from_sources({})
        with self.assertRaisesRegex(ValueError, "templates or selectors"):
            ln.DiagramRender.from_sources(
                {"main.typ": b"hello"},
                config=ln.RenderSettings(
                    selectors=ln.DrawingSelectors(edge=lambda edge: None)
                ),
            )


if __name__ == "__main__":
    unittest.main()
