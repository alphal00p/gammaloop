"""All Community graph owners honor the same advanced rendering settings."""

import json
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from tempfile import TemporaryDirectory

# HepKit diagrams and tensor networks exist only in the Symbolica Community
# wheel. FeynKit's Rust tests run this file without the suite's support module.
try:
    from symbolica import E
except ModuleNotFoundError:
    raise unittest.SkipTest("needs the Symbolica Community wheel") from None
from symbolica.community import graph, hepkit, tensor


class CommunityRenderingTests(unittest.TestCase):
    def test_diagrams_paint_configured_flow_and_momentum_marks(self):
        diagram = (
            hepkit.Model.qed()
            .process(["e-", "e+"], ["a", "a"])
            .generate_diagrams(loops=0, progress=None)
            .diagrams[0]
        )
        flow = graph.Mark("stealth", length="0.3cm", fill=graph.Color("#1a2b3c"))
        momentum = graph.Mark.combine(
            graph.Mark("bar", stroke=graph.Color("#4d5e6f")),
            "1pt",
            graph.Mark("triangle", fill=graph.Color("#7a8b9c")),
        )
        # `Mark` values and their public dictionaries configure the heads of
        # native scenes and of scenes whose other settings need Typst.
        for flow_arrow, momentum_arrow in (
            (flow, momentum),
            (flow.to_dict(), momentum.to_dict()),
        ):
            style = graph.GraphStyleOptions(
                edge_style={
                    "flow-arrow": flow_arrow,
                    "momentum-arrow": momentum_arrow,
                }
            )
            for route, layouts in (
                ("native", graph.INHERIT),
                ("typst", graph.LayoutSettings(impred_steps=1).then(impred_steps=1)),
            ):
                with self.subTest(route=route, data=type(flow_arrow).__name__):
                    settings = graph.RenderSettings(style=style, layouts=layouts)
                    svg = diagram.render(config=settings, momenta=True).to_svg()
                    for paint in ("#1a2b3c", "#4d5e6f", "#7a8b9c"):
                        self.assertIn(paint, svg)

    def test_advanced_styles_and_snapshots_across_owners(self):
        diagram = (
            hepkit.Model.phi3()
            .process(["phi"], ["phi", "phi"])
            .generate_diagrams(loops=0, progress=None)
            .diagrams[0]
        )
        network = tensor.TensorNetwork(E("x+2"))
        for owner in (diagram, network):
            with self.subTest(owner=type(owner).__name__):
                settings = graph.RenderSettings(
                    layouts=graph.LayoutSettings(impred_steps=1).then(impred_steps=1),
                    drawing=graph.DrawOptions(
                        node_fill=graph.Color("#123456"),
                        edge_stroke=graph.Stroke(
                            paint=graph.Color("#246813"), thickness=graph.Length.pt(2)
                        ),
                    ),
                )
                result = owner.render(config=settings)
                self.assertIs(type(result), graph.DiagramRender)
                settings.drawing = graph.DrawOptions(node_fill=graph.Color("#abcdef"))
                svg = result.to_svg()
                self.assertIn("#123456", svg)
                self.assertIn("#246813", svg)
                self.assertNotIn("#abcdef", svg)
                self.assertIn(
                    "prefers-color-scheme:dark"
                    if owner is diagram
                    else "color-scheme:light dark",
                    svg,
                )
                self.assertEqual(result.to_svg(), svg)
                self.assertNotIn("#linnet-", svg)
                interactive = result.to_svg(interactive=True)
                self.assertIn("#123456", interactive)
                details = [
                    json.loads(element.attrib["data-linnet-detail"])
                    for element in ET.fromstring(interactive).iter()
                    if element.get("data-linnet-kind") == "edge"
                ]
                self.assertTrue(details)
                if owner is diagram:
                    self.assertTrue(
                        any(item.get("particle") == "phi" for item in details)
                    )

    def test_all_owners_preserve_multipage_snapshots_and_exports(self):
        diagrams = (
            hepkit.Model.phi3()
            .process(["phi", "phi"], ["phi", "phi"])
            .generate_diagrams(loops=0, progress=None)
            .diagrams
        )
        diagram = diagrams[0]
        generic = diagram.to_graph()
        owners = (
            (generic, 2),
            (generic.full_subgraph(), 2),
            (diagram, 2),
            (tensor.TensorNetwork(E("x+2")), 2),
            (hepkit.Amplitude(diagrams[:2]), 4),
        )
        self.assertFalse(hasattr(hepkit, "AmplitudeRender"))
        with TemporaryDirectory() as directory:
            root = Path(directory)
            template = root / "pages.typ"
            template.write_text(
                "#let render(config) = {\n"
                "  set page(width: 80pt, height: 40pt, margin: 2pt, fill: none)\n"
                "  [First #pagebreak() Second]\n"
                "}\n"
            )
            settings = graph.RenderSettings(template=template)
            results = [
                (owner.render(config=settings), pages) for owner, pages in owners
            ]
            template.unlink()
            for index, (result, count) in enumerate(results):
                with self.subTest(owner=type(owners[index][0]).__name__):
                    self.assertIs(type(result), graph.DiagramRender)
                    self.assertEqual(len(result.to_svg_pages()), count)
                    self.assertEqual(result.to_html().count("<svg"), count)
                    self.assertTrue(result.typst_source)
                    self.assertIsNone(result._repr_svg_())
                    self.assertEqual(result._mime_(), ("text/html", result.to_html()))
                    with self.assertRaisesRegex(ValueError, f"{count} pages"):
                        result.to_svg()
                    output = root / str(index)
                    self.assertEqual(
                        len(result.save_pages(output, format="png")), count
                    )
                    result.save(output / "all.pdf")
                    self.assertTrue(
                        (output / "all.pdf").read_bytes().startswith(b"%PDF")
                    )
            self.assertEqual(len(results[-1][0].diagrams), 2)
            self.assertTrue(
                all(
                    type(child) is graph.DiagramRender
                    for child in results[-1][0].diagrams
                )
            )

    def test_advanced_momentum_arrows_and_live_selector_views(self):
        diagram = (
            hepkit.Model.phi3()
            .process(["phi"], ["phi", "phi"])
            .generate_diagrams(loops=0, progress=None)
            .diagrams[0]
        )
        seen = []

        def select(edge):
            self.assertIs(type(edge), graph.EdgeView)
            self.assertIn("particle", edge.data)
            seen.append(edge.index)
            return graph.EdgeDrawing(
                style={"stroke": graph.Stroke(paint=graph.Color("#123456"))}
            )

        result = diagram.render(
            config=graph.RenderSettings(
                layouts=graph.LayoutSettings(impred_steps=1).then(impred_steps=1),
                selectors=graph.DrawingSelectors(edge=select),
            ),
            style=hepkit.DiagramStyle(show_momentum=True, momentum_arrows=True),
        )
        self.assertEqual(len(seen), len(diagram.edges))
        self.assertIn("#123456", result.to_svg())
        result.to_svg()
        self.assertEqual(len(seen), len(diagram.edges))


if __name__ == "__main__":
    unittest.main()
