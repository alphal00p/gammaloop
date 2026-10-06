"""All Community graph owners honor the same advanced rendering settings."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import xml.etree.ElementTree as ET

from symbolica import E
from symbolica.community import graph, hepkit, tensor


class CommunityRenderingTests(unittest.TestCase):
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
                details = [
                    json.loads(element.attrib["data-linnet-detail"])
                    for element in ET.fromstring(svg).iter()
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
            self.assertIs(type(edge), graph.Edge)
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
