"""Native graph views expose MathML labels without compiling a document."""

import unittest
import xml.etree.ElementTree as ET
from rendering_metrics import external_metrics_snapshot

from symbolica.community import graph as ln


class RenderSourcesTests(unittest.TestCase):
    def test_full_document_entry_points_are_removed(self):
        self.assertFalse(hasattr(ln.DiagramRender, "from_sources"))
        for name, value in (
            ("template", "custom.typ"),
            ("source_root", "."),
            ("template_options", {}),
        ):
            with (
                self.subTest(name=name),
                self.assertRaisesRegex(TypeError, "unknown RenderSettings option"),
            ):
                ln.RenderSettings(**{name: value})

    def test_labels_and_titles_use_mathml_without_packages(self):
        graph = ln.build(ln.node("node", label=ln.MathSymbol("p", subscript=1)))
        rendered = graph.render(config=ln.RenderSettings(title=ln.TextLabel("A title")))
        self.assertFalse(hasattr(rendered, "typst_source"))
        request = rendered.measurement_request()
        self.assertEqual(len(request["labels"]), 2)
        self.assertIn("A title", request["labels"][0]["mathml"])
        self.assertIn("<msub>", request["labels"][1]["mathml"])
        root = ET.fromstring(external_metrics_snapshot(rendered).to_svg())
        self.assertEqual(root.get("data-linnet-renderer"), "native")
        self.assertEqual(
            len(root.findall(".//{http://www.w3.org/1998/Math/MathML}math")), 2
        )

    def test_executable_style_is_rejected_without_typst_fallback(self):
        graph = ln.build(ln.node("node"))
        config = ln.RenderSettings(
            style=ln.GraphStyleOptions(
                node_style=ln.TypstModule.file("unused.typ").function("style")
            )
        )
        with self.assertRaisesRegex(ValueError, "native SVG renderer"):
            graph.render(config=config)


if __name__ == "__main__":
    unittest.main()
