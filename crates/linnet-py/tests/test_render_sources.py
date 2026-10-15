"""The Python renderer accepts native graphs and typesets only their labels."""

import unittest
import xml.etree.ElementTree as ET
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

    def test_labels_and_titles_are_typeset_without_packages(self):
        graph = ln.build(ln.node("node", label=ln.MathSymbol("p", subscript=1)))
        rendered = graph.render(config=ln.RenderSettings(title=ln.TextLabel("A title")))
        self.assertFalse(hasattr(rendered, "typst_source"))
        root = ET.fromstring(rendered.to_svg())
        self.assertEqual(root.get("data-linnet-renderer"), "native")
        self.assertTrue(any(element.tag.endswith("path") for element in root.iter()))

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
