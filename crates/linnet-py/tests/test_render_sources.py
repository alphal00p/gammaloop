"""Authored physics documents share typed configuration and source snapshots."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import linnet as ln


class RenderSourcesTests(unittest.TestCase):
    def test_configuration_and_imported_modules_are_snapshotted(self):
        with TemporaryDirectory() as directory:
            module_path = Path(directory) / "labels.typ"
            module_path.write_text("#let label = [original label]\n")
            config = ln.RenderConfig(
                title=ln.TypstModule.file(module_path).content("label"),
                layouts=ln.LayoutOptions(seed=13, steps=0),
                template_options={"marker": "original"},
            )
            prepared = ln.PreparedRender.from_sources(
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

    def test_invalid_configuration_and_missing_entrypoint_are_rejected(self):
        with self.assertRaises(TypeError):
            ln.PreparedRender.from_sources({"main.typ": b"hello"}, config={})
        with self.assertRaisesRegex(ValueError, "main.typ"):
            ln.PreparedRender.from_sources({})
        with self.assertRaisesRegex(ValueError, "templates or selectors"):
            ln.PreparedRender.from_sources(
                {"main.typ": b"hello"},
                config=ln.RenderConfig(
                    selectors=ln.DrawingSelectors(edge=lambda edge: None)
                ),
            )


if __name__ == "__main__":
    unittest.main()
