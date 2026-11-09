"""A missing upstream example must fail even when local discovery is consistent."""

import importlib.util
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

spec = importlib.util.spec_from_file_location(
    "notebook_export_coverage",
    Path(__file__).resolve().parents[1] / "examples/export_wasm.py",
)
exporter = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = exporter
spec.loader.exec_module(exporter)


class CanonicalCoverageTest(unittest.TestCase):
    def test_four_loop_omission_is_not_hidden_by_local_discovery(self):
        # Reproduce the actual omission: neither source nor registry entry exists.
        # Merely comparing registered notebooks with local Marimo files passes.
        notebooks = tuple(
            notebook
            for notebook in exporter.NOTEBOOKS
            if notebook.ready_value != "05_four_loop_numerator_marimo"
        )
        sources = [(notebook, notebook.source.resolve()) for notebook in notebooks]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            examples = root / "crates/linnet-py/examples"
            examples.mkdir(parents=True)
            for notebook, source in sources:
                target = examples / notebook.filename
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)
            with (
                patch.object(exporter, "EXAMPLES_DIR", examples),
                patch.object(exporter, "NOTEBOOKS", notebooks),
                self.assertRaisesRegex(ValueError, "Missing canonical Community"),
            ):
                exporter.main(["--output", str(root / "export")])


if __name__ == "__main__":
    unittest.main()
