import json
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("measure-notebook-cells.py")
HEADER = "import marimo\napp = marimo.App()\n"


class CellTimingTests(unittest.TestCase):
    def collect(self, cells, name="example.py", expect_error=False):
        with tempfile.TemporaryDirectory() as directory:
            notebook = Path(directory) / name
            notebook.parent.mkdir(parents=True, exist_ok=True)
            notebook.write_text(HEADER + textwrap.dedent(cells))
            output = Path(directory) / "result.json"
            process = subprocess.run(
                [sys.executable, SCRIPT, "--notebook", notebook, "--output", output],
                capture_output=True,
                check=False,
                text=True,
            )
            if expect_error:
                self.assertNotEqual(process.returncode, 0)
                self.assertFalse(output.exists())
                return process.stderr
            self.assertEqual(process.returncode, 0, process.stderr)
            return json.loads(output.read_text())

    def test_real_setup_and_awaited_cells_in_dependency_order(self):
        receipt = self.collect(
            """
            with app.setup:
                import asyncio

            @app.cell
            def _(x):
                assert x == 17
                return

            @app.cell
            async def producer():
                await asyncio.sleep(0.015)
                x = 17
                return (x,)

            @app.cell(disabled=True)
            def disabled():
                raise RuntimeError("disabled")
            """
        )
        self.assertEqual(set(receipt["cells"]), {"setup", "uses.x", "named.producer"})
        self.assertGreaterEqual(receipt["cells"]["named.producer"], 0.015)
        self.assertTrue(all(seconds > 0 for seconds in receipt["cells"].values()))

    def test_errors_and_stopped_leaf_cannot_produce_success(self):
        self.assertIn(
            "cell error",
            self.collect(
                "@app.cell\ndef _():\n    raise ValueError('cell error')\n",
                expect_error=True,
            ),
        )
        self.assertIn(
            "did not finish",
            self.collect(
                "@app.cell\ndef _():\n    import marimo as mo\n    mo.stop(True)\n",
                expect_error=True,
            ),
        )

    def test_submitted_defaults_drive_dependent_cell(self):
        receipt = self.collect(
            """
            @app.cell
            def _():
                import marimo as mo
                integration_settings = mo.ui.dictionary({
                    "size": mo.ui.dropdown([1, 4], value=4)
                }).form()
                return (integration_settings,)

            @app.cell
            def _(integration_settings):
                assert integration_settings.value == {"size": 4}
                return
            """,
            "examples/notebooks/feynkit/08_numerical_cff_ltd_marimo.py",
        )
        self.assertEqual(receipt["settings"], {"size": 4})
        self.assertEqual(receipt["scenario"], "submit-integration-defaults")
        self.assertEqual(len(receipt["cells"]), 2)

    def test_identifiers_survive_reordering_and_implementation_edits(self):
        first = "@app.cell\ndef _():\n    x = 1\n    return (x,)\n"
        second = "@app.cell\ndef _(x):\n    x + 1\n    return\n"
        third = "@app.cell\ndef _(x):\n    x * 3\n    return\n"
        baseline = self.collect(first + second + third)
        revised = self.collect(third + first.replace("x = 1", "x = 2") + second)
        self.assertEqual(baseline["labels"], revised["labels"])
        self.assertEqual(set(baseline["cells"]), set(revised["cells"]))
        self.assertEqual(len(baseline["cells"]), 3)

    def test_background_failure_in_layout_cell_fails_capture(self):
        error = self.collect(
            """
            @app.cell
            def _():
                import marimo as mo
                from ipywidgets import Widget as LayoutView
                return LayoutView, mo

            @app.cell
            def _(LayoutView, mo):
                def _animate():
                    LayoutView()
                    raise ValueError("layout worker failed")
                mo.Thread(target=_animate, daemon=True).start()
                return
            """,
            "crates/linnet-py/examples/layout_stream.py",
            expect_error=True,
        )
        self.assertIn("layout worker failed", error)
        self.assertIn("Layout animation thread failed", error)


if __name__ == "__main__":
    unittest.main()
