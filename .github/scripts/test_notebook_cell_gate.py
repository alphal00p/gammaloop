import io
import json
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import notebook_cell_gate as gate

NOTEBOOK = gate.NOTEBOOKS / "chapter/example.py"
APP = "from marimo import App as Notebook\napp = Notebook()\n"


class NotebookGateTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.roots = {s: self.directory / s for s in gate.SIDES}
        self.output = self.directory / "result.json"
        self.summary = self.directory / "summary.md"
        for root in self.roots.values():
            path = root / NOTEBOOK
            path.parent.mkdir(parents=True)
            path.write_text(APP)

    def invoke(self, collect, *, bootstrap=False):
        args = ["gate", "--output", str(self.output), "--summary", str(self.summary)]
        for side, root in self.roots.items():
            args.extend([f"--{side}-root", str(root)])
            if side == "candidate" or not bootstrap:
                args.extend([f"--{side}-python", sys.executable])
        with (
            patch.object(sys, "argv", args),
            patch.object(gate, "collect", side_effect=collect),
            redirect_stdout(io.StringIO()),
            redirect_stderr(io.StringIO()),
        ):
            code = gate.main()
        return code, json.loads(self.output.read_text())

    @staticmethod
    def sample(cells):
        return {
            "cells": cells,
            "labels": {key: key.title() for key in cells},
            "scenario": "default",
            "settings": {},
        }

    def test_discovery_and_independent_confirmation(self):
        calls = []

        def capture(root, python, path):
            calls.append(root.name)
            candidate = root.name == "candidate"
            cells = {"work": 1.5 if candidate else 1.0, "steady": 1.0}
            cells["new" if candidate else "retired"] = 0.01
            return self.sample(cells)

        code, report = self.invoke(capture)
        entry = report["notebooks"][str(NOTEBOOK)]
        self.assertEqual((code, report["status"]), (1, "regression"))
        self.assertEqual(
            calls,
            list(
                ("baseline", "candidate", "candidate", "baseline") * 2
                + ("baseline", "candidate", "candidate", "baseline")
                + ("baseline", "candidate", "candidate", "baseline") * 2
            ),
        )
        self.assertEqual(entry["cells"]["work"]["confirmation"]["slower_pairs"], 5)
        self.assertEqual(
            {key: record["status"] for key, record in entry["cells"].items()},
            {
                "work": "regression",
                "steady": "passed",
                "new": "added",
                "retired": "removed",
            },
        )
        self.assertIn("Work</code> — regression", self.summary.read_text())
        self.assertEqual(len(entry["phases"]["discovery"]), 5)

    def test_tiny_noisy_cells_pass_without_confirmation(self):
        calls = []

        def capture(root, *_):
            calls.append(root.name)
            duration = 0.0001 if root.name == "baseline" else len(calls) * 0.0001
            return self.sample({"tiny": duration})

        code, report = self.invoke(capture)
        self.assertEqual((code, report["status"], len(calls)), (0, "passed", 10))

    def test_nonreplicating_slowdown_is_inconclusive(self):
        calls = []

        def capture(root, *_):
            calls.append(root.name)
            slow = root.name == "candidate" and len(calls) <= 10
            return self.sample({"work": 2.0 if slow else 1.0})

        code, report = self.invoke(capture)
        self.assertEqual((code, report["status"], len(calls)), (2, "inconclusive", 20))

    def test_noisy_discovery_gets_one_independent_confirmation(self):
        noisy = [0.7, 1.3, 1.0, 1.1, 0.9]
        for confirmation, expected in [([1.0] * 5, "passed"), (noisy, "inconclusive")]:
            with self.subTest(expected=expected):
                candidate = iter(noisy + confirmation)

                def capture(root, *_, values=candidate):
                    elapsed = next(values) if root.name == "candidate" else 1.0
                    return self.sample({"work": elapsed})

                code, report = self.invoke(capture)
                cell = report["notebooks"][str(NOTEBOOK)]["cells"]["work"]
                self.assertEqual(code, 0 if expected == "passed" else 2)
                self.assertEqual(
                    (cell["discovery"]["status"], cell["status"]), (2, expected)
                )
                self.assertEqual(cell["confirmation"]["pair_count"], 5)
                if expected == "passed":
                    self.assertEqual(cell["combined"]["pair_count"], 10)

    def test_mismatched_form_settings_do_not_produce_comparisons(self):
        def capture(root, *_):
            result = self.sample({"work": 0.5})
            result["settings"]["samples"] = 100 if root.name == "baseline" else 200
            return result

        code, report = self.invoke(capture)
        self.assertEqual((code, report["status"]), (2, "error"))
        entry = report["notebooks"][str(NOTEBOOK)]
        self.assertIn("different scenarios or settings", entry["error"])
        self.assertEqual(entry["scenarios"]["baseline"]["settings"], {"samples": 100})
        self.assertEqual(entry["cells"], {})

    def test_collection_error_and_invalid_coverage_fail(self):
        for failure in [
            subprocess.CalledProcessError(1, ["collector"]),
            self.sample({}),
            self.sample({"cell": 0.0}),
            self.sample({"cell": float("nan")}),
        ]:
            with self.subTest(failure=failure):

                def capture(*_, failure=failure):
                    if isinstance(failure, Exception):
                        raise failure
                    return failure

                code, report = self.invoke(capture)
                self.assertEqual((code, report["status"]), (2, "error"))
        calls = []

        def change(root, *_):
            calls.append(root.name)
            key = "wrong" if calls.count(root.name) > 1 else "expected"
            return self.sample({key: 0.1})

        code, report = self.invoke(change)
        self.assertEqual(code, 2)
        self.assertIn("coverage changed", report["notebooks"][str(NOTEBOOK)]["error"])

    def test_discovery_finds_nested_apps_and_excludes_helpers(self):
        root = self.roots["candidate"]
        helper = root / gate.NOTEBOOKS / "helper.py"
        helper.write_text('# marimo.App()\n"marimo.App()"\n')
        other = root / gate.NOTEBOOKS / "yet_another.py"
        other.write_text("import marimo as mo\napp = mo.App()\n")
        crate = root / "crates/widget/examples/showcase.py"
        crate.parent.mkdir(parents=True)
        crate.write_text(APP)
        self.assertEqual(
            gate.notebooks(root),
            {NOTEBOOK, other.relative_to(root), crate.relative_to(root)},
        )

    def test_full_suite_requires_candidate_notebooks(self):
        (self.roots["candidate"] / NOTEBOOK).unlink()
        code, report = self.invoke(lambda *_: self.fail("Unexpected measurement"))
        self.assertEqual((code, report["status"]), (2, "error"))
        self.assertIn("no candidate Marimo notebooks", report["error"])

    def test_unmatched_notebook_is_observed_without_claiming_a_comparison(self):
        (self.roots["baseline"] / NOTEBOOK).unlink()
        calls = []

        def capture(root, *_):
            calls.append(root.name)
            return self.sample({"cell": 1.0})

        code, report = self.invoke(capture)
        self.assertEqual(
            (code, report["status"], calls), (0, "uncompared", ["candidate"])
        )
        self.assertEqual(report["notebooks"][str(NOTEBOOK)]["status"], "added")

    def test_bootstrap_requires_host_introduction_and_executes_candidate(self):
        manifest = gate.NOTEBOOKS / "symbolica-host/Cargo.toml"
        host = self.roots["candidate"] / manifest
        host.parent.mkdir()
        host.touch()
        calls = []

        def capture(root, *_):
            calls.append(root.name)
            return self.sample({"cell": 0.1})

        code, report = self.invoke(capture, bootstrap=True)
        self.assertEqual(
            (code, report["status"], calls), (0, "baseline-unavailable", ["candidate"])
        )
        self.assertEqual(report["baseline_inventory"], [str(NOTEBOOK)])
        self.assertEqual(
            report["notebooks"][str(NOTEBOOK)]["phases"]["observation"],
            {"candidate": {"cell": 0.1}},
        )
        code, report = self.invoke(lambda *_: self.sample({}), bootstrap=True)
        self.assertEqual((code, report["status"]), (2, "error"))

        baseline_manifest = self.roots["baseline"] / manifest
        baseline_manifest.parent.mkdir()
        baseline_manifest.touch()
        code, report = self.invoke(capture, bootstrap=True)
        self.assertEqual((code, report["status"]), (2, "error"))
        self.assertIn("--baseline-python", report["error"])
        self.assertEqual(calls, ["candidate"])


if __name__ == "__main__":
    unittest.main()
