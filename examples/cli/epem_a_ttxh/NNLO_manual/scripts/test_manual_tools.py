#!/usr/bin/env python3
"""Focused tests for the manual NNLO audit utilities."""

from __future__ import annotations

import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

import psutil

SCRIPTS = Path(__file__).resolve().parent


def load_script(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


analyze_approach = load_script("analyze_approach")
collect_results = load_script("collect_results")
normalize_dot = load_script("normalize_dot")


def approach_document(slope: float, *, asymmetric: bool = False) -> dict:
    points = []
    index = 0
    for sign in (-1.0, 1.0):
        branch_slope = slope + (1.0 if asymmetric and sign > 0 else 0.0)
        for exponent in range(2, 27):
            parameter = sign * 10.0 ** (-exponent / 4)
            magnitude = abs(parameter) ** branch_slope
            points.append(
                {
                    "index": index,
                    "axis_index": 0,
                    "t": parameter,
                    "status": "evaluated",
                    "evaluation": {
                        "total_weight": {"re": magnitude, "im": 0.0},
                        "metadata": {"is_nan": False},
                    },
                }
            )
            index += 1
    return {"schema_version": 3, "points": points}


class DotNormalizationTests(unittest.TestCase):
    def test_removes_only_requested_attributes_and_preserves_multiline_toml(self):
        source = """digraph GL001{
    threshold_counterterms = "
schema_version = 1
expression = \\"dir=not-an-attribute\\"
";
    0 -> 1 [id=0, pin="1,2", lmb_rep="k1 + p1, p2", dir=back, lmb_id="0", particle="g"];
}
"""
        output, counts = normalize_dot.normalize_text(source)
        self.assertEqual(counts, {"pin": 1, "dir": 1, "lmb_rep": 1})
        self.assertIn('expression = \\"dir=not-an-attribute\\"', output)
        self.assertIn('lmb_id="0"', output)
        self.assertIn('particle="g"', output)
        self.assertNotRegex(output.split('";\n', 1)[1], normalize_dot.RESIDUAL)

    def test_single_line_embedded_document_does_not_mask_following_edges(self):
        source = """digraph G{
 threshold_counterterms="schema_version = 1\\ncuts = []";
 0 -> 1 [id=0, pin="1,2", dir=none, lmb_id="0"];
}
"""
        output, counts = normalize_dot.normalize_text(source)
        self.assertEqual(counts["pin"], 1)
        self.assertEqual(counts["dir"], 1)
        self.assertIn("threshold_counterterms", output)


class ApproachFitTests(unittest.TestCase):
    def fit(self, document: dict, rank: float = 3.0):
        points = analyze_approach.evaluated_points(document, 0)
        return analyze_approach.fit_branches(points, [8, 12, 16])

    def test_synthetic_passing_and_failing_powers(self):
        passing, errors = self.fit(approach_document(-1.0))
        self.assertFalse(errors)
        self.assertTrue(all(math.isclose(fit.p, 1.0) for fit in passing))
        failing, errors = self.fit(approach_document(-4.0))
        self.assertFalse(errors)
        self.assertTrue(all(math.isclose(fit.p, 4.0) for fit in failing))

    def test_asymmetric_branches_expand_p_envelope(self):
        fits, errors = self.fit(approach_document(-1.0, asymmetric=True))
        self.assertFalse(errors)
        self.assertGreater(max(fit.p for fit in fits) - min(fit.p for fit in fits), 0.5)

    def test_nonfinite_and_missing_points_are_rejected(self):
        document = approach_document(-1.0)
        for point in document["points"]:
            point["evaluation"]["metadata"]["is_nan"] = True
        fits, errors = self.fit(document)
        self.assertFalse(fits)
        self.assertTrue(errors)

    def test_threshold_activity_is_read_from_event_decomposition(self):
        document = approach_document(0.0)
        self.assertFalse(analyze_approach.has_active_threshold_counterterm(document, 0))
        document["points"][0]["evaluation"]["threshold_counterterm_summary"] = {
            "counterterm_sum": {"re": 1.0, "im": 0.0}
        }
        self.assertTrue(analyze_approach.has_active_threshold_counterterm(document, 0))

    def test_threshold_activity_uses_current_additional_contribution_schema(self):
        document = approach_document(0.0)
        document["points"][0]["evaluation"]["additional_contribution_sums"] = {
            "original": {"re": 1.0, "im": 0.0},
            "threshold_counterterm_0": {"re": 0.0, "im": 2.0},
        }
        self.assertTrue(analyze_approach.has_active_threshold_counterterm(document, 0))


class CollectionTests(unittest.TestCase):
    def test_summary_contains_every_graph_and_live_progress(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            state = root / "state_GL297"
            state.mkdir()
            state.joinpath("result_record.json").write_text(
                json.dumps(
                    {
                        "threshold_edited": "yes",
                        "all_limits": "inconclusive",
                        "observation": "pilot | still running",
                        "generation": {"time_seconds": 12.5},
                        "integration": {"converged": False},
                    }
                )
            )
            rendered = collect_results.summary_markdown(root)
            self.assertIn("Progress: **1/71**", rendered)
            self.assertEqual(rendered.count("\n| GL"), 71)
            self.assertIn("pilot \\| still running", rendered)

    def test_collects_internal_timings_and_absolute_convergence(self):
        with tempfile.TemporaryDirectory() as temporary:
            state = Path(temporary)
            (state / "generation_summary.json").write_text(
                json.dumps(
                    {
                        "peak_ram_bytes": 123,
                        "reports": [
                            {"stats": {"total_time": {"secs": 2, "nanos": 500_000_000}}}
                        ],
                    }
                )
            )
            workspace = state / "integration_workspace"
            workspace.mkdir()
            workspace.joinpath("integration_result.json").write_text(
                json.dumps(
                    {
                        "slots": [
                            {
                                "integral": {
                                    "neval": 100,
                                    "result": {"re": 1.0, "im": 2.0},
                                    "error": {"re": 0.1, "im": 0.2},
                                },
                                "integration_statistics": {"nan_percentage": 0.0},
                                "max_weight_info": [],
                                "absolute": {
                                    "integral": {
                                        "neval": 100,
                                        "result": {"re": 3.0, "im": 4.0},
                                        "error": {"re": 0.3, "im": 0.4},
                                    },
                                    "table_results": [
                                        {
                                            "component": "|re|",
                                            "value": 3.0,
                                            "relative_error_percent": 10.0,
                                        },
                                        {
                                            "component": "|im|",
                                            "value": 4.0,
                                            "relative_error_percent": 10.0,
                                        },
                                    ],
                                    "max_weight_info": [],
                                },
                            }
                        ]
                    }
                )
            )
            log = state / "integration.log"
            log.write_text("last update: 42.5 µs /sample/core\n")
            generation = collect_results.generation_data(state)
            integration = collect_results.integration_data(state, log)
            self.assertEqual(generation["time_seconds"], 2.5)
            self.assertEqual(integration["runtime_per_sample_per_core"], "42.5µs")
            self.assertTrue(integration["converged"])


class GuardTests(unittest.TestCase):
    def test_timeout_cleans_nested_process_tree(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            record = root / "record.json"
            child_code = (
                "import subprocess,time; "
                f"subprocess.Popen([{sys.executable!r}, '-c', "
                "'import time; time.sleep(60)']); time.sleep(60)"
            )
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPTS / "run_guarded.py"),
                    "--job",
                    "selftest",
                    "--log",
                    str(root / "guard.log"),
                    "--record",
                    str(record),
                    "--disk-path",
                    str(root),
                    "--disk-floor-gib",
                    "0",
                    "--min-available-ram-gib",
                    "0",
                    "--timeout-seconds",
                    "0.5",
                    "--interrupt-grace-seconds",
                    "0.5",
                    "--terminate-grace-seconds",
                    "0.5",
                    "--timeout-success",
                    "--",
                    sys.executable,
                    "-c",
                    child_code,
                ],
                check=False,
                timeout=8,
            )
            self.assertEqual(completed.returncode, 0)
            data = json.loads(record.read_text())
            self.assertEqual(data["exit_reason"], "timeout")
            time.sleep(0.2)
            for identity in data["known_processes"]:
                try:
                    process = psutil.Process(identity["pid"])
                except psutil.NoSuchProcess:
                    continue
                self.assertNotAlmostEqual(
                    process.create_time(), identity["create_time"], places=3
                )


if __name__ == "__main__":
    unittest.main()
