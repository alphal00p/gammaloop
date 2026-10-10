import copy
import json
import math
import subprocess
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path

from performance_gate import PerformanceGate

BASELINE = "a" * 40
CANDIDATE = "b" * 40


def measurements(generation=10.0, warm=10e-6, *, confirmation=False):
    def phase(run_id, first):
        return {
            "run_id": run_id,
            "pairs": [
                {
                    "order": first
                    if index % 2 == 0
                    else "BA"
                    if first == "AB"
                    else "AB",
                    "baseline": {
                        "generation_seconds": 10.0,
                        "warm_runtime_seconds_per_sample": 10e-6,
                    },
                    "candidate": {
                        "generation_seconds": generation,
                        "warm_runtime_seconds_per_sample": warm,
                    },
                    "baseline_output": [1.25, -0.4, 0.1, 0.0],
                    "candidate_output": [1.25, -0.4, 0.1, 0.0],
                }
                for index in range(5)
            ],
        }

    identity = {
        "fixture_sha256": "1" * 64,
        "corpus_sha256": "2" * 64,
        "config_sha256": "3" * 64,
    }
    data = {
        "schema_version": 1,
        "baseline": {"commit": BASELINE, **identity},
        "candidate": {"commit": CANDIDATE, **identity},
        "required_metrics": ["generation_seconds", "warm_runtime_seconds_per_sample"],
        "correctness_tolerance": {"absolute": 1e-12, "relative": 1e-9},
        "discovery": phase(1, "AB"),
    }
    if confirmation:
        data["confirmation"] = phase(2, "BA")
    return data


class PerformanceGateTests(unittest.TestCase):
    def setUp(self):
        self.gate = PerformanceGate(BASELINE, CANDIDATE)

    def test_clean_measurement_and_speedup_pass_without_confirmation(self):
        for data in (measurements(), measurements(8.0, 7e-6)):
            with self.subTest(data=data):
                self.assertEqual(self.gate.analyze(data)["status"], 0)

    def test_regression_needs_independent_confirmation(self):
        report = self.gate.analyze(measurements(13.0, 13e-6))
        self.assertEqual(report["status"], 2)
        self.assertEqual(report["confirmation_required"], 1)
        self.assertEqual(report["invalid_input"], 0)
        confirmed = self.gate.analyze(measurements(13.0, 13e-6, confirmation=True))
        self.assertEqual(confirmed["status"], 1)

    def test_all_three_metric_thresholds_have_relative_and_absolute_floors(self):
        # Large percentage changes on tiny operations must not fail the gate.
        for metric, before, after in (
            ("generation_seconds", 0.1, 0.2),
            ("load_seconds", 0.01, 0.02),
            ("warm_runtime_seconds_per_sample", 1e-7, 2e-7),
        ):
            data = measurements(confirmation=True)
            data["required_metrics"] = [metric]
            for phase in (data["discovery"], data["confirmation"]):
                for pair in phase["pairs"]:
                    pair["baseline"] = {metric: before}
                    pair["candidate"] = {metric: after}
            with self.subTest(metric=metric):
                self.assertEqual(self.gate.analyze(data)["status"], 0)
        data = measurements(confirmation=True)
        data["required_metrics"] = ["load_seconds"]
        for phase in (data["discovery"], data["confirmation"]):
            for pair in phase["pairs"]:
                pair["baseline"] = {"load_seconds": 1.0}
                pair["candidate"] = {"load_seconds": 1.3}
        self.assertEqual(self.gate.analyze(data)["status"], 1)

    def test_large_absolute_delta_below_relative_threshold_passes(self):
        self.assertEqual(self.gate.analyze(measurements(11.0))["status"], 0)

    def test_single_slow_outlier_does_not_gate(self):
        data = measurements()
        data["discovery"]["pairs"][-1]["candidate"]["generation_seconds"] = 100.0
        self.assertEqual(self.gate.analyze(data)["status"], 0)

    def test_noisy_measurements_do_not_silently_pass(self):
        data = measurements(confirmation=True)
        for phase in (data["discovery"], data["confirmation"]):
            for pair, value in zip(phase["pairs"], [7.0, 15.0, 10.0, 15.0, 7.0]):
                pair["candidate"]["generation_seconds"] = value
        report = self.gate.analyze(data)
        self.assertEqual(report["status"], 2)
        self.assertEqual(report["invalid_input"], 0)

    def test_uncertainty_crossing_threshold_is_inconclusive(self):
        data = measurements()
        for pair, value in zip(
            data["discovery"]["pairs"], [11.9, 11.8, 12.0, 11.8, 12.1]
        ):
            pair["candidate"]["generation_seconds"] = value
        self.assertEqual(self.gate.analyze(data)["status"], 2)

    def test_at_least_four_pairs_must_be_slower(self):
        data = measurements(confirmation=True)
        for phase in (data["discovery"], data["confirmation"]):
            for pair, value in zip(phase["pairs"], [9.0, 9.0, 12.5, 12.6, 12.7]):
                pair["candidate"]["generation_seconds"] = value
        self.assertEqual(self.gate.analyze(data)["status"], 2)

    def test_confirmation_must_reproduce_discovery(self):
        data = measurements(13.0, confirmation=True)
        for pair in data["confirmation"]["pairs"]:
            pair["candidate"]["generation_seconds"] = 10.0
        self.assertEqual(self.gate.analyze(data)["status"], 2)

    def test_incomplete_phase_and_missing_metrics_are_invalid(self):
        cases = []
        data = measurements()
        data["discovery"]["pairs"].pop()
        cases.append(data)
        data = measurements()
        del data["discovery"]["pairs"][2]["candidate"]["generation_seconds"]
        cases.append(data)
        data = measurements()
        data["required_metrics"] = []
        cases.append(data)
        data = measurements()
        data["required_metrics"].append("generation_seconds")
        cases.append(data)
        for data in cases:
            with self.subTest(data=data):
                report = self.gate.analyze(data)
                self.assertEqual(report["status"], 2)
                self.assertEqual(report["invalid_input"], 1)

    def test_invalid_nonfinite_and_overflowing_timings_are_inconclusive(self):
        for value in (math.nan, math.inf, -math.inf, 0.0, -1.0, True, "1.0", 10**400):
            data = measurements()
            data["discovery"]["pairs"][0]["candidate"]["generation_seconds"] = value
            with self.subTest(value=value):
                report = self.gate.analyze(data)
                self.assertEqual(report["status"], 2)
                self.assertEqual(report["invalid_input"], 1)
                json.dumps(report, allow_nan=False)
        data = measurements()
        pair = data["discovery"]["pairs"][0]
        pair["baseline"]["generation_seconds"] = 1e-308
        pair["candidate"]["generation_seconds"] = 1e308
        self.assertEqual(self.gate.analyze(data)["invalid_input"], 1)

    def test_correctness_mismatch_fails_without_timing_confirmation(self):
        data = measurements(8.0, 7e-6)
        data["discovery"]["pairs"][0]["candidate_output"][0] = 1.3
        report = self.gate.analyze(data)
        self.assertEqual(report["status"], 1)
        self.assertEqual(report["correctness_mismatches"], 1)

    def test_numerical_tolerance_allows_roundoff(self):
        data = measurements()
        data["discovery"]["pairs"][0]["candidate_output"][0] += 1e-11
        self.assertEqual(self.gate.analyze(data)["status"], 0)

    def test_scalar_tolerance_detects_changes_smaller_than_physical_absolute_floor(
        self,
    ):
        data = measurements()
        manifest = (
            Path(__file__).resolve().parents[2] / "benchmarks/performance/suite.toml"
        )
        data["correctness_tolerance"] = tomllib.loads(manifest.read_text())["scalar"][
            "correctness_tolerance"
        ]
        for pair in data["discovery"]["pairs"]:
            pair["baseline_output"] = [8e-12, 4e-14]
            pair["candidate_output"] = [8.001e-12, 4e-14]
        self.assertEqual(self.gate.analyze(data)["status"], 1)
        for pair in data["discovery"]["pairs"]:
            pair["candidate_output"] = [8e-12 * (1 + 1e-11), 4e-14]
        self.assertEqual(self.gate.analyze(data)["status"], 0)

    def test_zero_physical_reference_is_invalid_and_zero_candidate_fails(self):
        data = measurements()
        pair = data["discovery"]["pairs"][0]
        pair["baseline_output"] = [0.0] * 4
        pair["candidate_output"] = [0.0] * 4
        report = self.gate.analyze(data)
        self.assertEqual(report["status"], 2)
        self.assertEqual(report["invalid_input"], 1)
        pair["baseline_output"] = [1e-15, 0.0, 0.0, 0.0]
        report = self.gate.analyze(data)
        self.assertEqual(report["status"], 1)
        self.assertEqual(report["correctness_mismatches"], 1)

    def test_invalid_correctness_data_does_not_pass(self):
        for value in ([], [math.nan], [1.0], "bad"):
            data = measurements()
            data["discovery"]["pairs"][0]["candidate_output"] = value
            with self.subTest(value=value):
                self.assertEqual(self.gate.analyze(data)["invalid_input"], 1)
        data = measurements()
        data["correctness_tolerance"]["absolute"] = -1.0
        self.assertEqual(self.gate.analyze(data)["invalid_input"], 1)

    def test_pin_and_shared_workload_identity_must_match(self):
        for key, value in (("commit", "c" * 40), ("corpus_sha256", "4" * 64)):
            data = measurements()
            data["candidate"][key] = value
            with self.subTest(key=key):
                report = self.gate.analyze(data)
                self.assertEqual(report["status"], 2)
                self.assertEqual(report["identity_mismatch"], 1)

    def test_pair_order_and_confirmation_independence_are_required(self):
        data = measurements()
        data["discovery"]["pairs"][1]["order"] = "AB"
        self.assertEqual(self.gate.analyze(data)["invalid_input"], 1)
        data = measurements(confirmation=True)
        data["confirmation"]["run_id"] = 1
        self.assertEqual(self.gate.analyze(data)["invalid_input"], 1)
        data = measurements(confirmation=True)
        data["confirmation"] = copy.deepcopy(data["discovery"])
        data["confirmation"]["run_id"] = 2
        self.assertEqual(self.gate.analyze(data)["invalid_input"], 1)

    def test_unknown_schema_or_fields_are_invalid(self):
        data = measurements()
        data["schema_version"] = 2
        self.assertEqual(self.gate.analyze(data)["invalid_input"], 1)
        data = measurements()
        data["debug_output"] = "This must never reach the report"
        self.assertEqual(self.gate.analyze(data)["invalid_input"], 1)

    def test_report_values_are_numeric_only(self):
        def check(value):
            if isinstance(value, dict):
                for child in value.values():
                    check(child)
            else:
                self.assertIn(type(value), (int, float))
                self.assertTrue(math.isfinite(value))

        check(self.gate.analyze(measurements(13.0, confirmation=True)))

    def test_cli_exit_codes_and_durable_report(self):
        script = Path(__file__).with_name("performance_gate.py")
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            source, output = folder / "input.json", folder / "reports" / "report.json"
            for data, expected in (
                (measurements(), 0),
                (measurements(13.0, confirmation=True), 1),
                (measurements(13.0), 2),
            ):
                source.write_text(json.dumps(data))
                process = subprocess.run(
                    [
                        sys.executable,
                        str(script),
                        str(source),
                        "--baseline-commit",
                        BASELINE,
                        "--candidate-commit",
                        CANDIDATE,
                        "--report",
                        str(output),
                    ],
                    capture_output=True,
                    text=True,
                    check=False,
                )
                with self.subTest(expected=expected):
                    self.assertEqual(process.returncode, expected, process.stderr)
                    self.assertEqual(json.loads(output.read_text())["status"], expected)
                    self.assertNotIn(BASELINE, output.read_text())
                    self.assertIn("Performance gate:", process.stderr)
                    self.assertIn("generation_seconds discovery: 10 ->", process.stderr)
                    self.assertIn("policy >20% and >0.5", process.stderr)
                    self.assertIn("paired relative noise", process.stderr)
            source.write_text("invalid json")
            process = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    str(source),
                    "--baseline-commit",
                    BASELINE,
                    "--candidate-commit",
                    CANDIDATE,
                    "--report",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(process.returncode, 2)
            self.assertEqual(json.loads(output.read_text())["invalid_input"], 1)


if __name__ == "__main__":
    unittest.main()
