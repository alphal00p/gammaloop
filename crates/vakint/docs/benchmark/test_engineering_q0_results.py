"""Synthetic receipt contracts; no physics, cache payload or executable is opened."""

import copy
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import merge_physical_capture_results as adapter


class EngineeringQ0Contracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="synthetic-q0-engineering-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.observations = {}
        for name in ("before", "after_off", "after_on"):
            source = ("a" if name == "before" else "b") * 40
            report = {
                "schema": "vakint.qcd-ladder.benchmark.v1",
                "source_commit": source,
                "compiled_cohort": {"main_revision": source},
                "reduction_completed": True,
                "retained_result_available": True,
                "performance_claim": False,
                "projection": "zero-momentum-metric",
                "input_mode": "prepared",
                "fresh_native_preparation": False,
                "preparation_seconds": None,
                "input_sha256": "c" * 64,
                "result_sha256": ("d" if name == "before" else "e") * 64,
                "result_bytes": 12,
                "guards": {"count": 3},
                "families": {"test": "synthetic"},
                "cache_identity": {"test": True},
                "cache_identity_after": {"test": True},
                "catalogue_closure": {
                    "closed": True,
                    "unknown_count": 0,
                    "master_count": 1,
                },
                "stages": {
                    "saved_rule_application": [1.0],
                    "retained_catalogue_exact_combination": [0.25],
                },
                "reduction_seconds": 10.0 if name == "before" else 2.0,
                "end_to_end_seconds": 12.0 if name == "before" else 3.0,
                "end_to_end_scope": "verified prepared-input loading through cold result serialization; original preparation excluded",
                "cache_hit_diagnostic_enabled": name == "after_on",
            }
            if name == "before":
                report["stages"]["retained_outer_exact_combination"] = [0.5]
            execution = {
                "schema": "vakint.final-instrumented-q0-replay-execution.v1",
                "status": "PASS",
                "source_head": source,
                "collector_exit_code": 0,
                "collector_wait_is_authoritative": True,
                "errors": [],
                "signals_sent": [],
                "elapsed_seconds": report["end_to_end_seconds"] + 1,
            }
            verified = {
                "schema": "vakint.final-instrumented-q0-verification.v1"
                if name == "before"
                else "vakint.coefficient-local-q0-verification.v1",
                "status": "PASS",
                "collector_source_head": source,
                "whole_ordered_guard_Atoms_equal": True,
                "families_equal": True,
                "cache_identity_unchanged": True,
                "guard_count": 3,
                "whole_canonical_input_sha256": report["input_sha256"],
                "whole_canonical_result_sha256": report["result_sha256"],
                "cache_identity": report["cache_identity"],
                "measured_public_reduction_seconds": report["reduction_seconds"],
                "pure_ibp_seconds": None,
                "HEAVY_file_performance_claim": False,
                "numerical_evaluations": 0,
                "master_values_evaluated": 0,
                "IBP_runs_in_verifier": 0,
                "stages": {
                    phase: {"samples": len(values), "total_seconds": sum(values)}
                    for phase, values in report["stages"].items()
                },
                "whole_retained_Atom_equal": True,
                "whole_native_baseline_equal": True,
                "whole_raw_result_equal": True,
                "full58_scalar_coefficient_parity": True,
                "original_native_baseline_authenticated": True,
                "Float_numeric_atoms_rejected": True,
                "master_identity_union_count": 1,
                "full_scalar_coefficient_count": 2,
            }
            records = {
                "benchmark": report,
                "execution": execution,
                "verification": verified,
            }
            if name != "before":
                proof = {
                    "schema": "vakint.complete-master-product-parity.v1",
                    "status": "PASS",
                    "master_count": 1,
                    "product_union_count": 2,
                    "full_scalar_proofs": [
                        {"product": product, "exact_difference_zero": True}
                        for product in ("SYNTHETIC_A", "SYNTHETIC_B")
                    ],
                    "global_together_or_full_amplitude_difference": False,
                    "unit_mass_or_coupling_replacement": False,
                    "numerical_sampling": False,
                }
                records["coefficient_proof"] = proof
                verified["full58_proof"] = self.write(f"{name}/proof.json", proof)
            self.observations[name] = records

    def write(self, name, value):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
        return adapter.evidence(path)

    def index(self):
        observations = copy.deepcopy(self.observations)
        for name, records in observations.items():
            if name != "before":
                records["verification"]["full58_proof"] = self.write(
                    f"{name}/proof.json", records["coefficient_proof"]
                )
        self.write(
            "index.json",
            {
                "schema": "vakint.q0-engineering-evidence-index.v1",
                "observations": {
                    name: {
                        kind: self.write(f"{name}/{kind}.json", value)
                        for kind, value in records.items()
                    }
                    for name, records in observations.items()
                },
            },
        )
        return self.root / "index.json"

    def test_scope_and_timing_boundaries_remain_separate(self):
        run = adapter.engineering_comparison(self.index())
        self.assertFalse(run["performance_eligible"])
        self.assertFalse(run["matched_serialized_repetitions"])
        self.assertFalse(run["cache_only_speedup_claim"])
        self.assertEqual(
            run["observed_before_over_after_off"],
            {
                "public_reduction": 5.0,
                "prepared_input_to_serialization": 4.0,
            },
        )
        self.assertEqual(
            run["observations"]["after_off"]["supervised_wall_seconds"], 4.0
        )
        self.assertIsNone(run["observations"]["after_off"]["pure_ibp_seconds"])

    def test_rejects_incomplete_or_duplicate_product_proof(self):
        proof = self.observations["after_off"]["coefficient_proof"]
        original = copy.deepcopy(proof["full_scalar_proofs"])
        for incomplete in (True, False):
            with self.subTest(incomplete=incomplete):
                rows = copy.deepcopy(original)
                if incomplete:
                    rows.pop()
                else:
                    rows[1]["product"] = rows[0]["product"]
                proof["full_scalar_proofs"] = rows
                with self.assertRaisesRegex(ValueError, "coefficient proof"):
                    adapter.engineering_comparison(self.index())

    def test_rejects_different_diagnostic_result(self):
        self.observations["after_on"]["benchmark"]["result_bytes"] += 1
        with self.assertRaisesRegex(ValueError, "OFF/ON results"):
            adapter.engineering_comparison(self.index())

    def test_rejects_phase_drift(self):
        self.observations["after_off"]["benchmark"]["stages"][
            "saved_rule_application"
        ] = [2.0]
        with self.assertRaisesRegex(ValueError, "phase differs"):
            adapter.engineering_comparison(self.index())

    def test_rejects_missing_displayed_phase(self):
        for kind in ("benchmark", "verification"):
            del self.observations["before"][kind]["stages"][
                "retained_outer_exact_combination"
            ]
        with self.assertRaisesRegex(ValueError, "phase inventory"):
            adapter.engineering_comparison(self.index())

    def test_rejects_different_timer_scope(self):
        self.observations["after_off"]["benchmark"]["end_to_end_scope"] = (
            "DIFFERENT SYNTHETIC TIMER"
        )
        with self.assertRaisesRegex(ValueError, "timer scope"):
            adapter.engineering_comparison(self.index())

    def test_rejects_changed_receipt(self):
        index = self.index()
        (self.root / "after_off/benchmark.json").write_text("{}")
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            adapter.engineering_comparison(index)

    def test_cli_preserves_history_and_physical_selection(self):
        base = {
            "runs": 48,
            "raw_runs": [{"synthetic_sample": sample} for sample in range(48)],
            "selected_reduction_run": "historical-four-loop",
            "selected_five_loop_run": "historical-five-loop",
            "selected_physical_capture_run": "historical-physical",
            "physical_capture_runs": {
                "historical-physical": {"performance_eligible": False}
            },
        }
        self.write("base.json", base)
        output = self.root / "output.json"
        subprocess.run(
            [
                sys.executable,
                adapter.__file__,
                "--base",
                str(self.root / "base.json"),
                "--engineering",
                str(self.index()),
                "--output",
                str(output),
                "--select-current",
            ],
            check=True,
            capture_output=True,
        )
        merged = json.loads(output.read_text())
        for key, value in base.items():
            self.assertEqual(merged[key], value)
        selected = merged["selected_engineering_optimization_run"]
        self.assertIn(selected, merged["engineering_optimization_runs"])


if __name__ == "__main__":
    unittest.main()
