"""Synthetic JSON contracts; referenced cache files deliberately do not exist."""

import copy
import json
import subprocess
import sys
import unittest
from pathlib import Path

import merge_physical_capture_results as adapter
import test_engineering_q0_results as q0_fixtures


class CacheFileEngineeringContracts(unittest.TestCase):
    def setUp(self):
        fixture = q0_fixtures.EngineeringQ0Contracts()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.fixture = fixture
        self.root = fixture.root
        for observation in fixture.observations.values():
            observation["benchmark"]["executable_sha256"] = "f" * 64
        prior = adapter.engineering_comparison(fixture.index())
        self.baseline_id = "synthetic-coefficient-engineering"
        self.base = {
            "runs": 48,
            "raw_runs": [{"synthetic_sample": sample} for sample in range(48)],
            "selected_reduction_run": "historical-four-loop",
            "selected_five_loop_run": "historical-five-loop",
            "selected_physical_capture_run": "historical-physical",
            "selected_engineering_optimization_run": self.baseline_id,
            "physical_capture_runs": {
                "historical-physical": {"performance_eligible": False}
            },
            "engineering_optimization_runs": {self.baseline_id: prior},
        }
        report = copy.deepcopy(
            prior["observations"]["after_off"]["records"]["benchmark"]["raw_receipt"]
        )
        report.update(
            warm_requested=False,
            reduction_seconds=1.5,
            end_to_end_seconds=4.5,
            data_loading_seconds=0.5,
            cache_loading_seconds=2.0,
            cache_identity={"cache_sha256": "6" * 64},
            cache_identity_after={"cache_sha256": "6" * 64},
            running_cache_read_delta=False,
            saved_application_cache={
                "native": {"cache_hits": 2, "rule_applications": 0},
                "persistent": {"captured_records": 0, "route_hits": 2, "misses": 0},
            },
        )
        verified = copy.deepcopy(
            prior["observations"]["after_off"]["records"]["verification"]["raw_receipt"]
        )
        verified["running_capture_durably_closed"] = True
        execution = copy.deepcopy(
            prior["observations"]["after_off"]["records"]["execution"]["raw_receipt"]
        )
        execution["elapsed_seconds"] = 5.5
        heavy = {
            "schema": "vakint.full-q0-shared-file-engineering-final-handoff.v1",
            "status": "PASS_SINGLE_COLD_FULL_Q0_ENGINEERING_AND_FULL_SCIENTIFIC_VERIFICATION",
            "performance_claim": False,
            "full_actual_benchmark_metadata": report,
            "full_actual_scientific_verification_metadata": verified,
            "full_actual_physical_completion_metadata": execution,
            "Native_work": copy.deepcopy(report["saved_application_cache"]),
            "clock_summary_seconds": {
                "public_call": 1.5,
                "collector_E2E": 4.5,
                "data_loading": 0.5,
                "Native_FILE_loading_auth_decode_index": 2.0,
                "inclusive_controller_elapsed_after_preflight": 5.5,
                "decoder_only": None,
                "pure_IBP": None,
            },
            "legacy_FILE_contract": {
                "externally_supplied_typed_FileIdentity": False,
                "Native_framing_and_semantic_binding_checked_on_load": True,
                "actual_before": report["cache_identity"],
                "actual_after": report["cache_identity_after"],
            },
            "clock_boundaries": {
                "wholeFILE_SHA_scan_count": 2,
                "SHA_durations_independently_measured": False,
                "controller_minus_collector_E2E_is_not_pure_SHA": True,
            },
        }
        candidates = []
        pilots = []
        for label, additional, compressed_bytes in (("A", 0, 64), ("B", 2, 80)):
            cache = {
                "path": str(self.root / f"absent-cache-{label}.rrcache.gz"),
                "bytes": compressed_bytes,
                "sha256": label.lower() * 64,
            }
            count = 54 + additional
            candidate = {
                "schema": "vakint.physical-route-budget-candidate-result.v1"
                if label == "A"
                else "vakint.physical-route-budget-additions-candidate-result.v1",
                "physical_route_limit": 54,
                "additional_keys"
                if label == "A"
                else "additional_apply_keys": additional,
                "native_export_accepted": True,
                "exact_record_and_ordered_guard_equality": True,
                "final_held_out_used": False,
                "cache": cache,
                "actual_compressed_bytes": compressed_bytes,
                "actual_decoded_bytes": compressed_bytes * 2,
                "native_counts": {
                    "selected_records": count,
                    "reimported_records": count,
                    "complete_records": count,
                    "missing_keys": 0,
                    "partial_records": 0,
                },
                "phase_seconds": {
                    phase: 0.25
                    for phase in (
                        "load_saved_application",
                        "encode_selected_completed_cache",
                        "write_selected_cache",
                        "reimport_selected_cache",
                        "compare_selected_native_visitor",
                    )
                },
            }
            capture_count = 4 if label == "A" else 3
            pilots.append(
                {
                    "label": label,
                    "cache_file": cache,
                    "parity": {
                        "full_canonical_value_bytes_equal": True,
                        "all_guards_bytes_equal": True,
                        "canonical_master_families_equal": True,
                        "catalogue_closed": True,
                        "frozen_input_owner_roots_equal": True,
                        "fixed_reduction_settings_equal": True,
                        "unknown_count": 0,
                    },
                    "closed_running_capture": {
                        "report": {
                            "closed": True,
                            "pending_records": 0,
                            "pending_bytes": 0,
                            "writer_error": None,
                            "committed_records": capture_count,
                            "accepted_records": capture_count,
                        },
                        "all_new_completed_nodes_preserved": True,
                        "read_delta": False,
                    },
                    "actual_cache_statistics": {
                        "persistent": {"captured_records": capture_count}
                    },
                    "timing": {
                        "evaluation_seconds": 8.0,
                        "data_loading_seconds": 2.0,
                        "cache_loading_seconds": 1.0,
                        "saved_application_including_owner_seconds": 6.0,
                        "stages": {"saved_rule_application": [2.0, 4.0]},
                        "native_timing_scope": {"pure_ibp_seconds": None},
                    },
                }
            )
            candidates.append(candidate)
        light_baseline = {
            "schema": "vakint.five-loop-cache-measurement.v1",
            "case_id": "synthetic-tuning-case",
            "identity": {
                "executable_sha256": "1" * 64,
                "workload_sha256": "2" * 64,
            },
            "success": True,
            "worker_success": True,
            "evaluation_seconds": 6.0,
            "saved_application_including_owner_seconds": 5.0,
            "ibp_excluding_owner_seconds": 4.5,
            "cache_statistics": {"native": {"rule_applications": 7}},
            "native_timing_scope": {"pure_ibp_seconds": None},
            "stages": {"saved_rule_application": [5.0]},
        }
        self.documents = {
            "heavy": heavy,
            "light": {
                "schema": "vakint.light-utility-engineering-pilots.v1",
                "status": "PASS",
                "single_math_worker": True,
                "A_B_full_canonical_result_and_all_guards_equal": True,
                "case_id": light_baseline["case_id"],
                "full_workload_sha256": "2" * 64,
                "compiled_executable_sha256": "1" * 64,
                "baseline": {
                    "evaluation_seconds": 6.0,
                    "saved_application_including_owner_seconds": 5.0,
                    "ibp_excluding_owner_seconds": 4.5,
                    "pure_ibp_seconds": None,
                    "new_native_rules": 7,
                    "stages": light_baseline["stages"],
                },
                "pilots": pilots,
            },
            "route_budget": {
                "schema": "vakint.physical-route-budget-pilot-completion.v1",
                "final_held_out_used": False,
            },
            "additions_budget": {
                "schema": "vakint.physical-route-budget-additions-pilot-completion.v1",
                "final_performance_selection": False,
                "final_held_out_used": False,
            },
            "candidate_a": candidates[0],
            "candidate_b": candidates[1],
            "light_baseline": light_baseline,
        }

    def write(self, name, value):
        reference = self.fixture.write(name, value)
        reference["bytes"] = Path(reference["path"]).stat().st_size
        return reference

    def index(self):
        documents = copy.deepcopy(self.documents)
        candidate_a = self.write("candidate-a.json", documents["candidate_a"])
        candidate_b = self.write("candidate-b.json", documents["candidate_b"])
        documents["route_budget"]["accepted_result"] = candidate_a
        documents["additions_budget"]["accepted_candidates"] = [
            {"cache": documents["candidate_b"]["cache"], "result": candidate_b}
        ]
        documents["light"]["prior_measurement"] = self.write(
            "light-baseline.json", documents["light_baseline"]
        )
        self.write(
            "file-index.json",
            {
                "schema": "vakint.cache-file-engineering-evidence-index.v1",
                "baseline_engineering_run": self.baseline_id,
                "records": {
                    name: self.write(f"{name}.json", documents[name])
                    for name in ("heavy", "light", "route_budget", "additions_budget")
                },
            },
        )
        return self.root / "file-index.json"

    def test_accepts_metadata_without_opening_cache_and_preserves_boundaries(self):
        original = copy.deepcopy(self.base)
        run = adapter.cache_file_engineering(self.index(), self.base)
        self.assertEqual(self.base, original)
        self.assertEqual(run["status"], "verified-engineering-observations")
        self.assertFalse(run["performance_eligible"])
        self.assertFalse(run["matched_serialized_repetitions"])
        self.assertFalse(run["final_light_selection"])
        self.assertEqual(
            run["heavy_observed_change_percent"],
            {"public_call": -25.0, "collector_E2E": 50.0},
        )
        heavy = run["records"]["heavy"]["raw_receipt"]
        self.assertEqual(heavy["clock_summary_seconds"]["collector_E2E"], 4.5)
        self.assertEqual(
            heavy["clock_summary_seconds"][
                "inclusive_controller_elapsed_after_preflight"
            ],
            5.5,
        )
        self.assertIsNone(heavy["clock_summary_seconds"]["decoder_only"])
        self.assertIsNone(heavy["clock_summary_seconds"]["pure_IBP"])
        for pilot in run["records"]["light"]["raw_receipt"]["pilots"]:
            self.assertFalse(Path(pilot["cache_file"]["path"]).exists())

    def test_rejects_changed_source_metadata(self):
        index = self.index()
        path = self.root / "heavy.json"
        path.write_text(path.read_text().replace("vakint.", "xakint.", 1))
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            adapter.cache_file_engineering(index, self.base)

    def test_rejects_wrong_baseline_or_unaccepted_completion(self):
        heavy = self.documents["heavy"]
        mutations = (
            (heavy["full_actual_benchmark_metadata"], "result_sha256", "0" * 64),
            (
                heavy["full_actual_physical_completion_metadata"],
                "collector_wait_is_authoritative",
                False,
            ),
        )
        for target, key, value in mutations:
            with self.subTest(field=key):
                original = target[key]
                target[key] = value
                with self.assertRaisesRegex(ValueError, "full-parity baseline"):
                    adapter.cache_file_engineering(self.index(), self.base)
                target[key] = original

    def test_rejects_mismatched_measured_phase(self):
        self.documents["heavy"]["full_actual_scientific_verification_metadata"][
            "stages"
        ]["saved_rule_application"]["total_seconds"] += 1
        with self.assertRaisesRegex(ValueError, "phase differs"):
            adapter.cache_file_engineering(self.index(), self.base)

    def test_rejects_file_capture_or_counter_summary_drift(self):
        heavy = self.documents["heavy"]
        verified = heavy["full_actual_scientific_verification_metadata"]
        for target, field, value in (
            (verified, "master_identity_union_count", 2),
            (verified, "full_scalar_coefficient_count", 3),
            (verified, "running_capture_durably_closed", False),
            (heavy["Native_work"]["native"], "cache_hits", 3),
            (heavy["full_actual_benchmark_metadata"], "running_cache_read_delta", True),
        ):
            with self.subTest(field=field):
                original = target[field]
                target[field] = value
                with self.assertRaisesRegex(ValueError, "full-parity baseline"):
                    adapter.cache_file_engineering(self.index(), self.base)
                target[field] = original

    def test_rejects_relabelled_clocks_and_measured_zero_for_unknown_phases(self):
        clocks = self.documents["heavy"]["clock_summary_seconds"]
        for field, value in (
            ("collector_E2E", 5.5),
            ("decoder_only", 0),
            ("pure_IBP", 0),
        ):
            with self.subTest(field=field):
                original = clocks[field]
                clocks[field] = value
                with self.assertRaisesRegex(ValueError, "clock summaries"):
                    adapter.cache_file_engineering(self.index(), self.base)
                clocks[field] = original

    def test_rejects_cache_above_95mb_even_when_metadata_agrees(self):
        candidate = self.documents["candidate_a"]
        candidate["actual_compressed_bytes"] = 95_000_001
        candidate["cache"]["bytes"] = 95_000_001
        with self.assertRaisesRegex(ValueError, "count/byte/reimport"):
            adapter.cache_file_engineering(self.index(), self.base)

    def test_rejects_incomplete_native_reimport(self):
        for field in ("complete_records", "reimported_records", "missing_keys"):
            with self.subTest(field=field):
                counts = self.documents["candidate_b"]["native_counts"]
                original = counts[field]
                counts[field] += 1
                with self.assertRaisesRegex(ValueError, "count/byte/reimport"):
                    adapter.cache_file_engineering(self.index(), self.base)
                counts[field] = original

    def test_rejects_missing_displayed_export_phase(self):
        del self.documents["candidate_b"]["phase_seconds"]["write_selected_cache"]
        with self.assertRaisesRegex(ValueError, "phase"):
            adapter.cache_file_engineering(self.index(), self.base)

    def test_rejects_drift_from_original_light_measurement(self):
        baseline = self.documents["light_baseline"]
        for target, field, value in (
            (baseline, "evaluation_seconds", 7.0),
            (baseline, "case_id", "different-synthetic-case"),
            (baseline["identity"], "executable_sha256", "0" * 64),
            (baseline["cache_statistics"]["native"], "rule_applications", 8),
        ):
            with self.subTest(field=field):
                original = target[field]
                target[field] = value
                with self.assertRaisesRegex(ValueError, "historical baseline"):
                    adapter.cache_file_engineering(self.index(), self.base)
                target[field] = original

    def test_rejects_incomplete_exact_parity(self):
        parity = self.documents["light"]["pilots"][0]["parity"]
        parity["all_guards_bytes_equal"] = False
        with self.assertRaisesRegex(ValueError, "parity or durable capture"):
            adapter.cache_file_engineering(self.index(), self.base)

    def test_rejects_unclosed_pending_or_inconsistent_capture(self):
        capture = self.documents["light"]["pilots"][1]["closed_running_capture"]
        for field, value in (
            ("closed", False),
            ("pending_records", 1),
            ("writer_error", "synthetic write failure"),
            ("committed_records", 2),
        ):
            with self.subTest(field=field):
                original = capture["report"][field]
                capture["report"][field] = value
                with self.assertRaisesRegex(ValueError, "parity or durable capture"):
                    adapter.cache_file_engineering(self.index(), self.base)
                capture["report"][field] = original

    def test_rejects_light_saved_call_drift(self):
        self.documents["light"]["pilots"][0]["timing"]["stages"][
            "saved_rule_application"
        ] = [1.0, 4.0]
        with self.assertRaisesRegex(ValueError, "saved-call samples"):
            adapter.cache_file_engineering(self.index(), self.base)

    def test_cli_appends_without_replacing_history_or_prior_selections(self):
        self.write("base.json", self.base)
        output = self.root / "appended.json"
        subprocess.run(
            [
                sys.executable,
                adapter.__file__,
                "--base",
                str(self.root / "base.json"),
                "--file-engineering",
                str(self.index()),
                "--output",
                str(output),
                "--select-current",
            ],
            check=True,
            capture_output=True,
        )
        merged = json.loads(output.read_text())
        for key, value in self.base.items():
            if key != "engineering_optimization_runs":
                self.assertEqual(merged[key], value)
        runs = merged["engineering_optimization_runs"]
        self.assertEqual(
            runs[self.baseline_id],
            self.base["engineering_optimization_runs"][self.baseline_id],
        )
        selected = merged["selected_cache_file_engineering_run"]
        self.assertEqual(set(runs), {self.baseline_id, selected})
        self.assertTrue(selected.startswith("cache-file-engineering-"))
        self.assertFalse(runs[selected]["performance_eligible"])


if __name__ == "__main__":
    unittest.main()
