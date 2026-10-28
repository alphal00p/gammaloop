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


class DeliveryCheckpointContracts(unittest.TestCase):
    def setUp(self):
        self.fixture = CacheFileEngineeringContracts()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.base = copy.deepcopy(self.fixture.base)
        self.base["engineering_optimization_runs"]["synthetic-file"] = (
            adapter.cache_file_engineering(self.fixture.index(), self.base)
        )
        self.base["physical_capture_runs"]["historical-physical"][
            "observed_seconds"
        ] = {"producer_reduction": 100.0}
        self.documents = {
            "light_freeze": {
                "schema": "vakint.light-cache.user-directed-freeze.v1",
                "performance_claim": False,
                "future_tuning_allowed": False,
                "all_heldout_outcomes_used_for_selection": False,
                "compressed_bytes": 80,
                "decoded_bytes": 160,
                "sha256": "b" * 64,
                "selected_records": 56,
            },
            "authorization": {
                "schema": "vakint.q0-delivery-authorization.v1",
                "new_uncached_repetitions_waived": True,
                "required_physical_examples": ["zero-momentum"],
                "user_prompt_verbatim": "Synthetic user waiver fixture.",
            },
            "monitor": {
                "schema": "vakint.q0-monitor-responsiveness-actual-final-handoff.v1",
                "status": "PASS_ACTUAL_ONE_CACHED_Q0_MONITOR_AND_FULL55_58_146678_SCIENCE",
                "science": {
                    "all_original_mass_coupling_normalization_preserved": True,
                    "input_result_guard_and_family_identity_authenticated": True,
                },
                "monitor": {
                    "fault_count": 0,
                    "all_actual_worker_births_match": True,
                    "maximum_conservative_freshness_bound_seconds": 2.0,
                    "required_limit_seconds": 30,
                },
                "actual_OS_resources": {"signals_sent": []},
            },
            "heavy_sources": {
                "schema": "vakint.heavy-final-union-source-manifest.v1",
                "status": "complete",
                "pending_sources": [],
                "source_count": 2,
                "total_stored_record_references": 3,
                "sources": [
                    {"stored_record_references": 1},
                    {"stored_record_references": 2},
                ],
            },
        }

    def index(self):
        return self.fixture.write(
            "delivery-index.json",
            {
                "schema": "vakint.q0-delivery-evidence-index.v1",
                "cache_file_engineering_run": "synthetic-file",
                "historical_physical_run": "historical-physical",
                "records": {
                    name: self.fixture.write(name + ".json", value)
                    for name, value in self.documents.items()
                },
            },
        )["path"]

    def test_keeps_functional_evidence_separate_and_final_sizes_pending(self):
        run = adapter.delivery_checkpoint(self.index(), self.base)
        self.assertFalse(run["performance_eligible"])
        self.assertFalse(run["cache_performance_winner"])
        self.assertEqual(run["historical_uncached_seconds"], 100.0)
        self.assertTrue(all(value is None for value in run["pending"].values()))

    def test_rejects_cache_freeze_drift(self):
        self.documents["light_freeze"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "frozen LIGHT differs"):
            adapter.delivery_checkpoint(self.index(), self.base)

    def test_rejects_selection_using_heldout_or_claiming_performance(self):
        for key in ("all_heldout_outcomes_used_for_selection", "performance_claim"):
            with self.subTest(key=key):
                self.documents["light_freeze"][key] = True
                with self.assertRaisesRegex(ValueError, "non-performance freeze"):
                    adapter.delivery_checkpoint(self.index(), self.base)
                self.documents["light_freeze"][key] = False

    def test_rejects_absent_waiver(self):
        self.documents["authorization"]["new_uncached_repetitions_waived"] = False
        with self.assertRaisesRegex(ValueError, "waiver or physical scope"):
            adapter.delivery_checkpoint(self.index(), self.base)

    def test_rejects_unaccepted_monitor_and_stale_heartbeat(self):
        for key, value in (
            ("fault_count", 1),
            ("all_actual_worker_births_match", False),
            ("maximum_conservative_freshness_bound_seconds", 31.0),
        ):
            with self.subTest(key=key):
                original = self.documents["monitor"]["monitor"][key]
                self.documents["monitor"]["monitor"][key] = value
                with self.assertRaisesRegex(ValueError, "monitor"):
                    adapter.delivery_checkpoint(self.index(), self.base)
                self.documents["monitor"]["monitor"][key] = original

    def test_rejects_unclosed_heavy_sources_or_wrong_reference_total(self):
        for key, value in (("pending_sources", ["open"]), ("source_count", 3)):
            with self.subTest(key=key):
                original = self.documents["heavy_sources"][key]
                self.documents["heavy_sources"][key] = value
                with self.assertRaisesRegex(ValueError, "source inventory"):
                    adapter.delivery_checkpoint(self.index(), self.base)
                self.documents["heavy_sources"][key] = original

    def test_cli_preserves_every_historical_field_and_selection(self):
        self.fixture.write("delivery-base.json", self.base)
        output = self.root / "delivery-appended.json"
        subprocess.run(
            [
                sys.executable,
                adapter.__file__,
                "--base",
                str(self.root / "delivery-base.json"),
                "--delivery-checkpoint",
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
            self.assertEqual(merged[key], value)
        selected = merged["selected_delivery_checkpoint"]
        self.assertTrue(selected.startswith("q0-delivery-"))
        self.assertFalse(
            merged["delivery_checkpoints"][selected]["performance_eligible"]
        )

    def add_publication(self):
        self.documents.update(
            publication={
                "schema": "vakint.external-data.serial-publication.v1",
                "status": "PASS",
                "source": "1" * 40,
                "version": "v0.1.2",
                "immutable_release": "https://github.com/alphal00p/vakint-external-data/releases/tag/v0.1.2",
                "latest_alias": "https://github.com/alphal00p/vakint-external-data/releases/tag/data-latest",
                "stages": [
                    {"stage": "source-push", "exit_code": 0},
                    {"stage": "immutable-then-latest", "exit_code": 0},
                ],
            },
            release_assets={
                "schema": "vakint.frozen-light-full-release-pack-preflight.v1",
                "status": "PASS_TWO_INDEPENDENT_FULL_PACKS_AND_LOCAL_PUBLISHER_PREFLIGHT",
                "source_commit": "1" * 40,
                "checksums_byte_identical": True,
                "local_publisher_checked_assets_both_passed": True,
                "all_assets_below_decimal2GB": True,
                "assets": {
                    "manifest.json": {"bytes": 10, "sha256": "2" * 64},
                    "vakint-external-data-v0.1.2-5L.tar.gz": {
                        "bytes": 20,
                        "sha256": "3" * 64,
                    },
                },
            },
        )

    def test_publication_does_not_establish_download_or_offline_reuse(self):
        self.add_publication()
        run = adapter.delivery_checkpoint(self.index(), self.base)
        self.assertEqual(run["publication"]["version"], "v0.1.2")
        self.assertNotIn("published_release_url", run["pending"])
        self.assertIsNone(run["pending"]["actual_download_and_offline_reuse"])
        self.assertIsNone(run["pending"]["final_heavy_unique_records"])

    def test_rejects_reversed_publication_order_and_unbound_asset_source(self):
        self.add_publication()
        self.documents["publication"]["stages"].reverse()
        with self.assertRaisesRegex(ValueError, "ordered stages"):
            adapter.delivery_checkpoint(self.index(), self.base)
        self.documents["publication"]["stages"].reverse()
        self.documents["release_assets"]["source_commit"] = "0" * 40
        with self.assertRaisesRegex(ValueError, "source-bound release assets"):
            adapter.delivery_checkpoint(self.index(), self.base)

    def test_rejects_partial_final_heavy_or_q0_proof_groups(self):
        for key in ("heavy_export", "final_q0_benchmark"):
            with self.subTest(key=key):
                self.documents[key] = {"synthetic_incomplete_group": True}
                with self.assertRaisesRegex(
                    ValueError, "unknown delivery evidence index"
                ):
                    adapter.delivery_checkpoint(self.index(), self.base)
                del self.documents[key]


class StrictClosedCaptureContracts(unittest.TestCase):
    @staticmethod
    def capture(**native_fields):
        native = {
            "schema": "rustred.running-application-cache.v1",
            "closed": True,
            "writer_error": None,
            "last_error": None,
            "capture_metadata_exhausted": False,
            "generation": 1,
            "accepted_records": 2,
            "committed_records": 2,
            "pending_records": 0,
            "pending_bytes": 0,
            "capture_skipped_records": 0,
            **native_fields,
        }
        return {
            "terminal": True,
            "flush_complete": True,
            "flush_error": None,
            "snapshot_error": None,
            "reduction_error": None,
            "report": native,
            "flush_report": copy.deepcopy(native),
            "close": {"error": None, "report": copy.deepcopy(native)},
            "snapshot": {"schema": "rustred.running-application-cache-snapshot.v1"},
        }

    def test_clean_empty_replay_and_nonempty_capture_are_valid(self):
        for count in (0, 2):
            capture = self.capture(accepted_records=count, committed_records=count)
            self.assertEqual(
                adapter.strict_closed_capture(capture)["committed_records"], count
            )

    def test_rejects_uncommitted_skipped_exhausted_or_pending_records(self):
        for key, value in (
            ("accepted_records", 3),
            ("capture_skipped_records", 1),
            ("capture_metadata_exhausted", True),
            ("pending_records", 1),
            ("pending_bytes", 1),
            ("closed", False),
        ):
            # All report copies agree, so the strict invariant must reject this.
            with (
                self.subTest(key=key),
                self.assertRaisesRegex(ValueError, "clean flush/snapshot/close"),
            ):
                adapter.strict_closed_capture(self.capture(**{key: value}))

    def test_rejects_writer_and_last_errors(self):
        for key in ("writer_error", "last_error"):
            with (
                self.subTest(key=key),
                self.assertRaisesRegex(ValueError, "clean flush/snapshot/close"),
            ):
                adapter.strict_closed_capture(self.capture(**{key: "synthetic fault"}))

    def test_rejects_flush_snapshot_reduction_and_close_errors(self):
        for key in ("flush_error", "snapshot_error", "reduction_error", "close"):
            with self.subTest(key=key):
                capture = self.capture()
                if key == "close":
                    capture["close"]["error"] = "synthetic close fault"
                else:
                    capture[key] = "synthetic fault"
                with self.assertRaisesRegex(ValueError, "clean flush/snapshot/close"):
                    adapter.strict_closed_capture(capture)

    def test_missing_error_field_is_not_an_explicit_success(self):
        capture = self.capture()
        del capture["snapshot_error"]
        with self.assertRaisesRegex(ValueError, "clean flush/snapshot/close"):
            adapter.strict_closed_capture(capture)

    def test_missing_reduction_error_requires_successful_physical_execution(self):
        capture = self.capture()
        del capture["reduction_error"]
        with self.assertRaisesRegex(ValueError, "clean flush/snapshot/close"):
            adapter.strict_closed_capture(capture)
        execution = {
            "schema": "vakint.final-instrumented-q0-replay-execution.v1",
            "status": "PASS",
            "collector_exit_code": 0,
            "collector_wait_is_authoritative": True,
            "errors": [],
            "signals_sent": [],
        }
        self.assertEqual(
            adapter.strict_closed_capture(capture, physical_execution=execution),
            capture["report"],
        )
        with self.assertRaisesRegex(ValueError, "clean flush/snapshot/close"):
            adapter.strict_closed_capture(
                capture, physical_execution={**execution, "collector_exit_code": 1}
            )
        capture["reduction_error"] = "synthetic explicit reduction failure"
        with self.assertRaisesRegex(ValueError, "clean flush/snapshot/close"):
            adapter.strict_closed_capture(capture, physical_execution=execution)


class FinalCacheCompletionContracts(unittest.TestCase):
    """Synthetic bounded receipts; absent cache/physics paths must never be opened."""

    @staticmethod
    def record(value, sha="5" * 64):
        return {"evidence": {"sha256": sha}, "raw_receipt": value}

    def setUp(self):
        self.fixture = q0_fixtures.EngineeringQ0Contracts()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        rows = [
            {
                "id": name,
                "directory": f"/absent/synthetic-{name}",
                "generation": 1,
                "native_index_identity": {"sha256": name},
                "stored_record_references": count,
            }
            for name, count in (("one", 1), ("two", 2))
        ]
        native_origin = {"pid": 77, "start_ticks": 111, "uid": 1}
        controller_origin = {"pid": 88, "start_ticks": 222, "uid": 1}
        identity = {
            "compressed_bytes": 10,
            "decoded_bytes": 20,
            "compressed_blake3": [1] * 32,
            "decoded_blake3": [2] * 32,
        }
        artifact = {
            "path": "/absent/synthetic-heavy.rrcache.gz",
            "sha256": "f" * 64,
            "compressed_bytes": 10,
            "decoded_bytes": 20,
            "typed_Native_identity": identity,
        }
        census = {
            "schema": "vakint.final-thirteen-source-census-postflight.v1",
            "status": "PASS_ACTUAL_FULL13_CENSUS_ORIGINAL_PIDFDS_TERMINAL",
            "source_count": 2,
            "input_all_version_references": 3,
            "unique_records": 2,
            "identical_duplicates": 1,
            "full_row_conflicts": 0,
            "original12_source_rows_and_seals_unchanged_superset": True,
            "Native_seconds": 1.0,
            "controller_elapsed_seconds": 1.5,
        }
        exported = {
            "schema": "vakint.final-thirteen-source-export-postflight.v1",
            "status": "PASS_ACTUAL_ONE_FINAL13_EXPORT_ORIGINAL_PIDFDS_TERMINAL",
            "source_manifest": {"sha256": "8" * 64},
            "actual_accepted_census": {"sha256": "7" * 64},
            "request": {"sha256": "6" * 64},
            "source_count": 2,
            "all13_input_references": 3,
            "actual_unique_records": 2,
            "identical_duplicates": 1,
            "actual_writer_statistics": {"records": 2, "unique_coefficients": 4},
            "all13_source_seals_unchanged_and_equal_actual_census": True,
            "original12_source_rows_and_index_seals_unchanged_superset": True,
            "all_census_record_duplicate_coefficient_payload_and_row_spool_metrics_match": True,
            "artifact": artifact,
            "Native_seconds": 2.0,
            "controller_elapsed_seconds": 2.5,
        }
        imported = {
            "schema": "vakint.heavy-native-phase.v1",
            "phase": "reimport",
            "status": "PASS",
            "performance_claim": False,
            "request_sha256": "6" * 64,
            "helper_executable": {"sha256": "4" * 64},
            "pid": native_origin["pid"],
            "source_seals": [
                {"source_id": row["id"], **{k: v for k, v in row.items() if k != "id"}}
                for row in rows
            ],
            "native_revision": "c" * 40,
            "symbolica_revision": "d" * 40,
            "main_revision": "0" * 40,
            "seconds": 3.0,
            "details": {
                "input_records": 3,
                "all_exact_record_values_and_ordered_guards_equal": True,
                "state_imports_shared": True,
                "cache_path": artifact["path"],
                "cache_sha256": artifact["sha256"],
                "cache_identity": identity,
                "source_visits": [
                    {
                        "source_id": row["id"],
                        "emitted_references": row["stored_record_references"],
                    }
                    for row in rows
                ],
            },
        }
        controller = {
            "schema": "vakint.heavy-finalization-controller-phase.v1",
            "phase": "reimport",
            "status": "PASS",
            "request_sha256": "6" * 64,
            "child_returncode": 0,
            "signals_sent": [],
            "driver_receipt_sha256": "5" * 64,
            "helper_sha256": "4" * 64,
            "child": native_origin,
            "controller": controller_origin,
            "elapsed_seconds": 3.5,
        }
        terminal = {
            "schema": "vakint.final-heavy-phase-original-pidfds-terminal.v1",
            "phase": "reimport",
            "all_original_PIDFDs_terminal": True,
            "signals_sent": [],
            "terminal": {
                name: {"original_PIDFD_terminal": True, "original": origin}
                for name, origin in (
                    ("helper", native_origin),
                    ("controller", controller_origin),
                )
            },
        }
        self.heavy_records = {
            "heavy_sources": self.record(
                {
                    "source_count": 2,
                    "total_stored_record_references": 3,
                    "sources": rows,
                },
                "8" * 64,
            ),
            "heavy_census": self.record(census, "7" * 64),
            "heavy_export": self.record(exported),
            "heavy_reimport": self.record(imported),
            "heavy_reimport_controller": self.record(controller),
            "heavy_reimport_terminals": self.record(terminal),
        }
        self.heavy = adapter.final_heavy_cache(self.heavy_records)
        observation = copy.deepcopy(self.fixture.observations["after_off"])
        report, execution, verified = (
            observation[key] for key in ("benchmark", "execution", "verification")
        )
        report.update(
            warm_requested=False,
            executable_sha256="9" * 64,
            input_loading_seconds=0.1,
            data_loading_seconds=0.2,
            cache_loading_seconds=0.3,
            cache_identity={"cache_sha256": artifact["sha256"]},
            cache_identity_after={"cache_sha256": artifact["sha256"]},
        )
        report["compiled_cohort"].update(
            native_revision="c" * 40, kernel_revision="d" * 40
        )
        worker, parent, guard = (
            {"pid": pid, "start_ticks": pid * 10, "uid": 1} for pid in (101, 102, 103)
        )
        execution.update(
            worker=worker,
            controller=parent,
            guard_registration={
                "guard": {key: guard[key] for key in ("pid", "start_ticks")}
            },
            controller_inclusive_elapsed_seconds=4.0,
        )
        capture = StrictClosedCaptureContracts.capture(
            accepted_records=0, committed_records=0
        )
        del capture[
            "reduction_error"
        ]  # Physical collector schema; execution proves success.
        capture["path"] = "/absent/synthetic-empty-capture"
        capture["read_delta"] = False
        verified.update(
            schema="vakint.shared-file-physical-q0-engineering-verification.v1",
            cache_identity=report["cache_identity"],
            file_attachment={
                "path": artifact["path"],
                "known_sha256": artifact["sha256"],
            },
            running_cache_path=capture["path"],
            running_cache_read_delta=False,
            running_capture_durably_closed=True,
            actual_native_work={"rule_applications": 0},
            actual_persistent_work_and_captures={"route_hits": 1, "misses": 0},
        )
        report["running_cache_read_delta"] = False
        report["saved_application_cache"] = {
            "native": verified["actual_native_work"],
            "persistent": verified["actual_persistent_work_and_captures"],
        }
        monitor = {
            "schema": "vakint.q0-monitor-responsiveness.v1",
            "status": "PASS",
            "fault_count": 0,
            "performance_claim": False,
            "actual_worker_terminal": True,
            "worker": worker,
            "max_conservative_freshness_bound_seconds": 1.0,
        }
        q0_terminals = {
            "schema": "vakint.q0-monitor-original-three-pidfds-terminal.v1",
            "all_three_original_PIDFDs_terminal": True,
            "signals_sent": [],
            "registration": {
                "originals": copy.deepcopy(
                    {"worker": worker, "controller": parent, "guard": guard}
                )
            },
            "terminal": {
                name: {"original_PIDFD_terminal": True, "original": origin}
                for name, origin in (
                    ("worker", worker),
                    ("controller", parent),
                    ("guard", guard),
                )
            },
        }
        self.q0_records = {
            "final_q0_benchmark": self.record(report),
            "final_q0_execution": self.record(execution),
            "final_q0_verification": self.record(verified),
            "final_q0_monitor": self.record(monitor),
            "final_q0_terminals": self.record(q0_terminals),
            "final_q0_capture": self.record(capture),
            "final_q0_coefficient_proof": {
                "evidence": verified["full58_proof"],
                "raw_receipt": observation["coefficient_proof"],
            },
        }
        self.reference = copy.deepcopy(report)

    def test_final_file_observation_is_single_and_never_a_matched_median(self):
        result = adapter.final_cached_q0(self.q0_records, self.heavy, self.reference)
        self.assertEqual(result["sample_count"], 1)
        self.assertFalse(result["performance_eligible"])
        self.assertIsNone(result["pure_ibp_seconds"])

    def test_rejects_false_original_native_guard_equality(self):
        self.heavy_records["heavy_reimport"]["raw_receipt"]["details"][
            "all_exact_record_values_and_ordered_guards_equal"
        ] = False
        with self.assertRaisesRegex(ValueError, "ordered-guard round-trip"):
            adapter.final_heavy_cache(self.heavy_records)

    def test_rejects_final_file_identity_or_collector_cohort_mismatch(self):
        self.heavy["artifact"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "science/file identity"):
            adapter.final_cached_q0(self.q0_records, self.heavy, self.reference)
        self.heavy["artifact"]["sha256"] = "f" * 64
        self.q0_records["final_q0_benchmark"]["raw_receipt"]["compiled_cohort"][
            "native_revision"
        ] = "0" * 40
        with self.assertRaisesRegex(ValueError, "saved-input/source observation"):
            adapter.final_cached_q0(self.q0_records, self.heavy, self.reference)

    def test_rejects_missing_full_scalar_parity(self):
        self.q0_records["final_q0_verification"]["raw_receipt"][
            "full58_scalar_coefficient_parity"
        ] = False
        with self.assertRaisesRegex(ValueError, "coefficient proof"):
            adapter.final_cached_q0(self.q0_records, self.heavy, self.reference)

    def test_guard_uid_must_match_authenticated_original_registration(self):
        self.q0_records["final_q0_terminals"]["raw_receipt"]["terminal"]["guard"][
            "original"
        ]["uid"] = 99
        with self.assertRaisesRegex(ValueError, "actual original terminals"):
            adapter.final_cached_q0(self.q0_records, self.heavy, self.reference)


if __name__ == "__main__":
    unittest.main()
