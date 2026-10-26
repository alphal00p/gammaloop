"""Synthetic file-contract tests only; no reduction or performance measurement."""

import copy
import hashlib
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import merge_five_loop_results as merge


class ReportContractTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(
            prefix="synthetic-benchmark-report-"
        )
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.binary = self.root / "TEST-NOT-AN-EXECUTABLE"
        self.binary.write_text("SYNTHETIC TEST FIXTURE ONLY")
        self.cache = self.root / "TEST-NOT-A-CACHE"
        self.cache.write_text("SYNTHETIC TEST FIXTURE ONLY")
        self.stamps = {
            "VAKINT_CACHE_BUILD_GIT_REVISION": "1" * 40,
            "VAKINT_CACHE_BUILD_RUSTRED_REVISION": "2" * 40,
            "VAKINT_CACHE_BUILD_SYMBOLICA_REVISION": "3" * 40,
            "VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256": "4" * 64,
            "VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256": "5" * 64,
            "VAKINT_CACHE_BUILD_PROFILE": "release",
        }
        self.artifact = self.write(
            "artifact.json",
            {
                "test_fixture": True,
                "source_drift": [],
                "source_unchanged": True,
                "stamps": self.stamps,
                "actual_codegen": {"optimization": "3", "fat_lto": True, "test": False},
                "binary": {
                    "sha256": merge.digest(self.binary),
                    "immutable_path": str(self.binary),
                },
            },
        )
        self.threads = {
            key: "1"
            for key in ("RAYON_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
        }
        self.records = []
        self.window_rows = []

    def write(self, name, value):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, indent=2))
        return path

    def sample(self, variant, index, case="quadratic-four-rung-uv"):
        directory = self.root / f"{case}/{variant}-{index}"
        directory.mkdir(parents=True)
        for name, content in (
            ("result.symbolica", "TEST_EXACT_RESULT"),
            ("vakint-input.symbolica", "TEST_INPUT"),
            ("nonzero-conditions.symbolica", "TEST_GUARD\n"),
            ("unknown-masters.json", "[]"),
        ):
            (directory / name).write_text(content)
        scientific = {
            name: "e" * 64
            for name in (
                "physical.rs",
                "projector.rs",
                "vacuum.rs",
                "ladder.dot",
                "qcd.json",
                "phi3.json",
            )
        }
        sources = {**scientific, "mod.rs": "f" * 64}
        conditions = ["TEST_GUARD"] if case == "quadratic-four-rung-uv" else []
        contributions = 21 if case == "quadratic-four-rung-uv" else 1
        report = {
            "test_fixture": True,
            "schema": "vakint.qcd-ladder.benchmark.v1",
            "projection": merge.PRIMARY[case],
            "coefficient_representation": "sparse",
            "executable_sha256": merge.digest(self.binary),
            "source_dirty": False,
            "source_commit": self.stamps["VAKINT_CACHE_BUILD_GIT_REVISION"],
            "fresh_native_preparation": True,
            "retained_result_available": True,
            "input_sha256": merge.digest(directory / "vakint-input.symbolica"),
            "fixture_sources_sha256": sources,
            "consumer_fixture_sources_sha256": sources,
            "input_bytes": len("TEST_INPUT"),
            "contributions": contributions,
            "scalar_terms": 30,
            "prepared_nonzero_conditions": conditions,
            "data_manifest_sha256": "7" * 64,
            "result_sha256": merge.digest(directory / "result.symbolica"),
            "guards": {
                "count": 1,
                "file": "nonzero-conditions.symbolica",
                "sha256": merge.digest(directory / "nonzero-conditions.symbolica"),
                "includes_preparation_conditions": True,
            },
            "catalogue_closure": {
                "closed": True,
                "unknown_count": 0,
                "master_count": 1,
                "unknown_file": "unknown-masters.json",
                "unknown_sha256": merge.digest(directory / "unknown-masters.json"),
            },
            "reduction_completed": True,
            "persistent_cache": str(self.cache) if variant == "cached" else None,
            "persistent_cache_attachments_before": {
                str(loop): False for loop in range(1, 6)
            },
            "persistent_cache_attachments_during_evaluation": {
                str(loop): loop == 5 and variant == "cached" for loop in range(1, 6)
            },
            "preparation_seconds": 2 + index,
            "data_loading_seconds": 3 + index,
            "cache_loading_seconds": 0 if variant == "uncached" else 1,
            "reduction_seconds": 10 + index if variant == "uncached" else 5 + index,
            "end_to_end_seconds": 25 + index,
            "stages": {"saved_rule_application": [2, 3]},
            "preparation": {"graph_generation_seconds": 0.5},
        }
        path = self.write(
            str(directory.relative_to(self.root) / "benchmark.json"), report
        )
        launch = self.write(
            str(directory.relative_to(self.root) / "launch.json"),
            {
                "test_fixture": True,
                "final_performance_eligible": True,
                "fresh_process": True,
                "binary_sha256": merge.digest(self.binary),
                "case_id": case,
                "variant": variant,
                "sample": index,
                "cpu_affinity": "TEST_CPU",
                "thread_environment": self.threads,
                "cache_sha256": merge.digest(self.cache)
                if variant == "cached"
                else None,
            },
        )
        generator = self.write(
            f"{case}-generator.json",
            {
                "test_fixture": True,
                "schema": "vakint.current-quadratic-preparation.scientific-parity.v1"
                if contributions == 21
                else "vakint.current-q0-preparation.scientific-parity.v1",
                "status": "PASS_ALL21_FRESH_CURRENT_DEPENDENCY_RAW_INPUT_PARITY"
                if contributions == 21
                else "PASS_Q0_FRESH_CURRENT_DEPENDENCY_RAW_INPUT_PARITY",
                "collector_exit_code": 0,
                "source_drift": [],
                "projection": merge.PRIMARY[case],
                "contributions": contributions,
                "input": {
                    "sha256": report["input_sha256"],
                    "bytes": report["input_bytes"],
                    "exact_canonical_bytes_equal_accepted_input": True,
                },
                "projector_guards": {
                    "sha256": merge.hashlib.sha256(
                        "\n".join(conditions).encode()
                    ).hexdigest(),
                    "values": conditions,
                    "bytes_equal_accepted": True,
                },
                "compiled_fixture_source_hashes": sources,
                "source_commit": "9" * 40,
                "immutable_elf_sha256": "8" * 64,
            },
        )
        launch_doc = merge.read(launch)
        launch_doc["generator_acceptance"] = merge.evidence(generator)
        launch.write_text(json.dumps(launch_doc))
        refs = {
            "record": merge.evidence(path),
            "launcher": merge.evidence(launch),
            "artifact": merge.evidence(self.artifact),
            "case_id": case,
        }
        self.records.append(refs)
        start = 100 + len(self.window_rows) * 100
        self.window_rows.append(
            {
                "report_sha256": merge.digest(path),
                "launcher_sha256": merge.digest(launch),
                "artifact_sha256": merge.digest(self.artifact),
                "case_id": case,
                "variant": variant,
                "sample": index,
                "started_unix_seconds": start,
                "finished_unix_seconds": start + 50,
                "exit_code": 0,
                "fresh_process": True,
            }
        )
        return path

    def paired(self, count=3, case="quadratic-four-rung-uv"):
        for index in range(count):
            for variant in ("uncached", "cached"):
                self.sample(variant, index, case)

    def index(self, qualified=True, window_edit=None):
        window = {
            "schema": "vakint.serialized-performance-window.v1",
            "test_fixture": True,
            "final_performance_eligible": True,
            "accepted": True,
            "serialized": True,
            "all_processes_drained": True,
            "source_unchanged": True,
            "started_unix_seconds": 0,
            "finished_unix_seconds": 10000,
            "cpu_affinity": "TEST_CPU",
            "thread_environment": self.threads,
            "samples": self.window_rows,
        }
        if window_edit:
            window_edit(window)
        window_path = self.write("window.json", window)
        rows = copy.deepcopy(self.records)
        if qualified:
            for row in rows:
                row["window"] = merge.evidence(window_path)
        return self.write(
            "index.json",
            {"schema": "vakint.five-loop-report-evidence.v1", "records": rows},
        )

    def evaluate(self, **kwargs):
        return merge.primary_records(merge.ReportEvidence(self.index(**kwargs)))

    def test_deferred_primary_preserves_functional_records_without_measurements(self):
        self.paired()
        self.paired(case="zero-momentum-diagnostic")
        provenance = merge.ReportEvidence(self.index())
        result = merge.primary_records(
            provenance,
            {
                "required": ["zero-momentum-diagnostic"],
                "deferred": ["quadratic-four-rung-uv"],
            },
        )
        deferred = result["quadratic-four-rung-uv"]
        self.assertEqual(deferred["status"], "deferred")
        self.assertEqual(len(deferred["functional_records"]), 6)
        self.assertTrue(
            all(
                not record["performance_eligible"]
                for record in deferred["functional_records"]
            )
        )
        for field in (
            "uncached",
            "cached",
            "catalogue_closure",
            "exact_output_and_guard_parity",
        ):
            self.assertIsNone(deferred[field])
        self.assertEqual(result["zero-momentum-diagnostic"]["status"], "measured")
        self.assertEqual(provenance.used, set(provenance.records))
        # The same historical evidence retains its original meaning without a scope.
        self.assertEqual(
            self.evaluate()["quadratic-four-rung-uv"]["status"], "measured"
        )

    def test_deferral_does_not_bypass_exact_retained_payload_checks(self):
        report = self.sample("uncached", 0)
        (report.parent / "result.symbolica").write_text("changed retained result")
        with self.assertRaisesRegex(ValueError, "exact result changed"):
            merge.primary_records(
                merge.ReportEvidence(self.index(False)),
                {
                    "required": ["zero-momentum-diagnostic"],
                    "deferred": ["quadratic-four-rung-uv"],
                },
            )

    def test_explicit_scope_does_not_invent_primary_completion(self):
        result = merge.primary_records(
            merge.ReportEvidence(self.index(False)),
            {
                "required": ["zero-momentum-diagnostic"],
                "deferred": ["quadratic-four-rung-uv"],
            },
        )
        self.assertEqual(result["quadratic-four-rung-uv"]["status"], "deferred")
        self.assertEqual(result["zero-momentum-diagnostic"]["status"], "pending")

    def change_report(self, index, edit):
        path = Path(self.records[index]["record"]["path"])
        report = merge.read(path)
        edit(report)
        path.write_text(json.dumps(report))
        self.records[index]["record"] = merge.evidence(path)
        self.window_rows[index]["report_sha256"] = merge.digest(path)

    def test_three_matched_pairs_and_unmeasured_metrics(self):
        for case in merge.PRIMARY:
            self.paired(case=case)
        for item in self.evaluate().values():
            self.assertEqual(item["status"], "measured")
            self.assertEqual(
                item["uncached"]["metrics"]["reduction_seconds"]["median_seconds"], 11
            )
            self.assertEqual(
                item["cached"]["metrics"]["reduction_seconds"]["min_seconds"], 5
            )
            self.assertEqual(
                item["cached"]["metrics"]["reduction_seconds"]["max_seconds"], 7
            )
            for key in (
                "pure_ibp_seconds",
                "terminal_substitution_seconds",
                "ibp_excluding_owner_seconds",
            ):
                self.assertIsNone(item["cached"]["metrics"][key])
            self.assertIsNone(item["cached"]["stages"]["frontend_canonicalization"])

    def test_fresh_release_without_window_stays_functional(self):
        self.paired()
        item = self.evaluate(qualified=False)["quadratic-four-rung-uv"]
        self.assertEqual(item["status"], "functional-only")
        self.assertIsNone(item["cached"])
        self.assertFalse(
            any(row["performance_eligible"] for row in item["functional_records"])
        )

    def test_concurrent_functional_launcher_cannot_be_promoted(self):
        self.paired()
        path = Path(self.records[0]["launcher"]["path"])
        launcher = merge.read(path)
        launcher["functional_only_not_isolated_three_run_performance"] = True
        path.write_text(json.dumps(launcher))
        self.records[0]["launcher"] = merge.evidence(path)
        self.window_rows[0]["launcher_sha256"] = merge.digest(path)
        with self.assertRaisesRegex(ValueError, "functional launcher"):
            self.evaluate()

    def test_two_pairs_do_not_establish_performance(self):
        self.paired(count=2)
        with self.assertRaisesRegex(ValueError, "three matched"):
            self.evaluate()

    def test_exact_result_guard_and_input_files_are_required(self):
        self.paired()
        for name in (
            "result.symbolica",
            "nonzero-conditions.symbolica",
            "vakint-input.symbolica",
            "unknown-masters.json",
        ):
            with self.subTest(name=name):
                file = Path(self.records[0]["record"]["path"]).parent / name
                original = file.read_bytes()
                file.write_text("CORRUPTED TEST FIXTURE")
                with self.assertRaises((ValueError, json.JSONDecodeError)):
                    self.evaluate()
                file.write_bytes(original)

    def test_missing_sidecar_fails(self):
        self.paired()
        (Path(self.records[0]["record"]["path"]).parent / "result.symbolica").unlink()
        with self.assertRaises(FileNotFoundError):
            self.evaluate()

    def test_coherent_guard_change_is_still_parity_failure(self):
        self.paired()
        file = (
            Path(self.records[0]["record"]["path"]).parent
            / "nonzero-conditions.symbolica"
        )
        file.write_text("DIFFERENT_TEST_GUARD\n")
        self.change_report(0, lambda r: r["guards"].update(sha256=merge.digest(file)))
        with self.assertRaisesRegex(ValueError, "parity differs"):
            self.evaluate()

    def test_missing_preparation_guards_rejected(self):
        self.paired()
        self.change_report(
            0, lambda r: r["guards"].update(includes_preparation_conditions=False)
        )
        with self.assertRaisesRegex(ValueError, "complete primary guards"):
            self.evaluate()

    def test_retained_unknown_is_functional_not_closed_performance(self):
        self.paired()
        file = Path(self.records[0]["record"]["path"]).parent / "unknown-masters.json"
        file.write_text('["TEST_FAMILY::Unknown(1)"]')

        def edit(report):
            report["catalogue_closure"].update(
                closed=False, unknown_count=1, unknown_sha256=merge.digest(file)
            )
            report["reduction_completed"] = False

        self.change_report(0, edit)
        item = self.evaluate(qualified=False)["quadratic-four-rung-uv"]
        self.assertEqual(item["functional_records"][0]["status"], "retained-incomplete")
        with self.assertRaisesRegex(ValueError, "catalogue is incomplete"):
            self.evaluate()

    def test_window_requires_acceptance_drain_and_nonoverlap(self):
        self.paired()
        for key in (
            "accepted",
            "serialized",
            "all_processes_drained",
            "source_unchanged",
            "final_performance_eligible",
        ):
            with (
                self.subTest(key=key),
                self.assertRaisesRegex(ValueError, "serialized performance window"),
            ):
                self.evaluate(window_edit=lambda w, key=key: w.update({key: False}))
        with self.assertRaisesRegex(ValueError, "overlap"):
            self.evaluate(
                window_edit=lambda w: w["samples"][1].update(started_unix_seconds=110)
            )

    def test_window_binds_sample_artifact_and_launcher(self):
        self.paired()
        for key in (
            "artifact_sha256",
            "launcher_sha256",
            "case_id",
            "variant",
            "sample",
        ):
            with (
                self.subTest(key=key),
                self.assertRaisesRegex(ValueError, "binding differs"),
            ):
                self.evaluate(
                    window_edit=lambda w, key=key: w["samples"][0].update(
                        {key: "WRONG_TEST_ID"}
                    )
                )

    def test_timing_rejects_nonfinite_negative_and_boolean(self):
        self.paired()
        for value in (float("nan"), float("inf"), -1, True):
            with self.subTest(value=value):
                self.change_report(
                    0, lambda r, value=value: r.update(reduction_seconds=value)
                )
                with self.assertRaisesRegex(ValueError, "timing"):
                    self.evaluate()

    def test_baseline_and_representation_must_match(self):
        self.paired()
        self.change_report(
            0, lambda r: r["persistent_cache_attachments_before"].update({"1": True})
        )
        with self.assertRaisesRegex(ValueError, "cache policy"):
            self.evaluate()
        self.change_report(
            0, lambda r: r["persistent_cache_attachments_before"].update({"1": False})
        )
        self.change_report(
            0, lambda r: r.update(coefficient_representation="factorized")
        )
        with self.assertRaisesRegex(ValueError, "parity differs"):
            self.evaluate()

    def test_archived_q0_retains_its_original_cohort(self):
        artifact = {
            "native": {"base_revision": "a" * 40},
            "symbolica_revision": "b" * 40,
            "native_profile": {"opt_level": 2},
            "driver_source_sha256": "c" * 64,
            "input_sha256": "d" * 64,
            "data": {"revision": "e" * 40},
            "binary_sha256": "f" * 64,
        }
        artifact_path = self.write("old-q0-artifact.json", artifact)
        launch = self.write(
            "old-q0-launch.json",
            {
                "schema": "vakint.q0-full-weighted-launch.v1",
                "bindings": [
                    {"sha256": merge.digest(artifact_path)},
                    {"sha256": "f" * 64},
                ],
            },
        )
        report = {
            "schema": "vakint.q0-full-weighted-native-diagnostic.v2",
            "compiled_cohort": artifact,
            "input_sha256": "d" * 64,
            "status": "TEST_FUNCTIONAL_ONLY",
        }
        path = self.write("old-q0-report.json", report)
        self.records.append(
            {
                "case_id": "zero-momentum-diagnostic",
                "record": merge.evidence(path),
                "launcher": merge.evidence(launch),
                "artifact": merge.evidence(artifact_path),
            }
        )
        item = self.evaluate(qualified=False)["zero-momentum-diagnostic"]
        self.assertEqual(item["functional_records"][0]["raw_report"], report)
        self.assertIsNone(item["uncached"])
        with self.assertRaisesRegex(ValueError, "exploratory native q0"):
            self.evaluate()

    def test_separately_accepted_overlapping_windows_cannot_be_combined(self):
        self.paired()
        index = self.index()
        entries = merge.read(index)
        window = merge.read(self.root / "window.json")
        window["test_fixture_note"] = "SECOND WINDOW SAME INTERVALS"
        second = self.write("second-window.json", window)
        entries["records"][0]["window"] = merge.evidence(second)
        index.write_text(json.dumps(entries))
        with self.assertRaisesRegex(ValueError, "parity differs"):
            merge.primary_records(merge.ReportEvidence(index))

    def test_multiple_driver_artifact_shape_uses_same_binary_identity(self):
        self.paired()
        original = merge.read(self.artifact)
        del original["binary"]
        del original["actual_codegen"]
        original["artifacts"] = [
            {
                "sha256": merge.digest(self.binary),
                "normal_cli": True,
                "immutable_path": str(self.binary),
                "profile": {"opt_level": "3", "test": False},
                "actual_rustc": {"lto": "fat"},
            }
        ]
        self.artifact.write_text(json.dumps(original))
        for record, row in zip(self.records, self.window_rows):
            record["artifact"] = merge.evidence(self.artifact)
            row["artifact_sha256"] = merge.digest(self.artifact)
        self.assertEqual(
            self.evaluate()["quadratic-four-rung-uv"]["status"], "measured"
        )

    def synthetic_suite(self):
        # Minimal synthetic tuning suite. These are file-shape controls only.
        cohort = dict(
            zip(
                (
                    "main_revision",
                    "native_revision",
                    "kernel_revision",
                    "main_source_snapshot_sha256",
                    "native_source_snapshot_sha256",
                ),
                list(self.stamps.values())[:5],
            )
        )
        identity = {
            "suite": "tuning",
            "repetitions": 3,
            "cache_sha256": merge.digest(self.cache),
            "data_manifest_sha256": "7" * 64,
            "workload_sha256": "8" * 64,
            "source_census_sha256": "9" * 64,
            "owner_proof": {"compiled_cohort": cohort},
            "build_git_revision": "1" * 40,
            "build_profile": "release",
            "rayon_num_threads": "1",
            "cpu_affinity": "TEST_CPU",
            "load_prepared_caches": False,
            "executable_sha256": merge.digest(self.binary),
        }
        build = {
            **identity,
            "build_identity": {"build": {"git_revision": "1" * 40}},
            "preparation": {"options": {"single_file": {"binding": "TEST_BINDING"}}},
        }
        samples = []
        closure = {
            "closed": False,
            "unknown_count": 1,
            "unknown_masters": ["TEST_UNKNOWN"],
        }
        for index in range(3):
            for variant in ("uncached", "cached"):
                self.sample(variant, index)
                directory = self.root / f"case-000-{index:02}-{variant}"
                directory.mkdir()
                (directory / "value.symbolica").write_text("TEST_RESULT")
                (directory / "guards.json").write_text('["TEST_GUARD"]')
                (directory / "unknown-masters.json").write_text('["TEST_UNKNOWN"]')
                sample = {
                    "case_id": "quadratic-four-rung-uv",
                    "variant": variant,
                    "sample": index,
                    "identity": identity,
                    "success": True,
                    "worker_success": True,
                    "fresh_process": True,
                    "native_keys_authenticated": True,
                    "owner_proof": identity["owner_proof"],
                    "input_sha256": "f" * 64,
                    "supervision": {
                        "success": True,
                        "exit_code": 0,
                        "deadline_exhausted": False,
                    },
                    "frozen_protected_roots": [[1] * 15],
                    "actual_lowered_roots": [[1] * 15],
                    "frozen_lowered_roots": [[1] * 15],
                    "persistent_cache_attachments_before": {
                        str(i): False for i in range(1, 6)
                    },
                    "persistent_cache_attachments_during_evaluation": {
                        str(i): i == 5 and variant == "cached" for i in range(1, 6)
                    },
                    "cache_statistics": {
                        "binding": "TEST_BINDING",
                        "persistent": {
                            "requested_root_hits": 0,
                            "intermediate_hits": 1,
                        },
                    },
                    "catalogue_closure": {
                        **closure,
                        "unknown_sha256": merge.digest(
                            directory / "unknown-masters.json"
                        ),
                    },
                    "value_sha256": merge.digest(directory / "value.symbolica"),
                    "guards_sha256": merge.digest(directory / "guards.json"),
                    "stages": {},
                }
                sample.update(
                    {key: 2 if variant == "uncached" else 1 for key in merge.METRICS}
                )
                observations = {
                    "child_pid": 2000 + len(samples),
                    "supervisor_pid": 1000,
                    "started_unix_seconds": self.window_rows[-1][
                        "started_unix_seconds"
                    ],
                    "finished_unix_seconds": self.window_rows[-1][
                        "finished_unix_seconds"
                    ],
                    "interval_scope": "supervisor pre-spawn setup through supervisor-observed worker completion",
                }
                sample["supervision"].update(observations)
                launcher_path = Path(self.records[-1]["launcher"]["path"])
                launcher = merge.read(launcher_path)
                launcher.update(observations)
                launcher_path.write_text(json.dumps(launcher))
                self.records[-1]["launcher"] = merge.evidence(launcher_path)
                self.window_rows[-1].update(observations)
                self.window_rows[-1]["launcher_sha256"] = merge.digest(launcher_path)
                saved = self.write(
                    str(directory.relative_to(self.root) / "measurement.json"), sample
                )
                self.records[-1]["record"] = merge.evidence(saved)
                self.window_rows[-1]["report_sha256"] = merge.digest(saved)
                samples.append(sample)
        case = {
            "id": "quadratic-four-rung-uv",
            "successful": True,
            "exact_output_and_guard_parity": True,
            "native_keys_authenticated": True,
            "catalogue_closure": closure,
            "requested_root_hits": 0,
            "intermediate_hits": 3,
            "backing_reuse": True,
            "uncached": merge.summarize([2] * 3),
            "cached": merge.summarize([1] * 3),
            "uncached_over_cached_ratio": 2,
        }
        suite = {
            "schema": "vakint.five-loop-cache-suite.v1",
            "identity": identity,
            "complete": True,
            "provisional": False,
            "native_keys_authenticated": True,
            "execution_provenance_complete": True,
            "selection_allowed": True,
            "measurements": samples,
            "cases": [case],
        }
        path = self.write("suite.json", suite)
        build_path = self.write("suite-build.json", build)
        return (
            path,
            suite,
            {
                "raw_report": build,
                "evidence": merge.evidence(build_path),
                "cache_file": merge.evidence(self.cache),
            },
            closure,
        )

    def test_suite_checks_closure_cohort_and_window(self):
        path, suite, build, closure = self.synthetic_suite()
        evidence = merge.ReportEvidence(self.index())
        accepted = merge.suite_record(path, "tuning", build, evidence)
        self.assertFalse(accepted["rows"][0]["catalogue_closure"]["closed"])
        suite["cases"][0]["catalogue_closure"] = {
            "closed": True,
            "unknown_count": 0,
            "unknown_masters": [],
        }
        path.write_text(json.dumps(suite))
        with self.assertRaisesRegex(ValueError, "closure differs"):
            merge.suite_record(path, "tuning", build, evidence)
        del suite["cases"][0]["catalogue_closure"]
        path.write_text(json.dumps(suite))
        with self.assertRaisesRegex(ValueError, "report catalogue closure"):
            merge.suite_record(path, "tuning", build, evidence)
        suite["cases"][0]["catalogue_closure"] = closure
        path.write_text(json.dumps(suite))
        for bundle in evidence.records.values():
            bundle["artifact"]["stamps"]["VAKINT_CACHE_BUILD_RUSTRED_REVISION"] = (
                "a" * 40
            )
        with self.assertRaisesRegex(ValueError, "compiled cohort differs"):
            merge.suite_record(path, "tuning", build, evidence)

    def suite_with_hit_counts(self, roots=0, others=1):
        path, suite, producer, _ = self.synthetic_suite()
        suite["cases"][0].update(
            requested_root_hits=3 * roots,
            intermediate_hits=3 * others,
            backing_reuse="intermediate-only",
        )
        for sample, reference, window in zip(
            suite["measurements"], self.records, self.window_rows
        ):
            if sample["variant"] == "cached":
                sample["cache_statistics"]["persistent"].update(
                    requested_root_hits=roots, intermediate_hits=others
                )
            sample_path = Path(reference["record"]["path"])
            sample_path.write_text(json.dumps(sample))
            reference["record"] = merge.evidence(sample_path)
            window["report_sha256"] = merge.digest(sample_path)
        path.write_text(json.dumps(suite))
        return path, suite, producer

    def test_zero_route_root_hits_do_not_prove_intermediate_only_reuse(self):
        path, suite, producer = self.suite_with_hit_counts()
        result = merge.suite_record(
            path, "tuning", producer, merge.ReportEvidence(self.index())
        )
        row = result["rows"][0]
        self.assertEqual(row["requested_root_hits"], 0)
        self.assertEqual(row["intermediate_hits"], 3)
        self.assertEqual(row["backing_reuse"], "intermediate-only")
        self.assertEqual(result["raw_suite"], suite)
        self.assertEqual(
            row["cache_hit_scope"],
            {
                "some_exact_requested_root_reuse_proven": False,
                "other_hits": "unclassified",
            },
        )

    def test_positive_route_root_hits_prove_some_but_not_all_root_reuse(self):
        path, _, producer = self.suite_with_hit_counts(roots=1)
        result = merge.suite_record(
            path, "tuning", producer, merge.ReportEvidence(self.index())
        )
        self.assertEqual(
            result["rows"][0]["cache_hit_scope"],
            {
                "some_exact_requested_root_reuse_proven": True,
                "other_hits": "unclassified",
            },
        )

    def test_no_hits_remain_no_observed_reuse_despite_reported_label(self):
        path, _, producer = self.suite_with_hit_counts(others=0)
        result = merge.suite_record(
            path, "tuning", producer, merge.ReportEvidence(self.index())
        )
        self.assertEqual(
            result["rows"][0]["cache_hit_scope"],
            {
                "some_exact_requested_root_reuse_proven": False,
                "other_hits": "none-observed",
            },
        )

    def test_unknown_revision_or_claim_cannot_authenticate_finer_hit_scope(self):
        path, suite, producer = self.suite_with_hit_counts()
        # The fixture's Native revision is deliberately synthetic, not an
        # allowlisted release. Unvalidated future-looking labels are not proof.
        suite["cases"][0].update(
            cache_hit_scope={"other_hits": "proper-descendants-only"},
            initial_frontier_authenticated=True,
            backing_reuse="proper-descendants-only",
        )
        path.write_text(json.dumps(suite))
        result = merge.suite_record(
            path, "tuning", producer, merge.ReportEvidence(self.index())
        )
        self.assertEqual(result["raw_suite"], suite)
        self.assertEqual(
            result["rows"][0]["cache_hit_scope"]["other_hits"], "unclassified"
        )

    def test_noninteger_root_counter_cannot_attest_exact_root_reuse(self):
        path, suite, producer = self.suite_with_hit_counts(roots=1)
        suite["cases"][0]["requested_root_hits"] = 3.0
        path.write_text(json.dumps(suite))
        with self.assertRaisesRegex(ValueError, "cache-hit aggregate differs"):
            merge.suite_record(
                path, "tuning", producer, merge.ReportEvidence(self.index())
            )

    def test_suite_window_cannot_invent_supervisor_timestamps_or_process_ids(self):
        path, _, producer, _ = self.synthetic_suite()
        for field in (
            "child_pid",
            "supervisor_pid",
            "started_unix_seconds",
            "finished_unix_seconds",
        ):
            with self.subTest(field=field):
                index = self.index(
                    window_edit=lambda window, field=field: window["samples"][0].update(
                        {field: window["samples"][0][field] + 1}
                    )
                )
                with self.assertRaisesRegex(
                    ValueError, "actual supervisor observations"
                ):
                    merge.suite_record(
                        path, "tuning", producer, merge.ReportEvidence(index)
                    )
                self.window_rows[0][field] -= 1

    def test_suite_without_recorded_supervisor_facts_cannot_qualify(self):
        path, suite, producer, _ = self.synthetic_suite()
        first = suite["measurements"][0]
        del first["supervision"]["child_pid"]
        sample_path = Path(self.records[0]["record"]["path"])
        sample_path.write_text(json.dumps(first))
        self.records[0]["record"] = merge.evidence(sample_path)
        self.window_rows[0]["report_sha256"] = merge.digest(sample_path)
        path.write_text(json.dumps(suite))
        with self.assertRaisesRegex(ValueError, "supervisor process identity"):
            merge.suite_record(
                path, "tuning", producer, merge.ReportEvidence(self.index())
            )

    def cross_cohort_fixture(self):
        binding = {
            "program_sha256": [1] * 32,
            "family_fingerprint": "TEST_FAMILY",
            "physical_arity": 15,
        }
        old_cohort = {field: self.stamps[stamp] for field, stamp in merge.COHORT_STAMPS}
        old_cohort["semantic_binding"] = binding
        new_cohort = {
            **old_cohort,
            "main_revision": "a" * 40,
            "main_source_snapshot_sha256": "b" * 64,
        }
        artifact = merge.read(self.artifact)
        artifact["stamps"].update(
            {stamp: new_cohort[field] for field, stamp in merge.COHORT_STAMPS}
        )
        self.artifact.write_text(json.dumps(artifact))

        def root(index):
            return {
                "id": f"TEST_ROOT_{index}",
                "powers": [index] + [1] * 14,
                "public_input": None if index == 0 else f"TEST_INPUT_{index}",
                "rationale": "SYNTHETIC FILE CONTRACT ONLY",
                "protected_roots": [],
                "lowered_roots": [],
            }

        selection = {
            "schema": "vakint.five-loop-cache-workload.v1",
            "data_manifest_sha256": "7" * 64,
            "source_census_sha256": "9" * 64,
            "training": [root(0)],
            "tuning": [root(1)],
            "held_out": [root(i) for i in range(2, 14)],
        }
        refs, workloads, proofs = {}, {}, {}
        for role, cohort in (("producer", old_cohort), ("consumer", new_cohort)):
            proof_path = self.write(
                f"{role}-proof.json",
                {
                    "schema": "vakint.public-lowering-owner-map-proof.v1",
                    "compiled_cohort": cohort,
                    "semantic_binding": binding,
                    "data_manifest_sha256": "7" * 64,
                    "measurement_executable_sha256": merge.digest(self.binary),
                },
            )
            proof_identity = {**merge.evidence(proof_path), "compiled_cohort": cohort}
            workload = {**selection, "owner_proof": merge.evidence(proof_path)}
            workload_path = self.write(f"{role}-workload.json", workload)
            refs[role] = {
                "owner_proof": merge.evidence(proof_path),
                "compiled_cohort": merge.evidence(
                    self.write(f"{role}-cohort.json", cohort)
                ),
                "workload" if role == "producer" else "final_workload": merge.evidence(
                    workload_path
                ),
            }
            workloads[role], proofs[role] = workload_path, proof_identity
        identity = {
            "cache_sha256": merge.digest(self.cache),
            "data_manifest_sha256": "7" * 64,
            "source_census_sha256": "9" * 64,
            "workload_sha256": merge.digest(workloads["consumer"]),
            "owner_proof": proofs["consumer"],
            "build_git_revision": new_cohort["main_revision"],
            "executable_sha256": merge.digest(self.binary),
        }
        build = {
            **identity,
            "workload_sha256": merge.digest(workloads["producer"]),
            "owner_proof": proofs["producer"],
            "build_identity": {"build": {"git_revision": "1" * 40}},
            "preparation": {"options": {"single_file": {"binding": binding}}},
        }
        build_path = self.write("cross-producer-build.json", build)
        producer = {
            "raw_report": build,
            "evidence": merge.evidence(build_path),
            "cache_file": merge.evidence(self.cache),
        }
        refs["producer"].update(
            cache=merge.evidence(self.cache), build_report=merge.evidence(build_path)
        )
        refs["consumer"].update(
            collector_artifacts=merge.evidence(self.artifact),
            qualification_workload=merge.evidence(
                self.write(
                    "bootstrap-workload.json", {**selection, "owner_proof": None}
                )
            ),
        )
        for name in (
            "protected_union",
            "training_route_qualification",
            "alias_exclusion",
        ):
            refs["consumer"][name] = merge.evidence(
                self.write(name + ".json", {"TEST_RECEIPT_ONLY": name})
            )
        accepted = {
            "schema": "vakint.cache-consumer-requalification.v1",
            "accepted": True,
            "test_fixture": True,
            **refs,
            "semantic_binding": binding,
            "immutable_selection": merge.selected_inputs(selection),
            "checks": {
                key: True
                for key in (
                    "canonical_consumer_proof",
                    "unchanged_selection",
                    "twelve_distinct_held_out_sets",
                    "training_raw_and_nonzero_owner_disjoint",
                    "cache_import_validated",
                )
            },
        }
        acceptance = self.write("requalification.json", accepted)
        index = self.index(qualified=False)
        document = merge.read(index)
        document["cache_requalification"] = merge.evidence(acceptance)
        index.write_text(json.dumps(document))
        return identity, producer, acceptance, index

    def refresh_acceptance(self, path, index, document):
        path.write_text(json.dumps(document))
        evidence = merge.read(index)
        evidence["cache_requalification"] = merge.evidence(path)
        index.write_text(json.dumps(evidence))

    def test_new_consumer_suite_imports_original_cache_with_accepted_bridge(self):
        suite_path, suite, _, _ = self.synthetic_suite()
        identity, producer, accepted_path, _ = self.cross_cohort_fixture()
        suite["identity"].update(identity)
        binding = producer["raw_report"]["preparation"]["options"]["single_file"][
            "binding"
        ]
        for sample, reference, window in zip(
            suite["measurements"], self.records, self.window_rows
        ):
            sample["identity"] = suite["identity"]
            sample["owner_proof"] = identity["owner_proof"]
            sample["cache_statistics"]["binding"] = binding
            sample_path = Path(reference["record"]["path"])
            sample_path.write_text(json.dumps(sample))
            reference["record"] = merge.evidence(sample_path)
            reference["artifact"] = merge.evidence(self.artifact)
            window["report_sha256"] = merge.digest(sample_path)
            window["artifact_sha256"] = merge.digest(self.artifact)
        suite_path.write_text(json.dumps(suite))
        index = self.index()
        document = merge.read(index)
        document["cache_requalification"] = merge.evidence(accepted_path)
        index.write_text(json.dumps(document))
        result = merge.suite_record(
            suite_path, "tuning", producer, merge.ReportEvidence(index)
        )
        self.assertEqual(
            result["cache_requalification"]["evidence"], merge.evidence(accepted_path)
        )
        self.assertEqual(result["identity"]["build_git_revision"], "a" * 40)
        self.assertEqual(
            producer["raw_report"]["build_identity"]["build"]["git_revision"], "1" * 40
        )
        self.assertFalse(result["rows"][0]["catalogue_closure"]["closed"])

    def test_cross_cohort_accepted_without_relabeling_producer_or_parsing_cache(self):
        identity, producer, path, index = self.cross_cohort_fixture()
        old = copy.deepcopy(producer)
        result = merge.ReportEvidence(index).cache_qualification(identity, producer)
        self.assertEqual(result["evidence"], merge.evidence(path))
        self.assertNotEqual(
            identity["build_git_revision"],
            producer["raw_report"]["build_identity"]["build"]["git_revision"],
        )
        self.assertNotEqual(
            identity["owner_proof"], producer["raw_report"]["owner_proof"]
        )
        self.assertEqual(producer, old)
        with self.assertRaises(json.JSONDecodeError):
            merge.read(
                self.cache
            )  # Successful validation above hashed these opaque bytes.

    def test_cross_cohort_requires_explicit_accepted_qualification(self):
        identity, producer, path, index = self.cross_cohort_fixture()
        original = merge.read(path)
        for value in (False, None, "true"):
            document = {**original, "accepted": value}
            self.refresh_acceptance(path, index, document)
            with self.assertRaisesRegex(ValueError, "no accepted"):
                merge.ReportEvidence(index).cache_qualification(identity, producer)
        self.refresh_acceptance(path, index, original)
        evidence = merge.read(index)
        del evidence["cache_requalification"]
        index.write_text(json.dumps(evidence))
        with self.assertRaisesRegex(ValueError, "missing hashed evidence"):
            merge.ReportEvidence(index).cache_qualification(identity, producer)

    def test_cross_cohort_cannot_relabel_original_producer_revision(self):
        identity, producer, path, index = self.cross_cohort_fixture()
        document = merge.read(path)
        producer["raw_report"]["build_identity"]["build"]["git_revision"] = "d" * 40
        build_path = Path(producer["evidence"]["path"])
        build_path.write_text(json.dumps(producer["raw_report"]))
        producer["evidence"] = merge.evidence(build_path)
        document["producer"]["build_report"] = merge.evidence(build_path)
        self.refresh_acceptance(path, index, document)
        with self.assertRaisesRegex(ValueError, "recorded revision differs"):
            merge.ReportEvidence(index).cache_qualification(identity, producer)

    def test_cross_cohort_tampered_qualification_receipt_rejected(self):
        identity, producer, path, index = self.cross_cohort_fixture()
        document = merge.read(path)
        Path(document["consumer"]["alias_exclusion"]["path"]).write_text("{}")
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            merge.ReportEvidence(index).cache_qualification(identity, producer)

    def test_cross_cohort_changed_selection_rejected_even_with_new_hash(self):
        identity, producer, path, index = self.cross_cohort_fixture()
        document = merge.read(path)
        workload_path = Path(document["consumer"]["final_workload"]["path"])
        workload = merge.read(workload_path)
        workload["held_out"][0]["powers"][0] = 999
        workload_path.write_text(json.dumps(workload))
        document["consumer"]["final_workload"] = merge.evidence(workload_path)
        identity["workload_sha256"] = merge.digest(workload_path)
        self.refresh_acceptance(path, index, document)
        with self.assertRaisesRegex(ValueError, "frozen selection"):
            merge.ReportEvidence(index).cache_qualification(identity, producer)

    def test_cross_cohort_binding_stamps_data_and_alias_verdict_stay_strict(self):
        identity, producer, path, index = self.cross_cohort_fixture()
        original = merge.read(path)
        changed = copy.deepcopy(original)
        changed["semantic_binding"]["program_sha256"][0] = 2
        self.refresh_acceptance(path, index, changed)
        with self.assertRaisesRegex(ValueError, "native binding"):
            merge.ReportEvidence(index).cache_qualification(identity, producer)
        changed = copy.deepcopy(original)
        changed["checks"]["training_raw_and_nonzero_owner_disjoint"] = False
        self.refresh_acceptance(path, index, changed)
        with self.assertRaisesRegex(ValueError, "no accepted"):
            merge.ReportEvidence(index).cache_qualification(identity, producer)
        self.refresh_acceptance(path, index, original)
        with self.assertRaisesRegex(ValueError, "data or census"):
            merge.ReportEvidence(index).cache_qualification(
                {**identity, "data_manifest_sha256": "a" * 64}, producer
            )
        artifact = merge.read(self.artifact)
        artifact["stamps"]["VAKINT_CACHE_BUILD_SYMBOLICA_REVISION"] = "c" * 40
        self.artifact.write_text(json.dumps(artifact))
        changed = copy.deepcopy(original)
        changed["consumer"]["collector_artifacts"] = merge.evidence(self.artifact)
        self.refresh_acceptance(path, index, changed)
        with self.assertRaisesRegex(ValueError, "consumer artifact/stamps"):
            merge.ReportEvidence(index).cache_qualification(identity, producer)

    def prepared_fixture(self, case="quadratic-four-rung-uv"):
        selected = [i for i, row in enumerate(self.records) if row["case_id"] == case]
        first = merge.read(self.records[selected[0]]["record"]["path"])
        generator_ref = merge.read(self.records[selected[0]]["launcher"]["path"])[
            "generator_acceptance"
        ]
        generator = merge.read(generator_ref["path"])
        directory = self.root / ("SYNTHETIC_PREPARED_" + case)
        directory.mkdir()
        input_path = directory / "vakint-input.symbolica"
        input_path.write_text("TEST_INPUT")
        guard_path = directory / "projector-conditions.symbolica"
        guard_path.write_text("\n".join(first["prepared_nonzero_conditions"]))
        timings = {
            name: 1.0
            for name in (
                "graph_generation_seconds",
                "ingress_seconds",
                "scalar_contraction_seconds",
                "taylor_seconds",
                "angular_seconds",
                "total_seconds",
            )
        }
        manifest = {
            "schema": "vakint.qcd-ladder.prepared-input.v1",
            "definition": "massless-qcd-four-rung-native-uv-v1",
            "projection": merge.PRIMARY[case],
            "scientific_sources_sha256": {
                name: value
                for name, value in first["fixture_sources_sha256"].items()
                if name != "mod.rs"
            },
            "scientific_assembly_sha256": "7" * 64,
            "producer": {
                "fixture_sources_sha256": first["fixture_sources_sha256"],
                "source_commit": generator["source_commit"],
                "executable_sha256": generator["immutable_elf_sha256"],
            },
            "input": {
                "bytes": input_path.stat().st_size,
                "sha256": merge.digest(input_path),
            },
            "nonzero_conditions_file": {
                "bytes": guard_path.stat().st_size,
                "sha256": merge.digest(guard_path),
            },
            "nonzero_conditions": first["prepared_nonzero_conditions"],
            "contributions": first["contributions"],
            "scalar_terms": first["scalar_terms"],
            "preparation_seconds": 999.0,
            "timings": timings,
        }
        path = self.write(
            str(directory.relative_to(self.root) / "prepared-input.json"), manifest
        )
        for index in selected:

            def edit(report):
                report.update(
                    input_mode="prepared",
                    fresh_native_preparation=False,
                    input_loading_seconds=0.25,
                    preparation=None,
                    preparation_seconds=None,
                    prepared_input={
                        "manifest_path": str(path),
                        "manifest_sha256": merge.digest(path),
                        "manifest": manifest,
                    },
                )
                report["consumer_fixture_sources_sha256"] = {
                    **report["fixture_sources_sha256"],
                    "mod.rs": "b" * 64,
                }

            self.change_report(index, edit)
            (
                Path(self.records[index]["record"]["path"]).parent
                / "vakint-input.symbolica"
            ).unlink()
            launch_path = Path(self.records[index]["launcher"]["path"])
            launcher = merge.read(launch_path)
            launcher.update(
                input_mode="prepared", prepared_input_manifest_sha256=merge.digest(path)
            )
            launch_path.write_text(json.dumps(launcher))
            self.records[index]["launcher"] = merge.evidence(launch_path)
            self.window_rows[index]["launcher_sha256"] = merge.digest(launch_path)
        return path

    def test_prepared_scope_loads_once_per_sample_without_charging_historical_preparation(
        self,
    ):
        self.paired()
        manifest = self.prepared_fixture()
        row = self.evaluate()["quadratic-four-rung-uv"]
        self.assertEqual(row["input_mode"], "prepared")
        self.assertEqual(
            row["prepared_input"]["manifest_sha256"], merge.digest(manifest)
        )
        self.assertEqual(row["prepared_input"]["manifest"]["preparation_seconds"], 999)
        for variant in ("cached", "uncached"):
            self.assertIsNone(row[variant]["metrics"]["preparation_seconds"])
            self.assertEqual(
                row[variant]["metrics"]["input_loading_seconds"]["median_seconds"], 0.25
            )
            self.assertEqual(
                row[variant]["metrics"]["end_to_end_seconds"]["median_seconds"], 26
            )
            self.assertTrue(
                all(value is None for value in row[variant]["preparation"].values())
            )
            self.assertIsNone(row[variant]["metrics"]["pure_ibp_seconds"])

    def test_prepared_historical_work_cannot_be_reported_as_current(self):
        self.paired()
        self.prepared_fixture()
        self.change_report(0, lambda row: row.update(preparation_seconds=999))
        with self.assertRaisesRegex(ValueError, "historical preparation"):
            self.evaluate()

    def test_prepared_payload_and_projector_guards_are_rehashed(self):
        self.paired()
        path = self.prepared_fixture()
        for name in ("vakint-input.symbolica", "projector-conditions.symbolica"):
            with self.subTest(name=name):
                payload = path.parent / name
                original = payload.read_bytes()
                payload.write_text("CORRUPTED_TEST_PAYLOAD")
                with self.assertRaisesRegex(ValueError, "payload changed"):
                    self.evaluate()
                payload.write_bytes(original)

    def test_prepared_modes_cannot_mix_fresh_and_reused_input(self):
        self.paired()
        path = self.prepared_fixture()
        manifest = merge.read(path)
        self.change_report(
            0,
            lambda row: row.update(
                input_mode="fresh",
                fresh_native_preparation=True,
                preparation_seconds=manifest["preparation_seconds"],
                preparation=manifest["timings"],
                input_loading_seconds=None,
            ),
        )
        launch_path = Path(self.records[0]["launcher"]["path"])
        launcher = merge.read(launch_path)
        launcher["input_mode"] = "fresh"
        launch_path.write_text(json.dumps(launcher))
        self.records[0]["launcher"] = merge.evidence(launch_path)
        self.window_rows[0]["launcher_sha256"] = merge.digest(launch_path)
        with self.assertRaisesRegex(ValueError, "parity differs"):
            self.evaluate()

    def test_mutually_equal_forged_prepared_input_is_not_the_accepted_physical_fixture(
        self,
    ):
        self.paired()
        path = self.prepared_fixture()
        manifest = merge.read(path)
        payload = path.parent / "vakint-input.symbolica"
        payload.write_text("DIFFERENT_CANONICAL_LOOKING_TEST_INPUT")
        manifest["input"] = {
            "bytes": payload.stat().st_size,
            "sha256": merge.digest(payload),
        }
        path.write_text(json.dumps(manifest))
        for index in range(len(self.records)):
            self.change_report(
                index,
                lambda row: row.update(
                    input_sha256=manifest["input"]["sha256"],
                    input_bytes=manifest["input"]["bytes"],
                    prepared_input={
                        "manifest_path": str(path),
                        "manifest_sha256": merge.digest(path),
                        "manifest": manifest,
                    },
                ),
            )
            launch_path = Path(self.records[index]["launcher"]["path"])
            launcher = merge.read(launch_path)
            launcher["prepared_input_manifest_sha256"] = merge.digest(path)
            launch_path.write_text(json.dumps(launcher))
            self.records[index]["launcher"] = merge.evidence(launch_path)
            self.window_rows[index]["launcher_sha256"] = merge.digest(launch_path)
        with self.assertRaisesRegex(ValueError, "accepted complete physical fixture"):
            self.evaluate()

    def test_prepared_manifest_and_generator_must_be_pinned_by_launcher(self):
        self.paired()
        self.prepared_fixture()
        launch_path = Path(self.records[0]["launcher"]["path"])
        original = merge.read(launch_path)
        for field in ("prepared_input_manifest_sha256", "generator_acceptance"):
            with self.subTest(field=field):
                launcher = {
                    key: value for key, value in original.items() if key != field
                }
                launch_path.write_text(json.dumps(launcher))
                self.records[0]["launcher"] = merge.evidence(launch_path)
                self.window_rows[0]["launcher_sha256"] = merge.digest(launch_path)
                with self.assertRaises(ValueError):
                    self.evaluate()

    def test_changed_hashed_artifact_rejected(self):
        self.paired()
        index = self.index()
        self.artifact.write_text("{}")
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            merge.ReportEvidence(index)

    def strict_build(self):
        return self.write(
            "build.json",
            {
                "schema": "vakint.five-loop-cache-build.v2",
                "phase": "final",
                "complete": True,
                "publishable": True,
                "test_fixture": True,
                "elapsed_seconds": 12,
                "preparation": {
                    "status": "complete",
                    "total_roots": 1,
                    "completed_roots": 1,
                    "native_profile": {"strict": True, "unresolved_policy": "reject"},
                },
                "supervision": {"success": True, "hard_deadline_exhausted": False},
                "cache_sha256": merge.digest(self.cache),
                "compressed_bytes": self.cache.stat().st_size,
                "data_manifest_sha256": "7" * 64,
                "within_50_mib_target": True,
            },
        )

    def selected_pilot_fixture(self):
        identity, producer, bridge_path, index_path = self.cross_cohort_fixture()
        strict = merge.read(self.strict_build())
        build = copy.deepcopy(producer["raw_report"])
        build.update(
            {key: value for key, value in strict.items() if key != "preparation"}
        )
        build["phase"] = "pilot"
        build["preparation"].update(strict["preparation"])
        build_path = Path(producer["evidence"]["path"])
        build_path.write_text(json.dumps(build))
        bridge = merge.read(bridge_path)
        bridge["producer"]["build_report"] = merge.evidence(build_path)
        self.refresh_acceptance(bridge_path, index_path, bridge)
        candidate_path = self.write(
            "candidate-acceptance.json",
            {
                "schema": "vakint.strict-complete-small-candidate-acceptance.v1",
                "state": "STRICT_COMPLETE_CANDIDATE_NOT_FINAL_SELECTION",
                "utc": "2026-10-23T12:00:00Z",
                "test_fixture": True,
                "reports": {"build": merge.evidence(build_path)},
                "cache_file": {
                    **merge.evidence(self.cache),
                    "compressed_bytes": self.cache.stat().st_size,
                },
                "workload_sha256": build["workload_sha256"],
                "selected_training_roots": build["preparation"]["total_roots"],
                "completed_strict_roots": build["preparation"]["completed_roots"],
            },
        )
        selection_path = self.write(
            "final-selection.json",
            {
                "schema": "vakint.cache-final-selection.v1",
                "accepted": True,
                "decision": "retain_complete_candidate",
                "reason": "SYNTHETIC explicit final selection of this complete unchanged candidate",
                "decided_utc": "2026-10-24T12:00:00Z",
                "producer_phase": "pilot",
                "test_fixture": True,
                "cache": merge.evidence(self.cache),
                "build_report": merge.evidence(build_path),
                "producer_workload": bridge["producer"]["workload"],
                "consumer_workload": bridge["consumer"]["final_workload"],
                "candidate_acceptance": merge.evidence(candidate_path),
                "cache_requalification": merge.evidence(bridge_path),
            },
        )
        index = merge.read(index_path)
        index["final_selection"] = merge.evidence(selection_path)
        index_path.write_text(json.dumps(index))
        return {
            "identity": identity,
            "build": build_path,
            "bridge": bridge_path,
            "index": index_path,
            "candidate": candidate_path,
            "selection": selection_path,
        }

    def selected_pilot(self, fixture):
        return merge.build_record(
            fixture["build"], self.cache, merge.ReportEvidence(fixture["index"])
        )

    def change_pilot_selection(self, fixture, edit):
        selection = merge.read(fixture["selection"])
        edit(selection)
        fixture["selection"].write_text(json.dumps(selection))
        index = merge.read(fixture["index"])
        index["final_selection"] = merge.evidence(fixture["selection"])
        fixture["index"].write_text(json.dumps(index))

    def test_explicit_final_selection_preserves_complete_pilot_provenance(self):
        fixture = self.selected_pilot_fixture()
        originals = {
            name: fixture[name].read_bytes()
            for name in ("build", "bridge", "candidate", "selection")
        }
        result = self.selected_pilot(fixture)
        self.assertEqual(result["raw_report"], json.loads(originals["build"]))
        self.assertEqual(result["raw_report"]["phase"], "pilot")
        self.assertEqual(
            result["final_selection"]["evidence"], merge.evidence(fixture["selection"])
        )
        self.assertEqual(result["cache_file"], merge.evidence(self.cache))
        for name, contents in originals.items():
            self.assertEqual(fixture[name].read_bytes(), contents)

    def test_unselected_pilot_stays_ineligible_with_or_without_provenance(self):
        fixture = self.selected_pilot_fixture()
        with self.assertRaises(ValueError):
            merge.build_record(fixture["build"], self.cache)
        index = merge.read(fixture["index"])
        del index["final_selection"]
        fixture["index"].write_text(json.dumps(index))
        with self.assertRaises(ValueError):
            self.selected_pilot(fixture)

    def test_original_final_build_needs_no_pilot_selection(self):
        path = self.strict_build()
        original = path.read_bytes()
        record = merge.build_record(path, self.cache)
        self.assertEqual(record["raw_report"], json.loads(original))
        self.assertEqual(record["raw_report"]["phase"], "final")
        self.assertIsNone(record.get("final_selection"))

    def test_final_selection_decision_flags_reason_and_utc_are_explicit(self):
        fixture = self.selected_pilot_fixture()
        original = merge.read(fixture["selection"])
        for field, value in (
            ("accepted", False),
            ("accepted", "true"),
            ("decision", "automatic"),
            ("producer_phase", "final"),
            ("reason", ""),
            ("reason", "   "),
            ("decided_utc", None),
            ("decided_utc", "not-a-date"),
            ("decided_utc", "2026-10-24T12:00:00"),
            ("decided_utc", "2026-10-24T12:00:00+02:00"),
        ):
            with self.subTest(field=field, value=value):
                self.change_pilot_selection(
                    fixture,
                    lambda selection, original=original: (
                        selection.clear(),
                        selection.update(copy.deepcopy(original)),
                    ),
                )
                self.change_pilot_selection(
                    fixture,
                    lambda selection, field=field, value=value: selection.update(
                        {field: value}
                    ),
                )
                with self.assertRaises((ValueError, TypeError)):
                    self.selected_pilot(fixture)

    def test_final_selection_reference_bytes_cannot_change_silently(self):
        fixture = self.selected_pilot_fixture()
        paths = [fixture[name] for name in ("selection", "candidate", "bridge")] + [
            self.cache
        ]
        for path in paths:
            original = path.read_bytes()
            with self.subTest(path=path.name):
                path.write_bytes(b"ALTERED")
                with self.assertRaises(ValueError):
                    self.selected_pilot(fixture)
                path.write_bytes(original)

    def test_final_selection_cannot_point_to_different_cache_or_workloads(self):
        fixture = self.selected_pilot_fixture()
        original = merge.read(fixture["selection"])
        for field in ("cache", "producer_workload", "consumer_workload"):
            replacement = self.write(
                f"different-{field}.json", {"SYNTHETIC_DIFFERENT": field}
            )
            with self.subTest(field=field):
                self.change_pilot_selection(
                    fixture,
                    lambda selection, original=original: (
                        selection.clear(),
                        selection.update(copy.deepcopy(original)),
                    ),
                )
                self.change_pilot_selection(
                    fixture,
                    lambda selection, field=field, replacement=replacement: (
                        selection.update({field: merge.evidence(replacement)})
                    ),
                )
                with self.assertRaises((ValueError, KeyError)):
                    self.selected_pilot(fixture)

    def test_final_selection_and_index_must_pin_the_same_accepted_bridge(self):
        fixture = self.selected_pilot_fixture()
        bridge = merge.read(fixture["bridge"])
        alternate = self.write(
            "different-bridge.json", {**bridge, "different_receipt": True}
        )
        self.change_pilot_selection(
            fixture,
            lambda selection: selection.update(
                cache_requalification=merge.evidence(alternate)
            ),
        )
        with self.assertRaises(ValueError):
            self.selected_pilot(fixture)

    def test_selected_pilot_bridge_acceptance_cannot_be_revoked_or_omitted(self):
        fixture = self.selected_pilot_fixture()
        bridge = merge.read(fixture["bridge"])
        bridge["accepted"] = False
        self.refresh_acceptance(fixture["bridge"], fixture["index"], bridge)
        self.change_pilot_selection(
            fixture,
            lambda selection: selection.update(
                cache_requalification=merge.evidence(fixture["bridge"])
            ),
        )
        with self.assertRaisesRegex(ValueError, "no accepted requalification"):
            self.selected_pilot(fixture)
        index = merge.read(fixture["index"])
        del index["cache_requalification"]
        fixture["index"].write_text(json.dumps(index))
        with self.assertRaises((ValueError, KeyError)):
            self.selected_pilot(fixture)

    def test_pilot_never_uses_same_cohort_shortcut_without_accepted_bridge(self):
        fixture = self.selected_pilot_fixture()
        build = merge.read(fixture["build"])
        identity = {
            **build,
            "build_git_revision": build["build_identity"]["build"]["git_revision"],
        }
        producer = {
            "raw_report": build,
            "evidence": merge.evidence(fixture["build"]),
            "cache_file": merge.evidence(self.cache),
        }
        index = merge.read(fixture["index"])
        del index["cache_requalification"]
        fixture["index"].write_text(json.dumps(index))
        with self.assertRaisesRegex(ValueError, "missing hashed evidence"):
            merge.ReportEvidence(fixture["index"]).cache_qualification(
                identity, producer
            )

    def test_selected_pilot_still_rejects_partial_failed_or_nonstrict_builds(self):
        fixture = self.selected_pilot_fixture()
        original = fixture["build"].read_bytes()
        for change in (
            lambda build: build.update(complete=False),
            lambda build: build.update(publishable=False),
            lambda build: build["preparation"].update(status="partial"),
            lambda build: build["preparation"].update(completed_roots=0),
            lambda build: build["preparation"].update(failure="synthetic failure"),
            lambda build: build["supervision"].update(success=False),
            lambda build: build["supervision"].update(hard_deadline_exhausted=True),
            lambda build: build["preparation"]["native_profile"].update(strict=False),
            lambda build: build["preparation"]["native_profile"].update(
                unresolved_policy="preserve"
            ),
        ):
            with self.subTest(change=change):
                build = json.loads(original)
                change(build)
                fixture["build"].write_text(json.dumps(build))
                with self.assertRaises(ValueError):
                    self.selected_pilot(fixture)

    def test_candidate_acceptance_counts_and_original_build_binding_are_checked(self):
        fixture = self.selected_pilot_fixture()
        original = fixture["candidate"].read_bytes()
        for change in (
            lambda candidate: candidate.update(schema="other-schema"),
            lambda candidate: candidate.update(state="INCOMPLETE"),
            lambda candidate: candidate.update(selected_training_roots=2),
            lambda candidate: candidate.update(completed_strict_roots=0),
            lambda candidate: candidate.update(workload_sha256="0" * 64),
            lambda candidate: candidate["reports"]["build"].update(sha256="0" * 64),
        ):
            with self.subTest(change=change):
                candidate = json.loads(original)
                change(candidate)
                fixture["candidate"].write_text(json.dumps(candidate))
                self.change_pilot_selection(
                    fixture,
                    lambda selection: selection.update(
                        candidate_acceptance=merge.evidence(fixture["candidate"])
                    ),
                )
                with self.assertRaises(ValueError):
                    self.selected_pilot(fixture)

    def test_selected_pilot_append_preserves_raw_phase_and_all_history(self):
        fixture = self.selected_pilot_fixture()
        historical = {
            "selected_reduction_run": "OLD_4L",
            "reduction_runs": {"OLD_4L": {"raw": [1, 2]}},
            "selected_five_loop_run": "OLD_5L",
            "five_loop_runs": {"OLD_5L": {"phase": "pilot", "raw": [3, 4]}},
        }
        base = self.write("pilot-history.json", historical)
        output = self.root / "appended-selected-pilot.json"
        original_build = fixture["build"].read_bytes()
        args = [
            "adapter",
            "--base",
            str(base),
            "--build",
            str(fixture["build"]),
            "--cache",
            str(self.cache),
            "--evidence",
            str(fixture["index"]),
            "--output",
            str(output),
        ]
        with patch("sys.argv", args):
            merge.main()
        result = merge.read(output)
        self.assertEqual(merge.read(base), historical)
        self.assertEqual(
            result["selected_reduction_run"], historical["selected_reduction_run"]
        )
        self.assertEqual(
            result["selected_five_loop_run"], historical["selected_five_loop_run"]
        )
        self.assertEqual(result["reduction_runs"], historical["reduction_runs"])
        self.assertEqual(
            result["five_loop_runs"]["OLD_5L"], historical["five_loop_runs"]["OLD_5L"]
        )
        new = next(
            value for key, value in result["five_loop_runs"].items() if key != "OLD_5L"
        )
        self.assertEqual(new["cache"]["raw_report"], json.loads(original_build))
        self.assertEqual(new["cache"]["raw_report"]["phase"], "pilot")
        self.assertIsNotNone(new["cache"]["final_selection"])
        self.assertEqual(fixture["build"].read_bytes(), original_build)

    def test_primary_performance_requires_matching_completed_cache(self):
        self.paired()
        base = self.write("base.json", {"test_fixture": True})
        output = self.root / "absent.json"
        args = [
            "adapter",
            "--base",
            str(base),
            "--evidence",
            str(self.index()),
            "--output",
            str(output),
        ]
        with (
            patch("sys.argv", args),
            self.assertRaisesRegex(ValueError, "completed strict cache"),
        ):
            merge.main()
        self.assertFalse(output.exists())
        build = self.strict_build()
        document = merge.read(build)
        document["data_manifest_sha256"] = "a" * 64
        build.write_text(json.dumps(document))
        args += ["--build", str(build), "--cache", str(self.cache)]
        with (
            patch("sys.argv", args),
            self.assertRaisesRegex(ValueError, "cache/data differs"),
        ):
            merge.main()
        self.assertFalse(output.exists())
        document["data_manifest_sha256"] = "7" * 64
        document["preparation"]["completed_roots"] = 0
        build.write_text(json.dumps(document))
        with (
            patch("sys.argv", args),
            self.assertRaisesRegex(ValueError, "roots are incomplete"),
        ):
            merge.main()
        self.assertFalse(output.exists())

    def test_append_preserves_all_historical_values_and_rejects_overwrite(self):
        self.paired()
        base_value = {
            "test_fixture": True,
            "historical": [1, {"raw": "unchanged"}],
            "selected_reduction_run": "old4L",
            "reduction_runs": {"old4L": {"data": 12}},
            "five_loop_runs": {"old5L": {"raw": "unchanged"}},
        }
        base = self.write("base.json", base_value)
        output = self.root / "new-report.json"
        args = [
            "adapter",
            "--base",
            str(base),
            "--evidence",
            str(self.index()),
            "--build",
            str(self.strict_build()),
            "--cache",
            str(self.cache),
            "--output",
            str(output),
            "--select-current",
        ]
        with patch("sys.argv", args):
            merge.main()
        result = merge.read(output)
        for key, value in base_value.items():
            if key == "five_loop_runs":
                self.assertEqual(result[key]["old5L"], value["old5L"])
            else:
                self.assertEqual(result[key], value)
        self.assertEqual(merge.read(base), base_value)
        with patch("sys.argv", args), self.assertRaisesRegex(ValueError, "new path"):
            merge.main()
        args[2] = str(output)
        args[args.index("--output") + 1] = str(self.root / "duplicate.json")
        with (
            patch("sys.argv", args),
            self.assertRaisesRegex(ValueError, "already appended"),
        ):
            merge.main()
        self.assertFalse((self.root / "duplicate.json").exists())

    def warm_pairs(self, count=3):
        """Synthetic files only; the two calls are nested in one process record."""

        def external(path):
            return {**merge.evidence(path), "bytes": path.stat().st_size}

        def sidecar(path):
            return {
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": merge.digest(path),
            }

        self.warm_external = external
        self.warm_sidecar = sidecar
        dot = "SYNTHETIC CONNECTED FIVE-LOOP GRAPH"
        powers = [1, 2, 1, 1, 1, 1, 1, 1, 1]
        manifests = {
            name: external(self.write("data/" + name, {"test_fixture": name}))
            for name in ("manifest.json", *(f"{i}L/manifest.json" for i in range(1, 6)))
        }
        application = external(
            self.write("application.json", {"test_fixture": "native application"})
        )
        old = self.write(
            "science.json",
            {
                "test_fixture": True,
                "pinched_dot": dot,
                "pinched_powers": powers,
                "pinched_loop_blocks": [5],
                "exact_equality": True,
                "data_manifest_sha256": manifests["manifest.json"]["sha256"],
            },
        )
        fixture = self.write(
            "fixture.json",
            {
                "test_fixture": True,
                "source_report": str(old),
                "source_report_sha256": merge.digest(old),
                "dot_field": "pinched_dot",
                "dot_sha256": merge.hashlib.sha256(dot.encode()).hexdigest(),
                "dot_bytes": len(dot),
                "powers_native_order": powers,
                "historical_pinched_loop_blocks": [5],
            },
        )
        sources = {
            "crates/vakint/examples/five_loop_saved_application.rs": "a" * 64,
            "crates/vakint/examples/saved_application/dot.rs": "b" * 64,
            "crates/vakint/tests/fixtures/phi_3_4.json": "c" * 64,
            "Cargo.lock": "d" * 64,
        }
        snapshot = self.write(
            "source.json",
            {
                "head": self.stamps["VAKINT_CACHE_BUILD_GIT_REVISION"],
                "source_sha256": sources,
            },
        )
        self.stamps["VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256"] = merge.digest(
            snapshot
        )
        self.stamps["VAKINT_CACHE_BUILD_RUSTC"] = "rustc TEST"
        artifact = merge.read(self.artifact)
        artifact["stamps"] = self.stamps
        self.artifact.write_text(json.dumps(artifact))
        cache = {"persistent": {str(i): {"attached": i < 5} for i in range(1, 6)}}
        for n in range(count):
            directory = self.root / f"pair-{n}"
            directory.mkdir()
            dot_path = directory / "input.dot"
            dot_path.write_text(dot)
            input_path = directory / "input.symbolica"
            input_path.write_text("SYNTHETIC SAME INPUT")
            inputs = {
                "schema": "vakint.five-loop-connected-timed-pair.v1",
                "test_fixture": True,
                "powers_native_order": powers,
                "mass_squared": "1",
                "numerator": "1",
                "dot": sidecar(dot_path),
                "dot_source": external(dot_path),
                "input": sidecar(input_path),
                "data_manifests": manifests,
                "application": application,
                "build_rustc": self.stamps["VAKINT_CACHE_BUILD_RUSTC"],
                "build_profile": "release",
                "build_git_revision": self.stamps["VAKINT_CACHE_BUILD_GIT_REVISION"],
                "build_rustred_revision": self.stamps[
                    "VAKINT_CACHE_BUILD_RUSTRED_REVISION"
                ],
                "build_symbolica_revision": self.stamps[
                    "VAKINT_CACHE_BUILD_SYMBOLICA_REVISION"
                ],
                "build_main_source_snapshot_sha256": merge.digest(snapshot),
                "build_native_source_snapshot_sha256": self.stamps[
                    "VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256"
                ],
                "compiled_example_sha256": sources[
                    "crates/vakint/examples/five_loop_saved_application.rs"
                ],
                "compiled_dot_adapter_sha256": sources[
                    "crates/vakint/examples/saved_application/dot.rs"
                ],
                "model_sha256": sources["crates/vakint/tests/fixtures/phi_3_4.json"],
                "compiled_cargo_lock_sha256": sources["Cargo.lock"],
                "executable": external(self.binary),
                "cache_statistics_before_cold": cache,
                "settings_debug": "TEST SAME SETTINGS",
            }
            input_manifest = self.write(f"pair-{n}/input.json", inputs)
            phases = {}
            for phase, duration in (("cold", 10 + n), ("warm", 2 + n)):
                value = directory / (phase + ".symbolica")
                value.write_text("SYNTHETIC SAME RESULT")
                guards = self.write(f"pair-{n}/{phase}-guards.json", ["G1", "G2"])
                families = self.write(
                    f"pair-{n}/{phase}-families.json", [["F1", {"TEST": "family"}]]
                )
                path = self.write(
                    f"pair-{n}/{phase}.json",
                    {
                        "phase": phase,
                        "completed": True,
                        "whole_public_application_seconds": duration,
                        "value": sidecar(value),
                        "ordered_nonzero_conditions": sidecar(guards),
                        "families": sidecar(families),
                        "nonzero_condition_count": 2,
                        "family_count": 1,
                        "cache_statistics_after": cache,
                    },
                )
                phases[phase + "_record"] = sidecar(path)
            pair = self.write(
                f"pair-{n}/pair.json",
                {
                    "schema": "vakint.five-loop-connected-timed-pair.v1",
                    "test_fixture": True,
                    "completed": True,
                    "exact_value_parity": True,
                    "ordered_full_guard_parity": True,
                    "family_identity_parity": True,
                    "additional_demo_checks_run": False,
                    **phases,
                },
            )
            observed = {
                "child_pid": 1000 + n,
                "supervisor_pid": 900,
                "child_start_ticks": 100 + n,
                "started_unix_seconds": 100 + 100 * n,
                "finished_unix_seconds": 150 + 100 * n,
                "interval_scope": "supervisor pre-spawn setup through supervisor-observed worker completion",
            }
            launcher = self.write(
                f"pair-{n}/launcher.json",
                {
                    "test_fixture": True,
                    "final_performance_eligible": True,
                    "fresh_process": True,
                    "binary_sha256": merge.digest(self.binary),
                    "case_id": merge.WARM_REUSE_CASE,
                    "variant": "cold-warm-pair",
                    "sample": n,
                    "cpu_affinity": "TEST_CPU",
                    "thread_environment": self.threads,
                    **observed,
                    "input_manifest": merge.evidence(input_manifest),
                    "scientific_fixture": merge.evidence(fixture),
                    "source_snapshot": merge.evidence(snapshot),
                },
            )
            self.records.append(
                {
                    "record": merge.evidence(pair),
                    "launcher": merge.evidence(launcher),
                    "artifact": merge.evidence(self.artifact),
                    "case_id": merge.WARM_REUSE_CASE,
                }
            )
            self.window_rows.append(
                {
                    "report_sha256": merge.digest(pair),
                    "launcher_sha256": merge.digest(launcher),
                    "artifact_sha256": merge.digest(self.artifact),
                    "case_id": merge.WARM_REUSE_CASE,
                    "variant": "cold-warm-pair",
                    "sample": n,
                    "exit_code": 0,
                    "fresh_process": True,
                    **observed,
                }
            )

    def warm_evaluate(self, **kwargs):
        return merge.warm_reuse_records(merge.ReportEvidence(self.index(**kwargs)))

    def warm_refresh(self, n=0):
        directory = self.root / f"pair-{n}"
        pair = merge.read(directory / "pair.json")
        for phase in ("cold", "warm"):
            pair[phase + "_record"] = self.warm_sidecar(directory / (phase + ".json"))
        (directory / "pair.json").write_text(json.dumps(pair))
        launcher = merge.read(directory / "launcher.json")
        launcher["input_manifest"] = merge.evidence(directory / "input.json")
        (directory / "launcher.json").write_text(json.dumps(launcher))
        self.records[n]["record"] = merge.evidence(directory / "pair.json")
        self.records[n]["launcher"] = merge.evidence(directory / "launcher.json")
        self.window_rows[n]["report_sha256"] = merge.digest(directory / "pair.json")
        self.window_rows[n]["launcher_sha256"] = merge.digest(
            directory / "launcher.json"
        )

    def test_warm_reuse_three_fresh_pairs_have_separate_nested_medians(self):
        self.warm_pairs()
        result = self.warm_evaluate()
        self.assertEqual(result["status"], "measured-paired-reuse")
        self.assertEqual(result["first_call"]["samples_seconds"], [10, 11, 12])
        self.assertEqual(result["warm_call"]["samples_seconds"], [2, 3, 4])
        self.assertEqual(result["warm_over_first"], 3 / 11)
        self.assertIsNone(result["catalogue_closure"])
        self.assertIsNone(result["pure_ibp_seconds"])
        self.assertIsNone(result["owner_preparation_seconds"])
        self.assertEqual(len(result["functional_records"]), 3)
        self.assertTrue(
            all(x["performance_eligible"] for x in result["functional_records"])
        )

    def test_warm_reuse_without_accepted_window_stays_functional(self):
        self.warm_pairs()
        result = self.warm_evaluate(qualified=False)
        self.assertEqual(result["status"], "functional-only")
        self.assertIsNone(result["first_call"])
        self.assertIsNone(result["warm_call"])
        with self.assertRaises(ValueError):
            self.warm_evaluate(window_edit=lambda w: w.update(accepted=False))

    def test_warm_reuse_requires_three_independent_process_pairs(self):
        self.warm_pairs(2)
        with self.assertRaises(ValueError):
            self.warm_evaluate()
        self.records.clear()
        self.window_rows.clear()
        for n in range(2):
            import shutil

            shutil.rmtree(self.root / f"pair-{n}")
        self.warm_pairs()
        launcher = self.root / "pair-1/launcher.json"
        value = merge.read(launcher)
        value.update(child_pid=1000, child_start_ticks=100)
        launcher.write_text(json.dumps(value))
        self.warm_refresh(1)
        self.window_rows[1].update(child_pid=1000, child_start_ticks=100)
        with self.assertRaisesRegex(ValueError, "reuses a process"):
            self.warm_evaluate()

    def test_warm_reuse_checks_actual_launcher_window_process_facts(self):
        self.warm_pairs()
        for field in (
            "child_pid",
            "supervisor_pid",
            "child_start_ticks",
            "started_unix_seconds",
            "finished_unix_seconds",
            "interval_scope",
        ):
            with self.subTest(field=field):
                old = self.window_rows[0][field]
                self.window_rows[0][field] = (
                    old + 1 if isinstance(old, (int, float)) else "WRONG"
                )
                with self.assertRaises(ValueError):
                    self.warm_evaluate()
                self.window_rows[0][field] = old

    def test_warm_reuse_checks_every_cache_attachment_snapshot(self):
        self.warm_pairs()
        directory = self.root / "pair-0"
        for file, field in (
            ("input.json", "cache_statistics_before_cold"),
            ("cold.json", "cache_statistics_after"),
            ("warm.json", "cache_statistics_after"),
        ):
            path = directory / file
            original = path.read_bytes()
            for loop in ("1", "5"):
                with self.subTest(file=file, loop=loop):
                    value = json.loads(original)
                    value[field]["persistent"][loop]["attached"] = loop == "5"
                    path.write_text(json.dumps(value))
                    self.warm_refresh()
                    with self.assertRaisesRegex(ValueError, "attachments differ"):
                        self.warm_evaluate()
            path.write_bytes(original)
            self.warm_refresh()

    def test_warm_reuse_rehashes_input_and_complete_phase_files(self):
        self.warm_pairs()
        directory = self.root / "pair-0"
        for name in (
            "input.json",
            "input.dot",
            "input.symbolica",
            "cold.json",
            "warm.json",
            "cold.symbolica",
            "warm.symbolica",
            "cold-guards.json",
            "warm-guards.json",
            "cold-families.json",
            "warm-families.json",
        ):
            path = directory / name
            original = path.read_bytes()
            with self.subTest(file=name):
                path.write_bytes(b"CORRUPTED")
                with self.assertRaises((ValueError, KeyError)):
                    self.warm_evaluate()
                path.write_bytes(original)

    def test_warm_reuse_does_not_trust_rehashed_ordered_parity_flags(self):
        self.warm_pairs()
        directory = self.root / "pair-0"
        for field, file, value in (
            ("ordered_nonzero_conditions", "warm-guards.json", ["G2", "G1"]),
            ("families", "warm-families.json", [["DIFFERENT", {}]]),
        ):
            path = directory / file
            old = path.read_bytes()
            phase = directory / "warm.json"
            old_phase = phase.read_bytes()
            with self.subTest(field=field):
                path.write_text(json.dumps(value))
                record = json.loads(old_phase)
                record[field] = self.warm_sidecar(path)
                phase.write_text(json.dumps(record))
                self.warm_refresh()
                with self.assertRaisesRegex(ValueError, "parity failed"):
                    self.warm_evaluate()
            path.write_bytes(old)
            phase.write_bytes(old_phase)
            self.warm_refresh()

    def test_warm_reuse_rejects_incomplete_calls_or_extra_demo_work(self):
        self.warm_pairs()
        directory = self.root / "pair-0"
        for file, field, value in (
            ("pair.json", "completed", False),
            ("pair.json", "additional_demo_checks_run", True),
            ("cold.json", "completed", False),
            ("warm.json", "whole_public_application_seconds", -1),
        ):
            path = directory / file
            old = path.read_bytes()
            with self.subTest(file=file, field=field):
                document = json.loads(old)
                document[field] = value
                path.write_text(json.dumps(document))
                self.warm_refresh()
                with self.assertRaises(ValueError):
                    self.warm_evaluate()
            path.write_bytes(old)
            self.warm_refresh()

    def test_warm_reuse_rejects_changed_science_source_data_and_settings(self):
        self.warm_pairs()
        path = self.root / "pair-0/input.json"
        old = path.read_bytes()
        for field, value in (
            ("powers_native_order", [1] * 9),
            ("mass_squared", "2"),
            ("compiled_example_sha256", "e" * 64),
            ("build_git_revision", "a" * 40),
            ("settings_debug", "DIFFERENT SETTINGS"),
        ):
            with self.subTest(field=field):
                document = json.loads(old)
                document[field] = value
                path.write_text(json.dumps(document))
                self.warm_refresh()
                with self.assertRaises(ValueError):
                    self.warm_evaluate()
        path.write_bytes(old)
        self.warm_refresh()
        data = self.root / "data/manifest.json"
        data.write_text("CHANGED")
        with self.assertRaises(ValueError):
            self.warm_evaluate()

    def test_warm_reuse_append_preserves_every_old_cold_result(self):
        self.warm_pairs()
        base_value = {
            "test_fixture": True,
            "summary": [{"cold_seconds": 123}],
            "reduction_runs": {"old": {"value": 42}},
        }
        base = self.write("base.json", base_value)
        output = self.root / "new-warm.json"
        args = [
            "adapter",
            "--base",
            str(base),
            "--evidence",
            str(self.index()),
            "--output",
            str(output),
            "--select-current",
        ]
        with patch("sys.argv", args):
            merge.main()
        merged = merge.read(output)
        for key, value in base_value.items():
            self.assertEqual(merged[key], value)
        current = merged["five_loop_runs"][merged["selected_five_loop_run"]]
        self.assertEqual(current["warm_reuse"]["status"], "measured-paired-reuse")
        self.assertTrue(
            all(x["status"] == "pending" for x in current["primary"].values())
        )
        self.assertIsNone(current["held_out"])
        self.assertIsNone(current["tuning"])
        self.assertIsNone(current["cache"])

    def test_warm_reuse_requires_same_outputs_across_independent_pairs(self):
        self.warm_pairs()
        directory = self.root / "pair-1"
        for phase in ("cold", "warm"):
            value = directory / (phase + ".symbolica")
            value.write_text("COHERENTLY DIFFERENT SECOND PAIR")
            path = directory / (phase + ".json")
            record = merge.read(path)
            record["value"] = self.warm_sidecar(value)
            path.write_text(json.dumps(record))
        self.warm_refresh(1)
        with self.assertRaisesRegex(ValueError, "identical-context"):
            self.warm_evaluate()

    def test_warm_reuse_checks_connected_science_and_payload_lengths(self):
        self.warm_pairs()
        science = self.root / "science.json"
        original = science.read_bytes()
        value = json.loads(original)
        value["pinched_loop_blocks"] = [1, 4]
        science.write_text(json.dumps(value))
        fixture = self.root / "fixture.json"
        authority = merge.read(fixture)
        authority["source_report_sha256"] = merge.digest(science)
        fixture.write_text(json.dumps(authority))
        for n in range(3):
            launch = self.root / f"pair-{n}/launcher.json"
            record = merge.read(launch)
            record["scientific_fixture"] = merge.evidence(fixture)
            launch.write_text(json.dumps(record))
            self.warm_refresh(n)
        with self.assertRaisesRegex(ValueError, "genuinely connected"):
            self.warm_evaluate()
        science.write_bytes(original)
        authority["source_report_sha256"] = merge.digest(science)
        fixture.write_text(json.dumps(authority))
        for n in range(3):
            launch = self.root / f"pair-{n}/launcher.json"
            record = merge.read(launch)
            record["scientific_fixture"] = merge.evidence(fixture)
            launch.write_text(json.dumps(record))
            self.warm_refresh(n)
        phase = self.root / "pair-0/cold.json"
        record = merge.read(phase)
        record["value"]["bytes"] += 1
        phase.write_text(json.dumps(record))
        self.warm_refresh()
        with self.assertRaisesRegex(ValueError, "length differs"):
            self.warm_evaluate()

    def test_warm_reuse_uses_existing_release_artifact_gate(self):
        self.warm_pairs()
        artifact = merge.read(self.artifact)
        artifact["actual_codegen"]["fat_lto"] = False
        self.artifact.write_text(json.dumps(artifact))
        for record, row in zip(self.records, self.window_rows):
            record["artifact"] = merge.evidence(self.artifact)
            row["artifact_sha256"] = merge.digest(self.artifact)
        with self.assertRaisesRegex(ValueError, "release executable"):
            self.warm_evaluate()


class CumulativeCacheContractTests(unittest.TestCase):
    """Opaque synthetic indexes and receipts only; no native codec or physics."""

    setUp = ReportContractTests.setUp
    write = ReportContractTests.write
    index = ReportContractTests.index
    sample = ReportContractTests.sample
    paired = ReportContractTests.paired
    synthetic_suite = ReportContractTests.synthetic_suite
    strict_build = ReportContractTests.strict_build

    def cumulative_fixture(self, completed=False):
        cache = self.root / "native-pack"
        cache.mkdir()
        index = cache / "cache.rrindex"
        index.write_text("OPAQUE SYNTHETIC INDEX; NOT A NATIVE CACHE")
        source_workload = self.write(
            "source-workload.json",
            {
                "schema": "vakint.five-loop-cache-workload.v1",
                "training": [{"id": "SOURCE_TRAINING", "powers": [2] * 15}],
            },
        )
        source_producer = self.write(
            "source-producer.json",
            {
                "workload_sha256": merge.digest(source_workload),
                "schema": "vakint.five-loop-cache-build.v2",
                "complete": False,
                "preparation": {
                    "status": "failed",
                    "completed_roots": 1,
                    "total_roots": 3,
                    "native_profile": {"strict": True, "unresolved_policy": "reject"},
                },
                "functional_acceptance": completed,
                "process_zero_acceptance": completed,
                "collector_exit_code": 0 if completed else 1,
            },
        )
        native = {
            "schema": "rustred.running-application-cache.v1",
            "generation": 1,
            "accepted_records": 5,
            "committed_records": 5,
            "pending_records": 0,
            "pending_bytes": 0,
            "writer_error": None,
            "last_error": None,
            "closed": completed,
        }
        capture = self.write(
            "source-capture.json",
            {"report": native, "flush_complete": completed, "close": {"error": None}},
        )
        snapshot = self.write(
            "source-snapshot.json",
            {
                "schema": "rustred.running-application-cache-snapshot.v1",
                "delta": {
                    "generation": 1,
                    "manifest_digest": "NATIVE-DIGEST-NOT-SHA256",
                    "manifest_path": str(index),
                    "manifest_bytes": index.stat().st_size,
                },
            },
        )
        barrier = self.write(
            "source-barrier.json",
            {
                "schema": "vakint.native-cache-snapshot-barrier.v1",
                "accepted": True,
                "committed_snapshot": True,
                "producer_completed": completed,
                "producer": merge.evidence(source_producer),
                "capture": merge.evidence(capture),
                "snapshot": merge.evidence(snapshot),
                "writer_fault": None,
                "generation": 1,
                "manifest_digest": "NATIVE-DIGEST-NOT-SHA256",
                **{
                    k: native[k]
                    for k in (
                        "accepted_records",
                        "committed_records",
                        "pending_records",
                        "pending_bytes",
                    )
                },
            },
        )
        source = {
            "id": "physical-q0",
            "kind": "running-capture",
            "scope": "completed-example" if completed else "committed-prefix",
            "producer": merge.evidence(source_producer),
            "capture": merge.evidence(capture),
            "snapshot": merge.evidence(snapshot),
            "barrier": merge.evidence(barrier),
            "index": merge.evidence(index),
        }
        report = {
            "schema": "vakint.cumulative-native-cache-report.v1",
            "cache_sha256": None,
            "cache_native_index_identity": merge.native_cache_identity(cache),
            "semantic_binding": {
                "family_fingerprint": "TEST_FAMILY",
                "physical_arity": 15,
                "artifact_digest": [1] * 32,
                "layer": "saved-application-route-apply-v1",
            },
            "data_manifest_sha256": "7" * 64,
            "sources": [source],
            "inventory": {
                "compressed_bytes": 20_000_000_000,
                "stored_record_references": 7,
                "physical_chunks": 3,
            },
        }
        validation = self.write(
            "native-validation.json",
            {
                "schema": "vakint.cumulative-native-cache-validation.v1",
                "accepted": True,
                **{
                    key: report[key]
                    for key in (
                        "cache_native_index_identity",
                        "semantic_binding",
                        "data_manifest_sha256",
                        "sources",
                        "inventory",
                    )
                },
                "checks": {
                    key: True
                    for key in (
                        "immutable_index_snapshot",
                        "authenticated_chunks",
                        "strict_completed_records_only",
                        "authenticated_bindings",
                        "committed_source_snapshots",
                    )
                },
            },
        )
        report["validation"] = merge.evidence(validation)
        path = self.write("cumulative-report.json", report)
        return path, cache

    def revise(self, report_path, name, mutate):
        """Rebind enclosing test evidence, so semantic checks see the changed facts."""
        report = merge.read(report_path)
        source = report["sources"][0]
        target = Path(source[name]["path"])
        doc = merge.read(target)
        mutate(doc)
        target.write_text(json.dumps(doc))
        source[name] = merge.evidence(target)
        if name in ("producer", "capture", "snapshot"):
            barrier_path = Path(source["barrier"]["path"])
            barrier = merge.read(barrier_path)
            barrier[name] = merge.evidence(target)
            barrier_path.write_text(json.dumps(barrier))
            source["barrier"] = merge.evidence(barrier_path)
        self.refresh_validation(report_path, report)

    def refresh_validation(self, path, report):
        target = Path(report["validation"]["path"])
        validation = merge.read(target)
        for key in (
            "cache_native_index_identity",
            "semantic_binding",
            "data_manifest_sha256",
            "sources",
            "inventory",
        ):
            validation[key] = report[key]
        target.write_text(json.dumps(validation))
        report["validation"] = merge.evidence(target)
        path.write_text(json.dumps(report))

    def test_prefix_stays_incomplete_without_aggregate_size_cap(self):
        path, cache = self.cumulative_fixture()
        result = merge.build_record(path, cache)
        self.assertFalse(result["physical_examples_complete"])
        self.assertFalse(result["performance_eligible"])
        self.assertNotIn("preparation", result["raw_report"])
        self.assertEqual(
            result["raw_report"]["inventory"]["stored_record_references"], 7
        )
        self.assertFalse(
            result["sources"][0]["receipts"]["producer"]["raw_receipt"]["complete"]
        )
        self.assertEqual(
            result["raw_report"]["inventory"]["compressed_bytes"], 20_000_000_000
        )

    def test_failed_workload_cannot_become_completed_with_acceptance_booleans(self):
        path, cache = self.cumulative_fixture(completed=True)
        with self.assertRaisesRegex(ValueError, "failed finite workload"):
            merge.build_record(path, cache)

    def test_partial_source_cannot_be_relabelled_completed(self):
        path, cache = self.cumulative_fixture()
        report = merge.read(path)
        report["sources"][0]["scope"] = "completed-example"
        self.refresh_validation(path, report)
        with self.assertRaisesRegex(ValueError, "partial source"):
            merge.build_record(path, cache)

    def test_snapshot_barrier_binds_actual_generation_and_counters(self):
        for name, change in (
            ("snapshot", lambda d: d["delta"].update(generation=2)),
            ("capture", lambda d: d["report"].update(committed_records=4)),
            ("barrier", lambda d: d.update(committed_snapshot=False)),
        ):
            with (
                self.subTest(name=name),
                tempfile.TemporaryDirectory(dir=self.root) as unused,
            ):
                old_root, self.root = self.root, Path(unused)
                try:
                    path, cache = self.cumulative_fixture()
                    self.revise(path, name, change)
                    with self.assertRaisesRegex(
                        ValueError, "snapshot barrier|native capture"
                    ):
                        merge.build_record(path, cache)
                finally:
                    self.root = old_root

    def test_changed_native_child_index_is_rejected(self):
        path, cache = self.cumulative_fixture()
        (cache / "cache.rrindex").write_text("CHANGED OPAQUE BYTES")
        with self.assertRaisesRegex(ValueError, "indexes changed"):
            merge.build_record(path, cache)

    def test_wrong_native_binding_is_rejected_even_with_same_indexes(self):
        path, cache = self.cumulative_fixture()
        report = merge.read(path)
        report["semantic_binding"]["physical_arity"] = 14
        path.write_text(json.dumps(report))
        with self.assertRaisesRegex(ValueError, "binding/source inventory"):
            merge.build_record(path, cache)

    def test_native_index_identity_is_location_independent_and_unambiguous(self):
        _, cache = self.cumulative_fixture()
        first = merge.native_cache_identity(cache)
        second = self.root / "other-pack"
        second.mkdir()
        (second / "cache.rrindex").write_bytes((cache / "cache.rrindex").read_bytes())
        self.assertEqual(first, merge.native_cache_identity(second))
        (second / "bundle.rrindex").write_text("AMBIGUOUS")
        with self.assertRaisesRegex(ValueError, "exactly one"):
            merge.native_cache_identity(second)
        (second / "bundle.rrindex").unlink()
        (second / "cache.rrindex").unlink()
        (second / "cache.rrindex").symlink_to(cache / "cache.rrindex")
        with self.assertRaisesRegex(ValueError, "regular file"):
            merge.native_cache_identity(second)

    def test_native_bundle_hashes_all_sorted_child_indexes(self):
        cache = self.root / "native-bundle"
        cache.mkdir()
        (cache / "bundle.rrindex").write_text("OPAQUE BUNDLE")
        for name in ("f" * 64, "0" * 64):
            child = cache / "bindings" / name
            child.mkdir(parents=True)
            (child / "cache.rrindex").write_text(name)
        identity = merge.native_cache_identity(cache)
        self.assertEqual(list(identity["indexes"]), sorted(identity["indexes"]))
        before = identity["sha256"]
        (cache / "bindings" / ("0" * 64) / "cache.rrindex").write_text("CHANGED")
        self.assertNotEqual(before, merge.native_cache_identity(cache)["sha256"])

    def test_old_completed_file_build_has_no_obsolete_size_gate(self):
        path = self.strict_build()
        report = merge.read(path)
        report["within_50_mib_target"] = False
        path.write_text(json.dumps(report))
        self.assertEqual(merge.build_record(path, self.cache)["raw_report"], report)

    def test_cumulative_append_preserves_all_historical_measurements(self):
        path, cache = self.cumulative_fixture()
        historical = {
            "selected_reduction_run": "OLD",
            "reduction_runs": {str(i): {"raw": [i]} for i in range(48)},
            "selected_five_loop_run": "EXISTING",
            "five_loop_runs": {"EXISTING": {"raw": "UNCHANGED"}},
        }
        base = self.write("historical.json", historical)
        original = base.read_bytes()
        output = self.root / "new.json"
        with patch(
            "sys.argv",
            [
                "adapter",
                "--base",
                str(base),
                "--build",
                str(path),
                "--cache",
                str(cache),
                "--evidence",
                str(self.index(False)),
                "--output",
                str(output),
            ],
        ):
            merge.main()
        self.assertEqual(base.read_bytes(), original)
        new = merge.read(output)
        self.assertEqual(new["reduction_runs"], historical["reduction_runs"])
        self.assertEqual(new["selected_reduction_run"], "OLD")
        self.assertEqual(new["selected_five_loop_run"], "EXISTING")
        self.assertEqual(
            new["five_loop_runs"]["EXISTING"], historical["five_loop_runs"]["EXISTING"]
        )

    def test_scoped_cumulative_append_binds_override_and_preserves_48_measurements(
        self,
    ):
        path, cache = self.cumulative_fixture()
        scope = {
            "required": ["zero-momentum-diagnostic"],
            "deferred": ["quadratic-four-rung-uv"],
        }
        authorization = self.write(
            "physical-scope.json",
            {
                "schema": "vakint.benchmark-physical-scope.v1",
                "physical_scope": scope,
                "reason": "Synthetic user scope: quadratic deferred until after delivery.",
            },
        )
        historical = {
            "reduction_runs": {str(i): {"raw": [i, i + 0.25]} for i in range(48)},
            "selected_reduction_run": "OLD",
            "five_loop_runs": {"old": {"raw": "unchanged"}},
            "selected_five_loop_run": "old",
        }
        base = self.write("historical-scoped.json", historical)
        original = base.read_bytes()
        index = self.index(False)
        identifiers = []
        for scoped in (False, True):
            output = self.root / f"scoped-{scoped}.json"
            command = [
                "adapter",
                "--base",
                str(base),
                "--build",
                str(path),
                "--cache",
                str(cache),
                "--evidence",
                str(index),
                "--output",
                str(output),
            ]
            if scoped:
                command += ["--physical-scope", str(authorization)]
            with patch("sys.argv", command):
                merge.main()
            new = merge.read(output)
            self.assertEqual(new["reduction_runs"], historical["reduction_runs"])
            self.assertEqual(new["selected_reduction_run"], "OLD")
            self.assertEqual(new["selected_five_loop_run"], "old")
            self.assertEqual(
                new["five_loop_runs"]["old"], historical["five_loop_runs"]["old"]
            )
            identifier = next(key for key in new["five_loop_runs"] if key != "old")
            identifiers.append(identifier)
            run = new["five_loop_runs"][identifier]
            # A preserved prefix never becomes a completed q0 run through scope alone.
            self.assertFalse(run["cache"]["physical_examples_complete"])
            if scoped:
                self.assertEqual(run["physical_scope"], scope)
                self.assertEqual(run["cache"]["physical_scope"], scope)
                self.assertEqual(
                    run["physical_scope_authorization"]["evidence"],
                    merge.evidence(authorization),
                )
                self.assertIn(merge.digest(authorization), run["notes"][0])
                self.assertEqual(
                    run["primary"]["quadratic-four-rung-uv"]["status"], "deferred"
                )
            else:
                self.assertNotIn("physical_scope", run)
                self.assertNotIn("physical_scope", run["cache"])
        self.assertNotEqual(*identifiers)
        self.assertEqual(base.read_bytes(), original)

    def public_import_fixture(self, producer, cohort):
        """Actual receipt schemas with opaque tiny test files; no import is run."""
        prefix = "synthetic-public-import"
        root = self.root / prefix
        output = root / "run1/imported-native-bundle"
        native = producer["raw_report"]["cache_native_index_identity"]
        files = {"cache/" + name: item for name, item in native["indexes"].items()}
        for name in native["indexes"]:
            target = output / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes((Path(producer["cache_directory"]) / name).read_bytes())
        total = producer["raw_report"]["inventory"]["compressed_bytes"]
        count = producer["raw_report"]["inventory"]["physical_chunks"]
        for index in range(count):
            files[f"cache/chunks/TEST-{index}.rrchunk.gz"] = {
                "bytes": total // count + (total % count if index == 0 else 0),
                "sha256": "f" * 64,
            }
        descriptor = self.write(
            prefix + "/descriptor.json",
            {
                "schema": "vakint.nativepack-transport.v1",
                "loop_count": 5,
                "binding": producer["raw_report"]["semantic_binding"],
                "generation": 1,
                "stored_record_references": producer["raw_report"]["inventory"][
                    "stored_record_references"
                ],
                "files": files,
                "parts": [{"archive_bytes": 42}],
            },
        )
        core = self.write(prefix + "/TEST-NOT-A-CORE", {"TEST_ONLY": "CORE"})
        wheel = self.write(prefix + "/TEST-NOT-A-WHEEL", {"TEST_ONLY": "WHEEL"})
        sources = {
            "vakint": cohort["main_revision"],
            "rustred": cohort["native_revision"],
            "official_kernel": cohort["kernel_revision"],
        }
        data_path = self.root / "public-import-data"
        packet = self.write(
            prefix + "/packet.json",
            {
                "schema": "vakint.public-native-large-prefix-import-source-packet.v1",
                "unit": "SYNTHETIC-TEST-UNIT",
                "whole_seconds": 3600,
                "command": ["TEST-NOT-EXECUTED", "worker.py"],
                "public_sources": sources,
                "descriptor": merge.evidence(descriptor),
                "fresh_destination": str(output),
                "frozen_data": str(data_path),
                "core": merge.evidence(core),
                "wheel": merge.evidence(wheel),
                "pinned_inputs": [
                    merge.evidence(descriptor),
                    merge.evidence(core),
                    merge.evidence(wheel),
                    merge.evidence(data_path / "manifest.json"),
                ],
            },
        )
        go = self.write(
            prefix + "/root-go.json",
            {
                "decision": "GO_ONCE_EXACT_PUBLIC_IMPORT",
                "packet_sha256": merge.digest(packet),
                "unit": "SYNTHETIC-TEST-UNIT",
                "whole_seconds": 3600,
            },
        )
        indexes = {
            name: item for name, item in files.items() if name.endswith(".rrindex")
        }
        namespaces = sum(name.startswith("cache/bindings/") for name in indexes)
        imported = {
            "format": "multipart_manifest",
            "output": str(output),
            "sha256": merge.digest(descriptor),
            "archive_bytes": descriptor.stat().st_size + 42,
            "extracted_bytes": sum(item["bytes"] for item in files.values()),
            "entries": len(files),
            "namespaces": namespaces,
            "copied_bases": namespaces,
            "copied_chunks": count,
            "preserved_record_references": producer["raw_report"]["inventory"][
                "stored_record_references"
            ],
            "compressed_bytes": total,
            "index_bytes": sum(item["bytes"] for item in indexes.values()),
            "transfer_buffer_bytes": 65536,
        }
        host = {
            "sources": sources,
            "core_sha256": merge.digest(core),
            "wheel_sha256": merge.digest(wheel),
            "mapped_core_paths": [str(core)],
        }
        host_path = self.write(prefix + "/run1/actual-host.json", host)
        state = {str(loop): {"attached": False} for loop in range(1, 6)}
        state["5"]["saved_application"] = {
            "binding": producer["raw_report"]["semantic_binding"],
            "read_backing": False,
            "native": {"rule_applications": 0},
            "persistent": {"lookups": 0, "captured_records": 0},
        }
        loaded = {
            "path": str(data_path),
            "version": merge.read(data_path / "manifest.json")["version"],
            "load_prepared_caches": False,
            "loops": [
                {"loop_count": loop, "cache_loaded": False} for loop in range(1, 6)
            ],
        }
        worker = self.write(
            prefix + "/run1/worker-result.json",
            {
                "schema": "vakint.public-large-native-import.v1",
                "state": "PASS",
                "public_sources": sources,
                "actual_host": str(host_path),
                "data_report": loaded,
                "import_report": imported,
                "exact_native_files_verified": len(files),
                "refused_repeat_preserved_full_output_and_cache_state": True,
                "cache_attachment_performed": False,
                "integral_reduction_or_numerical_master_evaluation": False,
                "cache_coefficient_lookup_or_capture": False,
                "repeat_refusal": {
                    "type": "ValueError",
                    "message": "destination already exists",
                },
                "final_attachment_state": state,
            },
        )
        self.write(prefix + "/run1/public-import-report.json", imported)
        self.write(prefix + "/run1/data-load-report.json", loaded)
        self.write(prefix + "/run1/cache-state-before-repeat.json", state)
        terminal_doc = {
            "schema": "vakint.public-native-import-controller.v1",
            "state": "PASS",
            "worker_returncode": 0,
            "packet_sha256": merge.digest(packet),
            "root_go_sha256": merge.digest(go),
            "worker_result": merge.evidence(worker),
            "whole_elapsed_seconds": 1.0,
            "whole_seconds": 3600,
            "output": str(output),
        }
        terminal = self.write(prefix + "/run1/terminal-outcome.json", terminal_doc)
        admission_doc = {
            "packet_sha256": merge.digest(packet),
            "frozen_inputs_and_host_hashes_verified": True,
        }
        admission = self.write(prefix + "/run1/admission.json", admission_doc)
        invocation = self.write(
            prefix + "/run1/invocation-start.json",
            {
                "packet_sha256": merge.digest(packet),
                "command": ["TEST-NOT-EXECUTED", "worker.py"],
                "pid": 22,
                "controller_pid": 11,
            },
        )
        launch = self.write(
            prefix + "/actual-launch-receipt.json",
            {
                "schema": "vakint.public-large-native-import.actual-launch.v1",
                "unit": "SYNTHETIC-TEST-UNIT",
                "service": {
                    "InvocationID": "SYNTHETIC-TEST-INVOCATION",
                    "MainPID": "11",
                },
                "processes": [
                    {"pid": 11, "start_ticks": 110},
                    {"pid": 22, "start_ticks": 220},
                ],
                "admission": admission_doc,
                "actual_host": host,
                "artifact_sha256": {
                    str(path): merge.digest(path)
                    for path in (host_path, admission, invocation)
                },
            },
        )
        journal = self.write(
            prefix + "/actual-invocation-journal.json",
            [
                {
                    "USER_UNIT": "SYNTHETIC-TEST-UNIT",
                    "USER_INVOCATION_ID": "SYNTHETIC-TEST-INVOCATION",
                }
            ],
        )
        observation = self.write(
            prefix + "/actual-terminal-observation.json",
            {
                "schema": "vakint.public-native-import.actual-terminal-observation.v1",
                "pidfd_completion_observed": True,
                "saved_launch_and_journal_identity_joined": True,
                "public_import_terminal_accepted": True,
                "frozen_inputs_unchanged": True,
                "current_service": {
                    "ActiveState": "inactive",
                    "Result": "success",
                    "ExecMainStatus": "0",
                    "InvocationID": "",
                },
                "unit": "SYNTHETIC-TEST-UNIT",
                "invocation": "SYNTHETIC-TEST-INVOCATION",
                "controller_terminal": terminal_doc,
                "artifact_sha256": {
                    str(path): merge.digest(path)
                    for path in (packet, go, launch, journal, terminal, worker)
                },
            },
        )
        return merge.evidence(observation)

    def native_suite_fixture(self, *, data_manifest=None):
        self.records, self.window_rows = [], []
        path, suite, _, _ = self.synthetic_suite()
        suite["cases"][0]["id"] = "native-tuning-synthetic"
        report_path, cache = self.cumulative_fixture()
        data_path = data_manifest or self.write(
            "public-import-data/manifest.json", {"version": "v-synthetic-import-only"}
        )
        source_report = merge.read(report_path)
        source_report["data_manifest_sha256"] = merge.digest(data_path)
        self.refresh_validation(report_path, source_report)
        producer = merge.build_record(report_path, cache)
        build = producer["raw_report"]
        identity = suite["identity"]
        identity["data_manifest_sha256"] = build["data_manifest_sha256"]
        identity["cache_sha256"] = None
        identity["cache_native_index_identity"] = build["cache_native_index_identity"]
        binding = build["semantic_binding"]
        cohort = identity["owner_proof"]["compiled_cohort"]
        cohort["semantic_binding"] = binding
        workload = {
            "source_census_sha256": identity["source_census_sha256"],
            "data_manifest_sha256": identity["data_manifest_sha256"],
            "training": [],
            "tuning": [
                {
                    "id": suite["cases"][0]["id"],
                    "powers": [1] * 15,
                    "public_input": "SYNTHETIC INPUT",
                }
            ],
            "held_out": [],
        }
        proof_path = self.write(
            "consumer-proof.json",
            {
                "compiled_cohort": cohort,
                "semantic_binding": binding,
                "data_manifest_sha256": build["data_manifest_sha256"],
                "measurement_executable_sha256": merge.digest(self.binary),
            },
        )
        identity["owner_proof"] = {
            **merge.evidence(proof_path),
            "compiled_cohort": cohort,
        }
        workload["owner_proof"] = merge.evidence(proof_path)
        workload_path = self.write("consumer-workload.json", workload)
        identity["workload_sha256"] = merge.digest(workload_path)
        accepted = {
            "schema": "vakint.cumulative-cache-consumer-requalification.v1",
            "accepted": True,
            "producer": merge.evidence(report_path),
            "cache_native_index_identity": build["cache_native_index_identity"],
            "semantic_binding": binding,
            "data_manifest_sha256": build["data_manifest_sha256"],
            "checks": {
                key: True
                for key in (
                    "canonical_consumer_proof",
                    "cache_import_validated",
                    "source_capture_inventory_authenticated",
                    "unchanged_selection",
                    "training_raw_and_nonzero_owner_disjoint",
                    "twelve_distinct_held_out_sets",
                )
            },
            "consumer": {
                "collector_artifacts": merge.evidence(self.artifact),
                "compiled_cohort": merge.evidence(
                    self.write("consumer-cohort.json", cohort)
                ),
                "final_workload": merge.evidence(workload_path),
                "owner_proof": merge.evidence(proof_path),
                **{
                    name: merge.evidence(
                        self.write(name + ".json", {"TEST_FILE_CONTRACT_ONLY": name})
                    )
                    for name in (
                        "protected_union",
                        "training_route_qualification",
                        "alias_exclusion",
                    )
                },
            },
            "immutable_selection": merge.selected_inputs(workload),
            "case_exposure": {suite["cases"][0]["id"]: "tuning"},
        }
        sources = build["sources"]
        accepted["source_starting_keys"] = merge.evidence(
            self.write(
                "source-starting-keys.json",
                {
                    "schema": "vakint.cumulative-source-starting-keys.v1",
                    "sources": [
                        {
                            "id": sources[0]["id"],
                            "producer": sources[0]["producer"],
                            "kind": "selected-workload",
                            "inputs": merge.evidence(
                                self.root / "source-workload.json"
                            ),
                        }
                    ],
                },
            )
        )
        protected = self.write(
            "protected_union.json",
            {
                "schema": "vakint.protected-raw-owner-union.v1",
                "owner_proof": merge.evidence(proof_path),
                "compiled_cohort": accepted["consumer"]["compiled_cohort"],
                "raw_lowered_keys": [],
                "direct_owner_keys": [],
                "declared_case_powers": [],
                "forbidden_direct_target_or_owner_keys": [],
            },
        )
        routes = self.write(
            "training_route_qualification.json",
            {
                "schema": "vakint.source-candidate-owner-qualification.v2",
                "all_routes_succeeded": True,
                "selected_rows_complete": True,
                "compiled_cohort": cohort,
                "semantic_binding": binding,
                "rows": [
                    {
                        "raw_powers": [2] * 15,
                        "route_success": True,
                        "route_error": None,
                        "owner_endpoints": [{"coefficient": "1", "powers": [3] * 15}],
                    }
                ],
            },
        )
        aliases = self.write(
            "alias_exclusion.json",
            {
                "schema": "vakint.source-candidate-external-protected-alias-exclusion.v1",
                "qualification": merge.evidence(routes),
                "protected_union": merge.evidence(protected),
                "rows": [
                    {
                        "powers": [2] * 15,
                        "route_complete_nonempty": True,
                        "raw_protected_intersection": False,
                        "owner_protected_intersections": [],
                        "eligible_after_external_alias_exclusion": True,
                    }
                ],
            },
        )
        accepted["consumer"].update(
            protected_union=merge.evidence(protected),
            training_route_qualification=merge.evidence(routes),
            alias_exclusion=merge.evidence(aliases),
        )
        accepted["cache_import"] = self.public_import_fixture(producer, cohort)
        acceptance = self.write("native-consumer-acceptance.json", accepted)
        for sample, reference, window in zip(
            suite["measurements"], self.records, self.window_rows
        ):
            sample["case_id"] = suite["cases"][0]["id"]
            reference["case_id"] = sample["case_id"]
            window["case_id"] = sample["case_id"]
            sample["identity"] = identity
            sample["owner_proof"] = identity["owner_proof"]
            sample["cache_statistics"]["binding"] = binding
            sample_path = Path(reference["record"]["path"])
            sample_path.write_text(json.dumps(sample))
            reference["record"] = merge.evidence(sample_path)
            window["report_sha256"] = merge.digest(sample_path)
            launcher_path = Path(reference["launcher"]["path"])
            launcher = merge.read(launcher_path)
            launcher["case_id"] = sample["case_id"]
            launcher["cache_sha256"] = None
            if sample["variant"] == "cached":
                launcher["cache_native_index_identity"] = build[
                    "cache_native_index_identity"
                ]
            launcher_path.write_text(json.dumps(launcher))
            reference["launcher"] = merge.evidence(launcher_path)
            window["launcher_sha256"] = merge.digest(launcher_path)
        path.write_text(json.dumps(suite))
        index = self.index()
        doc = merge.read(index)
        doc["cache_requalification"] = merge.evidence(acceptance)
        index.write_text(json.dumps(doc))
        return path, suite, producer, index, acceptance

    def test_native_suite_requires_real_repeated_serialized_samples(self):
        path, _, producer, index, _ = self.native_suite_fixture()
        result = merge.suite_record(
            path, "tuning", producer, merge.ReportEvidence(index)
        )
        self.assertEqual(result["rows"][0]["cache_exposure"], "tuning")
        self.assertEqual(result["rows"][0]["intermediate_hits"], 3)
        index_doc = merge.read(index)
        del index_doc["records"][0]["window"]
        index.write_text(json.dumps(index_doc))
        with self.assertRaises((KeyError, ValueError)):
            merge.suite_record(path, "tuning", producer, merge.ReportEvidence(index))

    def test_native_prefix_alone_cannot_supply_consumer_qualification(self):
        path, _, producer, index, _ = self.native_suite_fixture()
        doc = merge.read(index)
        del doc["cache_requalification"]
        index.write_text(json.dumps(doc))
        with self.assertRaisesRegex(ValueError, "hashed evidence reference"):
            merge.suite_record(path, "tuning", producer, merge.ReportEvidence(index))

    def test_native_consumer_binding_and_cohort_are_not_relabelled(self):
        path, _, producer, index, acceptance = self.native_suite_fixture()
        accepted = merge.read(acceptance)
        accepted["semantic_binding"]["physical_arity"] = 14
        acceptance.write_text(json.dumps(accepted))
        doc = merge.read(index)
        doc["cache_requalification"] = merge.evidence(acceptance)
        index.write_text(json.dumps(doc))
        with self.assertRaisesRegex(ValueError, "another producer/snapshot/binding"):
            merge.suite_record(path, "tuning", producer, merge.ReportEvidence(index))

    def test_native_single_sample_or_reused_process_cannot_be_performance(self):
        path, suite, producer, index, _ = self.native_suite_fixture()
        suite["identity"]["repetitions"] = 1
        path.write_text(json.dumps(suite))
        with self.assertRaisesRegex(ValueError, "three repetitions"):
            merge.suite_record(path, "tuning", producer, merge.ReportEvidence(index))

    def rebind_acceptance(self, index, acceptance, mutate):
        doc = merge.read(acceptance)
        mutate(doc)
        acceptance.write_text(json.dumps(doc))
        evidence = merge.read(index)
        evidence["cache_requalification"] = merge.evidence(acceptance)
        index.write_text(json.dumps(evidence))

    def test_native_changed_consumer_cohort_fails(self):
        path, _, producer, index, acceptance = self.native_suite_fixture()
        doc = merge.read(acceptance)
        cohort_path = Path(doc["consumer"]["compiled_cohort"]["path"])
        cohort = merge.read(cohort_path)
        cohort["main_revision"] = "a" * 40
        cohort_path.write_text(json.dumps(cohort))
        self.rebind_acceptance(
            index,
            acceptance,
            lambda d: d["consumer"].update(compiled_cohort=merge.evidence(cohort_path)),
        )
        with self.assertRaisesRegex(ValueError, "proof/cohort/artifact"):
            merge.suite_record(path, "tuning", producer, merge.ReportEvidence(index))

    def test_native_guard_payload_and_fresh_process_are_enforced(self):
        path, suite, producer, index, _ = self.native_suite_fixture()
        guards = path.parent / "case-000-00-cached" / "guards.json"
        original = guards.read_bytes()
        guards.write_text('["DIFFERENT GUARD"]')
        with self.assertRaisesRegex(ValueError, "result/guard sidecar"):
            merge.suite_record(path, "tuning", producer, merge.ReportEvidence(index))
        guards.write_bytes(original)
        suite["measurements"][0]["fresh_process"] = False
        path.write_text(json.dumps(suite))
        with self.assertRaisesRegex(ValueError, "unsuccessful or mismatched sample"):
            merge.suite_record(path, "tuning", producer, merge.ReportEvidence(index))

    def test_native_seen_inputs_are_labelled_not_promoted_to_held_out(self):
        path, suite, producer, index, acceptance = self.native_suite_fixture()
        case_id = suite["cases"][0]["id"]
        self.rebind_acceptance(
            index,
            acceptance,
            lambda d: d["case_exposure"].update({case_id: "seen-starting-input"}),
        )
        result = merge.suite_record(
            path, "tuning", producer, merge.ReportEvidence(index)
        )
        self.assertEqual(result["rows"][0]["cache_exposure"], "seen-starting-input")
        accepted = merge.read(acceptance)
        workload_path = Path(accepted["consumer"]["final_workload"]["path"])
        workload = merge.read(workload_path)
        workload["held_out"] = workload.pop("tuning")
        workload["tuning"] = []
        workload_path.write_text(json.dumps(workload))
        identity = suite["identity"]
        identity["workload_sha256"] = merge.digest(workload_path)
        self.rebind_acceptance(
            index,
            acceptance,
            lambda d: (
                d["consumer"].update(final_workload=merge.evidence(workload_path)),
                d.update(immutable_selection=merge.selected_inputs(workload)),
            ),
        )
        with self.assertRaisesRegex(ValueError, "seen starting input"):
            merge.ReportEvidence(index).cache_qualification(identity, producer)

    def test_native_explicit_training_key_cannot_be_held_out_under_another_name(self):
        _, suite, producer, index, acceptance = self.native_suite_fixture()
        accepted = merge.read(acceptance)
        workload_path = Path(accepted["consumer"]["final_workload"]["path"])
        workload = merge.read(workload_path)
        workload["training"] = [
            {"id": "TRAINING", "powers": [1] * 15, "public_input": None}
        ]
        workload["held_out"] = [
            {"id": "HIDDEN_ALIAS", "powers": [1] * 15, "public_input": None}
        ]
        workload_path.write_text(json.dumps(workload))
        suite["identity"]["workload_sha256"] = merge.digest(workload_path)
        self.rebind_acceptance(
            index,
            acceptance,
            lambda d: (
                d["consumer"].update(final_workload=merge.evidence(workload_path)),
                d.update(immutable_selection=merge.selected_inputs(workload)),
                d["case_exposure"].update(
                    TRAINING="seen-starting-input", HIDDEN_ALIAS="held-out"
                ),
            ),
        )
        with self.assertRaisesRegex(ValueError, "starting training key"):
            merge.ReportEvidence(index).cache_qualification(suite["identity"], producer)

    def test_native_launcher_must_name_the_actual_snapshot(self):
        path, _, producer, index, _ = self.native_suite_fixture()
        doc = merge.read(index)
        for record in doc["records"]:
            launcher_path = Path(record["launcher"]["path"])
            launcher = merge.read(launcher_path)
            if launcher["variant"] == "cached":
                launcher["cache_native_index_identity"]["sha256"] = "f" * 64
                launcher_path.write_text(json.dumps(launcher))
                record["launcher"] = merge.evidence(launcher_path)
                break
        index.write_text(json.dumps(doc))
        with self.assertRaisesRegex(ValueError, "launcher native cache"):
            merge.suite_record(path, "tuning", producer, merge.ReportEvidence(index))

    def preserved_prefix_fixture(self):
        path, cache = self.cumulative_fixture()
        report = merge.read(path)
        source = report["sources"][0]
        index = cache / "cache.rrindex"
        metadata = self.write(
            "actual-shape-union-metadata.json",
            {
                "immutable_indices": [
                    {"role": "old", "path": str(index), "bytes": index.stat().st_size}
                ],
                "metadata_only": True,
                "run_concat_verified": True,
                "generation": 1,
                "compressed_bytes": report["inventory"]["compressed_bytes"],
                "preserved_record_references": 7,
            },
        )
        terminal = self.write(
            "terminal-outcome.json",
            {
                "schema": "vakint.real-opaque-union-prefix-diagnostic.v1",
                "state": "PASS",
                "independent_verification": {
                    "schema": "vakint.independent-opaque-union-prefix-preservation.v1",
                    "state": "PASS",
                    "every_captured_source_and_output_chunk_sha256_and_length_equal": True,
                    "all_chunk_files_stable_during_reads": True,
                    "native_exact_run_concatenation_verified": True,
                    "completed_q0_or_final_cache_claimed": False,
                    "semantic_qualification_or_cache_speedup_claimed": False,
                },
            },
        )
        binding = self.write(
            "actual-shape-binding.json", {"original_build_receipt": source["producer"]}
        )
        proposal = self.write(
            "actual-shape-proposal.json", {"binding": merge.evidence(binding)}
        )
        proof = self.write(
            "actual-shape-prefix-proof.json",
            {
                "schema": "vakint.postflight-real-union-prefix-proof-qualified.v1",
                "data_preservation_state": "PASS",
                "original_runtime_receipts_unchanged": True,
                "proposal": merge.evidence(proposal),
                "output": str(cache),
                "generation": 1,
                **report["inventory"],
                "artifacts": [
                    merge.evidence(metadata),
                    merge.evidence(terminal),
                    merge.evidence(index),
                ],
            },
        )
        source.update(
            kind="native-index-prefix",
            prefix_role="old",
            capture=merge.evidence(metadata),
            snapshot=merge.evidence(index),
            barrier=merge.evidence(proof),
        )
        source.pop("index")
        report["validation"] = merge.evidence(proof)
        path.write_text(json.dumps(report))
        return path, cache, proof, terminal

    def test_actual_union_prefix_shape_does_not_become_physical_completion(self):
        path, cache, _, _ = self.preserved_prefix_fixture()
        result = merge.build_record(path, cache)
        self.assertEqual(result["status"], "preserved-committed-prefix")
        self.assertFalse(result["physical_examples_complete"])
        self.assertFalse(result["performance_eligible"])
        self.assertFalse(
            result["sources"][0]["receipts"]["producer"]["raw_receipt"]["complete"]
        )

    def test_actual_union_prefix_cannot_claim_semantics_or_final_q0(self):
        path, cache, proof_path, terminal = self.preserved_prefix_fixture()
        bad = merge.read(terminal)
        bad["independent_verification"]["completed_q0_or_final_cache_claimed"] = True
        terminal.write_text(json.dumps(bad))
        proof = merge.read(proof_path)
        proof["artifacts"] = [
            merge.evidence(terminal) if Path(ref["path"]) == terminal else ref
            for ref in proof["artifacts"]
        ]
        proof_path.write_text(json.dumps(proof))
        report = merge.read(path)
        report["validation"] = merge.evidence(proof_path)
        report["sources"][0]["barrier"] = merge.evidence(proof_path)
        path.write_text(json.dumps(report))
        with self.assertRaisesRegex(ValueError, "limited scope"):
            merge.build_record(path, cache)

    def completed_source_fixture(self, job):
        from test_merge_physical_capture_results import PhysicalCaptureContracts

        physical = PhysicalCaptureContracts()
        physical.setUp()
        self.addCleanup(physical.doCleanups)
        config = physical.producer if job == "q0" else physical.quadratic_fixture()
        config_ref = (
            physical.producer_config
            if job == "q0"
            else physical.file("quadratic.final.json", config)
        )
        done = (
            physical.producer_done
            if job == "q0"
            else physical.run_fixture(config, config_ref, 13, 30.0, True)
        )
        root = Path(config["output_root"])
        cache_path = root / "run" / "running-cache.json"
        captured = merge.read(cache_path)
        captured["report"].update(
            schema="rustred.running-application-cache.v1", generation=3
        )
        cache_ref = physical.file(cache_path.relative_to(physical.root), captured)
        completion_path = Path(done["path"])
        completion = merge.read(completion_path)
        completion["terminal_artifacts"]["running-cache.json"] = {
            **cache_ref,
            "exists": True,
        }
        completion["outputs"]["running_cache"]["report"] = cache_ref
        completion_ref = physical.file(
            completion_path.relative_to(physical.root), completion
        )
        path, cache = self.cumulative_fixture()
        report = merge.read(path)
        snapshot_path = root / "run" / "snapshot.json"
        snapshot = merge.read(snapshot_path)
        native = captured["report"]
        source = report["sources"][0]
        source.update(
            scope="completed-example",
            producer=completion_ref,
            capture=cache_ref,
            snapshot=merge.evidence(snapshot_path),
            index=merge.evidence(snapshot["delta"]["manifest_path"]),
            config=config_ref,
        )
        barrier = self.write(
            "completed-source-barrier.json",
            {
                "schema": "vakint.native-cache-snapshot-barrier.v1",
                "accepted": True,
                "committed_snapshot": True,
                "producer_completed": True,
                "producer": source["producer"],
                "capture": source["capture"],
                "snapshot": source["snapshot"],
                "writer_fault": None,
                "generation": 3,
                "manifest_digest": "a" * 64,
                **{
                    key: native[key]
                    for key in (
                        "accepted_records",
                        "committed_records",
                        "pending_records",
                        "pending_bytes",
                    )
                },
            },
        )
        source["barrier"] = merge.evidence(barrier)
        self.refresh_validation(path, report)
        return (
            path,
            cache,
            physical,
            completion,
            completion_path,
            report,
            source,
            barrier,
        )

    def test_completed_source_reuses_full_existing_physical_validation(self):
        for job, case in (
            ("q0", "zero-momentum-diagnostic"),
            ("quadratic", "quadratic-four-rung-uv"),
        ):
            with (
                self.subTest(job=job),
                tempfile.TemporaryDirectory(dir=self.root) as temporary,
            ):
                original_root, self.root = self.root, Path(temporary)
                try:
                    (
                        path,
                        cache,
                        physical,
                        completion,
                        completion_path,
                        report,
                        source,
                        barrier,
                    ) = self.completed_source_fixture(job)
                    completed = merge.build_record(path, cache)
                    self.assertTrue(completed["all_source_runs_completed"])
                    self.assertEqual(completed["completed_physical_examples"], [case])
                    self.assertFalse(completed["physical_examples_complete"])
                    completion["functional_acceptance"] = False
                    physical.file(
                        completion_path.relative_to(physical.root), completion
                    )
                    source["producer"] = merge.evidence(completion_path)
                    barrier_doc = merge.read(barrier)
                    barrier_doc["producer"] = source["producer"]
                    barrier.write_text(json.dumps(barrier_doc))
                    source["barrier"] = merge.evidence(barrier)
                    self.refresh_validation(path, report)
                    with self.assertRaisesRegex(ValueError, "functional completion"):
                        merge.build_record(path, cache)
                finally:
                    self.root = original_root

    def test_required_completion_comes_only_from_actual_matching_physical_receipts(
        self,
    ):
        scope = {
            "required": ["zero-momentum-diagnostic"],
            "deferred": ["quadratic-four-rung-uv"],
        }
        for job, complete in (("q0", True), ("quadratic", False)):
            with (
                self.subTest(job=job),
                tempfile.TemporaryDirectory(dir=self.root) as temporary,
            ):
                original_root, self.root = self.root, Path(temporary)
                try:
                    path, cache, *_ = self.completed_source_fixture(job)
                    scoped = merge.build_record(path, cache, physical_scope=scope)
                    self.assertIs(scoped["physical_examples_complete"], complete)
                    self.assertEqual(scoped["physical_scope"], scope)
                    self.assertEqual(
                        scoped["completed_physical_examples"],
                        [
                            "zero-momentum-diagnostic"
                            if job == "q0"
                            else "quadratic-four-rung-uv"
                        ],
                    )
                    self.assertFalse(
                        merge.build_record(path, cache)["physical_examples_complete"]
                    )
                finally:
                    self.root = original_root

    def test_native_omitted_source_keys_and_owner_aliases_are_rejected(self):
        path, _, producer, index, acceptance = self.native_suite_fixture()
        accepted = merge.read(acceptance)
        routes_path = Path(accepted["consumer"]["training_route_qualification"]["path"])
        routes = merge.read(routes_path)
        routes["rows"] = []
        routes_path.write_text(json.dumps(routes))
        self.rebind_acceptance(
            index,
            acceptance,
            lambda d: d["consumer"].update(
                training_route_qualification=merge.evidence(routes_path)
            ),
        )
        with self.assertRaisesRegex(ValueError, "omits source starting keys"):
            merge.suite_record(path, "tuning", producer, merge.ReportEvidence(index))

    def test_byte_preserved_prefix_cannot_supply_primary_performance(self):
        report_path, cache, _, _ = self.preserved_prefix_fixture()
        native = merge.native_cache_identity(cache)
        self.paired()
        for ref, window in zip(self.records, self.window_rows):
            report_path_sample = Path(ref["record"]["path"])
            launcher_path = Path(ref["launcher"]["path"])
            report, launcher = merge.read(report_path_sample), merge.read(launcher_path)
            if launcher["variant"] == "cached":
                report.update(
                    persistent_cache=str(cache), cache_native_index_identity=native
                )
                launcher.update(cache_sha256=None, cache_native_index_identity=native)
            report_path_sample.write_text(json.dumps(report))
            launcher_path.write_text(json.dumps(launcher))
            ref.update(
                record=merge.evidence(report_path_sample),
                launcher=merge.evidence(launcher_path),
            )
            window.update(
                report_sha256=merge.digest(report_path_sample),
                launcher_sha256=merge.digest(launcher_path),
            )
        output = self.root / "must-not-exist.json"
        base = self.write(
            "base-empty.json", {"reduction_runs": {"historical": {"unchanged": True}}}
        )
        with (
            patch(
                "sys.argv",
                [
                    "adapter",
                    "--base",
                    str(base),
                    "--build",
                    str(report_path),
                    "--cache",
                    str(cache),
                    "--evidence",
                    str(self.index()),
                    "--output",
                    str(output),
                ],
            ),
            self.assertRaisesRegex(ValueError, "not prefix preservation"),
        ):
            merge.main()
        self.assertFalse(output.exists())

    def test_primary_and_suite_keep_their_distinct_elf_identities(self):
        _, suite, producer, index, acceptance = self.native_suite_fixture()
        primary = self.root / "SYNTHETIC-PRIMARY-ELF"
        primary.write_text("SYNTHETIC ONLY; DISTINCT FROM MEASUREMENT ELF")
        artifact = merge.read(self.artifact)
        artifact["artifacts"] = [
            artifact.pop("binary"),
            {"sha256": merge.digest(primary), "immutable_path": str(primary)},
        ]
        self.artifact.write_text(json.dumps(artifact))
        self.rebind_acceptance(
            index,
            acceptance,
            lambda d: d["consumer"].update(
                collector_artifacts=merge.evidence(self.artifact)
            ),
        )
        index_doc = merge.read(index)
        for ref in index_doc["records"]:
            ref["artifact"] = merge.evidence(self.artifact)
        index.write_text(json.dumps(index_doc))
        provenance = merge.ReportEvidence(index)
        identity = copy.deepcopy(suite["identity"])
        identity["executable_sha256"] = merge.digest(primary)
        self.assertIsNotNone(
            provenance.cache_qualification(identity, producer, primary=True)
        )
        with self.assertRaisesRegex(ValueError, "proof/cohort/artifact"):
            provenance.cache_qualification(identity, producer)

    def test_union_roles_cannot_hide_or_relabel_actual_sources(self):
        for role in ("output", "q0", "omit-q0"):
            with (
                self.subTest(role=role),
                tempfile.TemporaryDirectory(dir=self.root) as temporary,
            ):
                original_root, self.root = self.root, Path(temporary)
                try:
                    path, cache, proof_path, _ = self.preserved_prefix_fixture()
                    report = merge.read(path)
                    source = report["sources"][0]
                    metadata_path = Path(source["capture"]["path"])
                    metadata = merge.read(metadata_path)
                    if role == "omit-q0":
                        metadata["immutable_indices"].append(
                            {**metadata["immutable_indices"][0], "role": "q0"}
                        )
                    else:
                        metadata["immutable_indices"][0]["role"] = role
                        source["prefix_role"] = role
                    metadata_path.write_text(json.dumps(metadata))
                    proof = merge.read(proof_path)
                    proof["artifacts"] = [
                        merge.evidence(metadata_path)
                        if Path(ref["path"]) == metadata_path
                        else ref
                        for ref in proof["artifacts"]
                    ]
                    proof_path.write_text(json.dumps(proof))
                    source.update(
                        capture=merge.evidence(metadata_path),
                        barrier=merge.evidence(proof_path),
                    )
                    report["validation"] = merge.evidence(proof_path)
                    path.write_text(json.dumps(report))
                    with self.assertRaisesRegex(
                        ValueError,
                        "not a source prefix|old producer|every captured source",
                    ):
                        merge.build_record(path, cache)
                finally:
                    self.root = original_root

    def test_native_suite_cli_appends_without_rewriting_history(self):
        suite_path, _, producer, index, _ = self.native_suite_fixture()
        base = self.write(
            "native-suite-base.json",
            {
                "selected_reduction_run": "HISTORICAL",
                "reduction_runs": {"HISTORICAL": {"raw": [1, 2, 3]}},
            },
        )
        historical = base.read_bytes()
        output = self.root / "native-suite-output.json"
        with patch(
            "sys.argv",
            [
                "adapter",
                "--base",
                str(base),
                "--build",
                producer["evidence"]["path"],
                "--cache",
                producer["cache_directory"],
                "--tuning",
                str(suite_path),
                "--evidence",
                str(index),
                "--output",
                str(output),
            ],
        ):
            merge.main()
        result = merge.read(output)
        self.assertEqual(base.read_bytes(), historical)
        self.assertEqual(result["reduction_runs"], merge.read(base)["reduction_runs"])
        report = next(iter(result["five_loop_runs"].values()))
        self.assertEqual(report["tuning"]["rows"][0]["status"], "measured")
        self.assertFalse(report["cache"]["physical_examples_complete"])

    def alter_public_import(self, index, acceptance, relative, mutate):
        """Rebind every outer hash so tests exercise parsed facts, not stale SHA."""
        accepted = merge.read(acceptance)
        observation_path = Path(accepted["cache_import"]["path"])
        root = observation_path.parent
        target = root / relative
        value = merge.read(target)
        mutate(value)
        target.write_text(json.dumps(value))
        packet_path, go_path = root / "packet.json", root / "root-go.json"
        packet, go = merge.read(packet_path), merge.read(go_path)
        if relative == "descriptor.json":
            packet["descriptor"] = merge.evidence(target)
            packet_path.write_text(json.dumps(packet))
        go["packet_sha256"] = merge.digest(packet_path)
        go_path.write_text(json.dumps(go))
        admission_path, invocation_path = (
            root / "run1/admission.json",
            root / "run1/invocation-start.json",
        )
        admission, invocation = merge.read(admission_path), merge.read(invocation_path)
        admission["packet_sha256"] = invocation["packet_sha256"] = merge.digest(
            packet_path
        )
        admission_path.write_text(json.dumps(admission))
        invocation_path.write_text(json.dumps(invocation))
        worker_path, terminal_path = (
            root / "run1/worker-result.json",
            root / "run1/terminal-outcome.json",
        )
        worker, terminal = merge.read(worker_path), merge.read(terminal_path)
        for name, document in (
            ("public-import-report.json", worker["import_report"]),
            ("data-load-report.json", worker["data_report"]),
            ("cache-state-before-repeat.json", worker["final_attachment_state"]),
        ):
            (root / "run1" / name).write_text(json.dumps(document))
        terminal.update(
            packet_sha256=merge.digest(packet_path),
            root_go_sha256=merge.digest(go_path),
            worker_result=merge.evidence(worker_path),
        )
        terminal_path.write_text(json.dumps(terminal))
        launch_path = root / "actual-launch-receipt.json"
        launch = merge.read(launch_path)
        launch["admission"] = admission
        launch["actual_host"] = merge.read(root / "run1/actual-host.json")
        launch["artifact_sha256"] = {
            name: merge.digest(name) for name in launch["artifact_sha256"]
        }
        launch_path.write_text(json.dumps(launch))
        observation = merge.read(observation_path)
        observation["controller_terminal"] = terminal
        observation["artifact_sha256"] = {
            name: merge.digest(name) for name in observation["artifact_sha256"]
        }
        observation_path.write_text(json.dumps(observation))
        self.rebind_acceptance(
            index,
            acceptance,
            lambda d: d.update(cache_import=merge.evidence(observation_path)),
        )

    def test_import_check_boolean_alone_cannot_qualify(self):
        _, suite, producer, index, acceptance = self.native_suite_fixture()
        self.rebind_acceptance(index, acceptance, lambda d: d.pop("cache_import"))
        with self.assertRaisesRegex(ValueError, "hashed evidence reference"):
            merge.ReportEvidence(index).cache_qualification(suite["identity"], producer)

    def test_public_import_derives_chunk_counts_from_descriptor(self):
        _, suite, producer, index, acceptance = self.native_suite_fixture()
        imported = merge.ReportEvidence(index).cache_qualification(
            suite["identity"], producer
        )["public_import"]
        self.assertEqual(imported["raw_import_report"]["copied_chunks"], 3)
        self.assertEqual(imported["raw_import_report"]["entries"], 4)
        receipt = merge.read(acceptance)["cache_import"]
        observed = merge.read(receipt["path"])
        worker = observed["controller_terminal"]["worker_result"]
        self.assertNotIn(
            "all_both_indexes_and_18639_chunks_length_sha256_equal",
            merge.read(worker["path"]),
        )

    def test_import_requires_actual_successful_worker_and_controller(self):
        for relative, mutation in (
            ("run1/worker-result.json", lambda d: d.update(state="FAIL")),
            ("run1/terminal-outcome.json", lambda d: d.update(worker_returncode=1)),
            (
                "run1/terminal-outcome.json",
                lambda d: d.update(whole_elapsed_seconds=3600),
            ),
        ):
            with (
                self.subTest(relative=relative),
                tempfile.TemporaryDirectory(dir=self.root) as temporary,
            ):
                original_root, self.root = self.root, Path(temporary)
                try:
                    _, suite, producer, index, acceptance = self.native_suite_fixture()
                    self.alter_public_import(index, acceptance, relative, mutation)
                    with self.assertRaisesRegex(
                        ValueError, "public import|actual public import"
                    ):
                        merge.ReportEvidence(index).cache_qualification(
                            suite["identity"], producer
                        )
                finally:
                    self.root = original_root

    def test_import_descriptor_and_report_must_match_exact_snapshot(self):
        changes = (
            ("descriptor.json", lambda d: d["binding"].update(physical_arity=14)),
            (
                "descriptor.json",
                lambda d: d["files"]["cache/cache.rrindex"].update(sha256="0" * 64),
            ),
            ("run1/worker-result.json", lambda d: d["import_report"].update(entries=1)),
            (
                "run1/worker-result.json",
                lambda d: d.update(exact_native_files_verified=1),
            ),
        )
        for relative, mutation in changes:
            with (
                self.subTest(relative=relative),
                tempfile.TemporaryDirectory(dir=self.root) as temporary,
            ):
                original_root, self.root = self.root, Path(temporary)
                try:
                    _, suite, producer, index, acceptance = self.native_suite_fixture()
                    self.alter_public_import(index, acceptance, relative, mutation)
                    with self.assertRaisesRegex(
                        ValueError,
                        "descriptor names another|actual public import report",
                    ):
                        merge.ReportEvidence(index).cache_qualification(
                            suite["identity"], producer
                        )
                finally:
                    self.root = original_root

    def test_import_host_and_data_provenance_cannot_be_substituted(self):
        changes = (
            (
                "../public-import-data/manifest.json",
                lambda d: d.update(version="WRONG-DATA-MANIFEST"),
            ),
            ("run1/actual-host.json", lambda d: d.update(core_sha256="0" * 64)),
            (
                "run1/worker-result.json",
                lambda d: d["data_report"].update(version="WRONG-VERSION"),
            ),
            (
                "run1/worker-result.json",
                lambda d: d["final_attachment_state"]["5"]["saved_application"][
                    "binding"
                ].update(physical_arity=14),
            ),
            (
                "run1/worker-result.json",
                lambda d: d["final_attachment_state"]["5"]["saved_application"][
                    "persistent"
                ].update(lookups=1),
            ),
            (
                "run1/worker-result.json",
                lambda d: d.update(cache_attachment_performed=True),
            ),
        )
        for relative, mutation in changes:
            with (
                self.subTest(relative=relative),
                tempfile.TemporaryDirectory(dir=self.root) as temporary,
            ):
                original_root, self.root = self.root, Path(temporary)
                try:
                    _, suite, producer, index, acceptance = self.native_suite_fixture()
                    self.alter_public_import(index, acceptance, relative, mutation)
                    with self.assertRaisesRegex(
                        ValueError,
                        "host artifact|data/binding/attachment|actual public import report|another data manifest",
                    ):
                        merge.ReportEvidence(index).cache_qualification(
                            suite["identity"], producer
                        )
                finally:
                    self.root = original_root

    def test_import_observer_cannot_hide_replacement_invocation(self):
        _, suite, producer, index, acceptance = self.native_suite_fixture()
        self.alter_public_import(
            index,
            acceptance,
            "actual-terminal-observation.json",
            lambda d: d["current_service"].update(InvocationID="OTHER-INVOCATION"),
        )
        with self.assertRaisesRegex(ValueError, "launch/journal identity"):
            merge.ReportEvidence(index).cache_qualification(suite["identity"], producer)


class MetadataRequalificationTests(unittest.TestCase):
    """Synthetic manifest/opaque-payload parity; no native import or evaluation."""

    def cache_fixture(self, count, directory, files):
        return {"path": "cache/prepared.gz"}

    def setUp(self):
        self.base = CumulativeCacheContractTests()
        self.base.setUp()
        self.addCleanup(self.base.doCleanups)
        self.root, self.write = self.base.root, self.base.write
        self.source = self.root / "public-import-data"
        root = {
            "schema": "vakint.external-data.v1",
            "version": "v0.1.2",
            "loop_packs": [],
            "provenance_assets": [],
        }
        for count in range(1, 6):
            directory = self.source / f"{count}L"
            files = {}
            for name in (
                "ibps/program.bin",
                "masters.json",
                "graphs/master.dot",
                "cache/prepared.gz",
            ):
                payload = directory / name
                payload.parent.mkdir(parents=True, exist_ok=True)
                payload.write_text(f"SYNTHETIC {count}L {name}")
                files[name] = {
                    "sha256": merge.digest(payload),
                    "bytes": payload.stat().st_size,
                }
            manifest = self.write(
                f"public-import-data/{count}L/manifest.json",
                {
                    "schema": "vakint.external-data.loop.v1",
                    "version": "v0.1.2",
                    "loop_count": count,
                    "normalization": "FMFTandMATAD",
                    "coverage": {"predicate": "SYNTHETIC"},
                    "cache": self.cache_fixture(count, directory, files),
                    "masters": [{"id": "TEST"}],
                    "ibps": {"path": "ibps/program.bin"},
                    "files": files,
                },
            )
            root["loop_packs"].append(
                {
                    "loop_count": count,
                    "manifest": f"{count}L/manifest.json",
                    "manifest_sha256": merge.digest(manifest),
                    "archive": f"vakint-external-data-v0.1.2-{count}L.tar.gz",
                    "archive_sha256": None,
                    "archive_bytes": None,
                }
            )
        self.old_manifest = self.write("public-import-data/manifest.json", root)
        _, self.suite, self.producer, self.index, self.acceptance = (
            self.base.native_suite_fixture(data_manifest=self.old_manifest)
        )
        self.accepted = merge.read(self.acceptance)
        proof = merge.read(self.accepted["consumer"]["owner_proof"]["path"])
        proof.update(
            schema="vakint.public-lowering-owner-map-proof.v1",
            application_sha256="a" * 64,
            kernel_revision=proof["compiled_cohort"]["kernel_revision"],
            mapper_executable_sha256=merge.digest(self.base.binary),
            mapper_git_revision=proof["compiled_cohort"]["main_revision"],
            mapper_profile="release",
            mapper_rustc="SYNTHETIC",
            mapper_source_sha256="c" * 64,
            native_options="SYNTHETIC",
            rustred_revision=proof["compiled_cohort"]["native_revision"],
            rows=[
                {
                    "id": "native-tuning-synthetic",
                    "logical_cohort": "tuning",
                    "public_input_source_sha256": "d" * 64,
                    "input_sha256": "f" * 64,
                    "actual_lowered_roots": [[1] * 15],
                    "protected_owner_roots": [[1] * 15],
                    "protected_root_set_sha256": "e" * 64,
                    "application_binding": proof["semantic_binding"],
                    "complete_owner_transport": True,
                    "reduction_succeeded": True,
                    "mapping_error": None,
                    "native_transports": [
                        {
                            "powers": [1] * 15,
                            "nonzero_conditions": ["TEST_GUARD"],
                            "owner_endpoints": [
                                {"powers": [1] * 15, "coefficient": "1"}
                            ],
                        }
                    ],
                }
            ],
        )
        proof["loop_manifest_sha256"] = root["loop_packs"][4]["manifest_sha256"]
        proof["compiled_cohort"].update(
            data_manifest_sha256=merge.digest(self.old_manifest),
            loop_manifest_sha256=proof["loop_manifest_sha256"],
        )
        workload = merge.read(self.accepted["consumer"]["final_workload"]["path"])
        for index in range(12):
            key = [index + 4] * 15
            case_id = f"frozen-held-out-{index}"
            workload["held_out"].append(
                {
                    "id": case_id,
                    "powers": key,
                    "public_input": f"SYNTHETIC_HELD_OUT_{index}",
                    "lowered_roots": [key],
                    "protected_roots": [key],
                }
            )
            row = copy.deepcopy(proof["rows"][0])
            row.update(
                id=case_id,
                logical_cohort="held-out",
                actual_lowered_roots=[key],
                protected_owner_roots=[key],
                native_transports=[
                    {
                        "powers": key,
                        "nonzero_conditions": ["TEST_GUARD"],
                        "owner_endpoints": [{"powers": key, "coefficient": "1"}],
                    }
                ],
            )
            proof["rows"].append(row)
            self.accepted["case_exposure"][case_id] = "held-out"
        self.write_workers(proof, "frozen")
        self.frozen_proof = self.write("frozen-owner-proof.json", proof)
        workload["owner_proof"] = merge.evidence(self.frozen_proof)
        self.accepted["immutable_selection"] = merge.selected_inputs(workload)
        workload["tuning"][0].update(
            lowered_roots=[[1] * 15], protected_roots=[[1] * 15]
        )
        self.frozen_workload = self.write("frozen-workload.json", workload)
        self.target = self.root / "published-data"
        shutil.copytree(self.source, self.target)
        packet = merge.read(
            Path(self.accepted["cache_import"]["path"]).parent / "packet.json"
        )
        descriptor = Path(packet["descriptor"]["path"])
        selected = {
            "loop_count": 5,
            "manifest": "5L/cache/running.json",
            "manifest_sha256": merge.digest(descriptor),
            "manifest_bytes": descriptor.stat().st_size,
        }
        (self.target / selected["manifest"]).write_bytes(descriptor.read_bytes())
        root.update(version="v0.1.3", runtime_cache_snapshots=[selected])
        for pack in root["loop_packs"]:
            path = self.target / pack["manifest"]
            loop = merge.read(path)
            loop["version"] = "v0.1.3"
            if pack["loop_count"] == 5:
                loop["files"]["cache/running.json"] = {
                    "sha256": selected["manifest_sha256"],
                    "bytes": selected["manifest_bytes"],
                }
            path.write_text(json.dumps(loop))
            pack.update(
                manifest_sha256=merge.digest(path),
                archive=f"vakint-external-data-v0.1.3-{pack['loop_count']}L.tar.gz",
                archive_sha256="a" * 64,
                archive_bytes=42,
            )
        self.new_manifest = self.write("published-data/manifest.json", root)
        self.relation_path = self.write(
            "metadata-relation.json",
            {
                "schema": "vakint.external-data-metadata-requalification.v1",
                "source_manifest": merge.evidence(self.old_manifest),
                "consumer_manifest": merge.evidence(self.new_manifest),
                "runtime_cache_snapshot": selected,
                "frozen_workload": merge.evidence(self.frozen_workload),
                "frozen_owner_proof": merge.evidence(self.frozen_proof),
                "frozen_collector_artifacts": self.accepted["consumer"][
                    "collector_artifacts"
                ],
            },
        )
        self.consumer_proof = self.write("consumer-proof.json", proof)
        self.consumer_workload = self.write("consumer-workload.json", workload)
        self.refresh()

    def write_workers(self, proof, role):
        for index, row in enumerate(proof["rows"]):
            worker = self.write(
                f"{role}-workers/{index}.json",
                {
                    "schema": "vakint.five-loop-cache-measurement.v1",
                    "success": True,
                    "worker_success": True,
                    "fresh_process": True,
                    "variant": "uncached",
                    "sample": 0,
                    "supervision": {
                        "success": True,
                        "exit_code": 0,
                        "deadline_exhausted": False,
                    },
                    "case_id": row["id"],
                    "input_sha256": row["input_sha256"],
                    "actual_lowered_roots": row["actual_lowered_roots"],
                    "identity": {
                        "data_manifest_sha256": proof["data_manifest_sha256"],
                        "build_git_revision": proof["compiled_cohort"]["main_revision"],
                        "executable_sha256": proof["measurement_executable_sha256"],
                    },
                },
            )
            row.update(receipt_path=str(worker), receipt_sha256=merge.digest(worker))

    def refresh(self):
        root = merge.read(self.new_manifest)
        for pack in root["loop_packs"]:
            pack["manifest_sha256"] = merge.digest(self.target / pack["manifest"])
        self.new_manifest.write_text(json.dumps(root))
        relation = merge.read(self.relation_path)
        relation["consumer_manifest"] = merge.evidence(self.new_manifest)
        self.relation_path.write_text(json.dumps(relation))
        proof = merge.read(self.consumer_proof)
        proof["data_manifest_sha256"] = merge.digest(self.new_manifest)
        proof["loop_manifest_sha256"] = root["loop_packs"][4]["manifest_sha256"]
        proof["compiled_cohort"].update(
            data_manifest_sha256=proof["data_manifest_sha256"],
            loop_manifest_sha256=proof["loop_manifest_sha256"],
        )
        self.write_workers(proof, "consumer")
        self.consumer_proof.write_text(json.dumps(proof))
        cohort_path = Path(self.accepted["consumer"]["compiled_cohort"]["path"])
        cohort_path.write_text(json.dumps(proof["compiled_cohort"]))
        self.accepted["consumer"]["compiled_cohort"] = merge.evidence(cohort_path)
        route_path = Path(
            self.accepted["consumer"]["training_route_qualification"]["path"]
        )
        route = merge.read(route_path)
        route["compiled_cohort"] = proof["compiled_cohort"]
        route_path.write_text(json.dumps(route))
        self.accepted["consumer"]["training_route_qualification"] = merge.evidence(
            route_path
        )
        workload = merge.read(self.consumer_workload)
        workload.update(
            data_manifest_sha256=merge.digest(self.new_manifest),
            owner_proof=merge.evidence(self.consumer_proof),
        )
        self.consumer_workload.write_text(json.dumps(workload))
        self.accepted["consumer"].update(
            owner_proof=merge.evidence(self.consumer_proof),
            final_workload=merge.evidence(self.consumer_workload),
        )
        union_path = Path(self.accepted["consumer"]["protected_union"]["path"])
        union = merge.read(union_path)
        union["owner_proof"] = merge.evidence(self.consumer_proof)
        union["compiled_cohort"] = merge.evidence(cohort_path)
        protected = [row["powers"] for row in workload["held_out"]]
        union.update(
            raw_lowered_keys=protected,
            direct_owner_keys=protected,
            declared_case_powers=protected,
            forbidden_direct_target_or_owner_keys=protected,
        )
        union_path.write_text(json.dumps(union))
        self.accepted["consumer"]["protected_union"] = merge.evidence(union_path)
        alias_path = Path(self.accepted["consumer"]["alias_exclusion"]["path"])
        alias = merge.read(alias_path)
        alias["protected_union"] = merge.evidence(union_path)
        alias["qualification"] = merge.evidence(route_path)
        alias_path.write_text(json.dumps(alias))
        self.accepted["consumer"]["alias_exclusion"] = merge.evidence(alias_path)
        self.accepted["data_metadata_requalification"] = merge.evidence(
            self.relation_path
        )
        self.suite["identity"].update(
            data_manifest_sha256=merge.digest(self.new_manifest),
            workload_sha256=merge.digest(self.consumer_workload),
            owner_proof={
                **merge.evidence(self.consumer_proof),
                "compiled_cohort": proof["compiled_cohort"],
            },
        )
        self.base.rebind_acceptance(
            self.index, self.acceptance, lambda doc: doc.update(self.accepted)
        )

    def qualify(self):
        return merge.ReportEvidence(self.index).cache_qualification(
            self.suite["identity"], self.producer
        )

    def test_metadata_only_release_keeps_source_import_identity(self):
        original = copy.deepcopy(self.producer["raw_report"])
        qualified = self.qualify()
        relation = qualified["data_metadata_requalification"]
        self.assertEqual(relation["verified_original_files"], 20)
        self.assertEqual(len(merge.read(self.frozen_workload)["held_out"]), 12)
        self.assertEqual(len(merge.read(self.frozen_proof)["rows"]), 13)
        self.assertNotEqual(
            relation["source_manifest"]["sha256"],
            relation["consumer_manifest"]["sha256"],
        )
        self.assertEqual(self.producer["raw_report"], original)
        self.assertEqual(
            merge.read(qualified["public_import"]["packet"]["path"])["frozen_data"],
            str(self.source),
        )

    def test_changed_dataset_needs_explicit_relation(self):
        self.base.rebind_acceptance(
            self.index,
            self.acceptance,
            lambda doc: doc.pop("data_metadata_requalification"),
        )
        with self.assertRaisesRegex(ValueError, "explicit metadata requalification"):
            self.qualify()

    def change_loop(self, change, count=5):
        path = self.target / f"{count}L/manifest.json"
        loop = merge.read(path)
        change(loop)
        path.write_text(json.dumps(loop))
        self.refresh()

    def test_scientific_manifest_changes_are_rejected(self):
        for key in ("normalization", "coverage", "cache", "masters", "ibps"):
            with self.subTest(key=key):
                case = MetadataRequalificationTests()
                case.setUp()
                self.addCleanup(case.doCleanups)
                case.change_loop(
                    lambda doc, key=key: doc.update({key: "ILLICIT CHANGE"})
                )
                with self.assertRaisesRegex(ValueError, "scientific semantics"):
                    case.qualify()

    def test_original_payload_bytes_cannot_change(self):
        for name in (
            "ibps/program.bin",
            "masters.json",
            "graphs/master.dot",
            "cache/prepared.gz",
        ):
            with self.subTest(name=name):
                case = MetadataRequalificationTests()
                case.setUp()
                self.addCleanup(case.doCleanups)
                (case.target / "5L" / name).write_text("ILLICIT PAYLOAD CHANGE")
                with self.assertRaisesRegex(ValueError, "payload bytes"):
                    case.qualify()

    def test_rehashed_payload_replacement_is_not_metadata(self):
        payload = self.target / "5L/ibps/program.bin"
        payload.write_text("ILLICIT BUT REHASHED PROGRAM")
        self.change_loop(
            lambda doc: doc["files"].update(
                {
                    "ibps/program.bin": {
                        "sha256": merge.digest(payload),
                        "bytes": payload.stat().st_size,
                    }
                }
            )
        )
        with self.assertRaisesRegex(ValueError, "scientific file inventory"):
            self.qualify()

    def test_hidden_payloads_in_either_or_both_trees_are_rejected(self):
        for both, declared in ((False, False), (True, False), (False, True)):
            with self.subTest(both=both, declared=declared):
                case = MetadataRequalificationTests()
                case.setUp()
                self.addCleanup(case.doCleanups)
                hidden = case.target / "5L/hidden.bin"
                hidden.write_text("UNDECLARED SYNTHETIC PAYLOAD")
                if both:
                    (case.source / "5L/hidden.bin").write_bytes(hidden.read_bytes())
                if declared:
                    case.change_loop(
                        lambda doc, hidden=hidden: doc["files"].update(
                            {
                                "hidden.bin": {
                                    "sha256": merge.digest(hidden),
                                    "bytes": hidden.stat().st_size,
                                }
                            }
                        )
                    )
                with self.assertRaisesRegex(
                    ValueError, "hidden payload|scientific file inventory"
                ):
                    case.qualify()

    def test_added_descriptor_requires_original_binding_and_index(self):
        for field in ("binding", "index"):
            with self.subTest(field=field):
                case = MetadataRequalificationTests()
                case.setUp()
                self.addCleanup(case.doCleanups)
                descriptor_path = case.target / "5L/cache/running.json"
                descriptor = merge.read(descriptor_path)
                if field == "binding":
                    descriptor["binding"]["physical_arity"] = 14
                else:
                    descriptor["files"]["cache/cache.rrindex"]["sha256"] = "0" * 64
                descriptor_path.write_text(json.dumps(descriptor))
                root = merge.read(case.new_manifest)
                selected = root["runtime_cache_snapshots"][0]
                selected.update(
                    manifest_sha256=merge.digest(descriptor_path),
                    manifest_bytes=descriptor_path.stat().st_size,
                )
                case.new_manifest.write_text(json.dumps(root))
                relation = merge.read(case.relation_path)
                relation["runtime_cache_snapshot"] = selected
                case.relation_path.write_text(json.dumps(relation))
                case.change_loop(
                    lambda doc, selected=selected: doc["files"].update(
                        {
                            "cache/running.json": {
                                "sha256": selected["manifest_sha256"],
                                "bytes": selected["manifest_bytes"],
                            }
                        }
                    )
                )
                with self.assertRaisesRegex(ValueError, "cache descriptor differs"):
                    case.qualify()

    def test_frozen_case_and_owner_keys_cannot_be_transformed(self):
        for field in ("powers", "lowered_roots", "protected_roots"):
            with self.subTest(field=field):
                case = MetadataRequalificationTests()
                case.setUp()
                self.addCleanup(case.doCleanups)
                workload = merge.read(case.consumer_workload)
                row = workload["tuning"][0]
                if field == "powers":
                    row[field][0] += 1
                else:
                    row[field][0][0] += 1
                case.consumer_workload.write_text(json.dumps(workload))
                case.accepted["immutable_selection"] = merge.selected_inputs(workload)
                case.refresh()
                with self.assertRaisesRegex(
                    ValueError, "transformed frozen workload/owner keys"
                ):
                    case.qualify()

    def test_native_transport_guards_and_executable_provenance_stay_exact(self):
        for field in (
            "native_transports",
            "mapper_executable_sha256",
            "compiled_cohort",
            "reduction_succeeded",
        ):
            with self.subTest(field=field):
                case = MetadataRequalificationTests()
                case.setUp()
                self.addCleanup(case.doCleanups)
                frozen = merge.read(case.frozen_proof)
                if field == "native_transports":
                    frozen["rows"][0][field][0]["nonzero_conditions"] = [
                        "CHANGED_GUARD"
                    ]
                elif field == "compiled_cohort":
                    frozen[field]["main_revision"] = "a" * 40
                elif field == "reduction_succeeded":
                    frozen["rows"][0][field] = False
                else:
                    frozen[field] = "0" * 64
                case.frozen_proof.write_text(json.dumps(frozen))
                workload = merge.read(case.frozen_workload)
                workload["owner_proof"] = merge.evidence(case.frozen_proof)
                case.frozen_workload.write_text(json.dumps(workload))
                relation = merge.read(case.relation_path)
                relation.update(
                    frozen_owner_proof=merge.evidence(case.frozen_proof),
                    frozen_workload=merge.evidence(case.frozen_workload),
                )
                case.relation_path.write_text(json.dumps(relation))
                case.refresh()
                with self.assertRaisesRegex(
                    ValueError,
                    "proof/artifact/cohort|native owner transports/guards",
                ):
                    case.qualify()

    def test_other_root_changes_and_extra_snapshot_refs_are_rejected(self):
        root = merge.read(self.new_manifest)
        root["provenance_assets"] = [{"archive": "ILLICIT NEW PAYLOAD"}]
        self.new_manifest.write_text(json.dumps(root))
        self.refresh()
        with self.assertRaisesRegex(ValueError, "release semantics/provenance"):
            self.qualify()
        root["provenance_assets"] = []
        root["runtime_cache_snapshots"].append(root["runtime_cache_snapshots"][0])
        self.new_manifest.write_text(json.dumps(root))
        self.refresh()
        with self.assertRaisesRegex(ValueError, "exactly the selected"):
            self.qualify()

    def test_duplicate_keys_and_symlinked_payloads_are_rejected(self):
        path = self.target / "5L/manifest.json"
        path.write_text('{"version":"v9.9.9",' + path.read_text()[1:])
        self.refresh()
        with self.assertRaisesRegex(ValueError, "duplicate metadata JSON"):
            self.qualify()
        case = MetadataRequalificationTests()
        case.setUp()
        self.addCleanup(case.doCleanups)
        payload = case.target / "5L/masters.json"
        payload.unlink()
        payload.symlink_to(case.source / "5L/masters.json")
        with self.assertRaisesRegex(ValueError, "symlink"):
            case.qualify()

    @unittest.skipUnless(hasattr(os, "mkfifo"), "requires POSIX named pipes")
    def test_special_payload_is_rejected_without_reading_it(self):
        os.mkfifo(self.target / "5L/hidden-pipe")
        with self.assertRaisesRegex(ValueError, "special file"):
            self.qualify()

    def test_frozen_proof_cannot_be_rebound_from_another_dataset(self):
        proof = merge.read(self.frozen_proof)
        proof["data_manifest_sha256"] = "b" * 64
        self.frozen_proof.write_text(json.dumps(proof))
        workload = merge.read(self.frozen_workload)
        workload.update(
            data_manifest_sha256="b" * 64, owner_proof=merge.evidence(self.frozen_proof)
        )
        self.frozen_workload.write_text(json.dumps(workload))
        relation = merge.read(self.relation_path)
        relation.update(
            frozen_owner_proof=merge.evidence(self.frozen_proof),
            frozen_workload=merge.evidence(self.frozen_workload),
        )
        self.relation_path.write_text(json.dumps(relation))
        self.refresh()
        with self.assertRaisesRegex(ValueError, "proof names another dataset"):
            self.qualify()

    def test_metadata_relation_keeps_distinct_authenticated_code_cohorts(self):
        proof = merge.read(self.frozen_proof)
        proof["compiled_cohort"]["main_revision"] = "a" * 40
        proof["mapper_git_revision"] = "a" * 40
        binary = self.root / "OLD-SYNTHETIC-COLLECTOR"
        binary.write_text("SYNTHETIC OLD COLLECTOR; NOT EXECUTABLE")
        proof["measurement_executable_sha256"] = proof["mapper_executable_sha256"] = (
            merge.digest(binary)
        )
        self.write_workers(proof, "distinct-frozen")
        self.frozen_proof.write_text(json.dumps(proof))
        workload = merge.read(self.frozen_workload)
        workload["owner_proof"] = merge.evidence(self.frozen_proof)
        self.frozen_workload.write_text(json.dumps(workload))
        artifact = merge.read(self.base.artifact)
        artifact["stamps"]["VAKINT_CACHE_BUILD_GIT_REVISION"] = "a" * 40
        artifact["binary"] = {
            "sha256": merge.digest(binary),
            "immutable_path": str(binary),
        }
        frozen_artifact = self.write("old-collector-artifact.json", artifact)
        relation = merge.read(self.relation_path)
        relation.update(
            frozen_owner_proof=merge.evidence(self.frozen_proof),
            frozen_workload=merge.evidence(self.frozen_workload),
            frozen_collector_artifacts=merge.evidence(frozen_artifact),
        )
        self.relation_path.write_text(json.dumps(relation))
        self.refresh()
        qualified = self.qualify()["data_metadata_requalification"]
        self.assertEqual(
            qualified["frozen_collector_artifacts"]["sha256"],
            merge.digest(frozen_artifact),
        )
        self.assertNotEqual(
            merge.read(self.frozen_proof)["compiled_cohort"]["main_revision"],
            self.suite["identity"]["build_git_revision"],
        )

    def test_new_data_proof_cannot_reuse_old_worker_receipts(self):
        qualified = self.qualify()
        docs = {
            name: merge.read(self.accepted["consumer"][name]["path"])
            for name in ("final_workload", "owner_proof", "collector_artifacts")
        }
        docs["owner_proof"]["rows"] = merge.read(self.frozen_proof)["rows"]
        with self.assertRaisesRegex(
            ValueError, "worker receipt belongs to another dataset"
        ):
            merge.ReportEvidence(self.index).data_metadata_requalification(
                merge.evidence(self.relation_path),
                self.suite["identity"],
                self.producer,
                self.accepted,
                docs,
                qualified["public_import"],
                self.acceptance.parent,
            )

    def test_all_fresh_samples_keep_new_dataset_identity(self):
        suite = self.suite
        index_doc = merge.read(self.index)
        window_path = self.root / "window.json"
        window = merge.read(window_path)
        for sample, reference, observed in zip(
            suite["measurements"], index_doc["records"], window["samples"]
        ):
            sample["owner_proof"] = suite["identity"]["owner_proof"]
            sample_path = Path(reference["record"]["path"])
            sample_path.write_text(json.dumps(sample))
            reference["record"] = merge.evidence(sample_path)
            observed["report_sha256"] = merge.digest(sample_path)
        window_path.write_text(json.dumps(window))
        for reference in index_doc["records"]:
            reference["window"] = merge.evidence(window_path)
        self.index.write_text(json.dumps(index_doc))
        suite_path = self.root / "suite.json"
        suite_path.write_text(json.dumps(suite))
        measured = merge.suite_record(
            suite_path, "tuning", self.producer, merge.ReportEvidence(self.index)
        )
        self.assertEqual(
            measured["identity"]["data_manifest_sha256"],
            merge.digest(self.new_manifest),
        )
        first = suite["measurements"][0]
        first["identity"] = copy.deepcopy(first["identity"])
        first["identity"]["data_manifest_sha256"] = merge.digest(self.old_manifest)
        suite_path.write_text(json.dumps(suite))
        with self.assertRaisesRegex(ValueError, "unsuccessful or mismatched sample"):
            merge.suite_record(
                suite_path, "tuning", self.producer, merge.ReportEvidence(self.index)
            )

    def test_worker_success_flags_cannot_hide_failed_original_process(self):
        proof = merge.read(self.frozen_proof)
        row = proof["rows"][0]
        worker_path = Path(row["receipt_path"])
        worker = merge.read(worker_path)
        worker["supervision"]["exit_code"] = 1
        worker_path.write_text(json.dumps(worker))
        row["receipt_sha256"] = merge.digest(worker_path)
        self.frozen_proof.write_text(json.dumps(proof))
        workload = merge.read(self.frozen_workload)
        workload["owner_proof"] = merge.evidence(self.frozen_proof)
        self.frozen_workload.write_text(json.dumps(workload))
        relation = merge.read(self.relation_path)
        relation.update(
            frozen_owner_proof=merge.evidence(self.frozen_proof),
            frozen_workload=merge.evidence(self.frozen_workload),
        )
        self.relation_path.write_text(json.dumps(relation))
        self.refresh()
        with self.assertRaisesRegex(ValueError, "worker receipt belongs"):
            self.qualify()


class LightMetadataRequalificationTests(unittest.TestCase):
    """Opaque synthetic LIGHT receipts; no native cache export or evaluation."""

    def setUp(self):
        self.binding = {
            "family_fingerprint": "TEST_FAMILY",
            "physical_arity": 15,
            "artifact_digest": [1] * 32,
            "layer": "saved-application-route-apply-v1",
        }
        # Compose the Native fixture, without inheriting its snapshot-only tests.
        self.fixture = MetadataRequalificationTests()
        self.fixture.cache_fixture = self.cache_fixture
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.target = self.fixture.target
        self.cache = self.target / "5L/cache/prepared.rrcache.gz"
        self.metadata = self.target / "5L/cache.json"
        self.loop_path = self.target / "5L/manifest.json"
        old_root = merge.read(self.fixture.new_manifest)
        selected = old_root.pop("runtime_cache_snapshots")[0]
        (self.target / selected["manifest"]).unlink()
        loop = merge.read(self.loop_path)
        del loop["files"][selected["manifest"].removeprefix("5L/")]
        self.put(self.fixture.new_manifest, old_root)
        self.put(self.loop_path, loop)
        self.cache.write_bytes(b"OPAQUE SYNTHETIC LIGHT CACHE; NOT A NATIVE FILE")
        producer = self.fixture.producer["raw_report"]
        self.assertEqual(producer["semantic_binding"], self.binding)
        self.source_record = {
            "producer": self.provenance(
                "producer.json",
                Path(self.fixture.producer["evidence"]["path"]).read_bytes(),
            ),
            "data_manifest_sha256": merge.digest(self.fixture.old_manifest),
            "cache_native_index_identity": producer["cache_native_index_identity"],
            "inventory": self.provenance(
                "source-inventory.json",
                {"synthetic": True, "inventory": producer["inventory"]},
            ),
        }
        self.selection = {
            "keys": self.provenance(
                "selected-keys.json", {"synthetic": True, "keys": ["ROUTE", "APPLY"]}
            ),
            "unique_keys": 2,
            "usage_profile": self.provenance(
                "usage-profile.json",
                {"synthetic": True, "tuning": 1, "physical-example": 3},
            ),
            "allowed_input_classes": ["tuning", "physical-example"],
            "final_held_out_used": False,
        }
        self.identity = {
            "compressed_bytes": self.cache.stat().st_size,
            "compressed_blake3": [4] * 32,
            "decoded_bytes": 512,
            "decoded_blake3": [5] * 32,
        }
        self.decoded_sha256 = "f" * 64
        exporter, exporter_artifact, exporter_sha = self.tool("exporter", "8")
        reader, reader_artifact, reader_sha = self.tool("reader", "9")
        self.validation_fields = {
            "selected_records": 2,
            "reimported_records": 2,
            "complete_records": 2,
            "missing_keys": 0,
            "partial_records": 0,
            "exact_record_and_ordered_guard_equality": True,
            "declared_terminal_validation": True,
            "source_unchanged": True,
            "codec_state_format": 6,
        }
        self.raw_path = self.target / "provenance/native-validation.json"
        self.put(
            self.raw_path,
            {
                "schema": "vakint.strict-light-cache-native-validation.v1",
                "success": True,
                "semantic_binding": self.binding,
                "cache_sha256": merge.digest(self.cache),
                "cache_identity": self.identity,
                "decoded_sha256": self.decoded_sha256,
                "selected_keys_sha256": self.selection["keys"]["sha256"],
                "source_inventory_sha256": self.source_record["inventory"]["sha256"],
                "source_cache_native_index_identity": self.source_record[
                    "cache_native_index_identity"
                ],
                "usage_profile_sha256": self.selection["usage_profile"]["sha256"],
                "compiled_cohort": reader,
                "collector_artifacts": reader_artifact,
                "executable_sha256": reader_sha,
                **self.validation_fields,
            },
        )
        self.put(
            self.metadata,
            {
                "schema": "vakint.external-data.saved-cache.v1",
                "strict": True,
                "coverage_claim": False,
                "selection": {
                    "kind": "usage-selected-strict-snapshot",
                    "source": self.source_record,
                    **self.selection,
                    "validation": self.ref(self.raw_path),
                },
            },
        )
        # A synthetic file-contract witness, not an executed source admission.
        self.admission_path = self.target / "provenance/selection-admission.json"
        self.put(
            self.admission_path,
            {
                "schema": "vakint.light-cache-selection-admission.v1",
                "accepted": True,
                "cache_metadata_sha256": merge.digest(self.metadata),
                "selected_keys_sha256": self.selection["keys"]["sha256"],
                "usage_profile_sha256": self.selection["usage_profile"]["sha256"],
                "source_inventory_sha256": self.source_record["inventory"]["sha256"],
                "native_validation_sha256": merge.digest(self.raw_path),
                "semantic_binding": self.binding,
                "source_cache_native_index_identity": self.source_record[
                    "cache_native_index_identity"
                ],
                "checks": {
                    key: True
                    for key in (
                        "authenticated_eligible_sources",
                        "dispatch_counts_reconciled",
                        "positive_frequency_selected_subset",
                        "final_held_out_excluded",
                    )
                },
                "selected_unique_keys": 2,
                "eligible_profile_unique_keys": 3,
            },
        )
        self.export_path = self.target / "provenance/light-export.json"
        self.put(
            self.export_path,
            {
                "schema": "vakint.strict-light-cache-export-validation.v1",
                "accepted": True,
                "source": self.source_record,
                "semantic_binding": self.binding,
                "selection": self.selection,
                "output": {
                    "cache": self.ref(self.cache, size=True),
                    "metadata": self.ref(self.metadata, size=True),
                    "identity": self.identity,
                    "decoded_sha256": self.decoded_sha256,
                },
                "compiled_cohort": exporter,
                "collector_artifacts": exporter_artifact,
                "executable_sha256": exporter_sha,
                "validation": {
                    "receipt": self.ref(self.raw_path),
                    **self.validation_fields,
                },
                "selection_admission": self.ref(self.admission_path),
            },
        )
        loop = merge.read(self.loop_path)
        loop["cache"].update(identity=self.identity, decoded_sha256=self.decoded_sha256)
        self.put(self.loop_path, loop)
        relation = merge.read(self.fixture.relation_path)
        relation.pop("runtime_cache_snapshot")
        relation["prepared_cache_replacement"] = {
            "loop_count": 5,
            "export_validation": merge.evidence(self.export_path),
        }
        self.put(self.fixture.relation_path, relation)
        self.refresh()

    @staticmethod
    def put(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, indent=2))

    def ref(self, path, *, size=False):
        result = {
            "path": path.relative_to(self.target).as_posix(),
            "sha256": merge.digest(path),
        }
        if size:
            result["bytes"] = path.stat().st_size
        return result

    def provenance(self, name, document):
        path = self.target / "provenance" / name
        if isinstance(document, bytes):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(document)
        else:
            self.put(path, document)
        return self.ref(path)

    def cache_fixture(self, count, directory, files):
        if count != 5:
            return {"path": "cache/prepared.gz"}
        prepared = directory / "cache/prepared.rrcache.gz"
        (directory / "cache/prepared.gz").rename(prepared)
        files["cache/prepared.rrcache.gz"] = files.pop("cache/prepared.gz")
        metadata = directory / "cache.json"
        self.put(metadata, {"synthetic": "ORIGINAL PREPARED CACHE METADATA"})
        files["cache.json"] = {
            "sha256": merge.digest(metadata),
            "bytes": metadata.stat().st_size,
        }
        return {
            "available": True,
            "path": "cache/prepared.rrcache.gz",
            "metadata": "cache.json",
            "format": "rustred.saved-application-cache-file.v1",
            "binding": self.binding,
            "identity": {
                "compressed_bytes": prepared.stat().st_size,
                "compressed_blake3": [1] * 32,
                "decoded_bytes": 256,
                "decoded_blake3": [2] * 32,
            },
            "decoded_sha256": "c" * 64,
        }

    def tool(self, role, digit):
        cohort = copy.deepcopy(merge.read(self.fixture.frozen_proof)["compiled_cohort"])
        cohort["main_revision"] = digit * 40
        binary = self.target / "provenance" / f"{role}.bin"
        binary.parent.mkdir(parents=True, exist_ok=True)
        binary.write_text(f"SYNTHETIC {role}; NOT AN EXECUTABLE")
        artifact = copy.deepcopy(merge.read(self.fixture.base.artifact))
        artifact["stamps"]["VAKINT_CACHE_BUILD_GIT_REVISION"] = cohort["main_revision"]
        artifact["binary"] = {
            "sha256": merge.digest(binary),
            "immutable_path": str(binary),
        }
        return (
            cohort,
            self.provenance(f"{role}-artifacts.json", artifact),
            merge.digest(binary),
        )

    def refresh(self, *, admission_links=True):
        """Rebind changed opaque files without repairing their semantic contents."""
        loop = merge.read(self.loop_path)
        for path in (self.cache, self.metadata):
            loop["files"][path.relative_to(self.target / "5L").as_posix()] = {
                "sha256": merge.digest(path),
                "bytes": path.stat().st_size,
            }
        self.put(self.loop_path, loop)
        exported = merge.read(self.export_path)
        exported["output"]["metadata"] = self.ref(self.metadata, size=True)
        if admission_links:
            admission = merge.read(self.admission_path)
            admission.update(
                cache_metadata_sha256=merge.digest(self.metadata),
                selected_keys_sha256=exported["selection"]["keys"]["sha256"],
                usage_profile_sha256=exported["selection"]["usage_profile"]["sha256"],
                source_inventory_sha256=exported["source"]["inventory"]["sha256"],
                native_validation_sha256=exported["validation"]["receipt"]["sha256"],
                semantic_binding=exported["semantic_binding"],
                source_cache_native_index_identity=exported["source"][
                    "cache_native_index_identity"
                ],
            )
            self.put(self.admission_path, admission)
        exported["selection_admission"] = self.ref(self.admission_path)
        self.put(self.export_path, exported)
        relation = merge.read(self.fixture.relation_path)
        if "prepared_cache_replacement" in relation:
            relation["prepared_cache_replacement"]["export_validation"] = (
                merge.evidence(self.export_path)
            )
        self.put(self.fixture.relation_path, relation)
        self.fixture.refresh()

    def alter_export(self, change):
        document = merge.read(self.export_path)
        change(document)
        self.put(self.export_path, document)
        self.refresh()

    def alter_raw(self, change):
        document = merge.read(self.raw_path)
        change(document)
        self.put(self.raw_path, document)
        reference = self.ref(self.raw_path)
        metadata = merge.read(self.metadata)
        metadata["selection"]["validation"] = reference
        self.put(self.metadata, metadata)
        exported = merge.read(self.export_path)
        exported["validation"]["receipt"] = reference
        self.put(self.export_path, exported)
        self.refresh()

    def fresh(self):
        case = LightMetadataRequalificationTests()
        case.setUp()
        self.addCleanup(case.doCleanups)
        return case

    def qualify(self):
        return self.fixture.qualify()["data_metadata_requalification"]

    def test_light_replacement_preserves_science_and_distinct_tool_provenance(self):
        producer = copy.deepcopy(self.fixture.producer["raw_report"])
        qualified = self.qualify()
        self.assertEqual(self.fixture.producer["raw_report"], producer)
        self.assertEqual(
            qualified["source_manifest"], merge.evidence(self.fixture.old_manifest)
        )
        self.assertEqual(
            qualified["consumer_manifest"], merge.evidence(self.fixture.new_manifest)
        )
        self.assertNotEqual(
            merge.read(self.export_path)["compiled_cohort"],
            merge.read(self.raw_path)["compiled_cohort"],
        )
        self.assertEqual(
            merge.read(self.fixture.new_manifest).get("runtime_cache_snapshots", []),
            [],
        )

    def test_light_requires_exactly_one_relation_mode(self):
        for change in ("missing", "both", "wrong-loop", "extra-field"):
            with self.subTest(change=change):
                case = self.fresh()
                relation = merge.read(case.fixture.relation_path)
                if change == "missing":
                    relation.pop("prepared_cache_replacement")
                elif change == "both":
                    relation["runtime_cache_snapshot"] = {"loop_count": 5}
                elif change == "wrong-loop":
                    relation["prepared_cache_replacement"]["loop_count"] = 4
                else:
                    relation["prepared_cache_replacement"]["extra"] = True
                case.put(case.fixture.relation_path, relation)
                case.fixture.refresh()
                with self.assertRaises(ValueError):
                    case.qualify()

    def test_light_cannot_be_shadowed_by_a_runtime_snapshot(self):
        for key in ("source_manifest", "consumer_manifest"):
            with self.subTest(root=key):
                case = self.fresh()
                if key == "source_manifest":
                    path = case.fixture.old_manifest
                else:
                    path = case.fixture.new_manifest
                root = merge.read(path)
                root["runtime_cache_snapshots"] = [{"loop_count": 5}]
                case.put(path, root)
                case.fixture.refresh()
                with self.assertRaises(ValueError):
                    case.qualify()

    def test_light_declared_size_cap_and_native_identity_are_checked(self):
        for field, value in (
            ("bytes", 0),
            ("bytes", True),
            ("bytes", 95_000_001),
            ("bytes", 1),
            ("sha256", "0" * 64),
        ):
            with self.subTest(field=field, value=value):
                case = self.fresh()
                case.alter_export(
                    lambda doc, field=field, value=value: doc["output"]["cache"].update(
                        {field: value}
                    )
                )
                with self.assertRaises(ValueError):
                    case.qualify()
        for field, value in (
            ("compressed_bytes", 1),
            ("compressed_blake3", "not-a-digest"),
            ("decoded_bytes", 0),
            ("decoded_blake3", "0" * 64),
        ):
            with self.subTest(identity=field):
                case = self.fresh()
                case.alter_export(
                    lambda doc, field=field, value=value: doc["output"][
                        "identity"
                    ].update({field: value})
                )
                with self.assertRaises(ValueError):
                    case.qualify()

    def test_light_rejects_oversized_actual_file_before_hashing_it(self):
        with self.cache.open("r+b") as stream:
            stream.truncate(95_000_001)
        exported = merge.read(self.export_path)
        exported["output"]["cache"]["bytes"] = 95_000_001
        exported["output"]["identity"]["compressed_bytes"] = 95_000_001
        self.put(self.export_path, exported)
        relation = merge.read(self.fixture.relation_path)
        relation["prepared_cache_replacement"]["export_validation"] = merge.evidence(
            self.export_path
        )
        self.put(self.fixture.relation_path, relation)
        self.fixture.refresh()
        original_digest = merge.digest

        def bounded_digest(path):
            self.assertNotEqual(Path(path).resolve(), self.cache.resolve())
            return original_digest(path)

        with (
            patch.object(merge, "digest", bounded_digest),
            self.assertRaisesRegex(ValueError, "95,000,000"),
        ):
            self.qualify()

    def test_light_accepts_the_inclusive_compressed_byte_cap(self):
        # A sparse opaque fixture exercises the real filesystem size boundary.
        # It is never imported or presented as an actual Native cache.
        with self.cache.open("r+b") as stream:
            stream.truncate(95_000_000)
        identity = {**self.identity, "compressed_bytes": 95_000_000}
        cache_reference = self.ref(self.cache, size=True)
        loop = merge.read(self.loop_path)
        loop["cache"]["identity"] = identity
        self.put(self.loop_path, loop)
        self.alter_raw(
            lambda doc: doc.update(
                cache_sha256=cache_reference["sha256"], cache_identity=identity
            )
        )
        self.alter_export(
            lambda doc: doc["output"].update(cache=cache_reference, identity=identity)
        )
        self.qualify()

    def test_light_requires_complete_matching_native_validation(self):
        for field, value in (
            ("selected_records", 0),
            ("selected_records", True),
            ("reimported_records", 1),
            ("complete_records", 1),
            ("missing_keys", 1),
            ("partial_records", 1),
            ("exact_record_and_ordered_guard_equality", False),
            ("declared_terminal_validation", False),
            ("source_unchanged", False),
            ("codec_state_format", 5),
        ):
            with self.subTest(field=field):
                case = self.fresh()
                case.alter_export(
                    lambda doc, field=field, value=value: doc["validation"].update(
                        {field: value}
                    )
                )
                case.alter_raw(
                    lambda doc, field=field, value=value: doc.update({field: value})
                )
                with self.assertRaises(ValueError):
                    case.qualify()

    def test_light_raw_receipt_cannot_be_forged_by_wrapper_success(self):
        for field, value in (
            ("schema", "OTHER"),
            ("success", False),
            ("cache_sha256", "0" * 64),
            ("decoded_sha256", "0" * 64),
            ("selected_keys_sha256", "0" * 64),
            ("source_inventory_sha256", "0" * 64),
            ("usage_profile_sha256", "0" * 64),
            ("executable_sha256", "0" * 64),
        ):
            with self.subTest(field=field):
                case = self.fresh()
                case.alter_raw(
                    lambda doc, field=field, value=value: doc.update({field: value})
                )
                with self.assertRaises(ValueError):
                    case.qualify()
        case = self.fresh()
        other_source = case.root / "OTHER-SYNTHETIC-SOURCE"
        other_source.mkdir()
        (other_source / "cache.rrindex").write_text("SYNTHETIC DIFFERENT SOURCE INDEX")
        other_identity = merge.native_cache_identity(other_source)
        case.alter_raw(
            lambda doc: doc.update(source_cache_native_index_identity=other_identity)
        )
        with self.assertRaises(ValueError):
            case.qualify()
        self.raw_path.write_text("TRUNCATED JSON")
        with self.assertRaises(ValueError):
            self.qualify()

    def test_light_native_binding_keeps_json_integer_types(self):
        for location in (
            "raw-binding",
            "raw-identity",
            "loop-binding",
            "loop-identity",
        ):
            with self.subTest(location=location):
                case = self.fresh()

                def change(document, location=location):
                    if location == "raw-binding":
                        document["semantic_binding"]["artifact_digest"][0] = True
                    elif location == "raw-identity":
                        document["cache_identity"]["compressed_blake3"][0] = 4.0
                    elif location == "loop-binding":
                        document["cache"]["binding"]["artifact_digest"][0] = True
                    else:
                        identity = document["cache"]["identity"]
                        identity["compressed_bytes"] = float(
                            identity["compressed_bytes"]
                        )

                if location.startswith("raw-"):
                    case.alter_raw(change)
                else:
                    case.fixture.change_loop(change)
                with self.assertRaises(ValueError):
                    case.qualify()

    def test_light_selection_admission_is_required_and_binds_actual_sources(self):
        mutations = [
            lambda d: d.update(schema="OTHER"),
            lambda d: d.update(accepted=False),
            lambda d: d["semantic_binding"].update(physical_arity=14),
            lambda d: d["source_cache_native_index_identity"].update(sha256="0" * 64),
        ]
        for field in (
            "cache_metadata_sha256",
            "selected_keys_sha256",
            "usage_profile_sha256",
            "source_inventory_sha256",
            "native_validation_sha256",
        ):
            mutations.append(lambda d, field=field: d.update({field: "0" * 64}))
        for field, value in (
            ("selected_unique_keys", 0),
            ("selected_unique_keys", True),
            ("selected_unique_keys", 3),
            ("eligible_profile_unique_keys", 0),
            ("eligible_profile_unique_keys", True),
            ("eligible_profile_unique_keys", 1),
        ):
            mutations.append(
                lambda d, field=field, value=value: d.update({field: value})
            )
        for field in (
            "authenticated_eligible_sources",
            "dispatch_counts_reconciled",
            "positive_frequency_selected_subset",
            "final_held_out_excluded",
        ):
            mutations.append(lambda d, field=field: d["checks"].update({field: False}))
        for index, change in enumerate(mutations):
            with self.subTest(mutation=index):
                case = self.fresh()
                admission = merge.read(case.admission_path)
                change(admission)
                case.put(case.admission_path, admission)
                case.refresh(admission_links=False)
                with self.assertRaises(ValueError):
                    case.qualify()
        exported = merge.read(self.export_path)
        exported.pop("selection_admission")
        self.put(self.export_path, exported)
        relation = merge.read(self.fixture.relation_path)
        relation["prepared_cache_replacement"]["export_validation"] = merge.evidence(
            self.export_path
        )
        self.put(self.fixture.relation_path, relation)
        self.fixture.refresh()
        with self.assertRaises(ValueError):
            self.qualify()

    def test_light_source_binding_selection_and_exporter_are_authenticated(self):
        mutations = (
            lambda d: d.update(accepted=False),
            lambda d: d["semantic_binding"].update(physical_arity=14),
            lambda d: d["source"].update(data_manifest_sha256="0" * 64),
            lambda d: d["source"]["cache_native_index_identity"].update(
                sha256="0" * 64
            ),
            lambda d: d["source"]["producer"].update(sha256="0" * 64),
            lambda d: d["selection"].update(unique_keys=1),
            lambda d: d["selection"].update(final_held_out_used=True),
            lambda d: d["selection"].update(allowed_input_classes=["held_out"]),
            lambda d: d["selection"].update(
                allowed_input_classes=["training", "physical-example"]
            ),
            lambda d: d["selection"].update(
                allowed_input_classes=["tuning", "physical-example", "physical-example"]
            ),
            lambda d: d.update(executable_sha256="0" * 64),
            lambda d: d["compiled_cohort"].update(main_revision="0" * 40),
        )
        for index, change in enumerate(mutations):
            with self.subTest(mutation=index):
                case = self.fresh()
                case.alter_export(change)
                with self.assertRaises(ValueError):
                    case.qualify()

    def test_light_only_allows_the_two_existing_five_loop_payloads(self):
        for count, name in (
            (1, "cache/prepared.gz"),
            (4, "cache/prepared.gz"),
            (5, "ibps/program.bin"),
            (5, "masters.json"),
            (5, "graphs/master.dot"),
            (5, "cache/extra.gz"),
        ):
            with self.subTest(loop=count, payload=name):
                case = self.fresh()
                path = case.target / f"{count}L" / name
                path.write_text("ILLICIT BUT REHASHED SYNTHETIC PAYLOAD")
                case.fixture.change_loop(
                    lambda doc, name=name, path=path: doc["files"].update(
                        {
                            name: {
                                "sha256": merge.digest(path),
                                "bytes": path.stat().st_size,
                            }
                        }
                    ),
                    count=count,
                )
                with self.assertRaises(ValueError):
                    case.qualify()

    def test_light_preserves_scientific_and_immutable_cache_fields(self):
        for field in ("normalization", "coverage", "masters", "ibps"):
            with self.subTest(scientific=field):
                case = self.fresh()
                case.fixture.change_loop(
                    lambda doc, field=field: doc.update({field: "ILLICIT"})
                )
                with self.assertRaises(ValueError):
                    case.qualify()
        for field, value in (
            ("available", False),
            ("available", 1),
            ("format", "other-format"),
            ("binding", {"other": True}),
            ("path", "cache/other.gz"),
            ("metadata", "other.json"),
            ("extra", True),
        ):
            with self.subTest(cache=field, value=value):
                case = self.fresh()
                case.fixture.change_loop(
                    lambda doc, field=field, value=value: doc["cache"].update(
                        {field: value}
                    )
                )
                with self.assertRaises(ValueError):
                    case.qualify()

    def test_light_metadata_must_match_export_and_keep_narrow_schema(self):
        mutations = (
            lambda d: d.update(strict=False),
            lambda d: d.update(strict=1),
            lambda d: d.update(coverage_claim=True),
            lambda d: d.update(coverage_claim=0),
            lambda d: d.update(extra=True),
            lambda d: d["selection"].update(kind="OTHER"),
            lambda d: d["selection"].update(unique_keys=3),
            lambda d: d["selection"].update(final_held_out_used=True),
            lambda d: d["selection"]["keys"].update(sha256="0" * 64),
            lambda d: d["selection"]["validation"].update(sha256="0" * 64),
            lambda d: d["selection"]["source"].update(data_manifest_sha256="0" * 64),
        )
        for index, change in enumerate(mutations):
            with self.subTest(mutation=index):
                case = self.fresh()
                document = merge.read(case.metadata)
                change(document)
                case.put(case.metadata, document)
                case.refresh()
                with self.assertRaises(ValueError):
                    case.qualify()

    def test_light_metadata_references_are_portable_root_provenance_paths(self):
        for kind in ("absolute", "parent", "dot", "parent-segment", "outside"):
            with self.subTest(kind=kind):
                case = self.fresh()
                keys = case.target / "provenance/selected-keys.json"
                if kind == "absolute":
                    path = str(keys)
                elif kind == "parent":
                    path = "../published-data/provenance/selected-keys.json"
                elif kind == "dot":
                    path = "./provenance/selected-keys.json"
                elif kind == "parent-segment":
                    path = "provenance/../provenance/selected-keys.json"
                else:
                    path = "selected-keys.json"
                    (case.target / path).write_bytes(keys.read_bytes())
                candidate = Path(path)
                candidate = (
                    candidate if candidate.is_absolute() else case.target / candidate
                )
                self.assertEqual(merge.digest(candidate), merge.digest(keys))
                document = merge.read(case.metadata)
                document["selection"]["keys"]["path"] = path
                case.put(case.metadata, document)
                case.refresh()
                with self.assertRaises(ValueError):
                    case.qualify()

    def test_light_duplicate_keys_symlinks_and_hidden_payloads_stay_rejected(self):
        for kind in ("duplicate", "symlink", "hidden"):
            with self.subTest(kind=kind):
                case = self.fresh()
                if kind == "duplicate":
                    case.metadata.write_text(
                        '{"strict":false,' + case.metadata.read_text()[1:]
                    )
                    case.refresh()
                elif kind == "symlink":
                    case.cache.unlink()
                    case.cache.symlink_to(
                        case.fixture.source / "5L/cache/prepared.rrcache.gz"
                    )
                else:
                    (case.target / "5L/hidden.bin").write_text("HIDDEN")
                with self.assertRaises(ValueError):
                    case.qualify()

    def test_light_release_provenance_and_frozen_workload_remain_exact(self):
        for change in ("provenance", "version", "archive", "workload"):
            with self.subTest(change=change):
                case = self.fresh()
                if change == "workload":
                    workload = merge.read(case.fixture.consumer_workload)
                    workload["held_out"][0]["powers"][0] += 1
                    case.put(case.fixture.consumer_workload, workload)
                else:
                    root = merge.read(case.fixture.new_manifest)
                    if change == "provenance":
                        root["provenance_assets"] = [{"archive": "ILLICIT"}]
                    elif change == "version":
                        root["version"] = "not-a-release"
                    else:
                        root["loop_packs"][4]["archive"] = "wrong.tar.gz"
                    case.put(case.fixture.new_manifest, root)
                case.fixture.refresh()
                with self.assertRaises(ValueError):
                    case.qualify()


class QualificationProjectionTests(unittest.TestCase):
    """Synthetic metadata only: no actual native success or acceptance is issued."""

    def setUp(self):
        self.make_fixture()

    @staticmethod
    def key(index):
        return [index] + [1] * 14

    def make_fixture(self):
        # Reuse the existing strict bridge fixture without inheriting/rerunning
        # its sixty tests. Every payload remains in a synthetic temporary tree.
        self.base = ReportContractTests()
        self.base.setUp()
        self.addCleanup(self.base.doCleanups)
        self.identity, self.producer, self.accepted_path, self.index = (
            self.base.cross_cohort_fixture()
        )
        self.accepted = merge.read(self.accepted_path)
        self.binding = self.accepted["semantic_binding"]
        self.docs, self.paths = {}, {}
        self.source_files = []
        source_hashes = {}
        for filename in (
            "measure_cache.rs",
            "build_cache/workload.rs",
            "build_cache/owner_proof.rs",
        ):
            name = "crates/vakint/examples/" + filename
            path = self.base.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("SYNTHETIC WORKER SOURCE ONLY: " + filename)
            self.source_files.append(merge.evidence(path))
            source_hashes[name] = merge.digest(path)
        cohort = merge.read(self.accepted["consumer"]["compiled_cohort"]["path"])
        self.put(
            "snapshot",
            {
                "test_fixture": True,
                "head": cohort["main_revision"],
                "source_sha256": source_hashes,
            },
        )
        self.save("snapshot")
        cohort["main_source_snapshot_sha256"] = merge.digest(self.paths["snapshot"])
        self.put("cohort", cohort)
        self.put("artifact", merge.read(self.base.artifact))
        self.paths["artifact"] = self.base.artifact
        self.docs["artifact"]["stamps"].update(
            {stamp: cohort[field] for field, stamp in merge.COHORT_STAMPS}
        )
        self.docs["artifact"]["main_source_snapshot"] = self.reference("snapshot")

        def training(index):
            return {
                "id": f"SYNTHETIC_TRAINING_{index}",
                "powers": self.key(10 + index),
                "public_input": None,
                "rationale": "SYNTHETIC UNEXECUTED TRAINING METADATA",
                "lowered_roots": [],
                "protected_roots": [],
            }

        public = []
        for index in range(13):
            public.append(
                {
                    "id": f"SYNTHETIC_PUBLIC_{index}",
                    "powers": self.key(100 + index),
                    "public_input": f"SYNTHETIC PUBLIC INPUT {index}",
                    "rationale": "SYNTHETIC FILE CONTRACT ONLY",
                    "lowered_roots": [self.key(1000 + index)],
                    "protected_roots": [self.key(2000 + index)]
                    + ([self.key(3000 + index)] if index < 12 else []),
                }
            )
        selection = {
            "schema": "vakint.five-loop-cache-workload.v1",
            "data_manifest_sha256": "7" * 64,
            "source_census_sha256": "9" * 64,
            "training": [training(0), training(1)],
            "tuning": copy.deepcopy(public[:1]),
            "held_out": copy.deepcopy(public[1:]),
        }
        self.put("target", copy.deepcopy(selection))
        self.put("producer_workload", copy.deepcopy(selection))
        source = copy.deepcopy(selection)
        source["training"].append(training(2))
        source["owner_proof"] = None
        for root in source["tuning"] + source["held_out"]:
            root["lowered_roots"], root["protected_roots"] = [], []
        self.put("qualification", source)
        self.put("proof", merge.read(self.accepted["consumer"]["owner_proof"]["path"]))
        self.docs["proof"]["compiled_cohort"] = cohort
        self.docs["proof"]["rows"] = []
        details = []
        for index, root in enumerate(public):
            worker_path = (
                self.base.root / f"projection-current-{index}/measurement.json"
            )
            old_path = (
                self.base.root / f"projection-historical-{index}/measurement.json"
            )
            for path in (worker_path, old_path):
                path.parent.mkdir()
                (path.parent / "value.symbolica").write_text(
                    f"SYNTHETIC EXACT VALUE {index}"
                )
                (path.parent / "guards.json").write_text(
                    json.dumps([f"SYNTHETIC COMPLETE GUARD {index}"])
                )
                (path.parent / "unknown-masters.json").write_text("[]")
            identity = {
                "data_manifest_sha256": "7" * 64,
                "source_census_sha256": "9" * 64,
                "build_git_revision": cohort["main_revision"],
                "executable_sha256": merge.digest(self.base.binary),
                "cache_sha256": None,
                **{
                    field: source_hashes["crates/vakint/examples/" + filename]
                    for field, filename in (
                        ("compiled_driver_sha256", "measure_cache.rs"),
                        ("compiled_workload_schema_sha256", "build_cache/workload.rs"),
                        ("compiled_owner_proof_sha256", "build_cache/owner_proof.rs"),
                    )
                },
            }
            worker = {
                "test_fixture": True,
                "success": True,
                "worker_success": True,
                "fresh_process": True,
                "variant": "uncached",
                "sample": 0,
                "case_id": root["id"],
                "declared_common_powers": root["powers"],
                "identity": identity,
                "actual_lowered_roots": root["lowered_roots"],
                "supervision": {
                    "success": True,
                    "exit_code": 0,
                    "deadline_exhausted": False,
                },
                "persistent_cache_attachments_before": {
                    str(loop): False for loop in range(1, 6)
                },
                "persistent_cache_attachments_during_evaluation": {
                    str(loop): False for loop in range(1, 6)
                },
                "cache_statistics": {"binding": self.binding},
                "value_sha256": merge.digest(worker_path.with_name("value.symbolica")),
                "guards_sha256": merge.digest(worker_path.with_name("guards.json")),
                "catalogue_closure": {
                    "closed": True,
                    "unknown_count": 0,
                    "unknown_sha256": merge.digest(
                        worker_path.with_name("unknown-masters.json")
                    ),
                },
                "evaluation_seconds": 12345,
            }
            self.put(f"worker{index}", worker, worker_path)
            self.put(f"historical{index}", copy.deepcopy(worker), old_path)
            self.docs["proof"]["rows"].append(
                {
                    "id": root["id"],
                    "logical_cohort": "tuning" if index == 0 else "held-out",
                    "public_input_source_sha256": hashlib.sha256(
                        root["public_input"].encode()
                    ).hexdigest(),
                    "application_binding": self.binding,
                    "reduction_succeeded": True,
                    "complete_owner_transport": True,
                    "mapping_error": None,
                    "actual_lowered_roots": root["lowered_roots"],
                    "protected_owner_roots": root["protected_roots"],
                }
            )
            details.append(
                {
                    "id": root["id"],
                    "role": "tuning" if index == 0 else "held_out",
                    "_fixture_worker": index,
                }
            )
        forbidden = sorted(
            {
                tuple(key)
                for root in public
                for key in [
                    root["powers"],
                    *root["lowered_roots"],
                    *root["protected_roots"],
                ]
            }
        )
        self.assertEqual(len(forbidden), 51)
        self.put(
            "union",
            {"forbidden_direct_target_or_owner_keys": forbidden, "forbidden_count": 51},
        )
        self.put(
            "routes",
            {
                "semantic_binding": self.binding,
                "compiled_cohort": cohort,
                "source_census_sha256": "9" * 64,
                "rows": [
                    {
                        "raw_powers": training(index)["powers"],
                        "route_success": True,
                        "route_error": None,
                        "owner_endpoints": [
                            {"powers": self.key(7000 + index), "coefficient": "1"}
                        ],
                    }
                    for index in range(3)
                ],
            },
        )
        self.put(
            "aliases",
            {
                "rows": [
                    {
                        "powers": training(index)["powers"],
                        "eligible_after_external_alias_exclusion": True,
                        "raw_protected_intersection": False,
                        "owner_protected_intersections": [],
                    }
                    for index in range(3)
                ]
            },
        )
        self.put(
            "validator",
            {
                "valid": True,
                "public_cases_checked": 13,
                "data_manifest_sha256": "7" * 64,
                "reduction_or_cache_import_performed": False,
                "identity": {"compiled_cohort": cohort},
            },
        )
        cached_path = self.base.root / "projection-cached/worker.json"
        cached_path.parent.mkdir()
        for name in ("value.symbolica", "guards.json"):
            cached_path.with_name(name).write_bytes(
                self.paths["worker0"].with_name(name).read_bytes()
            )
        cached = copy.deepcopy(self.docs["worker0"])
        cached.update(variant="cached", native_keys_authenticated=True)
        cached["identity"] = {
            "build_git_revision": cohort["main_revision"],
            "data_manifest_sha256": "7" * 64,
            "source_census_sha256": "9" * 64,
            "executable_sha256": merge.digest(self.base.binary),
            "cache_sha256": merge.digest(self.base.cache),
        }
        cached["persistent_cache_attachments_during_evaluation"]["5"] = True
        cached["cache_statistics"]["file"] = {"records": 3}
        self.put("cached", cached, cached_path)
        self.put("import_supervision", {"exit_code": 0, "deadline_exhausted": False})
        self.put(
            "import",
            {
                "checks": {
                    key: True
                    for key in (
                        "compressed_and_decoded_identity",
                        "single_member_and_framing",
                        "native_binding_and_every_record_validation",
                        "fresh_cached_case_exact_value_and_all_guards",
                        "catalogue_closure",
                    )
                },
                "actual_file_statistics": {"records": 3},
                "producer_records": 3,
            },
        )
        self.put(
            "failed",
            {
                "complete": False,
                "accepted": False,
                "counts": {
                    "functional_success": 10,
                    "deadline_failed": 1,
                    "not_started_after_failure": 2,
                    "started_native_workers": 11,
                    "planned_public_cases": 13,
                },
                "failed_case": {"exit_code": 124, "worker_json_exists": False},
                "mapper_launched": False,
                "performance_eligible": False,
            },
        )
        self.put(
            "detail",
            {
                "schema": "vakint.qualification-workload-projection-evidence.v1",
                "fresh_training_reductions_executed": 0,
                "full_case_key_value_and_guard_evidence": details,
            },
        )
        self.put("native_reviewed_roles", {})
        self.put("packet", {})
        self.put(
            "review",
            {
                "schema": "vakint.independent-fallback2-semantic-projection-review.v1",
                "verdict": "SOUND_FOR_NATIVE_COMPATIBILITY_WITH_EXPLICIT_AUTHENTICATED_PROJECTION_NOT_A_NEW_F022_SUITE_PASS",
                "source_files": self.source_files,
            },
        )
        self.put(
            "scope",
            {"schema": "vakint.independent-projection-validator-scope-addendum.v1"},
        )
        self.put(
            "relation",
            {
                "schema": "vakint.qualification-workload-projection.v1",
                "relation": "unexecuted-training-prefix-pruning",
            },
        )
        self.accepted["checks"]["unchanged_selection"] = False
        self.rebind()

    def put(self, name, value, path=None):
        self.docs[name] = value
        self.paths[name] = path or self.base.root / f"projection-{name}.json"

    def save(self, name):
        self.paths[name].write_text(json.dumps(self.docs[name], indent=2))

    def reference(self, name):
        return merge.evidence(self.paths[name])

    def rebind(self):
        """Refresh only file/cohort links, leaving scientific payload claims editable."""
        for name in ("snapshot", "cohort", "artifact", "qualification"):
            self.save(name)
        source_hash = merge.digest(self.paths["qualification"])
        for index in range(13):
            self.docs[f"worker{index}"]["identity"]["workload_sha256"] = source_hash
            self.save(f"worker{index}")
            self.save(f"historical{index}")
        for row in self.docs["proof"]["rows"]:
            index = int(row["id"].rsplit("_", 1)[1])
            row["receipt_sha256"] = merge.digest(self.paths[f"worker{index}"])
            row["workload_sha256"] = source_hash
        self.save("proof")
        self.docs["target"]["owner_proof"] = self.reference("proof")
        self.docs["producer_workload"]["owner_proof"] = self.accepted["producer"][
            "owner_proof"
        ]
        self.save("target")
        self.save("producer_workload")
        target_hash = merge.digest(self.paths["target"])
        self.docs["cached"]["identity"].update(
            workload_sha256=target_hash,
            owner_proof={
                **self.reference("proof"),
                "compiled_cohort": self.docs["cohort"],
            },
        )
        self.docs["cached"]["owner_proof"] = {
            **self.reference("proof"),
            "compiled_cohort": self.docs["cohort"],
        }
        self.save("cached")
        self.save("import_supervision")
        self.docs["import"].update(
            actual_cached_worker=self.reference("cached"),
            supervision=self.reference("import_supervision"),
        )
        self.save("import")
        self.docs["validator"].update(workload_sha256=target_hash)
        self.docs["validator"]["identity"]["sha256"] = merge.digest(self.paths["proof"])
        self.save("validator")
        self.docs["union"]["owner_proof"] = self.reference("proof")
        self.docs["union"]["compiled_cohort"] = self.reference("cohort")
        self.save("union")
        self.save("routes")
        self.docs["aliases"].update(
            qualification=self.reference("routes"),
            protected_union=self.reference("union"),
            canonical_cohort_sha256=merge.digest(self.paths["cohort"]),
        )
        self.save("aliases")
        self.accepted["producer"]["workload"] = self.reference("producer_workload")
        self.accepted["consumer"].update(
            collector_artifacts=self.reference("artifact"),
            compiled_cohort=self.reference("cohort"),
            qualification_workload=self.reference("qualification"),
            final_workload=self.reference("target"),
            owner_proof=self.reference("proof"),
            protected_union=self.reference("union"),
            training_route_qualification=self.reference("routes"),
            alias_exclusion=self.reference("aliases"),
        )
        build = self.producer["raw_report"]
        build["workload_sha256"] = merge.digest(self.paths["producer_workload"])
        path = Path(self.producer["evidence"]["path"])
        path.write_text(json.dumps(build))
        self.producer["evidence"] = merge.evidence(path)
        self.accepted["producer"]["build_report"] = merge.evidence(path)
        self.identity.update(
            workload_sha256=target_hash,
            owner_proof={
                **self.reference("proof"),
                "compiled_cohort": self.docs["cohort"],
            },
        )
        self.accepted["immutable_selection"] = merge.selected_inputs(
            self.docs["target"]
        )
        detail = self.docs["detail"]
        detail.update(
            source_fresh_qualification_workload=self.reference("qualification"),
            compatibility_target_workload=self.reference("target"),
            producer_workload=self.reference("producer_workload"),
            canonical_proof=self.reference("proof"),
            compiled_cohort=self.reference("cohort"),
            actual_existing_validator_two_root_result=self.reference("validator"),
            strict_import_equality_validation=self.reference("import"),
            two_root_cached_public_case=self.reference("cached"),
        )
        for row in detail["full_case_key_value_and_guard_evidence"]:
            index = row["_fixture_worker"]
            row.update(
                actual_qualification_receipt=self.reference(f"worker{index}"),
                original_producer_receipt=self.reference(f"historical{index}"),
            )
        self.save("detail")
        self.docs["native_reviewed_roles"].update(
            producer=copy.deepcopy(self.accepted["producer"]),
            consumer=copy.deepcopy(self.accepted["consumer"]),
        )
        self.save("native_reviewed_roles")
        self.docs["packet"].update(
            exact_semantic_projection_evidence=self.reference("detail"),
            native_compatibility_draft=self.reference("native_reviewed_roles"),
        )
        self.save("packet")
        self.save("failed")
        self.docs["review"].update(
            evidence_packet=self.reference("packet"),
            new_attempt_preserved={"postflight": self.reference("failed")},
        )
        self.save("review")
        self.docs["scope"]["review"] = self.reference("review")
        self.save("scope")
        self.docs["relation"].update(
            evidence=self.reference("detail"),
            independent_review=self.reference("review"),
            validator_scope=self.reference("scope"),
            incomplete_attempt=self.reference("failed"),
        )
        self.save("relation")
        self.accepted["qualification_projection"] = self.reference("relation")
        self.base.refresh_acceptance(self.accepted_path, self.index, self.accepted)

    def evaluate(self):
        return merge.ReportEvidence(self.index).cache_qualification(
            self.identity, self.producer
        )

    def assert_rejected(self):
        with self.assertRaises(ValueError):
            self.evaluate()

    def test_explicit_projection_keeps_source_and_target_roles_distinct(self):
        result = self.evaluate()
        projection = result["qualification_projection"]
        self.assertIn("Native compatibility only", projection["scope"])
        self.assertEqual(set(projection), {"evidence", "raw_relation", "scope"})
        self.assertNotEqual(
            result["raw_receipt"]["consumer"]["qualification_workload"],
            result["raw_receipt"]["consumer"]["final_workload"],
        )

    def test_missing_projection_preserves_the_strict_default(self):
        del self.accepted["qualification_projection"]
        self.accepted["checks"]["unchanged_selection"] = True
        self.base.refresh_acceptance(self.accepted_path, self.index, self.accepted)
        self.assert_rejected()

    def test_projection_requires_truthful_false_full_selection_flag(self):
        self.accepted["checks"]["unchanged_selection"] = True
        self.rebind()
        self.assert_rejected()

    def test_training_metadata_must_be_unexecuted(self):
        self.docs["detail"]["fresh_training_reductions_executed"] = 1
        self.rebind()
        self.assert_rejected()

    def test_projection_requires_the_specific_independent_review_contract(self):
        for change in ("relation", "review", "scope"):
            with self.subTest(change=change):
                self.make_fixture()
                if change == "relation":
                    self.docs["relation"]["relation"] = "SYNTHETIC UNKNOWN PROJECTION"
                elif change == "review":
                    self.docs["review"]["verdict"] = "SYNTHETIC UNREVIEWED"
                else:
                    self.docs["scope"]["schema"] = "SYNTHETIC UNKNOWN SCOPE"
                self.rebind()
                self.assert_rejected()

    def test_training_must_be_exact_two_of_three_prefix(self):
        for change in ("source_order", "new_target", "nonprefix", "missing_source"):
            with self.subTest(change=change):
                self.make_fixture()
                if change == "source_order":
                    self.docs["qualification"]["training"].reverse()
                elif change == "new_target":
                    for name in ("target", "producer_workload"):
                        self.docs[name]["training"][0]["powers"] = self.key(9999)
                elif change == "nonprefix":
                    for name in ("target", "producer_workload"):
                        self.docs[name]["training"] = copy.deepcopy(
                            self.docs["qualification"]["training"][1:]
                        )
                else:
                    self.docs["qualification"]["training"].pop()
                self.rebind()
                self.assert_rejected()

    def test_public_case_identity_role_and_order_are_preserved(self):
        for field in ("id", "powers", "public_input", "rationale", "order", "missing"):
            with self.subTest(field=field):
                self.make_fixture()
                roots = self.docs["qualification"]["held_out"]
                if field == "order":
                    roots.reverse()
                elif field == "missing":
                    roots.pop()
                elif field == "powers":
                    roots[0][field] = self.key(999)
                else:
                    roots[0][field] += " ALTERED"
                self.rebind()
                self.assert_rejected()

    def test_proof_and_review_details_cannot_omit_or_reorder_cases(self):
        for change in (
            "missing_proof",
            "duplicate_proof",
            "reordered_details",
            "missing_details",
        ):
            with self.subTest(change=change):
                self.make_fixture()
                if change == "missing_proof":
                    self.docs["proof"]["rows"].pop()
                elif change == "duplicate_proof":
                    self.docs["proof"]["rows"].append(
                        copy.deepcopy(self.docs["proof"]["rows"][0])
                    )
                elif change == "reordered_details":
                    self.docs["detail"][
                        "full_case_key_value_and_guard_evidence"
                    ].reverse()
                else:
                    self.docs["detail"]["full_case_key_value_and_guard_evidence"].pop()
                self.rebind()
                self.assert_rejected()

    def test_partial_failed_or_backed_source_worker_is_rejected(self):
        for change in (
            "success",
            "worker_success",
            "exit",
            "deadline",
            "backing",
            "variant",
        ):
            with self.subTest(change=change):
                self.make_fixture()
                worker = self.docs["worker4"]
                if change in ("success", "worker_success"):
                    worker[change] = False
                elif change == "exit":
                    worker["supervision"]["exit_code"] = 124
                elif change == "deadline":
                    worker["supervision"]["deadline_exhausted"] = True
                elif change == "backing":
                    worker["persistent_cache_attachments_during_evaluation"]["5"] = True
                else:
                    worker["variant"] = "cached"
                self.rebind()
                self.assert_rejected()

    def test_source_worker_cannot_change_case_or_cohort_bindings(self):
        for field in (
            "case_id",
            "declared_common_powers",
            "build_git_revision",
            "data_manifest_sha256",
            "executable_sha256",
        ):
            with self.subTest(field=field):
                self.make_fixture()
                if field == "case_id":
                    self.docs["worker4"][field] = "SYNTHETIC_OTHER"
                elif field == "declared_common_powers":
                    self.docs["worker4"][field] = self.key(9999)
                else:
                    self.docs["worker4"]["identity"][field] = "e" * 64
                self.rebind()
                self.assert_rejected()

    def test_reviewed_worker_source_bindings_are_required(self):
        for change in (
            "missing_review_source",
            "wrong_compiled_source",
            "wrong_snapshot",
        ):
            with self.subTest(change=change):
                self.make_fixture()
                if change == "missing_review_source":
                    self.docs["review"]["source_files"].pop()
                elif change == "wrong_compiled_source":
                    self.docs["worker2"]["identity"]["compiled_driver_sha256"] = (
                        "e" * 64
                    )
                else:
                    self.docs["snapshot"]["source_sha256"][
                        "crates/vakint/examples/measure_cache.rs"
                    ] = "e" * 64
                self.rebind()
                self.assert_rejected()

    def test_rehashed_changed_full_result_or_guards_fail_historical_parity(self):
        for filename, field in (
            ("value.symbolica", "value_sha256"),
            ("guards.json", "guards_sha256"),
        ):
            with self.subTest(filename=filename):
                self.make_fixture()
                path = self.paths["worker3"].with_name(filename)
                path.write_text("SYNTHETIC CHANGED COMPLETE PAYLOAD")
                self.docs["worker3"][field] = merge.digest(path)
                self.rebind()
                self.assert_rejected()

    def test_complete_source_sidecars_are_rehashed(self):
        self.paths["worker3"].with_name("guards.json").unlink()
        with self.assertRaises(FileNotFoundError):
            self.evaluate()

    def test_authenticated_raw_owner_sets_and_union_are_exact(self):
        for change in ("raw", "owner", "union", "mapping_error"):
            with self.subTest(change=change):
                self.make_fixture()
                if change == "raw":
                    self.docs["proof"]["rows"][3]["actual_lowered_roots"] = [
                        self.key(9999)
                    ]
                elif change == "owner":
                    self.docs["proof"]["rows"][3]["protected_owner_roots"] = [
                        self.key(9999)
                    ]
                elif change == "union":
                    self.docs["union"]["forbidden_direct_target_or_owner_keys"].pop()
                else:
                    self.docs["proof"]["rows"][3]["mapping_error"] = "SYNTHETIC FAILURE"
                self.rebind()
                self.assert_rejected()

    def test_selected_training_cannot_alias_a_protected_raw_key(self):
        protected = self.docs["target"]["tuning"][0]["powers"]
        for name in ("qualification", "target", "producer_workload"):
            self.docs[name]["training"][0]["powers"] = protected
        self.docs["routes"]["rows"][0]["raw_powers"] = protected
        self.docs["aliases"]["rows"][0]["powers"] = protected
        self.rebind()
        self.assert_rejected()

    def test_twelve_held_out_owner_sets_must_remain_distinct(self):
        first = self.docs["target"]["held_out"][0]["protected_roots"]
        self.docs["target"]["held_out"][1]["protected_roots"] = first
        self.docs["producer_workload"]["held_out"][1]["protected_roots"] = first
        self.docs["proof"]["rows"][2]["protected_owner_roots"] = first
        self.rebind()
        self.assert_rejected()

    def test_direct_routes_must_bind_the_source_census(self):
        self.docs["routes"]["source_census_sha256"] = "e" * 64
        self.rebind()
        self.assert_rejected()

    def test_selected_training_cannot_alias_any_nonzero_direct_owner(self):
        self.docs["routes"]["rows"][0]["owner_endpoints"][0]["powers"] = self.docs[
            "target"
        ]["held_out"][5]["protected_roots"][0]
        self.rebind()
        self.assert_rejected()

    def test_direct_routes_require_nonzero_coefficients_and_eligibility(self):
        for change in ("zero", "empty", "route_error", "alias"):
            with self.subTest(change=change):
                self.make_fixture()
                if change == "zero":
                    self.docs["routes"]["rows"][0]["owner_endpoints"][0][
                        "coefficient"
                    ] = "0"
                elif change == "empty":
                    self.docs["routes"]["rows"][0]["owner_endpoints"] = []
                elif change == "route_error":
                    self.docs["routes"]["rows"][0]["route_error"] = "SYNTHETIC FAILURE"
                else:
                    self.docs["aliases"]["rows"][0][
                        "eligible_after_external_alias_exclusion"
                    ] = False
                self.rebind()
                self.assert_rejected()

    def test_current_target_native_api_validation_remains_required(self):
        for field, value in (
            ("valid", False),
            ("public_cases_checked", 12),
            ("data_manifest_sha256", "e" * 64),
            ("reduction_or_cache_import_performed", True),
        ):
            with self.subTest(field=field):
                self.make_fixture()
                self.docs["validator"][field] = value
                self.rebind()
                self.assert_rejected()

    def test_missing_current_native_api_or_import_reference_is_rejected(self):
        for field in (
            "actual_existing_validator_two_root_result",
            "strict_import_equality_validation",
        ):
            with self.subTest(field=field):
                self.make_fixture()
                self.docs["detail"][field] = None
                self.save("detail")
                # Rehash the containing authority chain without repairing the missing API ref.
                self.docs["packet"]["exact_semantic_projection_evidence"] = (
                    self.reference("detail")
                )
                self.save("packet")
                self.docs["review"]["evidence_packet"] = self.reference("packet")
                self.save("review")
                self.docs["scope"]["review"] = self.reference("review")
                self.save("scope")
                self.docs["relation"].update(
                    evidence=self.reference("detail"),
                    independent_review=self.reference("review"),
                    validator_scope=self.reference("scope"),
                )
                self.save("relation")
                self.accepted["qualification_projection"] = self.reference("relation")
                self.base.refresh_acceptance(
                    self.accepted_path, self.index, self.accepted
                )
                self.assert_rejected()

    def test_cached_import_requires_actual_success_source_and_backing(self):
        for change in (
            "success",
            "auth",
            "cohort",
            "data",
            "census",
            "cache",
            "backing",
            "import_check",
            "exit",
        ):
            with self.subTest(change=change):
                self.make_fixture()
                if change == "success":
                    self.docs["cached"]["success"] = False
                elif change == "auth":
                    self.docs["cached"]["native_keys_authenticated"] = False
                elif change == "cohort":
                    self.docs["cached"]["identity"]["build_git_revision"] = "e" * 40
                elif change == "data":
                    self.docs["cached"]["identity"]["data_manifest_sha256"] = "e" * 64
                elif change == "census":
                    self.docs["cached"]["identity"]["source_census_sha256"] = "e" * 64
                elif change == "cache":
                    self.docs["cached"]["identity"]["cache_sha256"] = "e" * 64
                elif change == "backing":
                    self.docs["cached"][
                        "persistent_cache_attachments_during_evaluation"
                    ]["4"] = True
                elif change == "import_check":
                    self.docs["import"]["checks"][
                        "native_binding_and_every_record_validation"
                    ] = False
                else:
                    self.docs["import_supervision"]["exit_code"] = 124
                self.rebind()
                self.assert_rejected()

    def test_original_producer_identity_cannot_be_relabelled(self):
        self.producer["raw_report"]["build_identity"]["build"]["git_revision"] = (
            "e" * 40
        )
        self.rebind()
        self.assert_rejected()

    def test_failed_target_attempt_is_visible_and_never_promoted(self):
        for change in ("complete", "accepted", "counts", "worker", "mapper"):
            with self.subTest(change=change):
                self.make_fixture()
                if change in ("complete", "accepted"):
                    self.docs["failed"][change] = True
                elif change == "counts":
                    self.docs["failed"]["counts"]["functional_success"] = 13
                elif change == "worker":
                    self.docs["failed"]["failed_case"]["worker_json_exists"] = True
                else:
                    self.docs["failed"]["mapper_launched"] = True
                self.rebind()
                self.assert_rejected()

    def test_projection_does_not_transfer_source_timings_or_create_a_window(self):
        result = self.evaluate()
        self.assertNotIn("evaluation_seconds", result["qualification_projection"])
        self.assertNotIn("metrics", result["qualification_projection"])
        functional = {
            "case_id": "SYNTHETIC_PUBLIC_0",
            "launcher": {"final_performance_eligible": False},
            "artifact": self.docs["artifact"],
        }
        with self.assertRaisesRegex(ValueError, "functional launcher"):
            merge.ReportEvidence(self.index).qualify(
                functional,
                merge.digest(self.base.binary),
                "SYNTHETIC_PUBLIC_0",
                "uncached",
                0,
            )


if __name__ == "__main__":
    unittest.main()
