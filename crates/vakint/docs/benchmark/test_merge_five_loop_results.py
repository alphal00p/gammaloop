"""Synthetic file-contract tests only; no reduction or performance measurement."""

import copy
import json
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
            "fixture_sources_sha256": {"TEST_SOURCE": "6" * 64},
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


if __name__ == "__main__":
    unittest.main()
