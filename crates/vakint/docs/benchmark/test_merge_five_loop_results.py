"""Synthetic file-contract tests only; no reduction or performance measurement."""

import copy
import hashlib
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
