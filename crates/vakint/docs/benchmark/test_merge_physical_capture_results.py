"""Synthetic file contracts only: no native reduction, cache import or benchmarks."""

import copy
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import merge_physical_capture_results as adapter


class PhysicalCaptureContracts(unittest.TestCase):
    def setUp(self):
        # TemporaryDirectory honours the caller's development TMPDIR.
        self.temp = tempfile.TemporaryDirectory(prefix="synthetic-physical-capture-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.binary = self.file("collector", b"\x7fELFsynthetic-not-executable")
        self.data = self.file("data.json", {"version": "v-test"})
        self.input = self.file("prepared/vakint-input.symbolica", b"fixture-input")
        self.guards = self.file("prepared/projector-conditions.symbolica", b"")
        self.manifest_data = {
            "schema": "vakint.qcd-ladder.prepared-input.v1",
            "projection": "zero-momentum-metric",
            "input": self.input,
            "nonzero_conditions_file": self.guards,
            "contributions": 1,
        }
        self.manifest = self.file("prepared/prepared-input.json", self.manifest_data)
        self.acceptance = self.file("science.json", {"synthetic": True})
        self.source_revision = "c" * 40
        self.native_tree = "b" * 40
        self.source_freeze = self.file(
            "source-freeze.json",
            {"head": self.source_revision, "native_tree": self.native_tree},
        )
        self.artifact = self.file(
            "artifact.json",
            {
                "state": "PASS_OPTIMIZED_CAPTURE_PRODUCER_BUILD_AND_TINY_GATE",
                "elf": {
                    "sha256": self.binary["sha256"],
                    "native_tree": self.native_tree,
                    "source_freeze_sha256": self.source_freeze["sha256"],
                    "profile": {
                        "opt_level": "3",
                        "test": False,
                        "debug_assertions": False,
                    },
                },
                "source_head": self.source_revision,
                "source_freeze": self.source_freeze,
                "actual_rustc_codegen": {"codegen_flags": ["opt-level=3", "lto=fat"]},
            },
        )
        prepared = {
            "manifest": self.manifest,
            "manifest_projection": "zero-momentum-metric",
            "contributions": 1,
            "input": self.input,
            "guards": self.guards,
            "scientific_acceptance": self.acceptance,
            "provenance_addendum": self.acceptance,
        }
        cohort = {
            "executable": self.binary,
            "artifact_receipt": self.artifact,
            "vakint_revision": self.source_revision,
            "data_manifest": self.data,
            "data_identity": "data-source",
            "data_version": "v-test",
            "dependency_receipts": [],
        }
        self.producer = {
            "schema": "vakint.inactive-capture-physical-user-unit.v1",
            "job_id": "q0",
            "projection": "zero-momentum",
            "performance_eligible": False,
            "output_root": str(self.root / "producer"),
            "prepared_input": prepared,
            "validated_cohort": cohort,
            "running_cache_directory": str(self.root / "cache"),
            "capture_only_producer": True,
            "required_predecessor": None,
            "running_cache_arguments": ["--running-cache-capture-only"],
            "cpus": [1],
        }
        self.producer_config = self.file("q0.final.json", self.producer)
        self.producer_done = self.run_fixture(
            self.producer, self.producer_config, 11, 20.0, True
        )
        self.reuse = copy.deepcopy(self.producer)
        self.reuse.update(
            output_root=str(self.root / "reuse"),
            capture_only_producer=False,
            running_cache_arguments=[],
            required_predecessor=self.producer_done,
            seen_reuse_contract={
                "fresh_process": True,
                "same_elf_input_settings_data_affinity": True,
                "producer_config_sha256": self.producer_config["sha256"],
                "producer_completion_sha256": self.producer_done["sha256"],
            },
        )
        self.reuse_config = self.file("q0.seen-reuse.final.json", self.reuse)
        self.reuse_done = self.run_fixture(
            self.reuse, self.reuse_config, 12, 2.0, False
        )
        self.comparison = {
            "state": "PASS_ONE_FRESH_SEEN_REUSE_EXACT_PARITY",
            "family_map_equal": True,
            "seen_not_heldout": True,
            "fresh_process": True,
            "native_snapshot_unchanged": True,
            "newly_captured_records": 0,
            "pure_ibp_seconds": None,
            "producer_completion_sha256": self.producer_done["sha256"],
            "reuse_completion_sha256": self.reuse_done["sha256"],
            "exact_complete_payload_sha256_equal": dict.fromkeys(
                [
                    "result.symbolica",
                    "nonzero-conditions.symbolica",
                    "unknown-masters.json",
                ],
                True,
            ),
            "cached_persistent_statistics": {
                "lookups": 4,
                "route_hits": 1,
                "apply_hits": 3,
            },
            "producer_reduction_seconds": 20.0,
            "cached_reduction_seconds": 2.0,
            "producer_wall_seconds": 25.0,
            "reuse_wall_seconds": 7.0,
            "raw_reduction_speedup": 10.0,
            "cache_preserved": self.producer["running_cache_directory"],
        }
        self.file("seen-reuse-comparison.json", self.comparison)

    def file(self, name, data):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(data, bytes):
            path.write_bytes(data)
        else:
            path.write_text(json.dumps(data, indent=2) + "\n")
        return {
            "path": str(path),
            "sha256": adapter.digest(path),
            "bytes": path.stat().st_size,
        }

    def run_fixture(self, config, config_ref, pid, duration, producer):
        root = Path(config["output_root"])
        prepared = config["prepared_input"]
        prefix = root.relative_to(self.root) / "run"
        result = self.file(prefix / "result.symbolica", b"full-fixture-result")
        guards = self.file(prefix / "nonzero-conditions.symbolica", b"D-4\n")
        unknown = self.file(prefix / "unknown-masters.json", [])
        snapshot = {
            "schema": "rustred.running-application-cache-snapshot.v1",
            "directory": config["running_cache_directory"],
            "delta": {"generation": 3, "manifest_digest": "a" * 64},
        }
        index = self.file(prefix / "snapshot-index", b"synthetic-index-not-a-cache")
        snapshot["delta"].update(
            manifest_path=index["path"], manifest_bytes=index["bytes"]
        )
        descriptor = self.file(prefix / "snapshot.json", snapshot)
        native = {
            "pending_records": 0,
            "pending_bytes": 0,
            "accepted_records": 8 if producer else 0,
            "committed_records": 8 if producer else 0,
            "writer_error": None,
            "closed": True,
            "last_error": None,
            "capture_metadata_exhausted": False,
        }
        cache = {
            "report": native,
            "snapshot": snapshot,
            "close": {"error": None},
            "flush_complete": True,
        }
        cache_ref = self.file(prefix / "running-cache.json", cache)
        fields = {
            "flush_complete": True,
            "writer_fault": None,
            **{
                k: native[k]
                for k in (
                    "pending_records",
                    "pending_bytes",
                    "accepted_records",
                    "committed_records",
                )
            },
            **{
                key: snapshot["delta"][key] for key in ("generation", "manifest_digest")
            },
        }
        report = {
            "schema": "vakint.qcd-ladder.benchmark.v1",
            "retained_result_available": True,
            "reduction_completed": True,
            "source_commit": self.source_revision,
            "executable_sha256": self.binary["sha256"],
            "data_manifest_sha256": self.data["sha256"],
            "input_mode": "prepared",
            "coefficient_representation": "sparse",
            "fresh_native_preparation": False,
            "preparation_seconds": None,
            "projection": prepared["manifest_projection"],
            "contributions": prepared["contributions"],
            "prepared_input": {
                "manifest_sha256": prepared["manifest"]["sha256"],
                "manifest": adapter.read(prepared["manifest"]["path"]),
            },
            "input_sha256": prepared["input"]["sha256"],
            "input_bytes": prepared["input"]["bytes"],
            "result_sha256": result["sha256"],
            "result_bytes": result["bytes"],
            "guards": {
                "file": "nonzero-conditions.symbolica",
                "sha256": guards["sha256"],
                "count": 1,
                "includes_preparation_conditions": True,
            },
            "catalogue_closure": {
                "closed": True,
                "unknown_count": 0,
                "unknown_sha256": unknown["sha256"],
            },
            "persistent_cache": None,
            "running_cache_attached": True,
            "persistent_cache_attachments_before": {str(n): False for n in range(1, 6)},
            "persistent_cache_attachments_during_evaluation": {
                str(n): not producer and n == 5 for n in range(1, 6)
            },
            "data": {
                "load_prepared_caches": False,
                "version": "v-test",
                "loops": [
                    {"loop_count": n, "cache_loaded": False, "source": "local"}
                    for n in range(1, 6)
                ],
            },
            "reduction_seconds": duration,
            "conventions": {"fixture-settings": True},
            "families": {"family": "identity"},
            "running_cache_capture_only": producer,
            "running_cache_read_backing": not producer,
            "saved_application_cache": {
                "read_backing": not producer,
                "persistent": {
                    "lookups": 0 if producer else 4,
                    "route_hits": 0 if producer else 1,
                    "apply_hits": 0 if producer else 3,
                },
            },
        }
        report_ref = self.file(prefix / "benchmark.json", report)
        completion = {
            "job_id": config["job_id"],
            "performance_eligible": False,
            "functional_acceptance": True,
            "process_zero_acceptance": True,
            "collector_exit_code": 0,
            "config_sha256": config_ref["sha256"],
            "pure_ibp_seconds": None,
            "native_child_restarted_or_signalled": False,
            "collector_wall_seconds_through_actual_wait": duration + 5,
            "collector": {"pid": pid, "start_ticks": pid * 10},
            "post_executable": self.binary,
            "post_artifact_receipt": self.artifact,
            "post_selected_data": {
                "manifest": self.data,
                "revision": "data-source",
                "version": "v-test",
            },
            "terminal_artifacts": {
                Path(ref["path"]).name: {**ref, "exists": True}
                for ref in (report_ref, result, guards, unknown, cache_ref)
            },
            "outputs": {
                "benchmark_sha256": report_ref["sha256"],
                "result_sha256": result["sha256"],
                "result_bytes": result["bytes"],
                "guards_sha256": guards["sha256"],
                "guards_count": 1,
                "unknown_sha256": unknown["sha256"],
                "unknown_count": 0,
                "running_cache": {
                    "report": cache_ref,
                    "native_fields": fields,
                    "snapshot_descriptor": descriptor,
                    "snapshot_index": index,
                },
            },
        }
        self.file(
            root.relative_to(self.root) / "launch.json",
            {
                "job_id": config["job_id"],
                "config_sha256": config_ref["sha256"],
                "collector": completion["collector"],
                "argv": [
                    self.binary["path"],
                    "--projection",
                    config["projection"],
                ],
            },
        )
        return self.file(root.relative_to(self.root) / "completion.json", completion)

    def quadratic_fixture(self):
        config = copy.deepcopy(self.producer)
        config.update(
            job_id="quadratic",
            projection="quadratic",
            output_root=str(self.root / "quadratic"),
        )
        prepared = config["prepared_input"]
        prepared.update(
            manifest_projection="quadratic-transverse",
            contributions=21,
            input=self.file(
                "quadratic-prepared/vakint-input.symbolica", b"full-quadratic-fixture"
            ),
            guards=self.file(
                "quadratic-prepared/projector-conditions.symbolica",
                b"d-1\nd-2\nd-3\nd-4\nd-5\n",
            ),
        )
        prepared["manifest"] = self.file(
            "quadratic-prepared/prepared-input.json",
            {
                **self.manifest_data,
                "projection": prepared["manifest_projection"],
                "contributions": prepared["contributions"],
                "input": prepared["input"],
                "nonzero_conditions_file": prepared["guards"],
            },
        )
        return config

    def refresh_bindings(self):
        """Rehash synthetic receipts so semantic failures cannot hide behind stale hashes."""
        for mode, filename in (
            ("producer", "q0.final.json"),
            ("reuse", "q0.seen-reuse.final.json"),
        ):
            config = adapter.read(self.root / filename)
            if mode == "reuse":
                config["required_predecessor"] = self.producer_done
                config["seen_reuse_contract"].update(
                    producer_config_sha256=adapter.digest(self.root / "q0.final.json"),
                    producer_completion_sha256=self.producer_done["sha256"],
                )
            config_ref = self.file(filename, config)
            prefix = Path(mode) / "run"
            completion = adapter.read(self.root / mode / "completion.json")
            completion["config_sha256"] = config_ref["sha256"]
            launch_path = Path(mode) / "launch.json"
            launch = adapter.read(self.root / launch_path)
            launch["config_sha256"] = config_ref["sha256"]
            launch["collector"] = completion["collector"]
            self.file(launch_path, launch)
            report = adapter.read(self.root / prefix / "benchmark.json")
            cache = adapter.read(self.root / prefix / "running-cache.json")
            bound = completion["outputs"]["running_cache"]
            bound["snapshot_descriptor"] = self.file(
                prefix / "snapshot.json", cache["snapshot"]
            )
            for name, item in completion["terminal_artifacts"].items():
                item.update(
                    sha256=adapter.digest(self.root / prefix / name),
                    bytes=(self.root / prefix / name).stat().st_size,
                )
            completion["outputs"]["benchmark_sha256"] = adapter.digest(
                self.root / prefix / "benchmark.json"
            )
            bound["report"] = {
                **bound["report"],
                "sha256": adapter.digest(self.root / prefix / "running-cache.json"),
                "bytes": (self.root / prefix / "running-cache.json").stat().st_size,
            }
            bound["native_fields"] = {
                "flush_complete": cache["flush_complete"],
                "writer_fault": cache["report"]["writer_error"],
                **{
                    key: cache["report"][key]
                    for key in (
                        "pending_records",
                        "pending_bytes",
                        "accepted_records",
                        "committed_records",
                    )
                },
                **{
                    key: cache["snapshot"]["delta"][key]
                    for key in ("generation", "manifest_digest")
                },
            }
            completed = self.file(Path(mode) / "completion.json", completion)
            if mode == "producer":
                self.producer_done = completed
            else:
                self.reuse_done = completed
            self.comparison[mode + "_completion_sha256"] = completed["sha256"]
            self.comparison[
                "producer_reduction_seconds"
                if mode == "producer"
                else "cached_reduction_seconds"
            ] = report["reduction_seconds"]
        self.file("seen-reuse-comparison.json", self.comparison)

    def test_accepts_one_functional_observation_without_performance_aggregates(self):
        run = adapter.physical_capture(self.root)
        self.assertEqual(run["status"], "functional-seen-physical-reuse")
        self.assertFalse(run["performance_eligible"])
        self.assertTrue(run["seen_not_heldout"])
        self.assertEqual(run["sample_count"], 1)
        self.assertEqual(run["observed_seconds"]["cached_reduction"], 2.0)
        self.assertEqual(run["raw_comparison"], self.comparison)
        self.assertFalse(
            any("median" in key or "ratio" in key or "speedup" in key for key in run)
        )
        self.assertIsNone(run["pure_ibp_seconds"])

    def test_completed_run_accepts_quadratic_only_with_explicit_job(self):
        config = self.quadratic_fixture()
        config_ref = self.file("quadratic.final.json", config)
        self.run_fixture(config, config_ref, 13, 30.0, True)
        run = adapter.completed_run(config_ref["path"], expected_job_id="quadratic")
        self.assertEqual(run["config"]["job_id"], "quadratic")
        self.assertEqual(run["report"]["projection"], "quadratic-transverse")
        self.assertEqual(run["report"]["contributions"], 21)
        self.assertFalse(run["completion"]["performance_eligible"])
        with self.assertRaisesRegex(ValueError, "functional completion"):
            adapter.completed_run(config_ref["path"])

    def test_completed_run_rejects_unknown_expected_jobs(self):
        for job in ("", "other", "Q0", "zero-momentum", None):
            with (
                self.subTest(job=job),
                self.assertRaisesRegex(ValueError, "unknown physical job"),
            ):
                adapter.completed_run(self.producer_config["path"], expected_job_id=job)

    def test_completed_run_rejects_job_and_projection_mismatches(self):
        for job in ("q0", "quadratic"):
            original = self.producer if job == "q0" else self.quadratic_fixture()
            other = "quadratic" if job == "q0" else "q0"
            for field, value in (
                ("job_id", other),
                ("job_id", "unrecognized"),
                ("projection", "quadratic" if job == "q0" else "zero-momentum"),
                ("projection", "unrecognized"),
                ("prepared_input.manifest_projection", "quadratic-transverse"),
                ("prepared_input.manifest_projection", "zero-momentum-metric"),
                ("prepared_input.manifest_projection", "unrecognized"),
            ):
                config = copy.deepcopy(original)
                target = config
                keys = field.split(".")
                for key in keys[:-1]:
                    target = target[key]
                if target[keys[-1]] == value:
                    continue
                target[keys[-1]] = value
                with self.subTest(job=job, field=field, value=value):
                    config_ref = self.file("candidate.final.json", config)
                    self.run_fixture(config, config_ref, 13, 30.0, True)
                    with self.assertRaisesRegex(
                        ValueError, "functional completion|job/projection identity"
                    ):
                        adapter.completed_run(config_ref["path"], expected_job_id=job)
            config_ref = self.file("candidate.final.json", original)
            completion_ref = self.run_fixture(original, config_ref, 13, 30.0, True)
            completion = adapter.read(completion_ref["path"])
            completion["job_id"] = other
            self.file(completion_ref["path"], completion)
            with (
                self.subTest(job=job, field="completion.job_id"),
                self.assertRaisesRegex(ValueError, "functional completion"),
            ):
                adapter.completed_run(config_ref["path"], expected_job_id=job)

    def test_quadratic_rejects_coherently_substituted_zero_momentum_fixture(self):
        config = copy.deepcopy(self.producer)
        config["job_id"] = "quadratic"
        config_ref = self.file("quadratic.final.json", config)
        self.run_fixture(config, config_ref, 13, 30.0, True)
        with self.assertRaisesRegex(ValueError, "job/projection identity"):
            adapter.completed_run(config_ref["path"], expected_job_id="quadratic")

    def test_completed_run_rejects_missing_and_mismatched_launches(self):
        for job in ("q0", "quadratic"):
            config = self.producer if job == "q0" else self.quadratic_fixture()
            config_ref = self.file("candidate.final.json", config)
            self.run_fixture(config, config_ref, 13, 30.0, True)
            launch_path = Path(config["output_root"]) / "launch.json"
            original = adapter.read(launch_path)
            launch_path.unlink()
            with (
                self.subTest(job=job, field="missing launch"),
                self.assertRaises(FileNotFoundError),
            ):
                adapter.completed_run(config_ref["path"], expected_job_id=job)
            binary = self.binary["path"]
            projection = config["projection"]
            for field, value in (
                ("job_id", "quadratic" if job == "q0" else "q0"),
                ("config_sha256", "0" * 64),
                ("collector.pid", 14),
                ("collector.start_ticks", 131),
                ("argv", []),
                ("argv", None),
                ("argv", [binary]),
                ("argv", [binary, "--projection"]),
                ("argv", [binary, "--projection", "unrecognized"]),
                (
                    "argv",
                    [
                        binary,
                        "--projection",
                        "quadratic" if job == "q0" else "zero-momentum",
                    ],
                ),
                (
                    "argv",
                    [binary, "--projection", projection, "--projection", projection],
                ),
                (
                    "argv",
                    [binary, "--projection", projection, "--projection=" + projection],
                ),
                ("argv", [binary, "--projection=" + projection]),
                ("argv", [binary, "--projection", False]),
                (
                    "argv",
                    [str(self.root / "other-collector"), "--projection", projection],
                ),
            ):
                with self.subTest(job=job, field=field, value=value):
                    launch = copy.deepcopy(original)
                    target = launch
                    keys = field.split(".")
                    for key in keys[:-1]:
                        target = target[key]
                    target[keys[-1]] = value
                    self.file(launch_path, launch)
                    with self.assertRaisesRegex(ValueError, "physical launch"):
                        adapter.completed_run(config_ref["path"], expected_job_id=job)

    def test_quadratic_keeps_prepared_manifest_and_report_projection_checks(self):
        config = self.quadratic_fixture()
        manifest = adapter.read(config["prepared_input"]["manifest"]["path"])
        manifest["projection"] = "zero-momentum-metric"
        config["prepared_input"]["manifest"] = self.file(
            config["prepared_input"]["manifest"]["path"], manifest
        )
        config_ref = self.file("quadratic.final.json", config)
        self.run_fixture(config, config_ref, 13, 30.0, True)
        with self.assertRaisesRegex(ValueError, "prepared manifest differs"):
            adapter.completed_run(config_ref["path"], expected_job_id="quadratic")

        config = self.quadratic_fixture()
        config_ref = self.file("quadratic.final.json", config)
        completion_ref = self.run_fixture(config, config_ref, 13, 30.0, True)
        root = Path(config["output_root"])
        report_path = root / "run/benchmark.json"
        report = adapter.read(report_path)
        report["projection"] = "zero-momentum-metric"
        report_ref = self.file(report_path, report)
        completion = adapter.read(completion_ref["path"])
        completion["outputs"]["benchmark_sha256"] = report_ref["sha256"]
        completion["terminal_artifacts"]["benchmark.json"].update(report_ref)
        self.file(completion_ref["path"], completion)
        with self.assertRaisesRegex(ValueError, "prepared-input identity differs"):
            adapter.completed_run(config_ref["path"], expected_job_id="quadratic")

    def test_rejects_missing_completion(self):
        Path(self.reuse_done["path"]).unlink()
        with self.assertRaises(FileNotFoundError):
            adapter.physical_capture(self.root)

    def test_rejects_corrupted_complete_payload_and_executable(self):
        for path in (
            self.root / "reuse/run/result.symbolica",
            Path(self.binary["path"]),
        ):
            with self.subTest(path=path):
                original = path.read_bytes()
                path.write_bytes(b"corrupt")
                with self.assertRaises(ValueError):
                    adapter.physical_capture(self.root)
                path.write_bytes(original)

    def test_rejects_changed_comparison_claims_and_seconds(self):
        mutations = {
            "seen_not_heldout": False,
            "fresh_process": False,
            "native_snapshot_unchanged": False,
            "newly_captured_records": 1,
            "pure_ibp_seconds": 1,
            "cached_reduction_seconds": 3,
            "producer_completion_sha256": "0" * 64,
            "cached_persistent_statistics": {"lookups": 0},
            "family_map_equal": False,
        }
        for key, value in mutations.items():
            with self.subTest(key=key):
                self.file("seen-reuse-comparison.json", {**self.comparison, key: value})
                with self.assertRaises(ValueError):
                    adapter.physical_capture(self.root)
        self.file("seen-reuse-comparison.json", self.comparison)

    def test_rejects_configuration_drift(self):
        self.file("q0.seen-reuse.final.json", {**self.reuse, "cpus": [2]})
        with self.assertRaises(ValueError):
            adapter.physical_capture(self.root)

    def test_rejects_all_misses_and_invalid_hit_counts_with_coherent_receipts(self):
        path = Path("reuse/run/benchmark.json")
        original = adapter.read(self.root / path)
        for route_hits, apply_hits in (
            (0, 0),
            (False, 1),
            (1, True),
            (1.0, 0),
            (0, -1),
        ):
            with self.subTest(route_hits=route_hits, apply_hits=apply_hits):
                report = copy.deepcopy(original)
                persistent = report["saved_application_cache"]["persistent"]
                persistent.update(route_hits=route_hits, apply_hits=apply_hits)
                self.assertGreater(persistent["lookups"], 0)
                self.comparison["cached_persistent_statistics"] = persistent.copy()
                self.file(path, report)
                self.refresh_bindings()
                with self.assertRaisesRegex(
                    ValueError, "successful persistent cache hits"
                ):
                    adapter.physical_capture(self.root)

    def test_rejects_missing_malformed_and_unbound_provenance(self):
        original = adapter.read(Path(self.artifact["path"]))
        mutations = (
            ("elf.native_tree", None),
            ("elf.native_tree", "a" * 39),
            ("elf.native_tree", "z" * 40),
            ("elf.native_tree", "a" * 40),
            ("elf.source_freeze_sha256", None),
            ("elf.source_freeze_sha256", "a" * 63),
            ("elf.source_freeze_sha256", 42),
            ("elf.source_freeze_sha256", "a" * 64),
            ("source_head", "short"),
            ("source_freeze", None),
        )
        for field, value in mutations:
            with self.subTest(field=field, value=value):
                artifact = copy.deepcopy(original)
                target = artifact
                parts = field.split(".")
                for key in parts[:-1]:
                    target = target[key]
                if value is None:
                    del target[parts[-1]]
                else:
                    target[parts[-1]] = value
                artifact_ref = self.file("artifact.json", artifact)
                for mode, config_name in (
                    ("producer", "q0.final.json"),
                    ("reuse", "q0.seen-reuse.final.json"),
                ):
                    config = adapter.read(self.root / config_name)
                    config["validated_cohort"]["artifact_receipt"] = artifact_ref
                    self.file(config_name, config)
                    completion_path = Path(mode) / "completion.json"
                    completion = adapter.read(self.root / completion_path)
                    completion["post_artifact_receipt"] = artifact_ref
                    self.file(completion_path, completion)
                self.refresh_bindings()
                with self.assertRaisesRegex(ValueError, "provenance"):
                    adapter.physical_capture(self.root)

    def test_rejects_same_process_even_with_different_observer_metadata(self):
        path = Path("reuse/completion.json")
        completion = adapter.read(self.root / path)
        completion["collector"] = {
            "pid": 11,
            "start_ticks": 110,
            "observer_extra": True,
        }
        self.file(path, completion)
        self.refresh_bindings()
        with self.assertRaisesRegex(ValueError, "fresh process"):
            adapter.physical_capture(self.root)

    def test_rejects_bad_durable_close_with_coherent_receipt_hashes(self):
        path = Path("reuse/run/running-cache.json")
        original = adapter.read(self.root / path)
        for field, value in (
            ("flush_complete", False),
            ("close.error", "failed"),
            ("report.pending_records", 1),
            ("report.pending_bytes", 1),
            ("report.writer_error", "failed"),
            ("report.closed", False),
        ):
            with self.subTest(field=field):
                cache = copy.deepcopy(original)
                target = cache
                parts = field.split(".")
                for key in parts[:-1]:
                    target = target[key]
                target[parts[-1]] = value
                self.file(path, cache)
                self.refresh_bindings()
                with self.assertRaisesRegex(ValueError, "durable checkpoint"):
                    adapter.physical_capture(self.root)
        self.file(path, original)
        self.refresh_bindings()

    def test_rejects_changed_snapshot_with_coherent_receipt_hashes(self):
        path = Path("reuse/run/running-cache.json")
        cache = adapter.read(self.root / path)
        cache["snapshot"]["delta"]["generation"] += 1
        self.file(path, cache)
        self.refresh_bindings()
        with self.assertRaisesRegex(ValueError, "changed cache snapshot"):
            adapter.physical_capture(self.root)

    def test_rejects_performance_upgrade_with_coherent_receipt_hashes(self):
        path = Path("reuse/completion.json")
        completion = adapter.read(self.root / path)
        completion["performance_eligible"] = True
        self.file(path, completion)
        self.refresh_bindings()
        with self.assertRaisesRegex(ValueError, "functional completion"):
            adapter.physical_capture(self.root)

    def test_append_preserves_history_and_selection_and_refuses_overwrite(self):
        base = {
            "selected_reduction_run": "four-old",
            "selected_five_loop_run": "five-old",
            "five_loop_runs": {
                "five-old": {"primary": "untouched", "held_out": "untouched"}
            },
            "physical_capture_runs": {"previous": {"raw": True}},
            "selected_physical_capture_run": "previous",
        }
        source = self.file("base.json", base)
        command = [
            sys.executable,
            str(Path(adapter.__file__)),
            "--base",
            source["path"],
            "--campaign",
            str(self.root),
            "--output",
            str(self.root / "new.json"),
        ]
        subprocess.run(command, check=True, capture_output=True)
        output = adapter.read(self.root / "new.json")
        self.assertEqual(output["selected_physical_capture_run"], "previous")
        for key in base.keys() - {"physical_capture_runs"}:
            self.assertEqual(output[key], base[key])
        self.assertEqual(
            output["physical_capture_runs"]["previous"],
            base["physical_capture_runs"]["previous"],
        )
        self.assertEqual(adapter.read(source["path"]), base)
        self.assertNotEqual(
            subprocess.run(command, check=False, capture_output=True).returncode, 0
        )
        command[-1] = str(self.root / "selected.json")
        subprocess.run(command + ["--select-current"], check=True, capture_output=True)
        selected = adapter.read(self.root / "selected.json")
        self.assertNotEqual(selected["selected_physical_capture_run"], "previous")
        self.assertEqual(selected["selected_five_loop_run"], "five-old")


if __name__ == "__main__":
    unittest.main()
