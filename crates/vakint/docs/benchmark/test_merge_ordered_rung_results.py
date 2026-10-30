"""Synthetic file-contract tests; these do not run or time a reduction."""

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import merge_ordered_rung_results as merge


class OrderedReportContractTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="synthetic-ordered-report-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.base = {
            "historical": {"samples": list(range(48)), "scope": "unchanged"},
            "selected_reduction_run": "original",
            "reduction_runs": {"original": {"raw": [1.25, 2.5]}},
        }
        binary = self.root / "SYNTHETIC-NOT-AN-EXECUTABLE"
        binary.write_text("SYNTHETIC")
        self.producer = {
            "commit": "1" * 40,
            "cargo_lock_sha256": "2" * 64,
            "debug_assertions": False,
            "build_profile": "release",
            "executable_sha256": merge.digest(binary),
            "allowed_cpus": "160",
            "comparison_process": "one engine in this process",
        }
        self.parity = self.pair("sample-1")

    def write(self, path, document):
        path = self.root / path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(document))
        return path

    def update_report(self, name, change):
        parity = merge.read(self.parity)
        path = Path(parity["reports"][name]["path"])
        report = merge.read(path)
        change(report)
        path.write_text(json.dumps(report))
        parity["reports"][name] = merge.evidence(path)
        self.parity.write_text(json.dumps(parity))

    def add_numerical(self):
        reference = {
            "seconds": 1,
            "point": {"mUV": 1, "mu_r_squared": 1, "UFO_G": 1},
            "epsilon_range": "SYNTHETIC poles through finite",
            "selection": {
                "requested_decimal_precision": 300,
                "surviving_dependencies": [],
            },
            "coupling_definitions": {"G": "1"},
            "coefficients": [{"epsilon_power": -1}, {"epsilon_power": 0}],
        }
        for name in merge.read(self.parity)["reports"]:
            self.update_report(
                name, lambda row: row.update(numerical_reference=reference)
            )
        parity = merge.read(self.parity)
        parity["numerical_enclosures_agree"] = True
        self.parity.write_text(json.dumps(parity))

    def pair(self, name):
        reports, retained = {}, {}
        for engine, strategy in (("dfs", "DepthFirst"), ("ordered", "Ordered")):
            directory = self.root / name / engine
            directory.mkdir(parents=True)
            payload = directory / "retained.symbolica"
            payload.write_text(f"SYNTHETIC NON-PHYSICS {engine}")
            retained[engine] = merge.digest(payload)
            report = {
                "schema": "vakint.ordered-rung-comparison.v1",
                "loop_count": 2,
                "gluon_rungs": 1,
                "projection": "zero-momentum-metric",
                "strategy": strategy,
                "persistent_cache": False,
                "monitoring_enabled": False,
                "dfs_request_memoization": "SYNTHETIC normal per-run memoization",
                "cache_policy": {
                    "prepared_cache_loaded": False,
                    "previous_run_cache_restored": False,
                    "depth_first_intrinsic_memoization": "SYNTHETIC normal per-run memoization",
                    "ordered_result_cache": "not supported",
                },
                "worker_policy": {"owner_preparation_workers": 1},
                "reduction_completed": True,
                "catalogue_closure": {"closed": True, "master_count": 1, "unknown": []},
                "input_sha256": "3" * 64,
                "evaluated_input_sha256": "4" * 64,
                "dimension_binding": {
                    "scope": "SYNTHETIC",
                    "source": "d",
                    "target": "4-2epsilon",
                },
                "data": {"source": "SYNTHETIC DATA"},
                "normalization": "MSbar",
                "loading_seconds": 2,
                "retained_reduction_seconds": 10 if engine == "dfs" else 5,
                "reduction_math_seconds": 4 if engine == "dfs" else 1,
                "reduction_math_timing": {
                    "scope_id": merge.MATH_SCOPE,
                    "complete": True,
                    "gross_seconds": 5 if engine == "dfs" else 2,
                    "measured_owner_preparation_seconds": 1,
                    "interval_count": 2 if engine == "dfs" else 1,
                    "incompleteness_reason": None,
                },
                "completed_seconds": 12 if engine == "dfs" else 7,
                "retained": {
                    "file": payload.name,
                    "bytes": payload.stat().st_size,
                    "sha256": retained[engine],
                },
                "memory_after_retained_reduction": {"process_peak_rss_bytes": None},
                "numerical_reference": None,
                "producer": self.producer,
            }
            report_path = self.write(f"{name}/{engine}/report.json", report)
            reports[engine] = merge.evidence(report_path)
        return self.write(
            f"{name}/parity.json",
            {
                "schema": "vakint.ordered-rung-parity.v1",
                "exact_retained_equality": True,
                "numerical_enclosures_agree": None,
                "input_sha256": "3" * 64,
                "reports": reports,
                "retained": retained,
                "equality_scope": "SYNTHETIC OWNER ASSERTION ONLY",
            },
        )

    def triplet(self, name):
        path = self.pair(name)
        parity = merge.read(path)
        parity["schema"] = "vakint.ordered-rung-parity.v2"
        parity["intrinsic_memoization"] = merge.MEMOIZATION.copy()
        for arm in ("dfs", "ordered"):
            report_path = Path(parity["reports"][arm]["path"])
            raw = merge.read(report_path)
            raw.pop("dfs_request_memoization")
            raw["intrinsic_memoization"] = merge.MEMOIZATION[arm]
            raw["cache_policy"]["depth_first_intrinsic_memoization"] = (
                merge.MEMOIZATION[arm]
            )
            report_path.write_text(json.dumps(raw))
            parity["reports"][arm] = merge.evidence(report_path)
        directory = self.root / name / "dfs-no-memo"
        directory.mkdir()
        payload = directory / "retained.symbolica"
        payload.write_text("SYNTHETIC NO-MEMO FILE CONTRACT; NOT A REDUCTION")
        raw = copy.deepcopy(merge.read(Path(parity["reports"]["dfs"]["path"])))
        raw.update(
            strategy="DepthFirstWithoutMemoization",
            intrinsic_memoization=False,
            reduction_math_seconds=2,
        )
        raw["cache_policy"]["depth_first_intrinsic_memoization"] = False
        raw["reduction_math_timing"].update(gross_seconds=3)
        raw["retained"] = {
            "file": payload.name,
            "bytes": payload.stat().st_size,
            "sha256": merge.digest(payload),
        }
        report_path = self.write(f"{name}/dfs-no-memo/report.json", raw)
        parity["reports"]["dfs_no_memo"] = merge.evidence(report_path)
        parity["retained"]["dfs_no_memo"] = merge.digest(payload)
        path.write_text(json.dumps(parity))
        return path

    def protocol(self, parity=None):
        parity = merge.read(parity or self.parity)
        binary = self.root / "SYNTHETIC-NOT-AN-EXECUTABLE"
        binary.write_text("SYNTHETIC")
        executable = merge.evidence(binary)
        sources = {
            "main_revision": "1" * 40,
            "main_source_tree_sha256": "4" * 64,
            "native_revision": "5" * 40,
            "native_source_tree_sha256": "6" * 64,
            "kernel_revision": "7" * 40,
            "cargo_lock_sha256": "2" * 64,
        }
        log = self.root / "SYNTHETIC-BUILD-LOG"
        log.write_text("SYNTHETIC FILE CONTRACT ONLY")
        build = self.write(
            "build-receipt.json",
            {
                "schema": "vakint.ordered-rung-build.v1",
                "success": True,
                "build_profile": "release",
                "profile_overrides": False,
                "command": ["cargo", "build", "--locked", "--release"],
                "build_log": merge.evidence(log),
                "executable": executable,
                "sources": sources,
            },
        )
        threads = {
            key: "1"
            for key in ("RAYON_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
        }
        samples = {}
        for index, name in enumerate(parity["reports"]):
            report = Path(parity["reports"][name]["path"])
            launch = self.write(
                f"guard/{name}/launch.json",
                {
                    "schema": "vakint.ordered-rung-process.v1",
                    "pid": 100 + index,
                    "start_ticks": 1000,
                    "guard_armed_before_native_exec": True,
                    "executable": executable,
                    "cpu_affinity": [160],
                    "thread_environment": threads,
                    "utc": f"2026-10-29T12:0{index}:00+00:00",
                    "command": [
                        str(binary),
                        "--engine",
                        name.replace("_", "-"),
                        "--loops",
                        "2",
                        "--projection",
                        "zero-momentum",
                        "--build-profile",
                        "release",
                        "--no-monitor",
                        "--work-dir",
                        str(report.parent.parent),
                    ],
                },
            )
            completion = self.write(
                f"guard/{name}/completion.json",
                {
                    "utc": f"2026-10-29T12:0{index}:30+00:00",
                    "success": True,
                    "exit_code": 0,
                    "failure": None,
                    "signals_sent_to_original_pidfd": [],
                },
            )
            samples[name] = {
                "report_sha256": parity["reports"][name]["sha256"],
                "pid": 100 + index,
                "start_ticks": 1000,
                "exit_code": 0,
                "executable_sha256": executable["sha256"],
                "launch": merge.evidence(launch),
                "completion": merge.evidence(completion),
            }
        return self.write(
            "protocol.json",
            {
                "schema": "vakint.ordered-rung-protocol.v2"
                if "dfs_no_memo" in parity["reports"]
                else "vakint.ordered-rung-protocol.v1",
                "fresh_processes": True,
                "serialized_execution": True,
                "build_profile": "release",
                "profile_overrides": False,
                "cpu_affinity": "160",
                "thread_environment": threads,
                "executable": executable,
                "sources": sources,
                "build_receipt": merge.evidence(build),
                "samples": [samples],
            },
        )

    def asymmetric(self):
        anchor = self.triplet("asymmetric-anchor")
        receipts = [anchor]
        protocol_path = self.protocol(anchor)
        protocol = merge.read(protocol_path)
        protocol["schema"] = "vakint.ordered-rung-protocol.v3"
        protocol["sample_counts"] = {"dfs": 1, "dfs_no_memo": 1, "ordered": 3}
        initial = protocol["samples"][0]["ordered"]
        for index in (1, 2):
            receipt = merge.read(anchor)
            source = Path(receipt["reports"]["ordered"]["path"])
            directory = self.root / f"asymmetric-ordered-{index}" / "ordered"
            directory.mkdir(parents=True)
            (directory / "retained.symbolica").write_bytes(
                (source.parent / "retained.symbolica").read_bytes()
            )
            report = merge.read(source)
            report["reduction_math_seconds"] = index + 1
            report["reduction_math_timing"]["gross_seconds"] = index + 2
            new_report = self.write(
                str(directory.relative_to(self.root) / "report.json"), report
            )
            receipt["reports"]["ordered"] = merge.evidence(new_report)
            receipt_path = self.write(f"asymmetric-parity-{index}.json", receipt)
            receipts.append(receipt_path)
            launch = merge.read(Path(initial["launch"]["path"]))
            launch.update(
                pid=200 + index,
                start_ticks=2000 + index,
                utc=f"2026-10-29T12:0{index + 2}:00+00:00",
            )
            launch["command"][launch["command"].index("--work-dir") + 1] = str(
                directory.parent
            )
            launch_path = self.write(f"extra-guard-{index}/launch.json", launch)
            completion = merge.read(Path(initial["completion"]["path"]))
            completion["utc"] = f"2026-10-29T12:0{index + 2}:30+00:00"
            completion_path = self.write(
                f"extra-guard-{index}/completion.json", completion
            )
            protocol["samples"].append(
                {
                    "ordered": {
                        **initial,
                        "report_sha256": merge.digest(new_report),
                        "pid": 200 + index,
                        "start_ticks": 2000 + index,
                        "launch": merge.evidence(launch_path),
                        "completion": merge.evidence(completion_path),
                    }
                }
            )
        protocol_path.write_text(json.dumps(protocol))
        return receipts, protocol_path

    def historical_pair(self):
        pair = merge.read(self.parity)
        ordered_path = Path(pair["reports"]["ordered"]["path"])
        ordered = merge.read(ordered_path)
        ordered.update(loop_count=5, gluon_rungs=4, normalization="MSbar")
        ordered_path.write_text(json.dumps(ordered))
        payload = self.root / "SYNTHETIC-HISTORICAL-RESULT"
        payload.write_text("SYNTHETIC HISTORICAL PHYSICS-FREE FIXTURE")
        closure = {"closed": True, "master_count": 55, "unknown_count": 0}
        baseline = self.write(
            "historical-report.json",
            {
                "schema": "vakint.qcd-ladder.benchmark.v1",
                "reduction_completed": True,
                "retained_result_available": True,
                "catalogue_closure": closure,
                "projection": "zero-momentum-metric",
                "input_sha256": "3" * 64,
                "result_sha256": merge.digest(payload),
                "result_bytes": payload.stat().st_size,
                "reduction_seconds": 100,
                "source_commit": "8" * 40,
            },
        )
        acceptance = self.write(
            "historical-acceptance.json",
            {
                "schema": "vakint.physical-q0.cold-completion-parent-acceptance.v1",
                "status": "ACCEPTED_FUNCTIONAL_COLD_COMPLETE_CACHE_CLOSED",
                "catalogue_closure": closure,
                "collector_reduction_seconds": 100,
                "measurement_scope": "SYNTHETIC FILE CONTRACT ONLY",
                "origin": {"data_manifest_sha256": "9" * 64},
                "output_artifacts_verified": {
                    name: {"sha256": merge.digest(path), "bytes": path.stat().st_size}
                    for name, path in (
                        ("benchmark.json", baseline),
                        ("result.symbolica", payload),
                    )
                },
            },
        )
        definition = self.write(
            "historical-definition.json",
            {
                "schema": "vakint.ordered-rungs.legacy-definition-identity.v1",
                "status": "PASS",
                "normalization": "pySecDec",
                "historical_scalar_input": {"sha256": "3" * 64},
                "historical_result_sha256": merge.digest(payload),
            },
        )
        converted = self.root / "SYNTHETIC-CONVERTED-RESULT"
        converted.write_text("SYNTHETIC CONVERSION, NOT A MATHEMATICAL EXPRESSION")
        factor = "SYNTHETIC_N_TARGET_OVER_SOURCE"
        return self.write(
            "historical-parity.json",
            {
                "schema": "vakint.ordered-rung-historical-parity.v1",
                "exact_retained_equality": True,
                "input_sha256": "3" * 64,
                "historical": {
                    "acceptance": merge.evidence(acceptance),
                    "report": merge.evidence(baseline),
                    "retained": merge.evidence(payload),
                    "definition": merge.evidence(definition),
                },
                "ordered_report": merge.evidence(ordered_path),
                "normalization_conversion": {
                    "schema": "vakint.ordered-rung-normalization-conversion.v1",
                    "source_normalization": "pySecDec",
                    "target_normalization": "MSbar",
                    "api": "VakintSettings::get_integral_normalization_factor_atom",
                    "ratio_identity": "N_target/N_source",
                    "loop_count": 5,
                    "dimension_binding": ordered["dimension_binding"],
                    "epsilon_symbol": "SYNTHETIC_EPSILON",
                    "mu_r_sq_symbol": "SYNTHETIC_MU_R_SQ",
                    "source_retained_sha256": merge.digest(payload),
                    "target_report_sha256": merge.digest(ordered_path),
                    "factor": {
                        "canonical_expression": factor,
                        "sha256": hashlib.sha256(factor.encode()).hexdigest(),
                    },
                    "converted_retained": {
                        **merge.evidence(converted),
                        "bytes": converted.stat().st_size,
                    },
                    "exact_conversion_verified": True,
                    "applied_before_laurent_truncation": True,
                },
                "retained": {
                    "historical": merge.digest(payload),
                    "historical_converted": merge.digest(converted),
                    "ordered": ordered["retained"]["sha256"],
                },
            },
        )

    def test_append_preserves_history_and_null_stages(self):
        original = copy.deepcopy(self.base)
        result = merge.append(self.base, "test", [self.parity])
        self.assertEqual(self.base, original)
        for key in original:
            self.assertEqual(result[key], original[key])
        row = result["ordered_rung_runs"]["test"]["summary"][0]
        self.assertIsNone(row["ordered_over_dfs"])
        self.assertFalse(row["performance_eligible"])
        self.assertEqual(row["sample_count"], 1)
        self.assertIsNone(row["metrics"]["ordered"]["pure_ibp_seconds"])
        self.assertEqual(row["metrics"]["ordered"]["ordered_progress"], [None])

    def test_bound_fresh_protocol_qualifies_observed_ratio(self):
        result = merge.append(self.base, "test", [self.parity], self.protocol())
        row = result["ordered_rung_runs"]["test"]["summary"][0]
        self.assertTrue(row["performance_eligible"])
        self.assertEqual(row["ordered_over_dfs"], 0.25)
        self.assertEqual(
            row["metrics"]["dfs"]["public_retained_seconds"]["median_seconds"], 10
        )
        self.assertEqual(
            row["metrics"]["dfs"]["reduction_math_seconds"]["median_seconds"], 4
        )

    def test_changed_report_is_rejected(self):
        raw = merge.read(self.parity)
        path = Path(raw["reports"]["ordered"]["path"])
        path.write_text(path.read_text() + " ")
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            merge.append(self.base, "test", [self.parity])

    def test_missing_math_measurement_never_falls_back_to_public_timer(self):
        for name in ("dfs", "ordered"):
            self.update_report(
                name,
                lambda row: (
                    row.pop("reduction_math_seconds"),
                    row.pop("reduction_math_timing"),
                ),
            )
        result = merge.append(self.base, "test", [self.parity])
        row = result["ordered_rung_runs"]["test"]["summary"][0]
        self.assertIsNone(row["metrics"]["dfs"]["reduction_math_seconds"])
        self.assertIsNone(row["ordered_over_dfs"])
        with self.assertRaisesRegex(ValueError, "complete reduction-math telemetry"):
            merge.append(self.base, "test", [self.parity], self.protocol())

    def test_monitoring_is_diagnostic_only_for_wall_time_comparison(self):
        for name in ("dfs", "ordered"):
            self.update_report(name, lambda row: row.update(monitoring_enabled=True))
        result = merge.append(self.base, "test", [self.parity])
        self.assertIsNone(
            result["ordered_rung_runs"]["test"]["summary"][0]["ordered_over_dfs"]
        )
        with self.assertRaisesRegex(ValueError, "requires monitoring disabled"):
            merge.append(self.base, "test", [self.parity], self.protocol())

    def test_incomplete_math_telemetry_preserves_physics_without_net_time(self):
        self.update_report(
            "ordered",
            lambda row: row.update(
                reduction_math_seconds=None,
                reduction_math_timing={
                    "scope_id": merge.MATH_SCOPE,
                    "complete": False,
                    "gross_seconds": None,
                    "measured_owner_preparation_seconds": None,
                    "interval_count": 0,
                    "incompleteness_reason": "SYNTHETIC feature disabled",
                },
            ),
        )
        result = merge.append(self.base, "test", [self.parity])
        row = result["ordered_rung_runs"]["test"]["summary"][0]
        self.assertTrue(row["exact_retained_equality"])
        self.assertIsNone(row["metrics"]["ordered"]["reduction_math_seconds"])
        with self.assertRaisesRegex(ValueError, "complete reduction-math telemetry"):
            merge.append(self.base, "test", [self.parity], self.protocol())

    def test_complete_math_deduction_scope_and_workers_are_checked(self):
        self.update_report(
            "ordered",
            lambda row: row["reduction_math_timing"].update(
                measured_owner_preparation_seconds=3
            ),
        )
        with self.assertRaisesRegex(ValueError, "invalid owner-preparation deduction"):
            merge.append(self.base, "test", [self.parity])
        self.parity = self.pair("timing-scope")
        self.update_report(
            "ordered",
            lambda row: row["reduction_math_timing"].update(
                scope_id="SYNTHETIC_OTHER_SCOPE"
            ),
        )
        with self.assertRaisesRegex(ValueError, "invalid reduction-math scope"):
            merge.append(self.base, "test", [self.parity])
        self.parity = self.pair("timing-workers")
        self.update_report(
            "ordered",
            lambda row: row["worker_policy"].update(owner_preparation_workers=2),
        )
        with self.assertRaisesRegex(ValueError, "worker policy"):
            merge.append(self.base, "test", [self.parity])

    def test_changed_retained_payload_is_rejected(self):
        raw = merge.read(self.parity)
        path = Path(raw["reports"]["ordered"]["path"]).parent / "retained.symbolica"
        path.write_text("SYNTHETIC CHANGED RESULT")
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            merge.append(self.base, "test", [self.parity])

    def test_incomplete_or_unclosed_reduction_is_rejected(self):
        self.update_report("ordered", lambda row: row.update(reduction_completed=False))
        with self.assertRaisesRegex(ValueError, "incomplete"):
            merge.append(self.base, "test", [self.parity])
        self.parity = self.pair("unclosed")
        self.update_report(
            "ordered",
            lambda row: row["catalogue_closure"].update(unknown=["SYNTHETIC"]),
        )
        with self.assertRaisesRegex(ValueError, "catalogue closure"):
            merge.append(self.base, "test", [self.parity])

    def test_evaluated_input_and_dimension_binding_are_compared(self):
        self.update_report(
            "ordered", lambda row: row.update(evaluated_input_sha256="9" * 64)
        )
        with self.assertRaisesRegex(ValueError, "evaluated_input_sha256 differs"):
            merge.append(self.base, "test", [self.parity])
        self.parity = self.pair("dimension")
        self.update_report(
            "ordered",
            lambda row: row["dimension_binding"].update(target="SYNTHETIC_OTHER"),
        )
        with self.assertRaisesRegex(ValueError, "dimension_binding differs"):
            merge.append(self.base, "test", [self.parity])

    def test_paired_numerical_points_and_orders_are_compared(self):
        self.add_numerical()
        merge.append(self.base, "test", [self.parity])
        self.update_report(
            "ordered", lambda row: row["numerical_reference"]["point"].update(mUV=2)
        )
        with self.assertRaisesRegex(ValueError, "numerical point differs"):
            merge.append(self.base, "test", [self.parity])
        self.parity = self.pair("orders")
        self.add_numerical()
        self.update_report(
            "ordered",
            lambda row: row["numerical_reference"].update(
                coefficients=[{"epsilon_power": 0}]
            ),
        )
        with self.assertRaisesRegex(ValueError, "Laurent orders differ"):
            merge.append(self.base, "test", [self.parity])

    def test_development_profile_is_diagnostic_with_nullable_self_hash(self):
        for name in ("dfs", "ordered"):
            self.update_report(
                name,
                lambda row: row["producer"].update(
                    build_profile="dev", debug_assertions=True, executable_sha256=None
                ),
            )
        result = merge.append(self.base, "test", [self.parity])
        self.assertIsNone(
            result["ordered_rung_runs"]["test"]["summary"][0]["ordered_over_dfs"]
        )
        with self.assertRaisesRegex(ValueError, "fresh-process release cohort"):
            merge.append(self.base, "test", [self.parity], self.protocol())

    def test_equivalent_exact_source_inventory_may_differ_but_cap_must_match(self):
        self.add_numerical()
        self.update_report(
            "ordered",
            lambda row: row["numerical_reference"]["selection"].update(
                surviving_dependencies=[
                    {
                        "identity": "SYNTHETIC e",
                        "available_precision": "ArbitraryPrecision",
                    }
                ]
            ),
        )
        merge.append(self.base, "test", [self.parity])
        self.update_report(
            "ordered",
            lambda row: row["numerical_reference"]["selection"].update(
                surviving_dependencies=[
                    {
                        "identity": "SYNTHETIC finite",
                        "available_precision": "FixedPrecision(300)",
                    }
                ]
            ),
        )
        with self.assertRaisesRegex(ValueError, "effective source cap differs"):
            merge.append(self.base, "test", [self.parity])

    def test_cached_run_and_failed_parity_are_rejected(self):
        raw = merge.read(self.parity)
        raw["exact_retained_equality"] = False
        self.parity.write_text(json.dumps(raw))
        with self.assertRaisesRegex(ValueError, "successful exact parity"):
            merge.append(self.base, "test", [self.parity])

        self.parity = self.pair("cached-sample")
        raw = merge.read(self.parity)
        path = Path(raw["reports"]["ordered"]["path"])
        report = merge.read(path)
        report["persistent_cache"] = True
        path.write_text(json.dumps(report))
        raw["reports"]["ordered"] = merge.evidence(path)
        self.parity.write_text(json.dumps(raw))
        with self.assertRaisesRegex(ValueError, "attaches a cache"):
            merge.append(self.base, "test", [self.parity])

    def test_same_process_cannot_be_qualified_as_fresh(self):
        path = self.protocol()
        raw = merge.read(path)
        raw["samples"][0]["ordered"]["pid"] = raw["samples"][0]["dfs"]["pid"]
        path.write_text(json.dumps(raw))
        with self.assertRaisesRegex(ValueError, "reuses a worker"):
            merge.append(self.base, "test", [self.parity], path)

    def test_bad_build_receipt_never_downgrades_to_diagnostic(self):
        path = self.protocol()
        protocol = merge.read(path)
        build_path = Path(protocol["build_receipt"]["path"])
        build = merge.read(build_path)
        build["success"] = False
        build_path.write_text(json.dumps(build))
        protocol["build_receipt"] = merge.evidence(build_path)
        path.write_text(json.dumps(protocol))
        with self.assertRaisesRegex(ValueError, "build receipt"):
            merge.append(self.base, "test", [self.parity], path)

    def test_guard_engine_and_serialized_windows_are_checked(self):
        path = self.protocol()
        protocol = merge.read(path)
        worker = protocol["samples"][0]["ordered"]
        launch_path = Path(worker["launch"]["path"])
        launch = merge.read(launch_path)
        launch["command"][launch["command"].index("--engine") + 1] = "both"
        launch_path.write_text(json.dumps(launch))
        worker["launch"] = merge.evidence(launch_path)
        path.write_text(json.dumps(protocol))
        with self.assertRaisesRegex(ValueError, "guard command"):
            merge.append(self.base, "test", [self.parity], path)
        path = self.protocol()
        protocol = merge.read(path)
        worker = protocol["samples"][0]["ordered"]
        launch_path = Path(worker["launch"]["path"])
        launch = merge.read(launch_path)
        launch["utc"] = "2026-10-29T12:00:10+00:00"
        launch_path.write_text(json.dumps(launch))
        worker["launch"] = merge.evidence(launch_path)
        path.write_text(json.dumps(protocol))
        with self.assertRaisesRegex(ValueError, "windows overlap"):
            merge.append(self.base, "test", [self.parity], path)

    def test_repeated_samples_keep_cohort_and_compute_real_median(self):
        second = self.pair("sample-2")
        result = merge.append(self.base, "test", [self.parity, second])
        row = result["ordered_rung_runs"]["test"]["summary"][0]
        self.assertEqual(row["sample_count"], 2)
        self.assertEqual(
            row["metrics"]["dfs"]["public_retained_seconds"]["samples"], [10, 10]
        )
        with self.assertRaisesRegex(ValueError, "duplicate parity receipt"):
            merge.append(self.base, "test", [self.parity, self.parity])

    def test_historical_reference_stays_separate_with_no_ratio(self):
        path = self.historical_pair()
        result = merge.append(self.base, "test", [], historical_receipts=[path])
        run = result["ordered_rung_runs"]["test"]
        self.assertEqual(run["summary"], [])
        reference = run["historical_comparisons"][0]
        self.assertIsNone(reference["ordered_over_dfs"])
        self.assertFalse(reference["performance_eligible"])
        self.assertIsNone(reference["historical"]["reduction_math_seconds"])
        self.assertEqual(reference["historical"]["public_retained_seconds"], 100)
        self.assertEqual(result["historical"], self.base["historical"])

    def test_new_pysecdec_reports_are_rejected(self):
        self.update_report("ordered", lambda row: row.update(normalization="pySecDec"))
        with self.assertRaisesRegex(ValueError, "require MSbar"):
            merge.append(self.base, "test", [self.parity])

    def test_historical_conversion_proof_and_factor_digest_are_required(self):
        path = self.historical_pair()
        raw = merge.read(path)
        proof = raw.pop("normalization_conversion")
        path.write_text(json.dumps(raw))
        with self.assertRaisesRegex(ValueError, "conversion proof"):
            merge.append(self.base, "test", [], historical_receipts=[path])
        raw["normalization_conversion"] = proof
        proof["factor"]["canonical_expression"] = "SYNTHETIC CHANGED FACTOR"
        path.write_text(json.dumps(raw))
        with self.assertRaisesRegex(ValueError, "factor identity differs"):
            merge.append(self.base, "test", [], historical_receipts=[path])

    def test_triplet_has_three_actual_arms_and_two_net_math_ratios(self):
        self.parity = self.triplet("triplet")
        loaded = merge.load_triplet(self.parity)
        self.assertEqual(tuple(loaded["reports"]), ("dfs", "dfs_no_memo", "ordered"))
        result = merge.append(self.base, "three", [self.parity], self.protocol())
        run = result["ordered_rung_runs"]["three"]
        self.assertEqual(run["schema"], "vakint.ordered-rung-report.v2")
        row = run["summary"][0]
        self.assertEqual(row["ordered_over_dfs"], 0.25)
        self.assertEqual(row["ordered_over_dfs_no_memo"], 0.5)
        self.assertFalse(row["intrinsic_memoization"]["dfs_no_memo"])
        self.assertEqual(
            row["metrics"]["dfs_no_memo"]["reduction_math_seconds"]["median_seconds"], 2
        )

    def test_triplet_cannot_label_an_ordinary_dfs_as_no_memo(self):
        self.parity = self.triplet("wrong-policy")
        self.update_report("dfs_no_memo", lambda row: row.update(strategy="DepthFirst"))
        with self.assertRaisesRegex(ValueError, "strategy"):
            merge.load_triplet(self.parity)
        self.parity = self.triplet("wrong-bool")
        self.update_report(
            "dfs_no_memo", lambda row: row.update(intrinsic_memoization=True)
        )
        with self.assertRaisesRegex(ValueError, "memoization"):
            merge.load_triplet(self.parity)

    def test_triplet_requires_the_no_memo_arm_and_keeps_v1_history(self):
        original = merge.append(self.base, "original-pair", [self.parity])
        self.parity = self.triplet("new-three")
        result = merge.append(original, "new-three", [self.parity])
        self.assertEqual(
            result["ordered_rung_runs"]["original-pair"],
            original["ordered_rung_runs"]["original-pair"],
        )
        raw = merge.read(self.parity)
        del raw["reports"]["dfs_no_memo"]
        self.parity.write_text(json.dumps(raw))
        with self.assertRaisesRegex(ValueError, "omits or adds an arm"):
            merge.load_triplet(self.parity)

    def test_append_does_not_requalify_an_existing_historical_v1_run(self):
        protocol = self.protocol()
        original = merge.append(self.base, "old-qualified", [self.parity], protocol)
        binary = Path(merge.read(protocol)["executable"]["path"])
        binary.write_text("SYNTHETIC LATER BUILD OVERWROTE THE ORIGINAL TARGET PATH")
        self.parity = self.triplet("later-triplet")
        result = merge.append(original, "new-diagnostic", [self.parity])
        self.assertTrue(
            result["ordered_rung_runs"]["old-qualified"]["summary"][0][
                "performance_eligible"
            ]
        )
        self.assertEqual(
            result["ordered_rung_runs"]["old-qualified"],
            original["ordered_rung_runs"]["old-qualified"],
        )

    def test_asymmetric_counts_use_one_real_baseline_and_three_ordered_values(self):
        receipts, protocol = self.asymmetric()
        result = merge.append(self.base, "asymmetric", receipts, protocol)
        run = result["ordered_rung_runs"]["asymmetric"]
        self.assertEqual(run["schema"], "vakint.ordered-rung-report.v3")
        row = run["summary"][0]
        self.assertEqual(
            row["sample_counts"], {"dfs": 1, "dfs_no_memo": 1, "ordered": 3}
        )
        self.assertIsNone(row["sample_count"])
        self.assertEqual(
            row["metrics"]["dfs"]["reduction_math_seconds"]["samples"], [4]
        )
        self.assertEqual(
            row["metrics"]["dfs_no_memo"]["reduction_math_seconds"]["samples"], [2]
        )
        self.assertEqual(
            row["metrics"]["ordered"]["reduction_math_seconds"]["samples"], [1, 2, 3]
        )
        self.assertEqual(row["ordered_over_dfs"], 0.5)
        self.assertEqual(row["ordered_over_dfs_no_memo"], 1)

    def test_repeated_baseline_references_need_explicit_asymmetric_protocol(self):
        receipts, _ = self.asymmetric()
        with self.assertRaisesRegex(ValueError, "explicit asymmetric protocol"):
            merge.append(self.base, "wrong", receipts)

    def test_asymmetric_protocol_rejects_fake_counts_or_extra_baseline_worker(self):
        receipts, path = self.asymmetric()
        protocol = merge.read(path)
        protocol["sample_counts"]["dfs"] = 3
        path.write_text(json.dumps(protocol))
        with self.assertRaisesRegex(ValueError, "one real baseline"):
            merge.append(self.base, "wrong", receipts, path)
        protocol["sample_counts"]["dfs"] = 1
        protocol["samples"][1]["dfs"] = protocol["samples"][0]["dfs"]
        path.write_text(json.dumps(protocol))
        with self.assertRaisesRegex(ValueError, "omits or adds an arm"):
            merge.append(self.base, "wrong", receipts, path)

    def diagnostic_receipts(self):
        parity = merge.read(self.parity)
        output = {}
        for arm in ("dfs", "ordered"):
            report_path = Path(parity["reports"][arm]["path"])
            raw = merge.read(report_path)
            raw.update(loop_count=4, gluon_rungs=3)
            raw["catalogue_closure"]["unknown_count"] = 0
            report_path.write_text(json.dumps(raw))
            completion = self.write(
                f"diagnostic-{arm}-completion.json",
                {
                    "success": True,
                    "exit_code": 0,
                    "failure": None,
                    "signals_sent_to_original_pidfd": [],
                    "elapsed_seconds": 20,
                },
            )
            output[arm] = self.write(
                f"diagnostic-{arm}.json",
                {
                    "schema": "vakint.ordered-rung-functional-probe-receipt.v1",
                    "status": "PASS-functional-catalogue-closure-only",
                    "success": True,
                    "main_revision": "1" * 40,
                    "report": merge.evidence(report_path),
                    "guard_completion": merge.evidence(completion),
                },
            )
        return output

    def test_completed_diagnostics_do_not_infer_parity_or_ratios(self):
        result = merge.append_diagnostic(
            self.base,
            "pending",
            self.diagnostic_receipts(),
            "SYNTHETIC pending exact comparison",
        )
        row = result["ordered_rung_diagnostics"]["pending"]
        self.assertIsNone(row["exact_retained_equality"])
        self.assertIsNone(row["numerical_validation"])
        self.assertIsNone(row["ordered_over_dfs"])
        self.assertFalse(row["performance_eligible"])
        self.assertEqual(row["sample_counts"]["dfs_no_memo"], 0)
        self.assertEqual(result["historical"], self.base["historical"])

    def test_diagnostic_changed_payload_or_failed_worker_is_rejected(self):
        observations = self.diagnostic_receipts()
        receipt = merge.read(observations["ordered"])
        report_path = Path(receipt["report"]["path"])
        (report_path.parent / "retained.symbolica").write_text("SYNTHETIC CORRUPTION")
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            merge.append_diagnostic(self.base, "bad", observations, "SYNTHETIC")

    def test_functional_ci_diagnostic_requires_named_pass_proof(self):
        self.add_numerical()
        observations = self.diagnostic_receipts()
        records = {arm: merge.read(path) for arm, path in observations.items()}
        summary = self.write(
            "functional-ci.json",
            {
                "schema": "vakint.four-loop-functional-ci-case-completion.v1",
                "status": "PASS",
                "performance_eligible": False,
                "test": "four_loop_zero_momentum_rung_engine_parity",
                "source_revision": "1" * 40,
                "completion": records["dfs"]["guard_completion"],
                "reports": {arm: {"report": records[arm]["report"]} for arm in records},
            },
        )
        row = merge.append_diagnostic(
            self.base,
            "functional",
            {arm: summary for arm in observations},
            "SYNTHETIC",
            functional_ci=summary,
        )["ordered_rung_diagnostics"]["functional"]
        self.assertTrue(row["exact_retained_equality"])
        self.assertTrue(row["numerical_validation"])
        self.assertFalse(row["performance_eligible"])
        self.assertIsNone(row["ordered_over_dfs"])

    def test_diagnostic_failed_worker_cannot_supply_observed_timings(self):
        observations = self.diagnostic_receipts()
        receipt = merge.read(observations["ordered"])
        path = Path(receipt["guard_completion"]["path"])
        completion = merge.read(path)
        completion.update(success=False, exit_code=1)
        path.write_text(json.dumps(completion))
        receipt["guard_completion"] = merge.evidence(path)
        observations["ordered"].write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError, "did not complete successfully"):
            merge.append_diagnostic(self.base, "failed", observations, "SYNTHETIC")

    def exact_diagnostic_fixture(self):
        pending = merge.append_diagnostic(
            self.base, "pending", self.diagnostic_receipts(), "SYNTHETIC pending"
        )
        reports = pending["ordered_rung_diagnostics"]["pending"]["observations"]
        proof = self.write(
            "saved-proof/proof/functional-pair.json",
            {
                "schema": "vakint.ordered-rung-functional-pair.v1",
                "exact_retained_equality": True,
                "functional_only": True,
                "performance_eligible": False,
                "new_evaluations": 0,
                "retain_only_explicitly_requested": True,
                "numerical_enclosures_agree": None,
                "numerical_check": "not checked",
                "input_sha256": "3" * 64,
                "reports": {arm: reports[arm]["report"]["evidence"] for arm in reports},
                "retained": {
                    arm: reports[arm]["report"]["raw"]["retained"]["sha256"]
                    for arm in reports
                },
            },
        )
        self.write(
            "saved-proof/guard/launch.json",
            {
                "schema": "vakint.ordered-rung-process.v1",
                "guard_armed_before_native_exec": True,
                "command": [
                    "SYNTHETIC_COMPARATOR",
                    "--compare-two-reports",
                    *[
                        str(Path(reports[arm]["report"]["evidence"]["path"]).parent)
                        for arm in ("dfs", "ordered")
                    ],
                    "--no-numerical-reference",
                ],
            },
        )
        self.write(
            "saved-proof/guard/completion.json",
            {
                "success": True,
                "exit_code": 0,
                "failure": None,
                "signals_sent_to_original_pidfd": [],
            },
        )
        return pending, proof

    def test_exact_diagnostic_preserves_pending_entry_and_no_numeric_claim(self):
        pending, proof = self.exact_diagnostic_fixture()
        result = merge.append_exact_diagnostic(
            pending, "exact", "pending", proof, "SYNTHETIC exact only"
        )
        self.assertEqual(
            result["ordered_rung_diagnostics"]["pending"],
            pending["ordered_rung_diagnostics"]["pending"],
        )
        row = result["ordered_rung_diagnostics"]["exact"]
        self.assertTrue(row["exact_retained_equality"])
        self.assertIsNone(row["numerical_validation"])
        self.assertIsNone(row["ordered_over_dfs"])
        self.assertEqual(
            row["observations"],
            pending["ordered_rung_diagnostics"]["pending"]["observations"],
        )

    def test_exact_diagnostic_rejects_unbound_output_or_unchecked_numeric_claim(self):
        pending, proof = self.exact_diagnostic_fixture()
        raw = merge.read(proof)
        raw["retained"]["ordered"] = "9" * 64
        proof.write_text(json.dumps(raw))
        with self.assertRaisesRegex(ValueError, "does not bind"):
            merge.append_exact_diagnostic(pending, "bad", "pending", proof, "SYNTHETIC")
        raw["retained"]["ordered"] = pending["ordered_rung_diagnostics"]["pending"][
            "observations"
        ]["ordered"]["report"]["raw"]["retained"]["sha256"]
        raw["numerical_enclosures_agree"] = True
        proof.write_text(json.dumps(raw))
        with self.assertRaisesRegex(ValueError, "exact-only"):
            merge.append_exact_diagnostic(pending, "bad", "pending", proof, "SYNTHETIC")

    def test_current_numeric_checkpoint_is_distinct_from_numeric_pair_comparison(self):
        pending, proof = self.exact_diagnostic_fixture()
        current = pending["ordered_rung_diagnostics"]["pending"]["observations"][
            "ordered"
        ]
        numerical = {
            "coefficients": [
                {
                    "epsilon_power": 0,
                    "real_significant_digits": None,
                    "imaginary_significant_digits": 213,
                    "real_error_quality": "Estimated",
                    "imaginary_error_quality": "Estimated",
                }
            ],
            "warnings": ["SYNTHETIC limiting sources"],
            "selection": {"requested_decimal_precision": 215},
        }
        current["report"]["raw"]["numerical_reference"] = numerical
        current["receipt"]["raw"].update(
            schema="vakint.ordered-rung-functional-diagnostic-completion.v1",
            numerical_precision_by_epsilon=numerical["coefficients"],
        )
        report_path = Path(current["report"]["evidence"]["path"])
        report_path.write_text(json.dumps(current["report"]["raw"]))
        current["report"]["evidence"] = merge.evidence(report_path)
        proof_raw = merge.read(proof)
        proof_raw["reports"]["ordered"] = merge.evidence(report_path)
        proof.write_text(json.dumps(proof_raw))
        receipt_path = Path(current["receipt"]["evidence"]["path"])
        receipt_path.write_text(json.dumps(current["receipt"]["raw"]))
        current["receipt"]["evidence"] = merge.evidence(receipt_path)
        result = merge.append_exact_diagnostic(
            pending, "current", "pending", proof, "SYNTHETIC", current_numerical=True
        )
        checkpoint = result["ordered_rung_diagnostics"]["current"]
        self.assertTrue(checkpoint["numerical_validation"])
        self.assertIsNone(checkpoint["numerical_pair_comparison"])
        self.assertEqual(
            checkpoint["observations"]["ordered"]["report"]["raw"][
                "numerical_reference"
            ]["warnings"],
            numerical["warnings"],
        )
        current["receipt"]["raw"]["numerical_precision_by_epsilon"] = []
        receipt_path.write_text(json.dumps(current["receipt"]["raw"]))
        current["receipt"]["evidence"] = merge.evidence(receipt_path)
        with self.assertRaisesRegex(ValueError, "precision evidence differs"):
            merge.append_exact_diagnostic(
                pending, "bad", "pending", proof, "SYNTHETIC", current_numerical=True
            )

    def test_historical_changed_acceptance_or_timing_is_rejected(self):
        path = self.historical_pair()
        raw = merge.read(path)
        acceptance_path = Path(raw["historical"]["acceptance"]["path"])
        acceptance = merge.read(acceptance_path)
        acceptance["collector_reduction_seconds"] = 101
        acceptance_path.write_text(json.dumps(acceptance))
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            merge.append(self.base, "test", [], historical_receipts=[path])
        raw["historical"]["acceptance"] = merge.evidence(acceptance_path)
        path.write_text(json.dumps(raw))
        with self.assertRaisesRegex(ValueError, "accepted timing differs"):
            merge.append(self.base, "test", [], historical_receipts=[path])


if __name__ == "__main__":
    unittest.main()
