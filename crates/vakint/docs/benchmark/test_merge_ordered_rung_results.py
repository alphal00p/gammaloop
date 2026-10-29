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
        for name in ("dfs", "ordered"):
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
        for index, name in enumerate(("dfs", "ordered")):
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
                        name,
                        "--loops",
                        "2",
                        "--projection",
                        "zero-momentum",
                        "--build-profile",
                        "release",
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
                "schema": "vakint.ordered-rung-protocol.v1",
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
        self.assertEqual(row["ordered_over_dfs"], 0.5)

    def test_changed_report_is_rejected(self):
        raw = merge.read(self.parity)
        path = Path(raw["reports"]["ordered"]["path"])
        path.write_text(path.read_text() + " ")
        with self.assertRaisesRegex(ValueError, "changed evidence"):
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
