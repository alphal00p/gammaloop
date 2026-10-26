"""Synthetic opaque-union receipts only; no Native codec, payload or reductions."""

import copy
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import compose_cumulative_cache_report as composer
import merge_five_loop_results as reports
import test_completed_capture_source as capture_sources
import test_strict_preparation_source as finite_sources


class ComposeUnionSourcesTests(unittest.TestCase):
    def write(self, path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
        return reports.evidence(path)

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="synthetic-union-composition-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.native_root = self.root / "actual-native-output"
        self.cache = self.native_root / "cache"
        self.cache.mkdir(parents=True)
        self.binding = {
            "artifact_digest": [1] * 32,
            "family_fingerprint": "SYNTHETIC",
            "physical_arity": 15,
            "layer": "saved-application-route-apply-v1",
        }
        self.proof_root = self.root / "proof-root"
        self.data_manifest = self.root / "report-context-manifest.json"
        self.write(
            self.data_manifest,
            {"schema": "vakint.external-data.v1", "version": "v99.0.0"},
        )
        self.finite = finite_sources.StrictPreparationSourceTests()
        self.finite.setUp()
        self.addCleanup(self.finite.doCleanups)
        self.source_paths, raw_sources, snapshots, censuses = [], [], [], []
        source_rows = [
            (self.finite.source, self.finite.cache, 3, 4),
            *[self.capture_source(position) for position in (1, 2)],
        ]
        for position, (source, directory, generation, records) in enumerate(
            source_rows
        ):
            self.source_paths.append(self.root / f"logical-source-{position}.json")
            self.write(self.source_paths[-1], source)
            identity = reports.native_cache_identity(directory)
            index = directory / "cache.rrindex"
            census = {
                "generation": generation,
                "runs": 1,
                "chunk_references": 2,
                "stored_record_references": records,
                "referenced_compressed_bytes": 100,
                "manifest_bytes": index.stat().st_size,
            }
            censuses.append(census)
            raw_sources.append(
                {
                    "id": source["id"],
                    "kind": "pack",
                    "directory": str(directory),
                    "cache_native_index_identity": identity,
                    # Unknown prior statistics stay null; the actual census supplies them.
                    "expected": {
                        **census,
                        "runs": None,
                        "chunk_references": None,
                        "referenced_compressed_bytes": None,
                    },
                }
            )
            frozen = self.native_root / f"snapshots/source-{position}.rrindex"
            frozen.parent.mkdir(exist_ok=True)
            frozen.write_bytes(index.read_bytes())
            snapshots.append(
                {
                    "manifest": f"snapshots/source-{position}.rrindex",
                    "generation": generation,
                    "manifest_bytes": census["manifest_bytes"],
                    "manifest_blake3": [position] * 32,
                    "manifest_sha256": reports.digest(frozen),
                    "refresh_allowed": False,
                }
            )
        child = self.cache / "bindings" / ("b" * 64) / "cache.rrindex"
        child.parent.mkdir(parents=True)
        child.write_text("SYNTHETIC UNION CHILD INDEX; NOT NATIVE DATA")
        (self.cache / "bundle.rrindex").write_text("SYNTHETIC UNION BUNDLE INDEX")
        identity = reports.native_cache_identity(self.cache)
        totals = {
            key: sum(row[key] for row in censuses)
            for key in (
                "runs",
                "chunk_references",
                "stored_record_references",
                "referenced_compressed_bytes",
            )
        }
        cohort = {
            field: str(position + 1) * (64 if field.endswith("_sha256") else 40)
            for position, (field, _) in enumerate(reports.COHORT_STAMPS)
        }
        cohort["measurement_executable_sha256"] = "9" * 64
        cohort["stamps"] = {
            stamp: cohort[field] for field, stamp in reports.COHORT_STAMPS
        }
        artifact_name = "provenance/native-exporter-v4-census/collector-artifacts.json"
        artifact = self.write(
            self.proof_root / artifact_name,
            {
                "compiled_cohort": cohort,
                "stamps": cohort["stamps"],
                "binary": {"sha256": "9" * 64},
            },
        )
        index_bytes = sum(row["bytes"] for row in identity["indexes"].values())
        self.raw = {
            "schema": "vakint.native-opaque-union-validation.v1",
            "success": True,
            "request_sha256": "a" * 64,
            "semantic_binding": self.binding,
            "sources_newest_first": raw_sources,
            "captured_snapshots": snapshots,
            "source_censuses": censuses,
            "source_indexes_unchanged": True,
            "frozen_readers_used": True,
            "cache_native_index_identity": identity,
            "output_census": {
                "generation": 1,
                **totals,
                "manifest_bytes": child.stat().st_size,
            },
            "source_reference_totals": totals,
            "copy_statistics": {
                "namespaces": 1,
                "copied_bases": 3,
                "unique_physical_chunks": 4,
                "physical_compressed_bytes": 200,
                "index_bytes": index_bytes,
                "preserved_record_references": totals["stored_record_references"],
                "transfer_buffer_bytes": 65536,
                "charged_disk_bytes": 200 + index_bytes,
            },
            "distinct_keys": None,
            "limits": {
                "max_disk_bytes": 10000,
                "max_manifest_bytes": 10000,
                "max_directory_owned_bytes": 1000000,
            },
            "compiled_cohort": cohort,
            "executable_sha256": "9" * 64,
            "collector_artifacts": {**artifact, "path": artifact_name},
            "resource_guard_admission": {
                "schema": "vakint.resource-guard-admission.v1",
                "admitted": True,
                "worker": {
                    "pid": 11,
                    "start_ticks": 12,
                    "uid": 13,
                    "executable_sha256": "9" * 64,
                },
                "guard": {"pid": 21, "start_ticks": 22, "uid": 13},
                "launch_spec_sha256": "6" * 64,
                "guard_config_sha256": "7" * 64,
                "slot": "SYNTHETIC",
                "admission_utc": "SYNTHETIC",
                "MemAvailable_bytes": 100000000,
            },
            "compiled_codec_state_format": 6,
            "coefficient_context_created": False,
            "state_or_coefficient_imports": 0,
            "IBP_applications": 0,
        }
        self.raw_path = self.native_root / "raw-validation.json"
        self.counter = 0

    def capture_source(self, position):
        root = self.root / f"capture-{position}"
        directory = root / "native-pack"
        directory.mkdir(parents=True)
        index = directory / "cache.rrindex"
        index.write_text(f"SYNTHETIC CAPTURE INDEX {position}; NOT NATIVE DATA")
        generation, records = 4 + position, 5 + position
        producer = self.write(
            root / "worker.json",
            {
                "schema": "vakint.five-loop-cache-measurement.v1",
                "success": True,
                "case_id": f"original-tuning-{position}",
                "provisional": True,
                "native_keys_authenticated": False,
            },
        )
        native = {
            "schema": "rustred.running-application-cache.v1",
            "generation": generation,
            "accepted_records": records,
            "committed_records": records,
            "pending_records": 0,
            "pending_bytes": 0,
            "writer_error": None,
            "last_error": None,
            "closed": False,
        }
        capture = self.write(root / "running-cache.json", {"report": native})
        snapshot = self.write(
            root / "snapshot.json",
            {
                "schema": "rustred.running-application-cache-snapshot.v1",
                "delta": {
                    "generation": generation,
                    "manifest_digest": bytes([position] * 32).hex(),
                    "manifest_path": str(index),
                    "manifest_bytes": index.stat().st_size,
                },
            },
        )
        barrier = self.write(
            root / "barrier.json",
            {
                "schema": "vakint.native-cache-snapshot-barrier.v1",
                "accepted": True,
                "committed_snapshot": True,
                "producer_completed": False,
                "producer": producer,
                "capture": capture,
                "snapshot": snapshot,
                "writer_fault": None,
                "generation": generation,
                "manifest_digest": bytes([position] * 32).hex(),
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
        source = {
            "id": f"synthetic-capture-{position}",
            "kind": "running-capture",
            "scope": "committed-prefix",
            "producer": producer,
            "capture": capture,
            "snapshot": snapshot,
            "barrier": barrier,
            "index": reports.evidence(index),
        }
        return source, directory, generation, records

    def compose(self, *, paths=None, cache=None, proof_root=None):
        self.write(self.raw_path, self.raw)
        self.counter += 1
        return composer.compose_union(
            self.source_paths if paths is None else paths,
            self.raw_path,
            self.cache if cache is None else cache,
            self.root / f"report-{self.counter}",
            data_manifest=self.data_manifest,
            proof_root=self.proof_root if proof_root is None else proof_root,
        )

    def test_actual_census_supplies_null_expectations_and_preserves_originals(self):
        original = [path.read_bytes() for path in self.source_paths]
        value = self.compose()
        report = value["raw_report"]
        self.assertEqual(
            report["inventory"],
            {
                "compressed_bytes": 200,
                "stored_record_references": 17,
                "physical_chunks": 4,
            },
        )
        self.assertEqual(
            report["data_manifest_sha256"], reports.digest(self.data_manifest)
        )
        self.assertNotIn("data_manifest_sha256", self.raw["compiled_cohort"])
        self.assertEqual(
            value["sources"][0]["original_producer_cohort"]["data_manifest_sha256"],
            "d" * 64,
        )
        self.assertFalse(value["all_source_runs_completed"])
        self.assertFalse(value["physical_examples_complete"])
        self.assertFalse(value["performance_eligible"])
        self.assertEqual(value["completed_physical_examples"], [])
        self.assertEqual([path.read_bytes() for path in self.source_paths], original)
        self.assertEqual(
            value["validation"]["raw_receipt"]["native_union_validation"],
            reports.evidence(self.raw_path),
        )

    def test_source_ids_order_and_exact_snapshots_are_required(self):
        for mutation in (
            "omitted",
            "reordered",
            "wrong-index",
            "wrong-generation",
            "wrong-count",
            "wrong-blake3",
            "unknown-source",
        ):
            with self.subTest(mutation=mutation):
                original = copy.deepcopy(self.raw)
                paths = self.source_paths
                if mutation == "omitted":
                    paths = paths[:-1]
                elif mutation == "reordered":
                    paths = list(reversed(paths))
                elif mutation == "wrong-index":
                    self.raw["captured_snapshots"][1]["manifest_sha256"] = "f" * 64
                elif mutation == "wrong-generation":
                    self.raw["sources_newest_first"][1]["expected"]["generation"] += 1
                elif mutation == "wrong-count":
                    self.raw["source_censuses"][1]["stored_record_references"] += 1
                elif mutation == "wrong-blake3":
                    self.raw["captured_snapshots"][1]["manifest_blake3"] = [42] * 32
                else:
                    self.raw["sources_newest_first"][1]["id"] = "unknown"
                with self.assertRaises(ValueError):
                    self.compose(paths=paths)
                self.raw = original

    def test_completed_capture_forms_compose_without_physical_promotion(self):
        for position, kind in enumerate(("cf1", "monitor"), start=1):
            captured = capture_sources.CompletedCaptureSourceTests()
            captured.setUp()
            self.addCleanup(captured.doCleanups)
            captured.binding = copy.deepcopy(self.binding)
            captured.worker["cache_statistics"]["binding"] = copy.deepcopy(self.binding)
            captured.snapshot["binding"] = {
                **copy.deepcopy(self.binding),
                "artifact_digest": bytes(self.binding["artifact_digest"]).hex(),
            }
            captured.index.write_text(f"SYNTHETIC DISTINCT CLOSED {kind} INDEX")
            captured.snapshot["delta"]["manifest_bytes"] = captured.index.stat().st_size
            captured.snapshot["delta"]["manifest_digest"] = bytes([position] * 32).hex()
            captured.capture["snapshot"] = copy.deepcopy(captured.snapshot)
            source, _ = captured.arguments(kind)
            source["id"] = f"synthetic-closed-{kind}"
            self.write(self.source_paths[position], source)
            actual = self.raw["source_censuses"][position]
            actual.update(
                generation=captured.native["generation"],
                stored_record_references=captured.native["committed_records"],
                manifest_bytes=captured.index.stat().st_size,
            )
            self.raw["sources_newest_first"][position].update(
                id=source["id"],
                directory=str(captured.capture_directory),
                cache_native_index_identity=reports.native_cache_identity(
                    captured.capture_directory
                ),
                expected={
                    **actual,
                    "runs": None,
                    "chunk_references": None,
                    "referenced_compressed_bytes": None,
                },
            )
            frozen = self.native_root / f"snapshots/source-{position}.rrindex"
            frozen.write_bytes(captured.index.read_bytes())
            self.raw["captured_snapshots"][position].update(
                generation=actual["generation"],
                manifest_bytes=actual["manifest_bytes"],
                manifest_sha256=reports.digest(frozen),
            )
        total = sum(
            row["stored_record_references"] for row in self.raw["source_censuses"]
        )
        self.raw["source_reference_totals"]["stored_record_references"] = total
        self.raw["output_census"]["stored_record_references"] = total
        self.raw["copy_statistics"]["preserved_record_references"] = total
        value = self.compose()
        self.assertTrue(value["all_source_runs_completed"])
        self.assertEqual(value["completed_physical_examples"], [])
        self.assertFalse(value["physical_examples_complete"])
        self.assertFalse(value["performance_eligible"])
        self.assertEqual(total, 14)
        for source in value["sources"][1:]:
            original = source["receipts"]["producer"]["raw_receipt"]
            self.assertTrue(original["provisional"])
            self.assertFalse(original["native_keys_authenticated"])

    def test_bad_native_copy_arithmetic_or_inferred_nulls_refuse(self):
        mutations = (
            lambda r: r["copy_statistics"].update(copied_bases=2),
            lambda r: r["copy_statistics"].update(namespaces=True),
            lambda r: r["copy_statistics"].update(charged_disk_bytes=1),
            lambda r: r["copy_statistics"].update(unique_physical_chunks=7),
            lambda r: r["copy_statistics"].update(physical_compressed_bytes=301),
            lambda r: r["source_reference_totals"].update(stored_record_references=99),
            lambda r: r["output_census"].update(generation=2),
            lambda r: r["source_censuses"][1].update(runs=None),
            lambda r: r["sources_newest_first"][1]["expected"].pop("runs"),
            lambda r: r["sources_newest_first"][1]["expected"].update(
                manifest_bytes=True
            ),
            lambda r: r["captured_snapshots"][1].update(refresh_allowed=True),
        )
        for mutate in mutations:
            with self.subTest(mutation=mutate):
                original = copy.deepcopy(self.raw)
                mutate(self.raw)
                with self.assertRaises(ValueError):
                    self.compose()
                self.raw = original

    def test_actual_output_location_and_frozen_file_are_bound(self):
        with self.assertRaisesRegex(ValueError, "output directory"):
            self.compose(cache=self.finite.cache)
        (self.native_root / "snapshots/source-1.rrindex").write_text("CHANGED SNAPSHOT")
        with self.assertRaisesRegex(ValueError, "snapshot bytes"):
            self.compose()

    def test_own_artifact_admission_and_opaque_scope_are_required(self):
        mutations = (
            lambda r: r.update(success=False),
            lambda r: r.update(source_indexes_unchanged=False),
            lambda r: r.update(frozen_readers_used=False),
            lambda r: r.update(distinct_keys=17),
            lambda r: r.update(IBP_applications=True),
            lambda r: r.update(coefficient_context_created=True),
            lambda r: r["compiled_cohort"].update(main_revision="f" * 40),
            lambda r: r["resource_guard_admission"].update(admitted=False),
            lambda r: r["resource_guard_admission"]["worker"].update(
                executable_sha256="f" * 64
            ),
            lambda r: r["semantic_binding"].update(family_fingerprint="FOREIGN"),
        )
        for mutate in mutations:
            with self.subTest(mutation=mutate):
                original = copy.deepcopy(self.raw)
                mutate(self.raw)
                with self.assertRaises(ValueError):
                    self.compose()
                self.raw = original
        with self.assertRaisesRegex(ValueError, "regular file"):
            self.compose(proof_root=self.root)

    def test_native_output_cannot_alias_a_source(self):
        self.raw["sources_newest_first"][1]["directory"] = str(self.cache)
        self.raw["sources_newest_first"][1]["kind"] = "bundle"
        self.raw["sources_newest_first"][1]["cache_native_index_identity"] = self.raw[
            "cache_native_index_identity"
        ]
        with self.assertRaisesRegex(ValueError, "overlaps the Native output"):
            self.compose()

    def test_existing_source_validator_rejects_forged_completion(self):
        source = reports.read(self.source_paths[1])
        source["scope"] = "completed-capture"
        self.write(self.source_paths[1], source)
        with self.assertRaises((ValueError, KeyError)):
            self.compose()
        self.assertFalse((self.root / f"report-{self.counter}").exists())
        self.assertEqual(list(self.root.glob(".report-*")), [])

    def test_cli_requires_explicit_union_context(self):
        command = [
            sys.executable,
            str(Path(composer.__file__).resolve()),
            "--source",
            str(self.source_paths[0]),
            "--source",
            str(self.source_paths[1]),
            "--cache",
            str(self.cache),
            "--output-directory",
            str(self.root / "unused-cli"),
        ]
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 2)
        self.assertIn("require --union-validation", result.stderr)
        self.assertFalse((self.root / "unused-cli").exists())


if __name__ == "__main__":
    unittest.main()
