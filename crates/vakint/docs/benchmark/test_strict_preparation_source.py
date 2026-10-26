"""Synthetic provenance contracts only; no Native codec, imports or reductions."""

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import compose_cumulative_cache_report as composer
import merge_five_loop_results as reports


class StrictPreparationSourceTests(unittest.TestCase):
    def write(self, name, value):
        path = self.root / name
        path.write_text(json.dumps(value))
        return reports.evidence(path)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="synthetic-strict-source-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.cache = self.root / "opaque-native"
        self.cache.mkdir()
        index = self.cache / "cache.rrindex"
        index.write_text("SYNTHETIC OPAQUE INDEX; NOT A NATIVE CACHE")
        index.chmod(0o444)
        identity = reports.native_cache_identity(self.cache)
        binding = {
            "artifact_digest": [1] * 32,
            "family_fingerprint": "SYNTHETIC",
            "physical_arity": 15,
            "layer": "saved-application-route-apply-v1",
        }
        cohort = {
            field: str(i + 1) * (64 if field.endswith("_sha256") else 40)
            for i, (field, _) in enumerate(reports.COHORT_STAMPS)
        }
        artifact = {
            "stamps": {stamp: cohort[field] for field, stamp in reports.COHORT_STAMPS},
            "binary": {"sha256": "9" * 64},
        }
        original_cohort = {
            "main_revision": "a" * 40,
            "application_sha256": "b" * 64,
            "data_manifest_sha256": "d" * 64,
        }
        single = {
            "records": 4,
            "decoded_sha256": "c" * 64,
            "identity": {"test": True},
            "checkpoint_path": "/historical/original",
            "binding": binding,
        }
        producer = {
            "schema": "vakint.five-loop-cache-build.v2",
            "complete": True,
            "publishable": True,
            "data_manifest_sha256": "d" * 64,
            "cache_sha256": "e" * 64,
            "build_identity": {
                "data_manifest_sha256": "d" * 64,
                "build": {"git_revision": "a" * 40},
            },
            "preparation": {
                "status": "complete",
                "total_roots": 2,
                "completed_roots": 2,
                "native_profile": {"strict": True, "unresolved_policy": "reject"},
                "options": {"single_file": single},
            },
            "supervision": {"success": True, "hard_deadline_exhausted": False},
            "owner_proof": {"compiled_cohort": original_cohort},
        }
        self.source = {
            "id": "synthetic-finite-preparation",
            "kind": "strict-preparation",
            "scope": "completed-preparation",
        }
        self.source["producer"] = self.write("producer.json", producer)
        self.source["acceptance"] = self.write(
            "acceptance.json",
            {
                "schema": "vakint.strict-complete-small-candidate-acceptance.v1",
                "state": "STRICT_COMPLETE_CANDIDATE_NOT_FINAL_SELECTION",
                "reports": {"build": self.source["producer"]},
                "selected_training_roots": 2,
                "completed_strict_roots": 2,
                "record_count": 4,
                "cache_file": {
                    "sha256": "e" * 64,
                    "decoded_sha256": "c" * 64,
                    "native_identity": single["identity"],
                },
            },
        )
        preservation = {
            "schema": "vakint.accepted-small-native-cache-preservation.v1",
            "state": "PASS_STORAGE_PRESERVATION",
            "original_bundle": "/historical/original",
            "accepted_gzip": {"sha256": "e" * 64},
            "stored_record_references": 4,
            "generation": 3,
            "compressed_chunk_bytes": 100,
            "chunks": 2,
            "indices": [
                {
                    "relative_path": name,
                    "sha256": row["sha256"],
                    "destination": {"bytes": row["bytes"]},
                }
                for name, row in identity["indexes"].items()
            ],
        }
        preservation.update(
            {
                key: True
                for key in (
                    "all_source_and_destination_compressed_bytes_equal",
                    "all_gzip_members_crc_lengths_hashes_and_header_bindings_valid",
                    "all_original_sources_stable",
                    "all_output_inodes_independent",
                    "read_only_preserved_bundle",
                )
            }
        )
        self.source["preservation"] = self.write("preservation.json", preservation)
        self.source["inventory"] = self.write(
            "inventory.json", {"cache_native_index_identity": identity}
        )
        self.source["native_artifacts"] = self.write("artifact.json", artifact)
        self.worker = {"pid": 101, "start_ticks": 1001, "uid": 1124}
        native = {
            "schema": "vakint.strict-cache-source-parity-native-validation.v1",
            "success": True,
            "source_cache_native_index_identity": identity,
            "source_inventory_sha256": self.source["inventory"]["sha256"],
            "semantic_binding": binding,
            "application_sha256": "b" * 64,
            "legacy_cache_sha256": "e" * 64,
            "native_unique_records": 4,
            "legacy_validated_records": 4,
            "source_native_stored_record_references": 4,
            "source_native_generation": 3,
            "missing_keys": 0,
            "partial_records": 0,
            "exact_record_and_ordered_guard_equality": True,
            "declared_terminal_validation": True,
            "source_unchanged": True,
            "functional_IBP_comparison": False,
            "new_reduction_or_capture": False,
            "codec_state_format": 6,
            "collector_artifacts": self.source["native_artifacts"],
            "compiled_cohort": cohort,
            "executable_sha256": "9" * 64,
            "resource_guard_admission": {
                "schema": "vakint.resource-guard-admission.v1",
                "admitted": True,
                "launch_spec_sha256": "6" * 64,
                "guard_config_sha256": "7" * 64,
                "slot": "synthetic-source-parity",
                "worker": {**self.worker, "executable_sha256": "9" * 64},
            },
        }
        self.source["native_validation"] = self.write("native.json", native)
        self.source["native_supervision"] = self.write(
            "supervision.json",
            {
                "schema": "vakint.guarded-small-native-gzip-source-parity.v1",
                "accepted": True,
                "worker_exit_code": 0,
                "errors": [],
                "raw_validation": self.source["native_validation"],
                "worker": self.worker,
                "admission_receipt": native["resource_guard_admission"],
                "registration": {"signals_sent": []},
                "guard_after_worker_exit": {"ActiveState": "inactive"},
            },
        )

    def replace(self, name, mutate):
        path = Path(self.source[name]["path"])
        value = reports.read(path)
        mutate(value)
        self.source[name] = self.write(path.name, value)
        # Rebind only actual enclosing evidence, so tests exercise semantic checks.
        if name == "native_validation":
            self.replace(
                "native_supervision",
                lambda row: row.update(raw_validation=self.source[name]),
            )

    def test_separate_cohorts_and_completed_finite_scope(self):
        value = reports.strict_preparation_source(self.source, self.root)
        self.assertNotEqual(
            value["original_producer_cohort"]["main_revision"],
            value["current_validation_cohort"]["main_revision"],
        )
        self.assertEqual(value["scope"], "completed-preparation")
        self.assertEqual(value["inventory_counts"]["stored_record_references"], 4)

    def test_selected_subset_cannot_qualify_complete_source(self):
        self.replace(
            "native_validation",
            lambda row: row.update(
                schema="vakint.strict-light-cache-native-validation.v1",
                selected_records=2,
            ),
        )
        with self.assertRaisesRegex(ValueError, "full bound"):
            reports.strict_preparation_source(self.source, self.root)

    def test_incomplete_producer_refuses_without_final_selection_promotion(self):
        self.replace("producer", lambda row: row.update(complete=False))
        with self.assertRaisesRegex(ValueError, "completed strict"):
            reports.strict_preparation_source(self.source, self.root)

    def test_original_producer_data_and_revision_cannot_be_relabelled(self):
        for field, value in (
            ("data_manifest_sha256", "f" * 64),
            ("main_revision", "f" * 40),
        ):
            with self.subTest(field=field):
                original = reports.read(self.source["producer"]["path"])
                changed = copy.deepcopy(original)
                changed["owner_proof"]["compiled_cohort"][field] = value
                self.source["producer"] = self.write("changed-producer.json", changed)
                with self.assertRaisesRegex(ValueError, "completed strict"):
                    reports.strict_preparation_source(self.source, self.root)
                self.source["producer"] = reports.evidence(self.root / "producer.json")

    def test_false_ordered_guard_proof_and_bool_counts_refuse(self):
        for update in (
            {"exact_record_and_ordered_guard_equality": False},
            {"missing_keys": False},
        ):
            with self.subTest(update=update):
                before = copy.deepcopy(self.source)
                path = Path(self.source["native_validation"]["path"])
                data = path.read_bytes()
                self.replace("native_validation", lambda row: row.update(update))
                with self.assertRaisesRegex(ValueError, "full bound"):
                    reports.strict_preparation_source(self.source, self.root)
                path.write_bytes(data)
                self.source = before

    def test_different_source_index_binding_and_application_refuse(self):
        for key, value in (
            ("source_inventory_sha256", "0" * 64),
            ("application_sha256", "0" * 64),
            ("source_cache_native_index_identity", {}),
        ):
            with self.subTest(field=key):
                native = reports.read(self.source["native_validation"]["path"])
                changed = copy.deepcopy(native)
                changed[key] = value
                self.source["native_validation"] = self.write(
                    "changed-native.json", changed
                )
                with self.assertRaisesRegex(ValueError, "full bound"):
                    reports.strict_preparation_source(self.source, self.root)
                self.source["native_validation"] = reports.evidence(
                    self.root / "native.json"
                )

    def test_current_reader_cohort_must_match_its_own_artifact(self):
        self.replace(
            "native_validation",
            lambda row: row["compiled_cohort"].update(main_revision="f" * 40),
        )
        with self.assertRaisesRegex(ValueError, "artifact"):
            reports.strict_preparation_source(self.source, self.root)

    def test_supervised_worker_cannot_be_substituted(self):
        self.replace(
            "native_supervision", lambda row: row["worker"].update(start_ticks=1002)
        )
        with self.assertRaisesRegex(ValueError, "guarded validation"):
            reports.strict_preparation_source(self.source, self.root)

    def test_supervised_launch_guard_tuple_cannot_be_crossed(self):
        self.replace(
            "native_validation",
            lambda row: row["resource_guard_admission"].update(
                guard_config_sha256="0" * 64
            ),
        )
        with self.assertRaisesRegex(ValueError, "guarded validation"):
            reports.strict_preparation_source(self.source, self.root)

    def test_compose_actual_index_report_preserves_original_data_identity(self):
        source = self.root / "source.json"
        source.write_text(json.dumps(self.source))
        value = composer.compose(source, self.cache, self.root / "output")
        self.assertEqual(value["raw_report"]["data_manifest_sha256"], "d" * 64)
        self.assertTrue(value["all_source_runs_completed"])
        self.assertFalse(value["performance_eligible"])
        self.assertEqual(value["completed_physical_examples"], [])

    def test_direct_consumer_rejects_rebound_other_snapshot_or_original_data(self):
        source = self.root / "source.json"
        source.write_text(json.dumps(self.source))
        composer.compose(source, self.cache, self.root / "output")
        report_path = self.root / "output/cumulative-report.json"
        original = reports.read(report_path)
        for change in ("index", "inventory", "data"):
            report = copy.deepcopy(original)
            if change == "index":
                identity = report["cache_native_index_identity"]
                identity["indexes"]["cache.rrindex"]["sha256"] = "f" * 64
                body = {
                    key: value for key, value in identity.items() if key != "sha256"
                }
                identity["sha256"] = hashlib.sha256(
                    json.dumps(body, sort_keys=True, separators=(",", ":")).encode()
                ).hexdigest()
            elif change == "inventory":
                report["inventory"]["stored_record_references"] += 1
            else:
                report["data_manifest_sha256"] = "f" * 64
            validation_path = self.root / "output/native-validation.json"
            validation = reports.read(validation_path)
            for key in (
                "cache_native_index_identity",
                "inventory",
                "data_manifest_sha256",
            ):
                validation[key] = report[key]
            validation_path.write_text(json.dumps(validation))
            report["validation"] = reports.evidence(validation_path)
            report_path.write_text(json.dumps(report))
            with (
                self.subTest(change=change),
                self.assertRaisesRegex(ValueError, "one-source finite report"),
            ):
                reports.cumulative_cache_record(
                    report_path,
                    None,
                    source_identity=report["cache_native_index_identity"],
                )

    def test_changed_actual_index_refuses_before_output_creation(self):
        index = self.cache / "cache.rrindex"
        index.chmod(0o644)
        index.write_text("DIFFERENT SYNTHETIC INDEX")
        index.chmod(0o444)
        source = self.root / "source.json"
        source.write_text(json.dumps(self.source))
        with self.assertRaisesRegex(ValueError, "actual preserved Native indexes"):
            composer.compose(source, self.cache, self.root / "output")
        self.assertFalse((self.root / "output").exists())


if __name__ == "__main__":
    unittest.main()
