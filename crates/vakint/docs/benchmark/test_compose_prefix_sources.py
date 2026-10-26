"""Synthetic historical prefix proofs and opaque bytes; no Native cache reads."""

import copy
import unittest
from pathlib import Path

import merge_five_loop_results as reports
import test_compose_union_sources as union_sources


class ComposePrefixSourcesTests(unittest.TestCase):
    def setUp(self):
        self.base = union_sources.ComposeUnionSourcesTests()
        self.base.setUp()
        self.addCleanup(self.base.doCleanups)
        self.root = self.base.root
        self.raw = self.base.raw
        self.sources = [reports.read(path) for path in self.base.source_paths]
        self.old_producer = {
            "schema": "vakint.five-loop-cache-build.v2",
            "complete": False,
            "publishable": False,
            "data_manifest_sha256": "e" * 64,
            "preparation": {
                "status": "failed",
                "completed_roots": 1,
                "total_roots": 3,
                "native_profile": {"strict": True, "unresolved_policy": "reject"},
            },
        }
        self.q0_producer = {
            "schema": "vakint.inactive-capture-physical-user-unit.v1",
            "running_cache_directory": self.raw["sources_newest_first"][2]["directory"],
        }
        rows = []
        self.summary = {
            "schema": "vakint.native-cache-real-union-proof.v1",
            "native_union_publication_success": True,
            # Original publication preceded the independent accepted barrier.
            "status": "PENDING_INDEPENDENT_INVENTORY",
            "binding": copy.deepcopy(self.base.binding),
        }
        for position, role in ((1, "old"), (2, "q0")):
            index = self.sources[position]["index"]
            frozen = self.raw["captured_snapshots"][position]
            census = self.raw["source_censuses"][position]
            rows.append(
                {
                    "role": role,
                    "path": index["path"],
                    "bytes": frozen["manifest_bytes"],
                    "native_blake3": bytes(frozen["manifest_blake3"]).hex(),
                    "identity": index["sha256"],
                }
            )
            self.summary[f"{role}_census"] = {
                **copy.deepcopy(census),
                "diagnostics_enabled": False,
                "extra_historical_diagnostic": {"source": role, "seconds": 0.5},
            }
            self.summary[f"{role}_snapshot"] = {
                "chunk_root": self.raw["sources_newest_first"][position]["directory"],
                "generation": census["generation"],
                "manifest_bytes": census["manifest_bytes"],
                "manifest_digest": copy.deepcopy(frozen["manifest_blake3"]),
                "manifest_digest_algorithm": "blake3",
                "manifest_path": index["path"],
                "refresh_allowed": False,
            }
            self.sources[position] = {
                "id": self.sources[position]["id"],
                "kind": "native-index-prefix",
                "scope": "committed-prefix",
                "prefix_role": role,
                "snapshot": copy.deepcopy(index),
            }
        self.metadata = {
            "immutable_indices": rows,
            "metadata_only": True,
            "run_concat_verified": True,
            # These describe the whole historical union, not either input prefix.
            "generation": 1,
            "compressed_bytes": 200,
            "preserved_record_references": 13,
        }
        self.terminal = {
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
        }
        self.publish_history()

    def write(self, name, value):
        return self.base.write(self.root / name, value)

    def publish_history(self, *, include_summary=True, embed_metadata=True):
        old = self.write("historical/failed-producer.json", self.old_producer)
        q0 = self.write("historical/q0-config.json", self.q0_producer)
        binding = self.write("historical/binding.json", {"original_build_receipt": old})
        proposal = self.write(
            "historical/proposal.json",
            {
                "binding": binding,
                "inputs_newest_first": [
                    {"role": "q0", "path": self.q0_producer["running_cache_directory"]}
                ],
            },
        )
        metadata = self.write("historical/metadata-proof.json", self.metadata)
        terminal = self.write("historical/terminal-outcome.json", self.terminal)
        if embed_metadata:
            self.summary["native_structural_metadata_proof"] = copy.deepcopy(
                self.metadata
            )
        summary = self.write("historical/native-union-summary.json", self.summary)
        artifacts = [
            metadata,
            terminal,
            self.sources[1]["snapshot"],
            self.sources[2]["snapshot"],
        ]
        if include_summary:
            artifacts.append(summary)
        self.proof = {
            "schema": "vakint.postflight-real-union-prefix-proof-qualified.v1",
            "data_preservation_state": "PASS",
            "original_runtime_receipts_unchanged": True,
            "proposal": proposal,
            "output": str(self.root / "historical-native-output"),
            "generation": self.metadata["generation"],
            "compressed_bytes": self.metadata["compressed_bytes"],
            "stored_record_references": self.metadata["preserved_record_references"],
            "physical_chunks": 4,
            "artifacts": artifacts,
        }
        barrier = self.write("historical/barrier.json", self.proof)
        for position, producer in ((1, old), (2, q0)):
            self.sources[position].update(
                producer=producer, capture=metadata, barrier=barrier
            )
        self.write_sources()

    def write_sources(self):
        for path, source in zip(self.base.source_paths, self.sources):
            self.base.write(path, source)

    def recalculate_raw_totals(self):
        for key in (
            "runs",
            "chunk_references",
            "stored_record_references",
            "referenced_compressed_bytes",
        ):
            total = sum(row[key] for row in self.raw["source_censuses"])
            self.raw["source_reference_totals"][key] = total
            self.raw["output_census"][key] = total
        copied = self.raw["copy_statistics"]
        copied["copied_bases"] = len(self.raw["sources_newest_first"])
        copied["preserved_record_references"] = self.raw["source_reference_totals"][
            "stored_record_references"
        ]

    def reject(self, *, paths=None, message=None):
        with self.assertRaisesRegex(ValueError, message or ".*"):
            self.base.compose(paths=paths)
        self.assertFalse((self.root / f"report-{self.base.counter}").exists())
        self.assertEqual(list(self.root.glob(".report-*")), [])

    def test_old_and_q0_prefixes_compose_with_failed_original_unchanged(self):
        originals = [path.read_bytes() for path in self.base.source_paths]
        failed_path = self.sources[1]["producer"]["path"]
        failed_bytes = Path(failed_path).read_bytes()
        value = self.base.compose()
        self.assertFalse(value["all_source_runs_completed"])
        self.assertFalse(value["physical_examples_complete"])
        self.assertFalse(value["performance_eligible"])
        self.assertEqual(value["completed_physical_examples"], [])
        self.assertEqual(
            value["raw_report"]["inventory"]["stored_record_references"], 17
        )
        self.assertEqual(
            [path.read_bytes() for path in self.base.source_paths], originals
        )
        self.assertEqual(Path(failed_path).read_bytes(), failed_bytes)
        self.assertFalse(
            value["sources"][1]["receipts"]["producer"]["raw_receipt"]["complete"]
        )
        for source, role in zip(value["sources"][1:], ("old", "q0")):
            self.assertEqual(source["scope"], "committed-prefix")
            native = source["receipts"]["native_summary"]
            self.assertEqual(
                native["evidence"],
                reports.evidence(self.root / "historical/native-union-summary.json"),
            )
            self.assertEqual(
                native["raw_receipt"][f"{role}_census"]["generation"],
                self.summary[f"{role}_snapshot"]["generation"],
            )
            self.assertNotEqual(
                self.metadata["generation"],
                self.summary[f"{role}_census"]["generation"],
            )

    def test_missing_q0_cannot_be_excused_by_new_native_union_success(self):
        for key in ("sources_newest_first", "captured_snapshots", "source_censuses"):
            self.raw[key] = self.raw[key][:2]
        self.recalculate_raw_totals()
        self.reject(
            paths=self.base.source_paths[:2], message="every captured source prefix"
        )

    def test_larger_later_q0_capture_does_not_prove_historical_prefix_containment(self):
        later, directory, generation, records = self.base.capture_source(3)
        later_path = self.root / "later-full-q0-source.json"
        self.base.write(later_path, later)
        original = self.summary["q0_census"]
        self.assertGreater(generation, original["generation"])
        self.assertGreater(records, original["stored_record_references"])
        index = directory / "cache.rrindex"
        census = {
            "generation": generation,
            "runs": 1,
            "chunk_references": 2,
            "stored_record_references": records,
            "referenced_compressed_bytes": 100,
            "manifest_bytes": index.stat().st_size,
        }
        self.raw["source_censuses"][2] = census
        self.raw["sources_newest_first"][2] = {
            "id": later["id"],
            "kind": "pack",
            "directory": str(directory),
            "cache_native_index_identity": reports.native_cache_identity(directory),
            "expected": copy.deepcopy(census),
        }
        frozen = self.base.native_root / "snapshots/source-2.rrindex"
        frozen.write_bytes(index.read_bytes())
        self.raw["captured_snapshots"][2].update(
            generation=generation,
            manifest_bytes=census["manifest_bytes"],
            manifest_blake3=[3] * 32,
            manifest_sha256=reports.digest(frozen),
        )
        self.recalculate_raw_totals()
        self.reject(
            paths=[*self.base.source_paths[:2], later_path],
            message="every captured source prefix",
        )

    def test_duplicate_or_foreign_roles_are_rejected(self):
        for role in ("old", "output", "invented-prefix"):
            with self.subTest(role=role):
                self.sources[2]["prefix_role"] = role
                self.write_sources()
                self.reject()
        self.sources[2]["prefix_role"] = "q0"
        self.write_sources()
        foreign = {**self.proof, "foreign_proof_marker": True}
        self.sources[2]["barrier"] = self.write(
            "historical/foreign-barrier.json", foreign
        )
        self.write_sources()
        self.reject(message="every captured source prefix")

    def test_prefix_cannot_be_relabelled_completed(self):
        for scope in (
            "completed-preparation",
            "completed-example",
            "completed-capture",
        ):
            with self.subTest(scope=scope):
                self.sources[1]["scope"] = scope
                self.write_sources()
                self.reject()

    def test_failed_whole_job_and_strict_policy_remain_authoritative(self):
        original = copy.deepcopy(self.old_producer)
        for mutate in (
            lambda p: p.update(complete=True),
            lambda p: p["preparation"].update(status="complete"),
            lambda p: p["preparation"]["native_profile"].update(strict=False),
            lambda p: p["preparation"]["native_profile"].update(
                unresolved_policy="retain"
            ),
        ):
            with self.subTest(mutation=mutate):
                self.old_producer = copy.deepcopy(original)
                mutate(self.old_producer)
                self.publish_history()
                self.reject()

    def test_actual_native_role_census_must_match_authenticated_historical_projection(
        self,
    ):
        original = copy.deepcopy(self.raw)
        for field in (
            "runs",
            "chunk_references",
            "stored_record_references",
            "referenced_compressed_bytes",
        ):
            with self.subTest(field=field):
                self.base.raw = self.raw = copy.deepcopy(original)
                self.raw["source_censuses"][1][field] += 1
                expected = self.raw["sources_newest_first"][1]["expected"]
                if expected[field] is not None:
                    expected[field] += 1
                self.recalculate_raw_totals()
                self.reject()

    def test_snapshot_sha_bytes_and_native_digest_cannot_be_substituted(self):
        original = copy.deepcopy(self.raw)
        for field, value in (
            ("manifest_sha256", "f" * 64),
            ("manifest_bytes", 1),
            ("manifest_blake3", [99] * 32),
        ):
            with self.subTest(field=field):
                self.base.raw = self.raw = copy.deepcopy(original)
                self.raw["captured_snapshots"][1][field] = value
                self.reject()

    def test_summary_and_nested_metadata_are_actual_authenticated_artifacts(self):
        self.publish_history(include_summary=False)
        self.reject()
        self.publish_history()
        self.summary["native_structural_metadata_proof"]["metadata_only"] = False
        self.publish_history(embed_metadata=False)
        self.reject()
        self.summary["native_union_publication_success"] = False
        self.publish_history()
        self.reject()

    def test_historical_census_and_snapshot_metadata_are_typed_and_bound(self):
        original = copy.deepcopy(self.summary)
        for mutate in (
            lambda s: s["old_census"].update(stored_record_references=True),
            lambda s: s["old_census"].update(runs=1.0),
            lambda s: s["old_snapshot"].update(refresh_allowed=True),
            lambda s: s["old_snapshot"].update(manifest_digest_algorithm="sha256"),
            lambda s: s["old_snapshot"].update(manifest_bytes=1),
            lambda s: s["old_snapshot"].update(manifest_digest=[98] * 32),
            lambda s: s["binding"].update(family_fingerprint="FOREIGN"),
        ):
            with self.subTest(mutation=mutate):
                self.summary = copy.deepcopy(original)
                mutate(self.summary)
                self.publish_history()
                self.reject()

    def test_prefix_source_cannot_alias_native_output(self):
        source = self.raw["sources_newest_first"][1]
        source.update(
            directory=str(self.base.cache),
            kind="bundle",
            cache_native_index_identity=self.raw["cache_native_index_identity"],
        )
        self.reject(message="overlaps the Native output")


if __name__ == "__main__":
    unittest.main()
