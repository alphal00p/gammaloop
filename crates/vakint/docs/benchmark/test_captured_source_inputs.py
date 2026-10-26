"""Synthetic captured-input provenance only; no Native codec or reductions."""

import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import merge_five_loop_results as reports


class CapturedSourceInputTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="synthetic-captured-inputs-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.final_data_manifest_sha256 = "a" * 64
        self.binding = {
            "artifact_digest": [1] * 32,
            "family_fingerprint": "SYNTHETIC-FAMILY",
            "physical_arity": 15,
            "layer": "saved-application-route-apply-v1",
        }
        self.cohort = {
            "main_revision": "1" * 40,
            "native_revision": "2" * 40,
            "kernel_revision": "3" * 40,
            "semantic_binding": self.binding,
        }
        self.tuning_keys = [[2] + [1] * 13 + [i] for i in range(1, 8)]
        self.q0_keys = [[3] + [1] * 13 + [i] for i in range(1, 4)]
        self.held_out = [9] * 15
        public_input = "SYNTHETIC public tuning input with α"
        self.workload = {
            "schema": "vakint.five-loop-cache-workload.v1",
            "data_manifest_sha256": "a" * 64,
            "source_census_sha256": "c" * 64,
            "training": [],
            "tuning": [
                {
                    "id": "synthetic-tuning",
                    "public_input": public_input,
                    "powers": [1] * 15,
                    "lowered_roots": self.tuning_keys,
                }
            ],
            "held_out": [],
        }
        workload = self.write("workload.json", self.workload)
        self.worker = {
            "schema": "vakint.five-loop-cache-measurement.v1",
            "success": True,
            "worker_success": True,
            "error": None,
            "case_id": "synthetic-tuning",
            "provisional": True,
            "native_keys_authenticated": False,
            "input_sha256": "b" * 64,
            "actual_lowered_roots": self.tuning_keys,
            "cache_statistics": {"binding": copy.deepcopy(self.binding)},
            "identity": {
                "workload_sha256": workload["sha256"],
                "data_manifest_sha256": "a" * 64,
                "source_census_sha256": "c" * 64,
                "build_git_revision": "f" * 40,
                "suite": "tuning",
                "only_case": "synthetic-tuning",
                "provisional": True,
            },
        }
        self.proof_row = {
            "id": "synthetic-tuning",
            "logical_cohort": "tuning",
            "public_input_source_sha256": hashlib.sha256(
                public_input.encode()
            ).hexdigest(),
            "input_sha256": self.worker["input_sha256"],
            "actual_lowered_roots": copy.deepcopy(self.tuning_keys),
            "application_binding": copy.deepcopy(self.binding),
            "reduction_succeeded": True,
            "complete_owner_transport": True,
            "mapping_error": None,
        }
        self.proof = {
            "schema": "vakint.public-lowering-owner-map-proof.v1",
            "compiled_cohort": copy.deepcopy(self.cohort),
            "semantic_binding": copy.deepcopy(self.binding),
            "rows": [self.proof_row],
        }
        scalar = self.root / "vakint-input.symbolica"
        scalar.write_bytes(b"SYNTHETIC scalar input; not a Native Atom")
        scalar_identity = {
            "sha256": reports.digest(scalar),
            "bytes": scalar.stat().st_size,
        }
        self.manifest = {
            "schema": "vakint.qcd-ladder.prepared-input.v1",
            "definition": "massless-qcd-four-rung-native-uv-v1",
            "projection": "zero-momentum-metric",
            "input": copy.deepcopy(scalar_identity),
        }
        self.config = {
            "schema": "vakint.inactive-capture-physical-user-unit.v1",
            "job_id": "q0",
            "projection": "zero-momentum",
            "prepared_input": {
                "directory": str(self.root),
                "manifest": self.write("prepared-input.json", self.manifest),
                "manifest_projection": "zero-momentum-metric",
                "input": copy.deepcopy(scalar_identity),
            },
        }
        self.census = {
            "schema": "vakint.qcd-ladder.campaign-census.v1",
            "source_sha256": scalar_identity["sha256"],
            "source_bytes": scalar_identity["bytes"],
            "projection": "zero-momentum-metric",
            "family_fingerprint": self.binding["family_fingerprint"],
            "combined_roots": [{"powers": key} for key in self.q0_keys],
        }
        self.inventory_rows = {
            "tuning": {
                "id": "captured-tuning",
                "kind": "captured-tuning-workload",
                "case_id": "synthetic-tuning",
            },
            "q0": {"id": "synthetic-physical-q0", "kind": "physical-census"},
        }
        self.sources = {
            name: {
                "id": row["id"],
                "kind": "running-capture",
                "scope": "committed-prefix",
            }
            for name, row in self.inventory_rows.items()
        }
        self.routes = {
            "schema": "vakint.source-candidate-owner-qualification.v2",
            "all_routes_succeeded": True,
            "selected_rows_complete": True,
            "compiled_cohort": copy.deepcopy(self.cohort),
            "semantic_binding": copy.deepcopy(self.binding),
            "rows": [],
        }
        self.aliases = {
            "schema": "vakint.source-candidate-external-protected-alias-exclusion.v1",
            "rows": [],
        }
        self.protected = {
            "schema": "vakint.protected-raw-owner-union.v1",
            "raw_lowered_keys": [self.held_out],
            "direct_owner_keys": [],
            "declared_case_powers": [self.held_out],
            "forbidden_direct_target_or_owner_keys": [self.held_out],
        }

    def write(self, name, value):
        path = self.root / name
        path.write_text(json.dumps(value, indent=2))
        return reports.evidence(path)

    def route_keys(self, *names):
        keys = {"tuning": self.tuning_keys, "q0": self.q0_keys}
        self.routes["rows"] = [
            {
                "raw_powers": key,
                "source_row_index": index,
                "census_root_index": index,
                "route_success": True,
                "route_error": None,
                "owner_endpoints": [{"coefficient": "1", "powers": [5] + key[1:]}],
            }
            for index, key in enumerate(key for name in names for key in keys[name])
        ]
        self.aliases["rows"] = [
            {
                "powers": row["raw_powers"],
                "route_complete_nonempty": True,
                "raw_protected_intersection": False,
                "owner_protected_intersections": [],
                "eligible_after_external_alias_exclusion": True,
            }
            for row in self.routes["rows"]
        ]

    def arguments(self, *names):
        # Rehash enclosing evidence after mutations so semantic checks are tested.
        workload = self.write("workload.json", self.workload)
        worker = self.write("worker.json", self.worker)
        self.config["prepared_input"]["manifest"] = self.write(
            "prepared-input.json", self.manifest
        )
        config = self.write("config.json", self.config)
        census = self.write("census.json", self.census)
        origins = {"tuning": worker, "q0": config}
        rows, sources = [], []
        for name in names:
            rows.append(
                {
                    **self.inventory_rows[name],
                    "producer": origins[name],
                    "inputs": workload if name == "tuning" else census,
                }
            )
            sources.append(
                {
                    **self.sources[name],
                    "producer": origins[name],
                    **(
                        {"config": config}
                        if name == "q0"
                        and self.sources[name]["scope"] == "completed-example"
                        else {}
                    ),
                }
            )
        proof = self.write("proof.json", self.proof)
        cohort = self.write("cohort.json", self.cohort)
        self.protected.update(owner_proof=proof, compiled_cohort=cohort)
        protected = self.write("protected.json", self.protected)
        routes = self.write("routes.json", self.routes)
        self.aliases.update(qualification=routes, protected_union=protected)
        inventory = self.write(
            "inventory.json",
            {"schema": "vakint.cumulative-source-starting-keys.v1", "sources": rows},
        )
        accepted = {
            "source_starting_keys": inventory,
            "consumer": {
                "owner_proof": proof,
                "compiled_cohort": cohort,
                "training_route_qualification": routes,
                "protected_union": protected,
            },
        }
        docs = {
            "compiled_cohort": self.cohort,
            "owner_proof": self.proof,
            "training_route_qualification": self.routes,
            "protected_union": self.protected,
            "alias_exclusion": self.aliases,
            "final_workload": {
                "data_manifest_sha256": self.final_data_manifest_sha256,
                "held_out": [
                    {
                        "powers": self.held_out,
                        "protected_roots": [self.held_out],
                        "lowered_roots": [self.held_out],
                    }
                ],
            },
        }
        producer = {
            "raw_report": {
                "sources": sources,
                "semantic_binding": self.binding,
                "data_manifest_sha256": self.final_data_manifest_sha256,
            }
        }
        return accepted, docs, producer, self.root

    def check(self, *names):
        reports.ReportEvidence.cumulative_source_exclusion(
            None, *self.arguments(*names)
        )

    def reject(self, *names):
        with self.assertRaises(ValueError):
            self.check(*names)

    def test_tuning_keeps_historical_flags_and_all_seven_raw_roots(self):
        self.route_keys("tuning")
        original = copy.deepcopy(self.worker)
        self.check("tuning")
        self.assertEqual(self.worker, original)
        self.assertTrue(self.worker["provisional"])
        self.assertFalse(self.worker["native_keys_authenticated"])
        self.assertNotEqual(
            self.worker["identity"]["build_git_revision"], self.cohort["main_revision"]
        )

    def test_physical_census_accepts_authenticated_pathless_prepared_input(self):
        self.route_keys("q0")
        self.assertNotIn("path", self.config["prepared_input"]["input"])
        self.check("q0")
        self.assertEqual(self.sources["q0"]["scope"], "committed-prefix")

    def test_mixed_sources_keep_the_complete_starting_union(self):
        self.route_keys("tuning", "q0")
        self.check("tuning", "q0")

    def test_mixed_source_preserves_original_snapshot_with_current_semantic_binding(
        self,
    ):
        self.route_keys("tuning", "q0")
        self.final_data_manifest_sha256 = "d" * 64
        self.cohort["data_manifest_sha256"] = self.final_data_manifest_sha256
        self.proof.update(
            data_manifest_sha256=self.final_data_manifest_sha256,
            compiled_cohort=copy.deepcopy(self.cohort),
        )
        self.routes["compiled_cohort"] = copy.deepcopy(self.cohort)
        original_worker, original_workload = copy.deepcopy((self.worker, self.workload))
        self.check("tuning", "q0")
        self.assertEqual(self.worker, original_worker)
        self.assertEqual(self.workload, original_workload)
        self.assertEqual(self.worker["identity"]["data_manifest_sha256"], "a" * 64)
        self.assertNotEqual(
            self.worker["identity"]["data_manifest_sha256"],
            self.final_data_manifest_sha256,
        )
        self.worker["identity"]["data_manifest_sha256"] = "e" * 64
        self.reject("tuning", "q0")
        self.worker["identity"]["data_manifest_sha256"] = "a" * 64
        self.worker["cache_statistics"]["binding"]["family_fingerprint"] = "FOREIGN"
        self.reject("tuning", "q0")

    def test_equal_malformed_original_source_snapshot_hashes_are_rejected(self):
        self.route_keys("tuning")
        for value in (None, "short", True, "A" * 64):
            with self.subTest(source_data_manifest=value):
                self.workload["data_manifest_sha256"] = value
                self.worker["identity"]["data_manifest_sha256"] = value
                self.worker["identity"]["workload_sha256"] = self.write(
                    "workload.json", self.workload
                )["sha256"]
                self.reject("tuning")

    def test_historical_native_census_remains_an_unknown_source_kind(self):
        self.route_keys("q0")
        self.inventory_rows["q0"]["kind"] = "q0-native-census"
        self.census.update(
            schema="vakint.qcd-ladder.q0-native-census.v1",
            family_fingerprint=self.binding["family_fingerprint"],
            projection="zero-momentum-metric",
            ibp_application=False,
            route_search_complete=False,
        )
        with self.assertRaisesRegex(ValueError, "unknown starting-key source"):
            self.check("q0")

    def test_tuning_case_status_and_normalized_input_must_match(self):
        for field, value in (
            ("schema", "unknown"),
            ("case_id", "different-case"),
            ("success", False),
            ("success", 1),
            ("error", "original failure"),
            ("failure", "original failure"),
            ("resource_failure", "original resource refusal"),
            ("input_sha256", "0" * 64),
            ("provisional", False),
            ("provisional", "true"),
            ("native_keys_authenticated", True),
            ("native_keys_authenticated", "false"),
        ):
            with self.subTest(field=field, value=value):
                original = copy.deepcopy(self.worker)
                self.worker[field] = value
                self.route_keys("tuning")
                self.reject("tuning")
                self.worker = original

    def test_tuning_workload_hash_and_case_are_bound(self):
        self.route_keys("tuning")
        self.worker["identity"]["workload_sha256"] = "0" * 64
        self.reject("tuning")
        self.worker["identity"]["workload_sha256"] = self.write(
            "workload.json", self.workload
        )["sha256"]
        self.inventory_rows["tuning"]["case_id"] = "different-case"
        self.reject("tuning")

    def test_tuning_original_identity_cannot_be_relabelled(self):
        self.route_keys("tuning")
        for field, value in (
            ("suite", "held-out"),
            ("only_case", "other-case"),
            ("provisional", False),
            ("provisional", 1),
            ("source_census_sha256", "0" * 64),
            ("data_manifest_sha256", "0" * 64),
        ):
            with self.subTest(field=field):
                original = self.worker["identity"][field]
                self.worker["identity"][field] = value
                self.reject("tuning")
                self.worker["identity"][field] = original
        self.worker["cache_statistics"]["binding"] = {}
        self.reject("tuning")

    def test_optional_worker_and_supervision_status_must_not_contradict_success(self):
        self.route_keys("tuning")
        self.worker.pop("worker_success")
        self.check("tuning")
        supervision = {"success": True, "deadline_exhausted": False, "exit_code": 0}
        self.worker.update(worker_success=True, supervision=supervision)
        self.check("tuning")
        for flag in (False, 1, "true", None):
            with self.subTest(worker_success=flag):
                self.worker["worker_success"] = flag
                self.reject("tuning")
        self.worker["worker_success"] = True
        for wrapper in (
            None,
            [],
            {**supervision, "success": False},
            {**supervision, "success": 1},
            {**supervision, "deadline_exhausted": True},
            {**supervision, "deadline_exhausted": 0},
            {**supervision, "exit_code": 1},
            {**supervision, "exit_code": False},
            {**supervision, "exit_code": 0.0},
            {"success": True, "deadline_exhausted": False},
        ):
            with self.subTest(supervision=wrapper):
                self.worker["supervision"] = wrapper
                self.reject("tuning")

    def test_raw_roots_cannot_be_replaced_by_a_combined_subset(self):
        self.route_keys("tuning")
        for target in (self.worker, self.proof_row):
            original = target["actual_lowered_roots"]
            for replacement in (self.tuning_keys[:5], list(reversed(self.tuning_keys))):
                target["actual_lowered_roots"] = replacement
                self.reject("tuning")
            target["actual_lowered_roots"] = original

    def test_workload_raw_root_changes_cannot_hide_behind_new_file_hashes(self):
        self.route_keys("tuning")
        for keys in (self.tuning_keys[:5], list(reversed(self.tuning_keys))):
            self.workload["tuning"][0]["lowered_roots"] = keys
            self.worker["identity"]["workload_sha256"] = self.write(
                "workload.json", self.workload
            )["sha256"]
            self.reject("tuning")

    def test_tuning_case_and_current_proof_cannot_be_missing_or_ambiguous(self):
        self.route_keys("tuning")
        for rows in ([], [self.proof_row, copy.deepcopy(self.proof_row)]):
            self.proof["rows"] = rows
            self.reject("tuning")
        self.proof["rows"] = [self.proof_row]
        case = self.workload["tuning"][0]
        self.workload["tuning"].append(copy.deepcopy(case))
        self.worker["identity"]["workload_sha256"] = self.write(
            "workload.json", self.workload
        )["sha256"]
        self.reject("tuning")

    def test_current_tuning_proof_requires_exact_source_transport_and_cohort(self):
        self.route_keys("tuning")
        for field, value in (
            ("logical_cohort", "held-out"),
            ("public_input_source_sha256", "0" * 64),
            ("complete_owner_transport", False),
            ("complete_owner_transport", 1),
            ("reduction_succeeded", False),
            ("mapping_error", "transport failed"),
            ("application_binding", {}),
        ):
            with self.subTest(field=field):
                original = self.proof_row[field]
                self.proof_row[field] = value
                self.reject("tuning")
                self.proof_row[field] = original
        self.proof["compiled_cohort"]["main_revision"] = "e" * 40
        self.reject("tuning")

    def test_physical_census_projection_family_and_input_are_bound(self):
        self.route_keys("q0")
        for field, value in (
            ("schema", "unknown"),
            ("projection", "quadratic-transverse"),
            ("family_fingerprint", "DIFFERENT"),
            ("source_sha256", "0" * 64),
            ("source_bytes", 1),
            ("source_bytes", True),
        ):
            with self.subTest(field=field):
                original = self.census[field]
                self.census[field] = value
                self.reject("q0")
                self.census[field] = original

    def test_completed_source_resolves_its_own_hashed_config(self):
        self.route_keys("q0")
        self.sources["q0"]["scope"] = "completed-example"
        self.check("q0")

    def test_physical_config_manifest_and_scalar_cannot_be_substituted(self):
        self.route_keys("q0")
        for target in (self.config["prepared_input"]["input"], self.manifest["input"]):
            for field, value in (("sha256", "0" * 64), ("bytes", 1), ("bytes", True)):
                with self.subTest(field=field):
                    original = target[field]
                    target[field] = value
                    self.reject("q0")
                    target[field] = original
        (self.root / "vakint-input.symbolica").write_text("DIFFERENT INPUT")
        self.reject("q0")

    def test_physical_manifest_definition_projection_and_directory_are_bound(self):
        self.route_keys("q0")
        for field, value in (
            ("schema", "unknown"),
            ("definition", "DIFFERENT"),
            ("projection", "quadratic-transverse"),
        ):
            with self.subTest(field=field):
                original = self.manifest[field]
                self.manifest[field] = value
                self.reject("q0")
                self.manifest[field] = original
        self.config["prepared_input"]["directory"] = "relative-input"
        self.reject("q0")

    def test_physical_root_vectors_require_fifteen_actual_integers(self):
        self.route_keys("q0")
        for key in ([1] * 14, [True] + [1] * 14, [1.0] + [1] * 14, ["1"] + [1] * 14):
            with self.subTest(key=key):
                original = self.census["combined_roots"]
                self.census["combined_roots"] = [{"powers": key}, *original[1:]]
                self.reject("q0")
                self.census["combined_roots"] = original

    def test_changed_hashed_source_evidence_is_rejected(self):
        for name, filename in (
            ("tuning", "workload.json"),
            ("tuning", "worker.json"),
            ("q0", "census.json"),
            ("q0", "prepared-input.json"),
        ):
            with self.subTest(source=name, filename=filename):
                self.route_keys(name)
                arguments = self.arguments(name)
                (self.root / filename).write_text("{}")
                with self.assertRaisesRegex(ValueError, "changed evidence"):
                    reports.ReportEvidence.cumulative_source_exclusion(None, *arguments)

    def test_inventory_must_name_every_actual_source_and_its_producer(self):
        self.route_keys("tuning", "q0")
        for mutation in ("omitted-source", "other-producer"):
            with self.subTest(mutation=mutation):
                accepted, docs, producer, base = self.arguments("tuning", "q0")
                inventory = reports.read(accepted["source_starting_keys"]["path"])
                if mutation == "omitted-source":
                    inventory["sources"].pop()
                else:
                    inventory["sources"][0]["producer"]["sha256"] = "0" * 64
                accepted["source_starting_keys"] = self.write(
                    "inventory.json", inventory
                )
                with self.assertRaises(ValueError):
                    reports.ReportEvidence.cumulative_source_exclusion(
                        None, accepted, docs, producer, base
                    )

    def test_missing_route_or_alias_cannot_hide_an_authenticated_starting_key(self):
        for name in ("tuning", "q0"):
            for document in (self.routes, self.aliases):
                with self.subTest(source=name, document=document["schema"]):
                    self.route_keys(name)
                    document["rows"].pop()
                    self.reject(name)

    def test_protected_raw_keys_and_nonzero_owners_are_rejected(self):
        for name in ("tuning", "q0"):
            for owner in (False, True):
                with self.subTest(source=name, owner=owner):
                    self.route_keys(name)
                    row = self.routes["rows"][0]
                    protected = (
                        row["owner_endpoints"][0]["powers"]
                        if owner
                        else row["raw_powers"]
                    )
                    self.protected["raw_lowered_keys"] = [self.held_out, protected]
                    self.protected["forbidden_direct_target_or_owner_keys"] = [
                        self.held_out,
                        protected,
                    ]
                    self.reject(name)

    def test_route_and_alias_success_flags_do_not_override_failure(self):
        for name in ("tuning", "q0"):
            for field, value in (
                ("route_complete_nonempty", False),
                ("raw_protected_intersection", True),
                ("owner_protected_intersections", [[7] * 15]),
                ("eligible_after_external_alias_exclusion", False),
            ):
                with self.subTest(source=name, field=field):
                    self.route_keys(name)
                    self.aliases["rows"][0][field] = value
                    self.reject(name)
        self.route_keys("tuning")
        self.routes["compiled_cohort"]["main_revision"] = "e" * 40
        self.reject("tuning")
        self.routes["compiled_cohort"] = copy.deepcopy(self.cohort)
        self.routes["semantic_binding"] = {}
        self.reject("tuning")


if __name__ == "__main__":
    unittest.main()
