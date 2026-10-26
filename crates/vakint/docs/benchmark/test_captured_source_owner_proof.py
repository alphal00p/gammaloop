"""Synthetic source/control proof joins only; no Native or binary reads."""

import copy
import unittest

import merge_five_loop_results as reports
import test_captured_source_inputs as captured_inputs


class CapturedSourceOwnerProofTests(unittest.TestCase):
    write = captured_inputs.CapturedSourceInputTests.write
    route_keys = captured_inputs.CapturedSourceInputTests.route_keys

    def setUp(self):
        captured_inputs.CapturedSourceInputTests.setUp(self)
        self.route_keys("tuning")
        self.cohort.update(
            main_source_snapshot_sha256="4" * 64,
            native_source_snapshot_sha256="5" * 64,
            data_manifest_sha256="a" * 64,
            loop_manifest_sha256="d" * 64,
            application_sha256="e" * 64,
            mapper_executable_sha256="6" * 64,
            measurement_executable_sha256="7" * 64,
        )
        headers = {
            key: self.cohort[key]
            for key in (
                "data_manifest_sha256",
                "loop_manifest_sha256",
                "application_sha256",
                "mapper_executable_sha256",
                "measurement_executable_sha256",
            )
        }
        headers.update(
            kernel_revision=self.cohort["kernel_revision"],
            rustred_revision=self.cohort["native_revision"],
            mapper_git_revision=self.cohort["main_revision"],
        )
        self.source_row = copy.deepcopy(self.proof_row)
        self.source_row["source_census_sha256"] = "c" * 64
        self.source_proof = {
            **copy.deepcopy(self.proof),
            **headers,
            "compiled_cohort": copy.deepcopy(self.cohort),
            "rows": [self.source_row],
        }
        self.source_cohort = copy.deepcopy(self.cohort)
        self.qualification_workload = {
            **copy.deepcopy(self.workload),
            "held_out": [{"id": "old-qualification-control", "powers": [4] * 15}],
        }
        own_workload = self.write(
            "source-qualification-workload.json", self.qualification_workload
        )
        self.qualification_worker = {
            **copy.deepcopy(self.worker),
            "fresh_process": True,
            "variant": "uncached",
            "sample": 0,
            "supervision": {
                "success": True,
                "exit_code": 0,
                "deadline_exhausted": False,
            },
            "identity": {
                **copy.deepcopy(self.worker["identity"]),
                "workload_sha256": own_workload["sha256"],
                "build_git_revision": self.cohort["main_revision"],
                "executable_sha256": self.cohort["measurement_executable_sha256"],
            },
        }
        new_row = {
            **copy.deepcopy(self.proof_row),
            "id": "new-unseen-tuning",
            "public_input_source_sha256": "8" * 64,
            "input_sha256": "9" * 64,
            "actual_lowered_roots": [[4] + [1] * 14],
        }
        self.proof.update(
            **headers,
            compiled_cohort=copy.deepcopy(self.cohort),
            rows=[new_row],
        )
        self.routes["compiled_cohort"] = copy.deepcopy(self.cohort)
        self.artifact = {
            "stamps": {
                stamp: self.cohort[field] for field, stamp in reports.COHORT_STAMPS
            },
            "artifacts": [
                {"name": "authenticate_owner_roots", "sha256": "6" * 64},
                {"name": "measure_cache", "sha256": "7" * 64},
            ],
        }
        self.use_source_proof = True

    def arguments(self):
        accepted, docs, producer, base = (
            captured_inputs.CapturedSourceInputTests.arguments(self, "tuning")
        )
        own_cohort = self.write("source-cohort.json", self.source_cohort)
        own_workload = self.write(
            "source-qualification-workload.json", self.qualification_workload
        )
        own_worker = self.write(
            "source-qualification-worker.json", self.qualification_worker
        )
        self.source_row.update(
            receipt_path=own_worker["path"],
            receipt_sha256=own_worker["sha256"],
            workload_path=own_workload["path"],
            workload_sha256=own_workload["sha256"],
        )
        self.source_proof.update(
            cohort_receipt_path=own_cohort["path"],
            cohort_receipt_sha256=own_cohort["sha256"],
        )
        source_proof = self.write("source-owner-proof.json", self.source_proof)
        inventory = reports.read(accepted["source_starting_keys"]["path"])
        if self.use_source_proof:
            inventory["sources"][0]["source_owner_proof"] = source_proof
        accepted["source_starting_keys"] = self.write("inventory.json", inventory)
        accepted["consumer"]["collector_artifacts"] = self.write(
            "current-artifacts.json", self.artifact
        )
        docs["collector_artifacts"] = self.artifact
        docs["final_workload"]["tuning"] = [
            {"id": "new-unseen-tuning", "public_input": "SYNTHETIC NEW INPUT"}
        ]
        return accepted, docs, producer, base

    def check(self):
        reports.ReportEvidence.cumulative_source_exclusion(None, *self.arguments())

    def reject(self):
        with self.assertRaises(ValueError):
            self.check()

    def rebind_qualification_workload(self):
        self.qualification_worker["identity"]["workload_sha256"] = self.write(
            "source-qualification-workload.json", self.qualification_workload
        )["sha256"]

    def test_separate_source_proof_preserves_old_input_without_current_eligibility(
        self,
    ):
        original = copy.deepcopy((self.worker, self.workload, self.proof))
        self.check()
        self.assertEqual((self.worker, self.workload, self.proof), original)
        self.assertEqual(
            [row["id"] for row in self.proof["rows"]], ["new-unseen-tuning"]
        )
        self.assertNotEqual(
            self.worker["identity"]["workload_sha256"],
            self.qualification_worker["identity"]["workload_sha256"],
        )
        self.assertNotEqual(
            self.worker["identity"]["build_git_revision"],
            self.qualification_worker["identity"]["build_git_revision"],
        )
        self.assertTrue(self.worker["provisional"])
        self.assertFalse(self.worker["native_keys_authenticated"])

    def test_no_source_reference_cannot_recover_an_ineligible_old_case(self):
        self.use_source_proof = False
        self.reject()

    def test_source_identity_is_not_rebased_to_current_metadata_headers(self):
        # This tests the input join, not the separate metadata-requalification gate.
        self.proof["data_manifest_sha256"] = "9" * 64
        self.proof["loop_manifest_sha256"] = "8" * 64
        self.final_data_manifest_sha256 = "9" * 64
        self.check()
        self.assertEqual(self.source_proof["data_manifest_sha256"], "a" * 64)
        self.assertEqual(self.source_proof["loop_manifest_sha256"], "d" * 64)
        self.assertEqual(self.worker["identity"]["data_manifest_sha256"], "a" * 64)

    def test_source_proof_hash_and_own_cohort_hash_are_authenticated(self):
        for filename in ("source-owner-proof.json", "source-cohort.json"):
            with self.subTest(filename=filename):
                arguments = self.arguments()
                (self.root / filename).write_text("{}")
                with self.assertRaisesRegex(ValueError, "changed evidence"):
                    reports.ReportEvidence.cumulative_source_exclusion(None, *arguments)

    def test_source_proof_case_binding_and_raw_roots_are_exact(self):
        for field, value in (
            ("id", "other-case"),
            ("logical_cohort", "held-out"),
            ("public_input_source_sha256", "0" * 64),
            ("input_sha256", "0" * 64),
            ("actual_lowered_roots", self.tuning_keys[:5]),
            ("actual_lowered_roots", list(reversed(self.tuning_keys))),
            ("application_binding", {}),
            ("complete_owner_transport", False),
            ("reduction_succeeded", False),
            ("mapping_error", "failed"),
            ("source_census_sha256", "0" * 64),
        ):
            with self.subTest(field=field):
                original = self.source_row[field]
                self.source_row[field] = value
                self.reject()
                self.source_row[field] = original

    def test_source_header_and_full_own_cohort_cannot_be_relabelled(self):
        for field in (
            "data_manifest_sha256",
            "loop_manifest_sha256",
            "application_sha256",
            "mapper_executable_sha256",
            "measurement_executable_sha256",
            "kernel_revision",
            "rustred_revision",
            "mapper_git_revision",
        ):
            with self.subTest(field=field):
                original = self.source_proof[field]
                self.source_proof[field] = "0" * len(original)
                self.reject()
                self.source_proof[field] = original
        self.source_proof["compiled_cohort"]["main_source_snapshot_sha256"] = "0" * 64
        self.reject()

    def test_same_binding_does_not_authorize_a_different_proof_cohort(self):
        self.source_cohort["main_revision"] = "f" * 40
        self.source_proof["compiled_cohort"] = copy.deepcopy(self.source_cohort)
        self.source_proof["mapper_git_revision"] = "f" * 40
        self.reject()

    def test_source_proof_requires_both_current_authenticated_binaries(self):
        original = self.artifact["artifacts"]
        for artifacts in (original[:1], original[1:]):
            with self.subTest(artifacts=artifacts):
                self.artifact["artifacts"] = artifacts
                self.reject()
        self.artifact["artifacts"] = original

    def test_source_qualification_worker_is_complete_and_bound_to_its_own_proof(self):
        for field, value in (
            ("schema", "unknown"),
            ("case_id", "other-case"),
            ("success", False),
            ("worker_success", False),
            ("fresh_process", False),
            ("variant", "cached"),
            ("sample", True),
            ("input_sha256", "0" * 64),
            ("actual_lowered_roots", self.tuning_keys[:5]),
        ):
            with self.subTest(field=field):
                original = self.qualification_worker[field]
                self.qualification_worker[field] = value
                self.reject()
                self.qualification_worker[field] = original
        for field, value in (
            ("data_manifest_sha256", "0" * 64),
            ("build_git_revision", "f" * 40),
            ("executable_sha256", "0" * 64),
            ("workload_sha256", "0" * 64),
            ("source_census_sha256", "0" * 64),
        ):
            with self.subTest(identity_field=field):
                original = self.qualification_worker["identity"][field]
                self.qualification_worker["identity"][field] = value
                self.reject()
                self.qualification_worker["identity"][field] = original

    def test_own_qualification_workload_is_independent_but_exact_for_the_source(self):
        for field, value in (
            ("public_input", "SYNTHETIC DIFFERENT INPUT"),
            ("lowered_roots", self.tuning_keys[:5]),
            ("id", "other-case"),
        ):
            with self.subTest(field=field):
                case = self.qualification_workload["tuning"][0]
                original = case[field]
                case[field] = value
                self.rebind_qualification_workload()
                self.reject()
                case[field] = original

    def test_own_workload_schema_data_and_census_are_bound(self):
        for field, value in (
            ("schema", "unknown"),
            ("data_manifest_sha256", "0" * 64),
            ("source_census_sha256", "0" * 64),
        ):
            with self.subTest(field=field):
                original = self.qualification_workload[field]
                self.qualification_workload[field] = value
                self.rebind_qualification_workload()
                self.reject()
                self.qualification_workload[field] = original

    def test_current_routes_and_exclusions_still_cover_every_source_raw_key(self):
        self.routes["rows"].pop()
        self.reject()
        self.route_keys("tuning")
        self.aliases["rows"].pop()
        self.reject()
        self.route_keys("tuning")
        self.routes["compiled_cohort"]["main_revision"] = "f" * 40
        self.reject()

    def test_source_proof_cannot_override_current_protected_raw_or_owner_keys(self):
        for owner in (False, True):
            with self.subTest(owner=owner):
                self.route_keys("tuning")
                row = self.routes["rows"][0]
                forbidden = (
                    row["owner_endpoints"][0]["powers"] if owner else row["raw_powers"]
                )
                self.protected["raw_lowered_keys"] = [self.held_out, forbidden]
                self.protected["forbidden_direct_target_or_owner_keys"] = [
                    self.held_out,
                    forbidden,
                ]
                self.reject()


if __name__ == "__main__":
    unittest.main()
