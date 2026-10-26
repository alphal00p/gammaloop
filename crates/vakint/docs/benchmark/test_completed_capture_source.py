"""Synthetic closure receipts and opaque index bytes only; no Native cache I/O."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

import merge_five_loop_results as reports


class CompletedCaptureSourceTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="synthetic-capture-closure-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.capture_directory = self.root / "capture"
        self.capture_directory.mkdir()
        self.index = self.capture_directory / "cache.rrindex"
        self.index.write_bytes(b"SYNTHETIC OPAQUE INDEX; NOT A NATIVE CACHE")
        self.value = self.root / "value.symbolica"
        self.value.write_bytes(b"SYNTHETIC COMPLETE VALUE")
        self.guards = self.root / "guards.json"
        self.guards.write_text('["SYNTHETIC GUARD A", "SYNTHETIC GUARD B"]')
        self.binding = {
            "artifact_digest": [171] * 32,
            "family_fingerprint": "SYNTHETIC-FAMILY",
            "physical_arity": 15,
            "layer": "saved-application-route-apply-v1",
        }
        self.native = {
            "schema": "rustred.running-application-cache.v1",
            "generation": 7,
            "accepted_records": 5,
            "committed_records": 5,
            "pending_records": 0,
            "pending_bytes": 0,
            "closed": True,
            "writer_error": None,
            "last_error": None,
        }
        self.snapshot = {
            "schema": "rustred.running-application-cache-snapshot.v1",
            "binding": {
                **copy.deepcopy(self.binding),
                "artifact_digest": bytes(self.binding["artifact_digest"]).hex(),
            },
            "delta": {
                "generation": 7,
                "manifest_digest": "SYNTHETIC-NATIVE-MANIFEST-DIGEST",
                "manifest_path": str(self.index),
                "manifest_bytes": self.index.stat().st_size,
            },
        }
        self.capture = {
            "path": str(self.capture_directory),
            "terminal": True,
            "flush_complete": True,
            "flush_error": None,
            "snapshot_error": None,
            "reduction_error": None,
            "report": self.native,
            "flush_report": {**self.native, "closed": False},
            "close": {"error": None, "report": copy.deepcopy(self.native)},
            "snapshot": copy.deepcopy(self.snapshot),
            "snapshot_descriptor": str(self.root / "snapshot.json"),
            "snapshot_sidecar": str(self.index),
        }
        self.workload = self.write(
            "workload.json",
            {
                "schema": "vakint.five-loop-cache-workload.v1",
                "data_manifest_sha256": "a" * 64,
                "source_census_sha256": "c" * 64,
                "training": [],
                "tuning": [{"id": "synthetic-tuning", "public_input": "SYNTHETIC"}],
                "held_out": [],
            },
        )
        self.worker = {
            "schema": "vakint.five-loop-cache-measurement.v1",
            "success": True,
            "provisional": True,
            "native_keys_authenticated": False,
            "case_id": "synthetic-tuning",
            "input_sha256": "b" * 64,
            "value_sha256": reports.digest(self.value),
            "guards_sha256": reports.digest(self.guards),
            "cache_statistics": {"binding": copy.deepcopy(self.binding)},
            "identity": {
                "suite": "tuning",
                "only_case": "synthetic-tuning",
                "provisional": True,
                "data_manifest_sha256": "a" * 64,
                "workload_sha256": self.workload["sha256"],
                "build_git_revision": "f" * 40,
            },
        }
        self.terminal = {
            "schema": "vakint.large-native-cache.tuning-pilot.v1",
            "state": "PASS",
            "worker_returncode": 0,
            "case_id": "synthetic-tuning",
            "value_bytes_equal": True,
            "complete_sorted_unique_guard_set_equal": True,
        }
        self.scope = {
            "provisional": True,
            "native_keys_authenticated": False,
            "held_out": False,
            "performance_qualification": False,
            "final_cache_complete": False,
        }

    def write(self, name, document):
        path = self.root / name
        path.write_text(json.dumps(document, indent=2))
        return reports.evidence(path)

    def arguments(
        self, kind="cf1", *, producer=None, completion=None, result=None, retained=None
    ):
        # Changes precede enclosing hashes, so negative cases test semantic joins.
        self.worker["running_cache"] = copy.deepcopy(self.capture)
        if producer:
            producer(self.worker)
        refs = {
            "producer": self.write("worker.json", self.worker),
            "capture": self.write("capture.json", self.capture),
            "snapshot": self.write("snapshot.json", self.snapshot),
            "index": reports.evidence(self.index),
        }
        documents = {
            "producer": self.worker,
            "capture": self.capture,
            "snapshot": self.snapshot,
        }
        if kind == "cf1":
            terminal = {**self.terminal, "worker_receipt": refs["producer"]}
            terminal_ref = self.write("terminal.json", terminal)
            closure = {
                "schema": "vakint.large-native-cache.functional-tuning-receipt.v1",
                "status": "PASS_EXPOSED_TUNING_ONLY",
                "actual_worker": refs["producer"],
                "actual_terminal": terminal_ref,
                "workload": self.workload,
                "new_closed_capture": {
                    "path": self.capture["path"],
                    **self.native,
                    "error": None,
                    "snapshot_descriptor": refs["snapshot"],
                    "cache_index": {
                        **refs["index"],
                        "bytes": self.index.stat().st_size,
                    },
                },
                "functional_parity": {
                    "value_sha256": self.worker["value_sha256"],
                    "complete_exported_guard_set_sha256": self.worker["guards_sha256"],
                },
                "scope": copy.deepcopy(self.scope),
            }
        else:
            evidence = {
                "native_report": copy.deepcopy(self.native),
                "flush_report": copy.deepcopy(self.capture["flush_report"]),
                "snapshot": copy.deepcopy(self.snapshot),
                "retained_files": {
                    ref["path"]: {
                        "sha256": ref["sha256"],
                        "bytes": Path(ref["path"]).stat().st_size,
                    }
                    for name, ref in refs.items()
                    if name in ("capture", "snapshot", "index")
                },
            }
            if retained:
                retained(evidence)
            refs["capture_evidence"] = self.write("capture-evidence.json", evidence)
            observed = {
                "status": "PASS",
                "current_identity": copy.deepcopy(self.worker["identity"]),
                **{
                    key: self.worker[key]
                    for key in ("input_sha256", "value_sha256", "guards_sha256")
                },
                "running_capture_evidence_sha256": refs["capture_evidence"]["sha256"],
                "running_accepted_records": self.native["accepted_records"],
                "running_committed_records": self.native["committed_records"],
                "performance_claim": False,
                "full_value_and_guards_byte_parity": True,
                "runs": [
                    {
                        "mode": "monitored-q",
                        "returncode": 0,
                        "timed_out": False,
                        "errors": [],
                    }
                ],
            }
            if result:
                result(observed)
            refs["result"] = self.write("result.json", observed)
            closure = {
                "status": "PASS_ACTUAL_ONE_REGULAR_CLI_Q_AND_RUNNING_CAPTURE",
                "result_sha256": refs["result"]["sha256"],
                "capture_evidence_sha256": refs["capture_evidence"]["sha256"],
                "unit_terminal_properties": {
                    "ActiveState": "inactive",
                    "SubState": "dead",
                    "Result": "success",
                    "ExecMainStatus": "0",
                },
                "owned_processes_absent": {"collector": True, "worker": True},
            }
        if completion:
            completion(closure)
        refs["completion"] = self.write("completion.json", closure)
        barrier = {
            "schema": "vakint.native-cache-snapshot-barrier.v1",
            "accepted": True,
            "committed_snapshot": True,
            "producer_completed": True,
            "writer_fault": None,
            "manifest_digest": self.snapshot["delta"]["manifest_digest"],
            **{
                key: self.native[key]
                for key in (
                    "generation",
                    "accepted_records",
                    "committed_records",
                    "pending_records",
                    "pending_bytes",
                )
            },
            **{
                key: refs[key]
                for key in ("producer", "capture", "snapshot", "completion")
            },
        }
        refs["barrier"] = self.write("barrier.json", barrier)
        documents["barrier"] = barrier
        source = {
            "id": "synthetic-tuning-capture",
            "kind": "running-capture",
            "scope": "completed-capture",
            **refs,
        }
        docs = {
            name: {
                "evidence": ref,
                **({"raw_receipt": documents[name]} if name in documents else {}),
            }
            for name, ref in refs.items()
            if name in ("producer", "capture", "snapshot", "index", "barrier")
        }
        return source, docs

    def check(self, kind="cf1", **changes):
        source, docs = self.arguments(kind, **changes)
        return reports.completed_capture_source(
            source, docs, self.root, semantic_binding=self.binding
        )

    def reject(self, kind="cf1", **changes):
        with self.assertRaises(ValueError):
            self.check(kind, **changes)

    def test_both_actual_receipt_forms_complete_only_the_capture(self):
        for kind in ("cf1", "monitor"):
            with self.subTest(kind=kind):
                accepted = self.check(kind)
                self.assertEqual(accepted["scope"], "completed-capture")
                self.assertEqual(
                    accepted["original_producer_identity"], self.worker["identity"]
                )
                self.assertTrue(
                    accepted["receipts"]["producer"]["raw_receipt"]["provisional"]
                )
                self.assertFalse(
                    accepted["receipts"]["producer"]["raw_receipt"][
                        "native_keys_authenticated"
                    ]
                )
                self.assertIn("no physical example", accepted["scope_note"])
                self.assertIn("completion", accepted["receipts"])

    def test_completed_capture_does_not_complete_physical_examples_or_performance(self):
        for kind in ("cf1", "monitor"):
            with self.subTest(kind=kind):
                source, _ = self.arguments(kind)
                report = {
                    "schema": reports.CUMULATIVE_SCHEMA,
                    "cache_sha256": None,
                    "cache_native_index_identity": reports.native_cache_identity(
                        self.capture_directory
                    ),
                    "semantic_binding": self.binding,
                    "data_manifest_sha256": "a" * 64,
                    "sources": [source],
                    "inventory": {
                        "compressed_bytes": 100,
                        "stored_record_references": 5,
                        "physical_chunks": 1,
                    },
                }
                validation = {
                    **report,
                    "schema": "vakint.cumulative-native-cache-validation.v1",
                    "accepted": True,
                    "checks": {
                        key: True
                        for key in (
                            "immutable_index_snapshot",
                            "authenticated_chunks",
                            "strict_completed_records_only",
                            "authenticated_bindings",
                            "committed_source_snapshots",
                        )
                    },
                }
                report["validation"] = self.write("validation.json", validation)
                reference = self.write("cumulative.json", report)
                accepted = reports.cumulative_cache_record(
                    reference["path"], self.capture_directory
                )
                self.assertTrue(accepted["all_source_runs_completed"])
                self.assertEqual(accepted["completed_physical_examples"], [])
                self.assertFalse(accepted["physical_examples_complete"])
                self.assertFalse(accepted["performance_eligible"])

    def test_source_scope_and_original_capture_must_not_be_substituted(self):
        for field, value in (
            ("kind", "strict-preparation"),
            ("scope", "completed-example"),
            ("scope", "committed-prefix"),
        ):
            with self.subTest(field=field):
                source, docs = self.arguments()
                source[field] = value
                with self.assertRaisesRegex(ValueError, "its own source scope"):
                    reports.completed_capture_source(
                        source, docs, self.root, semantic_binding=self.binding
                    )
        self.reject(
            producer=lambda row: row["running_cache"].update(path="SYNTHETIC OTHER RUN")
        )
        self.reject(completion=lambda row: row.update(status="FAIL"))
        self.reject("monitor", completion=lambda row: row.update(status="FAIL"))

    def test_common_flush_close_errors_and_terminal_flags_refuse(self):
        for field, value in (
            ("terminal", False),
            ("flush_complete", False),
            ("flush_error", "failed"),
            ("snapshot_error", "failed"),
            ("reduction_error", "failed"),
        ):
            with self.subTest(field=field):
                original = self.capture[field]
                self.capture[field] = value
                self.reject()
                self.capture[field] = original
        self.capture["close"]["error"] = "close failed"
        self.reject("monitor")

    def test_native_closed_counters_and_errors_are_not_summary_flags(self):
        for field, value in (
            ("closed", False),
            ("accepted_records", 6),
            ("committed_records", True),
            ("pending_records", 1),
            ("pending_bytes", 1),
            ("generation", 0),
            ("writer_error", "failed"),
            ("last_error", "failed"),
        ):
            with self.subTest(field=field):
                original = self.native[field]
                self.native[field] = value
                self.capture["close"]["report"] = copy.deepcopy(self.native)
                self.capture["flush_report"] = {**self.native, "closed": False}
                self.reject()
                self.native[field] = original
        self.native["accepted_records"] = self.native["committed_records"] = 0
        self.capture["close"]["report"] = copy.deepcopy(self.native)
        self.capture["flush_report"] = {**self.native, "closed": False}
        self.reject("monitor")

    def test_flush_close_and_snapshot_must_belong_to_this_capture(self):
        for target, field, value in (
            (self.capture["close"]["report"], "generation", 8),
            (self.capture["flush_report"], "committed_records", 4),
            (self.capture["snapshot"]["delta"], "manifest_digest", "OTHER"),
            (self.capture, "snapshot_descriptor", str(self.root / "other.json")),
            (self.capture, "snapshot_sidecar", str(self.root / "other.rrindex")),
        ):
            with self.subTest(field=field):
                original = target[field]
                target[field] = value
                self.reject()
                target[field] = original

    def test_snapshot_binding_uses_exact_native_hex_digest_and_semantics(self):
        original = copy.deepcopy(self.snapshot["binding"])
        for field, value in (
            ("family_fingerprint", "OTHER-FAMILY"),
            ("artifact_digest", "0" * 64),
            ("artifact_digest", original["artifact_digest"].upper()),
            ("artifact_digest", self.binding["artifact_digest"]),
            ("physical_arity", 14),
            ("physical_arity", 15.0),
            ("layer", "OTHER-LAYER"),
        ):
            with self.subTest(field=field, value=value):
                invalid = {**original, field: value}
                self.snapshot["binding"] = invalid
                self.capture["snapshot"]["binding"] = copy.deepcopy(invalid)
                self.reject()
                self.reject("monitor")
        self.snapshot.pop("binding")
        self.capture["snapshot"].pop("binding")
        self.reject()

    def test_failed_or_promoted_original_measurement_refuses(self):
        for field, value in (
            ("success", False),
            ("provisional", False),
            ("native_keys_authenticated", True),
            ("error", "failed"),
        ):
            with self.subTest(field=field):
                original = self.worker.get(field)
                self.worker[field] = value
                self.reject()
                self.worker[field] = original
        self.worker["cache_statistics"]["binding"] = {}
        self.reject("monitor")

    def test_full_value_and_guards_are_rehashed(self):
        for path in (self.value, self.guards):
            with self.subTest(path=path.name):
                original = path.read_bytes()
                path.write_bytes(b"CHANGED FULL PAYLOAD")
                self.reject()
                self.reject("monitor")
                path.write_bytes(original)

    def test_cf1_requires_successful_terminal_and_exact_worker(self):
        for field, value in (
            ("state", "FAIL"),
            ("worker_returncode", 1),
            ("worker_returncode", False),
            ("case_id", "OTHER"),
            ("value_bytes_equal", False),
            ("complete_sorted_unique_guard_set_equal", False),
        ):
            with self.subTest(field=field):
                original = self.terminal[field]
                self.terminal[field] = value
                self.reject()
                self.terminal[field] = original
        other = self.write("other-worker.json", {"synthetic": "another run"})
        self.reject(completion=lambda row: row.update(actual_worker=other))
        foreign_worker = self.write("foreign-worker.json", copy.deepcopy(self.worker))
        foreign_terminal = self.write(
            "foreign-terminal.json",
            {**self.terminal, "worker_receipt": foreign_worker},
        )
        self.reject(completion=lambda row: row.update(actual_terminal=foreign_terminal))

    def test_cf1_cannot_promote_scope_or_change_parity(self):
        for field in (
            "native_keys_authenticated",
            "held_out",
            "performance_qualification",
            "final_cache_complete",
        ):
            with self.subTest(field=field):
                self.reject(
                    completion=lambda row, field=field: row["scope"].update(
                        {field: True}
                    )
                )
        self.reject(completion=lambda row: row["scope"].update(provisional=False))
        for field in ("value_sha256", "complete_exported_guard_set_sha256"):
            self.reject(
                completion=lambda row, field=field: row["functional_parity"].update(
                    {field: "0" * 64}
                )
            )

    def test_cf1_closed_snapshot_and_workload_cannot_be_substituted(self):
        workload = reports.read(self.workload["path"])
        workload["tuning"][0]["public_input"] = "SYNTHETIC OTHER INPUT"
        other = self.write("other-workload.json", workload)
        self.reject(completion=lambda row: row.update(workload=other))
        other_snapshot = self.write("other-snapshot.json", copy.deepcopy(self.snapshot))
        self.reject(
            completion=lambda row: row["new_closed_capture"].update(
                snapshot_descriptor=other_snapshot
            )
        )
        for field, value in (
            ("closed", False),
            ("generation", 8),
            ("committed_records", 4),
            ("pending_records", 1),
            ("error", "failed"),
            ("path", str(self.root / "other")),
        ):
            with self.subTest(field=field):
                self.reject(
                    completion=lambda row, field=field, value=value: row[
                        "new_closed_capture"
                    ].update({field: value})
                )
        for field, value in (
            ("sha256", "0" * 64),
            ("bytes", True),
            ("path", str(self.root / "other" / "cache.rrindex")),
        ):
            self.reject(
                completion=lambda row, field=field, value=value: row[
                    "new_closed_capture"
                ]["cache_index"].update({field: value})
            )

    def test_monitor_raw_result_identity_counts_and_run_must_succeed(self):
        for field, value in (
            ("status", "FAIL"),
            ("current_identity", {}),
            ("input_sha256", "0" * 64),
            ("value_sha256", "0" * 64),
            ("guards_sha256", "0" * 64),
            ("running_accepted_records", 6),
            ("running_committed_records", True),
            ("performance_claim", True),
            ("full_value_and_guards_byte_parity", False),
        ):
            with self.subTest(field=field):
                self.reject(
                    "monitor",
                    result=lambda row, field=field, value=value: row.update(
                        {field: value}
                    ),
                )
        for field, value in (
            ("mode", "unmonitored"),
            ("returncode", 1),
            ("returncode", False),
            ("timed_out", True),
            ("errors", ["failed"]),
        ):
            self.reject(
                "monitor",
                result=lambda row, field=field, value=value: row["runs"][0].update(
                    {field: value}
                ),
            )
        self.reject("monitor", result=lambda row: row.update(runs=[]))

    def test_monitor_completion_and_retained_files_must_match_actual_source(self):
        for field in ("result_sha256", "capture_evidence_sha256"):
            self.reject(
                "monitor",
                completion=lambda row, field=field: row.update({field: "0" * 64}),
            )
        for field in ("native_report", "flush_report", "snapshot"):
            self.reject(
                "monitor", retained=lambda row, field=field: row.update({field: {}})
            )
        self.reject("monitor", retained=lambda row: row.update(retained_files={}))
        for filename in ("capture.json", "snapshot.json", "capture/cache.rrindex"):
            self.reject(
                "monitor",
                retained=lambda row, filename=filename: row["retained_files"].pop(
                    str(self.root / filename)
                ),
            )
        self.reject(
            "monitor",
            retained=lambda row: row["retained_files"][str(self.index)].update(
                bytes=True
            ),
        )
        self.reject(
            "monitor",
            retained=lambda row: row["retained_files"][str(self.index)].update(
                sha256="0" * 64
            ),
        )

    def test_monitor_process_and_unit_outcome_cannot_be_overridden(self):
        for field, value in (
            ("ActiveState", "active"),
            ("SubState", "running"),
            ("Result", "failed"),
            ("ExecMainStatus", "1"),
        ):
            self.reject(
                "monitor",
                completion=lambda row, field=field, value=value: row[
                    "unit_terminal_properties"
                ].update({field: value}),
            )
        for absent in (
            {},
            {"collector": True, "worker": False},
            {"collector": True, "worker": 1},
        ):
            self.reject(
                "monitor",
                completion=lambda row, absent=absent: row.update(
                    owned_processes_absent=absent
                ),
            )

    def test_completion_hash_and_barrier_cannot_name_another_receipt(self):
        source, docs = self.arguments()
        other = self.write("other-completion.json", {"status": "SYNTHETIC OTHER"})
        docs["barrier"]["raw_receipt"]["completion"] = other
        with self.assertRaisesRegex(ValueError, "another capture completion"):
            reports.completed_capture_source(
                source, docs, self.root, semantic_binding=self.binding
            )
        Path(source["completion"]["path"]).write_text("{}")
        with self.assertRaisesRegex(ValueError, "changed evidence"):
            reports.completed_capture_source(
                source, docs, self.root, semantic_binding=self.binding
            )


if __name__ == "__main__":
    unittest.main()
