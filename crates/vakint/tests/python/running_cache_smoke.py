"""Running-cache Python bindings on one coherently built native host.

VAKINT_TEST_EXTERNAL_DATA selects the authenticated local campaign. One engine
loads it once, and the stored P65 graph supplies the existing five-tadpole input
with powers [1,1,1,1,5]. No graph generation or numerical master evaluation runs.
All cache outputs belong to this test's temporary directory.
"""

import hashlib
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

if os.environ.get("VAKINT_TEST_HOST") == "standalone":
    import vakint_hepkit_standalone as core
    from vakint_hepkit_standalone import (
        CacheFileIdentity,
        CacheFileOptions,
        E,
        MasterSymbol,
        RunningCacheControl,
        Vakint,
        VakintEvaluationMethod,
        integral_from_diagram,
    )
else:
    from symbolica import E, core
    from symbolica.community.hepkit.vakint import (
        CacheFileIdentity,
        CacheFileOptions,
        MasterSymbol,
        RunningCacheControl,
        Vakint,
        VakintEvaluationMethod,
        integral_from_diagram,
    )


def file_hashes(directory):
    return {
        str(path.relative_to(directory)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in directory.rglob("*")
        if path.is_file()
    }


class RunningCacheSmoke(unittest.TestCase):
    def test_closed_capture_exports_one_native_file_with_lazy_shared_loading(self):
        self.obj.detach_cache(5)
        reference = self.obj.evaluate_rustred_symbolic(self.input)
        delta = self.directory / "single-file-delta"
        output = self.directory / "single-file.rrcache.gz"
        control = self.obj.attach_running_cache(5, delta, read_backing=False)
        self.addCleanup(self.close_control, control)
        self.addCleanup(self.obj.detach_cache, 5)
        observed = self.obj.evaluate_rustred_symbolic(self.input)
        self.assert_complete_native_result(observed, reference)
        with self.assertRaisesRegex(ValueError, "explicitly closed"):
            self.obj.save_cache(5, str(output), mode="full-snapshot")
        self.assertFalse(output.exists())
        self.assert_durable(control.close(), closed=True)
        original = file_hashes(delta)
        report = self.obj.save_cache(5, str(output), mode="full-snapshot")
        self.assertTrue(report["complete"])
        self.assertGreater(report["single_file"]["statistics"]["records"], 0)
        identity = report["single_file"]["identity"]
        decoded = self.directory / "single-file-decoded"
        decoded.mkdir()
        options = CacheFileOptions(
            max_compressed_bytes=identity["compressed_bytes"],
            max_decoded_bytes=identity["decoded_bytes"],
            decoded_directory=str(decoded),
        )
        self.obj.load_cache(5, str(output), file_options=options)
        trusted = self.obj.cache_statistics(5)["saved_application"]["file"]
        self.assertEqual(trusted["state_imports"], 0)
        self.assertEqual(trusted["decoded_coefficients"], 0)
        typed_identity = CacheFileIdentity.from_document(identity)
        self.assertEqual(typed_identity.to_document(), identity)
        self.obj.load_cache(
            5, str(output), file_options=options, file_identity=typed_identity
        )
        before = self.obj.cache_statistics(5)["saved_application"]["file"]
        self.assertEqual(before["state_imports"], 0)
        self.assertEqual(before["decoded_coefficients"], 0)
        for field in (
            "compressed_bytes",
            "decoded_bytes",
            "compressed_blake3",
            "decoded_blake3",
        ):
            invalid = typed_identity.to_document()
            if field.endswith("bytes"):
                invalid[field] -= 1
            else:
                invalid[field][0] ^= 1
            with self.subTest(field=field):
                with self.assertRaisesRegex(ValueError, "mismatch"):
                    self.obj.load_cache(
                        5,
                        str(output),
                        file_options=options,
                        file_identity=CacheFileIdentity.from_document(invalid),
                    )
                self.assertEqual(
                    self.obj.cache_statistics(5)["saved_application"]["file"], before
                )
        cached = self.obj.evaluate_rustred_symbolic(self.input)
        self.assert_complete_native_result(cached, reference)
        after = self.assert_backing_hit()["saved_application"]["file"]
        self.assertEqual(after["state_imports"], 1)
        self.assertGreater(after["decoded_coefficients"], 0)
        self.assertEqual(file_hashes(delta), original)
        self.assertEqual(list(decoded.iterdir()), [])

    @classmethod
    def setUpClass(cls):
        if sys.flags.optimize or not core.is_licensed():
            raise RuntimeError(
                "running-cache tests need normal Python and a licensed host"
            )
        data = Path(os.environ["VAKINT_TEST_EXTERNAL_DATA"])
        cls.scratch = tempfile.TemporaryDirectory(prefix="vakint-running-python-")
        cls.addClassCleanup(cls.scratch.cleanup)
        cls.directory = Path(cls.scratch.name)
        cls.obj = Vakint(
            evaluation_order=[
                VakintEvaluationMethod.new_rustred_method(master_mode="retain")
            ],
            epsilon_symbol=E("ep", default_namespace="vakint"),
            integral_normalization_factor="FMFTandMATAD",
            use_dot_product_notation=True,
            form_exe_path="/running-cache-python-test-must-not-use-form",
            temporary_directory=str(cls.directory / "native"),
        )
        loaded = cls.obj.load_official_external_data(
            5, local_artifacts_path=data, load_prepared_caches=False
        )
        assert loaded["load_prepared_caches"] is False
        assert all(
            not cls.obj.cache_statistics(loop)["attached"] for loop in range(1, 6)
        )
        manifest = json.loads((data / "5L/manifest.json").read_text())
        metadata = next(
            master for master in manifest["masters"] if master["name"] == "P65"
        )
        assert metadata["powers"] == [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1]
        assert metadata["graph_definition"]["graph_edge_to_family_slot"] == [
            [0, 4],
            [1, 7],
            [2, 8],
            [3, 11],
            [4, 14],
        ]
        master = cls.obj.get_master(MasterSymbol(metadata["master_id"]))
        assert master.loop_count == 5 and not master.is_computed
        graph = master.graph
        assert graph.loop_count == 5
        assert graph.numerator_expression().to_expression() == E("1")
        propagators = graph.propagator_family()
        family = graph.integral_family(kinematics=propagators.kinematics)
        assert family.denominators[:5] == propagators.denominators
        powers = [1, 1, 1, 1, 5] + [0] * (len(family.denominators) - 5)
        cls.input = integral_from_diagram(
            graph, family, E("1"), powers=powers
        ).to_expression()

    def assert_report(self, report, *, closed=False):
        self.assertEqual(report["schema"], "rustred.running-application-cache.v1")
        self.assertEqual(report["writer_mode"], "thread")
        self.assertIs(report["closed"], closed)
        self.assertIsNone(report["last_error"])
        self.assertIsNone(report["writer_error"])
        for key in (
            "accepted_records",
            "committed_records",
            "pending_records",
            "pending_bytes",
            "generation",
        ):
            self.assertIs(type(report[key]), int)
            self.assertGreaterEqual(report[key], 0)

    def assert_durable(self, report, *, closed=False):
        self.assert_report(report, closed=closed)
        self.assertEqual(report["accepted_records"], report["committed_records"])
        self.assertEqual(report["pending_records"], 0)
        self.assertEqual(report["pending_bytes"], 0)

    def assert_snapshot(self, snapshot, sidecar, generation):
        self.assertEqual(
            snapshot["schema"], "rustred.running-application-cache-snapshot.v1"
        )
        self.assertEqual(Path(snapshot["directory"]), self.directory / "delta")
        self.assertEqual(snapshot["binding"]["physical_arity"], 15)
        self.assertTrue(snapshot["binding"]["family_fingerprint"])
        self.assertEqual(Path(snapshot["delta"]["manifest_path"]), sidecar)
        self.assertEqual(snapshot["delta"]["generation"], generation)
        self.assertEqual(sidecar.stat().st_size, snapshot["delta"]["manifest_bytes"])
        self.assertGreater(sidecar.stat().st_size, 0)
        self.assertEqual(len(snapshot["delta"]["manifest_digest"]), 64)

    def assert_complete_native_result(self, observed, reference):
        # Preserve the public sequence order: no sorting, sets or CAS comparison.
        self.assertIs(type(observed[0]), core.Expression)
        self.assertIs(type(observed[1]), list)
        self.assertIs(type(reference[1]), list)
        self.assertEqual(observed[0], reference[0])
        self.assertEqual(observed[1], reference[1])
        self.assertEqual(observed[2], reference[2])
        self.assertTrue(reference[2])

    def close_control(self, control):
        closed = False
        try:
            closed = control.report()["closed"]
        finally:
            if not closed:
                control.close()

    def assert_latest_monitor(self, path):
        document = json.loads(path.read_text())
        self.assertEqual(document["schema"], "rustred.saved-application-progress.v1")
        self.assertEqual(document["observation_scope"], "caller")
        self.assertEqual(document["activity"], "caller")
        self.assertIs(document["request_finished"], True)
        self.assertEqual(document["status"], "succeeded")
        self.assertEqual(document["core"]["status"], "succeeded")
        self.assertEqual(document["monitor"]["view"]["status"], "Succeeded")
        self.assertIs(document["closure_claim"], False)
        self.assertIs(document["global_remaining_work_known"], False)
        self.assertIsNone(document["eta_seconds"])

    def assert_backing_hit(self):
        stats = self.obj.cache_statistics(5)
        self.assertIs(stats["attached"], True)
        counts = stats["saved_application"]["persistent"]
        self.assertGreater(counts["route_hits"] + counts["apply_hits"], 0)
        self.assertTrue(
            all(not self.obj.cache_statistics(loop)["attached"] for loop in range(1, 5))
        )
        return stats

    def test_control_lifecycle_pinned_snapshot_and_native_monitor_parity(self):
        delta = self.directory / "delta"
        missing = self.directory / "missing"
        with self.assertRaises(ValueError):
            self.obj.load_running_cache(5, missing)
        self.assertFalse(missing.exists())
        with self.assertRaises(ValueError):
            self.obj.attach_running_cache(1, missing)
        with self.assertRaises(ValueError):
            self.obj.attach_running_cache(
                5, missing, checkpoint_interval_seconds=float("nan")
            )
        self.assertFalse(missing.exists())

        reference = self.obj.evaluate_rustred_symbolic(self.input)
        ordinary = self.obj.evaluate(self.input)
        self.assertNotEqual(ordinary, E("0"))
        self.assertEqual(ordinary, reference[0])
        self.obj.clear_cache(5)
        control = self.obj.attach_running_cache(5, delta, read_backing=False)
        # unittest runs every cleanup even if detachment/another cleanup fails;
        # keep the original test exception and attempt each writer's close.
        self.addCleanup(self.close_control, control)
        self.addCleanup(self.obj.detach_cache, 5)
        self.assertIs(type(control), RunningCacheControl)
        self.assert_report(control.report())
        path = self.directory / "capture-monitor.json"
        path.write_text('{"stale":true}\n')
        observed = self.obj.evaluate_rustred_symbolic(
            self.input, progress_snapshot=path, progress_interval_seconds=0
        )
        self.assert_complete_native_result(observed, reference)
        self.assert_latest_monitor(path)
        self.assertNotIn("stale", json.loads(path.read_text()))
        durable = self.obj.save_running_cache(5)
        self.assert_durable(durable)
        self.assertGreater(durable["accepted_records"], 0)
        self.assertEqual(control.flush()["generation"], durable["generation"])
        statistics = self.obj.cache_statistics(5)["saved_application"]
        self.assertEqual(statistics["persistent"]["lookups"], 0)
        self.assertIs(statistics["running"]["read_backing"], False)
        self.assertEqual(
            statistics["running"]["native"]["generation"], durable["generation"]
        )

        sidecar = self.directory / "control-frozen.rrindex"
        snapshot = control.snapshot(sidecar)
        self.assert_snapshot(snapshot, sidecar, durable["generation"])
        second_sidecar = self.directory / "engine-frozen.rrindex"
        second = self.obj.freeze_running_cache(5, second_sidecar)
        self.assert_snapshot(second, second_sidecar, durable["generation"])
        self.assertEqual(snapshot["binding"], second["binding"])
        self.assertEqual(
            snapshot["delta"]["manifest_digest"], second["delta"]["manifest_digest"]
        )
        # The returned object is the public descriptor. The .rrindex file is
        # Native binary metadata, not JSON and not a second Python cache format.
        descriptor = json.loads(json.dumps(snapshot))
        self.obj.detach_cache(5)
        self.assertIs(control.report()["closed"], False)
        self.assert_durable(control.flush())
        self.assert_durable(control.close(), closed=True)
        self.assertIs(control.report()["closed"], True)
        with self.assertRaises(ValueError):
            control.flush()

        reopened = self.obj.load_running_cache(5, delta)
        self.addCleanup(self.close_control, reopened)
        self.addCleanup(self.obj.detach_cache, 5)
        self.assertIs(type(reopened), RunningCacheControl)
        self.obj.clear_cache(5)
        warm = self.obj.evaluate_rustred_symbolic(
            self.input,
            monitor="json",
            progress_snapshot=path,
            progress_interval_seconds=0,
        )
        self.assert_complete_native_result(warm, reference)
        self.assert_latest_monitor(path)
        self.assert_backing_hit()
        resumed = reopened.flush()
        self.assert_durable(resumed)
        self.assertEqual(resumed["accepted_records"], 0)
        self.assertEqual(resumed["generation"], durable["generation"])
        self.obj.detach_cache(5)
        self.assert_durable(reopened.close(), closed=True)

        frozen_files = file_hashes(delta)
        self.obj.load_cache_snapshot(5, snapshots=[descriptor])
        self.obj.clear_cache(5)
        frozen = self.obj.evaluate_rustred_symbolic(self.input)
        self.assert_complete_native_result(frozen, reference)
        self.assertIsNone(self.assert_backing_hit()["saved_application"]["running"])
        self.assertEqual(self.obj.evaluate(self.input), ordinary)
        with self.assertRaises(ValueError):
            self.obj.save_running_cache(5)
        self.obj.detach_cache(5)
        self.assertEqual(file_hashes(delta), frozen_files)


if __name__ == "__main__":
    unittest.main(verbosity=2)
