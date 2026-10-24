"""Run against ONE coherently built Symbolica/HepKit extension, not an embedded kernel.

Set TMPDIR to an owned scratch directory. VAKINT_CACHE_SMOKE_KEEP=1 retains
packs/reports for auditing. These bounded one-loop controls compute no new
master values. This script is not proof of four-loop runtime completion.
"""

import hashlib
import os
import shutil
import signal
import sys
import tarfile
import tempfile
import threading
import unittest
from pathlib import Path

if os.environ.get("VAKINT_TEST_HOST") == "standalone":
    from vakint_hepkit_standalone import E, Vakint, VakintEvaluationMethod
else:
    from symbolica import E
    from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod

ROOTS = [(0, [power]) for power in range(5)]


def engine(**overrides):
    options = {
        "evaluation_order": [
            VakintEvaluationMethod.new_rustred_method(master_mode="retain")
        ],
        "form_exe_path": "/native-cache-python-test-must-not-use-form",
        "integral_normalization_factor": "FMFTandMATAD",
        "use_dot_product_notation": True,
    }
    options.update(overrides)
    return Vakint(**options)


def integral(power, mass=4):
    return E(f"topo(I1L({mass},{power}))", default_namespace="vakint")


class ProgressFailure(RuntimeError):
    pass


class NativeCacheSmoke(unittest.TestCase):
    def setUp(self):
        self.scratch = Path(tempfile.mkdtemp(prefix="vakint-cache-python-"))

    def tearDown(self):
        if os.environ.get("VAKINT_CACHE_SMOKE_KEEP") == "1":
            print(f"retained_cache_fixture={self.scratch}", flush=True)
        else:
            shutil.rmtree(self.scratch)

    def prepare(self, obj, path, **kwargs):
        return obj.prepare_cache(
            1, str(path), workload=ROOTS, time_limit_seconds=30, **kwargs
        )

    def same(self, left, right):
        self.assertEqual((left - right).expand().together(), E("0"))

    @unittest.skipUnless(
        sys.platform in ("linux", "darwin"), "native cache transfer platform"
    )
    def test_checksum_import_is_explicit_and_preserves_native_records(self):
        obj = engine()
        prepared = self.scratch / "prepared-transfer"
        self.assertEqual(self.prepare(obj, prepared)["status"], "complete")
        archive = self.scratch / "cache.tar.gz"
        with tarfile.open(archive, "w:gz", format=tarfile.USTAR_FORMAT) as bundle:
            bundle.add(prepared, arcname="cache")
        digest = hashlib.sha256(archive.read_bytes()).hexdigest()
        imported = self.scratch / "imported"
        report = obj.import_cache(1, str(archive), str(imported), sha256=digest)
        self.assertEqual(report["sha256"], digest)
        self.assertEqual(report["preserved_record_references"], 5)
        self.assertFalse(obj.cache_statistics(1)["attached"])
        self.assertFalse((imported / "preparation.rrcheckpoint").exists())
        obj.load_cache(1, str(imported))
        self.same(
            obj.evaluate_integral(integral(4)), engine().evaluate_integral(integral(4))
        )
        self.assertEqual(obj.cache_statistics(1)["native"]["backing_hits"], 1)
        for checksum in ["not-a-sha256", "0" * 64]:
            failed = self.scratch / "failed-import"
            with self.assertRaises(ValueError):
                obj.import_cache(1, str(archive), str(failed), sha256=checksum)
            self.assertFalse(failed.exists())
        with self.assertRaises(ValueError):
            obj.import_cache(1, str(archive), str(imported), sha256=digest)

    def test_native_hit_full_snapshot_and_reentrant_progress(self):
        obj = engine()
        events = []

        def progress(report):
            # The callback can query the same wrapper without holding a native
            # cache lock across the Python call.
            self.assertFalse(obj.cache_statistics(1)["attached"])
            events.append(report["status"])

        prepared = self.scratch / "prepared"
        report = self.prepare(obj, prepared, progress=progress)
        self.assertEqual(report["status"], "complete")
        self.assertEqual(report["completed_roots"], 5)
        self.assertEqual(report["options"]["batch_records"], 64)
        self.assertEqual(events[0], "running")
        self.assertEqual(events[-1], "complete")
        obj.load_cache(1, str(prepared), cache_policy="evict-aged-lfu")
        self.assertEqual(obj.cache_statistics(1)["disk"]["chunks_read"], 0)
        reference = engine()
        cached = obj.evaluate_integral(integral(4))
        self.same(cached, reference.evaluate_integral(integral(4)))
        self.same(cached, obj.evaluate_integral(integral(4)))
        first = obj.cache_statistics(1)
        self.assertEqual(first["native"]["rule_applications"], 0)
        self.assertEqual(first["native"]["backing_hits"], 1)
        self.assertEqual(first["native"]["ram_hits"], 1)
        held_out = obj.evaluate_integral(integral(5))
        snapshot = self.scratch / "snapshot"
        saved = obj.save_cache(1, str(snapshot), mode="full-snapshot")
        self.assertTrue(saved["complete"])
        self.assertEqual(saved["base_copy"]["copied_bases"], 1)
        restored = engine()
        restored.load_cache(1, str(snapshot))
        self.same(held_out, restored.evaluate_integral(integral(5)))
        self.assertEqual(restored.cache_statistics(1)["native"]["rule_applications"], 0)
        # Root2 was never loaded into the source session: full copy preserves it.
        self.same(
            restored.evaluate_integral(integral(2)),
            reference.evaluate_integral(integral(2)),
        )
        self.assertEqual(
            restored.cache_statistics(1)["native"]["completed_reductions"], 0
        )
        restored.clear_cache(1)
        self.assertEqual(restored.cache_statistics(1)["native"]["cached_integrals"], 0)
        self.same(cached, restored.evaluate_integral(integral(4)))
        # Unit-mass cache identity does not contain display epsilon, precision
        # or integral normalization. A new presentation keeps the same pack.
        alternate = engine(
            epsilon_symbol=E("ep_cache_smoke"),
            run_time_decimal_precision=48,
            integral_normalization_factor="MSbar",
        )
        alternate_reference = engine(
            epsilon_symbol=E("ep_cache_smoke"),
            run_time_decimal_precision=48,
            integral_normalization_factor="MSbar",
        )
        alternate.load_cache(1, str(snapshot))
        self.same(
            alternate.evaluate_integral(integral(4, mass=9)),
            alternate_reference.evaluate_integral(integral(4, mass=9)),
        )
        self.assertEqual(
            alternate.cache_statistics(1)["native"]["rule_applications"], 0
        )

    def test_false_progress_and_cancel_predicate_checkpoint_then_resume(self):
        for name, options in [
            ("progress", {"progress": lambda _report: False}),
            ("cancel", {"cancel": lambda: True}),
        ]:
            obj = engine()
            output = self.scratch / name
            partial = self.prepare(obj, output, **options)
            self.assertEqual(partial["status"], "cancelled")
            self.assertEqual(partial["completed_roots"], 0)
            resumed = self.prepare(obj, output)
            self.assertEqual(resumed["status"], "complete")
            self.assertEqual(resumed["completed_roots"], 5)

    def test_detach_preserves_pack_and_resumes_uncached_evaluation(self):
        obj = engine()
        prepared = self.scratch / "detach"
        self.assertEqual(self.prepare(obj, prepared)["status"], "complete")

        def pack_digest():
            return {
                str(path.relative_to(prepared)): hashlib.sha256(
                    path.read_bytes()
                ).hexdigest()
                for path in prepared.rglob("*")
                if path.is_file()
            }

        original_pack = pack_digest()
        obj.load_cache(1, str(prepared))
        cached = obj.evaluate_integral(integral(4))
        self.assertTrue(obj.cache_statistics(1)["attached"])
        self.assertEqual(obj.cache_statistics(1)["native"]["backing_hits"], 1)
        obj.detach_cache(1)
        detached = obj.cache_statistics(1)
        self.assertFalse(detached["attached"])
        self.assertIsNone(detached["disk"])
        self.assertEqual(detached["namespaces"], [])
        self.assertEqual(detached["native"]["cached_integrals"], 0)
        self.assertEqual(detached["native"]["cached_coefficient_terms"], 0)
        self.assertEqual(detached["native"]["cached_coefficient_bytes"], 0)
        self.assertTrue(prepared.is_dir())
        self.assertEqual(pack_digest(), original_pack)
        self.same(cached, obj.evaluate_integral(integral(4)))
        # Detach removes the observed backing backend. Ordinary uncached work
        # uses a request-local reducer outside these attached-backend counters.
        after = obj.cache_statistics(1)
        self.assertFalse(after["attached"])
        self.assertIsNone(after["disk"])
        self.assertEqual(after["namespaces"], [])
        self.assertEqual(after["native"]["cached_integrals"], 0)
        self.assertEqual(after["native"]["cached_coefficient_terms"], 0)
        self.assertEqual(after["native"]["cached_coefficient_bytes"], 0)
        self.assertEqual(after["native"]["backing_hits"], 0)
        self.assertEqual(pack_digest(), original_pack)
        obj.load_cache(1, str(prepared))
        self.same(cached, obj.evaluate_integral(integral(4)))
        self.assertEqual(obj.cache_statistics(1)["native"]["backing_hits"], 1)
        self.assertEqual(pack_digest(), original_pack)

    def test_callback_exception_keeps_original_and_durable_prefix(self):
        obj = engine()
        output = self.scratch / "exception"
        calls = 0

        def progress(_report):
            nonlocal calls
            calls += 1
            if calls == 3:
                raise ProgressFailure("deliberate progress failure")

        with self.assertRaisesRegex(
            ProgressFailure, "deliberate progress failure"
        ) as caught:
            self.prepare(obj, output, progress=progress, progress_interval_seconds=0)
        report = caught.exception.cache_report
        self.assertEqual(report["status"], "cancelled")
        self.assertEqual(report["completed_roots"], 2)
        self.assertIsNone(report["checkpoint_failure"])
        self.assertEqual(calls, 3)
        resumed = self.prepare(obj, output)
        self.assertEqual(resumed["status"], "complete")
        self.assertEqual(resumed["reused_roots"], 2)

    def test_keyboard_interrupt_is_rethrown_after_checkpoint(self):
        obj = engine()
        output = self.scratch / "interrupt"
        calls = 0

        def progress(_report):
            nonlocal calls
            calls += 1
            if calls == 3:
                # Exercise an actual Python signal at a supported cooperative
                # boundary, not an asynchronous interruption inside native CAS.
                signal.raise_signal(signal.SIGINT)

        previous = signal.signal(signal.SIGINT, signal.default_int_handler)
        try:
            with self.assertRaises(KeyboardInterrupt) as caught:
                self.prepare(
                    obj, output, progress=progress, progress_interval_seconds=0
                )
        finally:
            signal.signal(signal.SIGINT, previous)
        self.assertEqual(caught.exception.cache_report["completed_roots"], 2)
        self.assertIsNone(caught.exception.cache_report["checkpoint_failure"])
        self.assertEqual(self.prepare(obj, output)["reused_roots"], 2)

    def test_checkpoint_failure_is_reported_without_hiding_callback_error(self):
        obj = engine()
        output = self.scratch / "checkpoint-failure"
        calls = 0

        def progress(_report):
            nonlocal calls
            calls += 1
            if calls == 3:
                checkpoint = output / "preparation.rrcheckpoint"
                checkpoint.rename(output / "prior-checkpoint")
                checkpoint.mkdir()
                raise ProgressFailure("callback failed before checkpoint publication")

        with self.assertRaises(ProgressFailure) as caught:
            self.prepare(obj, output, progress=progress, progress_interval_seconds=0)
        self.assertIsNotNone(caught.exception.cache_report["checkpoint_failure"])
        self.assertEqual(caught.exception.cache_report["completed_roots"], 2)
        self.assertTrue((output / "prior-checkpoint").is_file())

    def test_partial_native_overrides_preserve_active_ram_and_work(self):
        prepared = self.scratch / "prepared"
        obj = engine()
        self.assertEqual(self.prepare(obj, prepared)["status"], "complete")
        obj.load_cache(1, str(prepared), ram_bytes=1)
        obj.load_cache(1, str(prepared), cache_terms=500)
        with self.assertRaises(ValueError):
            obj.evaluate_integral(integral(4))
        # Change RAM explicitly, then an independent term override must not
        # reset the active zero-rule allowance to ordinary native defaults.
        obj.load_cache(1, str(prepared), ram_bytes=1_000_000, max_rule_applications=0)
        obj.load_cache(1, str(prepared), cache_terms=500)
        obj.evaluate_integral(integral(4))
        with self.assertRaises(ValueError):
            obj.evaluate_integral(integral(5))

    def test_explicit_batch_target_changes_layout_not_reduction(self):
        for batch_records, expected_chunks in ((2, 3), (4097, 1)):
            with self.subTest(batch_records=batch_records):
                obj = engine()
                output = self.scratch / f"batches-{batch_records}"
                report = self.prepare(obj, output, batch_records=batch_records)
                self.assertEqual(report["status"], "complete")
                self.assertEqual(report["completed_roots"], 5)
                self.assertEqual(report["options"]["batch_records"], batch_records)
                obj.load_cache(1, str(output))
                disk = obj.cache_statistics(1)["disk"]
                self.assertEqual(disk["records"], 5)
                self.assertEqual(disk["chunks"], expected_chunks)
                self.assertEqual(disk["chunks_read"], 0)
                reference = engine()
                for power in (1, 2, 4):
                    self.same(
                        obj.evaluate_integral(integral(power)),
                        reference.evaluate_integral(integral(power)),
                    )
                self.assertEqual(
                    obj.cache_statistics(1)["native"]["rule_applications"], 0
                )

    def test_large_existing_chunk_loads_and_resumes_with_smaller_batch_target(self):
        output = self.scratch / "old-large-chunk"
        roots = [(0, [power]) for power in range(67)]
        calls = 0

        def progress(report):
            nonlocal calls
            if report["status"] == "running":
                calls += 1
                return calls <= 65

        partial = engine().prepare_cache(
            1,
            str(output),
            workload=roots,
            batch_records=4096,
            time_limit_seconds=30,
            progress=progress,
            progress_interval_seconds=0,
        )
        self.assertEqual(partial["status"], "cancelled")
        self.assertEqual(partial["completed_roots"], 65)
        self.assertIsNone(partial["checkpoint_failure"])
        self.assertEqual(partial["options"]["batch_records"], 4096)

        fresh = engine()
        fresh.load_cache(1, str(output))
        disk = fresh.cache_statistics(1)["disk"]
        self.assertEqual(disk["chunks"], 1)
        self.assertEqual(disk["records"], 65)
        self.assertEqual(disk["chunks_read"], 0)
        reference = engine()
        self.same(
            fresh.evaluate_integral(integral(64)),
            reference.evaluate_integral(integral(64)),
        )
        self.assertEqual(fresh.cache_statistics(1)["native"]["rule_applications"], 0)

        resumed = engine().prepare_cache(
            1, str(output), workload=roots, time_limit_seconds=30
        )
        self.assertEqual(resumed["status"], "complete")
        self.assertEqual(resumed["options"]["batch_records"], 64)
        self.assertEqual(resumed["completed_roots"], 67)
        self.assertEqual(resumed["reused_roots"], 65)
        self.assertEqual(resumed["newly_stored_roots"], 2)
        fresh.load_cache(1, str(output))
        self.same(
            fresh.evaluate_integral(integral(66)),
            reference.evaluate_integral(integral(66)),
        )
        self.assertEqual(fresh.cache_statistics(1)["native"]["rule_applications"], 0)

    def test_native_work_releases_gil(self):
        obj = engine()
        ready = threading.Event()
        gate = threading.Lock()
        gate.acquire()
        returned = threading.Event()
        observed = threading.Event()

        def worker():
            ready.set()
            with gate:
                if not returned.is_set():
                    observed.set()

        thread = threading.Thread(target=worker)
        thread.start()
        self.assertTrue(ready.wait(timeout=5), "GIL test worker did not start")
        previous = sys.getswitchinterval()
        # Prevent ordinary Python bytecode scheduling just before/after the
        # native call from masquerading as native GIL release.
        sys.setswitchinterval(1000)
        try:
            # release() cannot block acquiring another lock before the call.
            gate.release()
            report = obj.prepare_cache(
                1,
                str(self.scratch / "gil"),
                workload=[(0, [power]) for power in range(33)],
                time_limit_seconds=30,
            )
            returned.set()
        finally:
            returned.set()
            sys.setswitchinterval(previous)
            thread.join(timeout=5)
        self.assertFalse(thread.is_alive())
        self.assertEqual(report["status"], "complete")
        self.assertTrue(
            observed.is_set(), "another Python thread did not run during native work"
        )

    def test_controls_fail_before_unrequested_io(self):
        obj = engine()
        path = self.scratch / "invalid"
        for loop_count in (0, 6):
            with self.assertRaises(ValueError):
                obj.prepare_cache(loop_count, str(path))
        with self.assertRaisesRegex(ValueError, "explicit selected workload"):
            obj.prepare_cache(5, str(path))
        for options in (
            {"cache_policy": "unbounded"},
            {"workload": ROOTS, "extra_dots": 1},
            {"progress": object()},
            {"progress_interval_seconds": float("nan")},
            {"batch_records": 0},
        ):
            with self.assertRaises((TypeError, ValueError)):
                obj.prepare_cache(1, str(path), **options)
        self.assertFalse(path.exists())
        checkpoint = self.scratch / "invalid-checkpoint"
        with self.assertRaisesRegex(ValueError, "only used by the five-loop"):
            obj.prepare_cache(1, str(path), checkpoint_path=checkpoint)
        self.assertFalse(checkpoint.exists())
        self.assertFalse(path.exists())
        with self.assertRaises(ValueError):
            obj.save_cache(1, str(path), mode="automatic")
        self.assertFalse(path.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
