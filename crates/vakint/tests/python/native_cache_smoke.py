"""Run against ONE coherently built Symbolica/HepKit extension, not an embedded kernel.

Set TMPDIR to an owned scratch directory. VAKINT_CACHE_SMOKE_KEEP=1 retains
packs/reports for auditing. These bounded one-loop controls compute no new
master values. This script is not proof of four-loop runtime completion.
"""

import os
import shutil
import signal
import sys
import tempfile
import threading
import unittest
from pathlib import Path

from symbolica import E
from symbolica.community.hepkit.vakint import Vakint, VakintEvaluationMethod

ROOTS = [(0, [power]) for power in range(5)]


def engine(**overrides):
    options = {
        "evaluation_order": [
            VakintEvaluationMethod.new_rustred_method(substitute_masters=False)
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
        for loop_count in (0, 5):
            with self.assertRaises(ValueError):
                obj.prepare_cache(loop_count, str(path))
        for options in (
            {"cache_policy": "unbounded"},
            {"workload": ROOTS, "extra_dots": 1},
            {"progress": object()},
            {"progress_interval_seconds": float("nan")},
        ):
            with self.assertRaises((TypeError, ValueError)):
                obj.prepare_cache(1, str(path), **options)
        self.assertFalse(path.exists())
        with self.assertRaises(ValueError):
            obj.save_cache(1, str(path), mode="automatic")
        self.assertFalse(path.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
