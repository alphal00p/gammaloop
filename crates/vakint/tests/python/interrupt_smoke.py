"""Ctrl-C, progress callbacks and the GIL during long native calls.

Run with the installed native host: ``python -m unittest interrupt_smoke``.
Signals are sent with ``os.kill`` from a timer thread while the main thread is
inside a native call; ``KeyboardInterrupt`` must arrive within a fraction of a
second, and the same engine must afterwards return the recorded result.

Classes:
- ``LowerLoop``: the vendored 1L-3L data only (the gallery's 3L quadratic
  ladder input).
- ``DeepInput``: deeply nested input evaluates on the request thread as on
  the calling thread (child processes).
- ``RestrictedSymbolica``: the same in a child process without a Symbolica
  license (one Symbolica thread: the work stays on the calling thread).
- ``FourLoop``: ``VAKINT_TEST_EXTERNAL_DATA`` (a data checkout with 4L) and
  the 4L quadratic ladder input (``VAKINT_TEST_LADDER_INPUTS``, a directory
  with ``4L-quadratic.vakint-input.symbolica``); about 22 s per reduction on one
  core. ``VAKINT_INTERRUPT_TEST_CORES`` (default ``1,4``) selects core counts.
- ``FiveLoop``: opt-in ``VAKINT_INTERRUPT_TEST_5L=1`` (about 15 GB): the shipped
  5L zero-momentum benchmark, interrupted during loading and reduction.

``VAKINT_INTERRUPT_TEST_REPORT`` names a JSON file for the measured latencies.
"""

import hashlib
import json
import os
import resource
import signal
import sys
import threading
import time
import unittest
from pathlib import Path

from symbolica import E
from symbolica.community.hepkit.vakint import (
    Vakint,
    VakintEvaluationMethod,
    VakintProgress,
)

HERE = Path(__file__).resolve().parent
LADDER = HERE.parent.parent / "examples" / "notebooks" / "data" / "ladder"
CASES = json.loads((LADDER / "cases.json").read_text())["cases"]
DATA = os.environ.get("VAKINT_TEST_EXTERNAL_DATA")
INPUTS = os.environ.get("VAKINT_TEST_LADDER_INPUTS")
REPORT = os.environ.get("VAKINT_INTERRUPT_TEST_REPORT")
CORES = [
    int(value)
    for value in os.environ.get("VAKINT_INTERRUPT_TEST_CORES", "1,4").split(",")
    if value.strip()
]
# Latency budget of a Ctrl-C in a checkpointed stage; above HARD it fails.
TARGET_SECONDS = 1.0
HARD_SECONDS = 2.0
MEASUREMENTS = []


def record(**fields):
    MEASUREMENTS.append(fields)
    print(json.dumps(fields), file=sys.stderr, flush=True)
    if REPORT:
        Path(REPORT).write_text(json.dumps(MEASUREMENTS, indent=1) + "\n")


def sha256(text):
    return hashlib.sha256(text.encode()).hexdigest()


def engine(loop_count, **overrides):
    """The gallery's ladder engine (retained masters, MSbar)."""
    options = {
        "evaluation_order": [
            VakintEvaluationMethod.new_rustred_method(master_mode="retain")
        ],
        "run_time_decimal_precision": 32,
        "integral_normalization_factor": "MSbar",
        "form_exe_path": "/interrupt-test-must-not-use-form",
        "use_dot_product_notation": True,
        "number_of_terms_in_epsilon_expansion": loop_count + 1,
        "ordered_cache": "off",
    }
    options.update(overrides)
    return Vakint(**options)


def ladder_input(key, directory):
    case = CASES[key]
    payload = (Path(directory) / case["input"]["file"]).read_bytes()
    assert hashlib.sha256(payload).hexdigest() == case["input"]["sha256"]
    source = E(payload.decode(), default_namespace="qcd_ladder")
    binding = case["dimension_binding"]
    evaluated = source.replace(E(binding["source"]), E(binding["target"]))
    assert sha256(evaluated.to_canonical_string()) == case["evaluated_input_sha256"]
    return evaluated, case["retained_sha256"]


def thread_names():
    names = []
    for task in Path("/proc/self/task").iterdir():
        try:
            names.append((task / "comm").read_text().strip())
        except OSError:
            pass
    return names


def rss_bytes():
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    return 0


def interrupt_after(seconds, call):
    """Run ``call()`` on this (main) thread; SIGINT after ``seconds``.

    Returns (outcome, latency seconds from the signal to the exception, call
    seconds). ``outcome`` is "interrupted" or "finished" (the call returned
    before the signal)."""
    sent = []
    done = threading.Event()

    def fire():
        if not done.wait(seconds):
            sent.append(time.monotonic())
            os.kill(os.getpid(), signal.SIGINT)

    timer = threading.Thread(target=fire, daemon=True)
    started = time.monotonic()
    caught = None
    timer.start()
    try:
        try:
            call()
        except KeyboardInterrupt:
            caught = time.monotonic()
        done.set()
        timer.join()
        if sent and caught is None:
            # A signal that arrived after the native call surfaces here.
            time.sleep(0.5)
    except KeyboardInterrupt:
        pass
    finally:
        done.set()
        timer.join()
    if caught is None:
        return "finished", None, time.monotonic() - started
    return "interrupted", caught - sent[0], caught - started


class GilTicker:
    """A Python thread that ticks every 10 ms; counts ticks per second."""

    def __init__(self):
        self.ticks = []
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self.run, daemon=True)

    def run(self):
        while not self.stop.is_set():
            self.ticks.append(time.monotonic())
            time.sleep(0.01)

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *_):
        self.stop.set()
        self.thread.join()

    def worst_rate(self, start, end):
        """Fewest ticks in any whole second of [start, end)."""
        rates = []
        second = start
        while second + 1.0 <= end:
            rates.append(sum(second <= t < second + 1.0 for t in self.ticks))
            second += 1.0
        return min(rates) if rates else None


class LowerLoop(unittest.TestCase):
    """Vendored data only: the 3L quadratic ladder (ten terms)."""

    @classmethod
    def setUpClass(cls):
        cls.input, cls.expected = ladder_input("3L-quadratic", LADDER)

    def slow_engine(self):
        # Exact arithmetic on one core: the longest lower-loop reduction.
        return engine(3, n_cores=1, ordered_weight_arithmetic="exact")

    def test_sigint_interrupts_and_engine_returns_the_recorded_result(self):
        obj = self.slow_engine()
        reference = obj.evaluate(self.input, progress=None)
        self.assertEqual(sha256(reference.to_canonical_string()), self.expected)
        # The first call also prepares lower-loop data: time a warm one.
        started = time.monotonic()
        obj.evaluate(self.input, progress=None)
        duration = time.monotonic() - started
        outcome, latency, _ = interrupt_after(
            duration / 2, lambda: obj.evaluate(self.input, progress=None)
        )
        record(test="3L sigint", outcome=outcome, latency=latency, duration=duration)
        self.assertEqual(outcome, "interrupted")
        self.assertLess(latency, HARD_SECONDS)
        again = obj.evaluate(self.input, progress=None)
        self.assertEqual(sha256(again.to_canonical_string()), self.expected)

    def test_progress_callback_receives_snapshots_and_final_status(self):
        obj = self.slow_engine()
        snapshots = []
        result = obj.evaluate(
            self.input, progress=snapshots.append, progress_interval_seconds=0.1
        )
        self.assertEqual(sha256(result.to_canonical_string()), self.expected)
        self.assertGreaterEqual(len(snapshots), 2, snapshots)
        self.assertTrue(all(isinstance(s, VakintProgress) for s in snapshots))
        self.assertEqual(snapshots[-1].status, "succeeded")
        self.assertEqual(snapshots[-1].to_dict()["status"], "Succeeded")
        self.assertTrue(all(s.status == "running" for s in snapshots[:-1]))
        self.assertTrue(all(s.summary for s in snapshots))
        elapsed = [s.elapsed_seconds for s in snapshots]
        self.assertEqual(elapsed, sorted(elapsed))
        self.check_phase_order(snapshots)
        record(
            test="3L progress",
            snapshots=len(snapshots),
            stages=sorted({s.stage for s in snapshots}),
            last=snapshots[-1].summary,
        )

    def check_phase_order(self, snapshots):
        # Within one native request the Ordered schedule only moves forward.
        last = {}
        for snapshot in snapshots:
            model = snapshot.to_dict()
            ordered = model.get("ordered")
            if not ordered:
                continue
            request = ordered.get("request")
            index = ordered["phase"]["index"]
            self.assertGreaterEqual(index, last.get(request, 0), snapshot.summary)
            last[request] = index
            self.assertEqual(snapshot.phase, (index, ordered["phase"]["total"]))

    def test_callback_exception_cancels_the_call_and_is_reraised(self):
        obj = self.slow_engine()

        class Stop(Exception):
            pass

        calls = []

        def callback(snapshot):
            calls.append(time.monotonic())
            if len(calls) == 2:
                raise Stop("stop now")

        started = time.monotonic()
        with self.assertRaises(Stop):
            obj.evaluate(self.input, progress=callback, progress_interval_seconds=0)
        left = time.monotonic()
        self.assertEqual(len(calls), 2)
        record(test="3L callback stop", latency=left - calls[-1], total=left - started)
        self.assertLess(left - calls[-1], HARD_SECONDS)
        self.assertEqual(
            sha256(obj.evaluate(self.input, progress=None).to_canonical_string()),
            self.expected,
        )

    def test_progress_must_be_auto_callable_or_none(self):
        obj = engine(1, evaluation_order=[])
        for value in ("tui", 3, True):
            with self.subTest(value=value), self.assertRaises(TypeError):
                obj.evaluate(E("1"), progress=value)

    def test_auto_outside_marimo_starts_no_monitor(self):
        self.assertNotIn("marimo", sys.modules)
        obj = self.slow_engine()
        seen = set()
        stop = threading.Event()

        def sample():
            while not stop.is_set():
                seen.update(thread_names())
                time.sleep(0.005)

        sampler = threading.Thread(target=sample, daemon=True)
        sampler.start()
        try:
            obj.evaluate(self.input)  # progress="auto" (the default)
        finally:
            stop.set()
            sampler.join()
        self.assertIn("vakint-request", seen)
        self.assertNotIn("vakint-monitor", seen)
        # A requested monitor does start one (the sampler sees it).
        seen.clear()
        stop.clear()
        sampler = threading.Thread(target=sample, daemon=True)
        sampler.start()
        try:
            obj.evaluate(self.input, progress=lambda snapshot: None)
        finally:
            stop.set()
            sampler.join()
        self.assertIn("vakint-monitor", seen)

    def test_other_threads_run_on_the_caller_unless_progress_is_consumed(self):
        # Python delivers no signal to this thread (as a marimo app kernel):
        # without progress nothing can cancel, so no request thread starts.
        obj = self.slow_engine()
        outcome = {}

        def call():
            # The work's CPU time lands on this thread when it runs here.
            cpu, wall = time.thread_time(), time.monotonic()
            value = obj.evaluate(self.input, progress=None)
            outcome["cpu"], outcome["wall"] = (
                time.thread_time() - cpu,
                time.monotonic() - wall,
            )
            outcome["sha"] = sha256(value.to_canonical_string())

        worker = threading.Thread(target=call)
        worker.start()
        worker.join()
        self.assertEqual(outcome["sha"], self.expected)
        self.assertGreater(outcome["cpu"], 0.5 * outcome["wall"], outcome)
        # On the main thread (a cancellable call) the work runs elsewhere.
        cpu, wall = time.thread_time(), time.monotonic()
        obj.evaluate(self.input, progress=None)
        self.assertLess(time.thread_time() - cpu, 0.5 * (time.monotonic() - wall))
        self.assertEqual(obj.interrupt, "fine")
        # "off": not even the main thread installs a cancellation scope.
        off = engine(3, n_cores=1, ordered_weight_arithmetic="exact", interrupt="off")
        self.assertEqual(off.interrupt, "off")
        value = off.evaluate(self.input, progress=None)
        self.assertEqual(sha256(value.to_canonical_string()), self.expected)
        with self.assertRaises(ValueError):
            engine(3, interrupt="sometimes")


DEEP_CHILD = """
import json, sys
sys.path.insert(0, sys.argv[1])
import interrupt_smoke as t
from symbolica import E, S
mode, depth, interrupt = sys.argv[2], int(sys.argv[3]), sys.argv[4]
obj = t.engine(1, n_cores=1, verify_numerator_identification=False, interrupt=interrupt)
if mode == "sum":
    muvsq = E("muvsq", default_namespace="vakint")
    coefficient = E("dot(k(1),k(1))", default_namespace="vakint")
    for _ in range(depth):
        coefficient = 1 + muvsq * coefficient
    expression = coefficient * E("topo(I1L(muvsq,4))", default_namespace="vakint")
else:
    g = S("deepinput::g")
    coefficient = E("deepinput::x")
    for _ in range(depth):
        coefficient = g(coefficient)
    expression = coefficient * E("topo(I1L(muvsq,2))", default_namespace="vakint")
value = obj.evaluate(expression, progress=None)
print(json.dumps({"sha": t.sha256(value.to_canonical_string())}))
"""


class DeepInput(unittest.TestCase):
    """Deeply nested input that the calling thread's stack (8 MiB on Python's
    main thread) evaluates also evaluates on the request thread: Symbolica's
    normalization recurses through it. Each case runs in a child process, so a
    stack overflow is a failure here rather than a dead test runner."""

    def evaluate(self, mode, depth, interrupt):
        import subprocess

        child = subprocess.run(
            [sys.executable, "-c", DEEP_CHILD, str(HERE), mode, str(depth), interrupt],
            check=False,
            capture_output=True,
            text=True,
            timeout=900,
        )
        self.assertEqual(
            child.returncode,
            0,
            f"{mode} depth {depth} ({interrupt}): rc {child.returncode}\n"
            + child.stderr[-1500:],
        )
        return json.loads(child.stdout.strip().splitlines()[-1])["sha"]

    def test_nested_numerator_and_function_coefficient(self):
        # Before the request thread got the caller's stack, the nested sum
        # crashed the process at depth 400 and the nested function at 4000.
        for mode, depth in (("sum", 800), ("function", 4000)):
            with self.subTest(mode=mode, depth=depth):
                on_request_thread = self.evaluate(mode, depth, "fine")
                on_caller = self.evaluate(mode, depth, "off")
                self.assertEqual(on_request_thread, on_caller)
                record(test=f"deep {mode} {depth}", ok=True)


RESTRICTED_CHILD = """
import json, sys, time
sys.path.insert(0, sys.argv[1])
import interrupt_smoke as t
from symbolica import core
assert not core.is_licensed()
evaluated, expected = t.ladder_input("3L-quadratic", t.LADDER)
obj = t.engine(3, n_cores=1, ordered_weight_arithmetic="exact")
reference = obj.evaluate(evaluated)
# The first call also prepares lower-loop data: time a warm one.
started = time.monotonic()
obj.evaluate(evaluated)
duration = time.monotonic() - started
outcome, latency, _ = t.interrupt_after(duration / 2, lambda: obj.evaluate(evaluated))
snapshots = []
again = obj.evaluate(evaluated, progress=snapshots.append)
print(json.dumps({
    "duration": duration, "outcome": outcome, "latency": latency,
    "reference_ok": t.sha256(reference.to_canonical_string()) == expected,
    "again_ok": t.sha256(again.to_canonical_string()) == expected,
    "final_status": snapshots[-1].status if snapshots else None,
}))
"""


class RestrictedSymbolica(unittest.TestCase):
    """Unlicensed Symbolica allows one Symbolica thread per user: the work
    stays on the calling thread and Ctrl-C still interrupts it."""

    def test_unlicensed_calls_stay_on_the_caller_and_stay_interruptible(self):
        import subprocess

        environment = {
            key: value
            for key, value in os.environ.items()
            if not key.startswith("SYMBOLICA_LICENSE")
        }
        child = subprocess.run(
            [sys.executable, "-c", RESTRICTED_CHILD, str(HERE)],
            check=False,
            env=environment,
            capture_output=True,
            text=True,
            timeout=600,
        )
        if "thread allowance is in use" in child.stdout and child.returncode:
            self.skipTest(
                "another unlicensed Symbolica process of this user is running"
            )
        self.assertEqual(
            child.returncode, 0, child.stdout[-2000:] + child.stderr[-2000:]
        )
        result = json.loads(child.stdout.strip().splitlines()[-1])
        record(test="restricted 3L sigint", **result)
        self.assertTrue(result["reference_ok"] and result["again_ok"], result)
        self.assertEqual(result["final_status"], "succeeded")
        self.assertEqual(result["outcome"], "interrupted", result)
        self.assertLess(result["latency"], HARD_SECONDS)


@unittest.skipUnless(
    DATA and INPUTS, "VAKINT_TEST_EXTERNAL_DATA and VAKINT_TEST_LADDER_INPUTS needed"
)
class FourLoop(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.input, cls.expected = ladder_input("4L-quadratic", INPUTS)

    def loaded(self, n_cores, **options):
        obj = engine(4, n_cores=n_cores, **options)
        with GilTicker() as ticker:
            started = time.monotonic()
            obj.load_official_external_data(
                4, version="latest", local_artifacts_path=DATA, progress=None
            )
            ended = time.monotonic()
        record(
            test="4L load",
            n_cores=n_cores,
            seconds=ended - started,
            worst_ticks_per_second=ticker.worst_rate(started, ended),
        )
        worst = ticker.worst_rate(started, ended)
        if worst is not None:
            self.assertGreaterEqual(worst, 50)
        return obj

    def test_sigint_latency_and_reuse(self):
        for n_cores in CORES:
            with self.subTest(n_cores=n_cores):
                obj = self.loaded(n_cores)
                baseline_threads = len(thread_names())
                started = time.monotonic()
                with GilTicker() as ticker:
                    reference = obj.evaluate(self.input, progress=None)
                    ended = time.monotonic()
                duration = ended - started
                self.assertEqual(sha256(reference.to_canonical_string()), self.expected)
                worst = ticker.worst_rate(started, ended)
                record(
                    test="4L evaluate",
                    n_cores=n_cores,
                    seconds=duration,
                    worst_ticks_per_second=worst,
                )
                if worst is not None:
                    self.assertGreaterEqual(worst, 50)
                rss = []
                for after in (2.0, 8.0, 15.0):
                    if after >= duration - 1.0:
                        continue
                    outcome, latency, _ = interrupt_after(
                        after, lambda obj=obj: obj.evaluate(self.input, progress=None)
                    )
                    record(
                        test="4L sigint",
                        n_cores=n_cores,
                        after=after,
                        outcome=outcome,
                        latency=latency,
                    )
                    self.assertEqual(outcome, "interrupted")
                    self.assertLess(latency, HARD_SECONDS)
                    again = obj.evaluate(self.input, progress=None)
                    self.assertEqual(sha256(again.to_canonical_string()), self.expected)
                    deadline = time.monotonic() + 2.0
                    while (
                        len(thread_names()) > baseline_threads
                        and time.monotonic() < deadline
                    ):
                        time.sleep(0.05)
                    self.assertLessEqual(len(thread_names()), baseline_threads)
                    rss.append(rss_bytes())
                record(test="4L rss after reruns", n_cores=n_cores, rss=rss)
                if len(rss) >= 2:
                    self.assertLessEqual(rss[-1], rss[0] * 1.05)

    def test_sigint_during_a_load_keeps_the_previous_data(self):
        obj = engine(4, n_cores=max(CORES))
        snapshots = []
        outcome, latency, _ = interrupt_after(
            3.0,
            lambda: obj.load_official_external_data(
                4,
                version="latest",
                local_artifacts_path=DATA,
                progress=snapshots.append,
            ),
        )
        record(
            test="4L load sigint",
            outcome=outcome,
            latency=latency,
            final=snapshots[-1].status if snapshots else None,
        )
        self.assertEqual(outcome, "interrupted")
        # Coarse checkpoints: between data files and program load steps.
        self.assertLess(latency, 5.0)
        self.assertEqual(snapshots[-1].status, "interrupted")
        self.assertEqual(snapshots[-1].stage, "data_load")
        # Without a progress consumer the load stays on this thread and its
        # own checkpoints consume the SIGINT: same stop, same exception.
        outcome, latency, _ = interrupt_after(
            3.0,
            lambda: obj.load_official_external_data(
                4, version="latest", local_artifacts_path=DATA, progress=None
            ),
        )
        record(test="4L load sigint (no progress)", outcome=outcome, latency=latency)
        self.assertEqual(outcome, "interrupted")
        self.assertLess(latency, 5.0)
        # The engine kept its previous (vendored 1L-3L) data and loads again.
        three, _ = ladder_input("3L-quadratic", LADDER)
        self.assertNotEqual(obj.evaluate(three, progress=None), E("0"))
        # A progress callback may use the engine while it loads (the load
        # publishes its data only when it finishes).
        seen_cores = []
        obj.load_official_external_data(
            4,
            version="latest",
            local_artifacts_path=DATA,
            progress=lambda snapshot: seen_cores.append(obj.n_cores),
        )
        self.assertTrue(seen_cores)
        self.assertEqual(set(seen_cores), {max(CORES)})
        result = obj.evaluate(self.input, progress=None)
        self.assertEqual(sha256(result.to_canonical_string()), self.expected)

    def test_progress_snapshots_and_cache_integrity_after_interrupt(self):
        # An in-memory read-write step cache: an interrupted request must not
        # leave records that change a later result.
        obj = self.loaded(max(CORES), ordered_cache="read-write")
        outcome, latency, _ = interrupt_after(
            8.0, lambda: obj.evaluate(self.input, progress=None)
        )
        record(test="4L sigint (read-write cache)", outcome=outcome, latency=latency)
        self.assertEqual(outcome, "interrupted")
        snapshots = []
        result = obj.evaluate(
            self.input, progress=snapshots.append, progress_interval_seconds=0.5
        )
        self.assertEqual(sha256(result.to_canonical_string()), self.expected)
        self.assertGreaterEqual(len(snapshots), 3)
        self.assertEqual(snapshots[-1].status, "succeeded")
        LowerLoop.check_phase_order(self, snapshots)
        record(
            test="4L progress",
            snapshots=len(snapshots),
            stages=sorted({s.stage for s in snapshots}),
            phases=max((s.phase[1] for s in snapshots if s.phase), default=None),
            last=snapshots[-1].summary,
        )


@unittest.skipUnless(
    DATA and os.environ.get("VAKINT_INTERRUPT_TEST_5L") == "1",
    "opt-in: VAKINT_INTERRUPT_TEST_5L=1 and VAKINT_TEST_EXTERNAL_DATA (about 15 GB)",
)
class FiveLoop(unittest.TestCase):
    def test_sigint_during_load_and_reduction(self):
        obj = engine(5, n_cores=max(CORES), ordered_cache="read-write")
        # Interrupt the load, then load for real on the same engine.
        outcome, latency, seconds = interrupt_after(
            20.0,
            lambda: obj.load_official_external_data(
                5, version="latest", local_artifacts_path=DATA, progress=None
            ),
        )
        # Loads stop at Vakint's coarse checkpoints: the native five-loop
        # application load step finishes first (recorded, not asserted).
        record(test="5L load sigint", outcome=outcome, latency=latency, seconds=seconds)
        self.assertEqual(outcome, "interrupted")
        started = time.monotonic()
        obj.load_official_external_data(
            5, version="latest", local_artifacts_path=DATA, progress=None
        )
        record(test="5L load", seconds=time.monotonic() - started)
        benchmark = Path(DATA) / "5L" / "benchmarks" / "zero-momentum"
        payload = (benchmark / "vakint-input.symbolica").read_text()
        source = E(payload, default_namespace="qcd_ladder")
        evaluated = source.replace(E("rustred::{}::d"), E("-2*vakint::{}::ε+4"))
        for after in [
            float(value)
            for value in os.environ.get(
                "VAKINT_INTERRUPT_TEST_5L_AFTER", "30,90,150"
            ).split(",")
        ]:
            snapshots = []
            outcome, latency, _ = interrupt_after(
                after,
                lambda snapshots=snapshots: obj.evaluate(
                    evaluated, progress=snapshots.append, progress_interval_seconds=1.0
                ),
            )
            stage = snapshots[-1].summary if snapshots else None
            record(
                test="5L sigint",
                after=after,
                outcome=outcome,
                latency=latency,
                stage=stage,
                peak_rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
            )
            if outcome == "interrupted" and latency >= HARD_SECONDS:
                # A checkpoint-free stretch (e.g. the modular lift): report it.
                print(f"5L: {latency:.1f} s at {stage}", file=sys.stderr)
        started = time.monotonic()
        result = obj.evaluate(evaluated, progress=None)
        record(
            test="5L evaluate after interrupts",
            seconds=time.monotonic() - started,
            retained_sha256=sha256(result.to_canonical_string()),
        )


if __name__ == "__main__":
    unittest.main()
