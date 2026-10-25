"""Public monitor API regressions for one coherently built native host.

The 5L class loads VAKINT_TEST_EXTERNAL_DATA once. Its fast input is the pinched
five-tadpole base of the existing single_positive_line Rust fixture. The explicit
PTY case additionally needs a host built with progress-tui and
VAKINT_MONITOR_TUI_TEST=1; component tests alone do not prove that public call.
"""

import json
import os
import select
import selectors
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

if os.name == "posix":
    import fcntl
    import pty
    import termios

from symbolica import E, core
from symbolica.community.hepkit.vakint import (
    MasterSymbol,
    Vakint,
    VakintEvaluationMethod,
    integral_from_diagram,
)


def engine(**overrides):
    options = {
        "evaluation_order": [
            VakintEvaluationMethod.new_rustred_method(master_mode="retain")
        ],
        "epsilon_symbol": "ep",
        "form_exe_path": "/monitor-python-test-must-not-use-form",
        "integral_normalization_factor": "FMFTandMATAD",
        "use_dot_product_notation": True,
    }
    options.update(overrides)
    return Vakint(**options)


def one_loop_zero():
    return E("topo(I1L(1,2))-(1-ep)*topo(I1L(1,1))", default_namespace="vakint")


def latest(path, status):
    # json.loads rejects concatenated frames as well as partial JSON writes.
    document = json.loads(path.read_text())
    assert document["schema"] == "rustred.saved-application-progress.v1"
    assert document["event"] == "application_progress"
    assert document["observation_scope"] == "caller"
    assert document["activity"] == "caller"
    assert document["request_finished"] is True
    assert document["status"] == status
    assert document["closure_claim"] is False
    assert document["global_remaining_work_known"] is False
    assert document["eta_seconds"] is None
    assert document["percent_complete"] is None
    assert document["arithmetic_internals_known"] is False
    return document


def public_pty_child(output):
    # Popen starts a new session in its native spawn path. Attach this child's
    # inherited fresh PTY here, without an unsafe Python preexec callback.
    fcntl.ioctl(0, termios.TIOCSCTTY, 0)
    assert not sys.flags.optimize and core.is_licensed()
    obj = engine()
    expression = E("topo(I1L(1,2))", default_namespace="vakint")
    reference = obj.evaluate(expression)
    assert reference != E("0")
    print("PUBLIC_MONITOR_READY", flush=True)
    # The controller queues q before the actual public call. Do not consume it:
    # the real terminal worker must handle it, even for a short reduction.
    assert select.select([sys.stdin], [], [], 10)[0]
    observed = obj.evaluate(expression, monitor="tui", progress_interval_seconds=0)
    assert observed == reference
    os.set_blocking(sys.stdin.fileno(), False)
    try:
        remaining = os.read(sys.stdin.fileno(), 4096)
    except BlockingIOError:
        remaining = b""
    assert b"q" not in remaining, (
        "terminal worker did not consume the queued detach key"
    )
    output.write_text(
        json.dumps({"exact_native_value_parity": True, "queued_q_consumed": True})
        + "\n"
    )
    print("PUBLIC_MONITOR_DONE", flush=True)


class PublicMonitorPty(unittest.TestCase):
    @unittest.skipUnless(
        os.name == "posix" and os.environ.get("VAKINT_MONITOR_TUI_TEST") == "1",
        "explicit progress-tui native host and a fresh controlling PTY required",
    )
    def test_q_detach_preserves_the_actual_public_evaluate_result(self):
        with tempfile.TemporaryDirectory(
            prefix="vakint-public-monitor-pty-"
        ) as scratch:
            output = Path(scratch) / "result.json"
            master, slave = pty.openpty()
            termios.tcsetwinsize(slave, (36, 110))
            before = termios.tcgetattr(slave)
            os.set_blocking(master, False)
            selector = selectors.DefaultSelector()
            selector.register(master, selectors.EVENT_READ)
            environment = os.environ.copy()
            environment["TERM"] = "xterm-256color"
            process = None
            transcript = bytearray()
            sent = False
            try:
                process = subprocess.Popen(
                    [sys.executable, "-B", __file__, "--public-pty-child", str(output)],
                    stdin=slave,
                    stdout=slave,
                    stderr=slave,
                    env=environment,
                    start_new_session=True,
                    close_fds=True,
                )
                deadline = time.monotonic() + 30
                while process.poll() is None:
                    self.assertLess(
                        time.monotonic(), deadline, "public PTY child timed out"
                    )
                    for _, _ in selector.select(0.05):
                        try:
                            transcript.extend(os.read(master, 65536))
                        except (BlockingIOError, OSError):
                            pass
                    if not sent and b"PUBLIC_MONITOR_READY" in transcript:
                        # Newline makes the canonical input readable before raw
                        # mode starts; q itself remains queued for the real UI.
                        os.write(master, b"q\n")
                        sent = True
                while True:
                    try:
                        chunk = os.read(master, 65536)
                    except (BlockingIOError, OSError):
                        break
                    if not chunk:
                        break
                    transcript.extend(chunk)
                self.assertEqual(
                    process.returncode, 0, transcript.decode(errors="replace")
                )
                self.assertTrue(sent)
                self.assertIn(b"PUBLIC_MONITOR_DONE", transcript)
                self.assertEqual(
                    json.loads(output.read_text()),
                    {
                        "exact_native_value_parity": True,
                        "queued_q_consumed": True,
                    },
                )
                self.assertEqual(termios.tcgetattr(slave), before)
                self.assertEqual(transcript.count(b"\x1b[?1049h"), 1)
                self.assertEqual(transcript.count(b"\x1b[?1049l"), 1)
                self.assertIn(b"\x1b[?25l", transcript)
                self.assertIn(b"\x1b[?25h", transcript)
            finally:
                if process is not None and process.poll() is None:
                    process.kill()
                    process.wait()
                termios.tcsetattr(slave, termios.TCSANOW, before)
                selector.close()
                os.close(master)
                os.close(slave)


class ReductionMonitorSmoke(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if sys.flags.optimize or not core.is_licensed():
            raise RuntimeError(
                "monitor regressions need unoptimized Python and a licensed host"
            )
        cls.scratch = tempfile.TemporaryDirectory(prefix="vakint-monitor-python-")
        cls.directory = Path(cls.scratch.name)
        cls.obj = engine(temporary_directory=str(cls.directory / "native"))
        data = Path(os.environ["VAKINT_TEST_EXTERNAL_DATA"])
        report = cls.obj.load_official_external_data(
            5, local_artifacts_path=data, load_prepared_caches=False
        )
        assert report["load_prepared_caches"] is False
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
        # The pinched side of single_positive_line has the same five tadpoles
        # and powers [1,1,1,1,5]; auxiliary family coordinates remain pinched.
        powers = [1, 1, 1, 1, 5] + [0] * (len(family.denominators) - 5)
        cls.input = integral_from_diagram(
            graph, family, E("1"), powers=powers
        ).to_expression()

    @classmethod
    def tearDownClass(cls):
        cls.scratch.cleanup()

    def setUp(self):
        self.case = tempfile.TemporaryDirectory(prefix="case-", dir=self.directory)
        self.path = Path(self.case.name) / "latest.json"

    def tearDown(self):
        self.case.cleanup()

    def test_default_off_creates_no_snapshot(self):
        reference = self.obj.evaluate(self.input)
        observed = self.obj.evaluate(self.input, monitor="off")
        self.assertEqual(observed, reference)
        self.assertNotEqual(reference, E("0"))
        self.assertEqual(list(self.path.parent.iterdir()), [])

    def test_path_alone_preserves_usual_retained_value_and_replaces_latest(self):
        reference = self.obj.evaluate(self.input)
        self.path.write_text('{"stale":true}\n')
        observed = self.obj.evaluate(
            self.input, progress_snapshot=self.path, progress_interval_seconds=0
        )
        self.assertEqual(observed, reference)
        document = latest(self.path, "succeeded")
        self.assertEqual(document["core"]["status"], "succeeded")
        self.assertNotIn("stale", document)
        self.assertEqual(set(self.path.parent.iterdir()), {self.path})

    def test_special_symbolic_preserves_whole_value_ordered_guards_and_families(self):
        reference = self.obj.evaluate_rustred_symbolic(self.input)
        observed = self.obj.evaluate_rustred_symbolic(
            self.input,
            monitor="json",
            progress_snapshot=self.path,
            progress_interval_seconds=0,
        )
        self.assertEqual(observed[0], reference[0])
        self.assertEqual(observed[1], reference[1])
        self.assertEqual(observed[2], reference[2])
        self.assertTrue(reference[1])
        self.assertTrue(reference[2])
        self.assertEqual(latest(self.path, "succeeded")["core"]["status"], "succeeded")

    def test_real_snapshot_fault_preserves_native_value(self):
        self.path.mkdir()
        expression = E("topo(I1L(1,2))", default_namespace="vakint")
        reference = self.obj.evaluate(expression)
        observed = self.obj.evaluate(
            expression, progress_snapshot=self.path, progress_interval_seconds=0
        )
        self.assertEqual(observed, reference)
        self.assertNotEqual(reference, E("0"))
        self.assertTrue(self.path.is_dir())
        self.assertEqual(list(self.path.iterdir()), [])
        self.assertEqual(set(self.path.parent.iterdir()), {self.path})


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] == "--public-pty-child":
        public_pty_child(Path(sys.argv[2]))
    else:
        unittest.main()
