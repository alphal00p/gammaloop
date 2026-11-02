"""Gallery launch, cancellation and display contracts.

No native host runs. Most tests use fake commands; the job-process tests start
short ``sleep`` commands to check cancellation, deadlines and orphan cleanup.
"""

import ast
import contextlib
import fcntl
import hashlib
import importlib.util
import io
import json
import os
import re
import signal
import subprocess
import sys
import tempfile
import time
import unittest
import uuid
from pathlib import Path
from types import MappingProxyType, SimpleNamespace
from unittest.mock import Mock, patch

import benchmark_job
import gallery_jobs
from shared_gallery import Gallery


class GalleryBenchmarkTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="gallery-contract-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.work = self.root / "work"
        self.prepared = self.root / "prepared"
        self.environment = {
            "VAKINT_GALLERY_BENCHMARK": str(self.root / "NOT-AN-EXECUTABLE"),
            "VAKINT_GALLERY_WORK": str(self.work),
            "VAKINT_GALLERY_ARTIFACTS": str(self.root / "artifacts"),
            "VAKINT_GALLERY_PREPARED_INPUT_ROOT": str(self.prepared),
        }
        environment_patch = patch.dict(
            os.environ, {**self.environment, "VAKINT_GALLERY_ALLOW_5L": "1"}, clear=True
        )
        environment_patch.start()
        self.addCleanup(environment_patch.stop)
        for projection in ("quadratic", "zero-momentum"):
            bundle = self.prepared / projection
            bundle.mkdir(parents=True)
            # The native collector validates content. This layer checks that a
            # configured bundle exists before launching that collector.
            (bundle / "prepared-input.json").write_text("{}")

    def test_prepared_input_requires_server_configuration(self):
        for value in (None, ""):
            with self.subTest(value=value):
                if value is None:
                    os.environ.pop("VAKINT_GALLERY_PREPARED_INPUT_ROOT", None)
                else:
                    os.environ["VAKINT_GALLERY_PREPARED_INPUT_ROOT"] = value
                with self.assertRaisesRegex(
                    ValueError, "configured|configuration|ROOT"
                ):
                    Gallery.prepared_input_directory("quadratic")

    def test_prepared_input_rejects_unknown_projection(self):
        for projection in ("../quadratic", "five-loop", ""):
            with (
                self.subTest(projection=projection),
                self.assertRaisesRegex(ValueError, "projection"),
            ):
                Gallery.prepared_input_directory(projection)

    def test_prepared_input_reports_missing_directory(self):
        bundle = self.prepared / "quadratic"
        (bundle / "prepared-input.json").unlink()
        bundle.rmdir()
        with self.assertRaises(ValueError) as caught:
            Gallery.prepared_input_directory("quadratic")
        self.assertIn(str(bundle), str(caught.exception))

    def test_prepared_input_reports_missing_manifest(self):
        manifest = self.prepared / "quadratic" / "prepared-input.json"
        manifest.unlink()
        with self.assertRaises(ValueError) as caught:
            Gallery.prepared_input_directory("quadratic")
        self.assertIn(str(manifest), str(caught.exception))

    def test_prepared_input_resolves_each_projection(self):
        for projection in ("quadratic", "zero-momentum"):
            with self.subTest(projection=projection):
                self.assertEqual(
                    Path(Gallery.prepared_input_directory(projection)),
                    (self.prepared / projection).resolve(),
                )

    def test_availability_requires_configuration_without_reading_bundle(self):
        self.assertTrue(Gallery.benchmark_available())
        # Availability is a cheap configuration check. Content validation occurs
        # on an explicit run, rather than on every reactive page refresh.
        (self.prepared / "quadratic" / "prepared-input.json").write_text("invalid")
        with patch.object(Path, "read_text", side_effect=AssertionError("eager read")):
            self.assertTrue(Gallery.benchmark_available())
        for name in self.environment:
            with self.subTest(name=name), patch.dict(os.environ, {name: ""}):
                self.assertFalse(Gallery.benchmark_available())

    def test_launch_records_projection_bundle_and_background_worker(self):
        projection = "zero-momentum"
        inherited = []

        def spawn(*args, inherited=inherited, **kwargs):
            (descriptor,) = kwargs["pass_fds"]
            os.fstat(descriptor)
            inherited.append(descriptor)

        with patch("gallery_jobs.subprocess.Popen", side_effect=spawn) as popen:
            work = Path(Gallery.start_benchmark(projection))
        command = json.loads((work / "request.json").read_text())["command"]
        self.assertEqual(command[0], self.environment["VAKINT_GALLERY_BENCHMARK"])
        self.assertEqual(command[command.index("--projection") + 1], projection)
        self.assertEqual(
            command[command.index("--prepared-input") + 1],
            str((self.prepared / projection).resolve()),
        )
        self.assertEqual(command[command.index("--work-dir") + 1], str(work))
        self.assertEqual(command.count("--progress-snapshot"), 1)
        self.assertEqual(
            command[command.index("--progress-snapshot") + 1],
            str(work / "monitor.json"),
        )
        self.assertNotIn("--monitor", command)
        self.assertEqual(
            command[command.index("--artifacts") + 1],
            self.environment["VAKINT_GALLERY_ARTIFACTS"],
        )
        # Ordered reduction reads and writes no persistent cache.
        self.assertFalse(
            any(
                argument.startswith(("--cache", "--running-cache"))
                for argument in command
            )
        )
        self.assertTrue(work.is_dir())
        self.assertEqual(
            json.loads((work / "status.json").read_text()),
            {"state": "starting"},
        )
        popen.assert_called_once()
        arguments = popen.call_args.args[0]
        self.assertEqual(arguments[0], sys.executable)
        self.assertEqual(Path(arguments[1]), Gallery.directory / "benchmark_job.py")
        self.assertEqual(arguments[2:], [str(work), str(inherited[0])])
        self.assertEqual(popen.call_args.kwargs["stderr"], subprocess.STDOUT)
        self.assertNotIn("timeout", popen.call_args.kwargs)
        # One process group per job, so Cancel signals exactly this job.
        self.assertTrue(popen.call_args.kwargs["start_new_session"])
        self.assertEqual(popen.call_args.kwargs["stdin"], subprocess.DEVNULL)
        request = json.loads((work / "request.json").read_text())
        self.assertEqual(request["kind"], "five-loop")
        self.assertGreater(request["max_seconds"], 0)
        cores = request["n_cores"]
        self.assertLessEqual(cores, Gallery.max_cores())
        self.assertEqual(request["environment"]["VAKINT_N_CORES"], str(cores))
        self.assertEqual(request["environment"]["RAYON_NUM_THREADS"], str(cores))
        with self.assertRaises(OSError):
            os.fstat(inherited[0])

    def test_five_loop_run_requires_the_operator_opt_in(self):
        for value in (None, "", "0", "yes"):
            with self.subTest(value=value):
                if value is None:
                    os.environ.pop("VAKINT_GALLERY_ALLOW_5L", None)
                else:
                    os.environ["VAKINT_GALLERY_ALLOW_5L"] = value
                self.assertTrue(Gallery.benchmark_available())
                self.assertFalse(Gallery.five_loop_runs_enabled())
                with (
                    patch("gallery_jobs.subprocess.Popen") as popen,
                    self.assertRaisesRegex(ValueError, "explore-only"),
                ):
                    Gallery.start_benchmark("zero-momentum")
                popen.assert_not_called()
                self.assertFalse(self.work.exists())
        os.environ["VAKINT_GALLERY_ALLOW_5L"] = "1"
        self.assertTrue(Gallery.five_loop_runs_enabled())

    def test_quadratic_launch_is_deferred_without_starting_a_worker(self):
        with (
            patch("gallery_jobs.subprocess.Popen") as popen,
            self.assertRaisesRegex(ValueError, "quadratic reduction is deferred"),
        ):
            Gallery.start_benchmark("quadratic")
        popen.assert_not_called()
        self.assertFalse(self.work.exists())

    def test_fresh_output_directories_share_only_prepared_input(self):
        with patch("gallery_jobs.subprocess.Popen"):
            first = Path(Gallery.start_benchmark("zero-momentum"))
            second = Path(Gallery.start_benchmark("zero-momentum"))
        self.assertNotEqual(first, second)
        for work in (first, second):
            self.assertEqual(work.parent, self.work.resolve())
            self.assertTrue((work / "request.json").is_file())
            self.assertNotEqual(work, self.prepared / "zero-momentum")
            command = json.loads((work / "request.json").read_text())["command"]
            self.assertEqual(command[command.index("--work-dir") + 1], str(work))

    def test_held_lock_prevents_another_benchmark(self):
        self.work.mkdir()
        with (self.work / "benchmark.lock").open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with (
                patch("gallery_jobs.subprocess.Popen") as popen,
                self.assertRaisesRegex(RuntimeError, "already running"),
            ):
                Gallery.start_benchmark("zero-momentum")
            popen.assert_not_called()
        self.assertEqual(list(self.work.iterdir()), [self.work / "benchmark.lock"])

    def test_launch_failure_closes_lock_without_retry(self):
        with (
            patch(
                "gallery_jobs.subprocess.Popen",
                side_effect=OSError("synthetic launch error"),
            ) as popen,
            self.assertRaisesRegex(OSError, "synthetic launch error"),
        ):
            Gallery.start_benchmark("zero-momentum")
        popen.assert_called_once()
        (descriptor,) = popen.call_args.kwargs["pass_fds"]
        with self.assertRaises(OSError):
            os.fstat(descriptor)
        with (self.work / "benchmark.lock").open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)


class GalleryWorkerTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="gallery-worker-contract-")
        self.addCleanup(temporary.cleanup)
        self.work = Path(temporary.name)
        self.command = ["NOT-AN-EXECUTABLE", "--prepared-input", "synthetic-bundle"]
        (self.work / "request.json").write_text(json.dumps({"command": self.command}))

    def run_worker(self, *, result=0, error=None, output="", timeouts=0):
        descriptor = os.open(
            self.work / "benchmark.lock", os.O_CREAT | os.O_RDWR, 0o600
        )
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        clock = SimpleNamespace(elapsed=0.0)
        self.running_statuses = []
        self.wait_timeouts = []
        self.child_context_entered = False
        self.child_context_exited = False
        owner = self

        class Child:
            pid = 123456789

            def __enter__(self):
                owner.child_context_entered = True
                return self

            def __exit__(self, *_):
                owner.child_context_exited = True
                return False

            def wait(self, *, timeout):
                owner.assertEqual(timeout, 5)
                os.fstat(descriptor)
                owner.running_statuses.append(
                    json.loads((owner.work / "status.json").read_text())
                )
                owner.wait_timeouts.append(timeout)
                if len(owner.wait_timeouts) <= timeouts:
                    clock.elapsed += timeout
                    raise subprocess.TimeoutExpired(owner.command, timeout)
                return result

        def spawn(command, **kwargs):
            self.assertEqual(command, self.command)
            self.assertNotIn("timeout", kwargs)
            self.assertNotIn("check", kwargs)
            self.assertEqual(kwargs["stderr"], subprocess.STDOUT)
            os.fstat(descriptor)
            kwargs["stdout"].write(output)
            if error is not None:
                raise error
            return Child()

        def resources(pid):
            self.assertEqual(pid, Child.pid)
            return {"sampled_unix_seconds": 1000 + clock.elapsed}

        with (
            patch.object(
                sys, "argv", ["benchmark_job.py", str(self.work), str(descriptor)]
            ),
            patch("benchmark_job.monotonic", side_effect=lambda: clock.elapsed),
            patch("benchmark_job.time", side_effect=lambda: 1000 + clock.elapsed),
            patch("benchmark_job.worker_resources", side_effect=resources),
            patch("benchmark_job.subprocess.Popen", side_effect=spawn) as popen,
        ):
            benchmark_job.main()
        popen.assert_called_once()
        self.assertEqual(self.child_context_entered, error is None)
        self.assertEqual(self.child_context_exited, error is None)
        with self.assertRaises(OSError):
            os.fstat(descriptor)
        with (self.work / "benchmark.lock").open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return json.loads((self.work / "status.json").read_text())

    def test_success_has_no_default_deadline_and_preserves_receipt(self):
        receipt = {"fresh_native_preparation": False, "preparation_seconds": None}
        (self.work / "benchmark.json").write_text(json.dumps(receipt))
        status = self.run_worker(result=0)
        self.assertEqual(status["state"], "complete")
        self.assertEqual(status["receipt"], receipt)
        self.assertGreaterEqual(status["elapsed_seconds"], 0)

    def test_heartbeat_stays_fresh_before_native_progress_without_a_job_deadline(self):
        status = self.run_worker(timeouts=7)
        self.assertEqual(status["state"], "complete")
        self.assertEqual(status["elapsed_seconds"], 35)
        self.assertEqual(self.wait_timeouts, [5] * 8)
        self.assertEqual(
            [row["elapsed_seconds"] for row in self.running_statuses],
            [0, 5, 10, 15, 20, 25, 30, 35],
        )
        for row in self.running_statuses:
            self.assertEqual(row["state"], "running")
            self.assertEqual(row["updated_unix_seconds"], 1000 + row["elapsed_seconds"])
            self.assertEqual(row["native_progress"]["state"], "pending")
            self.assertEqual(row["worker_pid"], 123456789)
            self.assertEqual(
                row["resources"]["sampled_unix_seconds"], row["updated_unix_seconds"]
            )

    def test_heartbeat_does_not_make_a_stale_native_snapshot_look_new(self):
        monitor = self.work / "monitor.json"
        raw = {"phase": "saved-rule-application", "terms": 12}
        monitor.write_text(json.dumps(raw))
        os.utime(monitor, (995, 995))
        status = self.run_worker(timeouts=2)
        self.assertEqual(status["state"], "complete")
        self.assertEqual(
            [row["native_progress"]["age_seconds"] for row in self.running_statuses],
            [5, 10, 15],
        )
        for row in self.running_statuses:
            self.assertEqual(row["native_progress"]["state"], "available")
            self.assertEqual(row["native_progress"]["snapshot"], raw)
        self.assertEqual(monitor.read_text(), json.dumps(raw))

    def test_missing_or_bad_native_monitor_does_not_fail_a_successful_child(self):
        monitor = self.work / "monitor.json"
        for label, payload, expected in (
            ("missing", None, "pending"),
            ("malformed", "{", "unavailable"),
            ("nonobject", "[]", "unavailable"),
            ("oversize", json.dumps({"value": "x" * 100}), "unavailable"),
        ):
            with self.subTest(case=label):
                if payload is None:
                    monitor.unlink(missing_ok=True)
                else:
                    monitor.write_text(payload)
                with patch.object(Gallery, "display_receipt_bytes", 64):
                    status = self.run_worker(timeouts=1)
                self.assertEqual(status["state"], "complete")
                self.assertEqual(status["native_progress"]["state"], expected)
                self.assertEqual(self.wait_timeouts, [5, 5])

    def test_deep_native_monitor_does_not_interrupt_the_child_or_change_the_artifact(
        self,
    ):
        monitor = self.work / "monitor.json"
        raw = '{"layers":' + "[" * 2000 + "0" + "]" * 2000 + "}"
        monitor.write_text(raw)
        decode_json = json.loads

        def decode(payload, *args, **kwargs):
            # Decoder nesting limits differ across supported Python versions.
            if payload == raw.encode():
                raise RecursionError("synthetic decoder nesting limit")
            return decode_json(payload, *args, **kwargs)

        with (
            patch.object(Gallery, "display_receipt_bytes", 10_000),
            patch("benchmark_job.json.loads", side_effect=decode),
        ):
            status = self.run_worker(timeouts=1)
        self.assertEqual(status["state"], "complete")
        self.assertEqual(status["native_progress"]["state"], "unavailable")
        self.assertEqual(self.wait_timeouts, [5, 5])
        self.assertEqual(monitor.read_text(), raw)

    def test_nonzero_exit_retains_failure_and_log(self):
        prefix = self.work / "partial-output"
        prefix.write_bytes(b"synthetic durable prefix")
        status = self.run_worker(result=17, output="synthetic native failure")
        self.assertEqual(status["state"], "failed")
        self.assertEqual(status["exit_code"], 17)
        self.assertIn("synthetic native failure", status["message"])
        self.assertEqual(prefix.read_bytes(), b"synthetic durable prefix")
        self.assertEqual(
            (self.work / "run.log").read_text(), "synthetic native failure"
        )

    def test_spawn_error_is_reported_and_releases_descriptor(self):
        status = self.run_worker(error=OSError("synthetic worker error"))
        self.assertEqual(status["state"], "failed")
        self.assertIn("synthetic worker error", status["message"])


class GalleryWorkerObservationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="gallery-observation-contract-")
        self.addCleanup(temporary.cleanup)
        self.work = Path(temporary.name)

    def test_native_progress_uses_bounded_display_and_file_age(self):
        raw = {"phase": "working", "long_field": "x" * 3000}
        monitor = self.work / "monitor.json"
        monitor.write_text(json.dumps(raw))
        os.utime(monitor, (100, 100))
        with patch("benchmark_job.time", return_value=107):
            progress = benchmark_job.read_native_progress(self.work)
        self.assertEqual(progress["state"], "available")
        self.assertEqual(progress["artifact"], str(monitor))
        self.assertEqual(progress["age_seconds"], 7)
        self.assertEqual(progress["snapshot"], Gallery.display_receipt(raw))
        self.assertLess(len(progress["snapshot"]["long_field"]), len(raw["long_field"]))
        self.assertEqual(json.loads(monitor.read_text()), raw)

    def test_unreadable_native_progress_is_best_effort(self):
        with patch.object(
            Path, "open", side_effect=PermissionError("synthetic denial")
        ):
            progress = benchmark_job.read_native_progress(self.work)
        self.assertEqual(progress["state"], "unavailable")
        self.assertIn("synthetic denial", progress["message"])

    def test_worker_resources_samples_the_requested_pid_and_machine_memory(self):
        proc = self.work / "proc"
        (proc / "321").mkdir(parents=True)
        (proc / "321/status").write_text("Name: synthetic\nVmRSS: 7 kB\nVmHWM: 11 kB\n")
        (proc / "meminfo").write_text("MemTotal: 1000 kB\nMemAvailable: 750 kB\n")
        with patch("benchmark_job.time", return_value=1234):
            resources = benchmark_job.worker_resources(321, proc=proc)
        self.assertEqual(resources["worker_rss_bytes"], 7 * 1024)
        self.assertEqual(resources["worker_peak_rss_bytes"], 11 * 1024)
        self.assertEqual(resources["machine_available_bytes"], 750 * 1024)
        self.assertEqual(resources["machine_total_bytes"], 1000 * 1024)
        self.assertEqual(resources["sampled_unix_seconds"], 1234)
        self.assertIn("/proc", resources["source"])
        self.assertIn("unavailable fields are null", resources["source"])

    def test_missing_proc_resources_are_null_and_do_not_raise(self):
        with patch("benchmark_job.time", return_value=1234):
            resources = benchmark_job.worker_resources(321, proc=self.work / "absent")
        for field in (
            "worker_rss_bytes",
            "worker_peak_rss_bytes",
            "machine_available_bytes",
            "machine_total_bytes",
        ):
            self.assertIsNone(resources[field])
        self.assertEqual(resources["sampled_unix_seconds"], 1234)

    def test_bad_worker_memory_does_not_hide_available_machine_memory(self):
        proc = self.work / "proc"
        (proc / "321").mkdir(parents=True)
        (proc / "321/status").write_text("VmRSS: invalid kB\n")
        (proc / "meminfo").write_text("MemAvailable: 750 kB\n")
        resources = benchmark_job.worker_resources(321, proc=proc)
        self.assertIsNone(resources["worker_rss_bytes"])
        self.assertIsNone(resources["worker_peak_rss_bytes"])
        self.assertEqual(resources["machine_available_bytes"], 750 * 1024)
        self.assertIsNone(resources["machine_total_bytes"])


class GalleryReceiptTests(unittest.TestCase):
    def test_recorded_advisory_does_not_become_an_unknown_or_performance_claim(self):
        mo = SimpleNamespace(
            Html=lambda value: value,
            vstack=lambda values: values,
            accordion=lambda values: values,
            json=lambda value: value,
        )
        with patch.dict(
            os.environ,
            {"VAKINT_GALLERY_EVIDENCE": str(Path(__file__).parent / "data")},
        ):
            record = Gallery.receipt("warning-controls")
            card = Gallery.evidence_card(mo, "warning-controls", "Recorded advisory")
        self.assertFalse(record["performance_claim"])
        self.assertFalse(record["global_campaign_coverage_claim"])
        self.assertFalse(record["unknown_residual_example"])
        self.assertFalse(record["cancelled_residual_example"])
        self.assertEqual(record["outcome"]["unknown_count"], 0)
        rendered = repr(card[:3])
        self.assertIn("reduction continued", rendered)
        self.assertIn("outside the recorded shifted-index entry predicates", rendered)
        self.assertIn("no unknown integrals", rendered)
        self.assertIn("no performance, held-out or speedup claim", rendered)
        self.assertNotIn("Native reduction · seconds", rendered)

    def test_running_snapshot_never_promotes_progress_to_reduction_timing(self):
        mo = SimpleNamespace(
            Html=lambda value: value,
            vstack=lambda values: values,
            accordion=lambda values: values,
            json=lambda value: value,
        )
        receipt = {
            "schema": "vakint.gallery.physical-status.v1",
            "state": "running_at_observation",
            "observed_utc": "2026-10-25T22:55:23Z",
            "running_cache": {"accepted_records": 106699, "committed_records": 106687},
            "reduction_completed": False,
            "final_reduction_seconds": None,
            "speedup": None,
        }
        with (
            patch.object(Gallery, "receipt", return_value=receipt),
            patch.object(Gallery, "evidence_links", return_value="download"),
        ):
            card = Gallery.evidence_card(mo, "synthetic", "Producer snapshot")
        rendered = repr(card[0]) + repr(card[1])
        self.assertIn("106,699", rendered)
        self.assertIn("106,687", rendered)
        self.assertIn("acknowledged durable", rendered)
        self.assertIn("not a completed reduction or a measured speedup", rendered)
        self.assertIn("not a held-out benchmark", rendered)
        self.assertNotIn("Native reduction · seconds", rendered)

    def test_reused_input_null_preparation_is_not_a_current_timing(self):
        tables = []

        def table(rows, **kwargs):
            tables.append(rows)
            return rows

        mo = SimpleNamespace(
            Html=lambda value: value,
            vstack=lambda values: values,
            accordion=lambda values: values,
            json=lambda value: value,
            ui=SimpleNamespace(table=table),
        )
        receipt = {
            "schema": "vakint.qcd-ladder.benchmark.v1",
            "input_mode": "prepared",
            "fresh_native_preparation": False,
            "preparation_seconds": None,
            "preparation": None,
            "input_loading_seconds": 0.125,
            "reduction_seconds": 1.25,
            "end_to_end_seconds": 1.5,
            "reduction_completed": True,
        }
        with (
            patch.object(Gallery, "receipt", return_value=receipt),
            patch.object(Gallery, "evidence_links", return_value="download"),
        ):
            card = Gallery.evidence_card(mo, "synthetic", "Synthetic benchmark")
        self.assertIn("1.250", repr(card))
        self.assertTrue(tables)
        for row in tables[0]:
            self.assertIsNotNone(row["Seconds"])
            self.assertNotIn("Native trace contraction", row["Stage"])
        self.assertIsNone(receipt["preparation_seconds"])


class GalleryClosurePresentationTests(unittest.TestCase):
    def setUp(self):
        self.path = Gallery.directory / "04_pinches_and_cache.py"
        self.tree = ast.parse(self.path.read_text())
        self.cell = next(
            cell
            for cell in self.tree.body
            if isinstance(cell, ast.FunctionDef)
            and any(
                isinstance(node, ast.Name)
                and isinstance(node.ctx, ast.Store)
                and node.id == "closure_view"
                for node in ast.walk(cell)
            )
        )

    def test_tab_selection_does_not_invalidate_native_widget_owner(self):
        self.assertNotIn("closure_view", [arg.arg for arg in self.cell.args.args])
        self.assertFalse(
            any(
                isinstance(node, ast.Attribute)
                and isinstance(node.value, ast.Name)
                and node.value.id == "closure_view"
                for node in ast.walk(self.tree)
            )
        )

    def test_all_native_panels_are_prepared_before_tab_selection(self):
        self.cell.decorator_list = []
        namespace = {}
        # Only the display cell runs, with fake HTML and no native host.
        exec(  # noqa: S102
            compile(
                ast.Module(body=[self.cell], type_ignores=[]), str(self.path), "exec"
            ),
            namespace,
        )
        render = Mock(side_effect=lambda value: ("html", value))
        static_html = Mock(side_effect=lambda value: ("static-html", value))
        mo = SimpleNamespace(
            as_html=render,
            Html=static_html,
            md=lambda value: value,
            vstack=lambda values: values,
            hstack=lambda values: values,
            ui=SimpleNamespace(
                tabs=lambda contents: SimpleNamespace(
                    contents=contents, value="The identity"
                )
            ),
        )
        arguments = {arg.arg: object() for arg in self.cell.args.args}
        arguments.update(
            mo=mo,
            Gallery=SimpleNamespace(
                metrics=lambda *_: "metrics", native_formula=Gallery.native_formula
            ),
            closure_seconds=0.1,
        )
        formula_names = (
            "cancelling_numerator",
            "source_integral",
            "target_integral",
            "source_result",
            "target_result",
            "difference",
        )
        formula_html = {name: f"<math><mi>{name}</mi></math>" for name in formula_names}
        # The native printer's omission metadata must never be stripped.
        formula_html["source_result"] = (
            '<math><mtext data-symbolica-marker="0">⋯</mtext></math>'
            "<p>Expression preview · items 1–25 of 30</p>"
        )
        for name in formula_names:
            arguments[name] = SimpleNamespace(
                _repr_html_=Mock(return_value=formula_html[name]),
                _display_=Mock(side_effect=AssertionError("Interactive formula")),
            )
        (tabs,) = namespace["_"](**arguments)
        self.assertEqual(static_html.call_count, len(formula_names))
        for name in formula_names:
            arguments[name]._repr_html_.assert_called_once_with()
            arguments[name]._display_.assert_not_called()
            static_html.assert_any_call(
                '<div class="vg-native-formula">' + formula_html[name] + "</div>"
            )
        self.assertEqual(render.call_count, 2)
        render.assert_any_call(arguments["source_graph"])
        render.assert_any_call(arguments["target_graph"])
        self.assertEqual(
            list(tabs.contents),
            ["The identity", "Retained-master results", "Independent native graphs"],
        )
        for label in tabs.contents:
            tabs.value = label
            self.assertTrue(tabs.contents[tabs.value])
        self.assertEqual(render.call_count, 2)
        self.assertEqual(static_html.call_count, len(formula_names))


class GalleryNativeFormulaTests(unittest.TestCase):
    def test_missing_native_html_is_an_explicit_error(self):
        for payload in (None, "", "  ", b"not a native HTML string"):
            with self.subTest(payload=payload):
                mo = SimpleNamespace(Html=Mock())
                expression = SimpleNamespace(_repr_html_=Mock(return_value=payload))
                with self.assertRaisesRegex(ValueError, "native HTML"):
                    Gallery.native_formula(mo, expression)
                mo.Html.assert_not_called()

    def test_all_four_tensor_panes_use_the_native_static_protocol(self):
        path = Gallery.directory / "01_graphs_and_numerators.py"
        tree = ast.parse(path.read_text())
        cell = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef)
            and {arg.arg for arg in node.args.args}
            == {
                "Gallery",
                "contraction_seconds",
                "mo",
                "numerator",
                "projected",
                "projector",
                "reduced",
            }
        )
        cell.decorator_list = []
        namespace = {}
        exec(  # noqa: S102 - Execute only this presentation cell with fake native values.
            compile(ast.Module(body=[cell], type_ignores=[]), str(path), "exec"),
            namespace,
        )
        tabs = Mock(side_effect=lambda contents: SimpleNamespace(contents=contents))
        mo = SimpleNamespace(
            Html=Mock(side_effect=lambda text: text),
            vstack=lambda values: values,
            ui=SimpleNamespace(tabs=tabs),
        )
        names = ("reduced", "projector", "numerator", "projected")
        values = {
            name: SimpleNamespace(
                _repr_html_=Mock(return_value=f"<math><mi>{name}</mi></math>"),
                _display_=Mock(side_effect=AssertionError("Interactive widget")),
            )
            for name in names
        }
        namespace["_"](
            Gallery=SimpleNamespace(
                native_formula=Gallery.native_formula,
                metrics=lambda *_: "metrics",
                note=lambda *_: "note",
            ),
            mo=mo,
            contraction_seconds=0.1,
            **values,
        )
        contents = tabs.call_args.args[0]
        self.assertEqual(
            list(contents),
            [
                "Scalar dot products",
                "Projector",
                "Tensor numerator",
                "Projected tensor",
            ],
        )
        for name, payload in zip(names, contents.values(), strict=True):
            values[name]._repr_html_.assert_called_once_with()
            values[name]._display_.assert_not_called()
            self.assertEqual(
                payload,
                f'<div class="vg-native-formula"><math><mi>{name}</mi></math></div>',
            )

    def test_physical_run_cell_exposes_only_zero_momentum(self):
        tree = ast.parse((Gallery.directory / "03_five_loop_uv.py").read_text())
        calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == "Gallery"
            and node.func.attr == "start_benchmark"
        ]
        self.assertEqual(len(calls), 1)
        self.assertIsInstance(calls[0].args[0], ast.Constant)
        self.assertEqual(calls[0].args[0].value, "zero-momentum")
        self.assertFalse(
            any(
                isinstance(node, ast.Name) and node.id == "benchmark_choice"
                for node in ast.walk(tree)
            )
        )


class GalleryFooterTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="gallery-footer-contract-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        environment = patch.dict(
            os.environ, {"VAKINT_GALLERY_EVIDENCE": str(self.root)}
        )
        environment.start()
        self.addCleanup(environment.stop)
        self.mo = SimpleNamespace(Html=lambda value: value)

    def evidence(self, runtime=None):
        receipts = {}
        if runtime is not None:
            payload = json.dumps(runtime).encode()
            (self.root / "runtime.json").write_bytes(payload)
            receipts["runtime"] = {
                "path": "runtime.json",
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
        (self.root / "manifest.json").write_text(
            json.dumps({"schema": 1, "receipts": receipts})
        )

    def test_footer_displays_authenticated_runtime_fields(self):
        self.evidence(
            {"display": [["Compiled Vakint", "native-test"], ["UI", "ui-test"]]}
        )
        footer = Gallery.footer(self.mo)
        self.assertIn("Compiled Vakint: native-test", footer)
        self.assertIn("UI: ui-test", footer)
        self.assertIn("External data and provenance</a>", footer)

    def test_footer_without_runtime_keeps_existing_links(self):
        self.evidence()
        footer = Gallery.footer(self.mo)
        self.assertIn("Native Symbolica expressions", footer)
        self.assertIn("External data and provenance</a></footer>", footer)
        self.assertNotIn("Compiled Vakint:", footer)

    def test_footer_escapes_runtime_labels_and_values(self):
        self.evidence(
            {"display": [["UI<script>", '</footer><script>alert("x")</script>']]}
        )
        footer = Gallery.footer(self.mo)
        self.assertNotIn("<script>", footer)
        self.assertIn("UI&lt;script&gt;:", footer)
        self.assertIn("&lt;/footer&gt;", footer)
        self.assertIn("alert(&quot;x&quot;)", footer)


class GalleryEngineTests(unittest.TestCase):
    """The engine forwards bounded worker cores and the Ordered arithmetic."""

    def engine(self, environment, **arguments):
        vakint = Mock(name="Vakint")
        module = SimpleNamespace(Vakint=vakint, VakintEvaluationMethod=Mock())
        with (
            patch.dict(os.environ, environment, clear=True),
            patch.dict(sys.modules, {"symbolica.community.hepkit.vakint": module}),
            patch("shared_gallery.os.sched_getaffinity", return_value=set(range(64))),
        ):
            Gallery.engine(**arguments)
        return vakint.call_args.kwargs

    def test_defaults_bound_the_cores_and_leave_the_arithmetic_to_vakint(self):
        arguments = self.engine({})
        self.assertEqual(arguments["n_cores"], Gallery.default_max_cores)
        self.assertIsNone(arguments["ordered_weight_arithmetic"])
        self.assertEqual(arguments["integral_normalization_factor"], "FMFTandMATAD")

    def test_server_environment_selects_cores_and_arithmetic(self):
        arguments = self.engine(
            {
                "VAKINT_GALLERY_N_CORES": "4",
                "VAKINT_GALLERY_ARITHMETIC": "modular_sampled",
            }
        )
        self.assertEqual(arguments["n_cores"], 4)
        self.assertEqual(arguments["ordered_weight_arithmetic"], "modular_sampled")

    def test_explicit_arguments_override_the_environment(self):
        arguments = self.engine(
            {
                "VAKINT_GALLERY_N_CORES": "4",
                "VAKINT_GALLERY_ARITHMETIC": "modular_sampled",
            },
            n_cores=1,
            ordered_weight_arithmetic="exact",
        )
        self.assertEqual(arguments["n_cores"], 1)
        self.assertEqual(arguments["ordered_weight_arithmetic"], "exact")

    def test_requested_cores_never_exceed_the_server_ceiling(self):
        for environment, requested, expected in (
            ({}, 64, Gallery.default_max_cores),
            ({"VAKINT_GALLERY_MAX_CORES": "2"}, 8, 2),
            (
                {"VAKINT_GALLERY_MAX_CORES": "2", "VAKINT_GALLERY_N_CORES": "16"},
                None,
                2,
            ),
            ({"VAKINT_GALLERY_MAX_CORES": "invalid"}, 64, Gallery.default_max_cores),
            ({}, 0, 1),
        ):
            with self.subTest(environment=environment, requested=requested):
                arguments = self.engine(environment, n_cores=requested)
                self.assertEqual(arguments["n_cores"], expected)

    def test_ceiling_respects_the_cpus_this_server_may_use(self):
        with (
            patch.dict(os.environ, {"VAKINT_GALLERY_MAX_CORES": "16"}, clear=True),
            patch("shared_gallery.os.sched_getaffinity", return_value={3, 5}),
        ):
            self.assertEqual(Gallery.max_cores(), 2)
            self.assertEqual(Gallery.job_cpus(2), [3, 5])

    def test_options_pass_through(self):
        arguments = self.engine(
            {}, normalization="MSbar", use_dot_product_notation=True
        )
        self.assertEqual(arguments["integral_normalization_factor"], "MSbar")
        self.assertTrue(arguments["use_dot_product_notation"])

    def test_execution_names_the_arithmetic_and_cores(self):
        engine = SimpleNamespace(n_cores=1, ordered_weight_arithmetic="exact")
        self.assertEqual(Gallery.execution(engine), "exact · 1 core")
        engine = SimpleNamespace(n_cores=8, ordered_weight_arithmetic="modular_sampled")
        self.assertEqual(Gallery.execution(engine), "modular_sampled · 8 cores")

    def test_official_data_never_installs_into_the_working_directory(self):
        engine = Mock()
        with tempfile.TemporaryDirectory() as cache:
            with patch.dict(
                os.environ, {"VAKINT_GALLERY_DATA_CACHE": cache}, clear=True
            ):
                Gallery.load_data(engine, 4)
            arguments = engine.load_official_external_data.call_args
            self.assertEqual(arguments.args, (4,))
            self.assertEqual(
                arguments.kwargs["local_artifacts_path"], str(Path(cache).resolve())
            )
            with patch.dict(
                os.environ, {"VAKINT_GALLERY_ARTIFACTS": "/data/checkout"}, clear=True
            ):
                self.assertEqual(Gallery.data_directory(), Path("/data/checkout"))

    def test_five_loop_catalogue_load_needs_the_operator_opt_in(self):
        engine = Mock()
        with tempfile.TemporaryDirectory() as cache:
            environment = {"VAKINT_GALLERY_DATA_CACHE": cache}
            with (
                patch.dict(os.environ, environment, clear=True),
                self.assertRaisesRegex(ValueError, "VAKINT_GALLERY_ALLOW_5L"),
            ):
                Gallery.load_data(engine, 5)
            engine.load_official_external_data.assert_not_called()
            with patch.dict(
                os.environ, {**environment, "VAKINT_GALLERY_ALLOW_5L": "1"}, clear=True
            ):
                Gallery.load_data(engine, 5)
            self.assertEqual(engine.load_official_external_data.call_args.args, (5,))


class FakeProc:
    """A synthetic /proc: {pid: (state, pgid, start_ticks)}."""

    def __init__(self, root, processes):
        self.root = Path(root)
        for pid, (state, pgid, ticks) in processes.items():
            directory = self.root / str(pid)
            directory.mkdir(parents=True)
            # comm with spaces and parentheses must not shift the fields.
            fields = [state, "1", str(pgid), str(pgid)] + ["0"] * 15 + [str(ticks)]
            (directory / "stat").write_text(f"{pid} (odd (name) x) " + " ".join(fields))


class GalleryJobGuardTests(unittest.TestCase):
    """Signals reach only this job's processes, never a recycled PID."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="gallery-job-guards-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.work = self.root / "work" / "job"
        self.work.mkdir(parents=True)
        (self.work / "status.json").write_text(json.dumps({"state": "running"}))

    def proc(self, processes):
        return FakeProc(
            self.root / f"proc-{len(list(self.root.iterdir()))}", processes
        ).root

    def job(self, pid=4242, ticks=1000):
        (self.work / "job.json").write_text(
            json.dumps({"pid": pid, "pgid": pid, "start_ticks": ticks})
        )

    def test_stat_parsing_survives_odd_process_names(self):
        proc = self.proc({7: ("S", 7, 99), 8: ("Z", 7, 120), 9: ("R", 3, 130)})
        self.assertEqual(gallery_jobs.process_start_ticks(7, proc), 99)
        self.assertIsNone(gallery_jobs.process_start_ticks(8, proc), "zombie")
        self.assertIsNone(gallery_jobs.process_start_ticks(10, proc), "absent")
        self.assertEqual(gallery_jobs.group_members(7, proc), {7: 99})

    def test_live_job_signals_its_whole_group(self):
        self.job()
        proc = self.proc({4242: ("S", 4242, 1000), 4243: ("R", 4242, 1001)})
        with (
            patch("gallery_jobs.os.killpg") as killpg,
            patch("gallery_jobs.os.kill") as kill,
        ):
            report = gallery_jobs.cancel(self.work, proc=proc, sleep=lambda _: None)
        self.assertEqual(report["signalled"], [4242, 4243])
        self.assertEqual(killpg.call_args_list[0].args, (4242, signal.SIGTERM))
        self.assertEqual(killpg.call_args_list[-1].args, (4242, signal.SIGKILL))
        kill.assert_not_called()
        status = json.loads((self.work / "status.json").read_text())
        self.assertEqual(status["state"], "cancelled")
        self.assertEqual(
            json.loads((self.work / "cancel-report.json").read_text())["signalled"],
            [4242, 4243],
        )

    def test_recycled_pid_is_never_signalled(self):
        self.job(ticks=1000)
        proc = self.proc({4242: ("S", 4242, 5555), 4243: ("S", 4242, 5556)})
        with (
            patch("gallery_jobs.os.killpg") as killpg,
            patch("gallery_jobs.os.kill") as kill,
        ):
            report = gallery_jobs.cancel(self.work, proc=proc, sleep=lambda _: None)
        killpg.assert_not_called()
        kill.assert_not_called()
        self.assertEqual(report["signalled"], [])
        self.assertFalse(gallery_jobs.is_running(self.work, proc))

    def test_orphans_of_a_dead_supervisor_are_signalled_individually(self):
        self.job(ticks=1000)
        # Supervisor gone; one member started after it (ours), one before (not).
        proc = self.proc({4243: ("S", 4242, 1001), 4244: ("S", 4242, 900)})
        with (
            patch("gallery_jobs.os.killpg") as killpg,
            patch("gallery_jobs.os.kill") as kill,
        ):
            report = gallery_jobs.cancel(self.work, proc=proc, sleep=lambda _: None)
        killpg.assert_not_called()
        self.assertEqual(report["signalled"], [4243])
        self.assertEqual(
            {call.args for call in kill.call_args_list},
            {(4243, signal.SIGTERM), (4243, signal.SIGKILL)},
        )

    def test_own_or_init_process_group_is_refused(self):
        for pgid in (1, 0, os.getpgrp()):
            with self.subTest(pgid=pgid):
                self.job(pid=pgid, ticks=1)
                proc = self.proc({pgid: ("S", pgid, 1)} if pgid > 0 else {})
                with (
                    patch("gallery_jobs.os.killpg") as killpg,
                    patch("gallery_jobs.os.kill") as kill,
                ):
                    gallery_jobs.cancel(self.work, proc=proc, sleep=lambda _: None)
                killpg.assert_not_called()
                kill.assert_not_called()

    def test_dead_supervisor_without_final_status_is_lost(self):
        self.job()
        proc = self.proc({})
        self.assertEqual(gallery_jobs.status(self.work, proc)["state"], "lost")
        (self.work / "status.json").write_text(json.dumps({"state": "complete"}))
        self.assertEqual(gallery_jobs.status(self.work, proc)["state"], "complete")

    def test_active_job_is_found_after_a_reload(self):
        self.job()
        proc = self.proc({4242: ("S", 4242, 1000)})
        (self.work / "request.json").write_text('{"kind": "ladder"}')
        self.assertEqual(
            gallery_jobs.active_job(self.work.parent, proc), str(self.work)
        )
        self.assertEqual(
            gallery_jobs.active_job(self.work.parent, proc, kind="ladder"),
            str(self.work),
        )
        self.assertIsNone(
            gallery_jobs.active_job(self.work.parent, proc, kind="five-loop")
        )
        proc = self.proc({4242: ("S", 4242, 77)})
        self.assertIsNone(gallery_jobs.active_job(self.work.parent, proc))


class GalleryJobProcessTests(unittest.TestCase):
    """Real supervisor processes running a short ``sleep`` command."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="gallery-job-process-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.sleeper = [sys.executable, "-c", "import time; time.sleep(120)"]

    def launch(self, **request):
        return Path(
            gallery_jobs.launch(
                self.root, "sleep", {"command": self.sleeper, "nice": 5, **request}
            )
        )

    def wait_for(self, work, predicate, limit=20.0):
        deadline = time.monotonic() + limit
        while time.monotonic() < deadline:
            status = json.loads((work / "status.json").read_text())
            if predicate(status):
                return status
            time.sleep(0.05)
        self.fail(f"timed out; last status {status}")

    def assert_job_gone(self, work):
        job = json.loads((work / "job.json").read_text())
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline:
            gallery_jobs._reap()
            if not gallery_jobs.group_members(job["pgid"]):
                break
            time.sleep(0.05)
        self.assertEqual(gallery_jobs.group_members(job["pgid"]), {})
        # The work root's lock is free again: a new job may start.
        with (self.root / gallery_jobs.LOCK_NAME).open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def test_cancel_stops_supervisor_and_command(self):
        work = self.launch(environment={"GALLERY_TEST": "1"})
        running = self.wait_for(work, lambda status: status.get("worker_pid"))
        child = running["worker_pid"]
        self.assertTrue(gallery_jobs.is_running(work))
        self.assertEqual(gallery_jobs.active_job(self.root), str(work))
        with self.assertRaisesRegex(RuntimeError, "already running"):
            self.launch()
        report = gallery_jobs.cancel(work)
        self.assertEqual(report["state"], "cancelled")
        self.assertIn(child, report["signalled"])
        self.assertEqual(report["survivors"], [])
        self.assertIsNone(gallery_jobs.process_start_ticks(child))
        self.assertEqual(gallery_jobs.status(work)["state"], "cancelled")
        self.assert_job_gone(work)
        self.assertIsNone(gallery_jobs.active_job(self.root))

    @unittest.skipUnless(importlib.util.find_spec("marimo"), "needs marimo")
    def test_headless_page_never_cancels_another_sessions_job(self):
        work = self.launch(kind="ladder")
        self.wait_for(work, lambda status: status.get("worker_pid"))
        environment = {**os.environ, "VAKINT_GALLERY_WORK": str(self.root)}
        page = Gallery.directory / "05_ladder_reductions.py"
        for arguments in (
            ["--case", "3L-zero-momentum", "--n-cores", "1", "--cancel-after", "0"],
            ["--case", "no-such-case", "--cancel-after", "0"],
        ):
            with self.subTest(arguments=arguments):
                finished = subprocess.run(
                    [sys.executable, str(page), *arguments],
                    capture_output=True,
                    text=True,
                    env=environment,
                    timeout=120,
                    check=False,
                )
                self.assertEqual(finished.returncode, 2, finished.stderr)
                summary = json.loads(finished.stdout)
                self.assertEqual(summary["state"], "not started")
                self.assertFalse(summary["ok"])
                self.assertTrue(finished.stderr.strip())
        self.assertTrue(gallery_jobs.is_running(work))
        self.assertFalse((work / "cancel-requested.json").exists())
        self.assertEqual(len(gallery_jobs.jobs(self.root)), 1)
        gallery_jobs.cancel(work)
        self.assert_job_gone(work)

    def test_deadline_stops_the_command(self):
        work = self.launch(max_seconds=0.5)
        status = self.wait_for(
            work, lambda status: status["state"] in gallery_jobs.FINAL_STATES
        )
        self.assertEqual(status["state"], "timed_out")
        self.assertIn("deadline", status["message"])
        self.assert_job_gone(work)

    def test_command_dies_with_a_killed_supervisor(self):
        work = self.launch()
        running = self.wait_for(work, lambda status: status.get("worker_pid"))
        job = json.loads((work / "job.json").read_text())
        os.kill(job["pid"], signal.SIGKILL)
        self.assert_job_gone(work)
        self.assertIsNone(gallery_jobs.process_start_ticks(running["worker_pid"]))
        self.assertEqual(gallery_jobs.status(work)["state"], "lost")

    def test_command_gets_bounded_cpus_and_environment(self):
        allowed = sorted(os.sched_getaffinity(0))[:1]
        probe = (
            "import json, os, sys; json.dump({'cpus': sorted(os.sched_getaffinity(0)),"
            " 'nice': os.nice(0), 'value': os.environ.get('GALLERY_TEST')},"
            " open(sys.argv[1], 'w'))"
        )
        output = self.root / "probe.json"
        work = Path(
            gallery_jobs.launch(
                self.root,
                "probe",
                {
                    "command": [sys.executable, "-c", probe, str(output)],
                    "cpus": allowed,
                    "nice": 3,
                    "environment": {"GALLERY_TEST": "bounded"},
                },
            )
        )
        status = self.wait_for(
            work, lambda status: status["state"] in gallery_jobs.FINAL_STATES
        )
        self.assertEqual(status["state"], "complete")
        observed = json.loads(output.read_text())
        self.assertEqual(observed["cpus"], allowed)
        self.assertGreaterEqual(observed["nice"], 3)
        self.assertEqual(observed["value"], "bounded")
        self.assert_job_gone(work)


class GalleryLadderTests(unittest.TestCase):
    """The 3L/4L cases: shipped inputs, bounded launches and recorded hashes."""

    # Retained-result sha256 prefixes of the rung fixtures (recorded with vakint
    # data v0.1.2, reproduced with v0.2.0).
    recorded = {
        "3L-zero-momentum": "2f48233c",
        "3L-quadratic": "5acbeebe",
        "4L-zero-momentum": "b93080b3",
        "4L-quadratic": "6986f629",
    }

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="gallery-ladder-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        environment = patch.dict(
            os.environ,
            {
                "VAKINT_GALLERY_WORK": str(self.root / "work"),
                "VAKINT_GALLERY_DATA_CACHE": str(self.root / "data"),
            },
            clear=True,
        )
        environment.start()
        self.addCleanup(environment.stop)

    def test_manifest_inputs_and_recorded_results_are_consistent(self):
        cases = Gallery.ladder_cases()
        self.assertEqual(set(cases), set(self.recorded))
        resources = Gallery.directory.parents[1] / "tests/resources/ordered_rungs"
        for key, case in cases.items():
            with self.subTest(case=key):
                self.assertTrue(case["retained_sha256"].startswith(self.recorded[key]))
                self.assertEqual(case["loop_count"], int(key[0]))
                if case["input"]["shipped"]:
                    payload = Gallery.ladder_input(key).read_bytes()
                    self.assertEqual(len(payload), case["input"]["bytes"])
                    self.assertEqual(
                        hashlib.sha256(payload).hexdigest(), case["input"]["sha256"]
                    )
                reference = resources / f"{key}.json"
                if reference.is_file():
                    reference = json.loads(reference.read_text())
                    self.assertEqual(
                        reference["source_input_sha256"], case["input"]["sha256"]
                    )
                    self.assertEqual(
                        reference["evaluated_input_sha256"],
                        case["evaluated_input_sha256"],
                    )
                    self.assertEqual(
                        reference["dimension_binding"]["target"],
                        case["dimension_binding"]["target"],
                    )
        defaults = [key for key, case in cases.items() if case["default"]]
        self.assertEqual(
            defaults, ["3L-zero-momentum", "3L-quadratic", "4L-zero-momentum"]
        )

    def test_unshipped_input_needs_server_configuration(self):
        self.assertFalse(Gallery.ladder_case_available("4L-quadratic"))
        with (
            patch("gallery_jobs.subprocess.Popen") as popen,
            self.assertRaisesRegex(ValueError, "VAKINT_GALLERY_LADDER_INPUTS"),
        ):
            Gallery.start_reduction("4L-quadratic")
        popen.assert_not_called()
        inputs = self.root / "inputs"
        inputs.mkdir()
        (inputs / "4L-quadratic.vakint-input.symbolica").write_text("x")
        with patch.dict(os.environ, {"VAKINT_GALLERY_LADDER_INPUTS": str(inputs)}):
            self.assertTrue(Gallery.ladder_case_available("4L-quadratic"))

    def test_invalid_requests_launch_nothing(self):
        for arguments, message in (
            ({"case": "5L-zero-momentum"}, "Unknown ladder case"),
            ({"case": "../3L-zero-momentum"}, "Unknown ladder case"),
            ({"case": "3L-zero-momentum", "arithmetic": "dfs"}, "arithmetic"),
            ({"case": "3L-zero-momentum", "max_seconds": 0}, "deadline"),
        ):
            with (
                self.subTest(arguments=arguments),
                patch("gallery_jobs.subprocess.Popen") as popen,
                self.assertRaisesRegex(ValueError, message),
            ):
                Gallery.start_reduction(arguments.pop("case"), **arguments)
            popen.assert_not_called()
        self.assertFalse((self.root / "work").exists())

    def test_launch_is_bounded_and_runs_the_reduce_worker(self):
        with (
            patch.dict(
                os.environ,
                {"VAKINT_GALLERY_MAX_CORES": "2", "VAKINT_GALLERY_MAX_SECONDS": "90"},
            ),
            patch("gallery_jobs.subprocess.Popen") as popen,
        ):
            work = Path(
                Gallery.start_reduction(
                    "4L-zero-momentum", n_cores=64, max_seconds=10_000
                )
            )
        request = json.loads((work / "request.json").read_text())
        cores = min(2, len(os.sched_getaffinity(0)))
        self.assertEqual(request["n_cores"], cores)
        self.assertEqual(request["environment"]["VAKINT_N_CORES"], str(cores))
        self.assertEqual(len(request["cpus"]), cores)
        self.assertEqual(request["max_seconds"], 90)
        self.assertEqual(request["kind"], "ladder")
        self.assertEqual(
            request["command"],
            [sys.executable, str(Gallery.directory / "reduce_worker.py"), str(work)],
        )
        spec = request["spec"]
        self.assertEqual(spec["case"], "4L-zero-momentum")
        self.assertEqual(spec["arithmetic"], "modular_polynomial")
        self.assertEqual(spec["n_cores"], cores)
        # 4L data installs into the gallery cache, never the working directory.
        self.assertEqual(spec["data_directory"], str((self.root / "data").resolve()))
        self.assertTrue(spec["ordered_cache_directory"].startswith(str(self.root)))
        self.assertTrue(popen.call_args.kwargs["start_new_session"])
        self.assertEqual(work.parent, (self.root / "work").resolve())

    def test_three_loop_cases_use_the_built_in_registry(self):
        with patch("gallery_jobs.subprocess.Popen"):
            work = Path(
                Gallery.start_reduction("3L-zero-momentum", ordered_cache=False)
            )
        spec = json.loads((work / "request.json").read_text())["spec"]
        self.assertIsNone(spec["data_directory"])
        self.assertIsNone(spec["ordered_cache_directory"])

    def test_ladder_dot_reproduces_the_five_loop_fixture(self):
        pattern = r'(\w+) -> (\w+) \[particle="(\w+)"'
        fixture = Gallery.directory.parents[1] / "tests/fixtures/qcd_ladder/ladder.dot"
        if fixture.is_file():
            self.assertEqual(
                set(re.findall(pattern, Gallery.ladder_dot(5))),
                set(re.findall(pattern, fixture.read_text())),
            )
        for loops in range(1, 6):
            edges = re.findall(pattern, Gallery.ladder_dot(loops))
            self.assertEqual(len(edges), 3 * loops + 1)
            self.assertEqual(sum(kind == "g" for *_, kind in edges), loops + 1)
        with self.assertRaises(ValueError):
            Gallery.ladder_dot(0)

    def test_job_status_overlays_live_stage_and_progress(self):
        work = self.root / "work" / "job"
        work.mkdir(parents=True)
        pid = os.getpid()
        (work / "job.json").write_text(
            json.dumps(
                {
                    "pid": pid,
                    "pgid": pid,
                    "start_ticks": gallery_jobs.process_start_ticks(pid),
                }
            )
        )
        (work / "status.json").write_text(
            json.dumps({"state": "running", "worker_pid": pid, "stage": None})
        )
        (work / "stage.json").write_text(json.dumps({"stage": "Ordered reduction"}))
        (work / "monitor.json").write_text(json.dumps({"caller_phase": "reduction"}))
        status = Gallery.job_status(work)
        self.assertEqual(status["stage"]["stage"], "Ordered reduction")
        self.assertEqual(status["native_progress"]["state"], "available")
        self.assertIsNotNone(status["resources"]["worker_rss_bytes"])
        (work / "status.json").write_text(json.dumps({"state": "complete"}))
        self.assertEqual(Gallery.job_status(work), {"state": "complete"})
        with self.assertRaisesRegex(ValueError, "escapes"):
            Gallery.job_status(self.root)

    def test_wait_job_cancels_like_the_button(self):
        states = iter(["running", "running", "cancelled"])
        with (
            patch.object(
                Gallery, "job_status", side_effect=lambda _: {"state": next(states)}
            ),
            patch.object(
                Gallery, "cancel_job", return_value={"signalled": [1], "survivors": []}
            ) as cancel,
        ):
            status = Gallery.wait_job("job", cancel_after=0, poll_seconds=0)
        cancel.assert_called_once_with("job")
        self.assertEqual(status["state"], "cancelled")
        self.assertEqual(status["cancel_report"]["survivors"], [])

    def test_wait_job_cancels_its_job_on_ctrl_c(self):
        def status(_):
            if not cancel.called:
                raise KeyboardInterrupt
            return {"state": "cancelled"}

        with (
            patch.object(Gallery, "job_status", side_effect=status),
            patch.object(
                Gallery, "cancel_job", return_value={"signalled": [1], "survivors": []}
            ) as cancel,
        ):
            with self.assertRaises(KeyboardInterrupt):
                Gallery.wait_job("job", poll_seconds=0)
            cancel.assert_not_called()
            status = Gallery.wait_job("job", poll_seconds=0, cancel_on_interrupt=True)
        cancel.assert_called_once_with("job")
        self.assertEqual(status["state"], "cancelled")
        self.assertTrue(status["interrupted"])

    # The ladder page's default run form.
    form = MappingProxyType(
        {
            "case": "3L-zero-momentum",
            "arithmetic": "modular_polynomial",
            "n_cores": 2,
            "max_seconds": 600,
            "ordered_cache": True,
        }
    )

    def test_headless_request_takes_the_command_line_and_rejects_typos(self):
        request = Gallery.headless_request(
            {"case": "4L-zero-momentum", "n-cores": 4, "cancel-after": 3}, self.form
        )
        self.assertEqual(
            request,
            {
                "case": "4L-zero-momentum",
                "n_cores": 4,
                "arithmetic": "modular_polynomial",
                "max_seconds": 600.0,
                "ordered_cache": True,
                "cancel_after": 3.0,
            },
        )
        request = Gallery.headless_request(
            {"case": "3L-quadratic", "arithmetic": "exact", "max-seconds": 90},
            self.form,
        )
        self.assertEqual(request["arithmetic"], "exact")
        self.assertEqual(request["max_seconds"], 90.0)
        self.assertIsNone(request["cancel_after"])
        for arguments, message in (
            ({}, "Name one case"),
            ({"case": ""}, "Name one case"),
            ({"case": "3L-quadratic", "ncores": 2}, "Unknown option --ncores"),
            ({"case": "3L-quadratic", "n-cores": 2.5}, "whole number"),
            ({"case": "3L-quadratic", "n-cores": [1, 2]}, "n-cores"),
            ({"case": "3L-quadratic", "max-seconds": "soon"}, "max-seconds"),
            ({"case": "3L-quadratic", "cancel-after": -1}, "non-negative"),
        ):
            with (
                self.subTest(arguments=arguments),
                self.assertRaisesRegex(ValueError, message),
            ):
                Gallery.headless_request(arguments, self.form)

    def finished_job(self, receipt=None, **status):
        work = Gallery.work_root() / f"3L-zero-momentum-{uuid.uuid4().hex[:12]}"
        work.mkdir(parents=True)
        (work / "request.json").write_text(json.dumps({"kind": "ladder"}))
        (work / "status.json").write_text(json.dumps(status))
        if receipt is not None:
            (work / "benchmark.json").write_text(json.dumps(receipt))
        return str(work)

    def receipt(self, matches=True, closed=True):
        return {
            "schema": "vakint.gallery.ladder-run.v1",
            "case": "3L-zero-momentum",
            "execution": "modular_polynomial · 2 cores",
            "retained": {
                "sha256": "2f48233c" if matches else "deadbeef",
                "recorded_sha256": "2f48233c",
                "matches_recorded": matches,
            },
            "catalogue_closure": {"master_count": 14, "closed": closed},
            "seconds": {"Ordered reduction": 0.2, "total": 2.1, "other": 1},
        }

    def test_headless_summary_passes_only_a_reproduced_or_requested_outcome(self):
        plain = {"case": "3L-zero-momentum", "cancel_after": None}
        cancelling = {**plain, "cancel_after": 3.0}
        clean = {"signalled": [10, 11], "survivors": []}
        cases = (
            ("reproduced", self.receipt(), {"state": "complete"}, plain, 0),
            (
                "other hash",
                self.receipt(matches=False),
                {"state": "complete"},
                plain,
                1,
            ),
            ("not closed", self.receipt(closed=False), {"state": "complete"}, plain, 1),
            ("no receipt", None, {"state": "complete"}, plain, 1),
            (
                "failed",
                {"schema": "vakint.gallery.ladder-run.v1", "error": "ValueError: sha"},
                {"state": "failed"},
                plain,
                1,
            ),
            ("timed out", None, {"state": "timed_out"}, plain, 1),
            ("lost", None, {"state": "lost"}, plain, 1),
            (
                "cancelled as asked",
                None,
                {"state": "cancelled", "cancel_report": clean},
                cancelling,
                0,
            ),
            (
                "cancelled, not asked",
                None,
                {"state": "cancelled", "cancel_report": clean},
                plain,
                1,
            ),
            (
                "survivor",
                None,
                {
                    "state": "cancelled",
                    "cancel_report": {"signalled": [10], "survivors": [10]},
                },
                cancelling,
                1,
            ),
            (
                "ctrl-c",
                None,
                {"state": "cancelled", "cancel_report": clean, "interrupted": True},
                plain,
                130,
            ),
        )
        for name, receipt, status, request, code in cases:
            with self.subTest(name):
                work = self.finished_job(receipt, **status)
                summary, exit_code = Gallery.headless_summary(work, status, request)
                self.assertEqual(exit_code, code)
                self.assertEqual(summary["ok"], code == 0)
                self.assertEqual(summary["job"], work)
                self.assertEqual(summary["state"], status["state"])
        summary, _ = Gallery.headless_summary(
            self.finished_job(self.receipt(), state="complete"),
            {"state": "complete"},
            plain,
        )
        self.assertEqual(summary["retained_sha256"], "2f48233c")
        self.assertTrue(summary["matches_recorded"])
        self.assertEqual(
            summary["catalogue_closure"], {"master_count": 14, "closed": True}
        )
        self.assertEqual(summary["seconds"], {"Ordered reduction": 0.2, "total": 2.1})

    def test_headless_report_prints_one_line_and_fails_a_different_hash(self):
        request = {"case": "3L-zero-momentum", "cancel_after": None}
        work = self.finished_job(self.receipt(), state="complete")
        output, errors = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
            Gallery.headless_report(work, {"state": "complete"}, request)
        self.assertTrue(json.loads(output.getvalue())["ok"])
        self.assertEqual(output.getvalue().count("\n"), 1)
        self.assertEqual(errors.getvalue(), "")
        work = self.finished_job(self.receipt(matches=False), state="complete")
        output, errors = io.StringIO(), io.StringIO()
        with (
            contextlib.redirect_stdout(output),
            contextlib.redirect_stderr(errors),
            self.assertRaises(SystemExit) as stop,
        ):
            Gallery.headless_report(work, {"state": "complete"}, request)
        self.assertEqual(stop.exception.code, 1)
        self.assertFalse(json.loads(output.getvalue())["matches_recorded"])
        self.assertIn("differs from the recorded fixture", errors.getvalue())
        output = io.StringIO()
        with (
            contextlib.redirect_stdout(output),
            contextlib.redirect_stderr(io.StringIO()),
            self.assertRaises(SystemExit) as stop,
        ):
            Gallery.headless_exit({"case": "nope"}, "Unknown ladder case: nope")
        self.assertEqual(stop.exception.code, 2)
        self.assertEqual(
            json.loads(output.getvalue()),
            {
                "ok": False,
                "state": "not started",
                "case": "nope",
                "error": "Unknown ladder case: nope",
            },
        )

    def test_ladder_job_reopens_only_ladder_jobs_of_this_work_root(self):
        work = self.finished_job(self.receipt(), state="complete")
        name = Path(work).name
        self.assertEqual(Gallery.ladder_job(name), work)
        other = Gallery.work_root() / "zero-momentum-0123"
        other.mkdir()
        (other / "request.json").write_text(json.dumps({"kind": "five-loop"}))
        (other / "status.json").write_text(json.dumps({"state": "complete"}))
        outside = self.root / "outside"
        outside.mkdir()
        (outside / "request.json").write_text(json.dumps({"kind": "ladder"}))
        (outside / "status.json").write_text("{}")
        for value in (
            None,
            "",
            ["a", "b"],
            other.name,
            "../outside",
            "..",
            str(outside),
            "missing-job",
        ):
            with self.subTest(value=value):
                self.assertIsNone(Gallery.ladder_job(value))

    def test_progress_rows_summarize_the_native_snapshot(self):
        status = {
            "state": "running",
            "elapsed_seconds": 12.5,
            "stage": {"stage": "Ordered reduction"},
            "resources": {"worker_rss_bytes": 2 * 1024**3},
            "native_progress": {
                "snapshot": {
                    "monitor": {"phase": "ordered_reduction"},
                    "ordered_reduction": {
                        "completed_phases": 3,
                        "total_phases": 9,
                        "phase_kind": "Routing",
                        "frontier_keys": 1234,
                        "processed_keys": 5,
                    },
                }
            },
        }
        rows = dict(Gallery.progress_rows(status))
        self.assertEqual(rows["Stage"], "Ordered reduction")
        self.assertEqual(rows["Ordered phases"], "3 / 9 (Routing)")
        self.assertEqual(rows["Frontier keys"], "1,234")
        self.assertEqual(rows["Worker RSS"], "2.00 GiB")
        self.assertEqual(dict(Gallery.progress_rows({}))["State"], "unknown")


class GalleryExploreOnlyTests(unittest.TestCase):
    """No page reduces at five loops, or contracts 4L/5L numerators, by itself."""

    def cell(self, notebook, predicate):
        path = Gallery.directory / notebook
        tree = ast.parse(path.read_text())
        cell = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and predicate(node)
        )
        cell.decorator_list = []
        namespace = {}
        exec(  # noqa: S102 - Execute only this cell with fake values.
            compile(ast.Module(body=[cell], type_ignores=[]), str(path), "exec"),
            namespace,
        )
        return namespace["_"], {arg.arg for arg in cell.args.args}

    class Stop(Exception):
        pass

    def fake_mo(self, mode="run", arguments=None, query=None):
        def stop(condition, output=None):
            if condition:
                raise self.Stop(output)

        query = {} if query is None else query
        return SimpleNamespace(
            stop=stop,
            app_meta=lambda: SimpleNamespace(mode=mode),
            cli_args=lambda: dict(arguments or {}),
            query_params=lambda: query,
            callout=lambda *args, **kwargs: ("callout", args),
            md=lambda text: text,
        )

    def test_five_loop_numerator_is_never_contracted_in_the_page(self):
        run, arguments = self.cell(
            "01_graphs_and_numerators.py",
            lambda node: any(
                isinstance(call, ast.Attribute) and call.attr == "apply_feynman_rules"
                for call in ast.walk(node)
            ),
        )
        self.assertIn("graph_choice", arguments)
        for loops in (5, 6):
            diagram = Mock()
            with self.subTest(loops=loops), self.assertRaises(self.Stop):
                run(
                    S=Mock(),
                    TensorExpression=Mock(),
                    contract_button=SimpleNamespace(value=True),
                    diagram=diagram,
                    graph_choice=SimpleNamespace(value=loops),
                    mo=self.fake_mo(),
                    model=Mock(),
                    perf_counter=lambda: 0.0,
                )
            diagram.apply_feynman_rules.assert_not_called()

    def test_five_loop_page_offers_no_enabled_run_without_opt_in(self):
        source = (Gallery.directory / "03_five_loop_uv.py").read_text()
        self.assertIn("disabled=not _enabled", source)
        self.assertIn("Gallery.five_loop_runs_enabled()", source)
        heavy = {
            "evaluate",
            "evaluate_integral",
            "evaluate_rustred_symbolic",
            "load_data",
            "load_official_external_data",
            "load_rustred_application",
            "apply_feynman_rules",
        }
        calls = [
            node.func.attr
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr in heavy
        ]
        self.assertEqual(calls, [])

    def test_ladder_page_starts_nothing_until_asked(self):
        run, arguments = self.cell(
            "05_ladder_reductions.py",
            lambda node: any(
                isinstance(call, ast.Attribute) and call.attr == "start_reduction"
                for call in ast.walk(node)
            ),
        )
        gallery = SimpleNamespace(start_reduction=Mock(return_value="job"))
        form = SimpleNamespace(
            value={
                "case": "3L-zero-momentum",
                "arithmetic": "exact",
                "n_cores": 2,
                "max_seconds": 60,
                "ordered_cache": False,
            }
        )
        headless, request, started = run(
            Gallery=gallery,
            mo=self.fake_mo(),
            run_button=SimpleNamespace(value=False),
            run_form=form,
        )
        self.assertEqual((headless, request, started), (False, None, None))
        gallery.start_reduction.assert_not_called()
        headless, request, started = run(
            Gallery=gallery,
            mo=self.fake_mo(),
            run_button=SimpleNamespace(value=True),
            run_form=form,
        )
        self.assertEqual(started, "job")
        self.assertIsNone(request["cancel_after"])
        gallery.start_reduction.assert_called_once_with(
            "3L-zero-momentum",
            n_cores=2,
            arithmetic="exact",
            max_seconds=60.0,
            ordered_cache=False,
        )

    def ladder_cells(self):
        run, _ = self.cell(
            "05_ladder_reductions.py",
            lambda node: any(
                isinstance(call, ast.Attribute) and call.attr == "start_reduction"
                for call in ast.walk(node)
            ),
        )
        watch, _ = self.cell(
            "05_ladder_reductions.py",
            lambda node: any(
                isinstance(target, ast.Name) and target.id == "watched_job"
                for assign in ast.walk(node)
                if isinstance(assign, ast.Assign)
                for target in assign.targets
            ),
        )
        return run, watch

    def test_refused_headless_start_exits_and_never_touches_another_job(self):
        run, watch = self.ladder_cells()
        form = SimpleNamespace(value=GalleryLadderTests.form)
        for arguments, error in (
            ({"case": "3L-zero-momentum"}, gallery_jobs.JobBusyError("/work/other")),
            ({"case": "3L-zero-momentum"}, ValueError("Prepared input is missing")),
            ({"case": "no-such-case"}, ValueError("Unknown ladder case")),
            ({}, None),
        ):
            gallery = SimpleNamespace(
                headless_request=Gallery.headless_request,
                headless_exit=Gallery.headless_exit,
                start_reduction=Mock(side_effect=error),
                active_job=Mock(return_value="/work/other"),
                cancel_job=Mock(),
                wait_job=Mock(),
            )
            output = io.StringIO()
            with (
                self.subTest(arguments=arguments),
                contextlib.redirect_stdout(output),
                contextlib.redirect_stderr(io.StringIO()),
                self.assertRaises(SystemExit) as stop,
            ):
                run(
                    Gallery=gallery,
                    mo=self.fake_mo("script", arguments),
                    run_button=SimpleNamespace(value=False),
                    run_form=form,
                )
            self.assertEqual(stop.exception.code, 2)
            summary = json.loads(output.getvalue())
            self.assertFalse(summary["ok"])
            self.assertEqual(summary["state"], "not started")
            gallery.active_job.assert_not_called()
            gallery.cancel_job.assert_not_called()
            gallery.wait_job.assert_not_called()
        # Script mode watches only the job it started, never a running one.
        gallery = SimpleNamespace(active_job=Mock(return_value="/work/other"))
        shown = Mock(return_value=None)
        self.assertEqual(
            watch(
                Gallery=gallery,
                get_shown_job=shown,
                headless=True,
                mo=self.fake_mo("script"),
                set_shown_job=Mock(),
                started_job=None,
            ),
            (None,),
        )
        gallery.active_job.assert_not_called()

    def test_ladder_page_keeps_showing_its_last_run(self):
        _, watch = self.ladder_cells()
        gallery = SimpleNamespace(active_job=Mock(return_value=None))
        shown = {"job": None}
        query = {}

        def show(job=None):
            return shown["job"]

        def set_shown(job):
            shown["job"] = job

        def watched(started):
            return watch(
                Gallery=gallery,
                get_shown_job=show,
                headless=False,
                mo=self.fake_mo(query=query),
                set_shown_job=set_shown,
                started_job=started,
            )[0]

        self.assertIsNone(watched(None))
        self.assertEqual(
            watched("/work/3L-zero-momentum-1"), "/work/3L-zero-momentum-1"
        )
        self.assertEqual(query, {"job": "3L-zero-momentum-1"})
        # A control changes after the run ended: the run cell re-runs without
        # starting anything and nothing is running any more.
        self.assertEqual(watched(None), "/work/3L-zero-momentum-1")
        gallery.active_job.assert_called_once()
        # A fresh page (no last run) re-attaches to a run still going.
        shown["job"] = None
        gallery.active_job.return_value = "/work/4L-zero-momentum-2"
        self.assertEqual(watched(None), "/work/4L-zero-momentum-2")
        self.assertEqual(shown["job"], "/work/4L-zero-momentum-2")
        self.assertEqual(query, {"job": "4L-zero-momentum-2"})

    def test_inspection_page_offers_the_five_loop_load_only_with_opt_in(self):
        choice, _ = self.cell(
            "02_constants_and_masters.py",
            lambda node: any(
                isinstance(call, ast.Attribute) and call.attr == "radio"
                for call in ast.walk(node)
            ),
        )
        for opt_in, loops in (({}, [4]), ({"VAKINT_GALLERY_ALLOW_5L": "1"}, [4, 5])):
            radio = Mock()
            mo = SimpleNamespace(
                ui=SimpleNamespace(radio=radio, run_button=Mock()),
                vstack=Mock(),
                md=lambda text: text,
                Html=lambda value: value,
            )
            with (
                self.subTest(opt_in=opt_in),
                tempfile.TemporaryDirectory() as cache,
                patch.dict(
                    os.environ,
                    {"VAKINT_GALLERY_DATA_CACHE": cache, **opt_in},
                    clear=True,
                ),
            ):
                choice(Gallery=Gallery, mo=mo)
            self.assertEqual(list(radio.call_args.args[0].values()), loops)

    def test_optional_missing_evidence_renders_nothing(self):
        mo = SimpleNamespace(md=lambda text: ("md", text), Html=lambda value: value)
        with patch.object(Gallery, "receipt", return_value=None):
            self.assertEqual(
                Gallery.evidence_card(mo, "closure", "Suite", optional=True),
                ("md", ""),
            )
            self.assertIn(
                "awaiting validated measurements",
                Gallery.evidence_card(mo, "closure", "Suite"),
            )

    def test_five_loop_plan_reads_only_the_selection_index(self):
        with tempfile.TemporaryDirectory() as root:
            ibps = Path(root) / "5L/ibps"
            ibps.mkdir(parents=True)
            (Path(root) / "5L/manifest.json").write_text('{"version": "v9"}')
            (ibps / "selection.json").write_text(
                json.dumps(
                    {
                        "owners": [
                            {"mask": "0111", "bytes": 1024**2},
                            {"mask": "0011", "bytes": 3 * 1024**2},
                            {"mask": "1011", "bytes": 1024**2},
                        ],
                        "route_record_count": 5,
                        "total_bundle_bytes": 5 * 1024**2,
                    }
                )
            )
            plan = Gallery.five_loop_plan(root)
        self.assertEqual(plan["version"], "v9")
        self.assertEqual(plan["owners"], 3)
        self.assertEqual(plan["route_records"], 5)
        self.assertEqual(
            [(row["owner support (lines)"], row["owners"]) for row in plan["supports"]],
            [(3, 2), (2, 1)],
        )
        with tempfile.TemporaryDirectory() as empty:
            self.assertIsNone(Gallery.five_loop_plan(empty))


if __name__ == "__main__":
    unittest.main()
