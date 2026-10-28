"""Gallery launch and display contracts; no native host or subprocesses run."""

import ast
import fcntl
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import benchmark_job
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
            "VAKINT_GALLERY_CACHE": str(self.root / "NOT-A-CACHE"),
        }
        environment_patch = patch.dict(os.environ, self.environment, clear=True)
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
        for name in self.environment.keys() - {"VAKINT_GALLERY_CACHE"}:
            with self.subTest(name=name), patch.dict(os.environ, {name: ""}):
                self.assertFalse(Gallery.benchmark_available())

    def test_launch_records_projection_bundle_and_background_worker(self):
        projection = "zero-momentum"
        inherited = []

        def spawn(*args, inherited=inherited, **kwargs):
            (descriptor,) = kwargs["pass_fds"]
            os.fstat(descriptor)
            inherited.append(descriptor)

        with patch("shared_gallery.subprocess.Popen", side_effect=spawn) as popen:
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
        self.assertNotIn("--cache", command)
        self.assertIn("--running-cache-capture-only", command)
        self.assertNotIn("--running-cache-no-read-delta", command)
        self.assertEqual(
            command[command.index("--running-cache") + 1],
            str(work / "running-cache"),
        )
        self.assertTrue(work.is_dir())
        self.assertFalse((work / "running-cache").exists())
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
        with self.assertRaises(OSError):
            os.fstat(inherited[0])

    def test_quadratic_launch_is_deferred_without_starting_a_worker(self):
        with (
            patch("shared_gallery.subprocess.Popen") as popen,
            self.assertRaisesRegex(ValueError, "quadratic reduction is deferred"),
        ):
            Gallery.start_benchmark("quadratic")
        popen.assert_not_called()
        self.assertFalse(self.work.exists())

    def test_cached_launch_uses_configured_cache(self):
        with patch("shared_gallery.subprocess.Popen"):
            work = Path(Gallery.start_benchmark("zero-momentum", cached=True))
        command = json.loads((work / "request.json").read_text())["command"]
        self.assertEqual(
            command[command.index("--cache") + 1],
            self.environment["VAKINT_GALLERY_CACHE"],
        )
        self.assertEqual(
            command[command.index("--running-cache") + 1],
            str(work / "running-cache"),
        )
        self.assertNotIn("--running-cache-capture-only", command)
        self.assertNotIn("--running-cache-no-read-delta", command)
        self.assertEqual(
            command[command.index("--progress-snapshot") + 1],
            str(work / "monitor.json"),
        )

    def test_cached_launch_accepts_an_immutable_native_directory(self):
        cache = Path(self.environment["VAKINT_GALLERY_CACHE"])
        cache.mkdir()
        index = cache / "bundle.rrindex"
        index.write_bytes(b"synthetic opaque native index")
        before = index.read_bytes()
        with patch("shared_gallery.subprocess.Popen"):
            work = Path(Gallery.start_benchmark("zero-momentum", cached=True))
        command = json.loads((work / "request.json").read_text())["command"]
        self.assertEqual(command[command.index("--cache") + 1], str(cache))
        self.assertEqual(index.read_bytes(), before)
        self.assertEqual(list(cache.iterdir()), [index])

    def test_missing_cache_does_not_launch_uncached_instead(self):
        os.environ.pop("VAKINT_GALLERY_CACHE")
        with (
            patch("shared_gallery.subprocess.Popen") as popen,
            self.assertRaisesRegex(ValueError, "cache"),
        ):
            Gallery.start_benchmark("zero-momentum", cached=True)
        popen.assert_not_called()

    def test_fresh_output_directories_share_only_prepared_input(self):
        with patch("shared_gallery.subprocess.Popen"):
            first = Path(Gallery.start_benchmark("zero-momentum"))
            second = Path(Gallery.start_benchmark("zero-momentum"))
        self.assertNotEqual(first, second)
        for work in (first, second):
            self.assertEqual(work.parent, self.work.resolve())
            self.assertTrue((work / "request.json").is_file())
            self.assertNotEqual(work, self.prepared / "zero-momentum")
            command = json.loads((work / "request.json").read_text())["command"]
            self.assertEqual(
                command[command.index("--running-cache") + 1],
                str(work / "running-cache"),
            )

    def test_held_lock_prevents_another_benchmark(self):
        self.work.mkdir()
        with (self.work / "benchmark.lock").open("a+") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with (
                patch("shared_gallery.subprocess.Popen") as popen,
                self.assertRaisesRegex(RuntimeError, "already running"),
            ):
                Gallery.start_benchmark("zero-momentum")
            popen.assert_not_called()
        self.assertEqual(list(self.work.iterdir()), [self.work / "benchmark.lock"])

    def test_launch_failure_closes_lock_without_retry(self):
        with (
            patch(
                "shared_gallery.subprocess.Popen",
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
        capture = self.work / "running-cache"
        capture.mkdir()
        prefix = capture / "committed-prefix"
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


if __name__ == "__main__":
    unittest.main()
