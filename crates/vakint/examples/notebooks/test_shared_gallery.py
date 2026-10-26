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
        for projection in ("quadratic", "zero-momentum"):
            with self.subTest(projection=projection):
                inherited = []

                def spawn(*args, inherited=inherited, **kwargs):
                    (descriptor,) = kwargs["pass_fds"]
                    os.fstat(descriptor)
                    inherited.append(descriptor)

                with patch(
                    "shared_gallery.subprocess.Popen", side_effect=spawn
                ) as popen:
                    work = Path(Gallery.start_benchmark(projection))
                command = json.loads((work / "request.json").read_text())["command"]
                self.assertEqual(
                    command[0], self.environment["VAKINT_GALLERY_BENCHMARK"]
                )
                self.assertEqual(command[command.index("--projection") + 1], projection)
                self.assertEqual(
                    command[command.index("--prepared-input") + 1],
                    str((self.prepared / projection).resolve()),
                )
                self.assertEqual(command[command.index("--work-dir") + 1], str(work))
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
                self.assertEqual(
                    Path(arguments[1]), Gallery.directory / "benchmark_job.py"
                )
                self.assertEqual(arguments[2:], [str(work), str(inherited[0])])
                self.assertEqual(popen.call_args.kwargs["stderr"], subprocess.STDOUT)
                self.assertNotIn("timeout", popen.call_args.kwargs)
                with self.assertRaises(OSError):
                    os.fstat(inherited[0])

    def test_cached_launch_uses_configured_cache(self):
        with patch("shared_gallery.subprocess.Popen"):
            work = Path(Gallery.start_benchmark("quadratic", cached=True))
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

    def test_cached_launch_accepts_an_immutable_native_directory(self):
        cache = Path(self.environment["VAKINT_GALLERY_CACHE"])
        cache.mkdir()
        index = cache / "bundle.rrindex"
        index.write_bytes(b"synthetic opaque native index")
        before = index.read_bytes()
        with patch("shared_gallery.subprocess.Popen"):
            work = Path(Gallery.start_benchmark("quadratic", cached=True))
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
            Gallery.start_benchmark("quadratic", cached=True)
        popen.assert_not_called()

    def test_fresh_output_directories_share_only_prepared_input(self):
        with patch("shared_gallery.subprocess.Popen"):
            first = Path(Gallery.start_benchmark("quadratic"))
            second = Path(Gallery.start_benchmark("quadratic"))
        self.assertNotEqual(first, second)
        for work in (first, second):
            self.assertEqual(work.parent, self.work.resolve())
            self.assertTrue((work / "request.json").is_file())
            self.assertNotEqual(work, self.prepared / "quadratic")
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
                Gallery.start_benchmark("quadratic")
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
            Gallery.start_benchmark("quadratic")
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

    def run_worker(self, *, result=None, error=None, output=""):
        descriptor = os.open(
            self.work / "benchmark.lock", os.O_CREAT | os.O_RDWR, 0o600
        )

        def run(command, **kwargs):
            self.assertEqual(command, self.command)
            self.assertNotIn("timeout", kwargs)
            self.assertFalse(kwargs["check"])
            self.assertEqual(kwargs["stderr"], subprocess.STDOUT)
            os.fstat(descriptor)
            kwargs["stdout"].write(output)
            if error is not None:
                raise error
            return SimpleNamespace(returncode=result)

        with (
            patch.object(
                sys, "argv", ["benchmark_job.py", str(self.work), str(descriptor)]
            ),
            patch("benchmark_job.subprocess.run", side_effect=run) as subprocess_run,
        ):
            benchmark_job.main()
        subprocess_run.assert_called_once()
        with self.assertRaises(OSError):
            os.fstat(descriptor)
        return json.loads((self.work / "status.json").read_text())

    def test_success_has_no_default_deadline_and_preserves_receipt(self):
        receipt = {"fresh_native_preparation": False, "preparation_seconds": None}
        (self.work / "benchmark.json").write_text(json.dumps(receipt))
        status = self.run_worker(result=0)
        self.assertEqual(status["state"], "complete")
        self.assertEqual(status["receipt"], receipt)
        self.assertGreaterEqual(status["elapsed_seconds"], 0)

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

    def test_spawn_error_is_reported_and_releases_descriptor(self):
        status = self.run_worker(error=OSError("synthetic worker error"))
        self.assertEqual(status["state"], "failed")
        self.assertIn("synthetic worker error", status["message"])


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
            Gallery=SimpleNamespace(metrics=lambda *_: "metrics"),
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
