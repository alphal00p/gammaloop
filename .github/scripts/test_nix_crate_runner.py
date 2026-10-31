"""Exercise crate selection against real Nix archives (requires SYMBOLICA_LICENSE)."""

import importlib.util
import json
import os
import shutil
import signal
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import MagicMock, patch
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[2]
FIRST_TEST = "collinear_backtracking_cubic_keeps_distinct_parameter_crossings"
SECOND_TEST = "coil_pattern_adds_longitudinal_and_lateral_offsets"


class BrowserSelectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location(
            "browser_check", ROOT / "nix/browser-check.py"
        )
        cls.browser = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.browser)

    def report(self, directory, package, cases):
        # Real Nextest JUnit uses the binary ID for both testsuite name and
        # testcase classname, with the Rust test name in testcase name.
        root = ET.Element(
            "testsuites",
            name="nextest-run",
            tests=str(len(cases)),
            skipped="0",
            failures="0",
            errors="0",
        )
        for binary in dict.fromkeys(binary for binary, _ in cases):
            names = [name for candidate, name in cases if candidate == binary]
            suite = ET.SubElement(
                root,
                "testsuite",
                name=binary,
                tests=str(len(names)),
                skipped="0",
                failures="0",
                errors="0",
            )
            for name in names:
                ET.SubElement(suite, "testcase", name=name, classname=binary)
        destination = Path(directory) / f"{package}.xml"
        ET.ElementTree(root).write(destination, encoding="utf-8", xml_declaration=True)
        return destination

    def test_filtered_feynkit_clients(self):
        cff = (
            "feynkit-py::cff_display",
            "cff_explorer_keeps_native_families_and_graph_identity",
        )
        combined = ("feynkit-py::combined_community", "browser_display_fixtures")
        for cases, scripts in (
            ([cff], ["cff_display_browser.py"]),
            (
                [combined],
                [
                    "ltd_display_browser.py",
                    "energy_scoped_display_browser.py",
                    "generalized_energy_display_browser.py",
                ],
            ),
            (
                [cff, combined],
                [
                    "cff_display_browser.py",
                    "ltd_display_browser.py",
                    "energy_scoped_display_browser.py",
                    "generalized_energy_display_browser.py",
                ],
            ),
            ([], []),
        ):
            with self.subTest(cases=cases), tempfile.TemporaryDirectory() as directory:
                report = self.report(directory, "feynkit-py", cases)
                arguments = [
                    "browser-check.py",
                    "--package",
                    "feynkit-py",
                    "--report",
                    str(report),
                    "--browser-report",
                    str(Path(directory) / "browser"),
                    "--",
                    "cargo",
                    "nextest",
                    "run",
                ]
                environment = {
                    "CFF_DISPLAY_TEST_OUTPUT": "cff-browser",
                    "FEYNKIT_BROWSER_FIXTURES": "feynkit-browser",
                    "LINNET_TEST_CHROMIUM": "/pinned/chromium",
                }
                native = MagicMock()
                native.wait.return_value = 0
                native.poll.return_value = 0
                with (
                    patch.object(self.browser.sys, "argv", arguments),
                    patch.dict(os.environ, environment),
                    patch.object(self.browser.signal, "signal"),
                    patch.object(
                        self.browser.subprocess, "Popen", return_value=native
                    ) as launch,
                    patch.object(
                        self.browser.subprocess,
                        "run",
                    ) as run,
                ):
                    self.assertEqual(self.browser.main(), 0)
                self.assertEqual(launch.call_count, 1 + len(scripts))
                self.assertEqual(
                    launch.call_args_list[0].args[0], ["cargo", "nextest", "run"]
                )
                self.assertEqual(
                    [Path(call.args[0][1]).name for call in launch.call_args_list[1:]],
                    scripts,
                )
                run.assert_not_called()

    def test_full_tensor_and_mathml_identity(self):
        for package, cases, expected in (
            (
                "spynso3",
                [("spynso3::tensor_paging", "direct_tensor_paging")],
                ("tensor",),
            ),
            (
                "linnet-render-py",
                [("linnet-render-py::mathml_view", "mathml_view_browser_producer")],
                ("mathml",),
            ),
            (
                "linnet-render-py",
                [("linnet-render-py::mathml_view", "mathml_view_outlined_exports")],
                (),
            ),
            (
                "feynkit-py",
                [
                    (
                        "feynkit-py::combined_community",
                        "cff_explorer_keeps_native_families_and_graph_identity",
                    )
                ],
                (),
            ),
        ):
            with (
                self.subTest(package=package, cases=cases),
                tempfile.TemporaryDirectory() as directory,
            ):
                report = self.report(directory, package, cases)
                self.assertEqual(
                    self.browser.selected_browser_checks(package, report), expected
                )

    def test_skipped_producer_does_not_select_client(self):
        with tempfile.TemporaryDirectory() as directory:
            report = self.report(
                directory,
                "spynso3",
                [("spynso3::tensor_paging", "direct_tensor_paging")],
            )
            tree = ET.parse(report)
            tree.getroot().set("skipped", "1")
            tree.find("testsuite").set("skipped", "1")
            ET.SubElement(tree.find("testsuite/testcase"), "skipped")
            tree.write(report)
            self.assertEqual(
                self.browser.selected_browser_checks("spynso3", report), ()
            )

    def test_mathml_client_requires_started_and_executed_producer(self):
        producer = ("linnet-render-py::mathml_view", "mathml_view_browser_producer")
        for executed, started in ((False, False), (True, True), (True, False)):
            with (
                self.subTest(executed=executed, started=started),
                tempfile.TemporaryDirectory() as directory,
            ):
                report = self.report(
                    directory, "linnet-render-py", [producer] if executed else []
                )
                root = Path(directory) / "runtime"
                destination = Path(directory) / "browser"
                arguments = [
                    "browser-check.py",
                    "--package",
                    "linnet-render-py",
                    "--report",
                    str(report),
                    "--browser-report",
                    str(destination),
                    "--",
                    "cargo",
                    "nextest",
                    "run",
                ]
                native = MagicMock()
                native.poll.side_effect = [None, 0] if started else [0, 0]
                native.wait.return_value = 0
                client = MagicMock()
                client.poll.return_value = 0
                client.wait.return_value = 0

                def start_native(
                    command,
                    *_args,
                    started=started,
                    root=root,
                    native=native,
                    client=client,
                    **_kwargs,
                ):
                    if command[0] != "cargo":
                        return client
                    if started:
                        (root / "ci").mkdir()
                        (root / "ci/ready.json").write_text("{}", encoding="utf-8")
                    return native

                with (
                    patch.object(self.browser.sys, "argv", arguments),
                    patch.dict(os.environ, {"LINNET_MATHML_BROWSER_ROOT": str(root)}),
                    patch.object(self.browser.signal, "signal"),
                    patch("builtins.print"),
                    patch.object(
                        self.browser.subprocess, "Popen", side_effect=start_native
                    ) as launch,
                    patch.object(self.browser.subprocess, "run") as validators,
                ):
                    if executed != started:
                        with self.assertRaisesRegex(RuntimeError, "evidence disagree"):
                            self.browser.main()
                    else:
                        self.assertEqual(self.browser.main(), 0)
                self.assertEqual(launch.call_count, 1 + int(started))
                self.assertEqual(
                    launch.call_args_list[0].args[0], ["cargo", "nextest", "run"]
                )
                validators.assert_not_called()
                self.assertEqual(destination.exists(), started)

    def test_missing_and_malformed_evidence_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            missing = Path(directory) / "missing.xml"
            with self.assertRaises(FileNotFoundError):
                self.browser.selected_browser_checks("feynkit-py", missing)
            for content in (
                "<testsuites",
                '<testsuites name="nextest-run"/>',
                '<testsuites name="nextest-run" tests="1" skipped="0" failures="0" errors="0"/>',
                '<testsuite name="nextest-run" tests="0" skipped="0" failures="0" errors="0"/>',
            ):
                with self.subTest(content=content):
                    report = Path(directory) / "malformed.xml"
                    report.write_text(content, encoding="utf-8")
                    with self.assertRaises((ET.ParseError, KeyError, ValueError)):
                        self.browser.selected_browser_checks("feynkit-py", report)

    def test_no_tests_status_is_preserved_without_report_or_clients(self):
        arguments = [
            "browser-check.py",
            "--package",
            "feynkit-py",
            "--report",
            "/missing/report.xml",
            "--browser-report",
            "/unused/browser",
            "--",
            "cargo",
            "nextest",
            "run",
        ]
        native = MagicMock()
        native.wait.return_value = 4
        native.poll.return_value = 4
        with (
            patch.object(self.browser.sys, "argv", arguments),
            patch.object(self.browser.signal, "signal"),
            patch.object(
                self.browser.subprocess, "Popen", return_value=native
            ) as launch,
            patch.object(
                self.browser.subprocess,
                "run",
            ) as run,
        ):
            self.assertEqual(self.browser.main(), 4)
        launch.assert_called_once()
        native.wait.assert_called()
        run.assert_not_called()

    def test_stress_rejection_precedes_native_execution_and_side_effects(self):
        for package in ("feynkit-py", "spynso3", "linnet-render-py"):
            for flags in (
                ["--stress-count", "1"],
                ["--stress-count=2"],
                ["--stress-count=infinite"],
                ["--stress-duration", "1s"],
                ["--stress-duration=1s"],
            ):
                arguments = [
                    "browser-check.py",
                    "--package",
                    package,
                    "--report",
                    "/unused/report.xml",
                    "--browser-report",
                    "/unused/browser",
                    "--",
                    "cargo",
                    "nextest",
                    "run",
                    *flags,
                ]
                with (
                    self.subTest(package=package, flags=flags),
                    patch.object(self.browser.sys, "argv", arguments),
                    patch.object(self.browser.sys, "stderr") as stderr,
                    patch.object(self.browser, "run_owned") as native,
                    patch.object(self.browser.subprocess, "Popen") as launch,
                    patch.object(self.browser.subprocess, "run") as client,
                    patch.object(Path, "mkdir") as mkdir,
                ):
                    with self.assertRaises(SystemExit) as rejected:
                        self.browser.main()
                    self.assertEqual(rejected.exception.code, 2)
                    self.assertIn(
                        "does not support Nextest stress-count or stress-duration",
                        "".join(call.args[0] for call in stderr.write.call_args_list),
                    )
                native.assert_not_called()
                launch.assert_not_called()
                client.assert_not_called()
                mkdir.assert_not_called()

    def test_stress_tokens_after_native_separator_are_filters(self):
        command = [
            "cargo",
            "nextest",
            "run",
            "--",
            "--skip",
            "--stress-count",
            "--stress-duration=1s",
        ]
        arguments = [
            "browser-check.py",
            "--package",
            "feynkit-py",
            "--report",
            "/missing/report.xml",
            "--browser-report",
            "/unused/browser",
            "--",
            *command,
        ]
        with (
            patch.object(self.browser.sys, "argv", arguments),
            patch.dict(
                os.environ, {"NEXTEST_STRESS_CURRENT": "0", "NEXTEST_STRESS_TOTAL": "2"}
            ),
            patch.object(self.browser, "run_owned", return_value=(4, False)) as native,
        ):
            self.assertEqual(self.browser.main(), 4)
        native.assert_called_once_with(command, "feynkit-py", Path("/unused/browser"))

    def test_native_cancellation_forwards_repeats_and_reaps_test_groups(self):
        for package in ("feynkit-py", "spynso3", "linnet-render-py"):
            with (
                self.subTest(package=package),
                tempfile.TemporaryDirectory() as directory,
            ):
                handlers = {}
                state = {"running": True, "tests_running": True, "first_signal": False}
                native = MagicMock()

                def first_cancel(handlers=handlers, state=state):
                    if not state["first_signal"]:
                        state["first_signal"] = True
                        handlers[signal.SIGTERM](signal.SIGTERM, None)

                def poll(package=package, state=state, first_cancel=first_cancel):
                    if package == "linnet-render-py" and state["running"]:
                        first_cancel()
                    return None if state["running"] else 143

                def wait(
                    *args,
                    handlers=handlers,
                    state=state,
                    first_cancel=first_cancel,
                    **kwargs,
                ):
                    self.assertEqual((args, kwargs), ((), {}))
                    if state["running"]:
                        first_cancel()
                        handlers[signal.SIGINT](signal.SIGINT, None)
                        self.assertTrue(state["tests_running"])
                        # Only Nextest's completed wait certifies its test groups
                        # are gone; receiving a signal alone is insufficient.
                        state["tests_running"] = False
                        state["running"] = False
                    return 143

                native.poll.side_effect = poll
                native.wait.side_effect = wait
                previous = object()
                with (
                    patch.dict(
                        os.environ,
                        {
                            "LINNET_MATHML_BROWSER_ROOT": str(
                                Path(directory) / "runtime"
                            )
                        },
                    ),
                    patch.object(
                        self.browser.signal, "getsignal", return_value=previous
                    ),
                    patch.object(
                        self.browser.signal,
                        "signal",
                        side_effect=lambda sig, handler, handlers=handlers: (
                            handlers.__setitem__(sig, handler)
                        ),
                    ),
                    patch.object(
                        self.browser.subprocess, "Popen", return_value=native
                    ) as launch,
                    patch("builtins.print"),
                    self.assertRaises(SystemExit) as stopped,
                ):
                    self.browser.run_owned(
                        ["cargo", "nextest", "run"],
                        package,
                        Path(directory) / "browser",
                    )
                self.assertEqual(stopped.exception.code, 143)
                self.assertFalse(state["running"])
                self.assertFalse(state["tests_running"])
                self.assertEqual(
                    [call.args[0] for call in native.send_signal.call_args_list],
                    [signal.SIGTERM, signal.SIGINT],
                )
                native.kill.assert_not_called()
                native.terminate.assert_not_called()
                native.wait.assert_called()
                self.assertTrue(launch.call_args.kwargs["start_new_session"])
                self.assertEqual(
                    handlers, {signal.SIGTERM: previous, signal.SIGINT: previous}
                )

    def test_main_validator_cancellation_forwards_and_reaps_before_exit(self):
        for package, binary, producer, script in (
            (
                "feynkit-py",
                "cff_display",
                "cff_explorer_keeps_native_families_and_graph_identity",
                "cff_display_browser.py",
            ),
            (
                "feynkit-py",
                "combined_community",
                "browser_display_fixtures",
                "ltd_display_browser.py",
            ),
            (
                "spynso3",
                "tensor_paging",
                "direct_tensor_paging",
                "browser_math_font.py",
            ),
        ):
            with (
                self.subTest(package=package, binary=binary),
                tempfile.TemporaryDirectory() as directory,
            ):
                report = self.report(
                    directory, package, [(f"{package}::{binary}", producer)]
                )
                arguments = [
                    "browser-check.py",
                    "--package",
                    package,
                    "--report",
                    str(report),
                    "--browser-report",
                    str(Path(directory) / "browser"),
                    "--",
                    "cargo",
                    "nextest",
                    "run",
                ]
                native = MagicMock()
                native.wait.return_value = 0
                native.poll.return_value = 0
                validator = MagicMock()
                handlers = {}
                state = {"running": True}

                def wait(
                    *args, handlers=handlers, state=state, native=native, **kwargs
                ):
                    self.assertEqual((args, kwargs), ((), {}))
                    native.wait.assert_called()
                    if state["running"]:
                        handlers[signal.SIGTERM](signal.SIGTERM, None)
                        self.assertTrue(state["running"])
                        state["running"] = False
                    return -signal.SIGTERM

                validator.wait.side_effect = wait
                validator.poll.side_effect = lambda state=state: (
                    None if state["running"] else -signal.SIGTERM
                )
                previous = object()
                environment = {
                    "CFF_DISPLAY_TEST_OUTPUT": "cff-browser",
                    "FEYNKIT_BROWSER_FIXTURES": "feynkit-browser",
                    "SPENSO_MATHML_BROWSER_FIXTURES": "tensor-browser.json",
                    "LINNET_TEST_CHROMIUM": "/pinned/chromium",
                    "SYMBOLICA_DATA_DIRECTORY": "/pinned/render-data",
                }
                with (
                    patch.object(self.browser.sys, "argv", arguments),
                    patch.dict(os.environ, environment),
                    patch.object(
                        self.browser.signal, "getsignal", return_value=previous
                    ),
                    patch.object(
                        self.browser.signal,
                        "signal",
                        side_effect=lambda sig, handler, handlers=handlers: (
                            handlers.__setitem__(sig, handler)
                        ),
                    ),
                    patch.object(
                        self.browser.subprocess,
                        "Popen",
                        side_effect=[native, validator],
                    ) as launch,
                    patch.object(self.browser.subprocess, "run") as shortcut,
                    self.assertRaises(SystemExit) as stopped,
                ):
                    self.browser.main()
                self.assertEqual(stopped.exception.code, 143)
                self.assertFalse(state["running"])
                self.assertEqual(launch.call_count, 2)
                self.assertEqual(
                    launch.call_args_list[0].args[0], ["cargo", "nextest", "run"]
                )
                self.assertEqual(Path(launch.call_args_list[1].args[0][1]).name, script)
                validator.send_signal.assert_called_once_with(signal.SIGTERM)
                validator.wait.assert_called()
                validator.kill.assert_not_called()
                native.send_signal.assert_not_called()
                shortcut.assert_not_called()
                self.assertEqual(
                    handlers, {signal.SIGTERM: previous, signal.SIGINT: previous}
                )


def nix_environment(env=None):
    result = dict(os.environ if env is None else env)
    for key in (
        "SYMBOLICA_LICENSE",
        "SYMBOLICA_LICENSE_SIGNED",
        "GAMMALOOP_NIX_TEST_ARGS",
    ):
        result.pop(key, None)
    return result


@contextmanager
def flake_source(root):
    supplied = os.environ.get("GAMMALOOP_FLAKE_SOURCE")
    if supplied:
        source = Path(supplied).resolve()
        if not source.is_dir() or not (source / "flake.nix").is_file():
            raise ValueError(
                "GAMMALOOP_FLAKE_SOURCE must name a flake source directory"
            )
        yield source
        return
    # The test harness snapshots tracked live files explicitly; runtime commands
    # use the shell bootstrap and never import a Python source manager.
    command = (
        ["jj", "file", "list", "--no-pager", "--color", "never", "-T", 'path ++ "\\0"']
        if (root / ".jj").is_dir()
        else ["git", "ls-files", "--cached", "-z"]
    )
    paths = subprocess.check_output(command, cwd=root).split(b"\0")
    with tempfile.TemporaryDirectory(prefix="gammaloop-runner-test-") as temporary:
        source = Path(temporary)
        for name in dict.fromkeys(paths):
            if not name:
                continue
            relative = Path(os.fsdecode(name))
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError(f"Invalid tracked path: {relative}")
            if relative.parts[0] in ("target", ".venv"):
                continue
            if any((root / parent).is_symlink() for parent in relative.parents):
                raise ValueError(
                    f"Tracked path traverses a directory symlink: {relative}"
                )
            original = root / relative
            if not original.exists() and not original.is_symlink():
                continue
            destination = source / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original, destination, follow_symlinks=False)
        yield source


class CrateRunnerTests(unittest.TestCase):
    def command(self, *args, succeeds=True, env=None):
        if args[0] in ("nix", "nix-store"):
            env = nix_environment(env)
        result = subprocess.run(
            args, cwd=ROOT, env=env, text=True, capture_output=True, check=False
        )
        if succeeds:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stderr)
        return result

    def reports(self, result):
        paths = [
            Path(line.removeprefix("Test reports: "))
            for line in result.stderr.splitlines()
            if line.startswith("Test reports: ")
        ]
        self.assertEqual(len(paths), 1, result.stderr)
        output = paths[0]
        self.assertTrue(output.is_dir())
        self.assertEqual(
            output.resolve().parent, (ROOT / "target/test-reports").resolve()
        )
        return output

    def run_crate(self, *args):
        result = self.command("just", "test-nix", "kurvst", *args)
        output = self.reports(result)
        report = ET.parse(output / "linnet" / "kurvst.xml")
        names = {test.attrib["name"] for test in report.iter("testcase")}
        self.assertEqual(report.getroot().attrib["failures"], "0")
        self.assertEqual(report.getroot().attrib["errors"], "0")
        return output, names

    def derivation(self, attribute, *, env=None, succeeds=True):
        with flake_source(ROOT) as source:
            result = self.command(
                "nix",
                "eval",
                "--impure",
                "--accept-flake-config",
                "--no-write-lock-file",
                "--raw",
                f"path:{quote(str(source), safe='/')}#{attribute}.drvPath",
                env=env,
                succeeds=succeeds,
            )
        return result.stdout.strip() if succeeds else result

    def archive_inputs(self, derivation):
        """Follow runtime/archive edges without inspecting derivation environments."""
        pending = [derivation]
        visited = set()
        archives = set()
        while pending:
            path = pending.pop()
            if path in visited:
                continue
            visited.add(path)
            references = self.command("nix-store", "--query", "--references", path)
            for reference in references.stdout.splitlines():
                if not reference.endswith(".drv"):
                    continue
                name = Path(reference).name.split("-", 1)[1]
                if name.startswith("gammaloop-nextest-binaries-"):
                    archives.add(reference)
                elif name.startswith(("gammaloop-nextest-", "run-nextest-")):
                    pending.append(reference)
        return archives

    def test_filter_union_and_distinct_runtime_reports(self):
        output, names = self.run_crate(
            "-E", f"test({FIRST_TEST})", "--filterset", f"test({SECOND_TEST})"
        )
        self.assertEqual(
            names,
            {f"curve_api::tests::{FIRST_TEST}", f"curve_api::tests::{SECOND_TEST}"},
        )
        other_output, other_names = self.run_crate("--", FIRST_TEST)
        self.assertEqual(other_names, {f"curve_api::tests::{FIRST_TEST}"})
        self.assertNotEqual(output, other_output)
        self.assertTrue((output / "linnet/kurvst.xml").is_file())

    def test_selected_runner_reuses_ci_archive(self):
        # Reports belong to runtime invocations, so inspect the runner's pure
        # derivation graph separately rather than asking for a report deriver.
        env = nix_environment()
        env.update(
            GAMMALOOP_NIX_TEST_PACKAGES=json.dumps(["kurvst"]),
            GAMMALOOP_NIX_TEST_PROFILE="ci-optim",
            GAMMALOOP_NIX_TEST_LIST="0",
        )
        runner = self.derivation("test-selected", env=env)
        archives = self.archive_inputs(runner)
        self.assertEqual(len(archives), 1)
        self.assertIn(
            "-gammaloop-nextest-binaries-linnet-kurvst-", next(iter(archives))
        )

        system = self.command(
            "nix", "eval", "--impure", "--raw", "--expr", "builtins.currentSystem"
        ).stdout.strip()
        ci_archive = self.derivation(
            f"checks.{system}.gammaloop-nextest-binaries-linnet", env=env
        )
        self.assertTrue(archives <= self.archive_inputs(ci_archive))
        env["GAMMALOOP_NIX_TEST_LIST"] = "1"
        list_runner = self.derivation("test-selected", env=env)
        self.assertNotEqual(runner, list_runner)
        self.assertEqual(archives, self.archive_inputs(list_runner))

    def test_feature_validation_is_evaluation_only(self):
        env = nix_environment()
        env.update(
            GAMMALOOP_NIX_TEST_PACKAGES='["kurvst"]',
            GAMMALOOP_NIX_TEST_PROFILE="ci-optim",
            GAMMALOOP_NIX_TEST_LIST="0",
        )
        original = self.derivation("test-selected", env=env)
        for features in (
            '["missing-feature"]',
            '["kurvst/missing-feature"]',
            "{}",
            "[1]",
        ):
            with self.subTest(features=features):
                env["GAMMALOOP_NIX_TEST_FEATURES"] = features
                failure = self.derivation("test-selected", env=env, succeeds=False)
                self.assertIn("feature", failure.stderr.lower())
        env["GAMMALOOP_NIX_TEST_FEATURES"] = "[]"
        self.assertEqual(original, self.derivation("test-selected", env=env))

    def test_supported_feature_groups_preserve_archive_selection(self):
        env = nix_environment()
        env.update(
            GAMMALOOP_NIX_TEST_PACKAGES='["spenso"]',
            GAMMALOOP_NIX_TEST_PROFILE="ci-optim",
            GAMMALOOP_NIX_TEST_LIST="0",
            GAMMALOOP_NIX_TEST_FEATURES="[]",
        )
        original = self.derivation("test-selected", env=env)
        for features in (
            ["shadowing"],
            ["spenso/shadowing"],
            ["shadowing,python_stubgen"],
        ):
            with self.subTest(features=features):
                env["GAMMALOOP_NIX_TEST_FEATURES"] = json.dumps(features)
                self.assertEqual(original, self.derivation("test-selected", env=env))
        env["GAMMALOOP_NIX_TEST_PACKAGES"] = '["kurvst", "spenso"]'
        env["GAMMALOOP_NIX_TEST_FEATURES"] = '["spenso/shadowing"]'
        self.derivation("test-selected", env=env)
        env["GAMMALOOP_NIX_TEST_FEATURES"] = '["shadowing"]'
        self.derivation("test-selected", env=env, succeeds=False)

    def test_target_scope_intersects_user_filter_union(self):
        _, names = self.run_crate(
            "--lib", "-E", f"test({FIRST_TEST})", "-E", f"test({SECOND_TEST})"
        )
        self.assertEqual(
            names,
            {f"curve_api::tests::{FIRST_TEST}", f"curve_api::tests::{SECOND_TEST}"},
        )
        # test_runs is declared in this workspace but absent from kurvst's
        # archive. A genuinely unknown name is a filter error, not an empty set.
        for selector in (("--test", "test_runs"), ("--test=test_runs",)):
            with self.subTest(selector=selector):
                _, names = self.run_crate(*selector, "--no-tests=pass", "-E", "all()")
                self.assertEqual(names, set())
        _, names = self.run_crate("--test=test_runs", "--lib", "--", FIRST_TEST)
        self.assertEqual(names, {f"curve_api::tests::{FIRST_TEST}"})

    def test_unknown_target_name_remains_a_filter_error(self):
        for selector in (("--test", "missing-binary"), ("--test=missing-binary",)):
            with self.subTest(selector=selector):
                result = self.command(
                    "just",
                    "test-nix",
                    "kurvst",
                    *selector,
                    "--no-tests=pass",
                    succeeds=False,
                )
                self.assertEqual(result.returncode, 94, result.stderr)
                self.assertIn("operator didn't match any binary names", result.stderr)

    def test_list_inventory_is_fresh_and_contextual(self):
        for selector in (FIRST_TEST, SECOND_TEST):
            result = self.command(
                "bash",
                "scripts/develop.sh",
                "--backend",
                "nix",
                "list",
                "-p",
                "kurvst",
                "--lib",
                "--",
                selector,
            )
            inventory = json.loads(result.stdout)
            reports = self.reports(result)
            self.assertEqual(set(inventory), {"linnet"})
            self.assertEqual(set(inventory["linnet"]), {"kurvst"})
            self.assertEqual(
                inventory["linnet"]["kurvst"],
                json.loads((reports / "linnet/kurvst.json").read_text()),
            )
            suites = inventory["linnet"]["kurvst"]["rust-suites"]
            matched = {
                name
                for suite in suites.values()
                for name, case in suite["testcases"].items()
                if case["filter-match"]["status"] == "matches"
            }
            self.assertEqual(matched, {f"curve_api::tests::{selector}"})

    def test_list_inventory_merges_runtime_contexts(self):
        result = self.command(
            "bash",
            "scripts/develop.sh",
            "--backend",
            "nix",
            "list",
            "-p",
            "gammaloop-integration-tests",
            "-E",
            "none()",
        )
        inventory = json.loads(result.stdout)
        reports = self.reports(result)
        self.assertEqual(set(inventory), {"integration", "python-api"})
        for context in inventory:
            self.assertEqual(set(inventory[context]), {"gammaloop-integration-tests"})
            self.assertEqual(
                inventory[context]["gammaloop-integration-tests"],
                json.loads(
                    (reports / context / "gammaloop-integration-tests.json").read_text()
                ),
            )

    def test_empty_selection_fails(self):
        result = self.command(
            "bash",
            "scripts/develop.sh",
            "--backend",
            "nix",
            "test",
            "-p",
            "kurvst",
            "-E",
            "none()",
            succeeds=False,
        )
        self.assertEqual(result.returncode, 4, result.stderr)
        self.assertIn("No tests selected in crate selected", result.stderr)
        self.assertTrue((self.reports(result) / "linnet/no-tests").is_file())

    def test_empty_selection_can_be_allowed(self):
        for policy in ("pass", "warn"):
            with self.subTest(policy=policy):
                _, names = self.run_crate("-E", "none()", f"--no-tests={policy}")
                self.assertEqual(names, set())

    def test_profile_override_and_substring(self):
        _, names = self.run_crate(
            "-P", "ci_slow", "--ignore-default-filter", "--", FIRST_TEST
        )
        self.assertEqual(names, {f"curve_api::tests::{FIRST_TEST}"})

    def test_invalid_filter_fails(self):
        result = self.command(
            "bash",
            "scripts/develop.sh",
            "--backend",
            "nix",
            "test",
            "-p",
            "kurvst",
            "-E",
            "(",
            succeeds=False,
        )
        self.assertNotEqual(result.returncode, 4, result.stderr)
        self.reports(result)

    def test_short_filter_and_runtime_options(self):
        _, names = self.run_crate(
            f"-E=test({FIRST_TEST})", "--max-fail", "2", "--final-status-level", "pass"
        )
        self.assertEqual(names, {f"curve_api::tests::{FIRST_TEST}"})

    def test_filter_clusters_are_rejected(self):
        result = self.command(
            "just", "test-nix", "kurvst", "-vE", "none()", succeeds=False
        )
        self.assertIn("Pass -E and -P separately", result.stderr)

    def test_unknown_crate_fails(self):
        result = self.command(
            "just", "test-nix", "missing-workspace-crate", succeeds=False
        )
        self.assertIn("Unknown test selection. Available crates:", result.stderr)

    def test_split_runtime_contexts(self):
        result = self.command(
            "just",
            "test-nix",
            "gammaloop-integration-tests",
            "--test",
            "test_python_api",
            "-E",
            "test(python_settings_expose_ecm_relative_stability_tolerances)",
        )
        output = self.reports(result)
        self.assertTrue((output / "integration" / "no-tests").is_file())
        report = ET.parse(output / "python-api" / "gammaloop-integration-tests.xml")
        self.assertEqual(report.getroot().attrib["tests"], "1")
        self.assertEqual(report.getroot().attrib["failures"], "0")

    def test_docs_runtime_fixtures(self):
        result = self.command(
            "just",
            "test-nix",
            "alphal00p-docs-builder",
            "-E",
            "test(checked_in_registry_is_valid)",
        )
        output = self.reports(result)
        report = ET.parse(output / "docs" / "alphal00p-docs-builder.xml")
        self.assertEqual(report.getroot().attrib["tests"], "1")
        self.assertEqual(report.getroot().attrib["failures"], "0")

    def test_insta_snapshot_fixture(self):
        result = self.command(
            "just",
            "test-nix",
            "idenso",
            "-E",
            "test(=reference_cases::tests::reference_case_simplification_snapshots)",
        )
        report = ET.parse(self.reports(result) / "spenso/idenso.xml")
        self.assertEqual(report.getroot().attrib["tests"], "1")
        self.assertEqual(report.getroot().attrib["failures"], "0")


if __name__ == "__main__":
    unittest.main()
