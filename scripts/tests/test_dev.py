"""Focused frontend tests; Nix execution is always mocked."""

import fcntl
import importlib.util
import io
import json
import os
import shutil
import stat
import subprocess
import tempfile
import threading
import unittest
from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    "gammaloop_dev", Path(__file__).resolve().parents[1] / "dev.py"
)
dev = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dev)

CATALOG = {
    "packages": ["feynkit-core", "gammalooprs", "gammaloop-api", "gammaloop-tests"],
    "groups": [
        {
            "name": "core",
            "packages": ["feynkit-core", "gammalooprs", "gammaloop-api"],
            "default": True,
        },
        {"name": "integration", "packages": ["gammaloop-tests"], "default": True},
        {"name": "python-api", "packages": ["gammaloop-api"], "default": False},
    ],
    "cargoProfile": "ci-optim",
}
SOURCE = "/nix/store/00000000000000000000000000000000-source"


class TemporaryTestCase(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="gammaloop-dev-test-")
        self.addCleanup(temporary.cleanup)
        self.temporary = Path(temporary.name)

    def make_result(
        self, name, *, contexts=None, cases=1, status=0, packages=(), logs=False
    ):
        wrapper = self.temporary / name
        result = wrapper / "result"
        result.mkdir(parents=True)
        (wrapper / "development-state.json").write_text(
            json.dumps(
                {"version": 1, "contexts": contexts or {"core": {"opaque": name}}}
            )
        )
        (result / "exit-status").write_text(f"{status}\n")
        (result / "tests.xml").write_text(
            "<testsuite>" + '<testcase name="example"/>' * cases + "</testsuite>"
        )
        for package in packages:
            (result / package).write_text("usable product")
        if logs:
            (result / "logs").mkdir()
            (result / "logs" / "trace.jsonl").write_text('{"event":"fixture"}\n')
        return wrapper


class SelectionTests(unittest.TestCase):
    def requests(self, *arguments, environment=None):
        return dev.requests_for(
            dev.parse_args(list(arguments)),
            CATALOG,
            {} if environment is None else environment,
        )

    def test_default_and_name_only_selections_keep_all_default_groups(self):
        for arguments in (
            ("test",),
            ("test", "a test name"),
            ("test", "-E", "test(alpha)"),
        ):
            with self.subTest(arguments=arguments):
                requests = self.requests(*arguments)
                self.assertEqual(
                    [request["group"] for request in requests], ["core", "integration"]
                )
                self.assertEqual(
                    [request["packages"] for request in requests],
                    [
                        group["packages"]
                        for group in CATALOG["groups"]
                        if group["default"]
                    ],
                )
                self.assertTrue(
                    all(request["noTests"] == "pass" for request in requests)
                )

    def test_packages_split_across_groups_without_duplicates(self):
        requests = self.requests(
            "test", "-p", "gammaloop-tests", "-p", "gammalooprs", "-p", "gammalooprs"
        )
        self.assertEqual(
            [(request["group"], request["packages"]) for request in requests],
            [("core", ["gammalooprs"]), ("integration", ["gammaloop-tests"])],
        )
        self.assertEqual(
            self.requests("test", "--group", "python-api")[0]["group"], "python-api"
        )
        for arguments in (
            ("test", "--group", "missing"),
            ("test", "-p", "missing"),
            ("test", "--group", "python-api", "-p", "gammalooprs"),
        ):
            with self.subTest(arguments=arguments), self.assertRaises(dev.DevError):
                self.requests(*arguments)

    def test_filter_union_is_intersected_with_binary_union_and_lib(self):
        request = self.requests(
            "test",
            "-E",
            "test(alpha)",
            "-E",
            "test(beta)",
            "--test",
            "test_one",
            "--test",
            "test_two",
            "--lib",
            "--profile",
            "ci_full",
        )[0]
        self.assertEqual(
            request["filter"],
            "((test(alpha)) or (test(beta))) and "
            "((kind(test) and binary(=test_one)) or "
            "(kind(test) and binary(=test_two)) or kind(lib))",
        )
        self.assertEqual(request["profile"], "ci_full")
        self.assertEqual(self.requests("test", "--lib")[0]["filter"], "(kind(lib))")

    def test_static_requests_have_no_test_only_inputs(self):
        for command in ("build", "check", "clippy"):
            with self.subTest(command=command):
                self.assertEqual(
                    self.requests(command, environment={"INSTA_UPDATE": "always"}),
                    [{"command": command, "packages": CATALOG["packages"]}],
                )
                self.assertEqual(
                    self.requests(command, "-p", "feynkit-core"),
                    [{"command": command, "packages": ["feynkit-core"]}],
                )

    def test_python_build_requires_only_api_package(self):
        expected = [
            {"command": "build", "packages": ["gammaloop-api"], "variant": "python"}
        ]
        self.assertEqual(self.requests("build", "--python"), expected)
        self.assertEqual(
            self.requests("build", "--python", "-p", "gammaloop-api"), expected
        )
        with self.assertRaises(dev.DevError):
            self.requests("build", "--python", "-p", "gammalooprs")

    def test_compilation_flags_are_rejected_on_both_sides_of_separator(self):
        flags = (
            ["--release"],
            ["--features", "python_api"],
            ["--all-features"],
            ["--no-default-features"],
            ["--target", "wasm32-unknown-unknown"],
            ["--config", "other.toml"],
            ["--profile", "release"],
        )
        for flag in flags:
            for separator in ([], ["--"]):
                with (
                    self.subTest(flag=flag, separator=separator),
                    redirect_stderr(io.StringIO()),
                    self.assertRaises((SystemExit, dev.DevError)),
                ):
                    dev.parse_args(["test", *separator, *flag])
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            dev.parse_args(["build", "--", "--nocapture"])

    def test_runtime_argv_remains_literal_and_owned_options_are_rejected(self):
        literals = [
            "a name",
            "'quoted'",
            '"double"',
            r"a\b",
            "$(touch marker)",
            "`touch marker`",
            "line\nbreak",
        ]
        args = dev.parse_args(
            ["test", "--retries", "2", "--", "--skip", literals[0], *literals[1:]]
        )
        self.assertEqual(
            args.nextest_args, ["--retries", "2", "--", "--skip", *literals]
        )
        dev.validate_runtime(args.nextest_args)
        for arguments in (
            ["--", "-p", "gammalooprs"],
            ["--no-tests=pass"],
            ["--retries"],
            ["--fail-fast=yes"],
        ):
            with self.subTest(arguments=arguments), self.assertRaises(dev.DevError):
                dev.validate_runtime(arguments)

    def test_runtime_environment_is_explicit_and_license_free(self):
        environment = {
            "GL_DISPLAY_FILTER": "off",
            "GL_LOGFILE_FILTER": "[{#uv}]=debug",
            "GL_ALL_LOG_FILTER": "warn",
            "RUST_MIN_STACK": "8388608",
            "RAYON_NUM_THREADS": "3",
            "INSTA_UPDATE": "no",
            "GL_TEST_LOG_DIR": "/caller/logs",
            "SYMBOLICA_LICENSE": "fixture-secret",
            "RUSTFLAGS": "-C target-cpu=native",
            "UNRELATED": "ignored",
        }
        request = self.requests(
            "test", "--memory-limit-gb", "2.5", environment=environment
        )[0]
        self.assertEqual(
            request["runtimeEnv"],
            {
                name: environment[name]
                for name in (
                    "GL_DISPLAY_FILTER",
                    "GL_LOGFILE_FILTER",
                    "GL_ALL_LOG_FILTER",
                    "RUST_MIN_STACK",
                    "RAYON_NUM_THREADS",
                    "INSTA_UPDATE",
                )
            },
        )
        self.assertTrue(request["collectLogs"])
        self.assertEqual(request["memoryLimitGb"], 2.5)
        self.assertNotIn("fixture-secret", json.dumps(request))
        self.assertTrue(self.requests("test", "--collect-logs")[0]["collectLogs"])

    def test_invalid_runtime_environment_and_memory_limit_are_rejected(self):
        for environment in (
            {"GL_DISPLAY_FILTER": ""},
            {"GL_LOGFILE_FILTER": "bad\0value"},
            {"RUST_MIN_STACK": "0"},
            {"RAYON_NUM_THREADS": "1.5"},
            {"INSTA_UPDATE": "always"},
        ):
            with self.subTest(environment=environment), self.assertRaises(dev.DevError):
                self.requests("test", environment=environment)
        for value in ("0", "-1", "nan", "inf"):
            with (
                self.subTest(value=value),
                redirect_stderr(io.StringIO()),
                self.assertRaises(SystemExit),
            ):
                dev.parse_args(["test", "--memory-limit-gb", value])

    def test_inherited_threads_yield_to_explicit_nextest_or_libtest_threads(self):
        environment = {"NEXTEST_TEST_THREADS": "7"}
        self.assertEqual(
            self.requests("test", environment=environment)[0]["nextestArgs"],
            ["--test-threads", "7"],
        )
        for arguments in (
            ["--test-threads", "2"],
            ["--test-threads=2"],
            ["-j", "2"],
            ["-j2"],
            ["--", "--test-threads", "2"],
        ):
            with self.subTest(arguments=arguments):
                expected = dev.parse_args(["test", *arguments]).nextest_args
                self.assertEqual(
                    self.requests("test", *arguments, environment=environment)[0][
                        "nextestArgs"
                    ],
                    expected,
                )
        with self.assertRaises(dev.DevError):
            self.requests("test", environment={"NEXTEST_TEST_THREADS": "0"})

    def test_rerun_and_rebuild_add_one_fresh_test_nonce(self):
        first = self.requests("test", "--rerun")
        second = self.requests("test", "--rebuild")
        self.assertRegex(first[0]["rerun"], r"^[0-9a-f]{32}$")
        self.assertEqual({request["rerun"] for request in first}, {first[0]["rerun"]})
        self.assertNotEqual(first[0]["rerun"], second[0]["rerun"])
        self.assertNotIn("rerun", self.requests("test")[0])


class SnapshotTests(TemporaryTestCase):
    def setUp(self):
        super().setUp()
        self.root = self.temporary / "worktree"
        self.root.mkdir()
        self.git("init", "--quiet")
        (self.root / ".gitignore").write_text("target/\n.venv/\n.jj/\n")
        (self.root / "dirty.rs").write_text("staged old source\n")
        (self.root / "deleted.rs").write_text("deleted source\n")
        (self.root / "executable").write_text("#!/bin/sh\n")
        (self.root / "executable").chmod(0o644)
        (self.root / "link").symlink_to("deleted.rs")
        self.git("add", "--", ".")
        (self.root / "executable").chmod(0o755)
        (self.root / "dirty.rs").write_text("current dirty source\n")
        (self.root / "deleted.rs").unlink()
        (self.root / "new\nsource.rs").write_text("new source\n")
        (self.root / "link").unlink()
        (self.root / "link").symlink_to("dirty.rs")
        for directory in ("target", ".venv"):
            (self.root / directory).mkdir()
            (self.root / directory / "large-ignored").write_bytes(b"ignored" * 300_000)
        self.snapshot = dev.SourceSnapshot(self.root)

    def git(self, *arguments):
        return subprocess.check_output(
            ["git", "--no-optional-locks", *arguments],
            cwd=self.root,
            env=os.environ
            | {"GIT_CONFIG_NOSYSTEM": "1", "GIT_CONFIG_GLOBAL": os.devnull},
            stderr=subprocess.PIPE,
        )

    def assert_current_snapshot(self, snapshot):
        index = self.root / ".git" / "index"
        before = (index.read_bytes(), index.stat().st_mtime_ns)
        destination = self.temporary / "captured"
        destination.mkdir()
        snapshot.copy(destination)
        self.assertEqual(
            (destination / "dirty.rs").read_text(), "current dirty source\n"
        )
        self.assertEqual((destination / "new\nsource.rs").read_text(), "new source\n")
        self.assertFalse((destination / "deleted.rs").exists())
        self.assertEqual(os.readlink(destination / "link"), "dirty.rs")
        self.assertEqual(
            stat.S_IMODE((destination / "executable").stat().st_mode), 0o755
        )
        self.assertEqual(stat.S_IMODE((destination / "dirty.rs").stat().st_mode), 0o644)
        self.assertEqual(
            {path.name for path in destination.iterdir()},
            {".gitignore", "dirty.rs", "executable", "link", "new\nsource.rs"},
        )
        self.assertEqual((index.read_bytes(), index.stat().st_mtime_ns), before)
        self.assertEqual(self.git("show", ":dirty.rs"), b"staged old source\n")

    def test_git_capture_uses_dirty_files_without_touching_index(self):
        self.assert_current_snapshot(self.snapshot)

    def test_jj_listing_uses_nonmutating_flags_and_an_empty_git_index(self):
        tracked = self.git("ls-files", "--cached", "-z")
        (self.root / ".jj").mkdir()
        (self.root / ".jj" / "large-ignored").write_bytes(b"ignored" * 300_000)
        check_output = subprocess.check_output
        indexes = []
        jj_commands = []

        def output(command, **kwargs):
            if command[0] == "jj":
                jj_commands.append(command)
                self.assertIn("--ignore-working-copy", command)
                self.assertIn("--at-op=@", command)
                if "list" in command:
                    return tracked
                self.assertEqual(command[-2:], ["git", "root"])
                return os.fsencode(self.root / ".git") + b"\n"
            self.assertEqual(command[0], "git")
            if "ls-files" not in command:
                return check_output(command, **kwargs)
            index = Path(kwargs["env"]["GIT_INDEX_FILE"])
            self.assertFalse(index.exists())
            self.assertTrue(index.parent.is_dir())
            indexes.append(index)
            return check_output(command, **kwargs)

        with patch.object(dev.subprocess, "check_output", side_effect=output):
            self.assert_current_snapshot(self.snapshot)
        self.assertEqual(len(jj_commands), 4)
        self.assertEqual(len(indexes), 2)
        self.assertTrue(all(not index.parent.exists() for index in indexes))

    def test_directory_replaced_by_symlink_captures_link_without_reading_target(self):
        directory = self.root / "nested"
        directory.mkdir()
        (directory / "secret.rs").write_text("tracked source")
        self.git("add", "--", "nested/secret.rs")
        shutil.rmtree(directory)
        outside = self.temporary / "outside"
        outside.mkdir()
        (outside / "secret.rs").write_text("must not read")
        directory.symlink_to(outside, target_is_directory=True)
        destination = self.temporary / "captured"
        destination.mkdir()
        copied_inodes = []
        copy = shutil.copyfileobj

        def copy_source(source, target):
            copied_inodes.append(os.fstat(source.fileno()).st_ino)
            copy(source, target)

        with patch.object(dev.shutil, "copyfileobj", side_effect=copy_source):
            self.snapshot.copy(destination)
        self.assertTrue((destination / "nested").is_symlink())
        self.assertEqual(os.readlink(destination / "nested"), str(outside))
        self.assertNotIn((outside / "secret.rs").stat().st_ino, copied_inodes)
        self.assertEqual((outside / "secret.rs").read_text(), "must not read")

    def test_directory_replaced_by_file_omits_old_tracked_descendants(self):
        directory = self.root / "nested"
        directory.mkdir()
        (directory / "child.rs").write_text("old child")
        self.git("add", "--", "nested/child.rs")
        shutil.rmtree(directory)
        directory.write_text("new parent file")
        destination = self.temporary / "captured"
        destination.mkdir()
        self.snapshot.copy(destination)
        self.assertTrue((destination / "nested").is_file())
        self.assertEqual((destination / "nested").read_text(), "new parent file")

    def test_snapshot_revalidates_files_after_later_files_are_copied(self):
        copy = shutil.copyfileobj
        changed = False

        def changing_copy(source, target):
            nonlocal changed
            copy(source, target)
            if Path(target.name).name == "executable":
                (self.root / "dirty.rs").write_text(
                    "changed after its copy completed\n"
                )
                changed = True

        destination = self.temporary / "captured"
        destination.mkdir()
        with (
            patch.object(dev.shutil, "copyfileobj", side_effect=changing_copy),
            self.assertRaisesRegex(dev.DevError, "source changed"),
        ):
            self.snapshot.copy(destination)
        self.assertTrue(changed)

    def test_snapshot_rejects_change_during_file_copy(self):
        copy = shutil.copyfileobj

        def changing_copy(source, target):
            copy(source, target)
            if Path(target.name).name == "dirty.rs":
                (self.root / "dirty.rs").write_text("source changed while copying\n")

        destination = self.temporary / "captured"
        destination.mkdir()
        with (
            patch.object(dev.shutil, "copyfileobj", side_effect=changing_copy),
            self.assertRaisesRegex(dev.DevError, "source changed"),
        ):
            self.snapshot.copy(destination)

    def test_tracked_file_replaced_by_directory_captures_new_child(self):
        (self.root / "dirty.rs").unlink()
        (self.root / "dirty.rs").mkdir()
        (self.root / "dirty.rs" / "new.rs").write_text("replacement source\n")
        destination = self.temporary / "captured"
        destination.mkdir()
        self.snapshot.copy(destination)
        self.assertEqual(
            (destination / "dirty.rs" / "new.rs").read_text(), "replacement source\n"
        )

    def test_equivalent_worktrees_use_identical_store_input_and_fixed_basename(self):
        other = self.temporary / "differently-named-worktree"
        shutil.copytree(self.root, other, symlinks=True)
        check_output = subprocess.check_output
        captures = []
        temporary_sources = []

        def output(command, **kwargs):
            if command[0] != "nix":
                return check_output(command, **kwargs)
            self.assertEqual(
                command[:5], ["nix", "store", "add-path", "--name", "source"]
            )
            source = Path(command[5])
            self.assertEqual(source.name, "source")
            temporary_sources.append(source)
            captures.append(
                {
                    str(path.relative_to(source)): (
                        ("symlink", os.readlink(path))
                        if path.is_symlink()
                        else (
                            "file",
                            stat.S_IMODE(path.stat().st_mode),
                            path.read_bytes(),
                        )
                    )
                    for path in source.rglob("*")
                    if path.is_file() or path.is_symlink()
                }
            )
            return SOURCE + "\n"

        with patch.object(dev.subprocess, "check_output", side_effect=output):
            self.assertEqual(self.snapshot.store(), SOURCE)
            self.assertEqual(dev.SourceSnapshot(other).store(), SOURCE)
        self.assertEqual(captures[0], captures[1])
        self.assertTrue(all(not source.exists() for source in temporary_sources))


class TransportTests(TemporaryTestCase):
    def test_catalog_and_default_system_use_pinned_development_namespace(self):
        with patch.object(
            dev.subprocess,
            "check_output",
            side_effect=["x86_64-linux\n", json.dumps(CATALOG)],
        ) as output:
            development = dev.NixDevelopment(self.temporary, SOURCE)
            self.assertEqual(development.catalog(), CATALOG)
        self.assertEqual(
            output.call_args_list[0].args[0][-2:], ["--expr", "builtins.currentSystem"]
        )
        command = output.call_args_list[1].args[0]
        self.assertIn("--no-update-lock-file", command)
        self.assertEqual(command[-1], SOURCE + "#lib.x86_64-linux.development.catalog")
        with self.assertRaises(dev.DevError):
            dev.NixDevelopment(self.temporary, SOURCE, 'linux"; abort "injected')

    def test_json_file_preserves_literal_arguments_and_large_previous_state(self):
        marker = self.temporary / "shell-must-not-run"
        request = {
            "command": "test",
            "packages": ["gammalooprs"],
            "group": "core",
            "nextestArgs": [
                "a b",
                "'single'",
                '"double"',
                r"back\slash",
                f"$(touch {marker})",
                f"`touch {marker}`",
                "line\nbreak",
            ],
            "previousState": {"version": 1, "contexts": {"opaque": "x" * 200_000}},
        }
        files = []

        def run(command, **kwargs):
            path = Path(kwargs["env"]["GAMMALOOP_DEV_REQUEST_FILE"])
            files.append(path)
            text = path.read_text()
            self.assertEqual(json.loads(text), request)
            self.assertNotIn("fixture-license", text)
            self.assertFalse(kwargs.get("shell", False))
            self.assertTrue(kwargs["check"])
            self.assertNotIn("--apply", command)
            self.assertEqual(
                command[command.index("--expr") + 1], dev.REQUEST_BUILD_EXPRESSION
            )
            self.assertIn(
                'builtins.getEnv "GAMMALOOP_DEV_REQUEST_FILE"',
                dev.REQUEST_BUILD_EXPRESSION,
            )
            self.assertEqual(kwargs["env"]["GAMMALOOP_DEV_SOURCE"], SOURCE)
            self.assertEqual(kwargs["env"]["GAMMALOOP_DEV_SYSTEM"], "x86_64-linux")
            self.assertIn(
                'builtins.getEnv "GAMMALOOP_DEV_SOURCE"', dev.REQUEST_BUILD_EXPRESSION
            )
            self.assertIn(
                'builtins.getEnv "GAMMALOOP_DEV_SYSTEM"', dev.REQUEST_BUILD_EXPRESSION
            )
            self.assertNotIn(str(path), command)
            self.assertTrue(
                all("x" * 200_000 not in value for value in kwargs["env"].values())
            )
            self.assertNotIn("fixture-license", " ".join(command))
            return subprocess.CompletedProcess(command, 0)

        development = dev.NixDevelopment(self.temporary, SOURCE, "x86_64-linux")
        with (
            patch.dict(
                os.environ, {"SYMBOLICA_LICENSE": "fixture-license"}, clear=True
            ),
            patch.object(dev.subprocess, "run", side_effect=run) as mocked,
        ):
            development.run(request, self.temporary / "out link")
            development.run(request)
        linked, dry = [call.args[0] for call in mocked.call_args_list]
        self.assertEqual(linked[-2:], ["--out-link", str(self.temporary / "out link")])
        self.assertEqual(dry[-2:], ["--dry-run", "--no-link"])
        self.assertEqual(linked[:-2], dry[:-2])
        self.assertFalse(marker.exists())
        self.assertTrue(all(not path.exists() for path in files))

    def test_request_tempfile_is_removed_after_build_failure(self):
        paths = []

        def fail(command, **kwargs):
            paths.append(Path(kwargs["env"]["GAMMALOOP_DEV_REQUEST_FILE"]))
            self.assertTrue(paths[-1].is_file())
            raise subprocess.CalledProcessError(1, command)

        with (
            patch.object(dev.subprocess, "run", side_effect=fail),
            self.assertRaises(subprocess.CalledProcessError),
        ):
            dev.NixDevelopment(self.temporary, SOURCE, "x86_64-linux").run(
                {"command": "check"}
            )
        self.assertFalse(paths[0].exists())


class OutcomeTests(TemporaryTestCase):
    def test_junit_counts_namespaced_testcases_including_skips_not_summary(self):
        result = self.make_result("junit") / "result"
        (result / "tests.xml").write_text(
            '<testsuites xmlns="urn:junit" tests="999"><testsuite tests="999">'
            '<testcase name="pass"/><testcase name="skip"><skipped/></testcase>'
            "</testsuite></testsuites>"
        )
        (result / "more.xml").write_text(
            '<testsuite tests="999"><testcase/></testsuite>'
        )
        self.assertEqual(dev.junit_count(result), 3)

    def test_zero_requires_explicit_valid_xml_and_failures_are_rejected(self):
        result = self.make_result("empty", cases=0) / "result"
        self.assertEqual(dev.junit_count(result), 0)
        for xml in (
            "<broken",
            "<testsuite><testcase><failure/></testcase></testsuite>",
            '<testsuite xmlns="urn:junit"><error/></testsuite>',
        ):
            with self.subTest(xml=xml):
                (result / "tests.xml").write_text(xml)
                with self.assertRaises(dev.DevError):
                    dev.junit_count(result)
        (result / "tests.xml").unlink()
        with self.assertRaisesRegex(dev.DevError, "missing JUnit"):
            dev.junit_count(result)
        with self.assertRaisesRegex(dev.DevError, "missing JUnit"):
            dev.test_outcome([result], "pass")

    def test_empty_overall_policy_and_empty_plus_matched_groups(self):
        empty = self.make_result("empty", cases=0, status=4) / "result"
        matched = self.make_result("matched") / "result"
        for policy, expected in (("fail", 4), ("warn", 0), ("pass", 0)):
            with self.subTest(policy=policy), redirect_stderr(io.StringIO()):
                self.assertEqual(dev.test_outcome([empty], policy), expected)
                self.assertEqual(dev.test_outcome([empty, matched], policy), 0)

    def test_real_exit_status_wins_even_with_positive_cases_and_no_tests_pass(self):
        failure = self.make_result("failure", cases=1, status=101) / "result"
        empty = self.make_result("empty", cases=0, status=4) / "result"
        with redirect_stderr(io.StringIO()):
            self.assertEqual(dev.test_outcome([empty, failure], "pass"), 101)

    def test_exit_status_is_required_and_validated(self):
        result = self.make_result("status") / "result"
        for status in ("", "-1", "256", "0 trailing", "success"):
            with self.subTest(status=status):
                (result / "exit-status").write_text(status)
                with self.assertRaisesRegex(dev.DevError, "invalid.*exit status"):
                    dev.test_outcome([result], "pass")
        (result / "exit-status").unlink()
        with self.assertRaisesRegex(dev.DevError, "missing.*exit status"):
            dev.test_outcome([result], "pass")


class RegistryTests(TemporaryTestCase):
    def setUp(self):
        super().setUp()
        self.directory = self.temporary / "cache"
        self.registry = dev.CacheRegistry(self.directory, keep=2)

    def publish(self, name, *, context="core", source=SOURCE, successful=True):
        result = self.make_result(name, contexts={context: {"opaque": name}})
        with self.registry.generation() as generation:
            output = generation / "result-0"
            output.symlink_to(result, target_is_directory=True)
            self.registry.publish(generation, source, [output], successful=successful)
        return generation

    def test_registry_write_is_atomic_and_failed_replace_preserves_index(self):
        self.publish("first")
        index = self.directory / "state-index.json"
        before = index.read_bytes()
        with (
            patch.object(
                dev.os, "replace", side_effect=OSError("injected replace failure")
            ),
            self.assertRaisesRegex(OSError, "injected"),
        ):
            self.registry.write({"version": 1, "generations": []})
        self.assertEqual(index.read_bytes(), before)
        self.assertEqual(list(self.directory.glob(".index-*")), [])
        self.assertEqual(len(self.registry.read()["generations"]), 1)

    def test_registry_rejects_untrusted_generation_paths_and_invalid_state(self):
        for record in (
            {"id": "../unrelated", "state": {"version": 1, "contexts": {}}},
            {"id": "a" * 32, "state": {"version": 1, "contexts": []}},
        ):
            with self.subTest(record=record):
                self.registry.write({"version": 1, "generations": [record]})
                with self.assertRaises(dev.DevError):
                    self.registry.read()

    def test_sequential_retention_preserves_unregistered_directories(self):
        foreign = self.directory / "runs" / "unregistered"
        foreign.mkdir()
        (foreign / "keep").write_text("someone else's files")
        generations = [self.publish(f"result-{index}") for index in range(4)]
        self.assertEqual(
            [record["id"] for record in self.registry.read()["generations"]],
            [path.name for path in generations[-2:]],
        )
        self.assertTrue(all(not path.exists() for path in generations[:-2]))
        self.assertEqual((foreign / "keep").read_text(), "someone else's files")

    def test_contexts_retain_independent_recent_generations(self):
        self.registry = dev.CacheRegistry(self.directory, keep=1)
        core = self.publish("core-old")
        integration = self.publish("integration", context="integration")
        latest = self.publish("core-new")
        self.assertFalse(core.exists())
        self.assertTrue(integration.exists())
        self.assertTrue(latest.exists())
        with self.registry.previous() as previous:
            self.assertEqual(
                previous["contexts"],
                {
                    "core": {"opaque": "core-new"},
                    "integration": {"opaque": "integration"},
                },
            )

    def test_failed_results_do_not_displace_successful_context_history(self):
        self.registry = dev.CacheRegistry(self.directory, keep=1)
        successful = self.publish("successful")
        failed_old = self.publish("failed-old", successful=False)
        failed_new = self.publish("failed-new", successful=False)
        self.assertTrue(successful.exists())
        self.assertFalse(failed_old.exists())
        self.assertTrue(failed_new.exists())
        self.assertEqual(
            {
                record["id"]: record["successful"]
                for record in self.registry.read()["generations"]
            },
            {successful.name: True, failed_new.name: False},
        )

    def test_global_maximum_bounds_independent_context_families(self):
        self.registry = dev.CacheRegistry(self.directory, keep=2, maximum=2)
        generations = [
            self.publish(f"result-{index}", context=f"context-{index}")
            for index in range(3)
        ]
        self.assertFalse(generations[0].exists())
        self.assertTrue(all(path.exists() for path in generations[1:]))
        self.assertEqual(len(self.registry.read()["generations"]), 2)

    def test_max_age_expires_old_contexts_after_active_seed_lease_ends(self):
        self.registry = dev.CacheRegistry(self.directory, keep=3, max_age_days=1)
        with patch.object(dev.time, "time", return_value=1_000_000) as clock:
            expired = self.publish("expired", context="old-context")
            with self.registry.previous():
                clock.return_value += 86401
                fresh = self.publish("fresh", context="fresh-context")
                self.assertTrue(expired.exists())
            newest = self.publish("newest", context="new-context")
        self.assertFalse(expired.exists())
        self.assertTrue(fresh.exists())
        self.assertTrue(newest.exists())
        self.assertEqual(
            {record["id"] for record in self.registry.read()["generations"]},
            {fresh.name, newest.name},
        )

    def test_active_seed_lease_skips_pruning_and_registry_lock_is_released(self):
        self.registry = dev.CacheRegistry(self.directory, keep=1)
        first = self.publish("first")
        with self.registry.previous() as previous:
            with (self.directory / ".registry.lock").open("a+b") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            second = self.publish("second")
            self.assertTrue(first.is_dir())
            self.assertEqual(previous["contexts"], {"core": {"opaque": "first"}})
        third = self.publish("third")
        self.assertFalse(first.exists())
        self.assertFalse(second.exists())
        self.assertTrue(third.exists())

    def test_concurrent_publications_merge_disjoint_state_without_build_lock(self):
        barrier = threading.Barrier(2)
        results = [
            self.make_result(name, contexts={name: {"opaque": name}})
            for name in ("core", "integration")
        ]

        def publish(result):
            registry = dev.CacheRegistry(self.directory, keep=2)
            with registry.generation() as generation:
                output = generation / "result-0"
                output.symlink_to(result, target_is_directory=True)
                with registry.previous():
                    # Both simulated builds must enter before either publishes.
                    barrier.wait(timeout=5)
                registry.publish(generation, SOURCE, [output])

        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(publish, result) for result in results]
            for future in futures:
                future.result(timeout=10)
        self.assertEqual(len(self.registry.read()["generations"]), 2)
        with self.registry.previous() as previous:
            self.assertEqual(
                previous["contexts"],
                {"core": {"opaque": "core"}, "integration": {"opaque": "integration"}},
            )

    def test_warm_identical_outputs_keep_existing_root(self):
        result = self.make_result("immutable")
        with patch.object(dev.time, "time", return_value=1_000_000) as clock:
            with self.registry.generation() as first:
                output = first / "result-0"
                output.symlink_to(result, target_is_directory=True)
                self.registry.publish(first, SOURCE, [output])
            before = output.lstat()
            before_record = self.registry.read()["generations"][0]
            clock.return_value += 100
            with self.registry.generation() as duplicate:
                duplicate_output = duplicate / "result-0"
                duplicate_output.symlink_to(result, target_is_directory=True)
                self.registry.publish(duplicate, SOURCE, [duplicate_output])
        self.assertFalse(duplicate.exists())
        self.assertEqual(
            [record["id"] for record in self.registry.read()["generations"]],
            [first.name],
        )
        after = output.lstat()
        self.assertEqual(
            (after.st_ino, after.st_mtime_ns), (before.st_ino, before.st_mtime_ns)
        )
        self.assertEqual(output.resolve(), result)
        after_record = self.registry.read()["generations"][0]
        for key in ("id", "source", "state", "outputs"):
            self.assertEqual(after_record[key], before_record[key])
        self.assertGreater(after_record["lastUsed"], before_record["lastUsed"])

    def test_warm_dedup_refreshes_used_root_and_prunes_expired_unleased_overflow(self):
        self.registry = dev.CacheRegistry(
            self.directory, keep=1, maximum=1, max_age_days=1
        )
        with patch.object(dev.time, "time", return_value=1_000_000) as clock:
            first = self.publish("leased-first")
            with self.registry.previous():
                second = self.publish("warm-second")
                self.assertTrue(first.exists())
            self.assertEqual(len(self.registry.read()["generations"]), 2)
            clock.return_value += 86401
            with self.registry.generation() as duplicate:
                output = duplicate / "result-0"
                output.symlink_to(
                    self.temporary / "warm-second", target_is_directory=True
                )
                self.registry.publish(duplicate, SOURCE, [output])
        self.assertFalse(first.exists())
        self.assertFalse(duplicate.exists())
        self.assertTrue(second.exists())
        self.assertEqual(
            [record["id"] for record in self.registry.read()["generations"]],
            [second.name],
        )

    def test_unregistered_failed_invocation_cleans_only_its_own_generation(self):
        retained = self.publish("retained")
        with (
            self.assertRaisesRegex(RuntimeError, "build failed"),
            self.registry.generation() as abandoned,
        ):
            (abandoned / "partial").write_text("unfinished")
            raise RuntimeError("build failed")
        self.assertFalse(abandoned.exists())
        self.assertTrue(retained.exists())

    def test_private_and_provisioned_shared_modes_survive_restrictive_umask(self):
        shared = self.temporary / "shared-cache"
        shared.mkdir(mode=0o700)
        shared.chmod(0o2770)
        foreign = shared / "unrelated"
        foreign.mkdir(mode=0o751)
        foreign.chmod(0o751)
        (foreign / "keep").write_text("untouched")
        (foreign / "keep").chmod(0o604)
        result = self.make_result("mode-result")
        old_umask = os.umask(0o077)
        try:
            for directory, directory_mode, file_mode in (
                (self.temporary / "private-cache", 0o700, 0o600),
                (shared, 0o2770, 0o660),
            ):
                with self.subTest(directory=directory):
                    registry = dev.CacheRegistry(directory)
                    with registry.generation() as generation:
                        output = generation / "result-0"
                        output.symlink_to(result, target_is_directory=True)
                        registry.publish(generation, SOURCE, [output])
                    for path in (directory, directory / "runs", generation):
                        self.assertEqual(
                            stat.S_IMODE(path.stat().st_mode), directory_mode
                        )
                    for path in (
                        directory / ".registry.lock",
                        directory / "state-index.json",
                        generation / ".lease",
                    ):
                        self.assertEqual(stat.S_IMODE(path.stat().st_mode), file_mode)
        finally:
            os.umask(old_umask)
        self.assertEqual(stat.S_IMODE(foreign.stat().st_mode), 0o751)
        self.assertEqual(stat.S_IMODE((foreign / "keep").stat().st_mode), 0o604)


class ExecuteTests(TemporaryTestCase):
    def setUp(self):
        super().setUp()
        self.cache = self.temporary / "cache"
        self.outputs = {}
        self.calls = []
        self.enterContext(
            patch.dict(os.environ, {"GAMMALOOP_DEV_CACHE": str(self.cache)}, clear=True)
        )
        self.enterContext(
            patch.object(dev.SourceSnapshot, "store", return_value=SOURCE)
        )
        self.development = self.enterContext(
            patch.object(dev, "NixDevelopment", autospec=True)
        ).return_value
        self.development.catalog.return_value = CATALOG
        self.development.run.side_effect = self.fake_run

    def fake_run(self, request, out_link=None):
        self.calls.append((request, out_link))
        with (self.cache / ".registry.lock").open("a+b") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if out_link is not None:
            out_link.symlink_to(
                self.outputs[request.get("group", request["command"])],
                target_is_directory=True,
            )

    def execute(self, *arguments):
        stdout, stderr = io.StringIO(), io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            status = dev.execute(dev.parse_args(list(arguments)), self.temporary)
        return status, stdout.getvalue(), stderr.getvalue()

    def test_build_prints_immutable_product_and_warm_run_reuses_root_and_seed(self):
        result = self.make_result("build-output", packages=["gammaloop-api"])
        self.outputs["build"] = result
        first = self.execute("build", "-p", "gammaloop-api")
        registry = dev.CacheRegistry(self.cache)
        first_id = registry.read()["generations"][0]["id"]
        second = self.execute("build", "-p", "gammaloop-api")
        self.assertEqual(first[0], 0)
        self.assertEqual(second[0], 0)
        self.assertEqual(first[1], str(result / "result" / "gammaloop-api") + "\n")
        self.assertEqual(second[1], first[1])
        self.assertEqual(self.calls[0][0]["previousState"], {})
        self.assertEqual(
            self.calls[1][0]["previousState"],
            json.loads((result / "development-state.json").read_text()),
        )
        self.assertEqual(
            [record["id"] for record in registry.read()["generations"]], [first_id]
        )
        self.assertEqual(len(list((self.cache / "runs").iterdir())), 1)

    def test_dry_run_passes_seed_without_publishing_roots(self):
        result = self.make_result("seed")
        registry = dev.CacheRegistry(self.cache)
        with registry.generation() as generation:
            output = generation / "result-0"
            output.symlink_to(result, target_is_directory=True)
            registry.publish(generation, SOURCE, [output])
        before = registry.read()
        status, stdout, _ = self.execute("test", "--dry-run")
        self.assertEqual(status, 0)
        self.assertEqual(stdout, "")
        self.assertTrue(all(link is None for _, link in self.calls))
        self.assertTrue(
            all(
                request["previousState"]["contexts"] == {"core": {"opaque": "seed"}}
                for request, _ in self.calls
            )
        )
        self.assertEqual(registry.read(), before)
        self.assertEqual(len(list((self.cache / "runs").iterdir())), 1)

    def test_failed_materialized_tests_return_real_status_and_keep_diagnostics(self):
        self.outputs["core"] = self.make_result("failed-tests", status=101, logs=True)
        status, _, _ = self.execute("test", "-p", "gammalooprs", "--no-tests", "pass")
        self.assertEqual(status, 101)
        record = dev.CacheRegistry(self.cache).read()["generations"][0]
        self.assertFalse(record["successful"])
        root = self.cache / "runs" / record["id"] / "result-0"
        self.assertTrue(root.is_dir())
        self.assertEqual((root / "result" / "exit-status").read_text(), "101\n")
        self.assertTrue((root / "result" / "logs" / "trace.jsonl").is_file())

    def test_invalid_outcome_retains_unsuccessful_result_before_raising(self):
        self.outputs["core"] = self.make_result("invalid-tests")
        (self.outputs["core"] / "result" / "exit-status").unlink()
        with self.assertRaisesRegex(dev.DevError, "missing.*exit status"):
            self.execute("test", "-p", "gammalooprs")
        record = dev.CacheRegistry(self.cache).read()["generations"][0]
        self.assertFalse(record["successful"])
        self.assertTrue((self.cache / "runs" / record["id"] / "result-0").is_dir())

    def test_later_nix_failure_keeps_completed_failed_group_and_host_diagnostics(self):
        self.outputs["core"] = self.make_result(
            "completed-failed", status=101, logs=True
        )
        logs = self.temporary / "host-logs"

        def run(request, out_link=None):
            if request["group"] == "integration":
                raise subprocess.CalledProcessError(73, ["nix", "build"])
            return self.fake_run(request, out_link)

        self.development.run.side_effect = run
        with (
            patch.dict(os.environ, {"GL_TEST_LOG_DIR": str(logs)}),
            self.assertRaises(subprocess.CalledProcessError) as error,
        ):
            self.execute("test")
        self.assertEqual(error.exception.returncode, 73)
        registry = dev.CacheRegistry(self.cache)
        record = registry.read()["generations"][0]
        self.assertFalse(record["successful"])
        result = self.cache / "runs" / record["id"] / "result-0" / "result"
        self.assertEqual((result / "exit-status").read_text(), "101\n")
        self.assertTrue((result / "logs" / "trace.jsonl").is_file())
        self.assertEqual(
            (logs / record["id"] / "core").resolve(),
            self.outputs["core"] / "result" / "logs",
        )
        with registry.previous(successful_only=True) as previous:
            self.assertEqual(previous, {})

    def test_missing_or_malformed_junit_keeps_materialized_diagnostics(self):
        for name in ("missing-report", "malformed-report"):
            with self.subTest(name=name):
                wrapper = self.make_result(name, logs=True)
                self.outputs["core"] = wrapper
                report = wrapper / "result" / "tests.xml"
                if name == "missing-report":
                    report.unlink()
                else:
                    report.write_text("<not-valid")
                with self.assertRaisesRegex(dev.DevError, "JUnit"):
                    self.execute("test", "-p", "gammalooprs")
                record = dev.CacheRegistry(self.cache).read()["generations"][-1]
                self.assertFalse(record["successful"])
                root = self.cache / "runs" / record["id"] / "result-0"
                self.assertEqual(root.resolve(), wrapper)
                self.assertTrue((root / "result" / "logs" / "trace.jsonl").is_file())

    def test_build_exception_cleans_unregistered_invocation(self):
        registry = dev.CacheRegistry(self.cache)
        foreign = self.cache / "runs" / "someone-elses-run"
        foreign.mkdir()
        (foreign / "keep").write_text("preserved")
        self.development.run.side_effect = subprocess.CalledProcessError(
            1, ["nix", "build"]
        )
        with self.assertRaises(subprocess.CalledProcessError):
            self.execute("check", "-p", "gammalooprs")
        self.assertEqual(list((self.cache / "runs").iterdir()), [foreign])
        self.assertEqual(registry.read()["generations"], [])

    def test_host_log_links_are_unique_by_generation_and_group(self):
        logs = self.temporary / "host-logs"
        logs.mkdir()
        (logs / "core").write_text("preexisting caller file")
        self.outputs = {
            "core": self.make_result("core-tests", cases=0, status=4, logs=True),
            "integration": self.make_result("integration-tests", logs=True),
        }
        with patch.dict(os.environ, {"GL_TEST_LOG_DIR": str(logs)}):
            self.assertEqual(self.execute("test")[0], 0)
            self.assertEqual(self.execute("test")[0], 0)
        self.assertEqual((logs / "core").read_text(), "preexisting caller file")
        generations = [path for path in logs.iterdir() if path.is_dir()]
        self.assertEqual(len(generations), 2)
        for generation in generations:
            self.assertRegex(generation.name, r"^[0-9a-f]{32}$")
            for group, result in self.outputs.items():
                self.assertTrue((generation / group).is_symlink())
                self.assertEqual(
                    (generation / group).resolve(), result / "result" / "logs"
                )
        self.assertTrue(all(request["collectLogs"] for request, _ in self.calls))
        self.assertTrue(
            all(
                "GL_TEST_LOG_DIR" not in request["runtimeEnv"]
                for request, _ in self.calls
            )
        )


if __name__ == "__main__":
    unittest.main()
