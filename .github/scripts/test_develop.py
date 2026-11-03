"""Black-box contracts for the native shell entry point and its Nix helper."""

import json
import os
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
LITERAL = ["a b; $(touch never)", "`never`", "'quoted'\nsecond line", "--backend=cargo"]


class BootstrapTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        (self.root / ".config").mkdir()
        shutil.copyfile(
            REPOSITORY / ".config/build-profile", self.root / ".config/build-profile"
        )
        self.bin = self.root / "bin"
        self.bin.mkdir()
        scripts = self.root / "scripts"
        scripts.mkdir()
        for name in ("develop.sh", "develop_nix.sh"):
            shutil.copy2(REPOSITORY / "scripts" / name, scripts / name)
        self.script = scripts / "develop.sh"
        self.bash = shutil.which("bash")
        self.log = self.root / "commands.jsonl"
        self.tools = self.root / "pinned"
        (self.tools / "bin").mkdir(parents=True)
        self.output = self.root / "nix-output"
        (self.output / "bin").mkdir(parents=True)
        self.home = self.root / "home"
        self.home.mkdir()
        self.env = {
            "PATH": str(self.bin),
            "HOME": str(self.home),
            "TEST_LOG": str(self.log),
            "TEST_TOOLS": str(self.tools),
            "TEST_OUTPUT": str(self.output),
        }
        # No host Python, jq, or Nix on PATH. Stubs use the harness interpreter
        # absolutely; production shell commands cannot discover it by name.
        for name in ("dirname", "mktemp", "rm", "cat"):
            (self.bin / name).symlink_to(shutil.which(name))
        for name in (
            "bash",
            "cat",
            "chmod",
            "cmp",
            "cp",
            "dirname",
            "find",
            "flock",
            "jq",
            "mkdir",
            "mktemp",
            "mv",
            "rm",
            "sort",
            "tar",
            "tr",
            "uname",
        ):
            target = shutil.which(name)
            if target:
                (self.tools / "bin" / name).symlink_to(target)
        self.stub("cargo")

    def stub(self, name, body="", *, directory=None):
        executable = (directory or self.bin) / name
        executable.write_text(
            f"#!{sys.executable}\n"
            "import json, os, sys, time, signal, fcntl\n"
            "from pathlib import Path\n"
            "record = dict(program=Path(sys.argv[0]).name, args=sys.argv[1:], "
            "cwd=os.getcwd(), root=os.environ.get('GAMMALOOP_WORKSPACE_ROOT'), "
            "license='SYMBOLICA_LICENSE' in os.environ, "
            "signed='SYMBOLICA_LICENSE_SIGNED' in os.environ, "
            "stale_args='GAMMALOOP_NIX_TEST_ARGS' in os.environ, "
            "env={k:v for k,v in os.environ.items() if k.startswith('GAMMALOOP_') "
            "and k != 'GAMMALOOP_NIX_TEST_ARGS' or k in "
            "('NO_SYMBOLICA_OEM_LICENSE','RUST_MIN_STACK','RAYON_NUM_THREADS',"
            "'NEXTEST_TEST_THREADS')})\n"
            "with open(os.environ['TEST_LOG'], 'a') as log:\n"
            "    log.write(json.dumps(record) + '\\n')\n"
            + body
            + "\nsys.exit(int(os.environ.get('TEST_EXIT', '0')))\n"
        )
        executable.chmod(0o755)
        return executable

    def records(self):
        if not self.log.exists():
            return []
        return [json.loads(line) for line in self.log.read_text().splitlines()]

    def invoke(self, arguments):
        self.log.unlink(missing_ok=True)
        result = subprocess.run(
            [self.bash, str(self.script), *arguments],
            cwd=self.bin,
            env=self.env,
            text=True,
            capture_output=True,
            check=False,
            timeout=20,
        )
        records = self.records()
        for record in records:
            self.assertEqual(record["cwd"], str(self.root))
            self.assertEqual(record["root"], str(self.root))
        return result, records

    def cargo(self, arguments, *, action="test"):
        result, records = self.invoke(["--backend=cargo", action, *arguments])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([r["program"] for r in records], ["cargo"])
        return records[0]

    def enable_nix(self, *, vcs="git"):
        self.assertTrue((self.tools / "bin/jq").exists(), "run in nix develop")
        (self.root / (".jj" if vcs == "jj" else ".git")).mkdir()
        self.env["TEST_FILES"] = json.dumps(["tracked.txt", "deleted.txt"])
        (self.root / "tracked.txt").write_text("live uncommitted content")
        (self.root / "flake.lock").write_text("{}")
        self.stub(
            vcs,
            "sys.stdout.buffer.write(b''.join(p.encode() + b'\\0' for p in json.loads(os.environ['TEST_FILES'])))\nsys.exit(0)\n",
        )
        self.stub(
            "nix",
            "args = sys.argv[1:]\n"
            "if args[0] == 'develop': sys.exit(0)\n"
            "if '--expr' not in args and os.environ.get('TEST_FAIL_BEFORE_LINK'):\n"
            "    sys.exit(int(os.environ['TEST_BUILD_EXIT']))\n"
            "output = Path(os.environ['TEST_TOOLS'] if '--expr' in args else os.environ['TEST_OUTPUT'])\n"
            "if '--out-link' in args:\n"
            "    link = Path(args[args.index('--out-link') + 1])\n"
            "    link.parent.mkdir(parents=True, exist_ok=True)\n"
            "    link.unlink(missing_ok=True)\n"
            "    link.symlink_to(output)\n"
            "if '--expr' not in args:\n"
            "    source = Path(os.environ['GAMMALOOP_FLAKE_SOURCE'])\n"
            "    snapshot = {str(p.relative_to(source)): "
            "({'link':os.readlink(p)} if p.is_symlink() else p.read_text()) "
            "for p in source.rglob('*') if p.is_symlink() or p.is_file()}\n"
            "    Path(os.environ['TEST_LOG'] + '.snapshot').write_text(json.dumps(snapshot))\n"
            "    sys.exit(int(os.environ.get('TEST_BUILD_EXIT', '0'))) if os.environ.get('TEST_BUILD_EXIT') else None\n"
            "print(output)\n"
            "sys.exit(0)\n",
        )
        self.stub(
            "nix-store",
            "args = sys.argv[1:]\n"
            "root = Path(args[args.index('--add-root') + 1])\n"
            "if os.environ.get('TEST_ROOT_EXIT') and '/bin/' in str(root):\n"
            "    sys.exit(int(os.environ['TEST_ROOT_EXIT']))\n"
            "if '/bin/' in str(root):\n"
            "    lock = root.parents[2] / ('.install-' + root.parent.name + '.lock')\n"
            "    with lock.open('a') as handle:\n"
            "        try: fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)\n"
            "        except BlockingIOError: pass\n"
            "        else: raise AssertionError('binary root committed without install lock')\n"
            "    time.sleep(float(os.environ.get('TEST_ROOT_DELAY', '0')))\n"
            "root.parent.mkdir(parents=True, exist_ok=True)\n"
            "root.unlink(missing_ok=True)\n"
            "root.symlink_to(args[args.index('--realise') + 1])\n"
            "sys.exit(0)\n",
        )
        self.stub(
            "nix-test-selected",
            "reports = Path(sys.argv[1])\n"
            "assert reports.is_dir()\n"
            "for call in map(json.loads, Path(os.environ['TEST_LOG']).read_text().splitlines()):\n"
            "    if call['program'] == 'nix' and '--out-link' in call['args']:\n"
            "        assert Path(call['args'][call['args'].index('--out-link') + 1]).exists()\n"
            "assert Path(os.environ['GAMMALOOP_FLAKE_SOURCE']).is_dir()\n"
            "(reports / 'runtime-result').write_text(os.environ.get('TEST_EXIT', '0'))\n"
            "if os.environ.get('TEST_WAIT'):\n"
            "    Path(os.environ['TEST_WAIT']).touch()\n"
            "    time.sleep(60)\n",
            directory=self.output / "bin",
        )

    def assert_cleaned(self, records):
        for record in records:
            source = record["env"].get("GAMMALOOP_FLAKE_SOURCE")
            if source:
                self.assertFalse(Path(source).exists(), source)
            if record["program"] == "nix" and "--expr" in record["args"]:
                link = Path(record["args"][record["args"].index("--out-link") + 1])
                self.assertFalse(link.parent.exists(), str(link.parent))
        self.assertEqual(list((self.root / "target/nix").glob(".pending-*")), [])

    def test_all_cargo_actions_need_neither_python_nor_jq_nor_nix(self):
        for action in (
            "build",
            "test",
            "check",
            "clippy",
            "doc",
            "run",
            "list",
            "fmt",
            "clean",
        ):
            with self.subTest(action=action):
                args = ["--workspace"] if action in ("test", "list") else []
                result, records = self.invoke([action, *args])
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual([r["program"] for r in records], ["cargo"])
                self.assertEqual(result.stderr, "Backend: cargo\n")
                expected = (
                    ["nextest", "run" if action == "test" else "list"]
                    if action in ("test", "list")
                    else [action]
                )
                self.assertEqual(records[0]["args"][: len(expected)], expected)

    def test_shared_profile_file_controls_cargo_and_nix_defaults(self):
        (self.root / ".config/build-profile").write_text("dev\n")
        result, records = self.invoke(["build"])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            records[0]["args"][records[0]["args"].index("--profile") + 1], "dev"
        )
        self.enable_nix()
        result, records = self.invoke(["check"])
        self.assertEqual(result.returncode, 0, result.stderr)
        request = next(
            record
            for record in records
            if record["program"] == "nix" and "--expr" not in record["args"]
        )
        self.assertEqual(
            json.loads(request["env"]["GAMMALOOP_DEVELOP_SPEC"])["profile"], "dev"
        )

    def test_manual_profile_is_shared_by_cargo_actions(self):
        for profile in ("dev", "release"):
            for action in ("build", "check", "clippy", "doc", "run", "test", "list"):
                with self.subTest(profile=profile, action=action):
                    self.env["GAMMALOOP_BUILD_PROFILE"] = profile
                    arguments = ["--workspace"] if action in ("test", "list") else []
                    result, records = self.invoke([action, *arguments])
                    self.assertEqual(result.returncode, 0, result.stderr)
                    flag = (
                        "--cargo-profile" if action in ("test", "list") else "--profile"
                    )
                    self.assertEqual(
                        records[0]["args"][records[0]["args"].index(flag) + 1], profile
                    )

    def test_just_profile_selection_reaches_native_and_python_recipes(self):
        (self.bin / "bash").symlink_to(self.bash)
        (self.bin / "sh").symlink_to(shutil.which("sh"))
        self.stub("python3")
        just_root = self.root / "recipes"
        for name in (
            "justfile",
            "just/ci.just",
            "just/tests.just",
            "just/drawing.just",
            "crates/linnet/Justfile",
            "tydenso/Justfile",
        ):
            target = just_root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            # A Nix build sandbox has no /usr/bin/env. Preserve the actual
            # recipes while resolving their Bash interpreter explicitly.
            target.write_text(
                (REPOSITORY / name)
                .read_text()
                .replace("#!/usr/bin/env bash", f"#!{self.bash}")
            )
        for name in (
            ".config/build-profile",
            "scripts/develop.sh",
            "scripts/develop_nix.sh",
        ):
            target = just_root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(REPOSITORY / name, target)
        (just_root / ".config/build-profile").write_text("dev\n")
        result = subprocess.run(
            [
                shutil.which("just"),
                "--justfile",
                str(just_root / "justfile"),
                "--working-directory",
                str(self.root),
                "build",
                "kurvst",
            ],
            env=self.env,
            text=True,
            capture_output=True,
            check=False,
            timeout=20,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        record = self.records()[-1]
        self.assertEqual(record["args"][record["args"].index("--profile") + 1], "dev")
        for profile in ("dev", "release"):
            for recipe in (
                "build",
                "build-cli",
                "build-cli-no-pyo3",
                "build-all",
                "doc",
                "clippy",
                "test",
                "python-wheel",
                "build-api",
            ):
                with self.subTest(profile=profile, recipe=recipe):
                    self.log.unlink(missing_ok=True)
                    arguments = ["kurvst"] if recipe in ("build", "test") else []
                    result = subprocess.run(
                        [
                            shutil.which("just"),
                            "--justfile",
                            str(just_root / "justfile"),
                            "--working-directory",
                            str(self.root),
                            f"profile={profile}",
                            recipe,
                            *arguments,
                        ],
                        env=self.env,
                        text=True,
                        capture_output=True,
                        check=False,
                        timeout=20,
                    )
                    self.assertEqual(result.returncode, 0, result.stderr)
                    record = self.records()[-1]
                    self.assertEqual(record["env"]["GAMMALOOP_BUILD_PROFILE"], profile)
                    if record["program"] == "cargo":
                        flag = "--cargo-profile" if recipe == "test" else "--profile"
                        self.assertEqual(
                            record["args"][record["args"].index(flag) + 1], profile
                        )
                    else:
                        self.assertIn("gammaloop", record["args"])

            self.log.unlink(missing_ok=True)
            result = subprocess.run(
                [
                    shutil.which("just"),
                    "--justfile",
                    str(just_root / "justfile"),
                    f"profile={profile}",
                    "linnet",
                    "build",
                ],
                env=self.env,
                text=True,
                capture_output=True,
                check=False,
                timeout=20,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            record = self.records()[-1]
            self.assertEqual(
                record["args"][record["args"].index("--profile") + 1], profile
            )

    def test_nix_requests_use_the_shared_profile(self):
        self.enable_nix()
        for profile in ("dev", "release"):
            self.env["GAMMALOOP_BUILD_PROFILE"] = profile
            for action in ("check", "test"):
                with self.subTest(profile=profile, action=action):
                    result, records = self.invoke([action, "-p", "kurvst"])
                    self.assertEqual(result.returncode, 0, result.stderr)
                    request = next(
                        record
                        for record in records
                        if record["program"] == "nix" and "--expr" not in record["args"]
                    )
                    if action == "test":
                        self.assertEqual(
                            request["env"]["GAMMALOOP_NIX_TEST_PROFILE"], profile
                        )
                    else:
                        self.assertEqual(
                            json.loads(request["env"]["GAMMALOOP_DEVELOP_SPEC"])[
                                "profile"
                            ],
                            profile,
                        )

    def test_explicit_cargo_never_uses_present_nix(self):
        self.stub("nix", "raise AssertionError('explicit Cargo consulted Nix')\n")
        for prefix in (["--backend", "cargo"], ["--backend=cargo"]):
            result, records = self.invoke([*prefix, "check", *LITERAL])
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual([r["program"] for r in records], ["cargo"])
            self.assertEqual(records[0]["args"][-len(LITERAL) :], LITERAL)

    def test_backend_prefix_errors_precede_tools(self):
        self.stub("nix")
        for prefix in (
            ["--backend"],
            ["--backend="],
            ["--backend=bad"],
            ["--backend", "bad"],
            ["--backend=auto", "--backend=cargo"],
            ["--backend", "cargo", "--backend", "nix"],
        ):
            with self.subTest(prefix=prefix):
                result, records = self.invoke(
                    prefix
                    if len(prefix) == 1 and prefix[0] == "--backend"
                    else [*prefix, "check"]
                )
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertEqual(records, [])

    def test_explicit_missing_nix_does_not_fall_back(self):
        result, records = self.invoke(["--backend=nix", "check"])
        self.assertEqual(result.returncode, 2)
        self.assertEqual(records, [])
        self.assertIn("not installed", result.stderr)

    def test_native_defaults_only_when_absent_and_literal_tail(self):
        for action in ("build", "check", "clippy", "doc", "run"):
            with self.subTest(action=action):
                command = self.cargo([], action=action)["args"]
                self.assertEqual(
                    command,
                    [
                        action,
                        "--locked",
                        "--profile",
                        "ci-optim",
                        "-p",
                        "gammaloop-api",
                    ],
                )
                explicit = ["-p", "spenso", "--profile=dev", "--", *LITERAL]
                command = self.cargo(explicit, action=action)["args"]
                self.assertEqual(command, [action, "--locked", *explicit])
                command = self.cargo(["--release", "--workspace"], action=action)[
                    "args"
                ]
                self.assertNotIn("--profile", command)
                self.assertNotIn("gammaloop-api", command)
        self.assertFalse((self.root / "never").exists())

    def test_options_after_separator_do_not_suppress_defaults(self):
        command = self.cargo(
            ["--", "--release", "--package=other", *LITERAL], action="run"
        )["args"]
        self.assertEqual(
            command[:6],
            ["run", "--locked", "--profile", "ci-optim", "-p", "gammaloop-api"],
        )
        self.assertEqual(command[6:], ["--", "--release", "--package=other", *LITERAL])

    def test_one_workspace_nextest_call_and_complete_report_on_success_and_failure(
        self,
    ):
        self.stub(
            "cargo",
            "report = Path('target/nextest/ci_gammaloop/junit.xml')\nreport.parent.mkdir(parents=True, exist_ok=True)\nreport.write_text(os.environ['TEST_JUNIT'])\n",
        )
        report = self.root / "target/nextest/ci_gammaloop/junit.xml"
        self.env["TEST_JUNIT"] = (
            '<testsuites tests="2"><testsuite name="spenso"/><testsuite name="kurvst"/></testsuites>'
        )
        for status in (0, 31):
            with self.subTest(status=status):
                self.env["TEST_EXIT"] = str(status)
                result, records = self.invoke(["test", "--workspace"])
                self.assertEqual(result.returncode, status, result.stderr)
                self.assertEqual(len(records), 1)
                command = records[0]["args"]
                self.assertEqual(
                    command[:6],
                    [
                        "nextest",
                        "run",
                        "--locked",
                        "--cargo-profile",
                        "ci-optim",
                        "--workspace",
                    ],
                )
                self.assertNotIn("--features", command)
                self.assertEqual(records[0]["env"]["NO_SYMBOLICA_OEM_LICENSE"], "1")
                self.assertEqual(report.read_text(), self.env["TEST_JUNIT"])

    def test_nextest_preserves_native_selectors_features_and_literal_arguments(self):
        selectors = [
            "--lib",
            "--test",
            "tensor",
            "--features",
            "new-member/extra spenso/new-mode",
            "--all-features",
            "--no-default-features",
        ]
        tail = ["--", "--workspace", "--test", "--profile", *LITERAL]
        record = self.cargo(["-p", "new-member", "-p", "spenso", *selectors, *tail])
        command = record["args"]
        self.assertEqual(command[5:9], ["-p", "new-member", "-p", "spenso"])
        for item in selectors:
            self.assertIn(item, command)
        self.assertEqual(command[-len(tail) :], tail)
        self.assertNotIn("python_api", command)
        self.assertFalse((self.root / "never").exists())

    def test_nextest_release_and_runtime_profile_are_independent(self):
        record = self.cargo(
            [
                "--workspace",
                "--cargo-profile",
                "dev",
                "--release",
                "--profile=ci_full",
                "--no-tests=warn",
                "-E",
                LITERAL[0],
            ]
        )
        command = record["args"]
        self.assertEqual(command[command.index("--cargo-profile") + 1], "release")
        self.assertEqual(command[command.index("--profile") + 1], "ci_full")
        self.assertEqual(command[command.index("--no-tests") + 1], "warn")
        self.assertIn(LITERAL[0], command[command.index("-E") + 1])

    def test_ci_slow_environment_defaults_preserve_explicit_threads(self):
        for threads, expected in ((None, "16"), ("3", "3")):
            with self.subTest(threads=threads):
                self.env["NEXTEST_TEST_THREADS"] = "16"
                if threads:
                    self.env["RAYON_NUM_THREADS"] = threads
                record = self.cargo(["--workspace", "-Pci_slow"])
                self.assertEqual(record["env"]["RUST_MIN_STACK"], "134217728")
                self.assertEqual(record["env"]["RAYON_NUM_THREADS"], expected)
                self.assertEqual(record["env"]["GAMMALOOP_TESTS_NO_CLEAN_STATE"], "1")

    def test_cargo_failure_status_never_retries(self):
        self.env["TEST_EXIT"] = "41"
        for action in ("build", "test", "list", "run", "fmt", "clean"):
            args = ["--workspace"] if action in ("test", "list") else []
            result, records = self.invoke([action, *args])
            self.assertEqual(result.returncode, 41, result.stderr)
            self.assertEqual(len(records), 1)

    def test_list_preserves_nextest_json_stdout(self):
        inventory = {"test-count": 2, "rust-suites": {"spenso": {}, "kurvst": {}}}
        self.env["TEST_INVENTORY"] = json.dumps(inventory)
        self.stub("cargo", "print(os.environ['TEST_INVENTORY'])\n")
        result, records = self.invoke(["list", "--workspace"])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout), inventory)
        command = records[0]["args"]
        self.assertEqual(command[:2], ["nextest", "list"])
        self.assertEqual(command[command.index("--message-format") + 1], "json")
        self.assertNotIn("--no-tests", command)

    def test_runtime_license_precedence_and_no_secret_output(self):
        config = self.home / ".config/symbolica/license"
        config.parent.mkdir(parents=True)
        config.write_text("config-sentinel")
        self.stub(
            "cargo",
            "assert os.environ['SYMBOLICA_LICENSE'] == os.environ['TEST_EXPECT_LICENSE']\n",
        )
        for exported, signed, expected in (
            ("export-sentinel", "signed-sentinel", "export-sentinel"),
            ("", "signed-sentinel", "signed-sentinel"),
            ("", "", "config-sentinel"),
        ):
            self.env.update(
                SYMBOLICA_LICENSE=exported,
                SYMBOLICA_LICENSE_SIGNED=signed,
                TEST_EXPECT_LICENSE=expected,
            )
            result, records = self.invoke(["test", "--workspace"])
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(records[0]["license"])
            self.assertNotIn(
                "sentinel", result.stdout + result.stderr + self.log.read_text()
            )
        alternate = self.root / "xdg/symbolica/license"
        alternate.parent.mkdir(parents=True)
        alternate.write_text("xdg-sentinel")
        self.env.update(
            XDG_CONFIG_HOME=str(alternate.parent.parent),
            TEST_EXPECT_LICENSE="xdg-sentinel",
        )
        result, _ = self.invoke(["run"])
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_nix_bootstrap_and_build_without_host_python_or_license(self):
        self.enable_nix()
        self.env.update(
            SYMBOLICA_LICENSE="private-sentinel",
            SYMBOLICA_LICENSE_SIGNED="signed-sentinel",
            GAMMALOOP_NIX_TEST_ARGS='["stale"]',
        )
        result, records = self.invoke(
            ["--backend=auto", "test", "-p", "kurvst", "--", *LITERAL]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(records[0]["program"], "nix")
        expression = records[0]["args"][records[0]["args"].index("--expr") + 1]
        self.assertIn("flake.lock", expression)
        self.assertNotIn("python", expression.lower())
        self.assertEqual(records[-1]["program"], "nix-test-selected")
        self.assertEqual(records[-1]["args"][-len(LITERAL) :], LITERAL)
        self.assertTrue(records[-1]["license"])
        for record in records:
            if record["program"] in ("nix", "nix-store"):
                self.assertFalse(record["license"])
                self.assertFalse(record["signed"])
                self.assertFalse(record["stale_args"])
        self.assertNotIn(
            "sentinel", result.stdout + result.stderr + self.log.read_text()
        )
        self.assert_cleaned(records)

    def test_nix_failure_never_falls_back_and_cleans_source_and_roots(self):
        self.enable_nix()
        self.env["TEST_BUILD_EXIT"] = "47"
        result, records = self.invoke(["check"])
        self.assertEqual(result.returncode, 47, result.stderr)
        self.assertNotIn("cargo", [r["program"] for r in records])
        self.assert_cleaned(records)

    def test_failed_static_commands_preserve_previous_permanent_roots(self):
        self.enable_nix()
        previous = self.root / "previous-output"
        previous.mkdir()
        self.env.update(TEST_BUILD_EXIT="47", TEST_FAIL_BEFORE_LINK="1")
        for action in ("check", "clippy", "doc"):
            with self.subTest(action=action):
                permanent = self.root / f"target/nix/{action}-gammaloop-api-ci-optim"
                permanent.parent.mkdir(parents=True, exist_ok=True)
                permanent.symlink_to(previous)
                result, records = self.invoke([action])
                self.assertEqual(result.returncode, 47, result.stderr)
                self.assertTrue(permanent.is_symlink())
                self.assertEqual(permanent.resolve(), previous)
                self.assert_cleaned(records)

    def test_python_helpers_can_resolve_their_bundled_uv_without_host_tools(self):
        self.enable_nix()
        self.stub("uv", directory=self.output / "bin")
        self.stub(
            "python3",
            "import shutil, subprocess\n"
            "assert shutil.which('uv') == str(Path(sys.argv[0]).with_name('uv'))\n"
            "subprocess.run(['uv', '--version'], check=True)\n",
            directory=self.output / "bin",
        )
        for arguments in (("python", "gammaloop", "--install"), ("upload-preview",)):
            with self.subTest(arguments=arguments):
                result, records = self.invoke(arguments)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(
                    [record["program"] for record in records][-2:], ["python3", "uv"]
                )
                self.assert_cleaned(records)

    def test_failed_bootstrap_cleans_temporary_tools(self):
        self.stub("nix", "sys.exit(39)\n")
        result, records = self.invoke(["check"])
        self.assertEqual(result.returncode, 39, result.stderr)
        self.assertEqual(len(records), 1)
        self.assert_cleaned(records)

    def test_nix_runtime_reports_and_protection_survive_each_status(self):
        self.enable_nix()
        other = self.root / "target/nix/.pending-tests-other"
        other.parent.mkdir(parents=True)
        other.symlink_to(self.root / "other-output")
        reports = []
        for status in (0, 4, 31):
            self.env["TEST_EXIT"] = str(status)
            result, records = self.invoke(["test", "-p", "kurvst"])
            self.assertEqual(result.returncode, status, result.stderr)
            report = Path(records[-1]["args"][0])
            reports.append(report)
            self.assertEqual((report / "runtime-result").read_text(), str(status))
            self.assertIn(str(report), result.stderr)
            self.assertTrue(other.is_symlink())
            self.assertEqual(list(other.parent.glob(".pending-*")), [other])
            other.unlink()
            self.assert_cleaned(records)
            other.symlink_to(self.root / "other-output")
        self.assertEqual(len(set(reports)), 3)

    def test_nix_selection_forwards_targets_and_filters_to_runner(self):
        self.enable_nix()
        forwarded = [
            "--test",
            "tensor",
            "--lib",
            "-P",
            "ci_slow",
            "-E",
            LITERAL[0],
            "--",
            *LITERAL,
        ]
        result, records = self.invoke(
            [
                "test",
                "-p",
                "spenso",
                "--cargo-profile=dev-optim",
                "--features",
                "spenso/shadowing",
                *forwarded,
            ]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        build = next(
            r for r in records if r["program"] == "nix" and "--expr" not in r["args"]
        )
        self.assertEqual(
            json.loads(build["env"]["GAMMALOOP_NIX_TEST_PACKAGES"]), ["spenso"]
        )
        self.assertEqual(
            json.loads(build["env"]["GAMMALOOP_NIX_TEST_FEATURES"]),
            ["spenso/shadowing"],
        )
        self.assertEqual(build["env"]["GAMMALOOP_NIX_TEST_PROFILE"], "dev-optim")
        self.assertEqual(records[-1]["args"][1:], forwarded)

    def test_nix_workspace_uses_default_selection_without_generated_metadata(self):
        self.enable_nix()
        result, records = self.invoke(["test", "--workspace", "--release"])
        self.assertEqual(result.returncode, 0, result.stderr)
        build = next(
            r for r in records if r["program"] == "nix" and "--expr" not in r["args"]
        )
        self.assertNotIn("GAMMALOOP_NIX_TEST_PACKAGES", build["env"])
        self.assertEqual(build["env"]["GAMMALOOP_NIX_TEST_PROFILE"], "release")

    def test_snapshots_capture_live_files_exclude_outputs_and_preserve_symlinks(self):
        for vcs in ("git", "jj"):
            with self.subTest(vcs=vcs):
                # Separate fixture roots make both discovery routes independent.
                self.enable_nix(vcs=vcs)
                (self.root / "dir").mkdir(exist_ok=True)
                (self.root / "dir/space\nname").write_text("live newline filename")
                (self.root / "--option-like-name").write_text("literal leading dash")
                (self.root / "link").symlink_to("tracked.txt")
                for directory in ("target", ".venv"):
                    (self.root / directory).mkdir(exist_ok=True)
                    (self.root / directory / "ignored").write_text("excluded")
                self.env["TEST_FILES"] = json.dumps(
                    [
                        "tracked.txt",
                        "deleted.txt",
                        "dir/space\nname",
                        "--option-like-name",
                        "link",
                        "target/ignored",
                        ".venv/ignored",
                    ]
                )
                result, records = self.invoke(["check"])
                self.assertEqual(result.returncode, 0, result.stderr)
                snapshot = json.loads(Path(str(self.log) + ".snapshot").read_text())
                self.assertEqual(
                    snapshot,
                    {
                        "tracked.txt": "live uncommitted content",
                        "dir/space\nname": "live newline filename",
                        "--option-like-name": "literal leading dash",
                        "link": {"link": "tracked.txt"},
                    },
                )
                self.assert_cleaned(records)
                (self.root / "link").unlink()
                (self.root / (".jj" if vcs == "jj" else ".git")).rmdir()

    def test_snapshot_rejects_traversal_and_directory_symlinks(self):
        self.enable_nix()
        outside = self.root / "outside"
        outside.mkdir()
        (outside / "file").write_text("must not copy")
        (self.root / "alias").symlink_to(outside, target_is_directory=True)
        for name in ("../escape", "/absolute", "alias/file"):
            self.env["TEST_FILES"] = json.dumps([name])
            result, records = self.invoke(["check"])
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(sum(r["program"] == "nix" for r in records), 1)
            self.assert_cleaned(records)

    def test_git_snapshot_uses_index_paths_and_live_worktree_contents(self):
        self.enable_nix()
        (self.root / ".git").rmdir()
        (self.bin / "git").unlink()
        git = shutil.which("git")
        self.assertIsNotNone(git, "Git is a declared test input")
        (self.bin / "git").symlink_to(git)

        def git_run(*arguments):
            subprocess.run(
                [git, *arguments],
                cwd=self.root,
                env=self.env,
                capture_output=True,
                check=True,
                timeout=10,
            )

        git_run("init", "--quiet")
        (self.root / "deleted.txt").write_text("deleted after staging")
        git_run("add", "tracked.txt", "deleted.txt")
        (self.root / "tracked.txt").write_text("edited after staging")
        (self.root / "deleted.txt").unlink()
        (self.root / "untracked.txt").write_text("not a source input")
        result, records = self.invoke(["check"])
        self.assertEqual(result.returncode, 0, result.stderr)
        snapshot = json.loads(Path(str(self.log) + ".snapshot").read_text())
        self.assertEqual(snapshot, {"tracked.txt": "edited after staging"})
        self.assert_cleaned(records)

    def test_nix_invalid_selection_fails_before_artifact_build(self):
        self.enable_nix()
        for arguments in (
            ["test"],
            ["test", "--workspace", "-p", "spenso"],
            ["build", "-p", "spenso", "-p", "kurvst"],
        ):
            with self.subTest(arguments=arguments):
                result, records = self.invoke(arguments)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertEqual(sum(r["program"] == "nix" for r in records), 1)
                self.assert_cleaned(records)

    def test_nix_native_spec_preserves_features_profile_and_literal_clippy_tail(self):
        self.enable_nix()
        tail = ["--", "-D", "warnings", *LITERAL]
        result, records = self.invoke(
            [
                "clippy",
                "-p",
                "spenso",
                "--profile",
                "dev",
                "--release",
                "--features",
                "z,a",
                "-F",
                "a,b",
                "--no-default-features",
                *tail,
            ]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        build = next(
            r for r in records if r["program"] == "nix" and "--expr" not in r["args"]
        )
        self.assertTrue(build["args"][-1].endswith("#development-clippy"))
        self.assertEqual(
            json.loads(build["env"]["GAMMALOOP_DEVELOP_SPEC"]),
            {
                "package": "spenso",
                "profile": "release",
                "features": ["a", "b", "z"],
                "default_features": False,
                "args": tail,
            },
        )
        self.assert_cleaned(records)

    def test_workspace_doc_open_uses_the_default_package_index(self):
        self.enable_nix()
        self.stub("xdg-open")
        result, records = self.invoke(["doc", "--workspace", "--open"])
        self.assertEqual(result.returncode, 0, result.stderr)
        opened = [record for record in records if record["program"] == "xdg-open"]
        self.assertEqual(len(opened), 1)
        self.assertEqual(
            opened[0]["args"],
            [str(self.output / "share/doc/gammaloop_api/index.html")],
        )

    def test_fmt_and_clean_use_nix_tools_without_native_artifact_build(self):
        self.enable_nix()
        for action in ("fmt", "clean"):
            result, records = self.invoke([action, *LITERAL])
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual([r["program"] for r in records], ["nix", "git", "nix"])
            self.assertEqual(records[-1]["args"][0], "develop")
            self.assertEqual(records[-1]["args"][-len(LITERAL) :], LITERAL)
            self.assert_cleaned(records)

    def binary(self, name="gammaloop", content="new executable", output=None):
        binary = (output or self.output) / "bin" / name
        binary.parent.mkdir(parents=True, exist_ok=True)
        binary.write_text(content)
        binary.chmod(0o555)
        os.utime(binary, (1, 1))
        return binary

    def test_binary_install_profiles_mtime_permissions_and_read_only_replacement(self):
        self.enable_nix()
        binary = self.binary()
        for profile, directory in (
            ("dev", "debug"),
            ("release", "release"),
            ("dev-optim", "dev-optim"),
            ("ci-optim", "ci-optim"),
        ):
            destination = self.root / "target" / directory / "gammaloop"
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text("old executable")
            destination.chmod(0o555)
            started = time.time() - 1
            with destination.open("rb") as old:
                result, records = self.invoke(["build", "--profile", profile])
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(old.read(), b"old executable")
            self.assertFalse(destination.is_symlink())
            self.assertEqual(destination.read_bytes(), binary.read_bytes())
            self.assertEqual(stat.S_IMODE(destination.stat().st_mode), 0o755)
            self.assertGreaterEqual(destination.stat().st_mtime, started)
            self.assertEqual(binary.stat().st_mtime, 1)
            self.assertEqual(stat.S_IMODE(binary.stat().st_mode), 0o555)
            self.assertEqual(
                (self.root / "target/nix/bin" / profile / "gammaloop").resolve(),
                self.output,
            )
            self.assert_cleaned(records)

    def test_partial_build_retains_independent_binary_roots(self):
        self.enable_nix()
        self.binary("first", "first")
        self.binary("second", "second")
        result, _ = self.invoke(["build", "--profile=dev"])
        self.assertEqual(result.returncode, 0, result.stderr)
        replacement = self.root / "replacement"
        self.binary("second", "replacement", output=replacement)
        self.env["TEST_OUTPUT"] = str(replacement)
        result, records = self.invoke(["build", "--profile=dev", "--bin", "second"])
        self.assertEqual(result.returncode, 0, result.stderr)
        roots = self.root / "target/nix/bin/dev"
        self.assertEqual((roots / "first").resolve(), self.output)
        self.assertEqual((roots / "second").resolve(), replacement)
        self.assertEqual((self.root / "target/debug/first").read_text(), "first")
        self.assertEqual((self.root / "target/debug/second").read_text(), "replacement")
        self.assert_cleaned(records)

    def test_failed_binary_root_retains_pending_and_previous_closures(self):
        self.enable_nix()
        self.binary()
        previous = self.root / "previous"
        previous.mkdir()
        permanent = self.root / "target/nix/bin/dev/gammaloop"
        permanent.parent.mkdir(parents=True)
        permanent.symlink_to(previous)
        self.env["TEST_ROOT_EXIT"] = "23"
        result, records = self.invoke(["build", "--profile=dev"])
        self.assertEqual(result.returncode, 23, result.stderr)
        pending = list((self.root / "target/nix").glob(".pending-*"))
        self.assertEqual(len(pending), 1)
        self.assertEqual(pending[0].resolve(), self.output)
        self.assertEqual(permanent.resolve(), previous)
        self.assertIn(str(pending[0]), result.stderr)
        self.assertEqual(
            list((self.root / "target/debug").iterdir()),
            [self.root / "target/debug/gammaloop"],
        )
        pending[0].unlink()
        self.assert_cleaned(records)

    def test_failed_binary_copy_preserves_old_binary_and_root(self):
        self.enable_nix()
        self.binary()
        previous = self.root / "previous"
        previous.mkdir()
        permanent = self.root / "target/nix/bin/dev/gammaloop"
        permanent.parent.mkdir(parents=True)
        permanent.symlink_to(previous)
        destination = self.root / "target/debug/gammaloop"
        destination.parent.mkdir(parents=True)
        destination.write_text("old executable")
        # The artifact copy fails after creating a partial temporary file.
        # The real entry point must remove that file.
        (self.tools / "bin/cp").unlink()
        cp = shutil.which("cp")
        self.stub(
            "cp",
            "if sys.argv[1] == '--':\n"
            "    Path(sys.argv[-1]).write_text('partial')\n"
            "    sys.exit(29)\n"
            f"os.execv({cp!r}, [{cp!r}, *sys.argv[1:]])\n",
            directory=self.tools / "bin",
        )
        result, records = self.invoke(["build", "--profile=dev"])
        self.assertEqual(result.returncode, 29, result.stderr)
        self.assertEqual(destination.read_text(), "old executable")
        self.assertEqual(permanent.resolve(), previous)
        self.assertEqual(list(destination.parent.iterdir()), [destination])
        self.assert_cleaned(records)

    def test_run_preserves_literal_args_license_and_pending_until_runtime_finishes(
        self,
    ):
        self.enable_nix()
        (self.output / "bin/nix-test-selected").unlink()
        self.stub(
            "gammaloop",
            "assert os.environ['SYMBOLICA_LICENSE'] == 'runtime-sentinel'\nassert Path(os.environ['GAMMALOOP_FLAKE_SOURCE']).is_dir()\nassert list(Path('target/nix').glob('.pending-*'))\n",
            directory=self.output / "bin",
        )
        self.env.update(SYMBOLICA_LICENSE="runtime-sentinel", TEST_EXIT="19")
        result, records = self.invoke(["run", "--profile=dev", "--", *LITERAL])
        self.assertEqual(result.returncode, 19, result.stderr)
        self.assertEqual(records[-1]["program"], "gammaloop")
        self.assertEqual(records[-1]["args"], LITERAL)
        self.assertTrue(records[-1]["license"])
        self.assert_cleaned(records)

    def test_interruption_cleans_owned_snapshot_bootstrap_and_pending_root(self):
        self.enable_nix()
        ready = self.root / "runtime-ready"
        self.env["TEST_WAIT"] = str(ready)
        process = subprocess.Popen(
            [self.bash, str(self.script), "test", "-p", "kurvst"],
            cwd=self.root,
            env=self.env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        try:
            deadline = time.monotonic() + 10
            while (
                not ready.exists()
                and process.poll() is None
                and time.monotonic() < deadline
            ):
                time.sleep(0.02)
            self.assertTrue(ready.exists(), "runtime never became ready")
            os.killpg(process.pid, signal.SIGTERM)
            _, stderr = process.communicate(timeout=10)
            self.assertIn(
                process.returncode, (-signal.SIGTERM, 128 + signal.SIGTERM), stderr
            )
            self.assert_cleaned(self.records())
        finally:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL)
            process.communicate(timeout=10)

    def test_concurrent_installs_leave_matching_binary_and_closure(self):
        self.enable_nix()
        (self.output / "bin/nix-test-selected").unlink()
        self.binary(content="first")
        replacement = self.root / "replacement"
        self.binary(content="second", output=replacement)
        processes = []
        try:
            for output in (self.output, replacement):
                env = dict(self.env, TEST_OUTPUT=str(output), TEST_ROOT_DELAY="0.1")
                processes.append(
                    subprocess.Popen(
                        [self.bash, str(self.script), "build", "--profile=dev"],
                        cwd=self.root,
                        env=env,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                    )
                )
            for process in processes:
                _, stderr = process.communicate(timeout=20)
                self.assertEqual(process.returncode, 0, stderr)
            root = self.root / "target/nix/bin/dev/gammaloop"
            destination = self.root / "target/debug/gammaloop"
            self.assertEqual(
                destination.read_bytes(), (root / "bin/gammaloop").read_bytes()
            )
            self.assert_cleaned(self.records())
        finally:
            for process in processes:
                if process.poll() is None:
                    process.kill()
                process.communicate(timeout=10)

    def test_launcher_profile_selection_overrides_ties_and_debug_fallback(self):
        launcher = self.root / "gammaloop"
        shutil.copyfile(REPOSITORY / "gammaloop", launcher)
        (self.bin / "tr").symlink_to(shutil.which("tr"))
        binaries = {}
        for offset, profile in enumerate(("release", "dev-optim", "ci-optim", "debug")):
            binary = self.root / "target" / profile / "gammaloop"
            binary.parent.mkdir(parents=True)
            binary.write_text(
                f'#!{self.bash}\nprintf "{profile}\\n"\nprintf "<%s>\\n" "$@"\n'
            )
            binary.chmod(0o755)
            os.utime(binary, (100 + offset, 100 + offset))
            binaries[profile] = binary

        def selected(arguments, expected):
            result = subprocess.run(
                [self.bash, str(launcher), *arguments, LITERAL[0]],
                cwd=self.root,
                env=self.env,
                capture_output=True,
                check=False,
                timeout=10,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, f"{expected}\n<{LITERAL[0]}>\n".encode())

        selected([], "ci-optim")
        for profile in binaries:
            selected([f"--{profile}"], profile)
        os.utime(binaries["release"], (102, 102))
        selected([], "release")
        for profile in ("release", "dev-optim", "ci-optim"):
            binaries[profile].unlink()
        selected([], "debug")

    def test_launcher_selects_fresh_nix_install_over_older_release(self):
        self.enable_nix()
        launcher = self.root / "gammaloop"
        shutil.copyfile(REPOSITORY / "gammaloop", launcher)
        (self.bin / "tr").symlink_to(shutil.which("tr"))
        binary = self.binary(
            content=f'#!{self.bash}\nprintf "nix-dev\\n"\nprintf "<%s>\\n" "$@"\n'
        )
        release = self.root / "target/release/gammaloop"
        release.parent.mkdir(parents=True)
        release.write_text(
            f'#!{self.bash}\nprintf "old-release\\n"\nprintf "<%s>\\n" "$@"\n'
        )
        release.chmod(0o755)
        old = time.time() - 3600
        os.utime(release, (old, old))
        destination = self.root / "target/dev-optim/gammaloop"
        destination.parent.mkdir(parents=True)
        destination.symlink_to(binary)

        def launched():
            return subprocess.run(
                [self.bash, str(launcher), LITERAL[0]],
                cwd=self.root,
                env=self.env,
                capture_output=True,
                check=True,
                timeout=10,
            ).stdout

        self.assertEqual(launched(), f"old-release\n<{LITERAL[0]}>\n".encode())
        result, records = self.invoke(["build", "--profile=dev-optim"])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(launched(), f"nix-dev\n<{LITERAL[0]}>\n".encode())
        self.assertFalse(destination.is_symlink())
        self.assertEqual(binary.stat().st_mtime, 1)
        self.assert_cleaned(records)


if __name__ == "__main__":
    unittest.main()
