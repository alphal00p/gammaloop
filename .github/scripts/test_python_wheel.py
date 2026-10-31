"""Native development command contracts, without builds or license access."""

import base64
import contextlib
import csv
import fcntl
import hashlib
import importlib.util
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from unittest import mock
from urllib.parse import quote

REPOSITORY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY / "scripts"))
MODULE_SPEC = importlib.util.spec_from_file_location(
    "python_wheel", REPOSITORY / "scripts/python_wheel.py"
)
wheel = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(wheel)


class WheelResourcesTests(unittest.TestCase):
    def test_missing_font_provider_does_not_copy_destination_onto_itself(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)
            (destination / "LICENSE").write_text("existing notice")
            with mock.patch.object(wheel, "run", return_value='{"packages": []}'):
                wheel.wheel_render_data({"cargo": "cargo"}, "Cargo.toml", destination)
            self.assertEqual((destination / "LICENSE").read_text(), "existing notice")

    def test_native_repair_updates_record_without_resource_initializer(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "example.whl"
            native = root / "native.so"
            native.write_bytes(b"repaired native payload")
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("gammaloop/native.so", b"original")
                archive.writestr("gammaloop-0.0.0.dist-info/RECORD", "")
            wheel.wheel_data_directory(
                path, None, "gammaloop", {"gammaloop/native.so": native}
            )
            with zipfile.ZipFile(path) as archive:
                self.assertEqual(
                    archive.read("gammaloop/native.so"), native.read_bytes()
                )
                record = archive.read("gammaloop-0.0.0.dist-info/RECORD").decode()
                digest = (
                    base64.urlsafe_b64encode(
                        hashlib.sha256(native.read_bytes()).digest()
                    )
                    .rstrip(b"=")
                    .decode()
                )
                self.assertIn("sha256=" + digest, record)

    def test_packaged_resource_initializer_and_record_remain_valid(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "example.whl"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr(
                    "gammaloop/__init__.py",
                    '"""Package docs."""\nfrom __future__ import annotations\nvalue = 7\n',
                )
                archive.writestr("gammaloop/native.so", b"native payload")
                archive.writestr("gammaloop-0.0.0.dist-info/RECORD", "")
            wheel.wheel_data_directory(
                path, Path("/nix/store/example/data"), "gammaloop._gammaloop"
            )
            with zipfile.ZipFile(path) as archive:
                source = archive.read("gammaloop/__init__.py").decode()
                compile(source, "__init__.py", "exec")
                self.assertIn('setdefault("SYMBOLICA_DATA_DIRECTORY",', source)
                records = csv.reader(
                    io.StringIO(
                        archive.read("gammaloop-0.0.0.dist-info/RECORD").decode()
                    )
                )
                for name, digest, size in records:
                    if not digest:
                        self.assertTrue(name.endswith("/RECORD"))
                        continue
                    content = archive.read(name)
                    self.assertEqual(int(size), len(content))
                    self.assertEqual(
                        digest,
                        "sha256="
                        + base64.urlsafe_b64encode(hashlib.sha256(content).digest())
                        .rstrip(b"=")
                        .decode(),
                    )


class WheelEntryPointTests(unittest.TestCase):
    def test_help_describes_the_packaging_options(self):
        output = io.StringIO()
        with (
            contextlib.redirect_stdout(output),
            self.assertRaises(SystemExit) as raised,
        ):
            wheel.main(["--backend", "nix", "gammaloop", "--help"])
        self.assertEqual(raised.exception.code, 0)
        self.assertIn("--install", output.getvalue())
        self.assertIn("--abi3", output.getvalue())
        self.assertIn("--python", output.getvalue())

    def test_duplicate_backend_cannot_override_the_shell_choice(self):
        for arguments in (
            ["--backend", "nix", "gammaloop", "--backend", "cargo"],
            ["--backend=cargo", "--backend=nix", "community"],
        ):
            with (
                self.subTest(arguments=arguments),
                mock.patch.object(wheel, "python_package") as package,
                contextlib.redirect_stderr(io.StringIO()),
                self.assertRaises(SystemExit) as raised,
            ):
                wheel.main(arguments)
            self.assertEqual(raised.exception.code, 2)
            package.assert_not_called()


class LicenseEnvironmentTests(unittest.TestCase):
    def test_empty_or_missing_config_home_uses_home_config(self):
        with tempfile.TemporaryDirectory() as temporary:
            home = Path(temporary)
            license_file = home / ".config/symbolica/license"
            license_file.parent.mkdir(parents=True)
            license_file.write_text("  home-license-sentinel\n")
            for env in ({}, {"XDG_CONFIG_HOME": ""}):
                with (
                    self.subTest(env=env),
                    mock.patch.dict(os.environ, env, clear=True),
                    mock.patch.object(Path, "home", return_value=home) as resolve_home,
                ):
                    self.assertEqual(
                        wheel.license_environment()["SYMBOLICA_LICENSE"],
                        "home-license-sentinel",
                    )
                    resolve_home.assert_called_once_with()
                    self.assertEqual(dict(os.environ), env)

    def test_explicit_config_home_does_not_resolve_home(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = Path(temporary)
            license_file = config / "symbolica/license"
            license_file.parent.mkdir()
            license_file.write_text("  config-license-sentinel\n")
            with (
                mock.patch.dict(
                    os.environ, {"XDG_CONFIG_HOME": str(config)}, clear=True
                ),
                mock.patch.object(
                    Path, "home", side_effect=AssertionError("home was resolved")
                ),
            ):
                self.assertEqual(
                    wheel.license_environment()["SYMBOLICA_LICENSE"],
                    "config-license-sentinel",
                )

    def test_exported_signed_and_config_license_precedence(self):
        with tempfile.TemporaryDirectory() as temporary:
            config = Path(temporary)
            license_file = config / "symbolica/license"
            license_file.parent.mkdir()
            license_file.write_text("  config-license-sentinel\n")
            for exported, signed, expected in (
                (
                    "exported-license-sentinel",
                    "signed-license-sentinel",
                    "exported-license-sentinel",
                ),
                ("", "signed-license-sentinel", "signed-license-sentinel"),
                (None, "signed-license-sentinel", "signed-license-sentinel"),
                (None, "", "config-license-sentinel"),
                (None, None, "config-license-sentinel"),
            ):
                with self.subTest(exported=exported, signed=signed):
                    env = {"XDG_CONFIG_HOME": str(config), "PRESERVE_THIS": "value"}
                    if exported is not None:
                        env["SYMBOLICA_LICENSE"] = exported
                    if signed is not None:
                        env["SYMBOLICA_LICENSE_SIGNED"] = signed
                    with mock.patch.dict(os.environ, env, clear=True):
                        loaded = wheel.license_environment()
                        self.assertEqual(loaded["SYMBOLICA_LICENSE"], expected)
                        self.assertEqual(loaded["PRESERVE_THIS"], "value")
                        self.assertEqual(dict(os.environ), env)

    def test_exported_and_signed_licenses_never_read_config(self):
        for env in (
            {"SYMBOLICA_LICENSE": "exported-license-sentinel"},
            {"SYMBOLICA_LICENSE_SIGNED": "signed-license-sentinel"},
        ):
            with (
                self.subTest(env=env),
                mock.patch.dict(os.environ, env, clear=True),
                mock.patch.object(
                    Path, "read_text", side_effect=AssertionError("license file read")
                ),
            ):
                self.assertEqual(
                    wheel.license_environment()["SYMBOLICA_LICENSE"],
                    next(iter(env.values())),
                )


class WheelFixture(unittest.TestCase):
    def setUp(self):
        self.root = Path(self.enterContext(tempfile.TemporaryDirectory())).resolve()
        (self.root / ".config").mkdir()
        shutil.copyfile(
            REPOSITORY / ".config/build-profile", self.root / ".config/build-profile"
        )
        self.output = self.root / "nix-output"
        self.enterContext(mock.patch.object(wheel, "ROOT", self.root))
        self.enterContext(mock.patch.dict(os.environ, {}, clear=True))
        self.stdout = self.enterContext(contextlib.redirect_stdout(io.StringIO()))
        self.stderr = self.enterContext(contextlib.redirect_stderr(io.StringIO()))
        self.source = self.root / "filtered-source"
        self.source.mkdir()

        (self.source / "flake.nix").touch()
        os.environ["GAMMALOOP_FLAKE_SOURCE"] = str(self.source)

        # The real helper can read a user's license file. Never call it here.
        self.license_environment = self.enterContext(
            mock.patch.object(
                wheel, "license_environment", return_value={"PATH": "/tools"}
            )
        )
        self.which = self.enterContext(
            mock.patch.object(
                wheel.shutil,
                "which",
                side_effect=lambda name: "/tools/nix" if name == "nix" else None,
            )
        )
        self.process = self.enterContext(
            mock.patch.object(
                wheel.subprocess,
                "run",
                return_value=subprocess.CompletedProcess([], 0, f"{self.output}\n"),
            )
        )
        self.process.side_effect = self.nix_process

    def nix_process(self, command, **kwargs):
        if command[:2] == ["nix", "build"] and "--out-link" in command:
            root = Path(command[command.index("--out-link") + 1])
            root.unlink(missing_ok=True)
            root.symlink_to(self.output)
        elif command[0] == "nix-store":
            root = Path(command[command.index("--add-root") + 1])
            output = Path(command[command.index("--realise") + 1])
            root.unlink(missing_ok=True)
            root.symlink_to(output)
        return self.process.return_value

    def nix_request(self, *, call_count=1):
        self.assertEqual(self.process.call_count, call_count)
        invocation = self.process.call_args_list[0]
        command = invocation.args[0]
        self.assertIsInstance(command, list)
        self.assertFalse(invocation.kwargs.get("shell", False))
        self.assertEqual(
            command[:5],
            ["nix", "build", "--impure", "--print-out-paths", "--no-write-lock-file"],
        )
        return command, invocation.kwargs["env"]


class SnapshotTests(WheelFixture):
    def test_source_path_is_encoded_and_borrowed(self):
        source = self.root / "source #with ?literal %characters"
        source.mkdir()
        (source / "flake.nix").touch()
        with mock.patch.dict(os.environ, {"GAMMALOOP_FLAKE_SOURCE": str(source)}):
            wheel.nix_build("development-wheel")
        command, _ = self.nix_request()
        self.assertEqual(
            command[-1], f"path:{quote(str(source), safe='/')}#development-wheel"
        )
        self.assertTrue((source / "flake.nix").is_file())

    def test_backend_is_required_and_rust_actions_are_rejected(self):
        for arguments in (
            ["gammaloop"],
            ["--backend", "auto", "gammaloop"],
            ["--backend", "cargo", "test"],
        ):
            with self.subTest(arguments=arguments), self.assertRaises(SystemExit):
                wheel.main(arguments)
            self.process.assert_not_called()

    def test_missing_or_invalid_source_never_builds(self):
        for source in ("", str(self.root / "missing"), str(self.root)):
            with (
                self.subTest(source=source),
                mock.patch.dict(os.environ, {"GAMMALOOP_FLAKE_SOURCE": source}),
            ):
                with self.assertRaisesRegex(ValueError, "GAMMALOOP_FLAKE_SOURCE"):
                    wheel.nix_build("development-wheel")
                self.process.assert_not_called()

    def test_nix_build_scrubs_runtime_secrets_and_stale_requests(self):
        with mock.patch.dict(
            os.environ,
            {
                "SYMBOLICA_LICENSE": "sentinel",
                "SYMBOLICA_LICENSE_SIGNED": "sentinel",
                "GAMMALOOP_NIX_TEST_ARGS": "stale",
                "GAMMALOOP_DEVELOP_SPEC": "stale",
            },
        ):
            wheel.nix_build("development-wheel", spec={"package": "wheel"})
            command, env = self.nix_request()
            for key in (
                "SYMBOLICA_LICENSE",
                "SYMBOLICA_LICENSE_SIGNED",
                "GAMMALOOP_NIX_TEST_ARGS",
            ):
                self.assertNotIn(key, env)
            self.assertEqual(
                json.loads(env["GAMMALOOP_DEVELOP_SPEC"]), {"package": "wheel"}
            )
            self.assertEqual(command[-1], f"path:{self.source}#development-wheel")
            self.assertTrue((self.source / "flake.nix").is_file())


class PythonEnvironmentTests(WheelFixture):
    def setUp(self):
        super().setUp()
        self.details = {
            "version": [3, 13],
            "implementation": "cpython",
            "free_threaded": False,
            "prefix": str(self.root / "selected-venv"),
            "base": "/system/python",
        }
        self.process.return_value.stdout = json.dumps(self.details)

    def test_explicit_virtualenv_conda_and_local_interpreters_have_defined_precedence(
        self,
    ):
        local = self.root / ".venv/bin/python"
        local.parent.mkdir(parents=True)
        local.touch()
        for requested, environment, expected in (
            (
                "/chosen/python",
                {"VIRTUAL_ENV": "/active", "CONDA_PREFIX": "/conda"},
                "/chosen/python",
            ),
            (
                None,
                {"VIRTUAL_ENV": "/active", "CONDA_PREFIX": "/conda"},
                "/active/bin/python",
            ),
            (None, {"CONDA_PREFIX": "/conda"}, "/conda/bin/python"),
            (None, {}, str(local)),
        ):
            with self.subTest(requested=requested, environment=environment):
                with mock.patch.dict(os.environ, environment, clear=True):
                    executable, details = wheel.python_environment(
                        requested, install=True
                    )
                self.assertEqual(executable, expected)
                self.assertEqual(details, self.details)
                self.assertEqual(self.process.call_args.args[0][:2], [expected, "-c"])
                self.assertFalse(self.process.call_args.kwargs.get("shell", False))

    def test_named_interpreter_is_resolved_on_path(self):
        self.which.side_effect = lambda name: (
            "/custom/python3.12" if name == "python3.12" else None
        )
        executable, _ = wheel.python_environment("python3.12")
        self.assertEqual(executable, "/custom/python3.12")

    def test_without_environment_wheel_build_uses_current_python(self):
        executable, _ = wheel.python_environment()
        self.assertEqual(executable, wheel.sys.executable)

    def test_system_interpreter_is_allowed_for_build_but_not_install(self):
        self.details["prefix"] = self.details["base"]
        self.process.return_value.stdout = json.dumps(self.details)
        wheel.python_environment("/system/python", install=False)
        with self.assertRaisesRegex(ValueError, "virtual environment"):
            wheel.python_environment("/system/python", install=True)

    def test_conda_environment_can_install_with_equal_base_prefix(self):
        self.details["prefix"] = self.details["base"] = "/conda/env"
        self.process.return_value.stdout = json.dumps(self.details)
        with mock.patch.dict(os.environ, {"CONDA_PREFIX": "/conda/env"}):
            executable, _ = wheel.python_environment(install=True)
        self.assertEqual(executable, "/conda/env/bin/python")

    def test_unrelated_conda_environment_does_not_authorize_system_install(self):
        self.details["prefix"] = self.details["base"]
        self.process.return_value.stdout = json.dumps(self.details)
        with (
            mock.patch.dict(os.environ, {"CONDA_PREFIX": "/conda/env"}),
            self.assertRaisesRegex(ValueError, "virtual environment|Conda|conda"),
        ):
            wheel.python_environment("/system/python", install=True)


class PythonPackageTests(WheelFixture):
    def setUp(self):
        super().setUp()
        self.venv = self.root / "selected env; $(never) `never`"
        self.venv.mkdir()
        self.interpreter = str(self.venv / "bin/python")
        self.details = {
            "version": [3, 12],
            "implementation": "cpython",
            "free_threaded": False,
            "prefix": str(self.venv),
            "base": "/system/python",
        }
        self.environment = self.enterContext(
            mock.patch.object(
                wheel,
                "python_environment",
                return_value=(self.interpreter, self.details),
            )
        )
        self.wheel = self.output / "dist/gammaloop-0.0.0-cp312-cp312-linux_x86_64.whl"
        self.wheel.parent.mkdir(parents=True)
        self.wheel.write_bytes(b"fixture wheel")
        self.native_payload()
        self.process.side_effect = self.nix_process

    def native_payload(self):
        for name in ("gammaloop/_native.so", "runtime.dylib", "README.txt"):
            path = self.output / "lib" / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"fixture native payload; never actually loaded")

    def assert_locked(self, path):
        with (
            path.open("a") as lock,
            self.assertRaises(BlockingIOError),
        ):
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def nix_process(self, command, **kwargs):
        self.assertIsInstance(command, list)
        self.assertFalse(kwargs.get("shell", False))
        if command[:2] == ["nix", "build"]:
            root = Path(command[command.index("--out-link") + 1])
            root.unlink(missing_ok=True)
            root.symlink_to(self.output)
        elif command[:3] == ["uv", "pip", "install"]:
            self.assert_locked(self.venv / ".gammaloop-install.lock")
            self.assertEqual(Path(command[-1]), self.wheel)
            self.assertTrue(self.wheel.is_file())
        elif command[:2] == [self.interpreter, "-c"]:
            self.assert_locked(self.venv / ".gammaloop-install.lock")
            if "ctypes.CDLL" in command[2]:
                self.assertEqual(
                    command[3:],
                    [
                        str(self.output / "lib/gammaloop/_native.so"),
                        str(self.output / "lib/runtime.dylib"),
                    ],
                )
                self.assertTrue(all(Path(path).is_file() for path in command[3:]))
            else:
                self.assertIn(
                    command[2],
                    (
                        "from gammaloop import GammaLoopAPI, Expression, Model, FeynmanDiagram, TensorExpression",
                        "import symbolica.community.hepkit, symbolica.community.graph, symbolica.community.tensor, symbolica.community.hepkit.gammaloop",
                    ),
                )
                self.assertEqual(len(command), 3)
        elif command[0] == "nix-store":
            root = Path(command[command.index("--add-root") + 1])
            output = Path(command[command.index("--realise") + 1])
            if root.parent == self.venv / ".nix-python":
                self.assert_locked(self.venv / ".gammaloop-install.lock")
            else:
                self.assert_locked(root.parent / ".gammaloop-wheel.lock")
            root.unlink(missing_ok=True)
            root.symlink_to(output)
        else:
            self.fail(f"Unexpected process: {command}")
        return subprocess.CompletedProcess(command, 0, f"{self.output}\n")

    def previous_root(self, name="gammaloop"):
        previous = self.root / f"previous-{name}"
        previous.mkdir()
        root = self.venv / ".nix-python" / name
        root.parent.mkdir(parents=True, exist_ok=True)
        root.symlink_to(previous)
        return root, previous

    def test_nix_python_variants_profiles_and_interpreter(self):
        for (
            distribution,
            extra,
            package,
            profile,
            features,
            defaults,
            manifest,
            pyproject,
        ) in (
            (
                "gammaloop",
                [],
                "gammaloop-api",
                "ci-optim",
                ["pyo3-extension-module"],
                False,
                "crates/gammaloop-api/Cargo.toml",
                ".",
            ),
            (
                "gammaloop",
                ["--abi3"],
                "gammaloop-api",
                "ci-optim",
                ["pyo3-extension-module", "python_abi"],
                False,
                "crates/gammaloop-api/Cargo.toml",
                ".",
            ),
            (
                "community",
                [],
                "alphal00p-notebook-host",
                "ci-optim",
                [],
                True,
                "examples/notebooks/symbolica-host/Cargo.toml",
                "examples/notebooks/symbolica-host",
            ),
        ):
            with self.subTest(distribution=distribution, extra=extra):
                self.process.reset_mock()
                wheel.main(
                    [
                        "--backend",
                        "nix",
                        distribution,
                        "--python",
                        self.interpreter,
                        *extra,
                    ]
                )
                command, env = self.nix_request(call_count=2)
                self.assertEqual(command[-1], f"path:{self.source}#development-wheel")
                self.assertEqual(
                    json.loads(env["GAMMALOOP_DEVELOP_SPEC"]),
                    {
                        "package": package,
                        "manifest": manifest,
                        "pyproject": pyproject,
                        "profile": profile,
                        "features": features,
                        "default_features": defaults,
                        "python": "python312",
                        "args": [],
                    },
                )
                self.environment.assert_called_with(self.interpreter, False)
                name = "symbolica" if distribution == "community" else "gammaloop"
                output_dir = self.root / "target/wheels"
                pending = Path(command[6])
                self.assertEqual(pending.parent, output_dir)
                self.assertRegex(pending.name, rf"^\.pending-{name}-[0-9a-f]{{32}}$")
                self.assertFalse(pending.is_symlink())
                export_root = output_dir / f".nix-{self.wheel.name}"
                self.assertEqual(
                    self.process.call_args_list[1].args[0],
                    [
                        "nix-store",
                        "--realise",
                        str(self.output),
                        "--add-root",
                        str(export_root),
                        "--indirect",
                    ],
                )
                self.assertTrue(export_root.is_symlink())
                self.assertEqual(export_root.resolve(), self.output)
                self.assertEqual(
                    (self.root / "target/wheels" / self.wheel.name).read_bytes(),
                    b"fixture wheel",
                )

    def test_host_rejects_unsupported_default_profile_before_build(self):
        with self.assertRaisesRegex(ValueError, "no dev-optim profile"):
            wheel.main(["--backend", "nix", "community", "--profile", "dev-optim"])
        self.process.assert_not_called()

    def test_non_native_distributions_are_rejected_before_build(self):
        for distribution in ("pyodide", "browser", "tydenso"):
            with self.subTest(distribution=distribution), self.assertRaises(SystemExit):
                wheel.main(["--backend", "nix", distribution])
        self.process.assert_not_called()

    def test_nix_rejects_pypy_free_threaded_and_debug_python(self):
        for implementation, free_threaded, debug in (
            ("pypy", False, False),
            ("cpython", True, False),
            ("cpython", False, True),
        ):
            with self.subTest(
                implementation=implementation, free_threaded=free_threaded, debug=debug
            ):
                self.details.update(
                    implementation=implementation,
                    free_threaded=free_threaded,
                    debug=debug,
                )
                with self.assertRaisesRegex(ValueError, "standard CPython"):
                    wheel.main(["--backend", "nix", "gammaloop"])
                self.process.assert_not_called()

    def test_nix_rejects_mismatched_python_architecture_including_abi3(self):
        self.details["machine"] = "x86_64"
        for extra in ([], ["--abi3"]):
            with (
                self.subTest(extra=extra),
                mock.patch.object(wheel.platform, "machine", return_value="aarch64"),
                self.assertRaisesRegex(ValueError, "architecture.*Nix host"),
            ):
                wheel.main(["--backend", "nix", "gammaloop", "--install", *extra])
            self.process.assert_not_called()

    def test_nix_accepts_matching_python_architecture_aliases(self):
        for host, interpreter in (
            ("aarch64", "arm64"),
            ("arm64", "aarch64"),
            ("x86_64", "AMD64"),
            ("AMD64", "x86_64"),
        ):
            with (
                self.subTest(host=host, interpreter=interpreter),
                mock.patch.object(wheel.platform, "machine", return_value=host),
            ):
                self.process.reset_mock()
                self.details.update(machine=interpreter, debug=False)
                wheel.main(["--backend", "nix", "gammaloop"])
                self.nix_request(call_count=2)

    def test_failed_nix_wheel_never_runs_maturin_or_installer(self):
        root, previous = self.previous_root()
        for rooted in (False, True):
            with self.subTest(rooted=rooted):
                self.process.reset_mock()

                def completed(command, rooted=rooted, **kwargs):
                    if rooted:
                        self.nix_process(command, **kwargs)
                    return subprocess.CompletedProcess(command, 43, "")

                self.process.side_effect = completed
                with self.assertRaises(subprocess.CalledProcessError) as raised:
                    wheel.main(["--backend", "nix", "gammaloop", "--install"])
                self.assertEqual(raised.exception.returncode, 43)
                self.process.assert_called_once()
                self.assertEqual(self.process.call_args.args[0][0], "nix")
                self.assertEqual(root.resolve(), previous)
                pending = Path(self.process.call_args.args[0][6])
                self.assertFalse(pending.is_symlink())
                self.assertFalse(pending.exists())

    def test_failed_or_interrupted_build_removes_only_its_pending_root(self):
        environment_root, previous = self.previous_root()
        for install in (False, True):
            directory = (
                environment_root.parent if install else self.root / "target/wheels"
            )
            directory.mkdir(parents=True, exist_ok=True)
            retained = directory / ".pending-gammaloop-earlier-install"
            retained.symlink_to(previous)
            for failure in (
                OSError("Nix failed"),
                KeyboardInterrupt(),
                SystemExit(130),
            ):
                with self.subTest(install=install, failure=type(failure).__name__):
                    self.process.reset_mock()

                    def interrupted(command, failure=failure, **kwargs):
                        self.nix_process(command, **kwargs)
                        raise failure

                    self.process.side_effect = interrupted
                    with self.assertRaises(type(failure)):
                        wheel.main(
                            [
                                "--backend",
                                "nix",
                                "gammaloop",
                                *(["--install"] if install else []),
                            ]
                        )
                    self.process.assert_called_once()
                    pending = Path(self.process.call_args.args[0][6])
                    self.assertEqual(pending.parent, directory)
                    self.assertFalse(pending.is_symlink())
                    self.assertFalse(pending.exists())
                    self.assertEqual(environment_root.resolve(), previous)
                    self.assertEqual(retained.resolve(), previous)
                    self.assertEqual(list(directory.glob(".pending-*")), [retained])

    def test_install_selects_environment_package_and_persistent_gc_root(self):
        for distribution, name, no_deps in (
            ("gammaloop", "gammaloop", []),
            ("community", "symbolica", ["--no-deps"]),
        ):
            with self.subTest(distribution=distribution):
                self.process.reset_mock()
                output_dir = self.root / "wheel output; $(never)"
                root = self.venv / ".nix-python" / name
                wheel.main(
                    [
                        "--backend",
                        "nix",
                        distribution,
                        "--install",
                        "--out",
                        str(output_dir),
                        *no_deps,
                    ]
                )
                build, export, loader, install, imported, commit = (
                    self.process.call_args_list
                )
                pending = Path(build.args[0][6])
                self.assertEqual(pending.parent, root.parent)
                self.assertRegex(pending.name, rf"^\.pending-{name}-[0-9a-f]{{32}}$")
                self.assertEqual(
                    build.args[0],
                    [
                        "nix",
                        "build",
                        "--impure",
                        "--print-out-paths",
                        "--no-write-lock-file",
                        "--out-link",
                        str(pending),
                        "--accept-flake-config",
                        f"path:{self.source}#development-wheel",
                    ],
                )
                self.assertEqual(
                    install.args[0],
                    [
                        "uv",
                        "pip",
                        "install",
                        "--python",
                        self.interpreter,
                        "--reinstall-package",
                        name,
                        *no_deps,
                        str(self.wheel),
                    ],
                )
                self.environment.assert_called_with(None, True)
                self.assertFalse(install.kwargs.get("shell", False))
                self.assertEqual(loader.args[0][:2], [self.interpreter, "-c"])
                self.assertIn("ctypes.CDLL", loader.args[0][2])
                self.assertEqual(
                    imported.args[0],
                    [
                        self.interpreter,
                        "-c",
                        "import symbolica.community.hepkit, symbolica.community.graph, symbolica.community.tensor, symbolica.community.hepkit.gammaloop"
                        if distribution == "community"
                        else "from gammaloop import GammaLoopAPI, Expression, Model, FeynmanDiagram, TensorExpression",
                    ],
                )
                export_root = output_dir / f".nix-{self.wheel.name}"
                self.assertEqual(
                    export.args[0],
                    [
                        "nix-store",
                        "--realise",
                        str(self.output),
                        "--add-root",
                        str(export_root),
                        "--indirect",
                    ],
                )
                self.assertTrue(export_root.is_symlink())
                self.assertEqual(export_root.resolve(), self.output)
                self.assertEqual(
                    commit.args[0],
                    [
                        "nix-store",
                        "--realise",
                        str(self.output),
                        "--add-root",
                        str(root),
                        "--indirect",
                    ],
                )
                self.assertTrue(root.is_symlink())
                self.assertEqual(root.resolve(), self.output)
                self.assertFalse(pending.is_symlink())

    def test_repeated_export_replaces_read_only_wheel_with_writable_copy(self):
        output_dir = self.root / "repeated wheels"
        output_dir.mkdir()
        wheel_name = self.wheel.name
        exported = output_dir / wheel_name
        exported.write_bytes(b"previous read-only export")
        export_root = output_dir / f".nix-{wheel_name}"
        environment_root = self.venv / ".nix-python/gammaloop"
        temporary_paths = []
        copyfile = shutil.copyfile

        def copy(source, temporary):
            self.assert_locked(output_dir / ".gammaloop-wheel.lock")
            self.assertEqual(source, self.wheel)
            self.assertEqual(temporary.parent, output_dir)
            self.assertRegex(temporary.name, r"^\.[0-9a-f]{32}\.whl$")
            self.assertFalse(temporary.exists())
            temporary_paths.append(temporary)
            return copyfile(source, temporary)

        with mock.patch.object(wheel.shutil, "copyfile", side_effect=copy):
            for index, contents in enumerate((b"first rebuild", b"second rebuild")):
                with self.subTest(rebuild=index):
                    self.process.reset_mock()
                    self.output = self.root / f"rebuilt-output-{index}"
                    self.wheel = self.output / "dist" / wheel_name
                    self.wheel.parent.mkdir(parents=True)
                    self.wheel.write_bytes(contents)
                    self.native_payload()
                    self.wheel.chmod(0o444)
                    exported.chmod(0o444)
                    previous = exported.read_bytes()
                    with exported.open("rb") as reader:
                        wheel.main(
                            [
                                "--backend",
                                "nix",
                                "gammaloop",
                                "--install",
                                "--out",
                                str(output_dir),
                            ]
                        )
                        # An existing reader keeps the old inode across replacement.
                        self.assertEqual(reader.read(), previous)
                    self.assertEqual(exported.read_bytes(), contents)
                    self.assertTrue(exported.stat().st_mode & 0o200)
                    self.assertFalse(self.wheel.stat().st_mode & 0o200)
                    self.assertEqual(export_root.resolve(), self.output)
                    self.assertEqual(environment_root.resolve(), self.output)
                    self.assertEqual(
                        [call.args[0][0] for call in self.process.call_args_list],
                        [
                            "nix",
                            "nix-store",
                            self.interpreter,
                            "uv",
                            self.interpreter,
                            "nix-store",
                        ],
                    )
                    self.assertEqual(
                        set(output_dir.iterdir()),
                        {exported, export_root, output_dir / ".gammaloop-wheel.lock"},
                    )
        self.assertEqual(len(temporary_paths), 2)
        self.assertNotEqual(*temporary_paths)
        self.assertEqual(list(environment_root.parent.iterdir()), [environment_root])

    def test_different_wheel_tags_keep_independent_export_roots(self):
        output_dir = self.root / "tagged wheels"
        roots = {}
        pending_roots = []
        for minor in (12, 13):
            with self.subTest(python=minor):
                self.process.reset_mock()
                self.details["version"] = [3, minor]
                self.output = self.root / f"python3{minor}-output"
                name = f"gammaloop-0.0.0-cp3{minor}-cp3{minor}-linux_x86_64.whl"
                self.wheel = self.output / "dist" / name
                self.wheel.parent.mkdir(parents=True)
                self.wheel.write_bytes(name.encode())
                self.wheel.chmod(0o444)
                wheel.main(
                    [
                        "--backend",
                        "nix",
                        "gammaloop",
                        "--out",
                        str(output_dir),
                    ]
                )
                command, _ = self.nix_request(call_count=2)
                pending = Path(command[6])
                self.assertEqual(pending.parent, output_dir)
                pending_roots.append(pending)
                roots[output_dir / f".nix-{name}"] = self.output
                self.assertEqual((output_dir / name).read_bytes(), name.encode())
                for root, output in roots.items():
                    self.assertTrue(root.is_symlink())
                    self.assertEqual(root.resolve(), output)
                self.assertEqual(list(output_dir.glob(".pending-*")), [])
        self.assertNotEqual(*pending_roots)
        self.assertFalse((self.venv / ".nix-python").exists())

    def test_replaced_export_cannot_change_install_wheel_or_environment_root(self):
        output_dir = self.root / "target/wheels"
        exported = output_dir / self.wheel.name
        export_root = output_dir / f".nix-{self.wheel.name}"
        competing_output = self.root / "competing-output"
        competing_wheel = competing_output / "dist" / self.wheel.name
        competing_wheel.parent.mkdir(parents=True)
        competing_wheel.write_bytes(b"competing wheel")

        def completed(command, **kwargs):
            if command[0] == "uv":
                # Another export can finish between our export and installation.
                with (output_dir / ".gammaloop-wheel.lock").open("a") as lock:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    exported.write_bytes(competing_wheel.read_bytes())
                    export_root.unlink()
                    export_root.symlink_to(competing_output)
            return self.nix_process(command, **kwargs)

        self.process.side_effect = completed
        wheel.main(["--backend", "nix", "gammaloop", "--install"])
        self.assertEqual(self.process.call_args_list[3].args[0][-1], str(self.wheel))
        self.assertEqual(self.wheel.read_bytes(), b"fixture wheel")
        self.assertEqual(exported.read_bytes(), b"competing wheel")
        self.assertEqual(export_root.resolve(), competing_output)
        environment_root = self.venv / ".nix-python/gammaloop"
        self.assertEqual(environment_root.resolve(), self.output)
        self.assertEqual(list(environment_root.parent.iterdir()), [environment_root])

    def test_failed_export_root_commit_retains_pending_and_previous_roots(self):
        environment_root, previous = self.previous_root()
        output_dir = self.root / "target/wheels"
        output_dir.mkdir(parents=True)
        export_root = output_dir / f".nix-{self.wheel.name}"
        export_root.symlink_to(previous)

        def completed(command, **kwargs):
            if command[0] == "nix-store":
                self.assertEqual(command[4], str(export_root))
                self.assert_locked(output_dir / ".gammaloop-wheel.lock")
                return subprocess.CompletedProcess(command, 31, "")
            return self.nix_process(command, **kwargs)

        self.process.side_effect = completed
        for extra in ([], ["--install"]):
            with self.subTest(extra=extra):
                self.process.reset_mock()
                with self.assertRaises(subprocess.CalledProcessError) as raised:
                    wheel.main(["--backend", "nix", "gammaloop", *extra])
                self.assertEqual(raised.exception.returncode, 31)
                self.assertEqual(
                    [call.args[0][0] for call in self.process.call_args_list],
                    ["nix", "nix-store"],
                )
                pending = Path(self.process.call_args_list[0].args[0][6])
                self.assertTrue(pending.is_symlink())
                self.assertEqual(pending.resolve(), self.output)
                self.assertEqual(export_root.resolve(), previous)
                self.assertEqual(environment_root.resolve(), previous)

    def test_invalid_nix_wheel_count_does_not_install(self):
        root, previous = self.previous_root()
        self.wheel.unlink()
        for count in (0, 2):
            with self.subTest(count=count):
                self.process.reset_mock()
                for index in range(count):
                    (self.output / "dist" / f"wheel{index}.whl").touch()
                with self.assertRaisesRegex(ValueError, "exactly one wheel"):
                    wheel.main(["--backend", "nix", "gammaloop", "--install"])
                self.assertEqual(self.process.call_count, 1)
                self.assertEqual(root.resolve(), previous)
                pending = Path(self.process.call_args.args[0][6])
                self.assertTrue(pending.is_symlink())
                self.assertEqual(pending.resolve(), self.output)

    def test_maturin_builds_exact_variant_and_installs_only_new_wheel(self):
        self.which.side_effect = None
        self.which.return_value = None
        output_dir = self.root / "fallback wheels"
        output_dir.mkdir()
        (output_dir / "unrelated-old.whl").write_bytes(b"old wheel")
        produced_wheel = output_dir / self.wheel.name

        def completed(command, **kwargs):
            if command[0] == "maturin":
                produced_wheel.write_bytes(b"new wheel")
            elif command[0] == "uv":
                self.assert_locked(self.venv / ".gammaloop-install.lock")
            return subprocess.CompletedProcess(command, 0, "")

        self.process.side_effect = completed
        wheel.main(
            [
                "--backend",
                "cargo",
                "gammaloop",
                "--abi3",
                "--profile",
                "release",
                "--install",
                "--no-deps",
                "--python",
                self.interpreter,
                "--out",
                str(output_dir),
            ]
        )
        build, install, imported = self.process.call_args_list
        self.assertEqual(
            imported.args[0],
            [
                self.interpreter,
                "-c",
                "from gammaloop import GammaLoopAPI, Expression, Model, FeynmanDiagram, TensorExpression",
            ],
        )
        self.assertEqual(
            build.args[0],
            [
                "maturin",
                "build",
                "--locked",
                "--manifest-path",
                "crates/gammaloop-api/Cargo.toml",
                "--profile",
                "release",
                "--interpreter",
                self.interpreter,
                "--out",
                str(output_dir),
                "--no-default-features",
                "--features",
                "pyo3-extension-module,python_abi",
            ],
        )
        self.assertEqual(
            install.args[0],
            [
                "uv",
                "pip",
                "install",
                "--python",
                self.interpreter,
                "--reinstall-package",
                "gammaloop",
                "--no-deps",
                str(produced_wheel),
            ],
        )
        self.assertFalse((self.venv / ".nix-python").exists())

    def test_profile_overrides_apply_to_both_distributions(self):
        (self.root / ".config/build-profile").write_text("dev\n")
        for distribution in ("gammaloop", "community"):
            for configured, arguments, expected in (
                (None, [], "dev"),
                ("dev", [], "dev"),
                ("release", [], "release"),
                ("release", ["--profile", "dev"], "dev"),
                ("dev", ["--release"], "release"),
                ("release", ["--profile", "dev", "--release"], "release"),
            ):
                with self.subTest(
                    distribution=distribution,
                    configured=configured,
                    arguments=arguments,
                ):
                    self.process.reset_mock()
                    if configured is None:
                        os.environ.pop("GAMMALOOP_BUILD_PROFILE", None)
                    else:
                        os.environ["GAMMALOOP_BUILD_PROFILE"] = configured
                    wheel.main(["--backend", "nix", distribution, *arguments])
                    _, env = self.nix_request(call_count=2)
                    self.assertEqual(
                        json.loads(env["GAMMALOOP_DEVELOP_SPEC"])["profile"], expected
                    )

    def test_maturin_host_uses_the_shared_profile_and_its_own_manifest(self):
        output_dir = self.root / "host wheels"

        def completed(command, **kwargs):
            (output_dir / "symbolica-fixture.whl").write_bytes(b"host wheel")
            return subprocess.CompletedProcess(command, 0, "")

        self.process.side_effect = completed
        wheel.main(["--backend", "cargo", "community", "--out", str(output_dir)])
        self.process.assert_called_once()
        self.assertEqual(
            self.process.call_args.args[0],
            [
                "maturin",
                "build",
                "--locked",
                "--manifest-path",
                "examples/notebooks/symbolica-host/Cargo.toml",
                "--profile",
                "ci-optim",
                "--interpreter",
                self.interpreter,
                "--out",
                str(output_dir),
            ],
        )

    def test_maturin_success_without_new_wheel_cannot_install_stale_wheel(self):
        self.process.side_effect = None
        output_dir = self.root / "target/wheels"
        output_dir.mkdir(parents=True)
        (output_dir / "stale.whl").touch()
        with self.assertRaisesRegex(ValueError, "exactly one new wheel"):
            wheel.main(["--backend", "cargo", "gammaloop", "--install"])
        self.process.assert_called_once()

    def test_installer_failure_preserves_previous_and_pending_roots_without_rebuild(
        self,
    ):
        root, previous = self.previous_root()

        def completed(command, **kwargs):
            self.assertEqual(root.resolve(), previous)
            result = self.nix_process(command, **kwargs)
            if command[0] == "uv":
                result.returncode = 19
            return result

        self.process.side_effect = completed
        with self.assertRaises(subprocess.CalledProcessError) as raised:
            wheel.main(["--backend", "nix", "gammaloop", "--install"])
        self.assertEqual(raised.exception.returncode, 19)
        self.assertEqual(
            [call.args[0][0] for call in self.process.call_args_list],
            ["nix", "nix-store", self.interpreter, "uv"],
        )
        pending = Path(self.process.call_args_list[0].args[0][6])
        self.assertTrue(root.is_symlink())
        self.assertEqual(root.resolve(), previous)
        self.assertTrue(pending.is_symlink())
        self.assertEqual(pending.resolve(), self.output)
        export_root = self.root / "target/wheels" / f".nix-{self.wheel.name}"
        self.assertEqual(export_root.resolve(), self.output)

    def test_loader_failure_preserves_packages_and_previous_and_pending_roots(self):
        for distribution, name in (
            ("gammaloop", "gammaloop"),
            ("community", "symbolica"),
        ):
            with self.subTest(distribution=distribution):
                self.process.reset_mock()
                root, previous = self.previous_root(name)

                def completed(command, root=root, previous=previous, **kwargs):
                    self.assertEqual(root.resolve(), previous)
                    result = self.nix_process(command, **kwargs)
                    if command[0] == self.interpreter:
                        self.assertIn("ctypes.CDLL", command[2])
                        result.returncode = 17
                    return result

                self.process.side_effect = completed
                with self.assertRaises(subprocess.CalledProcessError) as raised:
                    wheel.main(["--backend", "nix", distribution, "--install"])
                self.assertEqual(raised.exception.returncode, 17)
                self.assertEqual(
                    [call.args[0][0] for call in self.process.call_args_list],
                    ["nix", "nix-store", self.interpreter],
                )
                pending = Path(self.process.call_args_list[0].args[0][6])
                self.assertTrue(root.is_symlink())
                self.assertEqual(root.resolve(), previous)
                self.assertTrue(pending.is_symlink())
                self.assertEqual(pending.resolve(), self.output)
                export_root = self.root / "target/wheels" / f".nix-{self.wheel.name}"
                self.assertEqual(export_root.resolve(), self.output)

    def test_missing_native_payload_rejects_install_before_replacing_packages(self):
        root, previous = self.previous_root()
        (self.output / "lib/gammaloop/_native.so").unlink()
        (self.output / "lib/runtime.dylib").unlink()
        for missing_directory in (False, True):
            with self.subTest(missing_directory=missing_directory):
                self.process.reset_mock()
                if missing_directory:
                    shutil.rmtree(self.output / "lib")
                with self.assertRaisesRegex(ValueError, "no retained native payload"):
                    wheel.main(["--backend", "nix", "gammaloop", "--install"])
                self.assertEqual(
                    [call.args[0][0] for call in self.process.call_args_list],
                    ["nix", "nix-store"],
                )
                pending = Path(self.process.call_args_list[0].args[0][6])
                self.assertTrue(pending.is_symlink())
                self.assertEqual(pending.resolve(), self.output)
                self.assertEqual(root.resolve(), previous)

    def test_import_failure_keeps_pending_root_without_committing_environment_root(
        self,
    ):
        for distribution, name in (
            ("gammaloop", "gammaloop"),
            ("community", "symbolica"),
        ):
            with self.subTest(distribution=distribution):
                self.process.reset_mock()
                root, previous = self.previous_root(name)

                def completed(command, root=root, previous=previous, **kwargs):
                    self.assertEqual(root.resolve(), previous)
                    result = self.nix_process(command, **kwargs)
                    if (
                        command[0] == self.interpreter
                        and "ctypes.CDLL" not in command[2]
                    ):
                        result.returncode = 29
                    return result

                self.process.side_effect = completed
                with self.assertRaises(subprocess.CalledProcessError) as raised:
                    wheel.main(["--backend", "nix", distribution, "--install"])
                self.assertEqual(raised.exception.returncode, 29)
                self.assertEqual(
                    [call.args[0][0] for call in self.process.call_args_list],
                    ["nix", "nix-store", self.interpreter, "uv", self.interpreter],
                )
                pending = Path(self.process.call_args_list[0].args[0][6])
                self.assertTrue(root.is_symlink())
                self.assertEqual(root.resolve(), previous)
                self.assertTrue(pending.is_symlink())
                self.assertEqual(pending.resolve(), self.output)
                export_root = self.root / "target/wheels" / f".nix-{self.wheel.name}"
                self.assertEqual(export_root.resolve(), self.output)

    def test_failed_root_commit_preserves_previous_and_pending_roots(self):
        root, previous = self.previous_root()

        def completed(command, **kwargs):
            self.assertEqual(root.resolve(), previous)
            if command[0] == "nix-store" and command[4] == str(root):
                self.assert_locked(self.venv / ".gammaloop-install.lock")
                return subprocess.CompletedProcess(command, 23, "")
            return self.nix_process(command, **kwargs)

        self.process.side_effect = completed
        with self.assertRaises(subprocess.CalledProcessError) as raised:
            wheel.main(["--backend", "nix", "gammaloop", "--install"])
        self.assertEqual(raised.exception.returncode, 23)
        self.assertEqual(
            [call.args[0][0] for call in self.process.call_args_list],
            ["nix", "nix-store", self.interpreter, "uv", self.interpreter, "nix-store"],
        )
        pending = Path(self.process.call_args_list[0].args[0][6])
        self.assertTrue(root.is_symlink())
        self.assertEqual(root.resolve(), previous)
        self.assertTrue(pending.is_symlink())
        self.assertEqual(pending.resolve(), self.output)
        export_root = self.root / "target/wheels" / f".nix-{self.wheel.name}"
        self.assertEqual(export_root.resolve(), self.output)

    def test_success_replaces_previous_root_before_removing_pending_protection(self):
        root, previous = self.previous_root()

        def completed(command, **kwargs):
            self.assertEqual(root.resolve(), previous)
            if command[0] in ("uv", "nix-store", self.interpreter):
                pending = Path(self.process.call_args_list[0].args[0][6])
                self.assertTrue(pending.is_symlink())
                self.assertEqual(pending.resolve(), self.output)
            return self.nix_process(command, **kwargs)

        self.process.side_effect = completed
        wheel.main(["--backend", "nix", "gammaloop", "--install"])
        self.assertEqual(
            [call.args[0][0] for call in self.process.call_args_list],
            ["nix", "nix-store", self.interpreter, "uv", self.interpreter, "nix-store"],
        )
        self.assertTrue(root.is_symlink())
        self.assertEqual(root.resolve(), self.output)
        self.assertEqual(list(root.parent.iterdir()), [root])
        export_root = self.root / "target/wheels" / f".nix-{self.wheel.name}"
        self.assertEqual(export_root.resolve(), self.output)
        with (self.venv / ".gammaloop-install.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)

    def test_repeated_attempts_use_distinct_pending_roots(self):
        root, previous = self.previous_root()

        def completed(command, **kwargs):
            result = self.nix_process(command, **kwargs)
            if command[0] == "uv":
                result.returncode = 19
            return result

        self.process.side_effect = completed
        pending_roots = []
        for _ in range(2):
            self.process.reset_mock()
            with self.assertRaises(subprocess.CalledProcessError):
                wheel.main(["--backend", "nix", "gammaloop", "--install"])
            pending_roots.append(Path(self.process.call_args_list[0].args[0][6]))
            self.assertEqual(root.resolve(), previous)
        self.assertNotEqual(*pending_roots)

        self.process.reset_mock()
        self.process.side_effect = self.nix_process
        wheel.main(["--backend", "nix", "gammaloop", "--install"])
        successful_pending = Path(self.process.call_args_list[0].args[0][6])
        self.assertNotIn(successful_pending, pending_roots)
        self.assertFalse(successful_pending.is_symlink())
        self.assertEqual(root.resolve(), self.output)
        for pending in pending_roots:
            self.assertTrue(pending.is_symlink())
            self.assertEqual(pending.resolve(), self.output)


if __name__ == "__main__":
    unittest.main()
