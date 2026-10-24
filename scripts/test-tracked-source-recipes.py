#!/usr/bin/env python3
"""Run the public Just entry points against Nix/exporter process fixtures."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parent.parent
SOURCE = "path:/nix/store/public-fixture-source"


class RecipeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="tracked-recipe-fixture-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        (self.root / "bin").mkdir()
        (self.root / "scripts").mkdir()
        bash = shutil.which("bash")
        self.assertIsNotNone(bash, "Recipe tests require the pinned Bash")
        # Nix sandboxes lack /usr/bin/env. Pin only the copied recipes' interpreter;
        # their bodies and argument forwarding remain the actual repository code.
        recipes = (REPOSITORY / "just/ci.just").read_text()
        recipes = recipes.replace("#!/usr/bin/env bash", f"#!{bash}")
        (self.root / "justfile").write_text(
            "mod linnet 'crates/linnet/Justfile'\n" + recipes
        )
        (self.root / "crates/linnet").mkdir(parents=True)
        shutil.copyfile(
            REPOSITORY / "crates/linnet/Justfile", self.root / "crates/linnet/Justfile"
        )
        (self.root / "netrc").write_text("public fixture, no network\n")
        self.just = shutil.which("just")
        self.assertIsNotNone(self.just, "Recipe tests require the pinned Just")
        self.env = dict(
            os.environ,
            PATH=str(self.root / "bin") + os.pathsep + os.environ["PATH"],
            TMPDIR=str(self.root),
            FIXTURE_ROOT=str(self.root),
            NIXCI_NETRC=str(self.root / "netrc"),
            SYMBOLICA_LICENSE="public-explicit",
            SYMBOLICA_LICENSE_SIGNED="public-signed",
        )
        self.env.pop("CI_FLAKE_SOURCE", None)
        self.env.pop("PYTHONOPTIMIZE", None)
        (self.root / "bin/python3").symlink_to(sys.executable)
        (self.root / "scripts/tracked_source.py").write_text(
            "import os, sys\nfrom pathlib import Path\n"
            'root = Path(os.environ["FIXTURE_ROOT"])\n'
            'if "--flake-source" in sys.argv:\n'
            '    source = sys.argv[sys.argv.index("--flake-source") + 1]\n'
            f"    if source != {SOURCE!r}: sys.exit(29)\n"
            "else:\n"
            '    with (root / "exports").open("a") as out: out.write("exported\\n")\n'
            f"print({SOURCE!r})\n"
        )
        self.executable(
            "bin/nix",
            """
import json, os, sys
from pathlib import Path
assert "SYMBOLICA_LICENSE" not in os.environ
assert "SYMBOLICA_LICENSE_SIGNED" not in os.environ
root = Path(os.environ["FIXTURE_ROOT"])
args = sys.argv[1:]
with (root / "commands").open("a") as out:
    out.write(json.dumps(args) + "\\n")
if args[0] == "eval":
    print("fixture-system")
elif args[0] == "run":
    assert args[1] == "path:/nix/store/public-fixture-source#ci-runtime-license-boundary"
    assert args[args.index("--flake-source") + 1] == "path:/nix/store/public-fixture-source"
    sys.exit(int(os.environ.get("FIXTURE_BOUNDARY_STATUS", "0")))
elif args[0] == "build":
    assert "--expr" in args or any(arg.startswith("path:/nix/store/public-fixture-source#") for arg in args)
    print(root / "launcher")
elif args[0] == "print-dev-env":
    assert args[1] == "path:/nix/store/public-fixture-source"
    print("export CI_FIXTURE_SHELL=prepared")
    print("shellHook='export CI_FIXTURE_HOOK_COUNT=$(( ${CI_FIXTURE_HOOK_COUNT:-0} + 1 ))'")
    print('eval "${shellHook:-}"')
elif args[0] == "flake":
    assert args[1:3] == ["check", "path:/nix/store/public-fixture-source"]
elif args[0] == "copy":
    assert sys.stdin.read().strip() == str(root / "launcher")
else:
    raise AssertionError(args)
""",
        )
        self.executable(
            "launcher/bin/fixture-runtime",
            """
import json, os, sys
from pathlib import Path
assert os.environ["SYMBOLICA_LICENSE"] == "public-explicit"
assert os.environ["SYMBOLICA_LICENSE_SIGNED"] == "public-signed"
root = Path(os.environ["FIXTURE_ROOT"])
(root / "runtime").write_text(json.dumps(sys.argv[1:]))
with (root / "runtime_calls").open("a") as out:
    out.write(Path(sys.argv[0]).name + "\\n")
sys.exit(int(os.environ.get("FIXTURE_RUNTIME_STATUS", "0")))
""",
        )
        for package in (
            "nix-ci-check-gammaloop-doctest",
            "nix-ci-check-gammaloop-nextest",
            "nix-ci-check-alphal00p-docs",
        ):
            (self.root / "launcher/bin" / package).symlink_to("fixture-runtime")
        self.executable(
            "bin/check-development-shell",
            """
import json, os, sys
from pathlib import Path
assert os.environ["CI_FIXTURE_SHELL"] == "prepared"
assert os.environ["CI_FIXTURE_HOOK_COUNT"] == "1"
assert os.environ["SYMBOLICA_LICENSE"] == "public-explicit"
assert os.environ["SYMBOLICA_LICENSE_SIGNED"] == "public-signed"
(Path(os.environ["FIXTURE_ROOT"]) / "runtime").write_text(json.dumps(sys.argv[1:]))
""",
        )

    def executable(self, relative, body):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"#!{sys.executable}\n" + body)
        path.chmod(0o755)

    def recipe(self, *arguments, **environment):
        return subprocess.run(
            [self.just, "--justfile", str(self.root / "justfile"), *arguments],
            cwd=self.root,
            env=self.env | environment,
            capture_output=True,
            text=True,
        )

    def commands(self):
        return [
            json.loads(line)
            for line in (self.root / "commands").read_text().splitlines()
        ]

    def assert_success(self, result):
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_linnet_wasm_uses_the_tracked_source_wrapper(self):
        self.assert_success(self.recipe("linnet::wasm"))
        self.assertEqual(
            self.commands(), [["build", SOURCE + "#linnest-wasm", "--print-build-logs"]]
        )
        self.assertEqual((self.root / "exports").read_text().splitlines(), ["exported"])

    def test_guard_reuses_one_export_for_build_and_probe(self):
        self.assert_success(
            self.recipe("ci-runtime-check", CI_FLAKE_SOURCE="stale-source")
        )
        self.assertEqual((self.root / "exports").read_text().splitlines(), ["exported"])
        self.assertEqual(
            [args[0] for args in self.commands()], ["eval", "build", "run"]
        )

    def test_guard_reuses_only_an_explicit_snapshot(self):
        self.assert_success(
            self.recipe("ci-runtime-check", SOURCE, CI_FLAKE_SOURCE="stale-source")
        )
        self.assertFalse((self.root / "exports").exists())
        self.assertEqual(
            [args[0] for args in self.commands()], ["eval", "build", "run"]
        )

    def test_ci_checks_share_one_tracked_source_and_execute_all_runtime_checks(self):
        # A top-level check must snapshot the current checkout even if its
        # caller inherited an earlier source; only child recipes reuse it.
        self.assert_success(self.recipe("ci-checks", CI_FLAKE_SOURCE="stale-source"))
        self.assertEqual((self.root / "exports").read_text().splitlines(), ["exported"])
        commands = self.commands()
        self.assertEqual([args[0] for args in commands[:3]], ["eval", "build", "run"])
        self.assertEqual(
            [
                arg
                for args in commands
                if args[0] == "build"
                for arg in args
                if "#" in arg
            ],
            [
                SOURCE + "#checks.fixture-system.gammaloop-ci-runtime",
                SOURCE + "#checks.fixture-system.gammaloop-clippy",
                SOURCE + "#checks.fixture-system.gammaloop-fmt",
                SOURCE + "#checks.fixture-system.gammaloop-guppy-workspace-graph",
                SOURCE + "#nix-ci-check-gammaloop-doctest",
                SOURCE + "#nix-ci-check-gammaloop-nextest",
                SOURCE + "#nix-ci-check-alphal00p-docs",
            ],
        )
        self.assertEqual(
            (self.root / "runtime_calls").read_text().splitlines(),
            [
                "nix-ci-check-gammaloop-doctest",
                "nix-ci-check-gammaloop-nextest",
                "nix-ci-check-alphal00p-docs",
            ],
        )

    def test_ci_checks_stop_before_static_and_runtime_jobs_on_probe_failure(self):
        result = self.recipe("ci-checks", FIXTURE_BOUNDARY_STATUS="37")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(
            [args[0] for args in self.commands()], ["eval", "build", "run"]
        )
        self.assertFalse((self.root / "runtime_calls").exists())

    def test_reused_source_is_validated_before_nix(self):
        self.assert_success(self.recipe("nix-build", "fixture", CI_FLAKE_SOURCE=SOURCE))
        self.assertFalse((self.root / "exports").exists())
        before = self.commands()
        result = self.recipe("nix-build", "fixture", CI_FLAKE_SOURCE=".")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.commands(), before)

    def test_boundary_failure_prevents_upload_and_producer_build(self):
        result = self.recipe("ci-checks-and-upload", FIXTURE_BOUNDARY_STATUS="37")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual([args[0] for args in self.commands()], ["run"])
        self.assertNotIn("CI upload completed", result.stdout)

    def test_passing_boundary_precedes_build_and_upload(self):
        result = self.recipe("ci-checks-and-upload")
        self.assert_success(result)
        self.assertEqual(
            [args[0] for args in self.commands()], ["run", "build", "copy"]
        )
        self.assertEqual((self.root / "exports").read_text().splitlines(), ["exported"])
        self.assertIn("CI upload completed", result.stdout)

    def test_runtime_keeps_keys_out_of_build_and_preserves_arguments(self):
        self.assert_success(
            self.recipe(
                "nix-run", "fixture-runtime", "--filter", "test(=public fixture)"
            )
        )
        self.assertEqual(
            json.loads((self.root / "runtime").read_text()),
            ["--filter", "test(=public fixture)"],
        )
        self.assertEqual([args[0] for args in self.commands()], ["build"])

    def test_runtime_failure_is_not_replaced_by_build_success(self):
        result = self.recipe("nix-run", "fixture-runtime", FIXTURE_RUNTIME_STATUS="37")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("exit code 37", result.stderr)
        self.assertTrue((self.root / "runtime").exists())

    def test_development_environment_restores_keys_only_at_runtime(self):
        self.assert_success(
            self.recipe(
                "nix-develop", "check-development-shell", "two words", "test(=fixture)"
            )
        )
        self.assertEqual([args[0] for args in self.commands()], ["print-dev-env"])
        self.assertEqual(
            json.loads((self.root / "runtime").read_text()),
            ["two words", "test(=fixture)"],
        )

    def test_pure_build_and_flake_check_use_tracked_source(self):
        self.assert_success(self.recipe("nix-build", "fixture", "--no-link"))
        self.assert_success(self.recipe("nix-flake-check", "--no-build"))
        self.assertEqual(
            self.commands(),
            [
                ["build", SOURCE + "#fixture", "--no-link"],
                ["flake", "check", SOURCE, "--no-build"],
            ],
        )


if __name__ == "__main__":
    unittest.main()
