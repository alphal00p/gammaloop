"""Exercise the public launcher after Nix and explicit Cargo builds."""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="gammaloop launcher ")
        self.root = Path(self.temporary.name)
        shutil.copy2(Path(__file__).resolve().parents[2] / "gammaloop", self.root)

    def tearDown(self):
        self.temporary.cleanup()

    def executable(self, path, label, modified):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("#!/bin/sh\nprintf '%s\\n' '" + label + '\' "$@"\n')
        path.chmod(0o755)
        os.utime(path, (modified, modified))

    def cargo(self, profile, modified=100):
        self.executable(self.root / "target" / profile / "gammaloop", profile, modified)

    def nix(self, modified=200, broken=False):
        package = self.root / "store product"
        if not broken:
            self.executable(package / "bin" / "gammaloop", "nix", 1)
        directory = self.root / "target" / "nix-dev"
        directory.mkdir(parents=True)
        (directory / "gammaloop").symlink_to(package)
        stamp = directory / "gammaloop.built"
        stamp.touch()
        os.utime(stamp, (modified, modified))

    def run_launcher(self, *args):
        return subprocess.run(
            [str(self.root / "gammaloop"), *args],
            capture_output=True,
            text=True,
            check=False,
        )

    def test_clean_target_uses_nix_product_and_preserves_arguments(self):
        self.nix()
        result = self.run_launcher("--help", "argument with spaces", "$(literal)")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            result.stdout.splitlines(),
            ["nix", "--help", "argument with spaces", "$(literal)"],
        )

    def test_fresh_nix_selection_overrides_stale_cargo_binary(self):
        self.cargo("release")
        self.nix()
        self.assertEqual(self.run_launcher().stdout.strip(), "nix")

    def test_later_cargo_build_wins_implicit_selection(self):
        self.nix()
        self.cargo("release", modified=300)
        self.cargo("dev-optim", modified=400)
        self.assertEqual(self.run_launcher().stdout.strip(), "dev-optim")

    def test_explicit_profile_overrides_newer_nix_selection(self):
        self.cargo("release")
        self.nix()
        self.assertEqual(self.run_launcher("--release").stdout.strip(), "release")

    def test_broken_current_nix_product_does_not_run_stale_cargo(self):
        self.cargo("release")
        self.nix(broken=True)
        result = self.run_launcher()
        self.assertEqual(result.returncode, 1)
        self.assertEqual(result.stdout, "")
        self.assertIn("target/nix-dev/gammaloop/bin/gammaloop", result.stderr)

    def test_cargo_only_selection_remains_available(self):
        self.cargo("release")
        self.cargo("dev-optim", modified=200)
        self.assertEqual(self.run_launcher().stdout.strip(), "dev-optim")


if __name__ == "__main__":
    unittest.main()
