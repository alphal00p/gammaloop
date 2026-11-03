"""Exercise GammaLoop's shell integration against the dyndrv worker contracts."""

import argparse
import os
import shutil
import sys
import tempfile
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[2]
arguments = argparse.ArgumentParser(add_help=False)
arguments.add_argument(
    "--bootstrap", type=Path, default=os.environ.get("GAMMALOOP_DYNDRV_BOOTSTRAP")
)
bootstrap, _ = arguments.parse_known_args()
if bootstrap.bootstrap is None:
    raise SystemExit("Use just test-dyndrv, or supply --bootstrap")
sys.path.insert(0, str(bootstrap.bootstrap / "worker-tests"))
from test_integration import Harness as WorkerHarness

TOOLING = (
    "scripts/develop.sh",
    "scripts/develop_dyndrv.sh",
    "nix/dyndrv.nix",
    "nix/patches/dyndrv-native-run-path.patch",
    "nix/dyndrv-bootstrap.nix",
    "flake.lock",
    ".config/build-profile",
)


class Harness(WorkerHarness):
    repository = REPOSITORY
    tooling = TOOLING
    immutable_names = tuple(name for name in TOOLING if name.startswith("nix/")) + (
        "flake.lock",
    )

    def __init__(self, args):
        if not os.environ.get("SYMBOLICA_LICENSE"):
            license_file = (
                Path(os.environ.get("XDG_CONFIG_HOME", str(Path.home() / ".config")))
                / "symbolica/license"
            )
            license_value = os.environ.get("SYMBOLICA_LICENSE_SIGNED", "")
            if not license_value and license_file.is_file():
                license_value = license_file.read_text().rstrip("\n")
            if license_value:
                os.environ["SYMBOLICA_LICENSE"] = license_value
        super().__init__(args)

    @staticmethod
    def parser(description, directory, **options):
        parser = WorkerHarness.parser(description, directory, **options)
        parser.set_defaults(
            directory=REPOSITORY / "target/dyndrv-validation" / directory
        )
        return parser

    @staticmethod
    def copy_tooling(destination):
        WorkerHarness.copy_tooling(destination)
        for name in TOOLING:
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(REPOSITORY / name, target)

    def invocation(self, workspace, action, arguments, backend="dyndrv"):
        if action in {"test", "list"}:
            arguments = (*arguments, "--profile", "ci_dyndrv")
        return [
            "bash",
            workspace / "scripts/develop.sh",
            f"--backend={backend}",
            action,
            *arguments,
        ]

    def prepare(self):
        super().prepare()
        self.verify_vendor_provider()

    def verify_vendor_provider(self):
        if self.spec.get("vendor_lockfile_only") is not True:
            return
        with (
            self.step("vendor provider ignores caller configuration and overlays"),
            tempfile.TemporaryDirectory(dir=self.run_directory) as temporary,
        ):
            home = Path(temporary)
            environment = {
                name: value
                for name, value in self.environment.items()
                if name not in self.spec.get("private_environment_names", [])
            } | {
                "HOME": str(home),
                "NIX_PATH": "",
                "NIXPKGS_CONFIG": "",
                "CARGO_DYNDRV_CARGO_LOCK": str(self.a / "Cargo.lock"),
            }
            command = [
                self.spec["nix"],
                "eval",
                "--impure",
                "--raw",
                "--expr",
                f"({self.spec['vendor_expression']}).drvPath",
            ]
            expected, _, _ = self.command(
                "vendor-provider-default", command, env=environment
            )
            config = home / ".config/nixpkgs"
            config.mkdir(parents=True)
            for name in ("config.nix", "overlays.nix"):
                (config / name).write_text('throw "caller configuration was read"\n')
            environment["NIXPKGS_CONFIG"] = str(config / "config.nix")
            environment["NIX_PATH"] = f"nixpkgs-overlays={config / 'overlays.nix'}"
            actual, _, _ = self.command(
                "vendor-provider-poisoned", command, env=environment
            )
            self.require(
                actual == expected,
                "Caller configuration changed the vendor provider",
            )


if __name__ == "__main__":
    sys.exit(Harness.run_cli(Harness.parser(__doc__, "integration", timeout=120)))
