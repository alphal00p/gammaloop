#!/usr/bin/env python3
"""Export only tracked working-tree files before Nix sees a checkout."""

import argparse
import json
import io
import os
import shutil
import shlex
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


LICENSE_VARIABLES = ("SYMBOLICA_LICENSE", "SYMBOLICA_LICENSE_SIGNED")


def without_credentials(env=None):
    """Keep user configuration and strip only runtime license variables."""
    return {
        key: value
        for key, value in (os.environ if env is None else env).items()
        if key not in LICENSE_VARIABLES
    }


def command_error(error):
    """Keep actionable command diagnostics without echoing inherited keys."""
    if isinstance(error, subprocess.CalledProcessError):
        message = f"{shlex.join(error.cmd)} exited with status {error.returncode}"
        if error.stderr:
            message += "\n" + error.stderr.decode(errors="replace").strip()
    else:
        message = str(error)
    for key in LICENSE_VARIABLES:
        if value := os.environ.get(key):
            message = message.replace(value, "[redacted]")
    return message


def run(command, env):
    return subprocess.run(command, env=env, check=True, capture_output=True).stdout


def copy_tracked_files(root, destination, names):
    """Copy only listed entries, without following workspace symlinks."""
    for name in names:
        relative = Path(name)
        if relative.is_absolute() or not relative.parts or ".." in relative.parts:
            raise ValueError(f"Tracked source contains an invalid path: {name!r}")
        source = root / relative
        if source.parent.resolve(strict=True) != source.parent:
            raise ValueError(f"Tracked source has a symlink parent: {name!r}")
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        mode = source.lstat().st_mode
        if stat.S_ISLNK(mode):
            link = os.readlink(source)
            relative_target = Path(os.path.normpath(relative.parent / link))
            if (
                relative_target.is_absolute()
                or ".." in relative_target.parts
                or not source.resolve().is_relative_to(root)
            ):
                raise ValueError(f"Tracked source has an escaping symlink: {name!r}")
            target.symlink_to(link)
        elif stat.S_ISREG(mode):
            # Refuse a leaf swapped to a symlink after lstat. Never recursively
            # copy a tracked directory and accidentally admit ignored children.
            descriptor = os.open(source, os.O_RDONLY | os.O_NOFOLLOW)
            with os.fdopen(descriptor, "rb") as incoming, target.open("wb") as outgoing:
                shutil.copyfileobj(incoming, outgoing)
            target.chmod(stat.S_IMODE(mode))
        else:
            raise ValueError(f"Tracked source contains a non-file entry: {name!r}")
    if not (destination / "flake.nix").is_file():
        raise ValueError("Tracked source has no flake.nix")


def safe_flake_source(root, explicit, env):
    """Freeze tracked files before Nix can import any workspace directory."""
    env = without_credentials(env)
    if explicit is not None:
        source = Path(explicit.removeprefix("path:"))
        if not source.is_absolute():
            raise ValueError("--flake-source requires an existing Nix-store snapshot")
        source = source.resolve(strict=True)
        if (
            not source.is_relative_to("/nix/store")
            or source == Path("/nix/store")
            or not (source / "flake.nix").is_file()
        ):
            raise ValueError(
                "--flake-source requires an existing Nix-store flake snapshot"
            )
        return "path:" + str(source)
    root = root.resolve(strict=True)
    if (root / ".jj").exists():
        # JSON handles spaces, newlines and quotes without jj's display quoting.
        # Snapshot with the user's VCS ignore configuration still in effect.
        names = [
            json.loads(line)
            for line in run(
                ["jj", "--no-pager", "file", "list", "-T", 'json(path) ++ "\\n"'], env
            ).splitlines()
        ]
    else:
        git_root = Path(
            run(["git", "rev-parse", "--show-toplevel"], env).decode().strip()
        ).resolve(strict=True)
        if git_root != root:
            raise ValueError("Run from the Git worktree root or supply --flake-source")
        # Include dirty tracked files and staged additions, never untracked or
        # ignored files. Missing tracked entries fail the export before Nix runs.
        names = [
            os.fsdecode(name)
            for name in run(["git", "ls-files", "--cached", "-z"], env).split(b"\0")
            if name
        ]
    with tempfile.TemporaryDirectory(prefix="gammaloop-tracked-source-") as tmp:
        temporary = Path(tmp)
        source = temporary / "source"
        source.mkdir()
        copy_tracked_files(root, source, sorted(set(names)))
        stored = (
            run(
                ["nix", "store", "add-path", "--name", "source", str(source)],
                env,
            )
            .decode()
            .strip()
        )
    return "path:" + stored


class SourceSelectionTests(unittest.TestCase):
    """Check the exact export admitted to Nix, including ignored credentials."""

    def test_jj_export_excludes_ignored_credentials_and_preserves_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / ".jj").mkdir()
            (root / ".envrc").write_text("ignored fixture credential\n")
            (root / "config").mkdir()
            (root / "config/license").write_text("ignored nested credential\n")
            names = ["flake.nix", 'name with "quotes"\nand spaces', "link"]
            (root / names[0]).write_text("{}\n")
            (root / names[1]).write_text("tracked edit\n")
            (root / "link").symlink_to(names[1])
            calls = []

            def execute(command, env):
                calls.append(command)
                if command[0] == "jj":
                    self.assertEqual(command[-1], 'json(path) ++ "\\n"')
                    return "\n".join(map(json.dumps, names)).encode()
                self.assertEqual(
                    command[:5], ["nix", "store", "add-path", "--name", "source"]
                )
                source = Path(command[-1])
                self.assertEqual(
                    sorted(p.name for p in source.iterdir()), sorted(names)
                )
                self.assertEqual((source / names[1]).read_text(), "tracked edit\n")
                self.assertTrue((source / "link").is_symlink())
                self.assertNotIn("SYMBOLICA_LICENSE", env)
                self.assertNotIn("SYMBOLICA_LICENSE_SIGNED", env)
                self.assertEqual(env["HOME"], os.environ["HOME"])
                self.assertEqual(env["XDG_CONFIG_HOME"], "/fixture/config")
                self.assertEqual(env["NIX_CONFIG"], "fixture-public-config")
                return b"/nix/store/fixture-source\n"

            with (
                patch(__name__ + ".run", side_effect=execute),
                patch.dict(
                    os.environ,
                    XDG_CONFIG_HOME="/fixture/config",
                    NIX_CONFIG="fixture-public-config",
                    SYMBOLICA_LICENSE="fixture",
                    SYMBOLICA_LICENSE_SIGNED="fixture",
                ),
            ):
                self.assertEqual(
                    safe_flake_source(root, None, without_credentials()),
                    "path:/nix/store/fixture-source",
                )
            self.assertEqual(len(calls), 2)
            self.assertFalse(
                Path(calls[-1][-1]).exists(), "Export scratch survived cleanup"
            )

    def test_git_export_uses_nul_delimited_tracked_entries(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / ".git").write_text("gitdir: fixture\n")
            (root / "flake.nix").write_text("{}\n")
            (root / ".env").write_text("untracked fixture credential\n")
            with patch(
                __name__ + ".run",
                side_effect=[
                    str(root).encode(),
                    b"flake.nix\0",
                    b"/nix/store/fixture-source\n",
                ],
            ) as execute:
                self.assertEqual(
                    safe_flake_source(root, None, {}), "path:/nix/store/fixture-source"
                )
            self.assertEqual(
                execute.call_args_list[1].args[0], ["git", "ls-files", "--cached", "-z"]
            )

    def test_vcs_failure_never_imports_a_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for is_jj in (False, True):
                if is_jj:
                    (root / ".jj").mkdir()
                with (
                    self.subTest(is_jj=is_jj),
                    patch(
                        __name__ + ".run",
                        side_effect=subprocess.CalledProcessError(1, ["vcs"]),
                    ) as execute,
                ):
                    with self.assertRaises(subprocess.CalledProcessError):
                        safe_flake_source(root, None, {})
                    self.assertEqual(execute.call_count, 1)
            (root / ".jj").rmdir()
            with patch(
                __name__ + ".run", return_value=str(root.parent).encode()
            ) as execute:
                with self.assertRaises(ValueError):
                    safe_flake_source(root, None, {})
                self.assertEqual(execute.call_count, 1)

    def test_invalid_tracked_entries_fail_before_store_import(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / ".jj").mkdir()
            (root / "flake.nix").write_text("{}\n")
            (root / "directory").mkdir()
            (root / "escape").symlink_to(root.parent, target_is_directory=True)
            (root / "alias").symlink_to(root, target_is_directory=True)
            (root / "absolute").symlink_to(root / "flake.nix")
            for name in (
                "missing",
                "directory",
                "escape",
                "absolute",
                "alias/flake.nix",
                "../flake.nix",
                str(root / "flake.nix"),
            ):
                with (
                    self.subTest(name=name),
                    patch(
                        __name__ + ".run", return_value=json.dumps(name).encode()
                    ) as execute,
                ):
                    with self.assertRaises((ValueError, OSError)):
                        safe_flake_source(root, None, {})
                    self.assertEqual(execute.call_count, 1)

    def test_explicit_checkout_paths_and_symlink_escapes_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "flake.nix").write_text("{}\n")
            link = root / "alias"
            link.symlink_to(root, target_is_directory=True)
            for source in (
                str(root),
                "path:" + str(root),
                str(link),
                ".",
                "path:.",
                "/nix/store/../../tmp",
                "/nix/store-fixture",
                "/nix/store",
            ):
                with self.subTest(source=source), patch(__name__ + ".run") as execute:
                    with self.assertRaises((ValueError, FileNotFoundError)):
                        safe_flake_source(root, source, {})
                    execute.assert_not_called()

    def test_existing_store_snapshot_is_accepted(self):
        source = Path("/nix/store/" + "a" * 32 + "-source")
        with (
            patch.object(Path, "resolve", return_value=source),
            patch.object(Path, "is_file", return_value=True),
            patch(__name__ + ".run") as execute,
        ):
            self.assertEqual(
                safe_flake_source(Path.cwd(), "path:" + str(source), {}),
                "path:" + str(source),
            )
            execute.assert_not_called()

    def test_export_errors_include_command_stderr_and_remediation(self):
        error = subprocess.CalledProcessError(
            17,
            ["nix", "store", "add-path"],
            stderr=b"fixture: enable nix-command; public-key-fixture",
        )
        with (
            patch.dict(os.environ, SYMBOLICA_LICENSE="public-key-fixture"),
            patch(__name__ + ".safe_flake_source", side_effect=error),
            patch.object(sys, "argv", ["tracked_source.py", "--print-source"]),
            patch.object(sys, "stderr", new_callable=io.StringIO) as diagnostic,
            self.assertRaises(SystemExit) as failed,
        ):
            main()
        self.assertEqual(failed.exception.code, 2)
        message = diagnostic.getvalue()
        self.assertIn("nix store add-path exited with status 17", message)
        self.assertIn("fixture: enable nix-command", message)
        self.assertIn("--flake-source", message)
        self.assertNotIn("public-key-fixture", message)

    def test_only_credential_environment_entries_are_removed(self):
        environment = {
            "HOME": "/public/home",
            "XDG_CONFIG_HOME": "/public/config",
            "NIX_USER_CONF_FILES": "/public/nix.conf",
            "NIX_CONFIG": "experimental-features = nix-command flakes",
            "SYMBOLICA_LICENSE": "public-explicit",
            "SYMBOLICA_LICENSE_SIGNED": "public-signed",
        }
        expected = {
            key: value
            for key, value in environment.items()
            if key not in LICENSE_VARIABLES
        }
        self.assertEqual(without_credentials(environment), expected)
        self.assertEqual(environment["SYMBOLICA_LICENSE"], "public-explicit")


def main():
    if sys.argv[1:] == ["--self-test"]:
        unittest.main(argv=[sys.argv[0]])
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--print-source",
        action="store_true",
        required=True,
        help="Export tracked Git/jj source to the Nix store without flake evaluation",
    )
    parser.add_argument(
        "--flake-source", help="Reuse an existing Nix-store flake snapshot"
    )
    args = parser.parse_args()
    try:
        print(safe_flake_source(Path.cwd(), args.flake_source, without_credentials()))
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.error(
            "Cannot export tracked source: "
            + command_error(error)
            + "\nRun at the repository root with Python 3, Git/jj and Nix available, "
            "or supply --flake-source with an existing Nix-store flake snapshot."
        )


if __name__ == "__main__":
    main()
