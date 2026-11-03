#!/usr/bin/env python3
"""Inspect publication coverage without realizing selected outputs or running tests.

The shell backend supplies a live tracked-source snapshot. This module only
performs offline evaluation and local store queries; IFD is disabled.
"""

import argparse
import json
import os
import platform
import subprocess
import sys
from pathlib import Path
from urllib.parse import quote


class PublicationPreview:
    def __init__(self, root):
        self.root = root
        # Ignore caller-injected Nix settings, including credential paths/hooks.
        self.env = os.environ.copy()
        for key in (
            "SYMBOLICA_LICENSE",
            "SYMBOLICA_LICENSE_SIGNED",
            "GAMMALOOP_NIX_TEST_ARGS",
            "GAMMALOOP_DEVELOP_SPEC",
            "NIX_CONFIG",
            "NIX_USER_CONF_FILES",
            "NIXCI_NETRC",
            "NIX_REMOTE",
        ):
            self.env.pop(key, None)
        self.nix = [
            "nix",
            "--extra-experimental-features",
            "nix-command flakes",
            "--offline",
            "--option",
            "allow-import-from-derivation",
            "false",
            "--option",
            "substituters",
            "",
            "--option",
            "netrc-file",
            "/dev/null",
            "--option",
            "post-build-hook",
            "",
        ]

    def command(self, arguments, *, valid=False):
        result = subprocess.run(
            arguments,
            cwd=self.root,
            env=self.env,
            text=True,
            capture_output=True,
            check=False,
        )
        if valid:
            return result.returncode == 0
        if result.returncode:
            raise ValueError(f"Read-only query failed: {result.stderr.strip()}")
        return result.stdout

    def inspect(self, scope):
        machine = {"arm64": "aarch64", "AMD64": "x86_64"}.get(
            platform.machine(), platform.machine()
        )
        system = f"{machine}-{platform.system().lower()}"
        supplied = os.environ.get("GAMMALOOP_FLAKE_SOURCE")
        if not supplied:
            raise ValueError(
                "Publication preview requires GAMMALOOP_FLAKE_SOURCE from the shell backend"
            )
        source = Path(supplied).resolve()
        if not source.is_dir() or not (source / "flake.nix").is_file():
            raise ValueError("GAMMALOOP_FLAKE_SOURCE must name a flake snapshot")
        flake = f"path:{quote(str(source), safe='/')}"
        expression = f"""let
          f = builtins.getFlake {json.dumps(flake)};
          m = f.lib.{json.dumps(system)}.cachePublication;
          attrs = builtins.attrNames (builtins.listToAttrs
            (map (name: {{ inherit name; value = true; }}) (m.local ++ m.github)));
          resolve = attr: builtins.foldl' (value: name: value.${{name}}) f
            (builtins.filter builtins.isString (builtins.split "[.]" attr));
        in {{ metadata = m; outputs = builtins.listToAttrs (map (name:
          let d = resolve name; in {{ inherit name; value = map
            (output: {{ inherit output; path = builtins.unsafeDiscardStringContext d.${{output}}.outPath; }}) d.outputs;
          }}) attrs); }}"""
        data = json.loads(
            self.command(
                self.nix
                + [
                    "eval",
                    "--impure",
                    "--no-write-lock-file",
                    "--json",
                    "--expr",
                    expression,
                ]
            )
        )
        metadata, outputs = data["metadata"], data["outputs"]
        source_dependent = set(metadata["sourceRevisionDependent"])
        validity = {
            row["path"]: self.command(
                ["nix-store", "--check-validity", row["path"]], valid=True
            )
            for attr, rows in outputs.items()
            if attr not in source_dependent
            for row in rows
        }
        closures = {}
        print(
            "Read-only coverage; launcher realization does not establish passing tests."
        )
        for name in ("local", "github"):
            if name == "local" and scope == "github":
                continue
            excluded = [attr for attr in metadata[name] if attr in source_dependent]
            rows = [
                (attr, row)
                for attr in metadata[name]
                if attr not in source_dependent
                for row in outputs[attr]
            ]
            realized = sorted({row["path"] for _, row in rows if validity[row["path"]]})
            missing = sum(not validity[row["path"]] for _, row in rows)
            info = (
                json.loads(
                    self.command(
                        self.nix + ["path-info", "--json", "--recursive", *realized]
                    )
                )
                if realized
                else {}
            )
            if isinstance(info, list):
                info = {row["path"]: row for row in info}
            if any(row is None for row in info.values()):
                raise ValueError(
                    "Store paths disappeared during preview; retry the read-only query"
                )
            closures[name] = set(info)
            if scope not in (name, "both"):
                continue
            print(
                f"{name}: {len(rows) - missing}/{len(rows)} comparable root outputs realized; "
                f"{'PARTIAL coverage' if missing else 'complete comparable coverage'}"
            )
            for attr, row in rows:
                print(
                    f"  {attr}.{row['output']} {row['path']} "
                    f"{'realized' if validity[row['path']] else 'missing'}"
                )
            for attr in excluded:
                for row in outputs[attr]:
                    print(
                        f"  {attr}.{row['output']} {row['path']} "
                        "not comparable across source identities"
                    )
            if excluded:
                print(
                    "  Entire publication size unknown: source-dependent roots excluded"
                )
            size = sum(row["narSize"] for row in info.values())
            print(
                f"  Comparable realized-root subtotal: {len(info)} deduplicated closure paths; {size} raw NAR bytes"
            )
            if name == "github":
                print(
                    f"  Publication cap: < {metadata['githubMaxNarBytes']} raw NAR bytes (unchanged)"
                )
        github_missing = any(
            attr in source_dependent for attr in metadata["github"]
        ) or any(
            not validity[row["path"]]
            for attr in metadata["github"]
            if attr not in source_dependent
            for row in outputs[attr]
        )
        print(
            "Selected CI producers outside realized GitHub requisites "
            f"({'full set unknown: missing GitHub roots' if github_missing else 'confirmed absent from full set'}):"
        )
        for attr in metadata["producers"]:
            if attr in source_dependent:
                continue
            for row in outputs[attr]:
                if row["path"] not in closures["github"]:
                    print(
                        f"  {attr}.{row['output']} {row['path']} "
                        f"({'realized' if validity[row['path']] else 'missing'})"
                    )


def preview(root, arguments):
    parser = argparse.ArgumentParser(prog="ci-upload-preview", allow_abbrev=False)
    parser.add_argument(
        "scope", nargs="?", choices=("local", "github", "both"), default="both"
    )
    options = parser.parse_args(arguments)
    PublicationPreview(root).inspect(options.scope)


if __name__ == "__main__":
    try:
        root = Path(
            os.environ.get(
                "GAMMALOOP_WORKSPACE_ROOT", Path(__file__).resolve().parents[1]
            )
        ).resolve()
        preview(root, sys.argv[1:])
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(2)
