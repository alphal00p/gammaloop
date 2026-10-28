#!/usr/bin/env python3
"""Probe only new pure inputs with dummy credentials, never a user's key."""

import argparse
import io
import json
import os
import secrets
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tracked_source import (
    LICENSE_VARIABLES,
    command_error,
    run,
    safe_flake_source,
    without_credentials,
)


REQUIRED_ROOTS = (
    "alphal00p-docs-build",
    "ci-runtime-environment",
    "ci-doctest-artifacts",
    "alphal00p-docs-cargo-artifacts",
    "cargoCheckArtifacts",
)


def probe_environment(home):
    env = dict(
        without_credentials(), HOME=str(home), XDG_CONFIG_HOME=str(home / "config")
    )
    env["NIX_CONFIG"] = (
        env.get("NIX_CONFIG", "") + "\nextra-experimental-features = nix-command flakes"
    )
    return env


def checked_environment(derivation):
    """Reject evaluated license fields and requests for impure build credentials."""
    for container in ("env", "structuredAttrs"):
        fields = derivation.get(container, {})
        if set(LICENSE_VARIABLES).intersection(fields):
            raise ValueError("License derivation field")
        impure = fields.get("impureEnvVars", [])
        if isinstance(impure, str):
            impure = impure.split()
        if set(LICENSE_VARIABLES).intersection(impure):
            raise ValueError("License in impureEnvVars")


def check_boundary(source_ref, report_dir):
    # The caller supplies an immutable tracked source. Never import a checkout.
    with tempfile.TemporaryDirectory(prefix="gammaloop-license-metadata-") as home:
        env = probe_environment(Path(home))
        metadata = json.loads(
            run(
                [
                    "nix",
                    "flake",
                    "metadata",
                    "--json",
                    "--no-write-lock-file",
                    source_ref,
                ],
                env,
            )
        )
    source = metadata["path"]
    expr = """let
      f = builtins.getFlake @FLAKE@;
      p = f.packages.x86_64-linux;
      c = f.checks.x86_64-linux;
      selected = builtins.filter (name:
        builtins.match "nix-ci-check-.*" name != null
        || builtins.elem name [@REQUIRED_ROOTS@])
        (builtins.attrNames p);
      archives = builtins.filter (name: builtins.match "gammaloop-nextest-binaries.*" name != null)
        (builtins.attrNames c);
      identity = d: { drv = d.drvPath; out = d.outPath; };
    in if f ? hydraJobsMain then throw "Legacy licensed Hydra jobs are exposed"
    else if p ? alphal00p-docs-pages || p ? alphal00p-docs-snapshot-fixture
      then throw "Legacy licensed documentation packages are exposed"
    else if c ? gammaloop-doctest || c ? gammaloop-nextest || c ? alphal00p-docs
      then throw "Legacy licensed test checks are exposed"
    else builtins.listToAttrs
      ((map (name: {inherit name; value = identity p.${name};}) selected)
      ++ (map (name: {inherit name; value = identity c.${name};}) archives))
    """.replace("@FLAKE@", json.dumps("path:" + source)).replace(
        "@REQUIRED_ROOTS@", " ".join(map(json.dumps, REQUIRED_ROOTS))
    )
    # Fresh public markers cannot match the guard's own packaged self-tests.
    # A uses the prescribed unsigned-only environment; B supplies both forms.
    # Neither may affect identity or appear in recursively inspected derivations.
    markers = tuple(
        f"PR3_BOUNDARY_DUMMY_{secrets.token_hex(16)}_{suffix}" for suffix in ("A", "B")
    )
    identities = []
    with tempfile.TemporaryDirectory(prefix="gammaloop-license-probe-") as home:
        for index, marker in enumerate(markers):
            env = probe_environment(Path(home))
            env["SYMBOLICA_LICENSE"] = marker
            if index == 1:
                env["SYMBOLICA_LICENSE_SIGNED"] = marker
            identity = json.loads(
                run(["nix", "eval", "--impure", "--json", "--expr", expr], env)
            )
            if not isinstance(identity, dict) or not identity:
                raise ValueError("No producer or launcher identities were inspected")
            missing_roots = set(REQUIRED_ROOTS).difference(identity)
            if missing_roots:
                raise ValueError(
                    "Missing required producer or launcher identities: "
                    + ", ".join(sorted(missing_roots))
                )
            identities.append(identity)
            drvs = sorted({entry["drv"] for entry in identity.values()})
            serialized = run(["nix", "derivation", "show", "--recursive", *drvs], env)
            if any(value.encode() in serialized for value in markers):
                raise ValueError("Dummy marker in derivation metadata")
            data = json.loads(serialized)
            # Nix 2.34 wraps derivations in {version, derivations}; older
            # versions return a mapping directly.
            derivations = data.get("derivations", data)
            if not isinstance(derivations, dict) or not derivations:
                raise ValueError("No derivation metadata was inspected")
            for derivation in derivations.values():
                checked_environment(derivation)
            print(
                f"Dummy {marker[-1]}: {len(identity)} roots, "
                f"{len(derivations)} marker-free derivations",
                flush=True,
            )
    if identities[0] != identities[1]:
        raise ValueError("Credential changed producer or launcher identity")
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "identities.json").write_text(
        json.dumps(identities[0], indent=2) + "\n"
    )
    print(
        "All A/B identities match; both credential forms are outside derivation metadata."
    )


class BoundaryTests(unittest.TestCase):
    def test_environment_is_private_before_metadata(self):
        with patch.dict(
            os.environ,
            SYMBOLICA_LICENSE="fixture-explicit",
            SYMBOLICA_LICENSE_SIGNED="fixture-signed",
        ):
            env = probe_environment(Path("/private-fixture"))
        self.assertNotIn("SYMBOLICA_LICENSE", env)
        self.assertNotIn("SYMBOLICA_LICENSE_SIGNED", env)
        self.assertEqual(env["HOME"], "/private-fixture")
        self.assertEqual(env["XDG_CONFIG_HOME"], "/private-fixture/config")

    def test_derivation_license_fields_and_impure_requests_are_rejected(self):
        for variable in LICENSE_VARIABLES:
            for container in ("env", "structuredAttrs"):
                for fields in (
                    {variable: "public-fixture"},
                    {"impureEnvVars": f"NIX_PATH {variable}"},
                    {"impureEnvVars": ["NIX_PATH", variable]},
                ):
                    with self.subTest(
                        variable=variable, container=container, fields=fields
                    ):
                        with self.assertRaises(ValueError):
                            checked_environment({container: fields})
            with self.assertRaises(ValueError):
                checked_environment(
                    {
                        "env": {"impureEnvVars": variable},
                        "structuredAttrs": {"impureEnvVars": []},
                    }
                )

    def test_runtime_script_references_and_unrelated_impure_variables_are_allowed(self):
        checked_environment(
            {
                "env": {
                    "text": 'test -n "$SYMBOLICA_LICENSE"',
                    "impureEnvVars": ["NIX_PATH"],
                }
            }
        )

    def test_flake_source_is_required_before_any_command(self):
        with (
            patch.object(sys, "argv", ["probe", "--report-dir", "/unused"]),
            patch.object(sys, "stderr", new_callable=io.StringIO) as diagnostic,
            patch(__name__ + ".safe_flake_source") as select,
            patch(__name__ + ".run") as execute,
            self.assertRaises(SystemExit) as failed,
        ):
            main()
        self.assertEqual(failed.exception.code, 2)
        self.assertIn("--flake-source", diagnostic.getvalue())
        select.assert_not_called()
        execute.assert_not_called()

    def test_raw_checkout_is_rejected_before_nix(self):
        with (
            patch.object(
                sys, "argv", ["probe", "--report-dir", "/unused", "--flake-source", "."]
            ),
            patch.object(sys, "stderr", new_callable=io.StringIO) as diagnostic,
            patch(__name__ + ".run") as execute,
            self.assertRaises(SystemExit) as failed,
        ):
            main()
        self.assertEqual(failed.exception.code, 2)
        self.assertIn("Nix-store snapshot", diagnostic.getvalue())
        execute.assert_not_called()

    def test_ab_identity_and_metadata_gates(self):
        identity = {
            name: {"drv": "/nix/store/fixture.drv", "out": "/nix/store/out"}
            for name in REQUIRED_ROOTS
        }
        for failure, message in (
            (None, None),
            ("empty_identity", "No producer or launcher identities"),
            *(
                (name, "Missing required producer or launcher identities: " + name)
                for name in REQUIRED_ROOTS
            ),
            ("marker", "Dummy marker in derivation metadata"),
            ("empty_metadata", "No derivation metadata"),
            ("license_field", "License derivation field"),
            ("drv", "Credential changed producer or launcher identity"),
            ("out", "Credential changed producer or launcher identity"),
        ):
            calls = []

            def execute(command, env):
                calls.append((command, env))
                if command[1] == "flake":
                    for variable in LICENSE_VARIABLES:
                        self.assertNotIn(variable, env)
                    return b'{"path": "/nix/store/fixture-source"}'
                marker = env["SYMBOLICA_LICENSE"]
                self.assertTrue(marker.startswith("PR3_BOUNDARY_DUMMY_"))
                if marker.endswith("_A"):
                    self.assertNotIn("SYMBOLICA_LICENSE_SIGNED", env)
                else:
                    self.assertTrue(marker.endswith("_B"))
                    self.assertEqual(env["SYMBOLICA_LICENSE_SIGNED"], marker)
                if command[1] == "eval":
                    if failure == "empty_identity":
                        return b"{}"
                    if failure in REQUIRED_ROOTS:
                        return json.dumps(
                            {
                                name: entry
                                for name, entry in identity.items()
                                if name != failure
                            }
                        ).encode()
                    if failure in ("drv", "out") and len(calls) == 4:
                        changed = {
                            name: dict(entry) for name, entry in identity.items()
                        }
                        changed[REQUIRED_ROOTS[0]][failure] = "/nix/store/changed"
                        return json.dumps(changed).encode()
                    return json.dumps(identity).encode()
                self.assertEqual(command[1:4], ["derivation", "show", "--recursive"])
                fields = {"text": 'test -n "$SYMBOLICA_LICENSE"'}
                if failure == "marker":
                    fields["text"] = env[LICENSE_VARIABLES[0]]
                elif failure == "license_field":
                    fields[LICENSE_VARIABLES[0]] = ""
                derivations = {"/nix/store/fixture.drv": {"env": fields}}
                if failure == "empty_metadata":
                    derivations = {}
                if len(calls) == 5:
                    derivations = {"version": 1, "derivations": derivations}
                return json.dumps(derivations).encode()

            with (
                self.subTest(failure=failure),
                tempfile.TemporaryDirectory() as tmp,
                patch(__name__ + ".run", side_effect=execute),
                patch.object(sys, "stdout", new_callable=io.StringIO),
            ):
                report = Path(tmp) / "report"
                if message:
                    with self.assertRaisesRegex(ValueError, message):
                        check_boundary("path:/nix/store/fixture-source", report)
                    self.assertFalse(report.exists())
                else:
                    check_boundary("path:/nix/store/fixture-source", report)
                    self.assertEqual(len(calls), 5)
                    self.assertNotEqual(
                        calls[1][1]["SYMBOLICA_LICENSE"],
                        calls[3][1]["SYMBOLICA_LICENSE"],
                    )
                    self.assertEqual(
                        json.loads((report / "identities.json").read_text()), identity
                    )

    def test_command_failures_include_redacted_command_and_stderr(self):
        for command in (
            ["nix", "flake", "metadata"],
            ["nix", "eval"],
            ["nix", "derivation", "show"],
        ):
            error = subprocess.CalledProcessError(
                17,
                [*command, "fixture-explicit"],
                stderr=b"fixture diagnostic: fixture-explicit fixture-signed",
            )

            def execute(actual, env):
                if actual[: len(command)] == command:
                    raise error
                if actual[1] == "flake":
                    return b'{"path": "/nix/store/fixture-source"}'
                self.assertEqual(actual[1], "eval")
                return json.dumps(
                    {name: {"drv": "/nix/store/fixture.drv"} for name in REQUIRED_ROOTS}
                ).encode()

            with (
                self.subTest(command=command),
                patch.dict(
                    os.environ,
                    SYMBOLICA_LICENSE="fixture-explicit",
                    SYMBOLICA_LICENSE_SIGNED="fixture-signed",
                ),
                patch.object(
                    sys,
                    "argv",
                    ["probe", "--report-dir", "/unused", "--flake-source", "/source"],
                ),
                patch(__name__ + ".safe_flake_source", return_value="path:/source"),
                patch(__name__ + ".run", side_effect=execute),
                patch.object(sys, "stderr", new_callable=io.StringIO) as diagnostic,
                self.assertRaises(SystemExit) as failed,
            ):
                main()
            self.assertEqual(failed.exception.code, 2)
            message = diagnostic.getvalue()
            self.assertIn(" ".join(command), message)
            self.assertIn("exited with status 17", message)
            self.assertIn("fixture diagnostic:", message)
            self.assertIn("[redacted]", message)
            self.assertNotIn("fixture-explicit", message)
            self.assertNotIn("fixture-signed", message)


def main():
    if sys.argv[1:] == ["--self-test"]:
        unittest.main(argv=[sys.argv[0]])
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument(
        "--flake-source",
        required=True,
        help="Existing immutable Nix-store snapshot of tracked source",
    )
    args = parser.parse_args()
    try:
        source_ref = safe_flake_source(
            Path.cwd(), args.flake_source, without_credentials()
        )
        check_boundary(source_ref, args.report_dir)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.error("License boundary probe failed: " + command_error(error))


if __name__ == "__main__":
    main()
