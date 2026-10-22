#!/usr/bin/env python3
"""Probe only new pure inputs with dummy credentials, never a user's key."""

import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path


def run(command, env):
    return subprocess.run(command, env=env, check=True, capture_output=True).stdout


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--outputs-file", type=Path)
    parser.add_argument(
        "--flake-source", help="Existing immutable flake snapshot to validate"
    )
    args = parser.parse_args()
    args.report_dir.mkdir(parents=True, exist_ok=True)
    root = Path.cwd()
    # Freeze the source once: independent contributors may edit this workspace
    # while the four intentionally impure probes run.
    metadata = json.loads(
        run(
            [
                "nix",
                "flake",
                "metadata",
                "--json",
                "--no-write-lock-file",
                args.flake_source or "path:" + str(root),
            ],
            os.environ,
        )
    )
    source = metadata["path"]
    expr = """let
      f = builtins.getFlake @FLAKE@;
      p = f.packages.x86_64-linux;
      c = f.checks.x86_64-linux;
      selected = builtins.filter (name:
        builtins.match "nix-ci-check-.*" name != null
        || builtins.elem name ["alphal00p-docs-build" "ci-runtime-environment"
          "ci-doctest-artifacts" "alphal00p-docs-cargo-artifacts" "cargoCheckArtifacts"])
        (builtins.attrNames p);
      archives = builtins.filter (name: builtins.match "gammaloop-nextest-binaries.*" name != null)
        (builtins.attrNames c);
      identity = d: { drv = d.drvPath; out = d.outPath; };
    in assert !(p ? alphal00p-docs-pages || p ? alphal00p-docs-snapshot-fixture);
    assert !(c ? gammaloop-doctest || c ? gammaloop-nextest || c ? alphal00p-docs);
    builtins.listToAttrs
      ((map (name: {inherit name; value = identity p.${name};}) selected)
      ++ (map (name: {inherit name; value = identity c.${name};}) archives))
    """.replace("@FLAKE@", json.dumps("path:" + source))
    markers = ("PR3_BOUNDARY_DUMMY_A", "PR3_BOUNDARY_DUMMY_B")
    identities = []
    with tempfile.TemporaryDirectory(prefix="gammaloop-license-probe-") as home:
        for variable in ("SYMBOLICA_LICENSE", "SYMBOLICA_LICENSE_SIGNED"):
            for marker in markers:
                env = dict(os.environ, HOME=home, XDG_CONFIG_HOME=home + "/config")
                for key in ("SYMBOLICA_LICENSE", "SYMBOLICA_LICENSE_SIGNED"):
                    env.pop(key, None)
                env[variable] = marker
                identity = json.loads(
                    run(["nix", "eval", "--impure", "--json", "--expr", expr], env)
                )
                identities.append(identity)
                drvs = sorted({entry["drv"] for entry in identity.values()})
                serialized = run(
                    ["nix", "derivation", "show", "--recursive", *drvs], env
                )
                assert all(value.encode() not in serialized for value in markers), (
                    "Dummy marker in derivation metadata"
                )
                data = json.loads(serialized)
                # Nix 2.34 wraps derivations in {version, derivations}; older
                # versions return a mapping directly.
                derivations = data.get("derivations", data)
                assert isinstance(derivations, dict) and derivations
                scripts_checked = 0
                for drv in derivations.values():
                    fields = drv.get("env", {}) | drv.get("structuredAttrs", {})
                    assert not any(
                        key in fields
                        for key in ("SYMBOLICA_LICENSE", "SYMBOLICA_LICENSE_SIGNED")
                    ), "License derivation field"
                    text = fields.get("text", "")
                    if (
                        "runtime_name=" in text
                        or "ci-runtime-environment/environment" in text
                    ):
                        subprocess.run(
                            ["bash", "-n"],
                            input=text.encode(),
                            check=True,
                            capture_output=True,
                        )
                        scripts_checked += 1
                assert scripts_checked >= 20, (
                    "Generated runtime scripts were not inspected"
                )
                receipt = {
                    "variable": variable,
                    "source": source,
                    "source_nar_hash": metadata["locked"]["narHash"],
                    "schema_keys": list(data),
                    "generated_scripts_syntax_checked": scripts_checked,
                    "serialized_sha256": hashlib.sha256(serialized).hexdigest(),
                    "dummy": marker,
                    "roots": len(identity),
                    "derivations": len(derivations),
                    "marker_free": True,
                    "license_environment_fields": False,
                }
                (args.report_dir / f"{variable}-{marker[-1]}.json").write_text(
                    json.dumps(receipt, indent=2) + "\n"
                )
                print(
                    f"{variable} dummy {marker[-1]}: {len(identity)} stable roots, {len(derivations)} marker-free derivations",
                    flush=True,
                )
    assert all(identity == identities[0] for identity in identities), (
        "Credential changed producer or launcher identity"
    )
    (args.report_dir / "identities.json").write_text(
        json.dumps(identities[0], indent=2) + "\n"
    )
    if args.outputs_file:
        outputs = args.outputs_file.read_text().splitlines()
        clean_env = {
            key: value
            for key, value in os.environ.items()
            if key not in ("SYMBOLICA_LICENSE", "SYMBOLICA_LICENSE_SIGNED")
        }
        paths = (
            run(["nix-store", "--query", "--requisites", *outputs], clean_env)
            .decode()
            .splitlines()
        )
        scanned_files = scanned_bytes = 0
        for path in map(Path, paths):
            candidates = path.rglob("*") if path.is_dir() else [path]
            for candidate in candidates:
                if candidate.is_symlink() or not candidate.is_file():
                    continue
                with candidate.open("rb") as handle:
                    tail = b""
                    while block := handle.read(1024 * 1024):
                        assert all(
                            marker.encode() not in tail + block for marker in markers
                        ), "Dummy marker in selected output closure"
                        scanned_bytes += len(block)
                        tail = block[-64:]
                scanned_files += 1
        (args.report_dir / "output-scan.json").write_text(
            json.dumps(
                {
                    "outputs": outputs,
                    "closure_paths": len(paths),
                    "files": scanned_files,
                    "bytes": scanned_bytes,
                    "marker_free": True,
                },
                indent=2,
            )
            + "\n"
        )
    print(
        "All A/B identities match; both credential forms are outside derivation metadata."
    )


if __name__ == "__main__":
    main()
