#!/usr/bin/env python3
"""Exercise the packaged lifecycle with value-free fixture receipts."""

import argparse
import json
import os
import subprocess
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("launcher", type=Path)
    parser.add_argument("--report-dir", type=Path, required=True)
    args = parser.parse_args()
    args.report_dir.mkdir(parents=True, exist_ok=True)
    launcher = args.launcher.resolve()
    markers = ("PR3_FIXTURE_EXPLICIT", "PR3_FIXTURE_SIGNED", "PR3_FIXTURE_FILE")
    receipts = []
    with tempfile.TemporaryDirectory(prefix="gammaloop-runtime-contract-") as tmp:
        root = Path(tmp)
        config = root / "config"
        (config / "symbolica").mkdir(parents=True)
        license_file = config / "symbolica/license"
        license_file.write_text(markers[2] + "\n")
        counter = root / "counter"
        for name, explicit, signed, expected, status in (
            ("explicit", markers[0], markers[1], markers[0], 0),
            ("cached-second-run", markers[0], markers[1], markers[0], 0),
            ("signed", "", markers[1], markers[1], 0),
            ("file", "", "", markers[2], 0),
            ("failure", markers[0], "", markers[0], 37),
            ("missing", "", "", "", 1),
        ):
            if name == "missing":
                license_file.unlink()
            env = dict(
                os.environ,
                HOME=str(root / "caller-home"),
                XDG_CONFIG_HOME=str(config),
                CI_RUNTIME_COUNTER=str(counter),
                CI_RUNTIME_EXPECTED=expected,
                CI_RUNTIME_EXIT=str(status),
                CI_RUNTIME_OUTPUT_EXPECTED="",
                NIX_BUILD_CORES="2",
            )
            for key, value in (
                ("SYMBOLICA_LICENSE", explicit),
                ("SYMBOLICA_LICENSE_SIGNED", signed),
            ):
                env.pop(key, None)
                if value:
                    env[key] = value
            report = args.report_dir / name
            result = subprocess.run(
                [str(launcher), "--report-dir", str(report)],
                env=env,
                capture_output=True,
                check=False,
            )
            assert result.returncode == status, (
                f"{name}: exit {result.returncode}, expected {status}"
            )
            assert all(
                marker.encode() not in result.stdout + result.stderr
                for marker in markers
            ), "Credential in output"
            if name != "missing":
                assert (report / "fixture.txt").read_text() == "retained\n"
                assert not Path(
                    (report / "scratch-path.txt").read_text().strip()
                ).exists(), "Scratch directory survived cleanup"
                assert all(
                    marker.encode() not in (report / "runtime.log").read_bytes()
                    for marker in markers
                )
            else:
                assert (
                    b"Missing SYMBOLICA_LICENSE environment variable" in result.stderr
                )
            receipts.append(
                {
                    "case": name,
                    "exit": result.returncode,
                    "report_retained": name != "missing",
                }
            )
        external = root / "external"
        external.mkdir()
        (external / "keep.txt").write_text("user data\n")
        link = root / "output-link"
        link.symlink_to(external, target_is_directory=True)
        for name, output in (
            ("output-symlink", str(link)),
            ("output-parent-component", str(root / "unused" / ".." / "output-link")),
        ):
            env.update(
                SYMBOLICA_LICENSE=markers[0],
                CI_RUNTIME_EXPECTED=markers[0],
                CI_RUNTIME_EXIT="0",
                CI_RUNTIME_OUTPUT_EXPECTED=output,
            )
            result = subprocess.run(
                [
                    str(launcher),
                    "--output",
                    output,
                    "--report-dir",
                    str(args.report_dir / name),
                ],
                env=env,
                capture_output=True,
                check=False,
            )
            assert result.returncode == 0, (
                "Lifecycle canonicalized output before its owner"
            )
            assert (external / "keep.txt").read_text() == "user data\n"
            assert link.is_symlink()
            receipts.append({"case": name, "exit": 0, "lexical_output_preserved": True})
        assert counter.read_text().splitlines() == ["executed"] * 7, (
            "Cached launcher skipped a child execution"
        )
    (args.report_dir / "receipt.json").write_text(json.dumps(receipts, indent=2) + "\n")
    print(
        "Runtime contract passed: explicit/signed/file precedence, missing key, exit 37, reports and repeated cached execution."
    )


if __name__ == "__main__":
    main()
