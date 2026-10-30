"""Pinned report/static checks, without a license in evaluations or builds."""

import hashlib
import importlib.util
import os
from pathlib import Path
import subprocess

SPEC = importlib.util.spec_from_file_location(
    "certificate", Path(__file__).with_name("validate.py")
)
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)
OUT = Path(
    os.environ.get(
        "FAMILY_ROUND2_OUT",
        str(V.WORKSPACE / "target/oracle-numerator-family-round2-final"),
    )
)


def main():
    env = dict(os.environ)
    env.pop("SYMBOLICA_LICENSE", None)
    env.pop("SYMBOLICA_LICENSE_SIGNED", None)
    processes = []
    for procdir in Path("/proc").iterdir():
        if not procdir.name.isdecimal():
            continue
        try:
            argv = (procdir / "cmdline").read_bytes().decode().split("\0")
        except (OSError, UnicodeDecodeError):
            continue
        if not argv or not (
            argv[0].endswith("python3") or argv[0].endswith("oracle_family_exact")
        ):
            continue
        if any(
            str(V.WORKSPACE) in a or str(V.HIST_WORKSPACE) in a for a in argv
        ) and any(
            a.endswith("/validate.py")
            or a.endswith("/controls.py")
            or a.endswith("/replay.py")
            or a.endswith("/oracle_family_exact")
            for a in argv
        ):
            processes.append({"pid": int(procdir.name), "argv": argv})
    V.require(not processes, "process-guard", str(processes))
    V.write(OUT / "process-guard.json", {"owned_diagnostic_processes": processes})
    commands = [
        [
            "/nix/store/2szzba1xi60y8i3qcyxwvgq2i6w9iy5s-typst-0.15.1/bin/typst",
            "compile",
            "--root",
            ".",
            "--creation-timestamp",
            "0",
            "docs/architecture/generation-oracle-family-diagnosis-round2.typ",
            str(OUT / "report.pdf"),
        ],
        [
            "/nix/store/0zh3g7fz7h3nh34cnzda1dbvg0cgwfsy-ruff-0.16.6/bin/ruff",
            "format",
            "--check",
            *[str(p) for p in sorted(Path(__file__).parent.glob("*.py"))],
        ],
        [str(V.WORKSPACE / "target/debug/alphal00p-docs-builder"), "check"],
    ]
    rows = []
    for index, argv in enumerate(commands):
        proc = subprocess.run(
            argv, cwd=V.WORKSPACE, env=env, capture_output=True, text=True
        )
        log = OUT / f"static-check-{index}.log"
        log.write_text(proc.stdout + proc.stderr)
        row = {
            "argv": argv,
            "exit": proc.returncode,
            "log": str(log),
            "license_variables_absent": True,
        }
        rows.append(row)
        if index < 2:
            V.require(proc.returncode == 0, "static-check", str(argv))
        else:
            V.require(
                proc.returncode == 1
                and "shrinking legacy inventory" in proc.stderr
                and proc.stderr.count(
                    '"docs/architecture/evidence/generation-integration/mark-geometry/'
                )
                == 8,
                "inherited-docs",
                proc.stderr,
            )
            row["classification"] = (
                "exactly eight inherited mark-geometry HTML inventory paths"
            )
    metadata = V.load(
        V.WORKSPACE / "target/oracle-numerator-family-round2/cargo-metadata.json"
    )
    generator = next(
        p for p in metadata["packages"] if p["name"] == "feynkit-generator"
    )
    V.require(
        not any(t["name"].startswith("oracle_family_") for t in generator["targets"]),
        "cargo-targets",
        "diagnostic example remains discovered",
    )
    V.write(
        OUT / "static-checks.json",
        {
            "commands": rows,
            "ordinary_diagnostic_cargo_targets": [],
            "pdf_sha256": hashlib.sha256((OUT / "report.pdf").read_bytes()).hexdigest(),
            "metadata_path": str(
                V.WORKSPACE
                / "target/oracle-numerator-family-round2/cargo-metadata.json"
            ),
        },
    )


if __name__ == "__main__":
    main()
