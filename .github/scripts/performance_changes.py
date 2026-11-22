"""Select code changes that need a fresh performance comparison."""

import argparse
import fnmatch
import os
import re
import subprocess

# fnmatch's '*' includes '/', so nested crate sources and examples are covered.
SHARED_INPUTS = (
    "Cargo.toml",
    "Cargo.lock",
    "build.rs",
    "rust-toolchain*",
    ".cargo/*",
    "crates/*.rs",
    "crates/*/Cargo.toml",
    "crates/*/form_src/*",
    "crates/*/templates/*",
    "crates/*/assets/*",
    "assets/*",
    "flake.nix",
    "flake.lock",
    "nix/*",
    ".github/scripts/performance_changes.py",
    ".github/workflows/performance-main.yml",
    ".github/workflows/nixci-readiness.yml",
)
CLI_INPUTS = SHARED_INPUTS + (
    "tests/Cargo.toml",
    "benchmarks/performance/*",
    ".github/scripts/run_performance.py",
    ".github/scripts/performance_metadata.py",
    ".github/workflows/performance.yml",
)
NOTEBOOK_INPUTS = SHARED_INPUTS + (
    "examples/*.py",
    "examples/notebooks/*",
    "crates/*/examples/*",
    "crates/*.py",
    "crates/*/src/*",
    "crates/*/data/*",
    "crates/linnet-render-py/*",
    "crates/*/typst/*",
    "vendor/*",
    "docs/api/python/spynso3.pyi",
    "benchmarks/performance/gammabench.rev",
    "scripts/measure-notebook-cells.py",
    "scripts/test_measure_notebook_cells.py",
    ".github/scripts/notebook_cell_gate.py",
    ".github/actions/nix-ci-cache/*",
    ".github/workflows/notebook-performance.yml",
)


def requires_measurement(paths, patterns=CLI_INPUTS):
    return any(
        fnmatch.fnmatchcase(path, pattern) for path in paths for pattern in patterns
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline")
    parser.add_argument("candidate")
    args = parser.parse_args()
    if any(not re.fullmatch(r"[0-9a-f]{40}", rev) for rev in vars(args).values()):
        parser.error("Use pinned complete commit hashes")
    paths = (
        subprocess.run(
            [
                "git",
                "diff",
                "--name-only",
                "--no-renames",
                "-z",
                args.baseline,
                args.candidate,
            ],
            check=True,
            capture_output=True,
        )
        .stdout.decode()
        .split("\0")
    )
    required = str(requires_measurement(paths)).lower()
    notebooks = str(requires_measurement(paths, NOTEBOOK_INPUTS)).lower()
    with open(os.environ["GITHUB_OUTPUT"], "a") as output:
        output.write(f"run-performance={required}\n")
        output.write(f"run-notebooks={notebooks}\n")
    print(f"Fresh performance comparison required: {required}")
    print(f"Notebook cell comparison required: {notebooks}")


if __name__ == "__main__":
    main()
