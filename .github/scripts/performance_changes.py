"""Select code changes that need a fresh performance comparison."""

import argparse
import fnmatch
import os
import re
import subprocess

PATTERNS = (
    "Cargo.toml",
    "Cargo.lock",
    "build.rs",
    "rust-toolchain*",
    ".cargo/*",
    "crates/*.rs",
    "crates/*/Cargo.toml",
    "tests/Cargo.toml",
    "crates/*/form_src/*",
    "crates/*/templates/*",
    "crates/*/assets/*",
    "assets/*",
    "flake.nix",
    "flake.lock",
    "nix/*",
    "benchmarks/performance/*",
    ".github/scripts/run_performance.py",
    ".github/scripts/performance_gate.py",
    ".github/scripts/performance_changes.py",
    ".github/workflows/performance.yml",
    ".github/workflows/nixci-readiness.yml",
)


def requires_measurement(paths):
    # fnmatch's '*' includes '/': a crate's nested Rust files and manifests
    # remain covered without depending on shell globstar settings.
    return any(
        fnmatch.fnmatchcase(path, pattern) for path in paths for pattern in PATTERNS
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
    with open(os.environ["GITHUB_OUTPUT"], "a") as output:
        output.write(f"run-performance={required}\n")
    print(f"Fresh performance comparison required: {required}")


if __name__ == "__main__":
    main()
