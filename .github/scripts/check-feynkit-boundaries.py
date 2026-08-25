#!/usr/bin/env python3
"""Reject dependencies from standalone FeynKit crates into GammaLoop."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
FORBIDDEN_PACKAGES = {"gammaloop-api", "gammalooprs"}
FORBIDDEN_SOURCE_REFERENCES = (
    "gammaloop-api",
    "gammaloop_api",
    "gammalooprs",
)


def cargo_metadata() -> dict[str, object]:
    command = [
        "cargo",
        "metadata",
        "--all-features",
        "--format-version=1",
        "--locked",
    ]
    result = subprocess.run(
        command,
        cwd=REPOSITORY_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode:
        sys.stderr.write(result.stderr)
        raise SystemExit(result.returncode)
    return json.loads(result.stdout)


def is_feynkit(package: dict[str, object]) -> bool:
    name = str(package["name"])
    return name == "feynkit" or name.startswith("feynkit-")


def dependency_violations(metadata: dict[str, object]) -> list[str]:
    packages = {package["id"]: package for package in metadata["packages"]}
    nodes = {node["id"]: node for node in metadata["resolve"]["nodes"]}
    roots = [package for package in packages.values() if is_feynkit(package)]
    if not roots:
        return ["cargo metadata did not contain any FeynKit packages"]

    violations = []
    for root in roots:
        pending = [(root["id"], [root["name"]])]
        visited = set()
        while pending:
            package_id, path = pending.pop()
            if package_id in visited:
                continue
            visited.add(package_id)
            for dependency_id in nodes[package_id]["dependencies"]:
                dependency = packages[dependency_id]
                dependency_path = [*path, dependency["name"]]
                if dependency["name"] in FORBIDDEN_PACKAGES:
                    violations.append(" -> ".join(dependency_path))
                else:
                    pending.append((dependency_id, dependency_path))
    return sorted(set(violations))


def source_violations(metadata: dict[str, object]) -> list[str]:
    violations = []
    for package in metadata["packages"]:
        if not is_feynkit(package):
            continue
        package_root = Path(package["manifest_path"]).parent
        for source in package_root.rglob("*.rs"):
            text = source.read_text(encoding="utf-8")
            for reference in FORBIDDEN_SOURCE_REFERENCES:
                if reference in text:
                    relative = source.relative_to(REPOSITORY_ROOT)
                    violations.append(f"{relative}: references {reference}")
    return sorted(violations)


def main() -> None:
    metadata = cargo_metadata()
    dependency_errors = dependency_violations(metadata)
    source_errors = source_violations(metadata)
    if dependency_errors or source_errors:
        print(
            "FeynKit must remain independent of GammaLoop application crates.",
            file=sys.stderr,
        )
        for violation in [*dependency_errors, *source_errors]:
            print(f"  - {violation}", file=sys.stderr)
        raise SystemExit(1)
    print("FeynKit dependency boundary is valid.")


if __name__ == "__main__":
    main()
