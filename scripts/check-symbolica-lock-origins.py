#!/usr/bin/env python3
"""Validate each workspace's declared Symbolica family against its own Cargo lock."""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from pathlib import Path
from urllib.parse import parse_qsl, urlsplit, urlunsplit

FAMILY = ("symbolica", "numerica", "graphica")
CRATES_IO = "registry+https://github.com/rust-lang/crates.io-index"
COMPONENT = r"(?:0|[1-9][0-9]*)"
VERSION = rf"{COMPONENT}\.{COMPONENT}\.{COMPONENT}"


class InvalidFamily(ValueError):
    """The lock does not implement the workspace's declared family."""


def version_tuple(value: str) -> tuple[int, int, int]:
    if not re.fullmatch(VERSION, value):
        raise InvalidFamily(f"unsupported version syntax {value!r}; expected x.y.z")
    return tuple(map(int, value.split(".")))


def permits(requirement: str, resolved: str) -> bool:
    """The exact/bare/caret full versions used by these workspace manifests."""
    match = re.fullmatch(rf"\s*([=^]?)\s*({VERSION})\s*", requirement)
    if match is None:
        raise InvalidFamily(f"unsupported Cargo version requirement {requirement!r}")
    operator, minimum = match.groups()
    lower, actual = version_tuple(minimum), version_tuple(resolved)
    if operator == "=":
        return actual == lower
    major, minor, patch = lower
    upper = (
        (major + 1, 0, 0)
        if major
        else (0, minor + 1, 0)
        if minor
        else (0, 0, patch + 1)
    )
    return lower <= actual < upper


def git_origin(spec: dict) -> tuple[str, str]:
    if not isinstance(spec, dict):
        raise InvalidFamily(f"unsupported Git patch declaration: {spec!r}")
    unsupported = set(spec) - {"git", "rev", "version"}
    if unsupported or not isinstance(spec.get("git"), str):
        raise InvalidFamily(f"unsupported Git patch declaration: {spec!r}")
    revision = spec.get("rev", "")
    if not isinstance(revision, str) or not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise InvalidFamily(
            "Symbolica Git patches require an immutable full 40-hex rev"
        )
    parsed = urlsplit(spec["git"])
    if parsed.scheme != "https" or not parsed.netloc or parsed.query or parsed.fragment:
        raise InvalidFamily(f"unsupported Git URL {spec['git']!r}")
    return spec["git"].rstrip("/").removesuffix(".git"), revision


def locked_origin(source: str) -> tuple[str, str] | str:
    if source == CRATES_IO:
        return source
    if not isinstance(source, str) or not source.startswith("git+"):
        raise InvalidFamily(f"unsupported lock source {source!r}")
    parsed = urlsplit(source[4:])
    selectors = parse_qsl(parsed.query, keep_blank_values=True)
    if len(selectors) != 1 or selectors[0][0] != "rev":
        raise InvalidFamily(
            f"lock Git source needs exactly one rev selector: {source!r}"
        )
    revision = selectors[0][1]
    if not re.fullmatch(r"[0-9a-f]{40}", revision) or parsed.fragment != revision:
        raise InvalidFamily(f"lock selector and resolved commit differ: {source!r}")
    url = urlunsplit((parsed.scheme, parsed.netloc, parsed.path, "", ""))
    return git_origin({"git": url, "rev": revision})


def validate(manifest: dict, lock: dict) -> tuple[str, str]:
    patches = manifest.get("patch", {}).get("crates-io", {})
    selected = [name for name in FAMILY if name in patches]
    if selected and len(selected) != len(FAMILY):
        raise InvalidFamily(
            "crates-io patches must declare the complete Symbolica family"
        )
    expected = CRATES_IO
    if selected:
        origins = {git_origin(patches[name]) for name in FAMILY}
        if len(origins) != 1:
            raise InvalidFamily("Symbolica family patches have different Git origins")
        expected = origins.pop()
    packages = {}
    for name in FAMILY:
        entries = [
            entry for entry in lock.get("package", []) if entry.get("name") == name
        ]
        if len(entries) != 1:
            raise InvalidFamily(
                f"Cargo.lock must contain exactly one {name}; found {len(entries)}"
            )
        packages[name] = entries[0]
    versions = {entry["version"] for entry in packages.values()}
    if len(versions) != 1:
        raise InvalidFamily("Symbolica family resolves multiple versions")
    version = versions.pop()
    version_tuple(version)
    for name, entry in packages.items():
        if locked_origin(entry.get("source")) != expected:
            raise InvalidFamily(
                f"{name} lock source differs from this workspace's declaration"
            )
    # Crates.io patches cannot authorize replacing a direct Git/path dependency.
    # Unsupported declarations fail explicitly rather than approximating Cargo.
    for owner in (manifest, manifest.get("workspace", {})):
        for table in ("dependencies", "build-dependencies", "dev-dependencies"):
            for name, declared in owner.get(table, {}).items():
                actual_name = (
                    declared.get("package", name)
                    if isinstance(declared, dict)
                    else name
                )
                if actual_name not in FAMILY:
                    continue
                spec = {"version": declared} if isinstance(declared, str) else declared
                if not isinstance(spec, dict) or set(spec) & {
                    "git",
                    "path",
                    "registry",
                    "workspace",
                }:
                    raise InvalidFamily(
                        f"unsupported direct {name} dependency declaration: {declared!r}"
                    )
                requirement = spec.get("version")
                if not isinstance(requirement, str) or not permits(
                    requirement, version
                ):
                    raise InvalidFamily(
                        f"{name} resolved {version} violates declared version {requirement!r}"
                    )
    for name in selected:
        requirement = patches[name].get("version")
        if requirement is not None and not permits(requirement, version):
            raise InvalidFamily(
                f"{name} patch version does not permit resolved {version}"
            )
    return version, str(expected)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspaces", nargs="+", type=Path)
    args = parser.parse_args()
    for workspace in args.workspaces:
        try:
            with (workspace / "Cargo.toml").open("rb") as stream:
                manifest = tomllib.load(stream)
            with (workspace / "Cargo.lock").open("rb") as stream:
                lock = tomllib.load(stream)
            version, origin = validate(manifest, lock)
        except (OSError, ValueError, KeyError, TypeError) as error:
            print(f"{workspace}: {error}", file=sys.stderr)
            return 1
        print(f"{workspace}: coherent Symbolica family {version} from {origin}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
