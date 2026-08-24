#!/usr/bin/env python3
"""Build and publish the repository's independently citable Zenodo archives.

This module deliberately isolates Zenodo's documented legacy deposition API.  It
is the only API which currently supports the fully automated draft/new-version,
community and publish workflow needed here, and keeping it behind ``ZenodoClient``
makes a future records-API migration local.
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime as dt
import gzip
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import tarfile
import time
import tomllib
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Iterable, Mapping
from pathlib import Path, PurePosixPath
from typing import Any, Protocol

PRODUCTION_API = "https://zenodo.org/api"
SANDBOX_API = "https://sandbox.zenodo.org/api"
ALLOWED_APIS = {PRODUCTION_API, SANDBOX_API}
SEMVER = re.compile(
    r"(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)"
    r"(?:-([0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*))?"
    r"(?:\+([0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*))?"
)
SAFE_SUFFIX = re.compile(r"[0-9A-Za-z.-]+")
SHA256 = re.compile(r"[0-9a-f]{64}")


class ReleaseError(RuntimeError):
    """A release invariant was not satisfied."""


class CollisionError(ReleaseError):
    """An immutable version already exists with different contents."""


@dataclasses.dataclass(frozen=True)
class Family:
    name: str
    title: str
    version_source: str
    tag: str
    marker: str
    concept_recid: str | None
    metadata_file: str
    citation_file: str
    full_repository: bool
    archive_paths: tuple[str, ...]
    crates: tuple[str, ...]
    include_sdist: bool
    legacy_concept_doi: str | None = None
    production_baseline_version: str | None = None


@dataclasses.dataclass(frozen=True)
class Registry:
    repository: str
    community: str
    grant: str
    families: dict[str, Family]


@dataclasses.dataclass(frozen=True)
class ArchiveEntry:
    name: str
    data: bytes
    mode: int = 0o644


@dataclasses.dataclass(frozen=True)
class Response:
    status: int
    headers: Mapping[str, str]
    body: bytes

    def json(self) -> Any:
        if not self.body:
            return None
        return json.loads(self.body)


class Transport(Protocol):
    def request(
        self, method: str, url: str, headers: Mapping[str, str], body: bytes | None
    ) -> Response: ...


class UrllibTransport:
    """HTTP transport which never delegates redirect handling to urllib."""

    class _NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, req, fp, code, msg, headers, newurl):
            return None

    def __init__(self, timeout: int = 60):
        self.timeout = timeout
        self.opener = urllib.request.build_opener(self._NoRedirect)

    def request(
        self, method: str, url: str, headers: Mapping[str, str], body: bytes | None
    ) -> Response:
        request = urllib.request.Request(
            url, data=body, headers=dict(headers), method=method
        )
        try:
            with self.opener.open(request, timeout=self.timeout) as response:
                return Response(
                    response.status, dict(response.headers), response.read()
                )
        except urllib.error.HTTPError as error:
            return Response(error.code, dict(error.headers), error.read())


def normalize_version(value: str) -> str:
    """Return SemVer from ``0.17.0``, ``v0.17.0`` or ``linnet-v0.17.0``."""

    candidate = value.strip()
    match = re.fullmatch(
        r"(?:(?:[A-Za-z][0-9A-Za-z-]*-v)|v)?(" + SEMVER.pattern + r")", candidate
    )
    if not match:
        raise ReleaseError(f"invalid release version: {value!r}")
    return match.group(1)


def _stable_version(value: str) -> str:
    version = normalize_version(value)
    if "-" in version or "+" in version:
        raise ReleaseError(f"production releases require stable SemVer: {value!r}")
    return version


def load_registry(path: Path) -> Registry:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReleaseError(f"cannot read registry {path}: {error}") from error
    if not isinstance(document, dict) or document.get("schema") != 1:
        raise ReleaseError(".zenodo/records.json must have schema 1")
    if set(document) != {"schema", "repository", "community", "grant", "families"}:
        raise ReleaseError(".zenodo/records.json fields differ from schema 1")
    for field in ("repository", "community", "grant"):
        if not isinstance(document.get(field), str) or not document[field]:
            raise ReleaseError(f"registry {field!r} must be a non-empty string")
    raw_families = document.get("families")
    if not isinstance(raw_families, dict) or not raw_families:
        raise ReleaseError("registry families must be a non-empty object")
    families: dict[str, Family] = {}
    required = {
        "title": str,
        "version_source": str,
        "tag": str,
        "marker": str,
        "metadata_file": str,
        "citation_file": str,
        "full_repository": bool,
        "archive_paths": list,
        "crates": list,
        "include_sdist": bool,
    }
    for name, raw in raw_families.items():
        if not isinstance(name, str) or not isinstance(raw, dict):
            raise ReleaseError("family names and definitions must be objects")
        unknown = (
            set(raw)
            - set(required)
            - {
                "concept_recid",
                "legacy_concept_doi",
                "production_baseline_version",
            }
        )
        if unknown:
            raise ReleaseError(
                f"family {name}: unknown registry fields {sorted(unknown)}"
            )
        if not re.fullmatch(r"[a-z0-9]+(?:-[a-z0-9]+)*", name):
            raise ReleaseError(f"invalid family name: {name!r}")
        for field, kind in required.items():
            if not isinstance(raw.get(field), kind):
                raise ReleaseError(f"family {name}: {field} must be {kind.__name__}")
        if raw["version_source"] != "release" and not raw["version_source"].endswith(
            "Cargo.toml"
        ):
            raise ReleaseError(f"family {name}: invalid version_source")
        for field in ("metadata_file", "citation_file"):
            _safe_repo_path(raw[field])
        if raw["version_source"] != "release":
            _safe_repo_path(raw["version_source"])
        concept = raw.get("concept_recid")
        if concept is not None and not isinstance(concept, str):
            raise ReleaseError(f"family {name}: concept_recid must be null or a string")
        for field in ("archive_paths", "crates"):
            if not all(isinstance(value, str) and value for value in raw[field]):
                raise ReleaseError(
                    f"family {name}: {field} must contain non-empty strings"
                )
        for path_value in raw["archive_paths"]:
            _safe_repo_path(path_value)
        if not all(re.fullmatch(r"[A-Za-z0-9_-]+", crate) for crate in raw["crates"]):
            raise ReleaseError(f"family {name}: invalid crate name")
        marker = raw["marker"]
        if urllib.parse.urlsplit(marker).scheme != "https":
            raise ReleaseError(f"family {name}: marker must be an HTTPS URL")
        families[name] = Family(
            name=name,
            title=raw["title"],
            version_source=raw["version_source"],
            tag=raw["tag"],
            marker=marker,
            concept_recid=concept,
            metadata_file=raw["metadata_file"],
            citation_file=raw["citation_file"],
            full_repository=raw["full_repository"],
            archive_paths=tuple(raw["archive_paths"]),
            crates=tuple(raw["crates"]),
            include_sdist=raw["include_sdist"],
            legacy_concept_doi=raw.get("legacy_concept_doi"),
            production_baseline_version=raw.get("production_baseline_version"),
        )
        for field in ("legacy_concept_doi", "production_baseline_version"):
            value = raw.get(field)
            if value is not None and not isinstance(value, str):
                raise ReleaseError(f"family {name}: {field} must be null or a string")
    return Registry(
        document["repository"], document["community"], document["grant"], families
    )


def _git(repo: Path, *arguments: str, input_bytes: bytes | None = None) -> bytes:
    result = subprocess.run(
        ["git", *arguments],
        cwd=repo,
        input=input_bytes,
        capture_output=True,
        check=False,
    )
    if result.returncode:
        raise ReleaseError(
            f"git {' '.join(arguments)} failed: {result.stderr.decode(errors='replace').strip()}"
        )
    return result.stdout


def _commit(repo: Path, ref: str) -> str:
    return _git(repo, "rev-parse", "--verify", f"{ref}^{{commit}}").decode().strip()


def _git_text(repo: Path, commit: str, path: str) -> str:
    return _git(repo, "show", f"{commit}:{path}").decode("utf-8")


def _tag_date(repo: Path, tag: str, commit: str, synthetic: bool) -> tuple[str, int]:
    if synthetic:
        timestamp = int(_git(repo, "show", "-s", "--format=%ct", commit).decode())
    else:
        kind = _git(repo, "cat-file", "-t", f"refs/tags/{tag}").decode().strip()
        if kind != "tag":
            raise ReleaseError(f"release tag {tag!r} must be annotated")
        timestamp = int(
            _git(
                repo, "for-each-ref", "--format=%(taggerdate:unix)", f"refs/tags/{tag}"
            )
            .decode()
            .strip()
        )
    date = dt.datetime.fromtimestamp(timestamp, tz=dt.UTC).date().isoformat()
    return date, timestamp


def _manifest_version(repo: Path, commit: str, source: str) -> str:
    try:
        package = tomllib.loads(_git_text(repo, commit, source))["package"]
        return _stable_version(package["version"])
    except (KeyError, tomllib.TOMLDecodeError) as error:
        raise ReleaseError(
            f"cannot read package version from {source} at {commit}"
        ) from error


def _safe_repo_path(value: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or not path.parts
        or any(part in ("", ".", "..") for part in path.parts)
    ):
        raise ReleaseError(f"unsafe repository path: {value!r}")
    return path


def _git_tree(repo: Path, commit: str) -> dict[str, tuple[bytes, int, str | None]]:
    """Read one immutable Git tree, materializing safe in-repository symlinks."""

    raw_archive = _git(repo, "archive", "--format=tar", commit)
    raw: dict[str, tuple[bytes, int, str | None]] = {}
    with tarfile.open(fileobj=io.BytesIO(raw_archive), mode="r:") as source:
        for member in source:
            name = str(_safe_repo_path(member.name.rstrip("/")))
            if member.isdir():
                continue
            if member.isfile():
                stream = source.extractfile(member)
                if stream is None:
                    raise ReleaseError(f"cannot read {name} from Git archive")
                raw[name] = (stream.read(), member.mode, None)
            elif member.issym():
                raw[name] = (b"", 0o644, member.linkname)
            else:
                raise ReleaseError(f"unsupported Git entry type: {name}")

    def resolve(name: str, trail: frozenset[str]) -> list[tuple[str, bytes, int]]:
        if name in trail:
            raise ReleaseError(f"symlink cycle at {name}")
        data, mode, link = raw[name]
        if link is None:
            return [(name, data, 0o755 if mode & 0o111 else 0o644)]
        link_path = PurePosixPath(link)
        if link_path.is_absolute():
            raise ReleaseError(f"absolute symlink target at {name}: {link}")
        stack = list(PurePosixPath(name).parent.parts)
        for part in link_path.parts:
            if part in ("", "."):
                continue
            if part == "..":
                if not stack:
                    raise ReleaseError(f"symlink escapes repository at {name}: {link}")
                stack.pop()
            else:
                stack.append(part)
        target = "/".join(stack)
        if target in raw:
            resolved = resolve(target, trail | {name})
            return [
                (name + path[len(target) :], body, bits)
                for path, body, bits in resolved
            ]
        prefix = target + "/"
        descendants = sorted(path for path in raw if path.startswith(prefix))
        if not descendants:
            raise ReleaseError(f"dangling symlink at {name}: {link}")
        materialized: list[tuple[str, bytes, int]] = []
        for descendant in descendants:
            for path, body, bits in resolve(descendant, trail | {name}):
                materialized.append((name + path[len(target) :], body, bits))
        return materialized

    tree: dict[str, tuple[bytes, int, str | None]] = {}
    for name in sorted(raw):
        for path, data, mode in resolve(name, frozenset()):
            tree[path] = (data, mode, None)
    return tree


def _selected_source_entries(
    tree: Mapping[str, tuple[bytes, int, str | None]], family: Family
) -> list[ArchiveEntry]:
    paths = tuple(
        str(_safe_repo_path(path)).rstrip("/") for path in family.archive_paths
    )
    selected = []
    for name, (data, mode, _) in tree.items():
        if family.full_repository or any(
            name == path or name.startswith(path + "/") for path in paths
        ):
            selected.append(ArchiveEntry(name, data, mode))
    if not selected:
        raise ReleaseError(
            f"family {family.name}: archive paths select no source files"
        )
    return selected


def read_release_assets(directory: Path) -> dict[str, tuple[Path, str]]:
    checksum_file = directory / "SHA256SUMS"
    try:
        lines = checksum_file.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise ReleaseError(f"cannot read {checksum_file}: {error}") from error
    expected: dict[str, str] = {}
    for line in lines:
        match = re.fullmatch(r"([0-9a-f]{64}) [ *]([^/\\]+)", line)
        if not match or match.group(2) == "SHA256SUMS":
            raise ReleaseError(f"invalid SHA256SUMS line: {line!r}")
        digest, name = match.groups()
        if name in expected:
            raise ReleaseError(f"duplicate SHA256SUMS entry: {name}")
        expected[name] = digest
    actual = {
        path.name
        for path in directory.iterdir()
        if path.is_file() and path.name != "SHA256SUMS"
    }
    if actual != set(expected):
        raise ReleaseError(
            f"release assets and SHA256SUMS differ: missing={sorted(set(expected) - actual)}, "
            f"unlisted={sorted(actual - set(expected))}"
        )
    assets: dict[str, tuple[Path, str]] = {}
    for name, expected_digest in expected.items():
        path = directory / name
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != expected_digest:
            raise ReleaseError(f"SHA-256 mismatch for {name}")
        assets[name] = (path, digest)
    return assets


def _crate_versions(
    tree: Mapping[str, tuple[bytes, int, str | None]], crates: Iterable[str]
) -> dict[str, str]:
    manifests: dict[str, str] = {}
    for path, (data, _, _) in tree.items():
        if not path.startswith("crates/") or not path.endswith("/Cargo.toml"):
            continue
        try:
            package = tomllib.loads(data.decode())["package"]
        except (KeyError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
            raise ReleaseError(f"cannot parse {path} from release tree") from error
        if isinstance(package.get("name"), str) and isinstance(
            package.get("version"), str
        ):
            manifests[package["name"]] = _stable_version(package["version"])
    requested = tuple(crates)
    missing = sorted(set(requested) - set(manifests))
    if missing:
        raise ReleaseError(f"release tree has no manifest for crates: {missing}")
    return {name: manifests[name] for name in requested}


def _asset_entries(
    family: Family,
    component_versions: Mapping[str, str],
    family_base_version: str,
    assets: Mapping[str, tuple[Path, str]],
) -> tuple[list[ArchiveEntry], list[dict[str, str]]]:
    selected: list[tuple[str, Path, str]] = []
    for crate in family.crates:
        name = f"{crate}-{component_versions[crate]}.crate"
        if name not in assets:
            raise ReleaseError(f"family {family.name}: missing canonical crate {name}")
        path, digest = assets[name]
        selected.append((name, path, digest))
    if family.include_sdist:
        candidates = [
            (name, value)
            for name, value in assets.items()
            if name.endswith(f"-{family_base_version}.tar.gz")
            and not name.endswith(".crate")
        ]
        if len(candidates) != 1:
            raise ReleaseError(
                f"family {family.name}: expected one {family_base_version} Python sdist, found "
                f"{[name for name, _ in candidates]}"
            )
        name, (path, digest) = candidates[0]
        selected.append((name, path, digest))
    entries = [ArchiveEntry(name, path.read_bytes()) for name, path, _ in selected]
    manifest = [
        {"name": name, "sha256": digest, "source": path.name}
        for name, path, digest in selected
    ]
    return entries, manifest


def render_citation(
    template: str,
    version: str,
    publication_date: str,
    version_doi: str | None = None,
    concept_doi: str | None = None,
) -> str:
    top_level = {
        line.split(":", 1)[0]
        for line in template.splitlines()
        if line and not line[0].isspace() and ":" in line
    }
    if "version" in top_level or "date-released" in top_level:
        raise ReleaseError("checked-in CFF must remain versionless")
    result = (
        template.rstrip()
        + f'\nversion: "{version}"\ndate-released: "{publication_date}"\n'
    )
    additions = []
    for doi, description in (
        (version_doi, "DOI for this software version"),
        (concept_doi, "Concept DOI for all software versions"),
    ):
        if doi and doi not in template and doi not in "".join(additions):
            additions.extend(
                [
                    "  - type: doi\n",
                    f'    value: "{doi}"\n',
                    f"    description: {description}\n",
                ]
            )
    if additions:
        if "identifiers" in top_level:
            lines = result.splitlines(keepends=True)
            index = (
                next(
                    i for i, line in enumerate(lines) if line.startswith("identifiers:")
                )
                + 1
            )
            while index < len(lines) and (
                not lines[index].strip() or lines[index][0].isspace()
            ):
                index += 1
            lines[index:index] = additions
            result = "".join(lines)
        else:
            result += "identifiers:\n" + "".join(additions)
    return result


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    ).encode()


def _with_generated_entries(
    entries: Iterable[ArchiveEntry],
    root: str,
    citation: str,
    provenance: Mapping[str, Any],
) -> list[ArchiveEntry]:
    replacements = {
        f"{root}/CITATION.cff": citation.encode(),
        f"{root}/PROVENANCE.json": _json_bytes(provenance),
    }
    kept = [
        entry
        for entry in entries
        if entry.name not in replacements and not entry.name.endswith("/SHA256SUMS")
    ]
    kept.extend(ArchiveEntry(name, data) for name, data in replacements.items())
    checksums = "".join(
        f"{hashlib.sha256(entry.data).hexdigest()}  {entry.name.removeprefix(root + '/')}\n"
        for entry in sorted(kept, key=lambda item: item.name)
    )
    kept.append(ArchiveEntry(f"{root}/SHA256SUMS", checksums.encode()))
    return kept


def write_archive(path: Path, entries: Iterable[ArchiveEntry], timestamp: int) -> str:
    ordered = sorted(entries, key=lambda entry: entry.name)
    names = [entry.name for entry in ordered]
    if len(names) != len(set(names)):
        raise ReleaseError(f"duplicate path in archive {path.name}")
    output = io.BytesIO()
    with (
        gzip.GzipFile(
            fileobj=output, mode="wb", filename="", mtime=0, compresslevel=9
        ) as zipped,
        tarfile.open(fileobj=zipped, mode="w", format=tarfile.PAX_FORMAT) as archive,
    ):
        for entry in ordered:
            _safe_repo_path(entry.name)
            info = tarfile.TarInfo(entry.name)
            info.size = len(entry.data)
            info.mode = entry.mode
            info.mtime = timestamp
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            archive.addfile(info, io.BytesIO(entry.data))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(output.getvalue())
    return hashlib.sha256(output.getvalue()).hexdigest()


def read_archive(path: Path) -> list[ArchiveEntry]:
    entries = []
    with tarfile.open(path, "r:gz") as archive:
        for member in archive:
            if not member.isfile():
                raise ReleaseError(
                    f"prepared archive contains non-file entry: {member.name}"
                )
            stream = archive.extractfile(member)
            if stream is None:
                raise ReleaseError(f"cannot read prepared archive entry: {member.name}")
            entries.append(ArchiveEntry(member.name, stream.read(), member.mode))
    return entries


def _metadata_template(
    root: Path, family: Family, registry: Registry
) -> dict[str, Any]:
    try:
        metadata = json.loads((root / family.metadata_file).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReleaseError(
            f"cannot read metadata for {family.name}: {error}"
        ) from error
    if metadata.get("title") != family.title:
        raise ReleaseError(
            f"family {family.name}: metadata title differs from registry"
        )
    for field, expected in (
        ("upload_type", "software"),
        ("access_right", "open"),
        ("license", "mit"),
    ):
        if str(metadata.get(field, "")).lower() != expected:
            raise ReleaseError(
                f"family {family.name}: metadata {field} must be {expected!r}"
            )
    if {item.get("identifier") for item in metadata.get("communities", [])} != {
        registry.community
    }:
        raise ReleaseError(
            f"family {family.name}: metadata community differs from registry"
        )
    if {item.get("id") for item in metadata.get("grants", [])} != {registry.grant}:
        raise ReleaseError(
            f"family {family.name}: metadata grant differs from registry"
        )
    for mutable in ("version", "publication_date", "prereserve_doi"):
        if mutable in metadata:
            raise ReleaseError(
                f"family {family.name}: template contains mutable field {mutable}"
            )
    return metadata


def _relation(identifier: str, relation: str) -> dict[str, str]:
    return {"identifier": identifier, "relation": relation, "resource_type": "software"}


def merge_relations(
    existing: Iterable[Mapping[str, Any]],
    marker: str,
    family_name: str,
    concept_dois: Mapping[str, str],
    version_dois: Mapping[str, str],
) -> list[dict[str, Any]]:
    relations = [dict(relation) for relation in existing]
    additions = [_relation(marker, "isAlternateIdentifier")]
    umbrella = "gammaloop-monorepo"
    if family_name == umbrella:
        additions.extend(
            _relation(doi, "hasPart")
            for name, doi in version_dois.items()
            if name != umbrella
        )
    elif umbrella in concept_dois:
        additions.append(_relation(concept_dois[umbrella], "isPartOf"))
    seen = {(item.get("identifier"), item.get("relation")) for item in relations}
    for addition in additions:
        key = (addition["identifier"], addition["relation"])
        if key not in seen:
            relations.append(addition)
            seen.add(key)
    return relations


def build_metadata(
    template: Mapping[str, Any],
    family: Family,
    metadata_version: str,
    publication_date: str,
    concept_dois: Mapping[str, str],
    version_dois: Mapping[str, str],
) -> dict[str, Any]:
    metadata = json.loads(json.dumps(template))
    metadata["version"] = metadata_version
    metadata["publication_date"] = publication_date
    metadata["related_identifiers"] = merge_relations(
        metadata.get("related_identifiers", []),
        family.marker,
        family.name,
        concept_dois,
        version_dois,
    )
    return metadata


def metadata_for_environment(
    metadata: Mapping[str, Any], production: bool
) -> dict[str, Any]:
    result = json.loads(json.dumps(metadata))
    if not production:
        result["related_identifiers"] = [
            relation
            for relation in result.get("related_identifiers", [])
            if not str(relation.get("identifier", "")).startswith("10.5281/zenodo.")
        ]
    return result


def _prepare(args: argparse.Namespace) -> None:
    root = Path(args.repository_root).resolve()
    registry_path = (root / args.registry).resolve()
    if not registry_path.is_relative_to(root):
        raise ReleaseError("registry must be inside the repository")
    registry = load_registry(registry_path)
    base_release_version = _stable_version(args.release_version)
    release_commit = _commit(root, args.release_commit)
    registry_relative = registry_path.relative_to(root).as_posix()
    if (
        _git(root, "show", f"{release_commit}:{registry_relative}")
        != registry_path.read_bytes()
    ):
        raise ReleaseError("registry must exactly match the immutable release commit")
    if args.version_suffix:
        if not SAFE_SUFFIX.fullmatch(args.version_suffix):
            raise ReleaseError("--version-suffix must match [0-9A-Za-z.-]+")
        suffix = f"+sandbox.{args.version_suffix}"
        normalize_version(base_release_version + suffix)
    else:
        suffix = ""
        if _commit(root, "HEAD") != release_commit:
            raise ReleaseError(
                "production prepare requires HEAD to equal --release-commit"
            )
        if _commit(root, f"refs/tags/{args.release_tag}") != release_commit:
            raise ReleaseError(
                "production release tag does not resolve to --release-commit"
            )
    release_date, release_timestamp = _tag_date(
        root, args.release_tag, release_commit, bool(args.version_suffix)
    )
    assets = read_release_assets(Path(args.assets_dir).resolve())
    tree_cache: dict[str, dict[str, tuple[bytes, int, str | None]]] = {}
    output = Path(args.output_dir).resolve()
    plan_path = Path(args.plan).resolve()
    if plan_path.parent != output:
        raise ReleaseError("--plan must be a file directly inside --output-dir")
    output.mkdir(parents=True, exist_ok=True)
    plan: dict[str, Any] = {
        "schema": 1,
        "repository": registry.repository,
        "registry_sha256": hashlib.sha256(registry_path.read_bytes()).hexdigest(),
        "release": {
            "version": base_release_version + suffix,
            "base_version": base_release_version,
            "tag": args.release_tag,
            "commit": release_commit,
            "publication_date": release_date,
            "source_date_epoch": release_timestamp,
            "sandbox_suffix": args.version_suffix,
        },
        "families": {},
    }
    planned_names: set[str] = set()
    for name, family in registry.families.items():
        base_version = (
            base_release_version
            if family.version_source == "release"
            else _manifest_version(root, release_commit, family.version_source)
        )
        version = base_version + suffix
        configured_tag = family.tag.format(version=base_version)
        if suffix:
            tag, commit = args.release_tag, release_commit
            publication_date, timestamp = release_date, release_timestamp
        else:
            tag = configured_tag
            commit = _commit(root, f"refs/tags/{tag}")
            if family.version_source != "release":
                tagged_version = _manifest_version(root, commit, family.version_source)
                if tagged_version != base_version:
                    raise ReleaseError(
                        f"family {name}: {tag} contains {tagged_version}, expected {base_version}"
                    )
            publication_date, timestamp = _tag_date(root, tag, commit, False)
        if commit not in tree_cache:
            tree_cache[commit] = _git_tree(root, commit)
        root_name = f"{name}-source-{version}"
        base_filename = f"prepared-{root_name}.tar.gz"
        filename = f"{root_name}.tar.gz"
        if filename in planned_names:
            raise ReleaseError(f"duplicate planned Zenodo filename: {filename}")
        planned_names.add(filename)
        source_entries = [
            ArchiveEntry(f"{root_name}/{entry.name}", entry.data, entry.mode)
            for entry in _selected_source_entries(tree_cache[commit], family)
        ]
        component_versions = _crate_versions(tree_cache[commit], family.crates)
        asset_entries, asset_manifest = _asset_entries(
            family, component_versions, base_version, assets
        )
        canonical_dir = output / "canonical"
        canonical_dir.mkdir(exist_ok=True)
        for entry in asset_entries:
            canonical_path = canonical_dir / entry.name
            if canonical_path.exists() and canonical_path.read_bytes() != entry.data:
                raise ReleaseError(f"canonical asset collision: {entry.name}")
            canonical_path.write_bytes(entry.data)
        for item in asset_manifest:
            item["prepared_path"] = f"canonical/{item['source']}"
        license_data = _git(root, "show", f"{commit}:LICENSE.md")
        source_entries = [
            entry for entry in source_entries if entry.name != f"{root_name}/LICENSE.md"
        ]
        source_entries.append(ArchiveEntry(f"{root_name}/LICENSE.md", license_data))
        citation_template = _git_text(root, release_commit, family.citation_file)
        metadata_bytes = _git(root, "show", f"{release_commit}:{family.metadata_file}")
        current_metadata = (root / family.metadata_file).read_bytes()
        current_citation = (root / family.citation_file).read_bytes()
        if (
            metadata_bytes != current_metadata
            or citation_template.encode() != current_citation
        ):
            raise ReleaseError(
                f"family {name}: metadata and citation templates must match {release_commit}"
            )
        _metadata_template(root, family, registry)
        provenance = {
            "schema": 1,
            "repository": registry.repository,
            "family": name,
            "version": version,
            "base_version": base_version,
            "tag": tag,
            "commit": commit,
            "publication_date": publication_date,
            "source_date_epoch": timestamp,
            "marker": family.marker,
            "component_versions": component_versions,
            "assets": asset_manifest,
            "concept_doi": None,
            "version_doi": None,
        }
        citation = render_citation(citation_template, version, publication_date)
        prepared_entries = _with_generated_entries(
            source_entries, root_name, citation, provenance
        )
        base_path = output / base_filename
        digest = write_archive(base_path, prepared_entries, timestamp)
        plan["families"][name] = {
            "title": family.title,
            "base_version": base_version,
            "version": version,
            "metadata_version": configured_tag + suffix,
            "tag": tag,
            "commit": commit,
            "publication_date": publication_date,
            "source_date_epoch": timestamp,
            "root": root_name,
            "base_archive": base_filename,
            "base_sha256": digest,
            "archive": filename,
            "assets": asset_manifest,
            "component_versions": component_versions,
            "metadata_sha256": hashlib.sha256(metadata_bytes).hexdigest(),
            "citation_sha256": hashlib.sha256(current_citation).hexdigest(),
            "template_commit": release_commit,
        }
    plan_path.parent.mkdir(parents=True, exist_ok=True)
    plan_path.write_bytes(_json_bytes(plan))
    print(
        json.dumps(
            {"plan": str(plan_path), "families": list(plan["families"])}, indent=2
        )
    )


def has_marker(deposition: Mapping[str, Any], marker: str) -> bool:
    return any(
        relation.get("identifier") == marker
        and str(relation.get("relation", "")).lower() == "isalternateidentifier"
        for relation in deposition.get("metadata", {}).get("related_identifiers", [])
        if isinstance(relation, dict)
    )


def select_lineage(
    depositions: Iterable[Mapping[str, Any]],
    marker: str,
    concept_recid: str | None = None,
) -> list[dict[str, Any]]:
    records = [dict(record) for record in depositions]
    marked = [record for record in records if has_marker(record, marker)]
    concepts = {str(record.get("conceptrecid")) for record in marked}
    concepts.discard("None")
    if len(concepts) > 1:
        raise CollisionError(
            f"marker {marker} is attached to multiple concepts: {sorted(concepts)}"
        )
    if concepts and concept_recid is not None and concept_recid not in concepts:
        raise CollisionError(
            f"marker {marker} resolves to concept {next(iter(concepts))}, not configured "
            f"concept {concept_recid}"
        )
    selected_concept = next(iter(concepts), concept_recid)
    if selected_concept is None:
        return []
    return sorted(
        (
            record
            for record in records
            if str(record.get("conceptrecid")) == selected_concept
        ),
        key=lambda record: int(record.get("id", 0)),
    )


def _published(record: Mapping[str, Any]) -> bool:
    return bool(record.get("submitted")) or record.get("state") == "done"


def _record_version(record: Mapping[str, Any]) -> str | None:
    value = record.get("metadata", {}).get("version")
    if not isinstance(value, str):
        return None
    try:
        return normalize_version(value)
    except ReleaseError:
        return None


def _concept_doi(record: Mapping[str, Any], prefix: str = "10.5281") -> str:
    value = record.get("conceptdoi")
    if isinstance(value, str) and value:
        return value
    recid = record.get("conceptrecid")
    if recid is None:
        raise ReleaseError("Zenodo response omitted conceptrecid")
    return f"{prefix}/zenodo.{recid}"


def _version_doi(record: Mapping[str, Any], prefix: str = "10.5281") -> str:
    for value in (
        record.get("doi"),
        record.get("metadata", {}).get("doi"),
        record.get("metadata", {}).get("prereserve_doi", {}).get("doi"),
    ):
        if isinstance(value, str) and value:
            return value
    recid = record.get("record_id") or record.get("id")
    if recid is None:
        raise ReleaseError("Zenodo response omitted reserved DOI")
    return f"{prefix}/zenodo.{recid}"


def _file_identity(record: Mapping[str, Any]) -> dict[str, str]:
    result = {}
    for item in record.get("files", []):
        name = item.get("filename") or item.get("key")
        checksum = str(item.get("checksum", ""))
        checksum = checksum.removeprefix("md5:")
        if name:
            result[str(name)] = checksum
    return result


def assert_same_version_files(
    record: Mapping[str, Any], version: str, expected: Mapping[str, str]
) -> None:
    if _record_version(record) != normalize_version(version):
        raise CollisionError("record version does not match collision check")
    actual = _file_identity(record)
    if actual != dict(expected):
        raise CollisionError(
            f"Zenodo version {version} exists with different files: expected={dict(expected)}, actual={actual}"
        )


def assert_stable_metadata(
    record: Mapping[str, Any], expected: Mapping[str, Any]
) -> None:
    actual = record.get("metadata", {})
    if _record_version(record) != normalize_version(str(expected.get("version", ""))):
        raise CollisionError("Zenodo metadata has a different normalized version")
    for field in (
        "publication_date",
        "title",
        "creators",
        "license",
        "grants",
        "communities",
    ):
        if actual.get(field) != expected.get(field):
            raise CollisionError(
                f"Zenodo metadata field {field!r} differs from the release plan"
            )
    expected_relations = {
        (item.get("identifier"), str(item.get("relation", "")).lower())
        for item in expected.get("related_identifiers", [])
    }
    actual_relations = {
        (item.get("identifier"), str(item.get("relation", "")).lower())
        for item in actual.get("related_identifiers", [])
    }
    if not expected_relations.issubset(actual_relations):
        raise CollisionError(
            f"Zenodo metadata is missing required relations: "
            f"{sorted(expected_relations - actual_relations)}"
        )


def is_legacy_linnet_match(
    family: Family, record: Mapping[str, Any], target_version: str, production: bool
) -> bool:
    return (
        family.name == "linnet"
        and production
        and family.concept_recid is not None
        and str(record.get("conceptrecid")) == family.concept_recid
        and _record_version(record) == normalize_version(target_version)
        and not has_marker(record, family.marker)
    )


def _version_order(value: str) -> tuple[int, int, int]:
    normalized = normalize_version(value)
    core = normalized.split("-", 1)[0].split("+", 1)[0]
    return tuple(int(part) for part in core.split("."))  # type: ignore[return-value]


class ZenodoClient:
    def __init__(self, token: str, api_url: str, transport: Transport | None = None):
        self.api_url = api_url.rstrip("/")
        if self.api_url not in ALLOWED_APIS:
            raise ReleaseError(f"unsupported Zenodo API URL: {api_url!r}")
        if not token:
            raise ReleaseError("ZENODO_TOKEN is required")
        self.token = token
        self.transport = transport or UrllibTransport()
        self.origin = urllib.parse.urlsplit(self.api_url)[:2]

    @property
    def production(self) -> bool:
        return self.api_url == PRODUCTION_API

    @property
    def doi_prefix(self) -> str:
        return "10.5281" if self.production else "10.5072"

    def concept_doi(self, record: Mapping[str, Any]) -> str:
        return _concept_doi(record, self.doi_prefix)

    def version_doi(self, record: Mapping[str, Any]) -> str:
        return _version_doi(record, self.doi_prefix)

    def validate_plan_environment(self, plan: Mapping[str, Any]) -> None:
        sandbox_suffix = plan["release"]["sandbox_suffix"]
        if self.production and sandbox_suffix is not None:
            raise ReleaseError("refusing to publish a sandbox plan to production")
        if not self.production and sandbox_suffix is None:
            raise ReleaseError("refusing to publish an unsuffixed plan to the sandbox")

    def _trusted_url(self, url: str) -> str:
        parsed = urllib.parse.urlsplit(url)
        if parsed.scheme != "https" or (parsed.scheme, parsed.netloc) != self.origin:
            raise ReleaseError(f"Zenodo returned an untrusted authenticated URL: {url}")
        return url

    def _request(
        self,
        method: str,
        url: str,
        payload: Any = None,
        raw: bytes | None = None,
        expected: tuple[int, ...] = (200, 201, 202),
    ) -> Any:
        url = self._trusted_url(url)
        headers = {
            "Authorization": f"Bearer {self.token}",
            "Accept": "application/json",
            "User-Agent": "gammaloop-zenodo-release/1",
        }
        body = raw
        if payload is not None:
            body = json.dumps(payload).encode()
            headers["Content-Type"] = "application/json"
        if raw is not None:
            headers["Content-Type"] = "application/octet-stream"
        attempts = 4 if method in {"GET", "PUT", "DELETE"} else 1
        for attempt in range(attempts):
            try:
                response = self.transport.request(method, url, headers, body)
            except (OSError, TimeoutError) as error:
                if attempt + 1 == attempts:
                    raise ReleaseError(
                        f"Zenodo {method} {url} failed: {error}"
                    ) from error
                time.sleep(2**attempt)
                continue
            if response.status in expected:
                return response.json()
            if (
                response.status not in {429, 500, 502, 503, 504}
                or attempt + 1 == attempts
            ):
                detail = response.body.decode(errors="replace")[:1000]
                raise ReleaseError(
                    f"Zenodo {method} {url} returned {response.status}: {detail}"
                )
            time.sleep(2**attempt)
        raise AssertionError("unreachable")

    def list_depositions(self) -> list[dict[str, Any]]:
        records = []
        page = 1
        while True:
            url = f"{self.api_url}/deposit/depositions?all_versions=true&size=100&page={page}"
            batch = self._request("GET", url)
            if not isinstance(batch, list):
                raise ReleaseError("Zenodo deposition list was not an array")
            records.extend(batch)
            if len(batch) < 100:
                return records
            page += 1

    def get(self, identifier_or_url: str | int) -> dict[str, Any]:
        url = (
            str(identifier_or_url)
            if str(identifier_or_url).startswith("https://")
            else f"{self.api_url}/deposit/depositions/{identifier_or_url}"
        )
        return self._request("GET", url)

    def create(self, metadata: Mapping[str, Any]) -> dict[str, Any]:
        return self._request(
            "POST", f"{self.api_url}/deposit/depositions", {"metadata": dict(metadata)}
        )

    def update(
        self, identifier: str | int, metadata: Mapping[str, Any]
    ) -> dict[str, Any]:
        return self._request(
            "PUT",
            f"{self.api_url}/deposit/depositions/{identifier}",
            {"metadata": dict(metadata)},
        )

    def new_version(self, latest: Mapping[str, Any]) -> dict[str, Any]:
        latest_id = latest["id"]
        existing_url = latest.get("links", {}).get("latest_draft")
        if existing_url:
            existing = self.get(self._trusted_url(existing_url))
            if not _published(existing):
                return existing
        try:
            response = self._request(
                "POST",
                f"{self.api_url}/deposit/depositions/{latest_id}/actions/newversion",
            )
        except ReleaseError:
            refreshed = self.get(latest_id)
            draft_url = refreshed.get("links", {}).get("latest_draft")
            if not draft_url:
                raise
            draft = self.get(self._trusted_url(draft_url))
            if _published(draft):
                raise
            return draft
        draft_url = response.get("links", {}).get("latest_draft")
        if not draft_url:
            refreshed = self.get(latest_id)
            draft_url = refreshed.get("links", {}).get("latest_draft")
        if not draft_url:
            raise ReleaseError("newversion did not expose links.latest_draft")
        return self.get(self._trusted_url(draft_url))

    def replace_files(
        self, draft: Mapping[str, Any], files: Mapping[str, Path]
    ) -> dict[str, Any]:
        for item in draft.get("files", []):
            url = item.get("links", {}).get("self")
            if not url and item.get("id") is not None:
                url = f"{self.api_url}/deposit/depositions/{draft['id']}/files/{item['id']}"
            if not url:
                raise ReleaseError("inherited draft file has no deletion link")
            self._request("DELETE", url, expected=(200, 204, 404))
        bucket = draft.get("links", {}).get("bucket")
        if not bucket:
            raise ReleaseError("draft has no file bucket")
        bucket = self._trusted_url(bucket).rstrip("/")
        for name, path in sorted(files.items()):
            filename = urllib.parse.quote(_safe_plan_name(name), safe="")
            self._request("PUT", f"{bucket}/{filename}", raw=path.read_bytes())
        return self.get(draft["id"])

    def publish(self, identifier: str | int) -> dict[str, Any]:
        try:
            self._request(
                "POST",
                f"{self.api_url}/deposit/depositions/{identifier}/actions/publish",
            )
        except ReleaseError:
            if not _published(self.get(identifier)):
                raise
        for attempt in range(12):
            record = self.get(identifier)
            if _published(record) and record.get("files"):
                return record
            if attempt != 11:
                time.sleep(5)
        raise ReleaseError(
            f"Zenodo deposition {identifier} was not integrated after publication"
        )


def _load_plan(
    path: Path, root: Path, registry_path: Path
) -> tuple[dict[str, Any], Registry]:
    try:
        plan = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ReleaseError(f"cannot read plan {path}: {error}") from error
    if not isinstance(plan, dict) or plan.get("schema") != 1:
        raise ReleaseError("unsupported Zenodo plan schema")
    registry = load_registry(registry_path)
    if hashlib.sha256(registry_path.read_bytes()).hexdigest() != plan.get(
        "registry_sha256"
    ):
        raise ReleaseError("registry differs from prepared plan")
    if plan.get("repository") != registry.repository:
        raise ReleaseError("plan repository differs from registry")
    families = plan.get("families")
    if not isinstance(families, dict) or set(families) != set(registry.families):
        raise ReleaseError("plan family set differs from registry")
    release = plan.get("release")
    if not isinstance(release, dict) or not re.fullmatch(
        r"[0-9a-f]{40}", str(release.get("commit", ""))
    ):
        raise ReleaseError("plan release commit is invalid")
    sandbox_suffix = release.get("sandbox_suffix")
    if sandbox_suffix is not None and (
        not isinstance(sandbox_suffix, str) or not SAFE_SUFFIX.fullmatch(sandbox_suffix)
    ):
        raise ReleaseError("plan sandbox suffix is invalid")
    plan_root = path.parent.resolve()
    for name, prepared in families.items():
        if not isinstance(prepared, dict):
            raise ReleaseError(f"family {name}: invalid plan entry")
        for field in ("base_archive", "archive", "root"):
            if not isinstance(prepared.get(field), str):
                raise ReleaseError(
                    f"family {name}: plan field {field} must be a string"
                )
            _safe_plan_name(prepared[field])
        assets = prepared.get("assets")
        if not isinstance(assets, list):
            raise ReleaseError(f"family {name}: plan assets must be an array")
        asset_names = []
        for asset in assets:
            if not isinstance(asset, dict):
                raise ReleaseError(f"family {name}: invalid planned asset")
            asset_names.append(_safe_plan_name(asset.get("source", "")))
            relative = PurePosixPath(asset.get("prepared_path", ""))
            if (
                relative.is_absolute()
                or not relative.parts
                or any(part in ("", ".", "..") for part in relative.parts)
            ):
                raise ReleaseError(f"family {name}: unsafe prepared asset path")
            candidate = (plan_root / Path(*relative.parts)).resolve()
            if not candidate.is_relative_to(plan_root) or not candidate.is_file():
                raise ReleaseError(
                    f"family {name}: prepared asset is outside the plan directory"
                )
            if not SHA256.fullmatch(str(asset.get("sha256", ""))):
                raise ReleaseError(f"family {name}: invalid planned asset SHA-256")
        if len(asset_names) != len(set(asset_names)):
            raise ReleaseError(f"family {name}: duplicate planned asset filename")
        base = (plan_root / prepared["base_archive"]).resolve()
        if not base.is_relative_to(plan_root) or not base.is_file():
            raise ReleaseError(
                f"family {name}: prepared archive is outside the plan directory"
            )
        base_digest = prepared.get("base_sha256")
        if not isinstance(base_digest, str) or not SHA256.fullmatch(base_digest):
            raise ReleaseError(f"family {name}: invalid prepared archive SHA-256")
        if hashlib.sha256(base.read_bytes()).hexdigest() != base_digest:
            raise ReleaseError(f"family {name}: prepared archive changed")
        if prepared.get("template_commit") != release["commit"]:
            raise ReleaseError(
                f"family {name}: template commit differs from release commit"
            )
        family = registry.families[name]
        for field, filename in (
            ("metadata_sha256", family.metadata_file),
            ("citation_sha256", family.citation_file),
        ):
            digest = prepared.get(field)
            if not isinstance(digest, str) or not SHA256.fullmatch(digest):
                raise ReleaseError(f"family {name}: invalid {field}")
            if hashlib.sha256((root / filename).read_bytes()).hexdigest() != digest:
                raise ReleaseError(
                    f"family {name}: {filename} differs from prepared plan"
                )
    return plan, registry


def _seed_metadata(
    root: Path,
    registry: Registry,
    family: Family,
    prepared: Mapping[str, Any],
    production: bool,
) -> dict[str, Any]:
    metadata = _metadata_template(root, family, registry)
    metadata["version"] = prepared["metadata_version"]
    metadata["publication_date"] = prepared["publication_date"]
    metadata["prereserve_doi"] = True
    metadata["related_identifiers"] = merge_relations(
        metadata.get("related_identifiers", []), family.marker, family.name, {}, {}
    )
    return metadata_for_environment(metadata, production)


def _existing_draft(lineage: Iterable[Mapping[str, Any]]) -> dict[str, Any] | None:
    drafts = [dict(record) for record in lineage if not _published(record)]
    if len(drafts) > 1:
        raise CollisionError("Zenodo concept has multiple unpublished drafts")
    return drafts[0] if drafts else None


def _matching_published(
    lineage: Iterable[Mapping[str, Any]], target_version: str
) -> list[dict[str, Any]]:
    records = [dict(record) for record in lineage]
    target = normalize_version(target_version)
    matching = [
        record
        for record in records
        if _published(record) and _record_version(record) == target
    ]
    if len(matching) > 1:
        raise CollisionError(f"duplicate published Zenodo version {target}")
    if matching and any(not _published(record) for record in records):
        raise CollisionError(f"published Zenodo version {target} has a leftover draft")
    return matching


def _production_preflight(
    registry: Registry,
    lineages: Mapping[str, list[dict[str, Any]]],
) -> None:
    missing = []
    for name, family in registry.families.items():
        baseline = family.production_baseline_version
        if not baseline:
            raise ReleaseError(
                f"family {name}: production_baseline_version is required"
            )
        published = [record for record in lineages[name] if _published(record)]
        versions = {_record_version(record) for record in published}
        if normalize_version(baseline) not in versions:
            missing.append(f"{name} {baseline}")
    if missing:
        raise ReleaseError(
            "production migration baselines are missing; initialize them with the dedicated "
            f"backfill workflow first: {', '.join(missing)}"
        )


def _safe_plan_name(value: str) -> str:
    if value != Path(value).name or value in ("", ".", ".."):
        raise ReleaseError(f"unsafe payload filename in plan: {value!r}")
    return value


def _finalize_payload(
    plan_path: Path,
    prepared: Mapping[str, Any],
    family: Family,
    version_doi: str,
    concept_doi: str,
    citation_template: str,
) -> dict[str, Path]:
    base_name = _safe_plan_name(prepared["base_archive"])
    archive_name = _safe_plan_name(prepared["archive"])
    base = plan_path.parent / base_name
    digest = hashlib.sha256(base.read_bytes()).hexdigest()
    if digest != prepared["base_sha256"]:
        raise ReleaseError(f"prepared archive changed: {base.name}")
    entries = read_archive(base)
    provenance_entry = next(
        entry
        for entry in entries
        if entry.name == f"{prepared['root']}/PROVENANCE.json"
    )
    provenance = json.loads(provenance_entry.data)
    provenance["concept_doi"] = concept_doi
    provenance["version_doi"] = version_doi
    citation = render_citation(
        citation_template,
        prepared["version"],
        prepared["publication_date"],
        version_doi,
        concept_doi,
    )
    final_entries = _with_generated_entries(
        entries, prepared["root"], citation, provenance
    )
    payload_dir = plan_path.parent / "payload" / family.name
    payload_dir.mkdir(parents=True, exist_ok=True)
    files = {archive_name: payload_dir / archive_name}
    write_archive(files[archive_name], final_entries, prepared["source_date_epoch"])
    files["CITATION.cff"] = payload_dir / "CITATION.cff"
    files["CITATION.cff"].write_text(citation, encoding="utf-8")
    files["LICENSE.md"] = payload_dir / "LICENSE.md"
    license_entry = next(
        entry
        for entry in final_entries
        if entry.name == f"{prepared['root']}/LICENSE.md"
    )
    files["LICENSE.md"].write_bytes(license_entry.data)
    files["PROVENANCE.json"] = payload_dir / "PROVENANCE.json"
    files["PROVENANCE.json"].write_bytes(_json_bytes(provenance))
    for asset in prepared["assets"]:
        name = _safe_plan_name(asset["source"])
        source = plan_path.parent / asset["prepared_path"]
        if hashlib.sha256(source.read_bytes()).hexdigest() != asset["sha256"]:
            raise ReleaseError(f"prepared canonical asset changed: {name}")
        files[name] = payload_dir / name
        files[name].write_bytes(source.read_bytes())
    if len(files) != 4 + len(prepared["assets"]):
        raise ReleaseError(f"family {family.name}: duplicate final payload filename")
    checksum_text = "".join(
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {name}\n"
        for name, path in sorted(files.items())
    )
    files["SHA256SUMS"] = payload_dir / "SHA256SUMS"
    files["SHA256SUMS"].write_text(checksum_text, encoding="utf-8")
    return files


def _md5_files(files: Mapping[str, Path]) -> dict[str, str]:
    return {
        name: hashlib.md5(path.read_bytes()).hexdigest() for name, path in files.items()
    }


def _publish(args: argparse.Namespace) -> None:
    root = Path(args.repository_root).resolve()
    plan_path = Path(args.plan).resolve()
    registry_path = (root / args.registry).resolve()
    plan, registry = _load_plan(plan_path, root, registry_path)
    client = ZenodoClient(args.token, args.api_url, args.transport)
    client.validate_plan_environment(plan)
    all_records = client.list_depositions()
    lineages = {
        name: select_lineage(
            all_records,
            family.marker,
            family.concept_recid if client.production else None,
        )
        for name, family in registry.families.items()
    }
    if client.production:
        _production_preflight(registry, lineages)

    staged: dict[str, dict[str, Any]] = {}
    legacy_skips: set[str] = set()
    for name, family in registry.families.items():
        prepared = plan["families"][name]
        target = normalize_version(prepared["version"])
        lineage = lineages[name]
        matching = _matching_published(lineage, target)
        if matching:
            record = matching[0]
            if is_legacy_linnet_match(family, record, target, client.production):
                legacy_skips.add(name)
            staged[name] = record
            continue
        if client.production:
            latest_published = next(
                (record for record in reversed(lineage) if _published(record)), None
            )
            if latest_published is not None:
                latest_version = _record_version(latest_published)
                if latest_version is None:
                    raise CollisionError(
                        f"family {name}: latest record has no valid version"
                    )
                if _version_order(target) < _version_order(latest_version):
                    raise CollisionError(
                        f"family {name}: target {target} regresses behind {latest_version}"
                    )
        draft = _existing_draft(lineage)
        if draft is not None:
            draft_version = _record_version(draft)
            latest_published = next(
                (record for record in reversed(lineage) if _published(record)), None
            )
            previous = _record_version(latest_published) if latest_published else None
            if draft_version not in {target, previous}:
                raise CollisionError(
                    f"family {name}: unrelated draft version {draft_version} blocks {target}"
                )
        elif lineage:
            latest = next(
                (record for record in reversed(lineage) if _published(record)), None
            )
            if latest is None:
                raise CollisionError(
                    f"family {name}: lineage has no usable published record"
                )
            draft = client.new_version(latest)
        else:
            if client.production:
                raise ReleaseError(
                    f"refusing to create a fresh production concept for {name}"
                )
            draft = client.create(
                _seed_metadata(root, registry, family, prepared, client.production)
            )
        staged[name] = client.update(
            draft["id"],
            _seed_metadata(root, registry, family, prepared, client.production),
        )

    concept_dois = {name: client.concept_doi(record) for name, record in staged.items()}
    version_dois = {name: client.version_doi(record) for name, record in staged.items()}
    finalized: dict[str, dict[str, Path]] = {}
    for name, family in registry.families.items():
        record = staged[name]
        prepared = plan["families"][name]
        citation_template = _git_text(
            root, prepared["template_commit"], family.citation_file
        )
        finalized[name] = _finalize_payload(
            plan_path,
            prepared,
            family,
            version_dois[name],
            concept_dois[name],
            citation_template,
        )

    release_data: dict[str, tuple[dict[str, Path], dict[str, str], dict[str, Any]]] = {}
    for name, family in registry.families.items():
        prepared = plan["families"][name]
        files = finalized[name]
        expected = _md5_files(files)
        expected_metadata = metadata_for_environment(
            build_metadata(
                _metadata_template(root, family, registry),
                family,
                prepared["metadata_version"],
                prepared["publication_date"],
                concept_dois,
                version_dois,
            ),
            client.production,
        )
        release_data[name] = (files, expected, expected_metadata)

    # Detect every known immutable collision before publishing the first draft.
    published_results: dict[str, dict[str, Any]] = {}
    for name in registry.families:
        record = staged[name]
        if not _published(record):
            continue
        _, expected, expected_metadata = release_data[name]
        if name not in legacy_skips:
            assert_same_version_files(
                record, plan["families"][name]["version"], expected
            )
            assert_stable_metadata(record, expected_metadata)
        published_results[name] = record

    for name in registry.families:
        record = staged[name]
        if _published(record):
            continue
        prepared = plan["families"][name]
        files, expected, expected_metadata = release_data[name]
        draft = client.update(record["id"], expected_metadata)
        uploaded = client.replace_files(draft, files)
        if _file_identity(uploaded) != expected:
            raise ReleaseError(f"family {name}: uploaded MD5 verification failed")
        result = client.publish(record["id"])
        if not _published(result) or _record_version(result) != normalize_version(
            prepared["version"]
        ):
            raise ReleaseError(f"family {name}: published record verification failed")
        if _file_identity(result) != expected:
            raise ReleaseError(f"family {name}: published file verification failed")
        assert_stable_metadata(result, expected_metadata)
        published_results[name] = result
    print(
        json.dumps(
            {
                name: {
                    "concept_doi": client.concept_doi(record),
                    "version_doi": client.version_doi(record),
                    "legacy_skip": name in legacy_skips,
                }
                for name, record in published_results.items()
            },
            indent=2,
            sort_keys=True,
        )
    )


def _audit(args: argparse.Namespace) -> None:
    root = Path(args.repository_root).resolve()
    plan_path = Path(args.plan).resolve()
    plan, registry = _load_plan(plan_path, root, (root / args.registry).resolve())
    client = ZenodoClient(args.token, args.api_url, args.transport)
    client.validate_plan_environment(plan)
    records = client.list_depositions()
    report: dict[str, Any] = {}
    failures = []
    discovered: dict[
        str, tuple[Family, list[dict[str, Any]], list[dict[str, Any]], bool]
    ] = {}
    for name, family in registry.families.items():
        lineage = select_lineage(
            records, family.marker, family.concept_recid if client.production else None
        )
        target = normalize_version(plan["families"][name]["version"])
        matches = [
            record
            for record in lineage
            if _published(record) and _record_version(record) == target
        ]
        legacy = bool(
            len(matches) == 1
            and is_legacy_linnet_match(family, matches[0], target, client.production)
        )
        discovered[name] = (family, lineage, matches, legacy)
    complete = all(len(matches) == 1 for _, _, matches, _ in discovered.values())
    concept_dois = (
        {
            name: client.concept_doi(matches[0])
            for name, (_, _, matches, _) in discovered.items()
        }
        if complete
        else {}
    )
    version_dois = (
        {
            name: client.version_doi(matches[0])
            for name, (_, _, matches, _) in discovered.items()
        }
        if complete
        else {}
    )
    for name, (family, lineage, matches, legacy) in discovered.items():
        status = "published" if len(matches) == 1 else "missing"
        if len(matches) > 1:
            status = "collision"
        if len(matches) != 1:
            failures.append(name)
        elif any(not _published(record) for record in lineage):
            status = "leftover-draft"
            failures.append(name)
        elif not legacy:
            payload_dir = plan_path.parent / "payload" / name
            paths = (
                {path.name: path for path in payload_dir.iterdir()}
                if payload_dir.is_dir()
                else {}
            )
            if not paths:
                status = "payload-missing"
                failures.append(name)
            else:
                if _file_identity(matches[0]) != _md5_files(paths):
                    status = "file-mismatch"
                    failures.append(name)
                elif complete:
                    expected_metadata = metadata_for_environment(
                        build_metadata(
                            _metadata_template(root, family, registry),
                            family,
                            plan["families"][name]["metadata_version"],
                            plan["families"][name]["publication_date"],
                            concept_dois,
                            version_dois,
                        ),
                        client.production,
                    )
                    try:
                        assert_stable_metadata(matches[0], expected_metadata)
                    except CollisionError:
                        status = "metadata-mismatch"
                        failures.append(name)
        report[name] = {
            "status": status,
            "legacy": legacy,
            "concept_doi": client.concept_doi(matches[0])
            if len(matches) == 1
            else None,
            "version_doi": client.version_doi(matches[0])
            if len(matches) == 1
            else None,
            "drafts": [
                record.get("id") for record in lineage if not _published(record)
            ],
        }
    print(json.dumps(report, indent=2, sort_keys=True))
    if failures:
        raise ReleaseError(f"Zenodo audit failed for: {', '.join(failures)}")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--repository-root", default=".")
    common.add_argument("--registry", default=".zenodo/records.json")
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare", parents=[common])
    prepare.add_argument("--release-version", required=True)
    prepare.add_argument("--release-tag", required=True)
    prepare.add_argument("--release-commit", required=True)
    prepare.add_argument("--assets-dir", required=True)
    prepare.add_argument("--output-dir", required=True)
    prepare.add_argument("--plan", required=True)
    prepare.add_argument("--version-suffix")
    prepare.set_defaults(handler=_prepare)
    for command, handler in (("publish", _publish), ("audit", _audit)):
        operation = subparsers.add_parser(command, parents=[common])
        operation.add_argument("--plan", required=True)
        operation.add_argument(
            "--api-url",
            default=os.environ.get("ZENODO_API_URL", PRODUCTION_API),
            choices=sorted(ALLOWED_APIS),
        )
        operation.add_argument("--token", default=os.environ.get("ZENODO_TOKEN"))
        operation.set_defaults(handler=handler, transport=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    try:
        arguments = _parser().parse_args(argv)
        arguments.handler(arguments)
    except (ReleaseError, OSError, ValueError) as error:
        print(f"zenodo: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
