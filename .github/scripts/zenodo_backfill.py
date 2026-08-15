"""Prepare, publish, and audit the one-time historical Zenodo migration.

Preparation is deliberately token-free.  Publishing is deliberately limited to
Zenodo production and requires an exact typed confirmation.  Historical inputs
are accepted only when they are pinned by ``.zenodo/backfill.json``.
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime as dt
import email.parser
import gzip
import hashlib
import io
import json
import os
import re
import sys
import tarfile
import tomllib
import urllib.parse
from collections.abc import Iterable, Mapping
from itertools import pairwise
from pathlib import Path, PurePosixPath
from typing import Any

import zenodo

CONFIRMATION = "PUBLISH HISTORICAL ZENODO BACKFILL"
EXPECTED_FAMILY_ORDER = (
    "gammaloop-monorepo",
    "gammaloop",
    "vakint",
    "spenso",
    "idenso",
    "linnet",
)
EXPECTED_VERSIONS = {
    "gammaloop-monorepo": ("0.2.0", "0.2.1", "0.3.0", "0.3.3", "0.3.4"),
    "gammaloop": (
        "0.0.1",
        "0.1.0",
        "0.1.1",
        "0.2.0",
        "0.2.1",
        "0.3.0",
        "0.3.1",
        "0.3.2",
        "0.3.3",
        "0.3.4",
    ),
    "vakint": ("0.1.0", "0.1.2"),
    "spenso": ("0.6.0",),
    "idenso": ("0.3.0",),
    "linnet": ("0.17.0",),
}
EXPECTED_COMPONENT_KEYS = {
    "gammaloop-monorepo": {
        version: ("gammalooprs", "gammaloop-python")
        for version in EXPECTED_VERSIONS["gammaloop-monorepo"]
    },
    "gammaloop": {
        "0.0.1": ("gammaloop-python",),
        "0.1.0": ("gammaloop-python",),
        "0.1.1": ("gammaloop-python",),
        "0.2.0": ("gammalooprs", "gammaloop-python"),
        "0.2.1": ("gammalooprs", "gammaloop-python"),
        "0.3.0": ("gammalooprs", "gammaloop-python"),
        "0.3.1": ("gammaloop-python",),
        "0.3.2": ("gammaloop-python",),
        "0.3.3": ("gammalooprs", "gammaloop-python"),
        "0.3.4": ("gammalooprs", "gammaloop-python"),
    },
    "vakint": {"0.1.0": ("vakint",), "0.1.2": ("vakint",)},
    "spenso": {"0.6.0": ("spenso", "spenso-macros", "spenso-hep-lib", "spynso3")},
    "idenso": {"0.3.0": ("idenso",)},
    "linnet": {"0.17.0": ("linnet", "clinnet", "linnest")},
}
EXPECTED_ARTIFACT_KINDS = {
    "gammaloop-monorepo": {
        version: ("git-archive",) for version in EXPECTED_VERSIONS["gammaloop-monorepo"]
    },
    "gammaloop": {
        "0.0.1": ("pypi-sdist",),
        "0.1.0": ("pypi-sdist",),
        "0.1.1": ("pypi-sdist",),
        "0.2.0": ("git-archive", "pypi-sdist"),
        "0.2.1": ("git-archive",),
        "0.3.0": ("git-archive", "pypi-sdist"),
        "0.3.1": ("pypi-sdist",),
        "0.3.2": ("pypi-sdist",),
        "0.3.3": ("git-archive", "pypi-sdist"),
        "0.3.4": ("git-archive",),
    },
    "vakint": {"0.1.0": ("crates-io",), "0.1.2": ("crates-io",)},
    "spenso": {"0.6.0": ("crates-io",) * 4},
    "idenso": {"0.3.0": ("crates-io",)},
    "linnet": {"0.17.0": ()},
}
EXPECTED_HOSTS = {
    "crates_io": ("crates.io", "static.crates.io"),
    "pypi": ("files.pythonhosted.org",),
    "zenodo": ("zenodo.org",),
    "github_api": ("api.github.com",),
}
HEX40 = re.compile(r"[0-9a-f]{40}")
HEX32 = re.compile(r"[0-9a-f]{32}")


@dataclasses.dataclass(frozen=True)
class HistoricalVersion:
    family: str
    version: str
    publication_date: str
    component_versions: dict[str, str]
    artifacts: tuple[dict[str, Any], ...]
    adoption: dict[str, Any] | None


@dataclasses.dataclass(frozen=True)
class HistoricalFamily:
    name: str
    registry_family: str
    lineage: str
    versions: tuple[HistoricalVersion, ...]
    production_only: bool


@dataclasses.dataclass(frozen=True)
class BackfillManifest:
    path: Path
    sha256: str
    repository: str
    hosts: dict[str, tuple[str, ...]]
    families: tuple[HistoricalFamily, ...]

    @staticmethod
    def _object(value: Any, label: str) -> dict[str, Any]:
        if not isinstance(value, dict):
            raise zenodo.ReleaseError(f"{label} must be an object")
        return value

    @staticmethod
    def _keys(
        value: Mapping[str, Any], required: set[str], optional: set[str], label: str
    ) -> None:
        keys = set(value)
        missing = required - keys
        unknown = keys - required - optional
        if missing or unknown:
            raise zenodo.ReleaseError(
                f"{label} schema differs: missing={sorted(missing)}, unknown={sorted(unknown)}"
            )

    @staticmethod
    def _text(value: Any, label: str) -> str:
        if not isinstance(value, str) or not value:
            raise zenodo.ReleaseError(f"{label} must be a non-empty string")
        return value

    @staticmethod
    def _date(value: Any, label: str) -> str:
        text = BackfillManifest._text(value, label)
        try:
            parsed = dt.date.fromisoformat(text)
        except ValueError as error:
            raise zenodo.ReleaseError(f"{label} must be an ISO date") from error
        if parsed.isoformat() != text:
            raise zenodo.ReleaseError(f"{label} must be a canonical ISO date")
        return text

    @staticmethod
    def _safe_filename(value: Any, label: str) -> str:
        name = BackfillManifest._text(value, label)
        if name != Path(name).name or name in (".", ".."):
            raise zenodo.ReleaseError(f"{label} is not a safe filename")
        if name.lower().endswith(".whl"):
            raise zenodo.ReleaseError("wheel files are forbidden in Zenodo backfill")
        return name

    @classmethod
    def load(cls, path: Path, registry: zenodo.Registry) -> BackfillManifest:
        try:
            data = path.read_bytes()
            document = json.loads(data)
        except (OSError, json.JSONDecodeError) as error:
            raise zenodo.ReleaseError(
                f"cannot read backfill manifest {path}: {error}"
            ) from error
        document = cls._object(document, "backfill manifest")
        cls._keys(
            document,
            {
                "schema",
                "repository",
                "redirect_host_allowlists",
                "relation_policy",
                "families",
            },
            set(),
            "backfill manifest",
        )
        if document["schema"] != 1:
            raise zenodo.ReleaseError(".zenodo/backfill.json must have schema 1")
        repository = cls._text(document["repository"], "backfill repository")
        if repository != registry.repository:
            raise zenodo.ReleaseError("backfill and registry repositories differ")
        hosts_raw = cls._object(
            document["redirect_host_allowlists"], "redirect allowlists"
        )
        if set(hosts_raw) != set(EXPECTED_HOSTS):
            raise zenodo.ReleaseError(
                "backfill redirect allowlist categories differ from policy"
            )
        hosts: dict[str, tuple[str, ...]] = {}
        for name, expected in EXPECTED_HOSTS.items():
            values = hosts_raw[name]
            if not isinstance(values, list) or tuple(values) != expected:
                raise zenodo.ReleaseError(
                    f"backfill {name} redirect allowlist differs from policy"
                )
            hosts[name] = tuple(values)
        policy = cls._object(document["relation_policy"], "relation policy")
        cls._keys(
            policy,
            {
                "umbrella_family",
                "historical_component_relation",
                "historical_component_target",
                "reason",
            },
            set(),
            "relation policy",
        )
        expected_policy = {
            "umbrella_family": "gammaloop-monorepo",
            "historical_component_relation": "isPartOf",
            "historical_component_target": "umbrella-concept",
        }
        if any(policy.get(key) != value for key, value in expected_policy.items()):
            raise zenodo.ReleaseError(
                "historical relation policy differs from the approved policy"
            )
        cls._text(policy.get("reason"), "relation policy reason")

        families_raw = cls._object(document["families"], "backfill families")
        if tuple(families_raw) != EXPECTED_FAMILY_ORDER:
            raise zenodo.ReleaseError(
                "backfill family ordering differs from the approved migration"
            )
        if set(families_raw) != set(registry.families):
            raise zenodo.ReleaseError("backfill and registry family sets differ")
        families = tuple(
            cls._family(name, cls._object(raw, f"family {name}"), registry)
            for name, raw in families_raw.items()
        )
        if sum(len(family.versions) for family in families) != 20:
            raise zenodo.ReleaseError(
                "historical migration must contain exactly 20 events"
            )
        return cls(
            path.resolve(),
            hashlib.sha256(data).hexdigest(),
            repository,
            hosts,
            families,
        )

    @classmethod
    def _family(
        cls, name: str, raw: dict[str, Any], registry: zenodo.Registry
    ) -> HistoricalFamily:
        cls._keys(
            raw,
            {"registry_family", "lineage", "versions"},
            {
                "historical_layout_reason",
                "legacy_concept_doi",
                "legacy_relation",
                "production_only",
            },
            f"family {name}",
        )
        if raw["registry_family"] != name:
            raise zenodo.ReleaseError(
                f"family {name}: registry_family must equal its key"
            )
        lineage = raw["lineage"]
        expected_lineage = "continue-production" if name == "linnet" else "new"
        if lineage != expected_lineage:
            raise zenodo.ReleaseError(
                f"family {name}: lineage must be {expected_lineage}"
            )
        production_only = raw.get("production_only", False)
        if not isinstance(production_only, bool) or production_only != (
            name == "linnet"
        ):
            raise zenodo.ReleaseError(f"family {name}: invalid production_only policy")
        legacy_doi = raw.get("legacy_concept_doi")
        legacy_relation = raw.get("legacy_relation")
        if name in {"spenso", "idenso"}:
            if (
                legacy_doi != "10.5281/zenodo.15913113"
                or legacy_relation != "isDerivedFrom"
            ):
                raise zenodo.ReleaseError(
                    f"family {name}: invalid legacy lineage relation"
                )
        elif legacy_doi is not None or legacy_relation is not None:
            raise zenodo.ReleaseError(
                f"family {name}: unexpected legacy lineage relation"
            )
        if registry.families[name].legacy_concept_doi != legacy_doi:
            raise zenodo.ReleaseError(
                f"family {name}: legacy concept differs from registry"
            )
        if "historical_layout_reason" in raw:
            cls._text(raw["historical_layout_reason"], f"family {name} layout reason")
        versions_raw = raw["versions"]
        if not isinstance(versions_raw, list) or not versions_raw:
            raise zenodo.ReleaseError(
                f"family {name}: versions must be a non-empty array"
            )
        versions = tuple(
            cls._version(name, index, item) for index, item in enumerate(versions_raw)
        )
        if tuple(item.version for item in versions) != EXPECTED_VERSIONS[name]:
            raise zenodo.ReleaseError(
                f"family {name}: historical version list differs from policy"
            )
        ordered = [zenodo._version_order(item.version) for item in versions]
        if any(left >= right for left, right in pairwise(ordered)):
            raise zenodo.ReleaseError(
                f"family {name}: versions are not strictly increasing"
            )
        dates = [item.publication_date for item in versions]
        if dates != sorted(dates):
            raise zenodo.ReleaseError(
                f"family {name}: publication dates are not ordered"
            )
        baseline = registry.families[name].production_baseline_version
        if (
            baseline is None
            or zenodo.normalize_version(baseline) != versions[-1].version
        ):
            raise zenodo.ReleaseError(
                f"family {name}: final backfill version is not its baseline"
            )
        return HistoricalFamily(
            name, raw["registry_family"], lineage, versions, production_only
        )

    @classmethod
    def _version(cls, family: str, index: int, value: Any) -> HistoricalVersion:
        raw = cls._object(value, f"family {family} version {index}")
        required = {"version", "publication_date", "component_versions"}
        optional = {
            "artifacts",
            "mode",
            "adoption",
            "absent_components",
            "historical_tag_policy",
            "source_version_mismatch",
            "no_repository_tag_reason",
        }
        cls._keys(raw, required, optional, f"family {family} version {index}")
        version = zenodo._stable_version(
            cls._text(raw["version"], "historical version")
        )
        publication_date = cls._date(
            raw["publication_date"], "historical publication_date"
        )
        components = cls._object(
            raw["component_versions"], "historical component_versions"
        )
        if not components or not all(
            isinstance(key, str) and key for key in components
        ):
            raise zenodo.ReleaseError("historical component_versions must not be empty")
        normalized_components = {
            key: zenodo._stable_version(cls._text(item, f"component {key} version"))
            for key, item in components.items()
        }
        if tuple(normalized_components) != EXPECTED_COMPONENT_KEYS[family][version]:
            raise zenodo.ReleaseError(
                f"family {family} {version}: component set or ordering differs from policy"
            )
        adoption = raw.get("adoption")
        if family == "linnet":
            if raw.get("mode") != "adopt" or "artifacts" in raw:
                raise zenodo.ReleaseError(
                    "Linnet baseline must be an adoption without artifacts"
                )
            adoption = cls._adoption(cls._object(adoption, "Linnet adoption"), version)
            absent = cls._object(
                raw.get("absent_components"), "Linnet absent_components"
            )
            if set(absent) != {"linnet-py"} or not all(
                isinstance(reason, str) and reason for reason in absent.values()
            ):
                raise zenodo.ReleaseError(
                    "Linnet adoption must explain the absent linnet-py package"
                )
            artifacts: tuple[dict[str, Any], ...] = ()
        else:
            if "mode" in raw or adoption is not None or "absent_components" in raw:
                raise zenodo.ReleaseError(
                    f"family {family}: unexpected adoption fields"
                )
            artifacts_raw = raw.get("artifacts")
            if not isinstance(artifacts_raw, list) or not artifacts_raw:
                raise zenodo.ReleaseError(
                    f"family {family}: artifacts must be non-empty"
                )
            artifacts = tuple(
                cls._artifact(
                    family,
                    version,
                    normalized_components,
                    cls._object(item, "historical artifact"),
                )
                for item in artifacts_raw
            )
            filenames = [item["filename"] for item in artifacts]
            if len(filenames) != len(set(filenames)):
                raise zenodo.ReleaseError(
                    f"family {family} {version}: duplicate artifact filename"
                )
        if (
            tuple(item["kind"] for item in artifacts)
            != EXPECTED_ARTIFACT_KINDS[family][version]
        ):
            raise zenodo.ReleaseError(
                f"family {family} {version}: artifact set or ordering differs from policy"
            )
        crate_names = tuple(
            item["crate"] for item in artifacts if item["kind"] == "crates-io"
        )
        if crate_names and crate_names != tuple(normalized_components):
            raise zenodo.ReleaseError(
                f"family {family} {version}: crate archive order differs from components"
            )
        for artifact in artifacts:
            if artifact["kind"] != "crates-io":
                continue
            try:
                published = dt.datetime.fromisoformat(artifact["published_at"])
            except ValueError as error:
                raise zenodo.ReleaseError(
                    "crate published_at must be an ISO timestamp"
                ) from error
            if published.date() > dt.date.fromisoformat(publication_date):
                raise zenodo.ReleaseError(
                    "crate was published after the historical family record"
                )
        has_git = any(item["kind"] == "git-archive" for item in artifacts)
        if has_git:
            policy = cls._object(
                raw.get("historical_tag_policy"), "historical tag policy"
            )
            cls._keys(policy, {"lightweight", "reason"}, set(), "historical tag policy")
            if policy["lightweight"] is not True:
                raise zenodo.ReleaseError(
                    "historical Git tag policy must explicitly allow lightweight tags"
                )
            cls._text(policy["reason"], "historical tag reason")
        elif family != "linnet":
            cls._text(
                raw.get("no_repository_tag_reason"), "missing repository tag reason"
            )
        mismatch = raw.get("source_version_mismatch")
        if mismatch is not None:
            mismatch = cls._object(mismatch, "source_version_mismatch")
            cls._keys(
                mismatch,
                {"record_version", "embedded_versions", "reason"},
                set(),
                "source_version_mismatch",
            )
            if mismatch["record_version"] != version or not isinstance(
                mismatch["embedded_versions"], dict
            ):
                raise zenodo.ReleaseError(
                    "source_version_mismatch does not describe this version"
                )
            cls._text(mismatch["reason"], "source version mismatch reason")
        return HistoricalVersion(
            family,
            version,
            publication_date,
            normalized_components,
            artifacts,
            adoption,
        )

    @classmethod
    def _artifact(
        cls,
        family: str,
        release_version: str,
        components: Mapping[str, str],
        raw: dict[str, Any],
    ) -> dict[str, Any]:
        kind = raw.get("kind")
        shared = {"kind", "filename", "sha256"}
        if kind == "git-archive":
            if family not in {"gammaloop-monorepo", "gammaloop"}:
                raise zenodo.ReleaseError(
                    f"family {family}: Git archives are not approved"
                )
            cls._keys(
                raw,
                {
                    "kind",
                    "filename",
                    "sha256",
                    "repository",
                    "ref",
                    "tag_object_type",
                    "commit",
                    "commit_date",
                    "source_date_epoch",
                    "archive_paths",
                    "archive_prefix",
                },
                set(),
                "git archive",
            )
            if raw["repository"] != "https://github.com/alphal00p/gammaloop":
                raise zenodo.ReleaseError(
                    "git archive repository differs from the approved repository"
                )
            if raw["ref"] != f"refs/tags/v{release_version}":
                raise zenodo.ReleaseError(
                    "git archive ref differs from the historical record version"
                )
            if raw["tag_object_type"] != "commit" or not HEX40.fullmatch(
                str(raw["commit"])
            ):
                raise zenodo.ReleaseError(
                    "historical Git tags must be pinned lightweight tags"
                )
            if (
                not isinstance(raw["source_date_epoch"], int)
                or raw["source_date_epoch"] <= 0
            ):
                raise zenodo.ReleaseError(
                    "git archive source_date_epoch must be positive"
                )
            try:
                timestamp = int(
                    dt.datetime.fromisoformat(raw["commit_date"]).timestamp()
                )
            except (TypeError, ValueError) as error:
                raise zenodo.ReleaseError(
                    "git archive commit_date must be an ISO timestamp"
                ) from error
            if timestamp != raw["source_date_epoch"]:
                raise zenodo.ReleaseError(
                    "git archive commit_date and source_date_epoch differ"
                )
            paths = raw["archive_paths"]
            if paths != ["."]:
                raise zenodo.ReleaseError(
                    "historical Git archives must contain the pinned full tree"
                )
            for path in paths:
                if path != ".":
                    zenodo._safe_repo_path(cls._text(path, "git archive path"))
            prefix = cls._text(raw["archive_prefix"], "git archive prefix")
            if not prefix.endswith("/"):
                raise zenodo.ReleaseError("git archive prefix must end in /")
            zenodo._safe_repo_path(prefix.rstrip("/"))
            expected_stem = (
                f"gammaloop-monorepo-{release_version}"
                if family == "gammaloop-monorepo"
                else f"gammaloop-source-{release_version}"
            )
            expected_prefix = (
                f"gammaloop-monorepo-{release_version}/"
                if family == "gammaloop-monorepo"
                else f"gammaloop-{release_version}/"
            )
            if (
                raw["filename"] != f"{expected_stem}.tar.gz"
                or prefix != expected_prefix
            ):
                raise zenodo.ReleaseError(
                    "historical Git archive filename or prefix differs from policy"
                )
        elif kind == "crates-io":
            if family not in {"vakint", "spenso", "idenso"}:
                raise zenodo.ReleaseError(
                    f"family {family}: crates.io artifacts are not approved"
                )
            cls._keys(
                raw,
                shared | {"crate", "version", "url", "published_at", "cargo_vcs_info"},
                set(),
                "crates.io artifact",
            )
            crate = cls._text(raw["crate"], "crate name")
            version = zenodo._stable_version(cls._text(raw["version"], "crate version"))
            if raw["filename"] != f"{crate}-{version}.crate":
                raise zenodo.ReleaseError(
                    "crates.io artifact filename differs from crate/version"
                )
            if components.get(crate) != version:
                raise zenodo.ReleaseError(
                    "crates.io artifact differs from event component versions"
                )
            if (
                raw["url"]
                != f"https://crates.io/api/v1/crates/{crate}/{version}/download"
            ):
                raise zenodo.ReleaseError("crates.io artifact URL is not canonical")
            cls._text(raw["published_at"], "crate published_at")
            vcs = cls._object(raw["cargo_vcs_info"], "cargo_vcs_info")
            cls._keys(vcs, {"sha1", "path_in_vcs"}, set(), "cargo_vcs_info")
            if not HEX40.fullmatch(str(vcs["sha1"])) or not isinstance(
                vcs["path_in_vcs"], str
            ):
                raise zenodo.ReleaseError("invalid cargo_vcs_info pin")
            if vcs["path_in_vcs"]:
                zenodo._safe_repo_path(vcs["path_in_vcs"])
        elif kind == "pypi-sdist":
            if family != "gammaloop":
                raise zenodo.ReleaseError(
                    f"family {family}: PyPI artifacts are not approved"
                )
            cls._keys(
                raw,
                shared | {"project", "url"},
                set(),
                "PyPI source distribution",
            )
            project = cls._text(raw["project"], "PyPI project")
            if (
                project != family
                or components.get("gammaloop-python") != release_version
                or raw["filename"] != f"{project}-{release_version}.tar.gz"
            ):
                raise zenodo.ReleaseError(
                    "PyPI source distribution identity differs from release"
                )
            parsed = urllib.parse.urlsplit(str(raw["url"]))
            if parsed.scheme != "https" or parsed.hostname != "files.pythonhosted.org":
                raise zenodo.ReleaseError(
                    "PyPI source distribution URL is not canonical"
                )
        else:
            raise zenodo.ReleaseError(f"unsupported historical artifact kind: {kind!r}")
        cls._safe_filename(raw.get("filename"), "historical artifact filename")
        digest = raw.get("sha256")
        if not zenodo.SHA256.fullmatch(str(digest)):
            raise zenodo.ReleaseError(
                "historical artifact SHA-256 is missing or invalid"
            )
        if any(
            str(value).lower().endswith(".whl")
            for value in raw.values()
            if isinstance(value, str)
        ):
            raise zenodo.ReleaseError("wheel files are forbidden in Zenodo backfill")
        return json.loads(json.dumps(raw))

    @classmethod
    def _adoption(cls, raw: dict[str, Any], version: str) -> dict[str, Any]:
        cls._keys(
            raw,
            {
                "concept_recid",
                "concept_doi",
                "record_id",
                "record_doi",
                "published_at",
                "title",
                "metadata_version",
                "normalized_version",
                "license",
                "creators",
                "related_tag_url",
                "tag_ref",
                "tag_object",
                "tag_object_type",
                "commit",
                "files",
                "reason",
            },
            set(),
            "Linnet adoption",
        )
        exact = {
            "concept_recid": "15494393",
            "concept_doi": "10.5281/zenodo.15494393",
            "record_id": "18429583",
            "record_doi": "10.5281/zenodo.18429583",
            "normalized_version": version,
            "tag_object_type": "tag",
            "related_tag_url": "https://github.com/alphal00p/linnet/tree/linnet-v0.17.0",
            "tag_ref": "refs/tags/linnet-v0.17.0",
        }
        if any(str(raw.get(key)) != expected for key, expected in exact.items()):
            raise zenodo.ReleaseError(
                "Linnet adoption identity differs from the approved record"
            )
        if zenodo.normalize_version(str(raw["metadata_version"])) != version:
            raise zenodo.ReleaseError(
                "Linnet adopted metadata version differs from baseline"
            )
        for field in (
            "published_at",
            "title",
            "license",
            "related_tag_url",
            "tag_ref",
            "reason",
        ):
            cls._text(raw[field], f"Linnet adoption {field}")
        if not HEX40.fullmatch(str(raw["tag_object"])) or not HEX40.fullmatch(
            str(raw["commit"])
        ):
            raise zenodo.ReleaseError("Linnet adoption Git objects are not pinned")
        if not isinstance(raw["creators"], list) or not all(
            isinstance(name, str) and name for name in raw["creators"]
        ):
            raise zenodo.ReleaseError("Linnet adoption creators are invalid")
        files = raw["files"]
        if not isinstance(files, list) or len(files) != 1:
            raise zenodo.ReleaseError("Linnet adoption must pin exactly one file")
        item = cls._object(files[0], "Linnet adoption file")
        cls._keys(
            item, {"filename", "url", "size", "md5", "sha256"}, set(), "adoption file"
        )
        if not isinstance(item["size"], int) or item["size"] <= 0:
            raise zenodo.ReleaseError("Linnet adoption file size is invalid")
        if not HEX32.fullmatch(str(item["md5"])) or not zenodo.SHA256.fullmatch(
            str(item["sha256"])
        ):
            raise zenodo.ReleaseError("Linnet adoption file checksums are invalid")
        parsed = urllib.parse.urlsplit(str(item["url"]))
        if parsed.scheme != "https" or parsed.hostname != "zenodo.org":
            raise zenodo.ReleaseError(
                "Linnet adoption file URL is not on Zenodo production"
            )
        cls._text(item["filename"], "Linnet adoption filename")
        return json.loads(json.dumps(raw))

    def events(self) -> Iterable[HistoricalVersion]:
        for family in self.families:
            yield from family.versions

    def family(self, name: str) -> HistoricalFamily:
        return next(family for family in self.families if family.name == name)


class Downloader:
    """Fetch immutable public inputs while enforcing every redirect origin."""

    def __init__(self, transport: zenodo.Transport | None = None):
        self.transport = transport or zenodo.UrllibTransport()

    @staticmethod
    def _trusted(url: str, hosts: tuple[str, ...]) -> str:
        parsed = urllib.parse.urlsplit(url)
        if (
            parsed.scheme != "https"
            or parsed.hostname not in hosts
            or parsed.username is not None
            or parsed.password is not None
            or parsed.port is not None
            or parsed.fragment
        ):
            raise zenodo.ReleaseError(f"untrusted historical artifact URL: {url}")
        return url

    def fetch(self, url: str, hosts: tuple[str, ...], accept: str = "*/*") -> bytes:
        current = self._trusted(url, hosts)
        headers = {
            "Accept": accept,
            "User-Agent": "gammaloop-zenodo-backfill/1",
        }
        for _ in range(6):
            response = self.transport.request("GET", current, headers, None)
            if response.status == 200:
                return response.body
            if response.status not in {301, 302, 303, 307, 308}:
                detail = response.body.decode(errors="replace")[:500]
                raise zenodo.ReleaseError(
                    f"historical download {current} returned {response.status}: {detail}"
                )
            location = response.headers.get("Location") or response.headers.get(
                "location"
            )
            if not location:
                raise zenodo.ReleaseError(
                    "historical download redirect omitted Location"
                )
            current = self._trusted(urllib.parse.urljoin(current, location), hosts)
        raise zenodo.ReleaseError("historical download exceeded five redirects")


class BackfillPreparer:
    def __init__(
        self,
        root: Path,
        manifest: BackfillManifest,
        registry: zenodo.Registry,
        downloader: Downloader | None = None,
    ):
        self.root = root.resolve()
        self.manifest = manifest
        self.registry = registry
        self.downloader = downloader or Downloader()

    def prepare(self, output: Path, plan_path: Path) -> dict[str, Any]:
        output = output.resolve()
        plan_path = plan_path.resolve()
        if plan_path.parent != output:
            raise zenodo.ReleaseError("--plan must be directly inside --output-dir")
        output.mkdir(parents=True, exist_ok=True)
        template_hashes = self._template_hashes()
        events = []
        generated_hashes: dict[str, str] = {}
        for version in self.manifest.events():
            event_dir = output / "canonical" / version.family / version.version
            event_dir.mkdir(parents=True, exist_ok=True)
            assets = []
            if version.adoption is not None:
                assets.append(self._prepare_adoption(version, event_dir, output))
            else:
                for artifact in version.artifacts:
                    prepared = self._prepare_artifact(
                        version, artifact, event_dir, output
                    )
                    assets.append(prepared)
                    if artifact["kind"] == "git-archive":
                        key = (
                            f"{version.family}/{version.version}/{artifact['filename']}"
                        )
                        generated_hashes[key] = prepared["sha256"]
            events.append(
                {
                    "family": version.family,
                    "version": version.version,
                    "metadata_version": self.registry.families[
                        version.family
                    ].tag.format(version=version.version),
                    "publication_date": version.publication_date,
                    "component_versions": version.component_versions,
                    "mode": "adopt" if version.adoption is not None else "publish",
                    "assets": assets,
                }
            )
        report = output / "git-archive-sha256.json"
        report.write_bytes(zenodo._json_bytes(generated_hashes))
        plan = {
            "schema": 1,
            "kind": "gammaloop-zenodo-historical-backfill",
            "repository": self.registry.repository,
            "registry_sha256": hashlib.sha256(
                (self.root / ".zenodo/records.json").read_bytes()
            ).hexdigest(),
            "backfill_manifest_sha256": self.manifest.sha256,
            "template_sha256": template_hashes,
            "events": events,
            "git_archive_hashes": generated_hashes,
        }
        plan_path.parent.mkdir(parents=True, exist_ok=True)
        plan_path.write_bytes(zenodo._json_bytes(plan))
        print(
            json.dumps(
                {"plan": str(plan_path), "git_archive_hashes": generated_hashes},
                indent=2,
            )
        )
        return plan

    def _template_hashes(self) -> dict[str, str]:
        paths = {"LICENSE.md"}
        for family in self.registry.families.values():
            paths.add(family.metadata_file)
            paths.add(family.citation_file)
        result = {}
        for value in sorted(paths):
            path = self.root / value
            if not path.is_file():
                raise zenodo.ReleaseError(f"backfill template is missing: {value}")
            result[value] = hashlib.sha256(path.read_bytes()).hexdigest()
        return result

    def _prepare_artifact(
        self,
        version: HistoricalVersion,
        artifact: dict[str, Any],
        event_dir: Path,
        output: Path,
    ) -> dict[str, Any]:
        destination = event_dir / artifact["filename"]
        if artifact["kind"] == "git-archive":
            digest = self._git_archive(artifact, destination)
            if digest != artifact["sha256"]:
                raise zenodo.ReleaseError(
                    f"generated Git archive differs: {destination.name}"
                )
        else:
            category = "crates_io" if artifact["kind"] == "crates-io" else "pypi"
            data = self._download_pinned(
                artifact["url"],
                artifact["sha256"],
                self.manifest.hosts[category],
                destination,
            )
            digest = hashlib.sha256(data).hexdigest()
            if artifact["kind"] == "crates-io":
                self._verify_crate(data, artifact)
            else:
                self._verify_sdist(data, artifact)
        return {
            "kind": artifact["kind"],
            "filename": artifact["filename"],
            "prepared_path": destination.relative_to(output).as_posix(),
            "sha256": digest,
        }

    def _download_pinned(
        self, url: str, expected: str, hosts: tuple[str, ...], destination: Path
    ) -> bytes:
        if (
            destination.is_file()
            and hashlib.sha256(destination.read_bytes()).hexdigest() == expected
        ):
            return destination.read_bytes()
        data = self.downloader.fetch(url, hosts)
        if hashlib.sha256(data).hexdigest() != expected:
            raise zenodo.ReleaseError(
                f"SHA-256 mismatch for historical download {destination.name}"
            )
        destination.write_bytes(data)
        return data

    def _git_archive(self, artifact: Mapping[str, Any], destination: Path) -> str:
        ref = artifact["ref"]
        if (
            zenodo._git(self.root, "cat-file", "-t", ref).decode().strip()
            != artifact["tag_object_type"]
        ):
            raise zenodo.ReleaseError(f"historical ref changed object type: {ref}")
        commit = zenodo._commit(self.root, ref)
        if commit != artifact["commit"]:
            raise zenodo.ReleaseError(f"historical ref changed commit: {ref}")
        timestamp = int(
            zenodo._git(self.root, "show", "-s", "--format=%ct", commit).decode()
        )
        commit_date = (
            zenodo._git(self.root, "show", "-s", "--format=%cI", commit)
            .decode()
            .strip()
        )
        if (
            timestamp != artifact["source_date_epoch"]
            or commit_date != artifact["commit_date"]
        ):
            raise zenodo.ReleaseError(f"historical commit metadata changed: {ref}")
        prefix = artifact["archive_prefix"]
        raw = zenodo._git(
            self.root,
            "archive",
            "--format=tar",
            f"--prefix={prefix}",
            commit,
            *artifact["archive_paths"],
        )
        self._verify_git_tar(raw, prefix)
        output = io.BytesIO()
        with gzip.GzipFile(
            fileobj=output, mode="wb", filename="", mtime=0, compresslevel=9
        ) as zipped:
            zipped.write(raw)
        destination.write_bytes(output.getvalue())
        return hashlib.sha256(output.getvalue()).hexdigest()

    @staticmethod
    def _verify_git_tar(data: bytes, prefix: str) -> None:
        prefix_path = PurePosixPath(prefix.rstrip("/"))
        try:
            with tarfile.open(fileobj=io.BytesIO(data), mode="r:") as archive:
                members = archive.getmembers()
                if not members:
                    raise zenodo.ReleaseError("historical Git archive is empty")
                for member in members:
                    path = PurePosixPath(member.name.rstrip("/"))
                    zenodo._safe_repo_path(str(path))
                    if path.parts[: len(prefix_path.parts)] != prefix_path.parts:
                        raise zenodo.ReleaseError(
                            f"historical Git archive path escapes prefix: {member.name}"
                        )
                    if member.isfile() or member.isdir():
                        continue
                    if not member.issym():
                        raise zenodo.ReleaseError(
                            f"unsupported historical Git entry: {member.name}"
                        )
                    target = PurePosixPath(member.linkname)
                    if target.is_absolute():
                        raise zenodo.ReleaseError(
                            f"absolute historical Git symlink: {member.name}"
                        )
                    stack = list(path.parent.parts)
                    for part in target.parts:
                        if part in ("", "."):
                            continue
                        if part == "..":
                            if len(stack) <= len(prefix_path.parts):
                                raise zenodo.ReleaseError(
                                    f"historical Git symlink escapes archive: {member.name}"
                                )
                            stack.pop()
                        else:
                            stack.append(part)
        except tarfile.TarError as error:
            raise zenodo.ReleaseError(
                f"invalid historical Git archive: {error}"
            ) from error

    @staticmethod
    def _safe_tar(data: bytes, label: str) -> list[tuple[tarfile.TarInfo, bytes]]:
        result = []
        try:
            with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
                for member in archive:
                    path = PurePosixPath(member.name)
                    if path.is_absolute() or any(
                        part in ("", ".", "..") for part in path.parts
                    ):
                        raise zenodo.ReleaseError(
                            f"unsafe path in {label}: {member.name}"
                        )
                    if member.isdir():
                        continue
                    if not member.isfile():
                        raise zenodo.ReleaseError(
                            f"unsupported entry in {label}: {member.name}"
                        )
                    stream = archive.extractfile(member)
                    if stream is None:
                        raise zenodo.ReleaseError(
                            f"cannot read {member.name} from {label}"
                        )
                    result.append((member, stream.read()))
        except tarfile.TarError as error:
            raise zenodo.ReleaseError(
                f"invalid tar archive {label}: {error}"
            ) from error
        return result

    @classmethod
    def _verify_crate(cls, data: bytes, artifact: Mapping[str, Any]) -> None:
        entries = cls._safe_tar(data, artifact["filename"])
        root = f"{artifact['crate']}-{artifact['version']}"
        expected = f"{root}/.cargo_vcs_info.json"
        matches = [body for member, body in entries if member.name == expected]
        if len(matches) != 1:
            raise zenodo.ReleaseError(
                f"{artifact['filename']} has no unique cargo_vcs_info"
            )
        try:
            vcs = json.loads(matches[0])
        except json.JSONDecodeError as error:
            raise zenodo.ReleaseError(
                f"{artifact['filename']} has invalid cargo_vcs_info"
            ) from error
        expected_vcs = artifact["cargo_vcs_info"]
        actual = {
            "sha1": vcs.get("git", {}).get("sha1"),
            "path_in_vcs": vcs.get("path_in_vcs", ""),
        }
        if actual != expected_vcs:
            raise zenodo.ReleaseError(
                f"{artifact['filename']} cargo_vcs_info differs from manifest"
            )
        manifests = [
            body for member, body in entries if member.name == f"{root}/Cargo.toml"
        ]
        if len(manifests) != 1:
            raise zenodo.ReleaseError(
                f"{artifact['filename']} has no unique Cargo.toml"
            )
        try:
            package = tomllib.loads(manifests[0].decode())["package"]
        except (KeyError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
            raise zenodo.ReleaseError(
                f"{artifact['filename']} has an invalid Cargo.toml"
            ) from error
        if (
            package.get("name") != artifact["crate"]
            or package.get("version") != artifact["version"]
        ):
            raise zenodo.ReleaseError(
                f"{artifact['filename']} Cargo.toml differs from manifest"
            )

    @classmethod
    def _verify_sdist(cls, data: bytes, artifact: Mapping[str, Any]) -> None:
        entries = cls._safe_tar(data, artifact["filename"])
        matches = [
            body for member, body in entries if member.name.endswith("/PKG-INFO")
        ]
        if len(matches) != 1:
            raise zenodo.ReleaseError(f"{artifact['filename']} has no unique PKG-INFO")
        metadata = email.parser.BytesParser().parsebytes(matches[0])
        normalized = re.sub(r"[-_.]+", "-", str(metadata.get("Name", "")).lower())
        project = re.sub(r"[-_.]+", "-", artifact["project"].lower())
        expected_version = (
            artifact["filename"]
            .removesuffix(".tar.gz")
            .removeprefix(artifact["project"] + "-")
        )
        if normalized != project or metadata.get("Version") != expected_version:
            raise zenodo.ReleaseError(
                f"{artifact['filename']} PKG-INFO differs from manifest"
            )

    def _prepare_adoption(
        self, version: HistoricalVersion, event_dir: Path, output: Path
    ) -> dict[str, Any]:
        assert version.adoption is not None
        adoption = version.adoption
        record_url = f"https://zenodo.org/api/records/{adoption['record_id']}"
        record_data = self.downloader.fetch(
            record_url, self.manifest.hosts["zenodo"], "application/json"
        )
        try:
            record = json.loads(record_data)
        except json.JSONDecodeError as error:
            raise zenodo.ReleaseError(
                "Linnet adoption record returned invalid JSON"
            ) from error
        self._verify_adoption_record(record, adoption)
        ref_url = (
            "https://api.github.com/repos/alphal00p/linnet/git/ref/tags/linnet-v0.17.0"
        )
        tag_url = (
            "https://api.github.com/repos/alphal00p/linnet/git/tags/"
            f"{adoption['tag_object']}"
        )
        try:
            ref_record = json.loads(
                self.downloader.fetch(
                    ref_url,
                    self.manifest.hosts["github_api"],
                    "application/vnd.github+json",
                )
            )
            tag_record = json.loads(
                self.downloader.fetch(
                    tag_url,
                    self.manifest.hosts["github_api"],
                    "application/vnd.github+json",
                )
            )
        except json.JSONDecodeError as error:
            raise zenodo.ReleaseError(
                "Linnet Git tag API returned invalid JSON"
            ) from error
        self._verify_adoption_git(ref_record, tag_record, adoption)
        item = adoption["files"][0]
        destination = event_dir / Path(item["filename"]).name
        data = self._download_pinned(
            item["url"], item["sha256"], self.manifest.hosts["zenodo"], destination
        )
        if len(data) != item["size"] or hashlib.md5(data).hexdigest() != item["md5"]:
            raise zenodo.ReleaseError(
                "Linnet adoption file size or MD5 differs from manifest"
            )
        return {
            "kind": "zenodo-adoption",
            "filename": item["filename"],
            "prepared_path": destination.relative_to(output).as_posix(),
            "sha256": item["sha256"],
            "md5": item["md5"],
            "size": item["size"],
        }

    @staticmethod
    def _verify_adoption_record(
        record: Mapping[str, Any], adoption: Mapping[str, Any]
    ) -> None:
        metadata = record.get("metadata", {})
        license_value = metadata.get("license", {})
        if isinstance(license_value, dict):
            license_value = license_value.get("id")
        actual = {
            "concept_recid": str(record.get("conceptrecid")),
            "concept_doi": record.get("conceptdoi"),
            "record_id": str(record.get("id")),
            "record_doi": record.get("doi"),
            "published_at": record.get("created"),
            "title": metadata.get("title"),
            "metadata_version": metadata.get("version"),
            "license": license_value,
            "creators": [item.get("name") for item in metadata.get("creators", [])],
        }
        expected = {key: adoption[key] for key in actual}
        if actual != expected:
            raise zenodo.ReleaseError(
                "public Linnet adoption metadata differs from manifest"
            )
        relations = {
            item.get("identifier") for item in metadata.get("related_identifiers", [])
        }
        if adoption["related_tag_url"] not in relations:
            raise zenodo.ReleaseError(
                "public Linnet adoption record lost its tag relation"
            )
        files = {
            item.get("key") or item.get("filename"): (
                int(item.get("size", -1)),
                str(item.get("checksum", "")).removeprefix("md5:"),
            )
            for item in record.get("files", [])
        }
        expected_files = {
            item["filename"]: (item["size"], item["md5"]) for item in adoption["files"]
        }
        if files != expected_files:
            raise zenodo.ReleaseError(
                "public Linnet adoption files differ from manifest"
            )

    @staticmethod
    def _verify_adoption_git(
        ref_record: Mapping[str, Any],
        tag_record: Mapping[str, Any],
        adoption: Mapping[str, Any],
    ) -> None:
        if (
            ref_record.get("ref") != adoption["tag_ref"]
            or ref_record.get("object", {}).get("type") != adoption["tag_object_type"]
            or ref_record.get("object", {}).get("sha") != adoption["tag_object"]
            or tag_record.get("sha") != adoption["tag_object"]
            or tag_record.get("object", {}).get("type") != "commit"
            or tag_record.get("object", {}).get("sha") != adoption["commit"]
        ):
            raise zenodo.ReleaseError("public Linnet Git tag differs from manifest")


class PreparedBackfill:
    def __init__(
        self,
        root: Path,
        plan_path: Path,
        manifest: BackfillManifest,
        registry: zenodo.Registry,
    ):
        self.root = root.resolve()
        self.plan_path = plan_path.resolve()
        self.manifest = manifest
        self.registry = registry
        try:
            self.plan = json.loads(self.plan_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise zenodo.ReleaseError(
                f"cannot read prepared backfill plan: {error}"
            ) from error
        self._validate()

    def _validate(self) -> None:
        expected_keys = {
            "schema",
            "kind",
            "repository",
            "registry_sha256",
            "backfill_manifest_sha256",
            "template_sha256",
            "events",
            "git_archive_hashes",
        }
        if not isinstance(self.plan, dict) or set(self.plan) != expected_keys:
            raise zenodo.ReleaseError("prepared backfill plan schema differs")
        if (
            self.plan["schema"] != 1
            or self.plan["kind"] != "gammaloop-zenodo-historical-backfill"
        ):
            raise zenodo.ReleaseError("unsupported prepared backfill plan")
        if self.plan["repository"] != self.registry.repository:
            raise zenodo.ReleaseError("prepared backfill repository differs")
        registry_path = self.root / ".zenodo/records.json"
        if (
            hashlib.sha256(registry_path.read_bytes()).hexdigest()
            != self.plan["registry_sha256"]
        ):
            raise zenodo.ReleaseError(
                "Zenodo registry changed after backfill preparation"
            )
        if self.manifest.sha256 != self.plan["backfill_manifest_sha256"]:
            raise zenodo.ReleaseError("historical manifest changed after preparation")
        expected_templates = {"LICENSE.md"}
        for family in self.registry.families.values():
            expected_templates.update((family.metadata_file, family.citation_file))
        template_hashes = self.plan["template_sha256"]
        if (
            not isinstance(template_hashes, dict)
            or set(template_hashes) != expected_templates
        ):
            raise zenodo.ReleaseError(
                "prepared template hash map differs from registry"
            )
        for path, digest in template_hashes.items():
            if not zenodo.SHA256.fullmatch(str(digest)):
                raise zenodo.ReleaseError(f"invalid prepared template digest: {path}")
            candidate = self._resolved(self.root, path, "template")
            if (
                not candidate.is_file()
                or hashlib.sha256(candidate.read_bytes()).hexdigest() != digest
            ):
                raise zenodo.ReleaseError(
                    f"backfill template changed after preparation: {path}"
                )
        expected_events = [
            (item.family, item.version) for item in self.manifest.events()
        ]
        events = self.plan.get("events")
        if (
            not isinstance(events, list)
            or not all(isinstance(item, dict) for item in events)
            or [(item.get("family"), item.get("version")) for item in events]
            != expected_events
        ):
            raise zenodo.ReleaseError(
                "prepared event ordering differs from historical manifest"
            )
        planned_git_hashes: dict[str, str] = {}
        for event, historical in zip(events, self.manifest.events(), strict=True):
            expected_event_keys = {
                "family",
                "version",
                "metadata_version",
                "publication_date",
                "component_versions",
                "mode",
                "assets",
            }
            if set(event) != expected_event_keys:
                raise zenodo.ReleaseError("prepared event schema differs")
            family = self.registry.families[historical.family]
            expected_metadata_version = family.tag.format(version=historical.version)
            if event["metadata_version"] != expected_metadata_version:
                raise zenodo.ReleaseError(
                    "prepared metadata version differs from registry"
                )
            if event["publication_date"] != historical.publication_date:
                raise zenodo.ReleaseError(
                    "prepared publication date differs from historical manifest"
                )
            if event["component_versions"] != historical.component_versions:
                raise zenodo.ReleaseError(
                    "prepared component versions differ from historical manifest"
                )
            expected_mode = "adopt" if historical.adoption is not None else "publish"
            if event["mode"] != expected_mode:
                raise zenodo.ReleaseError(
                    "prepared event mode differs from historical manifest"
                )
            assets = event["assets"]
            sources = (
                historical.adoption["files"]
                if historical.adoption is not None
                else historical.artifacts
            )
            if not isinstance(assets, list) or len(assets) != len(sources):
                raise zenodo.ReleaseError(
                    "prepared asset count differs from historical manifest"
                )
            for asset, source in zip(assets, sources, strict=True):
                expected_asset_keys = {"kind", "filename", "prepared_path", "sha256"}
                if historical.adoption is not None:
                    expected_asset_keys.update(("md5", "size"))
                if not isinstance(asset, dict) or set(asset) != expected_asset_keys:
                    raise zenodo.ReleaseError("prepared asset schema is invalid")
                if not zenodo.SHA256.fullmatch(str(asset.get("sha256", ""))):
                    raise zenodo.ReleaseError("prepared asset SHA-256 is invalid")
                expected_kind = (
                    "zenodo-adoption"
                    if historical.adoption is not None
                    else source["kind"]
                )
                if (
                    asset["kind"] != expected_kind
                    or asset["filename"] != source["filename"]
                ):
                    raise zenodo.ReleaseError(
                        "prepared asset identity differs from historical manifest"
                    )
                if asset["sha256"] != source["sha256"]:
                    raise zenodo.ReleaseError(
                        "prepared asset digest differs from historical manifest"
                    )
                if historical.adoption is not None and (
                    asset["md5"] != source["md5"] or asset["size"] != source["size"]
                ):
                    raise zenodo.ReleaseError(
                        "prepared adoption file identity differs from manifest"
                    )
                basename = Path(source["filename"]).name
                expected_path = (
                    f"canonical/{historical.family}/{historical.version}/{basename}"
                )
                if asset["prepared_path"] != expected_path:
                    raise zenodo.ReleaseError(
                        "prepared asset path differs from canonical layout"
                    )
                path = self._resolved(
                    self.plan_path.parent, asset["prepared_path"], "asset"
                )
                if (
                    not path.is_file()
                    or hashlib.sha256(path.read_bytes()).hexdigest() != asset["sha256"]
                ):
                    raise zenodo.ReleaseError(
                        f"prepared historical asset changed: {path.name}"
                    )
                if str(asset.get("filename", "")).lower().endswith(".whl"):
                    raise zenodo.ReleaseError(
                        "wheel files are forbidden in Zenodo backfill"
                    )
                if expected_kind == "git-archive":
                    key = (
                        f"{historical.family}/{historical.version}/{source['filename']}"
                    )
                    planned_git_hashes[key] = asset["sha256"]
        git_hashes = self.plan["git_archive_hashes"]
        if not isinstance(git_hashes, dict) or git_hashes != planned_git_hashes:
            raise zenodo.ReleaseError(
                "prepared Git archive hash map differs from event assets"
            )
        if not all(
            isinstance(key, str) and zenodo.SHA256.fullmatch(str(digest))
            for key, digest in git_hashes.items()
        ):
            raise zenodo.ReleaseError("prepared Git archive hash map is invalid")

    @staticmethod
    def _safe_relative(value: Any) -> Path:
        if not isinstance(value, str):
            raise zenodo.ReleaseError("prepared path must be a string")
        path = PurePosixPath(value)
        if (
            path.is_absolute()
            or not path.parts
            or any(part in ("", ".", "..") for part in path.parts)
        ):
            raise zenodo.ReleaseError(f"unsafe prepared path: {value!r}")
        return Path(*path.parts)

    @classmethod
    def _resolved(cls, base: Path, value: Any, label: str) -> Path:
        base = base.resolve()
        try:
            candidate = (base / cls._safe_relative(value)).resolve(strict=True)
        except OSError as error:
            raise zenodo.ReleaseError(
                f"prepared {label} path is missing: {value!r}"
            ) from error
        if not candidate.is_relative_to(base):
            raise zenodo.ReleaseError(
                f"prepared {label} path escapes its trusted directory"
            )
        return candidate

    def events(self) -> Iterable[tuple[dict[str, Any], HistoricalVersion]]:
        return zip(self.plan["events"], self.manifest.events(), strict=True)

    def payload(
        self,
        event: Mapping[str, Any],
        historical: HistoricalVersion,
        version_doi: str,
        concept_doi: str,
    ) -> dict[str, Path]:
        if historical.adoption is not None:
            raise zenodo.ReleaseError("adopted records do not have a generated payload")
        family = self.registry.families[historical.family]
        directory = (
            self.plan_path.parent / "payload" / historical.family / historical.version
        )
        directory.mkdir(parents=True, exist_ok=True)
        files: dict[str, Path] = {}
        for asset in event["assets"]:
            name = BackfillManifest._safe_filename(
                asset["filename"], "prepared filename"
            )
            source = self._resolved(
                self.plan_path.parent, asset["prepared_path"], "asset"
            )
            target = directory / name
            target.write_bytes(source.read_bytes())
            files[name] = target
        citation = zenodo.render_citation(
            (self.root / family.citation_file).read_text(encoding="utf-8"),
            historical.version,
            historical.publication_date,
            version_doi,
            concept_doi,
        )
        generated = {
            "CITATION.cff": citation.encode(),
            "LICENSE.md": (self.root / "LICENSE.md").read_bytes(),
            "PROVENANCE.json": zenodo._json_bytes(
                {
                    "schema": 1,
                    "historical_backfill": True,
                    "repository": self.registry.repository,
                    "family": historical.family,
                    "version": historical.version,
                    "metadata_version": event["metadata_version"],
                    "publication_date": historical.publication_date,
                    "component_versions": historical.component_versions,
                    "manifest_sha256": self.manifest.sha256,
                    "concept_doi": concept_doi,
                    "version_doi": version_doi,
                    "sources": [
                        {
                            **source,
                            "sha256": asset["sha256"],
                        }
                        for source, asset in zip(
                            historical.artifacts, event["assets"], strict=True
                        )
                    ],
                }
            ),
        }
        for name, data in generated.items():
            if name in files:
                raise zenodo.ReleaseError(
                    f"historical payload filename collision: {name}"
                )
            files[name] = directory / name
            files[name].write_bytes(data)
        checksums = "".join(
            f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {name}\n"
            for name, path in sorted(files.items())
        )
        files["SHA256SUMS"] = directory / "SHA256SUMS"
        files["SHA256SUMS"].write_text(checksums, encoding="utf-8")
        return files


class BackfillPublisher:
    def __init__(self, prepared: PreparedBackfill, client: zenodo.ZenodoClient):
        if not client.production:
            raise zenodo.ReleaseError("historical backfill is production-only")
        self.prepared = prepared
        self.client = client
        self.registry = prepared.registry

    @staticmethod
    def _md5(files: Mapping[str, Path]) -> dict[str, str]:
        return {
            name: hashlib.md5(path.read_bytes()).hexdigest()
            for name, path in files.items()
        }

    def _metadata(
        self,
        event: Mapping[str, Any],
        historical: HistoricalVersion,
        umbrella_concept: str | None,
    ) -> dict[str, Any]:
        family = self.registry.families[historical.family]
        metadata = zenodo._metadata_template(self.prepared.root, family, self.registry)
        relations = [dict(item) for item in metadata.get("related_identifiers", [])]
        additions = [
            {
                "identifier": family.marker,
                "relation": "isAlternateIdentifier",
                "resource_type": "software",
            }
        ]
        if historical.family != "gammaloop-monorepo":
            if umbrella_concept is None:
                raise zenodo.ReleaseError(
                    "umbrella concept DOI must exist before component backfill"
                )
            additions.append(
                {
                    "identifier": umbrella_concept,
                    "relation": "isPartOf",
                    "resource_type": "software",
                }
            )
        seen = {
            (item.get("identifier"), str(item.get("relation", "")).lower())
            for item in relations
        }
        for addition in additions:
            key = (addition["identifier"], addition["relation"].lower())
            if key not in seen:
                relations.append(addition)
                seen.add(key)
        metadata["version"] = event["metadata_version"]
        metadata["publication_date"] = historical.publication_date
        metadata["related_identifiers"] = relations
        return metadata

    @staticmethod
    def _assert_metadata(
        record: Mapping[str, Any], expected: Mapping[str, Any]
    ) -> None:
        zenodo.assert_stable_metadata(record, expected)

    @staticmethod
    def _lineage_state(
        lineage: Iterable[Mapping[str, Any]], target: str
    ) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
        records = [dict(item) for item in lineage]
        drafts = [item for item in records if not zenodo._published(item)]
        if len(drafts) > 1:
            raise zenodo.CollisionError("historical lineage has multiple drafts")
        matches = [
            item
            for item in records
            if zenodo._published(item) and zenodo._record_version(item) == target
        ]
        if len(matches) > 1:
            raise zenodo.CollisionError(f"historical version {target} is duplicated")
        later = [
            item
            for item in records
            if zenodo._published(item)
            and zenodo._record_version(item) is not None
            and zenodo._version_order(zenodo._record_version(item))
            > zenodo._version_order(target)
        ]
        if not matches and later:
            raise zenodo.CollisionError(
                f"historical version {target} is missing before a later version"
            )
        return (matches[0] if matches else None, drafts[0] if drafts else None)

    def _preflight(self, records: list[dict[str, Any]]) -> None:
        """Validate every immutable collision before the first remote mutation."""

        umbrella_concept: str | None = None
        incomplete_family: str | None = None
        for family_spec in self.prepared.manifest.families:
            if family_spec.name == "linnet":
                self._adopt(records, family_spec)
                continue
            family = self.registry.families[family_spec.name]
            lineage = zenodo.select_lineage(
                records, family.marker, family.concept_recid
            )
            drafts = [item for item in lineage if not zenodo._published(item)]
            if len(drafts) > 1:
                raise zenodo.CollisionError(
                    f"family {family.name}: historical lineage has multiple drafts"
                )
            published = [item for item in lineage if zenodo._published(item)]
            by_version: dict[str, dict[str, Any]] = {}
            for record in published:
                version = zenodo._record_version(record)
                if version is None:
                    raise zenodo.CollisionError(
                        f"family {family.name}: published record has no valid version"
                    )
                if version in by_version:
                    raise zenodo.CollisionError(
                        f"family {family.name}: published version {version} is duplicated"
                    )
                by_version[version] = record
            expected = [version.version for version in family_spec.versions]
            present = [version for version in expected if version in by_version]
            if present != expected[: len(present)]:
                raise zenodo.CollisionError(
                    f"family {family.name}: historical published versions are not a prefix"
                )
            complete = len(present) == len(expected)
            unexpected = [version for version in by_version if version not in expected]
            if unexpected and not complete:
                raise zenodo.CollisionError(
                    f"family {family.name}: unexpected versions precede migration completion: "
                    f"{sorted(unexpected, key=zenodo._version_order)}"
                )
            if complete and any(
                zenodo._version_order(version) <= zenodo._version_order(expected[-1])
                for version in unexpected
            ):
                raise zenodo.CollisionError(
                    f"family {family.name}: unexpected version contaminates historical lineage"
                )
            if drafts and complete:
                raise zenodo.CollisionError(
                    f"family {family.name}: completed history has a leftover draft"
                )
            if drafts:
                previous = present[-1] if present else None
                target = expected[len(present)]
                draft_version = zenodo._record_version(drafts[0])
                if draft_version not in {None, previous, target}:
                    raise zenodo.CollisionError(
                        f"family {family.name}: unrelated draft {draft_version} blocks {target}"
                    )
            if incomplete_family is not None and lineage:
                raise zenodo.CollisionError(
                    f"family {family.name}: lineage exists after incomplete {incomplete_family}"
                )
            if not complete and incomplete_family is None:
                incomplete_family = family.name

            concepts = {str(record.get("conceptrecid")) for record in lineage}
            concepts.discard("None")
            if len(concepts) > 1:
                raise zenodo.CollisionError(
                    f"family {family.name}: versions span concepts"
                )
            if family_spec.name == "gammaloop-monorepo" and lineage:
                umbrella_concept = self.client.concept_doi(lineage[0])
            elif (
                family_spec.name != "gammaloop-monorepo"
                and lineage
                and umbrella_concept is None
            ):
                raise zenodo.CollisionError(
                    f"family {family.name}: component lineage exists before umbrella concept"
                )

            historical_by_version = {
                item.version: item for item in family_spec.versions
            }
            event_by_version = {
                historical.version: event
                for event, historical in self.prepared.events()
                if historical.family == family.name
            }
            for version in present:
                record = by_version[version]
                historical = historical_by_version[version]
                event = event_by_version[version]
                concept_doi = self.client.concept_doi(record)
                files = self.prepared.payload(
                    event,
                    historical,
                    self.client.version_doi(record),
                    concept_doi,
                )
                zenodo.assert_same_version_files(record, version, self._md5(files))
                self._assert_metadata(
                    record, self._metadata(event, historical, umbrella_concept)
                )

    def publish(self) -> dict[str, Any]:
        records = self.client.list_depositions()
        self._preflight(records)
        report: dict[str, Any] = {}
        umbrella_concept: str | None = None
        for family_spec in self.prepared.manifest.families:
            family = self.registry.families[family_spec.name]
            if family_spec.name == "linnet":
                self._adopt(records, family_spec)
                adoption = family_spec.versions[0].adoption
                assert adoption is not None
                report[family_spec.name] = {
                    "adopted": True,
                    "concept_doi": adoption["concept_doi"],
                    "concept_recid": adoption["concept_recid"],
                    "versions": [
                        {
                            "version": "0.17.0",
                            "record_id": adoption["record_id"],
                            "version_doi": adoption["record_doi"],
                        }
                    ],
                }
                continue
            lineage = zenodo.select_lineage(
                records, family.marker, family.concept_recid
            )
            published_versions = []
            family_concept: str | None = None
            family_concept_recid: str | None = None
            for event, historical in (
                pair
                for pair in self.prepared.events()
                if pair[1].family == family_spec.name
            ):
                target = historical.version
                existing, draft = self._lineage_state(lineage, target)
                if existing is not None:
                    record = existing
                else:
                    published = [item for item in lineage if zenodo._published(item)]
                    if draft is None:
                        if published:
                            latest = max(
                                published,
                                key=lambda item: zenodo._version_order(
                                    zenodo._record_version(item) or "0.0.0"
                                ),
                            )
                            draft = self.client.new_version(latest)
                        else:
                            if family.concept_recid is not None:
                                raise zenodo.CollisionError(
                                    f"configured concept for {family.name} is not accessible"
                                )
                            seed = self._metadata(event, historical, umbrella_concept)
                            seed["prereserve_doi"] = True
                            draft = self.client.create(seed)
                    draft_version = zenodo._record_version(draft)
                    previous_versions = {
                        zenodo._record_version(item)
                        for item in lineage
                        if zenodo._published(item)
                    }
                    if draft_version not in {target, *previous_versions, None}:
                        raise zenodo.CollisionError(
                            f"unrelated draft {draft_version} blocks historical version {target}"
                        )
                    seed = self._metadata(event, historical, umbrella_concept)
                    seed["prereserve_doi"] = True
                    record = self.client.update(draft["id"], seed)
                concept_doi = self.client.concept_doi(record)
                version_doi = self.client.version_doi(record)
                if family_spec.name == "gammaloop-monorepo":
                    if umbrella_concept is not None and umbrella_concept != concept_doi:
                        raise zenodo.CollisionError(
                            "umbrella versions span multiple concepts"
                        )
                    umbrella_concept = concept_doi
                if family_concept is not None and family_concept != concept_doi:
                    raise zenodo.CollisionError(
                        f"family {family.name}: versions span concepts"
                    )
                family_concept = concept_doi
                family_concept_recid = str(record.get("conceptrecid"))
                files = self.prepared.payload(
                    event, historical, version_doi, concept_doi
                )
                expected_files = self._md5(files)
                expected_metadata = self._metadata(event, historical, umbrella_concept)
                if zenodo._published(record):
                    zenodo.assert_same_version_files(record, target, expected_files)
                    self._assert_metadata(record, expected_metadata)
                else:
                    record = self.client.update(record["id"], expected_metadata)
                    uploaded = self.client.replace_files(record, files)
                    if zenodo._file_identity(uploaded) != expected_files:
                        raise zenodo.ReleaseError(
                            f"family {family.name}: uploaded files differ"
                        )
                    record = self.client.publish(record["id"])
                    zenodo.assert_same_version_files(record, target, expected_files)
                    self._assert_metadata(record, expected_metadata)
                    lineage = [
                        item
                        for item in lineage
                        if str(item.get("id")) != str(record.get("id"))
                    ]
                    lineage.append(record)
                    records.append(record)
                published_versions.append(
                    {
                        "version": target,
                        "record_id": str(record.get("record_id") or record.get("id")),
                        "version_doi": self.client.version_doi(record),
                    }
                )
            if any(not zenodo._published(item) for item in lineage):
                raise zenodo.CollisionError(
                    f"family {family.name}: leftover historical draft"
                )
            report[family_spec.name] = {
                "adopted": False,
                "concept_doi": family_concept,
                "concept_recid": family_concept_recid,
                "versions": published_versions,
            }
        return report

    def _adopt(self, records: list[dict[str, Any]], family: HistoricalFamily) -> None:
        historical = family.versions[0]
        assert historical.adoption is not None
        adoption = historical.adoption
        configured = self.registry.families[family.name].concept_recid
        if configured != adoption["concept_recid"]:
            raise zenodo.CollisionError(
                "Linnet registry concept differs from adopted concept"
            )
        lineage = zenodo.select_lineage(
            records, self.registry.families[family.name].marker, configured
        )
        matches = [
            item for item in lineage if str(item.get("id")) == adoption["record_id"]
        ]
        if len(matches) != 1 or not zenodo._published(matches[0]):
            raise zenodo.ReleaseError(
                "adopted Linnet production record is not accessible"
            )
        record = matches[0]
        BackfillPreparer._verify_adoption_record(record, adoption)
        if (
            self.client.concept_doi(record) != adoption["concept_doi"]
            or self.client.version_doi(record) != adoption["record_doi"]
            or zenodo._record_version(record) != historical.version
        ):
            raise zenodo.CollisionError("adopted Linnet deposition identity differs")
        expected = {item["filename"]: item["md5"] for item in adoption["files"]}
        if zenodo._file_identity(record) != expected:
            raise zenodo.CollisionError("adopted Linnet deposition files differ")
        if any(not zenodo._published(item) for item in lineage):
            raise zenodo.CollisionError("adopted Linnet lineage has a leftover draft")

    def audit(self) -> dict[str, Any]:
        records = self.client.list_depositions()
        report: dict[str, Any] = {}
        umbrella_family = self.registry.families["gammaloop-monorepo"]
        umbrella_lineage = zenodo.select_lineage(
            records, umbrella_family.marker, umbrella_family.concept_recid
        )
        umbrella_published = [
            item for item in umbrella_lineage if zenodo._published(item)
        ]
        if not umbrella_published:
            raise zenodo.ReleaseError("umbrella historical lineage is missing")
        umbrella_concept = self.client.concept_doi(umbrella_published[0])
        for family_spec in self.prepared.manifest.families:
            if family_spec.name == "linnet":
                self._adopt(records, family_spec)
                adoption = family_spec.versions[0].adoption
                assert adoption is not None
                report[family_spec.name] = {
                    "status": "adopted",
                    "baseline": "0.17.0",
                    "concept_doi": adoption["concept_doi"],
                    "concept_recid": adoption["concept_recid"],
                    "versions": [
                        {
                            "version": "0.17.0",
                            "record_id": adoption["record_id"],
                            "version_doi": adoption["record_doi"],
                        }
                    ],
                }
                continue
            family = self.registry.families[family_spec.name]
            lineage = zenodo.select_lineage(
                records, family.marker, family.concept_recid
            )
            if any(not zenodo._published(item) for item in lineage):
                raise zenodo.CollisionError(
                    f"family {family.name}: leftover historical draft"
                )
            version_records = []
            for event, historical in (
                pair
                for pair in self.prepared.events()
                if pair[1].family == family_spec.name
            ):
                matches = [
                    item
                    for item in lineage
                    if zenodo._published(item)
                    and zenodo._record_version(item) == historical.version
                ]
                if len(matches) != 1:
                    raise zenodo.ReleaseError(
                        f"family {family.name}: expected one published {historical.version}"
                    )
                record = matches[0]
                files = self.prepared.payload(
                    event,
                    historical,
                    self.client.version_doi(record),
                    self.client.concept_doi(record),
                )
                zenodo.assert_same_version_files(
                    record, historical.version, self._md5(files)
                )
                self._assert_metadata(
                    record, self._metadata(event, historical, umbrella_concept)
                )
                version_records.append(
                    {
                        "version": historical.version,
                        "record_id": str(record.get("record_id") or record.get("id")),
                        "version_doi": self.client.version_doi(record),
                    }
                )
            report[family_spec.name] = {
                "status": "published",
                "baseline": family_spec.versions[-1].version,
                "concept_doi": self.client.concept_doi(lineage[0]),
                "concept_recid": str(lineage[0].get("conceptrecid")),
                "versions": version_records,
            }
        return report


def _load(args: argparse.Namespace) -> tuple[Path, zenodo.Registry, BackfillManifest]:
    root = Path(args.repository_root).resolve()
    registry = zenodo.load_registry(root / args.registry)
    manifest = BackfillManifest.load(root / args.manifest, registry)
    return root, registry, manifest


def _prepare(args: argparse.Namespace) -> None:
    root, registry, manifest = _load(args)
    BackfillPreparer(root, manifest, registry).prepare(
        Path(args.output_dir), Path(args.plan)
    )


def _operation(args: argparse.Namespace, audit: bool) -> None:
    if args.confirmation != CONFIRMATION:
        raise zenodo.ReleaseError(f"confirmation must exactly equal {CONFIRMATION!r}")
    root, registry, manifest = _load(args)
    prepared = PreparedBackfill(root, Path(args.plan), manifest, registry)
    client = zenodo.ZenodoClient(args.token, zenodo.PRODUCTION_API)
    publisher = BackfillPublisher(prepared, client)
    report = publisher.audit() if audit else publisher.publish()
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        report_path = Path(args.report).resolve()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--repository-root", default=".")
    common.add_argument("--registry", default=".zenodo/records.json")
    common.add_argument("--manifest", default=".zenodo/backfill.json")
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare", parents=[common])
    prepare.add_argument("--output-dir", required=True)
    prepare.add_argument("--plan", required=True)
    prepare.set_defaults(handler=_prepare)
    for name, audit in (("publish", False), ("audit", True)):
        operation = commands.add_parser(name, parents=[common])
        operation.add_argument("--plan", required=True)
        operation.add_argument("--confirmation", required=True)
        operation.add_argument("--report")
        operation.add_argument("--token", default=os.environ.get("ZENODO_TOKEN"))
        operation.set_defaults(
            handler=lambda args, audit=audit: _operation(args, audit)
        )
    return parser


def main(argv: list[str] | None = None) -> int:
    try:
        arguments = _parser().parse_args(argv)
        arguments.handler(arguments)
    except (zenodo.ReleaseError, OSError, ValueError) as error:
        print(f"zenodo-backfill: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
