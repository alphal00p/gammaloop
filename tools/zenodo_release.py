#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tomllib
import urllib.error
import urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SUPPORTED_CRATES = {
    "spenso": ROOT / "crates" / "spenso",
    "idenso": ROOT / "crates" / "idenso",
    "vakint": ROOT / "crates" / "vakint",
    "gammalooprs": ROOT / "crates" / "gammalooprs",
    "linnet": ROOT / "crates" / "linnet",
}
TAG_RE = re.compile(
    r"^(?P<crate>[a-z0-9][a-z0-9_-]*)-v(?P<version>[0-9][0-9A-Za-z.+-]*)$"
)


class CffError(RuntimeError):
    pass


def parse_scalar(raw: str) -> str:
    raw = raw.strip()
    if not raw:
        return ""
    if raw.startswith('"') and raw.endswith('"'):
        return json.loads(raw)
    if raw.startswith("'") and raw.endswith("'"):
        return raw[1:-1].replace("''", "'")
    return raw


def parse_cff(path: Path) -> dict[str, object]:
    data: dict[str, object] = {}
    lines = path.read_text(encoding="utf-8").splitlines()
    i = 0

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            i += 1
            continue
        if line.startswith(" "):
            raise CffError(f"Unexpected indentation in {path}: {line!r}")

        key, sep, rest = line.partition(":")
        if not sep:
            raise CffError(f"Malformed line in {path}: {line!r}")
        key = key.strip()
        rest = rest.strip()
        i += 1

        if rest:
            data[key] = parse_scalar(rest)
            continue

        if key == "authors":
            authors: list[dict[str, str]] = []
            current: dict[str, str] | None = None
            while i < len(lines):
                nested = lines[i]
                nested_stripped = nested.strip()
                if not nested_stripped or nested_stripped.startswith("#"):
                    i += 1
                    continue
                if not nested.startswith("  "):
                    break
                if nested.startswith("  - "):
                    if current is not None:
                        authors.append(current)
                    current = {}
                    tail = nested[4:].strip()
                    if tail:
                        nkey, nsep, nrest = tail.partition(":")
                        if not nsep:
                            raise CffError(
                                f"Malformed author entry in {path}: {nested!r}"
                            )
                        current[nkey.strip()] = parse_scalar(nrest)
                elif nested.startswith("    "):
                    if current is None:
                        raise CffError(
                            f"Unexpected author continuation in {path}: {nested!r}"
                        )
                    nkey, nsep, nrest = nested[4:].partition(":")
                    if not nsep:
                        raise CffError(f"Malformed author field in {path}: {nested!r}")
                    current[nkey.strip()] = parse_scalar(nrest)
                else:
                    raise CffError(
                        f"Unsupported nested structure in {path}: {nested!r}"
                    )
                i += 1
            if current is not None:
                authors.append(current)
            data[key] = authors
            continue

        if key == "keywords":
            keywords: list[str] = []
            while i < len(lines):
                nested = lines[i]
                nested_stripped = nested.strip()
                if not nested_stripped or nested_stripped.startswith("#"):
                    i += 1
                    continue
                if not nested.startswith("  "):
                    break
                if not nested.startswith("  - "):
                    raise CffError(f"Malformed keyword entry in {path}: {nested!r}")
                keywords.append(parse_scalar(nested[4:]))
                i += 1
            data[key] = keywords
            continue

        raise CffError(f"Unsupported nested key in {path}: {key}")

    return data


def load_manifest(crate_dir: Path) -> dict[str, object]:
    manifest_path = crate_dir / "Cargo.toml"
    with manifest_path.open("rb") as handle:
        manifest = tomllib.load(handle)
    return manifest["package"]


def resolve_release(tag: str) -> dict[str, object]:
    match = TAG_RE.fullmatch(tag)
    if not match:
        raise SystemExit(f"Tag {tag!r} does not match '<crate>-v<version>'")

    crate = match.group("crate")
    version = match.group("version")
    crate_dir = SUPPORTED_CRATES.get(crate)
    if crate_dir is None:
        supported = ", ".join(sorted(SUPPORTED_CRATES))
        raise SystemExit(f"Unsupported crate {crate!r}; supported crates: {supported}")

    package = load_manifest(crate_dir)
    manifest_version = str(package["version"])
    if manifest_version != version:
        raise SystemExit(
            f"Tag version {version} does not match {crate_dir / 'Cargo.toml'} version {manifest_version}"
        )

    cff_path = crate_dir / "CITATION.cff"
    citation = parse_cff(cff_path)
    cff_version = citation.get("version")
    if cff_version is not None and str(cff_version) != version:
        raise SystemExit(
            f"CFF version {cff_version} does not match tag version {version} for {crate}"
        )

    package_name = str(package["name"])
    archive_path = ROOT / "target" / "package" / f"{package_name}-{version}.crate"

    return {
        "crate": crate,
        "crate_dir": str(crate_dir.relative_to(ROOT)),
        "manifest_path": str((crate_dir / "Cargo.toml").relative_to(ROOT)),
        "citation_path": str(cff_path.relative_to(ROOT)),
        "package": package_name,
        "version": version,
        "archive_path": str(archive_path.relative_to(ROOT)),
        "citation": citation,
        "manifest": package,
    }


def build_zenodo_payload(release: dict[str, object]) -> dict[str, object]:
    citation = release["citation"]
    manifest = release["manifest"]
    authors = citation.get("authors", [])
    if not authors:
        raise SystemExit(f"No authors found in {release['citation_path']}")

    creators: list[dict[str, str]] = []
    for author in authors:
        family = str(author.get("family-names", "")).strip()
        given = str(author.get("given-names", "")).strip()
        if not family or not given:
            raise SystemExit(
                f"Each author in {release['citation_path']} must have family-names and given-names"
            )
        creator = {"name": f"{family}, {given}"}
        if author.get("affiliation"):
            creator["affiliation"] = str(author["affiliation"])
        if author.get("orcid"):
            creator["orcid"] = str(author["orcid"])
        creators.append(creator)

    description = str(
        citation.get("abstract") or manifest.get("description") or release["crate"]
    )
    metadata: dict[str, object] = {
        "title": str(citation.get("title", release["crate"])),
        "upload_type": "software",
        "description": description,
        "creators": creators,
        "version": str(release["version"]),
        "publication_date": date.today().isoformat(),
        "notes": (
            f"Release tag: {release['crate']}-v{release['version']}\n"
            f"Source repository: {citation.get('repository-code', '')}"
        ).strip(),
    }
    if citation.get("license"):
        metadata["license"] = str(citation["license"])
    if citation.get("keywords"):
        metadata["keywords"] = list(citation["keywords"])
    return {"metadata": metadata}


def write_github_output(tag: str, output_path: Path) -> None:
    release = resolve_release(tag)
    with output_path.open("a", encoding="utf-8") as handle:
        for key in (
            "crate",
            "crate_dir",
            "manifest_path",
            "citation_path",
            "package",
            "version",
            "archive_path",
        ):
            handle.write(f"{key}={release[key]}\n")


def metadata_command(tag: str) -> int:
    release = resolve_release(tag)
    payload = build_zenodo_payload(release)
    json.dump(payload, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


def request_json(
    method: str, url: str, token: str, payload: dict[str, object] | None = None
) -> dict[str, object]:
    data = None
    headers = {"Authorization": f"Bearer {token}"}
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(request) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(
            f"Zenodo API request failed ({exc.code} {exc.reason}): {body}"
        ) from exc


def upload_archive(bucket_url: str, archive_path: Path, token: str) -> None:
    request = urllib.request.Request(
        f"{bucket_url.rstrip('/')}/{archive_path.name}",
        data=archive_path.read_bytes(),
        headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/octet-stream",
        },
        method="PUT",
    )
    try:
        with urllib.request.urlopen(request):
            return
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(
            f"Zenodo upload failed ({exc.code} {exc.reason}): {body}"
        ) from exc


def publish_command(tag: str, archive: Path) -> int:
    release = resolve_release(tag)
    archive_path = archive if archive.is_absolute() else ROOT / archive
    if not archive_path.exists():
        raise SystemExit(f"Archive not found: {archive_path}")

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        raise SystemExit("ZENODO_TOKEN is required")
    api_root = (os.environ.get("ZENODO_API_URL") or "https://zenodo.org/api").rstrip(
        "/"
    )

    deposition = request_json(
        "POST", f"{api_root}/deposit/depositions", token, payload={}
    )
    bucket_url = deposition["links"]["bucket"]
    deposition_id = deposition["id"]

    upload_archive(bucket_url, archive_path, token)
    payload = build_zenodo_payload(release)
    request_json(
        "PUT", f"{api_root}/deposit/depositions/{deposition_id}", token, payload=payload
    )
    published = request_json(
        "POST",
        f"{api_root}/deposit/depositions/{deposition_id}/actions/publish",
        token,
        payload={},
    )

    result = {
        "deposition_id": deposition_id,
        "record_html": published.get("links", {}).get("record_html"),
        "doi": published.get("metadata", {}).get("prereserve_doi", {}).get("doi"),
    }
    json.dump(result, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Resolve crate release metadata and publish tagged crate archives to Zenodo."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    github_output_parser = subparsers.add_parser("github-output")
    github_output_parser.add_argument("--tag", required=True)
    github_output_parser.add_argument("--output", required=True)

    metadata_parser = subparsers.add_parser("metadata")
    metadata_parser.add_argument("--tag", required=True)

    publish_parser = subparsers.add_parser("publish")
    publish_parser.add_argument("--tag", required=True)
    publish_parser.add_argument("--archive", required=True)

    args = parser.parse_args()

    if args.command == "github-output":
        write_github_output(args.tag, Path(args.output))
        return 0
    if args.command == "metadata":
        return metadata_command(args.tag)
    if args.command == "publish":
        return publish_command(args.tag, Path(args.archive))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
