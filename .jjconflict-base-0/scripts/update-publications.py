#!/usr/bin/env python3
"""Refresh the checked-in publication cache from INSPIRE HEP.

The list of scholarly authors is owned by ``docs/portal.toml``.  Every
``[[people]]`` entry with ``publications = true`` must provide a stable portal
``id`` and an ``inspire_bai``.  The updater queries each BAI separately so a
publication shared by multiple members can be deduplicated while retaining all
matching people for client-side filtering.

Only public bibliographic metadata is cached.  The update is assembled fully
in memory and written atomically, so network or validation failures cannot
leave a partial file behind.
"""

from __future__ import annotations

import argparse
import html
import json
import os
import re
import sys
import tempfile
import time
from dataclasses import dataclass
from datetime import date, datetime, timezone
from html.parser import HTMLParser
from http.client import HTTPException
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urljoin, urlsplit
from urllib.request import Request, urlopen

try:
    import tomllib
except ModuleNotFoundError as exc:  # pragma: no cover - depends on interpreter
    raise SystemExit("update-publications.py requires Python 3.11 or newer") from exc


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PORTAL = REPOSITORY_ROOT / "docs/portal.toml"
DEFAULT_OUTPUT = REPOSITORY_ROOT / "docs/data/publications.json"
INSPIRE_ORIGIN = "https://inspirehep.net"
INSPIRE_API = f"{INSPIRE_ORIGIN}/api/literature"
INSPIRE_AUTHORS_API = f"{INSPIRE_ORIGIN}/api/authors"
SCHEMA_VERSION = 1
PAGE_SIZE = 250
U64_MAX = 2**64 - 1
USER_AGENT = (
    "alphal00p-publications/1 "
    "(+https://github.com/alphal00p/gammaloop; public documentation cache)"
)
FIELDS = (
    "control_number",
    "titles.title",
    "authors.full_name",
    "earliest_date",
    "citation_count",
    "document_type",
    "publication_info",
    "arxiv_eprints.value",
    "dois.value",
)
BAI_RE = re.compile(r"^[A-Za-z0-9_.-]+\.\d+$")
PORTAL_ID_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
SUPERSCRIPT = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")
SUBSCRIPT = str.maketrans(
    {
        "0": "₀",
        "1": "₁",
        "2": "₂",
        "3": "₃",
        "4": "₄",
        "5": "₅",
        "6": "₆",
        "7": "₇",
        "8": "₈",
        "9": "₉",
        "+": "₊",
        "-": "₋",
        "a": "ₐ",
        "e": "ₑ",
        "h": "ₕ",
        "i": "ᵢ",
        "j": "ⱼ",
        "k": "ₖ",
        "l": "ₗ",
        "m": "ₘ",
        "n": "ₙ",
        "o": "ₒ",
        "p": "ₚ",
        "r": "ᵣ",
        "s": "ₛ",
        "t": "ₜ",
        "u": "ᵤ",
        "v": "ᵥ",
        "x": "ₓ",
        "H": "ₕ",
    }
)


class UpdateError(RuntimeError):
    """A validation or remote-data error that must not replace the cache."""


@dataclass(frozen=True)
class ScholarlyPerson:
    id: str
    name: str
    inspire_bai: str
    inspire_recid: int


class MarkupTextExtractor(HTMLParser):
    """Collect visible text from the small MathML fragments in INSPIRE titles."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.parts: list[str] = []

    def handle_data(self, data: str) -> None:
        self.parts.append(data)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Refresh docs/data/publications.json from INSPIRE author BAIs."
    )
    parser.add_argument(
        "--portal",
        type=Path,
        default=DEFAULT_PORTAL,
        help=f"people metadata TOML (default: {DEFAULT_PORTAL})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"generated JSON cache (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="timeout in seconds for each INSPIRE request (default: 30)",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="retries for rate limits and transient failures (default: 3)",
    )
    parser.add_argument(
        "--updated",
        help="override the cache date (YYYY-MM-DD), primarily for reproducible builds",
    )
    parser.add_argument(
        "--allow-shrink",
        action="store_true",
        help="allow cached records, authors, or person associations to disappear",
    )
    parser.add_argument(
        "--keep-existing-on-error",
        action="store_true",
        help="leave an existing cache untouched and exit successfully on update errors",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the generated cache instead of writing it",
    )
    args = parser.parse_args(argv)
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    if args.retries < 0:
        parser.error("--retries cannot be negative")
    if args.updated:
        try:
            parsed_updated = date.fromisoformat(args.updated)
        except ValueError:
            parser.error("--updated must be a real calendar date in YYYY-MM-DD form")
        if parsed_updated.isoformat() != args.updated:
            parser.error("--updated must use zero-padded YYYY-MM-DD form")
    return args


def require_nonempty_string(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise UpdateError(f"{context} must be a non-empty string")
    return value.strip()


def load_scholarly_people(path: Path) -> list[ScholarlyPerson]:
    try:
        with path.open("rb") as source:
            portal = tomllib.load(source)
    except (OSError, UnicodeError, tomllib.TOMLDecodeError) as exc:
        raise UpdateError(f"failed to read {path}: {exc}") from exc

    entries = portal.get("people")
    if not isinstance(entries, list):
        raise UpdateError(f"{path} must define one or more [[people]] entries")

    people: list[ScholarlyPerson] = []
    all_ids: set[str] = set()
    all_bais: set[str] = set()
    for index, entry in enumerate(entries, start=1):
        context = f"{path} [[people]] entry {index}"
        if not isinstance(entry, dict):
            raise UpdateError(f"{context} must be a table")

        enabled = entry.get("publications", False)
        if not isinstance(enabled, bool):
            raise UpdateError(f"{context} publications must be true or false")
        if not enabled:
            continue

        person_id = require_nonempty_string(entry.get("id"), f"{context} id")
        name = require_nonempty_string(entry.get("name"), f"{context} name")
        bai = require_nonempty_string(
            entry.get("inspire_bai"), f"{context} inspire_bai"
        )
        recid = parse_unsigned(
            entry.get("inspire_recid"), f"{context} inspire_recid"
        )
        if not PORTAL_ID_RE.fullmatch(person_id):
            raise UpdateError(
                f"{context} id {person_id!r} must be a lowercase route segment"
            )
        if not BAI_RE.fullmatch(bai):
            raise UpdateError(
                f"{context} inspire_bai {bai!r} does not look like an INSPIRE BAI"
            )
        if person_id in all_ids:
            raise UpdateError(f"duplicate scholarly person id {person_id!r}")
        if bai in all_bais:
            raise UpdateError(f"duplicate INSPIRE BAI {bai!r}")
        all_ids.add(person_id)
        all_bais.add(bai)
        people.append(ScholarlyPerson(person_id, name, bai, recid))

    if not people:
        raise UpdateError(
            f"{path} enables no publication authors; add id, publications = true, "
            "and inspire_bai to the relevant [[people]] entries"
        )
    return people


def query_url(bai: str) -> str:
    query = urlencode(
        {
            "sort": "mostrecent",
            "size": PAGE_SIZE,
            "page": 1,
            "q": f"a {bai}",
            "fields": ",".join(FIELDS),
        }
    )
    return f"{INSPIRE_API}?{query}"


def aggregate_api_url(people: list[ScholarlyPerson]) -> str:
    expression = " or ".join(f"a {person.inspire_bai}" for person in people)
    return f"{INSPIRE_API}?{urlencode({'sort': 'mostrecent', 'q': f'({expression})'})}"


def validate_api_url(url: str) -> str:
    try:
        parsed = urlsplit(url)
    except ValueError as exc:
        raise UpdateError(f"invalid INSPIRE API URL: {url}") from exc
    if parsed.scheme != "https" or parsed.hostname != "inspirehep.net":
        raise UpdateError(f"refusing unexpected INSPIRE pagination URL: {url}")
    if not (
        parsed.path.startswith("/api/literature")
        or parsed.path.startswith("/api/authors/")
    ):
        raise UpdateError(f"refusing unexpected INSPIRE API path: {url}")
    return url


def retry_delay(error: HTTPError | None, attempt: int) -> float:
    if error is not None:
        retry_after = error.headers.get("Retry-After")
        if retry_after:
            try:
                return max(0.0, min(float(retry_after), 30.0))
            except ValueError:
                pass
    return min(2.0**attempt, 15.0)


def fetch_json(url: str, *, timeout: float, retries: int) -> dict[str, Any]:
    validate_api_url(url)
    last_error: Exception | None = None
    for attempt in range(retries + 1):
        try:
            request = Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": USER_AGENT,
                },
            )
            with urlopen(request, timeout=timeout) as response:
                payload = json.load(response)
            if not isinstance(payload, dict):
                raise UpdateError(f"INSPIRE returned a non-object response for {url}")
            return payload
        except HTTPError as exc:
            last_error = exc
            transient = exc.code == 429 or 500 <= exc.code < 600
            if not transient or attempt == retries:
                break
            delay = retry_delay(exc, attempt)
        except (
            URLError,
            TimeoutError,
            HTTPException,
            UnicodeError,
            json.JSONDecodeError,
        ) as exc:
            last_error = exc
            if attempt == retries:
                break
            delay = retry_delay(None, attempt)
        except ValueError as exc:
            raise UpdateError(f"invalid INSPIRE request or response for {url}: {exc}") from exc
        print(
            f"warning: transient INSPIRE error; retrying in {delay:g}s: {last_error}",
            file=sys.stderr,
        )
        time.sleep(delay)

    raise UpdateError(f"failed to fetch {url}: {last_error}")


def clean_text(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    normalized = " ".join(value.split())
    return normalized or None


def strip_markup(value: str) -> str:
    if "<" not in value:
        return html.unescape(value)
    parser = MarkupTextExtractor()
    try:
        parser.feed(value)
        parser.close()
    except ValueError as exc:
        raise UpdateError(f"failed to parse markup in INSPIRE title: {exc}") from exc
    return html.unescape("".join(parser.parts))


def script_text(value: str, table: dict[int, str]) -> str:
    translated = value.translate(table)
    return translated if len(translated) == len(value) else value


def normalize_tex_fragment(fragment: str) -> str:
    fragment = re.sub(
        r"\\(?:overline|bar)\{([^{}]+)\}",
        lambda match: f"{match.group(1)}\u0305",
        fragment,
    )
    commands = {
        r"\alpha": "α",
        r"\gamma": "γ",
        r"\ln": "ln",
        r"\pm": "±",
        r"\sim": "∼",
    }
    for command, replacement in commands.items():
        fragment = fragment.replace(command, replacement)
    fragment = re.sub(
        r"\^\{([^{}]+)\}",
        lambda match: script_text(match.group(1), SUPERSCRIPT),
        fragment,
    )
    fragment = re.sub(
        r"_\{([^{}]+)\}",
        lambda match: script_text(match.group(1), SUBSCRIPT),
        fragment,
    )
    fragment = re.sub(
        r"\^([0-9+\-*∗])",
        lambda match: script_text(match.group(1), SUPERSCRIPT),
        fragment,
    )
    fragment = re.sub(
        r"_([0-9A-Za-z])",
        lambda match: script_text(match.group(1), SUBSCRIPT),
        fragment,
    )
    return fragment.replace("{", "").replace("}", "")


def clean_title(value: Any) -> str | None:
    if not isinstance(value, str):
        return None
    title = strip_markup(value)
    parts = re.split(r"(\$[^$]*\$)", title)
    title = "".join(
        normalize_tex_fragment(part[1:-1])
        if len(part) >= 2 and part.startswith("$") and part.endswith("$")
        else part
        for part in parts
    )
    return clean_text(title)


def first_value(entries: Any, key: str) -> str | None:
    if not isinstance(entries, list):
        return None
    for entry in entries:
        if isinstance(entry, dict):
            value = clean_text(entry.get(key))
            if value:
                return value
    return None


def first_title(entries: Any) -> str | None:
    if not isinstance(entries, list):
        return None
    for entry in entries:
        if isinstance(entry, dict):
            title = clean_title(entry.get("title"))
            if title:
                return title
    return None


def parse_unsigned(value: Any, context: str, *, default: int | None = None) -> int:
    if value is None and default is not None:
        return default
    if not isinstance(value, int) or isinstance(value, bool):
        raise UpdateError(f"{context} must be a non-negative integer")
    if value < 0 or value > U64_MAX:
        raise UpdateError(f"{context} must be a non-negative integer")
    return value


def parse_partial_date(value: str) -> tuple[str, int] | None:
    candidate = value.strip()
    try:
        if re.fullmatch(r"\d{4}", candidate):
            parsed = date(int(candidate), 1, 1)
        elif re.fullmatch(r"\d{4}-\d{2}", candidate):
            year, month = (int(part) for part in candidate.split("-"))
            parsed = date(year, month, 1)
        elif re.fullmatch(r"\d{4}-\d{2}-\d{2}", candidate):
            parsed = date.fromisoformat(candidate)
        else:
            return None
    except ValueError:
        return None
    return candidate, parsed.year


def publication_date(metadata: dict[str, Any], hit: dict[str, Any], record_id: int) -> tuple[str, int]:
    candidates: list[Any] = [metadata.get("earliest_date")]
    publication_info = metadata.get("publication_info")
    if isinstance(publication_info, list):
        candidates.extend(
            entry.get("year") for entry in publication_info if isinstance(entry, dict)
        )
    candidates.append(hit.get("created"))

    for candidate in candidates:
        if isinstance(candidate, int):
            candidate = str(candidate)
        if not isinstance(candidate, str):
            continue
        parsed = parse_partial_date(candidate.strip()[:10])
        if parsed:
            return parsed
    raise UpdateError(f"INSPIRE literature {record_id} has no usable publication date")


def publication_authors(metadata: dict[str, Any], record_id: int) -> list[str]:
    raw_authors = metadata.get("authors")
    if not isinstance(raw_authors, list):
        raise UpdateError(f"INSPIRE literature {record_id} has no author list")
    authors: list[str] = []
    seen: set[str] = set()
    for entry in raw_authors:
        if not isinstance(entry, dict):
            continue
        name = clean_text(entry.get("full_name"))
        if name and name not in seen:
            authors.append(name)
            seen.add(name)
    if not authors:
        raise UpdateError(f"INSPIRE literature {record_id} has an empty author list")
    return authors


def publication_types(metadata: dict[str, Any]) -> list[str]:
    raw_types = metadata.get("document_type")
    if not isinstance(raw_types, list):
        return ["other"]
    types = sorted(
        {
            normalized
            for value in raw_types
            if (normalized := clean_text(value)) is not None
        }
    )
    return types or ["other"]


def publication_venue(metadata: dict[str, Any]) -> str | None:
    publication_info = metadata.get("publication_info")
    if not isinstance(publication_info, list):
        return None
    for entry in publication_info:
        if not isinstance(entry, dict):
            continue
        journal = clean_text(entry.get("journal_title"))
        if not journal:
            continue
        parts = [journal]
        volume = clean_text(entry.get("journal_volume"))
        if volume:
            parts.append(volume)
        year = entry.get("year")
        if isinstance(year, int) and not isinstance(year, bool):
            parts.append(f"({year})")
        article = clean_text(entry.get("artid")) or clean_text(entry.get("page_start"))
        if article:
            parts.append(article)
        return " ".join(parts)
    return None


def normalize_hit(hit: Any) -> dict[str, Any]:
    if not isinstance(hit, dict):
        raise UpdateError("INSPIRE search returned a non-object literature hit")
    metadata = hit.get("metadata")
    if not isinstance(metadata, dict):
        raise UpdateError("INSPIRE literature hit has no metadata object")

    record_id = parse_unsigned(
        metadata.get("control_number", hit.get("id")), "INSPIRE control number"
    )
    if record_id == 0:
        raise UpdateError("INSPIRE control number cannot be zero")
    title = first_title(metadata.get("titles"))
    if not title:
        raise UpdateError(f"INSPIRE literature {record_id} has no title")
    date, year = publication_date(metadata, hit, record_id)
    if year > 65535:
        raise UpdateError(f"INSPIRE literature {record_id} year does not fit u16")

    normalized: dict[str, Any] = {
        "id": record_id,
        "title": title,
        "date": date,
        "year": year,
        "authors": publication_authors(metadata, record_id),
        "people": [],
    }
    venue = publication_venue(metadata)
    if venue:
        normalized["venue"] = venue
    doi = first_value(metadata.get("dois"), "value")
    if doi:
        normalized["doi"] = doi
    arxiv = first_value(metadata.get("arxiv_eprints"), "value")
    if arxiv:
        normalized["arxiv"] = arxiv
    normalized.update(
        {
            "citations": parse_unsigned(
                metadata.get("citation_count"),
                f"INSPIRE literature {record_id} citation count",
                default=0,
            ),
            "types": publication_types(metadata),
            "url": f"{INSPIRE_ORIGIN}/literature/{record_id}",
            "bibtex_url": f"{INSPIRE_API}/{record_id}?format=bibtex",
        }
    )
    return normalized


def validate_person_identity(
    person: ScholarlyPerson, *, timeout: float, retries: int
) -> None:
    url = f"{INSPIRE_AUTHORS_API}/{person.inspire_recid}"
    payload = fetch_json(url, timeout=timeout, retries=retries)
    metadata = payload.get("metadata")
    if not isinstance(metadata, dict):
        raise UpdateError(
            f"INSPIRE author {person.inspire_recid} has no metadata object"
        )
    control_number = parse_unsigned(
        metadata.get("control_number"),
        f"INSPIRE author control number for {person.name}",
    )
    if control_number != person.inspire_recid:
        raise UpdateError(
            f"INSPIRE author record {person.inspire_recid} returned control number "
            f"{control_number}"
        )
    raw_ids = metadata.get("ids")
    if not isinstance(raw_ids, list):
        raise UpdateError(f"INSPIRE author {person.inspire_recid} has no identifier list")
    bais = {
        value
        for entry in raw_ids
        if isinstance(entry, dict)
        and clean_text(entry.get("schema")) == "INSPIRE BAI"
        and (value := clean_text(entry.get("value"))) is not None
    }
    if person.inspire_bai not in bais:
        found = ", ".join(sorted(bais)) or "none"
        raise UpdateError(
            f"portal BAI {person.inspire_bai!r} does not belong to INSPIRE author "
            f"{person.inspire_recid} ({person.name}); record BAIs: {found}"
        )


def fetch_person_publications(
    person: ScholarlyPerson, *, timeout: float, retries: int
) -> list[dict[str, Any]]:
    url: str | None = query_url(person.inspire_bai)
    visited: set[str] = set()
    records: dict[int, dict[str, Any]] = {}
    expected_total: int | None = None

    while url:
        url = validate_api_url(url)
        if url in visited:
            raise UpdateError(f"INSPIRE pagination loop while querying {person.inspire_bai}")
        visited.add(url)
        payload = fetch_json(url, timeout=timeout, retries=retries)
        hits_container = payload.get("hits")
        if not isinstance(hits_container, dict):
            raise UpdateError(f"INSPIRE query for {person.inspire_bai} has no hits object")
        total = parse_unsigned(
            hits_container.get("total"), f"INSPIRE total for {person.inspire_bai}"
        )
        if expected_total is None:
            expected_total = total
        elif total != expected_total:
            raise UpdateError(
                f"INSPIRE total changed during pagination for {person.inspire_bai}: "
                f"{expected_total} to {total}"
            )

        hits = hits_container.get("hits")
        if not isinstance(hits, list):
            raise UpdateError(f"INSPIRE query for {person.inspire_bai} has no hit list")
        for hit in hits:
            record = normalize_hit(hit)
            records[record["id"]] = record

        links = payload.get("links")
        next_url = links.get("next") if isinstance(links, dict) else None
        if next_url is not None and not isinstance(next_url, str):
            raise UpdateError(f"INSPIRE next link for {person.inspire_bai} is not a URL")
        if next_url:
            try:
                url = urljoin(url, next_url)
            except ValueError as exc:
                raise UpdateError(
                    f"invalid INSPIRE next link for {person.inspire_bai}: {next_url}"
                ) from exc
        else:
            url = None

    if expected_total is None or expected_total == 0:
        raise UpdateError(
            f"INSPIRE BAI {person.inspire_bai!r} returned no publications; "
            "verify docs/portal.toml before replacing the cache"
        )
    if len(records) != expected_total:
        raise UpdateError(
            f"INSPIRE query for {person.inspire_bai} returned {len(records)} unique records "
            f"but advertised {expected_total}"
        )
    print(
        f"{person.name}: {len(records)} INSPIRE publication(s)", file=sys.stderr
    )
    return list(records.values())


def merge_publication(
    publications: dict[int, dict[str, Any]],
    incoming: dict[str, Any],
    person_id: str,
) -> None:
    record_id = incoming["id"]
    existing = publications.get(record_id)
    if existing is None:
        incoming["people"] = [person_id]
        publications[record_id] = incoming
        return

    if person_id not in existing["people"]:
        existing["people"].append(person_id)
    existing_stable = {
        key: value
        for key, value in existing.items()
        if key not in {"people", "citations"}
    }
    incoming_stable = {
        key: value
        for key, value in incoming.items()
        if key not in {"people", "citations"}
    }
    if existing_stable != incoming_stable:
        raise UpdateError(
            f"INSPIRE literature {record_id} changed between author queries; "
            "keeping the previous cache rather than constructing a hybrid record"
        )
    # Citation counts are explicitly volatile and do not affect record identity.
    existing["citations"] = max(existing["citations"], incoming["citations"])


def generate_cache(
    people: list[ScholarlyPerson], *, timeout: float, retries: int
) -> dict[str, Any]:
    publications: dict[int, dict[str, Any]] = {}
    person_order = {person.id: index for index, person in enumerate(people)}
    for person in people:
        validate_person_identity(person, timeout=timeout, retries=retries)
        for publication in fetch_person_publications(
            person, timeout=timeout, retries=retries
        ):
            merge_publication(publications, publication, person.id)

    ordered_publications = list(publications.values())
    for publication in ordered_publications:
        publication["people"].sort(key=person_order.__getitem__)
    ordered_publications.sort(
        key=lambda publication: (publication["date"], publication["id"]),
        reverse=True,
    )
    return {
        "schema": SCHEMA_VERSION,
        "source": INSPIRE_ORIGIN,
        "updated": "",  # Filled after comparison with the existing cache.
        "api_url": aggregate_api_url(people),
        "authors": [
            {
                "id": person.id,
                "name": person.name,
                "inspire_bai": person.inspire_bai,
            }
            for person in people
        ],
        "publications": ordered_publications,
    }


def load_existing(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        with path.open(encoding="utf-8") as source:
            existing = json.load(source)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise UpdateError(f"failed to read existing cache {path}: {exc}") from exc
    if not isinstance(existing, dict) or not isinstance(
        existing.get("publications"), list
    ):
        raise UpdateError(f"existing cache {path} has an invalid shape")
    return existing


def publication_ids(cache: dict[str, Any], context: str) -> set[int]:
    ids: set[int] = set()
    for index, publication in enumerate(cache.get("publications", []), start=1):
        if not isinstance(publication, dict):
            raise UpdateError(f"{context} publication {index} is not an object")
        record_id = parse_unsigned(publication.get("id"), f"{context} publication id")
        if record_id in ids:
            raise UpdateError(f"{context} contains duplicate publication {record_id}")
        ids.add(record_id)
    return ids


def cache_authors(cache: dict[str, Any], context: str) -> dict[str, str]:
    raw_authors = cache.get("authors")
    if not isinstance(raw_authors, list):
        raise UpdateError(f"{context} has no author list")
    authors: dict[str, str] = {}
    for index, author in enumerate(raw_authors, start=1):
        if not isinstance(author, dict):
            raise UpdateError(f"{context} author {index} is not an object")
        person_id = require_nonempty_string(
            author.get("id"), f"{context} author {index} id"
        )
        bai = require_nonempty_string(
            author.get("inspire_bai"), f"{context} author {person_id} inspire_bai"
        )
        if person_id in authors:
            raise UpdateError(f"{context} contains duplicate author {person_id}")
        authors[person_id] = bai
    return authors


def publication_people(cache: dict[str, Any], context: str) -> dict[int, set[str]]:
    associations: dict[int, set[str]] = {}
    for index, publication in enumerate(cache.get("publications", []), start=1):
        if not isinstance(publication, dict):
            raise UpdateError(f"{context} publication {index} is not an object")
        record_id = parse_unsigned(publication.get("id"), f"{context} publication id")
        raw_people = publication.get("people")
        if not isinstance(raw_people, list):
            raise UpdateError(f"{context} publication {record_id} has no people list")
        people = {
            require_nonempty_string(
                person_id, f"{context} publication {record_id} person id"
            )
            for person_id in raw_people
        }
        if len(people) != len(raw_people):
            raise UpdateError(
                f"{context} publication {record_id} contains duplicate people"
            )
        associations[record_id] = people
    return associations


def guard_against_shrink(
    existing: dict[str, Any] | None,
    generated: dict[str, Any],
    *,
    allow_shrink: bool,
) -> None:
    if existing is None or allow_shrink:
        return
    old_ids = publication_ids(existing, "existing cache")
    new_ids = publication_ids(generated, "generated cache")
    missing = sorted(old_ids - new_ids)
    if missing:
        preview = ", ".join(str(record_id) for record_id in missing[:10])
        suffix = " …" if len(missing) > 10 else ""
        raise UpdateError(
            f"refusing suspicious cache shrink: {len(missing)} existing publication(s) "
            f"disappeared ({preview}{suffix}); verify the BAI queries or rerun with "
            "--allow-shrink after review"
        )
    old_authors = cache_authors(existing, "existing cache")
    new_authors = cache_authors(generated, "generated cache")
    removed_or_changed_authors = sorted(
        person_id
        for person_id, bai in old_authors.items()
        if new_authors.get(person_id) != bai
    )
    if removed_or_changed_authors:
        raise UpdateError(
            "refusing suspicious cache shrink: scholarly author metadata disappeared "
            f"or changed for {', '.join(removed_or_changed_authors)}; rerun with "
            "--allow-shrink after review"
        )
    old_people = publication_people(existing, "existing cache")
    new_people = publication_people(generated, "generated cache")
    lost_matches = [
        (record_id, sorted(people - new_people.get(record_id, set())))
        for record_id, people in old_people.items()
        if people - new_people.get(record_id, set())
    ]
    if lost_matches:
        preview = ", ".join(
            f"{record_id}: {','.join(people)}"
            for record_id, people in lost_matches[:5]
        )
        suffix = " …" if len(lost_matches) > 5 else ""
        raise UpdateError(
            "refusing suspicious cache shrink: publication/person associations "
            f"disappeared ({preview}{suffix}); rerun with --allow-shrink after review"
        )


def payload_without_updated(cache: dict[str, Any]) -> dict[str, Any]:
    comparable = dict(cache)
    comparable.pop("updated", None)
    return comparable


def choose_updated(
    generated: dict[str, Any],
    existing: dict[str, Any] | None,
    override: str | None,
) -> str:
    if override:
        return override
    if existing is not None and payload_without_updated(existing) == payload_without_updated(
        generated
    ):
        prior = existing.get("updated")
        if isinstance(prior, str):
            try:
                if date.fromisoformat(prior).isoformat() == prior:
                    return prior
            except ValueError:
                pass
    return datetime.now(timezone.utc).date().isoformat()


def serialize(cache: dict[str, Any]) -> str:
    return json.dumps(cache, ensure_ascii=False, indent=2) + "\n"


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent, text=True
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as target:
            target.write(content)
            target.flush()
            os.fsync(target.fileno())
        os.chmod(temporary, 0o644)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def update(args: argparse.Namespace) -> int:
    people = load_scholarly_people(args.portal)
    existing = load_existing(args.output)
    generated = generate_cache(people, timeout=args.timeout, retries=args.retries)
    guard_against_shrink(
        existing, generated, allow_shrink=args.allow_shrink
    )
    generated["updated"] = choose_updated(generated, existing, args.updated)
    content = serialize(generated)
    if args.dry_run:
        sys.stdout.write(content)
    else:
        atomic_write(args.output, content)
        print(
            f"updated {args.output} with {len(generated['publications'])} "
            "unique publication(s)",
            file=sys.stderr,
        )
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        return update(args)
    except (UpdateError, OSError) as exc:
        if args.keep_existing_on_error and args.output.is_file():
            print(
                f"warning: publication update failed; keeping {args.output}: {exc}",
                file=sys.stderr,
            )
            return 0
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
