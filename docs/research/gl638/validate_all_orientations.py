#!/usr/bin/env python3
"""Verify a GL638 unrestricted-generation artifact without loading its state.

This is an acceptance verifier, not a generator.  It reads the small
``generation.json`` provenance record, the saved global settings, and the
orientation prefix of the (potentially multi-gigabyte) standalone JSON.  The
standalone file is streamed only until the GL638 orientation array has been
decoded.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import tomllib
from typing import Any

EXPECTED_ORIENTATIONS = 936
ORIENTATION_WIDTH = 15
MAX_PREFIX_BYTES = 64 * 1024 * 1024
ORIENTATION_RE = re.compile(
    rb'"graph_name"\s*:\s*"GL638"\s*,\s*"orientations"\s*:\s*'
)


class ValidationError(RuntimeError):
    """An artifact failed a required acceptance check."""


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fail(message: str) -> None:
    raise ValidationError(message)


def require_file(path: Path, label: str) -> Path:
    if not path.is_file():
        fail(f"{label} does not exist or is not a file: {path}")
    return path


def read_orientations(path: Path) -> list[list[int]]:
    """Decode the GL638 orientation array from a standalone JSON prefix."""
    decoder = json.JSONDecoder()
    buffered = b""
    with path.open("rb") as stream:
        while len(buffered) <= MAX_PREFIX_BYTES:
            chunk = stream.read(1024 * 1024)
            if not chunk:
                break
            buffered += chunk
            match = ORIENTATION_RE.search(buffered)
            if match is None:
                continue
            try:
                text = buffered[match.end() :].decode("utf-8")
                value, _ = decoder.raw_decode(text)
            except (UnicodeDecodeError, json.JSONDecodeError):
                # The array may end in the next bounded chunk.
                continue
            if not isinstance(value, list):
                fail(f"GL638 orientations value is not an array in {path}")
            if not all(isinstance(row, list) for row in value):
                fail(f"GL638 orientation keys are not arrays in {path}")
            return value
    fail(
        f"could not decode the GL638 orientation array within {MAX_PREFIX_BYTES} "
        f"bytes of {path}; refusing summary-only validation"
    )


def read_orientation_file(path: Path) -> list[list[int]]:
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"cannot read orientation-key file {path}: {exc}")
    if not isinstance(value, list) or not all(isinstance(row, list) for row in value):
        fail(f"orientation-key file must contain an array of arrays: {path}")
    return value


def check_orientation_keys(
    orientations: list[list[int]],
    artifact: dict[str, Any],
    expected_count: int,
) -> dict[str, Any]:
    if len(orientations) != expected_count:
        fail(f"standalone has {len(orientations)} orientations, expected {expected_count}")
    keys = [tuple(row) for row in orientations]
    if len(set(keys)) != len(keys):
        fail("standalone orientation keys are not unique")
    if not all(len(row) == ORIENTATION_WIDTH for row in orientations):
        fail(f"orientation keys must have width {ORIENTATION_WIDTH}")
    if not all(type(value) is int for row in orientations for value in row):
        fail("orientation keys contain non-integer values")
    if not all(value in (-1, 0, 1) for row in orientations for value in row):
        fail("orientation keys contain values outside {-1, 0, 1}")
    for field in ("orientation_count", "unique_orientation_count"):
        recorded = artifact.get(field)
        if recorded != len(orientations):
            fail(f"artifact {field}={recorded!r} disagrees with standalone key count")
    return {
        "orientation_count": len(orientations),
        "unique_orientation_count": len(set(keys)),
        "orientation_width": ORIENTATION_WIDTH,
        "first_key": list(orientations[0]),
        "last_key": list(orientations[-1]),
    }


def load_toml(path: Path, label: str) -> dict[str, Any]:
    try:
        return tomllib.loads(path.read_text())
    except (OSError, tomllib.TOMLDecodeError) as exc:
        fail(f"cannot parse {label} {path}: {exc}")


def get_nested(mapping: dict[str, Any], *keys: str, label: str) -> Any:
    value: Any = mapping
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            fail(f"{label} is missing [{'.'.join(keys)}]")
        value = value[key]
    return value


def check_provenance(
    artifact: dict[str, Any], card: Path, graph: Path, binary: Path
) -> dict[str, str]:
    for key, path in (("card_sha256", card), ("graph_sha256", graph), ("binary_sha256", binary)):
        recorded = artifact.get(key)
        if not isinstance(recorded, str):
            fail(f"artifact has no {key}; provenance is incomplete")
        actual = digest(path)
        if actual != recorded:
            fail(f"stale {key}: artifact={recorded}, current={actual}, path={path}")
    return {
        "card_sha256": artifact["card_sha256"],
        "graph_sha256": artifact["graph_sha256"],
        "binary_sha256": artifact["binary_sha256"],
    }


def check_saved_settings(artifact: dict[str, Any], state: Path) -> dict[str, Any]:
    settings = require_file(state / "global_settings.toml", "saved global settings")
    actual_hash = digest(settings)
    recorded_settings = artifact.get("saved_settings", {})
    recorded = (
        recorded_settings.get("global_settings.toml")
        if isinstance(recorded_settings, dict)
        else None
    )
    if isinstance(recorded, dict) and isinstance(recorded.get("sha256"), str):
        if actual_hash != recorded["sha256"]:
            fail(
                "saved global settings hash mismatch: "
                f"artifact={recorded['sha256']}, current={actual_hash}"
            )
        settings_attestation = "artifact hash matched"
    else:
        # Older generation-summary.json files omitted this digest. The
        # structural checks below still run, but the result is summary-only
        # provenance rather than a byte-for-byte state attestation.
        settings_attestation = "artifact omitted settings hash; structural state check only"
    current = load_toml(settings, "saved global settings")
    generation = get_nested(current, "global", "generation", label="saved settings")
    uv = get_nested(generation, "uv", label="saved settings")
    thresholds = get_nested(generation, "threshold_subtraction", label="saved settings")
    if generation.get("explicit_orientation_sum_only") is not True:
        fail("generation was not configured with explicit_orientation_sum_only=true")
    orientation_pattern = generation.get("orientation_pattern")
    if orientation_pattern not in ({}, None):
        fail(f"saved orientation pattern is selected, not unrestricted: {orientation_pattern!r}")
    if uv.get("generate_integrated") is not True or uv.get("subtract_uv") is not True:
        fail("saved settings do not enable both integrated UV generation and UV subtraction")
    if uv.get("final_integrand") != "ThreeD":
        fail(f"saved local UV route is not ThreeD: {uv.get('final_integrand')!r}")
    if thresholds.get("enable_thresholds") is not True:
        fail("saved threshold subtraction is disabled")
    return {
        "explicit_orientation_sum_only": True,
        "orientation_pattern": orientation_pattern,
        "integrated_uv": True,
        "subtract_uv": True,
        "local_uv_route": uv["final_integrand"],
        "thresholds": True,
        "saved_settings_sha256": actual_hash,
        "saved_settings_attestation": settings_attestation,
    }


def check_metadata(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    require_file(path, "metadata comparison")
    value = json.loads(path.read_text())
    if value.get("variants_equal") is not True or value.get("components_equal") is not True:
        fail(f"metadata comparison is not equal: {path}")
    if value.get("generated_variants") != 19 or value.get("active_variants") != 19:
        fail(f"metadata comparison does not contain all 19 variants: {path}")
    return {key: value.get(key) for key in ("generated_variants", "active_variants")}


def check_current_commit(artifact: dict[str, Any], scope: str, repo: Path) -> str:
    recorded = artifact.get("source_commit")
    if scope == "archive":
        return "archive provenance (source commit not required)"
    if not isinstance(recorded, str) or not recorded:
        fail("current validation requires artifact source_commit; use --scope archive for old evidence")
    try:
        current = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
        ).strip()
        dirty = subprocess.check_output(
            ["git", "-C", str(repo), "status", "--porcelain", "--untracked-files=no"], text=True
        ).strip()
    except OSError as exc:
        fail(f"cannot inspect repository commit: {exc}")
    if current != recorded:
        fail(f"stale source commit: artifact={recorded}, current={current}")
    if dirty:
        fail("current validation requires a clean tracked worktree")
    return f"current checkout at {current}"


def validate(
    artifact_path: Path,
    state: Path,
    card: Path,
    graph: Path,
    binary: Path,
    standalone: Path | None,
    orientation_file: Path | None,
    metadata: Path | None,
    scope: str,
    expected_count: int,
    repo: Path,
) -> dict[str, Any]:
    require_file(artifact_path, "generation artifact")
    try:
        artifact = json.loads(artifact_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"cannot read generation artifact {artifact_path}: {exc}")
    if artifact.get("returncode") != 0 or artifact.get("timed_out") is True:
        fail("generation did not finish successfully")
    state = require_file(state / "global_settings.toml", "state global settings").parent
    card = require_file(card, "run card")
    graph = require_file(graph, "graph DOT")
    binary = require_file(binary, "gammaloop binary")
    provenance = check_provenance(artifact, card, graph, binary)
    source = check_current_commit(artifact, scope, repo)
    settings = check_saved_settings(artifact, state)
    if standalone is not None:
        standalone = require_file(standalone, "standalone archive")
        orientations = read_orientations(standalone)
        orientation_source = str(standalone)
    elif orientation_file is not None:
        orientation_file = require_file(orientation_file, "orientation-key file")
        orientations = read_orientation_file(orientation_file)
        orientation_source = str(orientation_file)
    else:
        fail("provide --standalone or --orientations; summary counts alone are not accepted")
    orientation_report = check_orientation_keys(orientations, artifact, expected_count)
    metadata_report = check_metadata(metadata)
    report = {
        "status": "ok",
        "scope": scope,
        "provenance": source,
        "artifact": str(artifact_path),
        "state": str(state),
        "orientation_source": orientation_source,
        "expected_orientation_count": expected_count,
        **orientation_report,
        **settings,
        **provenance,
    }
    if metadata_report is not None:
        report["metadata"] = metadata_report
    return report


def self_test() -> None:
    with tempfile.TemporaryDirectory(prefix="gl638-verifier-") as raw:
        root = Path(raw)
        card, graph, binary = (root / name for name in ("card.toml", "graph.dot", "binary"))
        card.write_text("[cli_settings.global.generation]\nexplicit_orientation_sum_only=true\n")
        graph.write_text("digraph GL638 {}\n")
        binary.write_bytes(b"synthetic-binary")
        state = root / "state"
        state.mkdir()
        settings = {
            "global": {
                "generation": {
                    "explicit_orientation_sum_only": True,
                    "orientation_pattern": {},
                    "uv": {"generate_integrated": True, "subtract_uv": True, "final_integrand": "ThreeD"},
                    "threshold_subtraction": {"enable_thresholds": True},
                }
            }
        }
        settings_path = state / "global_settings.toml"
        settings_path.write_text(
            "[global.generation]\nexplicit_orientation_sum_only=true\n\n"
            "[global.generation.orientation_pattern]\n\n"
            "[global.generation.uv]\ngenerate_integrated=true\nsubtract_uv=true\nfinal_integrand=\"ThreeD\"\n\n"
            "[global.generation.threshold_subtraction]\nenable_thresholds=true\n"
        )
        orientations = [[1 if (i >> bit) & 1 else -1 for bit in range(15)] for i in range(EXPECTED_ORIENTATIONS)]
        (root / "orientations.json").write_text(json.dumps(orientations))
        artifact = {
            "returncode": 0,
            "timed_out": False,
            "orientation_count": EXPECTED_ORIENTATIONS,
            "unique_orientation_count": EXPECTED_ORIENTATIONS,
            "card_sha256": digest(card),
            "graph_sha256": digest(graph),
            "binary_sha256": digest(binary),
            "saved_settings": {"global_settings.toml": {"sha256": digest(settings_path)}},
        }
        artifact_path = root / "generation.json"
        artifact_path.write_text(json.dumps(artifact))
        report = validate(
            artifact_path, state, card, graph, binary, None, root / "orientations.json", None,
            "archive", EXPECTED_ORIENTATIONS, root
        )
        assert report["status"] == "ok"
        graph.write_text("stale\n")
        try:
            validate(
                artifact_path, state, card, graph, binary, None, root / "orientations.json", None,
                "archive", EXPECTED_ORIENTATIONS, root
            )
        except ValidationError:
            return
        fail("self-test did not reject stale graph hash")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, help="generation.json provenance artifact")
    parser.add_argument("--state", type=Path, help="saved state directory")
    parser.add_argument("--card", type=Path, help="run card used by generation")
    parser.add_argument("--graph", type=Path, help="GL638 DOT used by generation")
    parser.add_argument("--binary", type=Path, help="gammaloop binary used by generation")
    parser.add_argument("--standalone", type=Path, help="standalone JSON containing source orientations")
    parser.add_argument("--orientations", type=Path, help="small extracted orientation-key JSON")
    parser.add_argument("--metadata", type=Path, help="metadata-comparison.json (optional)")
    parser.add_argument("--repo", type=Path, default=Path.cwd(), help="repository for source-commit validation")
    parser.add_argument("--expected-count", type=int, default=EXPECTED_ORIENTATIONS)
    parser.add_argument("--scope", choices=("current", "archive"), default="current")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        if args.self_test:
            self_test()
            print(json.dumps({"status": "ok", "self_test": True}))
            return 0
        required = {"artifact": args.artifact, "state": args.state, "card": args.card, "graph": args.graph, "binary": args.binary}
        missing = [name for name, value in required.items() if value is None]
        if missing:
            parser.error("missing required options: " + ", ".join("--" + name for name in missing))
        report = validate(
            args.artifact, args.state, args.card, args.graph, args.binary,
            args.standalone, args.orientations, args.metadata, args.scope,
            args.expected_count, args.repo
        )
    except ValidationError as exc:
        print(json.dumps({"status": "failed", "error": str(exc)}), file=sys.stderr)
        return 1
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
