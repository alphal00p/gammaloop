#!/usr/bin/env python3
"""Summarize 3Drep evaluate manifests in a workspace."""

from __future__ import annotations

import argparse
import json
import math
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any

try:
    from prettytable import PrettyTable
except ImportError as exc:  # pragma: no cover - user environment guard
    raise SystemExit(
        "Missing dependency: prettytable. Install it with `python3 -m pip install prettytable`."
    ) from exc


RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"
RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
MAGENTA = "\033[35m"
CYAN = "\033[36m"
DETAILS_MAX_WIDTH = 40


def color(text: object, code: str, enabled: bool) -> str:
    value = str(text)
    if not enabled:
        return value
    return f"{code}{value}{RESET}"


def human_sig(value: float | None, scale: float, suffix: str) -> str:
    if value is None or not math.isfinite(value):
        return "-"
    return f"{value * scale:.3g} {suffix}"


def with_ratio(
    value: float | None,
    baseline: float | None,
    scale: float,
    suffix: str,
    use_color: bool,
) -> str:
    formatted = human_sig(value, scale, suffix)
    if (
        value is None
        or baseline is None
        or not math.isfinite(value)
        or not math.isfinite(baseline)
        or baseline == 0.0
    ):
        return formatted
    ratio = value / baseline
    if not math.isfinite(ratio):
        return formatted
    ratio_color = GREEN if ratio <= 1.05 else YELLOW if ratio <= 2.0 else RED
    return f"{formatted} ({color(f'{ratio:.3g}x', ratio_color, use_color)})"


def normalize_representation(manifest: dict[str, Any], path: Path) -> str:
    evaluation = manifest.get("evaluation", {})
    representation = evaluation.get("representation")
    if representation:
        return str(representation).upper()
    lowered_parts = [part.lower() for part in path.parts]
    if any(part.startswith("ltd_") for part in lowered_parts):
        return "LTD"
    if any(part.startswith("cff_") for part in lowered_parts):
        return "CFF"
    return "N/A"


def normalize_precision(manifest: dict[str, Any]) -> str:
    evaluation = manifest.get("evaluation", {})
    settings = manifest.get("settings", {})
    return str(evaluation.get("runtime_precision") or settings.get("runtime_precision") or "-")


def normalize_evaluator_mode(manifest: dict[str, Any]) -> str:
    settings = manifest.get("settings", {})
    evaluation = manifest.get("evaluation", {})
    backend = str(
        settings.get("evaluator_backend") or evaluation.get("evaluator_backend") or "eager"
    )
    force_eager = bool(settings.get("force_eager", False))
    compiled_available = bool(settings.get("compiled_backend_available", False))
    if force_eager or backend == "eager" or not compiled_available:
        return "eager"
    if backend in {"assembly", "symjit", "cpp", "c++"}:
        normalized = "c++" if backend == "cpp" else backend
        return f"compiled {normalized}"
    return backend


def numerator_kind(manifest: dict[str, Any], path: Path) -> str:
    if manifest.get("numerator_only"):
        return "only_numerator"
    if any("no_numerator" in part.lower() for part in path.parts):
        return "no_numerator"
    return "full"


def sort_index(value: str, order: list[str]) -> tuple[int, str]:
    try:
        return (order.index(value), value)
    except ValueError:
        return (len(order), value)


def graph_sort_key(name: str) -> tuple[int, str]:
    preferred = ["Box", "DoubleBox", "TripleBox", "QuadrupleBox"]
    return sort_index(name, preferred)


def detail_fields(manifest: dict[str, Any], path: Path, workspace: Path) -> dict[str, str]:
    settings = manifest.get("settings", {})
    rel_path = path.relative_to(workspace) if path.is_relative_to(workspace) else path
    return {
        "process": str(manifest.get("process_id", "-")),
        "integrand": str(manifest.get("integrand_name", "-")),
        "graph_id": str(manifest.get("graph_id", "-")),
        "norm": str(settings.get("numerator_samples_normalization", "-")),
        "method": str(settings.get("evaluator_method", "-")),
        "seed": str(settings.get("seed", "-")),
        "scale": str(settings.get("scale", "-")),
        "path": str(rel_path.parent),
    }


def wrap_text(text: str) -> str:
    if not text:
        return text
    return "\n".join(
        textwrap.wrap(
            text,
            width=DETAILS_MAX_WIDTH,
            break_long_words=True,
            break_on_hyphens=False,
        )
    )


@dataclass(frozen=True)
class Entry:
    graph_name: str
    representation: str
    numerator_kind: str
    precision: str
    evaluator_mode: str
    build_seconds: float | None
    sample_seconds: float | None
    details: dict[str, str]
    manifest_basename: str

    @property
    def category(self) -> tuple[str, str, str, str, str]:
        return (
            self.graph_name,
            self.representation,
            self.numerator_kind,
            self.precision,
            self.evaluator_mode,
        )

    @property
    def baseline_key(self) -> tuple[str, str, str]:
        return (self.graph_name, self.representation, self.numerator_kind)

    @property
    def baseline_priority(self) -> tuple[int, tuple[Any, ...]]:
        precision_ok = self.precision == "Double"
        candidates = [
            precision_ok and self.evaluator_mode == "compiled assembly",
            precision_ok and self.evaluator_mode == "eager",
        ]
        for index, matches in enumerate(candidates):
            if matches:
                return (index, self.sort_key)
        return (len(candidates), self.sort_key)

    @property
    def sort_key(self) -> tuple[Any, ...]:
        return (
            graph_sort_key(self.graph_name),
            sort_index(self.representation, ["LTD", "CFF", "N/A"]),
            sort_index(self.numerator_kind, ["no_numerator", "only_numerator", "full"]),
            sort_index(self.precision, ["Double", "Quad", "ArbPrec"]),
            sort_index(
                self.evaluator_mode,
                ["eager", "compiled assembly", "compiled symjit"],
            ),
            tuple(self.details.items()),
            self.manifest_basename,
        )

    @property
    def dedup_key(self) -> tuple[Any, ...]:
        return (
            self.category,
            tuple(self.details.items()),
        )


def read_manifest(path: Path, workspace: Path) -> Entry | None:
    try:
        manifest = json.loads(path.read_text())
    except Exception as exc:
        print(f"warning: could not read {path}: {exc}", file=sys.stderr)
        return None

    if "evaluation" not in manifest or "graph_name" not in manifest:
        return None

    evaluation = manifest.get("evaluation", {})
    profile = evaluation.get("profile") or {}
    sample_seconds = (
        profile.get("timing_per_sample_seconds")
        or evaluation.get("sample_evaluation_timing_seconds")
    )
    build_seconds = evaluation.get("evaluator_build_timing_seconds")
    graph_name = str(manifest.get("graph_name") or "-")
    return Entry(
        graph_name=graph_name,
        representation=normalize_representation(manifest, path),
        numerator_kind=numerator_kind(manifest, path),
        precision=normalize_precision(manifest),
        evaluator_mode=normalize_evaluator_mode(manifest),
        build_seconds=float(build_seconds) if build_seconds is not None else None,
        sample_seconds=float(sample_seconds) if sample_seconds is not None else None,
        details=detail_fields(manifest, path, workspace),
        manifest_basename=path.name,
    )


def deduplicate_entries(entries: list[Entry]) -> list[Entry]:
    selected: dict[tuple[Any, ...], Entry] = {}
    for entry in sorted(entries, key=lambda item: item.sort_key):
        existing = selected.get(entry.dedup_key)
        if existing is None or (
            existing.manifest_basename == "evaluate_manifest.json"
            and entry.manifest_basename != "evaluate_manifest.json"
        ):
            selected[entry.dedup_key] = entry
    return sorted(selected.values(), key=lambda item: item.sort_key)


def collapse_duplicate_categories(entries: list[Entry]) -> list[Entry]:
    selected: dict[tuple[str, str, str, str, str], Entry] = {}
    for entry in sorted(entries, key=lambda item: item.sort_key):
        selected.setdefault(entry.category, entry)
    return sorted(selected.values(), key=lambda item: item.sort_key)


def display_sort_key(entry: Entry, baselines: dict[tuple[str, str, str], Entry]) -> tuple[Any, ...]:
    baseline = baselines.get(entry.baseline_key)
    return (
        graph_sort_key(entry.graph_name),
        sort_index(entry.representation, ["LTD", "CFF", "N/A"]),
        sort_index(entry.numerator_kind, ["no_numerator", "only_numerator", "full"]),
        0 if entry == baseline else 1,
        sort_index(entry.precision, ["Double", "Quad", "ArbPrec"]),
        sort_index(
            entry.evaluator_mode,
            ["compiled assembly", "eager", "compiled symjit"],
        ),
        tuple(entry.details.items()),
        entry.manifest_basename,
    )


def collect_entries(workspace: Path) -> list[Entry]:
    entries = []
    for path in sorted(workspace.rglob("*.json")):
        entry = read_manifest(path, workspace)
        if entry is not None:
            entries.append(entry)
    return deduplicate_entries(entries)


def baseline_map(entries: list[Entry]) -> dict[tuple[str, str, str], Entry]:
    baselines = {}
    for baseline_key in sorted({entry.baseline_key for entry in entries}):
        candidates = [
            entry
            for entry in entries
            if entry.baseline_key == baseline_key
        ]
        ranked = [
            entry
            for entry in sorted(candidates, key=lambda candidate: candidate.baseline_priority)
            if entry.baseline_priority[0] < 2
        ]
        if ranked:
            baselines[baseline_key] = ranked[0]
    return baselines


def format_details(entry: Entry, reference: Entry, use_color: bool) -> str:
    if entry == reference:
        return ""
    lines = []
    for key, value in entry.details.items():
        if value == reference.details.get(key):
            continue
        wrapped = wrap_text(f"{key}={value}")
        lines.append(color(wrapped, RED, use_color))
    return "\n".join(lines)


def render_table(entries: list[Entry], use_color: bool, show_duplicates: bool) -> str:
    baselines = baseline_map(entries)
    reference_by_category: dict[tuple[str, str, str, str, str], Entry] = {}
    for entry in entries:
        reference_by_category.setdefault(entry.category, entry)
    duplicate_categories = {
        category
        for category, reference in reference_by_category.items()
        if any(other.category == category and other != reference for other in entries)
    }
    display_entries = entries if show_duplicates else collapse_duplicate_categories(entries)
    display_entries = sorted(display_entries, key=lambda entry: display_sort_key(entry, baselines))
    include_details = show_duplicates and bool(duplicate_categories)

    table = PrettyTable()
    fields = [
        "graph_name",
        "representation",
        "kind",
        "precision",
        "evaluator mode",
        "evaluator build time",
        "evaluation / sample",
    ]
    if include_details:
        fields.append("details")
    fields.append("manifest")
    table.field_names = fields
    table.align = "l"
    table.title = "3Drep Evaluate Manifest Summary"

    previous_graph = None
    for entry in display_entries:
        if previous_graph is not None and previous_graph != entry.graph_name:
            table.add_divider()
        previous_graph = entry.graph_name

        baseline = baselines.get(entry.baseline_key)
        build_baseline = baseline.build_seconds if baseline else None
        sample_baseline = baseline.sample_seconds if baseline else None
        row = [
            color(entry.graph_name, BLUE, use_color),
            color(entry.representation, MAGENTA if entry.representation == "CFF" else CYAN, use_color),
            color(entry.numerator_kind, YELLOW, use_color),
            color(entry.precision, GREEN, use_color),
            color(entry.evaluator_mode, BOLD, use_color),
            with_ratio(entry.build_seconds, build_baseline, 1.0, "s", use_color),
            with_ratio(entry.sample_seconds, sample_baseline, 1_000_000.0, "µs", use_color),
        ]
        if include_details:
            reference = reference_by_category[entry.category]
            details = format_details(entry, reference, use_color)
            row.append(details)
        manifest_color = GREEN if entry == baseline else BLUE
        row.append(color(entry.manifest_basename, manifest_color, use_color))
        table.add_row(row)

    return table.get_string()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize all 3Drep evaluate_manifest.json files under a workspace."
    )
    parser.add_argument("workspace", type=Path, help="3Drep workspace directory")
    parser.add_argument(
        "--no-color",
        action="store_true",
        help="Disable ANSI colors in the table.",
    )
    parser.add_argument(
        "--show-duplicates",
        action="store_true",
        help="Show entries that share the same displayed category but differ in details.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = args.workspace.expanduser().resolve()
    if not workspace.exists() or not workspace.is_dir():
        print(f"error: workspace does not exist or is not a directory: {workspace}", file=sys.stderr)
        return 2

    entries = collect_entries(workspace)
    if not entries:
        print(f"No evaluate_manifest.json files found under {workspace}")
        return 1

    print(
        render_table(
            entries,
            use_color=not args.no_color,
            show_duplicates=args.show_duplicates,
        )
    )
    print(
        color(
            f"\nCollected {len(entries)} evaluate manifest(s) from {workspace}",
            DIM,
            not args.no_color,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
