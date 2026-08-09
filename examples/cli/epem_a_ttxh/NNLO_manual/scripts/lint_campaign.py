#!/usr/bin/env python3
"""Validate the reproducibility contract of the manual NNLO graph campaign."""

from __future__ import annotations

import argparse
import re
import tomllib
from pathlib import Path

from collect_results import GRAPH_IDS
from normalize_dot import normalize_text

PROCESS_SPEC = (
    "e+ e- > t t~ h | e+ e- g t t~ h d d~ ghG ghG~ a QCD^2==4 QED^2==6 [{{4}} QCD=2]"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "root",
        type=Path,
        nargs="?",
        default=Path(__file__).resolve().parent.parent,
    )
    parser.add_argument("--allow-incomplete", action="store_true")
    return parser.parse_args()


def check_dot(path: Path) -> list[str]:
    errors = []
    raw = path.read_text()
    _, removed = normalize_text(raw)
    if any(removed.values()):
        errors.append(f"{path}: contains generated/drawing attributes {removed}")
    if len(re.findall(r"\blmb_id\s*=", raw)) != 4:
        errors.append(f"{path}: does not contain exactly four lmb_id attributes")
    if "threshold_counterterms" in raw:
        if not re.search(r'threshold_counterterms\s*=\s*"\r?\n', raw):
            errors.append(
                f"{path}: threshold_counterterms is not a multiline DOT string"
            )
        if re.search(r"(?m)^\s*cuts\s*=\s*\[\s*\]\s*$", raw):
            errors.append(f"{path}: explicitly stores an empty threshold document")
    return errors


def command_text(card: dict) -> str:
    return "\n".join(
        command
        for block in card.get("command_blocks", [])
        for command in block.get("commands", [])
    )


def check_card(path: Path, graph_id: str) -> list[str]:
    errors = []
    try:
        card = tomllib.loads(path.read_text())
    except tomllib.TOMLDecodeError as error:
        return [f"{path}: invalid TOML: {error}"]
    commands = command_text(card)
    expected_import = (
        f"import graphs ./examples/cli/epem_a_ttxh/NNLO_manual/graphs/{graph_id}.dot "
        f"--process-spec '{PROCESS_SPEC}' -p epem_a_tth -i NNLO -o"
    )
    if expected_import not in commands:
        errors.append(f"{path}: missing exact process-aware graph import")
    for prohibited in ("--orientation-id", "force_cuts", "orientation_pattern"):
        if prohibited in commands:
            errors.append(f"{path}: command contains prohibited {prohibited}")
    settings = card.get("cli_settings", {}).get("global", {})
    generation = settings.get("generation", {})
    if "force_cuts" in generation or "orientation_pattern" in generation:
        errors.append(f"{path}: filters cuts or orientations in generation settings")
    evaluator = generation.get("evaluator", {})
    expected_evaluator = {
        "do_algebra": False,
        "iterative_orientation_optimization": False,
        "store_atom": True,
        "compile": False,
        "summed": False,
        "summed_function_map": False,
    }
    for key, expected in expected_evaluator.items():
        if evaluator.get(key) != expected:
            errors.append(f"{path}: generation.evaluator.{key} != {expected!r}")
    cores = settings.get("n_cores", {})
    if cores.get("generate") != 4 or cores.get("compile") != 0:
        errors.append(f"{path}: requires generate=4 and compile=0 cores")
    runtime_general = card.get("default_runtime_settings", {}).get("general", {})
    if runtime_general.get("evaluator_method") != "SingleParametric":
        errors.append(f"{path}: runtime evaluator is not SingleParametric")
    state_folder = card.get("cli_settings", {}).get("state", {}).get("folder")
    expected_state = f"./examples/cli/epem_a_ttxh/NNLO_manual/states/state_{graph_id}"
    if state_folder != expected_state:
        errors.append(f"{path}: state folder is not {expected_state}")
    block_names = {block.get("name") for block in card.get("command_blocks", [])}
    if "all_approaches" not in block_names:
        errors.append(f"{path}: has no all_approaches block")
    if not {"integrate", "integration"} & block_names:
        errors.append(f"{path}: has no integration block")
    return errors


def main() -> int:
    args = parse_args()
    expected = set(GRAPH_IDS)
    found_dots = {path.stem for path in (args.root / "graphs").glob("GL*.dot")}
    found_cards = {
        path.stem.removeprefix("run_")
        for path in (args.root / "run_cards").glob("run_GL*.toml")
    }
    found_results = {path.stem for path in (args.root / "results").glob("GL*.md")}
    errors = []
    for label, found in (
        ("DOT", found_dots),
        ("card", found_cards),
        ("result", found_results),
    ):
        extras = sorted(found - expected)
        missing = sorted(expected - found)
        if extras:
            errors.append(f"unexpected {label} IDs: {' '.join(extras)}")
        if missing and not args.allow_incomplete:
            errors.append(f"missing {label} IDs: {' '.join(missing)}")
    for graph_id in sorted(found_dots & expected):
        errors.extend(check_dot(args.root / "graphs" / f"{graph_id}.dot"))
    for graph_id in sorted(found_cards & expected):
        errors.extend(
            check_card(args.root / "run_cards" / f"run_{graph_id}.toml", graph_id)
        )
    summary = args.root / "results" / "summary.md"
    if (
        not summary.exists()
        or sum(line.startswith("| GL") for line in summary.read_text().splitlines())
        != 71
    ):
        errors.append(f"{summary}: must contain exactly 71 graph rows")
    if errors:
        print("\n".join(errors))
        return 1
    print(
        f"validated {len(found_dots)} DOTs, {len(found_cards)} cards, "
        f"and {len(found_results)} result files"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
