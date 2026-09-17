#!/usr/bin/env python3
"""Compare tth@NNLO integration workspaces, including GL134, GL297 and GL638."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from colorama import Fore, Style, init
from prettytable import PrettyTable, TableStyle

HERE = Path(__file__).resolve().parent


def _defaults(workspace_root: Path | None) -> dict[str, dict[str, Path]]:
    roots = (
        [workspace_root]
        if workspace_root is not None
        else [Path.cwd(), HERE, HERE.parent / "NNLO_experiment"]
    )
    groups: dict[str, dict[str, Path]] = {}
    seen = set()
    for root in roots:
        for directory in (
            root,
            root / "workspaces",
            *sorted(root.glob("*/workspaces")),
        ):
            directory = directory.resolve()
            if directory in seen:
                continue
            seen.add(directory)
            local: dict[str, dict[str, Path]] = {}
            for manifest in sorted(directory.glob("*/manifest.json")):
                match = re.fullmatch(
                    r"(GL\d+)_(advanced_sampling(?:_with_optimized_lmbs|_max_weight)?|optimized_lmbs)(?:_workspace)?",
                    manifest.parent.name,
                    re.IGNORECASE,
                )
                if match is None:
                    continue  # Explicit comparison names exclude fixed-grid pilots and backups.
                graph, strategy = match[1].upper(), match[2].lower()
                paths = local.setdefault(graph, {})
                if strategy in paths:
                    raise ValueError(
                        f"Ambiguous {graph} {strategy}: {paths[strategy]} and {manifest.parent}. "
                        "Use explicit --advanced-workspace / --optimized-workspace overrides."
                    )
                paths[strategy] = manifest.parent
            for graph, paths in local.items():
                groups.setdefault(graph, paths)  # Keep a comparison in one directory.
    return groups


def latest_snapshot(workspace: Path) -> tuple[int, dict[str, Any]] | None:
    files = list((workspace / "results").glob("integration_result_iter_*.json"))
    if not files:
        return None
    path = max(files, key=lambda p: int(p.stem.rsplit("_", 1)[-1]))
    with path.open(encoding="utf-8") as stream:
        document = json.load(stream)
    slots = document.get("slots", [])
    if not slots:
        raise ValueError(f"snapshot {path} contains no integration slots")
    return int(path.stem.rsplit("_", 1)[-1]), slots[0]


def maximum(slot: dict[str, Any], component: str, absolute: bool = False) -> float:
    entries = (
        slot["absolute"]["max_weight_info"] if absolute else slot["max_weight_info"]
    )
    values = [
        abs(float(item["max_eval"]))
        for item in entries
        if item["component"] == component
    ]
    return max(values, default=0.0)


def format_value(value: float, error: float | None = None) -> str:
    if error is None:
        return f"{value:.5e}"
    return f"{value:.5e} ± {error:.5e}"


def metrics(slot: dict[str, Any]) -> dict[str, str]:
    signed = slot["integral"]
    absolute = slot["absolute"]["integral"]
    stats = slot["integration_statistics"]
    return {
        "samples": f"{signed['neval']:,}",
        "Re": format_value(signed["result"]["re"], signed["error"]["re"]),
        "|Re|": format_value(absolute["result"]["re"], absolute["error"]["re"]),
        "Im": format_value(signed["result"]["im"], signed["error"]["im"]),
        "|Im|": format_value(absolute["result"]["im"], absolute["error"]["im"]),
        "max |Re|": format_value(maximum(slot, "re")),
        "max |Im|": format_value(maximum(slot, "im")),
        "max abs |Re|": format_value(maximum(slot, "re", True)),
        "max abs |Im|": format_value(maximum(slot, "im", True)),
        "Re impact": f"{slot['table_results'][0]['max_weight_impact']:.3%}",
        "Im impact": f"{slot['table_results'][1]['max_weight_impact']:.3%}",
        "|Re| impact": f"{slot['absolute']['table_results'][0]['max_weight_impact']:.3%}",
        "|Im| impact": f"{slot['absolute']['table_results'][1]['max_weight_impact']:.3%}",
        "f64": f"{stats['f64_percentage']:.3f}%",
        "f128": f"{stats['f128_percentage']:.3f}%",
        "ArbPrec": f"{stats['arb_percentage']:.3f}%",
        "unstable": f"{stats['nan_or_unstable_percentage']:.3f}%",
        "time/sample": f"{stats['average_total_time_seconds']:.4f} s",
        "sampling time": f"{stats['average_parameterization_time_seconds']:.4f} s",
        "evaluator time": f"{stats['average_evaluator_time_seconds']:.4f} s",
        "integrand overhead": f"{stats['average_integrand_time_seconds'] - stats['average_evaluator_time_seconds']:.4f} s",
    }


_ITERATION_RE = re.compile(
    r"Iteration #\s*(\d+)\s*\(\s*completed\s*\).*?# samples total =\s*([0-9.]+)([KMG]?)"
)
_OBSERVABLE_RE = re.compile(
    r"│\s*All\s+(re|\|re\||im|\|im\|)\s+│\s*([+-]?\d+(?:\.\d+)?)\((\d+(?:\.\d+)?)\)e([+-]?\d+)\s+.*?[0-9.]+%\s+│\s*[0-9.]+\s+([0-9.]+e[+-]?\d+)"
)
_ABS_OBSERVABLE_RE = re.compile(
    r"│\s*(?:All\s+)?(\|re\||\|im\|)\s+│\s*([+-]?\d+(?:\.\d+)?)\((\d+(?:\.\d+)?)\)e([+-]?\d+)"
)
_MAX_RE = re.compile(
    r"│\s*epem_a_tth@NNLO\s+(re|im)\s+\[([+-])\]\s*│\s*([+-]?\d+(?:\.\d+)?)e([+-]?\d+)"
)
_TIMING_RE = re.compile(
    r"│\s+timing\s+│\s+total\s+([0-9.]+)\s+(ns|µs|ms|s)\s+│\s+param\s+([0-9.]+)\s+(ns|µs|ms|s)\s+│\s+itg\s+([0-9.]+)\s+(ns|µs|ms|s)\s+│\s+evaluators\s+([0-9.]+)\s+(ns|µs|ms|s)"
)
_EVAL_RE = re.compile(
    r"│\s+evals\s+│\s+f64\s+([0-9.]+)%\s+│\s+f128\s+([0-9.]+)%\s+│\s+arb\s+([0-9.]+)%\s+│\s+nans\+unstable\s+([0-9.]+)%"
)


def _compact_value(value: str, uncertainty: str, exponent: str) -> tuple[float, float]:
    coefficient = float(value)
    error = float(uncertainty)
    if "." not in uncertainty:
        error /= 10 ** len(value.partition(".")[2])
    scale = 10 ** int(exponent)
    return coefficient * scale, error * scale


def live_metrics(*log_paths: Path) -> tuple[int, dict[str, str]] | None:
    log_path = max(
        (path for path in log_paths if path.exists()),
        key=lambda path: path.stat().st_mtime,
        default=None,
    )
    if log_path is None:
        return None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    iterations = list(_ITERATION_RE.finditer(text))
    if not iterations:
        return None
    iteration = iterations[-1]
    if text.rfind("Streaming integration updates disabled") > iteration.end():
        return None
    samples = (
        float(iteration.group(2))
        * {"": 1, "K": 1e3, "M": 1e6, "G": 1e9}[iteration.group(3)]
    )
    tail = text[iteration.end() :]
    values: dict[str, tuple[float, float]] = {}
    impacts: dict[str, float] = {}
    for match in _OBSERVABLE_RE.finditer(tail):
        values[match.group(1)] = _compact_value(
            match.group(2), match.group(3), match.group(4)
        )
        impacts[match.group(1)] = float(match.group(5))
    for match in _ABS_OBSERVABLE_RE.finditer(tail):
        values[match.group(1)] = _compact_value(
            match.group(2), match.group(3), match.group(4)
        )
    timing = _TIMING_RE.search(tail)
    evaluation = _EVAL_RE.search(tail)
    maxima: dict[str, float] = {}
    for match in _MAX_RE.finditer(tail):
        maxima[match.group(1)] = max(
            maxima.get(match.group(1), 0.0),
            abs(float(f"{match.group(3)}e{match.group(4)}")),
        )
    if (
        not {"re", "|re|", "im", "|im|"} <= values.keys()
        or timing is None
        or evaluation is None
    ):
        return None
    total, parameterization, integrand, evaluator = (
        float(value) * {"ns": 1e-9, "µs": 1e-6, "ms": 1e-3, "s": 1}[unit]
        for value, unit in zip(timing.groups()[::2], timing.groups()[1::2])
    )
    f64, f128, arb, unstable = map(float, evaluation.groups())
    unavailable = "n/a (not logged live)"
    return int(iteration.group(1)), {
        "samples": f"{int(samples):,}",
        "Re": format_value(*values["re"]),
        "|Re|": format_value(*values["|re|"]),
        "Im": format_value(*values["im"]),
        "|Im|": format_value(*values["|im|"]),
        "max |Re|": format_value(maxima["re"]) if "re" in maxima else unavailable,
        "max |Im|": format_value(maxima["im"]) if "im" in maxima else unavailable,
        "max abs |Re|": unavailable,
        "max abs |Im|": unavailable,
        "Re impact": f"{impacts['re']:.3%}",
        "Im impact": f"{impacts['im']:.3%}",
        "|Re| impact": unavailable,
        "|Im| impact": unavailable,
        "f64": f"{f64:.3f}%",
        "f128": f"{f128:.3f}%",
        "ArbPrec": f"{arb:.3f}%",
        "unstable": f"{unstable:.3f}%",
        "time/sample": f"{total:.4f} s",
        "sampling time": f"{parameterization:.4f} s",
        "evaluator time": f"{evaluator:.4f} s",
        "integrand overhead": f"{integrand - evaluator:.4f} s",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--graph", default="gl638", help="graph name, e.g. gl134, gl297 or gl638"
    )
    parser.add_argument(
        "--list-available-graphs",
        action="store_true",
        help="list discovered comparison workspaces and exit",
    )
    parser.add_argument(
        "--workspace-root",
        type=Path,
        help="search only this directory (or its workspaces/ children); defaults to the current directory and example directories beside this script",
    )
    parser.add_argument(
        "--live",
        action="store_true",
        help="read the latest completed summary from integration.log or *.long.log",
    )
    parser.add_argument("--advanced-workspace", type=Path)
    parser.add_argument("--optimized-workspace", type=Path)
    parser.add_argument(
        "--max-weight-workspace",
        type=Path,
        help="GL638 advanced workspace with the additional outlier-focused channel",
    )
    parser.add_argument(
        "--augmented-workspace",
        type=Path,
        help="GL638 advanced workspace with optimized-LMB channels",
    )
    args = parser.parse_args()
    init(autoreset=True)
    try:
        available = (
            {}
            if args.advanced_workspace
            and args.optimized_workspace
            and not args.list_available_graphs
            else _defaults(args.workspace_root)
        )
    except ValueError as error:
        parser.error(str(error))
    if args.list_available_graphs:
        table = PrettyTable(["Graph", "Comparison directory", "Available setups"])
        table.set_style(TableStyle.SINGLE_BORDER)
        table.align = "l"
        for graph, paths in sorted(available.items()):
            table.add_row(
                [graph, str(next(iter(paths.values())).parent), ", ".join(paths)]
            )
        print(
            table
            if available
            else "No comparison workspaces found. Use --workspace-root PATH."
        )
        return
    graph = args.graph.upper()
    paths = available.get(graph, {}).copy()
    for strategy, override in (
        ("advanced_sampling", args.advanced_workspace),
        ("optimized_lmbs", args.optimized_workspace),
        ("advanced_sampling_with_optimized_lmbs", args.augmented_workspace),
        ("advanced_sampling_max_weight", args.max_weight_workspace),
    ):
        if override is not None:
            paths[strategy] = override.resolve()
    if not paths:
        parser.error(
            f"No comparison workspaces found for {graph}; use --list-available-graphs or --workspace-root PATH."
        )
    columns = []
    for strategy, title, colour in (
        ("advanced_sampling", "Advanced sampling", Fore.CYAN),
        ("optimized_lmbs", "Optimized-LMB only", Fore.YELLOW),
        (
            "advanced_sampling_with_optimized_lmbs",
            "Advanced + optimized-LMB",
            Fore.MAGENTA,
        ),
        (
            "advanced_sampling_max_weight",
            "Advanced + max-weight channel",
            Fore.GREEN,
        ),
    ):
        workspace = paths.get(strategy)
        if workspace is None:
            continue
        result = None
        source = "iteration dump"
        if workspace is not None and args.live:
            base = (
                workspace.parent.parent
                if workspace.parent.name == "workspaces"
                else workspace.parent
            )
            state_strategy = strategy.replace("_with_", "_")
            result = live_metrics(
                workspace / "integration.log",
                base / f"gammaloop_state_{graph}_{state_strategy}" / "integration.log",
                base / "states" / f"{graph}_{strategy}" / "integration.log",
                base / f"{graph.lower()}_{strategy}.long.log",
            )
            if result:
                source = "live log"
        if result is None and workspace is not None:
            snapshot = latest_snapshot(workspace)
            if snapshot is not None:
                result = snapshot[0], metrics(snapshot[1])
        iteration, data = result or (None, None)
        if data is not None:
            data["iteration"] = str(iteration)
            data["source"] = source
        columns.append((title, data, colour, iteration))
    table = PrettyTable()
    table.set_style(TableStyle.SINGLE_BORDER)
    table.field_names = ["Observable / statistic"] + [column[0] for column in columns]
    table.align["Observable / statistic"] = "l"
    for title, _, _, _ in columns:
        table.align[title] = "r"
    mode = "live summary" if args.live else "latest iteration"
    iteration_title = "  |  ".join(
        f"{colour}{title} #{iteration or 'pending'}{Style.RESET_ALL}"
        for title, _, colour, iteration in columns
    )
    table.title = f"{graph} {mode}  {iteration_title}"
    rows = [
        ("Completed iterations", "iteration"),
        ("Data source", "source"),
        ("Samples", "samples"),
        ("Re", "Re"),
        ("|Re|", "|Re|"),
        ("Im", "Im"),
        ("|Im|", "|Im|"),
        ("Maximum |Re| weight", "max |Re|"),
        ("Maximum |Im| weight", "max |Im|"),
        ("Maximum absolute-Re weight", "max abs |Re|"),
        ("Maximum absolute-Im weight", "max abs |Im|"),
        ("Max-weight impact on Re", "Re impact"),
        ("Max-weight impact on |Re|", "|Re| impact"),
        ("Max-weight impact on Im", "Im impact"),
        ("Max-weight impact on |Im|", "|Im| impact"),
        ("f64 evaluations", "f64"),
        ("f128 evaluations", "f128"),
        ("ArbPrec evaluations", "ArbPrec"),
        ("Unstable or NaN", "unstable"),
        ("Time per sample", "time/sample"),
        ("  Sampling maps", "sampling time"),
        ("  Integrand overhead", "integrand overhead"),
        ("  Evaluators", "evaluator time"),
    ]
    for label, key in rows:
        if key == "unstable":
            values = []
            for _, data, _, _ in columns:
                value = data[key] if data else "pending (no completed iteration)"
                values.append(
                    (Fore.GREEN if value.startswith("0.000") else Fore.RED)
                    + value
                    + Style.RESET_ALL
                )
        else:
            values = []
            for _, data, colour, _ in columns:
                value = data[key] if data else "pending (no completed iteration)"
                values.append(colour + value + Style.RESET_ALL)
        table.add_row([label] + values)
    print(table)
    for strategy, workspace in paths.items():
        print(f"{Style.DIM}{strategy}: {workspace}{Style.RESET_ALL}")


if __name__ == "__main__":
    main()
