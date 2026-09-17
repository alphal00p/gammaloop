#!/usr/bin/env python3
"""Compare completed iterations for the GL297/GL638 tth@NNLO runs."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from colorama import Fore, Style, init
from prettytable import PrettyTable, TableStyle


HERE = Path(__file__).resolve().parent


def _defaults(graph: str) -> tuple[Path, Path, Path | None, Path, Path | None, str]:
    graph = graph.upper()
    if graph not in {"GL297", "GL638"}:
        raise ValueError("--graph must be gl297 or gl638")
    stem = graph
    advanced = HERE / "workspaces" / f"{stem}_advanced_sampling_workspace"
    optimized = HERE / "workspaces" / f"{stem}_optimized_lmbs_workspace"
    advanced_state = HERE / f"gammaloop_state_{stem}_advanced_sampling"
    if graph == "GL638":
        augmented = (
            HERE
            / "workspaces"
            / "GL638_advanced_sampling_with_optimized_lmbs_workspace"
        )
        augmented_state = (
            HERE / "gammaloop_state_GL638_advanced_sampling_optimized_lmbs"
        )
    else:
        augmented = augmented_state = None
    return advanced, optimized, augmented, advanced_state, augmented_state, stem


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
    r"│\s*All\s+(re|\|re\||im|\|im\|)\s+│\s*([+-]?\d+(?:\.\d+)?)\((\d+(?:\.\d+)?)\)e([+-]\d+)\s+.*?[0-9.]+%\s+│\s*[0-9.]+\s+([0-9.]+e[+-]\d+)"
)
_ABS_OBSERVABLE_RE = re.compile(
    r"│\s*(?:All\s+)?(\|re\||\|im\|)\s+│\s*([+-]?\d+(?:\.\d+)?)\((\d+(?:\.\d+)?)\)e([+-]\d+)"
)
_MAX_RE = re.compile(
    r"│\s*epem_a_tth@NNLO\s+(re|im)\s+\[([+-])\]\s*│\s*([+-]?\d+(?:\.\d+)?)e([+-]\d+)"
)
_TIMING_RE = re.compile(
    r"│\s+timing\s+│\s+total\s+([0-9.]+)\s+ms\s+│\s+param\s+([0-9.]+)\s+ms\s+│\s+itg\s+([0-9.]+)\s+ms\s+│\s+evaluators\s+([0-9.]+)\s+ms"
)
_EVAL_RE = re.compile(
    r"│\s+evals\s+│\s+f64\s+([0-9.]+)%\s+│\s+f128\s+([0-9.]+)%\s+│\s+arb\s+([0-9.]+)%\s+│\s+nans\+unstable\s+([0-9.]+)%"
)


def _compact_value(value: str, uncertainty: str, exponent: str) -> tuple[float, float]:
    coefficient = float(value)
    error = float(uncertainty)
    if "." not in uncertainty:
        error /= 10 ** len(value.split(".", 1)[1])
    scale = 10 ** int(exponent)
    return coefficient * scale, error * scale


def live_metrics(log_path: Path) -> tuple[int, dict[str, str]] | None:
    if not log_path.exists():
        return None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    iterations = list(_ITERATION_RE.finditer(text))
    if not iterations:
        return None
    iteration = iterations[-1]
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
            maxima.get(match.group(1), 0.0), float(f"{match.group(3)}e{match.group(4)}")
        )
    if (
        not {"re", "|re|", "im", "|im|"} <= values.keys()
        or timing is None
        or evaluation is None
    ):
        return None
    total_ms, parameterization_ms, integrand_ms, evaluator_ms = map(
        float, timing.groups()
    )
    f64, f128, arb, unstable = map(float, evaluation.groups())
    unavailable = "n/a (not logged live)"
    return int(iteration.group(1)), {
        "samples": f"{int(samples):,}",
        "Re": format_value(*values["re"]),
        "|Re|": format_value(*values["|re|"]),
        "Im": format_value(*values["im"]),
        "|Im|": format_value(*values["|im|"]),
        "max |Re|": format_value(maxima.get("re", 0.0)),
        "max |Im|": format_value(maxima.get("im", 0.0)),
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
        "time/sample": f"{total_ms / 1000:.4f} s",
        "sampling time": f"{parameterization_ms / 1000:.4f} s",
        "evaluator time": f"{evaluator_ms / 1000:.4f} s",
        "integrand overhead": f"{(integrand_ms - evaluator_ms) / 1000:.4f} s",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--graph", default="gl638", help="graph selector: gl297 or gl638"
    )
    parser.add_argument(
        "--live",
        action="store_true",
        help="read the latest completed summary from integration.log",
    )
    parser.add_argument("--advanced-workspace", type=Path)
    parser.add_argument("--optimized-workspace", type=Path)
    parser.add_argument(
        "--augmented-workspace",
        type=Path,
        help="GL638 advanced workspace with optimized-LMB channels",
    )
    args = parser.parse_args()
    (
        advanced_default,
        optimized_default,
        augmented_default,
        advanced_state,
        augmented_state,
        graph,
    ) = _defaults(args.graph)
    advanced_workspace = args.advanced_workspace or advanced_default
    optimized_workspace = args.optimized_workspace or optimized_default
    augmented_workspace = args.augmented_workspace or augmented_default
    init(autoreset=True)
    if args.live:
        advanced_live = live_metrics(advanced_state / "integration.log")
        optimized_live = live_metrics(
            advanced_state.parent
            / f"gammaloop_state_{graph}_optimized_lmbs"
            / "integration.log"
        )
        augmented_live = (
            live_metrics(augmented_state / "integration.log")
            if augmented_state
            else None
        )
        advanced_iter, advanced = advanced_live or (None, None)
        optimized_iter, optimized = optimized_live or (None, None)
        augmented_iter, augmented = augmented_live or (None, None)
    else:
        advanced_snapshot = latest_snapshot(advanced_workspace)
        optimized_snapshot = latest_snapshot(optimized_workspace)
        advanced_iter = advanced_snapshot[0] if advanced_snapshot else None
        optimized_iter = optimized_snapshot[0] if optimized_snapshot else None
        advanced = metrics(advanced_snapshot[1]) if advanced_snapshot else None
        optimized = metrics(optimized_snapshot[1]) if optimized_snapshot else None
        augmented_snapshot = (
            latest_snapshot(augmented_workspace) if augmented_workspace else None
        )
        augmented_iter = augmented_snapshot[0] if augmented_snapshot else None
        augmented = metrics(augmented_snapshot[1]) if augmented_snapshot else None
    table = PrettyTable()
    table.set_style(TableStyle.SINGLE_BORDER)
    columns = [
        ("Advanced sampling", advanced, Fore.CYAN, advanced_iter),
        ("Optimized-LMB only", optimized, Fore.YELLOW, optimized_iter),
    ]
    if graph == "GL638":
        columns.append(
            ("Advanced + optimized-LMB", augmented, Fore.MAGENTA, augmented_iter)
        )
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


if __name__ == "__main__":
    main()
