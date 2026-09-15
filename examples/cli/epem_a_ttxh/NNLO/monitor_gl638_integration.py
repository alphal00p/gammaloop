#!/usr/bin/env python3
"""Compare the latest GL638 iteration snapshots from both sampling setups."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from colorama import Fore, Style, init
from prettytable import PrettyTable, TableStyle


HERE = Path(__file__).resolve().parent
DEFAULT_ADVANCED = HERE / "workspaces" / "GL638_advanced_sampling_workspace"
DEFAULT_OPTIMIZED = HERE / "workspaces" / "GL638_optimized_lmbs_workspace"


def latest_snapshot(workspace: Path) -> tuple[int, dict[str, Any]] | None:
    files = sorted((workspace / "results").glob("integration_result_iter_*.json"))
    if not files:
        return None
    path = files[-1]
    with path.open(encoding="utf-8") as stream:
        document = json.load(stream)
    slots = document.get("slots", [])
    if not slots:
        raise ValueError(f"snapshot {path} contains no integration slots")
    return int(path.stem.rsplit("_", 1)[-1]), slots[0]


def maximum(slot: dict[str, Any], component: str, absolute: bool = False) -> float:
    entries = slot["absolute"]["max_weight_info"] if absolute else slot["max_weight_info"]
    values = [abs(float(item["max_eval"])) for item in entries if item["component"] == component]
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
    r"│\s*All\s+(re|\|re\||im|\|im\|)\s+│\s*"
    r"([+-]?\d+(?:\.\d+)?)\((\d+(?:\.\d+)?)\)e([+-]?\d+)"
    r"\s+.*?[0-9.]+%\s+│\s*[0-9.]+\s+([0-9.]+e[+-]?\d+)"
)
_ABS_OBSERVABLE_RE = re.compile(
    r"│\s*(?:All\s+)?(\|re\||\|im\|)\s+│\s*"
    r"([+-]?\d+(?:\.\d+)?)\((\d+(?:\.\d+)?)\)e([+-]?\d+)"
)
_MAX_RE = re.compile(
    r"│\s*epem_a_tth@NNLO\s+(re|im)\s+\[([+-])\]\s*│\s*"
    r"([+-]?\d+(?:\.\d+)?)e([+-]?\d+)"
)
_TIMING_RE = re.compile(
    r"│\s+timing\s+│\s+total\s+([0-9.]+)\s+ms\s+│\s+param\s+"
    r"([0-9.]+)\s+ms\s+│\s+itg\s+([0-9.]+)\s+ms\s+│\s+evaluators\s+"
    r"([0-9.]+)\s+ms"
)
_EVAL_RE = re.compile(
    r"│\s+evals\s+│\s+f64\s+([0-9.]+)%\s+│\s+f128\s+"
    r"([0-9.]+)%\s+│\s+arb\s+([0-9.]+)%\s+│\s+nans\+unstable\s+"
    r"([0-9.]+)%"
)


def _scientific_value(coefficient: str, exponent: str) -> float:
    return float(f"{coefficient}e{exponent}")


def _compact_value(value: str, uncertainty: str, exponent: str) -> tuple[float, float]:
    coefficient = float(value)
    uncertainty_coefficient = float(uncertainty)
    if "." not in uncertainty:
        uncertainty_coefficient /= 10 ** len(value.split(".", 1)[1])
    scale = 10 ** int(exponent)
    return coefficient * scale, uncertainty_coefficient * scale


def live_metrics(log_path: Path) -> tuple[int, dict[str, str]] | None:
    """Parse the latest completed iteration table emitted to an integration log.

    The non-TTY renderer does not emit partial ratatui updates. Consequently,
    ``--live`` means the latest completed ordinary iteration summary already
    flushed to the log, while a currently running iteration remains visible as
    that last completed summary.
    """
    if not log_path.exists():
        return None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    iterations = list(_ITERATION_RE.finditer(text))
    if not iterations:
        return None
    iteration = iterations[-1]
    latest_start = text.rfind("Streaming integration updates disabled")
    if latest_start > iteration.end():
        return None
    samples_text, suffix = iteration.group(2), iteration.group(3)
    samples = float(samples_text) * {"": 1, "K": 1e3, "M": 1e6, "G": 1e9}[suffix]
    end = iterations[-1].end()
    tail = text[end:]
    values: dict[str, tuple[float, float]] = {}
    impacts: dict[str, float] = {}
    for match in _OBSERVABLE_RE.finditer(tail):
        label = match.group(1)
        values[label] = _compact_value(
            match.group(2), match.group(3), match.group(4)
        )
        impacts[label] = float(match.group(5))
    for match in _ABS_OBSERVABLE_RE.finditer(tail):
        values[match.group(1)] = _compact_value(
            match.group(2), match.group(3), match.group(4)
        )
    timing = _TIMING_RE.search(tail)
    evaluation = _EVAL_RE.search(tail)
    maxima: dict[str, float] = {}
    for match in _MAX_RE.finditer(tail):
        key = match.group(1)
        maxima[key] = max(
            maxima.get(key, 0.0), _scientific_value(match.group(3), match.group(4))
        )
    if not {"re", "|re|", "im", "|im|"} <= values.keys() or timing is None or evaluation is None:
        return None
    total_ms, parameterization_ms, integrand_ms, evaluator_ms = map(float, timing.groups())
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
        "--live",
        action="store_true",
        help="read the latest completed ordinary iteration summaries from integration.log",
    )
    parser.add_argument("--advanced-workspace", type=Path, default=DEFAULT_ADVANCED)
    parser.add_argument("--optimized-workspace", type=Path, default=DEFAULT_OPTIMIZED)
    args = parser.parse_args()
    init(autoreset=True)

    if args.live:
        advanced_live = live_metrics(
            HERE / "gammaloop_state_GL638_advanced_sampling" / "integration.log"
        )
        optimized_live = live_metrics(
            HERE / "gammaloop_state_GL638_optimized_lmbs" / "integration.log"
        )
        advanced_iter = advanced_live[0] if advanced_live else None
        optimized_iter = optimized_live[0] if optimized_live else None
        advanced = advanced_live[1] if advanced_live else None
        optimized = optimized_live[1] if optimized_live else None
    else:
        advanced_snapshot = latest_snapshot(args.advanced_workspace)
        optimized_snapshot = latest_snapshot(args.optimized_workspace)
        advanced_iter = advanced_snapshot[0] if advanced_snapshot else None
        optimized_iter = optimized_snapshot[0] if optimized_snapshot else None
        advanced = metrics(advanced_snapshot[1]) if advanced_snapshot else None
        optimized = metrics(optimized_snapshot[1]) if optimized_snapshot else None

    table = PrettyTable()
    table.set_style(TableStyle.SINGLE_BORDER)
    table.field_names = [
        "Observable / statistic",
        "Advanced sampling",
        "Optimized-LMB only",
    ]
    table.align["Observable / statistic"] = "l"
    table.align["Advanced sampling"] = "r"
    table.align["Optimized-LMB only"] = "r"
    table.title = (
        f"GL638 {'live summary' if args.live else 'latest iteration'}  {Fore.CYAN}advanced #{advanced_iter or 'pending'}{Style.RESET_ALL}  |  "
        f"{Fore.YELLOW}optimized #{optimized_iter or 'pending'}{Style.RESET_ALL}"
    )
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
        left = advanced[key] if advanced else "pending (no completed iteration)"
        right = optimized[key] if optimized else "pending (no completed iteration)"
        if key == "unstable":
            left = (Fore.GREEN if left.startswith("0.000") else Fore.RED) + left + Style.RESET_ALL
            right = (Fore.GREEN if right.startswith("0.000") else Fore.RED) + right + Style.RESET_ALL
        else:
            left = Fore.CYAN + left + Style.RESET_ALL
            right = Fore.YELLOW + right + Style.RESET_ALL
        table.add_row([label, left, right])
    print(table)


if __name__ == "__main__":
    main()
