#!/usr/bin/env python3
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys
import time
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore[no-redef]


ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
DEFAULT_OUTPUT = "outputs/bench_result.json"
DEFAULT_CONFIG = "bench_config.toml"
DEFAULT_GRAPHS = ["G1"]
DEFAULT_PROFILES = ["minimal"]
DEFAULT_MAX_GEN_RAM = "15GB"
DEFAULT_MAX_GEN_TIME = "120s"
DEFAULT_BINARY_PROFILE = "release"
DEFAULT_TARGET_PROFILE_TIME = "10s"
USE_COLOUR = True
RESET = "\033[0m"
COLORS = {
    "bold": "\033[1m",
    "dim": "\033[2m",
    "red": "\033[31m",
    "green": "\033[32m",
    "yellow": "\033[33m",
    "blue": "\033[34m",
    "magenta": "\033[35m",
    "cyan": "\033[36m",
}


def color(text: Any, *styles: str) -> str:
    value = str(text)
    if not USE_COLOUR:
        return value
    prefix = "".join(COLORS[style] for style in styles if style in COLORS)
    return f"{prefix}{value}{RESET}" if prefix else value


def visible_len(text: str) -> int:
    return len(ANSI_RE.sub("", text))


def pad(text: str, width: int, align: str = "left") -> str:
    extra = max(width - visible_len(text), 0)
    if align == "right":
        return " " * extra + text
    if align == "center":
        left = extra // 2
        return " " * left + text + " " * (extra - left)
    return text + " " * extra


def parse_duration_seconds(value: str) -> float:
    raw = str(value).strip().lower()
    match = re.fullmatch(r"([0-9]*\.?[0-9]+)\s*([a-z]*)", raw)
    if not match:
        raise argparse.ArgumentTypeError(f"invalid duration: {value}")
    number = float(match.group(1))
    unit = match.group(2) or "s"
    factors = {
        "s": 1.0,
        "sec": 1.0,
        "secs": 1.0,
        "second": 1.0,
        "seconds": 1.0,
        "m": 60.0,
        "min": 60.0,
        "mins": 60.0,
        "minute": 60.0,
        "minutes": 60.0,
        "h": 3600.0,
        "hr": 3600.0,
        "hour": 3600.0,
        "hours": 3600.0,
    }
    if unit not in factors:
        raise argparse.ArgumentTypeError(f"invalid duration unit in {value}")
    return number * factors[unit]


def parse_size_bytes(value: str) -> int:
    raw = str(value).strip().lower().replace(" ", "")
    match = re.fullmatch(r"([0-9]*\.?[0-9]+)([kmgt]?i?b?|)", raw)
    if not match:
        raise argparse.ArgumentTypeError(f"invalid size: {value}")
    number = float(match.group(1))
    unit = match.group(2) or "gb"
    factors = {
        "b": 1,
        "kb": 1000,
        "k": 1000,
        "mb": 1000**2,
        "m": 1000**2,
        "gb": 1000**3,
        "g": 1000**3,
        "tb": 1000**4,
        "t": 1000**4,
        "kib": 1024,
        "mib": 1024**2,
        "gib": 1024**3,
        "tib": 1024**4,
    }
    if unit not in factors:
        raise argparse.ArgumentTypeError(f"invalid size unit in {value}")
    return int(number * factors[unit])


def format_sig(value: float, sig: int = 3) -> str:
    if value is None or not math.isfinite(value):
        return "N/A"
    if value == 0:
        return "0"
    magnitude = math.floor(math.log10(abs(value)))
    if magnitude >= sig or magnitude < -2:
        return f"{value:.{sig - 1}e}"
    decimals = max(sig - magnitude - 1, 0)
    return f"{value:.{decimals}f}"


def format_with_uncertainty(value: float, error: float | None) -> str:
    if error is None or error <= 0:
        return format_sig(value)
    deviation = abs(error)

    def decimals_for(positive_value: float, offset: int) -> int:
        decimals = int(offset - math.log10(positive_value))
        if decimals > 0 and positive_value * 10.0**decimals >= [0.5, 9.5, 99.5][offset]:
            decimals -= 1
        return max(decimals, 0)

    if math.isnan(value) or math.isnan(deviation):
        return f"{value:e} ± {deviation:e}"
    if math.isinf(deviation):
        return f"{value:e} ± inf"
    if not math.isfinite(value):
        return f"{value:e} ± {deviation:e}"
    if value == 0.0 and not (1e-4 <= deviation < 1e5):
        if deviation == 0.0:
            return "0(0)"
        mantissa, exponent = f"{deviation:.1e}".split("e")
        return f"0.0({mantissa})e{exponent}"
    if value == 0.0:
        if deviation >= 9.95:
            return f"0({deviation:.0f})"
        if deviation >= 0.995:
            return f"0.0({deviation:.1f})"
        decimals = decimals_for(deviation, 2)
        return f"{value:.{decimals}f}({deviation * 10.0**decimals:.0f})"
    if deviation == 0.0:
        mantissa, exponent = f"{value:e}".split("e")
        return f"{mantissa}(0)e{exponent}" if exponent != "0" else f"{mantissa}(0)"
    if deviation > 1e4 * abs(value):
        return f"{value:.1e} ± {deviation:.2e}"
    if abs(value) >= 1e6 or abs(value) < 1e-5:
        exponent = math.floor(math.log10(abs(value)))
        scale = 10.0**exponent
        mantissa = format_with_uncertainty(value / scale, deviation / scale)
        _, exponent_text = f"{scale:.0e}".split("e")
        return f"{mantissa}e{exponent_text}"
    if deviation >= 9.95:
        if abs(value) >= 9.5:
            return f"{value:.0f}({deviation:.0f})"
        decimals = decimals_for(abs(value), 1)
        return f"{value:.{decimals}f}({deviation:.{decimals}f})"
    if deviation >= 0.995:
        if abs(value) >= 0.95:
            return f"{value:.1f}({deviation:.1f})"
        decimals = decimals_for(abs(value), 1)
        return f"{value:.{decimals}f}({deviation:.{decimals}f})"

    decimals = max(decimals_for(abs(value), 1), decimals_for(deviation, 2))
    return f"{value:.{decimals}f}({deviation * 10.0**decimals:.0f})"


def format_seconds(seconds: float | None, stderr: float | None = None) -> str:
    if seconds is None or not math.isfinite(seconds):
        return color("N/A", "red")
    abs_value = abs(seconds)
    units = [
        (1e-9, "ns"),
        (1e-6, "µs"),
        (1e-3, "ms"),
        (1.0, "s"),
        (60.0, "min"),
        (3600.0, "h"),
    ]
    scale, unit = 1.0, "s"
    for candidate_scale, candidate_unit in units:
        scaled = abs_value / candidate_scale if candidate_scale else abs_value
        if 0.1 <= scaled < 1000:
            scale, unit = candidate_scale, candidate_unit
            break
    else:
        if abs_value >= 3600:
            scale, unit = 3600.0, "h"
        elif abs_value < 1e-9 and abs_value > 0:
            scale, unit = 1e-9, "ns"
    scaled_value = seconds / scale
    scaled_error = stderr / scale if stderr is not None else None
    return f"{format_with_uncertainty(scaled_value, scaled_error)} {unit}"


def format_bytes(byte_count: float | None) -> str:
    if byte_count is None or not math.isfinite(byte_count):
        return color("N/A", "red")
    abs_value = abs(byte_count)
    units = [(1.0, "B"), (1000.0, "KB"), (1000.0**2, "MB"), (1000.0**3, "GB")]
    scale, unit = units[0]
    for candidate_scale, candidate_unit in units:
        scaled = abs_value / candidate_scale
        if scaled < 1000:
            scale, unit = candidate_scale, candidate_unit
            break
    return f"{format_sig(byte_count / scale)} {unit}"


def format_limit_seconds(seconds: float | None) -> str:
    return format_seconds(seconds) if seconds is not None else color("no cap", "dim")


def format_multiplier(value: float | None, reference: float | None) -> str:
    if value is None or reference is None or reference <= 0 or not math.isfinite(value):
        return color("N/A", "red")
    multiplier = value / reference
    rendered = f"x{format_sig(multiplier)}"
    return color(rendered, "red" if multiplier > 2.0 else "green")


def format_with_optional_multiplier(
    rendered: str, value: float | None, reference: float | None
) -> str:
    if (
        value is None
        or not math.isfinite(value)
        or reference is None
        or reference <= 0
        or not math.isfinite(reference)
    ):
        return rendered
    return f"{rendered} ({format_multiplier(value, reference)})"


def optional_alpha(graph: dict[str, Any]) -> dict[str, Any] | None:
    alpha = graph.get("alphaloop")
    return alpha if isinstance(alpha, dict) else None


def alpha_float(alpha: dict[str, Any] | None, key: str) -> float | None:
    if alpha is None or key not in alpha:
        return None
    try:
        value = float(alpha[key])
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def alpha_peak_rss_bytes(alpha: dict[str, Any] | None) -> float | None:
    if alpha is None:
        return None
    for key, scale in (
        ("peak_rss_bytes", 1.0),
        ("peak_rss_mb", 1_000_000.0),
        ("rss_mb", 1_000_000.0),
        ("peak_ram_mb", 1_000_000.0),
        ("peak_rss_gb", 1_000_000_000.0),
        ("rss_gb", 1_000_000_000.0),
    ):
        value = alpha_float(alpha, key)
        if value is not None:
            return value * scale
    return None


def alpha_disk_bytes(alpha: dict[str, Any] | None) -> float | None:
    disk_mb = alpha_float(alpha, "disk_mb")
    return disk_mb * 1_000_000.0 if disk_mb is not None else None


def alpha_eval_seconds(alpha: dict[str, Any] | None) -> float | None:
    t_eval_ms = alpha_float(alpha, "t_eval_ms")
    return t_eval_ms / 1000.0 if t_eval_ms is not None else None


def format_optional_reference(
    value: float | None, formatter: Any, missing: str = "N/A"
) -> str:
    return formatter(value) if value is not None else color(missing, "dim")


def is_status_table_cell(cell: str) -> bool:
    plain = ANSI_RE.sub("", cell).strip()
    return plain in {"N/A", "Error", "Ctrl+C", "Aborted"} or plain.startswith("> ")


def result_cell_alignment(category: str, subcolumn_index: int, cell: str) -> str:
    if is_status_table_cell(cell):
        return "center"
    return "left" if category != "Δ target" and subcolumn_index == 0 else "right"


def metadata_float(
    metadata: dict[str, Any], key: str, default: float | int | str | None = None
) -> float | None:
    value = metadata.get(key, default)
    if value is None:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def command_status_cell(
    run: dict[str, Any], stage: str, metadata: dict[str, Any]
) -> str:
    stage_result = run.get(stage) or {}
    status = stage_result.get("status") or run.get("status")
    abort_reason = str(
        stage_result.get("abort_reason") or run.get("abort_reason") or ""
    )
    returncode = stage_result.get("returncode")

    if status == "interrupted" or "Ctrl+C" in abort_reason:
        return color("Ctrl+C", "red")
    if stage == "generation" and status == "aborted":
        if "RAM" in abort_reason:
            cap = metadata_float(
                metadata, "max_gen_ram_bytes", parse_size_bytes(DEFAULT_MAX_GEN_RAM)
            )
            return color(f"> {format_bytes(cap)}", "red")
        if "exceeded" in abort_reason:
            cap = metadata_float(
                metadata, "max_gen_time_s", parse_duration_seconds(DEFAULT_MAX_GEN_TIME)
            )
            return color(f"> {format_seconds(cap)}", "red")
    if status == "error" or (returncode is not None and returncode != 0):
        return color("Error", "red")
    if status == "aborted":
        return color("Aborted", "red")
    return color("N/A", "red")


def generation_unavailable_cell(
    run: dict[str, Any], metadata: dict[str, Any]
) -> str | None:
    generation = run.get("generation") or {}
    if generation.get("status") == "ok" and generation.get("returncode") == 0:
        return None
    if generation:
        return command_status_cell(run, "generation", metadata)
    if run.get("status") == "error":
        return color("Error", "red")
    return color("N/A", "red")


def benchmark_unavailable_cell(
    run: dict[str, Any], metadata: dict[str, Any]
) -> str | None:
    generation_cell = generation_unavailable_cell(run, metadata)
    if generation_cell is not None:
        return generation_cell
    benchmark = run.get("benchmark") or {}
    if run.get("status") == "ok" and benchmark.get("status") == "ok":
        return None
    if benchmark:
        return command_status_cell(run, "benchmark", metadata)
    if run.get("status") == "error":
        return color("Error", "red")
    return color("N/A", "red")


def table(
    title: str,
    headers: list[str],
    rows: list[list[str]],
    align: list[str] | None = None,
) -> str:
    align = align or ["left"] * len(headers)
    widths = [
        max([visible_len(headers[i])] + [visible_len(row[i]) for row in rows])
        for i in range(len(headers))
    ]
    top = "┌" + "┬".join("─" * (width + 2) for width in widths) + "┐"
    sep = "├" + "┼".join("─" * (width + 2) for width in widths) + "┤"
    bottom = "└" + "┴".join("─" * (width + 2) for width in widths) + "┘"
    lines = [color(title, "bold", "cyan"), top]
    lines.append(
        "│"
        + "│".join(
            f" {pad(color(headers[i], 'bold', 'blue'), widths[i], align[i])} "
            for i in range(len(headers))
        )
        + "│"
    )
    lines.append(sep)
    for row in rows:
        lines.append(
            "│"
            + "│".join(
                f" {pad(row[i], widths[i], align[i])} " for i in range(len(headers))
            )
            + "│"
        )
    lines.append(bottom)
    return "\n".join(lines)


def grouped_result_table(
    graphs: list[str],
    profiles: list[str],
    results: dict[str, dict[str, Any]],
    graph_configs: dict[str, Any],
    use_evaluator_timings: bool,
    compare_eval: bool,
    binary_profile: str,
    metadata: dict[str, Any],
) -> str:
    categories = ["t_gen", "peak_rss", "disk_size", "t_eval"]
    if compare_eval:
        categories.append("Δ target")
    category_subheaders = {
        category: profiles if category == "Δ target" else ["AlphaLoop", *profiles]
        for category in categories
    }
    category_labels = {
        category: (
            "t_eval " + ("evaluator" if use_evaluator_timings else "total")
            if category == "t_eval"
            else "peak RSS"
            if category == "peak_rss"
            else category
        )
        for category in categories
    }
    subcolumn_gap = "  "

    row_values: dict[str, dict[str, list[str]]] = {}
    for graph_name in graphs:
        graph_result = results.get(graph_name, {})
        alpha = optional_alpha(graph_configs[graph_name])
        row_values[graph_name] = {}
        for category in categories:
            cells: list[str] = []
            if category == "t_gen":
                reference = alpha_float(alpha, "t_gen_s")
                cells.append(format_optional_reference(reference, format_seconds))
                for profile in profiles:
                    run = graph_result.get(profile, {})
                    unavailable = generation_unavailable_cell(run, metadata)
                    if unavailable is not None:
                        cells.append(unavailable)
                        continue
                    seconds = run.get("generation", {}).get("seconds")
                    if seconds is None:
                        cells.append(color("N/A", "red"))
                    else:
                        cells.append(
                            format_with_optional_multiplier(
                                format_seconds(seconds), seconds, reference
                            )
                        )
            elif category == "peak_rss":
                reference = alpha_peak_rss_bytes(alpha)
                cells.append(format_optional_reference(reference, format_bytes))
                for profile in profiles:
                    run = graph_result.get(profile, {})
                    unavailable = generation_unavailable_cell(run, metadata)
                    if unavailable is not None:
                        cells.append(unavailable)
                        continue
                    peak_rss = run.get("generation", {}).get("peak_rss_bytes")
                    if peak_rss is None:
                        cells.append(color("N/A", "red"))
                    else:
                        cells.append(
                            format_with_optional_multiplier(
                                format_bytes(float(peak_rss)),
                                float(peak_rss),
                                reference,
                            )
                        )
            elif category == "disk_size":
                disk_mb = alpha_float(alpha, "disk_mb")
                reference = disk_mb * 1_000_000.0 if disk_mb is not None else None
                cells.append(format_optional_reference(reference, format_bytes))
                for profile in profiles:
                    run = graph_result.get(profile, {})
                    unavailable = generation_unavailable_cell(run, metadata)
                    if unavailable is not None:
                        cells.append(unavailable)
                        continue
                    size = run.get("disk_size_bytes")
                    if size is None:
                        cells.append(color("N/A", "red"))
                    else:
                        cells.append(
                            format_with_optional_multiplier(
                                format_bytes(size), size, reference
                            )
                        )
            elif category == "t_eval":
                t_eval_ms = alpha_float(alpha, "t_eval_ms")
                reference = t_eval_ms / 1000.0 if t_eval_ms is not None else None
                cells.append(format_optional_reference(reference, format_seconds))
                timing_key = "evaluator" if use_evaluator_timings else "Total"
                for profile in profiles:
                    run = graph_result.get(profile, {})
                    unavailable = benchmark_unavailable_cell(run, metadata)
                    if unavailable is not None:
                        cells.append(unavailable)
                        continue
                    benchmark = run.get("benchmark") or {}
                    timings = benchmark.get("timings", {})
                    timing = timings.get(timing_key)
                    if timing is None and timing_key == "evaluator":
                        timing = timings.get("evaluators")
                    if timing is None and timing_key == "Total":
                        timing = timings.get("total")
                    if not timing:
                        cells.append(color("N/A", "red"))
                    else:
                        mean = timing.get("mean_s")
                        stderr = timing.get("stderr_s")
                        cells.append(
                            format_with_optional_multiplier(
                                format_seconds(mean, stderr), mean, reference
                            )
                        )
            else:
                for profile in profiles:
                    run = graph_result.get(profile, {})
                    unavailable = benchmark_unavailable_cell(run, metadata)
                    if unavailable is not None:
                        cells.append(unavailable)
                        continue
                    benchmark = run.get("benchmark") or {}
                    delta = benchmark.get("delta_target")
                    if delta is None:
                        cells.append(color("N/A", "red"))
                    else:
                        cells.append(
                            color(
                                format_sig(delta), "green" if delta < 1.0e-5 else "red"
                            )
                        )
            row_values[graph_name][category] = cells

    graph_width = max(
        [len("graph")] + [visible_len(color(name, "bold")) for name in graphs]
    )
    subwidths: dict[str, list[int]] = {}
    group_widths: dict[str, int] = {}
    for category in categories:
        subheaders = category_subheaders[category]
        widths = []
        for index, subheader in enumerate(subheaders):
            widths.append(
                max(
                    visible_len(subheader),
                    *[
                        visible_len(row_values[graph_name][category][index])
                        for graph_name in graphs
                    ],
                )
            )
        subwidths[category] = widths
        subcolumns_width = sum(widths) + len(subcolumn_gap) * max(len(widths) - 1, 0)
        group_widths[category] = max(
            subcolumns_width,
            visible_len(category_labels[category]),
        )

    group_sizes = [graph_width, *[group_widths[category] for category in categories]]
    top = "┌" + "┬".join("─" * (width + 2) for width in group_sizes) + "┐"
    sep = "├" + "┼".join("─" * (width + 2) for width in group_sizes) + "┤"
    bottom = "└" + "┴".join("─" * (width + 2) for width in group_sizes) + "┘"

    header_1 = [" graph ".ljust(graph_width + 2)]
    for category in categories:
        header_1.append(
            " "
            + pad(
                color(category_labels[category], "bold", "blue"),
                group_widths[category],
                "center",
            )
            + " "
        )
    header_2 = [" " + pad("", graph_width) + " "]
    for category in categories:
        subheaders = category_subheaders[category]
        pieces = [
            pad(
                color(subheaders[i], "bold", "magenta"),
                subwidths[category][i],
                "center",
            )
            for i in range(len(subheaders))
        ]
        header_2.append(
            " "
            + pad(subcolumn_gap.join(pieces), group_widths[category], "center")
            + " "
        )

    lines = [
        f"{color('Results', 'bold', 'cyan')} ({binary_profile} build)",
        top,
        "│" + "│".join(header_1) + "│",
        "│" + "│".join(header_2) + "│",
        sep,
    ]
    for graph_name in graphs:
        groups = [" " + pad(color(graph_name, "bold"), graph_width) + " "]
        for category in categories:
            subheaders = category_subheaders[category]
            pieces = [
                pad(
                    row_values[graph_name][category][i],
                    subwidths[category][i],
                    result_cell_alignment(
                        category, i, row_values[graph_name][category][i]
                    ),
                )
                for i in range(len(subheaders))
            ]
            groups.append(
                " "
                + pad(
                    subcolumn_gap.join(pieces),
                    group_widths[category],
                    "right" if category == "Δ target" else "left",
                )
                + " "
            )
        lines.append("│" + "│".join(groups) + "│")
    lines.append(bottom)
    return "\n".join(lines)


def load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def placeholder_defaults_from_card(path: Path) -> dict[str, str]:
    try:
        content = path.read_text(encoding="utf-8")
    except OSError:
        return {}
    defaults: dict[str, str] = {}
    for match in re.finditer(r"\$\(([A-Za-z_][A-Za-z0-9_]*):([^)]*)\)", content):
        defaults.setdefault(match.group(1), match.group(2))
    return defaults


def graph_point(graph: dict[str, Any]) -> str:
    return " ".join(format(value, ".17g") for value in graph["x"])


def process_tree_rss_bytes(pid: int) -> int:
    try:
        output = subprocess.check_output(["ps", "-axo", "pid=,ppid=,rss="], text=True)
    except (OSError, subprocess.SubprocessError):
        return 0
    children: dict[int, list[int]] = {}
    rss_by_pid: dict[int, int] = {}
    for line in output.splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            child_pid, parent_pid, rss_kb = map(int, parts)
        except ValueError:
            continue
        children.setdefault(parent_pid, []).append(child_pid)
        rss_by_pid[child_pid] = rss_kb * 1024
    stack = [pid]
    total = 0
    seen: set[int] = set()
    while stack:
        current = stack.pop()
        if current in seen:
            continue
        seen.add(current)
        total += rss_by_pid.get(current, 0)
        stack.extend(children.get(current, []))
    return total


class InterruptState:
    def __init__(self, repeat_window_seconds: float = 2.0) -> None:
        self.repeat_window_seconds = repeat_window_seconds
        self.last_interrupt = 0.0
        self.abort_all = False

    def register(self) -> bool:
        now = time.monotonic()
        repeated = now - self.last_interrupt <= self.repeat_window_seconds
        self.last_interrupt = now
        if repeated:
            self.abort_all = True
        return self.abort_all


def clear_live_line() -> None:
    if USE_COLOUR:
        print("\r\033[K", end="", file=sys.stderr, flush=True)
    else:
        print("\r" + " " * 180 + "\r", end="", file=sys.stderr, flush=True)


def print_process_header(
    graph_name: str,
    profile_name: str,
    stage: str,
    max_seconds: float | None,
    max_rss_bytes: int | None,
    command_file: Path,
    output_file: Path,
    state_dir: Path,
    cwd: Path,
) -> None:
    cap_text = (
        f"time ≤ {format_limit_seconds(max_seconds)}, RAM ≤ {format_bytes(float(max_rss_bytes))}"
        if max_seconds is not None and max_rss_bytes is not None
        else "generation caps disabled for this run"
    )
    rows = [
        f"{color(stage.upper(), 'bold', 'magenta')} {color(graph_name, 'bold', 'green')} / {color(profile_name, 'bold', 'blue')}",
        f"sandbox: {cap_text}",
        f"command: {color(display_path(command_file, cwd), 'dim')}",
        f"output: {color(display_path(output_file, cwd), 'dim')}",
        f"state: {color(display_path(state_dir, cwd), 'dim')}",
    ]
    width = max(78, *(visible_len(row) for row in rows))
    lines = [
        "",
        color("╭" + "─" * (width + 2) + "╮", "cyan"),
        *[
            color("│", "cyan") + " " + pad(row, width) + " " + color("│", "cyan")
            for row in rows
        ],
        color("╰" + "─" * (width + 2) + "╯", "cyan"),
    ]
    print("\n".join(lines), file=sys.stderr)


def live_process_line(
    graph_name: str,
    profile_name: str,
    stage: str,
    elapsed: float,
    rss_bytes: int,
    max_seconds: float | None,
    max_rss_bytes: int | None,
) -> str:
    elapsed_color = "green"
    memory_color = "green"
    if max_seconds is not None and elapsed > 0.8 * max_seconds:
        elapsed_color = "yellow"
    if max_seconds is not None and elapsed > max_seconds:
        elapsed_color = "red"
    if max_rss_bytes is not None and rss_bytes > 0.8 * max_rss_bytes:
        memory_color = "yellow"
    if max_rss_bytes is not None and rss_bytes > max_rss_bytes:
        memory_color = "red"
    elapsed_part = f"{color(format_seconds(elapsed), elapsed_color)} / {format_limit_seconds(max_seconds)}"
    memory_part = (
        f"{color(format_bytes(float(rss_bytes)), memory_color)} / "
        f"{format_bytes(float(max_rss_bytes)) if max_rss_bytes is not None else color('no cap', 'dim')}"
    )
    line_prefix = "\r\033[K" if USE_COLOUR else "\r"
    return (
        f"{line_prefix}{color('⟳', 'cyan')} "
        f"{color(graph_name, 'bold')} {color(profile_name, 'blue')} {stage} | "
        f"elapsed {elapsed_part} | RSS {memory_part}"
    )


def terminate_process_group(proc: subprocess.Popen[str]) -> None:
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return


def display_path(path: Path, cwd: Path) -> str:
    try:
        relative = path.resolve().relative_to(cwd.resolve())
    except ValueError:
        return str(path)
    return f"./{relative}" if str(relative) != "." else "."


def write_command_script(
    command_file: Path, output_file: Path, cwd: Path, cmd: list[str]
) -> None:
    command_file.parent.mkdir(parents=True, exist_ok=True)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    command_line = " ".join(shlex.quote(part) for part in cmd)
    output_target = display_path(output_file, cwd)
    command_file.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f"cd {shlex.quote(str(cwd))}\n\n"
        f"{command_line} > {shlex.quote(output_target)} 2>&1\n",
        encoding="utf-8",
    )
    command_file.chmod(0o755)


def run_monitored(
    cmd: list[str],
    cwd: Path,
    graph_name: str,
    profile_name: str,
    stage: str,
    interrupt_state: InterruptState,
    command_file: Path,
    output_file: Path,
    state_dir: Path,
    max_seconds: float | None = None,
    max_rss_bytes: int | None = None,
) -> dict[str, Any]:
    write_command_script(command_file, output_file, cwd, cmd)
    print_process_header(
        graph_name,
        profile_name,
        stage,
        max_seconds,
        max_rss_bytes,
        command_file,
        output_file,
        state_dir,
        cwd,
    )
    start = time.monotonic()
    peak_rss = 0
    abort_reason = None
    interrupted = False
    abort_all = False
    with output_file.open("w", encoding="utf-8") as output_handle:
        proc = subprocess.Popen(
            cmd,
            cwd=cwd,
            stdout=output_handle,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=True,
        )
        try:
            while proc.poll() is None:
                elapsed = time.monotonic() - start
                peak_rss = max(peak_rss, process_tree_rss_bytes(proc.pid))
                print(
                    live_process_line(
                        graph_name,
                        profile_name,
                        stage,
                        elapsed,
                        peak_rss,
                        max_seconds,
                        max_rss_bytes,
                    ),
                    end="",
                    file=sys.stderr,
                    flush=True,
                )
                if max_seconds is not None and elapsed > max_seconds:
                    abort_reason = f"{stage} exceeded {max_seconds:.3g}s"
                    break
                if max_rss_bytes is not None and peak_rss > max_rss_bytes:
                    abort_reason = (
                        f"{stage} exceeded {format_bytes(float(max_rss_bytes))} RAM"
                    )
                    break
                time.sleep(0.25)
        except KeyboardInterrupt:
            interrupted = True
            abort_all = interrupt_state.register()
            abort_reason = (
                "received repeated Ctrl+C; aborting remaining runs"
                if abort_all
                else "received Ctrl+C; aborting current graph/profile"
            )
        if abort_reason:
            terminate_process_group(proc)
            try:
                proc.wait(timeout=1)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                proc.wait()
        else:
            proc.wait()
    elapsed = time.monotonic() - start
    peak_rss = max(peak_rss, process_tree_rss_bytes(proc.pid))
    clear_live_line()
    status = "aborted" if abort_reason else ("ok" if proc.returncode == 0 else "failed")
    status_color = "green" if status == "ok" else "red"
    print(
        f"{color('✓' if status == 'ok' else '✗', status_color)} "
        f"{graph_name}/{profile_name} {stage}: {color(status, status_color)} "
        f"in {format_seconds(elapsed)}, peak RSS {format_bytes(float(peak_rss))}"
        + (f" ({abort_reason})" if abort_reason else ""),
        file=sys.stderr,
    )
    try:
        output = output_file.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        output = f"<could not read output log {output_file}: {exc}>"
    return {
        "returncode": proc.returncode,
        "seconds": elapsed,
        "peak_rss_bytes": peak_rss,
        "output": output,
        "stdout": output,
        "stderr": "",
        "command_file": str(command_file),
        "output_log": str(output_file),
        "aborted": abort_reason is not None,
        "abort_reason": abort_reason,
        "interrupted": interrupted,
        "abort_all": abort_all,
    }


def shell_join_define(key: str, value: Any) -> str:
    return f"-D {shlex.quote(f'{key}={value}')}"


def run_block_command(block: str, defines: dict[str, Any]) -> str:
    parts = ["run", block]
    for key, value in defines.items():
        if value is None:
            continue
        parts.append(shell_join_define(key, value))
    return " ".join(parts)


def gammaloop_command_sequence(block: str, defines: dict[str, Any]) -> str:
    return f"{run_block_command(block, defines)}; quit -o"


def state_size_bytes(state_dir: Path) -> int:
    if not state_dir.exists():
        return 0
    total = 0
    for path in state_dir.rglob("*"):
        if path.is_file():
            try:
                total += path.stat().st_size
            except OSError:
                pass
    return total


def sanitize_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def print_error_log_notice(path: Path, cwd: Path) -> None:
    relative_path = display_path(path, cwd)
    print(
        f"{color('↳ error log', 'bold', 'red')}: {color(relative_path, 'bold', 'yellow')}",
        file=sys.stderr,
    )


def write_error_log(
    path: Path,
    title: str,
    cmd: list[str],
    result: dict[str, Any],
    cwd: Path | None = None,
) -> None:
    body = [
        title,
        "",
        "Command file:",
        result.get("command_file", ""),
        "",
        "Output log:",
        result.get("output_log", ""),
        "",
        "Expanded command:",
        " ".join(shlex.quote(part) for part in cmd),
        "",
        f"Return code: {result.get('returncode')}",
        f"Elapsed: {result.get('seconds')}",
    ]
    if result.get("abort_reason"):
        body.append(f"Abort reason: {result['abort_reason']}")
    body.extend(
        [
            "",
            "OUTPUT:",
            result.get("output", result.get("stdout", "")),
        ]
    )
    path.write_text("\n".join(str(part) for part in body), encoding="utf-8")
    if cwd is not None:
        print_error_log_notice(path, cwd)


def command_prefix(repo_root: Path, binary_profile: str) -> list[str]:
    return [str(repo_root / "gammaloop"), f"--{binary_profile}"]


def benchmark_json_timings(
    payload: dict[str, Any],
) -> dict[str, dict[str, float | None]]:
    timings: dict[str, dict[str, float | None]] = {}
    for row in payload.get("summary", []):
        category = row.get("category")
        if not category:
            continue
        raw_category = str(category)
        normalized_category = {
            "evaluators": "evaluator",
            "evaluator": "evaluator",
            "event processing": "event_processing",
            "other/overhead": "other",
            "total": "Total",
        }.get(raw_category.lower(), raw_category)
        timing = {
            "mean_s": row.get("mean_seconds_per_sample"),
            "stderr_s": row.get("standard_error_seconds_per_sample"),
            "percentage_of_total": row.get("percentage_of_total"),
        }
        timings[raw_category] = timing
        timings[normalized_category] = timing
    return timings


def complex_from_any(value: Any) -> complex | None:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return complex(float(value), 0.0)
    if isinstance(value, list) and len(value) >= 2:
        return complex(float(value[0]), float(value[1]))
    if isinstance(value, dict):
        return complex(float(value.get("re", 0.0)), float(value.get("im", 0.0)))
    return None


def relative_delta(value: Any, target: Any) -> float | None:
    observed = complex_from_any(value)
    expected = complex_from_any(target)
    if observed is None or expected is None:
        return None
    denominator = max(abs(expected), 1.0e-300)
    return abs(observed - expected) / denominator


def graph_defines(graph_name: str, graph: dict[str, Any]) -> dict[str, str]:
    defines = {
        "model": str(graph["model"]),
        "graph_path": str(graph["path"]),
        "process": str(graph.get("process", "bench")),
        "integrand": str(graph.get("integrand", "default")),
        "process_spec": str(graph["process_spec"]),
        "point": graph_point(graph),
    }
    for index, value in enumerate(graph["x"]):
        defines[f"x{index}"] = format(value, ".17g")
    return defines


def run_case(
    graph_name: str,
    graph: dict[str, Any],
    profile_name: str,
    profile: dict[str, Any],
    args: argparse.Namespace,
    bench_dir: Path,
    repo_root: Path,
    interrupt_state: InterruptState,
) -> dict[str, Any]:
    safe_case = f"{sanitize_name(graph_name)}_{sanitize_name(profile_name)}"
    outputs_dir = bench_dir / "outputs"
    state_dir = outputs_dir / safe_case
    error_log = outputs_dir / f"{safe_case}_errors.log"
    json_output = outputs_dir / f"{safe_case}_inspect.json"
    generate_command_file = outputs_dir / f"{safe_case}_generate_command.sh"
    generate_output_log = outputs_dir / f"{safe_case}_generate_output.log"
    inspect_command_file = outputs_dir / f"{safe_case}_inspect_command.sh"
    inspect_output_log = outputs_dir / f"{safe_case}_inspect_output.log"
    state_arg = display_path(state_dir, bench_dir)

    base_defines = graph_defines(graph_name, graph)
    profile_defines = {
        key: str(value) for key, value in profile.get("defines", {}).items()
    }
    generate_defines = {**base_defines, **profile_defines}
    bench_defines = {
        **generate_defines,
        "x_dim": str(len(graph["x"])),
        "bench_duration": format(args.target_profile_time, ".17g"),
        "json_output": str(json_output),
    }

    generate_command = gammaloop_command_sequence("generate", generate_defines)
    generate_cmd = [
        *command_prefix(repo_root, args.profile),
        "-s",
        state_arg,
        "--clean-state",
        "bench.toml",
        "run",
        "-c",
        generate_command,
    ]
    generation = run_monitored(
        generate_cmd,
        bench_dir,
        graph_name,
        profile_name,
        "generation",
        interrupt_state,
        generate_command_file,
        generate_output_log,
        state_dir,
        max_seconds=args.max_gen_time,
        max_rss_bytes=args.max_gen_ram,
    )
    result: dict[str, Any] = {
        "status": "ok",
        "state_dir": str(state_dir),
        "error_log": None,
        "generation": {
            "status": "ok",
            "seconds": generation["seconds"],
            "peak_rss_bytes": generation["peak_rss_bytes"],
            "returncode": generation["returncode"],
            "command_file": str(generate_command_file),
            "output_log": str(generate_output_log),
        },
        "disk_size_bytes": None,
        "benchmark": None,
    }
    if generation["aborted"] or generation["returncode"] != 0:
        status = (
            "interrupted"
            if generation.get("interrupted")
            else ("aborted" if generation["aborted"] else "error")
        )
        result["status"] = status
        result["generation"]["status"] = status
        result["generation"]["abort_reason"] = generation.get("abort_reason")
        result["abort_all_requested"] = generation.get("abort_all", False)
        result["error_log"] = str(error_log)
        write_error_log(
            error_log,
            f"{graph_name}/{profile_name} generation {status}",
            generate_cmd,
            generation,
            bench_dir,
        )
        result["disk_size_bytes"] = state_size_bytes(state_dir)
        return result

    result["disk_size_bytes"] = state_size_bytes(state_dir)

    bench_command = gammaloop_command_sequence("bench", bench_defines)
    bench_cmd = [
        *command_prefix(repo_root, args.profile),
        "-s",
        state_arg,
        "bench.toml",
        "run",
        "-c",
        bench_command,
    ]
    bench = run_monitored(
        bench_cmd,
        bench_dir,
        graph_name,
        profile_name,
        "inspect",
        interrupt_state,
        inspect_command_file,
        inspect_output_log,
        state_dir,
    )
    if bench.get("aborted"):
        result["status"] = "interrupted" if bench.get("interrupted") else "aborted"
        result["error_log"] = str(error_log)
        result["abort_all_requested"] = bench.get("abort_all", False)
        write_error_log(
            error_log,
            f"{graph_name}/{profile_name} benchmark {result['status']}",
            bench_cmd,
            bench,
            bench_dir,
        )
        result["benchmark"] = {
            "status": result["status"],
            "seconds": bench["seconds"],
            "returncode": bench["returncode"],
            "json_output": str(json_output),
            "command_file": str(inspect_command_file),
            "output_log": str(inspect_output_log),
            "abort_reason": bench.get("abort_reason"),
        }
        return result
    if bench["returncode"] != 0:
        result["status"] = "error"
        result["error_log"] = str(error_log)
        write_error_log(
            error_log,
            f"{graph_name}/{profile_name} benchmark error",
            bench_cmd,
            bench,
            bench_dir,
        )
        result["benchmark"] = {
            "status": "error",
            "seconds": bench["seconds"],
            "returncode": bench["returncode"],
            "json_output": str(json_output),
            "command_file": str(inspect_command_file),
            "output_log": str(inspect_output_log),
        }
        return result

    try:
        payload = json.loads(json_output.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        bench["stderr"] += (
            f"\nCould not read inspect JSON output {json_output}: {exc}\n"
        )
        result["status"] = "error"
        result["error_log"] = str(error_log)
        write_error_log(
            error_log,
            f"{graph_name}/{profile_name} benchmark JSON error",
            bench_cmd,
            bench,
            bench_dir,
        )
        result["benchmark"] = {
            "status": "error",
            "seconds": bench["seconds"],
            "returncode": bench["returncode"],
            "json_output": str(json_output),
            "command_file": str(inspect_command_file),
            "output_log": str(inspect_output_log),
        }
        return result

    target = graph.get("target_weight")
    result["benchmark"] = {
        "status": "ok",
        "seconds": bench["seconds"],
        "returncode": bench["returncode"],
        "json_output": str(json_output),
        "command_file": str(inspect_command_file),
        "output_log": str(inspect_output_log),
        "timings": benchmark_json_timings(payload),
        "displayed_result": payload.get("displayed_result"),
        "delta_target": relative_delta(payload.get("displayed_result"), target)
        if target is not None
        else None,
        "total_samples": payload.get("total_samples"),
        "n_batches": payload.get("n_batches"),
        "minimal_integrand": payload.get("minimal_integrand"),
    }
    return result


def render_profiles(config: dict[str, Any], profile_names: list[str]) -> str:
    profiles = config["profiles"]
    card_defaults = placeholder_defaults_from_card(
        Path(__file__).resolve().parent / "bench.toml"
    )
    keys = sorted(
        {key for name in profile_names for key in profiles[name].get("defines", {})}
    )
    rows = []
    for key in keys:
        default = card_defaults.get(key)
        values = [
            str(profiles[name].get("defines", {}).get(key, ""))
            for name in profile_names
        ]
        if default is not None and all(value == default for value in values):
            continue
        rows.append(
            [
                color(key, "cyan"),
                *[
                    color(value, "dim")
                    if default is not None and value == default
                    else value
                    for value in values
                ],
            ]
        )
    if not rows:
        rows.append(
            [
                color("(all selected values match bench.toml defaults)", "dim"),
                *["" for _ in profile_names],
            ]
        )
    return table("Run Profiles", ["define", *profile_names], rows)


def process_characteristics(process_spec: str) -> dict[str, str]:
    final_state = re.search(r">\s*(.*?)\s*\|", process_spec)
    qed_squared = re.search(r"QED\^2==(\d+)", process_spec)
    qcd_squared = re.search(r"QCD\^2==(\d+)", process_spec)
    supergraph_loop = re.search(
        r"\[\s*\{\{\s*(\d+)\s*\}\}\s*QCD=(\d+)\s*\]", process_spec
    )
    amplitude_loop = re.search(r"\[\s*\{\s*(\d+)\s*\}\s*QCD=(\d+)\s*\]", process_spec)
    loop_match = supergraph_loop or amplitude_loop
    loop_kind = "SG" if supergraph_loop else "amp"
    return {
        "final_state": final_state.group(1) if final_state else "?",
        "qed_squared": qed_squared.group(1) if qed_squared else "?",
        "qcd_squared": qcd_squared.group(1) if qcd_squared else "?",
        "qcd_order": loop_match.group(2) if loop_match else "?",
        "loop_selector": f"{loop_match.group(1)} {loop_kind} loops"
        if loop_match
        else "?",
    }


def render_graphs(config: dict[str, Any], graph_names: list[str]) -> str:
    graphs = config["graphs"]
    rows = []
    fields = [
        ("path", lambda g: str(g["path"]).removeprefix("./")),
        ("model", lambda g: g["model"]),
        (
            "process / integrand",
            lambda g: f"{g.get('process', 'bench')} / {g.get('integrand', 'default')}",
        ),
        (
            "final state",
            lambda g: process_characteristics(g["process_spec"])["final_state"],
        ),
        ("QED²", lambda g: process_characteristics(g["process_spec"])["qed_squared"]),
        ("QCD²", lambda g: process_characteristics(g["process_spec"])["qcd_squared"]),
        (
            "QCD order",
            lambda g: process_characteristics(g["process_spec"])["qcd_order"],
        ),
        (
            "loop selector",
            lambda g: process_characteristics(g["process_spec"])["loop_selector"],
        ),
        ("x dim", lambda g: str(len(g["x"]))),
        (
            "AlphaLoop t_gen",
            lambda g: format_optional_reference(
                alpha_float(optional_alpha(g), "t_gen_s"), format_seconds
            ),
        ),
        (
            "AlphaLoop peak RSS",
            lambda g: format_optional_reference(
                alpha_peak_rss_bytes(optional_alpha(g)), format_bytes
            ),
        ),
        (
            "AlphaLoop disk",
            lambda g: format_optional_reference(
                alpha_disk_bytes(optional_alpha(g)), format_bytes
            ),
        ),
        (
            "AlphaLoop t_eval",
            lambda g: format_optional_reference(
                alpha_eval_seconds(optional_alpha(g)), format_seconds
            ),
        ),
        ("target", lambda g: "yes" if "target_weight" in g else color("no", "yellow")),
    ]
    for label, getter in fields:
        rows.append(
            [color(label, "cyan"), *[str(getter(graphs[name])) for name in graph_names]]
        )
    return table("Graphs", ["field", *graph_names], rows)


def render_report(
    payload: dict[str, Any], use_evaluator_timings: bool, compare_eval: bool
) -> str:
    profile_names = payload["selected_profiles"]
    graph_names = payload["selected_graphs"]
    binary_profile = payload.get("metadata", {}).get(
        "binary_profile", DEFAULT_BINARY_PROFILE
    )
    sections = [
        color("GammaLoop Benchmark Comparison", "bold", "green"),
        render_profiles(payload["config"], profile_names),
        render_graphs(payload["config"], graph_names),
        grouped_result_table(
            graph_names,
            profile_names,
            payload["results"],
            payload["config"]["graphs"],
            use_evaluator_timings,
            compare_eval,
            binary_profile,
            payload.get("metadata", {}),
        ),
    ]
    return "\n\n".join(sections)


def default_config_choices() -> tuple[str, str]:
    config_path = Path(__file__).resolve().parent / "bench_config.toml"
    try:
        config = load_toml(config_path)
    except (OSError, tomllib.TOMLDecodeError):
        return "(could not read default config)", "(could not read default config)"
    graphs = ", ".join(sorted(config.get("graphs", {}))) or "(none)"
    profiles = ", ".join(sorted(config.get("profiles", {}))) or "(none)"
    return graphs, profiles


def parse_args() -> argparse.Namespace:
    graph_choices, profile_choices = default_config_choices()
    parser = argparse.ArgumentParser(
        description="Benchmark GammaLoop example graphs against AlphaLoop reference timings.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config", default=DEFAULT_CONFIG, help="benchmark configuration TOML"
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help="JSON result path",
    )
    parser.add_argument(
        "--rerun",
        action="store_true",
        help="ignore an existing output JSON file and overwrite it with a fresh run",
    )
    parser.add_argument(
        "--graphs",
        nargs="+",
        default=DEFAULT_GRAPHS,
        help=f"graph names to run, or all; valid names in default config: {graph_choices}",
    )
    parser.add_argument(
        "--profiles",
        nargs="+",
        default=DEFAULT_PROFILES,
        help=f"run profile names to run, or all; valid names in default config: {profile_choices}",
    )
    parser.add_argument(
        "--max-gen-RAM",
        "--max-gen-ram",
        dest="max_gen_ram",
        type=parse_size_bytes,
        default=DEFAULT_MAX_GEN_RAM,
        help="generation RAM cap",
    )
    parser.add_argument(
        "--max-gen-time",
        type=parse_duration_seconds,
        default=DEFAULT_MAX_GEN_TIME,
        help="generation time cap",
    )
    parser.add_argument(
        "--profile",
        choices=["release", "debug", "dev-optim"],
        default=DEFAULT_BINARY_PROFILE,
        help="GammaLoop binary profile",
    )
    parser.add_argument(
        "--compare-eval",
        action="store_true",
        help="show relative distance to configured target weights",
    )
    parser.add_argument(
        "--use-evaluator-timings",
        action="store_true",
        help="compare the inspect evaluator row instead of Total",
    )
    parser.add_argument(
        "--target-profile-time",
        type=parse_duration_seconds,
        default=DEFAULT_TARGET_PROFILE_TIME,
        help="target runtime for inspect --bench",
    )
    parser.add_argument(
        "--no-colour",
        "--no-color",
        dest="no_colour",
        action="store_true",
        help="disable ANSI colours in progress messages and rendered tables",
    )
    return parser.parse_args()


def validate_selection(
    config: dict[str, Any], graphs: list[str], profiles: list[str]
) -> None:
    missing_graphs = [name for name in graphs if name not in config.get("graphs", {})]
    missing_profiles = [
        name for name in profiles if name not in config.get("profiles", {})
    ]
    if missing_graphs:
        raise SystemExit(f"Unknown graph(s): {', '.join(missing_graphs)}")
    if missing_profiles:
        raise SystemExit(f"Unknown profile(s): {', '.join(missing_profiles)}")


def expand_all_selection(
    requested: list[str], available: dict[str, Any], label: str
) -> list[str]:
    if any(name.lower() == "all" for name in requested):
        expanded = list(available)
        if not expanded:
            raise SystemExit(f"No {label}s are defined in the configuration.")
        return expanded
    return requested


def main() -> int:
    global USE_COLOUR
    args = parse_args()
    USE_COLOUR = not args.no_colour
    invocation_cwd = Path.cwd()
    bench_dir = Path(__file__).resolve().parent
    repo_root = bench_dir.parents[2]
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = (
            bench_dir / output_path
            if args.output == DEFAULT_OUTPUT
            else invocation_cwd / output_path
        )

    if output_path.exists() and not args.rerun:
        payload = json.loads(output_path.read_text(encoding="utf-8"))
        compare_eval = args.compare_eval or bool(
            payload.get("metadata", {}).get("compare_eval")
        )
        print(render_report(payload, args.use_evaluator_timings, compare_eval))
        return 0

    config_path = Path(args.config)
    if not config_path.is_absolute():
        candidate = invocation_cwd / config_path
        config_path = candidate if candidate.exists() else bench_dir / config_path
    config = load_toml(config_path)
    args.graphs = expand_all_selection(args.graphs, config.get("graphs", {}), "graph")
    args.profiles = expand_all_selection(
        args.profiles, config.get("profiles", {}), "profile"
    )
    validate_selection(config, args.graphs, args.profiles)

    payload: dict[str, Any] = {
        "metadata": {
            "created_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "config_path": str(config_path),
            "binary_profile": args.profile,
            "max_gen_ram_bytes": args.max_gen_ram,
            "max_gen_time_s": args.max_gen_time,
            "target_profile_time_s": args.target_profile_time,
            "compare_eval": args.compare_eval,
            "rerun": args.rerun,
        },
        "selected_graphs": args.graphs,
        "selected_profiles": args.profiles,
        "config": config,
        "results": {},
    }

    interrupt_state = InterruptState()
    stop_remaining = False
    for graph_name in args.graphs:
        payload["results"][graph_name] = {}
        for profile_name in args.profiles:
            if stop_remaining:
                break
            print(
                color(f"running {graph_name}/{profile_name}", "bold", "cyan"),
                file=sys.stderr,
            )
            try:
                payload["results"][graph_name][profile_name] = run_case(
                    graph_name,
                    config["graphs"][graph_name],
                    profile_name,
                    config["profiles"][profile_name],
                    args,
                    bench_dir,
                    repo_root,
                    interrupt_state,
                )
            except KeyboardInterrupt:
                if interrupt_state.register():
                    stop_remaining = True
                    break
                payload["results"][graph_name][profile_name] = {
                    "status": "interrupted",
                    "generation": {"status": "interrupted"},
                    "disk_size_bytes": None,
                    "benchmark": None,
                    "error_log": None,
                }
            if payload["results"][graph_name][profile_name].get("abort_all_requested"):
                stop_remaining = True
                break
        if stop_remaining:
            payload["metadata"]["interrupted"] = True
            break

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(render_report(payload, args.use_evaluator_timings, args.compare_eval))
    print(color(f"\nSaved JSON results to {output_path}", "green"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
