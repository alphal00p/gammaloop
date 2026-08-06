#!/usr/bin/env python3
"""Plot JSON output produced by `gammaloop approach`."""

from __future__ import annotations

import argparse
import json
import math
import re
import shlex
import sys
import textwrap
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


@dataclass
class ResultFile:
    path: Path
    data: dict[str, Any]


@dataclass
class Series:
    label: str
    t_values: list[float]
    values: list[float]


class NonFiniteWarnings:
    def __init__(self) -> None:
        self._items: OrderedDict[tuple[Any, ...], dict[str, Any]] = OrderedDict()

    def add(
        self,
        result: ResultFile,
        point: dict[str, Any],
        label: str,
        component: str,
        raw_value: Any,
    ) -> None:
        key = (
            result.path,
            int(point.get("axis_index", -1)),
            int(point.get("axis_point_index", point.get("index", -1))),
            repr(point.get("t")),
            tuple(repr(value) for value in point.get("point", [])),
        )
        item = self._items.setdefault(
            key,
            {
                "result": result,
                "point": point,
                "fields": OrderedDict(),
            },
        )
        item["fields"][f"{label}.{component}"] = raw_value_label(raw_value)

    def emit(self) -> None:
        for item in self._items.values():
            result = item["result"]
            point = item["point"]
            axis_index = int(point.get("axis_index", -1))
            axis_point_index = int(
                point.get("axis_point_index", point.get("index", -1))
            )
            t_value = point.get("t")
            coordinate = pretty_exact_vector(point.get("point", []))
            print("plot_approach_result.py: warning", file=sys.stderr)
            print(
                "  Replaced non-finite/null value(s) with 0.0 for plotting.",
                file=sys.stderr,
            )
            print(f"  Source: {display_path(result.path)}", file=sys.stderr)
            print("  Sample:", file=sys.stderr)
            print(f"    axis: {axis_index + 1}", file=sys.stderr)
            print(f"    axis point: {axis_point_index + 1}", file=sys.stderr)
            print(f"    t: {exact_number(t_value)}", file=sys.stderr)
            print("    x:", file=sys.stderr)
            for line in textwrap.wrap(
                coordinate,
                width=100,
                break_long_words=False,
            ):
                print(f"      {line}", file=sys.stderr)
            print("  Affected fields:", file=sys.stderr)
            for line in grouped_field_lines(item["fields"]):
                print(f"    - {line}", file=sys.stderr)
            print("  Inspect command:", file=sys.stderr)
            for line in inspect_command_lines(result, point):
                print(f"    {line}", file=sys.stderr)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot one or more gammaloop approach result JSON files."
    )
    parser.add_argument(
        "results", nargs="+", type=Path, help="Approach result JSON file(s)"
    )
    parser.add_argument("--output", required=True, type=Path, help="Output PDF path")
    parser.add_argument(
        "--component",
        choices=("real", "imag", "abs"),
        default="abs",
        help="Complex component to plot",
    )
    parser.add_argument(
        "--t-branch",
        choices=("both", "positive", "negative"),
        default="both",
        help="Plot both signed approach branches or only one branch",
    )
    parser.add_argument(
        "--x-log-scale",
        action="store_true",
        help="Plot |t| on a logarithmic x-axis; requires --t-branch positive or negative",
    )
    parser.add_argument(
        "--fit-log-slope",
        action="store_true",
        help="Append the fitted log-log power over the visible x-range to each legend label",
    )
    y_scale_group = parser.add_mutually_exclusive_group()
    y_scale_group.add_argument(
        "--y-log-scale",
        dest="y_log_scale",
        action="store_true",
        help="Use a logarithmic y-axis; real/imag values are plotted by absolute value and negative values are dashed (default)",
    )
    y_scale_group.add_argument(
        "--linear-y-scale",
        dest="y_log_scale",
        action="store_false",
        help="Use a signed linear y-axis instead of the default logarithmic y-axis",
    )
    parser.set_defaults(y_log_scale=True)
    parser.add_argument(
        "--include-contribution",
        action="append",
        default=[],
        help="Regex for labels to include; can be repeated (default: total_weight only)",
    )
    parser.add_argument(
        "--exclude-contribution",
        action="append",
        default=[],
        help="Regex for labels to exclude; can be repeated",
    )
    parser.add_argument(
        "--combine-plots",
        action="store_true",
        help="Overlay compatible result files on the same plot",
    )
    parser.add_argument(
        "--combine-axes",
        action="store_true",
        help="Overlay all approach axes on the same plot",
    )
    parser.add_argument(
        "--sum-lmb-samples-per-cut",
        action="store_true",
        help="Aggregate contribution curves over all lmb_sample_id values within each cut",
    )
    parser.add_argument(
        "--threshold-decomposition",
        choices=("summary", "components", "all"),
        help="Add threshold-counterterm summary curves, individual weighted components, or both",
    )
    parser.add_argument(
        "--title",
        help="Replace the generated plot title",
    )
    parser.add_argument(
        "--hide-info-box",
        action="store_true",
        help="Omit the process, integrand, kinematics, and source summary box",
    )
    parser.add_argument(
        "--x-range",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        help="Visible x-axis range to show",
    )
    parser.add_argument(
        "--y-range",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        help="Visible y-axis range to show; with logarithmic y-scale both values must be positive",
    )
    parser.add_argument(
        "--axis-id",
        action="append",
        default=[],
        metavar="ID[,ID...]",
        help="Only plot the displayed approach axis ID(s); can be repeated or comma-separated",
    )
    return parser.parse_args()


def validate_range_option(
    option_name: str, values: list[float] | None, *, log_scale: bool = False
) -> None:
    if values is None:
        return
    lower, upper = values
    if not math.isfinite(lower) or not math.isfinite(upper):
        raise ValueError(f"{option_name} values must be finite.")
    if lower >= upper:
        raise ValueError(
            f"{option_name} expects MIN < MAX, got {lower:.17g} >= {upper:.17g}."
        )
    if log_scale and lower <= 0.0:
        raise ValueError(
            f"{option_name} lower bound must be positive with logarithmic scale, got {lower:.17g}."
        )


def validate_args(args: argparse.Namespace) -> None:
    validate_range_option("--x-range", args.x_range, log_scale=args.x_log_scale)
    validate_range_option("--y-range", args.y_range, log_scale=args.y_log_scale)
    if args.x_log_scale and args.t_branch == "both":
        raise ValueError(
            "--x-log-scale requires --t-branch positive or --t-branch negative."
        )
    if args.fit_log_slope and not (args.x_log_scale and args.y_log_scale):
        raise ValueError("--fit-log-slope requires logarithmic x and y axes.")
    args.selected_axis_indices = parse_selected_axis_indices(args.axis_id)


def parse_selected_axis_indices(raw_axis_ids: list[str]) -> set[int] | None:
    if not raw_axis_ids:
        return None
    selected = set()
    for raw_axis_id in raw_axis_ids:
        for token in raw_axis_id.split(","):
            token = token.strip()
            if not token:
                continue
            try:
                axis_id = int(token)
            except ValueError as error:
                raise ValueError(
                    f"--axis-id expects positive integer IDs, got {token!r}."
                ) from error
            if axis_id < 1:
                raise ValueError(
                    f"--axis-id uses displayed one-based axis IDs; got {axis_id}."
                )
            selected.add(axis_id - 1)
    if not selected:
        raise ValueError("--axis-id was supplied but no axis IDs were parsed.")
    return selected


def load_result(path: Path) -> ResultFile:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if "points" not in data:
        raise ValueError(
            f"{path} is not an approach result JSON file: missing 'points'"
        )
    return ResultFile(path=path, data=data)


def compile_patterns(patterns: Iterable[str]) -> list[re.Pattern[str]]:
    return [re.compile(pattern) for pattern in patterns]


def label_is_selected(
    label: str,
    include_patterns: list[re.Pattern[str]],
    exclude_patterns: list[re.Pattern[str]],
) -> bool:
    if include_patterns and not any(
        pattern.search(label) for pattern in include_patterns
    ):
        return False
    return not any(pattern.search(label) for pattern in exclude_patterns)


def axis_indices(
    result: ResultFile, selected_axis_indices: set[int] | None = None
) -> list[int]:
    indices = sorted(
        {int(point["axis_index"]) for point in result.data.get("points", [])}
    )
    if selected_axis_indices is None:
        return indices
    return [axis_index for axis_index in indices if axis_index in selected_axis_indices]


def axis_label(result: ResultFile, axis_index: int) -> str:
    axes = result.data.get("axes", [])
    if 0 <= axis_index < len(axes):
        vector = ", ".join(f"{float(component):.4g}" for component in axes[axis_index])
        return f"axis {axis_index + 1}: [{vector}]"
    return f"axis {axis_index + 1}"


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(Path.cwd().resolve()))
    except ValueError:
        return str(path)


def raw_value_label(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, float) and math.isnan(value):
        return "NaN"
    if isinstance(value, float) and math.isinf(value):
        return "+inf" if value > 0.0 else "-inf"
    return repr(value)


def exact_number(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return raw_value_label(value)
    return f"{number:.17g}"


def exact_vector(values: list[Any]) -> str:
    return "[" + ",".join(exact_number(value) for value in values) + "]"


def pretty_exact_vector(values: list[Any]) -> str:
    return "[" + ", ".join(exact_number(value) for value in values) + "]"


def grouped_field_lines(fields: OrderedDict[str, str]) -> list[str]:
    grouped: OrderedDict[str, OrderedDict[str, str]] = OrderedDict()
    for field, value in fields.items():
        if "." in field:
            label, component = field.rsplit(".", 1)
        else:
            label = field
            component = "value"
        grouped.setdefault(label, OrderedDict())[component] = value
    return [
        f"{label}: "
        + ", ".join(f"{component}={value}" for component, value in components.items())
        for label, components in grouped.items()
    ]


def sanitized_float(
    value: Any,
    warnings: NonFiniteWarnings,
    result: ResultFile,
    point: dict[str, Any],
    label: str,
    component: str,
) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        warnings.add(result, point, label, component, value)
        return 0.0
    if not math.isfinite(number):
        warnings.add(result, point, label, component, value)
        return 0.0
    return number


def sanitized_complex_component(
    value: dict[str, Any],
    component: str,
    warnings: NonFiniteWarnings,
    result: ResultFile,
    point: dict[str, Any],
    label: str,
) -> float:
    re_value = sanitized_float(
        value.get("re", 0.0), warnings, result, point, label, "re"
    )
    im_value = sanitized_float(
        value.get("im", 0.0), warnings, result, point, label, "im"
    )
    if component == "real":
        return re_value
    if component == "imag":
        return im_value
    return math.hypot(re_value, im_value)


def value_or_none(
    value: Any,
    component: str,
    warnings: NonFiniteWarnings,
    result: ResultFile,
    point: dict[str, Any],
    label: str,
) -> float | None:
    if not isinstance(value, dict):
        return None
    return sanitized_complex_component(value, component, warnings, result, point, label)


def inspect_command_segments(
    result: ResultFile, point: dict[str, Any]
) -> list[list[str]]:
    process = result.data.get("process", {}).get("name")
    integrand = result.data.get("integrand", {}).get("name")
    command_head = ["inspect"]
    if process is not None:
        command_head.extend(["-p", str(process)])
    if integrand is not None:
        command_head.extend(["-i", str(integrand)])
    segments = [command_head, ["-x", exact_vector(point.get("point", [])).strip("[]")]]
    if result.data.get("space") == "momentum":
        segments.append(["--momentum-space"])
    approach_command = result.data.get("command", {})
    discrete_dim = approach_command.get("discrete_dim", [])
    if discrete_dim:
        segments.append(
            ["--discrete-dim", ",".join(str(value) for value in discrete_dim)]
        )
    graph_id = approach_command.get("graph_id")
    if graph_id is not None:
        segments.append(["--graph-id", str(graph_id)])
    orientation_id = approach_command.get("orientation_id")
    if orientation_id is not None:
        segments.append(["--orientation-id", str(orientation_id)])
    return segments


def inspect_command(result: ResultFile, point: dict[str, Any]) -> str:
    return " ".join(
        shlex.quote(part)
        for segment in inspect_command_segments(result, point)
        for part in segment
    )


def inspect_command_lines(result: ResultFile, point: dict[str, Any]) -> list[str]:
    lines = []
    segments = inspect_command_segments(result, point)
    for index, segment in enumerate(segments):
        line = " ".join(shlex.quote(part) for part in segment)
        if index + 1 < len(segments):
            line = f"{line} \\"
        if index > 0:
            line = f"  {line}"
        lines.append(line)
    return lines


def contribution_group_label(
    contribution: dict[str, Any], sum_lmb_samples: bool
) -> str:
    graph = contribution.get("graph_name")
    if graph is None:
        graph = f"#{contribution.get('graph_id')}"
    parts = [
        str(contribution.get("contribution", "contribution")),
        str(graph),
        f"cut{contribution.get('cut_id')}",
    ]
    orientation_id = contribution.get("orientation_id")
    if orientation_id is not None:
        parts.append(f"ori{orientation_id}")
    if sum_lmb_samples:
        parts.append("lmb=sum")
    else:
        lmb_sample_id = contribution.get("lmb_sample_id")
        if lmb_sample_id is not None:
            parts.append(f"lmb{lmb_sample_id}")
    return " ".join(parts)


def complex_value(
    value: Any,
    warnings: NonFiniteWarnings,
    result: ResultFile,
    point: dict[str, Any],
    label: str,
) -> complex:
    if not isinstance(value, dict):
        return 0j
    return complex(
        sanitized_float(value.get("re", 0.0), warnings, result, point, label, "re"),
        sanitized_float(value.get("im", 0.0), warnings, result, point, label, "im"),
    )


def threshold_component_catalog(result: ResultFile) -> dict[tuple[int, int], str]:
    catalog: dict[tuple[int, int], str] = {}
    entries = result.data.get("integrand", {}).get("threshold_counterterms", [])
    if not isinstance(entries, list):
        return catalog
    kind_labels = {
        "local": "L",
        "integrated": "I",
        "local_local": "LL",
        "local_integrated": "LI",
        "integrated_local": "IL",
        "integrated_integrated": "II",
    }
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        graph_id = int(entry.get("graph_id", -1))
        registry = entry.get("registry", entry)
        if not isinstance(registry, dict):
            continue
        graph_name = str(registry.get("graph_name", f"#{graph_id}"))
        variants = {
            int(variant.get("variant_id", -1)): str(variant.get("name", "?"))
            for variant in registry.get("variants", [])
            if isinstance(variant, dict)
        }
        for component in registry.get("components", []):
            if not isinstance(component, dict):
                continue
            component_id = int(component.get("component_id", -1))
            kind = str(component.get("kind", "component"))
            kind_label = kind_labels.get(kind, kind)
            variant_label = "×".join(
                f"v{variant_id}:{variants.get(int(variant_id), '?')}"
                for variant_id in component.get("variant_ids", [])
            )
            cut_group_id = component.get("cut_group_id")
            cut_label = "" if cut_group_id is None else f" cg{cut_group_id}"
            catalog[(graph_id, component_id)] = (
                f"threshold:component c{component_id} {kind_label} {variant_label} "
                f"[{graph_name}{cut_label}]"
            )
    return catalog


def threshold_decomposition_values(
    result: ResultFile,
    point: dict[str, Any],
    mode: str,
    warnings: NonFiniteWarnings,
) -> list[tuple[str, complex]]:
    evaluation = point.get("evaluation", {})
    events = evaluation.get("events", []) if isinstance(evaluation, dict) else []
    original = 0j
    counterterms = 0j
    components: dict[str, complex] = defaultdict(complex)
    catalog = threshold_component_catalog(result)
    found = False
    for event in events:
        if not isinstance(event, dict):
            continue
        decomposition = event.get("threshold_counterterms")
        if not isinstance(decomposition, dict):
            continue
        found = True
        original += complex_value(
            decomposition.get("original"), warnings, result, point, "threshold:original"
        )
        graph_id = int(event.get("graph_id", -1))
        for component in decomposition.get("components", []):
            if not isinstance(component, dict):
                continue
            component_id = int(component.get("component_id", -1))
            label = catalog.get(
                (graph_id, component_id),
                f"threshold:component c{component_id} [graph#{graph_id}]",
            )
            weighted = complex_value(
                component.get("weighted"), warnings, result, point, label
            )
            counterterms += weighted
            components[label] += weighted
    if not found:
        return []

    values: list[tuple[str, complex]] = []
    if mode in ("summary", "all"):
        values.extend(
            [
                ("threshold:original", original),
                ("threshold:counterterm_sum", counterterms),
                ("threshold:decomposition_total", original + counterterms),
            ]
        )
    if mode in ("components", "all"):
        values.extend(sorted(components.items()))
    return values


def selected_t_value(
    point: dict[str, Any], branch: str, x_log_scale: bool
) -> float | None:
    value = float(point["t"])
    if branch == "positive" and value <= 0.0:
        return None
    if branch == "negative" and value >= 0.0:
        return None
    return abs(value) if x_log_scale else value


def collect_series(
    result: ResultFile,
    axis_index: int,
    component: str,
    include_patterns: list[re.Pattern[str]],
    exclude_patterns: list[re.Pattern[str]],
    sum_lmb_samples: bool,
    warnings: NonFiniteWarnings,
    t_branch: str,
    x_log_scale: bool,
    threshold_decomposition: str | None,
    label_prefix: str = "",
) -> list[Series]:
    points = [
        point
        for point in result.data.get("points", [])
        if int(point.get("axis_index", -1)) == axis_index
    ]
    points.sort(
        key=lambda point: int(point.get("axis_point_index", point.get("index", 0)))
    )

    series_values: OrderedDict[str, list[tuple[float, float]]] = OrderedDict()
    for point in points:
        if point.get("status") != "evaluated":
            continue
        evaluation = point.get("evaluation")
        if not isinstance(evaluation, dict):
            continue
        t_value = selected_t_value(point, t_branch, x_log_scale)
        if t_value is None:
            continue

        base_values: list[tuple[str, dict[str, Any]]] = [
            ("total_weight", evaluation.get("total_weight", {})),
            ("event_weight_sum", evaluation.get("event_weight_sum", {})),
        ]
        for weight_name, weight_value in sorted(
            evaluation.get("additional_weight_sums", {}).items()
        ):
            base_values.append((f"additional:{weight_name}", weight_value))

        for label, serialized_value in base_values:
            if not label_is_selected(label, include_patterns, exclude_patterns):
                continue
            y_value = value_or_none(
                serialized_value, component, warnings, result, point, label
            )
            if y_value is not None and math.isfinite(y_value):
                series_values.setdefault(f"{label_prefix}{label}", []).append(
                    (t_value, y_value)
                )

        if threshold_decomposition is not None:
            for label, value in threshold_decomposition_values(
                result, point, threshold_decomposition, warnings
            ):
                if not label_is_selected(label, include_patterns, exclude_patterns):
                    continue
                if component == "real":
                    y_value = value.real
                elif component == "imag":
                    y_value = value.imag
                else:
                    y_value = abs(value)
                if math.isfinite(y_value):
                    series_values.setdefault(f"{label_prefix}{label}", []).append(
                        (t_value, y_value)
                    )

        contribution_accumulator: dict[str, complex] = defaultdict(complex)
        for contribution in evaluation.get("contributions", []):
            label = contribution_group_label(contribution, sum_lmb_samples)
            full_label = f"contribution:{label}"
            if not label_is_selected(full_label, include_patterns, exclude_patterns):
                continue
            weight = contribution.get("weight", {})
            if not isinstance(weight, dict):
                weight = {}
            contribution_accumulator[label] += complex(
                sanitized_float(
                    weight.get("re", 0.0), warnings, result, point, full_label, "re"
                ),
                sanitized_float(
                    weight.get("im", 0.0), warnings, result, point, full_label, "im"
                ),
            )

        for label, value in sorted(contribution_accumulator.items()):
            full_label = f"{label_prefix}contribution:{label}"
            if component == "real":
                y_value = value.real
            elif component == "imag":
                y_value = value.imag
            else:
                y_value = abs(value)
            if math.isfinite(y_value):
                series_values.setdefault(full_label, []).append((t_value, y_value))

    series = []
    for label, values in series_values.items():
        values.sort(key=lambda item: item[0])
        series.append(
            Series(
                label=label,
                t_values=[item[0] for item in values],
                values=[item[1] for item in values],
            )
        )
    return series


def finite_abs(value: float) -> float:
    if value == 0.0:
        return math.nan
    return abs(value)


def log_slope(series: Series, x_range: list[float] | None = None) -> float | None:
    points = [
        (math.log(x), math.log(abs(y)))
        for x, y in zip(series.t_values, series.values)
        if x > 0.0
        and y != 0.0
        and math.isfinite(x)
        and math.isfinite(y)
        and (x_range is None or x_range[0] <= x <= x_range[1])
    ]
    if len(points) < 2:
        return None
    mean_x = sum(x for x, _ in points) / len(points)
    mean_y = sum(y for _, y in points) / len(points)
    denominator = sum((x - mean_x) ** 2 for x, _ in points)
    if denominator == 0.0:
        return None
    return sum((x - mean_x) * (y - mean_y) for x, y in points) / denominator


def plot_series(
    ax: plt.Axes,
    series: Series,
    component: str,
    y_log_scale: bool,
    fit_log_slope: bool,
    x_range: list[float] | None,
) -> None:
    if not series.t_values:
        return
    label = series.label
    if fit_log_slope and (slope := log_slope(series, x_range)) is not None:
        label = f"{label} (p={slope:+.2f})"
    if not y_log_scale or component == "abs":
        y_values = [
            finite_abs(value) if y_log_scale else value for value in series.values
        ]
        ax.plot(series.t_values, y_values, linewidth=1.4, label=label)
        return

    positive = [
        finite_abs(value) if value >= 0.0 else math.nan for value in series.values
    ]
    negative = [
        finite_abs(value) if value < 0.0 else math.nan for value in series.values
    ]
    positive_line = None
    if any(math.isfinite(value) for value in positive):
        (positive_line,) = ax.plot(
            series.t_values,
            positive,
            linewidth=1.4,
            label=label,
        )
    if any(math.isfinite(value) for value in negative):
        color = positive_line.get_color() if positive_line is not None else None
        ax.plot(
            series.t_values,
            negative,
            linewidth=1.4,
            linestyle="--",
            color=color,
            label=f"{label} (negative)",
        )


def compact_vector(values: list[Any], limit: int = 8) -> str:
    formatted = [f"{float(value):.4g}" for value in values]
    if len(formatted) > limit:
        formatted = formatted[:limit] + ["..."]
    return "[" + ", ".join(formatted) + "]"


def info_box_text(result: ResultFile) -> str:
    process = result.data.get("process", {})
    integrand = result.data.get("integrand", {})
    spacing = result.data.get("spacing", {})
    return "\n".join(
        [
            f"process: {process.get('name', process.get('id', '?'))}",
            f"integrand: {integrand.get('name', '?')} ({integrand.get('kind', '?')})",
            f"space: {result.data.get('space', '?')}",
            f"base: {compact_vector(result.data.get('base_point', []))}",
            f"spacing: {spacing.get('kind', '?')} n={spacing.get('n_points', '?')}",
            f"source: {result.path.name}",
        ]
    )


def decorate_axes(
    ax: plt.Axes,
    title: str,
    result: ResultFile,
    component: str,
    y_log_scale: bool,
    x_log_scale: bool,
    hide_info_box: bool,
) -> None:
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("|approach parameter t|" if x_log_scale else "approach parameter t")
    if x_log_scale:
        ax.set_xscale("log")
    ylabel = "weight magnitude" if component == "abs" else f"{component} weight"
    if y_log_scale:
        if component != "abs":
            ylabel = f"|{ylabel}|"
        ax.set_yscale("log")
    ax.set_ylabel(ylabel)
    if not x_log_scale:
        ax.axvline(0.0, color="0.35", linewidth=0.9, alpha=0.7)
    ax.grid(True, which="major", color="0.88", linewidth=0.8)
    ax.grid(True, which="minor", color="0.94", linewidth=0.5)
    if not hide_info_box:
        ax.text(
            0.012,
            0.988,
            info_box_text(result),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.5,
            bbox={
                "boxstyle": "round,pad=0.36",
                "facecolor": "white",
                "edgecolor": "0.7",
                "alpha": 0.92,
            },
        )


def apply_axis_ranges(ax: plt.Axes, args: argparse.Namespace) -> None:
    if args.x_range is not None:
        ax.set_xlim(args.x_range[0], args.x_range[1])
    if args.y_range is not None:
        ax.set_ylim(args.y_range[0], args.y_range[1])


def draw_page(
    pdf: PdfPages,
    result_axis_pairs: list[tuple[ResultFile, int]],
    args: argparse.Namespace,
    include_patterns: list[re.Pattern[str]],
    exclude_patterns: list[re.Pattern[str]],
    warnings: NonFiniteWarnings,
) -> None:
    fig, ax = plt.subplots(figsize=(13.5, 8.2), constrained_layout=False)
    has_series = False
    multi_result = len({pair[0].path for pair in result_axis_pairs}) > 1
    multi_axis = len({pair[1] for pair in result_axis_pairs}) > 1

    for result, axis_index in result_axis_pairs:
        prefix_parts = []
        if multi_result or args.combine_plots:
            prefix_parts.append(result.path.stem)
        if multi_axis or args.combine_axes:
            prefix_parts.append(f"axis {axis_index + 1}")
        label_prefix = (" / ".join(prefix_parts) + " / ") if prefix_parts else ""
        series = collect_series(
            result,
            axis_index,
            args.component,
            include_patterns,
            exclude_patterns,
            args.sum_lmb_samples_per_cut,
            warnings,
            args.t_branch,
            args.x_log_scale,
            args.threshold_decomposition,
            label_prefix=label_prefix,
        )
        for item in series:
            has_series = True
            plot_series(
                ax,
                item,
                args.component,
                args.y_log_scale,
                args.fit_log_slope,
                args.x_range,
            )

    first_result, first_axis = result_axis_pairs[0]
    title = args.title or axis_label(first_result, first_axis)
    if args.title is None and len(result_axis_pairs) > 1:
        title = "combined approach curves"
    decorate_axes(
        ax,
        title,
        first_result,
        args.component,
        args.y_log_scale,
        args.x_log_scale,
        args.hide_info_box,
    )
    apply_axis_ranges(ax, args)
    if has_series:
        _, labels = ax.get_legend_handles_labels()
        if len(labels) > 8:
            legend_columns = min(5, max(2, math.ceil(len(labels) / 10)))
            ax.legend(
                loc="upper center",
                bbox_to_anchor=(0.5, -0.14),
                fontsize=5.8,
                frameon=True,
                ncols=legend_columns,
                borderaxespad=0.0,
            )
            fig.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.28)
        else:
            ax.legend(loc="upper right", fontsize=7.0, frameon=True, ncols=1)
            fig.tight_layout()
    else:
        ax.text(
            0.5,
            0.5,
            "No evaluated series matched the requested filters.",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
            color="0.35",
        )
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def page_groups(
    results: list[ResultFile], args: argparse.Namespace
) -> list[list[tuple[ResultFile, int]]]:
    if args.combine_plots and args.combine_axes:
        pairs = [
            (result, axis_index)
            for result in results
            for axis_index in axis_indices(result, args.selected_axis_indices)
        ]
        return [pairs] if pairs else []

    if args.combine_plots:
        all_axes = sorted(
            {
                axis_index
                for result in results
                for axis_index in axis_indices(result, args.selected_axis_indices)
            }
        )
        return [
            [
                (result, axis_index)
                for result in results
                if axis_index in axis_indices(result, args.selected_axis_indices)
            ]
            for axis_index in all_axes
        ]

    if args.combine_axes:
        return [
            [
                (result, axis_index)
                for axis_index in axis_indices(result, args.selected_axis_indices)
            ]
            for result in results
            if axis_indices(result, args.selected_axis_indices)
        ]

    return [
        [(result, axis_index)]
        for result in results
        for axis_index in axis_indices(result, args.selected_axis_indices)
    ]


def main() -> int:
    args = parse_args()
    try:
        validate_args(args)
        results = [load_result(path) for path in args.results]
        warnings = NonFiniteWarnings()
        include_contributions = list(args.include_contribution)
        if not include_contributions:
            include_contributions.append(r"^total_weight$")
            if args.threshold_decomposition is not None:
                include_contributions.append(r"^threshold:")
        include_patterns = compile_patterns(include_contributions)
        exclude_patterns = compile_patterns(args.exclude_contribution)
        groups = page_groups(results, args)
        if not groups:
            raise ValueError(
                "No axes found in the supplied approach result JSON files."
            )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with PdfPages(args.output) as pdf:
            for group in groups:
                draw_page(
                    pdf, group, args, include_patterns, exclude_patterns, warnings
                )
        warnings.emit()
        print(f"plot_approach_result.py: PDF created at {display_path(args.output)}")
    except Exception as error:  # noqa: BLE001 - this script should print clean CLI errors.
        print(f"plot_approach_result.py: error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
