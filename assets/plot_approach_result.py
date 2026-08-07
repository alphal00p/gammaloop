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
from matplotlib.lines import Line2D


THRESHOLD_KIND_LABELS = {
    "local": "L",
    "integrated": "I",
    "local_local": "LL",
    "local_integrated": "LI",
    "integrated_local": "IL",
    "integrated_integrated": "II",
}
THRESHOLD_KIND_ORDER = {kind: index for index, kind in enumerate(THRESHOLD_KIND_LABELS)}


@dataclass
class ResultFile:
    path: Path
    data: dict[str, Any]
    label: str | None = None


@dataclass
class Series:
    label: str
    t_values: list[float]
    values: list[float]


@dataclass(frozen=True)
class ThresholdVariant:
    graph_id: int
    graph_name: str
    variant_id: int
    name: str
    cut_group_id: int | None
    associations: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...]
    side: str
    subspace: tuple[int, ...]

    def cut_label(self) -> str:
        cuts = sorted({cut for cut, _ in self.associations})
        return "/".join(f"c({','.join(map(str, cut))})" for cut in cuts) or "c(?)"

    def compact_label(self) -> str:
        thresholds = sorted({threshold for _, threshold in self.associations})
        eta = (
            "/".join(f"η({','.join(map(str, threshold))})" for threshold in thresholds)
            or "η(?)"
        )
        side = self.side[:1].upper() or "?"
        subspace = ",".join(map(str, self.subspace)) or "max"
        return f"v{self.variant_id} {self.name} · {side} {eta} · S[{subspace}]"


@dataclass(frozen=True)
class ThresholdComponent:
    graph_id: int
    component_id: int
    cut_group_id: int | None
    kind: str
    variant_ids: tuple[int, ...]

    @property
    def key(self) -> tuple[int, int]:
        return self.graph_id, self.component_id


@dataclass
class ThresholdRegistry:
    variants: dict[tuple[int, int], ThresholdVariant]
    components: dict[tuple[int, int], ThresholdComponent]

    @classmethod
    def from_result(cls, result: ResultFile) -> ThresholdRegistry:
        variants = {}
        components = {}
        entries = result.data.get("integrand", {}).get("threshold_counterterms", [])
        if not isinstance(entries, list):
            return cls(variants, components)
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            graph_id = int(entry.get("graph_id", -1))
            registry = entry.get("registry", entry)
            if not isinstance(registry, dict):
                continue
            graph_name = str(registry.get("graph_name", f"#{graph_id}"))
            for raw_variant in registry.get("variants", []):
                if not isinstance(raw_variant, dict):
                    continue
                variant_id = int(raw_variant.get("variant_id", -1))
                associations = []
                for association in raw_variant.get("associations", []):
                    if not isinstance(association, dict):
                        continue
                    associations.append(
                        (
                            tuple(
                                int(edge) for edge in association.get("cut_edges", [])
                            ),
                            tuple(
                                int(edge)
                                for edge in association.get("threshold_edges", [])
                            ),
                        )
                    )
                raw_cut_group_id = raw_variant.get("cut_group_id", -1)
                variants[(graph_id, variant_id)] = ThresholdVariant(
                    graph_id=graph_id,
                    graph_name=graph_name,
                    variant_id=variant_id,
                    name=str(raw_variant.get("name", "?")),
                    cut_group_id=(
                        None if raw_cut_group_id is None else int(raw_cut_group_id)
                    ),
                    associations=tuple(associations),
                    side=str(raw_variant.get("side", "?")),
                    subspace=tuple(
                        int(edge) for edge in raw_variant.get("resolved_subspace") or []
                    ),
                )
            for raw_component in registry.get("components", []):
                if not isinstance(raw_component, dict):
                    continue
                raw_cut_group_id = raw_component.get("cut_group_id", -1)
                component = ThresholdComponent(
                    graph_id=graph_id,
                    component_id=int(raw_component.get("component_id", -1)),
                    cut_group_id=(
                        None if raw_cut_group_id is None else int(raw_cut_group_id)
                    ),
                    kind=str(raw_component.get("kind", "component")),
                    variant_ids=tuple(
                        int(variant_id)
                        for variant_id in raw_component.get("variant_ids", [])
                    ),
                )
                components[component.key] = component
        return cls(variants, components)

    def variant(self, graph_id: int, variant_id: int) -> ThresholdVariant:
        return self.variants.get(
            (graph_id, variant_id),
            ThresholdVariant(
                graph_id,
                f"#{graph_id}",
                variant_id,
                "?",
                -1,
                (),
                "?",
                (),
            ),
        )

    def cut_group_label(self, graph_id: int, cut_group_id: int | None) -> str:
        group_variants = [
            variant
            for (candidate_graph_id, _), variant in self.variants.items()
            if candidate_graph_id == graph_id and variant.cut_group_id == cut_group_id
        ]
        graph_name = group_variants[0].graph_name if group_variants else f"#{graph_id}"
        cuts = sorted(
            {cut for variant in group_variants for cut, _ in variant.associations}
        )
        cut_label = "/".join(f"c({','.join(map(str, cut))})" for cut in cuts) or "c(?)"
        if cut_group_id is None:
            return f"{graph_name} · {cut_label}"
        return f"{graph_name} · cg{cut_group_id} · {cut_label}"

    def component_label(self, component: ThresholdComponent) -> str:
        variants = "×".join(
            self.variant(component.graph_id, variant_id).compact_label()
            for variant_id in component.variant_ids
        )
        kind = THRESHOLD_KIND_LABELS.get(component.kind, component.kind)
        return (
            f"threshold:component c{component.component_id} {kind} {variants} "
            f"[{self.cut_group_label(component.graph_id, component.cut_group_id)}]"
        )

    def observed_components(
        self, result: ResultFile, axis_index: int
    ) -> set[tuple[int, int]]:
        observed = set()
        for point in result.data.get("points", []):
            if (
                int(point.get("axis_index", -1)) != axis_index
                or point.get("status") != "evaluated"
            ):
                continue
            evaluation = point.get("evaluation", {})
            for event in evaluation.get("events", []):
                if not isinstance(event, dict):
                    continue
                graph_id = int(event.get("graph_id", -1))
                decomposition = event.get("threshold_counterterms", {})
                if not isinstance(decomposition, dict):
                    continue
                for raw_component in decomposition.get("components", []):
                    if not isinstance(raw_component, dict):
                        continue
                    key = graph_id, int(raw_component.get("component_id", -1))
                    if key in self.components:
                        observed.add(key)
        return observed

    def single_facets(self, observed: set[tuple[int, int]]) -> list[ReportFacet]:
        groups: dict[tuple[int, int | None], dict[int, list[tuple[int, int]]]] = (
            defaultdict(lambda: defaultdict(list))
        )
        for key in observed:
            component = self.components[key]
            if len(component.variant_ids) != 1 or component.kind not in (
                "local",
                "integrated",
            ):
                continue
            groups[(component.graph_id, component.cut_group_id)][
                component.variant_ids[0]
            ].append(key)

        facets = []
        for (graph_id, cut_group_id), variants in sorted(
            groups.items(),
            key=lambda item: (
                item[0][0],
                item[0][1] is None,
                item[0][1] if item[0][1] is not None else -1,
            ),
        ):
            entries = []
            for variant_id, component_keys in sorted(variants.items()):
                keys_by_kind: dict[str, list[tuple[int, int]]] = defaultdict(list)
                for component_key in component_keys:
                    keys_by_kind[self.components[component_key].kind].append(
                        component_key
                    )
                for kind, kind_component_keys in sorted(
                    keys_by_kind.items(),
                    key=lambda item: THRESHOLD_KIND_ORDER.get(item[0], 99),
                ):
                    entries.append(
                        (
                            f"{self.variant(graph_id, variant_id).compact_label()} · {kind}",
                            tuple(sorted(kind_component_keys)),
                        )
                    )
            chunks = [entries[index : index + 8] for index in range(0, len(entries), 8)]
            for chunk_index, chunk in enumerate(chunks):
                title = self.cut_group_label(graph_id, cut_group_id)
                if len(chunks) > 1:
                    title += f" · {chunk_index + 1}/{len(chunks)}"
                facets.append(ReportFacet(title, OrderedDict(chunk)))
        return facets

    def pair_facets(self, observed: set[tuple[int, int]]) -> list[ReportFacet]:
        groups: dict[
            tuple[int, int | None, tuple[int, ...]],
            dict[str, list[tuple[int, int]]],
        ] = defaultdict(lambda: defaultdict(list))
        for key in observed:
            component = self.components[key]
            if len(component.variant_ids) != 2:
                continue
            groups[(component.graph_id, component.cut_group_id, component.variant_ids)][
                component.kind
            ].append(key)

        facets = []
        for (graph_id, cut_group_id, variant_ids), kinds in sorted(
            groups.items(),
            key=lambda item: (
                item[0][0],
                item[0][1] is None,
                item[0][1] if item[0][1] is not None else -1,
                item[0][2],
            ),
        ):
            variants = [
                self.variant(graph_id, variant_id).compact_label()
                for variant_id in variant_ids
            ]
            component_groups = OrderedDict(
                (
                    THRESHOLD_KIND_LABELS.get(kind, kind),
                    tuple(sorted(component_keys)),
                )
                for kind, component_keys in sorted(
                    kinds.items(),
                    key=lambda item: THRESHOLD_KIND_ORDER.get(item[0], 99),
                )
            )
            facets.append(
                ReportFacet(
                    "\n".join(
                        [
                            self.cut_group_label(graph_id, cut_group_id),
                            variants[0],
                            f"× {variants[1]}",
                        ]
                    ),
                    component_groups,
                )
            )
        return facets


@dataclass
class ReportFacet:
    title: str
    component_groups: OrderedDict[str, tuple[tuple[int, int], ...]]
    summary_labels: tuple[str, ...] = ()


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
        help="Plot |t| on a logarithmic x-axis",
    )
    parser.add_argument(
        "--branch-layout",
        choices=("auto", "overlay", "split"),
        default="auto",
        help=(
            "Layout for both logarithmic t branches: auto selects mirrored split "
            "panels, overlay draws both branches on one panel"
        ),
    )
    parser.add_argument(
        "--result-label",
        action="append",
        default=[],
        help=(
            "Legend label for the corresponding input result; repeat once per "
            "result file"
        ),
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
        "--threshold-report",
        action="store_true",
        help=(
            "Create a curated multipage threshold-counterterm closure, variant, "
            "and iterated-component report"
        ),
    )
    parser.add_argument(
        "--threshold-report-section",
        action="append",
        choices=("summary", "singles", "pairs"),
        default=[],
        help="Threshold report section to emit; can be repeated (default: all)",
    )
    parser.add_argument(
        "--facets-per-page",
        type=int,
        default=3,
        help="Maximum one-sided cut-group facets per threshold-report page (default: 3)",
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
    if args.branch_layout == "split" and not (
        args.x_log_scale and args.t_branch == "both"
    ):
        raise ValueError(
            "--branch-layout split requires --x-log-scale --t-branch both."
        )
    if args.fit_log_slope and not (args.x_log_scale and args.y_log_scale):
        raise ValueError("--fit-log-slope requires logarithmic x and y axes.")
    if args.result_label and len(args.result_label) != len(args.results):
        raise ValueError(
            "--result-label must be repeated exactly once per input result file "
            f"({len(args.results)} expected, {len(args.result_label)} supplied)."
        )
    if any(not label.strip() for label in args.result_label):
        raise ValueError("--result-label values must be nonblank.")
    if args.threshold_report and args.threshold_decomposition is not None:
        raise ValueError(
            "--threshold-report and --threshold-decomposition are mutually exclusive."
        )
    if args.threshold_report_section and not args.threshold_report:
        raise ValueError("--threshold-report-section requires --threshold-report.")
    if args.facets_per_page < 1:
        raise ValueError("--facets-per-page must be a positive integer.")
    args.split_branches = (
        args.x_log_scale
        and args.t_branch == "both"
        and args.branch_layout in ("auto", "split")
    )
    args.threshold_report_sections = list(
        dict.fromkeys(args.threshold_report_section or ("summary", "singles", "pairs"))
    )
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


def load_result(path: Path, label: str | None = None) -> ResultFile:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if "points" not in data:
        raise ValueError(
            f"{path} is not an approach result JSON file: missing 'points'"
        )
    return ResultFile(path=path, data=data, label=label)


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
    registry = ThresholdRegistry.from_result(result)
    return {
        key: registry.component_label(component)
        for key, component in registry.components.items()
    }


def threshold_decomposition_values(
    result: ResultFile,
    point: dict[str, Any],
    mode: str,
    warnings: NonFiniteWarnings,
) -> list[tuple[str, complex]]:
    catalog = threshold_component_catalog(result)
    found, original, components_by_id = threshold_point_values(
        result, point, warnings, catalog
    )
    if not found:
        return []

    components: dict[str, complex] = defaultdict(complex)
    for key, value in components_by_id.items():
        components[
            catalog.get(key, f"threshold:component c{key[1]} [graph#{key[0]}]")
        ] += value
    counterterms = sum(components_by_id.values(), 0j)
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


def threshold_point_values(
    result: ResultFile,
    point: dict[str, Any],
    warnings: NonFiniteWarnings,
    catalog: dict[tuple[int, int], str],
) -> tuple[bool, complex, dict[tuple[int, int], complex]]:
    evaluation = point.get("evaluation", {})
    events = evaluation.get("events", []) if isinstance(evaluation, dict) else []
    original = 0j
    components: dict[tuple[int, int], complex] = defaultdict(complex)
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
            key = graph_id, component_id
            label = catalog.get(
                key,
                f"threshold:component c{component_id} [graph#{graph_id}]",
            )
            weighted = complex_value(
                component.get("weighted"), warnings, result, point, label
            )
            components[key] += weighted
    return found, original, components


def selected_t_value(
    point: dict[str, Any], branch: str, x_log_scale: bool
) -> float | None:
    value = float(point["t"])
    if branch == "positive" and value <= 0.0:
        return None
    if branch == "negative" and value >= 0.0:
        return None
    return abs(value) if x_log_scale else value


def displayed_component(value: complex, component: str) -> float:
    if component == "real":
        return value.real
    if component == "imag":
        return value.imag
    return abs(value)


def collect_threshold_summary_series(
    result: ResultFile,
    axis_index: int,
    branch: str,
    args: argparse.Namespace,
    warnings: NonFiniteWarnings,
    selected_labels: tuple[str, ...],
) -> list[Series]:
    labels = ("event sum", "original O", "Σ CT", "reconstructed O+ΣCT", "closure Δ")
    values: OrderedDict[str, list[tuple[float, float]]] = OrderedDict(
        (label, []) for label in labels
    )
    catalog = threshold_component_catalog(result)
    for point in result.data.get("points", []):
        if (
            int(point.get("axis_index", -1)) != axis_index
            or point.get("status") != "evaluated"
        ):
            continue
        t_value = selected_t_value(point, branch, args.x_log_scale)
        if t_value is None:
            continue
        found, original, components = threshold_point_values(
            result, point, warnings, catalog
        )
        if not found:
            continue
        counterterms = sum(components.values(), 0j)
        reconstructed = original + counterterms
        event_sum = complex_value(
            point.get("evaluation", {}).get("event_weight_sum"),
            warnings,
            result,
            point,
            "event_weight_sum",
        )
        for label, value in zip(
            labels,
            (
                event_sum,
                original,
                counterterms,
                reconstructed,
                event_sum - reconstructed,
            ),
        ):
            values[label].append((t_value, displayed_component(value, args.component)))
    return [
        Series(
            label,
            [t_value for t_value, _ in samples],
            [value for _, value in samples],
        )
        for label, samples in values.items()
        if samples and (not selected_labels or label in selected_labels)
    ]


def collect_threshold_facet_series(
    result: ResultFile,
    axis_index: int,
    branch: str,
    facet: ReportFacet,
    component: str,
    args: argparse.Namespace,
    warnings: NonFiniteWarnings,
) -> list[Series]:
    values: OrderedDict[str, list[tuple[float, float]]] = OrderedDict(
        (label, []) for label in facet.component_groups
    )
    catalog = threshold_component_catalog(result)
    for point in result.data.get("points", []):
        if (
            int(point.get("axis_index", -1)) != axis_index
            or point.get("status") != "evaluated"
        ):
            continue
        t_value = selected_t_value(point, branch, args.x_log_scale)
        if t_value is None:
            continue
        _, _, components = threshold_point_values(result, point, warnings, catalog)
        for label, component_keys in facet.component_groups.items():
            value = sum((components.get(key, 0j) for key in component_keys), 0j)
            values[label].append((t_value, displayed_component(value, component)))
    return [
        Series(
            label,
            [t_value for t_value, _ in samples],
            [value for _, value in samples],
        )
        for label, samples in values.items()
        if samples
    ]


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
) -> bool:
    if not series.t_values:
        return False
    label = series.label
    if fit_log_slope and (slope := log_slope(series, x_range)) is not None:
        label = f"{label} (p={slope:+.2f})"
    if not y_log_scale or component == "abs":
        y_values = [
            finite_abs(value) if y_log_scale else value for value in series.values
        ]
        ax.plot(series.t_values, y_values, linewidth=1.4, label=label)
        return False

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
            label="_nolegend_" if positive_line is not None else label,
        )
        return True
    return False


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


def decorate_axis(
    ax: plt.Axes,
    result: ResultFile,
    component: str,
    y_log_scale: bool,
    x_log_scale: bool,
    hide_info_box: bool,
    branch: str | None = None,
    show_ylabel: bool = True,
    show_xlabel: bool = True,
) -> None:
    if show_xlabel:
        if x_log_scale and branch in ("negative", "positive"):
            sign = "< 0" if branch == "negative" else "> 0"
            ax.set_xlabel(f"|approach parameter t|  (t {sign})")
        else:
            ax.set_xlabel(
                "|approach parameter t|" if x_log_scale else "approach parameter t"
            )
    if x_log_scale:
        ax.set_xscale("log")
    ylabel = "weight magnitude" if component == "abs" else f"{component} weight"
    if y_log_scale:
        if component != "abs":
            ylabel = f"|{ylabel}|"
        ax.set_yscale("log")
    if show_ylabel:
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


def apply_axis_ranges(
    ax: plt.Axes, args: argparse.Namespace, *, reverse_x: bool = False
) -> None:
    if args.x_range is not None:
        limits = args.x_range[::-1] if reverse_x else args.x_range
        ax.set_xlim(limits[0], limits[1])
    elif reverse_x:
        ax.invert_xaxis()
    if args.y_range is not None:
        ax.set_ylim(args.y_range[0], args.y_range[1])


def branches_on_axis(args: argparse.Namespace) -> list[tuple[str, str]]:
    if args.t_branch != "both":
        return [(args.t_branch, "")]
    if args.x_log_scale:
        return [("negative", "t<0"), ("positive", "t>0")]
    return [("both", "")]


def shared_figure_legend(
    fig: plt.Figure, axes: Iterable[plt.Axes], has_negative_values: bool
) -> float:
    entries: OrderedDict[str, Any] = OrderedDict()
    for ax in axes:
        handles, labels = ax.get_legend_handles_labels()
        for handle, label in zip(handles, labels):
            if label and label != "_nolegend_":
                entries.setdefault(label, handle)
    if has_negative_values:
        entries["dashed = negative value"] = Line2D(
            [0], [0], color="0.25", linewidth=1.4, linestyle="--"
        )
    if not entries:
        return 0.08
    columns = min(5, max(1, len(entries)))
    rows = math.ceil(len(entries) / columns)
    fig.legend(
        list(entries.values()),
        list(entries),
        loc="lower center",
        bbox_to_anchor=(0.5, 0.018),
        fontsize=6.2 if len(entries) > 8 else 7.0,
        frameon=True,
        ncols=columns,
    )
    return min(0.34, 0.1 + 0.045 * rows)


def result_label_prefix(
    result: ResultFile,
    axis_index: int,
    *,
    multi_result: bool,
    multi_axis: bool,
) -> str:
    parts = []
    if result.label is not None or multi_result:
        parts.append(result.label or result.path.stem)
    if multi_axis:
        parts.append(f"axis {axis_index + 1}")
    return (" / ".join(parts) + " / ") if parts else ""


def draw_page(
    pdf: PdfPages,
    result_axis_pairs: list[tuple[ResultFile, int]],
    args: argparse.Namespace,
    include_patterns: list[re.Pattern[str]],
    exclude_patterns: list[re.Pattern[str]],
    warnings: NonFiniteWarnings,
) -> None:
    column_count = 2 if args.split_branches else 1
    fig, raw_axes = plt.subplots(
        1,
        column_count,
        figsize=(13.5, 8.2),
        sharey=column_count == 2,
        squeeze=False,
    )
    axes = list(raw_axes[0])
    has_series = False
    has_negative_values = False
    multi_result = len({pair[0].path for pair in result_axis_pairs}) > 1
    multi_axis = len({pair[1] for pair in result_axis_pairs}) > 1

    for column, ax in enumerate(axes):
        branch_specs = (
            ([("negative", "")] if column == 0 else [("positive", "")])
            if args.split_branches
            else branches_on_axis(args)
        )
        axis_has_series = False
        for branch, branch_suffix in branch_specs:
            for result, axis_index in result_axis_pairs:
                label_prefix = result_label_prefix(
                    result,
                    axis_index,
                    multi_result=multi_result or args.combine_plots,
                    multi_axis=multi_axis or args.combine_axes,
                )
                series = collect_series(
                    result,
                    axis_index,
                    args.component,
                    include_patterns,
                    exclude_patterns,
                    args.sum_lmb_samples_per_cut,
                    warnings,
                    branch,
                    args.x_log_scale,
                    args.threshold_decomposition,
                    label_prefix=label_prefix,
                )
                for item in series:
                    if branch_suffix:
                        item.label = f"{item.label} · {branch_suffix}"
                    has_series = axis_has_series = True
                    has_negative_values |= plot_series(
                        ax,
                        item,
                        args.component,
                        args.y_log_scale,
                        args.fit_log_slope,
                        args.x_range,
                    )
        if not axis_has_series:
            ax.text(
                0.5,
                0.5,
                "No evaluated series matched the requested filters.",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=11,
                color="0.35",
            )

    first_result, first_axis = result_axis_pairs[0]
    title = args.title or axis_label(first_result, first_axis)
    if args.title is None and len(result_axis_pairs) > 1:
        title = "combined approach curves"
    fig.suptitle(title, fontsize=12, fontweight="bold")
    for column, ax in enumerate(axes):
        branch = ("negative", "positive")[column] if args.split_branches else None
        if args.split_branches:
            ax.set_title("t < 0" if branch == "negative" else "t > 0", fontsize=9)
        decorate_axis(
            ax,
            first_result,
            args.component,
            args.y_log_scale,
            args.x_log_scale,
            args.hide_info_box or column > 0,
            branch=branch,
            show_ylabel=column == 0,
        )
        apply_axis_ranges(ax, args, reverse_x=branch == "negative")
    bottom = (
        shared_figure_legend(fig, axes, has_negative_values) if has_series else 0.08
    )
    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        top=0.9,
        bottom=bottom,
        wspace=0.06 if args.split_branches else 0.2,
    )
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def draw_threshold_report_page(
    pdf: PdfPages,
    result: ResultFile,
    axis_index: int,
    facets: list[ReportFacet],
    args: argparse.Namespace,
    warnings: NonFiniteWarnings,
    *,
    title: str,
    component: str,
    summary: bool = False,
) -> None:
    row_count = len(facets)
    column_count = 2 if args.split_branches else 1
    fig, axes = plt.subplots(
        row_count,
        column_count,
        figsize=(13.5, 8.2),
        sharey="row" if column_count == 2 else False,
        squeeze=False,
    )
    has_negative_values = False
    flat_axes = [ax for row in axes for ax in row]
    for row, facet in enumerate(facets):
        for column in range(column_count):
            ax = axes[row][column]
            branch_specs = (
                ([("negative", "")] if column == 0 else [("positive", "")])
                if args.split_branches
                else branches_on_axis(args)
            )
            has_series = False
            for branch, branch_suffix in branch_specs:
                series = (
                    collect_threshold_summary_series(
                        result,
                        axis_index,
                        branch,
                        args,
                        warnings,
                        facet.summary_labels,
                    )
                    if summary
                    else collect_threshold_facet_series(
                        result,
                        axis_index,
                        branch,
                        facet,
                        component,
                        args,
                        warnings,
                    )
                )
                for item in series:
                    if branch_suffix:
                        item.label = f"{item.label} · {branch_suffix}"
                    has_series = True
                    has_negative_values |= plot_series(
                        ax,
                        item,
                        component,
                        args.y_log_scale,
                        args.fit_log_slope,
                        args.x_range,
                    )
            if not has_series:
                ax.text(
                    0.5,
                    0.5,
                    "No observed components.",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    color="0.35",
                )
            branch = ("negative", "positive")[column] if args.split_branches else None
            axis_title = facet.title
            if args.split_branches:
                axis_title += "\n" + ("t < 0" if branch == "negative" else "t > 0")
            ax.set_title(
                axis_title,
                fontsize=8.2,
                fontweight="semibold",
                linespacing=1.12,
            )
            decorate_axis(
                ax,
                result,
                component,
                args.y_log_scale,
                args.x_log_scale,
                args.hide_info_box or row > 0 or column > 0,
                branch=branch,
                show_ylabel=column == 0,
                show_xlabel=row + 1 == row_count,
            )
            apply_axis_ranges(ax, args, reverse_x=branch == "negative")

    fig.suptitle(title, fontsize=12, fontweight="bold")
    bottom = shared_figure_legend(fig, flat_axes, has_negative_values)
    panel_title_lines = max(facet.title.count("\n") + 1 for facet in facets)
    if args.split_branches:
        panel_title_lines += 1
    extra_title_lines = max(0, panel_title_lines - 2)
    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        top=0.88 - 0.018 * extra_title_lines,
        bottom=bottom,
        hspace=0.52 + 0.12 * extra_title_lines,
        wspace=0.06 if args.split_branches else 0.2,
    )
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def threshold_report_title(
    result: ResultFile, axis_index: int, args: argparse.Namespace, section: str
) -> str:
    if args.title:
        return f"{args.title} · {section}"
    result_label = result.label or result.path.stem
    return f"{result_label} · axis {axis_index + 1} · {section}"


def draw_threshold_report(
    pdf: PdfPages,
    results: list[ResultFile],
    args: argparse.Namespace,
    warnings: NonFiniteWarnings,
) -> int:
    page_count = 0
    for result in results:
        registry = ThresholdRegistry.from_result(result)
        for axis_index in axis_indices(result, args.selected_axis_indices):
            observed = registry.observed_components(result, axis_index)
            if not observed:
                continue
            if "summary" in args.threshold_report_sections:
                summary_facets = [
                    ReportFacet(
                        "event decomposition",
                        OrderedDict(),
                        ("event sum", "original O", "Σ CT", "reconstructed O+ΣCT"),
                    ),
                    ReportFacet(
                        "closure residual",
                        OrderedDict(),
                        ("closure Δ",),
                    ),
                ]
                draw_threshold_report_page(
                    pdf,
                    result,
                    axis_index,
                    summary_facets,
                    args,
                    warnings,
                    title=threshold_report_title(result, axis_index, args, "summary"),
                    component=args.component,
                    summary=True,
                )
                page_count += 1

            if "singles" in args.threshold_report_sections:
                single_facets = registry.single_facets(observed)
                chunks = [
                    single_facets[index : index + args.facets_per_page]
                    for index in range(0, len(single_facets), args.facets_per_page)
                ]
                for chunk_index, facets in enumerate(chunks):
                    suffix = (
                        f" {chunk_index + 1}/{len(chunks)}" if len(chunks) > 1 else ""
                    )
                    draw_threshold_report_page(
                        pdf,
                        result,
                        axis_index,
                        facets,
                        args,
                        warnings,
                        title=threshold_report_title(
                            result,
                            axis_index,
                            args,
                            f"variants{suffix}",
                        ),
                        component=args.component,
                    )
                    page_count += 1

            if "pairs" in args.threshold_report_sections:
                pair_facets = registry.pair_facets(observed)
                chunks = [
                    pair_facets[index : index + 2]
                    for index in range(0, len(pair_facets), 2)
                ]
                for chunk_index, facets in enumerate(chunks):
                    suffix = (
                        f" {chunk_index + 1}/{len(chunks)}" if len(chunks) > 1 else ""
                    )
                    draw_threshold_report_page(
                        pdf,
                        result,
                        axis_index,
                        facets,
                        args,
                        warnings,
                        title=threshold_report_title(
                            result,
                            axis_index,
                            args,
                            f"pairs{suffix}",
                        ),
                        component=args.component,
                    )
                    page_count += 1
    return page_count


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
        labels = args.result_label or [None] * len(args.results)
        results = [
            load_result(path, label) for path, label in zip(args.results, labels)
        ]
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
            if args.threshold_report:
                if draw_threshold_report(pdf, results, args, warnings) == 0:
                    raise ValueError(
                        "No observed threshold-counterterm report sections were found."
                    )
            else:
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
