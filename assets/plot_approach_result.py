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


def compact_edge_sets(
    edge_sets: Iterable[tuple[int, ...]], head: str, *, limit: int = 2
) -> str:
    unique = sorted(set(edge_sets))
    labels = [f"{head}({','.join(map(str, edges))})" for edges in unique[:limit]]
    if len(unique) > limit:
        labels.append(f"…(+{len(unique) - limit})")
    return "/".join(labels) or f"{head}(?)"


@dataclass
class ResultFile:
    path: Path
    data: dict[str, Any]
    label: str | None = None
    supplemental_series: list[SupplementalSeries] | None = None


@dataclass
class Series:
    label: str
    t_values: list[float]
    values: list[float]


@dataclass(frozen=True)
class SupplementalSeries:
    label: str
    axis_index: int
    values: dict[int, Any]


@dataclass(frozen=True)
class LogFit:
    slope: float
    r_squared: float
    point_count: int


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
        side = {"left": "left", "right": "right", "amplitude": "amp"}.get(
            self.side, self.side or "?"
        )
        subspace = ",".join(map(str, self.subspace)) or "max"
        return f"v{self.variant_id} {self.name} · {side} {eta} · S[{subspace}]"

    def panel_label(self) -> str:
        thresholds = (threshold for _, threshold in self.associations)
        eta = compact_edge_sets(thresholds, "η")
        side = {"left": "L", "right": "R", "amplitude": "A"}.get(
            self.side, self.side[:1].upper() or "?"
        )
        subspace = ",".join(map(str, self.subspace)) or "max"
        return f"v{self.variant_id}:{self.name} {side} {eta} S[{subspace}]"


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
            summary = evaluation.get("threshold_counterterm_summary", {})
            if isinstance(summary, dict):
                for raw_component in summary.get("components", []):
                    if not isinstance(raw_component, dict):
                        continue
                    key = (
                        int(raw_component.get("graph_id", -1)),
                        int(raw_component.get("component_id", -1)),
                    )
                    if key in self.components:
                        observed.add(key)
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
            for variant_id, component_keys in sorted(variants.items()):
                keys_by_kind: dict[str, list[tuple[int, int]]] = defaultdict(list)
                for component_key in component_keys:
                    keys_by_kind[self.components[component_key].kind].append(
                        component_key
                    )
                component_groups = OrderedDict(
                    (
                        THRESHOLD_KIND_LABELS.get(kind, kind),
                        tuple(sorted(kind_component_keys)),
                    )
                    for kind, kind_component_keys in sorted(
                        keys_by_kind.items(),
                        key=lambda item: THRESHOLD_KIND_ORDER.get(item[0], 99),
                    )
                )
                facets.append(
                    ReportFacet(
                        "\n".join(
                            [
                                self.cut_group_label(graph_id, cut_group_id),
                                self.variant(graph_id, variant_id).compact_label(),
                            ]
                        ),
                        component_groups,
                    )
                )
        return facets

    def multiplier_facets(self, observed: set[tuple[int, int]]) -> list[ReportFacet]:
        facets = []
        for key in sorted(observed):
            component = self.components[key]
            variants = [
                self.variant(component.graph_id, variant_id)
                for variant_id in component.variant_ids
            ]
            graph_name = (
                variants[0].graph_name if variants else f"#{component.graph_id}"
            )
            cut_group = (
                "cg–"
                if component.cut_group_id is None
                else f"cg{component.cut_group_id}"
            )
            kind = THRESHOLD_KIND_LABELS.get(component.kind, component.kind)
            facets.append(
                ReportFacet(
                    "\n".join(
                        [
                            f"{graph_name} · {cut_group} · c{component.component_id} {kind}",
                            " × ".join(variant.panel_label() for variant in variants),
                        ]
                    ),
                    OrderedDict((("multipliers", (key,)),)),
                )
            )
        return facets

    def filter_components(
        self,
        observed: set[tuple[int, int]],
        include_patterns: list[re.Pattern[str]],
        exclude_patterns: list[re.Pattern[str]],
    ) -> set[tuple[int, int]]:
        return {
            key
            for key in observed
            if label_is_selected(
                self.component_label(self.components[key]),
                include_patterns,
                exclude_patterns,
            )
        }

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
    y_label: str | None = None


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

    def emit(self, limit: int) -> None:
        items = list(self._items.values())
        for item in items[:limit]:
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
                "  Rendered non-finite/null value(s) as gaps instead of numerical zeros.",
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
        omitted = len(items) - limit
        if omitted > 0:
            print(
                "plot_approach_result.py: warning: "
                f"{omitted} additional gap sample(s) omitted; "
                "raise --max-gap-warnings to inspect more.",
                file=sys.stderr,
            )


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
        help=(
            "Use logarithmic |t| scaling; with both branches, the automatic "
            "layout uses one signed symmetric-log axis"
        ),
    )
    parser.add_argument(
        "--branch-layout",
        choices=("auto", "signed", "overlay", "split"),
        default="auto",
        help=(
            "Layout for both logarithmic t branches: auto selects one signed "
            "symmetric-log axis, overlay draws both against |t|, and split uses "
            "mirrored panels"
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
        help=(
            "Append a log-log power fit over the smallest visible approach "
            "parameters on each signed branch"
        ),
    )
    parser.add_argument(
        "--fit-points",
        type=int,
        help="Number of smallest visible points to use per log-log fit (default: 4)",
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
        choices=("summary", "singles", "pairs", "multipliers"),
        default=[],
        help=(
            "Threshold report section to emit; can be repeated "
            "(default: summary, singles, and pairs)"
        ),
    )
    parser.add_argument(
        "--threshold-quantity",
        action="append",
        choices=("weighted", "bare"),
        default=[],
        help=(
            "Counterterm quantity for variant and pair report pages; repeat to "
            "show both (default: weighted)"
        ),
    )
    parser.add_argument(
        "--multiplier-log-scale",
        action="store_true",
        help="Use a logarithmic y-axis for the multiplier report section",
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
        "--y-label",
        help="Replace the generated y-axis label",
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
        help=(
            "Visible x-axis range to show; signed symmetric-log plots use "
            "the signed limits directly"
        ),
    )
    parser.add_argument(
        "--x-symlog-linthresh",
        type=float,
        help=(
            "Width of the linear region around zero for a signed symmetric-log "
            "x-axis (default: spacing.min_abs_t)"
        ),
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
    parser.add_argument(
        "--series-json",
        action="append",
        default=[],
        type=Path,
        help=(
            "Supplemental scalar-series JSON keyed by approach point index; can "
            "be repeated and currently requires exactly one approach result"
        ),
    )
    parser.add_argument(
        "--max-gap-warnings",
        type=int,
        default=5,
        help=(
            "Maximum detailed non-finite/skipped-sample gap warnings to print "
            "(default: 5; 0 prints only the omitted count)"
        ),
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
    signed_branches = (
        args.x_log_scale
        and args.t_branch == "both"
        and args.branch_layout in ("auto", "signed")
    )
    requested_report_sections = set(
        args.threshold_report_section or ("summary", "singles", "pairs")
    )
    report_uses_weight_scale = bool(
        requested_report_sections.intersection(("summary", "singles", "pairs"))
    )
    report_uses_multiplier_scale = "multipliers" in requested_report_sections
    y_range_uses_log_scale = (
        args.y_log_scale
        if not args.threshold_report
        else (
            (report_uses_weight_scale and args.y_log_scale)
            or (report_uses_multiplier_scale and args.multiplier_log_scale)
        )
    )
    validate_range_option(
        "--x-range", args.x_range, log_scale=args.x_log_scale and not signed_branches
    )
    if (
        signed_branches
        and args.x_range is not None
        and not (args.x_range[0] < 0.0 < args.x_range[1])
    ):
        raise ValueError(
            "A signed symmetric-log --x-range must straddle zero; use overlay or "
            "split layout for a positive |t| range."
        )
    validate_range_option("--y-range", args.y_range, log_scale=y_range_uses_log_scale)
    if args.branch_layout in ("signed", "split") and not (
        args.x_log_scale and args.t_branch == "both"
    ):
        raise ValueError(
            f"--branch-layout {args.branch_layout} requires "
            "--x-log-scale --t-branch both."
        )
    if args.x_symlog_linthresh is not None:
        if not math.isfinite(args.x_symlog_linthresh) or args.x_symlog_linthresh <= 0:
            raise ValueError("--x-symlog-linthresh must be finite and positive.")
        if not signed_branches:
            raise ValueError(
                "--x-symlog-linthresh requires a signed symmetric-log layout."
            )
    fit_has_log_y_axis = (
        args.y_log_scale
        if not args.threshold_report
        else (
            (report_uses_weight_scale and args.y_log_scale)
            or (report_uses_multiplier_scale and args.multiplier_log_scale)
        )
    )
    if args.fit_log_slope and not (args.x_log_scale and fit_has_log_y_axis):
        raise ValueError("--fit-log-slope requires logarithmic x and y axes.")
    if args.fit_points is not None and not args.fit_log_slope:
        raise ValueError("--fit-points requires --fit-log-slope.")
    if args.fit_points is not None and args.fit_points < 2:
        raise ValueError("--fit-points must be at least 2.")
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
    if args.threshold_quantity and not args.threshold_report:
        raise ValueError("--threshold-quantity requires --threshold-report.")
    if args.multiplier_log_scale and not args.threshold_report:
        raise ValueError("--multiplier-log-scale requires --threshold-report.")
    if args.facets_per_page < 1:
        raise ValueError("--facets-per-page must be a positive integer.")
    if args.max_gap_warnings < 0:
        raise ValueError("--max-gap-warnings must be non-negative.")
    if args.series_json and len(args.results) != 1:
        raise ValueError(
            "--series-json currently requires exactly one approach result."
        )
    if args.series_json and args.threshold_report:
        raise ValueError("--series-json is only available for ordinary approach plots.")
    args.split_branches = (
        args.x_log_scale and args.t_branch == "both" and args.branch_layout == "split"
    )
    args.signed_branches = signed_branches
    args.threshold_report_sections = list(
        dict.fromkeys(args.threshold_report_section or ("summary", "singles", "pairs"))
    )
    args.threshold_quantities = list(
        dict.fromkeys(args.threshold_quantity or ("weighted",))
    )
    args.fit_points = args.fit_points or 4
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
    if data.get("schema_version") != 3:
        raise ValueError(
            f"{path} uses approach schema_version {data.get('schema_version')!r}; "
            "plot_approach_result.py requires schema_version 3"
        )
    parameter_name = data.get("parameter_name")
    if not isinstance(parameter_name, str) or not parameter_name.strip():
        raise ValueError(f"{path} has an invalid or missing 'parameter_name'.")
    axes = data.get("axes")
    if not isinstance(axes, list) or not axes:
        raise ValueError(f"{path} has an invalid or empty 'axes' list.")
    for axis_index, axis in enumerate(axes):
        if not isinstance(axis, dict):
            raise ValueError(f"{path} axis {axis_index + 1} is not an object.")
        axis_name = axis.get("label")
        vector = axis.get("vector")
        if not isinstance(axis_name, str) or not axis_name.strip():
            raise ValueError(f"{path} axis {axis_index + 1} has an invalid label.")
        if (
            not isinstance(vector, list)
            or not vector
            or any(
                not isinstance(component, (int, float))
                or not math.isfinite(float(component))
                for component in vector
            )
        ):
            raise ValueError(f"{path} axis {axis_index + 1} has an invalid vector.")
    spacing = data.get("spacing")
    if not isinstance(spacing, dict) or "max_abs_t" not in spacing:
        raise ValueError(f"{path} has no schema-v3 spacing.max_abs_t value.")
    try:
        max_abs_t = float(spacing["max_abs_t"])
    except (TypeError, ValueError) as error:
        raise ValueError(f"{path} has an invalid spacing.max_abs_t value.") from error
    if not math.isfinite(max_abs_t) or max_abs_t <= 0.0:
        raise ValueError(f"{path} has an invalid spacing.max_abs_t value.")
    if not isinstance(data.get("points"), list):
        raise ValueError(
            f"{path} is not an approach result JSON file: missing 'points'."
        )
    point_indices = set()
    for point_index, point in enumerate(data["points"]):
        if not isinstance(point, dict):
            raise ValueError(f"{path} point {point_index + 1} is not an object.")
        stored_point_index = point.get("index")
        if (
            not isinstance(stored_point_index, int)
            or stored_point_index in point_indices
        ):
            raise ValueError(
                f"{path} point {point_index + 1} has an invalid or duplicate index."
            )
        point_indices.add(stored_point_index)
        axis_index = point.get("axis_index")
        if not isinstance(axis_index, int) or not 0 <= axis_index < len(axes):
            raise ValueError(
                f"{path} point {point_index + 1} has an invalid axis_index."
            )
        if point.get("status") != "evaluated":
            continue
        evaluation = point.get("evaluation")
        if not isinstance(evaluation, dict):
            raise ValueError(
                f"{path} evaluated point {point_index + 1} has no evaluation."
            )
        if not isinstance(evaluation.get("additional_contribution_sums"), dict):
            raise ValueError(
                f"{path} evaluated point {point_index + 1} has no "
                "additional_contribution_sums map."
            )
        events = evaluation.get("events")
        if not isinstance(events, list):
            raise ValueError(
                f"{path} evaluated point {point_index + 1} has no events list."
            )
        for event_index, event in enumerate(events):
            if not isinstance(event, dict) or not isinstance(
                event.get("display_only_weights"), dict
            ):
                raise ValueError(
                    f"{path} point {point_index + 1} event {event_index + 1} has no "
                    "schema-v3 display_only_weights map."
                )
        summary = evaluation.get("threshold_counterterm_summary")
        if summary is None:
            continue
        if not isinstance(summary, dict) or not all(
            field in summary
            for field in (
                "original",
                "counterterm_sum",
                "reconstructed_total",
                "decomposed_event_sum",
                "closure",
                "components",
            )
        ):
            raise ValueError(
                f"{path} point {point_index + 1} has an invalid "
                "threshold_counterterm_summary."
            )
    return ResultFile(path=path, data=data, label=label, supplemental_series=[])


def load_supplemental_series(paths: list[Path], result: ResultFile) -> None:
    known_points = {
        int(point["index"]): point for point in result.data.get("points", [])
    }
    known_axes = set(axis_indices(result))
    seen = set()
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        if data.get("schema_version") != 1 or not isinstance(data.get("series"), list):
            raise ValueError(
                f"{path} must use supplemental-series schema_version 1 with a 'series' list."
            )
        for series_index, raw_series in enumerate(data["series"]):
            if not isinstance(raw_series, dict):
                raise ValueError(f"{path} series {series_index + 1} is not an object.")
            label = raw_series.get("label")
            axis_index = raw_series.get("axis_index")
            raw_values = raw_series.get("values")
            if not isinstance(label, str) or not label.strip():
                raise ValueError(
                    f"{path} series {series_index + 1} has an invalid label."
                )
            if not isinstance(axis_index, int) or axis_index not in known_axes:
                raise ValueError(
                    f"{path} series {series_index + 1} has unknown axis_index {axis_index!r}."
                )
            if not isinstance(raw_values, list):
                raise ValueError(
                    f"{path} series {series_index + 1} has no values list."
                )
            key = axis_index, label
            if key in seen:
                raise ValueError(
                    f"duplicate supplemental series label {label!r} for axis {axis_index + 1}."
                )
            seen.add(key)
            values = {}
            for value_index, raw_value in enumerate(raw_values):
                if not isinstance(raw_value, dict):
                    raise ValueError(
                        f"{path} series {series_index + 1} value {value_index + 1} is not an object."
                    )
                point_index = raw_value.get("point_index")
                if not isinstance(point_index, int) or point_index not in known_points:
                    raise ValueError(
                        f"{path} series {series_index + 1} references unknown point_index "
                        f"{point_index!r}."
                    )
                if int(known_points[point_index].get("axis_index", -1)) != axis_index:
                    raise ValueError(
                        f"{path} series {series_index + 1} point {point_index} belongs to "
                        f"a different approach axis."
                    )
                if point_index in values:
                    raise ValueError(
                        f"{path} series {series_index + 1} repeats point_index {point_index}."
                    )
                values[point_index] = raw_value.get("value")
            result.supplemental_series.append(
                SupplementalSeries(f"series:{label}", axis_index, values)
            )


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
        axis = axes[axis_index]
        return str(axis["label"])
    return f"axis {axis_index + 1}"


def axis_vector(result: ResultFile, axis_index: int) -> list[Any]:
    axes = result.data.get("axes", [])
    if 0 <= axis_index < len(axes):
        return axes[axis_index]["vector"]
    return []


def parameter_name(result: ResultFile) -> str:
    return str(result.data["parameter_name"])


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
        return math.nan
    if not math.isfinite(number):
        warnings.add(result, point, label, component, value)
        return math.nan
    return number


def sanitized_complex_component(
    value: dict[str, Any],
    component: str,
    warnings: NonFiniteWarnings,
    result: ResultFile,
    point: dict[str, Any],
    label: str,
) -> float:
    re_value = sanitized_float(value.get("re"), warnings, result, point, label, "re")
    im_value = sanitized_float(value.get("im"), warnings, result, point, label, "im")
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
        warnings.add(result, point, label, "value", value)
        return complex(math.nan, math.nan)
    return complex(
        sanitized_float(value.get("re"), warnings, result, point, label, "re"),
        sanitized_float(value.get("im"), warnings, result, point, label, "im"),
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
    summary = threshold_summary(point)
    if summary is None:
        return []
    original = complex_value(
        summary.get("original"), warnings, result, point, "threshold:original"
    )
    counterterms = complex_value(
        summary.get("counterterm_sum"),
        warnings,
        result,
        point,
        "threshold:counterterm_sum",
    )
    reconstructed = complex_value(
        summary.get("reconstructed_total"),
        warnings,
        result,
        point,
        "threshold:decomposition_total",
    )
    components: dict[str, complex] = defaultdict(complex)
    for raw_component in summary.get("components", []):
        if not isinstance(raw_component, dict):
            continue
        key = (
            int(raw_component.get("graph_id", -1)),
            int(raw_component.get("component_id", -1)),
        )
        components[
            catalog.get(key, f"threshold:component c{key[1]} [graph#{key[0]}]")
        ] += complex_value(
            raw_component.get("weighted"),
            warnings,
            result,
            point,
            catalog.get(key, f"threshold:component c{key[1]} [graph#{key[0]}]"),
        )
    values: list[tuple[str, complex]] = []
    if mode in ("summary", "all"):
        values.extend(
            [
                ("threshold:original", original),
                ("threshold:counterterm_sum", counterterms),
                ("threshold:decomposition_total", reconstructed),
            ]
        )
    if mode in ("components", "all"):
        values.extend(sorted(components.items()))
    return values


def threshold_summary(point: dict[str, Any]) -> dict[str, Any] | None:
    evaluation = point.get("evaluation")
    if not isinstance(evaluation, dict):
        return None
    summary = evaluation.get("threshold_counterterm_summary")
    return summary if isinstance(summary, dict) else None


def axis_has_threshold_summary(result: ResultFile, axis_index: int) -> bool:
    return any(
        int(point.get("axis_index", -1)) == axis_index
        and point.get("status") == "evaluated"
        and threshold_summary(point) is not None
        for point in result.data.get("points", [])
    )


def threshold_component_values(
    result: ResultFile,
    point: dict[str, Any],
    quantity: str,
    warnings: NonFiniteWarnings,
    catalog: dict[tuple[int, int], str],
    selected_keys: set[tuple[int, int]],
) -> dict[tuple[int, int], complex]:
    if quantity == "weighted":
        summary = threshold_summary(point)
        if summary is None:
            return {}
        components: dict[tuple[int, int], complex] = defaultdict(complex)
        for raw_component in summary.get("components", []):
            if not isinstance(raw_component, dict):
                continue
            key = (
                int(raw_component.get("graph_id", -1)),
                int(raw_component.get("component_id", -1)),
            )
            if key not in selected_keys:
                continue
            components[key] += complex_value(
                raw_component.get("weighted"),
                warnings,
                result,
                point,
                catalog.get(
                    key,
                    f"threshold:component c{key[1]} [graph#{key[0]}]",
                ),
            )
        return components

    evaluation = point.get("evaluation", {})
    events = evaluation.get("events", []) if isinstance(evaluation, dict) else []
    components: dict[tuple[int, int], complex] = defaultdict(complex)
    missing_bare: set[tuple[int, int]] = set()
    observed_counts: dict[tuple[int, int], int] = defaultdict(int)
    for event in events:
        if not isinstance(event, dict):
            continue
        decomposition = event.get("threshold_counterterms")
        if not isinstance(decomposition, dict):
            continue
        graph_id = int(event.get("graph_id", -1))
        for component in decomposition.get("components", []):
            if not isinstance(component, dict):
                continue
            component_id = int(component.get("component_id", -1))
            key = graph_id, component_id
            if key not in selected_keys:
                continue
            observed_counts[key] += 1
            label = catalog.get(
                key,
                f"threshold:component c{component_id} [graph#{graph_id}]",
            )
            bare = component.get("bare")
            if component.get("evaluation_skipped") is True or not isinstance(
                bare, dict
            ):
                missing_bare.add(key)
                warnings.add(
                    result,
                    point,
                    label,
                    "bare",
                    "missing because the component evaluation was skipped",
                )
                continue
            components[key] += complex_value(bare, warnings, result, point, label)
    summary = threshold_summary(point)
    if summary is not None:
        for raw_component in summary.get("components", []):
            if not isinstance(raw_component, dict):
                continue
            key = (
                int(raw_component.get("graph_id", -1)),
                int(raw_component.get("component_id", -1)),
            )
            if key not in selected_keys:
                continue
            expected_count = int(raw_component.get("occurrence_count", 0))
            skipped_count = int(raw_component.get("skipped_count", 0))
            if skipped_count == 0 and observed_counts[key] >= expected_count:
                continue
            missing_bare.add(key)
            label = catalog.get(
                key,
                f"threshold:component c{key[1]} [graph#{key[0]}]",
            )
            warnings.add(
                result,
                point,
                label,
                "bare",
                (
                    "missing because the component evaluation was skipped"
                    if skipped_count
                    else "missing detailed runtime occurrence"
                ),
            )
    for key in missing_bare:
        components[key] = complex(math.nan, math.nan)
    return components


def selected_t_value(
    point: dict[str, Any], branch: str, x_log_scale: bool, signed_x: bool = False
) -> float | None:
    value = float(point["t"])
    if branch == "positive" and value <= 0.0:
        return None
    if branch == "negative" and value >= 0.0:
        return None
    return abs(value) if x_log_scale and not signed_x else value


def selected_evaluated_points(
    result: ResultFile,
    axis_index: int,
    branch: str,
    x_log_scale: bool,
    signed_x: bool = False,
) -> list[tuple[dict[str, Any], float]]:
    points = []
    for point in result.data.get("points", []):
        if (
            int(point.get("axis_index", -1)) != axis_index
            or point.get("status") != "evaluated"
        ):
            continue
        t_value = selected_t_value(point, branch, x_log_scale, signed_x)
        if t_value is not None:
            points.append((point, t_value))
    points.sort(key=lambda item: item[1])
    return points


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
    labels = (
        "decomposed event sum",
        "original O",
        "Σ weighted CT",
        "reconstructed O+ΣCT",
        "relative closure |Δ|/scale",
    )
    values: OrderedDict[str, list[tuple[float, float]]] = OrderedDict(
        (label, []) for label in labels
    )
    for point, t_value in selected_evaluated_points(
        result, axis_index, branch, args.x_log_scale, args.signed_branches
    ):
        summary = threshold_summary(point)
        if summary is None:
            for label in labels:
                values[label].append((t_value, math.nan))
            continue
        event_sum = complex_value(
            summary.get("decomposed_event_sum"),
            warnings,
            result,
            point,
            "threshold:decomposed_event_sum",
        )
        original = complex_value(
            summary.get("original"), warnings, result, point, "threshold:original"
        )
        counterterms = complex_value(
            summary.get("counterterm_sum"),
            warnings,
            result,
            point,
            "threshold:counterterm_sum",
        )
        reconstructed = complex_value(
            summary.get("reconstructed_total"),
            warnings,
            result,
            point,
            "threshold:reconstructed_total",
        )
        closure = complex_value(
            summary.get("closure"),
            warnings,
            result,
            point,
            "threshold:closure",
        )
        displayed = tuple(
            displayed_component(value, args.component)
            for value in (event_sum, original, counterterms, reconstructed)
        )
        closure_component = displayed_component(closure, args.component)
        scale = max(
            abs(displayed[0]),
            abs(displayed[1]) + abs(displayed[2]),
            sys.float_info.min,
        )
        relative_closure = abs(closure_component) / scale
        for label, value in zip(labels, (*displayed, relative_closure)):
            values[label].append((t_value, value))
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
    quantity: str,
    args: argparse.Namespace,
    warnings: NonFiniteWarnings,
) -> list[Series]:
    selected_keys = {
        key
        for component_keys in facet.component_groups.values()
        for key in component_keys
    }
    values: OrderedDict[str, list[tuple[float, float]]] = OrderedDict(
        (label, []) for label in facet.component_groups
    )
    catalog = threshold_component_catalog(result)
    for point, t_value in selected_evaluated_points(
        result, axis_index, branch, args.x_log_scale, args.signed_branches
    ):
        components = threshold_component_values(
            result, point, quantity, warnings, catalog, selected_keys
        )
        for label, component_keys in facet.component_groups.items():
            component_values = [components.get(key) for key in component_keys]
            if any(value is None for value in component_values):
                value = complex(math.nan, math.nan)
            else:
                value = sum(
                    (value for value in component_values if value is not None), 0j
                )
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


def runtime_occurrence_label(event: dict[str, Any], component: dict[str, Any]) -> str:
    occurrence = component.get("occurrence", {})
    occurrence_kind = occurrence.get("kind", "?")
    event_bits = [
        f"eg{event.get('event_group_index', '?')}/e{event.get('event_index', '?')}",
        f"c{event.get('cut_id', '?')}",
    ]
    if event.get("orientation_id") is not None:
        event_bits.append(f"o{event['orientation_id']}")
    if event.get("lmb_sample_id") is not None:
        event_bits.append(f"b{event['lmb_sample_id']}")
    if occurrence_kind == "local_unitarity":
        overlap_groups = ",".join(
            str(value) for value in occurrence.get("overlap_groups", [])
        )
        event_bits.extend(
            [
                f"og[{overlap_groups}]",
                f"L{occurrence.get('left_threshold_order', '-')}",
                f"R{occurrence.get('right_threshold_order', '-')}",
                f"u{occurrence.get('lu_cut_order', '-')}",
            ]
        )
    elif occurrence_kind == "amplitude":
        event_bits.extend(
            [
                f"es{occurrence.get('raised_esurface_id', '?')}",
                f"og{occurrence.get('overlap_group', '?')}",
            ]
        )
    else:
        event_bits.append(str(occurrence_kind))
    return " ".join(event_bits)


def collect_multiplier_facet_series(
    result: ResultFile,
    axis_index: int,
    branch: str,
    facet: ReportFacet,
    args: argparse.Namespace,
    warnings: NonFiniteWarnings,
    registry: ThresholdRegistry,
) -> list[Series]:
    component_keys = {key for keys in facet.component_groups.values() for key in keys}
    points = selected_evaluated_points(
        result, axis_index, branch, args.x_log_scale, args.signed_branches
    )
    values: OrderedDict[str, dict[int, float]] = OrderedDict()
    for point, _ in points:
        point_index = int(point["index"])
        for event in point.get("evaluation", {}).get("events", []):
            if not isinstance(event, dict):
                continue
            graph_id = int(event.get("graph_id", -1))
            decomposition = event.get("threshold_counterterms")
            if not isinstance(decomposition, dict):
                continue
            for raw_component in decomposition.get("components", []):
                if not isinstance(raw_component, dict):
                    continue
                component_id = int(raw_component.get("component_id", -1))
                key = graph_id, component_id
                if key not in component_keys:
                    continue
                metadata = registry.components[key]
                occurrence = runtime_occurrence_label(event, raw_component)
                multiplier_values = raw_component.get("multiplier_values", [])
                if not isinstance(multiplier_values, list):
                    multiplier_values = []
                for index, raw_value in enumerate(multiplier_values):
                    variant_id = (
                        metadata.variant_ids[index]
                        if index < len(metadata.variant_ids)
                        else -1
                    )
                    variant = registry.variant(graph_id, variant_id)
                    label = f"m{index}:f(v{variant_id} {variant.name}) · {occurrence}"
                    value = sanitized_float(
                        raw_value,
                        warnings,
                        result,
                        point,
                        label,
                        "value",
                    )
                    if point_index in values.setdefault(label, {}):
                        raise ValueError(
                            "multiplier runtime identity collision for "
                            f"{label!r} at point {point_index}"
                        )
                    values[label][point_index] = value
                effective_label = f"f_eff · {occurrence}"
                effective = sanitized_float(
                    raw_component.get("effective_multiplier"),
                    warnings,
                    result,
                    point,
                    effective_label,
                    "value",
                )
                if point_index in values.setdefault(effective_label, {}):
                    raise ValueError(
                        "effective-multiplier runtime identity collision for "
                        f"{effective_label!r} at point {point_index}"
                    )
                values[effective_label][point_index] = effective
    series = []
    for label, samples in values.items():
        y_values = []
        for point, _ in points:
            point_index = int(point["index"])
            if point_index not in samples:
                warnings.add(
                    result,
                    point,
                    label,
                    "value",
                    "missing runtime occurrence",
                )
            y_values.append(samples.get(point_index, math.nan))
        series.append(
            Series(
                label,
                [t_value for _, t_value in points],
                y_values,
            )
        )
    return series


def collect_supplemental_series(
    result: ResultFile,
    axis_index: int,
    branch: str,
    x_log_scale: bool,
    signed_x: bool,
    include_patterns: list[re.Pattern[str]],
    exclude_patterns: list[re.Pattern[str]],
    warnings: NonFiniteWarnings,
) -> list[Series]:
    supplemental = result.supplemental_series or []
    points = selected_evaluated_points(
        result, axis_index, branch, x_log_scale, signed_x
    )
    collected = []
    for external in supplemental:
        if external.axis_index != axis_index or not label_is_selected(
            external.label, include_patterns, exclude_patterns
        ):
            continue
        samples = []
        for point, t_value in points:
            point_index = int(point["index"])
            raw_value = external.values.get(point_index)
            samples.append(
                (
                    t_value,
                    sanitized_float(
                        raw_value,
                        warnings,
                        result,
                        point,
                        external.label,
                        "value",
                    ),
                )
            )
        collected.append(
            Series(
                external.label,
                [t_value for t_value, _ in samples],
                [value for _, value in samples],
            )
        )
    return collected


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
    signed_x: bool,
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
        t_value = selected_t_value(point, t_branch, x_log_scale, signed_x)
        if t_value is None:
            continue

        base_values: list[tuple[str, dict[str, Any]]] = [
            ("total_weight", evaluation.get("total_weight", {})),
            ("event_weight_sum", evaluation.get("event_weight_sum", {})),
        ]
        for weight_name, weight_value in sorted(
            evaluation.get("additional_contribution_sums", {}).items()
        ):
            base_values.append((f"additional:{weight_name}", weight_value))

        for label, serialized_value in base_values:
            if not label_is_selected(label, include_patterns, exclude_patterns):
                continue
            y_value = value_or_none(
                serialized_value, component, warnings, result, point, label
            )
            if y_value is not None:
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
    for supplemental in collect_supplemental_series(
        result,
        axis_index,
        t_branch,
        x_log_scale,
        signed_x,
        include_patterns,
        exclude_patterns,
        warnings,
    ):
        supplemental.label = f"{label_prefix}{supplemental.label}"
        series.append(supplemental)
    return series


def finite_abs(value: float) -> float:
    return abs(value)


def log_fit(
    series: Series,
    fit_points: int,
    x_range: list[float] | None = None,
    signed_x: bool = False,
) -> LogFit | None:
    points = [
        (math.log(abs(x)), math.log(abs(y)))
        for x, y in zip(series.t_values, series.values)
        if x != 0.0
        and y != 0.0
        and math.isfinite(x)
        and math.isfinite(y)
        and (
            x_range is None
            or (
                x_range[0] <= x <= x_range[1]
                if signed_x
                else x_range[0] <= abs(x) <= x_range[1]
            )
        )
    ]
    points.sort(key=lambda point: point[0])
    points = points[:fit_points]
    if len(points) < 2:
        return None
    mean_x = sum(x for x, _ in points) / len(points)
    mean_y = sum(y for _, y in points) / len(points)
    denominator = sum((x - mean_x) ** 2 for x, _ in points)
    if denominator == 0.0:
        return None
    slope = sum((x - mean_x) * (y - mean_y) for x, y in points) / denominator
    intercept = mean_y - slope * mean_x
    residual_sum = sum((y - (intercept + slope * x)) ** 2 for x, y in points)
    total_sum = sum((y - mean_y) ** 2 for _, y in points)
    r_squared = 1.0 if total_sum == 0.0 else 1.0 - residual_sum / total_sum
    return LogFit(slope, r_squared, len(points))


def plot_series(
    ax: plt.Axes,
    series: Series,
    component: str,
    y_log_scale: bool,
    fit_log_slope: bool,
    fit_points: int,
    x_range: list[float] | None,
    signed_x: bool,
) -> bool:
    if not series.t_values:
        return False
    label = series.label
    if (
        fit_log_slope
        and (fit := log_fit(series, fit_points, x_range, signed_x)) is not None
    ):
        label = (
            f"{label} (p={fit.slope:+.2f}, R^2={fit.r_squared:.3f}, "
            f"n={fit.point_count})"
        )
    marker_every = max(1, math.ceil(len(series.t_values) / 30))
    if not y_log_scale or component == "abs":
        y_values = [
            finite_abs(value) if y_log_scale else value for value in series.values
        ]
        ax.plot(
            series.t_values,
            y_values,
            linewidth=1.4,
            marker="o",
            markersize=2.5,
            markevery=marker_every,
            label=label,
        )
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
            marker="o",
            markersize=2.5,
            markevery=marker_every,
            label=label,
        )
    if any(math.isfinite(value) for value in negative):
        color = positive_line.get_color() if positive_line is not None else None
        ax.plot(
            series.t_values,
            negative,
            linewidth=1.4,
            linestyle="--",
            marker="o",
            markersize=2.5,
            markevery=marker_every,
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
            (
                f"spacing: {spacing.get('kind', '?')} n={spacing.get('n_points', '?')} "
                f"|{parameter_name(result)}|≤{float(spacing.get('max_abs_t')):.4g}"
            ),
            f"source: {result.path.name}",
        ]
    )


def signed_log_bounds(
    results: Iterable[ResultFile], args: argparse.Namespace
) -> tuple[float, float | None]:
    linthresh_candidates = []
    extent_candidates = []
    for result in results:
        spacing = result.data.get("spacing", {})
        min_abs_t = spacing.get("min_abs_t")
        if (
            isinstance(min_abs_t, (int, float))
            and math.isfinite(float(min_abs_t))
            and float(min_abs_t) > 0.0
        ):
            linthresh_candidates.append(float(min_abs_t))
        else:
            linthresh_candidates.extend(
                abs(float(point["t"]))
                for point in result.data.get("points", [])
                if math.isfinite(float(point.get("t", 0.0)))
                and float(point.get("t", 0.0)) != 0.0
            )
        extent_candidates.append(float(spacing["max_abs_t"]))
    if args.x_symlog_linthresh is not None:
        linthresh = args.x_symlog_linthresh
    elif linthresh_candidates:
        linthresh = min(linthresh_candidates)
    else:
        raise ValueError(
            "Signed logarithmic plotting requires a positive spacing.min_abs_t "
            "or at least one non-zero approach point."
        )
    return linthresh, None if args.x_range is not None else max(extent_candidates)


def decorate_axis(
    ax: plt.Axes,
    result: ResultFile,
    component: str,
    y_log_scale: bool,
    x_log_scale: bool,
    hide_info_box: bool,
    y_label: str | None = None,
    branch: str | None = None,
    signed_x: bool = False,
    signed_linthresh: float | None = None,
    show_ylabel: bool = True,
    show_xlabel: bool = True,
) -> None:
    parameter = parameter_name(result)
    if show_xlabel:
        if signed_x:
            ax.set_xlabel(parameter)
        elif x_log_scale and branch in ("negative", "positive"):
            sign = "< 0" if branch == "negative" else "> 0"
            ax.set_xlabel(f"|{parameter}|  ({parameter} {sign})")
        else:
            ax.set_xlabel(f"|{parameter}|" if x_log_scale else parameter)
    else:
        ax.tick_params(axis="x", labelbottom=False)
    has_nonzero_x = any(
        math.isfinite(float(value)) and float(value) != 0.0
        for line in ax.lines
        for value in line.get_xdata()
    )
    if signed_x and has_nonzero_x:
        if signed_linthresh is None:
            raise ValueError("Signed logarithmic plotting requires a linthresh.")
        ax.set_xscale("symlog", linthresh=signed_linthresh)
    elif x_log_scale and has_nonzero_x:
        ax.set_xscale("log")
    ylabel = y_label or (
        "weight magnitude" if component == "abs" else f"{component} weight"
    )
    if y_log_scale and any(
        math.isfinite(float(value)) and float(value) > 0.0
        for line in ax.lines
        for value in line.get_ydata()
    ):
        if component != "abs" and y_label is None:
            ylabel = f"|{ylabel}|"
        ax.set_yscale("log")
    if show_ylabel:
        ax.set_ylabel(ylabel)
    if not x_log_scale or signed_x:
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
    ax: plt.Axes,
    args: argparse.Namespace,
    *,
    reverse_x: bool = False,
    signed_extent: float | None = None,
) -> None:
    if signed_extent is not None:
        ax.set_xlim(-signed_extent, signed_extent)
    elif args.x_range is not None:
        limits = args.x_range[::-1] if reverse_x else args.x_range
        ax.set_xlim(limits[0], limits[1])
    elif reverse_x:
        ax.invert_xaxis()
    if args.y_range is not None:
        ax.set_ylim(args.y_range[0], args.y_range[1])


def branches_on_axis(
    args: argparse.Namespace, result: ResultFile
) -> list[tuple[str, str]]:
    if args.t_branch != "both":
        return [(args.t_branch, "")]
    if args.x_log_scale:
        parameter = parameter_name(result)
        return [("negative", f"{parameter}<0"), ("positive", f"{parameter}>0")]
    return [("both", "")]


def add_figure_footer(fig: plt.Figure, result: ResultFile, axis_index: int) -> None:
    fig.text(
        0.01,
        0.006,
        (
            f"source: {result.path.name}  ·  axis {axis_index + 1}: "
            f"{axis_label(result, axis_index)}"
        ),
        ha="left",
        va="bottom",
        fontsize=6.2,
        color="0.35",
    )


def wrap_display_text(value: str, width: int) -> str:
    return "\n".join(
        wrapped
        for line in value.splitlines() or [value]
        for wrapped in (
            textwrap.wrap(
                line,
                width=width,
                break_long_words=True,
                break_on_hyphens=False,
            )
            or [""]
        )
    )


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
    font_size = 6.0 if len(entries) > 8 else 6.8
    label_width = 42 if len(entries) > 8 else 48
    wrapped_entries = OrderedDict(
        (wrap_display_text(label, label_width), handle)
        for label, handle in entries.items()
    )
    widest_line = max(
        len(line) for label in wrapped_entries for line in label.splitlines()
    )
    figure_width = float(fig.get_size_inches()[0])
    approximate_character_capacity = figure_width * 72.0 / (font_size * 0.58)
    columns = min(
        len(wrapped_entries),
        4,
        max(1, int(approximate_character_capacity / (widest_line + 8))),
    )
    rows = math.ceil(len(entries) / columns)
    wrapped_lines = max(label.count("\n") + 1 for label in wrapped_entries)
    fig.legend(
        list(wrapped_entries.values()),
        list(wrapped_entries),
        loc="lower center",
        bbox_to_anchor=(0.5, 0.022),
        fontsize=font_size,
        frameon=True,
        ncols=columns,
        columnspacing=0.8,
        handletextpad=0.45,
        borderaxespad=0.2,
    )
    return min(0.48, 0.085 + 0.031 * rows * wrapped_lines)


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
        parts.append(f"axis {axis_index + 1}: {axis_label(result, axis_index)}")
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
    signed_linthresh, signed_extent = (
        signed_log_bounds((result for result, _ in result_axis_pairs), args)
        if args.signed_branches
        else (None, None)
    )
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
    has_nonzero_values = False
    multi_result = len({pair[0].path for pair in result_axis_pairs}) > 1
    multi_axis = len({pair[1] for pair in result_axis_pairs}) > 1

    for column, ax in enumerate(axes):
        branch_specs = (
            ([("negative", "")] if column == 0 else [("positive", "")])
            if args.split_branches
            else branches_on_axis(args, result_axis_pairs[0][0])
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
                    args.signed_branches,
                    args.threshold_decomposition,
                    label_prefix=label_prefix,
                )
                for item in series:
                    if branch_suffix:
                        item.label = f"{item.label} · {branch_suffix}"
                    has_series = axis_has_series = True
                    has_nonzero_values |= any(
                        math.isfinite(value) and value != 0.0 for value in item.values
                    )
                    has_negative_values |= plot_series(
                        ax,
                        item,
                        args.component,
                        args.y_log_scale,
                        args.fit_log_slope,
                        args.fit_points,
                        args.x_range,
                        args.signed_branches,
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
    effective_y_log_scale = args.y_log_scale and has_nonzero_values
    title = args.title or axis_label(first_result, first_axis)
    if args.title is None and len(result_axis_pairs) > 1:
        title = "combined approach curves"
    fig.suptitle(title, fontsize=12, fontweight="bold")
    for column, ax in enumerate(axes):
        branch = ("negative", "positive")[column] if args.split_branches else None
        if args.split_branches:
            parameter = parameter_name(first_result)
            ax.set_title(
                f"{parameter} < 0" if branch == "negative" else f"{parameter} > 0",
                fontsize=9,
            )
        decorate_axis(
            ax,
            first_result,
            args.component,
            effective_y_log_scale,
            args.x_log_scale,
            args.hide_info_box or column > 0,
            y_label=args.y_label,
            branch=branch,
            signed_x=args.signed_branches,
            signed_linthresh=signed_linthresh,
            show_ylabel=column == 0,
        )
        apply_axis_ranges(
            ax,
            args,
            reverse_x=branch == "negative",
            signed_extent=signed_extent,
        )
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
    add_figure_footer(fig, first_result, first_axis)
    pdf.savefig(fig)
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
    quantity: str = "weighted",
    multiplier: bool = False,
) -> None:
    row_count = len(facets)
    column_count = 2 if args.split_branches else 1
    signed_linthresh, signed_extent = (
        signed_log_bounds((result,), args) if args.signed_branches else (None, None)
    )
    fig, axes = plt.subplots(
        row_count,
        column_count,
        figsize=(13.5, 8.2),
        sharey="row" if column_count == 2 else False,
        squeeze=False,
    )
    has_negative_values = False
    rows_with_nonzero_values = [False] * row_count
    y_log_scale = args.multiplier_log_scale if multiplier else args.y_log_scale
    registry = ThresholdRegistry.from_result(result)
    flat_axes = [ax for row in axes for ax in row]
    panel_title_width = 50 if column_count == 2 else 96
    panel_titles = [
        wrap_display_text(facet.title, panel_title_width) for facet in facets
    ]
    for row, facet in enumerate(facets):
        for column in range(column_count):
            ax = axes[row][column]
            branch_specs = (
                ([("negative", "")] if column == 0 else [("positive", "")])
                if args.split_branches
                else branches_on_axis(args, result)
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
                    else collect_multiplier_facet_series(
                        result,
                        axis_index,
                        branch,
                        facet,
                        args,
                        warnings,
                        registry,
                    )
                    if multiplier
                    else collect_threshold_facet_series(
                        result,
                        axis_index,
                        branch,
                        facet,
                        component,
                        quantity,
                        args,
                        warnings,
                    )
                )
                for item in series:
                    if branch_suffix:
                        item.label = f"{item.label} · {branch_suffix}"
                    has_series = True
                    rows_with_nonzero_values[row] |= any(
                        math.isfinite(value) and value != 0.0 for value in item.values
                    )
                    has_negative_values |= plot_series(
                        ax,
                        item,
                        component,
                        y_log_scale,
                        args.fit_log_slope and y_log_scale and args.x_log_scale,
                        args.fit_points,
                        args.x_range,
                        args.signed_branches,
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
            axis_title = panel_titles[row]
            if args.split_branches:
                parameter = parameter_name(result)
                axis_title += "\n" + (
                    f"{parameter} < 0" if branch == "negative" else f"{parameter} > 0"
                )
            ax.set_title(
                axis_title,
                fontsize=7.6 if axis_title.count("\n") >= 3 else 8.2,
                fontweight="semibold",
                linespacing=1.12,
            )
            decorate_axis(
                ax,
                result,
                component,
                y_log_scale,
                args.x_log_scale,
                args.hide_info_box or row > 0 or column > 0,
                y_label=(
                    args.y_label
                    or facet.y_label
                    or ("multiplier value" if multiplier else None)
                    or ("unmultiplied CT" if quantity == "bare" else None)
                ),
                branch=branch,
                signed_x=args.signed_branches,
                signed_linthresh=signed_linthresh,
                show_ylabel=column == 0,
                show_xlabel=row + 1 == row_count,
            )
            apply_axis_ranges(
                ax,
                args,
                reverse_x=branch == "negative",
                signed_extent=signed_extent,
            )

    if y_log_scale:
        for row, has_nonzero_values in enumerate(rows_with_nonzero_values):
            if has_nonzero_values:
                continue
            for column in range(column_count):
                axes[row][column].set_yscale("linear")

    fig.suptitle(
        wrap_display_text(title, 112),
        fontsize=12,
        fontweight="bold",
    )
    bottom = shared_figure_legend(fig, flat_axes, has_negative_values)
    panel_title_lines = max(panel_title.count("\n") + 1 for panel_title in panel_titles)
    if args.split_branches:
        panel_title_lines += 1
    extra_title_lines = max(0, panel_title_lines - 2)
    top = max(bottom + 0.2, 0.89 - 0.02 * extra_title_lines)
    fig.subplots_adjust(
        # Leave room for long scientific-notation tick labels together with the
        # shared quantity label on component pages.
        left=0.1,
        right=0.98,
        top=top,
        bottom=bottom,
        hspace=0.5 + 0.18 * extra_title_lines,
        wspace=0.06 if args.split_branches else 0.2,
    )
    add_figure_footer(fig, result, axis_index)
    pdf.savefig(fig)
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
    include_patterns: list[re.Pattern[str]],
    exclude_patterns: list[re.Pattern[str]],
) -> int:
    page_count = 0
    for result in results:
        registry = ThresholdRegistry.from_result(result)
        for axis_index in axis_indices(result, args.selected_axis_indices):
            all_observed = registry.observed_components(result, axis_index)
            observed = registry.filter_components(
                all_observed, include_patterns, exclude_patterns
            )
            if (
                "summary" in args.threshold_report_sections
                and axis_has_threshold_summary(result, axis_index)
            ):
                summary_facets = [
                    ReportFacet(
                        "event decomposition",
                        OrderedDict(),
                        (
                            "decomposed event sum",
                            "original O",
                            "Σ weighted CT",
                            "reconstructed O+ΣCT",
                        ),
                        (
                            "weight magnitude"
                            if args.component == "abs"
                            else f"{args.component} weight"
                        ),
                    ),
                    ReportFacet(
                        "relative closure residual",
                        OrderedDict(),
                        ("relative closure |Δ|/scale",),
                        "relative residual",
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

            if "singles" in args.threshold_report_sections and observed:
                single_facets = registry.single_facets(observed)
                for quantity in args.threshold_quantities:
                    chunks = [
                        single_facets[index : index + args.facets_per_page]
                        for index in range(0, len(single_facets), args.facets_per_page)
                    ]
                    for chunk_index, facets in enumerate(chunks):
                        suffix = (
                            f" {chunk_index + 1}/{len(chunks)}"
                            if len(chunks) > 1
                            else ""
                        )
                        quantity_label = (
                            "weighted" if quantity == "weighted" else "unmultiplied CT"
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
                                f"variants · {quantity_label}{suffix}",
                            ),
                            component=args.component,
                            quantity=quantity,
                        )
                        page_count += 1

            if "pairs" in args.threshold_report_sections and observed:
                pair_facets = registry.pair_facets(observed)
                for quantity in args.threshold_quantities:
                    chunks = [
                        pair_facets[index : index + 2]
                        for index in range(0, len(pair_facets), 2)
                    ]
                    for chunk_index, facets in enumerate(chunks):
                        suffix = (
                            f" {chunk_index + 1}/{len(chunks)}"
                            if len(chunks) > 1
                            else ""
                        )
                        quantity_label = (
                            "weighted" if quantity == "weighted" else "unmultiplied CT"
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
                                f"pairs · {quantity_label}{suffix}",
                            ),
                            component=args.component,
                            quantity=quantity,
                        )
                        page_count += 1

            if "multipliers" in args.threshold_report_sections and observed:
                multiplier_facets = registry.multiplier_facets(observed)
                chunks = [
                    multiplier_facets[index : index + 2]
                    for index in range(0, len(multiplier_facets), 2)
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
                            f"multipliers{suffix}",
                        ),
                        component="real",
                        multiplier=True,
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
        if args.series_json:
            load_supplemental_series(args.series_json, results[0])
        warnings = NonFiniteWarnings()
        include_contributions = list(args.include_contribution)
        if not include_contributions:
            include_contributions.append(r"^total_weight$")
            if args.threshold_decomposition is not None:
                include_contributions.append(r"^threshold:")
        include_patterns = compile_patterns(include_contributions)
        exclude_patterns = compile_patterns(args.exclude_contribution)
        report_include_patterns = compile_patterns(args.include_contribution)
        groups = page_groups(results, args)
        if not groups:
            raise ValueError(
                "No axes found in the supplied approach result JSON files."
            )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with PdfPages(args.output) as pdf:
            if args.threshold_report:
                if (
                    draw_threshold_report(
                        pdf,
                        results,
                        args,
                        warnings,
                        report_include_patterns,
                        exclude_patterns,
                    )
                    == 0
                ):
                    raise ValueError(
                        "No requested threshold-counterterm data was found. Summary "
                        "pages require evaluation.threshold_counterterm_summary; "
                        "component pages require observed registry components matching "
                        "the requested filters."
                    )
            else:
                for group in groups:
                    draw_page(
                        pdf, group, args, include_patterns, exclude_patterns, warnings
                    )
        warnings.emit(args.max_gap_warnings)
        print(f"plot_approach_result.py: PDF created at {display_path(args.output)}")
    except Exception as error:  # noqa: BLE001 - this script should print clean CLI errors.
        print(f"plot_approach_result.py: error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
