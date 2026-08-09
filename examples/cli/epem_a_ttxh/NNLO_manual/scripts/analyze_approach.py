#!/usr/bin/env python3
"""Fit signed GammaLoop approach branches and classify an IR limit."""

from __future__ import annotations

import argparse
import json
import math
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass
class Fit:
    branch: str
    window: int
    points: int
    slope: float
    p: float
    r_squared: float
    min_abs_lambda: float
    max_abs_lambda: float


def parse_windows(raw: str) -> list[int]:
    windows = sorted({int(value) for value in raw.split(",")})
    if not windows or windows[0] < 2:
        raise argparse.ArgumentTypeError("fit windows must contain integers >= 2")
    return windows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--axis-index", type=int, default=0)
    parser.add_argument("--rank", type=float, required=True)
    parser.add_argument(
        "--windows", type=parse_windows, default=parse_windows("8,12,16")
    )
    parser.add_argument("--minimum-r-squared", type=float, default=0.995)
    parser.add_argument("--maximum-p-span", type=float, default=0.5)
    parser.add_argument("--margin", type=float, default=0.25)
    parser.add_argument("--threshold-only", action="store_true")
    parser.add_argument("--series-json", type=Path)
    parser.add_argument("--no-kinematic-check", action="store_true")
    parser.add_argument("--kinematic-slope-tolerance", type=float, default=0.05)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def linear_fit(xs: Iterable[float], ys: Iterable[float]) -> tuple[float, float]:
    x_values = list(xs)
    y_values = list(ys)
    count = len(x_values)
    if count < 2 or count != len(y_values):
        raise ValueError(
            "a linear fit requires equally sized inputs with at least two points"
        )
    mean_x = sum(x_values) / count
    mean_y = sum(y_values) / count
    variance = sum((value - mean_x) ** 2 for value in x_values)
    if variance == 0.0:
        raise ValueError("cannot fit identical x values")
    slope = (
        sum(
            (x_value - mean_x) * (y_value - mean_y)
            for x_value, y_value in zip(x_values, y_values)
        )
        / variance
    )
    intercept = mean_y - slope * mean_x
    residual = sum(
        (y_value - (intercept + slope * x_value)) ** 2
        for x_value, y_value in zip(x_values, y_values)
    )
    total = sum((value - mean_y) ** 2 for value in y_values)
    r_squared = 1.0 if total == 0.0 and residual == 0.0 else 1.0 - residual / total
    return slope, r_squared


def evaluated_points(raw: dict, axis_index: int) -> list[dict]:
    selected = []
    for point in raw.get("points", []):
        if point.get("axis_index") != axis_index or point.get("status") != "evaluated":
            continue
        t = point.get("t")
        evaluation = point.get("evaluation") or {}
        total = evaluation.get("total_weight") or {}
        metadata = evaluation.get("metadata") or {}
        re_value = total.get("re")
        im_value = total.get("im")
        if metadata.get("is_nan"):
            continue
        if not all(
            isinstance(value, (int, float)) and math.isfinite(value)
            for value in (t, re_value, im_value)
        ):
            continue
        magnitude = math.hypot(re_value, im_value)
        if t == 0.0 or magnitude <= 0.0 or not math.isfinite(magnitude):
            continue
        selected.append(
            {
                "index": point.get("index"),
                "t": float(t),
                "value": magnitude,
            }
        )
    return selected


def has_active_threshold_counterterm(raw: dict, axis_index: int) -> bool:
    for point in raw.get("points", []):
        if point.get("axis_index") != axis_index or point.get("status") != "evaluated":
            continue
        summary = (point.get("evaluation") or {}).get(
            "threshold_counterterm_summary"
        ) or {}
        counterterm_sum = summary.get("counterterm_sum") or {}
        values = (counterterm_sum.get("re"), counterterm_sum.get("im"))
        if (
            all(
                isinstance(value, (int, float)) and math.isfinite(value)
                for value in values
            )
            and math.hypot(*values) > 0.0
        ):
            return True
    return False


def fit_branches(points: list[dict], windows: list[int]) -> tuple[list[Fit], list[str]]:
    fits: list[Fit] = []
    errors: list[str] = []
    for branch, predicate in (
        ("negative", lambda t: t < 0),
        ("positive", lambda t: t > 0),
    ):
        branch_points = sorted(
            (point for point in points if predicate(point["t"])),
            key=lambda point: abs(point["t"]),
        )
        for window in windows:
            if len(branch_points) < window:
                errors.append(
                    f"{branch} branch has {len(branch_points)} usable points, fewer than {window}"
                )
                continue
            selected = branch_points[:window]
            slope, r_squared = linear_fit(
                (math.log(abs(point["t"])) for point in selected),
                (math.log(point["value"]) for point in selected),
            )
            fits.append(
                Fit(
                    branch=branch,
                    window=window,
                    points=len(selected),
                    slope=slope,
                    p=-slope,
                    r_squared=r_squared,
                    min_abs_lambda=min(abs(point["t"]) for point in selected),
                    max_abs_lambda=max(abs(point["t"]) for point in selected),
                )
            )
    return fits, errors


def kinematic_fits(
    series_path: Path | None,
    raw_points: list[dict],
    axis_index: int,
    windows: list[int],
) -> tuple[list[dict], list[str]]:
    if series_path is None:
        return [], ["no kinematic series supplied"]
    raw_series = json.loads(series_path.read_text())
    t_by_index = {point["index"]: point["t"] for point in raw_points}
    output: list[dict] = []
    errors: list[str] = []
    for series in raw_series.get("series", []):
        if series.get("axis_index", 0) != axis_index:
            continue
        samples = []
        for value in series.get("values", []):
            t = t_by_index.get(value.get("point_index"))
            distance = value.get("value")
            if (
                t is None
                or t == 0.0
                or not isinstance(distance, (int, float))
                or not math.isfinite(distance)
                or distance <= 0.0
            ):
                continue
            samples.append({"t": t, "value": float(distance)})
        fits, fit_errors = fit_branches(samples, windows)
        output.append(
            {
                "label": series.get("label", "unnamed"),
                "fits": [asdict(fit) for fit in fits],
            }
        )
        errors.extend(
            f"{series.get('label', 'unnamed')}: {error}" for error in fit_errors
        )
    if not output:
        errors.append("kinematic series contains no data for the selected axis")
    return output, errors


def main() -> int:
    args = parse_args()
    raw = json.loads(args.input.read_text())
    points = evaluated_points(raw, args.axis_index)
    fits, errors = fit_branches(points, args.windows)
    kinematics, kinematic_errors = kinematic_fits(
        args.series_json, points, args.axis_index, args.windows
    )
    if not args.no_kinematic_check:
        errors.extend(kinematic_errors)

    p_values = [fit.p for fit in fits]
    p_min = min(p_values) if p_values else None
    p_max = max(p_values) if p_values else None
    p_span = p_max - p_min if p_values else None
    poor_fits = [fit for fit in fits if fit.r_squared < args.minimum_r_squared]
    if poor_fits:
        errors.append(
            f"{len(poor_fits)} total-weight fits have R^2 below {args.minimum_r_squared}"
        )
    if p_span is not None and p_span > args.maximum_p_span:
        errors.append(
            f"total-weight p envelope span {p_span:.6g} exceeds {args.maximum_p_span}"
        )
    active_threshold_counterterm = has_active_threshold_counterterm(
        raw, args.axis_index
    )
    if args.threshold_only and not active_threshold_counterterm:
        errors.append(
            "threshold-only path has no nonzero threshold-counterterm contribution"
        )

    kinematic_fit_rows = [fit for series in kinematics for fit in series["fits"]]
    if not args.no_kinematic_check and any(
        fit["r_squared"] < args.minimum_r_squared
        or abs(fit["slope"] - 1.0) > args.kinematic_slope_tolerance
        for fit in kinematic_fit_rows
    ):
        errors.append("one or more target distances do not scale linearly with lambda")

    classification = "inconclusive"
    if not errors and p_min is not None and p_max is not None:
        if args.threshold_only:
            if p_max <= args.margin:
                classification = "pass"
            elif p_min > args.margin:
                classification = "fail"
        elif p_max <= args.rank - args.margin:
            classification = "pass"
        elif p_min >= args.rank + args.margin:
            classification = "fail"

    result = {
        "schema_version": 1,
        "source": str(args.input),
        "axis_index": args.axis_index,
        "rank": args.rank,
        "threshold_only": args.threshold_only,
        "active_threshold_counterterm": active_threshold_counterterm,
        "classification": classification,
        "p_envelope": {"min": p_min, "max": p_max, "span": p_span},
        "fits": [asdict(fit) for fit in fits],
        "kinematic_fits": kinematics,
        "errors": errors,
        "criteria": {
            "windows": args.windows,
            "minimum_r_squared": args.minimum_r_squared,
            "maximum_p_span": args.maximum_p_span,
            "margin": args.margin,
            "kinematic_slope_tolerance": args.kinematic_slope_tolerance,
        },
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(serialized, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    return 0 if classification != "fail" else 2


if __name__ == "__main__":
    raise SystemExit(main())
