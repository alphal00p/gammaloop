#!/usr/bin/env python3
"""Scan approach lines through the latest high-stat imaginary max-weight points."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
from gammaloop import GammaLoopAPI

from scan_imag_maxweight_approaches import (
    OUTPUT_DIR,
    STATE_DIR,
    Direction,
    Target,
    add_arb_peak_checks,
    apply_runtime_settings,
    build_line_points,
    evaluate_rows,
    plot_target,
    write_inspect_spotchecks,
    write_json_tsv,
)
from inspect_imag_maxweight_thresholds import (
    ProbePoint,
    evaluate_probe,
    threshold_maps,
    write_markdown,
    write_tsv as write_threshold_tsv,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_WORKSPACE = (
    ROOT
    / "examples"
    / "cli"
    / "qq_hhh_2L"
    / "workspaces"
    / "integrate_xi3000_common_radial_lmb_mc_imag_10pct"
)
COORD_RE = re.compile(
    r"graph:\s*(?P<graph>\d+),\s*LMB channel:\s*(?P<channel>\d+),\s*xs:\s*\[(?P<xs>.*?)\]",
    re.S,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Wait for a high-stat integration iteration, extract its im[+]/im[-] "
            "max weights, and run fixed-bin approach scans plus threshold "
            "decomposition."
        )
    )
    parser.add_argument(
        "--workspace",
        type=Path,
        default=DEFAULT_WORKSPACE,
        help="Integration workspace to watch for result JSON files.",
    )
    parser.add_argument(
        "--result-json",
        type=Path,
        default=None,
        help="Use a specific integration_result_iter_*.json instead of waiting.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR / "highstat_current",
        help="Directory for scan and threshold-decomposition outputs.",
    )
    parser.add_argument(
        "--state-dir",
        type=Path,
        default=STATE_DIR,
        help="GammaLoop state directory to load read-only.",
    )
    parser.add_argument(
        "--points-per-line",
        type=int,
        default=101,
        help="Odd number of approach points per line, including the peak.",
    )
    parser.add_argument(
        "--poll-seconds",
        type=float,
        default=30.0,
        help="Polling interval while waiting for an integration result.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=0.0,
        help="Optional timeout while waiting. Zero means wait indefinitely.",
    )
    parser.add_argument(
        "--skip-arb-peak-checks",
        action="store_true",
        help="Skip arb checks at the exact max-weight points.",
    )
    parser.add_argument(
        "--skip-inspect-spotchecks",
        action="store_true",
        help="Skip CLI inspect spotcheck logs.",
    )
    parser.add_argument(
        "--skip-threshold-decomposition",
        action="store_true",
        help="Skip event-level additional-weight decomposition at the peaks.",
    )
    parser.add_argument(
        "--variant-label",
        default=None,
        help="Optional label stored in the output metadata for width scans.",
    )
    parser.add_argument(
        "--sliver-width",
        type=float,
        default=None,
        help="Override subtraction.local_ct_settings.uv_localisation.sliver_width.",
    )
    parser.add_argument(
        "--gaussian-width",
        type=float,
        default=None,
        help="Override subtraction.local_ct_settings.uv_localisation.gaussian_width.",
    )
    parser.add_argument(
        "--h-sigma",
        type=float,
        default=None,
        help=(
            "Override subtraction.integrated_ct_settings.range."
            "h_function_settings.sigma only. This intentionally does not "
            "touch the top-level LU h_function.sigma."
        ),
    )
    parser.add_argument(
        "--xi-scale",
        type=float,
        default=3.0,
        help="Use xi = xi_scale * (p1 + p2) in Python API evaluations.",
    )
    return parser.parse_args()


def wait_for_latest_result(workspace: Path, poll_seconds: float, timeout_seconds: float) -> Path:
    results_dir = workspace / "results"
    start = time.monotonic()
    while True:
        candidates = sorted(results_dir.glob("integration_result_iter_*.json"))
        if candidates:
            latest = candidates[-1]
            # Avoid reading a file still being written.
            size_before = latest.stat().st_size
            time.sleep(0.5)
            size_after = latest.stat().st_size
            if size_before == size_after and size_after > 0:
                return latest
        if timeout_seconds and time.monotonic() - start > timeout_seconds:
            raise TimeoutError(f"No integration result appeared under {results_dir}")
        print(f"waiting for integration result under {results_dir}", flush=True)
        time.sleep(poll_seconds)


def parse_coordinates(text: str) -> tuple[int, int, tuple[float, float, float, float, float, float]]:
    match = COORD_RE.search(text)
    if not match:
        raise ValueError(f"Cannot parse max-weight coordinates: {text}")
    xs = tuple(float(value) for value in match.group("xs").split())
    if len(xs) != 6:
        raise ValueError(f"Expected six x coordinates, got {len(xs)} from {text}")
    return int(match.group("graph")), int(match.group("channel")), xs  # type: ignore[return-value]


def nearest_boundary_axes(x_peak: tuple[float, ...]) -> tuple[int, int]:
    distances = sorted((min(x, 1.0 - x), index) for index, x in enumerate(x_peak))
    return distances[0][1], distances[1][1]


def targets_from_result(path: Path) -> list[Target]:
    data = json.loads(path.read_text())
    if not data.get("slots"):
        raise ValueError(f"No slots found in {path}")
    slot = data["slots"][0]
    iteration = path.stem.removeprefix("integration_result_iter_")
    targets: list[Target] = []
    for entry in slot.get("max_weight_info", []):
        if entry.get("component") != "im":
            continue
        graph, channel, xs = parse_coordinates(str(entry["coordinates"]))
        sign = str(entry["sign"])
        label_sign = "plus" if sign == "+" else "minus"
        label = f"highstat_iter_{iteration}_imag_{label_sign}"
        targets.append(
            Target(
                label=label,
                reported_component=f"im[{sign}]",
                reported_value=float(entry["max_eval"]),
                graph=graph,
                lmb_channel=channel,
                x_peak=xs,
                axis_indices=nearest_boundary_axes(xs),
            )
        )
    if len(targets) != 2:
        raise ValueError(f"Expected im[+]/im[-] max points in {path}, found {len(targets)}")
    targets.sort(key=lambda target: target.reported_component)
    return targets


def directions_for_target(target: Target) -> list[Direction]:
    directions = []
    for axis_index in target.axis_indices:
        vector = np.zeros(6, dtype=float)
        vector[axis_index] = 1.0
        directions.append(Direction(f"axis_x{axis_index}", vector))
    diagonal = np.array(target.x_peak, dtype=float) - 0.5
    norm = float(np.linalg.norm(diagonal))
    if norm == 0.0:
        diagonal = np.ones(6, dtype=float)
        norm = float(np.linalg.norm(diagonal))
    directions.append(Direction("center_diagonal", diagonal / norm))
    return directions


def run_approach_scan(
    targets: list[Target],
    output_dir: Path,
    state_dir: Path,
    points_per_line: int,
    skip_arb_peak_checks: bool,
    skip_inspect_spotchecks: bool,
    overrides: dict[str, float],
    xi_scale: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    groups: list[tuple[Target, Direction, list[dict[str, Any]]]] = []
    for target in targets:
        for direction in directions_for_target(target):
            groups.append((target, direction, build_line_points(target, direction, points_per_line)))

    api = GammaLoopAPI(state_folder=state_dir, read_only_state=True)
    apply_runtime_settings(api, xi_scale=xi_scale)
    apply_width_overrides(api, overrides)
    for index, (target, direction, group_rows) in enumerate(groups, start=1):
        print(
            f"approach scan {index}/{len(groups)}: {target.label} {direction.label}",
            flush=True,
        )
        evaluate_rows(api, group_rows, batch_size=1)
        rows.extend(group_rows)

    if not skip_arb_peak_checks:
        api = GammaLoopAPI(state_folder=state_dir, read_only_state=True)
        apply_runtime_settings(api, xi_scale=xi_scale)
        apply_width_overrides(api, overrides)
        add_arb_peak_checks(api, rows)

    write_json_tsv(rows, output_dir)
    for target in targets:
        plot_target(rows, target, output_dir)

    if not skip_inspect_spotchecks:
        write_inspect_spotchecks(targets, output_dir, state_dir, ROOT / "target" / "dev-optim" / "gammaloop")
    return rows


def run_threshold_decomposition(
    targets: list[Target], output_dir: Path, state_dir: Path, xi_scale: float
) -> None:
    report_dir = output_dir / "threshold_decomposition"
    report_dir.mkdir(parents=True, exist_ok=True)
    api = GammaLoopAPI(state_folder=state_dir, read_only_state=True)
    apply_runtime_settings(api, xi_scale=xi_scale)
    api.run(
        "set process -p qq_hhh_2L -i subtracted kv "
        "general.generate_events=true "
        "general.store_additional_weights_in_event=true "
        "general.enable_cache=false "
        "integrator.integrated_phase=\"imag\""
    )
    group_maps = threshold_maps(api)
    reports = []
    for index, target in enumerate(targets, start=1):
        print(f"threshold decomposition {index}/{len(targets)}: {target.label}", flush=True)
        probe = ProbePoint(
            label=f"{target.label}_exact_highstat_max",
            target=target.label,
            graph_group=target.graph,
            lmb_channel=target.lmb_channel,
            x=target.x_peak,
            source=f"{target.reported_component} {target.reported_value:+.8e}",
        )
        reports.append(evaluate_probe(api, group_maps, probe))
    (report_dir / "threshold_decomposition.json").write_text(json.dumps(reports, indent=2) + "\n")
    write_threshold_tsv(reports, report_dir / "threshold_decomposition.tsv")
    write_markdown(reports, report_dir / "threshold_decomposition.md")


def apply_width_overrides(api: GammaLoopAPI, overrides: dict[str, float]) -> None:
    if not overrides:
        return
    assignments = []
    if "sliver_width" in overrides:
        assignments.append(
            "subtraction.local_ct_settings.uv_localisation.sliver_width="
            f"{overrides['sliver_width']:.17g}"
        )
    if "gaussian_width" in overrides:
        assignments.append(
            "subtraction.local_ct_settings.uv_localisation.gaussian_width="
            f"{overrides['gaussian_width']:.17g}"
        )
    if "h_sigma" in overrides:
        h_sigma = overrides["h_sigma"]
        api.run(
            f"""set process -p qq_hhh_2L -i subtracted string '
[subtraction.integrated_ct_settings.range]
type = "infinite"

[subtraction.integrated_ct_settings.range.h_function_settings]
function = "poly_exponential"
sigma = {h_sigma:.17g}
enabled_dampening = true
'"""
        )
    if assignments:
        api.run("set process -p qq_hhh_2L -i subtracted kv " + " ".join(assignments))


def main() -> None:
    args = parse_args()
    result_json = (
        args.result_json
        if args.result_json is not None
        else wait_for_latest_result(args.workspace, args.poll_seconds, args.timeout_seconds)
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    copied_result = args.output_dir / result_json.name
    if result_json.resolve() != copied_result.resolve():
        shutil.copy2(result_json, copied_result)

    targets = targets_from_result(result_json)
    overrides = {
        key: value
        for key, value in {
            "sliver_width": args.sliver_width,
            "gaussian_width": args.gaussian_width,
            "h_sigma": args.h_sigma,
        }.items()
        if value is not None
    }
    (args.output_dir / "highstat_imag_max_points.json").write_text(
        json.dumps(
            {
                "variant_label": args.variant_label,
                "width_overrides": overrides,
                "xi_scale": args.xi_scale,
                "targets": [
                    {
                        "label": target.label,
                        "reported_component": target.reported_component,
                        "reported_value": target.reported_value,
                        "graph": target.graph,
                        "lmb_channel": target.lmb_channel,
                        "x_peak": list(target.x_peak),
                        "axis_indices": list(target.axis_indices),
                    }
                    for target in targets
                ],
            },
            indent=2,
        )
        + "\n"
    )

    run_approach_scan(
        targets,
        args.output_dir,
        args.state_dir,
        args.points_per_line,
        args.skip_arb_peak_checks,
        args.skip_inspect_spotchecks,
        overrides,
        args.xi_scale,
    )
    if not args.skip_threshold_decomposition:
        run_threshold_decomposition(targets, args.output_dir, args.state_dir, args.xi_scale)
    print(f"Wrote current high-stat max-weight diagnostics to {args.output_dir}", flush=True)


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
