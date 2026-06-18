#!/usr/bin/env python3
"""Scan fixed-bin x-space line approaches through qq_hhh_2L imaginary max-weight points."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from gammaloop import GammaLoopAPI


ROOT = Path(__file__).resolve().parents[3]
EXAMPLE_DIR = ROOT / "examples" / "cli" / "qq_hhh_2L"
RUN_CARD = EXAMPLE_DIR / "qq_hhh_2L.toml"
STATE_DIR = EXAMPLE_DIR / "state"
OUTPUT_DIR = EXAMPLE_DIR / "maxweight_approach_scans"
GAMMALOOP = ROOT / "target" / "dev-optim" / "gammaloop"

LOWER = 1.0e-15
UPPER = 1.0 - 1.0e-15


@dataclass(frozen=True)
class Target:
    label: str
    reported_component: str
    reported_value: float
    graph: int
    lmb_channel: int
    x_peak: tuple[float, float, float, float, float, float]
    axis_indices: tuple[int, int]


@dataclass(frozen=True)
class Direction:
    label: str
    vector: np.ndarray


TARGETS: dict[str, Target] = {
    "imag_plus_global": Target(
        label="imag_plus_global",
        reported_component="im[+]",
        reported_value=3.5795025591653502e6,
        graph=0,
        lmb_channel=11,
        x_peak=(
            0.27543703182272550,
            0.025975374287740258,
            0.50782398548203456,
            0.68942356972268137,
            0.80460957897792340,
            0.52443570508794690,
        ),
        axis_indices=(1, 4),
    ),
    "imag_minus_global": Target(
        label="imag_minus_global",
        reported_component="im[-]",
        reported_value=-7.7414299360872700e5,
        graph=1,
        lmb_channel=9,
        x_peak=(
            0.23262693392408407,
            0.00000034231896575351528,
            0.49485992364488252,
            0.99442214604999712,
            0.33719614224884459,
            0.78136672187917255,
        ),
        axis_indices=(1, 3),
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Evaluate line scans through the qq_hhh_2L imaginary max-weight "
            "points using the GammaLoop Python API."
        )
    )
    parser.add_argument(
        "--points-per-line",
        type=int,
        default=201,
        help="Odd number of samples per line, including the peak point.",
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        choices=sorted(TARGETS),
        default=sorted(TARGETS),
        help="Subset of targets to scan.",
    )
    parser.add_argument(
        "--directions",
        nargs="+",
        default=None,
        help=(
            "Optional subset of direction labels. Labels are axis_x1, axis_x3, "
            "axis_x4, and center_diagonal where available for the chosen target."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="Directory for JSON, TSV, logs, and PNG plots.",
    )
    parser.add_argument(
        "--skip-inspect-spotchecks",
        action="store_true",
        help="Skip CLI inspect peak spot checks.",
    )
    parser.add_argument(
        "--skip-arb-peak-checks",
        action="store_true",
        help="Skip arb-precision Python checks at the exact max-weight points.",
    )
    parser.add_argument(
        "--state-dir",
        type=Path,
        default=STATE_DIR,
        help="GammaLoop state directory to load read-only.",
    )
    parser.add_argument(
        "--run-card",
        type=Path,
        default=RUN_CARD,
        help="Run card providing runtime_common and sampling_lmb_mc blocks.",
    )
    parser.add_argument(
        "--gammaloop",
        type=Path,
        default=GAMMALOOP,
        help="GammaLoop CLI binary for inspect spot checks.",
    )
    parser.add_argument(
        "--evaluation-batch-size",
        type=int,
        default=1,
        help=(
            "Number of points per structured Python API evaluation call. Keep "
            "the default 1 for SingleParametric diagnostics; larger mixed "
            "batches can trip a debug ParamCache assertion."
        ),
    )
    return parser.parse_args()


def directions_for_target(target: Target) -> list[Direction]:
    directions: list[Direction] = []
    for axis_index in target.axis_indices:
        vector = np.zeros(6, dtype=float)
        vector[axis_index] = 1.0
        directions.append(Direction(f"axis_x{axis_index}", vector))

    diagonal = np.array(target.x_peak, dtype=float) - 0.5
    norm = float(np.linalg.norm(diagonal))
    if norm == 0.0:
        raise ValueError(f"Cannot build center diagonal for {target.label}: zero direction")
    directions.append(Direction("center_diagonal", diagonal / norm))
    return directions


def allowed_t_range(x_peak: np.ndarray, direction: np.ndarray) -> tuple[float, float]:
    t_min = -math.inf
    t_max = math.inf
    for x_i, d_i in zip(x_peak, direction, strict=True):
        if d_i > 0.0:
            t_min = max(t_min, (LOWER - x_i) / d_i)
            t_max = min(t_max, (UPPER - x_i) / d_i)
        elif d_i < 0.0:
            t_min = max(t_min, (UPPER - x_i) / d_i)
            t_max = min(t_max, (LOWER - x_i) / d_i)
    if not (t_min < 0.0 < t_max):
        raise ValueError(f"Invalid t range {t_min}, {t_max}")
    return t_min, t_max


def build_line_points(
    target: Target, direction: Direction, points_per_line: int
) -> list[dict[str, Any]]:
    if points_per_line < 3 or points_per_line % 2 == 0:
        raise ValueError("--points-per-line must be an odd integer >= 3")

    x_peak = np.array(target.x_peak, dtype=float)
    t_min, t_max = allowed_t_range(x_peak, direction.vector)
    side_count = (points_per_line - 1) // 2
    fractions = np.geomspace(1.0e-9, 1.0, side_count)

    samples: list[tuple[float, float]] = []
    samples.extend((t_min * fraction, -float(fraction)) for fraction in fractions[::-1])
    samples.append((0.0, 0.0))
    samples.extend((t_max * fraction, float(fraction)) for fraction in fractions)

    rows = []
    for t_value, signed_fraction in samples:
        x = np.clip(x_peak + t_value * direction.vector, LOWER, UPPER)
        rows.append(
            {
                "target": target.label,
                "reported_component": target.reported_component,
                "reported_value": target.reported_value,
                "graph": target.graph,
                "lmb_channel": target.lmb_channel,
                "direction": direction.label,
                "t": float(t_value),
                "signed_side_fraction": signed_fraction,
                "x": [float(value) for value in x],
            }
        )
    return rows


def apply_runtime_settings(api: GammaLoopAPI, xi_scale: float = 3.0) -> None:
    xi_energy = 1000.0 * xi_scale
    api.run("set process -p qq_hhh_2L -i subtracted defaults")
    api.run("set process -p qq_hhh_2L -i subtracted kv subtraction.disable_threshold_subtraction=false")
    api.run(
        """set process -p qq_hhh_2L -i subtracted string '
[kinematics.externals.data]
helicities = [1, -1, 0, 0, 0, 0, 0, 0, 0]
'"""
    )
    api.run(
        f"""set process -p qq_hhh_2L -i subtracted string '
[kinematics.externals.data]
momenta = [
    [5.0000000000000000e+02, 0.0000000000000000e+00, 0.0000000000000000e+00, 5.0000000000000000e+02],
    [5.0000000000000000e+02, 0.0000000000000000e+00, 0.0000000000000000e+00, -5.0000000000000000e+02],
    [{xi_energy:.17e}, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00],
    [{xi_energy:.17e}, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00],
    [{xi_energy:.17e}, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00],
    [{xi_energy:.17e}, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00],
    [4.3855556622469447e+02, 1.5533220018353779e+02, 3.4801603965135871e+02, -1.7737736157184119e+02],
    [3.5636963749219223e+02, -1.6802389008511000e+01, -3.1872911024360047e+02, 9.7487191636880979e+01],
    "dependent"
]
'"""
    )
    api.run(
        """set process -p qq_hhh_2L -i subtracted string '
[general]
generate_events = false
store_additional_weights_in_event = false
evaluator_method = "SingleParametric"
enable_cache = false
debug_cache = false
disable_flux_factor = false
integral_unit = "none"
'"""
    )
    api.run(
        """set process -p qq_hhh_2L -i subtracted string '
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = true
lmb_channels = "monte_carlo"
lmb_channel_weight = "inverse_jacobian"
coordinate_system = "spherical_common_radial"
mapping = "linear"
b = 1.0
'"""
    )
    api.run(
        "set process -p qq_hhh_2L -i subtracted kv "
        "general.generate_events=false "
        "general.store_additional_weights_in_event=false "
        "integrator.integrated_phase=\"imag\""
    )


def stability_summary(stability_results: Any) -> list[dict[str, Any]]:
    if not stability_results:
        return []
    return [
        {
            "precision": result.precision,
            "estimated_relative_accuracy": result.estimated_relative_accuracy,
            "status": result.status,
            "sample_count": result.sample_count,
            "total_time_seconds": result.total_time_seconds,
        }
        for result in stability_results
    ]


def fill_row_from_sample(row: dict[str, Any], sample: Any) -> None:
    value = sample.integrand_result
    raw_re = float(value.real)
    raw_im = float(value.imag)
    jacobian = sample.parameterization_jacobian
    integrator_weight = float(sample.integrator_weight)
    full_factor = (1.0 if jacobian is None else float(jacobian)) * integrator_weight
    re_val = raw_re * full_factor
    im_val = raw_im * full_factor
    row["raw_re"] = raw_re
    row["raw_im"] = raw_im
    row["raw_abs"] = math.hypot(raw_re, raw_im)
    row["re"] = re_val
    row["im"] = im_val
    row["abs"] = math.hypot(re_val, im_val)
    row["parameterization_jacobian"] = sample.parameterization_jacobian
    row["integrator_weight"] = integrator_weight
    row["weight_convention"] = "parameterization_jacobian_x_unit_fixed_bin_sampling"
    row["is_nan"] = sample.is_nan
    row["stability"] = stability_summary(sample.stability_results)
    row["arb_raw_re"] = None
    row["arb_raw_im"] = None
    row["arb_raw_abs"] = None
    row["arb_re"] = None
    row["arb_im"] = None
    row["arb_abs"] = None
    row["arb_rel_delta"] = None


def evaluate_rows(api: GammaLoopAPI, rows: list[dict[str, Any]], batch_size: int) -> None:
    if batch_size < 1:
        raise ValueError("--evaluation-batch-size must be >= 1")
    try:
        if batch_size == 1:
            for index, row in enumerate(rows, start=1):
                result = api.evaluate_sample(
                    row["x"],
                    integrand_name="subtracted",
                    return_events=False,
                    discrete_dim=[int(row["graph"]), int(row["lmb_channel"])],
                )
                fill_row_from_sample(row, result)
                if index % 50 == 0:
                    print(f"evaluated {index}/{len(rows)} points", flush=True)
            return

        for start in range(0, len(rows), batch_size):
            chunk = rows[start : start + batch_size]
            points = np.array([row["x"] for row in chunk], dtype=float)
            discrete_dims = np.array(
                [[int(row["graph"]), int(row["lmb_channel"])] for row in chunk],
                dtype=np.uintp,
            )
            result = api.evaluate_samples(
                points,
                integrand_name="subtracted",
                discrete_dims=discrete_dims,
                return_events=False,
            )
            samples = result.samples
            if len(samples) != len(chunk):
                raise RuntimeError(f"Expected {len(chunk)} samples, got {len(samples)}")
            for row, sample in zip(chunk, samples, strict=True):
                fill_row_from_sample(row, sample)
            print(
                f"evaluated {min(start + batch_size, len(rows))}/{len(rows)} points",
                flush=True,
            )
    except Exception as exc:
        raise RuntimeError(
            "Structured evaluation failed. If this state does not contain the "
            "generated non-compiled qq_hhh_2L integrand, regenerate it with:\n"
            "./target/dev-optim/gammaloop --clean-state "
            "examples/cli/qq_hhh_2L/qq_hhh_2L.toml quit\n"
            "./target/dev-optim/gammaloop -o -s ./examples/cli/qq_hhh_2L/state "
            "run generate_diagrams generate_integrands_single_parametric"
        ) from exc

def add_arb_peak_checks(api: GammaLoopAPI, rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    peak_rows = [row for row in rows if row["direction"] == "axis_x1" and row["t"] == 0.0]
    summaries: dict[str, dict[str, Any]] = {}
    for row in peak_rows:
        result = api.evaluate_sample(
            row["x"],
            integrand_name="subtracted",
            use_arb_prec=True,
            return_events=False,
            discrete_dim=[int(row["graph"]), int(row["lmb_channel"])],
        )
        value = result.integrand_result
        arb_raw_re = float(value.real)
        arb_raw_im = float(value.imag)
        arb_jacobian = result.parameterization_jacobian
        arb_factor = (1.0 if arb_jacobian is None else float(arb_jacobian)) * float(
            result.integrator_weight
        )
        arb_re = arb_raw_re * arb_factor
        arb_im = arb_raw_im * arb_factor
        arb_abs = math.hypot(arb_re, arb_im)
        default_abs = float(row["abs"])
        rel_delta = math.hypot(arb_re - float(row["re"]), arb_im - float(row["im"])) / max(
            default_abs, 1.0e-300
        )
        summaries[str(row["target"])] = {
            "arb_raw_re": arb_raw_re,
            "arb_raw_im": arb_raw_im,
            "arb_raw_abs": math.hypot(arb_raw_re, arb_raw_im),
            "arb_re": arb_re,
            "arb_im": arb_im,
            "arb_abs": arb_abs,
            "arb_rel_delta": rel_delta,
        }

    for row in rows:
        summary = summaries.get(str(row["target"]))
        if summary and row["t"] == 0.0:
            row.update(summary)
    return summaries


def write_json_tsv(rows: list[dict[str, Any]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "imag_maxweight_approach_scan.json"
    tsv_path = output_dir / "imag_maxweight_approach_scan.tsv"
    json_path.write_text(json.dumps(rows, indent=2) + "\n")

    keys = [
        "target",
        "reported_component",
        "reported_value",
        "graph",
        "lmb_channel",
        "direction",
        "t",
        "signed_side_fraction",
        "x0",
        "x1",
        "x2",
        "x3",
        "x4",
        "x5",
        "raw_re",
        "raw_im",
        "raw_abs",
        "re",
        "im",
        "abs",
        "parameterization_jacobian",
        "integrator_weight",
        "weight_convention",
        "is_nan",
        "stability_summary",
        "arb_raw_re",
        "arb_raw_im",
        "arb_raw_abs",
        "arb_re",
        "arb_im",
        "arb_abs",
        "arb_rel_delta",
    ]
    with tsv_path.open("w") as handle:
        handle.write("\t".join(keys) + "\n")
        for row in rows:
            flat = dict(row)
            for index, value in enumerate(row["x"]):
                flat[f"x{index}"] = value
            flat["stability_summary"] = ";".join(
                f"{entry['precision']}:{entry['status']}:{entry['estimated_relative_accuracy']}"
                for entry in row.get("stability", [])
            )
            handle.write(
                "\t".join(
                    "" if flat.get(key) is None else str(flat.get(key, "")) for key in keys
                )
                + "\n"
            )


def plot_target(rows: list[dict[str, Any]], target: Target, output_dir: Path) -> None:
    target_rows = [row for row in rows if row["target"] == target.label]
    if not target_rows:
        return

    max_abs = max(float(row["abs"]) for row in target_rows)
    linthresh_y = max(max_abs * 1.0e-8, 1.0e-30)
    peak_arb = next(
        (
            row
            for row in target_rows
            if row["t"] == 0.0 and row.get("arb_rel_delta") is not None
        ),
        None,
    )
    arb_note = ""
    if peak_arb and float(peak_arb["arb_rel_delta"]) > 1.0e-6:
        arb_note = f" (arb peak rel delta {float(peak_arb['arb_rel_delta']):.2e})"

    fig, (ax_comp, ax_abs) = plt.subplots(2, 1, figsize=(8.5, 7.0), sharex=True)
    for direction in sorted({str(row["direction"]) for row in target_rows}):
        subset = [row for row in target_rows if row["direction"] == direction]
        subset.sort(key=lambda row: float(row["signed_side_fraction"]))
        x_axis = [float(row["signed_side_fraction"]) for row in subset]
        re_vals = [float(row["re"]) for row in subset]
        im_vals = [float(row["im"]) for row in subset]
        abs_vals = [max(float(row["abs"]), 1.0e-300) for row in subset]

        ax_comp.plot(x_axis, re_vals, "-", linewidth=1.1, label=f"{direction} re")
        ax_comp.plot(x_axis, im_vals, "--", linewidth=1.1, label=f"{direction} im")
        ax_abs.plot(x_axis, abs_vals, "-", linewidth=1.2, label=direction)

    for axis in (ax_comp, ax_abs):
        axis.axvline(0.0, color="black", linewidth=0.8, alpha=0.55)
        axis.set_xscale("symlog", linthresh=1.0e-6)
        axis.grid(True, which="both", alpha=0.25)

    ax_comp.set_yscale("symlog", linthresh=linthresh_y)
    ax_abs.set_yscale("log")
    ax_comp.set_ylabel("fixed-bin weight component")
    ax_abs.set_ylabel("|fixed-bin weight|")
    ax_abs.set_xlabel("signed fraction of available x-space line segment")
    ax_comp.legend(ncols=2, fontsize="small")
    ax_abs.legend(fontsize="small")
    fig.suptitle(
        f"{target.label}: graph {target.graph}, LMB channel {target.lmb_channel}, "
        f"reported {target.reported_component} {target.reported_value:.6e}{arb_note}\n"
        "Jacobian-weighted, fixed discrete bin, unit sampling weight"
    )
    fig.tight_layout()
    fig.savefig(output_dir / f"{target.label}_approaches.png", dpi=180)
    plt.close(fig)


def write_inspect_spotchecks(
    selected_targets: list[Target],
    output_dir: Path,
    state_dir: Path,
    gammaloop: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    card_path = output_dir / "imag_maxweight_approach_inspect_spotchecks.toml"
    log_path = output_dir / "imag_maxweight_approach_inspect_spotchecks.log"
    commands = [
        "set process -p qq_hhh_2L -i subtracted defaults",
        "set process -p qq_hhh_2L -i subtracted kv subtraction.disable_threshold_subtraction=false",
        """set process -p qq_hhh_2L -i subtracted string '
[kinematics.externals.data]
helicities = [1, -1, 0, 0, 0, 0, 0, 0, 0]
'""",
        """set process -p qq_hhh_2L -i subtracted string '
[kinematics.externals.data]
momenta = [
    [5.0000000000000000e+02, 0.0000000000000000e+00, 0.0000000000000000e+00, 5.0000000000000000e+02],
    [5.0000000000000000e+02, 0.0000000000000000e+00, 0.0000000000000000e+00, -5.0000000000000000e+02],
    [3.0000000000000000e+03, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00],
    [3.0000000000000000e+03, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00],
    [3.0000000000000000e+03, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00],
    [3.0000000000000000e+03, 0.0000000000000000e+00, 0.0000000000000000e+00, 0.0000000000000000e+00],
    [4.3855556622469447e+02, 1.5533220018353779e+02, 3.4801603965135871e+02, -1.7737736157184119e+02],
    [3.5636963749219223e+02, -1.6802389008511000e+01, -3.1872911024360047e+02, 9.7487191636880979e+01],
    "dependent"
]
'""",
        """set process -p qq_hhh_2L -i subtracted string '
[general]
generate_events = true
store_additional_weights_in_event = true
evaluator_method = "SingleParametric"
enable_cache = false
debug_cache = false
disable_flux_factor = false
integral_unit = "none"
'""",
        """set process -p qq_hhh_2L -i subtracted string '
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = true
lmb_channels = "monte_carlo"
lmb_channel_weight = "inverse_jacobian"
coordinate_system = "spherical_common_radial"
mapping = "linear"
b = 1.0
'""",
        (
            "set process -p qq_hhh_2L -i subtracted kv "
            "general.generate_events=true "
            "general.store_additional_weights_in_event=true "
            "integrator.integrated_phase=\"imag\""
        ),
    ]
    for target in selected_targets:
        commands.append(
            "inspect -p qq_hhh_2L -i subtracted "
            f"--discrete-dim {target.graph} {target.lmb_channel} -f -x "
            + " ".join(f"{value:.17g}" for value in target.x_peak)
        )
    body = ["[[command_blocks]]", 'name = "imag_maxweight_spotchecks"', "commands = ["]
    for command in commands:
        body.append("  " + json.dumps(command) + ",")
    body.append("]")
    card_path.write_text("\n".join(body) + "\n")

    result = subprocess.run(
        [
            str(gammaloop),
            "-n",
            "-s",
            str(state_dir),
            "--read-only-state",
            str(card_path),
            "run",
            "imag_maxweight_spotchecks",
        ],
        cwd=ROOT,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_path.write_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(
            f"GammaLoop inspect spot checks failed with status {result.returncode}; see {log_path}"
        )


def main() -> None:
    args = parse_args()
    selected_targets = [TARGETS[label] for label in args.targets]
    row_groups: list[tuple[Target, Direction, list[dict[str, Any]]]] = []
    for target in selected_targets:
        directions = directions_for_target(target)
        if args.directions:
            requested = set(args.directions)
            directions = [direction for direction in directions if direction.label in requested]
            missing = requested.difference(direction.label for direction in directions)
            if missing:
                raise ValueError(f"Requested directions not available for {target.label}: {missing}")
        for direction in directions:
            row_groups.append(
                (target, direction, build_line_points(target, direction, args.points_per_line))
            )

    rows: list[dict[str, Any]] = []
    for group_index, (target, direction, group_rows) in enumerate(row_groups, start=1):
        print(
            f"evaluating group {group_index}/{len(row_groups)}: "
            f"{target.label} {direction.label}",
            flush=True,
        )
        api = GammaLoopAPI(
            state_folder=args.state_dir,
            read_only_state=True,
        )
        apply_runtime_settings(api)
        evaluate_rows(api, group_rows, args.evaluation_batch_size)
        rows.extend(group_rows)

    if not args.skip_arb_peak_checks:
        api = GammaLoopAPI(
            state_folder=args.state_dir,
            read_only_state=True,
        )
        apply_runtime_settings(api)
        add_arb_peak_checks(api, rows)

    write_json_tsv(rows, args.output_dir)
    for target in selected_targets:
        plot_target(rows, target, args.output_dir)

    if not args.skip_inspect_spotchecks:
        write_inspect_spotchecks(
            selected_targets,
            args.output_dir,
            args.state_dir,
            args.gammaloop,
        )

    print(f"Wrote {len(rows)} scan rows to {args.output_dir}")


if __name__ == "__main__":
    main()
