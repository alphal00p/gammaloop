#!/usr/bin/env python3
"""Low-stat xi scan for qq_hhh_2L variance and max-weight diagnostics."""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
EXAMPLE_DIR = ROOT / "examples" / "cli" / "qq_hhh_2L"
STATE_DIR = EXAMPLE_DIR / "state"
GAMMALOOP = ROOT / "target" / "dev-optim" / "gammaloop"
OUTPUT_DIR = EXAMPLE_DIR / "maxweight_approach_scans" / "xi_lowstat_scan"
WORKSPACE_ROOT = EXAMPLE_DIR / "workspaces" / "xi_lowstat_scan"
COORD_RE = re.compile(
    r"graph:\s*(?P<graph>\d+),\s*LMB channel:\s*(?P<channel>\d+),\s*xs:\s*\[(?P<xs>.*?)\]",
    re.S,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run one low-stat integration per xi=t*(p1+p2) and summarize "
            "variance/max-weight diagnostics."
        )
    )
    parser.add_argument("--state-dir", type=Path, default=STATE_DIR)
    parser.add_argument("--gammaloop", type=Path, default=GAMMALOOP)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--workspace-root", type=Path, default=WORKSPACE_ROOT)
    parser.add_argument("--points", type=int, default=50_000)
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--seed", type=int, default=271000)
    parser.add_argument("--t-values", type=float, nargs="+", default=[0.5, 1.0, 2.0, 3.0])
    parser.add_argument("--sliver-width", type=float, default=None)
    parser.add_argument("--gaussian-width", type=float, default=None)
    parser.add_argument("--h-sigma", type=float, default=None)
    return parser.parse_args()


def format_momenta_block(xi_scale: float) -> str:
    xi_energy = 1000.0 * xi_scale
    return f"""set process -p qq_hhh_2L -i subtracted string '
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


def build_commands(
    xi_scale: float,
    points: int,
    seed: int,
    workspace: Path,
    batch_size: int,
    cores: int,
    sliver_width: float | None,
    gaussian_width: float | None,
    h_sigma: float | None,
) -> list[str]:
    commands = [
        "set process -p qq_hhh_2L -i subtracted defaults",
        "set process -p qq_hhh_2L -i subtracted kv subtraction.disable_threshold_subtraction=false",
        """set process -p qq_hhh_2L -i subtracted string '
[kinematics.externals.data]
helicities = [1, -1, 0, 0, 0, 0, 0, 0, 0]
'""",
        format_momenta_block(xi_scale),
        """set process -p qq_hhh_2L -i subtracted string '
[general]
generate_events = false
store_additional_weights_in_event = false
evaluator_method = "SingleParametric"
enable_cache = true
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
    ]
    assignments = [
        f"integrator.n_start={points}",
        "integrator.n_increase=0",
        f"integrator.n_max={points}",
        f"integrator.seed={seed}",
        'integrator.integrated_phase="imag"',
        "integrator.target_relative_accuracy=1.0e-1",
    ]
    if sliver_width is not None:
        assignments.append(
            f"subtraction.local_ct_settings.uv_localisation.sliver_width={sliver_width:.17g}"
        )
    if gaussian_width is not None:
        assignments.append(
            f"subtraction.local_ct_settings.uv_localisation.gaussian_width={gaussian_width:.17g}"
        )
    commands.append("set process -p qq_hhh_2L -i subtracted kv " + " ".join(assignments))
    if h_sigma is not None:
        commands.append(
            f"""set process -p qq_hhh_2L -i subtracted string '
[subtraction.integrated_ct_settings.range]
type = "infinite"

[subtraction.integrated_ct_settings.range.h_function_settings]
sigma = {h_sigma:.17g}
'"""
        )
    commands.append(
        f"""integrate -p qq_hhh_2L -i subtracted --n-cores {cores} --batch-size {batch_size}
            --workspace-path {workspace}
            --renderer tabled
            --show-phase both
            --show-max-weight-info
            --show-max-weight-info-for-discrete-bins
            --show-top-discrete-grid
            --show-discrete-contributions-sum
            --write-results-for-each-iteration
            --restart
        """
    )
    return commands


def parse_coordinates(coordinates: str | None) -> dict[str, Any] | None:
    if coordinates is None:
        return None
    match = COORD_RE.search(coordinates)
    if not match:
        return {"raw": coordinates}
    xs = [float(value) for value in match.group("xs").split()]
    return {
        "graph": int(match.group("graph")),
        "lmb_channel": int(match.group("channel")),
        "x": xs,
        "raw": coordinates,
    }


def slot_to_summary(xi_scale: float, workspace: Path, slot: dict[str, Any]) -> dict[str, Any]:
    integral = slot["integral"]
    stats = slot["integration_statistics"]
    table_results = []
    for result in slot["table_results"]:
        table_results.append(
            {
                "component": result["component"],
                "value": result["value"],
                "error": result["error"],
                "relative_error_percent": result["relative_error_percent"],
                "chi_sq_per_dof": result["chi_sq_per_dof"],
                "max_weight_impact": result["max_weight_impact"],
            }
        )
    max_weights = []
    for entry in slot["max_weight_info"]:
        max_weights.append(
            {
                "component": entry["component"],
                "sign": entry["sign"],
                "max_eval": entry["max_eval"],
                "coordinates": parse_coordinates(entry["coordinates"]),
            }
        )
    return {
        "xi_scale": xi_scale,
        "xi": [1000.0 * xi_scale, 0.0, 0.0, 0.0],
        "workspace": str(workspace),
        "neval": integral["neval"],
        "result": {"re": integral["result"]["re"], "im": integral["result"]["im"]},
        "error": {"re": integral["error"]["re"], "im": integral["error"]["im"]},
        "relative_error_percent": {
            "re": abs(integral["error"]["re"] / integral["result"]["re"]) * 100.0
            if integral["result"]["re"] != 0.0
            else math.inf,
            "im": abs(integral["error"]["im"] / integral["result"]["im"]) * 100.0
            if integral["result"]["im"] != 0.0
            else math.inf,
        },
        "chisq": {"re": integral["real_chisq"], "im": integral["im_chisq"]},
        "table_results": table_results,
        "integration_statistics": {
            "num_evals": stats["num_evals"],
            "average_total_time_seconds": stats["average_total_time_seconds"],
            "f64_percentage": stats["f64_percentage"],
            "f128_percentage": stats["f128_percentage"],
            "arb_percentage": stats["arb_percentage"],
            "nan_percentage": stats["nan_percentage"],
            "nan_or_unstable_percentage": stats["nan_or_unstable_percentage"],
        },
        "max_weight_info": max_weights,
    }


def write_outputs(summaries: list[dict[str, Any]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "xi_lowstat_summary.json").write_text(json.dumps(summaries, indent=2) + "\n")

    keys = [
        "xi_scale",
        "sliver_width",
        "gaussian_width",
        "h_sigma",
        "neval",
        "re",
        "re_error",
        "re_relative_error_percent",
        "im",
        "im_error",
        "im_relative_error_percent",
        "re_max_weight_impact",
        "im_max_weight_impact",
        "max_re_plus",
        "max_re_minus",
        "max_im_plus",
        "max_im_minus",
        "f128_percent",
        "arb_percent",
        "nan_or_unstable_percent",
        "workspace",
    ]
    with (output_dir / "xi_lowstat_summary.tsv").open("w") as handle:
        handle.write("\t".join(keys) + "\n")
        for summary in summaries:
            table = {row["component"]: row for row in summary["table_results"]}
            maxima = {
                f"{entry['component']}_{'plus' if entry['sign'] == '+' else 'minus'}": entry[
                    "max_eval"
                ]
                for entry in summary["max_weight_info"]
            }
            row = {
                "xi_scale": summary["xi_scale"],
                "sliver_width": summary["runtime_overrides"]["sliver_width"],
                "gaussian_width": summary["runtime_overrides"]["gaussian_width"],
                "h_sigma": summary["runtime_overrides"]["h_sigma"],
                "neval": summary["neval"],
                "re": summary["result"]["re"],
                "re_error": summary["error"]["re"],
                "re_relative_error_percent": summary["relative_error_percent"]["re"],
                "im": summary["result"]["im"],
                "im_error": summary["error"]["im"],
                "im_relative_error_percent": summary["relative_error_percent"]["im"],
                "re_max_weight_impact": table.get("re", {}).get("max_weight_impact"),
                "im_max_weight_impact": table.get("im", {}).get("max_weight_impact"),
                "max_re_plus": maxima.get("re_plus"),
                "max_re_minus": maxima.get("re_minus"),
                "max_im_plus": maxima.get("im_plus"),
                "max_im_minus": maxima.get("im_minus"),
                "f128_percent": summary["integration_statistics"]["f128_percentage"],
                "arb_percent": summary["integration_statistics"]["arb_percentage"],
                "nan_or_unstable_percent": summary["integration_statistics"][
                    "nan_or_unstable_percentage"
                ],
                "workspace": summary["workspace"],
            }
            handle.write("\t".join("" if row.get(key) is None else str(row[key]) for key in keys) + "\n")


def write_command_card(card_path: Path, commands: list[str]) -> None:
    body = ["[[command_blocks]]", 'name = "xi_lowstat"', "commands = ["]
    for command in commands:
        body.append("  " + json.dumps(command) + ",")
    body.append("]")
    card_path.write_text("\n".join(body) + "\n")


def run_cli(args: argparse.Namespace, card_path: Path, log_path: Path) -> None:
    result = subprocess.run(
        [
            str(args.gammaloop),
            "-n",
            "-s",
            str(args.state_dir),
            "--read-only-state",
            str(card_path),
            "run",
            "xi_lowstat",
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
            f"GammaLoop xi low-stat command failed with status {result.returncode}; see {log_path}"
        )


def load_integration_slot(workspace: Path) -> dict[str, Any]:
    result_path = workspace / "results" / "integration_result_iter_0001.json"
    data = json.loads(result_path.read_text())
    slots = data.get("slots", [])
    if len(slots) != 1:
        raise RuntimeError(f"Expected one integration slot in {result_path}, got {len(slots)}")
    return slots[0]


def run_scan(args: argparse.Namespace) -> list[dict[str, Any]]:
    summaries = []
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for index, xi_scale in enumerate(args.t_values, start=1):
        workspace = args.workspace_root / f"xi_t{xi_scale:g}_{args.points}"
        card_path = args.output_dir / f"xi_t{xi_scale:g}_{args.points}.toml"
        log_path = args.output_dir / f"xi_t{xi_scale:g}_{args.points}.log"
        print(
            f"[{index}/{len(args.t_values)}] integrating xi=(p1+p2)*{xi_scale:g} "
            f"with {args.points} samples, workspace {workspace}",
            flush=True,
        )
        commands = build_commands(
            xi_scale=xi_scale,
            points=args.points,
            seed=args.seed + index - 1,
            workspace=workspace,
            batch_size=args.batch_size,
            cores=args.cores,
            sliver_width=args.sliver_width,
            gaussian_width=args.gaussian_width,
            h_sigma=args.h_sigma,
        )
        write_command_card(card_path, commands)
        run_cli(args, card_path, log_path)
        summary = slot_to_summary(xi_scale, workspace, load_integration_slot(workspace))
        summary["runtime_overrides"] = {
            "sliver_width": args.sliver_width,
            "gaussian_width": args.gaussian_width,
            "h_sigma": args.h_sigma,
        }
        summaries.append(summary)
        write_outputs(summaries, args.output_dir)
    return summaries


def main() -> int:
    args = parse_args()
    summaries = run_scan(args)
    write_outputs(summaries, args.output_dir)
    print(f"Wrote xi scan summary to {args.output_dir}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
