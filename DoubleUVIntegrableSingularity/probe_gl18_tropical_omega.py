#!/usr/bin/env python3
"""Probe GL18 tropical target-omega values with short integrations."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "current_maxweight_diagnostics"
GAMMALOOP = ROOT.parent / "gammaloop"
OMEGAS = (2.0, 4.0)


def omega_label(omega: float) -> str:
    return str(omega).replace(".", "p").replace("-", "m")


def card_for_omega(omega: float) -> Path:
    label = omega_label(omega)
    card = Path(f"/tmp/aa_aa_2l_gl18_tropical_omega_{label}.toml")
    state = Path(f"/tmp/gammaloop_aa_aa_2l_gl18_tropical_omega_{label}_state")
    workspace = Path(f"/tmp/aa_aa_2l_gl18_tropical_omega_{label}_workspace")
    body = f"""[cli_settings]
[cli_settings.state]
name = "aa_aa_2L_GL18_tropical_omega_{label}"
folder = "{state}"

[cli_settings.global]
display_directive = "info,gammalooprs[{{#generation, #profile, #uv, #graph, #term, #summary}}]=debug"
logfile_directive = "off"

[cli_settings.global.generation.uv]
softct = false
generate_integrated = true
add_integrated_uv_with_n_loops = 1
subtract_uv = true
only_integrated = false
pole_part = false
add_sigma = false
keep_sigma = true
inner_products = true
use_legacy = true
cached = false
renormalization_schemes = []

[cli_settings.global.generation.evaluator]
iterative_orientation_optimization = false
store_atom = false
compile = false
summed = false
summed_function_map = true
horner_iterations = 5
cpe_iterations = 5

[cli_settings.global.generation.threshold_subtraction]
enable_thresholds = false
check_esurface_at_generation = false

[cli_settings.global.generation.tropical_subgraph_table]
panic_on_fail = true
disable_tropical_generation = false
target_omega = {omega}

[cli_settings.global.n_cores]
feyngen = 10
generate = 10
compile = 10
integrate = 20

[[command_blocks]]
name = "generate_gl18"
commands = [
  "import model sm-default.json",
  "set global kv global.generation.override_lmb_heuristics=true",
  "remove processes -p aa_aa",
  "import graphs ./examples/cli/aa_aa/2L/graphs/GL18.dot -p aa_aa -i 2L -o",
  "generate existing -p aa_aa -i 2L",
  "save state -o",
]

[[command_blocks]]
name = "set_model_parameters"
commands = [
  "set model MT=173.0",
  "set model WT=0.0",
  "set model ymt=173.0",
  "set model aS=0.118",
  "set model aEWM1=137.036",
  "set model Gf=1.166390e-05",
  "set model MZ=91.188",
]

[[command_blocks]]
name = "set_kinematics_a"
commands = [
  \"\"\"set process string '
    [kinematics.externals]
    type = "constant"
    [kinematics.externals.data]
    momenta = [
        [ 86.5, 0.0, 0.0, 86.5 ],
        [ 86.5, 0.0, 0.0, -86.5 ],
        [ 86.5, 0.0, -79.27855952273603, -34.6 ],
        "dependent"
    ]
    '\"\"\",
]

[[command_blocks]]
name = "integrate_gl18_tropical"
commands = [
  "run set_model_parameters",
  "run set_kinematics_a",
  "set process -p aa_aa kv kinematics.externals.data.helicities=[+1,-1,+1,-1] general.mu_r=91.188 general.m_uv=9.188",
  "set process -p aa_aa kv integrator.n_increase=0 integrator.n_start=100000 integrator.n_max=100000 integrator.seed=261101 integrator.integrated_phase=\\"imag\\" integrator.target_relative_accuracy=0.01",
  \"\"\"set process -p aa_aa string '
        [sampling]
        graphs = "monte_carlo"
        orientations = "summed"
        lmb_multichanneling = false
        coordinate_system = "tropical"
    '\"\"\",
  \"\"\"integrate -p aa_aa -i 2L
             --n-cores 20
             --workspace-path {workspace}
             --renderer tabled
             --batch-size 10000
             --show-phase both
             --show-max-weight-info
             --show-top-discrete-grid
             --show-discrete-contributions-sum
             --write-results-for-each-iteration
             --restart
    \"\"\",
]
"""
    card.write_text(body)
    return card


def run_omega(omega: float) -> dict[str, object]:
    label = omega_label(omega)
    card = card_for_omega(omega)
    state = Path(f"/tmp/gammaloop_aa_aa_2l_gl18_tropical_omega_{label}_state")
    workspace = Path(f"/tmp/aa_aa_2l_gl18_tropical_omega_{label}_workspace")
    shutil.rmtree(state, ignore_errors=True)
    shutil.rmtree(workspace, ignore_errors=True)
    log_path = DATA / f"gl18_tropical_omega_{label}_100k.log"
    result_path = DATA / f"gl18_tropical_omega_{label}_100k_result.json"
    DATA.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [
            str(GAMMALOOP),
            "--dev-optim",
            str(card),
            "run",
            "generate_gl18",
            "integrate_gl18_tropical",
        ],
        cwd=ROOT.parent,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_path.write_text(result.stdout)
    if result.returncode != 0:
        raise SystemExit(f"omega={omega} failed with status {result.returncode}; see {log_path}")

    state_result = workspace / "integration_result.json"
    if not state_result.exists():
        raise SystemExit(f"omega={omega} did not produce {state_result}; see {log_path}")
    shutil.copy2(state_result, result_path)
    return json.loads(result_path.read_text())


def main() -> None:
    summaries = []
    for omega in OMEGAS:
        result = run_omega(omega)
        slot = result["slots"][0]
        im_row = next(row for row in slot["table_results"] if row["component"] == "im")
        summaries.append(
            {
                "target_omega": omega,
                "im": slot["integral"]["result"]["im"],
                "im_error": slot["integral"]["error"]["im"],
                "im_relative_error_percent": im_row["relative_error_percent"],
                "im_max_weight_impact": im_row["max_weight_impact"],
                "avg_time": slot["integration_statistics"]["average_total_time_seconds"],
            }
        )
    summary_path = DATA / "gl18_tropical_target_omega_100k_summary.json"
    summary_path.write_text(json.dumps(summaries, indent=2))
    print(summary_path.read_text())


if __name__ == "__main__":
    main()
