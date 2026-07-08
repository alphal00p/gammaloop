#!/usr/bin/env python3
"""Run small import/generate/integrate checks for grouped 3L aa_aa DOT graphs."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
DEFAULT_GRAPHS = (
    ROOT / "examples/cli/aa_aa/3L/graphs_grouped/processes/amplitudes/aa_aa/3L"
)
DEFAULT_BINARY = ROOT / "target/release/gammaloop"
DEFAULT_OUTPUT = ROOT / "examples/cli/aa_aa/3L/grouped_diagram_smokes"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--graphs-dir", type=Path, default=DEFAULT_GRAPHS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument(
        "--only", nargs="*", default=None, help="Graph names such as GL000"
    )
    parser.add_argument("--n-start", type=int, default=20)
    parser.add_argument("--n-increase", type=int, default=0)
    parser.add_argument("--n-max", type=int, default=None)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--cores", type=int, default=20)
    parser.add_argument("--generate-cores", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=1)
    parser.add_argument("--phase", choices=["real", "imag", "both"], default="imag")
    parser.add_argument("--target-relative-accuracy", type=float, default=1.0)
    parser.add_argument("--print-full", action="store_true")
    return parser.parse_args()


def read_status_kib(pid: int) -> tuple[int, int]:
    rss = 0
    hwm = 0
    try:
        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                rss = int(line.split()[1])
            elif line.startswith("VmHWM:"):
                hwm = int(line.split()[1])
    except FileNotFoundError:
        pass
    return rss, hwm


def toml_card(
    dot_path: Path, state_dir: Path, workspace: Path, args: argparse.Namespace
) -> str:
    n_max = args.n_max if args.n_max is not None else args.n_start
    return f'''
commands = ["run smoke"]

[cli_settings]
override_state = false

[cli_settings.state]
name = "aa_aa_3L_grouped_single"
folder = "{state_dir}"

[cli_settings.global]
display_directive = "info"
logfile_directive = "info"

[cli_settings.global.generation.evaluator]
iterative_orientation_optimization = false
store_atom = false
compile = false
summed = false
summed_function_map = false
direct_translation = true
horner_iterations = 5
cpe_iterations = 5

[cli_settings.global.generation.threshold_subtraction]
enable_thresholds = false
check_esurface_at_generation = false

[cli_settings.global.generation.uv]
softct = false
generate_integrated = false
subtract_uv = false
final_integrand = "ThreeD"
pole_part = false
add_marker = false
keep_marker = true
inner_products = true
orchestrator = "legacy_dag_forest"

[cli_settings.global.generation.tropical_subgraph_table]
panic_on_fail = true
disable_tropical_generation = true

[cli_settings.global.n_cores]
feyngen = 1
generate = {args.generate_cores}
compile = 1
integrate = {args.cores}

[[command_blocks]]
name = "smoke"
commands = [
  "import model sm-default.json",
  "import graphs {dot_path} -p aa_aa -i 3L -o",
  "generate existing -p aa_aa -i 3L",
  "run set_model_parameters",
  "run set_kinematics_a",
  "set process -p aa_aa kv kinematics.externals.data.helicities=[+1,-1,+1,-1]",
  "set process -p aa_aa kv general.mu_r=91.188 general.m_uv=9.188 general.evaluator_method=\\"SingleParametric\\"",
  "set process -p aa_aa kv integrator.n_increase={args.n_increase} integrator.n_start={args.n_start} integrator.n_max={n_max} integrator.seed=1337 integrator.integrated_phase=\\"{args.phase}\\" integrator.target_relative_accuracy={args.target_relative_accuracy}",
  """set process -p aa_aa string '
        [sampling]
        graphs = "monte_carlo"
        orientations = "monte_carlo"
        lmb_multichanneling = true
        lmb_channels = "monte_carlo"
        lmb_channel_weight = "inverse_jacobian"
        coordinate_system = "spherical"
        mapping = "linear"
        b = 1.0
    '""",
  """integrate -p aa_aa -i 3L
             --n-cores {args.cores}
             --workspace-path {workspace}
             --renderer tabled
             --batch-size {args.batch_size}
             --show-phase both
             --show-max-weight-info
             --show-top-discrete-grid
             --show-discrete-contributions-sum
             --write-results-for-each-iteration
             --restart
    """,
  "quit",
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
  """set process string '
    [kinematics.externals]
    type = "constant"
    [kinematics.externals.data]
    momenta = [
        [ 86.5, 0.0, 0.0, 86.5 ],
        [ 86.5, 0.0, 0.0, -86.5 ],
        [ 86.5, 0.0, -79.27855952273603, -34.6 ],
        "dependent"
    ]
    '""",
]

[default_runtime_settings.kinematics]
e_cm = 300.0

[default_runtime_settings.kinematics.externals]
type = "constant"

[default_runtime_settings.kinematics.externals.data]
momenta = [
  [86.5, 0.0, 0.0, 86.5],
  [86.5, 0.0, 0.0, -86.5],
  [86.5, 0.0, -79.27855952273603, -34.6],
  "dependent",
]
helicities = [+1, -1, +1, -1]

[default_runtime_settings.general]
evaluator_method = "SingleParametric"
enable_cache = false
debug_cache = false
generate_events = false
store_additional_weights_in_event = false
integral_unit = "none"
mu_r = 91.188
m_uv = 9.188

[default_runtime_settings.subtraction]
disable_threshold_subtraction = true

[default_runtime_settings.sampling]
graphs = "monte_carlo"
orientations = "monte_carlo"
lmb_multichanneling = true
lmb_channels = "monte_carlo"
lmb_channel_weight = "inverse_jacobian"
coordinate_system = "spherical"
mapping = "linear"
b = 1.0

[default_runtime_settings.integrator]
n_start = {args.n_start}
n_increase = {args.n_increase}
n_max = {n_max}
integrated_phase = "{args.phase}"
target_relative_accuracy = {args.target_relative_accuracy}
seed = 1337
'''


ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")
PEAK_RE = re.compile(r"peak RAM\s+([0-9.]+)\s+([A-Za-z]+)")


def summarize_log(text: str) -> dict[str, object]:
    clean = ANSI_RE.sub("", text)
    summary: dict[str, object] = {}
    if "Integrand generation summary" in clean:
        summary["generated"] = True
    if "Integrating using" in clean:
        summary["started_integration"] = True
    if "Finished" in clean or "Completed" in clean:
        summary["has_finished_marker"] = True
    peak = PEAK_RE.search(clean)
    if peak:
        summary["generation_peak_ram"] = f"{peak.group(1)} {peak.group(2)}"
    if "Max weight" in clean or "max weight" in clean:
        summary["has_max_weight_info"] = True
    if re.search(r"\bNaN\b", clean):
        summary["mentions_nan"] = True
    if " ERROR " in clean or "The application panicked" in clean:
        summary["mentions_error"] = True
    return summary


def summarize_result(workspace: Path) -> dict[str, object]:
    result_path = workspace / "integration_result.json"
    if not result_path.exists():
        return {"result_emitted": False}
    data = json.loads(result_path.read_text())
    slot = data["slots"][0]
    integral = slot["integral"]
    stats = slot.get("integration_statistics", {})
    summary: dict[str, object] = {
        "result_emitted": True,
        "neval": integral.get("neval"),
        "result_re": integral["result"]["re"],
        "result_im": integral["result"]["im"],
        "error_re": integral["error"]["re"],
        "error_im": integral["error"]["im"],
        "avg_total_time_s": stats.get("average_total_time_seconds"),
        "avg_integrand_time_s": stats.get("average_integrand_time_seconds"),
        "avg_evaluator_time_s": stats.get("average_evaluator_time_seconds"),
        "nan_or_unstable_percent": stats.get("nan_or_unstable_percentage"),
    }
    for table_result in slot.get("table_results", []):
        component = table_result.get("component")
        if component in {"re", "im"}:
            summary[f"rel_error_{component}_percent"] = table_result.get(
                "relative_error_percent"
            )
            summary[f"max_weight_impact_{component}"] = table_result.get(
                "max_weight_impact"
            )
    for max_weight in slot.get("max_weight_info", []):
        component = max_weight.get("component")
        sign = max_weight.get("sign")
        if component in {"re", "im"} and sign in {"+", "-"}:
            prefix = f"max_weight_{component}_{'plus' if sign == '+' else 'minus'}"
            summary[prefix] = max_weight.get("max_eval")
            summary[f"{prefix}_coordinates"] = max_weight.get("coordinates")
    return summary


def run_one(dot_path: Path, args: argparse.Namespace) -> dict[str, object]:
    graph = dot_path.stem
    work_root = args.output_dir / graph
    state_dir = work_root / "state"
    workspace = work_root / "workspace"
    log_path = work_root / "run.log"
    card_path = work_root / "run.toml"

    if work_root.exists():
        shutil.rmtree(work_root)
    work_root.mkdir(parents=True)
    card_path.write_text(
        toml_card(dot_path.resolve(), state_dir.resolve(), workspace.resolve(), args)
    )

    cmd = [str(args.binary), "--no-save-state", "--clean-state", str(card_path)]
    start = time.monotonic()
    max_rss = 0
    max_hwm = 0
    timed_out = False
    with log_path.open("w") as log_file:
        proc = subprocess.Popen(
            cmd,
            cwd=ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        deadline = start + args.timeout
        while True:
            line = proc.stdout.readline()
            if line:
                log_file.write(line)
                log_file.flush()
            rss, hwm = read_status_kib(proc.pid)
            max_rss = max(max_rss, rss)
            max_hwm = max(max_hwm, hwm)
            if proc.poll() is not None:
                for rest in proc.stdout:
                    log_file.write(rest)
                break
            if time.monotonic() > deadline:
                timed_out = True
                proc.terminate()
                try:
                    proc.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                break
            if not line:
                time.sleep(0.2)

    duration = time.monotonic() - start
    text = log_path.read_text(errors="replace")
    result = {
        "graph": graph,
        "returncode": proc.returncode,
        "timed_out": timed_out,
        "duration_s": round(duration, 3),
        "max_rss_gib": round(max_rss / 1024 / 1024, 3),
        "max_hwm_gib": round(max_hwm / 1024 / 1024, 3),
        "log": str(log_path),
        "workspace": str(workspace),
    }
    result.update(summarize_log(text))
    result.update(summarize_result(workspace))
    return result


def main() -> int:
    args = parse_args()
    dots = sorted(args.graphs_dir.glob("GL*.dot"))
    if args.only:
        requested = set(args.only)
        dots = [dot for dot in dots if dot.stem in requested]
    if args.limit is not None:
        dots = dots[: args.limit]
    if not dots:
        print("No DOT graphs selected.", file=sys.stderr)
        return 2

    args.output_dir.mkdir(parents=True, exist_ok=True)
    jsonl_path = args.output_dir / "summary.jsonl"
    with jsonl_path.open("a") as jsonl:
        for index, dot in enumerate(dots, 1):
            print(f"[{index}/{len(dots)}] {dot.stem}", flush=True)
            result = run_one(dot, args)
            if args.print_full:
                print(json.dumps(result, sort_keys=True), flush=True)
            else:
                print(
                    "  "
                    f"rc={result.get('returncode')} emitted={result.get('result_emitted')} "
                    f"wall={result.get('duration_s')}s "
                    f"rss={result.get('max_hwm_gib')}GiB "
                    f"avg={result.get('avg_total_time_s')}s "
                    f"re={result.get('result_re')} +/- {result.get('error_re')} "
                    f"rel_re={result.get('rel_error_re_percent')}% "
                    f"im={result.get('result_im')} +/- {result.get('error_im')} "
                    f"rel_im={result.get('rel_error_im_percent')}%",
                    flush=True,
                )
            jsonl.write(json.dumps(result, sort_keys=True) + "\n")
            jsonl.flush()
    print(f"Wrote {jsonl_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
