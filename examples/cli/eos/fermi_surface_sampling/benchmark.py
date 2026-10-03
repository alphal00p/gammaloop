"""Measure thermal vacuum convergence and timing with explicit Fermi channels."""

import argparse
import hashlib
import json
import math
import os
import platform
import random
import resource
import subprocess
import time
import tomllib
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--binary", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--generate", action="store_true")
parser.add_argument(
    "--example", choices=["sunrise", "tennis_ball", "mercedes"], default="sunrise"
)
temperature = parser.add_mutually_exclusive_group()
temperature.add_argument("--betas", type=float, nargs="+", default=[1, 10, 100, 1000])
temperature.add_argument("--zero-temperature", action="store_true")
parser.add_argument("--seeds", type=int, nargs="+", default=[1337, 2843, 4517])
parser.add_argument("--samples", type=int, default=100_000)
parser.add_argument("--iterations", type=int, default=10)
parser.add_argument("--cpu", type=int)
precision = parser.add_mutually_exclusive_group()
precision.add_argument(
    "--arb-only",
    action="store_true",
    help="Use Arb for every sample in a separate precision control",
)
precision.add_argument(
    "--skip-quad",
    action="store_true",
    help="Retain Double and Arb, skipping Quad in a separate precision control",
)
parser.add_argument(
    "--modes",
    nargs="+",
    choices=[
        "ordinary",
        "fermi",
        "mixed",
        "optimized",
        "optimized_fermi",
        "optimized_fermi_repeated",
    ],
)
args = parser.parse_args()
if args.samples <= 0 or args.iterations <= 0 or args.samples % args.iterations:
    parser.error("samples must be positive and divisible by iterations")
if any(not math.isfinite(beta) or beta <= 0 for beta in args.betas):
    parser.error(
        "beta must be finite and positive; use --zero-temperature for exact T=0"
    )
if args.zero_temperature:
    args.betas = [None]

here = Path(__file__).resolve().parent
repo = here.parents[3]
example = here if args.example == "sunrise" else here / "three_loop" / args.example
process = "fermi_sunrise" if args.example == "sunrise" else f"quark_{args.example}"
integrand = "two_loop" if args.example == "sunrise" else "3L"
fermi_channel = "fermi_pair" if args.example == "sunrise" else "fermi_shells"
if args.modes is None:
    args.modes = (
        ["ordinary", "fermi", "mixed", "optimized", "optimized_fermi"]
        if args.example == "sunrise" and not args.zero_temperature
        else ["optimized", "optimized_fermi"]
    )
    if args.zero_temperature and args.example == "tennis_ball":
        args.modes.append("optimized_fermi_repeated")
if args.example != "sunrise" and set(args.modes) - {
    "optimized",
    "optimized_fermi",
    "optimized_fermi_repeated",
}:
    parser.error("three-loop examples require optimized LMB coverage")
if args.example != "tennis_ball" and "optimized_fermi_repeated" in args.modes:
    parser.error("the repeated-momentum diagnostic is specific to tennis_ball")
sampling_card = example / "sampling.toml"
sampling_overlay = sampling_card.read_text() if sampling_card.exists() else ""
binary = args.binary.resolve()
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
state = output / "state"
generation_source = (example / "generate.toml").read_text()
if args.zero_temperature:
    assert generation_source.count('mode = "thermodynamic_equilibrium"') == 1
    generation_source = generation_source.replace(
        'mode = "thermodynamic_equilibrium"', 'mode = "zero_temperature_equilibrium"'
    )
    # This normalized profile localizes T=0 distribution terms independently of
    # the optional importance-sampling Fermi channel.
    generation_source += (
        '\n[default_runtime_settings.h_function]\nfunction = "poly_exponential"\n'
        "sigma = 1.0\npower = 0\n"
    )
stability_overlay = ""
if args.arb_only or args.skip_quad:
    for level in tomllib.loads(generation_source)["default_runtime_settings"][
        "stability"
    ]["levels"]:
        if level["precision"] == "Arb" or (
            args.skip_quad and level["precision"] == "Double"
        ):
            stability_overlay += (
                "\n[[stability.levels]]\n"
                + "\n".join(
                    f"{key} = {json.dumps(value)}" for key, value in level.items()
                )
                + "\n"
            )
prefix = [] if args.cpu is None else ["taskset", "-c", str(args.cpu)]
environment = dict(os.environ, RAYON_NUM_THREADS="1")
revision = subprocess.check_output(
    ["git", "rev-parse", "HEAD"], cwd=repo, text=True
).strip()
with binary.open("rb") as source:
    binary_sha256 = hashlib.file_digest(source, "sha256").hexdigest()
metadata = {
    "example": args.example,
    "zero_temperature": args.zero_temperature,
    "arb_only": args.arb_only,
    "skip_quad": args.skip_quad,
    "revision": revision,
    "binary": str(binary),
    "binary_sha256": binary_sha256,
    "host": platform.node(),
    "cpu_affinity": args.cpu,
    "samples": args.samples,
    "iterations": args.iterations,
    "workers": 1,
    "generation_card_sha256": hashlib.sha256(generation_source.encode()).hexdigest(),
    "sampling_overlay": sampling_overlay,
    "reference": (
        "User-supplied exact-T=0 targets (2026-10-03)"
        if args.zero_temperature
        else "https://arxiv.org/pdf/hep-ph/0305183 (3.15), one massless flavor term"
        if args.example == "sunrise"
        else "User-supplied thermal targets: run cards at beta=1,10; table at beta=100,1000"
    ),
}
if args.generate:
    if (state / "state_manifest.toml").exists():
        parser.error("generation requires a fresh output directory")
    generation_card = output / "generate.toml"
    generation_card.write_text(generation_source)
    started = time.perf_counter()
    with (output / "generate.log").open("w") as log:
        subprocess.run(
            prefix
            + [
                str(binary),
                "-s",
                str(state),
                "-n",
                str(generation_card),
                "quit",
                "-n",
            ],
            cwd=repo,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    metadata["generation_wall_seconds"] = time.perf_counter() - started
    (output / "generation.json").write_text(json.dumps(metadata, indent=2) + "\n")
else:
    generation = json.loads((output / "generation.json").read_text())
    for key in ("revision", "binary_sha256", "generation_card_sha256"):
        if generation[key] != metadata[key]:
            parser.error(
                f"saved state has a different {key}; generate in a fresh output directory"
            )
    metadata["generation_wall_seconds"] = generation["generation_wall_seconds"]

selections = {
    "ordinary": ["lmb(0,1)"],
    "fermi": [fermi_channel],
    "mixed": ["lmb(0,1)", fermi_channel],
    "optimized": ["auto:optimized_lmb"],
    "optimized_fermi": ["auto:optimized_lmb", fermi_channel],
    "optimized_fermi_repeated": ["auto:optimized_lmb", "fermi_repeated"],
}
jobs = [
    (beta, seed, mode)
    for beta in args.betas
    for seed in args.seeds
    for mode in args.modes
]
random.Random(20261002).shuffle(jobs)
records = []
for beta, seed, mode in jobs:
    temperature_label = "T0" if args.zero_temperature else f"beta_{beta:g}"
    run = output / f"{temperature_label}_{mode}_seed_{seed}_n_{args.samples}"
    run.mkdir(exist_ok=False)
    # SU(3), alpha_s=1, mu_q=muB/3=1. The pure-gluon term is absent.
    target = (
        {
            "sunrise": -0.0161257672165997446,
            "tennis_ball": 0.011121480775908030,
            "mercedes": -0.030797946764053006,
        }[args.example]
        if args.zero_temperature
        else -(
            5 * math.pi / (18 * beta**4)
            + 1 / (math.pi * beta**2)
            + 1 / (2 * math.pi**3)
        )
        if args.example == "sunrise"
        else {
            "tennis_ball": {
                1: 0.7113504235805803,
                10: 0.01398321747318332,
                100: 0.011149085893291288,
                1000: 0.011121756598238203,
            },
            "mercedes": {
                1: -1.739558840398423,
                10: -0.04062795126001039,
                100: -0.030915837358603345,
                1000: -0.030799358173484585,
            },
        }[args.example].get(beta)
    )
    overlay = run / "overlay.toml"
    overlay.write_text(
        f"[general]\ninverse_temperature = {1.0 if args.zero_temperature else beta}\n"
        f"[integrator]\nseed = {seed}\nn_start = {args.samples // args.iterations}\n"
        f"n_increase = 0\nn_max = {args.samples}\n"
        f"[sampling]\ndefault_channel_selection = {json.dumps(selections[mode])}\n"
        'power = 2.0\nsampling_channel_weight = "map_density"\n'
        + sampling_overlay
        + stability_overlay
    )
    commands = [
        f"set process -p {process} -i {integrand} file '{overlay}'",
        f"display integrand -p {process} -i {integrand} --show-sampling pretty",
        (
            f"integrate -p {process} -i {integrand} -c 1 --workspace-path '{run / 'workspace'}' "
            + (f"--target 0,{target} " if target is not None else "")
            + "--no-stream-iterations --no-stream-updates --renderer tabled"
        ),
    ]
    card = run / "run.toml"
    card.write_text("commands = " + json.dumps(commands, indent=2) + "\n")
    started = time.perf_counter()
    cpu_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    with (run / "run.log").open("w") as log:
        subprocess.run(
            prefix
            + [
                str(binary),
                "-s",
                str(state),
                "--read-only-state",
                "-n",
                str(card),
                "quit",
                "-n",
            ],
            cwd=repo,
            env=environment,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    wall_seconds = time.perf_counter() - started
    cpu_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = json.loads((run / "workspace" / "integration_result.json").read_text())
    (slot,) = result["slots"]
    integral = slot["integral"]
    assert integral["neval"] == args.samples, integral
    value, error = integral["result"]["im"], integral["error"]["im"]
    assert math.isfinite(value) and math.isfinite(error) and error > 0, integral
    assert abs(integral["result"]["re"]) <= 2 * integral["error"]["re"] + 1e-12, (
        integral
    )
    assert 0 <= integral["error"]["re"] <= 1e-12, integral
    record = {
        "beta": beta,
        "seed": seed,
        "mode": mode,
        "target": target,
        "value": value,
        "error": error,
        "pull": (value - target) / error if target is not None else None,
        "relative_error": error / abs(target) if target is not None else None,
        "wall_seconds": wall_seconds,
        "cpu_seconds": cpu_after.ru_utime
        + cpu_after.ru_stime
        - cpu_before.ru_utime
        - cpu_before.ru_stime,
        "variance_times_wall": error**2 * wall_seconds,
        "result": result,
    }
    (run / "record.json").write_text(json.dumps(record, indent=2) + "\n")
    records.append(record)
    (output / "results.json").write_text(
        json.dumps({"metadata": metadata, "runs": records}, indent=2) + "\n"
    )
    print(
        f"{temperature_label} {mode:16s} seed={seed}: {value:.9g} +/- {error:.3g}; "
        f"{wall_seconds:.2f}s; pull={record['pull']}",
        flush=True,
    )
