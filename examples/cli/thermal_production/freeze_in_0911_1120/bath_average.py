"""Fold GammaLoop cross sections with the Maxwell-Boltzmann incoming bath.

Run from any directory with Python + NumPy. All native integrations use the
adjacent bath_average.toml; the generated per-node run commands are saved too.
The thermal column deliberately tests the current equilibrium cut, not a
claim that this computes the paper's out-of-equilibrium production observable.
"""

import argparse
import json
import math
import subprocess
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--order", type=int, default=8)
parser.add_argument("--samples", type=int, default=100000)
parser.add_argument("--work-dir", type=Path)
parser.add_argument("--skip-run", action="store_true")
args = parser.parse_args()
here = Path(__file__).resolve().parent
repo = here.parents[3]
work = (args.work_dir or Path(f"/tmp/gammaloop-freeze-in-bath-n{args.order}")).resolve()
work.mkdir(parents=True, exist_ok=True)

# MX = T = lambda = 1. Let a = sqrt(s) = 1 + u, q0 = a + v.
# gamma = 3/(32*pi^4) int ds s int dq0 sqrt(q0^2-s) exp(-q0) sigma(s,q0).
# The factor 3 counts the three distinct incoming bath-species pairs.
# Integrate u with Laguerre(alpha=0) and v with Laguerre(alpha=1/2).
u, wu = np.polynomial.laguerre.laggauss(args.order)
k = np.arange(args.order)
jacobi = np.diag(2 * k + 1.5)
off_diagonal = np.sqrt(k[1:] * (k[1:] + 0.5))
jacobi += np.diag(off_diagonal, 1) + np.diag(off_diagonal, -1)
v, eigenvectors = np.linalg.eigh(jacobi)
wv = math.gamma(1.5) * eigenvectors[0] ** 2
nodes = []
commands = []
# Native integration uses seed + iteration; separate nodes by more than their
# maximum possible iteration count so their random streams cannot overlap.
seed_stride = args.samples // 5000 + 2
for ui, wi in zip(u, wu):
    a = 1 + float(ui)
    s = a * a
    for vj, wj in zip(v, wv):
        q0 = a + float(vj)
        q = math.sqrt(float(vj) * (2 * a + float(vj)))
        e1 = (q0 + q) / 2
        e2 = s / (4 * e1)
        index = len(nodes)
        coefficient = (
            3
            * math.exp(-1)
            / (16 * math.pi**4)
            * wi
            * wj
            * a**3
            * math.sqrt(2 * a + vj)
        )
        sigma = (s - 1) / (16 * math.pi * s * s)
        # Angle-averaged 1 + n_B(E_B3) + n_B(E_X), independently of GammaLoop.
        width = (s - 1) * q / s
        b_min = (s - 1) * (q0 - q) / (2 * s)
        x_min = ((s + 1) * q0 - (s - 1) * q) / (2 * s)
        enhancement = 1.0
        for energy_min in (b_min, x_min):
            enhancement += (
                math.log(-math.expm1(-energy_min - width))
                - math.log(-math.expm1(-energy_min))
            ) / width
        nodes.append(
            {
                "index": index,
                "sqrt_s": a,
                "q0": q0,
                "q": q,
                "coefficient": float(coefficient),
                "vacuum_exact": sigma,
                "equilibrium_cut_exact": sigma * enhancement,
            }
        )
        commands.append(
            f"run sample -D ECM={a:.17g} -D E1={e1:.17g} -D E2={e2:.17g} -D PZ2={-e2:.17g} -D SEED={1337 + seed_stride * index} -D SAMPLES={args.samples} -D WORKSPACE={work}/node_{index:04d}"
        )

(work / "nodes.json").write_text(json.dumps(nodes, indent=2) + "\n")
# This extra block makes the exact executed node sequence directly replayable.
base_card = (here / "bath_average.toml").read_text()
card = work / "executed_run_card.toml"
card.write_text(
    base_card
    + '\n[[command_blocks]]\nname = "bath_nodes"\ncommands = [\n'
    + ",\n".join("  " + json.dumps(c) for c in commands)
    + "\n]\n"
)
if not args.skip_run:
    command = [
        str(repo / "gammaloop"),
        "-s",
        str(work / "state"),
        "-n",
        "-L",
        "off",
        "-l",
        "warn",
        str(card),
        "run",
        "generate",
        "bath_nodes",
    ]
    (work / "invocation.json").write_text(json.dumps(command, indent=2) + "\n")
    print(
        f"Integrating {len(nodes)} bath nodes, {args.samples} samples each; log: {work / 'run.log'}",
        flush=True,
    )
    with (work / "run.log").open("w") as log:
        subprocess.run(
            command, cwd=repo, stdout=log, stderr=subprocess.STDOUT, check=True
        )

# Independent integral representation K1(1) = int_0^infty exp(-cosh(t))*cosh(t) dt.
t, wt = np.polynomial.legendre.leggauss(160)
t = 4 * (t + 1)
k1 = float(np.sum(4 * wt * np.exp(-np.cosh(t)) * np.cosh(t)))
reference = 3 * k1 / (128 * math.pi**5)
summary = {
    "order": args.order,
    "samples_per_node": args.samples,
    "bath_nodes": len(nodes),
    "units": "MX=T=lambda=1; rates in MX^4",
    "paper_equation": "0911.1120 (6.24)",
    "paper_rate": reference,
    "results": {},
}
for process, exact_key in [
    ("vacuum", "vacuum_exact"),
    ("thermal", "equilibrium_cut_exact"),
]:
    total = variance = quadrature_exact = 0.0
    worst_pull = 0.0
    unstable = 0
    for node in nodes:
        result = json.loads(
            (work / f"node_{node['index']:04d}" / "integration_result.json").read_text()
        )
        slot = next(slot for slot in result["slots"] if slot["process"] == process)
        integral = slot["integral"]
        value, error = integral["result"]["re"], integral["error"]["re"]
        weight = node["coefficient"]
        total += weight * value
        variance += (weight * error) ** 2
        quadrature_exact += weight * node[exact_key]
        if error > 0:
            worst_pull = max(worst_pull, abs(value - node[exact_key]) / error)
        unstable += slot["integration_statistics"]["nan_or_unstable_percentage"] > 0
        node[process] = {"result": value, "error": error, "samples": integral["neval"]}
    summary["results"][process] = {
        "rate": total,
        "monte_carlo_error": math.sqrt(variance),
        "analytic_same_quadrature": quadrature_exact,
        "ratio_to_paper": total / reference,
        "worst_node_pull": worst_pull,
        "nodes_with_unstable_samples": unstable,
    }
(work / "node_results.json").write_text(json.dumps(nodes, indent=2) + "\n")
(work / "comparison.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
