"""Compare simple-cut mass limits with the exact raised self-energy graph.

Only the second serial H mass changes; bubble and spectator masses stay fixed.
All comparisons are bare pointwise integrands, not UV-renormalized rates.
"""

import argparse
import json
import re
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--work", type=Path, default=Path("/tmp/gammaloop-self-energy-scan")
)
parser.add_argument(
    "--output", type=Path, default=Path(__file__).with_name("mass_limit.json")
)
args = parser.parse_args()
here = Path(__file__).resolve().parent
root = here.parents[3]
args.work.mkdir(parents=True, exist_ok=True)
source = (here / "split_mass.dot").read_text()
rows = []
for step in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    for sign in [-1, 1]:
        mass = 0.6 + sign * step
        graph = args.work / f"split_{mass:.6f}.dot"
        graph.write_text(source.replace('mass="0.601"', f'mass="{mass:.17g}"'))
        for mode in ["vacuum", "thermodynamic_equilibrium"]:
            prefix = args.work / f"{mode}_{mass:.6f}"
            command = [
                str(root / "gammaloop"),
                "--read-only-state",
                "-s",
                str(prefix) + "-state",
                "-n",
                "-L",
                "off",
                str(here / "split_mass.toml"),
                "run",
                "generate",
                "inspect",
                "-D",
                "mode=" + mode,
                "-D",
                "graph=" + str(graph),
                "-D",
                "output=" + str(prefix),
            ]
            log_path = Path(str(prefix) + ".log")
            with log_path.open("w") as log:
                process = subprocess.run(
                    command, cwd=root, stdout=log, stderr=subprocess.STDOUT, check=False
                )
            row = {
                "step": step,
                "sign": sign,
                "split_mass": mass,
                "mode": mode,
                "exit_code": process.returncode,
            }
            if process.returncode:
                row["error"] = re.sub(r"\x1b\[[0-9;]*m", "", log_path.read_text())
            else:
                evaluation = json.loads(
                    Path(str(prefix) + "-inspect.json").read_text()
                )["evaluation"]
                row["integrand"] = evaluation["integrand_result"]
                row["metadata"] = evaluation["evaluation_metadata"]
            rows.append(row)
            print(
                json.dumps({k: v for k, v in row.items() if k != "metadata"}),
                flush=True,
            )

limits = []
for mode in ["vacuum", "thermodynamic_equilibrium"]:
    for step in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
        pair = [row for row in rows if row["mode"] == mode and row["step"] == step]
        if all(row["exit_code"] == 0 for row in pair):
            limits.append(
                {
                    "mode": mode,
                    "step": step,
                    "symmetric_mean": sum(row["integrand"]["re"] for row in pair) / 2,
                }
            )
args.output.write_text(
    json.dumps(
        {
            "method": "Symmetric split-mass limit of the same graph with simple cut poles; this is a native consistency oracle, not an independent integrated literature rate.",
            "parameters": {"M": 1.6, "H": 0.6, "L": 0.25, "T": 0.5},
            "point": [0.4, 0.2, 0.1, -0.3, 0.5, 0.2],
            "precision": "Arb via inspect -f",
            "rows": rows,
            "limits": limits,
        },
        indent=2,
        allow_nan=False,
    )
    + "\n"
)
