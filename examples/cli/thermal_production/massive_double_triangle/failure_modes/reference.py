"""Independent finite rates at vacuum and thermal branch points.

Reuse the published Eq. (58) cut reference from ../above_threshold/reference.py.
Only endpoint coordinates, cancellation-safe algebra, and the positive-mass
input guard change. No native GammaLoop result or expression enters this scan.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--masses",
    nargs="+",
    type=float,
    default=[
        0.249,
        0.25,
        0.251,
        0.499,
        0.5,
        0.501,
        0.949,
        0.95,
        0.951,
        1.199,
        1.2,
        1.201,
        1.449,
        1.45,
        1.451,
    ],
)
parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
args = parser.parse_args()
here = Path(__file__).resolve().parent
base = here.parent / "above_threshold" / "reference.py"
source = base.read_text()
# Keep the reference formulas shared: this wrapper adds quadrature charts for
# exact endpoints and narrow thermal boundary layers, rather than a second oracle.
replacements = [
    ("if T <= 0 or M <= 2 * L:", "if T <= 0 or M <= 0:"),
    (
        "use a positive temperature and an external mass above 2 m_L=0.5",
        "use positive temperature and external mass",
    ),
    (
        "for order in (32, 64, 128, 256, 512):",
        "for order in (64, 128, 256, 512, 1024):",
    ),
    ("z = u\n                dz = wu", "z = 1-u*u\n                dz = 2*u*wu"),
    (
        "if disc[0] > 0:",
        "if M == 2*side and left == 0 and right == 1:\n                # Exact branch point: evaluate (z-1)=-u^2 without subtraction.\n                other=4*H*H/(4*ext*ext-M*M)\n                disc=(4*ext*ext-M*M)*(-u*u)*(z-other)\n            if disc[0] > 0:",
    ),
    (
        "e = side + T * (u / (1 - u)) ** 2",
        "energy_shift = T * (u / (1 - u)) ** 2\n        e = side + energy_shift",
    ),
    (
        "p = np.sqrt(e * e - side * side)\n        n =",
        "p = np.sqrt(energy_shift*(2*side+energy_shift))\n        n =",
        1,
    ),
    (
        "log_minus = np.log(np.abs((a0 + M * e + b) / (a0 + M * e - b)))",
        "log_minus = np.log1p(2*b/(a0+M*e-b))",
    ),
    (
        "log_plus = np.log(np.abs((a0 - M * e + b) / (a0 - M * e - b)))",
        "log_plus = np.log1p(2*b/(a0-M*e-b))",
    ),
    (
        "de * n * log_plus / (M * M - 2 * M * e)",
        "de * n * log_plus / (-2*M*(side-M/2+energy_shift))",
    ),
    (
        "e = M / 2 + d + T * u / (1 - u)\n            de = T / (1 - u) ** 2 * wu",
        "e = M / 2 + d + T*(u/(1-u))**2\n            de = 2*T*u/(1-u)**3*wu",
    ),
    (
        "angular = np.log(\n            np.abs((a0 + b - d) * (a0 - b + d) / ((a0 + b + d) * (a0 - b - d)))\n        ) / (2 * d * b)",
        "angular = (np.log1p(2*b/(a0-b-d))-np.log1p(2*b/(a0-b+d)))/(2*d*b)",
    ),
    ("if maximum <= L:", "if M <= 2*H+L:"),
    (
        "f = (threshold + T * u / (1 - u))[:, None]\n            df = (T / (1 - u) ** 2 * wu)[:, None]",
        "f = (threshold + T*(u/(1-u))**2)[:, None]\n            df = (2*T*u/(1-u)**3*wu)[:, None]",
    ),
]
for original, replacement, *count in replacements:
    if source.count(original) != 1:
        raise RuntimeError(
            f"Shared reference changed; expected one occurrence of {original!r}"
        )
    source = source.replace(original, replacement, *count)

results = []
with tempfile.TemporaryDirectory(prefix="gammaloop-threshold-reference-") as work:
    output = Path(work) / "point.json"
    script = Path(work) / "reference_endpoint.py"
    script.write_text(source)
    for mass in args.masses:
        subprocess.run(
            [sys.executable, str(script), "--mass", str(mass), "--output", str(output)],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        data = json.loads(output.read_text())
        result = {
            "M": mass,
            "convergence": [
                {
                    key: row[key]
                    for key in (
                        "order",
                        "Gamma_vacuum",
                        "Gamma_virtual_thermal",
                        "Gamma_real_thermal",
                        "Gamma_inclusive_thermal",
                    )
                }
                for row in data["convergence"]
            ],
        }
        results.append(result)
        print(json.dumps({"M": mass, **result["convergence"][-1]}), flush=True)
args.output.write_text(
    json.dumps(
        {
            "source": "https://arxiv.org/pdf/hep-ph/0406109, Eqs. (4), (58), (59)",
            "method": "Independent numerical evaluation of the published arbitrary-mass integral representation; these rates are not published tabulated numbers or native GammaLoop outputs.",
            "shared_reference": "../above_threshold/reference.py",
            "endpoint_changes": [
                "Allow any positive M; closed two-particle cuts contribute zero, while crossed thermal cuts may remain open.",
                "Square the closed vacuum-triangle endpoint; at M=2H factor its vanishing discriminant directly as z-1=-u^2.",
                "Evaluate side-line momenta, threshold denominators and angular logarithms without subtracting nearby values.",
                "Square virtual-PV tail and crossed-real energy coordinates to resolve narrow boundary layers.",
                "Close the three-particle decay using M<=2H+L, including the exact zero-phase-space endpoint.",
            ],
            "parameters": {
                "H": 0.6,
                "L": 0.25,
                "T": 0.5,
                "symmetry": 0.5,
                "cubic_couplings": 1.0,
                "chemical_potentials": 0.0,
            },
            "normalization": "Pi_R=-symmetry*I_R; Gamma_rho=symmetry*Im(I_R)/M; production source=n_B(M)*Gamma_rho.",
            "stationary_masses": {
                "thermal_H_minus_H_plus_L": 0.25,
                "vacuum_LL": 0.5,
                "thermal_H_plus_H_minus_L": 0.95,
                "vacuum_HH": 1.2,
                "vacuum_HHL": 1.45,
            },
            "results": results,
        },
        indent=2,
        allow_nan=False,
    )
    + "\n"
)
