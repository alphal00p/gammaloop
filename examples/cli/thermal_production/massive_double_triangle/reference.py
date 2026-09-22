"""Independent massive scalar thermal double-triangle cuts from hep-ph/0406109, Eq. (58).

No GammaLoop expression is used. Analytically take the retarded discontinuity
and integrate each cut delta function before doing Gauss-Legendre quadrature.
"""

import argparse
import json
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--temperature", type=float, default=0.5)
parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
args = parser.parse_args()
if args.temperature <= 0:
    parser.error(
        "--temperature must be positive; every run also computes the vacuum limit"
    )

M, H, L, T, symmetry = 1.0, 0.6, 0.25, args.temperature, 0.5
results = []
for order in (32, 64, 128, 256, 512):
    x, w = np.polynomial.legendre.leggauss(order)
    u, wu = (x + 1) / 2, w / 2
    # Positive Euclidean HHH triangle, all three external invariants subthreshold.
    a, b = u[:, None], (1 - u[:, None]) * u[None, :]
    c = 1 - a - b
    triangle_vacuum = np.sum(
        wu[:, None]
        * wu[None, :]
        * (1 - a)
        / (H * H - M * M * a * b - L * L * c * (a + b))
    ) / (16 * np.pi**2)
    eq = M / 2
    q = np.sqrt(eq * eq - L * L)
    p = (H * u / (1 - u))[:, None]
    dp = (H / (1 - u) ** 2 * wu)[:, None]
    ep = np.sqrt(p * p + H * H)
    er = np.sqrt(p * p + q * q - 2 * p * q * x[None, :] + H * H)
    np_ = np.exp(-ep / T) / (-np.expm1(-ep / T))
    nr = np.exp(-er / T) / (-np.expm1(-er / T))
    # Thermal part of the virtual B-cut residue: n(si E)-n(si E)|T=0=si n(E).
    summed = np.zeros_like(er)
    summed_vac = np.zeros_like(er)
    for s1 in (-1, 1):
        for s3 in (-1, 1):
            for s5 in (-1, 1):
                A = M + (s3 - s1) * ep
                C = s5 * er + eq - s1 * ep
                D = C - A
                summed += (
                    s1
                    * s3
                    * s5
                    * (s1 * np_ / (A * C) - s3 * np_ / (A * D) + s5 * nr / (C * D))
                )
                summed_vac += (
                    s1
                    * s3
                    * s5
                    * (-(s1 < 0) / (A * C) + (s3 < 0) / (A * D) - (s5 < 0) / (C * D))
                )
    triangle_thermal = np.sum(dp * w[None, :] * p * p / (ep * ep * er) * summed) / (
        32 * np.pi**2
    )
    triangle_vacuum_energy = np.sum(
        dp * w[None, :] * p * p / (ep * ep * er) * summed_vac
    ) / (32 * np.pi**2)
    enhancement = 1 + 2 * np.exp(-eq / T) / (-np.expm1(-eq / T))
    virtual_vacuum = q / eq / (16 * np.pi) * triangle_vacuum
    virtual_thermal = (
        enhancement * q / eq / (16 * np.pi) * (triangle_vacuum + triangle_thermal)
    )
    # Real pole U: M=e5+e2-e3. Sum the two uncut energy signs analytically.
    # Three kinematically allowed crossed cuts; the all-outgoing cut is closed.
    real_parts = []
    for s2, s3, s5 in ((1, 1, 1), (1, -1, -1), (-1, -1, 1)):
        threshold = (
            (M * M + L * L) / (2 * M)
            if s2 == 1
            else (4 * H * H - M * M - L * L) / (2 * M)
        )
        if s2 == 1:
            eq = (threshold + T * u / (1 - u))[:, None]
            deq = (T / (1 - u) ** 2 * wu)[:, None]
        else:
            # Square the coordinate to remove the pair-production endpoint.
            eq = (max(threshold, L) + T * (u / (1 - u)) ** 2)[:, None]
            deq = (2 * T * u / (1 - u) ** 3 * wu)[:, None]
        q = np.sqrt(eq * eq - L * L)
        aq = M - s2 * eq
        t = M * M + L * L - 2 * M * s2 * eq
        spread = q * np.sqrt(1 - 4 * H * H / t)
        if s2 == 1:
            ep_min = (-aq * s3 + spread) / 2
            ep = ep_min + T * u[None, :] / (1 - u[None, :])
            dep = (T / (1 - u) ** 2 * wu)[None, :]
        else:
            ep_min = (aq - spread) / 2
            ep = ep_min + spread * u[None, :]
            dep = spread * wu[None, :]
        er = s5 * (M - s2 * eq + s3 * ep)
        n2p = np.exp(-eq / T) / (-np.expm1(-eq / T))
        n3p = np.exp(-ep / T) / (-np.expm1(-ep / T))
        n5p = np.exp(-er / T) / (-np.expm1(-er / T))
        # Explicit positive gain-minus-loss weights avoid subtracting unit terms.
        if s2 == 1 and s3 == 1:
            numerator = n3p * (1 + n2p + n5p) - n2p * n5p
        elif s2 == 1:
            numerator = n5p * (1 + n2p + n3p) - n2p * n3p
        else:
            numerator = n2p * (1 + n3p + n5p) - n3p * n5p
        integrand = (
            s2
            * s3
            * s5
            * numerator
            / ((M * M + 2 * M * s3 * ep) * (M * M - 2 * M * s2 * eq))
        )
        real_parts.append(float(np.sum(deq * dep * integrand) / (64 * np.pi**3)))
    # Reflection maps U-cut to V-cut because m1=m3 and m2=m4.
    real_thermal = 2 * sum(real_parts)
    row = {
        "order": order,
        "triangle_vacuum": float(triangle_vacuum),
        "triangle_vacuum_energy": float(triangle_vacuum_energy),
        "triangle_thermal": float(triangle_thermal),
        "virtual_ImI_vacuum": float(virtual_vacuum),
        "virtual_ImI_thermal": float(virtual_thermal),
        "real_U_parts_ImI": real_parts,
        "real_ImI_thermal": float(real_thermal),
        "inclusive_ImI_thermal": float(virtual_thermal + real_thermal),
        "Gamma_vacuum": float(symmetry / M * virtual_vacuum),
        "Gamma_virtual_thermal": float(symmetry / M * virtual_thermal),
        "Gamma_real_thermal": float(symmetry / M * real_thermal),
        "Gamma_inclusive_thermal": float(
            symmetry / M * (virtual_thermal + real_thermal)
        ),
        "production_source_thermal": float(
            symmetry / M * (virtual_thermal + real_thermal) / np.expm1(M / T)
        ),
        "real_virtual_cancellation_fraction": float(-real_thermal / virtual_thermal),
    }
    print(json.dumps(row), flush=True)
    results.append(row)
args.output.write_text(
    json.dumps(
        {
            "source": "https://arxiv.org/pdf/hep-ph/0406109, Eqs. (4), (58), (59)",
            "method": "Independent cut quadrature of the published arbitrary-mass integral representation; these benchmark values are not a published numerical table.",
            "conventions": {
                "Pi_R": "-symmetry * I_R",
                "Gamma_rho": "symmetry * Im(I_R) / M",
                "production_source": "n_B(M) * Gamma_rho = df_X/dt at f_X=0",
                "real_U_signs": [[1, 1, 1], [1, -1, -1], [-1, -1, 1]],
                "real_U_sign_order": ["s2", "s3", "s5"],
                "real_U_absorbed_line": [3, 5, 2],
            },
            "parameters": {
                "M": M,
                "H": H,
                "L": L,
                "T": T,
                "symmetry": symmetry,
                "cubic_couplings": 1.0,
                "chemical_potentials": 0.0,
            },
            "convergence": results,
        },
        indent=2,
        allow_nan=False,
    )
    + "\n"
)
