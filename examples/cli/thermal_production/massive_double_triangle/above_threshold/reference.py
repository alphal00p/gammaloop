"""Independent retarded massive double-triangle cuts, including virtual principal values.

Evaluate the massive integral representation in hep-ph/0406109, Eq. (58).
No GammaLoop expression is used; the discontinuity and angular integrals are
performed analytically before quadrature. The default opens every vacuum cut.
"""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mass", type=float, default=1.6)
parser.add_argument("--temperature", type=float, default=0.5)
parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
args = parser.parse_args()
M, H, L, T, kappa = args.mass, 0.6, 0.25, args.temperature, 0.5
if T <= 0 or M <= 2 * L:
    parser.error("use a positive temperature and an external mass above 2 m_L=0.5")
rows = []
for order in (32, 64, 128, 256, 512):
    x, w = np.polynomial.legendre.leggauss(order)
    u, wu = (x + 1) / 2, w / 2
    triangles = []
    for label, side, ext in [("B", H, L), ("A", L, H)]:
        if M <= 2 * ext:
            continue
        q = np.sqrt(M * M / 4 - ext * ext)
        # Vacuum Feynman parameters u=x+y: integrate the other parameter exactly.
        # The discriminant can change sign; split there and remove its square root.
        roots = np.roots(
            [4 * ext * ext - M * M, 4 * (side * side - H * H - ext * ext), 4 * H * H]
        )
        inside = [
            float(r.real) for r in roots if abs(r.imag) < 1e-12 and 0 < r.real < 1
        ]
        intervals = [0] + sorted(inside) + [1]
        c_vac = 0.0
        c_vac_im = 0.0
        for left, right in itertools.pairwise(intervals):
            if left == 0 and right == 1:
                z = u
                dz = wu
            elif left == 0:
                z = right * (1 - u * u)
                dz = 2 * right * u * wu
            else:
                z = left + (right - left) * u * u
                dz = 2 * (right - left) * u * wu
            a = side * side * z + H * H * (1 - z) - ext * ext * z * (1 - z)
            disc = 4 * a - M * M * z * z
            if inside:
                # Evaluate the vanishing factor directly in the endpoint variable;
                # subtracting z-ucrit loses precision after the square map.
                critical = inside[0]
                other = next(r.real for r in roots if abs(r.real - critical) > 1e-10)
                delta = -right * u * u if left == 0 else (right - left) * u * u
                disc = (4 * ext * ext - M * M) * delta * (z - other)
            if disc[0] > 0:
                f = 4 / (M * np.sqrt(disc)) * np.arctan(M * z / np.sqrt(disc))
                c_vac += np.sum(dz * f) / (16 * np.pi**2)
            else:
                r = np.sqrt(-disc)
                f = 2 / (M * r) * np.log(np.abs((M * z - r) / (M * z + r)))
                c_vac += np.sum(dz * f) / (16 * np.pi**2)
                c_vac_im += np.sum(dz * 2 * np.pi / (M * r)) / (16 * np.pi**2)
        # Thermal triangle: put each internal line on shell, angular-integrate,
        # then take the principal value of the remaining side-line energy pole.
        a0 = ext * ext + side * side - H * H
        e = side + T * (u / (1 - u)) ** 2
        de = 2 * T * u / (1 - u) ** 3 * wu
        p = np.sqrt(e * e - side * side)
        n = np.exp(-e / T) / (-np.expm1(-e / T))
        b = 2 * q * p
        log_minus = np.log(np.abs((a0 + M * e + b) / (a0 + M * e - b)))
        side_minus = np.sum(de * n * log_minus / (M * M + 2 * M * e)) / (
            8 * np.pi**2 * q
        )
        if M / 2 <= side:
            log_plus = np.log(np.abs((a0 - M * e + b) / (a0 - M * e - b)))
            side_plus = np.sum(de * n * log_plus / (M * M - 2 * M * e)) / (
                8 * np.pi**2 * q
            )
        else:
            # Pair equal distances around the pole E=M/2. This is exactly PV,
            # and leaves an ordinary exponentially convergent upper tail.
            d = M / 2 - side
            v = d * (1 - u * u)
            dv = 2 * d * u * wu
            pair = []
            for sign in (1, -1):
                e = M / 2 + sign * v
                p = np.sqrt(e * e - side * side)
                n = np.exp(-e / T) / (-np.expm1(-e / T))
                b = 2 * q * p
                pair.append(n * np.log(np.abs((a0 - M * e + b) / (a0 - M * e - b))))
            pv = np.sum(dv * (pair[0] - pair[1]) / v)
            e = M / 2 + d + T * u / (1 - u)
            de = T / (1 - u) ** 2 * wu
            p = np.sqrt(e * e - side * side)
            n = np.exp(-e / T) / (-np.expm1(-e / T))
            b = 2 * q * p
            f = n * np.log(np.abs((a0 - M * e + b) / (a0 - M * e - b)))
            pv += np.sum(de * f / (e - M / 2))
            side_plus = -pv / (16 * np.pi**2 * M * q)
        # The central-line thermal cut has no threshold pole because each HHL
        # vertex is stable; angular integration leaves a smooth logarithm.
        e = H + T * (u / (1 - u)) ** 2
        de = 2 * T * u / (1 - u) ** 3 * wu
        p = np.sqrt(e * e - H * H)
        n = np.exp(-e / T) / (-np.expm1(-e / T))
        b = 2 * q * p
        a0 = ext * ext + H * H - side * side
        d = M * e
        angular = np.log(
            np.abs((a0 + b - d) * (a0 - b + d) / ((a0 + b + d) * (a0 - b - d)))
        ) / (2 * d * b)
        central = np.sum(de * p * n * angular) / (4 * np.pi**2)
        c_thermal = side_plus + side_minus + central
        enhancement = 1 + 2 * np.exp(-M / (2 * T)) / (-np.expm1(-M / (2 * T)))
        im_vac = 2 * q / M * c_vac / (16 * np.pi)
        im_thermal = enhancement * 2 * q / M * (c_vac + c_thermal) / (16 * np.pi)
        triangles.append(
            {
                "cut": label,
                "triangle_vacuum_real": float(c_vac),
                "triangle_vacuum_imag": float(c_vac_im),
                "triangle_thermal": float(c_thermal),
                "triangle_thermal_side_plus": float(side_plus),
                "triangle_thermal_side_minus": float(side_minus),
                "triangle_thermal_central": float(central),
                "ImI_vacuum": float(im_vac),
                "ImI_thermal": float(im_thermal),
                "Gamma_vacuum": float(kappa / M * im_vac),
                "Gamma_thermal": float(kappa / M * im_thermal),
            }
        )
    real_parts = []
    real_vac = 0.0
    for s2, s3, s5 in ((1, 1, 1), (1, -1, -1), (-1, -1, 1), (1, -1, 1)):
        decay = s2 == 1 and s3 == -1 and s5 == 1
        if decay:
            maximum = (M * M + L * L - 4 * H * H) / (2 * M)
            if maximum <= L:
                continue
            f = (L + (maximum - L) * np.sin(np.pi * u / 2) ** 2)[:, None]
            df = ((maximum - L) * np.pi / 2 * np.sin(np.pi * u) * wu)[:, None]
        elif s2 == 1:
            threshold = (M * M + L * L) / (2 * M)
            f = (threshold + T * u / (1 - u))[:, None]
            df = (T / (1 - u) ** 2 * wu)[:, None]
        else:
            threshold = max(L, (4 * H * H - M * M - L * L) / (2 * M))
            f = (threshold + T * (u / (1 - u)) ** 2)[:, None]
            df = (2 * T * u / (1 - u) ** 3 * wu)[:, None]
        q = np.sqrt(f * f - L * L)
        aq = M - s2 * f
        t = M * M + L * L - 2 * M * s2 * f
        spread = q * np.sqrt(1 - 4 * H * H / t)
        if s2 == 1 and not decay:
            ep_min = (-aq * s3 + spread) / 2
            e = ep_min + T * u[None, :] / (1 - u[None, :])
            de = (T / (1 - u) ** 2 * wu)[None, :]
        else:
            ep_min = (aq - spread) / 2
            e = ep_min + spread * u[None, :]
            de = spread * wu[None, :]
        r = s5 * (M - s2 * f + s3 * e)
        n2 = np.exp(-f / T) / (-np.expm1(-f / T))
        n3 = np.exp(-e / T) / (-np.expm1(-e / T))
        n5 = np.exp(-r / T) / (-np.expm1(-r / T))
        if decay:
            numerator = -(1 + n2 + n3 + n5 + n2 * n3 + n2 * n5 + n3 * n5)
        elif s2 == 1 and s3 == 1:
            numerator = n3 * (1 + n2 + n5) - n2 * n5
        elif s2 == 1:
            numerator = n5 * (1 + n2 + n3) - n2 * n3
        else:
            numerator = n2 * (1 + n3 + n5) - n3 * n5
        factor = s2 * s3 * s5 / ((M * M + 2 * M * s3 * e) * (M * M - 2 * M * s2 * f))
        value = float(np.sum(df * de * factor * numerator) / (64 * np.pi**3))
        real_parts.append({"signs": [s2, s3, s5], "decay": decay, "ImI_U": value})
        if decay:
            real_vac = float(-np.sum(df * de * factor) / (64 * np.pi**3))
    v_vac = sum(t["ImI_vacuum"] for t in triangles)
    v_thermal = sum(t["ImI_thermal"] for t in triangles)
    r_vac = 2 * real_vac
    r_thermal = 2 * sum(t["ImI_U"] for t in real_parts)
    row = {
        "order": order,
        "triangles": triangles,
        "real_U_parts": real_parts,
        "Gamma_virtual_vacuum": kappa / M * v_vac,
        "Gamma_real_vacuum": kappa / M * r_vac,
        "Gamma_vacuum": kappa / M * (v_vac + r_vac),
        "Gamma_virtual_thermal": kappa / M * v_thermal,
        "Gamma_real_thermal": kappa / M * r_thermal,
        "Gamma_inclusive_thermal": kappa / M * (v_thermal + r_thermal),
        "production_source_thermal": float(
            kappa / M * (v_thermal + r_thermal) / np.expm1(M / T)
        ),
    }
    rows.append(row)
    print(json.dumps(row), flush=True)
args.output.write_text(
    json.dumps(
        {
            "source": "https://arxiv.org/pdf/hep-ph/0406109, Eqs. (4), (58), (59)",
            "vacuum_cross_check_source": "https://arxiv.org/pdf/hep-ph/9501201, Eqs. (4)-(12)",
            "method": "Independent evaluation of the published arbitrary-mass integral representation, with analytic cut deltas and angular integrals, explicit virtual principal values, and Feynman-parameter vacuum triangles; these masses are not a published numerical table.",
            "conventions": {
                "Pi_R": "-symmetry * I_R",
                "Gamma_rho": "symmetry * Im(I_R) / M",
                "production_source": "n_B(M) * Gamma_rho = df_X/dt at f_X=0",
                "A_cut": "HH virtual; remaining triangle has internal masses LLH",
                "B_cut": "LL virtual; remaining triangle has internal masses HHH",
                "real_U_sign_order": ["s2", "s3", "s5"],
                "real_reflection": "I_V = I_U",
            },
            "parameters": {
                "M": M,
                "H": H,
                "L": L,
                "T": T,
                "symmetry": kappa,
                "cubic_couplings": 1.0,
                "chemical_potentials": 0.0,
            },
            "convergence": rows,
        },
        indent=2,
        allow_nan=False,
    )
    + "\n"
)
