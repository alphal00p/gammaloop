"""Vacuum MSbar insertion reference from hep-ph/9409388, Eqs. (3), (13), (15)."""

import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

# Only the serial propagator mass is differentiated; the bubble masses stay fixed.
M, H, L, mu_r = 1.6, 0.6, 0.25, 1.6
s, a, b = M**2, H**2, H**2
root_lambda = np.sqrt((s - a - b) ** 2 - 4 * a * b)
phi = root_lambda / (16 * np.pi * s)
phi_prime = (a - s - b) / (16 * np.pi * s * root_lambda)
t_min, t_max = (H + L) ** 2, (M - H) ** 2
rows = []
for order in [16, 32, 64, 128, 256]:
    nodes, weights = leggauss(order)
    x, wx = (nodes + 1) / 2, weights / 2
    delta = x * H**2 + (1 - x) * L**2 - x * (1 - x) * a
    bubble = np.dot(wx, -np.log(delta / mu_r**2)) / (16 * np.pi**2)
    bubble_prime = np.dot(wx, x * (1 - x) / delta) / (16 * np.pi**2)
    virtual = -(phi_prime * bubble + phi * bubble_prime) / M

    # A sine-squared map removes both three-body square-root endpoints.
    angle, wa = (nodes + 1) * np.pi / 4, weights * np.pi / 4
    t = t_min + (t_max - t_min) * np.sin(angle) ** 2
    dt = 2 * (t_max - t_min) * np.sin(angle) * np.cos(angle)
    spectral_bubble = np.sqrt((t - t_min) * (t - (H - L) ** 2)) / (16 * np.pi**2 * t)
    outer_phase = np.sqrt((s - t - b) ** 2 - 4 * t * b) / (16 * np.pi * s)
    real = np.dot(wa, dt * spectral_bubble * outer_phase / (t - a) ** 2) / M
    rows.append(
        {
            "order": order,
            "bubble_msbar": float(bubble),
            "bubble_derivative": float(bubble_prime),
            "gamma_virtual_hh": float(virtual),
            "gamma_real_hhl": float(real),
            "gamma_inclusive": float(virtual + real),
        }
    )

# The two-particle term also follows from the nondegenerate divided difference.
# This separately verifies the derivative sign and the limit of two simple cuts.
split_rows = []
for step in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    shifted_masses = a + np.array([-step, step])
    shifted_delta = (
        x[:, None] * H**2
        + (1 - x[:, None]) * L**2
        - x[:, None] * (1 - x[:, None]) * shifted_masses
    )
    shifted_bubbles = wx @ -np.log(shifted_delta / mu_r**2) / (16 * np.pi**2)
    shifted_phases = np.sqrt((s - shifted_masses - b) ** 2 - 4 * shifted_masses * b) / (
        16 * np.pi * s
    )
    two_simple_cuts = shifted_phases * shifted_bubbles
    split_rows.append(
        {
            "mass_squared_step": step,
            "gamma_virtual_hh": float(
                (two_simple_cuts[0] - two_simple_cuts[1]) / (2 * step * M)
            ),
        }
    )

result = {
    "source": "https://arxiv.org/pdf/hep-ph/9409388",
    "source_equations": [3, 13, 15],
    "reference_kind": "Independent numerical evaluation of published integral relations; not a published numerical table.",
    "parameters": {"M": M, "H": H, "L": L, "mu_r": mu_r, "T": 0.0},
    "normalization": "Labelled graph factor 1, unit cubic couplings; Gamma_rho = Im I_R/M = -Im Pi_R/M for positive Euclidean five-propagator master I.",
    "renormalization": "MSbar bubble B_R(a) = -int_0^1 dx log((x H^2 +(1-x) L^2 -x(1-x)a)/mu_r^2)/(16 pi^2).",
    "formulas": {
        "phi": "sqrt(lambda(M^2,a,H^2))/(16 pi M^2)",
        "virtual": "-d/da [phi(a,H^2) B_R(a)]/M, at a=H^2; inner bubble masses held fixed",
        "real": "int_(H+L)^2^((M-H)^2) dt sqrt(lambda(t,H^2,L^2))/(16 pi^2 t) phi(t,H^2)/((t-H^2)^2 M)",
    },
    "quadrature": rows,
    "virtual_divided_difference_check": split_rows,
    "gamma_derivative_wrt_log_mu_r_squared": float(-phi_prime / (16 * np.pi**2 * M)),
    "mpmath_50_digit_check": {
        "gamma_virtual_hh": "0.000061180533582880222503890749715318155803584096",
        "gamma_real_hhl": "0.00000793099916099000896111455226614401311312677752",
        "gamma_inclusive": "0.0000691115327438702314650053019814621689167108736",
    },
}
output = Path(__file__).with_suffix(".json")
output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(rows[-1], indent=2))
