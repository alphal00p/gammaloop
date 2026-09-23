"""Independent MSbar spectral quadrature for the scalar self-energy insertion.

Bubble ingredients follow hep-ph/0302098, section III; the raised-line spectral
reduction is derived in reference_derivation.md. No GammaLoop output is read.
"""

import argparse
import json
from itertools import pairwise
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--temperature", type=float, default=0.5)
parser.add_argument("--mu-r", type=float, default=1.6)
parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
args = parser.parse_args()
if args.temperature <= 0 or args.mu_r <= 0:
    parser.error("temperature and mu-r must be positive")

M, H, L, T, mu_r = 1.6, 0.6, 0.25, args.temperature, args.mu_r
convergence = []
for order in (32, 64, 128, 256, 512):
    x, w = np.polynomial.legendre.leggauss(order)
    u, wu = (x + 1) / 2, w / 2
    derivatives = []
    # The bubble is analytic near the stable H shell. Complex-step derivatives
    # change only the squared mass of the repeated line, retaining fixed H,L
    # masses inside the inserted bubble and on the spectator propagator.
    for step in (1e-4, 1e-8, 1e-20):
        a = H * H + 1j * step
        ea = (M * M + a - H * H) / (2 * M)
        eb = M - ea
        q = np.sqrt(ea * ea - a)
        bubble_vacuum = -np.dot(
            wu, np.log((u * H * H + (1 - u) * L * L - u * (1 - u) * a) / mu_r**2)
        ) / (16 * np.pi**2)
        k = T * u / (1 - u)
        dk = T * wu / (1 - u) ** 2
        bubble_thermal = 0j
        for mi, mj in ((H, L), (L, H)):
            energy = np.sqrt(k * k + mi * mi)
            occupation = np.exp(-energy / T) / (-np.expm1(-energy / T))
            delta = mj * mj - mi * mi - a
            logs = np.log(
                (delta + 2 * ea * energy + 2 * q * k)
                / (delta + 2 * ea * energy - 2 * q * k)
            ) + np.log(
                (delta - 2 * ea * energy + 2 * q * k)
                / (delta - 2 * ea * energy - 2 * q * k)
            )
            bubble_thermal += np.sum(dk * k * occupation / energy * logs) / (
                16 * np.pi**2 * q
            )
        phase_space = q / (8 * np.pi * M)
        enhancement = 1 + 1 / np.expm1(ea / T) + 1 / np.expm1(eb / T)
        pole_vacuum = phase_space * bubble_vacuum
        pole_thermal = phase_space * enhancement * (bubble_vacuum + bubble_thermal)
        derivatives.append(
            {
                "complex_mass_squared_step": step,
                "Gamma_HH_vacuum": float(-pole_vacuum.imag / step / M),
                "Gamma_HH_thermal": float(-pole_thermal.imag / step / M),
            }
        )

    # Split at every change of bubble spectral support. A sine-squared map
    # removes square-root endpoint behavior on the finite intervals.
    endpoints = [
        H,
        (M * M + H * H - (H + L) ** 2) / (2 * M),
        (M * M + H * H - (H - L) ** 2) / (2 * M),
        (M * M + H * H) / (2 * M),
        M,
        np.inf,
    ]
    continuum = []
    for lower, upper in pairwise(endpoints):
        if np.isinf(upper):
            energy = lower + T * u / (1 - u)
            de = T * wu / (1 - u) ** 2
        else:
            energy = lower + (upper - lower) * np.sin(np.pi * u / 2) ** 2
            de = (upper - lower) * np.pi / 2 * np.sin(np.pi * u) * wu
        q = np.sqrt(energy * energy - H * H)
        vacuum, thermal = 0.0, 0.0
        for eta in (-1, 1):
            omega = M + eta * energy
            omega_abs = np.abs(omega)
            s = omega * omega - q * q
            a_pair = s + H * H - L * L
            discriminant = a_pair * a_pair - 4 * s * H * H
            # spectral_integral = 16*pi*q*Im B(|omega|,q). The angular
            # delta function leaves energy integrals with elementary Bose
            # antiderivative T*log(1-exp(-E/T)).
            spectral_integral = np.zeros_like(energy)
            spectral_vacuum = np.zeros_like(energy)
            mask = s > (H + L) ** 2
            width = q[mask] * np.sqrt(discriminant[mask]) / s[mask]
            lo = (
                omega_abs[mask] * a_pair[mask] - q[mask] * np.sqrt(discriminant[mask])
            ) / (2 * s[mask])
            hi = lo + width
            spectral_integral[mask] = width + T * (
                np.log(-np.expm1(-hi / T))
                - np.log(-np.expm1(-lo / T))
                + np.log(-np.expm1(-(omega_abs[mask] - lo) / T))
                - np.log(-np.expm1(-(omega_abs[mask] - hi) / T))
            )
            spectral_vacuum[mask] = width
            # Timelike Landau cut: H is heavier than L and H -> Q + L.
            mask = (s > 0) & (s < (H - L) ** 2)
            lo = (
                omega_abs[mask] * a_pair[mask] - q[mask] * np.sqrt(discriminant[mask])
            ) / (2 * s[mask])
            hi = (
                omega_abs[mask] * a_pair[mask] + q[mask] * np.sqrt(discriminant[mask])
            ) / (2 * s[mask])
            spectral_integral[mask] = T * (
                np.log(-np.expm1(-(hi - omega_abs[mask]) / T))
                - np.log(-np.expm1(-(lo - omega_abs[mask]) / T))
                - np.log(-np.expm1(-hi / T))
                + np.log(-np.expm1(-lo / T))
            )
            # Spacelike Landau cut: both energy-sign assignments contribute.
            mask = s < 0
            lo_plus = (
                omega_abs[mask] * a_pair[mask] - q[mask] * np.sqrt(discriminant[mask])
            ) / (2 * s[mask])
            lo_minus = (
                -omega_abs[mask] * a_pair[mask] - q[mask] * np.sqrt(discriminant[mask])
            ) / (2 * s[mask])
            spectral_integral[mask] = T * (
                -np.log(-np.expm1(-(lo_plus - omega_abs[mask]) / T))
                + np.log(-np.expm1(-lo_plus / T))
                - np.log(-np.expm1(-lo_minus / T))
                + np.log(-np.expm1(-(lo_minus + omega_abs[mask]) / T))
            )
            ne = np.exp(-energy / T) / (-np.expm1(-energy / T))
            no_positive = np.exp(-omega_abs / T) / (-np.expm1(-omega_abs / T))
            no = np.where(omega > 0, no_positive, -1 - no_positive)
            weight = 1 + no + ne if eta == -1 else ne - no
            weight_vacuum = (omega > 0).astype(float) if eta == -1 else 0.0
            prefactor = de * np.sign(omega) / (64 * np.pi**3 * (s - H * H) ** 2)
            thermal += float(np.sum(prefactor * spectral_integral * weight))
            vacuum += float(np.sum(prefactor * spectral_vacuum * weight_vacuum))
        continuum.append(
            {
                "spectator_energy_min": lower,
                "spectator_energy_max": None if np.isinf(upper) else upper,
                "Gamma_HHL_vacuum": vacuum / M,
                "Gamma_HHL_thermal": thermal / M,
            }
        )
    hh_vacuum = derivatives[-1]["Gamma_HH_vacuum"]
    hh_thermal = derivatives[-1]["Gamma_HH_thermal"]
    hhl_vacuum = sum(row["Gamma_HHL_vacuum"] for row in continuum)
    hhl_thermal = sum(row["Gamma_HHL_thermal"] for row in continuum)
    row = {
        "order": order,
        "bubble_vacuum_at_H_shell": float(bubble_vacuum.real),
        "bubble_thermal_at_H_shell": float(bubble_thermal.real),
        "complex_step_convergence": derivatives,
        "continuum_intervals": continuum,
        "Gamma_HH_vacuum": hh_vacuum,
        "Gamma_HHL_vacuum": hhl_vacuum,
        "Gamma_inclusive_vacuum": hh_vacuum + hhl_vacuum,
        "Gamma_HH_thermal": hh_thermal,
        "Gamma_HHL_thermal": hhl_thermal,
        "Gamma_inclusive_thermal": hh_thermal + hhl_thermal,
        "production_source_thermal": (hh_thermal + hhl_thermal) / np.expm1(M / T),
    }
    print(json.dumps(row), flush=True)
    convergence.append(row)

args.output.write_text(
    json.dumps(
        {
            "sources": [
                "https://arxiv.org/pdf/hep-ph/0302098, section III (thermal bubble)",
                "https://arxiv.org/pdf/hep-ph/9409388, Eq. (3), Eqs. (13)-(15) (vacuum derivative and dispersion)",
            ],
            "method": "Independent spectral convolution and repeated-mass derivative using published bubble ingredients; these numerical values are not a published table.",
            "conventions": {
                "Pi_R": "-I_R, labelled graph overall factor 1",
                "Gamma_rho": "Im(I_R)/M",
                "production_source": "n_B(M)*Gamma_rho",
                "renormalization": "MSbar at mu_r; vacuum contribution retained",
                "derivative": "Only the squared mass of the repeated line varies; bubble and spectator masses stay fixed.",
            },
            "parameters": {"M": M, "H": H, "L": L, "T": T, "mu_r": mu_r},
            "convergence": convergence,
        },
        indent=2,
        allow_nan=False,
    )
    + "\n"
)
