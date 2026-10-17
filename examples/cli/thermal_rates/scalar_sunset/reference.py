#!/usr/bin/env python3
"""Independent massless three-body phase-space oracle; requires NumPy."""

import argparse
import json
from math import pi, sqrt
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


def bose(energy, temperature):
    """Use decaying exponentials to keep energetic tails finite."""
    x = np.asarray(energy) / temperature
    return np.exp(-x) / -np.expm1(-x)


def pair_occupation(energy, momentum, invariant, temperature):
    """Angular average of n_B(E2)n_B(E3) for a massless two-body pair."""
    high = (energy + momentum) / 2
    low = invariant / (2 * (energy + momentum))
    logarithm = np.log(-np.expm1(-high / temperature)) - np.log(
        -np.expm1(-low / temperature)
    )
    correction = np.divide(
        2 * temperature * logarithm,
        momentum,
        out=np.array(2 * bose(energy / 2, temperature) + np.zeros_like(momentum)),
        where=momentum > 1e-9 * energy,
    )
    return bose(energy, temperature) * (1 + correction)


def phase_space(momentum, temperature, mass, coupling, order):
    """Factorize one bath leg and the other two legs, without LU expressions."""
    nodes, weights = leggauss(order)
    cosine, angular_weight = nodes[:, None], weights[:, None]
    unit, radial_weight = (nodes[None, :] + 1) / 2, weights[None, :] / 2
    energy = sqrt(mass * mass + momentum * momentum)

    # P -> k + Q: Q^2 >= 0 sets the spectator's upper endpoint.
    maximum = mass * mass / (2 * (energy - momentum * cosine))
    endpoint = (1 - unit) ** 4
    k = maximum * (1 - endpoint)
    decay_jacobian = maximum * 4 * (1 - unit) ** 3
    pair_energy = energy - k
    pair_momentum = np.sqrt(momentum**2 + k**2 - 2 * momentum * k * cosine)
    invariant = mass * mass * endpoint
    decay = np.sum(
        angular_weight * radial_weight * decay_jacobian * k * bose(k, temperature)
        * pair_occupation(pair_energy, pair_momentum, invariant, temperature)
    )

    # P + k -> Q: integrate the positive half-line in the bath frame.
    k = temperature * unit / (1 - unit)
    jacobian = temperature / (1 - unit) ** 2
    pair_energy = energy + k
    pair_momentum = np.sqrt(momentum**2 + k**2 + 2 * momentum * k * cosine)
    invariant = mass * mass + 2 * k * (energy - momentum * cosine)
    scattering = np.sum(
        angular_weight * radial_weight * jacobian * k * (1 + bose(k, temperature))
        * pair_occupation(pair_energy, pair_momentum, invariant, temperature)
    )
    prefactor = coupling**2 / (128 * pi**3 * energy)
    source_decay = float(prefactor * decay)
    source_scattering = float(3 * prefactor * scattering)
    source = source_decay + source_scattering
    vacuum = coupling**2 * mass**2 / (512 * pi**3 * energy)
    return {
        "order": order,
        "momentum": momentum,
        "energy": energy,
        "source_inverse_decay": source_decay,
        "source_scattering_three_channels": source_scattering,
        "source": source,
        "spectral_rate": source / float(bose(energy, temperature)),
        "rate_per_d3p": source / (2 * pi) ** 3,
        "vacuum_spectral_rate": vacuum,
        "mb_source_inverse_decay": vacuum * np.exp(-energy / temperature),
        "mb_source_scattering_three_channels": (
            3 * coupling**2 * temperature**2 * np.exp(-energy / temperature)
            / (64 * pi**3 * energy)
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--momenta", nargs="+", type=float, default=[0, 1, 3])
    parser.add_argument("--orders", nargs="+", type=int, default=[32, 64, 128, 256])
    parser.add_argument("--temperature", type=float, default=1)
    parser.add_argument("--mass", type=float, default=1)
    parser.add_argument("--coupling", type=float, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.temperature <= 0 or args.mass <= 0 or min(args.orders) < 2:
        parser.error("positive temperature and mass, and orders >= 2 required")
    if min(args.momenta) < 0:
        parser.error("momentum magnitudes must be nonnegative")
    output = {
        "observable": "Gamma_rho = source/n_B(E); source = df_X/dt at f_X=0",
        "method": "independent two-body factorization and Gauss-Legendre quadrature",
        "temperature": args.temperature,
        "mass": args.mass,
        "coupling": args.coupling,
        "chemical_potentials": [0, 0, 0, 0],
        "results": [
            phase_space(p, args.temperature, args.mass, args.coupling, order)
            for p in args.momenta for order in args.orders
        ],
    }
    encoded = json.dumps(output, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
