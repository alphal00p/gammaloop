#!/usr/bin/env python3
"""Independent massive scalar-bubble rates; optional quadrature requires NumPy."""

import argparse
import json
from dataclasses import asdict, dataclass
from math import exp, expm1, isfinite, log1p, pi, sqrt
from pathlib import Path


@dataclass(frozen=True)
class ScalarBubble:
    """Distinct real daughters of -g X B1 B2, with zero chemical potentials."""

    energy: float
    momentum: float
    mass1: float
    mass2: float
    temperature: float
    coupling: float = 1.0

    @property
    def invariant(self):
        return self.energy**2 - self.momentum**2

    @property
    def domain(self):
        if self.invariant < 0:
            return "spacelike"
        if self.invariant > (self.mass1 + self.mass2) ** 2:
            return "decay"
        if self.invariant < (self.mass1 - self.mass2) ** 2:
            return "absorption"
        return "closed"

    def bose(self, energy):
        scaled = energy / self.temperature
        return exp(-scaled) / -expm1(-scaled)

    def occupation_integral(self, lower, upper):
        """Integral of n_B(e) de, evaluated with decaying exponentials."""
        return self.temperature * log1p(
            self.bose(lower) * -expm1(-(upper - lower) / self.temperature)
        )

    def energy_intervals(self):
        """Bath-frame daughter energies after the two-body angular integral."""
        s = self.invariant
        if self.domain == "closed":
            return []
        if self.domain == "spacelike":
            intervals = []
            for incoming, outgoing in [
                (self.mass1, self.mass2),
                (self.mass2, self.mass1),
            ]:
                difference = s + incoming**2 - outgoing**2
                lower = (
                    self.energy * difference
                    + self.momentum * sqrt(difference**2 - 4 * s * incoming**2)
                ) / (-2 * s)
                intervals.append((lower, float("inf")))
            return intervals
        discriminant = sqrt(
            (s - (self.mass1 + self.mass2) ** 2) * (s - (self.mass1 - self.mass2) ** 2)
        )
        half_width = self.momentum * discriminant / (2 * s)
        if self.domain == "decay":
            center = self.energy * (s + self.mass1**2 - self.mass2**2) / (2 * s)
        else:
            center = self.energy * (abs(self.mass1**2 - self.mass2**2) - s) / (2 * s)
        return [(center - half_width, center + half_width)]

    def rates(self):
        vacuum, thermal = 0.0, 0.0
        if self.domain != "closed":
            if self.momentum == 0:
                s = self.invariant
                discriminant = sqrt(
                    (s - (self.mass1 + self.mass2) ** 2)
                    * (s - (self.mass1 - self.mass2) ** 2)
                )
                prefactor = self.coupling**2 * discriminant / (16 * pi * self.energy**3)
                lower, _ = self.energy_intervals()[0]
                if self.domain == "decay":
                    vacuum = prefactor
                    thermal = prefactor * (
                        self.bose(lower) + self.bose(self.energy - lower)
                    )
                else:
                    thermal = prefactor * (
                        self.bose(lower) - self.bose(self.energy + lower)
                    )
            else:
                prefactor = self.coupling**2 / (16 * pi * self.energy * self.momentum)
                for lower, upper in self.energy_intervals():
                    if self.domain == "decay":
                        vacuum = prefactor * (upper - lower)
                        thermal += prefactor * (
                            self.occupation_integral(lower, upper)
                            + self.occupation_integral(
                                self.energy - upper, self.energy - lower
                            )
                        )
                    elif self.domain == "absorption":
                        thermal += prefactor * (
                            self.occupation_integral(lower, upper)
                            - self.occupation_integral(
                                lower + self.energy, upper + self.energy
                            )
                        )
                    else:
                        thermal += prefactor * self.occupation_integral(
                            lower, lower + self.energy
                        )
        return {
            "spectral_rate": vacuum + thermal,
            "vacuum_spectral_rate": vacuum,
            "thermal_excess": thermal,
        }

    def rest_mass_derivatives(self):
        """Differentiate daughter mass squared at fixed external energy and T."""
        if self.momentum != 0:
            return None
        if self.domain == "closed":
            return None  # Endpoints may have singular one-sided derivatives.
        s, a, b = self.invariant, self.mass1**2, self.mass2**2
        discriminant = sqrt(
            (s - (self.mass1 + self.mass2) ** 2) * (s - (self.mass1 - self.mass2) ** 2)
        )
        prefactor = self.coupling**2 / (16 * pi * self.energy**3)
        energy, _ = self.energy_intervals()[0]
        n1 = self.bose(energy)
        n2 = self.bose(
            self.energy - energy if self.domain == "decay" else self.energy + energy
        )
        if self.domain == "decay":
            weight = n1 + n2
            derivative_weight = (-n1 * (1 + n1) + n2 * (1 + n2)) / (
                2 * self.energy * self.temperature
            )
        else:
            weight = n1 - n2
            derivative_weight = (n1 * (1 + n1) - n2 * (1 + n2)) / (
                2 * self.energy * self.temperature
            )
            if a > b:
                derivative_weight = -derivative_weight
        result = {}
        for name, derivative_discriminant, derivative_occupation in [
            ("mass1_squared", (a - s - b) / discriminant, derivative_weight),
            ("mass2_squared", (b - s - a) / discriminant, -derivative_weight),
        ]:
            vacuum = (
                prefactor * derivative_discriminant if self.domain == "decay" else 0.0
            )
            thermal = prefactor * (
                derivative_discriminant * weight + discriminant * derivative_occupation
            )
            result[name] = {
                "spectral_rate": vacuum + thermal,
                "vacuum_spectral_rate": vacuum,
                "thermal_excess": thermal,
            }
        return result

    def quadrature(self, order):
        """Integrate the occupation weights directly, independently of logarithms."""
        if self.domain == "closed" or self.momentum == 0:
            return self.rates()["spectral_rate"]
        from numpy.polynomial.legendre import leggauss

        nodes, weights = leggauss(order)
        total = 0.0
        for lower, upper in self.energy_intervals():
            for node, weight in zip(nodes, weights):
                unit = (float(node) + 1) / 2
                if self.domain == "spacelike":
                    energy = lower + self.temperature * unit / (1 - unit)
                    jacobian = self.temperature / (1 - unit) ** 2
                else:
                    energy = lower + (upper - lower) * unit
                    jacobian = upper - lower
                if self.domain == "decay":
                    occupation = 1 + self.bose(energy) + self.bose(self.energy - energy)
                else:
                    occupation = self.bose(energy) - self.bose(self.energy + energy)
                total += float(weight) * jacobian * occupation / 2
        return self.coupling**2 * total / (16 * pi * self.energy * self.momentum)


CASES = {
    "rest_decay": ScalarBubble(2, 0, 0.4, 0.7, 0.6),
    "boosted_decay": ScalarBubble(2.5, 1.5, 0.4, 0.7, 0.6),
    "closed": ScalarBubble(1, 0, 0.6, 0.8, 0.6),
    "absorption": ScalarBubble(1, 0, 0.3, 1.8, 0.6),
    "spacelike": ScalarBubble(1, 2, 0.5, 0.5, 1),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CASES)
    parser.add_argument("--energy", type=float)
    parser.add_argument("--momentum", type=float)
    parser.add_argument("--masses", nargs=2, type=float)
    parser.add_argument("--temperature", type=float)
    parser.add_argument("--coupling", type=float, default=1)
    parser.add_argument("--orders", nargs="+", type=int, default=[])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    custom = [args.energy, args.momentum, args.masses, args.temperature]
    if any(value is not None for value in custom):
        if args.case or any(value is None for value in custom):
            parser.error(
                "custom points require --energy, --momentum, --masses and --temperature, without --case"
            )
        cases = {
            "custom": ScalarBubble(
                args.energy,
                args.momentum,
                *args.masses,
                args.temperature,
                args.coupling,
            )
        }
    else:
        selected = {args.case: CASES[args.case]} if args.case else CASES
        cases = {
            name: ScalarBubble(**(asdict(point) | {"coupling": args.coupling}))
            for name, point in selected.items()
        }
    if any(order < 2 for order in args.orders):
        parser.error("quadrature orders must be at least two")
    results = []
    for name, point in cases.items():
        if (
            not all(isfinite(value) for value in asdict(point).values())
            or min(point.energy, point.mass1, point.mass2, point.temperature) <= 0
            or point.momentum < 0
            or point.invariant == 0
        ):
            parser.error(
                "require finite parameters, positive E, masses and T, p >= 0, and E² != p²"
            )
        results.append(
            {
                "case": name,
                **asdict(point),
                "invariant": point.invariant,
                "domain": point.domain,
                **point.rates(),
                "rest_mass_squared_derivatives": point.rest_mass_derivatives(),
                "quadrature": [
                    {"order": order, "spectral_rate": point.quadrature(order)}
                    for order in args.orders
                ],
            }
        )
    output = {
        "observable": "Gamma_rho = -Im Pi_R(P)/P0, with distinct scalar daughters",
        "method": "closed two-body phase-space formulas; optional occupation quadrature",
        "chemical_potentials": [0, 0, 0],
        "results": results,
    }
    encoded = json.dumps(output, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
