#!/usr/bin/env python3
"""Dependency-free quadrature checks for LU_H_MATCHED_SAMPLING.md.

This integrates one-dimensional model densities, not a GammaLoop integrand.
Run: python3 docs/research/advanced_sampling/lu_h_profile_checks.py
"""

import json
import math


def integral(function, intervals=16000):
    low, high = -7.0, 7.0
    step = (high - low) / intervals
    total = function(low) + function(high)
    for index in range(1, intervals):
        total += (4 if index % 2 else 2) * function(low + index * step)
    return total * step / 3


def peak(function):
    # Locate a bracket on a fixed grid, then maximize by golden-section search.
    points = [-7.0 + index * 0.002 for index in range(7001)]
    best = max(range(len(points)), key=lambda index: function(points[index]))
    low, high = points[max(0, best - 1)], points[min(best + 1, len(points) - 1)]
    ratio = (math.sqrt(5) - 1) / 2
    for _ in range(80):
        left, right = high - ratio * (high - low), low + ratio * (high - low)
        if function(left) > function(right):
            high = right
        else:
            low = left
    return function((low + high) / 2)


def profile(family, power):
    exponent = 2 if family == "poly_exponential" else 1
    raw = lambda y: math.exp(-power * y + 2 - 2 * math.cosh(exponent * y))
    normalization = integral(lambda y: math.exp(y) * raw(y))
    h = lambda y: raw(y) / normalization
    mode = math.asinh((1 - power) / (2 * exponent)) / exponent
    shape = 2 * exponent * math.sqrt(math.cosh(exponent * mode))
    density = lambda y: shape / (4 * math.exp(y) * math.cosh(shape * (y - mode) / 2) ** 2)
    return h, density, {"log_median": mode, "shape": shape, "h_normalization": normalization}


def row(name, h, density, **extra):
    value = lambda y: math.exp(y) * h(y) ** 2 / density(y) if h(y) else 0.0
    coarse = integral(value, 8000)
    fine = integral(value, 16000)
    return dict(
        proposal=name,
        second_moment=fine,
        quadrature_refinement_difference=abs(fine - coarse),
        peak_h_over_q=peak(lambda y: h(y) / density(y) if h(y) else 0.0),
        **extra,
    )


def main():
    h, density, parameters = profile("poly_exponential", 0)
    comparisons = [row("h_matched", h, h), row("log_logistic", h, density, **parameters)]
    for root_over_beta in [0.1, 1.0, 10.0, 100.0]:
        comparisons.append(row(
            "rational_raw_radius", h,
            lambda y, scale=root_over_beta: scale / (scale + math.exp(y)) ** 2,
            root_over_beta=root_over_beta,
        ))
    powers = [0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16]
    families = []
    for family in ["poly_exponential", "poly_left_right_exponential"]:
        for power in powers:
            hi, qi, params = profile(family, power)
            families.append(row("log_logistic", hi, qi, family=family, power=power, **params))
    cdf = lambda z: 0.5 * (math.erfc(1 / z - z) - math.exp(4) * math.erfc(z + 1 / z))
    residuals = []
    for z in [0.4, 0.7, 1.0, 1.7, 3.0]:
        step = z * 1e-5
        derivative = (cdf(z + step) - cdf(z - step)) / (2 * step)
        residuals.append(abs(derivative / h(math.log(z)) - 1))
    print(json.dumps(dict(
        method="Composite Simpson in log(t) on [-7,7], sigma=1; 8000/16000 intervals",
        scope="One-dimensional simple-residue model; no GammaLoop integration",
        default_comparisons=comparisons,
        all_current_polynomial_profiles=families,
        p0_cdf_derivative_max_relative_residual=max(residuals),
    ), indent=2))


if __name__ == "__main__":
    main()
