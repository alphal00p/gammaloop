#!/usr/bin/env python3
"""Standalone proxy for the GL00 double-UV relative-momentum ridge.

The proxy is intentionally small, but it keeps the structural element that
matters for the GL00 max-weight events in an LMB such as (e5,e10):

    q = v - u + Q

where u and v are sampled hard LMB momenta and Q is a fixed external shift.  A
local 3D UV approximation of the relative edge drops Q,

    q_UV = v - u,

so the locally subtracted relative-edge factor contains

    1 / E(q) - 1 / E(q_UV).

For generic double-UV rays, q and q_UV both scale with the common UV radius and
the difference falls.  If the rays are exactly collinear with equal norm, q_UV
stays small while u and v become hard; the same expression grows with the
common radius.  At fixed common radius, the q_UV -> 0 behavior is the
integrable 3D ridge 1 / |q_UV|.

The optional integrated-CT switch only multiplies the same relative-edge
structure by a smooth numerator-like factor, mimicking the fact that the nested
integrated self-energy has a p-slash + m identity Lorentz structure.  It must
not change the geometric singularity.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Callable, NamedTuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"
DATA = ROOT / "data" / "mock_gl00_proxy"

M_LOOP = 1.0
M_REL_PHYS = 0.35
M_REL_UV = 1.0e-6
M_REL_UV_4D = 9.188
Q_SHIFT = np.array([0.73, -0.21, 0.44])


def norm(vec: np.ndarray) -> float:
    return float(np.linalg.norm(vec))


def unit(vec: np.ndarray) -> np.ndarray:
    length = norm(vec)
    if length == 0.0:
        raise ValueError("cannot normalize a zero vector")
    return vec / length


N_HARD = unit(np.array([0.0933320938357391, 0.7051483309767351, -0.7028904264399691]))
N_PERP = unit(np.cross(N_HARD, np.array([0.0, 0.0, 1.0])))
N_RANDOM = unit(np.array([-0.43, 0.18, 0.885]))


class MCResult(NamedTuple):
    mode: str
    include_integrated_ct: bool
    domain: str
    channel: str
    sample_count: int
    estimate: float
    standard_error: float
    relative_error: float
    sample_variance: float
    variance_reduction_vs_uniform: float | None
    max_abs_weight: float
    p99_abs_weight: float
    p999_abs_weight: float


class ChannelPartitionRow(NamedTuple):
    path: str
    radius: float
    relative_radius: float
    raw_abs_weight: float
    bad_channel_probability_ose_alpha1: float
    bad_channel_probability_ose_alpha3: float
    bad_channel_probability_ose_alpha4: float
    bad_channel_probability_jacobian_linear: float
    bad_channel_probability_jacobian_log: float
    bad_channel_weighted_abs_ose_alpha1: float
    bad_channel_weighted_abs_ose_alpha3: float
    bad_channel_weighted_abs_ose_alpha4: float
    bad_channel_weighted_abs_jacobian_linear: float
    bad_channel_weighted_abs_jacobian_log: float


def energy(vec: np.ndarray, mass: float) -> float:
    return math.sqrt(float(np.dot(vec, vec)) + mass * mass)


def common_prefactor(u: np.ndarray, v: np.ndarray) -> float:
    """Jacobian-like hard prefactor left after the other smooth factors."""

    eu = energy(u, M_LOOP)
    ev = energy(v, M_LOOP)
    common = 0.5 * (eu + ev)
    return common**3 / (eu * ev)


def common_prefactor_4d_cff(u: np.ndarray, v: np.ndarray) -> float:
    """Hard prefactor for the 4D-first vacuum-CFF proxy.

    The extra power compensates the additional hard CFF energy denominator in
    `vacuum_cff_relative_factor`, so generic shifted-vs-unshifted rays still
    fall like 1/R while the unresolved relative ray grows like R.
    """

    eu = energy(u, M_LOOP)
    ev = energy(v, M_LOOP)
    common = 0.5 * (eu + ev)
    return common**4 / (eu * ev)


def numerator_factor(
    u: np.ndarray, v: np.ndarray, include_integrated_ct: bool
) -> float:
    """Smooth numerator proxy; the integrated CT only changes this coefficient."""

    eu = energy(u, M_LOOP)
    ev = energy(v, M_LOOP)
    cos_uv = float(np.dot(u, v)) / (eu * ev)
    original_numerator = 1.0 + 0.15 * cos_uv

    if not include_integrated_ct:
        return original_numerator

    integrated_self_energy = 0.42 * (1.0 - 0.08 * cos_uv)
    return original_numerator + integrated_self_energy


def original_integrand(
    u: np.ndarray, v: np.ndarray, include_integrated_ct: bool
) -> float:
    q = v - u + Q_SHIFT
    return (
        numerator_factor(u, v, include_integrated_ct)
        * common_prefactor(u, v)
        / energy(q, M_REL_PHYS)
    )


def local_double_uv_ct(
    u: np.ndarray, v: np.ndarray, include_integrated_ct: bool
) -> float:
    q_uv = v - u
    return (
        numerator_factor(u, v, include_integrated_ct)
        * common_prefactor(u, v)
        / energy(q_uv, M_REL_UV)
    )


def subtracted_integrand(
    u: np.ndarray, v: np.ndarray, include_integrated_ct: bool
) -> float:
    return original_integrand(u, v, include_integrated_ct) - local_double_uv_ct(
        u, v, include_integrated_ct
    )


def vacuum_cff_relative_factor(
    rel_spatial: np.ndarray, hard_radius: float, mass: float
) -> float:
    """Minimal CFF factor produced by a massive 4D vacuum relative propagator."""

    rel_energy = energy(rel_spatial, mass)
    return 1.0 / (rel_energy * (2.0 * hard_radius + rel_energy))


def original_integrand_4d_cff(
    u: np.ndarray, v: np.ndarray, include_integrated_ct: bool
) -> float:
    hard_radius = 0.5 * (energy(u, M_LOOP) + energy(v, M_LOOP))
    return (
        numerator_factor(u, v, include_integrated_ct)
        * common_prefactor_4d_cff(u, v)
        * vacuum_cff_relative_factor(v - u + Q_SHIFT, hard_radius, M_REL_PHYS)
    )


def local_double_uv_ct_4d_cff(
    u: np.ndarray, v: np.ndarray, include_integrated_ct: bool
) -> float:
    hard_radius = 0.5 * (energy(u, M_LOOP) + energy(v, M_LOOP))
    return (
        numerator_factor(u, v, include_integrated_ct)
        * common_prefactor_4d_cff(u, v)
        * vacuum_cff_relative_factor(v - u, hard_radius, M_REL_UV_4D)
    )


def subtracted_integrand_4d_cff(
    u: np.ndarray, v: np.ndarray, include_integrated_ct: bool
) -> float:
    return original_integrand_4d_cff(
        u, v, include_integrated_ct
    ) - local_double_uv_ct_4d_cff(u, v, include_integrated_ct)


def energy_array(vec: np.ndarray, mass: float) -> np.ndarray:
    return np.sqrt(np.sum(vec * vec, axis=-1) + mass * mass)


def numerator_factor_array(
    u: np.ndarray, v: np.ndarray, include_integrated_ct: bool
) -> np.ndarray:
    eu = energy_array(u, M_LOOP)
    ev = energy_array(v, M_LOOP)
    cos_uv = np.sum(u * v, axis=-1) / (eu * ev)
    original_numerator = 1.0 + 0.15 * cos_uv
    if not include_integrated_ct:
        return original_numerator
    integrated_self_energy = 0.42 * (1.0 - 0.08 * cos_uv)
    return original_numerator + integrated_self_energy


def common_prefactor_array(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    eu = energy_array(u, M_LOOP)
    ev = energy_array(v, M_LOOP)
    common = 0.5 * (eu + ev)
    return common**3 / (eu * ev)


def common_prefactor_4d_cff_array(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    eu = energy_array(u, M_LOOP)
    ev = energy_array(v, M_LOOP)
    common = 0.5 * (eu + ev)
    return common**4 / (eu * ev)


def vacuum_cff_relative_factor_array(
    rel_spatial: np.ndarray, hard_radius: np.ndarray, mass: float
) -> np.ndarray:
    rel_energy = energy_array(rel_spatial, mass)
    return 1.0 / (rel_energy * (2.0 * hard_radius + rel_energy))


def subtracted_integrand_array(
    common_momentum: np.ndarray,
    relative_momentum: np.ndarray,
    include_integrated_ct: bool,
    mode: str,
) -> np.ndarray:
    u = common_momentum - 0.5 * relative_momentum
    v = common_momentum + 0.5 * relative_momentum
    numerator = numerator_factor_array(u, v, include_integrated_ct)

    if mode == "3d_cff_expansion":
        prefactor = common_prefactor_array(u, v)
        original = prefactor / energy_array(relative_momentum + Q_SHIFT, M_REL_PHYS)
        local = prefactor / energy_array(relative_momentum, M_REL_UV)
    elif mode == "4d_first_then_cff":
        hard_radius = 0.5 * (energy_array(u, M_LOOP) + energy_array(v, M_LOOP))
        prefactor = common_prefactor_4d_cff_array(u, v)
        original = prefactor * vacuum_cff_relative_factor_array(
            relative_momentum + Q_SHIFT, hard_radius, M_REL_PHYS
        )
        local = prefactor * vacuum_cff_relative_factor_array(
            relative_momentum, hard_radius, M_REL_UV_4D
        )
    else:
        raise ValueError(f"unknown mode {mode!r}")

    return numerator * (original - local)


def rotate_from_hard(theta: float) -> np.ndarray:
    return math.cos(theta) * N_HARD + math.sin(theta) * N_PERP


def ray_locked_equal(radius: float) -> tuple[np.ndarray, np.ndarray]:
    return radius * N_HARD, radius * N_HARD


def ray_parallel_unequal(radius: float) -> tuple[np.ndarray, np.ndarray]:
    return radius * N_HARD, 1.04 * radius * N_HARD


def ray_small_angle(radius: float) -> tuple[np.ndarray, np.ndarray]:
    return radius * N_HARD, radius * rotate_from_hard(math.radians(1.0))


def ray_wider_angle(radius: float) -> tuple[np.ndarray, np.ndarray]:
    return radius * N_HARD, radius * rotate_from_hard(math.radians(8.0))


def ray_random(radius: float) -> tuple[np.ndarray, np.ndarray]:
    return radius * N_HARD, radius * N_RANDOM


RAY_CASES: list[tuple[str, Callable[[float], tuple[np.ndarray, np.ndarray]]]] = [
    ("locked equal", ray_locked_equal),
    ("parallel, 4% radius mismatch", ray_parallel_unequal),
    ("equal radius, 1 deg opening", ray_small_angle),
    ("equal radius, 8 deg opening", ray_wider_angle),
    ("generic directions", ray_random),
]


def fit_slope(x: np.ndarray, y: np.ndarray, min_x: float) -> float:
    mask = (x >= min_x) & np.isfinite(y) & (np.abs(y) > 0.0)
    log_x = np.log(x[mask])
    log_y = np.log(np.abs(y[mask]))
    slope, _ = np.polyfit(log_x, log_y, 1)
    return float(slope)


def stable_logsumexp_pair(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    reference = np.maximum(first, second)
    return reference + np.log(np.exp(first - reference) + np.exp(second - reference))


def spherical_map_log_jacobian(
    radius: np.ndarray, mapping: str, e_cm: float = 91.188, b: float = 1.0
) -> np.ndarray:
    safe_radius = np.maximum(radius, np.finfo(float).tiny)
    if mapping == "linear":
        log_dr_dx = 2.0 * np.log(e_cm * b + safe_radius) - math.log(e_cm * b)
    elif mapping == "log":
        z = safe_radius / e_cm
        log_a_plus_b = np.empty_like(z)
        small_z = z < 50.0
        log_a_plus_b[small_z] = np.log(np.expm1(z[small_z]) + b)
        log_a_plus_b[~small_z] = z[~small_z] + np.log1p(
            (b - 1.0) * np.exp(-z[~small_z])
        )
        log_dr_dx = math.log(e_cm / b) + 2.0 * log_a_plus_b - z
    else:
        raise ValueError(f"unknown mapping {mapping!r}")
    return math.log(4.0 * math.pi) + 2.0 * np.log(safe_radius) + log_dr_dx


def bad_channel_probability_from_log_weights(
    direct_log_weight: np.ndarray, relative_log_weight: np.ndarray
) -> np.ndarray:
    return np.exp(
        direct_log_weight
        - stable_logsumexp_pair(direct_log_weight, relative_log_weight)
    )


def ose_bad_channel_probability(
    u: np.ndarray, v: np.ndarray, common_momentum: np.ndarray, alpha: float
) -> np.ndarray:
    relative_momentum = v - u
    direct_log_weight = -alpha * (
        np.log(energy_array(u, M_LOOP)) + np.log(energy_array(v, M_LOOP))
    )
    relative_log_weight = -alpha * (
        np.log(energy_array(common_momentum, M_LOOP))
        + np.log(energy_array(relative_momentum, M_REL_UV))
    )
    return bad_channel_probability_from_log_weights(
        direct_log_weight, relative_log_weight
    )


def jacobian_bad_channel_probability(
    u: np.ndarray, v: np.ndarray, common_momentum: np.ndarray, mapping: str
) -> np.ndarray:
    relative_momentum = v - u
    direct_log_weight = -(
        spherical_map_log_jacobian(np.linalg.norm(u, axis=-1), mapping)
        + spherical_map_log_jacobian(np.linalg.norm(v, axis=-1), mapping)
    )
    relative_log_weight = -(
        spherical_map_log_jacobian(np.linalg.norm(common_momentum, axis=-1), mapping)
        + spherical_map_log_jacobian(
            np.linalg.norm(relative_momentum, axis=-1), mapping
        )
    )
    return bad_channel_probability_from_log_weights(
        direct_log_weight, relative_log_weight
    )


def radial_scan(include_integrated_ct: bool, mode: str) -> list[dict[str, float | str]]:
    radii = np.logspace(2.0, 8.0, 220)
    rows: list[dict[str, float | str]] = []
    for label, ray in RAY_CASES:
        for radius in radii:
            u, v = ray(float(radius))
            if mode == "3d_cff_expansion":
                original = original_integrand(u, v, include_integrated_ct)
                local_ct = local_double_uv_ct(u, v, include_integrated_ct)
                subtracted = subtracted_integrand(u, v, include_integrated_ct)
            elif mode == "4d_first_then_cff":
                original = original_integrand_4d_cff(u, v, include_integrated_ct)
                local_ct = local_double_uv_ct_4d_cff(u, v, include_integrated_ct)
                subtracted = subtracted_integrand_4d_cff(u, v, include_integrated_ct)
            else:
                raise ValueError(f"unknown mode {mode!r}")
            rows.append(
                {
                    "mode": mode,
                    "include_integrated_ct": str(include_integrated_ct),
                    "case": label,
                    "radius": float(radius),
                    "original": original,
                    "local_double_uv_ct": local_ct,
                    "subtracted": subtracted,
                    "relative_norm": norm(v - u),
                }
            )
    return rows


def summarize_radial(
    rows: list[dict[str, float | str]],
) -> list[dict[str, float | str]]:
    summary = []
    modes = sorted({str(row["mode"]) for row in rows})
    cases = sorted({str(row["case"]) for row in rows})
    for mode in modes:
        for include in ["False", "True"]:
            for case in cases:
                selected = [
                    row
                    for row in rows
                    if row["mode"] == mode
                    and row["include_integrated_ct"] == include
                    and row["case"] == case
                ]
                radii = np.array([float(row["radius"]) for row in selected])
                sub = np.array([float(row["subtracted"]) for row in selected])
                orig = np.array([float(row["original"]) for row in selected])
                ct = np.array([float(row["local_double_uv_ct"]) for row in selected])
                summary.append(
                    {
                        "mode": mode,
                        "include_integrated_ct": include,
                        "case": case,
                        "subtracted_tail_slope": fit_slope(radii, sub, 1.0e5),
                        "original_tail_slope": fit_slope(radii, orig, 1.0e5),
                        "local_ct_tail_slope": fit_slope(radii, ct, 1.0e5),
                        "tail_relative_norm_over_radius": float(
                            np.mean(
                                [
                                    float(row["relative_norm"]) / float(row["radius"])
                                    for row in selected
                                    if float(row["radius"]) >= 1.0e5
                                ]
                            )
                        ),
                    }
                )
    return summary


def write_radial_outputs(
    rows: list[dict[str, float | str]], summary: list[dict[str, float | str]]
) -> None:
    DATA.mkdir(parents=True, exist_ok=True)

    with (DATA / "radial_scan.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "include_integrated_ct",
                "mode",
                "case",
                "radius",
                "original",
                "local_double_uv_ct",
                "subtracted",
                "relative_norm",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    (DATA / "summary.json").write_text(json.dumps(summary, indent=2))


def plot_radial(
    rows: list[dict[str, float | str]], summary: list[dict[str, float | str]]
) -> None:
    FIGURES.mkdir(exist_ok=True)

    summary_by_key = {
        (str(row["mode"]), str(row["include_integrated_ct"]), str(row["case"])): float(
            row["subtracted_tail_slope"]
        )
        for row in summary
    }

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.5), sharex=True, sharey=True)
    for row_index, mode in enumerate(["3d_cff_expansion", "4d_first_then_cff"]):
        for axis, include in zip(axes[row_index], ["False", "True"], strict=True):
            for case, _ in RAY_CASES:
                selected = [
                    row
                    for row in rows
                    if row["mode"] == mode
                    and row["include_integrated_ct"] == include
                    and row["case"] == case
                ]
                radii = np.array([float(row["radius"]) for row in selected])
                values = np.abs(
                    np.array([float(row["subtracted"]) for row in selected])
                )
                slope = summary_by_key[(mode, include, case)]
                axis.loglog(radii, values, label=f"{case}, slope {slope:.2f}")

            title = (
                "without integrated CT" if include == "False" else "with integrated CT"
            )
            mode_title = (
                "3D expansion" if mode == "3d_cff_expansion" else "4D first, then CFF"
            )
            axis.set_title(f"{mode_title}; {title}")
            axis.set_xlabel("common UV radius R")
            axis.grid(True, which="both", alpha=0.25)
            axis.legend(fontsize=6.5)

    axes[0, 0].set_ylabel("|original - local double-UV CT|")
    axes[1, 0].set_ylabel("|original - local double-UV CT|")
    fig.suptitle(
        "Mock GL00 proxy: the locked equal-norm ray survives 4D-first expansion"
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "mock_gl00_proxy_radial.png", dpi=180)
    plt.close(fig)


def plot_components(rows: list[dict[str, float | str]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), sharex=True)
    for axis, case in zip(
        axes, ["locked equal", "equal radius, 1 deg opening"], strict=True
    ):
        selected = [
            row
            for row in rows
            if row["mode"] == "3d_cff_expansion"
            and row["include_integrated_ct"] == "False"
            and row["case"] == case
        ]
        radii = np.array([float(row["radius"]) for row in selected])
        for key, label in [
            ("original", "original"),
            ("local_double_uv_ct", "local UV CT"),
            ("subtracted", "subtracted"),
        ]:
            values = np.abs(np.array([float(row[key]) for row in selected]))
            axis.loglog(radii, values, label=label)
        axis.set_title(case)
        axis.set_xlabel("common UV radius R")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=8)

    axes[0].set_ylabel("absolute weight")
    fig.suptitle("Original and local-CT pieces in the no-integrated-CT proxy")
    fig.tight_layout()
    fig.savefig(FIGURES / "mock_gl00_proxy_components.png", dpi=180)
    plt.close(fig)


def relative_ridge_scan(
    include_integrated_ct: bool, mode: str
) -> list[dict[str, float | str]]:
    radius = 1.0e6
    deltas = np.logspace(-5.0, 1.0, 240)
    rows = []
    u = radius * N_HARD
    for delta in deltas:
        v = u + delta * N_PERP
        if mode == "3d_cff_expansion":
            weight = subtracted_integrand(u, v, include_integrated_ct)
        elif mode == "4d_first_then_cff":
            weight = subtracted_integrand_4d_cff(u, v, include_integrated_ct)
        else:
            raise ValueError(f"unknown mode {mode!r}")
        rows.append(
            {
                "mode": mode,
                "include_integrated_ct": str(include_integrated_ct),
                "delta": float(delta),
                "abs_subtracted": abs(weight),
                "measure_weighted_delta2_abs_subtracted": float(
                    delta * delta * abs(weight)
                ),
            }
        )
    return rows


def write_and_plot_relative_ridge(
    rows: list[dict[str, float | str]],
) -> dict[str, float]:
    with (DATA / "relative_ridge_scan.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "include_integrated_ct",
                "mode",
                "delta",
                "abs_subtracted",
                "measure_weighted_delta2_abs_subtracted",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    summary = {}
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0), sharex=True)
    for mode, axis_row in zip(
        ["3d_cff_expansion", "4d_first_then_cff"], axes, strict=True
    ):
        for include, style in [("False", "-"), ("True", "--")]:
            selected = [
                row
                for row in rows
                if row["mode"] == mode and row["include_integrated_ct"] == include
            ]
            delta = np.array([float(row["delta"]) for row in selected])
            weight = np.array([float(row["abs_subtracted"]) for row in selected])
            measured = np.array(
                [
                    float(row["measure_weighted_delta2_abs_subtracted"])
                    for row in selected
                ]
            )
            mask = (delta >= 1.0e-4) & (delta <= 1.0e-2)
            ridge_slope = float(
                np.polyfit(np.log(delta[mask]), np.log(weight[mask]), 1)[0]
            )
            summary[f"ridge_slope_{mode}_include_integrated_ct_{include}"] = ridge_slope

            label = (
                "with integrated CT" if include == "True" else "without integrated CT"
            )
            axis_row[0].loglog(
                delta, weight, style, label=f"{label}, slope {ridge_slope:.2f}"
            )
            axis_row[1].loglog(delta, measured, style, label=label)

        mode_title = (
            "3D expansion" if mode == "3d_cff_expansion" else "4D first, then CFF"
        )
        axis_row[0].set_title(f"{mode_title}: relative ridge")
        axis_row[1].set_title(f"{mode_title}: with 3D relative measure")
        axis_row[0].set_ylabel("|subtracted proxy|")
        axis_row[1].set_ylabel(r"$\delta^2 |$subtracted proxy$|$")
        for axis in axis_row:
            axis.set_xlabel(r"relative momentum size $\delta=|v-u|$")
            axis.grid(True, which="both", alpha=0.25)
            axis.legend(fontsize=8)
    fig.suptitle("Local integrability of the relative-momentum ridge")
    fig.tight_layout()
    fig.savefig(FIGURES / "mock_gl00_proxy_relative_ridge.png", dpi=180)
    plt.close(fig)

    return summary


def random_unit_vectors(rng: np.random.Generator, sample_count: int) -> np.ndarray:
    z = rng.uniform(-1.0, 1.0, sample_count)
    phi = rng.uniform(0.0, 2.0 * math.pi, sample_count)
    transverse = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    return np.column_stack((transverse * np.cos(phi), transverse * np.sin(phi), z))


def sample_uniform_relative_ball(
    rng: np.random.Generator, sample_count: int, q_max: float
) -> tuple[np.ndarray, np.ndarray]:
    radius = q_max * rng.random(sample_count) ** (1.0 / 3.0)
    direction = random_unit_vectors(rng, sample_count)
    q = radius[:, None] * direction
    density = np.full(sample_count, 3.0 / (4.0 * math.pi * q_max**3))
    return q, density


def sample_ridge_relative_ball(
    rng: np.random.Generator, sample_count: int, q_max: float
) -> tuple[np.ndarray, np.ndarray]:
    radius = q_max * np.sqrt(rng.random(sample_count))
    direction = random_unit_vectors(rng, sample_count)
    q = radius[:, None] * direction
    safe_radius = np.maximum(radius, np.finfo(float).tiny)
    density = 1.0 / (2.0 * math.pi * q_max**2 * safe_radius)
    return q, density


def uniform_relative_density(q_max: float) -> float:
    return 3.0 / (4.0 * math.pi * q_max**3)


def ridge_relative_density(relative_momentum: np.ndarray, q_max: float) -> np.ndarray:
    radius = np.linalg.norm(relative_momentum, axis=-1)
    safe_radius = np.maximum(radius, np.finfo(float).tiny)
    return 1.0 / (2.0 * math.pi * q_max**2 * safe_radius)


def sample_relative_momentum(
    rng: np.random.Generator,
    sample_count: int,
    q_max: float,
    channel: str,
    ridge_probability: float,
) -> tuple[np.ndarray, np.ndarray]:
    if channel == "uniform":
        q, _ = sample_uniform_relative_ball(rng, sample_count, q_max)
    elif channel == "ridge":
        q, _ = sample_ridge_relative_ball(rng, sample_count, q_max)
    elif channel == "multichannel":
        choose_ridge = rng.random(sample_count) < ridge_probability
        q = np.empty((sample_count, 3))
        if np.any(choose_ridge):
            q[choose_ridge], _ = sample_ridge_relative_ball(
                rng, int(np.sum(choose_ridge)), q_max
            )
        if np.any(~choose_ridge):
            q[~choose_ridge], _ = sample_uniform_relative_ball(
                rng, int(np.sum(~choose_ridge)), q_max
            )
    else:
        raise ValueError(f"unknown channel {channel!r}")

    uniform_density = uniform_relative_density(q_max)
    ridge_density = ridge_relative_density(q, q_max)
    if channel == "uniform":
        sampling_density = np.full(sample_count, uniform_density)
    elif channel == "ridge":
        sampling_density = ridge_density
    else:
        sampling_density = (
            1.0 - ridge_probability
        ) * uniform_density + ridge_probability * ridge_density
    return q, sampling_density


def sample_common_momentum(
    rng: np.random.Generator,
    sample_count: int,
    domain: str,
    hard_radius: float,
    radius_min: float,
    radius_max: float,
) -> tuple[np.ndarray, float]:
    if domain == "fixed_hard_radius":
        radii = np.full(sample_count, hard_radius)
        log_radius_volume = 1.0
    elif domain == "log_hard_radius":
        log_min = math.log(radius_min)
        log_max = math.log(radius_max)
        radii = np.exp(rng.uniform(log_min, log_max, sample_count))
        log_radius_volume = log_max - log_min
    else:
        raise ValueError(f"unknown domain {domain!r}")

    return radii[:, None] * N_HARD[None, :], log_radius_volume


def run_mc_integral(
    mode: str,
    include_integrated_ct: bool,
    domain: str,
    channel: str,
    sample_count: int,
    seed: int,
    q_max: float,
    hard_radius: float,
    radius_min: float,
    radius_max: float,
    ridge_probability: float,
) -> MCResult:
    rng = np.random.default_rng(seed)
    common_momentum, log_radius_volume = sample_common_momentum(
        rng, sample_count, domain, hard_radius, radius_min, radius_max
    )
    relative_momentum, q_density = sample_relative_momentum(
        rng, sample_count, q_max, channel, ridge_probability
    )
    values = subtracted_integrand_array(
        common_momentum, relative_momentum, include_integrated_ct, mode
    )
    weights = log_radius_volume * values / q_density
    estimate = float(np.mean(weights))
    sample_variance = float(np.var(weights, ddof=1))
    standard_error = math.sqrt(sample_variance / sample_count)
    relative_error = standard_error / max(abs(estimate), np.finfo(float).tiny)
    abs_weights = np.abs(weights)
    return MCResult(
        mode=mode,
        include_integrated_ct=include_integrated_ct,
        domain=domain,
        channel=channel,
        sample_count=sample_count,
        estimate=estimate,
        standard_error=standard_error,
        relative_error=relative_error,
        sample_variance=sample_variance,
        variance_reduction_vs_uniform=None,
        max_abs_weight=float(np.max(abs_weights)),
        p99_abs_weight=float(np.quantile(abs_weights, 0.99)),
        p999_abs_weight=float(np.quantile(abs_weights, 0.999)),
    )


def mc_result_to_dict(result: MCResult) -> dict[str, float | int | str | bool | None]:
    return result._asdict()


def run_mc_variance_experiment() -> list[MCResult]:
    sample_count = 250_000
    q_max = 0.2
    hard_radius = 1.0e6
    radius_min = 1.0e3
    radius_max = 1.0e7
    ridge_probability = 0.75
    seed = 20260614
    channels = ["uniform", "ridge", "multichannel"]
    domains = ["fixed_hard_radius", "log_hard_radius"]
    modes = ["3d_cff_expansion", "4d_first_then_cff"]

    results: list[MCResult] = []
    for mode in modes:
        for include_integrated_ct in [False, True]:
            for domain in domains:
                domain_results = [
                    run_mc_integral(
                        mode=mode,
                        include_integrated_ct=include_integrated_ct,
                        domain=domain,
                        channel=channel,
                        sample_count=sample_count,
                        seed=seed + 17 * len(results) + channels.index(channel),
                        q_max=q_max,
                        hard_radius=hard_radius,
                        radius_min=radius_min,
                        radius_max=radius_max,
                        ridge_probability=ridge_probability,
                    )
                    for channel in channels
                ]
                uniform_variance = next(
                    result.sample_variance
                    for result in domain_results
                    if result.channel == "uniform"
                )
                for result in domain_results:
                    reduction = (
                        None
                        if result.channel == "uniform"
                        else uniform_variance / result.sample_variance
                    )
                    results.append(
                        result._replace(variance_reduction_vs_uniform=reduction)
                    )
    return results


def write_and_plot_mc_variance(results: list[MCResult]) -> None:
    rows = [mc_result_to_dict(result) for result in results]
    (DATA / "mc_variance_summary.json").write_text(json.dumps(rows, indent=2))
    with (DATA / "mc_variance_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    plot_rows = [
        result
        for result in results
        if result.mode == "3d_cff_expansion"
        and result.domain in {"fixed_hard_radius", "log_hard_radius"}
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.5), sharey=True)
    for axis, domain in zip(
        axes, ["fixed_hard_radius", "log_hard_radius"], strict=True
    ):
        selected = [
            result
            for result in plot_rows
            if result.domain == domain and not result.include_integrated_ct
        ]
        names = [result.channel for result in selected]
        rel_errors = [result.relative_error for result in selected]
        reductions = [
            1.0
            if result.variance_reduction_vs_uniform is None
            else result.variance_reduction_vs_uniform
            for result in selected
        ]
        bars = axis.bar(names, rel_errors, color=["#7a7a7a", "#2f7fbf", "#42a45f"])
        axis.set_yscale("log")
        axis.set_title(domain.replace("_", " "))
        axis.set_ylabel("relative MC error")
        axis.grid(True, axis="y", which="both", alpha=0.25)
        for bar, reduction in zip(bars, reductions, strict=True):
            axis.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                f"{reduction:.1f}x var",
                ha="center",
                va="bottom",
                fontsize=8,
            )
    fig.suptitle("Variance reduction from a relative-momentum channel")
    fig.tight_layout()
    fig.savefig(FIGURES / "mock_gl00_proxy_mc_variance.png", dpi=180)
    plt.close(fig)


def build_channel_partition_rows() -> list[ChannelPartitionRow]:
    rows: list[ChannelPartitionRow] = []

    path_points: list[tuple[str, np.ndarray, np.ndarray]] = []
    fixed_q = 0.2
    for radius in np.logspace(2.0, 7.0, 160):
        common = radius * N_HARD
        q = fixed_q * N_PERP
        path_points.append(("fixed relative q", common, q))

    fixed_fraction = 1.0e-3
    for radius in np.logspace(2.0, 7.0, 160):
        common = radius * N_HARD
        q = fixed_fraction * radius * N_PERP
        path_points.append(("fixed q/R", common, q))

    fixed_radius = 1.0e6
    for q_radius in np.logspace(-5.0, 1.0, 160):
        common = fixed_radius * N_HARD
        q = q_radius * N_PERP
        path_points.append(("shrinking relative q", common, q))

    for path, common, q in path_points:
        common_batch = common[None, :]
        q_batch = q[None, :]
        u = common_batch - 0.5 * q_batch
        v = common_batch + 0.5 * q_batch
        raw_abs_weight = float(
            abs(
                subtracted_integrand_array(
                    common_batch, q_batch, False, "3d_cff_expansion"
                )[0]
            )
        )
        p_ose_1 = float(ose_bad_channel_probability(u, v, common_batch, 1.0)[0])
        p_ose_3 = float(ose_bad_channel_probability(u, v, common_batch, 3.0)[0])
        p_ose_4 = float(ose_bad_channel_probability(u, v, common_batch, 4.0)[0])
        p_jac_linear = float(
            jacobian_bad_channel_probability(u, v, common_batch, "linear")[0]
        )
        p_jac_log = float(
            jacobian_bad_channel_probability(u, v, common_batch, "log")[0]
        )
        rows.append(
            ChannelPartitionRow(
                path=path,
                radius=float(np.linalg.norm(common)),
                relative_radius=float(np.linalg.norm(q)),
                raw_abs_weight=raw_abs_weight,
                bad_channel_probability_ose_alpha1=p_ose_1,
                bad_channel_probability_ose_alpha3=p_ose_3,
                bad_channel_probability_ose_alpha4=p_ose_4,
                bad_channel_probability_jacobian_linear=p_jac_linear,
                bad_channel_probability_jacobian_log=p_jac_log,
                bad_channel_weighted_abs_ose_alpha1=raw_abs_weight * p_ose_1,
                bad_channel_weighted_abs_ose_alpha3=raw_abs_weight * p_ose_3,
                bad_channel_weighted_abs_ose_alpha4=raw_abs_weight * p_ose_4,
                bad_channel_weighted_abs_jacobian_linear=raw_abs_weight * p_jac_linear,
                bad_channel_weighted_abs_jacobian_log=raw_abs_weight * p_jac_log,
            )
        )
    return rows


def summarize_channel_partitions(
    rows: list[ChannelPartitionRow],
) -> list[dict[str, float | str]]:
    summary: list[dict[str, float | str]] = []
    paths = sorted({row.path for row in rows})
    observables = [
        "raw_abs_weight",
        "bad_channel_weighted_abs_ose_alpha1",
        "bad_channel_weighted_abs_ose_alpha3",
        "bad_channel_weighted_abs_ose_alpha4",
        "bad_channel_weighted_abs_jacobian_linear",
    ]

    for path in paths:
        selected = [row for row in rows if row.path == path]
        if path == "shrinking relative q":
            x = np.array([row.relative_radius for row in selected])
            min_x = 1.0e-4
            fit_label = "relative_radius"
        else:
            x = np.array([row.radius for row in selected])
            min_x = 1.0e5
            fit_label = "radius"
        for observable in observables:
            y = np.array([float(getattr(row, observable)) for row in selected])
            summary.append(
                {
                    "path": path,
                    "observable": observable,
                    "fit_variable": fit_label,
                    "tail_slope": fit_slope(x, y, min_x),
                }
            )
    return summary


def write_and_plot_channel_partitions(rows: list[ChannelPartitionRow]) -> None:
    summary = summarize_channel_partitions(rows)
    row_dicts = [row._asdict() for row in rows]
    (DATA / "channel_partition_scan.json").write_text(
        json.dumps({"rows": row_dicts, "summary": summary}, indent=2)
    )
    with (DATA / "channel_partition_scan.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row_dicts[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(row_dicts)
    (DATA / "channel_partition_summary.json").write_text(json.dumps(summary, indent=2))

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6))
    plot_specs = [
        ("fixed relative q", "radius", "common radius R"),
        ("fixed q/R", "radius", "common radius R"),
        ("shrinking relative q", "relative_radius", r"relative radius $|q|$"),
    ]
    series = [
        ("raw_abs_weight", "no partition", "#4d4d4d"),
        ("bad_channel_weighted_abs_ose_alpha3", "OSE alpha=3", "#2f7fbf"),
        ("bad_channel_weighted_abs_ose_alpha4", "OSE alpha=4", "#42a45f"),
        (
            "bad_channel_weighted_abs_jacobian_linear",
            "inverse J, linear map",
            "#ba4a00",
        ),
    ]
    for axis, (path, x_name, x_label) in zip(axes, plot_specs, strict=True):
        selected = [row for row in rows if row.path == path]
        x = np.array([float(getattr(row, x_name)) for row in selected])
        for y_name, label, color in series:
            y = np.array([float(getattr(row, y_name)) for row in selected])
            safe_y = np.maximum(y, np.finfo(float).tiny)
            axis.loglog(x, safe_y, label=label, color=color)
        axis.set_title(path)
        axis.set_xlabel(x_label)
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=8)
    axes[0].set_ylabel("bad direct-channel contribution")
    fig.suptitle("Current OSE partition versus inverse-Jacobian channel balance")
    fig.tight_layout()
    fig.savefig(FIGURES / "mock_gl00_proxy_channel_partitions.png", dpi=180)
    plt.close(fig)


def main() -> None:
    rows = (
        radial_scan(False, "3d_cff_expansion")
        + radial_scan(True, "3d_cff_expansion")
        + radial_scan(False, "4d_first_then_cff")
        + radial_scan(True, "4d_first_then_cff")
    )
    summary = summarize_radial(rows)
    write_radial_outputs(rows, summary)
    plot_radial(rows, summary)
    plot_components(rows)

    ridge_rows = (
        relative_ridge_scan(False, "3d_cff_expansion")
        + relative_ridge_scan(True, "3d_cff_expansion")
        + relative_ridge_scan(False, "4d_first_then_cff")
        + relative_ridge_scan(True, "4d_first_then_cff")
    )
    ridge_summary = write_and_plot_relative_ridge(ridge_rows)
    mc_results = run_mc_variance_experiment()
    write_and_plot_mc_variance(mc_results)
    channel_partition_rows = build_channel_partition_rows()
    write_and_plot_channel_partitions(channel_partition_rows)

    combined_summary = {
        "model": {
            "m_loop": M_LOOP,
            "m_rel_phys": M_REL_PHYS,
            "m_rel_uv": M_REL_UV,
            "m_rel_uv_4d": M_REL_UV_4D,
            "q_shift": Q_SHIFT.tolist(),
            "description": "3D mode expands CFF on-shell energy directly; 4D-first mode uses a massive vacuum relative propagator before CFF.",
        },
        "radial_slopes": summary,
        "relative_ridge": ridge_summary,
        "mc_variance": [mc_result_to_dict(result) for result in mc_results],
        "channel_partition": {
            "summary": summarize_channel_partitions(channel_partition_rows),
            "description": "Two-channel comparison between the bad direct (u,v) interpretation and a relative-momentum (P,q) interpretation.",
        },
        "figures": [
            str(FIGURES / "mock_gl00_proxy_radial.png"),
            str(FIGURES / "mock_gl00_proxy_components.png"),
            str(FIGURES / "mock_gl00_proxy_relative_ridge.png"),
            str(FIGURES / "mock_gl00_proxy_mc_variance.png"),
            str(FIGURES / "mock_gl00_proxy_channel_partitions.png"),
        ],
    }
    (DATA / "combined_summary.json").write_text(json.dumps(combined_summary, indent=2))

    print(json.dumps(combined_summary, indent=2))


if __name__ == "__main__":
    main()
