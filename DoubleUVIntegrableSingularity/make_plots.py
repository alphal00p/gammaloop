#!/usr/bin/env python3
"""Generate toy-model figures for the GL00 double-UV report."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
FIGURES = ROOT / "figures"
DATA = ROOT / "data"


def toy_weight(r: np.ndarray, epsilon: float, q0: float = 1.0) -> np.ndarray:
    """Integrator-weight toy model after the radial jacobians are included."""

    return r**2 / np.sqrt((epsilon * r) ** 2 + q0**2)


def make_toy_radial_scaling() -> None:
    r = np.logspace(0.0, 5.0, 400)
    epsilons = [1.0, 0.1, 0.03, 0.01, 0.0]

    plt.figure(figsize=(6.4, 4.3))
    for eps in epsilons:
        label = "epsilon = 0" if eps == 0.0 else f"epsilon = {eps:g}"
        plt.loglog(r, toy_weight(r, eps), label=label)

    ref = toy_weight(np.array([1.0e3]), 0.03)[0]
    plt.loglog(r, ref * (r / 1.0e3), "k--", lw=1.0, label="slope 1")
    plt.loglog(r, ref * (r / 1.0e3) ** 2, "k:", lw=1.0, label="slope 2")
    plt.xlabel(r"common UV radius $R$")
    plt.ylabel(r"toy integrator weight $W_{\rm toy}$")
    plt.title(r"Radial growth from a small relative edge $q=v-u+Q$")
    plt.legend(fontsize=8, ncol=2)
    plt.grid(True, which="both", alpha=0.25)
    plt.tight_layout()
    plt.savefig(FIGURES / "toy_radial_scaling.png", dpi=180)
    plt.close()


def make_toy_relative_ridge() -> None:
    r = 1000.0
    q0 = 1.0
    rho = np.linspace(0.92, 1.08, 320)
    theta = np.linspace(0.0, 0.14, 280)
    rho_grid, theta_grid = np.meshgrid(rho, theta)
    qhat = np.sqrt(
        1.0 + rho_grid**2 - 2.0 * rho_grid * np.cos(theta_grid)
    )
    weight_over_r = r / np.sqrt((r * qhat) ** 2 + q0**2)

    plt.figure(figsize=(6.4, 4.5))
    levels = np.linspace(-2.0, 3.0, 26)
    contour = plt.contourf(
        rho_grid,
        theta_grid,
        np.log10(weight_over_r),
        levels=levels,
        cmap="viridis",
        extend="both",
    )
    plt.colorbar(contour, label=r"$\log_{10}(W_{\rm toy}/R)$")
    plt.xlabel(r"radius ratio $\rho=|v|/|u|$")
    plt.ylabel(r"opening angle $\theta(u,v)$")
    plt.title(r"Relative-momentum ridge at equal radius and collinear UV rays")
    plt.tight_layout()
    plt.savefig(FIGURES / "toy_relative_ridge.png", dpi=180)
    plt.close()


def make_gl00_compact_fit_plot() -> None:
    path = DATA / "gl00_uv_ray_weight.tsv"
    rows = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append({key: float(value) for key, value in row.items()})

    lam = np.array([row["lambda"] for row in rows])
    off = np.array([row["abs_lmb_off"] for row in rows])
    on = np.array([row["abs_lmb_on"] for row in rows])

    mask = lam >= 1.0
    off_slope, off_intercept = np.polyfit(np.log(lam[mask]), np.log(off[mask]), 1)
    on_slope, on_intercept = np.polyfit(np.log(lam[mask]), np.log(on[mask]), 1)

    plt.figure(figsize=(6.4, 4.3))
    plt.loglog(lam, off, label=f"LMB MC off, slope {off_slope:.3f}")
    plt.loglog(lam, on, label=f"LMB MC on, slope {on_slope:.3f}")
    plt.loglog(
        lam[mask],
        np.exp(off_intercept) * lam[mask] ** off_slope,
        "k--",
        lw=1.0,
        label="large-radius fit",
    )
    plt.xlabel(r"radial scale multiplier $\lambda$")
    plt.ylabel(r"$|$GammaLoop integrator weight$|$")
    plt.title("GL00 approach along the high-weight double-UV ray")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(FIGURES / "gl00_uv_ray_fit_compact.png", dpi=180)
    plt.close()


def write_toy_data() -> None:
    path = DATA / "toy_radial_scaling.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["R", "epsilon", "W_toy"])
        for eps in [1.0, 0.1, 0.03, 0.01, 0.0]:
            for r in np.logspace(0.0, 5.0, 80):
                writer.writerow([f"{r:.16e}", f"{eps:.16e}", f"{toy_weight(np.array([r]), eps)[0]:.16e}"])


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    make_toy_radial_scaling()
    make_toy_relative_ridge()
    make_gl00_compact_fit_plot()
    write_toy_data()


if __name__ == "__main__":
    main()
