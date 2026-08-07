#!/usr/bin/env python3
"""Build the GL638 generic-component and correlated-cure PDF reports."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import fitz
import matplotlib.pyplot as plt


ANALYSIS_DIR = Path(__file__).resolve().parent
TOP_MASS = 173.0
LAMBDA_PER_T = 1.0e-2


def repository_root() -> Path:
    for parent in ANALYSIS_DIR.parents:
        if (parent / "assets" / "plot_approach_result.py").is_file():
            return parent
    raise RuntimeError("could not locate assets/plot_approach_result.py")


def run_threshold_report(
    root: Path,
    result: Path,
    output: Path,
    title: str,
    minimum_t: float,
) -> None:
    subprocess.run(
        [
            sys.executable,
            str(root / "assets" / "plot_approach_result.py"),
            str(result),
            "--output",
            str(output),
            "--threshold-report",
            "--component",
            "real",
            "--t-branch",
            "both",
            "--x-log-scale",
            "--branch-layout",
            "auto",
            "--hide-info-box",
            "--x-range",
            str(minimum_t),
            "1",
            "--title",
            title,
        ],
        cwd=root,
        check=True,
    )


def log_slope(x_values: list[float], y_values: list[float]) -> float:
    points = [
        (math.log(x), math.log(y))
        for x, y in zip(x_values, y_values)
        if x > 0.0 and y > 0.0 and math.isfinite(x) and math.isfinite(y)
    ]
    mean_x = sum(x for x, _ in points) / len(points)
    mean_y = sum(y for _, y in points) / len(points)
    denominator = sum((x - mean_x) ** 2 for x, _ in points)
    return sum((x - mean_x) * (y - mean_y) for x, y in points) / denominator


def kinematic_distances(point: list[float]) -> tuple[float, float, float]:
    k0 = point[0:3]
    k3 = point[9:12]
    q13 = [left - right for left, right in zip(k0, k3)]

    def norm(vector: list[float]) -> float:
        return math.sqrt(sum(component**2 for component in vector))

    energy0 = math.sqrt(norm(k0) ** 2 + TOP_MASS**2)
    energy3 = math.sqrt(norm(k3) ** 2 + TOP_MASS**2)
    return (
        norm(q13),
        abs(2.0 * energy3 - 1000.0),
        abs(energy0 + energy3 + norm(q13) - 1000.0),
    )


def draw_kinematic_page(result_path: Path, output: Path) -> None:
    data = json.loads(result_path.read_text())
    branches: dict[str, list[tuple[float, tuple[float, float, float]]]] = {
        "negative": [],
        "positive": [],
    }
    for point in data["points"]:
        if point.get("status") != "evaluated" or point["t"] == 0.0:
            continue
        branch = "negative" if point["t"] < 0.0 else "positive"
        branches[branch].append(
            (abs(float(point["t"]) * LAMBDA_PER_T), kinematic_distances(point["point"]))
        )

    labels = (
        r"$|q_{13}|$",
        r"$|\eta(\mathrm{eset}(5,10))|$",
        r"$|\eta(\mathrm{eset}(5,12,13))|$",
    )
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 8.2), sharey=True)
    for axis, branch in zip(axes, ("negative", "positive")):
        samples = sorted(branches[branch])
        lambdas = [sample[0] for sample in samples]
        for quantity, label in enumerate(labels):
            values = [sample[1][quantity] for sample in samples]
            slope = log_slope(lambdas, values)
            axis.plot(
                lambdas,
                values,
                linewidth=1.7,
                marker="o",
                label=f"{label}  p={slope:.3f}",
            )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.grid(True, which="major", color="0.86")
        axis.grid(True, which="minor", color="0.94")
        axis.set_xlabel(
            r"$|\lambda|$"
            + (r"  ($\lambda<0$)" if branch == "negative" else r"  ($\lambda>0$)")
        )
        axis.set_title(
            r"approach from $\lambda<0$"
            if branch == "negative"
            else r"approach from $\lambda>0$"
        )
        axis.legend(loc="best", fontsize=8.5, frameon=True)
    axes[0].set_xlim(max(branches["negative"])[0], min(branches["negative"])[0])
    axes[1].set_xlim(min(branches["positive"])[0], max(branches["positive"])[0])
    axes[0].set_ylabel("kinematic distance [GeV]")
    fig.suptitle(
        "GL638 cut (2,4,12): common q13 soft and two-threshold approach",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.035,
        r"$K_0=K^*+\lambda a$,  $K_3=K^*+\lambda(a+d)$,  "
        r"$q_{13}=-\lambda d$; all three fitted powers are unity.",
        ha="center",
        fontsize=9.5,
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.9, bottom=0.13, wspace=0.06)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def merge_pdfs(inputs: list[Path], output: Path) -> None:
    merged = fitz.open()
    try:
        for input_path in inputs:
            with fitz.open(input_path) as document:
                merged.insert_pdf(document)
        merged.save(output)
    finally:
        merged.close()


def main() -> None:
    root = repository_root()
    with tempfile.TemporaryDirectory(prefix="gl638_approach_") as temporary:
        temporary_dir = Path(temporary)

        run_threshold_report(
            root,
            ANALYSIS_DIR / "component_approach.json",
            ANALYSIS_DIR / "multiplier_decomposition.pdf",
            "GL638: expanded threshold-counterterm components",
            1.0e-3,
        )

        correlated_json = ANALYSIS_DIR / "q13_correlated.json"
        kinematics = temporary_dir / "q13_kinematics.pdf"
        correlated_report = temporary_dir / "q13_threshold_report.pdf"
        draw_kinematic_page(correlated_json, kinematics)
        run_threshold_report(
            root,
            correlated_json,
            correlated_report,
            "GL638 cut (2,4,12): q13 correlated cure",
            1.0e-4,
        )
        merge_pdfs(
            [kinematics, correlated_report],
            ANALYSIS_DIR / "q13_correlated_cure.pdf",
        )


if __name__ == "__main__":
    main()
