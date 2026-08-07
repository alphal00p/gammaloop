#!/usr/bin/env python3
"""Build the double-sided GL638 correlated soft/threshold cure report."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import fitz


HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[6]
PLOTTER = ROOT / "assets" / "plot_approach_result.py"
TOP_MASS = 173.0


def run_plotter(*arguments: str) -> None:
    subprocess.run(
        [sys.executable, str(PLOTTER), *arguments],
        cwd=ROOT,
        check=True,
    )


def merge_pdfs(output: Path, inputs: list[Path]) -> None:
    merged = fitz.open()
    try:
        for path in inputs:
            with fitz.open(path) as source:
                merged.insert_pdf(source)
        merged.save(output)
    finally:
        merged.close()


def subtract(left: list[float], right: list[float]) -> list[float]:
    return [a - b for a, b in zip(left, right)]


def norm(vector: list[float]) -> float:
    return math.sqrt(sum(component * component for component in vector))


def energy(vector: list[float]) -> float:
    return math.sqrt(norm(vector) ** 2 + TOP_MASS**2)


def write_sidecar(result_path: Path, output_path: Path) -> None:
    result = json.loads(result_path.read_text())
    if result.get("schema_version") != 3:
        raise ValueError(f"{result_path} is not fresh approach schema v3 output")
    series: dict[str, list[dict[str, float | int]]] = {}
    for record in result["points"]:
        if record.get("status") != "evaluated":
            continue
        point = record["point"]
        k0, k3 = point[:3], point[9:12]
        q13 = subtract(k0, k3)
        values = {
            "|q13|": norm(q13),
            "|eta(5,10)|": abs(2.0 * energy(k3) - 1000.0),
            "|eta(5,12,13)|": abs(energy(k0) + energy(k3) + norm(q13) - 1000.0),
        }
        for label, value in values.items():
            series.setdefault(label, []).append(
                {"point_index": record["index"], "value": value}
            )
    output_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "series": [
                    {"label": label, "axis_index": 0, "values": values}
                    for label, values in series.items()
                ],
            },
            indent=2,
        )
    )


def main() -> None:
    threshold_off = HERE / "threshold_off.json"
    legacy = HERE / "legacy.json"
    ir_safe = HERE / "ir_safe.json"
    for path in (threshold_off, legacy, ir_safe):
        if not path.is_file():
            raise FileNotFoundError(f"Run the analysis card first; missing {path}")

    variant_filter = r"^(?=.*c\(2,4,12\))(?=.*(?:intrinsic_1l|embedded_2l)).*$"
    multiplier_filter = (
        r"^(?=.*c\(2,4,12\))(?=.*(?:intrinsic_1l|embedded_2l))"
        r"(?=.*η\(5,10\)).*$"
    )
    common_view = [
        "--component",
        "abs",
        "--x-log-scale",
        "--branch-layout",
        "split",
        "--x-range",
        "1e-6",
        "1e-2",
    ]
    common_log = [*common_view, "--fit-log-slope", "--fit-points", "4"]
    with tempfile.TemporaryDirectory(prefix="gl638_correlated_plots_") as temporary:
        temp = Path(temporary)
        sidecar = temp / "kinematics.json"
        kinematics_pdf = temp / "kinematics.pdf"
        comparison_pdf = temp / "comparison.pdf"
        summary_pdf = temp / "summary.pdf"
        components_pdf = temp / "components.pdf"
        multipliers_pdf = temp / "multipliers.pdf"
        write_sidecar(ir_safe, sidecar)

        run_plotter(
            str(ir_safe),
            "--series-json",
            str(sidecar),
            "--include-contribution",
            r"^series:",
            "--output",
            str(kinematics_pdf),
            *common_log,
            "--y-label",
            "kinematic distance [GeV]",
            "--title",
            "GL638: simultaneous q13-soft and two-threshold approach",
        )
        run_plotter(
            str(threshold_off),
            str(legacy),
            str(ir_safe),
            "--output",
            str(comparison_pdf),
            "--result-label",
            "threshold subtraction off",
            "--result-label",
            "legacy maximal subspaces",
            "--result-label",
            "IR-safe directives",
            "--combine-plots",
            *common_log,
            "--title",
            "GL638: correlated soft + both-threshold cure",
        )
        run_plotter(
            str(ir_safe),
            "--output",
            str(summary_pdf),
            "--threshold-report",
            "--threshold-report-section",
            "summary",
            "--include-contribution",
            variant_filter,
            *common_view,
            "--title",
            "GL638 cut (2,4,12): correlated cure decomposition",
        )
        run_plotter(
            str(ir_safe),
            "--output",
            str(components_pdf),
            "--threshold-report",
            "--threshold-report-section",
            "singles",
            "--threshold-report-section",
            "pairs",
            "--threshold-quantity",
            "weighted",
            "--threshold-quantity",
            "bare",
            "--include-contribution",
            variant_filter,
            "--facets-per-page",
            "2",
            *common_log,
            "--title",
            "GL638 cut (2,4,12): correlated cure decomposition",
        )
        run_plotter(
            str(ir_safe),
            "--output",
            str(multipliers_pdf),
            "--threshold-report",
            "--threshold-report-section",
            "multipliers",
            "--include-contribution",
            multiplier_filter,
            "--x-log-scale",
            "--branch-layout",
            "split",
            "--x-range",
            "1e-6",
            "1e-2",
            "--title",
            "GL638 cut (2,4,12): correlated multiplier contexts",
        )

        output = HERE / "correlated_soft_threshold_cure.pdf"
        merge_pdfs(
            output,
            [
                kinematics_pdf,
                comparison_pdf,
                summary_pdf,
                components_pdf,
                multipliers_pdf,
            ],
        )
    print(f"created {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
