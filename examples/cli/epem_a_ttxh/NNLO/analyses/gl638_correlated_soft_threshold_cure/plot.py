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

    shared_ct_filter = r"^(?=.*η\(8,12,14\))(?=.*(?:shared_1l|default)).*$"
    duplicated_single_filter = r"^threshold:component c(?:8|9|10|11) "
    duplicated_pair_filter = r"^threshold:component c(?:30|31|32|33|38|39|40|41) "
    common_axis = [
        "--x-log-scale",
        "--branch-layout",
        "signed",
        "--x-range",
        "-0.01",
        "1e-2",
        "--x-symlog-linthresh",
        "1e-6",
    ]
    common_log = [*common_axis, "--fit-log-slope", "--fit-points", "12"]
    with tempfile.TemporaryDirectory(prefix="gl638_correlated_plots_") as temporary:
        temp = Path(temporary)
        sidecar = temp / "kinematics.json"
        kinematics_pdf = temp / "kinematics.pdf"
        comparison_pdf = temp / "comparison.pdf"
        cut_cancellation_pdf = temp / "cut_cancellation.pdf"
        ct_cancellation_pdf = temp / "ct_cancellation.pdf"
        shared_components_pdf = temp / "shared_components.pdf"
        duplicated_singles_pdf = temp / "duplicated_singles.pdf"
        duplicated_pairs_pdf = temp / "duplicated_pairs.pdf"
        write_sidecar(ir_safe, sidecar)

        run_plotter(
            str(ir_safe),
            "--series-json",
            str(sidecar),
            "--include-contribution",
            r"^series:",
            "--output",
            str(kinematics_pdf),
            "--component",
            "abs",
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
            "--include-contribution",
            r"^total_weight$",
            "--combine-plots",
            "--component",
            "abs",
            *common_log,
            "--title",
            "GL638: correlated soft + both-threshold cure",
        )
        run_plotter(
            str(ir_safe),
            "--output",
            str(cut_cancellation_pdf),
            "--include-contribution",
            r"^total_weight$",
            "--include-contribution",
            r"^contribution:event_weight GL638 cut(?:0|1|3) ori0$",
            "--component",
            "imag",
            *common_log,
            "--y-label",
            "|imaginary cut weight|",
            "--title",
            "GL638: leading LU cut cancellation",
        )
        run_plotter(
            str(ir_safe),
            "--output",
            str(ct_cancellation_pdf),
            "--include-contribution",
            r"^additional:threshold_counterterm_0$",
            "--include-contribution",
            r"^contribution:threshold_counterterm_0 GL638 cut(?:0|1|3) ori0$",
            "--component",
            "imag",
            *common_log,
            "--y-label",
            "|imaginary threshold CT weight|",
            "--title",
            "GL638: cross-cut threshold-counterterm cancellation",
        )
        run_plotter(
            str(ir_safe),
            "--output",
            str(shared_components_pdf),
            "--threshold-report",
            "--threshold-report-section",
            "singles",
            "--threshold-quantity",
            "weighted",
            "--include-contribution",
            shared_ct_filter,
            "--facets-per-page",
            "3",
            "--component",
            "abs",
            *common_log,
            "--title",
            "GL638: shared one-loop threshold CTs across cancelling cuts",
        )
        run_plotter(
            str(ir_safe),
            "--output",
            str(duplicated_singles_pdf),
            "--threshold-report",
            "--threshold-report-section",
            "singles",
            "--threshold-quantity",
            "bare",
            "--include-contribution",
            duplicated_single_filter,
            "--facets-per-page",
            "2",
            "--component",
            "abs",
            "--hide-info-box",
            "--max-gap-warnings",
            "1",
            *common_log,
            "--title",
            "GL638: duplicated η(7,8) variants before multiplier weighting",
        )
        run_plotter(
            str(ir_safe),
            "--output",
            str(duplicated_pairs_pdf),
            "--threshold-report",
            "--threshold-report-section",
            "pairs",
            "--threshold-quantity",
            "bare",
            "--include-contribution",
            duplicated_pair_filter,
            "--component",
            "abs",
            "--hide-info-box",
            "--max-gap-warnings",
            "1",
            *common_axis,
            "--title",
            "GL638: duplicated variants paired with η(5,10) before weighting",
        )

        output = HERE / "correlated_soft_threshold_cure.pdf"
        merge_pdfs(
            output,
            [
                kinematics_pdf,
                comparison_pdf,
                cut_cancellation_pdf,
                ct_cancellation_pdf,
                shared_components_pdf,
                duplicated_singles_pdf,
                duplicated_pairs_pdf,
            ],
        )
    print(f"created {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
