#!/usr/bin/env python3
"""Build the curated GL297 approach PDFs from fresh schema-v3 CLI output."""

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


def vector_add(*vectors: list[float]) -> list[float]:
    return [sum(components) for components in zip(*vectors)]


def vector_scale(vector: list[float], factor: float) -> list[float]:
    return [factor * component for component in vector]


def norm(vector: list[float]) -> float:
    return math.sqrt(sum(component * component for component in vector))


def energy(vector: list[float]) -> float:
    return math.sqrt(norm(vector) ** 2 + TOP_MASS**2)


def kinematic_values(direction: str, point: list[float]) -> dict[str, float]:
    k0, k1, k2, k3 = [point[index : index + 3] for index in range(0, 12, 3)]
    if direction == "e12":
        q = vector_add(k0, k1, vector_scale(k2, -1.0), vector_scale(k3, -1.0))
        p3 = vector_add(k0, k1, vector_scale(k2, -1.0))
        return {
            "|q12|": norm(q),
            "|eta(3,9)|": abs(2.0 * energy(p3) - 1000.0),
            "|eta(3,11,12)|": abs(energy(p3) + energy(k3) + norm(q) - 1000.0),
        }
    q = vector_add(k3, vector_scale(k1, -1.0))
    return {
        "|q13|": norm(q),
        "|eta(5,7)|": abs(2.0 * energy(k1) - 1000.0),
        "|eta(5,11,13)|": abs(energy(k1) + energy(k3) + norm(q) - 1000.0),
    }


def write_sidecar(direction: str, result_path: Path, output_path: Path) -> None:
    result = json.loads(result_path.read_text())
    if result.get("schema_version") != 3:
        raise ValueError(f"{result_path} is not fresh approach schema v3 output")
    series: dict[str, list[dict[str, float | int]]] = {}
    for record in result["points"]:
        if record.get("status") != "evaluated":
            continue
        for label, value in kinematic_values(direction, record["point"]).items():
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


def build_direction(direction: str) -> Path:
    threshold_off = HERE / f"{direction}_threshold_off.json"
    legacy = HERE / f"{direction}_legacy.json"
    forced = HERE / f"{direction}_forced.json"
    for path in (threshold_off, legacy, forced):
        if not path.is_file():
            raise FileNotFoundError(f"Run the analysis card first; missing {path}")

    cut = r"2,6,7" if direction == "e12" else r"2,4,9"
    component_filter = rf"^(?=.*c\({cut}\))(?=.*forced_1l).*$"
    title_name = direction.upper()
    with tempfile.TemporaryDirectory(prefix=f"gl297_{direction}_plots_") as temporary:
        temp = Path(temporary)
        sidecar = temp / "kinematics.json"
        kinematics_pdf = temp / "kinematics.pdf"
        comparison_pdf = temp / "comparison.pdf"
        summary_pdf = temp / "summary.pdf"
        components_pdf = temp / "components.pdf"
        write_sidecar(direction, forced, sidecar)

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
        run_plotter(
            str(forced),
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
            f"GL297 {title_name}: simultaneous soft and threshold distances",
        )
        run_plotter(
            str(threshold_off),
            str(legacy),
            str(forced),
            "--output",
            str(comparison_pdf),
            "--result-label",
            "threshold subtraction off",
            "--result-label",
            "legacy maximal subspace",
            "--result-label",
            "forced one-loop subspace",
            "--combine-plots",
            *common_log,
            "--title",
            f"GL297 {title_name}: correlated soft + threshold approach",
        )
        run_plotter(
            str(forced),
            "--output",
            str(summary_pdf),
            "--threshold-report",
            "--threshold-report-section",
            "summary",
            "--include-contribution",
            component_filter,
            *common_view,
            "--title",
            f"GL297 {title_name}: forced-subspace decomposition",
        )
        run_plotter(
            str(forced),
            "--output",
            str(components_pdf),
            "--threshold-report",
            "--threshold-report-section",
            "singles",
            "--threshold-report-section",
            "pairs",
            "--threshold-quantity",
            "weighted",
            "--include-contribution",
            component_filter,
            "--facets-per-page",
            "2",
            *common_log,
            "--title",
            f"GL297 {title_name}: forced-subspace decomposition",
        )

        output = HERE / f"{direction}_correlated_soft_threshold.pdf"
        merge_pdfs(
            output,
            [kinematics_pdf, comparison_pdf, summary_pdf, components_pdf],
        )
    return output


def main() -> None:
    for direction in ("e12", "e13"):
        print(f"created {build_direction(direction).relative_to(ROOT)}")


if __name__ == "__main__":
    main()
