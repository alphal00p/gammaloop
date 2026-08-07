#!/usr/bin/env python3
"""Build the finite-path GL638 multiplier/decomposition report."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import fitz


HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[6]
PLOTTER = ROOT / "assets" / "plot_approach_result.py"
RESULT = HERE / "finite_multiplier_path.json"


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


def main() -> None:
    if not RESULT.is_file():
        raise FileNotFoundError(f"Run the analysis card first; missing {RESULT}")
    if json.loads(RESULT.read_text()).get("schema_version") != 3:
        raise ValueError(f"{RESULT} is not fresh approach schema v3 output")

    multiplier_filter = r"^threshold:component c(?:30|38) "
    common_axis = [
        "--x-log-scale",
        "--branch-layout",
        "signed",
        "--x-range",
        "-1",
        "1",
        "--x-symlog-linthresh",
        "1e-4",
    ]
    with tempfile.TemporaryDirectory(prefix="gl638_multiplier_plots_") as temporary:
        temp = Path(temporary)
        overview_pdf = temp / "overview.pdf"
        multipliers_pdf = temp / "multipliers.pdf"

        run_plotter(
            str(RESULT),
            "--output",
            str(overview_pdf),
            "--component",
            "real",
            "--linear-y-scale",
            *common_axis,
            "--title",
            "GL638: finite multiplier-deformation total",
        )
        run_plotter(
            str(RESULT),
            "--output",
            str(multipliers_pdf),
            "--threshold-report",
            "--threshold-report-section",
            "multipliers",
            "--include-contribution",
            multiplier_filter,
            "--facets-per-page",
            "2",
            *common_axis,
            "--title",
            "GL638 cut (2,4,12): star-root selector for local-local pairs",
        )

        output = HERE / "multiplier_decomposition.pdf"
        merge_pdfs(output, [overview_pdf, multipliers_pdf])
    print(f"created {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
