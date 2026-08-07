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

    variant_filter = r"^(?=.*c\(2,4,12\))(?=.*(?:intrinsic_1l|embedded_2l)).*$"
    multiplier_filter = (
        r"^(?=.*c\(2,4,12\))(?=.*(?:intrinsic_1l|embedded_2l))"
        r"(?=.*η\(5,10\)).*$"
    )
    with tempfile.TemporaryDirectory(prefix="gl638_multiplier_plots_") as temporary:
        temp = Path(temporary)
        overview_pdf = temp / "overview.pdf"
        components_pdf = temp / "components.pdf"
        multipliers_pdf = temp / "multipliers.pdf"

        run_plotter(
            str(RESULT),
            "--output",
            str(overview_pdf),
            "--component",
            "real",
            "--linear-y-scale",
            "--title",
            "GL638: finite multiplier-deformation total",
        )
        run_plotter(
            str(RESULT),
            "--output",
            str(components_pdf),
            "--threshold-report",
            "--threshold-report-section",
            "summary",
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
            "--component",
            "abs",
            "--title",
            "GL638 cut (2,4,12): finite expanded-component decomposition",
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
            "--title",
            "GL638 cut (2,4,12): occurrence-resolved multipliers",
        )

        output = HERE / "multiplier_decomposition.pdf"
        merge_pdfs(output, [overview_pdf, components_pdf, multipliers_pdf])
    print(f"created {output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
