#!/usr/bin/env python3
"""Build the GL297 correlated-cure comparison and component reports."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import fitz


ANALYSIS_DIR = Path(__file__).resolve().parent


def repository_root() -> Path:
    for parent in ANALYSIS_DIR.parents:
        if (parent / "assets" / "plot_approach_result.py").is_file():
            return parent
    raise RuntimeError("could not locate assets/plot_approach_result.py")


def run_plotter(root: Path, arguments: list[str]) -> None:
    subprocess.run(
        [sys.executable, str(root / "assets" / "plot_approach_result.py"), *arguments],
        cwd=root,
        check=True,
    )


def merge_pdfs(inputs: list[Path], output: Path) -> None:
    merged = fitz.open()
    try:
        for input_path in inputs:
            with fitz.open(input_path) as document:
                merged.insert_pdf(document)
        merged.save(output)
    finally:
        merged.close()


def plot(direction: str, soft_edge: str) -> None:
    root = repository_root()
    with tempfile.TemporaryDirectory(prefix=f"gl297_{direction}_") as temporary:
        temporary_dir = Path(temporary)
        overview = temporary_dir / "overview.pdf"
        components = temporary_dir / "components.pdf"
        common = [
            "--t-branch",
            "both",
            "--x-log-scale",
            "--branch-layout",
            "auto",
            "--y-log-scale",
            "--hide-info-box",
            "--x-range",
            "0.0001",
            "1",
        ]

        run_plotter(
            root,
            [
                *(
                    str(ANALYSIS_DIR / f"{direction}_{case}.json")
                    for case in ("threshold_off", "legacy", "forced")
                ),
                "--output",
                str(overview),
                "--result-label",
                "threshold CT off",
                "--result-label",
                "legacy maximal-subspace CT",
                "--result-label",
                "forced one-loop CT",
                "--combine-plots",
                "--component",
                "real",
                "--fit-log-slope",
                *common,
                "--title",
                f"GL297: correlated {soft_edge} soft/threshold cure",
            ],
        )
        run_plotter(
            root,
            [
                str(ANALYSIS_DIR / f"{direction}_forced.json"),
                "--output",
                str(components),
                "--threshold-report",
                "--component",
                "real",
                "--facets-per-page",
                "3",
                *common,
                "--title",
                f"GL297: forced one-loop components for {soft_edge}",
            ],
        )
        merge_pdfs([overview, components], ANALYSIS_DIR / f"{direction}_cure.pdf")


def main() -> None:
    plot("e12", "e12")
    plot("e13", "e13")


if __name__ == "__main__":
    main()
