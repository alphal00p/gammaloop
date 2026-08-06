#!/usr/bin/env python3
"""Plot the GL638 event-level threshold decomposition smoke trajectory."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ANALYSIS_DIR = Path(__file__).resolve().parent


def repository_root() -> Path:
    for parent in ANALYSIS_DIR.parents:
        if (parent / "assets" / "plot_approach_result.py").is_file():
            return parent
    raise RuntimeError("could not locate assets/plot_approach_result.py")


def main() -> None:
    root = repository_root()
    subprocess.run(
        [
            sys.executable,
            str(root / "assets" / "plot_approach_result.py"),
            str(ANALYSIS_DIR / "component_approach.json"),
            "--output",
            str(ANALYSIS_DIR / "multiplier_decomposition.pdf"),
            "--component",
            "abs",
            "--threshold-decomposition",
            "summary",
            "--exclude-contribution",
            "^threshold:decomposition_total$",
            "--t-branch",
            "positive",
            "--x-log-scale",
            "--y-log-scale",
            "--hide-info-box",
            "--x-range",
            "0.001",
            "1",
            "--title",
            "GL638: finite expanded threshold-counterterm decomposition",
        ],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
