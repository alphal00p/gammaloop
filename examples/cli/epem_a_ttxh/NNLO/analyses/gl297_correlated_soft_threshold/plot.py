#!/usr/bin/env python3
"""Plot the GL297 bare, legacy, and forced correlated approaches."""

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


def plot(direction: str, soft_edge: str) -> None:
    root = repository_root()
    command = [
        sys.executable,
        str(root / "assets" / "plot_approach_result.py"),
        *(
            str(ANALYSIS_DIR / f"{direction}_{case}.json")
            for case in ("bare", "legacy", "forced")
        ),
        "--output",
        str(ANALYSIS_DIR / f"{direction}_cure.pdf"),
        "--combine-plots",
        "--component",
        "real",
        "--t-branch",
        "positive",
        "--x-log-scale",
        "--y-log-scale",
        "--fit-log-slope",
        "--hide-info-box",
        "--x-range",
        "0.0001",
        "1",
        "--title",
        f"GL297: correlated soft/threshold limit for {soft_edge}",
    ]
    subprocess.run(command, cwd=root, check=True)


def main() -> None:
    plot("e12", "e12")
    plot("e13", "e13")


if __name__ == "__main__":
    main()
