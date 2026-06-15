#!/usr/bin/env python3
"""Inspect radial approaches around the current all-graph max-weight points.

The script intentionally runs against a copied, read-only state so it does not
touch any long-running integration workspace.  It builds a temporary command
card, evaluates the selected graph bins with the same sampling settings as the
all-graph run, and writes compact TSV/JSON/PNG diagnostics.
"""

from __future__ import annotations

import json
import math
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent
STATE = Path("/tmp/gammaloop_aa_aa_2l_maxweight_inspect_state")
CARD = Path("/tmp/aa_aa_2l_current_maxweight_scan.toml")
LOG = Path("/tmp/aa_aa_2l_current_maxweight_scan.log")
DATA = ROOT / "data" / "current_maxweight_diagnostics"
FIGURES = ROOT / "figures"
GAMMALOOP = ROOT.parent / "gammaloop"


@dataclass(frozen=True)
class ScanPoint:
    label: str
    graph: int
    xs: tuple[float, float, float, float, float, float]


@dataclass(frozen=True)
class ScanMode:
    label: str
    sampling: str


MODES = [
    ScanMode(
        "no_lmb",
        '''set process -p aa_aa string '
        [sampling]
        graphs = "monte_carlo"
        orientations = "summed"
        lmb_multichanneling = false
        lmb_channels = "summed"
        coordinate_system = "spherical"
        mapping = "linear"
        b = 1.0
    ' ''',
    ),
    ScanMode(
        "ose",
        '''set process -p aa_aa string '
        [sampling]
        graphs = "monte_carlo"
        orientations = "summed"
        lmb_multichanneling = true
        lmb_channels = "summed"
        lmb_channel_weight = "ose"
        coordinate_system = "spherical"
        mapping = "linear"
        b = 1.0
    ' ''',
    ),
    ScanMode(
        "inverse_jacobian",
        '''set process -p aa_aa string '
        [sampling]
        graphs = "monte_carlo"
        orientations = "summed"
        lmb_multichanneling = true
        lmb_channels = "summed"
        lmb_channel_weight = "inverse_jacobian"
        coordinate_system = "spherical"
        mapping = "linear"
        b = 1.0
    ' ''',
    ),
    ScanMode(
        "inverse_jacobian_log",
        '''set process -p aa_aa string '
        [sampling]
        graphs = "monte_carlo"
        orientations = "summed"
        lmb_multichanneling = true
        lmb_channels = "summed"
        lmb_channel_weight = "inverse_jacobian"
        coordinate_system = "spherical"
        mapping = "log"
        b = 1.0
    ' ''',
    ),
]


POINTS = [
    ScanPoint(
        "GL00_im_plus",
        0,
        (
            0.9993972417847613,
            0.40776395552741984,
            0.5435242330138342,
            0.9999052822352779,
            0.4204277904161551,
            0.3417750999292709,
        ),
    ),
    ScanPoint(
        "GL03_re_minus_im_plus",
        2,
        (
            0.99985130679222345,
            0.7524997097227335,
            0.36093249964439938,
            0.99996982475085883,
            0.6721874876322349,
            0.6252914825413105,
        ),
    ),
    ScanPoint(
        "GL05_re_plus",
        3,
        (
            0.9992159405808112,
            0.6868444838619521,
            0.6429740315399803,
            0.9990699795198853,
            0.11418044330903651,
            0.3953319419540728,
        ),
    ),
    ScanPoint(
        "GL06_re_minus",
        4,
        (
            0.9995215227215548,
            0.5787928959151282,
            0.22059957686624845,
            0.9999299728315765,
            0.8989299227500787,
            0.42030419319463375,
        ),
    ),
    ScanPoint(
        "GL18_bins128_re_minus_im_plus",
        10,
        (
            0.99972605511944601,
            0.59527066218796598,
            0.028121655386487698,
            0.99929042279752922,
            0.99432932400174268,
            0.86220910818747232,
        ),
    ),
    ScanPoint(
        "GL14_bins128_im_minus",
        7,
        (
            0.99928560786268694,
            0.69209294959349577,
            0.28659450247430496,
            0.99985395518723508,
            0.11911699688086141,
            0.050153678028659998,
        ),
    ),
    ScanPoint(
        "GL06_prob1000_re_plus_im_plus",
        4,
        (
            0.99988492079262081,
            0.44226073563398893,
            0.29422178156572731,
            0.99995583397035226,
            0.25056458851596020,
            0.069123913846451557,
        ),
    ),
    ScanPoint(
        "GL17_im_minus",
        9,
        (
            0.9997425961322831,
            0.1305730414150327,
            0.8833503679259446,
            0.9999938714344689,
            0.44900669230836104,
            0.9348993796532749,
        ),
    ),
]

SCALES = (8.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.125, 0.0625)


def rescale_radial(xs: tuple[float, ...], scale: float) -> tuple[float, ...]:
    out = list(xs)
    for i in (0, 3):
        out[i] = 1.0 - scale * (1.0 - xs[i])
    return tuple(out)


def mapped_radius(x: float) -> float:
    return x / (1.0 - x)


def build_card() -> dict[str, tuple[str, float, ScanPoint, tuple[float, ...]]]:
    records: dict[str, tuple[str, float, ScanPoint, tuple[float, ...]]] = {}
    commands = [
        "run set_model_parameters",
        "run set_kinematics_a",
        "set process -p aa_aa kv kinematics.externals.data.helicities=[+1,-1,+1,-1] general.mu_r=91.188 general.m_uv=9.188",
    ]

    for mode in MODES:
        commands.append(mode.sampling)
        for point in POINTS:
            for scale in SCALES:
                xs = rescale_radial(point.xs, scale)
                key = f"{mode.label}_{point.label}_scale_{scale:g}"
                records[key] = (mode.label, scale, point, xs)
                commands.append(
                    "inspect -p aa_aa -i 2L -d "
                    f"{point.graph} -x " + " ".join(f"{x:.17g}" for x in xs)
                )

    body = ["[[command_blocks]]", 'name = "scan_current_maxweights"', "commands = ["]
    for command in commands:
        body.append("  " + json.dumps(command) + ",")
    body.append("]")
    CARD.write_text("\n".join(body) + "\n")
    return records


def run_scan() -> str:
    result = subprocess.run(
        [
            str(GAMMALOOP),
            "--dev-optim",
            "-s",
            str(STATE),
            "--read-only-state",
            str(CARD),
            "run",
            "scan_current_maxweights",
        ],
        cwd=ROOT.parent,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    LOG.write_text(result.stdout)
    if result.returncode != 0:
        raise SystemExit(f"gammaloop inspect scan failed with status {result.returncode}; see {LOG}")
    return result.stdout


def parse_log(
    text: str, records: dict[str, tuple[str, float, ScanPoint, tuple[float, ...]]]
) -> list[dict[str, object]]:
    command_re = re.compile(r"Running command inspect .* -d (\d+) -x ([^\n]+)")
    value_re = re.compile(
        r"The evaluation of integrand '2L' is:\n\n"
        r"\( ([+-][0-9.eE+-]+), ([+-][0-9.eE+-]+) i\)"
    )
    jac_re = re.compile(r"Parameterization jacobian for this point: ([+-][0-9.eE+-]+)")
    chunks = re.split(r"(?=\d{4}-\d{2}-\d{2}T[^\n]*Running command inspect)", text)

    rows = []
    record_items = iter(records.items())
    for chunk in chunks:
        command_match = command_re.search(chunk)
        value_match = value_re.search(chunk)
        jac_match = jac_re.search(chunk)
        if not (command_match and value_match and jac_match):
            continue
        key, (mode, scale, point, xs) = next(record_items)
        r0 = mapped_radius(xs[0])
        r1 = mapped_radius(xs[3])
        rows.append(
            {
                "key": key,
                "mode": mode,
                "label": point.label,
                "graph": point.graph,
                "scale": scale,
                "x0": xs[0],
                "x3": xs[3],
                "radius0": r0,
                "radius1": r1,
                "mean_radius": 0.5 * (r0 + r1),
                "re": float(value_match.group(1)),
                "im": float(value_match.group(2)),
                "abs": math.hypot(float(value_match.group(1)), float(value_match.group(2))),
                "jacobian": float(jac_match.group(1)),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, object]]) -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    (DATA / "current_maxweight_radial_scan.json").write_text(json.dumps(rows, indent=2))
    with (DATA / "current_maxweight_radial_scan.tsv").open("w") as handle:
        keys = [
            "mode",
            "label",
            "graph",
            "scale",
            "x0",
            "x3",
            "mean_radius",
            "re",
            "im",
            "abs",
            "jacobian",
        ]
        handle.write("\t".join(keys) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[key]) for key in keys) + "\n")

    for label in sorted({str(row["label"]) for row in rows}):
        plt.figure(figsize=(6.4, 4.2))
        for mode in [scan_mode.label for scan_mode in MODES]:
            subset = [
                row for row in rows if row["label"] == label and row["mode"] == mode
            ]
            subset.sort(key=lambda row: float(row["mean_radius"]))
            r = [float(row["mean_radius"]) for row in subset]
            abs_val = [float(row["abs"]) for row in subset]
            plt.loglog(r, abs_val, "o-", label=mode)
        plt.xlabel("mean mapped radial coordinate")
        plt.ylabel("inspect value after continuous jacobian")
        plt.title(f"Current max-weight radial scan: {label}")
        plt.grid(True, which="both", alpha=0.25)
        plt.legend()
        plt.tight_layout()
        plt.savefig(FIGURES / f"current_maxweight_{label}_radial_scan.png", dpi=180)
        plt.close()


def main() -> None:
    records = build_card()
    text = run_scan()
    rows = parse_log(text, records)
    write_outputs(rows)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
