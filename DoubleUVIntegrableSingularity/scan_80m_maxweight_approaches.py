#!/usr/bin/env python3
"""Inspect radial approaches around the 80M all-graph max-weight points."""

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
CARD = Path("/tmp/aa_aa_2l_80m_maxweight_scan.toml")
LOG = Path("/tmp/aa_aa_2l_80m_maxweight_scan.log")
DATA = ROOT / "data" / "aa_aa_2l_all_graphs_inverse_jacobian_sumlmb_80M"
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
        "GL18_re_plus_80M",
        10,
        (
            0.9999738238928433,
            0.6761869937119277,
            0.479166862842833,
            0.9999638043698932,
            0.017088467296021152,
            0.8935175085829773,
        ),
    ),
    ScanPoint(
        "GL01_re_minus_80M",
        1,
        (
            0.9997436851028048,
            0.42277449086096408,
            0.3044967933923729,
            0.9998203392091457,
            0.47458938787899707,
            0.47185501006813324,
        ),
    ),
    ScanPoint(
        "GL03_im_plus_80M",
        2,
        (
            0.99985130679222345,
            0.7524997097227335,
            0.36093249964439938,
            0.9999698247508588,
            0.6721874876322349,
            0.6252914825413105,
        ),
    ),
    ScanPoint(
        "GL17_im_minus_80M",
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

    body = ["[[command_blocks]]", 'name = "scan_80m_maxweights"', "commands = ["]
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
            "scan_80m_maxweights",
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


def slope(rows: list[dict[str, object]]) -> float:
    ordered = sorted(rows, key=lambda row: float(row["mean_radius"]))
    x = [math.log(float(row["mean_radius"])) for row in ordered[-4:]]
    y = [math.log(float(row["abs"])) for row in ordered[-4:]]
    x_mean = sum(x) / len(x)
    y_mean = sum(y) / len(y)
    numerator = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y))
    denominator = sum((xi - x_mean) ** 2 for xi in x)
    return numerator / denominator


def write_outputs(rows: list[dict[str, object]]) -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    (DATA / "maxweight_80m_radial_scan.json").write_text(json.dumps(rows, indent=2))
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
    with (DATA / "maxweight_80m_radial_scan.tsv").open("w") as handle:
        handle.write("\t".join(keys) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[key]) for key in keys) + "\n")

    summary = []
    for label in sorted({str(row["label"]) for row in rows}):
        plt.figure(figsize=(6.4, 4.2))
        for mode in [scan_mode.label for scan_mode in MODES]:
            subset = [
                row for row in rows if row["label"] == label and row["mode"] == mode
            ]
            subset.sort(key=lambda row: float(row["mean_radius"]))
            summary.append(
                {
                    "label": label,
                    "mode": mode,
                    "slope_last4": slope(subset),
                    "abs_at_recorded_point": next(
                        float(row["abs"]) for row in subset if float(row["scale"]) == 1.0
                    ),
                    "abs_at_deepest_point": float(subset[-1]["abs"]),
                }
            )
            r = [float(row["mean_radius"]) for row in subset]
            abs_val = [float(row["abs"]) for row in subset]
            plt.loglog(r, abs_val, "o-", label=mode)
        plt.xlabel("mean mapped radial coordinate")
        plt.ylabel("inspect value after continuous jacobian")
        plt.title(f"80M max-weight radial scan: {label}")
        plt.grid(True, which="both", alpha=0.25)
        plt.legend()
        plt.tight_layout()
        plt.savefig(FIGURES / f"maxweight_80m_{label}_radial_scan.png", dpi=180)
        plt.close()
    (DATA / "maxweight_80m_radial_scan_summary.json").write_text(
        json.dumps(summary, indent=2)
    )


def main() -> None:
    records = build_card()
    text = run_scan()
    rows = parse_log(text, records)
    write_outputs(rows)
    summary = json.loads((DATA / "maxweight_80m_radial_scan_summary.json").read_text())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
