#!/usr/bin/env python3
"""Compare independent and relative-spherical densities on the same 80M rays."""

from __future__ import annotations

import json
import math
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent
STATE = Path("/tmp/gammaloop_aa_aa_2l_maxweight_inspect_state")
CARD = Path("/tmp/aa_aa_2l_80m_relative_same_physical_scan.toml")
LOG = Path("/tmp/aa_aa_2l_80m_relative_same_physical_scan.log")
DATA = ROOT / "data" / "aa_aa_2l_all_graphs_inverse_jacobian_sumlmb_80M"
GAMMALOOP = ROOT.parent / "gammaloop"
E_CM = 300.0


@dataclass(frozen=True)
class ScanPoint:
    label: str
    graph: int
    xs: tuple[float, float, float, float, float, float]


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

SAMPLING_INDEPENDENT = '''set process -p aa_aa string '
        [sampling]
        graphs = "monte_carlo"
        orientations = "summed"
        lmb_multichanneling = true
        lmb_channels = "summed"
        lmb_channel_weight = "inverse_jacobian"
        coordinate_system = "spherical"
        mapping = "linear"
        b = 1.0
    ' '''

SAMPLING_RELATIVE = '''set process -p aa_aa string '
        [sampling]
        graphs = "monte_carlo"
        orientations = "summed"
        lmb_multichanneling = true
        lmb_channels = "summed"
        lmb_channel_weight = "inverse_jacobian"
        coordinate_system = "relative_spherical"
        mapping = "linear"
        b = 1.0
    ' '''


def rescale_radial(xs: tuple[float, ...], scale: float) -> tuple[float, ...]:
    out = list(xs)
    for i in (0, 3):
        out[i] = 1.0 - scale * (1.0 - xs[i])
    return tuple(out)


def spherical_to_vec(xs: tuple[float, float, float]) -> tuple[float, float, float]:
    radius = E_CM * xs[0] / (1.0 - xs[0])
    phi = 2.0 * math.pi * xs[1]
    cos_theta = -1.0 + 2.0 * xs[2]
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    return (
        radius * sin_theta * math.cos(phi),
        radius * sin_theta * math.sin(phi),
        radius * cos_theta,
    )


def vec_to_spherical(vec: tuple[float, float, float]) -> tuple[float, float, float]:
    x, y, z = vec
    radius = math.sqrt(x * x + y * y + z * z)
    if radius == 0.0:
        return (0.0, 0.0, 0.5)
    x0 = radius / (E_CM + radius)
    phi = math.atan2(y, x)
    if phi < 0.0:
        phi += 2.0 * math.pi
    x1 = phi / (2.0 * math.pi)
    x2 = 0.5 * (1.0 + z / radius)
    return (x0, x1, x2)


def to_relative_spherical_xs(xs: tuple[float, ...]) -> tuple[float, ...]:
    first = spherical_to_vec((xs[0], xs[1], xs[2]))
    second = spherical_to_vec((xs[3], xs[4], xs[5]))
    common = tuple(0.5 * (a + b) for a, b in zip(first, second))
    relative = tuple(a - b for a, b in zip(first, second))
    return vec_to_spherical(common) + vec_to_spherical(relative)


def build_card() -> dict[str, tuple[str, float, ScanPoint, tuple[float, ...]]]:
    records: dict[str, tuple[str, float, ScanPoint, tuple[float, ...]]] = {}
    commands = [
        "run set_model_parameters",
        "run set_kinematics_a",
        "set process -p aa_aa kv kinematics.externals.data.helicities=[+1,-1,+1,-1] general.mu_r=91.188 general.m_uv=9.188",
        SAMPLING_INDEPENDENT,
    ]

    for point in POINTS:
        for scale in SCALES:
            xs = rescale_radial(point.xs, scale)
            key = f"independent_{point.label}_scale_{scale:g}"
            records[key] = ("independent", scale, point, xs)
            commands.append(
                "inspect -p aa_aa -i 2L -d "
                f"{point.graph} -x " + " ".join(f"{x:.17g}" for x in xs)
            )

    commands.append(SAMPLING_RELATIVE)
    for point in POINTS:
        for scale in SCALES:
            independent_xs = rescale_radial(point.xs, scale)
            xs = to_relative_spherical_xs(independent_xs)
            key = f"relative_same_physical_{point.label}_scale_{scale:g}"
            records[key] = ("relative_same_physical", scale, point, xs)
            commands.append(
                "inspect -p aa_aa -i 2L -d "
                f"{point.graph} -x " + " ".join(f"{x:.17g}" for x in xs)
            )

    body = [
        "[[command_blocks]]",
        'name = "scan_80m_relative_same_physical"',
        "commands = [",
    ]
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
            "scan_80m_relative_same_physical",
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
        rows.append(
            {
                "key": key,
                "mode": mode,
                "label": point.label,
                "graph": point.graph,
                "scale": scale,
                "xs": xs,
                "re": float(value_match.group(1)),
                "im": float(value_match.group(2)),
                "abs": math.hypot(float(value_match.group(1)), float(value_match.group(2))),
                "jacobian": float(jac_match.group(1)),
            }
        )
    return rows


def slope(rows: list[dict[str, object]]) -> float:
    ordered = sorted(rows, key=lambda row: float(row["scale"]), reverse=True)
    tail = ordered[-4:]
    radii = [1.0 / float(row["scale"]) for row in tail]
    values = [float(row["abs"]) for row in tail]
    x = [math.log(radius) for radius in radii]
    y = [math.log(value) for value in values]
    x_mean = sum(x) / len(x)
    y_mean = sum(y) / len(y)
    numerator = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y))
    denominator = sum((xi - x_mean) ** 2 for xi in x)
    return numerator / denominator


def write_outputs(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    DATA.mkdir(parents=True, exist_ok=True)
    (DATA / "maxweight_80m_relative_same_physical_scan.json").write_text(
        json.dumps(rows, indent=2)
    )
    keys = ["mode", "label", "graph", "scale", "re", "im", "abs", "jacobian"]
    with (DATA / "maxweight_80m_relative_same_physical_scan.tsv").open("w") as handle:
        handle.write("\t".join(keys) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[key]) for key in keys) + "\n")

    summary = []
    for label in sorted({str(row["label"]) for row in rows}):
        for mode in ("independent", "relative_same_physical"):
            subset = [
                row for row in rows if row["label"] == label and row["mode"] == mode
            ]
            summary.append(
                {
                    "label": label,
                    "mode": mode,
                    "slope_last4": slope(subset),
                    "abs_at_recorded_point": next(
                        float(row["abs"]) for row in subset if float(row["scale"]) == 1.0
                    ),
                    "abs_at_deepest_point": next(
                        float(row["abs"]) for row in subset if float(row["scale"]) == 0.0625
                    ),
                }
            )
    (DATA / "maxweight_80m_relative_same_physical_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    return summary


def main() -> None:
    records = build_card()
    rows = parse_log(run_scan(), records)
    print(json.dumps(write_outputs(rows), indent=2))


if __name__ == "__main__":
    main()
