#!/usr/bin/env python3
"""Probe hard-radial max-weight approaches with the common-radial spherical map."""

from __future__ import annotations

import json
import math
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
STATE = Path("/tmp/gammaloop_aa_aa_2l_all_graphs_invj_sumlmb_state")
CARD = Path("/tmp/aa_aa_2l_common_radial_scan.toml")
LOG = Path("/tmp/aa_aa_2l_common_radial_scan.log")
DATA = ROOT / "data" / "current_maxweight_diagnostics"
GAMMALOOP = ROOT.parent / "gammaloop"
E_CM = 300.0

TARGETS = (
    {
        "label": "GL01_re_minus_80M",
        "graph": 1,
        "channels": 11,
        "xs": (
            0.9997436851028048,
            0.4227744908609641,
            0.3044967933923729,
            0.9998203392091457,
            0.47458938787899707,
            0.47185501006813324,
        ),
    },
    {
        "label": "GL03_im_plus_80M",
        "graph": 2,
        "channels": 11,
        "xs": (
            0.9998513067922234,
            0.7524997097227335,
            0.3609324996443994,
            0.9999698247508588,
            0.6721874876322349,
            0.6252914825413105,
        ),
    },
    {
        "label": "GL17_im_minus_80M",
        "graph": 9,
        "channels": 14,
        "xs": (
            0.9997425961322831,
            0.1305730414150327,
            0.8833503679259446,
            0.9999938714344689,
            0.44900669230836104,
            0.9348993796532749,
        ),
    },
    {
        "label": "GL18_re_plus_80M",
        "graph": 10,
        "channels": 14,
        "xs": (
            0.9999738238928433,
            0.6761869937119277,
            0.479166862842833,
            0.9999638043698932,
            0.017088467296021152,
            0.8935175085829773,
        ),
    },
)
SCALES = (8.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.125, 0.0625)


def radius_from_x(x: float) -> float:
    return E_CM * x / (1.0 - x)


def x_from_radius(radius: float) -> float:
    return radius / (E_CM + radius)


def rescale_independent_radial(xs: tuple[float, ...], scale: float) -> tuple[float, ...]:
    out = list(xs)
    for index in (0, 3):
        out[index] = 1.0 - scale * (1.0 - xs[index])
    return tuple(out)


def independent_to_common_radial(xs: tuple[float, ...]) -> tuple[float, ...]:
    r0 = radius_from_x(xs[0])
    r1 = radius_from_x(xs[3])
    common = r0 + r1
    return (
        x_from_radius(common),
        r0 / common,
        xs[1],
        xs[2],
        xs[4],
        xs[5],
    )


def mean_radius(common_xs: tuple[float, ...]) -> float:
    return radius_from_x(common_xs[0])


def slope(rows: list[dict[str, object]]) -> float:
    ordered = sorted(rows, key=lambda row: float(row["common_radius"]))
    tail = ordered[-4:]
    logs_r = [math.log(float(row["common_radius"])) for row in tail]
    logs_w = [math.log(float(row["abs"])) for row in tail]
    mean_r = sum(logs_r) / len(logs_r)
    mean_w = sum(logs_w) / len(logs_w)
    numerator = sum((r - mean_r) * (w - mean_w) for r, w in zip(logs_r, logs_w))
    denominator = sum((r - mean_r) ** 2 for r in logs_r)
    return numerator / denominator


def sampling_block(lmb_channels: str) -> str:
    return f"""set process -p aa_aa string '
        [sampling]
        graphs = "monte_carlo"
        orientations = "summed"
        lmb_multichanneling = true
        lmb_channels = "{lmb_channels}"
        lmb_channel_weight = "inverse_jacobian"
        coordinate_system = "spherical_common_radial"
        mapping = "linear"
        b = 1.0
    ' """


def build_card() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    commands = [
        "run set_model_parameters",
        "run set_kinematics_a",
        "set process -p aa_aa kv kinematics.externals.data.helicities=[+1,-1,+1,-1] general.mu_r=91.188 general.m_uv=9.188",
        sampling_block("summed"),
    ]

    for target in TARGETS:
        for scale in SCALES:
            common_xs = independent_to_common_radial(
                rescale_independent_radial(tuple(target["xs"]), scale)
            )
            records.append(
                {
                    "label": target["label"],
                    "graph": target["graph"],
                    "mode": "summed",
                    "channel": None,
                    "scale": scale,
                    "xs": common_xs,
                }
            )
            commands.append(
                "inspect -p aa_aa -i 2L -d "
                f"{target['graph']} -x " + " ".join(f"{x:.17g}" for x in common_xs)
            )

    commands.append(sampling_block("monte_carlo"))
    for target in TARGETS:
        for channel in range(int(target["channels"])):
            for scale in SCALES:
                common_xs = independent_to_common_radial(
                    rescale_independent_radial(tuple(target["xs"]), scale)
                )
                records.append(
                    {
                        "label": target["label"],
                        "graph": target["graph"],
                        "mode": "selected_channel",
                        "channel": channel,
                        "scale": scale,
                        "xs": common_xs,
                    }
                )
                commands.append(
                    "inspect -p aa_aa -i 2L --discrete-dim "
                    f"{target['graph']} {channel} -x "
                    + " ".join(f"{x:.17g}" for x in common_xs)
                )

    body = ["[[command_blocks]]", 'name = "scan_common_radial"', "commands = ["]
    for command in commands:
        body.append("  " + json.dumps(command) + ",")
    body.append("]")
    CARD.write_text("\n".join(body) + "\n")
    return records


def parse_log(text: str, records: list[dict[str, object]]) -> list[dict[str, object]]:
    value_re = re.compile(
        r"The evaluation of integrand '2L' is:\n\n"
        r"\( ([+-][0-9.eE+-]+), ([+-][0-9.eE+-]+) i\)"
    )
    jac_re = re.compile(r"Parameterization jacobian for this point: ([+-][0-9.eE+-]+)")
    chunks = re.split(r"(?=\d{4}-\d{2}-\d{2}T[^\n]*Running command inspect)", text)

    rows = []
    record_iter = iter(records)
    for chunk in chunks:
        value_match = value_re.search(chunk)
        jac_match = jac_re.search(chunk)
        if not (value_match and jac_match):
            continue
        record = next(record_iter)
        xs = tuple(float(x) for x in record["xs"])
        re_val = float(value_match.group(1))
        im_val = float(value_match.group(2))
        rows.append(
            {
                "mode": record["mode"],
                "channel": record["channel"],
                "label": record["label"],
                "graph": record["graph"],
                "scale": record["scale"],
                "common_radius": mean_radius(xs),
                "fraction": xs[1],
                "re": re_val,
                "im": im_val,
                "abs": math.hypot(re_val, im_val),
                "jacobian": float(jac_match.group(1)),
            }
        )
    return rows


def write_outputs(rows: list[dict[str, object]]) -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    json_path = DATA / "common_radial_ridge_scan.json"
    tsv_path = DATA / "common_radial_ridge_scan.tsv"
    summary_path = DATA / "common_radial_ridge_scan_summary.json"

    json_path.write_text(json.dumps(rows, indent=2))
    keys = (
        "mode",
        "label",
        "channel",
        "graph",
        "scale",
        "common_radius",
        "fraction",
        "re",
        "im",
        "abs",
        "jacobian",
    )
    with tsv_path.open("w") as handle:
        handle.write("\t".join(keys) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[key]) for key in keys) + "\n")

    summary = []
    for target in TARGETS:
        label = target["label"]
        summed = [
            row for row in rows if row["mode"] == "summed" and row["label"] == label
        ]
        summary.append(
            {
                "label": label,
                "mode": "summed",
                "channel": None,
                "slope_last4": slope(summed),
                "abs_at_recorded_point": next(
                    row["abs"] for row in summed if float(row["scale"]) == 1.0
                ),
                "abs_at_deepest_point": sorted(
                    summed, key=lambda row: float(row["common_radius"])
                )[-1]["abs"],
            }
        )
        channel_count = next(int(target["channels"]) for target in TARGETS if target["label"] == label)
        for channel in range(channel_count):
            subset = [
                row
                for row in rows
                if row["mode"] == "selected_channel"
                and row["label"] == label
                and row["channel"] == channel
            ]
            summary.append(
                {
                    "label": label,
                    "mode": "selected_channel",
                    "channel": channel,
                    "slope_last4": slope(subset),
                    "abs_at_recorded_point": next(
                        row["abs"] for row in subset if float(row["scale"]) == 1.0
                    ),
                    "abs_at_deepest_point": sorted(
                        subset, key=lambda row: float(row["common_radius"])
                    )[-1]["abs"],
                }
            )
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"wrote {json_path}")
    print(f"wrote {summary_path}")


def main() -> None:
    records = build_card()
    with LOG.open("w") as handle:
        subprocess.run(
            [
                str(GAMMALOOP),
                "-s",
                str(STATE),
                "--read-only-state",
                str(CARD),
                "run",
                "scan_common_radial",
            ],
            check=True,
            stdout=handle,
            stderr=subprocess.STDOUT,
            cwd=ROOT.parent,
        )
    rows = parse_log(LOG.read_text(errors="replace"), records)
    if len(rows) != len(records):
        raise RuntimeError(f"Parsed {len(rows)} rows for {len(records)} inspect records")
    write_outputs(rows)


if __name__ == "__main__":
    main()
