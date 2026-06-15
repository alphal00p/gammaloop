#!/usr/bin/env python3
"""Compare ridge scaling channel by channel for representative hard-edge points."""

from __future__ import annotations

import json
import math
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
STATE = Path("/tmp/gammaloop_aa_aa_2l_maxweight_inspect_state")
CARD = Path("/tmp/aa_aa_2l_gl18_per_lmb_channel_scan.toml")
LOG = Path("/tmp/aa_aa_2l_gl18_per_lmb_channel_scan.log")
DATA = ROOT / "data" / "current_maxweight_diagnostics"
GAMMALOOP = ROOT.parent / "gammaloop"

CHANNELS = tuple(range(14))
TARGETS = (
    {
        "label": "GL06_prob1000",
        "graph": 4,
        "xs": (
            0.9998849207926208,
            0.44226073563398893,
            0.2942217815657273,
            0.9999558339703523,
            0.2505645885159602,
            0.06912391384645156,
        ),
    },
    {
        "label": "GL14_bins128_im_minus",
        "graph": 7,
        "xs": (
            0.9992856078626869,
            0.6920929495934958,
            0.28659450247430496,
            0.9998539551872351,
            0.11911699688086141,
            0.05015367802866,
        ),
    },
    {
        "label": "GL18_re_plus_80M",
        "graph": 10,
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
SCALES = (16.0, 8.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.125, 0.0625)


def rescale_radial(xs: tuple[float, ...], scale: float) -> tuple[float, ...]:
    out = list(xs)
    for index in (0, 3):
        out[index] = 1.0 - scale * (1.0 - xs[index])
    return tuple(out)


def mapped_radius(x: float) -> float:
    return x / (1.0 - x)


def slope(rows: list[dict[str, object]]) -> float:
    ordered = sorted(rows, key=lambda row: float(row["mean_radius"]))
    tail = ordered[-4:]
    logs_r = [math.log(float(row["mean_radius"])) for row in tail]
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
        coordinate_system = "spherical"
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
            xs = rescale_radial(tuple(target["xs"]), scale)
            records.append(
                {
                    "label": target["label"],
                    "graph": target["graph"],
                    "mode": "summed",
                    "channel": None,
                    "scale": scale,
                    "xs": xs,
                }
            )
            commands.append(
                "inspect -p aa_aa -i 2L -d "
                f"{target['graph']} -x " + " ".join(f"{x:.17g}" for x in xs)
            )

    commands.append(sampling_block("monte_carlo"))
    for target in TARGETS:
        for channel in CHANNELS:
            for scale in SCALES:
                xs = rescale_radial(tuple(target["xs"]), scale)
                records.append(
                    {
                        "label": target["label"],
                        "graph": target["graph"],
                        "mode": "selected_channel",
                        "channel": channel,
                        "scale": scale,
                        "xs": xs,
                    }
                )
                commands.append(
                    "inspect -p aa_aa -i 2L --discrete-dim "
                    f"{target['graph']} {channel} -x "
                    + " ".join(f"{x:.17g}" for x in xs)
                )

    body = ["[[command_blocks]]", 'name = "scan_gl18_per_lmb_channel"', "commands = ["]
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
        radius0 = mapped_radius(xs[0])
        radius1 = mapped_radius(xs[3])
        re_val = float(value_match.group(1))
        im_val = float(value_match.group(2))
        row = {
            "mode": record["mode"],
            "channel": record["channel"],
            "label": record["label"],
            "graph": record["graph"],
            "scale": record["scale"],
            "x0": xs[0],
            "x3": xs[3],
            "radius0": radius0,
            "radius1": radius1,
            "mean_radius": 0.5 * (radius0 + radius1),
            "re": re_val,
            "im": im_val,
            "abs": math.hypot(re_val, im_val),
            "jacobian": float(jac_match.group(1)),
        }
        rows.append(row)
    return rows


def write_outputs(rows: list[dict[str, object]]) -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    json_path = DATA / "hard_edge_ridges_per_lmb_channel_radial_scan.json"
    tsv_path = DATA / "hard_edge_ridges_per_lmb_channel_radial_scan.tsv"
    summary_path = DATA / "hard_edge_ridges_per_lmb_channel_radial_scan_summary.json"

    json_path.write_text(json.dumps(rows, indent=2))
    keys = (
        "mode",
        "label",
        "channel",
        "graph",
        "scale",
        "mean_radius",
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
                    summed, key=lambda row: float(row["mean_radius"])
                )[-1]["abs"],
            }
        )
        for channel in CHANNELS:
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
                        subset, key=lambda row: float(row["mean_radius"])
                    )[-1]["abs"],
                }
            )
    summary_path.write_text(json.dumps(summary, indent=2))


def main() -> None:
    records = build_card()
    result = subprocess.run(
        [
            str(GAMMALOOP),
            "--dev-optim",
            "-s",
            str(STATE),
            "--read-only-state",
            str(CARD),
            "run",
            "scan_gl18_per_lmb_channel",
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

    rows = parse_log(result.stdout, records)
    write_outputs(rows)
    print(
        (DATA / "hard_edge_ridges_per_lmb_channel_radial_scan_summary.json").read_text()
    )


if __name__ == "__main__":
    main()
