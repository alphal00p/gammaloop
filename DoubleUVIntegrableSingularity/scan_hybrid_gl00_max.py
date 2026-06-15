#!/usr/bin/env python3
"""Scan the GL00 max-weight point found by the hybrid radial sampler."""

from __future__ import annotations

import json
import math
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
STATE = Path("/tmp/gammaloop_aa_aa_2l_all_graphs_invj_sumlmb_state")
CARD = Path("/tmp/aa_aa_2l_hybrid_gl00_max_scan.toml")
LOG = Path("/tmp/aa_aa_2l_hybrid_gl00_max_scan.log")
DATA = ROOT / "data" / "current_maxweight_diagnostics"
GAMMALOOP = ROOT.parent / "gammaloop"
E_CM = 300.0
TARGET = {
    "label": "GL00_hybrid_10M_im_plus",
    "graph": 0,
    "channels": 11,
    "xs": (
        0.9999999760943541,
        0.6826389519199705,
        0.2131752810821759,
        0.820828867381572,
        0.22209636435284774,
        0.9059183610506123,
    ),
}
SCALES = (16.0, 8.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.125, 0.0625)


def branch_radius_x(xs: tuple[float, ...]) -> float:
    return 2.0 * xs[0] - 1.0


def common_radius(xs: tuple[float, ...]) -> float:
    x = branch_radius_x(xs)
    return E_CM * x / (1.0 - x)


def rescale_common_branch(xs: tuple[float, ...], scale: float) -> tuple[float, ...]:
    branch_x = branch_radius_x(xs)
    branch_x = 1.0 - scale * (1.0 - branch_x)
    out = list(xs)
    out[0] = 0.5 + 0.5 * branch_x
    return tuple(out)


def slope(rows: list[dict[str, object]]) -> float:
    positive_rows = [row for row in rows if float(row["abs"]) > 0.0]
    tail = sorted(positive_rows, key=lambda row: float(row["common_radius"]))[-4:]
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
        coordinate_system = "spherical_product_common_radial"
        mapping = "linear"
        b = 1.0
    ' """


def build_card() -> list[dict[str, object]]:
    records = []
    commands = [
        "run set_model_parameters",
        "run set_kinematics_a",
        "set process -p aa_aa kv kinematics.externals.data.helicities=[+1,-1,+1,-1] general.mu_r=91.188 general.m_uv=9.188",
        sampling_block("summed"),
    ]
    for scale in SCALES:
        xs = rescale_common_branch(tuple(TARGET["xs"]), scale)
        records.append({"mode": "summed", "channel": None, "scale": scale, "xs": xs})
        commands.append(
            "inspect -p aa_aa -i 2L -d "
            f"{TARGET['graph']} -x " + " ".join(f"{x:.17g}" for x in xs)
        )

    commands.append(sampling_block("monte_carlo"))
    for channel in range(int(TARGET["channels"])):
        for scale in SCALES:
            xs = rescale_common_branch(tuple(TARGET["xs"]), scale)
            records.append(
                {"mode": "selected_channel", "channel": channel, "scale": scale, "xs": xs}
            )
            commands.append(
                "inspect -p aa_aa -i 2L --discrete-dim "
                f"{TARGET['graph']} {channel} -x "
                + " ".join(f"{x:.17g}" for x in xs)
            )

    body = ["[[command_blocks]]", 'name = "scan_hybrid_gl00_max"', "commands = ["]
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
                "label": TARGET["label"],
                "graph": TARGET["graph"],
                "mode": record["mode"],
                "channel": record["channel"],
                "scale": record["scale"],
                "common_radius": common_radius(xs),
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
    json_path = DATA / "hybrid_gl00_max_radial_scan.json"
    summary_path = DATA / "hybrid_gl00_max_radial_scan_summary.json"
    json_path.write_text(json.dumps(rows, indent=2))

    summary = []
    summed = [row for row in rows if row["mode"] == "summed"]
    summary.append(
        {
            "label": TARGET["label"],
            "mode": "summed",
            "channel": None,
            "slope_last4": slope(summed),
            "abs_at_recorded_point": next(row["abs"] for row in summed if row["scale"] == 1.0),
            "abs_at_deepest_point": sorted(summed, key=lambda row: row["common_radius"])[-1][
                "abs"
            ],
        }
    )
    for channel in range(int(TARGET["channels"])):
        subset = [
            row
            for row in rows
            if row["mode"] == "selected_channel" and row["channel"] == channel
        ]
        summary.append(
            {
                "label": TARGET["label"],
                "mode": "selected_channel",
                "channel": channel,
                "slope_last4": slope(subset),
                "abs_at_recorded_point": next(
                    row["abs"] for row in subset if row["scale"] == 1.0
                ),
                "abs_at_deepest_point": sorted(subset, key=lambda row: row["common_radius"])[-1][
                    "abs"
                ],
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
                "scan_hybrid_gl00_max",
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
