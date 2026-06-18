#!/usr/bin/env python3
"""Run targeted inspect calls for qq_hhh_2L NaN debugging."""

from __future__ import annotations

import json
import re
import subprocess
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "examples" / "cli" / "qq_hhh_2L"
STATE = EXAMPLE / "state"
LOG_DIR = EXAMPLE / "logs" / "nan_point_matrix"
GAMMALOOP = ROOT / "target" / "dev-optim" / "gammaloop"
ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")

POINTS = {
    "lmb_summed_max": [
        6.5263450898169484e-01,
        5.7196634612326147e-01,
        1.4361695614139025e-03,
        4.8286308473039141e-01,
        6.3503950391523345e-01,
        8.7877566375154159e-01,
    ],
    "unit_mid": [5.0e-01] * 6,
    "random_bulk": [3.7e-01, 4.2e-01, 6.1e-01, 2.8e-01, 7.4e-01, 5.3e-01],
    "known_variance": [
        3.1390634220043045e-01,
        6.1945042268059569e-01,
        9.0716850422920681e-01,
        7.1401078676936824e-01,
        4.9683778437982462e-01,
        7.8198648800276493e-02,
    ],
    "smoke_max": [
        5.4906782482100000e-01,
        6.9976306919957410e-01,
        5.3783697990307630e-01,
        2.0101240866401980e-01,
        6.6644143106342780e-01,
        3.1764002852414430e-01,
    ],
}

GROUPS = [0, 1]
LMB_CHANNELS = [0, 1, 2]


@dataclass
class InspectSummary:
    point: str
    group: int
    lmb_channel: int
    returncode: int
    metadata_nan: str
    stability_status: str
    eval_result: str
    max_exponent: int | None
    nan_line_count: int
    ct19_max_exponent: int | None
    ct28_max_exponent: int | None
    original_max_exponent: int | None
    log: str
    fatal_error: str


def run_command(args: list[str], log_path: Path | None = None) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(
        args,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if log_path is not None:
        log_path.write_text(result.stdout)
    return result


def max_exp(text: str) -> int | None:
    exps = [int(match.group(1)) for match in re.finditer(r"e([+-]?\d+)", text)]
    if not exps:
        return None
    return max(exps)


def first_match(text: str, pattern: str) -> str:
    match = re.search(pattern, text, flags=re.DOTALL)
    return match.group(1).strip() if match else ""


def clean_output(text: str) -> str:
    return ANSI_RE.sub("", text)


def table_value(text: str, name: str) -> str:
    match = re.search(rf"│\s*{re.escape(name)}\s*│\s*([^│]+?)\s*│", text)
    return match.group(1).strip() if match else ""


def lines_with(text: str, needle: str) -> str:
    return "\n".join(line for line in text.splitlines() if needle in line)


def parse_summary(
    point: str,
    group: int,
    lmb_channel: int,
    log_path: Path,
    result: subprocess.CompletedProcess[str],
) -> InspectSummary:
    text = clean_output(result.stdout)
    eval_result = first_match(
        text,
        r"The evaluation of integrand 'subtracted' is:\s*\n\s*\((.*?)\)\s*\n",
    )
    stability_status = first_match(text, r"│\s*(Stable\([^│]+?\)|Unstable\([^│]+?\))\s*│")
    fatal_error = first_match(text, r"\n\s*0:\s*(.*?)\n\nLocation:")
    return InspectSummary(
        point=point,
        group=group,
        lmb_channel=lmb_channel,
        returncode=result.returncode,
        metadata_nan=table_value(text, "nan"),
        stability_status=stability_status,
        eval_result=eval_result,
        max_exponent=max_exp(text),
        nan_line_count=sum(1 for line in text.splitlines() if "NaN" in line),
        ct19_max_exponent=max_exp(lines_with(text, "threshold_counterterm:19:0")),
        ct28_max_exponent=max_exp(lines_with(text, "threshold_counterterm:28:0")),
        original_max_exponent=max_exp(lines_with(text, "original")),
        log=str(log_path.relative_to(ROOT)),
        fatal_error=fatal_error,
    )


def write_outputs(log_dir: Path, summaries: list[InspectSummary]) -> None:
    rows = [summary.__dict__ for summary in summaries]
    (log_dir / "summary.json").write_text(json.dumps(rows, indent=2))

    columns = [
        "point",
        "group",
        "lmb_channel",
        "returncode",
        "metadata_nan",
        "stability_status",
        "max_exponent",
        "nan_line_count",
        "ct19_max_exponent",
        "ct28_max_exponent",
        "original_max_exponent",
        "eval_result",
        "log",
        "fatal_error",
    ]
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join("" if row[column] is None else str(row[column]) for column in columns))
    (log_dir / "summary.tsv").write_text("\n".join(lines) + "\n")


def parse_csv_ints(value: str) -> list[int]:
    return [int(item.strip()) for item in value.split(",") if item.strip()]


def parse_args() -> tuple[Path, Path, list[str], list[int], list[int], str]:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--state", type=Path, default=STATE)
    parser.add_argument("--log-dir", type=Path, default=LOG_DIR)
    parser.add_argument("--points", default=",".join(POINTS))
    parser.add_argument("--groups", default=",".join(str(group) for group in GROUPS))
    parser.add_argument(
        "--lmb-channels",
        default=",".join(str(channel) for channel in LMB_CHANNELS),
    )
    parser.add_argument("--runtime-block", default="runtime_common")
    args = parser.parse_args()

    point_names = [item.strip() for item in args.points.split(",") if item.strip()]
    unknown_points = [point for point in point_names if point not in POINTS]
    if unknown_points:
        parser.error(f"unknown points: {', '.join(unknown_points)}")

    state = args.state if args.state.is_absolute() else ROOT / args.state
    log_dir = args.log_dir if args.log_dir.is_absolute() else ROOT / args.log_dir

    return (
        state,
        log_dir,
        point_names,
        parse_csv_ints(args.groups),
        parse_csv_ints(args.lmb_channels),
        args.runtime_block,
    )


def main() -> int:
    state, log_dir, point_names, groups, lmb_channels, runtime_block = parse_args()
    log_dir.mkdir(parents=True, exist_ok=True)

    setup = run_command(
        [
            str(GAMMALOOP),
            "-o",
            "-s",
            str(state),
            "run",
            runtime_block,
            "sampling_lmb_mc",
        ],
        log_dir / "00_setup_lmb_mc.log",
    )
    if setup.returncode != 0:
        print(setup.stdout)
        return setup.returncode

    summaries: list[InspectSummary] = []
    for point in point_names:
        xs = POINTS[point]
        for group in groups:
            for lmb_channel in lmb_channels:
                log_path = log_dir / f"{point}_group{group}_lmb{lmb_channel}.log"
                result = run_command(
                    [
                        str(GAMMALOOP),
                        "-n",
                        "-s",
                        str(state),
                        "inspect",
                        "-p",
                        "qq_hhh_2L",
                        "-i",
                        "subtracted",
                        "--discrete-dim",
                        str(group),
                        str(lmb_channel),
                        "-f",
                        "-x",
                        *[f"{value:.17e}" for value in xs],
                    ],
                    log_path,
                )
                summary = parse_summary(point, group, lmb_channel, log_path, result)
                summaries.append(summary)
                print(
                    f"{point:14s} group={group} lmb={lmb_channel} "
                    f"nan={summary.metadata_nan or '?':5s} "
                    f"status={summary.stability_status or '?':20s} "
                    f"max_e={summary.max_exponent} nan_lines={summary.nan_line_count} "
                    f"ct19_e={summary.ct19_max_exponent} ct28_e={summary.ct28_max_exponent} "
                    f"fatal={summary.fatal_error[:72] if summary.fatal_error else '-'}"
                )

    write_outputs(log_dir, summaries)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
