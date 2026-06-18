#!/usr/bin/env python3
"""Inspect the near-endpoint NaN point found in qq_hhh_2L integration."""

from __future__ import annotations

import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
EXAMPLE = ROOT / "examples" / "cli" / "qq_hhh_2L"
STATE = EXAMPLE / "state"
LOG_DIR = EXAMPLE / "logs" / "nan_endpoint_scan"
GAMMALOOP = ROOT / "target" / "dev-optim" / "gammaloop"

BASE_XS = [
    6.5263450898169484e-01,
    5.7196634612326147e-01,
    1.4361695614139025e-03,
    4.8286308473039141e-01,
    6.3503950391523345e-01,
    8.7877566375154159e-01,
]
X2_VALUES = [
    1.0e-1,
    3.0e-2,
    1.0e-2,
    3.0e-3,
    BASE_XS[2],
    5.0e-4,
    1.0e-4,
]
LMB_CHANNELS = [0, 1, 2]


@dataclass
class InspectSummary:
    label: str
    lmb_channel: int
    x2: float
    returncode: int
    eval_result: str
    metadata_nan: str
    stability_status: str
    max_exponent: int | None
    nan_line_count: int
    ct19_max_exponent: int | None
    ct28_max_exponent: int | None
    log: str


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
    match = re.search(pattern, text)
    return match.group(1).strip() if match else ""


def parse_summary(label: str, lmb_channel: int, x2: float, log_path: Path, result: subprocess.CompletedProcess[str]) -> InspectSummary:
    text = result.stdout
    eval_result = first_match(
        text,
        r"The evaluation of integrand 'subtracted' is:\n\n\((.*?)\)\n",
    )
    metadata_nan = first_match(text, r"nan\s+.*?\|\s+(true|false)\s+\|")
    stability_status = first_match(text, r"\|\s+(Stable\([^|]+?\)|Unstable\([^|]+?\))\s+\|")
    ct19_lines = "\n".join(line for line in text.splitlines() if "threshold_counterterm:19:0" in line)
    ct28_lines = "\n".join(line for line in text.splitlines() if "threshold_counterterm:28:0" in line)
    return InspectSummary(
        label=label,
        lmb_channel=lmb_channel,
        x2=x2,
        returncode=result.returncode,
        eval_result=eval_result,
        metadata_nan=metadata_nan,
        stability_status=stability_status,
        max_exponent=max_exp(text),
        nan_line_count=sum(1 for line in text.splitlines() if "NaN" in line),
        ct19_max_exponent=max_exp(ct19_lines),
        ct28_max_exponent=max_exp(ct28_lines),
        log=str(log_path.relative_to(ROOT)),
    )


def write_outputs(summaries: list[InspectSummary]) -> None:
    rows = [summary.__dict__ for summary in summaries]
    (LOG_DIR / "nan_endpoint_scan_summary.json").write_text(json.dumps(rows, indent=2))

    columns = [
        "label",
        "lmb_channel",
        "x2",
        "returncode",
        "metadata_nan",
        "stability_status",
        "max_exponent",
        "nan_line_count",
        "ct19_max_exponent",
        "ct28_max_exponent",
        "eval_result",
        "log",
    ]
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join("" if row[column] is None else str(row[column]) for column in columns))
    (LOG_DIR / "nan_endpoint_scan_summary.tsv").write_text("\n".join(lines) + "\n")


def main() -> int:
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    setup = run_command(
        [
            str(GAMMALOOP),
            "-o",
            "-s",
            str(STATE),
            "run",
            "runtime_common",
            "sampling_lmb_mc",
        ],
        LOG_DIR / "00_setup_lmb_mc.log",
    )
    if setup.returncode != 0:
        print(setup.stdout)
        return setup.returncode

    summaries: list[InspectSummary] = []
    for x2 in X2_VALUES:
        xs = list(BASE_XS)
        xs[2] = x2
        for lmb_channel in LMB_CHANNELS:
            label = f"x2_{x2:.16e}_lmb{lmb_channel}".replace("+", "")
            log_path = LOG_DIR / f"{label}.log"
            args = [
                str(GAMMALOOP),
                "-n",
                "-s",
                str(STATE),
                "inspect",
                "-p",
                "qq_hhh_2L",
                "-i",
                "subtracted",
                "--discrete-dim",
                "1",
                str(lmb_channel),
                "-f",
                "-x",
                *[f"{value:.17e}" for value in xs],
            ]
            result = run_command(args, log_path)
            summaries.append(parse_summary(label, lmb_channel, x2, log_path, result))

    write_outputs(summaries)
    for summary in summaries:
        print(
            f"lmb={summary.lmb_channel} x2={summary.x2:.6e} "
            f"nan={summary.metadata_nan or '?'} status={summary.stability_status or '?'} "
            f"max_e={summary.max_exponent} nan_lines={summary.nan_line_count} "
            f"ct19_e={summary.ct19_max_exponent} ct28_e={summary.ct28_max_exponent}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
