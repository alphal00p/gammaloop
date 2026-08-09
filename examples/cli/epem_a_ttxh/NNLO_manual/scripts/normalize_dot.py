#!/usr/bin/env python3
"""Remove drawing and generated momentum-representation attributes from DOT files."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

REMOVED_KEYS = ("pin", "dir", "lmb_rep")
ATTRIBUTE = re.compile(
    r"(?<![A-Za-z0-9_])(?P<key>pin|dir|lmb_rep)\s*=\s*"
    r'(?P<value>"(?:\\.|[^"\\])*"|<[^>]*>|[^,\]\s]+)'
    r"(?P<trailing>\s*,?\s*)"
)
RESIDUAL = re.compile(r"(?<![A-Za-z0-9_])(?:pin|dir|lmb_rep)\s*=")


def unescaped_quote_count(line: str) -> int:
    count = 0
    preceding_backslashes = 0
    for character in line:
        if character == "\\":
            preceding_backslashes += 1
            continue
        if character == '"' and preceding_backslashes % 2 == 0:
            count += 1
        preceding_backslashes = 0
    return count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--output-dir", type=Path)
    output.add_argument("--output", type=Path)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    if args.output is not None and len(args.inputs) != 1:
        parser.error("--output requires exactly one input")
    return args


def normalize_text(raw: str) -> tuple[str, dict[str, int]]:
    counts = {key: 0 for key in REMOVED_KEYS}
    output: list[str] = []
    in_embedded_string = False
    for line in raw.splitlines(keepends=True):
        if "threshold_counterterms" in line and '"' in line:
            output.append(line)
            in_embedded_string = unescaped_quote_count(line) % 2 == 1
            continue
        if in_embedded_string:
            output.append(line)
            if unescaped_quote_count(line) % 2 == 1:
                in_embedded_string = False
            continue
        if "[" not in line or "]" not in line:
            output.append(line)
            continue

        def remove(match: re.Match[str]) -> str:
            counts[match.group("key")] += 1
            return " "

        normalized = ATTRIBUTE.sub(remove, line)
        normalized = re.sub(r"\[\s+,", "[", normalized)
        normalized = re.sub(r",\s*\]", "]", normalized)
        normalized = re.sub(r" {2,}", " ", normalized)
        normalized = normalized.replace("\t ", "\t")
        if RESIDUAL.search(normalized):
            raise ValueError("a removable DOT attribute remains after normalization")
        output.append(normalized)
    normalized_raw = "".join(output)
    return normalized_raw, counts


def target_for(args: argparse.Namespace, source: Path) -> Path:
    if args.output is not None:
        return args.output
    assert args.output_dir is not None
    return args.output_dir / source.name


def main() -> int:
    args = parse_args()
    summaries = []
    for source in args.inputs:
        raw = source.read_text()
        normalized, counts = normalize_text(raw)
        target = target_for(args, source)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(normalized)
        summaries.append(
            {
                "source": str(source),
                "output": str(target),
                "removed": counts,
                "lmb_id_count": len(re.findall(r"\blmb_id\s*=", normalized)),
                "threshold_counterterms": "threshold_counterterms" in normalized,
            }
        )
    if args.summary is not None:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(json.dumps(summaries, indent=2) + "\n")
    else:
        print(json.dumps(summaries, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
