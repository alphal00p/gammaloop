#!/usr/bin/env python3
"""Export JUnit timing metadata without copying captured test output."""

import argparse
import csv
import xml.etree.ElementTree as ET
from decimal import Decimal, InvalidOperation
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("reports", type=Path, nargs="+")
    args = parser.parse_args()

    rows = []
    try:
        args.output.unlink(missing_ok=True)
        for report in args.reports:
            for case in ET.parse(report).iter("testcase"):
                name = case.get("name")
                if not name:
                    raise ValueError("Missing testcase name")
                seconds = case.get("time", "")
                if seconds:
                    duration = Decimal(seconds)
                    if not duration.is_finite() or duration < 0:
                        raise ValueError("Invalid testcase duration")
                    seconds = str(duration)
                status = "passed"
                if case.find("error") is not None:
                    status = "error"
                elif case.find("failure") is not None:
                    status = "failed"
                elif case.find("skipped") is not None:
                    status = "skipped"
                rows.append(
                    [args.check, case.get("classname", ""), name, status, seconds]
                )
    except (ET.ParseError, OSError, InvalidOperation, ValueError):
        parser.error("Could not read valid JUnit testcase timing metadata")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(["check", "test_binary", "test_name", "status", "seconds"])
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
