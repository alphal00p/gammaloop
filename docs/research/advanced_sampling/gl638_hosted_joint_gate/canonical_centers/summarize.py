#!/usr/bin/env python3
"""Summarize archived reports and native traces; never load a physics state."""

import gzip
import json
import re
from decimal import Decimal, localcontext
from pathlib import Path

base = Path(__file__).resolve().parent
data = json.loads(
    gzip.decompress((base / "failure-diagnostic-progress.json.gz").read_bytes())
)
old = json.loads(
    gzip.decompress(
        (
            base / "../performance_1f/diagnostic/failure-diagnostic-progress.json.gz"
        ).read_bytes()
    )
)
lines = gzip.decompress((base / "run.log.gz").read_bytes()).decode().splitlines()
audit = {"branches": [], "worker_controls": [], "native_trace": {}}

with localcontext() as context:
    context.prec = 350
    for branch in data["physical_matrix"]:
        rows = [row["result"]["evaluation"] for row in branch["rows"]]
        ordinary = next(
            row["result"]["evaluation"]
            for row in branch["rows"]
            if not row["forced_arb"]
        )
        arb = next(
            row["result"]["evaluation"] for row in branch["rows"] if row["forced_arb"]
        )
        previous = next(
            row for row in old["physical_matrix"] if row["name"] == branch["name"]
        )
        previous_arb = next(
            row["result"]["evaluation"] for row in previous["rows"] if row["forced_arb"]
        )
        scale = sum(
            Decimal(arb["integrand_result"][phase]) ** 2 for phase in ("re", "im")
        ).sqrt()
        differences = {
            phase: Decimal(ordinary["integrand_result"][phase])
            - Decimal(arb["integrand_result"][phase])
            for phase in ("re", "im")
        }
        prior_differences = {
            phase: Decimal(arb["integrand_result"][phase])
            - Decimal(previous_arb["integrand_result"][phase])
            for phase in ("re", "im")
        }
        audit["branches"].append(
            {
                "name": branch["name"],
                "same_baseline_anchor": branch["same_baseline_anchor"],
                "same_anchor_after_evaluations": branch[
                    "same_anchor_after_evaluations"
                ],
                "rows": [
                    {
                        "precision": row["precision"],
                        "valid": row["valid"],
                        "is_nan": row["evaluation_metadata"]["is_nan"],
                        "cut_ids": [
                            event["cut_info"]["cut_id"] for event in row["events"]
                        ],
                        "stability": row["evaluation_metadata"]["stability_results"],
                        "canonical_physical_preparation_time": row[
                            "evaluation_metadata"
                        ]["canonical_physical_preparation_time"],
                    }
                    for row in rows
                ],
                "ordinary_vs_arb": {
                    "norm_relative_error": str(
                        sum(value * value for value in differences.values()).sqrt()
                        / scale
                    ),
                    "components": {
                        phase: {
                            "absolute_error": str(abs(value)),
                            "norm_scaled_error": str(abs(value) / scale),
                            "relative_error": str(
                                abs(value)
                                / abs(Decimal(arb["integrand_result"][phase]))
                            ),
                        }
                        for phase, value in differences.items()
                    },
                    "scope": "Complex-norm stability does not certify relative accuracy of a negligible component; CT-off real part is consistent with zero on the complex scale.",
                },
                "previous_1f_forced_arb_identity": {
                    "total_exact": arb["integrand_result"]
                    == previous_arb["integrand_result"],
                    "six_cut_weights_exact": [
                        event["weight"] for event in arb["events"]
                    ]
                    == [event["weight"] for event in previous_arb["events"]],
                    "norm_relative_change": str(
                        sum(
                            value * value for value in prior_differences.values()
                        ).sqrt()
                        / scale
                    ),
                    "public_anchor_exact": branch["sampling"]["anchor"]
                    == previous["sampling"]["anchor"],
                },
                "reference": {
                    "valid": branch["reference"]["valid"],
                    "reporting_precision": branch["reference"]["reporting_precision"],
                    "metadata": branch["reference"]["evaluation"][
                        "evaluation_metadata"
                    ],
                },
            }
        )
    euler, exact_z = data["physical_matrix"][:2]
    audit["Euler_vs_Pi2Z_returned_identity"] = [
        {
            "precision": left["result"]["evaluation"]["precision"],
            "total_exact": left["result"]["evaluation"]["integrand_result"]
            == right["result"]["evaluation"]["integrand_result"],
            "six_cut_weights_exact": [
                event["weight"] for event in left["result"]["evaluation"]["events"]
            ]
            == [event["weight"] for event in right["result"]["evaluation"]["events"]],
        }
        for left, right in zip(euler["rows"], exact_z["rows"])
    ]
    for worker in data["component_workers"]:
        audit["worker_controls"].append(
            {
                "setup": worker["setup"],
                "sources": [
                    {
                        "name": source["name"],
                        "same_baseline_anchor": source["same_baseline_anchor"],
                        "valid": source["physical"]["evaluation"]["valid"],
                        "precision": source["physical"]["evaluation"]["precision"],
                    }
                    for source in worker["sources"]
                ],
            }
        )
    for branch in ("ct_on_euler", "ct_on_exact_z"):
        start = next(
            i
            for i, line in enumerate(lines)
            if f"failure diagnostic: {branch}, forced_arb=true" in line
        )
        end = next(
            i for i in range(start + 1, len(lines)) if "failure diagnostic:" in lines[i]
        )
        rotations = [
            i for i in range(start, end) if "Evaluating rotation:" in lines[i]
        ][:2]
        records = []
        for begin, stop in zip(rotations, rotations[1:] + [end]):
            rays = {}
            for i in range(begin, stop):
                match = re.search(
                    r"LU (left|right) evaluator input: cut_group_id=0, (?:left|right)_threshold_id=(\d+).*?, r=([^,]+), rstar=(\S+)",
                    lines[i],
                )
                if match:
                    side, index, radius, star = match.groups()
                    rays[f"{side}{index}"] = {"r": radius, "rstar": star, "line": i + 1}
            records.append({"rotation_line": begin + 1, "rays": rays})
        left, right = records
        audit["native_trace"][branch] = {
            "records": records,
            "second_minus_identity": {
                key: {
                    field: str(
                        Decimal(right["rays"][key][field]) - Decimal(value[field])
                    )
                    for field in ("r", "rstar")
                }
                for key, value in left["rays"].items()
            },
        }

audit["center_scope"] = (
    "center_evidence.json is reproduced with ../performance_1f/extract_center_evidence.py. Its center tables now contain unrotated canonical binary64 values. Native r/rstar/alpha test consumption covariance; the display log omits file.active_center vectors, so no direct native-vector equality is inferred from those tables."
)
audit["timing_scope"] = (
    "Traced one-worker diagnostics with RAYON_NUM_THREADS=20; canonical physical preparation is a subset of P. No runtime-budget or scaling claim; detached hosted leaf errors are not completed timings and successful detached costs are nonadditive."
)
(base / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
