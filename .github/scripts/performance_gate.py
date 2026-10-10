#!/usr/bin/env python3
"""Compare paired performance measurements; emit only numeric report values.

Schema 1 has baseline/candidate commit and fixture/corpus/config SHA256 metadata,
required_metrics, correctness_tolerance {absolute, relative}, and discovery.
Each discovery/confirmation has a positive run_id and five alternating AB/BA
pairs. A pair has baseline/candidate metric dictionaries and baseline_output /
candidate_output numeric vectors. Confirmation starts in the opposite order and
has a distinct run_id. Generation/load are seconds; warm runtime is seconds per
sample. CLI commit arguments bind measurements to the intended revisions.

Exit codes: 0 clean, 1 confirmed regression or correctness mismatch, 2 invalid,
incomplete, noisy, or awaiting independent confirmation. Timing thresholds are
owned here; benchmark configuration owns workload and correctness tolerances.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import statistics
import sys
from pathlib import Path
from typing import ClassVar


class InvalidMeasurements(ValueError):
    pass


class PerformanceGate:
    THRESHOLDS: ClassVar[dict[str, tuple[float, float]]] = {
        "generation_seconds": (0.20, 0.5),
        "load_seconds": (0.20, 0.1),
        "warm_runtime_seconds_per_sample": (0.15, 1e-6),
    }
    # Normal-consistent MAD is a robust estimate of paired dispersion. The
    # three-noise margin is a screening heuristic, not a confidence interval.
    MAD_SCALE = 1.4826
    NOISE_MULTIPLIER = 3

    def __init__(self, baseline_commit: str, candidate_commit: str):
        self.baseline_commit = baseline_commit
        self.candidate_commit = candidate_commit

    @staticmethod
    def number(value: object) -> float:
        if type(value) not in (int, float):
            raise InvalidMeasurements("Expected a finite number")
        try:
            number = float(value)
        except OverflowError as error:
            raise InvalidMeasurements("Number is out of range") from error
        if not math.isfinite(number):
            raise InvalidMeasurements("Expected a finite number")
        return number

    @staticmethod
    def fields(value: object, required: set[str], optional: set[str] = frozenset()):
        if not isinstance(value, dict) or not required <= value.keys():
            raise InvalidMeasurements("Missing required fields")
        if value.keys() - required - optional:
            raise InvalidMeasurements("Unexpected fields")
        return value

    @staticmethod
    def fingerprint(value: object, length: int):
        if not isinstance(value, str) or not re.fullmatch(
            f"[0-9a-f]{{{length}}}", value
        ):
            raise InvalidMeasurements("Invalid identity fingerprint")

    def validate_identity(self, data: dict, report: dict):
        keys = {"commit", "fixture_sha256", "corpus_sha256", "config_sha256"}
        baseline = self.fields(data["baseline"], keys)
        candidate = self.fields(data["candidate"], keys)
        for metadata in (baseline, candidate):
            self.fingerprint(metadata["commit"], 40)
            for key in keys - {"commit"}:
                self.fingerprint(metadata[key], 64)
        if (
            baseline["commit"] != self.baseline_commit
            or candidate["commit"] != self.candidate_commit
            or any(baseline[key] != candidate[key] for key in keys - {"commit"})
        ):
            report["identity_mismatch"] = 1
            raise InvalidMeasurements("Measurement identities do not match")

    def validate_phase(self, phase: object, metrics: list[str], tolerance: dict):
        phase = self.fields(phase, {"run_id", "pairs"})
        if type(phase["run_id"]) is not int or phase["run_id"] <= 0:
            raise InvalidMeasurements("Invalid run identifier")
        pairs = phase["pairs"]
        if not isinstance(pairs, list) or len(pairs) != 5:
            raise InvalidMeasurements("Exactly five pairs are required")
        results = {metric: [] for metric in metrics}
        mismatches = 0
        previous_order = None
        output_size = None
        for pair in pairs:
            pair = self.fields(
                pair,
                {
                    "order",
                    "baseline",
                    "candidate",
                    "baseline_output",
                    "candidate_output",
                },
            )
            order = pair["order"]
            if order not in ("AB", "BA") or order == previous_order:
                raise InvalidMeasurements("Pair order must alternate AB and BA")
            previous_order = order
            baseline = self.fields(pair["baseline"], set(metrics))
            candidate = self.fields(pair["candidate"], set(metrics))
            for metric in metrics:
                before = self.number(baseline[metric])
                after = self.number(candidate[metric])
                if before <= 0 or after <= 0:
                    raise InvalidMeasurements("Timings must be positive")
                results[metric].append((before, after))
            before_output, after_output = (
                pair["baseline_output"],
                pair["candidate_output"],
            )
            if (
                not isinstance(before_output, list)
                or not isinstance(after_output, list)
                or not before_output
                or len(before_output) != len(after_output)
                or (output_size is not None and len(before_output) != output_size)
            ):
                raise InvalidMeasurements(
                    "Correctness output vectors must match in size"
                )
            output_size = len(before_output)
            before_output = [self.number(value) for value in before_output]
            after_output = [self.number(value) for value in after_output]
            if not any(before_output):
                raise InvalidMeasurements(
                    "A physical correctness reference must be nonzero"
                )
            pair_mismatches = sum(
                not math.isclose(
                    before,
                    after,
                    rel_tol=tolerance["relative"],
                    abs_tol=tolerance["absolute"],
                )
                for before, after in zip(before_output, after_output)
            )
            # Even an absolute tolerance must not accept a completely zeroed
            # physical candidate. The producer also checks each complex point.
            mismatches += max(pair_mismatches, int(not any(after_output)))
        return results, mismatches, output_size

    def summarize(self, samples: list[tuple[float, float]], metric: str) -> dict:
        relative_limit, absolute_limit = self.THRESHOLDS[metric]
        deltas = [self.number(after - before) for before, after in samples]
        relative = [self.number(after / before - 1) for before, after in samples]
        delta = statistics.median(deltas)
        effect = statistics.median(relative)
        noise = self.number(
            self.MAD_SCALE
            * statistics.median(abs(value - effect) for value in relative)
        )
        absolute_noise = self.number(
            self.MAD_SCALE * statistics.median(abs(value - delta) for value in deltas)
        )
        margin = self.number(self.NOISE_MULTIPLIER * noise)
        absolute_margin = self.number(self.NOISE_MULTIPLIER * absolute_noise)
        slower = sum(after > before for before, after in samples)
        exceeds_limits = effect > relative_limit and delta > absolute_limit
        qualified = exceeds_limits and slower >= 4 and effect > margin
        # A noisy phase cannot pass merely because its median happens to fall
        # below a threshold. Also flag a near-threshold uncertainty interval.
        uncertain = (
            margin >= relative_limit
            or (exceeds_limits and not qualified)
            or (
                not qualified
                and effect + margin > relative_limit
                and delta + absolute_margin > absolute_limit
            )
        )
        return {
            "status": 1 if qualified else 2 if uncertain else 0,
            "baseline_median": statistics.median(before for before, _ in samples),
            "candidate_median": statistics.median(after for _, after in samples),
            "median_absolute_delta": delta,
            "median_relative_delta": effect,
            "paired_noise": noise,
            "paired_absolute_noise": absolute_noise,
            "slower_pairs": slower,
            "pair_count": len(samples),
            "relative_threshold": relative_limit,
            "absolute_threshold": absolute_limit,
        }

    def analyze(self, data: object) -> dict:
        report = {
            "schema_version": 1,
            "status": 2,
            "confirmation_required": 0,
            "invalid_input": 0,
            "identity_mismatch": 0,
            "correctness_mismatches": 0,
            "metrics": {},
        }
        try:
            data = self.fields(
                data,
                {
                    "schema_version",
                    "baseline",
                    "candidate",
                    "required_metrics",
                    "correctness_tolerance",
                    "discovery",
                },
                {"confirmation"},
            )
            if type(data["schema_version"]) is not int or data["schema_version"] != 1:
                raise InvalidMeasurements("Unsupported schema version")
            self.validate_identity(data, report)
            metrics = data["required_metrics"]
            if (
                not isinstance(metrics, list)
                or not metrics
                or any(not isinstance(metric, str) for metric in metrics)
                or len(set(metrics)) != len(metrics)
                or not set(metrics) <= self.THRESHOLDS.keys()
            ):
                raise InvalidMeasurements("Invalid required metrics")
            tolerance = self.fields(
                data["correctness_tolerance"], {"absolute", "relative"}
            )
            tolerance = {key: self.number(value) for key, value in tolerance.items()}
            if any(value < 0 for value in tolerance.values()):
                raise InvalidMeasurements("Correctness tolerances cannot be negative")
            discovery, mismatch, output_size = self.validate_phase(
                data["discovery"], metrics, tolerance
            )
            report["correctness_mismatches"] = mismatch
            confirmation = None
            if "confirmation" in data:
                confirmation, mismatch, confirmation_size = self.validate_phase(
                    data["confirmation"], metrics, tolerance
                )
                if (
                    data["confirmation"]["run_id"] == data["discovery"]["run_id"]
                    or data["confirmation"]["pairs"][0]["order"]
                    == data["discovery"]["pairs"][0]["order"]
                    or output_size != confirmation_size
                ):
                    raise InvalidMeasurements(
                        "Confirmation must be an independent matched run"
                    )
                report["correctness_mismatches"] += mismatch
            statuses = []
            for metric in metrics:
                first = self.summarize(discovery[metric], metric)
                entry = {"discovery": first}
                if confirmation is not None:
                    second = self.summarize(confirmation[metric], metric)
                    entry["confirmation"] = second
                    status = (
                        1
                        if first["status"] == second["status"] == 1
                        else 0
                        if first["status"] == second["status"] == 0
                        else 2
                    )
                else:
                    status = 0 if first["status"] == 0 else 2
                    report["confirmation_required"] |= first["status"] == 1
                entry["status"] = status
                report["metrics"][metric] = entry
                statuses.append(status)
            report["status"] = (
                1
                if report["correctness_mismatches"] or 1 in statuses
                else 2
                if 2 in statuses
                else 0
            )
        except (InvalidMeasurements, ArithmeticError):
            report["invalid_input"] = 1
            report["status"] = 2
            report["metrics"] = {}
        return report

    @staticmethod
    def print_report(report: dict) -> None:
        """Describe an analyzed numeric report without exposing measurement inputs."""
        outcome = {0: "pass", 1: "confirmed regression", 2: "inconclusive"}[
            report["status"]
        ]
        if report["correctness_mismatches"]:
            outcome = "correctness failure"
        print(f"Performance gate: {outcome}", file=sys.stderr)
        phase_labels = {0: "within policy", 1: "regression signal", 2: "inconclusive"}
        for metric, entry in report["metrics"].items():
            for phase in ("discovery", "confirmation"):
                if phase not in entry:
                    continue
                summary = entry[phase]
                percentage = summary["median_relative_delta"] * 100
                change = (
                    f"{percentage:+.2f}%"
                    if math.isfinite(percentage)
                    else f"{summary['median_relative_delta']:+.3g} relative change"
                )
                print(
                    f"  {metric} {phase}: {summary['baseline_median']:.6g} -> "
                    f"{summary['candidate_median']:.6g} ({change}); "
                    f"policy >{summary['relative_threshold']:.0%} and "
                    f">{summary['absolute_threshold']:.6g}; "
                    f"paired relative noise {summary['paired_noise']:.3g}; "
                    f"{phase_labels[summary['status']]}",
                    file=sys.stderr,
                )
        if report["confirmation_required"]:
            print("  Independent confirmation required", file=sys.stderr)
        if report["invalid_input"]:
            print("  Invalid or incomplete measurements", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--baseline-commit", required=True)
    parser.add_argument("--candidate-commit", required=True)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args(argv)
    gate = PerformanceGate(args.baseline_commit, args.candidate_commit)
    try:
        with args.input.open() as stream:
            data = json.load(stream)
    except (OSError, ValueError):
        data = None
    report = gate.analyze(data)
    try:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    except OSError:
        print("Could not write numeric performance report", file=sys.stderr)
        return 2
    gate.print_report(report)
    return report["status"]


if __name__ == "__main__":
    sys.exit(main())
