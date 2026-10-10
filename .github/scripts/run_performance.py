#!/usr/bin/env python3
"""Measure paired fresh generation and warmed fixed-point CLI benchmarks."""

import argparse
import hashlib
import json
import math
import os
import re
import shlex
import signal
import statistics
import subprocess
import sys
import time
import tomllib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from string import Template

from performance_gate import PerformanceGate


def positive_number(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{label} must be a number")
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{label} must be positive and finite")
    return value


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def check_deadline(deadline):
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise ValueError("The performance harness exceeded its 240-second budget")
    return remaining


class Fixture:
    @staticmethod
    def cases(suite):
        cases = {
            f"{name}-{route}" if route else name: (workload, route)
            for name, workload in suite["workloads"].items()
            for route in (suite["routes"] if workload["family"] == "scalar" else [None])
        }
        if not cases or any(
            not re.fullmatch(r"[A-Za-z0-9_-]+", name) for name in cases
        ):
            raise ValueError("Invalid performance case names")
        return cases

    @staticmethod
    def source(input_root, name):
        path = Path(name)
        source = (input_root / path).resolve()
        if (
            path.is_absolute()
            or ".." in path.parts
            or not source.is_relative_to(input_root)
        ):
            raise ValueError("Fixture inputs must stay in the shared repository")
        return source

    @classmethod
    def from_suite(cls, suite, name, input_root):
        workload, route = cls.cases(suite)[name]
        family = suite[workload["family"]]
        config = suite["settings"] | {
            key: value for key, value in family.items() if key not in ("card", "corpus")
        }
        corpus = family["corpus"] | {
            "process": workload["process"],
            "graph_name": workload["graph_name"],
        }
        directory = Path("benchmarks/performance")
        card = cls.source(input_root, str(directory / family["card"])).read_text()
        if route:
            graph = str(directory / workload["graph"])
            corpus["inputs"] = [graph, *corpus["inputs"]]
            values = suite["routes"][route] | {
                "process": corpus["process"],
                "graph": graph,
            }
            card = Template(card).substitute(
                {
                    key: str(value).lower() if type(value) is bool else value
                    for key, value in values.items()
                }
            )
        return cls(card, corpus, config, input_root)

    def __init__(self, card, corpus, config, input_root):
        self.card = card
        self.input_root = input_root.resolve()
        self.corpus = corpus
        self.config = config
        fixture_hash = hashlib.sha256(card.encode())
        inputs = self.corpus["inputs"]
        if not isinstance(inputs, list) or any(not isinstance(p, str) for p in inputs):
            raise ValueError("Fixture inputs must be relative paths")
        if len(set(inputs)) != len(inputs):
            raise ValueError("Fixture inputs must be distinct")
        for input_file in inputs:
            source = self.source(self.input_root, input_file)
            fixture_hash.update(b"\0" + input_file.encode() + b"\0")
            fixture_hash.update(source.read_bytes())
        self.hashes = {
            "fixture_sha256": fixture_hash.hexdigest(),
            "corpus_sha256": hashlib.sha256(
                json.dumps(corpus, sort_keys=True, allow_nan=False).encode()
            ).hexdigest(),
            "config_sha256": hashlib.sha256(
                json.dumps(config, sort_keys=True, allow_nan=False).encode()
            ).hexdigest(),
        }
        card = tomllib.loads(card)
        if card["default_runtime_settings"]["general"]["enable_cache"] is not False:
            raise ValueError("The production benchmark must disable point caching")
        if self.config["schema_version"] != 1 or self.config["pairs"] != 5:
            raise ValueError("The performance gate requires schema 1 and five pairs")
        if self.config["required_metrics"] != [
            "generation_seconds",
            "warm_runtime_seconds_per_sample",
        ]:
            raise ValueError("Unexpected required performance metrics")
        self.timeout = positive_number(
            self.config["fixture_timeout_seconds"], "Fixture timeout"
        )
        if self.timeout > 240:
            raise ValueError("The fixture timeout must not exceed 240 seconds")
        positive_number(self.config["bench_duration_seconds"], "Bench duration")
        batches = self.config["bench_batches"]
        if isinstance(batches, bool) or not isinstance(batches, int) or batches < 2:
            raise ValueError("Bench requires at least two batches")
        for label in ("process", "integrand", "graph_name"):
            if not re.fullmatch(r"[A-Za-z0-9_]+", self.corpus[label]):
                raise ValueError(f"Invalid corpus {label}")
        if type(self.corpus["process_id"]) is not int or self.corpus["process_id"] < 0:
            raise ValueError("Invalid process id")
        dimension = self.corpus["continuous_dimension"]
        if type(dimension) is not int or dimension <= 0 or dimension % 3:
            raise ValueError("The corpus requires three coordinates per loop")
        dims = self.corpus["discrete_dim"]
        if not isinstance(dims, list) or any(type(v) is not int or v < 0 for v in dims):
            raise ValueError("Invalid discrete dimensions")
        points = self.corpus["points"]
        if len(points) < 2 or any(len(p) != dimension for p in points):
            raise ValueError(
                "The corpus requires multiple points of the declared dimension"
            )
        for point in points:
            if any(
                isinstance(v, bool)
                or not isinstance(v, (int, float))
                or not math.isfinite(v)
                or not 0 < v < 1
                for v in point
            ):
                raise ValueError("Invalid x-space point")
        if len({tuple(p) for p in points}) != len(points):
            raise ValueError("The fixed-point corpus must contain distinct points")

    def benchmark_card(self, directory):
        commands = []
        for index, point in enumerate(self.corpus["points"]):
            point_args = " ".join(format(v, ".17g") for v in point)
            dims = self.corpus["discrete_dim"]
            selector = "-d " + " ".join(str(v) for v in dims) if dims else ""
            path = directory / f"bench-{index}.json"
            commands.append(
                f"bench -p {self.corpus['process']} -i {self.corpus['integrand']} "
                f"-x {point_args} {selector} "
                f"--duration {self.config['bench_duration_seconds']} "
                f"--n-batches {self.config['bench_batches']} "
                f"--json-output {shlex.quote(str(path))}"
            )
        commands.append("quit --no-save-state")
        card = directory / "benchmark.toml"
        card.write_text("commands = " + json.dumps(commands, ensure_ascii=False) + "\n")
        return card


def run_cli(command, repository, log_path, deadline):
    check_deadline(deadline)
    environment = os.environ.copy()
    environment.update({"RAYON_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"})
    with log_path.open("w") as log:
        process = subprocess.Popen(
            command,
            cwd=repository,
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        # POSIX wait(timeout=...) polls with sleeps of up to 50 ms. Wait for
        # the actual exit in a worker and enforce the shared deadline separately.
        with ThreadPoolExecutor(max_workers=1) as waiter:
            try:
                completion = waiter.submit(process.wait)
                status = completion.result(timeout=check_deadline(deadline))
            except TimeoutError:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    completion.result(timeout=5)
                except TimeoutError:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    completion.result()
                raise ValueError(
                    "The performance harness exceeded its 240-second budget"
                )
            finally:
                # A successful parent must not leave compiler or evaluator children.
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
    if status:
        raise ValueError(f"CLI exited with status {status}; private log: {log_path}")


def benchmark_result(path, fixture, index):
    value = json.loads(path.read_text())
    if value["point"] != fixture.corpus["points"][index]:
        raise ValueError("Bench did not evaluate the requested corpus point")
    if value["discrete_dim"] != fixture.corpus["discrete_dim"]:
        raise ValueError("Bench changed the requested discrete dimensions")
    if (
        value["process_id"] != fixture.corpus["process_id"]
        or value["integrand_name"] != fixture.corpus["integrand"]
    ):
        raise ValueError("Bench changed the requested integrand")
    if value["target_seconds"] != fixture.config["bench_duration_seconds"]:
        raise ValueError("Bench changed the requested target duration")
    if value["warmup_samples"] != 10:
        raise ValueError("Bench did not complete the expected ten warmup samples")
    positive_number(value["warmup_seconds"], "Bench warmup duration")
    if value["minimal_integrand"] is not False or value["momentum_space"] is not False:
        raise ValueError("Bench must evaluate the production x-space path")
    if value["n_batches"] != fixture.config["bench_batches"]:
        raise ValueError("Bench did not complete all requested batches")
    if (
        type(value["total_samples"]) is not int
        or value["total_samples"] < value["n_batches"]
    ):
        raise ValueError("Bench did not complete a valid sample count")
    batches = value["batches"]
    if len(batches) != fixture.config["bench_batches"]:
        raise ValueError("Bench did not record all requested batches")
    batch_totals = [
        positive_number(batch["total"], "Bench batch Total") for batch in batches
    ]
    output = value["displayed_result"]
    numbers = [output["re"], output["im"]]
    if any(type(v) not in (int, float) or not math.isfinite(v) for v in numbers):
        raise ValueError("Bench returned a nonfinite physical result")
    if math.hypot(*numbers) == 0:
        raise ValueError("The physical benchmark must not evaluate to zero")
    totals = [row for row in value["summary"] if row["category"] == "Total"]
    if len(totals) != 1:
        raise ValueError("Bench must provide exactly one Total summary")
    timing = positive_number(totals[0]["mean_seconds_per_sample"], "Bench Total")
    if not math.isclose(timing, statistics.mean(batch_totals), rel_tol=1e-12):
        raise ValueError("Bench Total disagrees with its per-sample batch timings")
    return timing, numbers, value["total_samples"]


def measure(binary, repository, directory, fixture, deadline):
    directory.mkdir(parents=True, exist_ok=False)
    state = directory / "state"
    generation_card = directory / "generation.toml"
    generation_card.write_text(fixture.card)
    started = time.monotonic()
    run_cli(
        [str(binary), "--clean-state", "-s", str(state), str(generation_card)],
        repository,
        directory / "generation.log",
        deadline,
    )
    generation_seconds = time.monotonic() - started
    summaries = list(state.rglob("generation_summary.json"))
    if len(summaries) != 1:
        raise ValueError("Fresh generation must save one integrand summary")
    summary = json.loads(summaries[0].read_text())
    reports = summary["reports"]
    if len(reports) != 1 or reports[0]["graph_name"] != fixture.corpus["graph_name"]:
        raise ValueError("Fresh generation must build the requested graph")
    stats = reports[0]["stats"]
    count = stats["evaluator_count"]
    if type(count) is not int or count <= 0:
        raise ValueError("Fresh generation did not build any evaluators")
    generation_timings = {}
    for stage in (
        "total_time",
        "evaluator_spenso_time",
        "evaluator_symbolica_time",
        "evaluator_compile_time",
    ):
        duration = stats[stage]
        seconds, nanos = duration["secs"], duration["nanos"]
        if (
            type(seconds) is not int
            or not 0 <= seconds < 2**64
            or type(nanos) is not int
            or not 0 <= nanos < 10**9
        ):
            raise ValueError(f"Invalid native generation duration: {stage}")
        generation_timings[f"generation_{stage.removesuffix('_time')}_seconds"] = (
            seconds + nanos / 10**9
        )
    benchmark_card = fixture.benchmark_card(directory)
    run_cli(
        [str(binary), "--read-only-state", "-n", "-s", str(state), str(benchmark_card)],
        repository,
        directory / "benchmark.log",
        deadline,
    )
    timings, outputs, sample_counts = [], [], []
    for index in range(len(fixture.corpus["points"])):
        check_deadline(deadline)
        timing, output, samples = benchmark_result(
            directory / f"bench-{index}.json", fixture, index
        )
        timings.append(timing)
        outputs.extend(output)
        sample_counts.append(samples)
    check_deadline(deadline)
    return (
        {
            "generation_seconds": generation_seconds,
            "warm_runtime_seconds_per_sample": statistics.mean(timings),
        },
        outputs,
        {
            "samples_per_point": sample_counts,
            "evaluator_count": count,
            "generation_peak_ram_bytes": summary["peak_ram_bytes"],
            **generation_timings,
        },
    )


def measure_phase(
    arguments,
    fixture,
    phase,
    run_id,
    first_order,
    details,
    deadline,
    measurement,
    checkpoint,
):
    result = {"run_id": run_id, "pairs": []}
    measurement[phase] = result
    checkpoint()
    for index in range(fixture.config["pairs"]):
        order = first_order if index % 2 == 0 else first_order[::-1]
        pair = {"order": order}
        result["pairs"].append(pair)
        for side in order:
            label = "baseline" if side == "A" else "candidate"
            metrics, outputs, metadata = measure(
                getattr(arguments, label + "_binary"),
                arguments.input_root,
                arguments.work_dir / phase / str(index) / label,
                fixture,
                deadline,
            )
            pair[label] = metrics
            pair[label + "_output"] = outputs
            details.append({"phase": phase, "pair": index, "side": label, **metadata})
            checkpoint()
    return result


def run_case(args, suite, name):
    measurement = {}
    gate = PerformanceGate(args.baseline_commit, args.candidate_commit)
    output_ready = False
    try:
        args.output_dir.mkdir(parents=True, exist_ok=False)
        output_ready = True
        fixture = Fixture.from_suite(suite, name, args.input_root)
        measurement = {
            "schema_version": 1,
            "required_metrics": fixture.config["required_metrics"],
            "correctness_tolerance": fixture.config["correctness_tolerance"],
            **{
                side: {"commit": getattr(args, side + "_commit"), **fixture.hashes}
                for side in ("baseline", "candidate")
            },
        }
        details = []
        input_path = args.output_dir / "measurements.json"
        report_path = args.output_dir / "report.json"

        def checkpoint():
            check_deadline(deadline)
            write_json(input_path, measurement)
            write_json(args.output_dir / "measurement_metadata.json", details)

        # One deadline covers every discovery trial and any confirmation.
        # Executables are built before invoking this measurement harness.
        deadline = time.monotonic() + fixture.timeout
        measure_phase(
            args,
            fixture,
            "discovery",
            1,
            "AB",
            details,
            deadline,
            measurement,
            checkpoint,
        )
        check_deadline(deadline)
        report = gate.analyze(measurement)
        check_deadline(deadline)
        if report["status"] == 2 and report["confirmation_required"] == 1:
            measure_phase(
                args,
                fixture,
                "confirmation",
                2,
                "BA",
                details,
                deadline,
                measurement,
                checkpoint,
            )
            check_deadline(deadline)
            report = gate.analyze(measurement)
            check_deadline(deadline)
        write_json(report_path, report)
        gate.print_report(report)
        return report["status"]
    except (OSError, ValueError, KeyError, TypeError) as error:
        if output_ready:
            # The comparator owns numeric failure reports too. An incomplete
            # phase cannot become either a pass or a confirmed regression.
            write_json(args.output_dir / "measurements.json", measurement)
            report = gate.analyze(None)
            write_json(args.output_dir / "report.json", report)
            gate.print_report(report)
        print(f"Performance measurement failed: {error}", file=sys.stderr)
        return 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for side in ("baseline", "candidate"):
        parser.add_argument(f"--{side}-binary", type=Path, required=True)
        parser.add_argument(f"--{side}-commit", required=True)
    parser.add_argument(
        "--input-root",
        type=Path,
        required=True,
        help="Shared candidate repository supplying both revisions' inputs",
    )
    parser.add_argument(
        "--case", help="Run one manifest case; otherwise run the complete suite"
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        required=True,
        help="Private state/log directory, never uploaded",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Metadata and numeric report directory",
    )
    args = parser.parse_args()
    try:
        for side in ("baseline", "candidate"):
            if not re.fullmatch(r"[0-9a-f]{40}", getattr(args, side + "_commit")):
                raise ValueError("Commit identities must be complete hexadecimal SHA1s")
            binary = getattr(args, side + "_binary").resolve()
            if not binary.is_file() or not os.access(binary, os.X_OK):
                raise ValueError("Missing executable")
            setattr(args, side + "_binary", binary)
        args.input_root = args.input_root.resolve()
        if not args.input_root.is_dir():
            raise ValueError("Missing shared input repository")
        args.work_dir = args.work_dir.resolve()
        args.output_dir = args.output_dir.resolve()
        if args.work_dir.exists():
            raise ValueError(
                "Private work directory must be new; no generated state may be reused"
            )
        if args.work_dir.is_relative_to(
            args.output_dir
        ) or args.output_dir.is_relative_to(args.work_dir):
            raise ValueError(
                "Private state/logs must be separate from uploaded metadata"
            )
        if args.output_dir.exists():
            raise ValueError(
                "Uploaded metadata directory must be new; stale reports cannot be reused"
            )
        suite = tomllib.loads(
            (args.input_root / "benchmarks/performance/suite.toml").read_text()
        )
        cases = Fixture.cases(suite)
        names = [args.case] if args.case else list(cases)
        if any(name not in cases for name in names):
            raise ValueError("Unknown performance case")
        args.work_dir.mkdir(parents=True, exist_ok=False)
        args.output_dir.mkdir(parents=True, exist_ok=False)
        status = 0
        for name in names:
            print(f"Measuring {name}", flush=True)
            case_args = argparse.Namespace(
                **(
                    vars(args)
                    | {
                        "work_dir": args.work_dir / name,
                        "output_dir": args.output_dir / name,
                    }
                )
            )
            status = max(status, run_case(case_args, suite, name))
        return status
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"Performance suite failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
