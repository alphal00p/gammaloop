"""Record the existing native tensor benchmarks without mixing build and run clocks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import resource
import shutil
import subprocess
import time
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_lines(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def run(command: list[str], directory: Path) -> dict:
    directory.mkdir(parents=True)
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter_ns()
    with (
        (directory / "stdout.log").open("w") as stdout,
        (directory / "stderr.log").open("w") as stderr,
    ):
        process = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
    elapsed = time.perf_counter_ns() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt = {
        "command": command,
        "exit_code": process.returncode,
        "process_wall_ns": elapsed,
        "child_cpu_ns": round(
            (after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime) * 1e9
        ),
        "cumulative_child_peak_rss_kib": after.ru_maxrss,
        "stdout_sha256": sha256(directory / "stdout.log"),
        "stderr_sha256": sha256(directory / "stderr.log"),
    }
    (directory / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binaries", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--form", type=Path, required=True)
    parser.add_argument("--profile", default="dev-optim")
    parser.add_argument("--samples", type=int, default=3)
    parser.add_argument("--target-ms", type=float, default=8)
    parser.add_argument("--form-max-length", type=int, default=12)
    args = parser.parse_args()
    if args.samples < 1 or args.target_ms <= 0:
        parser.error("samples and target-ms must be positive")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    args.source = args.source.resolve()
    binaries = {
        name: (args.binaries / name).resolve()
        for name in (
            "metric_contraction_benchmark",
            "contraction_phase_benchmark",
            "trace_form",
        )
    }
    report = {
        "profile": args.profile,
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "cpu_affinity": sorted(os.sched_getaffinity(0)),
        "settings": {"samples": args.samples, "target_ms": args.target_ms},
        "binary_sha256": {name: sha256(path) for name, path in binaries.items()},
        "lockfile_sha256": sha256(args.source / "Cargo.lock"),
        "form_sha256": sha256(args.form),
        "boundaries": {
            "native_rows": "Existing benchmark method clocks; see example source.",
            "process_receipts": "Include process startup, input setup, checks and output.",
            "form_rows": "Process wall clocks include parsing, sorting and scalar checks.",
            "memory": "Cumulative child maximum, not isolated per-operation allocations.",
            "compilation": "Must complete separately before invoking this runner.",
            "diagnostic": (
                "Separate instrumented call after primary clocks; thread-local wall "
                "time, nested inclusive totals overlap. outermost_ns avoids "
                "same-phase recursive double counting. Absent on frozen baseline."
            ),
        },
    }
    metric = args.output / "metric"
    report["metric_run"] = run(
        [
            str(binaries["metric_contraction_benchmark"]),
            str(metric / "fixtures"),
            str(args.samples),
            str(args.target_ms),
        ],
        metric,
    )
    report["metric_rows"] = json_lines(metric / "stdout.log")
    report["phase_runs"] = []
    report["phase_rows"] = []
    # Keep terminal-rerun assertions and their failures in a separate complete run.
    full = args.output / "full_phases"
    report["full_phase_run"] = run(
        [
            str(binaries["contraction_phase_benchmark"]),
            str(metric / "fixtures"),
            str(full / "outputs"),
            str(args.samples),
            str(args.target_ms),
        ],
        full,
    )
    # Each original source gets a receipt, including refused inputs. A refusal
    # cannot hide the remaining cases or replace the failed complete cohort.
    for source in sorted((metric / "fixtures").glob("*.gamma_input")):
        case = args.output / "phases" / source.stem
        inputs = case / "inputs"
        inputs.mkdir(parents=True)
        shutil.copyfile(source, inputs / source.name)
        result = run(
            [
                str(binaries["contraction_phase_benchmark"]),
                str(inputs),
                str(case / "outputs"),
                str(args.samples),
                str(args.target_ms),
            ],
            case / "run",
        )
        result.update(case=source.stem, input_sha256=sha256(source))
        report["phase_runs"].append(result)
        report["phase_rows"].extend(
            {"source_case": source.stem, **row}
            for row in json_lines(case / "run" / "stdout.log")
        )
    form = args.output / "form"
    report["form_run"] = run(
        [
            str(binaries["trace_form"]),
            str(args.form.resolve()),
            str(args.form_max_length),
        ],
        form,
    )
    with (form / "stdout.log").open() as source:
        report["form_rows"] = list(csv.DictReader(source))
    for line in (form / "stderr.log").read_text().splitlines():
        if line.startswith("FORM programs and output: "):
            shutil.copytree(Path(line.split(": ", 1)[1]), form / "programs")
    (args.output / "results.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
