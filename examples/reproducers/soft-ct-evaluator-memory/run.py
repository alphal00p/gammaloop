#!/usr/bin/env python3
"""Compare evaluator construction in separate, bounded processes.

Usage: run.py BINARY INPUT_JSON OUTPUT_DIRECTORY [--perf /path/to/perf]
The input is the existing evaluator_symbolica_input_dump payload. Only common
pair elimination differs between runs; input and binary hashes are recorded.
"""

import argparse
import datetime
import hashlib
import json
import os
import signal
import subprocess
import time
from decimal import Decimal, localcontext
from pathlib import Path


def identify(path):
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def write(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def stop(process, sig=signal.SIGTERM):
    if process.poll() is None:
        os.killpg(process.pid, sig)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()


def run(args, mode, identity):
    folder = args.output / mode
    folder.mkdir()
    command = [str(args.binary), str(args.input), mode]
    receipt = {
        "command": command,
        "identity": identity,
        "mode": mode,
        "started_utc": datetime.datetime.now(datetime.UTC).isoformat(),
        "wall_limit_seconds": args.seconds,
        "rss_limit_gib": args.gib,
        "sample_interval_seconds": 0.2,
        "profiled": args.perf is not None,
    }
    started = time.monotonic()
    receipt["started_monotonic_seconds"] = started
    peak = hwm = 0
    reason = None
    error = None
    profiler = None
    stage = None
    events = []
    env = dict(os.environ, RAYON_NUM_THREADS="1", RUST_MIN_STACK=str(128 * 1024**2))
    with (
        (folder / "console.jsonl").open("w") as log,
        (folder / "console.jsonl").open() as reader,
        (folder / "resources.jsonl").open("w") as samples,
        (folder / "perf.log").open("w") as perf_log,
    ):
        process = subprocess.Popen(
            command,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            start_new_session=True,
        )
        receipt["pid"] = process.pid
        write(folder / "active.json", receipt)
        print(
            json.dumps({"stage": "started", "mode": mode, "pid": process.pid}),
            flush=True,
        )
        if args.perf:
            perf_command = [
                str(args.perf),
                "record",
                "-e",
                "cycles:u",
                "-F",
                "49",
                "--call-graph",
                "dwarf,8192",
                "-o",
                str(folder / "perf.data"),
                "-p",
                str(process.pid),
            ]
            receipt["perf_command"] = perf_command
            profiler = subprocess.Popen(
                perf_command,
                stdout=perf_log,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        pending = ""
        try:
            while process.poll() is None:
                elapsed = time.monotonic() - started
                lines = (pending + reader.read()).split("\n")
                pending = lines.pop()
                for line in lines:
                    try:
                        event = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    events.append(event)
                    stage = event.get("stage", stage)
                try:
                    status = dict(
                        line.split(":", 1)
                        for line in Path(f"/proc/{process.pid}/status")
                        .read_text()
                        .splitlines()
                        if ":" in line
                    )
                    rss = int(status.get("VmRSS", "0 kB").split()[0])
                    hwm = max(hwm, int(status.get("VmHWM", "0 kB").split()[0]))
                except FileNotFoundError:
                    rss = 0
                peak = max(peak, rss)
                samples.write(
                    json.dumps(
                        {
                            "elapsed_seconds": elapsed,
                            "rss_kib": rss,
                            "stage": stage,
                        }
                    )
                    + "\n"
                )
                samples.flush()
                if rss > args.gib * 1024**2:
                    reason = "RSS limit"
                    break
                if elapsed > args.seconds:
                    reason = "wall limit"
                    break
                try:
                    process.wait(timeout=0.2)
                except subprocess.TimeoutExpired:
                    pass
        except (KeyboardInterrupt, OSError, ValueError) as exc:
            reason = "supervisor interrupted or failed"
            error = repr(exc)
        finally:
            stop(process)
            if profiler:
                stop(profiler, signal.SIGINT)
    # Read again after exit to include the final event and check complete lines.
    events = []
    for line in (folder / "console.jsonl").read_text().splitlines():
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError:
            pass
    receipt.update(
        elapsed_seconds=time.monotonic() - started,
        exit_code=process.returncode,
        termination_reason=reason,
        error=error,
        sampled_peak_rss_kib=peak,
        peak_process_hwm_kib=hwm,
        events=events,
        profiler_exit_code=profiler.returncode if profiler else None,
        input_unchanged=identify(args.input) == identity["input"],
        binary_unchanged=identify(args.binary) == identity["binary"],
    )
    write(folder / "receipt.json", receipt)
    print(
        json.dumps(
            {
                k: receipt[k]
                for k in (
                    "mode",
                    "elapsed_seconds",
                    "exit_code",
                    "termination_reason",
                    "sampled_peak_rss_kib",
                    "peak_process_hwm_kib",
                )
            }
        ),
        flush=True,
    )
    return receipt


def compare(results):
    points = [
        {
            e["point_index"]: e
            for e in result["events"]
            if e.get("stage") == "numeric_point"
        }
        for result in results
    ]
    if len(points) != 2 or any(set(p) != {0, 1, 2} for p in points):
        return {"passed": False, "reason": "incomplete numeric checks"}
    residuals = []
    with localcontext() as context:
        context.prec = 350
        for index in range(3):
            left, right = (p[index] for p in points)
            if not all(p["finite"] and p["nonzero"] for p in (left, right)):
                return {"passed": False, "reason": "nonfinite or all-zero output"}
            for a, b in zip(left["outputs"], right["outputs"], strict=True):
                a, b = ([Decimal(x) for x in pair] for pair in (a, b))
                scale = max(abs(x) for x in a + b)
                error = max(abs(x - y) for x, y in zip(a, b, strict=True))
                residuals.append(error / scale if scale else Decimal(0))
    maximum = max(residuals)
    return {
        "passed": maximum < Decimal("1e-250"),
        "maximum_relative_residual": str(maximum),
        "relative_tolerance": "1e-250",
        "normalization": "complex infinity norm, symmetric between modes",
        "point_count": 3,
        "precision_bits": 1000,
        "scope": "synthetic complex parameter vectors, not physical momenta",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--perf", type=Path)
    parser.add_argument("--seconds", type=int, default=900)
    parser.add_argument("--gib", type=int, default=100)
    args = parser.parse_args()

    def interrupted(signum, _frame):
        raise KeyboardInterrupt(f"signal {signum}")

    signal.signal(signal.SIGTERM, interrupted)
    args.binary, args.input, args.output = (
        path.resolve() for path in (args.binary, args.input, args.output)
    )
    args.output.mkdir(parents=True, exist_ok=False)
    identity = {
        "input": identify(args.input),
        "binary": identify(args.binary),
        "supervisor": identify(Path(__file__).resolve()),
    }
    write(args.output / "identity.json", identity)
    results = []
    for mode in ["off", "default"]:
        results.append(run(args, mode, identity))
        write(args.output / "results.json", results)
        if results[-1]["error"]:
            break
    comparison = compare(results)
    write(args.output / "comparison.json", comparison)
    print(json.dumps(comparison), flush=True)
    return int(
        not comparison["passed"]
        or any(r["exit_code"] != 0 or r["termination_reason"] for r in results)
    )


if __name__ == "__main__":
    raise SystemExit(main())
