#!/usr/bin/env python3
"""Generate author targets independently of the native Figure15 scan."""

import argparse
import concurrent.futures
import decimal
import json
import math
import os
import pathlib
import queue
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--workers", type=int, default=6)
parser.add_argument(
    "--cpus",
    default="236:242",
    help="comma-separated CPUs or half-open ranges, e.g.236:242",
)
parser.add_argument("--deadline-epoch", type=float)
parser.add_argument("--output-directory", type=pathlib.Path, required=True)
parser.add_argument("--cache-directory", type=pathlib.Path)
parser.add_argument("--binary", type=pathlib.Path)
parser.add_argument(
    "--pilot-targets", type=pathlib.Path, default=ROOT / "pilot-targets.json"
)
parser.add_argument("--cli", type=pathlib.Path, default=ROOT / "reference_cli.py")
parser.add_argument("--dry-run", action="store_true")
args = parser.parse_args()
cpus = []
try:
    for item in args.cpus.split(","):
        if ":" in item:
            lower, upper = map(int, item.split(":"))
            cpus.extend(range(lower, upper))
        else:
            cpus.append(int(item))
except ValueError:
    parser.error("invalid --cpus")
if (
    len(set(cpus)) != len(cpus)
    or not 1 <= args.workers <= 6
    or args.workers > len(cpus)
):
    parser.error("require1..6 workers and at least that many distinct CPUs")
if not set(cpus).issubset(os.sched_getaffinity(0)):
    parser.error("requested CPUs are outside inherited affinity")
if args.deadline_epoch is not None and not math.isfinite(args.deadline_epoch):
    parser.error("deadline must be finite")
output = args.output_directory.resolve()
output.mkdir(parents=True, exist_ok=True)
cache = (args.cache_directory or output / "cache").resolve()
points = {
    (E, k)
    for k in (0.0, 0.5, 1.0, 1.5)
    for E in (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
    if E != k
}
for k in (0.5, 1.0, 1.5):
    for sign in (-1, 1):
        E = float(decimal.Decimal(str(k)) + sign * decimal.Decimal(".05"))
        points.add((E, k))
points = sorted(points, key=lambda point: (point[1], point[0]))
assert len(points) == 35
rows = {}
statuses = {f"{E:g},{k:g}": {"E": E, "k": k, "status": "pending"} for E, k in points}
started = time.time()


def atomic(path, value):
    temporary = path.with_suffix(f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def retain(E, k, value, source):
    for component in ("T", "L"):
        lo = value.get("rho" + component + "_lo")
        nlo = value.get("rho" + component + "_nlo")
        if (
            lo is not None
            and nlo is not None
            and math.isfinite(lo)
            and math.isfinite(nlo)
        ):
            rows[(E, k, component)] = {
                "E": E,
                "k": k,
                "component": component,
                "rho_lo": lo,
                "rho_nlo": nlo,
                "reference_tolerance": value["tolerance"],
                "rho_g3": lo + 9 * nlo,
                "source": source,
                "method": value.get("method", "author assembly"),
                "source_archive_sha256": value.get(
                    "source_archive_sha256",
                    "207a2cb3d39933a57a5c1a8592d664c7dc9d5750eb8cd10801e98ac9560446fb",
                ),
            }


if args.pilot_targets.exists():
    for value in json.loads(args.pilot_targets.read_text()):
        E, k = float(value["E"]), float(value["k"])
        if (E, k) in points:
            retain(E, k, value, "pilot")


def persist(phase):
    atomic(
        output / "targets.json",
        [rows[key] for key in sorted(rows, key=lambda key: (key[1], key[0], key[2]))],
    )
    atomic(
        output / "status.json",
        {
            "phase": phase,
            "started_epoch": started,
            "updated_epoch": time.time(),
            "deadline_epoch": args.deadline_epoch,
            "workers": args.workers,
            "cpus": cpus[: args.workers],
            "points": len(points),
            "available_target_rows": len(rows),
            "jobs": list(statuses.values()),
            "process_group": os.getpgrp(),
            "cache_directory": str(cache),
        },
    )


persist("prepared" if args.dry_run else "running")
if args.dry_run:
    sys.exit(0)
leases = queue.Queue()
for cpu in cpus[: args.workers]:
    leases.put(cpu)
thread_environment = os.environ.copy()
for name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
):
    thread_environment[name] = "1"


def evaluate(point):
    E, k = point
    cpu = leases.get()
    begin = time.time()
    label = f"E{E:g}-k{k:g}"
    result_path = output / f"{label}.json"
    timeout = (
        1800.0
        if args.deadline_epoch is None
        else min(1800.0, args.deadline_epoch - begin)
    )
    if timeout <= 0:
        leases.put(cpu)
        return {
            "E": E,
            "k": k,
            "status": "deadline",
            "cpu": cpu,
            "started_epoch": begin,
            "seconds": 0.0,
        }
    command = [
        sys.executable,
        str(args.cli),
        "--energy",
        str(E),
        "--momentum",
        str(k),
        "--tolerance",
        "1e-5",
        "--basis",
        "both",
        "--timeout",
        str(timeout),
        "--cache-directory",
        str(cache),
        "--cpu",
        str(cpu),
        "--output",
        str(result_path),
    ]
    if args.binary:
        command.extend(["--binary", str(args.binary)])
    try:
        # Inherit the controller process group; hard cutoff covers every child.
        with (
            (output / f"{label}.stdout").open("w") as stdout,
            (output / f"{label}.stderr").open("w") as stderr,
        ):
            process = subprocess.Popen(
                command, stdout=stdout, stderr=stderr, env=thread_environment
            )
            try:
                process.wait(timeout=timeout + 5)
            except BaseException:
                # The CLI handles SIGTERM by killing and reaping its source child.
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                raise
        value = json.loads(result_path.read_text()) if result_path.exists() else {}
        status = (
            "complete"
            if process.returncode == 0 and value.get("status") == "ok"
            else "failed"
        )
        return {
            "E": E,
            "k": k,
            "status": status,
            "cpu": cpu,
            "started_epoch": begin,
            "seconds": time.time() - begin,
            "returncode": process.returncode,
            "command": command,
            "result": value,
        }
    except (OSError, ValueError, subprocess.TimeoutExpired) as error:
        return {
            "E": E,
            "k": k,
            "status": "failed",
            "cpu": cpu,
            "started_epoch": begin,
            "seconds": time.time() - begin,
            "command": command,
            "error": str(error),
        }
    finally:
        leases.put(cpu)


with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
    futures = {executor.submit(evaluate, point): point for point in points}
    for future in concurrent.futures.as_completed(futures):
        record = future.result()
        E, k = record["E"], record["k"]
        if record["status"] == "complete":
            retain(E, k, record["result"], "reference-grid")
        statuses[f"{E:g},{k:g}"] = {
            key: value for key, value in record.items() if key != "result"
        }
        with (output / "jobs.jsonl").open("a") as log:
            log.write(json.dumps(record) + "\n")
        persist("running")
persist("complete")
