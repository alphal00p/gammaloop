"""Native profile window for repeated fresh reductions of the matched quark input."""

import argparse
import gc
import hashlib
import json
import runpy
import statistics
import sys
import time
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--mode", choices=["aliases", "no-aliases"], required=True)
p.add_argument("--iterations", type=int, required=True)
p.add_argument("--save", type=Path, required=True)
p.add_argument("--hep-examples", type=Path, required=True)
p.add_argument("--compare", type=Path, required=True)
p.add_argument("--python-root", type=Path)
p.add_argument("--control", type=Path)
p.add_argument("--ack", type=Path)
p.add_argument("--materialize", action="store_true")
args = p.parse_args()
if args.python_root:
    sys.path.insert(0, str(args.python_root))
sys.argv = [
    "ladder.py",
    "quark",
    "--hep-examples",
    str(args.hep_examples),
    "--save",
    str(args.save.with_suffix(".warmup.bin")),
    "--compare",
    str(args.compare),
    "--symbol-indices",
    "--max-steps-per-domain",
    "128",
]
if args.mode == "no-aliases":
    sys.argv.append("--no-aliases")
ns = runpy.run_path(str(Path(__file__).with_name("ladder.py")))
sources = [ns["sp"].TensorExpression(ns["raw"]) for _ in range(args.iterations)]
results = []
times = []
control = args.control.open("w") if args.control else None
ack = args.ack.open("rb", buffering=0) if args.ack else None


def command(cmd):
    if control:
        control.write(cmd + "\n")
        control.flush()
        reply = ack.read(64)
        assert reply.strip(b"\0\r\n ") == b"ack", repr(reply)


gc.collect()
command("enable")
start = time.perf_counter()
cpu = time.process_time()
for source in sources:
    t = time.perf_counter()
    result = source.simplify_algebra(**ns["settings"])
    times.append(time.perf_counter() - t)
    results.append(result)
elapsed = time.perf_counter() - start
cpu = time.process_time() - cpu
command("disable")
# Retain results until after profiling, excluding destruction from the window.
assert all(r.reduction_status == ns["sp"].ReductionStatus.Complete for r in results)
counts = {len(r.aliases) for r in results}
assert counts == ({0} if args.mode == "no-aliases" else {576}), counts
assert results[-1].expand().to_expression() == ns["scalar"]
record = {
    "mode": args.mode,
    "iterations": args.iterations,
    "wall": elapsed,
    "cpu": cpu,
    "seconds_per_reduction": times,
    "median": statistics.median(times),
    "mean": statistics.mean(times),
    "aliases": sorted(counts),
    "complete": True,
    "exact_last_comparison": True,
    "extension": str(ns["core"].__file__),
    "extension_sha256": hashlib.sha256(
        Path(ns["core"].__file__).read_bytes()
    ).hexdigest(),
    "profiled": control is not None,
    "settings": ns["settings"],
}
if args.materialize:
    materialization = []
    reruns = []
    for result in results:
        start = time.perf_counter()
        materialized = result.expand().to_expression()
        materialization.append(time.perf_counter() - start)
        assert materialized == ns["scalar"]
        start = time.perf_counter()
        repeated = result.simplify_algebra(**ns["settings"])
        reruns.append(time.perf_counter() - start)
        assert repeated.reduction_status == ns["sp"].ReductionStatus.Complete
    record["materialization_seconds"] = materialization
    record["completed_rerun_seconds"] = reruns
    record["mean_reduction_plus_materialization"] = statistics.mean(
        times
    ) + statistics.mean(materialization)
args.save.write_text(json.dumps(record, indent=2) + "\n")
print("PROFILE_RESULT", json.dumps(record), flush=True)
