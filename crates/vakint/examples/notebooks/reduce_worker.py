"""Reduce one prepared ladder input; the command of a gallery ``ladder`` job.

Run by the job supervisor (``benchmark_job.py``) as
``reduce_worker.py WORK``: it reads ``WORK/request.json``'s ``spec``,
verifies the shipped input against its recorded sha256, reduces it with the
Ordered engine to retained masters, and writes

- ``retained.symbolica``: the result's canonical text;
- ``benchmark.json``: timings, data, closure, Ordered step-cache statistics
  and the comparison with the recorded fixture hash;
- ``stage.json``: the current coarse stage, for the notebook's progress card;
- ``monitor.json``: the native progress snapshot (Vakint's ``progress_snapshot``).

The settings are those of the rung test (tests/ordered_reduction_rungs.rs):
MSbar, retained masters, dot-product notation and L+1 epsilon orders, so the
retained sha256 must equal the recorded one for every arithmetic and n_cores.
"""

import hashlib
import json
import os
import resource
import sys
import traceback
from pathlib import Path
from time import perf_counter, process_time, time

from shared_gallery import Gallery


def sha256(text):
    return hashlib.sha256(text.encode()).hexdigest()


class Stages:
    """Wall-clock and CPU seconds per stage. CPU time (all threads) is immune
    to the process being paused, e.g. by a host memory guard."""

    def __init__(self, work):
        self.work = work
        self.seconds = {}
        self.cpu_seconds = {}
        self.name = None
        self.started = None
        self.cpu_started = None

    def enter(self, name):
        now, cpu = perf_counter(), process_time()
        if self.name is not None:
            self.seconds[self.name] = now - self.started
            self.cpu_seconds[self.name] = cpu - self.cpu_started
        self.name, self.started, self.cpu_started = name, now, cpu
        temporary = self.work / "stage.tmp"
        temporary.write_text(json.dumps({"stage": name, "since_unix_seconds": time()}))
        temporary.replace(self.work / "stage.json")


def catalogue_closure(result, engine, E):
    """Catalogue masters present in the result, and whether any other remain."""
    present = []
    for master in engine.master_symbols():
        expression = master.expression
        if next(iter(result.match(expression)), None) is not None:
            present.append(expression)
    rest = result
    for expression in present:
        rest = rest.replace(expression, E("0"))
    unknown = next(iter(rest.match(E("vakint::RustRedMaster(x__)"))), None)
    return {"closed": unknown is None, "master_count": len(present)}


def reduce(work, spec, stages):
    cases = Gallery.ladder_cases()
    case = cases[spec["case"]]
    loops = case["loop_count"]
    receipt = {
        "schema": "vakint.gallery.ladder-run.v1",
        "case": spec["case"],
        "label": case["label"],
        "loop_count": loops,
        "projection": case["projection"],
        "n_cores": spec["n_cores"],
        "ordered_weight_arithmetic": spec["arithmetic"],
        "strategy": "Ordered",
        "normalization": "MSbar",
        "process": {"pid": os.getpid(), "cpus": sorted(os.sched_getaffinity(0))},
    }

    stages.enter("reading the prepared input")
    path = Gallery.ladder_input(spec["case"])
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != case["input"]["sha256"]:
        raise ValueError(f"{path.name}: sha256 differs from the recorded input")
    receipt["input"] = {
        "file": path.name,
        "bytes": len(payload),
        "sha256_verified": True,
    }

    stages.enter("importing the native host")
    from symbolica import E

    stages.enter("creating the engine")
    cache_directory = spec.get("ordered_cache_directory")
    engine = Gallery.engine(
        mode="retain",
        n_cores=spec["n_cores"],
        ordered_weight_arithmetic=spec["arithmetic"],
        normalization="MSbar",
        use_dot_product_notation=True,
        number_of_terms_in_epsilon_expansion=loops + 1,
        max_elapsed_seconds=float(spec["max_seconds"]),
        ordered_cache="read-write",
        ordered_cache_directory=cache_directory,
    )
    receipt["execution"] = Gallery.execution(engine)

    data_directory = spec.get("data_directory")
    if data_directory:
        stages.enter(f"loading {loops}L external data")
        loaded = engine.load_official_external_data(
            loops, version="latest", local_artifacts_path=data_directory
        )
        receipt["data"] = {"version": loaded["version"], "path": str(loaded["path"])}
    else:
        receipt["data"] = {"source": "built-in 1L-3L registry"}

    stages.enter("parsing and binding the dimension")
    source = E(payload.decode(), default_namespace="qcd_ladder")
    if sha256(source.to_canonical_string()) != case["input"]["sha256"]:
        raise ValueError("The parsed input does not reproduce its canonical text")
    binding = case["dimension_binding"]
    evaluated = source.replace(E(binding["source"]), E(binding["target"]))
    if sha256(evaluated.to_canonical_string()) != case["evaluated_input_sha256"]:
        raise ValueError("The dimension-bound input differs from the recorded one")
    receipt["evaluated_input_sha256_verified"] = True

    stages.enter("Ordered reduction")
    result = engine.evaluate(
        evaluated,
        progress_snapshot=str(work / "monitor.json"),
        progress_interval_seconds=0.5,
    )

    stages.enter("writing the retained result")
    text = result.to_canonical_string()
    (work / "retained.symbolica").write_text(text)
    digest = sha256(text)
    receipt["retained"] = {
        "file": "retained.symbolica",
        "bytes": len(text.encode()),
        "sha256": digest,
        "recorded_sha256": case["retained_sha256"],
        "matches_recorded": digest == case["retained_sha256"],
    }

    stages.enter("checking catalogue closure")
    receipt["catalogue_closure"] = catalogue_closure(result, engine, E)
    if cache_directory:
        receipt["ordered_cache"] = {
            "directory": cache_directory,
            "statistics": engine.ordered_cache_statistics(),
        }
    return receipt


def main():
    work = Path(sys.argv[1]).resolve()
    stages = Stages(work)
    started = perf_counter()
    spec = json.loads((work / "request.json").read_text())["spec"]
    try:
        receipt = reduce(work, spec, stages)
        status = 0
    except Exception as error:  # noqa: BLE001 - every failure is reported.
        traceback.print_exc()
        receipt = {
            "schema": "vakint.gallery.ladder-run.v1",
            "case": spec.get("case"),
            "error": f"{type(error).__name__}: {error}",
            "failed_stage": stages.name,
        }
        status = 1
    stages.enter("done" if status == 0 else f"failed while {stages.name}")
    seconds = {
        key: round(value, 4)
        for key, value in stages.seconds.items()
        if not key.startswith("failed while")
    }
    receipt["seconds"] = {**seconds, "total": round(perf_counter() - started, 4)}
    receipt["cpu_seconds"] = {
        **{key: round(value, 4) for key, value in stages.cpu_seconds.items()},
        "total": round(process_time(), 4),
    }
    receipt["peak_rss_bytes"] = (
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
    )
    temporary = work / "benchmark.tmp"
    temporary.write_text(json.dumps(receipt, indent=1, default=str))
    temporary.replace(work / "benchmark.json")
    return status


if __name__ == "__main__":
    sys.exit(main())
