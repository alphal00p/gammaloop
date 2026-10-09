#!/usr/bin/env python3
"""Instruction and allocation gates for Feynman-diagram generation.

Runs the workloads of crates/feynkit-generator/examples/generation_workloads.rs
under valgrind and attributes their instructions (callgrind Ir) to the generation
stages. Each workload runs in two modes:

  generation  cancelled at the first progress report after enumeration: stage-1
              topologies, stage-2 dressing, graph materialization and certificates,
              and feynkit's certificate sort, exactly as feynkit runs them;
  full        the whole pipeline, adding interaction filters, numerators and
              grouping.

Commands:

  measure [--out DIR] [--jobs N] [--runs N] [--modes generation,full]
          [--allocations] [--label TEXT] [CASE ...]
      Build the example in the `generation-ir` profile (release without LTO, 16
      codegen units and line tables, defined on the command line), run each case
      (default: the reference workloads) in each mode `--runs` times under
      callgrind, optionally once under DHAT for allocation counts, and write
      DIR/summary.json. The generation mode of the reference workloads takes
      minutes; the full pipeline of gg->gg at three loops runs for about an hour
      and a half under callgrind and twice that under DHAT.
  summarize FILE ...
      Print the stage table of callgrind outputs or of summary.json files.
  compare OLD.json NEW.json
      Print new/old ratios of the total and of every stage, per workload and mode.

Ir varies between runs of one binary: hash maps take random keys at run time (the
ahash maps of linnet's stage-2 feasibility checks dominate), which moves the totals
of dressing-heavy workloads by up to two percent (the photon self-energy at five
loops) and the dressing stage by several percent. Totals and stages are therefore
medians over the runs, with the totals of every run kept; compare medians of
several runs of builds with the same toolchain and flags. The [min, max] of a
baseline's runs is not a tolerance band: a run of the unchanged build can fall
outside it, so gates judge the ratio of the medians. A stage is the inclusive cost of
the functions its pattern matches, summed over every source file their inlined
code comes from; calls between matched functions count once, and recursive calls
of a function to itself are ignored. A stage whose functions are all inlined away
is reported as missing, never as zero. Work run on a rayon pool is attributed on
the worker thread, so the stages are not nested in `Generator::generate`; `total`
is the whole program, including the model build.
"""

import argparse
import collections
import concurrent.futures
import json
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROFILE = "generation-ir"
PROFILE_CONFIG = [
    f'profile.{PROFILE}.inherits="release"',
    f"profile.{PROFILE}.lto=false",
    f"profile.{PROFILE}.codegen-units=16",
    f'profile.{PROFILE}.debug="line-tables-only"',
]

LINNET = r"<linnet::half_edge::generate"
FEYNKIT = r"<feynkit_generator::generation"
# (key, description, pattern of the function names)
STAGES = [
    (
        "model",
        "model build",
        (
            r"^<feynkit_model::model::Model>::"
            r"(standard_model|qcd|qed|electroweak|qcd_qed|yang_mills|phi3|phi4|phi_3_4)$"
        ),
    ),
    (
        "enumeration",
        "linnet generation (Diagrams::next)",
        rf"^{LINNET}::Diagrams<.*> as core::iter::traits::iterator::Iterator>::next$",
    ),
    ("topologies", "  stage 1: topologies", rf"^{LINNET}::topology::Search>::next$"),
    (
        "topology_canonize",
        "    canonical augmentation",
        rf"^{LINNET}::topology::Search>::canonize$",
    ),
    ("dressing_start", "  stage 2: start", rf"^{LINNET}::decorate::Decorator>::start$"),
    ("dressing", "  stage 2: dressing", rf"^{LINNET}::decorate::Decorator>::next$"),
    (
        "filter",
        "  filter (process filters)",
        rf"^{LINNET}::view::Pruner<.*> as {LINNET[1:]}::view::Prune>::keeps$",
    ),
    (
        "materialize",
        "  HedgeGraph materialization",
        r"^<core::cell::once::OnceCell<linnet::half_edge::HedgeGraph<.*>>>::try_init",
    ),
    (
        "certify",
        "  certificates",
        r"^<linnet::half_edge::canon::Workspace>::certify::<",
    ),
    (
        "sort",
        "certificate sort",
        rf"^<\[{LINNET[1:]}::Diagram<.*>\]>::sort_by::<{FEYNKIT}::Generator>::generate",
    ),
    (
        "drop",
        "drop of the generated diagrams",
        r"^core::ptr::drop_glue::<alloc::vec::Vec<linnet::half_edge::generate::Diagram<",
    ),
    (
        "interaction_filters",
        "interaction filters",
        rf"^{FEYNKIT}::Generator>::diagram_input$",
    ),
    (
        "canonical_frame",
        "  numerator-frame canonicalization",
        rf"^{FEYNKIT}::ResolvedProcess>::canonicalize_numerator_graph_pass$",
    ),
    ("numerators", "numerators (to_diagram)", rf"^{FEYNKIT}::Generator>::to_diagram$"),
    ("grouping", "grouping", r"^feynkit_generator::grouping::group_diagrams$"),
]


def parse_callgrind(path):
    """Self costs per function and inclusive costs per call edge, aggregated by name."""
    names = {}
    self_cost = collections.Counter()
    edges = collections.Counter()
    function = callee = None
    call = False
    total = None

    def name(text):
        match = re.match(r"\((\d+)\)(?: (.*))?$", text)
        if not match:
            return text
        if match.group(2) is not None:
            names[match.group(1)] = match.group(2)
        return names[match.group(1)]

    with open(path) as lines:
        for line in lines:
            if not line.strip():
                continue
            first = line[0]
            if first.isdigit() or first in "+-*":
                parts = line.split()
                cost = int(parts[1]) if len(parts) > 1 else 0
                if call:
                    edges[(function, callee)] += cost
                    call = False
                else:
                    self_cost[function] += cost
            elif line.startswith("fn="):
                function = name(line[3:].rstrip("\n"))
            elif line.startswith("cfn="):
                callee = name(line[4:].rstrip("\n"))
            elif line.startswith("calls="):
                call = True
            elif line.startswith(("summary:", "totals:")):
                total = int(line.split()[1])
    return total, self_cost, edges


def stage_costs(path):
    total, self_cost, edges = parse_callgrind(path)
    inclusive = collections.Counter(self_cost)
    for (caller, callee), cost in edges.items():
        if caller != callee:
            inclusive[caller] += cost
    functions = set(self_cost) | {callee for _, callee in edges}
    stages = {}
    for key, _, pattern in STAGES:
        regex = re.compile(pattern)
        matched = {f for f in functions if regex.search(f)}
        if not matched:
            stages[key] = None
            continue
        internal = sum(
            cost
            for (caller, callee), cost in edges.items()
            if caller != callee and caller in matched and callee in matched
        )
        stages[key] = sum(inclusive[f] for f in matched) - internal
    return {"total": total, "stages": stages}


def table(rows):
    """rows: list of (label, value-or-None, total)."""
    out = []
    for label, value, total in rows:
        if value is None:
            out.append(f"  {label:40} {'missing':>17}")
        else:
            out.append(f"  {label:40} {value:>17,} {100 * value / total:6.2f}%")
    return "\n".join(out)


def describe(measurement):
    """The stage table of one callgrind output, or of the medians of several runs."""
    total = measurement["total"]
    label = "total (median)" if "runs" in measurement else "total"
    rows = [(label, total, total)]
    rows += [
        (description, measurement["stages"].get(key), total)
        for key, description, _ in STAGES
    ]
    text = table(rows)
    if "runs" in measurement:
        runs = ", ".join(f"{run:,}" for run in measurement["runs"])
        text += f"\n  {'total of every run':40} {runs}"
    if "allocations" in measurement:
        allocations = measurement["allocations"]
        text += (
            f"\n  {'allocations':40} {allocations['blocks']:>17,} blocks, "
            f"{allocations['bytes']:,} bytes, peak {allocations['peak_bytes']:,} bytes"
        )
    return text


def cargo_target():
    metadata = subprocess.run(
        ["cargo", "metadata", "--format-version", "1", "--no-deps"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return Path(json.loads(metadata.stdout)["target_directory"])


def build():
    command = [
        "cargo",
        "build",
        "-p",
        "feynkit-generator",
        "--example",
        "generation_workloads",
    ]
    command += ["--profile", PROFILE]
    for setting in PROFILE_CONFIG:
        command += ["--config", setting]
    subprocess.run(command, cwd=ROOT, check=True)
    return cargo_target() / PROFILE / "examples" / "generation_workloads"


def valgrind(tool, out_file, log, workload):
    """Run a workload command under a valgrind tool writing `out_file`, and keep its log."""
    process = subprocess.run(
        ["valgrind", f"--tool={tool}", f"--{tool}-out-file={out_file}", *workload],
        capture_output=True,
        text=True,
        check=False,
    )
    log.write_text(process.stdout + process.stderr)
    if process.returncode != 0:
        raise RuntimeError(f"{' '.join(workload)} failed under {tool}, see {log}")
    return process


def count(text):
    return int(text.replace(",", ""))


def callgrind(binary, out, mode, case, index):
    stem = f"{mode}.{case}.{index}"
    output = out / f"{stem}.callgrind.out"
    log = out / f"{stem}.callgrind.log"
    workload = [str(binary), mode, case]
    process = valgrind("callgrind", output, log, workload)
    measurement = stage_costs(output)
    result = re.search(r"^RESULT (.*)$", process.stdout, re.MULTILINE).group(1)
    measurement["result"] = dict(item.split("=", 1) for item in result.split())
    return measurement


def dhat(binary, out, mode, case):
    stem = f"{mode}.{case}"
    output = out / f"{stem}.dhat.json"
    log = out / f"{stem}.dhat.log"
    process = valgrind("dhat", output, log, [str(binary), mode, case])
    total = re.search(r"Total:\s+([\d,]+) bytes in ([\d,]+) blocks", process.stderr)
    peak = re.search(r"At t-gmax:\s+([\d,]+) bytes", process.stderr)
    return {
        "bytes": count(total.group(1)),
        "blocks": count(total.group(2)),
        "peak_bytes": count(peak.group(1)),
    }


def median(runs):
    """The medians of the totals and stages of several runs of one workload."""
    results = {json.dumps(run["result"], sort_keys=True) for run in runs}
    if len(results) != 1:
        raise RuntimeError(f"runs report different counts: {results}")
    stages = {}
    for key, _, _ in STAGES:
        values = [run["stages"][key] for run in runs]
        stages[key] = None if None in values else statistics.median_low(values)
    totals = [run["total"] for run in runs]
    return {
        "total": statistics.median_low(totals),
        "runs": totals,
        "stages": stages,
        "result": runs[0]["result"],
    }


def measure(arguments):
    for tool in ("cargo", "valgrind"):
        if shutil.which(tool) is None:
            sys.exit(f"{tool} is required")
    binary = build()
    cases = (
        arguments.cases
        or subprocess.run(
            [str(binary), "list", "reference"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.split()
    )
    modes = arguments.modes.split(",")
    out = Path(arguments.out)
    out.mkdir(parents=True, exist_ok=True)

    def version(tool):
        return subprocess.run(
            [tool, "--version"], check=True, capture_output=True, text=True
        ).stdout.strip()

    summary = {
        "label": arguments.label,
        "toolchain": version("rustc"),
        "valgrind": version("valgrind"),
        "profile": PROFILE_CONFIG,
        "workloads": {case: {} for case in cases},
    }
    workloads = [(mode, case) for mode in modes for case in cases]
    with concurrent.futures.ThreadPoolExecutor(arguments.jobs) as pool:
        runs = {
            workload: [
                pool.submit(callgrind, binary, out, *workload, index)
                for index in range(arguments.runs)
            ]
            for workload in workloads
        }
        allocations = {
            workload: pool.submit(dhat, binary, out, *workload)
            for workload in workloads
            if arguments.allocations
        }
        for mode, case in workloads:
            measurement = median([run.result() for run in runs[(mode, case)]])
            if arguments.allocations:
                measurement["allocations"] = allocations[(mode, case)].result()
            summary["workloads"][case][mode] = measurement
            print(f"{case} {mode}\n{describe(measurement)}", flush=True)
    path = out / "summary.json"
    path.write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
    print(f"wrote {path}")


def summarize(arguments):
    for file in arguments.files:
        if file.endswith(".json"):
            summary = json.loads(Path(file).read_text())
            print(f"{file}: {summary.get('label')}")
            for case, modes in summary["workloads"].items():
                for mode, measurement in modes.items():
                    print(f"{case} {mode}\n{describe(measurement)}")
        else:
            print(f"{file}\n{describe(stage_costs(file))}")


def compare(arguments):
    old, new = (
        json.loads(Path(file).read_text()) for file in (arguments.old, arguments.new)
    )
    print(f"new/old: {new.get('label')} / {old.get('label')}")
    for case, modes in new["workloads"].items():
        for mode, measurement in modes.items():
            before = old["workloads"].get(case, {}).get(mode)
            if before is None:
                continue
            print(f"{case} {mode}")
            rows = [("total (median)", before["total"], measurement["total"])]
            rows += [
                (description, before["stages"].get(key), measurement["stages"].get(key))
                for key, description, _ in STAGES
            ]
            if "allocations" in before and "allocations" in measurement:
                rows += [
                    (
                        f"allocations ({unit})",
                        before["allocations"][unit],
                        measurement["allocations"][unit],
                    )
                    for unit in ("blocks", "bytes")
                ]
            for label, a, b in rows:
                if a is None or b is None:
                    print(f"  {label:40} {'missing':>17}")
                else:
                    ratio = b / a if a else float("nan")
                    print(f"  {label:40} {a:>17,} -> {b:>17,}  {ratio:.4f}")
            spans = [
                f"[{min(m['runs']):,}, {max(m['runs']):,}]"
                for m in (before, measurement)
                if "runs" in m
            ]
            if len(spans) == 2:
                print(f"  {'total range over the runs':40} {spans[0]} -> {spans[1]}")
            if before["result"] != measurement["result"]:
                print(
                    f"  RESULT differs: {before['result']} -> {measurement['result']}"
                )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("measure")
    command.add_argument("cases", nargs="*")
    command.add_argument("--out", default=str(ROOT / "target" / "generation_ir"))
    command.add_argument("--jobs", type=int, default=16)
    command.add_argument("--runs", type=int, default=3)
    command.add_argument("--modes", default="generation,full")
    command.add_argument("--allocations", action="store_true")
    command.add_argument("--label", default="")
    command.set_defaults(handler=measure)
    command = commands.add_parser("summarize")
    command.add_argument("files", nargs="+")
    command.set_defaults(handler=summarize)
    command = commands.add_parser("compare")
    command.add_argument("old")
    command.add_argument("new")
    command.set_defaults(handler=compare)
    arguments = parser.parse_args()
    arguments.handler(arguments)


if __name__ == "__main__":
    main()
