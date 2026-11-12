#!/usr/bin/env python3
"""Measure the benchmark set of docs/benchmark.typ (resumable, one fresh process per sample).

usage: run_benchmarks.py --binary BIN --data-root ROOT --out-dir OUT [--cases cases.json]
                         [--only ID,...] [--skip-heavy | --heavy-only] [--repetitions N]
                         [--heavy-repetitions N] [--cpus LIST] [--form PATH] [--assemble-only]
                         [--heavy-n-cores N] [--heavy-cpus LIST] [--draft] [--allow-shared-l3]

Layout of OUT: identity.json, inputs/<case>.json (prepared inputs; warmup-<L>L.json: the warm-up
integrals), raw/<case>/<cell>-r<k>.json (one measurement each; a failed sample keeps its error),
raw/<case>/stale/ (samples of another binary, input, data root or FORM), logs/, results.json.
Resumable: a sample is kept only while its key (binary SHA-256, input and warm-up SHA-256, data-root
manifest SHA-256, FORM SHA-256, n_cores, protocol) matches the current run; any other sample is moved
to stale/ and measured again, and the assembly ignores it.

Every sample runs `taskset -c CPU BIN measure --warmup inputs/warmup-<L>L.json ...` with one worker
thread (RAYON/OMP/MKL = 1, vakint n_cores 1) on a CPU of --cpus (default: the CPUs this process may
use; one CPU per L3 domain unless --allow-shared-l3). Samples of one CPU never overlap; several CPUs
run samples concurrently; multi-core heavy samples (--heavy-n-cores > 1) run alone afterwards on
--heavy-cpus (default: all CPUs this process may use). Contention evidence per sample: busy fraction
of the SMT siblings, run-queue delay and foreign CPU time on the timed CPUs (/proc/schedstat minus the
process tree's own CPU time from wait4), run-queue delay of the harness threads around the timed call.
A contended sample (busy SMT sibling, foreign CPU time on the timed CPU, run-queue delay of the timed call;
thresholds --max-sibling-busy, --max-foreign-cpu, --max-run-delay) is measured again up to --contention-retries
times; the last attempt is kept; the assembly flags samples with the thresholds of the assembling run.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import math
import os
import queue
import shutil
import subprocess
import sys
import threading
import time
from decimal import Decimal, InvalidOperation, getcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
CRATE = HERE.parent.parent.parent  # crates/vakint
getcontext().prec = 80

BACKENDS = {  # vakint applicability (lib.rs EvaluationMethod::supports)
    "matad": (1, 3, "MATAD covers one to three loops"),
    "alphaloop": (1, 3, "alphaLoop covers one to three loops"),
    "fmft": (4, 4, "FMFT covers four loops only"),
}
CELLS_RUSTRED = [("rustred", "reduction", "cold"), ("rustred", "reduction", "warm"), ("rustred", "full", "default")]
# Version of the measurement protocol; samples of another protocol are stale.
# 2: untimed warm-up (twice, another small integral of the same loop count) before every timed call.
PROTOCOL = 2


def now():
    return dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256(path: Path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=1, ensure_ascii=False) + "\n")
    tmp.replace(path)


def expand_cpus(text):
    out = []
    for part in text.replace(" ", "").split(","):
        if not part:
            continue
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b or a) + 1))
    return out


def sibling(cpu):
    try:
        text = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list").read_text().strip()
        return [c for c in expand_cpus(text) if c != cpu]
    except OSError:
        return []


def l3_domain(cpu):
    try:
        return Path(f"/sys/devices/system/cpu/cpu{cpu}/cache/index3/shared_cpu_list").read_text().strip()
    except OSError:
        return None


def schedstat(cpus):
    """Per CPU (rq_cpu_time ns, run_delay ns) of /proc/schedstat (sched_info: tasks other than idle)."""
    out = {}
    try:
        for line in Path("/proc/schedstat").read_text().splitlines():
            if line.startswith("cpu"):
                f = line.split()
                c = int(f[0][3:])
                if c in cpus and len(f) >= 9:
                    out[c] = (int(f[7]), int(f[8]))
    except OSError:
        pass
    return out


def proc_cpu_seconds(pid):
    """CPU seconds of process PID with its threads and reaped children (/proc/PID/stat utime + stime + cutime + cstime)."""
    try:
        f = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
        return sum(int(x) for x in f[11:15]) / os.sysconf("SC_CLK_TCK")
    except (OSError, IndexError, ValueError):
        return None


class CallWindow:
    """Multi-core samples: polls the timed CPUs' run time (/proc/schedstat), the measured process's CPU time and this
    poller's own CPU time every 0.25 s, plus whether the output file exists (the harness writes it right after the timed
    call, which is the last call of the process). foreign() = other tasks' share of the timed CPUs during the window
    [output written - call_seconds, output written], i.e. during the timed call, not during the single-threaded data load."""

    def __init__(self, pid, cpus, done_path):
        self.pid, self.cpus, self.done_path, self.rows = pid, set(cpus), done_path, []
        self.stop = threading.Event()
        self.thread = threading.Thread(target=self.loop, daemon=True)
        self.thread.start()

    def loop(self):
        while not self.stop.is_set():
            q = schedstat(self.cpus)
            own = proc_cpu_seconds(self.pid)
            self.rows.append((time.time(), sum(v[0] for v in q.values()) * 1e-9, own, time.thread_time(),
                              self.done_path.exists()))
            self.stop.wait(0.25)

    def finish(self):
        self.stop.set()
        self.thread.join()

    def foreign(self, call_seconds):
        rows = [r for r in self.rows if r[2] is not None]
        end = next((r for r in rows if r[4]), None)
        if end is None or not call_seconds:
            return None
        before = [r for r in rows if r[0] <= end[0] - call_seconds]
        if not before:
            return None
        start = before[-1]
        dt = end[0] - start[0]
        if dt <= 0:
            return None
        other = (end[1] - start[1]) - (end[2] - start[2]) - (end[3] - start[3])
        return {"foreign_cpu_fraction": round(max(0.0, other) / dt, 5), "window_seconds": round(dt, 3),
                "polls": len(rows)}


def cpu_times(cpus):
    out = {}
    for line in Path("/proc/stat").read_text().splitlines():
        if line.startswith("cpu") and line[3:4].isdigit():
            name, *vals = line.split()
            c = int(name[3:])
            if c in cpus:
                vals = list(map(int, vals))
                idle = vals[3] + (vals[4] if len(vals) > 4 else 0)
                out[c] = (sum(vals), idle)
    return out


def cells_for(case):
    cells = list(CELLS_RUSTRED)
    for backend, (lo, hi, _) in BACKENDS.items():
        if lo <= case["loops"] <= hi:
            cells += [(backend, "reduction", "-"), (backend, "full", "-")]
    return cells


def cell_name(cell):
    backend, scope, cache = cell
    return f"{backend}.{scope}" + (f".{cache}" if cache != "-" else "")


# ------------------------------------------------------------- preparation --

def prepare(args, case, out):
    target = out / "inputs" / f"{case['id']}.json"
    spec = json.dumps(case["source"], sort_keys=True)
    if target.exists():
        rec = json.loads(target.read_text())
        if rec.get("case_source") == spec:
            return target
    src = case["source"]
    log = out / "logs" / f"prepare-{case['id']}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    tmp = out / "inputs" / f"{case['id']}.prep.json"
    started = time.time()
    if src["kind"] == "rung":
        work = out / "prep" / f"L{case['loops']}-{src['projection']}"
        shipped = args.data_root / f"{case['loops']}L" / "benchmarks" / src["projection"]
        if (shipped / "prepared-input.json").exists() and not work.exists():
            # the data root ships this prepared physical input (authenticated again by the fixture loader)
            shutil.copytree(shipped, work)
        cmd = [args.binary, "prepare-rung", "--loops", str(case["loops"]), "--projection", src["projection"],
               "--work-dir", str(work), "--out", str(tmp)]
    elif src["kind"] == "graph":
        master = src.get("master")
        manifest = json.loads((args.data_root / f"{case['loops']}L" / "manifest.json").read_text())
        if "master_powers" in src:  # robust to the renaming of a new basis: find the master by its powers
            found = [m["name"] for m in manifest["masters"] if m["powers"] == src["master_powers"]]
            if len(found) != 1:
                print(f"  prepare {case['id']}: no unique master with powers {src['master_powers']}", flush=True)
                return None
            master = found[0]
        elif src.get("select") == "fewest-lines":  # smallest master not used by a case (warm-up)
            used = {tuple(c["source"].get("master_powers", ())) for c in args.all_cases}
            pool = [m for m in manifest["masters"] if tuple(m["powers"]) not in used]
            pool.sort(key=lambda m: (sum(1 for p in m["powers"] if p > 0), sum(abs(p) for p in m["powers"]), m["name"]))
            if not pool:
                print(f"  prepare {case['id']}: no master to select", flush=True)
                return None
            master = pool[0]["name"]
        cmd = [args.binary, "prepare-graph", "--data-root", str(args.data_root), "--loops", str(case["loops"]),
               "--master", master, "--out", str(tmp)]
        if "powers" in src:
            cmd += ["--powers", ",".join(map(str, src["powers"]))]
    elif src["kind"] in ("expression", "fixture"):
        if src["kind"] == "fixture":
            rows = json.loads((CRATE / src["file"]).read_text())
            text = next(r for r in rows if r["case"] == src["case"])["input"]
        else:
            text = src["input"]
        rec = {"schema": "vakint.benchmark-input.v1", "kind": src["kind"], "loop_count": case["loops"],
               "input": text, "input_sha256": hashlib.sha256(text.encode()).hexdigest()}
        write_json(tmp, rec)
        cmd = None
    else:
        raise SystemExit(f"unknown source kind {src['kind']}")
    if cmd:
        env = dict(os.environ, RAYON_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        with open(log, "w") as f:
            rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env).returncode
        if rc != 0:
            print(f"  prepare {case['id']} FAILED rc={rc}; see {log}", flush=True)
            return None
    rec = json.loads(tmp.read_text())
    rec["case_source"] = spec
    if case.get("label"):
        rec["label"] = case["label"]
    rec["prepare_wall_seconds"] = time.time() - started
    write_json(target, rec)
    tmp.unlink(missing_ok=True)
    print(f"  prepared {case['id']} in {rec['prepare_wall_seconds']:.1f} s", flush=True)
    return target


# ------------------------------------------------------------- measurement --

def sample_n_cores(args, case, backend):
    return args.heavy_n_cores if (case.get("heavy") and backend == "rustred") else 1


def sample_key(args, case, cell):
    """Identity of a valid sample: a sample with another key is stale (other binary, input, data, FORM)."""
    backend = cell[0]
    return {"protocol": PROTOCOL, "executable_sha256": args.binary_sha256,
            "input_sha256": args.input_sha.get(case["id"]),
            "warmup_input_sha256": args.input_sha.get(f"warmup-{case['loops']}L"),
            "root_manifest_sha256": args.root_manifest_sha256,
            "form_sha256": args.form_sha256 if backend != "rustred" else None,
            "n_cores": sample_n_cores(args, case, backend)}


def sample_is_current(args, case, cell, rec):
    return (rec.get("harness") or {}).get("key") == sample_key(args, case, cell)


# Below this foreign CPU fraction no other task ran on the timed CPU: a run-queue wait of the timed call is then the
# process tree's own (vakint waiting while its FORM child runs on the same pinned CPU), not contention.
SELF_WAIT_FOREIGN = 0.002


def contention_score(args, meta):
    """Reasons a sample counts as contended (empty: clean). The CPU-level run-queue delay is recorded but is no criterion:
    it includes a backend's own processes waiting for each other (vakint and its FORM child, ~3 % on an idle CPU); the
    timed call's run-queue wait counts only when other tasks ran on the timed CPU (SELF_WAIT_FOREIGN).
    Multi-core samples: the foreign CPU time and the threads' run-queue wait are per timed CPU (the recorded fractions
    are sums over the timed CPUs / threads, divided here by their number), so that the thresholds mean the same share
    of a CPU for one core and for eight; the foreign CPU time is the one during the timed call where it was recorded
    (call_window), else over the whole process."""
    why = []
    n = max(1, len(meta.get("timed_cpus") or ()))
    busy = [v for v in (meta.get("sibling_busy_fraction") or {}).values() if v is not None]
    if busy and max(busy) > args.max_sibling_busy:
        why.append(f"SMT sibling {max(busy):.0%} busy")
    foreign = meta.get("foreign_cpu_fraction")
    window = meta.get("call_window") or {}
    during = n > 1 and window.get("foreign_cpu_fraction") is not None
    if during:  # multi-core: during the timed call (the process's single-threaded phases leave CPUs to other tasks)
        foreign = window["foreign_cpu_fraction"]
    foreign = None if foreign is None else foreign / n
    if (foreign or 0) > args.max_foreign_cpu:
        why.append(f"other tasks ran {foreign:.1%} of the wall time on the timed CPU" if n == 1 else
                   f"other tasks used {foreign:.1%} of each of the {n} timed CPUs"
                   + (" during the timed call" if during else ""))
    call = meta.get("call_run_delay_fraction")
    call = None if call is None else call / n
    if call is not None and call > args.max_run_delay and (foreign is None or foreign >= SELF_WAIT_FOREIGN):
        why.append(f"timed call waited {call:.1%} in the run queue" + ("" if n == 1 else " per thread"))
    return why


def measure(args, case, cell, rep, cpu, out, input_file):
    """One sample; measured again (up to --contention-retries times, multi-core heavy samples up to
    --heavy-contention-retries times) while it was contended, keeping the earlier attempts' times and reasons."""
    attempts = []
    retries = args.heavy_contention_retries if sample_n_cores(args, case, cell[0]) > 1 else args.contention_retries
    for attempt in range(retries + 1):
        meta = measure_once(args, case, cell, rep, cpu, out, input_file, attempts)
        if meta is None or meta["exit_code"] != 0:
            return
        why = meta["contention"]
        if not why:
            return
        if attempt < retries:
            target = out / "raw" / case["id"] / f"{cell_name(cell)}-r{rep}.json"
            attempts.append({"call_seconds": json.loads(target.read_text()).get("call_seconds"), "contention": why})
            target.unlink()
            print(f"  {case['id']:<14} {cell_name(cell):<26} r{rep}: {'; '.join(why)}, re-measuring", flush=True)


def measure_once(args, case, cell, rep, cpu, out, input_file, attempts):
    backend, scope, cache = cell
    name = cell_name(cell)
    target = out / "raw" / case["id"] / f"{name}-r{rep}.json"
    if target.exists():
        try:
            rec = json.loads(target.read_text())
        except (OSError, json.JSONDecodeError):
            rec = {}
        if sample_is_current(args, case, cell, rec):
            return None
        stale = target.parent / "stale"
        stale.mkdir(exist_ok=True)
        target.replace(stale / f"{name}-r{rep}.{int(time.time())}.json")
        print(f"  {case['id']:<14} {name:<26} r{rep}: stale sample (other binary/input/data/protocol) moved to stale/", flush=True)
    log = out / "logs" / case["id"] / f"{name}-r{rep}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    work = out / "raw" / case["id"] / f"{name}-r{rep}.part.json"
    n_cores = sample_n_cores(args, case, backend)
    timed = [cpu] if n_cores == 1 else expand_cpus(args.heavy_cpus)
    cmd = ["taskset", "-c", ",".join(map(str, timed)),
           args.binary, "measure", "--data-root", str(args.data_root), "--input", str(input_file),
           "--backend", backend, "--scope", scope, "--n-cores", str(n_cores), "--out", str(work)]
    warm = out / "inputs" / f"warmup-{case['loops']}L.json"
    if warm.exists():
        cmd += ["--warmup", str(warm)]
    if cache != "-":
        cmd += ["--cache", cache]
    if backend != "rustred":
        cmd += ["--form", args.form]
    env = dict(os.environ, RAYON_NUM_THREADS=str(n_cores), OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               SYMBOLICA_HIDE_BANNER="1")
    sib = sorted({s for c in timed for s in sibling(c)} - set(timed))
    load0 = os.getloadavg()[0]
    log_file = open(log, "w")
    t0 = cpu_times(timed + sib)
    time.sleep(0.002)  # block once: this thread's CPU time so far is accounted before the snapshot
    own0 = time.thread_time()
    q0 = schedstat(set(timed))
    started, start_utc = time.time(), now()
    proc = subprocess.Popen(cmd, stdout=log_file, stderr=subprocess.STDOUT, env=env)
    window = CallWindow(proc.pid, timed, work) if n_cores > 1 else None
    _, status, ru = os.wait4(proc.pid, 0)
    if window:
        window.finish()
    proc.returncode = os.waitstatus_to_exitcode(status)
    wall = time.time() - started
    q1 = schedstat(set(timed))
    harness_cpu = time.thread_time() - own0  # this worker thread (pinned to its CPU) inside the window
    t1 = cpu_times(timed + sib)
    log_file.close()
    returncode = proc.returncode

    def busy(c):
        if c not in t0 or c not in t1:
            return None
        total = t1[c][0] - t0[c][0]
        idle = t1[c][1] - t0[c][1]
        return round(1 - idle / total, 4) if total > 0 else None

    own_cpu = ru.ru_utime + ru.ru_stime  # the process tree (wait4: the child and its reaped descendants)
    rq_cpu = sum(q1[c][0] - q0[c][0] for c in timed if c in q0 and c in q1) * 1e-9 - harness_cpu
    rq_delay = sum(q1[c][1] - q0[c][1] for c in timed if c in q0 and c in q1) * 1e-9
    meta = {"case": case["id"], "cell": name, "repetition": rep, "cpu": cpu, "timed_cpus": timed, "siblings": sib,
            "l3_domain": l3_domain(cpu),
            "sibling_busy_fraction": {str(c): busy(c) for c in sib}, "cpu_busy_fraction": busy(cpu),
            "process_cpu_seconds": round(own_cpu, 4),
            "cpu_run_delay_fraction": round(rq_delay / wall, 5) if wall > 0 and q0 else None,
            "foreign_cpu_fraction": round(max(0.0, rq_cpu - own_cpu) / wall, 5) if wall > 0 and q0 else None,
            "load1_start": load0, "load1_end": os.getloadavg()[0], "started_utc": start_utc,
            "wall_seconds": wall, "exit_code": returncode, "command": cmd, "n_cores": n_cores,
            "key": sample_key(args, case, cell), "contended_attempts": list(attempts)}
    if returncode == 0 and work.exists():
        rec = json.loads(work.read_text())
        meta["call_run_delay_fraction"] = (rec.get("call_usage") or {}).get("run_delay_fraction")
        if window:  # other tasks on the timed CPUs during the timed call (sum over the CPUs, as foreign_cpu_fraction)
            meta["call_window"] = window.foreign(rec.get("call_seconds"))
        meta["contention"] = contention_score(args, meta)
        meta["contended"] = bool(meta["contention"])
        rec["harness"] = meta
        write_json(target, rec)
        work.unlink()
        flag = f"  [contended: {'; '.join(meta['contention'])}]" if meta["contended"] else ""
        print(f"  {case['id']:<14} {name:<26} r{rep} cpu{cpu}: {rec['call_seconds']:.4g} s{flag}", flush=True)
    else:
        meta["contention"] = []
        tail = log.read_text(errors="replace")[-4000:] if log.exists() else ""
        write_json(target, {"schema": "vakint.benchmark-measurement.v1", "failed": True, "error_tail": tail,
                            "harness": meta, "backend": backend, "scope": scope, "cache": cache})
        work.unlink(missing_ok=True)
        print(f"  {case['id']:<14} {name:<26} r{rep} cpu{cpu}: FAILED rc={returncode}", flush=True)
    shutil.rmtree(work.with_suffix(".work"), ignore_errors=True)  # FORM scratch of the harness
    return meta


# ---------------------------------------------------------------- assembly --

def dec(text):
    try:
        return Decimal(text)
    except (InvalidOperation, TypeError):
        return None


def stats(values):
    v = sorted(values)
    n = len(v)
    med = (v[(n - 1) // 2] + v[n // 2]) / 2
    return {"n": n, "min": v[0], "median": med, "max": v[-1], "samples": values}


args_precision = 32  # measure --precision default


def agreement(records):
    """Digits of agreement of every 'full' result with RustRed's, per epsilon order."""
    def coeffs(rec):
        num = rec.get("numerical") or {}
        out = {}
        for c in num.get("coefficients", []):
            digits = [c.get(k) for k in ("real_significant_digits", "imaginary_significant_digits")
                      if c.get(k) is not None]
            out[c["epsilon_power"]] = (dec(c["real"]), dec(c["imaginary"]), min(digits + [99]))
        return out
    ref = records.get("rustred")
    if not ref:
        return {"status": "n/a", "reason": "no RustRed numerical result"}
    rc = coeffs(ref)
    result = {"reference": "rustred", "compared": {}}
    worst = None
    for backend, rec in records.items():
        if backend == "rustred":
            continue
        oc = coeffs(rec)
        if not oc:
            result["compared"][backend] = {"status": "no numerical result",
                                           "error": (rec.get("numerical") or {}).get("error")}
            continue
        per = {}
        ok = True
        for p in sorted(set(rc) | set(oc)):
            a, b = rc.get(p), oc.get(p)
            if a is None or b is None:
                per[str(p)] = {"digits": None, "note": "order missing in one result"}
                continue
            num = abs(a[0] - b[0]) + abs(a[1] - b[1])
            den = max(abs(a[0]) + abs(a[1]), abs(b[0]) + abs(b[1]))
            # digits the less precise result carries (requested precision when unreported)
            expected = float(min(a[2], b[2], args_precision))
            if den == 0:
                per[str(p)] = {"digits": None, "note": "both exactly zero"}
                continue
            digits = expected if num == 0 else min(float(-(num / den).log10()), expected)
            per[str(p)] = {"digits": round(digits, 2), "expected": expected}
            if digits < max(expected - 4, 6):
                ok = False
        ds = [v["digits"] for v in per.values() if v.get("digits") is not None]
        md = min(ds) if ds else None
        result["compared"][backend] = {"status": "agree" if ok and ds else "DISAGREE", "min_digits": md,
                                       "orders": per}
        if md is not None:
            worst = md if worst is None else min(worst, md)
    statuses = [v.get("status") for v in result["compared"].values()]
    result["status"] = ("n/a" if not statuses else "DISAGREE" if "DISAGREE" in statuses else
                        "agree" if all(s == "agree" for s in statuses) else "partial")
    result["min_digits"] = worst
    return result


def summarize_error(tail):
    """Short, attributed form of a failed call's error (who rejected what), for the report."""
    import re
    lines = [l for l in (tail or "").splitlines() if l.strip()]
    last = lines[-1] if lines else "no output"
    m = re.search(r"NumericalPrecision\(\"coefficient multiplying (.*?): numerical precision: (real|imaginary) component has no "
                  r"reliable significant decimal digit; midpoint ([-+0-9.eE]+), absolute uncertainty ([-+0-9.eE]+)", last)
    if m:
        coeff, comp, mid, err = m.groups()
        g = lambda x: f"{float(x):.1e}".replace("-", "−")
        zero = abs(float(mid)) <= float(err)
        return {"source": "vakint precision guard",
                "text": f"vakint's precision guard rejects the backend's result: the {comp} part of the coefficient of "
                        f"{coeff} is {g(mid)} ± {g(err)}" + (" (consistent with zero, no reliable digit)" if zero else
                                                              " (no reliable digit)")}
    if "MissingMasterOrder {" in last:
        found = re.findall(r'MissingMasterOrder \{ master: "L(\d+)_[^"]*", first_missing_order: (-?\d+), last_needed_order: '
                           r'(?:Some\((-?\d+)\)|None), first_affected_result_order: (-?\d+)', last)
        if found:
            by_loop = {}
            for loops, first, last_needed, _ in found:
                lo, hi = by_loop.get(loops, (None, None))
                need = int(last_needed or first)
                by_loop[loops] = (int(first) if lo is None else min(lo, int(first)), need if hi is None else max(hi, need))
            affected = min(int(a) for *_, a in found)
            need = ", ".join(f"{L}L to ε^{hi}" for L, (_, hi) in sorted(by_loop.items(), reverse=True))
            return {"source": "vakint (fails closed)",
                    "text": f"vakint fails closed (MissingMasterOrders): {len(found)} masters lack Laurent orders this result needs "
                            f"(the reduction coefficients have poles at d = 4; deepest orders needed: {need}), so the result fails "
                            f"from ε^{affected} on".replace("^-", "^−")}
    return {"source": "backend", "text": last[:240]}


def source_equivalent(repo, build_commit):
    """Receipt rule 1 on the tools checkout: the oldest commit of HEAD's first-parent line whose source outside the neutral
    paths equals the measured build's (a reader can resolve it when the build commit itself was rebased away)."""
    if not build_commit:
        return None
    spec = importlib.util.spec_from_file_location("bench_receipt", HERE / "receipt.py")
    rc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rc)
    want = rc.source_digest(repo, build_commit)
    if want is None:
        return None
    found = None
    for c in (rc.git(repo, "rev-list", "--first-parent", "-300", "HEAD") or "").split():
        if rc.source_digest(repo, c) != want:
            break
        found = c
    return {"commit": found, "source_digest": want,
            "build_in_history": rc.git(repo, "merge-base", "--is-ancestor", build_commit, "HEAD") is not None}


def assemble(args, cases_doc, out, identity, pipeline):
    """OUT/annotations.json (optional): {"cells": {"<case> <cell>": "<why the cell was not timed in an idle window>"}},
    e.g. a long single run that overlapped other work; the PDF marks such cells like contended ones."""
    ann_path = out / "annotations.json"
    annotations = json.loads(ann_path.read_text()).get("cells", {}) if ann_path.exists() else {}
    cases = []
    stale_total = 0
    one_time = {}  # loops -> backend -> lists
    for case in cases_doc["cases"]:
        row = {k: case.get(k) for k in ("id", "loops", "title", "subtitle", "heavy")}
        row["kind"] = case["source"]["kind"]
        inp = out / "inputs" / f"{case['id']}.json"
        if inp.exists():
            rec = json.loads(inp.read_text())
            row["input"] = {k: rec.get(k) for k in ("input_sha256", "fixture_scalar_terms", "preparation_seconds",
                                                    "prepare_wall_seconds", "projection", "master", "powers")}
            row["input"]["bytes"] = len(rec.get("input", ""))
        cells = {}
        full_records = {}
        for cell in [c for c in cells_for(case)]:
            name = cell_name(cell)
            samples = sorted((out / "raw" / case["id"]).glob(f"{name}-r*.json")) if (out / "raw" / case["id"]).exists() else []
            recs = [json.loads(p.read_text()) for p in samples]
            current = [r for r in recs if sample_is_current(args, case, cell, r)]
            stale_total += len(recs) - len(current)
            good = [r for r in current if not r.get("failed")]
            bad = [r for r in current if r.get("failed")]
            if good:
                c = {"status": "ok", **stats([r["call_seconds"] for r in good])}
                c["loader_median"] = stats([r["loader_seconds"] for r in good])["median"]
                setups = [r["warmup"]["one_time_setup_seconds"] for r in good if isinstance(r.get("warmup"), dict)]
                if setups:
                    c["one_time_setup_median"] = stats(setups)["median"]
                    c["warmup_first_median"] = stats([r["warmup"]["first_seconds"] for r in good])["median"]
                    slot = one_time.setdefault(str(case["loops"]), {}).setdefault(name, {"setup": [], "loader": []})
                    slot["setup"] += setups
                    if cell[0] == "rustred":
                        slot["loader"] += [r["loader_seconds"] for r in good]
                why = [contention_score(args, r["harness"]) for r in good]
                c["contended"] = sum(1 for w in why if w)
                # the reported median rests on a contended sample (odd n: the middle one; even n: either middle one)
                vals = sorted((r["call_seconds"], bool(w)) for r, w in zip(good, why))
                mids = vals[len(vals) // 2:len(vals) // 2 + 1] if len(vals) % 2 else vals[len(vals) // 2 - 1:len(vals) // 2 + 1]
                c["median_contended"] = any(w for _, w in mids)
                c["contention"] = sorted({x for r, w in zip(good, why) for x in w
                                          if (r["call_seconds"], True) in mids or not c["median_contended"]})[:6]
                if f"{case['id']} {name}" in annotations:
                    c["overlap"] = annotations[f"{case['id']} {name}"]
                math_s = [r["reduction_capture"]["reduction_math_seconds"] for r in good
                          if isinstance(r.get("reduction_capture"), dict) and r["reduction_capture"].get("reduction_math_seconds")]
                if math_s:
                    c["reduction_math_median"] = stats(math_s)["median"]
                c["peak_rss_bytes"] = max((r.get("peak_rss_bytes") or 0) for r in good)
                c["masters_retained"] = good[0].get("retained_rustred_masters")
                c["output_sha256"] = sorted({r["output_sha256"] for r in good})
                c["n_cores"] = max(r.get("n_cores") or 1 for r in good)
                c["cpus"] = sorted({x for r in good for x in r["harness"].get("timed_cpus", [r["harness"]["cpu"]])})
                c["l3_domains_per_sample"] = max(len({l3_domain(x) for x in r["harness"].get("timed_cpus", [r["harness"]["cpu"]])})
                                                 for r in good)
                c["sibling_busy_max"] = max((v for r in good for v in r["harness"]["sibling_busy_fraction"].values()
                                             if v is not None), default=None)
                c["load1_max"] = max(r["harness"]["load1_start"] for r in good)
                if len(c["output_sha256"]) > 1:
                    c["warning"] = "repetitions returned different expressions"
                if cell[1] == "full":
                    full_records[cell[0]] = good[0]
            elif bad:
                tail = bad[-1].get("error_tail", "")
                h = bad[-1].get("harness") or {}
                full_log = out / "logs" / case["id"] / f"{name}-r{h.get('repetition')}.log"
                if full_log.exists():  # the whole error (a MissingMasterOrders list is longer than the stored tail)
                    tail = full_log.read_text(errors="replace")
                lines = [l for l in tail.splitlines() if l.strip()]
                c = {"status": "failed", "error": lines[-1][:300] if lines else "no output", "n": len(bad),
                     "error_summary": summarize_error(tail)}
            elif case.get("heavy") and args.skip_heavy:
                c = {"status": "pending", "reason": "heavy cell, measured in its own run"}
            else:
                c = {"status": "pending"}
            cells[name] = c
        for backend, (lo, hi, reason) in BACKENDS.items():
            if not lo <= case["loops"] <= hi:
                cells[f"{backend}.reduction"] = {"status": "n/a", "reason": reason}
                cells[f"{backend}.full"] = {"status": "n/a", "reason": reason}
        row["cells"] = cells
        row["agreement"] = agreement(full_records)
        if "rustred" in full_records:
            row["rustred_numerical"] = full_records["rustred"].get("numerical")
        cases.append(row)
    sib, loads, used, delays, foreign, n_samples, n_contended, l3, calls = [], [], set(), [], [], 0, 0, set(), []
    by_case = {c["id"]: c for c in cases_doc["cases"]}
    for raw in (out / "raw").glob("*/*.json"):
        try:
            rec = json.loads(raw.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        h = rec.get("harness") or {}
        case = by_case.get(h.get("case"))
        cell = next((c for c in cells_for(case) if cell_name(c) == h.get("cell")), None) if case else None
        if not cell or not sample_is_current(args, case, cell, rec) or rec.get("failed"):
            continue
        n_samples += 1
        n_contended += bool(contention_score(args, h))
        used |= set(h.get("timed_cpus") or [h["cpu"]])
        if h.get("l3_domain"):
            l3.add(h["l3_domain"])
        sib += [v for v in (h.get("sibling_busy_fraction") or {}).values() if v is not None]
        n_t = max(1, len(h.get("timed_cpus") or ()))  # per timed CPU (contention_score)
        if h.get("cpu_run_delay_fraction") is not None:
            delays.append(h["cpu_run_delay_fraction"] / n_t)
        if h.get("foreign_cpu_fraction") is not None:
            foreign.append(h["foreign_cpu_fraction"] / n_t)
        if h.get("call_run_delay_fraction") is not None:
            calls.append(h["call_run_delay_fraction"] / n_t)
        if h.get("load1_start") is not None:
            loads.append(h["load1_start"])

    def q(v, f):
        v = sorted(v)
        return round(v[min(len(v) - 1, int(f * len(v)))], 4) if v else None
    contention = {"samples": n_samples, "contended_samples": n_contended,
                  "sibling_busy_median": q(sib, 0.5), "sibling_busy_max": q(sib, 1.0),
                  "cpu_run_delay_median": q(delays, 0.5), "cpu_run_delay_max": q(delays, 1.0),
                  "call_run_delay_median": q(calls, 0.5), "call_run_delay_max": q(calls, 1.0),
                  "foreign_cpu_median": q(foreign, 0.5), "foreign_cpu_max": q(foreign, 1.0),
                  "l3_domains": len(l3),
                  "thresholds": {"sibling_busy": args.max_sibling_busy, "run_delay": args.max_run_delay,
                                 "foreign_cpu": args.max_foreign_cpu}}
    if loads:
        contention["load1_max"] = round(max(loads), 1)
    if used:
        identity.setdefault("host", {})["cpus_used"] = sorted(used)
    # One-time cost per process (excluded from every cell): RustRed data load, first-call setup of every backend.
    once = {L: {b: {"setup_median": stats(v["setup"])["median"] if v["setup"] else None,
                    "loader_median": stats(v["loader"])["median"] if v["loader"] else None}
                for b, v in per.items()} for L, per in one_time.items()}
    ok_cells = [(c, k, v) for c in cases for k, v in c["cells"].items() if v["status"] == "ok"]
    timing = {"multi_core_cells": [f"{c['id']} {k}" for c, k, v in ok_cells if (v.get("n_cores") or 1) > 1],
              "single_run_cells": [f"{c['id']} {k}" for c, k, v in ok_cells if v["n"] == 1],
              "min_repetitions": min((v["n"] for _, _, v in ok_cells), default=None),
              "max_repetitions": max((v["n"] for _, _, v in ok_cells), default=None)}
    timing["all_single_core"] = not timing["multi_core_cells"]
    result = {
        "schema": "vakint.benchmark-results.v1",
        "contention": contention,
        "one_time_per_process": once,
        "timing": timing,
        "stale_samples_ignored": stale_total,
        "protocol": PROTOCOL,
        "draft": bool(args.draft),
        "assembled_utc": now(),
        "identity": identity,
        "pipeline": pipeline,
        "groups": cases_doc["groups"],
        "cases": cases,
    }
    write_json(out / "results.json", result)
    return result


def host_info(cpus, form):
    cpuinfo = Path("/proc/cpuinfo").read_text()
    model = next((l.split(":", 1)[1].strip() for l in cpuinfo.splitlines() if l.startswith("model name")), None)
    try:
        fv = subprocess.run([form, "-v"], capture_output=True, text=True, timeout=60).stdout.strip().splitlines()[0]
        # drop the build date ("FORM 5.0.0 (Jan  1 1980, v5.0.0) 64-bits": Nix epoch artefact)
        import re
        fv = re.sub(r"\s*\(\s*\)", "", re.sub(r"[A-Z][a-z]{2}\s+\d{1,2}\s+\d{4},?\s*", "", fv))
    except Exception as e:  # FORM is optional for RustRed-only runs
        fv = f"unavailable: {e}"
    return {"cpu_model": model, "logical_cpus": os.cpu_count(), "cpus_used": cpus,
            "siblings": {str(c): sibling(c) for c in cpus}, "load1": os.getloadavg()[0],
            "hostname": os.uname().nodename, "form_version": fv, "form_path": shutil.which(form) or form,
            "form_sha256": sha256(Path(shutil.which(form))) if shutil.which(form) else None}


def git_info(path):
    def g(*a):
        try:
            r = subprocess.run(["git", "-C", str(path), *a], capture_output=True, text=True)
            return r.stdout.strip() if r.returncode == 0 else None
        except Exception:
            return None
    if g("rev-parse", "--is-inside-work-tree") != "true":
        return None
    return {"commit": g("rev-parse", "HEAD"), "branch": g("rev-parse", "--abbrev-ref", "HEAD"),
            "describe": g("describe", "--tags", "--always", "--dirty"),
            "remote": g("config", "--get", "remote.origin.url"),
            "dirty_files": len([l for l in (g("status", "--porcelain") or "").splitlines() if l.strip()])}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--cases", type=Path, default=HERE.parent / "cases.json")
    ap.add_argument("--only")
    ap.add_argument("--skip-heavy", action="store_true")
    ap.add_argument("--heavy-only", action="store_true")
    ap.add_argument("--repetitions", type=int, default=3)
    ap.add_argument("--heavy-repetitions", type=int, default=1)
    ap.add_argument("--heavy-n-cores", type=int, default=1)
    ap.add_argument("--heavy-cpus", help="CPU list for multi-core heavy RustRed samples")
    ap.add_argument("--cpus")
    ap.add_argument("--form", default="form")
    ap.add_argument("--assemble-only", action="store_true")
    ap.add_argument("--prepare-only", action="store_true", help="prepare inputs only (no measurement, no assembly)")
    ap.add_argument("--cells", help="measure only these cells, e.g. rustred.reduction.cold,fmft.full")
    ap.add_argument("--max-sibling-busy", type=float, default=0.25,
                    help="re-measure a sample whose SMT sibling was busier than this fraction (default 0.25)")
    ap.add_argument("--max-run-delay", type=float, default=0.03,
                    help="re-measure a sample whose timed call waited in the run queue longer than this fraction of its time")
    ap.add_argument("--max-foreign-cpu", type=float, default=0.03,
                    help="re-measure a sample when other tasks used more than this fraction of the timed CPU")
    ap.add_argument("--contention-retries", type=int, default=2)
    ap.add_argument("--heavy-contention-retries", type=int, default=0,
                    help="re-measure a contended multi-core sample this many times (default 0: they are long)")
    ap.add_argument("--allow-shared-l3", action="store_true", help="allow concurrent workers that share an L3 cache")
    ap.add_argument("--draft", action="store_true")
    args = ap.parse_args()
    args.binary = str(Path(args.binary).resolve())
    args.data_root = args.data_root.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    cases_doc = json.loads(args.cases.read_text())
    args.all_cases = cases_doc["cases"]
    cpus = expand_cpus(args.cpus) if args.cpus else sorted(os.sched_getaffinity(0))
    if args.heavy_n_cores > 1 and not args.heavy_cpus:
        args.heavy_cpus = ",".join(map(str, sorted(os.sched_getaffinity(0))))
    if args.heavy_n_cores > 1 and len(expand_cpus(args.heavy_cpus)) < args.heavy_n_cores:
        raise SystemExit(f"--heavy-cpus has fewer than --heavy-n-cores = {args.heavy_n_cores} CPUs")
    args.binary_sha256 = sha256(Path(args.binary))
    form_path = shutil.which(args.form)
    args.form_sha256 = sha256(Path(form_path)) if form_path else None
    started = time.time()
    pipeline_file = out / "pipeline.json"
    pipeline = json.loads(pipeline_file.read_text()) if pipeline_file.exists() else {"runs": []}

    ident_path = out / "identity.json"
    if not ident_path.exists() or not args.assemble_only:
        subprocess.run([args.binary, "identity", "--data-root", str(args.data_root), "--out", str(ident_path)], check=True)
    identity = json.loads(ident_path.read_text())
    if identity.get("executable_sha256") != args.binary_sha256 or Path(identity.get("data_root", "")) != args.data_root:
        subprocess.run([args.binary, "identity", "--data-root", str(args.data_root), "--out", str(ident_path)], check=True)
        identity = json.loads(ident_path.read_text())
    args.root_manifest_sha256 = identity.get("root_manifest_sha256")
    identity["host"] = host_info(cpus, args.form)
    # the measured build (embedded at build time), not the checkout these tools run from
    identity["source"] = identity.get("build") or {}
    identity["tools_checkout"] = git_info(CRATE)
    identity["source_equivalent"] = source_equivalent(CRATE, identity["source"].get("commit"))
    identity["data_git"] = git_info(args.data_root)
    identity["cases_sha256"] = sha256(args.cases)
    identity["tools_sha256"] = {p.name: sha256(p) for p in sorted(HERE.glob("*.py"))}

    run = {"started_utc": now(), "argv": sys.argv, "cpus": cpus}
    # input hashes of the prepared inputs (sample keys); warm-up integrals per loop count
    warmups = [{"id": f"warmup-{L}L", "loops": int(L), "label": w.get("label"), "source": w["source"]}
               for L, w in cases_doc.get("warmups", {}).items()]
    args.input_sha = {}
    for p in (out / "inputs").glob("*.json") if (out / "inputs").exists() else []:
        try:
            args.input_sha[p.stem] = json.loads(p.read_text()).get("input_sha256")
        except (OSError, json.JSONDecodeError):
            pass
    if not args.assemble_only:
        selected = cases_doc["cases"]
        if args.only:
            keep = set(args.only.split(","))
            selected = [c for c in selected if c["id"] in keep]
        if args.skip_heavy:
            selected = [c for c in selected if not c.get("heavy")]
        if args.heavy_only:
            selected = [c for c in selected if c.get("heavy")]
        t = time.time()
        inputs = {}
        for case in [w for w in warmups if w["loops"] in {c["loops"] for c in selected}] + selected:
            p = prepare(args, case, out)
            if p:
                inputs[case["id"]] = p
                args.input_sha[case["id"]] = json.loads(p.read_text()).get("input_sha256")
        run["prepare_seconds"] = time.time() - t
        if args.prepare_only:
            print(f"prepared {len(inputs)} inputs in {run['prepare_seconds']:.1f} s", flush=True)
            return
        jobs, heavy_jobs = queue.Queue(), []
        for case in selected:
            if case["id"] not in inputs:
                continue
            reps = args.heavy_repetitions if case.get("heavy") else args.repetitions
            for rep in range(reps):
                cells = cells_for(case)
                if args.cells:
                    wanted = set(args.cells.split(","))
                    cells = [c for c in cells if cell_name(c) in wanted]
                if rep % 2 == 1:
                    cells = list(reversed(cells))
                for cell in cells:
                    if sample_n_cores(args, case, cell[0]) > 1:
                        heavy_jobs.append((case, cell, rep))
                    else:
                        jobs.put((case, cell, rep))
        t = time.time()
        workers = cpus[:max(1, jobs.qsize())]
        domains = {}
        for c in workers:
            domains.setdefault(l3_domain(c), []).append(c)
        shared = [v for v in domains.values() if len(v) > 1]
        if shared and not args.allow_shared_l3:
            raise SystemExit(f"worker CPUs share an L3 cache: {shared}; give one CPU per L3 domain or --allow-shared-l3")

        def worker(cpu):
            try:  # this thread's own work stays on its CPU (never on another worker's timed CPU)
                os.sched_setaffinity(threading.get_native_id(), {cpu})
            except OSError:
                pass
            while True:
                try:
                    case, cell, rep = jobs.get_nowait()
                except queue.Empty:
                    return
                measure(args, case, cell, rep, cpu, out, inputs[case["id"]])

        threads = [threading.Thread(target=worker, args=(c,)) for c in workers]
        for th in threads:
            th.start()
        for th in threads:
            th.join()
        for case, cell, rep in heavy_jobs:  # multi-core samples alone, one after the other
            measure(args, case, cell, rep, expand_cpus(args.heavy_cpus)[0], out, inputs[case["id"]])
        run["measure_seconds"] = time.time() - t
    run["finished_utc"] = now()
    run["wall_seconds"] = time.time() - started
    pipeline["runs"].append(run)
    write_json(pipeline_file, pipeline)
    result = assemble(args, cases_doc, out, identity, pipeline)
    n_ok = sum(1 for c in result["cases"] for v in c["cells"].values() if v["status"] == "ok")
    n_all = sum(1 for c in result["cases"] for v in c["cells"].values() if v["status"] != "n/a")
    print(f"assembled {out / 'results.json'}: {n_ok}/{n_all} cells measured", flush=True)


if __name__ == "__main__":
    main()
