#!/usr/bin/env python3
"""Measure the benchmark set of docs/benchmark.typ (resumable, one fresh process per sample).

usage: run_benchmarks.py --binary BIN --data-root ROOT --out-dir OUT [--cases cases.json]
                         [--only ID,...] [--skip-heavy | --heavy-only] [--repetitions N]
                         [--heavy-repetitions N] [--cpus LIST] [--form PATH] [--assemble-only]
                         [--heavy-n-cores N] [--draft]

Layout of OUT: identity.json, inputs/<case>.json (prepared inputs), raw/<case>/<cell>-r<k>.json
(one measurement each; a failed sample keeps its error), logs/, results.json (assembled).
A sample that already exists is never re-measured, so the command resumes after an interruption.

Every sample runs `taskset -c CPU BIN measure ...` with one worker thread (RAYON/OMP/MKL = 1 and
vakint n_cores 1) on a CPU of --cpus (default: the CPUs this process may use). Samples of one CPU
never overlap; several CPUs run samples concurrently. The record keeps the CPU, the 1-minute load,
and the busy fraction of the CPU's SMT sibling during the sample (contention evidence).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
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
        cmd = [args.binary, "prepare-graph", "--data-root", str(args.data_root), "--loops", str(case["loops"]),
               "--master", src["master"], "--out", str(tmp)]
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
    rec["prepare_wall_seconds"] = time.time() - started
    write_json(target, rec)
    tmp.unlink(missing_ok=True)
    print(f"  prepared {case['id']} in {rec['prepare_wall_seconds']:.1f} s", flush=True)
    return target


# ------------------------------------------------------------- measurement --

def measure(args, case, cell, rep, cpu, out, input_file):
    """One sample; re-measured (up to --contention-retries times) while the SMT sibling of the
    timed CPU was busier than --max-sibling-busy, keeping the earlier attempts in the record."""
    attempts = []
    for attempt in range(args.contention_retries + 1):
        meta = measure_once(args, case, cell, rep, cpu, out, input_file, attempts)
        if meta is None:
            return
        busy = [v for v in meta["sibling_busy_fraction"].values() if v is not None]
        if not busy or max(busy) <= args.max_sibling_busy or meta["exit_code"] != 0:
            return
        if attempt < args.contention_retries:
            target = out / "raw" / case["id"] / f"{cell_name(cell)}-r{rep}.json"
            attempts.append(json.loads(target.read_text()).get("call_seconds"))
            target.unlink()
            print(f"  {case['id']:<14} {cell_name(cell):<26} r{rep}: sibling busy {max(busy):.2f}, re-measuring", flush=True)


def measure_once(args, case, cell, rep, cpu, out, input_file, attempts):
    backend, scope, cache = cell
    name = cell_name(cell)
    target = out / "raw" / case["id"] / f"{name}-r{rep}.json"
    if target.exists():
        return None
    log = out / "logs" / case["id"] / f"{name}-r{rep}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    work = out / "raw" / case["id"] / f"{name}-r{rep}.part.json"
    n_cores = args.heavy_n_cores if (case.get("heavy") and backend == "rustred") else 1
    cmd = ["taskset", "-c", str(cpu) if n_cores == 1 else args.heavy_cpus or str(cpu),
           args.binary, "measure", "--data-root", str(args.data_root), "--input", str(input_file),
           "--backend", backend, "--scope", scope, "--n-cores", str(n_cores), "--out", str(work)]
    if cache != "-":
        cmd += ["--cache", cache]
    if backend != "rustred":
        cmd += ["--form", args.form]
    threads = str(n_cores)
    env = dict(os.environ, RAYON_NUM_THREADS=threads, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               SYMBOLICA_HIDE_BANNER="1")
    sib = sibling(cpu)
    t0 = cpu_times([cpu] + sib)
    load0 = os.getloadavg()[0]
    started, start_utc = time.time(), now()
    with open(log, "w") as f:
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, env=env)
    wall = time.time() - started
    t1 = cpu_times([cpu] + sib)

    def busy(c):
        if c not in t0 or c not in t1:
            return None
        total = t1[c][0] - t0[c][0]
        idle = t1[c][1] - t0[c][1]
        return round(1 - idle / total, 4) if total > 0 else None

    meta = {"case": case["id"], "cell": name, "repetition": rep, "cpu": cpu, "siblings": sib,
            "sibling_busy_fraction": {str(c): busy(c) for c in sib}, "cpu_busy_fraction": busy(cpu),
            "load1_start": load0, "load1_end": os.getloadavg()[0], "started_utc": start_utc,
            "wall_seconds": wall, "exit_code": proc.returncode, "command": cmd, "n_cores": n_cores,
            "contended_attempt_seconds": list(attempts)}
    if proc.returncode == 0 and work.exists():
        rec = json.loads(work.read_text())
        rec["harness"] = meta
        write_json(target, rec)
        work.unlink()
        print(f"  {case['id']:<14} {name:<26} r{rep} cpu{cpu}: {rec['call_seconds']:.4g} s", flush=True)
    else:
        tail = log.read_text(errors="replace")[-4000:] if log.exists() else ""
        write_json(target, {"schema": "vakint.benchmark-measurement.v1", "failed": True, "error_tail": tail,
                            "harness": meta, "backend": backend, "scope": scope, "cache": cache})
        print(f"  {case['id']:<14} {name:<26} r{rep} cpu{cpu}: FAILED rc={proc.returncode}", flush=True)
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


def assemble(args, cases_doc, out, identity, pipeline):
    cases = []
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
            good = [r for r in recs if not r.get("failed")]
            bad = [r for r in recs if r.get("failed")]
            if good:
                c = {"status": "ok", **stats([r["call_seconds"] for r in good])}
                c["loader_median"] = stats([r["loader_seconds"] for r in good])["median"]
                math_s = [r["reduction_capture"]["reduction_math_seconds"] for r in good
                          if isinstance(r.get("reduction_capture"), dict) and r["reduction_capture"].get("reduction_math_seconds")]
                if math_s:
                    c["reduction_math_median"] = stats(math_s)["median"]
                c["peak_rss_bytes"] = max((r.get("peak_rss_bytes") or 0) for r in good)
                c["masters_retained"] = good[0].get("retained_rustred_masters")
                c["output_sha256"] = sorted({r["output_sha256"] for r in good})
                c["n_cores"] = good[0].get("n_cores")
                c["cpus"] = sorted({r["harness"]["cpu"] for r in good})
                c["sibling_busy_max"] = max((v for r in good for v in r["harness"]["sibling_busy_fraction"].values()
                                             if v is not None), default=None)
                c["load1_max"] = max(r["harness"]["load1_start"] for r in good)
                if len(c["output_sha256"]) > 1:
                    c["warning"] = "repetitions returned different expressions"
                if cell[1] == "full":
                    full_records[cell[0]] = good[0]
            elif bad:
                tail = bad[-1].get("error_tail", "")
                lines = [l for l in tail.splitlines() if l.strip()]
                c = {"status": "failed", "error": lines[-1][:300] if lines else "no output"}
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
    sib, loads, used = [], [], set()
    for raw in (out / "raw").glob("*/*.json"):
        try:
            h = json.loads(raw.read_text()).get("harness") or {}
        except (OSError, json.JSONDecodeError):
            continue
        if h.get("cpu") is not None:
            used.add(h["cpu"])
            if h.get("n_cores", 1) > 1 and h.get("command"):
                used |= set(expand_cpus(h["command"][2]))
        sib += [v for v in (h.get("sibling_busy_fraction") or {}).values() if v is not None]
        if h.get("load1_start") is not None:
            loads.append(h["load1_start"])
    contention = {}
    if sib:
        sib.sort()
        contention = {"sibling_busy_median": round(sib[len(sib) // 2], 3), "sibling_busy_max": round(sib[-1], 3)}
    if loads:
        contention["load1_max"] = round(max(loads), 1)
    if used:
        identity.setdefault("host", {})["cpus_used"] = sorted(used)
    result = {
        "schema": "vakint.benchmark-results.v1",
        "contention": contention,
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
    except Exception as e:  # FORM is optional for RustRed-only runs
        fv = f"unavailable: {e}"
    return {"cpu_model": model, "logical_cpus": os.cpu_count(), "cpus_used": cpus,
            "siblings": {str(c): sibling(c) for c in cpus}, "load1": os.getloadavg()[0],
            "hostname": os.uname().nodename, "form_version": fv, "form_path": shutil.which(form) or form,
            "form_sha256": sha256(Path(shutil.which(form))) if shutil.which(form) else None}


def git_info(path):
    def g(*a):
        try:
            return subprocess.run(["git", "-C", str(path), *a], capture_output=True, text=True).stdout.strip()
        except Exception:
            return None
    return {"commit": g("rev-parse", "HEAD"), "branch": g("rev-parse", "--abbrev-ref", "HEAD"),
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
    ap.add_argument("--contention-retries", type=int, default=2)
    ap.add_argument("--draft", action="store_true")
    args = ap.parse_args()
    args.binary = str(Path(args.binary).resolve())
    args.data_root = args.data_root.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    cases_doc = json.loads(args.cases.read_text())
    cpus = expand_cpus(args.cpus) if args.cpus else sorted(os.sched_getaffinity(0))
    started = time.time()
    pipeline_file = out / "pipeline.json"
    pipeline = json.loads(pipeline_file.read_text()) if pipeline_file.exists() else {"runs": []}

    ident_path = out / "identity.json"
    if not ident_path.exists() or not args.assemble_only:
        subprocess.run([args.binary, "identity", "--data-root", str(args.data_root), "--out", str(ident_path)], check=True)
    identity = json.loads(ident_path.read_text())
    identity["host"] = host_info(cpus, args.form)
    identity["source"] = git_info(CRATE)
    identity["cases_sha256"] = sha256(args.cases)
    identity["tools_sha256"] = {p.name: sha256(p) for p in sorted(HERE.glob("*.py"))}

    run = {"started_utc": now(), "argv": sys.argv, "cpus": cpus}
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
        for case in selected:
            p = prepare(args, case, out)
            if p:
                inputs[case["id"]] = p
        run["prepare_seconds"] = time.time() - t
        if args.prepare_only:
            print(f"prepared {len(inputs)} inputs in {run['prepare_seconds']:.1f} s", flush=True)
            return
        jobs = queue.Queue()
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
                    jobs.put((case, cell, rep))
        t = time.time()

        def worker(cpu):
            while True:
                try:
                    case, cell, rep = jobs.get_nowait()
                except queue.Empty:
                    return
                measure(args, case, cell, rep, cpu, out, inputs[case["id"]])

        threads = [threading.Thread(target=worker, args=(c,)) for c in cpus]
        for th in threads:
            th.start()
        for th in threads:
            th.join()
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
