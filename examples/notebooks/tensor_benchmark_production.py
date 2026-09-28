"""Paired production workloads for the tensor consolidation benchmark driver."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import shutil
import signal
import subprocess
from pathlib import Path
from statistics import median
from time import perf_counter_ns, time_ns

CASES = ("GL000", "GL053", "GL148", "uv-scalars")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(command, directory, output, cpu, timeout, memory_limit):
    """Measure one fresh process, including startup, without blocking on log lines."""

    def alarm(_signal, _frame):
        raise TimeoutError

    command = [
        shutil.which("prlimit"),
        f"--as={memory_limit}:{memory_limit}",
        "--",
        shutil.which("taskset"),
        "-c",
        str(cpu),
        *command,
    ]
    timed_out = False
    log_path = output / "run.log"
    with log_path.open("wb") as log:
        started = perf_counter_ns()
        process = subprocess.Popen(
            command,
            cwd=directory,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        previous = signal.signal(signal.SIGALRM, alarm)
        signal.setitimer(signal.ITIMER_REAL, timeout)
        try:
            _, status, usage = os.wait4(process.pid, 0)
        except TimeoutError:
            timed_out = True
            os.killpg(process.pid, signal.SIGKILL)
            _, status, usage = os.wait4(process.pid, 0)
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, previous)
        wall_ns = perf_counter_ns() - started
        process.returncode = os.waitstatus_to_exitcode(status)
    return {
        "command": command,
        "exit_code": process.returncode,
        "timed_out": timed_out,
        "wall_ns": wall_ns,
        "cpu_ns": round((usage.ru_utime + usage.ru_stime) * 1e9),
        "user_seconds": usage.ru_utime,
        "system_seconds": usage.ru_stime,
        "max_rss_kib": usage.ru_maxrss,
        "log": str(log_path),
        "log_sha256": sha256(log_path),
    }


def benchmark(args):
    if args.cpu is None or args.cpu not in os.sched_getaffinity(0):
        raise ValueError("Production comparisons require an available --cpu")
    if args.rounds < 1 or args.timeout < 1:
        raise ValueError("Rounds and timeout must be positive")
    fixture_root = args.fixture_root.resolve()
    runner = fixture_root / "examples/cli/aa_aa/3L/run_grouped_diagram_smokes.py"
    spec = importlib.util.spec_from_file_location("grouped_smoke", runner)
    grouped = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(grouped)
    cases = args.cases.split(",") if args.cases else list(CASES)
    if (
        not cases
        or any(case not in CASES for case in cases)
        or len(set(cases)) != len(cases)
    ):
        raise ValueError(f"Select distinct production cases from {CASES}")
    builds = {}
    inputs = {str(runner): sha256(runner), str(Path(__file__)): sha256(__file__)}
    for tool in ("prlimit", "taskset"):
        path = shutil.which(tool)
        if path is None:
            raise FileNotFoundError(f"Production measurements require {tool}")
        inputs[path] = sha256(path)
    for entry in args.production_build or []:
        label, path = entry.split("=", 1)
        if label in builds:
            raise ValueError(f"Duplicate build label: {label}")
        path = Path(path).resolve()
        results_path = path / "build-results.json"
        results = json.loads(results_path.read_text())
        binaries = {}
        for kind in ("test-lib", "cli"):
            row = next(row for row in results if row["label"] == kind)
            if row["exit_code"] or len(row["artifacts"]) != 1:
                raise ValueError(
                    f"Build {label}/{kind} is not a successful single artifact"
                )
            artifact = row["artifacts"][0]
            inputs[artifact["path"]] = artifact["sha256"]
            binaries[kind] = artifact
        manifest_path = path / "source-manifest.json"
        manifest = json.loads(manifest_path.read_text())
        inputs[str(manifest_path)] = sha256(manifest_path)
        inputs[str(results_path)] = sha256(results_path)
        inputs.update(
            {str(path / "source" / p): h for p, h in manifest["files"].items()}
        )
        builds[label] = {
            "path": str(path),
            "binaries": binaries,
            "source_manifest_sha256": sha256(manifest_path),
        }
    if len(builds) != 2:
        raise ValueError(
            "Supply exactly two --production-build LABEL=BUILD_DIRECTORY values"
        )
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(output)
    raw = output.with_suffix("")
    raw.mkdir(parents=True, exist_ok=False)
    fixture_args = argparse.Namespace(
        n_max=20,
        n_start=20,
        n_increase=0,
        cores=1,
        generate_cores=1,
        batch_size=1,
        phase="imag",
        target_relative_accuracy=1.0,
    )
    graphs = {}
    provenance = fixture_root / "examples/notebooks/fixtures/production/provenance.json"
    fixture_provenance = json.loads(provenance.read_text())
    inputs[str(provenance)] = sha256(provenance)
    for case in cases:
        if case == "uv-scalars":
            continue
        path = fixture_root / "examples/notebooks/fixtures/production" / f"{case}.dot"
        graphs[case] = path
        inputs[str(path)] = sha256(path)
        if inputs[str(path)] != fixture_provenance["records"][case]["prepared_sha256"]:
            raise ValueError(f"Canonical production fixture changed: {case}")
    labels = list(builds)
    schedule = [("warmup", 0, case, label) for case in cases for label in labels]
    schedule += [
        ("measurement", round_, case, label)
        for round_ in range(args.rounds)
        for case in cases
        for label in (labels if round_ % 2 == 0 else labels[::-1])
    ]
    report = {
        "started_ns": time_ns(),
        "completed": False,
        "builds": builds,
        "inputs": inputs,
        "fixture_provenance": fixture_provenance,
        "schedule": schedule,
        "records": [],
        "protocol": {
            "cpu": args.cpu,
            "timeout_seconds": args.timeout,
            "memory_limit_bytes": 64 * 1024**3,
            "rounds": args.rounds,
            "warmups_per_case_and_build": 1,
            "clock": "complete fresh command; startup included; preparation and checks excluded",
            "cases": cases,
            "samples": 20,
            "seed": 1337,
            "generation_cores": 1,
            "integration_cores": 1,
            "comparison_relative_tolerance": 1e-9,
            "comparison_absolute_tolerance": 1e-12,
            "scope": "individual three-loop aa→aa graphs; scalar UV HedgePoset profile",
        },
    }

    def save():
        output.write_text(json.dumps(report, indent=2) + "\n")

    def guard():
        changed = [
            path for path, expected in inputs.items() if sha256(path) != expected
        ]
        if changed:
            raise RuntimeError(f"Benchmark inputs changed: {changed}")

    save()
    guard()
    for sequence, (phase, round_, case, label) in enumerate(schedule):
        work = raw / f"{sequence:03d}-{phase}-{label}-{case}"
        work.mkdir()
        build = builds[label]
        if case == "uv-scalars":
            binary = build["binaries"]["test-lib"]["path"]
            command = [
                binary,
                "--exact",
                "uv::tests::scalars_profile_new",
                "--nocapture",
                "--test-threads",
                "1",
            ]
        else:
            card = work / "run.toml"
            card.write_text(
                grouped.toml_card(
                    graphs[case], work / "state", work / "workspace", fixture_args
                )
            )
            binary = build["binaries"]["cli"]["path"]
            command = [binary, "--no-save-state", "--clean-state", str(card)]
        guard()
        row = run(
            command,
            Path(build["path"]) / "source",
            work,
            args.cpu,
            args.timeout,
            report["protocol"]["memory_limit_bytes"],
        )
        row.update(phase=phase, round=round_, case=case, label=label)
        guard()
        text = Path(row["log"]).read_text(errors="replace")
        if case == "uv-scalars":
            passed = "1 passed; 0 failed" in text
        else:
            row["log_summary"] = grouped.summarize_log(text)
            result = row["result"] = grouped.summarize_result(work / "workspace")
            passed = (
                row["log_summary"].get("generated", False)
                and result.get("result_emitted", False)
                and result.get("neval") == 20
                and all(
                    math.isfinite(float(result[key]))
                    for key in ("result_re", "result_im", "error_re", "error_im")
                )
            )
            result_path = work / "workspace/integration_result.json"
            if result_path.exists():
                row["integration_result_sha256"] = sha256(result_path)
        row["checks_pass"] = passed and row["exit_code"] == 0 and not row["timed_out"]
        report["records"].append(row)
        save()
        print(
            f"{phase} {case} {label}: exit={row['exit_code']} checks={row['checks_pass']}",
            flush=True,
        )
        if not row["checks_pass"]:
            raise RuntimeError(
                f"Production qualification failed; complete record retained in {output}"
            )
    comparisons = []
    for round_ in range(args.rounds):
        for case in cases:
            rows = {
                row["label"]: row
                for row in report["records"]
                if row["phase"] == "measurement"
                and row["round"] == round_
                and row["case"] == case
            }
            equal = True
            if case != "uv-scalars":
                left, right = (rows[label]["result"] for label in labels)
                equal = all(
                    math.isclose(
                        float(left[key]), float(right[key]), rel_tol=1e-9, abs_tol=1e-12
                    )
                    for key in ("result_re", "result_im", "error_re", "error_im")
                )
                coordinates = {
                    key
                    for key in left.keys() | right.keys()
                    if key.endswith("_coordinates")
                }
                equal &= bool(coordinates) and all(
                    left.get(key) == right.get(key) for key in coordinates
                )
            comparisons.append({"case": case, "round": round_, "pass": equal})
    report["comparisons"] = comparisons
    report["summary"] = {
        case: {
            label: {
                clock.replace("_ns", "_ms"): median(
                    row[clock] / 1e6
                    for row in report["records"]
                    if row["phase"] == "measurement"
                    and row["case"] == case
                    and row["label"] == label
                )
                for clock in ("wall_ns", "cpu_ns")
            }
            for label in labels
        }
        for case in cases
    }
    report.update(
        completed=True,
        source_drift=False,
        checks_pass=all(row["pass"] for row in comparisons),
    )
    save()
    if not report["checks_pass"]:
        raise AssertionError(
            f"Production comparison failed; all rows retained in {output}"
        )
    return report
