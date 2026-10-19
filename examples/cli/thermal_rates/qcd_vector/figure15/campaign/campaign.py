#!/usr/bin/env python3
"""Prepare or run the bounded Figure 15 scan. Running requires an explicit 'run'."""

import argparse
import collections
import csv
import datetime as dt
import math
import os
import shlex
import signal
import subprocess
import sys
import time
import traceback
from pathlib import Path

from common import (
    append_event,
    file_hash,
    group_members,
    process_identity,
    read_json,
    same_process,
    set_process_affinity,
    signal_process,
    state_identity,
    utc_now,
    write_json,
)

HERE = Path(__file__).resolve().parent
COMMIT = "116767bd3ed30491a2b15559d228171d20082abe"


def points():
    result = [
        (energy / 4, momentum / 2)
        for energy in range(1, 9)
        for momentum in range(4)
        if energy / 4 != momentum / 2
    ]
    result += [
        (momentum + offset, momentum)
        for momentum in (0.5, 1.0, 1.5)
        for offset in (-0.05, 0.05)
    ]
    assert len(result) == len(set(result)) == 35
    return result


def point_name(energy, momentum):
    return f"E{energy:.2f}_k{momentum:.2f}".replace(".", "p")


def next_quota(samples):
    for quota in (100000, 250000, 500000):
        if samples < quota:
            return quota
    return ((samples // 500000) + 1) * 500000


def prepare(args):
    root = args.run_root.resolve()
    repository = args.repository.resolve()
    binary = (
        args.binary.resolve()
        if args.binary
        else repository / "target/dev-optim/gammaloop"
    )
    direct_state, rest_state = args.direct_state.resolve(), args.rest_state.resolve()
    available = sorted(os.sched_getaffinity(0))
    cpus = (
        [int(value) for value in args.cpus.split(",")]
        if args.cpus
        else available[64:124]
    )
    if not cpus and len(available) >= args.worker_limit:
        cpus = available[: args.worker_limit]
    if (
        len(cpus) < args.worker_limit
        or len(set(cpus)) != len(cpus)
        or not set(cpus) <= set(available)
    ):
        raise ValueError(
            "Need distinct available CPUs covering the worker limit; specify --cpus explicitly"
        )
    cpus = cpus[: args.worker_limit]
    selected = points()
    if args.test_point:
        selected = [tuple(map(float, value.split(","))) for value in args.test_point]
    if args.test_point_limit:
        selected = selected[: args.test_point_limit]
    jobs = []
    for index, (energy, momentum) in enumerate(selected):
        name = point_name(energy, momentum)
        directory = root / "jobs" / name
        rest = momentum == 0
        source = rest_state if rest else direct_state
        processes = ["vector_nlo"] if rest else ["vector_T_nlo", "vector_L_nlo"]
        seed = args.seed + index * 104729
        settings = []
        for process in processes:
            settings.append(
                f"set process -p {process} -i NLO kv "
                f"kinematics.externals.data.momenta=[[{energy},0.0,0.0,{momentum}]] "
                f"kinematics.e_cm={max(math.sqrt(abs(energy * energy - momentum * momentum)), 0.25):.17g} "
                f"general.generate_events=false integrator.n_start={args.iteration_samples} "
                f"integrator.n_increase=0 integrator.n_max=1000000000 integrator.seed={seed} "
                "integrator.integrated_phase=real integrator.target_relative_accuracy=1e-12"
            )
        selectors = " ".join(f"-p {process} -i NLO" for process in processes)
        integrate = (
            f"integrate {selectors} --n-cores {args.workers_per_job} --batch-timing 0.5 "
            f"--workspace-path {directory}/workspace --renderer tabled "
            "--no-stream-iterations --no-stream-updates --write-results-for-each-iteration"
        )
        commands = settings[:]
        for epoch in range(args.max_cycles):
            commands.append(
                f"! {shlex.quote(sys.executable)} {shlex.quote(str(HERE / 'gate.py'))} "
                f"{shlex.quote(str(root))} {shlex.quote(name)} {epoch}"
            )
            commands.append(integrate)
        # This final callback lets the controller label and terminate a cycle-limited session deliberately.
        commands.append(
            f"! {shlex.quote(sys.executable)} {shlex.quote(str(HERE / 'gate.py'))} "
            f"{shlex.quote(str(root))} {shlex.quote(name)} {args.max_cycles}"
        )
        command = [
            str(binary),
            "-n",
            "--read-only-state",
            "-s",
            str(source),
            "run",
            "-c",
            "; ".join(commands),
        ]
        if max(map(len, command)) >= 120000:
            raise ValueError("CLI argument too long; reduce --max-cycles")
        jobs.append(
            {
                "name": name,
                "k0_over_T": energy,
                "k_over_T": momentum,
                "T": 1,
                "e": 1,
                "g": 1,
                "Nc": 3,
                "Qd": -1 / 3,
                "seed": seed,
                "source_state": str(source),
                "directory": str(directory),
                "processes": processes,
                "slot_keys": [f"{process}@NLO" for process in processes],
                "output_scale": 3 * energy if rest else 9 * energy,
                "rest_isotropy": rest,
                "command": command,
                "workers": args.workers_per_job,
            }
        )
    reference_command = (
        shlex.split(args.reference_command) if args.reference_command else None
    )
    plan = {
        "source_commit": args.source_commit,
        "repository": str(repository),
        "binary": str(binary),
        "binary_sha256": file_hash(binary),
        "source_states": [state_identity(path) for path in (direct_state, rest_state)],
        "controller_sources": [
            {"path": str(path), "sha256": file_hash(path)}
            for path in sorted(HERE.glob("*.py"))
        ],
        "run_root": str(root),
        "python": sys.executable,
        "point_count": len(jobs),
        "jobs": jobs,
        "worker_limit": args.worker_limit,
        "workers_per_job": args.workers_per_job,
        "cpu_affinity": cpus,
        "relative_error": args.relative_error,
        "absolute_floor": args.absolute_floor,
        "minimum_samples": args.minimum_samples,
        "wall_seconds": args.wall_seconds,
        "graceful_margin_seconds": args.graceful_margin_seconds,
        "quantum_seconds": args.quantum_seconds,
        "iteration_samples": args.iteration_samples,
        "max_cycles": args.max_cycles,
        "reference_command": reference_command,
        "reserved_workers": args.reserved_workers,
        "reference_targets_json": str(args.reference_targets_json.resolve())
        if args.reference_targets_json
        else None,
        "production_settings": not args.allow_test_settings,
        "stopping_criterion": "error <= max(relative_error * abs(reference if present else estimate), absolute_floor), "
        "minimum training samples, zero final discarded or unstable samples; reference residual never controls stopping",
        "checkpoint_policy": "Completed-iteration native snapshots only; no metadata edits or on-disk state reloads",
        "prepared_utc": utc_now(),
    }
    return plan


def reference_targets(path):
    payload = read_json(path, []) if path else []
    if isinstance(payload, dict):
        payload = payload.get("rows", payload.get("references", []))
    result = {}
    for row in payload:
        try:
            energy = float(row.get("k0_over_T", row.get("E")))
            momentum = float(row.get("k_over_T", row.get("k")))
            component = row["component"].upper()
            value = float(row.get("rho_nlo", row.get("value")))
            if math.isfinite(value):
                result[(round(energy, 9), round(momentum, 9), component)] = value
        except (TypeError, ValueError, KeyError):
            continue
    return result


def checkpoint(job, plan, targets):
    path = Path(job["directory"]) / "workspace/integration_result.json"
    raw = read_json(path)
    if not raw:
        for archive in sorted(
            path.parent.glob("results/integration_result_iter_*.json"), reverse=True
        ):
            raw = read_json(archive)
            if raw:
                break
    if not raw or not raw.get("slots"):
        return {
            "samples": 0,
            "precision_reached": False,
            "invalid": False,
            "rows": [],
            "result_path": str(path),
        }
    rows, precision, provisional_precision, invalid = [], True, True, False
    counts = []
    slots = {slot["key"]: slot for slot in raw["slots"]}
    for key in job["slot_keys"]:
        slot = slots.get(key)
        if not slot:
            return {
                "samples": 0,
                "precision_reached": False,
                "invalid": True,
                "rows": [],
                "result_path": str(path),
                "error": "missing native slot",
            }
        integral, stats = slot["integral"], slot.get("integration_statistics", {})
        samples, scale = int(integral["neval"]), job["output_scale"]
        counts.append(samples)
        value = integral["result"]["re"] * scale
        error = abs(integral["error"]["re"] * scale)
        finite = math.isfinite(value) and math.isfinite(error) and error >= 0
        clean = (
            bool(stats)
            and stats.get("nan_percentage", 1) == 0
            and stats.get("nan_or_unstable_percentage", 1) == 0
        )
        invalid |= not finite or not clean
        components = (
            ("T", "L")
            if job["rest_isotropy"]
            else ("T",)
            if key == "vector_T_nlo@NLO"
            else ("L",)
        )
        for component in components:
            target = targets.get(
                (round(job["k0_over_T"], 9), round(job["k_over_T"], 9), component)
            )
            accepted_reference = job.get("acceptance_reference", {}).get(component)
            precision_reference = (
                accepted_reference if accepted_reference is not None else target
            )
            norm = (
                abs(precision_reference)
                if precision_reference is not None
                else abs(value)
            )
            absolute_goal = max(plan["relative_error"] * norm, plan["absolute_floor"])
            provisional = (
                finite
                and clean
                and samples >= plan["minimum_samples"]
                and absolute_goal > 0
                and error <= absolute_goal
            )
            reached = provisional and precision_reference is not None
            provisional_precision &= provisional
            precision &= reached
            table = next(
                (
                    item
                    for item in slot.get("table_results", [])
                    if item["component"] == "re"
                ),
                {},
            )
            rows.append(
                {
                    "name": job["name"],
                    "component": component,
                    "order": "nlo",
                    "k0_over_T": job["k0_over_T"],
                    "k_over_T": job["k_over_T"],
                    "samples": samples,
                    "value": value if finite else None,
                    "error": error if finite else None,
                    "reference": target,
                    "pull": (value - target) / error
                    if finite and target is not None and error
                    else None,
                    "relative_error": abs(error / value) if finite and value else None,
                    "relative_error_to_reference": abs(error / target)
                    if finite and target
                    else None,
                    "precision_norm_source": "frozen acceptance reference"
                    if accepted_reference is not None
                    else "independent reference"
                    if target is not None
                    else "current estimate; provisional only",
                    "precision_reference": precision_reference,
                    "independent_reference_available": target is not None,
                    "target_absolute_error": absolute_goal
                    if math.isfinite(absolute_goal)
                    else None,
                    "absolute_floor": plan["absolute_floor"],
                    "precision_reached": reached,
                    "provisional_precision_reached": provisional,
                    "final_stability_clean": clean,
                    "final_nonfinite_or_unstable_count": round(
                        samples * stats.get("nan_or_unstable_percentage", 0) / 100
                    ),
                    "statistics": stats,
                    "chi_sq_per_dof": table.get("chi_sq_per_dof"),
                    "max_weight_impact": table.get("max_weight_impact"),
                    "native_slot": key,
                    "native_multiplier": scale,
                    "result_path": str(path),
                    "rest_isotropy": job["rest_isotropy"],
                }
            )
    if len(set(counts)) != 1:
        invalid = True
    return {
        "samples": min(counts),
        "precision_reached": precision and not invalid,
        "provisional_precision_reached": provisional_precision and not invalid,
        "invalid": invalid,
        "rows": rows,
        "result_path": str(path),
    }


def report(root, plan, state):
    targets = reference_targets(plan.get("reference_targets_json"))
    rows = []
    for job in state["jobs"]:
        saved = checkpoint(job, plan, targets)
        job["samples"] = saved["samples"]
        job["precision_reached"] = saved["precision_reached"]
        job["final_stability_clean"] = bool(saved["rows"]) and not saved["invalid"]
        for row in saved["rows"]:
            row["status"] = job["status"]
            current = job.get("current_allocation")
            row["active_wall_seconds"] = job.get("active_wall_seconds", 0) + (
                max(0, time.monotonic() - current["reserved_monotonic"])
                if current
                else 0
            )
            row["workers"] = plan["workers_per_job"]
            rows.append(row)
    state["updated_utc"] = utc_now()
    state["summary_rows"] = len(rows)
    write_json(root / "manifest.json", state)
    write_json(root / "summary.json", rows)
    columns = [
        "name",
        "component",
        "k0_over_T",
        "k_over_T",
        "samples",
        "value",
        "error",
        "reference",
        "pull",
        "relative_error",
        "relative_error_to_reference",
        "precision_norm_source",
        "precision_reference",
        "independent_reference_available",
        "target_absolute_error",
        "precision_reached",
        "provisional_precision_reached",
        "final_stability_clean",
        "final_nonfinite_or_unstable_count",
        "active_wall_seconds",
        "status",
        "result_path",
    ]
    temporary = root / "summary.csv.tmp"
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, root / "summary.csv")
    return rows


def run(plan, start_times=None):
    root = Path(plan["run_root"])
    if (root / "campaign.json").exists():
        raise ValueError(
            "A campaign already exists at this path; use a new run root. No automatic restart/reload."
        )
    identity = process_identity(os.getpid())
    if identity["process_group"] != os.getpid() or identity["session"] != os.getpid():
        raise ValueError(
            "Launch controller with start_new_session=True so its complete process group can be bounded safely"
        )
    root.mkdir(parents=True, exist_ok=True)
    started, monotonic = start_times or (time.time(), time.monotonic())
    campaign = {
        "started_utc": dt.datetime.fromtimestamp(started, dt.UTC).isoformat(),
        "started_epoch": started,
        "started_monotonic": monotonic,
        "controller_identity": identity,
        "compute_process_group": os.getpgrp(),
        "hard_deadline_utc": dt.datetime.fromtimestamp(
            started + plan["wall_seconds"], dt.UTC
        ).isoformat(),
        "hard_deadline_monotonic": monotonic + plan["wall_seconds"],
        "soft_deadline_monotonic": monotonic
        + plan["wall_seconds"]
        - plan["graceful_margin_seconds"],
    }
    write_json(root / "campaign.json", campaign)
    write_json(root / "plan.json", plan)
    watchdog_log = (root / "watchdog.log").open("w")
    watchdog = subprocess.Popen(
        [sys.executable, str(HERE / "watchdog.py"), str(root)],
        stdout=watchdog_log,
        stderr=subprocess.STDOUT,
        start_new_session=True,
    )
    for _ in range(100):
        registration = read_json(root / "watchdog.json", {})
        if (
            registration.get("status") == "armed"
            and registration.get("target_group") == os.getpgrp()
        ):
            break
        if watchdog.poll() is not None:
            raise RuntimeError(
                "Independent watchdog failed before any compute was launched"
            )
        time.sleep(0.05)
    else:
        raise RuntimeError(
            "Independent watchdog registration timed out; no compute launched"
        )
    append_event(
        root,
        "watchdog_verified",
        watchdog_pid=watchdog.pid,
        hard_deadline_utc=campaign["hard_deadline_utc"],
    )
    jobs = [
        {key: value for key, value in job.items() if key != "command"}
        for job in plan["jobs"]
    ]
    for job in jobs:
        job.update(status="pending", samples=0, active_wall_seconds=0, epochs=[])
        directory = Path(job["directory"])
        directory.mkdir(parents=True, exist_ok=True)
    state = {
        "stage": "full_scan",
        "source_commit": COMMIT,
        "started_utc": campaign["started_utc"],
        "hard_deadline_utc": campaign["hard_deadline_utc"],
        "compute_process_group": os.getpgrp(),
        "watchdog_pid": watchdog.pid,
        "watchdog_verified": True,
        "worker_limit": plan["worker_limit"],
        "relative_error": plan["relative_error"],
        "absolute_floor": plan["absolute_floor"],
        "minimum_samples": plan["minimum_samples"],
        "jobs": jobs,
        "full_scan_started": True,
        "metadata_fingerprints_modified": False,
        "completion_notification": "Local JSON/CSV artifacts only; no automatic chat notification configured",
    }
    write_json(root / "reserved_workers.json", {"workers": plan["reserved_workers"]})
    references = None
    reference_log = None
    if plan["reference_command"]:
        reserved_cpus = plan["cpu_affinity"][-6:]
        reference_log = (root / "reference.log").open("w")
        references = subprocess.Popen(
            [
                "taskset",
                "-c",
                ",".join(map(str, reserved_cpus)),
                *plan["reference_command"],
            ],
            cwd=plan["repository"],
            stdout=reference_log,
            stderr=subprocess.STDOUT,
        )
        state["reference_helper"] = {
            "command": plan["reference_command"],
            "pid": references.pid,
            "cpu_affinity": reserved_cpus,
            "started_utc": utc_now(),
            "status": "running",
        }
        write_json(
            root / "reserved_workers.json",
            {"workers": max(6, plan["reserved_workers"])},
        )
    handles, processes, active = {}, {}, {}
    queue = collections.deque(job["name"] for job in jobs)
    by_name = {job["name"]: job for job in jobs}
    commands = {job["name"]: job["command"] for job in plan["jobs"]}
    cpu_slots = [
        plan["cpu_affinity"][index : index + plan["workers_per_job"]]
        for index in range(0, plan["worker_limit"], plan["workers_per_job"])
    ]
    shutting_down, last_report, sticky_reference_reservation, cleanup_started = (
        False,
        0,
        0,
        None,
    )

    def release(name, now):
        allocation = active.pop(name, None)
        if not allocation:
            return
        duration = now - allocation["reserved_monotonic"]
        job = by_name[name]
        job["active_wall_seconds"] += duration
        allocation["ended_utc"], allocation["wall_seconds"] = utc_now(), duration
        job["epochs"].append(allocation)
        job.pop("current_allocation", None)
        append_event(
            root,
            "allocation_released",
            name=name,
            epoch=allocation["epoch"],
            wall_seconds=duration,
        )

    def stop_job(job, status, reason):
        job.update(status=status, stop_reason=reason)
        write_json(
            Path(job["directory"]) / "stop.json", {"reason": reason, "utc": utc_now()}
        )
        if job["name"] in active:
            allocation = active[job["name"]]
            if not allocation.get("interrupt_sent"):
                signal_process(job.get("identity"), signal.SIGINT)
                allocation.update(
                    interrupt_sent=True,
                    interrupted_monotonic=time.monotonic(),
                    interrupt_reason=reason,
                )

    while True:
        now = time.monotonic()
        if (
            watchdog.poll() is not None
            and not (root / "all-compute-exited.json").exists()
        ):
            # Losing the guard cannot leave a 12-hour budget unenforced.
            try:
                os.write(
                    2,
                    b"Independent watchdog died; killing the complete compute process group.\n",
                )
            finally:
                os.killpg(os.getpgrp(), signal.SIGKILL)
        if references is not None and references.poll() is not None:
            state["reference_helper"].update(
                status="completed" if references.returncode == 0 else "failed",
                exit_code=references.returncode,
                ended_utc=utc_now(),
            )
            reference_log.close()
            if references.returncode != 0:
                # A failed helper may leave quadrature children in the guarded group.
                # Preserve their six-worker reservation until the global hard cutoff.
                sticky_reference_reservation = 6
                state["reference_helper"]["worker_reservation_retained"] = True
            references = None
            write_json(
                root / "reserved_workers.json",
                {"workers": sticky_reference_reservation},
            )
            append_event(root, "reference_helper_finished", **state["reference_helper"])
        if now >= campaign["soft_deadline_monotonic"] and not shutting_down:
            shutting_down = True
            write_json(
                root / "stop.json",
                {"reason": "original_wall_time_budget", "utc": utc_now()},
            )
            queue.clear()
            for job in jobs:
                if job["status"] not in (
                    "precision_reached",
                    "failed",
                    "invalid_samples",
                    "cycle_limit",
                ):
                    stop_job(job, "budget_stopped", "original_wall_time_budget")
            if references is not None:
                references.terminate()
            append_event(root, "graceful_deadline_started")
        targets = reference_targets(plan.get("reference_targets_json"))
        for job in jobs:
            name, directory = job["name"], Path(job["directory"])
            process = processes.get(name)
            if process is None:
                continue
            code = process.poll()
            if code is not None:
                release(name, now)
                if job["status"] in ("running", "parked", "starting"):
                    job.update(status="failed", stop_reason="unexpected_cli_exit")
                job.update(exit_code=code, ended_utc=utc_now())
                handles.pop(name).close()
                processes.pop(name)
                append_event(
                    root,
                    "native_exited",
                    name=name,
                    status=job["status"],
                    exit_code=code,
                )
                continue
            saved = checkpoint(job, plan, targets)
            job["samples"] = saved["samples"]
            gate = read_json(directory / "gate.json", {})
            allocation = active.get(name)
            returned = bool(
                allocation
                and gate.get("epoch", -1) > allocation["epoch"]
                and gate.get("state") in ("ready", "stopping")
            )
            if returned:
                release(name, now)
                job["next_epoch"] = gate["epoch"]
                if job["status"] in ("running", "starting"):
                    job["status"] = "parked"
            if saved["invalid"] and job["status"] not in ("failed", "invalid_samples"):
                stop_job(
                    job, "invalid_samples", "nonfinite_or_unstable_completed_samples"
                )
            if saved["precision_reached"] and job["status"] in (
                "running",
                "parked",
                "starting",
            ):
                job["acceptance_reference"] = {
                    row["component"]: row["precision_reference"]
                    for row in saved["rows"]
                }
                stop_job(job, "precision_reached", "isolated_nlo_precision_goal")
            if gate.get("state") in ("ready", "stopping"):
                if (
                    job["status"] not in ("running", "parked", "starting")
                    or shutting_down
                ):
                    write_json(
                        directory / "stop.json",
                        {
                            "reason": job.get("stop_reason", "campaign_stopping"),
                            "utc": utc_now(),
                        },
                    )
                elif gate["epoch"] >= plan["max_cycles"]:
                    stop_job(job, "cycle_limit", "maximum_persistent_session_cycles")
                elif name not in active and name not in queue:
                    job.update(status="parked", next_epoch=gate["epoch"])
                    queue.append(name)
            allocation = active.get(name)
            if not allocation:
                continue
            if (
                job["status"] == "starting"
                and gate.get("epoch") == allocation["epoch"]
                and gate.get("state") == "ready"
            ):
                allocation["started_monotonic"] = now
                job["status"] = "running"
                write_json(
                    directory / "grant.json",
                    {
                        "epoch": allocation["epoch"],
                        "native_start_ticks": job["identity"]["start_ticks"],
                        "quota": allocation["quota"],
                        "cpu_affinity": allocation["cpu_affinity"],
                        "utc": utc_now(),
                    },
                )
                append_event(
                    root,
                    "epoch_granted",
                    name=name,
                    epoch=allocation["epoch"],
                    quota=allocation["quota"],
                )
            if job["status"] == "running" and not allocation.get("interrupt_sent"):
                quota = saved["samples"] >= allocation["quota"]
                quantum = (
                    now - allocation["started_monotonic"] >= plan["quantum_seconds"]
                )
                if quota or quantum:
                    # A PID-only signal cannot kill a gate callback that won a completion race.
                    signal_process(job["identity"], signal.SIGINT)
                    allocation.update(
                        interrupt_sent=True,
                        interrupted_monotonic=now,
                        interrupt_reason="checkpoint_quota"
                        if quota
                        else "wall_quantum",
                    )
                    append_event(
                        root,
                        "quota_interrupt",
                        name=name,
                        epoch=allocation["epoch"],
                        samples=saved["samples"],
                        reason=allocation["interrupt_reason"],
                    )
            if (
                allocation.get("interrupt_sent")
                and now - allocation["interrupted_monotonic"] > 60
            ):
                signal_process(job.get("identity"), signal.SIGTERM)
                job.update(
                    status="failed" if not shutting_down else "budget_stopped",
                    stop_reason="interrupt_unresponsive_60_seconds",
                )
        reserved = read_json(root / "reserved_workers.json", {}).get("workers", 0)
        reserved = max(
            int(reserved),
            6 if references is not None else 0,
            sticky_reference_reservation,
        )
        capacity = max(0, (plan["worker_limit"] - reserved) // plan["workers_per_job"])
        state["reserved_workers"], state["active_workers"] = (
            reserved,
            len(active) * plan["workers_per_job"],
        )
        while queue and len(active) < capacity and not shutting_down:
            name = queue.popleft()
            job = by_name[name]
            if job["status"] not in ("pending", "parked"):
                continue
            used = {tuple(item["cpu_affinity"]) for item in active.values()}
            affinity = next(
                (slot for slot in cpu_slots[:capacity] if tuple(slot) not in used), None
            )
            if affinity is None:
                break
            epoch = job.get("next_epoch", 0)
            allocation = {
                "epoch": epoch,
                "cpu_affinity": affinity,
                "workers": plan["workers_per_job"],
                "quota": next_quota(job["samples"]),
                "started_utc": utc_now(),
                "started_monotonic": now,
                "reserved_monotonic": now,
                "initial_samples": job["samples"],
                "interrupt_sent": False,
            }
            if job["status"] == "pending":
                log = (Path(job["directory"]) / "run.log").open("w")
                env = os.environ.copy()
                env.update(
                    RAYON_NUM_THREADS=str(plan["workers_per_job"]),
                    OMP_NUM_THREADS="1",
                    OPENBLAS_NUM_THREADS="1",
                    MKL_NUM_THREADS="1",
                    GL_DISPLAY_FILTER="info",
                    GL_LOGFILE_FILTER="off",
                )
                process = subprocess.Popen(
                    ["taskset", "-c", ",".join(map(str, affinity)), *commands[name]],
                    cwd=plan["repository"],
                    env=env,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                )
                processes[name], handles[name] = process, log
                job.update(
                    identity=process_identity(process.pid),
                    pid=process.pid,
                    started_utc=utc_now(),
                    status="starting",
                )
                if (
                    not job["identity"]
                    or job["identity"]["process_group"] != os.getpgrp()
                ):
                    raise RuntimeError(
                        "Native child escaped the guarded compute process group"
                    )
                write_json(
                    Path(job["directory"]) / "session.json",
                    {
                        "identity": job["identity"],
                        "name": name,
                        "compute_process_group": os.getpgrp(),
                        "utc": utc_now(),
                    },
                )
                append_event(
                    root,
                    "native_started",
                    name=name,
                    pid=process.pid,
                    cpu_affinity=affinity,
                )
            else:
                set_process_affinity(job["identity"], affinity)
                job["status"] = "running"
                write_json(
                    Path(job["directory"]) / "grant.json",
                    {
                        "epoch": epoch,
                        "native_start_ticks": job["identity"]["start_ticks"],
                        "quota": allocation["quota"],
                        "cpu_affinity": affinity,
                        "utc": utc_now(),
                    },
                )
                append_event(
                    root,
                    "epoch_granted",
                    name=name,
                    epoch=epoch,
                    quota=allocation["quota"],
                )
            active[name] = allocation
            job["current_allocation"], job["cpu_affinity"] = allocation, affinity
        state["active_allocations"] = active
        if now - last_report >= 10:
            report(root, plan, state)
            last_report = now
        pending = any(
            job["status"] in ("pending", "running", "parked", "starting")
            for job in jobs
        )
        if not processes and (shutting_down or not pending) and references is None:
            members = [
                item
                for item in group_members(os.getpgrp())
                if item["pid"] != os.getpid() and item["state"] != "Z"
            ]
            if not members:
                break
            cleanup_started = cleanup_started or now
            state["status"] = "cleaning_guarded_children"
            state["remaining_compute_group_members"] = members
            for member in members:
                signal_process(
                    member,
                    signal.SIGKILL if now - cleanup_started > 2 else signal.SIGTERM,
                )
        time.sleep(0.25)
    state.update(
        finished_utc=utc_now(),
        status="budget_finished"
        if shutting_down
        else "completed_with_reference_failure"
        if sticky_reference_reservation
        else "completed",
        active_workers=0,
    )
    state["remaining_compute_group_members"] = []
    report(root, plan, state)
    # All children have been reaped. The independent guard checks group membership before disarming.
    write_json(
        root / "all-compute-exited.json",
        {"utc": utc_now(), "controller_pid": os.getpid()},
    )
    append_event(root, "campaign_finished", status=state["status"])
    print(f"Campaign {state['status']}; results: {root / 'summary.csv'}", flush=True)
    return 0


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("plan", "validate", "run", "report"))
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--repository", type=Path, default=Path.cwd())
    parser.add_argument(
        "--binary",
        type=Path,
        help="Immutable CLI path; default REPOSITORY/target/dev-optim/gammaloop",
    )
    parser.add_argument(
        "--direct-state",
        type=Path,
        help="Saved direct transverse/longitudinal projector state",
    )
    parser.add_argument(
        "--rest-state",
        type=Path,
        help="Saved Feynman trace state for exact rest-frame isotropy",
    )
    parser.add_argument("--source-commit", default=COMMIT)
    parser.add_argument(
        "--cpus",
        help="Comma-separated CPU IDs; production default is affinity positions 64:124",
    )
    parser.add_argument(
        "--reference-command",
        help="Quoted argv string, launched without a shell in the compute group",
    )
    parser.add_argument("--reference-targets-json", type=Path)
    parser.add_argument("--reserved-workers", type=int, default=0)
    parser.add_argument("--absolute-floor", type=float, default=0)
    parser.add_argument("--relative-error", type=float, default=0.005)
    parser.add_argument("--minimum-samples", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=191007552)
    parser.add_argument("--worker-limit", type=int, default=60)
    parser.add_argument("--workers-per-job", type=int, default=6)
    parser.add_argument("--wall-seconds", type=float, default=43200)
    parser.add_argument("--graceful-margin-seconds", type=float, default=120)
    parser.add_argument("--quantum-seconds", type=float, default=600)
    parser.add_argument("--iteration-samples", type=int, default=2000)
    parser.add_argument("--max-cycles", type=int, default=128)
    parser.add_argument("--test-point-limit", type=int, default=0)
    parser.add_argument(
        "--test-point",
        action="append",
        help="Tiny actual-process smoke point E,k; requires --allow-test-settings",
    )
    parser.add_argument("--allow-test-settings", action="store_true")
    args = parser.parse_args()
    if args.action != "report" and (
        args.direct_state is None or args.rest_state is None
    ):
        parser.error(
            "--direct-state and --rest-state are required outside the read-only report action"
        )
    if (
        args.worker_limit <= 0
        or args.workers_per_job <= 0
        or args.worker_limit % args.workers_per_job
    ):
        parser.error("worker limit must be divisible by positive workers-per-job")
    if args.absolute_floor < 0 or args.relative_error <= 0 or args.minimum_samples <= 0:
        parser.error("invalid accuracy criterion")
    if not 0 <= args.reserved_workers <= args.worker_limit:
        parser.error("invalid reserved worker allocation")
    if (
        args.wall_seconds <= args.graceful_margin_seconds
        or args.quantum_seconds <= 0
        or args.iteration_samples <= 0
    ):
        parser.error("invalid wall/iteration settings")
    if not args.allow_test_settings and (
        args.wall_seconds != 43200
        or args.worker_limit != 60
        or args.workers_per_job != 6
        or args.test_point_limit
        or args.test_point
        or args.minimum_samples < 50000
    ):
        parser.error("nonproduction settings require --allow-test-settings")
    return args


def main():
    args = arguments()
    # Source hashing and launch setup consume the same wall budget as integration.
    start_times = (time.time(), time.monotonic()) if args.action == "run" else None
    if args.action == "report":
        root = args.run_root.resolve()
        plan, state = read_json(root / "plan.json"), read_json(root / "manifest.json")
        if not plan or not state:
            raise SystemExit("Existing plan and manifest are required")
        report(root, plan, state)
        return 0
    plan = prepare(args)
    if args.action == "validate":
        print(
            f"Validated {plan['point_count']} points, {plan['worker_limit']} workers, "
            f"{plan['wall_seconds'] / 3600:g} hours; no processes launched or campaign files written."
        )
        print(
            f"Criterion: {plan['stopping_criterion']}; absolute floor={plan['absolute_floor']}; minimum N={plan['minimum_samples']}"
        )
        return 0
    if args.action == "plan":
        write_json(args.run_root / "prepared-plan.json", plan)
        print(
            f"Prepared {args.run_root / 'prepared-plan.json'}; no processes launched."
        )
        return 0
    try:
        return run(plan, start_times)
    except BaseException as error:
        # The independent guard remains alive; fail closed immediately if this controller fails.
        campaign = read_json(args.run_root / "campaign.json", {})
        identity = campaign.get("controller_identity")
        if identity and identity["pid"] == os.getpid() and same_process(identity):
            try:
                write_json(
                    args.run_root / "controller-failure.json",
                    {
                        "utc": utc_now(),
                        "error": repr(error),
                        "traceback": traceback.format_exc(),
                    },
                )
                write_json(
                    args.run_root / "stop.json",
                    {"reason": "controller_failed", "utc": utc_now()},
                )
                print(traceback.format_exc(), flush=True)
            finally:
                os.killpg(os.getpgrp(), signal.SIGKILL)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
