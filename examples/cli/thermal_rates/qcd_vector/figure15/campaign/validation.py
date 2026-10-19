#!/usr/bin/env python3
"""Bounded Python-only checks, including a forced group cutoff of harmless sleep processes."""

import argparse
import datetime as dt
import json
import os
import signal
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import campaign
from common import (
    file_hash,
    group_members,
    process_identity,
    read_json,
    utc_now,
    write_json,
)


def main(output):
    checks = []
    for path in Path(__file__).resolve().parent.glob("*.py"):
        compile(path.read_text(), str(path), "exec")
    checks.append("all Python sources compile")
    assert len(campaign.points()) == 35 and len(set(campaign.points())) == 35
    assert all(
        energy > 0 and energy != momentum for energy, momentum in campaign.points()
    )
    assert [campaign.next_quota(value) for value in (0, 100000, 250000, 500000)] == [
        100000,
        250000,
        500000,
        1000000,
    ]
    checks.append("35 unique finite points and quota transitions")
    with tempfile.TemporaryDirectory(prefix="fig15-python-validation-") as temporary:
        root = Path(temporary)
        job = {
            "directory": str(root / "point"),
            "name": "synthetic",
            "slot_keys": ["vector_nlo@NLO"],
            "output_scale": 1,
            "rest_isotropy": True,
            "k0_over_T": 1,
            "k_over_T": 0,
            "status": "running",
            "active_wall_seconds": 0,
            "epochs": [],
        }
        slot = {
            "key": "vector_nlo@NLO",
            "integral": {
                "neval": 100000,
                "result": {"re": 1.1, "im": 0},
                "error": {"re": 0.004, "im": 0},
            },
            "integration_statistics": {
                "nan_percentage": 0,
                "nan_or_unstable_percentage": 0,
            },
        }
        write_json(
            Path(job["directory"]) / "workspace/integration_result.json",
            {"slots": [slot]},
        )
        plan = {
            "relative_error": 0.005,
            "absolute_floor": 0,
            "minimum_samples": 100000,
            "reference_targets_json": None,
            "workers_per_job": 2,
        }
        targets = {(1, 0, "T"): 1, (1, 0, "L"): 1}
        saved = campaign.checkpoint(job, plan, targets)
        assert saved["precision_reached"] and saved["rows"][0]["pull"] > 20
        unknown = campaign.checkpoint(job, plan, {})
        assert (
            not unknown["precision_reached"]
            and unknown["provisional_precision_reached"]
        )
        slot["integration_statistics"]["nan_or_unstable_percentage"] = 0.001
        write_json(
            Path(job["directory"]) / "workspace/integration_result.json",
            {"slots": [slot]},
        )
        assert not campaign.checkpoint(job, plan, targets)["precision_reached"]
        slot["integration_statistics"]["nan_or_unstable_percentage"] = 0
        slot["integral"]["neval"] = 10000
        write_json(
            Path(job["directory"]) / "workspace/integration_result.json",
            {"slots": [slot]},
        )
        assert not campaign.checkpoint(job, plan, targets)["precision_reached"]
        checks.append(
            "reference-normalized precision independent of pull; unknown reference, invalid sample and insufficient training rejected"
        )
        started, monotonic = time.time(), time.monotonic()
        stub = (
            "import json, pathlib, subprocess, sys, time; "
            "child=subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)']); "
            "pathlib.Path(sys.argv[1]).write_text(json.dumps({'pid':child.pid})); time.sleep(60)"
        )
        controller = subprocess.Popen(
            [sys.executable, "-c", stub, str(root / "stub-child.json")],
            start_new_session=True,
        )
        watchdog = None
        try:
            identity = process_identity(controller.pid)
            assert identity["process_group"] == identity["session"] == controller.pid
            assert controller.pid != os.getpgrp()
            deadline = monotonic + 5
            write_json(
                root / "campaign.json",
                {
                    "controller_identity": identity,
                    "compute_process_group": controller.pid,
                    "hard_deadline_monotonic": deadline,
                    "hard_deadline_utc": dt.datetime.fromtimestamp(
                        started + 5, dt.UTC
                    ).isoformat(),
                },
            )
            job["current_allocation"] = {"reserved_monotonic": monotonic, "epoch": 0}
            write_json(root / "plan.json", plan)
            write_json(root / "manifest.json", {"jobs": [job], "status": "running"})
            watchdog = subprocess.Popen(
                [
                    sys.executable,
                    str(Path(__file__).with_name("watchdog.py")),
                    str(root),
                ],
                start_new_session=True,
            )
            assert controller.wait(timeout=12) == -signal.SIGKILL
            assert watchdog.wait(timeout=12) == 0
            child_pid = read_json(root / "stub-child.json")["pid"]
            child = process_identity(child_pid)
            assert not child or child["state"] == "Z"
            assert not [
                member
                for member in group_members(controller.pid)
                if member["state"] != "Z"
            ]
            guard = read_json(root / "watchdog.json")
            final = read_json(root / "manifest.json")
            assert (
                guard["status"] == "hard_cutoff"
                and final["status"] == "budget_finished"
            )
            assert (
                final["finalized_by_independent_watchdog"]
                and final["jobs"][0]["status"] == "budget_stopped"
            )
            assert "current_allocation" not in final["jobs"][0]
            elapsed = final["jobs"][0]["active_wall_seconds"]
            campaign.report(root, plan, final)
            assert (
                read_json(root / "manifest.json")["jobs"][0]["active_wall_seconds"]
                == elapsed
            )
            observed = dt.datetime.fromisoformat(guard["cutoff_utc"]).timestamp()
            assert -0.05 <= observed - (started + 5) <= 1
            cutoff = {
                "controller_pid": controller.pid,
                "child_pid": child_pid,
                "controller_exit_code": controller.returncode,
                "target_group": controller.pid,
                "watchdog_pid": watchdog.pid,
                "all_group_compute_dead": True,
                "registered_deadline_utc": guard["hard_deadline_utc"],
                "observed_cutoff_utc": guard["cutoff_utc"],
                "deadline_delay_seconds": observed - (started + 5),
                "final_manifest_status": final["status"],
                "live_elapsed_frozen_after_cutoff": True,
            }
            checks.append(
                "forced five-second deadline kills controller and child process group and finalizes stable JSON/CSV"
            )
        finally:
            if controller.poll() is None:
                os.killpg(controller.pid, signal.SIGKILL)
                controller.wait()
            if watchdog is not None and watchdog.poll() is None:
                watchdog.kill()
                watchdog.wait()
    certificate = {
        "status": "passed",
        "completed_utc": utc_now(),
        "checks": checks,
        "forced_hard_cutoff": cutoff,
        "benchmark_processes_launched": 0,
        "source_hashes": {
            path.name: file_hash(path)
            for path in Path(__file__).resolve().parent.glob("*.py")
        },
    }
    write_json(output, certificate)
    print(json.dumps(certificate, indent=2))
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).with_name("validation.json")
    )
    arguments = parser.parse_args()
    raise SystemExit(main(arguments.output))
