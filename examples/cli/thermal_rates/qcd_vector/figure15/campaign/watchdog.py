#!/usr/bin/env python3
"""Independent monotonic deadline guard for the entire compute process group."""

import argparse
import os
import signal
import time
from pathlib import Path

from common import (
    group_members,
    process_identity,
    read_json,
    same_process,
    utc_now,
    write_json,
)


def final_report(root, reason):
    # Computation has stopped. Rendering a local JSON/CSV summary may finish after the compute deadline.
    from campaign import report

    plan, state = read_json(root / "plan.json"), read_json(root / "manifest.json")
    if not plan or not state:
        return
    campaign = read_json(root / "campaign.json", {})
    ended = (
        min(time.monotonic(), campaign.get("hard_deadline_monotonic", float("inf")))
        if reason == "hard_cutoff"
        else time.monotonic()
    )
    state.update(
        status="budget_finished" if reason == "hard_cutoff" else reason,
        stop_reason=reason,
        finished_utc=utc_now(),
        active_workers=0,
        finalized_by_independent_watchdog=True,
    )
    for job in state["jobs"]:
        allocation = job.pop("current_allocation", None)
        if allocation:
            duration = max(0, ended - allocation["reserved_monotonic"])
            job["active_wall_seconds"] = job.get("active_wall_seconds", 0) + duration
            allocation.update(
                ended_utc=utc_now(), wall_seconds=duration, ended_by_watchdog=True
            )
            job.setdefault("epochs", []).append(allocation)
        if job.get("status") in ("pending", "starting", "running", "parked"):
            job.update(
                status="budget_stopped"
                if reason == "hard_cutoff"
                else "controller_failed",
                stop_reason=reason,
                ended_utc=utc_now(),
            )
    state["active_allocations"] = {}
    report(root, plan, state)


def run(root):
    campaign = read_json(root / "campaign.json")
    if not campaign:
        raise SystemExit("Campaign registration missing")
    controller = campaign["controller_identity"]
    group = campaign["compute_process_group"]
    if os.getpgrp() == group or not same_process(controller):
        raise SystemExit(
            "Watchdog must be independent of the verified live compute process group"
        )
    if controller["process_group"] != group or controller["session"] != group:
        raise SystemExit("Controller must own its session and process group")
    registration = {
        "status": "armed",
        "identity": process_identity(os.getpid()),
        "target_group": group,
        "controller_identity": controller,
        "hard_deadline_utc": campaign["hard_deadline_utc"],
        "hard_deadline_monotonic": campaign["hard_deadline_monotonic"],
        "registered_utc": utc_now(),
    }
    write_json(root / "watchdog.json", registration)
    while True:
        remaining = campaign["hard_deadline_monotonic"] - time.monotonic()
        if remaining <= 0:
            registration.update(status="hard_cutoff", cutoff_utc=utc_now())
            try:
                os.killpg(group, signal.SIGKILL)
            except ProcessLookupError:
                pass
            # The hard stop must not depend on successful filesystem writes.
            write_json(root / "watchdog.json", registration)
            for _ in range(50):
                if not [item for item in group_members(group) if item["state"] != "Z"]:
                    break
                time.sleep(0.1)
            final_report(root, "hard_cutoff")
            return 0
        if (root / "all-compute-exited.json").exists():
            # Never trust a completion marker if a compute child remains in the group.
            members = [
                identity
                for identity in group_members(group)
                if identity["pid"] != controller["pid"] and identity["state"] != "Z"
            ]
            if not members:
                registration.update(
                    status="disarmed_after_all_compute_exited", ended_utc=utc_now()
                )
                write_json(root / "watchdog.json", registration)
                return 0
        if not same_process(controller) and not group_members(group):
            # An empty old group cannot spawn future children; avoid targeting a reused PGID.
            registration.update(status="group_already_empty", ended_utc=utc_now())
            write_json(root / "watchdog.json", registration)
            if not (root / "all-compute-exited.json").exists():
                final_report(root, "controller_failed")
            return 0
        time.sleep(min(0.2, remaining))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    arguments = parser.parse_args()
    raise SystemExit(run(arguments.root.resolve()))
