#!/usr/bin/env python3
"""Park a persistent CLI between integrations until its matching epoch is granted."""

import argparse
import os
import signal
import time
from pathlib import Path

from common import (
    read_json,
    same_process,
    signal_process,
    utc_now,
    write_json,
)


def run(root, name, epoch):
    job = root / "jobs" / name
    registration = None
    for _ in range(100):
        registration = read_json(job / "session.json")
        if registration:
            break
        time.sleep(0.05)
    if not registration or not same_process(registration["identity"]):
        return 2
    identity = registration["identity"]
    write_json(
        job / "gate.json",
        {
            "epoch": epoch,
            "state": "ready",
            "gate_pid": os.getpid(),
            "native_pid": identity["pid"],
            "utc": utc_now(),
        },
    )
    while same_process(identity):
        campaign = read_json(root / "campaign.json", {})
        if (
            (job / "stop.json").exists()
            or (root / "stop.json").exists()
            or (
                campaign.get("soft_deadline_monotonic") is not None
                and time.monotonic() >= campaign["soft_deadline_monotonic"]
            )
        ):
            write_json(
                job / "gate.json",
                {
                    "epoch": epoch,
                    "state": "stopping",
                    "gate_pid": os.getpid(),
                    "native_pid": identity["pid"],
                    "utc": utc_now(),
                },
            )
            # The shell's immediate parent is not necessarily the CLI; use its registered identity.
            signal_process(identity, signal.SIGTERM)
            return 0
        grant = read_json(job / "grant.json", {})
        if (
            grant.get("epoch") == epoch
            and grant.get("native_start_ticks") == identity["start_ticks"]
        ):
            write_json(
                job / "gate.json",
                {
                    "epoch": epoch,
                    "state": "granted",
                    "gate_pid": os.getpid(),
                    "native_pid": identity["pid"],
                    "utc": utc_now(),
                    "quota": grant.get("quota"),
                },
            )
            return 0
        time.sleep(0.1)
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("name")
    parser.add_argument("epoch", type=int)
    arguments = parser.parse_args()
    raise SystemExit(run(arguments.root.resolve(), arguments.name, arguments.epoch))
