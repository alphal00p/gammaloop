#!/usr/bin/env python3
"""Run a command under process-tree RAM, disk, concurrency, and time guards."""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import queue
import shutil
import signal
import subprocess
import sys
import threading
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TextIO

import psutil

GIB = 1024**3
POLL_SECONDS = 0.25


@dataclass
class ProcessIdentity:
    pid: int
    create_time: float


@dataclass
class GuardRecord:
    job: str
    command: list[str]
    cwd: str
    started_at: float
    finished_at: float | None = None
    duration_seconds: float | None = None
    root_pid: int | None = None
    process_group: int | None = None
    returncode: int | None = None
    exit_reason: str = "starting"
    peak_rss_bytes: int = 0
    minimum_disk_free_bytes: int | None = None
    known_processes: list[dict[str, int | float]] | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", required=True)
    parser.add_argument("--cwd", type=Path, default=Path.cwd())
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--ram-gib", type=float, default=15.0)
    parser.add_argument("--disk-floor-gib", type=float, default=12.0)
    parser.add_argument("--disk-path", type=Path, default=Path.cwd())
    parser.add_argument("--min-available-ram-gib", type=float, default=6.0)
    parser.add_argument("--timeout-seconds", type=float)
    parser.add_argument("--deadline-epoch", type=float)
    parser.add_argument("--interrupt-grace-seconds", type=float, default=20.0)
    parser.add_argument("--terminate-grace-seconds", type=float, default=10.0)
    parser.add_argument("--slots-dir", type=Path)
    parser.add_argument("--max-heavy-jobs", type=int, default=2)
    parser.add_argument("--timeout-success", action="store_true")
    parser.add_argument("--cleanup-record", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command[:1] == ["--"]:
        args.command = args.command[1:]
    if args.cleanup_record is None and not args.command:
        parser.error("a command after '--' is required")
    for name in (
        "ram_gib",
        "disk_floor_gib",
        "min_available_ram_gib",
        "interrupt_grace_seconds",
        "terminate_grace_seconds",
    ):
        if getattr(args, name) < 0:
            parser.error(f"--{name.replace('_', '-')} must be non-negative")
    if args.max_heavy_jobs < 1:
        parser.error("--max-heavy-jobs must be at least 1")
    return args


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def process_identity(process: psutil.Process) -> ProcessIdentity | None:
    try:
        return ProcessIdentity(process.pid, process.create_time())
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return None


def matching_process(identity: ProcessIdentity) -> psutil.Process | None:
    try:
        process = psutil.Process(identity.pid)
        if abs(process.create_time() - identity.create_time) > 1.0e-3:
            return None
        return process
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return None


def tree_snapshot(root_pid: int) -> tuple[int, list[ProcessIdentity]]:
    try:
        root = psutil.Process(root_pid)
        processes = [root, *root.children(recursive=True)]
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return 0, []
    rss = 0
    identities: list[ProcessIdentity] = []
    for process in processes:
        identity = process_identity(process)
        if identity is None:
            continue
        identities.append(identity)
        try:
            rss += process.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return rss, identities


def signal_tree(
    process_group: int | None,
    known: dict[tuple[int, float], ProcessIdentity],
    requested_signal: signal.Signals,
) -> None:
    if process_group is not None:
        try:
            os.killpg(process_group, requested_signal)
        except ProcessLookupError:
            pass
        except PermissionError:
            pass
    for identity in known.values():
        process = matching_process(identity)
        if process is None:
            continue
        try:
            process.send_signal(requested_signal)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass


def wait_for_exit(
    known: dict[tuple[int, float], ProcessIdentity], timeout: float
) -> list[ProcessIdentity]:
    deadline = time.monotonic() + timeout
    survivors = list(known.values())
    while survivors and time.monotonic() < deadline:
        survivors = [identity for identity in survivors if matching_process(identity)]
        if survivors:
            time.sleep(POLL_SECONDS)
    return survivors


def stop_tree(
    process_group: int | None,
    known: dict[tuple[int, float], ProcessIdentity],
    interrupt_grace: float,
    terminate_grace: float,
) -> list[ProcessIdentity]:
    signal_tree(process_group, known, signal.SIGINT)
    survivors = wait_for_exit(known, interrupt_grace)
    if survivors:
        survivor_map = {(item.pid, item.create_time): item for item in survivors}
        signal_tree(process_group, survivor_map, signal.SIGTERM)
        survivors = wait_for_exit(survivor_map, terminate_grace)
    if survivors:
        survivor_map = {(item.pid, item.create_time): item for item in survivors}
        signal_tree(process_group, survivor_map, signal.SIGKILL)
        survivors = wait_for_exit(survivor_map, 2.0)
    return survivors


class SlotLease:
    def __init__(self, directory: Path | None, count: int):
        self.directory = directory
        self.count = count
        self.file: TextIO | None = None

    def acquire(self, deadline: float | None) -> None:
        if self.directory is None:
            return
        self.directory.mkdir(parents=True, exist_ok=True)
        while True:
            for index in range(self.count):
                candidate = (self.directory / f"slot_{index}.lock").open("a+")
                try:
                    fcntl.flock(candidate.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    candidate.close()
                    continue
                candidate.seek(0)
                candidate.truncate()
                candidate.write(f"pid={os.getpid()} acquired={time.time():.6f}\n")
                candidate.flush()
                self.file = candidate
                return
            if deadline is not None and time.time() >= deadline:
                raise TimeoutError(
                    "deadline expired while waiting for a heavy-job slot"
                )
            time.sleep(POLL_SECONDS)

    def release(self) -> None:
        if self.file is not None:
            fcntl.flock(self.file.fileno(), fcntl.LOCK_UN)
            self.file.close()
            self.file = None


def cleanup_from_record(args: argparse.Namespace) -> int:
    raw = json.loads(args.cleanup_record.read_text())
    identities = {
        (int(item["pid"]), float(item["create_time"])): ProcessIdentity(
            int(item["pid"]), float(item["create_time"])
        )
        for item in raw.get("known_processes", [])
    }
    root_pid = raw.get("root_pid")
    if root_pid is not None:
        try:
            root = psutil.Process(int(root_pid))
            identity = process_identity(root)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            identity = None
        if identity is not None and abs(identity.create_time - raw["started_at"]) < 5.0:
            identities[(identity.pid, identity.create_time)] = identity
            try:
                descendants = root.children(recursive=True)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                descendants = []
            for descendant in descendants:
                child_identity = process_identity(descendant)
                if child_identity is not None:
                    identities[(child_identity.pid, child_identity.create_time)] = (
                        child_identity
                    )
    survivors = stop_tree(
        raw.get("process_group"),
        identities,
        args.interrupt_grace_seconds,
        args.terminate_grace_seconds,
    )
    if survivors:
        print(
            "surviving matching processes: "
            + ", ".join(str(item.pid) for item in survivors),
            file=sys.stderr,
        )
        return 1
    return 0


def reader(stdout: TextIO, messages: queue.Queue[str | None]) -> None:
    try:
        for line in stdout:
            messages.put(line)
    finally:
        messages.put(None)


def drain_messages(
    messages: queue.Queue[str | None], log: TextIO, reader_done: bool
) -> bool:
    while True:
        try:
            line = messages.get_nowait()
        except queue.Empty:
            break
        if line is None:
            reader_done = True
            continue
        log.write(line)
        log.flush()
        sys.stdout.write(line)
        sys.stdout.flush()
    return reader_done


def run(args: argparse.Namespace) -> int:
    started_at = time.time()
    record = GuardRecord(
        job=args.job,
        command=args.command,
        cwd=str(args.cwd.resolve()),
        started_at=started_at,
    )
    atomic_json(args.record, asdict(record))
    absolute_deadline = args.deadline_epoch
    if args.timeout_seconds is not None:
        timeout_deadline = started_at + args.timeout_seconds
        absolute_deadline = (
            timeout_deadline
            if absolute_deadline is None
            else min(absolute_deadline, timeout_deadline)
        )

    lease = SlotLease(args.slots_dir, args.max_heavy_jobs)
    try:
        lease.acquire(absolute_deadline)
    except TimeoutError as error:
        record.exit_reason = "slot_deadline"
        record.finished_at = time.time()
        record.duration_seconds = record.finished_at - started_at
        atomic_json(args.record, asdict(record))
        print(error, file=sys.stderr)
        return 124

    process: subprocess.Popen[str] | None = None
    known: dict[tuple[int, float], ProcessIdentity] = {}
    requested_stop: str | None = None
    previous_handlers: dict[signal.Signals, object] = {}

    def request_stop(reason: str):
        def handler(_signum: int, _frame: object) -> None:
            nonlocal requested_stop
            requested_stop = requested_stop or reason

        return handler

    try:
        if psutil.virtual_memory().available < args.min_available_ram_gib * GIB:
            record.exit_reason = "insufficient_host_ram_before_launch"
            returncode = 75
            return returncode
        free = shutil.disk_usage(args.disk_path).free
        record.minimum_disk_free_bytes = free
        if free < args.disk_floor_gib * GIB:
            record.exit_reason = "disk_floor_before_launch"
            return 75

        args.log.parent.mkdir(parents=True, exist_ok=True)
        environment = os.environ.copy()
        environment.setdefault("MallocNanoZone", "1")
        with args.log.open("w", buffering=1) as log:
            process = subprocess.Popen(
                args.command,
                cwd=args.cwd,
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                start_new_session=True,
            )
            record.root_pid = process.pid
            record.process_group = process.pid
            atomic_json(args.record, asdict(record))
            assert process.stdout is not None
            messages: queue.Queue[str | None] = queue.Queue()
            reader_thread = threading.Thread(
                target=reader, args=(process.stdout, messages), daemon=True
            )
            reader_thread.start()
            for caught_signal in (signal.SIGINT, signal.SIGTERM):
                previous_handlers[caught_signal] = signal.getsignal(caught_signal)
                signal.signal(
                    caught_signal,
                    request_stop(
                        "keyboard_interrupt"
                        if caught_signal == signal.SIGINT
                        else "terminated"
                    ),
                )

            reader_done = False
            last_record_update = 0.0
            while process.poll() is None:
                reader_done = drain_messages(messages, log, reader_done)
                rss, identities = tree_snapshot(process.pid)
                for identity in identities:
                    known[(identity.pid, identity.create_time)] = identity
                record.peak_rss_bytes = max(record.peak_rss_bytes, rss)
                free = shutil.disk_usage(args.disk_path).free
                record.minimum_disk_free_bytes = min(
                    record.minimum_disk_free_bytes or free, free
                )
                now = time.time()
                if now - last_record_update >= 1.0:
                    record.known_processes = [
                        asdict(identity)
                        for identity in sorted(
                            known.values(), key=lambda item: item.pid
                        )
                    ]
                    atomic_json(args.record, asdict(record))
                    last_record_update = now
                if requested_stop is None and rss > args.ram_gib * GIB:
                    requested_stop = "ram_limit"
                if requested_stop is None and free < args.disk_floor_gib * GIB:
                    requested_stop = "disk_floor"
                if (
                    requested_stop is None
                    and psutil.virtual_memory().available
                    < args.min_available_ram_gib * GIB
                ):
                    requested_stop = "host_ram_reserve"
                if (
                    requested_stop is None
                    and absolute_deadline is not None
                    and now >= absolute_deadline
                ):
                    requested_stop = "timeout"
                if requested_stop is not None:
                    stop_tree(
                        record.process_group,
                        known,
                        args.interrupt_grace_seconds,
                        args.terminate_grace_seconds,
                    )
                    break
                time.sleep(POLL_SECONDS)

            try:
                process.wait(timeout=2.0)
            except subprocess.TimeoutExpired:
                requested_stop = requested_stop or "unresponsive"
                stop_tree(
                    record.process_group,
                    known,
                    args.interrupt_grace_seconds,
                    args.terminate_grace_seconds,
                )
                process.wait(timeout=2.0)
            reader_thread.join(timeout=1.0)
            drain_messages(messages, log, reader_done)
            record.returncode = process.returncode
            record.exit_reason = requested_stop or "completed"
    finally:
        for caught_signal, previous in previous_handlers.items():
            signal.signal(caught_signal, previous)
        if process is not None:
            _, identities = tree_snapshot(process.pid)
            for identity in identities:
                known[(identity.pid, identity.create_time)] = identity
        record.known_processes = [
            asdict(identity)
            for identity in sorted(known.values(), key=lambda item: item.pid)
        ]
        record.finished_at = time.time()
        record.duration_seconds = record.finished_at - started_at
        atomic_json(args.record, asdict(record))
        lease.release()

    if record.exit_reason == "completed":
        return record.returncode or 0
    if record.exit_reason == "timeout" and args.timeout_success:
        return 0
    if record.exit_reason in {"keyboard_interrupt", "terminated"}:
        return 130
    if record.exit_reason in {
        "disk_floor",
        "disk_floor_before_launch",
        "host_ram_reserve",
        "insufficient_host_ram_before_launch",
    }:
        return 75
    if record.exit_reason == "ram_limit":
        return 137
    return 124


def main() -> int:
    args = parse_args()
    if args.cleanup_record is not None:
        return cleanup_from_record(args)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
