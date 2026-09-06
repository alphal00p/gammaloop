#!/usr/bin/env python3
"""Guard one pipeline with a process-tree memory cap (macOS footprint, Linux RSS)."""

import argparse
import ctypes
import fcntl
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path


class RusageInfoV2(ctypes.Structure):
    # Darwin sys/resource.h: the versioned libproc ABI, including its UUID prefix.
    _fields_ = [("ri_uuid", ctypes.c_uint8 * 16)] + [
        (name, ctypes.c_uint64)
        for name in (
            "ri_user_time",
            "ri_system_time",
            "ri_pkg_idle_wkups",
            "ri_interrupt_wkups",
            "ri_pageins",
            "ri_wired_size",
            "ri_resident_size",
            "ri_phys_footprint",
            "ri_proc_start_abstime",
            "ri_proc_exit_abstime",
            "ri_child_user_time",
            "ri_child_system_time",
            "ri_child_pkg_idle_wkups",
            "ri_child_interrupt_wkups",
            "ri_child_pageins",
            "ri_child_elapsed_abstime",
            "ri_diskio_bytesread",
            "ri_diskio_byteswritten",
        )
    ]


def snapshot(known, groups, proc_pid_rusage):
    rows = subprocess.check_output(
        ["ps", "-axo", "pid=,ppid=,pgid=,rss="], text=True, timeout=2
    )
    processes = {}
    for row in rows.splitlines():
        pid, parent, group, rss = map(int, row.split())
        processes[pid] = (parent, group, rss * 1024)
    if not processes:
        raise ValueError("ps returned no processes")
    known.intersection_update(processes)
    groups.intersection_update(row[1] for row in processes.values())
    known.update(pid for pid, (_, group, _) in processes.items() if group in groups)
    while True:
        descendants = {
            pid for pid, (parent, _, _) in processes.items() if parent in known
        }
        if descendants <= known:
            break
        known.update(descendants)
    live = {pid: processes[pid] for pid in known if pid in processes}
    groups.update(row[1] for row in live.values())
    # Count only this tree, without adding whole-machine anonymous, wired or
    # compressed memory from unrelated applications. Darwin's footprint retains
    # its compressed/paged-out charge; Linux's ps RSS counts resident pages only.
    if proc_pid_rusage is not None:
        for pid, (parent, group, _) in list(live.items()):
            usage = RusageInfoV2()
            if proc_pid_rusage(pid, 2, ctypes.byref(usage)) != 0:  # RUSAGE_INFO_V2
                error = ctypes.get_errno()
                status = subprocess.run(
                    ["ps", "-p", str(pid), "-o", "stat="],
                    capture_output=True,
                    text=True,
                    timeout=2,
                )
                # A process may exit between enumeration and measurement.
                # Omit only confirmed exits/zombies; any live failure is fatal.
                if not status.stderr.strip() and (
                    (status.returncode == 1 and not status.stdout.strip())
                    or (
                        status.returncode == 0 and status.stdout.strip().startswith("Z")
                    )
                ):
                    del live[pid]
                    known.discard(pid)
                    continue
                raise OSError(error, f"cannot measure physical footprint for PID {pid}")
            live[pid] = (parent, group, usage.ri_phys_footprint)
    return live, sum(row[2] for row in live.values())


def stop_tree(known, groups, proc_pid_rusage):
    # Nextest can put tests in separate process groups. Track descendants as
    # well as the initial group so the entire pipeline is stopped at the cap.
    if not known and not groups:
        return
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    stopped = set()
    while True:
        try:
            live, _ = snapshot(known, groups, proc_pid_rusage)
        except (OSError, ValueError, KeyError, subprocess.SubprocessError):
            # Even failed monitoring must stop every process last observed. The
            # initial child also owns a process group because it starts a session.
            live = {pid: (0, 0, 0) for pid in known}
            break
        active_groups = groups - {os.getpgrp(), 0, 1}
        for group in active_groups:
            try:
                os.killpg(group, signal.SIGSTOP)
            except ProcessLookupError:
                pass
        for pid in live:
            try:
                os.kill(pid, signal.SIGSTOP)
            except ProcessLookupError:
                pass
        if set(live) <= stopped:
            break
        stopped.update(live)
        # Freeze every newly discovered child before rescanning. This closes
        # the spawn/setsid race between the first snapshot and termination.
    for group in groups - {os.getpgrp(), 0, 1}:
        try:
            os.killpg(group, signal.SIGKILL)
        except ProcessLookupError:
            pass
    for pid in live:
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    # Successful cleanup consumes these identities. A second call from finally
    # must not rediscover or signal their process groups after the child is reaped.
    # Leave them tracked if a signal fails with any other error, including EPERM.
    known.clear()
    groups.clear()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument(
        "--limit-gb",
        type=float,
        default=30.0,
        help="process-tree memory cap in decimal GB (macOS footprint, Linux RSS)",
    )
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or not 0 < args.limit_gb <= 30:
        parser.error("provide a command and a positive limit no greater than 30 GB")
    args.log.parent.mkdir(parents=True, exist_ok=True)
    lock = Path(__file__).resolve().with_name("heavy-run.lock").open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print(
            "RAM watchdog: another guarded pipeline is already active", file=sys.stderr
        )
        return 125
    known = set()
    groups = set()
    child = None
    proc_pid_rusage = None
    metric = "phys_footprint" if sys.platform == "darwin" else "rss"
    started = time.monotonic()
    peak_tree = 0
    with args.log.open("a", buffering=1) as log:

        def record(event, **values):
            log.write(
                json.dumps(
                    {"event": event, "elapsed_s": time.monotonic() - started, **values}
                )
                + "\n"
            )

        def interrupted(signum, _frame):
            raise InterruptedError(f"signal {signum}")

        signal.signal(signal.SIGINT, interrupted)
        signal.signal(signal.SIGTERM, interrupted)
        record(
            "start",
            command=command,
            tree_limit_bytes=int(args.limit_gb * 1e9),
            memory_metric=metric,
        )
        try:
            if sys.platform == "darwin":
                libproc = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
                proc_pid_rusage = libproc.proc_pid_rusage
                proc_pid_rusage.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
                proc_pid_rusage.restype = ctypes.c_int
                usage = RusageInfoV2()
                if proc_pid_rusage(os.getpid(), 2, ctypes.byref(usage)) != 0:
                    raise OSError(
                        ctypes.get_errno(), "physical-footprint preflight failed"
                    )
            elif not sys.platform.startswith("linux"):
                raise OSError(f"unsupported memory accounting platform: {sys.platform}")
            env = os.environ.copy()
            env.setdefault("CARGO_BUILD_JOBS", "2")
            env.setdefault("NEXTEST_TEST_THREADS", "2")
            child = subprocess.Popen(command, start_new_session=True, env=env)
            known.add(child.pid)
            groups.add(child.pid)
            last_record = -5.0
            while True:
                live, tree_bytes = snapshot(known, groups, proc_pid_rusage)
                peak_tree = max(peak_tree, tree_bytes)
                elapsed = time.monotonic() - started
                if tree_bytes >= args.limit_gb * 1e9:
                    record(
                        "memory_limit",
                        pids=list(live),
                        tree_bytes=tree_bytes,
                        peak_tree_bytes=peak_tree,
                    )
                    stop_tree(known, groups, proc_pid_rusage)
                    child.wait(timeout=5)
                    print(
                        "RAM watchdog: stopped the entire pipeline; "
                        f"process-tree {metric} {tree_bytes / 1e9:.2f}/{args.limit_gb:.2f} GB "
                        "(decimal GB); this is not a test verdict",
                        file=sys.stderr,
                    )
                    return 137
                if elapsed - last_record >= 5:
                    record(
                        "sample",
                        pids=list(live),
                        tree_bytes=tree_bytes,
                    )
                    last_record = elapsed
                result = child.poll()
                if result is not None:
                    record(
                        "complete",
                        returncode=result,
                        peak_tree_bytes=peak_tree,
                    )
                    return result if result >= 0 else 128 - result
                time.sleep(0.25)
        except (
            OSError,
            AttributeError,
            ValueError,
            KeyError,
            subprocess.SubprocessError,
            InterruptedError,
        ) as error:
            record("watchdog_error", error=str(error))
            print(
                f"RAM watchdog: stopping because monitoring failed: {error}",
                file=sys.stderr,
            )
            return 125
        finally:
            if child is not None:
                stop_tree(known, groups, proc_pid_rusage)
                child.wait(timeout=5)


if __name__ == "__main__":
    sys.exit(main())
