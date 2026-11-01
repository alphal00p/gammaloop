"""Background jobs of the gallery: launch, observe and cancel one heavy run.

A job is a directory ``WORK/<name>-<id>/`` holding

- ``request.json``: the command, its environment, CPUs, niceness and deadline;
- ``job.json``: the supervisor's PID, process group and ``/proc`` start time;
- ``status.json``: written by the supervisor (``benchmark_job.py``);
- ``worker.log`` / ``run.log``: supervisor and command output.

The supervisor starts in a new session, so the job is one process group that
the notebook can signal without touching anything else. It holds the work
root's lock for its whole life and hands the same descriptor to the command,
so the lock is released only once every process of the job has exited. The
command also dies with its supervisor (``PR_SET_PDEATHSIG``): no worker
outlives the job. Cancel sends SIGTERM to the group, SIGKILL after a grace
period, and checks that no member survives. Recorded start times guard every
signal against PID reuse. Jobs survive a page reload; any later session can
list them, re-attach and cancel.

Pure standard library: no marimo, Symbolica or native host import.
"""

from __future__ import annotations

import fcntl
import json
import os
import signal
import subprocess
import sys
import time
import uuid
from pathlib import Path

LOCK_NAME = "benchmark.lock"
SUPERVISOR = Path(__file__).resolve().parent / "benchmark_job.py"
FINAL_STATES = frozenset({"complete", "failed", "cancelled", "timed_out", "lost"})
PROC = Path("/proc")

# Supervisors launched by this process, reaped on observation (no zombies).
_LAUNCHED: dict[int, subprocess.Popen] = {}


class JobBusyError(RuntimeError):
    """Another job holds the work root's lock."""

    def __init__(self, active):
        self.active = active
        where = f" ({active})" if active else ""
        super().__init__(
            "A gallery job is already running"
            + where
            + "; cancel it or wait for it to complete"
        )


def _stat_fields(pid, proc=PROC):
    """Fields 3.. of /proc/<pid>/stat (state first), or None if unreadable."""
    try:
        text = (Path(proc) / str(pid) / "stat").read_text()
    except (OSError, ValueError):
        return None
    try:
        return text[text.rindex(")") + 2 :].split()
    except ValueError:
        return None


def process_start_ticks(pid, proc=PROC):
    """The process start time (stat field 22), or None if gone or a zombie."""
    fields = _stat_fields(pid, proc)
    if fields is None or len(fields) < 20 or fields[0] in ("Z", "X"):
        return None
    return int(fields[19])


def group_members(pgid, proc=PROC):
    """Live (non-zombie) processes of a process group: ``{pid: start_ticks}``."""
    members = {}
    try:
        entries = list(Path(proc).iterdir())
    except OSError:
        return members
    for entry in entries:
        if not entry.name.isdigit():
            continue
        fields = _stat_fields(entry.name, proc)
        if fields is None or len(fields) < 20 or fields[0] in ("Z", "X"):
            continue
        if int(fields[2]) == pgid:
            members[int(entry.name)] = int(fields[19])
    return members


def _write_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value))
    temporary.replace(path)


def read_json(path):
    return json.loads(Path(path).read_text())


def _reap():
    for pid, process in list(_LAUNCHED.items()):
        if process.poll() is not None:
            del _LAUNCHED[pid]


def launch(root, name, request, *, supervisor=SUPERVISOR):
    """Start one serialized job; return its directory.

    ``request`` (a dict, or a callable receiving the new job directory) holds
    at least ``command``; the supervisor applies its optional ``environment``,
    ``cpus``, ``nice`` and ``max_seconds``. Raises :class:`JobBusyError` while
    another job holds the lock.
    """
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    lock = (root / LOCK_NAME).open("a+")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        raise JobBusyError(active_job(root)) from None
    try:
        work = root / f"{name}-{uuid.uuid4().hex[:12]}"
        work.mkdir(parents=True)
        if callable(request):
            request = request(work)
        (work / "request.json").write_text(json.dumps(request))
        (work / "status.json").write_text(json.dumps({"state": "starting"}))
        with (work / "worker.log").open("w") as log:
            process = subprocess.Popen(
                [sys.executable, str(supervisor), str(work), str(lock.fileno())],
                stdin=subprocess.DEVNULL,
                stdout=log,
                stderr=subprocess.STDOUT,
                pass_fds=(lock.fileno(),),
                # One process group per job: cancel signals exactly this job.
                start_new_session=True,
            )
        pid = getattr(process, "pid", None)
        if isinstance(pid, int):
            _LAUNCHED[pid] = process
            _write_json(
                work / "job.json",
                {
                    "pid": pid,
                    "pgid": pid,
                    "start_ticks": process_start_ticks(pid),
                    "started_unix_seconds": time.time(),
                },
            )
    finally:
        # The supervisor owns the inherited descriptor until the job ends.
        lock.close()
    return str(work)


def _verified_members(job, proc=PROC):
    """Processes of this job, never a recycled PID or process group."""
    pgid = job.get("pgid")
    ticks = job.get("start_ticks")
    if not isinstance(pgid, int) or pgid <= 1 or pgid == os.getpgrp():
        return {}, False
    leader = process_start_ticks(pgid, proc)
    members = group_members(pgid, proc)
    if leader is not None:
        if ticks is not None and leader == ticks:
            return members, True
        # The PID now names another process; its group is not ours.
        return {}, False
    if ticks is None:
        return {}, False
    # The supervisor is gone; Linux keeps the group ID reserved while members
    # remain. Only processes started after the job can belong to it.
    return {pid: start for pid, start in members.items() if start >= ticks}, False


def is_running(work, proc=PROC):
    """True while the job's supervisor (verified by start time) is alive."""
    _reap()
    try:
        job = read_json(Path(work) / "job.json")
    except (OSError, ValueError):
        return False
    pid = job.get("pid")
    return (
        isinstance(pid, int)
        and job.get("start_ticks") is not None
        and process_start_ticks(pid, proc) == job["start_ticks"]
    )


def status(work, proc=PROC):
    """The supervisor's status; ``lost`` if it died without a final one."""
    work = Path(work)
    current = read_json(work / "status.json")
    if current.get("state") in FINAL_STATES or not (work / "job.json").exists():
        return current
    if is_running(work, proc):
        return current
    # It may have written its final status just now.
    current = read_json(work / "status.json")
    if current.get("state") in FINAL_STATES:
        return current
    return {
        **current,
        "state": "lost",
        "message": "The job's supervisor exited without recording a final status",
    }


def cancel(work, *, grace_seconds=5.0, proc=PROC, sleep=time.sleep):
    """Stop every process of a job; return ``{"state", "signalled", "survivors"}``.

    SIGTERM first (the supervisor records ``cancelled``), SIGKILL after
    ``grace_seconds``. Survivors are listed only if even SIGKILL failed.
    """
    work = Path(work)
    try:
        job = read_json(work / "job.json")
    except (OSError, ValueError):
        return {
            "state": status(work, proc).get("state"),
            "signalled": [],
            "survivors": [],
        }
    members, whole_group = _verified_members(job, proc)
    signalled = sorted(members)
    if members:
        _write_json(work / "cancel-requested.json", {"unix_seconds": time.time()})
    for signal_number in (signal.SIGTERM, signal.SIGKILL):
        if not members:
            break
        try:
            if whole_group:
                os.killpg(job["pgid"], signal_number)
            else:
                for pid in members:
                    os.kill(pid, signal_number)
        except ProcessLookupError:
            pass
        waited = 0.0
        limit = grace_seconds if signal_number == signal.SIGTERM else 2.0
        while True:
            _reap()
            members, whole_group = _verified_members(job, proc)
            if not members or waited >= limit:
                break
            sleep(0.1)
            waited += 0.1
    _reap()
    current = read_json(work / "status.json")
    if signalled and current.get("state") not in FINAL_STATES:
        # SIGKILL left the supervisor no chance to record the cancellation.
        current = {
            **current,
            "state": "cancelled",
            "message": "Cancelled from the notebook; stopped by signal",
        }
        _write_json(work / "status.json", current)
    report = {
        "state": current.get("state"),
        "signalled": signalled,
        "survivors": sorted(members),
    }
    if signalled:
        _write_json(work / "cancel-report.json", report)
    return report


def jobs(root):
    """Known jobs of a work root, newest first."""
    root = Path(root)
    if not root.is_dir():
        return []
    found = [path for path in root.iterdir() if (path / "status.json").is_file()]
    return [str(path) for path in sorted(found, key=_mtime, reverse=True)]


def _mtime(path):
    try:
        return (path / "request.json").stat().st_mtime
    except OSError:
        return 0.0


def job_kind(work):
    """The ``kind`` recorded in a job's request, or None."""
    try:
        return read_json(Path(work) / "request.json").get("kind")
    except (OSError, ValueError, AttributeError):
        return None


def active_job(root, proc=PROC, kind=None):
    """The running job of a work root (of one ``kind``), e.g. after a reload."""
    for work in jobs(root):
        if (kind is None or job_kind(work) == kind) and is_running(work, proc):
            return work
    return None


def tail(path, limit=4000):
    """The end of a log file, decoded leniently; empty if absent."""
    try:
        with Path(path).open("rb") as stream:
            stream.seek(max(0, os.fstat(stream.fileno()).st_size - limit))
            return stream.read().decode(errors="replace")
    except OSError:
        return ""
