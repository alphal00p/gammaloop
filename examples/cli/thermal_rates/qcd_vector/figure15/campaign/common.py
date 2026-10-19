"""File and process identities shared by the local scan controller and guards."""

import datetime as dt
import hashlib
import json
import os
import signal
import tempfile
from pathlib import Path


def utc_now():
    return dt.datetime.now(dt.UTC).isoformat()


def read_json(path, default=None):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError):
        return default


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w") as handle:
            json.dump(value, handle, indent=2, allow_nan=False)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def append_event(root, event, **fields):
    record = dict(utc=utc_now(), event=event, **fields)
    with (Path(root) / "events.jsonl").open("a") as handle:
        handle.write(json.dumps(record, allow_nan=False) + "\n")
    return record


def process_identity(pid):
    try:
        raw = Path(f"/proc/{int(pid)}/stat").read_text()
        fields = raw[raw.rfind(")") + 2 :].split()
        return {
            "pid": int(pid),
            "state": fields[0],
            "parent_pid": int(fields[1]),
            "process_group": int(fields[2]),
            "session": int(fields[3]),
            "start_ticks": int(fields[19]),
        }
    except (OSError, ValueError, IndexError):
        return None


def same_process(identity):
    current = process_identity(identity["pid"]) if identity else None
    return bool(current and current["start_ticks"] == identity["start_ticks"])


def signal_process(identity, number=signal.SIGTERM):
    if not same_process(identity):
        return False
    try:
        os.kill(identity["pid"], number)
        return True
    except ProcessLookupError:
        return False


def group_members(group):
    result = []
    for directory in Path("/proc").iterdir():
        if directory.name.isdigit():
            identity = process_identity(int(directory.name))
            if identity and identity["process_group"] == group:
                result.append(identity)
    return result


def set_process_affinity(identity, cpus):
    if not same_process(identity):
        raise ProcessLookupError(identity)
    # Parked sessions have no integration pool, but may retain other idle threads.
    for task in Path(f"/proc/{identity['pid']}/task").iterdir():
        try:
            os.sched_setaffinity(int(task.name), cpus)
        except ProcessLookupError:
            pass


def file_hash(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def state_identity(directory):
    root = Path(directory)
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and "logs" not in path.relative_to(root).parts:
            files.append(
                {
                    "path": str(path.relative_to(root)),
                    "bytes": path.stat().st_size,
                    "sha256": file_hash(path),
                }
            )
    digest = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
    return {"directory": str(root), "sha256": digest, "files": files}
