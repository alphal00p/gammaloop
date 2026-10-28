"""Private background worker for the explicitly requested gallery benchmark."""

import json
import os
import subprocess
import sys
from pathlib import Path
from time import monotonic, time

from shared_gallery import Gallery


def read_native_progress(work):
    """Preview the native snapshot without treating file age as new progress."""
    path = work / "monitor.json"
    result = {"artifact": str(path)}
    try:
        with path.open("rb") as stream:
            payload = stream.read(Gallery.display_receipt_bytes + 1)
            modified = os.fstat(stream.fileno()).st_mtime
        if len(payload) > Gallery.display_receipt_bytes:
            raise ValueError("Native progress exceeds the browser display limit")
        document = json.loads(payload)
        if not isinstance(document, dict):
            raise TypeError("Native progress is not a JSON object")
        return {
            **result,
            "state": "available",
            "age_seconds": max(0, time() - modified),
            "snapshot": Gallery.display_receipt(document),
        }
    except FileNotFoundError:
        return {
            **result,
            "state": "pending",
            "message": "Waiting for the first native progress snapshot; loading may still be in progress",
        }
    except (OSError, ValueError, TypeError, RecursionError) as error:
        return {**result, "state": "unavailable", "message": str(error)}


def worker_resources(pid, proc=Path("/proc")):
    """Observe worker RAM during loading as well as native reduction phases."""
    result = {
        "worker_rss_bytes": None,
        "worker_peak_rss_bytes": None,
        "machine_available_bytes": None,
        "machine_total_bytes": None,
        "sampled_unix_seconds": time(),
        "source": "Linux /proc worker status and /proc/meminfo; unavailable fields are null",
    }
    for path, fields in (
        (
            proc / str(pid) / "status",
            {"VmRSS": "worker_rss_bytes", "VmHWM": "worker_peak_rss_bytes"},
        ),
        (
            proc / "meminfo",
            {
                "MemAvailable": "machine_available_bytes",
                "MemTotal": "machine_total_bytes",
            },
        ),
    ):
        try:
            for line in path.read_text().splitlines():
                key, _, value = line.partition(":")
                if key in fields:
                    amount, unit = value.split()
                    if unit == "kB":
                        result[fields[key]] = int(amount) * 1024
        except (OSError, ValueError):
            continue
    return result


def main():
    work = Path(sys.argv[1]).resolve()
    lock_fd = int(sys.argv[2])
    started = monotonic()

    def status(state, **fields):
        temporary = work / "status.tmp"
        temporary.write_text(
            json.dumps(
                {
                    "state": state,
                    "elapsed_seconds": monotonic() - started,
                    "updated_unix_seconds": time(),
                    **fields,
                }
            )
        )
        temporary.replace(work / "status.json")

    def receipt():
        path = work / "benchmark.json"
        if not path.exists():
            return None
        if path.stat().st_size > Gallery.display_receipt_bytes:
            return {
                "artifact": str(path),
                "bytes": path.stat().st_size,
                "display": "Full receipt exceeds the browser display limit; scientific data remains in the artifact",
            }
        return Gallery.display_receipt(json.loads(path.read_text()))

    def log_tail():
        with (work / "run.log").open("rb") as stream:
            stream.seek(max(0, os.fstat(stream.fileno()).st_size - 5000))
            return stream.read().decode(errors="replace")

    try:
        command = json.loads((work / "request.json").read_text())["command"]
        status("running", native_progress=read_native_progress(work))
        with (
            (work / "run.log").open("w") as log,
            subprocess.Popen(
                command,
                stdout=log,
                stderr=subprocess.STDOUT,
            ) as process,
        ):
            while True:
                observed = {
                    "worker_pid": process.pid,
                    "resources": worker_resources(process.pid),
                }
                status(
                    "running",
                    **observed,
                    native_progress=read_native_progress(work),
                )
                try:
                    # This wait bounds status latency, never the job's duration.
                    exit_code = process.wait(timeout=5)
                    break
                except subprocess.TimeoutExpired:
                    continue
        if exit_code:
            status(
                "failed",
                **observed,
                native_progress=read_native_progress(work),
                exit_code=exit_code,
                message=log_tail(),
                receipt=receipt(),
            )
        else:
            status(
                "complete",
                **observed,
                native_progress=read_native_progress(work),
                receipt=receipt(),
            )
    except (OSError, ValueError) as error:
        status("failed", message=str(error))
    finally:
        os.close(lock_fd)


if __name__ == "__main__":
    main()
