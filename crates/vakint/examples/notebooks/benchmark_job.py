"""Supervisor of one gallery job: run its command, publish status, stop on request.

Started by ``gallery_jobs.launch`` in a new session as
``benchmark_job.py WORK LOCKFD``. It runs ``request.json``'s command with the
requested environment, CPUs and niceness, hands it the work root's lock and
makes it die with this supervisor. It writes ``status.json`` every five
seconds. SIGTERM/SIGINT/SIGHUP (the notebook's Cancel) stop the command and
record ``cancelled``; an optional ``max_seconds`` deadline records
``timed_out``. A command still alive five seconds after a stop is killed.
"""

import json
import os
import signal
import subprocess
import sys
from pathlib import Path
from time import monotonic, time

from shared_gallery import Gallery

STOP_SIGNALS = (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)
KILL_AFTER_SECONDS = 5


def _parent_death_signal():
    """Return a fork-safe callable that makes the caller die with its parent."""
    if not sys.platform.startswith("linux"):
        return lambda: None
    import ctypes

    prctl = ctypes.CDLL(None, use_errno=True).prctl

    def arm():
        # PR_SET_PDEATHSIG = 1
        if prctl(1, int(signal.SIGKILL), 0, 0, 0) != 0:
            raise OSError(ctypes.get_errno(), "prctl(PR_SET_PDEATHSIG) failed")

    return arm


def child_setup(request, parent):
    """Pre-exec hook: die with the supervisor, then apply niceness and CPUs."""
    arm = _parent_death_signal()
    cpus = request.get("cpus")
    niceness = request.get("nice")

    def setup():
        arm()
        if os.getppid() != parent:
            # The supervisor died before the death signal was armed.
            os._exit(70)
        if niceness:
            os.nice(int(niceness))
        if cpus:
            os.sched_setaffinity(0, [int(cpu) for cpu in cpus])

    return setup


def child_environment(request):
    extra = request.get("environment")
    if not extra:
        return None
    return {**os.environ, **{str(k): str(v) for k, v in extra.items()}}


def read_stage(work):
    """The command's own coarse stage marker (``stage.json``), if any."""
    try:
        stage = json.loads((work / "stage.json").read_text())
    except (OSError, ValueError):
        return None
    return stage if isinstance(stage, dict) else None


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
    stop = {"reason": None, "at": None}
    running = {}

    def status(state, **fields):
        temporary = work / "status.tmp"
        temporary.write_text(
            json.dumps(
                {
                    "state": state,
                    "elapsed_seconds": monotonic() - started,
                    "updated_unix_seconds": time(),
                    "stage": read_stage(work),
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

    def request_stop(reason):
        if stop["reason"] is None:
            stop["reason"] = reason
            stop["at"] = monotonic()
        process = running.get("process")
        if process is not None:
            try:
                process.send_signal(signal.SIGTERM)
            except OSError:
                pass

    def on_signal(_signal, _frame):
        request_stop("cancelled")

    previous = {}
    try:
        for number in STOP_SIGNALS:
            previous[number] = signal.signal(number, on_signal)
    except ValueError:
        # Not the main thread (embedded use): rely on the deadline only.
        pass
    try:
        request = json.loads((work / "request.json").read_text())
        command = request["command"]
        max_seconds = request.get("max_seconds")
        status("running", native_progress=read_native_progress(work))
        with (
            (work / "run.log").open("w") as log,
            subprocess.Popen(
                command,
                stdin=subprocess.DEVNULL,
                stdout=log,
                stderr=subprocess.STDOUT,
                # The command keeps the lock until it exits, even if this
                # supervisor is killed first.
                pass_fds=(lock_fd,),
                env=child_environment(request),
                # preexec_fn is required for PR_SET_PDEATHSIG, which no other
                # Popen option provides; the supervisor runs no other threads.
                preexec_fn=child_setup(request, os.getpid()),  # noqa: PLW1509
            ) as process,
        ):
            running["process"] = process
            if stop["reason"] is not None:
                request_stop(stop["reason"])
            while True:
                observed = {
                    "worker_pid": process.pid,
                    "resources": worker_resources(process.pid),
                }
                status(
                    "stopping" if stop["reason"] else "running",
                    **observed,
                    native_progress=read_native_progress(work),
                )
                try:
                    # This wait bounds status latency, never the job's duration.
                    exit_code = process.wait(timeout=5)
                    break
                except subprocess.TimeoutExpired:
                    if stop["reason"] is not None:
                        if monotonic() - stop["at"] >= KILL_AFTER_SECONDS:
                            process.kill()
                    elif max_seconds is not None and monotonic() - started >= float(
                        max_seconds
                    ):
                        request_stop("timed_out")
                    continue
        final = {
            **observed,
            "native_progress": read_native_progress(work),
            "exit_code": exit_code,
        }
        if stop["reason"] == "cancelled":
            status("cancelled", **final, message="Cancelled on request")
        elif stop["reason"] == "timed_out":
            status(
                "timed_out",
                **final,
                message=f"Stopped at the {float(max_seconds):g} s deadline",
            )
        elif exit_code:
            status("failed", **final, message=log_tail(), receipt=receipt())
        else:
            status("complete", **final, receipt=receipt())
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        status("failed", message=str(error))
    finally:
        os.close(lock_fd)
        for number, handler in previous.items():
            signal.signal(number, handler)


if __name__ == "__main__":
    main()
