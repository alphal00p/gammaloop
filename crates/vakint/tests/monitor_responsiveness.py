"""Exercise synthetic monitor fixtures in a fresh, continuously drained real PTY.

No registry, saved data, cache, numerical evaluator or reduction is loaded. The
caller supplies the freshly built libtest ELF and its authenticated SHA256.
"""

from __future__ import annotations

import argparse
import codecs
import fcntl
import hashlib
import html
import itertools
import json
import os
import pty
import re
import select
import signal
import struct
import subprocess
import sys
import termios
import time
import unicodedata
from pathlib import Path

CASES = {
    "embedded-silence": "rustred_evaluation::monitor::runtime::responsiveness_tests::real_pty_silent_or_blocked_publisher_keeps_controls_live",
    "blocked-publisher": "rustred_evaluation::monitor::runtime::responsiveness_tests::real_pty_silent_or_blocked_publisher_keeps_controls_live",
    "viewer-silence": "rustred_evaluation::monitor::viewer::tests::real_pty_silent_worker_keeps_live_resources_and_detaches_without_control",
    "job-preparation": "rustred_evaluation::monitor::runtime::responsiveness_tests::real_pty_job_scope_covers_preparation_and_post_api_silence",
    "job-output": "rustred_evaluation::monitor::runtime::responsiveness_tests::real_pty_job_scope_covers_preparation_and_post_api_silence",
}
PHASES = {
    "embedded-silence": "synthetic_opaque_call",
    "blocked-publisher": "synthetic_opaque_call",
    "viewer-silence": "synthetic_paused_producer",
    "job-preparation": "prepared_input_loading",
    "job-output": "result_expression_serialization",
}
CSI = re.compile(r"\x1b\[([0-?]*)([ -/]*)([@-~])")


class Screen:
    """The cursor-addressed VT subset emitted by Crossterm (no external package)."""

    def __init__(self, width=240, height=80):
        self.width, self.height = width, height
        self.grid = [[" "] * width for _ in range(height)]
        self.styles = [[("#dddddd", "#000000", False)] * width for _ in range(height)]
        self.fg, self.bg, self.bold = "#dddddd", "#000000", False
        self.x = self.y = 0
        self.pending = ""
        self.decoder = codecs.getincrementaldecoder("utf-8")("replace")

    def feed(self, data):
        self.pending += self.decoder.decode(data)
        while self.pending:
            ch = self.pending[0]
            if ch == "\x1b":
                match = CSI.match(self.pending)
                if not match:
                    if self.pending.startswith("\x1b[") or len(self.pending) == 1:
                        return
                    self.pending = self.pending[2:]
                    continue
                self.pending = self.pending[match.end() :]
                args, _, op = match.groups()
                values = [int(v) if v else 0 for v in args.lstrip("?").split(";")]
                n = values[0] or 1
                if op in "Hf":
                    self.y = min((values[0] or 1) - 1, self.height - 1)
                    self.x = min(
                        (values[1] if len(values) > 1 else 1) - 1, self.width - 1
                    )
                elif op == "G":
                    self.x = min(n - 1, self.width - 1)
                elif op == "d":
                    self.y = min(n - 1, self.height - 1)
                elif op in "ABCD":
                    self.x = min(
                        max(self.x + (n if op == "C" else -n if op == "D" else 0), 0),
                        self.width - 1,
                    )
                    self.y = min(
                        max(self.y + (n if op == "B" else -n if op == "A" else 0), 0),
                        self.height - 1,
                    )
                elif op == "J" and values[0] in (2, 3):
                    self.grid = [[" "] * self.width for _ in range(self.height)]
                    self.styles = [
                        [(self.fg, self.bg, self.bold)] * self.width
                        for _ in range(self.height)
                    ]
                elif op == "K":
                    start, end = (
                        (0, self.width)
                        if values[0] == 2
                        else (0, self.x + 1)
                        if values[0] == 1
                        else (self.x, self.width)
                    )
                    self.grid[self.y][start:end] = [" "] * (end - start)
                    self.styles[self.y][start:end] = [(self.fg, self.bg, self.bold)] * (
                        end - start
                    )
                elif op == "m":
                    self.sgr(values)
                continue
            self.pending = self.pending[1:]
            if ch == "\r":
                self.x = 0
            elif ch == "\n":
                self.y = min(self.y + 1, self.height - 1)
            elif ch >= " " and not unicodedata.combining(ch):
                if self.x >= self.width:
                    self.x = 0
                    self.y = min(self.y + 1, self.height - 1)
                self.grid[self.y][self.x] = ch
                self.styles[self.y][self.x] = (self.fg, self.bg, self.bold)
                self.x += 2 if unicodedata.east_asian_width(ch) in "WF" else 1

    def text(self):
        return "\n".join("".join(row).rstrip() for row in self.grid)

    def sgr(self, values):
        palette = [
            "#000000",
            "#cd0000",
            "#00cd00",
            "#cdcd00",
            "#0000ee",
            "#cd00cd",
            "#00cdcd",
            "#e5e5e5",
            "#7f7f7f",
            "#ff0000",
            "#00ff00",
            "#ffff00",
            "#5c5cff",
            "#ff00ff",
            "#00ffff",
            "#ffffff",
        ]
        index = 0
        while index < len(values):
            value = values[index]
            if value == 0:
                self.fg, self.bg, self.bold = "#dddddd", "#000000", False
            elif value in (1, 22):
                self.bold = value == 1
            elif value in (39, 49):
                if value == 39:
                    self.fg = "#dddddd"
                else:
                    self.bg = "#000000"
            elif 30 <= value <= 37 or 90 <= value <= 97:
                self.fg = palette[value - (30 if value < 90 else 82)]
            elif 40 <= value <= 47 or 100 <= value <= 107:
                self.bg = palette[value - (40 if value < 100 else 92)]
            elif value in (38, 48) and index + 1 < len(values):
                color = None
                if values[index + 1] == 2 and index + 4 < len(values):
                    color = "#" + "".join(
                        f"{v:02x}" for v in values[index + 2 : index + 5]
                    )
                    index += 4
                elif values[index + 1] == 5 and index + 2 < len(values):
                    n = values[index + 2]
                    if n < 16:
                        color = palette[n]
                    elif n < 232:
                        n -= 16
                        levels = [0, 95, 135, 175, 215, 255]
                        color = "#" + "".join(
                            f"{levels[v]:02x}" for v in (n // 36, n // 6 % 6, n % 6)
                        )
                    else:
                        color = "#" + f"{8 + 10 * (n - 232):02x}" * 3
                    index += 2
                if color is not None:
                    if value == 38:
                        self.fg = color
                    else:
                        self.bg = color
            index += 1

    def svg(self):
        """Render only the captured cell grid and its actual ANSI SGR colors."""
        parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.width * 8}" height="{self.height * 18}" viewBox="0 0 {self.width * 8} {self.height * 18}">',
            "<metadata>Actual captured VT cells; ANSI true colors retained; default xterm black background, light foreground; monospace cell8x18.</metadata>",
            '<rect width="100%" height="100%" fill="#000000"/>',
            '<g font-family="DejaVu Sans Mono,Liberation Mono,monospace" font-size="14">',
        ]
        for y, row in enumerate(self.grid):
            for x, ch in enumerate(row):
                fg, bg, bold = self.styles[y][x]
                if bg != "#000000":
                    parts.append(
                        f'<rect x="{x * 8}" y="{y * 18}" width="8" height="18" fill="{bg}"/>'
                    )
                if ch != " ":
                    weight = ' font-weight="bold"' if bold else ""
                    parts.append(
                        f'<text x="{x * 8}" y="{y * 18 + 14}" fill="{fg}"{weight}>{html.escape(ch)}</text>'
                    )
        return "".join(parts) + "</g></svg>"


def identity(pid):
    path = Path("/proc") / str(pid)
    stat = (path / "stat").read_text()
    return {
        "pid": pid,
        "start_ticks": int(stat[stat.rfind(")") + 2 :].split()[19]),
        "uid": path.stat().st_uid,
    }


def available():
    return next(
        int(line.split()[1]) * 1024
        for line in Path("/proc/meminfo").read_text().splitlines()
        if line.startswith("MemAvailable:")
    )


def run_case(binary, scenario, output, cpus):
    case = output / scenario
    case.mkdir()
    snapshot = case / "snapshot.json"
    master, slave = pty.openpty()
    fcntl.ioctl(slave, termios.TIOCSWINSZ, struct.pack("HHHH", 80, 240, 0, 0))
    original = termios.tcgetattr(slave)
    env = os.environ.copy()
    env.update(
        TERM="xterm-256color",
        VAKINT_MONITOR_PTY_SCENARIO=scenario,
        VAKINT_MONITOR_PTY_SNAPSHOT=str(snapshot),
        OMP_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        RAYON_NUM_THREADS="1",
        TMPDIR=str(case),
    )
    env.pop("NO_COLOR", None)
    control_read = control_write = None
    if scenario.startswith("job-"):
        control_read, control_write = os.pipe()
        env["VAKINT_MONITOR_PTY_CONTROL_FD"] = str(control_write)
    command = [
        "taskset",
        "-c",
        cpus,
        str(binary),
        "--exact",
        CASES[scenario],
        "--ignored",
        "--nocapture",
        "--test-threads=1",
        "--color",
        "never",
    ]
    process = subprocess.Popen(
        [sys.executable, str(Path(__file__).resolve()), "--pty-exec", *command],
        stdin=slave,
        stdout=slave,
        stderr=slave,
        env=env,
        start_new_session=True,
        close_fds=True,
        pass_fds=(() if control_write is None else (control_write,)),
    )
    if control_write is not None:
        os.close(control_write)
    owned = identity(process.pid)
    pidfd = os.pidfd_open(process.pid)
    captured = bytearray()
    control = bytearray()
    screen = Screen()
    start = time.monotonic()
    ready = q_sent = detached = None
    actions = set()
    observations = []
    latest_publication = None
    before_detach = None
    result = {
        "scenario": scenario,
        "command": command,
        "process": owned,
        "floor_bytes": 110_000_000_000,
        "control_markers": "PTY stdout"
        if control_read is None
        else "dedicated inherited pipe",
    }
    try:
        while process.poll() is None:
            now = time.monotonic()
            if now - start > 50:
                raise RuntimeError(
                    "synthetic PTY fixture exceeded its diagnostic 50s deadline"
                )
            free = available()
            if free < 110_000_000_000:
                raise RuntimeError(
                    "owned synthetic fixture stopped below its 110GB host floor"
                )
            for descriptor in select.select(
                [master] + ([] if control_read is None else [control_read]),
                [],
                [],
                0.05,
            )[0]:
                try:
                    data = os.read(descriptor, 65536)
                except OSError:
                    data = b""
                if descriptor == master:
                    captured.extend(data)
                    screen.feed(data)
                else:
                    control.extend(data)
            messages = captured if control_read is None else control
            if ready is None and f"PTY_READY:{scenario}".encode() in messages:
                ready = time.monotonic()
            if ready is None:
                continue
            elapsed = now - ready
            for action, when, key in [
                ("help", 3, b"?"),
                ("close-help", 4, b"?"),
                ("freeze", 8, b" "),
            ]:
                if elapsed >= when and action not in actions:
                    os.write(master, key)
                    actions.add(action)
            if (
                scenario in ("embedded-silence", "job-preparation", "job-output")
                and snapshot.exists()
            ):
                document = json.loads(snapshot.read_text())
                monitor = document["monitor"]
                published = monitor["published_unix_seconds"]
                if published != latest_publication:
                    observations.append(
                        {
                            "publication": published,
                            "phase": monitor["phase"],
                            "before_q": q_sent is None,
                            "after_resume": detached is not None,
                            "event_age": monitor["event_age_seconds"],
                            "phase_elapsed": monitor["phase_elapsed_seconds"],
                            "sequence": document.get("sequence"),
                            "rss": monitor["resources"]["rss_bytes"],
                            "MemAvailable": monitor["resources"]["mem_available_bytes"],
                            "resource_age": monitor["resources"]["sample_age_seconds"],
                            "verified": monitor["resources"][
                                "process_identity_verified"
                            ],
                            "scope_kind": monitor.get("scope_kind"),
                            "api_status": monitor.get("api_status"),
                            "job_finished": monitor.get("job_finished"),
                            "status": monitor["status"],
                            "owner": monitor["request_sequence"],
                        }
                    )
                    latest_publication = published
            if elapsed >= 36 and q_sent is None:
                before_detach = screen.text()
                (case / "screen-before-detach.txt").write_text(before_detach)
                (case / "screen-before-detach.svg").write_text(screen.svg())
                assert (
                    PHASES[scenario].replace("_", " ") in before_detach
                    and "active" in before_detach
                )
                assert "display frozen" in before_detach
                assert (
                    "Process RSS" in before_detach and "MemAvailable" in before_detach
                )
                assert "indeterminate" in before_detach.lower()
                assert any(
                    float(age) >= 35
                    for age in re.findall(r"active (\d+\.\d+)s", before_detach)
                )
                if scenario == "blocked-publisher":
                    assert b"PTY_SINK_BLOCKED:blocked-publisher" in captured
                    assert b"PTY_SINK_RELEASED:blocked-publisher" not in captured
                q_sent = time.monotonic()
                os.write(master, b"q")
            if q_sent is not None and detached is None:
                if f"PTY_DETACHED:{scenario}".encode() in messages:
                    detached = time.monotonic()
                    assert detached - q_sent < 1, (
                        "detach or terminal restoration took >=1s"
                    )
                    assert termios.tcgetattr(slave) == original, (
                        "terminal state was not restored before producer resumed"
                    )
                    assert process.poll() is None, "paused caller unexpectedly finished"
                    os.write(master, b"continue\n")
                elif now - q_sent >= 1:
                    raise RuntimeError("no detach acknowledgment within 1s")
        while select.select([master], [], [], 0)[0]:
            try:
                data = os.read(master, 65536)
            except OSError:
                break
            if not data:
                break
            captured.extend(data)
        if control_read is not None:
            while select.select([control_read], [], [], 0)[0]:
                data = os.read(control_read, 65536)
                if not data:
                    break
                control.extend(data)
        messages = captured if control_read is None else control
        assert process.returncode == 0, captured.decode(errors="replace")[-5000:]
        assert ready is not None and detached is not None
        assert f"PTY_DONE:{scenario}".encode() in messages
        if scenario == "blocked-publisher":
            assert b"PTY_SINK_RELEASED:blocked-publisher" in captured
        if scenario in ("embedded-silence", "job-preparation", "job-output"):
            quiet_observations = [
                row
                for row in observations
                if row["before_q"] and row["phase"] == PHASES[scenario]
            ]
            assert len(quiet_observations) >= 30
            gaps = [
                b["publication"] - a["publication"]
                for a, b in itertools.pairwise(quiet_observations)
            ]
            assert max(gaps) < 2, "independent publication heartbeat stalled"
            assert (
                quiet_observations[-1]["event_age"] >= 35
                and quiet_observations[-1]["phase_elapsed"] >= 35
            )
            sequences = {
                row["sequence"]
                for row in quiet_observations
                if row["sequence"] is not None
            }
            assert len(sequences) == 1, "heartbeat fabricated native observer events"
            verified = [row for row in quiet_observations if row["verified"]]
            assert (
                verified
                and verified[0]["publication"] - quiet_observations[0]["publication"]
                < 2
            )
            assert all(
                row["rss"] and row["MemAvailable"] and row["resource_age"] < 1.5
                for row in verified
            )
            assert all(
                row["verified"]
                for row in quiet_observations
                if row["publication"] >= verified[0]["publication"]
            )
            if scenario.startswith("job-"):
                assert {row["owner"] for row in observations} == {
                    quiet_observations[0]["owner"]
                }
                assert all(
                    row["scope_kind"] == "job"
                    and row["status"] == "Running"
                    and row["job_finished"] is False
                    for row in quiet_observations
                )
                expected_api = (
                    "NotStarted" if scenario == "job-preparation" else "Succeeded"
                )
                assert all(
                    row["api_status"] == expected_api for row in quiet_observations
                )
                final = json.loads(snapshot.read_text())["monitor"]
                assert (
                    final["job_finished"] is True
                    and final["api_status"] == "Succeeded"
                    and final["status"] == "Succeeded"
                )
                assert (case / "completed-job.json").read_text() == "17\n"
        result.update(
            status="PASS",
            exit_code=process.returncode,
            quiet_seconds=q_sent - ready,
            detach_seconds=detached - q_sent,
            elapsed_seconds=time.monotonic() - start,
            observations=observations,
        )
    except BaseException as error:
        result.update(status="FAIL", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        try:
            if process.poll() is None:
                try:
                    still_owned = identity(process.pid) == owned
                except (FileNotFoundError, ProcessLookupError):
                    still_owned = False
                if still_owned:
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    try:
                        process.wait(timeout=1)
                    except subprocess.TimeoutExpired:
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                        process.wait()
        finally:
            try:
                result.update(
                    observations=observations,
                    quiet_phase_observations=[
                        row
                        for row in observations
                        if row["before_q"] and row["phase"] == PHASES[scenario]
                    ],
                    post_q_observations=[
                        row
                        for row in observations
                        if not row["before_q"] and not row["after_resume"]
                    ],
                    post_resume_observations=[
                        row for row in observations if row["after_resume"]
                    ],
                    quiet_seconds=(
                        q_sent - ready
                        if q_sent is not None and ready is not None
                        else None
                    ),
                    detach_seconds=(
                        detached - q_sent
                        if detached is not None and q_sent is not None
                        else None
                    ),
                )
                (case / "terminal.bin").write_bytes(captured)
                (case / "control.bin").write_bytes(control)
                (case / "result.json").write_text(json.dumps(result, indent=2) + "\n")
            finally:
                os.close(pidfd)
                os.close(master)
                os.close(slave)
                if control_read is not None:
                    os.close(control_read)
    return result


def run_blocked_stderr(binary, output, cpus):
    """Keep the actual stderr pipe full until scientific completion is reported."""
    case = output / "blocked-stderr"
    case.mkdir()
    read_fd, write_fd = os.pipe()
    os.set_blocking(write_fd, False)
    filled = 0
    while True:
        try:
            filled += os.write(write_fd, b"x" * 65536)
        except BlockingIOError:
            break
    assert filled == fcntl.fcntl(write_fd, fcntl.F_GETPIPE_SZ)
    os.set_blocking(write_fd, True)
    command = [
        "taskset",
        "-c",
        cpus,
        str(binary),
        "--exact",
        "rustred_evaluation::monitor::runtime::responsiveness_tests::blocked_stderr_and_enabled_logger_do_not_block_scientific_completion",
        "--ignored",
        "--nocapture",
        "--test-threads=1",
        "--color",
        "never",
    ]
    process = subprocess.Popen(
        command,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=write_fd,
        start_new_session=True,
    )
    os.close(write_fd)
    owned = identity(process.pid)
    stdout = bytearray()
    stderr = bytearray()
    started = time.monotonic()
    locked = completed = None
    result = {
        "command": command,
        "process": owned,
        "initial_full_pipe_bytes": filled,
        "stderr_undrained_until_scientific_completion": True,
    }
    try:
        assert process.stdout is not None and process.stdin is not None
        while completed is None:
            now = time.monotonic()
            if now - started > 8:
                raise RuntimeError("completion blocked behind an undrained stderr pipe")
            if available() < 110_000_000_000:
                raise RuntimeError("owned stderr fixture stopped below110GB host floor")
            if select.select([process.stdout], [], [], 0.05)[0]:
                chunk = os.read(process.stdout.fileno(), 65536)
                stdout.extend(chunk)
            if locked is None and b"STDERR_WRITER_LOCKED" in stdout:
                locked = now
            if b"STDERR_SCIENTIFIC_COMPLETE:17" in stdout:
                completed = now
            if process.poll() is not None:
                raise RuntimeError("stderr fixture exited before the completion marker")
        assert locked is not None and completed - locked < 1
        pending = struct.unpack(
            "I", fcntl.ioctl(read_fd, termios.FIONREAD, struct.pack("I", 0))
        )[0]
        assert pending == filled, "stderr was unexpectedly consumed before completion"
        # The logger and pending writer may run only after this deliberate drain.
        while len(stderr) < filled:
            stderr.extend(os.read(read_fd, filled - len(stderr)))
        process.stdin.write(b"continue\n")
        process.stdin.flush()
        while process.poll() is None:
            if time.monotonic() - started > 10:
                raise RuntimeError(
                    "stderr fixture did not finish after its sink drained"
                )
            for stream in select.select([read_fd, process.stdout], [], [], 0.05)[0]:
                chunk = os.read(
                    stream if isinstance(stream, int) else stream.fileno(), 65536
                )
                (stderr if isinstance(stream, int) else stdout).extend(chunk)
        for stream in [read_fd, process.stdout]:
            while select.select([stream], [], [], 0)[0]:
                chunk = os.read(
                    stream if isinstance(stream, int) else stream.fileno(), 65536
                )
                if not chunk:
                    break
                (stderr if isinstance(stream, int) else stdout).extend(chunk)
        assert process.returncode == 0, stdout.decode(errors="replace")
        assert b"STDERR_DONE" in stdout
        assert b"STDERR_LOGGER_CONFIRMED_AFTER_DRAIN" in stderr
        result.update(
            status="PASS",
            exit_code=process.returncode,
            locked_to_complete_seconds=completed - locked,
            elapsed_seconds=time.monotonic() - started,
        )
    except BaseException as error:
        result.update(status="FAIL", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        try:
            if process.poll() is None:
                try:
                    still_owned = identity(process.pid) == owned
                except (FileNotFoundError, ProcessLookupError):
                    still_owned = False
                if still_owned:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
        finally:
            (case / "stdout.bin").write_bytes(stdout)
            (case / "stderr.bin").write_bytes(stderr)
            (case / "result.json").write_text(json.dumps(result, indent=2) + "\n")
            os.close(read_fd)
            if process.stdout is not None:
                process.stdout.close()
            if process.stdin is not None:
                process.stdin.close()
    return result


def main():
    if sys.argv[1:2] == ["--pty-exec"]:
        # Popen already created a session; this fresh helper now claims its
        # terminal and execs the fixture without an unsafe post-fork callback.
        fcntl.ioctl(0, termios.TIOCSCTTY, 0)
        os.execvp(sys.argv[2], sys.argv[2:])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpus", default="200-203")
    args = parser.parse_args()
    assert args.binary.is_absolute() and args.output.is_absolute()
    assert hashlib.sha256(args.binary.read_bytes()).hexdigest() == args.binary_sha256
    args.output.mkdir()
    results = [
        run_case(args.binary, scenario, args.output, args.cpus) for scenario in CASES
    ]
    blocked_stderr = run_blocked_stderr(args.binary, args.output, args.cpus)
    (args.output / "result.json").write_text(
        json.dumps(
            {
                "schema": "vakint.synthetic-monitor-pty.v1",
                "status": "PASS",
                "binary": str(args.binary),
                "binary_sha256": args.binary_sha256,
                "cases": results,
                "blocked_stderr": blocked_stderr,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
