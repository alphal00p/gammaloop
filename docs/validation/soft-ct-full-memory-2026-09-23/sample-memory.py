from pathlib import Path
from collections import deque
import ctypes, datetime, json, os, resource, shlex, signal, subprocess, time

root = Path("/tmp/soft-ct-memory-2026-09-23")
binary = Path("/common/dev/gammaloop/lcnbr/target/dev-optim/gammaloop")
card = root / "orientation-58.toml"
override = root / "muv.json"
# Recreate the input override in the local run directory.
root.mkdir(parents=True, exist_ok=True)
override.write_text(json.dumps({'global': {'generation': {'uv': {'renormalization_prescription': {'log_divergent': 'MUV', 'massive_power_divergent': 'MUV', 'massless_power_divergent': 'MUV', 'overrides': []}}}}}, indent=2) + "\n")
(root / "bare.json").write_text(json.dumps({"global": {"generation": {"uv": {"subtract_uv": False}}}}, indent=2) + "\n")
log_filter = "off,[{generation,summary}]=debug,[{generation,cff,profile}]=debug,[{generation,uv,numerator,!dump}]=debug"
parent_pid = os.getpid()
libc = ctypes.CDLL(None, use_errno=True)
gdb = "/nix/store/2xrxnqpkwlr4smxwhg5x34s2nbs2rpi6-gdb-17.2/bin/gdb"
results = []
for label, scheme, selected_card, wall in [
    ("all-u-gdb", "U", root / "all.toml", 600),
    ("all-h-gdb", "H", root / "all.toml", 900),
]:
    snapshots = []
    thresholds = [3, 6, 10]
    folder = root / label
    folder.mkdir(exist_ok=False)
    commands = (
        f"set global file {override}; "
        if scheme == "U"
        else (f"set global file {root / 'bare.json'}; " if scheme == "bare" else "")
    ) + f"run generate; save state --path {folder}/state"
    argv = [
        str(binary),
        str(selected_card),
        "-n",
        "-s",
        str(folder / "state"),
        "-t",
        "trace",
        "run",
        "-c",
        commands,
    ]
    env = dict(
        os.environ,
        RUST_BACKTRACE="full",
        GL_DISPLAY_FILTER="off",
        GL_LOGFILE_FILTER=log_filter,
        GL_TEST_LOG_DIR=str(folder / "test-logs"),
    )
    env.pop("GL_ALL_LOG_FILTER", None)
    (folder / "command.json").write_text(
        json.dumps(
            {
                "argv": argv,
                "shell_equivalent": shlex.join(argv),
                "cwd": str(root),
                "scheme": scheme,
                "input_card": str(selected_card),
            },
            indent=2,
        )
        + "\n"
    )
    start = time.monotonic()
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    peak_tree = 0
    peak_single = 0
    cap = None
    cap_sample = None
    last_telemetry = 0
    print(f"START {label}", flush=True)
    with (
        (folder / "stdout.log").open("w") as log,
        (folder / "resources.jsonl").open("w") as telemetry,
    ):
        proc = subprocess.Popen(
            argv,
            cwd=root,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
            preexec_fn=lambda: libc.prctl(0x59616D61, parent_pid, 0, 0, 0),
        )
        while proc.poll() is None:
            pending = deque([proc.pid])
            seen = set()
            rss_tree = 0
            sample = []
            while pending:
                pid = pending.popleft()
                if pid in seen:
                    continue
                seen.add(pid)
                try:
                    fields = {
                        line.split(":", 1)[0]: line.split(":", 1)[1].strip()
                        for line in Path(f"/proc/{pid}/status").read_text().splitlines()
                        if ":" in line
                    }
                    rss = int(fields.get("VmRSS", "0 kB").split()[0])
                    hwm = int(fields.get("VmHWM", "0 kB").split()[0])
                    rss_tree += rss
                    peak_single = max(peak_single, hwm)
                    sample.append({"pid": pid, "rss_kib": rss, "hwm_kib": hwm})
                    pending.extend(
                        int(x)
                        for x in Path(f"/proc/{pid}/task/{pid}/children")
                        .read_text()
                        .split()
                    )
                except (FileNotFoundError, ProcessLookupError):
                    pass
            elapsed = time.monotonic() - start
            peak_tree = max(peak_tree, rss_tree)
            row = {
                "elapsed_seconds": elapsed,
                "rss_tree_kib": rss_tree,
                "processes": sample,
            }
            if (
                elapsed - last_telemetry >= 1
                or rss_tree > 12 * 1024 * 1024
                or elapsed >= wall
            ):
                telemetry.write(json.dumps(row) + "\n")
                telemetry.flush()
                last_telemetry = elapsed
            if thresholds and rss_tree >= thresholds[0] * 1024 * 1024:
                threshold = thresholds.pop(0)
                snapshot = folder / f"stack-{threshold}GiB.txt"
                snapshot_start = time.monotonic()
                gdb_argv = [
                    gdb,
                    "-nx",
                    "-nh",
                    "-batch",
                    "-iex",
                    "set debuginfod enabled off",
                    "-iex",
                    "set auto-load off",
                    "-iex",
                    "set pagination off",
                    "-iex",
                    "set print frame-arguments none",
                    "-p",
                    str(proc.pid),
                    "-ex",
                    "thread apply all bt 45",
                    "-ex",
                    "detach",
                ]
                with snapshot.open("w") as out:
                    try:
                        stack = subprocess.run(
                            gdb_argv,
                            stdout=out,
                            stderr=subprocess.STDOUT,
                            env=dict(env, DEBUGINFOD_URLS=""),
                            timeout=90,
                        )
                        stack_code = stack.returncode
                    except subprocess.TimeoutExpired:
                        stack_code = "timeout"
                snapshots.append(
                    {
                        "threshold_gib": threshold,
                        "elapsed_seconds": elapsed,
                        "rss_tree_kib": rss_tree,
                        "gdb_seconds": time.monotonic() - snapshot_start,
                        "exit_code": stack_code,
                        "file": str(snapshot),
                    }
                )
                print("SNAPSHOT " + label + " " + json.dumps(snapshots[-1]), flush=True)
            if rss_tree > 12 * 1024 * 1024 or elapsed >= wall:
                cap = (
                    "rss_tree_exceeded_12_GiB"
                    if rss_tree > 12 * 1024 * 1024
                    else f"wall_exceeded_{wall}_seconds"
                )
                cap_sample = row
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(proc.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    proc.wait()
                break
            time.sleep(0.2)
        exit_code = proc.wait()
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = {
        "snapshots": snapshots,
        "instrumented": True,
        "label": label,
        "scheme": scheme,
        "wall_limit_seconds": wall,
        "rss_limit_kib": 12 * 1024 * 1024,
        "started_utc": started,
        "finished_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - start,
        "exit_code": exit_code,
        "termination_reason": cap,
        "cap_sample": cap_sample,
        "peak_sampled_tree_rss_kib": peak_tree,
        "peak_observed_process_hwm_kib": peak_single,
        "user_seconds": after.ru_utime - before.ru_utime,
        "system_seconds": after.ru_stime - before.ru_stime,
        "argv": argv,
        "cwd": str(root),
        "state_saved": (folder / "state/run.toml").exists(),
    }
    (folder / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    results.append(result)
    print(json.dumps(result), flush=True)
    if label == "single-u" and exit_code != 0 and cap is None:
        break
(root / "generation-results.json").write_text(
    json.dumps({"results": results}, indent=2) + "\n"
)
