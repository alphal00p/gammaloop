# Guarding local test runs

`just test_LU_scalar_xs` runs all 166 scalar LU cross-section cases in release
mode with fail-fast, one Cargo worker, one test at a time and zero retries. A
watchdog covers compilation and every test descendant with a 30 GB process-tree
memory limit. Unrelated applications cannot trigger this limit.

To guard another command on macOS:

```sh
python3 bin/ram_watchdog.py --limit-gb 30 \
  --log target/ram_watchdog/run-unique/watchdog.jsonl -- COMMAND ARGUMENTS
```

The helper requires Python 3, `ps`, POSIX locks and signals. On macOS it uses
the native libproc `proc_pid_rusage` API with `RUSAGE_INFO_V2`; on Linux it uses
`ps` RSS in KiB. Failed measurement of a live owned process stops the command;
confirmed exits and zombies are omitted. There is no macOS RSS fallback.
The scalar recipe creates a fresh log directory automatically. For a manual
invocation, choose a fresh log path: the helper creates parents and appends JSONL
records, so reusing a filename combines multiple runs.

The limit uses decimal GB (1 GB = 10^9 bytes) and sums memory over tracked
descendants and process groups. On macOS, `ri_phys_footprint` includes the
process's charged compressed memory, including compressor-backed pages moved
to swap; a falling RSS therefore cannot hide that charge. This kernel accounting
metric is neither virtual address-space size nor a count of currently resident
pages. See Apple's [task-ledger definition](https://github.com/apple-oss-distributions/xnu/blob/main/osfmk/kern/task.c).

Linux retains an RSS-only limit: shared resident pages can be counted in more
than one process, and compressed or swapped-out pages are excluded. It does not
provide the macOS footprint guarantee. Whole-machine memory is not monitored.
Each start record names the metric; `tree_bytes` and `peak_tree_bytes` use it.

The default and maximum cap is 30 GB. The helper defaults `CARGO_BUILD_JOBS`
and `NEXTEST_TEST_THREADS` to 2 when unset; the
scalar recipe explicitly sets both to 1. Sampling occurs about every 0.25 seconds,
with ordinary log records every five seconds and observed peaks recorded at
completion or a cap stop. Allocations can briefly exceed a cap between samples.

The child starts a new session. Cleanup freezes newly discovered descendants
before killing the tracked process tree, including separate test process groups.
Normal completion also cleans up remaining descendants. SIGINT and SIGTERM to
the watchdog trigger this cleanup.

| Outcome | Exit status |
| --- | --- |
| Child exits normally | Child's status, including test failures |
| Child is killed by a signal | 128 plus the signal number |
| Memory cap reached | 137 |
| Lock contention, monitoring failure or interruption | 125 |

A resource stop or refusal is an incomplete run. The ignored lock beside the
script coordinates invocations within this checkout, including symlinked paths.
Use one watchdog around a pipeline; nesting another watchdog will contend for
that lock. Separate checkouts or copies of the script have independent locks.
