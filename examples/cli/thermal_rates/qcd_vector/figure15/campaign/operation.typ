= Figure 15 campaign controller

`campaign.py plan` prepares commands and source hashes without starting work.
`campaign.py validate` checks the 35-point matrix, source states and CPU allocation.
Both require explicit `--run-root`, `--direct-state` and `--rest-state` paths.
Use `--binary` to select an immutable executable and `--repository` for its working
directory. Exact executable, state and controller hashes are stored in `plan.json`.

Only `campaign.py run` starts numerical work. Launch it through
`subprocess.Popen(..., start_new_session=True)` with stdout/stderr captured.
Its separate-session watchdog guards the complete controller process group,
including native sessions, shell gates and an optional reference helper.
The 12-hour monotonic deadline includes source hashing and setup. Graceful
interruption starts two minutes earlier; the watchdog sends `SIGKILL` before
writing its cutoff record, then regenerates the final JSON/CSV summary.

Ten allocations each contain six disjoint CPUs. Explicit `--cpus` accepts sixty
comma-separated CPU IDs. A `--reference-command` argv string runs without a shell
on the final six CPUs and reserves six workers until clean helper completion.
The helper must reap its own subprocesses before returning zero. An abnormal
exit retains the reservation; verified remaining children are terminated before
the controller records completion. `reserved_workers.json` contains
`{"workers": 6}` and may be reduced to zero when an external reservation ends.

Each point keeps one native CLI alive through its entire lifetime. Gate callbacks
park between integrations, with explicit ready/granted epochs. Quota interruption
targets only the native PID; running session allocations remain reserved until
the gate acknowledges return. Completed iterations are restored by the same
process with its existing warmed integrands. No fingerprint metadata is changed.
Queued points get broad first coverage before later continuations. A 600-second
wall quantum supplements the 100k, 250k, 500k and subsequent sample quotas.

The controller alone accepts precision. Native integration uses a fixed tight
accuracy to prevent early returns before training. The production minimum is
100k samples, the relative goal is 0.5%, and the absolute floor defaults to zero.
Final acceptance requires an independent reference and zero final nonfinite or
unstable discards. Error relative to the current estimate remains provisional.
Reference residuals never control stopping. The reference normalization is frozen
when a point passes; later reports retain that exact acceptance criterion.

`--reference-targets-json` reads an atomically replaceable array of objects with
`k0_over_T`, `k_over_T`, `component` (`T` or `L`) and `rho_nlo`. `E`, `k` and
`value` are accepted field alternatives. Values are the unit-coupling NLO
coefficient divided by $T^2$, with electromagnetic charge removed.

`manifest.json` stores statuses, current allocations and completed epochs;
`summary.json` and `summary.csv` contain physical transverse/longitudinal results.
`events.jsonl` records scheduling, signals and process identities. Every native
workspace retains raw checkpoints, convergence archives and failed-sample evidence.
The run has local completion artifacts and does not configure a chat notification.

`plot.py RUN_ROOT` renders the saved results in Figure 15 conventions. Native
markers use the author analytic LO term plus nine times the independently
computed native NLO coefficient. Smooth comparison lines come from the
published PDF's vector paths in `../reference/published-curves.json`, relative
to the renderer's directory. `--published-curves FILE` overrides that path.
When staging the renderer in `RUN_ROOT/controller`, also stage the curve JSON
in `RUN_ROOT/reference`. It writes full-range and explicitly clipped comparison
plots, isolated-NLO residuals, `figure15-comparison.csv`, and the marker/source
audit `figure15-render.json`. Numerical residuals still use computed author
targets. Install Matplotlib in the Python environment used for rendering;
the recorded plots used version 3.10.8. A separate virtual environment is
sufficient:

```sh
python3 -m venv .venv-figure15
.venv-figure15/bin/python -m pip install 'matplotlib==3.10.8'
.venv-figure15/bin/python examples/cli/thermal_rates/qcd_vector/figure15/campaign/plot.py RUN_ROOT
```

`render_on_exit.py RUN_ROOT` waits on kernel process descriptors for both the
registered controller and watchdog before rendering. It performs no integration
and records its state in `completion-renderer.json`. Run it with the same
Matplotlib-enabled interpreter; it passes that interpreter to `plot.py`.

`validation.py` performs bounded Python checks and kills a harmless sleep-process
group at a five-second forced deadline. Its certificate is `validation.json`.
Tiny native smoke tests require `--allow-test-settings`; use a fresh root with
`--test-point 1,0.5 --worker-limit 2 --workers-per-job 2 --wall-seconds 150
--graceful-margin-seconds 5 --quantum-seconds 30 --iteration-samples 200`.
These tests are distinct from the full scan and do not resume an existing run.
