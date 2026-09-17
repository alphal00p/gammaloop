# Comparing tth NNLO integration workspaces

From the `NNLO` example directory, use the local environment:

```bash
.venv/bin/python monitor_tth_nnlo_integration.py --list-available-graphs
.venv/bin/python monitor_tth_nnlo_integration.py --graph gl134
.venv/bin/python monitor_tth_nnlo_integration.py --graph gl297
.venv/bin/python monitor_tth_nnlo_integration.py --graph gl638 --live
```

The monitor discovers workspaces in the current directory, its `workspaces/`
directory and immediate child `workspaces/` directories. It also checks the
directory containing the script and its sibling `NNLO_experiment`. These are
relative layout conventions; no host-specific path is embedded in the script.
Use `--workspace-root PATH` to search only a chosen common directory.

Default comparisons use manifest-bearing directories named
`GL<number>_advanced_sampling`, `GL<number>_optimized_lmbs`, and optionally
`GL<number>_advanced_sampling_with_optimized_lmbs`, with an optional `_workspace`
suffix. Fixed-grid pilots, timing runs and timestamped backups are excluded.
Each graph's automatically selected comparison stays in one directory; the
current directory takes precedence over the script-relative fallbacks.
The selected paths are printed below the table. Explicit
`--advanced-workspace`, `--optimized-workspace`, and `--augmented-workspace`
arguments can select workspaces with other names.

The GL134 default reads the two completed adaptive 3-million-point runs in
`../NNLO_experiment/workspaces/GL134_{advanced_sampling,optimized_lmbs}`.
They both use 600 GeV and train the imaginary part, which carries this graph's
nonzero result (approximately 1.994e-9 and 1.999e-9 respectively). Both real
and imaginary components and their absolute integrals are displayed.

Normal mode uses the highest numerically numbered completed-iteration JSON.
`--live` reads the most recently updated matching integration log when it
contains a complete summary. If it does not, the column falls back to its
iteration dump. Every column reports its iteration count and data source;
rounded log values remain distinguishable from full-precision dump values.
The monitor reports completed results, not whether an integration process is
currently alive.
