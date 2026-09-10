= Three-Loop aa to aa Performance Smoke

#quote(block: true)[
#strong[Status:] Archived, partially reproducible performance evidence

This record was captured on 2026-06-15 at evidence revision
`4fdbf430b29edd24b9e1292c07ab24a151427dc7`. The cards, runner, 339 ungrouped and 155 grouped DOT
files, and compact CSV summaries remain. The generated saved states, workspaces, and `/tmp` logs
named below do not. These are historical measurements, not current performance guarantees.
The retained cards and runner were later schema-migrated in revision
`21cefaeebe9de3f15be8efcf4f871aa6a1275954`; use the evidence revision above for the exact
capture inputs.

#strong[Resource warning:] The recorded whole-amplitude attempts reached about 400--570 GiB RSS.
Do not repeat those commands on a shared host without an explicit memory guard and resource
authorization.
]

This file records the first whole-amplitude three-loop `aa -> aa` smoke test using the lightweight `SingleParametric` runtime-evaluation path. Compilation, summed evaluators, and summed function maps are disabled in the card.

== Card and Artifacts

- Card: `examples/cli/aa_aa/3L/aa_aa.toml`
- DOT export: `examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL000.dot` through `GL338.dot`
- DOT count: 339 graph files
- Runtime evaluator: `SingleParametric`
- Function maps/compilation: disabled (`compile=false`, `summed=false`, `summed_function_map=false`)
- UV/threshold/tropical generation: disabled for this raw runtime smoke test
- Kinematics/helicity: kinematics A, helicity `+-+-`
- Final card integration smoke settings after this pass: `n_start=20`, `n_max=20`, `--batch-size 1`

== Diagram Generation

The card uses feyngen with

```text
generate amp a a > a a | a t t~ g ghG ghG~ QCD==4 QED==4 [{3}]
  --numerator-grouping no_grouping
  --symmetrize-initial-states=true --symmetrize-final-states=true
  -p aa_aa -i 3L --only-diagrams
```

Observed counts from `/tmp/aa_aa_3l_full_sequence_nolimit.log`:

#table(
  columns: 2,
  table.header([*Stage*], [*Graphs*]),
  [Symbolica generation], [6192],
  [After vetoed topologies], [4284],
  [After complete-graph filters], [1296],
  [After closed-fermion-chain analysis], [1296],
  [After external-state symmetrization], [339],
  [After canonization], [339],
  [After numerator-aware grouping (`no_grouping`)], [339],
)

The 339 graphs contain 189 isomorphically unique graph shapes. The sum of generated graph symmetry factors reported by feyngen is `-444`.

A previous attempt with numerator-aware scalar-rescaling grouping reached the same canonical 339-graph stage but was not kept as the default because it was much heavier. The retained 3L card uses `no_grouping` to keep this first all-graph generation pass tractable.

=== Release-Build Scalar-Rescaling Grouping Attempt

A follow-up card was added at `examples/cli/aa_aa/3L/aa_aa_grouped.toml` to isolate the one-time feyngen pass with

```text
--numerator-grouping group_identical_graphs_up_to_scalar_rescaling
```

The card saves DOTs to `examples/cli/aa_aa/3L/graphs_grouped` and would save the grouped state to `examples/cli/aa_aa/3L/gammaloop_state_grouped` if the feyngen grouping pass completed.

Command used with a freshly built release binary:

```bash
target/release/gammaloop --clean-state examples/cli/aa_aa/3L/aa_aa_grouped.toml run generate_diagrams
```

The release binary was built with:

```bash
nix develop --no-write-lock-file -c cargo build -p gammaloop-api --bin gammaloop --features ufo_support --release
```

Observed behavior from `/tmp/aa_aa_3l_grouped_release.log` and process monitoring:

#table(
  columns: 2,
  table.header([*Quantity*], [*Value*]),
  [Symbolica generation], [6192 graphs],
  [After vetoed topologies], [4284 graphs],
  [After complete-graph filters], [1296 graphs],
  [After closed-fermion-chain analysis], [1296 graphs],
  [After external-state symmetrization], [339 graphs],
  [After canonization], [339 graphs],
  [Time spent after canonization before manual stop], [about 80 min],
  [Peak observed RSS during grouping], [about 166 GiB],
  [CPU usage during grouping], [about 9-10 cores],
  [Grouped DOT files produced], [0],
  [Usable grouped saved state produced], [no],
)

The run was interrupted after about 80 minutes because it remained inside the scalar-rescaling numerator-aware grouping stage with no progress log after the 339 canonized graph line. A tiny startup state scaffold was removed afterward to avoid confusing it with a usable grouped saved state. This confirms that using the release binary does not make the scalar-rescaling grouping pass tractable for the full 3L `aa -> aa` graph set as currently configured.

== Integrand Generation

Command used for the no-cap run:

```bash
./gammaloop --clean-state -s /tmp/gammaloop_aa_aa_3l_perf_nolimit /tmp/aa_aa_3l_noauto.toml run generate_diagrams generate_integrands integrate_diagrams
```

Generation summary:

#table(
  columns: 2,
  table.header([*Quantity*], [*Value*]),
  [Wall time, `generate existing`], [about 29m 58s],
  [Peak RAM reported by GammaLoop], [9.59 GiB],
  [Generation cores], [10],
  [Per-graph timing sum], [291.11 min],
  [Mean graph generation time], [51.52 s],
  [Median graph generation time], [43.30 s],
  [Fastest graph], [119 ms],
  [Slowest graph], [3.15 min],
  [Graphs above 60 s], [112],
  [Graphs at least 120 s], [25],
  [Graphs below 1 s], [38],
)

Timing share over the 339 per-graph timings:

#table(
  columns: 3,
  table.header([*Phase*], [*Sum*], [*Share*]),
  [Expression build], [277.33 min], [95.27%],
  [Spenso], [6.32 min], [2.17%],
  [Symbolica eval], [7.46 min], [2.56%],
  [Compile], [0 ms], [0.00%],
)

Expression construction is the dominant generation cost. Spenso and Symbolica runtime evaluation are subdominant in this light `SingleParametric` setup.

== Whole-Amplitude Integration Smoke

Two all-graph integration attempts were made.

+ With a 100 GiB virtual-memory cap, integrand generation completed, but integration aborted before the first iteration with `memory allocation of 16896 bytes failed`. The process had reached the virtual-memory cap during integration setup and still held about 77 GiB RSS when stopped.
+ Without the cap, the same setup reached actual integration over 3 nested discrete grids, 339 graphs, and 9 continuous dimensions using Monte Carlo over graphs, orientations, and LMBs. Integration-time memory rose well beyond the capped limit: about 51 GiB RSS at integration setup, then above 570 GiB RSS before any 1000-sample iteration completed. The run was stopped manually after roughly 7 minutes of sampling to avoid unnecessary pressure on the shared machine.

No completed integration iteration was obtained from the 1000-sample, `--batch-size 100` smoke run, so this record contains no reliable MC uncertainty or ms/sample/core number for the full 339-graph amplitude. The capture-time performance conclusion was that all-graph raw 3L `SingleParametric` integration was memory dominated before it was statistics dominated. The retained card uses `n_start=20`, `n_max=20`, and `--batch-size 1` for a smaller whole-amplitude smoke attempt.

== Slowest Generation Graphs

#table(
  columns: 6,
  table.header([*Rank*], [*Graph*], [*Expression build*], [*Spenso*], [*Symbolica eval*], [*Total*]),
  [1], [GL296], [3m 0.9s], [3.95s], [4.08s], [3.15 min],
  [2], [GL299], [3m 2.0s], [2.86s], [2.54s], [3.12 min],
  [3], [GL191], [2m 30.2s], [2.63s], [2.88s], [2.60 min],
  [4], [GL169], [2m 28.5s], [3.19s], [3.17s], [2.58 min],
  [5], [GL179], [2m 30.1s], [2.18s], [2.49s], [2.58 min],
  [6], [GL213], [2m 24.6s], [2.25s], [2.29s], [2.49 min],
  [7], [GL192], [2m 21.0s], [2.40s], [2.49s], [2.43 min],
  [8], [GL216], [2m 17.9s], [3.15s], [2.44s], [2.39 min],
  [9], [GL219], [2m 18.6s], [2.30s], [2.24s], [2.39 min],
  [10], [GL226], [2m 17.3s], [2.56s], [2.63s], [2.37 min],
  [11], [GL215], [2m 16.7s], [2.49s], [2.46s], [2.36 min],
  [12], [GL178], [2m 16.4s], [2.41s], [2.72s], [2.36 min],
  [13], [GL170], [2m 15.1s], [2.27s], [2.78s], [2.34 min],
  [14], [GL186], [2m 14.8s], [2.69s], [2.59s], [2.33 min],
  [15], [GL220], [2m 15.2s], [2.46s], [2.33s], [2.33 min],
  [16], [GL187], [2m 14.0s], [2.75s], [2.59s], [2.32 min],
  [17], [GL214], [2m 12.3s], [2.75s], [2.31s], [2.29 min],
  [18], [GL198], [2m 10.6s], [2.63s], [2.34s], [2.26 min],
  [19], [GL197], [2m 8.2s], [2.84s], [2.26s], [2.22 min],
  [20], [GL309], [2m 6.0s], [2.78s], [2.09s], [2.18 min],
)

The complete 339-row timing matrix is preserved as
#link("../../examples/cli/aa_aa/3L/aa_aa_generation_timings.csv")[`aa_aa_generation_timings.csv`] rather than embedded in this page.

== Follow-Up

- Re-run the checked-in card with `n_start=20` and `--batch-size 1` to see whether the integration memory blow-up was primarily batch-size driven.
- Use the saved DOT files to build one graph at a time and establish per-graph runtime/memory behavior before returning to the full all-graph Monte Carlo setup.
- If the batch-size-one all-graph run still allocates hundreds of GiB before the first iteration, inspect the integration-state layout for graph/orientation/LMB Monte Carlo to identify what scales with all 339 graph evaluators.

== Grouped Rescaling Follow-Up

This follow-up used the release binary `target/release/gammaloop` built with

```bash
nix develop --no-write-lock-file -c cargo build -p gammaloop-api --bin gammaloop --features ufo_support --release
```

The grouped card is `examples/cli/aa_aa/3L/aa_aa_grouped.toml`. The default symbolic scalar-rescaling comparison remained impractical, so the successful one-time feyngen pass used the same scalar-rescaling grouping mode but with a single fully numerical numerator-comparison sample:

```text
--numerator-grouping group_identical_graphs_up_to_scalar_rescaling
--number-of-samples-for-numerator-comparisons 1
--fully-numerical-substitution-when-comparing-numerators=true
--compare-canonized-numerator=false
```

=== Grouped Diagram Generation

At capture time, the grouped feyngen pass completed and saved reusable DOT/state artifacts. Only
the DOT export is retained in this checkout:

#table(
  columns: 2,
  table.header([*Quantity*], [*Value*]),
  [Symbolica generation], [6192 graphs],
  [After vetoed topologies], [4284 graphs],
  [After complete-graph filters], [1296 graphs],
  [After external-state symmetrization], [339 graphs],
  [After canonization], [339 graphs],
  [After numerator-aware grouping], [155 graphs],
  [Isomorphically unique graphs], [163],
  [Color zeros], [58],
  [Grouped by scalar rescaling], [118],
  [Cancellations], [8],
  [Wall time for grouped feyngen stage], [about 7 min],
  [Peak observed RSS during numeric grouping], [about 37 GiB],
  [DOT export], [`examples/cli/aa_aa/3L/graphs_grouped`],
  [Saved state at capture time], [`examples/cli/aa_aa/3L/gammaloop_state_grouped` (not retained)],
)

The earlier release-build attempt with the default symbolic rescaling comparison remained useful as a negative control: it was stopped after about 80 min and about 166 GiB peak RSS without producing grouped DOTs/state.

=== Grouped Integrand Generation

Command:

```bash
target/release/gammaloop --override-state \
  -s examples/cli/aa_aa/3L/gammaloop_state_grouped \
  generate existing -p aa_aa -i 3L
```

Observed grouped integrand generation:

#table(
  columns: 2,
  table.header([*Quantity*], [*Value*]),
  [Graph groups], [155],
  [Wall time], [11m36s],
  [Peak RAM reported by GammaLoop], [9.15 GiB],
  [Generation cores], [10],
  [Saved integrand payload], [153.15 MiB],
  [Slowest generated graphs], [GL299, GL296, GL191, GL226, GL215],
)

=== All-Grouped-Amplitude Smoke

The grouped state was then integrated with kinematics A, helicity `+-+-`, `SingleParametric`, Monte Carlo over graphs/orientations/LMBs, and inverse-Jacobian LMB channel weights. The 20-sample all-graph smoke reached the integration loop over 155 graphs and 9 continuous dimensions, but it was stopped before the first iteration completed because RSS kept increasing:

#table(
  columns: 2,
  table.header([*Quantity*], [*Value*]),
  [Workspace], [`/tmp/aa_aa_3l_grouped_overall_smoke_mc_workspace`],
  [Requested samples], [20],
  [Batch size], [1],
  [RSS after setup/sampling start], [about 52 GiB],
  [RSS when stopped], [about 400 GiB],
  [Completed iterations], [0],
)

At capture time, grouping reduced the previous 339-graph all-amplitude setup, but the combined grouped evaluator still aggregated too much runtime state to be useful in the tested all-at-once `SingleParametric` workflow.

=== Per-Diagram Grouped Smoke Scan

A separate single-DOT runner was added at `examples/cli/aa_aa/3L/run_grouped_diagram_smokes.py`. It imports one grouped DOT at a time, generates the `SingleParametric` integrand, and runs a 20-sample integration with kinematics A and helicity `+-+-`. The compact per-graph table is stored in `examples/cli/aa_aa/3L/grouped_diagram_smokes_summary.csv`; phase-probe runs are stored in `examples/cli/aa_aa/3L/grouped_phase_probe_summary.csv`.

All 155 grouped diagrams completed and emitted integration results.

#table(
  columns: 5,
  table.header([*Quantity*], [*Min*], [*Median*], [*Mean*], [*Max*]),
  [End-to-end wall time per graph], [15.79 s], [61.26 s], [66.58 s], [169.54 s],
  [Peak RSS per graph], [2.88 GiB], [17.79 GiB], [18.11 GiB], [46.55 GiB],
  [Avg total runtime per sample], [0.176 ms], [0.402 ms], [0.447 ms], [1.627 ms],
  [Avg integrand runtime per sample], [0.102 ms], [0.267 ms], [0.314 ms], [1.417 ms],
  [NAN or unstable samples], [0.0%], [0.0%], [0.065%], [5.0%],
)

Top peak-RSS graphs:

#table(
  columns: 5,
  table.header([*Graph*], [*Wall s*], [*RSS GiB*], [*Avg ms/sample*], [*Unstable %*]),
  [GL309], [149.6], [46.6], [0.638], [0.0],
  [GL301], [163.7], [46.5], [1.627], [0.0],
  [GL303], [140.0], [46.5], [0.521], [0.0],
  [GL307], [135.3], [46.5], [0.490], [0.0],
  [GL244], [126.7], [35.9], [1.605], [0.0],
)

Top runtime-per-sample graphs:

#table(
  columns: 5,
  table.header([*Graph*], [*Wall s*], [*RSS GiB*], [*Avg ms/sample*], [*Unstable %*]),
  [GL301], [163.7], [46.5], [1.627], [0.0],
  [GL244], [126.7], [35.9], [1.605], [0.0],
  [GL296], [169.5], [25.3], [1.406], [5.0],
  [GL191], [137.3], [24.5], [1.278], [0.0],
  [GL000], [32.9], [9.5], [1.204], [0.0],
)

Only `GL296` and `GL299` reported nonzero instability in this 20-sample smoke scan, each with one unstable/NAN-class evaluation out of 20 samples. Both still emitted integration results.

=== Real vs Imaginary Training Probe

The user pointed out that the three-loop setup may need real-part rather than imaginary-part training. I ran separate high-stat single-graph probes using the same runner with explicit `integrator.integrated_phase` choices.

The first probes (`GL026`, `GL001`, `GL047`) were not reliable phase discriminators: at higher statistics, rare-weight variance dominated and neither real nor imaginary converged cleanly. A lower-max-weight graph, `GL105`, was pushed to 10M samples with real-part training:

#table(
  columns: 7,
  table.header([*Probe*], [*Graph*], [*Samples*], [*Re*], [*Re err*], [*Im*], [*Im err*]),
  [real-trained], [GL105], [10,000,000], [-5.31e-7], [1.19e-6], [6.84e-7], [1.37e-6],
)

This probe is compatible with zero in both components. I therefore did not restart the full per-diagram sweep with real-part training. The completed smoke scan should be read primarily as a runtime/memory/stability scan; its 20-sample central values are too noisy for a phase conclusion. Both real and imaginary columns are retained in the CSV summaries for follow-up.
