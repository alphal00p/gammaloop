= GL638 X1 auxiliary-radial pilot
<gl638-x1-auxiliary-radial-pilot>
This pilot does #strong[not] establish a reliable efficiency improvement. LU-h reduces the largest excursion seen with the matched ordinary-cut control, but that comparison is dominated by one seed; the other two seeds have larger LU-h errors. The remaining shape singularities are unresolved, and these results do not establish convergence, finite variance or bounded weights.

The frozen source is `3f9b1a28473f313d6242a556ada2e8dbcd6f6477`, built with `--profile dev-optim --locked -p gammaloop-api --bin gammaloop --features ufo_support`. Binary SHA256: `f700ee575153067ccc77894ea776a96f139105812a15014d094dd67ba7db6d61`. All controls and candidates use this binary and the fresh #link("GL638_CHECKPOINT06D_BASELINE.json")[checkpoint state];. All #strong[936 production orientations] are summed, with the original threshold metadata, full 3D local UV and integrated UV, and the complete physical cut/CT sum.

Four distinct proposals each use seeds 1337, 7331 and 424242, with 2,048 draws in one untrained iteration on 20 cores: #strong[24,576 draws] total. All twelve runs completed with no final NaN/unstable sample, no Arb rescue, and unchanged saved-state hashes. The whole batch, including one state load, five full-state point checks, settings changes and worker setup, took 586.75 seconds. The explicit six-LMB and `auto:optimized_lmb` catalogues resolve to the same basis IDs `[1,2,13,14,21,22]`; their point results agree exactly, so only one is integrated.

#figure(
  align(center)[#table(
    columns: 6,
    align: (auto,auto,right,right,right,right,),
    table.header([Proposal], [Phase], [Pooled signed estimate ± SE \[pb\]], [Pooled absolute estimate ± SE \[pb\]], [Largest absolute weight], [Mean recorded sample time \[ms\]],),
    table.hline(),
    [Six optimized LMBs], [Re], [4.1816e-4 ± 2.66e-4], [7.2896e-4 ± 2.66e-4], [1.455], [78.71],
    [Six optimized LMBs], [Im], [1.4520e-4 ± 3.24e-4], [1.1993e-3 ± 3.24e-4], [1.245], [78.71],
    [Six cut charts, ordinary power 1], [Re], [1.7420e-5 ± 2.33e-5], [5.0191e-5 ± 2.33e-5], [0.1334], [79.01],
    [Six cut charts, ordinary power 1], [Im], [-3.4981e-4 ± 3.17e-4], [3.7984e-4 ± 3.17e-4], [1.937], [79.01],
    [Six cut charts, LU-h], [Re], [-9.4737e-7 ± 1.99e-5], [5.8255e-5 ± 1.98e-5], [0.07951], [86.13],
    [Six cut charts, LU-h], [Im], [-1.7649e-4 ± 1.92e-4], [2.8316e-4 ± 1.92e-4], [1.153], [86.13],
    [LU-h plus optimized LMBs], [Re], [4.6928e-4 ± 4.76e-4], [7.6544e-4 ± 4.76e-4], [2.901], [76.94],
    [LU-h plus optimized LMBs], [Im], [5.4156e-4 ± 4.49e-4], [1.1794e-3 ± 4.49e-4], [2.483], [76.94],
  )]
  , kind: table
  )

The absolute monitors are `|Re I|` and `|Im I|` separately. Equal-size independent runs are pooled using within-run sums of squares inferred from the reported errors, plus between-seed mean differences. These empirical standard errors are not a finite-variance certificate. Per-sample time is the recorded evaluation statistic, not end-to-end wall time. The seeds are paired across proposals; no significance estimate neglecting their covariance is claimed.

The LU-h/cut-control imaginary error comparison is `5.669e-4 / 9.504e-4` for seed 1337, `9.997e-5 / 2.076e-5` for seed 7331, and `6.798e-6 / 3.901e-6` for seed 424242. Thus the smaller pooled error is not a consistent three-seed gain. Large seed dependence and the much smaller cut-only absolute estimates also prevent interpreting those estimates as converged. Surface-focused angular sampling and LMB products visit very different shape regions at this sample count.

All cut-only draws passed in Double precision. The LMB control used 0.83--1.46% Quad and the mixture 0.44--0.73% Quad. Some lower-precision attempts printed `process evaluation is nonfinite`; the final stability-stack results contain no invalid sample.

The six physical cuts are:

#figure(
  align(center)[#table(
    columns: 2,
    align: (right,auto,),
    table.header([CutId], [Energy edges],),
    table.hline(),
    [0], [`[2,6,12,13]`],
    [1], [`[2,6,10]`],
    [2], [`[2,6,7,13,14]`],
    [3], [`[2,4,12]`],
    [4], [`[2,4,10,13]`],
    [5], [`[2,4,7,14]`],
  )]
  , kind: table
  )

Each named cut chart uses `phase_space(cut(...))` with parent/subspace `[3,4,7,10]`, its actual numeric `on_cut` ID, and explicit `map_density`. The LU-h profile is the log-logistic approximation to physical `poly_exponential` h with sigma 1 and power 3, plus a 2% ordinary raw radial floor with beta 300 GeV. The comparison cut chart has the same 12-dimensional angular map and ordinary radial power 1. This matched control separates the radial change from the difference between a hyperspherical cut chart and products of 3D LMB spheres. Threshold localization settings are unchanged.

Two setup mistakes were identified and corrected before integration. The reduced o0 smoke only activates CutIds `[0,1,3]`, so its smoke-only selection was narrowed while all-orientation cards retained all six cuts. Also, `set process ... file` merges nested settings: omitting `radial_profile` under a reused name retained its prior value. Distinct `radial_cut_N` and `lu_cut_N` names prevent that carryover. The actual saved settings verify the selected ordinary control has no radial profile. Both earlier smoke artifacts are retained as provenance.

The portable #link("GL638_X1_PILOT.json")[JSON artifact] embeds all five runtime cards, resolved per-mode settings, source/binary/state/input hashes, original CLI argv, all 24 per-seed component results, precision statistics, and the largest-weight cube/discrete point for each run and component (24 records). Scratch paths are provenance; extract the cards and substitute fresh binary/state/workspace paths to reproduce. Maximum replay must account for the original discrete/grid probability: fixed-index `inspect` does not itself reproduce that outer integration weight. No maximum replay is included in this pilot.

The next GL638 step remains the generic conditional host/threshold chart, with normalized-reference, inverse-density, cut-preservation and native-precision gates before a direct-H pilot. Longer adaptive or equal-wall-time comparisons can then assess a candidate that addresses the shape singularity; this auxiliary-radial pilot does not replace those checks.
