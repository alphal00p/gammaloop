= Matched ordinary-U comparison against clean main

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

Clean main also exceeds the 12-GiB bound on the exact full-orientation ordinary-U
case. The previous investigation identified where the rebased stack allocates
memory, but did not establish a regression against main. This matched control
closes that gap: main stops after 122.143 seconds, versus 132.806 seconds for
ordinary U on the stack. These are times to a resource stop, not completed
generation times or a statistically established speed ratio.

== Inputs and measurements

Main is `4c84f3cdc68ac184f986c1cb2bc1e00446882950`, built unmodified in an
isolated Jujutsu workspace. The stack executable is the preserved build of
`2c2694036e93`. Both use the same locked dependencies, Symbolica 3.0.0,
`dev-optim` profile, cached Rust toolchain, and runtime user-license mode.
The working workspace's stack executable was restored after building main.

Both import the same absolute DOT fixture,
`paper_figure_b1_double_triangle_soft_ir`, with its physical numerator and
projector. The original command cards and MUV overlay are reused: direct 3D,
hedge-poset orchestration, ordinary U everywhere, integrated counterterms off,
thresholds off, evaluator compilation off, and one graph-generation worker.
Persisted selected-run global settings match except for isolated state paths;
runtime settings match exactly. Full runs request all explicit orientations.

#table(
  columns: (auto, auto, auto, auto, auto),
  [Case], [Main seconds], [Stack seconds], [Main/stack GiB], [Outcome],
  [Selected pattern], [4.413], [3.209], [0.758 / 0.757], [Both completed],
  [Orientation 58], [7.621], [7.420], [0.881 / 0.890], [Both completed],
  [Full ordinary U], [122.143], [132.806], [12.071 / 12.266], [Both hit memory bound],
)

Memory columns are observed process RSS high-water marks. The monitor enforces
a 12-GiB process-tree RSS bound with 0.2-second polling, so sampled peaks can
overshoot. Full generation has a 600-second wall bound. Main measurements are
new; stack timings come from the earlier uninstrumented rerun. The machine is
shared, and these are single measurements. No claim of a precise percentage
regression follows from their timing difference. Debugger runs are separate
and excluded from this table.

== Complete pipeline comparison

+ #strong[Shared input and raw CFF.] Same graph and settings; 102 native source
  maps on both revisions.
+ #strong[First possible differences: forest classification and preparation.]
  The stack classifies disconnected UV components separately and reports
  routing failures. Main additionally builds unused local-4D atoms even with
  integrated counterterms disabled; the direct stack skips that work.
+ #strong[Different Taylor coordinates and mass bookkeeping.]
  The stack chooses graph-canonical component coordinates, retains nested
  carriers and scopes physical/UV masses. These changes also affect ordinary U;
  selecting U does not reproduce main's implementation.
+ #strong[Shared Taylor-family engine.]
  `direct_3d/branches.rs` is unchanged. Both ordinary-U runs log 20 shared
  coefficient-family expansions, with maxima of 6 source families, 102 residue
  rows and 8 coefficient families per event. Maximum retained source-body
  sizes are 9,526 bytes on main and 9,931 on the stack; coefficient-body maxima
  are 17,319 and 18,599 bytes. These are per-event body totals, not complete
  integrand sizes or total live memory.
+ #strong[Different finalization order.]
  Main materializes the residue sum before metric/color/dot normalization.
  The stack normalizes branch roots and retained definitions before
  materialization. The normalization operations themselves are shared apart
  from removal of the stack's temporary mass-support tags. This is a possible
  representation difference, not an established cause of a regression.
+ #strong[Shared final sum and factor collection.]
  `Integrands::zip_add`, opaque-factor protection, the fixed-point
  `collect_factors()` call and Symbolica's collector are unchanged. Neither
  full-U run saves an evaluator state.

There is no newly introduced global `.expand()` in the ordinary-U production
path. Native Taylor coefficient extraction, color/tensor collection and finite
component contractions still perform algebra; absence of that method call is
not a claim of zero distribution. The detailed expansion-boundary audit remains
in #link("soft-ct-full-memory-2026-09-23.typ")[the memory investigation].

== Collector dimensions

Separate bounded debugger snapshots of clean main and the stack both stop at
Symbolica `collect.rs:1343`. Both read exactly 1,030 summand maps and 7,112
distinct factor keys, with current common exponent zero. Equal dimensions do
not establish symbolic equality of the two expressions, but they directly
identify the same allocation mechanism and the same map-entry count.
The collector inserts missing keys into each summand map even at zero exponent.
Completing this phase therefore entails 7,325,360 map entries for this one sum,
before accounting for owned-key copies and other live structures. This does
not require a globally expanded physical numerator.

== Scope of the conclusion

This exact full-U failure is already present on clean main. Soft H adds real
Taylor work: the stack logs 33 family expansions instead of 20, and the earlier
orientation-58 H control takes 33.486 seconds versus 7.420 seconds for U. That
additional work is separate from the memory cliff already reproduced without
soft subtraction. The full-H 448.838-second result is also censored by the
memory bound and cannot be interpreted as a completed-generation slowdown.

The comparison covers this graph and these settings. A different successful
main command card requires its own matched comparison. No production code,
tests, algebraic settings or CI settings were changed for this investigation.

== Evidence

The adjacent `soft-ct-main-comparison-2026-09-23/` directory retains input cards,
the MUV overlay, command receipts, result/resource records, compressed stage
logs, debugger snapshots and scripts, executable fingerprints, the build log,
the settings comparison and a source/method manifest. Full scratch states and
binaries remain under `/tmp/soft-ct-main-comparison-2026-09-23/`.
