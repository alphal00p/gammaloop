= Soft-counterterm rebase onto main and bottleneck reruns

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

The nine soft-counterterm changes have been rebased with Jujutsu onto
`main@origin` at `4c84f3cdc68a` (23 September 2026). The local `main` bookmark
was advanced to the same revision. All ancestors of the new validation change
are free of merge conflicts. No changes were pushed.

The replacement runtime license restored application startup. The selected
U/H controls, orientation-58 U/H bottleneck cases, and full bare control now
complete generation and save evaluators. Both full direct-3D U/H cases exceed the
12 GiB diagnostic memory limit before evaluator construction. Projected-4D H
reaches its 600-second wall bound, also before evaluator construction.

Core unit-test compilation requires mechanical API migrations; permission to
edit those test call sites was requested and remains pending. No failing
assertions were weakened.

== Rebase and implementation

The preserved stack runs from change `yssunply` through `lmmvovxt`, followed by
validation change `vnpolmzv`. Its former base was `c08b0d9f`; the saved
pre-rebase operation,
original change mapping, and
rebased change mapping
record the recoverable Jujutsu history. Rebase command:

```sh
jj rebase -s yssunply -o main@origin
```

Conflict resolutions combine the soft prescription and component mass ownership
with main's branch-family Taylor expansion and retained numerator definitions.
The scalar root receives the loop measure once; numerator definitions undergo
the corresponding momentum and mass deformation without another measure.
Positive-degree H still uses `S(X) + U(X-S(X))`; degree-zero H uses U.

Final normalization now reaches scalar roots and retained numerator definitions.
The added raw 4D U and S series use main's existing factor-preserving series
API. Analytic UV preparation retains main's scalar-product protection and also
protects denominator metadata before tensor parsing. Forest-provenance JSON
uses main's decoded DOT-attribute contract in both writer and reader.

Both sets of upstream and soft-CT tensor tests were retained during conflict
resolution. Conflicting architecture additions were transferred to the canonical
Typst sources. The new validation change records the remaining production API
repairs separately from the rebased stack.

== Validation

- `cargo fmt --all -- --check`: passed.
- `cargo check --locked --offline --profile dev-optim -p gammaloop-api --bin gammaloop`:
  passed after replacing a removed DOT-export helper.
- `cargo build --locked --offline --profile dev-optim -p gammaloop-api --bin gammaloop`:
  passed, including the final factor-preserving 4D series changes.
- `cargo check --locked --offline --profile dev-optim -p gammaloop-integration-tests --test uv`:
  passed. This is compilation, not executed numerical acceptance.
- Final Clippy for `gammalooprs` and `gammaloop-api` passed with one existing
  `too_many_arguments` warning in UV profiling; see the
  lint log.
- Core unit-test compilation failed with 48 diagnostics, including cascades from
  obsolete `emr_vec_index` and `simplify_final` calls, the branch-family Taylor
  signature, and newly optional Symbolica coefficients and spinney construction.
  The compiler log
  preserves every diagnostic. Updating these calls is pending user approval.

The initial offline dependency check needed one pinned upstream dependency that
was not cached. A subsequent locked online check populated the cache; later
builds and checks ran offline. The source manifest records Rust, input hashes,
revision, settings and resource limits. The executable was freshly built with
Symbolica 3.0.0; the previous campaign used Symbolica 2.2.0. Any later comparison
therefore measures the complete main update rather than isolating one change.

The default-enabled CI toggle was left enabled. Automatic approval review
rejected disabling it as an unauthorized persistent CI change; local compilation
and rebase work proceeded without changing that setting.

== Performance campaign

The isolated harness is `/tmp/soft-ct-rebase-2026-09-23/rerun-generation.py`.
It reuses the previous cards and prescription overlay, serially running:

+ The original selected-orientation U/H controls: 180 seconds each.
+ Source orientation 58 with U/H: 60 seconds each.
+ Full-orientation bare control: 180 seconds.
+ Full-orientation ordinary U and soft H in direct 3D: 600 seconds each.
+ Full-orientation soft H through projected 4D: 600 seconds.

Every process tree has a sampled 12 GiB RSS limit and 0.2-second polling.
Rayon uses eight threads, generation requests one core, and the Rust default
thread stack is 128 MiB. Each run saves its command, resource trace, stage logs,
and any successfully generated state in a fresh directory. Existing September
14 states and reports are preserved. Time or memory stops are diagnostic bounds,
not acceptance passes.

The first U control aborted before generation. Its
startup log and
receipt
are retained; the harness stopped instead of repeating the same startup failure
for every card. The failed scratch directory was renamed to
`startup-license-failure`, leaving the intended run names available for resumption.

The CLI registers its application key unless built with
`NO_SYMBOLICA_OEM_LICENSE`. It was rebuilt with that supported setting enabled
and the user-supplied runtime `SYMBOLICA_LICENSE` loaded from the ignored local
`.envrc` (permissions 0600). No license value is stored in this report or its
artifacts, and no application licensing checks were removed. The
original binary fingerprint
identifies the startup failure; the successful runtime-license build has a
separate fingerprint.

The source/environment manifest
identifies the tested source and settings. Timing and memory improvements do
not establish numerical UV/IR acceptance.

== Completed generation results

All eight diagnostic cases ran sequentially with the same cards and
resource monitor. The five successful cases saved fresh evaluators. The three
bounded stops saved no completed evaluator. Peak RSS is the largest sampled
process high-water mark; the process-tree limit can overshoot between polls.

#table(
  columns: (2fr, 1fr, 1fr, 1.5fr),
  table.header([Case], [Wall time (s)], [Peak RSS (GiB)], [Outcome]),
  [Selected U], [3.209], [0.757], [Saved evaluator],
  [Selected H], [4.213], [0.778], [Saved evaluator],
  [Orientation 58 U], [7.420], [0.890], [Saved evaluator],
  [Orientation 58 H], [33.486], [1.482], [Saved evaluator],
  [Full bare], [1.806], [0.722], [Saved evaluator],
  [Full U, direct 3D], [132.806], [12.266], [12 GiB stop],
  [Full H, direct 3D], [448.838], [12.139], [12 GiB stop],
  [Full H, projected 4D], [600.812], [6.012], [600 s stop],

)

The measurement matrix
links the wall/RSS outcomes with completed evaluator milestones. Per-case command
receipts, resource samples, stdout, and compressed complete stage traces are in
the same evidence directory. The cards, prescription overlays and bounded
runner are preserved there as well.

Orientation 58 previously stopped at 60 seconds for both prescriptions, with
6.810 GiB (U) and 11.439 GiB (H) observed high-water. Both now finish inside that
bound, at 0.890 and 1.482 GiB. The retained evaluator representation contains 56
reachable generated numerator definitions for U and 252 for H. Current finite
preprocessing takes 3.768/20.302 seconds and evaluator construction then takes
0.980/6.113 seconds. The old fully materialized scalar byte counts are not
directly comparable with the new scalar roots plus shared definitions.

The full bare control completes in 1.806 seconds, compared with 2.407 seconds
previously; its observed RSS is higher (0.722 GiB versus about 0.079 GiB).
The small selected controls also use more memory than before. These are single
runs on a shared machine across a complete main/Symbolica update, not repeated
throughput measurements or isolation of one implementation change.

Both full direct-3D stops have no evaluator-start event. Their last recorded
milestones are shared local Taylor completion followed by the global-numerator
message. The resource limit is reached later in the interval before evaluator
construction; this trace does not locate a particular allocation or hot function.
Previously full U reached evaluator preprocessing and stopped at 12 GiB after
412.854 seconds. Earlier termination at the same memory limit is not a speedup.

Projected H also has no evaluator-start event. Its final recorded milestone is
certified exact CFF assignment selection; the long subsequent interval has no
completed stage log. Its 600-second termination is a diagnostic wall stop, not
a Symbolica panic, numerical failure, or completed correctness verdict.

The earlier matched evidence is in
#link("soft-ct-expression-growth-2026-09-14.typ")[the expression-growth report].
The rerun removes the observed orientation-58 generation bottleneck within the
chosen bounds; full-graph generation remains unresolved.

== Numerical checks of the successful H states

The original selected-orientation H state and orientation-58 H state both pass
the CLI's UV verdict at maximum degree of divergence -0.9. Each covers 28 LMBs,
196 subsets and 196 corresponding orientation entries: 392 resolved, passing
fits per state, with no missing-fit waiver. The 25 scaling points run from
`1e8` through `1e12`, seed 1337, with `--selected-limits all --per-orientation`.
The profiles take 18.046 and 59.945 seconds respectively, after generation.

The numerical summary
and command receipts
record coverage and outcomes; complete profiles are retained compressed in the
same evidence directory. These checks reload fresh states read-only and write
reports outside those states. They do not regenerate evaluators, edit tests,
or establish acceptance for the full-orientation or projected-4D cases.

All benchmark and numerical-check processes have finished. The remaining
limitations are full-graph generation within the stated bounds and the pending
unit-test API migrations. No numerical or assertion thresholds were changed.
