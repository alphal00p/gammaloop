= FeynKit merge readiness and performance

This record follows the reconciliation of local FeynKit, remote FeynKit and
main, then the repairs required before merge. Runtime assertions, numerical
reference points, tolerances and the performance gate remain in force. A
successful producer build/upload is distinct from a successful runtime suite.

== Reproducible baseline

The main baseline is `ff533fcd91d80577e8372c398c17a5924613e73c`.
The reconciled comparison is `d3b91d31461e50d91dce4702d68c45294faca350`.
Both CLI executables were built with the normal Nix `gammaloop` package and
uploaded to the configured cache. The candidate supplies identical cards and
inputs to both binaries. Gammabench is pinned by
`benchmarks/performance/gammabench.rev`; every gate uses the declared corpus,
five paired discovery trials and, where required, five confirmation trials.
Generation includes a fresh CLI process; runtime uses warmed fixed points.
Timing comparisons wait for this task's own builds and tests to finish.

All thirteen numerical corpora agree. Eleven cases pass their performance
gates. Confirmed regressions are:

#table(
  columns: (2fr, 1fr, 1fr, 1fr),
  [Metric], [Main median], [Reconciled median], [Paired change],
  [GL1 SingleParametric runtime], [156.7 µs], [186.7 µs], [+19.1%],
  [Thermal GL1 SingleParametric runtime], [819.7 µs], [1023.5 µs], [+24.6%],
  [Thermal GL1 generation peak RSS], [119.6 MB], [194.9 MB], [+62.9%],
)

The thermal generation-time gate passes. GL1 guard generation time remains
inconclusive under the original noise policy. Other evaluator runtime modes
pass their gates. GL15 integrand construction is 45.3% faster. A separate
five-pair topology check measures 0.892 to 0.400 seconds, 55.2% faster, with the
same 66 labelled topologies and overall factors. This inventory is not a full
numerator equivalence certificate. The GL15 integrand benchmark freezes the
original physical topology because generator numbering changed; it retains
all original numerical points and tolerances and excludes enumeration timing.

The local review artifacts are under `target/reconciliation`: sealed
`performance-readiness-final-report` snapshots, `performance-readiness-summary.json`,
`topology-readiness-final-report`, and captured evaluator programs in
`program-capture`. Their binary and input hashes distinguish these measured
results from diagnostic timing probes. Private generated states remain local.

== Pipeline comparison and demonstrated costs

Both routes use the same graph, runtime inputs, selected evaluator settings,
orientation selectors and reference points. Their numerator representations
first differ at tensor preprocessing. FeynKit uses the checked Idenso tensor
carrier and shared identity planner before explicit Spenso component execution;
main's historical route performs its older specialized algebra preprocessing.
The downstream scalar compiler receives equivalent but differently structured
expressions. Runtime parameter lists and guard counts agree.

Three causes have been investigated separately:

+ Mixed physical `gamma0` traces remained unevaluated before finite component
  execution. The shared Idenso Dirac simplifier now feeds the physical time
  basis directly to its existing trace recurrence. The rule requires explicit
  four-dimensional spin and Lorentz representations and preserves formal
  dimensions and requests to leave traces unevaluated. It avoids a cycle through
  temporal-slash normalization. This is shared tensor-library behavior, not a
  GammaLoop-specific evaluator shortcut. Commit `fd1fef9f` contains the fix and
  exact trace, dimension-boundary and idempotence regressions.
+ Symbolic admission/reduction graphs overlapped in lifetime with the component
  graph. GammaLoop now consumes the reduced carrier and drops borrowed-input
  and expression owners before constructing the next graph. Validation still
  runs at each public boundary. Isolated diagnostic RSS fell substantially;
  only the final paired comparison above establishes the remaining memory cost.
+ In GL1, equivalent expression structure changes the compiler's sharing across
  conditionals. Replaying the captured main and candidate inputs through the
  same community compiler reproduces the main operation counts and the
  candidate difference. Candidate static instructions fall from 16,375 to
  14,981, but arithmetic outside the first guard increases. A diagnostic
  interpreter with the exact production selector IDs counts more executed
  operations for each candidate orientation. This supports branch locality as
  the remaining runtime cause; synthetic diagnostic inputs are not a physics
  certificate. Merely disabling common-pair elimination or selecting Horner
  optimization did not reduce the counted arithmetic cost, and no production
  optimization setting was changed to obtain a passing result.

The original GL1 discovery showed substantially larger runtime and generation
costs before the temporal-trace and lifetime fixes. Those reports remain
retained, including censored confirmation trials, but censored measurements are
not ranked as qualified performance results.

== Current repair work

The resumed merge-readiness work incorporates the on-shell agent's separately
named shared fixes: `3d166ca1` preserves dependent binomial relations in finite
tensor certificates, and `eab3d857` batches exact scalar sums during Spenso
network parsing. The latter avoids repeated copying of a growing scalar sum;
inexact coefficient domains retain their original accumulation order. These
changes require new same-scope measurements before any performance claim.
The on-shell-specific implementation remains a separate stack.

The broader initial core run selected 1,528 cases: 1,346 passed, 164 failed and
18 timed out, with 34 skipped by the existing configuration. It ran alongside
Nix compilation, so its timeouts are not an isolated baseline comparison.
Failures are retained individually in `core-failure-report.json`; the existing
failed snapshots and numerical expectations were not accepted as replacements.
The first repair batch makes missing projector, rule-leg slot and loop-basis
metadata explicit in fully authored synthetic graphs. Existing numerators,
edge/half-edge identities, declared bases, comments and test assertions are
preserved. Fixtures requiring model-rule generation or canonical cross-section
sewing are handled at their proper import boundary instead of weakening the
strict finalized-artifact parser.

The next scalar-compiler repair collects compact factors in each finite numerator
component and its scalar definitions before registration. Tensor contraction
creates new repeated factors after the earlier source-expression collection.
An identical-input compiler replay isolates the benefit to function bodies:
collecting only the outer root changes no operations. Compact collection changes
7,853 additions and 10,491 multiplications to 5,710 and 7,951, preserving all
299 inversions and 76 calls. The diagnostic interpreter counts 3,564--4,042
executed instructions across nine orientations, versus 4,682--5,249 before the
change and 3,187--3,847 on main. These instruction counts are diagnostic, not
elapsed-time gates. Scalar powers, lazy guards and numerator products remain
factorized; no compiler optimization setting or benchmark threshold changes.
The source captures, replay programs and per-orientation counts are retained in
`target/readiness-repair`. Normal CLI measurements remain required.

The runtime import repairs also preserve canonical cut and threshold inventories
for process inference after preprocessing, use the canonical `eN` edge names for
forced-cut selection, and exclude dummy edges from the spanning-forest guide as
well as the physical graph. Synthetic graphs explicitly declare their numerator,
projector, vertex slots and loop coordinates. Synthetic scalar cross sections
with two-valent insertions use canonical artifacts built from the original
incidence, with all two-particle cuts and external scalar phases retained.
