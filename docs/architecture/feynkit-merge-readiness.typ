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
Counting each n-ary addition/multiplication by its number of binary operations
and every other executed instruction once gives ranges of 3,556--4,344 on main,
5,072--5,814 before component collection, and 3,789--4,392 afterward. This second
diagnostic records operand cost as well as control flow; it is not a cycle model.
The source captures, replay programs and per-orientation counts are retained in
`target/readiness-repair`. Normal CLI measurements remain required.

The runtime import repairs infer processed cross sections from their cut owner,
retain the existing ownership transfer, use the canonical `eN` edge names for
forced-cut selection, and exclude dummy edges from the spanning-forest guide as
well as the physical graph. Synthetic graphs explicitly declare their numerator,
projector, vertex slots and loop coordinates. Synthetic scalar cross sections
with two-valent insertions use canonical artifacts built from the original
incidence, with all two-particle cuts and external scalar phases retained.

A Valgrind Massif diagnostic of thermal GL1 locates 86.7 MB of the 126.8 MB
peak useful heap in the initial Idenso candidate/region observation walk, plus
8.8 MB in its leaf hash sets. Identical subtrees appeared in several branches;
the completed scope map reused them during later replacements but the initial
walk retained separate observation trees through each parent's child list.
The repair interns completed observations in that same existing scope map
before attaching them to their parent. Every occurrence is still scanned and
counted; only immutable observation storage is shared. The nested-inventory
regression now requires one retained copy per distinct function, a stricter
bound than its former one-copy-per-occurrence expectation. No algebra budget,
identity rule, or validation check is skipped. The raw heap profile and decoded
peak tree are retained under `target/readiness-repair/massif-thermal-candidate`;
new runtime and RSS gates are required before claiming this fixes the regression.

The frozen two-loop ghost input also retained the pre-correction vertex signs.
The original main test input was imported and exported with the pinned main
CLI. Across the six diagrams used by this test, all five non-ghost vertices
agree exactly; all nineteen V_35 numerators are exact negatives of the frozen
ones. Only those signs are refreshed from main. Topology, momentum bases,
projectors and propagator phases remain unchanged, as do all independent RQFT
Laurent expectations. The source card, baseline output and vertex-by-vertex
certificate are in `target/readiness-repair/ghost-two-loop-main`.

Cross-section DOT exports now transport the runtime owner's existing cuts and
threshold candidates through the native FeynmanDiagram validator. The transport
retains every half-edge ID, edge orientation, rule slot, numerator fragment and
chosen loop coordinate, and checks all routing signatures before writing. It
neither applies Feynman rules nor enumerates replacement cuts. Runtime selections
continue to have one owner; the export view is temporary. Import uses the normal
model-aware artifact boundary and restores runtime metadata such as threshold
directives, names and mass overrides. Split diagnostic views retain their
existing explicit runtime-artifact contract.

The first diagnostic run with shared initial observations measures GL1 guard
peak RSS at 74.0 MB versus main's 84.0 MB. Thermal GL1 falls to 144.3 MB versus
main's 121.0 MB, compared with the earlier reconciled 194.9 MB. These are local
single-run diagnostics, not replacements for the pinned paired gate. Thermal
memory at the end of tensor parsing falls from roughly 174 MB to 93 MB.

Integrated counterterms had another migration mismatch: the old cleanup closed
momentum pairs before emitting scalar products, whereas the port disabled
rank-one contraction before a notation-only `to_dots` call. A captured d58 input
therefore multiplied independently integrated factors containing the same local
`uvind(0,2)` label. It produced `Q1(mu) Q5(mu)^2 Q7(mu)` and was correctly refused
by tensor admission. The integrated-output boundary now closes those internal
pairs before reinsertion, keeping open external ports and numerator factors.
The independent d58 and d70 Laurent references pass after this change. Cases
with remaining failures are retained separately in the local validation logs.

The next UV boundary comparison isolates an adjoint word whose last structure
constant becomes `f(a,b,b)` when its endpoint is substituted. It is exactly zero,
but retaining the word through materialization removed that factor's slots and
left a misleading rank mismatch against the other summands. The shared color
owner now removes this zero at either incident endpoint, with orientation and
nonzero controls. The captured before-materialization expressions and smallest
sum are retained in `target/readiness-repair/uv-materialization-inputs`.

Fifteen shared-library snapshots required only presentation updates. Every old
and new scalar/tensor expression parses to the same exact symbolic Atom (31
comparisons); DOT node/edge attributes, external labels and port inventories are
unchanged. Local rendering-port permutations are compared independently of
the node expressions. The proposal, exact proof and original source copies are
retained under `target/readiness-repair/shared-presentation-proposal`. Proposal
capture with forced snapshot continuation is not counted as validation; normal
runs keep exact snapshot assertions enabled.

After the endpoint-zero repair, d71 and d78 pass their unchanged RQFT Laurent
references. The d76 trace advances to scalar evaluation but reaches the strict
240-second limit while compilation is concurrent. Two scalar-evaluation stages
take about 97.94 seconds each and receive exactly the same numerator and
topology. That numerator retains repeated raw adjoint matrix pairs under
different internal labels. The color owner now applies their local Casimir
identity before materialization, with independent explicit-product and reversed
orientation tests. All 1,366 selected shared-library tests pass. The quiet d76
rerun passes the unchanged numerical reference in 33.295 seconds. Its 26 completed
scalar evaluations take 21.10 seconds in total, with the longest at 3.79 seconds;
the earlier timed-out run spent 203.75 seconds on only 20 completed evaluations.
The trace comparison is retained in
`target/readiness-repair/d76-evaluation-boundary/comparison.json`. This is a
diagnostic comparison with different concurrent load; the unchanged pinned
paired performance gates remain the final acceptance criterion.

Broader workspace discovery exposed two additional shared boundaries. Eager
opaque-word materialization contracted internal tensor ports but returned the
original, larger layout to the outer parser. Spenso now returns each realized
tensor together with its own layout; symbolic leaves retain their source layout.
A finite two-dimensional matrix test compares three parsing modes with distinct
endpoint vectors, exact component polynomials and no residual ports. The native
sampling/API reproducer passes after this repair.

Index canonicalization also erased duality before renaming and restored it by
function-name patterns. For symmetric metrics this could turn a closed color
cycle into an invalid non-self-dual square. Dual wrappers now stay in Symbolica's
canonical graph. The regression preserves a powered scalar spectator, checks
the exact color trace, label invariance and idempotence. All 1,380 selected
Spenso/Idenso checks and all 105 native generator unit tests pass, including the
Z-to-quark-pair reproducer. One canonical-string snapshot changes only by the
simultaneous swap of bound fundamental labels `d_3` and `d_4`; the exact
certificate is `target/readiness-repair/dual-canonicalization-reference-proof.json`.

The first full-workspace discovery run selected 5,562 tests. It reported 5,016
passes and 546 failures, with the existing 128 skips. Of those failures, 420 were
test-launch errors caused by a diagnostic rebuild replacing an executable while
the run was active; 126 were assertion/runtime failures. This is not a final
qualification run. Final validation must use frozen test artifacts and complete
its runtime checks before overlapping builds. The published `d3b91d31` integration
archive independently reproduces the Z grouping error and the extra ttH master;
its isolated control log is `target/readiness-repair/published-grouping-control.log`.

Finite-component sampling required a separate coefficient boundary in the native
generator. Charge-conjugate sewn self-energy words agree under five exact mass
probes but initially retain different nested scalar coefficient forms. The
sampling owner now collects exact complex-rational coefficients after tensor
contraction, protecting complete function calls and composite powers. It accepts
only strictly smaller intermediate forms. All 105 native generator unit tests
pass, including exact symbolic equality for that pair; the ttH selection test
then passes its unchanged two-master assertion. This repair is named
`feynkit-finite-sample-collection` (`ec4f5a19`). Normal paired generation and
evaluation timing gates still have to qualify the complete repair stack.

The remaining GL20 grouping discrepancy is explicitly retained. A read-only
export of the generated candidate state shows GL03 and GL06 with opposite
closed top-loop flows and weights -1 each. Pinned main GL38 has GL06's directed
topology and weight -2. Selecting only one candidate would change the physical
contribution. The comparison is recorded in
`target/readiness-repair/gl20-directed-topology-comparison.json`; the unique-master
assertion remains unchanged while the missing charge-conjugate grouping is
investigated.

The sunrise evaluator failure was localized to UV history arguments, not a
missing retained-scalar alias. Its complete source capture contains 21 declared
CT parameters and seven incomplete marker histories in outer roots; the child
histories are hidden in shared numerator definitions. The prepared repair makes
complete marker calls explicit numerator-family arguments so outer Approx and
prefix operations can transform them without distributing the numerator.
Validation of this repair is pending in the frozen workspace archive.

Documentation graph frames now reuse Linnest's SVG inspection decorator before
HTML link processing. This preserves drawing geometry and ordinary links while
turning native inspection links into accessible controls. Both the focused
inline-frame test and the two-build revision-provenance test pass (26.280 seconds
total); link validation remains enabled. Cargo Hakari verification and the
regenerated CI dependency graph include the new existing-library dependency.
