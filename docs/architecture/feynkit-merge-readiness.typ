= FeynKit merge readiness and performance

The qualification below belongs to checkpoint `6002674d`. A subsequent
user-requested migration retains `SymbolicTensor` through numerator processing;
its implementation and validation are recorded at the end of this document.
The earlier timing results do not qualify that new implementation.

This record follows the reconciliation of local FeynKit, remote FeynKit and
main, then the repairs required before merge. Runtime assertions, numerical
reference points, tolerances and the performance gate remain in force. A
successful producer build/upload is distinct from a successful runtime suite.

The implementation failures found by the frozen full-suite run are repaired:
all 24 affected acceptance tests and the related Vakint checks pass in a
70-test rerun. All 112 selected renormalization and UV integration tests also
pass, followed by all 73 evaluator and scalar-factor checks. All thirteen paired
performance gates pass with zero numerical mismatches. Ten audited
generation-bookkeeping snapshots still await user approval. Final CI preparation
and publication have not completed.

== Final implementation qualification, 17 October 2026

The normal Nix CLI at `6937daf944d7dcb2b24aaf7ee7a73a9aa5773bc3` passes every
declared performance gate against main `ff533fcd`. All production and test
sources match the reviewed implementation; later edits only update this record.
The original corpus, numerical tolerances, five-pair policy and timing limits
are unchanged. No build or test process from this task overlapped measurement.

#table(
  columns: (2fr, 1fr, 1fr, 1fr),
  [Metric], [Main median], [Repaired median], [Paired change],
  [GL1 generation], [3.9653 s], [4.5784 s], [+15.4%],
  [GL1 SingleParametric runtime], [157.75 µs], [135.28 µs], [−13.0%],
  [Thermal GL1 generation], [6.8730 s], [6.8564 s], [−0.55%],
  [Thermal GL1 SingleParametric runtime], [817.04 µs], [729.85 µs], [−10.0%],
  [Thermal generation peak RSS], [123.86 MB], [130.34 MB], [+4.9%],
)

The change column is the median of paired relative differences; it need not
equal the ratio of the separately reported medians. GL1 construction still has
an overhead of about 0.61 seconds. Its 15.4 percent paired increase now passes
the original 20 percent limit with the original noise allowance. Thermal
generation is essentially unchanged. The former runtime and memory regressions
pass their unchanged gates; this does not imply that every metric equals main.
All thirteen sealed reports, input hashes, raw trials and the consolidated
`performance-final-summary.json` are retained under `target/readiness-repair`.

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
  temporal-slash normalization. The concrete matrix already exists as
  `spenso_hep_lib::gamma0_weyl` and is registered in the component library.
  The missing behavior was the symbolic Idenso trace rule, which runs before
  components are constructed. Keeping both layers lets symbolic traces reduce
  early while component evaluation retains its independent matrix definition.
  Commit `fd1fef9f` contains the fix and
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
  certificate. At that earlier, uncollected-component checkpoint, disabling
  common-pair elimination or selecting Horner optimization did not reduce the
  counted arithmetic cost. The later thermal investigation below repeats the
  comparison after component collection and identifies a separate remaining
  conditional-sharing cost.

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

The GL20 comparison isolates opposite closed top-loop flows with weights -1
and -1, while main retains one master with weight -2. Both contributions must
survive. Finite sampled scalar coefficients now have an exact polynomial
identity certificate: conservative per-variable degree bounds determine a full
Cartesian grid of degree plus one points per variable, evaluated over exact
Gaussian rationals. The grid is capped at 256 points; unsupported functions,
floating-point coefficients, inverse/fractional powers and larger grids remain
unproved. Complete functions and powered numerator factors stay intact. A
polynomial that vanishes at the old probe points but not identically is an
explicit negative control. The restricted ttH NLO model now gives the unchanged
16 masters and 30 nonzero assignments, with 16 color-zero diagrams excluded.
The exact input model/options and captured samples are regression fixtures.
The comparison and former mismatch are retained in
`target/readiness-repair/gl20-directed-topology-comparison.json`.

The sewn incoming-state symmetry boundary also overcounted identical incoming
particles by an extra factorial. Raw cross-section external labels now remain
exact. Canonical comparison handles simultaneous permutations of the incoming
legs and their sewn partners, independently of left/right exchange. The
regression checks all eight symmetry configurations for a phi-four Born graph,
plus an asymmetric three-cycle and a one-sided negative control. Canonical-only
numerator grouping again permits zero numerical probes. All 110 native generator
unit tests pass after these changes; physical phase-space acceptance is checked
separately.

The sunrise evaluator failure was localized to UV history arguments, not a
missing retained-scalar alias. Its complete source capture contains 21 declared
CT parameters and seven incomplete marker histories in outer roots; the child
histories are hidden in shared numerator definitions. The prepared repair makes
complete marker calls explicit numerator-family arguments so outer Approx and
prefix operations can transform them without distributing the numerator.
The Taylor guard accepts constant complete marker arguments while still rejecting
nested retained families and dependence on the expansion variable. The focused
localizer test checks exact commuting Approx/Integrate/prefix transformations
and constant Taylor expansion. Both physical sunrise tests subsequently pass
with their original numerical assertions.

Documentation graph frames now reuse Linnest's SVG inspection decorator before
HTML link processing. This preserves drawing geometry and ordinary links while
turning native inspection links into accessible controls. Both the focused
inline-frame test and the two-build revision-provenance test pass (26.280 seconds
total); link validation remains enabled. Cargo Hakari verification and the
regenerated CI dependency graph include the new existing-library dependency.

A second full discovery run used frozen executables, eliminating the executable
replacement failures. It ran 5,566 tests: 5,523 passed, 41 failed and two reached
the unchanged timeout, with the existing 128 skips. This archive predates the
last grouping, initial-symmetry, model-serialization and UV-boundary repairs;
it is evidence for triage, not final acceptance. Its source identity and all 43
unique failures are retained with `workspace-repair2` in the task artifacts.

Model fingerprints must survive process-local symbol-registration order. Model
expression export now uses the existing canonical printer, preserving namespaces
and symbol semantics. A fresh-child-process regression registers Gamma/Gamma5 in
opposite orders and verifies identical complete Standard Model fingerprints and
JSON round trips. All 39 model tests pass. The existing canonical fixture
regenerator reproduces five graphs byte for byte except their model fingerprints;
three authored scalar collections also pass strict import after updating only
that derived field. The proof retains all graph numerators, cuts, thresholds,
slots and routing verbatim. No fingerprint check is bypassed.

The fixed-component bridge to Vakint now distinguishes a physical component
label from a summed Einstein index. A fourth power of the temporal component
therefore becomes four contractions with the fixed temporal basis, rather than
a dummy index occurring four times. Spatial contractions carry the Minkowski
metric sign. Independent formal-dimensional Wick expectations cover temporal,
spatial and mixed quadratic/quartic moments with an untouched powered spectator.
Both focused projection checks pass, including the existing on-shell temporal-dot
case. This GammaLoop bridge is separate from the shared Idenso gamma-zero trace
repair.

Imported edge mass overrides must precede construction of energy definitions.
The boundary capture shows the split-mass fixture retaining distinct `m1`/`m2`
on its edges while all registered energies still use the original model mass.
The import owner now installs edge metadata before building its parameter and
energy map, then compiles overridden-mass evaluators against the completed
parameter layout. The mapped compact-import regression checks the energy's
explicit mass dependence as well as its existing cut and routing assertions.
Both repeated-bubble mass-limit integrations subsequently pass, as do the
compact-import mass and energy regression checks.

The next focused physical run passes 24 of 30 checks, including both sunrise UV
cases, scalar Born phase-space comparisons, the one-through-four-loop unit
amplitudes, the massive-top inclusive NLO ratio, CLI state handling and the
updated LTD comparison. Six failures remain in that run: GL20 absorptive-sign
comparison, three label-dependent NLO checks, the pre-repair mass-override case,
and another incomplete authored diagnostic fixture. The NLO topology comparison
proves that old GL0 maps to native GL1 and old GL2 to native GL0 for both currents
and both quark masses, with unchanged weights. Selectors and edge-owned degree
bounds follow that map; all published graph cross sections and tolerances stay
unchanged. The mass-override and fixed-component unit checks both pass in the
subsequent focused run. Logs and source identities are retained with
`acceptance-repair3` and `import-mass-cache-repair`.

The exact binary64 six-photon point changes by 6 real-component ULPs and 34
imaginary-component ULPs after compact evaluator collection. Identical complete
input cards on main and the repaired CLI give exactly the same Arb-evaluated
output: `(-5.127347252278009e-10, 2.8273655564307803e-10)` at the public binary64
output boundary. The binary64 snapshot is refreshed for operation ordering and
the test now also enforces that independent main Arb result with exact equality.
Neither assertion uses a relaxed tolerance. Both cards, JSON results and the ULP
comparison remain under `target/readiness-repair/photon-point-*`.

The photon spin-sum zero was a sampling-chart failure before evaluator execution.
All explicit helicities also failed; event metadata reports an invalid interior
point for `E3+E4-Q0-Q1`, with value 58.1663 at zero radius. Native default routing
stored external attachments as `(source, target)` but selected the target while
its documented legacy convention requires the outgoing source. Correcting that
endpoint restores a physical cut momentum as the loop coordinate. A sewn-box
unit check verifies the spatial origin in the centre-of-mass frame and preserves
topology/numerator identity. All 110 generator unit tests pass, and the unchanged
photon helicity-average integration test passes after the repair. The corrected
native GL0 selection also passes both original electron-positron smoke references.
GL20 still reverses its absorptive component and remains an open acceptance
failure; its numerical reference is not changed.

The remaining d54 timeout is reproduced with no concurrent task builds under
the unchanged 240-second limit, and again on current code with boundary logging.
The trace records 146.47 seconds between canonicalization and completed tensor
projection stages, including individual intervals of 23.34, 18.68, 9.83, 65.62
and 14.62 seconds. These timings locate work to investigate, not a qualified
speedup. The factored backend inputs and full trace remain under
`target/readiness-repair/d54-tensor-inputs` and `d54-current-trace`.

The analytic timeout was localized further on a captured d54 input. Native tensor
projection took 66.35 seconds, but the shared Idenso phases accounted for less
than five seconds. CPU sampling identified repeated Symbolica replacements.
The existing dot converter replaced one occurrence, then scanned the entire
factorized numerator again. The isolated conversion took 54.500 seconds. A
single bottom-up traversal now performs it in 0.0467 seconds, retaining the
existing reservation of dummy indices and independent scopes for powered sums.
All 242 selected graph, generator and tensor-reduction checks pass, including
reciprocal powers, existing indices and a new factored-product round trip. The
complete d54 test passes its unchanged exact Laurent reference in 97.101 seconds,
under the original 240-second cap. The named repair is
`feynkit-linear-dot-conversion` (`102b9a12`). The measurements are diagnostic;
normal Nix CLI performance qualification remains a separate gate.

Numerical numerator comparison previously performed full symbolic Lorentz
contraction before substituting sample momenta. Its preparation pass now limits
structural work to color; finite-component execution still performs the complete
Lorentz calculation. The four-loop photon-current case falls from a timeout
above 240 seconds to 28.073 seconds; pinned main took 44.671 seconds. Its ten
extra retained masters were exactly vanishing weighted classes. Native grouping
now removes complete scalar-weight cancellations, including rational model
normalizations over a common denominator, while leaving graph numerators
factored. The original main regression for the FullSM GL51--56 rational class
is reproduced in the native grouping owner, with a nonzero control.

The remaining symbolic rescaling checks also distinguish scalar conjugation
from tensor conjugation. A conjugated CKM/model parameter is a valid scalar
ratio; conjugation of a tensor is not. The registered broadcasting conjugation
head now validates its argument instead of being rejected unconditionally.
Cross-sample comparisons normalize only their proven scalar ratios over a
common denominator. No source numerator, numerical-sample budget or tolerance
is relaxed. All 118 native generator unit tests pass after these repairs.

The final grouping boundary repair is `c1f0ad7f`. Isolating a color coefficient
can expose a closed fundamental metric loop; that isolated coefficient now
passes through the existing tensor contractor before finite samples are
compared. Canonical tensors collect only numerical content of retained sums,
so opposite and quarter-turned bodies share a scalar normalization. Completed
group weights use Symbolica's exact zero certificate with both sampling budgets
set to zero. This also removes the six charged-ghost/G0 classes whose weights
were `1-sqrt(1/2)*sqrt(2)`. The original vacuum counts are restored: 300 in both
canonical modes, 247 in symbolic sampled modes, 288 in fully numerical mode,
and 67 in the nonfactorizable subset. The four-loop count remains 161.
The fresh snapshot capture records these counts but is not a passing snapshot
validation; representative-dependent weights are still being reviewed.

The on-shell agent's validated nested-vector admission fix `e2d55fbf` is included
as `589480e5`. It normalizes intrinsic nested-vector contractions before tensor
interface admission, preserving opaque payloads and callback boundaries. Its
original validation record is `on-shell-review-validation-record` (`0893e362`);
that record certifies the agent's own OS stack, not this task's final revision.

== Conditional evaluator performance follow-up

The normal Nix producer at `c0e78d263f03f6cb16c6055230e12eaf36da347d`
predates the last grouping changes. Its paired GL1 runtime and memory metrics
pass; generation remains inconclusive near the gate (confirmation medians
3.982 seconds on main and 4.766 seconds on the candidate). Thermal GL1 still
fails SingleParametric runtime: 820.25 versus 988.24 microseconds, a paired
20.60 percent increase. All numerical corpus comparisons pass. These immutable
reports are in `target/readiness-repair/performance-provisional-*-report` and
must not be relabelled as final-revision results.

The two complete thermal pipelines produce the same 32 orientation selectors
and thermal function calls. Their scalar compiler inputs first differ in the
factorized numerator bodies. The candidate has fewer static instructions, but
common-pair elimination moves repeated arithmetic outside mutually exclusive
residue branches. A replay of the exact programs with their production
selectors executes 86,001 additions and 82,155 multiplications on the candidate,
against 43,619 and 64,865 on main. The replay uses synthetic nonzero scalar
inputs to count work; it is not a physics oracle.

Disabling common-pair elimination in that replay reduces the candidate totals
to 38,684 additions and 54,551 multiplications, with the same thermal calls.
A normal-Nix diagnostic using the existing setting measures 732 and 735
microseconds at the two thermal corpus points, versus 812 and 818 on main;
nonthermal GL1 measures 143 and 139 microseconds, versus 155 at both points.
Every pilot value passes the unchanged corpus tolerance. The proposed production
policy, named `feynkit-conditional-evaluator-cost` (`20a09717`), skips this
optimization by default only for guarded runtime residue
dispatch and preserves explicit user settings. Direct DAG sharing remains
enabled, including during merges of scalar output programs. All 73 selected
evaluator and scalar-factor tests pass. The normal Nix binary at
`f6d6b810cd15acde53003f0228bf5b7301dac492` passes twelve of the thirteen paired
workloads with no numerical mismatches. GL1 guard evaluation improves from
159.9614 to 136.8547 microseconds (14.46 percent), and thermal evaluation from
838.7351 to 744.5759 microseconds (10.91 percent). The remaining GL1 generation
gate is inconclusive: 3.9841 versus 4.7989 seconds, a paired increase of 19.85
percent. Its component-factor collection takes 426 milliseconds across eighteen
distinct bodies. Captured-argument discovery takes only 0.4 milliseconds and
post-algebra contraction 12.5 milliseconds; those stages do not explain the
remaining cost. The separate topology check preserves all 66 inventories over
five pairs and improves from 0.9089 to 0.4836 seconds (46.8 percent). Timing
thresholds and paired-run policy remain unchanged.

The user approved the two semantic generation-reference corrections: uniform
formal-Nc CP ratios, and ten canonical masters when symmetric polarizations
are explicitly enabled. A controlled toggle retains twelve masters when that
option is disabled. Five further reference edits change source IDs only; their
diagram counts, complete scalar-weight ledgers and evaluated totals agree
exactly before and after.

The user also approved retaining the original 49 scalar graph frames and the
GL20 event frame. Scalar point, edge-numerator, UV-profile and precision
assertions keep their original inputs and tolerances. A separate native
generation check compares all 49 colored graphs, local factors, overall weights,
and initial, physical and threshold cut incidence. GL20 still checks the unique
native master, its three-loop topology, all 255 loop bases, four cuts and four
active channels before importing the original directed frame for the complete
event oracle. No numerical event reference is changed. All five focused
fixture checks pass. The frozen full-source run at `52293751` completes 5,587
tests: 5,559 pass and 28 fail, with 128 existing exclusions and no timeout.

The remaining tensor failures exposed two narrower Vakint boundaries. A
projected external coordinate dot such as `dot(delta(cind(1),mink(D)),p(3))`
is a closed scalar, although its coordinate vector alone is not. The scalar
parser now uses Spenso's existing fixed-basis matcher only within that closed
dot. It still rejects open external ports, loop-dependent coordinate dots and
out-of-range components. Separately, generated dot dummy index 2 could collide
with a retained numeric source index 2 when both entered GammaLoop's restored
UV index chart. Vakint's existing allocator now reserves raw and represented
source slots before introducing dummy labels. Commit `7d1ec9b5` preserves
these fixes; all 22 targeted tests pass, including an independent formal-D
rank-one Wick projection with a retained powered spectator. The complete
acceptance rerun then exposed two adjacent boundaries: scalar integration
appended indices beside native vector representations, and a retained source
slot on a loop momentum could become a raw numeric backend index. Moving the
existing native-slot repair into the common dot converter and preserving those
source slots resolves both. Commit `9f59b7d9` passes all 70 selected tests,
including every one of the 24 previously failing acceptance cases and the full
Vakint tensor-reduction mode tests. The identical complete paired-index and
fixed-component GL00 inputs also generate successfully with main's binary.

Commit `5025e085` avoids repeated tensor scans within a certified scalar
collection scope and skips summed-base normalization when its lookup catalog
is empty. Guard, power and inverse ownership are unchanged. Replaying all
eighteen original finite components takes 0.1945 seconds instead of 0.2406
seconds, with exactly identical expressions and idempotent outputs. This is an
isolated diagnostic, not a replacement for the paired generation gate. The
profiling split attributes about 102 milliseconds each to factor protection
and numeric normalization, versus 5 milliseconds to Horner collection and
18 milliseconds to factor extraction. A temporary short-alias variant was
discarded; it did not improve the cost enough to justify changing expression
ordering.

The follow-up `2c3d6e3d` also skips occurrence-locality scans when the complete
scope contains no guards or numerator families and the current factor has no
tensor-bearing inverse. Under exactly those conditions the existing predicate
is identically false. The same eighteen replay expressions remain identical
and idempotent; isolated time falls to 0.1687 seconds. All 73 evaluator and
scalar-factor tests pass after this change. Whole-workspace Clippy, including
all targets with warnings denied, also passes. The 112-test broader UV run
passes with its original 240-second per-test limit; its durations overlap
compilation and are not performance-gate measurements.

The archived documentation diagnostic failure is a relocation artifact:
trybuild compares relative source paths against the relocated absolute path.
The same binary passes at the original workspace path without changing any
diagnostic reference. A missed scalar canonical fixture receives only the
deterministic model fingerprint, verified by strict import and an exact diff.
The photonic binary64 reference follows the changed compiler operation order;
fresh arbitrary-precision evaluation still equals main's exact stored value,
and both equality assertions remain enforced.

== Pending generation-reference review

Ten prepared snapshot updates remain unapplied. They retain every diagram
count and exact assertion. The audit is in
`target/readiness-repair/feyngen-boundary-proposal/remaining-reference-proposal.json`,
with the original and proposed strings and reproducible evidence for each case.

All sixteen closed-ghost colored masters retain exactly the same individual
weights. Two provenance entries of weight minus two become four entries of
weight minus one; the total remains minus 29/2. The asymmetric twelve-master
case also has the same complete scalar-weight ledger after exact Casimir
normalization. These changes do not alter a physical contribution.

The vacuum snapshots select different representatives of equivalent massless
current classes. A down-quark photon numerator is one quarter of the former
up-quark representative, accompanied by four times its weight. For the neutral
current, define `Pu = 9*cw^4 - 6*cw^2*sw^2 + 17*sw^4` and
`Pd = 9*cw^4 + 6*cw^2*sw^2 + 5*sw^4`. The new numerator is `Pd/Pu` times the
old one and its weight is `Pu/Pd` times the old weight. Fourteen exact symbolic
identities, the actual model couplings and the selected source diagrams certify
the normalization changes; numerical sampling is disabled in those proofs.

In the four-loop case all 161 colored topologies match, and 159 master weights
are identical. The other two representatives reverse one odd-length vector
fermion cycle. The independent SU(3) color contraction changes from 7/3 to
minus 2/3, while charge conjugation reverses the Dirac sign. The numerator ratio
is therefore 2/7, compensating the weight change from 18/7 to 9. Exact affine
momentum maps have determinant one and retain the sewn external pair. The raw
weight-sum change of 90/7 is fully accounted for by those two representatives;
their weighted physical contributions are unchanged.

== Retained symbolic tensors: implementation and API notes

The symbolic numerator should retain its tensor interface until finite-component
lowering. Inspecting an immutable expression is different from replacing its
owner with an `Atom`: the latter loses logical port order, the rank of zero,
and reusable validation facts. Scalar coefficients, source-model recipes whose
ports are not instantiated yet, and the final scalar evaluator have distinct
representation requirements.

The migration changes runtime edge/vertex fragments,
`SymbolicExpression<State>.expr`, FourD sectors and projection witnesses, and
retained numerator function bodies to `SymbolicTensor<PartialStructure>`. The existing
`FnMapEntry` becomes generic in its body type; scalar definitions still hold
`Atom`. Color/gamma processing, direct-3D numerator mapping, Taylor coefficients,
final normalization and evaluator component lowering retain typed bodies.
Opaque family calls in scalar dispatch roots must not be inferred as complete
scalar tensor expressions: their body catalog owns the hidden interfaces.

Against main `ff533fcd`, the concrete storage changes are:

#table(
  columns: (1.1fr, 1fr, 1.4fr),
  inset: 6pt,
  [Stage], [Main], [Retained-tensor branch],
  [Runtime edge/vertex numerator], [`Autogen<Atom>`],
    [`Autogen<SymbolicTensor<PartialStructure>>`],
  [Numerator processing state], [`SymbolicExpression<State>.expr: Atom`],
    [Same state machine, typed tensor payload],
  [Shared numerator definition], [`Arc<FnMapEntry>` with `rhs: Atom`],
    [`Arc<NumeratorDefinition>` with typed RHS; scalar `FnMapEntry<Atom>` remains separate],
  [Local-4D sector and source witness], [`FourDSector.atom`, `FourDTerm.numerator`,
    `CanonicalUvTerm.source_numerator`: `Atom`],
    [Each retains the tensor and its logical interface],
  [Component preparation], [Raw expression parsed into `ParsingNet`],
    [Consuming typed `into_network` returns `Canonicalized<ParsingNet>`],
  [Compiled evaluator], [`AliasedAtom`, scalar `FunctionMap`, Symbolica program],
    [Same scalar execution engine after explicit component lowering],
)

This migration targets repeated interface discovery and representation conversion
before compilation. It cannot by itself remove component contraction, scalar
factor collection, evaluator compilation or numerical execution. Retained
metadata also occupies memory; releasing owned planning caches at component
lowering limits their lifetime, but does not establish a measured memory saving.

Stage timings also need matching boundaries. Raw-root cooking and admission now
precede `preprocess_tensor`, so its cumulative `normalization_done` timer omits
that work. A separate `evaluator_stack_raw_root_admission` event records it.
Retained family bodies actually avoid their former re-admission. The enclosing
evaluator preprocessing interval and paired full-process generation time remain
the comparable totals; a smaller normalization timer alone is insufficient
evidence of a speedup. The exact timer map and outstanding traversals are in
`target/symbolic-tensor-flow/performance-boundaries.json`.

Concrete Rust API difficulties found during implementation:

- Symbolica's Rust `AtomCore` trait is sealed, so a downstream tensor carrier
  cannot implement it directly. Its checked `replace` takes
  `TensorRule` rather than the Atom replacement builder. Arbitrary scalar or
  kinematic substitutions currently need `expression()` followed by
  `with_checked_expression`. Transformations which change an interface need
  `with_transformed_expression`; `with_rewritten_expression` is reserved for
  established tensor identities and must not bypass arbitrary-substitution
  validation.
- There is no tensor Taylor-series API. Direct-3D currently keeps the existing
  scalar series engine internal and restores each coefficient through the
  retained source interface. A zero coefficient keeps its tensor rank. This
  boundary needs an ergonomic checked series API before eliminating the
  remaining application-side adaptation.
- Existing Atom helpers can hide index-transport assumptions. The Taylor
  helper's Dirac prepass used to decode every cooked symbol on return. Native
  input then acquired raw indices alongside retained encoded ones. The shared
  helper now preserves existing index and dimension symbols and decodes only
  symbols introduced by its own prepass; fully native input needs no decode
  traversal. Raw/native/mixed-encoding controls check the same exact Dirac trace
  identity and spectator factors, with explicit Taylor dependence inside the word.
- The old `expand_dots` returned `Atom`, even when only a closed scalar dot
  inside an open tensor changed. Its typed return must retain the untouched
  open-port identities. A checker that only rediscovers a rank cannot certify
  the order of anonymous ports; callback-induced rank changes remain errors.
- The old `to_network` specialization only supported one ordered-index and
  mixed-numeric tensor type. The new consuming `into_network` entry accepts the actual
  component tensor library, function library and parsing settings. It checks
  the resulting external ports against the retained interface and returns the
  root's `CanonicalLayout` with the network. Component coordinates, dimensions
  and generated shadows must use that same logical axis order; checking only
  the port multiset is insufficient. Shorthand
  expansion, fresh-index reservation and library materialization still belong
  to the component parser; this change alone does not eliminate that parse.
- GammaLoop's owner-bearing index syntax is more permissive than the native
  `AbstractIndex` syntax. Typed bodies retain reversible encoded indices;
  decoding occurs at component or explicit expression-export boundaries.
  A generic uncook operation cannot silently change the carrier's index type.
- Persistence previously had no carrier support. Archives must store payload
  and logical interface, including ranked zero and unresolved axis order.
  Decode remaps symbols and validates those parts once. Proof and planner
  caches are never serialized. Existing equality also observes the metric and
  composite flags, so those flags are retained to preserve exact round trips.
- Declared source interfaces need `from_declared_parts`; merely attaching
  `checked_parts` does not certify an arbitrary external payload against its
  declared ports. Runtime DOT exports add `num_interface` only when the
  expression cannot recover its logical interface, notably for ranked zero.
  Source scalar overrides remain scalars. Full state version 19 rejects older
  positional archives; standalone scalar evaluator layouts are unchanged.
- `with_lorentz_dimension` specifically changes four-dimensional Lorentz
  spaces. General dimension substitutions are currently less convenient and
  require explicit checked interface transformation.
- A no-op reduction can establish new completion facts. Mathematical equality
  alone must not discard its returned carrier when reusing an `Arc`;
  `same_reduction_state` distinguishes that case without traversing cached graphs.
- A mixed-space `from_signature` example currently constructs canonical
  argument order with different logical metadata, then fails its own admission
  check. The standalone reproducer is
  `target/symbolic-tensor-flow/idenso/from-signature-mixed-layout.json`.
  Production admission uses `infer` or the checked declared-interface owner;
  the constructor issue is recorded separately rather than bypassing validation.

The component handoff releases its owned planning caches before parsing the
component graph. Another retained definition may still share those caches;
only the full-process RSS measurement can establish the lifetime's actual cost.
Additional Vakint normalization is admitted and required to be scalar before
typed multiplication, so arbitrary configured factors cannot inherit an
unrelated tensor-interface certificate.

Evaluator family preparation reads retained logical ports even when the payload
is zero. A ranked zero keeps an empty sparse component table and a typed shadow;
only a closed scalar zero may collapse directly to scalar zero. Standalone open
zero and mismatched-port controls are rejected, while a matching contraction
and an ordinary scalar zero remain accepted.

The initial 530-case runtime pass exposed six failures at these new boundaries.
Three local-4D cases identified the Taylor prepass encoding error described
above. A temporal-component fixture used the raw shorthand `Q(edge,0)` and now
declares the same component as `Q(edge,cind(0))`. An owned-dot fixture used
literal mass 1: blanket mass tagging changed `cind(1)` but left the numerator's
minus-one coefficient unchanged. The fixture now retains `mass_scalar_1` as the
physical mass symbol and asserts that every edge still has numerical mass 1;
its original numerical comparison remains intact.

The sixth failure concerned the exact zero of four independent Taylor terms.
Typed composition retained equivalent scalar contractions as both `g(Q,Q)`
and `dot(Q,Q)`, which the coefficient comparison treated as different keys.
The existing `to_dots` conversion puts these keys in the same notation before
the original exact coefficient assertion. The formal Lorentz dimension remains
symbolic and graph numerators remain factorized. An independent replay removes
the literal common numerator from all four terms and proves that their scalar
denominator weights sum to zero with zero numerical or finite-field samples;
the frozen earlier binary passes the same original physical test. Evidence is
in `target/symbolic-tensor-flow/collective-scalar-proof.json`.

The on-shell agent's named scalar-power repair `c8ff5b8f` is included separately
as `589e447b` (`feynkit-shared-scalar-powers-20261017`). It lets generated
intrinsic scalar powers retain their existing certificate after an identity
kernel, without re-admitting the complete tensor. Both original regressions pass.

Canonical FeynKit diagrams and raw model recipes remain source-artifact
representations. Runtime graph construction admits their instantiated fragments
once; subsequent subgraph assembly composes the retained carriers through the
shared graph-selection fold.

Validation and performance measurements for this migration are pending. Existing
physics references, tolerances and performance thresholds remain unchanged.
Artifacts and logs are retained under `target/symbolic-tensor-flow`.
