= FeynKit merge readiness and performance

Production revision `fd9ed336` passes the normal Nix CLI build and all numerical
checks in the thirteen original CLI cases. Twelve cases pass their performance
gates; GL16 SingleParametric runtime remains inconclusive under the unchanged
noise policy. The separate topology comparison preserves all 66 entries in each
of five pairs. The production-equivalent test child `663ee049` passes all 1,054
selected Idenso tests and all 664 selected consumers, strict workspace Clippy,
formatting and the documentation registry check. Later documentation checks also
pass at `575606be`. A later quiet ghost follow-up of the `663ee049` consumer
binary passes all eight exact checks without the earlier large slowdown,
subject to the build-feature limits recorded below.

The complete 5,661-test workspace sweep, earlier quiet ghost comparisons and
successful 105-output local CI preparation belong to the earlier `e6a3bb46`
checkpoint.
Their source identities remain separate from the newer qualification. At this
documentation checkpoint, CI is enabled, but preparation of all 105 outputs for
the final exact source, cache upload and push/publication remain pending. The early `9c76b732` preparation was
intentionally stopped as superseded; it is neither active nor a completed pass.
The remaining inconclusive gate is not waived, and this is not a merge-ready
checkpoint.

This record follows the reconciliation of local FeynKit, remote FeynKit and
main, then the repairs required before merge. Exact assertions, numerical
tolerances, sample budgets and performance limits remain enforced. Individually
justified reference corrections are recorded below, including the historical
six-photon binary64 update checked independently with arbitrary precision.
A successful build or focused rerun is distinct from a complete runtime or
performance qualification.

== Latest shared fixes, 18 October 2026

The on-shell agent's two additional named shared fixes are imported separately:
`aa118f7852cef1ffc67743cc18f91d023c00764a` maps to
`b3b3d533876c07fe817dbcd6b99e0037c990b556` for native fixed-selector recognition,
and `3db4014f5fcc23256e427ddf8bf7e112b55bd675` maps to
`fd9ed3361517fd36bcfc6059eb626122f92924f0` for uniformly reversed contracted
`gamma0` sandwiches. No on-shell prescription ancestry is imported. The child
`663ee04900a1ffe358b94630f107a5f3a801335d` only disambiguates a test map's type.

Frozen `663ee049` passes all 1,054 selected Idenso tests: the exact earlier
1,048-test inventory plus six new controls, with 24 unchanged skips. Strict
workspace Clippy, formatting and the documentation registry check also pass.
All 664 selected GammaLoop and tensor consumers also pass in 635.552 seconds,
with the exact previous inventory retained, no failures or timeouts, and unchanged
source hashes. The normal Nix CLI build passes at `fd9ed336`; its source, exit
status and binary identity are retained in
`target/symbolic-tensor-flow/shared-selector-temporal-cli`. The first CI attempts
exposed the test-map type ambiguity and a missing compile-time fixture. The map
correction above preserves the exact assertion; `ea221df3` adds the 819-byte
captured temporal chain to Idenso's existing compile-time test assets. Receipts are under
`target/symbolic-tensor-flow/shared-selector-temporal-validation-2`, with import
mappings in `onshell-api-coordination.json`. The reviewed fixture-source digest
and documentation registry pass at `575606be`. At this documentation checkpoint,
final exact-source CI preparation remains pending; the newer CLI performance
results follow separately from the historical e6 full-workspace sweep.

== Current fd9 CLI qualification

The aggregate at
`target/symbolic-tensor-flow/performance-final-fd9ed336-aggregate-1/report.json`
verifies all thirteen fresh, predeclared case attempts against main
`ff533fcd`, with zero evidence errors and every numerical check passing. It
records twelve performance passes and one inconclusive case. No e6 trials are
substituted or pooled into these results.

#table(
  columns: (2fr, 1fr, 1fr, 1.2fr, 0.7fr),
  [Metric], [Main median], [fd9 median], [Median paired change], [Gate],
  [GL1 generation], [4.08535 s], [4.53491 s], [+10.7000%], [Pass],
  [GL1 SingleParametric runtime], [158.9323 µs], [137.6361 µs], [−14.3240%], [Pass],
  [GL1 generation peak RSS], [83.6362 MB], [99.2707 MB], [+16.9070%], [Pass],
  [Thermal generation], [7.13200 s], [6.45814 s], [−9.5860%], [Pass],
  [Thermal SingleParametric runtime], [825.4038 µs], [744.0008 µs], [−10.1271%], [Pass],
  [Thermal generation peak RSS], [112.7711 MB], [155.1483 MB], [+25.8318%], [Pass],
)

MB denotes decimal megabytes. Median paired changes are not ratios of the two
reported medians. Thermal RSS retains a measurable increase: the median paired
absolute difference is 31.04 MB. It passes the original combined relative and
absolute rule, including the 64 MiB floor; a passing gate does not mean unchanged
memory use. All GL1 and thermal gates pass without changing their limits.

GL16 alias SingleParametric runtime has a +8.3801 percent median paired change
and 2.6763 percent robust paired noise. The median is below the 15 percent
allowance, but the median plus three times noise crosses it, with the absolute
floor also exceeded. This produces the unchanged inconclusive status, not a
pass or a confirmed regression. The screen is not a statistical confidence
interval. Its other metrics pass.

The separate topology comparison preserves all 66 inventory entries in each of
five pairs. Fresh-process medians are 0.913906216 seconds on main and
0.962177605 seconds at fd9, with a +5.1210 percent median paired change. Timing
includes model import, generation and state save, and excludes integrand
construction. This is a supplemental inventory comparison, not a fourteenth
performance gate or a full numerator-equivalence certificate. Its receipt is
`target/symbolic-tensor-flow/performance-final-fd9ed336-topology-1/topology-report/report.json`.

== Current CPE policy controls

The SingleParametric compiler keeps the scoped `cpe_iterations = 0` policy for
conditional residue dispatch when the user has not specified a value. CPE means
common-pair elimination. Explicit settings remain honored, and Horner ordering,
ordinary common-instruction elimination and stack reuse still run.

An fd9 two-pair ABBA comparison of explicit zero versus one CPE iteration gives
mixed results. GL16 SingleParametric runtime changes by −5.87 percent in one
pair and +0.29 percent in the other; thermal changes by −0.48 and +0.44 percent.
All 48 physical checks pass and stored settings are verified. The original
runner exits with a post-check GL16 history-path error; the preserved history
was subsequently verified at its actual cross-section path. The failure receipt
is retained. This comparison neither establishes a safe repeatable gain nor
isolates the scoped production override: explicit zero also changes other modes
from their unspecified fixed-point policy. One iteration is a bound on each
construction or merge, not the native fixed-point default. Evidence and limits
are in `target/symbolic-tensor-flow/fd9-cpe-counterfactual/analysis.json`.

A separate fd9 ABBA control compares the production default with explicit
fixed-point CPE. All 48 physical checks and eight stored-settings checks pass.
Separate invariants also hold: the preoptimization source shapes,
non-SingleParametric operation vectors and other scalar-build vectors match
across policies. These invariants describe source/operation metadata, not
byte-identical serialized programs. On GL16,
fixed-point CPE makes SingleParametric runtime about 5.2 percent faster in both
pairs. On thermal GL1, it makes that runtime 32.7--38.5 percent slower and
evaluator construction 46.7--56.4 percent slower. These are bounded two-pair observations on the original
physical points, not a replacement gate or an exact decomposition of the
main-versus-branch difference. The complete receipt is
`target/symbolic-tensor-flow/fd9-cpe-fixedpoint-counterfactual/conclusions.json`.

Both programs have fewer static additions and multiplications with fixed-point
CPE, yet thermal evaluation becomes slower. Static operation counts therefore
do not establish executed cost. The measured tradeoff supports retaining the
existing scoped zero policy. It does not support removing that policy globally
or adding workload-specific heuristics from these two cases. No new production
policy or gate threshold is introduced by these diagnostics.

== Earlier frozen e6 qualification

`e6a3bb4664d6927c46b0e463eee0850fc6285fa4` passes 1,048 of 1,048 selected Idenso
tests, with all 24 existing skips preserved. Strict workspace Clippy completed
in 23.57 seconds and the documentation check passed. The exact d53 and d76
assertions pass in the repaired diagnostic capture. These results are recorded
in `internal-chain-tests-2.log`, `loop-square-clippy-1.log`,
`loop-square-docs-1.log` and `ghost-performance-audit/loop-square-stage-analysis.json`
under `target/symbolic-tensor-flow`.

The complete e6 workspace run passes all 5,661 tests with 128 unchanged test
skips, four slow cases and zero timeouts, in 3,176.327 seconds. Its source tree
remained unchanged through execution; the post-run diff is empty. Run
`237e15f9-d4a7-46e8-a39a-38c920b739a2`, the JUnit hash and complete result inventory
are recorded in `target/symbolic-tensor-flow/color-proof-full3-summary.json`.
The four slow cases are the jet-generation test and three bubble-UV wrappers,
also present in the earlier slow category. d75 is not an additional slow case:
it takes 108.600 seconds versus 115.907 seconds in the earlier `880ad11d` run.
This run overlapped normal Nix CLI and local CI artifact preparation. Its
durations establish correctness and timeout evidence, not isolated performance
qualification.

The normal e6 Nix CLI build also passes. Its binary SHA-256 is
`9423734bfb3f42b948053654295d4d8691fe7073038d4bad339ad094d0c70a31`, and it reports
`gammaLoop 0.1.0`. The exact source revision, build exit status, binary hash and
version receipts are retained in
`target/symbolic-tensor-flow/native-adjoint-loops-cli`.

Local CI preparation at e6 successfully realizes all 105 selected outputs,
covering static checks and compiled runtime producers. This preparation performs
no external cache upload and no Nix runtime-launcher execution; the complete
Cargo runtime suite above is separate evidence. Receipts are retained under
`target/symbolic-tensor-flow/e6-ci-cache-preparation`.

At e6, all thirteen original CLI cases completed: eleven passed and two were
inconclusive, with all numerical checks passing. The coverage index verifies all
evidence with zero errors. It preserves four original sealed reports and nine
predeclared continuations without pooling partial trials or selecting retries.
The initial interrupted invocation remains recorded with exit 143; no cause is
assigned here. The aggregate and exact per-case reports are indexed by
`target/symbolic-tensor-flow/performance-final-e6a3bb46-aggregate-1/report.json`.

At that earlier e6 checkpoint, GL1 passed every unchanged gate while retaining
measurable generation and memory costs:

#table(
  columns: (2fr, 1fr, 1fr, 1.2fr, 0.7fr),
  [Metric], [Main median], [e6 median], [Median paired change], [Gate],
  [GL1 generation], [3.95567 s], [4.46103 s], [+12.1132%], [Pass],
  [GL1 SingleParametric runtime], [157.9431 µs], [135.0273 µs], [−14.1263%], [Pass],
  [GL1 generation peak RSS], [82.14528 MB], [100.78208 MB], [+20.5537%], [Pass],
)

MB denotes decimal megabytes. Paired changes are medians of paired relative
changes, not ratios of the two reported medians. GL1's median RSS increases by
about 18.64 MB; it passes the original combined relative/absolute rule, including
the 64 MiB absolute floor. Passing does not mean unchanged memory use.

The e6 GL16 alias case was inconclusive for SingleParametric runtime: its
median paired change is +1.0652 percent with 6.3495 percent robust paired noise.
The e6 thermal GL1 case was inconclusive for SingleParametric runtime and
peak RSS. Its runtime median improves by 8.3747 percent, with 6.7367 percent
paired noise; RSS increases by 26.4654 percent (31.38 MB absolute), with
18.5406 percent paired noise. Thermal generation passes, improving by
11.8838 percent.

These metrics trigger the unchanged high-noise rule: three times relative
paired noise reaches the applicable relative allowance, and the median paired
absolute change plus three times absolute noise exceeds its absolute floor.
The runtime allowance is 15 percent with a 1 µs floor; RSS uses 20 percent and
64 MiB. Both conditions matter, including for thermal runtime's negative median.
This is a robust-noise screen, not a statistical confidence interval. The two
e6 cases were neither passes nor confirmed regressions. The later fd9 results
are recorded above without replacing this evidence or relaxing the policy.

Two predeclared fixed-state GL16 ABBA diagnostics retain the original four
methods, two points and numerical checks: one unbound, then one with both
binaries pinned to the same singleton CPU 192. All 64 point-result JSONs pass
and both saved states remain unchanged. Pinning does not remove SingleParametric
variance: its batch coefficient of variation remains about 2.65--5.96 percent,
while the other methods remain about 0.25--1.56 percent. This does not isolate an SMT, cache or
precision-path cause and does not replace either inconclusive gate.

The unbound slow episode stays on one CPU, with negligible observed leader
run-queue waiting and no cgroup throttle increments. Hardware counters report
full running coverage, but point-window instruction counts interpolate 100 ms
intervals; they are not exact per-batch or symbolic-program execution counts.
`RUSAGE_CHILDREN` under the profiler does not account for the payload here, and
its zero software migration/context-switch counts contradict `/proc`. Those
fields are not used to establish absence of work or scheduling effects. The
observations and limits are recorded in `observed-analysis.json`,
`singleton192-comparison.json` and `diagnostic-counter-quality.json` under
`target/symbolic-tensor-flow/gl16-final-variance`.

A separate original thermal-pipeline ABBA diagnostic passes all numerical
checks. One candidate run exposes a missed generation peak: external observed
RSS reaches 174,297,088 bytes while the native 100 ms sampler reports
136,749,056 bytes. The peak lies after normalization and before the tensor
boundary is finished, before finite component execution. This locates the stage,
not its internal allocation cause. Some same-process `VmHWM` observations also
decrease, so the final value is not treated as an exact lifetime bound. The
native sampling rule and every gate remain unchanged; evidence is in
`target/symbolic-tensor-flow/thermal-final-variance/observed-analysis.json`.

The scoped `cpe_iterations = 0` policy remains in place. Main lets common-pair
elimination repeat to a fixed point; on the branch's changed factorized sources,
sharing across mutually exclusive residue branches can shrink the static program
while making each selected residue execute more work. The override preserves
branch-local arithmetic, Horner ordering, ordinary DAG sharing and stack reuse,
and honors explicit user settings. This predates retained tensor storage: a
normal-Nix diagnostic at `c0e78d26` measured thermal runtime at 969/974 microseconds
with default CPE and 732/735 with zero, passing the original corpus comparisons.
Later scoped production passed its paired runtime gate. These are historical
`c0e78d26` measurements, not a default-versus-zero comparison at e6. The fd9
zero/one and fixed-point controls in the current CPE section above provide the
later measured tradeoff supporting the retained policy.

The separate five-pair topology comparison preserves all 66 inventory entries
in every pair. Fresh-process medians are 0.910630843 seconds on main and
0.950127492 seconds at e6, with a +5.23124 percent median paired change. Timing
includes model import, generation and state save; it excludes integrand
construction, DOT export and inventory verification. This supplemental check is
not one of the thirteen metric gates or a full numerator-equivalence certificate.
Its receipt is
`target/symbolic-tensor-flow/performance-final-e6a3bb46-topology-1/topology-report/report.json`.
These topology results qualify e6 only; current fd9 qualification and the
remaining final-source CI/upload/publication steps are recorded above.

The earlier full run of `880ad11dcce02724c178eac628e0f28dad6e17fa` finished in
3,175.357 seconds with 5,639 passes, five failures, five slow cases and zero
timeouts; 128 existing skips remain. Its compile-time source tree stayed
unchanged throughout execution. Run `ab964f08-8546-476b-9dc4-6ba4066eeb39` and its
JUnit identity are recorded in
`target/symbolic-tensor-flow/native-complete3-full-summary.json`.

`8f24069f3eb8dee3b71883d294e87a9e529a8168` passes all 126 focused source/interface
tests, strict workspace Clippy, formatting and the documentation check. This
covers the two Q3 validation-order failures and the subsequently identified
selector boundary. Q3 source validation now precedes generic tensor admission
without cancelling invalid source functions across separate payloads; the
original arity, edge-bound, dimension, component and cancellation guards remain.

The selector repair admits exact integrated/external vector selectors once
through the shared source codec. Their Rust builders are fallible; graph
promotion uses `map_representation_dimensions`. Raw compound dimensions and
native encoded dimensions can then identify the same vector without re-inferring
the retained numerator. Frozen failure/isolated repair evidence and canonical
controls preserve the unselected same-head vector, declared logical axes,
ranked zero and unchanged source.

The user approved all ten remaining bookkeeping updates, applied at
`79c954094b97adba6927cfc4054f5f434603ae68`. The three complete generation test
functions pass in 53.346 seconds, covering all ten exact assertions, rather than
stopping at the first snapshot in each function. This commit changes only those
reference strings; its production code is identical to `8f24069f`. The approval
scope and the limits of the supporting proofs are recorded below. These focused
results address all five failures of the earlier full run; they are not reported
as a new complete post-fix workspace run.

At the earlier `8f24069f` checkpoint, the normal Nix CLI build and its sealed
performance corpus against main completed. All thirteen cases preserve their
numerical references; twelve pass their
performance gates. The remaining GL16 alias case is inconclusive: its
single-parametric runtime has a 7.4 percent median paired increase, below the
15 percent limit, but its paired-noise screen crosses that limit. This requires
a further isolated comparison, not a relaxed threshold. Final CI preparation
and publication have not completed.

The same run measures GL1 generation at 3.996 to 4.569 seconds (14.1 percent
paired increase) and single-parametric evaluation at 159.49 to 135.81 microseconds
(15.5 percent paired decrease). Thermal GL1 generation improves from 7.040 to
6.371 seconds and its single-parametric evaluation from 834.27 to 743.19
microseconds. Peak memory increases from 81.37 to 99.12 MB for GL1 and from
117.92 to 145.91 MB for thermal GL1; both remain within the unchanged gates.
These costs remain part of the performance review even when a gate passes.

The separate topology check preserves all 66 diagram inventory entries in every
pair. Its fresh-process median time changes from 0.907 to 0.953 seconds, with a
5.6 percent median paired increase. This includes model import and saving the
generated state, excludes integrand construction, and is not a full numerator
equivalence certificate. All thirteen gate reports, binary identities and the
five topology pairs are retained in
`target/symbolic-tensor-flow/performance-source-boundary-8f24069f-1`.

== Integrated ghost regressions and repairs

Two isolated quiet pairs, in alternating `16787f71`/`8f24069f` order, pass the
unchanged exact d53 reference within the original 240-second deadline. Median
wall time rises from 19.515 to 166.810 seconds, about 8.5 times; median child CPU
time rises from 19.373 to 165.514 seconds. This is a confirmed computation cost,
not an explanation based solely on the earlier concurrent full-suite timings.
The relevant test macro, settings, 78-diagram source, Standard Model input and
exact reference are identical. Preserved binaries, hashes, environments and
trial records are under `target/symbolic-tensor-flow/ghost-performance-audit`.

The three bubble UV wrappers remain near their earlier 150-second whole-test
times. A separate two-pair d76 comparison also confirms extra computation:
`16787f71` takes 31.688/31.758 seconds and `8f24069f` takes 70.197/62.732 seconds.
All four exact assertions pass. Both ghost cases require repaired-build timing;
whole-suite durations and the old GL1 measurements do not replace those checks.

The d53 diagnostic traces have the same stage call counts. Coefficient sizes
grow before the unchanged Laurent-series engine, which accounts for about
75.6 percent of the sampled candidate CPU profile. For the matched
`current=16C4`, `given=sXg`, `reduced=DeO` context, the first observed growth is
at color simplification before Lorentz projection or scalar vacuum evaluation.
The earlier path closes the color traces into Casimirs and indices; the
unrepaired `8f24069f` path retains traces and structure constants. Three matched projected numerators
subsequently grow from 1,422/2,864/2,017 to 23,298/25,744/49,858 raw string bytes.
Small native controls now prove the signed explicit-product/adjoint-chain
equivalence with identical external ports. Each individual trace should give
the same nonzero Casimir coefficient, while their symmetric sum must be a
rank-two zero. The explicit forms reduce correctly; the chains remain unresolved
through both retained and freshly checked second passes. Opening only the small
chain restores the exact expected results. The full captured input reproduces
the same unreduced color inventory; missing process-local data in legacy cooked
labels limits exact text comparison of that larger replay.

The implemented repair belongs to Idenso's product-color owner. It exposes one
connected adjoint chain through the existing reserved dummy allocator and
commits the result only when an existing cross-factor identity reduces it.
Unrelated chains and spectator factors remain intact. The production change
passes all 1,040 existing selected Idenso tests. After correcting construction
of chain placeholders in the new test fixtures, all six selected adjoint-chain
tests pass, including the three new controls. These runs include the shared
on-shell fixes `428558eb`, `2c3d3260` and `4138873d`, imported as separate commits
`036dcfac`, `07b3cecc` and `5d06ad55`. These initial controls preceded the combined
checkpoints and quiet comparisons below. The comparison preserves numerator
factorization and the original assertions; disabling checks
or expanding the graph numerator is not a remedy.
The exact boundary pairs are retained under
`ghost-performance-audit/runs/current-d53-inputs-1/boundary-comparison`.

The combined `58ab7be2` checkpoint passes all 1,045 selected Idenso tests,
strict workspace Clippy and the documentation check. Its d53 and d76 numerical
assertions also pass. In that bounded stage capture, d53's largest cleaned
coefficient returns to the earlier 309,148 bytes; d76 still reaches 7,517,268
bytes. Those test durations overlap build work and are not timing qualification.

d76's remaining growth first appears at color simplification in
`current=1680`, `given=sLI`, `reduced=Dmi`. A two-matrix adjoint chain retains an
internal endpoint/side contraction where the baseline trace reduces immediately.
An exact native local replay proves
`chain(A,B; F(B),F(C)) = C_A g(A,C)`, with only A and C public. The baseline
`Tr(F(C) F(A))` has the opposite sign; the captured overall imaginary factors
compensate. The follow-up uses the existing two-f owner inside chain
simplification. A shared edge guard distinguishes summed indices from fixed
components; the latter are not evidence for a color-index sum. The matched
inputs, local proof and independent review are retained under
`ghost-performance-audit`; the native checks and completed quiet comparisons
are distinguished below.

The first full Idenso run of that follow-up passed 1,047 of 1,048 tests; its new
scalar-closure control exposed a normalized `f^2` form. The correction routes
both the separate and squared structure constants through the same guarded
contraction. The exact `8*C_A` assertion is retained and extended with weighted
closures and fixed-index controls. The combined e6 correction now passes all
1,048 selected tests and both exact ghost assertions.

In the repaired e6 stage capture, the largest cleaned coefficients are 309,148
bytes for d53 and 255,806 bytes for d76. The previously observed 18,059,032-byte
d53 coefficient at `8f24069f` and 7,517,268-byte d76 coefficient at `58ab7be2`
are absent. Stage call counts match the paired `16787f71` baseline:
47/36/11/36/36 for d53 and 33/26/7/26/26 for d76 across stages one through five.
All four exact assertions in that baseline/candidate capture pass. These are
expression-size and correctness observations from stage logging, not quiet
performance samples. The separate quiet measurements follow below. The
stage-only log does not establish per-call symbolic input equivalence or a
full-source performance certificate.

== Quiet ghost and scalar checks, 18 October 2026

Two alternating pairs per ghost case compare frozen `e6a3bb46` with the earlier
retained-tensor baseline `16787f71`, preserving the original 240-second deadline.
All eight exact assertions pass, with no task-owned build or test overlap. The
host is shared; these are quiet diagnostic comparisons, not a claim of exclusive
CPU access or a substitute for the thirteen-case sealed CLI corpus.

#table(
  columns: (1fr, 1fr, 1fr, 1.3fr),
  [Case], [16787f71 median wall], [e6a3bb46 median wall], [Median paired change],
  [d53], [19.80694 s], [20.55913 s], [+3.7994%],
  [d76], [32.59799 s], [18.91189 s], [−41.9735%],
)

The change column is the median of the paired relative differences, not the
ratio of the separately reported medians. The earlier 8.5-fold d53 increase is
no longer present; its remaining roughly 0.75-second (3.8 percent) overhead is
still reported. d76 improves by about 42 percent. These measurements do not
attribute the residual d53 cost to an unmeasured stage. Raw trials and the exact
source/fixture scope are recorded in
`target/symbolic-tensor-flow/ghost-performance-audit/loop-square-quiet-final-analysis.json`.

A later quiet follow-up uses the `663ee049` consumer binary, whose production
code is `fd9ed336`, against frozen `16787f71`. All eight exact checks pass under
the original 240-second deadline, with no timeouts. The nine relevant source,
input and reference groups are reverified byte-identical. Candidate d53 takes
20.971--21.890 seconds versus 20.087--20.131 on the baseline; candidate d76 takes
19.585--20.206 seconds versus 32.928--33.404. The earlier large slowdowns do not
return in these binaries. This is a bounded follow-up, not a sealed CLI gate or
a feature-identical causal timing comparison: 663 omits Idenso's optional
`reference-cases` hooks present in the saved 167/e6 builds, and 167 has broader
historical Python/portable-payload linkage. The hooks are dormant in these tests,
but their TLS/code-layout presence can still affect cost. Small differences are
therefore not attributed to the late shared fixes. The e6 comparison above
remains separate; the follow-up receipt is
`target/symbolic-tensor-flow/ghost-performance-audit/consumer663-quiet-analysis.json`.

A separate scalar GL01 check compares `880ad11d` with e6 using sequential ABBA
runs, byte-identical test/fixture inputs and the unchanged 240-second deadline.
All four original assertions pass. The median times are 11.0530 and 11.0754
seconds respectively, so the approximately 40 percent increase observed during
concurrent full runs is not reproduced. This check is distinct from the sealed
CLI gates. Its raw results and scope are recorded in
`target/symbolic-tensor-flow/scalar-gl01-duration-audit/quiet-analysis.json`.

== Measured GL1 preparation cost and proof reuse

The existing instrumentation at `8f24069f` attributes about 0.915 seconds of
additional GL1 work to tensor preparation, offset by about 0.270 seconds less in
scalar compilation. Within preparation, the nineteen normalization calls add
about 0.402 seconds. Work after normalization adds about 0.392 seconds, including
0.320 seconds of finite-component factor collection; those figures overlap and
must not be added together. Raw-root admission adds about 0.139 seconds, root
dot conversion 0.038 seconds, and family combination saves about 0.061 seconds.
Sampling supports substantial inference/validation cost but does not isolate
each validation substage. These diagnostic intervals explain the work performed;
the sealed paired timings above remain the performance comparison.

An unchanged `map_representation_dimensions` previously rebuilt checked parts,
discarding the carrier's completed reductions and cached observations. Its
follow-up keeps the mapper and callback checks, then uses the existing
`finish_parts` owner when both the payload and complete logical structure are
unchanged. This preserves established facts while still validating unchecked
inputs and finishing pending product normalization. A zero whose represented
dimension changes takes the normal changed-value path. This is a separate
optimization; its regression tests remain passing in the e6 1,048-test run,
while its contribution to GL1 timing still needs qualification. Evidence and
the validation plan are in
`target/symbolic-tensor-flow/dimension-noop-proposal`.

== Historical qualification at 6002674d, 17 October 2026

At that checkpoint, all 24 affected acceptance tests and related Vakint checks
passed in a 70-test rerun. All 112 selected renormalization and UV integration
tests passed, followed by all 73 evaluator and scalar-factor checks. The normal
Nix CLI at `6937daf944d7dcb2b24aaf7ee7a73a9aa5773bc3` passed every declared
performance gate against main `ff533fcd`, with zero numerical mismatches. Its
production and test sources matched the reviewed `6002674d` implementation.
The original corpus, numerical tolerances, five-pair policy and timing limits
are unchanged. No build or test process from this task overlapped measurement.

#table(
  columns: (2fr, 1fr, 1fr, 1fr),
  [Metric], [Main median], [6002674d median], [Paired change],
  [GL1 generation], [3.9653 s], [4.5784 s], [+15.4%],
  [GL1 SingleParametric runtime], [157.75 µs], [135.28 µs], [−13.0%],
  [Thermal GL1 generation], [6.8730 s], [6.8564 s], [−0.55%],
  [Thermal GL1 SingleParametric runtime], [817.04 µs], [729.85 µs], [−10.0%],
  [Thermal generation peak RSS], [123.86 MB], [130.34 MB], [+4.9%],
)

The change column is the median of paired relative differences; it need not
equal the ratio of the separately reported medians. At that checkpoint, GL1
construction had an overhead of about 0.61 seconds. Its 15.4 percent paired
increase passed the original 20 percent limit with the original noise allowance.
Thermal generation was essentially unchanged. The former runtime and memory
regressions passed their unchanged gates; this does not imply that every metric
equaled main or qualify the later migration.
All thirteen sealed reports, input hashes, raw trials and the consolidated
`performance-final-summary.json` are retained under `target/readiness-repair`.

== Historical reconciliation baseline

The main baseline is `ff533fcd91d80577e8372c398c17a5924613e73c`.
The reconciled comparison is `d3b91d31461e50d91dce4702d68c45294faca350`.
Both CLI executables were built with the normal Nix `gammaloop` package and
uploaded to the configured cache. The candidate supplies identical cards and
inputs to both binaries. Gammabench is pinned by
`benchmarks/performance/gammabench.rev`; every gate uses the declared corpus,
five paired discovery trials and, where required, five confirmation trials.
Generation includes a fresh CLI process; runtime uses warmed fixed points.
Timing comparisons wait for this task's own builds and tests to finish.

The frozen main uses the registry release of Symbolica 3.0.1; the later
`8f24069f` candidate uses community revision `c78cffc`. Sealed CLI timings
compare those complete binaries. Replaying both captured sources through one
community compiler instead isolates source structure and optimization settings;
it does not by itself separate their contribution from upstream compiler
changes. A measured runtime improvement remains valid even when that causal
split is unresolved. The separate ghost comparison from `16787f71` to
`8f24069f` uses the same community compiler revision on both sides.

All thirteen numerical corpora agreed in the initial reconciliation comparison.
Eleven cases passed their performance gates. Its confirmed regressions were:

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

The source comparison uses three immutable checkpoints: main `ff533fcd`, the
previously qualified FeynKit implementation `6002674d`, and the first retained-
tensor candidate `48b8acfa`. This distinguishes the earlier FeynKit integration
from the newer tensor-ownership migration. The full path/symbol evidence is in
`target/symbolic-tensor-flow/main-candidate-rust-flow.json`.

+ #strong[Generation and model.] Main's `ProcessDefinition` directly owns its
  PDG, filtering and grouping options; `diagram_generator` constructs a
  `SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>`. Its model uses
  `ArcParticle`, `ArcPropagator` and `ArcVertexRule`. The reconciled branch uses
  `GenerationProcess` and `GenerationOptions`; `feynkit_generator::Process`
  returns `GenerationResult` containing `Vec<FeynmanDiagram>` and a
  `GenerationReport`. The canonical model supplies `ParticleId`, `PropagatorId`
  and `VertexRuleId`. These are earlier architectural changes; tensor retention
  does not change the generator or enumerate diagrams again.
+ #strong[Runtime graph, routing and cuts.] Both routes use GammaLoop's
  `Graph`, containing `HedgeGraph<Edge, Vertex, HedgeData>`, a
  `LoopMomentumBasis` and cut/surface data. The loop basis still carries
  `loop_edges`, `ext_edges` and `EdgeVec<LoopExtSignature>`. The branch imports
  `FinalizedCut` and `FinalizedTopologyThresholdCandidate` records from the
  canonical diagram and delegates routing to `feynkit_graph::routing::MomentumRouting`.
  The external momentum rows and source frames must survive this conversion.
  The tensor migration changes fragment admission at this boundary; the LMB
  and cuts implementations remain unchanged from `6002674d`.
+ #strong[Numerators.] Main stores `Autogen<Atom>` on edges and vertices and an
  `Atom` in `SymbolicExpression<State>` inside `Numerator<State>`. The earlier
  FeynKit branch also stores Atoms, constructing temporary Idenso tensors for
  algebra and exporting their expressions afterward. The migration retains
  `SymbolicTensor<PartialStructure>` in those same owners. The Atom remains its
  payload; logical ports, zero rank and reduction facts now travel with it.
  `UltravioletGraph::numerator` becomes fallible because typed composition must
  validate the fragments it combines.
+ #strong[CFF and UV.] `AmplitudeGraph::preprocess` still orders CFF generation,
  UV construction, sampling channels and optional thresholds. `GraphThreeDSource`
  supplies `GeneratedThreeDExpression<Esurface, Hsurface>`;
  `OrientationExpression` carries exact `LinearEnergyExpr` maps and `CFFVariant`
  data. `OperationNode`, `ComputeStore` and `ComputeNode` schedule the hedge-poset
  forests. Main already has `Arc`-shared numerator definitions and separate
  scalar roots in `Integrands`; this sharing is not a new tensor-migration
  optimization. The migration makes `FnMapEntry<E = Atom>` generic and uses
  `NumeratorDefinition = FnMapEntry<SymbolicTensor<PartialStructure>>` for
  tensor bodies. Direct-3D mapping, `FourDSector`, `FourDTerm` and
  `CanonicalUvTerm` retain those bodies. The forest combinatorics, cut keys and
  residue maps stay with their existing owners. Vakint and scalar Laurent
  series remain explicit backend boundaries.
+ #strong[Finite components.] Main preprocesses an `AtomCore` input, builds a
  `ParsingNet`, and executes finite components. The earlier branch inserts
  temporary Idenso admission/algebra before that parse. The migration passes
  the retained tensor to `into_network`, which returns
  `Canonicalized<ParsingNet>`. Here `ParsingNet` is a Spenso `Network` with
  `NetworkStore<ParamTensor<ShadowedStructure<Aind>>, Atom>`, `ExplicitKey<Aind>`,
  `Symbol` and GammaLoop's `Aind`. Its `CanonicalLayout` maps component storage
  back to the declared logical axes. This lowering occurs during generation,
  before sample evaluation.
+ #strong[Scalar evaluators.] `EvaluatorStack` still contains single-parametric,
  iterative, summed-function-map and summed `GenericEvaluator` variants.
  `GenericEvaluator` consumes scalar `FnMapEntry<Atom>` definitions and
  Symbolica `ExpressionEvaluator` programs; `ParamBuilder<T: FloatLike>` supplies
  numeric values. The precision lanes and eager, SymJIT and external backends
  remain the same. Component factor collection and the guarded common-pair
  elimination policy are earlier performance repairs. Equivalent expression
  organization can still change compiler work and executed instruction counts.
+ #strong[Sampling and events.] `ProcessIntegrand` evaluates
  `MomentumSample<T>` / `BareMomentumSample<T>` with their momenta, Jacobian,
  orientation and parameterization branch. `SamplingMapContext` and proposal/
  host plans retain their existing responsibilities. The sampling subtree is
  unchanged from `6002674d`; the migration adds no per-sample tensor inference
  and does not change the Monte Carlo estimator or event accumulation.
+ #strong[Persistence.] Main's full-state format is version 15; the current
  full-state format is version 20, including retained native diagrams and global
  prefactors. Regenerate older full states. The payload stores expressions and
  logical interfaces, remaps symbols and validates them on load, and excludes
  planner/proof caches. Standalone amplitude version 12 and cross-section version
  15 remain unchanged because they store already
  lowered scalar evaluators. Canonical diagram/model recipes remain source
  artifacts; runtime tensor interfaces are retained at their admission boundary.

The model/generator split and stable IDs are architecture choices. Correct
routing, source provenance, tensor interfaces and scalar admission are necessary
invariants. Repeated expression transport, repeated validation of an already
certified value, and overlapping cached graphs are implementation costs to
measure and reduce. Tensor retention changes ownership; it does not require a
new CFF scheduler, sampler or numerical evaluator.

For the measured GL1 evaluator comparison, both routes use the same frozen graph,
runtime inputs, selected evaluator settings,
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
  interpreter uses the exact production selector IDs with synthetic nonzero
  scalar parameters and counts more interpreter operations for each candidate
  orientation. This supports a branch-locality hypothesis; it does not measure
  instruction execution on the physical corpus or establish the whole runtime
  cause. The separate CLI corpus comparisons and timing gates provide the
  numerical and runtime evidence. At that earlier, uncollected-component checkpoint, disabling
  common-pair elimination or selecting Horner optimization did not reduce the
  counted arithmetic cost. The later thermal investigation below repeats the
  comparison after component collection and identifies a separate remaining
  conditional-sharing cost.

The original GL1 discovery showed substantially larger runtime and generation
costs before the temporal-trace and lifetime fixes. Those reports remain
retained, including censored confirmation trials, but censored measurements are
not ranked as qualified performance results.

== Repairs preceding the retained-tensor migration

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
299 inversions and 76 calls. With synthetic nonzero scalar parameters and the
nine production selector assignments, the diagnostic interpreter takes
3,564--4,042 steps, versus 4,682--5,249 before the change and 3,187--3,847 on
main. These are synthetic program-replay counts, not instruction measurements
on the physical corpus or elapsed-time gates. The partial diagnostic interpreter
tolerates nonfinite synthetic intermediates, so these counts are heuristic work
estimates. Scalar powers, lazy guards and numerator products remain
factorized; no compiler optimization setting or benchmark threshold changes.
Counting each n-ary addition/multiplication by its number of binary operations
and every other instruction visited by that synthetic replay once gives ranges
of 3,556--4,344 on main, 5,072--5,814 before component collection, and
3,789--4,392 afterward. This second diagnostic records synthetic operand cost
and control flow; it is not a physical-corpus execution count or cycle model.
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
`feynkit-finite-sample-collection` (`ec4f5a19`). At that stage, the complete repair
stack still needed its subsequent normal paired generation and evaluation gates.

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

== Historical conditional evaluator performance follow-up

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
residue branches. A diagnostic replay of the exact programs uses all 32
production selector assignments with synthetic nonzero scalar parameters. It
visits 86,001 additions and 82,155 multiplications on the candidate, against
43,619 and 64,865 on main. These synthetic program counts are neither a physics
oracle nor measured instruction execution on the thermal corpus. The separate
normal-Nix corpus comparisons and paired runtime gates below remain the
correctness and performance evidence.

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

== Approved generation bookkeeping references

The ten updates approved by the user are applied at `79c95409`: two ghost
provenance ledgers, one asymmetric scalar-weight ledger, six vacuum
representative/normalization snapshots, and one four-loop representative-weight
snapshot. All diagram counts, exact snapshot assertions, numerical references,
tolerances and comparison sample budgets remain unchanged. The original/proposed
strings and supporting evidence are in
`target/readiness-repair/feyngen-boundary-proposal/remaining-reference-proposal.json`;
that historical proposal remains labelled as prepared, while the applied-site
manifest and successful three-function rerun are under
`target/symbolic-tensor-flow`.

For the ghost updates, all sixteen colored masters retain their individual
weights. Two provenance entries of weight minus two become four entries of
weight minus one; the total remains minus 29/2. This certifies the colored-master
inventory and scalar-weight bookkeeping. It is not an independent full-numerator
or kinematic equivalence proof for every ghost contribution. The separate
physical acceptance references remain enforced. The asymmetric twelve-master
case has the same complete scalar-weight ledger after exact Casimir normalization.

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
  [Integrated counterterm], [`Series<AtomField>`],
    [Typed pole and signed finite projections; scalar Laurent series is temporary],
  [Global graph numerator/projector], [`Atom`],
    [Typed global factors; direct native-to-runtime transport],
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
After family shadows are installed, finite dots inside scalar function
arguments are lowered before the component parser, which treats those
arguments as opaque. `evaluator_stack_raw_root_scalar_components` records this
additional root boundary separately.
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
- A tensor is not an ordered Atom key. Exact inventory comparisons use an
  explicit exported-expression key together with its ordered interface
  metadata; inventing `Ord` from the expression alone would ignore ranked zero
  and logical axes.
- There is no tensor Taylor-series API. Direct-3D currently keeps the existing
  scalar series engine internal and restores each coefficient through the
  retained source interface. A zero coefficient keeps its tensor rank. This
  boundary needs an ergonomic checked series API before eliminating the
  remaining application-side adaptation.
- Tensor composition is fallible: `multiply`, `try_add`, `try_sub` and `pow`
  must check compatible interfaces and index multiplicities. Replacing these
  calls with ordinary Atom operators would discard those checks. Component
  tensor containers have Atom-like forwarding APIs, but those APIs do not
  establish the symbolic carrier's contraction and logical-layout invariants.
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
  Internal symbolic materialization and the integrated residual-contraction
  observer now use `AbstractIndex` directly, removing their legacy decode,
  parse and recook detour while retaining all structural checks.
- Persistence previously had no carrier support. Archives must store payload
  and logical interface, including ranked zero and unresolved axis order.
  Decode remaps symbols and validates those parts once. Proof and planner
  caches are never serialized. Existing equality also observes the metric and
  composite flags, so those flags are retained to preserve exact round trips.
- Declared source interfaces need `from_declared_parts`; merely attaching
  `checked_parts` does not certify an arbitrary external payload against its
  declared ports. Runtime DOT exports add `num_interface` only when the
  expression cannot recover its logical interface, notably for ranked zero.
  Source scalar overrides remain scalars. Full state version 20 rejects older
  positional archives; standalone scalar evaluator layouts are unchanged.
- `with_lorentz_dimension` specifically changes four-dimensional Lorentz
  spaces. `map_representation_dimensions` now owns general dimension changes
  in both payloads and declared ports, including ranked zero. A compound raw
  dimension is not directly representable by `Dimension::try_from`; source
  selectors and graph promotion need their proper admission/mapping boundary.
- Custom vector heads must be registered before native admission. Head-name
  selectors require the qualified Symbolica name; a bare string can silently
  select no vector. The tensor README now registers its vectors in the parser's
  caller namespace and passes the registered symbol to `with_integrated_head`.
  Its isolated example check distinguishes the intended `K.K*P.P/D` result
  from an unreduced `(K.P)^2`, using the frozen e6 libraries. Evidence is in
  `target/symbolic-tensor-flow/tensor-readme-example`.
- A no-op reduction can establish new completion facts. Mathematical equality
  alone must not discard its returned carrier when reusing an `Arc`;
  `same_reduction_state` distinguishes that case without traversing cached graphs.
- Symbolic dimensions do not provide a total canonical storage order. A scalar
  substitution can change the syntax-derived order of explicit ports from
  `[D:1, 4:1]` to `[4:1, D:1]` without changing either port. Checked replacement
  now compares complete explicit port identities with their multiplicities and
  keeps the original logical layout. Anonymous ports retain their positional
  checks. The separate repair is `7f673da6`
  (`feynkit-explicit-port-interface-order`); all 1,028 selected Idenso tests pass,
  including asymmetric finite-component and invalid-interface controls.
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

The native audit found instantiated values still stored as Atoms before the
runtime boundary in the reconciled FeynKit implementation. These native
FeynKit owners are absent from main `ff533fcd`; their before/after comparison
is against the pre-migration FeynKit repair `6002674d`.
Live `DiagramVertex::numerator`, `DiagramEdge::numerator`,
the aggregate `FeynmanDiagram` numerator, its prefactor and projector now retain
the carrier. `NumeratorCache` keeps completed vertex/propagator tensors;
`PreparedNumerator` keeps exact, sampling-source and color-valued sampled
tensors. Finite Lorentz components and certified scalar ratios remain Atoms.
Model templates and actual DOT/JSON strings are the source boundaries. The
existing `NumeratorInterface` codec now belongs to `feynkit-graph` and is
shared with GammaLoop; it records unrecoverable logical ports explicitly.

Standalone `Amplitude::terms` and `SquaredAmplitude::expression` also retain
tensors. `TensorReductionTerm::tensor` includes foreign-representation
spectators and remains typed, while its generated orbit weight is scalar.
The reduction result retains a zero with the source interface when no terms
survive. GammaLoop transfers native diagram carriers directly instead of
exporting them through raw parse payloads. Its `GlobalPrefactor` now keeps both
the numerator and projector typed as well. The completed `880ad11d` run above
qualifies the covered native paths and records its remaining failures; the
earlier runtime-only checkpoints did not qualify these additions.

Three further API gaps required operations on their existing native owners.
`dirac_adjoint` must dualize logical slots and retain a zero's interface;
same-representation anonymous axes need the existing materialization/axis
ownership machinery before adjunction. General representation-dimension
mapping must update both indexed and compact metadata, including zero;
ordinary scalar substitution cannot transform an encoded compound color
dimension. Checked simultaneous symbol renaming must transport the explicit
interface, reject collisions and preserve callback validation. These are not
safe uses of an Atom replacement followed by fresh inference. Python transfers
reuse `TensorExpression::from_structured` to keep the same carrier.

Current production `fd9ed336` has twelve passing CLI cases and one inconclusive
GL16 SingleParametric case, with all numerical checks and evidence verified.
Its normal Nix CLI build and five-pair topology comparison are complete.
The `663ee049` test child passes 1,054 native and 664 consumer tests; its separate
quiet follow-up passes eight exact ghost checks with the feature limits above.
The complete 5,661-test sweep, earlier quiet ghost comparisons and completed
105-output preparation remain explicitly tied to e6. The CPE controls above support keeping
the scoped policy and do not replace the remaining inconclusive gate. At this
documentation checkpoint, final exact-source preparation of all 105 CI outputs,
cache upload and push/publication remain pending; the early 9c76 preparation was
stopped as superseded. Artifacts and logs are retained under `target/symbolic-tensor-flow`.

The first immutable candidate `48b8acfa` passed a fresh 533-case core sweep and
workspace Clippy. The broader 151-case physics/state/backend run then reported
129 passes, 21 failures and one timeout. This prevents treating the earlier core
pass as complete qualification. The failures exposed explicit-port order matching,
finite dots retained inside scalar function arguments, a legacy color collector
which still exported native indices, and a diagnostic parser's dummy-name
collision. The scalar-argument and diagnostic-collision repairs passed focused
production cases; the native port-order repair passed the full Idenso package.
The later results below supersede those application failures.

The follow-up frozen boundary candidate `65a1b835` passes all ten selected
application/regression controls, including projected GL00, the quark local-IR
profile, both UV wrappers, scalar-argument lowering and the two ghost cases.
The d76 test passes with numeric-only stage logging in 32.418 seconds and with
both log sinks disabled in 32.718 seconds on that same frozen binary. The
earlier 240-second failure belongs to a different candidate and remains an
unresolved intermittent observation; these passing controls do not justify
changing the deadline or declaring a measured performance improvement.

The end-to-end audit found two additional retained-representation gaps after
that checkpoint. `IntegratedCts` now retains typed pole and signed finite
projections instead of a persistent `Series<AtomField>`; its existing consumers
need only these projections. The same native series engine performs temporary
Laurent arithmetic and truncation. Completed projections retain the source
interface even when they vanish, and nested forests/localizers consume them
without re-admission.

Local-4D denominator decomposition now restores every completed numerator
through the original carrier, preserving its logical layout. The parser still
performs temporary expression arithmetic and re-inference. Selecting every
`den` with the native tensor collector can distribute positive denominator
polynomials which belong to the factorized numerator; it cannot replace the
existing rational-shell parser without an API extension. The regression keeps
the exact positive polynomial factor, open-port order and ranked-zero checks.

Frozen `16787f71` builds all 102 workspace test binaries. Its combined core,
physics, state and backend run passes 685 of 687 tests in 805.811 seconds,
with no timeouts; all existing tests in that selection pass. The two failures
are new regressions: one compares indexed contractions with equivalent dot
notation and now uses the existing exact normalization; the other exposed an
anonymous-axis witness bug. The checked rewrite fallback attempted tensor-only
inference on scalar product factors. It now uses the same structural filter as
ordinary product inference, while retaining complete result/scope/multiplicity
checks. The same frozen d76 test passes in 32.215 seconds. The two corrections
and the subsequent native-generator migration were later exercised by the
completed `880ad11d` run recorded above.

Earlier ghost d76 timing diagnostics were inconsistent. A rebuilt binary exceeded
240 seconds with tracing disabled yet passed in about 35 seconds with narrow
tracing. Neither observation establishes an unconditional algorithmic timing.
Scalar-only diagnostic events now isolate stages without expression printers;
the original reference and 240-second timeout remain unchanged. The incomplete
48b8 CLI build was stopped as superseded, with completed shared Nix outputs
retained. No performance result from the previous qualified implementation is
being relabeled as a measurement of this migration.

Before the `880ad11d` full run, the native migration and backend-interface
corrections passed workspace-wide type checking, strict Clippy with warnings
denied, formatting and the documentation check. The isolated Idenso run with
`reference-cases` passed all 1,029 selected cases in 33.492 seconds; the later
workspace run also included its four feature-gated bincode tests. The formerly
failing scalar rescaling guard now verifies exact nonsymmetric components and
rejects a changed anonymous matrix head; no ambiguous-axis control is removed.

Source review also identified unmeasured costs for the next profiling pass:
source export re-admits an expression to decide whether interface metadata may
be omitted; identity index transport can rebuild planning state; numerical
samples repeat some sample-independent color preparation; and diagram export
can compute the same metadata twice. These are candidates for measurement,
not an attribution of GL1's historical slowdown. The completed frozen workspace
run and normal CLI build do not replace fresh paired timing. The selector and
Q3 repairs have since passed their canonical focused checks recorded above.

Frozen native candidate `53312d8b` built all 102 workspace test binaries and
passed 1,519 of 1,535 selected native/interface tests. The sixteen failures
identified boundaries requiring correction: a malformed toy singlet
fixture adds a colored summand; graph artifact validation needs an exact scalar
coefficient certificate; several tests still expect raw Python expressions;
and the new backend controls expose encoded dimensions and invalid TOML output.
The subsequent repairs kept complete tensor admission, exact numerical oracles,
and interface assertions. The broader run of these same binaries was diagnostic;
source-dependent documentation tests read the live compile-time workspace and
needed a rerun there after its review digests settled. Archive remapping alone
does not isolate those tests.

The corrected `28f5259f` focused run passed 1,706 of 1,717 tests in 365.232
seconds. All five Standard Model acceptance cases, the GL20 event snapshot,
global tensor round trips, native spin-sum controls and the scalar-vacuum
reference passed. At that checkpoint, remaining repairs covered a toy propagator's
obsolete `Slash(P)` spelling, typed Python comparisons, compact block-probe dimensions, and the
exact malformed-input diagnostic at the reducer's external source boundary.
The new Laurent control now requests four coefficients so that its unchanged
epsilon-squared assertion is supported by the backend truncation. Exact rational
arithmetic proves the CP snapshot's only remaining change is the spelling of
three simultaneous factor signs; that cosmetic refresh retains the count,
source IDs and total. It is separate from the ten bookkeeping proposals later
approved and validated as recorded above.

The subsequent `fd27769a` native package run passed 501 of 502 tests. Its only
failure was presentation: encoded index symbols hid the original edge/vertex
labels from Python's raw-index formatter. The existing display owners now
resolve individual index payloads while keeping native identities as keys.
Temporary presentation uses existing portable print-head tags and adds no
callbacks or live-tensor rewrites. Both new controls and all 171 Python/tensor
wrapper tests pass, including the original plain, LaTeX and Typst provenance
assertions. Those display fixes were included in the completed `880ad11d`
application/workspace run.

The final ownership audit found two remaining uses of fresh interface inference
after transformations of a retained tensor. The evaluator's scalar color
invariant rewrite now uses `with_checked_expression`. Final specialization
uses `map_representation_dimensions` for both `D` and the evanescent dimension,
then checks the remaining scalar rewrite. Nonzero reversed-axis and ranked-zero
controls cover the distinction. The native reducer likewise compares decoded
dimension metadata while retaining the original encoded slots; it does not
decode and re-admit the complete tensor. The runtime fermion spin-sum recipe
now allocates its dummy through `AbstractIndex`, avoiding a legacy function-valued
index inside an otherwise native projector. These repairs and the formal-color
comparison correction were included in the completed `880ad11d` run.

The formal-color probe establishes that `Nc` survives typed composition,
algebra and canonical preparation. The exact canceled quotient contains an
encoded scalar dimension for `Nc^2 - 1`; the scalar-name validator rejected
that symbol, causing the existing numerical fallback to return its value at
`Nc = 3`. The repair decodes dimension payloads only in this scalar quotient,
then applies the unchanged strict scalar check and the requested symbolic or
concrete-color option. Encoded index coefficients remain rejected. It neither
decodes the retained tensor bodies nor changes fallback sampling policy.

The first complete native run at `53312d8b` finished with 5,344 passes and 292
failures out of 5,636 selected tests. Most failures stopped at the same
encoded-dimension or legacy spin-sum index boundaries; their original numerical checks were not
reached. The remaining classes include the focused source/API cases above,
four generation snapshot tests (including the formal-color regression), six
source-dependent documentation failures, and one global tensor JSON round trip.
The latter wrote display-formatted payloads beside canonical interface metadata;
the serializer now uses canonical source strings and retains exact parse-back
checks. The later `880ad11d` full run kept its actual compile-time source tree
frozen and reduced the failure set to the five cases subsequently addressed by
the focused source/interface and approved-reference runs at the start of this record. The earlier diagnostic failures remain recorded as failures of
their own checkpoint.
