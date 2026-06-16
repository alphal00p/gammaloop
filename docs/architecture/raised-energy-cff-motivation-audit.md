# Raised-energy CFF change motivation and benchmarks

Current evidence, physical inputs, measured costs and review dispositions · 11 September 2026

## Assessment

The changes do not all have the same justification. Several repair concrete wrong values, index capture, invalid-input handling or output contracts. Others add requested capabilities or simplify ownership. The measurements support grouped dense/mixed contraction, no-chain collection, direct alias registration and exact base caching on their reached workloads. They also expose cost without demonstrated benefit on the tested inputs: CFF compaction and unproductive assignment contenders deserve reconsideration. Powered scalar-dot freshening has a reproduced correctness requirement at the whole-numerator FORM boundary. No universally best contraction ordering is established.

The inventory covers **73 logical categories, all 579 commit/path records and 505 distinct paths** in the seven reconstructed commits. A category is a semantic group, not an assertion that each changed hunk has a separate causal reproducer. The complete inventory includes source revisions, concrete inputs, current conclusions and a path-to-owner map. A successful current test alone is not evidence that the old implementation failed.

Benchmark measurements and the 579-record path map are pinned to `d55a0f3e9d25e5c64bdb9ef9fb57749ecf680390` (functional parent `1a747039`); the original tensor baseline is `39561014`. Most A/B experiments restore one old method or disable one optimization on those dependencies. They are not whole historical-release comparisons. Current closeout dispositions are identified separately: inspection JSON uses a pair-list adapter, the three compact-parser tests exercise public execution, and powered-dot freshening remains required. Final reconstructed-stack validation is supplied in the [stack review](raised-energy-cff-stack-review.md); these benchmark counts do not certify it.

## Inputs and measurement contract

The requested physical input is an isolated two-loop gluon two-point double triangle: a three-gluon vertex at 0, gluons 0–1 and 0–2, and a closed top-quark loop 1→2→3→1. Vertices use the SM `V_36` and three `V_137` rules, with the fermion-loop sign and complete propagator/vertex numerators retained. An external Lorentz/color trace closes the tensor. This is one physical diagram, not a complete gauge-invariant amplitude or integrated cross section. The 140,005-byte captured contraction input includes its complete CFF expression; UV and thresholds are disabled for that contraction benchmark.

The ttH capture is the physical GL38 NLO numerator (the source's GL20 label), with UV and threshold terms. BNL supplies four scalar aliases; nine legacy delta-head names are migrated to the current owner without changing arguments, coefficients or factorization. The factorized rank-four stress input supplies large scalar sums specifically to reach automatic deferred broadcasting. Controlled 4D/8D DD/DS/SD/SS tensors isolate storage effects. The CFF-only double-triangle-derived cases use owner-compatible energy-polynomial witnesses, not the full spin/color numerator or certified Taylor sectors. Quark-bubble Vakint inputs retain actual tensor and color algebra.

Root contraction, foundation, shared CFF and core variants have one excluded warm-up followed by six balanced forward/reverse rounds. Vakint and Feyngen body timing instead use six interleaved rounds with unequal pair-leading order. Feyngen end-to-end uses two exact forward/reverse rounds after a warm-up, a resource-based design revision declared before its measurements. Those limitations constrain small-difference claims. Rust builds use `dev-optim` (GammaLoop optimization level 2, dependency package overrides at 3). Executables, inputs and patches are hash-pinned. Root, foundation, core and Vakint timings use one worker on CPU 8; shared CFF uses CPU 9. Feyngen body timing uses CPU 8. Feyngen end-to-end uses the same CPU affinity 8–11 with either one or four workers; a one-worker process may migrate within that mask. A shared lock serializes the timing families and heavy comparisons; builds were coordinated to finish before measurements. Other users run on the host, an AMD EPYC 9754; there is no exclusive-machine or fixed-frequency guarantee. Small differences with overlapping ranges are descriptive observations.

Stage times exclude parsing/output when named. Whole-child RSS includes initialization, parsing, contraction and output, and can include a launcher high-water floor. It is not isolated kernel allocation. Binary hashing streams bytes; an invalid preliminary whole-binary-hash RSS experiment is excluded. Process-tree watchdog caps are 8GB for the main workflow families and 15GB for foundation probes; builds have separate caps within the coordinated 30GB aggregate budget. Failed observations remain visible and are not retried. Snapshot updates are disabled. Resource-guard termination is distinct from a returned wrong value.

## Tensor scheduling

Cells show median network-execution milliseconds / median whole-child RSS in decimal MB. All six scheduling choices use the same complete input, alias threshold and execution mode. SmallestDegree is the existing strategy selected by the diagnostic switch, even though its underlying CLI argument names the default preset.


| Strategy | Double triangle | ttH | BNL |
| --- | --- | --- | --- |
| intermediate-cost | 94.564 / 34.7 | 393.707 / 123.8 | 7.810 / 149.2 |
| sparse-atom-aware | 80.164 / 37.8 | 289.702 / 251.3 | 7.871 / 489.0 |
| atom-aware | 80.223 / 37.8 | 300.307 / 251.0 | 7.545 / 487.5 |
| result-rank-only | 89.591 / 50.4 | 3032.963 / 2699.8 | 42.700 / 3830.7 |
| entry-aware | 83.464 / 37.8 | 303.174 / 251.3 | 7.440 / 490.3 |
| smallest-degree | 170.297 / 50.4 | 8GB guard; no result | 8GB guard; no result |


Whole diagnostic-process medians in seconds include parsing, complete alias resolution and output dumping. They are not full GammaLoop generation or integration times.


| Strategy | Double triangle s | ttH s | BNL s |
| --- | --- | --- | --- |
| intermediate-cost | 0.200 | 0.675 | 0.793 |
| sparse-atom-aware | 0.191 | 0.804 | 1.993 |
| atom-aware | 0.207 | 0.836 | 2.036 |
| result-rank-only | 0.212 | 7.962 | 19.526 |
| entry-aware | 0.204 | 0.808 | 1.977 |
| smallest-degree | 0.298 | 8GB guard; no result | 8GB guard; no result |


Intermediate-cost has the lowest whole-process median on ttH and BNL despite slower network execution than some alternatives. In BNL, alias resolution alone takes 0.281s with intermediate-cost versus approximately 1.44–1.46s with the atom/entry alternatives and 18.81s with result-rank-only. Choosing from the network timer alone would miss this measured downstream cost.

Intermediate-cost uses less whole-child memory on ttH and BNL, with some network-execution cost relative to faster presets. These measurements do not separate live intermediate allocation from parsing and output. Result-rank-only is especially costly on ttH and BNL. SmallestDegree crosses the 8GB guard during warm-up on both, at 8.33GB and 8.45GB process-tree RSS; those runs have no valid completion time or full result. Its twelve later observations are deliberately unattempted. This compares the provided choices, not every possible contraction tree, cost-weight configuration or threading strategy.

The scheduling/kernel matrix records 177 observations: 150 measured successes, 25 successful warm-ups and two guard-terminated warm-ups. All ten distinct successful output hashes map to passed complete-value comparisons at three shared exact algebraic points. One output prints 78 instances of the exactly integral literal 2.00000000000000; an explicit 2.0→2 conversion ledger covers that representation boundary while retaining the unchanged oracle and original dump. Representation differences are not correctness failures. Raw distributions retain minima, maxima and sample standard deviations.

## Multi-index contraction and scalar aliases

Dense–dense and mixed multi-index contraction already worked through the generic API. The change extends grouped symbolic accumulation to those storage combinations. It does not add previously missing contraction semantics. The physical g–g trace reaches dense contractions over two, three and four common axes; both versions also pass the independent storage-component oracle.


| Input | Old network ms | Current network ms | Old / current |
| --- | --- | --- | --- |
| gg-double-triangle | 140.115 | 94.564 | 1.48× |
| ttH-GL38 | 488.379 | 393.707 | 1.24× |
| BNL-alias-pressure | 8.136 | 7.810 | 1.04× |


The strongest physical kernel gains are on gg and ttH. BNL has a smaller approximately 4% median difference; raw paired samples and spread remain in the evidence. These are network-stage comparisons, not equivalent total-workflow speedups.


| Controlled input | Old ms | Current ms | Old / current |
| --- | --- | --- | --- |
| d4-density1-DD | 0.4494 | 0.3138 | 1.43× |
| d4-density1-DS | 0.3954 | 0.2860 | 1.38× |
| d4-density1-SD | 0.3803 | 0.2903 | 1.31× |
| d4-density1-SS | 0.2843 | 0.2739 | 1.04× |
| d4-density16-DD | 0.2328 | 0.0100 | 23.35× |
| d4-density16-DS | 0.0252 | 0.0130 | 1.94× |
| d4-density16-SD | 0.0265 | 0.0133 | 2.00× |
| d4-density16-SS | 0.0135 | 0.0134 | 1.01× |
| d8-density1-DD | 8.1510 | 4.1105 | 1.98× |
| d8-density1-DS | 8.1740 | 4.2139 | 1.94× |
| d8-density1-SD | 8.2313 | 4.1515 | 1.98× |
| d8-density1-SS | 4.1439 | 4.1941 | 0.99× |
| d8-density16-DD | 3.9245 | 0.0870 | 45.13× |
| d8-density16-DS | 0.3519 | 0.1266 | 2.78× |
| d8-density16-SD | 0.3997 | 0.1526 | 2.62× |
| d8-density16-SS | 0.1959 | 0.1372 | 1.43× |


Each controlled tensor has two common Euclidean axes and one free axis per operand. Density16 means one populated coordinate in sixteen, including deliberately zero-filled dense storage. Every output component is compared with an independent coordinate sum using expanded subtraction on finite components. This does not distribute a graph numerator. The 32 old/current storage smokes and all 224 timed/warm-up storage observations agree. Sparse–sparse is mostly unchanged; neither the old sparse fiber nor the new loop is a nonzero-only traversal. These dimensions/densities do not establish asymptotic performance for enormous sparse domains.

On BNL, direct registration changes the alias-registration stage from **4.669 ms to 0.079 ms** (58.8×), with four aliases and equivalent complete results. That stage is a small part of total runtime; the total-wall distributions do not establish an end-to-end speedup. The double triangle and ttH create no aliases, so their old-registration flag timings are inactive controls.

## Chain collection, factorization and powered dots


| Operation/input | Old µs | Current µs | Old / current |
| --- | --- | --- | --- |
| chain_physical | 14483.235 | 83.504 | 173.44× |
| chain_scalar | 31.408 | 0.239 | 131.19× |
| prune_nozero | 7430.241 | 360.509 | 20.61× |
| gl06 | 195.927 | 229.500 | 0.85× |
| nested | 54.603 | 115.079 | 0.47× |
| square | 39.461 | 55.606 | 0.71× |


The no-chain early return removes work while preserving the complete physical numerator. The no-zero result times full canonicalization under a scoped multi-file ablation; it cannot apportion the improvement solely to pruning. Separately, the old zero-pruning path distributes an unrelated spectator, violating the explicit factorization contract even though its expanded algebra is equal.

Powered scalar-dot freshening is required by the whole-numerator FORM route. For `(k·p + k·q)^2`, the correct vacuum average is `k^2*(p^2 + 2*p·q + q^2)/D`, with `D=4-2*eps`. Disabling only the final freshening pass gives `k^2*(p^2+q^2) + 2*k^2*p·q/D`: both diagonal terms lose `1/D`. A four-case public-API comparison gives correct complete values for both current modes and unfreshened projected mode; unfreshened whole-numerator mode is wrong. The [powered-dot follow-up](raised-energy-cff-powered-dot-freshening.md) retains the actual backend inputs and independent exact residual checks.

The square, nested and GL06 Rust round trips remain valid controls. They do not establish backend equivalence: the Rust Lorentz parser contracts a scalar base before its exponent, whereas FORM receives indexed copies. Freshening costs more on the measured conversions, but this counterexample justifies retaining it for correctness. No manual PySecDec integration is claimed.

## Deferred scalar-weight broadcasting


| Variant | Network ms | Whole wall s | Child RSS MB |
| --- | --- | --- | --- |
| deferred | 2.943 | 3.337 | 18.9 |
| old-eager | 465.080 | 4.409 | 550.4 |
| aliases | 1.560 | 3.758 | 189.3 |
| old-alias-registration | 1.566 | 3.676 | 189.8 |


The controlled factorized input is `(c0 q0⊗q0⊗q0⊗q0 + c1 q1⊗q1⊗q1⊗q1) : q2⊗q2⊗q2⊗q2`, where each coefficient sums 12,000 simple rational terms. The independent value is `c0(q0·q2)^4 + c1(q1·q2)^4`. Python Fraction arithmetic checks three complete values, including full alias definitions. It is a representative scalar-broadcast stress input, not a physical UV-locality certificate or an identity proof for all scalar parameters. Current versus old-eager with aliases disabled isolates automatic deferral; enabling aliases is a separate intervention. None of the ordinary physical captures reaches automatic deferral.

## Shared CFF generation

Generation medians in milliseconds; complete output serialization is outside the timer. The complete distributions and all five variants are retained in the evidence bundle.


| Input | Current | Eager layers | Eager fanout | No base cache | No compaction |
| --- | --- | --- | --- | --- | --- |
| gg affine | 0.993 | 1.001 | 1.002 | 1.003 | 1.002 |
| gg dotted | 6.831 | 6.926 | 6.897 | 7.024 | 6.996 |
| cubic box | 53.243 | 52.968 | 54.906 | 109.092 | 60.488 |
| quartic sunrise | 20.668 | 19.777 | 19.304 | 34.098 | 18.129 |
| two bubbles | 1.141 | 1.118 | 1.183 | 1.119 | 1.085 |
| gg rank-6 stress | 1433.976 | 1298.905 | 1297.506 | 2376.280 | 1156.223 |


All 210 observations complete (180 measured, 30 warm-ups), and every dumped expression matches its fully evaluated trace for the same variant. Exact base caching helps the reached box, sunrise and larger gg topology. On the larger gg diagnostic, disabling compaction is faster in every paired round; mean generation falls from 1.475s to 1.154s, with the same 22.4MB expression and approximately 44MB peak child RSS. This is a concrete simplification candidate for this path, not proof that compaction is unnecessary for products of multiple components.

Streaming reaches only one-component products in that larger stress case; the two-component bubble is tiny. Neither eager/streaming comparison demonstrates a resident-memory saving here. Source evidence records a GL16 30GB failure, but these isolated experiments do not reproduce or apportion that resource failure. No isolated Rational-versus-Atom or merged-versus-independent-occurrence benchmark is claimed.

The complete CFF checks include 285 numerical observations and 15 independent denominator-pinch oracles, with maximum pinch error 1.507e-15. The 28 fresh public contracts pass. Exact old-parser-method replay distinguishes `-2**2` (4→−4), `2**3**2` (64→512), and an oversized exponent (wrapped value→error). These are small public-API counterexamples, not physical benchmark inputs.

## Core UV caches and assignment scoring

The core experiments keep physical numerator owners, request-local capacities, exact reconstruction certificates and selected payload/assignment association intact. Total-cache disabling removes both payload and count reuse; count-only disabling retains completed winners. The one-proposal variant scores the first certified rank-ordered candidate but still constructs the proposal list. It is not a complete recreation of the old planner or union-envelope cache.


| Input / variant | Wall mean ± SD s | Wall median s | Peak child MB |
| --- | --- | --- | --- |
| GL04 / current | 15.601 ± 0.291 | 15.622 | 55.0 |
| GL04 / no-count-memo | 15.801 ± 0.392 | 15.643 | 54.3 |
| GL04 / no-exact-cache | 15.729 ± 0.576 | 15.552 | 54.5 |
| GL04 / one-proposal | 15.438 ± 0.394 | 15.416 | 54.4 |
| GL04-temporal-square / current | 5.746 ± 0.285 | 5.802 | 44.3 |
| GL04-temporal-square / no-count-memo | 5.766 ± 0.300 | 5.793 | 44.3 |
| GL04-temporal-square / no-exact-cache | 5.778 ± 0.313 | 5.829 | 44.3 |
| GL04-temporal-square / one-proposal | 5.615 ± 0.515 | 5.746 | 44.4 |


All 56 observations pass (48 measured, 8 warm-ups). The measured unit includes fresh import/generation, three inspections and computed forest export, including work repeated by that exporter. Mean differences of 0.3–2.3% are small relative to the 0.28–0.58s sample spread; these runs do not establish a convincing total-workflow speed or memory benefit.

Both GL04 inputs have exact equality of all three complete inspection JSON values and all nine computed forest/node DOT exports across variants. Their small real values are compared directly, without a unit-sized absolute tolerance floor. All alternative proposals lose: eleven GL04 groups have map counts `[14,20,20]`; seven temporal-square groups have `[32,56,44]`. Count-only disabling retains the same hit counts as current code, so these inputs do not demonstrate losing-count reuse. They exercise scoring overhead but do not justify its benefit. A partial GL16 trace does reach a winning third proposal: `[1200,1224,1104]` reduces the rank-first map count by 8%. All four GL16 variants time out at 180s before inspection or computed exports; none supplies full-value parity or a total-runtime comparison. The traced peak process-tree RSS is 545MB, so this is a time limit, not a reproduction of the source's 30GB failure. Matching diagnostic prefixes show far fewer generated contenders with memoization, and 52 paired parsed-source/bounds records match exactly. Cache hits omit the source and some canonical-key fields, so the traces are not an independent proof of every hit's key identity.

## Vakint and numerator matching

Both Vakint input modes and both forest owners produce the same complete nonzero Laurent expression for the physical gluon–top-quark bubble. Exact dummy-name normalization and rational cancellation check pole and finite terms; mode flags are verified at the actual backend boundary. The scalar self-energy is a control. All 56 measured/warm-up workflows complete and their full computed exports match the validated results. Scalar medians are 0.783–0.809s; quark-bubble medians are 1.382–1.405s, with 35MB/69MB peak child RSS. Paired monolithic/projected ratios are about 0.990 and 1.003 for the two owners, with overlapping ranges: no consistent speed advantage is established. Projected mode leads three of six paired rounds for quark/hedge and four of six for the other mode pairs; the schedule is interleaved, not exactly balanced.

For the full two-loop double triangle, projected processing exceeds 180 seconds for both owners, while whole-numerator processing reaches the 8GB guard for both. These are explicit coverage limits, not failed numerical comparisons. Both modes still require d-dimensional closure and adequate Laurent depth; choosing an input mode does not replace those correctness requirements.

The numerator-matching input is a two-loop three-gluon vertex family with top-quark loops: 36 candidates in 24 topologies. Twelve nonsingleton buckets contain two candidates each; all twelve compared pairs are rejected. The inputs carry five numerical and five polynomial samples, with no canonical numerator. Twenty repeated calls use each actual physical pair. This reaches early rejection but neither long bucket scans nor the canonical-numerator path. The smaller physical control additionally exercises a matching pair with complete result `Some(1)`.

The body experiment completes all 28 observations (24 measured, four warm-ups). Every pair key is unique; all complete comparison results agree between old and current methods, and every 36-graph DOT map agrees with the validated physical reference. Times below are medians per comparison call, averaged across the twelve physical pairs within each process. The timer includes black-box argument/result handling and a push into a preallocated result vector; setup, allocation and equality assertions are outside it. These microbenchmarks repeat the comparator, not the full generator.


| Comparison | Old ms/call | Current ms/call | Paired old/current ratio [range] |
| --- | --- | --- | --- |
| sign | 2.510 | 0.355 | 7.074 [6.875, 7.349] |
| scalar | 172.269 | 171.988 | 0.999 [0.984, 1.009] |


Ratios are computed within each paired round before taking their median; they need not equal the ratio of the two marginal time medians.

Sign-only rejection is faster on this reached workload; scalar-factor comparison shows no benefit beyond run-to-run variation. The six body rounds are interleaved with current first in four sign rounds and five scalar rounds, so small differences must not be presented as precise causal effects. Repetition-inclusive CLI medians are about 67s for sign and 110s for scalar: most of that work is numerator preparation, and twenty artificial comparator repetitions must not be confused with default generation cost.

The global-lock intervention retains current bucket storage, clone removal and zero-mask construction. It restores serialization around lookup/search/insertion rather than rebuilding the complete old architecture. The larger four-gluon case reaches 567 candidates but times out during preparation/grouping at 180 seconds; it supplies no completed graph parity or timing ratio. A 5,000-repetition calibration also times out, which is a limit of the artificially amplified body experiment, not a demonstrated timeout of default generation.

The end-to-end experiment accounts for 48 executed observations out of 48 planned: 48 successful, 32 measured, 0 failed and 0 omitted after a failed variant. Successful observations retain all 36 physical graph exports and agree directly with the validated default. Each variant has one excluded warm-up and two exact forward/reverse measured rounds; the repetition hook is disabled.

The cells below contain the range of the two complete CLI wall times in seconds. Preparation, lookup and grouping remain included; these are descriptive observations on a shared host.


| Mode / workers | Current s | Old comparator s | Global lock s | Both s |
| --- | --- | --- | --- | --- |
| sign / 1 workers | 66.002–66.507 | 66.206–67.265 | 66.208–66.614 | 66.516–66.776 |
| sign / 4 workers | 19.048–20.435 | 18.634–18.810 | 18.146–18.170 | 17.737–18.867 |
| scalar / 1 workers | 68.633–68.913 | 68.789–68.800 | 68.526–69.190 | 69.236–69.276 |
| scalar / 4 workers | 19.641–19.725 | 18.622–18.761 | 18.560–18.595 | 19.005–19.526 |



| Mode / workers | Old/current ratio | Global/current ratio | Both/current ratio |
| --- | --- | --- | --- |
| sign / 1 workers | 0.995–1.019 | 1.002–1.003 | 1.004–1.008 |
| sign / 4 workers | 0.920–0.978 | 0.888–0.954 | 0.868–0.991 |
| scalar / 1 workers | 0.998–1.002 | 0.994–1.008 | 1.005–1.009 |
| scalar / 4 workers | 0.948–0.951 | 0.943–0.945 | 0.968–0.990 |


No end-to-end gain is established. In the four-worker observations, current takes longer than each ablation in both rounds; these results must not be hidden by the faster isolated sign-comparison body. The small shared-host sample does not establish the cause or a general regression.

The numerator-preparation/grouping milestones, whole-child RSS and full samples are retained in the evidence. The milestone is an inclusive preparation/grouping duration, not isolated lock time. Two observations per variant cannot establish a precise speedup; worker scaling mostly reflects the complete parallel preparation workload. With at most two candidates per bucket, this experiment cannot certify the motivating large-bucket/global-lock claim.

## Reproduced defects and current dispositions

On the pinned benchmark source, the focused foundation suite has 22/22 passes. Reversing eight production files on current dependencies gives 12 passes and 10 failures: seven index/normalization/namespace/factorization failures and three newly supported compact-syntax probes. These are not ten independent physics bugs. The baseline regression build additionally enables Idenso serialization/reference-case features; the named arithmetic paths are unchanged, but the full build feature sets are not identical. The nested-sum Symbolica panic is a pinned-source counterexample; the fixed dependency remains in this local ablation. No failing product test is edited.

**Inspection JSON preserves structured additional-weight identifiers.** `GenericAdditionalWeightInfo.weights` uses `#[serde(with = "vectorize", bound(deserialize = "T: Deserialize<'de>"))]`, reusing the existing adapter to serialize `[key, value]` pairs, including `[]` for no weights. The explicit bound restores Serde's original deserialization requirement without strengthening the generic type. The motivating saved-GL04 probe isolates JSON serialization after successful numerical evaluation: the benchmark source rejects a structured threshold-counterterm map key with `key must be a string`. Current outward tests cover all four key variants and empty weights through JSON/bincode round trips, plus physical inspection with weight retention both enabled and disabled. They compare retained weights, event grouping and evaluation values; the CLI comparison does not claim every event metadata field. Final execution evidence is supplied in the [stack review](raised-energy-cff-stack-review.md).

A diagnostic boundary also needs care: reparsing Python Symbolica's canonical string for `(2+3i)/7*x` changes coefficient grouping. The plain printer round-trips in the reproducer. This is not evidence of a GammaLoop evaluation defect; the comparison harness retains live expressions and parses original Rust plain-string dumps, never intermediate canonical strings. The exact reproducer is retained.

The three compact-parser regressions in [parsing/test.rs](../../crates/spenso/src/network/parsing/test.rs) use public parsing, execution and complete scalar results. The weighted-vector test compares against an independent four-component Minkowski sum with signature `(+,-,-,-)` for weights `-2`, `a` and `a+b`. The two ambiguous arguments—two compact vectors, or a compact vector times an explicit tensor—retain their complete opaque function values. Exact expanded zero residuals check these values without fixing intermediate products, allocated slots or traversal order. Separate external-slot contracts remain. Exact factorization assertions remain appropriate where preservation of factorization is itself the intended public behavior.

## Rust patterns and KISS

| Pattern | Effectiveness supported by the audit |
| --- | --- |
| Shared `Rc<Cell<_>>` allocation plus explicit reservations | Shared allocation is established behavior; reservations fix a reproduced capture. Keep ownership local to one parse and do not claim this introduces sharing. |
| Borrowed input views and grouped iterator accumulation | Avoid repeated symbolic normalization/copying; dense/mixed measurements support the specialization. General sparse-coordinate traversal remains a limit. |
| Owned scalar aliases and direct registration | The alias already exists, so registering it avoids a redundant root replacement search. The four-alias BNL stage supports this simplification. |
| Explicit deferred-result alternatives | Preserve a factorized sum until terminal materialization when eager broadcasting is expensive. The added branch must earn its complexity on an activated workload. |
| Concrete `Rational` coefficients and boundary conversions | Excludes unsupported symbolic coefficient states and makes exact arithmetic intent clear. Previous rational Atoms were already exact; no rounding or speed claim follows from the type alone. |
| One expression assembler and fallible public conversions | Shared invariant ownership reduces duplicate construction logic. Small malformed-input counterexamples justify Result-based rejection before indexing. |
| Immutable assignment plans and typed occurrence/owner IDs | Keep capacity, routing, physical owner and selected expression together. Reconstruction certificates justify these distinctions; LMB coordinates do not replace physical ownership. |
| Exact-key caches and bounded heuristic search | Base caching pays on reached high-rank inputs. No benefit from losing-count memoization or extra contenders is demonstrated on the completed GL04 inputs; a smaller map is only a proxy for total cost. |
| Mode enums and builder-owned settings | Clarify supported alternatives and keep migration at an owner. Both Vakint modes are requested functionality; no measured universal advantage justifies presenting one as inherently superior. |
| Short bucket locks and staged exact/sample comparisons | The reached sign-rejection body is faster; scalar matching has no demonstrated gain. The physical family has only two-candidate buckets, so the large-bucket motivation remains unproved. Canonical-numerator matching is not reached, and sampled matches remain probabilistic. |
| Serde field adapter with an explicit deserialization bound | Pair-list serialization preserves structured identifiers without string parsing. Reusing vectorize keeps the fix at the owning field; its explicit bound retains the original generic deserialization requirement. |

Keep the demonstrated correctness fixes and required feature boundaries. Treat combined root/product compaction on the reached single-component path, unproductive proposal scoring and unproved large-domain sparse behavior as review findings with explicit evidence limits. Do not remove mathematically necessary ownership distinctions merely to reduce line count, or retain a performance mechanism only because current tests pass.

## Per-change motivation inventory

Each row states the current justification or measurement limit. The machine-readable inventory retains pinned benchmark evidence and the original path ownership, with separate `closeout_disposition` fields for the completed JSON/test changes and retained powered-dot requirement. Recorded counterexamples identify a concrete failure in retained source evidence; they are not relabelled as fresh executions. Capability, policy and consumer migrations need not invent a motivating bug.

## Commit 1: motivation rows

| Change | Evidence and conclusion |
| --- | --- |
| C1-01 — Reserve explicit compact-parser dummy names | Reproduced index capture: a compact contraction collides with an explicit Dummy(1000000). The source also records the physical ttH overcontraction. Reservation is the new behavior; shared Rc allocation already existed. |
| C1-02 — Reserve external and dual indices during canonical dummy naming | Reproduced external/dual index capture during canonical relabeling. Preserve external index identity; no timing claim is needed. |
| C1-03 — Ignore canonical names from canceled representation groups | Reproduced canonical-identity failure after one representation cancels. A canceled group must not consume names used by the surviving contraction. |
| C1-04 — Preserve completed contractions across nested tensor sums in pinned Symbolica | Pinned source records a nested-sum contraction-counter panic and a reduced expression. Current probes pass; the old Symbolica dependency itself was not rebuilt in this audit. |
| C1-05 — Prune antisymmetric zeros locally without distributing spectators | Reproduced unwanted distribution of a spectator factor. The no-zero fixture also runs faster under current full canonicalization; that timing does not isolate pruning from all naming changes. |
| C1-06 — Skip chain collection on expressions with no chains | Measured benefit: the captured physical numerator has no remaining chains. Complete output is unchanged, and the early return removes almost all collection work. |
| C1-07 — Allow one implicit compact axis beside explicit spectator slots | Supported syntax extension: one implicit axis can coexist with explicit spectator slots. Public probes check external indices; this is not evidence that dense/mixed contraction was broken. |
| C1-08 — Allow scalar-weighted compact vectors while retaining ambiguity checks | Supported syntax extension for scalar-weighted vectors, useful for momentum differences in a 3g vertex. Public parse/execute regressions now compare an independent weighted Minkowski sum and complete opaque values for ambiguous arguments. |
| C1-09 — Freshen independent copies of powered scalar dots | Reproduced correctness requirement: disabling only freshening breaks the whole-numerator FORM angular average of (k·p+k·q)^2, dropping 1/D from both diagonal terms. Both current modes and unfreshened projected mode match the complete independent oracle. Rust-only round trips miss this backend failure; the measured conversion cost does not justify removal. |
| C1-10 — Use numeric imaginary coefficients in Vakint normalization | Reproduced wrong symbolic normalization when a user registers the imaginary-name symbol. A numeric complex coefficient preserves the intended value independently of the symbol registry. |
| C1-11 — Preserve canonical user namespaces across FORM | Reproduced FORM parsing failure with canonical user namespaces. Complete tensor-reduction and adapter probes pass with the fix. |
| C1-12 — Generalize grouped Atom contraction across tensor storage kinds | Measured optimization. Old dense/mixed contractions already return correct values. Generalized grouping speeds the physical captures and most controlled storage cases; sparse–sparse gains are limited. |
| C1-13 — Score contraction candidates by intermediate symbolic cost | Measured tradeoff across all five exposed presets plus SmallestDegree. Intermediate-cost has the lowest whole diagnostic time and memory on ttH/BNL, although other presets have faster network steps. It is not a universal optimum. SmallestDegree crosses 8GB on two captures. |
| C1-14 — Register already-created scalar aliases without recompressing the root | Measured simplification on the four-alias BNL capture: direct registration avoids an unnecessary root search. The physical double triangle and ttH create no aliases and are inactive controls. |
| C1-15 — Defer large coefficient broadcasting while preserving ordinary terminal results | Measured on a factorized rank-four scalar-weighted tensor that reaches automatic deferral. The three ordinary captures do not reach it. Complete Fraction component oracles check both routes. |
| C1-16 — Select clang for macOS Rust unwinding | Pinned source records the macOS linker/panic-unwinding problem. Linux execution cannot reproduce that platform failure; no cross-platform speed claim is made. |
| C1-17 — Guard heavy process trees and keep local runtime outputs out of version control | Resource policy and tooling support. The audit exercises the process-tree guard and retains its terminations. Ignore patterns have no runtime performance justification. |
| C1-18 — Keep manual PySecDec integrations outside automated selections and require supported nextest | Explicit test-selection policy and supported nextest requirement. Manual PySecDec integrations remain outside this automatic audit; this is not a numerical algorithm change. |
| C1-19 — Pin compatible Python Symbolica and UFO loader | Pinned source records Python Symbolica import/fixed-point trouble and a UFO loader requirement. Dependency compatibility is the motivation; no independent physics speedup is claimed. |
| C1-20 — Regenerate shared dependency feature metadata | Required generated dependency metadata. It keeps feature resolution/builds consistent with the owning API changes; it is not an independent optimization. |
| C1-21 — Remove excess trait bounds and migrate sum-capable consumers | API simplification and necessary consumers of sum-capable terminal results. Removing unused bounds narrows conceptual requirements; no separate wrong-value bug is established. |
| C1-22 — Preserve physical numerator factors in diagnostic consumers | Required factorization-preserving diagnostics and outward assertions. Tests must compare complete values or an explicit factorization contract, not formatting or intermediate storage. |
| C1-23 — Regenerate current optional defaults and schema formatting | Generated schema/default alignment. Acceptance is agreement with serialized public settings; no contraction benchmark applies. |
| C1-24 — Clarify the scope of older large-expression observations | Documentation accuracy. Earlier memory observations do not certify a current speedup; this audit replaces broad claims with workload-specific evidence. |


## Commit 2: motivation rows

| Change | Evidence and conclusion |
| --- | --- |
| C2-01 — Shared generalized CFF engine and selectable representation | Requested raised-power capability. Shared CFF contracts and complete signed/pinched values provide acceptance; an unsupported parent input is not an old wrong-result bug. |
| C2-02 — Retain independent repeated propagator occurrences | Independent occurrence ownership is covered by current complete-value probes. No isolated merged-versus-independent performance comparison or old wrong-value counterexample is established here. |
| C2-03 — Stream component products and shorten intermediate lifetimes | Source records a GL16 30GB resource failure. Current streaming/eager comparisons show no speed or RSS advantage on the reached workloads; the larger stress products have only one component. |
| C2-04 — Reuse exact bases and compact completed products | Exact base caching helps the reached high-rank cases. Compaction adds cost without reducing the larger double-triangle-derived output; necessity for other multi-component inputs remains open. |
| C2-05 — Correct shared terminal/embedded contour normalization | Source records an embedded-versus-standalone signed contour mismatch exposed by the unchanged depth-three banana. Correct common normalization is a mathematical requirement, independent of speed. |
| C2-06 — Correct diagnostic numerator arithmetic | Fresh exact old-method replay: `-2**2` gives 4 instead of -4, `2**3**2` gives 64 instead of 512, and an oversized exponent wraps. Current public contracts reject/correct these. |
| C2-07 — Validate bounds before fast paths and reject malformed graph signatures | Source records accepted unknown bound 999 and a malformed-signature panic. Fresh current public probes reject both; this audit does not rebuild every old validation branch. |
| C2-08 — Additional input arithmetic/shape validation | Fresh current arithmetic, shape, scale and graph-balance contracts pass. Some guards address source-level boundary counterexamples rather than an observed physical failure. |
| C2-09 — Rational coefficient storage across shared engine and consumers | Rational storage expresses the coefficient domain and removes unsupported symbolic states. Rational-valued Atoms were already exact. No rounding bug or isolated storage speedup is claimed. |
| C2-10 — Shared expression assembly, surface copying and builder consolidation | KISS consolidation: one assembler owns interning, copying and finalization. Complete-value coverage supports behavior; there is no standalone measured benefit or old bug for each moved line. |
| C2-11 — Preserve zero values for filtered and infinite denominator trees | Concrete numeric/symbolic zero-contract mismatch: absent/Infinite denominator trees could panic or evaluate near 1e80. Fresh public zero/fusion checks pass. |


## Commit 3: motivation rows

| Change | Evidence and conclusion |
| --- | --- |
| C3-01 — Integrate shared raised-energy CFF into GammaLoop and evaluator persistence | Necessary GammaLoop integration of the new shared engine, with generated evaluators and persistence. Physical captures and exact scalar/source certificates exercise the boundary. |
| C3-02 — Exact owner-preserving local-4D UV source reconstruction | Source contains complete GL04 numerator and denominator reconstruction certificates isolating owner/routing errors. Equal denominator squares alone cannot certify an odd numerator sign. |
| C3-03 — Typed direct/local/projected UV routes and explicit selector-free sums | Requested direct/projected routes plus removal of stale state. Explicit selector truth tables are the behavioral contract; API cleanup has no separate performance claim. |
| C3-04 — Owner-specific integrated UV mass scaling | Source records the GL04 finite-localizer coefficient discrepancy and enclosing/disjoint controls. Owner-specific vacuum-mass scaling is a correctness fix, not a global mass rescaling. |
| C3-05 — Keep natural Taylor denominator topologies and factorized numerator ownership | Source records expensive common-denominator reconstruction and the natural Taylor-topology requirement. This audit does not isolate old-versus-current topology collection performance. |
| C3-06 — Factorized degree analysis and independent outer rank envelopes | Source records slow factorized degree/rank-envelope analysis. No isolated timing of old degree-analysis methods is established here; combined generation timing cannot apportion its benefit. |
| C3-07 — Bound exact-CFF cache keys to each request capacity | Request-local capacity prevents union envelopes from enlarging unrelated requests. Core cache-removal measurements keep capacities fixed, so they do not benchmark the old union-envelope policy. |
| C3-08 — Choose up to three certified derivative-energy assignments by generated map size | Core experiment compares complete scoring against the first certified proposal, including downstream cost. It retains proposal planning and certificates; fewer generated rows alone are not a speedup certificate. |
| C3-09 — Route child-soft carriers into parent hard chart | Source identifies GL08 retaining an unintegrated child momentum. Compatible parent-chart routing is the recorded fix; no fresh reduced old/current GL08 replay is claimed here. |
| C3-10 — Resolve symbolic constants before SymJIT compilation | Pinned source records unresolved pi^-3 producing infinity in SymJIT while the mapped evaluator is finite. Resolve constants at the compiler boundary; no speed motivation is needed. |
| C3-11 — Validate source integer arithmetic and raised residue selections | Concrete invalid-input and overflow boundary cases motivate fallible selection/conversion. No typical physical graph reaching the integer limits is claimed. |
| C3-12 — Uniform nonzero M evaluator parameter and finalized expression ownership | Uniform nonzero M and finalized expression ownership simplify the evaluator contract. This is an explicit API choice, not a demonstrated rounding fix or isolated optimization. |


## Commit 4: motivation rows

| Change | Evidence and conclusion |
| --- | --- |
| C4-01 — Command templates and invocation scope | Source contains real CLI old-fail controls for nested inline/boot template scope. Lexical inheritance and local overrides are observable behavior; graph timing is irrelevant. |
| C4-02 — Explicit output, state detection/versioning and history persistence | Source reproduces unsaved scalar-box 3D output writing into the state directory. Explicit write/version boundaries protect persistence. A Serde pair-list adapter now preserves structured additional-weight keys in inspect JSON; round-trip and physical inspection tests cover the public output contract. |
| C4-03 — 3D representation build/validate/export CLI | Requested 3D inspection/export workflow. Complete exported-function evaluations are the acceptance criterion; output paths and saved/reloaded functions are public contracts. |
| C4-04 — Benchmark command, timing reports and temporary minimal evaluation settings | Requested benchmark observability and temporary-settings restoration. Existing public success/error restoration checks support it; no old settings-loss reproduction or large-batch memory bound is claimed. |
| C4-05 — UV profiling and stable finite-range numerical reports | Concrete extreme-value counterexamples motivate nonzero arbitrary-precision bounds and finite UV sample validation. Small physical integrands cannot exercise every reporting exponent boundary. |
| C4-06 — Numerical-stability histograms, dashboard and integration checkpoint outputs | New statistics/checkpoint workflow and required consumers. Public batch, persistence and restoration contracts motivate it; no separate physical bug or reporting-overhead improvement is established. |
| C4-07 — Resource-aware test execution, Nix inputs and tracing | Recorded resource/build issues and explicit tooling policy. Platform-specific memory accounting is not inferred from Linux performance measurements. |
| C4-08 — Behavioral numerical tests and archive unconsumed scalar snapshots | Test-detection defects are concrete: scale-dependent tolerances accept a factor-two error, nonfinite predicates can pass, and unused snapshots assert nothing. Independent value oracles are the remedy. |
| C4-09 — Remove obsolete wrapped model and migrate examples/metadata | Cleanup of an unreferenced wrapped model plus necessary examples/schemas. Reference and live-model audits justify removal; no independent speedup is claimed. |
| C4-10 — Preserve precise stability digits before f64 narrowing and report inclusive bounds | Concrete precise-output errors: narrowing 1e-400 to f64 zero loses information; an inclusive overflow bin must not be labelled with a strict bound. Current boundary tests cover the contracts. |
| C4-11 — Process import/reference contracts and named-integrand provenance | Public import/provenance contracts. The audit itself encountered inferred-all-cuts versus explicit-process inputs, confirming why the process specification matters; this is not a core CFF mismatch. |


## Commit 5: motivation rows

| Change | Evidence and conclusion |
| --- | --- |
| C5-01 — Physical amplitude/LU and threshold phase normalization | Independent signed scalar/left-right contour counterexamples establish the phase requirement. A common wrong phase can survive strategy-equality checks; no speedup is needed. |
| C5-02 — Inverse-process sewing, complex couplings and explicit CP opt-in | Complex charged-current and scalar spin-matrix oracles reject the proposed extra-adjoint route; this is not a fresh reproduction of an original-main defect. Real-coupling gluon graphs cannot certify this contract. CP optimization remains an explicit assumption. |
| C5-03 — Consistent W/Z virtual propagators, covariant cut multiplets and physical event labels | Full covariant versus physical polarization sums and Ward identities motivate gauge sewing/model consistency. QCD-only benchmarks do not establish the electroweak contract. |
| C5-04 — Reject unsupported complex masses/custom pole rules and use ghost anticommutation | Unsupported real-energy/complex-mass and custom-pole boundaries require validation. Ghost loop/exchange signs share the anticommutation predicate; this is a supported-domain contract. |
| C5-05 — Feyngen accepted-numerator lookup optimization | Physical repeated-comparator measurements show a 7.07× median sign-rejection benefit and no scalar-matching benefit on twelve rejecting pairs. All 48 end-to-end observations preserve complete exports, but no workflow gain is established; current four-worker runs are slower than all ablations in both measured rounds. With only two rounds and two-candidate buckets, neither a precise speed effect nor long-bucket scaling is established. Canonical-numerator matching remains unexercised. The global mutex is a serialization ablation, not the complete old map/clone architecture. |
| C5-06 — Retain fixed-momentum internal bridges during initial-state subtraction | Recorded reproducible GL3 bridge bug: an internal fixed-momentum edge touches neither initial-cut endpoint, yet the old helper removes a vertex and destroys cut-side loop rank. |


## Commit 6: motivation rows

| Change | Evidence and conclusion |
| --- | --- |
| C6-01 — Projected tensor kernels/scalar coefficient preservation and whole-numerator alternative | Both Vakint input modes and forest owners agree on complete nonzero quark-bubble Laurent output. Timings show no consistent mode advantage there; full double-triangle UV exceeds declared limits. |
| C6-02 — Complete d-dimensional Dirac/metric closure before Laurent truncation | Independent d-dimensional Dirac/metric and scalar-master oracles expose finite evanescent terms lost by premature 4D specialization. Both input modes need the same closure boundary. |
| C6-03 — Account for coefficient/normalization epsilon poles and reject insufficient depth | Laurent order budgeting must add numerator-coefficient and normalization pole orders. Exact single/double-pole counterexamples justify the change without a performance claim. |
| C6-04 — Validate Lorentz domains, opaque tensor slots and scoped gamma5 errors | Source records odd-open/renamed-dummy boundary failures; gamma5 rejection is a scoped unsupported-domain policy. Cancellation cannot hide invalid tensor slots before validation. |
| C6-05 — Remove exactly-zero factorized scalar coefficient subtrees | Source records GL29/35/38/40/46 generation failures and a concrete exactly-zero scalar coefficient. Removing that zero while retaining nonzero spectator factors is the value/factorization contract. |
| C6-06 — Canonical backend serialization and post-restoration notation | Source gives corrupt FORM input (1+2i)*(a+b)*tensor losing parentheses, canonical-name corruption and invalid decimal-complex tokens. Adapter correctness, not speed, motivates the fixes. |
| C6-07 — Preserve signed numerical momenta, complex parameters and combined estimator | Signed real momenta, invalid complex inputs and transactional settings have concrete adapter contracts. Manual PySecDec integration evidence is not fresh automatic coverage in this audit. |
| C6-08 — Normalize completed analytic return before reinsertion and migrate forest owners | Source gives a sunrise analytic-return normalization discrepancy with equal exact Laurent algebra. Both forest owners and modes require the same reinsertion normalization; not a new integral value. |


## Commit 7: motivation rows

| Change | Evidence and conclusion |
| --- | --- |
| C7-01 — Current reports, Typst/PDF and durable source/validation provenance | Requested documentation and provenance. The report records actual commands, complete-result coverage, measured costs, current dispositions and unsupported motivation claims. No runtime speedup applies. |


## Pinned benchmark validation evidence

Counts describe the accepted command families in this audit. Warm-ups are retained separately from measured observations; repeated calls within one process are not independent test executions. The report does not combine these heterogeneous checks into a headline test total.

| Family | Fresh result | Complete-result evidence |
| --- | --- | --- |
| Foundation regressions | Current 22/22; scoped baseline 12 passes and 10 observed distinguishing failures | Seven defect/factorization cases and three supported-syntax probes; unchanged product tests |
| Foundation timing and storage | 308 successful observations: 264 measured, 44 warm-ups; 32 storage smokes pass | Independent component/value oracles outside the timed operation |
| Physical contraction matrix | 175 successful observations: 150 measured, 25 warm-ups; two guard terminations | All ten successful output hashes checked at three exact shared algebraic points |
| Deferred scalar-weight stress | 28 successful observations: 24 measured, four warm-ups | Three distinct full output hashes checked with the independent Fraction oracle |
| Shared CFF | 210/210 observations: 180 measured, 30 warm-ups; 28/28 public contracts | Full-expression trace linkage, 285 numerical observations and 15 independent pinch oracles |
| Core GL04 reconstruction | 56/56 observations: 48 measured, eight warm-ups; four GL16 timeouts | Three full inspections and nine complete computed forest/node exports |
| Vakint modes and owners | 56/56 observations: 48 measured, eight warm-ups; four full two-loop limits | Full nonzero Laurent oracle and direct linkage of every timed export |

Feyngen's body and end-to-end counts are reported with their timing tables and actual schedule. The reproduction bundle records formatting, locked checks and builds of the diagnostic executables, exact compiler/dependency hashes, commands and failures. These counts belong to the pinned benchmark sources and isolated interventions. Final reconstructed-stack commands and results are supplied in the [stack review](raised-energy-cff-stack-review.md), without reusing benchmark or earlier workspace-suite totals as fresh validation.

## Evidence and reproduction

The [complete inventory](raised-energy-cff-motivation-inventory.json) and [evidence bundle](raised-energy-cff-motivation-evidence.tar.gz) accompany this report. The bundle contains source pins, exact intervention patches, public-API probes, physical graph/model captures, input hashes, full timing manifests, per-output comparison receipts, resource failures and independent methodology review. It excludes executable binaries and license material. The reproducibility guide distinguishes valid final observations from discarded setup measurements and states how to provide a licensed local environment.

The validation counts in this audit come from its own pinned commands and receipts. No manual PySecDec integration is run, no full gauge-invariant two-loop gg acceptance is claimed, and no result demonstrates every possible contraction tree or arbitrary sparse domain. Unsupported performance necessity is a finding, not a silently completed benchmark.
