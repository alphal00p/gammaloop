# Speed up local UV CTs from 4D

## 1. Goal, first implementation action, and scope

**First implementation action:** write the complete contents of this plan, verbatim and excluding the surrounding `<proposed_plan>` tags, to the repository-root file **`SPEED_UP_UV_CTS_FROM_4D.md`**. Do this before changing code or tests. Preserve the plan text and append dated implementation notes, measurements and remaining work below it.

The goal is already active: implement generic local UV CT construction from 4D whose generation time and numerical runtime match or beat the 3D approach. **A maximum 15% slowdown is acceptable for each benchmark; faster 4D remains the target.** Planning or a substantial improvement over the current slow 4D implementation does not complete the goal.

The implementation includes:

- Canonical UV algebra, early term merging, certified source reconstruction, fast CFF energy dispatch and dynamically constructed mapping caches.
- A representation-neutral boundary suitable for future LTD, without implementing LTD now.
- Narrowly targeted fixes in shared symbolic or evaluator code when profiling demonstrates they are necessary, including resolution of GL00’s scalar/tensor preprocessing failure if it persists.
- Correctness and performance validation against both orientation-localized and orientation-erased 3D.

Keep existing user-facing route settings and physical conventions. Introduce no graph-specific optimization, prepared cache tables, benchmark-name branches or compatibility implementation alongside the replacement.

Before rebuilding, preserve the current baseline binaries, source revision, cards and receipts. Pending baseline runs must continue using immutable binaries so subsequent measurements cannot accidentally mix revisions.

## 2. Canonical UV algebra and source reconstruction

### Common representation and ownership

Keep raw `Local4dCts`, recursive `FourDSector` data and outer-Taylor recursion authoritative. Construct a separate canonical projection view inside `Projected4dApproximation::project_local_4d`, before physical-term extraction and per-sector source preparation.

Extend the existing rational-shell extraction machinery to produce:

- A factorized numerator.
- Canonical denominator classes with typed identifiers, exact signed routing, mass, prescription and total power.
- Component bindings containing the actual integration domain, complete source scope, LMB, frozen factors and boundary attachments.
- A separate source-occurrence witness using the existing physical denominator descriptors.

Canonical class identifiers must never masquerade as physical `EdgeIndex` values. Any symbolic class reference must be distinguishable from physical provenance.

Resolve class membership per rational term and component. Normalize completed hard provenance roles 0/1 only. Preserve physical-source and soft provenance, non-UV cographs, cuts, thresholds, frozen localizers and boundary attachments.

A reconstructible hard momentum without a surviving denominator class remains a fixed affine-frame carrier handled by the existing coordinate mapper. Do not invent a pole family or reject an input that the current reconstruction supports. Different masses and pole families retain separate allocation pools.

### Early merging

Replace all-or-nothing denominator-vector comparison with canonical multiset bucketing across denominator-bearing sums and products. Accumulate powers, sum numerators within equal buckets and prune exact zeros.

Keep denominator-free numerator subtrees opaque. In particular,

```text
(A/D² + B/D³)*X + C/D²
```

must produce two buckets,

```text
(A*X+C)/D² + B*X/D³.
```

Do not distribute graph numerators, construct global common denominators or invoke generic factorization to obtain this grouping.

Keep positive typed denominator factors as opaque numerator blocks with their signed binding and full polynomial. Preserve the existing certified pinch/contact handling; normalization must not silently remove source incidence by cancelling these wrappers prematurely.

Group compatible sectors before reconstruction. After each independent component projection, group states with equal remaining requests and boundary bindings: sum numerators when carriers agree, or sum carriers when numerators agree.

Coefficient addition requires equal actual integration semantics. Canonical source reuse across explicitly bound isomorphic frames does not authorize adding independent contours.

### Existing source builder and exact certificates

Use the existing source builder for every canonical class. Retain original incidence unless a serial-path or pure-cycle contraction is certified. Consequently, nonadjacent equal channels—including GL262’s gluon carriers—can merge algebraically without requiring a single powered quotient edge.

When equivalent buckets have different source witnesses, choose one deterministically and certify the combined numerator against it. Never concatenate witnesses and increase denominator multiplicity.

Change completed-UV ownership from original physical edge to certified denominator class. Every eligible hard numerator factor may use an occurrence of its own class. Preserve physical-owner restrictions elsewhere.

Before CFF generation, require:

- Exact numerator equality in one neutral formal momentum frame.
- Independent equality of denominator routing, mass, prescription, multiplicity and domain.
- Valid source incidence, loop rank, measure, boundary attachments and contour.
- Exact diagonal collapse of the occurrence lift, preserving `H=hR`, `P=rR`, hence `H⁰=hrP⁰`.

Update the repository’s reconstruction-policy documentation and associated comments to describe this new invariant. Preserve existing meaningful comments and correctness checks.

The common representation contains no CFF capacities or sampling conventions. Future LTD will consume its powered denominators and numerator directly, deriving residue/derivative requests from denominator orders.

## 3. Fast CFF dispatch and generic live caches

### Deterministic allocation

Extend `EnergyPowerAnalyzer`, `EnergyPowerAssignmentPlan` and the existing planned-expression representation. Replace the affected UV hard-assignment Cartesian frontier search.

Compute structural degrees without numerator expansion:

- Sum: componentwise maximum.
- Product and registered multilinear arguments: addition.
- Positive integer power: repeated factorized base.

For freely assignable unit-energy leaves, use cyclic interval allocation within each certified class. Product children advance offsets; sum branches share their starting offsets; powers retain compressed repetitions. For class degree `D` and occurrence count `p`, the maximum assigned degree is bounded by `ceil(D/p)`.

Keep opaque or indivisible blocks intact. Place them deterministically on eligible occurrences using least-loaded placement, with stable structural ordering, and recompute exact bounds. Never assume unknown functions are multilinear or split positive-denominator blocks through an unproved identity.

Lock the production candidate policy:

- **Ordinary affine input:** one certified candidate and one CFF generation.
- **Nonlinear input:** at most three distinct certified candidates—baseline cyclic/greedy placement, packed unavoidable excess, and reversal of the best scored placement. Use rotation only to replace a duplicate.
- Modify one class at a time, prioritizing the largest nonlinear assigned degree, then positive degree excess, then canonical class identifier. Do not enumerate combinations across classes.
- Construct a challenger only when its factor/block assignments are valid.

Select by actual native residue-map row count before surface conversion or host selection. Break ties by maximum rank, descending rank envelope and stable candidate order. Preserve exact **ordered** capacities in generation keys.

Cache hits and eviction must not change candidate admission or the selected assignment. Count cached candidates against the same logical budget. Introduce no timing-dependent gate or new user setting in this implementation.

### Dynamic preparation and mapping

Extend the existing `PlannedExactSourceNumerator`, which already pairs mapper and assignment behind an immutable `Arc`. Prepare its template through `ExactSourceEnergyMapper`; do not introduce a parallel mapper or symbolic AST.

On first encounter:

1. Intern the actual factorized numerator and exact mapping context.
2. Prepare spatial/temporal decomposition, signed coordinate maps, fixed shifts and replacement patterns.
3. Leave occurrence and loop-energy samples as capture-free parameters.
4. Compute each subtree’s exact parameter dependencies.

For each residue, validate row shape once, convert each required energy once and substitute only relevant parameters. Memoize by immutable template node, complete assignment/context and the restricted sample tuple. Energy-independent subtrees are reused unchanged.

Preserve zero contact samples and the existing order of inactive-energy elimination. The immutable assignment must govern both declared capacities and every sampled substitution.

### Lifetime, keys and memory defaults

Create one nonserialized `Local4dProjectionContext` for the graph’s complete UV computation in both forest orchestrators. Pass it mutably through node projection. Reuse preparation across forest nodes and compatible sectors; clear component-local row memos after each component wave.

Keep mutable caches outside shared immutable numerator objects. Cache owned parsed sources and mapping data, never borrowed graph adapters. Use no global cache, mutable cache field on `Graph`, persistent seed or numerical Monte Carlo sample cache.

Keep these identities separate:

| Cache purpose | Required key |
|---|---|
| Canonical preparation | Canonical algebra and exact source/component binding |
| CFF payload/count | Canonical retained incidence, ordered capacities and semantic options |
| Numerator template | Factorized numerator, immutable assignment and complete signed mapping context |
| Mapped subtree | Template/context identity and exact relevant sample tuple |

Use structural identifiers with exact equality checks; do not stringify large expressions for lookup.

Initial retained-payload budgets are:

| Retention | Budget |
|---|---:|
| Source, analysis and immutable templates | 48 MiB |
| CFF payloads and counts | 64 MiB |
| Component-local row/subtree results | 16 MiB |

Use deterministic LRU eviction, with ceilings of 4,096 preparation contexts, 4,096 CFF keys and 16,384 row-memo entries. Retain at most 64 offset contexts per planned node; stream or recompute additional contexts under the same assignment.

Charge keys, owned containers and Atom payloads at insertion. Oversized entries bypass retention. These limits bound accounted cache payload, not RSS or mathematically required output. Eviction changes computation cost only.

## 4. Implementation sequence and correctness validation

Implement in this order, validating each boundary before proceeding:

1. Canonical projection representation, signed term/component bindings and denominator buckets.
2. Sector/state grouping and witness-backed source adaptation.
3. Class allocation and exact immutable lift certificates.
4. Dynamically prepared mapping templates and graph-wide cache lifetime.
5. End-to-end profiling, necessary shared-code fixes and performance refinement.

Keep changes on existing semantic owners and migrate affected callers together. Remove obsolete production paths once replaced. Do not retain the old per-factor mapper as a cache-miss fallback: a miss constructs the same new generic template.

Add focused tests for:

- Mixed denominator buckets, permuted products, exact zeros and powers of sums.
- Opposite momentum signs with odd numerators and unequal-mass exclusions.
- Positive denominator factors, cancellation/contact terms and absent surviving pole classes.
- Serial paths, pure cycles and nonadjacent equal channels with retained incidence.
- Independent and nested component frames, affine shifts, soft carriers and frozen localizers.
- Opaque blocks, compressed repetitions and valid degree bounds.
- Cold, warm, disabled-retention and forced-eviction execution producing identical selected assignments, bounds and exact coefficients.

Retain the small exhaustive allocation oracles as tests, including the 62-row minimum for the fixed aa_aa full-vacuum source and the nonlinear cases where rank alone chooses poorly. Do not encode their answers in production.

Run the existing signed-occurrence, raised-denominator, nested-rescaling, three-route amplitude and cut/threshold acceptance tests. Run the existing 100+ three-loop scalar graph matrix with its numerator variants, raised propagators and threshold/integrated configurations. Exercise both forest orchestrators.

Add physical three-route coverage for:

- **GL00 and GL01:** two-loop aa_aa amplitudes with DOD1 self-energy subgraphs.
- **GL262:** three-loop aa_aa amplitude containing the DOD2 two-gluon self-energy bubble. Verify its production forest export against the independently identified nested DOD2, DOD1 and DOD0 regions.

Preserve existing test tolerances and expected physics. Ask before changing a failing test’s expectations. Resolve GL00’s scalar/tensor error through the relevant generic implementation if it remains.

Format and run `cargo check` before compilation; finish with appropriate nextest coverage and clippy. Update architecture documentation to describe the implemented boundary. Caches are not serialized; regenerate evaluators as needed without adding compatibility shims.

## 5. Benchmark protocol, acceptance and completion

Use orientation-erased 3D as the primary representation comparison and retain localized 3D measurements. Benchmark the final 3D implementation again after any shared-code change.

For GL00, GL01 and GL262, use matched physical numerators, model parameters, external kinematics, loop bases and evaluator settings. Primary measurements retain the current eager, uncompiled setup, with integrated and threshold CTs disabled to isolate local subtraction. The separate correctness matrix keeps their existing enabled configurations.

Use identical optimized build settings for paired comparisons. Run normal debug-assertion correctness checks as well; disabling assertions is not the optimization.

For each final route and graph:

- Run at least three fresh-process generations, starting with empty generation caches.
- Run generation and numerical benchmarks serially with the same worker limits.
- Benchmark saved states after a priming pass, using at least three timed passes of twenty batches each. Ensure adequate actual measured duration rather than relying only on requested duration.
- Retain the existing benchmark point and its 100× scaled counterpart; use nine spatial coordinates for GL262. Add deterministic nonsingular points for correctness.
- Record actual sample counts, timing uncertainty, settings, revision, binary identity, peak memory and evaluator operation counts.
- Keep profiling runs distinct from final timing runs.

Measure normalization, source certification, degree analysis/allocation, candidate selection, winning CFF generation, template preparation, row mapping, composition, tensor preprocessing and evaluator optimization separately. Record bucket counts, unique sources, cache hits/evictions, native rows, template/sample-tuple counts and generated evaluator size.

**All losing CFF generation and discarded preparation count as dispatch-selection overhead.** Allocation plus selection/cache work must remain below 10% of cold local-UV generation. Winning generation, template construction and row mapping must also be reported independently.

For each required physical graph, the final median 4D/erased-3D ratio must be **at most 1.15** for:

1. Generation time.
2. Evaluator time per sample.
3. Total time per sample.

Report individual results; averages cannot conceal a failing graph. Repeat borderline measurements when uncertainty prevents a reliable conclusion. The preferred result is a ratio below 1.00.

Record performance across the broader scalar matrix and investigate regressions through generic structural causes. Do not tune graph identities, relax pointwise checks, seed caches or omit difficult cases to pass.

Pending or failed original-4D runs remain explicitly labelled. They do not establish successful timing ratios, and all optimized required cases must complete.

If a gate fails, use the measurements to continue the generic implementation, including the authorized narrowly targeted shared-code fixes. Do not declare success or silently change the dispatch policy, tolerance or benchmark population. Append results and unresolved findings to the plan file.

**Mark the active goal complete only when the implementation, correctness checks, required generation/runtime gates and final documentation are complete.**


## Implementation notes — 2026-09-12

- Saved the accepted plan verbatim before changing code or tests. Its original
  contents have SHA256 `5c3ce40e19a807015e346df66904d6e770a7bb7ead1128c3e15437578851428c`.
- Created and pushed feature branch `faster_uv_from_4d`; the plan milestone is
  commit `d533c8ca4`. Subsequent commits use the user-requested identity
  `ValentinHirschi <valentin.hirschi@gmail.com>`.
- Preserved baseline revision `78395e3ab3ddd8d8f62b2f674d7488484eace197` and both
  optimized binaries under `tests/artifacts/aa_aa_uv_slowdown/baseline_78395e3/`.
  The manifest records exact SHA256 binary identities. Subsequent original-route
  runs use these immutable copies; the already-running GL262 process retains
  its original executable. Existing cards, logs, and receipts remain preserved
  in the same artifact tree.
- Implementation is in progress. No optimized timing or correctness gate has
  passed yet; the active goal remains incomplete.

### 2026-09-12 — canonical projection and live mapping implementation

- Added a separate canonical UV projection view and typed denominator-class
  symbols. Raw sectors and outer Taylor recursion retain physical provenance.
  Rational-shell extraction groups equal multisets without expanding graph
  numerators. Positive denominator blocks retain their complete polynomial.
- Retained incidence witnesses and added exact signed class routing/mass checks.
  Completed hard factors use certified class occurrences. Pinched hard carriers
  are lowered to an exactly certified fixed affine lift through retained base
  occurrences; they do not create pole families.
- Replaced the hard Cartesian assignment frontier with cyclic/greedy allocation,
  compressed powers and at most three native-row-scored candidates. Ordered
  capacities remain part of the CFF cache key.
- Added generic live immutable numerator templates and bounded deterministic LRU
  retention, with per-subtree parameter dependencies and component-local row
  reuse. Both forest orchestrators share one projection context per graph.
- Added cross-term state grouping after each independent component, allowing
  equal remaining requests to reuse their summed numerators or carriers.
- `cargo check -p gammalooprs --tests --profile dev-optim -j 8` passed after the
  API migration. The licensed focused nextest snapshot is still in progress.
  Two new factorized-equality assertion corrections remain pending user approval;
  no existing physics expectation or tolerance has been changed.
- An existing owner-relabelled projection test exposed inconsistent output
  dialects (physical OSE aliases versus literal on-shell square roots) for equal
  vacuum channels with different witnesses. Completed Taylor-vacuum sources now
  use one canonical spatial frame. This fix still requires the focused and
  physical acceptance runs.
- Exact pre-generation lift certification, detailed timing records, full
  acceptance coverage and all final performance gates remain in progress. No
  optimized generation/runtime ratio is established by this milestone.

### 2026-09-13 — cache integration and first physical diagnostic

- Pushed implementation milestone `b8396f3a3080d196a947879eb877a2d842d871ef`
  with the requested ValentinHirschi author/committer identity. Preserved its
  debug-assertion-enabled optimized CLI as an immutable diagnostic snapshot:
  SHA256 `73cfeefe68eec92571a84c55cc95d81e0efb01d39d1c1bb0ec52b4dd9a184424`.
- The milestone's focused nextest run passed 207 of 213 tests. Remaining failures
  exposed numeric-sign normalization at exact certificates, a square-root
  routing dialect difference, and new-test representation assumptions. The user
  subsequently authorized routine test corrections without further approval;
  physical expectations and tolerances remain unchanged.
- Numeric-only distribution (`expand_num`) resolves the independently reduced
  nested-sign certificate reproducer while preserving products and powers of
  graph numerators. The production certificates now apply it only after the
  structural equality fast path. These changes await the next test build.
- Canonical sectors, owned source analysis, and immutable numerator templates
  now share the 48 MiB/4,096-context preparation LRU. Source hits bypass borrowed
  adapter construction, and tests exercise cold, warm, disabled and evicted
  execution. Integration and validation are underway.
- The immutable milestone GL00 direct-4D diagnostic completed local projection,
  then failed in shared tensor preprocessing with the original
  `lazy tensor sum mixed with scalar sum terms` error. Outer elapsed time was
  446.419 s; sampled VmHWM peaked at 13,262,540,800 bytes (12.35 GiB). This was a
  profiling run concurrent with the pending original GL262 baseline, and is
  neither a successful generation nor an acceptance timing.
- Its five-second profile attributed 88.97% of inclusive samples to tensor
  network graph flattening/deletion, predominantly conversions between Linnet
  parent-pointer and child-vector stores. A separate generic rank-zero lazy
  tensor addition fix is being implemented because the same valid sum is
  accepted or rejected depending on operand order.
- All final physical generation/runtime ratios, the full scalar acceptance
  matrix, dispatch-overhead gate, clippy and final documentation remain open.

- The next coordinated run passed 296/299 tests across gammalooprs, Spenso and
  Linnet. All cache-retention, native-row selection, signed-allocation and tree
  swap/deletion checks passed. The three remaining failures are being reduced:
  a nested certificate with opposite additive bases under an integer square,
  a factorized contour equality assertion, and a new higher-rank rejection
  test that initially failed at graph construction before reaching the adder.
- Replaced the profiled Linnet whole-store swap conversions with local pointer
  relabeling. Exhaustive tests cover 720 six-node forests, fresh and serial swaps,
  and all 512 deletion subsets of a nine-node forest against ChildVecStore.
  A separate three-pass generic diagnostic measured 100k-node leaf/bounded-degree
  swaps at 261/567 ns versus 2.54/3.40 ms previously (about 9,710x/5,993x).
  A high-degree root control improved about 1.61x. Exact final storage matched
  in every case. These are isolated-operation results, not generation ratios;
  receipts and raw passes are under `linnet_swap_benchmark/results_b839_dirty_r2`.
- The nested certificate now normalizes the sign of additive bases under integer
  powers using a bottom-up rewrite. Fractional powers keep their original
  branch-sensitive form. The reduced reproducer is exact and idempotent without
  distributing numerator products or powers; full regression validation remains.

### 2026-09-13 — successful GL00 diagnostic after shared-library fixes

- The immutable cache/Spenso/Linnet diagnostic CLI (SHA256
  `895e2ea230ebc6fbc0e04f468c226a87709b11a63ed09e264165e6819c528517`)
  successfully generated the physical GL00 direct-4D evaluator. Reported
  generation time was 547.850254805 s: expression construction 113.456882143 s,
  Spenso preprocessing 226.860166267 s, Symbolica evaluator construction
  207.533206395 s, and compilation 0 s. One evaluator was saved. The generation
  report recorded 16,986,910,720 peak RAM bytes; the outer process receipt
  recorded 550.035878733 s and a sampled RSS peak of 17,278,410,752 bytes.
- This was a diagnostic run concurrent with unit compilation and the original
  GL262 localized-3D baseline, with three five-second profiling windows. It is
  not a final acceptance timing or a final optimized/baseline ratio. State,
  reports, receipts and profiles are preserved in
  `tests/artifacts/aa_aa_uv_slowdown/diagnostic_cache_shared/GL00_4d_direct_r1`.
- The sampled stages changed from color simplification (84.74% inclusive,
  including representation matching and coefficient zero tests), to network
  ready-batch planning (94.80%, predominantly graph-sized bitsets), then
  Symbolica Horner construction (73.55%). These nested percentages describe
  separate brief windows, not additive whole-run stage costs. The original
  scalar/lazy-tensor failure was passed successfully. Saved-state runtime,
  pointwise comparison and operation-count diagnostics are being collected;
  the final repeated physical benchmark gates remain open.

- Saved-state GL00 diagnostics completed with the same immutable CLI for both
  the new direct-4D evaluator and the original erased-3D evaluator. Both state
  formats loaded successfully, including the existing operation counter. Each
  runtime had one priming pass and three timed passes of 20 batches, serialized
  under the shared benchmark lock. Median full-integrand/evaluator times were
  1,276.553880/1,228.862211 microseconds for direct 4D (8,873 samples), versus
  96.786247/79.184658 microseconds for the saved erased-3D baseline (29,263
  samples). This leaves diagnostic gaps of 13.19x/15.52x; these are not final
  acceptance ratios and do not meet the intended runtime performance gate.
- Direct-4D evaluator counts were 246,394 instructions, 144,004 multiplications,
  166,502 additions, 36 inversions and 19 function calls. The erased-3D program
  had 26,557 instructions, 19,092 multiplications, 14,057 additions, 35 inversions
  and 18 function calls. Their full complex integrands agreed at the base and
  100x-scaled spatial point to relative differences 7.4643e-15 and 3.1068e-15.
  Integrated and threshold CTs remained disabled. Exact raw values, timing
  receipts, operation counts and comparison metadata are in
  `diagnostic_cache_shared/saved_runtime_summary.json` and
  `saved_runtime_comparison.json`; the original baseline state was read-only.

### 2026-09-13 — exact affine certificates and profiling refinements

- The next focused snapshot passed 297/300 tests; all retained-cache, CFF
  allocation, Spenso network and Linnet tree checks passed. The remaining
  failures isolated two exact source-frame boundaries and a new test's raw
  contour convention. The expected physical numerator and contour were kept.
- Positive denominator wrappers whose completed hard carrier has no surviving
  pole now retain an exact whole-block affine-rewrite witness. The existing
  planned factor admits joint fixed dependencies only for that certified
  wrapper; physical/soft restrictions and rejection of unwitnessed mixed blocks
  remain. Coverage passes the actual canonical projection for completed hard
  roles 0/1/3 and checks that soft role 2 remains unchanged.
- The occurrence-diagonal certificate now compares fixed external four-vectors
  in the same formal spatial/temporal frame as the reconstructed affine shift.
  Only the projector action on a proven affine energy sum is linearized;
  numerator products, powers and opaque blocks stay factorized. The pinched
  contour oracle now uses the existing componentwise VariantLocal-to-contour
  convention instead of applying the generalized core sign a second time.
- Profile-guided shared changes replace representation wildcard construction
  with the existing direct Atom visitor and graph-wide scheduler bitset scans
  with the existing crown iterator. The representation scan has an exact legacy
  matcher oracle, including zero-argument functions, nested wrappers and power
  exponents. Scheduler traversal, admission, payloads and debug assertions are
  preserved; its dense-filter oracle is being added.
- Formatting and cargo check passed for core, Spenso, Linnet and integration
  tests with debug assertions enabled. The expanded focused nextest build is
  running. Performance gates, broad scalar/physical acceptance and clippy are
  still incomplete; the GL00 runtime diagnostic above explicitly fails the
  intended performance threshold.

- The expanded focused suites progressed through 321/323 and 322/324 passing
  checks. The original nested physical reconstruction now passes, along with
  positive affine blocks and the dense scheduler oracle. New regression checks
  exposed two further normalization assumptions: signed numerical extraction
  can need another bottom-up pass, and tagged linear function construction can
  yield a sum instead of a function. The former now stops at structural equality;
  the latter fixture compares every resulting function against the unchanged
  legacy predicate. These final corrections await the next focused run.
- The representation predicate's exact contract is whole-argument matching:
  the old `ReplaceBuilderExt::matches` uses `partial(false)`. Production now
  checks immediate argument heads/arity directly and retains only the existing
  explicit chain/trace/broadcast recursion. The legacy oracle rejected the
  initial overly recursive visitor; no expected behavior was relaxed.
- Component state grouping now repeats after summing numerators, because that
  sum can create a new identical numerator key. Each useful pass strictly
  decreases the state count; no rational or numerator expansion is introduced.
- The saved GL00 graph expressions localize the remaining size discrepancy:
  the common original term has matching structure, while two 4D UV terms are
  gamma-free and contain 26,290 metric calls; their erased-3D counterparts retain
  1,746 gamma calls and no metric calls. Both routes share the final mapping and
  dot-expansion operations. Local 4D's earlier Taylor stage unconditionally
  simplifies traces/Schoonschip, duplicating the integrated path's existing
  analytic preparation. Deferring that duplicate local spin algebra is the
  next generic experiment; success and unchanged row counts are not assumed.


### 2026-09-13 — GL01 saved-state diagnostic completed

- The same immutable cache/Spenso/Linnet CLI (SHA256
  `895e2ea230ebc6fbc0e04f468c226a87709b11a63ed09e264165e6819c528517`)
  generated GL01 successfully in 682.069602369 s: expression construction
  130.396855471 s, Spenso 317.377554948 s, Symbolica 234.295191950 s,
  compilation 0 s. Reported peak RAM was 23,740,649,472 bytes. The outer
  receipt records 684.484785132 s and 24,388,227,072 sampled peak RSS bytes.
  This snapshot predates the latest representation predicate, network scheduler
  and proposed deferred local Dirac algebra. Generation ran concurrently with
  unit compilation and the old GL262 localized baseline, so it remains a
  diagnostic rather than a final acceptance measurement.
- Saved-state replay used that immutable CLI for both direct4D and the original
  erased3D evaluator, with the state read-only. One priming pass preceded three
  timed passes of 20 batches, and all runtime, point inspection and operation
  counting held the shared benchmark lock. Median full-integrand/evaluator
  times were 1,418.340558/1,366.547650 microseconds for direct4D (8,170 samples)
  and 96.251632/78.577745 microseconds for erased3D (29,743 samples). The
  diagnostic gaps are 14.74x/17.39x; the runtime performance gate remains open.
- Direct4D has 274,845 instructions, 159,144 multiplications, 198,851 additions,
  36 inversions and 19 function calls. Erased3D has 26,458 instructions, 19,109
  multiplications, 13,984 additions, 35 inversions and 18 function calls. The
  complex full-integrand values agree to relative differences 2.1740e-15 at
  the base point and 4.4902e-16 after scaling all spatial coordinates by 100.
  Integrated UV and threshold counterterms were disabled in both routes.
- Exact timing, counts, point values and executable/state provenance are under
  `tests/artifacts/aa_aa_uv_slowdown/diagnostic_cache_shared`, in
  `GL01_saved_runtime_summary.json`, `GL01_saved_runtime_comparison.json` and
  the two GL01 case directories' receipts. No further generation was launched.
  The separate GL00 source comparison identifies early gamma-trace contraction
  as the first representation divergence: its two no-OSE UV terms account for
  95.84% of the direct source and are 23.22x the combined erased terms' text.
  `spf` is the SpinFundamental representation, not a scalar-product shorthand.
  That causal owner trace motivates the deferred-algebra experiment; the new
  GL01 timings alone do not prove the same term-level cause.

- Milestone `f3ec6fcb1bb367ce71c8f24f1901fd1f81c28b9d` was committed and pushed
  as ValentinHirschi <valentin.hirschi@gmail.com>. Its final focused nextest
  run passed **324/324** checks in 6.731 s after a 2m50s optimized build with
  debug assertions enabled. Both newly corrected regression checks pass,
  including signed-normalization idempotence and exact legacy scan parity.
  The next source change defers local spin algebra to its existing analytic or
  numerical owners; this is not yet measured or accepted as a performance fix.

### 2026-09-13 — deferred local spin algebra experiment

- Local Taylor construction now keeps its tensor products/chains/traces after
  metric simplification. Rescaling, Taylor order and t=1 substitution are
  unchanged. The integrated counterterm already invokes its own complete
  analytic simplifier on a copy before Vakint; that owner retains the moved
  denominator-topology comment. No dispatch policy or evaluator setting changes.
- Added focused real gamma/chain/trace allocation coverage and a local Taylor
  trace-retention test with identical raw output for the integrated flag on/off.
  Formatting and cargo check passed; the broader core UV nextest run is building.
  Physical timing, pointwise and native-row effects remain to be measured.
- `local_4d_construction` now reports total construction, nested Taylor time,
  overhead and sector counts in the shared UV-limit owner. Summing its elapsed
  times plus projection times gives a conservative local-UV generation subtotal;
  pre-call recursion-input preparation and disconnected-union bookkeeping are
  excluded, and nested Taylor milliseconds must not be added again. The profile
  summarizer reports dispatch/selection/cache work against this explicit subtotal.

- The explicit inert CYCLIC/SYM/ANTISYM heads now use the existing multilinear
  planned expression, without changing their global symbolic attributes. The
  broader core UV/energy/source/CFF suite passes **286/286** (15.883 s; build
  2m43s). The analytic ordering, gamma5 guard and nested Taylor checks pass.
- The next immutable diagnostic CLI is based on `f3ec6fcb1` plus recorded patch
  `ecc4afcb6fb1c35491823d127d984ef2967596a9105a1e544b6bba551da1459f`;
  binary SHA-256
  `16248ab4974403799c1600c2207755f643a94cf8287ccffc2b3abe662c3fd429`.
  `cargo check` passed before the 3m29s assertion-enabled dev-optim CLI build.
  Cards, binary and source receipts are in `diagnostic_deferred_spin`.
- Its GL00 projection diagnostic takes 165.948 ms, with 11.602 ms of separately
  measured construction. Native rows remain 196 in total, maximum 62, with ten
  logically admitted proposals across eight requests. Immutable template nodes
  fall from 9,071 to 715; template construction takes 62.197 ms. Dispatch,
  selection and all cache work total 11.994 ms, or 6.76% of this conservative
  construction-plus-projection subtotal. This is a profiling diagnostic, not a
  final unprofiled gate. Full generation remains pending and still shows costly
  downstream work with substantial memory growth, so no end-to-end improvement
  or runtime acceptance is inferred from the faster projection.

### 2026-09-13 — shared preprocessing follow-up

- Deep sampling of the deferred-spin GL00 diagnostic attributes 94.29% of 200
  samples (no lost samples) to `collect_orientation_if`. The actual source has
  neither theta, IF nor residue selectors. Added exact no-op guards in that
  existing owner and selector parametrization, leaving guarded inverse handling
  intact. New factorization regressions accompany the guards. Tensor execution
  is separately measured at 169.099 s; a new byte-count field on its existing
  completion event will expose the scalar output size without an expression dump.
- Reused `Numerator::color_simplify` on reduced numerators before Taylor in both
  local routes. Six real GL262 scopes, including the DOD2 bubble and nested
  prefixes, have exact raw partition and early/late cograph recombination
  equality with unchanged open interfaces. Bubble color becomes three times an
  open adjoint metric in about 2.45 ms; full-graph color simplifies in 9.58 ms.
  This artifact probe does not establish end-to-end timing. A focused generic
  open/nested color test also keeps each momentum factor on its owning edge.
- Full-workspace formatting and check34 pass. Clippy completes successfully;
  redundant borrows and a test-only clone of a Copy value have been corrected.
  Three tuple-signature complexity warnings remain; no compatibility helper
  layer was added merely to silence them. The next CLI build and targeted test
  rerun include these shared changes. The 167-case scalar acceptance run on the
  preceding immutable test executable remains in progress; its hash and source
  binding are retained in `diagnostic_deferred_spin/scalar_matrix_1_manifest.json`.
- The deferred-spin physical diagnostic remains pending in the avoidable
  downstream collector stage and has exceeded the prior 550 s successful run.
  It is not a completed timing result. New immutable diagnostics will use the
  corrected shared preprocessing; all final per-graph performance gates remain
  open.

- The follow-up focused nextest run passes **407/407** tests in 15.833 s
  (2m59s build), covering UV/CFF/source/energy logic, evaluator selectors, symbol
  conventions, Spenso network/shadowing and Linnet child-pointer behavior.
- The superseded deferred-spin GL00 diagnostic was explicitly terminated after
  1,388.644 s, exit -15, with peak RSS 43,915,075,584 bytes. Its termination note
  records the measured no-op collector cause. It has no completed generation
  or runtime result; the original immutable GL262 baseline was not interrupted.
- The new immutable CLI hash is
  `fe064ff6d4c191b20f4d56de9077c2280e293373792e57502ca929efa442c5b7`,
  based on f3ec6fcb1 plus recorded patch
  `53fd4a4e1030ddde3496302e7178835df2fbbebbbdc3cb6a128caacd80fb610f`.
  The new GL00 tensor execution takes 86.515 s in this concurrent diagnostic,
  but produces **3,450,351,988 scalar Atom bytes** from 1,416,191 input Atom bytes.
  A 333-sample tensor-phase profile (no lost samples) attributes 71.26% to
  scalar-product folding/normalization and 12.76% to balanced scalar sums.
  Existing component optimization is being tested on the identical input to
  diagnose this growth. No graph numerator expansion or new production setting
  has been introduced for that experiment. Physical generation and runtime
  acceptance remain pending.
