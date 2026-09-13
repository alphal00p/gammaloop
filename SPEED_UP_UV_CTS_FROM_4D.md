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

### 2026-09-13: locating the remaining tensor growth

- Milestone `a1673901fa3417ce3b4eebc480d6edd658e87255` was committed and
  pushed as ValentinHirschi. The current immutable diagnostic CLI above contains
  its production code. These are concurrent diagnostic measurements, not final
  acceptance timings.
- Matched current-code GL00 generations completed: direct 4D **829.149287255 s**
  (expression 5.174840171 s, Spenso 90.901857385 s, Symbolica 733.072589699 s),
  erased 3D **19.740054989 s** (expression 6.642393406 s, Spenso 4.111505172 s,
  Symbolica 8.986156411 s). Both use the same immutable binary, physical inputs
  and settings except the route. Direct reported peak memory is 41,532,669,952
  bytes. Generation remains far outside the required gate.
- The shared tensor owner receives similarly sized inputs: direct 1,416,191
  bytes and erased 1,432,379 bytes. Erased produces 41,119,289 scalar bytes in
  2.776 s tensor execution; direct produces 3,450,351,988 bytes in 86.515 s,
  about 83.91 times the output size. Standalone replays reproduce both output
  sizes exactly. The original, unsubtracted term is identical at this boundary.
- Existing `HornerAtomComponents<4096>` reduces the direct scalar size by only
  160,846 bytes (0.00466%); it is not a useful remedy. A scalar-fold trace shows
  every multiplying fold has exactly two factors. The sampled late hotspot is
  copying already enormous scalar results, not multiplying many factors.
- Existing lazy tensor sums are active; all 21 logged nontrivial sums remain
  lazy and no distributed-term threshold forces materialization. The first
  large external-vector contraction receives 20 already completed rank-four
  tensors, 5,120 component entries and 246,405,632 bytes. The growth therefore
  precedes the external contractions. The full-vacuum contribution has 62
  mapped numerator copies in both routes, but direct groups them into 12
  carriers before attaching the external cograph. The next diagnostic tests
  closing the native residue summands before forming their open tensor sum,
  while preserving every incoming graph numerator's factorization. No
  production grouping policy has been changed on this evidence alone.

- The matched GL00 diagnostic completed with the same immutable `fe064ff6`
  executable for direct4D and a fresh erased3D generation. Generation took
  829.149 s versus 19.740 s (42.00x): expression work 5.175/6.642 s,
  Spenso 90.902/4.112 s, and Symbolica 733.073/8.986 s. Reported peak RAM
  was 41,532,669,952/688,590,848 bytes. These single generations overlapped
  other diagnostic work and remain observational, not final acceptance timings.
- The subsequent saved-state runtime commands and counters held the shared
  benchmark lock. Three primed passes of twenty batches each retained 7,200
  direct samples and 109,603 erased samples. Actual measured Total durations
  were 4.061/4.038/4.090 s for direct and 3.895/3.686/3.728 s for erased.
  The erased route's initial approximately one-second attempts were preserved
  but excluded; adaptive retries requested nineteen seconds and each exceeded
  the three-second actual-duration gate. Median Total time was
  **1691.152/103.050 microseconds (16.411x)** and evaluator time was
  **1633.378/91.313 microseconds (17.888x)**, direct/erased respectively.
  The original long-running GL262 baseline remained active, so these are
  diagnostic measurements even though runtime commands were serialized.
- Direct/erased instruction counts are **408,227/31,382 (13.008x)**, with
  307,911/23,192 multiplications and 250,738/17,309 additions. Both programs
  have 35 inversions, eighteen function calls, 214 inputs, one output and a
  62-orientation catalog. Full-amplitude plus local-UV values agree to relative
  complex norm differences **1.79e-15** at the base point and **2.94e-15**
  at the point with spatial coordinates scaled by 100. Integrated UV and
  threshold CTs are off; evaluation is eager and uncompiled in both routes.
- The exact current before-network inputs are 1,416,191 direct Atom bytes and
  1,432,379 erased Atom bytes. The erased route creates sixteen scalar aliases
  totalling 69,918 bytes (largest 5,941 bytes); direct creates none at the same
  4,096-byte threshold. Tensor execution yields 3,450,351,988/41,119,289 scalar
  Atom bytes in 86.515/2.776 s. This locates the large growth after the actual
  shared tensor boundary and motivates comparing the differing factorization
  of its inputs. It does not establish an algebraic error or a successful
  performance fix. Receipts, all timing attempts, counts, point values and
  exact binary/state provenance are in
  `tests/artifacts/aa_aa_uv_slowdown/diagnostic_early_color_and_guards/GL00_adaptive_runtime_summary.json`
  and its case directories. The saved erased reference is the fresh state in
  this same diagnostic directory, despite the replay output directory's
  `GL00_3d_erased_baseline_state` suffix.

- A bounded finite-tensor closure experiment now confirms that mechanism.
  It identifies the 21 native residue aggregation scopes in this captured GL00
  input and preserves all 298 complete numerator summands, including their
  internal sums and powers. Exact reverse reconstruction of every original
  scope passes. Attaching the existing external tensors before aggregating
  those summands changes input size from 1,416,191 to 1,547,061 bytes (+9.24%),
  costs 17.2 ms, and reduces tensor execution from 84.77 to **3.755 s**.
  Result extraction falls from 1.812 to 0.0189 s and scalar output from
  3,450,351,988 to **38,471,036 bytes**, slightly below erased 3D's
  41,119,289 bytes. The diagnostic's inspected gamma/polarization predicates
  are confined to ignored artifacts and must not enter production.
- The implementation owner for this refinement is the existing tensor network
  graph/execution boundary. Its pending Product/Sum structure retains the
  necessary contraction context after all UV mappings are complete. Closing
  tensor indices through that pending sum can preserve canonical UV grouping,
  selected immutable assignments, cache identities and the original input
  Atom. Selection must use generic tensor slots and reduce the open boundary;
  scalar factors, products of sums, powers and unknown functions must not be
  distributed. Native-row sidecars in projected UV coefficients would instead
  have to survive every prior component grouping and both hard/soft mappings,
  so they are not being introduced as a second representation.
- The existing `finite_part_quark_lo` analytic regression passes **1/1** in
  0.396 s, with dev-optim assertions enabled and the standard licensed runtime.
  Its separate target passed cargo check (23.30 s) before compilation (3m10s);
  the running scalar matrix's `test_runs` executable was not rebuilt. Binary
  provenance, nextest metadata, check/build logs and the execution receipt are
  in `tests/artifacts/aa_aa_uv_slowdown/finite_part_quark_lo_check_1`.
- Bounded source-only diagnostics reused the same immutable `fe064ff6` CLI and
  matched physical cards under the shared lock. GL01 reached the exact
  before-network event after 5.172 s, containing four terms and 1,420,147 Atom
  bytes; its payload was saved and the newly launched process deliberately
  terminated. GL262 reached the 180.028 s cap without that payload and was
  deliberately terminated too. Neither is a completed generation timing or a
  numerical failure. Receipts and the GL01 input are in
  `diagnostic_early_color_and_guards/GL01_4d_source_capture_r1` and
  `GL262_4d_source_capture_r1` under the measurement artifact directory.
- During the capped GL262 capture, a six-second 19 Hz sample of its worker and
  auxiliary resource-monitor thread collected 122 user-CPU samples with none
  lost. Of those samples, 90.98% include `simplify_final` → `simplify_color` →
  `collect_chains` → `collect_symbol`, and 88.52% include symbolic zero tests.
  The last logged forest node was topo index zero, DOD0, with no parents.
  This isolates an upstream color-collection bottleneck; it does not establish
  a DOD2-specific issue or describe the whole generation. The original GL262
  baseline remains active and was neither paused nor terminated by these runs.
- The scalar acceptance matrix completed **167/167 tests passing** in
  2,949.565 s. It used the immutable deferred-spin test executable
  `f7a5f864bd7ba23756b5e0fd9d9bd7e76b4a5d7a7aa85b8897b32907caacfd17`,
  with its source manifest recorded before the run. This validates the canonical
  projection, assignments and caches across the existing scalar graph,
  numerator, raised-propagator and integrated/threshold configurations on
  that snapshot. The pending finite-tensor preparation and final source still
  require a new matrix run; this result is not a final performance gate.
- The GL262 root profile identifies another missing early-color boundary:
  the untouched cograph numerator in final assembly. Both final assembly routes
  now call its existing color simplifier before any residue mapping repeats
  that numerator. The later color pass remains responsible for color indices
  closed by attached UV terms or projectors. Raw stored graph numerators remain
  unchanged; the new order will be checked in the next paired snapshot.


### 2026-09-13 — generic finite tensor preparation

- Added preparation on the existing Spenso contraction graph: a ready tensor
  whose complete exposed boundary meets one immediate sum is attached to each
  intact arm before that sum materializes. It uses tensor-store references and
  exact self-dual slot bindings, preserves scalar spectators and unrelated
  incidence, and never multiplies sums through other sums, powers or unknown
  functions. The original symbolic Atom and UV projection representation remain
  unchanged. Each move strictly decreases an existing sum's exposed rank.
- The first broad run completed 416/425 tests successfully. Eight failures came
  from accidentally overriding the repository filter for explicitly marked
  `::failing::` tests; these are recorded separately from the new regression.
  The new independent coordinate test exposed a real implementation defect:
  the sum-first operand order closed four boundaries, but the outside-first
  order closed one and left six dangling slots. Its expected values and counts
  were not changed.
- The cause was the existing `SmartEdgeVec::set_flow`: its source/sink setters
  ignored dangling Identity edges. It now uses the existing involution flip,
  with a regression checking dangling and paired flow, payload consistency,
  idempotence and exact round trips. The preparation also exposed stale slot
  bookkeeping: `NetworkGraph::delete` did not truncate its slot-order vector,
  and `merge_ops` bypassed that owner entirely. Both now use the aligned graph
  deletion path. Restoration additionally checks the exact residual slot count.
- The broader scalar fixture now always emits its existing paired setup-plus-
  generation timings, so the next matrix run can retain per-case observational
  comparisons. No physics expectation, numerator, route setting or tolerance
  changed. Final physical generation/runtime gates remain open.

- The corrected focused run passes **633/633 tests** in 58.093 s (3m15s
  compilation), including the unchanged independent coordinate oracle, both
  operand orders, four execution strategies, graph seam/slot regressions,
  UV/CFF/source/energy coverage, evaluator guards and the Idenso library suite.
  The preceding cargo check passes in 11.07 s. The source retains normal debug
  assertions. A newly built immutable CLI will provide the next physical
  measurements; neither these unit results nor the small structural replay
  establish the generation/runtime acceptance gates.

- Clippy completes successfully; the same three documented tuple-complexity
  warnings remain and the new tensor/graph changes introduce no warnings.
  The next immutable assertion-enabled CLI is based on a1673901 plus patch
  `ed44598ebfb02533a79109f9eccdd0118a713bfe56db7afa220a20578a3af175`,
  with binary SHA-256
  `40f7f703f770f6942f8fdf21cf5be9089a9dae0a8836ff80ed880926f12738fd`.
  Its 3m58s build, source patch, binary and six matched diagnostic cards are
  recorded in `tests/artifacts/aa_aa_uv_slowdown/diagnostic_tensor_boundaries`.


### 2026-09-13 — first physical generation/runtime improvement

- Milestone `d21ce5c92e7030bfd5cbf5a977b04ea6cedb04ef` was committed and
  pushed as ValentinHirschi. All 633 focused tests pass; the frozen `40f7f703`
  CLI contains this production implementation. The separately built physical
  and scalar integration executables are also copied and hashed in
  `diagnostic_tensor_boundaries/integration_binaries` before further edits.
- Production finite-tensor replay of the captured GL00 and GL01 inputs makes
  84 boundary moves each. Preparation takes 4.221/4.213 s, execution takes
  4.301/4.314 s, and scalar outputs are 39,489,750/39,252,026 bytes. The GL00
  output is reduced from 3,450,351,988 bytes without modifying its incoming
  symbolic Atom. Exact coordinate tests and graph seam certificates remain
  unchanged. Preparation itself is now substantial; batching compatible
  closing leaves of one sum is being evaluated as a generic refinement.
- Fresh matched GL00 generations complete in **18.840463931 s (4D)** and
  **19.603097009 s (erased 3D)**, ratio **0.9611**. Expression work takes
  4.712439665/6.481769429 s, Spenso takes 8.976910534/4.202542809 s, and
  Symbolica takes 5.151113732/8.918784771 s. Reported peak memory is
  651,010,048/832,045,056 bytes. These are single concurrent-build diagnostics,
  not the required three fresh-process final measurements.
- Matched saved-state runtime uses three adequate primed passes of twenty
  batches each. Median evaluator time is **79.3137/94.3283 microseconds**
  (4D/erased, ratio **0.840827**); median Total time is
  **90.6994/106.1245 microseconds** (ratio **0.854651**). Direct retains
  121,548 samples over actual measured pass durations 3.722/3.741/3.582 s;
  erased retains 104,935 samples over 3.615/3.757/3.839 s. Pointwise relative
  complex differences are **1.1943e-15** at the base point and **1.2644e-15**
  at the 100-times-scaled point. Direct instruction/multiplication/addition
  counts are 26,465/19,294/15,347 versus erased 31,382/23,192/17,309.
  Full state-content identities are verified before and after read-only reuse.
  The same frozen binary and matching physical cards serve both routes.
- Generation and runtime commands are serialized under the shared benchmark
  lock, but the original GL262 baseline and compilation remain concurrent, so
  these are diagnostic observations. GL01 and GL262 remain required, and no
  final acceptance gate is declared complete. Original GL262 is still pending
  after five hours without a saved result.
- Profiling now records per-call preparation/row-cache work in addition to
  cumulative lifetime counters, allowing exact overhead sums across contexts.
  The summarizer refuses to infer multi-context cost from cumulative maxima;
  old single-context profiles require an explicit scope assertion. Synthetic
  delta, uncertified cumulative and certified single-context controls pass.

### 2026-09-13 — first successful tensor-boundary physical comparison

- Frozen CLI `40f7f703f770f6942f8fdf21cf5be9089a9dae0a8836ff80ed880926f12738fd`
  contains the production changes committed in `d21ce5c92`. The exact precommit
  source binding is retained in `diagnostic_tensor_boundaries/source_manifest.json`.
  Standalone production preparation made 84 boundary moves on each captured
  GL00/GL01 input. Preparation took 4.221/4.213 s and tensor execution
  4.301/4.314 s, producing 39,489,750/39,252,026 scalar Atom bytes. These use
  the generic network method, with the original input Atom unchanged.
- Full GL00 direct generation now completes in **18.840464 s**, versus
  **19.603097 s** for a fresh erased3D generation with the same binary. Direct
  expression/Spenso/Symbolica times are 4.712440/8.976911/5.151114 s; erased
  times are 6.481769/4.202543/8.918785 s. Reported peak RAM is
  651,010,048/832,045,056 bytes. Both produced one evaluator.
- Serialized adaptive saved-state measurements give direct/erased median Total
  times **90.699/106.125 microseconds (ratio 0.854651)** and evaluator times
  **79.314/94.328 microseconds (ratio 0.840827)**. Three twenty-batch passes
  retained 121,548/104,935 samples. Actual measured durations are
  3.722/3.741/3.582 s for direct and 3.615/3.757/3.839 s for erased; initial
  short attempts remain preserved and excluded. Direct/erased instructions are
  26,465/31,382, multiplications 19,294/23,192, and additions 15,347/17,309.
- Base and 100x-point values agree at relative complex norm differences
  **1.1943e-15/1.2644e-15**. The full saved-state file inventories were hashed
  and remained unchanged through read-only runtime/counting. Results and all
  attempts are in `diagnostic_tensor_boundaries/GL00_adaptive_runtime_summary.json`.
  These are promising single-generation diagnostics under concurrent build and
  original-baseline load. They do not replace the required three fresh isolated
  generations or the GL01/GL262 correctness and performance gates.

- The matching GL01 diagnostic completes in **18.838888610/19.576989172 s**
  (4D/erased), with reported peak memory **626,102,272/803,491,840 bytes**.
  Expression/Spenso/Symbolica times are 4.691252835/8.911709881/5.235925894 s
  for direct and 6.425137982/4.235130774/8.916720416 s for erased.
  Median evaluator time is **79.6264/96.0689 microseconds** and Total time is
  **91.2902/107.9295 microseconds**. Three primed twenty-batch passes retain
  118,815/105,221 samples; their actual measured durations are
  3.364/3.812/3.700 s and 3.774/3.832/3.797 s. Both physical points agree,
  with relative complex differences **2.3019e-15** and **5.4482e-16**.
  Direct instruction/multiplication/addition counts are 26,169/19,029/15,158
  versus erased 31,999/23,657/17,588. The frozen binary, complete state hashes,
  counters and timing attempts are retained in the same diagnostic directory.
- The GL00/GL01 local construction-plus-projection subtotals are
  182.708/181.572 ms. Allocation, candidate certification/selection and cache
  work total 11.984/11.200 ms, or **6.56%/6.17%**. These include losing CFF
  generation once: 4.685/4.683 ms. Winning CFF generation takes
  20.870/20.991 ms. Both have ten admitted candidates across eight requests
  and 196 native rows in total. Summaries use the explicitly certified single
  graph-owned context of these legacy-forest cards. These are profiling
  diagnostics, separate from final unprofiled performance gates.
- The generic batch refinement now closes all eligible tensor leaves attached
  exclusively to the same Sum in one pass. It keeps that Sum's native arms,
  retains store references, and preserves the complete residual boundary.
  A partial-closure regression keeps two coupled sums separate in either
  operand order. The unchanged coordinate and graph invariants, plus the
  broader focused suite, pass **634/634 tests** in 58.225 s (3m12s build).
  Check39 passes before compilation. Physical replay and timing of the batch
  implementation remain pending.

- The batch replay improves GL00 preparation from **4.220910 to 1.068252 s**
  and GL01 from **4.212994 to 1.058098 s**, about four times faster in both.
  Tensor execution is 4.083161/4.162146 s and outputs are
  39,268,458/39,370,594 bytes. All 84 moved leaves and four root terms remain;
  the exact tests are unchanged. The combined preparation/execution stages
  improve by factors 1.654/1.634. These are single captured-input replays;
  complete generation and runtime of the batch revision remain to be measured.
  The immutable CLI hash is
  `d381b796830d86c8868b42e2d2b9e5ec781915956865040daca101b1ccf30dd8`,
  based on d21ce5c92 plus recorded patch
  `925dfab412485338691ac18d393e7156e3183b0e14d12496eeb8aefb0bda1e05`.
  The CLI build takes 4m00s; clippy again succeeds with the same three tuple
  warnings. Probe receipts are in `tensor_component_replay/prepared_nextest14`.
- GL262 with the preceding `40f7f703` CLI remains incomplete. It was explicitly
  terminated after 301.196892 s (exit -15, peak RSS 463,589,376 bytes), still
  inside the first forest node, with 3,038 production orientations. There is
  no evaluator, saved state or production forest export, hence no timing ratio
  or production-forest verification. The original baseline remains untouched.
- A fresh 124-sample GL262 profile (zero lost samples) locates 91.13% inclusive
  cost in final `collect_color` / representation collection, including 89.52%
  in symbolic zero tests. This is a different subpath from the earlier
  `collect_chains` bottleneck: early cograph simplification removed that work,
  but final color-tensor collection is still expensive. New existing-boundary
  debug events expose the exact cograph and mapped integrand before final
  color simplification. Their source-only capture will distinguish remaining
  scalar color invariants from open tensors before the next generic fix.

### 2026-09-13 — GL262 scalar color-invariant boundary

- Milestone `04501f7eba67f1387ca86d6a9e92fd223fae59e3` is committed and pushed
  as ValentinHirschi. Its immutable diagnostic CLI has SHA-256
  `7817960155b15016cbc869632b23ddf339eb4966157dce0cd4119f2cae4f33fb`;
  check40 and CLI build7 succeeded with debug assertions enabled.
- Its GL262 source-only run reached the requested pre-color boundary and was
  deliberately stopped after 26.1735 s, peak RSS 684,453,888 bytes. The cograph
  is 3,078 Atom bytes; the mapped input is 17,459,639 Atom bytes. Exact captured
  arguments show only `idx(2,cof(3))` and `cas(2,coad(8))`: one of each before
  mapping, 3,038 of each afterward. All color representation calls have arity
  one and are accounted for inside those scalar invariants; there are no open
  color tensor slots. Files, hashes and inventories are retained in
  `diagnostic_final_color_boundary/GL262_4d_source_capture_r1`.
- The existing color simplifier now applies its requested fundamental-dimension
  invariant substitution immediately after local color rewrites, before
  representation collection. This handles both preexisting and newly produced
  scalar invariants without distributing the accompanying numerator. The
  setting and final substitution remain intact. Focused correctness checks and
  a fresh full GL262 generation will validate this change; the source-only run
  is not a completed generation or a performance acceptance result.
- Check42 passes; the focused nextest suite passes **639/639 tests** in
  58.214 s after a 3m00s build. Five new color tests cover factorized spectators,
  generated and unsupported invariants, disabled trace/Fierz operations and
  idempotence. Clippy4 succeeds with the same three type-complexity warnings.
- CLI8 SHA-256
  `c43e77703ad0e23b0b06ba2c5a56606f90f77c0902e6d4e6205913e3330d8c37`
  builds in 58.85 s. Its first GL262 forest node completes in 55.392331 s;
  the next full-DOD0 projection takes 2.570277 s, and the DOD2 bubble projection
  takes 41.925 ms. This passes the previous first-node color bottleneck.
- The DOD2 bubble's final assembly then reaches a separate color-collection
  bottleneck. A fresh 116-sample profile (zero lost samples) attributes 96.55%
  to collection and 94.83% to statistical zero tests. Its 24,756,504-byte Atom
  contains genuine open adjoint slots: 19,344 two-argument `coad` calls, plus
  7,936 scalar adjoint invariants. Exact inputs and profile receipts are in
  `diagnostic_color_invariants/GL262_4d_direct_r1`. The run remains incomplete;
  the next investigation concerns factorized coefficient handling in the
  existing shared tensor collector, not another scalar-invariant assumption.

### 2026-09-13 — Factorized coefficients during tensor collection

- Milestone `c4acf58a0` is committed and pushed as ValentinHirschi. CLI8's
  DOD2 diagnostic was deliberately stopped once the new bottleneck was
  established: 412.881438 s wall time, exit -15, peak RSS 1,085,669,376 bytes.
  It produced no evaluator, saved state or production forest export.
- An artifact-only experiment on that exact mapped input protects maximal
  unselected composite coefficients before the existing polynomial collector.
  After the existing scalar invariant conversion, the input is 24,533,304 Atom
  bytes. Wrapping takes 552 ms, protection 173 ms, collection 992 ms and
  restoration 290 ms. Its 344 structurally interned aliases reduce the temporary
  expression to 1,043,079 bytes; the restored result is 23,982,068 bytes.
  Exact restoration before collection reproduces the full wrapped input in
  64 ms. The entire process, including parsing and canonical export, takes
  6.68 s. Source, binary/input hashes, receipts and seven tiny algebraic checks
  are retained in `color_coefficient_shielding`.
- Six tiny outputs exactly match the previous collector. The intentional
  boundary distinction is explicit: an opaque coefficient such as
  `(x+y)^2-x^2-2*x*y-y^2` remains factorized, whereas the previous collector
  could statistically recognize it as zero. Equal scalar terms still cancel
  on restoration; the tensor collector does not expand the numerator to prove
  scalar polynomial identities. Existing physical zero checks remain required.
- The shared `Collectable::collect_with_map` now uses this protection through
  the existing Symbolica `AliasedAtom` owner. It preserves complete selected
  tensors, avoids input-symbol capture, groups through the same collector, and
  restores all aliases before public tensor callbacks. No graph-specific rule,
  persistent table, setting or parallel mapper is introduced. Check43 passes
  in 11.86 s. Focused tests and complete physical generation remain pending.
- Check44 passes in 9.79 s; nextest16 passes **644/644 tests** in 58.558 s
  after a 1m12s build. Existing color polynomial/physical-zero snapshots and
  checks pass without expectation changes. Five new Spenso tests cover symbolic
  powers, repeated coefficients, input-symbol collisions, structural
  cancellation, opaque polynomial coefficients and complete callback payloads.
- CLI9 builds in 1m02s and is frozen with SHA-256
  `3e061d56ce138250b1c2cd707780abd296ec552ba6e7205c6f567316414f569a`,
  based on `c4acf58a0aca1f901d2cf231ef3d79a5abc283b0` plus recorded patch
  `df4ceedfd759fe4c8ec22a06f7581e2f0f824078be5102c93d7ebb7b61084b4c`.
  Fresh GL262 generation is running under the shared measurement lock; no
  generation/runtime gate is inferred from the captured-input improvement.
- Clippy5 succeeds with the same three type-complexity warnings. Check45 and
  integration build3 also pass (14.29 s / 3m24s). The new `uv` and `test_runs`
  binaries are preserved under `diagnostic_coefficient_protection/integration_binaries`
  with complete source/patch/binary receipts; these builds are not executions.
- The next CLI9 profile, on a DOD1 region rather than the DOD2 bubble, finds
  no collector or statistical-zero-test frames. Its sampled work is ordinary
  final substitutions (59.09%) and integrand addition/copying (22.73%). The
  DOD1 mapped expression is 72,021,365 Atom bytes. Forest traversal indices
  differ between runs, so stage comparisons use DOD and subgraph masks; the
  association remains separate from the still-pending production DOT export.
- Existing scalar-suite coverage is precisely 166 three-route cases across
  49 labels, GL00–GL48, plus one standalone sampling-scale test. This matches
  `test_LU_scalar_xs` and the historical suite inventory. It must not be reported
  as 167 distinct topologies. The complete existing suite remains required.

### 2026-09-13 — Route chain composition through the shared collector

- Milestone `89cda93ffce3ef4325e43e4143abd2d91b0b4958` is committed and pushed
  as ValentinHirschi. A later CLI9 sample, taken after the DOD1 pre-color event,
  differs from its earlier substitutions/addition sample: **94.07%** is in
  `Chain::collect_chains`' direct polynomial collector and **89.83%** in
  statistical zero tests (118 samples, zero lost). This caller bypassed the
  protected shared collector. The captured DOD1 input contains 14,608 chains
  and 29,216 color generators. CLI9 was deliberately stopped after 484.549 s
  without a complete evaluator or production forest export.
- Chain composition now collects its existing selected chain heads through
  `Collectable::collect_with_map`, unwraps them, and applies exactly the existing
  composition and normalization rules. This reuses coefficient protection
  without another engine, policy or setting. Check46 passes in 9.60 s.
- CLI10 builds in 59.30 s and is frozen with SHA-256
  `da042ba03c784fb3a66746c4d1acb6aae64e21f976655164822160c280cc6e85`,
  based on `89cda93ffce3ef4325e43e4143abd2d91b0b4958` plus patch
  `eb4d93a78761a4e017e521a7d43a25e8a03b2b34913dfa9ea5a603a7fb24f80b`.
  The next diagnostic allows 1,200 s for ordinary assembly/evaluator progress
  while retaining early termination on another conclusively profiled stall.
  Focused chain checks and completed physical timings remain pending.
- The original no-assertion GL262 localized-3D baseline was capped after
  **21,702.392314 s** (6h01m42s), exit -15, without a generation summary,
  evaluator or production forest export. Its loaded executable hash matched
  the immutable original `45df452b7c1069e51c1b39547572bd7add9244accda55b68cdbbf8565b54d911`.
  Final RSS was 3,213,508 KiB and peak RSS 4,787,520 KiB. Original cards,
  binary and termination receipts remain preserved; this is a lower bound on
  an incomplete attempt, not a completed baseline or a successful timing ratio.
  No persistent original-baseline process remains to compete with final runs.
- Check47 passes in 9.43 s; nextest17 passes **645/645 tests** in 58.450 s
  after a 49.79 s build. The added chain regression verifies ordered Dirac
  composition with an intact factorized/symbolic-power spectator, and unchanged
  Dirac payload when collecting the unrelated color representation. Existing
  snapshots and zero checks pass unchanged. Clippy6 succeeds in 21.88 s with
  the same three type-complexity warnings.

### 2026-09-13 — Complete-forest diagnostic and final validation preparation

- Milestone `eef89f3d22ecf285d59da6bd4cd108143d0d8b31` is committed and
  pushed as ValentinHirschi. An independent read-only review found no actionable
  defect in the scalar-invariant, coefficient-protection and chain-collection
  changes; this review does not replace the pending integration checks.
- CLI10 completes the previously obstructed DOD2 node in 94.634396 s, with
  71.499 ms in local projection. Its later DOD1 node completes in 309.192255 s,
  with 215.432 ms in local projection. A 106-sample transition profile has no
  collector or statistical-zero-test frames. A subsequent DOD0 term completes
  in 231.875341 s; its 123-sample profile attributes 79.67% to ordinary
  Schoonschip normalization, including 71.54% in dot normalization. Both
  profiles report zero lost samples. These are diagnostic stage measurements,
  not complete generation times or performance acceptance results.
- After six of eight terms completed and the seventh was progressing, the
  diagnostic's total allowance was extended explicitly from 1,200 to 2,400 s.
  The original limit and extension receipts remain part of the run record.
  Evaluator construction, pointwise evaluation and actual production forest
  export remain pending; no successful GL262 ratio is inferred.
- Check49 and integration build5 pass. Four immutable, assertion-enabled
  `dev-optim` test executables from clean revision `eef89f3d2` are frozen in
  `validation_candidate_eef89`, with exact binary sizes/hashes, compiler and
  lockfile identity in `final_build_manifest.json`. The complete harness
  dry-run selects 184 checks in six groups. No execution result is claimed by
  this build or dry-run. The verbatim plan prefix still has SHA-256
  `5c3ce40e19a807015e346df66904d6e770a7bb7ead1128c3e15437578851428c`.
- The physical harness now clears six inherited Spenso diagnostic/override
  variables and all three logging overrides; cards determine final logging.
  The correctness harness also clears the sparse-shadow override. Median gates
  are unchanged, but runtime repetition is additionally required when the
  envelope of all accepted pass means plus/minus twice their reported SEMs
  crosses 1.15. This is a conservative repeat trigger, not a confidence
  interval. Eleven synthetic controls, 27-card preparation, Ruff and syntax
  checks pass; results and exact policy are in `physical_protocol_validation`.
- Seven frozen-binary boundary checks (four UV composition, two API
  cut/threshold and one analytic renormalization check) are queued under the
  same measurement lock, behind CLI10. The full scalar and physical matrix
  and final repeated generation/runtime gates remain outstanding.
