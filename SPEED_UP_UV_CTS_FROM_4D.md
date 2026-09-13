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

### 2026-09-13 — Bulk final residue assembly and production-node inventory

- The two-point runtime harness is prepared and passes fifteen synthetic
  controls. Its 27 fresh generation cases yield 54 point-specific runtime rows;
  generation statistics count each state once. Base and 100-times-scaled
  points have separate priming, actual-duration requirements and performance
  gates. Existing diagnostic receipts are preserved. Policy, controls and
  prepared cards are in `physical_protocol_two_point_validation` and
  `optimized_physical_matrix_two_points`.
- CLI10 reached its explicitly authorized 2,400 s cap and exited -15 after
  2,400.264068 s, with peak RSS 4,167,876,608 bytes. It produced no evaluator,
  saved state or production forest export. The extension watcher encountered
  a process-disappearance race after SIGTERM; its `finally` cleanup resumed
  the supervisor correctly. Raw supervisor/watcher receipts retain that note.
- **Correction to the live progress inventory above:** the complete post-stop
  log contains nine completed nodes and ten starts. The previously described
  final node completed in 928.612672 s at 04:37:01.924 UTC, another completed
  in 2.302546 s, and the following projection was interrupted. The production
  traversal includes `current=18NfM, given=11vs, DOD0`, followed by
  `current=18OW0, given=18NfM`. The earlier eight-node assumption was inferred
  from the three initially identified regions, not certified by an export.
  Additional production nodes are under independent investigation. No test
  expectation or definitive forest count is changed on this incomplete record.
- A six-second final-assembly profile has 28 user-CPU samples, zero lost;
  92.86% is copying inside binary `Integrands::zip_add` / Symbolica addition.
  The existing native `Atom::add_many` performs one merge of normalized sums
  without distributing products. The existing `Integrands::zip_add` API now
  accepts a whole summand wave, validates every cut-key shape, and uses that
  bulk operation. Both projected assembly and direct selector materialization
  consume it; no alternate mapper or graph-specific rule is added.
- Check52 passes in 24.21 s after fixing two ordinary compile errors in checks
  50/51 (iterator reference depth and the infallible test parsing macro). Two
  new regressions check intact numerator powers, denominator buckets, exact
  cancellation and missing/extra cut keys. Focused nextest18 is running.
- The seven frozen `eef89` boundary checks pass: UV composition 4/4 in
  22.220 s, API cut/threshold 2/2 in 4.399 s, and analytic renormalization 1/1
  in 0.387 s. These binaries precede bulk assembly and do not certify that
  subsequent change. Large-Atom diagnostic logging and the interrupted final
  projection make CLI10's aggregate dispatch ratio indicative only; the final
  dispatch gate requires a separate complete cold profile with minimal payloads.
- Independent enumeration of all 1,023 nonempty internal-edge subsets confirms
  exactly four aggregate-admitted regions: the three connected regions and
  `U=[4,6,8,9,10,11,12,13]`, a two-loop disconnected union of the DOD -2
  top-quark hexagon and DOD +2 gluon bubble. The existing forest builder admits
  U by aggregate DOD zero, before canonical UV projection. Its component-parent
  rules produce ten nodes/nine arcs: empty once, bubble once, self-energy twice,
  U once and full graph five times. The two additional paths are bubble-to-U
  and bubble-to-U-to-full. Actual production export remains required. The new
  GL262 test's incomplete region assertion is corrected to include U; existing
  physics expectations, tolerances and forest semantics are unchanged.
- Nextest18 passes **647/647 tests** in 58.760 s after a 2m50s build. Check53
  for the matched CLI configuration passes in 8.89 s. CLI build11 is running;
  its prepared GL262 cards use an initial 2,400 s cap and exclude `dump` and
  `trace` tags while retaining stage counters.
- CLI10's matched completed-projection prefix has 125.921 ms construction plus
  19,650.092 ms projection; conservative dispatch/cache work is 741.048 ms
  (3.7472%, indicative only). It includes 100 preparations, 162 admitted
  candidates, 7,828 mapped native rows and 15 native generations. Template
  preparation takes 8.653 s, row mapping 1.947 s, component composition
  278.101 ms and sector composition 4.874 s; these are nested timings.
  The CFF cache records 247 hits/15 misses/no evictions. Row memos record
  115,364 hits/42,095 misses/28,777 evictions. Full and matched-prefix summaries
  remain preserved separately; the unfinished-call aggregate is not a gate.
- CLI11 build passes in 3m32s and its immutable executable has SHA-256
  `c76a89674c2d02df791890691cd975d801693072cce40012d0f2382af33060d9`.
  Its source revision, patch, compiler and settings are recorded under
  `diagnostic_batched_integrand_sum`. The GL262 diagnostic started at
  04:50:57.341817 UTC with a 2,400 s cap. Its initial log filter accidentally
  required explicitly false `dump`/`trace` fields, suppressing stage logs.
  This attempt continues with its original card and receipt preserved; a
  separate complete cold profile will use the corrected absent-field filter
  `[{#generation,#profile,!dump,!trace}]=debug`. No dispatch timing gate can
  be inferred from this uninstrumented attempt.
- Independent review finds no bulk-addition correctness blocker. Every full
  cut-key shape is validated, including zero summands; multiplication, powers
  and function arguments stay intact. Clippy7 passes with the same three
  tuple-complexity warnings in existing assignment/source preparation code.
  Formatting and integration check54 pass; the four updated integration
  executables will be rebuilt and frozen before broader execution.

### 2026-09-13 — Complete dispatch accounting and color payload isolation

- Bulk assembly is committed and pushed as `d5938769c`. Integration build6
  passes in 4m35s. Four assertion-enabled executables from that clean source
  are frozen in `validation_candidate_d593`; seven boundary checks are queued
  behind CLI11 under the common benchmark lock.
- CLI11's ten-minute stack sample has 107 user-CPU samples and zero lost.
  Color simplification accounts for 90.65% inclusive, with chain collection
  53.27%, chain normalization 26.17%, identity-chain collapse 18.69% and
  trace closure 17.76%. These overlap. The previous bulk-addition and
  polynomial-zero-test hotspots do not appear in this sample. The full
  generation remains pending; no overall speedup is inferred from one window.
- A shared color refinement uses the existing tensor collector to isolate
  complete color payloads before the existing fixed-point rewrite loop. The
  selector also retains generic chain/trace normalization, accepted untyped
  color heads and wrappers. Global metric normalization remains before and
  after that boundary for mixed tensors. The existing collector now presents
  positive integer powers as complete compressed monomials to callbacks;
  numerator powers are not expanded. New regressions cover generic chains,
  wrappers, mixed slots and complete powered callback payloads.
- Checks55/56 caught ordinary new-test parsing and AtomView conversion errors;
  check57 passes after those corrections. Focused nextest19 is running.
  Existing physics expectations and tolerances are unchanged.
- The previous dispatch subtotal omits unconditional raw physical-degree
  reporting and outer soft-routing preparation. Added timing events cover
  those boundaries and separate each raw CFF request into source rebuilding,
  degree/preparation, native generation and postprocessing. Outer selection
  records its complete elapsed time and winning proposal; all losing work
  remains charged. This changes telemetry only, not admission or selection.
- For a complete local-only run, the full dispatch denominator is the sum of
  completed nonroot forest-node times, which includes child projection and
  outer assembly but excludes root production CFF and subsequent amplitude
  evaluator work. The old construction-plus-projection denominator remains
  a separately labelled child-only diagnostic. Nested structural-certificate
  time is already part of allocation and must not be added again. Interrupted
  or dump-contaminated prefixes cannot certify the complete overhead gate.
- Nextest19 completes 650 tests in 57.472 s: 648 pass and two representation
  assertions fail. The new generic-chain fixture omitted the canonical cyclic
  wrapper. The existing three-loop color snapshot differs only by moving a
  numerical factor 16 into two coefficients (`16/128 = 1/8`), without changing
  their scalar polynomials. These routine assertions are corrected under the
  user's standing authorization; no physical expectation or tolerance changes.
- Further review identifies a factorization hazard in selecting all generic
  traces as collector variables: a power of a sum of color-free traces could
  distribute. The refinement separates factor-preserving generic chain/metric
  normalization from actual color collection, reusing the existing fixed-point
  reducer. A new exact no-op regression covers this boundary. Check58 passes
  in 3m17s, but its idenso-only package selection unnecessarily rebuilt a
  different dependency feature set; subsequent checks retain the full set.
- CLI11 advances from color rewriting to tensor construction. A second
  six-second profile has 129 samples and zero lost: 88.37% inclusive and
  78.29% self are in Linnet's `SmartEdgeVec::connect_identities`. Its sole-use
  edge-reference helper scans the entire involution after each join. The
  existing `Swap<EdgeIndex>` already updates the two exact owning half-edges
  in constant time, so that operation replaces the scan without adding a
  helper or changing merge flow. A new 48-case test compares every half-edge
  against the independent direct-involution operation, including paired and
  dangling moved edges and both merge-flow/source-order choices.
- The corrected color boundary now runs the same reducer on generic nodes
  separately, reaching a fixed point with metric normalization before and after
  actual-color collection. This also handles an identity chain whose resulting
  metric closes another chain. The new tests preserve a powered non-color trace
  sum both alone and alongside a color square, with trace evaluation on and off.
  Check59 exposed reserved-keyword parsing in the new fixtures; check60 passes
  in 1.98 s after that routine correction. Nextest20 is running with the new
  graph-join oracle included.
- Nextest20 completes 653 tests in 58.017 s: 652 pass, including the powered
  trace-sum and generic chain-closure regressions and the graph-join oracle.
  The remaining new fixture compares identical metrics with differently
  qualified index symbols; its input namespaces are corrected. Targeted
  nextest21 passes that check. Clippy8 passes in 26.88 s with the same three
  tuple-complexity warnings. The complete dispatch parser passes fifteen
  synthetic controls and documents exact stage/overlap accounting in
  `projection_profile_validation/README.md`.
- Targeted nextest22 also passes the metric-driven chain-closure fixture after
  index qualification. CLI check62 passes in 13.21 s; build12 is running.
- CLI11 reaches its unchanged 2,400 s cap. SIGTERM is requested at
  2,400.0830 s; final outer wall is 2,401.9509 s, exit -15, and peak RSS is
  20,023,398,400 bytes. It produces no saved state, evaluator or forest export.
  A separately recorded 30 GiB RSS guard was armed late in this diagnostic
  and exits without intervention. The earlier reference to a preexisting
  external diagnostic RAM guard was incorrect: the existing 30 GiB watchdog
  belonged to the correctness harness. Both supervisor and RSS-guard receipts
  remain preserved. This is another incomplete diagnostic, not a timing ratio.
- All seven frozen `d593` boundary checks pass: four UV composition checks in
  21.827 s, two API cut/threshold checks in 4.783 s and analytic renormalization
  in 0.393 s. This validates bulk residue assembly, before the subsequent
  color-payload and graph-join changes. Counter11 reproduces every recorded
  evaluator count in the known GL00 saved state exactly and leaves all state
  hashes unchanged; its validation runs after those checks under the same lock.
- CLI12 build passes in 3m30s. The immutable executable in
  `diagnostic_collected_color_payloads/immutable_binary/gammaloop` has SHA-256
  `bcc1f87e946ad77fed108fdb64704b638dd3d83f4793ac7f5c246c53063e6c03`.
  Its manifest binds revision `d5938769c`, source patch
  `c2b9f0b1efdffb9a5c6e96885ada136317d3c2ffed779ad3f40f2b0de0ae9ba6`,
  assertion-enabled build settings and compiler identity. The next GL262
  diagnostic uses the corrected minimal log filter, explicit cold-start
  attestation, a 2,400 s cap and a 30 GiB guard from launch. Its completion,
  actual forest export, pointwise checks and repeated timing gates remain pending.
- Integration check63 passes in 10.85 s and the four integration executables
  build in 3m50s from clean revision `9427afbef85325df90835172e5a951866493b6a7`.
  `validation_candidate_9427/final_build_manifest.json` binds their immutable
  copies, hashes, compiler and assertion-enabled settings. Fourteen UV,
  cut/threshold and analytic checks are queued behind CLI12 under the common
  benchmark lock; they are not yet reported as passed.
- CLI12 completes nine of ten forest nodes and reaches the final node at
  05:48:42 UTC. A six-second sample of its second DOD1 node contains 123
  samples with zero lost: 91.87% include final-integrand `expand_dots`, with
  network construction, Schoonschip normalization and execution below that
  boundary. The earlier color reducer and dense edge-reference scan are absent
  from this sample. The DOD1 node subsequently completes in 249.383 s.
  This diagnostic overlaps earlier compilation and includes sampling overhead;
  it is neither a final generation measurement nor a completed ratio.
- Read-only inspection finds that dot expansion currently constructs a metric
  network per replacement match. Symbolica offers a per-call RHS match cache,
  but its default size of zero does not establish repeated inputs in this run.
  Existing saved pre-color expressions are 58–231 MiB and precede the actual
  dot-mapping boundary. No speculative cache or replay of those large inputs is
  introduced on that evidence alone.
- All ten CLI12 forest nodes finish: 912.777859 s for the complete forest,
  including 19.830121 s for the root and 892.945160 s for nonroot nodes.
  Complete dispatch accounting is 6.150606 s, or 0.688800% of that nonroot
  denominator. Evaluator completion and the successful-run receipt are still
  required to certify the cold-run gate. The artifact parser passes 27 synthetic
  controls, including standalone capacity scopes and exact overlap rejection.
- The completed forest contains 27 bucket visits, 127 source requests/builds,
  127 template builds, 201 scored candidates and 8,494 selected native rows,
  each mapped once. The cold context makes 15 native CFF calls; all 15 canonical
  source/ordered-capacity projections are distinct, certifying 15 complete CFF
  keys once the cold successful-run receipt is available. Dropping capacities
  leaves 13 canonical parsed source descriptors. These are distinct from
  preparation and row misses: 257 preparation evictions and 60,769 row evictions
  prevent interpreting their miss counts as globally unique keys.
- Evaluator network parsing finishes in 210.181 s. The next six-second profile
  has 119 samples, zero lost, and 94.96% inclusive in operator merging's dense
  bit-vector union. Each operator island allocates a graph-sized mask and unions
  it into the deletion set. The replacement uses the existing disjoint group
  identities to mark only newly internal edges from their node crowns, then
  calls the existing node-identification API. It preserves every edge kind,
  original self-loop and intergroup connection. No new helper or symbolic
  representation is introduced. Independent review confirms the invariant;
  the regression compares complete graphs and slot orders with the existing
  dense helper, including idempotence, three islands and dangling/loop/crossing
  slots. Check64 passes in 10.15 s; focused nextest23 is running. CLI12 remains
  immutable and continues under its existing bounds.
- Nextest23 builds in 3m06s and passes all 654 focused checks in 57.988 s,
  including the new operator-merging oracle. Clippy9 passes in 22.50 s with
  the same three tuple-complexity warnings. CLI check65 passes in 11.09 s;
  build13 is running. Because operator merging is shared, paired 3D generation
  will use this new binary as well. The pending CLI12 result cannot establish
  a timing ratio for that later implementation.
- CLI13 builds in 3m24s and is frozen in
  `diagnostic_local_operator_merges/immutable_binary/gammaloop`, SHA-256
  `f9da9df2cee4077a121b2fb9d6c75e7c04a2a3a073562eae6733295be0d6b045`.
  Its manifest records base revision `9427afbef`, patch
  `920ad426ed129b441afec299ca15e97da82dc0210aaab04432e57270ad452b53`,
  compiler, flags, assertion-enabled profile and immutable byte count. The next
  diagnostic is queued after the existing frozen-run/check sequence; no new
  generation or runtime result is implied by this build.
- CLI12 reaches its 2,400 s cap during operator merging: termination is
  requested at 2,400.1950 s, final outer wall is 2,402.2633 s, exit -15,
  and peak RSS is 21,575,184,384 bytes. Its guard does not intervene. The
  forest is complete but evaluator preparation is incomplete, with no saved
  state, evaluator or production forest export. This remains an incomplete
  diagnostic; the full cold-run gate is not certified.
- All fourteen frozen `9427` boundary checks pass after CLI12 releases the
  lock: four UV composition checks in 21.827 s, seven cut/threshold checks in
  17.818 s, two API checks in 4.329 s and analytic renormalization in 0.411 s.
  These checks precede the operator-merging replacement.
- The operator-merging milestone is committed and pushed as
  `1c4131224ee10189c2a42659cf24d25701d695e2`, with ValentinHirschi as both
  author and committer. Counter13 reproduces the known GL00 evaluator counts
  without changing any state hash. CLI13 starts at 06:22:39.323946 UTC, with
  the same 2,400 s cap and its 30 GiB guard armed from launch.
- Integration check66 passes in 9.53 s and build8 passes in 3m24s. Four
  immutable executables from clean `1c4131224` are bound by
  `validation_candidate_1c413/final_build_manifest.json`; fourteen current
  boundary checks are queued behind CLI13. The 167-case scalar matrix, three
  physical route tests, final repeated generations and both runtime-point
  gates remain pending on this implementation.

### 2026-09-13 — Independent additive tensor preprocessing

- CLI13 completes all ten forest nodes in 913.929492 s: root 19.978591 s
  and nine nonroot nodes 893.948735 s. Full dispatch accounting is
  6.167046 s (0.68987% of nonroot generation), pending a successful complete
  generation receipt. It uses 13 canonical source descriptors, 15 distinct
  ordered-capacity source keys, 127 requests, 201 candidates and 8,494 mapped
  native rows. A full-ancestry comparison with CLI12 matches source/capacity,
  scored candidate, selected-row and outer-choice multisets; these counts alone
  do not certify exact coefficients or immutable assignment equality.
- CLI13 network parsing takes 213.163113 s, followed by 153.504 ms for 476
  aliases over 653 scalar terms (4,272,626 bytes). The next six-second profile
  has 97 samples with zero lost: 92.78% include boundary preparation and 85.57%
  include network extraction. Offline source resolution places 83 samples at
  the whole-network extraction of a ready sum boundary, not extraction of
  individual sum arms. The preceding operator-mask hotspot is absent from
  the >=1% report. CLI13 remains immutable under its original cap and guard.
- The next shared fix stays in EvaluatorStack::preprocess_atom: normalize once,
  contract each already-existing top-level additive term in its own network,
  restore its local aliases, and bulk-add scalar outputs before global evaluator
  optimization. No products, powers or nested sums are distributed. Independent
  review confirms that open tensor zeros remain invalid scalar outputs and
  aliases must resolve before combining independent stores. Three regressions
  cover the whole-network finite-component oracle, factorized scalar aliases,
  contraction-induced cancellation and open indices. Check67 passes in 10.05 s;
  focused nextest24 is building. A separate normalization interval and per-term
  indices preserve meaningful preprocessing timing accounting.
- Nextest24 builds in 3m25s and passes 37 of 38 selected checks. Its only
  failure is the newly added alias oracle comparing a tensor-library floating
  coefficient 4.0 structurally with integer 4. The assertion now requires exact
  zero after subtraction, with no tolerance and no factor expansion. This is
  a routine test representation correction authorized by the user's standing
  instruction. The whole-network oracle also covers the existing fallback and
  both contraction modes. Check68 passes in 8.16 s; broader nextest25 is running
  with the curated profile and zero retries. The initial nextest24 invocation
  used the default profile and retried its failure once; it is retained as a
  failed diagnostic, not clean coverage.
- CLI13 reaches its unchanged 2,400 s cap: termination is requested at
  2,400.179081 s, outer wall is 2,402.146918 s, exit -15, and peak RSS is
  21,493,309,440 bytes. Its 30 GiB guard does not intervene. Boundary
  preparation never completes, so no evaluator, saved state, production forest
  export or runtime measurement exists. The corresponding erased-3D card
  remains unrun. This is an incomplete diagnostic, not a successful timing ratio.
- The existing artifact parser now has a separate evaluator timing ledger;
  all 40 synthetic controls pass (27 previous controls plus 13 new ones).
  It accounts for normalization + per-term parsing + per-term execution +
  finalization/alias restoration/addition/destruction residual, validating term
  counts, ordering and complete successful execution. Boundary preparation is
  a subset of parsing; cumulative milestones are never summed as phase times.
  Historical unindexed parsing retains its normalization-inclusive meaning.
  The old parser is preserved and the new parser has SHA-256
  `6ca6934c6a367426ad0e5c5f91ea274c1d39f72a927c725dbd5789c4bcbedcb9`.
- All fourteen frozen `1c4131224` boundary checks pass: four UV composition
  checks in 22.184 s, seven cut/threshold checks in 18.088 s, two API checks
  in 4.427 s and analytic renormalization in 0.423 s. These validate the
  operator-merging revision, before additive preprocessing.
- Nextest25 builds in 2m51s and passes 306 of 307 broader UV, CFF, numerator
  and evaluator checks in 47.752 s. The same new alias fixture still fails:
  floating zero is not the canonical Atom::Zero used by is_zero. Comparing the
  complete original outputs confirms that only the two 4.0 versus 4 coefficient
  spellings differ. The independent trace oracle now uses the tensor library's
  explicit f64 coefficient representation and restores exact structural Atom
  equality; no physics value, tolerance or factorization changes. Focused
  nextest26 is rebuilding after check69. Clippy10 passes with the same three
  tuple-complexity warnings (2m39s including its wait for the build lock).
- Check69 passes in 7.73 s. Focused nextest26 builds in 2m27s and passes all
  five evaluator-preprocessing checks in 0.253 s, including the corrected exact
  alias oracle. Together with nextest25 this covers all 307 selected checks;
  the only intervening change was the new test's coefficient representation.
  No physical expectation or tolerance was changed. Final clippy11 is running.
- Clippy11 passes in 18.11 s with the same three tuple-complexity warnings.
- The additive-preprocessing milestone is committed and pushed as
  `1eb32575d9b16fee2910f3bb2c791fbae2eae5ba`, with ValentinHirschi as author
  and committer. CLI14 check passes in 9.441 s and build in 205.144 s. Its
  immutable clean-source binary has SHA-256
  `7aad642f9b065f89e8842d6943173f15aeeb17b24230779c7a050b75021829db`,
  1,575,509,560 bytes, in `diagnostic_additive_component_preprocessing`.
  Manifests record the empty source patch, exact flags, compiler and assertions.
- Counter14 links in 22.382 s with unchanged direct dependency hashes and is
  frozen with SHA-256
  `96038ccd0115814c67999f760dc71bbabde4dfa320ff5b9041d89bdb44e66b78`.
  Its known GL00 validation reproduces every evaluator-count row and preserves
  all saved-state hashes. CLI14's fresh GL262 direct run starts at
  07:19:22.549834 UTC with the same 2,400 s cap and 30 GiB guard from launch.
  It overlaps the following integration compilation and is a diagnostic,
  not a final timing run.
- Integration check70 passes in 10.32 s and build9 passes in 3m42s. Four
  immutable executables from clean `1eb32575d` are bound by
  `validation_candidate_1eb32/final_build_manifest.json`. The exact four,
  seven, two and one boundary-check groups are queued under the common lock
  after CLI14. The full scalar/physical matrix and repeated performance gates
  remain pending. CLI14 has completed six of ten forest nodes; evaluator
  preprocessing has not yet been reached.


### 2026-09-13 — Batched finite tensor boundary preparation

- CLI14 completes all ten GL262 forest nodes in 913.275 s, with the root
  taking 20.137 s. Evaluator normalization takes 4.654 s and finds 3,050
  top-level summands. The first summand takes 0.575 s for parsing/preparation
  (including 0.383 s for sixteen boundary closures) and 3.247 s for execution.
  The second summand remains in boundary preparation near the diagnostic cap;
  no completed evaluator or performance ratio is claimed.
- A separate six-second, 19 Hz user-CPU stack sample of CLI14's generation
  worker records 85 samples with none lost. Boundary preparation appears in
  98.82% of inclusive stacks; extraction in 42.35%, expression-tree building
  in 23.53%, and deletion in 17.65%. These overlapping stack percentages
  characterize this window, not the whole run. The immutable binary identity,
  thread, tool, commands and artifact hashes are recorded in
  `diagnostic_additive_component_preprocessing/GL262_4d_direct_r1/` under
  `perf_large_summand_manifest.json`.
- The existing finite-component boundary owner now selects disjoint outermost
  Sum scopes, including all eligible siblings, from one expression traversal.
  It appends small Product wrappers holding copies of ready leaf references,
  rewires exact half-edge IDs and removes obsolete endpoints once per wave.
  Original Sum shells, native arms, residual seams, tensor data and scalar
  spectators stay in place. Traced leaf loops preserve both axis positions.
  One operator merge follows the batch. The obsolete per-arm extraction path
  is removed and its meaningful comments are carried to the replacement.
- The existing independent coordinate oracle now also covers eight disjoint
  scopes with reused slot labels, four sibling targets, sums sharing a residual
  connection, nested waves and a closing leaf with an internal trace. Each
  checks exact finite-component values, graph consistency, store sharing and
  idempotence. No physical expectation, tolerance or graph numerator changes.
- Formatting and check71 pass (10.04 s). Nextest27 exits before compilation
  because `-j8` aliases test threads and duplicates the explicit serial setting;
  nextest28 uses CARGO_BUILD_JOBS=8 and a single test-thread argument. Results
  and the remaining generation/runtime gates are pending.

- CLI14 reaches its unchanged cap: SIGTERM at 2,400.101208 s, outer wall
  2,400.666167 s, exit -15 and peak RSS 6,287,552,512 bytes (5.86 GiB).
  The memory guard did not intervene. Term one remains unfinished; there is
  no saved evaluator, production forest export or runtime measurement.
- A source/listing audit confirms the full existing scalar LU matrix is 166
  three-route cases over 49 graph labels: 49 base, 90 quadratic numerator,
  ten quartic numerator and seventeen owned-dot/component probes. One extra
  sampling-scale check gives the harness's 167. Integrated UV and threshold
  counterterms are enabled in all three-route cases. Dedicated raised-line
  amplitude and cut acceptances are also selected separately. Numerator
  degree variants are not propagator-power variants, and these counts do not
  establish more than 100 distinct topologies. No authored case is omitted.

- Nextest28 builds in 3m09s and passes all 343 selected Spenso, local UV,
  energy-degree and evaluator-preprocessing unit checks in 37.620 s. The five
  added coordinate cases pass without expectation changes. Clippy12 passes
  in 1m01s including the build-lock wait, with the same three existing
  tuple-complexity warnings. Formatting and diff whitespace checks pass.
- All fourteen frozen `1eb32575d` boundary acceptances pass: UV composition
  four in 22.205 s, cut/threshold seven in 18.057 s, API two in 4.492 s and
  analytic renormalization one in 0.401 s. These certify the previous additive
  preprocessing revision, not yet this batched graph transformation. The
  complete scalar/physical matrix and repeated performance gates remain open.

- The batched finite-boundary milestone is committed and pushed as
  `ee23541cef755a8ca932e44cda22d507efd9be3e`, with ValentinHirschi as
  author and committer. CLI15 check passes in 9.153 s and build in 214.360 s.
  The clean-source immutable CLI has SHA-256
  `192e54fe3eb9a419ec1adfd5094c40cadeeeacb5c079f63da92bbb796b4fa4ad`,
  1,576,124,376 bytes. Counter15 links in 22.431 s with unchanged dependency
  hashes and reproduces every known GL00 evaluator count without changing
  saved-state hashes. Its SHA-256 is
  `030b667aabaa924f369551b880a1005f9c8f5c140ff81e204d4e2a839c8c279a`,
  1,574,466,232 bytes. Both are in `diagnostic_batched_boundary_scopes`.
- To establish a current representation reference, CLI15 runs GL262 erased
  3D first, then direct 4D, with the same 2,400 s cap and 30 GiB guard for
  each fresh process. Erased 3D starts at 08:09:02.017811 UTC. Its root takes
  19.918 s, a DOD1 node 127.555 s and the raw overall DOD0 node 96.318 s.
  This is still a partial diagnostic; no total generation or runtime ratio
  is inferred. It overlaps the following integration build, recorded explicitly.
- Integration check72 passes in 9.49 s and build10 in 4m36s. Four frozen
  executables from clean `ee23541ce` are bound by
  `validation_candidate_ee235/final_build_manifest.json`; the exact fourteen
  boundary acceptances are queued after the erased-3D run under the common
  lock. The scalar/physical matrix and final repeated timing gates are pending.
- The scalar-suite audit also distinguishes orchestrators: its 166 three-route
  cases use the default HedgePoset implementation. The selected physical
  route tests use Compare, while the focused core test
  `scalars_integrated_cts_compare_legacy_and_hedge_poset` exercises both for
  integrated counterterms. This does not claim the entire scalar matrix has
  been repeated under both orchestrators.

- A separate six-second, 19 Hz user-CPU sample of CLI15 erased 3D, during
  nested DOD1 forest-node assembly after its Taylor output, records 102 samples
  with none lost. All sampled stacks pass through final integrand color
  simplification; 93.14% pass through DotNormalizer, including 48.04% in repeated
  redundant-metric replacement and 44.12% in the asymmetric-vector replacement
  traversal. These inclusive percentages overlap and are not whole-run shares.
  `perf_nested_assembly_manifest.json` binds the worker, binary and artifacts.
  This identifies another shared normalization cost for review; the pending
  direct-4D run must establish whether it needs further refinement. No metric
  normalization rule or numerical expectation has been changed.

- The selective evaluator-input capture filter is validated against the existing
  production parser and static callsite gate with six controls. Only ordinary
  profile fields and the normalized-input dump are evaluated; UV dumps, network
  DOT, trace and missing-term events remain unevaluated. The extra direct card
  has SHA-256 `853dfcc627bf5646cb1f3a09d75b8b90170d0779f7bb51c1be4e7c2dc441443c`.
  The minimal card is preserved. A capture run is diagnostic only and cannot
  certify dispatch overhead or final timing gates.
- CLI15 erased 3D completes the DOD2 bubble in 169.318 s, another DOD0 node
  in 141.971 s, the next DOD0 node in 213.334 s and the nested DOD1 node
  in 483.142 s. Seven of ten nodes are complete at 08:29:54.473 UTC. These
  remain partial, profiled measurements; complete generation is still pending.
- The shared metric-normalization review preserves existing stage order. In
  particular, metric-power normalization precedes metric tracing; exchanging
  their order can change g(i,i)^2 from D to D^2 under the current conventions.
  A possible narrow improvement is conservative structural gating of the
  unchanged root-level replacement rules, retaining all stages and the
  redundant-metric fixed point. No replacement rule has been changed yet.


### 2026-09-13 — Structural guards for metric normalization

- The two profiled DotNormalizer stages now use the existing top-down
  replace_map traversal, conservative symbol-tag guards and their unchanged
  replacement rules restricted to the current tree node. The asymmetric-vector
  guard scans all direct arguments rather than assuming a final slot: symmetric
  function heads can reorder arguments. The redundant-metric guard admits both
  representation tags and dind wrappers without requiring registered reps or
  imposing additional arities. False-positive eligibility is harmless because
  the original matcher still decides the replacement.
- The redundant stage retains the existing repeat_map fixed point. Vector
  powers, metric powers and traces retain their original separate passes and
  order, including the original negative-power filters. There is no new AST,
  cache, helper type, graph-specific condition or rule simplification. All
  existing meaningful comments are retained; stale extraction wording in the
  evaluator and boundary test is updated to describe current preparation.
- A new exact differential test restores the original unrestricted replacement
  settings as a test-only oracle. It covers 455 wrapped cases, including nested
  vectors, bare rank1 fallback, dual slots, unregistered tagged representations,
  symmetric and Linear heads, unequal dimensions and integer/symbolic/rational
  powers. Explicit metric-trace ordering and negative-odd-power checks plus a
  large factorized coefficient witness preserve the relevant invariants without
  expanding graph numerators.
- Formatting and check73 pass (8.84 s). Nextest29 and clippy13 are running.
  CLI15 erased 3D and its queued fourteen boundary checks continue against
  their immutable ee23541ce binaries; the new normalization code is not mixed
  into those runs. No final performance gate is claimed yet.

- Nextest29 builds in 2m27s and passes all 427 selected Idenso, CFF generation,
  UV, energy-degree and evaluator-preprocessing tests in 46.411 s, including
  the 455-case exact differential oracle and the existing 62-row exhaustive
  allocation reference. Clippy13 passes in 2m45s including its build-lock wait,
  with the same three existing tuple-complexity warnings. No expectations or
  tolerances needed correction. Formatting and whitespace checks pass.
- An independent read-only benchmark-harness audit confirms the default
  27-generation, 54-point-row and fifteen-gate arithmetic, actual duration and
  sample-count reconstruction, priming and uncertainty handling. A narrow
  per-stage resume check is being added to the existing replay owner so an
  unrelated successful receipt cannot be reused with a different command,
  card or provenance. Final runs still require the full population, verified
  assertions-enabled build and fresh unprofiled output directories.
- The unlaunched CLI15 direct cards are preserved as prepared, not run. The
  next direct diagnostic will use the normalization revision after its frozen
  CLI/counter are ready; erased 3D will then be remeasured at that revision.
  No ratio will mix CLI15 and the following binary. The current erased run
  and its queued snapshot-specific boundary checks continue unchanged.


- The guarded-normalization milestone is committed and pushed as
  `e69f6153d8f991af85b3f83d479be9a750a9e6be`. CLI16 check passes in
  10.285 s and build in 144.180 s. Its immutable binary has SHA-256
  `65e0576e6a7e3d61c1af3eb21c15c1407cf9211e4c5e1ef8b61983214b84f485`,
  1,576,695,904 bytes. Counter16 links in 22.743 s with unchanged dependency
  hashes, reproduces every known GL00 evaluator count and preserves its state
  hashes. Counter SHA-256 is
  `45f3c8c51a57978920b888fa2b2ce4c78cde113fbbdc5fac1bf38cae2b4e94f6`,
  1,575,041,864 bytes. Artifacts are in `diagnostic_root_guarded_metrics`.
- CLI15 erased 3D completes all ten GL262 forest nodes in 2,328.309505 s;
  the final node takes 650.886710 s. Evaluator parsing begins at
  08:48:54.650 UTC, eight seconds before the unchanged cap. The run terminates
  with exit -15 at 2,400.536083 s and peak RSS 4,628,992,000 bytes. No saved
  evaluator, state, forest export or runtime measurement is produced. The
  historical CLI14-direct/CLI15-erased semantic forest paths all match, but
  their separate revisions, compilation overlap and capped outcomes preclude
  any final generation or runtime ratio.
- All fourteen frozen ee235 boundary checks pass: UV composition four in
  22.481 s, cut/threshold seven in 18.184 s, API two in 4.465 s and analytic
  one in 0.409 s. Guarded-normalization integration check74 passes in 9.93 s
  and build11 in 151.671653 s. Four immutable executables are bound in
  `validation_candidate_e69f/final_build_manifest.json`; their fourteen checks
  are queued after CLI16 direct under the measurement lock.
- CLI16 direct begins at 08:52:04.791648 UTC with the same 2,400 s cap and
  30 GiB guard. Its validated selective input capture is diagnostic only.
  Integration build11 and subsequent checks overlap this run explicitly;
  no final timing or dispatch gate is inferred from it. Nine forest nodes
  are complete by 09:02:12.373 UTC, with evaluator construction still pending.
- Successful benchmark-receipt reuse now requires exact command, card hash,
  complete provenance, worker limits and measurement-lock identity. Fifteen
  synthetic controls pass without launching subprocesses, including missing
  identity fields and preserved failed/interrupted receipts. Ruff and Python
  compilation pass. The ignored helper SHA-256 is
  `c700369d89e5669a2e455ff0aa3b3bdbf1c84af2536fdde54f953543019728fa`;
  the exact patch and controls are in `receipt_resume_validation`.

### 2026-09-13 — Audit of compressed powers and certificate ordering

- Independent source review confirms the separate cache identities, exact signed
  bindings, graph-wide lifetime, retention limits and deterministic candidate
  policy. It finds two requirements that need correction before final gates.
  First, the initial 64 remembered offset contexts miss a packed repetition
  tail whose transient fills that table. A half-million power with 64 eligible
  occurrences consequently materializes nearly all factors. Second, the complete
  mapper diagonal certificate currently follows native candidate generation,
  although canonical source equality, signed occurrence eligibility and bound
  checks precede it. Existing winning outputs fail closed, but this ordering
  does not meet the plan.
- Repeat allocation now uses one deterministic checkpoint with a doubling
  interval. It finds cycles after long transients without increasing retention
  or changing assignment choices. The existing large-power test now covers
  5, 64, 65, 128 and 257 occurrences, exact ordered capacities for both policies,
  independent mapped-degree analysis, diagonal collapse and an O(p) prepared
  leaf bound. Graph numerators remain factorized.
- Every exact CFF candidate now prepares or reuses its certified immutable
  mapping template before native generation. The winner carries that same
  template into row mapping; no second preparation occurs afterwards. Existing
  cache/oracle callers migrate to the graph projection context. Selection
  records winning and losing template build times separately, and excludes
  already-accounted template build/cache costs from its residual interval.
  All losing preparation remains dispatch overhead. The parser is being updated
  to reconcile both this ordering and preserved historical logs.
- Check75 catches one unmigrated paired-mapper test binding before compilation;
  that mechanical caller correction changes no assertion. Check76 passes in
  10.87 s. Further accounting/readability refinements, focused nextest coverage
  and clippy are pending. The frozen CLI16 run is unaffected by these edits.
  Complete physical/scalar acceptance, final repeated generation/runtime gates,
  production GL262 forest export and final documentation remain required.


- Check78 passes in 9.94 s. Nextest30 builds in 2m58s and passes all 289
  selected CFF, UV, energy-degree and exact-source reconstruction checks in
  47.146 s, including the expanded compressed-power cases. Independent review
  confirms the full diagonal certificate now precedes native generation, that
  retention cannot affect candidate admission, and that the timing partition
  charges discarded work once. The template mapper remains private: external
  test oracles obtain their independent mapper through the existing source API.
  Check79 passes in 9.18 s; nextest31 builds in 2m50s and passes the 39 affected
  source tests in 2.253 s. Clippy14 passes in 20.66 s with the same three existing
  tuple-complexity warnings. Formatting and diff whitespace checks pass.
- The profiling parser now reconciles one template/score pair per admitted
  candidate before selection, exact winning/losing build totals, and zero-row
  discarded winners. Fifty-six controls pass, with Ruff and Python compilation.
  Completed CLI14 and actual CLI16 prefixes reproduce the old parser's dispatch,
  count and evaluator semantics. Parser SHA-256:
  `aa73853fb322125c2f2a00906d1b9cbdd1049e1d0c7c67a3e9ee7349ef6b24e5`.
  Evidence is in `projection_profile_validation/candidate_template_schema_final`.
- All fourteen e69f boundary checks pass: UV composition 22.352 s, cut/threshold
  18.459 s, API 4.507 s and analytic renormalization 0.390 s. These remain
  snapshot-specific boundary results; the full scalar/physical matrix is pending.
- CLI16 completes its ten-node forest in 920.816712 s, with final node
  314.133087 s, then evaluator normalization in 4.655648 s. Its factorized input
  capture is 678,173,112 JSONL bytes and contains 3,050 top-level summands; the
  exact offset/hash manifest is `normalized_input_capture.json`. Term zero closes
  sixteen boundaries in 72.227 ms, parses/prepares in 269.323 ms and executes
  in 3.200217 s. Term one closes 15,904 boundaries in 33.481867 s and completes
  parsing/preparation in 61.401224 s, passing the previous stalled stage.
- Term-one execution then reaches the existing 30 GiB guard. It requests SIGTERM
  at 09:10:25.334 UTC after observing RSS 35,930,345,472 bytes; two-second sampling
  permits that overshoot. The driver reports peak RSS 35,785,789,440 bytes,
  elapsed 1,103.731767 s and exit -15. No completed evaluator, state, production
  forest export or runtime exists. These are separate observed peaks, not a
  claimed hard RSS ceiling. The next investigation replays the preserved input
  through existing tensor APIs under the common lock, rather than repeating
  the forest. The planned CLI16 erased generation is held; all final gates
  and matching-revision measurements remain open.


### 2026-09-13 — Compact ready-operation membership

- The exact CLI16 normalized Atom was replayed through the existing public
  tensor pipeline, without redoing its forest or expanding numerators. The
  generic largest-summand selector chose reparsed index 4 (127,033,485 Atom
  bytes); this is not asserted to be original CLI term one because canonical
  parser ordering can change. Its prepared network contains 4,317,091 nodes,
  13,662,997 half-edges and 779,720 ready Products, with leaf slot ranks at most
  three. The complete input still has 3,050 existing top-level summands.
- Parsing the selected network takes 60.033 s; 29,920 closed tensor boundaries
  take 48.831 s. At primary execution, RSS rises from 6.44 GiB to 30.16 GiB and
  the 0.2-second RSS supervisor terminates the process with exit -15. Its 740
  user-CPU profile samples end with 180 consecutive samples in ready-operation
  enumeration, including `BitVec::repeat` allocation frames. No initial batch
  completion or tensor execution is reached. Using pre-execution counts,
  retaining a dense graph mask for every ready Product would require about
  1.211 TiB; this is a scale estimate, since initial operator merging can change
  those counts, not an observed complete allocation.
- Raw stacks, sampled resources, receipt, exact input/binary/library hashes and
  analysis are retained in `captured_tensor_replay/largest_cli16`. Version two
  of the standalone probe reuses one empty mask during diagnostic inventory,
  removing its two approximately 79-second inventory overheads; its tensor
  smoke passes. Neither capped replay establishes acceptance timing.
- The shared operation payload now retains sorted visible half-edge lists.
  All batch admission, profiling counts, sequential/parallel consumers and
  partial rewrite callers migrate together. Existing node identification is
  generalized to an ordered iterator and membership predicate, so collapsing
  each operation also avoids constructing a dense temporary mask. Extraction
  alone constructs one complete-extent subset. Sorted half-edge order preserves
  node representative selection; hidden-edge filtering and self-edge marking
  retain their previous semantics. No contraction policy or physical expression
  changes. Focused checks, replay and broader acceptance remain pending.

- Check80 detects missing `Inclusion` imports in the migrated test modules before
  compilation; those routine imports are corrected without changing assertions.
  Check81 passes in 18.21 s. The first nextest command is rejected for an
  incorrectly joined CLI flag before compiling or running tests; the corrected
  nextest33 invocation builds in 3m58s and passes all 603 selected tests across
  Linnet, Spenso and GammaLoop in 56.927 s. This includes the 4,096-operation
  storage bound, dense/compact hidden-edge and representative-node oracles,
  existing tensor strategy checks and 289 UV/CFF/source/allocation checks.
- Clippy15 passes in 33.88 s with the same three existing tuple-complexity
  warnings. Formatting, diff whitespace and the verbatim plan-prefix hash pass.
  Independent review finds no change to visible incidence, operation admission,
  traversal order or the existing deferred-deletion lifetime. CLI17 and a new
  provenance-bound replay are prepared next; no successful large-summand result
  or final performance ratio is claimed at this milestone.

- CLI17 is frozen from clean `947c6acef7b85aeeac86e54e9b099be7ed4678eb`:
  check 17.41 s, build 259.40 s, binary SHA-256
  `79645f8c4296abf1dd8b15541a61f4e91e48b5f9b46ac52ca507cf5329eb1836`
  (1,576,768,728 bytes). The artifact freezer required metadata-only finalization
  because Cargo reports host and runtime Symbolica artifacts separately; exact
  owning dependency fingerprints resolve the runtime libraries, without another
  Cargo build. Receipts are in `diagnostic_compact_ready_operations`.
- The unchanged v2 probe source links against the five exact CLI17 runtime
  libraries as v3; metadata check 0.161 s, link 21.414 s, frozen binary SHA-256
  `6b944423fcd6330fcfd8dd30114021085205d2d382e68ae0f563c12027b96c0f`.
  Its unprofiled tensor smoke has the same 36 closed boundaries and 11,996,015-byte
  final scalar as v2; matching output size is not an exact equality proof.
- The first v3 largest replay is preserved as a supervisor failure: a valid JSON
  scalar profiling line was incorrectly assumed to be an object. It stopped
  before primary execution. The corrected runner passes mixed object/scalar/
  null/list controls and a profiled tensor smoke, clears the actual sparse-tensor
  override, and terminates/waits its owned child group and perf if monitoring
  raises. The repeated run has a fresh directory and the same immutable binary.
- `captured_tensor_replay_v3/largest_cli17_profile2` reproduces the same input
  and prepared inventory. It completes 29,920 boundaries in 44.510 s and reaches
  primary execution at 151.818 s. Its first actual plan contains 779,720 operations
  and 4,695,382 included half-edges, with global extent 13,662,997; first-wave
  execution takes 1.895 s. Compact membership passes the original unfinished
  planning failure and thirteen complete batches (0–12) execute. After batch 12,
  the live graph has about 25,000 nodes but the global half-edge extent remains
  unchanged until deferred deletion.
- A later memory increase reaches the 30 GiB guard at 421.513 s, observing
  32,364,167,168 bytes; exit -15, with primary execution still incomplete.
  The 10,800-sample profile contains real contraction work; its final samples
  primarily construct the next traversal/cache roots. This later failure is
  distinct from retaining a dense mask for every ready operation. The profile
  does not establish how much tensor-store payload is still live. A separately
  guarded 64 GiB public-API inventory replay is prepared to measure that, without
  changing production execution. Profiling and build overlap preclude acceptance
  timing claims from these diagnostics.
- Integration check82 passes (14.78 s; wrapper 14.832 s), build12 passes
  (5m10s; wrapper 310.928 s), and all four executables are frozen under
  `validation_candidate_947c6`. All fourteen boundary checks pass: UV composition
  21.988 s, cut/threshold 18.307 s, API 4.472 s and analytic renormalization
  0.391 s. These are nextest summary durations and snapshot-specific correctness
  results, not the complete 184-test final matrix.
- Counter17 is frozen with SHA-256
  `d8e69655dff3ce3f877cac55ec6e5cae5f42e84c06f6ffb6aa51d72bf80063f7`.
  Integration compilation had selected the core's `eval` dependency feature
  and an API without `ufo_support`; its API is not mixed with CLI17's core.
  A separately recorded check/build restores the exact CLI17 API/core pairing,
  verifies all five previously frozen shared-library hashes, and records the
  restored API hash. Its known GL00 counter validation passes, with unchanged
  state and identical reference operation-count rows.
- Full GL262 direct/erased generation cards are prepared but unlaunched. No
  completed GL262 evaluator/state, numerical timing, production forest export
  or final repeated generation/runtime ratio is yet available. Final scalar/
  physical coverage and all original acceptance gates remain open.


### 2026-09-13 — Reclaiming completed tensor intermediates

- The v4 public-API inventory replay completes the largest captured GL262
  summand using unchanged CLI17 production code (`947c6ace`). It executes all
  23 native batches, exits successfully in 588.166 s, and reaches a measured
  peak of 56,981,630,976 bytes (53.068 GiB) under its diagnostic 64 GiB guard.
  The final scalar is 3,389,407,083 Atom bytes, emitted at 565.051 s. This is
  one selected summand, not a completed graph evaluator or saved state.
- Terminal stored Atom payload is 28,439,588,507 bytes: 3,389,407,083 live and
  25,050,181,424 dead (88.082%). All 21,542,810,048 tensor Atom bytes are dead,
  across 3,852,060 tensor slots. Scalar slots contain 3,507,371,376 dead bytes
  and the live final scalar; this source creates no aliases. Counts include
  Atom payload per stored entry, excluding containers, concrete tensor data,
  graph memory and allocator overhead; they are not interchangeable with RSS.
- Small independent controls validate the inventory's complete leaf-reference
  walk and pin the original scalar definitions needed for alias restoration.
  The tensor control has 95,538,693 stored Atom bytes, of which 11,996,015
  remain live. An alias control retains its 4,817-byte definition and 10-byte
  handle. Exact identities and full stage records are retained under
  `captured_tensor_replay_v4`, including `inventory_controls.json` and
  `largest_cli17_64g/analysis.json`.
- The profile has 14,432 samples, including 4,248 cache-root frames and 4,552
  contraction-pair selection frames. These overlapping stack categories are
  diagnostic attribution, not additive generation timings or acceptance ratios.
- This evidence motivates a narrow shared-owner change: reclaim unused tensors
  after complete contraction waves, remapping all live tensor leaf variants.
  Scalar indices and alias definitions remain stable. Root alignment also uses
  the existing traversal's discovery order without constructing an extra
  child-vector forest, and drops the traversal before converting the graph's
  node store. The existing conversion is retained even when roots do not change,
  since it also normalizes sibling ordering. Implementation and validation of
  this next milestone are in progress; no new performance ratio is claimed.

- Implementation now extends the existing tensor-reference walker used for
  index shifting and parallel result rebasing. The owned `NetworkStore` uses
  move-only retention; borrowed overlays defer reclamation to their joined
  wave. All execution strategies reclaim only after their complete graph is
  restored, with a final pass after terminal materialization. Scalar vectors
  and alias definitions are unchanged. Tests exercise all six leaf variants,
  repeated/shared references, non-Clone tensor payloads, overlays, aliases,
  sequential/ref/extracted/parallel execution and self-loop-only completion.
- Check83 identified associated-type normalization at concrete `ExecuteOp`
  supertrait bounds. The affected concrete store and overlay callers now state
  their existing tensor/scalar equalities explicitly; check84 passes in
  11.15 s. No expectation or physical tolerance changed. The selected nextest
  suite is compiling; its result remains pending.

- Nextest34 passes all 608 selected Linnet, Spenso, CFF, UV, energy-degree and
  exact-source tests (3m14s build; 52.179 s tests; 467 skipped by selection,
  including 15 by the existing profile). Clippy16 passes in 24.81 s with the
  same three existing type-complexity warnings. Formatting and diff checks
  pass, and the original plan-prefix SHA-256 remains unchanged.
- The v5 reference probe is frozen against the exact CLI17 libraries, before
  rebuilding them. It can stream a portable scalar export and compare two
  imports exactly in one Symbolica state. Small controls with different symbol
  introduction orders produce different serialized file digests but equal
  imported Atoms; internal symbol-ID hashes are not used as an equality proof.
  The largest reference export is running separately under its existing guard.
  CLI18, its matched counter and the identical candidate probe are prepared
  for rebuilding after this milestone. Full graph generation, saved evaluators,
  pointwise physical checks and final timing gates remain pending.

- Milestone eighteen is committed and pushed as
  `d6e4c13d9b9b301a9e9201875446cbc199028d56`, with ValentinHirschi as author and
  committer. CLI18 is frozen from that clean revision: check 12.378 s, build
  217.802 s; SHA-256
  `22b6560abceb9e65bbaec193bf351e8807829c682385b53a7a421e07d675fc0c`,
  1,577,498,240 bytes. The exact API/core and six-library identities were
  captured before integration compilation changed dependency feature selection.
- Counter18 is frozen with SHA-256
  `1228640d14286b58a7a7d29e99371a8ca74c1e5b3aa5052979a4b9e516b5fa77`.
  Its known-GL00 validation preserves the saved state and reproduces the
  reference operation-count row. No API/core restoration build was needed.
- The unprofiled CLI17 reference export completes in 565.912 s with peak
  52.9368 GiB, reproducing the complete v4 scalar/store inventory. Its portable
  scalar file has 3,389,426,376 bytes and artifact SHA-256
  `66240aab1ee9115af7d8e0c4e5c5291b2cf74ef74a931efe07ce09db84e7fce0`.
  The file digest identifies the artifact, not a cross-process symbolic proof.
  This remains a single-summand correctness diagnostic; concurrent compilation
  and its export work exclude acceptance timing claims.
- The identical v5 probe source is linked against CLI18 (metadata 0.1653 s,
  link 20.2991 s); candidate binary SHA-256
  `d90666ac6278aa9fac0d79ebc9bf405b9800061aa912d5cae3b4ce4e2de11050`.
  Small tensor and alias controls compare exactly after importing both route
  exports into one Symbolica state. The tensor control has no retained tensors;
  its 11,996,015-byte final scalar and 48,309,278 dead scalar bytes are unchanged,
  while the previous 35,233,400 dead tensor bytes are removed. The largest
  candidate replay and its same-state comparison are in progress.
- Integration check85 passes (9.73 s; wrapper 9.788 s), and build13 passes
  (3m37s; wrapper 218.308 s). Four immutable binaries and complete receipts are
  frozen under `validation_candidate_d6e4c`; execution of the fourteen boundary
  checks and the complete final matrix is pending. Builds and replays use
  separate locks; all diagnostic executions and numerical benchmarks remain
  serialized by the common measurement lock.

- The largest CLI18 candidate replay completes successfully in 509.560 s, with
  measured peak RSS 28,444,950,528 bytes (26.4914 GiB), compared with the frozen
  CLI17 reference's 565.912 s and 52.9368 GiB. Its final tensor store is empty:
  precisely 21,542,810,048 dead tensor Atom bytes are removed. The live scalar
  remains 3,389,407,083 Atom bytes and the retained dead scalar payload remains
  3,507,371,376 bytes. This comparison concerns one captured summand and includes
  diagnostic preparation/export work; it is not a final generation ratio.
- Exact import comparison passes after importing both portable scalar files
  into one Symbolica state (153.407 s), without expanding graph numerators.
  Both exports also have identical bytes (`cmp --silent` succeeds). The complete
  replay, immutable-library and control evidence is recorded in
  `diagnostic_tensor_store_reclamation/replay_validation.json`.
- All fourteen d6e4c integration boundary tests pass: four UV composition checks,
  seven cut/threshold checks, two threshold API checks, and the analytic quark
  finite-part check. Nextest execution times are respectively 22.281, 18.070,
  4.423 and 0.383 s. All eight list/run receipts exit successfully and identify
  the same frozen assertion-enabled build. The full scalar/physical matrix
  remains pending.
- Full GL262 direct generation is now running from immutable CLI18, under the
  common measurement lock. Both direct and erased cards retain matching eager,
  uncompiled local-subtraction settings. Their diagnostic supervision permits
  64 GiB RSS and 7,200 s, allowing headroom above the completed isolated replay
  for the remaining summands and evaluator optimization. These are diagnostic
  guards, not changed acceptance criteria; previous capped failures remain
  preserved. The erased run will follow serially. No successful full GL262
  state, runtime ratio or completed performance gate is claimed yet.

- User steering: "Continue as planned, but you can increase the memory guard to
  100 GB for GL262, although still try to limit the memory load as much as you
  can, but also without introducing crazy complications. Generation and runtime
  optimizations are most important, but often memory optimization helps with
  those objective too." Future GL262 supervision may therefore use
  100,000,000,000 bytes. The currently running direct process stays uninterrupted
  under its already-active 64 GiB guard; its memory remains well below that
  limit. This permission changes diagnostic headroom only, not performance or
  pointwise acceptance thresholds. Subsequent final paired runs will use the
  same supervision limit.

- The existing process-tree watchdog now accepts an explicitly requested positive
  finite memory cap while preserving its 30 GB default and scalar recipe.
  The final correctness harness forwards an explicit `--memory-limit-gb` and
  records it in provenance; the physical GL262 group can use 100 GB without
  numerical-code changes or graph-specific branches. Eighteen tiny execution
  and argument controls pass, including 100 GB, invalid/nonfinite limits and
  overflow-safe byte conversion. Ruff, Python compilation and diff checks pass.
  Results are retained in `watchdog_memory_cap_validation/results.json`.
- A 20-second, 961-sample profile of the running CLI18 nested forest node places
  approximately 99% of that interval in final-integrand dot normalization,
  including 69% in `normalize_dots`. These are inclusive sampled intervals,
  not additive stages or whole-generation percentages. Remaining root-matcher
  guards are under read-only review; no numerical implementation change is
  made while the full generation result remains pending.

- CLI18 completes all ten GL262 forest nodes. Completed nonroot local-UV time is
  879.897 s; full dispatch/selection/cache accounting is 23.444 s (2.6644%).
  All losing-candidate preparation is included: 74 additional candidate
  templates take 17.164 s, bringing candidate-template builds to 201. Structural
  work remains 27 canonical buckets, 17 component waves, 127 source requests,
  201 scored candidates, 15 native CFF generations and 8,494 selected row
  mappings. Fifteen distinct source/ordered-capacity projections and fifteen
  native calls give matching uniqueness bounds. These completed-forest
  observations await a successful whole-generation receipt before certification.
- Full evaluator preprocessing passes the first two large summands: preparation
  and contraction respectively take 59.776/221.433 s and 47.406/309.390 s.
  Peak process memory reaches approximately 28.9 GB before the next summand.
  This establishes progress beyond the former tensor-allocation failures, not
  completed evaluator generation or numerical acceptance.
- At 11:38:52 UTC, the running direct diagnostic receives the authorized
  100,000,000,000-byte RSS guard without pausing or restarting generation. The
  replacement watcher is armed against the same verified executable/PID before
  the original watcher is suspended; the original supervisor retains ownership
  and its watcher resumes after generation ends. Original receipts are preserved,
  and `guard_handover.json` plus `rss_guard_100gb_receipt.json` record the change.
  The generation driver continues recording the complete-run peak independently.
  Both direct and prepared erased diagnostics now allow 100 GB and 7,200 s;
  executable, cards, assignment policy and physical tolerances are unchanged.

- At 12:02:36 UTC, all 3,050 GL262 summands finish tensor preprocessing. The
  complete Atom scope takes 2,325.890 s: normalization 4.643 s, parsing and
  preparation 339.637 s (including 144.301 s of boundary preparation), execution
  1,785.673 s, and final scalar composition/network cleanup/loop overhead
  195.938 s. Boundary time is a subset of preparation, not an additional term.
  The complete scope has no missing or out-of-order stages in
  `preprocessing_observed_profile.json`. Peak process RSS reaches
  68,568,248,320 bytes before evaluator optimization; full generation is ongoing.
- A separate 20-second, 894-sample contraction profile records inclusive stack
  shares of 34.70% in contraction, 21.68% in root caching and 12.44% in ready
  operation preparation. These overlapping sampled scopes are not additive
  whole-generation timings. Read-only follow-up identifies existing earlier
  graph deletion and alias/function-map APIs as conditional next experiments;
  neither is implemented without evidence from the pending full comparison.

### 2026-09-13 — Separating scalar preparation and evaluator construction

- The full CLI18 direct diagnostic subsequently hits the 100 GB guard at
  12:06:08 UTC, after completing all tensor preprocessing. Its final receipt
  records exit -15 and 3,484.833 s; the guard observes 101,899,558,912 bytes.
  The completed scalar expression occupies 19,982,016,710 Atom bytes. The last
  completed log stage is `evaluator_stack_parse_atoms_done`; no evaluator,
  saved state or production forest export is available. Receipt `complete`
  describes finished supervision, not successful generation. The final profile
  remains ineligible for acceptance. Earlier completed-stage observations and
  failed-run evidence are preserved.
- User steering: "Continue as planned, but maybe in order to separate evaluator
  construction efficiency and expression building, you could lower the evaluator
  building hyperparameters to a minimum, i.e. only 1 horner iteration and 5 CPE
  round only? You can also for now increase the RAM watchdog to 500 GB (if current
  free RAM allows it), but record separately the RAM usage before entering
  evaluator building."
- Future paired diagnostics use one Horner iteration and five CPE rounds.
  The running CLI18 erased diagnostic retains its original five/five settings.
  At 12:16:56 UTC, with 762,059,776,000 bytes available on the host, its guard is
  raised from 100 to 500 GB using the same verified-process handover, without
  pausing generation. The new guard is armed before the old watcher is suspended.
  Before-evaluator RSS and peak RSS will be recorded separately; a larger guard
  and lower diagnostic optimization settings do not establish performance gates.
- Pinned Symbolica source confirms that an empty `replace_multiple` skips
  recursive matching but still copies the complete Atom. The existing generic
  evaluator applies two such passes even when function-map replacements are
  disabled. The candidate removes these two copies, transfers completed scalar
  ownership after optional evaluator variants finish borrowing it, and releases
  unstored source expressions before dual/numeric program construction. Existing
  selector transformations and nonempty two-pass replacements are preserved.
  Additional timing milestones separate expression preparation, each Symbolica
  build and numeric-program conversion. Validation and measured impact are
  pending; the running CLI18 binary remains immutable.
- Format and `cargo check` pass for the ownership change. All 19 selected
  evaluator/symbol tests pass, including nested function-map definitions with
  both replacement settings, both source-retention settings and multiple
  outputs; existing tests retain orientation-selector and hyperdual coverage.
  Clippy passes with the same three existing type-complexity warnings. The
  first check caught a test-only comparison of `DualOrNot` without `PartialEq`;
  the test now uses its existing `unwrap_real` accessor. No numerical expected
  values or tolerances changed. These correctness checks overlap the immutable
  erased diagnostic and are not performance measurements.
- An immutable standalone CLI18 probe is prepared for the existing 3.389 GB
  portable scalar. Its model restriction, physical overrides, LMB, ordered
  parameters and function-map setup are independently checked against the full
  generation path. It replays constituent operations and can stop before
  optimization, after Symbolica construction or after precision conversions.
  It is one term, not the unavailable complete 19.982 GB scalar. RSS/VmHWM and
  exact Atom equality are recorded without expanding the numerator. No probe
  workload has run while the erased generation owns the common measurement
  lock; full-input completion and runtime validation remain required.
- The ownership change is committed and pushed as
  `8ff1bc018af3c0d2f704a90b282acc7baaba46a1`. Immutable CLI19 passes its own
  check/build and has SHA256
  `3e75fa2d9c356af35acd08c44dec5c7e4db9fecc3b1f9f125da61b5fd888d65b`.
  Counter19 passes metadata/link against the exact frozen API/core dependency
  fingerprints; its SHA256 is
  `97f05d15cf636d06ccd6edc948cea81b7fd8fe6b65dd65cdf0f445e4391ce576`.
  Integration check88 and build14 also pass, freezing all four current
  executables under `validation_candidate_8ff1b` before any execution. The
  counter control, small replay control and fourteen integration boundaries
  are queued serially after the running erased diagnostic; none is claimed
  successful from build completion alone.
- New paired CLI19 GL262 diagnostics retain all physical settings and use
  H1/CPE5 with 500 GB guards. The generic final benchmark harness now accepts
  an explicit positive `--horner-iterations`, retaining its historical default
  of five and recording effective settings plus both card hashes. Twelve
  lightweight argument/provenance/resume controls pass. All 27 required final
  cases are prepared with `--horner-iterations 1`; no final generation has run.
- Inspection of the pinned Symbolica direct-translation path separates the
  effect of H1: it skips the additional scheme search but still counts
  indeterminates, constructs a full Hornered expression, processes function-map
  bodies and linearizes them before CSE/CPE. The transformed expression remains
  live until linearization returns. External perf will distinguish these phases;
  verbose optimization logging is avoided because it adds large expression
  scans and temporary operation-count tables. The probe invokes the existing
  evaluator operation and introduces no numerator-expansion command.

### 2026-09-13 — Matched memory stages and typed scalar aliases

- The initially small erased-3D RSS was observed during forest construction,
  before tensor preprocessing. Its immutable CLI18 GL262 run later reached
  213,034,311,680 bytes VmHWM while restoring scalar aliases after the second
  contracted summand. This exceeds direct-4D's completed tensor-preprocessing
  peak of 68,568,248,320 bytes. Neither observation is a completed evaluator
  comparison; direct-4D was stopped by its former 100 GB guard during evaluator
  construction, while erased-3D remains in preprocessing under the approved
  500 GB guard. Erased-3D has 5,274 initial scalar-contraction summands versus
  direct-4D's 3,050. Counts alone do not determine expression size or complexity.
- The completed same-binary forest timings are 900.274 s for direct-4D and
  2,254.471 s for erased-3D (ratio 0.3993). These are diagnostic stage timings,
  not accepted full-generation measurements. Nine nonroot nodes take 879.897 s
  and 2,233.803 s respectively; the shared root takes about 20.4/20.7 s.
- A separate 20 s, 49 Hz sample of erased-3D after term 1 contraction contains
  573 samples. Scalar alias restoration accounts for about 84% inclusive in
  this interval, with roughly 46.6% self time hashing Atom bytes and 23.2% self
  time copying memory. The generic alias map hashes every subtree even though
  this owner's keys are typed scalar-store indices. These sampled percentages
  are interval-local, overlapping where inclusive, and are not whole-generation
  phase fractions. The raw profile and report are retained in the erased case.
- The existing `ScalarAliases::resolve_atom` now uses the existing typed index
  decoder and borrows registered definitions directly. Its exact substitution
  fixed point remains, including forward/nested definitions and handles exposed
  by normalization. A symbol-presence guard avoids the final unchanged Atom
  copy once all handles disappear. No graph-specific rule, new mapper, cache,
  numerator expansion or prepared lookup table is introduced. A dedicated
  evaluator milestone separates alias resolution time and input/output Atom
  bytes from contraction time. Exact old-resolver oracles cover factorized
  spectators, unregistered and malformed handles, nested arguments and stable
  self-references. Validation and isolated timing are pending.
- The first scalar-stage RSS monitor stopped after 13:00:33 UTC. Its original
  files are preserved with hashes. The existing monitor was restarted at
  13:10:14 UTC as a detached process after checking for duplicates; VmHWM still
  captures the intervening peak, but no exact historical stage RSS is invented.
  The guard remained active and generation was never paused.
- Format, Spenso check89 and GammaLoop check91 pass. Four targeted alias/store
  tests and nineteen evaluator/symbol tests pass (nextest36/37), including the
  exact generic-alias oracle. Clippy18 passes with the same three existing
  type-complexity warnings. Check90 first caught an unqualified timer name in
  the new profiling event; it was corrected before compilation. No physics
  expectations or tolerances changed. These checks overlap the immutable
  diagnostic and supply correctness evidence only.
- The CLI20 build wrapper and paired GL262 H1/CPE5 cards are prepared under
  `diagnostic_typed_scalar_aliases`; all physical card contents match CLI19
  after changing only output paths. The supervisor records available RAM after
  acquiring the common measurement lock and refuses launch if it is below
  the requested 500 GB guard. A small generic alias A/B and exact-equality
  probe is prepared separately; its compilation and execution remain pending.

### 2026-09-13 — Closed-tensor bulk sums and separate memory stages

- Milestone 23 is committed and pushed as
  `6efe16006e3ada28dd6362b222df949838f20c87`. Immutable CLI20 has SHA256
  `21eda4e1c2a21715cbff70823a10c2389bc31b217ed5cd5a50447efa22271aad`;
  its matching counter has SHA256
  `6dcac4b16fe9bc525a64f69bc5a93af251257d7692e0ddf3154ca2bd923fa793`.
  Integration check92/build15 pass and all four executables are frozen under
  `validation_candidate_6efe1`. The generic typed-alias probe also passes its
  metadata/link checks. These builds do not establish runtime correctness or
  performance; their workloads have not run yet.
- The immutable erased-3D CLI18 run subsequently reached 326,512,156,672 bytes
  VmHWM after term 4 contraction, still before evaluator construction. Its
  earlier low forest-stage RSS does not describe preprocessing memory. Direct
  4D's observed preprocessing peak remains 68,568,248,320 bytes. The approved
  500 GB guard and original two-hour cap remain in force for erased 3D.
- A separate sample during term 4 contraction contains 790 samples with no lost
  records. About 74.3% inclusive lies in the balanced scalar-sum path, including
  repeated Atom addition and copying; memory copying is about 24.65% self time.
  The raw profile is `tensor_sample_3.perf` in the erased case. This interval
  measures contraction, separately from the previous alias-resolution profile.
  Unresolved kernel samples are not attributed to paging without evidence.
- The existing bulk Atom-sum helper previously admitted literal scalar leaves
  only. Rank-zero tensor results and lazy closed sums instead reached pairwise
  scalar addition. It now uses the same existing scalar conversion as that
  fallback, retaining scales, sparse zeros and typed aliases. Owned results
  move into the existing bulk/streaming sum. Open tensors still reject the
  optimization before any scalar-store insertion. Dispatch thresholds and
  factorization remain unchanged; no new production helper or setting is added.
- The first exact test exposed explicit zero terms retained by the pinned
  streaming implementation as an unnormalized `9*0`. The sum boundary now
  omits explicit zero inputs and represents an empty streamed result by exact
  scalar zero. The test's equality expectation is unchanged; added cases cover
  all-zero and cancelling sums on both sides of the streaming threshold.
  The initial format91/check93 passed and nextest38 ran eleven passing tests
  plus this failure. The corrected implementation is being checked again.
- An isolated public-network A/B probe is prepared for rank-zero tensor sums.
  CLI20 dependencies were copied and hash-verified before later compilation.
  Its first metadata check caught an artifact import path; the corrected probe
  passes metadata/link with the same frozen dependencies. Execution has not
  started. It retains powers of sums and checks exact equality against bulk
  addition of the original factorized inputs outside the execution timer.
- Profiling receipts now validate and report typed-alias resolution separately
  as a subset of scalar finalization, avoiding double counting. Seventy-one
  parser controls pass; old logs retain an unobserved alias interval rather
  than claiming zero cost. Reprocessing the failed CLI18 direct run preserves
  its previous phase totals and failed eligibility.
- The final physical benchmark harness now supports the explicit 500 GB RSS
  guard and checks available memory after acquiring the measurement lock.
  A guard-triggered attempt is ineligible even when its child handles SIGTERM
  and exits zero. Twenty-nine isolated guard/status controls, twelve existing
  Horner-setting controls and fifteen receipt controls pass. Prepared physical
  cards remain byte-identical. All final generation/runtime gates remain open.
- Format92/check94 pass for the corrected sum boundary. Nextest40 passes all
  twelve selected network tests; nextest41 passes all nineteen selected
  evaluator/symbol tests. Clippy20 passes with the same three existing
  type-complexity warnings. Exact zero, scale, alias and open-tensor assertions
  remain intact. The original plan prefix is unchanged. These checks overlap
  the old diagnostic and are correctness evidence only.

### 2026-09-13 — Capped erased baseline and isolated alias/sum controls

- Milestone 24 is committed and pushed as
  `a8b4c9d96c3994f77bba80826e16715d53ddbfc9`, with ValentinHirschi as author
  and committer. Immutable CLI21 has SHA256
  `91c9e3389d9b7c3cd614890d24d46422cb45c4339efa3333db8858f958733fe2`;
  counter21 has SHA256
  `10ad4be73def5d438ccbde7533478172074d82d53307c4615f0a6615a8dc68cd`.
  Check and build pass in 10.233 s and 203.393 s. The counter's known-GL00
  control passes with exact count equality and an unchanged saved state.
  Paired diagnostics and all 27 final H1/CPE5 cases have finalized identities;
  preparation alone is not generation or performance evidence.
- The immutable CLI18 erased-3D GL262 run stopped at its original time limit:
  exit -15, outer wall 7,215.814457280 s, termination requested after
  7,200.069243180 s. Its observed peak was 375,509,811,200 bytes. The approved
  500 GB RSS guard did not trigger. All ten forest nodes completed in
  2,254.471267532 s; only seven of 5,274 tensor summands reached execution
  completion. No pre-evaluator boundary, evaluator, saved state, production
  forest export, counter output or runtime exists for this failed attempt.
- Through term 6 execution, the observed preprocessing prefix is
  4,857.316288205 s: normalization 7.799106124 s; parsing and preparation
  297.422270102 s; contraction 1,911.610710178 s; preceding scalar finalization
  and loop overhead 2,640.484201801 s. Boundary preparation is a 121.758081121 s
  subset of parsing. Legacy alias-resolution time cannot be isolated from the
  residual and remains explicitly unobserved. Prior receipts are preserved;
  separate `final_failed_audit.json` and `final_failed_profile.json` record the
  completed observations and failed eligibility.
- The queued counter19 and H1/CPE5 scalar-replay controls pass. The latter
  reproduces exact value 36 and all numeric-domain conversions. All fourteen
  queued integration boundaries pass on their frozen milestone-22 executables.
  This is historical validation of that revision, not a substitute for current
  integration execution.
- The generic alias probe passes exact same-frame equality for nested and
  alias-free inputs at both tested sizes. At scale 65,536, three fresh runs per
  resolver give nested medians 2.456507315 s (typed) and 2.815747797 s (old),
  ratio 0.8724, with the same 75,487,207-byte result. Alias-free medians are
  0.071083854 s and 0.220198779 s. These are isolated synthetic diagnostics,
  not physical generation or numerical-runtime ratios.
- Nested process VmHWM medians are 388,702,208 and 380,751,872 bytes, a 2.09%
  increase for the typed resolver. Periodic current-RSS maxima suggested a
  much larger gap because they missed the old resolver's brief peak; use the
  kernel high-water marks for this comparison. Both owners use the same
  rebuilding/normalizing replacement routine. Temporary Atom buffers and
  allocator retention can explain modest endpoint differences, but allocator
  traffic was not measured. The probe has matched isolated lifetimes; its
  two-result exact check runs separately and is excluded from timing ratios.
- All 26 public-network sum-probe processes pass exact equality. Small inputs
  expose streaming overhead: at 32 by 1,024 terms the CLI20/CLI21 medians are
  6.967245/12.255105 ms. At 128 by 4,096 terms (26,181,642 output bytes), they
  are 135.520324/147.611895 ms. At 256 by 4,096 terms (52,363,274 bytes), they
  are 289.677145/284.544955 ms, with overlapping ranges. These short, generic
  measurements neither prove a physical speedup nor justify selecting a path
  by graph identity. Their memory samples cover whole processes, not the timed
  execution alone; tiny cases finish between polls.
- Source inspection confirms that in-memory bulk addition borrows normalized
  additive inputs through merge cursors, whereas streaming owns and sorts
  individual terms before rebuilding the result. All tested sizes fit below
  the existing 4 GiB stream buffer. A total-byte criterion may avoid needless
  streaming for small sums, but it has not been implemented or validated. The
  untimed oracle's residual interval includes extraction, bulk addition,
  destruction, comparison and logging; it is not a pure native-addition timer.
- Integration check95/build16 pass in 10.017374494 s and 203.933597533 s;
  all four current executables are frozen in `validation_candidate_a8b4c`.
  Their full correctness execution remains pending. A broader current unit
  suite is now running separately from acceptance measurements.
- The immutable CLI21 direct-4D GL262 diagnostic started at 14:18:09 UTC with
  H1/CPE5, one worker and the 500 GB guard. The independent stage monitor was
  armed after verifying binary, card, process start and guard identities.
  Available memory before monitor launch was 829,036,617,728 bytes. It records
  the tensor/evaluator boundary and subsequent preparation, Symbolica-build
  and numeric-program stages with explicit observation lag. No completed
  generation or final gate is claimed while this run is active.
- Format93/check96 pass. Nextest42 passes 597 selected Linnet, Spenso, UV,
  CFF, allocation and evaluator tests. A disjoint selection of the remaining
  signed/affine source-mapping cases passes all 32 tests in nextest43: 629
  current tests in total, including the fixed-source 62-row oracle and exact
  reconstruction boundaries. No assertions, tolerances or physics changed.
  These correctness checks overlap the diagnostic and are not timing runs.
- The existing `save uv-forest` command requires a generated integrand even
  for topology-only export. Its internal constructor/serializer already covers
  both forest orchestrators and needs no new topology algorithm. Current trace
  evidence establishes GL262's ten paths/nine arcs and independent four-region
  classification, but does not substitute for the still-pending production DOT.
  The existing topology verifier can consume that export after generation.
