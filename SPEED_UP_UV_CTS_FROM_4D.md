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
