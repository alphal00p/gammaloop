= Symbolica expansion performance reproducer

This standalone package isolates expansion and polynomial-to-expression
conversion. It depends only on Symbolica main `06906976`
and the Rust standard library. No Idenso rewrite rules, Spenso slot parsing,
tensor networks or pattern matching execute in the timed operation.

The input is the memoized, factored Clifford pairing polynomial with `Tr(1)=4`.
All free-index metrics are algebraically independent. The `tensor` payload is
`g(mink(D,mu_i),mink(D,mu_j))`; `shallow` removes the representation wrapper;
`variables` replaces each metric by a scalar variable. These deliberately
different payloads expose the cost of carrying tensor syntax through expansion.

== Run

From this directory:

```sh
cargo build --release --locked
./target/release/expansion 12 tensor expand
./target/release/expansion 12 tensor poly
./target/release/expansion 14 variables expand
```

Use even lengths from 2 to 14; compare at least 8, 10, 12 and 14. On Linux, pin
each process with `taskset -c N` and repeat in fresh processes, alternating
methods. The driver uses optimization level 2 and dependencies level 3, matching
the measured workspace configuration. GMP and MPFR are required. `Cargo.lock`
pins a separate standalone dependency tree; the historical matrix retains the
identities of its original frozen workspace libraries in the benchmark record.

Each invocation emits five JSON records. Input construction, an untimed reference
expansion, and exact output comparison are excluded from `expansion_ns`.
Destruction is recorded separately in `drop_ns`. The output must equal the
ordinary expansion and contain $(n-1)!!$ terms. This checks the two expansion APIs
against each other; independent Clifford and FORM validation belongs to the
full trace benchmark. It is not a proof based solely on counting terms.

== Conversion only

Use `polynomial_emission` to isolate the current sparse trace emitter's slow
Symbolica operation. The existing `expansion ... poly` command times
`expand_via_poly`, which also constructs the polynomial. The focused binary
constructs a signed perfect-matching polynomial with integer coefficients and
byte exponents before timing, then calls only `to_expression()` in the
non-inlined `profile_run` function. It uses the same polynomial domain and
conversion API as Idenso's sparse emitter, with no Idenso or Spenso dependency.

From this directory:

```sh
cargo build --release --locked --bin polynomial_emission
./target/release/polynomial_emission 12 tensor 5
./target/release/polynomial_emission 14 tensor 5
./target/release/polynomial_emission 12 variables 5
```

Arguments are even length 2–14, payload (`tensor`, `shallow`, `variables`), and
sample count; defaults are `12 tensor 5`. An optional fourth argument,
`roundtrip`, forces inverse-conversion validation for larger inputs.
Twelve produces 10,395 monomials;
fourteen produces 135,135. The tensor-shaped payload is
`g(mink(4,idx(mu_i,0)),mink(4,idx(mu_j,0)))`. These are inert Symbolica functions,
representative of the current nested slot spelling; they have no Spenso tags or
custom tensor callbacks. The shallow and variable modes retain the polynomial
while reducing its atoms' payload.

Each JSON line reports `emission_ns` and `drop_ns`. Polynomial construction,
warmup and validation are outside both clocks. `emission_and_drop_ns` sums the
two measured intervals; validation lies between them, so it is not a contiguous
full-lifecycle measurement. `exact_equal` compares a timed result with the warmup
output from the same library. It is not an independent correctness certificate
for a modified Symbolica implementation.

For lengths up to ten, an exact inverse polynomial roundtrip also compares every
coefficient with the input using the same variable ordering. Larger runs skip
that separate expensive operation and report `roundtrip_checked: false` and
`roundtrip_equal: null`. Pass `roundtrip` to enable it explicitly. The initial
unconditional inverse check passed through twelve but timed out at fourteen
after 120 seconds before any timing record; it is not part of the conversion
measurement being optimized here.

For an instruction profile of conversion alone:

```sh
valgrind --tool=callgrind --collect-atstart=no \
  --toggle-collect='*polynomial_emission::profile_run*' \
  ./target/release/polynomial_emission 12 tensor 1
```

To work on a local Symbolica checkout, replace this standalone manifest's
`git` and `rev` dependency keys with `path = "/path/to/symbolica"`, retaining
the listed features. This package has its own workspace and lockfile; the main
FeynKit workspace dependency remains separate. Build both revisions with the
same release settings and compare several sizes and payloads. The earlier
isolated patches in `patches/` are not applied to the default pinned build.

The unpatched focused binary was measured on 2026-09-25, on CPU 10 of the same
shared host, with three processes per size/payload and five samples per process.
Median conversion-only times are:

#table(
  columns: 5,
  [Length], [Terms], [Tensor, ms], [Shallow, ms], [Variables, ms],
  [8], [105], [0.0670], [0.0639], [0.0610],
  [10], [945], [0.7389], [0.6585], [0.6412],
  [12], [10,395], [10.7384], [8.8399], [8.3021],
  [14], [135,135], [235.6720], [142.2364], [131.1368],
)

All 36 processes and 180 timed-output comparisons pass. The small-input inverse
roundtrips pass; larger inputs use the documented warmup comparison. The
Callgrind toggle command was also checked on length eight and collected only
the profiled conversion. Build, formatting and scoped Clippy pass. Raw samples,
source identities, the setup timeout and validation limits are retained under
`focused_polynomial_emission_mre` in `primitive_measurements.json`. These figures
are a fresh MRE baseline, not an attributed improvement over historical tables.

== Measured bottleneck

The companion #link("../../notebooks/tensor_contraction_parity.json")[benchmark
record], under `expansion_primitive_followup`, retains the original matrix driver,
its exact build command, raw samples and profile. That larger driver additionally
checks inverse substitution for every payload, including contracted words.
Its timings are separate from this minimal reproducer, which adds an untimed
reference expansion and may therefore have a different allocator state.

On 2026-09-24, three processes with five samples each on CPU 9 of a shared
EPYC 9754 host measured the following expansion construction times:

#table(
  columns: 4,
  [Free gammas], [Tensor, ms], [Shallow, ms], [Variables, ms],
  [12], [37.49], [35.41], [33.77],
  [14], [756.93], [635.30], [580.23],
)

At length 14, `expand_via_poly` with byte exponents reduces the tensor case to
679.91 ms. Scalar variables reach 414.22 ms, but substitution back to metrics
takes another 353.61 ms, before counting the forward substitution. These results
do not establish an end-to-end gain from replacing tensor atoms with variables.
All 1,080 inverse-substitution checks in the six-case matrix agree exactly.

The length-12 ordinary expansion records 384.8 million Callgrind instructions
with tensor payloads and 340.1 million with scalar variables, an 11.6% reduction.
Exclusive instruction counts attribute about 19.8% to normalization, 11.7% to
product/sum buffer extension, 8.3% to factor comparison, 8.3% to copying and 4.5%
to byte comparison. These are instruction shares, not wall-time percentages.
The profile excludes trace construction and contains no pattern matching.

For context, the separately measured production length-14 pipeline takes
789.14 ms including trace construction, expansion and destruction, against
FORM 5.0.0's 53.67 ms internal `tracen` plus sort CPU time. FORM and Idenso
return the same 135,135-term polynomial. The clocks and experiment lifecycles
differ; neither the factored-output timing nor an expansion-only prototype
establishes parity with FORM.

== Isolated primitive experiments

The companion #link("primitive_measurements.json")[primitive measurement record]
contains paired samples, executable driver sources, build identities and exact
output checks. The patches in `patches/` target Symbolica `06906976`; they are
experiments and are not applied to this repository's dependency.

The original emission patch changed the general fallback's coefficient order.
A later review found two counterexamples: a one-bit rounding difference for
53-bit coefficients with compound variable mappings, and reversed normalization
callback order for raw variable/coefficient expressions. The current patch
preserves the original fallback order and converts each coefficient once.
The historical primitive timings below describe the original patch; the
corrected full-host comparison is recorded separately in the gamma guide.
The two permanent regression tests pass on the pinned baseline and corrected
patch, and fail on the original patch. Run them with:

```sh
cargo test --release --locked --bin polynomial_emission
```

`poly-emission.patch` emits normalized square-free monomials directly when their
variables are distinct normalized symbols or function calls. Other variable
forms, aliases and powers retain the existing normalization. The separate
`poly-presence.patch` checks whether a variable occurs, rather than whether its
maximum degree is nonzero. This also fixes the isolated signed-exponent boundary
case `1+x^-1`, whose maximum degree is zero although `x` occurs.

Three rotated process rounds with seven samples each measure the following
conversion times for already constructed polynomials:

#table(
  columns: 4,
  [Free gammas], [Original, ms], [Direct emission, ms], [Plus presence check, ms],
  [8], [0.0625], [0.0349], [0.0331],
  [10], [0.7394], [0.3333], [0.3133],
  [12], [8.8232], [4.3299], [3.6278],
  [14], [198.9511], [122.0897], [104.1423],
)

All 504 measured outputs agree exactly. Separate boundary checks exercise
aliasing, compound variables, powers and the signed-exponent case. This is not
the full upstream test suite. A power-heavy control regresses by about 4%.

The effect on complete expansion is smaller. A separate paired experiment
loads the exact factored production outputs and includes expansion plus result
destruction: length-14 `expand_via_poly` improves from 700.8 to 605.7 ms, while
ordinary `expand` stays at about 670 ms. Length-12 polynomial expansion improves
47.6 to 42.0 ms but remains slower than ordinary expansion at 34.6 ms. All 126
API comparisons and the fingerprints of all 630 timing records agree exactly.
These isolated builds share dependency identities within each comparison;
their timings must not be subtracted from the earlier production lifecycle.

The addition experiment uses flat sorting when every input to `add_many` is a
single term, retaining original input order when equal terms merge. Its first
versions improve addition of the 135,135 final length-14 monomials from 152.3
to 37.1 ms, but regress existing-sum controls by approximately 4–12%. Updating
the heap's front cursor in place makes the largest all-sum control slower still;
that variant is rejected.

The subsequent variant keeps the full cursors in a vector and moves only their
comparison keys and positions through the heap. Three rotated process rounds
on CPU 6 measure:

#table(
  columns: 3,
  [Addition input], [Original, µs], [Compact heap and singleton sort, µs],
  [Free-12 final monomials], [4,489.63], [1,078.19],
  [Free-14 final monomials], [158,176.08], [35,664.37],
  [64 existing sums], [52.30], [51.86],
  [512 existing sums], [595.34], [572.46],
  [512 mixed sums and single terms], [139.96], [112.93],
)

None of the 22 measured case medians regresses against the original addition.
This is an isolated addition result, not a complete trace or expansion result.
The earlier 63 boundary checks compare textual output; the timing matrix has
198 untimed oracle checks and 1,386 timing records. A stronger certificate
compares raw Atom bytes for 1,017 case/permutation records, including 939 with
actual sum cursors. All input streams, first outputs and rerun outputs match
between versions. It covers rounded 11-, 53- and 80-bit coefficients, zero gaps,
nonfinite values, finite fields and rational-polynomial coefficients. Some
rounded and finite-field zeros change representation on rerun in both versions;
the certificate preserves that baseline behavior rather than assuming bytewise
idempotence.
All seven upstream normalization unit tests also pass using the same isolated
wrapper lock and dependency features. This is not the full Symbolica test suite.

A further trial merges already normalized product factors during ordinary
expansion. It passes 156 binary baseline comparisons but regresses length 12
from 35.0 to 42.3 ms and length 14 from 671.6 to 759.2 ms. The length-12
Callgrind total grows from 385.4 to 406.1 million instructions: normalization
falls from 72.9 to 8.8 million, but collecting factor lists and eagerly
multiplying unit coefficients costs more than it saves. This trial is rejected;
its complete sources, raw samples and profiles remain in the measurement record.

== Direct expanded trace construction

A separate constructor evaluates a copy of the Clifford word identities directly
into an integer-coefficient polynomial, avoiding the factored-expression-to-polynomial
conversion. It starts from an already classified word and includes the metric
table, dimension parsing, leaf enumeration, polynomial construction, Atom
emission and destruction. It excludes public gamma-AST recognition and cleanup.
The #link("../../notebooks/tensor_contraction_parity.json")[contraction record]
retains the sources, dependency fingerprints and all samples under
`direct_expanded_constructor_followup`.

For free length 14, the constructor takes 256.79 ms with the frozen production
library. In a separate matched dependency family, the original library takes
248.23 ms and the polynomial-emission patches reduce this to 157.61 ms.
Adding the earlier addition patch gives 159.82 ms: polynomial emission builds
an Add directly and does not use the optimized `add_many` path. All 432 untimed
mode comparisons pass across the twelve-case, four-library matrix.

Expanding every contracted word is not a suitable default. With actual Spenso
metric registration, the branching-eight result grows from 1,852 to 3,398 Atom
bytes and its rerun slows from 29.5 to 52.6 µs. The order-sensitive twelve-gamma
result instead shrinks from 133,818 to 52,401 bytes and its rerun improves from
1.927 to 0.764 ms. All thirty downstream comparison processes agree algebraically.
The constructor driver registers a symmetric metric only; the downstream check
uses Spenso's full attributes and normalization callback. A public integration
must therefore be measured with the real registration and retain scalar
spectators in their original factored form.

A subsequent copied-library trial shares the production word evaluator between
factored and sparse output algebras. The sparse representation stores and clones
lists of integer monomials at intermediate subwords, unlike the earlier leaf
emitter. Its public order-sensitive-twelve lifecycle, including expansion and
destruction, improves from 5.995 to 1.472 ms; the separate FORM reference is
0.733 ms. However, the default first call slows from 0.974 to 1.405 ms.
The rerun improves from 1.855 to 0.725 ms.
This trial is rejected; the retained implementation is unchanged.
Both copied libraries use matching direct-rustc flags and frozen production
dependencies. Cargo's incremental compilation and codegen partitioning are not
reproduced, so this is a matched public-API experiment, not a replacement for the
retained Cargo-build checkpoint.

The public driver uses bare traces, which return directly from the terminal
evaluator. Thus post-trace cleanup cannot explain the 1.405 ms first call.
A subsequent profile identifies intermediate sparse storage as a material cost:
the order-sensitive case increases from 9.34 to 14.06 million instructions and
from 3,232 to 11,609 allocations. Memo-result cloning costs 2.31 million
instructions, and sparse finalization costs 6.77 million. Earlier wall-time
regressions in unchanged controls do not consistently reproduce. Free lengths
twelve and fourteen retain identical allocation counts, with instruction counts
within 1%; duplicate coefficient construction in the interior-pair control adds
only 0.03% instructions. Its causal profile is retained under
`direct_expanded_public_profile` in the contraction record.

The representation-only follow-up under `direct_expanded_shared_nodes_trial`
stores handles in the existing factored-trace nodes and emits polynomial leaves
once. The matched order-sensitive-twelve full lifecycle takes 6.276 ms with
factored output, 1.570 ms with sparse lists and 1.099 ms with shared nodes.
Instructions decrease from 14.06 to 9.56 million and allocations from 11,609 to
3,026. All 52 trace cases and 87 boundary/callback files exactly match the
validated sparse-list trial. This remains a scratch-library experiment; it does
not change production or repeat the FORM measurement.

Before timing, the trial also fixes two confirmed boundaries: mixed
free slots and compact vectors can invoke callbacks producing non-polynomial
variable forms, and replacing binary addition with variadic addition can change
factorization. The final trial conservatively excludes mixed residual arguments
and preserves the original binary operations. Its source, controls and timings
are retained in the contraction record.

== Inserting an atomic factor into a normalized product

A subsequent isolated Symbolica trial avoids renormalizing a product when
expansion multiplies it by one distinct variable or function, optionally with an
exact coefficient. It checks factor ordering and coefficient domains, inserts
the factor into the existing byte representation, and leaves merging, powers
and unsupported coefficient cases to normal multiplication. Both sides of this
comparison already include the compact addition and polynomial-emission trials.

Three alternating process pairs give ordinary expansion times of 35.16 to
25.51 ms for the twelve-gamma output and 667.28 to 474.83 ms for fourteen gammas.
The twelve-gamma instruction count falls from 386.8 to 265.6 million. All 156
archived expansion boundary comparisons pass, including changed flags and
reruns; seven real fixtures agree across ordinary and polynomial routes.
Input construction and parsing are excluded, with destruction recorded
separately. These remain shared-host measurements of an isolated dependency.

The first version allocates scratch storage even when every summand is a
variable, so the insertion routine immediately declines it. Simple
function-times-sum controls regress by about 9–11%, with 268–369 extra
instructions. Allocating scratch only when a normalized product needs it
reduces those excess instructions to 13–89. The final isolated version passes
225 exact expansion comparisons, including changed flags and reruns, plus the
existing header-width test. Its three-process median of process medians is:

#table(
  columns: (auto, auto, auto),
  [Fixture], [Baseline expansion], [Insertion with lazy scratch],
  [Free eight], [0.163 ms], [0.115 ms],
  [Free twelve], [34.51 ms], [25.54 ms],
  [Free fourteen], [647.87 ms], [453.52 ms],
  [Interior pair, five], [17.36 ms], [11.76 ms],
  [Power control, twenty], [3.23 ms], [3.22 ms],
)

The twelve-gamma instruction count is 268.9 million, compared with the same
386.8-million baseline. Small wall controls remain mixed: for example, the
positive function-times-sum case takes 1.000 to 1.025 µs and the rational-prefix
case 1.247 to 1.335 µs. Thus this is not a universal no-regression result, and
the dependency remains unchanged in production. The `full_fn_cmp` feature
configuration has not been built. The source patch is
`patches/atomic-product-insertion.patch`; validation, raw timings, source and
identities are stored under `atomic_product_insertion_trial` in
`primitive_measurements.json`, with the final results in
`lazy_scratch_validated_trial`.

A historical explicit trace-expansion API trial kept default factorization
unchanged and used the shared-node emitter when expansion was requested. Its
standalone order-sensitive twelve-gamma result improved from 5.96 to 0.83 ms,
and fourteen free gammas from 805 to 356 ms. That prototype was not integrated:
with the untouched scalar spectator `(x+y)^8`, the free-twelve comparison
regressed from 68.49 to 99.84 ms.
The reference completes the factored expression before expanding its extracted
trace body; the proposed API expands inside the trace rewrite, exposing the
larger output to outer cleanup. Profiles show essentially unchanged trace-kernel
and expansion costs, with the increase in subsequent scans and cleanup.
The candidate, validation and these explicit timing boundaries are archived
under `explicit_expanded_trace_api_trial` in the contraction parity record.

The later integrated opt-in API addresses that outer-cleanup boundary. See the
#link("../../../docs/products/idenso/content/gamma-simplification.typ")[gamma
simplification documentation] and `expanded_trace_live_integration` in the
contraction parity record for its implementation and validation. The historical
measurements above do not describe the current API's performance.

== Corrected conversion in the full host

The corrected emission patch and unchanged presence patch were subsequently
measured in a copied full Idenso/Spynso host, preserving all other package
versions, features and dependency relationships. Against the saved host,
complete tracen calls improve 18.508 to 13.088 ms at length 12 and 363.202 to
265.865 ms at length 14. Dispatch, result wrapping and recombination with a
factored scalar spectator are included; source construction and checks are not.
These figures replace neither the earlier isolated primitive clocks nor the
unpatched production dependency. Trace4 and ladder controls show no reliable gain;
the raw ladder median is 0.7% slower.

Both new permanent regressions pass the pinned baseline and corrected patch,
and fail the original patch. The final host also passes the ten earlier boundary
tests, the Laurent case, 207 exact Python behavior records and 117 HEP component
checks. Scoped Clippy, formatting and the standalone Cargo tests pass.
The corrected conversion still occupies 76.61% of sampled tracen cycles,
including 38.31% in final Atom normalization; coefficient-list construction
accounts for 14.33% separately. These are sampled cycle shares, not wall-time
phases. The `corrected_emission_full_host_followup` entry in
#link("primitive_measurements.json")[the primitive record] retains the separate
build, profiles, corrected patch, original failing examples and full timings.
