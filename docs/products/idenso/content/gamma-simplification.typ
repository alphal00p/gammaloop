#import "../../shared.typ": source-link, product-link

#let gamma-simplification = [
= Gamma simplification notebook

Explore the current Dirac simplifier through executable identities, configurable chain
ordering, and examples of expressions that deliberately remain explicit. Each identity
card checks its expected result and verifies that a second pass leaves it unchanged.

#context if target() == "html" {
  html.elem("div", attrs: (
    class: "live-notebook",
    "data-notebook": "gamma_simplification",
    "aria-label": "Idenso gamma simplification notebook",
  ))[
    #html.elem("p", attrs: (class: "live-notebook-fallback"))[
      This notebook requires JavaScript and the documentation build's notebook assets.
      Open the #source-link("examples/notebooks/gamma_simplification.py", label: "Marimo source")
      to run it locally.
    ]
  ]
} else {
  [Open the #source-link("examples/notebooks/gamma_simplification.py", label: "Marimo source")
  to run the interactive notebook.]
}

== Capabilities and conventions

The examples cover contracted open chains, four-dimensional Chisholm reduction, symbolic
Lorentz dimension with independent spinor trace normalization, slashed momenta, ordinary
traces, canonical ordering, gamma-five, chiral projectors, the optional three-gamma epsilon
expansion, and charge conjugation. A small adjustable trace illustrates the growth in the
number of metric pairings.

The default repeated-pair ordering avoids unnecessary open-chain expansion. Canonical
ordering can expose cancellations between different gamma orders. Gamma-five rules remain
strictly four-dimensional: the simplifier does not choose a dimensional-regularization
prescription. Mixed dimensions and unsupported transposed words can remain unevaluated.
The notebook checks these boundaries as well as successful reductions.

See #product-link("idenso", page: "reference/form-color-dirac/", label: "Shipped color and Dirac rules")
for the conventions and #product-link("idenso", page: "guides/showcase/", label: "the tensor display showcase")
for additional notation controls.

== Schoonschip performance and FORM comparison

The performance section follows equivalent paired and alternating slash traces through
three routes: evaluate a free-index trace then contract momenta; contract the indexed
input with `schoonschip_net()` before taking the trace; or start from the compact slash
expression. Each route must match an independent scalar identity. The displayed speedup
compares the first two complete pipelines, including early conversion cost. The compact
input timing is reported separately. Changing the display layout alone does not change
the algebra or its cost.

Measurements use three timed runs after a warm-up, rotate execution order, and exclude
construction, assertions, and rendering. They characterize the installed build and machine;
they are not a release benchmark of the engines.

Closed traces now reuse adjacent contractions and the four-dimensional Chisholm
identities from open chains. A macro dispatches lengths 1–14 to generated short-trace
kernels, with separate ordinary and gamma-five tables and immediate zero for odd words.
The integer recipes follow the same three-gamma epsilon identity as open chains and
are initialized lazily. Standalone short traces with distinct explicit indices can
return their terminal metric products directly; surrounding contractions retain the
full tensor simplification pass.

Ten, twelve and fourteen free indices now give 693, 4,383 and 26,931 ordinary trace
terms, matching FORM 5.0.0 `trace4`; generic-dimensional pairing gives 945, 10,395 and
135,135. FORM uses a general four-dimensional reduction internally: the arity-specific
caches are our implementation choice. This agreement does not establish full FORM
parity, including its cross-trace Chisholm and symmetrization options. See the
#link("https://github.com/form-dev/form/blob/master/doc/manual/gamma.tex")[FORM Dirac-algebra manual]
and #link("https://github.com/form-dev/form/blob/v5.0.0/sources/opera.c")[FORM 5.0.0 source].

When FORM is installed locally, the notebook runs the selected compact scalar case with
both `trace4` and `tracen`, plus a fourteen-index trace that better exposes their cost
and output-size differences. It checks the compact results and the free-index results
contracted with paired momenta against exact scalar identities. Generated programs and
FORM output remain visible in the notebook. An optional table times the same fourteen-index
trace in the installed Idenso host; it is disabled by default because reconstructing
Python tensor metadata for 26,931 terms can take minutes. First and warm calls distinguish
lazy initialization when the table has not already been used. Native FORM process timings
include startup, parsing, sorting and verification; Idenso's timings are in process.
This establishes selected identities and performance examples, not full FORM correctness
or performance parity. Browser exports skip native execution.

The native Rust example measures all even lengths through fourteen, for distinct free
indices and paired or alternating momenta. Both engines must pass exact scalar checks;
FORM programs and logs are retained at the printed path. Its CSV reports first-call
and warm Idenso timings alongside FORM `trace4` and `tracen` process timings. Use the
same optimization profile when comparing revisions.

// docs-example: syntax
```sh
cargo run -p idenso --profile dev-optim --example trace_form -- /path/to/form
```

== Recorded native measurement

On 2026-09-23, with Rust 1.98.1, Symbolica 3.0.0, `dev-optim` and FORM 5.0.0 on
a shared AMD EPYC 9754 host, the example measured the following medians in milliseconds:

#table(
  columns: 5,
  [Case], [Old Idenso], [New Idenso], [FORM `trace4`], [FORM `tracen`],
  [10 free indices], [168.05], [1.12], [10.81], [10.63],
  [12 free indices], [Not measured], [8.00], [16.54], [22.57],
  [14 free indices], [Not measured], [63.61], [54.42], [190.45],
  [14 paired slashes], [Not measured], [0.75], [9.12], [8.74],
  [14 alternating slashes], [Not measured], [9.52], [9.58], [9.50],
)

The ten-index case improves by about 150 times in this controlled before/after
comparison. The first fourteen-index call takes 73.12 ms, including initialization.
FORM's fourteen-index pipeline is about 3.5 times faster with `trace4` than `tracen`.
The roughly 9 ms FORM process overhead hides its advantage on smaller traces,
creating an apparent crossover near length fourteen. These timings do not establish
engine parity. Full symbolic scalar checks passed for
every reported input. The #source-link("examples/notebooks/gamma_trace_measurements.json", label: "measurement record")
contains every length, the baseline provenance, timing boundaries and source hashes.

The separately measured Python notebook's default eight-factor paired example gives
590.16 ms for tracing first, 18.88 ms for contracting first, and 4.66 ms for compact
input: about 31 times faster with early contraction. This unoptimized Python-host
measurement is not compared directly with the optimized native Rust numbers above.

== Amortized FORM timing and output construction

A follow-up benchmark on the same host amortizes FORM startup over six batches of
independent traces per process. Complete process wall times include parsing, tracing,
sorting and disposal, but exclude scalar verification. The median of three processes,
divided by the trace count, compares with five warm in-process Idenso wall-time samples:

#table(
  columns: 4,
  [Free indices], [Idenso (ms)], [FORM (ms)], [Idenso / FORM],
  [8], [0.150], [0.055], [2.75×],
  [10], [1.100], [0.386], [2.85×],
  [12], [8.094], [2.887], [2.80×],
  [14], [61.770], [22.198], [2.78×],
)

FORM is already faster at every measured length. Shared-host variability matters:
its fourteen-index samples range from 19.38 to 25.47 ms per trace. The batch approach
also holds more simultaneous expressions in FORM than the Idenso call loop. These
results characterize this workload, rather than general engine parity.

An instrumented copy of the fourteen-index kernel attributes about 36.0 ms to
materializing and normalizing products, 23.5 ms to canonicalizing the sum, and
0.15 ms to constructing metric atoms. Instrumentation changes the loop slightly,
so these diagnostic phase medians do not sum exactly to the API median. Warm calls
exclude integer-recipe generation, and standalone free traces already bypass the
later fixed-point passes. Output construction is the next optimization target.
FORM's `Trace4Gen` writes packed metric-index records into reusable scratch storage;
our kernel constructs general Symbolica products and sums. This supports the
implementation explanation, without isolating allocation as the sole cause.

The #source-link("examples/notebooks/gamma_trace_scaling.json", label: "scaling measurement record")
contains raw samples, source hashes and timing boundaries. Reproduce the complete
API and FORM measurements with the native example; it saves generated FORM programs
and logs and prints every sample as CSV. The historical phase instrumentation is
separate from this reproducible API comparison.

// docs-example: syntax
```sh
cargo run -p idenso --profile dev-optim --example trace_scaling -- /path/to/form
```

== Replacement passes and terminal expressions

The ordinary standalone trace bypasses the fixed-point pipeline. Other inputs can
spend much more time traversing output that is already simplified. In a diagnostic
copy before the metric-first contractor below, a ten-gamma trace takes about 1.10 ms alone
but 150 ms with a scalar prefactor, which prevents that narrow shortcut. Most of
the additional time is in Schoonschip, chain collection and dot normalization.

A twelve-gamma axial trace takes about 2.56 s, with 2.36 s (92 percent) in epsilon
simplification. In each of two outer passes, the epsilon pass visits 44,248 nodes,
runs `term.schoonschip()` at each node, and applies zero epsilon identities.
Every nested Schoonschip call also leaves its argument unchanged. Those nested
calls consume about 1.96 s altogether; the actual epsilon-rule checks take only
about 7.4 ms. The second outer pass runs because the first replaced the trace with
a polynomial; equality is tested only after every stage runs again.

Bypassing epsilon processing in an isolated diagnostic reduces that case to about
0.21 s with exactly the production result. This does not justify skipping epsilon
contractions on general inputs. The first gamma rewrite alone takes 2.33 ms
and already equals the final production polynomial in this case. A single
fourteen-index scalar-prefactor probe takes 8.92 s, compared with a three-sample
median of 59 ms for the bare trace. Repeated processing of terminal expressions
is the first optimization target. A fourteen-factor alternating slash
trace is different again: eleven outer passes take 9.66 ms, including 3.76 ms in
Schoonschip, 2.02 ms in collection and 2.50 ms in gamma rewriting.

For the standalone ordinary kernel, borrowed factor views remove avoidable copies:
a separate five-sample experiment changes length fourteen from 54.3 to 48.3 ms,
about eleven percent. Symbolica still normalizes each product and merges the sum.
Sampling highlights product normalization, byte comparisons and heap operations.
The #source-link("examples/notebooks/gamma_trace_profile.json", label: "detailed profile record")
preserves raw timings, source hashes, diagnostic variants and their limits.
Production algebra was unchanged during this investigation.

== Metric-first contraction

The pattern Schoonschip pass now selects one metric with
`g(a_,b_)*remainder___`, then walks the captured remainder to find a compatible
explicit tensor slot. It substitutes that slot with the opposite metric endpoint
and removes the metric. The matcher captures the entire product at the current
tree level; it does not search for a metric and its tensor partner or partition the
remainder into subsets. A failed candidate does not prevent trying another metric.

Slot compatibility includes representation, dimension and duality. Substitution
touches direct tensor slots, keeping scalar parameters and index payloads opaque.
The shared borrowed-slot matcher compares dimension and index expressions exactly,
including compound dimensions such as `Nc^2-1` and large integer index labels.
Chain and trace bodies follow the existing chain-like setting. A substitution
does not cross a sum or power, while complete contractions inside those scopes
can simplify independently. Rewriting continues to a fixed point.

Before entering this contractor, the shared explicit-index scan checks whether
any index payload repeats. A negative result skips metric substitution and the
chain/vector contraction rules, which also require repeated explicit slots.
Dot normalization, ordinary vector compaction and epsilon identities retain
their own processing. The scan is conservative: repeated payloads with incompatible
representations still reach the contractor and are rejected there.

The native benchmark separates the raw contractor, the repeated-index scan,
the guarded contractor, the complete Schoonschip pass and gamma simplification.
Inputs include productive contractions, incompatible slots, long products,
chains, cyclic and symmetric traces, and axial outputs of lengths six through
twelve. HEP component tests compare exact tensor-network evaluations before and
after rewriting, including antisymmetry, dual representations and scalar
parameters sharing an index name.

On 2026-09-24, five-sample medians on one pinned CPU of the shared host give the
following isolated contractor timings in microseconds. These exclude the
repeated-index guard and dot/bracket normalization.

#table(
  columns: 4,
  [Input], [Tuple rules], [Metric first], [Speedup],
  [One metric and tensor], [26.75], [2.31], [11.6×],
  [Metric and chain body], [37.08], [4.15], [8.9×],
  [Metric and cyclic trace], [39.73], [4.69], [8.5×],
  [Eight connected metrics], [110.79], [38.24], [2.9×],
  [Miss with 48 spectator factors], [58.08], [43.08], [1.35×],
)

For the already simplified axial-twelve output, the raw failed pass improves from
58.11 to 23.28 ms. The guard reduces the new metric stage to 0.382 ms; putting the
same guard before the old rules takes 0.413 ms. Thus the large terminal-output
gain comes from avoiding work, separately from faster productive contractions.
Complete Schoonschip still takes 67.51 ms, down from 125.44 ms.

The complete gamma pipeline remains much slower than FORM. These timings are in
milliseconds; FORM uses one warmup process followed by three measured processes,
each containing six batches. Its wall time includes amortized startup, input
parsing, tracing, sorting and cleanup. Rust input construction is outside timing.

#table(
  columns: 4,
  [Axial length], [Previous full gamma], [Current full gamma], [FORM `trace4` wall],
  [6], [5.79], [4.78], [0.00634],
  [8], [45.42], [36.15], [0.02023],
  [10], [332.40], [254.64], [0.09993],
  [12], [2444.34], [1791.53], [0.61944],
)

For axial twelve, the full pipeline improves by 1.36 times, while FORM remains
about 2,892 times faster under these timing boundaries. FORM's internal
trace-and-sort CPU estimate is 0.593 ms. This change does not establish FORM parity;
the surrounding normalization and epsilon passes still dominate the full pipeline.

The #source-link("examples/notebooks/metric_contraction.json", label: "metric contraction measurements")
record all 25 cases, raw samples, build/source identities, output diagnostics,
FORM programs and validation results. All 24 HEP tests pass. The broad Spenso and
Idenso suite has 589 passing tests and two preexisting failures; the record also
identifies the unrelated HEP library-unit-test compile failure encountered by
the broad Clippy command. Clippy passes for the changed libraries, examples and
HEP integration tests.

// docs-example: syntax
```sh
cargo run -p idenso --profile dev-optim --example metric_contraction_benchmark -- /tmp/metric-benchmark
cargo nextest run -p spenso-hep-lib --test metric_contraction_validation --cargo-profile dev-optim
```

== Contraction order and intermediate terms

The following investigation records the behavior before the contraction-order
update described below.

For distinct indices, cold Idenso recipe generation and FORM emit the same raw
monomial counts before collection. Warm Idenso calls reuse the collected recipes.

#table(
  columns: 3,
  [Trace], [Generated in either engine], [Collected recipes / FORM output],
  [Axial 6 / 8 / 10 / 12], [6 / 33 / 180 / 1053], [6 / 33 / 180 / 1029],
  [Ordinary 8 / 10 / 12], [117 / 801 / 5139], [105 / 693 / 4383],
)

Axial twelve grows from one trace to 1,029 terms during the first gamma rewrite.
Every later stage, including the whole second outer pass, leaves it unchanged.
Its slow cleanup therefore does not come from generating extra trace branches.
These counters measure cumulative emissions, not peak resident terms or memory.

Contraction order does matter when indices repeat. Write `Tr` for an ordinary
four-dimensional trace, `a,b` for summed Lorentz indices, and `p0..p7` for distinct
slashed vectors. For `Tr[a,p0,b,p1,p2,a,p3,p4,b,p5,p6,p7]`, Idenso selects the
shortest cyclic pair, `a`, whose even interior branches. It eventually evaluates
three length-eight kernels, materializing 315 cached monomials. FORM prioritizes
the odd-interior Chisholm identity: contracting `b`, then `a`, gives
`4 Tr[p3,p4,p0,p2,p1,p5,p6,p7]`. One length-eight kernel emits 117 raw monomials,
collected into 105. This priority is explicit in
#link("https://github.com/form-dev/form/blob/v5.0.0/sources/opera.c#L828")[FORM's `Trace4Gen`].

With external metrics,
`g(a,c)*g(b,d)*Tr[a,p0,b,p1,p2,c,p3,p4,d,p5,p6,p7]`, the default Idenso prepass
leaves trace bodies opaque and materializes 4,383 recipes before contraction.
Enabling the existing chain-like Schoonschip setting before gamma simplification
contracts the metrics first and reduces that count to 315. FORM contracts the
metrics before tracing and then follows its odd-interior priority.

Five warm, interleaved samples on CPU 6 of the shared host measured 24.73 ms for
the repeated-index input versus 8.33 ms after manually applying the first odd
identity. Preparing that shorter input is outside timing; this is not a patched
selector benchmark. For the external-metric input, the default takes 656.42 ms
versus 26.64 ms including the chain-aware prepass. These diagnostic timings
preceded the production update below.

Exact HEP network evaluations agree on three nonzero integer assignments for the
original, both odd reductions, the external-metric and precontracted inputs, and
the repeated-pair and odd-first simplified outputs. The large default
external-metric output was not HEP-evaluated. All 24 instrumented outputs match
the linked production outputs exactly. The
#source-link("examples/notebooks/gamma_trace_intermediates.json", label: "intermediate-count record")
preserves counters, assignments, raw timings, source identities and diagnostic
drivers. Counts preserve factorization; outer sum arity alone hides nested sums.

== Applied contraction ordering

Gamma simplification now contracts metrics into chains and traces before every
rewrite pass. Keeping this active through the fixed point also handles metrics
generated by evaluating a neighboring trace. The existing repeated-index scan
guards both metric substitution and chain/vector rules; terminal expressions
without repeated explicit indices skip those operations without tensor parsing.

Cyclic rotation and Chisholm selection share one priority: adjacent pairs,
applicable odd interiors, applicable even interiors, then other repeated pairs.
Only explicit four-dimensional index pairs with gamma-only interiors get Chisholm
priority. An interval crossing gamma-five cannot hide an applicable interval on
the other side of the cyclic boundary. Open-chain and generic-dimensional
repeated-pair ordering retain their existing behavior.

The repeated-index example now has 105 arithmetic leaves instead of 315; its
external-metric version goes from 4,383 to 105. The axial crossing example goes
from 99 to 33. These counts sum nested additions and multiply product counts
without expanding or collecting the expression. They are output-size diagnostics,
not peak intermediate-memory measurements or equality checks.

The benchmark covers 35 inputs, including 16 complete gamma pipelines, with 12
matching FORM cases. Five warm Idenso samples use the same optimized driver and
one pinned CPU before and after the change. FORM wall times amortize startup,
parsing, tracing, sorting and cleanup over six batches in each of three measured
processes. All runs use the shared host; Idenso input creation is outside timing.

#table(
  columns: 4,
  [Input], [Before (ms)], [After (ms)], [FORM wall (ms)],
  [Repeated-index trace, 12], [19.52], [6.62], [0.06043],
  [External metrics, 12], [508.55], [6.74], [0.06144],
  [Axial crossing, 12], [105.96], [32.22], [0.02116],
  [Two-gamma interior, 8], [0.269], [0.286], [0.00491],
  [Cyclic adjacent pair, 12], [50.98], [51.16], [0.40266],
  [Free indices, 8], [0.148], [0.148], [0.05719],
  [Free indices, 12], [6.24], [6.32], [2.91050],
  [Axial distinct indices, 12], [1769.61], [1785.88], [0.62148],
)

The two main ordinary cases improve by 2.95 and 75.4 times; the axial crossing
case improves by 3.29 times. The small two-gamma-interior case costs about 17
microseconds more. Distinct-index axial twelve remains dominated by cleanup,
and FORM remains substantially faster on the contracted examples.

Ten new unit regressions cover cyclic priority, gamma-five boundaries, symbolic
dimensions, inert traces and idempotence after generated metrics. Three new HEP
regressions compare original and simplified networks on nonzero exact integer
assignments, including the production external-metric output. The
#source-link("examples/notebooks/gamma_trace_ordering.json", label: "implementation measurements")
preserve raw samples, all matching FORM programs, source identities, factorized
output counts and validation results. Reproduce the native measurements with the
metric contraction benchmark command above.

== HEP tensor-component validation

Four-dimensional identities allow different symbolic expressions for the same tensor.
Term counts are output-size diagnostics, never the criterion for equivalence.
The notebook contracts original gamma traces using the HEP library's explicit
four-by-four matrices, and contracts simplified metric/epsilon networks against the
same exact integer momenta. Each sample spans four dimensions; nonzero epsilon
contractions prevent gamma-five checks from passing for degenerate kinematics.
The conventions are the (+---) metric and epsilon component `epsilon(0,1,2,3) = -i`.

For lengths four, eight and ten, ordinary and gamma-five checks compare these HEP
values. When FORM is available, its ordinary `trace4` and `tracen` polynomials are
imported as metric networks and evaluated using the same library. The ten-factor
case demonstrates agreement despite 693 versus 945 terms. Finite component samples
provide regression evidence, not a symbolic proof for all tensors.

The dedicated Rust HEP tests extend this coverage to lengths one through fourteen,
several gamma-five positions, repeated slashes and cyclic contracted Lorentz indices.
They compare exact Gaussian-integer results. To keep the tests small, repeated scalar
metric/epsilon subnetworks reuse their HEP contraction values across the polynomial;
the original gamma-matrix network is contracted independently.

// docs-example: syntax
```sh
cargo nextest run -p spenso-hep-lib --test short_trace_validation --cargo-profile dev-optim
```

== Comparing the first axial rewrite directly with FORM

For gamma-five followed by twelve ordinary gammas with distinct Lorentz indices,
Idenso's first rewrite and FORM 5.0.0 `trace4` each give 1,029 terms. The outputs
share 741 identical terms and have 288 terms unique to each. Their raw symbolic
difference has 576 terms. This example makes the non-uniqueness concrete.

Mapping FORM's `d_` and `e_` to the HEP metric and epsilon, three full-rank integer
momentum assignments give identical exact values for the original gamma network,
Idenso's first rewrite and FORM's polynomial: `-96536i`, `468485024i` and
`114798012i`. The first rewrite also equals Idenso's full production output exactly.
The epsilon convention follows the
#link("https://form-dev.github.io/form-docs/master/manual/#a-few-notes-on-the-use-of-a-metric")[FORM manual's trace normalization];
no sign or factor is fitted to the twelve-gamma result.

A fresh optimized run measures 2.42 ms for the first rewrite and 0.654 ms amortized
FORM wall time, about 3.7 times faster for FORM. FORM process latency is 11.97 ms,
while its internal trace-and-sort CPU time is 0.624 ms. The first rewrite timing
excludes Idenso's later passes; the complete pipeline previously measured 2.56 s.
Both batch and in-process measurements are on a shared host, without CPU pinning.

The #source-link("examples/notebooks/gamma_trace_axial_form.json", label: "axial FORM comparison record")
contains the exact assignments, result values, raw timing samples, source hashes
and executable FORM programs. These finite component checks support equivalence
without requiring a common symbolic normal form.

== Run locally

Use an interpreter containing the combined Symbolica community host with Spenso and Idenso,
Marimo 0.24.0, and Typst 0.15.0. Installing ordinary Symbolica alone does not supply the
community extension.

// docs-example: syntax
```sh
just notebook gamma_simplification /path/to/python
```

Put `form` on `PATH`, or select its executable explicitly for native comparisons:

// docs-example: syntax
```sh
FORM_EXECUTABLE=/path/to/form just notebook gamma_simplification /path/to/python
```
]
