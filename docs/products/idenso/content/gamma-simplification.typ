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

`GammaSimplifySettings(expand_traces=True)` requests expanded evaluated trace
bodies in Python; Rust uses `with_expanded_traces()`. The default remains factored.
Surrounding scalar factors and independent trace boundaries are preserved, so
`S * trace` becomes `S * expanded_body` without distributing S. Disabling trace
evaluation also disables trace expansion. The ordinary-trace example checks the
result by expanding only the standalone body and then verifies an unchanged rerun.

See #product-link("idenso", page: "reference/form-color-dirac/", label: "Shipped color and Dirac rules")
for the conventions and #product-link("idenso", page: "guides/showcase/", label: "the tensor display showcase")
for additional notation controls.

== Schoonschip performance and FORM comparison

The notebook also includes a massless three-loop propagator with a fermionic
outer ring, followed by four-loop fermionic and gluonic ladders. The new cases
share a physical routing: each vertex conserves momentum and the two incoming
momenta on every internal edge are opposite. External vector indices are
contracted with a metric. This differs from the historical gluonic example's
independent external polarization vectors, so their timings and term counts
are separate benchmarks.

The #source-link("examples/notebooks/fermion_ladder.py", label: "shared benchmark driver")
uses the installed public tensor API and runs handwritten/generated native
FORM programs. Fermionic cases compare strict four-dimensional `trace4` with
symbolic-D `tracen`, both with spinor trace normalization four. The gluonic
case compares four-dimensional and symbolic-D metric algebra, without a
Dirac trace. All cases omit couplings, color, denominators and integration;
the closed-fermion-loop minus sign is included in the fermionic numerators.

Complete algebra clocks include tracing or vertex replacement, contraction,
momentum routing and final scalar polynomial collection. Input construction,
prepared routing rules, component checks and output are excluded and setup
is recorded separately. Independent expressions are evaluated in alternating
warm batches on the same CPU; FORM's internal CPU timer avoids process-launch
overhead. Idenso's Python dispatch and typed result wrapping remain included.
The raw records retain CPU and wall samples, source and executable hashes,
FORM diagnostics and exact checks. Different batch sizes imply different
numbers of simultaneously retained results.

The fermionic route keeps momenta compact during the trace and inserts their
routing into scalar products afterward. Three-loop four-dimensional reduction
uses simultaneous dot substitution and polynomial expansion. The larger cases
convert the scalar trace to a polynomial once, group coefficients of all dots
touching one momentum, substitute that momentum simultaneously, and collect
before proceeding. Opposite outer momentum pairs are eliminated first; the final
Atom is emitted once. This uses existing Symbolica primitives after typed gamma
simplification, without bypassing tensor-interface checks.
Expanded traces now also use direct polynomial emission for contracted
sixteen-factor words. This reuses the existing Clifford recurrence and
collects its integer coefficients before constructing the scalar Atom.
The coefficient and exponent row buffers have a checked 16 MiB limit;
exceeding it evaluates the same completed recipe with its already normalized
metrics. Free sixteen-factor traces and callback-sensitive mixed free-index
and compact-vector inputs retain their established evaluation path.
FORM collects after each momentum substitution;
expanding all substitutions before collecting was substantially slower.
These are measured schedules, not claims of globally optimal contraction
ordering. Each final polynomial is checked exactly against FORM, including
symbolic D, and D=4 must agree with strict four-dimensional reduction. Three
rational assignments independently compare original HEP tensor networks and
scalar results; fermionic checks additionally use exact matrix contractions
of the HEP library's gamma components.

On 2026-09-27, five alternating warm batches per mode on CPU 8 gave these
complete algebra median CPU times with FORM 5.0.0 and the optimized installed
host `4eace950`. Ratios divide the two engines' medians. The machine is shared;
these are workload-specific measurements.

#table(
  columns: 5,
  [Numerator], [Dimension], [Idenso, ms], [FORM, ms], [Ratio],
  [Three-loop fermion], [4D], [0.612], [0.356], [1.72×],
  [Three-loop fermion], [D], [2.879], [2.231], [1.29×],
  [Four-loop fermion], [4D], [6.523], [6.583], [0.99×],
  [Four-loop fermion], [D], [133.257], [71.500], [1.86×],
  [Four-loop gluon], [4D], [773.030], [778.000], [0.99×],
  [Four-loop gluon], [D], [1,410.728], [1,479.000], [0.95×],
)

The three-loop case meets the factor-three threshold in both dimensions;
the four-loop comparisons follow that gate. The gluonic route uses
tensor-safe replacement and local contraction against opaque vertices,
with different measured rung-closing orders in the two engines. The
remaining larger gap is the four-loop D-dimensional fermion trace: one
diagnostic call spends 88.6 ms tracing, 8.8 ms converting to a polynomial,
29.7 ms substituting momenta, and 7.4 ms emitting the final scalar. These phase
clocks are separate from the complete-call medians.

Seven alternating paired batches against the saved simultaneous-routing driver
show median paired CPU reductions of *27.5%* for three-loop D, *18.2%* for
four-loop 4D, and *53.9%* for four-loop D. Each changed case is faster in all
seven pairs. Three-loop 4D retains its previous route. Preparing the changed
cases' rules costs about 0.5 ms more in that cohort, outside the algebra clocks.
The #source-link("examples/notebooks/fermion_ladder_routing_comparison.json", label: "routing comparison")
retains every paired sample and the baseline driver source. The gluonic route
is unchanged; its fresh measurements are controls. Shared-host absolute times
vary between cohorts and are not substituted for paired changes.

The #source-link("examples/reproducers/symbolica-expansion/polynomial_replacement.rs", label: "standalone polynomial replacement reproducer")
isolates another Symbolica cost: its current replacement rebuilds the right-hand
side power for each monomial, then repeatedly adds to the growing polynomial.
Grouping coefficients and using sparse Horner evaluation avoids that repeated
work. The reproducer uses synthetic scalar polynomials and unchanged Symbolica;
its substitution timings are separate from the complete ladder results above.

All six exact FORM polynomials, three D-to-four specializations, and 72
HEP component comparisons pass. Full measurements are in the
#source-link("examples/notebooks/fermion_ladder_three_staged_timing.json", label: "three-loop fermionic record"),
#source-link("examples/notebooks/fermion_ladder_four_staged_timing.json", label: "four-loop fermionic record"), and
#source-link("examples/notebooks/gluon_ladder_four_staged_control.json", label: "four-loop gluonic control").
The earlier simultaneous-routing records remain available beside these files.

// docs-example: syntax
```sh
python examples/notebooks/fermion_ladder.py --form /path/to/form \
  --loops 3 --rounds 5 --cpu 8 --output /tmp/fermion-three.json
python examples/notebooks/fermion_ladder.py --form /path/to/form \
  --loops 4 --rounds 5 --cpu 8 --output /tmp/fermion-four.json
python examples/notebooks/fermion_ladder.py --form /path/to/form \
  --particle gluonic --loops 4 --rounds 5 --cpu 8 --output /tmp/gluon-four.json
```

For the historical gluonic fixture, the qualified checkpoint takes *0.858 s* in the original order
and *0.300 s* with early rung contractions. Both use tensor-safe replacement
and include final expansion. Paired gains over the saved build are 9.09% and
6.66%; current FORM comparisons and mixed trace controls are reported in
#link(<resumable-observation-public>)[the public qualification below]. These
results do not establish full FORM parity.

The performance section follows equivalent paired and alternating slash traces through
three routes: evaluate a free-index trace then contract momenta; contract the indexed
input with `schoonschip_net()` before taking the trace; or start from the compact slash
expression. Each route must match an independent scalar identity. The displayed speedup
compares the first two complete pipelines, including early conversion cost. The compact
input timing is reported separately. Changing the display layout alone does not change
the algebra or its cost.

These notebook trace examples use three timed runs after a warm-up, rotate execution
order, and exclude construction, assertions, and rendering. They characterize the
installed build and machine; they are not a release benchmark of the engines.

Closed traces now reuse adjacent contractions and the four-dimensional Chisholm
identities from open chains. A macro dispatches lengths 1–14 to generated short-trace
kernels, with separate ordinary and gamma-five tables and immediate zero for odd words.
The integer recipes follow the same three-gamma epsilon identity as open chains and
are initialized lazily. Standalone ordinary or single-gamma5 traces with distinct
explicit 4D Lorentz slots can return their terminal metric/epsilon output directly;
surrounding contractions retain the full tensor simplification pass.

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
was the first optimization target; the applied cleanup below addresses it. A
fourteen-factor alternating slash trace is different again: eleven outer passes take 9.66 ms, including 3.76 ms in
Schoonschip, 2.02 ms in collection and 2.50 ms in gamma rewriting.

For the standalone ordinary kernel, borrowed factor views remove avoidable copies:
a separate five-sample experiment changes length fourteen from 54.3 to 48.3 ms,
about eleven percent. Symbolica still normalizes each product and merges the sum.
Sampling highlights product normalization, byte comparisons and heap operations.
The #source-link("examples/notebooks/gamma_trace_profile.json", label: "detailed profile record")
preserves raw timings, source hashes, diagnostic variants and their limits.
Production algebra was unchanged during this investigation.

== Earlier metric-first contraction

This section records the first contractor update, before the shared slot
substitution and epsilon cleanup below. At this stage, Schoonschip selected one
metric with `g(a_,b_)*remainder___`, then walked the captured remainder to find a
compatible explicit tensor slot. It substituted that slot with the opposite
metric endpoint and removed the metric. The matcher captured the entire product
at the current tree level; it did not search for a metric and its tensor partner
or partition the remainder into subsets. A failed candidate did not prevent
trying another metric.

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
Complete Schoonschip at this stage takes 67.51 ms, down from 125.44 ms.

At this stage, the complete gamma pipeline remains much slower than FORM. These
timings are in milliseconds; FORM uses one warmup process followed by three measured processes,
each containing six batches. Its wall time includes amortized startup, input
parsing, tracing, sorting and cleanup. Rust input construction is outside timing.

#table(
  columns: 4,
  [Axial length], [Before metric first], [After metric first], [FORM `trace4` wall],
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
FORM programs and validation results from that stage. All 24 HEP tests passed.
The broad Spenso and Idenso suite had 589 passing tests and two preexisting failures; the record also
identifies the unrelated HEP library-unit-test compile failure encountered by
the broad Clippy command. Clippy passes for the changed libraries, examples and
HEP integration tests.

// docs-example: syntax
```sh
cargo run -p idenso --features reference-cases --profile dev-optim --example metric_contraction_benchmark -- /tmp/metric-benchmark
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

== Earlier contraction-order update

The next update made gamma simplification contract metrics into chains and
traces before every rewrite pass. Keeping this active through the fixed point also handles metrics
generated by evaluating a neighboring trace. The existing repeated-index scan
guards both metric substitution and chain/vector rules; terminal expressions
without repeated explicit indices skip those operations without tensor parsing.

Cyclic rotation and Chisholm selection share one priority: adjacent pairs,
applicable odd interiors, applicable even interiors, then other repeated pairs.
Only explicit four-dimensional index pairs with gamma-only interiors get Chisholm
priority. An interval crossing gamma-five cannot hide an applicable interval on
the other side of the cyclic boundary. Open-chain and generic-dimensional
repeated-pair ordering retain their existing behavior.

After this update, the repeated-index example has 105 arithmetic leaves instead
of 315; its external-metric version goes from 4,383 to 105. The axial crossing example goes
from 99 to 33. These counts sum nested additions and multiply product counts
without expanding or collecting the expression. They are output-size diagnostics,
not peak intermediate-memory measurements or equality checks.

These historical measurements precede the cleanup below. The benchmark covers
35 inputs, including 16 complete gamma pipelines, with 12 matching FORM cases.
Five warm Idenso samples use the same optimized driver and
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
microseconds more. At this stage, distinct-index axial twelve remains dominated
by cleanup, and FORM remains substantially faster on the contracted examples.

Ten new unit regressions cover cyclic priority, gamma-five boundaries, symbolic
dimensions, inert traces and idempotence after generated metrics. Three new HEP
regressions compare original and simplified networks on nonzero exact integer
assignments, including the production external-metric output. The
#source-link("examples/notebooks/gamma_trace_ordering.json", label: "implementation measurements")
preserve raw samples, all matching FORM programs, source identities, factorized
output counts and validation results. Reproduce the native measurements with the
metric contraction benchmark command above.

== Earlier shared slot substitution and epsilon cleanup

These measurements record the shared-slot update before the local metric
normalization below.

Metrics and tagged vectors now use the same slot substitution. The contractor
scans product factors directly, trying metric sources before tagged rank-one
vectors; finding a source no longer uses the metric/remainder pattern. It joins
one compatible explicit slot to the opposite metric endpoint or compact vector,
then removes the source factor. Representation, dimension and duality checks,
opaque scalar parameters, and the chain-like and rank-one settings are retained.
The repeated-index scan still guards this work, and substitutions preserve sums
and powers rather than distributing them.

Epsilon cleanup remains a separate `simplify_epsilon()` operation. When epsilon
is present, it retains one initial whole-expression Schoonschip pass. Candidate
selection directly distinguishes powers, products and epsilon function heads.
It no longer runs Schoonschip at every visited node: another whole-expression
pass runs only after an epsilon identity changed the expression, to consume the
generated metrics. Epsilon powers and pairs with disjoint explicit indices still
receive their identities even when the repeated-index predicate is false.

The fresh comparison uses the same benchmark driver against the saved baseline
and the new libraries: 50 inputs and 316 case/method pairs. All before/after
output snapshots are byte-identical. Each Rust timing is a median of five warm
samples, with adaptive repetitions targeting 8 ms per sample, pinned to CPU 6
of the shared host; input construction is outside timing. The following
Schoonschip timings include normalization and the repeated-index guard, in
microseconds.

#table(
  columns: 4,
  [Input], [Before], [After], [Speedup],
  [Tagged vector and tensor], [44.33], [19.51], [2.27×],
  [Tagged vector and trace], [62.62], [29.76], [2.10×],
  [Tagged compact metric], [49.34], [24.46], [2.02×],
)

Isolated epsilon cleanup improves from 16.401 to 1.418 ms for a pair with eight
distinct explicit indices (11.57 times), and from 827.076 to 58.705 ms on the
already simplified axial-twelve output (14.09 times). These measurements include
the retained initial Schoonschip pass; they do not change the trace kernel or
reduce its final 1,029 terms.

Full gamma measurements and a fresh FORM 5.0.0 run give the following
milliseconds. FORM uses one warmup process and three measured processes, with
six batches each, pinned to CPU 3. Its amortized wall time includes startup,
parsing, input normalization, `trace4`, sorting and cleanup.

#table(
  columns: 4,
  [Input], [Before], [After], [FORM wall],
  [Axial distinct indices, 6], [4.854], [0.637], [—],
  [Axial distinct indices, 8], [36.759], [4.895], [—],
  [Axial distinct indices, 10], [256.810], [35.781], [—],
  [Axial distinct indices, 12], [1805.957], [259.161], [0.61493],
  [Axial crossing, 12], [32.130], [3.103], [0.02139],
  [External metrics, 12], [6.846], [6.105], [0.06142],
  [Free indices, 8], [0.153], [0.148], [0.05617],
  [Free indices, 12], [6.759], [6.477], [2.87345],
)

Axial twelve improves by 6.97 times, and the axial crossing case by 10.36 times.
Ordinary standalone traces change little because they already bypass most of
this cleanup. FORM remains about 421 times faster for axial twelve under these
timing boundaries; this does not establish parity. The
#source-link("examples/notebooks/gamma_cleanup.json", label: "shared-slot and epsilon-cleanup measurements")
record all cases, raw samples, output comparisons, build/source identities and
the fresh FORM programs. A dash marks a case without a fresh matching FORM run.

Normalization remains costly on the unchanged axial-twelve output. In an
instrumented copy, one `normalize_dots()` call takes 28.45 ms, including 17.13 ms
in redundant-metric rules and 6.36 ms in metric-trace rules. The epsilon candidate
pass takes 1.09 ms and the repeated-index scan 0.398 ms; complete Schoonschip
takes 59.80 ms. These are isolated five-round warm medians on CPU 6, with outputs
checked against production. They do not add up to a decomposition of the full
gamma pipeline.

The HEP coverage retains one known limitation: an explicit vector product with
dual slots evaluates, but both the saved baseline and the new compact-dot output
discard the dual orientation, so the tensor-network parser rejects that output.
The diagnostic `tagged_vector_dots_preserve_dual_slot_orientation` remains
explicitly ignored with this reason; this case is not certified by the passing
component checks. The measurement record identifies the limitation separately
from the unchanged benchmark snapshots.

== Applied local metric normalization

The measurements in this section use crates.io Symbolica 3.0.0. The main-branch
comparison below retains this normalizer and changes the pinned dependency.

The previous axial-twelve profile spent 17.128 ms checking redundant metrics
and 6.361 ms checking metric traces inside one 28.448 ms `normalize_dots()`
call. Every result was unchanged. These stages repeatedly walked an already
normalized expression and attempted wildcard matches at eligible metric heads.

Intrinsic metric identities now run locally when Symbolica constructs `g`:
compatible repeated abstract slots yield their dimension, and a compatible
explicit slot with a tagged compact vector yields that vector component.
Scalar linearity still applies. `normalize_dots()` handles a root-level identity
directly, then uses a read-only visitor to check for remaining nested-vector or
integral-power identities. If none applies, it returns the existing expression
without rebuilding it. Otherwise, one rewriting walk uses borrowed slot
recognition instead of whole-expression replacement rules. Scalar parameters
and slot payloads remain opaque, and numerator sums stay factored.

The read-only check matters even after removing wildcard rules: assigning an
unchanged node in Symbolica's rewriting walker to prune its children marks that
path as changed and rebuilds its ancestors. The visitor can prune the same
opaque payloads without triggering reconstruction or normalization.

Recognition requires canonical slot shapes, exact dimensions, compatible
representations and variance, and a final structural argument for vectors.
Malformed arguments and dimension mismatches remain explicit. Fixed components
marked by `cind` or `find` are excluded from abstract-index traces and power
contractions. The existing convention for abstract powers is retained: negative
even powers normalize; negative odd and nonintegral powers remain explicit.

Construction also makes scalar composition consistent: `g(i,i) = D` implies
`g(i,i)^2 = D^2` and `g(i,i)^3 = D^3`. The former order applied the metric-power
rule before the trace rule and incorrectly gave `D` for the square. The
distinct-index identity `g(i,j)^2 = D` remains a contraction. Exact HEP
regressions compare against a separate tensor with the metric's components and
no normalization callback, so the trace-power check does not merely compare
two expressions that already simplified during construction.

The expanded benchmark retains identical source strings and separately times
parsing, `normalize_dots()` on constructed atoms, and parsing followed by
normalization. Parse-only includes constructor callbacks; parse-plus-normalize
shows their total cost instead of hiding work moved out of the later pass.
There are 73 inputs and 641 case/method pairs. Each Rust median uses five warm
samples, rotating method order and adapting repetitions toward 8 ms per sample,
on CPU 6 of the shared host. The following construction-inclusive measurements
are in microseconds.

#table(
  columns: 5,
  [Input], [Parse before], [Parse after], [Parse + normalize before], [Parse + normalize after],
  [Metric self-trace], [3.105], [3.822], [12.391], [3.860],
  [Metric and tensor], [4.880], [5.471], [16.160], [6.051],
  [Metric and compact vector], [3.713], [4.737], [12.107], [4.971],
  [Mismatched dimensions], [3.124], [3.633], [13.148], [4.057],
  [Compact dot with scalar factors], [5.751], [6.365], [11.037], [6.880],
)

Parsing is more expensive in these examples because local metric checks now
happen there. The combined operation improves even after including that cost.
Constructed input snapshots can change eagerly, while the recorded source
strings stay identical. The trace-power correction and stricter canonical
recognition also intentionally change some outputs; the record distinguishes
those cases from unchanged output comparisons.

On the unchanged axial-twelve output, isolated `normalize_dots()` improves from
27.794 to 0.858 ms (32.39 times), and complete Schoonschip from 59.577 to 5.265 ms
(11.32 times). Full gamma timings below exclude initial input construction but
include metrics constructed while evaluating the trace. Matching FORM 5.0.0
wall times include amortized startup, parsing, tracing, sorting and cleanup:
one warmup process, then three measured processes with six batches each on CPU 3.
All entries are milliseconds.

#table(
  columns: 4,
  [Input], [Before], [After], [FORM wall],
  [Axial distinct indices, 6], [0.647], [0.204], [—],
  [Axial distinct indices, 8], [4.952], [1.141], [—],
  [Axial distinct indices, 10], [36.379], [7.161], [—],
  [Axial distinct indices, 12], [264.302], [49.084], [0.62948],
  [Axial crossing, 12], [3.268], [1.643], [0.03008],
  [External metrics, 12], [6.291], [3.524], [0.06158],
  [Free indices, 8], [0.150], [0.162], [0.05679],
  [Free indices, 12], [6.622], [6.496], [2.91656],
)

Axial twelve improves by 5.38 times. Free eight becomes about 7.9 percent slower;
that shortcut already avoids cleanup, while constructing already-normal metrics
can incur the new checks. Free twelve changes little. FORM remains about 78
times faster for axial twelve under these timing boundaries. The
#source-link("examples/notebooks/dot_normalization.json", label: "local metric normalization record")
preserves all cases, raw samples, input/output snapshot hashes, classified changed
outputs, source identities and matching FORM programs. The benchmark command
recreates the complete snapshots. A dash marks a case without a fresh matching
FORM run.

Eight dot-normalization calls process the large axial-twelve output. Multiplying
the isolated per-call reduction by eight estimates 215.5 ms saved, close to the
215.2 ms reduction in the complete pipeline. This is a timing-based attribution
estimate, not an exact measurement of those calls inside the pipeline.

A separate final profile uses five warm samples on CPU 7 and checks that every
stage preserves the terminal expression. It measures 0.867 ms for dot
normalization, 5.330 ms for Schoonschip and 6.722 ms for epsilon cleanup. In that measurement, chain
collection was the largest stage: `chainify` takes 20.957 ms and
`collect_gamma_chains` takes 24.475 ms; simplifying the already-terminal
expression takes 39.275 ms. These stage measurements overlap and cannot be added
to reconstruct full gamma runtime.

The selected validation suites have 632 passing tests: 336 in Idenso, 261 in
Spenso with shadowing enabled, and 35 HEP component tests. Two previously recorded
tests still fail, and 23 are skipped, including the known compact dual-dot
variance diagnostic above. The record names those limitations and retains the
commands and logs. All three scoped Clippy checks pass with warnings denied.

== Pinned Symbolica main comparison

The active Rust workspaces and Python lock now pin Symbolica main at
`06906976bca24fefc5203aee699d90d62ebe08cd`. Its package version still reads
3.0.0, so the source revision distinguishes it from the registry release above.
The same benchmark driver covers 73 inputs and 641 case/method pairs. All 787
plain-text source, input and output snapshots match the optimized registry
baseline exactly.

The following full gamma medians are milliseconds. Both executables were
measured afresh, registry first and pinned main second, with five warm samples
on CPU 6. Our validation compiler group was paused for both runs and resumed
afterward. Other shared-host load remains uncontrolled; builds were not
interleaved.

#table(
  columns: 3,
  [Input], [Registry 3.0.0], [Pinned main],
  [Axial distinct indices, 8], [1.145], [1.091],
  [Axial distinct indices, 10], [7.172], [6.870],
  [Axial distinct indices, 12], [48.962], [46.879],
  [Axial crossing, 12], [1.645], [1.621],
  [External metrics, 12], [3.524], [3.628],
  [Free indices, 8], [0.160], [0.159],
  [Free indices, 12], [6.399], [6.726],
)

The differences are small and mixed. Axial-twelve dot normalization takes
0.889 versus 0.848 ms, and its full pipeline is about 4.3 percent faster.
Free twelve is about 5.1 percent slower; the external-metric case is about
3 percent slower. This sequential shared-host comparison does not establish a
broad upstream speedup.
The #source-link("examples/notebooks/symbolica_main_update.json", label: "Symbolica-main comparison record")
preserves both sets of raw samples, source and lock identities, build metadata,
snapshot hashes and validation results. FORM values in that record are reused
from the preceding comparison; FORM was not rerun for this dependency update.

This revision changes the binary atom export format from 0 to 1 and rejects
legacy format-0 payloads. The surrounding full state-export format also moves
from version 5 to 6, requiring compatible payload consumers. To migrate saved
atoms, export expressions as strings using the older version and parse them
using the new version. Exact agreement
of the printed benchmark outputs does not imply binary-format compatibility.

The standalone native notebook host and Tydenso `wasm32` dependency checks pass
at this revision. It includes the upstream fix for decoding a stored 64-bit
length on a 32-bit target. Tydenso also passes `preserve_indices = false` at its
Dirac-adjoint call site, preserving its previous behavior. Separate native and
Wasm runtime probes check binary round trips, malformed lengths and legacy
format rejection; exact cross-platform byte equality covers the zero fixture.
Fifteen existing Idenso snapshots were updated only for imaginary-unit display;
their inputs and algebra are unchanged. The selected native suites retain 632
passing tests, the same two known failures and 23 skipped tests. All three
scoped Clippy checks pass with warnings denied. These native results do not
certify the separate Tydenso/Tymbolica runtime payload exchange.

That runtime exchange remains blocked: the pinned atom-payload dependency
`18e04916` accepts outer export version 5 and rejects version 6 before import.
The record keeps the reproduced failure and an isolated candidate patch;
the candidate is not applied to the workspace dependency, and the tracked
plugin assets retain their previous build. In isolation, the candidate passes
22 payload tests against each of registry and pinned-main Symbolica, 21
Tydenso runtime checks, and two-way exchange with the rebuilt Tymbolica engine.

== Further contraction and cleanup improvements

The next optimization keeps Symbolica pinned at `06906976` and changes how
Idenso schedules contractions and cleanup. A product finishes its local
contractions before rebuilding enclosing expressions, while preserving the
metric-first order. Compact vectors are constructed only after a compatible
partner is found. Borrowed visitors avoid rebuilding unchanged slot payloads;
chain, bracket and gamma passes skip work when their required symbols are absent.
Schoonschip repeats post-contraction cleanup only when a contraction changed
its input; its outer fixed point still handles initial normalization changes.

The existing terminal-trace shortcut now also accepts a standalone 4D spin
trace with exactly one gamma5 and at most fourteen ordinary gammas, whose
Lorentz indices are distinct explicit `mink(4,index)` slots. Each kernel
monomial uses every index once and contains at most one epsilon, so it needs
no further contraction cleanup. Repeated indices, compact slashes and
spectators retain the full pipeline. The ordinary free-trace shortcut is
unchanged in scope.

Full gamma medians below are milliseconds. These are five warm samples on
CPU 6, using the reversed after-before confirmation run. FORM 5.0.0 was
remeasured on CPU 3. Its amortized wall timing includes input preparation and
process overhead; its separate internal timer covers `trace4` plus sorting.
Rust input preparation is outside timing, and output destruction is included.
Our builds were quiet, but unrelated shared-host load was uncontrolled.

#table(
  columns: 4,
  [Input], [Before], [After], [FORM wall],
  [Axial distinct indices, 8], [1.112], [0.119], [0.021],
  [Axial distinct indices, 10], [6.794], [0.467], [0.102],
  [Axial distinct indices, 12], [46.343], [2.512], [0.646],
  [External metrics, 12], [3.332], [0.930], [0.062],
  [Free indices, 8], [0.162], [0.166], [0.060],
  [Free indices, 12], [6.542], [6.776], [2.939],
)

Axial twelve improves by *18.45 times*. FORM remains *3.89 times faster*
on that case, so this does not establish performance parity. Ordinary free
eight and twelve are 2.3% and 3.6% slower in this confirmation run; both were
slightly faster in the first pair. Those small mixed changes are consistent
with shared-host variation, rather than an established gain.

Productive public Schoonschip calls improve too; the following medians are
microseconds from the same confirmation run.

#table(
  columns: 3,
  [Contraction], [Before (µs)], [After (µs)],
  [Vector into tensor], [5.115], [3.356],
  [Metric into epsilon], [7.758], [4.935],
  [Metric into chain], [9.294], [6.129],
  [Metric into trace], [9.838], [6.563],
  [Eight-metric chain], [26.562], [23.753],
)

The original before-after run is retained alongside the confirmation. It
measured axial twelve at *46.658 → 2.383 ms*. Unchanged parsing controls
exposed a localized timing disturbance in several small cases, which motivated
the reversed repeat; the confirmation resolves those apparent regressions.
No confirmation case above one microsecond regresses by more than 10%.

A separate public-stage profile covers *92 expression shapes and 860
case/method pairs*, including both original and already-terminal traces.
On the terminal axial-twelve output, chain collection takes *23.635 → 1.209 ms*,
Schoonschip *5.280 → 1.567 ms*, epsilon cleanup *6.168 → 2.537 ms*, and full
gamma simplification *37.017 → 7.038 ms*. These overlapping stage timings
cannot be added to reconstruct the original trace call. The shortcut explains
why simplifying the original axial trace is now faster than processing its
large terminal output through the general pipeline.

Both full comparisons cover *73 inputs and 641 case/method pairs*. All
*787 full-driver snapshots and 952 phase snapshots match exactly*. The fresh
FORM run covers fifteen cases; twelve full polynomials match the historical
records exactly. The additional axial-six/eight/ten cases have fresh
counts and timings, without a new FORM component certificate. Term counts
are never used as an equivalence test.

The current validation has *383 passing tests*: 344 Idenso tests and 39 HEP
component tests. The known Idenso parsing snapshot failure remains; 23 tests
are skipped. Spenso's separate suite was not rerun in this round. Both scoped
Clippy checks pass with warnings denied. The new HEP checks simplify explicit
axial traces before attaching spectators, then independently evaluate their
metric/epsilon tensors against the same vectors as the original gamma-matrix
network. Lengths four, eight and twelve, three full-rank samples, and several
gamma5 positions preserve exact Gaussian-integer values.

The first length-twelve oracle timed out while materializing a tensor with
twelve open Lorentz slots. Binding each vector inside its raw gamma before
HEP takes the matrix trace avoids that intermediate without using gamma
simplification in the oracle. The record retains the timeout and the final
passing checks.

The #source-link("examples/notebooks/contraction_performance.json", label: "contraction performance record") preserves both run orders, raw samples, source identities, FORM programs and validation evidence. Reproduce the current Rust measurements with:

// docs-example: syntax
```sh
cargo run --locked -p idenso --features reference-cases --profile dev-optim --example metric_contraction_benchmark -- /tmp/contraction-full 5 8
cargo run --locked -p idenso --profile dev-optim --example contraction_phase_benchmark -- /tmp/contraction-full /tmp/contraction-phases 5 8
```

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

== Generic-dimensional traces and contraction progress

The 2026-09-24 comparison uses symbolic Lorentz dimension $D$ and independent
spin trace normalization `Tr(1)=4`. Canonical ordinary traces now reduce summed
pairs inside the memoized pairing evaluator, including branches with two or
more interior gammas. Adjacent pairs contribute $D$ and one-gamma sandwiches
contribute $2-D$. Longer sandwiches use the Clifford recurrence, with the
existing four-dimensional odd-interior Chisholm identity when applicable.
Branches share subword results and construct metrics only when needed, avoiding
intermediate trace atoms and repeated general cleanup.

This direct route requires canonical slots, uniform dimension and at most two
occurrences of each explicit index. Fixed components are not summed. Compact
slashes also reuse adjacent squares and the repeated-slash identity inside the
evaluator. Scalar spectators retain their factorization. Gamma-five remains
strictly four-dimensional, and free 4D traces retain their factored recipes.

Warm times below are milliseconds with `dev-optim`, Symbolica main `06906976`
and FORM 5.0.0 on the shared EPYC 9754 host, pinned to CPU 9. Idenso measures
in-process wall time including output destruction; FORM measures internal
`tracen` plus sorting CPU time over independent expressions. Both exclude
parsing and startup. Idenso uses three processes per version, five warm batches
each, in alternating order; FORM uses three processes with three retained
batches each. Other host activity is uncontrolled, so small differences should
not be treated as established improvements.

#table(
  columns: 3,
  [Free gammas], [Idenso factored], [FORM expanded],
  [8], [0.145], [0.036],
  [10], [0.483], [0.337],
  [12], [1.997], [3.900],
  [14], [18.184], [53.667],
)

*Expanded-output performance is not at parity.* Factored Idenso results defer
substantial work. A separate five-sample experiment measures trace construction,
expansion and destruction of both outputs:

#table(
  columns: 5,
  [Generic-D input], [Trace, ms], [Expansion, ms], [Full lifecycle, ms], [FORM, ms],
  [8 gammas, branching pair], [0.058], [0.050], [0.109], [0.0134],
  [Pair with four interior gammas], [0.441], [1.295], [1.731], [0.173],
  [Pair with five interior gammas], [2.230], [19.728], [21.894], [2.550],
  [12 gammas, order-sensitive pairs], [1.064], [4.892], [5.895], [0.763],
  [14 free gammas], [39.047], [739.466], [789.138], [53.667],
)

The fourteen-gamma polynomial has 135,135 terms and takes about fifteen times
FORM's time including expansion. Each column is a median of its measured
quantity; the total includes destruction. This separate experiment has a
different allocator state from the first-call table. Diagnostic expansion is
restricted to standalone trace polynomials; graph numerators and spectators
retain their factorization. The earlier checkpoint's 776 ms is comparable, not
an attributed optimization. FORM's earlier length-ten samples ranged from
0.337 to 0.503 ms, another reason to avoid a parity claim there.

A separate six-case expansion profile isolates Symbolica's expression work.
At length 14, expansion alone takes 756.9 ms with tensor slots, 635.3 ms with
shallow metrics and 580.2 ms with scalar variables. Using polynomial expansion
reduces the variable case to 414.2 ms, but converting it back to metrics takes
another 353.6 ms, before forward substitution. All 1,080 inverse-substitution
checks pass. Removing tensor syntax therefore does not provide an end-to-end
improvement. The length-12 instruction profile attributes about 19.8% to
normalization, 11.7% to product/sum buffer extension, and 8.3% each to factor
comparison and copying. Pattern matching and slot parsing do not appear in
this expansion capture. These percentages count instructions, not wall time.

A direct-output prototype emits the final pairings without generic expansion.
Five samples on CPU 6 measure 266.0 ms including destruction for length 14,
with exact equality to the production expanded polynomial. A second prototype
builds an integer polynomial in 19.4 ms but spends 248.8 ms converting it to
an Atom. These are isolated experiments, not production improvements, and
still leave roughly a fivefold gap to FORM. Caching expanded subexpressions
regresses most smaller cases and was rejected. The
#source-link("examples/reproducers/symbolica-expansion/expansion.rs", label: "standalone Symbolica expansion reproducer")
and its pinned Cargo manifest isolate this bottleneck without Spenso or Idenso.

Direct branching contraction substantially improves difficult summed words.
The baseline already includes the preceding adjacent/one-gamma reductions;
remaining gammas are compact slashes.

#table(
  columns: 4,
  [Generic-D input], [Before, ms], [Current, ms], [FORM, ms],
  [8 gammas, branching two-gamma interior], [1.142], [0.053], [0.0134],
  [Pair with three interior gammas], [2.085], [0.096], [0.018],
  [Pair with four interior gammas], [18.148], [0.412], [0.173],
  [Pair with five interior gammas], [228.236], [1.851], [2.550],
  [8 gammas, crossing pairs], [1.876], [0.042], [0.0088],
  [10 gammas, crossing pairs], [16.863], [0.192], [0.070],
  [10 gammas, nested branching pairs], [8.794], [0.112], [0.0408],
  [12 gammas, order-sensitive pairs], [160.308], [0.964], [0.763],
)

The branching-eight case improves 21.5 times, with a remaining fourfold gap
before expansion. At the preceding checkpoint it used general rewriting and
was about 55 times slower than that checkpoint's FORM measurement. Direct
identities now cover it and longer interiors. Short compact words still have
appreciable overhead: paired fourteen slashes take 11.6 microseconds against
FORM's 0.8, and alternating slashes take 36.0 against 9.4 microseconds.
A four-dimensional five-interior control improves from 0.626 to 0.527 ms.

At the shared-scan checkpoint, unchanged generic-D reruns measure:

#table(
  columns: 3,
  [Free gammas], [Idenso, ms], [FORM, ms],
  [8], [0.052], [0.031],
  [10], [0.483], [0.277],
  [12], [5.340], [3.400],
  [14], [72.235], [46.667],
)

The original fourteen-gamma rerun took about 142 ms. The shared scan first
reduced it to 81.6 ms; the final paired comparison measures 83.3 to 72.2 ms.
This last check includes cached slot metadata, uses three alternating process
rounds, and preserves every measured output exactly. The earlier 40 ms walker
prototype is a separate experiment, not the production result.

The repeated-explicit-index walk also observes bracket, dot and Dirac-operation
candidates. It writes head flags only on a hit and avoids observing an unwrapped
representation head twice. The observer now finishes this same syntactic walk
after a repeated index, with further index bookkeeping disabled, so a metric
chain can establish that bracket and dot normalization have no work. The
boolean-only predicate still stops at its first hit. Compound opaque payloads
leave later passes enabled, and any rewrite invalidates the observations.
The fourteen-gamma result still has 2.68 million tree nodes to visit.

All 39 case/dimension combinations pass 480 exact coordinate assignments,
comparing previous, current and FORM outputs against independent Clifford
multiplication in four and six dimensions. All 21 symbolic-D cases additionally
match FORM's complete polynomial exactly. Tests cover interiors two through
five, nested/crossing pairs, repeated compact slashes, free slots and scalar
spectators; eliminated summed labels are absent. The final scoped optimization
preserves every output hash from that complete certificate.

Earlier contraction checkpoints provide context: full axial-twelve
simplification improved from 2.356 to 0.621 ms, ordinary free-twelve from 6.422
to 2.794 ms, and public Schoonschip on a 64-metric chain from 738.4 to
34.7 microseconds. FORM's separately measured metric-chain CPU time was
3.2 microseconds. A single metric hit was slightly slower, 2.90 to
3.11 microseconds. These are historical measurements, separate from the
branching comparison above.

A later isolated `SlotMatcher` comparison caches the representation tag and
variance-wrapper IDs once per matcher, removing repeated Symbolica initialization
probes when many distinct tensor heads miss the small cache. Public Schoonschip
times include normalization:

#table(
  columns: 4,
  [Input], [Before, µs], [Cached metadata, µs], [FORM, µs],
  [Metric over 8 tensor terms], [11.98], [10.51], [4.0],
  [Metric over 128 tensor terms], [176.65], [141.68], [66.7],
  [Vector over 8 tensor terms], [13.74], [12.06], [4.5],
  [Vector over 128 tensor terms], [196.85], [166.19], [73.3],
)

The paired Idenso measurements use CPU 7; the separate FORM reference uses CPU 3
and times contraction plus sorting after input setup. Metric chains remain about
33 microseconds for 64 metrics at that checkpoint. Accepted slot syntax is unchanged.

The subsequent metric-component change parses endpoints while building the
incidence map and interns exact representation/dimension pairs. Since each
metric has at most two neighbors, walking from both sides directly finds a
path's boundaries or closes a loop; DFS worklists are unnecessary. Private
contractor times on CPU 7 improve from 18.85 to 14.89 microseconds for a
64-metric path, 16.56 to 12.58 for a loop and 37.53 to 30.62 for disjoint paths.
All 221 exact-output cases pass, including shuffled labels, mixed dimensions,
duality, ambiguous incidences and factored spectators. These measurements
exclude public normalization and cannot be compared directly with FORM's
complete contraction time. Callgrind records 19.3% fewer instructions;
endpoint decoding is unchanged and now accounts for 37.3% of the private cost.
The follow-up also passes all 36 Schoonschip unit tests and 46 HEP integration
tests, with Clippy clean.

The paired public API comparison, including normalization, improves 64-metric
paths from 32.86 to 28.63 microseconds, loops from 30.13 to 25.34 and disjoint
paths from 64.17 to 55.64. All 41 outputs agree exactly and the dependency
identities match; sum controls remain effectively unchanged. These measurements
use 11 alternating process pairs on CPU 7. FORM was not rerun for this change:
its earlier 3.2-microsecond path reference still leaves an approximate ninefold
gap to the public contractor.

Completing the shared observer after its first repeated index subsequently
improves public 64-metric paths from 28.09 to 23.44 microseconds, loops from
25.46 to 20.58, and metrics over 128 tensor terms from 140.60 to 124.07.
Vectors over sums remain effectively unchanged. All 61 public outputs agree
exactly across 11 alternating process pairs on CPU 7. Actual fourteen-gamma
reruns measure 72.44 to 73.53 ms, with the paired variability interval including
no change. This optimization avoids unnecessary bracket and dot passes before
metric contraction; it does not reduce the free-trace scan. The scoped checks pass
12 parser/slot tests, 36 Schoonschip tests and 46 HEP integration tests, plus
Clippy. These are separate from the earlier broader validation below.

A later cleanup experiment reuses the contractor's final traversal to detect
whether its changed output needs bracket or dot normalization. Input observations
alone are insufficient because tensor callbacks can introduce either form.
The copied-library trial matches all 249 fixtures under four settings, but is
not retained: metrics over 128 terms show a plausible gain from about 119 to
107 microseconds, while paired runs are noisy and small unchanged controls add
measurable overhead. Trace controls show no regression. The contraction record
preserves the phase profile, callback counterexamples and full timing ranges.

Isolated Symbolica patches improve conversion of the length-14 polynomial back
to an expression from 199.0 to 104.1 ms. Full polynomial expansion improves
700.8 to 605.7 ms, while ordinary expansion remains about 670 ms. A separate
flat-sort addition prototype initially regresses existing-sum inputs. Keeping
the full merge cursors outside the heap resolves those measured regressions:
adding the final length-14 monomials takes 158.2 to 35.7 ms, and adding 512
existing sums takes 595.3 to 572.5 microseconds. All 1,017 raw input/output/rerun
cases and seven upstream normalization tests pass. These are isolated primitive
measurements; none of the patches changes the production dependency. The
#source-link("examples/reproducers/symbolica-expansion/performance.typ", label: "primitive experiment notes")
record their boundaries, regressions, exact checks and paired measurements.

A separate one-pass dimension/index accessor trial does not survive the public
comparison. Despite isolated decoding gains, the 64-metric path slows from
23.24 to 24.28 microseconds, with similar 4–5% regressions in loops and gamma
reruns. All 61 outputs remain exact, but the accessor change is rejected and
the preceding observer implementation is retained.

The final retained-build checkpoint refreshes three generic-D traces and their
FORM programs. Times below are warm milliseconds; factored tracing and the full
expanded lifecycle are separate measurements.

#table(
  columns: 4,
  [Input], [Idenso, factored], [Idenso, fully expanded], [FORM, fully expanded],
  [Branching 8], [0.0506], [0.1165], [0.0132],
  [Order-sensitive 12], [0.9564], [5.9803], [0.7333],
  [Free 14], [17.4870], [798.5513], [54.0000],
)

The full-output gaps remain approximately 8.8, 8.2 and 14.8 times. Rust uses
wall time on CPU 9, including trace construction, expansion and destruction;
FORM uses its internal tracing-plus-sort CPU timer on CPU 7. Parsing and process
startup are excluded. For free length fourteen, expansion takes 753.20 ms,
about 94% of the full lifecycle. The retained build's fourteen-gamma rerun takes 73.58 ms
against FORM's 46.67 ms. Three processes with five samples each pass all nine
comparisons with historical factored outputs, and all reruns are unchanged.
The record retains cold first calls separately. Different allocator and working
set states mean the factored-only timing must not be added to an independent
expansion timing to estimate the full lifecycle.

A direct expanded constructor avoids converting the factored trace back into a
polynomial. Starting from an already classified free fourteen-gamma word, it
takes 256.79 ms with the frozen production library. A separate matched-library
experiment improves 248.23 to 157.61 ms with the polynomial-emission patches.
These timings include metric construction, polynomial assembly, Atom emission
and destruction, but exclude public gamma recognition and cleanup. They do not
replace the complete pipeline timings above.

The direct representation also has mixed downstream effects. With actual Spenso
metric normalization, the order-sensitive twelve-gamma result shrinks from
133,818 to 52,401 Atom bytes, improving its rerun from 1.927 to 0.764 ms. The
branching-eight result instead grows from 1,852 to 3,398 bytes and its rerun
slows from 29.5 to 52.6 microseconds. All thirty downstream comparison processes
agree algebraically. These cases rule out expanding every contracted word by
default; scalar spectators retain their factorization in both routes.

A subsequent public-library prototype shares the word evaluator between factored
and sparse output algebras. Its order-sensitive-twelve full expanded lifecycle
improves from 5.995 to 1.472 ms, against the separate 0.733 ms FORM reference.
However, its default first call slows from 0.974 to 1.405 ms. The rerun improves
from 1.855 to 0.725 ms. The prototype
is rejected; the retained timings above remain the production checkpoint.
These bare traces return directly from the terminal evaluator, so an outer
cleanup pass cannot explain the first-call regression.

A subsequent frozen-library profile identifies the intermediate sparse lists as
a material cost: the order-sensitive case increases from 9.34 to 14.06 million
instructions and from 3,232 to 11,609 allocations. Cloning memoized results
accounts for 2.31 million instructions; sparse finalization accounts for 6.77
million. Earlier wall-time regressions in unselected controls do not consistently
reproduce: free lengths twelve and fourteen have identical allocation counts
and instruction counts within 1%. The duplicate coefficient construction in
the interior-pair control adds only 0.03% instructions. These measurements do
not support attributing those control fluctuations to that extra construction.

An isolated follow-up replaces those intermediate lists with handles into the
existing factored-trace node storage and emits polynomial leaves once. It keeps
the recurrence and selection rules unchanged. In the matched comparison, the
order-sensitive twelve-gamma public call takes 0.993 ms with factored output,
1.483 ms with sparse lists and 0.832 ms with shared nodes; their complete expanded
lifecycles take 6.276, 1.570 and 1.099 ms respectively. Shared nodes reduce
instructions to 9.56 million and allocations to 3,026. All 52 trace cases and
87 boundary/callback files match the validated sparse-list trial exactly.
These are scratch-library measurements; production was unchanged at that checkpoint. FORM
was not rerun for this follow-up, and the earlier FORM CPU timing is not a
matched wall-clock comparison.

The retained scalar-power admission fix removes another unnecessary cleanup
pass. A power requests dot normalization only when its base is a metric function;
rank-one vector bases are already detected by the existing function observer.
The scan still visits bases and exponents, and opaque metadata retains its
conservative fallback. Scalar factors such as `(x+y)^8` no longer force
Schoonschip to walk an otherwise simplified trace polynomial.

The matched copied-production checkpoint below uses the default factored API,
with `S=(x+y)^8` kept intact. Times are warm milliseconds, including returned
output destruction, from three alternating processes with five samples each.

#table(
  columns: 5,
  [Input], [Before], [After], [Rerun before], [Rerun after],
  [S × free trace 8], [0.3758], [0.2837], [0.1452], [0.0520],
  [S × free trace 10], [2.4488], [1.6274], [1.2786], [0.4716],
  [S × free trace 12], [23.3361], [14.6916], [14.4584], [5.2234],
)

For the twelve-gamma case, instructions fall from 299.24 to 192.20 million;
the unnecessary outer Schoonschip walk disappears. Small controls remain mixed:
a scalar variable changes from 0.296 to 0.325 microseconds, while a scalar sum
changes from 0.362 to 0.347. The order-sensitive twelve-gamma source changes
from 152.6 to 158.0 ms, with its rerun changing from 3.382 to 3.360 ms.
These are workload-specific gains, not a FORM-parity result.

All 122 power cases agree across dot normalization, Schoonschip, default gamma
and experimental expanded gamma: 488 exact outputs. Eleven production inputs
and their reruns also agree. The experimental API compatibility suite passes
52 trace cases, 64 boundary/callback/scope cases, 87 unchanged default-output
files and four integration tests. Its explicit expanded-trace API was
unintegrated at that checkpoint: expanding a trace inside a larger expression still changed the
cost of later cleanup. The retained power correction is independent of that API.
The retained change passes all 40 targeted repository tests and Clippy. The unit
tests explicitly register their parsed vector namespace and check a productive power identity; the initial namespace-only
fixture failure is retained in the record. The `retained_scalar_power_admission`
entry preserves sources, samples, instruction profiles and validation logs.

The next retained change inspects the whole rebuilt expression before the
remaining cleanup passes. A complete observation with no contraction,
normalization or Dirac/epsilon work can return immediately. Otherwise cleanup
runs as before; its observation is reused on the next iteration only while every
cleanup result is exactly unchanged. This includes outer factors and callback
results, and does not assume that an unchanged spectator cannot contract with
new trace metrics.

After the scalar-power fix, a five-round comparison of the final concise
production source gives the following default-call times in milliseconds:

#table(
  columns: 3,
  [Input], [Before], [After],
  [S × free trace 8], [0.2844], [0.2064],
  [S × free trace 10], [1.6479], [0.9773],
  [S × free trace 12], [14.3991], [7.2032],
  [S × order-sensitive trace 12], [154.4985], [154.2121],
)

The twelve-free-gamma source drops from 188.59 to 96.58 million instructions.
Its rerun is nearly unchanged at 5.379 to 5.179 ms. An earlier variant added a
full scan when the certificate failed; carrying an unchanged observation
removes that duplication. All 29 final default/disabled/compute fixture modes
and 87 boundary/callback outputs match their baselines. A separate regression
checks a cleanup callback that introduces a new trace after the observation.
All 29 targeted repository Dirac simplification tests and Clippy pass; this
Clippy command does not enable warnings-as-errors.

The remaining small costs are reported explicitly: disabled-trace free-eight
and free-twelve controls increase by 1.8% and 3.6%. An untouched Symbolica scalar
expansion control changes by 1.7%, with identical instruction counts and an
unchanged rerun. These results establish a large free-trace cleanup gain,
without claiming zero overhead in every mode or FORM parity.

The expanded-trace API was unintegrated at that checkpoint. With the preceding carryover
variant on both routes, the same complete expanded free-twelve output takes
about 50.9 ms by simplifying first and expanding only the trace body, versus
62.7 ms when expansion happens inside the experimental API. Scalar spectators
remain factored in both. The `retained_post_rewrite_completion` record keeps
all three candidate stages, the final untouched-control comparison and the
callback regression separately from that experimental API.

=== Traces with exact scalar spectators

The first retained scalar-context dispatch reused the existing terminal evaluator for
one ordinary trace multiplied by exact function-free scalar factors. It keeps
$S=(x+y)^8$ intact. Matched default-call times, including output destruction,
are milliseconds:

#table(
  columns: 3,
  [Input], [Before], [Retained],
  [S × order-sensitive trace 12, D], [159.596], [2.878],
  [S × alternating slashes 12, D], [1.636], [0.0569],
  [S × paired slashes 12, D], [0.1717], [0.0163],
  [S × free trace 12, 4D], [10.995], [5.493],
  [S × free trace 12, D], [7.487], [7.243],
  [S × free trace 10, D], [0.9855], [1.0198],
)

The order-sensitive source improves 55.5 times; its rerun improves from 3.372
to 1.807 ms. This does not establish FORM parity. A separate matched-output
lifecycle simplifies the product, expands only the trace body and restores S:

#table(
  columns: 4,
  [Trace body], [Before, ms], [Retained, ms], [FORM CPU, ms],
  [Order-sensitive 12, D], [165.193], [8.104], [0.7433],
  [Free 12, D], [55.822], [53.940], [3.9333],
)

These Rust measurements use three alternating process rounds and five samples
per process on CPU 9, with native opt-level 2 libraries and debug assertions.
FORM 5.0.0 measures internal `tracen` plus sorting CPU time after setup; Rust
measures wall time. The expanded bodies still take about 11 and 14 times the
FORM reference. Factored default results and expanded FORM output are distinct
amounts of work.

Admission is deliberately restrictive: canonical ordinary words, at most two
occurrences per explicit index, exact scalar metadata and callback-free compact
vectors. Axial words, external tensor factors, compound metadata, custom vector
normalizers and numeric domains outside exact rational/complex-rational
coefficients retain the established route. This matters: a rounded trace-unit control exposed coefficient
reassociation, and an indexed-vector callback exposed a different rewrite order.
The guards preserve both original first-call and rerun results exactly.

At this checkpoint, the initial shared scan gated one-shot product admission.
The whole rebuilt product then had to pass the existing completion certificate;
otherwise that observation entered outer cleanup. The following refinement
removes this redundant output scan under the same strict admission. Neither
the experimental expanded-trace API nor the isolated Symbolica patches were
integrated by this change.

Small costs remain visible. Generic-D free-ten is 3.5% slower, while an untouched
scalar-expansion control is 1.7% slower; free-ten instructions decrease 0.73%.
The remaining wall shift has no established cause. The earlier 9.5% instruction
regression on a 128-metric product is removed by the scan gate; its final public
time changes from 28.36 to 27.77 microseconds. Heavy failed admission is unchanged
at about 2.53 ms.

Validation passes 20 public cases, 84 boundary cases, 82 exact trace-body
polynomial comparisons, 81 default/disabled/compute modes and eight exact
numeric-dimension status checks. Fifteen admitted raw products equal their
cleaned results exactly. Independent HEP networks pass 24 case/sample rows,
providing 48 comparisons against the original gamma products. The combined
retained source passes all 149 selected repository Dirac/Schoonschip tests and
scoped Clippy, without warnings-as-errors. Two older assertions were updated to
compare exact trace-body polynomials while explicitly preserving the spectator;
no production change was needed for those test corrections.
The `retained_scalar_context_trace_dispatch` entry stores source/build identities,
raw samples, independent reviews and superseded rounded/admission variants.

=== Closed terminal products and compact vectors

The strict scalar-context terminal evaluator now returns its rebuilt product
directly. Its trace recurrence emits each free explicit slot once per monomial,
and the admitted exact scalar spectators cannot add a contraction partner.
Wildcard metadata can leave a diagonal metric unevaluated, but that eliminated
dummy has no distinct partner and the diagonal is inert under existing cleanup.
The admission rules, callback exclusions and shared recurrence are unchanged.

A matched comparison of both retained changes measures these default calls
in milliseconds, against the same build before either change:

#table(
  columns: 3,
  [Input], [Before], [Combined],
  [S × order-sensitive trace 12, D], [2.907], [0.964],
  [S × free trace 10, D], [0.970], [0.569],
  [S × free trace 12, D], [7.438], [2.373],
  [S × free trace 12, 4D], [5.151], [2.974],
)

Including expansion of only the trace body gives 8.160 → 6.186 ms for the
order-sensitive word and 52.138 → 46.590 ms for free twelve. FORM references
remain 0.7433 and 3.9333 ms CPU respectively; these are not parity results.
The Rust measurements use the same native opt-level 2/debug-assertions setup,
three alternating rounds, five samples per process on CPU 8 and output destruction.
The isolated closure trial reduces contextual free-twelve instructions by 61.2%.
Its bare control moved 7.1% slower despite 1.3% fewer instructions; the combined
run instead moves about 1% slower. No cause for those wall shifts is established.
The order-sensitive rerun improves 1.819 → 0.918 ms; generic free-twelve stays
near 5.1 ms. First-call closure itself does not accelerate reruns, because the
result no longer contains the admitted trace.

The shared candidate scan also stops requesting dot normalization for canonical
compact vectors. It uses the existing strict `SlotMatcher` grammar and retains
head observation and conservative treatment of malformed or opaque shapes.
In the combined matrix, a product of 64 compact vectors takes
22.35 → 17.02 microseconds through Schoonschip and 33.49 → 18.13 microseconds
through gamma simplification. Malformed-vector fallback grows
0.435 → 0.547 microseconds; wrapped compact fallback grows 0.896 → 1.194.
A trial shape prefilter recovers those costs
but repeats decoding on successful vectors and slows those paths; it is not
retained. Paired-slash reruns grow 2.465 → 2.858 microseconds and disabled
free-twelve calls grow 8.408 → 9.191. The untouched scalar expansion control
moves 2.5% slower. This is not a uniform improvement across every input.
Three-vertex network times remain about 9.6 ms for smallest-degree order and
9.2 ms for minimum-product-terms order.

Closure validation adds 392 metadata/symbol cases and 557 wildcard cases, with
exact raw, cleaned, contextual and rerun agreement. Independent HEP networks
again pass 48 comparisons with original gamma products per build. Compact-vector
validation preserves 294 first/rerun mode records, callback counts and 42 public
benchmark outputs. The `retained_terminal_product_closure` and
`retained_compact_vector_candidates` records preserve both candidates, controls,
instruction counts and raw timings. Combined validation passes all 152 selected
repository tests, scoped Clippy, 42 compact/gamma/gluon output pairs and 23
trace lifecycle controls. The `retained_closure_compact_combination` entry
records the combined build and measurements separately from isolated trials.

=== Expanded trace output: integrated API and measured checkpoints

The source now includes the Rust and Python expanded-trace setting. It emits
polynomial output from the shared trace recurrence and existing factored recipe;
scalar spectators never enter the polynomial backend. The default remains
factored. Unsupported sparse shapes keep the established arithmetic followed by
trace-local expansion, preserving callback and metadata handling. Independent
traces keep independent expansion boundaries.

The initial matched trial, before strict bare admission and literal-four folding,
measured these complete calls with $S=(x+y)^8$ kept factored, in milliseconds:

#table(
  columns: 4,
  [Trace body], [Simplify then expand body], [Prototype], [FORM CPU],
  [Order-sensitive 12, D], [6.248], [0.851], [0.7533],
  [Free 10, D], [3.475], [1.345], [0.3433],
  [Free 12, D], [47.868], [17.153], [3.9667],
  [Free 14, D], [871.254], [418.609], [53.0],
)

This narrows the order-sensitive gap substantially, but does not establish
general FORM parity. Rust reports wall time with output destruction; FORM
reports internal trace-plus-sort CPU time with the spectator excluded. The
matched Rust matrix covers eight fixtures, bare and scalar contexts, three
alternating process rounds and five samples per mode on CPU 10. Bare free-fourteen
takes 886.7 → 369.1 ms. Expanded free-twelve reruns with S take
18.29 → 11.40 ms; faster first-call emission does not make large reruns free.

Exact comparisons pass for defaults, expanded bodies, preserved spectators,
rounded boundaries, nested callbacks and metadata, including 392 strict leaf
and 557 wildcard cases. Nine refreshed FORM references retain exact sorted-word
and archived-output checks across 27 processes. Existing HEP certificates are
referenced as prior independent evidence, not presented as fresh evaluations.

The initial prototype's broad standalone mixed-six route grows
25.27 → 34.87 microseconds because it completes possible callback output before
expansion; the strict scalar-context route avoids that work. Default scalar
free-ten also shows a noisy 24% median increase. A separate default-only probe
finds 5.082 → 5.007 million instructions, but 523 → 556 microseconds in wall
medians, with paired ratios spanning 0.978–1.222. This does not show an instruction
increase of that size; it neither explains nor eliminates the wall-time concern.
A bounded follow-up tries the same strict certificate first for bare traces
only when expanded output is requested. Mixed-six improves 34.68 → 24.53
microseconds; the order-sensitive D case stays near 0.84 ms. Exact checks still
pass, with 72 additional timed processes. Rejected metadata pays another
0.534 microseconds and callback/axial controls grow about 4%. This follow-up
is now integrated; the default wall-time concern remains recorded. The
`expanded_trace_closure_rebase_trial` entry preserves the initial full matrix
and subsequent bounded trials separately.

The integrated sparse backend also folds literal exact dimension four into its
integer coefficients, using the existing four-dimensional predicate. Rounded,
wildcard and unsupported dimensions retain their previous routes. The matched
copied-library comparison isolates this change; values below are microseconds:

#table(
  columns: 4,
  [Expanded trace body], [Before], [Exact-four path], [FORM CPU],
  [Order-sensitive 12, 4D], [432.27], [163.82], [53.33],
  [Five interior gammas, 4D], [4,086.94], [1,177.43], [—],
  [Repeated compact 8, 4D], [17.42], [13.83], [1.40],
  [Branching 8, 4D], [23.35], [21.22], [2.00],
  [Alternating 14, 4D], [20.97], [20.89], [10.20],
  [Mixed branching 8, 4D], [24.38], [24.42], [—],
)

With S preserved, order-sensitive twelve takes 443.60 → 172.34 microseconds
and the five-interior case 4.145 → 1.441 ms. The generic-D control stays
827.18 → 827.70 microseconds bare and 843.37 → 847.19 with S. Rust uses
three alternating process pairs and five samples per mode on CPU 10, timing
public evaluation, final emission and output destruction. FORM measures internal
`trace4` plus sorting CPU after setup, with S excluded; a dash means no matching
reference in this checkpoint. FORM remains faster, and mixed cases and reruns
do not uniformly improve. These measurements do not establish general parity.

All 84 timed processes, the reused trace/callback/wildcard matrix, 29 additional
numeric-domain cases and five public Rust tests pass against the frozen
candidate. The integrated build passes 158 selected Idenso tests, four installed
Python expansion tests, the pipeline check, scoped Clippy and the generated-stub
check. Notebook algebra assertions pass, but rich rendering still fails at the
Symbolica export-version boundary described above.
The checks preserve exact trace-body equality, reruns and scalar factorization;
this trial adds no fresh HEP component evaluation.

Late external metrics also contract through factored trace sums. A compatible
tensor can pass through a sum when every branch can absorb it, including metrics
inside the sum with an epsilon outside. Tagged vectors follow the same route
when vector contraction is enabled; scalar spectators retain their factorization.
Independent HEP checks cover ordinary, axial and symbolic-D traces simplified
before an external metric is attached.

An earlier broad trace/scan validation checkpoint recorded 363 passing Idenso
tests and 46 passing HEP integration tests. Three Idenso processes initially hit
the host's open-file limit; the targeted retry passed. The existing tensor-display snapshot failure remains,
with 23 tests skipped across the two suites. Both scoped Clippy checks at that
earlier checkpoint passed with warnings denied. The Spenso follow-up passes
19 relevant tests; its remaining dual-wrapper filter test fails identically with the matcher change
removed. All final measured trace outputs match the independent certificate.

The #source-link("examples/notebooks/tensor_contraction_parity.json", label: "contraction parity progress record")
retains raw samples, generated FORM programs, the generic-D benchmark driver,
source identities and validation evidence.

== Planning contractions before distribution

The gluon-rule example also exposes a separate expansion boundary. The ordered
ladder benchmark expands each new vertex into its accumulator before contracting.
A partial network can instead keep each factor opaque, discover shared outer
slots and apply metric or vector substitutions before distributing an incident
sum. Scalar spectators remain outside that local distribution. The existing
symbolic tensor, product contraction strategy and slot contractor already own
these operations; this does not require a second parser or graph abstraction.

In a closed three-vertex subcase using vertices 1, 2 and 8, the initial
expand-first measurement takes about 5.56 ms and produces 64 scalar terms. The
baseline network route with local sum distribution takes 24.7–27.3 ms; its final
polynomial emission separately takes 0.30–0.48 ms. Routes that only retain
factorization are faster but still contain uncontracted indices. Initial
depth-one partial parsing takes 0.064 ms, less than 0.3% of the network route.
The timing phases use independently warmed calls and must not be added as an
exact full lifecycle. Input string parsing and shared initial dot normalization
are excluded from these route timings.

An instruction profile identifies repeated scalar recursion as the largest
measured cost. The smallest-degree route executes 209.2 million instructions,
of which 103.8 million are inside scalar-leaf simplification. A separate trace
counts 487 network parse/merge passes and 16 contraction entries: four consume
ports and twelve rebuild plain products. All 53 captured scalar leaves have no
explicit slots; dot normalization changes none,
but numeric-coefficient distribution changes 25, so skipping all scalar cleanup
would change the output. The minimum-product-terms route similarly spends
119.5 of 230.1 million instructions inside scalar-leaf simplification. These
inclusive figures include its nested parsing and must not be added to parser
costs. The expand-first control executes 53.8 million instructions.

Source factors are also parsed again for each branch, and substitutions can
invoke multiple cleanup passes and residual-index searches. The initial parse
alone does not explain the gap. Existing order scores use
operand bytes, top-level terms and graph degree; they do not model the asymmetric
cost of distributing one source sum and substituting into its opaque target.

The retained scalar-leaf shortcut reuses the strict slot recognizer to certify
an index-free expression without parser-owned syntax. It also requires dot
normalization to leave the original expression unchanged. Such a leaf only
needs the requested numeric-coefficient distribution; all other leaves enter
the existing parser with their original expression. Brackets, trace and chain
heads, broadcasts and malformed compact slots therefore keep their existing
scope and error behavior.

At the scalar-shortcut-only checkpoint, before residual-query deduplication,
matched copied-library measurements with five alternating process pairs and
eleven samples per process give the following smallest-degree timings:

#table(
  columns: 3,
  [Input], [Before], [After],
  [One closed vertex], [0.906 ms], [0.713 ms],
  [Two closed vertices], [4.197 ms], [3.149 ms],
  [Three closed vertices], [23.696 ms], [13.619 ms],
  [Two contracted sums], [0.224 ms], [0.190 ms],
)

The three-vertex minimum-product-terms route improves from 26.631 to 14.076 ms;
the expand-first control remains 5.50–5.53 ms. Thus the network route still
has substantial overhead. Instructions decrease by 42.5–42.7%. Parse/merge
passes fall from 487 to 77 and from 607 to 117 respectively, with the same four
port-consuming contractions and twelve plain-product entries. Direct sum
contraction now accounts for 34.5% of the
smallest-degree instruction profile. These clocks include result destruction
but exclude initial normalization and final polynomial emission. The copied
libraries share compiler flags and dependencies, without reproducing Cargo's
incremental codegen layout; the host has unrelated concurrent work.

All 3,200 before/after comparisons agree, including both compile-time expansion
modes, eight ordering strategies, traversal choices and runtime sum settings.
The comparisons include 512 unchanged errors. All 39 repository Schoonschip
tests pass under process-isolated nextest, including new coefficient, scope,
malformed-input and factored-spectator cases. The scalar checkpoint and raw
measurements are retained alongside the contraction profile in the parity record.

A subsequent retained change deduplicates identical residual-index queries before
running the existing matcher, preserving its grammar and search scope. In the
three-vertex case, median residual-scan instructions fall from 12.53 to 6.33
million with smallest-degree ordering and from 12.18 to 6.18 million with
minimum-product-terms ordering, about 49% in each case. Whole-route instruction
medians fall by 4.9% and 2.9%, respectively. All 3,712 exact output/error
comparisons agree. The paired wall run is invalid for a speedup claim because
unchanged controls became strongly bimodal under concurrent host load. The
scalar-only timings above therefore remain a separate checkpoint; no additional
wall-time improvement is claimed for query deduplication.

The next retained change reuses the same certificate at depth-one network entry,
after top-level sums have been split. It skips network reconstruction for a
certified scalar term only when dot normalization is unchanged; the existing
coefficient cleanup still runs. The matched checkpoint below starts from the
retained scalar-leaf, unique-query and scalar-power changes. It uses the same
five-process-pair protocol, with smallest-degree ordering:

#table(
  columns: 3,
  [Input], [Before entry shortcut], [After entry shortcut],
  [One closed vertex], [0.698 ms], [0.440 ms],
  [Two closed vertices], [2.986 ms], [2.793 ms],
  [Three closed vertices], [13.432 ms], [10.403 ms],
  [Two contracted sums], [0.188 ms], [0.129 ms],
)

For three vertices, minimum-product-terms ordering improves from 13.695 to
9.895 ms. The unchanged expand-first medians stay within 2%; individual outliers
remain in the raw samples. Instructions decrease by 19.9% and 23.9% for the two
orders. These route clocks include result destruction and exclude input parsing,
shared initial normalization and final polynomial emission. They are a separate
matched checkpoint, not additive to earlier timings or to the initialization-gate
experiment.

The parser counts explain the gain. The earlier counts decompose as
`77 = 1 + 4*6 + 2*26` and `117 = 1 + 6 + 6 + 6 + 26 + 2*36`: one initial
network, intermediate result recursion, and final scalar-term passes. The new
route needs `19 = 1 + 3*6` passes for either order. All four port-consuming
contractions and twelve plain-product entries remain, as do the factored
outputs. The matrix has 16,960 matching output/error records plus 192 matching
caught-panic statuses for an invalid nonuniform sum boundary at deeper parsing.
It covers both expansion modes, parser depths, SinglePass, compact roots,
callbacks and scalar spectators. All 41 repository Schoonschip tests and the
scoped Clippy check pass.

The retained follow-up removes one duplicate dot-normalization pass. After a
scalar term is certified unchanged, the same `apply` owner now sends that
normalized result directly to coefficient cleanup. Parser-bound inputs still
receive both their original pre-parse normalization and final cleanup. The
matched five-process-pair checkpoint starts after the entry shortcut:

#table(
  columns: 3,
  [Input and order], [Before], [After],
  [Three vertices, smallest degree], [10.212 ms], [10.022 ms],
  [Three vertices, minimum product terms], [10.033 ms], [9.485 ms],
  [16 compact dots], [18.60 µs], [16.35 µs],
  [64 compact dots], [68.62 µs], [59.51 µs],
)

The smallest-degree wall gain is modest. One contracted-sum control is 1.8%
slower with overlapping process ranges; its instruction median still decreases.
Unchanged expand-first controls remain within 0.7%. Three-vertex instruction
medians decrease by 4.0% and 6.5%; compact-dot controls, including a factored
scalar spectator, decrease by 10–22%. All 20,032 output/error records and 192
caught-panic statuses agree across modes, depths and callback boundaries.
All 149 selected Dirac and Schoonschip tests and the scoped Clippy check pass.
These timings keep the previous route boundary and frozen dependency family;
they exclude final polynomial emission and the separate initialization-gate
experiment and do not establish FORM parity.

A sum boundary also needs a clear invariant. Fast inference uses the first
branch; the depth-limited parser does not certify matching later branches.
Generated tensor rules can validate their common boundary once. Substitution
must then respect complete slots and structural occurrences, reusing the
restricted slot contractor, and update the affected boundary. A repeated-index
scan remains a cheap candidate check, not a boundary certificate or contraction
order. Outer slots establish available contractions; estimated branch growth
and copying costs guide their order without proving global optimality.

For example, an exterior metric multiplying a sum can rename its contracted
slot in every branch while the rest of the product remains factored. When the
metric is inside a three-gluon vertex sum, exposing one vertex branch reveals
which other factor can absorb that metric. The planner can distribute that
vertex alone, substitute into the selected neighbour, and update the affected
interfaces. It need not expand every opaque factor in the product. A useful
cost estimate includes the exposed branch count, the size of the neighbour to
copy, and how many contractions become available immediately. Hidden
cancellations prevent outer structure alone from certifying an optimal order.

All 72 exact HEP component comparisons pass across four inputs, three integer
assignments and six routes. Nine FORM programs also match the final scalar
polynomials exactly. For the three-vertex subcase, FORM order 1, 2, 8 passes
through 6, 34 and 64 terms in about 328 microseconds CPU; order 8, 1, 2 passes
through 9, 45 and 64 terms in about 270 microseconds. More intermediate terms
can therefore accompany a faster run. FORM includes substitutions and sorting
from opaque vertices; the Rust wall measurements start from normalized
six-term rules. These different timing boundaries do not establish a matched
speedup ratio. The fixtures and validation are retained under
`gluon_network_planning_followup` in the contraction parity record.

=== Residual searches and Python result construction

Literal residual slots now use exact containment instead of entering the pattern
matcher. Genuine wildcard patterns retain matching semantics. In a fresh complete
scalar-output measurement, the three-vertex minimum-product-terms route improves
9.735 → 9.301 ms; smallest-degree improves 10.095 → 9.902 ms. The unchanged
expand-first control stays near 5.42 ms. Fresh FORM CPU measurements for one,
two and three ordered vertices are 5, 32.5 and 325 microseconds, with 6, 34 and
64 final terms. Rust starts from substituted normalized rules and measures wall
time including final scalar emission; FORM includes vertex substitution and
sorting. The gap remains substantial and these are different timing boundaries.

The installed Python API has a separate result-construction cost. Reusing an
unchanged normalized sum avoids repeatedly merging its growing prefixes while
preserving the previous grouping whenever a child changes. Three interleaved
process pairs, each with five samples and $S=(x+y)^8$ kept factored, measure
complete expanded-trace calls in milliseconds:

#table(
  columns: 4,
  [Trace body], [Before], [Unchanged-sum reuse], [Speedup],
  [Free 6, D], [0.827], [0.842], [0.98×],
  [Free 10, D], [117.687], [98.915], [1.19×],
  [Free 12, D], [4,203.895], [1,438.199], [2.92×],
  [Free 8, 4D], [8.204], [7.957], [1.03×],
  [Repeated 8, 4D], [0.199], [0.200], [1.00×],
  [Repeated 8, D], [1.528], [1.529], [1.00×],
)

These are release Python-host wall times, separate from the native Rust and FORM
tables. Construction of inputs, reference results and rendering are excluded;
result disposal is included. Ambient builds use other CPUs. The unconditional
n-way sum alternative was rejected because it changed established grouping.
All exact output and spectator checks pass; 43 composition/settings tests pass.

After this change, instruction profiles of free 6, 10 and 12 still attribute
95.0%, 99.0% and 99.2% of complete calls to result wrapping. Explicit multiplicity
validation accounts for 39.2%, 52.0% and 56.1%; interface reinference is another
48.0%, 40.6% and 37.2%. These are instruction shares, not additive wall-time
phases. They identify repeated validation and inference as the next boundary to
improve; speeding initial partial parsing alone cannot remove this work.

For a transformation returning the identical Atom, the Python wrapper now
allocates a fresh result with the existing interface and metadata. Changed
results retain full validation. A separate three-process-pair comparison measures
complete reruns after the sum improvement, in milliseconds:

#table(
  columns: 4,
  [Expanded body], [Revalidate unchanged result], [Reuse metadata], [Speedup],
  [Free 6, D], [0.786], [0.00738], [106×],
  [Free 10, D], [95.504], [0.676], [141×],
  [Free 12, D], [1,374.800], [9.222], [149×],
  [Free 8, 4D], [7.594], [0.0594], [128×],
  [Repeated 8, 4D], [0.165], [0.00229], [72×],
  [Repeated 8, D], [1.393], [0.0126], [111×],
)

The tensor algebra still executes on each rerun; this change removes repeated
result validation after exact equality is established. It preserves names,
scalar key arguments, unresolved ports, tensor-zero interfaces and logical port
permutations. Installed tests cover these cases and fresh object identity.
This is a rerun improvement, not a first-call speedup or a comparison with FORM.

A subsequent change shares the existing slot matcher across candidate collection
and multiplicity counting. It preserves typed index and dimension parsing,
duality, additive branch scopes and opaque metadata. Complete first-call medians
improve 0.826 → 0.711 ms for free six, 95.508 → 71.342 ms for free ten and
1,440.647 → 1,014.794 ms for free twelve in symbolic D. Free eight in 4D improves
8.461 → 5.876 ms; repeated eight improves 0.207 → 0.165 ms in 4D and
1.474 → 1.186 ms in D. All six unchanged-rerun controls stay within 1.3%.
The same three-pair/five-sample wall protocol applies; one baseline process was
slower than the others, so instruction counts provide an independent check.

At that checkpoint, complete free-six/ten/twelve instruction counts decrease by
20.6%, 27.4% and 28.0%. Multiplicity validation decreases by about half. Interface
reinference accounts for 53–61% of those calls, and validation still scans
separately for each candidate index. The change passes 44 selected
composition/settings Rust tests, six installed Python checks and scoped Clippy.

The next change collects multiplicities for all explicit indices together.
Products add counts; sums take the maximum for each index across alternatives,
without distributing the expression. Counts stop at three, the first invalid
Einstein multiplicity. The existing typed slot matcher still determines which
occurrences are structural, and opaque scalar metadata contributes no counts.
This removes the separate traversal for each candidate index.

Syntax validation and interface inference also reuse fully explicit built-in
tensor leaves within one call. The bounded cache keeps complete, unmerged
interfaces, including repeated ports and their logical order. Compact ports,
generic tensor functions and fresh-index materialization remain outside this
cache. Enclosing sums, products, powers and placeholder scopes still receive
their existing validation.

Three fresh interleaved process pairs measure these two changes together. The
baseline already includes the shared matcher and unchanged-result improvements;
all entries below are complete Python first-call wall times in milliseconds:

#table(
  columns: 4,
  [Expanded body], [Before], [After], [Speedup],
  [Free 2, D], [0.0272], [0.0262], [1.04×],
  [Free 4, D], [0.1069], [0.1023], [1.04×],
  [Free 6, D], [0.662], [0.437], [1.51×],
  [Free 10, D], [70.408], [30.332], [2.32×],
  [Free 12, D], [987.202], [401.789], [2.46×],
  [Free 8, 4D], [5.824], [2.765], [2.11×],
  [Repeated 8, 4D], [0.163], [0.133], [1.23×],
  [Repeated 8, D], [1.179], [0.700], [1.68×],
)

Each process uses five samples, retains $S=(x+y)^8$ and includes result disposal.
All eight unchanged-rerun controls remain within 2%; free twelve takes 9.07 ms.
Exact output, rerun and spectator checks pass. The count summary agrees with the
old traversal on 445,868 target queries across 2,322 fixtures, with another
36,864 representation-key comparisons. Six installed Python checks and scoped
Clippy pass. Of 90 selected Rust tests, 89 pass; the remaining wrapped-index test
also fails on the frozen baseline because its equal-index metric has already
normalized to a scalar. Its expectation is unchanged.

Instruction counts independently decrease by 35.0%, 55.3% and 59.5% for free
six, ten and twelve. In the twelve-gamma call, interface inference still takes
39.8% of instructions and tensor-power lowering takes another 18.3%. Explicit
multiplicity validation is down to 10.6%. Product normalization consumes 16.7%
across two construction sites, overlapping the constructor costs; these owner
shares must not all be added. The remaining overhead is repeated traversal and
construction around the algebra, rather than the initial network parse.

The explicit-vector candidate-scan refinement was rejected: avoiding dot walks
introduced enough additional classification work to leave the large gluon
instruction totals essentially unchanged. Source-factor parse caching likewise
has only a 0.24–0.30% measured ceiling. Normalizing newly formed nested vectors
within scoped substitution is a more promising follow-up, but cannot remove
cleanup for new repeated indices, powers or callback-produced work.

The `retained_literal_residual_index_search`, `retained_unchanged_sum_construction`,
`retained_unchanged_python_result`, `retained_shared_python_slot_matcher` and
`retained_explicit_multiplicity_and_interface_reuse`
archive entries retain validation, raw timings, profiles, controls and build
identities, including the rejected candidate.

=== Remaining complete-call gap

The following inventory combines the latest available matched-output native
checkpoints; it is not a uniform rerun of every case. Rust measures warm wall
time through output construction and disposal. FORM measures internal algebra
and sorting CPU time. Setup, parsing and process startup are excluded, and
factored scalar spectators remain outside polynomial expansion. Values are
milliseconds; the ratios describe these stated timing boundaries.

#table(
  columns: 4,
  [Native workload], [Our wall], [FORM CPU], [Observed ratio],
  [Axial 12, 4D (older checkpoint)], [0.621], [0.627], [0.99×],
  [Contracted 12, 4D, expanded], [0.164], [0.0533], [3.07×],
  [Contracted 12, D, expanded], [0.828], [0.753], [1.10×],
  [Free 10, D, expanded], [1.255], [0.343], [3.66×],
  [Free 12, D, expanded], [16.707], [3.967], [4.21×],
  [Free 14, D, expanded], [369.061], [53.000], [6.96×],
  [64-metric chain], [0.0234], [0.0032], [7.33×],
  [Metric over 128 terms], [0.1241], [0.0667], [1.86×],
  [Vector over 128 terms], [0.1651], [0.0667], [2.48×],
  [Three gluon vertices, network, complete], [9.221], [0.325], [28.4×],
  [Three gluon vertices, expand-first], [5.255], [0.325], [16.2×],
)

The gluon rows are fresh live-source measurements. FORM includes ordered vertex
substitution from opaque functions; Rust starts from substituted normalized
rules. The isolated scoped-substitution/dot-fusion candidate lowers the complete
network call to 8.734 ms, a 5.3% improvement. It agrees on 2,192 public
output/error/rerun records but is not integrated. Full cleanup remains necessary
for newly formed powers and repeated indices. No general FORM parity follows
from the nearly equal older axial-twelve checkpoint.

Fresh installed-Python measurements have a substantially larger gap:

#table(
  columns: 4,
  [Expanded body], [Complete Python wall], [FORM body CPU], [Observed ratio],
  [Free 2, D], [0.0262], [0.000640], [41×],
  [Free 4, D], [0.1023], [0.001350], [76×],
  [Free 6, D], [0.4369], [0.005400], [81×],
  [Free 10, D], [30.332], [0.353333], [86×],
  [Free 12, D], [401.789], [3.966667], [101×],
  [Free 8, 4D], [2.765], [0.048000], [58×],
  [Repeated 8, 4D], [0.1327], [0.002000], [66×],
  [Repeated 8, D], [0.7000], [0.014000], [50×],
)

All 24 fresh FORM processes pass exact polynomial/reference checks. Python
includes the wrapper, interface inference and handling of $S=(x+y)^8$; FORM
excludes S. These are user-visible boundary comparisons, not equal-stage backend
speedups. Python result construction owns 97.25% of the gated free-twelve
instruction profile. The native free-twelve expanded profile instead attributes
80.8% to polynomial-to-Atom emission and only 2.3% to the Clifford recurrence.
This identifies distinct targets: eliminate redundant validation and rebuilding
in Spynso, and improve expression emission in Symbolica.

The full eight-vertex ladder was also refreshed through the notebook's actual
Python/TensorExpression route. One complete ordered reduction takes *46.474 s*;
FORM's three-process warm wall median is *0.742 s* (0.71 s internal CPU), an
observed wall ratio of *62.7×*. All intermediate counts and every coefficient
of the 9,652-term final polynomial agree. The last two vertex stages take
13.311 and 25.147 s. This is a full-route measurement, separate from the
three-vertex network experiment; the earlier 54.970 s figure is historical.
Input construction, term counting and validation are outside the Python clock,
while FORM process wall includes startup and input parsing. The single Python
run is not a multi-run median.

A separate run of the same frozen ladder implementation adds timers at the
public operation boundaries. It takes *42.697 s*, with exact equality to the
previous final expression; these measurements are not a retroactive partition
of the 46.474 s run.

#table(
  columns: 3,
  [Ladder operation], [Wall time, s], [Share],
  [Construct TensorExpression from expanded input], [32.083], [75.14%],
  [Schoonschip, including result reconstruction], [9.846], [23.06%],
  [Multiply and initially expand], [0.700], [1.64%],
  [Normalize dots, including result reconstruction], [0.0387], [0.091%],
  [Final expansion and previous-result disposal], [0.0275], [0.064%],
  [Vertex-rule substitution], [0.0011], [0.003%],
  [Extract expression and dispose wrapper], [0.0003], [0.001%],
)

Native sampling of a second full run retains 4,289 samples within these operation
boundaries, excluding setup and validation. Within initial construction, disjoint
owners account for 64.52% of sampled user cycles in interface inference/syntax
analysis, 19.66% in tensor-power lowering and product normalization, and 15.63%
in index-multiplicity, placeholder and nesting validation. Within Schoonschip,
52.89% belongs to the Idenso operation and 45.30% to `from_transformed_atom` result
reconstruction; 1.81% is other or unresolved. These are sampling estimates, not
exact wall-time phases. The profiled run takes 44.430 s and also matches the
saved final expression exactly.

The dominant cost is reconstructing tensor metadata after distributing the
product. The output slot interface can instead be determined before expansion:
validate the small substituted factors, compose their interfaces with the
accumulator, and retain that interface through expansion and contraction. For
a contraction inside an already composed expression, its external interface
is unchanged. When composing two tensors, the surviving ports follow from the
contraction pairs, without inspecting the expanded output terms.

This requires an internal construction path that carries established invariants.
Merely using the current public tensor arithmetic or `preserving_interface`
still performs global validation or inference. The network contraction kernel
already carries its result structure but drops it when returning the Atom.
Input parsing and contraction compatibility checks remain necessary; repeated
output validation is avoidable for operations that guarantee those invariants.
Such a path must preserve logical port order, fresh dummy identities, tensor-zero
shape and the existing normalization semantics. It must not implicitly trust
arbitrary user rewrites or custom normalization callbacks. No implementation
speedup from carrying interfaces is claimed by this investigation. The `full_gluon_ladder_breakdown`
archive entry retains the driver, phases, exact comparison and build identity.

The subsequent constructor improvement reuses Spenso's cached `SlotMatcher`
for both explicit slots and compact representations. Recursive inference borrows
Atom views, shares its matcher, and reads ordinary tensor-leaf ports directly.
It no longer creates temporary dummy indices and rebuilds those leaves merely
to infer their interface. Bounded reuse of small leaf interfaces retains fresh
open-port identities; tensor-power lowering also avoids rebuilding unchanged
subtrees. Callback-bearing or order-sensitive leaves retain their normalization
semantics. Raw input still receives all-summand compatibility validation: the
existing fast inference mode's first summand cannot establish validity of an
arbitrary user-supplied sum.

Constructor-only measurements use the same saved, already parsed ladder inputs,
with one warmup and three timed calls per case. All eight output Atoms and
interfaces agree exactly with the preceding build.

#table(
  columns: 4,
  [Input terms], [Before, ms], [After, ms], [Ratio],
  [6], [0.128], [0.095], [1.35×],
  [36], [1.230], [0.771], [1.60×],
  [204], [8.213], [5.140], [1.60×],
  [1,152], [57.788], [35.181], [1.64×],
  [6,503], [419.630], [241.470], [1.74×],
  [50,684], [3,921.555], [2,241.664], [1.75×],
  [92,341], [8,151.410], [4,451.431], [1.83×],
  [186,516], [18,916.916], [9,874.551], [1.92×],
)

The complete ladder's measured operation sum improves *47.506 → 25.340 s*
(*1.88×*). Construction accounts for *36.420 → 16.941 s*, while Schoonschip,
including result reconstruction, improves *10.291 → 7.648 s*. This harness saves
fixtures between operation clocks, so its phase sum is distinct from the earlier
contiguous 46.474 s run. Fresh FORM process-wall median is *0.727 s* over three
runs (0.70 s internal CPU); the observed remaining ratio is *34.8×* and does not
establish parity. These are shared-host observations, not interleaved process
pairs. All intermediate outputs and the FORM-certified 9,652-term final Atom
match exactly. The 68-case adversarial constructor corpus, eleven metric/trace
controls, and the four-stage typed route also retain their outputs and errors.
The `retained_fast_constructor_inference` archive entry records the sources,
builds, timings and validation. Carrying a proven interface through contraction
remains a separate opportunity; this change speeds inference itself.

Separate five-call constructor medians also improve the 32-metric chain
458.7 → 378.4 µs, compact dot sum 45.2 → 24.2 µs, free-six D-dimensional trace
78.5 → 67.7 µs, and repeated-eight 4D trace 96.6 → 81.3 µs. These small-case
measurements exclude contraction; their complete output/error checks pass.

A native profile at that checkpoint takes 25.399 s in timed operations,
again with exact final equality. Construction still takes 16.952 s. Of its
1,640 samples, interface inference and syntax analysis own 44.71% of sampled
user cycles, product normalization 20.34%, explicit multiplicity/placeholder/
nesting checks 25.27%, and tensor-power lowering 9.43%. These are disjoint
sampling estimates, not wall-time subdivisions. Schoonschip takes 7.614 s:
68.31% of its sampled cycles belong to Idenso and 29.64% to result reconstruction.
The remaining constructor cost is spread across repeated global analysis and
normalization, beyond the now cheaper leaf-port extraction.

The next retained refinement avoids rebuilding unchanged products and powers:
At that checkpoint, `StructuredAtom` kept its existing Atom storage until bracket normalization
changes a child. The existing bounded inference cache also retains scalar dots
whose operands already have reusable interfaces. Syntax validation reuses checks
for directly readable leaves, with the first occurrence still visiting every
child. Callback-sensitive materialization and all validation boundaries remain.

A fresh before/after run improves the full ladder operation sum
*25.160 → 17.740 s (1.42×)*. Its constructor sum falls *16.851 → 10.155 s* and
Schoonschip including result reconstruction falls *7.547 → 6.801 s*. Fresh FORM
takes *0.719 s* process wall (three-run median; 0.69 s internal CPU), leaving an
observed *24.7×* ratio with the same phase-sum caveat. All eight intermediate
Atoms/interfaces and the final 9,652-term result agree exactly.

#table(
  columns: 4,
  [Constructor input terms], [Before, ms], [After, ms], [Ratio],
  [6,503], [251.410], [173.175], [1.45×],
  [50,684], [2,344.216], [1,589.343], [1.47×],
  [92,341], [4,651.301], [3,209.359], [1.45×],
  [186,516], [10,361.964], [7,133.692], [1.45×],
)

These are new three-call medians, separate from the preceding checkpoint.
The small 32-metric constructor improves 378.2 → 267.9 µs; the two-metric and
compact-dot controls instead move 24.6 → 25.6 and 24.2 → 25.8 µs. The cache adds
some work where an expression offers little repetition; no universal gain is
claimed. Five-sample complete trace transformations give:

#table(
  columns: 4,
  [Expanded trace], [Python before, ms], [Python after, ms], [FORM body CPU, ms],
  [Free 2, D], [0.02384], [0.02019], [0.00064],
  [Free 4, D], [0.09027], [0.07876], [0.00130],
  [Free 6, D], [0.39843], [0.34175], [0.00540],
  [Free 10, D], [29.1296], [24.0113], [0.34333],
  [Free 12, D], [387.3884], [316.3761], [4.06667],
  [Free 8, 4D], [2.60757], [2.13313], [0.04800],
  [Repeated 8, 4D], [0.11895], [0.10306], [0.00200],
  [Repeated 8, D], [0.65066], [0.54286], [0.01400],
)

Python includes its factored scalar spectator and result construction; FORM
measures trace-body expansion and sorting without that spectator. Every Python
result matches its reference and the previous build; all 24 fresh FORM processes
pass exact polynomial checks. First calls improve 1.15–1.22×, while reruns are
mostly unchanged. The `retained_unchanged_products_and_scalar_interfaces` archive
entry retains timings, sources, the 68-case constructor comparison and validation.
The targeted Rust suite passes 96 of 97 tests; the remaining wrapped-index
assertion is the previously reproduced baseline failure. Clippy and formatting
pass.

Ordinary changed, nonzero typed expansion now reuses the ordered interface after
validating the smaller factored input. Its admission check uses the existing
slot and interface owners. Callback-sensitive leaves, exposed tensor powers,
unresolved ports and other expansion modes retain output validation; unchanged
and zero results avoid the added preflight. This does not remove validation from
arbitrary transformations or assume that a supplied interface is certified.

A baseline/candidate/baseline comparison of the complete typed ladder gives
*27.740 / 20.358 / 28.684 s*. Initial expansion falls *8.687–9.250 → 1.849 s*,
a *4.70–5.00×* improvement; the complete typed route improves *1.36–1.41×*.
All eight stage Atoms and their typed logical interfaces agree, as does the
FORM-certified 9,652-term result. All 105 expansion records and 68 constructor
cases retain their results and errors. The targeted Rust suite passes 99 of 100
tests, with the same pre-existing wrapped-index failure; Clippy and formatting
pass. This is one candidate run bracketed by two baseline runs on the same CPU,
not a statistical confidence interval.

Fresh FORM takes *0.751 s* process wall (three-run median, 0.72 s internal CPU),
so the typed operation sum remains *27.1×* larger, with the phase-sum boundary
limitations above. The raw control moves 21.999 → 18.253 s despite not using
the new shortcut; that variation is not claimed as an optimization gain.
The earlier 17.740 s raw result is a separate checkpoint. Trace controls remain
approximately unchanged; no trace-kernel gain is attributed to this change.

Typed multiplication still consumes 10.603 s in the candidate run. A separate
profile of the final multiplication, with exact expanded-result checks, takes
9.136 and 8.273 s over two calls. Of 1,689 samples within those calls, port
rewriting owns 90.30% of sampled user cycles, multiplicity validation 5.00%,
and product normalization 3.93%. Sum extension and copying dominate the sampled
instructions: identity index substitutions still reconstruct the large sum.
These are sampled-cycle estimates, not exact wall-time shares. Borrowing
unchanged rewrite branches is a subsequent opportunity, with callback and
per-summand port bookkeeping preserved.

The `retained_factored_expansion_interface` archive entry includes the frozen
source, installed release, paired measurements, profile and validation.
The concurrent broader tensor API audit changed the wrapper after this release
was built; those newer source changes are preserved and are not certified by
this checkpoint's timings or checks.

An executed counterexample clarifies the next interface-retention boundary.
A valid rank-one `g(a,b)*T(a)` can invoke a custom normalizer that turns `T(b)`
into scalar one after index substitution. The current wrapper rejects the lost
tensor syntax; blindly retaining its rank-one interface would be wrong. Also,
the separately valid scalar `p(a)*q(a)` and vector `r(a)` cannot be multiplied
without accounting for the resulting three occurrences of `a`. Reusing proven
interfaces therefore needs callback and internal-index information as well as
the external slot shape. The archive retains both executed reproducers.

Measured, unintegrated Symbolica candidates include an initialization-completion
query, normalized polynomial emission, compact-heap sum merging, and lazy scratch
allocation during expansion. They have separate correctness and performance
records; their gains cannot be multiplied together. The initialization gate
improved an earlier complete gluon call by about 10–11%. Polynomial emission
improved 8.823 → 3.628 ms for twelve and 198.950 → 104.140 ms for fourteen;
those are isolated emission times. The updated
#source-link("examples/reproducers/symbolica-expansion/performance.typ", label: "Symbolica primitive experiments")
and `end_to_end_gap_checkpoint` archive entry retain the full inventory,
counterexamples, measured regressions and timing limits.

=== Shared symbolic tensor ownership

Idenso's `SymbolicTensor<S>` owns the expression together with its structure.
Its default storage remains `OrderedStructure<LibraryRep, AbstractIndex>` for
explicit symbolic networks. Positional composition uses the same type with
`PartialStructure`: this retains the canonical-to-logical layout and
occurrence-local unresolved ports. The former Spynso `StructuredAtom` has been
removed. A bare ordered structure cannot replace its logical layout; unresolved
port identifiers are metadata, not serialized abstract indices.

Composition, contraction, index substitution, chain and trace assembly, fast
structure inference, and result-interface validation now live beside that
existing Idenso owner. Spynso converts Python arguments, dispatches operations,
translates errors, and wraps results and presentation metadata. This follows
the existing dependency direction: Idenso uses Spenso's structure machinery,
and Spynso uses Idenso. Mapping a symbolic tensor's structure retains its
expression and classification flags instead of discarding them into a shell.

Known interfaces survive algebra when the operation's invariants establish
their validity. A zero retains its declared tensor shape, and unresolved ports
remain distinct and ordered. Callback-sensitive rewrites check the resulting
interface, including the metric counterexample above. External slots alone
cannot establish internal index multiplicity, so necessary index checks remain.
Result validation observes the normalized expression without materializing
temporary indices or replaying user callbacks. Constructor inference retains
its existing materialization behavior; both use the same inference machinery.

Port rewriting borrows unchanged branches, including identity substitutions.
Every sum branch still independently accounts for the requested ports; an
index missing from one summand is an error. Genuine relabeling uses Symbolica's
bulk sum and product builders when arithmetic is exact and callbacks cannot be
reordered. Rounded coefficients and callback-sensitive expressions retain their
existing evaluation order.

Unchanged results consume the already-owned result Atom while retaining the
existing structure and classification flags. Python wraps those results
directly, avoiding a temporary symbolic tensor and another expression copy.

Lorentz-dimension changes also belong to this shared owner. Changing four to
six dimensions can join explicit indices that previously occupied different
spaces, so the mapped interface must merge those pairs and reject excess
occurrences. Callback-sensitive results are checked without invoking their
normalizers again. In particular, a callback turning the changed tensor into a
nonzero scalar cannot retain the old tensor interface. Typed zeros retain their
mapped shape; Python clears descriptor metadata when contraction changes rank.
The callback check compares encoded interfaces before and after the change;
the retained logical order can differ from canonical argument storage order.

The contraction and dot-normalization walkers capture predefined symbol and
tag handles once per call. Their traversal and normalization order stay the
same, including callbacks reached through unchanged branches when another
branch changes. Skipping those branches would change observable behavior.

The final frozen release is compared with the saved pre-consolidation binary,
which already contains borrowed port rewriting. On the same pinned CPU:

#table(
  columns: (2fr, 1fr, 1fr),
  inset: 5pt,
  table.header([Operation], [Saved baseline], [Consolidated]),
  [Complete typed eight-vertex ladder], [12.792 s], [9.929 s],
  [Relabel a 1,024-term sum], [17.354 ms], [8.925 ms],
  [Identity relabel, 1,024 terms], [6.873 ms], [6.951 ms],
  [Free trace4, length 8], [2.119 ms], [1.797 ms],
  [Free tracen, length 12], [309.781 ms], [259.863 ms],
  [Axial trace4, length 12], [35.318 ms], [31.264 ms],
)

The complete typed loop is *22.4% shorter* in the adjacent final comparison.
An earlier baseline/candidate/baseline comparison gave
11.528 / 9.412 / 11.442 s, an 18% reduction. These are observed pinned-CPU
runs, not a confidence interval. The contiguous clock includes every
substitution, construction, composition, expansion, Schoonschip call and dot
normalization; imports, fixture setup, serialization and exact checks are outside
it. Every final Atom equals the FORM-certified 9,652-term result.

Fresh FORM takes *1.194 s process wall* (three-run median), leaving about an
*8.3×* full-ladder gap. Its free length-12 tracen body takes 4.033 ms CPU and
axial length-12 trace4 takes about 1.06 ms CPU. Python trace timings include a
scalar spectator and result wrapping, which FORM's trace-body clocks exclude.
The trace kernels themselves were not changed by this consolidation.

Tiny controls are mixed: free length-two first calls move 19.965 → 20.881 µs.
The final no-op follow-up restores short reruns: in its own paired comparison,
length two improves 1.825 → 1.039 µs and repeated-index trace4
3.654 → 2.312 µs. The complete loop does not establish an additional gain from
that follow-up. Raw-route observations vary: contiguous runs give
18.049 → 19.401 s at the first consolidation checkpoint and 17.442 s for the
final build, while a separate phase comparison gives 20.177 → 18.932 s.
No raw-route gain or regression is established. Its measured construction
phase still takes 11.402 s, with Schoonschip including wrapping at 6.840 s.

The final typed profile assigns 54.87% of sampled cycles to Schoonschip/dot
cleanup, 11.23% to interface inference/validation, 7.24% to explicit index
multiplicity checks, 6.07% to product normalization, and 2.63% to port rewriting.
These are disjoint sample-attribution groups, not wall-time phases. In the
preceding consolidated profile, Symbolica's initialization guard alone owned
about 10% of sampled cycles exclusively; its inclusive share overlaps the
owners above. Repeated symbol-property lookup remains an upstream optimization
target, alongside contraction and dot normalization.

Validation includes 461 enabled Idenso tests, 11 Spenso structure/slot tests,
111 of 113 binding tests plus the reflected-division regression, 34 exact
composition controls, and 117 fresh HEP component comparisons. All 68
constructor outcomes match apart from three shared-owner error-wording changes.
Eight expansion cases now correctly retain compact callback tensors: the old
errors came solely from validation inserting temporary indices. Actual
callback-induced rank loss still fails. The old strict logical-order comparison
is retained as a failure: the prior API work intentionally preserves incoming
logical order instead of re-inferring it. Separate four-stage current-order assertions and
exact algebra/component checks pass.

Clippy and formatting pass. One Idenso snapshot failure was reproduced with the
original owner and excluded; 22 existing tests remain ignored. The two binding
failures also reproduce on the baseline: wrapped-index admission and the
renderer rejecting Symbolica format 6. Rendering compatibility is not repaired
by this change. The `shared_symbolic_tensor_consolidation` archive entry retains
both measured releases, source snapshots, all strict failures, additional
certificates, profiles, and reproduction scripts. The concurrent reflected
division dispatch fix is preserved and tested; these timing drivers do not call
that method. The validated final release is installed in the notebook environment;
the installed 68-case constructor corpus matches the measured release exactly.

=== Cached operation handles and checked dimension changes

The next release captures metric symbols and tensor tags once per contraction
or dot-normalization call. It also moves Lorentz-dimension rewriting and its
interface checks into the shared symbolic tensor. The latter fixes stale rank
metadata after user normalization, including mixed representations whose logical
order differs from their encoded argument order.

Two alternating process pairs compare this release with the same saved
pre-consolidation baseline. Medians of the process measurements are:

#table(
  columns: (2fr, 1fr, 1fr),
  inset: 5pt,
  table.header([Operation], [Saved baseline], [Current]),
  [Complete typed eight-vertex ladder], [12.124 s], [8.726 s],
  [Free trace4, length 8], [2.160 ms], [1.812 ms],
  [Free tracen, length 12], [315.654 ms], [266.506 ms],
  [Axial trace4, length 12], [35.095 ms], [31.398 ms],
)

The complete ladder is about *28% shorter*. Its two paired reductions are
24.7% and 31.0%; these observations are not a confidence interval. Fresh FORM
takes *0.728 s process wall*, leaving about a *12×* full-ladder gap. The exact
9,652-term result is unchanged. FORM trace-body CPU clocks remain separate from
Python's spectator and wrapping costs.

Against the immediately preceding consolidated release, three paired complete
typed runs improve by 4.3–5.9%. Their process medians are 9.329 → 8.906 s;
raw-route medians are 17.689 → 16.915 s, with paired reductions of 1.2–4.8%.
The first typed pair is much slower for both builds (15.847 → 15.124 s), and
is retained. Axial first-call measurements in this comparison are slightly
slower, 31.098 → 31.773 ms; no additional axial pipeline gain is established.
Short reruns are approximately unchanged.

The isolated native matrix covers 91 cases and 890 timing rows. Large
`normalize_dots` examples improve by 32–35%, while complete Schoonschip calls on
large sums improve by 4–6% and long metric chains remain roughly flat. Tiny
no-work calls pay an additional 30–55 ns for capturing handles. The unchanged
phase harness fails for both binaries because it labels a newly multiplied
late metric as an already terminal gamma result. Its failed assertions and
fixture diagnosis are retained; no phase gate was weakened.

All 7,504 native snapshot comparisons, 207 exact Python behavior records, and
117 HEP component comparisons pass. The final current-source Rust run passes
466 Idenso tests and 113 of 115 binding tests, with the same two binding
failures, 22 ignored tests, and one known snapshot exclusion. Clippy and scoped
formatting pass. The `shared_tensor_handle_and_dimension_followup` archive entry
retains source and binary identities, original samples, FORM programs, callback
regressions, and the remaining profile costs.

The final complete-ladder profile assigns 53.65% of sampled cycles to
Schoonschip/dot cleanup, 11.19% to interface inference/validation, 7.03% to
explicit multiplicity, 6.19% to product normalization, and 2.77% to port
rewriting. Initialization probes still account for 9.63% exclusively;
their 24.60% inclusive share overlaps the owners above. These are sampled
cycle shares, not additive wall-time measurements. Symbolica initialization
and normalized result construction remain relevant upstream targets.

=== Reusing closed-trace interfaces

Profiling the expanded length-12 tracen call exposed a second reconstruction
cost: 94.23% of sampled cycles were spent validating and wrapping its output.
The shared symbolic tensor now reuses the Dirac evaluator's terminal-word
recognition and checks the short source interface. Only the trace's actual
cyclic wrapper is exempted from the callback check. Compact vector attributes,
unresolved ports, nested trace metadata, and user callbacks remain checked.
The proof is attempted only when output bytes exceed twice the source bytes;
small results retain their cheaper output-validation path.

Two alternating process pairs against the preceding saved release give these
medians for complete public calls, with expanded trace output and a factored
`S=(x+y)^8` spectator:

#table(
  columns: (2fr, 1fr, 1fr, 1fr),
  inset: 5pt,
  table.header([Operation], [Before], [After], [Speedup]),
  [Free trace4, length 8], [1.745 ms], [0.230 ms], [7.58×],
  [Free tracen, length 6], [0.353 ms], [0.115 ms], [3.07×],
  [Free tracen, length 10], [22.552 ms], [1.541 ms], [14.63×],
  [Free tracen, length 12], [253.761 ms], [18.836 ms], [13.47×],
  [Repeated-index tracen, length 8], [0.454 ms], [0.150 ms], [3.02×],
  [Axial trace4, length 12], [30.977 ms], [9.205 ms], [3.37×],
)

Length-2 and length-4 tracen, repeated-index trace4, and large-result reruns
remain approximately flat. The complete typed ladder is also flat:
*8.830 → 8.826 s* across three alternating pairs. The raw-expression route is
slightly slower, *17.031 → 17.302 s*, with paired increases of 0.73–2.22%; no
ladder improvement is attributed to this trace change. All runs retain the
exact FORM-certified 9,652-term result. A separate endpoint-recognition
optimization was rejected after its 91-case matrix showed no reliable gain.

Fresh FORM takes *0.753 s process wall* for the ladder, leaving about *11.7×*
between complete typed routes. FORM's trace-body CPU times are *0.0473 ms*
for free trace4 length 8, *3.967 ms* for free tracen length 12, and *0.600 ms*
for axial length 12. Those body clocks exclude the spectator and Python
wrapping; both pipelines exclude source construction from their trace clocks.
The remaining measured trace gaps are approximately 4.9×, 4.7×, and 15.3×,
respectively, across these different timing boundaries. This is not FORM parity.

The final profiles assign only 2.29% of free tracen and 2.78% of axial sampled
cycles to result wrapping; over 96% now belongs to the core trace pipeline.
These cycle shares are separate from the unprofiled timings above. They confirm
that expanded-output interface reconstruction is no longer the dominant cost.
Within the core, free tracen's polynomial-to-Atom conversion
(`to_expression_with_map`) accounts for 73.14% inclusively. This is Rust-side
Symbolica result construction. This profile uses the pinned Symbolica
`06906976` without the isolated `poly-emission.patch` and `poly-presence.patch`
experiments described above. At this checkpoint their combined effect on the
complete tracen pipeline had not been measured; the 73.14% share is not a
post-patch result.
Axial epsilon cleanup accounts for 64.23% inclusively, including nested
Schoonschip scans; the length-12 trace kernel itself accounts for 14.73%.
These inclusive shares overlap other profile categories and must not be added.
At this checkpoint the axial scalar-spectator path excluded gamma5 from its
terminal shortcut, and epsilon cleanup invoked Schoonschip before checking for
epsilon-pair work.

The final build passes 207 exact Python behavior records and 117 HEP component
comparisons. Current-worktree Rust checks pass 469 Idenso tests and 114 of 115
binding tests, with the existing wrapped-index admission failure, 22 ignored
tests, and one known Idenso snapshot exclusion. The concurrent renderer update
clears the previous rendering failure. All three new source-proof regressions,
Clippy, and scoped formatting pass. The callback regression emits a large scalar
result and verifies rank-loss rejection without replaying its normalizer.
`shared_tensor_trace_interface_followup` in the parity archive retains both
releases, samples, FORM programs, source snapshots, profiles, and the rejected
endpoint trial.

== Axial scalar contexts and corrected conversion

The next comparison starts from the preceding `43942853` build. Two changes are
measured separately against it; no combined speedup is inferred.

The retained Idenso change admits free four-dimensional axial traces in the
existing exact-scalar context shortcut. It keeps the prior expanded result,
including under default settings. One gamma5, distinct explicit symbolic
indices, supported dimensions and length, and function-free exact spectators
are required. These words emit one epsilon per term and have no contracted
metric left to simplify. Repeated indices, compact vectors, rounded coefficients
and callback-bearing spectators retain the full pass. The shared tensor owner
continues to validate callback-sensitive interfaces.

The separate Symbolica experiment applies the corrected emission and variable
presence patches to a copied `06906976` dependency. The first emission prototype
changed floating-point multiplication and callback order in its general fallback.
Two permanent regression tests now catch those errors; the corrected patch keeps
the original fallback order. Production still uses unpatched Symbolica.

#table(
  columns: 4,
  table.header([Change and input], [Before], [After], [Speedup]),
  [Idenso: axial trace4, length 12], [9.164 ms], [1.369 ms], [6.69×],
  [Symbolica experiment: tracen, length 10], [1.525 ms], [1.114 ms], [1.37×],
  [Symbolica experiment: tracen, length 12], [18.508 ms], [13.088 ms], [1.41×],
  [Symbolica experiment: tracen, length 14], [363.202 ms], [265.865 ms], [1.37×],
)

The full host timings include dispatch, result construction and wrapping, with
`(x+y)^8` kept factored. Source construction and exact checks are untimed.
Length-14 reruns remain about 184 ms. Trace4 length 8 stays near 0.23 ms.
The complete typed ladder stays near 8.8 s; one axial-comparison pair has a
17.8% slowdown despite nearly flat other pairs, so no ladder gain is claimed.
The raw route remains near 17 s; the conversion experiment observes a small
0.7% median slowdown, with paired increases of 0.35–2.90%. All outputs agree
exactly, including the ladder's 9,652-term FORM-certified result.

Fresh FORM trace-body CPU medians are 4.167 ms for tracen length 12,
53.667 ms for length 14, and 0.592 ms for axial length 12. They exclude the
scalar spectator and Python wrapping; their remaining gaps are approximately
3.1× and 5.0× for the conversion experiment and 2.3× for the retained axial
change. The ladder's fresh FORM process median is 1.111 s, versus 0.753 s at
the preceding checkpoint with the same source. This shared-host variability
leaves an approximate 8–12× ladder gap, not an improvement in our contractor.

Here expression construction means emitting the already computed polynomial
(coefficients and exponents) as Symbolica product and sum Atoms, including
canonical normalization. The expanded tracen recurrence builds that polynomial
directly, without first constructing a Symbolica sum: monomial collection is
deferred until the final coefficient list. The polynomial is an implementation
choice, not a requirement of the trace identity. Constructing the tensor wrapper
and validating its interface are separate operations.

Even with the corrected Symbolica patches, conversion accounts for 76.61% of
sampled tracen cycles: final Atom normalization takes 38.31%, other conversion
work 38.30%, and coefficient-list construction separately takes 14.33%.
In the axial profile, the trace kernel takes 76.23% and result wrapping/source
interface validation 17.22%; the previous epsilon-cleanup pass is absent.
These are sampled CPU-cycle shares, not precise wall-time phases. No wrapper
sample appears in the 290-sample conversion capture; that does not prove zero
wrapper cost.

Both candidates pass 207 exact Python records and 117 HEP component checks.
The axial regression compares 234 word/settings cases at every gamma5 position
through length 14, plus 936 scalar-decorated results and their reruns, against
the full pass. The corrected conversion passes both new regression tests,
ten previous boundaries and the Laurent case on the actual host library;
the original patch fails both new regressions. The parity archive's
`shared_tensor_axial_and_conversion_followup` entry preserves the separate
releases, exact sources, FORM programs, raw timings, profiles and validation.

== Direct Atom emission without a polynomial

A matched three-way experiment keeps the trace recurrence and its shared recipe
nodes, then replaces only final emission. Each recipe leaf builds one product
with `Atom::mul_many`; a single `Atom::add_many` collects the complete result.
It removes the coefficient-list polynomial, dense exponent arrays and
`to_expression` call. The saved baseline, corrected conversion experiment and
direct emitter all start from `43942853`, so this comparison does not include
the separately retained axial shortcut.

#table(
  columns: 4,
  table.header([Expanded tracen input], [Original polynomial, ms],
    [Patched conversion, ms], [Direct Atoms, ms]),
  [Free length 10], [1.472], [1.031], [1.605],
  [Free length 12], [18.814], [13.125], [18.866],
  [Free length 14], [379.800], [280.832], [368.768],
  [Order-sensitive contracted 12], [0.931], [0.825], [3.356],
  [Interior-pair contracted case], [3.282], [2.349], [8.926],
)

These are medians across three rotating process rounds, timing the complete
Python transformation with the scalar spectator kept factored. Exact checks,
source construction and warmup are outside the clock. Removing the polynomial
is approximately flat at free length 12 and substantially slower for both
contracted cases. The direct emitter is not retained in production.

The complete typed ladder takes 8.833, 8.804 and 9.044 seconds respectively;
the raw route takes 17.226, 16.954 and 17.044 seconds. Shared-host ranges and
paired observations are retained rather than interpreting these controls as
an improvement in contraction. Trace4 length 8 remains around 0.23 ms.
All three builds agree exactly, including the ladder's 9,652-term result.

Fresh FORM trace-body CPU medians for the two contracted examples are
0.757 and 2.550 ms. The corrected conversion experiment is competitive on these
examples (0.825 and 2.349 ms for whole Python calls), while direct emission
widens the gap. These clocks have different boundaries: FORM excludes the
scalar spectator and Python wrapping. Free-length FORM references remain the
preceding checkpoint's 4.167 ms at length 12 and 53.667 ms at length 14;
they were not refreshed in this comparison. Broader FORM parity remains open.

In the direct free-length-12 profile, product normalization owns 44.38% of
sampled cycles, other product construction 6.59%, and the final bulk sum 33.42%.
Per-leaf factor-ID sorting takes 3.77%. For the contracted order-sensitive word,
the corresponding shares are 29.69%, 8.27%, 35.20% and 6.86%. These are disjoint
owners from 288 and 290 samples, not precise wall-time phases. The polynomial
conversion frames are gone; ordinary Atom construction and normalization now
dominate. The corrected polynomial patch can skip product normalization for
proven canonical square-free factors, while the general bulk constructors still
perform it.

An untimed diagnostic confirms why contraction makes this worse: the
order-sensitive word has 2,220 recipe leaves but only 315 distinct final
monomials; the interior-pair word has 5,670 leaves and 1,890 monomials.
Polynomial coefficient collection therefore avoids roughly seven and three
times as many Atom product constructions. Free words have no duplicate
pairings: length 12 has 10,395 leaves and 10,395 monomials. All five diagnostic
outputs equal their saved full-host outputs exactly.

The isolated emitter passes 112 Rust tests, 207 exact Python records, 117 HEP
component checks and both added contracted-word checks. Its polynomial-oracle
regressions cover repeated metrics, powers, cancellation, large exact
coefficients and callback-sensitive metadata. The current worktree also passes
all 591 enabled Idenso/Spynso library tests after refreshing four tensor-power
diagram snapshots and the scoped-index return-type assertion; 22 tests remain
ignored. Clippy and formatting pass. The parity archive entry
`shared_tensor_direct_atom_emission_experiment` preserves the source, frozen
builds, benchmark drivers, all process observations and validation.

== Ladder-focused contraction profile

The current-workspace checkpoint (`f51d1774` host, unpatched Symbolica
`06906976`) runs the complete typed eight-vertex ladder in *8.551 s* and the
raw-expression route in *16.408 s*. Fresh FORM takes *0.733 s process wall*,
leaving an approximately *11.7×* gap for the typed route. These medians cover
three alternating process pairs; exact checks and fixture setup stay outside
the clock. All eight stage outputs match the saved reference, including the
FORM-certified 9,652-term final polynomial. This checkpoint also contains
unrelated workspace/dependency updates, so its small improvement over the
preceding release is not attributed to a single contraction change.

A separate unprofiled diagnostic assigns 6.009 s to Schoonschip, 1.579 s to
typed composition and 1.004 s to initial expansion. Raw construction alone
takes 10.110 s. These diagnostic phases are distinct from the contiguous
benchmark above. In the raw-constructor native profile, repeated Symbolica
initialization probes account for 27.13% of sampled cycles. Repeated accesses
to public lazy symbol bundles cause these probes even after initialization.

An untimed count explains the importance of collecting terms before dot
cleanup. At stage eight, the initial dot pass leaves 184,152 terms. One
productive contraction pass performs 197,841 metric substitutions and 299,093
rank-one substitutions, collecting 54,846 terms. Its 85,702 nested-vector
occurrences then normalize to the final 9,652 terms. There are no successful
metric-component shortcuts in this stage; the main remaining work is vector
substitution and normalized expression construction.

A subsequent directed-pair inventory finds *no opposite nested orientations*
in either ladder order. The original final stage has 26 nested-pair types, all
also present as compact `g` products; the rung-closing final stage has 26 types,
22 also present as compact `g`. Thus new nested contractions fail to collect
with existing dots until cleanup. The current vector constructor registers tags
and printing without the normalization hook that would establish that shared
form immediately. The installed frontend confirms `p(q(rep)) != q(p(rep))`
before explicit cleanup, while both clean up to the same symmetric metric.

An isolated native prototype installs that identity on default vector heads,
sharing its implementation with explicit dot cleanup. It establishes
`p(q(rep)) = q(p(rep)) = g(p(rep), q(rep))` during construction. Three
alternating process pairs measure the complete Schoonschip call on the same
preconstructed stage fixtures:

#table(
  columns: 4,
  table.header([Order and stage], [Baseline, ms], [Constructor hook, ms], [Change]),
  [Original, 6], [832.771], [702.641], [−15.6%],
  [Original, 8], [2,197.455], [2,406.311], [+9.5%],
  [Rungs closed early, 6], [58.559], [51.046], [−12.8%],
  [Rungs closed early, 8], [646.315], [653.410], [+1.1%],
)

The largest stage regresses, so this prototype is not retained. All 72 timed
outputs, expanded-output and fixed-point gates agree with the saved references;
139 focused Idenso tests and four constructor tests pass. Explicit custom
normalizers retain ownership, and automatic normalization leaves wildcard
patterns literal. Normalization now precedes enclosing scalar/metadata
construction; rounded arithmetic is checked against direct canonical-metric
construction rather than the old deferred collection order. The prototype
also needs an owner-controlled registration boundary before integration.
These native results exclude frontend construction and establish no full-ladder
or FORM improvement. The archive preserves both the measured source and the
test-only follow-up under `native_vector_constructor_normalization_experiment`.

Sampling that same prototype explains the regression: constructor normalization
owns 15.6% of the original-order stage-eight samples, while subsequent dot
cleanup still owns 30.9% and the remaining slot contraction owns 43.5%.
In the rung-closing order those shares are 12.5%, 32.2% and 44.9%.
They are disjoint sampled costs, not stopwatch phases. Canonicalizing earlier
removes one representation mismatch, but still repeatedly constructs dots before
equal scalar products combine. The `native_vector_constructor_profile` archive
preserves the exact measured libraries and decoded stacks.

A subsequent isolated prototype caches successful canonical nested-vector
identities in a bounded thread-local cache. Three alternating process pairs
reduce original-order stage eight from *2.172 to 2.066 s* and rung-closing
stage eight from *651 to 572 ms*. Stage six improves by 15.2% and 16.5%,
respectively. All 36 timed outputs and the callback/fixed-point gates agree;
139 Idenso tests and six constructor/cache tests pass. This remains an isolated
native experiment: owner-controlled registration and complete-host validation
are required before integration. `native_vector_constructor_cache_experiment`
preserves both the first screen and the refined cache's separate measurements.

The underlying symmetric-function constructor also has avoidable work: it
copies and rebuilds arguments even when already ordered. An isolated Symbolica
patch retains ordered arguments and sorts borrowed views otherwise. A scalar-dot
shaped two-argument microbenchmark improves from 359 to 227 ns when ordered,
and 357 to 304 ns when reversed. These primitive clocks establish no ladder
gain. The standalone `symmetric_construction` binary and
`patches/symmetric-construction.patch` in
`examples/reproducers/symbolica-expansion` preserve the reproducer and patch;
`symbolica_symmetric_borrowed_arguments_experiment` retains exact differential
gates, source identities and clocks.

Emitting compact dots during contraction passes exact and callback checks but
is not retained: stage six improves 832.73 → 683.74 ms, while stage eight
regresses 2,192.48 → 2,243.86 ms in three paired native runs. Matched profiles
show function construction moving from cleanup into the contractor. Doing it
before terms combine repeats work. These native stage timings exclude Python
wrapping and do not establish a full-ladder gain.

The parity archive entry `current_ladder_profile_checkpoint` preserves the
full-host baseline, exact and HEP checks, FORM source, raw clocks and sampled
profiles. Further changes must be compared against that same frozen host.

The retained constructor changes reuse initialized symbol-bundle handles within
each syntax predicate and consult the existing successful-interface cache before
repeating the syntax walk. The cache still stores shapes only: callback-aware
result validation, logical port merging and index checks remain in place.
Against the same frozen host, three paired full runs give:

#table(
  columns: 3,
  table.header([Complete ladder route], [Baseline], [Handles and shape reuse]),
  [Typed], [8.489 s], [7.936 s],
  [Raw expression], [16.798 s], [14.802 s],
)

Every typed pair improves; its median falls 6.5%. The raw median falls 11.9%,
but the first raw pair regresses 20.4%, so its improvement is less consistent.
Fresh FORM takes 0.733 s process wall, leaving a *10.8×* typed gap. Separate
phase diagnostics reduce raw construction from 10.128 to 7.431 s and typed
composition from 1.553 to 1.264 s. These are individual diagnostic runs, not
additive estimates of the medians above. Trace4/tracen controls show small
mixed changes; no trace improvement is claimed.

All 207 exact comparisons and 117 HEP component checks pass. The live Idenso
and Spynso libraries pass 592 tests, with 22 default ignored; Clippy and
formatting pass with no source drift during the checks. The archive entry
`retained_ladder_constructor_handles_and_shape_cache` records the isolated
two-file source change, all paired observations, FORM runs and the rejected
early-dot experiment.

Closing each rung earlier reduces the frontier from rank five to rank three.
On the frozen baseline, changing only the order gives the following complete
loop/process medians:

#table(
  columns: 4,
  table.header([Vertex order], [Idenso], [FORM], [Ratio]),
  [1, 2, 3, 4, 5, 6, 7, 8], [8.553 s], [0.729 s], [11.7×],
  [1, 2, 8, 3, 7, 4, 6, 5], [2.593 s], [0.178 s], [14.5×],
  [5, 4, 6, 3, 7, 2, 8, 1], [2.569 s], [0.173 s], [14.9×],
)

Both engines use the same selected order and retain the same final polynomial.
FORM improves more, so a poor ordering does not explain the relative gap.
Python accumulates applied vertices while FORM carries the unapplied vertex
functions throughout; intermediate term counts need not agree. The reverse
rung-closing order reduces Python's largest expanded input from 186,516 to
60,314 terms. These are two explicit alternative orders, not measurements of
an automatic planner or a claim of global optimality.

A further three-pair comparison on that reverse order isolates the retained
constructor changes: *2.579 → 2.297 s*, with every pair improving. Fresh FORM
takes *0.171 s*, leaving a *13.4×* gap. The notebook now carries typed interfaces
between stages and defaults to this rung-closing order; its selector also
retains the original order and applies the choice to FORM. Both complete
notebook routes check their 9,652-term result against FORM's exported polynomial.
The validated build is installed in the notebook interpreter; saved baseline
environments remain unchanged.

The preceding validation build also fixes a pre-existing strict-syntax inconsistency:
`dind` accepts representation slots without treating an arbitrary wrapped tensor
as a slot. Its complete library suite passes *873 tests*, with 22 default ignored,
plus all-target Clippy and formatting. The 207 exact and 117 HEP checks and both
actual notebook orders pass again. Fresh medians are *8.224 s* in the original
order and *2.490 s* closing rungs early, versus FORM's *0.748 s* and *0.182 s*:
remaining gaps of *11.0×* and *13.7×*. Paired timing changes from the preceding
build are mixed, so this correctness fix is not claimed as a speedup. The entry
`final_strict_wrapper_ladder_validation` records this final build separately
from the isolated performance comparisons above.

The shared `SlotContraction` now collects closed metric/vector components before
constructing intermediate scalar products. Each monomial supplies a small graph
of explicit indices: metrics join nodes, vector occurrences terminate paths,
closed paths yield canonical dots, and metric cycles yield dimensions. Existing
dots and newly closed paths share one variable table. Exact rational
coefficients are collected through Symbolica's coefficient-list polynomial
constructor, followed by one expression emission; no general
expression-to-polynomial conversion is required.

Admission requires plain vectors and canonical self-dual representations, checks
each component's index multiplicity, and preserves the existing local pairing
of positive powers. Unknown functions, metadata-bearing vectors, custom
callbacks, unsupported powers and open components retain the established
normalization schedule. A checked 64 MiB bound limits the dense exponent
storage. Failed admission constructs no speculative callback-bearing result.
The early attempt uses the existing repeated-explicit-index flag, so scalar
reruns avoid another collector walk.

Six native processes compare the selected implementation with the frozen
baseline and an unguarded intermediate. Original-order stage eight improves
*2.189 → 0.832 s*, and rung-closing stage eight improves *646 → 241 ms*;
scalar reruns remain approximately *13.5 ms*. The open stage-six controls are
1–4.3% slower in this screen, and tiny rejected inputs retain sub-microsecond
admission overhead. These are complete native Schoonschip calls on saved stage
inputs, not complete-host or FORM timings. The archive entry
`guarded_pre_dot_closed_component_collection_native_screen` retains all clocks,
exact/callback checks and the proof that the retained code differs only by
removing an experiment-specific fixture test.

The selected collector's complete-host validation uses three alternating process
pairs against a frozen baseline with the same current frontend and dependency
snapshot. Each complete loop checks the full 9,652-term FORM-certified polynomial
outside its clock:

#table(
  columns: 4,
  table.header([Route and order], [Baseline, s], [Selected, s], [Median change]),
  [Typed, original], [7.832], [6.979], [−10.9%],
  [Typed, rungs closed early], [2.349], [2.160], [−8.0%],
  [Raw, original], [13.784], [12.500], [−9.3%],
  [Raw, rungs closed early], [4.020], [3.554], [−11.6%],
)

All original-order typed pairs and all raw pairs improve. Rung-closing typed
pairs range from 13.7% faster to 3.1% slower, so their median gain has more
variability. Fresh FORM process medians are *0.725 s* and *0.174 s* in the same
two orders, leaving typed gaps of *9.6×* and *12.4×*. Python clocks include all
eight reduction steps but exclude imports, fixture setup and final checks;
FORM clocks include process startup. These results establish improvement,
not FORM parity.

Separate instrumented runs attribute the remaining operation time as follows:

#table(
  columns: 3,
  table.header([Public operation], [Original order, s], [Rungs closed early, s]),
  [Schoonschip], [4.389 (66.6%)], [1.214 (62.2%)],
  [Typed composition], [1.243], [0.467],
  [Initial expansion], [0.933], [0.256],
)

Public intervals include applicable result wrapping, validation and old-value
disposal. These clocks do not separately resolve those inner owners and are not
additive estimates of the complete-loop medians. The remaining phases account
for less than 1% in each diagnostic.

Trace controls show mixed changes: free length-12 `tracen` takes
18.777 → 18.707 ms, free length-8 `trace4` 0.229 → 0.226 ms, and free length-10
`tracen` 1.481 → 1.551 ms (+4.7%). Axial length-12 takes 1.560 → 1.358 ms.
No general trace improvement is attributed to this collector; the archive
preserves all first-pass and rerun samples.

All 877 library tests, 207 exact comparisons, 117 HEP component checks and both
actual notebook orders pass, with the 22 existing ignored tests unchanged.
All-target Clippy and formatting pass. The final library source matches the
measured host; a benchmark-only lint expectation is recorded separately.
The validated collector is installed in the notebook interpreter. The vector
constructor hook/cache and the isolated Symbolica patch remain experiments.
`retained_closed_component_ladder_validation` preserves the complete samples,
FORM programs, source identities and validation records.

A subsequent profile of that validated build separates native contraction from
shared interface proof, validation and wrapping. Within public Schoonschip
intervals, those groups account for 76.3% / 21.2% in original order and
66.3% / 32.0% when closing rungs early; the remaining 2.6% / 1.8% is unresolved.
These are sampled cycle shares (427 / 117 samples), not additive wall-time
measurements. The smaller reverse-order sample is directional evidence.
`current_closed_collector_ladder_profile` records the phase boundaries and
classification, with no samples attributed to both groups.

FORM's retained output counts agree with this build at all eight stage
boundaries in both orders. Its generated-term counter measures a different
boundary from the already combined expanded input passed to Schoonschip:

#table(
  columns: 4,
  table.header([Original stage], [FORM generated], [Expanded input], [Both retain]),
  [6], [72,828], [50,684], [10,947],
  [7], [131,364], [92,341], [23,937],
  [8], [287,244], [186,516], [9,652],
)

The unapplied FORM vertices are opaque multiplicative spectators, so they do
not add terms at these boundaries. Equal counts alone do not establish equality
of every intermediate expression; the complete final polynomial is checked
independently. FORM applies operations term by term and sorts generated terms
in buffers, as described in chapters 4 and 16 of its
#link("https://www.nikhef.nl/~form/maindir/documentation/reference/man.pdf")[reference
manual]. Its generated count is therefore not a peak number of simultaneously
materialized terms. Our expanded-input count describes a complete Atom sum.
The counters do not quantify the cost of that materialization or of postponing
contractions until after expansion.

The collector now also reduces open components: a vector-to-free-index path
becomes an indexed vector, and a path between two free indices becomes a metric.
It retains the exact external slot syntax; the existing shared symbolic tensor
owns logical order, metadata and typed zeros. Callback-sensitive or unsupported
expressions keep the established checked route. This extends the same collector,
renamed from `closed` to `components`, without introducing another tensor type.
The shared callback certificate also reuses its metric symbol handle within one
walk, preserving lazy initialization and every metadata/callback check.

Three fresh alternating process pairs compare these changes with the preceding
closed-component build, using identical frontend and dependency snapshots:

#table(
  columns: 4,
  table.header([Route and order], [Before, s], [After, s], [Median change]),
  [Typed, original], [6.547], [4.784], [−26.9%],
  [Typed, rungs closed early], [1.880], [1.502], [−20.1%],
  [Raw, original], [12.535], [10.684], [−14.8%],
  [Raw, rungs closed early], [3.574], [3.445], [−3.6%],
)

Every paired ladder comparison improves, with all complete 9,652-term results
checked exactly. Fresh FORM process medians are 0.724 / 0.177 s, leaving typed
gaps of 6.6× / 8.5×. The timing boundaries remain those described above; these
gains must not be added to the isolated native collector or certificate gains.
All 880 library tests pass (22 existing ignored), as do all-target Clippy,
formatting, 207 exact comparisons, 117 HEP component checks, and 30 additional
interface cases with 46 HEP component checks. The added cases cover open paths,
logical order, metadata, typed zero, unresolved ports, missing branch indices
and callbacks that change rank.

Separate diagnostic runs of the new build spend 2.484 / 0.823 s in public
Schoonschip, 1.241 / 0.442 s in composition and 0.925 / 0.257 s in initial
expansion (original / rung-closing order). Thus contraction still accounts for
about 53% of these runs, while composition and initial expansion together account
for about 46%. These are measured complete operation intervals, not a breakdown
of the paired medians or of their improvement. Trace controls show no general
gain: free length-12 `tracen` first calls change 19.032 → 19.615 ms (+3.1%),
while their reruns remain 9.207 → 9.210 ms; other first-call trace/axial controls
change between −0.1% and +0.7% in this checkpoint.

The supplied pure Symbolica recipe exposes a larger scheduling difference. It
keeps all unapplied vertices opaque, substitutes contracted indices into their
arguments, and expands and contracts each replacement's right-hand side before
multiplication into the surrounding expression. The existing replacement cache
can then reuse the result for identical bound vertices. The notebook includes
this explicit alternative using plain symmetric linear `d` and tagged indices,
derived from the same routing and vertex rule as the typed example. The selected
order applies to both, and their complete final polynomials are compared after
notation conversion, including the existing FORM certificate.

A controlled ablation keeps that same plain `d` notation and contraction kernel
while changing the schedule. The normalized expanded sum at the final vertex
has the following sizes, before the next separate global contraction, if any:

#table(
  columns: 3,
  table.header([Schedule], [Original order], [Rungs closed early]),
  [Accumulator, expand before contracting], [186,516], [60,314],
  [Bind future vertices; contract after expansion], [31,724], [22,148],
  [Bind future vertices; contract the held RHS locally], [9,652], [9,652],
)

Across all eight steps, the largest expanded sums are 186,516 / 60,314 terms
for the accumulator and 36,281 / 20,570 for local outside-in contraction.
All six schedule/order combinations, with the replacement cache both enabled
and disabled, give exactly the same final 9,652-term polynomial. The local and
global outside-in variants have identical counts after each completed step;
their temporary expanded sums differ substantially. Thus completed-stage counts
alone missed a material difference in the amount of intermediate algebra.
These counts describe materialized normalized sums, not FORM's generated-term
stream. The recipe is specific to this ladder; it does not replace the generic
typed API's interface validation or callback semantics.

Three fresh processes per configuration, with term counting and final checks
outside the ablation clocks, give the following medians with cache size 1000:

#table(
  columns: 3,
  table.header([Scalar Symbolica schedule], [Original, s], [Rungs closed early, s]),
  [Expand-first accumulator], [2.573], [0.677],
  [Bind future vertices; global contraction], [1.438], [0.584],
  [Bind future vertices; local held contraction], [0.867], [0.320],
)

The supplied loop itself, with only the pinned API's equivalent level-argument
spelling changed, takes *0.874 / 0.324 s*. Its clocks retain the supplied
per-stage count/print operations and exclude setup. Fresh FORM process medians
are *0.724 / 0.173 s*, including process startup: gaps of *1.21× / 1.87×* for
this scalar recipe, not parity of the generic tensor API. Every final polynomial
agrees exactly. The first ablation round is systematically slower; all raw
observations are retained, and these median comparisons must not be read as
precise confidence intervals. Later local/global observations and all literal
user/FORM observations are stable within their recorded ranges.

Disabling the replacement cache changes local outside-in medians to
2.477 / 0.850 s. Bound vertices repeat across the surrounding sum, so caching
their locally reduced replacement matters. The accumulator has only one
vertex match per step; its initially different cache-on/off medians are
temporally confounded and do not establish a caching benefit there.
A bounded follow-up alternates adjacent cache-on/off accumulator pairs: disabling
the cache changes time by −1.0% to +2.7%, confirming no substantial benefit for
that schedule. The original observations remain archived unchanged.
The plain-metric accumulator also changes representation, contraction kernel
and interface work relative to the typed route; its residual timing difference
cannot be assigned entirely to validation. The combined implementation,
notebook checks, frozen sources and complete measurements are recorded under
`retained_open_component_ladder_validation`.

A bounded cache of successful nested-vector-to-dot rewrites also passed its
52 focused tests and exact/callback checks, but is not retained: stage six
regresses 3.0% while stage eight improves only 2.1%. Reusing those local results
does not remove global product and sum normalization. The entries
`rejected_ladder_nested_dot_memoization` and `ladder_order_comparison` preserve
the rejected source and the order comparison respectively. The next structural
experiment can reuse Spenso's existing interface-based pair selection; shape
alone cannot distinguish the two mirrored rung-closing orders, which have
identical frontier ranks but different term counts.

== Applying the schedule to native tensor contractions

The same experiment now uses Spenso's intrinsic metric, tagged vectors and
existing Rust Schoonschip implementation. All future vertices remain in the
expression. A held replacement expands and contracts its substituted body;
the ambient contraction then binds exposed indices into future vertices.
The routing and six-term vertex rule are exported from the notebook fixture.

Three counterbalanced fresh processes per configuration give:

#table(
  columns: 3,
  table.header([Native Atom schedule], [Original, s], [Rungs closed early, s]),
  [Expand-first accumulator], [2.531], [0.675],
  [Accumulator with local RHS contraction], [2.536], [0.673],
  [Bind future vertices; global contraction], [1.580], [0.539],
  [Bind future vertices; local held contraction], [1.376], [0.481],
)

The outside-in reduction is 45.6% and 28.7% shorter than the corresponding
accumulator. Adding local contraction to the accumulator alone does not help.
A separate paired comparison gives 1.377 / 0.469 s with replacement caching,
versus 2.829 / 0.948 s without it. All 36 measured final results equal the full
9,652-term FORM-certified polynomial and are unchanged by another Schoonschip
pass. Four closed two- and three-vertex subgraphs also pass independent FORM
polynomial and exact component comparisons.

These are native Atom clocks on the saved tensor implementation, not Python
tensor-interface clocks. They include all eight replacement/expansion/contraction
steps and ordinary intermediate lifetimes, and exclude setup, input parsing,
term counting, final checks and final-result destruction. The ambient pass is
full Schoonschip, broader than the supplied recipe's targeted vertex absorption.
CPU affinity is fixed; unrelated jobs on other cores mean the host is not
globally quiet. Raw samples and source identities are retained.

Separate diagnostic phases take 0.208 s for held replacement, 0.575 s for
ambient expansion and 0.627 s for ambient Schoonschip in the original order.
The local Schoonschip kernels account for only 2.25 ms across 295 cached RHS
evaluations. Within the ambient Schoonschip profile, approximately 53% of
samples belong to slot contraction, 28% to dot normalization and 18% to
candidate scanning. These sampled shares are not additive wall-time phases.

The notebook also executes this schedule through `TensorExpression`, retaining
its interface checks. Routing momenta are explicitly opaque scalar metadata;
the final three vertex arguments are ports. Direct compact tagged vectors in
those positions consume ports without creating external indices. Their remaining
logical port order is retained, and the result has no atomic stored-data
descriptor. Scalar wrappers retain their metadata opacity; arbitrary tensor
metadata and malformed vector arguments remain errors.

The shared inference implementation recognizes these consumed ports directly.
Scalar algebra can retain an established interface for plain, normalized leaves,
including leaves with opaque routing metadata. This proof does not change how
constructors handle callbacks. Callback-sensitive contraction results are still
observed and checked: a normalizer that turns `T(b)` into a scalar after
contracting `g(a,b)*T(a)` must not leave a stale rank-one interface. Network
validation also retains encoded open-index owners, while partial interfaces
continue to describe unresolved ports by logical position.

Profiling the first typed implementation exposed repeated successful interface
proofs across terms: they occupied about 35% of operation-filtered samples in
both orders. The inference owner now remembers successful normalized function
leaves within the operation, using its existing limits of 256 entries and
256 bytes per key. This cache stores no inferred port identities and admits no
failed proof. Constructor inference and callback-sensitive validation retain
their separate semantics. Index-occurrence collection also uses an equivalent
shorter head predicate, avoiding redundant composite classification without
changing which indices are counted.

The completed stages are still not identical to FORM's. In the original order,
the typed outside-in route retains 11,959 and 25,662 terms after vertices six
and seven, compared with FORM's 10,947 and 23,937. In the rung-closing order,
the penultimate stage retains 11,836 versus 10,516. These counts include opaque
future vertices; every route finishes at the same 9,652-term polynomial.
Reducing temporary expanded sums therefore does not establish equality of the
intermediate algorithms or eliminate all remaining algebra work.

== Outside-in schedule checkpoint

Three counterbalanced fresh processes per route and order measure the complete
eight-step loop, retaining ordinary intermediate lifetimes. Input setup, counting
and final checks are outside the clock. The preserved baseline and final host
use the same pinned Symbolica source; the comparison includes the shared tensor
changes and the selected schedule.

#table(
  columns: 3,
  table.header([Route], [Original, s], [Rungs closed early, s]),
  [Previous typed accumulator], [4.827], [1.528],
  [Final typed accumulator], [4.470], [1.363],
  [Typed outside-in before proof caching and shorter index checks], [4.964], [1.676],
  [Final typed outside-in], [2.909], [1.044],
  [Supplied scalar Symbolica recipe], [0.897], [0.330],
  [FORM, complete process wall], [0.723], [0.175],
)

The final typed route is *39.7% / 31.7% faster* than the preserved accumulator,
and *41.4% / 37.7% faster* than the first typed implementation of the same
schedule. All 24 typed measurements and six fresh supplied-recipe measurements
equal the full FORM-certified polynomial. FORM remains *4.02× / 5.96× faster*
than the typed route under these clock boundaries; the supplied scalar recipe
is *1.24× / 1.88×* slower than FORM. This is not general contraction parity.
The shared host has visible variation: final original-order samples range from
2.906 to 3.128 s, and early-rung samples from 1.039 to 1.065 s. These are three
process medians, not confidence intervals.

The final operation-filtered profile reduces successful proof scanning to about
1–2% of samples. Arbitrary replacement's result observation and index checks
now account for approximately 52% of the complete typed loop, using disjoint
native parent categories. This is the largest remaining typed overhead.
Within expansion, ordinary Symbolica normalization accounts for about 68% of
samples. Profile shares are diagnostic and cannot be added to the benchmark
wall times as separately measured phases.

Trace controls retain the same factored scalar spectator and requested expanded
output before and after. Across nine first-call cases, the final/baseline
differences range from −2.1% to +0.6%; there is no substantial trace improvement
in this update. Representative fresh comparisons are in milliseconds:

#table(
  columns: 4,
  table.header([Case], [Previous Python], [Final Python], [FORM trace + sort CPU]),
  [Free length-8 trace4], [0.228], [0.229], [0.0480],
  [Repeated-index length-8 trace4], [0.0944], [0.0929], [0.0020],
  [Free length-12 tracen], [19.581], [19.493], [3.8667],
  [Repeated-index length-8 tracen], [0.1496], [0.1471], [0.0136],
  [Axial length-12 trace4], [1.376], [1.374], [0.6000],
)

FORM trace clocks amortize internal tracing and sorting over repeated copies;
they exclude the Python boundary and scalar spectator. Both first calls and
unchanged-result reruns are retained in the record. Exact FORM polynomial
certificates, 207 unchanged Python behavior records, 117 exact HEP component
rows, 26 new public-boundary controls and two complete factored-ladder component
assignments pass. The live worktree passes 897 Rust tests, with the previous
22 skips retained, Clippy and formatting. The notebook executes all three routes
in both orders against fresh FORM. Raw clocks, source identities, failed
intermediate diagnostics and final validation are recorded under
`retained_outside_in_typed_ladder_validation` in the
#source-link("examples/notebooks/tensor_contraction_parity.json", label: "contraction record").

== Shared replacement validation

Arbitrary replacement must check the resulting tensor interface: a normalizer
can remove a port, and a replacement can introduce a dummy index that collides
with another factor. The shared `SymbolicTensor` now observes the result and
finishes it through one checked boundary. Tensor scope checks performed during
observation are not repeated during construction. The Python wrapper retains
descriptor metadata and wraps the checked value without constructing it again.

The explicit-index collector also reuses bounded summaries of repeated,
normalized function nodes, including scalar dots with no explicit indices.
Each invocation owns its cache. Products add counts, sums take the maximum
over their branches, and encoded open-index owners remain distinct. This walk
inspects syntax and does not execute normalization callbacks. Logical port
order, typed zeros, unresolved ports, and existing error precedence are
unchanged; callback-sensitive rank loss still fails validation.

Three counterbalanced fresh processes per route and vertex order compare the
saved outside-in checkpoint with these validation changes, using the same
notebook, dependencies and contraction schedules. Medians of the complete
eight-step loops are in seconds:

#table(
  columns: 5,
  table.header([Route], [Original, before], [Original, after],
    [Early rungs, before], [Early rungs, after]),
  [Typed accumulator], [4.448], [3.865], [1.390], [1.136],
  [Typed outside-in], [2.935], [2.220], [1.066], [0.846],
)

The outside-in route improves by *24.4% / 20.6%* and the accumulator by
*13.1% / 18.2%*. All 24 runs retain exact equality to the 9652-term FORM
polynomial, scalar rank and an unchanged contraction rerun. Counting and
equality checks remain outside the clock. The early-rung candidate samples
include 0.792, 0.846 and 0.847 s; the median retains that variation rather
than selecting the fastest run.

Fresh complete FORM process medians are 0.749 / 0.180 s; the supplied scalar
Symbolica recipe takes 0.887 / 0.330 s. The typed outside-in route is therefore
still *2.97× / 4.71× slower than FORM*. These clock boundaries retain the
distinction between the typed loop and FORM process startup. The changes reduce
validation work without changing intermediate term counts or the algebraic
schedule.

In separate operation-filtered profiles, explicit-index counting falls from
about 19% to 3% of sampled cycles. Result validation and finishing together
fall from approximately 53% / 52% to 38% / 34%. The main remaining categories
are interface observation and scope validation, Schoonschip contraction, and
Symbolica normalization. These are sampled CPU shares, not independently
measured wall-time phases. In particular, removing duplicate validation does
not remove the obligation to inspect arbitrary replacement results.

All 901 Rust tests pass, with the previous 22 skips retained, along with Clippy
and formatting. The unchanged 207 API controls, 117 HEP component rows,
26 public metadata/callback controls, all three notebook routes in both
orders, and the full factored ladder at two exact component assignments pass.

Trace controls do not show a uniform improvement. Free length-12 tracen changes
from 20.434 to 20.798 ms in the paired process medians. Axial length-12 shows a
regression: the original short-loop control is 1.386 to 1.554 ms, and a separate
longer-loop check confirms 1.392 to 1.580 ms. Its unchanged-result rerun also
increases, from 2.537 to 2.902 ms in that follow-up. The rerun takes the existing
identity shortcut and bypasses the new validator; the cause is not established.
The retained change prioritizes the ladder's reduction of hundreds of
milliseconds, while preserving this trace regression in the measurement record.

Source identities, complete clocks, profiles, checks and the axial follow-up
are recorded under `retained_replacement_validation_consolidation` in the
#source-link("examples/notebooks/tensor_contraction_parity.json", label: "contraction record").

== Tensor-factor replacement

`TensorExpression.replace_tensor(pattern, rhs, rhs_cache_size=1000)` applies
tensor identities to whole factors inside sums and products. It uses
Symbolica's matcher and supports held right-hand sides, including the local
expansion and contraction used by the ladder. Tensor arguments and scalar
metadata are opaque to this traversal.

The operation checks the actual matched factor and the reduced right-hand side
for compatible explicit ports. A replacement may consume internal contracted
indices, but cannot introduce extra occurrences or unrelated dummy indices.
These restrictions let it preserve the surrounding tensor's established
interface. This certifies tensor structure; the supplied rule still determines
the algebraic identity. Repeated right-hand sides and their local summaries share a bounded
cache; conditions still run at every match. The left-hand side is not rebuilt
for validation, so its normalizer is not replayed. A source-admission walk is
still required; the saving is avoiding inference of the complete result.

Unresolved ports, remaining user normalization or evaluation hooks, and tensor
powers that the interface proof cannot certify are rejected. RHS callbacks may
produce an admitted normalized expression; set `rhs_cache_size=0` when those
callbacks have side effects. Use `replace` for general Symbolica traversal and
its full-result interface checks, and `rename_indices` or `reindex` for deliberate
changes of index labels. Tensor-factor replacement retains logical port order and typed
zeros while checking callback-sensitive replacement results.

=== Complete ladder measurements

The notebook now uses `replace_tensor` for the typed outside-in route. Three
counterbalanced fresh processes per route and order compare the saved `b2ee`
build with the new `60e11d75` build, including both replacement methods in the
new build. Medians of the complete eight-step reduction are in seconds:

#table(
  columns: 3,
  table.header([Route], [Original order], [Early rungs]),
  [Saved typed generic replacement], [2.632], [1.138],
  [New build, generic replacement], [2.303], [0.830],
  [New build, tensor-safe replacement], [2.157], [0.686],
  [Supplied scalar Symbolica recipe], [0.962], [0.344],
  [FORM complete process], [0.791], [0.186],
)

Within the same build, the ratio of medians improves by *6.4% / 17.4%*. Every
paired round improves in both orders. The saved-baseline samples vary widely:
2.416–3.097 s in the original order and 0.814–1.325 s with early rung closure.
The table retains those medians; it does not establish a universal percentage
gain. Tensor-safe early-rung samples are 0.669, 0.686 and 0.718 s.

The typed route remains *2.73× / 3.68× slower than FORM*. Its clocks include
replacement, expansion and contraction, excluding setup, counting and exact
checks; FORM includes process startup. All 30 typed reduction loops retain the
same 9652-term scalar, exact FORM equality and unchanged reruns. Intermediate
term counts are unchanged: original-order and early-rung peaks remain 25662
and 11836, respectively.

Separate 199 Hz profiles attribute roughly 33–36% of sampled cycles to the
Schoonschip kernel, 26–27% to replacement traversal, matching and reconstruction
(including its intrinsic callback guard), and 17–22% to other Symbolica
normalization. Local replacement-signature observation accounts for 1.5–2.8%,
with another 4.6–5.2% in the source-interface proof. These do not include all
validation: the intrinsic guard contributes 2.2–2.8% within the replacement
category. These are disjoint sampled CPU categories from separate diagnostic
loops, not wall-time fractions. The profiles contain 371 / 136 operation
samples, so small categories have limited precision.

The controls show no general trace speedup. Free length-10 tracen increases
from 1.590 to 1.668 ms, with all three paired processes slower. Free length-12
tracen changes from 20.052 to 21.022 ms with mixed paired signs. Axial
length-12 trace4 is essentially unchanged at 1.549 to 1.552 ms. The reverse
accumulator median also increases by 3.5%. These observations remain in the
record; their causes are not established. The retained improvement targets
the complete ladder through the new operation.

Validation passes all 910 Rust tests, Clippy and formatting, the unchanged
207 API records and 117 HEP component rows, 32 new public replacement controls
and 26 existing controls. Actual notebook cells pass in both orders, and the
complete factored eight-vertex network agrees at two exact component
assignments. The public tests compare logical ports independently of inferred
tensor names and also require descriptor parity with general replacement.

Sources, commands, raw timings, profiles and checks are preserved in the local
measurement bundle at `/tmp/idenso-tensor-safe-replace-host`. The
`retained_tensor_safe_replacement` checkpoint is prepared for the
#source-link("examples/notebooks/tensor_contraction_parity.json", label: "contraction record");
appending it awaits explicit approval. The existing archive is unchanged.

=== Symbol classification and remaining contraction walks

The tensor classifier now borrows Spenso's initialized tag bundle once per
classification. It checks a symbol's symmetry attributes before comparing it
with a projector symbol. Exact symbol equality still distinguishes projectors
from ordinary tensors carrying the same attributes. Registered and imported
projectors retain their nested ports; unrelated attributed tensor leaves retain
their declared ports and scalar metadata.

Outside initializer reentry, public lazy-bundle access synchronizes with
Symbolica's registered initializers through a discarded `State::is_builtin`
lookup. In the saved
`60e11d75` profiles, those lookups account for about 5.9% / 6.6% of inclusive
sampled cycles. An isolated warm primitive measures 10.61 ns for a bundle
access and 0.657 ns through a held reference. These primitive clocks do not
predict a complete-ladder speedup. The change retains the initialization
barrier, reentry handling and license checks; no global readiness cache was
introduced. A concurrent reproducer confirms that merely reading an existing
symbol can return before registered initializers finish.

A fresh 24-process matrix compares the saved `60e11d75` build with the retained
`9c8b475a` classifier cleanup. Each row has three fresh-process samples per
order; the typed variants run in counterbalanced rounds. Complete reduction
medians are seconds:

#table(
  columns: 3,
  table.header([Route], [Original order], [Early rungs]),
  [Saved build, tensor-safe replacement], [1.764], [0.676],
  [Classifier cleanup, tensor-safe replacement], [1.765], [0.663],
  [Saved build, generic replacement], [2.309], [0.822],
  [Classifier cleanup, generic replacement], [2.344], [0.823],
  [Supplied scalar Symbolica recipe], [1.003], [0.349],
  [FORM complete process], [0.827], [0.198],
)

The matched tensor-safe medians change by *+0.08% / −1.98%*. Paired samples
have mixed signs, so this does not establish a speedup in both orders. The
cleanup is retained as a small removal of redundant initialization checks.
The generic controls change by +1.52% / +0.07%. Cross-session differences are
larger than this change: the lower original-order time relative to the earlier
table also occurs in the unchanged saved build and must not be attributed to
the classifier. All samples, including the slower first baseline samples,
remain in the record.

The retained tensor-safe route remains *2.14× / 3.35× slower than FORM* in this
session. Python clocks cover the eight reduction steps; FORM clocks include
startup. All 24 reductions preserve the exact scalar result and intermediate
counts. The unchanged notebook cells pass in both orders against fresh FORM;
912 Rust tests, Clippy, formatting, 207 API records, 117 HEP rows, 58 public
controls and the complete network component checks pass. The complete
measurement bundle is `/tmp/idenso-initialization-guard-host`.

Trace controls remain mixed. First-transformation medians below are milliseconds;
the Python method clock and FORM's batched body CPU clock have different
boundaries, unlike the complete-process FORM ladder clock above.

#table(
  columns: 4,
  table.header([Case], [Saved build], [Classifier cleanup], [FORM body CPU]),
  [Free 2, tracen], [0.02150], [0.02094], [0.00070],
  [Free 4, tracen], [0.07532], [0.07432], [0.00145],
  [Free 6, tracen], [0.12253], [0.12138], [0.00560],
  [Free 10, tracen], [1.65475], [1.63159], [0.35667],
  [Free 12, tracen], [20.72884], [20.95333], [4.43333],
  [Free 8, trace4], [0.23939], [0.24758], [0.04933],
  [Repeated 8, tracen], [0.16067], [0.17939], [0.01440],
  [Repeated 8, trace4], [0.09369], [0.09209], [0.00200],
  [Axial 12, trace4], [1.86560], [1.65081], [0.64800],
)

The repeated length-8 tracen median increases by 11.7%, with mixed paired
signs. Free length-12 tracen's rerun control increases from 9.369 to 10.567 ms,
also with mixed paired signs. Axial trace4's first median decreases by 11.5%,
but its first process pair is essentially flat. The free length-8 trace4 rerun
is slower in all three pairs (1.3–7.6%). These observations do not establish a
general trace improvement; all exact output checks still pass.

Separate untimed counters explain the contractor's remaining discovery work.
They use the actual inputs to stages 1–7 of both tensor-safe ladder orders and
match all saved stage outputs exactly:

#table(
  columns: 3,
  table.header([Counter], [Original order], [Early rungs]),
  [Actual substitutions], [82,042], [42,444],
  [Discarded first-hit substitutions], [7], [7],
  [Discovery product searches], [263,406], [79,313],
  [Discovery visitor nodes], [3,007,991], [1,000,922],
)

Each active stage performs one modifying map followed by one unsuccessful
discovery pass. Reusing the seven discarded first results would save little;
the final no-work traversal is the stronger candidate. Callback purity alone
does not justify skipping it: an unchanged factor may contain a metric product
inside scalar metadata. A locality certificate would also need to exclude
that work and account for normalization that exposes products. The production
contractor retains its fixed-point walk. The scalar eighth-stage control is
excluded from the table because the public pipeline skips contraction there.

The initialization source, concurrency checks and primitive samples are in
`/tmp/idenso-initialization-probe`; the contraction counters and exact stage
fixtures are in `/tmp/idenso-certified-contraction-audit/diagnostic`.

=== Completion-scan experiment (not retained)

An isolated candidate (`72384254`) extended the existing index observation with
the location of each explicit metric/vector source. When all sources belonged
to outer products, normalization was intrinsic, and no replacement exposed
new arithmetic, the contractor could omit its final unsuccessful traversal.
Hidden sources, callbacks, component emission and chain rewriting retained the
existing fixed-point search. The modifying pass and arithmetic order stayed
unchanged.

The initial eight-process screen improved both order medians by about 4%, but
the separate 24-process matrix did not establish a consistent gain. Complete
reduction medians are seconds; each typed row has three fresh-process samples:

#table(
  columns: 3,
  table.header([Route], [Original order], [Early rungs]),
  [Retained build, tensor-safe replacement], [1.994], [0.689],
  [Experimental completion shortcut], [1.756], [0.891],
  [Retained build, generic replacement], [2.596], [0.848],
  [Experimental build, generic replacement], [2.708], [0.813],
  [Supplied scalar Symbolica recipe], [0.959], [0.354],
  [FORM complete process], [0.792], [0.176],
)

Two early-rung candidate samples were 27% and 39% slower; the third was 8%
faster. A prespecified follow-up of four counterbalanced pairs added process
and thread CPU clocks around the unchanged reduction. All four pairs improved,
with wall medians 0.727 → 0.670 s and process CPU medians 0.718 → 0.662 s.
It also contained a slow 1.013 s baseline sample. Scheduler waiting was small
in that follow-up; these observations cannot retrospectively explain the
original slow candidate samples. The screen, matrix and diagnostic remain
separate, with no discarded samples.

Trace controls showed no general improvement: free length-12 tracen changed
from 19.170 to 20.255 ms and axial length-12 trace4 from 1.564 to 1.615 ms,
both with mixed paired signs. Free length-10 tracen increased by 1.25%, with
all three pairs slower.

Correctness passed 917 Rust tests, Clippy, formatting, 207 API records,
117 HEP rows, 58 public controls, 24 adversarial cases, all 16 saved stages
and fixed points, actual notebook routes against fresh FORM, and the complete
factored ladder component checks. Intermediate counts were unchanged.
The seven-file experimental delta was nevertheless restored to the saved
`9c8b475a` sources: the added complexity did not demonstrate a consistent
end-to-end benefit. The installed build and tensor-safe replacement are
unchanged. In this session the retained build remains 2.52× / 3.93× slower
than FORM. Sources, raw samples, checks and the decision are preserved under
`/tmp/idenso-contraction-completion-host`.

=== Reusing Symbolica's tensor-factor matcher

Tensor-factor replacement now keeps one lazily constructed
`AtomMatchIterator` and its `WrappedMatchStack` for the operation. It clears
bindings before each new factor and uses `next()` to require fully satisfied
conditions. This replaces the former per-factor tree iterator, whose configured
search visited only the root. A fixed function-name pattern also skips factors
with another name before tensor classification. Wildcard heads, alternatives
and optional arithmetic patterns keep Symbolica's matching semantics.

The shared Idenso implementation still observes the actual matched factor,
checks each normalized RHS against its interface and explicit-index counts,
and validates cached results against each actual target. Conditions run before
cache lookup; RHS callbacks retain their existing cache behavior. Lazy matcher
construction preserves the existing behavior when no eligible leaf is reached.
No interface proof, callback guard or Python conversion behavior was removed.

The retained `ac275b56` build is compared with saved `9c8b475a` using three
counterbalanced fresh processes per route and order. Complete eight-step
reduction medians are seconds; setup and exact checks are outside the typed
clock, while FORM includes process startup:

#table(
  columns: 3,
  table.header([Route], [Original order], [Early rungs]),
  [Saved build, tensor-safe replacement], [1.723], [0.655],
  [Reused matcher, tensor-safe replacement], [1.637], [0.604],
  [Saved build, generic replacement], [2.376], [0.816],
  [New build, generic replacement], [2.271], [0.825],
  [Supplied scalar Symbolica recipe], [0.955], [0.360],
  [FORM complete process], [0.811], [0.180],
)

Tensor-safe medians decrease by *4.97% / 7.76%*. All six matched process pairs
improve: 1.5–6.5% in the original order and 4.6–11.0% with early rungs; process
CPU clocks give similar changes. Within the new build, tensor-safe replacement
is *27.9% / 26.8% faster* than generic replacement. The generic original-order
control also improves by 4.45%, so the entire original-order build difference
cannot be attributed solely to matcher reuse. These small samples establish
this measured result, not a universal speedup.

The separate initial eight-process screen improved three pairs, with the fourth
4.76% slower. Its order medians changed by −4.56% / −1.03%; it is not pooled
with the full comparison. An earlier fixed-name-only candidate retired
4.1–5.1% fewer user instructions in four hardware-counter pairs, but wall and CPU
times remained mixed. In its slow early-rung pair, fewer instructions coincided
with more cycles per instruction and cache misses. That diagnostic does not
identify the stall source or explain earlier unprofiled slow samples. The
intermediate candidate was not installed separately.

Final output remains the same 9652-term scalar, exactly equal to FORM.
Intermediate counts are unchanged: original-order and early-rung peaks remain
25662 and 11836. The retained typed route is still *2.02× / 3.35× slower than
FORM*. This change reduces replacement overhead; it does not change contraction
order or the number of generated terms.

Trace controls show no general improvement. First-transformation medians below
are milliseconds. The Python method clock and FORM's batched body CPU clock
have different boundaries:

#table(
  columns: 4,
  table.header([Case], [Saved build], [New build], [FORM body CPU]),
  [Free 2, tracen], [0.02079], [0.02012], [0.00066],
  [Free 4, tracen], [0.07261], [0.07050], [0.00140],
  [Free 6, tracen], [0.12126], [0.11961], [0.00540],
  [Free 10, tracen], [1.57321], [1.62166], [0.36000],
  [Free 12, tracen], [20.05999], [20.39048], [4.20000],
  [Free 8, trace4], [0.24481], [0.23886], [0.04867],
  [Repeated 8, tracen], [0.15591], [0.15489], [0.01540],
  [Repeated 8, trace4], [0.09393], [0.09460], [0.00200],
  [Axial 12, trace4], [1.60134], [1.61769], [0.65200],
)

Axial trace4's first transformation is slower in all three pairs. Rerun medians
increase by 0.09–3.92%; free length-4 tracen and repeated length-8 trace4 reruns
are slower in all three pairs. These observations are retained without a causal
claim about code paths unchanged by this patch.

Validation passes 918 Rust tests, Clippy and formatting, 207 API records,
117 HEP rows, 84 public controls, the actual notebook's three routes in both
orders, and the complete factored ladder at two independent exact component
assignments. Regressions compare condition and RHS callback transcripts across
cache settings, including failed matches, backtracking, changing arities,
symmetric metric arguments, and lazy construction. The two initially invalid
test fixtures were corrected to use admitted tensors; production checks were
unchanged.

The previous profile still identifies contraction and normalization as major
remaining costs; no new profile is claimed for this build. A further concrete
allocation candidate is the RHS cache key: lookup currently clones bindings
even on hits and when the cache is full. A borrowed lookup could defer that
copy until insertion; it is not implemented in this measurement.

Sources, raw samples, exact checks and the final summary are preserved under
`/tmp/idenso-reused-matcher-host`; the fixed-name-only diagnostic is under
`/tmp/idenso-fixed-head-replacement-host`. The existing 85-entry benchmark
archive is unchanged.

=== Borrowing RHS cache keys

Tensor-factor replacement now looks up cached right-hand sides using the
borrowed match bindings. It copies the bindings only when inserting a new cache
entry. Hits, disabled caches and full caches avoid the former temporary vector
and copies of sequence-wildcard bindings. Matching, conditions, normalized RHS
checks and the check against each actual target are unchanged.

This remains a tensor-safe operation: whole factors are replaced underneath
sums and products, and local interface checks justify retaining the surrounding
logical ports. It avoids full-result structure inference, but still scans the
source to establish that local checks suffice. It certifies structure, not the
mathematical validity of a user-supplied identity. Remaining normalization hooks
that can change rank require the general checked route; an RHS callback may
return an admitted normalized expression, which is checked before reuse.

The saved `ac275b56` and retained `d9846eac` builds were compared in three
counterbalanced fresh processes per route and order. All complete ladder
samples below are seconds, in paired run order:

#table(
  columns: 3,
  table.header([Order], [Saved tensor-safe build], [Borrowed cache lookup]),
  [Original], [1.682, 1.782, 1.796], [1.572, 1.555, 1.758],
  [Early rungs], [0.705, 0.711, 6.698], [0.607, 1.127, 0.587],
)

All three original-order pairs improve by 2.1–12.7%; the medians change from
1.782 to 1.572 seconds. The early-rung samples do not establish a reliable
speedup. The 6.698-second saved-build sample consumed 5.287 seconds of process
CPU time, while the 1.127-second candidate sample consumed 1.116 seconds.
Scheduler delays alone therefore cannot explain the variation; no cause was
established. All samples are retained. The separate eight-process screen is
not pooled with this comparison; three pairs improved and one was 2.3% slower.

Fresh FORM process medians are 0.809 / 0.196 seconds and the supplied scalar
Symbolica recipe takes 0.950 / 0.350 seconds. The candidate's original-order
median remains 1.94 times FORM. Generic-replacement controls also vary, with
medians changing from 2.400 / 0.986 to 2.585 / 0.847 seconds. These results do
not establish a general build-wide gain. Final results remain exactly equal
to FORM's 9652-term scalar; intermediate counts and contraction order are
unchanged.

Nine trace4/tracen cases and their reruns remain exact. These controls do not
use the changed RHS cache, so their timing changes are not attributed to this
patch. Free length-12 tracen changes from 19.087 to 19.588 milliseconds, and
axial length-12 trace4 from 1.573 to 1.619 milliseconds. FORM's corresponding
batched body CPU times are 4.133 and 0.604 milliseconds; the Python method and
FORM body clocks have different boundaries.

Validation passes 918 Rust tests, Clippy, formatting, 207 API records,
117 HEP rows, 84 public controls, all three notebook routes in both orders,
and the complete factored ladder at two exact component assignments. The
cache regression now covers disabled, saturated and reusable caches with
sequence wildcards, including condition and RHS callback counts. Source and
measurement evidence is preserved under `/tmp/idenso-borrowed-rhs-cache-host`;
the existing 85-entry benchmark archive is unchanged.

=== Planning contractions into opaque tensors

The existing metric-component collector now admits plain tensor factors with
explicit ports and already bound compact vectors. Each term supplies an index
graph: metrics connect slots, vectors terminate paths, and tensor terminals
identify argument positions. The collector resolves those positions and interns
the resulting tensor argument lists across sum terms before constructing their
Symbolica functions. Repeated source factors reuse their admitted port plans
and metadata checks. Resolved explicit endpoints are cached across terms, while
incidence counts remain local to each term; tensor scratch storage retains its
capacity. The shared Idenso tensor inference still owns the leaf interface
proof; Python dispatch is unchanged.

For the 32-term `g(a,b_i)*T(b_i)` control, this reuses 33 resolved slots rather
than making 288 endpoint-parser calls inside the collector. The existing
interface proof still runs; these operation counts are not a timing result.

Admission excludes normalization hooks, unresolved ports, incompatible spaces,
unsupported powers and ambiguous index multiplicities. Scalar metadata is
retained only when the existing pending-work observer certifies that it needs
no further contraction or normalization. Failed admission constructs no
speculative tensor results. This preserves the callback-sensitive case where
relabeling a tensor makes its normalizer return a scalar.

Two tensor terminals can retain a dummy index. The implementation preserves the
ordered contractor's surviving label, including its whole-product threshold
for the three-metric shortcut. A regression compares 912 small path arrangements
with the unchanged ordered contractor, including two ports on one tensor and
disconnected metrics. Other regressions cover bound vectors, opaque metadata,
local power pairing, callback behavior, logical port order and typed zeros.

All fourteen saved inputs to ladder stages one through seven are admitted and
match the existing intermediate expressions exactly. The largest original-order
case plans 339 variables, including 313 distinct tensor results; its dense
exponents require 24.60 MB (23.46 MiB), within the existing 64 MiB bound. Closing
rungs early needs at most 86 variables and 3.54 MB (3.37 MiB) on these inputs.
The untimed source-included diagnostic uses the public full-subtree intrinsic
check in place of the private head predicate; both admit these fixtures. It
establishes workload coverage, not an end-to-end speedup.

The native metric benchmark now calls the actual private contractor from the
existing `reference-cases` library module. Its example is a thin launcher and
requires `--features reference-cases`. Inputs, options, output snapshots and
timing loops are preserved; the source-included copy of the contractor is
removed.

The retained `63ba01ba` build is compared with the saved `d9846eac` build, which
already has tensor-safe replacement and borrowed RHS cache lookup. Both use
the same pinned Symbolica revision. In the initial full matrix, three
counterbalanced process pairs per route and order give these median seconds:

#table(
  columns: 5,
  table.header([Order], [Saved tensor-safe], [New collector], [FORM process], [Scalar recipe]),
  [Original], [1.916], [1.812], [0.960], [1.032],
  [Early rungs], [0.695], [0.552], [0.211], [0.370],
)

All three original-order pairs improve. Two early-rung pairs improve, while
the third changes from 0.633 to 0.851 seconds. That slow sample uses 0.837
seconds of thread CPU, so scheduler waiting alone does not explain it. Its
busy SMT sibling is a possible confounder, not a demonstrated cause. Generic
replacement controls have medians 2.593 to 2.464 seconds and 0.902 to 0.798
seconds. FORM measures complete process wall time, whereas Python measures
the complete eight-stage reduction after imports and fixture setup.

To assess this variation, a further fixed block runs eight fresh process pairs
per order and eight metric pairs. It retains every sample and is not pooled
with either the initial screen or the three-pair matrix. All sixteen ladder
pairs improve in this larger block:

#table(
  columns: 4,
  table.header([Order], [Saved median, s], [New median, s], [Candidate/saved paired range]),
  [Original], [1.830], [1.627], [0.740–0.995],
  [Early rungs], [0.624], [0.508], [0.758–0.905],
)

The median reductions are 11.1% and 18.5%; paired CPU clocks show the same
direction. These are shared-host observations, with the earlier contrary
sample preserved. They do not establish general FORM parity. The final scalar
has the same 9652 terms, and the notebook's stage counts are unchanged.
This improvement reduces work per stage; it does not remove
the existing difference in intermediate counts relative to FORM.

Small metric workloads do not improve uniformly. The eight-pair block gives
the following microseconds. The first two columns time the Python contraction
method, including dispatch and result wrapping. The next column also parses
canonical source text and constructs the typed input. FORM's separate
three-process batches include input declarations, contraction and the initial
sort, since FORM can contract metrics while reading the input. These are
different clock boundaries, not a pure-kernel comparison.

#table(
  columns: 5,
  table.header([Case], [Saved method], [New method], [New parse + construct + contract], [FORM input + contract CPU]),
  [Chain 8], [19.604], [19.563], [128.447], [5.000],
  [Chain 64], [109.187], [110.385], [884.715], [44.000],
  [Loop 64], [104.390], [107.137], [864.242], [42.000],
  [Opaque sum 32], [123.202], [168.078], [809.917], [46.000],
)

The opaque sum is slower in all eight pairs, by 14.9–40.2%; its median grows
by 36.4%, or 44.9 microseconds. The retained implementation prioritizes the
complete ladder gain and leaves this small-sum planning overhead unresolved.
The other metric cases have mixed paired results. No-op reruns remain exact.

The full matrix also retains nine trace controls and their reruns. Below are
first-call method wall-time medians in milliseconds; FORM gives batched trace-body CPU
time. These clocks have different boundaries. Several small trace controls
regress, and the variation does not support a general trace speedup:

#table(
  columns: 4,
  table.header([Case], [Saved], [New], [FORM body CPU]),
  [Free 2, tracen], [0.02092], [0.02056], [0.00072],
  [Free 4, tracen], [0.07189], [0.07371], [0.00155],
  [Free 6, tracen], [0.12156], [0.11997], [0.00640],
  [Free 8, trace4], [0.27796], [0.38502], [0.05000],
  [Repeated 8, trace4], [0.09411], [0.15894], [0.00220],
  [Repeated 8, tracen], [0.15694], [0.26465], [0.01480],
  [Free 10, tracen], [1.72363], [1.60299], [0.35000],
  [Free 12, tracen], [23.99286], [21.36665], [4.16667],
  [Axial 12, trace4], [3.08722], [2.87527], [0.64000],
)

Reruns do not improve uniformly either. Free length-12 tracen reruns have
medians 10.234 to 16.594 milliseconds, a 62.1% increase; all three paired
ratios exceed one (1.115, 1.004 and 1.621). The lower first-call median for
that case must not be read as a general improvement to its lifecycle. These
observations retain the full matrix's host-variation qualification.

A separate 199 Hz profile of the retained build records 389 original-order
and 101 early-rung samples inside the timed phases. Its disjoint attribution
assigns about 31–34% to Schoonschip, 22–27% to other Symbolica normalization,
19–20% to replacement traversal, matching and reconstruction, and 6–8% to
the source-interface proof. Coefficient-list collection and conversion back
to Symbolica expressions remain visible contraction costs. These small
sample counts guide further investigation; profile clocks are not benchmark
measurements or an exact wall-time decomposition.

Validation passes 922 Rust tests with 22 existing skips, Clippy, formatting,
207 API records, 117 HEP rows and 149 targeted controls. The complete factored
ladder passes four component evaluations across two assignments. All three
notebook routes agree exactly with fresh FORM results in both orders, including
the inactive-cell behavior. Callback rank loss, typed zeros, missing branch
indices, unresolved ports and logical order retain their checked behavior.
Sources, every timing sample and validation evidence are saved under
`/tmp/idenso-opaque-components-cache-host`, with independent controls under
`/tmp/idenso-opaque-components-cache-review`. The earlier `b14f5d66` experiment
is preserved separately; the existing 85-entry benchmark archive is unchanged.

=== Sparse collection without the polynomial boundary

An isolated experiment replaces the final dense coefficient-list polynomial
with sparse monomial keys and exact rational coefficient accumulation. It
retains the existing admission checks and canonical variable aliasing, then
uses bulk Symbolica sum/product construction only for nonzero collected
monomials. Unlike the earlier direct trace emitter, it combines duplicates
before constructing product Atoms.

The experiment passes 12 focused tests, including the existing 912 metric-path
comparisons and new alias, cancellation, large-rational and exponent-limit
checks. A matched native screen covers all 16 saved ladder stage inputs and
eight metric/power controls, checking exact outputs and fixed points. Three
counterbalanced process pairs retain all 864 timing batches. The following
selected input-contraction medians are milliseconds; ranges include all three paired
candidate/baseline wall-time ratios:

#table(
  columns: 4,
  table.header([Contraction input], [Polynomial], [Sparse], [Paired ratio range]),
  [Original stage 4], [9.405], [11.353], [0.993–1.494],
  [Original stage 6], [121.655], [123.299], [0.912–1.103],
  [Original stage 7], [221.667], [192.401], [0.825–0.942],
  [Early-rung stage 4], [3.581], [4.176], [1.037–1.701],
  [Early-rung stage 5], [10.557], [12.823], [1.012–1.602],
  [Early-rung stage 6], [32.851], [35.739], [1.011–1.145],
  [Early-rung stage 7], [89.594], [99.283], [0.973–1.141],
)

The largest original-order input improves in every pair, but several
early-rung inputs regress in every pair. The opaque 32-term control is mixed,
with medians 76.1 and 77.0 microseconds. These clocks measure the native
contraction method on preconstructed inputs; parsing, exact-result checks and output
destruction are outside the measured batches. They are not complete-ladder
or Python method timings. The already scalar eighth-stage inputs do not
exercise coefficient collection.

This variant is not retained or installed, and no complete-host gain is
claimed. All observations, compiler identities, source patches and checks are
saved under `/tmp/idenso-component-sum-emission`. The installed implementation
and its FORM/trace comparisons above remain unchanged.

An independent shared-interface experiment reuses already validated explicit
self-dual ports for ordinary tensor leaves, avoiding construction of a
temporary interface merely to check the same ports again. Builtin validation,
unresolved ports, metadata and callback-sensitive paths retain their existing
behavior. All 18 focused tests pass, including 123 comparisons with the old
proof and its callback schedule. This experiment does not include the sparse
collector change.

Its separate three-pair native screen improves the opaque 32-term case from
161.8 to 138.2 microseconds for contraction plus the shared typed result
operation, with every pair faster. The proof alone improves from 76.2 to 60.8
microseconds. Larger typed inputs remain mixed: original stage 6 has medians
145.0 to 157.1 milliseconds, with two of three pairs slower; original stage 7
improves 6.3%, early-rung stage 6 improves 9.9%, and early-rung stage 7 improves
14.3% by ratio of medians. The original stage-6 proof itself is slower in all
three pairs. Large baseline outliers remain in the saved observations, and
these measurements do not establish an end-to-end ladder improvement. The
variant remains isolated and unselected while work prioritizes the complete
ladder. Sources, old/new proof comparisons and every sample are saved under
`/tmp/idenso-explicit-leaf-proof`.

=== Locating the first intermediate-term difference

An untimed comparison of the retained build isolates the first difference in
the rung-closing order, after vertices 5, 4 and 6. The complete expressions
agree exactly after vertices 5 and 4, with 6 and 34 terms respectively. The
first syntax difference is created at vertex 6. Both the tensor-safe route
and the supplied scalar outside-in recipe produce 74 terms; FORM produces 63.
Converting the typed result to the scalar recipe's notation gives exactly the
same 74-term expression.

Subtracting the FORM expression leaves 14 pairs of terms. Each pair has
equal-magnitude, opposite-sign coefficients and identical momenta and
remaining vertices, but uses a different dummy label for the connection between
vertices 3 and 7: `mu11` versus `mu6`.
Renaming that label in this particular expression gives exactly FORM's 63
terms. Independently, Symbolica's existing tensor canonicalization gives equal
63-term expressions for both inputs. This is an exact algebraic diagnostic,
not a measured speedup or a general prescription to replace index names.

Feeding the saved, fully expanded 143-term input to FORM produces exactly
Symbolica's 74-term expression. The notation conversion is checked by an exact
round trip before execution. Thus the difference precedes contraction of the
collected input; it arises while the vertex replacement is constructed and
reduced.

A minimal FORM probe exposes the relevant ordering: `d_(a,b)*t(a)*u(b)` becomes
`t(b)*u(b)`, whereas writing `d_(b,a)` instead leaves `t(a)*u(a)`. FORM can
contract using that written argument orientation. Both the scalar recipe's
symmetric `d` and Spenso's symmetric `g` erase the orientation on construction:
their two argument orders already have identical Symbolica bytes before any
contraction.

One controlled change ties this mechanism to the complete third stage.
Reversing only the two `d_(mu2,mu3)` factors in FORM's vertex-6 replacement
changes its 63-term result to 74 terms, exactly equal to Symbolica's result.
Every other program byte is preserved, and reversing the patch recovers the
original program. All results still canonicalize to the same 63-term tensor.
These are untimed correctness comparisons; they quantify no runtime gain.

The probes also confirm that FORM contracts into ordinary `CFunction`
arguments. Index declaration order does not change the ladder counts, and
moving the scalar recipe's contraction from the held RHS into the ambient
product does not recover the orientation already erased by construction.
This behavior does not imply general dummy canonicalization in FORM:
`t(a)*u(a)-t(b)*u(b)` remains a two-term expression in a separate control.

The distinction matters for coefficient collection: products with different
dummy spellings are different Symbolica Atoms, so collecting equal products
alone cannot merge them. Symmetric metric construction remains algebraically
correct. A contraction strategy that combines more terms must choose dummy
representatives consistently while preserving free indices, avoiding dummy
collisions and respecting opaque metadata. The current contraction policy
remains unchanged. The scripts, exported expressions and exact checks are
saved under `/tmp/idenso-form-cfunction-boundary`, including the complete
orientation A/B under `orientation-stage`.

The existing scalar `canonize_tensors` operation was also checked after every
vertex of both complete ladders. It recovers FORM's retained counts at every
stage, with exact agreement on the final 9,652-term polynomial:

#table(
  columns: 3,
  table.header([Order], [Current typed contractions], [Scalar canonicalization / FORM]),
  [Original], [6, 34, 192, 1084, 6069, 11959, 25662, 9652],
  [6, 34, 192, 1084, 6069, 10947, 23937, 9652],
  [Early-rung], [6, 34, 74, 396, 944, 5046, 11836, 9652],
  [6, 34, 63, 352, 829, 4586, 10516, 9652],
)

This test uses the supplied scalar recipe and its eleven explicit index
symbols. Final contraction and canonicalization fixed points pass in both
orders. It does not establish typed metadata handling or a runtime gain.
The complete checks and the first-stage typed canonicalization failure are
saved under `/tmp/idenso-existing-canonization`.

Three counterbalanced process pairs per order measure the complete scalar
recipe with and without this extra pass. All twelve runs preserve the exact
final result and fixed points. Median whole-loop wall times, in seconds, are:

#table(
  columns: 4,
  table.header([Order], [Scalar recipe], [With canonicalization], [Fresh FORM process]),
  [Original], [1.110805], [1.989898], [0.832198],
  [Early-rung], [0.346671], [0.681879], [0.186588],
)

Every pair is slower with canonicalization: wall-time ratios range from
1.394 to 2.126 in the original order and 1.271 to 2.273 in the early-rung
order. Process and thread clocks confirm the direction. These measurements
exclude scalar setup, final projection and exact checks; FORM's separate
three-run medians include process startup and parsing. The unchanged
`63ba01ba` host, frozen notebook recipe, all observations and CPU/SMT-sibling
snapshots are recorded under `benchmarks` in that experiment directory.
No samples are discarded. General canonicalization after every step is not
adopted: its cost exceeds the savings from fewer intermediate terms.

Two separate instrumented runs attribute 0.701 seconds of a 1.561-second
original-order loop and 0.270 seconds of a 0.573-second early-rung loop to
canonicalization itself, about 45% and 47%. These single diagnostic runs add
per-stage clocks and counts; they are not pooled with the primary timings
or used to estimate a paired speedup.

Two isolated native experiments also test borrowing unchanged tensor functions
and allocating argument arrays only when a port changes. The second variant
uses the existing component graph to skip identity writes, including pairs
without metric edges, and avoids reading an old slot when the replacement
is necessarily a vector. Both variants pass 12 focused tests, including the
912 existing metric-path comparisons, and exact output/fixed-point checks on
all 16 saved ladder stages plus eight controls.

Each variant has its own three-pair screen, retaining all 864 batch samples.
Neither demonstrates a broad gain on the large stages. For the refined variant,
original stage 7 has contraction-method medians 212.747 to 254.424 milliseconds,
with paired ratios 1.196, 1.240 and 1.043; early-rung stage 7 has medians 89.709
to 89.785 milliseconds. The already scalar early-rung stage-8 control also
varies from 13.959 to 21.398 milliseconds, so these observations do not establish
that the patch causes every slowdown. Neither variant is retained or used to
claim a complete-ladder gain. All samples, rerun controls and decisions remain
under `/tmp/idenso-lazy-tensor-arguments` and
`/tmp/idenso-lazy-tensor-arguments-refined`.

=== Keeping scalar metadata opaque during materialization

The typed canonicalization probe exposed a separate parser defect. Compact
vectors inside a `Scalar` function were lowered into fresh tensor factors.
For example, `T(F(g(p(rep),q(rep))),a)*U(a)` could become
`p(d)*q(d)*T(F(1),a)*U(a)`, although `F` is scalar metadata and may be
nonlinear. Different fresh dangling indices in sum branches also caused the
network-addition panic seen on the ladder's first six-term expression.

The shared Spenso shorthand materializer now stops at Scalar functions in
both detection and lowering. Scalar attributes take precedence over tensor,
rank-one and broadcast tags. Actual compact vectors in tensor argument
positions still lower normally, leaving adjacent metadata unchanged. The
fix adds four early guards in the existing owner, with regressions for
nonlinear metadata, nested wrappers, sums, powers, conflicting tags and the
allocation of exactly one dummy for a real bound vector. No tensor type or
Python implementation is added.

This source correction is separate from the installed `63ba01ba` host used
for the timings above. That host already provides `replace_tensor` and runs
the notebook's tensor-safe route; these measurements do not include the new
materializer guards.

Nine materializer tests pass. Five native pipeline checks also pass, including
the two actual six-term ladder inputs: all 21 routing metadata fields survive
literally in each input, the interface remains scalar, canonicalization reaches
a fixed point, and an independent dummy-label comparison confirms equality
after Schoonschip contraction. These are correctness checks, not timings.
That worktree checkpoint also passes all 925 library tests across Spenso, Idenso
and Spynso, with 22 existing skips, plus Clippy for all targets and formatting.
Source snapshots remain unchanged throughout those checks.

There is an independent remaining limitation in the general canonicalizer:
`T(F(mink(4,a)),mink(4,a))*U(mink(4,a))`, with `F` Scalar, constructs as a
scalar tensor but canonicalization reports that the index occurs too often.
Symbolica's recursive labeling sees the slot hidden in metadata as well as
the two real ports. The materializer correction does not claim to fix this
separate labeling boundary. Its reproducer is saved as
`checks/hidden-fullslot-summary.json` in the canonicalization experiment;
the materializer patch and regressions are under
`/tmp/idenso-materializer-scalar-opacity`.

=== Collecting equivalent internal connections

The shared Idenso component contractor now recognizes products that differ
only by their internal index names. Its key records the tensor heads, exact
noninternal arguments, ordered port positions, representation spaces,
connections, and exact non-tensor factors and powers. Free labels and scalar
metadata remain literal. Equal keys reuse one complete contracted product as
the coefficient-collection representative;
the operation does not generate or substitute fresh dummy names.

This deliberately covers uniquely distinguishable tensor factors. Repeated
indistinguishable factors retain their literal result, since identifying those
graphs requires more general labeling. Estimated stored-key and transient work
budgets are bounded separately. Exceeding either budget declines the entire
collection before materialization. If ordinary contraction or dot cleanup
makes a previously declined sum eligible, the public simplifier finishes that
collection before returning, preserving its fixed point.

The same collector caches successful compact-representation recognition across
terms. The complete representation Atom is the key; dimension, variance and
syntax validation still run on misses. This reuses the existing slot matcher
and component plan. Neither change introduces another tensor abstraction or
Python algebra implementation.

The two complete ladder orders now match FORM's retained term counts at every
stage:

#table(
  columns: 3,
  table.header([Order], [Saved collector], [Connection collection / FORM]),
  [Original], [6, 34, 192, 1084, 6069, 11959, 25662, 9652],
  [6, 34, 192, 1084, 6069, 10947, 23937, 9652],
  [Early-rung], [6, 34, 74, 396, 944, 5046, 11836, 9652],
  [6, 34, 63, 352, 829, 4586, 10516, 9652],
)

Regressions cover self-connections, multiple connections, distinct spaces,
free ports, opaque metadata, rational coefficients, local powers, compact
vector aliases, typed zero with reordered logical ports, budget exhaustion and
cleanup that changes eligibility. Existing callback-sensitive paths retain
their checks. The integrated source passes all 932 library tests across Spenso,
Idenso and Spynso, with 22 existing skips, all-target Clippy and scoped
formatting checks. The test helper also now qualifies complete function heads,
so registering Scalar `routing` tests that actual head rather than accidentally
rewriting the end of its name.

The isolated native comparisons explain the selection without establishing a
Python speedup. Compact-space caching improves complete-loop medians from
1.437 to 1.375 seconds in the original order and 0.492 to 0.458 seconds in the
early-rung order. Adding connection collection to caching has a small further
gain in a separate cohort: 1.336 to 1.315 seconds and 0.442 to 0.436 seconds.
Seven of eight original-order pairs and five of eight early-rung pairs improve
in the latter comparison. Connection collection alone was neutral or mixed;
fewer retained terms do not remove the cost of constructing their keys.
These are native typed operation loops in persistent processes, excluding
Python dispatch. All samples remain under
`/tmp/idenso-native-typed-compact-pipeline` and
`/tmp/idenso-native-typed-alpha-compact-pipeline`; the separate cohorts are not
pooled.

=== Complete Python comparison of the retained collector

Host `97833bbd` includes connection collection, compact-space caching and the
Scalar materializer correction above. It is compared with the unchanged
`63ba01ba` host using eight fresh process pairs per contraction order. Each
measured loop includes all eight public tensor-safe replacements, expansions,
contractions, interface checks and Python result wrapping. Imports, initial
construction, exact comparisons and final result destruction are outside the
clocks. All first samples and slower pairs are retained.

#table(
  columns: 5,
  table.header([Order], [Saved host (s)], [New host (s)],
    [Median paired change], [FORM process (s)]),
  [Original], [1.625594], [1.612381], [−3.48%], [0.814914],
  [Early-rung], [0.551142], [0.501827], [−9.00%], [0.189660],
)

The time columns are independent side medians; the change column is the median
of the eight candidate/baseline ratios. These are different statistics: the
original-order ratio of side medians improves only 0.81%. Seven original-order
pairs and six early-rung pairs improve in wall, process and thread time. The
slower wall-time pairs remain in the record: +19.09% in the original order,
and +52.13% and +28.25% in the early-rung order. The observations establish a
modest, mixed improvement rather than a uniform speedup.
Processes are pinned to CPU 8, but the shared host is not exclusive: recorded
activity on its SMT sibling ranges from zero to 97.4%. These observations do
not establish the cause of individual slower pairs; none are excluded.

FORM uses the same programs and contraction orders, with one declared warmup
and three measured fresh processes per order. Its wall clock includes startup,
parsing and teardown, unlike the Python operation loop. The displayed medians
put the new loop at 1.98 and 2.65 times the separate FORM process times; this is
not parity or a comparison of isolated contraction kernels.

Both fresh FORM polynomials equal the final 9,652-term result exactly. The
host also preserves 207 saved API cases and 117 HEP component rows. Four
independent finite-network evaluations compare the full factored ladder and
its scalar result under two exact momentum assignments. The actual notebook
passes all three routes in both orders, and 25 additional public controls
cover materializer opacity, connection identities, typed zeros and mixed
representation spaces. No correctness assertion was relaxed.

Three fresh process pairs also measure the unchanged trace first-call and
rerun controls. First-call medians are milliseconds:

#table(
  columns: 5,
  table.header([Trace case], [Saved host], [New host],
    [Median paired change], [FORM body CPU]),
  [Free 2, D], [0.02024], [0.02108], [+7.25%], [0.00066],
  [Free 4, D], [0.07223], [0.07918], [+11.14%], [0.00135],
  [Free 6, D], [0.12117], [0.12447], [+2.72%], [0.00540],
  [Free 10, D], [1.95912], [1.92722], [+22.89%], [0.35333],
  [Free 12, D], [23.81537], [22.30606], [−6.34%], [4.23333],
  [Free 8, 4D], [0.25163], [0.25116], [−0.19%], [0.04800],
  [Repeated 8, 4D], [0.09484], [0.09373], [−0.25%], [0.00200],
  [Repeated 8, D], [0.16779], [0.16449], [−0.64%], [0.01420],
  [Axial 12, 4D], [1.88230], [1.76866], [−9.08%], [0.62400],
)

Python includes the typed simplification method and factored scalar spectator.
FORM's internal timer covers only `trace4` or `tracen` plus sorting, excluding
that spectator, Python dispatch and interface handling. Its declared warmup
batches and amortization are unchanged. These unequal scopes must not be read
as isolated-kernel speed ratios.

The trace observations do not establish a general improvement. Free 6 in D is
slower in all three pairs. Free 10 in D has paired first-call ratios 1.229,
0.892 and 2.085, despite its slightly smaller independent side median; that
last pair's rerun ratio is also 1.811. Free 12 rerun medians are 9.804 and
9.769 milliseconds, and axial reruns are 3.053 and 3.068 milliseconds. All
individual measurements remain available; no trace regression is dismissed
because the source change targets the ladder.

Separate profiles of the new host reuse the existing whole-loop harness, at
199 Hz with DWARF call stacks. They retain 291 original-order samples and 94
early-rung samples. The existing disjoint owner classification attributes
27.6% / 23.7% of sampled cycles to Symbolica normalization outside other
classified owners, 27.1% / 31.3% to Schoonschip, 20.5% / 19.9% to certified
replacement, and 7.4% / 6.2% to interface proof. These single profiles locate
remaining work; especially the shorter profile is too small for precise
percentage comparisons or a speedup estimate.
Within the original-order normalization category, byte comparison contributes
6.1 percentage points of the whole loop, term merging 3.0 and sorting 2.6.
Those Symbolica operations remain useful optimization targets alongside the
shared tensor code.

The inclusive `normalization_is_intrinsic` check accounts for 9.5% / 7.9% of
sampled cycles, overlapping the owner groups above. Its observed callers are
replacement preflight and the outer contraction's typed-result validation,
with one original-order sample in local RHS validation. Caching or avoiding
redundant proofs is a concrete remaining opportunity, but merely trusting the
old interface would revive the callback-induced rank-loss counterexample.
The current implementation keeps those checks. Coefficient collection and
normalization also remain substantive costs even with FORM-matching retained
term counts.

The verified wheel is installed in `/tmp/spenso-paper-venv`, with the immutable
`63ba01ba` baseline preserved separately. The installed core equals the tested
candidate byte for byte; all 25 focused public smoke checks pass, and the 33
other installed dependency records remain unchanged.

The complete commands, source manifest, raw timings, CPU/SMT observations and
validation results are under `/tmp/idenso-component-optimization-host`.
The source snapshot has 765 inputs; all frozen identities remain unchanged
through the formal comparison. These results do not overwrite or pool the
earlier cohorts.

=== Reusing intrinsic normalization proofs

The next inference change caches successful callback-safety checks for whole
normalized function subtrees inside the existing `InterfaceInference` object.
The bounded cache stores at most 256 expressions of at most 256 bytes each.
It includes function arguments and scalar metadata in its keys; a shared head
alone is insufficient to certify the contents. Negative results and unfinished
normalization are not cached, and no persistent proof is attached to a mutable
`SymbolicTensor` payload.

The entire source must still pass the intrinsic-normalization check before
interface inference starts. That ordering prevents inference or materialization
of an earlier leaf from running before a later callback is discovered. Once
the source passes, the leaf proof can omit its duplicate metadata scan. The
algebra-only entry point retains its independent check. Callback-sensitive
results, unresolved ports, multiplicities and logical interface order retain
their existing validation, including the metric contraction whose callback
turns `T(b)` into a scalar.

Matched native experiments compare this change with the retained `97833bbd`
implementation using the same public Rust operations as the notebook. Eight
paired calls per order improve in every pair: the median paired wall-time
reductions are 7.31% in the original order and 6.60% with early rung closure.
Side medians are 1.315459 → 1.231453 s and 0.429324 → 0.404600 s. Process and
thread clocks agree in sign. These are persistent native workers, excluding
Python dispatch, and must not be pooled with fresh Python-process samples.
The alternative head-only cache gave smaller, mixed gains and is not included.

Both libraries reproduce every saved intermediate expression literally and
the same final FORM polynomial and fixed point. Separate timed FORM processes
check stage counts; the fresh full-polynomial checks belong to the Python host
validation. All observations and the distinct timing boundaries are retained
under `/tmp/idenso-native-typed-intrinsic-subtree-pipeline`.

The integrated source passes all 936 Rust tests, with 22 existing skips,
all-target Clippy and scoped formatting. Four added regressions cover repeated
subtrees, the global callback barrier, cache limits and raw/oversized inputs,
and equivalence with the previous two-phase proof under all inference policies.

=== Python measurement of proof reuse

The `ee756634` host differs from saved host `97833bbd` only in this inference
implementation. Eight fresh process pairs per order time the complete eight
public replacement, expansion and contraction steps, including Python wrapping.
Setup, initial construction, final equality checks and disposal remain outside
the clocks. The unchanged FORM programs use one declared warmup and three
measured fresh processes per order.

#table(
  columns: 5,
  table.header([Order], [Saved host (s)], [New host (s)],
    [Median paired change], [FORM process (s)]),
  [Original], [1.452938], [1.350243], [−7.15%], [0.856800],
  [Early-rung], [0.478071], [0.451845], [−3.70%], [0.177630],
)

The time columns are separate medians, whose ratios improve by 7.07% and 5.49%;
the paired column summarizes the eight individual ratios. Seven original-order
pairs and six early-rung pairs improve. The slower pairs are retained: +15.69%
original, and +5.98% / +47.11% early-rung. Process and thread clocks agree in
sign, so wall-clock preemption alone cannot explain those observations.
The shared machine was not reserved exclusively for these measurements.

The new Python loop is still 1.58 / 2.54 times the separate FORM process wall
time. FORM includes process startup, parsing and teardown; these are unequal
bounds, not an isolated-kernel comparison or parity claim. All 32 Python
results equal the same 9,652-term scalar exactly and are fixed points.
Retained intermediate term counts remain equal to FORM in both orders.

The candidate host also passes the unchanged 207 API cases, 117 HEP
component rows, four full factored-network evaluations, two fresh full FORM
polynomial comparisons, 58 existing public controls, 25 additional controls
and all three current notebook routes in both orders. No correctness assertion
or algebra fixture was changed. Source, core, driver and oracle identities
remain fixed throughout the measured cohort; the saved baseline stays intact.
Raw observations and the frozen protocol are under
`/tmp/idenso-intrinsic-proof-host`.

Three fresh process pairs measure the same trace first-call and rerun controls.
First-call side medians below are milliseconds; the paired change is the median
of the three candidate/baseline ratios.

#table(
  columns: 5,
  table.header([Trace case], [Saved host], [New host],
    [Median paired change], [FORM body CPU]),
  [Free 2, D], [0.02214], [0.02159], [−2.50%], [0.00066],
  [Free 4, D], [0.07397], [0.07440], [+2.68%], [0.00140],
  [Free 6, D], [0.12057], [0.13649], [+13.65%], [0.00560],
  [Free 10, D], [1.64681], [1.68647], [+0.73%], [0.37000],
  [Free 12, D], [20.67946], [20.71930], [+0.29%], [4.13333],
  [Free 8, 4D], [0.23977], [0.30042], [+24.10%], [0.04933],
  [Repeated 8, 4D], [0.09377], [0.11534], [+23.00%], [0.00200],
  [Repeated 8, D], [0.15485], [0.18306], [+19.39%], [0.01420],
  [Axial 12, 4D], [1.58340], [1.70693], [+4.29%], [0.64800],
)

These observations do not establish a general trace improvement. Free 6 in D,
free 8 in 4D and both repeated-index cases are slower in two of three first-call
pairs. Axial 12 is slower in all three. Repeated 8 reruns are also slower in all
three pairs in both dimensions, with paired median increases of 12.31% in 4D
and 1.62% in D. Free 8 in 4D reruns increase 18.04% by the paired median.
Free 12 in D rerun side medians are 9.402 → 9.439 ms; axial reruns are
2.935 → 2.899 ms, with a paired median change of −0.15%. Every inner sample
and process median remains in the record; the measurements do not identify
the cause of these regressions. The retained change prioritizes the complete
ladder, as requested.

FORM's trace timer covers only the trace body and sorting. Python additionally
handles the factored scalar spectator and public typed method. Those columns
are not measurements of identical stages. The fresh FORM outputs and all
Python trace results pass the unchanged exact checks.

Separate 199 Hz profiles retain 329 original-order and 84 early-rung operation
samples. The disjoint sampled cycle shares are 32.6% / 31.1% Schoonschip,
30.5% / 24.4% other Symbolica normalization, 16.8% / 20.0% certified replacement,
and 6.2% / 7.7% algebra/interface proof. The intrinsic callback scan contributes
2.5% / 3.3% inclusively and overlaps those groups. These small single-loop
profiles identify remaining work; their percentages are not paired speed
estimates. The early-rung record retains one empty stack and one unknown leaf.

The notebook still materializes its expanded expression before raw Schoonschip
contracts it. The separate network-local expansion routine is not on this call
path. Combining expansion with component collection must therefore receive the
unexpanded expression at the actual replacement/expansion boundary, while
preserving scalar spectators, powers, callbacks and tensor interfaces. That
combination is a remaining opportunity, not part of this measured change.

At that checkpoint, the verified `ee756634` wheel was installed in
`/tmp/spenso-paper-venv`. The immutable `97833bbd` baseline was preserved,
all 33 other dependency records were unchanged, and the installed host passed
the 25 focused public smoke cases.

=== Expansion needs early collection

An isolated native experiment feeds unexpanded factors into the existing
component collector. All intermediate tensor graphs agree with the saved
pipeline after canonicalization, and both orders produce the exact 9,652-term
FORM scalar and a fixed point. Nevertheless, the first implementation is
slower. Four counterbalanced fresh processes per route and order give these
complete-loop median seconds:

#table(
  columns: 4,
  table.header([Order], [Saved expanded route], [Fused intake],
    [Same new build, expanded route]),
  [Original], [1.34406], [3.00912], [1.35602],
  [Early rungs], [0.45983], [0.53285], [0.48617],
)

The fused route loses every pair against both controls. Relative to the saved
route, its median paired wall-time ratios are 2.251 and 1.175; process and
thread CPU clocks agree closely. All 24 results pass the exact checks, and no
sample is discarded. These are native screening measurements, separate from
the Python and FORM timings above. The complete fused interval includes its
required final typed expansion; fixture setup and final checks are outside.
The experiment is not installed.

Counters on the actual original-order inputs explain the largest regression.
At stage 6, the fused intake plans 65,456 monomials; normal expansion first
combines the same input into 18,247 terms. At stage 7 it plans 121,158 rows
with 333 variables, requiring 80,691,228 bytes of dense exponents. That exceeds
the existing 64 MiB bound, so the complete attempt is discarded. Ordinary
contraction then runs, followed by another attempt with 120,108 rows and 213
variables. That retry fits and succeeds. Normal expansion of the original
stage-7 input has only 34,196 terms.

Thus avoiding expanded Atom construction also requires preserving expansion's
early combination and cancellation of equal monomials. Performing tensor
graph work before that collection multiplies the work even when no limit is
hit; exceeding the bound adds fallback and repeated planning. The diagnostic
counters run outside the timing cohort. Sources, all observations, canonical
stage checks and the admission records are preserved under
`/tmp/idenso-fused-component-expansion`. The retained implementation and
notebook remain unchanged by this first experiment.

=== Collect before planning contractions

The selected implementation extends the existing component collector behind
`expand_contracted_sums=True`. It compiles admitted sums, products and positive
integer powers into a bounded syntax tape, interns their borrowed factors once,
and collects equal factor/exponent keys with exact rational coefficients before
building contraction graphs. Cancelled terms never reach incidence checks.
Common scalar spectators remain factored; scalar differences within the indexed
core are distributed. Unsupported inputs use the existing ordinary traversal.
The final polynomial emitter and its existing memory bound are unchanged.

Local held right-hand sides still use Symbolica expansion before contraction.
The notebook removes only the ambient expansion before Schoonschip, then
expands the complete final scalar inside the last measured step. Eliminating
both expansions was slower in the native screen. This distinction matters:
the small local expansion is useful, while constructing the complete indexed
expression can be avoided.

Four counterbalanced fresh native processes per route and order compare the
saved implementation, the new build with its old expanded schedule, fully
fused intake, and the selected schedule. Complete-loop side medians are:

#table(
  columns: 5,
  table.header([Order], [Saved expanded (s)], [New expanded (s)],
    [Fully fused (s)], [Selected (s)]),
  [Original], [1.25712], [1.27572], [1.37565], [1.15887],
  [Early rungs], [0.416078], [0.428726], [0.435531], [0.380686],
)

The selected schedule improves every pair against both expanded controls.
Median paired reductions against the saved implementation are 7.44% and 8.60%;
against the same new build they are 8.21% and 10.91%. Fully fused intake is
slower than the saved implementation in every pair. All 32 observations are
retained, CPU clocks agree in direction, and all final results equal the exact
9,652-term FORM polynomial. Intermediate results agree by tensor graph
canonicalization when factorization or dummy-index spelling differs. These
native measurements exclude Python dispatch and are separate from the public
Python comparison.

The shared Idenso implementation retains typed interfaces and callback checks;
Spynso adds no second contraction implementation. Fourteen regressions cover
the admitted factored inputs, exact coefficient collection, cancellation before
incidence validation, unsupported surviving terms, spectators, callbacks and
typed zero. An empty collected coefficient list returns zero before entering
Symbolica's polynomial constructor. The full affected Rust suite passes
950 tests, with 22 existing skips; Clippy and formatting pass.

=== Public Python measurement of collected intake

Saved host `ee756634` and candidate `58dd008a` share their dependency snapshot;
the candidate changes the four Schoonschip implementation/settings files.
Eight fresh process pairs per ladder order compare the previous expanded
schedule with the selected schedule above. Both include all eight public
typed replacements, local expansion and contractions; the candidate also
includes its final typed expansion. Initial construction, term counting,
equality checks and disposal are outside the clocks.

#table(
  columns: 5,
  table.header([Order], [Saved host (s)], [Selected host (s)],
    [Median paired change], [FORM process (s)]),
  [Original], [1.357979], [1.234922], [−10.66%], [0.856578],
  [Early rungs], [0.453248], [0.409100], [−8.64%], [0.179846],
)

Every wall, process and thread clock pair improves in both orders. The ratios
of the side medians improve by 9.06% and 9.74%; these differ from the medians
of paired ratios shown in the table. All 32 timed results equal the same
9,652-term scalar exactly and are fixed points. No observations are excluded.
The shared machine is not reserved exclusively for these measurements.

The remaining ratios to FORM are 1.44 and 2.27. FORM includes process startup,
parsing and teardown, while the Python clock covers the prepared public loop.
These unequal bounds do not establish kernel parity. Fresh untimed full FORM
polynomials certify both orders independently of the timed term-count checks.

The candidate passes 207 API cases, 117 HEP component rows, four network
evaluations comparing the original factored ladder and certified scalar at
two assignments, 58 public and tensor-safe controls, 25 additional controls,
and all three actual notebook routes in both orders.
The control records match the saved host, including callbacks, errors and
logical port ordering. Source, extension and driver identities stay unchanged
through the complete cohort.

The unchanged trace controls use three fresh process pairs per case. First-call
side medians below are milliseconds; FORM measures its internal trace-and-sort
CPU body, excluding the Python scalar spectator and public wrapper.

#table(
  columns: 5,
  table.header([Trace case], [Saved host], [Selected host],
    [Median paired change], [FORM body CPU]),
  [Free 2, D], [0.02065], [0.02002], [−4.38%], [0.00064],
  [Free 4, D], [0.07297], [0.07106], [−2.62%], [0.00140],
  [Free 6, D], [0.12761], [0.12138], [−4.36%], [0.00540],
  [Free 10, D], [1.63196], [1.64530], [+7.12%], [0.35333],
  [Free 12, D], [19.86224], [20.10018], [+2.23%], [4.13333],
  [Free 8, 4D], [0.24366], [0.23789], [−2.49%], [0.05067],
  [Repeated 8, 4D], [0.09477], [0.09354], [−1.30%], [0.00200],
  [Repeated 8, D], [0.15984], [0.15822], [−1.01%], [0.01440],
  [Axial 12, 4D], [1.57719], [1.61152], [+1.82%], [0.60400],
)

This is not a general trace improvement. Free 10 and 12 in D and axial 12 are
slower in two of three first-call pairs. Seven of nine rerun controls have
positive median paired changes: free 6/10/12 in D (+0.61% / +0.98% / +1.96%),
free 8 in 4D (+0.43%), repeated 8 in 4D/D (+5.11% / +9.41%), and axial 12
(+0.77%). Both repeated-index reruns are slower in all three pairs. Free 12 in
D rerun side medians are 9.279 → 10.332 ms; axial reruns are 2.878 → 2.900 ms.
All inner samples, process medians and slower observations are retained.
Several processes show broad slowdowns across cases, so these measurements
alone do not identify a code-level cause for the trace changes.

Separate profiles of one native selected-schedule loop per order retain
238 / 71 operation samples at 199 Hz. The disjoint sampled cycle shares are
26.9% / 25.4% Symbolica normalization, 18.1% / 19.6% component collection,
17.4% / 20.5% interface inference and proofs, and 16.4% / 15.7% other
Schoonschip work. Symbolica matching/replacement accounts for 7.6% / 5.9%,
polynomial operations for 6.8% / 4.5%, and expansion for 5.1% / 5.5%.
These shares point to several remaining costs rather than a single matcher
bottleneck. Inclusive collector and replacement stacks overlap these groups
and must not be added to them. There are no empty operation stacks and
8 / 1 unknown leaves. The small profiles exclude Python callback overhead;
they are diagnostic attribution, not another speed comparison or precise
wall-time breakdown. Attribution distinguishes the `normalize_dots` module
from Symbolica normalization; the corrected overlay and raw stacks are both
retained.

Untimed counters on the selected original-order route show the avoided work.
Stage 6 generates 43,824 branches, combines them into 19,291 distinct keys,
and sends only 18,247 nonzero terms to contraction planning. Stage 7 generates
81,099 branches, combines them into 35,250 keys, and plans 34,181 nonzero terms.
Both planned counts equal ordinary Symbolica expansion of that same input.
The stage-7 dense exponent estimate is 23,174,718 bytes with 339 variables,
within the existing 64 MiB limit. Both stages complete one successful ambient
collector call, with no rejected attempt or fallback/retry. These inputs use
the retained local RHS expansion, so their generated counts are distinct from
the earlier fully fused experiment.

The frozen protocol, complete observations and exact checks are under
`/tmp/idenso-fused-component-host`; native experiments remain separate under
`/tmp/idenso-fused-component-expansion`. The source-controlled historical
contraction archive is unchanged.

The validated `58dd008a` wheel is installed in `/tmp/spenso-paper-venv`.
The immutable `ee756634` baseline remains available, all 33 other dependency
records are unchanged, and the installed extension passes the 25 focused
public smoke cases.

=== Final scalar expansion choices

A caller-level audit of the same saved profiles places 23.09% / 23.24% of
sampled cycles inside final typed expansion, including 15.89% / 14.91% in
Symbolica normalization. Local RHS construction contributes only
2.13% / 2.77% inclusively. The audit also separates three / one generic helper
frames whose type names mention `TensorInferenceError` from concrete inference
frames. Concrete inference frames occur in 16.00% / 19.02% of sampled cycles;
the earlier broad subsystem grouping must not be read as that exact quantity.
All these inclusive caller shares overlap and retain the small-profile limits.

The last vertex replacement leaves 23,714 / 10,516 scalar terms containing
small factored dot sums. No component-collector call occurs at this last stage:
the repeated-index scan correctly finds no remaining contraction work. Final
expansion must therefore distribute actual scalar products, not merely sort an
already expanded polynomial. Expanding arbitrary scalar sums automatically in
Schoonschip would change its existing no-work behavior.

The existing public expansion choices were compared on these exact saved
inputs using installed host `58dd008a`. Three fresh processes per route and
order give the following median final-materialization seconds:

#table(
  columns: 4,
  table.header([Order], [Ordinary expansion], [`via_poly=True`],
    [Whole polynomial conversion]),
  [Original], [0.296172], [0.741699], [6.402327],
  [Early rungs], [0.096620], [0.240256], [2.602664],
)

The whole-conversion route includes both
`source.to_expression().to_polynomial().to_expression()` and construction of
the resulting `TensorExpression`; it does not bypass result validation.
Setup, parsing the saved input, equality checks and disposal are outside the
clocks. All 18 observations produce the exact 9,652-term scalar, preserve the
receiver and empty logical interface, and pass expansion and Schoonschip
fixed-point checks. No timeout, repeat or sample exclusion occurs.

Both alternatives lose every pair in both orders. Median paired wall-time
ratios are 2.50 / 2.52 for `via_poly=True` and 21.62 / 26.94 for whole
conversion; process and thread clocks agree. Ordinary expansion stays in the
notebook. The prepared complete-loop comparison is left unexecuted because
the only changed stage already loses consistently; no complete-loop gain is
claimed for either alternative.

The pinned Symbolica implementation expands each Add term separately in
`expand_via_poly`, then normalizes the full Atom sum. Whole polynomial
conversion instead accumulates converted Add terms through sequential
polynomial additions. The source suggests repeated accumulation as a cost to
isolate; these timings alone do not attribute the whole conversion slowdown
to that loop. The fixed protocol, all observations and rejected-route decision
are preserved under `/tmp/idenso-final-scalar-expansion`.

An isolated Symbolica-only diagnostic separates conversion from emission and
also tries balanced polynomial accumulation. It parses the same saved scalar
text with functions opaque; it does not register tensor tags or callbacks, so
its clocks are separate from the typed measurements above. Three fresh
processes per route and input give these raw-operation medians:

#table(
  columns: 4,
  table.header([Input], [Ordinary (s)], [Sequential polynomial (s)],
    [Balanced polynomial (s)]),
  [Original ladder scalar], [0.250168], [6.795221], [0.881246],
  [Early-rung ladder scalar], [0.077039], [2.561941], [0.296572],
  [256-term synthetic], [0.002533], [0.423715], [0.006718],
)

For the original ladder, sequential conversion takes 6.786 s and emission
9.32 ms; balanced conversion takes 0.871 s and emission 10.35 ms. The
early-rung equivalents are 2.551 s / 11.05 ms and 0.286 s / 10.82 ms.
Conversion dominates this whole-polynomial route on the pinned, unpatched
Symbolica dependency. Balancing improves that route but still loses every pair
against ordinary expansion: median paired ratios are 3.55, 3.85 and 2.47.
All 27 results equal ordinary expansion of the identical parsed input and
pass its fixed-point check; neither route is adopted in the notebook.

This diagnostic changes both accumulation order and the variable map supplied
to each term conversion. It therefore does not isolate addition-tree shape as
the sole cause. Its synthetic case also grows the variable basis with input
size. These limits and all observations are retained under
`/tmp/idenso-final-scalar-expansion/balanced-poly`.

The standalone
#link("../../../../examples/reproducers/symbolica-expansion/performance.typ")[scalar
expansion reproducer] makes these four routes available without tensor-library
dependencies. Its synthetic input keeps twenty variables fixed while changing
the number of terms; it also accepts a saved expression file. All four routes
pass exact ordinary-expansion and fixed-point checks at 20, 256, 1,024 and
4,096 input terms and on both saved ladder scalars. These 24 checks validate
the reproducer; their single-call clocks are not a new benchmark cohort.

=== Small admission experiments

Two isolated native candidates were compared against a fresh build of the
installed source and dependency family. The unchanged selected schedule
includes final typed expansion inside its clocks. Six permutations of the
three variants, for each ladder order, retain 36 fresh-process observations.
Every candidate stage agrees literally with the corresponding baseline, and
the final scalar has 9,652 terms. These native clocks are not pooled with the
public Python measurements.

Caching exact positive-power interface proofs gives a median paired wall-time
increase of 0.59% in both orders, with three of six pairs faster in each.
That candidate is rejected. Checking a function's arity and tensor tags
before walking its arguments for vector admission reduces paired medians by
1.53% / 2.04%; five / four of six wall pairs improve. This is a modest
candidate, with slower observations retained, and has not been installed.
The full protocol, source identities and results are under
`/tmp/idenso-selected-hybrid-native`.

=== Reuse scalar collection for final expansion

The retained implementation moves expansion dispatch into the shared
`SymbolicTensor<PartialStructure>` abstraction. Spynso converts arguments and
wraps the checked result. The default expansion mode can reuse the existing
bounded coefficient collector for scalar sums, while selected-variable and
`via_poly` operations retain their Symbolica dispatch.

This private intake admits exact rational arithmetic, ordinary scalar variables
and compatible compact metric dots. It checks every source leaf before
cancellation, retains normalized dot atoms literally, and skips all tensor
graph and index-contraction work. Explicit tensor indices, unresolved ports,
user normalization hooks, unsupported powers and resource-limit failures use
the existing expansion and callback-aware finisher. Successful admission
already proves that the established interface is preserved, so it avoids a
second proof scan. Unchanged results retain their flags and metadata; typed
zeros retain their interface. Schoonschip's scalar no-work behavior is unchanged.

The collector builds exact coefficients and exponent lists directly; it does
not call general expression-to-polynomial conversion on the input. Output
still uses the existing coefficient-list polynomial emitter. The improvement
comes from collecting equal monomials before constructing the resulting Atom
sum and from reusing the scalar-interface proof.

A separate diagnostic build confirms successful admission of both actual final
ladder inputs, with no declined attempt. The original order generates 141,158
virtual products, collects 10,121 distinct monomial keys and removes zero
coefficients to emit 9,652 terms. Early rungs generate 47,416 products, collect
10,027 keys and emit the same 9,652 terms. Both use 26 literal scalar leaves
and 501,904 bytes of dense exponents. These are collector counts before Atom
materialization, not FORM's intermediate term statistics.

Four counterbalanced fresh-process pairs per order compare the exact saved
native build with the candidate. A separate cohort measures only the final
typed expansion; these are distinct clock scopes, each with 16 observations.

#table(
  columns: 5,
  table.header([Order], [Complete loop, saved (s)], [Candidate (s)],
    [Final expansion, saved (ms)], [Candidate (ms)]),
  [Original], [1.228201], [1.064481], [275.157], [97.765],
  [Early rungs], [0.422248], [0.368230], [82.924], [44.285],
)

Every pair improves on wall, process and thread clocks in both scopes. Median
paired complete-loop wall changes are −17.38% / −8.89%; final-only changes
are −64.73% / −46.82%. The complete-loop observations vary substantially:
original pair improvements range from 11.73% to 21.54%, and early-rung
improvements from 4.48% to 29.79%. No observation is excluded or repeated.
All 32 timed results agree exactly; all 16 intermediate stage expressions
match the baseline literally, and both final scalars equal FORM and remain
fixed points. The diagnostic binary is never used for timing.

The implementation passes 97 focused tests and the complete affected
Idenso/Spenso/Spynso suite in the actual worktree: 955 tests, with 22 existing
skips. All-target Clippy and formatting pass. Five new tests cover admitted
scalar sums, refusal before cancellation, exact fallback, dispatch options,
interface order, typed zero and callback behavior. The existing Python
metadata/zero regression also covers shared expansion. The current locked
dependencies and unrelated model, API and display changes are preserved.
Native evidence is under `/tmp/idenso-scalar-expansion-intake` and
`/tmp/idenso-scalar-expansion-native-harness/final`.

=== Public measurement of shared scalar expansion

Saved host `58dd008a` and selected host `9fb006cb` use the same dependency
snapshot and the same contraction schedule. Both use tensor-safe replacement,
local RHS expansion and fused ambient contraction; both include final typed
expansion inside the loop clock. The selected implementation changes five
shared tensor/Schoonschip and Python-dispatch source files. Symbolica remains
the pristine pinned `06906976` revision; isolated emission patches are not
included in either host.

Eight counterbalanced fresh-process pairs per order give these wall medians:

#table(
  columns: 5,
  table.header([Order], [Saved host (s)], [Selected host (s)],
    [Median paired change], [FORM process (s)]),
  [Original], [1.308102], [1.105130], [−16.41%], [0.851635],
  [Early rungs], [0.435243], [0.365242], [−14.60%], [0.220752],
)

All sixteen pairs improve on wall, process and thread clocks. The ratios of
the side medians improve by 15.52% / 16.08%, distinct from the median paired
changes in the table. Original-order observations range from 1.232–1.606 s
for the saved host and 1.033–1.169 s for the selected host; early-rung ranges
are 0.407–0.453 s and 0.358–0.410 s. No slower observation is removed and no
outcome-dependent repeat is performed. The shared machine is not reserved.

Against the same fresh FORM processes, the original-order ratio falls from
1.54 to 1.30 and the early-rung ratio from 1.97 to 1.65. These are unequal
clock boundaries: Python measures the prepared public eight-vertex loop,
excluding initial construction, counting, validation and disposal; FORM
measures startup, parsing, computation and teardown. They quantify the
observed remaining gap, not equal-scope engine parity. Each FORM order keeps
three measured processes after one designated warmup; raw warmups are saved.

Already expanded results now return without generic expansion. A separate
three-pair control on the actual 9,652-term final scalar improves expansion
rerun medians from 3.888 to 0.759 ms and from 3.548 to 0.878 ms. All six pairs
improve on all three clocks, including a retained 1.653 ms original-order
candidate observation. These controls are separate from the complete loop.

The trace controls below keep three fresh paired process medians per case,
each using five calibrated inner wall batches. Times are milliseconds; the
FORM column is the median of its retained internal trace-and-sort CPU
batches. It excludes the Python scalar spectator and wrapper, so those
columns also have different boundaries.

#table(
  columns: 5,
  table.header([Trace case], [Saved host], [Selected host],
    [Median paired change], [FORM body CPU]),
  [Free 2, D], [0.02213], [0.02171], [−7.45%], [0.00084],
  [Free 4, D], [0.07386], [0.07494], [+1.46%], [0.00170],
  [Free 6, D], [0.13121], [0.12079], [−3.19%], [0.00580],
  [Free 10, D], [1.61352], [1.75478], [+10.67%], [0.42667],
  [Free 12, D], [21.74959], [19.98120], [−6.61%], [4.96667],
  [Free 8, 4D], [0.23674], [0.26900], [+1.76%], [0.05467],
  [Repeated 8, 4D], [0.09309], [0.09647], [+3.30%], [0.00240],
  [Repeated 8, D], [0.15499], [0.15959], [+2.81%], [0.01740],
  [Axial 12, 4D], [1.59002], [1.62686], [+1.47%], [0.74800],
)

There is no general trace improvement. Free 10 in D is slower in all three
first-call pairs. Free 8 in 4D reruns are also slower in all three pairs
(median paired +8.17%, side medians 0.06001 → 0.06981 ms). Free 12 in D and
axial reruns have positive paired medians of +0.84% and +2.64%; their side
medians are 9.714 → 9.796 ms and 2.920 → 2.954 ms. The other six rerun
controls have negative median paired changes. All trace samples and slower
observations remain in the report; this cohort does not isolate the cause of
those changes.

A subsequent source audit follows the free-10 D control through
`simplify_gamma(expand_traces=True)`. Its ten distinct indices use the existing
generic terminal trace evaluator; the scalar spectator remains factored and
the evaluator returns before Schoonschip component collection. The Python
wrapper uses the existing rewritten-result finisher. Neither the new public
expansion method nor scalar intake is reached, and the trace kernel and
finisher bodies are byte-identical in the saved and selected sources. This
excludes direct execution of the new dispatch as the explanation for the
observed regression; effects from code generation, layout or allocation have
not been isolated. The source comparison and route evidence are recorded in
`/tmp/idenso-contraction-next/free10-dispatch-audit.json`.

All measured ladder results equal the same scalar exactly. Separate untimed
validation passes 207 API records, 117 HEP component rows, four network
evaluations comparing the original factored ladder and scalar at two
assignments, fresh complete FORM polynomials in both orders, 83 public and
tensor-safe controls, all three actual notebook routes in both orders, and
six installed metadata tests. Logical port order, callbacks, errors, typed
zero and the preserved input are checked independently of the timing loops.

The frozen protocol and all observations are under
`/tmp/idenso-scalar-expansion-host`; `formal-summary.json` records the qualified
comparison. Source and extension identities remain unchanged throughout the
cohort. The source-controlled historical contraction archive is unchanged.

At the end of this cohort, the selected `9fb006cb` wheel was installed in
`/tmp/spenso-paper-venv` and passed the 25 focused public controls and six
metadata tests there. A fresh
before/after inventory preserves all 28 other local distributions and five
additional distributions exposed by the development environment. The saved
`58dd008a` host remains unchanged. Installation records are separate from the
pre-installation benchmark summary under `canonical-install`.

=== Remaining cost after shared expansion

Separate native profiles of the `9fb006cb` implementation retain 783 / 272
operation samples from four complete loops per order at 199 Hz. All eight
loops produce the exact FORM scalar and pass fixed-point checks. The disjoint
sampled-cycle shares are:

#table(
  columns: 3,
  table.header([Owner], [Original], [Early rungs]),
  [Component intake and collection], [28.24%], [29.51%],
  [Other Schoonschip work], [17.21%], [16.60%],
  [Interface inference and proofs], [17.84%], [18.90%],
  [Symbolica normalization], [16.17%], [17.22%],
  [Symbolica polynomial operations], [9.59%], [6.15%],
  [Symbolica matching/replacement], [5.46%], [6.73%],
  [Typed replacement bookkeeping], [3.65%], [2.30%],
  [Index multiplicity and composition], [1.61%], [2.21%],
)

Final shared expansion accounts for 9.63% / 13.37% inclusively. Most sampled
Symbolica normalization now occurs while rebuilding sums/products after
tensor-factor replacement: 11.18% / 10.31% of the whole loop. This is where
changed branches are materialized into normalized Atoms. The remaining work
is spread across coefficient collection, tensor scans/proofs and Atom
construction; matcher search alone does not explain the gap. The inclusive
replacement and expansion shares overlap the table and must not be added to
its disjoint owners.

Intrinsic-normalization scans account for 9.88% / 8.70% inclusively across all
callers. Of that, 3.58% / 4.13% of the whole loop is the final scalar intake's
admission check before expansion. It is not reinference of the emitted
result. Ambient collector admission, Schoonschip finishing and replacement
preflight account for the other observed intrinsic scans.

These are sampled native cycle shares, not wall-time fractions or another
speed comparison. Four loops in one process may reuse caches differently
from the fresh-process benchmark; Python callback crossings are absent.
Inlining and truncated stacks limit attribution, with ten / one unknown
leaf samples and no empty retained stacks. Classification uses the actual
callee, excluding misleading generic parameter names such as
`normalize_dots` and `TensorInferenceError`. The unchanged executable,
sources and loaded runtime dependencies are verified. Raw profiles and
corrected caller reports are under
`/tmp/idenso-scalar-expansion-native-harness/final/selected-scalar-profile/v2`;
`interpretation-overlay.json` records the corrected collector name and
distinguishes scalar-intake admission from result finishing.

=== Certify each input leaf once

The collector now proves intrinsic normalization while compiling each distinct
input leaf, using its existing literal-leaf table. Previously it walked the
whole input first, repeating the same function-head and metadata checks for
every occurrence. Factored contraction still checks each original opaque
function subtree before caching it. The stricter scalar grammar already
certifies every admitted compact dot head and operand. Root scalar factors
that stay factored retain the existing complete/intrinsic observer check.

Compilation of the entire input still precedes distribution and cancellation.
A canceled term cannot conceal a user callback, and unsupported powers or
resource exhaustion decline the intake before constructing results. Callback
declines still use the checked ordinary path, including the counterexample
where contracting a metric turns a tensor into a scalar. This change reuses
the same tensor abstraction, collector and emitter; it adds no public mode
or persistent cache.

A separate allocation experiment reused a temporary monomial-key buffer.
Both candidates were compared with a freshly built baseline from the selected
`9fb006cb` source/dependency family. Each scope retains six permutations of
the three variants for both orders: 36 complete-loop observations and 36
separate final-expansion observations. Native clocks exclude Python wrapping
and are not pooled with the public measurements above.

#table(
  columns: 5,
  table.header([Order], [Complete baseline (s)], [Leaf proof (s)],
    [Final baseline (ms)], [Leaf proof (ms)]),
  [Original], [1.042415], [0.978628], [101.135], [59.778],
  [Early rungs], [0.353807], [0.333711], [44.320], [32.292],
)

Median paired complete-loop changes are −5.59% / −5.32%, with five of six
pairs faster in each order. Process and thread clocks agree. Final-only
changes are −38.72% / −27.12%, with six / four of six pairs faster. All slower
observations are retained, including two early-rung final-expansion pairs
that regress slightly. The complete-loop screen and its selection criteria
were fixed before timing; there are no outcome-dependent repeats.

The temporary-key experiment fails those criteria: complete-loop paired
changes are +0.59% / −1.72%, with only two / three of six pairs faster.
Final-only paired changes are also positive (+0.29% / +1.24%). It is not
retained. Its reduction in temporary allocation opportunities does not
establish a runtime improvement, and the prepared combined candidate remains
unbuilt.

All 72 timed results pass exact and fixed-point checks. Untimed gates compare
all sixteen intermediate expressions literally with the baseline, both final
9,652-term scalars with FORM, and two saved scalar inputs plus three fixed-basis
synthetic inputs. The leaf-proof implementation also passes 100 focused tests,
Clippy and formatting. Its three new regressions use independent ordinary
expansion/Schoonschip result oracles, check canceled and hidden callbacks,
and preserve metric-induced rank-loss validation and typed zero. The selected
source is integrated and passes the complete affected Idenso/Spenso/Spynso
worktree suite: 958 tests, 22 existing skips, all-target Clippy and formatting.
The original locked dependencies and all 677 checked source inputs remain
unchanged through validation.
The fixed dependency snapshot, protocols and all observations are under
`/tmp/idenso-contraction-next`, with source and validation records under
`/tmp/idenso-scalar-admission-proof`.

==== Public comparison after leaf admission

The public host changes only the collector admission source relative to the
saved `9fb006cb` host, using the same pristine Symbolica `06906976` dependency.
Both sides perform the same complete typed ladder loop, including final
expansion. The prepared input, result checks, term counts and disposal are
outside the clocks. Eight fresh paired calls per order retain all wall,
process and thread observations. The selection criterion was fixed before
timing: negative median paired changes on all three clocks and at least five
of eight faster wall pairs in each order.

#table(
  columns: 6,
  table.header([Order], [Saved host (s)], [Leaf proof (s)],
    [Paired wall change], [FORM process (s)], [New / FORM]),
  [Original], [1.138677], [1.037415], [−6.57%], [0.829428], [1.251×],
  [Early rungs], [0.383812], [0.350274], [−12.71%], [0.197007], [1.778×],
)

Seven of eight wall pairs improve in each order, and process/thread changes
agree. The candidate therefore passes the retention criterion. Ratios of
the separate side medians give −8.89% / −8.74%; these are different statistics
from the median paired changes in the table. Against the same fresh FORM
observations, the saved host ratios are 1.373× / 1.948×. FORM's whole-process
wall clock includes startup, parsing and teardown, whereas Python's prepared
loop excludes input construction. These ratios describe the measured
boundaries, not equal-scope engine overhead or parity. Native screens and
historical public cohorts are kept separate.

The separate final-expansion fixed-point controls remain mixed. Original-order
side medians rise from 1.068 to 1.497 ms even though the median paired change
is −7.36%, with two of three pairs faster. Early-rung medians rise from 0.840
to 0.860 ms, with a +6.56% paired change and only one of three pairs faster.
All observations are retained; this is not evidence of a rerun improvement.

The trace controls use three fresh paired processes per case and five
calibrated wall batches per process. The following first-call medians are in
milliseconds. FORM reports internal trace-and-sort CPU time, excluding the
Python scalar spectator and wrapping, and therefore has a different scope.

#table(
  columns: 5,
  table.header([Trace case], [Saved host], [Leaf proof],
    [Paired change], [FORM body CPU]),
  [Free 2, D], [0.02048], [0.02032], [−1.58%], [0.00066],
  [Free 4, D], [0.07232], [0.07189], [−0.97%], [0.00135],
  [Free 6, D], [0.11864], [0.12108], [+1.61%], [0.00520],
  [Free 10, D], [1.69740], [1.56434], [−4.73%], [0.35000],
  [Free 12, D], [22.07120], [20.12071], [−8.50%], [3.90000],
  [Free 8, 4D], [0.24780], [0.23661], [−5.22%], [0.04933],
  [Repeated 8, 4D], [0.09323], [0.09308], [−0.61%], [0.00220],
  [Repeated 8, D], [0.15511], [0.15609], [−0.05%], [0.01460],
  [Axial 12, 4D], [1.77805], [1.68551], [−4.57%], [0.64400],
)

Free 6 in D is slower in all three first-call pairs. Rerun controls below also
retain the positive paired changes for free 4 and free 6 in D. Times are
milliseconds and remain independent of the first-call samples.

#table(
  columns: 4,
  table.header([Trace rerun], [Saved host], [Leaf proof], [Paired change]),
  [Free 2, D], [0.001101], [0.001062], [−3.58%],
  [Free 4, D], [0.001872], [0.001882], [+0.85%],
  [Free 6, D], [0.007474], [0.007455], [+0.81%],
  [Free 10, D], [0.683301], [0.682330], [−0.14%],
  [Free 12, D], [9.810806], [9.418374], [−5.93%],
  [Free 8, 4D], [0.059934], [0.059437], [−0.83%],
  [Repeated 8, 4D], [0.002364], [0.002341], [−0.98%],
  [Repeated 8, D], [0.012684], [0.012618], [−0.52%],
  [Axial 12, 4D], [3.013321], [3.003201], [−2.19%],
)

There is no general trace-speedup claim. In particular, the free-10 route
audited above does not execute the changed admission code, so this cohort's
improvement does not establish a causal benefit from that change. The earlier
`9fb006cb` sampled profiles identify the motivation for the optimization;
their percentages are not reused as measurements of the new implementation.
The separate retained-host profile below measures the remaining work.

Untimed public validation passes 207 API records, 117 HEP component rows,
four network evaluations comparing the original factored ladder and scalar
at two assignments, both fresh complete FORM polynomials, 83 public and
tensor-safe controls, all three notebook routes in both orders, and six
metadata tests. All timed results pass their exact checks. No source or
extension identity changes during the cohort, and no observations are
discarded or repeated based on their outcome.

The selected host is `d58aee8a`. The frozen protocol, raw observations and
qualified report are under `/tmp/idenso-scalar-admission-host`, with
`formal-summary.json` recording the completed comparison before installation.
The historical contraction archive remains unchanged.

The selected wheel is installed in `/tmp/spenso-paper-venv`. Its 25 focused
public control records exactly match the isolated candidate, and six installed
metadata tests pass. A fresh inventory verifies that all 28 other local
distributions and five additional development-environment distributions are
unchanged. The immutable `9fb006cb` baseline remains available. The separate
`canonical-install/summary.json` records the wheel, installed extension,
dependency inventory and checks.

=== Factor recognition and terminal storage experiments

Two isolated candidates test additional collector costs against the retained
`d58aee8a` source. The first caches immutable factor recognition by input ID
after exact coefficient cancellation. It reuses resolved metric/vector
endpoints and opaque tensor ports while retaining exponent checks, fresh
incidence and contraction planning per monomial. Optional recognition storage
is bounded; exhaustion repeats ordinary recognition. Original callback
admission still precedes distribution and cancellation.

Here the original vertex order is `[1,2,3,4,5,6,7,8]`; "early rungs"
(called `reverse` by the measurement driver) is `[5,4,6,3,7,2,8,1]`.
It already closes successive triangles, with expanded-vertex frontier ranks
`[2,3,2,3,2,3,2,0]`. The opposite traversal `[1,2,8,3,7,4,6,5]` was measured
historically on a growing-accumulator implementation, not on this current
all-opaque loop. Its old times do not establish a new improvement here;
`exact-order-audit.json` under the native measurement directory records that
distinction and the exact arrays.

The second replaces each node's terminal vector with two inline entries.
Accepted degree-at-most-two graphs consist of paths and cycles, which have
two or zero component terminals after free endpoints are included. Checked
insertion and ordered merging preserve the existing behavior. This removes
terminal allocations but makes each node larger; fewer allocations alone do
not establish a runtime benefit.

A fixed native comparison retains all six variant permutations twice per
order: 72 fresh processes, twelve pairs per candidate and order. All variants
include the final shared expansion. Before timing, retention requires median
paired ratios below one on wall, process and thread clocks in both orders,
plus at least eight of twelve faster wall pairs per order.

#table(
  columns: 4,
  table.header([Complete loop], [Baseline (s)], [Recognition (s)], [Inline terminals (s)]),
  [Original], [0.970406], [0.947283], [0.985330],
  [Early rungs], [0.346926], [0.317414], [0.342484],
)

#table(
  columns: 5,
  table.header([Candidate], [Original paired change], [Faster pairs],
    [Early-rung paired change], [Faster pairs]),
  [Recognition], [−3.33%], [7/12], [−6.93%], [9/12],
  [Inline terminals], [+2.11%], [4/12], [−3.47%], [7/12],
)

Neither candidate qualifies. Process and thread paired medians have the same
directions as wall time. No observations are excluded or repeated, including
the recognition candidate's slower ratios of 1.315 and 1.472 in the two
orders. The separate side medians in the first table and median paired changes
in the second are distinct statistics. These native measurements omit Python
wrapping and are not pooled with the earlier public/FORM cohort. No combined
candidate or new public host is built; production and the installed `d58aee8a`
host retain the previous implementation.

Both candidates pass all sixteen literal ladder-stage comparisons, both
complete 9,652-term FORM scalar checks, fixed points and five scalar controls.
Inline storage passes 100 focused tests; recognition passes 103. Both pass
Clippy, formatting and independent source review. One new recognition fixture
initially failed because constructing the diagonal metric already returned
the dimension. The corrected fixture explicitly checks that normalization and
uses a private plain metric to exercise the lower-level self-loop guard. The
entire production prefix is unchanged, and all 100 original tests are retained.

The frozen dependency family, protocol and observations are under
`/tmp/idenso-factor-recognition-native`. Candidate sources and validation are
under `/tmp/idenso-factor-recognition` and `/tmp/idenso-inline-terminals`.
The earlier sample counts motivated these experiments; they do not establish
the cause of the measured variation or a removable wall-time fraction.

=== Remaining work in the retained implementation

A fresh profile of the retained `d58aee8a` implementation records eight
complete loops per order at 499 Hz. All sixteen results match the complete
9,652-term FORM scalar, retain rank zero, and pass both simplification and
expansion fixed-point checks. Source and runtime identities remain unchanged.
There are 3,316 operation samples for the original order and 1,145 for early
rungs. These repeated-loop native processes differ from the fresh-process
public timing cohort; sampled cycle weights are not wall-time fractions.

#table(
  columns: 3,
  table.header([Disjoint sampled owner], [Original], [Early rungs]),
  [Component collector], [29.01%], [29.63%],
  [Candidate observation], [12.67%], [13.17%],
  [Replacement reconstruction normalization], [12.28%], [9.16%],
  [Interface inference and proofs], [10.86%], [12.73%],
  [Polynomial operations], [9.56%], [6.46%],
  [Symbolica matching/replacement], [6.90%], [8.36%],
  [Dot normalization traversal], [6.77%], [7.34%],
)

The table lists selected disjoint owners, not a complete partition. Nearly all
candidate observation comes from ambient Schoonschip dispatch; collection's
own tensor observation is small. The existing shared index walker can stop
its observer at the first repeated explicit index. Trying the collector at
that point could avoid observing the remaining expression on successful
admission. A declined attempt must complete observation before cleanup uses
callback, bracket or dot flags. An expression with no repeated index already
receives a complete scan, so this design must not add a second scalar scan.

Replacement reconstruction means Rust calls to Symbolica's bulk sum/product
constructors after tensor-leaf replacement, including term merging and
ordering. It does not mean conversion through Python or a new tensor-network
parse. Tensor-safe replacement already retains the proved external interface;
it still has to construct normalized algebra. Skipping that normalization
would need a separate proof for downstream admission and validation.

The intrinsic-normalization scan alone accounts for 3.08% / 2.64% inclusively,
which overlaps the table's owners. Persistent proof caching is therefore a
smaller target than candidate observation, and the publicly mutable tensor
expression forbids blindly retaining such a cache. Collector stack frames
identify attempted intake, not whether an `Option` result succeeds or declines.
No success rate or per-stage term count is inferred from these samples.

The frozen protocol, raw profiles, disjoint attribution and caller analysis
are under `/tmp/idenso-collector-next-analysis/d58-profile`; `final-report.json`
records their identities. This profile supplies optimization targets, not a
measured speedup for an unimplemented change.

=== Early observation: ladder gain and refusal cost

An isolated candidate exposes the existing shared walker's early-stop mode to
Schoonschip dispatch. Eligible expressions try collection at the first repeated
explicit index. Successful collection skips the remaining observations; a
declined attempt restarts full observation before cleanup. No-repeat inputs
finish in one scan. The change also removes the immediate retry of an unchanged
expanded input that the collector already declined; factored-to-expanded
fallback and later changed-normalization retries retain their behavior.

Both Spenso and Idenso are rebuilt against one frozen dependency family. The
complete native loop includes all eight tensor-safe replacements, contractions
and final shared expansion, excluding prepared input construction, checks and
final disposal. The frozen comparison has twelve counterbalanced pairs per
order, with no exclusions or outcome-dependent repeats. Its primary retention
criterion requires lower median paired ratios on wall, process and thread
clocks in both orders, and at least eight of twelve faster wall pairs per order.

#table(
  columns: 5,
  table.header([Complete native loop], [Baseline (s)], [Candidate (s)],
    [Paired wall change], [Faster pairs]),
  [Original], [0.955089], [0.914345], [−3.85%], [10/12],
  [Early rungs], [0.347077], [0.315647], [−8.19%], [12/12],
)

The candidate meets that primary criterion. Paired process/thread changes are
−3.79% for the original order and −8.12% for early rungs. All sixteen intermediate
expressions, complete FORM polynomials, scalar interfaces and fixed points
match. These native times are kept separate from the public Python/FORM cohort.

Separate controls expose the cost of restarting after refusal. Each fixture
uses four counterbalanced pairs of fresh processes, with 32 calls per batch;
these observations are not pooled with the ladder. A wrapper containing 1,024
scalar arguments before its repeated explicit index uses factored intake,
which preserves the subsequent expanded fallback. Its wall batch median rises
from 2.115 to 2.570 ms; the median paired change is *+19.64%*. Process/thread
paired changes are +18.15% / +18.16%, with no faster thread pair. A related late
callback control is mixed and near flat (wall −0.11%, thread +1.77%). The small
unchanged opaque-power refusal also has mixed clocks (wall −17.20%, process
+14.96%, thread −13.60%). All outputs and callback transcripts match, including
the callback-free fixed-point reruns.

This candidate remains isolated despite passing the primary ladder criterion.
The measured refusal regression motivates continuing observation of the
unvisited suffix in the same walk when collection declines. At the end of this
comparison, production and the installed host retained `d58aee8a` while that
follow-up was prepared.

Focused validation passes 234 tests, all-target Spenso/Idenso Clippy, the HEP
library compile check and formatting. The separate HEP library unit-test
target cannot compile in this configuration: five trait-bound errors at
`hep_lib_atom::<AbstractIndex, i32>()` reproduce on the saved baseline. Those
unit tests were not run; the successful library check is compile coverage.
One new admission-count fixture needed explicit registration of its tensor
heads. Its corrected assertion retains the exact one-scan requirement and
additionally checks that collection changes the expression; production is
identical to the measured artifact.

Sources and checks are under `/tmp/idenso-early-observation`. The frozen
48-process primary comparison and separate 24-process controls are under
`/tmp/idenso-early-observation-native`; the source-review records are under
`/tmp/idenso-collector-next-analysis`. Earlier failures and all slower
observations remain in their evidence bundles.

=== Continuing observation after refused collection

The follow-up extends the existing shared index observer with a first-hit
decision callback. It runs once when an explicit index repeats and chooses
whether to stop or observe the remaining nodes. With no repetition it is not
called. Continuing uses the existing traversal with further index bookkeeping
disabled; the observer sees the same nodes in the same order as a full scan.
Constant full-observation and stop decisions serve the existing callers.

Schoonschip attempts the existing collector in that decision. Success stops
observation and returns the collected expression. Refusal completes the
unvisited suffix before callback-sensitive cleanup uses its flags. The
collector still proves source admission before distribution or cancellation.
Observation completion remains distinct from certainty about opaque metadata.
Previously certified intrinsic cleanup can carry pending observations into
the next iteration, preserving the existing fixed-point behavior.

A reviewed edge case needs a final memoized decision when the first repetition
is a direct explicit-slot exponent of a root power. There is no ancestor left
to request observation of a suffix. The fix invokes the decision before any
true result returns, reusing a prior decision when one exists. Regressions
verify that the power survives construction, and check callback count and
ordering for both root and nested cases.

A separate native cohort uses twelve pairs per ladder order against the same
immutable baseline. It includes the final scalar expansion and retains all
observations. Both orders pass the previously specified primary criterion.

#table(
  columns: 5,
  table.header([Complete native loop], [Baseline (s)], [Candidate (s)],
    [Paired wall change], [Faster pairs]),
  [Original], [0.970110], [0.904641], [−6.06%], [12/12],
  [Early rungs], [0.336101], [0.328655], [−3.00%], [10/12],
)

Paired process and thread changes are −6.03% for the original order and −3.02%
for early rungs. All intermediate and final exact checks pass. These
measurements are not pooled with the earlier early-stop experiment or the
public Python/FORM cohort.

Separate controls retain the same 32-call boundary, now with eight pairs of
fresh processes per fixture (sixteen processes). Before measurement, the
late-repeat remedy criterion requires median paired ratios at most 1.05 on
all three clocks and at least six of eight wall and thread pairs at most 1.10.
This tolerance checks the measured regression; it does not establish parity or zero overhead.

The late-repeat control passes narrowly: median paired changes are *+3.29%
wall, +4.90% process and +4.39% thread*. Six wall pairs and seven thread pairs
are within 10%; only three wall and two thread pairs are faster. Its wall
batch medians are 1.737 and 1.794 ms. The largest paired increases, +53.94%
wall / +53.15% process / +41.34% thread, remain in the record without an
asserted cause or exclusion. The other controls favor the candidate in their
paired medians, with small-batch variability retained.

The corrected candidate passes 235 focused tests, all-target Spenso/Idenso
Clippy, the HEP library compile check and formatting. The separate baseline
HEP unit-test compilation limitation remains as described above. Native
qualification permitted the public-host comparison below. Its independently
audited result qualifies the nine reviewed files for retention in production;
the full worktree validation also passes.

The corrected sources and checks are under `/tmp/idenso-resumable-observation`,
with the final manifest in `revision2/source-manifest.json`. Frozen native
protocols, raw observations and `comparison-summary.json` are under
`/tmp/idenso-resumable-observation-native`. Independent source and 96-record
reviews are under `/tmp/idenso-collector-next-analysis`.

=== Public qualification of resumable observation
<resumable-observation-public>

The complete Python host rebuild includes all consumers of the changed Spenso
and Idenso code. The frozen comparison uses the saved `d58aee8a` host and the
new `4eace950` host, with eight counterbalanced pairs per order. Both sides use
`replace_tensor`, the same contraction schedule, and final typed expansion.
The measured boundary excludes prepared input construction, checks and disposal.
Before timing, selection requires lower median paired ratios on wall, process
and thread clocks in both orders, with at least six of eight faster wall pairs.

#table(
  columns: 5,
  table.header([Complete ladder], [Saved host (s)], [New host (s)],
    [Paired wall change], [FORM process (s)]),
  [Original], [0.951250], [0.858030], [−9.09%], [0.759784],
  [Early rungs], [0.320533], [0.300276], [−6.66%], [0.180466],
)

All eight pairs are faster on all three clocks in both orders. Paired process
and thread changes are −9.09% / −6.64% for original / early-rung ordering.
Side medians and median paired changes are distinct statistics. Every sample
and designated warmup is retained; there are no outcome-dependent repeats.
The early-rung order is `[5,4,6,3,7,2,8,1]`, called `reverse` by the driver.

Against the same fresh FORM denominator, the original-order ratio falls from
1.252 to *1.129*, and the early-rung ratio from 1.776 to *1.664*. FORM timings
cover the whole process, including startup, parsing and teardown; Python times
cover the prepared operation loop. These boundaries differ, so neither the
ratios nor the remaining gaps establish equal-scope engine parity. The frozen
native and earlier public cohorts remain separate.

The nine trace controls cover symbolic-dimensional `tracen`, ordinary 4D
`trace4`, repeated indices and an axial trace. First-call wall medians are
below; FORM reports its internal trace-and-sort CPU body and excludes Python
dispatch, wrapping and the axial scalar spectator.

#table(
  columns: 5,
  table.header([Case], [Saved (ms)], [New (ms)], [Paired change], [FORM body (ms)]),
  [Free 2, D], [0.020059], [0.020845], [+4.07%], [0.000640],
  [Free 4, D], [0.072113], [0.070318], [−1.37%], [0.001350],
  [Free 6, D], [0.123821], [0.121820], [−1.62%], [0.005400],
  [Free 10, D], [1.597726], [1.562062], [−2.24%], [0.343333],
  [Free 12, D], [19.593131], [19.338730], [−0.87%], [4.000000],
  [Free 8, 4D], [0.234022], [0.236598], [+1.92%], [0.048667],
  [Repeated 8, 4D], [0.093199], [0.092720], [−0.69%], [0.002000],
  [Repeated 8, D], [0.154918], [0.156764], [+0.40%], [0.014000],
  [Axial 12, 4D], [1.552262], [1.565161], [−0.57%], [0.596000],
)

Each trace route has three paired fresh processes, each retaining five
calibrated inner samples. First calls and unchanged-result reruns are separate.
The axial first-call side median rises while its median paired ratio falls;
neither statistic is substituted for the other. Rerun paired changes are
+4.32%, −28.38%, +1.08%, +1.31%, +0.35%, +0.39%, +5.34%, +1.40% and +1.91%
in the table's order. Six rerun cases are slower in all three pairs. These
controls support no general trace or rerun speedup claim.

The separate final-expansion fixed-point control has twelve processes. Paired
wall changes are +3.03% / +0.15%, with only one of three faster pairs in each
order; CPU changes are about +2.3% / −1%. The native late-repeat control's
remaining slowdown and outlier are retained above. Selection prioritizes the
complete ladder gain without treating these controls as improvements.

Actual-worktree validation passes *964 tests*, with 22 existing skips,
all-target Clippy, formatting and the HEP library compile check. The pre-existing
HEP unit-test compilation failure described above is not counted as a pass.
Public correctness additionally covers 207 API records, 117 HEP component
rows, four network evaluations over two momentum assignments, both complete
FORM polynomials, 83 public controls including tensor-safe replacement,
three notebook routes in both orders
and six metadata cases. The source and dependency inventories remain unchanged
through measurement.

The frozen protocol, all raw timings and validation records are under
`/tmp/idenso-resumable-observation-host`; `formal-summary.json` records the
qualified public result. Native profiles are separate from these wall timings.

The qualified wheel is installed in the notebook environment. A fresh inventory
checks all 28 other installed distributions; their versions and recorded
metadata remain unchanged. Twenty-five selected controls and six metadata cases
pass after installation. The saved baseline and isolated candidate remain
available. Installation evidence is recorded separately in
`canonical-install/summary.json` beneath the same host directory.

=== Remaining costs after resumable observation

A new native profile records eight complete loops per order, with all sixteen
exact FORM, scalar-interface and fixed-point checks passing. It retains 4,583
operation samples for the original order and 1,097 for early rungs at 499 Hz.
These repeated-loop sampled cycle weights are separate from fresh-process
Python wall timings and are not precise fractions of elapsed time.

#table(
  columns: 3,
  table.header([Disjoint sampled owner], [Original], [Early rungs]),
  [Component collection], [32.42%], [32.39%],
  [Other Schoonschip work], [16.22%], [16.93%],
  [Symbolica normalization], [14.28%], [15.85%],
  [Interface inference and proofs], [12.68%], [13.66%],
  [Polynomial operations], [10.14%], [6.41%],
  [Symbolica matching/replacement], [9.48%], [8.20%],
  [Typed replacement], [3.08%], [3.91%],
  [Index occurrence and composition], [1.55%], [2.31%],
)

Collection now executes inside the observer's first-repeat callback. Consequently
the inclusive scanner share, 47.41% / 41.99%, includes nested collection and
cannot be read as scanning cost. Resolved nested collector ancestry accounts
for 38.74% / 33.07%. A disjoint 7.76% / 8.01% remains attributable only to
observation or an inlined decision; it is not assigned entirely to scanning.
Dot normalization traversal accounts for 7.38% / 8.44%. These more detailed
caller views overlap the owner table and must not be added to it.

Reconstructing replaced sums/products accounts for 10.80% / 10.86% in Symbolica
normalization. Within collection, the disjoint method samples include virtual
distribution (7.66% / 7.06%), vector recognition (3.92% / 3.55%), monomial keys
(3.32% / 2.28%), tensor planning (2.50% / 3.53%) and input compilation
(2.20% / 2.68%). These measurements identify repeated work in contraction and
reconstruction as native optimization targets.

The finishing call `with_rewritten_expression` has an inclusive share of only
3.93% / 3.40%, mostly after ambient contraction. That bounds the resolved target
for reusing a collector's interface proof; it is not the whole inference share.
The samples do not distinguish initial successful collection from other
finishing calls, so they neither prove that all of this work is redundant nor
predict a speedup. The native comparison below tests a certificate from initial
collection while retaining validation after callbacks or earlier cleanup.

There are 63 / 5 unresolved leaf frames; inlining and stack resolution limit
finer attribution. The interrupted first launch stopped at its environment
check before any samples. A versioned profile verified the same six shared
library hashes under the current environment before and after recording;
all source and executable guards passed. Its frozen protocol, raw stacks and
`final-report.json` are under
`/tmp/idenso-collector-next-analysis/resumable-profile/v2`. Old and new profile
percentages are not subtracted to infer runtime gains.

=== Rejected refinement of dot candidates

An isolated follow-up made the observation pass distinguish ordinary indexed
vectors, such as `p(mink(4,a))`, from nested vectors and explicit vector powers.
Ordinary indexed vectors cannot trigger a nested-vector identity, so the change
can avoid an otherwise unchanged dot-normalization walk. Powers, opaque
payloads and callback-sensitive admission retain their checks. This also makes
some inert vector metadata eligible for existing component collection.

The candidate passes 215 focused Idenso/Spynso tests, all-target Clippy and
formatting. Compiler paths, dependency records and the four new test names bind
these checks to the candidate. Native checks preserve all sixteen literal
intermediate ladder expressions, both exact 9,652-term FORM results, scalar
interfaces and fixed points, plus five scalar-expansion cases.

Nevertheless, it fails the fixed native performance criterion. Twelve fresh
pairs per order compare the retained `4eace950` source with this one-file
candidate. The measured operation includes the complete reduction and final
typed expansion; setup, checks and disposal are outside the clocks. These
native measurements are separate from the public Python timings above.

#table(
  columns: 5,
  table.header([Order], [Retained median (s)], [Candidate median (s)],
    [Paired wall change], [Faster pairs]),
  [Original], [0.941312], [0.935361], [+0.22%], [6/12],
  [Early rungs], [0.328228], [0.337743], [+6.49%], [3/12],
)

Paired process/thread changes are +0.25% / +0.25% for the original order and
+6.39% / +6.39% for early rungs. Neither order satisfies the prespecified
requirement of lower median paired ratios on all three clocks and at least
eight of twelve faster wall pairs. The slightly lower original-order side
median is not substituted for its paired result.

Separate controls retain eight pairs per fixture, each containing 32 calls on
an unchanged input. Paired wall changes are −15.93% for one explicit vector
(six faster pairs) and −14.90% for a 256-vector sum sharing one free port
(eight faster pairs). The compact-dot control is mixed: −1.97% wall, +5.95%
process and −1.21% thread. Short batches retain substantial outliers, including
a 2.85-times wall ratio for one single-vector pair. These controls model the
native wrapper work but omit Python object and descriptor allocation.

The local no-op improvement does not justify the ladder regression. The
candidate is not integrated or installed; the retained public build and FORM
comparisons above are unchanged. All 48 ladder and 48 control observations,
including slower pairs, are retained without outcome-dependent repeats under
`/tmp/idenso-dot-candidate-flags-native`. Source and validation records are
under `/tmp/idenso-dot-candidate-flags`. The measurements do not by themselves
attribute the regression to a particular added classification step.

=== Further native comparisons

Three independent candidates were compared with the same retained `4eace950`
implementation. Each cohort contains twelve fresh pairs per ladder order and
includes the complete eight-step reduction plus final typed expansion. Setup,
correctness checks and disposal are outside the operation clocks. Every
candidate preserves all sixteen literal intermediate expressions and both
exact 9,652-term FORM results, scalar interfaces and fixed points. The existing
component certificate is reused through exact final equality; these native
experiments do not claim new HEP evaluations or new FORM timings.

The acceptance rule was fixed before measurement: both orders need lower
median paired ratios on wall, process and thread clocks, and at least eight of
twelve faster wall pairs. The following percentages are paired wall changes,
each from its own cohort; absolute times across cohorts are not compared.

#table(
  columns: 5,
  table.header([Change], [Original], [Early rungs],
    [Faster pairs, original / early], [Decision]),
  [Reuse verified interface], [−4.94%], [−2.87%], [10/12 / 10/12],
    [Rejected in Python],
  [Check vector head first], [+10.73%], [−1.93%], [4/12 / 7/12], [Rejected],
  [Stop dot preflight globally], [−2.65%], [+4.26%], [7/12 / 4/12], [Rejected],
)

*Interface reuse.* Initial successful component admission can certify that the
contraction preserves an established explicit interface. A shared
`SymbolicTensor` operation consumes this fact instead of inferring the result's
interface again. Admission after earlier cleanup cannot certify the original
input. Callback-bearing or unresolved cases retain checked reconstruction,
including the counterexample where contraction changes `T(b)` into a scalar.
Identity results reuse the computed Atom, and typed zeros retain their ports.
The candidate passes 215 focused Idenso/Spynso tests, all-target Clippy and
formatting in a private, source-bound build.

Native side wall medians are 0.928064 → 0.890104 seconds for the original order
and 0.324804 → 0.314515 seconds for early rungs. Paired process/thread changes
are −4.90% / −4.90% and −2.82% / −2.82%, with ten faster pairs on each clock.
The native main loop omits Python's initial structured-input clone. Separate
controls model that clone and identity wrapping, but still omit Python object
and descriptor allocation. Eight pairs of 32 unchanged-input calls show a
*+35.20%* paired wall change for a tiny scalar and *+3.30%* for the 9,652-term
scalar; all eight wall pairs are slower in both controls. The tiny scalar's
CPU records include a zero thread duration and large outliers, so they do not
support precise CPU percentages. All records remain included. The subsequent
public Python comparison below rejects this candidate.

*Vector admission.* Checking arity and the rank-one tag before asking the slot
matcher to recognize a vector preserves malformed-input and conflicting-tag
checks. It passes 160 focused tests in separate processes and Clippy for the
Idenso test target. Its original-order paired process/thread changes are both
+10.53%; early-rung changes are both −2.11%. Neither order meets the full rule.
The candidate remains unintegrated; the measurements do not establish the cause
of the regression.

*Dot preflight.* Symbolica's visitor can prune a subtree, but its active ancestor
loops still invoke the callback for later siblings. A private recursive search
can stop globally at the first rewrite while preserving traversal order, opaque
boundaries and the first rewrite's callback effects. The candidate passes 161
focused tests in separate processes and Clippy for the Idenso test target.
Its original-order paired process/thread changes are both −2.64%; early-rung
changes are both +4.24%. It also fails the complete-ladder rule.

Separate controls contain eight pairs of 32 calls for an early rewrite, a late
rewrite, no rewrite, and an opaque payload. Paired wall changes are −6.54%,
−7.64%, −12.66% and −9.11%, with six, six, eight and seven faster pairs. These
measure the whole raw Schoonschip operation, including observation and any
subsequent replacement and cleanup; they are not preflight-only timings. The
same original input is used for each call, with all outputs retained until the
clocks stop. These improvements do not justify the ladder regression.

The standalone `visitor_shortcircuit` binary in
`examples/reproducers/symbolica-expansion` demonstrates the visitor behavior
without Idenso. For a first-argument match followed by 8,192 arguments, callback
counts are 8,194 versus two. All twelve first/last/absent/opaque cases, formatting
and Clippy checks pass. These are visit counts, not a runtime or ladder-speed
claim; the accompanying `performance.typ` describes how to run it.

All 80 interface-reuse, 48 vector-admission and 112 preflight observations are
retained without exclusions or outcome-dependent repeats. Our builds and other
benchmarks were quiet during each cohort; the shared host was not exclusive.
Frozen protocols, raw records and independent arithmetic/guard reviews are in
`/tmp/idenso-typed-schoonschip-proof-native`,
`/tmp/idenso-vector-admission-order-native`,
`/tmp/idenso-dot-preflight-shortcircuit-native` and
`/tmp/idenso-collector-next-analysis`. None of these candidates is installed;
the retained public timings and FORM ratios above are unchanged.

=== Public comparison of interface reuse

The isolated `23a302a3` Python host passes 207 API records, 117 HEP component
rows, four complete network evaluations, both fresh FORM polynomials, 83 public
controls, all three notebook routes in both orders and six metadata cases.
Additional checks cover the actual Python no-op calls and the raw Symbolica
reference. Callback-induced rank loss remains rejected. The candidate therefore
passes correctness, but fails the performance requirement fixed before timing.

Eight fresh pairs per order compare complete reductions against the retained
`4eace950` host. Both orders require lower median paired ratios on all three
clocks and at least six of eight faster wall pairs.

#table(
  columns: 5,
  table.header([Order], [Retained median (s)], [Candidate median (s)],
    [Paired wall change], [Faster pairs]),
  [Original], [1.343208], [1.312954], [−0.61%], [4/8],
  [Early rungs], [0.399306], [0.426222], [+4.10%], [3/8],
)

Paired process/thread changes are −0.50% for the original order and +4.09% for
early rungs. Separate public no-op controls contain eight pairs per fixture:
4,096 calls on a tiny scalar or 32 calls on the 9,652-term scalar, with outputs
retained through the clocks. Paired wall changes are *+27.79%* and *+6.12%*,
with only one and two faster pairs respectively. Unlike the native main loop,
this public path clones its structured input before contraction; this source
difference is established, but its contribution to the measured changes has
not been isolated.

The supplied raw Symbolica recipe is also measured separately on the retained
host, with three fresh processes per order: medians are 1.081411 and 0.384681 s.
Fresh FORM process medians are 1.035341 and 0.249036 s. The raw recipe uses
specialized scalar rules, while FORM includes process startup and teardown;
neither comparison isolates the cost of tensor-interface checks. These cohorts
are not pooled with earlier timings or the primary paired comparison.

Trace controls retain all three process pairs and five inner batches per route.
Free length-12 tracen has side medians 32.468 → 23.502 ms, but only a −2.71%
median paired change; these statistics must not be substituted for each other.
Axial length-12 unchanged-result reruns have a +19.58% paired change. The trace
algorithm is unchanged, and these variable controls do not establish a general
trace improvement. Separate final-expansion fixed-point controls have paired
wall changes of +157.42% and +18.74% for original and early-rung inputs. Their
expansion implementation is unchanged; the measurements do not establish the
cause of these differences.

All eight measurement phases and source/dependency guards pass. Every sample
is retained without outcome-dependent repeats; our other workloads were quiet,
but external work continued on the shared host. Full clocks, trace controls,
FORM output and correctness records are under
`/tmp/idenso-typed-schoonschip-proof-host/benchmarks`. The candidate is neither
integrated nor installed. The retained tensor-safe replacement implementation
and its qualified results above remain current.

=== Passing owned operands to replacement builders

An independent two-line candidate passes the already-collected
`Vec<AtomOrView>` directly to Symbolica's bulk sum and product builders, rather
than creating iterators of borrowed views. Matching, callback order, unchanged
branch reuse, normalization and post-product validation remain identical.
The possible benefit is allocation reuse; it does not remove normalization.
All fifteen existing replacement regressions pass in separate processes, along
with Clippy for that Idenso test target and formatting. Native checks also
preserve all sixteen intermediate ladder expressions, both exact FORM results
and five scalar controls.

The same twelve-pair native criterion rejects this candidate. Original-order
paired wall/process/thread changes are +0.89% / +0.84% / +0.84%, with five of
twelve faster pairs. Early-rung changes are −6.39% / −6.30% / −6.30%, with
eight faster pairs. Both orders had to pass. Wall side medians are
1.037943 → 1.238934 s and 0.392143 → 0.365881 s; the substantial difference
between paired and side-median statistics remains in the record.

All 48 observations and source/runtime guards are retained without exclusions
or repeats under `/tmp/idenso-owned-replacement-builders-native`. These native
clocks include final expansion but omit Python wrapping. No fresh FORM timing
is claimed, and no cause is inferred from the variable results. The candidate
is neither integrated nor installed.

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
