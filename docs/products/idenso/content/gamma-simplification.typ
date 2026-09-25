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
cargo run --locked -p idenso --profile dev-optim --example metric_contraction_benchmark -- /tmp/contraction-full 5 8
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
`StructuredAtom` keeps its existing Atom storage until bracket normalization
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

Port rewriting borrows unchanged branches, including identity substitutions.
Every sum branch still independently accounts for the requested ports; an
index missing from one summand is an error. Genuine relabeling uses Symbolica's
bulk sum and product builders when arithmetic is exact and callbacks cannot be
reordered. Rounded coefficients and callback-sensitive expressions retain their
existing evaluation order.

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
