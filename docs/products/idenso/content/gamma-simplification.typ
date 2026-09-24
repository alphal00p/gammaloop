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
normalization, 5.330 ms for Schoonschip and 6.722 ms for epsilon cleanup. Chain
collection is now the largest measured stage: `chainify` takes 20.957 ms and
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
