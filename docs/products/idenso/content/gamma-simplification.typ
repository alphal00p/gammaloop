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
