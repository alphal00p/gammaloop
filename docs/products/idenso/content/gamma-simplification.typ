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

Idenso's ordinary closed traces currently use signed pairing recursion. FORM's `trace4`
adds trace-specific reductions, including a four-dimensional reduction for distinct
arguments. In a FORM 5.0.0 probe, ten free indices gave 693 terms with `trace4` and 945
with `tracen`; Idenso's generic recursion gives 945. These are different representations
of the same four-dimensional tensor, not evidence of an algebraic discrepancy. See the
#link("https://github.com/form-dev/form/blob/master/doc/manual/gamma.tex")[FORM Dirac-algebra manual].

When FORM is installed locally, the notebook runs the selected compact scalar case with
both `trace4` and `tracen`, plus a fourteen-index trace that better exposes their cost
and output-size differences. It checks the compact results and the free-index results
contracted with paired momenta against exact scalar identities. Generated programs and
FORM output remain visible in the notebook. Native process timings include startup,
parsing, sorting and verification, so they are not comparable as isolated kernels to
Idenso's in-process timing. This establishes selected identities and performance examples,
not full FORM correctness or performance parity. Browser exports skip native execution.

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
