#import "../../shared.typ": callout, boundary, product-link, source-link

#let overview = [
= Overview

Vakint matches single-scale vacuum integrals to a library of known topologies, canonicalizes
their momentum routing, reduces tensor numerators, and evaluates them analytically or
numerically. Expressions are represented with Symbolica and contain a numerator multiplied by
a `topo(...)` structure built from propagators.

#callout("Separate matching from evaluation", [
  Canonical topology matching is useful on its own and is much cheaper than a complete
  evaluation. Tensor reduction and analytic backends invoke FORM; numerical sector decomposition
  invokes pySecDec and FORM. Choose the narrowest operation and evaluation order needed for the
  task.
])

== Workflow

A complete calculation proceeds through explicit stages:

- create one `Vakint` engine and reuse it, because topology setup has a cost;
- parse or construct a Symbolica expression in Vakint's namespace;
- canonicalize the topology and momentum routing;
- tensor-reduce the numerator when required;
- evaluate the integral with the first supported backend in the configured order;
- substitute numerical parameters and external momenta when a numerical result is required.

`VakintExpression` is the structured sum of numerator/topology terms. A recognized topology can
be rendered in canonical long or short form. If unknown integrals are allowed, matching may
preserve them, but evaluation still needs a backend that supports the resulting topology.

== External tools and reproducibility

#boundary("Construction validates the selected backends", [
  The default evaluation order includes AlphaLoop, MATAD, FMFT, and pySecDec. Settings validation
  therefore probes FORM and pySecDec even if a later example only canonicalizes an expression.
  FORM 4.2.1 or newer and pySecDec 1.6.4 or newer are the minimum versions checked by the current
  source. For a pure matching test, use an empty evaluation order in Rust; in Python, provide an
  empty `evaluation_order` when constructing `Vakint`.
])

Analytic AlphaLoop, MATAD, and FMFT methods depend on FORM. The pySecDec method depends on both
FORM and pySecDec and requires numerical masses and, where applicable, external momenta. Record
tool versions, normalization convention, epsilon depth, decimal precision, backend order, and
temporary-output reuse policy with any result intended to be reproduced.

== Test tiers and platform notes

The ordinary fast suite excludes pySecDec work while retaining verbose diagnostics and generated
temporary files for inspection:

```sh
VAKINT_SKIP_PYSECDEC_TESTS=true RUST_BACKTRACE=full VAKINT_NO_CLEAN_TMP_DIR=T \
  RUST_LOG=DEBUG cargo test --package vakint --no-default-features
```

Pure parsing, matching, canonicalization, and settings tests belong in the normal CI tier.
FORM-backed examples require a provisioned FORM job. pySecDec, MATAD, and FMFT comparisons are
opt-in or scheduled because they depend on external installations and can be slow. On macOS,
a GNU-compiler link failure mentioning missing `__emul...` symbols is handled by building with
`EXTRA_MACOS_LIBS_FOR_GNU_GCC=T`, as consumed by Vakint's build script.

Vakint is used by #product-link("gammaloop", label: "GammaLoop") for vacuum-integral work, but it
owns its topology and evaluation conventions independently. Start with the
#source-link("crates/vakint/README.md", label: "Vakint README") and generated catalogs.
]
