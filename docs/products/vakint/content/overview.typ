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

== Choose a task

- To recognize and canonicalize an integral without external tools, follow the
  #link("tutorial/")[matching-only tutorial] and compare with the
  #link("reference/topologies/")[generated topology table].
- To tensor-reduce and evaluate with an explicit backend policy, use the
  #link("manual/evaluation/")[evaluation manual] with the exact
  #link("reference/rust/supported/vakint/#supported-evaluationorder")[`EvaluationOrder`]
  reference.
- To embed Vakint in Rust or `symbolica.community.vakint`, use the
  #link("manual/interfaces/")[interface guide] and generated Rust/Python signatures before
  selecting external tools.

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
  Vakint requires FORM 4.2.1 or newer and pySecDec 1.6.4 or newer. For matching without
  evaluation, use an empty evaluation order in Rust; in Python, provide an empty
  `evaluation_order` when constructing `Vakint`.
])

Analytic AlphaLoop, MATAD, and FMFT methods depend on FORM. The pySecDec method depends on both
FORM and pySecDec and requires numerical masses and, where applicable, external momenta. Record
tool versions, normalization convention, epsilon depth, decimal precision, backend order, and
temporary-output reuse policy with any result intended to be reproduced.

== Diagnostics and platform notes

Matching with an empty evaluation order needs neither FORM nor pySecDec. For backend evaluation,
check the configured executable paths first. Verbose logging and retained temporary files are
useful when an external tool fails:

// docs-example: syntax
```sh
VAKINT_NO_CLEAN_TMP_DIR=T RUST_LOG=DEBUG cargo run
```

On macOS, if a GNU-compiler link fails with missing `__emul...` symbols, retry the build with
`EXTRA_MACOS_LIBS_FOR_GNU_GCC=T` set.

Vakint is used by #product-link("gammaloop", label: "GammaLoop") for vacuum-integral work, but it
owns its topology and evaluation conventions independently. The
#source-link("crates/vakint/README.md", label: "Vakint README") contains additional Rust examples;
use the API and topology references in this manual for the selected version.
]
