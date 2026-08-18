#import "../../shared.typ": callout, developer-link, product-link

#let overview = [
= Overview

Idenso is the symbolic identity and simplification layer for tensors encoded in the form that
Spenso can parse from Symbolica expressions. It provides representation symbols, index
tooling, reversible “cooking” of subexpressions, selective expansion, and identities for Dirac,
metric, epsilon, and color algebra.

#callout("Symbolic, not component expansion", [
  Idenso rewrites abstract tensor expressions. Functions such as `expand_mink` distribute
  factorized terms that carry selected index families; they do not turn every tensor into a
  dense array of components. Concrete tensor storage and network execution belong to
  #product-link("spenso", label: "Spenso").
])

== Choose a task

- To initialize representations and verify one metric contraction, follow the
  #link("tutorial/")[controlled identity tutorial] and the
  #link("reference/python/idenso-community/#simplify-metrics")[Python function reference].
- To isolate dummy-index namespaces or cook a large expression, use the
  #link("guides/algebra/")[algebra guide] with the exact
  #link("reference/rust/supported/idenso/#supported-indextooling")[`IndexTooling`] and
  #link("reference/rust/supported/idenso/#supported-cookable")[`Cookable`] references.
- To audit a sign or normalization in a Dirac/color pass, use the source-backed
  #link("reference/form-color-dirac/")[FORM, Symbolica, and Spenso rule specification].

== A controlled rewrite pipeline

There is no universally correct “simplify everything” order. A robust workflow makes each
phase explicit:

- initialize the standard representations and tensor symbols;
- inspect dangling indices and normalize or wrap dummy-index namespaces when combining
  expressions;
- expand only the sectors needed by the next identity pass;
- simplify gamma, metric, or color structures as appropriate;
- canonicalize and compare results using the conventions of the consuming calculation.

Expansion can grow an expression quickly, so perform it as late and selectively as possible.
Wrapping indices is especially important before multiplying independently constructed
expressions: equal printed index names can otherwise acquire an unintended contraction.

== Representations and cooking

The representation layer defines spin-fundamental, color-fundamental, color-sextet, bispinor,
and color-adjoint types together with their duality. `initialize()` forces the related
Symbolica symbols and Spenso tags to be registered before parsing or rewriting expressions.

The `Cookable` API can replace selected function-like subexpressions or index payloads with
compact symbols and later reverse the operation. Use reversible settings when the original
expression must be recovered. A cooked symbol is an intermediate encoding, not a portable
serialized physics result unless the associated mapping and settings are retained.

Idenso is consumed by #product-link("gammaloop", label: "GammaLoop") and its tensor conventions
are shared with #product-link("spenso", label: "Spenso"). Continue to the
#link("reference/interfaces/")[interface guide] and generated Rust/Python references for
supported functions, signatures, and feature gates.

Contributors changing representation syntax, symbolic-network parsing, index ownership, or
rewrite order should also read the source-audited
#developer-link("idenso-architecture", "idenso-architecture.typ", "Idenso implementation architecture").
]
