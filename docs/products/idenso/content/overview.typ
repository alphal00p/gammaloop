#import "../../shared.typ": callout, boundary, product-link, source-link

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

#boundary("Specifications versus supported entry points", [
  The files under `docs/idenso/` describe detailed FORM/Symbolica conventions and rewrite
  investigations. Their representation names and supported rewrite boundaries were checked
  against the current registrations and tests while preparing this manual. Experimental
  performance observations remain developer notes and do not silently become compatibility
  promises.
])

Idenso is consumed by #product-link("gammaloop", label: "GammaLoop") and its tensor conventions
are shared with #product-link("spenso", label: "Spenso"). Start with the
#source-link("crates/idenso/README.md", label: "Idenso README") and generated API catalog.
]
