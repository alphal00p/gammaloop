#import "../../shared.typ": callout, product-link, source-link

#let algebra = [
= Syntax, indices, and algebra passes

Idenso consumes the Symbolica function form emitted for Spenso structures. Function heads carry
representation meaning and their arguments carry abstract indices or tensor data. Initialize
Idenso before parsing so Symbolica attributes, Spenso tags, and dual representations are all
registered in one deterministic order.

== Indices and cooking

Free indices describe the result; repeated compatible indices describe contractions. Before
combining independently built expressions, wrap or rename dummy-index namespaces so equal
printed names do not create accidental contractions. Cooking temporarily replaces selected
index payloads or function subexpressions with compact symbols. Retain the generated map and
settings whenever uncooking is required.

Cooking is especially useful before expensive canonicalization, but it changes what downstream
matchers can see. Uncook before applying an identity whose pattern depends on the hidden tensor
head or index structure.

== Metric and epsilon operations

Metric contraction raises, lowers, or identifies compatible Lorentz indices according to the
registered representation. Epsilon identities depend on dimension, ordering, and sign
conventions; expand them only when the next pass benefits from the larger expression. Keep
Minkowski expansion selective to avoid distributing unrelated scalar factors.

== Dirac and color algebra

Dirac passes simplify gamma chains, traces, slashes, and spinor-compatible contractions. Color
passes handle fundamental/adjoint deltas, generators, structure constants, and registered
group parameters. Apply one algebra family at a time and inspect the intermediate expression;
an all-at-once fixed-point loop can obscure which convention produced a sign or normalization.

#callout("Canonical does not mean physically equivalent by itself", [
  Canonical ordering makes structurally equivalent expressions comparable under the registered
  rules. On-shell relations, gauge choices, dimension-specific identities, and model parameter
  substitutions must still be requested explicitly.
])

== Schoonschip network parsing

The Schoonschip path parses large gamma-chain expressions into the same Spenso-compatible
network model before contraction. It resolves aliases and normalizes the expression before
constructing the network. Spenso then plans and executes contractions, while Idenso supplies
the representation and algebra rules. Avoid substituting scalar parameters too early: doing so
can substantially increase intermediate expression size even when the resulting tensor network
is unchanged.

Concrete syntax and rewrite cases are documented in
#source-link("docs/architecture/spenso-symbolica-syntax-and-rewrites.typ", label: "the Spenso/Symbolica syntax note")
and the rendered #link("manual/form-color-dirac/")[FORM color and Dirac specification].
The #source-link("docs/architecture/schoonschip-net-parsing.typ", label: "Schoonschip parsing guide")
shows how normalization, network construction, and contraction fit together.

Tensor storage and execution remain owned by #product-link("spenso", label: "Spenso").
]
