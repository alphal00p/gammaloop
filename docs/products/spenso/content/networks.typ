#import "../../shared.typ": callout, boundary, product-link, source-link

#let networks = [
= Tensor data, contraction, and network execution

Spenso keeps tensor structure and tensor data separate so the same structural operations work
with dense, sparse, symbolic, or parametric storage. A structure supplies ordered slots,
representations, dimensions, and abstract indices. Data supplies values and a storage-specific
access strategy. Contraction first establishes structural compatibility and only then chooses a
data operation.

== Dense, sparse, and parametric data

Dense storage is appropriate when most components are populated and a predictable contiguous
layout matters. Sparse storage is appropriate when zero structure dominates, but contraction
cost depends on index overlap and lookup behavior rather than only the number of stored values.
Parametric and symbolic tensors defer numeric substitution; keep their parameter environment
with the tensor when serializing or moving a calculation between stages.

Conversions must preserve the structure's slot order. A tensor with the same dimensions but a
different representation or variance is not interchangeable until an explicit structural
operation establishes that relationship.

== Contraction and planning

Pairwise contraction matches compatible dual slots, determines the output structure, and then
contracts the selected data backends. A tensor network generalizes this to many nodes connected
through abstract indices. Planning selects a contraction order; execution applies that plan to
the current data.

The cheapest-looking local pair is not always the cheapest complete plan. Intermediate tensor
dimensions, sparsity, symbolic expression growth, and reusable subexpressions can dominate.
Record the planner/configuration with performance results and invalidate a stored plan when the
network structure changes.

#callout("Plan structure, execute data", [
  A plan is meaningful for the network structure it was derived from. Replacing values without
  changing slots can reuse it; adding a tensor, changing a representation, or rewiring an index
  requires validation or replanning.
])

== Symbolica shadowing

With `shadowing`, Spenso recognizes registered tensor functions inside Symbolica expressions
and builds a network whose nodes shadow those subexpressions. Initialization of representations
and the tensor library must happen before parsing. Unknown functions remain symbolic rather
than silently acquiring tensor semantics.

Shadowing is useful for extracting a contraction problem while retaining a reversible link to
the source expression. Keep the replacement map until the computed result has been substituted
back. For rewrite identities and index cooking, use #product-link("idenso", label: "Idenso");
Spenso itself owns storage, contraction, and execution.

== Custom representations and the HEP library

`SimpleRepresentation` derives the repetitive representation plumbing for a user enum. The
derive input still defines domain semantics: duality, dimension, slot compatibility, and any
Symbolica tags must agree. Generated representations are ordinary Rust API and appear in the
exhaustive `spenso-macros` sidecar.

`spenso-hep-lib` registers concrete high-energy-physics structures such as gamma matrices and
projectors. Treat that library as shared domain data. GammaLoop links to it rather than copying
its definitions into the application manual.

#boundary("Graph and tensor ownership", [
  #product-link("linnet", label: "Linnet") owns connectivity, cuts, cycles, and graph mutation.
  Spenso owns which connected slots contract and how a plan is executed. Keeping that boundary
  explicit prevents a graph rewrite from being mistaken for a tensor-algebra identity.
])

Execution regression notes are maintained in
#source-link("docs/architecture/spenso-network-execution-test-baseline.md", label: "the network execution baseline").
]
