#import "../../shared.typ": callout, boundary, product-link

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

== Execute a small tensor/scalar network

This self-contained network computes `2 a + 3 b` for two equally structured tensors. The dummy
libraries make the execution boundary explicit: every tensor is already concrete and any
unresolved function key is an error.

// docs-example: compile
```rust
use spenso::{
    network::{
        ExecutionResult, Network, Sequential, SmallestDegree,
        library::{DummyKey, DummyLibrary, DummyLibraryTensor, panicing::ErroringLibrary},
        store::NetworkStore,
    },
    structure::{
        OrderedStructure,
        representation::{Euclidean, RepName},
    },
    tensors::data::DenseTensor,
};

fn main() {
    type Tensor = DenseTensor<f64, OrderedStructure<Euclidean>>;
    type Store = NetworkStore<Tensor, f64>;
    type Net = Network<Store, DummyKey, DummyKey>;
    type LibraryTensor = DummyLibraryTensor<Tensor>;

    let structure = OrderedStructure::new(vec![Euclidean {}.new_slot(2, 1)]).structure;
    let a = DenseTensor::from_data(vec![1.0, 2.0], structure.clone()).unwrap();
    let b = DenseTensor::from_data(vec![3.0, 4.0], structure).unwrap();
    let tensors = DummyLibrary::new();
    let functions = ErroringLibrary::new();

    let mut network =
        Net::from_scalar(2.0) * Net::from_tensor(a) + Net::from_scalar(3.0) * Net::from_tensor(b);
    network
        .execute::<Sequential, SmallestDegree, LibraryTensor, _, _>(&tensors, &functions)
        .unwrap();

    let ExecutionResult::Val(result) = network.result_tensor::<LibraryTensor, _>(&tensors).unwrap()
    else {
        panic!("expected a concrete tensor result");
    };
    assert_eq!(result.data, vec![11.0, 16.0]);
}
```

The output invariant is the component vector `[11.0, 16.0]`; the original rank-one structure
must also remain attached to the result. The documentation harness compiles the complete program;
run it to execute the assertion in a provisioned Symbolica environment. See the exact
#link("reference/rust/supported/spenso/#supported-network")[`Network`],
#link("reference/rust/supported/spenso/#supported-densetensor")[`DenseTensor`], and
#link("reference/rust/supported/spenso/#supported-contract")[pairwise contraction] references.

#callout("Interpret execution failures by layer", [
  `from_data` errors are rank/dimension or storage-length mismatches. An unresolved function-key
  error means the network contains a symbolic library node despite the concrete-only setup. A
  different value with the same shape points to scalar placement or node arithmetic; a different
  shape points to changed slots or unintended contraction, which requires replanning.
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
Symbolica tags must agree. The resulting representations implement the ordinary Rust traits
used by the rest of Spenso; the `spenso-macros` API reference lists the supported derive
attributes.

`spenso-hep-lib` registers concrete high-energy-physics structures such as gamma matrices and
projectors. Treat that library as shared domain data. GammaLoop uses the same definitions.

#boundary("Graph and tensor ownership", [
  #product-link("linnet", label: "Linnet") owns connectivity, cuts, cycles, and graph mutation.
  Spenso owns which connected slots contract and how a plan is executed. Keeping that boundary
  explicit prevents a graph rewrite from being mistaken for a tensor-algebra identity.
])

]
