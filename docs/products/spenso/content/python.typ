#import "../../shared.typ": boundary, callout, source-link

#let python = [
= Python tensor workflows

Spenso's Python interface is a native adapter over the same tensor structures, libraries, and
network executor as the Rust crates. Use the generated reference for exact signatures and this
guide for the object boundaries and execution sequence that those signatures do not explain.

== Availability and version boundary

#boundary("A Symbolica community module", [
  Import Spenso as `symbolica.community.spenso`; there is no standalone `spenso` wheel. The
  installed Symbolica assembly determines whether the module is present and which Spenso API it
  contains. A source checkout or generated `.pyi` file does not add the native module to an
  existing Python environment.
])

Check the environment before building a workflow:

// docs-example: syntax
```sh
python -c "import symbolica.community.spenso as spenso; print(spenso.__name__)"
```

When a provider's Symbolica package does not include Spenso, build it through the
#link("https://github.com/benruijl/symbolica-community")[community-module assembly]. Record the
Symbolica assembly version with reproducible results; it is a more useful Python compatibility
fact than the version of an unrelated local Rust checkout.

== Choose the right object

- `Representation` and `Slot` define dimensions, duality, and abstract indices.
- `TensorIndices` carries concrete abstract indices; `TensorStructure` describes a reusable
  shape whose indices can be assigned later.
- `Tensor` owns ordinary dense or sparse data. `LibraryTensor` is the named form registered in a
  `TensorLibrary` for reuse by symbolic networks.
- `TensorNetwork` owns an expression graph. `ExecutionMode` selects the rewrite strategy:
  one smallest-degree rewrite per step, scalar work only, or the general smallest-degree strategy.
  Use `n_steps`, not the mode, to bound how many execution steps run.
- `TensorEvaluator` substitutes repeated numerical parameter batches into a symbolic tensor;
  its compiled counterpart owns the generated native evaluator.

The generated #link("reference/python/spynso3/Representation/")[`Representation`],
#link("reference/python/spynso3/Tensor/")[`Tensor`], and
#link("reference/python/spynso3/TensorNetwork/")[`TensorNetwork`] entries are the exact
versioned contracts. Do not infer index compatibility from two equal dimensions: names,
representations, and duality remain part of the value.

== Construct and inspect concrete data

This complete source creates a named rank-two tensor, verifies one component, and converts its
storage in place. The documentation harness compiles the Python source without importing the
native module.

// docs-example: compile
```python
from symbolica.community.spenso import Representation, Tensor, TensorIndices

rep = Representation.euc(2)
i = rep("i")
j = rep("j")
indices = TensorIndices(i, j)
matrix = Tensor.dense(indices, [1.0, 0.0, 0.0, 1.0])

assert len(matrix) == 4
assert matrix[0, 0] == 1.0
matrix.to_sparse()
assert matrix[1, 1] == 1.0
```

`Tensor.dense` requires row-major data whose length is the product of the structure dimensions.
`Tensor.sparse` instead needs the element type and starts empty. `to_dense()` and `to_sparse()`
mutate the storage representation; they do not change slots or re-index the tensor. See the
#link("reference/python/spynso3/Tensor/#exports-tensor-dense-associatedfunction")[dense constructor] and
#link("reference/python/spynso3/Tensor/#exports-tensor-to-sparse-method")[conversion contract].

#callout("Diagnose structure before storage", [
  A constructor failure usually means the data length and dimensions disagree. An unexpected
  contraction or exterior product is instead an index/duality problem. Print the structure and
  slots before changing dense/sparse storage, because a storage conversion cannot repair a
  structural mismatch.
])

== Register data and execute a network

A symbolic network resolves named tensors through a library. Register the data first, construct
the expression with the same name and compatible slots, then execute and request the result:

// docs-example: compile
```python
from symbolica.community.spenso import (
    ExecutionMode,
    LibraryTensor,
    Representation,
    TensorLibrary,
    TensorName,
    TensorNetwork,
    TensorStructure,
)

rep = Representation.euc(2)
A = TensorName("A")
structure = TensorStructure(rep, rep, name=A)
library = TensorLibrary()
library.register(LibraryTensor.dense(structure, [1.0, 0.0, 0.0, 1.0]))

network = TensorNetwork(A(rep("i"), rep("j")), library=library)
network.execute(library=library, mode=ExecutionMode.All)
result = network.result_tensor(library=library)
assert len(result) == 4
```

`ExecutionMode.Single` selects the smallest-degree single-rewrite strategy, but execution still
continues while work remains unless it is bounded; use `n_steps=1` to inspect exactly one step.
`Scalar` processes scalar work while retaining tensor structure, and `All` attempts the complete
available execution. A successful `execute()` does not make `result_scalar()` appropriate for a tensor result; choose
`result_tensor()` or `result_scalar()` according to the remaining structure. Exact behavior is
linked from #link("reference/python/spynso3/TensorNetwork/#exports-tensornetwork-execute-method")[`execute`],
#link("reference/python/spynso3/TensorNetwork/#exports-tensornetwork-result-tensor-method")[`result_tensor`],
and #link("reference/python/spynso3/ExecutionMode/")[`ExecutionMode`].

#callout("Interpret a network failure by ownership", [
  A missing-library error means a symbolic tensor name has no registered data. A result-kind
  error means the graph still has tensor structure when a scalar was requested, or vice versa.
  A changed network structure invalidates assumptions about contraction order; replacing only
  values with the same registered structure does not.
])

== Repeated symbolic evaluation

Use `Tensor.evaluator()` when the tensor structure stays fixed and only symbolic parameters
change across batches. Supply constants, custom functions, and parameters explicitly. Call
`evaluate()` for real inputs and `evaluate_complex()` when coefficients or inputs are complex.
`compile()` creates source and a shared library, so treat its filenames, compiler, architecture,
and optimization level as reproducibility inputs rather than incidental arguments.

Control Symbolica-backed Rayon work with `SymbolicParallelism`: `Serial` keeps work on the
calling thread, `Auto` checks license capability and applies workload heuristics, and `Parallel`
forces Rayon without that safety choice. Configure the policy before benchmarking; otherwise a
threading-policy change can be mistaken for an algorithmic improvement.

For exact parameters and defaults, use the
#link("reference/python/spynso3/Tensor/#exports-tensor-evaluator-method")[evaluator reference] and
#link("reference/python/spynso3/set_symbolica_rayon_enabled/")[symbolic parallelism reference].
The implementation starts in
#source-link("crates/spynso3/src/lib.rs", label: "the Spenso Python adapter").
]
