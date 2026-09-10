= Spenso implementation architecture

#quote(block: true)[
#strong[Status:] Current implementation architecture, audited against the Spenso source on
2026-08-18.

This note covers the `spenso` Rust crate. `spenso-macros`, `spenso-hep-lib`, and `spynso3` are
separate packages: they provide derives, concrete physics tensors, and a Python adapter rather
than a second tensor engine.
]

== Boundaries

Spenso is an in-process, generic tensor and tensor-network library. Its implementation is split
into these layers:

- `structure`: representations, dimensions, abstract and concrete indices, slots, ordering,
  permutations, names, and structure merging;
- `tensors`: dense, sparse, symbolic, and—when Symbolica support is enabled—parametric payloads;
- `algebra`, `iterators`, and `contraction`: scalar promotion, indexed traversal, traces,
  exterior products, and single- or multi-index contraction kernels;
- `network`: a Linnet-backed expression graph, separate value store, tensor and function
  libraries, contraction planning, and sequential or parallel execution;
- feature-gated Symbolica integration in `shadowing`, `network::parsing`,
  `symbolic_parallelism`, and `symbolica_init`.

The crate does not own a global physics tensor catalogue. A caller supplies a `Library` for
external tensor keys and a `FunctionLibrary` for opaque functions. `spenso-hep-lib` is one
possible source of concrete tensors; Idenso is a separate identity-and-rewrite layer over
Spenso-compatible Symbolica expressions.

== Core representation

The structural kernel is `Slot<R, Aind>`, which pairs a `Representation<R>` with an abstract
index. A representation combines a representation name and a dimension. Its `RepName`
implementation owns duality, matching, ordering, and any negative-metric convention. Two slots
contract only when their dimensions and representation orientations match and their abstract
indices are equal after dualization.

`TensorStructure` exposes an ordered external slot sequence and the conversions derived from it:
shape, size, flat and expanded indices, duality, permutations, and repeated-index traces.
Concrete implementations choose different ownership policies: `OrderedStructure` owns an
ordered slot vector, `NamedStructure` adds a global name and arguments, `PermutedStructure`
records representation and index permutations, and `SmartShadowStructure` couples symbolic and
ordered views. `HasStructure` keeps this metadata separate from a tensor's payload.

The principal materialized payloads are:

- `DenseTensor<T, S>`, an element vector plus structure. Its checked constructor requires the
  vector length to equal the structure size, apart from the one-element scalar case;
- `SparseTensor<T, S>`, a map from `FlatIndex` to non-default entries, an explicit zero value,
  and the structure. Checked index-based construction converts validated concrete indices to
  flat positions;
- `DataTensor<T, S>`, a dense-or-sparse enum that preserves one structure and coefficient type
  across generic operations.

These fields are public, so the checked constructors are the invariant-preserving boundary;
callers constructing structs literally are responsible for keeping payload extent and structure
consistent.

== Direct contraction flow

`Contract` first delegates topology to `StructureContract::merge`. The merge identifies slots
that match the other structure's dual slots, removes contracted pairs, and returns both operand
subsets plus `MergeInfo`, which records clean operand order or an interleaved result.

```text
left structure + right structure
  -> merge dual slots and construct result structure
  -> zero matches: exterior product
  -> one match: single contraction
  -> multiple matches: multi contraction
  -> apply normal, reversed, or interleaved data kernel
  -> tensor with the merged external structure
```

Coefficient types participate through fallible multiplication, accumulation, subtraction, and
promotion traits. Representation-specific negative entries can therefore change accumulation
signs without changing the structural matching rule. Dense–dense, dense–sparse, and
sparse–dense `DataTensor` contractions produce dense output; sparse–sparse preserves sparse
output. Internal pairs on one tensor use `Trace::internal_contract` instead of the two-operand
path.

== Network representation and execution

`Network<S, LibKey, FunKey, Aind>` owns three independent pieces of state:

- `NetworkGraph` owns expression topology in a Linnet half-edge graph. Nodes are either leaves or
  `Sum`, `Neg`, `Product`, `Power`, and `Function` operations. Head edges encode the expression
  tree; slot edges encode contractions between dual slots;
- the generic store `S` owns local tensor and scalar values. The standard `NetworkStore<T, Sc>`
  uses stable vector positions during an execution and keeps optional scalar aliases alongside
  the source scalars;
- `NetworkState` summarizes whether the expression is a pure scalar, scalar, tensor, or
  self-dual tensor so network algebra can maintain compatibility.

A leaf refers to a local store position, a lazy tensor sum or scaled tensor, a scalar position,
or a caller-owned library key with its current indices. The graph never embeds the external
library value. Materialization occurs through the supplied `Library`, while opaque function
nodes delegate to the supplied `FunctionLibrary`.

Execution follows one semantic path with interchangeable policies:

```text
Network::execute
  -> merge adjacent expression operations
  -> trace slot self-loops on leaves
  -> find ready operation subgraphs
  -> materialize library leaves when needed
  -> execute scalar, sum, power, function, or product operation
  -> choose product pairs through ContractionStrategy
  -> replace each completed subgraph with a result leaf
  -> refresh NetworkState
```

`ExecutionStrategy` controls graph scheduling; `ContractionStrategy` controls pair selection
inside products. `Sequential`, bounded `Steps`, and extraction-based variants share the same
`ExecuteOp` boundary. `Parallel` evaluates independent ready operations against read-only base
storage plus per-worker overlays, then rebases newly appended store indexes before graph
replacement. It deliberately falls back to sequential execution when Symbolica's Rayon support
is disabled or the contraction strategy performs partial graph rewrites.

=== Symbolic parallelism and optimized sums

With `shadowing`, the process-wide `SymbolicParallelism` policy controls Rayon work involving
Symbolica atoms. `Auto` checks the Symbolica license once and uses a workload heuristic where an
operation provides one; `Serial` disables Rayon; `Parallel` forces Rayon and bypasses `Auto`'s
license safety check. Configure the policy before tensor work: changing it concurrently with an
active operation is unsupported.

The network `Parallel` execution strategy still falls back to `Sequential` when that policy
resolves to serial or when the selected contraction strategy performs partial graph rewrites.
Separately, `FastTensorSum`, `FastTensorSumContractible`, balanced scalar/tensor sums,
contraction profiles, and pair estimates provide optional optimized paths inside the same
network semantics. Unsupported shapes or unprofitable cases use the ordinary sum and
contraction paths. Network counters and timers are opt-in diagnostic instrumentation; they do
not change the execution result.

== Symbolica parsing boundary

The `shadowing` feature enables `network::parsing`. `NetworkParse` and the more general
`Network::try_from_view` lower a Symbolica `AtomView` into the same graph and store used by
programmatic construction. Dispatch recognizes addition, multiplication, integer powers, and
functions before falling back to scalar or opaque-tensor construction.

`ParseSettings` is semantic policy, not formatting: it controls scalar precontraction, first-term
structure discovery, recursion depth, shorthand expansion versus opaque leaves, composite scalar
handling, and how strictly a function must be tagged as a tensor. The tensor library resolves
known keys; `TensorFromExpression` owns construction of an opaque tensor leaf; the function
library owns supported opaque operations. Fresh dummy indices come from a parse-local
`ParseState`, so callers combining independently parsed expressions must still manage index
namespaces deliberately.

The detailed dispatch and shorthand behavior are recorded in the
#link("parsing-flow.typ")[Symbolica-to-network parsing flow]. Syntax and rewrite ownership across
Spenso and Idenso are recorded in the
#link("spenso-symbolica-syntax-and-rewrites.typ")[Symbolica syntax and rewrite boundary].

== Features and serialization

The core crate has no default feature declaration. `shadowing` enables Symbolica atoms,
parametric tensors, parsing, state-aware decoding, and Linnet's Symbolica adapter. `python`
only adds PyO3 extraction for Spenso's complex scalar type; the full Python tensor surface lives
in `spynso3`. `python_stubgen` adds stub metadata and implies `python`.

Structures, dense and sparse tensors, `DataTensor`, `NetworkGraph`, `NetworkStore`, and
`Network` derive Serde and bincode representations when their generic fields support them.
Under `shadowing`, selected types also decode through Symbolica's state map. These are type-level
encodings, not a versioned checkpoint protocol: Spenso provides no state directory, database,
schema migration, or cross-version compatibility promise. The caller owns byte storage, format
versioning, Symbolica-state availability, and library-key reconstruction.

== Ownership and error boundaries

Structure and direct-contraction failures are represented by `SlotError`, `StructureError`, and
`ContractionError`. Network construction and execution widen these into `TensorNetworkError`,
which also covers graph invariants, library lookup, function-library application, incompatible
sums, unsupported powers, incomplete contraction, and I/O or opaque errors. Parsing returns that
same network error type after converting coefficient and structure failures.

The store owns local results; graph leaves only own indexes into it. External libraries and
function libraries are borrowed for an execution and remain caller-owned. Consequently a
serialized graph containing `LibraryKey` leaves is not self-contained unless the same key
semantics and values are restored. A completed network can expose zero, one, a scalar, or a
tensor through its result accessors; callers should use those checked boundaries rather than
assuming the final leaf kind.

Some low-level APIs retain programmer-invariant panics—for example literal structs with invalid
shape, inconsistent ordered structures, or an integer power of a non-self-dual tensor. Public
fallible constructors, parsing entry points, and execution methods are the intended boundaries
for untrusted or dynamically inferred inputs.

== Maintained invariants

- Slot contraction requires equal abstract indices and matching dual representations at equal
  dimensions; matching names alone is insufficient.
- Tensor payload order and flat-index interpretation come from the attached structure.
- A network's head edges form expression dependencies, while slot edges carry tensor
  contraction topology; the two edge roles are not interchangeable.
- Local leaf indexes refer to the associated store. Joining stores or committing parallel
  overlays must rebase every affected tensor and scalar reference.
- A library leaf owns a key and requested index order, not the library tensor. Permutations are
  applied when the value crosses the library boundary.
- Ready-operation replacement preserves the boundary slot edges of the collapsed subgraph and
  finishes Linnet's deferred node identifications before later scheduling.
- Parse depth and shorthand settings can change graph granularity without changing the external
  tensor structure that opaque leaves report.

== Verification and related documentation

Unit and snapshot tests live beside structures, dense and sparse data, contraction kernels,
network graph and store operations, libraries, execution strategies, symbolic parsing,
parametric tensors, and serialization-sensitive types. The default typed boundary is exercised
with `cargo test -p spenso`; parsing and symbolic execution require
`cargo test -p spenso --features shadowing`. Rayon behavior also depends on the Symbolica
parallelism setting, so performance checks should record both the Cargo feature set and runtime
thread configuration.

For supported workflows, start with the
#link("../../../products/spenso/latest/tutorial/")[first-contraction tutorial], continue with the
#link("../../../products/spenso/latest/guides/networks/")[network guide], and use the
#link("../../../products/spenso/latest/reference/rust/spenso/")[native Spenso Rustdoc]
for exact public signatures. Python construction and execution are documented in the
#link("../../../products/spenso/latest/guides/python/")[Spynso3 workflow], whose ownership is
separate from the core crate described here.
