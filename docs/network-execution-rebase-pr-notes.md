# Network Execution Rebase PR Notes

This file tracks the code changes being preserved while rebasing the network
execution work onto `symbolic-network-simplify-v2`. It is source material for
the PR description: the focus is what the new code does and how it does it.
Commit labels are kept only for traceability during the rebase.

## Integrated Changes

### Direct Boundary Contractions Over Sums

Trace commit: `pmrsnlov`.

What this adds:

- Allows the network-aware Schoonschip path to apply dummy-index contractions
  directly at a product boundary such as `p(mu) * sum(...)` or
  `g(mu, nu) * sum(...)`.
- Avoids sending these simple boundary contractions through the generic
  multiplication-pattern machinery when the network analysis already found the
  replacement.
- Keeps expansion as a fallback, not the normal path, for cases where residual
  contracted slots still remain after the direct replacement.

How it works:

- The network contraction analysis produces a replacement map for the dummy
  index/vector/metric relation seen at the product boundary.
- The replacement map is applied transitively before deciding whether the term
  still needs distribution.
- The direct replacement target is compact-cleaned before falling back to
  distribution, so already-resolved contractions do not force larger algebra.
- The fallback path distributes only terms that still contain the contracted
  slots after the direct rewrite.
- `IDENSO_TRACE_FINISH_CONTRACTS` can be used to inspect which contractions are
  finished by this direct path.

Validation run:

- `cargo fmt`
- `cargo check --package idenso --test network_dumb --profile dev-optim`

### Full-Depth Parse And Execution Profiling

Trace commit: `orwppqxy`.

What this adds:

- Adds diagnostic coverage for the large Spenso text inputs without prior
  simplification, so parsing and execution can be measured independently.
- Adds network profiling counters and timers for the expensive orchestration
  paths seen in the large inputs.
- Establishes a CI/test baseline document for the network execution work.

How it works:

- Large-input tests strip Symbolica tags, initialize idenso, parse the raw
  expression into symbolic networks, and then exercise execution separately.
- `crates/spenso/src/network/profile.rs` exposes runtime-gated timers and
  counters.
- The parser records profile spans inside the current modular parser in
  `crates/spenso/src/network/parsing/`.
- Network execution records graph joins, graph extraction, graph state scans,
  store extensions, operation execution, ready-operation discovery, and merge
  operations.
- Profiling is gated by `SPENSO_NETWORK_PROFILE`. With it disabled, expensive
  collection is skipped, but hot paths still pay a small call/branch cost.

Validation run:

- `cargo fmt`
- `cargo check --package idenso --test large_spenso_inputs --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`
- `cargo check --package spenso --tests --profile dev-optim`

### Live Operation Descriptors

Trace commit: `mkswtoln`.

What this adds:

- Changes network execution from "extract an operation graph, execute it, splice
  back a replacement graph" toward "describe a ready operation in the live graph
  and return the replacement leaf".
- Introduces `NetworkOperation` as the owned descriptor for a ready operation.
- Changes `ExecuteOp` so an operation executor receives the borrowed live
  `NetworkGraph` plus the `NetworkOperation`.
- Sets up later execution improvements where contracted edges can be hidden in
  an ignored subgraph instead of immediately deleting and invalidating graph
  state.

How it works:

- Ready-operation discovery builds a `NetworkOperation` that records the
  operation kind and the live subgraph to be collapsed.
- `ExecuteOp::execute(...)` runs the operation against the live graph
  descriptor and returns a `NetworkLeaf<K, Aind>` for the result.
- The network collapses the operation subgraph by identifying its nodes into
  the returned replacement leaf.
- Self-edges created by identifying operation nodes are collected in an ignored
  subgraph and deleted after the batch/iteration boundary.
- This commit still had a temporary `clone_subgraph(...)` compatibility helper
  so older product-contraction strategies could keep using graph-shaped input
  until the next change removed that adapter.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package idenso --tests --profile dev-optim`

### Product Contraction Without Graph Cloning

Trace commit: `xwqrolxv`.

What this adds:

- Replaces product contraction by cloned subgraph construction with a compact
  product-local operand state.
- Makes `ContractionStrategy` operate on the borrowed live graph plus a
  `NetworkOperation`, and return the replacement leaf directly.
- Removes the temporary `clone_subgraph(...)` helper from the previous change.
- Removes the `replacement_leaf(...)` graph-to-leaf adapter because strategies
  now return the final leaf themselves.

How it works:

- `ProductContraction::from_operation(...)` reads the operation subgraph and
  builds a list of `ProductOperand<K, Aind>` values.
- Product operands point at local tensors, scalars, or library keys. Library
  operands retain their source node until they must be materialized.
- Product-local helpers perform scalar folding, library materialization,
  best-pair selection by matching degree, pair contraction, and final collapse.
- Built-in strategies such as `ContractScalars`, `SmallestDegree`,
  `SmallestDegreeIter`, `SingleSmallestDegree`, and `SingleLargestDegree` now
  use that product operand state directly.
- The idenso Schoonschip strategies use the same product state and add
  expression-size based ordering metrics for symbolic tensor products.
- The live graph is mutated only after a strategy has produced the replacement
  leaf for the whole ready operation.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package idenso --tests --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`

### Optional Rayon Ready-Batch Execution

Trace commit: `prszmwmv`.

What this adds:

- Adds an opt-in `Network::execute_parallel(...)` path for symbolic networks.
- Keeps the existing sequential execution path as the default.
- Executes independent ready operations from the same ready batch in parallel
  with Rayon.
- Avoids cloning the full network store per worker by giving workers overlay
  stores that borrow the shared base store and own only worker-created tensor
  and scalar additions.
- Adds diagnostic coverage comparing sequential and parallel symbolic execution
  on `symbolica_expression.txt`.

How it works:

- `NetworkStoreAccess` abstracts over both `NetworkStore<T, Sc>` and
  `NetworkStoreOverlay<'_, T, Sc>`.
- Contraction helpers use `NetworkStoreAccess` for tensor/scalar lookup and for
  pushing newly created tensors/scalars. This lets the same contraction strategy
  run in the owning store and in worker overlays.
- `NetworkStoreOverlay` indexes into the borrowed base vectors for existing
  tensors/scalars and into local vectors for additions created by the worker.
- The parallel loop first plans a ready-operation batch from the live graph.
- A single ready operation still runs sequentially to avoid Rayon overhead.
- For a multi-operation batch, each operation executes in a worker-local
  overlay. Each worker returns the replacement leaf plus its local tensor/scalar
  additions.
- The main thread appends each worker's additions to the owning store and
  remaps replacement leaf indices from overlay-local offsets into final store
  offsets.
- After remapping, the graph identifies each operation subgraph into its
  replacement leaf and hides merge-created self-edges in the ignored subgraph
  until the iteration boundary.
- The parallel path is rebased onto the current `ProductContraction<K, Aind>`
  implementation, so it uses the same product-local operand state as
  sequential execution rather than reintroducing graph clone/splice execution.
- The large-input diagnostics use the current parser API
  (`StructureFromAtom`, `ParseSettings`) and compare sequential and parallel
  symbolic execution results.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package idenso --tests --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`
