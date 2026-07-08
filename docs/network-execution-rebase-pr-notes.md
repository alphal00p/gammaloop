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
  directly at product boundaries such as `p(mu) * sum(...)` or
  `g(mu, nu) * sum(...)`.
- Avoids sending these simple boundary contractions through the generic
  multiplication-pattern machinery when the network analysis already found the
  replacement.
- Keeps expansion as a fallback, not the normal path, for terms that still have
  residual contracted slots after the direct replacement.

How it works:

- The network contraction analysis produces a replacement map for the
  dummy-index/vector/metric relation seen at the product boundary.
- The replacement map is applied transitively before deciding whether a term
  still needs distribution.
- The direct replacement target is compact-cleaned before falling back to
  distribution.
- `IDENSO_TRACE_FINISH_CONTRACTS` can be used to inspect which contractions are
  finished by this direct path.

Validation run:

- `cargo fmt`
- `cargo check --package gammalooprs --test large_spenso_actual --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`
- `cargo check --package idenso --test network_dumb --profile dev-optim`

### Full-Depth Parse And Execution Profiling

Trace commit: `orwppqxy`.

What this adds:

- Adds diagnostic coverage for the large Spenso text inputs without prior
  simplification, so parsing and execution can be measured independently.
- Adds network profiling counters and timers for expensive orchestration paths.
- Establishes a CI/test baseline document for the network execution work.

How it works:

- Large-input tests strip Symbolica tags, initialize idenso, parse the raw
  expression into symbolic networks, and exercise execution separately.
- `crates/spenso/src/network/profile.rs` exposes runtime-gated timers and
  counters.
- The parser records profile spans inside the current modular parser in
  `crates/spenso/src/network/parsing/`.
- Network execution records graph joins, graph extraction, graph state scans,
  store extensions, operation execution, ready-operation discovery, and merge
  operations.
- Profiling is gated by `SPENSO_NETWORK_PROFILE`. With profiling off, expensive
  collection is skipped, but hot paths still pay a small call/branch cost.

Validation run:

- `cargo fmt`
- `cargo check --package idenso --test large_spenso_inputs --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`
- `cargo check --package spenso --tests --profile dev-optim`

### Live Operation Descriptors

Trace commit: `mkswtoln`.

What this adds:

- Changes execution from "extract an operation graph, execute it, splice back a
  replacement graph" toward "describe a ready operation in the live graph and
  return the replacement leaf".
- Introduces `NetworkOperation` as the owned descriptor for a ready operation.
- Changes `ExecuteOp` so an operation executor receives the borrowed live
  `NetworkGraph` plus the `NetworkOperation`.

How it works:

- Ready-operation discovery builds a `NetworkOperation` that records the
  operation kind and live subgraph to collapse.
- `ExecuteOp::execute(...)` runs against the live graph descriptor and returns
  a `NetworkLeaf<K, Aind>`.
- The network identifies operation-subgraph nodes into the returned replacement
  leaf.
- Merge-created self-edges are collected in an ignored subgraph and deleted at
  the iteration boundary.
- This step still had a temporary `clone_subgraph(...)` compatibility helper,
  removed by the next product-contraction change.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package idenso --tests --profile dev-optim`

### Product Contraction Without Graph Cloning

Trace commit: `xwqrolxv`.

What this adds:

- Replaces product contraction by cloned subgraph construction with compact
  product-local operand state.
- Makes `ContractionStrategy` operate on the borrowed live graph plus a
  `NetworkOperation`, returning the replacement leaf directly.
- Removes the temporary `clone_subgraph(...)` helper and the
  `replacement_leaf(...)` graph-to-leaf adapter.

How it works:

- `ProductContraction::from_operation(...)` reads the operation subgraph and
  builds `ProductOperand<K, Aind>` values.
- Product operands point at local tensors, scalars, library keys, and later
  tensor sums. Library operands retain their source node until materialization.
- Product-local helpers perform scalar folding, library materialization,
  best-pair selection, pair contraction, and final collapse.
- Built-in strategies and idenso Schoonschip strategies use the product operand
  state directly.
- The live graph is mutated only after a strategy produces the replacement leaf
  for the ready operation.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package idenso --tests --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`

### Optional Rayon Ready-Batch Execution

Trace commit: `prszmwmv`.

What this adds:

- Adds an opt-in `Network::execute_parallel(...)` path for symbolic networks.
- Keeps sequential execution as the default path.
- Executes independent ready operations from the same ready batch in parallel
  with Rayon.
- Avoids cloning the full network store per worker by using overlay stores.
- Adds diagnostics comparing sequential and parallel symbolic execution on
  `symbolica_expression.txt`.

How it works:

- `NetworkStoreAccess` abstracts over both `NetworkStore<T, Sc>` and
  `NetworkStoreOverlay<'_, T, Sc>`.
- Worker overlays borrow base tensors/scalars and own only worker-created
  additions.
- Each worker returns a replacement leaf plus local tensor/scalar additions.
- The main thread appends additions, remaps replacement indices, and identifies
  each ready-operation subgraph into its replacement leaf.
- The parallel path uses the same `ProductContraction<K, Aind>` operand state
  as sequential execution.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package idenso --tests --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`

### Concrete Tensor Parse/Execution Diagnostics

Trace commit: `mlnonqyu`.

What this adds:

- Adds ignored GammaLoop diagnostics for parsing and executing
  `symbolica_expression.txt` through the evaluator `ParsingNet` path.
- Adds a Hornered variant of the same diagnostic.
- Adds a small concrete HEP tensor probe using non-symbolic `hep_lib` data.
- Adds profile-only logging around sequential ready-operation batches and slow
  product-pair contractions.

How it works:

- Diagnostics initialize GammaLoop, strip Symbolica tags, parse the expression,
  and build a `ParsingNet` using GammaLoop's tensor library.
- Network statistics count concrete tensors, parametric tensors, scalars,
  logical entries, stored entries, tensor orders, and common tensor names.
- Step diagnostics execute `Steps<1>` repeatedly and emit profile reports.
- The concrete HEP probe uses `spenso_hep_lib::hep_lib(...)`.
- Product contraction logs selected operand pairs and reports slow pair
  contractions around `ProductContraction::contract_pair(...)`.
- The old `parse_inner_products = false` intent is expressed through current
  opaque shorthand parsing with fast structure inference.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package gammalooprs --test large_spenso_actual --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`

### Contraction Order MWE And MinResultRank Diagnostic

Trace commit: `nlzmnxqt`.

What this adds:

- Adds `MinResultRank`, a diagnostic contraction strategy that chooses the next
  pair by minimizing the rank/order of the intermediate contraction result.
- Adds ignored gamma-ladder MWE diagnostics comparing `SmallestDegree` to
  `MinResultRank`.
- Adds an equivalence check that expands the scalar difference between both
  contraction orders and asserts it is zero.

How it works:

- `ProductContraction::best_result_rank_pair(...)` scores candidate tensor
  pairs by `left.order + right.order - 2 * matching_degree`.
- Candidate pairs with no matching slots are treated as bad choices.
- Ties prefer higher matching degree.
- `contract_one_by_result_rank(...)` materializes library operands, selects the
  best result-rank pair, logs profile information, contracts the pair, and
  folds scalars.
- The strategy is rebased onto `ProductContraction<K, Aind>`: `Aind` stays at
  the impl/strategy level, helper calls add only `Store`, and the result is
  `NetworkLeaf<K, Aind>`.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package gammalooprs --test large_spenso_actual --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`

### Opt-In Lazy Tensor Sums

Trace commit: `womukytz`.

What this adds:

- Adds `NetworkLeaf::TensorSum(Vec<usize>)` as an internal leaf for
  tensor-valued sums whose terms should remain separate for diagnostic
  execution.
- Keeps lazy tensor sums opt-in behind `SPENSO_NETWORK_LAZY_TENSOR_SUMS`; the
  default sum execution path still materializes tensor sums.
- Allows one-sided contractions such as `(T1 + ... + Tn) * U` to distribute
  through the lazy sum by contracting each term with `U`.
- Adds profiling for lazy tensor-sum creation, materialization, and distributed
  product contractions.

How it works:

- Sum execution first tries existing scalar and balanced tensor-sum fast paths.
- When lazy sums are enabled, the sum is tensor-valued, non-scalar, and has at
  least `MIN_LAZY_TENSOR_SUM_TERMS`, execution returns
  `NetworkLeaf::TensorSum(terms)` instead of adding all tensor entries.
- `ProductContraction` treats `TensorSum` as tensor-like: it reads structure
  from the first term, scales all terms under scalar multiplication, and expands
  the term list for pair contraction.
- If both product sides are lazy and the Cartesian product exceeds
  `MAX_LAZY_TENSOR_SUM_DISTRIBUTED_TERMS`, operands are materialized first.
- Negation maps over terms; functions and powers materialize tensor sums.
- Final result extraction rejects an unmaterialized `TensorSum`.
- The rebase keeps `LibraryKey { key, indices }` and carries `Aind` through
  lazy-sum helper signatures as `NetworkLeaf<K, Aind>`.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`

### Component Horner Diagnostics

Trace commit: `zlvxmnxm`.

What this adds:

- Adds diagnostics comparing raw and component-Hornered execution paths.
- Adds the type-level contraction setting `HornerAtomComponents<N>`, used as
  `MinResultRank<HornerAtomComponents<N>>`, to Horner generated tensor
  component atoms after contraction.

How it works:

- Parametric tensor contraction keeps the parsed graph unchanged, then maps the
  generated `Atom` entries through Horner collection when the contraction
  strategy requests it.
- The obsolete flat `crates/spenso/src/network/parsing.rs` stays deleted; only
  the modular parser remains.
- GammaLoop diagnostics use `ShorthandParsing::Opaque` with
  `StructureInferenceMode::Fast`; component Hornering is selected through the
  execution strategy type, not by mutating parse settings.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package gammalooprs --test large_spenso_actual --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`

### Staged Index Disconnection Diagnostic

Trace commit: `ntrwxyor`.

What this adds:

- Adds a diagnostic MWE that manually stages a contraction by disconnecting one
  boundary index before parsing.
- Executes the disconnected tensor-valued intermediate first, then reconnects
  it with an explicit metric in a second network execution.
- Adds `docs/architecture/spenso-large-expression-pathology.md`, summarizing
  the large-expression behavior and what the lazy-sum, component-Horner, and
  staged-disconnect experiments show.

How it works:

- The direct MWE parses and executes the fully connected expression.
- The staged variant parses an expression where one index is renamed so the
  corresponding contraction boundary is temporarily disconnected.
- The first stage executes to a tensor-valued `ActualTensor` intermediate.
- The second stage builds a network from that intermediate and multiplies it by
  a reconnecting metric network.
- The diagnostic compares the expanded scalar result against the direct path
  and asserts the difference is zero.
- The rebase keeps current parser imports for opaque shorthand settings and
  also imports `ShadowedStructure` for the `ActualTensor` alias.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package gammalooprs --test large_spenso_actual --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`

### Scalar-Preserving Tensor Sum Terms

Trace commit: `spwmxtox`.

What this adds:

- Extends lazy tensor sums from plain tensor indices to tensor terms that can
  carry an optional scalar factor.
- Keeps scalar factors outside numeric tensor addition and contraction for as
  long as possible.
- Adds fast-sum hooks and common-factor extraction hooks so parametric tensors
  can avoid materializing large scalar-expanded entries before contraction.

How it works:

- `TensorTerm { tensor, scalar }` records a local tensor plus an optional local
  scalar factor.
- `NetworkLeaf::TensorTerm` represents one scaled tensor term, and
  `NetworkLeaf::TensorTermSum` represents a sum whose terms may have different
  scalar factors.
- Product contraction now converts tensor-like operands into `TensorTerm`
  lists, combines scalar factors separately, contracts tensor payloads, and
  attaches the scalar product back to the resulting term.
- Pure unscaled term lists still use `TensorSum(Vec<usize>)`; scaled term lists
  use `TensorTermSum(Vec<TensorTerm>)`.
- Fallback paths materialize term sums only when an operation needs a concrete
  tensor, such as functions, powers, scalar result extraction, or over-large
  distributed products.
- The rebase keeps the generic leaf shape `NetworkLeaf<K, Aind>` and the
  structured library variant `LibraryKey { key, indices }`.

Validation run:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package idenso --tests --profile dev-optim`
- `cargo check --package gammalooprs --test large_spenso_actual --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`
