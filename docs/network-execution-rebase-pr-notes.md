# Network Execution Rebase PR Notes

This file tracks the concrete changes preserved while rebasing the network
execution work onto `symbolic-network-simplify-v2`. It is intended as source
material for the eventual PR description.

## Resolved Changes

### `pmrsnlov` - avoid expanding rewritten direct contraction target sums

Resolved commit: `479ee4e9`.

Actual changes preserved:

- Kept the direct contraction rewrite path for Schoonschip-style dummy
  contractions over expanded sums.
- Preserved transitive replacement application before deciding whether a term
  still needs distribution.
- Kept the compact-cleaning step on the direct replacement target before
  falling back to distribution.
- Preserved the fallback path that distributes only when residual slots remain.
- Kept the `IDENSO_TRACE_FINISH_CONTRACTS` trace helper in
  `schoonschip/utils.rs`.
- Updated the added `network_dumb` calls to the current
  `schoonschip_with_net::<false, AbstractIndex>(&settings)` signature.

Resolution choices:

- Kept the destination bounds that require `DummyAind + 'static`.
- Did not pull later local-boundary cleanup changes into this commit.

Verification:

- `cargo fmt`
- `cargo check --package idenso --test network_dumb --profile dev-optim`

### `orwppqxy` - profile raw full-depth Spenso parsing and establish baseline

Resolved commit: `ddeab59a`.

Actual changes preserved:

- Kept the large Spenso input diagnostic tests in
  `crates/idenso/tests/large_spenso_inputs.rs`.
- Kept the profiling module at `crates/spenso/src/network/profile.rs`.
- Ported parse profiling into the current modular parser at
  `crates/spenso/src/network/parsing/mod.rs`.
- Kept network/store/graph profiling counters and timers for parse, graph join,
  graph extraction, graph state scans, store extension, execution operations,
  and merge operations.
- Kept the test-baseline document at
  `docs/architecture/spenso-network-execution-test-baseline.md`.
- Removed the obsolete conflicted flat parser file
  `crates/spenso/src/network/parsing.rs`; the destination parser is now the
  `network/parsing/` module tree.

Resolution choices:

- Kept the newer parser semantics: `StructureFromAtom`,
  `TensorFromExpression`, `ShorthandParsing`, and
  `ShadowedStructure::<Aind>::parse(...)`.
- Kept `Network::n_mul` using `graph.state()` for tensor-state recomputation,
  instead of the older `dangling_indices().is_empty() => Scalar` shortcut.
- Kept `library_tensor` storing the explicit external indices on
  `NetworkLeaf::LibraryKey`.
- Fixed the large-input test import to use `StructureFromAtom` instead of the
  removed `Parse` trait.

Verification:

- `cargo fmt`
- `cargo check --package idenso --test large_spenso_inputs --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`
- `cargo check --package spenso --tests --profile dev-optim`

Note:

- Profiling is runtime-gated by `SPENSO_NETWORK_PROFILE`. With profiling off,
  expensive work is skipped, but hot paths still pay a small call/branch cost.
  A compile-time feature or macro layer would be needed for truly zero-cost
  disabled profiling.

### `mkswtoln` - execute operations through live graph operation descriptors

Resolved commit: `74698466`.

Actual changes preserved:

- Kept `NetworkOperation` as the owned descriptor for a ready operation.
- Changed `ExecuteOp` to execute against a borrowed live `NetworkGraph` plus a
  `NetworkOperation`, returning a replacement leaf instead of a full replacement
  graph.
- Added `clone_subgraph(...)` as the compatibility adapter needed by existing
  graph-based product contraction strategies.
- Collapsed operation subgraphs by identifying their nodes into a replacement
  leaf.
- Accumulated merge-created self-edges in an ignored subgraph before deleting
  them.
- Kept the destination `Aind`-aware network graph types and adjusted the new
  execution path to return `NetworkLeaf<K, Aind>`.

Resolution choices:

- Kept `find_all_ready_ops` returning
  `Vec<(Self, NetworkOp<FK>, Vec<NetworkLeaf<K, Aind>>)>`.
- Changed deferred node identification to accept `NetworkNode<K, FK, Aind>`.
- Preserved the profiling added by `orwppqxy`.

Verification:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package idenso --tests --profile dev-optim`

### `xwqrolxv` - contract products through live graph operation descriptors

Resolved in working change over commit `d15fa8b1`.

Actual changes preserved:

- Changed `ContractionStrategy` to receive the borrowed live
  `NetworkGraph` plus the ready `NetworkOperation`, and to return a
  replacement leaf directly.
- Replaced product contraction by graph cloning with a lightweight
  `ProductContraction` operand state.
- Removed the temporary `clone_subgraph(...)` compatibility path introduced in
  `mkswtoln`.
- Added product-local operations for scalar folding, library materialization,
  best-pair selection by matching degree, pair contraction, and final product
  collapse.
- Rebased the built-in product strategies (`ContractScalars`,
  `SmallestDegree`, `SmallestDegreeIter`, `SingleSmallestDegree`, and
  `SingleLargestDegree`) onto `ProductContraction`.
- Rebased the idenso Schoonschip ordering strategies onto
  `ProductContraction`, including expression-size based order scoring.
- Kept the destination `Aind`-aware `NetworkLeaf<K, Aind>` and
  `NetworkLeaf::LibraryKey { .. }` shapes throughout.

Resolution choices:

- Chose the `xwqrolxv` product operand-state implementation over the previous
  graph-extraction/splice implementation.
- Kept the live-operation descriptor API from `mkswtoln`.
- Removed the `replacement_leaf(...)` adapter because contraction strategies
  now return the replacement leaf directly.
- Kept profiling from the earlier profiling change intact.

Verification:

- `cargo fmt`
- `cargo check --package spenso --tests --profile dev-optim`
- `cargo check --package idenso --tests --profile dev-optim`
- `cargo check --package gammalooprs --tests --profile dev-optim`

## Current Conflict

No current conflict recorded in this file yet. Update this section after
advancing to the next conflicted commit.
