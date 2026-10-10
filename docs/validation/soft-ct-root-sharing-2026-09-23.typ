= Reducing the production expression before global factor collection

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

The full ordinary-U root is already 79,062,885 bytes (75.40 MiB) before the
protection pass in `FinalIntegrands::into_integrands`. After protection it is
85,217,462 bytes (81.27 MiB). The earlier 18,599-byte figure described retained
numerator bodies in one Taylor event; it was not the assembled root size.
The Symbolica allocation failure and excessive duplication in our input are
therefore separate, concrete opportunities for improvement.

== Direct measurements

A bounded debugger run reads the root and protected root from the preserved
`2c2694036e93` stack executable while stopped at Symbolica 3.0.0
`collect.rs:1343`. This uses the same three-loop graph and full-U card as the
previous investigations. No production code or algebraic settings change.

#table(
  columns: (auto, auto, auto),
  [Property], [Before protection], [After protection],
  [Packed root bytes], [79,062,885], [85,217,462],
  [Outer summands], [1,030], [1,030],
  [Direct factor occurrences], [9,160], [9,160],
  [Distinct literal direct factors], [799], [7,111],
  [Sum-factor occurrences], [886], [886],
  [Bytes in those sum factors], [78,769,329], [84,827,149],
)

Nested sum factors account for 99.6% of the unprotected root's direct-factor
bytes. The largest is 2,014,327 bytes and contains only four immediate summands:
substantial structure sits below the outer sum. Deduplicating only complete
outer factors would still leave 70,385,061 unique factor bytes.

Protection assigns 238,366 occurrence tags. The tags intentionally keep inverse
and numerator/selector occurrences local, but also make previously identical
surrounding sums distinct. This increases distinct literal direct factors from
799 to 7,111. Symbolica's recursively collected factor dictionary has 7,112
keys; literal direct-factor counting is a different measurement.

== Exact structural sharing

Interning identical encoded subtrees of the unprotected root finds only
19,530 distinct nodes and 75,575 edges between those nodes. An experimental
immutable-node DAG file stores the same structure in 586,129 bytes (0.56 MiB),
about 135 times smaller than the packed tree. This is node sharing, not gzip,
factorization, cancellation or expansion.

The saved DAG was independently decoded; reconstructing the original byte
stream gives exactly the original 79,062,885 bytes and SHA-256 digest. Both
representations exclude the separately stored symbol metadata and retained
numerator-definition catalog. This establishes repeated structure in this root.
It does not measure application RSS or prove that evaluator aliases can safely
reuse every node across tensor bindings and branch guards.

== Where to preserve the structure

The existing `DirectResidueBranches` and `Integrands` types retain numerator
bodies separately, but scalar roots remain owned `Atom` trees. Taylor projection
returns those trees through `series_preserving_factors`. Final normalization
maps both the shared bodies and each root. `materialize()` then sums residue
branches for each forest node; `orientation_parametric_exprs()` sums the forest
nodes, splits remaining momenta, and invokes global factor collection.

The evaluator already has `AliasedAtom`, joint numerator-family construction,
and reachable-definition pruning. They run after the failing collection call.
This gives existing owners to extend instead of introducing a second unrelated
expression system. The temporary protection function currently embeds its
whole argument; it is not an out-of-line shared definition.

A reduction should preserve repeated scalar subexpressions as explicit
references through Taylor, final normalization and signed-forest assembly,
then pass those definitions into the existing evaluator machinery. Algebra
should operate on small roots and shared bodies, rather than repeatedly
reconstructing their full trees. Keep Taylor-dependent quantities visible to
series operations, preserve tensor interfaces and dummy bindings, and retain
full residue-map keys and selector-local inverse evaluation. Cancellations
needed before the final momentum split must still see the complete signed sum.

No extra global `.expand()`, `together()` or indiscriminate denominator merging
is justified by these measurements. Batching the existing forest `zip_add`
would reduce repeated construction copies, but would not by itself shrink the
final root. The remaining investigation is to measure scalar-root growth at
successive Taylor and normalization boundaries and prove a scope-preserving
retained representation on the smallest nonzero orientation before full U/H.

== Evidence and limits

The structural metrics above and the adjacent
#link("soft-ct-root-sharing-2026-09-23/root-structure.dag")[lossless DAG]
are retained. Original packed root dumps, debugger inputs, logs and receipts
were local campaign artifacts and are not distributed. The host-specific
capture and audit scripts have been removed.

The DAG begins with the eight bytes `SCFDAG1\0`, then little-endian `u32` root
ID and node count. Each node stores `u32` prefix length and child count, its
literal prefix bytes, then the ordered `u32` child IDs. Recursively emitting a
node's prefix followed by its children reconstructs the original packed byte
stream. No Symbolica symbol state or numerator-definition catalogue is attached;
this is structural evidence, not an evaluator input.

The initial debugger probe exposed invalid optimized-frame fields for the
definition store; those values are excluded. Only the subsequently verified
root/protected byte lengths are reported.
