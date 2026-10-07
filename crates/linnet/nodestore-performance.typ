= Parent/child forest performance

Measured on `itphlies`, 7 October 2026, against the Linnet sources in
`1730ba8ba41b` and Symbolica Graphica `bbc2d73c70b4`. The implementation lives in
jj change `vpoqpqxtumps`, authored by lcnbr.

== Changes

- Cache the sibling boundary during parent/child preorder traversal. Subtree
  traversal now stops at its starting node instead of escaping into siblings.
- Construct parent/child storage directly from parent records, removing an
  intermediate allocation in the graph builder.
- Write a new edge's paired involution entries directly.
- Check active root records without walking identification chains.
- Flatten identification trees when `forget_identification_history()` explicitly
  discards history. Only the branches that identification attached below an
  active root since the previous cleanup are flattened; each discarded node record
  names one of them. Preserve half-edge IDs, payloads, involution, crown order,
  isolated nodes and active root identity; active node indices can still change
  through the existing root compaction. Cleanup scans the node records, as before,
  plus work proportional to the merged crowns, and subsequent root lookups need
  at most one parent hop. Half-edge swaps flatten the branches at both swapped
  positions, so cleanup still reaches every merged branch when a swap has made a
  record's label stale.

The parent/child flattening operation splices sibling rings in place and allocates
no temporary incidence vectors. Child-vector forests share the cleanup semantics:
reparenting appends, so they rebuild only the suffix of a root's child vector that
holds the merged branches, using one node bitset per cleanup, and rebuild a whole
crown when a swap meets a deep branch. Identification retains its hierarchy until
cleanup; Python defaults and rendering behavior are unchanged.

== Release measurements

Each size and store runs in a fresh process, pinned to one CPU. The timing harness
uses Rust 1.98.1, release optimization, thin LTO and one codegen unit. Medians cover
nine batches. Build and clone include destruction; traversal uses identical DFS
code for all three stores. Other workloads were active on the machine, so timings
below 100 microseconds are indicative only; the cleanup comparison below relies on
instruction counts instead.

A ring plus perfect matching has 1,000 vertices and 1,500 edges. Values are
microseconds per operation:

#table(
  columns: 5,
  [Operation], [Parent/child before], [Parent/child after], [Child-vector after], [Graphica],
  [Build], [102.22], [76.64], [290.70], [42.88],
  [DFS], [16.38], [14.45], [54.26], [7.04],
  [Adjacency scan], [11.85], [8.51], [51.93], [2.95],
  [Clone], [10.66], [10.57], [40.79], [33.47],
)

Parent/child construction improves by 25%, DFS by 12%, and adjacency scanning by
28% on this fixture. At 10,000 vertices construction improves from 986.66 to
751.46 microseconds and DFS from 168.94 to 147.26 microseconds.

Repeatedly identifying layers of a cycle produces 128 final vertices. At depth
32, the final graph has 4,096 edges, with parallel edges between adjacent vertices.
DFS after discarding history takes:

#table(
  columns: 2,
  [Store/state], [Microseconds],
  [Parent/child before], [376.65],
  [Parent/child after], [53.89],
  [Parent/child freshly rebuilt], [56.52],
  [Child-vector before], [502.16],
  [Child-vector after], [71.29],
  [Graphica freshly built], [12.49],
)

Cleanup removes the depth penalty: parent/child traversal improves about 7 times
and is comparable to rebuilding the same topology. Graphica has no corresponding
identification operation in this comparison, so its number measures traversal of
the resulting graph only.

== Batched cleanup

Batched tensor-network contraction identifies a batch of nodes and then discards
history once per batch. The `hub_merge` benchmarks merge every node of a square
grid into the first one in row order, discarding history after every 64 merges or
only once at the end. Flattening every active root at each cleanup, which an
earlier revision of this change did, costs linear time in all half-edges per
batch; flattening only the merged branches avoids that.

Instruction counts come from Callgrind on single runs of the whole merge loop,
including identification and node lookups; graph construction and cloning are
excluded. Values are thousands of instructions:

#table(
  columns: 5,
  [Workload], [Store], [Before], [Flatten all roots], [Flatten merged branches],
  [960 half-edges, batches of 64], [Parent/child], [87], [116], [86],
  [960 half-edges, batches of 64], [Child-vector], [90], [313], [157],
  [9,800 half-edges, batches of 64], [Parent/child], [2,640], [6,786], [1,672],
  [9,800 half-edges, batches of 64], [Child-vector], [2,781], [22,262], [2,419],
  [9,800 half-edges, one cleanup], [Parent/child], [620], [658], [750],
  [9,800 half-edges, one cleanup], [Child-vector], [612], [1,436], [1,357],
)

"Before" is `1730ba8ba41b`, which does not flatten. Parent/child allocation counts
are identical in all three versions. For child-vector forests at 9,800 half-edges,
batched merging allocates 51, 3,538 and 369 times respectively, and a single
cleanup 12, 2,537 and 28 times.

The batched cases at 9,800 half-edges take more than 100 microseconds per call, so
wall-clock ratios are meaningful there. Over five pinned processes of 31 calls
each, the median per-call time relative to "Before" is:

#table(
  columns: 3,
  [Store], [Flatten all roots], [Flatten merged branches],
  [Parent/child], [6.3 (6.2--6.3)], [0.78 (0.76--0.79)],
  [Child-vector], [7.2 (7.1--7.4)], [0.84 (0.83--0.85)],
)

Ranges divide the smallest and largest process median by the "Before" median.
At 9,800 half-edges batched merging is therefore faster than before flattening
was introduced, because later lookups stay one hop from their root. At 960
half-edges this holds for parent/child (0.99 times the "Before" instructions), but
child-vector forests still take 1.73 times the "Before" instructions and 45
instead of 11 allocations: rebuilding child vectors is not repaid by the few
lookups of so short a run. The difference is 66 thousand instructions, too small
for a wall-clock conclusion. The remaining cleanup cost also includes the existing
scan of all node records for discarded ones; making that proportional to the
merged nodes would require recording them in the forest.

A single cleanup after many merges does some work per discarded record instead of
one pass over every crown. For parent/child, the whole merge loop on the
9,800-half-edge grid takes 1.14 times the instructions of flattening all roots.
Discarding the history of the depth-32 layered cycle above takes 538 thousand
instructions, against 373 thousand for flattening all roots and 1,640 thousand
before. Child-vector forests take 0.95 and 0.83 times the instructions of
flattening all roots in these two cases. Both cases run below 100 microseconds per
call, so no wall-clock conclusion is drawn. Both strategies produce the same flat
crowns for merged nodes, so traversal after cleanup is unchanged.

Unlike flattening all roots, cleanup no longer flattens trees that only other
operations made deep, such as rerooting, unless a node was merged into them or a
half-edge swap touched them.

== Half-edge swaps

Discarded node records name their merged branches by half-edge label. A half-edge
swap (`permute`, `delete`, extraction) exchanges two labels without updating
those records, so a later cleanup would flatten the wrong branch. Each swap
therefore flattens the branches at both swapped positions; the walk up that this
needs replaces the root lookup each swap already did. Parent/child forests flatten
just those branches. Child-vector forests rebuild the whole crown once one of its
branches is deep, because inserting single branches into a long child vector would
take quadratic time.

Permuting the node store of the 9,800-half-edge grid with a fixed random
permutation takes these instruction counts, in thousands; "64 pending" first merges
64 nodes into a hub and discards history after the permutation:

#table(
  columns: 5,
  [Workload], [Store], [Before], [Flatten all roots], [Flatten merged branches],
  [No history], [Parent/child], [20,674], [20,674], [20,850],
  [No history], [Child-vector], [2,325], [2,325], [2,579],
  [64 pending, permutation], [Parent/child], [20,706], [20,706], [20,836],
  [64 pending, permutation], [Child-vector], [2,351], [2,351], [2,812],
  [64 pending, cleanup], [Parent/child], [111], [192], [53],
  [64 pending, cleanup], [Child-vector], [119], [220], [58],
)

These calls take more than 100 microseconds. Measured as above, wall-clock time
relative to "Before" is 0.98 (0.98--0.99) for parent/child and 1.10
(1.08--1.11) for child-vector without history, and 0.96 (0.96--0.97) and 1.14
(1.14--1.19) with 64 merges pending, including the cleanup. Flattening all roots
measures 1.00 to 1.02 in these cases, but can leave merged trees deep after a swap.

Child-vector swaps scan the parent's child vector, so they cost time linear in the
crown size once crowns are flat. After 2,400 merges into the hub and a cleanup,
the same permutation takes 1,568 million instructions with
flattening, against 50 million before; with the merges still pending, swap
flattening reaches the same cost (1,612 million). Parent/child forests take 20 to
24 million in all these cases.

== Remaining costs

Graphica remains about 1.8 times faster at construction and 2 times faster at DFS
on the 1,000-vertex sparse fixture. It is about 4.3 times faster at DFS on the
128-vertex, high-valence quotient. These measurements do not establish parity.

Retained requested heap storage is unchanged: 376,384 bytes for parent/child,
384,384 for child-vector, and 138,496 for Graphica on the sparse 1,000-vertex graph.
This includes vector capacity but excludes allocator metadata, stack usage and
process RSS. DFS uses nine allocations for parent/child and Graphica, versus 3,009
for child-vector; parent/child adjacency scanning allocates nothing.

Appending 32 edges immediately after a clone still costs about 102 microseconds
for parent/child versus 2 microseconds for Graphica at 1,000 vertices. Cloning and
final destruction are excluded, but growing full buffers is included. A separate
run first grows those buffers outside the timer; parent/child then takes roughly
5.5 microseconds. This demonstrates that growth dominates the cold case, while a
local insertion gap also remains. Allocation-sensitive append timings vary, and
no append speedup is claimed.

These are native topology measurements with integer payloads and no half-edge
payload. They do not measure Python overhead, generation, canonicalization or
rendering, and do not change the Community wheel.

== Validation and reproduction

The Linnet suite passes 264 tests. New coverage enumerates 720 small forests,
including shuffled IDs and nonroot subtree starts, and flattens the branch of
every node as well as all branches at once against the child-vector oracle,
including the single-branch paths that swaps use; it
also checks flattening, loops, parallel edges, dangling edges, isolated nodes,
payload preservation, subsequent mutation and extraction. Repeated identification
batches, with absorbed trees nested three deep and a hub moved below a fresh root,
check that every crown is flat after each cleanup, crown order is preserved and
untouched crowns keep their exact tree, including one deepened by rerooting. For every
absorbed node of an eight-node graph and every half-edge, swapping the absorbed
record's label away before cleanup still leaves no half-edge two hops below its
root. New paired-edge construction is compared
against connecting identity edges for every orientation. Spenso's 40 network
library tests pass.

The 58 benchmark cases pass a smoke run. Clippy passes for Linnet's library,
tests and benchmarks; existing vendored Clarabel warnings remain.

```sh
cargo nextest run -p linnet --no-default-features --features nodestore-vec
cargo test -p linnet --no-default-features --features nodestore-vec --bench nodestores -- --test
cargo clippy -p linnet --no-default-features --features nodestore-vec --lib --tests --benches
cargo bench -p linnet --no-default-features --features nodestore-vec --bench nodestores
cargo test -p spenso --lib network
```

The batched-cleanup instruction counts use `valgrind --tool=callgrind
--toggle-collect` on non-inlined merge, cleanup and permutation functions that
mirror the `hub_merge` and `permute` benchmarks, built against copies of the three
revisions.

The benchmarks instantiate both forest types explicitly. The parent/child-default
library also passes `cargo check --lib --no-default-features --features
nodestore-forest-parent-child`. The separate test configuration that combines a
forest default with archive support remains blocked by existing Vec-specific
field accesses in `parser/archive.rs`; the regular test configuration above still
exercises both forest implementations.

The comparison harness, baseline and optimized executables, raw CSV files,
compiler logs and revision/hash record are saved under
`/shared/codex/graph-scale/nodestore-bench/parent-child-performance/`; the
batched-cleanup harness, executable and raw counts are in its `batched-cleanup`
subdirectory, and the swap measurements in `batched-cleanup/swaps`. CI remains disabled and this change has not been pushed.
