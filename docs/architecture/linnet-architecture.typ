= Linnet implementation architecture

#quote(block: true)[
#strong[Status:] Current implementation architecture, audited against the Linnet source on
2026-10-13.

This note covers the `linnet` Rust crate. `clinnet` is a separate command-line renderer and
`linnest` owns Typst layout and drawing; neither is part of Linnet's graph core.
]

== Boundaries

Linnet is an in-process, generic graph library. Its main boundary is
`HedgeGraph<E, V, H, N>` in `crates/linnet/src/half_edge.rs`: callers choose edge data `E`, node
data `V`, per-half-edge data `H`, and a `Topology` implementation `N`. The crate is divided
into four layers:

- `half_edge`: topology, typed indexes, graph mutation, subgraph views, traversal, cuts, cycles,
  canonicalization, constrained graph generation, and directed-acyclic-graph algorithms;
- `parser`: DOT parsing and serialization, including graph-, node-, edge-, half-edge-, and port
  attributes;
- `permutation`, `union_find`, and `tree`: reusable index-reordering and connectivity machinery
  used by graph operations;
- feature-gated adapters: `half_edge::layout` under `drawing`, Symbolica interop under
  `symbolica`, and binary archive formats under `bincode` or `rkyv`.

Linnet does not invoke Graphviz, a browser, Typst, or another renderer. It reads and writes DOT;
rendering and document layout belong to downstream tools. The top-level `drawing` module is not
a supported graph-rendering pipeline; the implemented layout code is the feature-gated
`half_edge::layout` module.

== Core representation

`HedgeGraph` keeps three stores whose indexes must remain synchronized:

- `hedge_data: HedgeVec<H>` has exactly one payload per half-edge;
- `edge_store: SmartEdgeVec<E>` owns edge payloads and an `Involution<EdgeIndex>`;
- `node_store: NodeStore<N, V>` pairs the topology `N`, which maps each node to its incident
  half-edges, with one payload per node slot.

`Hedge`, `EdgeIndex`, and `NodeIndex` are distinct typed wrappers around compact indexes. They
are locations inside one graph, not durable identities. Operations that compact, extract, or
permute storage can change them.

The involution is the topological kernel. For every half-edge `h`, `inv(inv(h)) == h`. A paired
edge has one `Source` mapping with its `EdgeIndex` and one `Sink` mapping back to that source. A
dangling edge is an `Identity` fixed point and records its underlying `Flow`. `Orientation` is a
separate, superficial direction used by callers and serialization; it must not be confused with
the source/sink flow that maintains the pairing.

Owned, shared, and mutable involution iterators delegate traversal and size hints
to the underlying half-edge vector. They implement `ExactSizeIterator` and
`FusedIterator`, so their remaining length is exact and exhaustion is permanent.
These iterator guarantees preserve the store's half-edge order.

Subgraphs are views over half-edge sets. `SubGraphLike` combines subset membership with
edge-pair inclusion semantics, so an edge crossing a subset boundary becomes a split edge rather
than silently disappearing. Concrete forms such as `SuBitGraph`, `InternalSubGraph`,
`HedgeNode`, `Cycle`, and `OrientedCut` add the invariants needed by particular algorithms.

The `Topology` trait is the node contract: crowns, identification, compaction, partitions and
permutations, with numbering fixed by the shared tail-fill kernel so implementations agree
exactly. `NodeForest` is the topology in use: each node's crown is a ring of its half-edges that
ascends from the smallest one in every state, so crown order never depends on history. A
test-only reference topology states the same contract naively, numbering by the tail-fill rule
written out rather than by the kernel, and a differential test runs random edits on both.
`NodeStore` owns everything that does not depend on the topology: payloads, mapping, iteration,
building and the contract check.

Identification never renumbers: absorbed nodes keep their slot and payload until `compact`
drops them, so earlier `NodeIndex` values stay valid until then. `node_bound()` counts every
slot and is the size for dense side tables; `n_live()` counts the nodes that were not absorbed,
and node iteration yields only those. Every operation that renumbers reports a typed `Remap`:
`compact` maps absorbed nodes to their survivor, `delete_hedges`, `extract`, `extract_nodes` and
`contract_subgraph` report where the remaining half-edges and nodes went, and `permute_hedges`
and `permute_nodes` take bijections. Callers follow indices through these remaps instead of
mirroring the kernel. Partitions and permutations require compacted identifications, and
a `NodeForest` on its own can open nested checkpoints that roll back node additions,
attachments, anchors and identifications; its renumberings refuse while one is open. A
`NodeStore` lends its topology only by shared reference, so it never holds a checkpoint and its
`compact` cannot fail.

== Construction and data flow

Programmatic construction follows one path:

```text
HedgeGraphBuilder
  -> add nodes
  -> add paired or dangling edges to an Involution<E>
  -> build the NodeStore topology from each half-edge's node
  -> SmartEdgeVec::new separates edge payloads from the involution
  -> HedgeGraph<E, V, H, N>
```

Mutation methods update all three stores. Adding an edge extends half-edge data, edge storage,
and node incidence; `add_node` appends an empty-crown node in place. Neither renumbers existing
nodes, edges or half-edges. Extraction and deletion partition the same half-edge set in each store.
Complete-node extraction also moves explicitly selected empty-crown nodes, which have no
half-edge representation in that partition.
Joining appends the stores; sewing converts matched identity half-edges into pairs.
`HedgeGraph::check` verifies involutivity and that every `EdgeIndex` points to the same
`HedgePair` recorded by `SmartEdgeVec`. `NodeStore::check_nodes` checks the topology contract and
one payload per node slot.

Algorithms consume the graph together with an explicit subgraph when boundary behavior matters.
Traversal and cycle/cut routines operate on those views. Directed algorithms take a
`DirectionBasis`: `Underlying` follows the involution's source/sink roles, while `Superficial`
applies each edge's `Orientation` and omits superficially undirected edges. Dangling and
subgraph-split edges do not form complete directed arcs. Full topological sort canonicalizes
identification-history aliases with shared incidence, counts incoming arcs in the selected basis,
and rejects cycles; transitive closure and reduction first require that DAG check, then add or
delete complete half-edge pairs rather than one side of an edge.

Exhaustive circuit enumeration requires an explicit cycle-rank bound and rejects powersets that
the platform cannot represent before enumeration. Tadpole classification requires a non-empty
set of distinct live terminals, each with an incident half-edge. A node absorbed by a pending
identification is rejected; callers supply its live representative. It compares components against
every terminal-incident edge and returns deduplicated structural results. These guards report
malformed input as errors.

Connectivity counts use one union-find pass over the nodes touched by a subgraph's half-edges.
Dangling and split edges touch a node without joining components; internal edges determine the
component count and independent cycle rank. Spanning-forest enumeration lives in
`half_edge::spanning_forests`: it borrows the involution and rolls branch identifications back
on a private checkpointed topology.

`all_cuts` enumerates whole-node partitions in which every component on each side contains a
terminal of that side. Empty or overlapping terminal regions and live isolated nodes produce
no partitions. Dangling half-edges stay with their node; disconnected terminal-rooted regions
can have an empty boundary. Results are sorted by the source-side mask.

`half_edge::canon` computes canonical forms and automorphisms from caller-supplied node, edge,
and half-edge keys, respecting physical orientation. Payloads omitted from those keys do not
distinguish graphs; canonical forms are comparable within one canonizer version.
`half_edge::generate` separates a checked `Plan` built from `Theory`, `Legs`, and `Query` from
execution controls in `Run`. It generates topologies, then decorates them with fields and vertex
rules. Pruning must reject every refinement of a rejected state; final selection may be arbitrary.
Accepted diagrams retain compact data and materialize graphs and canonical forms on demand.
Parallel execution preserves the serial search order.

`Query.classes(Classes::Auto)` enables a conservative class search for directed fields that
occur only in single diagonal lines. Undirected fields and participants in mixed or multi-line
rules remain distinct. Exact automorphism orbits of the coloured rule-incidence graph define
the eligible classes; active budgets and caller distinctions can only split them. External
label classes retain their original particle pins.

A complete coarse assignment is refined by whole directed lines, with its stabilizer used
for the second symmetry quotient. Fine automaton checks retain original rule indices,
including equal-signature rules, and exact weights and automorphism orders are computed for
the resulting particle assignment. Hereditary dressing hooks also inspect refinement states.
`Classes::None` bypasses the quotient and serves as the direct-particle comparison route.

DOT interchange follows a separate adapter path:

```text
DOT text/file
  -> parser::syntax: lexing and recursive descent, flattened in one pass
  -> FlatGraph: attribute statements, nodes keyed by decoded identifier, single edges
  -> DotVertexData / DotEdgeData / DotHedgeData
  -> HedgeGraphBuilder
  -> explicit id permutations, when ids were supplied
  -> DotGraph
```

Invisible endpoint nodes represent dangling edges. `dir` controls superficial orientation, while
the endpoint side determines underlying flow. Serialization maps caller data into the DOT data
types and writes through a caller-provided `fmt::Write` or `io::Write`; DOT parsing separately
offers string- and file-based entry points. `DotGraph` requires exactly one graph and rejects
trailing input; `DotGraphSet` accepts multiple graphs. Invalid DOT is a `DotSyntaxError` with
its line and column; the in-tree reader limits subgraph nesting to 256 levels, because the
Typst plugin parses user-supplied DOT. Malformed direction, endpoint, and explicit-ID data
produce `DotParseError` without publishing a partially constructed graph. Later graph-attribute
statements take precedence, and parallel edges with equal endpoints retain their statement order.
Quoted values preserve literal backslashes and Graphviz escapes. The shared `DotIndent` formatter
keeps newlines and indentation inside quoted values intact in whole-graph and subgraph output.

== Feature and persistence boundaries

The default `serde` feature derives serialization for the generic topology types when their
payloads support it. `bincode` adds its encode/decode derives. `rkyv` adds archive validation
metadata and byte-oriented APIs for `DotGraph`; its unsafe raw-root accessor still requires the
caller to supply bytes for exactly the archived type. Callers remain responsible for storing and
versioning those bytes. DOT is the human-readable interchange format. Linnet has no state
directory, database, or background persistence service.

The `drawing` feature enables `cgmath` and the layout modules. The `symbolica` feature enables
the half-edge Symbolica adapter. These features extend data conversion or analysis; they do not
change the involution and node-incidence contracts.

== Maintained invariants

- An internal edge is always represented by two mutually paired half-edges; an external edge is
  an identity half-edge.
- `hedge_data`, the involution, edge records, and node incidence describe the same half-edge
  index space.
- Edge payload and `HedgePair` records agree with the `EdgeIndex` stored in the involution.
- A subgraph belongs to the graph whose half-edge extent it was created from; reusing its bitset
  with an unrelated or structurally modified graph is invalid.
- Mutations that reorder compact storage also reorder every dependent store. Callers must retain
  returned graphs, offsets, or permutations instead of assuming old numeric indexes survive.
- Algorithms that require a directed acyclic graph report `TopoError` or `TransitiveError`
  rather than producing a partial ordering for a cyclic graph.

== Verification and related documentation

Unit tests live beside the implementation and exercise involution and node-store consistency,
exact vector/forest conversion, extraction/join/sew operations, subgraphs and cuts, DOT round
trips, permutations, union-find, topological and transitive algorithms, and optional archived
views. Snapshot tests pin graph, cycle, ordering, and DOT representations. The default boundary
is exercised with `cargo test -p linnet`; the `rkyv` archive path requires that feature
explicitly.

For supported usage, start with the
#link("../../../products/linnet/latest/tutorial/")[first-graph tutorial], continue with the
#link("../../../products/linnet/latest/guides/algorithms/")[algorithms guide], and use the
#link("../../../products/linnet/latest/reference/rust/linnet/")[native Linnet Rustdoc]
for exact public signatures. The
#link("../../../products/linnet/latest/guides/clinnet/")[Clinnet guide] and
#link("../../../products/linnet/latest/guides/linnest/")[Linnest guide] document the downstream
rendering boundaries that this crate intentionally does not own.

== Dense edge joins used by tensor contraction

Joining graph edges updates moved payloads through their owning half-edges. Node identification and
edge deletion preserve slot order, flow, graph extent, and complete node membership. These operations
remain on the existing half-edge owner; no alternative graph representation is introduced.
