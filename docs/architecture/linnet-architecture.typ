= Linnet implementation architecture

#quote(block: true)[
#strong[Status:] Current implementation architecture, audited against the Linnet source on
2026-08-18.

This note covers the `linnet` Rust crate. `clinnet` is a separate command-line renderer and
`linnest` owns Typst layout and drawing; neither is part of Linnet's graph core.
]

== Boundaries

Linnet is an in-process, generic graph library. Its main boundary is
`HedgeGraph<E, V, H, S>` in `crates/linnet/src/half_edge.rs`: callers choose edge data `E`, node
data `V`, per-half-edge data `H`, and a `NodeStorage` implementation `S`. The crate is divided
into four layers:

- `half_edge`: topology, typed indexes, graph mutation, subgraph views, traversal, cuts, cycles,
  and directed-acyclic-graph algorithms;
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
- `node_store: S` owns node payloads and maps each node to its incident half-edges.

`Hedge`, `EdgeIndex`, and `NodeIndex` are distinct typed wrappers around compact indexes. They
are locations inside one graph, not durable identities. Operations that compact, extract, or
permute storage can change them.

The involution is the topological kernel. For every half-edge `h`, `inv(inv(h)) == h`. A paired
edge has one `Source` mapping with its `EdgeIndex` and one `Sink` mapping back to that source. A
dangling edge is an `Identity` fixed point and records its underlying `Flow`. `Orientation` is a
separate, superficial direction used by callers and serialization; it must not be confused with
the source/sink flow that maintains the pairing.

Subgraphs are views over half-edge sets. `SubGraphLike` combines subset membership with
edge-pair inclusion semantics, so an edge crossing a subset boundary becomes a split edge rather
than silently disappearing. Concrete forms such as `SuBitGraph`, `InternalSubGraph`,
`HedgeNode`, `Cycle`, and `OrientedCut` add the invariants needed by particular algorithms.

The default build resolves `DefaultNodeStore` to `NodeStorageVec`. The two forest-backed node
stores are mutually exclusive features. The distinction matters: vector-backed node
identification can invalidate previous `NodeIndex` values, whereas the forest-backed
implementations preserve identification history until it is explicitly forgotten.

== Construction and data flow

Programmatic construction follows one path:

```text
HedgeGraphBuilder
  -> add nodes
  -> add paired or dangling edges to an Involution<E>
  -> build NodeStorage from each node's incident half-edges
  -> SmartEdgeVec::new separates edge payloads from the involution
  -> HedgeGraph<E, V, H, S>
```

Mutation methods update all three stores. Adding an edge extends half-edge data, edge storage,
and node incidence. Extraction and deletion partition the same half-edge set in each store.
Joining extends the stores and then revalidates node incidence; sewing converts matched identity
half-edges into pairs. `HedgeGraph::check` verifies involutivity and that every `EdgeIndex` points
to the same `HedgePair` recorded by `SmartEdgeVec`. Node-store operations provide the additional
partition and incidence checks.

Algorithms consume the graph together with an explicit subgraph when boundary behavior matters.
Traversal and cycle/cut routines operate on those views. Topological sort counts incoming
half-edges by `Flow` and rejects cycles; transitive closure and reduction first require that DAG
check, then add or delete complete half-edge pairs rather than one side of an edge.

DOT interchange follows a separate adapter path:

```text
DOT text/file
  -> dot-parser AST
  -> SubGraphFreeGraph plus inherited attributes
  -> DotVertexData / DotEdgeData / DotHedgeData
  -> HedgeGraphBuilder
  -> explicit id permutations, when ids were supplied
  -> DotGraph
```

Invisible endpoint nodes represent dangling edges. `dir` controls superficial orientation, while
the endpoint side determines underlying flow. Serialization maps caller data into the DOT data
types and writes through a caller-provided `fmt::Write` or `io::Write`; DOT parsing separately
offers string- and file-based entry points.

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
extraction/join/sew operations, subgraphs and cuts, DOT round trips, permutations, union-find,
topological and transitive algorithms, and optional archived views. Snapshot tests pin graph,
cycle, ordering, and DOT representations. The default boundary is exercised with
`cargo test -p linnet`; the `rkyv` archive path requires that feature explicitly.

For supported usage, start with the
#link("../../../products/linnet/latest/tutorial/")[first-graph tutorial], continue with the
#link("../../../products/linnet/latest/guides/algorithms/")[algorithms guide], and use the
#link("../../../products/linnet/latest/reference/rust/supported/linnet/")[curated Rust reference]
for exact public signatures. The
#link("../../../products/linnet/latest/guides/clinnet/")[Clinnet guide] and
#link("../../../products/linnet/latest/guides/linnest/")[Linnest guide] document the downstream
rendering boundaries that this crate intentionally does not own.
