#import "../../shared.typ": callout, boundary

#let algorithms = [
= Subgraphs, traversal, cuts, and drawing

Linnet algorithms operate on a graph together with an explicit subgraph view. This lets one
algorithm work on the full graph, a selected region, or the boundary created by cutting paired
half-edges. Choose the subgraph representation by operation: bit-backed sets are efficient for
repeated set algebra, while borrowed or mapped views avoid copying when identity must remain
attached to the owning graph.

== Traversal and connectivity

A traversal needs a starting half-edge or node, the active subgraph, and a policy describing
which adjacent half-edges may be crossed. Results retain typed indexes into the original graph.
Breadth-first traversal is useful for distance layers and connectivity; depth-first traversal
is useful for trees, back edges, and decompositions. A traversal tree is a derived view, not a
new graph, so mutations must be coordinated with the indexes it contains.

Connectivity and traversal results depend on the chosen subgraph. A traversal tree retains its
parent relation, so callers can inspect immediate children, strict ancestors, and the fundamental
cycle closed by a non-tree internal edge without reconstructing an adjacency model. Boundary
queries distinguish half-edges already selected at a split from every boundary half-edge incident
to a touched node; dangling half-edges are available as a separate structural selection.

Linnet has two independent notions of direction. `DirectionBasis::Underlying` follows the
source/sink roles stored by the involution. This is the default in Python and is the basis used by
GammaLoop momentum routing. `DirectionBasis::Superficial` additionally applies an edge's
`Default`, `Reversed`, or `Undirected` orientation; a superficially undirected edge contributes no
directed arc. Reachability, topological order, transitive closure, and transitive reduction take
this choice explicitly. Neither mode rewrites stored flow or orientation.

Directed algorithms ignore dangling half-edges and edges split by the active subgraph boundary,
because neither represents a complete arc inside that selection. Topological ties use node-index
order for the initial ready set and deterministic edge-discovery order thereafter. Linnet does not
currently expose a dedicated articulation-point or biconnected-component API; derive those from
traversal data only when your application owns that algorithm and its tests.

== Cycles and cuts

`cycle_basis` returns an independent cycle vector together with the spanning forest used to
construct it. Use `cyclotomatic_number` when the numerical cycle rank itself is the result.
`all_cycles` expands a basis and can grow exponentially; prefer a basis when assigning loop
momenta. Native callers supply a maximum accepted cycle rank. Python callers supply
`max_results`, which bounds the complete non-empty basis-combination search space rather than
only the circuits that survive filtering. Both guards are checked before powerset construction.
Symmetric differences combine cycle bitsets without inventing new graph identities.

Cut enumeration separates required left and right node sets and returns the left region, cut
content, and right region. A cut edge is represented through its half-edges, preserving which
side owns each endpoint. Validate direction/flow after enumeration when a physical Cutkosky or
tensor-network interpretation depends on orientation.

#callout("Complexity is part of the API choice", [
  Enumerating every cycle or every admissible cut is inherently combinatorial. Filter the
  subgraph and endpoint sets first, and use basis or connectivity operations when enumeration
  is not required by the caller.
])

== Enumerate separating cuts

This complete example builds a square with a diagonal and asks for every cut separating opposite
nodes. It is the smallest graph that makes the distinction between a path and a cut family
visible without external data or drawing support.

// docs-example: run linnet-separating-cuts
```rust
use linnet::half_edge::{HedgeGraph, builder::HedgeGraphBuilder};

fn main() {
    let mut builder = HedgeGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());

    builder.add_edge(a, b, (), true);
    builder.add_edge(b, c, (), true);
    builder.add_edge(c, d, (), true);
    builder.add_edge(d, a, (), true);
    builder.add_edge(b, d, (), true);

    let graph: HedgeGraph<(), (), ()> = builder.build();
    graph.check().expect("the half-edge involution is valid");
    let cuts = graph.all_cuts_from_ids(&[a], &[c]);

    assert_eq!(cuts.len(), 4);
    println!("separating cuts: {}", cuts.len());
}
```

The expected invariant is four distinct left/cut/right partitions separating `a` from `c`.
Changing topology can change that count, while changing only node or edge payloads does not.
The native #link("reference/rust/linnet/half_edge/struct.HedgeGraph.html")[`HedgeGraph`
Rustdoc] and
#link("reference/rust/linnet/half_edge/builder/struct.HedgeGraphBuilder.html")[builder Rustdoc]
give the exact compiled type boundaries for this revision.

The Python binding delegates to that same native enumeration and returns named graph-bound
results rather than nested tuples:

```python
partitions = graph.all_cuts(["a"], ["c"])
for partition in partitions:
    source_side = partition.source_side
    boundary = partition.boundary
    cut_edges = partition.edges
    target_side = partition.target_side
```

The source and target sides are structural `Subgraph` views. `boundary.left` and
`boundary.right` select the two oriented half-edges of every cut edge. The result is eagerly
materialized and sorted by half-edge membership at the Python boundary; do not use all-cut
enumeration where one witness or a minimum-cardinality cut is the actual question.

`all_bonds(min_size=..., max_size=...)` is the related native enumeration of
inclusion-minimal cutsets whose two induced sides are connected. A bond size is its boundary
half-edge count; it is not a max-flow result. `find_bond` short-circuits after one witness in the
same inclusive size range. Both full enumerations can be combinatorial.

`tadpoles(externals)` identifies a validated terminal set structurally and returns the
native tadpole components in stable order. Terminals must resolve to distinct structural nodes;
empty or duplicate sets, out-of-range nodes, and zero-crown terminals are errors instead of panic
paths.

#callout("Interpret a cut mismatch structurally", [
  A failed `graph.check()` indicates a storage or involution invariant and must be resolved before
  enumeration. Overlapping, empty, isolated, or closed terminal groups are rejected. For validated
  groups, zero cuts means the native enumerator found no admissible separating partition, commonly
  because the groups occupy disconnected components. An unexpected nonzero count means the
  constructed edges or endpoint sets do not describe the graph you intended; inspect those before
  adding post-enumeration filters.
])

== Canonical forms and symmetry

`HedgeGraph::canonize_by` takes one key function each for node, edge, and half-edge payloads
(`canon::ignore` disregards one) and returns a `Canonization`: canonical half-edge, node, and edge
labels, a `CanonicalForm` that is equal exactly for isomorphic graphs, automorphism generators as
half-edge permutations with their node permutations, orbits, and the exact automorphism count.
`canonical_form_by` (and `canonical_form`, with every payload as its key) computes the canonical
form alone, without labels, generators or count. A canonical form holds words that encode the
canonically relabelled structure, with colours that rank into the sorted distinct node, edge and
half-edge keys, and those key tables; `into_parts` moves them out. Canonical forms are comparable
only between graphs canonized by the same version of Linnet.
Symmetry is half-edge symmetry: an automorphism permutes half-edges compatibly with the
involution, so `m` interchangeable parallel edges or dangling legs contribute `m!` and an
undirected self-loop with indistinguishable ends contributes `2`. Undirected edges ignore their
stored source/sink roles, while directed edges and legs keep their physical direction, including
`Reversed` orientations. `is_isomorphic_by` and `isomorphism_by` compare two graphs through their
canonical forms.

== Feynman-graph generation

`generate::GraphGenerator` enumerates the connected graphs built from vertex signatures and
external legs, one per isomorphism class. A signature is a multiset of typed stubs: a field with a
direction relative to the vertex (undirected for self-conjugate fields, incoming or outgoing for
fermions, ghosts and other charged fields). Two stubs pair when they carry the same field and are
both undirected, or one is outgoing and the other incoming. Legs are `(label, stub)` pairs that
become dangling half-edges keyed by their label; repeating a label makes those legs
interchangeable. Settings cover a loop range, a vertex bound (legs are not vertices), a bridge
limit (`0` for one-particle-irreducible graphs), self-loops and zero-flow bridges (a bridge with no
leg on one side), budgets (sums of integer signature weights such as coupling and counterterm
orders, and counts of vertices, propagators, chords and closed lines of some fields, each within a
range), a filter that prunes the search, and progress and abort callbacks. The filter
sees read-only views of the partial untyped topologies, the partial field assignments and the
complete graphs; it must be hereditary (reject every completion of a state it rejects), so the
output is exactly the complete graphs it accepts, and a validation mode searches without pruning
and fails when a rejected state completes to an accepted graph. Only paired edges count as bridges: a one-leg process keeps its graphs without
internal bridges, where Graphica counts the leg as a zero-flow bridge and drops them.

Each output carries the `HedgeGraph` (fields on edges, signature indices on nodes), its loop
order, the exact automorphism count and the canonizer's canonical form under the keys
`(signature, field, Option<label>)`, so a graph from elsewhere can be matched against generated
ones by canonizing it with the same keys. Signatures with identical stubs are distinct vertex rules.

Generation runs in two stages per loop order. Untyped topologies are built by McKay's canonical
augmentation, one vertex at a time: rule 1 keeps one multiplicity vector per orbit of the parent's
node automorphism group, and rule 2 accepts a child only if the new vertex is in the orbit of its
canonical deletion (the deletable vertex with the largest invariant, ties broken by a node
canonization). Every prune (completability, permanent bridges, zero-flow bridges) is a necessary
condition for completion, so each output appears exactly once; the filter prunes stage 1 and the
partial field assignments of stage 2 the same way. Once the counts force every remaining vertex to
be a pendant on an open stub, the pendants are hung all at once, one assignment per orbit. Three or
more legs with distinct labels and a common stub are placed as anonymous legs while the topology
is shaped and labelled afterwards, each label on one vertex per orbit of the group fixing the
previous labels, so many distinct legs cost few symmetric topologies instead of a search through
labelled dead ends. Twin vertices (equal corollas and multiplicities to every other vertex) are
exchanged by an automorphism, so rule-2 ties among twins and groups whose colour classes are twin
classes need no node canonization. Each complete topology is then
decorated with fields, orientations and vertex rules edge by edge, keeping one decoration per orbit
of the topology's automorphism group. A partial vertex is a state of an automaton over the
sub-multisets of the signatures, so an edge tries only the fields both of its ends can still take,
and forward checking drops a partial decoration before the filter sees it once an edge at the ends
of the last decorated one has no field left; a topology with an edge that no field fits at its
ends gets no decoration. Typed canonical augmentation (signatures as the corollas) was measured as the
alternative and visits about four times as many partial objects on QED and QCD processes.

Budgets prune the decorations in the same way. Every automaton state carries the least and the
most weight of the signatures containing it, so the sums of these bounds over the vertices bound
each weight of a partial decoration, and a decoration step changes them at its two ends only;
propagator and chord counts are bounded by the edges counted and the edges still undecorated, and
closed lines by the lines closed so far (a union-find over the decorated edges of their fields),
unless a signature lets lines branch. A partial decoration whose bounds miss a range is dropped
before the filter sees it, and the vertex rules are checked exactly; the weights are invariant
under the stabilizer of a labelling, so keeping the orbit-minimal rule tuples that meet them is
exact. Before any search, a loop order is skipped when no multiset of signatures has the degrees
its Euler relation needs and weights within the ranges, and a request is empty when its legs
carry a charge or a parity that no combination of vertices balances (the leg vector lies outside
the span of the signatures' flows, modulo a large prime, or of their stub parities, modulo two).
The emptiness test follows the reachable sums of degrees and weights rather than the vertex
partitions themselves: on SMEFTsim with an NP budget the partitions number tens of thousands
(2->3 at three loops: 46,697) while the sums stay a few thousand (2,533).

#callout("Streaming, deterministic order and parallelism", [
  Both stages are resumable depth-first searches with explicit stacks, so the iterator streams:
  each `next` searches only up to the next graph, and dropping it ends the search. Graphs come in
  the depth-first order of the search. With the `rayon` feature, `parallel(true)` runs a window of
  independent subtrees on the current thread pool, each at most a fixed number of graphs ahead of
  the output, splits the unsearched rest of the oldest subtrees while fewer subtrees than threads
  can advance, and returns their results in that same order, so the output does not depend on the
  number of threads. The first step runs only the oldest subtree, up to its first graph, and
  further steps buffer one graph per subtree until the first is returned, so the first graph does
  not wait for a step of the whole window.
])

== Mutation boundaries

Python insertion extends a graph with reusable `Node` and `Edge` descriptions or direct
`add_node` / `add_edge` arguments. Existing live element views and their drawing metadata remain
valid after appending nodes or edges. Structural selections and algorithm results are tied to the
complete topology revision and become stale even after insertion.
`split_edge` separates one paired edge into two dangling edges, while `connect` joins two exact
dangling half-edges. Extraction copies a selected region, excision separates a graph along its
boundary, sewing joins compatible dangling half-edges, and contraction identifies structure.
These transformations return replacement views or indices where needed and conservatively
invalidate older live element views. Reacquire elements by name or index after changing existing
structure. Detached Python payloads remain usable independently of graph topology.

Transitive reduction is an owning transformation over a directed acyclic graph. The Python method
returns a new graph, preserves arbitrary payload identity and copied drawing configuration, and
leaves views into the source graph live. Transitive closure remains a Rust-only operation until the
Python API has an explicit factory for every new edge and half-edge value it synthesizes.

== DOT and drawing

DOT parsing and emission are interchange/debugging boundaries. Attributes can carry graph,
node, edge, and half-edge data, but the typed Rust invariants are reconstructed by the parser;
handle parser errors before running algorithms. The `drawing` feature enables layout; rendering
belongs to a DOT consumer, Clinnet, or Linnest. Layout coordinates are presentation data and may
vary without changing topology, subgraph membership, or edge flow.

#boundary("Tensor contraction uses Spenso", [
  Linnet determines connectivity and graph transformations. Tensor compatibility, contraction
  cost, and execution belong to Spenso even when a Spenso network delegates its graph
  algorithms to Linnet.
])
]
