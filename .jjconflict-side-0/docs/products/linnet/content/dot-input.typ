#import "../../shared.typ": callout, product-link

#let dot-input = [
= DOT input and graph indices

Linnet reads DOT into a graph of nodes, edges, and half-edges. This guide explains how DOT
endpoints become graph records and how their indices are assigned. The same rules apply when
DOT enters through Python, Linnest, or GammaLoop; a consumer can subsequently transform or
reorder the parsed graph.

== Edges and half-edges

An internal edge has two half-edges, one attached to each endpoint. A dangling edge has one
half-edge attached to its real node. Represent its missing endpoint with an invisible DOT node:

// docs-example: syntax
```text
digraph example {
  ext [style=invis];
  ext -> a;
  a -> b;
  b -> ext;
}
```

This graph has two nodes, three edges, and four half-edges. The invisible `ext` is not a graph
node. `ext -> a` contributes one incoming half-edge at `a`; `b -> ext` contributes one outgoing
half-edge at `b`; `a -> b` contributes one half-edge at each node. A self-loop also contributes
two half-edges, both attached to the same node. Connecting two invisible endpoints is invalid.

One invisible node can represent several dangling endpoints: its incident edges remain separate
edges, not a connection through `ext`. Keep shared invisible nodes otherwise attribute-free,
since their attributes are merged into each attached dangling edge.

== Indices, not arbitrary labels

An edge attribute `id=2` requests edge index 2. A numeric endpoint port such as `a:3` requests
half-edge index 3 at `a`. These are positions in separate, densely indexed graph collections:

#table(
  columns: (auto, auto, 1fr),
  table.header([*Record*], [*DOT spelling*], [*Allowed indices*]),
  [Node], [`a [id=0]`], [From `0` through `V - 1`, where `V` is the number of visible nodes.],
  [Edge], [`a -> b [id=2]`], [From `0` through `E - 1`, where `E` is the total number of edges.],
  [Half-edge], [`a:3 -> b`], [From `0` through `H - 1`, where `H` is the total number of half-edges.],
)

With `I` internal edges and `D` dangling edges, `E = I + D` and `H = 2 I + D`. Count the complete
parsed graph, including dangling edges. An invisible endpoint contributes no extra half-edge;
put a dangling half-edge's numeric port on its real node.

Each explicit index must occur only once within its collection, across the entire graph.
`a:0 -> b:0` requests the same half-edge index twice and is rejected, even though the endpoints
are different nodes. Two edges with `id=0` also collide. In contrast, an edge index 0 and a
half-edge index 0 can coexist: they identify different kinds of records. The node-index
collection is independent of both.

Indices must be nonnegative and in range. On the three-edge, four-half-edge graph above,
`id=3` and a numeric port `a:4` are out of bounds. An explicit index does not create empty
records to extend the collection. Use ordinary attributes such as `name` or `label` for
application names or displayed text.

Node identifiers that parse as nonnegative integers, including quoted numeric identifiers,
also request that node index unless a node `id` attribute overrides it. Thus a one-node graph
containing only `7;` is out of bounds; `7 [id=0];` explicitly assigns its sole node to index 0.
Alphabetic identifiers such as `a` leave the index unspecified unless they have an `id`.

Write numeric half-edge ports without quotes: `a:3` or `a:3:n`. Named ports and compass points
such as `a:output:n` and `a:n` do not reserve half-edge indices. Currently, quoted numeric ports
such as `a:"3"` also remain port labels rather than fixing an index.

#callout("Numeric ports are global", [
  `a:3` means half-edge 3 of the whole graph, not the fourth connection at node `a`.
  It neither selects an edge index nor encodes a screen coordinate. The two halves of an
  internal edge can have unrelated indices; do not infer them from its edge index.
])

== Omission and inferred order

You can omit every index, supply every index, or fix only the records referenced elsewhere.
The parser reserves the explicitly requested positions, then fills the remaining positions
in ascending order. It handles all three collections independently, so specifying an edge's
`id` does not also specify either of its half-edge indices.

The order used to fill omitted indices comes from parsing the DOT endpoints, not simply from
the order of lines in the file:

- Nodes are traversed in lexicographic DOT identifier order.
- Edges are sorted by their source endpoint and then their target endpoint. An endpoint compares
  its DOT identifier spelling first, then any port or compass information. Ports therefore
  participate in sorting, even when they are not numeric indices. Equal endpoint pairs retain
  their order in the flattened DOT input, so statement order still matters for parallel edges
  with identical endpoints.
- Half-edges follow that edge traversal: an internal edge contributes its source half-edge
  before its sink half-edge; a dangling edge contributes only its real half-edge.

Each collection skips its explicit reservations when assigning omitted indices. The sorting
uses DOT spellings, not numeric magnitude, and happens before the explicit index permutations.

For example, this input deliberately puts `b -> ext` before `a -> b`:

// docs-example: syntax
```text
digraph indices {
  ext [style=invis];
  b -> ext;
  a:3 -> b [id=2];
  ext -> a:0;
}
```

The endpoint traversal is `a -> b`, `b -> ext`, then `ext -> a`. Edge 2 and half-edges 0 and 3
are already reserved. The resulting graph is:

#table(
  columns: (auto, auto, auto, auto),
  table.header([*Edge index*], [*Connection*], [*Source half-edge*], [*Sink half-edge*]),
  [`0`], [`b -> ext`], [`2` at `b`], [None],
  [`1`], [`ext -> a`], [None], [`0` at `a`],
  [`2`], [`a -> b`], [`3` at `a`], [`1` at `b`],
)

Changing an endpoint name, adding a port, or inserting an edge can change inferred indices.
Supply the indices that an external reference depends on, and inspect or export the parsed
graph before using numeric references. Graph edits after parsing can also change indices;
keep any index mappings returned by the operation.

== Chains, subgraphs, and repeated defaults

The parser expands edge chains and subgraph endpoints before counting and assigning indices.
For example, `a -> b -> c` creates two edges. Its trailing attributes apply to both, so
`a -> b -> c [id=0]` requests edge 0 twice and is rejected. Likewise, `a -> {b c} [id=0]`
creates a duplicate index. A numeric port in the middle of a chain is also repeated:
`a -> b:0 -> c` requests half-edge 0 on both incident edges at `b` and is rejected. Write
separate edge statements if these edges or half-edges need explicit indices, or omit the
indices and inspect the resulting order.

The same uniqueness rule applies when node or edge defaults supply `id` values. A subgraph does
not start a new index space: all visible nodes, expanded edges, and half-edges belong to the
containing graph. Parallel edges remain distinct records.

== Inspect and draw

Use the #link("quickstart/python/")[Python quickstart] to parse and inspect graph records,
or the #link("playground/")[live DOT playground] to edit a diagram in the browser.
For drawing an existing DOT collection, use the #link("guides/clinnet/")[Clinnet guide].

#product-link("gammaloop", page: "guides/dot-input/", label: "GammaLoop's DOT input guide")
explains which indices select external kinematics, identify numerator data, or appear in
momentum labels. Amplitude kinematics follow increasing dangling half-edge order. GammaLoop
uses the same order for initial drawing rows, with incoming and outgoing legs in separate
columns; its cross-section rows follow numeric `is_cut` order. Its example deliberately gives edges different indices from their
half-edges: changing drawing order does not renumber the edge-based $q$ labels. Cross-section
physics import can subsequently sew and renumber records, as described in that guide.
]
