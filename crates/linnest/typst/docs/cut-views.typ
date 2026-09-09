#set document(title: "Linnest directed cut views")
#set page(paper: "a4", margin: 20mm)
#set text(size: 10pt)

= Linnest directed cut views

Implemented design for `examples/xbox.typ`. Each view opens *one weighted
directed cut* of a shared master graph. `left` and `right` describe its two
drawing boundaries; repeated passages are hedge annotations, not separate named
cut objects. The general API tutorial is in `docs/manual.typ`; the momentum
algebra below is specific to this example.

== One graph, four openings

Ignoring layout-only compactification springs, the physical master is the
four-vertex, six-edge graph $K_4$:

#table(
  columns: 4,
  table.header([Edge], [Momentum], [Underlying source → sink], [Drawing orientation]),
  [D1], [$k-p_1$], [`c → d`], [default],
  [D2], [$k-p_2$], [`a → b`], [reversed],
  [D3], [$p_1$], [`c → a`], [reversed],
  [D4], [$p_2$], [`d → b`], [default],
  [D5], [$p_(12)-k$], [`a → d`], [default],
  [D6], [$k$], [`b → c`], [default],
)

`xbox` is already an opening, not the uncut master: its two D3 fragments are the
ends of one edge, as are its two D4 fragments. Conceptually identifying each
pair recovers the master; this is not a claim about data-preserving `graph.join`.

The following selectors refer to *original* half-edges. For example,
`D1.source` is the hedge at `c`, and `D6.sink` is the hedge at `c`. Use
`subgraph.select(g, source: (<D1>,), sink: (<D6>,))` to select that pair, rather
than relying on numeric hedge ordering. All windings are one except D6 in
`xbox-cut`, which has winding two on either or both selected halves.

#table(
  columns: (1.2fr, 2fr, 2fr, 1.5fr),
  table.header([View], [Left original hedges], [Right original hedges], [D1…D6 segment counts]),
  [`xbox`],
  [`D3.sink`, `D4.sink`],
  [`D3.source`, `D4.source`],
  [`1, 1, 2, 2, 1, 1`],
  [`xbox-opened`],
  [`D1.source`, `D6.sink`, `D4.sink`],
  [`D1.sink`, `D6.source`, `D4.source`],
  [`2, 1, 1, 2, 1, 2`],
  [`xbox-opened2`],
  [`D2.source`, `D3.sink`, `D6.sink`],
  [`D2.sink`, `D3.source`, `D6.source`],
  [`1, 2, 2, 1, 1, 2`],
  [`xbox-cut`],
  [`D1.source`, `D2.source`, `D6.sink` (winding 2)],
  [`D1.sink`, `D2.sink`, `D6.source` (winding 2)],
  [`2, 2, 1, 1, 1, 3`],
)

L and R mean the drawing's left and right boundaries, not structural source and
sink and not superficial fermion-arrow orientation. Count a seam passage
*R → L as positive*. An edge whose source hedge is on R and sink hedge is on L
therefore contributes its momentum positively; swapping the sides changes the
sign. Reversing a drawn fermion arrow does not change this convention.

Writing $q_e$ for the table's momentum and $epsilon_e$ for this sign, the weighted
cut carries $Q = sum_(e " cut") epsilon_e w_e q_e$. The four descriptions give
exactly the same total:

$ Q_("xbox") &= p_1 + p_2 \\
  Q_("xbox-opened") &= -(k-p_1) + k + p_2 = p_1 + p_2 \\
  Q_("xbox-opened2") &= -(k-p_2) + p_1 + k = p_1 + p_2 \\
  Q_("xbox-cut") &= -(k-p_1) - (k-p_2) + 2k = p_1 + p_2. $

The cuts need not be graph-theoretic separating partitions: cutting D3 alone
does not disconnect $K_4$. Terminal-cut or bond enumeration is not a prerequisite
for opening an explicitly selected edge set. The red Cutkosky line remains a
separate drawing annotation. Four blue shades and dash patterns identify the
original opening, the $p_1$ replacement, the $p_2$ replacement, and the combined
replacement. Each view has two matching vertical cuts through its left and
right free endpoints exactly, including the auxiliary endpoints of the
through-gluon. The endpoints on each side share their solved x coordinate.

On the first view, the $p_1$ Bezier crosses D1/D6 and the right dangling D4
($p_2$); the $p_2$ Bezier crosses D2/D6 and the left dangling D3 ($p_1$).
The combined pair crosses D6 twice in the same direction without crossing
external legs. Each alternative uses exactly the stroke of its resulting view.
All top/bottom endpoints share the red final-state cut's vertical extent and
arrive perpendicularly; the combined pair also meets the side boundaries
perpendicularly. These overlays never cross D5 and do not change topology or
layout.

== Winding two retains the middle segment

`graph.cut(g, left: left, right: right, boundary: patch)` is nonmutating.
Its selections must be disjoint involution partners on paired edges. A positive
integer `(winding: n)` annotation requests same-direction seam passages;
missing winding defaults to one, and values supplied on both sides must agree.
This is not a net winding that cancels back-and-forth crossings.

For D6 in `xbox-cut`, winding two generates `<D6.0>` at `b`, `<D6.1>` between
two auxiliary boundary nodes, and `<D6.2>` at `c`, ordered along underlying
`b → c` flow. Crossing 0 connects the right stub to the middle's left endpoint;
crossing 1 connects the middle's right endpoint to the left stub. Both seam
passages are R → L, although the middle segment runs L → R within the drawing.
A plain boolean union would forget the multiplicity and lose this segment.

The example annotates D6 with two `crossings` entries, each carrying a momentum
`group` and a `y` placement. Its local `cut-side` helper sets `winding` from the
entry count; `crossing` indexes the entries for placement and endpoint boxes.
These are view-specific annotations, not two named seam objects or an ordering
argument to `graph.cut`. Positive winding and underlying flow determine fragment
order; Linnest itself only interprets the `winding` field.

Native Linnet requires every hedge to belong to a node; an edge cannot have two
free endpoints without incident nodes, and an already dangling stub cannot be
split again. The implementation uses native splitting for the outer stubs and
appends disconnected middle edges on auxiliary nodes. Those nodes are
zero-sized and unpainted automatically, but remain queryable and positionable.

== Data, provenance, and placement contracts

- `subgraph.select`, `bits`, `label`, and `compass` return dictionaries, as do
  the selections in `graph.cycles` and `graph.forests`. *Raw selection bytes are
  no longer public inputs.* Pass complete objects to graph, layout, and drawing
  operations. Base62 labels preserve membership only, not annotations.
- `subgraph.with-data(g, s, data: value, hedge: updater)` stores subgraph-wide
  and per-hedge Typst annotations separately from graph data. Hedge callbacks
  see `edge-name`, structural `flow`, existing annotation `data`, and original
  half-edge `graph-data`. Content and callbacks stay on the Typst side.
- Ownership checks compare IDs, names, and source/sink incidence, not object
  identity, drawing orientation, placement, or user data. Layout does not
  invalidate a selection; changed topology requires selecting on the new graph.
- Cut origins remap graph, node, edge, and hedge sidecars without serialization.
  Node, edge, and hedge `origin` refer to the immediate input graph, including on
  repeated cuts. Edge origins contain `(edge, name, segment, winding)`; node and
  hedge origins are input IDs. Newly generated nodes have `origin: none`.
  Middle source hedges copy the input sink's data, and middle sink hedges copy
  the input source's, following their sides.
- Fragments use `<original.0>` through `<original.n>`. Unnamed edges use the
  native prefix `__linnest_cut_edge_ID`; generated name collisions are errors.
  Patch fragment-specific references such as `crossing-under: <D6.1>` rather
  than retaining a reference to a split master edge.
- Boundary nodes and dangling cut-edge records expose `boundary`, including
  current edge/node IDs, the creation input's side hedge, `side`, `crossing`,
  hedge annotation `data`, subgraph-wide `cut-data`, and fragment `origin` at
  creation. On subsequent cuts, `boundary.origin` and `boundary.hedge` retain
  creation provenance, as do `side`, `crossing`, `data`, and `cut-data`; only
  `boundary.edge` and `boundary.node` remap to current anchors. The `boundary`
  callback returns ordinary `graph.map`-style data/placement patches *only for
  newly created endpoints*, leaving inherited boundary placements and
  annotations intact.
- `graph.boundaries(view)` adds *current* endpoint `pos`: the dangling edge's
  position or the middle endpoint's auxiliary node position. Query after layout
  or inside `draw-after`; do not cache cut-time positions for the shaded boxes.
  `graph.nodes` includes auxiliary nodes, identifiable by non-`none` `boundary`.

Each view explicitly builds from the shared master items, opens its directed cut
with `graph.cut`, and applies geometry and depth pins with `graph.map`.
`diagram(g, ...)` receives that prepared graph and only styles, lays out, and
draws it, including the shared overlays. Boundary and hidden-node metadata
exclude auxiliary nodes from the cut-line bounds without recognizing node names.

Cut fragments reset old geometry rather than copying a paired edge's midpoint
pin to both new endpoints. Compactification springs remain layout-only builder
items, with their strength-independent length multiplier patched per view.
CeTZ's existing `on-layer(-1, ...)` and `rect` primitives provide background
boxes; no new shape API is needed.

Selection is still distinct from opening: `graph.nodes(g, subgraph: s)` returns
incident nodes, while selected edge records retain both available endpoints.
`draw(g, subgraph: s)` highlights half-edges; `layout(g, subgraph: s)` restricts
placement/optimization. Neither changes connectivity. Geometric `split-gap`,
`edge-halves`, and crossing gaps likewise do not create topological endpoints.

=== Unrelated join limitation

The existing `graph.join` wrapper still constructs its result with
`_empty-native-data()`. Earlier probes confirmed loss of node, edge, source, and
sink Typst payloads, including content-valued labels and momentum. The native
join also uses a left-biased edge-data merger. This is unrelated and unfixed:
the cut's explicit sidecar remapping does not make join its data-preserving
inverse. Replacing graph bytes alone cannot preserve sidecars after IDs or
cardinalities change.

== Sources and validation coverage

The public contracts live in `typst/src/graph.typ` and `typst/src/subgraph.typ`;
`typst/src/impl/graph.typ` and `typst/src/impl/subgraph.typ` own Typst-side
remapping and ownership checks. Native opening is `TypstGraph::cut_bytes` in
`crates/linnest/src/graph_api.rs`; automatic boundary-node styling and painting
are in `typst/src/impl/graph.typ` and `typst/src/impl/draw.typ`. Typst paths here
are relative to `crates/linnest`.

The earlier investigation queried all four diagrams to establish the physical
incidences, drawing orientations, and segment counts above; it also confirmed
the selection/query distinction and join data loss. Its native check and 13
focused Linnet tests are distinct from validation of the weighted-cut API.

Coverage is defined in the native weighted-cut tests and
`crates/clinnet/tests/resources/weighted-cut-behavior.typ`: winding and
orientation combinations, sidecar content/callbacks, origins, selection
ownership, repeated openings, current boundary positions, and invisible
auxiliary nodes. Check semantic records before comparing rendered pictures,
without overwriting the user's `xbox.pdf`.
