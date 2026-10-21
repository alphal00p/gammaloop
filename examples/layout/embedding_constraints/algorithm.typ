= Constrained embedding and optimal edge insertion

This experiment implements the embedding-constraint and edge-insertion algorithms
from Gutwenger, Klein and Mutzel,
#link("https://jgaa.info/index.php/jgaa/article/view/paper160")[Planarity Testing and Optimal Edge Insertion with Embedding Constraints],
JGAA 12(1), 2008. It supplies a combinatorial initialization for the diagram
layout experiment. The current playground constructs an EC seed with metric
coordinates, then relaxes it with ImPrEd. It stops there; native spring refinement
and annealing are excluded from this pipeline. Constraint projection remains
active during ImPrEd.

== Metric refinement forces

The playground uses ImPrEd's topology protection, exponent schedule and
certified route subdivision/contraction. Attraction acts on each segment
independently, with tension $T = A ell (ell / delta)^a$ for segment length
$ell$. A new route point adds a segment in series: under a fixed endpoint
pull, an edge can become longer as it gains points. On a straight section,
unequal adjacent lengths produce unequal tensions, moving their shared point
toward equal spacing.
The exponent decreases from $1$ to $0.4$ during refinement. The spacing
parameter $delta$ sets the force scale; it is not a prescribed rest length.

Each external tip and its intermediate points share a total point charge of
one; internal edge samples share a total charge of `0.35`. Physical vertices
keep unit charge. Points on the same flexible dangling connection do not repel
one another. Repulsion against other connections remains active, and flexible
internal edges retain their consecutive-point exemption. Subdividing an
external leg therefore adds compliance without adding total point charge or
introducing repulsion between points on that same connection.
Fixed charge conserves the total, not the spatial force field: moving or adding
samples can still change repulsion against other connections.

The constant external pull uses native topology normalization and active
coordinate groups. The playground exposes separate controls for its base
strength, balancing exponent, and distributed-attachment correction. The last
control increases demand from shared external X groups spread through the
network without increasing terminal-group or one-direction radial pull.
These force choices extend the paper's metric model; they do not change the EC
embedding algorithm, crossing guards, or certificates required to add or remove
route points.

=== Parallel-edge attraction balancing

Parallel physical edges contribute multiple attractive routes between the same
vertices. The experimental `parallel_attraction_balance` parameter applies
$A_e = A m_e^(-beta)$ to every segment of such an edge, where $m_e$ counts
internal edges sharing its two endpoints. Direction and intermediate drawing
points do not affect this count. Single edges, dangling edges and self-loops
keep their original attraction; charges, spacing and external pull are unchanged.

At `0`, the accepted unbalanced force is reproduced. At `1`, identical parallel
routes together exert the attraction of a single route. Intermediate values
partially compensate this multiplicity, allowing the parallel loop to grow
relative to adjoining single-edge bridges. This does not guarantee equal sizes:
repulsive forces and the directions of curved segments still matter.
It recognizes direct parallel edges, not arbitrary self-energy subgraphs with
additional physical vertices. In split cross-section drawings, a visible bridge
need not be a bridge of the unsplit physical graph.

== Constraint semantics

A constraint tree contains every edge incident to its vertex exactly once.
Its leaves are edge identities, not neighboring vertex identities. Thus parallel
edges remain distinct. A `group` node permits any permutation of its children,
a `mirror` node permits the listed child order or its reversal, and an `oriented`
node permits only the listed order. The complete vertex order is cyclic and
clockwise. These are the paper's gc-, mc- and oc-nodes, respectively.

Each internal node has at least two children. At a constrained vertex, no
incident edge may be omitted. Self-loops must first be subdivided into a
drawing-only cycle; subdivision does not change the physical graph.

The tree constrains incidences, not Cartesian coordinates. Incoming and outgoing
external groups can define an exterior scaffold. Equal X coordinates, X signs,
straight dangling legs, and any native fixed coordinates still need a metric
construction and independent validation.

The metric adapter uses temporary integer vertex identities during grid drawing
so NetworkX's canonical-ordering choices do not depend on Python's randomized
string hashes. Original graph and edge identities survive this conversion.

== Algorithm correspondence

- `ECExpansion` implements Section 4: replace a constrained vertex by its tree,
  and replace each mirror/oriented node of degree $d$ by a wheel with $2d$ rim
  vertices and spokes. Expansion edges cannot be crossed.
- `embed_expanded_graph` implements Algorithm 1: decompose the expansion into
  blocks, compute planar SPQR skeletons, and reject conflicting oriented hubs
  within an R skeleton. Independently orient the remaining R skeletons and
  assemble blocks through faces outside wheel interiors.
- `optimal_insertion` implements Algorithms 2--4: follow the allocation path in
  the SPQR tree, retain both exit-side costs and their backpointers, allow the
  paper's P permutations and permitted R reflections, then concatenate the
  block-path solutions. Local paths use dual BFS excluding protected edges.
- Mirrored local R cases are reconstructed and searched with off-path children
  retaining their valid orientation. This realizes the paper's starred-path
  transformation without reflecting oriented hubs in child expansions.
- Endpoint constraints use the expansion of the #emph[full] graph, with only the
  inserted physical edge removed. Its reserved leaf ports remain present.
- `planarize` retains a greedy maximal EC-planar subgraph, then inserts each
  remaining edge optimally. Existing crossings are expanded into mirror wheels
  for later insertions so their alternating incidences are preserved.

Each #emph[single insertion] minimizes its crossing count over all admissible
embeddings of the current graph. The sequence does not guarantee a globally
minimum crossing number. The Python implementation uses explicit graph copies
and reconstruction; it implements the paper's optimization but does not claim
its linear running-time bound. No enumeration of embeddings is used by the
solver. Exhaustive enumeration is confined to small independent test oracles.

== Native primitive and reproduction

The C++ helper uses OGDF's planar SPQR decomposition. OGDF's unconstrained
variable-edge inserter and its separate SyncPlan algorithm are not substituted
for the paper's EC algorithms. `native/build.py` pins dependency revisions and
verifies SHA256 checksums. Downloaded source, licenses, and build outputs stay
under `target/ec-planarity-native`; no third-party source is vendored here.

From the experiment directory, with `uv`, CMake and a C++17 compiler available:

```sh
cd examples/layout/embedding_constraints
uv run python native/build.py
uv run python -m unittest discover -v
uv run python -m unittest discover -s native -v
```

Set `EC_SPQR_BINARY` to use a prebuilt helper. The helper protocol is documented
in `native/protocol.typ`. The implementation and tests require no physics model
or Symbolica license. The metric adapter pins NetworkX 3.7 through `pyproject.toml`.
The topology fixtures in `diagrams.json` include all eight
existing diagram examples, including the nonplanar K3,3 and K5 examples.

The CLI input uses string identities:

```json
{
  "nodes": ["a", "b", "c"],
  "edges": {"ab": ["a", "b"], "bc": ["b", "c"], "ca": ["c", "a"]},
  "constraints": {"a": {"kind": "oriented", "children": ["ab", "ca"]}},
  "protected": []
}
```

```sh
uv run python cli.py problem.json --mode embed
uv run python cli.py problem.json --mode insert --insert ab
uv run python cli.py problem.json --mode planarize --output planarized.json
```

Embedding returns physical clockwise rotations. Insertion and planarization
return the expanded planar drawing graph, directed edge chains, crossing
provenance, and the original-to-expansion map. Crossing vertices are temporary
drawing objects, never physical interactions.
