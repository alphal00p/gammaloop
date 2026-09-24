= Constrained embedding and optimal edge insertion

This experiment implements the embedding-constraint and edge-insertion algorithms
from Gutwenger, Klein and Mutzel,
#link("https://jgaa.info/index.php/jgaa/article/view/paper160")[Planarity Testing and Optimal Edge Insertion with Embedding Constraints],
JGAA 12(1), 2008. It supplies a combinatorial initialization for the diagram
layout experiment. The current playground constructs an EC seed with metric
coordinates, then relaxes it with ImPrEd. It stops there; native spring refinement
and annealing are excluded from this pipeline. Constraint projection remains
active during ImPrEd.

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
