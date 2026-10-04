= Native constrained planarization protocol

The native engine implements Gutwenger, Klein, and Mutzel's embedding-constraint
and optimal edge-insertion algorithms using OGDF's SPQR decomposition. The
`initialize` operation builds a constrained planarization and its physical edge
routes; the `decompose` operation exposes the underlying SPQR primitive for the
Python reference implementation.

== Runtime and build

Direct Typst, Python rendering and CLI drawings call `graph_impred_seed` in
`linnest.wasm`: a Rust port of this engine (`crates/linnest/src/ec_seed`) with
the same request and result protocol, whose own SPQR decomposition replaces
OGDF. This native engine remains the reference: given its SPQR
decompositions, the port reproduces its seeds bit for bit. The standalone
`ec-layout.wasm` plugin imports only the two Typst plugin protocol functions.

`build.py` downloads checksum-verified pinned dependencies and builds a native
JSON-lines executable at `target/ec-planarity-native/ec-spqr`. It needs Python
3.11+, CMake, and a C++17 compiler. `build.py --wasm` also needs Emscripten and
Binaryen, and produces `ec-layout.wasm`. Build caches and downloaded sources stay
under `target/`. Set `EC_SPQR_BINARY` to override the executable used by the
Python reference implementation.

== Seed initialization

The WASM `initialize` export accepts JSON with `diagram`, optional `scale`
(default `2.4`), and optional `external_sides` (default `true`). The native
executable dispatches to initialization when the request contains `diagram`.
The diagram contains native graph topology, external flow and embedding
constraints; `graph_impred_diagram` prepares it from a Linnest graph snapshot.

The result contains vertex and route-point `positions`, physical `routes`,
`node_ids`, `external_ids`, `edge_endpoints`, explicit `crossings`, and a geometry
`report`. The constrained embedding preserves cyclic endpoint groups;
nonplanar inputs use individually optimal edge insertions. The insertion
sequence is a heuristic for the complete graph's crossing number. Artificial
crossing vertices have paired route-point identities belonging to their
physical edges. Downstream ImPrEd preserves this crossing configuration.

== SPQR decomposition

Each stdin line is one JSON request; exactly one JSON response is written to
stdout. EOF exits cleanly. A malformed request produces `{"error":"..."}` and
does not terminate the session. Diagnostics belong on stderr.

```json
{"nodes":["a","b","c"],"edges":[{"id":"ab","source":"a","target":"b"},{"id":"bc","source":"b","target":"c"},{"id":"ca","source":"c","target":"a"}],"root_edge":"ab"}
```

Node and original-edge IDs are unique strings. Parallel edges are supported;
self-loops must be handled before decomposition. `root_edge` is optional.
The caller separates biconnected blocks and handles blocks with fewer than three
edges. Nonplanar or non-biconnected inputs return their respective Boolean flags
and an empty `components` array.

The result contains `planar`, `biconnected`, `root`, `components`, and
`rotation_cw`. The latter maps original vertex IDs to original-edge IDs, in
clockwise cyclic order. Each component has:

- `id`: a string identifying the SPQR tree node;
- `type`: `S`, `P`, or `R` (OGDF omits Q nodes);
- `reference_edge`: a skeleton-edge ID, real at the root and virtual elsewhere;
- `nodes`: original graph vertex IDs;
- `edges`: records with `id`, original `source` and `target` vertex IDs,
  `real_edge` (original edge ID or `null`), and `twin` (either `null` or a record
  containing `component` and `edge` IDs);
- `rotation_cw`: original vertex IDs mapped to skeleton-edge IDs in clockwise
  cyclic order.

Skeleton-edge IDs are globally unique within a response: `component:local-edge`.
Virtual twins have the same two original endpoints, but their source/target
directions need not agree. Treat endpoint identity explicitly during splicing.
Clockwise output reverses OGDF's documented counterclockwise adjacency order.
It does not choose an absolute orientation satisfying an oriented constraint;
the constrained embedding layer makes that choice.

