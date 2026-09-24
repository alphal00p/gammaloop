= Native SPQR decomposition protocol

This helper supplies OGDF's exact planar SPQR decomposition to the implementation
of Gutwenger, Klein, and Mutzel's embedding-constraint algorithms. It does not
implement embedding constraints or optimal edge insertion itself.

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
the Python layer makes that choice.

`build.py` downloads checksum-verified pinned dependencies and builds a static
helper at `target/ec-planarity-native/ec-spqr`. It needs Python 3.11+, CMake,
and a C++17 compiler. Build caches and downloaded sources stay under `target/`.
Set `EC_SPQR_BINARY` to override this executable in the Python caller.
