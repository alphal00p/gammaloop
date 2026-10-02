= Four-photon top-loop self-energy routing diagnostic

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.

This note records the coordinate obstruction in the GL262 projected local-4D
soft-counterterm route and the input-only basis change used to pass it. Source
citations refer to commit `032969abd45c53e0a5bb10d883d7b1d3a4bb149b`.
The campaign is `/tmp/soft-ct-four-photon-top-gluon-se-2026-09-30`;
`routing-certificate.json` records the transformation, graph hashes, and point
mapping. No production code was changed for this routing control.

== Graph identities and original coordinates

The source graphs are
`examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL262.dot`
and the adjacent `GL256.dot`. Both have four external photons, a top-quark
loop, and a self-energy insertion between gluon edges 5 and 7. The bubble
edges 9 and 10 are gluons in GL262 and ghosts (`ghG`) in GL256. The two graphs
retain their own vertices, directed edges, projectors, and prefactors:
GL262 has `overall_factor_evaluated = -2`; GL256 has `-4`.
See each DOT file, lines 3–5, 7–14, and 23–32.

Write `K0 = Q4`, `K1 = Q9`, and `K2 = Q12` for the original basis
`[4, 9, 12]`; `P0`, `P1`, and `P2` are graph-external momenta. The GL262
routes recorded in its DOT, lines 23–32, give:

#table(
  columns: (auto, 1fr, 1fr),
  table.header([Edge], [Original GL262 route], [Aligned GL262 route]),
  [4], [`K0`], [`h - q`],
  [5], [`-K0 + K2 + P0`], [`q`],
  [6], [`K2 + P0`], [`h`],
  [7], [`K0 - K2 - P0`], [`-q`],
  [8], [`K2 + P0`], [`h`],
  [9], [`K1`], [`b`],
  [10], [`-K1 - K2 - P0 + K0`], [`-b - q`],
  [11], [`K2 + P0 + P1`], [`h + P1`],
  [12], [`K2`], [`h - P0`],
  [13], [`K2 + P0 + P1 - P2`], [`h + P1 - P2`],
)

== The observed guard

The original GL262 H projected-4D attempt exited after 1.118 s with:

#quote(block: true)[
  component 12nI cannot preserve required nested expansion coordinates
  \[EdgeIndex(5), EdgeIndex(7)\]: their canonical routes carry graph-external
  momentum, so promoting them would be an affine loop-momentum shift
]

Evidence: `runs/gl262-h-projected-4d/generation/console.log` and
`receipt.json`. The guard is
`crates/gammalooprs/src/uv/approx/local_4d.rs:1935–1953`.
It rejects a completed child's required expansion carrier when that carrier
is internal to the parent but has a nonzero external signature in the
canonical full-graph route. The original routes of edges 5 and 7 contain
`+P0` and `-P0`, respectively.

The component label `12nI` decodes to the top self-energy edge set
`[4, 5, 7, 9, 10]`. Its nested bubble is `[9, 10]`, whose two boundary
representatives are edges 5 and 7. This identification follows the graph
ports and component bitset; the summary log does not itself include the
full current/given forest context. The existing GL262 inventory check
records the DOD-2 gluon bubble and enclosing DOD-1 top self-energy:
`tests/tests/uv.rs:2190–2221`.

The guard protects the already constructed child's Taylor coordinates.
The code requires them to remain in the parent's basis and restricts that
selection to homogeneous canonical routes
(`local_4d.rs:1956–2014`). This failure precedes outer CFF reconstruction
and is not evidence of a CFF recursion or contour-sign fault.

== Aligned basis and invertibility

The replacement basis is `[5, 9, 6]`, with

```
q = -K0 + K2 + P0,    b = K1,    h = K2 + P0;
K0 = h - q,           K1 = b,    K2 = h - P0.
```

The linear part has determinant `-1`, hence absolute Jacobian one for
the loop integration measure. Removing edges `[5, 9, 6]` leaves the
connected seven-edge tree `[4, 7, 8, 10, 11, 12, 13]` on eight vertices,
so this is also a valid graph loop-momentum basis. The transformation is
affine only through fixed graph-external momenta; it changes coordinates,
not the physical graph.

Now the bubble boundary is `Q5 = q`, `Q7 = -q`, and the top self-energy
boundary is `Q6 = Q8 = h`. Both can be represented by canonical loop
carriers without an external shift. `internal_expansion_carriers` recognizes
equal routes up to sign and deduplicates the representatives
(`local_4d.rs:2108–2136`). Choosing `[5, 9, 12]` would fix the first
bubble obstruction but would leave the top boundary `K2 + P0` affine;
including edge 6 addresses both nested levels.

The aligned DOTs set `lmb_id = 0` on edge 5, retain `lmb_id = 1` on
edge 9, set `lmb_id = 2` on edge 6, and remove the old assignments on
edges 4 and 12. Stale `lmb_rep` attributes are removed. The parser reads
`lmb_id` (`graph/edge.rs:598–601`), cuts those edges and preserves their
requested ordering (`graph/parse/mod.rs:1015–1037,1282–1350`). In contrast,
`lmb_rep` is an ignored presentation attribute
(`graph/attribute_warnings.rs:60–70`). All these graph source paths are
relative to `crates/gammalooprs/src/`.

The saved certificate checks that all other graph data are unchanged.
Aligned graph SHA-256 values are:

- GL262: `1093d0fd54eeda7724998f6bbdc1bd8b2298c86c9d31e682bf98acb1ae27fd10`.
- GL256: `df1294df520559df2ec73e411230331bc4308e00e8d4e49d64eb8e4643049cd7`.

The first aligned GL262 attempt passed integrand construction in 208.660 s
and evaluator parsing in 62.387 s. It then reached the 12 GiB process-tree
RSS cap during Symbolica evaluator construction, after 347.435 s total.
Thus the coordinate obstruction is cleared, but that attempt supplies no
completed evaluator, pointwise comparison, or soft-scaling measurement.
See `runs/gl262-h-projected-4d-aligned/generation/receipt.json`.
Fresh larger-memory attempts have separate case directories and receipts.

== Comparison contract and ghost sign mapping

Primary route comparisons must use one common basis, the same external
data, subtraction prescription, orientation sum, and physical momentum
point. The certificate maps the old point coordinates to the aligned
ones. Reusing the same nine numerical components in both bases would
evaluate different physical momenta.

The invertible transformation establishes a coordinate correspondence for
the original graph. It does not by itself prove pointwise equality of
local counterterms regenerated in the two bases. A finite Taylor
projection depends on which hard and external coordinates are held fixed;
an external-momentum-dependent coordinate shift need not commute with
that projection. A cross-basis claim would require transforming the
actual U/H operators and their retained child coordinates consistently.
The implementation explicitly fixes those coordinates before Taylor
projection and reuses them for residues
(`local_4d.rs:2455–2464,2496–2532,2567–2576`). This note asserts neither
cross-basis local-counterterm equality nor equality of integrated
subtraction terms, which were disabled in these runs.

The gluon and ghost graph routings also cannot be identified by copying
the same value of the second loop coordinate. Edge 9 is oriented
`2 -> 3` in GL262 and `3 -> 2` in GL256. With each graph's own definition
`b = Q9`, their aligned edge-10 routes are respectively `-b - q` and
`b - q` (DOT lines 28–29). To compare the same physically directed bubble
line momenta, use `b_ghost = -b_gluon`; do not silently reverse an edge
while retaining its numerator or vertex assignment. Denominator evenness
does not erase an odd numerator or ghost momentum-flow sign. Preserve
the original graph-specific factors `-2` and `-4`; their ratio is not a
replacement for evaluating the two distinct numerators. A combined
gluon-plus-ghost soft test therefore needs this explicit point/ray mapping.

An explicit physical central soft ray for GL262 is `q -> lambda*q` with `b` and
`h` fixed, `lambda > 0`. In the original chart it is
`K0 -> K0 + (1-lambda)*q` with `K1,K2` fixed. Consequently both central
gluon lines soften together (`Q5 = lambda*q`, `Q7 = -lambda*q`) and
`Q10 = -b-lambda*q`; there is one independent soft momentum, not two.
This is the generic canonical ray used for the matched graph-pair
diagnostic. It is not the native bulk profiler's particular ray: that
profiler selected the basis `[5, 10, 13]` and holds its companion
coordinates fixed before routing back to `[5, 9, 6]`.
