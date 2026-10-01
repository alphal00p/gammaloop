= Exact single-orientation forest comparison

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<exact-single-orientation-forest-comparison>
The saved U and H cases have #strong[identical 16-node, 16-edge structural forests];, including the root, and exactly the same parent-key lists in computed exports. Both are localized direct 3D with the single generation pattern `(0,0,-,-,-,+,-,-,+,-)`, not explicit-sum generation. Their structure-only DOT files are byte-identical, SHA-256 `c9230e9e50a5ee3e61a03770853d335d108729849cd865b242345967206086ac`.

The central bubble `[e8,e9]` is literal zero on this selected orientation in both U and H, and every forest with that child prefix is zero. A second region `[e2,e5,e6,e7,e8,e9]` is also zero. In total #strong[10 of 16 node terms are literal zero] in each case. This is a material limit of this cheap control: it does not exercise a nonzero nested central soft child and cannot predict complete-orientation H complexity.

== Node expression sizes
<node-expression-sizes>
The existing direct-3D provenance schema contains routes, schemes, conceptual branches and materialization choices, but #strong[no packed-atom byte-size fields];. The table therefore measures UTF-8 byte length of each exported `full_num` expression as rendered by the existing Typst printer. These are exact serialization lengths, not Symbolica packed sizes, operation counts, evaluator sizes or resident memory. Different strings need not represent different mathematical functions. No graph numerator was expanded to compute these lengths.

#figure(
  align(center)[#table(
    columns: 5,
    align: (auto,auto,right,right,right,),
    table.header([Forest key], [Meaning], [U text bytes], [H text bytes], [H/U],),
    table.hline(),
    [`∅`], [Root], [1,183], [1,183], [1.000],
    [`{15zU}`], [`[e3,e4,e6,e7,e8,e9]`, ordinary logarithmic component], [1,802], [1,802], [1.000],
    [`{GS}`], [Outer fermion square `[e2,e3,e4,e5]`, ordinary logarithmic component], [1,246], [1,246], [1.000],
    [`{16C4}`], [Overall degree-two component `[e2..e9]`], [63,629], [103,727], [1.630],
    [`{15zU} · {16C4}`], [Ordinary logarithmic prefix, then overall component], [17,632], [20,683], [1.173],
    [`{GS} · {16C4}`], [Outer-square prefix, then overall component], [8,761], [1,821], [0.208],
    [Ten literal-zero nodes], [Each prints `0`], [10], [10], [1.000],
    [#strong[Total];], [Six nonzero and ten zero node terms], [#strong[94,263];], [#strong[130,472];], [#strong[1.384];],
  )]
  , kind: table
  )

Thirteen of sixteen node-expression strings are byte-identical between U and H. Only the three keys ending in the overall component shown above differ. The largest node is `{16C4}` in both cases. Its provenance changes from scheme MUV / `direct_u` to IR / `factorized_soft`, with conceptual coefficients `U:+1`, `S:+1`, `US:-1`. That records a factorized completed expression, not three additional structural forest nodes. The other nonzero primitive components remain MUV.

Subgraph keys were decoded using the repository\'s exact base-62 hedge-bitset encoding (`crates/linnet/src/half_edge/subgraph/subset.rs:502`) and matched to endpoint hedge IDs in the exported original graph. This is metadata decoding, not a numerical or symbolic approximation:

#figure(
  align(center)[#table(
    columns: 2,
    align: (auto,auto,),
    table.header([Key], [Internal edges],),
    table.hline(),
    [`11vs`], [8, 9],
    [`15zU`], [3, 4, 6, 7, 8, 9],
    [`168C`], [2, 5, 6, 7, 8, 9],
    [`GS`], [2, 3, 4, 5],
    [`16C4`], [2, 3, 4, 5, 6, 7, 8, 9],
    [`12CK`], [2, 3, 4, 5, 8, 9; disconnected union],
  )]
  , kind: table
  )

== Execution and limits
<execution-and-limits>
The existing `save uv-forest` help exposes only `--computed`, with no metadata-only or size-only option. Structure-only exports do not construct expressions. Computed ThreeD exports recompute the forest using the stored selected source orientations and necessarily write `full_num` text beside the atom-free provenance. These small expressions remain in the audit\'s `/tmp` directory and were not printed into the report or tool output. The computed export does not build a second 4D representation or a new numerical evaluator.

All four read-only commands exited zero under a 60-second / 2-GiB sampled-process-tree RSS guard. Structure-only U/H took 0.806/0.706 seconds; computed U/H took 1.108/1.209 seconds, with maximum sampled tree RSS below 51 MiB. These export costs are diagnostics, not generation benchmarks. Loaded states were not modified. No test or source code changed.

Evidence: commands and resource records;, all sixteen matched keys, exact lengths, parent relationships, component maps and complete atom-free provenance;, plus each case\'s structure/computed DOT files in the adjacent directories. The equality certificate concerns structural topology and recorded strings for this one selected orientation. It does not establish general H/U equivalence, strict soft convergence, or the complexity of the unfiltered source sum.
