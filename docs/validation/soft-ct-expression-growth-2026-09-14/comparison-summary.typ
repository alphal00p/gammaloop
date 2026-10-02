= Bounded same-graph counterterm complexity comparison

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<bounded-same-graph-counterterm-complexity-comparison>
Same imported three-loop Figure B.1 graph, original numerator/model/runtime, ThreeD final integrand, no integrated addback or thresholds. U clears both IR overrides and explicitly selects MUV for logarithmic, massive-power and massless-power prescriptions; H keeps both original IR overrides. Localized singles retain `(0,0,-,-,-,+,-,-,+,-)`; explicit all-orientation runs clear the filter because the existing API rejects explicit sums with a filter. Raw CFF source-map inventories include 102 maps even for the localized single-orientation runs. No numerical profile was requested.

Every primary cell has a 180-second wall bound and 12 GiB descendant-RSS bound, polled every 0.2 seconds, with TERM of its own process group on either cap. The separately authorized U extension uses 600 seconds with the same memory bound; the all-bare control uses 60 seconds and 2 GiB. RSS table uses the largest process VmHWM observed by the monitor; transient exit peaks may exceed the final sampled value.

#figure(
  align(center)[#table(
    columns: 10,
    align: (auto,auto,right,right,right,right,right,right,right,right,),
    table.header([Run], [Outcome], [Wall s], [Max observed process RSS MiB], [T projections], [Largest T output bytes], [Evaluator summands], [Scalar aliases], [Max resolved scalar term bytes], [Final scalar result bytes],),
    table.hline(),
    [single-u-logged], [completed], [6.018], [417.125], [20], [98873], [6], [4], [19925657], [25233258],
    [single-h-logged], [completed], [6.217], [531.605], [33], [140463], [6], [4], [27926697], [32898185],
    [all-u-explicit-logged], [wall\_exceeded\_180\_seconds], [180.191], [698.637], [846], [1409850], [not reached], [not reached], [not reached], [not reached],
    [all-h-explicit-logged], [wall\_exceeded\_180\_seconds], [180.299], [622.199], [1082], [1822911], [not reached], [not reached], [not reached], [not reached],
    [all-bare-explicit], [completed], [2.407], [80.988], [0], [---], [102], [0], [34540], [3315713],
    [all-u-explicit-600s-after-bare], [rss\_tree\_exceeded\_12\_GiB], [412.854], [12454.098], [1392], [1409850], [1074], [887], [773457405], [not reached],
  )]
  , kind: table
  )

The initial U/H attempts completed but lacked semantic trace files because the CLI -t argument is a suffix, not a path; their all-U continuation was explicitly cancelled for the logging repair. Those artifacts remain under single-u, single-h and all-u-explicit; the first 600-second U attempt was subsequently cancelled to insert the requested bare control, and is retained separately under all-u-explicit-600s; they are excluded from the stage comparison above. Working traces use -t trace and reside at each absolute state/logs/gammalog-trace.jsonl.

Exact commands, copied cards, overrides, resource samples and structured stage events are retained alongside this report. Scalar result bytes measure Symbolica atom serialization size, not peak RAM. Native source maps are not forest-node counts. The single-orientation ratios do not predict all-orientation behavior. The independent forest audit found 16 nodes in both localized singles, with 10 literal-zero terms; the central bubble prefix \[e8,e9\] and every descendant containing it are zero for this selected orientation, so this ratio does not test a nonzero central-soft branch. See forest-audit/report.md for the exact per-node inventory.

== Completed interpretation
<completed-interpretation>
The all-bare explicit control finishes with the same 102-source-map inventory in 2.407 seconds, 102 evaluator summands and 3,315,713 scalar-output bytes. Adding ordinary U on the same graph increases the evaluator inventory to 1074 summands; source assembly completes in 245.548 seconds after 1392 Taylor projections. Positive-degree soft operators are therefore not necessary for the observed large scalar materialization in this checkout.

The extended ordinary-U run stopped on its 12 GiB RSS bound, not its 600-second wall bound. The monitor observed 12,752,996 KiB (12.162204742 GiB) at 411.587407511 seconds against the limit of 12,582,912 KiB, an overshoot of 170,084 KiB (166.097656 MiB) between samples. TERM completed at 412.854073988 seconds with exit −15. The observed process HWM equals the sampled peak; the monitor polls every 0.2 seconds rather than imposing an address-space limit.

At cancellation, 290 of the 1074 scalar terms had completed alias resolution. The largest completed term occupied 773,457,405 serialized bytes; completed term sizes sum to 5,970,393,977 bytes. That cumulative sum is not the size of a final combined atom, because generation was stopped before the complete sum or saved state existed. The last completed log event is parsing of term index 290.

One concrete ordinary-U scalar boundary is term256: six scalar aliases cover 39,612 bytes (largest alias 8,499 bytes), one closed-sum boundary is identified, tensor execution leaves a 115,846,127-byte atom, and alias resolution produces 129,014,497 bytes. This locates a large growth step in the shared evaluator pipeline; the counters alone do not establish an algebraic correctness defect.

The full-H 180-second trace also stops before evaluator preprocessing, so no full-H/full-U scalar-size ratio is available from these bounded cells. The separate nonzero-child orientation audit is required to quantify the soft-specific increment; the zero-child single-orientation ratio cannot supply that result.

All controls use the current validation binary and its shared code. Clearing IR overrides selects ordinary U in that implementation; it does not restore a historical pre-soft implementation. A historical pre-soft performance regression remains unmeasured.

The ordinary-U prescription is reproducible from muv.json: log\_divergent=MUV, massive\_power\_divergent=MUV, massless\_power\_divergent=MUV, overrides=\[\]. The command applies this file before run generate. This follows the existing ordinary-prescription control in #link("../../../tests/tests/uv.rs:4911")[uv.rs];. Graph/card hashes, binary modification time, all absolute command paths and limits are retained in the input/command/result JSON files.
