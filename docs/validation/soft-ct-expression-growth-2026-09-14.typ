= Soft-counterterm expression-growth investigation --- 2026-09-14

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<soft-counterterm-expression-growth-investigation--2026-09-14>
#strong[Large scalar-expression amplification is reproducible with soft subtraction disabled in the current checkout.] The earlier 314 GiB failure does not establish that soft counterterms inherently require expressions of that size. A possible Symbolica length overflow explains a crash mechanism, not the upstream growth. The matched checks below isolate ordinary UV subtraction (U) from the soft prescription (H) on the #strong[same three-loop massless graph];. This is a supplement to the #link("soft-ct-2026-09-14.typ")[full validation report];, not a replacement for its failed full-orientation acceptance result.

== Inputs and pipeline
<inputs-and-pipeline>
The shared input is `paper_figure_b1_double_triangle_soft_ir`: one three-loop graph, eight internal edges, the same Standard Model numerator/projector, kinematics, loop-momentum basis, disabled thresholds/integrated counterterms, direct 3D generation, and the same evaluator settings. H keeps both original positive-degree IR overrides; U clears both overrides and explicitly selects MUV for all three divergence classes. The legacy `uv.softct` field is not a production switch. A bare control instead sets `subtract_uv=false` and removes the forest.

Ordered pipeline: stored parent CFF and exact residue maps → classified forest → component-owned numerator mapping/common OSE setup → #strong[first U/H difference: Taylor operation] → branch normalization → orientation/forest summation → finite tensor contraction and scalar-alias restoration → scalar evaluator construction. U applies U(X); positive-degree H already uses S(X)+U(X−S(X)). This changes coefficient algebra without adding U/S/US forest nodes. Logarithmic H takes U directly. See the #link("soft-ct-expression-growth-2026-09-14/current-pipeline-audit.typ")[current pipeline audit];.

All recorded raw-CFF generation events have #strong[102 native source-map entries];. The localized control retains the same audited single orientation in U and H. Full explicit mode retains all generated maps. An explicit single-orientation configuration is rejected by the existing settings API; it was not substituted silently for the localized control. No selector-local/full-explicit correctness equivalence is inferred from unequal inventories.

== Completed single-orientation control
<completed-single-orientation-control>
#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,right,),
    table.header([Measured quantity], [Ordinary U], [Soft H],),
    table.hline(),
    [Structural forest nodes / edges], [16 / 16], [16 / 16],
    [Exactly zero exported node expressions], [10], [10],
    [Taylor projection calls], [20], [33],
    [Largest completed Taylor output, packed atom bytes], [98,873], [140,463],
    [Top-level summands entering finite contraction], [6], [6],
    [Final preprocessed scalar, packed atom bytes], [25,233,258], [32,898,185],
  )]
  , kind: table
  )

The forest DOTs are byte-for-byte identical, and every node\'s parent keys agree. Thirteen of sixteen rendered node expressions are identical, including the ten literal zeros. Only three containing-component expressions differ. The total rendered node-expression text grows from 94,263 to 130,472 UTF-8 bytes (1.384×); this text-size metric is distinct from packed atom size and process RSS. The largest containing node grows 1.630×, another grows 1.173×, and the third shrinks to 0.208×. See exact forest comparison;.

The final scalar grows #strong[1.304×];, while the largest individual Taylor output grows 1.421×. This selected orientation does not reproduce the catastrophic growth. The central bubble \[e8,e9\], key `{11vs}`, is one of the ten literal-zero nodes, and every descendant containing this prefix is also zero. The only changed coefficients involve the overall \[e2…e9\] component `{16C4}`. This selected orientation does not exercise a nonzero nested central-soft-child coefficient and is a poor predictor of the full source inventory. It also demonstrates that finite contraction materially enlarges both representations; this observation alone does not identify a faulty contraction or prove lost factorization.

Both fresh stage-logged controls complete successfully in 6.018 s (U) and 6.217 s (H), with observed process high-water 417.125 and 531.605 MiB. Their raw traces and run receipts are retained: U trace;, H trace;, U receipt;, H receipt;. These are diagnostic generation runs on a shared machine, not repeated throughput benchmarks.

== Full-orientation controls
<full-orientation-controls>
The same full 102-map graph with `subtract_uv=false` completes in #strong[2.407 seconds];, with a preprocessed scalar of #strong[3,315,713 packed bytes] and observed process high-water of #strong[82,932 KiB];. It retains 102 independent scalar summands and needs no scalar aliases. The un-subtracted full parent CFF is therefore small and quick to process in this environment. See bare receipt and trace;.

The initial full-U and full-H diagnostic probes each stop at their 180-second wall bound #strong[before evaluator preprocessing];. Their largest observed individual Taylor outputs are 1,409,850 and 1,822,911 bytes respectively, with observed process high-water 715,404 and 637,132 KiB. These are incomplete trace prefixes at potentially different points in the forest, not a completed U/H size or speed ratio. Their last events are in final forest-node construction; no result from these bounded probes replaces the original full-H production panic.

The longer ordinary-U probe completes source assembly in #strong[245.548 seconds];, then enters evaluator preprocessing with #strong[1,074 summands];. It crosses the 12 GiB diagnostic memory bound and is terminated after #strong[412.854 seconds];, during tensor execution of zero-based summand 290. At that point 290 summands have completed, their sizes sum to #strong[5,970,393,977 packed bytes];, and the largest individual scalar term is #strong[773,457,405 bytes];. The cumulative size is not a completed combined atom, and no final U evaluator or numerical result was obtained.

The monitor sampled 12,752,996 KiB RSS against a 12,582,912 KiB bound: 170,084 KiB overshoot between 0.2-second polls. Its termination is a deliberate diagnostic memory stop, not a Symbolica panic or an OS OOM kill. See the final receipt;, compressed complete trace;, and #link("soft-ct-expression-growth-2026-09-14/comparison-summary.typ")[full measurement matrix];.

== Smaller reproducer: source orientation 58
<smaller-reproducer-source-orientation-58>
The existing orientation display independently lists `(0,0,+,+,-,-,+,+,+,-)` as source orientation 58. A copied card changes only the selected localized orientation; U/H then differ only in the same prescription overlay used above. This supplies a much smaller reproducer than generating every source map. See the #link("soft-ct-expression-growth-2026-09-14/orientation-58.toml")[card] and actual orientation inventory;.

Both U and H cross the initial 2 GiB RSS bound during #strong[finite tensor execution];, after parsing zero-based summand 5 and before its execution-complete event. Numeric summand indices are not stable coefficient identities across schemes or repeated runs; this locates the shared pipeline stage, not an identical U/H input term. All local Taylor projections had already completed. Repeating this one-orientation pair with 60-second/12-GiB bounds reaches the completed preprocessed scalar in both cases:

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,right,),
    table.header([Quantity, orientation 58], [Ordinary U], [Soft H],),
    table.hline(),
    [Taylor projection calls], [20], [33],
    [Largest individual Taylor output, packed bytes], [591,473], [803,832],
    [Top-level scalar summands], [14], [14],
    [Completed preprocessed scalar, packed bytes], [394,216,189], [541,813,987],
    [Finite preprocessing time], [14.672 s], [21.212 s],
    [Observed process high-water], [6.810 GiB], [11.439 GiB],
    [Final diagnostic outcome], [60 s wall bound], [60 s wall bound],
  )]
  , kind: table
  )

A separate bounded computed-forest export certifies the structural coverage of this orientation: U/H have the same 16 nodes and parent lists, with only two literal-zero coefficients. The central bubble `{11vs}`=\[e8,e9\] is represented by 7,854/8,815 UTF-8 bytes in U/H. The largest coefficient is its nested child→overall key `{11vs}·{16C4}`, at 372,890/592,162 rendered bytes. A non-literal-zero expression is evidence of a represented contribution, not a mathematical noncancellation certificate. See the exact forest comparison;.

These forest exports used scratch copies of the original small saved U/H states with #strong[only the saved generation-history orientation pattern changed];. The existing export API recomputes the forest from the stored full production CFF and those generation settings. It does not use the stale numerical evaluator for the computation. The copied states were diagnostic export inputs, never claimed to be newly generated orientation-58 evaluators; original states remained untouched. Both exports succeeded in 3.213/4.617 seconds below 66 MiB sampled RSS. Export receipts;.

The H scalar is #strong[1.374×] the U scalar at this same orientation. Both are hundreds of MB even though every individual local Taylor output is below 1 MB. These quantities compare different representation boundaries; their quotient is not a per-term growth factor. Both runs subsequently enter Symbolica evaluator construction but stop at their wall bound before construction finishes. The Symbolica build inputs are 133 bytes larger than the preprocessed scalars because intervening expression preparation changes the representation; those byte counts are deliberately kept separate.

The defaults used here are `do_algebra=false`, sequential finite tensor execution, `MinResultRank` contraction mode with the `IntermediateCost` ordering, and `compile=false`. The first reproduced large-memory boundary is inside `net.execute(...)` in #link("../../crates/gammalooprs/src/integrands/process/evaluators.rs:747")[preprocess\_atom];, before scalar-alias restoration for the failing summand. The later selector scan seen in the original long run consumes expressions that are already large. This does not identify the exact contraction pair or prove an algebraic error, lost factorization, or historical performance regression.

Exact narrowed receipts: U;, H;. Complete compressed stage traces: U;, H;.

The measured problem is therefore shared counterterm scalar materialization in the current tree, with an additional H cost. The exact historical introduction remains open. The smallest next implementation-level comparison is the factorized input/network for orientation 58, aligned by forest key and exact source map before comparing contraction steps. The nested central-child→overall key `{11vs}·{16C4}` is the largest exported coefficient and is a concrete starting point; matching a numeric summand index alone is insufficient. Preserve numerator factors and prove exact output equivalence for any proposed representation change; do not infer a repair from memory measurements alone.

== Historical comparison and exclusions
<historical-comparison-and-exclusions>
The ordinary two-loop Figure B1 control is a different graph and cannot establish a soft-only complexity multiplier. The exact three-loop fixture/card are absent from the pre-stack base, so prior acceptance on other fixtures does not constitute a matched historical measurement.

The pre-stack ordinary U routine already used hard-loop scaling and a Laurent series. The removed opaque OSE derivative helper belonged to the former soft S path. Current U has changed canonical routing, external-shift authority and scoped vacuum/expansion-mass representation; a historical old-U/current-U comparison would be needed to quantify those changes. The production evaluator constructor, preprocessing, selector scans and scalar-alias restoration already existed before this stack. Their presence in a hot stack does not identify a newly introduced evaluator regression. See the #link("soft-ct-expression-growth-2026-09-14/historical-pipeline-audit.typ")[historical source audit];.

== Execution and scope
<execution-and-scope>
No production source, failing test, assertion, or graph-numerator factorization was changed during this investigation. The two previously approved test API repairs remain in the parent validation change. Diagnostic cards, generated states, monitors and raw logs are isolated under `/tmp/soft-ct-complexity-2026-09-14`; explicitly saved states use absolute paths. Forest exports reload those states read-only and recompute only the selected forest, without rebuilding evaluators.

The first U/H timing runs completed, but an incorrect CLI trace filename suffix prevented stage logs from being written. They were retained as preliminary receipts. Their all-U continuation was deliberately cancelled after 73.83 seconds to repair logging; that cancellation is not a production failure. Fresh traced runs use `-t trace`, which writes `state/logs/gammalog-trace.jsonl`. Full-orientation probes have explicit wall/RSS bounds; a probe stopped at a diagnostic bound is not an acceptance verdict.

All diagnostic generation processes have stopped. Successful generation controls, deliberate resource stops, and preparation cancellations are kept distinct; none of these diagnostic probes is added to the original acceptance-test pass count. No new numerical UV/IR correctness claim follows from these size measurements. The original run checked its report links and raw artifacts; those artifacts are no longer distributed; the jj change contains only documentation/evidence and has no conflicted ancestors. The source/executable manifest records the exact tested binary.
