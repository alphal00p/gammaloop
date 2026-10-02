= Appendix B.1 supplementary numerical validation

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<appendix-b1-supplementary-numerical-validation>
#strong[Projected local4D fails the complete Appendix B.1 UV profile: 13 of 21 resolved summed limits fail. Both direct3D routes pass all 21 summed limits, and localized direct3D also passes all 294 orientation limits.] This is an actual route-dependent numerical failure, separate from the stale root-provenance assertion that stopped the Rust acceptance test before its full-forest physical profile. The acceptance test is still failed, and three-route numerical acceptance cannot be claimed.

All six supplementary CLI cases completed with written profiles and generated states. They ran sequentially after benchmark-state preparation, while the long Figure B1 test continued. No test or production source changed. The final throughput measurements ran after these supplementary runs finished.

== Result matrix and controls
<result-matrix-and-controls>
Every case covers one graph, seven LMBs and all 21 subsets at 25 scales from 1e8 to 1e12, seed 1337. Localized cases additionally cover fourteen orientations per subset (294 orientation limits). All returned fits resolve; there are no missing-fit failures. A fit passes if slope ≤ −0.9 and R² ≥ 0.99. “Bare” deliberately disables all subtractions and is expected to return a failed CLI numerical verdict.

#figure(
  align(center)[#table(
    columns: 9,
    align: (auto,auto,right,right,right,right,right,right,right,),
    table.header([Route], [Forest], [Summed resolved], [Summed pass], [Summed DOD fail], [Summed unstable], [Orientation resolved], [Orientation pass / DOD fail / unstable], [CLI exit],),
    table.hline(),
    [Localized direct3D], [full], [21/21], [21], [0], [0], [294/294], [294 / 0 / 0], [0],
    [Localized direct3D], [bare], [21/21], [4], [13], [4], [294/294], [116 / 100 / 78], [1],
    [Explicit direct3D], [full], [21/21], [21], [0], [0], [0/0], [0 / 0 / 0], [0],
    [Explicit direct3D], [bare], [21/21], [4], [13], [4], [0/0], [0 / 0 / 0], [1],
    [Projected local4D], [full], [21/21], [8], [6], [7], [0/0], [0 / 0 / 0], [1],
    [Projected local4D], [bare], [21/21], [4], [13], [4], [0/0], [0 / 0 / 0], [1],
  )]
  , kind: table
  )

All full/bare scale, graph, LMB, subset, initial-DOD and orientation-label identities match within each route. After removing orientation labels, all three route identities match for each forest. The original successful child-only profile also has exactly the same profile identity as the supplementary localized-full profile, restoring that previously unreached comparison. Each bare control has 13 summed DOD failures; localized bare additionally has 100 orientation DOD failures. Their expected control verdicts therefore pass. Bare failures are not new subtraction regressions.

For the full forest, localized summed slopes span −2.00000000007 to −0.999999997664 (minimum R²=1), and explicit summed slopes span −2.00000000007 to −0.999999673340 (minimum R²=0.999999999955). All 294 localized orientation fits pass; their worst slope is −0.999985336891 and minimum R² is 0.999999982930.

== All thirteen failed projected-local4D rays
<all-thirteen-failed-projected-local4d-rays>
The six d=1 rays are exactly the child-hard gamma1 rays: free e3/e4 and fixed e2/e5/e6. Their expected decay near −1 in both direct3D full forests becomes growth near +1 in projected local4D. The other seven failures are all-hard rays with poor fits, not missing fits. `fit_start` is the zero-based scale-sample index selected for the fitted suffix; full-range values remain in the JSON.

#figure(
  align(center)[#table(
    columns: 9,
    align: (right,auto,auto,right,right,right,right,auto,auto,),
    table.header([LMB], [Free edges], [Fixed edges], [Initial DOD], [Fitted slope], [R²], [Fit start], [Failure], [Precision],),
    table.hline(),
    [0], [e4], [e6], [1], [1.066173378682], [0.993304991849], [2], [dod\_exceeds\_threshold], [Arb retry],
    [0], [e4,e6], [none], [2], [0.021330714116], [0.265587531216], [0], [unstable\_fit], [Arb retry],
    [1], [e3], [e6], [1], [1.062121423459], [0.993670325982], [1], [dod\_exceeds\_threshold], [Arb retry],
    [1], [e3,e6], [none], [2], [-0.070124600344], [0.369586145495], [0], [unstable\_fit], [Arb retry],
    [2], [e4], [e5], [1], [1.053814450506], [0.991097018135], [0], [dod\_exceeds\_threshold], [Arb retry],
    [2], [e4,e5], [none], [2], [0.018591823112], [0.269096414752], [0], [unstable\_fit], [Arb retry],
    [3], [e3], [e5], [1], [1.060065637217], [0.993825694587], [1], [dod\_exceeds\_threshold], [Arb retry],
    [3], [e3,e5], [none], [2], [0.040102596508], [0.035965226455], [0], [unstable\_fit], [Arb retry],
    [4], [e3,e4], [none], [2], [-0.027668028427], [0.325187862607], [0], [unstable\_fit], [Arb retry],
    [5], [e4], [e2], [1], [1.032727845462], [0.997200553111], [1], [dod\_exceeds\_threshold], [Arb retry],
    [5], [e2,e4], [none], [2], [0.058571586454], [0.239353563478], [0], [unstable\_fit], [Arb retry],
    [6], [e3], [e2], [1], [1.025412539898], [0.998464164998], [2], [dod\_exceeds\_threshold], [Arb retry],
    [6], [e2,e3], [none], [2], [0.144878005584], [0.122698349306], [0], [unstable\_fit], [Arb retry],
  )]
  , kind: table
  )

All thirteen failures remain after the profiler\'s explicit arbitrary-precision retry (`used_arb_prec_retry=true`). The implementation\'s Arb type is `VarFloat<1000>` (#link("../../../crates/gammalooprs/src/utils/mod.rs:2640")[utils];). The profile retries a complete ray using `use_arb_prec=true` after an unacceptable initial fit (#link("../../../crates/gammalooprs/src/uv/profile.rs:2817")[profile];); that flag selects the Arb stability level (#link("../../../crates/gammalooprs/src/integrands/process/mod.rs:1589")[evaluation];). The JSON stores final fit values as f64, so it does not preserve full-precision complex samples. It records 20 of 21 projected summed rays retried, versus 21 of 21 localized summed rays and 13 of 21 explicit summed rays.

A concrete matched ray is free e4/fixed e6. The localized and explicit profile magnitudes are identical at the endpoints: 7.1649464317e−25 at scale 1e8 and 7.1649467015e−29 at 1e12. The projected magnitudes are 6.0682572410e−25 and 1.0966893086e−21. This exhibits the failing asymptotic behavior directly; it is not merely a marginal threshold crossing. These are numerical samples, not an algebraic certificate or identification of a particular faulty coefficient.

== Exact command and input equivalence
<exact-command-and-input-equivalence>
The supplementary cases use the original #link("../../../tests/resources/run_cards/paper_appendix_b1_nested_gluon_self_energy_ir.toml:1")[Appendix card];. The card imports the same SM graph, sets MT=173 and WT=0, Q=e\_cm=300, helicities \[1,−1\], mUV=muR=20, SingleParametric evaluation, hedge\_poset orchestration, disabled generated thresholds and disabled iterative orientation optimization. Its outer IR rule is external \[-21,21\], internal \[6,21\]; the other prescription defaults are MUV. Integrated counterterms are disabled.

The #link("../../../tests/tests/uv.rs:2510")[acceptance helper] sets `explicit_orientation_sum_only`, final integrand ThreeD and `local_uv_cts_from_expanded_4d_integrands`, requires unintegrated generation, and clears the orientation pattern for explicit sums. The original card\'s orientation pattern is already empty; every saved supplementary pattern is empty. Our CLI applies the same route values before running that same generation block. For bare controls it sets all three prescription defaults to Unsubtracted and clears `overrides=[]`, exactly as the helper does. There is no Compare-to-Hedge switch needed because this card already uses hedge\_poset.

The profile command follows #link("../../../tests/tests/uv.rs:3542")[cli\_uv\_profile\_pass\_fail];: original min/max exponents, point count and seed, appending `--selected-limits all`, and retaining `--per-orientation` only when explicit sums are disabled. Fixed UV directions/norms are unspecified, as in the Rust helper; both use the same seeded sampler. Neither sets `--use_f128`, analytic analysis, cut selection, graph filtering or fail-fast. The CLI writes a complete failed-limit report and then returns exit 1 on failure; the Rust helper directly invokes the profile handler to inspect that same report without turning the bare verdict into an early command error.

A concrete projected-full CLI command body is:

```text
set global kv global.generation.explicit_orientation_sum_only=true global.generation.uv.local_uv_cts_from_expanded_4d_integrands=true global.generation.uv.final_integrand=ThreeD global.generation.uv.generate_integrated=false; run generate; save state; profile ultra-violet -p paper_appendix_b1_nested_gluon_self_energy -i soft_ir --min-scaling 8 --max-scaling 12 --n-points 25 --seed 1337 --selected-limits all --output /tmp/soft-ct-validation-2026-09-14/appendix-supplementary-projected-full/profile
```

It is passed as an argv element to `/common/dev/gammaloop/lcnbr/target/dev-optim/gammaloop <original-card> -s <fresh-state> run -c <body>`. The complete six argv vectors are preserved in the plan;. State is explicitly saved before profiling so a numerical exit 1 still leaves an inspectable generated state.

Verification went beyond exit codes: all six have `amp.bin`, `integrand.bin` and actual `generation_summary.json`; saved global settings match route/final-dimension/unintegrated flags and full/bare prescriptions. All six generated `integrand/settings.toml` files are equal. All three full-forest `model_parameters.json` and runtime defaults are equal, and the original failed test\'s runtime defaults equal the supplementary defaults. The explicit-full and projected-full complete global settings are identical after changing only `local_uv_cts_from_expanded_4d_integrands`. Thus the strongest A/B is explicit direct3D versus projected local4D: both use an explicit full sum and the same numeric inputs, with only the construction route switched.

== First-boundary triage
<first-boundary-triage>
Shared stages are the graph/model/card, physical external state, UV prescription, unintegrated forest selection, imported canonical graph chart, explicit full-sum selector choice (for the strongest A/B), evaluator method, and seeded candidate-LMB/ray profile settings. Candidate ray momenta are passed through the same profiling evaluation entry point with the candidate LMB; source code explicitly casts to the requested precision before canonical routing. Full/bare profile identities and the explicit/projected runtime/settings comparisons rule out a changed ray or runtime selector as the explanation supported by these artifacts.

The first stage allowed to differ is local counterterm construction: direct3D applies the forest operations to the stored post-energy-integration production CFF; projected local4D constructs factorized four-dimensional Taylor sectors and reconstructs their projected CFF source before final forest/evaluator assembly. The completed evaluator expressions differ downstream. Therefore the shared numerical evaluator cannot be blamed without comparing its actual expression/parameter inputs, and no comparison here identifies CFF recursion, contour signs or residue aggregation as faulty. The provenance-root assertion only reads exported metadata and did not run in these CLI cases; it cannot cause this numerical discrepancy. Bare controls show the no-counterterm workload behaves consistently across route settings.

Projected logs contain 5 local-4D construction events, 5 local-projection completions, and 36 each of `source_reconstruction` and `exact_cff_assignment_selection`; the latter report a certified assignment. These success/timing records do #strong[not] supply a durable factor-preserving equality certificate comparing raw Taylor sectors, normalized reconstruction input and projected source coefficients for the failing forest terms. No such full equality has been established in this task. The #link("../../../CONTRIBUTING.typ#projected-local-4d-uv-reconstruction-stop-rule")[CONTRIBUTING reconstruction stop rule] therefore still applies before any deeper CFF investigation. Finite-precision arithmetic during construction, loss of cancellations, or a mathematical projection mismatch remain possibilities; none is yet proved.

Subsequent isolated-child diagnostics completed the proposed next comparison: both explicit direct3D and projected local4D pass all six gamma1-hard rays when only U1 is retained (12/12 selected fits, minimum R²=0.9999999999999941). Their other fifteen whole-profile limits fail as expected after gamma2/outer subtraction is removed; those whole-profile CLI exits are not failures of the six-ray control. The #link("appendix-child-diagnostic.typ")[isolated-child report] records matched settings and two-node structure-only forests. This narrows the discrepancy to behavior present in the complete forest but not reproduced by the isolated child on these rays. It does not certify U1 generally or distinguish gamma2, outer H2, their composition or numerical cancellation.

Further forced-Arb comparisons of saved explicit/projected full and bare states show agreement at the generic point and increasing full-forest disagreement on hard-scaled canonical points; all ten bare comparisons agree after f64 serialization. These generic rays differ from the original seeded candidate-LMB rays. The #link("appendix-point-comparisons.typ")[matched-point report] reports these numerical comparisons without claiming an error bound or exact reconstruction. The next symbolic comparison must preserve the complete raw Taylor sectors, denominator ownership and factorized numerators, and certify the source reconstruction at the first differing complete-forest operation before any lower-level CFF investigation. No such exact full-forest certificate has been established.

== Generation and elapsed time
<generation-and-elapsed-time>
These precise construction timings come from persisted `Duration` fields in generation\_summary.json, unlike rounded log table values. All six cases generate one evaluator, use one generation core and do no backend compilation. They ran as separate processes, so the cross-route retained-memory issue of the original multi-route Rust test does not apply. The monitor still samples process memory every 100 ms, can miss short spikes, and excludes child-process memory. Other validation was running on the host, so these single observations are not a controlled speed comparison.

#figure(
  align(center)[#table(
    columns: 8,
    align: (auto,auto,right,right,right,right,right,right,),
    table.header([Route], [Forest], [CLI wall s], [Generation s], [Expression s], [Spenso build s], [Symbolica build s], [Sampled process RAM GiB],),
    table.hline(),
    [Localized direct3D], [full], [85.818], [68.766098], [6.604145], [56.767360], [5.394594], [1.207],
    [Localized direct3D], [bare], [1.656], [0.098284], [0.057433], [0.029451], [0.011400], [0.040],
    [Explicit direct3D], [full], [67.838], [64.735452], [5.592514], [56.704930], [2.438008], [1.178],
    [Explicit direct3D], [bare], [0.740], [0.144885], [0.110499], [0.029629], [0.004757], [0.036],
    [Projected local4D], [full], [75.992], [72.650307], [7.597105], [63.896534], [1.156669], [6.246],
    [Projected local4D], [bare], [0.966], [0.147739], [0.089254], [0.050762], [0.007723], [0.039],
  )]
  , kind: table
  )

Spenso construction dominates all three full forests. Projected local4D used about 6.246 GiB sampled process memory versus about 1.2 GiB for the direct routes in these fresh processes, and its numerical profile failed. Generation timing or evaluator object count must not be interpreted as runtime throughput; the wall-time gap includes state serialization and profiling.

Artifacts: structured results with every fit and source path;, complete execution plan;, coordinator log;, projected full profile;, root-provenance failure diagnosis;. Original Rust acceptance status remains failed; supplementary direct3D physical coverage is positive and supplementary projected-local4D physical coverage is negative.
