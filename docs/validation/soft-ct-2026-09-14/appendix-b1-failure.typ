= Appendix B.1 acceptance failure diagnosis

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<appendix-b1-acceptance-failure-diagnosis>
The Rust test failure is a confirmed #strong[projection-provenance contract mismatch at the empty forest root];. The Rust acceptance test remains failed; its full-forest UV checks, matched bare control, and both later routes were not executed inside that test. Subsequent independent CLI diagnostics found a separate actual projected-local4D numerical failure: 13 of 21 summed UV fits fail, while both direct3D routes pass. See the supplementary numerical report;. The provenance mismatch must not be used to dismiss that later physical failure.

== Exact failure and first differing boundary
<exact-failure-and-first-differing-boundary>
The only printed route is `explicit_orientation_sum_only=false, project_local_4d=false`. It failed at #link("../../../tests/tests/uv.rs:6028")[uv.rs:6028];, `assert!(!paths.is_empty(), "Appendix-B.1 must retain at least the identity projection path")`, after 72.271 s according to nextest. The assertion expects an array containing at least one identity-path record (`steps=[]`) even for the root.

The complete pipeline is: import the Appendix graph/model and runtime settings → classify the UV wood → compute the same signed forest from the stored production CFF → assemble the local integrand/evaluator → independently export the same computed forest → serialize each node\'s local provenance → validate JSON provenance. Generation and export both succeeded. The first mismatching boundary is the provenance array: #link("../../../crates/gammalooprs/src/uv/approx/direct_3d/forest.rs:221")[Direct3dCts::projection\_paths] returns `Vec::new()` for `Root(_)`; the acceptance oracle requires a nonempty vector.

This identifies the failing node without an additional run: #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:879")[compatible\_topological\_order] explicitly moves/inserts the root first; #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:1420")[export\_node\_expressions] preserves that order, and the export collector preserves expression order. The root\'s stored production CFF is held independently in `Direct3dCts::Root`; its numerical `branches()` returns those branches. Omitting projection-history records cannot remove that root integrand. No numerator expansion, contour analysis or low-level CFF diagnosis is needed to explain this failure.

The projected-local4D route intentionally exposes an empty path list through `Local3DCts::projection_paths()` and the acceptance test skips path checks for that route. Explicit direct3D shares the root representation and would reach the same root-path assertion if earlier steps succeeded. Thus a route matrix follows directly from source: direct3D localized = fails observed; direct3D explicit = predicted same provenance failure; projected local4D = excludes this particular path assertion, otherwise untested here. This is not a numerical A/B disagreement between localization and explicit sums; selector algebra has not been implicated.

History check: `jj file show -r yssunply tests/tests/uv.rs` already contains the nonempty identity-path assertion at old lines 5388--5397, while the same revision owns `Root(_) => Vec::new()`. The contradiction therefore predates the later `ywptrwuz` adaptation. Current `jj annotate` attributes the assertion lines to `ywptrwuz` because its enclosing guard/indentation changed; that is not evidence of first introduction. The source establishes an existing stale/inconsistent acceptance contract, not a newly introduced rebase regression. No historical checkout was executed. Which public provenance contract to preserve still requires maintainer intent; there were no test or production edits.

There is a related deliberate rule for non-root zero sectors: they remain numerically stored but do not publish projection histories. The existing `zero_contributions_do_not_publish_projection_paths` unit test passed in this campaign (0.059 s). It should be considered when repairing the acceptance oracle; simply requiring every non-root term to publish a path can also conflict with zero-sector semantics. The present failure occurs at the root before this issue could arise.

== Checks actually completed
<checks-actually-completed>
+ Localized direct3D full-forest generation produced a 6-dimensional local integrand. Classification passed: exactly three connected components, outer H2 edges \[2,3,4,5,6\], d=2, identifier external \[-21,21\]/internal \[6,21\]; gamma1 U1 edges \[3,4\], d=1; gamma2 U0 edges \[2,3,5,6\], d=0.
+ The child-only control disabled outer H2 and gamma2, generated successfully, retained exactly gamma1 U1, and matched the full graph\'s canonical LMB route.
+ The child-only CLI UV profile completed 25 scales from 1e8 to 1e12 with seed 1337, all LMBs/subsets, automatic precision retries and per-orientation coverage. The saved JSON contains one graph, 7 LMBs, 21 summed rays and 294 orientation rays (14 distinct orientations per ray). Every summed and orientation fit resolved.
+ The six d=1 child-hard rays are exactly free e3/e4 crossed with fixed e2/e5/e6. All six summed slopes satisfy the explicit `< -0.9` assertion:

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,right,),
    table.header([Free], [Fixed], [Slope], [R²],),
    table.hline(),
    [e4], [e6], [-1.000000001195171], [1.0],
    [e3], [e6], [-1.000000001784386], [1.0],
    [e4], [e5], [-0.999999998525664], [1.0],
    [e3], [e5], [-1.000000000202399], [1.0],
    [e4], [e2], [-0.999999998976868], [1.0],
    [e3], [e2], [-1.000000000315168], [1.0],
  )]
  , kind: table
  )

#block[
#set enum(numbering: "1.", start: 5)
+ Structure-only forest export contained no constructed node atoms and exactly six nodes: root, three singleton operations and the two allowed child-to-outer chains. Computed export succeeded, had unique term identities, and retained six unique structural node keys. Its generic helper also parsed provenance for every term, verified 3D representation, nondecreasing child-to-parent topology order for published paths, consistent replay routes, and valid parent references.
+ The root\'s provenance parsed, had empty parent keys, native 3D representation and no legacy 4D component/branch payload. The next nonempty-path assertion failed.
]

The child-only profile is intentionally not required to remove the outer divergences. Reapplying the production verdict order (missing → R²\<0.99 → slope\>-0.9) to its saved fits gives 6 passed / 4 unstable / 11 DOD failures among summed rays and 140 passed / 110 unstable / 44 DOD failures among orientation rays. These remaining failures are expected for that partial control and must not be misreported as the complete H2 forest\'s numerical result. The acceptance test requires the six isolated gamma1-hard summed slopes to pass, and they did.

== Checks not reached by the original Rust test
<checks-not-reached-by-the-original-rust-test>
- The additional Appendix-specific root identity-path step assertion, all later Appendix-specific non-root checks, nonempty canonical/actual LMB fields, strictly increasing nested child-to-parent ordering, required nested-path existence, and positive-degree H2 conceptual branches U+S−US. The generic export helper had already checked nondecreasing published-path order, route consistency and parent references.
- Full H2-forest UV profile and its zero-failure assertion, including all localized orientations.
- Full/child profile identity comparison.
- Bare generation, bare UV nonvacuity/DOD-failure control, and full/bare profile identity comparison.
- Every generation, structural/provenance check and physical profile in the explicit-direct3D and projected-local4D loop iterations.

No saved full-forest `local_ct_uv_profile/uv_profile.json` or bare state exists in this failed acceptance directory. The computed forest export was retained in memory, not written as a standalone export artifact.

== Generation observations before failure
<generation-observations-before-failure>
The full localized forest summary shows expression build 5.68 s, Spenso evaluator construction 52.30 s, Symbolica evaluator construction 4.86 s, compile 0 ms, one evaluator and one core: approximately 62.84 s total displayed stages, sampled process RAM 1.18 GiB. Child-only generation shows 160 ms + 71 ms + 40 ms = approximately 271 ms, one evaluator/core, sampled RAM 1.08 GiB. The second RAM sample includes allocations retained from the first generation; it is not intrinsic child-only memory. The additional test time includes imports, the child profile and computed-forest replay.

== Supplementary numerical diagnostics
<supplementary-numerical-diagnostics>
Separately reported CLI diagnostics filled the missing full/bare numerical matrix without changing or bypassing the Rust acceptance assertion. Both direct3D full forests passed, while projected local4D failed 13/21 summed fits. All bare controls failed as expected. The original Rust test remains failed. The six diagnostics used fresh output states (three routes × full/bare), the original card and exact in-memory route/prescription overrides, then `run generate; save state; profile ultra-violet ...`. Save preceded the profile because a failing numerical verdict returns a nonzero CLI status after writing its report.

Every diagnostic uses `--min-scaling 8 --max-scaling 12 --n-points 25 --seed 1337 --selected-limits all`, with `--per-orientation` only for localized direct3D. Directions/norms remain unspecified so the profiler uses the same seed-driven rays as the acceptance helper. Preserve SingleParametric, MT=173, WT=0, Q=300, mUV=muR=20, disabled generated thresholds, disabled iterative orientation optimization and unintegrated counterterms. Use source defaults for precision, including automatic arbitrary-precision retry. Full profiles require zero failures and resolved fits; bare profiles require resolved DOD failures. Compare full/bare LMB/subset/orientation identities within each route and summed identities across routes.

The six supplementary runs completed at 10:08:28 UTC after benchmark-state preparation, alongside the long Figure B1 validation. Their generation wall times include shared-host contention and are not isolated throughput measurements. Later explicit/projected U1 child-only diagnostics pass all six intended child-hard rays in each route; see the #link("appendix-child-diagnostic.typ")[isolated-child report];. That numerical scope reduction does not certify exact reconstruction or identify a particular remaining forest operator as faulty.

Artifacts: failure log;, child profile;, generation log;.
