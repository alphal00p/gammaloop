= Appendix B.1: two ordinary children, outer H2 disabled

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<appendix-b1-two-ordinary-children-outer-h2-disabled>
UV slopes here describe the measure-weighted profile quantity `lambda^(3*n_hard) * |I(lambda)|`, not the raw integrand alone. The profiler applies this prefactor at #link("../../../crates/gammalooprs/src/uv/profile.rs:2854")[profile.rs:2854];. For these one-hard rays, a reported slope +2 corresponds to a raw integrand exponent about -1; a reported slope -1 corresponds to a raw exponent about -4. The acceptance threshold acts on the reported measure-weighted slope.

Restoring gamma2 U0 while keeping the outer H2 disabled gives a measure-adjusted UV degree of divergence of approximately +2 on all six gamma1-hard rays in #strong[both] explicit 3D and projected 4D. The UV profiler weights the integrand by `λ^(3 n_hard)` before fitting (#link("../../../crates/gammalooprs/src/uv/profile.rs:2854")[source];); these one-hard-loop fits therefore correspond to raw integrand exponents near −1. Their complete target-fit records, including sampled values, are identical across the routes. This partial forest therefore does not preserve the isolated U1 control's child-hard finiteness. It does not reproduce a projected-only discrepancy.

With the existing full-forest results, the numerical boundary is now: isolated U1 removes these subdivergences in both routes; adding gamma2 gives the same +2 measure-adjusted degree in both routes; adding the outer-containing forest contributions cancels that growth successfully in explicit 3D, while the full projected route still fails. The complete forest is where the route-dependent acceptance failure is observed in these comparisons. This does not locate a defect solely inside the outer H2 term: a smaller earlier difference could be hidden by the shared leading divergence, and agreement of the sampled partial-forest slopes/values does not prove equal coefficients or full functions. It is not an exact symbolic diagnosis or evidence that the incomplete two-child forest ought to pass these limits.

The controls use the original Appendix run card, unchanged MUV defaults, and only one override: `external=[-21,21], internal=[6,21]` is Unsubtracted. Relative to the isolated-U1 control, the sole physical change is removal of the `internal=[6]` Unsubtracted override, restoring gamma2. Source inventory identifies gamma1 as edges \[3,4\], U1; gamma2 as \[2,3,5,6\], U0; and outer H2 as \[2,3,4,5,6\]: #link("../../../tests/tests/uv.rs:5779")[uv.rs:5779];. Each saved structure-only DOT contains the root and two separate child nodes, with no outer node. The children overlap at edge e3 and neither edge set contains the other; they are neither disjoint nor nested. There is therefore no joint product branch for the two children. This agrees with the production factorized-join requirement of pairwise disjoint connected components: #link("../../../crates/gammalooprs/src/uv/hedge_poset.rs:101")[hedge\_poset.rs:101];. Neither graph numerators nor computed forest expressions were exported.

Both fresh states use complete orientation sums, local ThreeD integrands, no integrated counterterms, and the same seed 1337 with 25 scales from 1e8 through 1e12. Saved generation settings differ from the matching one-child/full references only in the intended override array and state-folder path. Default/runtime integrand settings, model and model parameters match. The complete graph/LMB/subset/scaling profile identities match across both routes and both reference variants. Commands ran in isolated /tmp working directories with explicit settings paths, explicit `save state -p` destinations, and autosave disabled.

#figure(
  align(center)[#table(
    columns: 6,
    align: (auto,auto,right,right,right,right,),
    table.header([Free child edge], [Fixed outer edge], [Two-child explicit slope], [Two-child projected slope], [Full explicit slope], [Full projected slope],),
    table.hline(),
    [e3], [e2], [2.000000000157], [2.000000000157], [-1.000000000351], [1.025412539898],
    [e3], [e5], [2.000000000164], [2.000000000164], [-1.000000001214], [1.060065637217],
    [e3], [e6], [1.999999999942], [1.999999999942], [-1.000001111766], [1.062121423459],
    [e4], [e2], [2.000000000000], [2.000000000000], [-1.000000669219], [1.032727845462],
    [e4], [e5], [2.000000000000], [2.000000000000], [-0.999999673340], [1.053814450506],
    [e4], [e6], [2.000000000000], [2.000000000000], [-0.999999997664], [1.066173378682],
  )]
  , kind: table
  )

All twelve two-child target fits have R²=1, start at sample index 0, and use arbitrary-precision retry. The original child-control slope bound `< -0.9` is not met by any of them. These are resolved growth measurements, not missing or unstable target fits. Exact values and arrays are preserved in the JSON artifact; their equality is equality of these recorded numerical observations, not a symbolic expression certificate.

Each CLI command exits 1 after writing all 21 fitted limits. Across the full two-child profile there are 17 slope failures and 4 unstable fits, with zero passing limits. Seven all-hard limits lie outside the intended target-ray comparison and remain unsubtracted without the outer component; the other non-target failures are likewise not evidence of a projected-only effect. The six target failures are reported explicitly and are not relabelled as expected passes.

#figure(
  align(center)[#table(
    columns: 5,
    align: (auto,right,right,auto,auto,),
    table.header([Route], [Command wall], [Peak child RSS], [Target rays], [Whole-profile verdict],),
    table.hline(),
    [explicit], [1.499963 s], [39.17 MiB], [6/6 resolved slope failures], [17 slope failures, 4 unstable],
    [projected], [1.666724 s], [57.17 MiB], [6/6 resolved slope failures], [17 slope failures, 4 unstable],
  )]
  , kind: table
  )

These command costs include startup, generation, saving, structure export and profiling while the large acceptance test remained active. They are not isolated throughput or speedup measurements. The completion gate was recreated after both processes exited, regardless of numerical verdict.

No low-level CFF/residue investigation followed this reduction. Any further projected reconstruction discrepancy analysis must retain raw post-Taylor factors and owner/provenance data, then certify exact numerator and denominator reconstruction before attributing it to CFF recursion or residue aggregation. The existing original Rust test outcomes remain unchanged.

Evidence: full diagnostic data and input checks;, exact plan;, runner;.

- explicit: profile;, command/resources;, log;, structure-only forest;.
- projected: profile;, command/resources;, log;, structure-only forest;.
