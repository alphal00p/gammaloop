= Appendix B.1: isolated U1 child control

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<appendix-b1-isolated-u1-child-control>
UV slopes here describe the measure-weighted profile quantity `lambda^(3*n_hard) * |I(lambda)|`, not the raw integrand alone. The profiler applies this prefactor at #link("../../../crates/gammalooprs/src/uv/profile.rs:2854")[profile.rs:2854];. For these one-hard rays, a reported slope +2 corresponds to a raw integrand exponent about -1; a reported slope -1 corresponds to a raw exponent about -4. The acceptance threshold acts on the reported measure-weighted slope.

UV slopes in this report include the `λ^(3 n_hard)` measure factor (#link("../../../crates/gammalooprs/src/uv/profile.rs:2854")[source];). For these one-hard-loop rays, a slope near −1 corresponds to a raw integrand exponent near −4.

The six gamma1-hard rays pass in both explicit 3D and projected 4D when only the ordinary U1 child counterterm is retained. The full projected-forest failure is not reproduced by this isolated child control at these rays. This is a numerical scope reduction; it does not certify the exact symbolic reconstruction or identify whether gamma2, the outer H2 terms, their composition, or numerical cancellation causes the full-forest discrepancy.

The diagnostic uses the original run card and the existing control at #link("../../../tests/tests/uv.rs:5797")[uv.rs:5797];. It changes the outer rule `external=[-21,21], internal=[6,21]` from IR to Unsubtracted and adds `external=[-21,21], internal=[6]` as Unsubtracted. The default MUV prescriptions remain. This is the source-prescribed `(1-U_child)I` control for child edges \[3,4\], degree 1. Each structure-only forest export contains exactly the root and one child node. No computed numerator export was requested.

Both routes use the complete orientation sum with an empty orientation filter, ThreeD final integrands, and integrated counterterms disabled. The standard CLI UV profile uses 25 logarithmically spaced scales from 1e8 to 1e12, seed 1337, and all limits. Parsed saved settings match the corresponding original full-route state except for the intended prescription overrides and state-folder path. Runtime settings, model, model parameters, and complete profile identities match. Across the two child routes the complete profile identities also match.

#figure(
  align(center)[#table(
    columns: 6,
    align: (auto,auto,right,right,right,right,),
    table.header([Free child edge], [Fixed outer edge], [Explicit child slope], [Projected child slope], [Explicit full slope], [Projected full slope],),
    table.hline(),
    [e3], [e2], [-0.999999978896], [-1.000000000315], [-1.000000000351], [1.025412539898],
    [e3], [e5], [-1.000000003945], [-1.000000004499], [-1.000000001214], [1.060065637217],
    [e3], [e6], [-0.999999969188], [-0.999999972317], [-1.000001111766], [1.062121423459],
    [e4], [e2], [-1.000000009441], [-0.999999998977], [-1.000000669219], [1.032727845462],
    [e4], [e5], [-0.999999998526], [-0.999999998526], [-0.999999673340], [1.053814450506],
    [e4], [e6], [-0.999999999336], [-1.000000001195], [-0.999999997664], [1.066173378682],
  )]
  , kind: table
  )

The source control requires exactly these six rays and slope \< -0.9: #link("../../../tests/tests/uv.rs:5910")[uv.rs:5910];. All twelve child-only fits also satisfy the standard R² \>= 0.99 stability threshold; the minimum R² is 0.9999999999999941. Explicit 3D uses arbitrary-precision retries for 1/6 child fits; projected 4D for 4/6. The explicit e4-hard/e6-fixed fit selects the tail beginning at point index 4 after an early outlier (its full-range R² is about 0.9584); the other eleven selected fits start at point index 0. Exact slopes, fit windows, full-range statistics and all sampled points are preserved in the JSON report.

Both CLI commands return exit code 1 because the remaining fifteen limits fail the whole-profile acceptance criterion. This is expected for the child-only diagnostic: gamma2 and overall UV subtraction were deliberately removed. Each profile resolves all 21 limits, with the six intended child-hard limits passing; these exit codes must not be counted as failures of the six-ray control or as successful full UV subtraction. The literal CLI summary `FAIL (6/21)` means six passing limits out of 21 in this output; the detailed failure list contains fifteen entries.

#figure(
  align(center)[#table(
    columns: 5,
    align: (auto,right,right,auto,auto,),
    table.header([Route], [Runner wall], [Peak child RSS], [Six-ray result], [Other 15 limits],),
    table.hline(),
    [explicit], [0.980248 s], [30.09 MiB], [6/6 pass], [4 unstable\_fit, 11 dod\_exceeds\_threshold],
    [projected], [0.964665 s], [39.20 MiB], [6/6 pass], [4 unstable\_fit, 11 dod\_exceeds\_threshold],
  )]
  , kind: table
  )

These command walls include startup, generation, save, structural export and profiling while the large Figure B1 test remained active. They are single observations, not isolated throughput or historical speedup measurements.

Current discrepancy boundary: the original full explicit route passes all six rays, the full projected route fails them, and both isolated U1 routes pass them with matched inputs. This excludes a reproduced failure of the isolated U1 child on the tested rays; it does not prove U1 correctness generally. The next symbolic comparison must retain raw post-Taylor sectors, original edge ownership/provenance and factorized numerators, and establish exact numerator/denominator reconstruction certificates before investigating low-level CFF recursion or residue aggregation. No such further investigation or production edit was performed here.

Evidence: machine-readable result and checks;, exact CLI plan;, runner;.

- explicit: profile;, log;, command/resources;, structure-only forest;.
- projected: profile;, log;, command/resources;, structure-only forest;.

These are supplementary CLI diagnostics; the original Rust acceptance and its assertions remain unchanged, with their original outcome preserved.
