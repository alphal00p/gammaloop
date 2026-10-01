= Completed regular numerical soft-CT validation

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<completed-regular-numerical-soft-ct-validation>
Extracted from 34 completed regular-profile JSON files under /tmp/soft-ct-validation-2026-09-14, for input revision efaf2da158c5b01bddbc5435592e7cfea7a45a91. No generators, profiles, builds or tests were run during extraction.

#strong[All 1,217 subtracted/MUV UV checks pass the implemented threshold:] 228 summed-limit records and 989 orientation-limit records. These include the separate single-orientation DGSE check and the ordinary-MUV disconnected control. All records have fitted slopes; #strong[zero results rely on certified-vanishing missing-fit exemptions];. The matched bare profiles deliberately contain failing UV limits, establishing non-vacuous controls.

The three routes are direct 3D with per-orientation analysis, direct 3D with an explicit orientation sum, and projected local 4D with an explicit orientation sum. The summed routes intentionally contain no individual orientation records; --- below does not denote a missing requested check.

Machine-readable extraction, including every limit and raw fit metadata;. Extraction script;.

== Acceptance interpretation
<acceptance-interpretation>
The test calls pass\_fail(-0.9): a fitted UV exponent must be at most -0.9 and finite R² must be at least 0.99. Resolved means a fit is present, exactly as counted by the Rust helper; a bare fit can be resolved yet fail for a slope above threshold or poor R². The schema nests summed fits at graphs\[\].lmbs\[\].subsets\[\].analysis.inspect\_level and individual fits at per\_orientation\_inspect\_entries\[\].analysis.

Vanishing is counted separately only when a fit is absent, allow\_vanishing\_missing\_fits=true, and the finite/positive-finite sample counts satisfy the production criterion. None of these regular profiles permits the exemption and none lacks a fit.

UV exponents are reported asymptotic scaling, not execution-time measurements.

== Subtracted UV results
<subtracted-uv-results>
Each resolved/total pair counts LMB/subset records, not distinct graph topologies. Slope ranges below are summed fits across selected LMBs.

#figure(
  align(center)[#table(
    columns: 7,
    align: (auto,auto,right,right,auto,right,right,),
    table.header([Fixture], [Route], [Summed resolved/total], [Orientation resolved/total], [Summed slope range], [Minimum summed R²], [UV failures],),
    table.hline(),
    [Consecutive scalar bubbles], [D3 per orientation;], [12/12], [48/48], [-4.000000000 to -2.000000000], [1.000000000], [0],
    [Consecutive scalar bubbles], [D3 summed;], [12/12], [---], [-4.000000000 to -2.000000000], [1.000000000], [0],
    [Consecutive scalar bubbles], [P4 summed;], [12/12], [---], [-4.000000000 to -2.000000000], [1.000000000], [0],
    [DGSE, one selected orientation], [D3 selected orientation;], [27/27], [27/27], [-3.000000088 to -0.999997565], [0.999999935], [0],
    [DGSE], [D3 per orientation;], [27/27], [810/810], [-2.000020527 to -0.999811725], [0.999999816], [0],
    [DGSE], [D3 summed;], [27/27], [---], [-2.000043970 to -0.999811792], [0.999999913], [0],
    [DGSE], [P4 summed;], [27/27], [---], [-2.000000024 to -0.998761370], [0.999970985], [0],
    [Disconnected scalar spectacles], [D3 per orientation;], [12/12], [48/48], [-4.000000000 to -2.000000000], [1.000000000], [0],
    [Disconnected scalar spectacles], [D3 summed;], [12/12], [---], [-4.000000000 to -2.000000000], [1.000000000], [0],
    [Disconnected scalar spectacles], [P4 summed;], [12/12], [---], [-4.000000000 to -2.000000000], [1.000000000], [0],
    [Massless gluon self-energy], [D3 per orientation;], [2/2], [4/4], [-1.000010844 to -1.000007259], [1.000000000], [0],
    [Massless gluon self-energy], [D3 summed;], [2/2], [---], [-1.000010844 to -1.000007259], [1.000000000], [0],
    [Massless gluon self-energy], [P4 summed;], [2/2], [---], [-1.000010844 to -1.000007259], [1.000000000], [0],
    [Massless quark self-energy], [D3 per orientation;], [2/2], [4/4], [-1.027403362 to -0.935352940], [0.992963692], [0],
    [Massless quark self-energy], [D3 summed;], [2/2], [---], [-1.027403362 to -0.935352940], [0.992963692], [0],
    [Massless quark self-energy], [P4 summed;], [2/2], [---], [-1.027403362 to -0.935352940], [0.992963692], [0],
  )]
  , kind: table
  )

DGSE resolves 27 summed rays in each route. Its direct per-orientation run covers 30 distinct orientation labels on every ray, giving 810 fitted orientation/ray records with zero failures. Individual slopes span -4.000000760 to -0.976733859, with minimum R² 0.993293954. The generation-filter check resolves precisely one orientation across 27 rays (27/27), also with zero failures.

Massless quark self-energy passes both LMBs. Its weakest summed slope is -0.935352940 with R² 0.992963692, inside the implemented bound of -0.9. Direct per-orientation quark minimum R² is 0.995936344. Gluon summed slopes are approximately -1.00001 with R² effectively one.

Disconnected and consecutive scalar fixtures give summed slopes approximately -2 when one loop is hard and -4 when both are hard. Each resolves 12 LMB/subset rays per route and 48 individual orientation/ray fits in the direct per-orientation route. Individual exponents are generally near -1 or -2; stronger summed exponents include cancellations across orientations.

== Matched bare controls
<matched-bare-controls>
These failures are expected control outcomes, not test failures. Controls use the same settings and LMB/subset/orientation inventory as their subtracted fixture.

#figure(
  align(center)[#table(
    columns: 6,
    align: (auto,auto,right,right,auto,auto,),
    table.header([Fixture], [Route], [Bare summed failures / total], [Bare orientation failures / total], [Failure reasons: summed; orientation], [Bare summed slope range],),
    table.hline(),
    [Consecutive scalar bubbles], [D3 per orientation;], [12/12], [48/48], [12 dod\_exceeds\_threshold; 48 dod\_exceeds\_threshold], [2.000000 to 4.000000],
    [Consecutive scalar bubbles], [D3 summed;], [12/12], [---], [12 dod\_exceeds\_threshold; not requested], [2.000000 to 4.000000],
    [Consecutive scalar bubbles], [P4 summed;], [12/12], [---], [12 dod\_exceeds\_threshold; not requested], [2.000000 to 4.000000],
    [DGSE], [D3 per orientation;], [17/27], [494/810], [8 dod\_exceeds\_threshold, 9 unstable\_fit; 224 dod\_exceeds\_threshold, 270 unstable\_fit], [-3.000000 to 1.000000],
    [DGSE], [D3 summed;], [17/27], [---], [8 dod\_exceeds\_threshold, 9 unstable\_fit; not requested], [-3.000000 to 1.000000],
    [DGSE], [P4 summed;], [17/27], [---], [8 dod\_exceeds\_threshold, 9 unstable\_fit; not requested], [-3.000000 to 1.000000],
    [Disconnected scalar spectacles], [D3 per orientation;], [12/12], [48/48], [12 dod\_exceeds\_threshold; 48 dod\_exceeds\_threshold], [2.000000 to 4.000000],
    [Disconnected scalar spectacles], [D3 summed;], [12/12], [---], [12 dod\_exceeds\_threshold; not requested], [2.000000 to 4.000000],
    [Disconnected scalar spectacles], [P4 summed;], [12/12], [---], [12 dod\_exceeds\_threshold; not requested], [2.000000 to 4.000000],
    [Massless gluon self-energy], [D3 per orientation;], [2/2], [4/4], [2 dod\_exceeds\_threshold; 4 dod\_exceeds\_threshold], [2.000000 to 2.000000],
    [Massless gluon self-energy], [D3 summed;], [2/2], [---], [2 dod\_exceeds\_threshold; not requested], [2.000000 to 2.000000],
    [Massless gluon self-energy], [P4 summed;], [2/2], [---], [2 dod\_exceeds\_threshold; not requested], [2.000000 to 2.000000],
    [Massless quark self-energy], [D3 per orientation;], [2/2], [4/4], [2 dod\_exceeds\_threshold; 4 dod\_exceeds\_threshold], [1.000000 to 1.000000],
    [Massless quark self-energy], [D3 summed;], [2/2], [---], [2 dod\_exceeds\_threshold; not requested], [1.000000 to 1.000000],
    [Massless quark self-energy], [P4 summed;], [2/2], [---], [2 dod\_exceeds\_threshold; not requested], [1.000000 to 1.000000],
  )]
  , kind: table
  )

Each DGSE bare summed run has 17 rejected fits: eight exceed the UV slope bound and nine have unstable fits. The individual-orientation bare run rejects 494/810: 224 exceed the slope bound and 270 have unstable fits. Calling all 494 records divergent slopes would be inaccurate. The eight/224 resolved slope failures establish the required bare controls independently of unstable fits.

Bare quark and gluon summed slopes are approximately +1 and +2. Scalar bare controls have slopes approximately +2 (one hard loop) and +4 (both hard), with every selected summed and individual ray exceeding the bound.

== Ordinary-MUV disconnected control
<ordinary-muv-disconnected-control>
Ordinary-MUV disconnected spectacles also passes 12/12 summed rays in each route and 48/48 individual orientation/ray records in the direct per-orientation route. Summed slopes range from approximately -4 to -2. This checks preservation of established ordinary UV behavior; these runs alone do not demonstrate a soft-specific gain.

- D3 per orientation ordinary-MUV profile
- D3 summed ordinary-MUV profile
- P4 summed ordinary-MUV profile

== DGSE soft limit: exact fitted values were not persisted
<dgse-soft-limit-exact-fitted-values-were-not-persisted>
The regular DGSE test runs matched S(e0) IR profiles with seed 1337 in all three routes, but leaves InfraRedProfile.output\_file=None. The handler returns fitted scalings in memory without writing soft\_profile.json. No standalone IR-fit JSON or retained scaling record was found in the DGSE state directories. Exact fitted soft slopes and their measured difference cannot be recovered from UV JSON and are not fabricated here.

The passing test certifies, on the identical recorded routed ray:

- Both fits are finite and have R² at least 0.9.
- Subtracted reported IR scaling is greater than 3.0, the fixture\'s bound for pointwise boundedness.
- Bare reported IR scaling is at most 3.0.
- Subtracted minus bare reported scaling is at least 1.0.

These are assertions verified by the run, not substituted measured values. Exact slopes would require profiling the preserved states with an explicit IR output file. Relevant code: #link("../../../tests/tests/uv.rs:4030")[test];, #link("../../../crates/gammaloop-api/src/commands/profile.rs:504")[optional IR output];.

== Precision and scope
<precision-and-scope>
The extraction JSON retains arbitrary-precision retry flags per fit. DGSE\'s subtracted per-orientation route retries eight of 27 summed fits and 16 of 810 individual fits. Every summed scalar disconnected/consecutive ray retries, as do quark and gluon summed fits. These are successful precision recoveries.

This report covers completed regular profile artifacts only. Mandatory top-bubble and Appendix acceptance, integration checks that emit no profile JSON, and runtime benchmarks belong to the main validation report.
