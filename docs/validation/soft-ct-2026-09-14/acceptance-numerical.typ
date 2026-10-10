= Mandatory acceptance numerical measurements

This is a historical run report. Its findings and summary measurements are retained; raw logs, JSON receipts, and JUnit reports are not distributed with the source. Historical artifact names and local paths below identify the original run, not files available in this checkout.
<mandatory-acceptance-numerical-measurements>
Extracted directly from isolated profile JSON reports and Nextest JUnit. Summary measurements, SHA-256 hashes and absolute provenance paths are retained in `/tmp/soft-ct-acceptance-numerical.json`; every compact UV fit record is retained separately in `/tmp/soft-ct-acceptance-numerical.fits.json`.

UV slopes are measure-adjusted: the profiler fits abs(evaluation) × lambda^(3#emph[n\_hard), so its reported exponent equals the raw evaluation exponent plus 3];n\_hard (uv/profile.rs:2854,3146). UV resolved fits pass iff slope \<= -0.9 and finite R² \>= 0.99. Missing fits count as certified vanishing only when allowed and finite\_samples \> positive\_finite\_samples with fewer than two positive samples (production `uv/profile.rs:2218,3179`). All UV rows below are standard 25-point windows 10^8→10^12 unless recorded otherwise in JSON.

== Execution outcomes
<execution-outcomes>
#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,auto,right,),
    table.header([Test (`slow::`)], [Outcome], [Testcase wall seconds],),
    table.hline(),
    [gamma\_star\_ddbar\_top\_bubble\_child\_only\_has\_two\_power\_soft\_improvement], [passed], [57.02],
    [gamma\_star\_ddbar\_top\_bubble\_has\_two\_power\_soft\_improvement], [passed], [198.364],
    [local\_nested\_soft\_fixture\_generates], [passed], [19.888],
    [paper\_appendix\_b1\_ir\_specialization\_generates\_and\_profiles], [failed], [72.271],
    [paper\_figure\_b1\_double\_triangle\_is\_a\_non\_vacuous\_uv\_fixture], [passed], [58.997],
    [paper\_figure\_b1\_massless\_bubble\_uses\_a\_complete\_soft\_wood], [failed], [6326.706],
  )]
  , kind: table
  )

== UV measured coverage
<uv-measured-coverage>
Each count cell is total/resolved/vanishing/failed. Bare failures are expected controls. These are profile artifacts, not an assertion that a failed containing test completed later stages.

#figure(
  align(center)[#table(
    columns: 7,
    align: (auto,auto,auto,auto,auto,auto,right,),
    table.header([Fixture], [Route], [Summed T/R/V/F], [Orientation T/R/V/F], [Summed slope range], [Orientation slope range], [Minimum R² (all resolved)],),
    table.hline(),
    [gamma\_star\_ddbar\_top\_bubble\_bare], [localized], [27/27/0/17], [810/810/0/494], [-3 to 2], [-3.00000001 to 2], [0],
    [gamma\_star\_ddbar\_top\_bubble\_bare], [explicit], [27/27/0/17], [0/0/0/0], [-3 to 2], [---], [0],
    [gamma\_star\_ddbar\_top\_bubble\_bare], [projected], [27/27/0/17], [0/0/0/0], [-3 to 2], [---], [0],
    [gamma\_star\_ddbar\_top\_bubble\_soft\_ir], [localized], [27/27/0/0], [810/810/0/0], [-3 to -0.999999994], [-3.00000002 to -0.999999934], [1],
    [gamma\_star\_ddbar\_top\_bubble\_soft\_ir], [explicit], [27/27/0/0], [0/0/0/0], [-3 to -0.999999994], [---], [1],
    [gamma\_star\_ddbar\_top\_bubble\_soft\_ir], [projected], [27/27/0/0], [0/0/0/0], [-3 to -0.999999994], [---], [1],
    [gamma\_star\_ddbar\_top\_bubble\_uv\_only], [localized], [27/27/0/0], [810/810/0/0], [-3 to -0.999999994], [-3.00000002 to -0.999692749], [0.9999706233],
    [gamma\_star\_ddbar\_top\_bubble\_uv\_only], [explicit], [27/27/0/0], [0/0/0/0], [-3 to -0.999999994], [---], [1],
    [gamma\_star\_ddbar\_top\_bubble\_uv\_only], [projected], [27/27/0/0], [0/0/0/0], [-3 to -0.999999994], [---], [1],
    [local\_nested\_soft\_top\_self\_energy\_bare], [localized], [21/21/0/13], [294/294/0/170], [-2 to 1.00000003], [-3 to 1.00000001], [1],
    [local\_nested\_soft\_top\_self\_energy\_bare], [explicit], [21/21/0/13], [0/0/0/0], [-2 to 1.00000003], [---], [1],
    [local\_nested\_soft\_top\_self\_energy\_bare], [projected], [21/21/0/13], [0/0/0/0], [-2 to 1.00000003], [---], [1],
    [local\_nested\_soft\_top\_self\_energy], [localized], [21/21/0/0], [294/294/0/0], [-2 to -0.999999557], [-4.00203489 to -0.999123635], [0.992388669],
    [local\_nested\_soft\_top\_self\_energy], [explicit], [21/21/0/0], [0/0/0/0], [-2 to -0.999999559], [---], [1],
    [local\_nested\_soft\_top\_self\_energy], [projected], [21/21/0/0], [0/0/0/0], [-2.00000063 to -0.99999956], [---], [1],
    [paper\_appendix\_b1\_child\_only], [localized], [21/21/0/15], [294/294/0/154], [-1 to 2], [-3 to 2.00000001], [0],
    [paper\_figure\_b1\_double\_triangle\_bare], [localized], [24/24/0/24], [432/432/0/312], [-1.42714797e-09 to 2], [-4 to 2], [0],
    [paper\_figure\_b1\_double\_triangle\_bare], [explicit], [24/24/0/24], [0/0/0/0], [-1.42714797e-09 to 2], [---], [0],
    [paper\_figure\_b1\_double\_triangle\_bare], [projected], [24/24/0/24], [0/0/0/0], [-1.42714797e-09 to 2], [---], [0],
    [paper\_figure\_b1\_double\_triangle], [localized], [24/24/0/0], [432/432/0/0], [-2 to -0.999991766], [-2 to -0.993269383], [0.9933226396],
    [paper\_figure\_b1\_double\_triangle], [explicit], [24/24/0/0], [0/0/0/0], [-2 to -0.996224197], [---], [0.995552061],
    [paper\_figure\_b1\_double\_triangle], [projected], [24/24/0/0], [0/0/0/0], [-2 to -0.999991766], [---], [0.9985439067],
    [paper\_figure\_b1\_double\_triangle\_soft\_ir\_bare], [localized], [196/196/0/62], [196/196/0/62], [-4.00000361 to 2], [-4.00000361 to 2], [0],
    [paper\_figure\_b1\_double\_triangle\_soft\_ir], [localized], [196/196/0/0], [196/196/0/0], [-6 to -0.999952054], [-6 to -0.999952054], [0.9996499398],
  )]
  , kind: table
  )

== Top-bubble soft comparison
<top-bubble-soft-comparison>
Profile `S(e6)`, 25 points 10^-2→10^-5, seed1337. The table reports measure-adjusted IR scaling s = b + 3\*n\_soft, where b is the raw exponent fitted by the IR profiler. The parser counts one soft vector in S(e6), so n\_soft=1 and b=s−3 (process/ir.rs:700,1768,1908; API profile.rs:524). These numbers must not be described as raw amplitude exponents. The Rust test requires sU\<−1, sH\>−0.5, improvement\>=1.5 and the specified R² bounds; those targeted improvement assertions pass. The generic CLI infrared verdict separately requires s\>0, and the negative summed H scalings fail it.

#figure(
  align(center)[#table(
    columns: 7,
    align: (auto,auto,right,right,right,right,right,),
    table.header([Fixture], [Route], [U reported scaling s], [H reported scaling s], [H−U], [U R²], [H R²],),
    table.hline(),
    [child-only], [localized], [-2.001104448869], [-0.000573378786], [2.000531070083], [0.999999970397], [0.999999985020],
    [child-only], [explicit], [-2.001104448869], [-0.000573378786], [2.000531070083], [0.999999970397], [0.999999985020],
    [child-only], [projected], [-2.001104448869], [-0.000573378786], [2.000531070083], [0.999999970397], [0.999999985020],
    [full], [localized], [-2.001104448869], [-0.000572806050], [2.000531642819], [0.999999970397], [0.999999985086],
    [full], [explicit], [-2.001104448869], [-0.000572806050], [2.000531642819], [0.999999970397], [0.999999985086],
    [full], [projected], [-2.001104448869], [-0.000572806050], [2.000531642819], [0.999999970397], [0.999999985086],
  )]
  , kind: table
  )

Each summed and orientation pair has identical LMB \[6,8\] and ray digest `94ac2449745d853c` across U/H (and bare where present). The full bare summed measure-adjusted scaling is -2.001114489716347, R²=0.9999999696759526 in all routes.

#figure(
  align(center)[#table(
    columns: 8,
    align: (auto,right,right,auto,auto,auto,auto,right,),
    table.header([Localized fixture], [Orientation count], [U scaling \<-1], [Bad-U reported scaling], [H reported scaling across all], [H reported scaling on bad U], [Improvement on bad U], [Minimum R² of bad-U/H fits],),
    table.hline(),
    [child-only], [30], [8], [-2.00142671 to -2.0003259], [-0.0016366931 to 4.03209512], [-0.0016366931 to -0.000138874746], [1.99979001 to 2.00018702], [0.999999818650],
    [full], [30], [8], [-2.00142671 to -2.0003259], [-0.00163667739 to 4.08074755], [-0.00163667739 to -0.000138574551], [1.99979003 to 2.00018732], [0.999999818657],
  )]
  , kind: table
  )

=== Raw exponents and generic strict-IR verdict
<raw-exponents-and-generic-strict-ir-verdict>
The raw exponent is b=s−3 for this one-soft-vector limit. The approximately two-power improvement is invariant under this shift; the raw H exponent is near the logarithmic boundary b=−3. These finite-window measurements establish the targeted relative improvement, but do not establish strict IR convergence or pointwise boundedness.

#figure(
  align(center)[#table(
    columns: 5,
    align: (auto,auto,right,right,auto,),
    table.header([Fixture], [Route], [Raw U exponent b], [Raw H exponent b], [Strict CLI U / H pass],),
    table.hline(),
    [child-only], [localized], [-5.001104448869], [-3.000573378786], [False / False],
    [child-only], [explicit], [-5.001104448869], [-3.000573378786], [False / False],
    [child-only], [projected], [-5.001104448869], [-3.000573378786], [False / False],
    [full], [localized], [-5.001104448869], [-3.000572806050], [False / False],
    [full], [explicit], [-5.001104448869], [-3.000572806050], [False / False],
    [full], [projected], [-5.001104448869], [-3.000572806050], [False / False],
  )]
  , kind: table
  )

Every saved top-bubble soft report has all\_passed=false (19/19 reports). All six summed H reports fail the strict s\>0 criterion, as do all six summed U and three summed bare reports. Each localized orientation report (child-only/full × U/H) has 22 strict passes and eight strict failures; the same eight orientations fail in U and H despite their approximately two-power improvement. No per-orientation bare report was requested.

DGSE uses exactly the same convention: its S(e0) also has n\_soft=1 and s=b+3. Its stronger targeted condition sH\>3 means bH\>0, whereas sBare\<=3 means bBare\<=0; improvement\>=1 is unchanged by subtracting three. The completed DGSE test verifies these bounds and R²\>=0.9, but exact DGSE soft fit values were not persisted.

The high R² requirement applies to the eight bad-U orientation pairs. All 30 H reported scalings meet the targeted lower bound of −0.5, but other orientations include low-quality power fits: minimum H R² across all 30 is 0.555708705283 (child-only) and 0.862094728884 (full). The eight enhancement-removal fits meet the strong R² criterion above.

== Bare-control failure breakdown
<bare-control-failure-breakdown>
Resolved here means a fit exists; a fit can still fail its quality criterion. This breakdown prevents unstable flat fits from being counted as evidence of a resolved divergence. Each bare control below contains actual excessive-slope failures, as required by the tests.

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,auto,auto,),
    table.header([Fixture], [Route], [Summed slope / unstable failures], [Orientation slope / unstable failures],),
    table.hline(),
    [gamma\_star\_ddbar\_top\_bubble\_bare], [localized], [8 / 9], [224 / 270],
    [gamma\_star\_ddbar\_top\_bubble\_bare], [explicit], [8 / 9], [0 / 0],
    [gamma\_star\_ddbar\_top\_bubble\_bare], [projected], [8 / 9], [0 / 0],
    [local\_nested\_soft\_top\_self\_energy\_bare], [localized], [13 / 0], [170 / 0],
    [local\_nested\_soft\_top\_self\_energy\_bare], [explicit], [13 / 0], [0 / 0],
    [local\_nested\_soft\_top\_self\_energy\_bare], [projected], [13 / 0], [0 / 0],
    [paper\_figure\_b1\_double\_triangle\_bare], [localized], [8 / 16], [144 / 168],
    [paper\_figure\_b1\_double\_triangle\_bare], [explicit], [8 / 16], [0 / 0],
    [paper\_figure\_b1\_double\_triangle\_bare], [projected], [8 / 16], [0 / 0],
    [paper\_figure\_b1\_double\_triangle\_soft\_ir\_bare], [localized], [28 / 34], [28 / 34],
  )]
  , kind: table
  )

== Appendix B.1 scope
<appendix-b1-scope>
The failed Appendix case has only a localized child-only UV artifact; no full-forest UV report, bare report, or explicit/projected UV report exists. The six initial-degree-one child-hard summed rays actually measured are:

#figure(
  align(center)[#table(
    columns: 4,
    align: (auto,auto,right,right,),
    table.header([Free], [Fixed], [Slope], [R²],),
    table.hline(),
    [\[4\]], [\[6\]], [-1.000000001195], [1.000000000000],
    [\[3\]], [\[6\]], [-1.000000001784], [1.000000000000],
    [\[4\]], [\[5\]], [-0.999999998526], [1.000000000000],
    [\[3\]], [\[5\]], [-1.000000000202], [1.000000000000],
    [\[4\]], [\[2\]], [-0.999999998977], [1.000000000000],
    [\[3\]], [\[2\]], [-1.000000000315], [1.000000000000],
  )]
  , kind: table
  )

The child-only profile is allowed to contain unremoved overall divergences; this control asserts only the six initial-degree-one rays. Full-forest numerical coverage cannot be inferred from this artifact.

== Final massless Figure B.1 scope and failure
<final-massless-figure-b1-scope-and-failure>
The case failed after 6326.706 seconds (105 minutes 26.706 seconds), with test exit code 101; it did not time out. The explicit direct-3D route panicked in bytes 1.11.1 with `advance out of bounds: the len is 1576152287 but advancing by 5871119583` (bytes/src/lib.rs:170), propagated by the fixture thread join at tests/tests/uv.rs:5549. The log omits the originating backtrace, so this records the observed failure boundary without identifying its upstream cause.

The completed localized route covers one graph, 28 LMBs and 196 summed hard-limit records, plus the same 196 records for the single retained orientation `00---+--+-|sigma(0)`. All 392 subtracted fits resolve and pass: measure-adjusted UV slopes range −5.99999999886739 to −0.9999520541631232, minimum R²=0.9996499398096927, with zero vanishing exemptions. Matched bare identities are exactly equal; each 196-record set has 62 failures, comprising 28 excessive-slope and 34 unstable-fit failures.

Only the localized H and bare UV reports exist. Explicit generation started but produced no UV report; the projected route never started and has no state or UV report. This failed test therefore does not validate either later route.

The command wrapper measured 6327.580892889 seconds wall, 6416.208603 user CPU seconds, 18080.821573 system CPU seconds and 329213524 KiB (313.962482 GiB) largest-child peak RSS. These command-wide resource measurements are distinct from per-stage or aggregate simultaneous memory measurements.
