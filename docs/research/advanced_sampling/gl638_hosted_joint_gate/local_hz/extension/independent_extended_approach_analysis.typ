= Independent local H/Z approach analysis
<independent-local-hz-approach-analysis>
Completed point files: 36; complete rows: 36; retained failures: 0. Native checks: 2040/2040 passed.

Independent Decimal380 postprocessing of retained native results only; slopes use actual forwarded R\_map, not historical WH-scaled radii. Full Y uses actual selected source/grid factors. Local two-direction finite-range evidence, no global bound/MC convergence claim.

Beta is −d log(abs(value))/d log(R); fits use the last three actual native forward radii. Complete source weights Y are in pb. qsum excludes discrete/grid probabilities, which are included once in Y.

#figure(
  align(center)[#table(
    columns: 8,
    align: (auto,auto,right,right,right,right,right,right,),
    table.header([Mode], [Direction], [beta F norm], [beta qsum], [beta Y norm], [smallest-R Re Y \[pb\]], [Im Y \[pb\]], [Y norm \[pb\]],),
    table.hline(),
    [direct\_h\_p2], [Hplus\_Zminus], [1.000003], [0.500247], [0.499755], [-922.83846], [2101.3624], [2295.0718],
    [direct\_h\_p2], [Hplus\_Zplus], [1.000002], [0.500247], [0.499755], [-54.060955], [2645.8246], [2646.3768],
    [joint\_hz\_plus\_lmb], [Hplus\_Zminus], [1.000003], [1.000000], [0.000003], [-0.0029120368], [0.0066308947], [0.0072421491],
    [joint\_hz\_plus\_lmb], [Hplus\_Zplus], [1.000002], [1.000000], [0.000002], [-0.00017059052], [0.0083489569], [0.0083506995],
    [optimized\_lmb], [Hplus\_Zminus], [1.000003], [-0.000000], [1.000003], [-585082.15], [1332269.6], [1455081.9],
    [optimized\_lmb], [Hplus\_Zplus], [1.000002], [0.000000], [1.000002], [-34274.793], [1677459.9], [1677810.1],
  )]
  , kind: table
  )

Actual native R\_map ranges \[GeV\]:

- direct\_h\_p2, Hplus\_Zminus: 0.251623841389 → 2.51623842234e-06.
- direct\_h\_p2, Hplus\_Zplus: 0.251623841388 → 2.51623840554e-06.
- joint\_hz\_plus\_lmb, Hplus\_Zminus: 0.251623841389 → 2.51623841472e-06.
- joint\_hz\_plus\_lmb, Hplus\_Zplus: 0.251623841388 → 2.51623841222e-06.
- optimized\_lmb, Hplus\_Zminus: 0.251623841389 → 2.51623852583e-06.
- optimized\_lmb, Hplus\_Zplus: 0.251623841388 → 2.51623836624e-06.

Largest actual-source radius displacement from the common native anchor: 4.415459E-8 relative; this is recorded, not treated as identical-point physics. Largest coordinate displacement from its intended binary64 raw target: 1.481737E-13 GeV. Different maps recover slightly different binary64 cubes; actual forwarded points are retained and are not claimed identical. Maximum same-map abs(J#emph[w];qsum−1): 5.435883E-293.

Original24point files remain hash-identical. The four new anchors extend each direction by two decades; fits combine original and extended actual forwarded radii. Raw binary64 momenta, common native comparison anchors and each candidate actual forward point are distinct and retained; no raw F / common-anchor q coincidence is assumed. Quad rounding allowance follows the actual106-bit backend; physical stability budgets remain unchanged. Maximum ordinary/Arb norm difference: 1.591925E-7. Maximum six-event sum / complete-estimator norm difference: 2.831740E-15. Maximum decomposition error relative to sum of term norms: 3.229585E-16.

Detailed check evidence and original-file hashes are in independent\_extended\_approach\_analysis.json. Event/decomposition weights are already normalized by the precise result owner, including outer grid weight; they are not multiplied again.
