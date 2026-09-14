# Independent local H/Z approach analysis

Completed point files: 84; complete rows: 84; retained failures: 0. Native checks: 4952/4952 passed.

Independent Decimal380 postprocessing of retained native results only; slopes use actual forwarded R_map, not historical WH-scaled radii. Full Y uses actual selected source/grid factors. Local two-direction finite-range evidence, no global bound/MC convergence claim.

Beta is −d log(abs(value))/d log(R); fits use the last three actual native forward radii. Complete source weights Y are in pb. qsum excludes discrete/grid probabilities, which are included once in Y.

| Mode | Direction | beta F norm | beta qsum | beta Y norm | smallest-R Re Y [pb] | Im Y [pb] | Y norm [pb] |
|---|---|---:|---:|---:|---:|---:|---:|
| cuts | Hplus_Zminus | 1.000003 | 0.000000 | 1.000002 | -4163059.1 | 9479552.8 | 10353404 |
| cuts | Hplus_Zplus | 1.000002 | 0.000000 | 1.000002 | -243876.85 | 11935701 | 11938192 |
| cuts_combined_joint | Hplus_Zminus | 1.000003 | 1.000000 | 0.000003 | -0.00036865746 | 0.00083945669 | 0.00091684014 |
| cuts_combined_joint | Hplus_Zplus | 1.000002 | 1.000000 | 0.000002 | -2.1596385e-05 | 0.0010569596 | 0.0010571802 |
| cuts_combined_joint_soft | Hplus_Zminus | 1.000003 | 1.000000 | 0.000003 | -0.00042132281 | 0.00095937907 | 0.0010478173 |
| cuts_combined_joint_soft | Hplus_Zplus | 1.000002 | 1.000000 | 0.000002 | -2.4681582e-05 | 0.0012079538 | 0.0012082059 |
| cuts_joint | Hplus_Zminus | 1.000003 | 1.000000 | 0.000003 | -0.0029120369 | 0.0066308948 | 0.0072421491 |
| cuts_joint | Hplus_Zplus | 1.000002 | 1.000000 | 0.000002 | -0.00017059052 | 0.0083489569 | 0.0083506996 |
| cuts_joint_soft | Hplus_Zminus | 1.000003 | 1.000000 | 0.000003 | -0.0033280421 | 0.0075781655 | 0.0082767419 |
| cuts_joint_soft | Hplus_Zplus | 1.000002 | 1.000000 | 0.000002 | -0.0001949606 | 0.0095416651 | 0.0095436566 |
| cuts_soft | Hplus_Zminus | 1.000003 | -0.000000 | 1.000003 | -2561667.9 | 5833082.1 | 6370791.9 |
| cuts_soft | Hplus_Zplus | 1.000002 | 0.000000 | 1.000002 | -150065.48 | 7344430.8 | 7345963.7 |
| optimized_lmb | Hplus_Zminus | 1.000003 | -0.000000 | 1.000003 | -585082.15 | 1332269.6 | 1455081.9 |
| optimized_lmb | Hplus_Zplus | 1.000002 | 0.000000 | 1.000002 | -34274.793 | 1677459.9 | 1677810.1 |

Actual native R_map ranges [GeV]:

- cuts, Hplus_Zminus: 0.251623841389 → 2.51623847522e-06.
- cuts, Hplus_Zplus: 0.251623841388 → 2.51623841647e-06.
- cuts_combined_joint, Hplus_Zminus: 0.251623841389 → 2.51623841472e-06.
- cuts_combined_joint, Hplus_Zplus: 0.251623841388 → 2.51623841222e-06.
- cuts_combined_joint_soft, Hplus_Zminus: 0.251623841389 → 2.51623841472e-06.
- cuts_combined_joint_soft, Hplus_Zplus: 0.251623841388 → 2.51623841222e-06.
- cuts_joint, Hplus_Zminus: 0.251623841389 → 2.51623841472e-06.
- cuts_joint, Hplus_Zplus: 0.251623841388 → 2.51623841222e-06.
- cuts_joint_soft, Hplus_Zminus: 0.251623841389 → 2.51623841472e-06.
- cuts_joint_soft, Hplus_Zplus: 0.251623841388 → 2.51623841222e-06.
- cuts_soft, Hplus_Zminus: 0.251623841389 → 2.51623847522e-06.
- cuts_soft, Hplus_Zplus: 0.251623841388 → 2.51623841647e-06.
- optimized_lmb, Hplus_Zminus: 0.251623841389 → 2.51623852583e-06.
- optimized_lmb, Hplus_Zplus: 0.251623841388 → 2.51623836624e-06.

Largest actual-source radius displacement from the common native anchor: 4.415459E-8 relative; this is recorded, not treated as identical-point physics.
Largest coordinate displacement from its intended binary64 raw target: 1.481737E-13 GeV. Different maps recover slightly different binary64 cubes; actual forwarded points are retained and are not claimed identical.
Maximum same-map abs(J*w*qsum−1): 1.567125E-291.

All 56 base point files match the retained hash manifest. The four new anchors extend each direction by two decades; fits combine original and extended actual forwarded radii. Raw binary64 momenta, common native comparison anchors and each candidate actual forward point are distinct and retained; no raw F / common-anchor q coincidence is assumed. Quad rounding allowance follows the actual106-bit backend; physical stability budgets remain unchanged.
Maximum ordinary/Arb norm difference: 2.034931E-7.
Maximum six-event sum / complete-estimator norm difference: 3.385079E-15.
Maximum decomposition error relative to sum of term norms: 4.755161E-16.

Detailed check evidence and original-file hashes are in independent_extended_approach_analysis.json. Event/decomposition weights are already normalized by the precise result owner, including outer grid weight; they are not multiplied again.
