# Independent local H/Z approach analysis

Completed point files: 24; complete rows: 24; retained failures: 0. Native checks: 1312/1312 passed.

Independent Decimal380 postprocessing of retained native results only; slopes use actual forwarded R_map, not historical WH-scaled radii. Full Y uses actual selected source/grid factors. Local two-direction finite-range evidence, no global bound/MC convergence claim.

Beta is −d log(abs(value))/d log(R); fits use the last three actual native forward radii. Complete source weights Y are in pb. qsum excludes discrete/grid probabilities, which are included once in Y.

| Mode | Direction | beta F norm | beta qsum | beta Y norm | smallest-R Re Y [pb] | Im Y [pb] | Y norm [pb] |
|---|---|---:|---:|---:|---:|---:|---:|
| direct_h_p2 | Hplus_Zminus | 1.000249 | 0.502482 | 0.497767 | -92.389295 | 210.37289 | 229.76626 |
| direct_h_p2 | Hplus_Zplus | 1.000227 | 0.502482 | 0.497745 | -5.4118981 | 264.88142 | 264.9367 |
| joint_hz_plus_lmb | Hplus_Zminus | 1.000249 | 0.999987 | 0.000262 | -0.0029120417 | 0.0066307965 | 0.0072420611 |
| joint_hz_plus_lmb | Hplus_Zplus | 1.000227 | 1.000001 | 0.000227 | -0.00017057909 | 0.0083488698 | 0.0083506122 |
| optimized_lmb | Hplus_Zminus | 1.000249 | -0.000041 | 1.000290 | -5850.8223 | 13322.478 | 14550.62 |
| optimized_lmb | Hplus_Zplus | 1.000227 | 0.000013 | 1.000214 | -342.72513 | 16774.433 | 16777.933 |

Actual native R_map ranges [GeV]:

- direct_h_p2, Hplus_Zminus: 0.251623841389 → 0.000251623841426.
- direct_h_p2, Hplus_Zplus: 0.251623841388 → 0.000251623841241.
- joint_hz_plus_lmb, Hplus_Zminus: 0.251623841389 → 0.000251623841472.
- joint_hz_plus_lmb, Hplus_Zplus: 0.251623841388 → 0.000251623841222.
- optimized_lmb, Hplus_Zminus: 0.251623841389 → 0.00025162384151.
- optimized_lmb, Hplus_Zplus: 0.251623841388 → 0.000251623841254.

Largest coordinate displacement from its intended binary64 raw target: 1.481737E-13 GeV. Different maps recover slightly different binary64 cubes; actual forwarded points are retained and are not claimed identical.
Maximum same-map abs(J*w*qsum−1): 2.912122E-294.

Maximum ordinary/Arb norm difference: 1.338696E-9.
Maximum six-event sum / complete-estimator norm difference: 2.831740E-15.
Maximum decomposition error relative to sum of term norms: 3.229585E-16.

Detailed check evidence and original-file hashes are in independent_approach_analysis.json. Event/decomposition weights are already normalized by the precise result owner, including outer grid weight; they are not multiplied again.
