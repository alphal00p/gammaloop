All 21 runs passed: 3 × 2048 samples per method, 20 workers; the 30-worker control is excluded.

Means and empirical errors below are in 10⁻⁴ pb; maxima are absolute signed-component maxima in pb.

| Method | Re ± SE | Abs Re ± SE | Im ± SE | Abs Im ± SE | Max Re / Im (pb) | Score |
|---|---:|---:|---:|---:|---:|---:|
| optimized_lmb | 1.4317 ± 1.927 | 6.1834 ± 1.925 | 1.0504 ± 2.465 | 8.9363 ± 2.463 | 1.04892 / 1.34675 | 596.6986 |
| cuts | 1.8994 ± 1.471 | 2.1214 ± 1.471 | 1.915 ± 1.699 | 2.3827 ± 1.699 | 0.875679 / 1.0342 | 347.7284 |
| cuts_joint | 0.19781 ± 0.5581 | 1.1609 ± 0.5579 | 0.35659 ± 0.2857 | 0.8743 ± 0.2855 | 0.224439 / 0.112123 | 50.04617 |
| cuts_soft | -4.5453 ± 6.468 | 10.479 ± 6.467 | 20.417 ± 19.47 | 23.437 ± 19.47 | 3.76756 / 11.9369 | 29160.24 |
| cuts_joint_soft | -5.4561 ± 7.389 | 11.7 ± 7.388 | 23.353 ± 22.25 | 26.571 ± 22.25 | 4.30578 / 13.6421 | 38085.04 |
| cuts_combined_joint | 0.18392 ± 0.558 | 1.1552 ± 0.5578 | 0.35645 ± 0.2853 | 0.88482 ± 0.2851 | 0.224439 / 0.112123 | 50.02642 |
| cuts_combined_joint_soft | -5.4704 ± 7.389 | 11.711 ± 7.388 | 23.36 ± 22.25 | 26.582 ± 22.25 | 4.30578 / 13.6421 | 38085.04 |

Selected a = optimized_lmb, b = cuts, c = cuts_combined_joint. The composed score is only 0.03947% below direct joint; this screen does not establish superiority between them.

167 independent checks passed. All 93 recorded input hashes still match; all 35 saved-state hashes are unchanged in each batch. Confirmation retains the frozen methods/cards/source/scales and uses independent seeds.
