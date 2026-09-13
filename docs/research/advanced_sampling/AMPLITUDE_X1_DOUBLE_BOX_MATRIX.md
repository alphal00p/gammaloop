# X1 double-box and boosted-kite sampling pilots

The C-only channel improves the double box’s real-component variance and maximum in this pilot, while the mixed channel helps the imaginary component. There is no single winner for both components.

Same frozen X1 source `3f9b1a28473f313d6242a556ada2e8dbcd6f6477` and binary SHA256 `f700ee575153067ccc77894ea776a96f139105812a15014d094dd67ba7db6d61` as the [kite pilot](AMPLITUDE_X1_KITE_MATRIX.md). All 98 double-box orientations are summed, threshold CTs are enabled, power 2, exact map-density partition, real training. Seeds 101/202/303 each have 10,000 draws in five batches on one core. All 12 integrations completed with no reported NaN/unstable points.

| Channels | Mean Re / Im | RMS error Re / Im | Max absolute weight Re / Im | Evaluation µs |
|---|---:|---:|---:|---:|
| optimized_lmb | -1.8515e-08 / -2.7573e-07 | 1.152e-08 / 1.746e-08 | 1.037e-04 / 1.137e-04 | 90.87 |
| all_lmb | -1.8515e-08 / -2.7573e-07 | 1.152e-08 / 1.746e-08 | 1.037e-04 / 1.137e-04 | 91.23 |
| surface | -2.5140e-08 / -2.7207e-07 | 6.776e-09 / 1.803e-08 | 2.144e-05 / 1.275e-04 | 85.37 |
| mixed | -2.3228e-08 / -2.6601e-07 | 1.145e-08 / 1.555e-08 | 1.408e-04 / 6.955e-05 | 95.41 |

C is `surface(5,6,10)` with full parent/subspace `[5,8]`. Both LMB selectors again yield identical samples/results. Relative to their baseline, C-only reported variance ratios are 0.346 (real) and 1.065 (imaginary); maximum ratios are 0.207 and 1.122. Mixed variance ratios are 0.988 and 0.793, with maximum ratios 1.357 and 0.612. Means are compatible at these short-run uncertainties. Three paired seeds do not establish asymptotic variance or bounded weights, and the real training objective does not optimize imaginary error.

| Reference channels | Normalization | Raw second moment |
|---|---:|---:|
| optimized_lmb | 0.98936452 | 8.93863347 |
| all_lmb | 0.98936452 | 8.93863347 |
| surface | 0.99472564 | 9.01221336 |
| mixed | 0.98962851 | 8.96729690 |

Reference target is 1 and second moment 9.1025, with exactly the same 8,192-point/channel shifted-Gaussian rule and 0.02/0.2 absolute gates as the kite. All pass. Per-channel normalization dispersion is retained; correlated contributions are not assigned a fabricated independent aggregate error. Moment dispersion was not collected. This is the frozen checkpoint’s f64 reference route, separate from the current native-reference development.

The existing boosted kite (`Q=(sqrt(26),0,0,1)`, invariant mass 5) also completed the four reference gates followed by all twelve 10,000-draw integrations. All 18 orientations are summed, with the same source, four profiles, seeds, power and one-core settings. Input states remained unchanged, and no run reported an invalid final sample. The modest boost keeps the selected full-rank C surface compatible with the checkpoint’s zero-center map; larger boosts requiring a new interior center are outside this evidence.

| Boosted-kite channels | Mean Re / Im | RMS error Re / Im | Max absolute weight Re / Im | Evaluation µs |
|---|---:|---:|---:|---:|
| optimized_lmb | -1.3698e-05 / 6.7005e-06 | 3.517e-07 / 1.062e-06 | 1.701e-03 / 1.158e-03 | 81.30 |
| all_lmb | -1.3698e-05 / 6.7005e-06 | 3.517e-07 / 1.062e-06 | 1.701e-03 / 1.158e-03 | 80.84 |
| surface | -1.4000e-05 / 7.7544e-06 | 8.759e-07 / 1.828e-06 | 5.027e-03 / 7.285e-03 | 74.64 |
| mixed | -1.3742e-05 / 5.7512e-06 | 3.017e-07 / 9.878e-07 | 8.137e-04 / 8.633e-04 | 83.84 |

C is again `surface(2,4,6)` with parent/subspace `[4,6]`. Mixed variance ratios are 0.736 (real) and 0.864 (imaginary), and its maxima are 52% and 25% smaller than LMB. Evaluation cost rises by about 3.1%. C-only variance ratios are 6.20 and 2.96, so it again performs worse. The imaginary mixed mean is lower than the LMB mean but still within these short-run reported uncertainties; these paired seeds do not support a formal independent-error significance test.

| Boosted-kite reference channels | Normalization | Raw second moment |
|---|---:|---:|
| optimized_lmb | 0.98292699 | 8.95521153 |
| all_lmb | 0.98292699 | 8.95521153 |
| surface | 0.98592521 | 8.96751469 |
| mixed | 0.98880972 | 9.00632004 |

All boosted reference checks pass the same targets and absolute tolerances above. Their frozen f64-path and dispersion limitations are identical. Together with the rest-kite pilot, this completes the bounded three-graph/frame inventory: 36 physical runs and twelve profile-specific reference checks. The mixed proposal helps the real component of both kite frames; the double box illustrates why a single full-rank surface does not guarantee simultaneous improvement in both components.

The [portable JSON](AMPLITUDE_X1_DOUBLE_BOX_MATRIX.json) retains the double-box study at the root and the complete `boosted_kite_matrix` study within it: source/binary/state/card hashes, original graph and generation inputs, orientation keys, driver source, per-run settings and argv, normalization/moment contributions, exact counts, timing/stability and maximum-weight points. Scratch paths are provenance; reproduce by extracting the embedded inputs and replacing those paths. No new topology or production source was generated for these measurements. Next compare the new proper-subspace channels on this same inventory after their implementation passes its gates, and use paired wall-time budgets to test whether the small fixed-count gains survive equal computational cost.
