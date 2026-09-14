# Matched GL638 sampling cost: 072c9999e to 1f83233dd

This report analyzes the two completed, untraced, optimized 20-worker runs. It executes no integrand. All 936 orientations and six physical cuts remain active. Every original nested source Sample, hard raw case, retained maximum, native total, precision, factor and validity flag matches between builds in both modes. Both runs preserve the same 35 saved-state file hashes.

Each mode has 960 representative calls (320 sources, three repetitions), 180 selected-hard calls, 180 bare-hard controls and 60 maximum replays. Each worker warms once before measurement. W is elapsed precise-call time; S is the existing sampling counter; P is disjoint actual physical-body time, including probes and retries. The conservative bound is sum(max(W−P,0))/sum(P). Grid/RNG generation is excluded. The ordinary configured stability stack and fresh-iteration max_eval=0 are used; forced-Arb controls are excluded from these ratios.

| Mode and set | Before S/P | After S/P | Before upper/P | After upper/P | Invalid after |
|---|---:|---:|---:|---:|---:|
| optimized_lmb: all_representative | 3.89% | 3.69% | 4.45% | 4.32% | 0 / 960 |
| optimized_lmb: valid_representative_only_diagnostic | 3.89% | 3.69% | 4.45% | 4.32% | 0 / 960 |
| optimized_lmb: selected_hard | 2.33% | 2.32% | 4.38% | 3.89% | 0 / 180 |
| optimized_lmb: bare_hard_control | 0.03% | 0.03% | 1.78% | 1.77% | 0 / 180 |
| optimized_lmb: representative_maxima | 6.22% | 6.14% | 6.52% | 6.42% | 0 / 60 |
| joint_hz_plus_lmb: all_representative | 40.93% | 39.48% | 41.25% | 39.76% | 3 / 960 |
| joint_hz_plus_lmb: valid_representative_only_diagnostic | 47.45% | 45.95% | 47.82% | 46.27% | 0 / 957 |
| joint_hz_plus_lmb: selected_hard | 18.65% | 18.08% | 18.87% | 18.33% | 0 / 180 |
| joint_hz_plus_lmb: bare_hard_control | 0.02% | 0.03% | 0.27% | 0.22% | 0 / 180 |
| joint_hz_plus_lmb: representative_maxima | 58.20% | 47.51% | 58.48% | 47.82% | 0 / 60 |

The valid-only row is a diagnostic. The primary totals retain all invalid calls; the joint configuration fails numerical and timing acceptance. Its three failures are repetitions of worker 12, draw 12, channel 0, with the same Quad and Arb rotation discrepancy near 2.1685e−10. Their slow physical rescues contribute 10.42 s to the new physical denominator; they must not be used to hide the 46.27% valid-call conservative overhead.

| Selected hard point, joint mode | Before upper/P | After upper/P | Native lane |
|---|---:|---:|---|
| hz_00_stack | 45.70% | 44.39% | Double |
| hz_03_stack | 48.80% | 43.33% | Double |
| soft22_ct_on_stack | 7.68% | 8.19% | Quad |

Joint representative S totals change from 27.9771 s to 27.8895 s (−0.31%). Mean channel-0 S changes from 56.96 to 55.82 ms; all other channel means remain near 22–28 ms. Representative upper/P p99 is 159.32%, with a 199.04% maximum. The maximum-replay aggregate improves to 47.82%, but its p99 is still 91.34%. These are ratios across matched calls, not integrator variance or bounded-weight evidence.

Baseline aggregate bounds remain below 10% (4.32% representative, 3.89% selected hard, 6.42% maximum), so no baseline optimization is warranted. This does not assert every individual baseline call is below 10%. The joint mode still requires substantial optimization. The worker clone ownership correction is important for actual Integrate, but both timing clients explicitly rewarm clones, so mutex sharing is not the explanation for this matched result. No compact/fallback visitation counts are claimed: proposal-policy records are not serialized.

Detailed metrics, hashes and retained-failure identities: [comparison_072_1f83233.json](comparison_072_1f83233.json). Input directories: `results-timing1` and `results-performance-timing1`. The latter uses the frozen optimized build of 1f83233dd79b5c7b8b09aef38ad5acdd5fde932c; the former uses 072c9999ee928e8e9979f92bd5fb16d4f17ca15f. This comparison makes no full GL638 normalization, integration-gain, CT-star cure or bounded-weight claim.
