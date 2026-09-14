# Final certificate timing comparison

This compares completed untraced 20-worker GL638 timings before and after reuse of radius-independent certificate data within one preparation. It analyzes retained reports only; no state or integrand is executed. Every original nested Sample, hard raw case, retained maximum, native total, precision, factor, validity flag and stability result matches exactly in both modes. Runtime settings, all 936 orientation identities, six cuts and all 35 state hashes are unchanged.

The before build is 1f83233dd79b5c7b8b09aef38ad5acdd5fde932c. The after build is optimized build 3, with the exact base/patch/build records retained in comparison_1f_certificate.json. Both use the actual configured ordinary stability stack, zero recorded previous maximum and the same 320 representative source draws repeated three times per mode. Forced-Arb diagnostic controls do not enter these timing denominators.

| Mode and set | Before S/P | After S/P | Before upper/P | After upper/P | After invalid/count |
|---|---:|---:|---:|---:|---:|
| optimized_lmb: all_representative | 3.69% | 4.16% | 4.32% | 4.59% | 0/960 |
| optimized_lmb: valid_representative_only_diagnostic | 3.69% | 4.16% | 4.32% | 4.59% | 0/960 |
| optimized_lmb: selected_hard | 2.32% | 2.41% | 3.89% | 4.76% | 0/180 |
| optimized_lmb: bare_hard_control | 0.03% | 0.02% | 1.77% | 1.00% | 0/180 |
| optimized_lmb: representative_maxima | 6.14% | 6.78% | 6.42% | 7.05% | 0/60 |
| joint_hz_plus_lmb: all_representative | 39.48% | 35.37% | 39.76% | 35.67% | 3/960 |
| joint_hz_plus_lmb: valid_representative_only_diagnostic | 45.95% | 41.45% | 46.27% | 41.79% | 0/957 |
| joint_hz_plus_lmb: selected_hard | 18.08% | 16.86% | 18.33% | 17.05% | 0/180 |
| joint_hz_plus_lmb: bare_hard_control | 0.03% | 0.02% | 0.22% | 0.21% | 0/180 |
| joint_hz_plus_lmb: representative_maxima | 47.51% | 47.77% | 47.82% | 48.10% | 0/60 |

S is the existing sampling counter; P is disjoint actual all 936 physical-body time including probes/retries; the conservative upper is sum(max(W−P,0))/sum(P), where W is elapsed precise-call time. Grid/RNG generation is excluded. Summed worker durations are not CPU time or elapsed batch wall time. Warmup is excluded. The valid-only row is a diagnostic; all failures remain in the primary totals.

| Joint selected hard point | Before upper/P | After upper/P | After mean S, ms |
|---|---:|---:|---:|
| hz_00_stack | 44.39% | 45.36% | 18.790 |
| hz_03_stack | 43.33% | 44.76% | 18.836 |
| soft22_ct_on_stack | 8.19% | 7.20% | 16.843 |

| Joint set | Before p99/max upper/P | After p99/max upper/P |
|---|---:|---:|
| all_representative | 159.32% / 199.04% | 146.78% / 159.36% |
| selected_hard | 73.86% / 81.20% | 54.23% / 54.93% |
| representative_maxima | 91.34% / 91.34% | 79.21% / 79.21% |

Joint representative S decreases from 27.8895 to 23.6353 seconds (15.25%). The unmodified baseline S also varies by −6.04%, so do not attribute all elapsed-time changes to the patch. The exact-output equality supports unchanged represented behavior; this single matched batch is not a universal speed guarantee. The joint configuration remains above the former 10% target. Per-point and tail results are retained rather than hidden by averaging.

The same worker 12/draw 12 joint source remains invalid in all three repetitions. Its Quad discrepancy is 2.1685013932030374e−10 and Arb discrepancy2.1685013932032444e−10, exactly unchanged. Those invalid calls contribute 10.2261 seconds to physical time; they cannot be used to certify cost or numerical acceptance. The independently traced center discrepancy remains the explanation under investigation; this certificate optimization does not change that physical center chooser.

As now requested by the user, performance tuning stops after this round regardless of these ratios. Next work prioritizes the canonical physical-center repair and actual-state Gaussian correctness, followed by matched equal-N ordinary/nonjoint/joint physics comparisons, signed and absolute integral errors, maxima and their origin, H/Z-corner scaling, and a defensible GL638 central value and uncertainty. This report claims neither normalization nor integration gain nor bounded weights.

Inputs: results-performance-timing1 and results-certificate-timing1. Exact hashes, call counts, per-channel means, per-point ratios and unchanged invalid-source identities are retained in [comparison_1f_certificate.json](comparison_1f_certificate.json).
