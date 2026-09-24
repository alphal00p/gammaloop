= Fixed GL638 Monte Carlo comparison
<fixed-gl638-monte-carlo-comparison>
Approved protocol, 2026-09-14. Analyze existing `Integrate` outputs; no new integration or accumulation engine. Physical validity and the actual-catalogue Gaussian acceptance precede the screen. Performance is reported, not optimized or used to block the physics experiment.

All seven proposals retain the same full 936-orientation GL638 graph contribution, six physical cuts, UV terms, threshold metadata, WH/WF, localization and stability stack. Proposal hosts never select physical cuts. The exact cards are in `cut-matched/cards_combined_seven_manifest.json`.

#figure(
  align(center)[#table(
    columns: 3,
    align: (auto,right,auto,),
    table.header([Proposal], [Channels], [Purpose],),
    table.hline(),
    [optimized\_lmb], [6], [Ordinary baseline a],
    [cuts], [6], [Six LU-h cut channels],
    [cuts\_joint], [7], [Same cuts plus direct H/Z joint channel],
    [cuts\_soft], [7], [Same cuts plus LMB(6,12,13,14)],
    [cuts\_joint\_soft], [8], [Same cuts, soft LMB and direct joint],
    [cuts\_combined\_joint], [7], [Same cuts plus one sequential cut1 LU-h/HZ channel],
    [cuts\_combined\_joint\_soft], [8], [Same cuts, soft LMB and sequential cut1 LU-h/HZ],
  )]
  , kind: table
  )

The last two test an actual conditional composition, not merely a mixture of cut and joint channels. They keep all six ordinary LU-h cut channels. Each joint/no-joint contrast retains exactly the same soft-LMB choice. Joint-only versus joint-plus-cut comparisons from earlier research are not these contrasts.

Screen all seven using seeds 11113, 22229, 33331 and 2048 samples per seed: 43,008 draws. Use one fresh untrained iteration per seed, with equal N, workers and batch settings. Record all signed Re/Im and absolute |Re|/|Im| means and errors, extrema, precision/failure counts, channel occupancies and time. Post-iteration PDFs describe an updated grid; they are not the PDFs that generated this first batch.

For each of the four observables, pool equal-N runs using the existing X2 analysis formula. With n samples per seed, k seeds and M=kn,

```
mu = sum(mu_i)/k
M2 = sum[n(n-1) error_i^2 + n(mu_i-mu)^2]
SE = sqrt(M2/[M(M-1)])
sample_variance = M SE^2
```

Also report the SE estimated from the k seed means, and individual seed means and errors. These are finite-sample empirical diagnostics, not a certificate of finite variance or tail convergence. Same seeds pair RNG schedules, not physical points. Matched differences use the paired seed means; their errors are not computed by blindly adding the two pooled error variances.

Freeze one common scale per phase: the baseline screen estimate of |Re| or |Im|. For each joint candidate define its score as the largest of its four sample variances divided by the corresponding phase scale squared. No signed mean is used as a normalization denominator. A zero baseline absolute scale prevents automatic selection; no numerical floor is introduced. Invalid, unstable, incomplete or ungated runs prevent eligibility and must remain visible.

Select c by the smallest score among the four joint candidates, with declared proposal order breaking exact ties. Set b to `cuts_soft` if c has the soft LMB, otherwise `cuts`; a remains `optimized_lmb`. This freezes a matched joint removal experiment. It does not independently pick a different soft policy for b. Publish all seven screen results, including losing or invalid candidates.

Before observing confirmation data, freeze a/b/c, card/source hashes and the screen scales. Confirm these three using fresh seeds 10007, 20011, 30013, 40009, 50021, 32768 samples each: 163840 per method, 491520 total. Do not select a new winner from confirmation or append samples selectively. Quote the preselected c estimate and empirical error as the primary current GL638 result, alongside a and its matched b, even if the anticipated improvement fails. Any additional equal-N block is a separately reported experiment.

Retain original `integration_state.bin` and decode extrema through its existing owner. Printed coordinate strings do not preserve complete Sample/grid weights. Replay each unique maximum with its original complete Sample, card and native stability stack; first-iteration max\_eval=0 is the original escalation context. Retain aliases when one source supplies multiple phase/sign extrema. Inspect forced Arb controls and native per-cut/threshold decomposition at the largest excursions. Returned precise event weights already include the outer grid factor; do not multiply it again. The statistics script never reconstructs Samples or evaluates physics.

Report maximum absolute weight, maximum single-sample contribution relative to the absolute estimate, and its fraction of the empirical second moment. Maxima growth across finite N is evidence about tails, not a boundedness proof. The separate direct H/Z scan levels off along two regular directions; star images, soft boundaries, tangencies and other sources of large weights remain open.
