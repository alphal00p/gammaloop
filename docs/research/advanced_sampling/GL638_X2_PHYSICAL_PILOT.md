# GL638 X2 direct-H physical pilot

The generic conditional H chart passes normalization, native inverse-density and physical replay checks with the full **936-orientation** GL638 state. This small untrained pilot gives encouraging imaginary-tail results for power 2, but does **not** establish a reliable global efficiency gain: real errors are essentially unchanged against the matched power-1 chart, the mixed catalogue does not improve, and rare excursions dominate several estimates. The direct-H proposal remains unbounded on the measured H/Z rays. It is neither the joint two-normal H/Z chart nor an affine-star pullback.

All results use frozen source `f2f64fb17dd5dab588eca298dc61f1100fa4962d`, optimized CLI SHA256 `a37b5b27e1fde706e585be9f9c5e10894152cb02aab962ef13ff2dd751f0ee8f`. Build/profile and exact library fingerprints are in the [machine artifact](GL638_X2_PHYSICAL_PILOT.json). The [fresh checkpoint state](GL638_CHECKPOINT06D_BASELINE.json), all six physical cuts, original threshold metadata, 3D local UV and integrated UV are retained. State hashes are unchanged. Unlike the earlier reference-only result, the raw replays and integrations here evaluate the complete physical orientation sum.

## Channels and unchanged physics

The named direct-H map is:

```toml
around = "then(complement(6,7,10),block(lmb(3),at_cut(cut(2,6,10),surface(2,4,12))))"
parent_lmb = [3,6,7,10]
subspace_lmb = [3]
on_cut = [1]
```

It samples H in the actual cut-1 LU kinematics; H need not be a CT of that cut. At fixed prepared a,t, `H=E(a+p)+E(p)-B`, `B=E(a+t)+E(t)`, with center `-a/2`. The generation-parent raw frame remains `[3,4,7,10]`. Power 1 is the otherwise identical chart without the singular H focus; power 2 supplies the signed-quadratic focus. With `mapping="linear"`, the ordinary complement/LMB maps ignore this power, making P1/P2 a focused comparison.

The four selections are optimized LMBs; H-P1 alone; H-P2 alone; and H-P2 plus optimized LMBs. Explicit `sampling_channel_weight="map_density"` is used throughout. The six ordinary bases, in canonical order, are `(6,12,13,14)`, `(4,12,13,14)`, `(6,7,13,14)`, `(4,7,13,14)`, `(6,10,13,14)`, `(4,10,13,14)`. The mixed catalogue has seven channels. Unique names prevent nested settings merges from retaining the wrong profile.

The physical multiplier prescription remains `WH=H²/[H²+(P Z/Q)²]` and its complementary sectors, evaluated at the appropriate CT stars, with the existing WF prescription and common-center/subspace grouping. Neither these functions nor threshold localization is altered. The physical LU h remains `poly_exponential`, sigma 1, power 3. This is a direct-H shape experiment; no additional LU-h radial profile is attached. Full settings and card text are embedded in the JSON.

## Normalization and physical gates

The Gaussian reference has width 300 GeV and center `[30,-20,10]` repeated four times; the expected raw second moment is 1,085,600 GeV². The production reference owner sums all channel contributions **per draw before squaring**. Its temporary summed-channel setting is restored to Monte Carlo for physics. The unchanged acceptance criteria are 6% normalization and 8% raw-moment error.

| Proposal | Halton points | Normalization | Raw-moment relative error | Criterion |
|---|---:|---:|---:|---|
| Six optimized LMBs | 8,192 | 1.245420 | +15.71% | Miss |
| Six optimized LMBs | 32,768 | 1.006498 | -2.92% | Pass |
| Direct H, power 1 | 8,192 | 1.033725 | +2.19% | Pass |
| H power 2 + six LMBs | 8,192 | 1.218556 | +22.79% | Miss |
| H power 2 + six LMBs | 32,768 | 0.984819 | -3.51% | Pass |
| Direct H, power 2 | 32,768 | 1.010322 | -3.47% | Pass |

The two 8,192-point misses were preserved and resolved by increasing the same quadrature to 32,768 points, without changing maps, widths or tolerances. All values were finite. The reported reference dispersion uses the IID formula on deterministic Halton values; it is not a randomized-QMC confidence interval. These gates test normalization, not the variance of a randomly selected channel.

Eight raw H/Z rays, the double-soft control with CTs on/off, and the sampling-invariance replays give **28 native physical evaluations**; four x-space mode checks also pass. All six cut events are retained. Ordinary-stack H/Z totals agree with forced Arb to at most `1.41e-9` relative in Re and `2.00e-10` in Im; raw totals and events are exactly unchanged when switching optimized LMBs to H-P2. The double-soft control rescues to Quad; its imaginary part agrees with Arb to about `4e-12`. CT-off Re is numerical zero, so relative error against zero is not meaningful. Full event sums match totals.

The initial forced-Arb CLI inspection stopped because its f64 event-reporting boundary could not represent a native event component `1.06664e-343`. The existing precise API successfully evaluated the same cases with events enabled, preserving all native total/event weights as decimal strings. No event was dropped and no physics/source workaround was applied. The initial failure is retained as reporting-range evidence, not a failed native physical evaluation.

For the fixed rays, `R=sqrt(H²+(|P0| Z/Q)²)`, with `P0=-254.52045717208898 GeV` and `Q=1000 GeV`. R is a threshold-normal energy scale in **cut-prepared** kinematics, not a soft-gluon radius or the map's radial coordinate. The archived raw points lie on that host shell (`t*=1`). The complete inverse density q includes complements and the conditional pullback, has units `GeV^-12`, and is expressed in the full raw generation-parent frame. Native Arb evaluates both F and q at exact promotions of the original binary64 points.

Using `|F|~R^(-beta_F)`, `q~R^(-beta_q)`, and `|F/q|~R^(-beta_W)`, least-squares fits to R = 0.02, 0.002, 0.0002 GeV give:

| Fixed direction | beta_F | beta_q | beta_W |
|---|---:|---:|---:|
| HplusZminus | 1.000249 | 0.502482 | 0.497767 |
| HplusZplus | 1.000227 | 0.502482 | 0.497745 |

The JSON also retains R=0.2 GeV and four-point fits. Thus the remaining approximately `R^-1/2` weight is measured using the **actual full physical sum**. At R=0.0002 GeV, the two rays have complex weight magnitudes `|F/q|=229.77` and `264.94 pb`, much larger than the random pilot's extrema below. The pilot has not explored this rare tail adequately. These two fixed-angle directions do not establish global variance or control tangencies, soft limits or affine-star images.

## Fixed-budget pilot

Each proposal uses seeds 1337, 7331 and 424242, 2,048 draws per seed, one fresh untrained iteration and 20 cores: **24,576 physical draws**. All twelve runs finish with zero final NaN/unstable samples and no Arb rescues. The one loaded CLI batch takes **565.70 seconds**, including setup and settings changes, excluding the surrounding file-hash audits. No proposal was enlarged after seeing its result.

| Proposal | Phase | Pooled signed estimate ± SE [pb] | Pooled absolute estimate ± SE [pb] | Largest absolute weight [pb] | Recorded total / map time [ms per draw] |
|---|---|---:|---:|---:|---:|
| Six optimized LMBs | Re | 0.00041816 ± 0.000266 | 0.00072896 ± 0.000266 | 1.4549 | 73.35 / 0.075 |
| Six optimized LMBs | Im | 0.0001452 ± 0.000324 | 0.0011993 ± 0.000324 | 1.2452 | 73.35 / 0.075 |
| Direct H, power 1 | Re | -2.2908e-06 ± 1.62e-05 | 8.0266e-05 ± 1.62e-05 | 0.071117 | 72.04 / 0.156 |
| Direct H, power 1 | Im | 0.00013017 ± 0.000147 | 0.00028061 ± 0.000147 | 0.89139 | 72.04 / 0.156 |
| Direct H, power 2 | Re | 1.0273e-05 ± 1.62e-05 | 8.1616e-05 ± 1.62e-05 | 0.067649 | 79.47 / 0.151 |
| Direct H, power 2 | Im | -4.9453e-06 ± 3.46e-05 | 0.0001628 ± 3.46e-05 | 0.15213 | 79.47 / 0.151 |
| H power 2 + six LMBs | Re | 0.00035843 ± 0.000293 | 0.00066464 ± 0.000293 | 1.6971 | 72.16 / 0.170 |
| H power 2 + six LMBs | Im | 0.0003572 ± 0.000295 | 0.00096122 ± 0.000295 | 1.4525 | 72.16 / 0.170 |

Absolute monitors are `|Re I|` and `|Im I|` separately. For equal per-seed size n=2048, pooled mean mu and reported seed means mu_i/errors e_i, the empirical error shown is

```
SE² = [sum_i n(n-1)e_i² + n sum_i(mu_i-mu)²] / [(3n)(3n-1)].
```

This pools within-seed sample second moments and between-seed mean differences; it is not the standard error estimated from just three seed means, nor a finite-variance certificate. Per-seed values/errors/timing are retained in JSON. Same seeds provide pairing, not common physical points across different maps; no covariance-blind significance claim is made. Recorded per-draw time is an evaluation statistic, not end-to-end wall time.

P1/P2 imaginary errors for the three seeds are respectively `2.223e-5 / 3.734e-5`, `4.388e-4 / 9.371e-5`, and `2.860e-5 / 2.506e-5`: the pooled improvement is largely driven by one P1 excursion and is not uniform across seeds. P2's real error barely changes. Quad usage increases from 2.47% for P1 to 10.07% for P2, versus 1.22% for optimized LMBs and 2.31% for the mixture; focusing closer to the surface therefore carries a rescue cost. Some lower-precision attempts report a nonfinite process value, but the completed stability stack leaves no invalid draw.

## Exact maximum replay and reproducibility

All **48 stored signed extrema** replay exactly, with relative difference zero. Each replay uses the original complete nested `Sample` from `workspace/state/integration_state.bin`, not merely printed coordinates or the final trained grid. The sample's outer total weight is applied once; nested weights already contain their children. Historical outer weights are exactly 1, 6 and 7 for direct-H, optimized LMB and mixed proposals, with continuous weight 1 and no orientation factor 936. Actual conditional probabilities are recorded.

Native map Jacobian, partition, selected factor and raw momenta are separate diagnostics: the current physical result already includes `J_c*w_c*F`, with reported unapplied J=1, and those factors must not be multiplied twice. All maximum replays select Double under the same stack. Their stored checkpoint lacks the original precision history, so this is an exact estimator/probability replay, not an independent forced-Arb audit of all 48 extrema. The independent Arb checks are the dedicated H/Z and soft controls above.

The JSON embeds cards, diagnostic/runner sources, exact link and CLI commands, source/binary/state/input hashes, all per-seed component results, normalized-reference misses and passes, native physical events, fixed-ray fits, and complete maximum Samples plus replay diagnostics. Scratch paths are provenance: extract the sources/cards and substitute fresh state/binary/workspace paths to reproduce. Existing public state, map, precise-evaluation and integration owners perform the work; there is no diagnostic sampling engine.

The next missing ingredients remain certified CT-star pullbacks and a generic joint two-normal H/Z chart. Longer adaptive/equal-wall-time comparisons should follow those geometry and normalization gates; this direct-H pilot does not replace them.
