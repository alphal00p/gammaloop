# Additional GL638 channel for the recovered real maximum

The existing 14-channel configuration combines the six optimized LMBs, six
Cutkosky-cut maps, the compact Cut-1 H/Z map, and the Cut-3 right-threshold map.
The new experiment retains these and adds one channel built entirely from the
existing conditional surface maps. It changes sampling only: all physical cuts,
threshold counterterms, subspaces, centers, and multiplier functions are retained.

## Why a further channel is useful

The independently replayed maximum has real weight
`-310.7628500264454 pb`. Its dominant terms are the Cut-1 integrated counterterms
for `A=[7,8]` and `Z=[3,10,13]`. The original sample lies outside the compact H/Z
chart's certified domain. At the Cut-1 prepared point (`sqrt(s)=600 GeV`),

| Quantity | Value (GeV) |
|---|---:|
| A = E7 + E8 - Q | 9.6302203511 |
| H = E2 + E4 + E12 - Q | 5.9596059591 |
| Z = E3 + E10 + E13 - Q | 2.7240141897 |
| Magnitude of q7 | 250.9646417154 |
| A-shell radius, sqrt((Q/2)^2 - MT^2) | 245.0938595722 |
| Radial distance from A | 5.8707821432 |

These are finite off-shell distances. The detailed decomposition and the
normal/1000-bit agreement are in
[the maximum replay](GL638_OPTIMIZED_LMB_SPIKE_20260917.md).

## Valid independent fibers

Use the parent LMB `[3,6,7,10]`, denoting its prepared vectors by `(p,a,s,v)`.
The host cut `[2,6,10]` depends only on `(a,v)`. A depends only on `s`, while H
and Z depend on `p` and the fixed host vectors. Therefore an ordered composition
of the host six-dimensional phase-space map, a three-dimensional A map, and a
three-dimensional H **or** Z map is valid. H and Z cannot be put in separate
independent `p` blocks because they share those three variables.

In particular, A is cylindrical in the native counterterm's six-dimensional
active space `[3,7]`; its sampling block is `[7]` with `[3]` treated separately.
A sampling chart need not reproduce the counterterm's projection subspace.
This change does not modify the native two-loop A counterterm.

The direct surface maps have normalized whole-fiber support, with the existing
ordinary-map fallback if a conditional surface does not exist. They avoid the
compact H/Z chart's support limitation. The LU profile is unchanged:
`{kind="lu_h", approximation="log_logistic", broad_fraction=0.02}`.

Three candidates are compared at the exact canonical maximum: host+A,
host+A+Z, and host+A+H. The latter two use no ordinary complement block because
the three blocks cover all twelve spatial coordinates.

## Selected channel and measured improvement

Choose **Cut 1 + A + Z**, whose two threshold fibers address the two dominant
real counterterms at this point:

```toml
[sampling.channel_definitions.GL638.cut1_left_7_8_and_Z]
around = "then(block(lmb(6,10),phase_space(cut(2,6,10))),block(lmb(7),at_cut(cut(2,6,10),surface(7,8))),block(lmb(3),at_cut(cut(2,6,10),surface(3,10,13))))"
parent_lmb = [3,6,7,10]
subspace_lmb = []
on_cut = [1]
radial_profile = { kind = "lu_h", approximation = "log_logistic", broad_fraction = 0.02 }
```

Append this name to the original selection. The catalogue then has fifteen
channels; the new channel's index is 14. The global radial power remains 2,
with `b=0.3`, so dimensional map scales continue to be relative to `Q`.

| Added channel | Density at the spike (GeV^-12) | Raw denominator gain | Fresh uniform-mixture gain | Inverse/forward coordinate error (GeV) |
|---|---:|---:|---:|---:|
| Cut 1 + A | 2.44241459e-37 | 1.46357 | 1.36599 | 5.68e-14 |
| **Cut 1 + A + Z** | **3.83903319e-36** | **8.28641** | **7.73398** | **5.68e-14** |
| Cut 1 + A + H | 1.26358201e-36 | 3.39825 | 3.17170 | 7.11e-14 |

All three support the exact original canonical point. The original fourteen
channels have summed density `5.268761564513515e-37 GeV^-12`; the original six
optimized-LMB channels have sum `4.072616472396053e-37 GeV^-12`. Compared with
a fresh uniform six-LMB mixture, the selected fifteen-channel mixture improves
the density at this point by **4.28806**.

The diagnostic reconstructs the canonical point using the original source
precision and requires exact agreement with the retained native replay before
evaluating inverse densities in Arb precision. It then rounds each inverse
coordinate to binary64 and checks its forward reconstruction. The errors above
compare binary64 momentum coordinates, so they describe that practical replay
accuracy rather than an exact native-arithmetic identity.

The chosen channel's inverse cube is:

```text
0.1251959013471031  0.6720992727031357  0.6683068027704575
0.3428091881658814  0.1868849890763291  0.7452067107805721
0.6413376292351303  0.6899232813031906  0.7296858295826409
0.4666995950359920  0.14642355954146338 0.8212867782367871
```

The retained unparameterized physical integrand predicts a fresh uniform-mixture
real weight of `-2.728128075 pb` for this point. This must **not** be presented as
a like-for-like reduction of the historical trained-grid weight `-310.76 pb`:
that old weight also contains its historical adaptive proposal correction.

A read-only production CLI replay of the chosen fifteen-channel runtime card
confirms this independently. Channel 14 with the cube above returns
`-0.18187520498453386 - 1.9882685737071093 i` with unit outer weight. Multiplying
by the fresh uniform-channel factor 15 gives
`-2.728128074768008 - 29.82402860560664 i`, agreeing with the raw-density
prediction to about `6.5e-14` relatively. The final reported result is finite.

The generation/backend display for the reused state identifies its actual
backend as **Eager**, with compilation disabled and no shared library. This
replay and the matching run card do not claim the historical SymJIT O3 setup.
The output, commands, and exact CLI binary hash are included in the evidence
archive. The authoritative backend evidence is the generation table in
[the complete CLI log](gl638_spike_channel_20260917/inspect.log.gz), recorded at
2026-09-17 10:15:06 UTC; it reports `compile enabled=false`,
`generation backend=eager`, `active f64 backend=eager`, and `.so size=0 B`.
A card's inactive compilation options do not override this frozen state.

[The exact invocation](gl638_spike_channel_20260917/inspect_argv.json) uses
`--read-only-state --no-save-state` and runs
[the retained command block](gl638_spike_channel_20260917/inspect_commands.toml):
display the loaded backend, apply `selected15.toml` in memory, then inspect
graph 0/channel 14. The process exited successfully. The
[production JSON](gl638_spike_channel_20260917/selected15_inspect.json) preserves
the unit outer weight and unit reported parameterization Jacobian; its precision
metadata is not retained. The
[verification record](gl638_spike_channel_20260917/production_replay_verification.json)
records the exact prediction comparison, completion time, and SHA-256 hashes.
Its componentwise relative discrepancies are `6.92e-14` (Re) and `6.54e-14`
(Im). These artifacts were checked byte-for-byte against the diagnostic outputs;
the compressed log was checked after decompression.

The [diagnostic source and exact outputs](gl638_spike_channel_20260917/) retain
the three-candidate runtime card, all inverse coordinates and support decisions,
the selected fifteen-channel card, the native CLI replay, and the comparison
summary. The diagnostic three-candidate catalogue has
seventeen channels; each table row above correctly adds **only that candidate**
to the original fourteen, rather than using the sum of all three candidates.

## Interpretation of the pointwise comparison

Let `S14` be the sum of the original fourteen unadapted map densities and `q`
the added channel density. The raw partition denominator increases by
`1+q/S14`. A fresh uniform channel mixture instead improves the density by
`(14/15)*(1+q/S14)`; the additional channel probability must be included.
Neither number is a prediction for the existing trained Havana grid, which has
different discrete probabilities and continuous densities. A new workspace is
required for the changed channel catalogue.

The probe records support, inverse cube coordinates, forward reconstruction,
and the full raw density at the *same physical point*. A short physical pilot
checks production evaluation before the long run; its outcome belongs to the
run report rather than this density-only measurement.

This is direct-threshold focusing. It is not an exact map of the H/Z equations
evaluated at the native A-counterterm star. The latter requires the projection
center and branch as part of a currently unavailable CT-star map. A density gain
at this maximum therefore does not establish globally bounded weights or finite
variance.

## Separate production comparison

The run card is [`gl638_advanced_sampling_max_weight.toml`](../gl638_advanced_sampling_max_weight.toml).
It uses a separate state `gammaloop_state_GL638_advanced_sampling_max_weight`
and a fresh `workspaces/GL638_advanced_sampling_max_weight_workspace`, both
under the NNLO example directory. The generated all-orientation state is copied
from the previous advanced-plus-optimized-LMB run; no generation setting or
physical counterterm changes. The old three workspaces remain intact.

Production keeps 50 workers, 81,920 samples per iteration, a 10-billion-sample
limit, seed 62001, both learning rates 0.25, real-part training and stability,
and per-iteration reporting of Re, Im, |Re| and |Im|. The energy is 600 GeV,
`mu_r=91.188 GeV`, `m_uv=50 GeV`, with `min_abs_wgt_for_escalation=1e-12`.
The frozen evaluator backend is Eager, as verified above; the card's O3 compile
setting is inactive while compilation is disabled.

A separate 128-point, four-worker pilot completed with no unstable or NaN
samples and wrote its iteration JSON. This checks startup, channel binding,
physical evaluation and output, rather than providing a statistical comparison.
The production workspace starts from a fresh grid, without that pilot's points.

From the repository root, with the existing state:

```bash
target/dev-optim/gammaloop -n --read-only-state \
  examples/cli/epem_a_ttxh/NNLO/gl638_advanced_sampling_max_weight.toml \
  run integrate -D renderer=tabled -c quit
```

To generate the state from scratch, use the card's `run all` workflow without
`--read-only-state`. From `examples/cli/epem_a_ttxh/NNLO`, the monitor discovers
the fourth column automatically once the new workspace manifest is written:

```bash
.venv/bin/python monitor_tth_nnlo_integration.py --graph gl638
```

Ordinary iteration summaries are also retained in
`gl638_advanced_sampling_max_weight.long.log`, which the same command can read
with `--live`.
