# GL638 hosted joint gate: first optimized measurements

The original H/Z joint channel runs on the complete GL638 integrand, but it
does **not** yet meet the sampling-cost or numerical-stability gates. These
measurements concern commit `072c9999ee928e8e9979f92bd5fb16d4f17ca15f` on
2026-09-14. They do not establish an integration gain or bounded weights.

## Fixed physical calculation and proposal

The saved state is `/common/dev/gl638_checkpoint_06d470409/all/state`:
`epem_a_tth@NNLO`, all 936 orientations, all six cuts, the existing threshold
metadata and WH/WF multipliers, direct 3D local UV with orientation localization,
and integrated UV. The 35 saved-file hashes were unchanged after both runs.
Cards differ only in sampling. The optimized-LMB control has six channels;
the candidate adds one joint channel to those same six:

```toml
around = "then(complement(6,7,10),block(lmb(3),at_cut(cut(2,6,10),intersect(surface(2,4,12),surface(3,10,13)))))"
parent_lmb = [3,6,7,10]
on_cut = [1]
```

The two original equations are H=`surface(2,4,12)` and
Z=`surface(3,10,13)`. The host is cut 1, `cut(2,6,10)`; selecting it does not
restrict the physical cut sum. The joint kernel uses alpha 1 and a trial maximum
radius of 300 GeV, subject to the actual certified compact-disk decision. It is
a direct joint chart; no CT-star target has been substituted. Ordinary LMB
siblings preserve full-volume coverage. The generation parent remains
`[3,4,7,10]`, distinct from this channel's requested parent.

One canonical Arb1000 map point, Jacobian and partition are retained through
physical rescue. The physical stack is Double/Quad/Arb with requested relative
accuracies `1e-6`/`1e-10`/`1e-12`, norm-based checking and the configured Euler
rotation `(0.1,0.2,0.3)`. No precision or acceptance tolerance was relaxed.

## Correctness established so far

The preceding generic gate passed 216 unique core tests, including the hosted
8192-point Gaussian/moment acceptance and the complete 18-orientation two-loop
amplitude fixture. These do not replace GL638 acceptance.

The GL638 smoke retained 16 native replay results and passed all 53 independent
physical checks. Baseline and candidate bare results agree exactly at matching
native precision. Ordinary versus forced-Arb totals satisfy the configured
accuracy requirements. Selected raw-input values agree with bare physics times
the independently evaluated canonical partition weight, including all six cut
events. Eight relative component discrepancies concern negligible real parts;
the retained norm-scaled component checks pass without loosening their budgets.

The smoke later stopped at a scratch-client channel lookup based on a display
label. The client now uses the existing catalogue kind query; both the failed
smoke and the corrected subsequent timing run are retained. Eight finite
Gaussian smoke draws are **not** a normalization test. Full GL638 normalization
and integration pilots have not run for this candidate.

## Warmed cost at 20 workers

The optimized build uses core/API optimization level 2 and dependency level 3.
Each worker is warmed separately. Per card, the measured set contains 320
distinct representative draws repeated three times, 180 selected hard calls,
180 bare hard calls and 60 representative-maximum replays. The original nested
Samples, probabilities and repeated outputs are retained. These are fresh
iteration-1 conditions (`max_eval=0`), not trained-grid histories.

Let W be elapsed precise-call time, S the sampling counter and P the disjoint
physical-body counter across actual probes and retries. The conservative
sampling bound is `sum(max(W-P,0))/sum(P)`; the table also reports `sum(S)/sum(P)`.
Grid/RNG work lies outside W. S includes foreign inverses and the partition;
forcing Arb for an accuracy check is not included in these timings. Aggregated
worker elapsed durations are not CPU time or batch wall time.

| Set | Optimized LMB: S/P; bound | Joint + LMB: S/P; bound | Joint status |
| --- | ---: | ---: | --- |
| Representative, all 960 calls | 3.89%; 4.45% | 40.93%; 41.25% | Three invalid calls retained |
| Representative, 957 valid joint calls | — | 47.45%; 47.82% | Diagnostic subset only |
| Selected hard points, combined | 2.33%; 4.38% | 18.65%; 18.87% | Valid; fails cost |
| Representative-maximum replay | 6.22%; 6.52% | 58.20%; 58.48% | Valid; fails cost |

The selected hard-point bound is 45.70% at hz00 and 48.80% at hz03. Soft22,
selected through an LMB sibling and requiring Quad physics, has a 7.68% bound;
its slower physical body dilutes the combined ratio. The baseline's corresponding
per-point averages are all below 10%. Per-call tails remain in the raw reports:
joint representative p99/max bounds are 149%/191%, and maximum-replay p99/max
are 553%. The baseline also has isolated timing tails; no tail was discarded.

Representative mean S is stable across repetitions (29.15, 29.16, 29.12 ms),
so lazy first-use warmup does not explain the joint cost. Both ordinary LMB
generators and the joint generator pay for the foreign joint inverse. The
current driver rewarms cloned workers, so a separately identified nested shared
evaluator mutex cannot explain this run, although its ownership must be fixed
before clone-only production integration. Channel identity alone does not prove
which internal compact/fallback branch was visited.

## Retained numerical failure

One source, worker 12/draw 12 on joint channel 0, fails identically on all three
repetitions. Its complete Sample is in the artifact directory. Double fails
before a physical probe. Quad and Arb each complete two probes, but report a
relative discrepancy of approximately `2.168501393203e-10`, above their
requested tolerances. The finite returned result is marked invalid; the slow
Arb retry must not be used to make the sampling-cost denominator look better.

The source point, Jacobian and partition remain fixed through the higher lanes,
and their native normal-accuracy checks pass. This excludes a proposal redraw
as the explanation. The following diagnostic varies CT activation, rotation
and the Gaussian reference independently at this exact source. A public
binary64 raw replay remains a secondary oracle, not an identical replacement
for the canonical Arb source.

## Same-source rotation diagnostic on frozen 072

The [diagnostic archive](gl638_hosted_joint_gate/diagnostics_072/run.json)
retains worker 12/draw 12's original nested Sample, including its binary64 cube
and outer weight 7. The public bridge checks promote those exact binary64
coordinates to Arb1000. The physical calls use the original Sample through
the existing canonical-source and stability owners. All 936 orientations and
six cuts remain active. No tolerance changed. Each physical case has three
ordinary calls and one forced-Arb call; each Gaussian case has one ordinary call.

| Target and second probe | Ordinary stack | Forced Arb |
| --- | --- | --- |
| CT on, Euler `(0.1,0.2,0.3)` | Invalid; Quad and Arb discrepancy `2.168501393203e-10` | Invalid; `2.1685013932032444e-10` |
| CT on, Pi/2 about z | Valid in Quad; `5.329944811107988e-14` | Valid; `5.329945071257791e-14` |
| Gaussian reference, Euler | Valid in Double; `5.001898468563083e-16` | Not requested |
| Gaussian reference, Pi/2 about z | Valid in Double; zero reported discrepancy | Not requested |

The metric is the existing norm discrepancy across the two probes. Gaussian
rows are single-point diagnostics, reported through the checked f64 reference
API; they do not establish normalization.

Native comparisons of the Arb raw point, Jacobian and every partition weight
are exactly equal across these settings and before/after evaluation. The
forced-Arb identity totals and all six identity cut weights also have identical
retained full native decimal values. The first differing boundary is therefore
the **second physical stability probe**. The Euler failure also reproduces in
both a warmed clone and a clone explicitly rewarmed; the valid joint control
(worker 0/draw 9) remains valid in both. These checks do not yet distinguish
threshold-center geometry from other physical rotation-dependent operations.

CT-off ran last and hit an empty-cache panic before returning a physical
result. The frozen 072 code collects counterterm representative samples even
when CT-off has deliberately left their cache empty. The source fix moves
that collection into the existing CT-enabled branch. Its generated conditional
cut regression has passed, including actual CT-off direct and selected calls,
per-cut original-weight comparisons, and rest/boosted parent frames. **Actual
GL638 CT-off remains unverified until the new build runs.** The archived exit-1
status and panic are preserved; they are not a CT-off numerical result.

The diagnostic also contains one-worker detached map timings. Direct hosted
channel forwards and joint inverses reject the missing runtime row context;
their elapsed times are unusable as completed component timings. Complete
bridge forwards, complete partitions and ordinary LMB inverses succeed, but
those operations overlap and are **nonadditive**. They provide no new 20-worker
cost or 10% budget claim.

`diagnostics_072/` stores the exact progress JSON and driver source compressed,
the exact build/library hashes, the production settings, log, and unchanged
35-file before/after state hashes. Its compact `run.json` preserves the run
error and provenance without duplicating the large repeated summary inventory.
Driver source SHA256 is
`7b5e474ea21b0dd20fa12181cd8b8dc635c7d5db0d17d6400d1e1c47b406123e`.
No binary is archived. The reproducible invocation was:

```sh
/tmp/gl638-hosted-joint-gate/drivers/gate-diagnostic-z-072 \
  /tmp/gl638-hosted-joint-gate/manifest.json \
  /tmp/gl638-hosted-joint-gate/results-failure-diagnostic2 \
  joint_hz_plus_lmb failure-diagnostic 1 1 3 1337
```

## Follow-up and reproducibility

Optimize the existing evaluator owner first: compile only the three active
joint Jacobian columns and evaluate scalar reconstruction through a scalar
program during inverses. Preserve the canonical precision, original equations,
all certificates and acceptance budgets. Separately ensure cloned warmed maps
own their mutable evaluator buffers. Profile certificate preparation only if
needed afterward. Stop runtime optimization once the matched warmed 10% gate
passes. Resolve the physical stability failure before ordinary integration.

The [artifact directory](gl638_hosted_joint_gate/) contains the exact timing
client, cards, build hashes, failure Sample, state hashes, compressed raw timing
and smoke reports, and independent Decimal audit. Gzip files contain the
unaltered source bytes; `artifact_hashes.json` records their uncompressed hashes.
The archived manifest is the input used for these runs and retains some
pre-execution status text; this report supplies the execution outcome. The
driver build record identifies the exact libraries rather than selecting
arbitrary cached rlibs. Reports are outside the saved state.

The completed timing invocation was:

```sh
RAYON_NUM_THREADS=20 /tmp/gl638-hosted-joint-gate/drivers/gate \
  /tmp/gl638-hosted-joint-gate/manifest.json \
  /tmp/gl638-hosted-joint-gate/results-timing1 \
  optimized_lmb,joint_hz_plus_lmb timing 20 16 3 1337
```

Use a new output directory for a replay. The ordinary full reference and pilot
stages remain pending. Further CT-star maps and automatic channel discovery
remain separate incomplete milestones in `ADVANCED_SAMPLING_PLAN.md`.
