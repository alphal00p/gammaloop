# GL638 hosted joint gate: first optimized measurements

The original H/Z joint channel runs on the complete GL638 integrand. Its
retained numerical failure now passes after the canonical physical-center repair
documented below. The original measurements at commit
`072c9999ee928e8e9979f92bd5fb16d4f17ca15f` failed both cost and stability gates;
those results remain archived. Performance tuning has stopped at the user's
request. No integration gain or bounded-weight claim follows from these tests.

The later seven-candidate cut-matched preflight exposed a foreign Cut1 inverse
failure at an ordinary-fallback point. Its exact reproduction and narrowly
scoped prepared-ray fix are recorded in
[PHASE_SPACE_RAY_CERTIFICATION.md](PHASE_SPACE_RAY_CERTIFICATION.md).
The rebuilt exact GL638 point replay passes with unchanged native forwards.
The complete seven-candidate preflight now also passes, as recorded below.

## Completed cut-matched acceptance

Optimized build5, with Rust source equal to commit `4d192b5be`, accepts all seven
unchanged candidate catalogues. Every mode passes the summed Gaussian
normalization and raw second-moment criteria at 32,768 indexed Halton points.
The targets are 1 and 1,085,600 GeV²; the fixed limits are 0.06 absolute and
0.08 relative. Earlier finite 8/8,192-point results remain in the reports;
the eight-point check is only a smoke/parallel-consistency test. Halton
dispersion is diagnostic, not a confidence interval.

| Catalogue | Channels | Normalization | Raw-moment relative error |
| --- | ---: | ---: | ---: |
| Optimized LMB | 6 | 1.006498 | −2.9234% |
| Six cuts | 6 | 1.031806 | +3.8081% |
| Six cuts + direct joint | 7 | 1.039626 | +4.6613% |
| Six cuts + soft LMB | 7 | 1.047922 | +5.8060% |
| Six cuts + direct joint + soft LMB | 8 | 1.051237 | +6.1107% |
| Six cuts + composed Cut1 LU-h→joint | 7 | 1.040412 | +4.6479% |
| Six cuts + composed Cut1 LU-h→joint + soft LMB | 8 | 1.048940 | +5.7802% |

All 56 retained physical rows pass the 221 independent comparison checks,
including ordinary/forced-Arb, unchanged bare physics across proposals,
six-cut event sums and the selected raw `bare × canonical partition` oracle.
Acceptance follows the declared complex-norm criterion. The audit retains
28 small-component relative misses; it does not promise separate relative
accuracy for every negligible real or imaginary component. All 35 saved-file
hashes, 936 orientations, six physical cuts and CT metadata remain unchanged.

The run completes in 4,255.332 s with 20 workers and peak process RSS
63,645,268 KiB. Part of it overlaps the local diagnostic runs, so these numbers
are not a sampling-cost or scaling comparison. The
[acceptance archive](gl638_hosted_joint_gate/cut_matched_acceptance/artifact_hashes.json)
retains exact reports, the reproduced physical audit, independent acceptance
summary and build/card/source links. The original failed preflight remains
preserved in the phase-space-ray archive.

The [cut-matched local study](LOCAL_CUT_MATCHED_HZ.md) separately establishes
plateaus on two regular hard directions. MC variance, maximum-weight origins
and a supported central estimate remain pending. Their
[predeclared protocol](gl638_hosted_joint_gate/cut_matched_acceptance/mc/MC_ANALYSIS_PROTOCOL.md)
screens all seven with three seeds and 2,048 samples per seed, then freezes
the ordinary baseline, matched non-joint control and selected joint candidate
before an independent five-seed, 32,768-sample confirmation.

The [execution declaration](gl638_hosted_joint_gate/cut_matched_acceptance/mc/SCREEN_BUILD5_COMMANDS.md)
fixes the 20/30-worker choice using throughput and memory before advanced
outcomes. The chosen baseline appears once in the screen; the other worker
control is retained but excluded. The archived multi-directory analyzer keeps
the existing pooling/selection formulas and authenticates unique mode/seed,
source, card, preflight, N and chosen workers. Its compatibility checks preserve
32 historical pooled means/errors and reject duplicate input directories.
That pre-results archive contains no pilot outcomes or inferred winner. The
separate [completed screen](GL638_MC_SCREEN.md) now retains all 21 stable runs
and 43,008 draws, freezes optimized LMB/cuts/composed joint for confirmation,
and preserves the excluded worker control. Its tiny direct/composed score
separation is not a demonstrated variance improvement between those variants.
Confirmation and physical maximum attribution remain separate follow-ups.

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
GL638 CT-off was still unverified at this stage.** The later successful 1f replay
is recorded below. The archived exit-1
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

## Matched measurements after evaluator specialization

The follow-up uses commit `1f83233dd79b5c7b8b09aef38ad5acdd5fde932c`.
It includes scalar inverse reconstruction, three active eager Jacobian columns,
independent nested worker evaluator buffers, and the CT-off representative-cache
fix. The relevant 137 core checks passed, including the full 18-orientation
amplitude fixture. Exact active-three/full-seventeen Jacobian comparisons pass;
separate scalar evaluation agrees within native rounding. No physical accuracy
budget, canonical precision, certificate or density was relaxed.

The matched untraced 20-worker run uses exactly the previous complete Samples,
hard raw cases and representative maxima. Native totals, final precisions,
factors and validity flags match the previous run in both modes. All 936
orientations and six cuts remain active; the same 35 saved-state hashes are
unchanged. The measurement scope and W/S/P definitions above remain unchanged.

| Set | Optimized LMB: S/P; bound | Joint + LMB: S/P; bound |
| --- | ---: | ---: |
| Representative, all 960 calls | 3.69%; 4.32% | 39.48%; 39.76% |
| Representative, 957 valid joint calls | — | 45.95%; 46.27% |
| Selected hard points, combined | 2.32%; 3.89% | 18.08%; 18.33% |
| Representative-maximum replay | 6.14%; 6.42% | 47.51%; 47.82% |

The valid-only row remains a diagnostic, not acceptance. The same source
(worker 12/draw 12) is invalid on all three repetitions. Its slow rescues add
10.42 seconds to P; including them must not conceal the valid-call cost. The
selected joint bounds are 44.39% at hz00 and 43.33% at hz03; the selected LMB
sibling at soft22 costs 8.19% and requires Quad physics. The representative
bound p99/max is 159%/199%; maximum-replay p99/max is 91.3%. No tail is removed.

Total representative S changes only from 27.9771 to 27.8895 seconds (−0.31%).
This improvement is insufficient. Both timing clients rewarm their worker
clones, so the clone-ownership fix cannot explain this comparison. Baseline
aggregate costs remain below 10%; no baseline optimization is needed. The
next bounded optimization reuses radius-independent directed certificate
ranges within one preparation. It must preserve every radius trial, predicate,
error guard and proposal decision. No new performance result for that change
is included here.

## Physical discrepancy localized to the center choice

The new CT-off code now runs successfully on actual GL638. Separate traced
diagnostics retain the same failed Sample and verify identical canonical raw
point, J and every partition weight across settings and before/after calls.
These diagnostics are not runtime-budget measurements.

| Target and probe | Ordinary configured stack | Forced Arb |
| --- | --- | --- |
| CT on, Euler `(0.1,0.2,0.3)` | Invalid; Quad/Arb discrepancy `2.168501393203e-10` | Invalid; `2.1685013932032444e-10` |
| CT on, Pi/2 about z | Valid Quad; `5.329944811107988e-14` | Valid; `5.329945071257791e-14` |
| CT off, Euler | Valid Quad; `8.814134297321171e-23` | Valid; `1.6901990237284813e-301` |

The more detailed CT trace identifies a first internal difference in physical
cut 3 (cut group 0), on its left side. Identity and Euler have the same complete
membership `[0,1]`, parent `LmbIndex(64)` and active index 2. Their independently
solved binary64 SOCP centers have norms differing by approximately
`−1.21607e−7 GeV` (`−3.40155e−10` relative). The subsequent native alpha solves
are accurate for these different rays: their residuals are of order `1e−299`.

| Quantity in the first Arb pair | Euler minus identity, GeV |
| --- | ---: |
| Left center norm | `−1.216073698e−7` |
| Left radius, both thresholds | `−3.109727606e−8` |
| Left threshold 0, r-star | `+7.176130159e−8` |
| Left threshold 1, r-star | `+4.065094032e−8` |
| Right radius and both r-star values | Exactly equal in retained native output |

The right center is zero. The extraction interprets printed center coordinates
as round-trippable binary64 values and r/r-star as native Arb decimals. Treating
the center printouts as exact decimal values instead changes the norm difference
by only about `7e−15 GeV`; this is immaterial to the finding. The raw trace,
line-numbered extraction and script are retained.

The center repair, validated in the final section below, chooses complete physical
overlaps once from each canonical point, independently of the sampling channel,
then promote the stored center bits exactly and rotate in native precision.
Native cuts retain their own complements, LU/alpha solves and raised packets.
Raw unselected, ordinary, summed and nonmaster graph rows require the same
physical convention for the same accepted cut set. Explicit channel-ID selectors
retain the real source annotation and may intentionally change that set; their
semantics must not be replaced by a fictitious channel or a numerical retry.
Unexplained native membership disagreement must fail without selecting a
replacement center. It is not an A-star map or an amplitude covariance claim.

## Follow-up archive and replay

The new files are in `gl638_hosted_joint_gate/performance_1f/`. The `timing/`
directory is untraced 20-worker cost evidence. `diagnostic/` and `ct_trace/`
contain separate traced one-worker runs; their timings are not additive
component costs or evidence for the 10% gate. In particular, detached hosted
component calls that reject a missing runtime row remain recorded as errors.

The archive retains exact raw compressed reports/logs, failure rows, settings,
state hashes, driver source/build record, compiler invocations and the source
patch proof tying optimized build 2 to commit `1f83233dd`. Large duplicate
summary reports are represented by their uncompressed hashes; no binary or
saved physics state is copied. `artifact_hashes.additions.json` follows the
existing archive's uncompressed-hash convention.

```sh
GL_DISPLAY_FILTER=off \
  /tmp/gl638-hosted-joint-gate/drivers/gate-performance \
  /tmp/gl638-hosted-joint-gate/manifest-performance-timing.json \
  /tmp/gl638-hosted-joint-gate/results-performance-timing1 \
  optimized_lmb,joint_hz_plus_lmb timing 20 16 3 1337

GL_DISPLAY_FILTER=off,gammalooprs::integrands::process=debug \
  /tmp/gl638-hosted-joint-gate/drivers/gate-performance \
  /tmp/gl638-hosted-joint-gate/manifest-performance.json \
  /tmp/gl638-hosted-joint-gate/results-performance-diagnostic1 \
  joint_hz_plus_lmb failure-diagnostic 1 1 3 1337

GL_DISPLAY_FILTER=off,gammalooprs::integrands::process=debug,gammalooprs::subtraction::lu_counterterm=debug \
  /tmp/gl638-hosted-joint-gate/drivers/gate-performance \
  /tmp/gl638-hosted-joint-gate/manifest-performance-timing.json \
  /tmp/gl638-hosted-joint-gate/results-ct-trace1 \
  joint_hz_plus_lmb failure-diagnostic 1 1 1 1337
```

The latest user direction is to complete the current certificate optimization
and then stop performance tuning even if its cost exceeds 10%. Retain honest
cost reporting, but prioritize the center repair, actual-state Gaussian
correctness, and matched physics comparisons: ordinary sampling, simpler
nonjoint advanced channels, and the best configuration including the joint
channel. Compare equal sample counts, signed and absolute integrals with Monte
Carlo errors, maximum weights and their origin, and H/Z-corner scaling. The
remaining target is the best defensible GL638 central value and uncertainty;
up to 30 cores and 300 GB are authorized when scaling is useful. No such new
physics comparison is claimed by this archive.

These reproduce the historical invocations; use fresh output directories for a
new run. The diagnostic stage returning successfully means its requested
matrix was retained, not that every physical row passed. Full GL638
normalization, integration gain and bounded weights remain unestablished.

## Final certificate round: performance work stopped

The [matched comparison](gl638_hosted_joint_gate/certificate_round/comparison_1f_certificate.md)
records optimized build 3: base `1f83233dd`, with the archived source patch
`61e9a9b07a858f152537dbc98d2778cbb05071558fdb8276b967d7427af4f06d`.
Radius-independent directed ranges are reused within one preparation. Radius
trials, predicate ordering, lower-bound arithmetic and underflow guards remain
unchanged. All eight focused tests, all-target checks and changed-line Clippy
checks pass. The newly added reuse fixture was corrected to include its known
admitted radius; no existing physical tolerance changed.

| Joint set | Previous conservative overhead | Final conservative overhead |
| --- | ---: | ---: |
| All 960 representative calls, including 3 invalid | 39.76% | 35.67% |
| 957 valid representative calls, diagnostic subset | 46.27% | 41.79% |
| Selected hard points, combined | 18.33% | 17.05% |
| Representative maxima | 47.82% | 48.10% |

Aggregate representative sampling time falls by 15.25% in this batch. The
unchanged baseline timing also varies; this is not a universal speed claim.
Every original Sample, hard case, retained maximum, native total, precision,
factor and stability outcome matches exactly. The same three physical failures
remain, and all 35 state hashes are unchanged. The archive retains individual
hard-point ratios, tails and the complete failed rows.

The user has explicitly accepted present performance. Further runtime tuning
stops despite exceeding the former 10% target. Work proceeds to the requested
matched physics comparisons after the physical-center validation below.
Neither bounded weights nor a variance improvement is established by this
performance experiment.

## Canonical physical centers: actual GL638 failure resolved

Optimized build 4 uses base `6e9bf401db75cb13ec20e18377a444bc2fcd4f33` plus
source patch `1d78e98b73e6a8dafaf6e446353c2a1607dc0724f8ea7d1c2d88c4962b81bba7`.
It retains complete physical overlaps from the canonical Arb source, then
promotes the binary64 center bits exactly before each native rotation. The
[implementation note](CANONICAL_CT_CENTERS.md) describes accepted-cut validation,
explicit channel selectors, native complements and raised packets. Fourteen
focused and 27 broader core tests pass, including the full hosted 8192-draw
fixture; checking and all-target Clippy pass with no changed-line warnings.

The same worker12/draw12 Sample now completes every requested physical branch.
Each has one ordinary call and one forced-Arb call, with all 936 orientations
and six retained cut weights. Tolerances, source cube and outer weight 7 are
unchanged. The values below are the actual final two-probe discrepancies.

| Target and second probe | Ordinary configured stack | Forced Arb |
| --- | --- | --- |
| CT on, Euler `(0.1,0.2,0.3)` | Valid Quad; `2.200848788466276e-21` | Valid; `4.2692671432208616e-290` |
| CT on, Pi/2 about z | Valid Quad; `2.346417869593957e-26` | Valid; zero |
| CT off, Euler | Valid Quad; `8.814134297321171e-23` | Valid; `1.6901990237284813e-301` |
| Gaussian reference, either Euler branch | Valid Double; `5.001898468563083e-16` | Not requested |
| Gaussian reference, Pi/2 about z | Valid Double; zero | Not requested |

All six physical results are finite and valid, with `is_nan=false`; Double
still fails before the physical body on the retained hard sample, then Quad
recovers. The valid control remains valid Double in both clone-only and
rewarmed workers. Gaussian results are checked f64 reports after native
evaluation, and these single points are not a normalization test.

Euler and Pi2Z return exactly equal native totals and all six cut weights within
each precision. The ordinary-versus-Arb complex-norm discrepancy is
`1.77734e-21` with CTs and `1.64143e-22` without. The tiny CT-off real part does
not have meaningful relative component accuracy; the archived audit retains
its absolute and norm-scaled discrepancies. Every canonical point, Jacobian and
partition comparison is exact across branches, workers and before/after calls.
CT-off forced-Arb totals and cut weights also equal the previous 1f result
exactly. The new CT-on identity value changes by `1.70996e-17` in complex norm
from the previous identity prescription; it is not claimed bit-identical.

The formerly differing physical cut3/group0 retains center
`(-96.7415997053086,234.15172244667613,252.2392543658604)`, membership `[0,1]`
and parent `LmbIndex(64)`. The table now displays the unrotated canonical center.
In the actual Arb Euler consumption, both left radii agree exactly and both
r-star differences are only `-4.8e-299 GeV`; the alpha differences are about
`-2.3e-302`, with energy residuals at most `1.92e-298`. Right radii/r-star values
agree exactly. Pi2Z radii, r-star and alpha agree exactly. The display trace
does not expose `file.active_center`, so equality of the stored center table
is not presented as a direct measurement of the rotated native vector.

The [archive](gl638_hosted_joint_gate/canonical_centers/run.json) stores the
unaltered progress JSON and trace compressed, all 35 unchanged state hashes,
cards, exact source patch, build/dependency provenance, gate logs and repeatable
analysis. Driver source is the same archived `7b5e474e…` source; no binary or
physics state is copied. Duplicate full summaries are represented by hashes.
The `canonical_physical_preparation_time` subset of P is about 3.1–3.3 ms in
these CT-on calls and exactly zero for CT-off and reference calls. This is a
traced one-worker diagnostic, not a cost or parallel-scaling measurement;
rejected detached hosted component calls remain unusable timing observations.

```sh
RAYON_NUM_THREADS=20 \
GL_DISPLAY_FILTER=off,gammalooprs::integrands::process=debug,gammalooprs::subtraction::lu_counterterm=debug \
  /tmp/gl638-hosted-joint-gate/drivers/gate-centers \
  /tmp/gl638-hosted-joint-gate/manifest-centers-diagnostic.json \
  /tmp/gl638-hosted-joint-gate/results-centers-diagnostic1 \
  joint_hz_plus_lmb failure-diagnostic 1 1 1 1337
```

The process exits zero and the complete matrix is valid. This resolves the
retained GL638 rotation failure; global numerical stability, CT-star coverage,
full-state normalization, H/Z boundedness and integration improvement still
require their own physical evidence.
