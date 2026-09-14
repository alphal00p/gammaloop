# GL638: compressed SymJIT and sampling timing

At the existing **√s = 1000 GeV**, **μr = mUV = 91.188 GeV** point,
compressed SymJIT O2 reached **95.415 samples/s with 20 actual workers**, versus
79.970 samples/s for Eager and 95.287 samples/s for uncompressed SymJIT O2 on
the same frozen samples. The compressed profile's measured increase over Eager
is 19.3%; its 0.13% throughput difference from uncompressed is too small to
establish a compression advantage from these runs. Sampling remains the main cost:
38.83 ms per sample versus 27.69 ms for the physical integrand. The earlier goal
of sampling below 10% of physical evaluation time is **not met**.

This benchmark retains the old seven-channel catalogue: six standalone Cutkosky
cut channels plus the composed Cut1 joint H/Z channel. It does not yet include
the additional soft LMB or Cut3→right[5,10] channel. The requested subsequent
600 GeV study, fixed μr = 91.188 GeV and mUV scan including 50 GeV is a separate
physics comparison; none of those changes enter these timings.

## Matched setup and source identity

Each profile replays the same 512 complete `Sample` records three times, at one
and then 20 workers: 1536 timed evaluations per worker-count setting. The deck
was generated once from the initial grid, with seed 1337 and without adaptation.
It preserves every outer/discrete sampling weight. All timed calls passed in
Double precision with the frozen ordinary stability stack; the benchmark does
not estimate rescue rates in more difficult corners.

The read-only state is `epem_a_tth@NNLO`, GL638, all 936 production orientation
keys, six cuts and 19 threshold variants, with the full original local 3D and
integrated UV subtraction. All 35 saved-state hashes agree before and after each
profile. Changing the backend occurs only in memory.

All four profiles use **build8**, base `33156a1c97f21b4c7569cc9f01837221b3e82a3f`
plus frozen source patch
`7e6d27b3673e92bb9a2c87b726b95db85c4adecd0fb2838529b773c2f8824799`.
The source preparation and full sampling weights remain **Arb1000**. Directed
joint support certification starts at 128 bits and escalates to 2048 only when
the same geometric predicate is inconclusive. Inclusive evaluator timers now
cover the actual wrappers, including stability probes/rescues and IFT-alpha
calls. Later fixed-Quad source-policy work is **not measured or validated** here.

The dependency is SymJIT **2.25.6** with the existing Symbolica fork
`a277a8dadbaae5a8d6605f2f381991a601e1f345`. The two complete standard SymJIT
configurations enable O2 and differ only in `compress = false` versus `true`,
while preserving `compact = false` and `direct = false`; compression and
compaction are different options. A fourth profile raises the compressed
configuration's stack limit. Quad/Arb evaluators remain Eager.

## Timings

The wall column is dispatch elapsed time for all 1536 evaluations. S and P are
inclusive means of worker-local elapsed durations per sample; E is a **subset of
P**, not an additional cost. These worker-local times cannot be divided directly
into dispatch throughput, particularly with concurrent execution and imbalance.

| Evaluator profile | Workers | Dispatch wall [s] | Samples/s | Sampling S [ms] | Physical P [ms] | Evaluators E [ms] | Outer overhead [ms] |
|---|---:|---:|---:|---:|---:|---:|---:|
| Eager | 1 | 126.004 | 12.190 | 42.049 | 38.829 | 34.992 | 0.443 |
| Eager | 20 | 19.207 | 79.970 | 39.754 | 55.993 | 47.958 | 1.880 |
| SymJIT O2 uncompressed, 1 MiB stack limit | 1 | 106.401 | 14.436 | 41.717 | 26.743 | 22.857 | 0.450 |
| SymJIT O2 uncompressed, 1 MiB stack limit | 20 | 16.120 | 95.287 | 41.905 | 33.525 | 27.592 | 1.097 |
| SymJIT O2 compressed, 1 MiB stack limit | 1 | 106.783 | 14.384 | 42.471 | 26.196 | 21.972 | 0.467 |
| SymJIT O2 compressed, 1 MiB stack limit | 20 | 16.098 | **95.415** | 38.832 | 27.691 | 22.294 | 0.957 |
| SymJIT O2 compressed, 16 MiB stack limit | 1 | 108.183 | 14.198 | 42.771 | 26.823 | 22.533 | 0.472 |
| SymJIT O2 compressed, 16 MiB stack limit | 20 | 18.461 | 83.205 | 55.371 | 41.472 | 32.484 | 2.835 |

For the fastest measured profile (compressed, normal stack limit, 20 workers),
the following components are **disjoint** and sum to 67.480 ms of worker time
per sample:

| Component | Mean [ms/sample] |
|---|---:|
| Canonical sampling preparation, C_S | 38.401 |
| Remaining sampling, S − C_S | 0.431 |
| Canonical physical preparation, C_P | 2.514 |
| All physical evaluator wrappers, E | 22.294 |
| Other physical work, P − C_P − E − events | 2.882 |
| Event processing | 0.000 |
| Outer wrapper work, worker elapsed − S − P | 0.957 |

The raw nanosecond durations and corresponding decomposition for every profile
and worker count are retained in the analysis. Subset inequalities are checked
before subtraction; no negative residual is silently clamped.

With 20 workers, Eager's evaluator time rises from 34.99 to 47.96 ms/sample;
uncompressed SymJIT rises from 22.86 to 27.59 ms/sample, while compressed SymJIT
remains near 22 ms/sample. This is compatible with the user's
cache-pressure concern, but these measurements do not isolate cache misses or
establish their cause. Both 1 MiB SymJIT profiles log ten SIMD requests downgraded
to scalar because of the stack limit. Raising the limit to 16 MiB removes those
warnings but is slower in this run. Neither configuration nor absence of warnings
establishes that every physical call executes SIMD instructions.

Load, activation and preparation are excluded from dispatch timing and retained
separately:

| Profile | State load [s] | Backend activation [s] | Map warmup [s] | Clone preparation/warmup: 1 / 20 workers [s] | Whole process [s] | Peak RSS [GiB] |
|---|---:|---:|---:|---:|---:|---:|
| Eager | 51.991 | <0.001 | 0.050 | 1.487 / 37.369 | 282.405 | 59.863 |
| Uncompressed, 1 MiB | 48.026 | 10.061 | 0.050 | 1.028 / 35.846 | 259.536 | 60.587 |
| Compressed, 1 MiB | 49.708 | 11.300 | 0.050 | 0.941 / 35.236 | 261.559 | 60.587 |
| Compressed, 16 MiB | 50.395 | 11.518 | 0.050 | 0.917 / 35.927 | 269.334 | 60.594 |

Each clone receives two untimed calls. Whole-process measurements also include
controls, auditing, serialization and destruction. Runs were sequential and were
the only owned GammaLoop workload; unrelated shared-host workloads remained.
CPU affinity and before/after process snapshots are archived. Three replays of
one deck are not repeated independent throughput trials: no throughput confidence
interval, linear 20-core scaling claim or MC variance improvement is inferred.

## Numerical checks and the real-part limitation

The complete Eager/uncompressed/compressed comparison passes **28,632 checks and
8704 complete-estimator comparisons**. The earlier Eager/compressed analysis
(19,090 checks and 5632 comparisons) and the separate 16 MiB comparison
(19,093 checks and 5632 comparisons) are also retained. These analyses reuse the
same reference and are not independent samples. Each profile first replays two
retained difficult native points in ordinary and forced-Arb modes. Canonical momenta, Jacobian and
partition agree exactly with their controls, and total/six-cut norms satisfy the
frozen physical acceptance criteria.

Passing those criteria does **not** certify relative accuracy of a nearly zero
real component. The complete comparison preserves twelve repeated real-relative
misses at source sample 64, six for each SymJIT profile: Eager gives −4.484×10⁻⁴⁰
and both uncompressed and compressed SymJIT give −1.183×10⁻³⁹, while
Im = −6.720×10⁻²⁵. The absolute difference,
7.349×10⁻⁴⁰, is below the configured complex-norm allowance 1.344×10⁻³⁰,
but the real component fails its own relative criterion. The native-control
audits similarly preserve four real-relative misses per profile for two
near-zero Cut2 contributions, counted in two comparisons each: ordinary residuals
−9.641×10⁻²⁰ and 5.420×10⁻²² versus Arb residuals about −1.723×10⁻³⁰⁶
and 2.783×10⁻³⁰⁸. Exact records and phase labels are retained.

The user requested ranking physics channels using Re and |Re| only. Therefore
these are qualified backend correctness checks under the **existing complex-norm
contract**, not a solution to the separate real-only stability issue. No new
cross-section estimate, bounded-weight claim or variance claim comes from this
fixed-deck benchmark.

## Source gates and retained evidence

Formatting and core/API checks passed. Focused tests covered adaptive joint
certification, conditional cut sampling, shared-group weights, evaluator handling
and timing accounting: 21 unique tests ultimately passed. The first compressed
run passed 20/21; its new timing fixture had inconsistent inclusive totals. A
routine fixture correction preserved the accounting assertions, then four tests
passed under compression-off (three overlap the earlier passing set). The failed
log is retained. Clippy reported 52 existing warnings and none with changed
primary locations. These gates predate later fixed-Quad changes.

The [runtime archive](gl638_hosted_joint_gate/runtime_build8/provenance.json)
contains 87 logical records, including frozen inputs, all result/control/progress
records, all three full analyses, configs, exact source/dependency/build records,
the two failed scratch-client link
attempts and successful third attempt, source gates, process logs and hash checks.
Original bytes are losslessly gzipped, with exact earlier artifacts reused by
hash; no executable, compiled library, saved-state payload or license value is
included. [artifact_hashes.json](gl638_hosted_joint_gate/runtime_build8/artifact_hashes.json)
maps each original logical file to its retained bytes and SHA256.
