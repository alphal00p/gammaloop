# Local 4D UV validation and performance

The implementation canonicalizes completed hard UV denominators before source
reconstruction, merges compatible forest terms and reuses dynamically prepared
numerator mappings. Its signed source certificates and CFF allocation are
described in [the architecture](architecture-current.md). The original plan and
dated implementation history are preserved in
[SPEED_UP_UV_CTS_FROM_4D.md](../../SPEED_UP_UV_CTS_FROM_4D.md).

On 2026-09-13 the user explicitly deferred GL262 and requested merge acceptance
for GL00/GL01, the scalar matrix and the remaining correctness checks. GL262's
test and diagnostic evidence remain available; it is not counted as passing.

The major physical GL00/GL01 slowdown is resolved at the measured revision:
complete generation is within 2.53% of erased 3D and numerical sampling is faster.
This does not resolve every broader scalar case. GL21 still has a measured
numerator-dependent regression with integrated and threshold counterterms enabled,
described separately below.

## Measurement boundary

The numerical implementation measured here is
`f938f9c1386eb2670a0a7f603f5eb0878086352a`. Initial merge cleanup changed tuple
type aliases, comments and test iteration syntax. A later correction to shared
Spenso odd tensor powers changes production execution and requires fresh
measurements. The results below remain evidence for the earlier numerical source;
all three GL00/GL01 routes will be regenerated and remeasured after the final
source passes validation. Existing stage logs do not record the tensor-power leaf
types and exponents needed to prove that the correction is unused here.
The earlier immutable CLI SHA-256 is
`b84c2d9e89f7a28fd1c1910db8b2d43cd1fda2db460fd2e5d6fa4661eca6c36a`.

Physical measurements use the same model, numerator, kinematics, LMB and eager,
uncompiled evaluator settings in all three routes. Horner iterations are 1,
common-pair-elimination rounds are 5, and worker limits are 1. Integrated and
threshold counterterms are disabled for these isolated local-UV measurements;
the scalar correctness matrix retains its enabled subtraction settings.
Debug assertions remain enabled. Each physical generation is guarded at 500 GB;
no guard intervention occurred. Numerical passes use fixed momentum-space points
and the existing minimal-integrand benchmark setting.

Each route has three fresh-process generations. Saved states are primed before
three accepted timing passes at each of the original and 100× points. Each pass
has twenty nonempty batches and at least three actual measured seconds.
Undersized attempts are retained but excluded. Ratios use route medians; observed
spread and twice the reported within-pass SEM supply conservative repeat triggers,
not confidence intervals. Profiling runs are separate from these timings.

The native generation summary distinguishes expression construction, Spenso
tensor preprocessing, Symbolica evaluator construction and optional compilation.
Subtracting Symbolica time alone retains tensor preprocessing and other graph
work. It must not be confused with the earlier expression-construction boundary
or with the UV forest alone. The native `evaluator_symbolica_time` interval
also includes function-map preparation and conversion to numerical programs;
it is broader than the literal Symbolica `.build()` call. Separate profiling
records the actual build-call intervals and RAM at both the post-Spenso and
build-entry boundaries. These profiling times do not replace unprofiled gates.

## Physical results at the measured revision

GL00 and GL01 are two-loop four-photon amplitudes with a top-quark loop and
a top/gluon self-energy insertion. Each input has six interaction vertices and
seven internal lines (six top-quark lines and one gluon), plus four external
photons. They differ in external-photon ordering and momentum routing. The
production forest exports each contain four computation nodes: the bare root,
the full-graph DOD0 subtraction, the self-energy DOD1 subtraction and the nested
self-energy/full-graph subtraction. Thus the count includes the bare state; it
is not four distinct UV-divergent regions. Inputs are
[`GL00.dot`](../../examples/cli/aa_aa/2L/graphs/GL00.dot) and
[`GL01.dot`](../../examples/cli/aa_aa/2L/graphs/GL01.dot).

GL00 and GL01 pass all ten generation/runtime gates against orientation-erased
3D. Each graph is assessed separately; the maximum accepted median ratio is 1.15.
The [benchmark receipt](local-4d-uv-benchmarks.json) records exact observations,
card/state identities, sample counts, operation counts and uncertainty envelopes.

| Graph | Generation 4D / erased | Evaluator base / scaled | Total sample base / scaled |
| --- | ---: | ---: | ---: |
| GL00 | 1.025273 | 0.795665 / 0.805925 | 0.818698 / 0.825173 |
| GL01 | 1.001063 | 0.807388 / 0.800332 | 0.830437 / 0.823653 |

All observed range and within-pass repeat-trigger envelopes remain below 1.15.
The phase split explains the remaining generation premium:

| Graph / stage | Localized 3D (s) | Erased 3D (s) | Direct 4D (s) | 4D / erased |
| --- | ---: | ---: | ---: | ---: |
| GL00 / Expression construction | 3.052140 | 3.065155 | 1.983849 | 0.647226 |
| GL00 / Spenso preprocessing | 1.868717 | 1.880194 | 3.366697 | 1.790612 |
| GL00 / Symbolica construction | 7.015245 | 1.352501 | 1.117818 | 0.826482 |
| GL00 / Outside Symbolica construction | 4.920857 | 4.955573 | 5.364433 | 1.082505 |
| GL00 / Complete graph generation | 11.925558 | 6.316841 | 6.476487 | 1.025273 |
| GL01 / Expression construction | 3.098538 | 3.140361 | 1.956529 | 0.623027 |
| GL01 / Spenso preprocessing | 1.914180 | 1.928317 | 3.392792 | 1.759457 |
| GL01 / Symbolica construction | 7.196926 | 1.364347 | 1.100288 | 0.806457 |
| GL01 / Outside Symbolica construction | 5.012718 | 5.070101 | 5.349730 | 1.055153 |
| GL01 / Complete graph generation | 12.209644 | 6.443168 | 6.450018 | 1.001063 |

Phase medians are independent and need not sum to the pipeline median. Expression
construction subtracts both tensor preprocessing and Symbolica construction from
native graph time; it includes more than the UV forest alone. Outside Symbolica
subtracts only Symbolica and retains tensor preprocessing and other graph work.
It is a difference of interval timers, not a timestamp at builder entry. Generated
evaluator compilation is disabled, with zero compilation time in every receipt.

Expression construction is 35.28% faster for GL00 and 37.70% faster for GL01.
Including tensor preprocessing but excluding Symbolica leaves 8.25% and 5.52%
premiums respectively. These measurements do not show every preparation stage
getting faster in 4D. Native graph time also excludes CLI startup and state
serialization: localized/erased/direct child-wall medians are 13.007/7.005/8.005 s
for GL00 and 14.006/8.005/8.005 s for GL01, with one-second supervisor polling.

| Graph / route | Evaluator μs/sample, base / scaled | Total μs/sample, base / scaled | Generation VmHWM (MB) | Original instructions |
| --- | ---: | ---: | ---: | ---: |
| GL00 / 3d_local | 77.895 / 77.150 | 89.350 / 88.533 | 520.507–573.981 | 87,768 |
| GL00 / 3d_erased | 75.702 / 75.789 | 87.057 / 87.307 | 297.927–355.844 | 25,840 |
| GL00 / 4d_direct | 60.233 / 61.080 | 71.274 / 72.044 | 394.109–408.158 | 20,210 |
| GL01 / 3d_local | 74.904 / 74.738 | 86.419 / 86.609 | 526.082–577.307 | 89,639 |
| GL01 / 3d_erased | 73.159 / 73.112 | 84.426 / 84.301 | 359.023–452.305 | 24,789 |
| GL01 / 4d_direct | 59.068 / 58.514 | 70.110 / 69.434 | 377.836–419.164 | 20,029 |

Memory uses decimal MB and observed process high-water marks; sampled RSS and
the graph-reported peak are separate fields in the receipt. The localized route
also retains a summed-function-map program (26,041 instructions for GL00; 24,987
for GL01). Retained programs are not additive per-sample operation counts. All
operation counts agree across the three independently generated states.

The 18 fresh generations yielded 108 accepted passes, 2,160 batches and
4,913,040 samples over 396.547288 actual timed seconds. GL00 contributed
2,388,838 samples / 195.916594 s with pass lengths 3.362306–3.906375 s; GL01
contributed 2,524,202 / 200.630695 s with pass lengths 3.340311–7.199305 s. All
163 child receipts succeeded with nonoverlapping recorded intervals. The 109
undersized attempts remain preserved and excluded. All 36 pointwise comparisons
pass the unchanged 1e−9 tolerance; worst relative complex-norm differences are
4.3294e−15 for GL00 and 4.4759e−15 for GL01.

## Dispatch and mapping profile

Two separate fresh CLI22 runs use the same physical cards and H1/CPE5 settings,
with minimal structured profiling and the existing exclusive measurement lock.
Their generation times are excluded from the preceding unprofiled gates. The
[dispatch receipt](local-4d-uv-dispatch-profile.json) certifies complete current
accounting, successful immutable receipts and cold projection caches.

| Graph | Complete nonroot UV (ms) | Allocation, selection and cache (ms) | Fraction | Winning hard CFF (ms) | Winning template (ms) | Losing template (ms) | Row mapping (ms) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| GL00 | 1641.227 | 65.248 | 3.9756% | 20.919 | 61.344 | 29.461 | 22.963 |
| GL01 | 1654.384 | 67.809 | 4.0988% | 20.867 | 64.230 | 31.169 | 22.862 |

Both fractions pass the strict 10% limit. The denominator sums the three nonroot
nodes of each four-node forest: Taylor construction, component projection and
outer-CT assembly. Root production CFF and later tensor/evaluator work are
excluded. Dispatch includes degree/allocation, certificates, all losing native
generation and template construction, raw physical-degree reports, cache work
and outer routing selection. Outer selection subtracts only the winner’s mandatory
native generation, source reconstruction and postprocessing. Nested timers are
charged once; losing preparation is never free.

Each run has eight source requests, ten admitted candidates and ten template
builds, producing 196 selected hard residue rows (at most 62 per request). Six
unique native keys are certified by equal independent lower/upper bounds. Final
cache snapshots record 12 CFF hits / 6 misses and 1,652 subtree hits / 602 misses,
with no LRU evictions and component-local row results cleared. Request/build
counts do not establish unique source bindings or full sample tuples.

## Correctness and merge checks

The frozen numerical source passes 630 focused unit tests. The final dedicated
integration selection contains 183 tests: 167 scalar checks, two physical
GL00/GL01 three-route comparisons and fourteen raised-propagator, UV-composition,
cut/threshold and analytic checks. All 183 passed on the frozen numerical source.

The scalar group has completed: all 167 tests passed. Its 166 single-run
setup/generation observations have a median direct/erased ratio of 0.814800;
131 are below 1.00 and eight exceed 1.15. The largest are GL21 (1.417311) and
GL17 (1.368180). These observations include evaluator construction and enabled
integrated/threshold subtraction; they are not measurements of expression
construction alone or repeated physical acceptance benchmarks. A separate frozen
GL21 base diagnostic passed and locates its representative cost. The no-numerator
UV/forest orchestration improves from 972.364 ms to 405.324 ms; the graph-owned
product numerator instead takes 5,944.078 ms in 4D versus 3,377.569 ms in erased
3D. New hard projection accounts for only 55.666 ms and outer routing selection
for 82.152 ms. The largest remaining gaps follow integrated-addback localization
and precede the next component or final assembly. The shared final-integrand
simplification and forest-aggregation boundary is therefore the next generic
optimization target; these logs do not identify one dominant internal function.
Numerator tensor preprocessing is also higher (1,005.779 versus 644.301 ms),
while its Symbolica stack construction is slightly lower (172.676 versus
176.267 ms). This is a single diagnostic with integrated/threshold terms enabled,
not a new physical gate. Exact intervals and limitations accompany the
[dispatch profile](local-4d-uv-dispatch-profile.json).

The scalar selection contains 166 route-comparison cases across 49 graph labels,
including quadratic/quartic numerator variants, plus one sampling-scale check.
It is not 167 distinct topologies. The physical comparisons use the Compare
orchestrator to exercise both forest implementations. Existing physics
expectations and tolerances are preserved.

The PR targets `codex/raised-energy-cff-reviewed`, its original base at
`78395e3ab3ddd8d8f62b2f674d7488484eace197`. The existing GitHub CI and Nix workflows
also run for this target. Formatting, workspace check/clippy, the normal CI
suite and the explicit slow-case selection must be green before readiness.

The initial full CI run exposed three failures among 2,137 tests:
`finite_part_ghost_2loop`, `se1l_uv` and `epem_a_bbx_amp_uv`. The first
constructs analytically integrated 4D renormalization terms and does not enter
the canonical 3D projection, CFF dispatch or evaluator stages. A matched
original-base comparison identified the first differing expression before
Vakint: analytic spin expansion changed `(p.k)^2` into `p^2 k^2` by reusing a
contraction index. The analytic owner now protects certified scalar products
through existing aliases while expanding open tensor contractions. A related
shared executor defect in odd tensor powers is also corrected.

The rebuilt focused selection passes all twenty analytic/physical checks,
including all three CI failures and the new powered-product/compact-metric
regression, with unchanged expectations and tolerances. The preceding broader
run passes all Spenso/Idenso checks, including the new five-leaf odd-power
coverage. The final source still requires fresh physical timings, the dedicated
183-test selection and green CI before readiness.

## Deferred Symbolica investigation

In the matched CLI18 diagnostic, GL262 forest construction took 900.274 s in
direct 4D and 2,254.471 s in erased 3D, a ratio of 0.3993. Direct tensor
preprocessing completed, but the erased run hit its time cap during
preprocessing. There is no completed paired measurement of the full pipeline.
The later direct run failed inside Symbolica evaluator construction with
compilation disabled; no complete GL262 evaluator/runtime result is claimed.

The separate nonsymmetric `f(x,y)` import regression is reproduced in the
[standalone Symbolica package](../../tests/artifacts/aa_aa_uv_slowdown/symbolica_evaluator_mre/README.md).
Its 2,480-byte input imports as `f(y,x)` when the reader registers the argument
symbols in reverse order. That confirmed import defect is distinct from the
large evaluator panic, which has no confirmed standalone reproducer.
