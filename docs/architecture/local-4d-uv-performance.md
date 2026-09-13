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

## Measurement boundary

The numerical implementation measured here is
`f938f9c1386eb2670a0a7f603f5eb0878086352a`. Subsequent merge cleanup changes tuple
type aliases, comments and test iteration syntax, without changing the production
algorithms. The immutable CLI SHA-256 is
`b84c2d9e89f7a28fd1c1910db8b2d43cd1fda2db460fd2e5d6fa4661eca6c36a`.

Physical measurements use the same model, numerator, kinematics, LMB and eager,
uncompiled evaluator settings in all three routes. Horner iterations are 1,
common-pair-elimination rounds are 5, and worker limits are 1. Integrated and
threshold counterterms are disabled for these isolated local-UV measurements;
the scalar correctness matrix retains its enabled subtraction settings.
Debug assertions remain enabled.

Each route has three fresh-process generations. Saved states are primed before
three accepted timing passes at each of the original and 100× points. Each pass
has twenty nonempty batches and at least three actual measured seconds.
Undersized attempts are retained but excluded. Ratios use route medians; observed
spread and twice the reported within-pass SEM supply conservative repeat triggers,
not confidence intervals. Profiling runs are separate from these timings.

The native generation summary distinguishes expression construction, Spenso
tensor preprocessing, Symbolica evaluator construction and optional compilation.
"Before Symbolica" includes tensor preprocessing. It must not be confused with
the earlier expression-construction boundary or with the UV forest alone.

## Final physical results

GL00 is complete. GL01's final repeat and the separate dispatch audit are pending.
This section will be finalized before the PR is marked ready.

| GL00 stage | Localized 3D (s) | Erased 3D (s) | Direct 4D (s) | 4D / erased |
| --- | ---: | ---: | ---: | ---: |
| Expression construction | 3.052140 | 3.065155 | 1.983849 | 0.647226 |
| Spenso preprocessing | 1.868717 | 1.880194 | 3.366697 | 1.790612 |
| Everything before Symbolica | 4.920857 | 4.955573 | 5.364433 | 1.082505 |
| Complete graph generation | 11.925558 | 6.316841 | 6.476487 | 1.025273 |

Stage medians are calculated independently and need not sum to the median of
the complete pipeline. The complete pipeline includes Symbolica construction,
with generated evaluator compilation disabled.

GL00's evaluator ratios are 0.795665 at the base point and 0.805925 at the scaled
point; total-sample ratios are 0.818698 and 0.825173. Every gate is below 1.15,
including the repeat-trigger envelopes. Expression construction is faster, while
the larger tensor-preprocessing cost leaves an 8.3% premium before Symbolica.
This is not a claim that every individual preparation stage is faster in 4D.

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
construction alone or repeated physical acceptance benchmarks. A separate
stage profile is pending to locate the representative GL21 cost.

The scalar selection contains 166 route-comparison cases across 49 graph labels,
including quadratic/quartic numerator variants, plus one sampling-scale check.
It is not 167 distinct topologies. The physical comparisons use the Compare
orchestrator to exercise both forest implementations. Existing physics
expectations and tolerances are preserved.

The PR targets `codex/raised-energy-cff-reviewed`, its original base at
`78395e3ab3ddd8d8f62b2f674d7488484eace197`. The existing GitHub CI and Nix workflows
also run for this target. Formatting, workspace check/clippy, the normal CI
suite and the explicit slow-case selection must be green before readiness.

The initial full CI run exposed a repeatable failure in
`finite_part_ghost_2loop`, on its closed-quark-loop contribution. This test
constructs analytically integrated 4D renormalization terms and does not enter
the canonical 3D projection, CFF dispatch or evaluator stages. Its unchanged
physics expectation remains the acceptance criterion; diagnosis and correction
are pending. Ubuntu's full suite has 2,134 passes and three failures: that test
and the integrated UV checks `se1l_uv` and `epem_a_bbx_amp_uv`. Their UV mass
invariance failures are also under investigation with unchanged expectations.

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
