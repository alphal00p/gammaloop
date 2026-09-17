# GL297: exact recovery and rescue of both interrupted integrations

The two GL297 jobs stopped because Clarabel returned an uncertified overlap
result. Both original failures have now been recovered from their unchanged
integration checkpoints and reproduced individually with the preserved old
binary. They are numerical SOCP failures, not evidence that the physical
integrand diverges at those points.

| Replay | Optimized LMB | Advanced sampling |
|---|---:|---:|
| Last original completed iteration | 715 | 649 |
| Original samples | 58,572,800 | 53,166,080 |
| Replayed iteration | 716 | 650 |
| Samples after replay | 58,654,720 | 53,248,000 |
| Failed samples recovered | 1 | 1 |
| Graph / sampling channel | 0 / 8 | 0 / 3 |
| Failing overlap surfaces | 1, 3 | 0, 1, 2, 3, 4, 5 |
| Original Clarabel status | AlmostPrimalInfeasible | AlmostPrimalInfeasible |
| Original solve iterations | 18 | 19 |
| Original primal residual | 0.1274421456440971 | 0.0010779539713490188 |
| Equivalent problem with all dimensional inputs divided by 2 | PrimalInfeasible | PrimalInfeasible |
| Rescaled solve iterations | 16 | 21 |

The optimized failure has cube coordinates

```text
0.33168777529733795 0.9182901767427527 0.8871151234013295
0.5237759320929373 0.4056451411943530 0.6190324006441688
0.6178135113081745 0.7667261103805785 0.9016624966708675
0.7373798919654603 0.2773975750381326 0.36987952654913625
```

The advanced failure has cube coordinates

```text
0.5156951665734136 0.5364636794716388 0.8318524038686087
0.5710377735171434 0.7546447565761167 0.8524813371503793
0.3525930855038039 0.014530525876098557 0.5224313652832420
0.6316763485194195 0.5120208610432367 0.5806413236969333
```

## Acceptance evidence and limits

With the first containment patch, both copied integrations completed their next
81,920-point iteration and wrote the complete failed sample to
`numerical_stability/failed_samples.jsonl`. The ordinary evaluation recorded
Double, Quad and Arb attempts, followed by `is_nan=true` and the original
SOCP diagnostic. Forced Arb inspection also retained the same numerical failure.
The conic preparation uses its canonical double-precision geometry, so simply
forcing the later integrand evaluation to Arb does not repair this failure.

The output API represents an invalid contribution by zero together with
`is_nan=true`; this is not a certified physical zero. Default minimal inspection
omits that metadata. The archived detailed inspections explicitly enable event
reporting, retain the invalid flag, and disable additional component weights.
The generic containment/failed-sample path remains necessary for exhausted
rescue attempts.

These two particular failures have a simple rescue. The archived raw conic
problems reproduce the original status, residuals and iteration counts exactly.
Both give a full `PrimalInfeasible` result after dividing **every** entry of the
right-hand side by 2. This includes all masses and spatial/energy shifts. With
zero quadratic objective and homogeneous cones, keeping A, q and the cones
unchanged is the same feasibility problem in different energy units. Any
returned primal coordinates must be multiplied by 2 before the existing strict
physical-interior test. A retry is refused if the rescaling does not round-trip
exactly in binary64, including problematic subnormal inputs.

The production change tries this once only when the original solve has neither
a physically certified interior center nor a full infeasibility certificate.
It preserves the same physical-center validation and never equates an
`AlmostPrimalInfeasible` status with an empty overlap. If neither attempt is
certified, both diagnostics accompany the typed numerical error and subsequent
NaN handling. Normalizing by 600 instead was unsuccessful in this experiment;
there is no claim that an arbitrary change of units always improves the solver.

The exact two matrices are in
`tests/resources/graphs/ir_safe_thresholds/GL297_socp_failures.json`. The regression
`gl297_recorded_uncertified_socp_problems_are_rescued_without_dropping_overlaps`
uses the same retry implementation as production and independently checks the
physical energy inequalities in the original units. A future solver may certify
the first attempt directly; the test does not demand the historical approximate
status.

**Final validation:** the production exact-matrix regression passes. Both full
physical points now evaluate successfully with `is_nan=false`, using ordinary
Double precision, and independently with forced Arb precision. Native traces
show the failed SOCP being retried at scale 2 and returning `PrimalInfeasible`;
the remaining overlaps and complete physical integrand are then evaluated.
These are complete sampling-weighted inspection values at unit outer Havana
weight, in picobarns, for the two distinct recovered points:

| Point | Ordinary Re | Ordinary Im | Relative Double/Arb difference, Re / Im |
|---|---:|---:|---:|
| Optimized | +1.0553291657152030e-05 | +9.7550512154723612e-06 | 1.493e-14 / 5.453e-14 |
| Advanced | -1.2246349837150626e-04 | +1.6310832812435497e-04 | 1.328e-15 / 3.656e-15 |

All four calls report two stable rotation probes. The final inspections,
commands and binary hash are in [rescued_points.json](gl297_socp_failures_20260917/rescued_points.json),
with compressed native traces alongside it. The earlier one-iteration replay
used the containment patch, before this rescue was added; the final validation
re-evaluates the individual recovered points and does not overwrite the original
long-run checkpoints.

## Reproduction and provenance

Run the exact-matrix regression from the repository development environment:

```bash
cargo test --profile dev-optim -p gammalooprs \
  gl297_recorded_uncertified_socp_problems_are_rescued_without_dropping_overlaps \
  -- --nocapture
```

For the optimized point, load its existing generated state without saving, apply
its saved integration settings, and inspect graph 0 / sampling channel 8:

```bash
target/dev-optim/gammaloop --read-only-state -n \
  -s examples/cli/epem_a_ttxh/NNLO/gammaloop_state_GL297_optimized_lmbs \
  run -c 'set process -p epem_a_tth -i NNLO file examples/cli/epem_a_ttxh/NNLO/workspaces/GL297_optimized_lmbs_workspace/integrands/epem_a_tth@NNLO/settings.toml; set process -p epem_a_tth -i NNLO kv general.generate_events=true general.store_additional_weights_in_event=false; inspect -p epem_a_tth -i NNLO -x 0.33168777529733795 0.9182901767427527 0.8871151234013295 0.5237759320929373 0.405645141194353 0.6190324006441688 0.6178135113081745 0.7667261103805785 0.9016624966708675 0.7373798919654603 0.2773975750381326 0.36987952654913625 -d 0 8; quit'
```

For advanced sampling, replace the state/workspace suffix with
`advanced_sampling`, substitute the second cube above, and use `-d 0 3`.
Add `-f` to force Arb. Native file tracing requires a disposable state copy,
since `--read-only-state` intentionally disables file logging.

The [machine artifact](gl297_socp_failures_20260917/replay_results.json) contains
both complete Havana samples and outer weights, detailed inspections, exact
commands, binary hashes, original manifests and state-copy hashes. Compressed
native traces and terminal logs are alongside it. The standalone
[same-library solver experiment](gl297_socp_failures_20260917/socp_rescale.rs)
and [its output](gl297_socp_failures_20260917/socp_rescale.log) isolate the conic
problem without regenerating either integrand.

All integration replays used disposable state/workspace copies under
`target/gl297_failure_diagnostics`, sequentially on 50 workers. Every one of the
36 non-log state files was SHA256-verified identical to its original. Original
workspaces were neither restarted nor overwritten. The original checkpoint
hashes are retained in the artifact.

A separate pre-existing resume-fingerprint nondeterminism initially blocked
both old and patched binaries: two identical one-point runs with the old binary
produced fingerprints `ad87249704a27250` and `883a7d44f6f81574`. For this diagnostic
only, a one-point run in a separate workspace established the current-process
fingerprint, which was copied into the **disposable replay manifest** in the
same CLI session. The original manifest and checkpoint bytes were preserved;
all model, generation, graph and checkpoint data remained unchanged. The exact
original SOCP residuals/iterations reproduced after this workaround. This was
an explicitly recorded diagnostic bypass, not a production resume repair.
Only the copied runtime sample cap was changed to stop after one iteration.
