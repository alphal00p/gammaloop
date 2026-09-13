# Sampling bridge compilation during warmup

Status: ownership design audited on 2026-09-13. The initial cached f64 bridge
milestone passed 154 library sampling tests, both saved-state/summed API tests,
core/API test checking and clippy. The native host extension now binds configured
precisions and lazy explicit requests from one catalogue/program cache; its
eight focused production gates, 169-test broader core suite, both saved-state/
summed API regressions and six event/precise integration API tests pass. No throughput
or physical variance improvement is claimed before measurement.
It extends [the implementation plan](../../../ADVANCED_SAMPLING_PLAN.md) and
the [amplitude benchmark study](AMPLITUDE_BENCHMARK_CANDIDATES.md).

## Scope and existing owners

Compile each currently supported graph's sampling bridge once per successful
warmup, then reuse it across draws, orientations, rotations, and precision
evaluations. The setup owns one resolved catalogue of canonical
`SamplingChannelId` values; native bridges retain bindings indexed by those
same entries and never independently select or enumerate channels.

The initial scope is existing LMB maps, full-rank amplitude surface maps, and
standalone full-parent `phase_space(cut(...))` maps. Their equations depend on
fixed graph routing, runtime settings, warmed masses, and fixed externals.
Their radial roots still depend on the sampled direction and are solved per
point. This does not implement conditional side-threshold compositions.

Relevant existing owners and call sites:

| Responsibility | Existing source |
| --- | --- |
| Nonpersistent, cloneable cache | `utils/mod.rs`: `RuntimeCache<T>` |
| Process and graph warmup | `ProcessIntegrand::warm_up`, `ProcessIntegrandImpl::warm_up`, `GraphTerm::warm_up` |
| Physical map construction | Fresh `GraphTerm::compile_sampling_bridge`, native `bind_sampling_bridge`; amplitude/cross-section implementations |
| Discrete forward map | `process/gammaloop_sample.rs`: `parameterize` |
| Summed and direct-momentum map use | `process/mod.rs`: common evaluation paths |
| Stateful Symbolica proxy | `sampling_partition.rs`: `SamplingScoreFunction` |
| Eager evaluator buffers | `sampling_evaluator.rs`: `SamplingExpressionEvaluator` |
| External numeric cache | `settings/runtime/kinematic.rs`: `Externals` |
| Worker cloning | `integrate/mod.rs`: worker setup and `evaluate_sample_list` |

Paths beginning with `process/` are relative to
`crates/gammalooprs/src/integrands/`; other paths in this table are relative to
`crates/gammalooprs/src/`.

## Construction and retrieval are different operations

Keep `compile_sampling_bridge(parameterization, e_cm, externals, orientation)`
as an explicit fresh constructor. Inspection and tests deliberately pass
different external vectors without mutating the process. It must use every
supplied argument and must neither return nor overwrite an unrelated cache.

Add cached retrieval at the existing graph owner, with no pretend override
arguments. It borrows the warmed bridge or returns an actionable missing-warmup
error. Do not silently rebuild on each draw or silently ignore changed settings.
Direct construction remains useful without establishing a warmed runtime epoch.

Replace the three production construction sites with this retrieval: discrete
forward sampling, summed-channel evaluation, and direct-momentum inverse
weighting. Where these paths currently resolve channel IDs again, retrieve IDs
from the same bridge. UI/grid inspection may retain explicit metadata resolution.
Only modes that actually use bridges need a runtime bridge; do not accidentally
activate unsupported map definitions in unrelated sampling modes.

The graph's existing `LmbMultiChannelingSetup` owns the catalogue, compiled
proxy programs and typed Double/Quad/Arb bridges in `RuntimeCache` fields,
instead of an `Arc` around the entire bridge. The hot path borrows it; the summed loop can finish an owned
forward-map result before mutably evaluating the graph. No per-draw bridge clone
is necessary. Existing immutable geometry callbacks may continue sharing `Arc`s.

## Transactional warmup and invalidation

1. Invalidate every graph's old sampling cache at process warmup entry, before
   runtime validation can fail. A failed warmup must not expose a previous epoch.
2. Run the existing model, graph, evaluator, and kinematic initialization.
   Amplitude warmup currently improves and caches externals **after** graph-term
   warmup. Compiling at the end of `GraphTerm::warm_up` would capture the wrong
   external data; compile after successful process-level initialization.
3. Construct all needed bridges into temporary storage, using the final improved
   externals and current real masses. Validate the complete set before publishing.
4. Publish the bridges only after all construction succeeds. No graph should
   retain a new bridge while another retains an old one after an error.

The existing mutable-settings invalidation boundary was expanded and renamed
from `invalidate_event_processing_runtime` to `invalidate_runtime_caches`,
retaining one settings owner. `ProcessIntegrand::get_mut_settings` invalidates
sampling caches before
returning its mutable reference. Graph filtering/replacement, routing changes,
and model changes require the same warmup lifecycle. Public direct Rust mutation
must be followed by warmup; untracked field mutation cannot safely coexist with
a promise of automatic cache coherence. Prefer existing setters over hashing the
entire runtime settings structure on every sample.

CLI `set` uses `get_mut_settings`; evaluate, integrate, approach, and profile
commands normally warm up before evaluation. Python sample evaluation delegates
to the same command owner. Global model changes are applied through the effective
model resolved before that warmup, so bindings must be rebuilt even if settings
are unchanged. Event-return overrides invalidate settings before warming and
again when restored: subsequent evaluation must follow the same warmup contract.

Audit direct Rust call sites during implementation. In particular,
`uv/profile.rs` clones a warmed integrand and changes `generate_events` through
`get_mut_settings` when projecting compatible cuts. Do not break this workflow:
either its raw-momentum route demonstrably needs no sampling bridge, or restore
warmup at the appropriate existing owner. A benign event setting is not evidence
that arbitrary mutable settings references can preserve the sampling cache.

### Existing stale-external-cache gap

`Externals::get_dependent_externals` returns its f64/f128 cache before inspecting
`Constant.momenta`; `improve_and_cache` calls that accessor before replacing the
caches. Direct in-place Rust edits of momenta can therefore survive even a later
warmup as stale numeric inputs. Reconstructing settings through TOML drops these
serde-skipped caches and does not have the same failure mode.

Clear the existing external caches at the appropriate mutation/recomputation
boundary before rebuilding improved externals. Merely invalidating the bridge
does not fix this upstream dependency. Add a regression changing a spatial
component while retaining `e_cm`, then verify fresh improved data and map roots.

## Program sharing and worker-local buffers

Existing integration workers clone the warmed integrand. `RuntimeCache` clones
its payload and serializes no payload, so saved states remain independent of
compiled programs and reload requires warmup. One bridge per graph is sufficient
for the presently supported orientation-independent maps; do not construct
hundreds of identical bridges for a graph's orientation catalogue.

The original `SamplingScoreFunction::Clone` shared its opaque callback `Arc`,
including an `Arc<Mutex<SamplingExpressionEvaluator>>` for Symbolica scores.
Worker clones therefore shared mutable evaluator buffers. The score owner now
retains a concrete compiled-expression variant whose clone copies the existing
evaluator into a fresh local mutex. Immutable Rust callbacks remain shared.
The partition borrows map/proxy evaluators at each point; it must not clone
their new worker-local buffers per draw. No reparse, recompilation or second
evaluator engine is introduced by worker cloning.

## Future conditional maps

Cache routing, channel definitions, expression/dual programs, derivative shape,
and immutable parameter bindings. Keep point-dependent `t*`, cut externals,
sampled complements, radial roots, overlap centers, branch status, and inverse
partition scores outside this cache. Use the existing prepared-cut/runtime-context
owners to bind each channel's current point, including foreign-channel inversion.

An absent/pinched/existing transition retains the canonical channel and selects
its valid normalized branch; it does not trigger enumeration or recompilation.
Specialize programs by equation and coordinate structure, not by orientation
number when orientation only changes parameters or masks. Rotation and precision
conversion of mapped samples do not justify copying an unchanged sampling program.

## Acceptance and performance evidence

- Compare fresh and warmed bridges on the generated six-dimensional kite and a
  genuine standalone cut map; retain Gaussian normalization and physical
  summed-versus-Monte-Carlo equivalence gates.
- Prove repeated draws and orientation/rotation evaluations reuse compilation,
  using program identity or focused instrumentation, not elapsed time alone.
- Fresh construction with external data B must differ where expected from warmed
  data A while leaving A's cache unchanged.
- Settings mutation must reject cached evaluation until warmup. Cover selection,
  radial profile, frame, externals, model masses, and filtered graph views.
- Check failed warmup clears old/partial bridges, in-place external mutation
clears every numeric cache, and saved/reloaded states warm to identical results.
- Check worker clones reuse programs without shared mutable proxy buffers;
  threaded proxy outputs must match serial outputs without compiling per point.
- Measure warmup cost, repeated-draw cost, memory, and 1/20-worker throughput with
  and without eager proxies before claiming improvement. These cache changes
  should not change densities, Jacobians, physical values, or estimator variance.
