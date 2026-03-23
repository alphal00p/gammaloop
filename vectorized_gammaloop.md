# Vectorized GammaLoop Plan

Status: implementation complete for the planned f64 batching scope.
Last updated: 2026-03-21

## Implementation status

Implemented:
- `runtime.integrator.gammaloop_batch_size` now controls entry-point rechunking for both `evaluate_samples_raw(...)` and `evaluate_momentum_configurations_raw(...)`.
- The f64 path keeps an exact scalar fast path for effective batch size `<= 1`.
- f64 batched evaluator calls now reach Symbolica through batch-evaluator entry points.
- Request-weighted `average_evaluator_batch_size` is surfaced through `EvaluationMetaData` and aggregated into integration statistics.
- Per-sample timing in the batched f64 path is normalized from explicit timing ledgers so the scalar interpretation of the timing fields is preserved.
- Amplitude and cross-section graph-term stacks both use the shared batched f64 substrate.
- `f128` and `ArbPrec` remain explicitly scalar.

Current intentional fallbacks:
- If runtime cache is enabled, the path stays scalar.
- If runtime `general.evaluator_method` is `Iterative`, the path stays scalar.
- Cross-section threshold subtraction remains in the same state as before: the existing runtime-enabled `todo!()` path is still not implemented, so batching does not change that behavior.

Validation:
- `cargo check`
- `EXTRA_MACOS_LIBS_FOR_GNU_GCC=T cargo test -p gammaloop-integration-tests --test test_runs test_batched_evaluation -- --nocapture`

Post-review fixes:
- Batched summed evaluation now sets `override_if` for `OwnedOrientations::All`, matching the scalar summed paths.
- Primary-rotation-only timing and batch-size accounting is restored; secondary stability rotations now contribute only to `total_timing`.
- Batched `total_timing` now includes stability checks, final result assembly, rotation preparation, and batch orchestration overhead that the scalar path already counted.
- Batch-output decoding no longer allocates one `Vec` per result chunk on the compiled f64 path.
- `evaluate_batch_single_f64(...)` now guards its single-output assumption and allocates compiled output buffers with the correct size.

Regression coverage added:
- `tests/tests/test_runs.rs:test_batched_evaluation`
  - standalone `gg_hhh` amplitude setup with physical kinematics and threshold counterterms enabled
  - standalone `e+ e- > d d~ g` differential-LU-style setup with generated events, selectors, observables, and additional weights
  - batched vs scalar comparisons for both raw API batch entry points, ignoring timing-only fields
  - explicit assertions that the evaluator batch monitor stays at `1.0` for scalar runs and exceeds `1.0` for batched runs

## Locked-in design choices

Confirmed from the latest review:
- Only the `f64` path is batched/vectorized.
- `f128` and `ArbPrec` stay explicitly scalar and remain the reference paths for side-effect and correctness checks.
- Amplitudes are in scope as well, not just cross sections.
- For now the only backend requirement is to call batch evaluators at the Symbolica boundary; no SIMD-specific backend work is required yet.
- The public `evaluate_samples_raw` call may receive a very large slice, but GammaLoop should internally re-batch it into smaller f64 sub-batches.
- The internal batch size becomes a runtime integrator setting exposed through `IntegratorSettings` as `gammaloop_batch_size`.
- Batch-of-one should stay very close to the current runtime, so allocation churn must be avoided and an exact scalar fast path is required.
- The batching layer should stay as explicit and local as possible, so the mostly linear evaluation stack remains readable and early exits remain effective.

## Updated assessment of the batching logic

The logic is sound.

What is correct:
- Re-batching a very large incoming slice into smaller internal sub-batches is the right control knob.
- It bounds the active frontier and keeps the batching state local to one `evaluate_samples_raw(...)` call.
- It still lets the hot evaluator sites accumulate enough requests to make backend vectorization worthwhile later.
- It preserves all current early-termination opportunities above the evaluator boundary.

Important caveats:
- `gammaloop_batch_size` bounds the number of active sample states, not the exact size of every concrete evaluator queue.
- Some queues can exceed it because of orientation expansion and threshold/cartesian branching.
- Sparse evaluator groups will still flush small; the occupancy wins will come mostly from the hot evaluator sites.

So the right interpretation is:
- `gammaloop_batch_size` caps the active frontier.
- Each concrete evaluator queue still needs its own flush threshold and stage-end flush points.

## Final deep-pass observations from the code

### 1. `gammaloop_batch_size` belongs in runtime integrator settings, not evaluator settings

Relevant code:
- `crates/gammalooprs/src/settings/mod.rs:47`
- `crates/gammalooprs/src/settings/runtime.rs:119`
- `crates/gammalooprs/src/processes/mod.rs:23`
- `crates/gammalooprs/src/settings/global.rs:29`

Observation:
- `IntegratorSettings` already owns the runtime controls for how incoming work is fed through the integrator.
- `processes::EvaluatorSettings` remains generation-scoped and contains compilation/optimization knobs.
- `gammaloop_batch_size` only controls the initial rechunking at the batch entry point; it does not describe low-level evaluator behavior and it does not need to cap downstream queue growth.

Planned resolution:
- Add `gammaloop_batch_size` to runtime `IntegratorSettings`.
- Keep it strictly as the entry-point rechunking limit used by `evaluate_samples_raw(...)` and the other batch raw entry points that should share the same policy.
- Do not reinterpret it later as a cap on concrete evaluator queues.
- Serialize it with the existing `SHOWDEFAULTS` helper path.

Why this matters:
- This keeps the setting aligned with what it actually controls and avoids coupling entry-point batching policy to evaluator-generation settings.

### 2. The batching layer should stay local and explicit

Relevant code:
- `crates/gammalooprs/src/integrands/process/mod.rs:592`
- `crates/gammalooprs/src/integrands/process/mod.rs:1704`
- `crates/gammalooprs/src/integrands/process/mod.rs:2536`

Observation:
- The stack is still mostly linear.
- The batch engine should therefore not become a generic opaque scheduler.

Conclusion:
- Keep stage-local batch collectors at evaluator boundaries.
- Let sample states continue through the same stages as today.
- Only the evaluator boundary becomes batched.

### 3. The f64 `ParamBuilder` already contains the right low-allocation template

Relevant code:
- `crates/gammalooprs/src/integrands/process/param_builder.rs:715`
- `crates/gammalooprs/src/integrands/process/param_builder.rs:790`
- `crates/gammalooprs/src/integrands/process/evaluators.rs:995`

Observation:
- The current f64 path mutates `self.values[value_index]` in place.
- `InputParams` already supports owned storage, but the f64 path currently uses borrowed storage.

Conclusion:
- The batched f64 path should use reusable queue arenas and parameter spans, not one `Vec` per request.
- The arena fill path should copy the prebuilt template slice into the next span and patch only the dynamic fields.
- Batch-of-one should bypass the arena entirely and use the current scalar path unchanged.

### 4. `EvaluationMetaData` changes are cross-cutting, not local

Relevant code:
- `crates/gammalooprs/src/integrands/evaluation.rs:594`
- `crates/gammalooprs/src/integrands/mod.rs:408`
- `crates/gammalooprs/src/integrands/mod.rs:568`
- `crates/gammalooprs/src/integrands/builtin/h_function.rs:194`
- `crates/gammalooprs/src/integrands/process/mod.rs:2694`
- `crates/gammalooprs/src/settings/runtime.rs:420`

Observation:
- `EvaluationMetaData` is manually constructed in several places.
- Integration summaries are built through `StatisticsCounter` and `IntegrationStatisticsSnapshot`.

Conclusion:
- Adding an evaluator-batch monitor is not just one new metadata field.
- It must also be threaded through:
  - `EvaluationMetaData::new_empty()`
  - manual metadata constructors
  - `StatisticsCounter`
  - `IntegrationStatisticsSnapshot`
  - any status-table rendering that should expose it.

Planned metric:
- Track the average batch size seen at the concrete evaluator call sites.
- Make it request-weighted, not flush-weighted.
- For one sample, compute:
  - `sum(batch_size_seen_by_each_logical_evaluator_request) / number_of_logical_evaluator_requests`
- For integration summaries, aggregate the same numerator/denominator across samples.

Why request-weighted is necessary:
- A plain average over flushes would misreport samples that contribute very different numbers of evaluator requests.

### 5. Timing normalization must stop using plain per-sample wall clock in the batched f64 path

Relevant code:
- `crates/gammalooprs/src/integrands/process/mod.rs:1815`
- `crates/gammalooprs/src/integrands/process/mod.rs:2198`
- `crates/gammalooprs/src/integrands/process/mod.rs:2550`
- `crates/gammalooprs/src/integrands/process/evaluators.rs:109`

Observation:
- Current timings are accumulated directly around scalar evaluator calls and around the scalar per-sample integrand evaluation.
- In the batched path, `start.elapsed()` for a whole sample would include waiting for peer samples in the same sub-batch.
- That would destroy the scalar interpretation of the timing fields.

Planned resolution:
- Introduce explicit per-sample timing ledgers for the batched f64 path.
- Sample-local scalar work stays charged directly to the owning sample.
- Shared batch flush time is apportioned back to participating samples by logical evaluator-request count, not duplicated onto every sample and not split evenly per sample.
- `integrand_evaluation_time`, `evaluator_evaluation_time`, and `total_timing` in the batched f64 path should be synthesized from that ledger instead of raw per-sample wall clock.

Normalization rule:
- If a batch flush covers `N` logical evaluator requests and takes wall time `T`, each request carries `T / N` of evaluator time.
- A sample gets the sum over its own logical requests.
- This preserves the current interpretation: if batching makes evaluator throughput 4x better, the reported average per-sample evaluator time also drops by about 4x.

### 6. `SingleOrAllOrientations` should not be stored across batching boundaries

Relevant code:
- `crates/gammalooprs/src/integrands/process/evaluators.rs:60`
- `crates/gammalooprs/src/integrands/process/evaluators.rs:502`

Observation:
- The orientation object is a borrowed view.

Conclusion:
- Store only owned orientation state in sample/cut state.
- Recreate the borrowed `SingleOrAllOrientations` view only when filling a concrete evaluator queue.

### 7. The batch monitor should instrument the low-level evaluator boundary, not the graph-level stage count

Relevant code:
- `crates/gammalooprs/src/integrands/process/evaluators.rs:502`
- `crates/gammalooprs/src/integrands/process/evaluators.rs:1106`

Observation:
- Different evaluator methods collapse work differently before reaching `GenericEvaluator`.
- In particular, `Iterative` and `SingleParametric` do not map orientations to low-level calls in the same way.

Conclusion:
- The batch-size monitor should be updated where a concrete `GenericEvaluator` batch flush happens.
- Measuring only batch size in terms of top-level samples would be misleading.

### 8. Cross section still has a staged continuation point, and amplitude remains the simpler second consumer

Relevant code:
- `crates/gammalooprs/src/integrands/process/cross_section_integrand.rs:667`
- `crates/gammalooprs/src/integrands/process/cross_section_integrand.rs:913`
- `crates/gammalooprs/src/integrands/process/amplitude/mod.rs:200`
- `crates/gammalooprs/src/subtraction/amplitude_counterterm.rs:132`

Observation:
- Cross section still needs staged sample state because pass 2 depends on pass 1 outputs and selector/event early exits sit in between.
- The currently wired cross-section threshold-subtraction call site is still `todo!()` when enabled.
- Amplitude is simpler, but it still has threshold-counterterm evaluator boundaries of the same general kind.

Conclusion:
- Cross section remains the first demanding consumer of the batch substrate.
- Amplitude should then reuse the same substrate directly.
- The queue API must be shaped so the currently inactive LU threshold-subtraction path can plug in later without another control-flow rewrite.

## Revised architecture

### A. Add `gammaloop_batch_size` to runtime integrator settings

Target structure:
- runtime `IntegratorSettings` owns `gammaloop_batch_size`.
- generation-side evaluator settings remain generation-only.

Serde requirement:
- `gammaloop_batch_size` must serialize with the standard `SHOWDEFAULTS` behavior through the existing helper.

### B. Keep an exact scalar fast path

Rule:
- If the effective internal f64 sub-batch size is `<= 1`, dispatch to the current scalar path unchanged.

Reason:
- This preserves current performance and keeps a direct reference path.

### C. Split large incoming sample slices into internal f64 sub-batches

At `evaluate_samples_raw`:
- accept the full incoming slice as today,
- internally process it in `gammaloop_batch_size` chunks for the f64 path,
- preserve output order exactly,
- finalize/process results in the same order as today.

### D. Use explicit stage-local batch collectors, not a global scheduler

Preferred shape:
- `evaluate_samples_raw` owns the active sample states for one internal sub-batch.
- Lower stages keep their current linear structure.
- Evaluator boundaries own local queues and flush them at natural stage boundaries or when full.

### E. Preallocate queue arenas and write params into spans

Per concrete evaluator queue:
- keep request metadata separate from raw parameter storage,
- reuse one contiguous parameter arena,
- reuse one contiguous output arena,
- avoid per-request heap allocation.

Required low-level additions:
- `GenericEvaluator` needs f64 batch helpers for both single-output and multi-output cases.
- The batch helpers should sit behind the same scalar fallback API.

### F. Make timing and monitoring first-class in the batch substrate

Each batch flush should report:
- wall time for the flush,
- number of logical evaluator requests in the flush,
- ownership of those requests by sample state.

Each sample ledger should accumulate:
- direct scalar stage time,
- apportioned evaluator flush time,
- request-weighted batch-size numerator and denominator.

Then `EvaluationMetaData` for the batched f64 path is assembled from the ledger instead of from one wall-clock span.

### G. Keep higher precisions entirely scalar

Policy:
- only f64 uses the batched substrate,
- escalated `f128` and `ArbPrec` jobs drop back to the existing scalar code path.

## Revised implementation plan

### Phase 1: settings and metric groundwork

Deliverables:
- add `gammaloop_batch_size` to runtime `IntegratorSettings`,
- decorate serialization with the existing `SHOWDEFAULTS` helper,
- extend `EvaluationMetaData`, `StatisticsCounter`, and `IntegrationStatisticsSnapshot` with evaluator-batch statistics.

### Phase 2: f64 batch substrate with scalar fast path

Deliverables:
- internal f64 sub-batching in the batch raw entry points,
- exact scalar fallback for `gammaloop_batch_size <= 1`,
- reusable queue arenas,
- f64 batch wrappers over Symbolica's batch evaluator API,
- request-weighted timing and batch-size accounting in the substrate.

### Phase 3: batched f64 orchestration through the existing process stack

Deliverables:
- batched f64 versions of:
  - `evaluate_from_source`,
  - `evaluate_stability_level`,
  - `evaluate_all_rotations`,
  - `evaluate_single`,
  - `evaluate_graph_term`,
- explicit handoff back to the untouched scalar path for escalated jobs.

### Phase 4: cross-section graph-term batching

Deliverables:
- `CrossSectionGraphTerm::evaluate_batch_f64(...)`, still reading like the current function,
- per-sample cut state containing the current staged continuation data,
- stage-local pass-1 and pass-2 evaluator queues,
- queue hooks shaped so LU threshold-subtraction jobs can be added later without another structural rewrite.

### Phase 5: amplitude graph-term batching

Deliverables:
- `AmplitudeGraphTerm::evaluate_batch_f64(...)`, reusing the same substrate,
- batched original integrand evaluator calls,
- batched amplitude threshold-counterterm evaluator calls.

### Phase 6: validation and regression checks

Correctness:
- compare batched f64 against scalar f64 with `batch_size = 1`,
- compare batched f64 against untouched scalar `f128` / `ArbPrec` escalation on targeted samples,
- verify event groups and observables are unchanged.

Performance:
- measure batch-of-one regression first,
- then measure hot evaluator occupancy and average batch size,
- only after that consider later backend vectorization work.

## Main weaknesses / culprits to keep in view

- The main settings subtlety is scope: `gammaloop_batch_size` must remain an entry-point integrator setting and must not silently turn into a downstream evaluator-queue cap.
- Timing normalization will be wrong if it is done per sample or per flush instead of per logical evaluator request.
- The new batch-size monitor will be misleading if it is flush-weighted or sample-count-weighted instead of request-weighted.
- `EvaluationMetaData` additions are cross-cutting because metadata is manually constructed in several places and summarized elsewhere.
- Cross-section threshold subtraction is still partially inactive at the call site, so the initial implementation should keep that path pluggable rather than pretending it is already live.
