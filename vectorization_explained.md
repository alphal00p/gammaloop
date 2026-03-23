# How GammaLoop Vectorization Works

This document explains how the f64 vectorized path was implemented in GammaLoop.

The goal is not just to list changed files, but to explain the design in the same order in which the code executes when a batched API entry point is called. The document is therefore organized as:

1. a high-level overview of the design philosophy,
2. a small mock-up that reproduces the same structural idea in simplified Rust,
3. a detailed walkthrough of the real implementation, with quoted code snippets and file locations.

Throughout this document, "vectorization" means: GammaLoop keeps several f64 samples alive at once, carries them through the normal evaluation pipeline, buffers work at evaluator boundaries, and eventually calls Symbolica's batch evaluator API with multiple inputs at once.

## 1. Overview

### 1.1 What problem was being solved?

GammaLoop already had public APIs that accepted batches:

- `evaluate_samples_raw(...)`
- `evaluate_momentum_configurations_raw(...)`

But the batch shape was destroyed immediately by scalar loops of the form:

```rust
for sample in samples {
    ...
}
```

This meant the rest of the evaluation pipeline was effectively scalar, even though the user had already provided a batch. As a result, the eventual evaluator call into Symbolica was also made sample by sample, which prevented the lower-level evaluator backend from seeing multiple inputs at once.

The refactor therefore had one central objective:

> keep several f64 samples alive together for as long as possible, and only collapse them when the algorithm truly becomes sample-specific.

In practice, the implementation pushes batching down to the evaluator boundary and only there issues a real batched evaluator call.

### 1.2 Design philosophy

The implementation does **not** try to rewrite GammaLoop into a SIMD-first algorithm.

Instead, it follows a more conservative design:

- keep the existing evaluation stack recognizable,
- keep the scalar f128 / arbitrary precision paths unchanged,
- keep early exits and sample-specific branching,
- only batch the work that naturally accumulates at evaluator call sites,
- materialize owned evaluator inputs so they can be buffered safely,
- flush those buffers when a stage is ready,
- scatter the results back into per-sample state and continue.

That is why the implementation feels like a "wavefront" or "queue-and-flush" design, but without introducing a global scheduler.

The main philosophy is:

1. keep the control flow local,
2. keep sample state explicit,
3. only batch where it actually pays off,
4. preserve the meaning of the timing and metadata that existed before.

This is important because GammaLoop has a lot of logic that cannot or should not be vectorized:

- parameterization may differ by sample,
- stability logic may escalate only some samples,
- selectors may reject only some events,
- threshold-subtraction related branches may activate only for some samples,
- cross-section evaluation has a true two-pass dependency: pass two can only be built after pass one has returned.

So the implementation keeps those branches, but it prevents them from forcing an early collapse back to purely scalar execution.

### 1.3 The core idea: buffer work, not full control flow

The easiest way to summarize the implementation is:

> each sample still walks through the normal GammaLoop pipeline, but whenever the code reaches a Symbolica evaluator boundary, the evaluator inputs are buffered into owned queues, then flushed in batch.

This gives a natural split:

- **sample-local work** remains local:
  - parameterization,
  - channel reinterpretation,
  - Newton solve,
  - event generation,
  - selector processing,
  - assembling pass-two parameters,
  - scattering results back into final event weights.

- **evaluator work** becomes batched:
  - original integrand evaluators,
  - summed / summed-function-map evaluator stacks,
  - amplitude counterterm prefactor evaluators,
  - amplitude counterterm esurface evaluators,
  - cross-section pass-one evaluators,
  - cross-section pass-two evaluators.

This is why the implementation introduces several "prepared request" structures: they are the owned payloads that cross the batching boundary.

### 1.4 Rebatching happens only at the API entry point

The user may call the API with a very large batch, for example 10,000 samples. The implementation does **not** keep all of those alive at once inside the batching logic. Instead, it rechunks them into smaller internal sub-batches controlled by:

- `runtime.integrator.gammaloop_batch_size`

So if:

- the API receives 10,000 samples,
- `gammaloop_batch_size = 100`,

GammaLoop processes 100-sample internal chunks.

This is important for two reasons:

1. it keeps the live buffers bounded,
2. it is still large enough that many evaluator queues naturally accumulate multiple requests and can benefit from batched backend execution.

It is also important that `gammaloop_batch_size` is **not** used as a hard cap on every internal queue. The queue sizes are allowed to grow according to the real branching structure. For example:

- a sample may fan out over multiple graphs,
- a graph may fan out over multiple channels,
- a cross-section cut may produce several pass-one or pass-two requests.

So the setting controls only the **entry-point rechunking**, not downstream queue lengths.

### 1.5 Why owned buffers were necessary

In the scalar path, many evaluator inputs are built on top of local or reused parameter buffers. That is fine when the evaluator is called immediately. It is not fine when the call is deferred.

If a sample prepares parameters, then the code moves on and prepares another sample's parameters into the same mutable storage, the first sample's evaluator input is gone.

So the batched path introduces a deliberate boundary:

- borrow-heavy, mutable construction on one side,
- owned, materialized evaluator inputs on the other side.

This is the role of `MaterializedInputParams<T>`:

- copy the current evaluator input into owned storage,
- retain orientation-related metadata,
- allow later modification such as orientation application or `override_if` toggling,
- survive until the queue is finally flushed.

This is one of the central enabling changes in the design.

### 1.6 The two-pass system in cross-section evaluation

The cross-section path is the trickiest part of the vectorization.

Its structure is not:

1. build parameters,
2. call evaluator,
3. done.

It is instead:

1. sample-local work determines whether a given threshold/cut branch is active,
2. pass one is evaluated,
3. pass one output is used to assemble pass-two inputs,
4. pass two is evaluated,
5. pass two output is combined with prefactors, counterterms, event weights, and accumulated cut results.

This means cross-section evaluation cannot simply queue "one request per sample".

It needs a staged queueing model:

- store per-sample continuation state,
- build `PassOneRequest`s,
- flush pass one,
- transform pass-one outputs into `PassTwoRequest`s,
- flush pass two,
- write results back into the stored sample state.

That is exactly what the implementation does.

The important point is that the control flow is still easy to recognize:

- the sample-local part of cut processing still happens in the cut loop,
- only the evaluator boundaries are batched,
- after a batched flush, execution resumes in normal GammaLoop logic.

### 1.7 How the evaluator batching is eventually realized

At the bottom of the stack, the actual "vectorized" call is not magical. It is a normal batch call:

1. flatten all queued inputs into one `flat_inputs` array,
2. allocate `flat_out`,
3. call the compiled evaluator's `evaluate_batch(...)`,
4. split `flat_out` back into per-request chunks,
5. decode those chunks into GammaLoop's `DualOrNot<Complex<...>>` representation.

This is performed in `GenericEvaluator::evaluate_batch_f64(...)`.

Above that layer:

- `EvaluatorStack::evaluate_batch_f64(...)` handles orientation expansion and evaluator-method-specific batching,
- graph-term batching feeds those evaluator batches,
- top-level batch orchestration feeds graph-term batches.

So the whole implementation is best seen as a staged preparation pipeline that gradually converts:

`API batch -> internal sample chunk -> per-graph work queues -> per-evaluator input queues -> Symbolica batch call`

### 1.8 Timing semantics were treated as part of the design

This was not only a correctness refactor. It also had to preserve the meaning of the existing timing fields.

If batching becomes more efficient, the timing reported **per sample** should improve accordingly. That means the metadata must not simply report wall-clock time of a shared batch as if each sample paid it in full.

The solution is a per-sample timing ledger:

- direct scalar work adds directly to the owning sample,
- shared evaluator flush time is divided over the logical requests in the flush,
- shared orchestration is apportioned back to the relevant owners,
- non-primary stability rotations contribute only to `total_timing`,
- primary-rotation evaluator work contributes to the evaluator and integrand timing metrics,
- the average evaluator batch size is tracked in the same request-weighted way.

This detail matters because the batching should not only be fast; its performance reporting should remain meaningful.

### 1.9 One small example

Suppose the user calls:

```text
evaluate_samples_raw(samples = 5 points, gammaloop_batch_size = 3)
```

Then GammaLoop internally processes:

- chunk A: samples `[0, 1, 2]`
- chunk B: samples `[3, 4]`

Inside chunk A:

- each sample is parameterized,
- each sample may fan out over graphs,
- the cross-section path may further fan out into pass-one and pass-two requests,
- all requests for the same evaluator stage are buffered together,
- the evaluator is called once with the full queue,
- results are scattered back to the owning samples.

So a sample is not "evaluated as a vector". Instead:

- its evaluator-boundary work participates in a vectorized batch,
- its sample-specific logic remains sample-specific.

That is the central structural idea of the implementation.

## 2. Simplified mock-up

The real GammaLoop code is complex because it has:

- rotations,
- precision escalation,
- multiple graph types,
- channels,
- orientations,
- events,
- threshold counterterms,
- dual numbers.

The simplified mock-up below removes all physics details but keeps the same structure:

- samples are received in batches,
- samples are rechunked,
- sample-local work builds owned requests,
- a two-pass system is used,
- evaluator calls are issued in batch,
- results are written back into per-sample state.

### 2.1 A scalar version

This is the shape GammaLoop used to have at a high level:

```rust
struct Sample {
    value: f64,
}

struct ResultValue {
    answer: f64,
}

struct Evaluator;

impl Evaluator {
    fn evaluate(&mut self, input: &[f64]) -> f64 {
        input.iter().sum()
    }
}

fn evaluate_samples_scalar(
    evaluator: &mut Evaluator,
    samples: &[Sample],
) -> Vec<ResultValue> {
    let mut out = Vec::with_capacity(samples.len());

    for sample in samples {
        let pass_one_input = vec![sample.value, 1.0];
        let pass_one = evaluator.evaluate(&pass_one_input);

        let pass_two_input = vec![pass_one, 2.0];
        let pass_two = evaluator.evaluate(&pass_two_input);

        out.push(ResultValue { answer: pass_two });
    }

    out
}
```

Even though `samples` is a slice, the batching disappears immediately.

### 2.2 A vectorized version with the same structure as GammaLoop

```rust
#[derive(Clone)]
struct Sample {
    value: f64,
}

#[derive(Default)]
struct SampleState {
    result: f64,
}

#[derive(Clone)]
struct PassOneRequest {
    sample_index: usize,
    input: Vec<f64>,
}

#[derive(Clone)]
struct PassTwoRequest {
    sample_index: usize,
    input: Vec<f64>,
}

struct Evaluator;

impl Evaluator {
    fn evaluate_batch(&mut self, inputs: &[Vec<f64>]) -> Vec<f64> {
        inputs.iter().map(|input| input.iter().sum()).collect()
    }
}

fn evaluate_samples_batched(
    evaluator: &mut Evaluator,
    samples: &[Sample],
    batch_size: usize,
) -> Vec<SampleState> {
    let mut final_results = Vec::with_capacity(samples.len());

    for chunk in samples.chunks(batch_size.max(1)) {
        let mut states = vec![SampleState::default(); chunk.len()];

        // Stage 1: sample-local work builds owned pass-one requests.
        let mut pass_one_requests = Vec::new();
        for (sample_index, sample) in chunk.iter().enumerate() {
            pass_one_requests.push(PassOneRequest {
                sample_index,
                input: vec![sample.value, 1.0],
            });
        }

        // Flush stage 1 in batch.
        let pass_one_inputs: Vec<Vec<f64>> = pass_one_requests
            .iter()
            .map(|req| req.input.clone())
            .collect();
        let pass_one_outputs = evaluator.evaluate_batch(&pass_one_inputs);

        // Stage 2: use pass-one outputs to build pass-two requests.
        let mut pass_two_requests = Vec::new();
        for (request, pass_one_output) in pass_one_requests.iter().zip(pass_one_outputs) {
            pass_two_requests.push(PassTwoRequest {
                sample_index: request.sample_index,
                input: vec![pass_one_output, 2.0],
            });
        }

        // Flush stage 2 in batch.
        let pass_two_inputs: Vec<Vec<f64>> = pass_two_requests
            .iter()
            .map(|req| req.input.clone())
            .collect();
        let pass_two_outputs = evaluator.evaluate_batch(&pass_two_inputs);

        // Scatter results back into sample-local state.
        for (request, output) in pass_two_requests.iter().zip(pass_two_outputs) {
            states[request.sample_index].result = output;
        }

        final_results.extend(states);
    }

    final_results
}
```

### 2.3 What this mock-up corresponds to in GammaLoop

- `samples.chunks(batch_size)` corresponds to internal rechunking by `gammaloop_batch_size`.
- `SampleState` corresponds to the real stored per-sample state such as:
  - `BatchedSourceState`,
  - `CrossSectionSampleState`,
  - `CutBatchState`.
- `PassOneRequest` and `PassTwoRequest` directly mirror the cross-section implementation.
- `evaluate_batch(...)` corresponds to `GenericEvaluator::evaluate_batch_f64(...)`.
- building `Vec<f64>` inputs corresponds to `MaterializedInputParams<f64>` and other owned parameter vectors.
- the final scatter back into `states[...]` is the same as the real code merging evaluator outputs back into graph results, event groups, cut totals, and metadata.

So the real implementation is more elaborate, but the structural pattern is exactly this one.

## 3. Detailed walkthrough of the real implementation

This section follows the order in which information flows when the public raw batch API is called.

### 3.1 Step 0: runtime setting that controls internal rechunking

Location: `crates/gammalooprs/src/settings/runtime.rs`

```rust
    #[serde(skip_serializing_if = "is_usize::<1>")]
    pub gammaloop_batch_size: usize,
```

```rust
            seed: 69,
            gammaloop_batch_size: 1,
            observables_output: ObservablesOutputSettings::default(),
```

This setting does not mean "force every evaluator queue to size N".

It only means:

- split the incoming public batch into internal sub-batches of this size,
- if the setting is `1`, behave as scalar,
- if the setting is larger than `1`, attempt batched f64 evaluation when the runtime mode supports it.

That design keeps the control knob simple and local.

### 3.2 Step 1: the public raw batch API decides whether to stay scalar or enter the batched f64 path

Location: `crates/gammalooprs/src/integrands/process/mod.rs`

```rust
        let mut results = Vec::with_capacity(samples.len());
        let batch_size = self.get_settings().integrator.gammaloop_batch_size.max(1);
        let batching_enabled = batch_size > 1
            && match self {
                ProcessIntegrand::Amplitude(integrand) => {
                    supports_batched_f64_evaluation(integrand)
                }
                ProcessIntegrand::CrossSection(integrand) => {
                    supports_batched_f64_evaluation(integrand)
                }
            };
        if !batching_enabled {
            for sample in samples {
                ...
            }
        } else {
            for batch in samples.chunks(batch_size) {
                let batch_results = match self {
                    ProcessIntegrand::Amplitude(integrand) => evaluate_sources_batched_f64(
                        integrand,
                        model,
                        &batch
                            .iter()
                            .map(|sample| (EvaluationSource::XSpace(sample), sample.get_weight()))
                            .collect::<Vec<_>>(),
                        use_arb_prec,
                        max_eval.clone(),
                    ),
                    ProcessIntegrand::CrossSection(integrand) => evaluate_sources_batched_f64(
                        integrand,
                        model,
                        &batch
                            .iter()
                            .map(|sample| (EvaluationSource::XSpace(sample), sample.get_weight()))
                            .collect::<Vec<_>>(),
                        use_arb_prec,
                        max_eval.clone(),
                    ),
                }?;
                ...
            }
        }
```

This is the first structural change:

- the old scalar path still exists,
- the new batched path is used only when the batch size is larger than one and the runtime mode supports batching.

The support check is intentionally explicit:

```rust
fn supports_batched_f64_evaluation<I: ProcessIntegrandImpl>(integrand: &I) -> bool {
    !integrand.get_settings().general.enable_cache
        && !matches!(
            integrand.get_settings().general.evaluator_method,
            evaluators::EvaluatorMethod::Iterative
        )
}
```

So the batching logic is not forced into paths where it would either be wrong or not worth it.

The same rebatching logic is also applied to `evaluate_momentum_configurations_raw(...)`.

### 3.3 Step 2: the batched path creates a per-sample state object instead of immediately finishing each sample

Location: `crates/gammalooprs/src/integrands/process/mod.rs`

```rust
fn evaluate_sources_batched_f64<'a, I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    sources: &[(EvaluationSource<'a>, F<f64>)],
    use_arb_prec: bool,
    max_eval: Complex<F<f64>>,
) -> Result<Vec<EvaluationResult>> {
    #[derive(Clone)]
    struct BatchedSourceState<'a> {
        source_index: usize,
        source: EvaluationSource<'a>,
        wgt: F<f64>,
        gammaloop_sample: GammaLoopSample<f64>,
        parameterization_time: Duration,
        stability_level: StabilityLevelSetting,
        total_levels: usize,
        loop_momenta_escalation: Option<LoopMomentaEscalationMetrics>,
    }
```

This `BatchedSourceState` is the first real queue-enabling structure.

Instead of:

- building a sample,
- evaluating it,
- immediately returning its final `EvaluationResult`,

the code now builds a batch of sample states that stay alive together.

The important design decision here is that only the f64-capable first stability level enters this path:

```rust
        if stability_level.precision != Precision::Double {
            results[source_index] = Some(evaluate_from_source(
                integrand,
                model,
                source,
                wgt,
                use_arb_prec,
                max_eval.clone(),
            )?);
            continue;
        }
```

So:

- f64 gets the batched substrate,
- escalated precision still uses the scalar reference implementation.

That was an explicit design choice, not an omission.

### 3.4 Step 3: timing and metadata are tracked per sample even though work is shared

Location: `crates/gammalooprs/src/integrands/process/mod.rs`

```rust
#[derive(Debug, Default, Clone)]
pub(crate) struct BatchedSampleTiming {
    total_seconds: f64,
    integrand_seconds: f64,
    evaluator_seconds: f64,
    evaluator_batch_size_sum: usize,
    evaluator_batch_size_count: usize,
}
```

```rust
pub(crate) fn record_shared_evaluator_flush(
    timings: &mut [BatchedSampleTiming],
    owner_indices: &[usize],
    shared_duration: Duration,
) {
    if owner_indices.is_empty() {
        return;
    }

    let batch_size = owner_indices.len();
    for &owner_index in owner_indices {
        timings[owner_index].record_shared_evaluator_flush(batch_size, shared_duration);
    }
}

pub(crate) fn record_shared_integrand_time(
    timings: &mut [BatchedSampleTiming],
    owner_indices: &[usize],
    shared_duration: Duration,
) {
    if owner_indices.is_empty() {
        return;
    }

    let shared_seconds = shared_duration.as_secs_f64() / owner_indices.len() as f64;
    for &owner_index in owner_indices {
        timings[owner_index].add_direct_integrand_seconds(shared_seconds);
    }
}
```

This timing ledger is not incidental bookkeeping. It is part of the batching design.

Why?

Because once work is shared, a simple wall-clock span does not mean the same thing anymore.

For example, if one evaluator flush handles 16 requests together:

- the total flush took one wall-clock duration,
- but 16 samples participated in it,
- and the metadata should reflect that shared cost per logical request.

The metadata structure was extended accordingly.

Location: `crates/gammalooprs/src/integrands/evaluation.rs`

```rust
pub struct EvaluationMetaData {
    pub total_timing: Duration,
    pub integrand_evaluation_time: Duration,
    pub evaluator_evaluation_time: Duration,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub average_evaluator_batch_size: Option<f64>,
    pub parameterization_time: Duration,
    pub event_processing_time: Duration,
    pub generated_event_count: usize,
    pub accepted_event_count: usize,
    pub relative_instability_error: Complex<F<f64>>,
    pub is_nan: bool,
    pub loop_momenta_escalation: Option<LoopMomentaEscalationMetrics>,
    pub stability_results: Vec<StabilityResult>,
    #[serde(skip_serializing)]
    pub(crate) evaluator_batch_size_sum: usize,
    #[serde(skip_serializing)]
    pub(crate) evaluator_batch_size_count: usize,
}
```

And the integration snapshot now exposes the aggregated monitor:

```rust
    pub(crate) fn snapshot(&self) -> IntegrationStatisticsSnapshot {
        IntegrationStatisticsSnapshot {
            num_evals: self.num_evals,
            average_total_time_seconds: self.get_avg_total_timing().as_secs_f64(),
            average_parameterization_time_seconds: self.get_avg_param_timing().as_secs_f64(),
            average_integrand_time_seconds: self.get_avg_integrand_timing().as_secs_f64(),
            average_evaluator_time_seconds: self.get_avg_evaluator_timing().as_secs_f64(),
            average_evaluator_batch_size: self.get_avg_evaluator_batch_size(),
            ...
        }
    }
```

### 3.5 Step 4: rotations are still handled, but only the primary rotation contributes to the batch metrics

Location: `crates/gammalooprs/src/integrands/process/mod.rs`

```rust
fn evaluate_all_rotations_batch_f64<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    gammaloop_samples: &[GammaLoopSample<f64>],
    timings: &mut [BatchedSampleTiming],
    record_rotated_results: bool,
) -> Result<(
    Vec<Vec<GraphEvaluationResult<f64>>>,
    usize,
    Vec<Vec<RotatedEvaluation>>,
)> {
    let rotations = integrand.get_rotations().cloned().collect_vec();
    let primary_rotation_index = rotations
        .iter()
        .position(Rotation::is_identity)
        .unwrap_or(0);
    ...
    for (rotation_index, rotation) in rotations.iter().enumerate() {
        let use_primary_timing = rotation_index == primary_rotation_index;
        let mut secondary_timings_storage = (!use_primary_timing)
            .then(|| vec![BatchedSampleTiming::default(); gammaloop_samples.len()]);
        let timing_target: &mut [BatchedSampleTiming] = match secondary_timings_storage.as_mut() {
            Some(secondary_timings) => secondary_timings.as_mut_slice(),
            None => timings,
        };
        ...
        let results =
            evaluate_single_batch_f64(integrand, model, &rotated_samples, rotation, timing_target)?;

        if let Some(secondary_timings) = secondary_timings_storage.as_ref() {
            for (timing, secondary_timing) in timings.iter_mut().zip(secondary_timings.iter()) {
                timing.add_total_only_seconds(secondary_timing.total_seconds);
            }
        }
        ...
    }
```

This is a subtle but important piece of the design.

Scalar GammaLoop already had the rule:

- the primary rotation is the one that contributes to the integrand/evaluator timing metrics,
- secondary rotations are used for stability logic but should not inflate the primary evaluator metrics.

The batched implementation preserves that rule by:

- using the real timing ledger for the primary rotation,
- using scratch ledgers for non-primary rotations,
- merging non-primary work back only into `total_timing`.

So the physics logic sees all rotations, but the performance metrics keep their scalar interpretation.

### 3.6 Step 5: at the graph boundary, sample-local work turns into queued graph requests

Location: `crates/gammalooprs/src/integrands/process/mod.rs`

```rust
pub(crate) struct GraphTermBatchRequest<'a> {
    pub sample_index: usize,
    pub sample: &'a MomentumSample<f64>,
    pub channel_id: Option<(ChannelIndex, F<f64>)>,
}
```

And the `GraphTerm` trait gained a batched method:

```rust
pub(crate) trait GraphTerm {
    fn evaluate<T: FloatLike>(...) -> Result<GraphEvaluationResult<T>>;
    ...
    fn evaluate_batch_f64(
        &mut self,
        requests: &[GraphTermBatchRequest<'_>],
        model: &Model,
        settings: &RuntimeSettings,
        event_processing_runtime: Option<&mut EventProcessingRuntime>,
        rotation: &Rotation,
        timings: &mut [BatchedSampleTiming],
    ) -> Result<Vec<GraphEvaluationResult<f64>>>;
}
```

This is one of the key signature changes in the implementation.

Before, the graph term API was essentially scalar:

- one sample in,
- one graph result out.

After the refactor, the f64 path also supports:

- many graph requests in,
- one graph result per request out.

The function that actually fills those graph requests is `evaluate_single_batch_f64(...)`.

Location: `crates/gammalooprs/src/integrands/process/mod.rs`

```rust
fn evaluate_single_batch_f64<I: ProcessIntegrandImpl>(
    integrand: &mut I,
    model: &Model,
    gammaloop_samples: &[GammaLoopSample<f64>],
    rotation: &Rotation,
    timings: &mut [BatchedSampleTiming],
) -> Result<Vec<GraphEvaluationResult<f64>>> {
    #[derive(Clone, Copy)]
    struct GraphTask<'a> {
        batch_index: usize,
        sample_index: usize,
        sample: &'a MomentumSample<f64>,
        channel_id: Option<(ChannelIndex, F<f64>)>,
        merge_into_group: bool,
        tropical_prefactor: Option<F<f64>>,
    }
    ...
```

This function does not yet batch evaluator calls directly. Its role is earlier:

- fan out each sample into graph-specific work,
- store that work in `tasks_by_graph`,
- then call `evaluate_graph_term_batch_f64(...)` once per graph with all requests for that graph.

That means the batching hierarchy is:

1. top-level batch of samples,
2. grouped per graph,
3. graph term internally grouped per evaluator stage.

The graph-level buffering is also where shared orchestration is timed:

```rust
        let owner_indices = tasks
            .iter()
            .map(|task| task.batch_index)
            .collect::<Vec<_>>();
        let request_build_start = Instant::now();
        let requests = tasks
            .iter()
            .map(|task| GraphTermBatchRequest {
                sample_index: task.sample_index,
                sample: task.sample,
                channel_id: task.channel_id,
            })
            .collect::<Vec<_>>();
        record_shared_integrand_time(timings, &owner_indices, request_build_start.elapsed());
```

And the results are later merged back into the correct owning sample:

```rust
        let merge_start = Instant::now();
        for (task, mut result) in tasks.into_iter().zip(results) {
            if let Some(prefactor) = task.tropical_prefactor {
                result.integrand_result *= Complex::new_re(prefactor);
            }

            if task.merge_into_group {
                for mut event_group in result.event_groups.drain(..) {
                    grouped_events[task.batch_index].append(&mut event_group);
                }
            }
            accumulators[task.batch_index].merge_in_place(result);
        }
        record_shared_integrand_time(timings, &owner_indices, merge_start.elapsed());
```

So the graph layer is where "sample batch" turns into "work queues keyed by graph".

### 3.7 Step 6: owned materialized inputs make deferred evaluator calls safe

Location: `crates/gammalooprs/src/integrands/process/evaluators.rs`

```rust
#[derive(Clone)]
pub struct MaterializedInputParams<T: FloatLike> {
    pub values: Vec<Complex<F<T>>>,
    pub orientations_start: usize,
    pub override_pos: usize,
    pub multiplicative_offset: usize,
}

impl<T: FloatLike> MaterializedInputParams<T> {
    pub fn from_input(mut input: InputParams<'_, T>) -> Self {
        Self {
            values: input.as_mut().to_vec(),
            orientations_start: input.orientations_start,
            override_pos: input.override_pos,
            multiplicative_offset: input.multiplicative_offset,
        }
    }

    pub fn as_slice(&self) -> &[Complex<F<T>>] {
        &self.values
    }

    pub fn apply_orientation<O: GraphOrientation>(&mut self, orientation: &O) {
        ...
    }

    pub fn set_override_if(&mut self, over_ride: bool) {
        ...
    }
}
```

This type is the exact boundary between:

- "I just built parameters and could evaluate immediately",
- "I want to defer evaluation and batch it with other requests later".

Without this owned copy, the queueing system would be unsafe or wrong because later parameter builds would overwrite earlier ones.

The extra metadata (`orientations_start`, `override_pos`, `multiplicative_offset`) is what lets the batched path still apply:

- orientation-specific edits,
- `override_if` changes for summed evaluator modes,

after the input has been materialized.

### 3.8 Step 7: `EvaluatorStack` batches inputs while preserving orientation semantics

Location: `crates/gammalooprs/src/integrands/process/evaluators.rs`

```rust
pub(crate) fn evaluate_batch_f64(
    &mut self,
    inputs: &[MaterializedInputParams<f64>],
    orientations: &[OwnedOrientations<OrientationID>],
    all_orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
    filter: &SubSet<OrientationID>,
    settings: &RuntimeSettings,
) -> Result<Vec<Vec<DualOrNot<Complex<F<f64>>>>>> {
    debug_assert_eq!(inputs.len(), orientations.len());
    match settings.general.evaluator_method {
        EvaluatorMethod::SingleParametric => {
            let mut expanded_inputs = Vec::new();
            let mut owners = Vec::new();
            for (owner, (input, orientation_mode)) in
                inputs.iter().zip(orientations).enumerate()
            {
                for (_, orientation) in
                    Self::orientation_iter(*orientation_mode, all_orientations, filter)
                {
                    let mut expanded = input.clone();
                    expanded.apply_orientation(orientation);
                    expanded_inputs.push(expanded);
                    owners.push(owner);
                }
            }
            ...
        }
        ...
    }
}
```

This layer is responsible for a critical GammaLoop-specific detail:

- one logical request does not always correspond to one physical evaluator call,
- because orientations may have to be expanded.

In the single-parametric path, the code:

- clones the materialized input,
- applies the orientation,
- pushes the physical evaluator input into `expanded_inputs`,
- remembers which logical owner it belongs to,
- later sums the results back into the correct owner slot.

The summed paths are slightly different and have their own batched helper:

```rust
fn evaluate_summed_like_batch_f64(
    evaluator: &mut GenericEvaluator,
    inputs: &[MaterializedInputParams<f64>],
    orientations: &[OwnedOrientations<OrientationID>],
    all_orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
    filter: &SubSet<OrientationID>,
) -> Result<Vec<Vec<DualOrNot<Complex<F<f64>>>>>> {
    let mut direct_inputs = Vec::new();
    let mut direct_owners = Vec::new();
    let mut expanded_inputs = Vec::new();
    let mut expanded_owners = Vec::new();

    for (owner, (input, orientation_mode)) in inputs.iter().zip(orientations).enumerate() {
        if matches!(orientation_mode, OwnedOrientations::All) {
            direct_owners.push(owner);
            let mut direct_input = input.clone();
            direct_input.set_override_if(true);
            direct_inputs.push(direct_input);
        } else {
            ...
        }
    }
    ...
}
```

That `set_override_if(true)` matters because the scalar summed path already relied on that behavior. The batched path had to preserve it explicitly.

### 3.9 Step 8: `GenericEvaluator` is where the actual batch call to Symbolica happens

Location: `crates/gammalooprs/src/integrands/process/evaluators.rs`

```rust
pub(crate) fn evaluate_batch_f64(
    &mut self,
    inputs: &[MaterializedInputParams<f64>],
) -> Vec<Vec<DualOrNot<Complex<F<f64>>>>> {
    if inputs.is_empty() {
        return Vec::new();
    }

    let out_size = self.compute_out_size();
    if let Some(compiled) = &mut self.f64_compiled {
        let mut flat_inputs =
            Vec::with_capacity(inputs.iter().map(|input| input.values.len()).sum());
        for input in inputs {
            flat_inputs.extend_from_slice(input.as_slice());
        }

        let mut flat_out = vec![Complex::default(); inputs.len() * out_size];
        compiled
            .evaluate_batch(inputs.len(), &flat_inputs, &mut flat_out)
            .expect("compiled evaluator batch evaluation failed");

        flat_out
            .chunks(out_size)
            .map(|chunk| self.decode_f64_output(chunk))
            .collect()
    } else {
        ...
    }
}
```

This is the bottom of the batching stack:

- all higher layers are about preparing and grouping requests,
- this function is where those queued requests finally become one batched evaluator call.

The flow is very direct:

1. flatten all inputs,
2. allocate the output vector with the right batched size,
3. call `compiled.evaluate_batch(...)`,
4. split the output back per logical request,
5. decode it into GammaLoop's expected representation.

There is also a special helper for evaluators that are expected to produce one complex number per request:

```rust
pub(crate) fn evaluate_batch_single_f64(
    &mut self,
    inputs: &[Vec<Complex<F<f64>>>],
) -> Vec<Complex<F<f64>>> {
    if inputs.is_empty() {
        return Vec::new();
    }

    let out_size = self.compute_out_size();
    assert_eq!(
        out_size, 1,
        "evaluate_batch_single_f64 requires evaluators with a single complex output"
    );

    let input_len = inputs[0].len();
    assert!(
        inputs.iter().all(|input| input.len() == input_len),
        "evaluate_batch_single_f64 requires all inputs to have the same arity"
    );

    if let Some(compiled) = &mut self.f64_compiled {
        let mut flat_inputs = Vec::with_capacity(inputs.len() * input_len);
        for input in inputs {
            flat_inputs.extend_from_slice(input);
        }

        let mut flat_out = vec![Complex::default(); inputs.len() * out_size];
        compiled
            .evaluate_batch(inputs.len(), &flat_inputs, &mut flat_out)
            .expect("compiled evaluator batch evaluation failed");
        flat_out.into_iter().step_by(out_size).collect()
    } else {
        ...
    }
}
```

This helper is used in places like:

- cross-section pass two,
- amplitude counterterm prefactor evaluation,

where the output shape is known to be a single complex value per request.

### 3.10 Step 9: amplitude batching is the simpler graph-term case

Location: `crates/gammalooprs/src/integrands/process/amplitude/mod.rs`

```rust
fn evaluate_batch_f64(
    &mut self,
    requests: &[GraphTermBatchRequest<'_>],
    model: &Model,
    settings: &RuntimeSettings,
    _event_processing_runtime: Option<&mut EventProcessingRuntime>,
    rotation: &Rotation,
    timings: &mut [BatchedSampleTiming],
) -> Result<Vec<GraphEvaluationResult<f64>>> {
    if requests.is_empty() {
        return Ok(Vec::new());
    }

    #[derive(Clone)]
    struct PreparedAmplitudeRequest {
        sample_index: usize,
        prefactor: F<f64>,
        orientations: OwnedOrientations<OrientationID>,
        input: MaterializedInputParams<f64>,
        sample: MomentumSample<f64>,
    }
    ...
```

This function illustrates the overall pattern in its simplest real form:

1. do sample-local preparation,
2. materialize evaluator inputs,
3. collect all prepared requests,
4. flush the main integrand evaluator in batch,
5. build and flush threshold-counterterm requests,
6. merge everything back per sample.

The core flush is:

```rust
        let inputs = prepared
            .iter()
            .map(|request| request.input.clone())
            .collect::<Vec<_>>();
        let orientations = prepared
            .iter()
            .map(|request| request.orientations)
            .collect::<Vec<_>>();
        let owners = prepared
            .iter()
            .map(|request| request.sample_index)
            .collect::<Vec<_>>();

        let start = std::time::Instant::now();
        let bare_results = self.original_integrand.evaluate_batch_f64(
            &inputs,
            &orientations,
            &self.orientations,
            &self.orientation_filter,
            settings,
        )?;
        record_shared_evaluator_flush(timings, &owners, start.elapsed());
```

The structure is very regular:

- gather inputs,
- gather orientation modes,
- gather owner indices,
- batch call,
- apportion timing back to owners.

That same pattern repeats throughout the implementation.

### 3.11 Step 10: amplitude counterterms are batched by grouping prepared requests by downstream evaluator

Location: `crates/gammalooprs/src/subtraction/amplitude_counterterm.rs`

```rust
pub(crate) fn evaluate_batch_f64(
    &mut self,
    requests: &[AmplitudeCountertermBatchRequest<'_>],
    graph: &Graph,
    model: &Model,
    esurfaces: &EsurfaceCollection,
    all_orientations: &TiVec<OrientationID, EdgeVec<Orientation>>,
    orientation_filter: &linnet::half_edge::subgraph::subset::SubSet<OrientationID>,
    rotation: &Rotation,
    settings: &RuntimeSettings,
    param_builder: &mut ParamBuilder<f64>,
    timings: &mut [BatchedSampleTiming],
) -> Vec<Complex<F<f64>>> {
    if requests.is_empty() {
        return Vec::new();
    }
    ...
    let mut prepared_requests = Vec::new();
    for (request_position, request) in requests.iter().enumerate() {
        ...
        let prepared_request = esurface_builder
            .solve_rstar()
            .rstar_samples()
            .prepare_batch_request(
                request_position,
                request.sample_index,
                param_builder,
                request.orientations,
                group_index,
            );
        prepared_requests.push(prepared_request);
        ...
    }
```

This is conceptually similar to amplitude proper, but with another level of grouping:

- first prepare all possible counterterm requests,
- then group those requests by prefactor evaluator,
- flush each prefactor evaluator in batch,
- then group requests by esurface evaluator,
- flush each esurface evaluator in batch,
- accumulate the results back into the owning sample.

The grouping code makes that explicit:

```rust
        let mut prefactor_groups = vec![Vec::<usize>::new(); self.overlap.overlap_groups.len()];
        for (prepared_index, request) in prepared_requests.iter().enumerate() {
            if let Some((group_index, _)) = request.prefactor_queue {
                prefactor_groups[group_index].push(prepared_index);
            }
        }
```

And later:

```rust
        let mut esurface_groups = vec![Vec::<usize>::new(); self.evaluators.len()];
        for (prepared_index, request) in prepared_requests.iter().enumerate() {
            esurface_groups[request.esurface_id.0].push(prepared_index);
        }
```

This is a good example of the overall design philosophy:

- do not force all branching to disappear,
- accept the branching,
- buffer the resulting requests,
- regroup them by the actual downstream evaluator,
- then batch.

### 3.12 Step 11: cross-section batching uses a true staged queue with pass one and pass two

Location: `crates/gammalooprs/src/integrands/process/cross_section_integrand.rs`

```rust
fn evaluate_batch_f64(
    &mut self,
    requests: &[GraphTermBatchRequest<'_>],
    model: &Model,
    settings: &RuntimeSettings,
    mut event_processing_runtime: Option<&mut EventProcessingRuntime>,
    _rotation: &Rotation,
    timings: &mut [BatchedSampleTiming],
) -> Result<Vec<GraphEvaluationResult<f64>>> {
    if requests.is_empty() {
        return Ok(Vec::new());
    }

    #[derive(Clone)]
    struct CrossSectionSampleState {
        sample_index: usize,
        momentum_sample: MomentumSample<f64>,
        orientations: OwnedOrientations<OrientationID>,
        cut_results: Vec<Complex<F<f64>>>,
        cut_threshold_counterterms: Vec<Complex<F<f64>>>,
        differential_result: GraphEvaluationResult<f64>,
        accepted_event_group: GenericEventGroup<f64>,
    }

    struct CutBatchState {
        accepted_event: Option<GenericEvent<f64>>,
        threshold_counterterm_weights: Vec<Complex<F<f64>>>,
        bare_cut_total: Complex<F<f64>>,
    }

    #[derive(Clone)]
    struct PassOneRequest {
        request_position: usize,
        sample_index: usize,
        prefactor: Complex<F<f64>>,
        esurface_derivatives: DualOrNot<F<f64>>,
        input: MaterializedInputParams<f64>,
    }

    #[derive(Clone)]
    struct PassTwoRequest {
        request_position: usize,
        sample_index: usize,
        prefactor: Complex<F<f64>>,
        ct_result: Complex<F<f64>>,
        params: Vec<Complex<F<f64>>>,
    }
```

This is the clearest expression of the queue-and-resume design.

There are three layers of state:

- `CrossSectionSampleState`: what lives for the whole sample,
- `CutBatchState`: what lives across one cut while waiting for pass-two completion,
- `PassOneRequest` / `PassTwoRequest`: what is actually sent to evaluators.

The flow inside one cut is:

1. do sample-local work,
2. maybe generate/process an event,
3. maybe skip the sample entirely for this cut,
4. otherwise build pass-one requests,
5. flush pass one,
6. use pass-one outputs to build pass-two requests,
7. flush pass two,
8. write everything back into the saved cut/sample state.

What matters for event semantics is that event construction is **not** moved behind the batched evaluator call. It still happens up front, before the pass-one queue is even filled, exactly because selector failure must be able to terminate the sample early.

Location: `crates/gammalooprs/src/integrands/process/cross_section_integrand.rs`

```rust
                let mut accepted_event = event;
                let mut selectors_pass = true;
                if let Some(evt) = accepted_event.as_mut() {
                    sample_state.differential_result.generated_event_count += 1;
                    let event_start = Instant::now();
                    if let Some(runtime) = event_processing_runtime.as_deref_mut() {
                        if runtime.has_selectors() || runtime.has_observables() {
                            selectors_pass = runtime.process_event(evt);
                        }
                    }
                    event_timing += event_start.elapsed();
                    if !selectors_pass {
                        accepted_event = None;
                    }
                }
                sample_state.differential_result.event_processing_time += event_timing;

                if !selectors_pass {
                    let zero = Complex::new_re(momentum_sample.zero());
                    for _ in 1..=max_occurance {
                        sample_state.cut_threshold_counterterms.push(zero.clone());
                        sample_state.cut_results.push(zero.clone());
                    }
                    timings[sample_state.sample_index].add_direct_integrand_time(start.elapsed());
                    continue;
                }

                cut_states[request_position] = Some(CutBatchState {
                    accepted_event,
                    threshold_counterterm_weights: Vec::with_capacity(max_occurance),
                    bare_cut_total: Complex::new_re(momentum_sample.zero()),
                });
```

This is structurally very important:

- event generation still happens before evaluator batching,
- selector processing still happens immediately on the event object,
- rejected events still terminate the cut early through `continue`,
- and only surviving events are stored in `CutBatchState`.

So the batching layer does **not** weaken or delay selector rejection. It only carries the already-accepted event payload through the pass-one / pass-two continuation.

The queue fill for pass one is:

```rust
                    let params = MaterializedInputParams::from_input(
                        <f64 as crate::integrands::process::GenericEvaluatorFloat>::get_parameters(
                            &mut self.param_builder,
                            (settings.general.enable_cache, settings.general.debug_cache),
                            &self.graph,
                            &rescaled_momenta,
                            hel,
                            &settings.additional_params(),
                            None,
                            None,
                            Some(&lu_params),
                        ),
                    );

                    pass_one_by_count[num_esurfaces - 1].push(PassOneRequest {
                        request_position,
                        sample_index: sample_state.sample_index,
                        prefactor,
                        esurface_derivatives,
                        input: params,
                    });
```

Then the pass-one flush:

```rust
                let inputs = pass_one_requests
                    .iter()
                    .map(|request| request.input.clone())
                    .collect::<Vec<_>>();
                let orientations = pass_one_requests
                    .iter()
                    .map(|request| sample_states[request.request_position].orientations)
                    .collect::<Vec<_>>();
                let owners = pass_one_requests
                    .iter()
                    .map(|request| request.sample_index)
                    .collect::<Vec<_>>();

                let eval_start = Instant::now();
                let values = self.integrand[raised_cut][num_esurfaces - 1].evaluate_batch_f64(
                    &inputs,
                    &orientations,
                    &self.orientations,
                    &self.orientation_filter,
                    settings,
                )?;
                record_shared_evaluator_flush(timings, &owners, eval_start.elapsed());
```

Then pass-one output is transformed into pass-two inputs:

```rust
                    let mut params_for_pass_two = vec![];
                    match result {
                        DualOrNot::Dual(dual_result) => {
                            params_for_pass_two
                                .extend_from_slice(&extract_t_derivatives_complex(dual_result));
                        }
                        DualOrNot::NonDual(non_dual_result) => {
                            params_for_pass_two.push(non_dual_result);
                        }
                    }

                    match &request.esurface_derivatives {
                        DualOrNot::Dual(dual_e_surface) => {
                            extract_t_derivatives(dual_e_surface.clone())[1..]
                                .iter()
                                .for_each(|value| {
                                    params_for_pass_two.push(Complex::new_re(value.clone()));
                                });
                        }
                        DualOrNot::NonDual(non_dual_e_surface) => {
                            params_for_pass_two.push(Complex::new_re(non_dual_e_surface.clone()));
                        }
                    }

                    pass_two_by_count[num_esurfaces - 1].push(PassTwoRequest {
                        request_position: request.request_position,
                        sample_index: request.sample_index,
                        prefactor: request.prefactor.clone(),
                        ct_result,
                        params: params_for_pass_two,
                    });
```

And finally pass two is flushed:

```rust
                let params = pass_two_requests
                    .iter()
                    .map(|request| request.params.clone())
                    .collect::<Vec<_>>();
                let owners = pass_two_requests
                    .iter()
                    .map(|request| request.sample_index)
                    .collect::<Vec<_>>();

                let eval_start = Instant::now();
                let values = self.raised_data.pass_two_evaluators[num_esurfaces - 1]
                    .evaluate_batch_single_f64(&params);
                record_shared_evaluator_flush(timings, &owners, eval_start.elapsed());
```

This is the heart of the "two-pass system" that the overview referred to.

Notice what is **not** done:

- the implementation does not try to merge pass one and pass two into one super-batch,
- because pass two genuinely depends on pass-one outputs.

Instead, the implementation stores enough continuation state to resume naturally after the pass-one flush.

That is the elegant part of the design.

It is also the reason the event path remains faithful to the scalar implementation: the event is already known, and possibly already rejected, before the batched evaluator work starts.

### 3.13 Step 12: finalization converts the batched internal state back into ordinary `EvaluationResult`s

Returning to `evaluate_sources_batched_f64(...)`, once all rotations and graph work are done, the batched function finishes each sample one by one:

```rust
        for ((state, mut graph_results), mut timing) in batched_states
            .into_iter()
            .zip(rotation_results.into_iter())
            .zip(timings.into_iter())
        {
            let stability_start = Instant::now();
            let results_for_stability = graph_results
                .iter()
                .map(|result| result.integrand_result.clone())
                .collect_vec();
            let (average_result, estimated_relative_accuracy, is_stable, _reason) =
                if integrand.get_settings().stability.check_on_norm {
                    stability_check_on_norm(...)
                } else {
                    stability_check(...)
                };
            timing.add_total_only_time(stability_start.elapsed());

            if !is_stable {
                results[state.source_index] = Some(evaluate_from_source(...)?);
                continue;
            }
            ...
            let mut evaluation_metadata = EvaluationMetaData::new_empty();
            timing.apply_to_metadata(&mut evaluation_metadata);
            evaluation_metadata.parameterization_time = state.parameterization_time;
            evaluation_metadata.event_processing_time = event_processing_time;
            evaluation_metadata.generated_event_count = graph_result.generated_event_count;
            evaluation_metadata.accepted_event_count = graph_result.accepted_event_count;
            ...
            results[state.source_index] = Some(EvaluationResult {
                integrand_result: nanless_result,
                parameterization_jacobian,
                integrator_weight: state.wgt,
                event_groups,
                evaluation_metadata,
            });
        }
```

However, it is worth being explicit about **what** is being finalized here.

At this point, the cross-section graph-term batching has already reconstructed the exact retained event list for each sample. The accepted event is resumed from `CutBatchState`, receives its final cut-dependent weight after pass two, and is pushed into the sample-local `accepted_event_group` only if it survived selectors.

Location: `crates/gammalooprs/src/integrands/process/cross_section_integrand.rs`

```rust
            for (request_position, sample_state) in sample_states.iter_mut().enumerate() {
                let Some(cut_state) = cut_states[request_position].take() else {
                    continue;
                };

                if let Some(mut event) = cut_state.accepted_event {
                    let threshold_counterterm_total =
                        cut_state.threshold_counterterm_weights.iter().fold(
                            Complex::new_re(sample_state.momentum_sample.zero()),
                            |acc, value| acc + value.clone(),
                        );
                    event.weight = cut_state.bare_cut_total.clone() + threshold_counterterm_total;

                    if settings.general.store_additional_weights_in_event {
                        event.additional_weights.weights.insert(
                            AdditionalWeightKey::Original,
                            cut_state.bare_cut_total.clone(),
                        );
                        for (subset_index, threshold_counterterm) in cut_state
                            .threshold_counterterm_weights
                            .into_iter()
                            .enumerate()
                        {
                            event.additional_weights.weights.insert(
                                AdditionalWeightKey::ThresholdCounterterm { subset_index },
                                threshold_counterterm,
                            );
                        }
                    }

                    sample_state.differential_result.accepted_event_count += 1;
                    if settings.should_buffer_generated_events() {
                        sample_state.accepted_event_group.push(event);
                    }
                }
            }
```

Only after that does the graph-term result move to the final per-sample assembly step, where the same full multiplicative factor logic as before is applied to the retained event groups:

Location: `crates/gammalooprs/src/integrands/process/mod.rs`

```rust
            let parameterization_jacobian = match state.source {
                EvaluationSource::XSpace(_) => {
                    Some(state.gammaloop_sample.get_default_sample().jacobian())
                }
                EvaluationSource::Momentum(_) => None,
            };
            let full_factor =
                full_event_multiplicative_factor(parameterization_jacobian, state.wgt);
            let mut event_groups = graph_result.event_groups;
            apply_full_event_multiplicative_factor(&mut event_groups, &full_factor);
            timing.add_total_only_time(finalization_start.elapsed());
```

So the event lifecycle through batching is:

1. generate the event sample-locally,
2. run selectors immediately,
3. if it fails, drop it immediately and skip any later event retention,
4. if it survives, store it inside the per-cut continuation state,
5. after pass two, attach the final cut weight and optional additional weights,
6. move it into the sample-local accepted event group,
7. during final sample assembly, apply the same overall multiplicative factor as in the scalar path,
8. return it inside the ordinary `EvaluationResult.event_groups`.

That is the precise sense in which batching preserves event construction "as before": batching only delays the moment when the **final numerical weight** becomes available; it does not delay selector rejection, and it does not change which events are retained.

This is where the implementation becomes "ordinary GammaLoop" again:

- each sample gets one final `EvaluationResult`,
- the batch disappears,
- the output ordering matches the input ordering,
- upstream callers do not need to know about the internal queueing.

That is another important design property:

> the batching is deep inside the execution pipeline, but the public API contract stays ordinary.

### 3.14 Step 13: validation checks the real behavior, not just compilation

Location: `tests/tests/test_runs.rs`

```rust
fn test_batched_evaluation() -> Result<()> {
    let gg_hhh_test_name = "test_batched_evaluation_gg_hhh";
    let gg_hhh_base = setup_gg_hhh_batched_evaluation_cli(gg_hhh_test_name)?;
    let gg_hhh_points = vec![default_xspace_point(&gg_hhh_base)?; 5];
    let mut gg_hhh_scalar = gg_hhh_base.clone();
    let mut gg_hhh_batched = gg_hhh_base;

    gg_hhh_scalar
        .run_command("set process -p gg_hhh -i 1L kv integrator.gammaloop_batch_size=1")?;
    gg_hhh_batched
        .run_command("set process -p gg_hhh -i 1L kv integrator.gammaloop_batch_size=3")?;

    let gg_hhh_scalar_result = evaluate_samples_api(&mut gg_hhh_scalar, &gg_hhh_points, false)?;
    let gg_hhh_batched_result = evaluate_samples_api(&mut gg_hhh_batched, &gg_hhh_points, false)?;

    assert_batch_evaluation_results_match_ignoring_timings(
        &gg_hhh_scalar_result,
        &gg_hhh_batched_result,
    );
    ...
}
```

This test is important because it checks more than "the code runs".

It checks:

- scalar and batched results match,
- the batch monitor is `1.0` in scalar mode,
- the batch monitor exceeds `1.0` in batched mode,
- both an amplitude and a differential cross-section style process are covered,
- generated event and observable behavior survives the refactor.

In other words, the implementation was validated at the level of:

- numerical agreement,
- event/observable behavior,
- and the new monitoring metadata.

## 4. Final mental model

If you want a one-paragraph summary of the whole implementation, this is the right one:

When `evaluate_samples_raw(...)` is called, GammaLoop first rechunks the public batch into internal f64-sized chunks. Each sample in a chunk is converted into persistent per-sample state. Those states are marched through the usual GammaLoop logic. Whenever the code reaches an evaluator boundary, it materializes owned inputs, buffers them into request queues, and flushes the appropriate evaluator in batch. In the cross-section path, this happens in two stages because pass two depends on pass-one outputs. After each flush, the results are scattered back into the stored sample state, and the normal algorithm continues. Timing and batch-size monitoring are accumulated per sample throughout the process, so the public metadata keeps the same meaning as before. Finally, the batched internal state is reassembled into the ordinary `EvaluationResult` objects returned by the existing API.

That is the essence of how GammaLoop's vectorization was achieved.
