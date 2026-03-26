# gammaLoop Current Architecture

## Scope
This document describes the current, implemented architecture of this repository.

The workspace is organized around two main Rust crates:
- `gammalooprs` (`crates/gammalooprs`): core physics/domain logic, graph processing, integrand construction, evaluation, and integration.
- `gammaloop-api` (`crates/gammaloop-api`): CLI, REPL, Python bindings, command parsing, and persisted state orchestration.

## High-Level Architecture
At a high level, gammaLoop uses a layered architecture with a stateful application shell.

```text
User / Automation
  -> CLI binary (`gammaloop`) and REPL
  -> Python module (`_gammaloop` via pyo3)
  -> run cards (`run.toml`)

Application Layer (`gammaloop-api`)
  -> clap command parsing
  -> command execution against mutable `State`
  -> state/load/save and logging control

Domain Layer (`gammalooprs`)
  -> model + parameter cards
  -> graph and process representations (amplitude/cross-section)
  -> preprocessing and CFF generation
  -> GL integrand construction + evaluator stacks
  -> differential event generation, selectors, and observables
  -> Monte Carlo integration and stability checks

Infrastructure
  -> filesystem persistence (`gammaloop_state/`)
  -> optional local scratch persistence (`.local/scratch/...`)
  -> Symbolica state import/export
  -> tracing/logging
```

## Main Components

### 1. Entry Points and Interfaces
- CLI entry point: `crates/gammaloop-api/src/cli.rs`.
- Main app orchestration: `crates/gammaloop-api/src/lib.rs` (`OneShot::load`, `OneShot::run`).
- Python module entry: `crates/gammaloop-api/src/python.rs` (`#[pymodule(name = "_gammaloop")]`).
- Python package shim: `crates/gammaloop-api/python/gammaloop/__init__.py`.

### 1.1 CLI architecture
- The CLI is state-centric: every session operates on one active state folder and one active in-memory `RunHistory`.
- Startup is orchestrated by `OneShot`:
  - resolve the state folder from CLI / boot card / defaults
  - load or create the state
  - load persisted settings and run history
  - optionally apply a boot card before entering the REPL or running a one-shot subcommand
- Command execution is centralized in `crates/gammaloop-api/src/session.rs` (`CliSession`):
  - top-level commands
  - boot-card replay
  - nested `run` execution
  - command-history recording policy
  - command-block-definition mode
- The REPL layer in `crates/gammaloop-api/src/repl.rs` is state-aware:
  - prompt text reflects the active state settings
  - tab completion is driven from the active state, active command blocks, and clap argument metadata

### 2. Application State and Command Model
- Central mutable app state: `crates/gammaloop-api/src/state.rs` (`State`).
- Command history and run cards: `CommandHistory`, `RunHistory`.
- Command dispatch: `crates/gammaloop-api/src/commands/mod.rs` with concrete subcommands in `commands/*`.

The command model is stateful by design: commands mutate a long-lived `State` that can be saved and resumed.

### 3. Domain Core (gammalooprs)
- Root module wiring: `crates/gammalooprs/src/lib.rs`.
- Model and parameters: `crates/gammalooprs/src/model/mod.rs`.
- Graph domain: `crates/gammalooprs/src/graph/mod.rs` and submodules.
- Process orchestration: `crates/gammalooprs/src/processes/mod.rs`, `crates/gammalooprs/src/processes/process.rs`.
- Amplitude and cross-section pipelines:
  - `crates/gammalooprs/src/processes/amplitude.rs`
  - `crates/gammalooprs/src/processes/cross_section.rs`
- Integrand abstraction and implementations:
  - `crates/gammalooprs/src/integrands/mod.rs`
  - `crates/gammalooprs/src/integrands/process/mod.rs`
  - `crates/gammalooprs/src/integrands/process/amplitude/mod.rs`
  - `crates/gammalooprs/src/integrands/process/cross_section_integrand.rs`
- Integration engine: `crates/gammalooprs/src/integrate/mod.rs`.
- Global/runtime settings: `crates/gammalooprs/src/settings/mod.rs`, `crates/gammalooprs/src/settings/global.rs`, `crates/gammalooprs/src/settings/runtime.rs`.
- Differential event/observable pipeline:
  - `crates/gammalooprs/src/observables/events.rs`
  - `crates/gammalooprs/src/observables/observables.rs`
  - `crates/gammalooprs/src/observables/clustering/*`
  - `crates/gammalooprs/src/integrands/evaluation.rs`

## Lifecycle and Data Flow

### 1. Startup
1. CLI parses command input (`OneShot::parse_env_with_capture`).
2. Initialization runs (`initialise()`), installs hooks and symbol registries.
3. State is loaded from state folder if available (`State::load`), otherwise default state is created.
4. CLI settings/runtime defaults are loaded and synchronized with tracing filters.

### 2. Process Generation Flow
1. `generate` command builds `ProcessDefinition` (from syntax or graph import).
2. `State::generate_integrand(s)` creates a generation thread pool.
3. `ProcessList::preprocess` delegates to amplitude/cross-section preprocessors.
4. `ProcessList::generate_integrands` builds `ProcessIntegrand` instances from preprocessed graphs.
5. Optional compile/export steps persist compiled evaluators and DOT/standalone outputs.

### 3. Evaluation and Integration Flow
1. Commands (`inspect`, `evaluate`, `integrate`) resolve process + integrand references.
2. Integrand is warmed up (`ProcessIntegrand::warm_up`) to initialize rotations and caches.
3. Sampling path parameterizes points and evaluates graph terms.
4. Stability checks may escalate precision (`f64 -> f128 -> arbitrary`) and rotate kinematics.
5. Cross-section graph evaluation returns a rich `GraphEvaluationResult<T>` rather
   than only a complex weight. This carries:
   - the graph contribution
   - grouped generated events
   - event-processing timing
   - generated / accepted event counts
6. Selectors are applied during event generation inside the cross-section graph
   evaluation path. Failed selector events are zeroed locally and do not survive
   into the final retained result.
7. Stability selection still compares only the complex graph weight, but the
   retained branch also carries the final grouped event payload.
8. The final `EvaluationResult` contains:
   - the stable `integrand_result`, before any parameterization Jacobian is
     applied
   - the top-level `parameterization_jacobian` when the sample came from
     x-space parameterization (`None` for direct momentum-space evaluation)
   - the separate `integrator_weight`, i.e. the Monte Carlo/grid weight only
   - grouped accepted events
   - event-processing timing in `evaluation_metadata`
   - generated / accepted event counts in `evaluation_metadata`
   - ordered per-level stability results, each with relative-accuracy and total
     time spent in that stability level
   - evaluation metadata
9. Observable filling happens only from the final stable `EvaluationResult`
   payload. Unstable branches do not contribute events or observables.
10. `havana_integrate` runs iterative Monte Carlo updates, merges worker-local
    observable accumulators, and writes integration artifacts.

### 3.1 Differential event-processing runtime

The differential pipeline is split into three layers:

- Declarative settings in `RuntimeSettings`:
  - `quantities`
  - `selectors`
  - `observables`
- A compiled runtime-only cache:
  - `EventProcessingRuntime`
  - owned by the process integrand
  - built in `warm_up()`
  - invalidated when runtime settings are mutated
- Per-sample / per-integration results:
  - `GraphEvaluationResult<T>`
  - `EvaluationResult`
  - observable accumulator state / snapshots

The runtime-only cache is intentionally not serialized in binary state dumps.
It is wrapped in a no-op `Encode` / default-empty `Decode` container and is
always rebuilt from settings at warm-up time.

### 3.2 Differential event model

Events are stored in precision-homogeneous containers:

- `GenericEvent<T>`
- `GenericEventGroup<T>`
- `GenericEventGroupList<T>`

This avoids enum-based mixed-precision layouts that would otherwise force
storage to be sized by the largest precision variant.

Each event carries:

- kinematics
- cut metadata
- one fully normalized event `weight`
- `additional_weights`

Observable-specific entry reweighting remains internal to the observable
runtime; it is not stored on the event.

`additional_weights` is a generic `BTreeMap` keyed by lightweight identifiers
such as:

- `FullMultiplicativeFactor`
- `Original`
- `ThresholdCounterterm { subset_index }`

This is populated only when
`settings.general.store_additional_weights_in_event = true`.

### 3.3 Event generation policy

Event generation is now governed by `settings.general.generate_events`:

- `false` (default):
  - events are generated only when selectors are present
  - selector-only events are discarded immediately after selection
  - selector+observable events may be generated temporarily to fill
    observables, but are cleared from returned sample results and integration
    event outputs afterwards
- `true`:
  - events are always generated
  - returned sample results keep the generated event groups even when no
    selector or observable is configured

This policy supersedes the earlier implicit "generate events whenever selectors
or observables are configured" rule.

Accepted LU events are retained as grouped event lists:

- one `EventGroup` per graph-group evaluation
- all accepted cuts from graphs sharing the same `group_id` are merged into the
  same retained event group
- each event stores its concrete `graph_id` and `cut_id`
- incoming PDGs are sourced from the initial-state cut edges

This preserves the correlation structure needed by downstream observables.

### 3.4 Observable architecture

Observable definitions are shared between selector and histogram use-cases:

- one quantity definition extracts zero or more `ObservableEntry<T>` values from
  an event
- selectors consume those entries on a single event
- histogram observables consume grouped events and fill per-group contributions

Histogram accumulation is based on:

- a histogram-level `sample_count`
- sparse per-bin sufficient statistics:
  - `entry_count`
  - `sum_weights`
  - `sum_weights_squared`
  - `mitigated_fill_count`
- `HistogramAccumulatorState` / `ObservableAccumulatorBundle` for mergeable
  runtime state
- `ObservableSnapshotBundle` for JSON / HwU output and API responses

Each Monte Carlo sample contributes exactly one statistical sample to a
histogram. Contributions from all retained event groups are summed first per
bin, and untouched bins are handled implicitly through the histogram-level
`sample_count` rather than explicit zero fills.

Histogram snapshots therefore remain fully mergeable and re-constructible into
live accumulator state. The Rust-facing histogram API exposes `merge(...)`,
`merge_in_place(...)`, and `rebin(...)`.

Current built-in quantities include:

- `particle` with `computation = scalar | count | pair`
- `jet` with `computation = scalar | count | pair`
- `afb`
- `cross_section`

Current pair quantities include `DeltaR`. Current scalar projections include
`E`, `CosTheta`, `PT`, `y`, `eta`, `Px`, `Py`, `Pz`, and `Mass`.

Object and pair quantities are ordered before `entry_selection` is applied.
The ordering contract is configurable per quantity through:

- `ordering = PT | Energy | AbsRapidity | Quantity`
- `order = Ascending | Descending`

Family defaults are:

- particle scalar quantities: `ordering = Quantity`
- jet scalar quantities: `ordering = PT`
- pair quantities: `ordering = Quantity`

`leading_only` and `nth_only` therefore operate on this explicit ordered list,
not on incidental source/event enumeration order.

Jet clustering is implemented natively in Rust with IRC-safe algorithms:

- `kt`
- `cambridge_aachen`
- `anti_kt`

and is validated against `fjcore` in tests.

## Persistence and Runtime Artifacts
Primary persisted state lives under `gammaloop_state/` (default):
- `state_manifest.toml` (state schema/version marker)
- `model.json`, `model_parameters.json`
- `symbolica_state.bin`
- `processes/` (amplitudes/cross_sections + integrands)
- `run.toml`
- `default_runtime_settings.toml`
- `global_settings.toml`
- `logs/`

For local experimentation, prefer an isolated path such as `.local/scratch/<run>/gammaloop_state` to keep repository root output minimal.

The persistence model is file-system based and intentionally human-editable for settings/run cards, mixed with binary artifacts for performance-heavy data.

### Persistence Compatibility Contract
- State format is versioned with `state_manifest.toml` (`version = 1` currently).
- `State::load` validates the manifest version and rejects states from newer binaries.
- If no manifest is present, load falls back to a legacy compatibility path (`version = 0`) and runs migration checks on required layout entries before loading.
- Process settings history now uses `settings_history.toml` consistently; loader still accepts legacy `settings_history.yaml` for backward compatibility and migration.

## Configuration Architecture
Configuration is split into:
- Global settings (`GlobalSettings`): generation, logging directives, parallelism.
- Runtime settings (`RuntimeSettings`): kinematics, integrator behavior, sampling, stability, subtraction, and differential event-processing configuration.

Differential runtime configuration currently includes:

- `general.generate_events`
- `general.store_additional_weights_in_event`
- named `quantities`
- named `selectors`
- named `observables`
- `integrator.observables_output`

This split is good: generation-time and evaluation-time concerns are separated and independently serializable.

## Concurrency Model
Concurrency is explicit and use-case scoped:
- Generation and compile thread pools use configurable thread counts.
- Integrator parallelism is controlled via runtime/global settings.
- Some loops over processes/integrands remain sequential at orchestration level while heavy operations inside are parallelized.

For the differential pipeline:

- each worker/integrand clone owns its own mutable observable accumulator state
- local worker results are merged back into the master integrand after chunk
  evaluation
- distributed batch mode can return either grouped events or serializable
  observable accumulator bundles
- observable snapshots are written only from the merged master state

## Testing Architecture
Testing is organized into layers:
- Unit tests in modules across `crates/gammalooprs/src/`.
- Integration test crate: `tests/` with reusable harness (`tests/src/lib.rs`).
- Additional Rust integration tests in `crates/gammalooprs/tests/`.
- Snapshot-heavy coverage for symbolic outputs and numerics.

Differential coverage is currently split into:

- `tests/tests/test_clustering.rs`
  - compares the native jet clustering implementation against `fjcore`
- `tests/tests/test_differential.rs`
  - validates grouped event propagation
  - selector behavior
  - observable filling and merge semantics
  - JSON / HwU snapshot output
  - API-facing result shaping using handcrafted fixtures
- `tests/tests/test_evaluation_api.rs`
  - real generated-process end-to-end coverage for Rust `evaluate_sample(s)`
  - graph/cut grouping and incoming-PDG propagation checks
  - selectors / observables / event-group retention
  - batch-global observable snapshots for `evaluate_samples`
  - per-sample evaluation metadata
  - x-space vs momentum-space Jacobian behavior
  - minimal-output mode
- `tests/tests/test_runs.rs`
  - real generated-process end-to-end coverage for `integrate`
  - integration-workspace observable output in JSON / HwU formats
  - repeated `save dot` overwrite behavior
- `tests/tests/test_python_api.rs`
  - subprocess-based end-to-end coverage for the public Python API
  - returned event grouping / ids / incoming PDGs
  - selectors / observables / returned metadata
  - minimal-output mode
  - x-space batch and momentum-space sample coverage
  - requires a separately prepared Python environment where
    `import gammaloop` already works



## Architectural Strengths
- Clear separation between application shell (`gammaloop-api`) and domain core (`gammalooprs`).
- Rich, serializable settings model with schema generation hooks.
- Strong process abstraction (`Process`, `ProcessCollection`, `ProcessIntegrand`) that supports amplitude and cross-section workflows.
- Robust stability pipeline with multi-precision escalation and rotation checks.
- Differential event-processing is now integrated into the core evaluation path
  without compromising the stable-weight selection logic.
- Observable accumulation is mergeable across rayon workers and distributed
  batch outputs.
