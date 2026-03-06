# gammaLoop Architecture

## Scope
This document describes the architecture of the current repository and proposes improvements that can reduce coupling, improve maintainability, and make runtime behavior easier to reason about.

The workspace is organized around two main Rust crates:
- `gammalooprs` (root crate): core physics/domain logic, graph processing, integrand construction, evaluation, and integration.
- `gammaloop-api` (workspace member): CLI, REPL, Python bindings, command parsing, and persisted state orchestration.

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
  -> Monte Carlo integration and stability checks

Infrastructure
  -> filesystem persistence (`gammaloop_state/`)
  -> Symbolica state import/export
  -> tracing/logging
```

## Main Components

### 1. Entry Points and Interfaces
- CLI entry point: `gammaloop-api/src/cli.rs`.
- Main app orchestration: `gammaloop-api/src/lib.rs` (`OneShot::load`, `OneShot::run`).
- Python module entry: `gammaloop-api/src/python.rs` (`#[pymodule] fn _gammaloop`).
- Python package shim: `gammaloop-api/python/gammaloop/__init__.py`.

### 2. Application State and Command Model
- Central mutable app state: `gammaloop-api/src/state.rs` (`State`).
- Command history and run cards: `CommandHistory`, `RunHistory`.
- Command dispatch: `gammaloop-api/src/commands/mod.rs` with concrete subcommands in `commands/*`.

The command model is stateful by design: commands mutate a long-lived `State` that can be saved and resumed.

### 3. Domain Core (gammalooprs)
- Root module wiring: `src/lib.rs`.
- Model and parameters: `src/model.rs`.
- Graph domain: `src/graph.rs` and submodules.
- Process orchestration: `src/processes/mod.rs`, `src/processes/process.rs`.
- Amplitude and cross-section pipelines:
  - `src/processes/amplitude.rs`
  - `src/processes/cross_section.rs`
- Integrand abstraction and implementations:
  - `src/integrands.rs`
  - `src/gammaloop_integrand/mod.rs`
  - `src/gammaloop_integrand/amplitude/mod.rs`
  - `src/gammaloop_integrand/cross_section_integrand.rs`
- Integration engine: `src/integrate.rs`.
- Global/runtime settings: `src/settings/mod.rs`, `src/settings/global.rs`, `src/settings/runtime.rs`.

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
4. `ProcessList::generate_integrands` builds `GLIntegrand` instances from preprocessed graphs.
5. Optional compile/export steps persist compiled evaluators and DOT/standalone outputs.

### 3. Evaluation and Integration Flow
1. Commands (`inspect`, `evaluate`, `integrate`) resolve process + integrand references.
2. Integrand is warmed up (`GLIntegrand::warm_up`) to initialize rotations and caches.
3. Sampling path parameterizes points and evaluates graph terms.
4. Stability checks may escalate precision (`f64 -> f128 -> arbitrary`) and rotate kinematics.
5. `havana_integrate` runs iterative Monte Carlo updates and writes integration artifacts.

## Persistence and Runtime Artifacts
Primary persisted state lives under `gammaloop_state/` (default):
- `model.json`, `model_parameters.json`
- `symbolica_state.bin`
- `processes/` (amplitudes/cross_sections + integrands)
- `run.toml`
- `default_runtime_settings.toml`
- `cli_settings.toml`
- `logs/`

The persistence model is file-system based and intentionally human-editable for settings/run cards, mixed with binary artifacts for performance-heavy data.

## Configuration Architecture
Configuration is split into:
- Global settings (`GlobalSettings`): generation, logging directives, parallelism.
- Runtime settings (`RuntimeSettings`): kinematics, integrator behavior, sampling, stability, subtraction.

This split is good: generation-time and evaluation-time concerns are separated and independently serializable.

## Concurrency Model
Concurrency is explicit and use-case scoped:
- Generation and compile thread pools use configurable thread counts.
- Integrator parallelism is controlled via runtime/global settings.
- Some loops over processes/integrands remain sequential at orchestration level while heavy operations inside are parallelized.

## Testing Architecture
Testing is organized into layers:
- Unit tests in modules across `src/`.
- Integration test crate: `integration-tests/` with reusable harness (`integration-tests/src/lib.rs`).
- Additional Rust integration tests in `tests/`.
- Snapshot-heavy coverage for symbolic outputs and numerics.

## Architectural Strengths
- Clear separation between application shell (`gammaloop-api`) and domain core (`gammalooprs`).
- Rich, serializable settings model with schema generation hooks.
- Strong process abstraction (`Process`, `ProcessCollection`, `GLIntegrand`) that supports amplitude and cross-section workflows.
- Robust stability pipeline with multi-precision escalation and rotation checks.

## Proposed Improvements

### P0: Fix persistence consistency and add state versioning
- Problem: process settings history is loaded from `settings_history.yaml` but saved as `settings_history.toml` in `src/processes/process.rs`.
- Improvement: enforce a single format and add an explicit state manifest version (`state_manifest.toml`), then gate load with migration checks.
- Impact: prevents silent partial loads and future format drift.

### P1: Introduce a dedicated application service layer
- Problem: clap command types currently execute domain logic directly and mutate `State` in-place.
- Improvement: add a service boundary (for example `AppService`) that owns use-cases (`generate`, `integrate`, `inspect`, `save`). Keep commands as parsing/DTO only.
- Impact: reduces coupling, improves testability, and makes Python/CLI behavior easier to keep equivalent.

### P1: Replace manual cache-ID protocol with typed cache sessions
- Problem: `GammaloopIntegrand` cache semantics rely on manual counter updates and conventions spread across many methods.
- Improvement: introduce a `CacheSession` / `ExternalConfigEpoch` abstraction that centralizes increment/revert semantics and can be validated at compile-time boundaries.
- Impact: fewer cache-invalidations bugs and cleaner evaluation code.

### P1: Unify amplitude/cross-section integrand pipeline shape
- Problem: both paths implement similar lifecycle steps (`preprocess`, `build`, `warm_up`, `compile`, `save/load`) with duplicate orchestration patterns.
- Improvement: extract shared pipeline traits/utilities for process and integrand lifecycle while retaining domain-specific internals.
- Impact: lower maintenance cost and more predictable feature rollout across both generation types.

### P2: Harden unsupported and partially-implemented runtime paths
- Problem: several user-facing paths still rely on `todo!()` or non-implemented branches (e.g. batch command, some event-processing and cross-section standalone behavior).
- Improvement: convert unsupported paths to explicit typed errors and feature flags; track with test coverage for error contracts.
- Impact: avoids panic-style failures in production workflows.

### P2: Remove panic-oriented error paths in sampling/inspect
- Problem: some paths still panic when settings are invalid instead of returning structured errors.
- Improvement: propagate `Result` consistently through inspect/integrand parameterization and command handlers.
- Impact: safer automation and clearer diagnostics in REPL/run cards.

### P3: Process-level orchestration parallelism
- Problem: outer loops over processes/integrands are mostly sequential even when independent.
- Improvement: parallelize process-level preprocessing/build/compile where deterministic ordering is not required.
- Impact: better throughput for multi-process workloads.

### P3: Improve observability for long-running pipelines
- Problem: tracing is already strong, but instrumentation is uneven across some high-cost transitions.
- Improvement: define a small set of canonical spans and metrics for `load -> preprocess -> build -> warm_up -> evaluate/integrate` and emit stable fields.
- Impact: easier profiling and regression triage.

## Suggested Implementation Order
1. Fix persistence format mismatch and introduce versioned manifest.
2. Extract `AppService` use-case boundary from command handlers.
3. Refactor cache handling into typed session semantics.
4. Consolidate amplitude/cross-section lifecycle scaffolding.
5. Replace remaining `todo!()`/panic paths in user-facing flows with typed errors.
