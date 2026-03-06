# gammaLoop Current Architecture

## Scope
This document describes the current, implemented architecture of this repository.

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
  -> optional local scratch persistence (`.local/scratch/...`)
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
- Model and parameters: `src/model/mod.rs`.
- Graph domain: `src/graph/mod.rs` and submodules.
- Process orchestration: `src/processes/mod.rs`, `src/processes/process.rs`.
- Amplitude and cross-section pipelines:
  - `src/processes/amplitude.rs`
  - `src/processes/cross_section.rs`
- Integrand abstraction and implementations:
  - `src/integrands/mod.rs`
  - `src/integrands/process/mod.rs`
  - `src/integrands/process/amplitude/mod.rs`
  - `src/integrands/process/cross_section_integrand.rs`
- Integration engine: `src/integrate/mod.rs`.
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
4. `ProcessList::generate_integrands` builds `ProcessIntegrand` instances from preprocessed graphs.
5. Optional compile/export steps persist compiled evaluators and DOT/standalone outputs.

### 3. Evaluation and Integration Flow
1. Commands (`inspect`, `evaluate`, `integrate`) resolve process + integrand references.
2. Integrand is warmed up (`ProcessIntegrand::warm_up`) to initialize rotations and caches.
3. Sampling path parameterizes points and evaluates graph terms.
4. Stability checks may escalate precision (`f64 -> f128 -> arbitrary`) and rotate kinematics.
5. `havana_integrate` runs iterative Monte Carlo updates and writes integration artifacts.

## Persistence and Runtime Artifacts
Primary persisted state lives under `gammaloop_state/` (default):
- `state_manifest.toml` (state schema/version marker)
- `model.json`, `model_parameters.json`
- `symbolica_state.bin`
- `processes/` (amplitudes/cross_sections + integrands)
- `run.toml`
- `default_runtime_settings.toml`
- `cli_settings.toml`
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
- Strong process abstraction (`Process`, `ProcessCollection`, `ProcessIntegrand`) that supports amplitude and cross-section workflows.
- Robust stability pipeline with multi-precision escalation and rotation checks.
