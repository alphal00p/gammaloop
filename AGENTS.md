# Repository Guidelines

## Project Structure & Module Organization
- `crates/gammalooprs/src/` holds the core Rust implementation.
- `crates/gammaloop-api/` is the workspace member for the Rust API + Python bindings (`crates/gammaloop-api/src` and `crates/gammaloop-api/python/gammaloop`).
- `tests/` contains the integration-test crate and shared fixtures in `tests/resources/`.
- `examples/` includes command cards and notebooks for end-to-end runs.
- `bin/` provides scripts like `compile.sh`, `run_tests.sh`, and the `gammaloop` entrypoints after build.
- `assets/` stores schemas/data and model files (`assets/models/`); Rust benches live under `crates/gammalooprs/benches/`.
- Architecture docs live in `docs/architecture/`:
  - `docs/architecture/architecture-current.md` for implemented architecture.
  - `docs/architecture/architecture-ideas.md` for roadmap/proposals.

## Build, Test, and Development Commands
- Build from source (full deps + binaries):
  ```bash
  ./bin/compile.sh
  ```
- Build Rust CLI in dev-optim:
  ```bash
  just build-cli
  ```
- Build Rust CLI with stable Python ABI:
  ```bash
  just build-cli-abi
  ```
- Build Python API bindings:
  ```bash
  just build-api
  ```
- Build Python API bindings with stable ABI:
  ```bash
  just build-api-abi
  ```
- Format and lint:
  ```bash
  just fmt
  just clippy
  ```
- If running into a linking issue complaining about missing `__emul...` symbol, then it may be a macos specific issue, that can be fixed by 
  building with the following environment variable change `EXTRA_MACOS_LIBS_FOR_GNU_GCC=T` (see impact of this in file `build.rs`.

## Coding Style & Naming Conventions
- Rust: use `rustfmt` formatting and address `clippy` warnings before merging.
- Naming: `snake_case` for functions/modules, `CamelCase` for types/traits.
- Python (bindings and helpers): 4-space indentation; `black` is listed in `pyproject.toml` dev deps.
- Aggressively use pub(crate) except when absolutely necessary to have fully public exposure (helps monitor unused components with clippy)

## Testing Guidelines
- Rust integration tests live in `tests/` (files like `test_runs.rs`).
- Use `cargo nextest` via `just test TEST_NAME` for targeted runs, or `./bin/run_tests.sh rust` for a full pass.
- Python tests use `pytest` (see `./bin/run_tests.sh python` or `python -m pytest` inside the module).

## Commit & Pull Request Guidelines
- Commit messages are short and descriptive, typically lowercase without scopes (e.g., "remove edge quotes").
- PRs should include: a clear summary, test command(s) run, and any example command card used to validate behavior.
- If changes touch CLI state, mention `gammaloop_state/` impacts or migration steps.

## Agent Workflow (jj)
- If the repo is colocated with jj metadata, prefer `jj` over `git`.
- Start new tasks with `jj new` and immediately `jj describe -m "<plan>"`.
- Use `jj log --no-pager` and `jj diff --no-pager` to avoid paging.
- Check the current change before editing; if the working copy has changes but no description, run `jj diff` then `jj describe -m "<description>"`.
- Update the current change description with `jj describe -m "<description>"` as the plan evolves.
- Use `jj diff` to see changed files.
- Reference: https://jj-vcs.github.io/jj/prerelease/cli-reference/

## Python API Features
- `python_api`: Basic Python API support (default for `gammaloop-api`)
- `python_abi`: Enables stable Python ABI (abi3) for cross-version compatibility
  - Use for distribution builds that need to work across Python versions
  - Adds some performance overhead compared to version-specific builds
  - Not enabled by default - opt in with `--features python_abi`

## Agent Instructions
- API stability is currently a non-goal for internal development: prioritize maintainability and structure over preserving external-facing module paths.
- Breaking Rust/Python API changes are acceptable when they simplify architecture; document major moves in the change description.
- Backward compatibility with old on-disk GammaLoop states is also currently a non-goal: it is acceptable to rename/remove legacy state files, layouts, and loaders when that simplifies the implementation.
- Do not keep compatibility fallbacks for old state formats/names unless the task explicitly asks for migration support.
- The settings serialization model relies on the `SHOWDEFAULTS` escape hatch in `gammalooprs::utils::serde_utils` when writing persisted state defaults and when building completion catalogs from serialized settings.
- When adding or changing settings fields, verify that `Default::default()` and the corresponding `skip_serializing_if` helper stay aligned; custom helpers must route through `show_defaults_helper(...)` (or `IsDefault::is_default`) or the field will disappear from saved `global_settings.toml` / `default_runtime_settings.toml` and from completion.

## Configuration & State Tips
- CLI runs create a `gammaloop_state/` directory by default; keep it out of commits unless intentionally sharing a reproducible state.
- Use `./bin/gammaloop -s <state_dir>` for isolated experiments.
