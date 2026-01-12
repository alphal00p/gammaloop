# Repository Guidelines

## Project Structure & Module Organization
- `src/` holds the core Rust implementation and CLI support code.
- `gammaloop-api/` is the workspace member for the Rust API + Python bindings (`gammaloop-api/src` and `gammaloop-api/python/gammaloop`).
- `tests/` contains Rust integration tests and shared fixtures in `tests/resources/`.
- `examples/` includes command cards and notebooks for end-to-end runs.
- `bin/` provides scripts like `compile.sh`, `run_tests.sh`, and the `gammaloop` entrypoints after build.
- `assets/`, `models/`, and `benches/` store schemas/data, model files, and benchmarks.

## Build, Test, and Development Commands
- Build from source (full deps + binaries):
  ```bash
  ./bin/compile.sh
  ```
- Build Rust CLI in dev-optim:
  ```bash
  just build-cli
  ```
- Build Python API bindings:
  ```bash
  just build-api
  ```
- Run an example workflow:
  ```bash
  ./bin/gammaloop examples/command_cards/simple_workflow.toml
  ```
- Run tests:
  ```bash
  ./bin/run_tests.sh rust
  ./bin/run_tests.sh python
  ```
- Format and lint:
  ```bash
  just fmt
  just clippy
  ```

## Coding Style & Naming Conventions
- Rust: use `rustfmt` formatting and address `clippy` warnings before merging.
- Naming: `snake_case` for functions/modules, `CamelCase` for types/traits.
- Python (bindings and helpers): 4-space indentation; `black` is listed in `pyproject.toml` dev deps.

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

## Agent Instructions

## Configuration & State Tips
- CLI runs create a `gammaloop_state/` directory by default; keep it out of commits unless intentionally sharing a reproducible state.
- Use `./bin/gammaloop -s <state_dir>` for isolated experiments.
