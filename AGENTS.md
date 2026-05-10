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
### Rust
- Always format and run cargo clippy warnings before finishing
- Naming: `snake_case` for functions/modules, `CamelCase` for types/traits.
- Prefer methods on types to bare functions
- Follow the clippy suggestions religiously
- Aggressively use pub(crate) except when absolutely necessary to have fully public exposure (helps monitor unused components with clippy)

### Python
- use the ruff formatter.

## Testing Guidelines
- Rust integration tests live in `tests/` (files like `test_runs.rs`).
- Use `cargo nextest` via `just test TEST_NAME` for targeted runs.

## Commit & Pull Request Guidelines
- Commit messages are short and descriptive, typically lowercase without scopes (e.g., "remove edge quotes"), but contain the full user prompt in the description
- If there is already a commit, prefer making a new one
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
- Always `cargo fmt` and `cargo check` before compiling to catch easy errors and ensure consistent formatting.
- API stability is currently a non-goal for internal development: prioritize maintainability and structure over preserving external-facing module paths.
- Breaking Rust/Python API changes are acceptable when they simplify architecture; document major moves in the change description.
- Backward compatibility with old on-disk GammaLoop states is also currently a non-goal: it is acceptable to rename/remove legacy state files, layouts, and loaders when that simplifies the implementation.
- Do not keep compatibility fallbacks for old state formats/names unless the task explicitly asks for migration support.
- Never be afraid to modify existing functions to obtain the most elegant and concise code as opposed to trying to be backward compatible or leave existing code untouched.
- When suitable, prefer implementing behavior as methods on existing structs/enums instead of free functions to improve discoverability.
- In generic or arbitrary-precision code, do not introduce constants or intermediate values through `from_f64`, `std::f64::consts::*`, or other lossy `f64` routes unless that exact location is an explicit `f64` boundary by design (for example persisted settings that are already `f64`, or histogram/output accumulation that is intentionally `f64`).
- `f64` values originating from user-supplied settings are an allowed boundary: converting those setting values into `F<T>` with `from_f64` is acceptable. Do not extend that exception to internally generated constants or intermediate values.
- When working with `F<T>` or other precision-generic numeric code, build constants from an in-scope representative value of the correct type using helpers such as `.zero()`, `.one()`, `.epsilon()`, `.PI()`, `.TAU()`, `.from_usize()`, `.from_isize()`, etc., so the active precision is preserved exactly.
- If a needed constant or operation cannot be expressed without a dubious precision-losing conversion, stop and ask instead of guessing.
- Do not explicitly bring `symbolica::domains::float::FloatLike` into scope in GammaLoop code. Prefer the GammaLoop wrappers in `crate::utils`, and if `F<T>` is missing an operation needed for precision-safe generic code, add or use a helper there instead. If it genuinely looks unavoidable to use Symbolica's trait directly, stop and ask first.
- The settings serialization model relies on the `SHOWDEFAULTS` escape hatch in `gammalooprs::utils::serde_utils` when writing persisted state defaults and when building completion catalogs from serialized settings.
- When adding or changing settings fields, verify that `Default::default()` and the corresponding `skip_serializing_if` helper stay aligned; custom helpers must route through `show_defaults_helper(...)` (or `IsDefault::is_default`) or the field will disappear from saved `global_settings.toml` / `default_runtime_settings.toml` and from completion.
- Treat saved-state detection as manifest-based only. Empty folders or folders containing only transient runtime artifacts such as `logs/` are scratch state, not legacy state, and must be treated as blank state.
- Do not introduce writes into the state folder outside explicit `save state` / `quit -o`, except for logfile tracing when that logger is actually enabled.
- Honor `--read-only-state` consistently. When it is enabled, do not write into the state folder and prefer cwd-based fallbacks for transient artifacts that would otherwise default into the state.
- Ignore unrelated untracked local artifacts by default (editor swap files, scratch docs, local example edits, profiling outputs, etc.) unless the task clearly requires them; do not stop work solely because such untracked files are present.

## Design / Refactor Preferences
- Prefer one semantic boundary per concept. Before adding a trait, helper, mode, or adapter, look for the existing abstraction that already owns the behavior and extend that instead.
- Collapse duplicate abstractions aggressively. If two names describe the same responsibility, merge them and update call sites rather than keeping compatibility shims.
- Prefer methods, trait methods, or impl-associated functions over bare free functions. A helper should live on the type or trait whose invariant it uses.
- Do not preserve internal API compatibility unless explicitly requested. Breaking internal Rust APIs is acceptable when it removes duplication or clarifies ownership.
- Avoid compatibility fallbacks and "old path plus new path" designs. Pick the clearer model and migrate usages.

## Implementation Style
- Keep changes minimal but complete: remove obsolete code paths, imports, docs, and call sites in the same change.
- Prefer concrete return types over sentinel-style APIs. Avoid `Option` where `None` means "nothing happened" if the caller always needs a usable result.
- Make transformations explicit and idempotent. If recursive rewriting is involved, include a fixed-point check or another clear loop guard.
- Avoid relying on fragile errors for control flow when a direct predicate or structural check would express the intent better.
- Do not add broad helper layers just to avoid touching call sites. Update the call sites when the API should change.

## Investigation Discipline
- Before fixing a suspected bug, find the existing code path that should own the behavior.
- When asked to identify failure modes, do not patch yet; report the concrete failing cases and where they enter the code.
- When asked for a plan, do not make code changes until the plan is accepted.
- If a requested change looks like only a rename, stop and check whether there is a deeper design change implied.

## Tests And Docs
- Put targeted unit tests next to the implementation they exercise. Use broader integration tests only for cross-module behavior.
- Add tests for edge cases that motivated the change, especially when collapsing duplicated logic.
- Documentation should explain conventions and algorithmic intent, not restate function names. Avoid tautological docs like "apply the convention"; say what the convention does.
- Keep docs synchronized with API refactors; remove references to deleted traits, functions, and modes.

## Differential LU / Event Processing Notes
- `differential_lu.md` is the detailed implementation log for the current differential LU stack; `docs/architecture/architecture-current.md` has the corresponding implemented-architecture summary. Keep both in sync when changing selectors, observables, event grouping, or sample-evaluation output.
- Event grouping semantics are by graph-group, not just by graph: if multiple LU graphs share the same `group_id`, all accepted cuts from all of those graphs belong in the same retained `EventGroup`.
- Histogram error propagation is based on Monte Carlo samples, not on individual observable entries. Treat the full retained event-group list from one evaluation as one statistical sample for each histogram.
- For histogram accumulation, do not explicitly zero-fill untouched bins. The current design is sparse:
  - histogram-level `sample_count`
  - per-bin raw stats `entry_count`, `sum_weights`, `sum_weights_squared`, `mitigated_fill_count`
  - underflow and overflow are full bins with the same raw-stat structure
- Histogram snapshots are intentionally raw-stat snapshots, not presentation-only views. They must remain mergeable and reconstructible back into live histogram state. Rust-side helpers already exist for `merge(...)`, `merge_in_place(...)`, `rebin(...)`, and reconstruction into accumulator state.
- `general.generate_events` controls whether events are retained in returned results; observables and selectors may still force temporary event generation internally when this is `false`.
- `general.store_additional_weights_in_event` controls whether `additional_weights` are stored on events. The map is a `BTreeMap` keyed by lightweight identifiers such as `Original`, `ThresholdCounterterm { subset_index }`, and `FullMultiplicativeFactor`.

## API / Python Interop Notes
- The Rust API has precise endpoints `evaluate_sample_precise` and `evaluate_samples_precise`; the Python API intentionally stays `f64`-only.
- Python API end-to-end tests are subprocess-based and assume the user already built/installed the Python extension in the active environment (`maturin develop` or `just build-api`). Do not try to embed Python into the Rust test binary for these tests.
- The top-level `.venv` in the repo is a reasonable default development environment for running the Python API examples/tests after building the extension there.
- If the extension module environment is an issue when compiling the Python api, your can try `build-api-abi`.

## Environment and Build Quirks
- `NO_SYMBOLICA_OEM_LICENSE` is read via `option_env!` in `gammalooprs` initialization, so it is a compile-time switch, not a pure runtime one. If you need to disable the OEM-license activation path, that variable must be present when building the relevant binary or Python extension.

## CLI / State Notes
- The run-card field is `command_blocks` (singular `command`), not `commands_blocks`. No compatibility alias is kept.
- Integration workspace resume must restore the exact state at the end of the last completed iteration. The authoritative resume data lives in the workspace manifest/internal state plus per-slot settings and internal observable checkpoints, not in user-facing `integration_result.json`.
- Treat workspace snapshots as completed-iteration checkpoints only: do not persist partial-iteration state, and keep user-facing latest snapshots plus optional `_iter_0001` history separate from the internal resume state.

## Configuration & State Tips
- CLI runs create a `gammaloop_state/` directory by default; keep it out of commits unless intentionally sharing a reproducible state.
- Use `./gammaloop -s <state_dir>` for isolated experiments.
