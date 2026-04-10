# symjit branch plan

## Scope

Introduce a compiled-backend selector that replaces `inline_asm`, make
`symjit` the default compiled backend when `global.generation.evaluator.compile`
is enabled, keep the on-disk state portable by persisting eager Symbolica
 evaluators, and rebuild runtime-only backends after generation/load.

This file is the working plan for the branch and should stay current as the
implementation evolves.

## Locked product decisions

1. `global.generation.evaluator.compile` remains the eager-vs-compiled switch.
2. `global.generation.compile` gets a new backend selector with exactly:
   - `c++`
   - `assembly`
   - `symjit`
3. `inline_asm` is removed entirely from the active settings model.
4. `symjit` is the default compiled backend.
5. Existing integrands are locked to generation-time backend metadata embedded
   in `integrand.bin`, not to current global settings.
6. Current globals only affect future generation, except for one narrow load
   fallback policy: if an external compiled artifact cannot be loaded, startup
   globals may authorize a fallback to symjit.
7. That fallback is allowed only when startup globals still say compiled mode
   with backend `symjit`.
8. For every integrand that actually falls back from external compiled loading
   to symjit, GammaLoop must emit one `warn!`.
9. `display integrands` detail metadata must show:
   - the generation-time compilation mode/options stored with the integrand
   - the currently active f64 backend flavor (`eager`, `c++`, `assembly`, or
     `symjit`)
10. No backward compatibility is required for old state or old settings keys.

## Key findings from the final review

### 1. Backend choice is still split today

- `crates/gammalooprs/src/processes/mod.rs`
  - `EvaluatorSettings.compile` is the eager-vs-compiled gate.
- `crates/gammalooprs/src/settings/global.rs`
  - `GammaloopCompileOptions.inline_asm` currently distinguishes ASM vs plain
    C++ export.
- `crates/gammaloop-api/src/commands/generate.rs`
  - external compilation is only triggered when
    `global_settings.generation.evaluator.compile` is true.

So the new selector belongs under `global.generation.compile`, but the resolved
frozen backend still has to come from both settings together.

### 2. The Symbolica bump is mandatory and gives a cleaner external-artifact type

At `aa51532f4b6be0aefec8658bbda9719333d88f6f`, Symbolica exposes both:
- `ExpressionEvaluator::jit_compile(...)`
- `CompiledCode<T>`

`CompiledCode<T>` is important because it already serializes only `(path,
function_name)` and loads lazily via `.load()`. That is a better persisted type
for GammaLoop than the currently eager-loading `CompiledComplexEvaluator` /
`CompiledComplexEvaluatorSpenso`.

### 3. Per-process settings history is still not precise enough

`settings_history.toml` in `crates/gammalooprs/src/processes/process.rs` is not
usable as authoritative per-integrand compile metadata because:
- one process can contain multiple integrands
- targeted regeneration can rebuild only one integrand
- that path does not keep per-integrand generation metadata synchronized

So all frozen compile metadata needed for restoration/display must live inside
`AmplitudeIntegrandData` / `CrossSectionIntegrandData`.

### 4. The frozen metadata cannot be just a bare backend enum anymore

The latest requirement on `display integrands` changes this materially.
Showing the generation-time compilation mode/options accurately after load
cannot rely on current globals or process-level settings history.

Therefore the embedded per-integrand metadata should contain only fields that
actually froze the generated evaluator state, but that set is larger than just
`Eager | Cpp | Assembly | Symjit`.

Current recommendation:
- store a resolved per-integrand generation backend descriptor like:
  - `Eager`
  - `Symjit`
  - `Cpp { optimization_level, fast_math, unsafe_math, compiler, custom }`
  - `Assembly { optimization_level, fast_math, unsafe_math, compiler, custom }`

Why this is consistent with the earlier constraint:
- those external compile options do affect the generated compiled artifacts
- they are therefore part of the frozen generation state for externally
  compiled integrands
- settings that do not freeze the generated integrand should still stay out of
  this embedded metadata

### 5. Runtime backend activation needs to be integrand-wide, not evaluator-local

`GenericEvaluator` is used throughout:
- main evaluator stacks
- LU counterterms
- amplitude counterterms
- overlap helpers
- cross-section helper evaluators

If backend activation is handled independently per `GenericEvaluator`, a partial
external-load failure could leave one integrand in a mixed runtime state.
That would conflict with the requested single active backend flavor shown by
`display integrands`.

So the runtime policy should be:
- frozen metadata is stored per integrand
- runtime activation is performed across all `GenericEvaluator`s owned by that
  integrand
- the integrand ends up with one active f64 flavor:
  - eager
  - c++
  - assembly
  - symjit
- if any external artifact load fails for a frozen `Cpp` / `Assembly`
  integrand and symjit fallback is authorized, the whole integrand switches to
  symjit

### 6. Two runtime-only caches fit the existing architecture cleanly

`RuntimeCache<T>` in `crates/gammalooprs/src/integrands/process/mod.rs` is
already the established pattern for non-serialized runtime rebuilds.
That pattern fits both:
- symjit evaluators
- loaded external `.so` evaluators

So the clean design is:
- persist eager/rational evaluators plus frozen generation backend metadata
- persist external compiled artifact descriptors, not loaded libraries
- keep loaded `.so` evaluators runtime-only
- keep symjit evaluators runtime-only

### 7. The load path still supports a narrow post-deserialize initialization pass

`crates/gammaloop-api/src/lib.rs` already knows the startup CLI/global settings
when it calls `State::load(...)`.
So the lightest architecture remains:
- `State::load(...)` does structural deserialization
- a post-load pass initializes f64 runtime backends for loaded integrands
- that pass uses embedded frozen metadata as the primary selector
- startup globals are consulted only for the opt-in symjit fallback gate

### 8. `display integrands` currently has no backend metadata path

Current detail rendering is split across:
- `crates/gammaloop-api/src/integrand_info.rs`
- `crates/gammaloop-api/src/commands/display.rs`
- `crates/gammaloop-api/src/python.rs`

The summary table currently only shows:
- kind
- graphs
- graph groups
- `.bin size`

So the display change needs explicit new fields on `IntegrandInfo`, plus the
Python mirror and schema surface if kept in sync.

## Final implementation plan

### Step 1: Bump Symbolica/graphica/numerica

Files:
- `Cargo.toml`
- `Cargo.lock`

Work:
- bump all three patched revs to
  `aa51532f4b6be0aefec8658bbda9719333d88f6f`
- validate the new Symbolica API surface used here:
  - `jit_compile`
  - `CompiledCode<T>`

### Step 2: Replace `inline_asm` with a compiled-backend enum

Primary files:
- `crates/gammalooprs/src/settings/global.rs`
- `crates/gammaloop-api/src/commands/set.rs`
- `assets/schemas/`

Work:
- add a new enum, likely `CompilationMode`, under
  `GammaloopCompileOptions`
- remove `inline_asm`
- make `symjit` the default backend
- keep `evaluator.compile` unchanged as the eager-vs-compiled switch
- fix `to_symbolica_compile_options()` and actually route external compilation
  through it
- keep `export_settings()` as the place that maps:
  - `c++` -> `InlineASM::None`
  - `assembly` -> default/native inline asm
  - `symjit` -> no external export/compile path
- update schema/tests and non-`OLD` TOML fixtures

### Step 3: Add embedded frozen generation backend metadata to integrands

Primary files:
- `crates/gammalooprs/src/integrands/process/amplitude/mod.rs`
- `crates/gammalooprs/src/integrands/process/cross_section/mod.rs`
- `crates/gammalooprs/src/processes/amplitude.rs`
- `crates/gammalooprs/src/processes/cross_section.rs`

Work:
- add a per-integrand embedded metadata type stored inside `integrand.data`
- keep it limited to fields that actually freeze the generated evaluator state
- set it during integrand construction from generation settings
- ensure `ProcessIntegrand::resume_fingerprint()` naturally includes it through
  `integrand.data`

Current recommended shape:
- `Eager`
- `Symjit`
- `Cpp { external_compile_options... }`
- `Assembly { external_compile_options... }`

### Step 4: Split persisted external artifacts from runtime-loaded evaluators

Primary file:
- `crates/gammalooprs/src/integrands/process/evaluators.rs`

Work:
- replace serialized `f64_compiled: Option<CompiledComplexEvaluatorSpenso>`
  with:
  - a serialized `CompiledCode<Complex<f64>>`-style artifact descriptor
  - a runtime-only loaded evaluator cache
- add a runtime-only symjit cache
- keep eager/rational/f128/arb evaluators persisted as before

This is required so load-time fallback can happen after decode instead of being
forced by eager `.so` loading inside bincode decode.

### Step 5: Make backend initialization explicit and uniform

Primary files:
- `crates/gammalooprs/src/integrands/process/evaluators.rs`
- the integrand/container types that own `GenericEvaluator`s

Work:
- add explicit runtime backend state so f64 dispatch is not based on optional
  field presence alone
- provide initialization helpers for:
  - eager
  - external compiled
  - symjit
- traverse every `GenericEvaluator` under an integrand uniformly
- keep `f128` / `arb` eager

Important rule:
- the runtime backend actually active for an integrand should be a single
  integrand-wide flavor, not a mixed per-evaluator accident

### Step 6: Keep generation behavior backend-aware

Primary touchpoints:
- `crates/gammaloop-api/src/commands/generate.rs`
- `crates/gammaloop-api/src/state.rs`
- `crates/gammalooprs/src/processes/*`
- `crates/gammalooprs/src/integrands/process/*`

Work:
- resolve the frozen backend during integrand generation from:
  - `evaluator.compile`
  - `global.generation.compile.compilation_mode`
- initialize runtime backend immediately after generation:
  - eager -> eager
  - symjit -> build symjit now
  - c++ / assembly -> keep eager until external compile step finishes, then load
    external artifacts into runtime cache
- only call `compile_integrands(...)` for external backends

### Step 7: Add a post-load backend activation pass

Primary touchpoints:
- `crates/gammaloop-api/src/lib.rs`
- `crates/gammaloop-api/src/state.rs`
- `crates/gammalooprs/src/processes/*`
- `crates/gammalooprs/src/integrands/process/*`

Work:
- keep `State::load(...)` focused on deserialization
- add a narrow post-load activation pass over loaded integrands
- use embedded frozen metadata as the primary selector
- use startup globals only for the opt-in fallback gate

Post-load runtime behavior:
- frozen `Eager` -> stay eager
- frozen `Symjit` -> rebuild symjit from eager/rational data
- frozen `Cpp` / `Assembly` -> try loading all external artifacts
- if any external artifact load fails and startup globals authorize symjit:
  - rebuild symjit across that integrand
  - emit one `warn!` naming the integrand and the reason for fallback
- otherwise fail load

### Step 8: Update `display integrands` metadata

Primary files:
- `crates/gammaloop-api/src/integrand_info.rs`
- `crates/gammaloop-api/src/commands/display.rs`
- `crates/gammaloop-api/src/python.rs`

Work:
- extend `IntegrandInfo` with generation-time compile metadata derived from the
  embedded per-integrand frozen state
- extend it with the currently active f64 backend flavor
- render those in the detail summary table shown by `display integrands -p ...
  -i ...`

Current recommendation for displayed rows:
- `compile enabled`
- `generation backend`
- `generation compile options` (only when external backend)
- `active f64 backend`
- existing size/count rows stay unchanged

### Step 9: Update fixtures, schemas, docs, and tests

Likely touchpoints:
- `tests/resources/run_cards/`
- `examples/cli/`
- `assets/schemas/`
- `docs/architecture/architecture-current.md`
- tests around settings merging / display / save-load

Coverage to add:
- settings serde/default tests for the new backend selector
- compile option plumbing test so external compile no longer silently uses
  `CompileOptions::default()`
- roundtrip test that runtime-only backend caches are skipped in serialization
- load test for frozen `Symjit`
- load test for frozen external backend success
- load test for frozen external backend failure with startup-symjit fallback
- display-info test covering generation metadata + active backend output

## Main risks to control

1. The dependency bump may expose unrelated compile fallout in Symbolica-facing
   code.
2. The frozen metadata must stay minimal enough to avoid dragging mutable global
   settings into `integrand.bin`, but rich enough to make display and restore
   semantics accurate.
3. The external-artifact split must not accidentally change normal eager /
   compiled evaluation semantics.
4. Integrand-wide backend activation is important; a mixed runtime backend state
   would make the display contract ambiguous and complicate debugging.
5. The startup fallback gate must stay narrow:
   `global.generation.evaluator.compile = true` and compiled backend selector =
   `symjit`.

## Current status

- [x] Repository, architecture, and load/save review
- [x] Confirmed Symbolica rev/API needed for symjit
- [x] Identified settings, evaluator, load, and display touchpoints
- [x] Refined the frozen-metadata design to satisfy both restore semantics and
      display requirements
- [x] Final plan review with user
- [x] Implement settings/backend/runtime-cache changes
- [x] Update schemas/docs/examples
- [x] Final validation sweep (`cargo fmt`, `cargo check`, `cargo clippy`,
      targeted tests)

## Implementation status

Implemented:
- bumped the patched `symbolica`, `graphica`, and `numerica` revisions to
  `aa51532f4b6be0aefec8658bbda9719333d88f6f`
- replaced `global.generation.compile.inline_asm` with
  `global.generation.compile.compilation_mode`
- kept `global.generation.evaluator.compile` as the eager-vs-compiled switch
- made `symjit` the default compiled backend
- embedded frozen per-integrand compilation metadata directly into
  `integrand.bin`
- switched persisted external compiled evaluators to lazy
  `CompiledCode<Complex<f64>>` artifacts with runtime-only loaded-library caches
- added runtime-only symjit caches and explicit active-f64-backend tracking
- restored integrand runtime backends after load from frozen metadata
- added the narrow startup-driven symjit fallback for failed external `.so`
  loads and log one `warn!` per affected integrand
- updated `display integrands` detail metadata with:
  - generation-time compile enabled/backend/options
  - actual active f64 backend
- updated non-`OLD` run cards/examples to `compilation_mode`
- regenerated the tracked schemas in `assets/schemas/`
- updated `docs/architecture/architecture-current.md`

Validation run so far:
- `cargo fmt`
- `cargo check`
- `cargo test -p gammalooprs compile_settings_ -- --nocapture`
- `cargo test -p gammaloop-api string_toml_supports_multiline_nested_tables -- --nocapture`
- `GAMMALOOP_SCHEMA_PATH=assets/schemas cargo run -p gammaloop-api -- -n save schema`
- `cargo clippy` (currently passes with pre-existing workspace warnings)

## Known open issue

- Validation on the `gg > hhh` 1L amplitude currently fails because the symjit
  backend does not yet reproduce the same value as the existing backend. This
  branch/PR remains work in progress until that discrepancy is understood and
  fixed.

### Current debug narrowing

- The mismatch is isolated to the `original` bare integrand contribution in the
  `gg > hhh` 1L inspect output; threshold counterterms match between assembly
  and symjit.
- Temporary debug instrumentation in
  `crates/gammalooprs/src/integrands/process/amplitude/mod.rs` shows that the
  discrepancy is already present in the original single-parametric evaluator.
- Eager and assembly match each other for that evaluator, while symjit differs
  by an almost exact global complex rescaling.
- For the inspected sample/orientations, the symjit/original-to-eager ratio is
  consistently about `97.4090910340024`, suggesting one missing or extra global
  multiplicative factor in the original symjit evaluator path.
- The summed-function-map path is not the source of the issue: eager summed and
  eager single-parametric match, and symjit summed and symjit single-parametric
  also match.
- The temporary debug code is intentionally local-only and should not be pushed
  as part of the PR.

## Temporary local debug edits to revert

- `crates/gammalooprs/src/integrands/process/amplitude/mod.rs` contains
  env-gated local debug instrumentation under `GL_SYMJIT_DEBUG_ORIGINAL` to
  compare eager/assembly/symjit results for the original bare integrand.
- The same file now also honors `GL_SYMJIT_DEBUG_INPUT_DUMP=<path>` and writes
  the live original-evaluator input vector as JSON `[[re, im], ...]` for
  standalone reproduction.
- These two debugging hooks are temporary and should be removed once the symjit
  mismatch is understood.

## Standalone debugging status

- The standalone rust-script loader generated by `save standalone --rust --json`
  now supports backend switching via `--backend` / `--compare-backends` and a
  custom evaluator input via `--input-json <path>`.
- The rust-script header now disables Symbolica default features so the
  generated standalone script does not pull in `mimalloc`; this avoids the
  macOS linker failure that previously blocked running the script directly.
- On the archive's built-in representative input, assembly and symjit agree to
  rounding noise for `GL15`, so that input is not a useful reproducer.
- On the exact inspect-sample input dumped from GammaLoop, the standalone JSON
  reproducer does show a backend discrepancy for the original `GL15`
  `summed_function_map` evaluator: eager/assembly reproduce the inspect value,
  while symjit returns zero.
- The current standalone reproducer command is:
  `rust-script standalone_evaluators_rust.rs standalone_evaluators.json --graph-name GL15 --compare-backends eager,symjit --method summed_function_map --input-json /tmp/gg_hhh_inspect_input.json`
- This gives a simpler path for minimization: edit the JSON expression payload
  and rerun the standalone script until the smallest backend-divergent example
  is isolated.

## Current minimal reproducer

- The standalone minimization reduced the failure to a dead-else `if`
  expression with no GammaLoop-specific function-map machinery:
  `if(o,0,s)`.
- Evaluated at the fixed point `o = 1`, `s = 1`, eager Symbolica returns
  exactly zero while symjit returns a large non-zero garbage value.
- The exact symjit garbage value is not stable across runs, which strengthens
  the suspicion that this is a dead-branch codegen/runtime bug rather than a
  deterministic arithmetic mismatch.
- The same minimal expression also keeps the assembly backend at zero, so the
  divergence remains isolated to symjit.
- A dedicated self-contained reproducer now lives in the top-level directory
  `MRE_symjit_vs_symbolica/`:
  - `MRE_symjit_vs_symbolica/run.rs`
  - `MRE_symjit_vs_symbolica/expression.txt`
  - `MRE_symjit_vs_symbolica/README.md`
- Verified commands:
  - `cd MRE_symjit_vs_symbolica && rust-script run.rs expression.txt`
  - `cd MRE_symjit_vs_symbolica && rust-script run.rs expression.txt --assembly`
- The rust script exits with status `1` when the mismatch is reproduced, so it
  can be used directly as a regression probe while reporting the eager,
  optional assembly, symjit, and diff values.

## Additional local regression work: Symbolica assembly on `epem_a_tth`

- I reproduced the separate assembly regression reported on
  `examples/cli/epem_a_ttxh/NLO/epem_a_tth_NLO.toml`.
- Inspect result with eager evaluation:
  `+1.1622074088447979e-38 -5.6387400781737620e-55i`.
- Inspect result with compiled assembly evaluation:
  `+0.0000000000000000e0 +0.0000000000000000e0i`.
- With temporary local debug instrumentation in
  `crates/gammalooprs/src/integrands/process/cross_section/mod.rs`, I narrowed
  the mismatch to the pass-one original contribution of graph `GL08`; the
  threshold counterterm stays zero.
- The same debug hook dumps the exact pass-one evaluator input to
  `/tmp/epem_pass1_input.json` and the pass-two input to
  `/tmp/epem_pass2_input.json` when these env vars are set:
  - `GL_SYMBOLICA_REGRESSION_DEBUG`
  - `GL_SYMBOLICA_REGRESSION_PASS_ONE_INPUT_DUMP=<path>`
  - `GL_SYMBOLICA_REGRESSION_PASS_TWO_INPUT_DUMP=<path>`

## `MRE_symbolica_regression`

- I created a top-level standalone reproducer in
  `MRE_symbolica_regression/`.
- Files:
  - `MRE_symbolica_regression/payload.json`
  - `MRE_symbolica_regression/generated_assembly.cpp`
  - `MRE_symbolica_regression/run.rs`
  - `MRE_symbolica_regression/README.md`
- `payload.json` is the pass-one
  `raised_cut[0].derivative[0].summed_function_map` evaluator extracted from
  `save standalone --json`, with the exact inspect input baked in.
- `generated_assembly.cpp` is the exact assembly backend source emitted by
  GammaLoop for that evaluator when generated with `compile=true` and
  `compilation_mode=assembly`.
- Running `cd MRE_symbolica_regression && rust-script run.rs` rebuilds the
  eager evaluator from `payload.json`, compiles `generated_assembly.cpp`, and
  compares the two results.
- Verified result:
  - eager: finite non-zero
  - generated assembly source: `NaN`
- Running `rust-script run.rs --fresh-compile` also shows that a fresh assembly
  build from `payload.json` itself stays finite, so the divergence is already
  present in the generated assembly source saved by GammaLoop.

## Saved-state reload bug

- I fixed the unrelated saved-state reload bug that existed on `main`.
- Root cause: process loading treated compile-artifact directories such as
  `cs_NLO/` as persisted cross sections and then tried to load a missing
  `cs.bin`.
- Fix: `Process::load_amplitude` and `Process::load_cross_section` now skip
  subdirectories that do not contain the corresponding persisted manifest file
  (`amp.bin` / `cs.bin`).
- This keeps process loading manifest-based and avoids mistaking compile-output
  scratch directories for saved amplitudes or cross sections.

## Generation summary table

- I added a generation-time summary table that is emitted once integrand
  generation completes.
- The summary is rendered from `crates/gammaloop-api/src/commands/generate.rs`
  with `tabled` and the existing CLI color scheme.
- Columns now include:
  - integrand
  - graph
  - total evaluator count built for that graph
  - expression-build time, with percentage of total graph generation time
  - evaluator-build time, with percentage of total graph generation time
  - compile time, with percentage of total graph generation time
- The timing model is:
  - preprocess/build-expression work contributes to total generation time
  - Symbolica evaluator construction contributes to evaluator-build time
  - backend compilation contributes to compile time
  - expression-build time is reported as the residual:
    `total - evaluator_build - evaluator_compile`
- Reporting covers:
  - amplitude graphs
  - cross-section graphs
  - symjit runtime JIT activation after generation
  - external compilation for `assembly` / `c++`
- I added non-persisted reporting types in
  `crates/gammalooprs/src/processes/generation_report.rs` and threaded them
  through:
  - `State`
  - `ProcessList`
  - `Process`
  - amplitude/cross-section generation paths
- For cross sections, pass-two evaluator construction is accounted for during
  preprocess, while the raised-cut and threshold evaluator stacks are accounted
  for during integrand construction.
- For amplitudes, the original stack and threshold-counterterm stacks are
  accounted for during integrand construction.

### Validation

- `cargo fmt`
- `cargo check`
- Verified the new summary table on the `gg_hhh` 1L amplitude with:
  - eager generation:
    `cargo run -p gammaloop-api --bin gammaloop -- --clean-state ./examples/cli/gg_hhh/1L/gg_hhh_1L.toml run -c "set global kv global.generation.evaluator.compile=false; run generate; generate existing -p 0 -i 1L; quit -n"`
  - symjit generation:
    `target/debug/gammaloop --clean-state ./examples/cli/gg_hhh/1L/gg_hhh_1L.toml run -c "set global kv global.generation.evaluator.compile=true; set global kv global.generation.compile.compilation_mode=symjit; run generate; generate existing -p 0 -i 1L; quit -n"`
- Observed the requested table with non-zero compile time in symjit mode.

### Clean commit preparation

- Removed the temporary local-only debug instrumentation from the tracked source
  files before preparing the next pushed commit.
- The clean commit scope is limited to the permanent changes:
  - post-generation timing summary table
  - per-graph evaluator-count and timing aggregation
  - saved-state reload fix that ignores artifact-only directories without
    `amp.bin` / `cs.bin`
- Debug-only local artifacts and MRE directories remain uncommitted.
