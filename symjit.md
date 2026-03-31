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
   to symjit, GammaLoop must emit one `info!`.
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

At `26e0b19f50135bc6bb66d70eebbe20033d49b9d1`, Symbolica exposes both:
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
  `26e0b19f50135bc6bb66d70eebbe20033d49b9d1`
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
  - emit one `info!` naming the integrand and the reason for fallback
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
  `26e0b19f50135bc6bb66d70eebbe20033d49b9d1`
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
  loads and log one `info!` per affected integrand
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
