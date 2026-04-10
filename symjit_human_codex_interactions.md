# symjit human <-> codex interactions

## Interaction 1

### User prompt

The user asked me to:
- work on branch `symjit` in the current `gammaloop_dev_branch` repo
- first inspect the codebase and reason about how to add a new
  `compilation_mode` setting with `eager`, `c++`, `assembly`, and `symjit`
- use Symbolica/symjit docs as needed
- write and maintain the implementation plan in `symjit.md`
- keep this interaction log updated through the session
- come back first with the plan for review before implementation

### Assistant final answer

I went through the settings layer, evaluator stack, generate/load flow, and the
current Symbolica sources, and I wrote the working branch plan in `symjit.md`.

Main findings:
- `compilation_mode` has to become the new source of truth because backend
  choice is currently split between
  `crates/gammalooprs/src/settings/global.rs` (`inline_asm`) and
  `crates/gammalooprs/src/processes/mod.rs` / 
  `crates/gammaloop-api/src/commands/generate.rs` (`evaluator.compile`)
- The current Symbolica pin in `Cargo.toml`
  (`650ba97bf3da7cf2ff5ada92875f92d5f71e7a31`) does not expose symjit yet, so
  this branch needs a Symbolica/graphica/numerica rev bump first
- GammaLoop already saves enough to rebuild symjit cleanly:
  `GenericEvaluator` stores the rational evaluator plus the eager evaluators in
  `crates/gammalooprs/src/integrands/process/evaluators.rs`
- The f64 dispatch path currently uses "compiled if present, otherwise eager",
  which is no longer enough; we need an explicit runtime backend selection so
  `eager`, `c++`/`assembly`, and `symjit` are chosen by mode, not by whichever
  optional backend happens to be populated
- State load currently does not receive global settings, so load-time symjit
  rebuild needs either:
  - threading the effective global settings into `State::load(...)`, or
  - an immediate post-load initialization pass

Implementation plan:
1. Bump the patched Symbolica revs so `jit_compile(...)` is available.
2. Add `CompilationMode` under `GammaloopCompileOptions`, default it to
   `symjit`, and make it the new backend selector.
3. Add a runtime-only symjit wrapper/cache in `GenericEvaluator`.
4. Make `GenericEvaluator` track the active f64 backend explicitly.
5. Initialize symjit right after evaluator construction and rebuild it again on
   state load.
6. Restrict `compile_integrands(...)` to the external backends (`c++` and
   `assembly`); eager and symjit should not trigger filesystem compilation.
7. Update schemas, tests, run cards, and examples to the new setting model.

Main review questions before I start implementation:
1. Do you want a legacy compatibility shim for
   `generation.evaluator.compile` and `generation.compile.inline_asm`, or should
   I hard-migrate the repo fixtures and treat old settings as a clean break?
2. Should I apply symjit uniformly to every runtime `GenericEvaluator` user
   (LU counterterms, overlap helpers, pass-two evaluators, etc.), or only to
   the main evaluator stacks? My recommendation is to do it uniformly at the
   `GenericEvaluator` level.

## Interaction 2

### User prompt

The user clarified that:
- the Symbolica pin should move to
  `26e0b19f50135bc6bb66d70eebbe20033d49b9d1`
- `compile=true/false` should remain the eager-vs-compiled switch
- the new setting should only absorb/deprecate `inline_asm`
- symjit should be applied uniformly across all evaluator stacks
- `to_symbolica_compile_options()` can be fixed or removed, whichever is
  cleaner
- the plan should be updated
- I should expand on the post-load symjit construction idea and come back with
  any remaining clarification points

### Assistant final answer

I updated `symjit.md` to reflect the revised structure:
- `global.generation.evaluator.compile` stays authoritative for eager vs
  compiled
- the new selector under `global.generation.compile` is only for compiled
  backends (`c++`, `assembly`, `symjit`)
- the Symbolica bump target is fixed to
  `26e0b19f50135bc6bb66d70eebbe20033d49b9d1`
- symjit is applied uniformly at the `GenericEvaluator` level

My proposed post-load solution is:
1. Keep `State::load(...)` focused on pure deserialization.
2. Add an explicit backend refresh pass, e.g.
   `State::refresh_evaluator_backends(&GlobalSettings)`.
3. Call it immediately after `State::load(...)` returns in
   `crates/gammaloop-api/src/lib.rs`, using the already-known startup global
   settings.
4. Reuse that same refresh entry point after generation, and possibly after
   compilation-related `set global ...` changes if we want existing evaluators
   to switch backend immediately.

The main additional subtlety I found is that externally compiled evaluators are
currently decoded eagerly: deserializing a saved compiled evaluator tries to
reload its shared library right away. That means a state can fail during load
before the post-load refresh pass has a chance to choose eager or symjit
instead. So I likely need to make external compiled loading lazy/tolerant, or
otherwise change that decode path, to keep the "deserialize first, choose
backend after load" design coherent.

The remaining clarification points are:
1. Should `compilation_mode = "eager"` still be accepted as a deprecated input
   alias, or should it disappear entirely now that `evaluator.compile` remains
   the eager-vs-compiled switch?
2. Should changing compilation settings through `set global ...` immediately
   reconfigure already-loaded/generated evaluators, or should those changes only
   apply to future load/generate actions?
3. Are you okay with me changing the external compiled-evaluator loading model
   to a lazier/tolerant form if that is the cleanest way to support automatic
   post-load backend selection?

## Interaction 3

### User prompt

The user clarified that:
- `compilation_mode = "eager"` should not exist anymore
- existing integrands are locked to the backend they were generated with
- current global compilation settings must affect only future generated
  integrands, never already-existing ones
- existing integrands loaded from disk should restore their backend choice from
  an integrand-side manifest, not from current globals
- only symjit should do an extra rebuild on load, because it is intentionally
  not serialized
- I should reassess the plan, inspect the code again, critique it, and come
  back with an updated review plan

### Assistant final answer

I re-audited the persistence and generation paths and updated `symjit.md`
accordingly.

The main correction is that process-level settings history is not sufficient for
backend restoration:
- `Process` only stores one `settings_history.toml`
- one process can contain multiple generated integrands
- targeted integrand regeneration can rebuild a single integrand without
  updating that process-level history

So if backend choice is truly locked per integrand, the authoritative metadata
must be persisted per integrand.

The revised plan is:
1. Bump Symbolica to `26e0b19f50135bc6bb66d70eebbe20033d49b9d1`.
2. Keep `global.generation.evaluator.compile` as the eager-vs-compiled switch.
3. Replace `inline_asm` with a compiled-backend selector
   (`c++`, `assembly`, `symjit`).
4. Persist a per-integrand frozen-backend manifest with each generated
   integrand.
5. Lock each generated integrand to that manifest at generation time.
6. Keep eager and external compiled integrands loading as they do now.
7. Rebuild only the symjit backend on load, using the saved eager/rational
   evaluator plus the integrand's saved manifest.
8. Do not let current global compilation settings reinterpret existing
   integrands.

The main extra critique point is that `AmplitudeIntegrand::save/load` and
`CrossSectionIntegrand::save/load` currently have no per-integrand generation
manifest at all, so that has to be introduced explicitly. After another code
review, the lightest place to store it appears to be inside
`AmplitudeIntegrandData` / `CrossSectionIntegrandData`, because those structs
already define `integrand.bin` and are already used by
`ProcessIntegrand::resume_fingerprint()`.

The remaining clarification points are:
1. Should the new per-integrand manifest store a resolved backend enum
   (`Eager`, `Cpp`, `Assembly`, `Symjit`) or a raw snapshot of
   `compile + GammaloopCompileOptions`? I recommend a resolved manifest.
2. Do you need any compatibility fallback for older saved states that do not
   have this new per-integrand manifest, or is it acceptable that the new
   manifest becomes the authoritative format going forward?
3. For `inline_asm`, should I keep a deprecated read-only parse shim and stop
   serializing it, or do you want a hard break immediately?

## Interaction 4

### User prompt

The user clarified that:
- the manifest should contain only fields that freeze the generated integrand
- there should be no backward compatibility
- `inline_asm` is fully gone and absorbed into the new setting
- non-`OLD` toml files should be updated to the new settings standard
- I should do one more deep review of the relevant code and come back with the
  updated plan and any remaining clarification needs

### Assistant final answer

I did another review of the save/load and evaluator ownership paths and updated
`symjit.md` again.

The main new design refinement is that the frozen backend metadata should
probably live as an embedded manifest inside `integrand.bin`, not as an extra
sidecar file. Concretely, the best place appears to be:
- `AmplitudeIntegrandData`
- `CrossSectionIntegrandData`

That is cleaner because:
- those are already the structs serialized into `integrand.bin`
- their load path already runs for every integrand
- `ProcessIntegrand::resume_fingerprint()` already hashes `integrand.data`, so
  the frozen backend would automatically become part of the resume identity

The updated plan is therefore:
1. Bump Symbolica to `26e0b19f50135bc6bb66d70eebbe20033d49b9d1`.
2. Keep `evaluator.compile` as the eager-vs-compiled switch.
3. Replace `inline_asm` with the new compiled-backend selector
   (`c++`, `assembly`, `symjit`) and update non-`OLD` tomls.
4. Persist only the resolved frozen backend per generated integrand, ideally as
   an embedded enum like `Eager | Cpp | Assembly | Symjit`.
5. Apply backend initialization uniformly across all `GenericEvaluator` owners.
6. Rebuild only symjit on load, based on that embedded frozen-backend manifest.
7. Keep eager and external compiled loading otherwise as close as possible to
   current behavior.

The main remaining risk I want to keep in mind is:
- if an integrand is locked to `Cpp` or `Assembly` but somehow its serialized
  compiled evaluator is missing/unloadable, the natural behavior is to treat
  that as an invalid/incomplete saved state rather than silently changing
  backend

At this point the only possible clarification I still see is whether you want
the frozen manifest stored as:
1. a resolved enum only, which I recommend, or
2. some larger structured manifest

My recommendation remains the resolved enum only, since you explicitly want to
store only what freezes the generated integrand.

## Interaction 5

### User prompt

The user clarified that:
- they are happy with embedding the frozen backend information directly inside
  the serialized integrand `.bin`
- if `.so` loading fails or the architecture is wrong, they want an automatic
  fallback to symjit, but only if the global settings at load time opt into
  symjit
- I should review that fallback design critically and return with an updated
  plan and comments

### Assistant final answer

After reviewing the current compiled-evaluator decode path, I think the symjit
fallback is supportable, but not by keeping the current "serialize a loaded
compiled evaluator directly" design.

The relevant issue is:
- `CompiledComplexEvaluator` in Symbolica decodes by immediately trying to load
  the shared library
- `CompiledComplexEvaluatorSpenso` currently inherits that eager behavior
- so if we keep the current serialized shape, a bad/missing/wrong-arch `.so`
  fails during bincode decode before GammaLoop has any chance to apply the
  fallback policy

So the clean implementation is likely:
1. Stop serializing the loaded external compiled evaluator directly.
2. Serialize only a lightweight external artifact descriptor, essentially
   `(path, function_name)`.
3. Keep the actually loaded `.so` evaluator in a runtime-only cache.
4. During integrand load, use the embedded frozen-backend enum as the normal
   selector.
5. If that frozen backend is `Cpp` or `Assembly`, explicitly try to load the
   external artifact.
6. If that load fails and startup globals explicitly opt into symjit, log an
   info and build symjit instead.
7. Otherwise fail the load as an invalid/incomplete saved state.

I think this is the cleanest way to satisfy your fallback requirement while
keeping the intended semantics:
- frozen backend remains embedded in the integrand and is not rewritten
- current globals are not used to reinterpret existing integrands in general
- current globals are used only as a narrow fallback permission when external
  compiled loading fails

My recommendation for the fallback gate is:
- require `global.generation.evaluator.compile = true`
- require the startup compiled backend selector to be `symjit`

That seems like the least ambiguous reading of "globals say symjit".

The only clarification I still see is whether you agree with that exact
fallback gate, i.e. whether `compile=false` should disable the fallback even if
the compiled-backend selector field still happens to say `symjit`.

## Interaction 6

### User prompt

The user clarified that:
- they are happy with embedding the extra frozen-backend information directly in
  the serialized integrand `.bin`
- they want `.so` load failure to fallback to symjit only if the global
  settings at load time say symjit
- I should review the plan accordingly in depth and come back with comments and
  clarification

### Assistant final answer

I think the fallback is viable, but only with one important structural change:
the external compiled evaluator can no longer stay serialized as an eagerly
loaded object.

The clean plan is:
1. Embed only the frozen backend enum directly into `AmplitudeIntegrandData` /
   `CrossSectionIntegrandData`.
2. Store external compiled backends as serialized artifact metadata
   `(path, function_name)`, not as an already-loaded evaluator.
3. Keep the actual loaded `.so` evaluator in a runtime-only cache.
4. During load, use the embedded frozen backend as the normal selector.
5. If that frozen backend is external compiled, explicitly try to load the
   artifact.
6. If loading fails and startup globals opt into symjit, log an info and build
   symjit instead.
7. Otherwise fail the load.

I also re-checked the load call chain and it is reasonably contained:
`lib.rs -> state.rs -> ProcessList::load -> Process::load_* ->
Amplitude/CrossSection::load -> *Integrand::load`.

So threading a narrow startup load policy down that chain is manageable. The
important point is that current globals are still not used as the normal backend
selector for existing integrands; they are only used as an opt-in permission for
the external-compiled-load failure fallback.

My recommended fallback gate remains:
- `global.generation.evaluator.compile = true`
- startup compiled backend selector = `symjit`

That is the least ambiguous reading of "globals say symjit".

## Interaction 7

### User prompt

The user clarified that:
- deferring actual `.so` loading into a runtime cache is acceptable given the
  post-load activation step
- GammaLoop must emit one `info!` for each integrand that actually falls back
  from external compiled loading to symjit
- `display integrands` metadata in the first table should show:
  - the generation-time compilation mode/options available at load time
  - the actual f64 backend currently active (`eager`, `c++`, `assembly`, or
    `symjit`)
- I should do one last global review of the plan and come back with any final
  clarifications if needed

### Assistant final answer

I did the last global pass through the settings, evaluator, load, and display
paths, and I updated `symjit.md` again.

The main correction from this review is that the embedded per-integrand frozen
metadata cannot stay just a bare backend enum if `display integrands` must show
accurate generation-time compile options after load. Since process-level
`settings_history.toml` is not precise per integrand, the frozen metadata should
contain only fields that truly froze the generated evaluator state, but for
external backends that includes the external compile options themselves.

So the current recommended embedded metadata shape is:
- `Eager`
- `Symjit`
- `Cpp { optimization_level, fast_math, unsafe_math, compiler, custom }`
- `Assembly { optimization_level, fast_math, unsafe_math, compiler, custom }`

I also think runtime activation must be handled integrand-wide, not separately
per `GenericEvaluator`. Otherwise a partial `.so` load failure could leave one
integrand in a mixed runtime state, which would make the requested single
"active f64 backend" display ambiguous. The cleaner rule is:
- frozen metadata is stored per integrand
- backend activation/fallback runs across all `GenericEvaluator`s under that
  integrand
- each integrand ends up with one active f64 flavor
- if any external artifact load fails and startup globals authorize symjit, the
  whole integrand falls back to symjit and one `info!` is emitted for that
  integrand

One other refinement is that Symbolica's newer `CompiledCode<T>` type is a
better persisted external-artifact representation than a home-grown path/function
pair, because it already serializes lazily and loads only when GammaLoop asks.

At this point I do not see another blocking clarification. The only non-blocking
implementation assumption I will use is that "the first table" means the detail
summary table rendered by `display integrands -p ... -i ...`, not the top-level
integrand list table.

## Interaction 8

### User prompt

The user instructed me to proceed with the implementation.

### Assistant final answer

I implemented the symjit branch changes and updated `symjit.md` to reflect the
finished plan and validation state.

Implemented:
- bumped the patched `symbolica`, `graphica`, and `numerica` revisions to
  `26e0b19f50135bc6bb66d70eebbe20033d49b9d1`
- replaced `global.generation.compile.inline_asm` with
  `global.generation.compile.compilation_mode`
- kept `global.generation.evaluator.compile` as the eager-vs-compiled switch
- made `symjit` the default compiled backend
- embedded frozen compilation metadata directly into the serialized integrand
  data
- switched persisted external compiled evaluators to lazy
  `CompiledCode<Complex<f64>>` artifacts with runtime-only loaded-library caches
- added runtime-only symjit caches and explicit active-f64-backend tracking
- restored runtime backends after load from frozen integrand metadata
- added the startup-gated symjit fallback for failed external `.so` loads and
  emit one `info!` per affected integrand
- updated `display integrands` detail metadata with generation-time compilation
  details and the active f64 backend
- updated non-`OLD` run cards/examples to `compilation_mode`
- regenerated `assets/schemas/*.json`
- updated `docs/architecture/architecture-current.md`

Validation run:
- `cargo fmt`
- `cargo check`
- `cargo test -p gammalooprs compile_settings_ -- --nocapture`
- `cargo test -p gammaloop-api string_toml_supports_multiline_nested_tables -- --nocapture`
- `GAMMALOOP_SCHEMA_PATH=assets/schemas cargo run -p gammaloop-api -- -n save schema`
- `cargo clippy`

`cargo clippy` completes successfully, but the workspace still has many
pre-existing warnings outside this task.

## Interaction 9

### User prompt

The user asked me to:
- prepare a short-title / long-body commit message
- commit the staged change
- push branch `symjit` to `origin`
- open a draft PR to `main` with `gh`
- make the PR body follow `symjit.md`, with a cleaner GitHub-facing summary
- explicitly note that validation on the `gg > hhh` 1L amplitude currently
  fails because symjit does not reproduce the same value yet

### Assistant final answer

I updated `symjit.md` with the open `gg > hhh` 1L validation issue, committed
the staged symjit work with a short title and detailed body, pushed `symjit` to
`origin`, and opened a draft PR to `main`.

The PR body summarizes:
- the new `compilation_mode` setting replacing `inline_asm`
- the frozen per-integrand compilation metadata
- the runtime symjit / external backend activation model
- the narrow symjit fallback on failed external artifact loading
- the `display integrands` metadata updates
- the validation commands run
- the current known discrepancy on the `gg > hhh` 1L amplitude, so the PR is
  clearly marked as work in progress

## Interaction 10

### User prompt

The user asked me to debug the symjit discrepancy on the `gg > hhh` 1L process
by comparing the rich `inspect` output for one fixed sample between:
- `global.generation.compile.compilation_mode=assembly`
- `global.generation.compile.compilation_mode=symjit`

They pointed out that the mismatch is isolated to the `original` contribution,
while the threshold counterterms match, and asked me to add temporary debug
instrumentation to narrow down the origin of that difference.

### Assistant final answer

I added temporary local-only debug instrumentation in
`crates/gammalooprs/src/integrands/process/amplitude/mod.rs` to compare the
original bare-integrand evaluator under eager, assembly, and symjit for the
same inspected sample and orientations.

The narrowing result is:
- the discrepancy is not in the threshold counterterms; they still match
- the discrepancy is not in the summed-function-map aggregation logic
- the discrepancy is already present in the original single-parametric bare
  evaluator
- eager and assembly match each other exactly for that evaluator
- symjit differs by an almost exact global multiplicative factor across all
  tested orientations

For the inspected sample, the symjit/eager ratio for the original contribution
is consistently about `97.4090910340024`, which strongly suggests one missing or
extra global multiplicative factor in the original symjit evaluator path.

So the current narrowed conclusion is:
- the bug is in symjit evaluation of the original single-parametric bare
  integrand for `gg > hhh` 1L
- it is not caused by backend restoration, load-time activation, threshold
  counterterms, or summed-function-map accumulation

The temporary debug changes are intentionally not for commit/push.

## Interaction 11

### User prompt

The user asked me to:
- keep track of the temporary debugging changes so they can be undone later
- improve `save standalone --json` so the generated standalone script can switch
  cleanly between Symbolica assembly compilation and symjit
- verify whether the discrepancy persists in standalone form
- use that improved standalone path to make future minimization easier, since
  the JSON payload can then be simplified into a small backend-divergent example

### Assistant final answer

I updated the amplitude standalone rust-script loader used by
`save standalone --rust --json` so it now supports:
- `--backend <eager|c++|assembly|symjit>`
- `--compare-backends <backend[,backend...]>`
- `--method <single_parametric|iterative|summed_function_map|summed>`
- `--graph-index` / `--graph-name`
- `--stack <original|ct:N>`
- `--input-json <path>` for an explicit evaluator input vector stored as JSON
  `[[re, im], ...]`

I also changed the generated rust-script cargo snippet to use
`symbolica` with `default-features = false`, so the standalone script no longer
pulls in `mimalloc`; this avoids the macOS linker failure that was blocking the
standalone execution path.

For the temporary local debugging state, I kept and documented two env-gated
hooks in `crates/gammalooprs/src/integrands/process/amplitude/mod.rs`:
- `GL_SYMJIT_DEBUG_ORIGINAL`
- `GL_SYMJIT_DEBUG_INPUT_DUMP=<path>`

Using those hooks, I dumped the exact live original-evaluator input vector for
the failing `inspect -x 0.1 0.2 0.3 -d 0 0` sample into
`/tmp/gg_hhh_inspect_input.json`.

Standalone verification result:
- on the archive's built-in representative input, assembly and symjit agree to
  rounding noise, so that input is not a useful reproducer
- on the exact dumped inspect input, the standalone JSON reproducer does show a
  backend discrepancy for `GL15` on the original `summed_function_map`
  evaluator
- in that standalone reproducer, eager and assembly match the inspect value,
  while symjit returns zero

So the standalone path is now useful for minimization: simplify the JSON
expression payload and rerun the standalone script until the smallest backend-
divergent expression is isolated.

## Interaction 12

### User prompt

The user clarified that the zero from the standalone reproducer was expected
because it already used a reduced path, then asked me to:
- keep track of the temporary debug edits so they can be reverted later
- create a top-level standalone reproducer in `MRE_symjit_vs_symbolica/`
- reduce the expression as aggressively as possible
- keep everything in one Rust script except for the separate input expression
  file
- add a short `README.md` explaining how to run it

### Assistant final answer

I reduced the standalone reproducer down to a dead-else `if` expression with no
GammaLoop-specific function-map dependency:

- expression: `if(o,0,s)`
- evaluation point: `o = 1`, `s = 1`

That is the current smallest expression I found that still reproduces the
symjit mismatch clearly:
- eager Symbolica returns `0`
- assembly also returns `0`
- symjit returns a large non-zero garbage value

I added a dedicated top-level reproducer in `MRE_symjit_vs_symbolica/`:
- `MRE_symjit_vs_symbolica/run.rs`
- `MRE_symjit_vs_symbolica/expression.txt`
- `MRE_symjit_vs_symbolica/README.md`

The script is self-contained apart from the expression file, reads the
expression from disk, evaluates it with eager and symjit, optionally compares
against Symbolica's compiled assembly backend via `--assembly`, prints the
results and diff, and exits with status `1` when the mismatch is reproduced.

I verified the reproducer with:
- `cd MRE_symjit_vs_symbolica && rust-script run.rs expression.txt`
- `cd MRE_symjit_vs_symbolica && rust-script run.rs expression.txt --assembly`

I also kept the existing temporary debug hooks in
`crates/gammalooprs/src/integrands/process/amplitude/mod.rs` documented for
later rollback:
- `GL_SYMJIT_DEBUG_ORIGINAL`
- `GL_SYMJIT_DEBUG_INPUT_DUMP=<path>`

## Interaction 13

### User prompt

The user asked me to continue after a pause and clarified that, for now, the
`epem_a_tth` assembly-regression reproducer does not need a smallest possible
expression. Using the full integrand expression from `save standalone --json`
is good enough.

### Assistant final answer

I continued from the existing regression-debugging state.

Current state:
- I had already reproduced the `epem_a_tth` eager-vs-assembly mismatch and the
  saved-state reload crash.
- I had narrowed the assembly regression to the pass-one original contribution
  of graph `GL08` and dumped the exact pass-one input.
- I was about to turn that into a top-level standalone reproducer and then fix
  the reload bug.

Planned next actions:
1. Build `MRE_symbolica_regression/` from the exported standalone payload and
   the dumped inspect input.
2. Fix the reload bug by making process loading ignore non-persisted compile
   artifact directories.
3. Update `symjit.md` / `symjit_human_codex_interactions.md` and rerun
   formatting and checks.

## Interaction 14

### User prompt

The user asked me to continue.

### Assistant final answer

I continued the implementation and validation.

What I added:
- A new top-level regression reproducer in `MRE_symbolica_regression/`:
  - `payload.json`
  - `generated_assembly.cpp`
  - `run.rs`
  - `README.md`
- The reproducer uses the pass-one
  `raised_cut[0].derivative[0].summed_function_map` evaluator extracted from
  `save standalone --json`, with the exact inspect input baked in.
- `run.rs` rebuilds the eager evaluator from `payload.json`, compiles the
  exact `generated_assembly.cpp` emitted by GammaLoop, and compares the two.
- Verified result:
  - eager evaluation is finite and non-zero
  - the generated assembly source evaluates to `NaN`
- `run.rs --fresh-compile` also shows that a fresh assembly build from the
  payload itself stays finite, so the regression is already present in the
  generated assembly source that GammaLoop writes.

I also fixed the unrelated saved-state reload bug:
- `crates/gammalooprs/src/processes/process.rs` now skips subdirectories that
  do not contain the persisted manifest file (`amp.bin` or `cs.bin`) when
  loading amplitudes / cross sections.
- This prevents compile-artifact directories such as `cs_NLO/` from being
  misidentified as persisted cross sections.

Validation performed:
- `cd MRE_symbolica_regression && rust-script run.rs`
- `cd MRE_symbolica_regression && rust-script run.rs --fresh-compile`
- `cargo fmt`
- `cargo check`
- reproduced the reload flow with:
  - create state: `cargo run -p gammaloop-api --bin gammaloop -- -s <tmp> --clean-state ./examples/cli/epem_a_ttxh/NLO/epem_a_tth_NLO.toml run -c "run generate_diagrams generate_integrands; quit -o;"`
  - reopen state: `cargo run -p gammaloop-api --bin gammaloop -- -s <tmp> ./examples/cli/epem_a_ttxh/NLO/epem_a_tth_NLO.toml run -c "quit -n"`
- The original `No such file or directory (os error 2)` crash no longer occurs.

## Interaction 15

### User prompt

The user asked me to enhance the closing message after integrand generation by
printing a colorized `tabled` summary with columns for:
- integrand
- graph
- total evaluator count
- cumulative input-expression build time with percentage of total graph time
- cumulative Symbolica evaluator-build time with percentage
- cumulative evaluator compile time with percentage

### Assistant final answer

I implemented the generation summary table and verified it.

What changed:
- Added generation-report plumbing in
  `crates/gammalooprs/src/processes/generation_report.rs` and threaded it
  through `State`, `ProcessList`, `Process`, and the amplitude/cross-section
  generation paths.
- The table is rendered from
  `crates/gammaloop-api/src/commands/generate.rs` using `tabled` with the
  existing CLI colors.
- The reported timing breakdown is:
  - expression build: total graph generation time minus evaluator build and
    compile time
  - evaluator build: Symbolica evaluator construction only
  - compile: symjit JIT activation or external `assembly` / `c++` compilation
- Evaluator counts are accumulated per generated graph across the relevant
  stacks for amplitudes and cross sections.

Validation performed:
- `cargo fmt`
- `cargo check`
- Verified the table on the `gg_hhh` 1L amplitude with:
  - eager mode:
    `cargo run -p gammaloop-api --bin gammaloop -- --clean-state ./examples/cli/gg_hhh/1L/gg_hhh_1L.toml run -c "set global kv global.generation.evaluator.compile=false; run generate; generate existing -p 0 -i 1L; quit -n"`
  - symjit mode:
    `target/debug/gammaloop --clean-state ./examples/cli/gg_hhh/1L/gg_hhh_1L.toml run -c "set global kv global.generation.evaluator.compile=true; set global kv global.generation.compile.compilation_mode=symjit; run generate; generate existing -p 0 -i 1L; quit -n"`
- The summary table appears as requested, and symjit mode shows a non-zero
  compile column.

## User prompt

Add all the changes that are *not* temporary debugging changes, commit them with a nice commit messages detailing the changes since last commit, and push them to the remote current active branch "symjit".

## Assistant final answer

I prepared a clean commit scope by removing the temporary debug instrumentation from the tracked source files, leaving only the permanent generation-summary and state-reload changes staged for commit. The debug-only local artifacts and MRE directories remain uncommitted.
