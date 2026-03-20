# Codex / Human Interaction Log

## 2026-03-20

### Human

````
# AGENTS.md instructions for /Users/vjhirsch/Documents/Work/gammaloop_vectorize

<INSTRUCTIONS>
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
- When exposing a function to within the crate only, prefer `pub(crate)` to more permissive `pub`.

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
- Always `cargo fmt` and `cargo check` before compiling to catch easy errors and ensure consistent formatting.
- API stability is currently a non-goal for internal development: prioritize maintainability and structure over preserving external-facing module paths.
- Breaking Rust/Python API changes are acceptable when they simplify architecture; document major moves in the change description.
- Backward compatibility with old on-disk GammaLoop states is also currently a non-goal: it is acceptable to rename/remove legacy state files, layouts, and loaders when that simplifies the implementation.
- Do not keep compatibility fallbacks for old state formats/names unless the task explicitly asks for migration support.
- Never be afraid to modify existing functions to obtain the most elegant and concise code as opposed to trying to be backward compatible or leave existing code untouched.
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

## Environment and Build Quirks
- `NO_SYMBOLICA_OEM_LICENSE` is read via `option_env!` in `gammalooprs` initialization, so it is a compile-time switch, not a pure runtime one. If you need to disable the OEM-license activation path, that variable must be present when building the relevant binary or Python extension.

## CLI / State Notes
- The run-card field is `command_blocks` (singular `command`), not `commands_blocks`. No compatibility alias is kept.
- Integration workspace resume must restore the exact state at the end of the last completed iteration. The authoritative resume data lives in the workspace manifest/internal state plus per-slot settings and internal observable checkpoints, not in user-facing `integration_result.json`.
- Treat workspace snapshots as completed-iteration checkpoints only: do not persist partial-iteration state, and keep user-facing latest snapshots plus optional `_iter_0001` history separate from the internal resume state.

## Configuration & State Tips
- CLI runs create a `gammaloop_state/` directory by default; keep it out of commits unless intentionally sharing a reproducible state.
- Use `./bin/gammaloop -s <state_dir>` for isolated experiments.


## Skills
A skill is a set of local instructions to follow that is stored in a `SKILL.md` file. Below is the list of skills that can be used. Each entry includes a name, description, and file path so you can open the source for full instructions when using a specific skill.
### Available skills
- openai-docs: Use when the user asks how to build with OpenAI products or APIs and needs up-to-date official documentation with citations, help choosing the latest model for a use case, or explicit GPT-5.4 upgrade and prompt-upgrade guidance; prioritize OpenAI docs MCP tools, use bundled references only as helper context, and restrict any fallback browsing to official OpenAI domains. (file: /Users/vjhirsch/.codex/skills/.system/openai-docs/SKILL.md)
- skill-creator: Guide for creating effective skills. This skill should be used when users want to create a new skill (or update an existing skill) that extends Codex's capabilities with specialized knowledge, workflows, or tool integrations. (file: /Users/vjhirsch/.codex/skills/.system/skill-creator/SKILL.md)
- skill-installer: Install Codex skills into $CODEX_HOME/skills from a curated list or a GitHub repo path. Use when a user asks to list installable skills, install a curated skill, or install a skill from another repo (including private repos). (file: /Users/vjhirsch/.codex/skills/.system/skill-installer/SKILL.md)
### How to use skills
- Discovery: The list above is the skills available in this session (name + description + file path). Skill bodies live on disk at the listed paths.
- Trigger rules: If the user names a skill (with $SkillName or plain text) OR the task clearly matches a skill's description shown above, you must use that skill for that turn. Multiple mentions mean use them all. Do not carry skills across turns unless re-mentioned.
- Missing/blocked: If a named skill isn't in the list or the path can't be read, say so briefly and continue with the best fallback.
- How to use a skill (progressive disclosure):
  1) After deciding to use a skill, open its `SKILL.md`. Read only enough to follow the workflow.
  2) When `SKILL.md` references relative paths (e.g., `scripts/foo.py`), resolve them relative to the skill directory listed above first, and only consider other paths if needed.
  3) If `SKILL.md` points to extra folders such as `references/`, load only the specific files needed for the request; don't bulk-load everything.
  4) If `scripts/` exist, prefer running or patching them instead of retyping large code blocks.
  5) If `assets/` or templates exist, reuse them instead of recreating from scratch.
- Coordination and sequencing:
  - If multiple skills apply, choose the minimal set that covers the request and state the order you'll use them.
  - Announce which skill(s) you're using and why (one short line). If you skip an obvious skill, say why.
- Context hygiene:
  - Keep context small: summarize long sections instead of pasting them; only load extra files when needed.
  - Avoid deep reference-chasing: prefer opening only files directly linked from `SKILL.md` unless you're blocked.
  - When variants exist (frameworks, providers, domains), pick only the relevant reference file(s) and note that choice.
- Safety and fallback: If a skill can't be applied cleanly (missing files, unclear instructions), state the issue, pick the next-best approach, and continue.
</INSTRUCTIONS><environment_context>
  <cwd>/Users/vjhirsch/Documents/Work/gammaloop_vectorize</cwd>
  <shell>zsh</shell>
  <current_date>2026-03-20</current_date>
  <timezone>Europe/Zurich</timezone>
</environment_context>
As you can see our `    pub fn evaluate_samples_raw(` API entry point in

[@mod.rs#L592:599](file:///Users/vjhirsch/Documents/Work/gammaloop_vectorize/crates/gammalooprs/src/integrands/process/mod.rs#L592:599)
<context ref="file:///Users/vjhirsch/Documents/Work/gammaloop_vectorize/crates/gammalooprs/src/integrands/process/mod.rs#L592:599">
    pub fn evaluate_samples_raw(
        &mut self,
        model: &Model,
        samples: &[Sample<F<f64>>],
        iter: usize,
        use_arb_prec: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<RawBatchEvaluationResult> {

</context> 

for our gammaloop code receives batch of samples, which will travel a long stack of our entire evaluation algorithm for calculating the corresponding integrand, events and observables for each of those samples.

Currently the vectorization is immediately undone, becaue we immediately do:

```
        for sample in samples {
        ...
```

The objective of this refactorization is to push through and retain the vectorization (over a batch size given through the settings, which can be different than the actual length of the samples vector given in input) all the way through the stacks, down to the call to symbolica's evaluator in:

[@evaluators.rs#L180:187](file:///Users/vjhirsch/Documents/Work/gammaloop_vectorize/crates/gammalooprs/src/integrands/process/evaluators.rs#L180:187)
<context ref="file:///Users/vjhirsch/Documents/Work/gammaloop_vectorize/crates/gammalooprs/src/integrands/process/evaluators.rs#L180:187">
    pub fn evaluate(&mut self, args: &[Complex<F<f64>>], out: &mut [Complex<F<f64>>]) {
        unsafe {
            self.0.evaluate(
                transmute::<&[Complex<F<f64>>], &[SymComplex<f64>]>(args),
                transmute::<&mut [Complex<F<f64>>], &mut [SymComplex<f64>]>(out),
            );
        }
    }

</context> 

which we sadly call individually for each sample.
However, since this evaluate call can internally benefit from SIMD, we'd like to call it with multiple inputs all at once.

The evaluation stack of gammaloop is mostly linear so it should not be too difficult.
However, one of the more tricky parts with branching is related to this evaluation stack

[@cross_section_integrand.rs#L666:1080](file:///Users/vjhirsch/Documents/Work/gammaloop_vectorize/crates/gammalooprs/src/integrands/process/cross_section_integrand.rs#L666:1080)
<context ref="file:///Users/vjhirsch/Documents/Work/gammaloop_vectorize/crates/gammalooprs/src/integrands/process/cross_section_integrand.rs#L666:1080">

    fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        model: &Model,
        settings: &RuntimeSettings,
        mut event_processing_runtime: Option<&mut EventProcessingRuntime>,
        _rotation: &Rotation,
        evaluation_metadata: &mut EvaluationMetaData,
        record_primary_timing: bool,
        channel_id: Option<(ChannelIndex, F<T>)>,
    ) -> Result<GraphEvaluationResult<T>> {
        let orientations =
            momentum_sample.orientations(&self.orientation_filter, &self.orientations);

        // let mut all_cut_result = Complex::new_re(momentum_sample.zero());
        let center = LoopMomenta::from_iter(vec![
            ThreeMomentum {
                px: momentum_sample.zero(),
                py: momentum_sample.zero(),
                pz: momentum_sample.zero(),
            };
            momentum_sample.loop_moms().len()
        ]);
        let masses = self.graph.get_real_mass_vector(model);
        let hel = settings.kinematics.externals.get_helicities();
        let mut cut_results = TiVec::<CutId, Complex<F<T>>>::new();
        let mut cut_threshold_counterterms = TiVec::<CutId, Complex<F<T>>>::new();
        let mut differential_result = GraphEvaluationResult::zero(momentum_sample.zero());
        let mut accepted_event_group = GenericEventGroup::default();

        let momentum_sample = if let Some((channel_id, _alpha)) = &channel_id {
            MomentumSample {
                sample: self.multi_channeling_setup.reinterpret_loop_momenta_impl(
                    *channel_id,
                    &momentum_sample.sample,
                    momentum_sample.sample.loop_mom_cache_id,
                ),
            }
        } else {
            momentum_sample.clone()
        };

        debug!("loop moms: {}", momentum_sample.loop_moms());

        for (raised_cut, raised_cut_group) in self.raised_data.raised_cut_groups.iter_enumerated() {
            let max_occurance = raised_cut_group.related_esurface_group.max_occurence;
            debug!("
 =====START EVALUTAION FOR CUT {}=====", raised_cut.0);
            let representative_esurface = &self.cut_esurface[raised_cut_group.cuts[0]];

            debug!("representative esurface: {:#?}", representative_esurface);

            let function = |t: &F<T>| {
                representative_esurface.compute_self_and_r_derivative(
                    t,
                    momentum_sample.loop_moms(),
                    &center,
                    momentum_sample.external_moms(),
                    &masses,
                    &self.graph.loop_momentum_basis,
                )
            };

            let (guess, _) = representative_esurface.get_radius_guess(
                momentum_sample.loop_moms(),
                momentum_sample.external_moms(),
                &self.graph.loop_momentum_basis,
            );

            let solution = newton_iteration_and_derivative(
                &guess,
                function,
                &F::from_f64(1.0),
                2000,
                &F::from_f64(settings.kinematics.e_cm),
            );

            debug!(
                "tolerance for newton solver: {}",
                F::from_f64(settings.kinematics.e_cm) * guess.epsilon()
            );

            debug!("solution: {:?}", solution);

            let mut event_timing = Duration::ZERO;
            let event = if settings.should_generate_events() {
                let start = Instant::now();
                let generated = self.generate_event_for_cut::<T>(
                    &model,
                    &solution,
                    &momentum_sample,
                    self.raised_data.raised_cut_groups[raised_cut].cuts[0],
                    &self.cuts[self.raised_data.raised_cut_groups[raised_cut].cuts[0]],
                )?;
                event_timing += start.elapsed();
                Some(generated)
            } else {
                None
            };

            let mut accepted_event = event;
            let mut selectors_pass = true;
            if let Some(evt) = accepted_event.as_mut() {
                differential_result.generated_event_count += 1;
                let start = Instant::now();
                if let Some(runtime) = event_processing_runtime.as_deref_mut() {
                    if runtime.has_selectors() || runtime.has_observables() {
                        selectors_pass = runtime.process_event(evt);
                    }
                }
                event_timing += start.elapsed();
                if !selectors_pass {
                    accepted_event = None;
                }
            }
            differential_result.event_processing_time += event_timing;

            if !selectors_pass {
                let zero = Complex::new_re(momentum_sample.zero());
                for _ in 1..=max_occurance {
                    cut_threshold_counterterms.push(zero.clone());
                    cut_results.push(zero.clone());
                }
                continue;
            }

            let mut bare_cut_total = Complex::new_re(momentum_sample.zero());
            let mut threshold_counterterm_weights = Vec::with_capacity(max_occurance);
            for num_esurfaces in 1..=max_occurance {
                let dual_shape = if num_esurfaces > 1 {
                    Some(HyperDual::<F<T>>::new(
                        self.raised_data.dual_shapes[num_esurfaces - 2].clone(),
                    ))
                } else {
                    None
                };

                let (tstar, h_function, esurface_derivatives, rescaled_momenta) =
                    if let Some(dual_shape) = dual_shape {
                        let dual_t_for_integrand =
                            dual_shape.variable(0, solution.solution.clone());
                        let dual_h_function =
                            h_dual(&dual_t_for_integrand, None, None, &settings.lu_h_function);
                        let dual_momenta_for_integrand = momentum_sample
                            .loop_moms()
                            .rescale_with_hyper_dual(&dual_t_for_integrand, None);

                        let dual_shape_for_esurface =
                            HyperDual::<F<T>>::new(shape_for_t_derivatives(num_esurfaces));

                        let dual_t_for_esurface =
                            dual_shape_for_esurface.variable(0, solution.solution.clone());
                        let dual_momenta_for_esurface = momentum_sample
                            .loop_moms()
                            .rescale_with_hyper_dual(&dual_t_for_esurface, None);

                        let dual_externals = momentum_sample
                            .external_moms()
                            .iter()
                            .map(|mom| FourMomentum {
                                temporal: Energy {
                                    value: new_constant(&dual_t_for_esurface, &mom.temporal.value),
                                },
                                spatial: ThreeMomentum {
                                    px: new_constant(&dual_t_for_esurface, &mom.spatial.px),
                                    py: new_constant(&dual_t_for_esurface, &mom.spatial.py),
                                    pz: new_constant(&dual_t_for_esurface, &mom.spatial.pz),
                                },
                            })
                            .collect();

                        let dual_e_surface = representative_esurface.compute_from_dual_momenta(
                            &self.graph.loop_momentum_basis,
                            &masses,
                            &dual_momenta_for_esurface,
                            &dual_externals,
                        );

                        let mut momentum_sample_with_duals = momentum_sample.clone();
                        momentum_sample_with_duals.sample.dual_loop_moms =
                            Some(dual_momenta_for_integrand);

                        (
                            DualOrNot::Dual(dual_t_for_integrand),
                            DualOrNot::Dual(dual_h_function),
                            DualOrNot::Dual(dual_e_surface),
                            momentum_sample_with_duals,
                        )
                    } else {
                        let h_function = h(&solution.solution, None, None, &settings.lu_h_function);
                        let rescaled_momenta = momentum_sample
                            .rescaled_loop_momenta(&solution.solution, Subspace::None);

                        (
                            DualOrNot::NonDual(solution.solution.clone()),
                            DualOrNot::NonDual(h_function),
                            DualOrNot::NonDual(solution.derivative_at_solution.clone()),
                            rescaled_momenta,
                        )
                    };

                debug!("tstar: {}", tstar);
                debug!("h(tstar): {}", h_function);
                debug!("esurface derivative at tstar: {}", esurface_derivatives);

                let lu_params = LUParams { h_function, tstar };

                let prefactor =
                    Complex::new_re(if let Some((_channel_index, _alpha)) = &channel_id {
                        if matches!(lu_params.tstar, DualOrNot::Dual(_)) {
                            panic!("multi channeling with duals not supported yet");
                        }

                        self.multi_channeling_setup.compute_prefactor_impl(
                            *_channel_index,
                            &rescaled_momenta,
                            model,
                            _alpha,
                        )
                    } else {
                        F::from_f64(1.0)
                    });

                let params = T::get_parameters(
                    &mut self.param_builder,
                    (settings.general.enable_cache, settings.general.debug_cache),
                    &self.graph,
                    &rescaled_momenta,
                    hel,
                    &settings.additional_params(),
                    None,
                    None,
                    Some(&lu_params),
                );

                let result = self.integrand[raised_cut][num_esurfaces - 1]
                    .evaluate(
                        params,
                        orientations,
                        settings,
                        evaluation_metadata,
                        record_primary_timing,
                    )
                    .expect("evaluation failed")
                    .pop()
                    .unwrap();

                let ct_result = if settings.subtraction.disable_threshold_subtraction {
                    Complex::new_re(momentum_sample.zero())
                } else {
                    todo!();
                };

                debug!("pass 1 result {}", result);

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

                match esurface_derivatives {
                    DualOrNot::Dual(dual_e_surface) => {
                        extract_t_derivatives(dual_e_surface)[1..]
                            .iter()
                            .for_each(|v| {
                                params_for_pass_two.push(Complex::new_re(v.clone()));
                            });
                    }
                    DualOrNot::NonDual(non_dual_e_surface) => {
                        params_for_pass_two.push(Complex::new_re(non_dual_e_surface));
                    }
                }

                let pass_two_evaluator =
                    &mut self.raised_data.pass_two_evaluators[num_esurfaces - 1];

                let pass_two_result = evaluate_evaluator_single(
                    pass_two_evaluator,
                    &params_for_pass_two,
                    evaluation_metadata,
                    record_primary_timing,
                );

                debug!("pass_two_result: {:+16e}", pass_two_result);

                cut_threshold_counterterms.push(ct_result.clone());
                //debug!("param builder for cut {}: 
{}", cut, self.param_builder);

                let bare_contribution = pass_two_result * prefactor; //   * Complex::new_im(-momentum_sample.one()).pow(subset.len() as u64),
                bare_cut_total += bare_contribution.clone();
                threshold_counterterm_weights.push(ct_result);
                cut_results.push(bare_contribution);
            }

            if let Some(mut event) = accepted_event {
                let threshold_counterterm_total = threshold_counterterm_weights
                    .iter()
                    .fold(Complex::new_re(momentum_sample.zero()), |acc, value| {
                        acc + value.clone()
                    });
                event.weight = bare_cut_total.clone() + threshold_counterterm_total;

                if settings.general.store_additional_weights_in_event {
                    event
                        .additional_weights
                        .weights
                        .insert(AdditionalWeightKey::Original, bare_cut_total);
                    for (subset_index, threshold_counterterm) in
                        threshold_counterterm_weights.into_iter().enumerate()
                    {
                        event.additional_weights.weights.insert(
                            AdditionalWeightKey::ThresholdCounterterm { subset_index },
                            threshold_counterterm,
                        );
                    }
                }

                differential_result.accepted_event_count += 1;
                if settings.should_buffer_generated_events() {
                    accepted_event_group.push(event);
                }
            }
        }

        let mut all_cut_result = Complex::new_re(momentum_sample.zero());
        for ((cut_id, result), ct_result) in cut_results
            .iter_enumerated()
            .zip(cut_threshold_counterterms.iter())
        {
            debug!("cut {} results:", cut_id);
            debug!("bare result: {:+16e}", result);
            debug!("threshold counterterm: {:+16e}", ct_result);
            all_cut_result += result + ct_result;
        }

        let flux_factor = if settings.general.disable_flux_factor {
            F::from_f64(1.0)
        } else {
            match momentum_sample.external_moms().len() {
                1 => {
                    momentum_sample.one()
                        / (F::from_f64(2.0)
                            * &momentum_sample
                                .external_moms()
                                .first()
                                .as_ref()
                                .unwrap()
                                .temporal
                                .value)
                }
                2 => {
                    let mom_1 = &(momentum_sample.external_moms()[ExternalIndex::from(0)]);
                    let mom_2 = &(momentum_sample.external_moms()[ExternalIndex::from(1)]);
                    let mass_factor = self
                        .graph
                        .initial_state_cut
                        .iter_edges(&self.graph)
                        .map(|(_, e)| e.data.mass.value(model, &self.param_builder).unwrap())
                        .fold(Complex::new_re(momentum_sample.one()), |acc, mass| {
                            acc * &mass * &mass
                        })
                        .re;

                    let f = F::from_f64(4.0) * (mom_1.dot(mom_2).square() - mass_factor).sqrt();

                    momentum_sample.one() / f
                        * if settings.general.use_picobarns {
                            F::from_ff64(PICOBARN_CONVERSION)
                        } else {
                            F::from_f64(1.0)
                        }
                }
                _ => unimplemented!(
                    "Flux factor for more than 3 or more incoming particles not implemented yet"
                ),
            }
        };

        let final_result = all_cut_result * flux_factor.clone();

        if settings.should_buffer_generated_events() {
            let flux_factor = Complex::new_re(flux_factor);
            for event in accepted_event_group.iter_mut() {
                event.weight *= flux_factor.clone();
                if !event.additional_weights.weights.is_empty() {
                    event.additional_weights.weights.insert(
                        AdditionalWeightKey::FullMultiplicativeFactor,
                        flux_factor.clone(),
                    );
                }
            }

            if !accepted_event_group.is_empty() {
                differential_result.event_groups.push(accepted_event_group);
            }
        }

        debug!(
            "{}",
            format!(
                "final result for graph: {}, {:+16e}",
                self.graph.name, final_result
            )
            .red()
        );

        differential_result.integrand_result = final_result;
        Ok(differential_result)
    }

</context> 

where multiple evaluators are being called per sample for each threshold CTs and using centers that were computed using a convex solver (that won't be vectorized), and wheter or not a threshold CT is active and have its evaluator called depends on the sample.
So for these parts one must store all paths in a buffer until it reaches the need for an evaluator to be called, and when all samples in `evaluate_samples` have been processed, all the inputs must be joined so that the vectorized calls to symbolica's evaluators can be issued.

There are many ways of achieving this, but ideally I'd like to avoid the use of any global thread-safe container and instead find a solution that involves extra buffered inputs being dragged along the vectorized chain perhaps.
Well, to be fair I'm not sure what would be the best approach but I would like you to deeply study the structure of the evaluation stack in gammaloop and come up with an elegant solution for implementing this with a clear identification of all the bottlenecks and your proposed solution for each.
(no implementation needed at this point, just planning for now).

Throughout the implementation and our discussions I would like your to:

a) Write your plan, and always keep it up to date in the top-level repository file "vectorized_gammaloop.md"

b) Record all our interactions (my prompts, copied verbatim, with *no edit at all from you*) as well as your responses (only what you emit back to me, not the result of intermediate reads you do, or code output your run, or tests outputs or search results, just the actual responses you send back to me here to read) in a top-level repository file "codex_human_interaction_log.md"

Ask questions if some additional input from me is needed to help you craft and efficient, elegant and accurate plan.
````

### Codex

````
I studied the stack from `crates/gammalooprs/src/integrands/process/mod.rs:592` through stability/rotation dispatch, graph-term evaluation, `EvaluatorStack`, `ParamBuilder`, `CrossSectionIntegrand`, and both threshold-counterterm stacks. The cleanest design is a local wavefront batch engine: keep all non-Symbolica work scalar per sample, but whenever a sample path reaches an evaluator boundary, materialize an owned parameter buffer, enqueue it into a batch keyed by evaluator identity, flush via Symbolica's batch API, scatter results back, and resume the sample path. That avoids any global container and lets the batch state travel through the call chain.

- Main bottlenecks:
  - `crates/gammalooprs/src/integrands/process/mod.rs:592` immediately scalarizes the batch.
  - `crates/gammalooprs/src/integrands/process/mod.rs:1704` and `crates/gammalooprs/src/integrands/process/mod.rs:2536` keep stability and rotation dispatch sample-local.
  - `crates/gammalooprs/src/integrands/process/param_builder.rs:702` returns borrowed shared parameter storage, so queued requests would overwrite each other.
  - `crates/gammalooprs/src/integrands/process/evaluators.rs:502` still loops orientations around immediate evaluator calls.
  - `crates/gammalooprs/src/integrands/process/cross_section_integrand.rs:667` has a true continuation point: pass 1 result is needed to build pass 2 inputs, with selector/event early exits in between.
  - `crates/gammalooprs/src/subtraction/lu_counterterm.rs:222` and `crates/gammalooprs/src/subtraction/amplitude_counterterm.rs:132` are the branching hotspots for threshold subtraction.

- Proposed architecture:
  - Add a local `BatchEvaluationSession` / wavefront driver under `evaluate_samples_raw`.
  - Introduce owned/materialized evaluator inputs so batch requests never borrow `self.param_builder.values`.
  - Add `GenericEvaluator::evaluate_batch_f64` / `evaluate_batch_single_f64` plus `EvaluatorStack::evaluate_batch`, grouping by evaluator id and output shape; use Symbolica's existing batch evaluator backend for compiled f64, with scalar fallback otherwise.
  - Make batched variants of `evaluate_from_source`, `evaluate_stability_level`, `evaluate_all_rotations`, `evaluate_single`, and graph-term `evaluate_batch`.
  - Implement `CrossSectionGraphTerm::evaluate_batch` as a staged loop over `(raised_cut, num_esurfaces)`, with per-sample cut-state continuations storing pending events, prefactors, derivatives, accumulators, and CT subjobs.
  - Then port amplitude onto the same substrate; it should be much simpler once the shared batch machinery exists.

- Recommended scope:
  - Keep parameterization, Newton solves, overlap detection, and convex/overlap solvers scalar for now; they just feed the batch queues.
  - Target double-precision batching first; keep `f128` / `ArbPrec` on the same orchestration shell but with scalar evaluator execution.
  - Preserve current event-group semantics by keeping event objects in sample-local state until the final weight is known.
  - Apportion batch timing back onto sample metadata instead of charging the whole batch wall-time to every participating sample.

I wrote the detailed plan to `vectorized_gammaloop.md` and logged this exchange in `codex_human_interaction_log.md`.

Two decisions would help before implementation:
1. Is phase 1 allowed to optimize only the `f64` raw path, while `f128` / `ArbPrec` stay scalar behind the same abstraction?
2. Do you want amplitude included in phase 1, or should I treat `CrossSectionIntegrand` as the primary target and port amplitude right after the shared substrate is stable?
3. Do you want me to explicitly target Symbolica's SIMD compiled evaluator as the optimized backend when available, with scalar compiled/eager fallback, or keep phase 1 backend-agnostic and wire SIMD in right after the control-flow refactor?
````


## 2026-03-20 (follow-up)

### Human

````
1. Yes, we only ever want to vectorize the f64 path, f128/ArbPrec can remain scalar. It is actually good if these paths remain explicitly scalar, so that we can always verify against them that the batching of f64 calls did not introduce side-effects. When not batching it (or equivalently supplying a batch of one) the code should have very close to the same performance as of now, so beware of allocation, the code you're touching is very runtime sensitive!

2. Do amplitudes as well.

3. For now just make sure to call batch evaluators, but do not try to make sure that the batch evaluators leverage SIMD. We have ideas and other plans for this; so long as the inputs of the evaluators are batched at their calling site, we are good.

A couple of design choices.

* Ideally the presence of the wavefront batch engine should disrupt as little as possible the understanding of the underlying mostly linear evaluation stack.
Because the vectorization only benefits the lowest-level calls in symbolica's evaluator, you should still retain and benefit from any early termination, just pooling evaluator calls as you go.

* The idea is that the batch size the top level "evaluate_samples" will be called for may well be large, (e.g. 10000), but then you'd re-batch it to smaller batches of user-defined size GAMMALOOP_BATCH_SIZE (so as to not overload / oversize the wavefront batching engine buffers). Then if GAMMALOOP_BATCH_SIZE is say about 100, then there are good chances that most evaluators will land with 16 or more inputs, meaning that they'd fully benefit from avx512.

Comment on the validity of this logic, along side with a revision of your plan and another deep pass at the code to fully grasp the subtelties that may arise in the implementation.

I'll then do another review of theplan
````

### Codex

````
I agree with your batching logic overall. Re-batching a very large incoming `evaluate_samples_raw(...)` slice into smaller internal f64 sub-batches is the right shape: it bounds the live frontier, keeps the batching layer local, and still gives the hot evaluators many chances to accumulate enough requests for future SIMD backends.

The two main caveats are:
- `GAMMALOOP_BATCH_SIZE` only bounds the number of active sample states, not the exact size of every evaluator queue. Some queues can exceed it because of orientation expansion and left/right threshold cartesian products.
- To keep batch-of-one performance close to current, we should preserve an exact scalar fast path and avoid per-request heap allocation in the f64 path.

I did another deep pass through the code and the main new subtleties I would call out are:
- `process_evaluation_result_runtime(...)` only runs after the per-sample `EvaluationResult` is finalized (`crates/gammalooprs/src/integrands/process/mod.rs:623`, `crates/gammalooprs/src/integrands/process/mod.rs:832`), so we can keep observables/result finalization exactly as today and only batch the internal evaluator work.
- The current f64 `ParamBuilder` already contains the right low-allocation template (`crates/gammalooprs/src/integrands/process/param_builder.rs:715`, `crates/gammalooprs/src/integrands/process/param_builder.rs:760`). In the batched path, we should reuse that by copying the template into preallocated queue spans and patching only the dynamic fields, not by allocating one `Vec` per request.
- `SingleOrAllOrientations` is borrowed (`crates/gammalooprs/src/momentum/sample.rs:765`), so we should not store it across batching boundaries. Store only owned orientation state in the sample/cut structs and recreate the borrowed view while filling a queue.
- Cross-section LU threshold subtraction is architecturally important, but the current active call site is still `todo!()` when subtraction is enabled (`crates/gammalooprs/src/integrands/process/cross_section_integrand.rs:913`). So phase 1 should batch the currently active pass 1 / pass 2 chain first, while designing the queue API so LU threshold jobs can plug in later without another refactor.
- The pinned Symbolica dependency already exposes the batch abstraction we need for now, even without any SIMD work. So the GammaLoop-side job is to make the call sites batch-aware, not to change the backend yet.

Because you want the presence of the batching layer to disrupt understanding as little as possible, I am revising the plan slightly: I no longer think the first implementation should look like a generic opaque wavefront scheduler. The better fit is to keep the current mostly linear control flow and thread explicit stage-local batch collectors through it. In other words: the sample states still march through the same stages as today, and only the evaluator boundaries become batched.

So the revised plan is:
- Keep an exact scalar fast path: if the effective internal sub-batch size is `<= 1`, dispatch to the existing scalar code unchanged.
- In `evaluate_samples_raw(...)`, split the large incoming slice into internal `GAMMALOOP_BATCH_SIZE` f64 sub-batches, preserve output order, and finalize/process results in the same order as today.
- Add f64 batch helpers on `GenericEvaluator` / `EvaluatorStack`, but use Symbolica's existing batch evaluator API rather than any SIMD-specific backend work.
- For the f64 batched path, use reusable queue arenas and write parameters into spans; do not allocate a fresh request vector for every evaluator call.
- Keep `f128` and `ArbPrec` entirely scalar and hand escalated jobs back to the existing scalar code.
- Implement `CrossSectionGraphTerm::evaluate_batch_f64(...)` so it still reads like the current function: build per-sample cut state, queue pass 1 requests, flush, queue pass 2 requests, flush, finalize events/weights.
- Then implement `AmplitudeGraphTerm::evaluate_batch_f64(...)` as the second consumer of exactly the same substrate.

I updated `vectorized_gammaloop.md` with that revised plan and appended this exchange to `codex_human_interaction_log.md`.

I think this is now aligned with your three decisions:
1. f64 only is batched.
2. amplitude is included.
3. the caller side becomes batch-aware now, while the backend remains batch-evaluator based and SIMD-agnostic for the moment.
````

## User
Ok, a couple more things:

* GAMMALOOP_BATCH_SIZE should be runtimesettings of the integrand, in the EvaluatorSettings struct. don't forget to properly decorate with serde if_default using the escape hatch "SHOW_DEFAULT".

* Add a monitor surfaced back through EvaluationMetaData of the average size of the batches when they reach the evaluator call sites.

* For timing measurements, make sure that things are normalized exactly as to keep their interpretation in the scalar case. i.e. the integrator timings should appear 4 times faster per sample if 4x SIMD were to be maximally efficient.

Revise your plan according to the above and do one last deep dive in the code to identify any culprit or weakness of your plan, and give me its update here!

## Assistant
I updated `vectorized_gammaloop.md`.

- Your batching logic still holds. Re-batching large incoming slices into smaller internal f64 sub-batches is the right control knob. The only caveat is that `GAMMALOOP_BATCH_SIZE` caps active sample states, not every concrete evaluator queue; hot queues can still grow beyond it through orientation and threshold branching.
- The main design correction is settings placement. The current `processes::EvaluatorSettings` is generation-scoped, while `evaluator_method` is runtime-scoped in `GeneralSettings`. Since `batch_size` is a runtime evaluation knob, my revised plan is to split generation and runtime evaluator settings cleanly, let runtime-side `EvaluatorSettings` own `method` and `batch_size`, and decorate `batch_size` with the existing `SHOWDEFAULTS` serde helper path.
- The batch monitor should live at the concrete `GenericEvaluator` flush boundary, not at the sample or graph-term level. I am planning a request-weighted `average_evaluator_batch_size` in `EvaluationMetaData`, then threading the same numerator/denominator through `StatisticsCounter` and `IntegrationStatisticsSnapshot` so it remains visible in aggregated integration stats.
- The largest implementation trap is timing. The current f64 path measures per-sample wall clock in `evaluate_from_source(...)` / `evaluate_single(...)`, and that becomes wrong once samples wait on shared flushes. The revised plan is to give each active sample a timing ledger: direct scalar work is charged directly, shared evaluator flush wall time is apportioned by logical evaluator-request count, and `integrand_evaluation_time`, `evaluator_evaluation_time`, and `total_timing` are synthesized from that ledger. That preserves scalar interpretation, so per-sample timings drop when batching improves throughput.
- The execution model stays the same otherwise: exact scalar fast path for effective batch size `<= 1`, explicit stage-local queues with preallocated parameter arenas, cross section as the first consumer, amplitude on the same substrate right after.

Main remaining weaknesses in the plan:
- The settings split is the main structural change; if we keep one mixed evaluator-settings struct, the config model will become harder to reason about than the batching logic.
- The batch-size monitor will be misleading if it is flush-weighted or sample-weighted instead of request-weighted.
- Cross-section threshold subtraction is still partially inactive at `crates/gammalooprs/src/integrands/process/cross_section_integrand.rs:913`, so the queue API has to stay pluggable there rather than assuming that path is already live.

I also appended this exchange to `codex_human_interaction_log.md`.

## User
Perfect, one last thing I obviously mispoke and `gammaloop_batch_size` should be part of the `runtime` settings of the integrator, not the generation-side settings of course, and it's ok that the size of the queues could grow beyond that. This only controls the re-batching at the entry point.
Update your plan one last time with the above and ask any question if some remains

## Assistant
I updated `vectorized_gammaloop.md` again.

- `gammaloop_batch_size` is now explicitly planned under runtime `IntegratorSettings`, not under any evaluator settings.
- I also tightened the wording so it is treated only as the entry-point rechunking limit. Concrete evaluator queues are allowed to grow beyond it, and the deeper batching logic should not reuse it as a queue cap.
- The request-weighted batch monitor and request-weighted timing normalization stay unchanged; those parts were already consistent with this correction.
- This simplifies the settings side of the plan: there is no need to move batch-size concerns into evaluator settings at all.

One remaining question:
- Do you want the same `runtime.integrator.gammaloop_batch_size` rechunking policy to apply to `evaluate_momentum_configurations_raw(...)` as well, or only to `evaluate_samples_raw(...)`?

I also appended this exchange to `codex_human_interaction_log.md`.

## 2026-03-20 (resume)

### Human

````
Are you still here? If so, say "yes" and continue
````

### Codex

````
yes

Implemented the batched f64 path and added standalone regression coverage.

- `runtime.integrator.gammaloop_batch_size` now rechunks both raw batch entry points.
- Batched f64 evaluator calls now flow through amplitude and cross-section stacks to Symbolica batch evaluators.
- `EvaluationMetaData` now carries request-weighted `average_evaluator_batch_size`; integration stats aggregate it.
- Batched timing is normalized per logical evaluator request so scalar interpretation is preserved.
- Added `test_batched_evaluation` in `tests/tests/test_runs.rs` covering standalone `gg_hhh` and `e+ e- > d d~ g`, including event/additional-weight surfacing and batch-monitor assertions.
- Updated `vectorized_gammaloop.md`.

Validation:
- `cargo check`
- `EXTRA_MACOS_LIBS_FOR_GNU_GCC=T cargo test -p gammaloop-integration-tests --test test_runs test_batched_evaluation -- --nocapture`

One caveat remains unchanged from the preexisting code: cross-section threshold subtraction is still on the existing `todo!()` path when enabled, so the new batching layer does not change that behavior.
````


## 2026-03-20 (backfilled after compaction)

### Human

````
1. Yes, it should apply to all batched inputs from our API which supports vectorizing downstream.

2. For testing, you should implement in `test_runs.rs` a test function `test_batched_evaluation` which uses a process like gg_hhh (with physical kinematics, so threshold counterterms are enabled) like obtained from

[@gg_hhh_1L.toml#L1:188](file:///Users/vjhirsch/Documents/Work/gammaloop_vectorize/examples/cli/gg_hhh/1L/gg_hhh_1L.toml#L1:188)
<context ref="file:///Users/vjhirsch/Documents/Work/gammaloop_vectorize/examples/cli/gg_hhh/1L/gg_hhh_1L.toml#L1:188">
commands = []

[cli_settings]
[cli_settings.state]
folder = "./examples/cli/gg_hhh/1L/state"
name = "gg_hhh_1l"

[cli_settings.global.generation.compile]
inline_asm = true
optimization_level = "O1"
fast_math = true
unsafe_math = true
custom = []

[cli_settings.global.generation.tropical_subgraph_table]
panic_on_fail = true
disable_tropical_generation = false

[cli_settings.global.generation.evaluator]
iterative_orientation_optimization = false
store_atom = false
compile = true
summed = false
summed_function_map = true

# Disable thresholds with
# [cli_settings.global.generation.threshold_subtraction]
# enable_thresholds = false
# [default_runtime_settings.subtraction]
# disable_threshold_subtraction = true

# Otherwise leave them enabled with
[cli_settings.global.generation.threshold_subtraction]
enable_thresholds = true
# The option below specifies the generation of threshold for the kinematics supplied
check_esurface_at_generation = true

[cli_settings.global.n_cores]
feyngen = 10
generate = 10
compile = 10
integrate = 10

[[command_blocks]]
name = "generate"
commands = [
  "import model sm-default.json",
  "set global kv global.n_cores.feyngen=10",

  "reset processes",
  """generate amp g g > h h h / u d c s b QED==3 [{1}]
           --only-diagrams
           --numerator-grouping only_detect_zeroes
           --select-graphs GL15
           --loop-momentum-bases GL15=8
           --global-prefactor-projector '1𝑖 * gammalooprs::ϵ(0,spenso::mink(4,gammalooprs::hedge(0)))
                                           * gammalooprs::ϵ(1,spenso::mink(4,gammalooprs::hedge(1)))
                                           * (1/8)*spenso::g(spenso::coad(8,gammalooprs::hedge(0)),spenso::coad(8,gammalooprs::hedge(1)))'
           -p gg_hhh
           -i 1L
    """,
  "generate",
  "save dot",

  "save state -o",
]

[[command_blocks]]
name = "integrate_euclidean"
commands = [
  "set model MT=1000.0",
  "set model WT=0.0",
  "set model ymt=1000.0",

  "set process -p gg_hhh -i 1L defaults",

  # Disable stability checks
  """set process -p gg_hhh -i 1L string '
    [stability]
    rotation_axis = [ ]
    '""",

  # Add stability checks
  # """set process -p gg_hhh -i 1L string '
  # [[stability.rotation_axis]]
  # type = "z"
  # [[stability.rotation_axis]]
  # type = "x"
  # [[stability.rotation_axis]]
  # type = "y"
  # '""",

  "set process -p gg_hhh -i 1L kv integrator.n_start=100000 integrator.n_max=1000000",
  # Target for "euclidean" point confirmed with MadLoop
  """integrate --process gg_hhh -i 1L --n-cores 10
        --target 6.56089133881216768e-4 -4.17078968913725420e-6
        --no-show-max-weights
        --restart
    """,
]

[[command_blocks]]
name = "integrate_physical"
commands = [
  "set process -p gg_hhh -i 1L defaults",

  # Disable stability checks
  """set process -p gg_hhh -i 1L string '
    [stability]
    rotation_axis = [ ]
    '""",

  # Target for the physical top mass of 173 GeV is also confirmed with MadLoop
  "set model MT=173.0",
  "set model WT=0.0",
  "set model ymt=173.0",
  "set process -p gg_hhh -i 1L kv integrator.n_start=100000 integrator.n_max=10000000",
  """integrate --process gg_hhh -i 1L --n-cores 10
        --target 1.89203604685291554e-6 -5.46603099163412881e-6
        --workspace-path ./GL_OUTPUT/gg_hhh/1L/state/integration_workspace
        --restart
    """,
]

[default_runtime_settings.general]
evaluator_method = "SummedFunctionMap"
enable_cache = false
debug_cache = false

[default_runtime_settings.kinematics.externals]
type = "constant"

[default_runtime_settings.kinematics.externals.data]
momenta = [
  [500.0, 0.0, 0.0, 500.0],
  [500.0, 0.0, 0.0, -500.0],
  [
    0.4385555662246945e3,
    0.1553322001835378e3,
    0.3480160396513587e3,
    -0.1773773615718412e3
  ],
  [
    0.3563696374921922e3,
    -0.1680238900851100e2,
    -0.3187291102436005e3,
    0.9748719163688098e2
  ],
  "dependent"
]
helicities = [1, 1, 0, 0, 0]

# Default "naive" parameterization
# [default_runtime_settings.sampling]
# type = "default"
# mode = "spherical"
# mapping = "linear"
# b = 10.0

## a) Tropical sampling run.
# [default_runtime_settings.sampling]
# type = "discrete_graph_sampling"
# [default_runtime_settings.sampling.sampling_type]
# subtype = "tropical"

## b) MC over LMB with explicit sum per sample point.
# [default_runtime_settings.sampling]
# type = "discrete_graph_sampling"
# [default_runtime_settings.sampling.sampling_type]
# subtype = "multi_channeling"

## c) MC over LMB with importance sampling of channels
[default_runtime_settings.sampling]
type = "discrete_graph_sampling"

[default_runtime_settings.sampling.sampling_type]
subtype = "discrete_multi_channeling"

## d) MC over LMB and single orientation per, both importance-sampled
# [default_runtime_settings.sampling]
# type = "discrete_graph_sampling"
# sample_orientations = true
# [default_runtime_settings.sampling.sampling_type]
# subtype = "discrete_multi_channeling"

[default_runtime_settings.subtraction]
disable_threshold_subtraction = false

</context> 

and also from `e+ e- > d d~ g` as in

[@run.toml#L1:131](file:///Users/vjhirsch/Documents/Work/gammaloop_vectorize/examples/api/rust/epem_a_ddxg_xs_LO/run.toml#L1:131)
<context ref="file:///Users/vjhirsch/Documents/Work/gammaloop_vectorize/examples/api/rust/epem_a_ddxg_xs_LO/run.toml#L1:131">
#:schema https://raw.githubusercontent.com/alphal00p/gammaloop/refs/heads/HEAD/assets/schemas/runhistory.json

# Differential LU example for the real two-graph
#   e+ e- > d d~ g | e- a d g QED^2==4 [{{2}} QCD=0]
# process used in the end-to-end differential tests.
#
# This card sets up:
# - event generation
# - additional partial weights in the generated events
# - one shared jet-pt quantity
# - one shared jet-count quantity
# - one shared down-quark energy quantity
# - one leading-jet selector
# - one leading-jet histogram observable
# - one jet-count histogram observable
# - one down/anti-down energy histogram observable
# - reusable display examples for quantities, observables, and selectors
# - a local state under ./state

commands = [
  "import model sm-default",
  "run generate",
  "run set_differential_setup",
  "save state -o",
  "save dot"
]

[[command_blocks]]
name = "generate"
commands = [
  "reset processes",
  "generate e+ e- > d d~ g | e- a d g QED^2==4 [{{2}} QCD=0] --numerator-grouping group_identical_graphs_up_to_sign --clear-existing-processes --only-diagrams",
  "generate",
]

[[command_blocks]]
name = "set_differential_setup"
commands = [
  """set process string '
[quantities.leading_jet_pt]
type = "jet_pt"
algorithm = "anti_kt"
dR = 0.4
min_jpt = 0.0

[quantities.jet_count]
type = "jet_count"
algorithm = "anti_kt"
dR = 0.4
min_jpt = 0.0

[quantities.down_energy]
type = "particle_scalar"
pdgs = [1, -1]
quantity = "E"

[selectors.leading_jet_pt_cut]
quantity = "leading_jet_pt"
selector = "value_range"
entry_selection = "leading_only"
min = 0.0

[observables.leading_jet_pt_hist]
quantity = "leading_jet_pt"
entry_selection = "leading_only"
phase = "real"
x_min = 0.0
x_max = 1000.0
n_bins = 8

[observables.down_energy_hist]
quantity = "down_energy"
entry_selection = "all"
phase = "real"
x_min = 0.0
x_max = 1000.0
n_bins = 8

[observables.jet_count_hist]
quantity = "jet_count"
entry_selection = "all"
phase = "real"
x_min = 0.0
x_max = 6.0
n_bins = 6
'"""
]

[[command_blocks]]
name = "display_named_settings_examples"
commands = [
  "display quantities -p #0",
  "display quantities -p #0 leading_jet_pt",
  "display observables -p #0",
  "display observables -p #0 leading_jet_pt_hist",
  "display selectors -p #0",
  "display selectors -p #0 leading_jet_pt_cut",
]

[cli_settings]
override_state = false
try_strings = true

[cli_settings.state]
folder = "./state"
name = "CLI_demo"

[cli_settings.global]
logfile_directive = "info,symbolica=off"
display_directive = "info,symbolica=off"

[default_runtime_settings.general]
generate_events = true
store_additional_weights_in_event = true

[default_runtime_settings.kinematics.externals]
type = "constant"

[default_runtime_settings.kinematics.externals.data]
momenta = [
  [32.0, 0.0, 0.0, 32.0],
  [32.0, 0.0, 0.0, -32.0],
]
helicities = [1, 1]

[default_runtime_settings.subtraction]
disable_threshold_subtraction = true

[default_runtime_settings.integrator.observables_output]
format = "json"

</context> 

which can test events surfacing in the presence of vectorization)

(but make sure the test is standalone and does not rely on those toml cards).

The test is then simply to compare that the outputs between the batched eval and the single eval are in agreement (bar timing).

You can go on with the implementation
````

### Human

````
I commited by mistake my changes to main, can you undo the commit and place it instead in a new branch called vectorize
````

### Codex

````
Done.

- Created `vectorize` at commit `b9d4daea`.
- Moved `main` back to `250e4f7c`.
- You are now on `vectorize`.
- Working tree is clean.
````

### Human

````
Create a PR to main, with the body of the PR taken from your vecorize.md plan but with an additional easy-to-read header file that reads nice when displayed as a PR body and which describes well your approach.
````

### Codex

````
Created PR `#49` to `main`: https://github.com/alphal00p/gammaloop/pull/49

- Title: `Vectorize f64 evaluator batching through the process stack`
- Body: starts with a short PR-oriented overview, then includes the full contents of `vectorized_gammaloop.md`
- Branch: `vectorize`

I also updated `codex_human_interaction_log.md` on `vectorize` and pushed that to the PR.
````


### Human

````
Pull github comments for this PR using gh and evaluate their relevance and address them if warranted.
You can ignore the comment regarding `codex_human_interaction_log.md`, stating that this will be removed post merge, same as for `vectorize_gammaloop.md`.
````

### Codex

````
Pulled the PR comments with `gh` and reviewed them.

- Ignored per your instruction: the `codex_human_interaction_log.md` comment.
- Relevant and fixed:
  - batched `Summed` / `SummedFunctionMap` evaluation now sets `override_if` for `OwnedOrientations::All`, matching the scalar path
  - `evaluate_batch_single_f64(...)` now enforces single-output use, requires consistent input arity, and sizes compiled output buffers correctly
  - compiled batch decoding now works from slices, so it no longer allocates one `Vec` per batch element
  - batched timing now includes batch orchestration, rotation preparation, stability checks, and final result assembly
  - integrand/evaluator timing and average batch-size accounting now stay on the primary stability rotation only; secondary rotations contribute only to `total_timing`

- I updated `vectorized_gammaloop.md` with the review follow-up.
- I pushed the fixes to `vectorize` at `686dfbab`, so PR `#49` is updated.
- I also posted a PR comment summarizing the addressed review points: https://github.com/alphal00p/gammaloop/pull/49#issuecomment-4101384496

Validation:
- `cargo check`
- `EXTRA_MACOS_LIBS_FOR_GNU_GCC=T cargo test -p gammaloop-integration-tests --test test_runs test_batched_evaluation -- --nocapture`
````
