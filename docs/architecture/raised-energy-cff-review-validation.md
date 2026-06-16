> Historical record for the revisions and dates stated below. Its test counts,
> API descriptions and continuation instructions do not certify or govern the
> reconstructed stack. See the [stack review](raised-energy-cff-stack-review.md)
> and [fresh validation record](raised-energy-cff-stack-review-validation.md),
> together with [current architecture](architecture-current.md) and
> [CONTRIBUTING.md](../../CONTRIBUTING.md). The original body is preserved.

# Raised-energy CFF review: executed validation

Evidence for `main` (`395610143`) through `raised_energy_cff_wip` (`91142139e`), reviewed 2026-09-07–08. Commands ran from the repository root through `nix develop .#default -c ...`, using Cargo 1.97.0 and nextest 0.9.140. Only review documents were changed. A runtime Symbolica license was supplied where needed; no credential is included here.

## Repository checks and tests

| Execution | Observed result | Local evidence |
| --- | --- | --- |
| Shared crate suite | 122 passed, 0 skipped; 30.044 s | `/tmp/raised-review-serial.log` |
| Earlier focused core suite | 68 passed, 628 outside selection; 11.127 s | `/tmp/raised-review-core-tests.log` |
| Earlier six core follow-ups | 6 passed, 690 outside selection; 1.507 s | `/tmp/raised-review-licensed-followup.log` |
| Broader core suite | Initially 206 passed and 4 stack-overflow aborts out of 210; 486 skipped | `/tmp/raised-review-broad-core.log` |
| Four stack-overflow cases rerun | All 4 passed with `RUST_MIN_STACK=67108864`; 13.735 s | `/tmp/raised-review-large-stack.log` |
| Focused API suite | 32 passed, 344 outside selection; 0.037 s | `/tmp/raised-review-api-targeted-tests.log` |
| Targeted integration tests | 7 passed, 252 outside selection; 5.541 s | `/tmp/raised-review-integration.log` |

These counts overlap: the broader core selection repeats earlier core tests, and the four-case rerun repeats its failures. They must not be added into a distinct-test total. The shared suite initially hit the unlicensed concurrency limit in a two-process attempt (12 aborts); its recorded serial rerun passed. Later licensed runs used four workers. The stack setting is an environmental accommodation; passing the rerun does not establish acceptable stack consumption under the default environment. No assertions were weakened to obtain these passes.

The following check/build commands all passed. Clippy emitted no diagnostic from the reviewed code; Cargo noted an upstream `proc-macro-error2` future-compatibility warning.

```sh
nix develop .#default -c cargo fmt --all -- --check
nix develop .#default -c cargo check -p three-dimensional-reps --all-features --locked
nix develop .#default -c cargo clippy -p three-dimensional-reps --all-features --all-targets --locked
nix develop .#default -c cargo check -p gammaloop-api --all-targets --locked
nix develop .#default -c cargo clippy -p gammaloop-api --all-targets --locked
nix develop .#default -c cargo build -p gammaloop-api --bin gammaloop --locked
```

Exact recorded test commands follow. The supplied license is an environment prerequisite for parallel Symbolica use; its value is deliberately omitted.

```sh
nix develop .#default -c cargo nextest run -p three-dimensional-reps --all-features --no-fail-fast --test-threads 1 --retries 0
nix develop .#default -c cargo nextest run -p gammalooprs --lib --locked -E 'test(numerator::energy_degree::tests) or test(uv::approx::local_3d::tests)' --no-fail-fast --test-threads 1 --retries 0
nix develop .#default -c cargo nextest run -p gammalooprs --lib --locked -P test_gammaloop -E 'test(cff::) or test(graph::three_d_source::) or test(graph::lmb::) or test(uv::approx::) or test(integrands::process::) or test(integrands::evaluation::)' --no-fail-fast --test-threads 4 --retries 0
RUST_MIN_STACK=67108864 nix develop .#default -c cargo nextest run -p gammalooprs --lib --locked -P test_gammaloop -E 'test(direct_nested_scalar_banana_replays_complete_cff_at_depth_three) or test(complete_self_energy_taylor_sum_matches_direct_3d_for_raised_lu_jets) or test(factorized_owned_dot_child_cff_matches_direct_3d_for_uncut_self_energy) or test(gl24_direct_3d_modes_preserve_orientation_selector_contracts)' --no-fail-fast --test-threads 4 --retries 0
nix develop .#default -c cargo nextest run -p gammaloop-api --lib --locked -E 'test(command_parser) | test(command_template) | test(commands::bench) | test(commands::profile) | test(completion_offers_run_define) | test(completion_filters_run_define)' --no-fail-fast
RUST_MIN_STACK=67108864 nix develop .#default -c cargo nextest run -p gammaloop-integration-tests --test test_runs --locked -P test_gammaloop -E 'test(test_3d_reps::) or test(repeated_masses::) or test(raised_cut_numerator_cancels_one_propagator_in_both_orientation_modes) or test(bench_cli_profiles_fixed_scalar_triangle_point_and_restores_settings) or test(integrate_writes_numerical_stability_histograms_for_scalar_triangle)' --test-threads 4 --retries 0 --no-fail-fast
```

The earlier six follow-ups used `-P test_gammaloop` and four workers. They covered the GL04 mapping certificate, powered rational identities, signed terminal contours, sampling-scale archive evaluation, deferred hyperdual derivatives, and SymJIT constants. Their summary is retained, but the exact original filter was not recovered; no replacement argv is presented as the historical command.

The full workspace suite, release scalar matrix, NLO numerical campaigns, and Python/standalone end-to-end campaigns were not rerun. The checks below are isolated reproductions, not additional repository tests or full physics campaigns.

## Real CLI reproductions

Use fresh paths when repeating these cases. These commands ran against the built `target/debug/gammaloop`.

**Nested inline variables.** The direct control exited 0; wrapping the same invocation in an inline run exited 1 with `Missing command-block variable 'level'`. Logs: `/tmp/raised-review-api-direct.log`, `/tmp/raised-review-api-nested.log`.

```sh
target/debug/gammaloop -s /tmp/raised-review-api-direct run -D level=warn -c 'set global kv global.display_directive=$(level)'
target/debug/gammaloop -s /tmp/raised-review-api-nested run -c 'run -D level=warn -c "set global kv global.display_directive=$(level)"'
```

**Boot-card variable inheritance.** Save these cards as `/tmp/raised-review-api-boot-direct.toml` and `/tmp/raised-review-api-boot-nested.toml`, respectively:

```toml
commands = ["run inner -D level=warn"]
[[command_blocks]]
name = "inner"
commands = ["set global kv global.display_directive=$(level)"]
```

```toml
commands = ["run outer -D level=warn"]
[[command_blocks]]
name = "inner"
commands = ["set global kv global.display_directive=$(level)"]
[[command_blocks]]
name = "outer"
commands = ["run inner"]
```

The direct card exited 0; the nested card exited 1 during boot validation with missing `level`, despite supplying it at invocation. Logs share the card basenames with `.log`.

```sh
target/debug/gammaloop -s /tmp/raised-review-api-boot-direct /tmp/raised-review-api-boot-direct.toml quit -n
target/debug/gammaloop -s /tmp/raised-review-api-boot-nested /tmp/raised-review-api-boot-nested.toml quit -n
```

**Unsaved state writes.** This exited 0, without `save state` and with `quit -n`, but created `threed_workspace/latest_oriented_expression_path.txt` and the per-graph `oriented_expression.json` under `/tmp/raised-review-api-3drep`. The repository's explicit-save state contract is the expected behavior; the control here is the command's explicit no-save exit. Log: `/tmp/raised-review-api-3drep.log`. Adapt the repository path when reproducing elsewhere.

```sh
target/debug/gammaloop -s /tmp/raised-review-api-3drep run -c 'import model scalars; import graphs /common/dev/gammaloop/higher-power-energies-review/tests/resources/graphs/scalar_box.dot -p box; 3Drep build -p box -i default -g 0 --no-pretty; quit -n'
```

## Compiled public-library reproductions

The following source programs were compiled and executed against built branch libraries, rather than mocked implementations. Save each block to its named `/tmp` file. These are the recorded effective compiler arguments (the original parser launcher selected its rlib with a shell glob); the bounds case was additionally recompiled during appendix preparation. Build artifacts have local hashes, so replace them with a matching built library when reproducing elsewhere.

```sh
nix develop .#default -c rustc --edition=2024 /tmp/raised-review-eval-parser.rs --extern three_dimensional_reps=target/debug/deps/libthree_dimensional_reps-397853b131b47fb5.rlib -L dependency=target/debug/deps -o /tmp/raised-review-eval-parser
/tmp/raised-review-eval-parser
nix develop .#default -c rustc --edition=2024 /tmp/raised-review-validator.rs -L dependency=target/debug/deps --extern three_dimensional_reps=target/debug/deps/libthree_dimensional_reps-da5aed3622a15d11.rlib -o /tmp/raised-review-validator
/tmp/raised-review-validator
nix develop .#default -c rustc --edition=2024 /tmp/raised-review-bounds.rs --extern three_dimensional_reps=target/debug/deps/libthree_dimensional_reps-da5aed3622a15d11.rlib -L dependency=target/debug/deps -o /tmp/raised-review-bounds-confirmed
nix develop .#default -c /tmp/raised-review-bounds-confirmed
```

**Diagnostic expression parser**, `/tmp/raised-review-eval-parser.rs`:

```rust
use std::collections::BTreeMap;
use three_dimensional_reps::{generate_3d_expression, Generate3DExpressionOptions, ParsedGraph};
use three_dimensional_reps::eval::{evaluate_expression, EvaluationInput};
fn main() {
    let graph = ParsedGraph { internal_edges: vec![], external_edges: vec![],
        initial_state_cut_edges: vec![], loop_names: vec![], external_names: vec![],
        node_name_to_internal: BTreeMap::new() };
    let generated = generate_3d_expression(&graph, &Generate3DExpressionOptions::default()).unwrap();
    let input = EvaluationInput { external_momenta: vec![], loop_spatial_momenta: vec![],
        masses: vec![], uniform_scale: None };
    for source in ["-2**2", "-(2**2)", "2**3**2", "2**4294967296"] {
        let result = evaluate_expression(&graph, &generated.expression, source, &input).unwrap();
        println!("{source} => {}", result.value);
    }
}
```

Observed: `-2**2 => 4`, `-(2**2) => -4`, `2**3**2 => 64`, `2**4294967296 => 1` (`/tmp/raised-review-eval-parser.log`). Conventional power precedence/association give −4 and 512 for the first/third forms; the last silently narrows the exponent. This is the public diagnostic evaluator, not compiled production integrand evaluation.

**Malformed graph validation**, `/tmp/raised-review-validator.rs`:

```rust
use std::collections::BTreeMap;
use three_dimensional_reps::{validate_parsed_graph, MomentumSignature, ParsedGraph};
use three_dimensional_reps::graph_io::ParsedGraphInternalEdge;
fn main() {
    let graph = ParsedGraph {
        internal_edges: vec![ParsedGraphInternalEdge {
            edge_id: 0, tail: 0, head: 1, label: "q".into(), mass_key: None,
            signature: MomentumSignature { loop_signature: vec![1], external_signature: vec![] },
            had_pow: false }],
        external_edges: vec![], initial_state_cut_edges: vec![],
        loop_names: vec![], external_names: vec![],
        node_name_to_internal: BTreeMap::from([("a".into(), 0), ("b".into(), 1)]) };
    let result = std::panic::catch_unwind(|| validate_parsed_graph(&graph));
    println!("validation_panicked={}", result.is_err());
}
```

Observed `validation_panicked=true`, with an out-of-bounds panic at `validator.rs:52`. The input is deliberately malformed; the diagnostic validator should return a structured rejection.

**Preserved-tree energy-bound validation**, `/tmp/raised-review-bounds.rs`:

```rust
use std::collections::BTreeMap;
use three_dimensional_reps::{generate_3d_expression, Generate3DExpressionOptions, MomentumSignature, ParsedGraph};
use three_dimensional_reps::graph_io::ParsedGraphInternalEdge;
fn main() {
    let graph = ParsedGraph {
        internal_edges: vec![ParsedGraphInternalEdge {
            edge_id: 0, tail: 0, head: 1, label: "tree".into(), mass_key: None,
            signature: MomentumSignature { loop_signature: vec![], external_signature: vec![1] },
            had_pow: false }],
        external_edges: vec![], initial_state_cut_edges: vec![],
        loop_names: vec![], external_names: vec!["p".into()],
        node_name_to_internal: BTreeMap::from([("a".into(), 0), ("b".into(), 1)]) };
    let options = Generate3DExpressionOptions {
        energy_degree_bounds: Some(vec![(999, 2)]),
        preserve_internal_edges_as_four_d_denominators: vec![0], ..Default::default() };
    match generate_3d_expression(&graph, &options) {
        Ok(generated) => println!("accepted invalid source bounds: {:?}", generated.source_energy_degree_bounds),
        Err(error) => println!("rejected invalid bounds: {error}"),
    }
}
```

Observed `accepted invalid source bounds: [(999, 2)]` (`/tmp/raised-review-bounds.log`). Edge 999 does not exist. This demonstrates an early-return validation inconsistency, not a wrong value for valid physics input.

## Executed numerical-assertion counterexamples

These isolate the test predicates. They show what the assertions accept; they do not claim that the current physics evaluators produced these invalid totals. Both programs were compiled using the actual branch `F<ArbPrec>` implementation:

```sh
nix develop .#default -c rustc --edition=2024 /tmp/raised-review-nan-repro.rs --extern gammalooprs=target/debug/deps/libgammalooprs-439f918cbac6a66b.rlib -L dependency=target/debug/deps -o /tmp/raised-review-nan-repro
nix develop .#default -c /tmp/raised-review-nan-repro
nix develop .#default -c rustc --edition=2024 /tmp/raised-review-predicates.rs --extern gammalooprs=target/debug/deps/libgammalooprs-439f918cbac6a66b.rlib -L dependency=target/debug/deps -o /tmp/raised-review-predicates
nix develop .#default -c /tmp/raised-review-predicates
```

Replace that rlib hash when building elsewhere. `/tmp/raised-review-nan-repro.rs`:

```rust
use gammalooprs::utils::{ArbPrec, F};
fn main() {
    let invalid = F::<ArbPrec>::from_f64(f64::NAN);
    let tolerance = F::<ArbPrec>::from_f64(1.0e-30);
    println!("nan > tolerance: {}", invalid > tolerance);
    println!("nan <= tolerance: {}", invalid <= tolerance);
}
```

Both comparisons printed `false`. Thus the failure aggregators at `cff/mod.rs:4418` and similar locations do not record a NaN error. `/tmp/raised-review-predicates.rs` exercises the f64 total predicate and the small-result Arb floor:

```rust
use gammalooprs::utils::{ArbPrec, F};
fn main() {
    for actual in [f64::NAN, f64::INFINITY] {
        let expected = 1.0_f64;
        let scale = (actual * actual).sqrt().max(expected).max(1.0);
        let tolerance = 1.0e-10 * scale;
        let distance = ((actual - expected) * (actual - expected)).sqrt();
        println!("total={actual}: failure recorded={}", distance > tolerance);
    }
    let actual = F::<ArbPrec>::from_f64(1.0e-18);
    let reference = F::<ArbPrec>::from_f64(2.0e-18);
    let distance = reference.clone() - actual;
    let relative = distance.clone() / reference;
    let floor_accepts = distance / F::<ArbPrec>::from_f64(1.0) <= F::<ArbPrec>::from_f64(1.0e-14);
    println!("small-result relative error={relative:e}; floor accepts={floor_accepts}");
}
```

Observed: `total=NaN: failure recorded=false`; `total=inf: failure recorded=false`; small-result relative error `0.5`, with `floor accepts=true`. These reproduce the acceptance conditions in `tests/tests/test_runs/utils.rs:635` and `scalar_3L_cross_section_inspects.rs:969`. Require finite values and affirmative successful comparisons; justify any absolute tolerance for certified zeros separately.
