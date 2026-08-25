#import "../../shared.typ": callout

#let quickstart-rust = [
= Using GammaLoop from Rust

Embed GammaLoop when a Rust program should own the state, command session, and typed evaluation
result directly. This example generates the same scalar bubble as the CLI and Python routes.

== Create a project

The ordinary crates.io command becomes available with the automated `0.3.4` release:

// docs-example: syntax
```sh
cargo new first-gammaloop
cd first-gammaloop
cargo add gammaloop-api@0.3.4 ndarray@0.17 eyre@0.6
```

Before that release, replace the `gammaloop-api` line with the tested release revision:

// docs-example: syntax
```sh
cargo add gammaloop-api \
  --git https://github.com/alphal00p/gammaloop.git \
  --rev eb1dab23c24345def88226ef10874f0dde6c4aa5
cargo add ndarray@0.17 eyre@0.6
```

== Drive one state and evaluate it

Replace `src/main.rs` with:

// docs-example: compile gammaloop-rust-quickstart
```rust
use eyre::Result;
use gammaloop_api::{
    commands::evaluate_samples::{evaluate_sample, EvaluateSamples},
    state::CommandHistory,
    StateLoadOption,
};
use ndarray::arr2;
use std::time::{SystemTime, UNIX_EPOCH};

fn main() -> Result<()> {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let state_folder = std::env::temp_dir().join(format!("gammaloop-quickstart-{nonce}"));
    let mut loaded = StateLoadOption {
        state_folder: Some(state_folder.clone()),
        ..StateLoadOption::default()
    }
    .load()?;

    for raw in [
        "import model scalars-default.json",
        "set global kv global.generation.uv.generate_integrated=false",
        concat!(
            "generate amp scalar_1 > scalar_1 [{1}] ",
            "--allowed-vertex-interactions V_3_SCALAR_122 ",
            "-p bubble -i one_loop"
        ),
    ] {
        let command = CommandHistory::from_raw_string(raw)?;
        let _ = loaded.cli_session().execute_command(command)?;
    }

    let point = arr2(&[[0.1, 0.2, 0.3]]);
    let result = evaluate_sample(
        &mut loaded.state,
        &EvaluateSamples {
            process_id: Some(0),
            integrand_name: Some("one_loop".into()),
            use_arb_prec: false,
            minimal_output: false,
            return_generated_events: None,
            momentum_space: false,
            points: point.view(),
            integrator_weights: None,
            discrete_dims: None,
            graph_names: None,
            orientations: None,
        },
    )?;

    println!("{result}");
    std::fs::remove_dir_all(state_folder)?;
    Ok(())
}
```

Compile and run with GammaLoop's public Symbolica build selector:

// docs-example: syntax
```sh
SYMBOLICA_OEM_LICENSE=SYMBOLICA_OEM_GAMMALOOP cargo run
```

#callout("The facade owns the supported embedding boundary", [
  Prefer `StateLoadOption`, `CliSession`, and the typed request objects from `gammaloop-api`.
  `gammalooprs` exposes lower-level implementation types, but their presence in Rustdoc is not a
  stability promise. The #link("reference/interfaces/")[interface guide] maps the supported
  workflow to the complete Rust reference.
])

The docs harness compiles this source against the current workspace. Running it additionally
requires the native toolchain and a Symbolica license mode appropriate for the user.
]
