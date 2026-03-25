#!/usr/bin/env rust-script
//! ```cargo
//! [dependencies]
//! color-eyre = "0.6"
//! ndarray = "0.17"
//! gammaloop-api = { path = "../../../../crates/gammaloop-api", default-features = false, features = ["cli"] }
//! gammalooprs = { path = "../../../../crates/gammalooprs" }
//!
//! [patch.crates-io]
//! #graphica = { git = "https://github.com/benruijl/symbolica", branch = "dev" }
//! #numerica = { git = "https://github.com/benruijl/symbolica", branch = "dev" }
//! #symbolica = { git = "https://github.com/benruijl/symbolica", branch = "dev" }
//! graphica = { git = "https://github.com/benruijl/symbolica", rev = "650ba97bf3da7cf2ff5ada92875f92d5f71e7a31" }
//! numerica = { git = "https://github.com/benruijl/symbolica", rev = "650ba97bf3da7cf2ff5ada92875f92d5f71e7a31" }
//! symbolica = { git = "https://github.com/benruijl/symbolica", rev = "650ba97bf3da7cf2ff5ada92875f92d5f71e7a31" }
//!
//! ```
//!
//! Run from the repository root with `rust-script`, for example:
//! `NO_SYMBOLICA_OEM_LICENSE=1 EXTRA_MACOS_LIBS_FOR_GNU_GCC=T rust-script --debug examples/api/rust/epem_a_ddxg_xs_LO/inspect_events.rs`
//!
//! If you change the OEM-license environment at build time, add `--force` so
//! `rust-script` does not reuse an incompatible cached build.

use std::{env, path::PathBuf};

use color_eyre::Result;
use gammaloop_api::{
    commands::evaluate_samples::{
        evaluate_sample, evaluate_sample_precise, EvaluateSamples, EvaluateSamplesPrecise,
    },
    state::CommandHistory,
    StateLoadOption,
};
use ndarray::{arr2, Array2};

fn example_dir() -> Result<PathBuf> {
    Ok(PathBuf::from(file!())
        .canonicalize()?
        .parent()
        .unwrap()
        .to_path_buf())
}

fn main() -> Result<()> {
    let example_dir = example_dir()?;
    let run_card = example_dir.join("run.toml");
    let state_dir = example_dir.join("state");

    let mut loaded = StateLoadOption {
        fresh_state: true,
        boot_commands_path: Some(run_card),
        state_folder: Some(state_dir),
        ..StateLoadOption::default()
    }
    .load()?;

    // The example card ships a small block of `display quantities`,
    // `display observables`, and `display selectors` commands. Run it through
    // the CLI session before switching to the lower-level evaluation APIs.
    {
        let mut session = loaded.cli_session();
        let _ = session.execute_command(CommandHistory::from_raw_string(
            "run display_named_settings_examples",
        )?)?;
    }

    let single_point: Array2<f64> = arr2(&[[0.17, 0.31, 0.53, 0.23, 0.41, 0.67]]);
    let single_result = evaluate_sample(
        &mut loaded.state,
        &EvaluateSamples {
            process_id: None,
            integrand_name: None,
            use_arb_prec: false,
            minimal_output: false,
            momentum_space: false,
            points: single_point.view(),
            integrator_weights: None,
            discrete_dims: None,
            graph_names: None,
            orientations: None,
        },
    )?;
    println!("\n== evaluate_sample ==\n");
    println!("{single_result}");

    // let momentum_point: Array2<f64> = arr2(&[[0.11, -0.07, 0.19, -0.13, 0.05, 0.29]]);
    // let momentum_result = evaluate_sample(
    //     &mut loaded.state,
    //     &EvaluateSamples {
    //         process_id: None,
    //         integrand_name: None,
    //         use_arb_prec: false,
    //         minimal_output: false,
    //         momentum_space: true,
    //         points: momentum_point.view(),
    //         integrator_weights: None,
    //         discrete_dims: None,
    //         graph_names: None,
    //         orientations: None,
    //     },
    // )?;
    // println!("\n== momentum-space evaluate_sample ==\n");
    // println!("{momentum_result}");

    // let precise_result = evaluate_sample_precise(
    //     &mut loaded.state,
    //     &EvaluateSamplesPrecise {
    //         process_id: None,
    //         integrand_name: None,
    //         use_arb_prec: true,
    //         minimal_output: false,
    //         momentum_space: false,
    //         points: single_point.view(),
    //         integrator_weights: None,
    //         discrete_dims: None,
    //         graph_names: None,
    //         orientations: None,
    //     },
    // )?;
    // println!("\n== evaluate_sample_precise ==\n");
    // println!("{precise_result}");

    // let batch_points: Array2<f64> = arr2(&[
    //     [0.17, 0.31, 0.53, 0.23, 0.41, 0.67],
    //     [0.11, 0.29, 0.47, 0.19, 0.37, 0.59],
    // ]);
    // let batch_results = EvaluateSamples {
    //     process_id: None,
    //     integrand_name: None,
    //     use_arb_prec: false,
    //     minimal_output: false,
    //     momentum_space: false,
    //     points: batch_points.view(),
    //     integrator_weights: None,
    //     discrete_dims: None,
    //     graph_names: None,
    //     orientations: None,
    // }
    // .run(&mut loaded.state)?;

    // println!("\n== evaluate_samples ==\n");
    // println!("{batch_results}");

    // let precise_batch_results = EvaluateSamplesPrecise {
    //     process_id: None,
    //     integrand_name: None,
    //     use_arb_prec: true,
    //     minimal_output: true,
    //     momentum_space: false,
    //     points: batch_points.view(),
    //     integrator_weights: None,
    //     discrete_dims: None,
    //     graph_names: None,
    //     orientations: None,
    // }
    // .run(&mut loaded.state)?;
    // println!("\n== evaluate_samples_precise ==\n");
    // println!("{precise_batch_results}");

    Ok(())
}
