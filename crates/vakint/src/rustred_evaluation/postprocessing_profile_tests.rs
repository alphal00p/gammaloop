//! Opt-in timing and digests of one complete public evaluation.
//!
//! Runs a `docs/benchmark` input record (`{"loop_count", "input"}`, as written
//! by `examples/benchmark_suite.rs prepare-*`) through the public `evaluate`
//! with the benchmark harness's settings (MSbar, `L + 1` ε terms, dot-product
//! notation, 32 requested digits, shipped step caches) and reports the wall
//! time, the SHA-256 of the result's canonical string and the numerical
//! Laurent coefficients at unit mass and scale, as the harness does. With
//! `RUST_LOG=vakint::postprocessing=debug` every post-reduction phase logs its
//! wall time, which separates the native reduction from the assembly after it.
//!
//! Env: `VAKINT_TEST_EXTERNAL_DATA` (data root), `VAKINT_PROFILE_INPUT` (input
//! record), `VAKINT_PROFILE_OUT` (JSON report), optional
//! `VAKINT_PROFILE_MODE` (`numerical` (default), `symbolic` or `retain`),
//! `VAKINT_N_CORES` (default 8) and `VAKINT_PROFILE_BACKEND` (`rustred`
//! (default) or `matad`: the harness's MATAD full call, with the FORM
//! executable `VAKINT_PROFILE_FORM`, default `form`, and its work files beside
//! the report).

use std::time::Instant;

use serde_json::json;
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore};

use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, MATADOptions,
    RustRedEvaluationOptions, RustRedMasterMode, Vakint, VakintSettings,
    utils::vakint_macros::{vk_parse, vk_symbol},
};

fn required(name: &str) -> String {
    std::env::var(name).unwrap_or_else(|_| panic!("{name} is required"))
}

#[test_log::test]
#[ignore = "diagnostic: one complete public evaluation of a benchmark input; needs VAKINT_TEST_EXTERNAL_DATA, VAKINT_PROFILE_INPUT and VAKINT_PROFILE_OUT"]
fn public_evaluation_postprocessing_profile() {
    let root = required("VAKINT_TEST_EXTERNAL_DATA");
    let record: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(required("VAKINT_PROFILE_INPUT")).unwrap())
            .unwrap();
    let loops = record["loop_count"].as_u64().unwrap() as usize;
    let mode = match std::env::var("VAKINT_PROFILE_MODE").as_deref() {
        Ok("retain") => RustRedMasterMode::Retain,
        Ok("symbolic") => RustRedMasterMode::Symbolic,
        _ => RustRedMasterMode::Numerical,
    };
    let n_cores = std::env::var("VAKINT_N_CORES")
        .ok()
        .map_or(8, |value| value.parse().unwrap());
    let matad = std::env::var("VAKINT_PROFILE_BACKEND").as_deref() == Ok("matad");
    let mut settings = VakintSettings {
        evaluation_order: if matad {
            EvaluationOrder::matad_only(Some(MATADOptions::default()))
        } else {
            EvaluationOrder(vec![EvaluationMethod::RustRed(RustRedEvaluationOptions {
                master_mode: mode,
            })])
        },
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        number_of_terms_in_epsilon_expansion: loops as i64 + 1,
        use_dot_product_notation: true,
        form_exe_path: if matad {
            std::env::var("VAKINT_PROFILE_FORM").unwrap_or_else(|_| "form".into())
        } else {
            "/profile-rustred-must-not-use-form".into()
        },
        temporary_directory: matad.then(|| format!("{}.work", required("VAKINT_PROFILE_OUT"))),
        run_time_decimal_precision: 32,
        ..Default::default()
    };
    settings.n_cores = n_cores;
    let started = Instant::now();
    let mut engine = Vakint::new().unwrap();
    engine.load_external_data(&root, loops).unwrap();
    let loading = started.elapsed().as_secs_f64();
    let input = vk_parse!(record["input"].as_str().unwrap()).unwrap();
    let started = Instant::now();
    let outcome = engine.evaluate(&settings, input.as_view());
    let seconds = started.elapsed().as_secs_f64();
    let mut report = json!({
        "schema": "vakint.postprocessing-profile.v1",
        "input_sha256": record["input_sha256"],
        "loop_count": loops,
        "mode": format!("{mode:?}"),
        "backend": if matad { "matad" } else { "rustred" },
        "n_cores": n_cores,
        "loading_seconds": loading,
        "call_seconds": seconds,
    });
    match &outcome {
        Ok(value) => {
            let canonical = value.to_canonical_string();
            report["output_sha256"] = json!(format!("{:x}", Sha256::digest(canonical.as_bytes())));
            report["output_bytes"] = json!(canonical.len());
            if canonical.len() <= 200_000 {
                report["output"] = json!(canonical);
            }
            if mode == RustRedMasterMode::Numerical {
                let mut evaluated = value.clone();
                for name in ["mUV", "mUVsq", "muvsq"] {
                    evaluated = evaluated
                        .replace(Atom::var(vk_symbol!(name)))
                        .with(Atom::one());
                }
                evaluated = evaluated
                    .replace(Atom::var(symbolica::symbol!("UFO::mUV")))
                    .with(Atom::one())
                    .replace(Atom::var(vk_symbol!(&settings.mu_r_sq_symbol)))
                    .with(Atom::one());
                report["numerical"] = match engine.numerical_evaluation(
                    &settings,
                    evaluated.as_view(),
                    &Default::default(),
                    &Default::default(),
                    None,
                ) {
                    Ok((result, _)) => json!({
                        "coefficients": result.coefficients().iter().map(|(power, c)| {
                            let p = c.precision();
                            json!({
                                "epsilon_power": power,
                                "real": c.value().re.to_string(),
                                "imaginary": c.value().im.to_string(),
                                "real_absolute_error": p.real().absolute_error().to_string(),
                                "imaginary_absolute_error": p.imaginary().absolute_error().to_string(),
                                "real_significant_digits": p.real().significant_digits(),
                                "imaginary_significant_digits": p.imaginary().significant_digits(),
                                "real_exact_zero": p.real().exact_zero(),
                                "imaginary_exact_zero": p.imaginary().exact_zero(),
                            })
                        }).collect::<Vec<_>>(),
                        "warnings": result.warnings().iter().map(ToString::to_string).collect::<Vec<_>>(),
                    }),
                    Err(error) => json!({"error": error.to_string()}),
                };
            }
        }
        Err(error) => report["error"] = json!(error.to_string()),
    }
    std::fs::write(
        required("VAKINT_PROFILE_OUT"),
        serde_json::to_vec_pretty(&report).unwrap(),
    )
    .unwrap();
    println!(
        "{}",
        json!({"call_seconds": seconds, "ok": outcome.is_ok()})
    );
}
