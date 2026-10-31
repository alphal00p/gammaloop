#![cfg(feature = "feynkit-ingress")]
//! Physical graph/tensor ingress, lower-loop engine parity and numerical references.

use ::vakint;
#[path = "fixtures/qcd_ladder/mod.rs"]
#[allow(dead_code)]
mod fixture;

use symbolica::atom::AtomCore;
use vakint::{OrderedWeightArithmetic, RustRedMonitorOptions};

fn golden_precision_target(golden_minimum: u64, selected_request: u64) -> u64 {
    golden_minimum.min(selected_request)
}

#[test]
fn golden_precision_target_respects_the_selected_source_cap() {
    assert_eq!(golden_precision_target(300, 300), 300);
    assert_eq!(golden_precision_target(300, 200), 200);
    assert_eq!(golden_precision_target(200, 300), 200);
}

/// Reduce one physical rung with exact and with modular Ordered weights. Both
/// must return the same retained expression, which must equal the historical
/// depth-first result frozen in `tests/resources/ordered_rungs` before
/// depth-first reduction was removed from Vakint, and match the recorded
/// high-precision numerical reference.
fn compare(loop_count: usize, projection: fixture::Projection) -> fixture::Result<()> {
    let temporary = tempfile::tempdir()?;
    let directory = if let Some(root) = std::env::var_os("VAKINT_ORDERED_RUN_EVIDENCE") {
        let root = std::path::PathBuf::from(root);
        let directory = root.join(format!("L{loop_count}-{projection:?}"));
        if ["exact", "modular"]
            .iter()
            .any(|arm| directory.join(arm).join("report.json").exists())
        {
            return Err(
                "preserve completed physical evidence and select a fresh evidence root".into(),
            );
        }
        std::fs::create_dir_all(root)?;
        // Fresh case directories keep graph preparation inside the measured CI
        // invocation and preserve interrupted evidence instead of reusing it.
        std::fs::create_dir(&directory)?;
        directory
    } else {
        temporary.path().to_owned()
    };
    let input_directory = directory.join("input");
    let prepared = fixture::comparison::prepared(&input_directory, projection, loop_count)?;
    assert!(!prepared.input.is_zero());
    assert_eq!(prepared.loop_count, loop_count);
    // Read again through the authenticated bundle boundary: repeated validation
    // must not regenerate graphs, traces, Taylor coefficients or routing.
    let loaded = fixture::LadderFixture::load_prepared_at_loop_count(
        &input_directory,
        projection,
        loop_count,
    )?;
    assert_eq!(prepared.input, loaded.input);
    let run = |name, arithmetic, numerical_reference| {
        fixture::comparison::run(
            &loaded,
            &directory.join(name),
            &fixture::comparison::RunOptions {
                arithmetic,
                artifacts: std::env::var_os("VAKINT_TEST_EXTERNAL_DATA").map(Into::into),
                monitor: RustRedMonitorOptions::default(),
                numerical_reference,
            },
        )
    };
    // The exact arm carries the numerical reference; the modular arm must
    // reproduce its retained expression exactly.
    let exact = run("exact", OrderedWeightArithmetic::exact(), true)?;
    let modular = run(
        "modular",
        fixture::comparison::modular_arithmetic_from_environment()?,
        false,
    )?;
    assert_eq!(
        exact.report["retained"]["sha256"], modular.report["retained"]["sha256"],
        "exact and modular Ordered weights retained different expressions"
    );
    // Same input, data and bindings; the modular arm has no numerical reference.
    let mut retained_only = fixture::comparison::CompletedRun {
        value: exact.value.clone(),
        report: exact.report.clone(),
    };
    retained_only.report["numerical_reference"] = serde_json::Value::Null;
    fixture::comparison::compare(&retained_only, &modular)?;
    for result in [&exact, &modular] {
        assert_eq!(result.report["strategy"], "Ordered");
        assert_eq!(result.report["intrinsic_memoization"], false);
        assert_eq!(
            result.report["cache_policy"]["prepared_cache_loaded"],
            false
        );
        assert_eq!(
            result.report["cache_policy"]["previous_run_cache_restored"],
            false
        );
        #[cfg(feature = "experimental-rustred")]
        {
            let timing = &result.report["reduction_math_timing"];
            assert_eq!(timing["scope_id"], "vakint.rustred.reduction-math.v1");
            assert_eq!(timing["complete"], true);
            assert!(timing["interval_count"].as_u64().unwrap_or(0) > 0);
            let math_seconds = result.report["reduction_math_seconds"]
                .as_f64()
                .ok_or("completed native reduction has no mathematical timing")?;
            assert!(math_seconds.is_finite() && math_seconds >= 0.0);
            assert!(
                math_seconds
                    <= result.report["retained_reduction_seconds"]
                        .as_f64()
                        .unwrap()
            );
        }
    }
    // The separately evaluated numerical reference is present, while the
    // completed capture above covers only the earlier retained call.
    assert!(exact.report["numerical_reference"]["seconds"].is_number());
    let case = format!(
        "{loop_count}L-{}",
        match projection {
            fixture::Projection::ZeroMomentumMetric => "zero-momentum",
            fixture::Projection::QuadraticTransverse => "quadratic",
        }
    );
    let resources =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/resources/ordered_rungs");
    let depth_first =
        std::fs::read_to_string(resources.join(format!("{case}.depth-first-retained.symbolica")))?;
    vakint::Vakint::initialize_vakint_symbols();
    let depth_first = symbolica::atom::Atom::parse(&depth_first, "vakint", Default::default())?;
    assert!(
        (fixture::comparison::canonical_retained(&exact.value)?
            - fixture::comparison::canonical_retained(&depth_first)?)
        .together()
        .expand()
        .is_zero(),
        "{case}: Ordered differs from the historical depth-first result"
    );
    let reference: serde_json::Value =
        serde_json::from_slice(&std::fs::read(resources.join(format!("{case}.json")))?)?;
    assert_eq!(reference["normalization"], "MSbar");
    assert_eq!(
        reference["source_input_sha256"],
        exact.report["input_sha256"]
    );
    assert_eq!(
        reference["evaluated_input_sha256"],
        exact.report["evaluated_input_sha256"]
    );
    assert_eq!(
        reference["dimension_binding"],
        exact.report["dimension_binding"]
    );
    if let Some(minimum_digits) = reference["minimum_achieved_decimal_precision"].as_u64() {
        let requested =
            exact.report["numerical_reference"]["selection"]["requested_decimal_precision"]
                .as_u64()
                .ok_or("completed reference has no selected precision request")?;
        let minimum_digits = golden_precision_target(minimum_digits, requested);
        for coefficient in exact.report["numerical_reference"]["coefficients"]
            .as_array()
            .ok_or("completed reference has no coefficient inventory")?
        {
            for component in ["real", "imaginary"] {
                if coefficient[format!("{component}_exact_zero")] != true {
                    assert!(
                        coefficient[format!("{component}_significant_digits")]
                            .as_u64()
                            .is_some_and(|digits| digits >= minimum_digits),
                        "epsilon^{} {component} falls below its validated {minimum_digits}-digit regression target",
                        coefficient["epsilon_power"],
                    );
                }
            }
        }
    }
    fixture::comparison::compare_numerical_references(
        &exact.report["numerical_reference"],
        &reference["numerical_reference"],
    )?;
    assert!(!exact.value.is_zero());
    Ok(())
}

#[test]
fn one_loop_zero_momentum_rung_matches_depth_first_and_reference() -> fixture::Result<()> {
    compare(1, fixture::Projection::ZeroMomentumMetric)
}

#[test]
fn one_loop_quadratic_rung_matches_depth_first_and_reference() -> fixture::Result<()> {
    compare(1, fixture::Projection::QuadraticTransverse)
}

#[test]
fn one_loop_rung_projections_obey_independent_dirac_and_tadpole_identity() -> fixture::Result<()> {
    use symbolica::atom::Atom;
    let directory = tempfile::tempdir()?;
    let run = |name, projection| -> fixture::Result<_> {
        let prepared = fixture::comparison::prepared(
            &directory.path().join(format!("input-{name}")),
            projection,
            1,
        )?;
        fixture::comparison::run(
            &prepared,
            &directory.path().join(name),
            &fixture::comparison::RunOptions {
                arithmetic: OrderedWeightArithmetic::exact(),
                artifacts: None,
                monitor: RustRedMonitorOptions::default(),
                numerical_reference: false,
            },
        )
    };
    let zero = run("zero", fixture::Projection::ZeroMomentumMetric)?;
    let quadratic = run("quadratic", fixture::Projection::QuadraticTransverse)?;
    // Independent Clifford trace and tadpole recurrence give, before their
    // common native coupling/phase factor, 8*d*(2-d)*I1 and
    // d*(d-2)/(12*M²)*I1. This checks the relative projector normalization,
    // external-momentum Taylor terms, color average and UV mass scaling.
    let mass_squared = Atom::var(vakint::vakint_symbol!("mUV")).pow(2);
    assert!(
        (fixture::comparison::canonical_retained(&zero.value)?
            + Atom::num(96)
                * mass_squared
                * fixture::comparison::canonical_retained(&quadratic.value)?)
        .together()
        .expand()
        .is_zero()
    );
    Ok(())
}

#[test]
fn two_loop_zero_momentum_rung_matches_depth_first_and_reference() -> fixture::Result<()> {
    compare(2, fixture::Projection::ZeroMomentumMetric)
}

#[test]
fn two_loop_quadratic_rung_matches_depth_first_and_reference() -> fixture::Result<()> {
    compare(2, fixture::Projection::QuadraticTransverse)
}

#[test]
fn three_loop_zero_momentum_rung_matches_depth_first_and_reference() -> fixture::Result<()> {
    compare(3, fixture::Projection::ZeroMomentumMetric)
}

#[test]
fn three_loop_quadratic_rung_matches_depth_first_and_reference() -> fixture::Result<()> {
    compare(3, fixture::Projection::QuadraticTransverse)
}

#[test]
fn four_loop_zero_momentum_rung_matches_depth_first_and_reference() -> fixture::Result<()> {
    compare(4, fixture::Projection::ZeroMomentumMetric)
}

#[test]
fn four_loop_quadratic_rung_matches_depth_first_and_reference() -> fixture::Result<()> {
    compare(4, fixture::Projection::QuadraticTransverse)
}

#[test]
#[ignore = "ordered-only full 5L reduction against the accepted saved zero-momentum result; requires explicit prepared input and durable evidence directories"]
fn five_loop_zero_momentum_ordered_matches_saved_depth_first_result() -> fixture::Result<()> {
    let completed = ordered_five_loop(fixture::Projection::ZeroMomentumMetric)?;
    let expected_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples/data/five_loop_q0/q0-retained.expr");
    let payload = std::fs::read(expected_path)?;
    let evidence = std::path::PathBuf::from(
        std::env::var_os("VAKINT_ORDERED_RUN_EVIDENCE")
            .ok_or("set VAKINT_ORDERED_RUN_EVIDENCE to the completed run directory")?,
    );
    fixture::comparison::historical::compare_saved_zero_momentum(&payload, &completed, &evidence)?;
    Ok(())
}

#[test]
#[ignore = "ordered-only quadratic 5L reduction; run only after ordered zero-momentum acceptance; requires explicit prepared input and durable evidence directories"]
fn five_loop_quadratic_ordered_retains_declared_masters() -> fixture::Result<()> {
    ordered_five_loop(fixture::Projection::QuadraticTransverse)?;
    Ok(())
}

fn ordered_five_loop(
    projection: fixture::Projection,
) -> fixture::Result<fixture::comparison::CompletedRun> {
    let input = std::path::PathBuf::from(
        std::env::var_os("VAKINT_QCD_LADDER_PREPARED_INPUT")
            .ok_or("set VAKINT_QCD_LADDER_PREPARED_INPUT to the authenticated 5L scalar bundle")?,
    );
    let evidence = std::path::PathBuf::from(
        std::env::var_os("VAKINT_ORDERED_RUN_EVIDENCE")
            .ok_or("set VAKINT_ORDERED_RUN_EVIDENCE to a fresh durable workspace directory")?,
    );
    if evidence.join("report.json").exists() {
        return Err(
            "preserve the completed physical run and select a new evidence directory".into(),
        );
    }
    let prepared = fixture::LadderFixture::load_prepared_at_loop_count(&input, projection, 5)?;
    if projection == fixture::Projection::ZeroMomentumMetric {
        use sha2::Digest;
        assert_eq!(
            format!(
                "{:x}",
                sha2::Sha256::digest(prepared.input.to_canonical_string().as_bytes())
            ),
            "9814c3b49a1b486b9790bbe5458d8bc582c3fb3c299ba6a1a0362db8ba71cb16",
            "the saved historical result belongs to a different physical input"
        );
    }
    fixture::comparison::run(
        &prepared,
        &evidence,
        &fixture::comparison::RunOptions {
            arithmetic: fixture::comparison::arithmetic_from_environment()?,
            artifacts: std::env::var_os("VAKINT_TEST_EXTERNAL_DATA").map(Into::into),
            monitor: RustRedMonitorOptions {
                progress_snapshot: Some(evidence.join("monitor.json")),
                ..Default::default()
            },
            numerical_reference: false,
        },
    )
}
