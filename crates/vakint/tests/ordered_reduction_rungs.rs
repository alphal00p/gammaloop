#![cfg(feature = "feynkit-ingress")]
//! Physical graph, tensor, ingress and exact/native numerical engine parity.

use ::vakint;
#[path = "fixtures/qcd_ladder/mod.rs"]
#[allow(dead_code)]
mod fixture;

use vakint::{RustRedMonitorOptions, RustRedReductionStrategy};

fn compare(loop_count: usize, projection: fixture::Projection) -> fixture::Result<()> {
    let directory = tempfile::tempdir()?;
    let input_directory = directory.path().join("input");
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
    let run = |name, strategy| {
        fixture::comparison::run(
            &loaded,
            &directory.path().join(name),
            &fixture::comparison::RunOptions {
                strategy,
                artifacts: std::env::var_os("VAKINT_TEST_EXTERNAL_DATA").map(Into::into),
                monitor: RustRedMonitorOptions::default(),
                numerical_reference: loop_count <= 4,
            },
        )
    };
    let dfs = run("dfs", RustRedReductionStrategy::DepthFirst)?;
    let ordered = run("ordered", RustRedReductionStrategy::Ordered)?;
    fixture::comparison::compare(&dfs, &ordered)?;
    #[cfg(feature = "experimental-rustred")]
    for result in [&dfs, &ordered] {
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
        // The separately evaluated numerical reference is present, while the
        // completed capture above covers only the earlier retained call.
        assert!(result.report["numerical_reference"]["seconds"].is_number());
    }
    let reference = match (loop_count, projection) {
        (1, fixture::Projection::ZeroMomentumMetric) => Some("1L-zero-momentum.json"),
        (1, fixture::Projection::QuadraticTransverse) => Some("1L-quadratic.json"),
        _ => None,
    };
    if let Some(reference) = reference {
        let reference: serde_json::Value = serde_json::from_slice(&std::fs::read(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/resources/ordered_rungs")
                .join(reference),
        )?)?;
        assert_eq!(reference["normalization"], "MSbar");
        assert_eq!(reference["source_input_sha256"], dfs.report["input_sha256"]);
        assert_eq!(
            reference["evaluated_input_sha256"],
            dfs.report["evaluated_input_sha256"]
        );
        assert_eq!(
            reference["dimension_binding"],
            dfs.report["dimension_binding"]
        );
        for result in [&dfs, &ordered] {
            if let Some(minimum_digits) = reference["minimum_achieved_decimal_precision"].as_u64() {
                for coefficient in result.report["numerical_reference"]["coefficients"]
                    .as_array()
                    .ok_or("completed reference has no coefficient inventory")?
                {
                    for component in ["real", "imaginary"] {
                        if coefficient[format!("{component}_exact_zero")] != true {
                            assert!(
                                coefficient[format!("{component}_significant_digits")]
                                    .as_u64()
                                    .is_some_and(|digits| digits >= minimum_digits),
                                "{} epsilon^{} {component} falls below its validated {minimum_digits}-digit regression target",
                                result.report["strategy"],
                                coefficient["epsilon_power"],
                            );
                        }
                    }
                }
            }
            fixture::comparison::compare_numerical_references(
                &result.report["numerical_reference"],
                &reference["numerical_reference"],
            )?;
        }
    }
    assert!(!dfs.value.is_zero());
    Ok(())
}

#[test]
fn one_loop_zero_momentum_rung_engine_parity() -> fixture::Result<()> {
    compare(1, fixture::Projection::ZeroMomentumMetric)
}

#[test]
fn one_loop_quadratic_rung_engine_parity() -> fixture::Result<()> {
    compare(1, fixture::Projection::QuadraticTransverse)
}

#[test]
fn one_loop_rung_projections_obey_independent_dirac_and_tadpole_identity() -> fixture::Result<()> {
    use symbolica::atom::{Atom, AtomCore};
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
                strategy: RustRedReductionStrategy::Ordered,
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
#[ignore = "sequential physical-fixture qualification: enable after the preceding loop gate passes and the complete comparison measures below fifteen minutes"]
fn two_loop_zero_momentum_rung_engine_parity() -> fixture::Result<()> {
    compare(2, fixture::Projection::ZeroMomentumMetric)
}

#[test]
#[ignore = "sequential physical-fixture qualification: enable after the preceding loop gate passes and the complete comparison measures below fifteen minutes"]
fn two_loop_quadratic_rung_engine_parity() -> fixture::Result<()> {
    compare(2, fixture::Projection::QuadraticTransverse)
}

#[test]
#[ignore = "sequential physical-fixture qualification: enable after the preceding loop gate passes and the complete comparison measures below fifteen minutes"]
fn three_loop_zero_momentum_rung_engine_parity() -> fixture::Result<()> {
    compare(3, fixture::Projection::ZeroMomentumMetric)
}

#[test]
#[ignore = "sequential physical-fixture qualification: enable after the preceding loop gate passes and the complete comparison measures below fifteen minutes"]
fn three_loop_quadratic_rung_engine_parity() -> fixture::Result<()> {
    compare(3, fixture::Projection::QuadraticTransverse)
}

#[test]
#[ignore = "sequential physical-fixture qualification: requires explicit 4L data and a complete comparison below fifteen minutes"]
fn four_loop_zero_momentum_rung_engine_parity() -> fixture::Result<()> {
    compare(4, fixture::Projection::ZeroMomentumMetric)
}

#[test]
#[ignore = "sequential physical-fixture qualification: requires explicit 4L data and a complete comparison below fifteen minutes"]
fn four_loop_quadratic_rung_engine_parity() -> fixture::Result<()> {
    compare(4, fixture::Projection::QuadraticTransverse)
}

#[test]
#[ignore = "ordered-only full 5L reduction against the accepted saved zero-momentum result; requires explicit prepared input and durable evidence directories"]
fn five_loop_zero_momentum_ordered_matches_saved_depth_first_result() -> fixture::Result<()> {
    let completed = ordered_five_loop(fixture::Projection::ZeroMomentumMetric)?;
    let expected_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples/data/five_loop_q0/q0-retained.expr");
    let payload = std::fs::read(expected_path)?;
    use sha2::Digest;
    assert_eq!(payload.len(), 6_683_821);
    assert_eq!(
        format!("{:x}", sha2::Sha256::digest(&payload)),
        "da909dc0b8f199a9ce46c27368bd2d1c88b277d72e765ab70f9a83680caa1212"
    );
    let expected =
        symbolica::atom::Atom::parse(&String::from_utf8(payload)?, "vakint", Default::default())?;
    let (expected, _) =
        fixture::comparison::convert_historical_retained(&expected, &completed.report)?;
    use symbolica::atom::AtomCore;
    assert!(
        (fixture::comparison::canonical_retained(&completed.value)?
            - fixture::comparison::canonical_retained(&expected)?)
        .together()
        .expand()
        .is_zero()
    );
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
        use symbolica::atom::AtomCore;
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
            strategy: RustRedReductionStrategy::Ordered,
            artifacts: std::env::var_os("VAKINT_TEST_EXTERNAL_DATA").map(Into::into),
            monitor: RustRedMonitorOptions {
                progress_snapshot: Some(evidence.join("monitor.json")),
                ..Default::default()
            },
            numerical_reference: false,
        },
    )
}
