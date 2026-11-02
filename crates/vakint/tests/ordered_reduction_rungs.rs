#![cfg(feature = "feynkit-ingress")]
//! Physical graph/tensor ingress, lower-loop engine parity and numerical references.

use ::vakint;
#[path = "fixtures/qcd_ladder/mod.rs"]
#[allow(dead_code)]
mod fixture;
#[path = "support/reference_table.rs"]
#[allow(dead_code)]
mod reference_table;

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
/// high-precision numerical reference. The rungs are the shipped-cache
/// workload, so they ignore shipped caches (they run cold and validate the
/// step computation) unless `VAKINT_ORDERED_CACHE_SHIPPED=on`, and they fail
/// when they replay step records saved in a `VAKINT_ORDERED_CACHE_DIR`
/// (`cache_policy.previous_run_cache_restored`).
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
        // Shipped caches are replayed only on request.
        assert!(
            result.report["cache_policy"]["prepared_cache_loaded"] == false
                || result.report["ordered_step_cache"]["shipped_setting"] == true
        );
        // Saved records of a cache directory are never replayed here: they
        // would skip the step computation this test validates.
        assert_eq!(
            result.report["cache_policy"]["previous_run_cache_restored"], false,
            "the rung tests validate cold: they replayed step records of an existing \
             VAKINT_ORDERED_CACHE_DIR cache (an earlier run's, or the exact arm's for the \
             modular arm); unset VAKINT_ORDERED_CACHE_DIR"
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
    print_rung_table(&case, &exact, &modular, &reference)
}

/// The physics of a rung case in words.
fn describe_rung(loop_count: u64, projection: &serde_json::Value) -> String {
    let projection = match projection.as_str() {
        Some("zero-momentum-metric") => "zero-momentum metric projection",
        Some("quadratic-transverse") => "quadratic transverse projection (complete q² coefficient)",
        _ => "projection",
    };
    format!(
        "{loop_count}L gluon self-energy (massless d-quark loop, {} gluon rung{}), {projection}",
        loop_count - 1,
        if loop_count == 2 { "" } else { "s" }
    )
}

/// Print the rung's Laurent coefficients next to the stored golden. The exact
/// checks above are the gate; this table repeats the numerical agreement in
/// readable form (it fails only if the values disagree beyond the golden's
/// stored digits).
fn print_rung_table(
    case: &str,
    exact: &fixture::comparison::CompletedRun,
    modular: &fixture::comparison::CompletedRun,
    reference: &serde_json::Value,
) -> fixture::Result<()> {
    let rows = |numerical: &serde_json::Value| -> fixture::Result<Vec<(i64, String, String)>> {
        numerical["coefficients"]
            .as_array()
            .ok_or("numerical reference has no coefficient inventory")?
            .iter()
            .map(|row| {
                let text = |field: &str| row[field].as_str().unwrap_or("0").to_owned();
                Ok((
                    row["epsilon_power"]
                        .as_i64()
                        .ok_or("numerical reference has an invalid epsilon order")?,
                    text("real"),
                    text("imaginary"),
                ))
            })
            .collect()
    };
    let values = rows(&exact.report["numerical_reference"])?
        .into_iter()
        .map(|(power, real, imaginary)| {
            Ok((
                power,
                reference_table::complex_from_decimal(&real, &imaginary)?,
            ))
        })
        .collect::<std::result::Result<Vec<_>, String>>()?;
    let digits = reference["minimum_achieved_decimal_precision"]
        .as_u64()
        .unwrap_or(64)
        .min(
            exact.report["numerical_reference"]["selection"]["requested_decimal_precision"]
                .as_u64()
                .unwrap_or(64),
        );
    let report = &exact.report;
    let seconds = |value: &serde_json::Value| value.as_f64().unwrap_or(f64::NAN);
    let table = reference_table::LaurentTable {
        title: describe_rung(
            report["loop_count"].as_u64().unwrap_or(0),
            &report["projection"],
        ),
        facts: vec![
            (
                "golden".into(),
                format!("tests/resources/ordered_rungs/{case}.json ({digits} digits)"),
            ),
            (
                "point".into(),
                "MSbar, mUV = μ_r² = G = 1, d = 4 - 2ε".into(),
            ),
            (
                "vakint".into(),
                format!(
                    "RustRed Ordered, exact weights: {} retained masters (catalogue closed), \
                     reduction {:.2} s, numerical evaluation {:.2} s ({} digits)",
                    report["catalogue_closure"]["master_count"],
                    seconds(&report["retained_reduction_seconds"]),
                    seconds(&report["numerical_reference"]["seconds"]),
                    report["numerical_reference"]["selection"]["requested_decimal_precision"],
                ),
            ),
            (
                "exact".into(),
                format!(
                    "{} weights retain the identical expression; it equals the frozen depth-first result",
                    modular.report["ordered_weight_arithmetic"]
                        .as_str()
                        .unwrap_or("modular")
                ),
            ),
        ],
        engine: "vakint".into(),
        values,
        // Both sides are complete inventories; an order only one lists is zero.
        orders: None,
        shown_digits: 40,
    };
    let golden = reference_table::Reference {
        source: "stored golden".into(),
        tolerance: 10f64.powi(-(digits.saturating_sub(10) as i32)),
        coefficients: rows(&reference["numerical_reference"])?,
    };
    table.check(&[golden])?;
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

/// Slow optional five-loop fixture tests (ignored by default). Invocation,
/// from the workspace root (release build; 5L zero takes minutes, 5L quadratic
/// about half an hour on 16 cores):
///
/// ```text
/// VAKINT_TEST_EXTERNAL_DATA=<vakint-external-data checkout> \
/// VAKINT_QCD_LADDER_PREPARED_INPUT=<authenticated 5L scalar bundle of the projection> \
/// VAKINT_ORDERED_RUN_EVIDENCE=<fresh evidence directory> \
/// [VAKINT_N_CORES=16] [VAKINT_ORDERED_WEIGHT_ARITHMETIC=modular_polynomial] \
/// cargo test --release -p vakint --features feynkit-ingress,experimental-rustred \
///     --test ordered_reduction_rungs <test name> -- --ignored --nocapture
/// ```
///
/// The runs are cold: shipped Ordered caches (whose workload includes these
/// inputs) are ignored unless `VAKINT_ORDERED_CACHE_SHIPPED=on`, which only
/// replays recorded steps (for timings; the report then sets
/// `cache_policy.prepared_cache_loaded`). A `VAKINT_ORDERED_CACHE_DIR` that
/// already holds records is replayed as well (for timings; the report then
/// sets `cache_policy.previous_run_cache_restored`). Validate code changes
/// cold.
///
/// The zero-momentum test compares with the accepted historical result shipped
/// in `examples/data/five_loop_q0`. The quadratic test compares the retained
/// expression and its nonzero conditions with saved digests, which hold only
/// for the five-loop pack they were produced with (see
/// [`FIVE_LOOP_QUADRATIC_FIXTURE`]); with another pack the test fails and says
/// how to regenerate them.
#[test]
#[ignore = "slow optional fixture: full ordered 5L zero-momentum reduction against the accepted saved result; needs VAKINT_TEST_EXTERNAL_DATA, VAKINT_QCD_LADDER_PREPARED_INPUT and a fresh VAKINT_ORDERED_RUN_EVIDENCE (see the invocation above)"]
fn five_loop_zero_momentum_ordered_matches_saved_depth_first_result() -> fixture::Result<()> {
    let completed = ordered_five_loop(fixture::Projection::ZeroMomentumMetric)?;
    let expected_path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples/data/five_loop_q0/q0-retained.expr");
    let payload = std::fs::read(expected_path)?;
    let evidence = std::path::PathBuf::from(
        std::env::var_os("VAKINT_ORDERED_RUN_EVIDENCE")
            .ok_or("set VAKINT_ORDERED_RUN_EVIDENCE to the completed run directory")?,
    );
    let certificate = fixture::comparison::historical::compare_saved_zero_momentum(
        &payload, &completed, &evidence,
    )?;
    print_five_loop_zero_momentum_summary(&completed, &certificate);
    Ok(())
}

/// Header lines shared by the five-loop summaries.
fn five_loop_facts(completed: &fixture::comparison::CompletedRun) -> Vec<(String, String)> {
    let report = &completed.report;
    vec![
        (
            "point".into(),
            "MSbar, mUV = μ_r² = G = 1, d = 4 - 2ε; five-loop masters stay symbolic".into(),
        ),
        (
            "vakint".into(),
            format!(
                "RustRed Ordered ({} weights, {} cores): data loading {:.1} s, reduction {:.1} s; \
                 retained {} bytes with {} nonzero conditions",
                report["ordered_weight_arithmetic"].as_str().unwrap_or("?"),
                report["execution"]["n_cores"],
                report["loading_seconds"].as_f64().unwrap_or(f64::NAN),
                report["retained_reduction_seconds"]
                    .as_f64()
                    .unwrap_or(f64::NAN),
                report["retained"]["bytes"],
                report["guards"]["count"],
            ),
        ),
    ]
}

/// `M3 = 5L master, family 3f9a01c2, powers (1,1,…)` lines.
fn master_legend(masters: &serde_json::Value) -> Vec<String> {
    masters
        .as_array()
        .map(|masters| {
            masters
                .iter()
                .map(|master| {
                    format!(
                        "{} = {}L master, family {}, powers ({})",
                        master["name"].as_str().unwrap_or("?"),
                        master["loops"],
                        master["family"].as_str().unwrap_or("?"),
                        master["powers"].as_str().unwrap_or("?"),
                    )
                })
                .collect()
        })
        .unwrap_or_default()
}

fn laurent_cell(value: &serde_json::Value) -> String {
    value.as_str().unwrap_or("0").to_owned()
}

/// Print every master-product coefficient of the five-loop zero-momentum
/// result next to the accepted historical one (already compared exactly).
fn print_five_loop_zero_momentum_summary(
    completed: &fixture::comparison::CompletedRun,
    certificate: &serde_json::Value,
) {
    let rows = certificate["products"]
        .as_array()
        .map(|products| {
            products
                .iter()
                .map(|product| reference_table::ProductRow {
                    label: product["label"].as_str().unwrap_or("?").to_owned(),
                    vakint: laurent_cell(&product["ordered_laurent"]),
                    reference: Some(format!(
                        "{}   (Δ = 0: equal rational functions of d)",
                        laurent_cell(&product["historical_laurent"])
                    )),
                })
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    let mut facts = five_loop_facts(completed);
    facts.push((
        "reference".into(),
        "accepted historical depth-first result (examples/data/five_loop_q0/q0-retained.expr, \
         pySecDec normalization converted exactly to MSbar)"
            .into(),
    ));
    reference_table::print_master_products(
        &describe_rung(5, &completed.report["projection"]),
        &facts,
        &rows,
        &master_legend(&certificate["masters"]),
        &format!(
            "PASS  all {} master-product coefficients ({} masters) equal the reference exactly; \
             comparison {:.1} s",
            certificate["master_product_count"],
            certificate["master_count"],
            certificate["comparison_seconds"]
                .as_f64()
                .unwrap_or(f64::NAN),
        ),
    );
}

/// The saved 5L quadratic result: digests of the canonical retained
/// expression and of its nonzero conditions, for one five-loop reduction
/// payload and one prepared input.
struct FiveLoopFixture {
    /// Data release the digests were produced with.
    data: &'static str,
    /// sha256 of the pack's `5L/ibps/application.json`, which pins every
    /// reduction payload of the saved application by BLAKE3.
    application_sha256: &'static str,
    /// [`five_loop_reduction_payload_sha256`] of that pack: the declared
    /// sha256 of every file its 5L manifest's `ibps.files` names. Derived
    /// files a pack may add or change (shipped Ordered caches, the route
    /// preparation add-on, documentation) do not enter it, so the gate holds
    /// for every pack with the same reduction.
    reduction_payload_sha256: &'static str,
    /// sha256 of the canonical prepared input.
    input_sha256: &'static str,
    /// sha256 of `retained.symbolica`.
    retained_sha256: &'static str,
    /// sha256 of `nonzero-conditions.symbolica`.
    nonzero_conditions_sha256: &'static str,
}

/// Produced with vakint-external-data v0.1.2 (commit f1626f67) and the prepared
/// quadratic bundle `prepared-input-quadratic-4bd7448-v1`; cross-checked by
/// independent `modular_sampled` and `modular_polynomial` runs (byte-identical
/// outputs). To regenerate it for another data release: run the test cold
/// (without `VAKINT_ORDERED_CACHE_SHIPPED=on` or a warm cache directory) with
/// the new data and a fresh evidence directory under two different weight
/// arithmetics (e.g. `modular_polynomial` and `modular_sampled`), check that
/// both `retained.symbolica` and `nonzero-conditions.symbolica` agree byte for
/// byte, then replace the release name and all five digests here (the test's
/// gate error prints the two payload digests of the data checkout; the others
/// are the `report.json` fields `input_sha256`, `retained.sha256` and
/// `guards.payload.sha256`). Data that only add or change derived files
/// (ordered caches, the route preparation add-on) keep the fixture.
const FIVE_LOOP_QUADRATIC_FIXTURE: FiveLoopFixture = FiveLoopFixture {
    data: "vakint-external-data v0.1.2 (f1626f67)",
    application_sha256: "976d574c4ad8b4021ac7a1914c624095a77f5975b98e243c7a20a9396a86ff03",
    reduction_payload_sha256: "aadd2d6b8521a8b0000119803ea00e96324acc1d824a80f43b4c880011613071",
    input_sha256: "276d17a136aa5292b82effe05354d0b7cb639966a5f3797f960bbdf9f301a1fa",
    retained_sha256: "093086472ecdc75fd7180f5ce08a758f822a3048a3661eb8cc60c06e07ce34a4",
    nonzero_conditions_sha256: "5a56684173706ff69d9991147666a47fe04fff88beb888d67f006e45041a7be3",
};

fn sha256_hex(bytes: &[u8]) -> String {
    use sha2::Digest;
    format!("{:x}", sha2::Sha256::digest(bytes))
}

/// The 5L reduction payload of a data checkout: sha256 over the lines
/// `<path> <declared sha256>\n`, sorted by path, of every file that its
/// `5L/manifest.json` lists in `ibps.files` (the saved application with all
/// its payloads, queries and selection); returns it with the sha256 of the
/// `ibps/application.json` bytes. Loading the data verifies every listed file
/// against its declared sha256.
fn five_loop_reduction_payload_sha256(data: &std::path::Path) -> fixture::Result<(String, String)> {
    let manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(data.join("5L/manifest.json"))?)?;
    let mut files = manifest["ibps"]["files"]
        .as_array()
        .ok_or("5L/manifest.json has no ibps.files")?
        .iter()
        .map(|path| {
            let path = path.as_str().ok_or("ibps.files holds a non-string")?;
            let sha256 = manifest["files"][path]["sha256"]
                .as_str()
                .ok_or_else(|| format!("5L/manifest.json does not list {path}"))?;
            Ok(format!("{path} {sha256}\n"))
        })
        .collect::<fixture::Result<Vec<_>>>()?;
    files.sort();
    Ok((
        sha256_hex(files.concat().as_bytes()),
        sha256_hex(&std::fs::read(data.join("5L/ibps/application.json"))?),
    ))
}

/// Derived files (shipped Ordered caches, the route preparation add-on) may
/// be added to the 5L pack without invalidating the quadratic fixture; any
/// change of a reduction payload invalidates it.
#[test]
fn five_loop_fixture_gate_pins_the_reduction_payload_only() -> fixture::Result<()> {
    let data = tempfile::tempdir()?;
    std::fs::create_dir_all(data.path().join("5L/ibps"))?;
    let application = b"{\"schema\": \"saved application\"}";
    std::fs::write(data.path().join("5L/ibps/application.json"), application)?;
    let mut manifest = serde_json::json!({
        "ibps": {"kind": "saved-application", "path": "ibps",
                 "files": ["ibps/application.json", "ibps/objects/a.rrbin.gz.part"]},
        "files": {
            "ibps/application.json": {"sha256": sha256_hex(application), "bytes": application.len()},
            "ibps/objects/a.rrbin.gz.part": {"sha256": "a".repeat(64), "bytes": 1},
            "coverage.json": {"sha256": "c".repeat(64), "bytes": 1},
        },
    });
    let write = |manifest: &serde_json::Value| {
        std::fs::write(
            data.path().join("5L/manifest.json"),
            serde_json::to_vec_pretty(manifest).unwrap(),
        )
    };
    write(&manifest)?;
    let gate = five_loop_reduction_payload_sha256(data.path())?;
    assert_eq!(gate.1, sha256_hex(application));
    for derived in [
        "route-preparation.rrbin",
        "ordered-cache.json",
        "ordered-cache/5L-x.0.0/ordered-step-cache.txt",
    ] {
        manifest["files"][derived] = serde_json::json!({"sha256": "d".repeat(64), "bytes": 1});
    }
    manifest["route_preparation"] = serde_json::json!({"path": "route-preparation.rrbin"});
    manifest["ordered_cache"] =
        serde_json::json!({"path": "ordered-cache", "metadata": "ordered-cache.json"});
    write(&manifest)?;
    assert_eq!(five_loop_reduction_payload_sha256(data.path())?, gate);
    manifest["files"]["ibps/objects/a.rrbin.gz.part"]["sha256"] = "b".repeat(64).into();
    write(&manifest)?;
    assert_ne!(five_loop_reduction_payload_sha256(data.path())?.0, gate.0);
    Ok(())
}

#[test]
#[ignore = "slow optional fixture: full ordered 5L quadratic reduction against the saved retained and nonzero-condition digests (data v0.1.2); needs VAKINT_TEST_EXTERNAL_DATA, VAKINT_QCD_LADDER_PREPARED_INPUT and a fresh VAKINT_ORDERED_RUN_EVIDENCE (see the invocation above)"]
fn five_loop_quadratic_ordered_retains_declared_masters() -> fixture::Result<()> {
    let fixture_data = &FIVE_LOOP_QUADRATIC_FIXTURE;
    let data = std::path::PathBuf::from(
        std::env::var_os("VAKINT_TEST_EXTERNAL_DATA")
            .ok_or("set VAKINT_TEST_EXTERNAL_DATA to the vakint-external-data checkout")?,
    );
    // Fail before the long run: an ignored test only runs on request, and a
    // silent skip would read as a pass.
    let (payload, application) = five_loop_reduction_payload_sha256(&data)?;
    if payload != fixture_data.reduction_payload_sha256
        || application != fixture_data.application_sha256
    {
        return Err(format!(
            "the saved digests belong to the 5L reduction payload of {} (payload sha256 {}, \
             ibps/application.json sha256 {}), but {} has payload sha256 {payload} and \
             ibps/application.json sha256 {application}. Run against that data, or regenerate \
             the fixture as documented at FIVE_LOOP_QUADRATIC_FIXTURE in \
             tests/ordered_reduction_rungs.rs.",
            fixture_data.data,
            fixture_data.reduction_payload_sha256,
            fixture_data.application_sha256,
            data.display()
        )
        .into());
    }
    let completed = ordered_five_loop(fixture::Projection::QuadraticTransverse)?;
    assert_eq!(
        completed.report["input_sha256"], fixture_data.input_sha256,
        "the saved digests belong to a different prepared quadratic input"
    );
    let evidence = std::path::PathBuf::from(
        std::env::var_os("VAKINT_ORDERED_RUN_EVIDENCE")
            .ok_or("set VAKINT_ORDERED_RUN_EVIDENCE to the completed run directory")?,
    );
    for (file, report, expected) in [
        (
            "retained.symbolica",
            &completed.report["retained"]["sha256"],
            fixture_data.retained_sha256,
        ),
        (
            "nonzero-conditions.symbolica",
            &completed.report["guards"]["payload"]["sha256"],
            fixture_data.nonzero_conditions_sha256,
        ),
    ] {
        assert_eq!(report, expected, "{file}: report digest");
        assert_eq!(
            sha256_hex(&std::fs::read(evidence.join(file))?),
            expected,
            "{file}: written payload"
        );
    }
    print_five_loop_quadratic_summary(&completed, fixture_data);
    Ok(())
}

/// Print the five-loop quadratic result: the saved-digest checks (its only
/// reference) and the leading Laurent terms of every master-product
/// coefficient. The summary is display only; the digests above are the gate.
fn print_five_loop_quadratic_summary(
    completed: &fixture::comparison::CompletedRun,
    fixture_data: &FiveLoopFixture,
) {
    let mut facts = five_loop_facts(completed);
    facts.push((
        "reference".into(),
        format!(
            "saved digests from {} (no independent five-loop quadratic result exists)",
            fixture_data.data
        ),
    ));
    for (file, digest) in [
        ("retained.symbolica", fixture_data.retained_sha256),
        (
            "nonzero-conditions.symbolica",
            fixture_data.nonzero_conditions_sha256,
        ),
    ] {
        facts.push((
            "digest".into(),
            format!("{file} sha256 {digest} = saved: ok"),
        ));
    }
    let (rows, legend, counts) =
        match fixture::comparison::historical::master_product_summary(completed) {
            Ok(summary) => (
                summary["products"]
                    .as_array()
                    .map(|products| {
                        products
                            .iter()
                            .map(|product| reference_table::ProductRow {
                                label: product["label"].as_str().unwrap_or("?").to_owned(),
                                vakint: laurent_cell(&product["laurent"]),
                                reference: None,
                            })
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default(),
                master_legend(&summary["masters"]),
                format!(
                    "{} master products of {} masters",
                    summary["master_product_count"], summary["master_count"]
                ),
            ),
            Err(error) => (
                Vec::new(),
                Vec::new(),
                format!("master-product summary unavailable: {error}"),
            ),
        };
    reference_table::print_master_products(
        &describe_rung(5, &completed.report["projection"]),
        &facts,
        &rows,
        &legend,
        &format!(
            "PASS  retained expression and nonzero conditions reproduce the saved digests; {counts}"
        ),
    );
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
