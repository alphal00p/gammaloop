use super::*;

fn numerical(center: &str) -> serde_json::Value {
    serde_json::json!({
        "point":{"mUV":1,"mu_r_squared":1,"UFO_G":1},"epsilon_range":"test exact coefficient",
        "coefficients":[{"epsilon_power":0,
            "real_exact_zero":false,"real_significant_digits":100,"real_exact_rational":center,"real_error_exact_rational":"0",
            "imaginary_exact_zero":true,"imaginary_exact_rational":"0","imaginary_error_exact_rational":"0"}]
    })
}

fn save(directory: &Path, strategy: &str, value: &str, numerical: serde_json::Value) -> Result<()> {
    Vakint::initialize_vakint_symbols();
    std::fs::create_dir_all(directory)?;
    let value = vakint_parse!(value)?.to_canonical_string();
    let memoization = strategy == "DepthFirst";
    let report = serde_json::json!({
        "schema":"vakint.ordered-rung-comparison.v1","persistent_cache":false,"reduction_completed":true,
        "catalogue_closure":{"closed":true,"unknown_count":0},"strategy":strategy,"intrinsic_memoization":memoization,
        "cache_policy":{"depth_first_intrinsic_memoization":memoization,"prepared_cache_loaded":false,"previous_run_cache_restored":false},
        "retained":{"file":"retained.symbolica","bytes":value.len(),"sha256":sha256(value.as_bytes())},
        "loop_count":4,"projection":"quadratic-transverse","normalization":"MSbar","data":{"version":"fixture"},
        "input_sha256":"fixture-input","evaluated_input_sha256":"fixture-lowered-input",
        "dimension_binding":{"source":"d","target":"4-2*epsilon"},"numerical_reference":numerical,
    });
    std::fs::write(directory.join("retained.symbolica"), value)?;
    std::fs::write(
        directory.join("report.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    Ok(())
}

fn change_report(directory: &Path, update: impl FnOnce(&mut serde_json::Value)) -> Result<()> {
    let path = directory.join("report.json");
    let mut report = serde_json::from_slice(&std::fs::read(&path)?)?;
    update(&mut report);
    std::fs::write(path, serde_json::to_vec_pretty(&report)?)?;
    Ok(())
}

#[test]
fn saved_pair_proves_nontrivial_exact_rational_parity_without_evaluation() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let dfs = directory.path().join("dfs");
    let ordered = directory.path().join("ordered");
    save(
        &dfs,
        "DepthFirst",
        "RustRedMaster(family,1)*(1/ε+1/(1-ε))",
        numerical("1"),
    )?;
    save(
        &ordered,
        "Ordered",
        "RustRedMaster(family,1)/(ε*(1-ε))",
        numerical("1"),
    )?;
    let destination = directory.path().join("proof.json");
    compare_two_reports(&dfs, &ordered, &destination, true)?;
    let proof: serde_json::Value = serde_json::from_slice(&std::fs::read(destination)?)?;
    assert_eq!(proof["schema"], "vakint.ordered-rung-functional-pair.v1");
    assert_eq!(proof["performance_eligible"], false);
    assert_eq!(proof["numerical_enclosures_agree"], true);
    assert_eq!(proof["new_evaluations"], 0);
    assert!(proof["reports"].get("dfs_no_memo").is_none());
    Ok(())
}

#[test]
fn asymmetric_numerical_reference_requires_explicit_retain_only_request() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let dfs = directory.path().join("dfs");
    let ordered = directory.path().join("ordered");
    save(
        &dfs,
        "DepthFirst",
        "RustRedMaster(family,1)",
        serde_json::Value::Null,
    )?;
    save(
        &ordered,
        "Ordered",
        "RustRedMaster(family,1)",
        numerical("1"),
    )?;
    let originals = [
        std::fs::read(dfs.join("report.json"))?,
        std::fs::read(ordered.join("report.json"))?,
    ];
    let destination = directory.path().join("proof.json");
    assert!(compare_two_reports(&dfs, &ordered, &destination, true).is_err());
    assert!(!destination.exists());
    compare_two_reports(&dfs, &ordered, &destination, false)?;
    let proof: serde_json::Value = serde_json::from_slice(&std::fs::read(destination)?)?;
    assert!(proof["numerical_enclosures_agree"].is_null());
    assert_eq!(proof["numerical_check"], "not checked");
    assert_eq!(proof["retain_only_explicitly_requested"], true);
    assert_eq!(std::fs::read(dfs.join("report.json"))?, originals[0]);
    assert_eq!(std::fs::read(ordered.join("report.json"))?, originals[1]);
    Ok(())
}

#[test]
fn pair_rejects_corruption_and_malformed_reports_before_proof() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let dfs = directory.path().join("dfs");
    let ordered = directory.path().join("ordered");
    save(&dfs, "DepthFirst", "1", serde_json::Value::Null)?;
    save(&ordered, "Ordered", "1", serde_json::Value::Null)?;
    let destination = directory.path().join("proof.json");
    std::fs::write(ordered.join("retained.symbolica"), "2")?;
    assert!(compare_two_reports(&dfs, &ordered, &destination, false).is_err());
    std::fs::write(ordered.join("report.json"), "{")?;
    assert!(compare_two_reports(&dfs, &ordered, &destination, false).is_err());
    assert!(!destination.exists());
    Ok(())
}

#[test]
fn two_absent_numerical_references_also_require_an_explicit_opt_out() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let dfs = directory.path().join("dfs");
    let ordered = directory.path().join("ordered");
    save(&dfs, "DepthFirst", "1", serde_json::Value::Null)?;
    save(&ordered, "Ordered", "1", serde_json::Value::Null)?;
    let destination = directory.path().join("proof.json");
    assert!(compare_two_reports(&dfs, &ordered, &destination, true).is_err());
    assert!(!destination.exists());
    compare_two_reports(&dfs, &ordered, &destination, false)?;
    Ok(())
}

#[test]
fn pair_preserves_strategy_cache_closure_and_metadata_checks() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let dfs = directory.path().join("dfs");
    let ordered = directory.path().join("ordered");
    save(&dfs, "DepthFirst", "1", serde_json::Value::Null)?;
    let destination = directory.path().join("proof.json");
    for (field, value) in [
        (
            "strategy",
            serde_json::json!("DepthFirstWithoutMemoization"),
        ),
        ("intrinsic_memoization", serde_json::json!(true)),
        ("reduction_completed", serde_json::json!(false)),
        ("input_sha256", serde_json::json!("different")),
        ("data", serde_json::json!({"version":"different"})),
        ("normalization", serde_json::json!("pySecDec")),
    ] {
        save(&ordered, "Ordered", "1", serde_json::Value::Null)?;
        change_report(&ordered, |report| report[field] = value)?;
        assert!(
            compare_two_reports(&dfs, &ordered, &destination, false).is_err(),
            "{field}"
        );
        assert!(!destination.exists());
    }
    save(&ordered, "Ordered", "1", serde_json::Value::Null)?;
    change_report(&ordered, |report| {
        report["cache_policy"]["prepared_cache_loaded"] = true.into()
    })?;
    assert!(compare_two_reports(&dfs, &ordered, &destination, false).is_err());
    Ok(())
}

#[test]
fn retain_only_never_ignores_an_exact_expression_mismatch() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let dfs = directory.path().join("dfs");
    let ordered = directory.path().join("ordered");
    save(
        &dfs,
        "DepthFirst",
        "RustRedMaster(family,1)/(1-ε)",
        serde_json::Value::Null,
    )?;
    save(
        &ordered,
        "Ordered",
        "RustRedMaster(family,1)",
        numerical("1"),
    )?;
    let destination = directory.path().join("proof.json");
    assert!(compare_two_reports(&dfs, &ordered, &destination, false).is_err());
    assert!(!destination.exists());
    Ok(())
}

#[test]
fn default_pair_checks_existing_numerics_and_three_report_contract_is_unchanged() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let dfs = directory.path().join("dfs");
    let no_memo = directory.path().join("no-memo");
    let ordered = directory.path().join("ordered");
    save(&dfs, "DepthFirst", "1", numerical("1"))?;
    save(
        &no_memo,
        "DepthFirstWithoutMemoization",
        "1",
        serde_json::Value::Null,
    )?;
    save(&ordered, "Ordered", "1", numerical("2"))?;
    let destination = directory.path().join("proof.json");
    assert!(compare_two_reports(&dfs, &ordered, &destination, true).is_err());
    assert!(compare_reports(&dfs, &no_memo, &ordered, &destination).is_err());
    assert!(!destination.exists());
    Ok(())
}
