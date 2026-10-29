use super::*;

fn report(loops: i64) -> Result<serde_json::Value> {
    super::super::Vakint::initialize_vakint_symbols();
    let dimension = super::super::super::dimension()?;
    let physical = vakint_parse!("4-2*ε")?;
    Ok(serde_json::json!({
        "loop_count":loops,"normalization":"MSbar",
        "epsilon_symbol":"vakint::ε","mu_r_sq_symbol":"vakint::mursq",
        "dimension_binding":{"source":dimension.to_canonical_string(),"target":physical.to_canonical_string()}
    }))
}

fn pair(old_core: Atom, new_core: Atom, loops: i64) -> Result<Comparison> {
    let report = report(loops)?;
    let epsilon = vakint_parse!("ε")?;
    let gamma = vakint_parse!("EulerGamma")?;
    let coupling = Atom::var(symbolica::symbol!("UFO::GC_11")).pow(2 * loops);
    let mass = vakint_parse!("mUV")?.pow(Atom::num(2) - Atom::num(2 * loops) * &epsilon);
    let old = old_core
        * &coupling
        * &mass
        * Atom::var(symbolica::atom::Symbol::E).pow(-Atom::num(loops) * &gamma * &epsilon);
    let new = new_core
        * coupling
        * mass
        * vakint_parse!("mursq")?.pow(Atom::num(loops) * &epsilon)
        * vakint_parse!("𝑖/(16*𝜋^2)")?.pow(loops);
    Comparison::checked(&old, &new, &report)
}

#[test]
fn complete_rational_products_agree_across_different_outer_factorizations() -> Result<()> {
    let _ = report(1)?;
    let old = vakint_parse!(
        "((1+ε)*RustRedMaster(family,1)+(1-ε)*RustRedMaster(family,1)*RustRedMaster(family,2))/(ε*(1-ε))"
    )?;
    let new = vakint_parse!(
        "RustRedMaster(family,1)*(1+ε)/(ε*(1-ε))+RustRedMaster(family,1)*RustRedMaster(family,2)/ε"
    )?;
    let result = pair(old, new, 1)?;
    assert_eq!(result.master_count, 2);
    assert_eq!(result.rows.len(), 2);
    assert!(result.rows.iter().all(|row| row["exactly_equal"] == true));
    Ok(())
}

#[test]
fn exact_poles_and_canceling_master_coefficients_do_not_require_laurent_truncation() -> Result<()> {
    let _ = report(1)?;
    let old = vakint_parse!(
        "RustRedMaster(family,1)*(1/ε+1/(1-ε))+RustRedMaster(family,2)*(ε^2/(1-ε)-1/(1-ε)+1+ε)"
    )?;
    let new = vakint_parse!("RustRedMaster(family,1)/(ε*(1-ε))")?;
    let result = pair(old, new, 1)?;
    assert_eq!(result.rows.len(), 1);
    assert_eq!(result.master_count, 1);
    Ok(())
}

#[test]
fn five_loop_normalization_phase_and_authenticated_coupling_expand_exactly() -> Result<()> {
    let _ = report(5)?;
    let core = vakint_parse!("RustRedMaster(family,1)*(2+ε)/(ε^6*(1-ε)^3)")?;
    let result = pair(core.clone(), core, 5)?;
    assert_eq!(result.master_count, 1);
    assert_eq!(result.rows.len(), 1);
    assert!(!result.factor.is_zero());
    Ok(())
}

#[test]
fn a_nonzero_extra_product_is_reported_explicitly() -> Result<()> {
    let _ = report(1)?;
    let old = vakint_parse!("RustRedMaster(family,1)/ε")?;
    let new = vakint_parse!("RustRedMaster(family,1)/ε+RustRedMaster(family,2)^2")?;
    let error = pair(old, new, 1)
        .err()
        .ok_or("unequal products were accepted")?;
    assert!(error.to_string().contains("RustRedMaster"));
    assert!(error.to_string().contains("differs"));
    Ok(())
}

#[test]
fn inverse_or_nested_master_dependencies_are_not_treated_as_scalar_coefficients() -> Result<()> {
    let _ = report(1)?;
    for source in [
        "RustRedMaster(family,1)^(-1)",
        "unregistered_function(RustRedMaster(family,1))",
    ] {
        let value = vakint_parse!(source)?;
        let error = pair(value.clone(), value, 1)
            .err()
            .ok_or("a nonpolynomial master dependency was accepted")?;
        assert!(error.to_string().contains("not a polynomial"));
    }
    Ok(())
}

#[test]
fn scalar_nonrational_functions_cannot_pass_as_full_rational_proofs() -> Result<()> {
    let _ = report(1)?;
    let value = vakint_parse!("RustRedMaster(family,1)*sin(ε)")?;
    assert!(pair(value.clone(), value, 1).is_err());
    Ok(())
}

#[test]
fn historical_hash_is_checked_before_parsing_or_creating_output() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let completed = CompletedRun {
        value: Atom::one(),
        report: serde_json::Value::Null,
    };
    let error =
        compare_saved_zero_momentum(b"not the archived result", &completed, directory.path())
            .err()
            .ok_or("an unauthenticated result was accepted")?;
    assert!(error.to_string().contains("SHA256"));
    assert_eq!(std::fs::read_dir(directory.path())?.count(), 0);
    Ok(())
}

#[test]
#[ignore = "full authenticated 6.68 MB saved-result classification; no reductions, caches or numerical master evaluation"]
fn saved_zero_momentum_reference_has_58_complete_rational_products() -> Result<()> {
    let report = report(5)?;
    let payload = std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/data/five_loop_q0/q0-retained.expr"),
    )?;
    assert_eq!(payload.len(), HISTORICAL_BYTES);
    assert_eq!(sha256(&payload), HISTORICAL_SHA256);
    let original = Atom::parse(std::str::from_utf8(&payload)?, "vakint", Default::default())?;
    let (converted, _) = convert_historical_retained(&original, &report)?;
    let classified = Comparison::checked(&original, &converted, &report)?;
    assert_eq!(classified.master_count, 55);
    assert_eq!(classified.rows.len(), 58);
    Ok(())
}
