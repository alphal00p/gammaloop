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
fn archived_builtin_euler_alias_cancels_current_normalization_exactly() -> Result<()> {
    let _ = report(5)?;
    let context = CoefficientContext::try_new(["d"])?;
    let epsilon = vakint_parse!("ε")?;
    let archived_euler = Atom::var(symbolica::transcendental::euler_gamma());
    let current_euler = vakint_parse!("EulerGamma")?;
    let exact = vakint_parse!("RustRedMaster(family,1)*(2+ε)/(ε^6*(1-ε)^3)")?;
    let archived = &exact
        * Atom::var(symbolica::atom::Symbol::E)
            .pow(Atom::num(5) * (current_euler - archived_euler) * &epsilon);
    let symbols = masters(&exact);
    let actual = coefficients(&archived, &symbols, &Atom::one(), &epsilon, &context, 1)?;
    let expected = coefficients(&exact, &symbols, &Atom::one(), &epsilon, &context, 1)?;
    assert_eq!(actual.len(), 1);
    assert_eq!(actual, expected);
    assert!(
        actual
            .values()
            .all(|coefficient| context.contains(coefficient))
    );
    Ok(())
}

#[test]
fn foreign_euler_homonyms_remain_outside_exact_coefficient_context() -> Result<()> {
    let _ = report(5)?;
    let context = CoefficientContext::try_new(["d"])?;
    let epsilon = vakint_parse!("ε")?;
    let current_euler = vakint_parse!("EulerGamma")?;
    let foreign_euler = Atom::var(symbolica::symbol!("historical_comparison_control::γ"));
    let master = vakint_parse!("RustRedMaster(family,1)")?;
    let value = &master
        * Atom::var(symbolica::atom::Symbol::E)
            .pow(Atom::num(5) * (current_euler - foreign_euler) * &epsilon);
    let error = coefficients(
        &value,
        &masters(&master),
        &Atom::one(),
        &epsilon,
        &context,
        1,
    )
    .err()
    .ok_or("an unrelated Euler-named symbol was silently identified")?;
    let diagnostic = error.to_string();
    assert!(diagnostic.contains("foreign variable map"), "{diagnostic}");
    assert!(
        diagnostic.contains("historical_comparison_control"),
        "{diagnostic}"
    );
    assert!(diagnostic.contains("RustRedMaster"), "{diagnostic}");
    Ok(())
}

#[test]
fn product_coefficients_and_their_first_error_do_not_depend_on_threads() -> Result<()> {
    let _ = report(5)?;
    let context = CoefficientContext::try_new(["d"])?;
    let epsilon = vakint_parse!("ε")?;
    let value = vakint_parse!(
        "RustRedMaster(family,1)*(2+ε)/(ε^6*(1-ε)^3)+RustRedMaster(family,2)^2*(1-ε)/ε\
         +RustRedMaster(family,1)*RustRedMaster(family,2)*(3+ε^2)/(2-ε)\
         +RustRedMaster(family,3)*(1/ε-1/(1+ε))+RustRedMaster(family,1)*ε^2"
    )?;
    let symbols = masters(&value);
    let sequential = coefficients(&value, &symbols, &Atom::one(), &epsilon, &context, 1)?;
    assert_eq!(sequential.len(), 4);
    for threads in [2, 3, 8] {
        assert_eq!(
            coefficients(&value, &symbols, &Atom::one(), &epsilon, &context, threads)?,
            sequential
        );
    }
    // Two invalid products: every thread count reports the first in collection order.
    let invalid =
        &value + vakint_parse!("RustRedMaster(family,4)*sin(ε)+RustRedMaster(family,5)*cos(ε)")?;
    let symbols = masters(&invalid);
    let expected = coefficients(&invalid, &symbols, &Atom::one(), &epsilon, &context, 1)
        .err()
        .ok_or("a nonrational product was accepted")?
        .to_string();
    for threads in [2, 3, 8] {
        let error = coefficients(
            &invalid,
            &symbols,
            &Atom::one(),
            &epsilon,
            &context,
            threads,
        )
        .err()
        .ok_or("a nonrational product was accepted")?;
        assert_eq!(error.to_string(), expected);
    }
    Ok(())
}

fn direct_rational(scalar: &Atom, context: &CoefficientContext) -> Coefficient {
    scalar
        .as_view()
        .try_to_rational_polynomial(&Q, &Z, Some(context.one().get_variables().clone()))
        .unwrap()
}

#[test]
fn grouped_rational_conversion_matches_term_by_term_conversion() -> Result<()> {
    let _ = report(5)?;
    let context = CoefficientContext::try_new(["d"])?;
    let variables = context.one().get_variables().clone();
    let d = context.parameter("d").unwrap().to_expression();
    for source in [
        // one shared denominator, several numerator powers and rational content
        "(3*x^2-x/2+7/3)*(1-2*x)^(-3)*(x+5)^(-1)",
        // several denominators, a bare negative power and a polynomial term
        "x^4/(x^2-3)+2/(x^2-3)^2-x/((x-1)*(x+1))+x^3/x^2+5",
        // exact cancellation between groups
        "1/(x-1)-1/(x+1)-2/((x-1)*(x+1))",
        // no denominators at all, and a single term
        "x^7-3*x+1",
        "x/(1-x)^5",
    ] {
        let scalar = vakint_parse!(source)?
            .replace(vakint_parse!("x")?)
            .with(d.clone())
            .expand();
        let grouped = rational_coefficient(scalar.as_view(), &variables)?;
        assert_eq!(grouped, direct_rational(&scalar, &context), "{source}");
        assert!(context.contains(&grouped));
    }
    // A foreign symbol still reaches the direct conversion's variable map.
    let foreign = vakint_parse!("x/(1-x)+y/(1-x)")?
        .replace(vakint_parse!("x")?)
        .with(d.clone())
        .expand();
    let grouped = rational_coefficient(foreign.as_view(), &variables)?;
    assert_eq!(grouped, direct_rational(&foreign, &context));
    assert!(!context.contains(&grouped));
    Ok(())
}

#[test]
#[ignore = "authenticated 6.68 MB saved result: grouped against term-by-term Q(d) conversion of all 58 products (several minutes)"]
fn saved_zero_momentum_grouped_conversion_matches_term_by_term_conversion() -> Result<()> {
    let report = report(5)?;
    let payload = std::fs::read(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/data/five_loop_q0/q0-retained.expr"),
    )?;
    assert_eq!(sha256(&payload), HISTORICAL_SHA256);
    let original = Atom::parse(std::str::from_utf8(&payload)?, "vakint", Default::default())?;
    let (converted, _) = convert_historical_retained(&original, &report)?;
    let value = unit_point(&converted, &report)?;
    let context = CoefficientContext::try_new(["d"])?;
    let variables = context.one().get_variables().clone();
    let epsilon = vakint_parse!("ε")?;
    let physical_epsilon =
        (Atom::num(4) - context.parameter("d").unwrap().to_expression()) / Atom::num(2);
    let common = vakint_parse!("𝑖/(16*𝜋^2)")?.pow(5);
    let symbols = masters(&value).into_iter().collect::<Vec<_>>();
    let rows = value.coefficient_list::<i32>(&symbols);
    assert_eq!(rows.len(), 58);
    let threads = super::super::fixture_threads();
    let mismatches = super::super::map_ordered(threads, &rows, |(_, scalar)| {
        let scalar = (scalar / &common)
            .replace(Atom::var(symbolica::transcendental::euler_gamma()))
            .with(Atom::var(vakint_symbol!("EulerGamma")));
        let scalar = canonical_retained(&scalar)
            .unwrap()
            .expand()
            .replace(epsilon.clone())
            .with(physical_epsilon.clone());
        let started = Instant::now();
        let grouped = rational_coefficient(scalar.as_view(), &variables).unwrap();
        let grouped_seconds = started.elapsed().as_secs_f64();
        let started = Instant::now();
        let direct = direct_rational(&scalar, &context);
        println!(
            "grouped {grouped_seconds:.3} s, term by term {:.3} s",
            started.elapsed().as_secs_f64()
        );
        grouped != direct
    });
    assert_eq!(mismatches.iter().filter(|mismatch| **mismatch).count(), 0);
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
