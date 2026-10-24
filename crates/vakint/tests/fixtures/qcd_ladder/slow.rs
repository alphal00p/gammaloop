//! Exact saved-application acceptance with fresh or verified saved scalar preparation.
use super::*;
use crate as vakint;
use symbolica::atom::AtomView;
#[path = "mod.rs"]
mod fixture;

#[test_log::test]
#[ignore = "full reduction: requires VAKINT_SAVED_APPLICATION; set VAKINT_QCD_LADDER_PREPARED_INPUT to reuse verified scalar preparation"]
fn qcd_ladder_quadratic_transverse_projection_and_saved_application() -> fixture::Result<()> {
    use crate::{EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions};
    let artifact = std::env::var_os("VAKINT_SAVED_APPLICATION").ok_or("set VAKINT_SAVED_APPLICATION to the explicitly selected saved publication; no default/download is permitted")?;
    let evidence = PathBuf::from(
        std::env::var_os("VAKINT_QCD_LADDER_EVIDENCE")
            .ok_or("set VAKINT_QCD_LADDER_EVIDENCE to a development workspace directory")?,
    );
    std::fs::create_dir_all(&evidence)?;
    let input_started = std::time::Instant::now();
    let force_fresh = match std::env::var_os("VAKINT_QCD_LADDER_FRESH_PREPARE") {
        None => false,
        Some(value) if value == "1" => true,
        Some(_) => return Err("VAKINT_QCD_LADDER_FRESH_PREPARE must be 1 or unset".into()),
    };
    let mut prepared_directory =
        std::env::var_os("VAKINT_QCD_LADDER_PREPARED_INPUT").map(PathBuf::from);
    if force_fresh && prepared_directory.is_some() {
        return Err("VAKINT_QCD_LADDER_FRESH_PREPARE cannot be combined with VAKINT_QCD_LADDER_PREPARED_INPUT".into());
    }
    if prepared_directory.is_none() && fixture::LadderFixture::has_prepared_input(&evidence)? {
        if force_fresh {
            return Err("fresh preparation requires a new VAKINT_QCD_LADDER_EVIDENCE directory; preserve the completed bundle".into());
        }
        prepared_directory = Some(evidence.clone());
    }
    let prepared_input = match &prepared_directory {
        Some(directory) => fixture::LadderFixture::load_prepared(
            directory,
            fixture::Projection::QuadraticTransverse,
        )?,
        None => {
            let prepared = fixture::prepare(&evidence, fixture::Projection::QuadraticTransverse)?;
            prepared.save_prepared(&evidence)?;
            prepared
        }
    };
    assert_eq!(prepared_input.contributions, 21);
    eprintln!(
        "qcd_ladder phase=input_ready input_mode={} seconds={:.6} scalar_terms={} contributions={}",
        if prepared_directory.is_some() {
            "prepared"
        } else {
            "fresh"
        },
        input_started.elapsed().as_secs_f64(),
        prepared_input.scalar_terms,
        prepared_input.contributions
    );
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: vakint::RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/qcd-ladder-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    };
    let vakint = Vakint::new()?.with_rustred_application(&mut settings, artifact)?;
    let application = vakint.saved_application(&settings)?;
    let mut expression = VakintExpression::try_from(prepared_input.input.as_view())?;
    expression.canonicalize(&settings, &vakint.topologies, false)?;
    let mut histogram = BTreeMap::<(i64, i64), usize>::new();
    let mut unsupported = BTreeSet::new();
    let mut lowered_count = 0usize;
    let mut report = String::from(
        "# unlowered quadratic source terms have A-r_mom=10; D_key below is a signed denominator-key degree\n# This census describes the reported native geometric witnesses only; a refusal here does not exclude another admitted exact witness. Full public scope-aware matching follows the census.\n# D_key coefficient_mass_squared_power_j terms; D_key-j=10 checks engineering mass dimension only, not individual-term UV degree at fixed MUV\n",
    );
    let muv = Atom::var(vk_symbol!("mUV"));
    for term in &expression.0 {
        let mut matched = vakint
            .topologies
            .match_topologies_to_user_input(term.integral.as_view(), false)?
            .ok_or("prepared ladder has no exact native topology witness")?;
        matched.apply_replacement_rules()?;
        let mut term = term.clone();
        term.apply_scalar_numerator_replacement_rules(&matched, &settings)?;
        let descriptor = select(&matched.canonical_topology)?;
        let input = MatchedScalarInput::prepare(
            term.numerator.as_view(),
            &matched,
            &descriptor.momenta,
            descriptor.physical_slots,
            15,
            &settings,
            ScalarNumeratorForm::PreserveProducts,
        )?;
        report.push_str(&format!(
            "witness topology={} target={:?}\n",
            matched.canonical_topology.get_integral().name,
            input.target.powers()
        ));
        let prepared = application
            .prepared
            .lock()
            .map_err(|error| error.to_string())?;
        let scalar = FamilyScalarNumeratorService::try_new(
            prepared.family(),
            S.dot,
            (1..=5).map(|axis| S.k.call(axis)).collect(),
            ScalarNumeratorLimits::default(),
        )?;
        let lowered =
            scalar.lower_weighted(&input.numerator, &input.target, &input.mass_squared)?;
        for term in lowered.terms() {
            lowered_count += 1;
            let degree = term.integral().powers().iter().sum::<i64>();
            let weight = term.coefficient().to_expression() * term.scalar_spectator();
            for (mass_power, coefficient) in weight.coefficient_list::<u16>(&[&muv]) {
                if coefficient.is_zero() {
                    continue;
                }
                let exponent = if mass_power.is_one() {
                    0
                } else if mass_power == muv {
                    1
                } else if let AtomView::Pow(power) = mass_power.as_view() {
                    let (base, exponent) = power.get_base_exp();
                    if base != muv.as_view() {
                        return Err("native mass collection returned another variable".into());
                    }
                    i64::try_from(exponent)?
                } else {
                    return Err("native mass collection did not return a monomial".into());
                };
                if exponent < 0 || exponent % 2 != 0 || degree - exponent / 2 != 10 {
                    return Err("stored key and explicit mass power violate engineering mass-dimension balance D_key-j=10".into());
                }
                *histogram.entry((degree, exponent / 2)).or_default() += 1;
            }
            if !prepared.routed_entry_is_admitted(term.integral())? {
                unsupported.insert(term.integral().clone());
            }
        }
    }
    for ((degree, mass_power), count) in &histogram {
        report.push_str(&format!("{degree} {mass_power} {count}\n"));
    }
    report.push_str(&format!(
        "lowered_terms={lowered_count}\nunsupported_distinct_keys={}\n",
        unsupported.len()
    ));
    for key in &unsupported {
        report.push_str(&format!("unsupported {:?}\n", key.powers()));
    }
    std::fs::write(evidence.join("saved-entry-admission.txt"), report)?;
    // The geometric census does not authorize/reject the physical input:
    // Public evaluation searches its bounded native witness alternatives.
    // Recorded coverage is advisory; individually safe rules may still reduce
    // a verified fallback, and surviving undeclared terminals are reported.
    let started = std::time::Instant::now();
    let first = vakint.evaluate_rustred_symbolic(&settings, prepared_input.input.as_view())?;
    let first_seconds = started.elapsed().as_secs_f64();
    let started = std::time::Instant::now();
    let warm = vakint.evaluate_rustred_symbolic(&settings, prepared_input.input.as_view())?;
    assert!((first.value - warm.value).together().expand().is_zero());
    assert_eq!(first.families, warm.families);
    assert_eq!(first.nonzero_conditions, warm.nonzero_conditions);
    // Preparation poles may cancel before Vakint sees its scalar input. They
    // remain explicit additional conditions of this result; evaluate() cannot
    // infer their provenance from an already simplified numerator.
    let d = rustred::algebra::CoefficientContext::try_new(["d"])?
        .parameter("d")
        .ok_or("missing native dimension")?
        .to_expression();
    let physical_d =
        Atom::num(4) - Atom::num(2) * Atom::var(vk_symbol!(settings.epsilon_symbol.as_str()));
    let prepared_conditions = prepared_input
        .nonzero_conditions
        .iter()
        .map(|condition| condition.replace(d.clone()).with(physical_d.clone()))
        .collect::<BTreeSet<_>>();
    assert!(
        prepared_conditions
            .iter()
            .all(|condition| !condition.is_zero())
    );
    let complete_conditions = first
        .nonzero_conditions
        .iter()
        .cloned()
        .chain(prepared_conditions.iter().cloned())
        .collect::<BTreeSet<_>>();
    let complete_warm_conditions = warm
        .nonzero_conditions
        .iter()
        .cloned()
        .chain(prepared_conditions)
        .collect::<BTreeSet<_>>();
    assert_eq!(complete_conditions, complete_warm_conditions);
    std::fs::write(
        evidence.join("guarded-result-conditions.symbolica"),
        complete_conditions
            .iter()
            .map(AtomCore::to_canonical_string)
            .collect::<Vec<_>>()
            .join("\n"),
    )?;
    eprintln!(
        "qcd_ladder phase=reduction_complete first_seconds={first_seconds:.6} warm_seconds={:.6} lowered_terms={lowered_count}",
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

#[test]
#[ignore = "native zero-momentum control; requires licensed tensor algebra and explicit evidence directory"]
fn qcd_ladder_zero_momentum_metric_retains_both_current_vertices() -> fixture::Result<()> {
    let evidence = PathBuf::from(
        std::env::var_os("VAKINT_QCD_LADDER_EVIDENCE")
            .ok_or("set VAKINT_QCD_LADDER_EVIDENCE to a development workspace directory")?,
    );
    let prepared = fixture::prepare(&evidence, fixture::Projection::ZeroMomentumMetric)?;
    assert!(matches!(
        prepared.projection,
        fixture::Projection::ZeroMomentumMetric
    ));
    assert_eq!(prepared.contributions, 1);
    assert!(!prepared.input.is_zero());
    assert!(
        prepared.nonzero_conditions.is_empty(),
        "the unaveraged metric contraction must not acquire transverse projector poles"
    );
    assert_eq!(prepared.timings.angular_seconds, 0.0);
    assert_eq!(prepared.timings.taylor_seconds, 0.0);
    Ok(())
}
