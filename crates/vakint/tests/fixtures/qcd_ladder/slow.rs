//! Full native preparation and exact saved-application acceptance. Never cached.
use super::*;
use crate as vakint;
use symbolica::atom::AtomView;
#[path = "mod.rs"]
mod fixture;

#[test_log::test]
#[ignore = "full slow pipeline: requires VAKINT_SAVED_APPLICATION and native licensed tensor algebra"]
fn qcd_ladder_quadratic_transverse_projection_and_saved_application() -> fixture::Result<()> {
    use crate::{EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions};
    let artifact = std::env::var_os("VAKINT_SAVED_APPLICATION").ok_or("set VAKINT_SAVED_APPLICATION to the explicitly selected saved publication; no default/download is permitted")?;
    let evidence = std::env::var_os("VAKINT_QCD_LADDER_EVIDENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join(format!("../../TMP/qcd-ladder-slow-{}", std::process::id()))
        });
    let prepared_input = fixture::prepare(&evidence)?;
    eprintln!(
        "qcd_ladder phase=preparation_complete seconds={:.6} scalar_terms={} contributions={}",
        prepared_input.preparation_seconds,
        prepared_input.scalar_terms,
        prepared_input.contributions
    );
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                substitute_masters: false,
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
        "# physical quadratic Taylor terms have A-R=10; D below is a stored shifted-denominator key, not a new physical starting request\n# This census describes the reported native geometric witnesses only; a refusal here does not exclude another admitted exact witness. Full public scope-aware matching follows the census.\n# D mass_squared_power_j terms; D-j=10 is checked exactly\n",
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
                    return Err("stored key and explicit mass power violate physical logarithmic homogeneity D-j=10".into());
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
    // public evaluation still searches its complete bounded native witness
    // alternatives and enforces the saved union before any rule application.
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
