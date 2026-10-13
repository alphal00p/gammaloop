//! Manual, ignored failure diagnosis; never selected by the full slow class.
//! Acceptance still reruns the complete native QCD preparation in slow.rs.
use super::*;
use symbolica::{atom::AtomView, parser::ParseSettings, prelude::Rational};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

// Metadata-only mirror of the existing conservative term-support estimate.
// It performs no symbolic algebra; saturation bounds the diagnostic itself.
fn supports(atom: AtomView<'_>) -> (u128, u128) {
    let mul = |a: (u128, u128), b: (u128, u128)| (a.0.saturating_mul(b.0), a.1.saturating_mul(b.1));
    match atom {
        AtomView::Mul(product) => product.iter().map(supports).fold((1, 1), mul),
        AtomView::Add(sum) => sum.iter().map(supports).fold((0u128, 1u128), |a, b| {
            (
                a.0.saturating_mul(b.1)
                    .saturating_add(b.0.saturating_mul(a.1)),
                a.1.saturating_mul(b.1),
            )
        }),
        AtomView::Pow(power) => {
            let Ok(exponent) = Rational::try_from(power.get_exp()) else {
                return (1, 1);
            };
            if !exponent.is_integer() {
                return (1, 1);
            }
            let Some(exponent) = exponent.numerator().to_i64() else {
                return (u128::MAX, u128::MAX);
            };
            let mut base = supports(power.get_base());
            if exponent < 0 {
                base = (base.1, base.0);
            }
            let mut count = exponent.unsigned_abs();
            let mut value = (1, 1);
            while count != 0 {
                if count & 1 != 0 {
                    value = mul(value, base);
                }
                count >>= 1;
                if count != 0 {
                    base = mul(base, base);
                }
            }
            value
        }
        _ => (1, 1),
    }
}

fn nodes(atom: AtomView<'_>) -> usize {
    1 + match atom {
        AtomView::Mul(product) => product.iter().map(nodes).sum(),
        AtomView::Add(sum) => sum.iter().map(nodes).sum(),
        AtomView::Pow(power) => nodes(power.get_base()) + nodes(power.get_exp()),
        AtomView::Fun(function) => function.iter().map(nodes).sum(),
        _ => 0,
    }
}

fn visible_poles(atom: AtomView<'_>, poles: &mut BTreeSet<Atom>) {
    match atom {
        AtomView::Pow(power) => {
            if Rational::try_from(power.get_exp()).is_ok_and(|exponent| exponent.is_negative())
                && !matches!(power.get_base(), AtomView::Num(_))
            {
                poles.insert(power.get_base().to_owned());
            }
            visible_poles(power.get_base(), poles);
            visible_poles(power.get_exp(), poles);
        }
        AtomView::Mul(product) => {
            for factor in product.iter() {
                visible_poles(factor, poles);
            }
        }
        AtomView::Add(sum) => {
            for term in sum.iter() {
                visible_poles(term, poles);
            }
        }
        AtomView::Fun(function) => {
            for argument in function.iter() {
                visible_poles(argument, poles);
            }
        }
        _ => {}
    }
}

// Flat native-expression capture for the core-only stage diagnostic. This is
// manual evidence, not a portable application format or a cached slow fixture.
fn capture_weighted_input(
    evidence: &std::path::Path,
    index: usize,
    input: &MatchedScalarInput,
    family: &rustred::family::IntegralFamily,
) -> Result<()> {
    let directory = evidence.join(format!("input-{index:02}"));
    std::fs::create_dir_all(&directory)?;
    for (name, value) in [
        ("numerator.symbolica", input.numerator.to_canonical_string()),
        ("mass.symbolica", input.mass_squared.to_canonical_string()),
        ("family-name.txt", family.name().to_owned()),
        ("family-fingerprint.txt", family.fingerprint().to_owned()),
        ("family-loops.txt", family.loop_momenta().join("\n")),
        (
            "family-parameters.txt",
            family.coefficient_context().parameter_names().join("\n"),
        ),
        (
            "family-dimension.symbolica",
            family.dimension().to_expression().to_canonical_string(),
        ),
        ("dot-head.txt", Atom::var(S.dot).to_canonical_string()),
        (
            "loop-labels.symbolica",
            (1..=family.loop_count())
                .map(|axis| S.k.call(axis).to_canonical_string())
                .collect::<Vec<_>>()
                .join("\n"),
        ),
        (
            "powers.txt",
            input
                .target
                .powers()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(" "),
        ),
        (
            "family-shifts.symbolica",
            family
                .power_shifts()
                .iter()
                .map(|value| value.to_expression().to_canonical_string())
                .collect::<Vec<_>>()
                .join("\n"),
        ),
        (
            "family-denominators.tsv",
            family
                .denominators()
                .iter()
                .map(|denominator| {
                    std::iter::once(denominator.constant())
                        .chain(denominator.coefficients())
                        .map(|value| value.to_expression().to_canonical_string())
                        .collect::<Vec<_>>()
                        .join("\t")
                })
                .collect::<Vec<_>>()
                .join("\n"),
        ),
    ] {
        std::fs::write(directory.join(name), value)?;
    }
    Ok(())
}

#[test]
#[ignore = "manual diagnostic only; requires explicitly saved prepared input and an external saved application"]
fn inspect_weighted_ladder_resource_failure() -> Result<()> {
    use crate::{EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions};
    let source = std::env::var_os("VAKINT_QCD_PREPARED_INPUT")
        .ok_or("set VAKINT_QCD_PREPARED_INPUT for diagnostic replay only")?;
    let evidence = PathBuf::from(
        std::env::var_os("VAKINT_QCD_WEIGHT_EVIDENCE").ok_or("set VAKINT_QCD_WEIGHT_EVIDENCE")?,
    );
    std::fs::create_dir_all(&evidence)?;
    let artifact =
        std::env::var_os("VAKINT_SAVED_APPLICATION").ok_or("set the explicit saved application")?;
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                substitute_masters: false,
            },
        )]),
        form_exe_path: "/diagnostic-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    };
    let vakint = Vakint::new()?.with_rustred_application(&mut settings, artifact)?;
    let source = Atom::parse(
        &std::fs::read_to_string(source)?,
        "qcd_weight_diagnostic",
        ParseSettings::default(),
    )?;
    let mut expression = VakintExpression::try_from(source.as_view())?;
    expression.canonicalize(&settings, &vakint.topologies, false)?;
    let application = vakint.saved_application(&settings)?;
    for (input_index, term) in expression.0.iter().enumerate() {
        let mut matched = vakint
            .topologies
            .match_topologies_to_user_input(term.integral.as_view(), false)?
            .ok_or("geometric witness missing")?;
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
        let failure =
            match scalar.lower_weighted(&input.numerator, &input.target, &input.mass_squared) {
                Ok(_) => continue,
                Err(error) => error.to_string(),
            };
        let lower = scalar.lower(&input.numerator.expand(), &input.target)?;
        let mut weights = BTreeMap::<IntegralKey, Atom>::new();
        let mut parts = BTreeMap::<IntegralKey, Vec<(Atom, Atom, u32)>>::new();
        for term in lower.terms() {
            let coefficient = term.coefficient().to_expression();
            let spectator = term.scalar_spectator().clone();
            let mass = term.common_mass_squared_power();
            parts.entry(term.integral().clone()).or_default().push((
                coefficient.clone(),
                spectator.clone(),
                mass,
            ));
            let value = coefficient * spectator * input.mass_squared.clone().pow(mass);
            *weights.entry(term.integral().clone()).or_insert(Atom::Zero) += value;
        }
        let mut report = format!(
            "input_index={input_index}\ntarget={:?}\noriginal_error={failure}\nlowered_components={}\ncoalesced_integrals={}\n",
            input.target.powers(),
            lower.terms().len(),
            weights.len()
        );
        for (weight_index, (integral, weight)) in weights.into_iter().enumerate() {
            let before = supports(weight.as_view());
            if before.0 <= 100_000 && before.1 <= 100_000 {
                continue;
            }
            let components = &parts[&integral];
            let denominators = components
                .iter()
                .map(|(coefficient, _, _)| {
                    let rational: Coefficient = coefficient
                        .try_to_rational_polynomial(
                            &symbolica::prelude::Q,
                            &symbolica::prelude::Z,
                            None,
                        )
                        .unwrap();
                    rational.denominator.to_expression()
                })
                .collect::<BTreeSet<_>>();
            // Pinned Symbolica collect_factors performs only exact syntactic
            // common-factor extraction; there is no polynomial expansion here.
            let factored = weight.collect_factors();
            let after = supports(factored.as_view());
            report.push_str(&format!("weight_index={weight_index}\nintegral={:?}\ncomponent_count={}\ndistinct_coefficient_denominators={}\nbefore_nodes={}\nbefore_support_bound={before:?}\nafter_nodes={}\nafter_support_bound={after:?}\n",integral.powers(),components.len(),denominators.len(),nodes(weight.as_view()),nodes(factored.as_view())));
            std::fs::write(
                evidence.join("weight-before.symbolica"),
                weight.to_canonical_string(),
            )?;
            std::fs::write(
                evidence.join("weight-collected.symbolica"),
                factored.to_canonical_string(),
            )?;
            std::fs::write(
                evidence.join("coefficient-denominators.symbolica"),
                denominators
                    .iter()
                    .map(AtomCore::to_canonical_string)
                    .collect::<Vec<_>>()
                    .join("\n"),
            )?;
            if after.0 <= 100_000 && after.1 <= 100_000 {
                let together = factored.together();
                report.push_str(&format!(
                    "native_collected_together_nodes={}\nnative_collected_together_support={:?}\n",
                    nodes(together.as_view()),
                    supports(together.as_view())
                ));
                std::fs::write(
                    evidence.join("weight-together.symbolica"),
                    together.to_canonical_string(),
                )?;
            }
            std::fs::write(evidence.join("weight-report.txt"), report)?;
            return Ok(());
        }
        std::fs::write(evidence.join("weight-report.txt"), report)?;
        return Err("failure was not the expected per-key rational support preflight; inspect preserved report".into());
    }
    Err("prepared diagnostic input no longer reproduces the weighted failure".into())
}

#[test]
#[ignore = "manual prepared-input validation only; full slow acceptance always rebuilds the QCD numerator"]
fn prepared_ladder_weighted_preflight_without_ibps() -> Result<()> {
    use crate::{EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions};
    let mut scalar_limits = ScalarNumeratorLimits::default();
    let manual_work_override = match std::env::var("VAKINT_QCD_DIAGNOSTIC_MAX_WORK") {
        Ok(value) => Some(value),
        Err(std::env::VarError::NotPresent) => None,
        Err(error) => return Err(error.into()),
    };
    if let Some(value) = &manual_work_override {
        scalar_limits.exact_algebra.max_term_operations = value.parse()?;
    }
    println!(
        "prepared_input_diagnostic scalar_max_term_operations={} manual_work_override={}",
        scalar_limits.exact_algebra.max_term_operations,
        manual_work_override.is_some()
    );
    let source = std::env::var_os("VAKINT_QCD_PREPARED_INPUT")
        .ok_or("set the saved prepared input for diagnostic replay only")?;
    let evidence = PathBuf::from(
        std::env::var_os("VAKINT_QCD_WEIGHT_EVIDENCE")
            .ok_or("set diagnostic evidence directory")?,
    );
    std::fs::create_dir_all(&evidence)?;
    let artifact =
        std::env::var_os("VAKINT_SAVED_APPLICATION").ok_or("set explicit saved application")?;
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                substitute_masters: false,
            },
        )]),
        form_exe_path: "/diagnostic-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    };
    let vakint = Vakint::new()?.with_rustred_application(&mut settings, artifact)?;
    let source = Atom::parse(
        &std::fs::read_to_string(source)?,
        "qcd_weight_diagnostic",
        ParseSettings::default(),
    )?;
    let mut expression = VakintExpression::try_from(source.as_view())?;
    expression.canonicalize(&settings, &vakint.topologies, false)?;
    assert_eq!(
        expression.0.len(),
        21,
        "the prepared ladder must retain all Taylor contributions"
    );
    let application = vakint.saved_application(&settings)?;
    let mut report = String::from(
        "# Diagnostic replay only; no fresh numerator preparation, no IBPs.\n# Geometric witnesses only: refusals do not rule out another exact presentation.\n",
    );
    report.push_str(&format!(
        "# scalar_max_term_operations={} manual_work_override={}\n# Contribution ordinals are run-local; compare exact recorded targets/numerators.\n",
        scalar_limits.exact_algebra.max_term_operations,
        manual_work_override.is_some()
    ));
    let mut count = 0usize;
    let mut admitted = 0usize;
    let mut guards = BTreeSet::new();
    let mut histogram = BTreeMap::<(i64, i64), usize>::new();
    let muv = Atom::var(vk_symbol!("mUV"));
    for (index, term) in expression.0.iter().enumerate() {
        let mut matched = vakint
            .topologies
            .match_topologies_to_user_input(term.integral.as_view(), false)?
            .ok_or("geometric witness missing")?;
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
            super::super::numerator::ScalarNumeratorForm::PreserveProducts,
        )?;
        let prepared = application
            .prepared
            .lock()
            .map_err(|error| error.to_string())?;
        let scalar = FamilyScalarNumeratorService::try_new(
            prepared.family(),
            S.dot,
            (1..=5).map(|axis| S.k.call(axis)).collect(),
            scalar_limits,
        )?;
        capture_weighted_input(&evidence, index, &input, prepared.family())?;
        let lowering =
            scalar.lower_weighted(&input.numerator, &input.target, &input.mass_squared)?;
        let mut original_poles = BTreeSet::new();
        visible_poles(input.numerator.as_view(), &mut original_poles);
        visible_poles(input.mass_squared.as_view(), &mut original_poles);
        for pole in original_poles {
            assert!(
                lowering.nonzero_conditions().contains(&pole),
                "original pole lost: {pole}"
            );
        }
        assert!(
            lowering
                .nonzero_conditions()
                .iter()
                .all(|condition| !condition.contains_symbol(S.k)),
            "a local denominator-coordinate label escaped in a guard"
        );
        report.push_str(&format!(
            "input={index} topology={} target={:?} lowered={}\n",
            matched.canonical_topology.get_integral().name,
            input.target.powers(),
            lowering.terms().len()
        ));
        guards.extend(lowering.nonzero_conditions().iter().cloned());
        for term in lowering.terms() {
            assert!(!term.scalar_spectator().contains_symbol(S.k));
            assert!(!term.coefficient().to_expression().contains_symbol(S.k));
            count += 1;
            let accepted = prepared.routed_entry_is_admitted(term.integral())?;
            admitted += usize::from(accepted);
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
                        return Err("unexpected mass monomial".into());
                    }
                    i64::try_from(exponent)?
                } else {
                    return Err("nonmonomial mass coefficient".into());
                };
                if exponent < 0 || exponent % 2 != 0 || degree - exponent / 2 != 10 {
                    return Err("physical D-j=10 homogeneity failed".into());
                }
                *histogram.entry((degree, exponent / 2)).or_default() += 1;
            }
            report.push_str(&format!(
                "key={:?} admitted={accepted}\n",
                term.integral().powers()
            ));
        }
        std::fs::write(evidence.join("prepared-preflight.txt"), &report)?;
    }
    report.push_str("# histogram stored_D explicit_mass_squared_power_j occurrences\n");
    for ((degree, mass), terms) in histogram {
        report.push_str(&format!("{degree} {mass} {terms}\n"));
    }
    report.push_str(&format!(
        "inputs={} lowered_occurrences={count} admitted_occurrences={admitted} guards={}\n",
        expression.0.len(),
        guards.len()
    ));
    std::fs::write(evidence.join("prepared-preflight.txt"), report)?;
    std::fs::write(
        evidence.join("prepared-preflight-guards.symbolica"),
        guards
            .iter()
            .map(AtomCore::to_canonical_string)
            .collect::<Vec<_>>()
            .join("\n"),
    )?;
    println!(
        "prepared_input_diagnostic inputs={} lowered_occurrences={count} admitted_occurrences={admitted} guards={}",
        expression.0.len(),
        guards.len()
    );
    Ok(())
}
