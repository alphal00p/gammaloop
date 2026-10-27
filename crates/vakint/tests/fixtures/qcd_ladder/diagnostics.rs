//! Manual, ignored failure diagnosis; never selected by the full slow class.
//! Acceptance still reruns the complete native QCD preparation in slow.rs.
use super::*;
use symbolica::{atom::AtomView, parser::ParseSettings, prelude::Rational};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CensusProjection {
    Quadratic,
    ZeroMomentum,
}

impl CensusProjection {
    fn from_requested(value: Option<&str>) -> Result<Self> {
        match value {
            None | Some("quadratic") => Ok(Self::Quadratic),
            Some("zero-momentum") => Ok(Self::ZeroMomentum),
            _ => Err("census projection must be quadratic or zero-momentum".into()),
        }
    }

    fn parameters(self) -> (&'static str, usize, i64) {
        match self {
            Self::Quadratic => ("quadratic-transverse", 21, 10),
            Self::ZeroMomentum => ("zero-momentum-metric", 1, 9),
        }
    }
}

#[test]
fn census_projection_configuration_preserves_physical_contracts() -> Result<()> {
    assert_eq!(
        CensusProjection::from_requested(None)?.parameters(),
        ("quadratic-transverse", 21, 10)
    );
    assert_eq!(
        CensusProjection::from_requested(Some("quadratic"))?,
        CensusProjection::Quadratic
    );
    assert_eq!(
        CensusProjection::from_requested(Some("zero-momentum"))?.parameters(),
        ("zero-momentum-metric", 1, 9)
    );
    for invalid in [
        "",
        "quadratic,zero-momentum",
        "quadratic zero-momentum",
        "zero-momentum-metric",
    ] {
        assert!(CensusProjection::from_requested(Some(invalid)).is_err());
    }
    Ok(())
}

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
                master_mode: crate::RustRedMasterMode::Retain,
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
    let scalar_limits = ScalarNumeratorLimits::default();
    println!("prepared_input_diagnostic scalar_limits={scalar_limits:?}");
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
                master_mode: crate::RustRedMasterMode::Retain,
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
        "# scalar_limits={scalar_limits:?}\n# Contribution ordinals are run-local; compare exact recorded targets/numerators.\n"
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
                    return Err("engineering mass-dimension balance D_key-j=10 failed".into());
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

/// One selected-witness handoff for the physical replay and frozen public corpus.
/// Native owner transport is inspected here; no saved IBP rule is applied.
struct CampaignCensus<'a> {
    engine: &'a Vakint,
    application: &'a SavedApplication,
    settings: &'a VakintSettings,
    binding: &'a serde_json::Value,
    evidence: &'a std::path::Path,
    physical: Option<CensusProjection>,
}

impl CampaignCensus<'_> {
    fn contribution(
        &self,
        index: usize,
        identity: &str,
        original: crate::VakintTerm,
    ) -> Result<serde_json::Value> {
        use serde_json::json;
        use std::time::Instant;
        let engine = self.engine;
        let application = self.application;
        let settings = self.settings;
        let binding = self.binding;
        let evidence = self.evidence;
        let mass = Atom::var(vk_symbol!("mUV"));
        let started = Instant::now();
        eprintln!("qcd_campaign_census contribution={index} phase=selected_witness_begin");
        let mut matched = engine
            .match_rustred_input(
                settings,
                original.integral.as_view(),
                original.numerator.as_view(),
            )?
            .ok_or("quadratic source has no verified geometric witness")?;
        let matching_seconds = started.elapsed().as_secs_f64();
        let selected_witness = matched.to_string();
        matched.apply_replacement_rules()?;
        let mut term = original;
        term.apply_scalar_numerator_replacement_rules(&matched, settings)?;
        let descriptor = select(&matched.canonical_topology)?;
        let input = MatchedScalarInput::prepare(
            term.numerator.as_view(),
            &matched,
            &descriptor.momenta,
            descriptor.physical_slots,
            15,
            settings,
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
            settings.rustred_scalar_numerator_limits,
        )?;
        capture_weighted_input(evidence, index, &input, prepared.family())?;
        let started = Instant::now();
        let lowered =
            scalar.lower_weighted(&input.numerator, &input.target, &input.mass_squared)?;
        let lowering_seconds = started.elapsed().as_secs_f64();
        let mut records = Vec::new();
        let mut batches = BTreeMap::<(Atom, i64), Vec<usize>>::new();
        let mut owner_images = Vec::new();
        let mut origin = None;
        let context = prepared.family().coefficient_context();
        let exact_limits = settings
            .rustred_application_options
            .reduction_limits
            .exact_algebra;
        for lowered_term in lowered.terms() {
            let key = lowered_term.integral();
            let weight = super::super::combine_exact(
                lowered_term.coefficient().to_expression() * lowered_term.scalar_spectator(),
            );
            let degree = key.powers().iter().sum::<i64>();
            let mut mass_powers = Vec::new();
            if let Some(physical) = self.physical {
                let (projection, _, engineering_degree_balance) = physical.parameters();
                for (monomial, coefficient) in weight.coefficient_list::<u16>(&[&mass]) {
                    if coefficient.is_zero() {
                        continue;
                    }
                    let exponent = if monomial.is_one() {
                        0
                    } else if monomial == mass {
                        1
                    } else if let AtomView::Pow(power) = monomial.as_view() {
                        if power.get_base() != mass.as_view() {
                            return Err("non-mass monomial".into());
                        }
                        i64::try_from(power.get_exp())?
                    } else {
                        return Err("nonmonomial mass dependence".into());
                    };
                    if exponent < 0
                        || exponent % 2 != 0
                        || degree - exponent / 2 != engineering_degree_balance
                    {
                        return Err(format!(
                        "mapped key violates exact {projection} engineering mass-dimension balance D_key-j={engineering_degree_balance}"
                    )
                    .into());
                    }
                    mass_powers.push(exponent / 2);
                }
            }
            let admitted = prepared.routed_entry_is_admitted(key)?;
            if !admitted {
                eprintln!(
                    "qcd_campaign_census uncovered contribution={index} family={} powers={:?} mass_squared_powers={mass_powers:?}",
                    prepared.family().name(),
                    key.powers()
                );
            }
            let mut native_images = Vec::new();
            let owner_route = if admitted || self.physical.is_none() {
                let routed = prepared.route_entry_to_owner(key)?;
                let mut images = Vec::new();
                for (endpoint, coefficient) in routed.terms() {
                    if coefficient.is_zero() {
                        continue;
                    }
                    let owner_weight =
                        super::super::combine_exact(&weight * coefficient.to_expression());
                    if self.physical.is_none() {
                        native_images.push((
                            endpoint.clone(),
                            context.try_mul(
                                lowered_term.coefficient(),
                                coefficient,
                                exact_limits,
                            )?,
                        ));
                    }
                    images.push(json!({"powers": endpoint.powers(), "weight": owner_weight.to_canonical_string(),
                        "transport_coefficient": coefficient.to_expression().to_canonical_string(),
                        "recorded_owner_admitted":prepared.entry_is_admitted(endpoint)?}));
                }
                Some(
                    json!({"images": images, "nonzero_conditions": routed.nonzero_conditions().iter()
                    .map(|coefficient| coefficient.to_expression().to_canonical_string()).collect::<Vec<_>>()}),
                )
            } else {
                None
            };
            let exponent = key
                .powers()
                .iter()
                .try_fold(0_i128, |sum, power| sum.checked_sub(i128::from(*power)))
                .ok_or("saved application mass exponent overflow")?;
            let exponent = i64::try_from(exponent)
                .map_err(|_| RustRedEvaluationError::MassExponentOverflow { exponent })?;
            if self.physical.is_none() {
                let anchor = *origin.get_or_insert(exponent);
                let channel = SavedApplication::mass_channel(
                    lowered_term.scalar_spectator(),
                    &input.mass_squared,
                    exponent,
                    anchor,
                    settings
                        .rustred_scalar_numerator_limits
                        .diagnostics_enabled
                        .then_some(
                            settings
                                .rustred_scalar_numerator_limits
                                .max_scalar_product_degree,
                        ),
                );
                batches.entry(channel).or_default().push(records.len());
            }
            owner_images.push(native_images);
            records.push(json!({"powers": key.powers(), "weight": weight.to_canonical_string(),
                "coefficient":lowered_term.coefficient().to_expression().to_canonical_string(),
                "scalar_spectator":lowered_term.scalar_spectator().to_canonical_string(),
                "original_mass_exponent":exponent,
                "D_key": degree,"r_idx":key.powers().iter().filter(|power|**power<0).map(|power|-power).sum::<i64>(),
                "omega_uv_at_d4_with_mass_fixed":20-2*degree,
                "admitted": admitted, "mass_squared_powers": mass_powers, "owner_route": owner_route}));
        }
        let mut requests = Vec::new();
        for ((spectator, exponent), indices) in batches {
            let ordinal = requests.len();
            let mut owner_terms = BTreeMap::<IntegralKey, Coefficient>::new();
            let mut transport_conditions = BTreeSet::new();
            for &row in &indices {
                for condition in records[row]["owner_route"]["nonzero_conditions"]
                    .as_array()
                    .ok_or("missing request transport conditions")?
                {
                    transport_conditions.insert(
                        condition
                            .as_str()
                            .ok_or("invalid transport condition")?
                            .to_owned(),
                    );
                }
                for (owner, coefficient) in &owner_images[row] {
                    if let std::collections::btree_map::Entry::Occupied(mut entry) =
                        owner_terms.entry(owner.clone())
                    {
                        let sum = context.try_add(entry.get(), coefficient, exact_limits)?;
                        if sum.is_zero() {
                            entry.remove();
                        } else {
                            *entry.get_mut() = sum;
                        }
                    } else if !coefficient.is_zero() {
                        owner_terms.insert(owner.clone(), coefficient.clone());
                    }
                }
            }
            let terms = indices.into_iter().map(|row| {
                records[row]["request_ordinal"] = json!(ordinal);
                json!({"source_root_row":row,"powers":lowered.terms()[row].integral().powers(),
                    "coefficient":lowered.terms()[row].coefficient().to_expression().to_canonical_string()})
            }).collect::<Vec<_>>();
            requests.push(json!({"ordinal":ordinal,"scalar_spectator":spectator.to_canonical_string(),
                "integer_mass_exponent":exponent,"terms":terms,
                "transport_nonzero_conditions":transport_conditions,
                "caller_nonzero_conditions":lowered.nonzero_conditions().iter().map(AtomCore::to_canonical_string).collect::<Vec<_>>(),
                "combined_owner_terms":owner_terms.into_iter().map(|(key,coefficient)|
                    json!({"powers":key.powers(),"coefficient":coefficient.to_expression().to_canonical_string(),"stage":"route"})).collect::<Vec<_>>(),
                "scope":"Exact production request before IBP; Native coefficients combined only within this batch"}));
        }
        let row = json!({"contribution": index, "source_identity":identity,"binding":binding,"topology": matched.canonical_topology.get_integral().name,
            "family_fingerprint": prepared.family().fingerprint(), "target": input.target.powers(),
            "selected_witness":selected_witness,"routed_numerator":input.numerator.to_canonical_string(),
            "request_batches":requests,
            "mass_squared": input.mass_squared.to_canonical_string(), "matching_seconds": matching_seconds,
            "lowering_seconds": lowering_seconds, "roots": records,
            "nonzero_conditions": lowered.nonzero_conditions().iter().map(AtomCore::to_canonical_string).collect::<Vec<_>>()});
        Ok(row)
    }
}

/// Authenticated inputs for the existing manual census, never reduction receipts.
#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct PublicCensusCorpus {
    schema: String,
    data_path: PathBuf,
    data_manifest_sha256: String,
    loop_manifest_sha256: String,
    application_sha256: String,
    source_inventory: CensusEvidence,
    cases: Vec<PublicCensusCase>,
}

#[derive(serde::Deserialize)]
#[serde(deny_unknown_fields)]
struct CensusEvidence {
    path: PathBuf,
    sha256: String,
}

impl CensusEvidence {
    fn bytes(&self) -> Result<Vec<u8>> {
        if !self.path.is_absolute() {
            return Err("census evidence paths must be absolute".into());
        }
        let bytes = std::fs::read(&self.path)?;
        if format!("{:x}", Sha256::digest(&bytes)) != self.sha256 {
            return Err(format!("census evidence changed: {}", self.path.display()).into());
        }
        Ok(bytes)
    }
}

#[derive(Clone, serde::Deserialize, serde::Serialize)]
#[serde(deny_unknown_fields)]
struct PublicCensusCase {
    id: String,
    role: String,
    public_input: String,
    public_input_sha256: String,
}

#[derive(serde::Deserialize)]
struct ExplicitSourceInventory {
    schema: String,
    raw_source_union_count: usize,
    direct_owner_source_union_count: usize,
    source_rows: Vec<ExplicitTrainingSource>,
}

#[derive(serde::Deserialize)]
struct ExplicitTrainingSource {
    source: String,
    source_raw_keys: Vec<Vec<i64>>,
    source_direct_owner_keys: Vec<Vec<i64>>,
}

#[derive(serde::Serialize)]
struct SourceQualification {
    source: String,
    raw_keys: BTreeSet<Vec<i64>>,
    direct_owner_keys: BTreeSet<Vec<i64>>,
    historical_direct_owner_keys: BTreeSet<Vec<i64>>,
    current_owner_images_match_historical: bool,
}

#[derive(serde::Serialize)]
struct OwnerTransportImage {
    powers: Vec<i64>,
    coefficient: String,
    recorded_owner_admitted: bool,
}

#[derive(serde::Serialize)]
struct SourceTransport {
    powers: Vec<i64>,
    recorded_coverage_admitted: bool,
    images: Vec<OwnerTransportImage>,
    nonzero_conditions: Vec<String>,
}

#[derive(serde::Serialize)]
struct CensusRoots {
    #[serde(rename = "actual_lowered_roots")]
    raw: BTreeSet<Vec<i64>>,
    #[serde(rename = "protected_owner_roots")]
    owners: BTreeSet<Vec<i64>>,
    #[serde(rename = "original_coverage_admitted")]
    covered: bool,
}

#[derive(serde::Serialize)]
struct CaseSourceOverlap {
    source: String,
    raw_overlap: Vec<Vec<i64>>,
    owner_overlap: Vec<Vec<i64>>,
}

#[derive(serde::Serialize)]
struct CasePairOverlap {
    left: String,
    right: String,
    raw_overlap: Vec<Vec<i64>>,
    owner_overlap: Vec<Vec<i64>>,
}

#[derive(serde::Serialize)]
struct CensusLowering {
    parsed_input: String,
    canonical_scalar_handoff: String,
    #[serde(flatten)]
    roots: CensusRoots,
    explicit_sources_disjoint: bool,
    nonempty_candidate: bool,
    source_overlaps: Vec<CaseSourceOverlap>,
    contributions: Vec<serde_json::Value>,
}

#[derive(serde::Serialize)]
struct CaseQualification {
    #[serde(flatten)]
    case: PublicCensusCase,
    family_fingerprint: String,
    application_binding: serde_json::Value,
    lowering_complete: bool,
    #[serde(flatten)]
    lowering: Option<CensusLowering>,
    error: Option<String>,
    native_before: String,
    native_after: String,
    native_work_unchanged: bool,
    reduction_succeeded: Option<bool>,
    ibp_application: bool,
    frozen: bool,
}

impl CaseQualification {
    fn qualified(&self) -> bool {
        self.lowering_complete
            && self.lowering.as_ref().is_some_and(|lowering| {
                lowering.roots.covered
                    && lowering.explicit_sources_disjoint
                    && lowering.nonempty_candidate
            })
    }
}

impl PublicCensusCorpus {
    fn validate(&self) -> Result<()> {
        if self.schema != "vakint.public-input-lowering-corpus.v1" || self.cases.is_empty() {
            return Err("unknown or empty public lowering corpus".into());
        }
        if self
            .cases
            .iter()
            .filter(|case| case.role == "held_out")
            .count()
            < 12
            || !self.cases.iter().any(|case| case.role == "tuning")
        {
            return Err("freeze at least twelve held-out inputs and one tuning input".into());
        }
        if !self.data_path.is_absolute() {
            return Err("public census requires an absolute local data checkout".into());
        }
        for hash in [
            &self.data_manifest_sha256,
            &self.loop_manifest_sha256,
            &self.application_sha256,
            &self.source_inventory.sha256,
        ] {
            if hash.len() != 64 || !hash.bytes().all(|byte| byte.is_ascii_hexdigit()) {
                return Err("invalid public census evidence digest".into());
            }
        }
        let mut ids = BTreeSet::new();
        let mut inputs = BTreeSet::new();
        for case in &self.cases {
            if case.id.is_empty()
                || !case
                    .id
                    .bytes()
                    .all(|c| c.is_ascii_alphanumeric() || b"_-".contains(&c))
                || !ids.insert(&case.id)
                || !matches!(case.role.as_str(), "tuning" | "held_out")
                || case.public_input.trim().is_empty()
                || !inputs.insert(&case.public_input)
                || format!("{:x}", Sha256::digest(case.public_input.as_bytes()))
                    != case.public_input_sha256
            {
                return Err("invalid, duplicate or changed frozen public input".into());
            }
        }
        Ok(())
    }

    /// A candidate inspection retains the full declaration and its role, but
    /// never certifies that the complete held-out suite has been qualified.
    fn selected_cases(&self, candidate: Option<&str>) -> Result<Vec<&PublicCensusCase>> {
        match candidate {
            None => Ok(self.cases.iter().collect()),
            Some(id) => self
                .cases
                .iter()
                .find(|case| case.id == id)
                .map(|case| vec![case])
                .ok_or_else(|| {
                    format!("candidate {id:?} is absent from the input declaration").into()
                }),
        }
    }

    fn validate_data(&self) -> Result<()> {
        for (file, expected) in [
            ("manifest.json", &self.data_manifest_sha256),
            ("5L/manifest.json", &self.loop_manifest_sha256),
            ("5L/ibps/application.json", &self.application_sha256),
        ] {
            if format!(
                "{:x}",
                Sha256::digest(std::fs::read(self.data_path.join(file))?)
            ) != *expected
            {
                return Err(format!("public census data binding changed: {file}").into());
            }
        }
        Ok(())
    }

    /// Keep request identities even when display weights cancel across batches.
    fn root_sets(contributions: &[serde_json::Value]) -> Result<CensusRoots> {
        let mut raw = BTreeSet::new();
        let mut owners = BTreeSet::new();
        let mut covered = true;
        for contribution in contributions {
            for root in contribution["roots"]
                .as_array()
                .ok_or("missing complete raw roots")?
            {
                let powers: Vec<i64> = serde_json::from_value(root["powers"].clone())?;
                if powers.len() != 15 {
                    return Err("invalid raw root arity".into());
                }
                raw.insert(powers);
                covered &= root["admitted"] == true;
                for image in root["owner_route"]["images"]
                    .as_array()
                    .ok_or("missing complete owner transport")?
                {
                    let powers: Vec<i64> = serde_json::from_value(image["powers"].clone())?;
                    if powers.len() != 15 {
                        return Err("invalid owner root arity".into());
                    }
                    owners.insert(powers);
                }
            }
        }
        Ok(CensusRoots {
            raw,
            owners,
            covered,
        })
    }

    fn write_document(
        directory: &std::path::Path,
        name: &str,
        value: &serde_json::Value,
    ) -> Result<()> {
        use std::io::Write;
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(directory.join(name))?;
        serde_json::to_writer_pretty(&mut file, value)?;
        file.write_all(b"\n")?;
        file.sync_all()?;
        std::fs::File::open(directory)?.sync_all()?;
        Ok(())
    }

    fn keys(keys: &[Vec<i64>]) -> Result<BTreeSet<Vec<i64>>> {
        let set = keys.iter().cloned().collect::<BTreeSet<_>>();
        if keys.iter().any(|key| key.len() != 15) || keys.len() != set.len() {
            return Err("census keys must be unique ordered fifteen-axis identities".into());
        }
        Ok(set)
    }

    fn compile_identity() -> Result<serde_json::Value> {
        let mut identity = serde_json::Map::new();
        for (name, value) in [
            (
                "base_main_revision",
                option_env!("VAKINT_LOWERING_BUILD_BASE_REVISION"),
            ),
            (
                "private_source_snapshot_sha256",
                option_env!("VAKINT_LOWERING_BUILD_SOURCE_SHA256"),
            ),
            (
                "native_revision",
                option_env!("VAKINT_LOWERING_BUILD_NATIVE_REVISION"),
            ),
            (
                "native_source_snapshot_sha256",
                option_env!("VAKINT_LOWERING_BUILD_NATIVE_SHA256"),
            ),
            (
                "kernel_revision",
                option_env!("VAKINT_LOWERING_BUILD_KERNEL_REVISION"),
            ),
            ("rustc", option_env!("VAKINT_LOWERING_BUILD_RUSTC")),
            ("profile", option_env!("VAKINT_LOWERING_BUILD_PROFILE")),
        ] {
            let value = value
                .filter(|value| !value.trim().is_empty())
                .ok_or_else(|| format!("missing census build identity: {name}"))?;
            let expected = if name.ends_with("sha256") {
                Some(64)
            } else if name.ends_with("revision") {
                Some(40)
            } else {
                None
            };
            if expected.is_some_and(|length| {
                value.len() != length || !value.bytes().all(|byte| byte.is_ascii_hexdigit())
            }) {
                return Err(format!("invalid census build identity: {name}").into());
            }
            identity.insert(name.to_owned(), value.into());
        }
        identity.insert(
            "diagnostic_source_sha256".into(),
            format!("{:x}", Sha256::digest(include_bytes!("diagnostics.rs"))).into(),
        );
        identity.insert(
            "executable_sha256".into(),
            format!(
                "{:x}",
                Sha256::digest(std::fs::read(std::env::current_exe()?)?)
            )
            .into(),
        );
        Ok(identity.into())
    }

    fn run(path: PathBuf) -> Result<()> {
        use crate::{
            EvaluationMethod, EvaluationOrder, OfficialExternalDataOptions,
            RustRedEvaluationOptions, RustRedMasterMode,
        };
        use serde_json::json;

        let bytes = std::fs::read(&path)?;
        let corpus: Self = serde_json::from_slice(&bytes)?;
        corpus.validate()?;
        let candidate = match std::env::var("VAKINT_PUBLIC_LOWERING_CASE") {
            Ok(id) => Some(id),
            Err(std::env::VarError::NotPresent) => None,
            Err(error) => return Err(error.into()),
        };
        let selected_cases = corpus.selected_cases(candidate.as_deref())?;
        let source_identity = Self::compile_identity()?;
        let inventory_bytes = corpus.source_inventory.bytes()?;
        let inventory: ExplicitSourceInventory = serde_json::from_slice(&inventory_bytes)?;
        if inventory.schema != "vakint.five-source-current-public-root-audit.v1" {
            return Err("unknown five-source exclusion inventory".into());
        }
        let sources = &inventory.source_rows;
        let expected = BTreeSet::from([
            "accepted-small-two-root",
            "old-completed-node-prefix",
            "q0-current-public-replay",
            "tuning-120",
            "q-monitor-1089",
        ]);
        let actual = sources
            .iter()
            .map(|source| source.source.as_str())
            .collect::<BTreeSet<_>>();
        if actual != expected || sources.len() != expected.len() {
            return Err(
                "all five explicit training source categories are required exactly once".into(),
            );
        }
        corpus.validate_data()?;
        let evidence = PathBuf::from(
            std::env::var_os("VAKINT_QCD_WEIGHT_EVIDENCE")
                .ok_or("set new census evidence directory")?,
        );
        if !evidence.is_absolute() {
            return Err("census evidence directory must be absolute".into());
        }
        std::fs::create_dir(&evidence)?;
        let mut settings = VakintSettings {
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions {
                    master_mode: RustRedMasterMode::Retain,
                },
            )]),
            form_exe_path: "/public-census-must-not-use-form".into(),
            use_dot_product_notation: true,
            ..Default::default()
        };
        settings.rustred_application_options.unresolved_policy =
            rustred::solver::UnresolvedPolicy::Preserve;
        settings.rustred_application_options.entry_coverage_policy =
            rustred_app::SavedApplicationEntryCoveragePolicy::Advisory;
        settings.rustred_application_options.cache_representation =
            rustred::solver::CandidateCacheRepresentation::Sparse;
        settings.rustred_application_options.preparation_workers = 1;
        let mut engine = Vakint::new()?;
        let load = engine.load_official_external_data(
            5,
            &OfficialExternalDataOptions {
                local_artifacts_path: Some(corpus.data_path.clone()),
                load_prepared_caches: false,
                ..Default::default()
            },
        )?;
        let cache_attachments = (1..=5)
            .map(|loops| -> Result<_> {
                Ok((loops, engine.cache_statistics(loops)?.has_read_backing()))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;
        if load
            .loops
            .iter()
            .map(|pack| pack.loop_count)
            .collect::<Vec<_>>()
            != [1, 2, 3, 4, 5]
            || load
                .loops
                .iter()
                .any(|pack| pack.cache_loaded || pack.source != "local")
            || cache_attachments.values().any(|attached| *attached)
        {
            return Err(
                "census must load only local cumulative data without prepared caches".into(),
            );
        }
        engine.validate_settings(&settings)?;
        let application = engine.saved_application(&settings)?;
        let binding = engine
            .cache_statistics(5)?
            .saved_application
            .ok_or("missing application binding")?
            .to_document()["binding"]
            .clone();
        let native_before = {
            let prepared = application
                .prepared
                .lock()
                .map_err(|error| error.to_string())?;
            let owners = prepared.owner_preparation_statistics();
            if owners.preparation_attempts != 0
                || owners.prepared_owners != 0
                || prepared.cache_statistics().rule_applications() != 0
            {
                return Err(
                    "census session unexpectedly prepared owner formulas or applied rules".into(),
                );
            }
            (
                owners,
                prepared.cache_statistics(),
                prepared.cache_work_statistics(),
            )
        };
        let mut mapped_sources = Vec::new();
        let mut training_raw = BTreeSet::new();
        let mut training_owners = BTreeSet::new();
        let mut historical_training_owners = BTreeSet::new();
        let mut training_transports = BTreeMap::<Vec<i64>, SourceTransport>::new();
        let mut sources_agree = true;
        for source in sources {
            let raw = Self::keys(&source.source_raw_keys)?;
            let historical_owners = Self::keys(&source.source_direct_owner_keys)?;
            let mut owners = BTreeSet::new();
            for key in &raw {
                if !training_transports.contains_key(key) {
                    let prepared = application
                        .prepared
                        .lock()
                        .map_err(|error| error.to_string())?;
                    let integral = IntegralKey::try_new(key.clone())?;
                    let mapped = prepared.route_entry_to_owner(&integral)?;
                    let images = mapped
                        .terms()
                        .iter()
                        .filter(|(_, coefficient)| !coefficient.is_zero())
                        .map(|(owner, coefficient)| -> Result<OwnerTransportImage> {
                            Ok(OwnerTransportImage {
                                powers: owner.powers().to_vec(),
                                coefficient: coefficient.to_expression().to_canonical_string(),
                                recorded_owner_admitted: prepared.entry_is_admitted(owner)?,
                            })
                        })
                        .collect::<Result<Vec<_>>>()?;
                    // This is routed_entry_is_admitted's original endpoint
                    // predicate after a successful verified route, without
                    // executing that same transport a second time.
                    training_transports.insert(
                        key.clone(),
                        SourceTransport {
                            powers: key.clone(),
                            recorded_coverage_admitted: images
                                .iter()
                                .all(|image| image.recorded_owner_admitted),
                            images,
                            nonzero_conditions: mapped
                                .nonzero_conditions()
                                .iter()
                                .map(|guard| guard.to_expression().to_canonical_string())
                                .collect(),
                        },
                    );
                }
                for image in &training_transports[key].images {
                    owners.insert(image.powers.clone());
                }
            }
            sources_agree &= owners == historical_owners;
            training_raw.extend(raw.iter().cloned());
            training_owners.extend(owners.iter().cloned());
            historical_training_owners.extend(historical_owners.iter().cloned());
            mapped_sources.push(SourceQualification {
                source: source.source.clone(),
                raw_keys: raw,
                current_owner_images_match_historical: owners == historical_owners,
                direct_owner_keys: owners,
                historical_direct_owner_keys: historical_owners,
            });
        }
        if training_raw.len() != inventory.raw_source_union_count
            || historical_training_owners.len() != inventory.direct_owner_source_union_count
        {
            return Err(
                "complete source inventories differ from authenticated union counts".into(),
            );
        }
        let family = application
            .prepared
            .lock()
            .map_err(|error| error.to_string())?
            .family()
            .fingerprint()
            .to_owned();
        let mut rows = Vec::new();
        for case in selected_cases {
            engine.clear_rustred_application_cache()?;
            let case_native_before = {
                let prepared = application
                    .prepared
                    .lock()
                    .map_err(|error| error.to_string())?;
                (
                    prepared.owner_preparation_statistics(),
                    prepared.cache_statistics(),
                    prepared.cache_work_statistics(),
                )
            };
            let case_dir = evidence.join(&case.id);
            std::fs::create_dir(&case_dir)?;
            let outcome = (|| -> Result<CensusLowering> {
                let input =
                    crate::vakint_parse!(&case.public_input).map_err(|error| error.to_string())?;
                let mut expression = VakintExpression::try_from(input.as_view())?;
                expression.canonicalize(&settings, &engine.topologies, false)?;
                expression.tensor_reduce(&engine, &settings)?;
                let handoff = Atom::from(expression);
                let terms = VakintExpression::split_integrals(handoff.as_view())?;
                let census = CampaignCensus {
                    engine: &engine,
                    application: &application,
                    settings: &settings,
                    binding: &binding,
                    evidence: &case_dir,
                    physical: None,
                };
                let mut contributions = Vec::new();
                for (index, term) in terms.into_iter().enumerate() {
                    let identity = format!(
                        "{:x}",
                        Sha256::digest(serde_json::to_vec(&(
                            term.integral.to_canonical_string(),
                            term.numerator.to_canonical_string()
                        ))?)
                    );
                    let row = census.contribution(index, &identity, term)?;
                    contributions.push(row);
                }
                let roots = Self::root_sets(&contributions)?;
                let overlaps = mapped_sources
                    .iter()
                    .map(|source| CaseSourceOverlap {
                        source: source.source.clone(),
                        raw_overlap: roots.raw.intersection(&source.raw_keys).cloned().collect(),
                        owner_overlap: roots
                            .owners
                            .intersection(&source.direct_owner_keys)
                            .cloned()
                            .collect(),
                    })
                    .collect();
                Ok(CensusLowering {
                    parsed_input: input.to_canonical_string(),
                    canonical_scalar_handoff: handoff.to_canonical_string(),
                    explicit_sources_disjoint: roots.raw.is_disjoint(&training_raw)
                        && roots.owners.is_disjoint(&training_owners),
                    nonempty_candidate: !roots.raw.is_empty() && !roots.owners.is_empty(),
                    roots,
                    source_overlaps: overlaps,
                    contributions,
                })
            })();
            let (lowering, error) = match outcome {
                Ok(lowering) => (Some(lowering), None),
                Err(error) => (None, Some(error.to_string())),
            };
            let case_native_after = {
                let prepared = application
                    .prepared
                    .lock()
                    .map_err(|error| error.to_string())?;
                (
                    prepared.owner_preparation_statistics(),
                    prepared.cache_statistics(),
                    prepared.cache_work_statistics(),
                )
            };
            let unchanged = case_native_before == case_native_after;
            let row = CaseQualification {
                case: case.clone(),
                family_fingerprint: family.clone(),
                application_binding: binding.clone(),
                lowering_complete: lowering.is_some() && unchanged,
                lowering,
                error: if unchanged {
                    error
                } else {
                    Some("lower-only case changed Native owner or cache work".into())
                },
                native_before: format!("{case_native_before:?}"),
                native_after: format!("{case_native_after:?}"),
                native_work_unchanged: unchanged,
                reduction_succeeded: None,
                ibp_application: false,
                frozen: false,
            };
            Self::write_document(
                &case_dir,
                "qualification.json",
                &serde_json::to_value(&row)?,
            )?;
            rows.push(row);
        }
        let native_after = {
            let prepared = application
                .prepared
                .lock()
                .map_err(|error| error.to_string())?;
            (
                prepared.owner_preparation_statistics(),
                prepared.cache_statistics(),
                prepared.cache_work_statistics(),
            )
        };
        if native_before != native_after
            || application
                .cache
                .lock()
                .map_err(|error| error.to_string())?
                .is_some()
        {
            return Err(
                "lower-only corpus changed native owner/cache work or attached capture".into(),
            );
        }
        let mut pairs = Vec::new();
        for (index, row) in rows.iter().enumerate() {
            let Some(lowering) = &row.lowering else {
                continue;
            };
            for other in rows.iter().skip(index + 1) {
                let Some(other_lowering) = &other.lowering else {
                    continue;
                };
                pairs.push(CasePairOverlap {
                    left: row.case.id.clone(),
                    right: other.case.id.clone(),
                    raw_overlap: lowering
                        .roots
                        .raw
                        .intersection(&other_lowering.roots.raw)
                        .cloned()
                        .collect(),
                    owner_overlap: lowering
                        .roots
                        .owners
                        .intersection(&other_lowering.roots.owners)
                        .cloned()
                        .collect(),
                });
            }
        }
        let complete = sources_agree
            && rows.iter().all(CaseQualification::qualified)
            && pairs
                .iter()
                .all(|pair| pair.raw_overlap.is_empty() && pair.owner_overlap.is_empty());
        let report = json!({"schema":if candidate.is_some() {"vakint.public-input-lowering-candidate-qualification.v1"} else {"vakint.public-input-lowering-owner-qualification.v1"},
            "qualification_complete":candidate.is_none()&&complete,"candidate_qualified":candidate.as_ref().map(|_|complete),
            "requested_candidate":candidate,"declared_case_count":corpus.cases.len(),
            "source_identity":source_identity,"input_corpus":{"path":path,"sha256":format!("{:x}",Sha256::digest(&bytes))},
            "source_inventory":{"path":corpus.source_inventory.path,"sha256":corpus.source_inventory.sha256},"data":load,
            "persistent_cache_attachments":cache_attachments,
            "data_manifest_sha256":corpus.data_manifest_sha256,"loop_manifest_sha256":corpus.loop_manifest_sha256,"application_sha256":corpus.application_sha256,
            "semantic_binding":binding,"family_fingerprint":family,"application_options":format!("{:?}",settings.rustred_application_options),
            "scalar_limits":format!("{:?}",settings.rustred_scalar_numerator_limits),"matching_limits":format!("{:?}",settings.rustred_momentum_match_limits),
            "native_before":format!("{native_before:?}"),"native_after":format!("{native_after:?}"),"native_work_unchanged":true,
            "coverage_authority":"Original recorded entry predicates; Advisory permits verified owner inspection only",
            "sources":mapped_sources,"source_transports":training_transports.into_values().collect::<Vec<_>>(),"source_owner_images_match_historical":sources_agree,
            "raw_source_union_count":training_raw.len(),"direct_owner_source_union_count":training_owners.len(),"cases":rows,"case_overlaps":pairs,
            "reduction_succeeded":null,"ibp_application":false,"measurement_claim_allowed":false,"cache_selection_claim":false,"frozen":false,
            "acceptance":"Requires independent coefficient/guard/exclusion postflight before freeze; incidental cache descendants are not excluded"});
        if std::fs::read(&path)? != bytes || corpus.source_inventory.bytes()? != inventory_bytes {
            return Err("corpus or source inventory changed during inspection".into());
        }
        corpus.validate_data()?;
        Self::write_document(&evidence, "public-corpus-qualification.json", &report)?;
        if !complete {
            return Err("public corpus has incomplete, uncovered, zero or overlapping candidates; see qualification report".into());
        }
        Ok(())
    }
}

#[test]
fn public_census_requires_complete_immutable_input_identities() -> Result<()> {
    use serde_json::{Value, json};
    let cases = (0..13).map(|index| {
        let input = format!("topo(prop(1,edge(1,1),k(1),1,{}))", index + 1);
        json!({"id":format!("case-{index}"),"role":if index==0 {"tuning"} else {"held_out"},
            "public_input_sha256":format!("{:x}",Sha256::digest(input.as_bytes())),"public_input":input})
    }).collect::<Vec<_>>();
    let document = json!({"schema":"vakint.public-input-lowering-corpus.v1","data_path":"/unused-local-checkout",
        "data_manifest_sha256":"0".repeat(64),"loop_manifest_sha256":"1".repeat(64),"application_sha256":"2".repeat(64),
        "source_inventory":{"path":"/unused-source-inventory","sha256":"3".repeat(64)},"cases":cases});
    let valid = |value: Value| -> bool {
        serde_json::from_value::<PublicCensusCorpus>(value)
            .is_ok_and(|corpus| corpus.validate().is_ok())
    };
    assert!(valid(document.clone()));
    let corpus: PublicCensusCorpus = serde_json::from_value(document.clone())?;
    assert_eq!(corpus.selected_cases(None)?.len(), 13);
    let selected = corpus.selected_cases(Some("case-12"))?;
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].id, "case-12");
    assert_eq!(selected[0].role, "held_out");
    assert!(corpus.selected_cases(Some("missing-case")).is_err());
    let mut changed = document.clone();
    changed["cases"].as_array_mut().unwrap().pop();
    assert!(!valid(changed), "eleven held-out inputs are insufficient");
    for (field, replacement) in [
        ("id", json!("case-1")),
        ("id", json!("../escape")),
        ("role", json!("training")),
        ("role", json!("held_out")),
        ("public_input", json!("changed")),
        ("public_input_sha256", json!("0".repeat(64))),
    ] {
        let mut changed = document.clone();
        changed["cases"][0][field] = replacement;
        assert!(!valid(changed), "invalid corpus field accepted: {field}");
    }
    let mut changed = document.clone();
    changed["reduction_succeeded"] = json!(true);
    assert!(
        !valid(changed),
        "lower-only inputs cannot masquerade as reduction receipts"
    );
    let mut changed = document;
    changed["cases"][0]["unknown"] = json!(true);
    assert!(!valid(changed));
    Ok(())
}

#[test]
fn public_census_conservative_roots_survive_cross_request_cancellation() -> Result<()> {
    use serde_json::json;
    let raw = vec![1_i64; 15];
    let owner = vec![2_i64; 15];
    let zero = vec![0_i64; 15];
    let contributions = [
        json!({"roots":[{"powers":raw,"weight":"1","request_ordinal":0,"admitted":true,
            "owner_route":{"images":[{"powers":owner,"weight":"1"}],"nonzero_conditions":["d-2"]}}]}),
        json!({"roots":[{"powers":raw,"weight":"-1","request_ordinal":1,"admitted":true,
            "owner_route":{"images":[{"powers":owner,"weight":"-1"}],"nonzero_conditions":["d-2"]}},
            {"powers":zero,"weight":"1","request_ordinal":2,"admitted":false,
            "owner_route":{"images":[],"nonzero_conditions":["d"]}}]}),
    ];
    let roots = PublicCensusCorpus::root_sets(&contributions)?;
    assert_eq!(roots.raw, BTreeSet::from([raw, zero]));
    assert_eq!(roots.owners, BTreeSet::from([owner]));
    assert!(!roots.covered);
    assert_eq!(
        contributions[1]["roots"][1]["owner_route"]["nonzero_conditions"],
        json!(["d"])
    );
    Ok(())
}

#[test]
fn public_census_checks_key_arity_and_preserves_public_symbol_namespaces() -> Result<()> {
    assert!(PublicCensusCorpus::keys(&[vec![1_i64; 14]]).is_err());
    assert!(PublicCensusCorpus::keys(&[vec![1_i64; 15], vec![1_i64; 15]]).is_err());
    assert_eq!(PublicCensusCorpus::keys(&[vec![1_i64; 15]])?.len(), 1);
    Vakint::initialize_vakint_symbols();
    let input = crate::vakint_parse!("dot(k(1),k(1))*topo(prop(1,edge(1,1),k(1),1,1))")
        .map_err(|error| error.to_string())?;
    let symbols = input.get_all_symbols(true);
    for expected in [S.dot, S.k, S.topo, S.prop, S.edge] {
        assert!(
            symbols.contains(&expected),
            "public input namespace changed"
        );
    }
    Ok(())
}

/// Replays an explicitly supplied prepared integrand for an early campaign
/// census. Fresh native preparation is still required by the slow acceptance.
#[test]
#[ignore = "manual selected-witness census; requires prepared input and explicit local external data"]
fn prepared_ladder_selected_campaign_census_without_ibps() -> Result<()> {
    if let Some(path) = std::env::var_os("VAKINT_PUBLIC_LOWERING_CORPUS") {
        return PublicCensusCorpus::run(path.into());
    }
    use crate::{
        EvaluationMethod, EvaluationOrder, OfficialExternalDataOptions, RustRedEvaluationOptions,
        RustRedMasterMode,
    };
    use serde_json::json;
    use std::time::Instant;

    let source_path = PathBuf::from(
        std::env::var_os("VAKINT_QCD_PREPARED_INPUT").ok_or("set VAKINT_QCD_PREPARED_INPUT")?,
    );
    let requested_projection = match std::env::var("VAKINT_QCD_CENSUS_PROJECTION") {
        Ok(value) => Some(value),
        Err(std::env::VarError::NotPresent) => None,
        Err(error) => return Err(error.into()),
    };
    let (projection, expected_groups, engineering_degree_balance) =
        CensusProjection::from_requested(requested_projection.as_deref())?.parameters();
    let data_path = PathBuf::from(
        std::env::var_os("VAKINT_EXTERNAL_DATA").ok_or("set VAKINT_EXTERNAL_DATA to a checkout")?,
    );
    let evidence = PathBuf::from(
        std::env::var_os("VAKINT_QCD_WEIGHT_EVIDENCE").ok_or("set VAKINT_QCD_WEIGHT_EVIDENCE")?,
    );
    std::fs::create_dir_all(&evidence)?;
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/campaign-census-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    };
    settings.rustred_application_options.unresolved_policy =
        rustred::solver::UnresolvedPolicy::Preserve;
    eprintln!("qcd_campaign_census phase=engine_construction_begin");
    let mut engine = Vakint::new()?;
    eprintln!("qcd_campaign_census phase=data_loading_begin");
    let started = Instant::now();
    let load = engine.load_official_external_data(
        5,
        &OfficialExternalDataOptions {
            local_artifacts_path: Some(data_path),
            ..Default::default()
        },
    )?;
    let application = engine.saved_application(&settings)?;
    let load_seconds = started.elapsed().as_secs_f64();
    eprintln!("qcd_campaign_census phase=data_loading_end seconds={load_seconds:.6}");
    let source = std::fs::read(&source_path)?;
    let source_sha256 = format!("{:x}", Sha256::digest(&source));
    let parsing_started = Instant::now();
    eprintln!(
        "qcd_campaign_census phase=input_parsing_begin bytes={}",
        source.len()
    );
    let input = Atom::parse(
        std::str::from_utf8(&source)?,
        "qcd_ladder",
        ParseSettings::default(),
    )?;
    eprintln!(
        "qcd_campaign_census phase=input_parsing_end seconds={:.6}",
        parsing_started.elapsed().as_secs_f64()
    );
    let started = Instant::now();
    eprintln!("qcd_campaign_census phase=canonicalization_begin");
    let mut expression = VakintExpression::try_from(input.as_view())?;
    expression.canonicalize(&settings, &engine.topologies, false)?;
    eprintln!(
        "qcd_campaign_census phase=canonicalization_end seconds={:.6}",
        started.elapsed().as_secs_f64()
    );
    let tensor_started = Instant::now();
    eprintln!("qcd_campaign_census phase=tensor_prepass_begin");
    expression.tensor_reduce(&engine, &settings)?;
    eprintln!(
        "qcd_campaign_census phase=tensor_prepass_end seconds={:.6}",
        tensor_started.elapsed().as_secs_f64()
    );
    let terms = VakintExpression::split_integrals(Atom::from(expression).as_view())?;
    let canonical_seconds = started.elapsed().as_secs_f64();
    assert_eq!(
        terms.len(),
        expected_groups,
        "all physical projection groups must survive"
    );
    let binding = json!({"schema":"vakint.qcd-ladder.census-checkpoint.v1",
        "input_sha256":source_sha256,
        "data_manifest_sha256":format!("{:x}",Sha256::digest(std::fs::read(load.path.join("manifest.json"))?)),
        "application_manifest_sha256":format!("{:x}",Sha256::digest(std::fs::read(load.path.join("5L/ibps/application.json"))?)),
        "scalar_limits":format!("{:?}",settings.rustred_scalar_numerator_limits),
        "application_options":format!("{:?}",settings.rustred_application_options)});
    // Fresh processes can enumerate the same sources in a different order.
    // Ordinals are presentation metadata, never checkpoint identities.
    let mut checkpoints = BTreeMap::new();
    for entry in std::fs::read_dir(&evidence)? {
        let path = entry?.path();
        if !path
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.starts_with("contribution-") && name.ends_with(".json"))
        {
            continue;
        }
        let bytes = std::fs::read(&path)?;
        let row: serde_json::Value = serde_json::from_slice(&bytes)?;
        let identity = row["source_identity"]
            .as_str()
            .ok_or("checkpoint has no exact source identity")?
            .to_owned();
        if checkpoints.insert(identity, (path, bytes, row)).is_some() {
            return Err("duplicate exact-source census checkpoints".into());
        }
    }
    let mut contributions = Vec::new();
    for (index, original) in terms.into_iter().enumerate() {
        let identity = format!(
            "{:x}",
            Sha256::digest(serde_json::to_vec(&(
                original.integral.to_canonical_string(),
                original.numerator.to_canonical_string()
            ))?)
        );
        eprintln!("qcd_campaign_census contribution={index} source_identity={identity}");
        if let Some((checkpoint, bytes, mut row)) = checkpoints.remove(&identity) {
            let started = Instant::now();
            for field in [
                "schema",
                "input_sha256",
                "data_manifest_sha256",
                "application_manifest_sha256",
            ] {
                if row["binding"][field] != binding[field] {
                    return Err(format!(
                        "census checkpoint {} has another semantic binding ({field})",
                        checkpoint.display()
                    )
                    .into());
                }
            }
            let prepared = application
                .prepared
                .lock()
                .map_err(|error| error.to_string())?;
            if row["family_fingerprint"] != json!(prepared.family().fingerprint()) {
                return Err("checkpoint family identity changed".into());
            }
            for root in row["roots"].as_array().ok_or("checkpoint has no roots")? {
                let key = IntegralKey::try_new(serde_json::from_value::<Vec<i64>>(
                    root["powers"].clone(),
                )?)?;
                if root["admitted"] != json!(prepared.routed_entry_is_admitted(&key)?) {
                    return Err("checkpoint admission disagrees with active predicates".into());
                }
            }
            row["imported_checkpoint"] = json!({"path": checkpoint, "sha256":format!("{:x}", Sha256::digest(&bytes)),
                "original_contribution": row["contribution"], "completed_under_prior_policy": row["binding"] != binding});
            row["contribution"] = json!(index);
            eprintln!(
                "qcd_campaign_census contribution={index} phase=checkpoint_reused seconds={:.6}",
                started.elapsed().as_secs_f64()
            );
            contributions.push(row);
            continue;
        }
        let row = CampaignCensus {
            engine: &engine,
            application: &application,
            settings: &settings,
            binding: &binding,
            evidence: &evidence,
            physical: Some(CensusProjection::from_requested(
                requested_projection.as_deref(),
            )?),
        }
        .contribution(index, &identity, original)?;
        std::fs::write(
            evidence.join(format!("contribution-{identity}.json")),
            serde_json::to_vec_pretty(&row)?,
        )?;
        contributions.push(row);
    }
    let mut roots = BTreeMap::<IntegralKey, (Atom, usize, BTreeSet<usize>)>::new();
    let mut owner_roots = BTreeMap::<IntegralKey, (Atom, usize, BTreeSet<usize>)>::new();
    for (index, contribution) in contributions.iter().enumerate() {
        for root in contribution["roots"]
            .as_array()
            .ok_or("checkpoint has no roots")?
        {
            let key =
                IntegralKey::try_new(serde_json::from_value::<Vec<i64>>(root["powers"].clone())?)?;
            let weight = Atom::parse(
                root["weight"].as_str().ok_or("checkpoint has no weight")?,
                "qcd_ladder",
                ParseSettings::default(),
            )?;
            let row = roots
                .entry(key)
                .or_insert_with(|| (Atom::Zero, 0, BTreeSet::new()));
            row.0 += weight;
            row.1 += 1;
            row.2.insert(index);
            if let Some(images) = root["owner_route"]["images"].as_array() {
                for image in images {
                    let key = IntegralKey::try_new(serde_json::from_value::<Vec<i64>>(
                        image["powers"].clone(),
                    )?)?;
                    let weight = Atom::parse(
                        image["weight"]
                            .as_str()
                            .ok_or("checkpoint owner has no weight")?,
                        "qcd_ladder",
                        ParseSettings::default(),
                    )?;
                    let row = owner_roots
                        .entry(key)
                        .or_insert_with(|| (Atom::Zero, 0, BTreeSet::new()));
                    row.0 += weight;
                    row.1 += 1;
                    row.2.insert(index);
                }
            }
        }
    }
    let prepared = application
        .prepared
        .lock()
        .map_err(|error| error.to_string())?;
    let mut combined = Vec::new();
    for (key, (weight, multiplicity, contribution_indices)) in roots {
        let weight = super::super::combine_exact(weight);
        if weight.is_zero() {
            continue;
        }
        combined.push(
            json!({"powers": key.powers(), "weight": weight.to_canonical_string(),
            "admitted": prepared.routed_entry_is_admitted(&key)?, "multiplicity": multiplicity,
            "contributions": contribution_indices}),
        );
    }
    let mut combined_owner_roots = Vec::new();
    for (key, (weight, multiplicity, contribution_indices)) in owner_roots {
        let weight = super::super::combine_exact(weight);
        if weight.is_zero() {
            continue;
        }
        prepared.validate_entry(&key)?;
        combined_owner_roots.push(
            json!({"powers": key.powers(), "weight": weight.to_canonical_string(),
            "multiplicity": multiplicity, "contributions": contribution_indices}),
        );
    }
    let report = json!({"schema": "vakint.qcd-ladder.campaign-census.v1", "diagnostic_replay_only": true,
        "projection":projection,"engineering_degree_balance":engineering_degree_balance,
        "degree_semantics":"D_key-j is engineering mass dimension, not the UV degree of an individual lowered term with MUV fixed; omega_UV(d=4)=20-2*D_key",
        "fresh_numerator_preparation": false, "ibp_application": false, "source": source_path,
        "source_sha256": source_sha256, "source_bytes": source.len(), "data": load,
        "family_fingerprint": prepared.family().fingerprint(), "load_seconds": load_seconds,
        "canonical_seconds": canonical_seconds, "scalar_limits": format!("{:?}", settings.rustred_scalar_numerator_limits),
        "application_options": format!("{:?}", settings.rustred_application_options),
        "contributions": contributions, "combined_roots": combined, "combined_owner_roots": combined_owner_roots});
    std::fs::write(
        evidence.join("campaign-census.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    Ok(())
}

/// A single exact term from the physical census, not the complete benchmark.
/// Starting-predicate membership and surviving catalogue closure are separate.
#[test]
#[ignore = "manual bounded uncovered-witness reduction; requires witness JSON and local data"]
fn physical_uncovered_witness_safe_reduction() -> Result<()> {
    use crate::{
        EvaluationMethod, EvaluationOrder, OfficialExternalDataOptions, RustRedEvaluationOptions,
        RustRedMasterMode,
    };
    use serde_json::json;
    use std::time::Instant;

    let witness_path = PathBuf::from(
        std::env::var_os("VAKINT_QCD_COVERAGE_WITNESS").ok_or("set VAKINT_QCD_COVERAGE_WITNESS")?,
    );
    let data_path =
        PathBuf::from(std::env::var_os("VAKINT_EXTERNAL_DATA").ok_or("set VAKINT_EXTERNAL_DATA")?);
    let evidence = PathBuf::from(
        std::env::var_os("VAKINT_QCD_WEIGHT_EVIDENCE").ok_or("set VAKINT_QCD_WEIGHT_EVIDENCE")?,
    );
    std::fs::create_dir_all(&evidence)?;
    let bytes = std::fs::read(&witness_path)?;
    let witness: serde_json::Value = serde_json::from_slice(&bytes)?;
    let key = IntegralKey::try_new(serde_json::from_value::<Vec<i64>>(
        witness["powers"].clone(),
    )?)?;
    let weight = Atom::parse(
        witness["weight_original"]
            .as_str()
            .ok_or("witness has no exact weight")?,
        "qcd_ladder",
        ParseSettings::default(),
    )?;
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/physical-witness-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    };
    settings.rustred_application_options.unresolved_policy =
        rustred::solver::UnresolvedPolicy::Preserve;
    let started = Instant::now();
    eprintln!("qcd_uncovered_witness phase=data_loading_begin");
    let mut engine = Vakint::new()?;
    let loaded = engine.load_official_external_data(
        5,
        &OfficialExternalDataOptions {
            local_artifacts_path: Some(data_path),
            ..Default::default()
        },
    )?;
    engine.detach_cache(5)?;
    let application = engine.saved_application(&settings)?;
    let admitted = application
        .prepared
        .lock()
        .map_err(|error| error.to_string())?
        .routed_entry_is_admitted(&key)?;
    if admitted {
        return Err("the selected witness is no longer outside the actual saved predicates".into());
    }
    let source = super::admission_census_tests::input_for_key(&key)?;
    if !source.numerator.is_one() {
        return Err("the concrete witness must have no reconstructed auxiliary numerator".into());
    }
    let mass_squared = Atom::var(vk_symbol!("mUV")).pow(2);
    let integral = source.integral.replace_map(|part, _, output| {
        if let AtomView::Fun(function) = part
            && function.get_symbol() == S.prop
        {
            let arguments = function.iter().collect::<Vec<_>>();
            if arguments.len() == 5 {
                **output = S.prop.call((
                    arguments[0],
                    arguments[1],
                    arguments[2],
                    &mass_squared,
                    arguments[4],
                ));
            }
        }
    });
    let input = weight * integral;
    std::fs::write(
        evidence.join("input.symbolica"),
        input.to_canonical_string(),
    )?;
    let mut report = json!({"schema":"vakint.qcd-ladder.uncovered-witness.v1",
        "complete_primary_benchmark":false,"witness":witness_path,
        "witness_sha256":format!("{:x}",Sha256::digest(&bytes)),"powers":key.powers(),
        "actual_saved_predicate_admitted":admitted,"physical_mass_squared":mass_squared.to_canonical_string(),
        "weight_original":witness["weight_original"],"application_options":format!("{:?}",settings.rustred_application_options),
        "scalar_limits":format!("{:?}",settings.rustred_scalar_numerator_limits),
        "data":loaded,"data_loading_seconds":started.elapsed().as_secs_f64(),"status":"reduction_started"});
    std::fs::write(
        evidence.join("result.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    eprintln!(
        "qcd_uncovered_witness phase=reduction_begin powers={:?}",
        key.powers()
    );
    let started = Instant::now();
    let outcome = engine.evaluate_rustred_symbolic(&settings, input.as_view());
    report["reduction_seconds"] = json!(started.elapsed().as_secs_f64());
    match outcome {
        Err(error) => {
            report["status"] = json!("operational_or_algebra_error");
            report["error"] = json!(error.to_string());
            std::fs::write(
                evidence.join("result.json"),
                serde_json::to_vec_pretty(&report)?,
            )?;
            return Err(error.into());
        }
        Ok(result) => {
            let known = engine
                .master_symbols()?
                .into_iter()
                .map(|value| value.as_atom().clone())
                .collect::<BTreeSet<_>>();
            let mut terminals = BTreeSet::new();
            result.value.visitor(&mut |part| {
                if matches!(part,AtomView::Fun(function) if function.get_symbol()==vk_symbol!("RustRedMaster")) {
                    terminals.insert(part.to_owned());
                    return false;
                }
                true
            });
            let unknown = terminals
                .difference(&known)
                .map(AtomCore::to_canonical_string)
                .collect::<Vec<_>>();
            let value = result.value.to_canonical_string();
            std::fs::write(evidence.join("result.symbolica"), &value)?;
            report["status"] = json!("retained_result");
            report["result_sha256"] = json!(format!("{:x}", Sha256::digest(value.as_bytes())));
            report["terminal_count"] = json!(terminals.len());
            report["unknown_terminals"] = json!(unknown);
            report["catalogue_closed"] = json!(unknown.is_empty());
            let guards = result
                .nonzero_conditions
                .iter()
                .map(AtomCore::to_canonical_string)
                .collect::<Vec<_>>()
                .join("\n");
            std::fs::write(evidence.join("nonzero-conditions.symbolica"), &guards)?;
            report["guards"] = json!({"count":result.nonzero_conditions.len(),
                "sha256":format!("{:x}",Sha256::digest(guards.as_bytes())),"file":"nonzero-conditions.symbolica"});
        }
    }
    std::fs::write(
        evidence.join("result.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    println!("{report}");
    Ok(())
}
