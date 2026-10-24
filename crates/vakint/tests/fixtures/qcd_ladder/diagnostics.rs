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

/// Replays an explicitly supplied prepared integrand for an early campaign
/// census. Fresh native preparation is still required by the slow acceptance.
#[test]
#[ignore = "manual selected-witness census; requires prepared input and explicit local external data"]
fn prepared_ladder_selected_campaign_census_without_ibps() -> Result<()> {
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
    let mass = Atom::var(vk_symbol!("mUV"));
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
        let started = Instant::now();
        eprintln!("qcd_campaign_census contribution={index} phase=selected_witness_begin");
        let mut matched = engine
            .match_rustred_input(
                &settings,
                original.integral.as_view(),
                original.numerator.as_view(),
            )?
            .ok_or("quadratic source has no verified geometric witness")?;
        let matching_seconds = started.elapsed().as_secs_f64();
        matched.apply_replacement_rules()?;
        let mut term = original;
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
            settings.rustred_scalar_numerator_limits,
        )?;
        capture_weighted_input(&evidence, index, &input, prepared.family())?;
        let started = Instant::now();
        let lowered =
            scalar.lower_weighted(&input.numerator, &input.target, &input.mass_squared)?;
        let lowering_seconds = started.elapsed().as_secs_f64();
        let mut records = Vec::new();
        for lowered_term in lowered.terms() {
            let key = lowered_term.integral();
            let weight = super::super::combine_exact(
                lowered_term.coefficient().to_expression() * lowered_term.scalar_spectator(),
            );
            let degree = key.powers().iter().sum::<i64>();
            let mut mass_powers = Vec::new();
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
            let admitted = prepared.routed_entry_is_admitted(key)?;
            if !admitted {
                eprintln!(
                    "qcd_campaign_census uncovered contribution={index} family={} powers={:?} mass_squared_powers={mass_powers:?}",
                    prepared.family().name(),
                    key.powers()
                );
            }
            let owner_route = if admitted {
                let routed = prepared.route_entry_to_owner(key)?;
                let mut images = Vec::new();
                for (endpoint, coefficient) in routed.terms() {
                    if coefficient.is_zero() {
                        continue;
                    }
                    let owner_weight =
                        super::super::combine_exact(&weight * coefficient.to_expression());
                    images.push(json!({"powers": endpoint.powers(), "weight": owner_weight.to_canonical_string(),
                        "transport_coefficient": coefficient.to_expression().to_canonical_string()}));
                }
                Some(
                    json!({"images": images, "nonzero_conditions": routed.nonzero_conditions().iter()
                    .map(|coefficient| coefficient.to_expression().to_canonical_string()).collect::<Vec<_>>()}),
                )
            } else {
                None
            };
            records.push(json!({"powers": key.powers(), "weight": weight.to_canonical_string(),
                "D_key": degree,"r_idx":key.powers().iter().filter(|power|**power<0).map(|power|-power).sum::<i64>(),
                "omega_uv_at_d4_with_mass_fixed":20-2*degree,
                "admitted": admitted, "mass_squared_powers": mass_powers, "owner_route": owner_route}));
        }
        let row = json!({"contribution": index, "source_identity":identity,"binding":binding,"topology": matched.canonical_topology.get_integral().name,
            "family_fingerprint": prepared.family().fingerprint(), "target": input.target.powers(),
            "mass_squared": input.mass_squared.to_canonical_string(), "matching_seconds": matching_seconds,
            "lowering_seconds": lowering_seconds, "roots": records,
            "nonzero_conditions": lowered.nonzero_conditions().iter().map(AtomCore::to_canonical_string).collect::<Vec<_>>()});
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
