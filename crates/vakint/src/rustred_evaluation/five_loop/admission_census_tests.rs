//! External-package admission only: never applies an IBP rule or prepares an
//! owner. The fixed inventory is a witness set, not an exhaustive point census.
use super::*;
use serde_json::{Value, json};
use std::time::Instant;

#[path = "../../topologies/five_loop_data.rs"]
#[allow(dead_code)]
mod census;

pub(super) struct ReconstructedInput {
    pub(super) parent_sector: u32,
    pub(super) integral: Atom,
    pub(super) numerator: Atom,
    pub(super) base: IntegralKey,
}

impl ReconstructedInput {
    pub(super) fn expression(&self) -> Atom {
        &self.numerator * &self.integral
    }
}

/// Reconstruct the actual common-basis integrand in a compiled defining
/// parent's original loop coordinates. All physical zero rows are retained,
/// and negative auxiliary slots enter the numerator rather than the graph.
pub(super) fn input_for_key(target: &IntegralKey) -> Result<ReconstructedInput, VakintError> {
    let invalid = |detail: String| VakintError::InvalidIntegralFormat(detail);
    if target.powers().len() != 15 {
        return Err(invalid(
            "the pinned admission fixture requires fifteen family slots".into(),
        ));
    }
    let parent = census::PARENTS
        .iter()
        .find(|parent| {
            target
                .powers()
                .iter()
                .enumerate()
                .all(|(slot, power)| *power <= 0 || parent.physical_slots.contains(&slot))
        })
        .ok_or_else(|| {
            invalid("no compiled defining parent contains this positive support".into())
        })?;
    let momenta = common_momentum_rows()
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .map(|(axis, coefficient)| Atom::num(*coefficient) * S.k.call(axis + 1))
                .sum::<Atom>()
        })
        .collect::<Vec<_>>();
    let base = IntegralKey::try_new(
        target
            .powers()
            .iter()
            .map(|power| (*power).max(0))
            .collect::<Vec<_>>(),
    )
    .map_err(|error| invalid(error.to_string()))?;
    let mut integral = Atom::one();
    for (physical, slot) in parent.physical_slots.iter().enumerate() {
        let signed_slot = (*slot + 1) as i64;
        let source = parent
            .vertices
            .iter()
            .position(|vertex| vertex.contains(&signed_slot))
            .ok_or_else(|| invalid("compiled parent omits a source incidence".into()))?;
        let target_vertex = parent
            .vertices
            .iter()
            .position(|vertex| vertex.contains(&-signed_slot))
            .ok_or_else(|| invalid("compiled parent omits a target incidence".into()))?;
        integral *= S.prop.call((
            physical + 1,
            S.edge.call((source + 1, target_vertex + 1)),
            &momenta[*slot],
            Atom::one(),
            Atom::num(base.powers()[*slot]),
        ));
    }
    let mut numerator = Atom::one();
    for (slot, power) in target.powers().iter().enumerate() {
        if *power < 0 {
            let exponent = u32::try_from(power.unsigned_abs()).map_err(|_| {
                invalid("fixture numerator exceeds native scalar exponent representation".into())
            })?;
            numerator *= (S.dot.call((&momenta[slot], &momenta[slot])) - Atom::one()).pow(exponent);
        }
    }
    Ok(ReconstructedInput {
        parent_sector: parent.sector,
        integral: S.topo.call(integral),
        numerator,
        base,
    })
}

fn check_identity(
    application: &SavedApplication,
    target: &IntegralKey,
    input: &ReconstructedInput,
) -> Result<(), String> {
    let prepared = application
        .prepared
        .lock()
        .map_err(|error| error.to_string())?;
    prepared
        .validate_entry(target)
        .map_err(|error| error.to_string())?;
    if !prepared
        .routed_entry_is_admitted(target)
        .map_err(|error| error.to_string())?
    {
        return Err(
            "original inventory key is not admitted by its native routed entry policy".into(),
        );
    }
    let scalar = FamilyScalarNumeratorService::try_new(
        prepared.family(),
        S.dot,
        (1..=5).map(|axis| S.k.call(axis)).collect(),
        Default::default(),
    )
    .map_err(|error| error.to_string())?;
    let lowered = scalar
        .lower_weighted(&input.numerator, &input.base, &Atom::one())
        .map_err(|error| error.to_string())?;
    if lowered.terms().len() != 1 {
        return Err(format!(
            "identity reconstruction produced {} keys instead of one",
            lowered.terms().len()
        ));
    }
    let term = &lowered.terms()[0];
    if term.integral() != target
        || term.coefficient() != &prepared.family().coefficient_context().one()
        || term.scalar_spectator() != &Atom::one()
        || term.common_mass_squared_power() != 0
    {
        return Err(format!(
            "identity reconstruction differs: key={:?}, coefficient={}, spectator={}, mass_power={}",
            term.integral().powers(),
            term.coefficient().to_expression(),
            term.scalar_spectator(),
            term.common_mass_squared_power()
        ));
    }
    Ok(())
}

#[test]
fn powered_reconstruction_keeps_compact_and_expanded_weighted_inputs_equivalent() {
    use rustred::algebra::CoefficientContext;
    use rustred::family::{AffineDenominator, IntegralFamily};

    Vakint::initialize_vakint_symbols();
    let context = CoefficientContext::try_new(["d"]).unwrap();
    let family = IntegralFamily::new(
        "admission_fixture_common_basis",
        (1..=5).map(|axis| format!("k{axis}")).collect(),
        vec![],
        context.clone(),
        context.parameter("d").unwrap(),
        common_momentum_rows()
            .iter()
            .map(|row| {
                AffineDenominator::new(
                    context.integer(-1),
                    (0..5)
                        .flat_map(|a| {
                            (a..5).map(move |b| row[a] * row[b] * if a == b { 1 } else { 2 })
                        })
                        .map(|value| context.integer(value))
                        .collect(),
                )
            })
            .collect(),
        vec![],
        vec![context.zero(); 15],
    )
    .unwrap();
    let target = IntegralKey::try_new([-2, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 7]).unwrap();
    let input = input_for_key(&target).unwrap();
    let original = input.numerator.clone();
    let scalar = FamilyScalarNumeratorService::try_new(
        &family,
        S.dot,
        (1..=5).map(|axis| S.k.call(axis)).collect(),
        Default::default(),
    )
    .unwrap();
    let expanded = input.numerator.expand();
    assert_ne!(original, expanded);
    for numerator in [&original, &expanded] {
        let lowered = scalar
            .lower_weighted(numerator, &input.base, &Atom::one())
            .unwrap();
        assert_eq!(lowered.terms().len(), 1);
        assert_eq!(lowered.terms()[0].integral(), &target);
        assert_eq!(lowered.terms()[0].coefficient(), &context.one());
        assert_eq!(lowered.terms()[0].scalar_spectator(), &Atom::one());
    }
    assert_eq!(
        input.numerator, original,
        "the frontend still receives the factored input"
    );
}

#[test_log::test]
#[ignore = "requires VAKINT_SAVED_APPLICATION and VAKINT_SAVED_APPLICATION_INVENTORY external generation-54 data"]
fn external_package_all_813_witnesses_pass_native_frontend_admission_without_reduction() {
    let path = std::env::var_os("VAKINT_SAVED_APPLICATION").expect("set external package path");
    let inventory = std::env::var_os("VAKINT_SAVED_APPLICATION_INVENTORY")
        .expect("set exact 813-witness inventory path");
    let document: Value = serde_json::from_slice(&std::fs::read(inventory).unwrap()).unwrap();
    assert_eq!(document["schema"], "rustred.saved-application-witnesses.v1");
    assert_eq!(document["arity"], 15);
    let rows = document["cases"].as_array().expect("inventory cases");
    assert_eq!(
        rows.len(),
        813,
        "this test must exercise the complete frozen witness inventory"
    );
    let mut settings = VakintSettings {
        evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
            crate::RustRedEvaluationOptions {
                substitute_masters: false,
            },
        )]),
        form_exe_path: "/admission-census-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    };
    if let Ok(workers) = std::env::var("VAKINT_SAVED_APPLICATION_PREPARATION_WORKERS") {
        settings.rustred_application_options.preparation_workers = workers
            .parse()
            .expect("positive native preparation worker count");
    }
    let native_operations_override =
        std::env::var("VAKINT_SAVED_APPLICATION_NATIVE_POLYNOMIAL_OPERATIONS")
            .ok()
            .map(|value| {
                value
                    .parse::<usize>()
                    .expect("native polynomial operation allowance must be an unsigned integer")
            });
    let transport_profile = std::env::var("VAKINT_SAVED_APPLICATION_TRANSPORT_PROFILE").ok();
    assert!(
        transport_profile.is_none() || native_operations_override.is_none(),
        "select either a complete transport profile or the single native-operation override"
    );
    if let Some(profile) = &transport_profile {
        assert_eq!(
            profile,
            test_profiles::MATERIALIZATION_PROFILE,
            "unknown test-only transport profile"
        );
        test_profiles::apply_materialization_profile(&mut settings);
    } else if let Some(operations) = native_operations_override {
        settings
            .rustred_application_options
            .trace_limits
            .expansion
            .max_native_polynomial_operations = operations;
    }
    println!(
        "{}",
        json!({"event":"admission_census_policy", "cases":rows.len(),
        "transport_profile":transport_profile,
        "application_options":format!("{:?}",settings.rustred_application_options),
        "scalar_limits":format!("{:?}",settings.rustred_scalar_numerator_limits),
        "matching_limits":format!("{:?}",settings.rustred_momentum_match_limits),
        "production_defaults_changed":false,"ibp_application":false})
    );
    let load_started = Instant::now();
    let vakint = Vakint::new()
        .unwrap()
        .with_rustred_application(&mut settings, PathBuf::from(path))
        .unwrap();
    let application = vakint.saved_application(&settings).unwrap();
    {
        let prepared = application.prepared.lock().unwrap();
        assert_eq!(prepared.metadata()["source_generation"], 54);
        assert_eq!(prepared.metadata()["source_query_count"], 804);
        assert_eq!(
            prepared.owner_preparation_statistics().preparation_attempts,
            0
        );
        assert_eq!(prepared.cache_statistics().cached_integrals(), 0);
    }
    println!(
        "{}",
        json!({"event":"admission_census_prepared", "seconds":load_started.elapsed().as_secs_f64(),
        "cases":rows.len(), "ibp_application":false, "exhaustive_point_enumeration":false,
        "native_polynomial_operations":settings.rustred_application_options.trace_limits.expansion.max_native_polynomial_operations,
        "native_polynomial_operations_override":native_operations_override.is_some(),
        "transport_profile":transport_profile})
    );
    let started = Instant::now();
    let mut failures = Vec::new();
    for row in rows {
        let id = row["id"].as_str().expect("inventory case ID");
        let target = IntegralKey::try_new(
            serde_json::from_value::<Vec<i64>>(row["powers"].clone()).unwrap(),
        )
        .unwrap();
        let case_started = Instant::now();
        let mut stage = "input_reconstruction";
        let mut parent = None;
        let mut diagnostics = Value::Null;
        let outcome = (|| -> Result<Vec<(IntegralKey, Atom)>, String> {
            let input = input_for_key(&target).map_err(|error| error.to_string())?;
            parent = Some(input.parent_sector);
            diagnostics = json!({"input":input.expression().to_string()});
            stage = "known_parent_identity_lowering";
            check_identity(&application, &target, &input)?;
            stage = "frontend_canonicalization";
            let mut expression = VakintExpression::try_from(input.expression().as_view())
                .map_err(|error| error.to_string())?;
            expression
                .canonicalize(&settings, &vakint.topologies, false)
                .map_err(|error| error.to_string())?;
            stage = "frontend_native_tensor_prepass";
            expression
                .tensor_reduce(&vakint, &settings)
                .map_err(|error| error.to_string())?;
            diagnostics["scalar_expression"] = json!(expression.to_string());
            let terms = VakintExpression::split_integrals(Atom::from(expression).as_view())
                .map_err(|error| error.to_string())?;
            let mut admitted = Vec::new();
            for term in terms {
                stage = "native_admitted_witness_search";
                let matched = vakint
                    .match_rustred_input(
                        &settings,
                        term.integral.as_view(),
                        term.numerator.as_view(),
                    )
                    .map_err(|error| error.to_string())?
                    .ok_or_else(|| "native matcher returned no topology witness".to_owned())?;
                // A different proved loop basis may legitimately permute or
                // expand this key. Check the selected verified witness rather
                // than demanding equality to the original presentation.
                stage = "selected_witness_admission";
                diagnostics["routing_witness"] = json!(matched.to_string());
                let selected = vakint
                    .rustred_matched_input_key(&settings, term.numerator.as_view(), &matched)
                    .map_err(|error| error.to_string())?
                    .ok_or_else(|| {
                        "selected native witness is outside the saved query union".to_owned()
                    })?;
                admitted.extend(selected);
            }
            if admitted.is_empty() {
                return Err(
                    "nonzero identity reconstruction left no admitted frontend entries".into(),
                );
            }
            Ok(admitted)
        })();
        match outcome {
            Ok(selected) => println!(
                "{}",
                json!({"event":"admission_census_passed", "id":id,
                "parent":parent, "seconds":case_started.elapsed().as_secs_f64(),
                "identity_reconstruction_exact":true, "selected_key_count":selected.len(),
                "selected_key_prefix":selected.iter().take(8).map(|(key, coefficient)| json!({"powers":key.powers(),"weight":coefficient.to_string()})).collect::<Vec<_>>(),
                "ibp_application":false})
            ),
            Err(error) => {
                println!(
                    "{}",
                    json!({"event":"admission_census_failed", "case":row,
                    "parent":parent, "stage":stage, "seconds":case_started.elapsed().as_secs_f64(),
                    "error":error, "diagnostics":diagnostics, "ibp_application":false})
                );
                failures.push(format!("{id} at {stage}: {error}"));
            }
        }
    }
    let prepared = application.prepared.lock().unwrap();
    let owners = prepared.owner_preparation_statistics();
    let cache = prepared.cache_statistics();
    assert_eq!(
        owners.preparation_attempts, 0,
        "admission must not prepare owner formulas"
    );
    assert_eq!(
        cache.cached_integrals(),
        0,
        "admission must not create coefficient cache entries"
    );
    assert_eq!(
        cache.rule_applications(),
        0,
        "admission must not apply IBP rules"
    );
    println!(
        "{}",
        json!({"event":"admission_census_finished", "checked":rows.len(),
        "passed":rows.len()-failures.len(), "failed":failures.len(), "seconds":started.elapsed().as_secs_f64(),
        "owner_preparation_attempts":owners.preparation_attempts,"cached_integrals":cache.cached_integrals(),
        "ibp_application":false,"exhaustive_point_enumeration":false,"closure_claim":false})
    );
    assert!(
        failures.is_empty(),
        "native frontend admission failures: {}",
        failures.join("; ")
    );
}
