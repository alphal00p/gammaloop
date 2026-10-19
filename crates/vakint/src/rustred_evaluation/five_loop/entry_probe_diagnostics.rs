//! Manual two-case eligibility diagnosis. No IBPs or new admission authority.
//! The rank-zero-order hypothesis below is retained as negative diagnostic
//! evidence: it fails on gen54's singleton high-rank candidate, as recorded in
//! five_loop_application.typ. It is not an acceptance gate or a slow-CI test.
use super::*;
use serde_json::json;

#[test_log::test]
#[ignore = "manual two-case resource-profile diagnosis; requires explicit profile and external package"]
fn coherent_transport_profile_checks_two_cases_without_reduction() {
    assert_eq!(
        std::env::var("VAKINT_ENTRY_PROBE_PROFILE").as_deref(),
        Ok(test_profiles::MATERIALIZATION_PROFILE),
        "explicitly select VAKINT_ENTRY_PROBE_PROFILE=materialize-512m"
    );
    let path = std::env::var_os("VAKINT_SAVED_APPLICATION")
        .expect("set explicit generation-54 application path");
    let mut settings = VakintSettings {
        evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
            crate::RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/entry-resource-profile-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    };
    let defaults = settings.rustred_application_options;
    test_profiles::apply_materialization_profile(&mut settings);
    assert_eq!(
        settings
            .rustred_application_options
            .trace_limits
            .expansion
            .exact_algebra,
        defaults.trace_limits.expansion.exact_algebra
    );
    assert_eq!(
        settings.rustred_application_options.reduction_limits,
        defaults.reduction_limits
    );
    assert_eq!(settings.rustred_scalar_numerator_limits, Default::default());
    println!(
        "{}",
        json!({"event":"entry_resource_profile_selected", "profile":test_profiles::MATERIALIZATION_PROFILE,
        "application_options":format!("{:?}", settings.rustred_application_options),
        "scalar_limits":format!("{:?}", settings.rustred_scalar_numerator_limits),
        "matching_limits":format!("{:?}", settings.rustred_momentum_match_limits),
        "production_defaults_changed":false,"changed_production_ordering":false,"ibp_application":false})
    );
    let started = std::time::Instant::now();
    let vakint = Vakint::new()
        .unwrap()
        .with_rustred_application(&mut settings, PathBuf::from(path))
        .unwrap();
    let application = vakint.saved_application(&settings).unwrap();
    println!(
        "{}",
        json!({"event":"entry_resource_profile_prepared", "load_seconds":started.elapsed().as_secs_f64()})
    );
    let mut failures = Vec::new();
    let cases = [("point-00352", 13, 10), ("point-00356", 14, 11)];
    for (id, power, rank) in cases {
        let target =
            IntegralKey::try_new([0, power, 1, 0, 0, 1, 0, 0, 0, 0, -rank, 1, 1, 1, 1]).unwrap();
        let started = std::time::Instant::now();
        println!(
            "{}",
            json!({"event":"entry_resource_profile_case_started", "id":id,"powers":target.powers()})
        );
        let result = (|| -> Result<Vec<(IntegralKey, Atom)>, VakintError> {
            assert!(
                application
                    .prepared
                    .lock()
                    .unwrap()
                    .routed_entry_is_admitted(&target)
                    .unwrap()
            );
            let source = admission_census_tests::input_for_key(&target)?;
            let mut expression = VakintExpression::try_from(source.expression().as_view())?;
            expression.canonicalize(&settings, &vakint.topologies, false)?;
            expression.tensor_reduce(&vakint, &settings)?;
            let terms = VakintExpression::split_integrals(Atom::from(expression).as_view())?;
            assert_eq!(terms.len(), 1);
            let term = &terms[0];
            // This is the actual production matching path, including its
            // original predicate order and deterministic minimum witness.
            let matched = vakint
                .match_rustred_input(&settings, term.integral.as_view(), term.numerator.as_view())?
                .expect("the finite corpus has a native matching witness");
            let selected = vakint
                .rustred_matched_input_key(&settings, term.numerator.as_view(), &matched)?
                .expect("the selected native witness must be admitted");
            assert!(!selected.is_empty());
            Ok(selected)
        })();
        match result {
            Ok(selected) => println!(
                "{}",
                json!({"event":"entry_resource_profile_case_passed","id":id,
                "seconds":started.elapsed().as_secs_f64(),"selected_key_count":selected.len(),
                "selected_keys":selected.iter().map(|(key,weight)|json!({"powers":key.powers(),"weight":weight.to_canonical_string()})).collect::<Vec<_>>(),"ibp_application":false})
            ),
            Err(error) => {
                println!(
                    "{}",
                    json!({"event":"entry_resource_profile_case_failed","id":id,
                    "seconds":started.elapsed().as_secs_f64(),"error":error.to_string(),"ibp_application":false})
                );
                // Record both cases, but fail the diagnostic after the census;
                // this does not turn an error into an admission decision.
                failures.push((id, error.to_string()));
            }
        }
    }
    let prepared = application.prepared.lock().unwrap();
    assert_eq!(
        prepared.owner_preparation_statistics().preparation_attempts,
        0
    );
    assert_eq!(prepared.cache_statistics().cached_integrals(), 0);
    assert_eq!(prepared.cache_statistics().rule_applications(), 0);
    println!(
        "{}",
        json!({"event":"entry_resource_profile_finished", "checked":cases.len(),
        "failed":failures.len(),"owner_preparation_attempts":0,"cached_integrals":0,"rule_applications":0})
    );
    assert!(
        failures.is_empty(),
        "resource-profile failures: {failures:?}"
    );
}

#[test_log::test]
#[ignore = "manual native entry-probe diagnosis; requires explicit external generation-54 package"]
fn coalesced_rank_zero_endpoint_rejects_expensive_bad_witnesses() {
    let path = std::env::var_os("VAKINT_SAVED_APPLICATION")
        .expect("set explicit generation-54 application path");
    let mut settings = VakintSettings {
        evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
            crate::RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/entry-probe-diagnostic-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    };
    // Reproduce the completed 813 diagnostic's explicit operation allowance.
    // All other entry-transport limits remain unchanged and are reported below.
    settings
        .rustred_application_options
        .trace_limits
        .expansion
        .max_native_polynomial_operations = 512_000_000;
    let started = std::time::Instant::now();
    let vakint = Vakint::new()
        .unwrap()
        .with_rustred_application(&mut settings, PathBuf::from(path))
        .unwrap();
    let application = vakint.saved_application(&settings).unwrap();
    println!(
        "{}",
        json!({"event":"entry_probe_diagnostic_prepared", "load_seconds":started.elapsed().as_secs_f64(),
        "trace_limits":format!("{:?}",settings.rustred_application_options.trace_limits),
        "ibp_application":false, "changed_production_ordering":false})
    );

    let mut total_old_errors = 0;
    for (id, power, rank) in [("point-00352", 13, 10), ("point-00356", 14, 11)] {
        let target =
            IntegralKey::try_new([0, power, 1, 0, 0, 1, 0, 0, 0, 0, -rank, 1, 1, 1, 1]).unwrap();
        let source = admission_census_tests::input_for_key(&target).unwrap();
        let mut expression = VakintExpression::try_from(source.expression().as_view()).unwrap();
        expression
            .canonicalize(&settings, &vakint.topologies, false)
            .unwrap();
        expression.tensor_reduce(&vakint, &settings).unwrap();
        let terms = VakintExpression::split_integrals(Atom::from(expression).as_view()).unwrap();
        assert_eq!(terms.len(), 1);
        let term = &terms[0];
        let mut candidates = 0;
        let mut old_errors = 0;
        let matched = vakint
            .topologies
            .match_topologies_to_user_input_with_acceptance(
                term.integral.as_view(),
                false,
                settings.rustred_momentum_match_limits,
                |raw| -> Result<Option<Vec<(IntegralKey, Atom)>>, VakintError> {
                    candidates += 1;
                    let candidate = candidates;
                    let mut matched = raw.clone();
                    matched.apply_replacement_rules()?;
                    let mut candidate_term = term.clone();
                    candidate_term.apply_scalar_numerator_replacement_rules(&matched, &settings)?;
                    let descriptor = select(&matched.canonical_topology)?;
                    let input = MatchedScalarInput::prepare(
                        candidate_term.numerator.as_view(),
                        &matched,
                        &descriptor.momenta,
                        descriptor.physical_slots,
                        15,
                        &settings,
                        ScalarNumeratorForm::PreserveProducts,
                    )?;
                    let prepared = application.prepared.lock().map_err(failure)?;
                    let service = FamilyScalarNumeratorService::try_new(
                        prepared.family(),
                        S.dot,
                        (1..=5).map(|axis| S.k.call(axis)).collect(),
                        settings.rustred_scalar_numerator_limits,
                    )
                    .map_err(failure)?;
                    let lowered = service
                        .lower_weighted(&input.numerator, &input.target, &input.mass_squared)
                        .map_err(failure)?;
                    println!("{}",json!({"event":"entry_probe_candidate","id":id,"candidate":candidate,
                        "target":input.target.powers(),"lowered_count":lowered.terms().len(),
                        "rank_zero_count":lowered.terms().iter().filter(|entry|entry.integral().powers().iter().all(|power|*power>=0)).count()}));
                    // This test changes only conjunction scheduling over the
                    // ALREADY coalesced native result. It invents no endpoint.
                    for entry in lowered
                        .terms()
                        .iter()
                        .filter(|entry| entry.integral().powers().iter().all(|power| *power >= 0))
                    {
                        assert!(!entry.coefficient().is_zero());
                        let admitted = prepared
                            .routed_entry_is_admitted(entry.integral())
                            .map_err(failure)?;
                        println!(
                            "{}",
                            json!({"event":"entry_probe_rank_zero", "id":id,"candidate":candidate,
                            "target":input.target.powers(),"lowered_count":lowered.terms().len(),
                            "powers":entry.integral().powers(),"coefficient":entry.coefficient().to_expression().to_canonical_string(),
                            "spectator":entry.scalar_spectator().to_canonical_string(),"admitted":admitted})
                        );
                        if !admitted {
                            // Independently record what the original order
                            // encounters. Its error is NOT the rejection proof:
                            // that is the exact false predicate immediately above.
                            for original in lowered.terms() {
                                match prepared.routed_entry_is_admitted(original.integral()) {
                                    Ok(true) => {}
                                    Ok(false) => {
                                        println!("{}",json!({"event":"entry_probe_original_order_outside","id":id,"candidate":candidate,"powers":original.integral().powers()}));
                                        break;
                                    }
                                    Err(error) => {
                                        old_errors += 1;
                                        println!("{}",json!({"event":"entry_probe_original_order_error","id":id,"candidate":candidate,"powers":original.integral().powers(),"error":error.to_string()}));
                                        break;
                                    }
                                }
                            }
                            return Ok(None);
                        }
                    }
                    // No cheap rejection: execute the complete unchanged
                    // predicate for every key. No error becomes false here.
                    let mut key = Vec::with_capacity(lowered.terms().len());
                    for entry in lowered.terms() {
                        let admitted = prepared
                            .routed_entry_is_admitted(entry.integral())
                            .map_err(|error| {
                                println!("{}",json!({"event":"entry_probe_unresolved_error","id":id,
                                    "candidate":candidate,"target":input.target.powers(),
                                    "powers":entry.integral().powers(),"error":error.to_string()}));
                                failure(error)
                            })?;
                        if !admitted {
                            return Ok(None);
                        }
                        key.push((
                            entry.integral().clone(),
                            (entry.coefficient().to_expression() * entry.scalar_spectator())
                                .together(),
                        ));
                    }
                    key.sort();
                    Ok(Some(key))
                },
            )
            .unwrap()
            .expect("an admitted alternative must remain available");
        let accepted = vakint
            .rustred_matched_input_key(&settings, term.numerator.as_view(), &matched)
            .unwrap()
            .expect("selected alternative must pass the ordinary complete predicate");
        assert!(!accepted.is_empty());
        total_old_errors += old_errors;
        println!(
            "{}",
            json!({"event":"entry_probe_diagnostic_case_finished","id":id,"candidates":candidates,"original_order_errors":old_errors,"selected_key_count":accepted.len(),"selected_keys":accepted.iter().map(|(key,weight)|json!({"powers":key.powers(),"weight":weight.to_canonical_string()})).collect::<Vec<_>>(),"ibp_application":false})
        );
    }
    assert!(
        total_old_errors > 0,
        "the diagnostic must reproduce an original-order failure"
    );
    let prepared = application.prepared.lock().unwrap();
    assert_eq!(
        prepared.owner_preparation_statistics().preparation_attempts,
        0
    );
    assert_eq!(prepared.cache_statistics().cached_integrals(), 0);
    assert_eq!(prepared.cache_statistics().rule_applications(), 0);
}
