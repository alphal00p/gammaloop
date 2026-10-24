//! Public retained/substituted behavior against the authenticated 640-master campaign.
use super::*;
use crate::{EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions};

fn engine() -> (Vakint, VakintSettings) {
    let path = super::super::cache_tests::external_path_through(5);
    let mut engine = Vakint::new().unwrap();
    engine.load_external_data(path, 5).unwrap();
    assert_eq!(
        engine
            .rustred_data()
            .unwrap()
            .masters
            .masters
            .values()
            .filter(|master| master.loop_count == 5)
            .count(),
        640,
        "these controls require the authenticated 640-master campaign"
    );
    let settings = VakintSettings {
        epsilon_symbol: "ep".into(),
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Retain,
            },
        )]),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        form_exe_path: "/five-loop-partial-test-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    };
    (engine, settings)
}

#[test_log::test]
fn outside_recorded_inputs_keep_safe_partial_reduction_semantics() {
    let (engine, settings) = engine();
    outside_recorded_degree(&engine, &settings);
    outside_recorded_rank(&engine, &settings);
    bare_application_registry_view(&engine, &settings);
}

fn bare_application_registry_view(engine: &Vakint, settings: &VakintSettings) {
    // Reuse this genuinely authenticated native program while removing only
    // the external master-data context. This isolates the read-only registry
    // view without repeating the expensive native import or changing rules.
    let mut bare = engine.clone();
    bare.rustred_data = None;
    let mut source = engine
        .rustred_application
        .lock()
        .unwrap()
        .as_ref()
        .unwrap()
        .clone();
    source.bare = true;
    bare.rustred_application = Arc::new(Mutex::new(Some(source)));
    let declarations = bare.announced_bare_masters().unwrap();
    let registry = bare.rustred_data().unwrap();
    let declared = declarations
        .keys()
        .find(|symbol| !registry.masters.masters.contains_key(*symbol))
        .unwrap()
        .clone();
    let requested = crate::MasterSymbol::try_from(declared.clone()).unwrap();
    assert!(bare.master_symbols().unwrap().contains(&requested));
    let error = bare.get_master(&requested, 16).unwrap_err();
    assert!(matches!(
        &error,
        VakintError::Inspection(crate::inspection::InspectionError::UnavailableGraph { .. })
    ));
    assert!(error.to_string().contains("matching external-data pack"));
    assert_eq!(
        bare.finalize_rustred(
            declared.clone(),
            settings,
            5,
            crate::RustRedMasterMode::Retain
        )
        .unwrap(),
        declared
    );
    for mode in [
        crate::RustRedMasterMode::Symbolic,
        crate::RustRedMasterMode::Numerical,
    ] {
        let error = bare
            .finalize_rustred(declared.clone(), settings, 5, mode)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("declared terminal has no computed coefficients"),
            "{error}"
        );
        assert!(!error.contains("undeclared residual master"), "{error}");
    }
    let error = bare
        .numerical_evaluation(
            settings,
            declared.as_view(),
            &Default::default(),
            &Default::default(),
            None,
        )
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("declared terminal has no computed coefficients"),
        "{error}"
    );
    assert!(engine.announced_bare_masters().unwrap().is_empty());
}

fn outside_recorded_degree(engine: &Vakint, settings: &VakintSettings) {
    // Five independent massive tadpoles, total degree 11. The exact entry
    // predicates, not this counting argument, establish that it is uncovered.
    let key = IntegralKey::try_new([0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 7]).unwrap();
    let application = engine.saved_application(settings).unwrap();
    assert!(
        !application
            .prepared
            .lock()
            .unwrap()
            .routed_entry_is_admitted(&key)
            .unwrap()
    );
    let input = super::admission_census_tests::input_for_key(&key)
        .unwrap()
        .expression();
    let result = engine
        .evaluate_rustred_symbolic(settings, input.as_view())
        .unwrap();
    assert!(!result.value.is_zero());
    let first = engine
        .evaluate_integral(
            settings,
            crate::vakint_parse!("topo(I1L(1,7))").unwrap().as_view(),
        )
        .unwrap();
    let other = engine
        .evaluate_integral(
            settings,
            crate::vakint_parse!("topo(I1L(1,1))").unwrap().as_view(),
        )
        .unwrap();
    let expected = first * other.pow(4);
    assert!(
        (&result.value - &expected).expand().together().is_zero(),
        "actual: {}; independent tadpole product: {}",
        result.value,
        expected
    );
    // For this factorized control the generic tadpole recurrence closes even
    // beyond the recorded starting domain. No new coverage is claimed.
    let data = engine.rustred_data().unwrap();
    result.value.visitor(&mut |part| {
        if matches!(part, AtomView::Fun(function) if function.get_symbol() == vk_symbol!("RustRedMaster")) {
            assert!(data.masters.masters.contains_key(&part.to_owned()), "unexpected residual {part}");
            return false;
        }
        true
    });
    let warm = engine
        .evaluate_rustred_symbolic(settings, input.as_view())
        .unwrap();
    assert_eq!(result.value, warm.value);
    assert_eq!(result.nonzero_conditions, warm.nonzero_conditions);
    assert!(
        !application
            .prepared
            .lock()
            .unwrap()
            .routed_entry_is_admitted(&key)
            .unwrap()
    );
}

fn outside_recorded_rank(engine: &Vakint, settings: &VakintSettings) {
    let key = IntegralKey::try_new([-7, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1]).unwrap();
    let application = engine.saved_application(settings).unwrap();
    assert!(
        !application
            .prepared
            .lock()
            .unwrap()
            .routed_entry_is_admitted(&key)
            .unwrap()
    );
    let input = super::admission_census_tests::input_for_key(&key)
        .unwrap()
        .expression();
    let result = engine
        .evaluate_rustred_symbolic(settings, input.as_view())
        .unwrap();
    assert!(!result.value.is_zero());
    let data = engine.rustred_data().unwrap();
    result.value.visitor(&mut |part| {
        if matches!(part, AtomView::Fun(function) if function.get_symbol() == vk_symbol!("RustRedMaster")) {
            assert!(data.masters.masters.contains_key(&part.to_owned()), "unexpected residual {part}");
            return false;
        }
        true
    });
    // This physical rank-seven request closes through the saved rules. Test
    // catalogue diagnostics independently with a qualified, absent terminal;
    // the native incomplete-owner tests exercise how such leaves are produced.
    let unknown = data.masters.symbol(
        &application.fingerprint,
        &[-8, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1],
    );
    assert!(!data.masters.masters.contains_key(&unknown));
    assert_eq!(
        engine
            .finalize_rustred(
                unknown.clone(),
                settings,
                5,
                crate::RustRedMasterMode::Retain
            )
            .unwrap(),
        unknown
    );
    for mode in [
        crate::RustRedMasterMode::Symbolic,
        crate::RustRedMasterMode::Numerical,
    ] {
        let error = engine
            .finalize_rustred(unknown.clone(), settings, 5, mode)
            .unwrap_err()
            .to_string();
        assert!(error.contains("undeclared residual master"), "{error}");
        assert!(error.contains(&data.masters.version), "{error}");
        assert!(error.contains(".dot"), "{error}");
        assert!(error.contains(&unknown.to_canonical_string()), "{error}");
        let zero = crate::vakint_parse!("ep^-2*((1-ep)^-1-1-ep)-(1-ep)^-1").unwrap();
        assert!(
            engine
                .finalize_rustred(zero * &unknown, settings, 5, mode)
                .unwrap()
                .is_zero()
        );
    }
    // Explicit scalar diagnostics stay errors after a successful warm native
    // request; preserving coverage misses must never suppress that policy.
    let mut limited = settings.clone();
    limited.rustred_scalar_numerator_limits.max_input_nodes = 0;
    limited.rustred_scalar_numerator_limits =
        limited.rustred_scalar_numerator_limits.with_diagnostics();
    let error = engine
        .evaluate_rustred_symbolic(&limited, input.as_view())
        .unwrap_err()
        .to_string();
    assert!(error.contains("limit"), "{error}");
}
