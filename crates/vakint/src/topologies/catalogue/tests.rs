use super::*;
use crate::topologies::TopologyContractions;
use crate::utils::vakint_macros::vk_parse;
use crate::{
    EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions, Vakint, VakintExpression,
    VakintSettings,
};

#[test]
fn absent_five_loop_catalogue_preserves_the_legacy_unknown_contract() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        allow_unknown_integrals: true,
        ..Default::default()
    };
    // The Unknown lane describes formal denominators without authenticating
    // their incidence as a vacuum graph. These are five independent tadpoles
    // with deliberately nonconserving, acyclic incidence metadata.
    let input = vk_parse!(
        "topo(prop(1,edge(1,2),k(1),m,1)*prop(2,edge(2,3),k(2),m,1)*prop(3,edge(3,4),k(3),m,1)*prop(4,edge(4,5),k(4),m,1)*prop(5,edge(5,6),k(5),m,1))"
    ).unwrap();
    let expected = input
        .replace(vk_parse!("topo(x_)").unwrap().to_pattern())
        .with(vk_parse!("topo(UNKNOWN(x_))").unwrap().to_pattern());
    assert_eq!(
        vakint
            .to_canonical(&settings, input.as_view(), false)
            .unwrap(),
        expected
    );
    let mut calls = 0;
    let matched = vakint
        .topologies
        .match_topologies_to_user_input_with_acceptance(
            input.as_view(),
            true,
            Default::default(),
            |_| {
                calls += 1;
                Ok(Some(()))
            },
        )
        .unwrap()
        .unwrap();
    assert!(matches!(matched.canonical_topology, Topology::Unknown(_)));
    assert_eq!(calls, 0);
    assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
    assert!(vakint.rustred_application.lock().unwrap().is_none());
    vakint
        .topologies
        .prepare_five_loop_catalogue_for_tests()
        .unwrap();
    assert!(matches!(
        vakint
            .topologies
            .match_topologies_to_user_input(input.as_view(), true),
        Err(VakintError::InvalidIntegralFormat(_))
    ));
}

#[test]
fn default_registry_and_lower_entrypoints_do_not_prepare_saved_catalogue() {
    let vakint = Vakint::new().unwrap();
    let before = vakint.topologies.to_string();
    assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
    assert!(
        !vakint
            .topologies
            .iter()
            .any(|t| matches!(t, Topology::FiveLoop(_)))
    );
    let lower = vk_parse!("topo(I1L(mass_squared,1))").unwrap();
    let settings = VakintSettings {
        rustred_artifact_path: Some("/nonexistent-unused-five-loop-application".into()),
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/lower-native-test-must-not-use-form".into(),
        ..Default::default()
    };
    let canonical = vakint
        .to_canonical(&settings, lower.as_view(), false)
        .unwrap();
    assert!(
        vakint
            .topologies
            .match_topologies_to_user_input(canonical.as_view(), false)
            .unwrap()
            .is_some()
    );
    assert_eq!(vakint.topologies.to_string(), before);
    assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
    assert!(vakint.rustred_application.lock().unwrap().is_none());
    let mut expression = VakintExpression::try_from(lower.as_view()).unwrap();
    expression
        .canonicalize(&settings, &vakint.topologies, false)
        .unwrap();
    let _ = vakint
        .evaluate_integral(&settings, lower.as_view())
        .unwrap();
    let _ = vakint.evaluate(&settings, lower.as_view()).unwrap();
    let _ = vakint.tensor_reduce(&settings, lower.as_view()).unwrap();
    assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
    assert!(vakint.rustred_application.lock().unwrap().is_none());
}

#[test]
fn failed_settings_and_explicit_loads_never_publish_catalogue() {
    let vakint = Vakint::new().unwrap();
    let mut settings = VakintSettings {
        rustred_artifact_path: Some("/nonexistent-lazy-five-loop-application".into()),
        ..Default::default()
    };
    let input = vk_parse!("topo(I5L_32745(1,1,1,1,1,1,1,1,1,1,1,1,1))").unwrap();
    assert!(
        vakint
            .to_canonical(&settings, input.as_view(), false)
            .is_err()
    );
    assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
    assert!(vakint.rustred_application.lock().unwrap().is_none());
    assert!(
        vakint
            .clone()
            .with_rustred_application(&mut settings, "/nonexistent-explicit-five-loop-application")
            .is_err()
    );
    assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
}

#[test]
fn explicit_native_census_preparation_is_shared_and_preserves_custom_registry() {
    let vakint = Vakint::new().unwrap();
    let clone_before = vakint.clone();
    let base = vakint
        .topologies
        .0
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    vakint
        .topologies
        .prepare_five_loop_catalogue_for_tests()
        .unwrap();
    let catalogue = vakint.topologies.five_loop_catalogue().unwrap();
    let statistics = vakint.topologies.five_loop_catalogue_statistics().unwrap();
    assert_eq!(statistics.topology_count, 67);
    for topology in vakint
        .topologies
        .iter()
        .filter(|topology| matches!(topology, Topology::FiveLoop(_)))
    {
        let short = topology.get_integral().short_expression.as_ref().unwrap();
        assert!(
            Topologies::needs_five_loop_catalogue(crate::symbols::S.topo.call(short).as_view())
                .unwrap()
        );
    }
    assert_eq!(
        vakint
            .topologies
            .0
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>(),
        base
    );
    assert!(matches!(
        vakint.topologies.iter().last(),
        Some(Topology::Unknown(_))
    ));
    vakint
        .topologies
        .prepare_five_loop_catalogue_for_tests()
        .unwrap();
    assert_eq!(
        vakint.topologies.five_loop_catalogue_statistics(),
        Some(statistics)
    );
    assert!(
        clone_before
            .topologies
            .five_loop_catalogue_statistics()
            .is_none()
    );
    clone_before
        .topologies
        .install_five_loop_catalogue(Arc::clone(catalogue));
    assert!(Arc::ptr_eq(
        catalogue,
        clone_before.topologies.five_loop_catalogue().unwrap()
    ));
    let clone_after = vakint.clone();
    assert!(Arc::ptr_eq(
        catalogue,
        clone_after.topologies.five_loop_catalogue().unwrap()
    ));

    let mut custom = Vakint::new().unwrap();
    custom.topologies = Topologies::generate_topologies_with_contractions(
        vk_parse!("topo(prop(1,edge(1,1),k(1),msq(1),pow(1)))")
            .unwrap()
            .as_view(),
        vk_parse!("custom_scalar(msq(1),pow(1))").unwrap().as_view(),
        TopologyContractions::Custom(vec![vec![]]),
        EvaluationOrder::empty(),
    )
    .unwrap();
    let custom_before = custom.topologies.to_string();
    custom
        .topologies
        .install_five_loop_catalogue(Arc::clone(catalogue));
    assert!(custom.topologies.to_string().starts_with(&custom_before));
    let input = vk_parse!("topo(custom_scalar(m2,3))").unwrap();
    assert!(
        custom
            .topologies
            .match_topologies_to_user_input(input.as_view(), false)
            .unwrap()
            .is_some()
    );
    assert_eq!(
        custom.topologies.five_loop_catalogue_statistics(),
        Some(statistics)
    );
}

#[test]
fn unrelated_unknown_lower_loops_do_not_open_a_configured_five_loop_artifact() {
    let vakint = Vakint::new().unwrap();
    let input = vk_parse!("topo(prop(1,edge(1,2),k(1),m1,1)*prop(2,edge(2,1),k(1),m2,1))").unwrap();
    // Unequal-mass one-loop bubble: genuinely outside the base registry, not
    // a known lower topology taking the existing early-return shortcut.
    assert!(
        vakint
            .topologies
            .match_topologies_to_user_input(input.as_view(), false)
            .unwrap()
            .is_none()
    );
    assert!(matches!(
        vakint
            .topologies
            .match_topologies_to_user_input(input.as_view(), true)
            .unwrap()
            .unwrap()
            .canonical_topology,
        Topology::Unknown(_)
    ));
    assert!(!Topologies::needs_five_loop_catalogue(input.as_view()).unwrap());
    let baseline = VakintSettings {
        allow_unknown_integrals: true,
        ..Default::default()
    };
    let expected = vakint
        .to_canonical(&baseline, input.as_view(), false)
        .unwrap();
    let configured = VakintSettings {
        rustred_artifact_path: Some("/nonexistent-unrelated-five-loop-package".into()),
        ..baseline
    };
    let actual = vakint
        .to_canonical(&configured, input.as_view(), false)
        .unwrap();
    assert_eq!(actual, expected);
    assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
    assert!(vakint.rustred_application.lock().unwrap().is_none());
}

#[test]
fn unmatched_actual_five_loops_still_load_before_unknown_fallback() {
    let vakint = Vakint::new().unwrap();
    let input = vk_parse!("topo(prop(1,edge(1,1),k(1),m1,1)*prop(2,edge(1,1),k(2),m2,1)*prop(3,edge(1,1),k(3),m3,1)*prop(4,edge(1,1),k(4),m4,1)*prop(5,edge(1,1),k(5),m5,1))").unwrap();
    // The unequal masses intentionally prevent the common-mass saved match.
    assert!(
        vakint
            .topologies
            .match_topologies_to_user_input(input.as_view(), false)
            .unwrap()
            .is_none()
    );
    assert!(Topologies::needs_five_loop_catalogue(input.as_view()).unwrap());
    let settings = VakintSettings {
        allow_unknown_integrals: true,
        rustred_artifact_path: Some("/nonexistent-actual-five-loop-package".into()),
        ..Default::default()
    };
    assert!(matches!(
        vakint.to_canonical(&settings, input.as_view(), false),
        Err(VakintError::RustRedEvaluation(
            crate::rustred_evaluation::RustRedEvaluationError::Reduction { .. }
        ))
    ));
    assert!(vakint.topologies.five_loop_catalogue_statistics().is_none());
    assert!(vakint.rustred_application.lock().unwrap().is_none());
}

#[test]
fn lazy_gate_uses_exact_symbols_and_native_loop_census_errors() {
    Vakint::initialize_vakint_symbols();
    for input in [
        "topo(I5L_not_a_registered_head(1))",
        "topo(foreign::I5L_32745(1))",
        "topo(prop(1,edge(1,1),0,m2,1))",
    ] {
        assert!(
            !Topologies::needs_five_loop_catalogue(vk_parse!(input).unwrap().as_view()).unwrap()
        );
    }
    let malformed =
        vk_parse!("topo(prop(1,edge(1,1),k(1),m1,1)*prop(2,edge(1,1),k(3),m2,1))").unwrap();
    assert!(matches!(
        Topologies::needs_five_loop_catalogue(malformed.as_view()),
        Err(VakintError::InvalidIntegralFormat(_))
    ));
}
