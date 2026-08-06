use std::collections::BTreeSet;

use crate::{
    dot,
    feyngen::GenerationType,
    graph::{
        FeynmanGraph, Graph, autogen::Autogen, parse::IntoGraph,
        threshold_counterterms::ThresholdCountertermSpec,
    },
    initialisation::test_initialise,
    integrands::process::{MomentumSpaceEvaluationInput, ProcessIntegrand},
    momentum::ThreeMomentum,
    processes::{Amplitude, CrossSection, ProcessDefinition},
    settings::{
        GlobalSettings, RuntimeSettings,
        global::{GenerationSettings, ThresholdSubtractionSettings},
        runtime::kinematic::KinematicsSettings,
    },
    utils::{F, load_generic_model},
};
use spenso::algebra::complex::Complex;
use symbolica::numerical_integration::Sample;

static NO_DIRECTIVE_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

fn generation_pool() -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .stack_size(256 * 1024 * 1024)
        .build()
        .unwrap()
}

fn with_explicit_empty_directives(mut graph: Graph) -> Graph {
    graph.threshold_counterterms = Autogen::explicit(ThresholdCountertermSpec::default());
    graph
}

fn assert_no_threshold_decomposition(result: &crate::integrands::evaluation::EvaluationResult) {
    let mut event_count = 0;
    for event in result
        .event_groups
        .iter()
        .flat_map(|event_group| event_group.iter())
    {
        event_count += 1;
        assert!(
            event.additional_weights.threshold_counterterms.is_none(),
            "legacy no-directive events must not allocate detailed threshold decomposition",
        );
    }
    assert!(
        event_count > 0,
        "the runtime comparison must exercise events"
    );
}

fn assert_event_metadata_equivalent(
    left: &crate::integrands::evaluation::EvaluationResult,
    right: &crate::integrands::evaluation::EvaluationResult,
) {
    assert_eq!(left.event_groups.len(), right.event_groups.len());
    for (left_group, right_group) in left.event_groups.iter().zip(right.event_groups.iter()) {
        assert_eq!(left_group.len(), right_group.len());
        for (left_event, right_event) in left_group.iter().zip(right_group.iter()) {
            assert_eq!(
                left_event.kinematic_configuration,
                right_event.kinematic_configuration
            );
            assert_eq!(left_event.weight, right_event.weight);
            assert_eq!(
                left_event.additional_weights.weights,
                right_event.additional_weights.weights,
            );
            assert_eq!(
                left_event.cut_info.particle_pdgs,
                right_event.cut_info.particle_pdgs,
            );
            assert_eq!(left_event.cut_info.cut_id, right_event.cut_info.cut_id);
            assert_eq!(left_event.cut_info.graph_id, right_event.cut_info.graph_id);
            assert_eq!(
                left_event.cut_info.graph_group_id,
                right_event.cut_info.graph_group_id,
            );
            assert_eq!(
                left_event.cut_info.orientation_id,
                right_event.cut_info.orientation_id,
            );
            assert_eq!(
                left_event.cut_info.lmb_channel_id,
                right_event.cut_info.lmb_channel_id,
            );
            assert_eq!(
                left_event.cut_info.lmb_channel_edge_ids,
                right_event.cut_info.lmb_channel_edge_ids,
            );
        }
    }
}

fn amplitude_graph() -> Graph {
    let graph: crate::processes::AmplitudeGraph = dot!(
        digraph no_directive_amplitude {
            // A massive external scalar above the two-massless-particle threshold ensures the
            // homogeneous threshold-counterterm lane is present in this equivalence regression.
            edge [particle=scalar_0]
            node [num=1]
            e [style=invis]
            e -> A:0 [id=3 particle=scalar_2]
            B:1 -> e [id=2 particle=scalar_2]
            A -> B [id=1]
            A -> B [id=0]
        },
        "scalars"
    )
    .unwrap();
    graph.graph
}

fn build_amplitude(
    graph: Graph,
    model: &crate::model::Model,
    generation: &GenerationSettings,
    runtime: &RuntimeSettings,
) -> (Amplitude, usize) {
    let mut amplitude = Amplitude::from_graph_list("no_directive_amplitude", vec![graph]).unwrap();
    let pool = generation_pool();
    amplitude
        .preprocess(model, generation, &runtime.into(), &pool)
        .unwrap();
    let reports = amplitude
        .build_integrand(
            model,
            "no_directive_amplitude",
            &GlobalSettings {
                generation: generation.clone(),
                ..Default::default()
            },
            runtime.into(),
            &pool,
        )
        .unwrap();
    let evaluator_count = reports
        .iter()
        .map(|report| report.stats.evaluator_count)
        .sum();
    amplitude
        .integrand
        .as_mut()
        .unwrap()
        .warm_up(model)
        .unwrap();
    (amplitude, evaluator_count)
}

fn assert_amplitude_legacy_storage(amplitude: &Amplitude) -> usize {
    assert!(
        amplitude.graphs[0]
            .derived_data
            .resolved_threshold_counterterms
            .as_ref()
            .is_some_and(|resolved| resolved.legacy_equivalent)
    );
    assert!(
        !amplitude.graphs[0]
            .derived_data
            .threshold_counterterms
            .is_empty(),
        "the amplitude comparison must exercise homogeneous threshold-counterterm storage",
    );
    let ProcessIntegrand::Amplitude(integrand) = amplitude.integrand.as_ref().unwrap() else {
        unreachable!("the scalar bubble is an amplitude fixture")
    };
    let counterterm = &integrand.data.graph_terms[0].threshold_counterterm;
    assert!(counterterm.legacy_equivalent);
    assert!(counterterm.threshold_multipliers.is_none());
    assert!(counterterm.variant_evaluators.is_empty());
    assert!(counterterm.variant_helper_evaluators.is_empty());
    assert!(counterterm.variant_subspaces.is_empty());
    assert!(counterterm.variant_metadata.is_empty());
    assert!(counterterm.metadata_registry.is_none());
    integrand
        .data
        .graph_terms
        .iter()
        .map(|term| term.generic_evaluator_count())
        .sum()
}

#[test]
fn amplitude_absent_and_explicit_empty_directives_are_runtime_equivalent() {
    std::thread::Builder::new()
        .name("amplitude-no-directive-equivalence".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| {
            let _guard = NO_DIRECTIVE_TEST_LOCK.lock().unwrap();
            test_initialise().unwrap();
            let model = load_generic_model("scalars");
            let graph = amplitude_graph();
            assert!(graph.threshold_counterterms.autogenerated);
            let mut runtime = RuntimeSettings {
                kinematics: KinematicsSettings::random(&graph, 42),
                ..Default::default()
            };
            runtime.general.generate_events = true;
            runtime.general.store_additional_weights_in_event = true;
            let generation = GenerationSettings {
                threshold_subtraction: ThresholdSubtractionSettings {
                    enable_thresholds: true,
                    assume_positive_external_energies: false,
                    check_esurface_at_generation: false,
                    ..Default::default()
                },
                ..Default::default()
            };

            let (mut absent, absent_report_count) =
                build_amplitude(graph.clone(), &model, &generation, &runtime);
            let explicit_graph = with_explicit_empty_directives(graph);
            assert!(!explicit_graph.threshold_counterterms.autogenerated);
            assert!(explicit_graph.threshold_counterterms.cuts.is_empty());
            let (mut explicit, explicit_report_count) =
                build_amplitude(explicit_graph, &model, &generation, &runtime);
            assert_eq!(absent_report_count, explicit_report_count);
            assert_eq!(
                assert_amplitude_legacy_storage(&absent),
                assert_amplitude_legacy_storage(&explicit),
            );

            let input = MomentumSpaceEvaluationInput {
                loop_momenta: vec![ThreeMomentum::new(F(0.35), F(-0.2), F(0.45))],
                integrator_weight: F(1.0),
                graph_id: Some(0),
                group_id: None,
                orientation: None,
                channel_id: None,
            };
            let absent_result = absent
                .integrand
                .as_mut()
                .unwrap()
                .evaluate_momentum_configuration(&model, &input, false)
                .unwrap();
            let explicit_result = explicit
                .integrand
                .as_mut()
                .unwrap()
                .evaluate_momentum_configuration(&model, &input, false)
                .unwrap();
            assert_eq!(
                absent_result.integrand_result,
                explicit_result.integrand_result
            );
            assert_no_threshold_decomposition(&absent_result);
            assert_no_threshold_decomposition(&explicit_result);
            assert_event_metadata_equivalent(&absent_result, &explicit_result);
        })
        .unwrap()
        .join()
        .unwrap();
}

fn cross_section_runtime() -> RuntimeSettings {
    let mut runtime: RuntimeSettings = toml::from_str(
        r#"
[kinematics]
e_cm = 1000.0

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
  [500.0, 0.0, 0.0, 500.0],
  [500.0, 0.0, 0.0, -500.0],
]
helicities = ["summed_averaged", "summed_averaged"]
"#,
    )
    .unwrap();
    runtime.general.generate_events = true;
    runtime.general.store_additional_weights_in_event = true;
    runtime
}

fn cross_section_generation() -> GenerationSettings {
    let mut generation = GenerationSettings::default();
    // This regression targets runtime-path equivalence; UV subtraction and integrated threshold
    // pieces have independent coverage and only add expensive evaluators here.
    generation.uv.softct = false;
    generation.uv.subtract_uv = false;
    generation.threshold_subtraction.disable_integrated_ct = true;
    generation
}

fn cross_section_graph() -> Graph {
    dot!(
        digraph no_directive_cross_section {
            graph [num=1, overall_factor=1]
            edge [particle=scalar_0]
            node [num=1]
            exte0 [style=invis, is_cut=0]
            exte1 [style=invis, is_cut=1]
            v1 -> exte1 [id=0]
            v1 -> exte0 [id=1]
            exte0 -> v2 [id=2]
            exte1 -> v2 [id=3]
            v2 -> a [id=4]
            a -> v1 [id=5]
            v2 -> b [id=6]
            b -> v1 [id=7]
        },
        "scalars"
    )
    .unwrap()
}

fn build_cross_section(
    graph: Graph,
    model: &crate::model::Model,
    generation: &GenerationSettings,
    runtime: &RuntimeSettings,
) -> (CrossSection, usize) {
    let definition = ProcessDefinition::from_graph_list(
        std::slice::from_ref(&graph),
        GenerationType::CrossSection,
        model,
    )
    .unwrap();
    let mut cross_section =
        CrossSection::from_graph_list("no_directive_cross_section".to_string(), vec![graph], model)
            .unwrap();
    let pool = generation_pool();
    cross_section
        .preprocess(model, &definition, generation, runtime.into(), &pool)
        .unwrap();
    let reports = cross_section
        .build_integrand(
            model,
            "no_directive_cross_section",
            &GlobalSettings {
                generation: generation.clone(),
                ..Default::default()
            },
            runtime.into(),
            &pool,
        )
        .unwrap();
    let evaluator_count = reports
        .iter()
        .map(|report| report.stats.evaluator_count)
        .sum();
    cross_section
        .integrand
        .as_mut()
        .unwrap()
        .warm_up(model)
        .unwrap();
    (cross_section, evaluator_count)
}

fn assert_cross_section_legacy_storage(cross_section: &CrossSection) -> usize {
    assert!(
        cross_section.supergraphs[0]
            .derived_data
            .resolved_threshold_counterterms
            .as_ref()
            .is_some_and(|resolved| resolved.legacy_equivalent)
    );
    let ProcessIntegrand::CrossSection(integrand) = cross_section.integrand.as_ref().unwrap()
    else {
        unreachable!("the forward scalar graph is a cross-section fixture")
    };
    let counterterm = &integrand.data.graph_terms[0].counterterm;
    assert!(counterterm.variant_subspaces.is_none());
    assert!(counterterm.metadata_registry.is_none());
    assert!(
        counterterm
            .evaluators
            .iter()
            .all(|evaluators| evaluators.threshold_multipliers.is_none())
    );
    counterterm
        .evaluators
        .iter()
        .map(|evaluators| evaluators.generic_compileable_evaluator_count())
        .sum()
}

fn evaluate_cross_section(
    cross_section: &mut CrossSection,
    model: &crate::model::Model,
) -> crate::integrands::evaluation::EvaluationResult {
    let n_dim = cross_section.supergraphs[0].graph.get_loop_number() * 3;
    let sample = Sample::Continuous(
        F(1.0),
        (0..n_dim)
            .map(|axis| F(0.12 + ((axis * 7 + 5) % 17) as f64 * 0.043))
            .collect(),
    );
    cross_section
        .integrand
        .as_mut()
        .unwrap()
        .evaluate_samples_raw(
            model,
            &[sample],
            0,
            false,
            false,
            Complex::new(F(0.0), F(0.0)),
        )
        .unwrap()
        .samples
        .remove(0)
}

#[test]
fn cross_section_absent_and_explicit_empty_directives_are_runtime_equivalent() {
    std::thread::Builder::new()
        .name("cross-section-no-directive-equivalence".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| {
            let _guard = NO_DIRECTIVE_TEST_LOCK.lock().unwrap();
            test_initialise().unwrap();
            let model = load_generic_model("scalars");
            let graph = cross_section_graph();
            assert!(graph.threshold_counterterms.autogenerated);
            let runtime = cross_section_runtime();
            let generation = cross_section_generation();
            assert!(
                generation.force_cuts.is_empty(),
                "the numerical comparison must retain every process-valid cut",
            );

            let (mut absent, absent_report_count) =
                build_cross_section(graph.clone(), &model, &generation, &runtime);
            let explicit_graph = with_explicit_empty_directives(graph);
            assert!(!explicit_graph.threshold_counterterms.autogenerated);
            assert!(explicit_graph.threshold_counterterms.cuts.is_empty());
            let (mut explicit, explicit_report_count) =
                build_cross_section(explicit_graph, &model, &generation, &runtime);
            assert!(!absent.supergraphs[0].cuts.is_empty());
            assert_eq!(
                absent.supergraphs[0].cuts.len(),
                explicit.supergraphs[0].cuts.len(),
            );
            for cross_section in [&absent, &explicit] {
                let graph = &cross_section.supergraphs[0];
                let grouped_cut_ids = graph
                    .derived_data
                    .cut_group_data
                    .cut_groups
                    .iter()
                    .flat_map(|group| group.cuts.iter().map(|cut_id| cut_id.0))
                    .collect::<BTreeSet<_>>();
                assert_eq!(
                    grouped_cut_ids,
                    (0..graph.cuts.len()).collect(),
                    "every discovered process-valid cut must remain in the numerical LU evaluation",
                );
            }
            assert_eq!(absent_report_count, explicit_report_count);
            assert_eq!(
                assert_cross_section_legacy_storage(&absent),
                assert_cross_section_legacy_storage(&explicit),
            );

            let absent_result = evaluate_cross_section(&mut absent, &model);
            let explicit_result = evaluate_cross_section(&mut explicit, &model);
            assert_eq!(
                absent_result.integrand_result,
                explicit_result.integrand_result
            );
            assert_no_threshold_decomposition(&absent_result);
            assert_no_threshold_decomposition(&explicit_result);
            assert_event_metadata_equivalent(&absent_result, &explicit_result);
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn cross_section_with_threshold_generation_disabled_builds_empty_legacy_lu_state() {
    std::thread::Builder::new()
        .name("cross-section-disabled-threshold-constructor".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| {
            let _guard = NO_DIRECTIVE_TEST_LOCK.lock().unwrap();
            test_initialise().unwrap();
            let model = load_generic_model("scalars");
            let mut generation = cross_section_generation();
            generation.threshold_subtraction.enable_thresholds = false;
            let mut runtime = cross_section_runtime();
            runtime.subtraction.disable_threshold_subtraction = true;
            let (cross_section, _) =
                build_cross_section(cross_section_graph(), &model, &generation, &runtime);

            let graph = &cross_section.supergraphs[0];
            assert!(graph.derived_data.resolved_threshold_counterterms.is_none());
            assert!(graph.derived_data.threshold_counterterms.is_empty());
            let ProcessIntegrand::CrossSection(integrand) =
                cross_section.integrand.as_ref().unwrap()
            else {
                unreachable!("the forward scalar graph is a cross-section fixture")
            };
            let counterterm = &integrand.data.graph_terms[0].counterterm;
            assert!(counterterm.evaluators.is_empty());
            assert!(counterterm.thresholds.is_empty());
            assert!(counterterm.variant_subspaces.is_none());
            assert!(counterterm.metadata_registry.is_none());
        })
        .unwrap()
        .join()
        .unwrap();
}
