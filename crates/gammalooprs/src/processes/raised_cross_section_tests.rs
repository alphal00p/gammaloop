use crate::integrands::process::EvaluationTarget;
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::Cursor,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use crate::{
    GammaLoopContextContainer,
    cff::CutCFFIndex,
    feyngen::GenerationType,
    graph::{
        Graph,
        autogen::Autogen,
        feynman_graph::FeynmanGraph,
        parse::IntoGraph,
        threshold_counterterms::{
            THRESHOLD_COUNTERTERM_SCHEMA_VERSION, ThresholdCountertermCut,
            ThresholdCountertermMultiplier, ThresholdCountertermSpec,
            ThresholdCountertermThreshold, ThresholdCountertermVariant,
        },
    },
    initialisation::test_initialise,
    integrands::{evaluation::EvaluationResult, process::ProcessIntegrand},
    model::Model,
    observables::ThresholdCountertermComponentOccurrence,
    processes::{
        CrossSection, CutGroupId, ProcessDefinition,
        threshold_counterterms::{
            ThresholdCountertermComponentKind, ThresholdCountertermMetadataRegistry,
            ThresholdCountertermSide, ThresholdCountertermVariantId,
        },
    },
    settings::{GlobalSettings, RuntimeSettings, global::GenerationSettings},
    utils::F,
};
use linnet::half_edge::involution::EdgeIndex;
use spenso::algebra::complex::Complex;
use symbolica::{numerical_integration::Sample, state::State};

const TRIPLE_DOTTED_BUBBLE: &str =
    include_str!("../../../../tests/resources/graphs/ir_safe_thresholds/triple_dotted_bubble.dot");
const SCALARS_2P_3P_MODEL: &str =
    include_str!("../../../../assets/models/json/scalars/scalars_2p_3p.json");
const RAISED_VARIANT_NAME: &str = "raised_order_two";

struct TemporaryDirectory(PathBuf);

impl TemporaryDirectory {
    fn new(prefix: &str) -> Self {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gammalooprs-{prefix}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&path).unwrap();
        Self(path)
    }
}

impl Drop for TemporaryDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn runtime_settings() -> RuntimeSettings {
    let mut runtime: RuntimeSettings = toml::from_str(
        r#"
[kinematics]
e_cm = 5.0

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [[5.0, 0.0, 0.0, 0.0]]
helicities = ["summed_averaged"]
"#,
    )
    .unwrap();
    runtime.general.generate_events = true;
    runtime.general.store_additional_weights_in_event = true;
    runtime
}

fn generation_settings() -> GenerationSettings {
    let mut settings = GenerationSettings::default();
    settings.uv.softct = false;
    settings.uv.subtract_uv = false;
    settings.threshold_subtraction.enable_thresholds = true;
    settings.threshold_subtraction.check_esurface_at_generation = false;
    settings.threshold_subtraction.skip_thresholds_that_are_cuts = false;
    settings
        .threshold_subtraction
        .assume_positive_external_energies = false;
    settings
}

fn generation_pool() -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .stack_size(256 * 1024 * 1024)
        .build()
        .unwrap()
}

fn preprocess_graph(
    graph: Graph,
    model: &Model,
    generation: &GenerationSettings,
    runtime: &RuntimeSettings,
) -> CrossSection {
    let definition = ProcessDefinition::from_graph_list(
        std::slice::from_ref(&graph),
        GenerationType::CrossSection,
        model,
    )
    .unwrap();
    let mut cross_section =
        CrossSection::from_graph_list("raised_cross_section".to_string(), vec![graph], model)
            .unwrap();
    cross_section
        .preprocess(
            model,
            &definition,
            generation,
            runtime.into(),
            &generation_pool(),
        )
        .unwrap();
    cross_section
}

fn assert_every_discovered_cut_is_grouped(cross_section: &CrossSection) {
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
        (0..graph.cuts.len()).collect::<BTreeSet<_>>(),
        "every discovered process-valid cut must remain in the numerical LU evaluation",
    );
}

/// Materialize a constant-multiplier directive for every physical association belonging to the
/// one raised threshold on each side of the middle cut group. Deriving these identities from the
/// legacy preprocessing keeps this synthetic regression independent of incidental edge ordering.
fn target_directives_from_legacy(cross_section: &CrossSection) -> ThresholdCountertermSpec {
    let graph = &cross_section.supergraphs[0];
    let resolved = graph
        .derived_data
        .resolved_threshold_counterterms
        .as_ref()
        .unwrap();
    assert!(resolved.legacy_equivalent);
    let target_groups = resolved
        .cross_section_cut_groups
        .iter()
        .filter(|group| !group.left.is_empty() && !group.right.is_empty())
        .collect::<Vec<_>>();
    assert_eq!(
        target_groups.len(),
        1,
        "the synthetic serial topology must have one cut group with CTs on both sides",
    );
    let target_group = target_groups[0];
    assert_eq!(target_group.left.len(), 1);
    assert_eq!(target_group.right.len(), 1);

    let mut cuts = BTreeMap::<Vec<EdgeIndex>, BTreeMap<Vec<EdgeIndex>, Vec<EdgeIndex>>>::new();
    for variant_id in target_group.left.iter().chain(&target_group.right) {
        let variant = &resolved.variants[*variant_id];
        assert_eq!(variant.raised_esurface_group.max_occurence, 2);
        assert_eq!(variant.subspace_loop_count, 1);
        assert!(
            variant
                .associations
                .iter()
                .all(|association| association.eligible)
        );
        for association in &variant.associations {
            cuts.entry(association.cut_edges.clone())
                .or_default()
                .insert(
                    association.threshold_edges.clone(),
                    association
                        .subspace
                        .get_lmb(graph.derived_data.lmbs.as_ref().unwrap())
                        .loop_edges
                        .raw
                        .clone(),
                );
        }
    }

    let counterterm = |parent_lmb| ThresholdCountertermVariant {
        name: Some(RAISED_VARIANT_NAME.to_string()),
        subspace: None,
        parent_lmb: Some(parent_lmb),
        group_id: None,
        disable: false,
        multiplier: Some(ThresholdCountertermMultiplier {
            expression: "2".to_string(),
            symmetrize: false,
            opaque_derivatives: true,
        }),
    };
    ThresholdCountertermSpec {
        schema_version: THRESHOLD_COUNTERTERM_SCHEMA_VERSION,
        cuts: cuts
            .into_iter()
            .map(|(edges, thresholds)| ThresholdCountertermCut {
                edges,
                thresholds: thresholds
                    .into_iter()
                    .map(|(edges, parent_lmb)| ThresholdCountertermThreshold {
                        edges,
                        counterterms: vec![counterterm(parent_lmb)],
                    })
                    .collect(),
            })
            .collect(),
    }
}

fn target_variants(
    cross_section: &CrossSection,
) -> (
    CutGroupId,
    ThresholdCountertermVariantId,
    ThresholdCountertermVariantId,
) {
    let resolved = cross_section.supergraphs[0]
        .derived_data
        .resolved_threshold_counterterms
        .as_ref()
        .unwrap();
    let targets =
        resolved
            .cross_section_cut_groups
            .iter_enumerated()
            .filter_map(|(cut_group_id, group)| {
                let left =
                    group.left.iter().copied().find(|variant_id| {
                        resolved.variants[*variant_id].name == RAISED_VARIANT_NAME
                    });
                let right =
                    group.right.iter().copied().find(|variant_id| {
                        resolved.variants[*variant_id].name == RAISED_VARIANT_NAME
                    });
                left.zip(right)
                    .map(|(left, right)| (cut_group_id, left, right))
            })
            .collect::<Vec<_>>();
    assert_eq!(targets.len(), 1);
    targets[0]
}

fn single_indices(side: ThresholdCountertermSide) -> BTreeSet<CutCFFIndex> {
    [1, 2]
        .into_iter()
        .flat_map(|threshold_order| {
            [1, 2].into_iter().map(move |lu_cut_order| CutCFFIndex {
                left_threshold_order: (side == ThresholdCountertermSide::Left)
                    .then_some(threshold_order),
                right_threshold_order: (side == ThresholdCountertermSide::Right)
                    .then_some(threshold_order),
                lu_cut_order: Some(lu_cut_order),
            })
        })
        .collect()
}

fn iterated_indices() -> BTreeSet<CutCFFIndex> {
    [1, 2]
        .into_iter()
        .flat_map(|left_threshold_order| {
            [1, 2].into_iter().flat_map(move |right_threshold_order| {
                [1, 2].into_iter().map(move |lu_cut_order| CutCFFIndex {
                    left_threshold_order: Some(left_threshold_order),
                    right_threshold_order: Some(right_threshold_order),
                    lu_cut_order: Some(lu_cut_order),
                })
            })
        })
        .collect()
}

fn component_ids_for_target(
    metadata: &ThresholdCountertermMetadataRegistry,
    cut_group_id: CutGroupId,
    left_variant_id: ThresholdCountertermVariantId,
    right_variant_id: ThresholdCountertermVariantId,
) -> BTreeSet<usize> {
    metadata
        .components
        .iter()
        .filter(|component| component.cut_group_id == Some(cut_group_id.0))
        .filter(|component| match component.variant_ids.as_slice() {
            [variant_id] => *variant_id == left_variant_id.0 || *variant_id == right_variant_id.0,
            [left, right] => *left == left_variant_id.0 && *right == right_variant_id.0,
            _ => false,
        })
        .map(|component| component.component_id)
        .collect()
}

fn assert_runtime_decomposition(
    results: &[EvaluationResult],
    metadata: &ThresholdCountertermMetadataRegistry,
    target_component_ids: &BTreeSet<usize>,
) -> BTreeMap<usize, BTreeSet<CutCFFIndex>> {
    let mut event_count = 0;
    let mut observed_target_components = BTreeSet::new();
    let mut observed_target_occurrences = BTreeMap::<usize, BTreeSet<CutCFFIndex>>::new();
    for result in results {
        assert!(!result.evaluation_metadata.is_nan);
        assert!(result.integrand_result.re.0.is_finite());
        assert!(result.integrand_result.im.0.is_finite());
        for event in result
            .event_groups
            .iter()
            .flat_map(|event_group| event_group.iter())
        {
            event_count += 1;
            let decomposition = event
                .additional_weights
                .threshold_counterterms
                .as_ref()
                .expect("the generalized path must record threshold components");
            assert_eq!(event.weight, decomposition.total());
            for component in &decomposition.components {
                let component_metadata = &metadata.components[component.component_id];
                assert_eq!(component_metadata.component_id, component.component_id);
                assert_eq!(
                    component_metadata.variant_ids.len(),
                    component.multiplier_values.len(),
                );
                assert!(component.weighted.re.0.is_finite());
                assert!(component.weighted.im.0.is_finite());
                if !target_component_ids.contains(&component.component_id) {
                    continue;
                }

                observed_target_components.insert(component.component_id);
                assert!(!component.evaluation_skipped);
                let bare = component.bare.as_ref().unwrap();
                assert!(bare.re.0.is_finite() && bare.im.0.is_finite());
                assert!(
                    component
                        .multiplier_values
                        .iter()
                        .all(|value| *value == F(2.0))
                );
                assert_eq!(
                    component.effective_multiplier,
                    F(if component_metadata.kind.variant_count() == 1 {
                        2.0
                    } else {
                        4.0
                    }),
                );
                let ThresholdCountertermComponentOccurrence::LocalUnitarity {
                    left_threshold_order,
                    right_threshold_order,
                    lu_cut_order,
                    ..
                } = component.occurrence
                else {
                    panic!("cross-section threshold components must have LU provenance")
                };
                observed_target_occurrences
                    .entry(component.component_id)
                    .or_default()
                    .insert(CutCFFIndex {
                        left_threshold_order,
                        right_threshold_order,
                        lu_cut_order,
                    });
            }
        }
    }
    assert!(event_count > 0, "the fixed samples must generate events");
    assert_eq!(&observed_target_components, target_component_ids);
    observed_target_occurrences
}

#[test]
fn raised_cross_section_rejects_different_groups_for_one_merged_variant() {
    std::thread::Builder::new()
        .name("raised-cross-section-group-conflict".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| {
            test_initialise().unwrap();
            let model = Model::from_str(SCALARS_2P_3P_MODEL.to_string(), "json").unwrap();
            let mut graph: Graph = TRIPLE_DOTTED_BUBBLE.into_graph(&model).unwrap();
            let generation = generation_settings();
            let runtime = runtime_settings();
            let baseline = preprocess_graph(graph.clone(), &model, &generation, &runtime);
            let resolved = baseline.supergraphs[0]
                .derived_data
                .resolved_threshold_counterterms
                .as_ref()
                .unwrap();
            let group = resolved
                .cross_section_cut_groups
                .iter()
                .find(|group| !group.left.is_empty() && !group.right.is_empty())
                .expect("the raised fixture must contain a middle cut group");
            let variant = &resolved.variants[group.left[0]];
            assert!(variant.associations.len() > 1);
            assert!(variant.raised_esurface_group.max_occurence > 1);
            let mut directives = target_directives_from_legacy(&baseline);
            for (index, association) in variant.associations.iter().enumerate() {
                let cut = directives
                    .cuts
                    .iter_mut()
                    .find(|cut| cut.edges == association.cut_edges)
                    .unwrap();
                let threshold = cut
                    .thresholds
                    .iter_mut()
                    .find(|threshold| threshold.edges == association.threshold_edges)
                    .unwrap();
                threshold.counterterms[0].group_id = Some(usize::from(index == 0));
            }
            graph.threshold_counterterms = Autogen::explicit(directives);
            let definition = ProcessDefinition::from_graph_list(
                std::slice::from_ref(&graph),
                GenerationType::CrossSection,
                &model,
            )
            .unwrap();
            let mut cross_section = CrossSection::from_graph_list(
                "raised_group_conflict".to_string(),
                vec![graph],
                &model,
            )
            .unwrap();
            let error = cross_section
                .preprocess(
                    &model,
                    &definition,
                    &generation,
                    (&runtime).into(),
                    &generation_pool(),
                )
                .expect_err("a merged raised variant cannot use different solve groups");
            let diagnostic = format!("{error:#}");
            for expected in [
                "one merged raised residue",
                "group_id=Some(0)",
                "group_id=Some(1)",
                "parent=",
                "signed_cycles=",
                "cut_edges=",
            ] {
                assert!(
                    diagnostic.contains(expected),
                    "missing {expected}: {diagnostic}"
                );
            }
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn generalized_raised_cross_section_covers_derivative_components_and_roundtrips() {
    std::thread::Builder::new()
        .name("raised-cross-section-directives".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| {
            test_initialise().unwrap();
            let model = Model::from_str(SCALARS_2P_3P_MODEL.to_string(), "json").unwrap();
            let graph: Graph = TRIPLE_DOTTED_BUBBLE.into_graph(&model).unwrap();
            let generation = generation_settings();
            assert!(
                generation.force_cuts.is_empty(),
                "the raised-order regression must retain every process-valid cut",
            );
            let runtime = runtime_settings();

            let baseline = preprocess_graph(graph.clone(), &model, &generation, &runtime);
            assert_every_discovered_cut_is_grouped(&baseline);
            assert_eq!(baseline.supergraphs[0].cuts.len(), 6);
            assert!(
                baseline.supergraphs[0]
                    .derived_data
                    .cut_group_data
                    .cut_groups
                    .iter()
                    .all(|group| group.related_esurface_group.max_occurence == 2)
            );
            let directives = target_directives_from_legacy(&baseline);

            let mut graph = graph;
            graph.threshold_counterterms = Autogen::explicit(directives);
            let mut cross_section = preprocess_graph(graph, &model, &generation, &runtime);
            assert_every_discovered_cut_is_grouped(&cross_section);
            let resolved = cross_section.supergraphs[0]
                .derived_data
                .resolved_threshold_counterterms
                .as_ref()
                .unwrap();
            assert!(!resolved.legacy_equivalent);
            let (cut_group_id, left_variant_id, right_variant_id) = target_variants(&cross_section);
            let left_variant = &resolved.variants[left_variant_id];
            let right_variant = &resolved.variants[right_variant_id];
            assert_eq!(left_variant.side, ThresholdCountertermSide::Left);
            assert_eq!(right_variant.side, ThresholdCountertermSide::Right);
            for variant in [left_variant, right_variant] {
                assert_eq!(variant.subspace_loop_count, 1);
                assert_eq!(variant.raised_esurface_group.max_occurence, 2);
                assert_eq!(variant.multiplier.as_ref().unwrap().expression, "2");
            }

            let generated = &cross_section.supergraphs[0]
                .derived_data
                .threshold_counterterms[cut_group_id];
            assert_eq!(
                generated
                    .left_variant_ids
                    .iter()
                    .copied()
                    .collect::<Vec<_>>(),
                [left_variant_id],
            );
            assert_eq!(
                generated
                    .right_variant_ids
                    .iter()
                    .copied()
                    .collect::<Vec<_>>(),
                [right_variant_id],
            );
            assert_eq!(generated.left_thresholds.first().unwrap().max_occurence, 2);
            assert_eq!(generated.right_thresholds.first().unwrap().max_occurence, 2);
            assert_eq!(
                generated
                    .left_atoms
                    .first()
                    .unwrap()
                    .integrands
                    .iter()
                    .map(|(index, _)| *index)
                    .collect::<BTreeSet<_>>(),
                single_indices(ThresholdCountertermSide::Left),
            );
            assert_eq!(
                generated
                    .right_atoms
                    .first()
                    .unwrap()
                    .integrands
                    .iter()
                    .map(|(index, _)| *index)
                    .collect::<BTreeSet<_>>(),
                single_indices(ThresholdCountertermSide::Right),
            );
            assert_eq!(generated.iterated.iter().count(), 1);
            assert_eq!(
                generated
                    .iterated
                    .iter()
                    .next()
                    .unwrap()
                    .integrands
                    .iter()
                    .map(|(index, _)| *index)
                    .collect::<BTreeSet<_>>(),
                iterated_indices(),
            );

            cross_section
                .build_integrand(
                    &model,
                    &ProcessDefinition::from_graph_list(
                        std::slice::from_ref(&cross_section.supergraphs[0].graph),
                        GenerationType::CrossSection,
                        &model,
                    )
                    .unwrap(),
                    &GlobalSettings {
                        generation: generation.clone(),
                        ..Default::default()
                    },
                    (&runtime).into(),
                    &generation_pool(),
                )
                .unwrap();
            let n_dim = cross_section.supergraphs[0].graph.get_loop_number() * 3;
            let samples = (0..2)
                .map(|sample_index| {
                    Sample::Continuous(
                        F(1.0),
                        (0..n_dim)
                            .map(|axis| {
                                F(0.12 + ((axis * 7 + sample_index * 5) % 17) as f64 * 0.043)
                            })
                            .collect(),
                    )
                })
                .collect::<Vec<_>>();

            let (metadata, before_save) = {
                let integrand = cross_section.integrand.as_mut().unwrap();
                let ProcessIntegrand::CrossSection(cross_section_integrand) = &*integrand else {
                    unreachable!("the synthetic topology is a cross section")
                };
                let term = &cross_section_integrand.data.graph_terms[0];
                let evaluators = &term.counterterm.evaluators[cut_group_id];
                let multipliers = evaluators.threshold_multipliers.as_ref().unwrap();
                assert_eq!(multipliers.evaluators().len(), 1);
                assert_eq!(multipliers.left_variants().len(), 1);
                assert_eq!(multipliers.right_variants().len(), 1);
                assert_eq!(
                    evaluators
                        .left_thresholds_evaluator
                        .first()
                        .unwrap()
                        .keys()
                        .copied()
                        .collect::<BTreeSet<_>>(),
                    single_indices(ThresholdCountertermSide::Left),
                );
                assert_eq!(
                    evaluators
                        .right_thresholds_evaluator
                        .first()
                        .unwrap()
                        .keys()
                        .copied()
                        .collect::<BTreeSet<_>>(),
                    single_indices(ThresholdCountertermSide::Right),
                );
                assert_eq!(
                    evaluators
                        .iterated_evaluator
                        .iter()
                        .next()
                        .unwrap()
                        .keys()
                        .copied()
                        .collect::<BTreeSet<_>>(),
                    iterated_indices(),
                );
                let metadata = term
                    .threshold_counterterm_metadata()
                    .expect("the generalized path must allocate metadata")
                    .clone();
                integrand.warm_up(&model).unwrap();
                let results = integrand
                    .evaluate_samples_raw(
                        EvaluationTarget::Physical(&model),
                        &samples,
                        0,
                        false,
                        false,
                        Complex::new(F(0.0), F(0.0)),
                    )
                    .unwrap()
                    .samples;
                (metadata, results)
            };
            let target_component_ids = component_ids_for_target(
                &metadata,
                cut_group_id,
                left_variant_id,
                right_variant_id,
            );
            assert_eq!(target_component_ids.len(), 8);
            assert_eq!(
                target_component_ids
                    .iter()
                    .map(|component_id| metadata.components[*component_id].kind)
                    .collect::<BTreeSet<_>>(),
                BTreeSet::from([
                    ThresholdCountertermComponentKind::Local,
                    ThresholdCountertermComponentKind::Integrated,
                    ThresholdCountertermComponentKind::LocalLocal,
                    ThresholdCountertermComponentKind::LocalIntegrated,
                    ThresholdCountertermComponentKind::IntegratedLocal,
                    ThresholdCountertermComponentKind::IntegratedIntegrated,
                ]),
            );
            let before_occurrences =
                assert_runtime_decomposition(&before_save, &metadata, &target_component_ids);
            // Sampling is a proposal for the complete raised-cut/CT sum. Its
            // Jacobian must multiply the result and every event only after
            // residue differentiation, independently of the chosen profile.
            {
                use crate::{
                    integrands::process::{MomentumSpaceEvaluationInput, SamplingChannelId},
                    momentum::ThreeMomentum,
                    settings::runtime::{
                        MultiChannelingSettings, SamplingChannelDefinition, SamplingRadialProfile,
                        SamplingSettings, StabilityLevelSetting,
                    },
                };
                let mut matched = cross_section.integrand.as_ref().unwrap().clone();
                let ProcessIntegrand::CrossSection(prepared) = &matched else {
                    unreachable!()
                };
                let term = &prepared.data.graph_terms[0];
                let graph_name = term.graph.name.clone();
                let parent_lmb = term
                    .graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .map(|edge| edge.0)
                    .collect::<Vec<_>>();
                let cut = term.cut_esurface.iter().next().unwrap();
                let edges = cut
                    .energies
                    .iter()
                    .map(|edge| edge.0.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                for scale in [None, Some(1.8)] {
                    let mut parameterization =
                        crate::settings::runtime::ParameterizationSettings::default();
                    parameterization.sampling_channels.default_channel_selection =
                        vec!["raised_cut".to_owned()];
                    parameterization
                        .sampling_channels
                        .channel_definitions
                        .entry(graph_name.clone())
                        .or_default()
                        .insert(
                            "raised_cut".to_owned(),
                            SamplingChannelDefinition {
                                around: format!("phase_space(cut({edges}))"),
                                parent_lmb: parent_lmb.clone(),
                                subspace_lmb: parent_lmb.clone(),
                                radial_profile: Some(SamplingRadialProfile {
                                    scale,
                                    ..Default::default()
                                }),
                                ..Default::default()
                            },
                        );
                    let settings = matched.get_mut_settings();
                    settings.sampling =
                        SamplingSettings::MultiChanneling(MultiChannelingSettings {
                            parameterization_settings: parameterization,
                            ..Default::default()
                        });
                    settings.stability.rotation_axis.clear();
                    settings.stability.levels = vec![StabilityLevelSetting::default_double()];
                    matched.warm_up(&model).unwrap();
                    let ProcessIntegrand::CrossSection(prepared) = &matched else {
                        unreachable!()
                    };
                    let bridge = prepared.data.graph_terms[0]
                        .multi_channeling_setup
                        .sampling_bridge::<f64>()
                        .unwrap()
                        .clone();
                    for sample in &samples {
                        let Sample::Continuous(_, coordinates) = sample else {
                            unreachable!()
                        };
                        let point = bridge
                            .forward(
                                SamplingChannelId(0),
                                &coordinates.iter().map(|x| x.0).collect::<Vec<_>>(),
                            )
                            .unwrap();
                        assert!(
                            point
                                .map
                                .diagnostics
                                .iter()
                                .any(|message| message.contains("max_occurrence=2")),
                            "{:?}",
                            point.map.diagnostics
                        );
                        let raw = matched
                            .evaluate_momentum_configuration(
                                &model,
                                &MomentumSpaceEvaluationInput {
                                    loop_momenta: point
                                        .raw_coordinates
                                        .chunks_exact(3)
                                        .map(|k| ThreeMomentum::new(F(k[0]), F(k[1]), F(k[2])))
                                        .collect(),
                                    integrator_weight: F(1.0),
                                    graph_id: Some(0),
                                    group_id: None,
                                    orientation: None,
                                    channel_id: None,
                                },
                                false,
                            )
                            .unwrap();
                        let mapped = matched
                            .evaluate_samples_raw(
                                EvaluationTarget::Physical(&model),
                                std::slice::from_ref(sample),
                                0,
                                false,
                                false,
                                Complex::new(F(0.0), F(0.0)),
                            )
                            .unwrap()
                            .samples
                            .remove(0);
                        let factor = point.selected_factor().unwrap();
                        assert_eq!(mapped.parameterization_jacobian, Some(F(1.0)));
                        for (left, right) in [
                            (
                                mapped.integrand_result.re.0,
                                raw.integrand_result.re.0 * factor,
                            ),
                            (
                                mapped.integrand_result.im.0,
                                raw.integrand_result.im.0 * factor,
                            ),
                        ] {
                            assert!(
                                (left - right).abs() <= 1.0e-8 * left.abs().max(right.abs()),
                                "{left} != {right}"
                            );
                        }
                        let raw_events = raw
                            .event_groups
                            .iter()
                            .flat_map(|group| group.iter())
                            .collect::<Vec<_>>();
                        let mapped_events = mapped
                            .event_groups
                            .iter()
                            .flat_map(|group| group.iter())
                            .collect::<Vec<_>>();
                        assert!(!raw_events.is_empty());
                        assert_eq!(raw_events.len(), mapped_events.len());
                        for (raw, mapped) in raw_events.iter().zip(mapped_events) {
                            for (left, right) in [
                                (mapped.weight.re.0, raw.weight.re.0 * factor),
                                (mapped.weight.im.0, raw.weight.im.0 * factor),
                            ] {
                                assert!(
                                    (left - right).abs() <= 1.0e-8 * left.abs().max(right.abs()),
                                    "event {left} != {right}"
                                );
                            }
                        }
                    }
                }
            }

            for component_id in &target_component_ids {
                let component = &metadata.components[*component_id];
                let expected = match component.kind {
                    ThresholdCountertermComponentKind::Local
                    | ThresholdCountertermComponentKind::Integrated => {
                        if component.variant_ids == [left_variant_id.0] {
                            single_indices(ThresholdCountertermSide::Left)
                        } else {
                            assert_eq!(component.variant_ids, [right_variant_id.0]);
                            single_indices(ThresholdCountertermSide::Right)
                        }
                    }
                    ThresholdCountertermComponentKind::LocalLocal
                    | ThresholdCountertermComponentKind::LocalIntegrated
                    | ThresholdCountertermComponentKind::IntegratedLocal
                    | ThresholdCountertermComponentKind::IntegratedIntegrated => iterated_indices(),
                };
                assert_eq!(before_occurrences[component_id], expected);
            }

            let save_root = TemporaryDirectory::new("raised-cross-section");
            cross_section.save(&save_root.0, true).unwrap();
            let mut state_bytes = Vec::new();
            State::export(&mut state_bytes).unwrap();
            let state_map = State::import(&mut Cursor::new(state_bytes), None).unwrap();
            let context = GammaLoopContextContainer {
                model: &model,
                state_map: &state_map,
            };
            let mut loaded =
                CrossSection::load(save_root.0.join("raised_cross_section"), context).unwrap();
            let loaded_integrand = loaded.integrand.as_mut().unwrap();
            let ProcessIntegrand::CrossSection(loaded_cross_section) = &*loaded_integrand else {
                unreachable!("the saved synthetic topology is a cross section")
            };
            assert_eq!(
                loaded_cross_section.data.graph_terms[0]
                    .threshold_counterterm_metadata()
                    .unwrap(),
                &metadata,
            );
            loaded_integrand.warm_up(&model).unwrap();
            let after_load = loaded_integrand
                .evaluate_samples_raw(
                    EvaluationTarget::Physical(&model),
                    &samples,
                    0,
                    false,
                    false,
                    Complex::new(F(0.0), F(0.0)),
                )
                .unwrap()
                .samples;
            assert_eq!(before_save.len(), after_load.len());
            for (before, after) in before_save.iter().zip(&after_load) {
                assert_eq!(before.integrand_result, after.integrand_result);
            }
            assert_eq!(
                assert_runtime_decomposition(&after_load, &metadata, &target_component_ids),
                before_occurrences,
            );
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn standalone_cut_sampling_compiles_from_production_cut_and_mass_data() {
    std::thread::Builder::new()
        .name("physical-cut-sampling".to_owned())
        .stack_size(64 * 1024 * 1024)
        .spawn(|| {
            use crate::integrands::process::{GraphTerm, SamplingChannelBridgeAcceptanceReport};
            use crate::settings::runtime::{ParameterizationSettings, SamplingChannelDefinition};

            test_initialise().unwrap();
            let model = Model::from_str(SCALARS_2P_3P_MODEL.to_string(), "json").unwrap();
            for (pdg, mass) in [(1001, 1.0), (1000, 0.0)] {
                let graph: Graph = r#"digraph cut_sampling_bubble {
                num=1; edge [pdg=PDG]; node [num=1];
                ext [style=invis, is_cut=0];
                ext -> a [id=0];
                a -> b [id=1, lmb_id=0];
                a -> b [id=2];
                b -> ext [id=3];
            }"#
                .replace("PDG", &pdg.to_string())
                .as_str()
                .into_graph(&model)
                .unwrap();
                let generation = generation_settings();
                let runtime = runtime_settings();
                let mut cross_section = preprocess_graph(graph, &model, &generation, &runtime);
                let definition = ProcessDefinition::from_graph_list(
                    std::slice::from_ref(&cross_section.supergraphs[0].graph),
                    GenerationType::CrossSection,
                    &model,
                )
                .unwrap();
                cross_section
                    .build_integrand(
                        &model,
                        &definition,
                        &GlobalSettings {
                            generation,
                            ..Default::default()
                        },
                        (&runtime).into(),
                        &generation_pool(),
                    )
                    .unwrap();
                let integrand = cross_section.integrand.as_mut().unwrap();
                integrand.warm_up(&model).unwrap();
                let ProcessIntegrand::CrossSection(integrand) = integrand else {
                    unreachable!()
                };
                let term = &integrand.data.graph_terms[0];
                let parent_lmb = term
                    .graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .map(|edge| edge.0)
                    .collect::<Vec<_>>();
                let (cut_id, cut_surface) = term.cut_esurface.iter_enumerated().next().unwrap();
                let cut_edges = cut_surface
                    .energies
                    .iter()
                    .map(|edge| edge.0.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                let mut parameterization = ParameterizationSettings::default();
                parameterization.sampling_channels.default_channel_selection =
                    vec!["physical_cut".into()];
                parameterization
                    .sampling_channels
                    .channel_definitions
                    .entry(term.graph.name.clone())
                    .or_default()
                    .insert(
                        "physical_cut".into(),
                        SamplingChannelDefinition {
                            around: format!("phase_space(cut({cut_edges}))"),
                            subspace_lmb: parent_lmb.clone(),
                            parent_lmb,
                            on_cut: vec![cut_id.0],
                            ..Default::default()
                        },
                    );
                let externals = runtime
                    .kinematics
                    .externals
                    .get_dependent_externals::<f64>(
                        crate::DependentMomentaConstructor::CrossSection,
                    )
                    .unwrap()
                    .iter()
                    .map(|momentum| {
                        [
                            momentum.temporal.value.0,
                            momentum.spatial.px.0,
                            momentum.spatial.py.0,
                            momentum.spatial.pz.0,
                        ]
                    })
                    .collect::<Vec<_>>();
                let bridge = term
                    .compile_sampling_bridge(&parameterization, &runtime, &externals, None)
                    .unwrap();
                assert_eq!(bridge.channels().len(), 1);
                assert_eq!(bridge.channels()[0].name, "physical_cut");
                // At the analytic two-particle threshold the LU scale is one.
                // This checks the graph's actual masses and cut equation, which
                // Gaussian normalization alone cannot distinguish from a chart
                // accidentally centered on some other radial shell.
                let masses = term.graph.get_real_mass_vector::<f64>(&model);
                assert!(
                    cut_surface
                        .energies
                        .iter()
                        .all(|edge| masses[*edge] == F(mass))
                );
                let threshold_radius = ((runtime.kinematics.e_cm / 2.0).powi(2)
                    - masses[cut_surface.energies[0]].0.powi(2))
                .sqrt();
                let shell_coordinate = threshold_radius
                    / (threshold_radius + runtime.kinematics.e_cm * parameterization.b);
                let shell = bridge
                    .forward(
                        crate::integrands::process::SamplingChannelId::from(0),
                        &[shell_coordinate, 0.31, 0.47],
                    )
                    .unwrap();
                let shell_loops = crate::momentum::sample::LoopMomenta::from_iter([
                    crate::momentum::ThreeMomentum::new(
                        F(shell.raw_coordinates[0]),
                        F(shell.raw_coordinates[1]),
                        F(shell.raw_coordinates[2]),
                    ),
                ]);
                let center = crate::momentum::sample::LoopMomenta::from_iter([
                    crate::momentum::ThreeMomentum::new(F(0.0), F(0.0), F(0.0)),
                ]);
                let prepared_externals = runtime
                    .kinematics
                    .externals
                    .get_dependent_externals::<f64>(
                        crate::DependentMomentaConstructor::CrossSection,
                    )
                    .unwrap();
                let (shell_value, derivative) = cut_surface.compute_self_and_r_derivative(
                    &F(1.0),
                    &shell_loops,
                    &center,
                    &prepared_externals,
                    &masses,
                    &term.graph.loop_momentum_basis,
                );
                assert!(shell_value.0.abs() < 1.0e-8, "{shell_value}");
                assert!(derivative.0 > 0.0);
                let report = SamplingChannelBridgeAcceptanceReport::normalized_gaussian(
                    &bridge,
                    4096,
                    1.0,
                    &vec![0.0; bridge.dimensions()],
                )
                .unwrap();
                assert!((report.normalization - 1.0).abs() < 0.02, "{report:?}");
                assert!(report.round_trip_residual_max < 1.0e-8, "{report:?}");

                // At fixed angular shape, simple-cut LU physics carries h(t),
                // whereas this proposal carries p(t|R). Their ratio is not
                // constant; multiplying the complete estimator by p/h must be.
                // The independent analytic density also checks that warmup uses
                // the actual runtime h family, width and every supported power.
                if mass == 1.0 {
                    use crate::settings::runtime::{
                        HFunction, HFunctionSettings, MultiChannelingSettings,
                        SamplingRadialProfile, SamplingSettings, StabilityLevelSetting,
                    };
                    for function in [
                        HFunction::Exponential,
                        HFunction::PolyExponential,
                        HFunction::PolyLeftRightExponential,
                    ] {
                        for power in [
                            None,
                            Some(0),
                            Some(1),
                            Some(3),
                            Some(4),
                            Some(6),
                            Some(7),
                            Some(9),
                            Some(10),
                            Some(12),
                            Some(13),
                            Some(15),
                            Some(16),
                        ] {
                            if function == HFunction::Exponential && power.is_some() {
                                continue;
                            }
                            for sigma in [0.6, 1.4] {
                                let h_settings = HFunctionSettings {
                                    function: function.clone(),
                                    sigma,
                                    power,
                                    ..Default::default()
                                };
                                let mut matched = ProcessIntegrand::CrossSection(integrand.clone());
                                let mut parameterization = parameterization.clone();
                                let profile = SamplingRadialProfile::default();
                                parameterization
                                    .sampling_channels
                                    .channel_definitions
                                    .get_mut(&term.graph.name)
                                    .unwrap()
                                    .get_mut("physical_cut")
                                    .unwrap()
                                    .radial_profile = Some(profile.clone());
                                let settings = matched.get_mut_settings();
                                settings.lu_h_function = h_settings.clone();
                                settings.sampling =
                                    SamplingSettings::MultiChanneling(MultiChannelingSettings {
                                        parameterization_settings: parameterization.clone(),
                                        ..Default::default()
                                    });
                                settings.stability.rotation_axis.clear();
                                settings.stability.levels =
                                    vec![StabilityLevelSetting::default_double()];
                                matched.warm_up(&model).unwrap();
                                let ProcessIntegrand::CrossSection(prepared) = &matched else {
                                    unreachable!()
                                };
                                let bridge = prepared.data.graph_terms[0]
                                    .multi_channeling_setup
                                    .sampling_bridge::<f64>()
                                    .unwrap()
                                    .clone();
                                let (scale, shape) = if function == HFunction::Exponential {
                                    (sigma / 2.0, 1.0)
                                } else {
                                    let exponent = if function == HFunction::PolyExponential {
                                        2.0
                                    } else {
                                        1.0
                                    };
                                    let shift =
                                        (1.0 - power.unwrap_or(0) as f64) / (2.0 * exponent);
                                    (
                                        sigma * (shift.asinh() / exponent).exp(),
                                        2.0 * exponent * (1.0 + shift * shift).powf(0.25),
                                    )
                                };
                                let mut ratios = Vec::new();
                                for u in [0.23, 0.47, 0.71] {
                                    let coordinates = [u, 0.31, 0.47];
                                    let point = bridge
                                        .forward(
                                            crate::integrands::process::SamplingChannelId(0),
                                            &coordinates,
                                        )
                                        .unwrap();
                                    let radius = point
                                        .raw_coordinates
                                        .iter()
                                        .map(|value| value * value)
                                        .sum::<f64>()
                                        .sqrt();
                                    let t = threshold_radius / radius;
                                    let z = t / scale;
                                    let a = threshold_radius
                                        / (runtime.kinematics.e_cm * parameterization.b);
                                    let density = (1.0 - profile.broad_fraction) * shape / scale
                                        * z.powf(shape - 1.0)
                                        / (1.0 + z.powf(shape)).powi(2)
                                        + profile.broad_fraction * a / (a + t).powi(2);
                                    let value = matched
                                        .evaluate_samples_raw(
                                            EvaluationTarget::Physical(&model),
                                            &[Sample::Continuous(
                                                F(1.0),
                                                coordinates.into_iter().map(F).collect(),
                                            )],
                                            0,
                                            false,
                                            false,
                                            Complex::new(F(0.0), F(0.0)),
                                        )
                                        .unwrap()
                                        .samples
                                        .remove(0);
                                    assert_eq!(value.parameterization_jacobian, Some(F(1.0)));
                                    let h = crate::utils::h(&F(t), None, None, &h_settings).0;
                                    ratios.push([
                                        value.integrand_result.re.0 * density / h,
                                        value.integrand_result.im.0 * density / h,
                                    ]);
                                }
                                let norm = ratios[0][0].hypot(ratios[0][1]);
                                assert!(
                                    norm > 0.0 && norm.is_finite(),
                                    "{h_settings:?}: {ratios:?}"
                                );
                                for ratio in &ratios[1..] {
                                    assert!(
                                        (ratio[0] - ratios[0][0]).hypot(ratio[1] - ratios[0][1])
                                            < 2.0e-7 * norm,
                                        "{h_settings:?}: {ratios:?}"
                                    );
                                }
                            }
                        }
                    }
                }

                // The same actual cut chart is used by the complete physical
                // source/stability route. The binary64 draw is mapped once in Arb,
                // so a collapsed Double map radius no longer requires a Quad
                // physical body. Compare the retained draw with Double-only
                // evaluation, including the complete physical event weights.
                {
                    use crate::{
                        integrands::evaluation::{PreciseEvaluationResult, StabilityStatus},
                        settings::runtime::{
                            MultiChannelingSettings, SamplingSettings, StabilityLevelSetting,
                        },
                    };
                    let mut runtime = ProcessIntegrand::CrossSection(integrand.clone());
                    let mut focused = parameterization.clone();
                    focused.power = 2.0;
                    let settings = runtime.get_mut_settings();
                    settings.sampling =
                        SamplingSettings::MultiChanneling(MultiChannelingSettings {
                            parameterization_settings: focused,
                            ..Default::default()
                        });
                    settings.stability.rotation_axis.clear();
                    settings.stability.levels = vec![
                        StabilityLevelSetting::default_double(),
                        StabilityLevelSetting::default_quad(),
                    ];
                    runtime.warm_up(&model).unwrap();
                    let source = Sample::Continuous(
                        F(1.0),
                        vec![F(shell_coordinate + 1.0e-10), F(0.31), F(0.47)],
                    );
                    let evaluated = runtime
                        .evaluate_sample_precise(
                            &source,
                            &model,
                            F(1.0),
                            false,
                            Complex::new_zero(),
                        )
                        .unwrap();
                    let PreciseEvaluationResult::Double(evaluated) = evaluated else {
                        panic!("the canonical cut draw must allow its regular Double physical body")
                    };
                    assert_eq!(evaluated.evaluation_metadata.stability_results.len(), 1);
                    assert!(!evaluated.evaluation_metadata.is_nan);
                    assert!(matches!(
                        evaluated.evaluation_metadata.stability_results[0].status,
                        StabilityStatus::Unknown
                    ));
                    assert!(
                        evaluated.integrand_result.re.0.is_finite()
                            && evaluated.integrand_result.im.0.is_finite()
                    );
                    assert_ne!(
                        evaluated.integrand_result,
                        Complex::new_re(evaluated.integrand_result.re.zero())
                    );
                    assert!(
                        evaluated
                            .parameterization_jacobian
                            .as_ref()
                            .is_some_and(|jacobian| jacobian == &jacobian.one())
                    );
                    runtime.get_mut_settings().stability.levels =
                        vec![StabilityLevelSetting::default_double()];
                    runtime.warm_up(&model).unwrap();
                    let direct = runtime
                        .evaluate_sample_precise(
                            &source,
                            &model,
                            F(1.0),
                            false,
                            Complex::new_zero(),
                        )
                        .unwrap();
                    let PreciseEvaluationResult::Double(direct) = direct else {
                        unreachable!()
                    };
                    assert!(!direct.evaluation_metadata.is_nan);
                    assert_eq!(evaluated.integrand_result, direct.integrand_result);
                    let evaluated_events = evaluated
                        .event_groups
                        .iter()
                        .flat_map(|group| group.iter())
                        .collect::<Vec<_>>();
                    let direct_events = direct
                        .event_groups
                        .iter()
                        .flat_map(|group| group.iter())
                        .collect::<Vec<_>>();
                    assert!(!evaluated_events.is_empty());
                    assert_eq!(evaluated_events.len(), direct_events.len());
                    for (evaluated, direct) in evaluated_events.iter().zip(direct_events) {
                        assert_eq!(evaluated.cut_info.sampling_channel_id, Some(0));
                        assert_eq!(evaluated.cut_info.sampling_channel_edge_ids, None);
                        assert_eq!(evaluated.weight, direct.weight);
                        assert_eq!(
                            evaluated.additional_weights.weights,
                            direct.additional_weights.weights
                        );
                    }
                }

                // User cut selectors are validated against the production CutId,
                // even though there is only one canonical channel in this test.
                parameterization
                    .sampling_channels
                    .channel_definitions
                    .get_mut(&term.graph.name)
                    .unwrap()
                    .get_mut("physical_cut")
                    .unwrap()
                    .on_cut = vec![usize::MAX];
                assert!(
                    term.compile_sampling_bridge(&parameterization, &runtime, &externals, None,)
                        .unwrap_err()
                        .to_string()
                        .contains("on_cut")
                );
            }
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn lu_h_sampling_preserves_normalization_and_raised_derivative_integrals() {
    use crate::{
        settings::runtime::{HFunction, HFunctionSettings},
        utils::{
            h, h_dual,
            hyperdual_utils::{extract_t_derivatives, simple_n_deriv_shape},
        },
    };
    use symbolica::domains::dual::HyperDual;

    test_initialise().unwrap();
    let shape = HyperDual::<F<f64>>::new(simple_n_deriv_shape(4));
    for function in [
        HFunction::Exponential,
        HFunction::PolyExponential,
        HFunction::PolyLeftRightExponential,
    ] {
        for power in [
            None,
            Some(0),
            Some(1),
            Some(3),
            Some(4),
            Some(6),
            Some(7),
            Some(9),
            Some(10),
            Some(12),
            Some(13),
            Some(15),
            Some(16),
        ] {
            if function == HFunction::Exponential && power.is_some() {
                continue;
            }
            for sigma in [0.4, 1.7] {
                let settings = HFunctionSettings {
                    function: function.clone(),
                    sigma,
                    power,
                    ..Default::default()
                };
                // Integrate in log(t/sigma), independently of the sampling map.
                // The exponential family has nonzero density at t=0 and needs
                // a longer lower tail; both polynomial families suppress it.
                let lower = if function == HFunction::Exponential {
                    -30.0
                } else {
                    -8.0
                };
                let steps = 8192;
                let step = (8.0 - lower) / steps as f64;
                let mut integrals = [0.0; 5];
                for index in 0..=steps {
                    let t = F(sigma * (lower + index as f64 * step).exp());
                    let derivatives =
                        extract_t_derivatives(h_dual(&shape.variable(0, t), None, None, &settings));
                    assert!(
                        (derivatives[0].0 - h(&t, None, None, &settings).0).abs() < 1.0e-12 / sigma
                    );
                    let quadrature = if index == 0 || index == steps {
                        1.0
                    } else if index % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    for (order, derivative) in derivatives.iter().enumerate() {
                        integrals[order] +=
                            quadrature * t.0.powi(order as i32 + 1) * derivative.0 * step / 3.0;
                    }
                }
                // Existing h tables contain rounded binary64 normalizations;
                // their agreement with one is independent of the IBP identity.
                assert!(
                    (integrals[0] - 1.0).abs() < 2.0e-8,
                    "{settings:?}: {integrals:?}"
                );
                let mut factorial = 1.0;
                for order in 1..=4 {
                    factorial *= -(order as f64);
                    assert!(
                        (integrals[order] / factorial - integrals[0]).abs() < 2.0e-8,
                        "{settings:?}, derivative {order}: {integrals:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn conditional_cut_sampling_preserves_both_sides_and_raised_sum() {
    std::thread::Builder::new()
        .name("conditional-cut-sampling".into())
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            use crate::{
                DependentMomentaConstructor,
                graph::LmbIndex,
                integrands::process::{
                    GaussianReferenceFunction, GraphTerm, MomentumSpaceEvaluationInput,
                    SamplingChannelBridge, SamplingChannelBridgeAcceptanceReport, SamplingChannelId,
                    sampling_maps::{SamplingEvaluationError, SamplingMapAffine},
                },
                momentum::{
                    ThreeMomentum,
                    sample::{LoopMomenta, SubspaceData},
                },
                settings::runtime::{
                    DiscreteGraphSamplingSettings, DiscreteGraphSamplingType, MultiChannelingSettings,
                    ParameterizationSettings, SamplingChannelDefinition, SamplingRadialProfile,
                    SamplingSettings, StabilityLevelSetting,
                },
                utils::{
                    ArbPrec,
                    newton_solver::{RadialRootDiagnostics, RadialRootIdentity},
                },
            };
            use itertools::Itertools;
            use symbolica::prelude::SingleFloat;
            use typed_index_collections::TiVec;

            test_initialise().unwrap();
            let model = Model::from_str(SCALARS_2P_3P_MODEL.to_string(), "json").unwrap();
            let graph = TRIPLE_DOTTED_BUBBLE.into_graph(&model).unwrap();
            let generation = generation_settings();
            let runtime = runtime_settings();
            let mut cross_section = preprocess_graph(graph, &model, &generation, &runtime);
            let definition = ProcessDefinition::from_graph_list(
                std::slice::from_ref(&cross_section.supergraphs[0].graph),
                GenerationType::CrossSection,
                &model,
            )
            .unwrap();
            cross_section
                .build_integrand(
                    &model,
                    &definition,
                    &GlobalSettings {
                        generation,
                        ..Default::default()
                    },
                    (&runtime).into(),
                    &generation_pool(),
                )
                .unwrap();
            let mut integrand = cross_section.integrand.take().unwrap();
            integrand.warm_up(&model).unwrap();

            for (parent, boost) in [(vec![1, 4, 7], 0.0), (vec![3, 6, 9], 1.0)] {
                // The second parent has L=Q-k, so BQ is nonzero in a boosted
                // frame. It catches a mistaken L_phys=t*L_raw rescaling.
                let mut runtime = runtime.clone();
                runtime.kinematics.externals = toml::from_str(&format!(
                    "type='constant'\n[data]\nmomenta=[[5.0,{boost},0.0,0.0]]\nhelicities=['summed_averaged']"
                )).unwrap();
                let externals = runtime
                    .kinematics
                    .externals
                    .get_dependent_externals::<f64>(DependentMomentaConstructor::CrossSection)
                    .unwrap();
                let external_arrays = externals
                    .iter()
                    .map(|p| {
                        [
                            p.temporal.value.0,
                            p.spatial.px.0,
                            p.spatial.py.0,
                            p.spatial.pz.0,
                        ]
                    })
                    .collect_vec();
                let ProcessIntegrand::CrossSection(prepared) = &integrand else {
                    unreachable!()
                };
                let term = &prepared.data.graph_terms[0];
                let (host_id, host) = term
                    .cut_esurface
                    .iter_enumerated()
                    .find(|(_, surface)| {
                        surface
                            .energies
                            .iter()
                            .map(|edge| edge.0)
                            .sorted()
                            .eq([4, 6])
                    })
                    .unwrap();
                assert!(term.cut_group_data.cut_groups.iter().any(|group| {
                    group.cuts.contains(&host_id) && group.related_esurface_group.max_occurence == 2
                }));
                let associations = &term.cut_threshold_associations[host_id];
                let left =
                    &term.topological_threshold_esurfaces[associations.left[0].topological_threshold_id];
                let right =
                    &term.topological_threshold_esurfaces[associations.right[0].topological_threshold_id];
                let edge_string = |surface: &crate::cff::esurface::Esurface| {
                    surface
                        .energies
                        .iter()
                        .map(|edge| edge.0)
                        .sorted()
                        .join(",")
                };
                let cut_edges = edge_string(host);
                let left_edges = edge_string(left);
                let right_edges = edge_string(right);
                let graph_name = term.graph.name.clone();
                let lmbs = TiVec::from(vec![
                    term.multi_channeling_setup
                        .sampling_parent_lmb(&parent)
                        .unwrap(),
                ]);
                let basis_id = LmbIndex::from(0);
                let lmb = &lmbs[basis_id];
                let frame = term
                    .multi_channeling_setup
                    .lmb_frame_map(lmb, &external_arrays)
                    .unwrap();
                let native_origin = frame.inverse(&[0.0; 9], &[]).unwrap().coordinates;
                assert_eq!(
                    native_origin,
                    vec![boost, 0.0, 0.0, boost, 0.0, 0.0, boost, 0.0, 0.0]
                );
                let masses = term.graph.get_real_mass_vector::<f64>(&model);
                let zero_loops =
                    LoopMomenta::from_iter((0..3).map(|_| ThreeMomentum::new(F(0.0), F(0.0), F(0.0))));
                let to_loops = |point: &[f64]| {
                    LoopMomenta::from_iter(
                        point
                            .chunks_exact(3)
                            .map(|v| ThreeMomentum::new(F(v[0]), F(v[1]), F(v[2]))),
                    )
                };
                let lu = |point: &[f64]| {
                    let loops = to_loops(point);
                    let (guess, _) =
                        host.get_radius_guess(&loops, &externals, &term.graph.loop_momentum_basis);
                    RadialRootDiagnostics::default()
                        .solve(
                            &RadialRootIdentity::new("conditional fixture host".into()),
                            &F(0.0),
                            &guess,
                            |t| {
                                host.compute_self_and_r_derivative(
                                    t,
                                    &loops,
                                    &zero_loops,
                                    &externals,
                                    &masses,
                                    &term.graph.loop_momentum_basis,
                                )
                            },
                            &F(1.0),
                            2000,
                            64,
                            &F(5.0),
                        )
                        .unwrap()
                        .solution
                        .0
                };

                let mut parameterization = ParameterizationSettings::default();
                parameterization.sampling_channels.default_channel_selection =
                    vec!["left_only".into(), "both".into()];
                parameterization.sampling_channels.channel_definitions.insert(graph_name.clone(), BTreeMap::from([
                    ("left_only".into(), SamplingChannelDefinition {
                        around: format!("then(block(lmb({}),phase_space(cut({cut_edges}))),block(lmb({}),left(surface({left_edges}))),lmb({}))", parent[1], parent[0], parent[2]),
                        parent_lmb: parent.clone(), on_cut: vec![host_id.0], ..Default::default()
                    }),
                    ("both".into(), SamplingChannelDefinition {
                        around: format!("then(block(lmb({}),phase_space(cut({cut_edges}))),block(lmb({}),left(surface({left_edges}))),block(lmb({}),right(surface({right_edges}))))", parent[1], parent[0], parent[2]),
                        parent_lmb: parent.clone(), on_cut: vec![host_id.0], ..Default::default()
                    }),
                ]));
                let bridge = term
                    .compile_sampling_bridge(&parameterization, &runtime, &external_arrays, None)
                    .unwrap();
                // An explicitly hosted graph surface remains a sampling
                // target without any host-side CT association. A side-qualified
                // spelling still enforces that side's physical identity.
                let mut direct_term = term.clone();
                direct_term.cut_threshold_associations[host_id].left.clear();
                let mut direct_settings = parameterization.clone();
                direct_settings.sampling_channels.default_channel_selection = vec!["direct".into()];
                direct_settings.sampling_channels.channel_definitions.get_mut(&graph_name).unwrap().insert("direct".into(), SamplingChannelDefinition {
                    around: format!("then(block(lmb({}),phase_space(cut({cut_edges}))),block(lmb({}),surface({left_edges})),lmb({}))", parent[1], parent[0], parent[2]),
                    parent_lmb: parent.clone(), on_cut: vec![host_id.0], ..Default::default()
                });
                direct_term
                    .compile_sampling_bridge(&direct_settings, &runtime, &external_arrays, None)
                    .unwrap();
                let error = direct_term
                    .compile_sampling_bridge(&parameterization, &runtime, &external_arrays, None)
                    .unwrap_err();
                assert!(
                    error.to_string().contains("target") && error.to_string().contains("host cut"),
                    "{error:?}"
                );
                let mut invalid = direct_settings.clone();
                invalid
                    .sampling_channels
                    .channel_definitions
                    .get_mut(&graph_name)
                    .unwrap()
                    .get_mut("direct")
                    .unwrap()
                    .around = format!(
                    "then(block(lmb({}),phase_space(cut({cut_edges}))),block(lmb({}),surface({left_edges})),lmb({}))",
                    parent[0], parent[1], parent[2]
                );
                let error = term
                    .compile_sampling_bridge(&invalid, &runtime, &external_arrays, None)
                    .unwrap_err();
                assert!(
                    error
                        .to_string()
                        .contains("phase-space cut depends on omitted/active coordinates"),
                    "{error:?}"
                );
                assert_eq!(bridge.dimensions(), 9);
                assert_eq!(bridge.channels().len(), 2);
                let coordinates = [0.31, 0.27, 0.61, 0.42, 0.33, 0.73, 0.57, 0.23, 0.67];
                let channel = SamplingChannelId(1);
                let mapped = bridge.forward(channel, &coordinates).unwrap();
                // Both side consumers share the same literal host-dependent
                // prior, even though the second sees additional host-null input.
                assert_eq!(mapped.prepared_lu_hosts.len(), 1);
                let retained = &mapped.prepared_lu_hosts[0];
                assert_eq!(retained.plan.parent_lmb, parent);
                assert_eq!(retained.plan.required_prior_lmb, vec![parent[1]]);
                assert_eq!(retained.source.generating_channel, channel);
                assert_eq!(retained.source.target_channel, channel);
                let mut changed_null = coordinates;
                changed_null[3] = 0.36;
                changed_null[7] = 0.41;
                let changed = bridge.forward(channel, &changed_null).unwrap();
                assert_eq!(changed.prepared_lu_hosts.len(), 1);
                assert_eq!(retained.prior, changed.prepared_lu_hosts[0].prior);
                assert_eq!(retained.solution.solution, changed.prepared_lu_hosts[0].solution.solution);
                for t in [F(0.7), F(1.0), F(1.3)] {
                    assert_eq!(retained.ray.evaluate(&t), changed.prepared_lu_hosts[0].ray.evaluate(&t));
                }
                let recovered_host = bridge.inverse(channel, &mapped.raw_coordinates).unwrap().unwrap();
                assert_eq!(recovered_host.prepared_lu_hosts.len(), 1);
                let inverse_native = frame.inverse(&mapped.raw_coordinates, &[]).unwrap().coordinates;
                assert_eq!(recovered_host.prepared_lu_hosts[0].prior, inverse_native[3..6]);
                // Inverse work never overwrites the initial forward source.
                assert_eq!(mapped.prepared_lu_hosts[0].prior, retained.prior);
                if boost != 0.0 {
                    // An ordinary raw prefix does not first add BQ. Its small
                    // native component therefore exposes actual subtraction
                    // roundoff through native -> master -> native, unlike a
                    // cut prefix already rounded by the same BQ translation.
                    let mut ordinary_prefix = parameterization.clone();
                    ordinary_prefix.sampling_channels.default_channel_selection = vec!["ordinary_host".into()];
                    ordinary_prefix.sampling_channels.channel_definitions.get_mut(&graph_name).unwrap().insert(
                        "ordinary_host".into(), SamplingChannelDefinition {
                            around: format!("then(lmb({}),block(lmb({}),at_cut(cut({cut_edges}),left(surface({left_edges})))),lmb({}))",
                                parent[1], parent[0], parent[2]),
                            parent_lmb: parent.clone(), on_cut: vec![host_id.0], ..Default::default()
                        });
                    let ordinary_bridge = term.compile_sampling_bridge(
                        &ordinary_prefix, &runtime, &external_arrays, None).unwrap();
                    let witness = (1..=8).map(|index| {
                        let mut cube = coordinates;
                        cube[0] = 0.025 * index as f64;
                        let forward = ordinary_bridge.forward(SamplingChannelId(0), &cube).unwrap();
                        let inverse = ordinary_bridge.inverse(SamplingChannelId(0), &forward.raw_coordinates)
                            .unwrap().unwrap();
                        (forward, inverse)
                    }).find(|(forward, inverse)|
                        forward.prepared_lu_hosts[0].prior != inverse.prepared_lu_hosts[0].prior)
                        .expect("boosted ordinary prefix must exhibit an actual native affine roundoff");
                    let (forward, inverse) = witness;
                    assert_eq!(forward.prepared_lu_hosts.len(), 1);
                    assert_eq!(inverse.prepared_lu_hosts.len(), 1);
                    let expected_prior = frame.inverse(&forward.raw_coordinates, &[]).unwrap().coordinates;
                    assert_eq!(inverse.prepared_lu_hosts[0].prior, expected_prior[3..6]);
                    assert_ne!(forward.prepared_lu_hosts[0].prior, inverse.prepared_lu_hosts[0].prior);
                    // Partition-time own inversion and a fresh direct inversion
                    // both evaluate q at the actual supplied point/context.
                    let inverse_log_density = inverse.map.inverse_jacobian.ln();
                    assert!((forward.partition.log_scores[0].unwrap() - inverse_log_density).abs() < 1.0e-12);
                }
                let tau = lu(&mapped.raw_coordinates);
                let physical_master = mapped.raw_coordinates.iter().map(|k| k * tau).collect_vec();
                let physical_native = frame.inverse(&physical_master, &[]).unwrap().coordinates;
                let master_loops = to_loops(&physical_master);
                let native_loops = to_loops(&physical_native);
                let host_value = host
                    .compute_self_and_r_derivative(
                        &F(1.0),
                        &master_loops,
                        &zero_loops,
                        &externals,
                        &masses,
                        &term.graph.loop_momentum_basis,
                    )
                    .0
                    .0;
                assert!(host_value.abs() < 2.0e-10, "{host_value}");
                for (surface, active_edge) in [(left, parent[0]), (right, parent[2])] {
                    let subspace = SubspaceData::new_from_parent_basis_edges(
                        &[EdgeIndex(active_edge)],
                        &term.graph.full_filter(),
                        basis_id,
                        &term.graph,
                        &lmbs,
                    )
                    .unwrap();
                    let global = surface
                        .compute_self_and_r_derivative(
                            &F(1.0),
                            &master_loops,
                            &zero_loops,
                            &externals,
                            &masses,
                            &term.graph.loop_momentum_basis,
                        )
                        .0
                        .0;
                    let native_global = surface
                        .compute_self_and_r_derivative(
                            &F(1.0),
                            &native_loops,
                            &zero_loops,
                            &externals,
                            &masses,
                            lmb,
                        )
                        .0
                        .0;
                    let subspace_value = surface
                        .compute_self_and_r_derivative_subspace(
                            &F(1.0),
                            &native_loops,
                            &zero_loops,
                            &externals,
                            &masses,
                            &subspace,
                            &lmbs,
                            &term.graph,
                        )
                        .0
                        .0;
                    let spatial: TiVec<crate::momentum::sample::ExternalIndex, _> =
                        externals.iter().map(|p| p.spatial).collect();
                    let energy_sum = |surface: &crate::cff::esurface::Esurface| {
                        surface
                            .energies
                            .iter()
                            .map(|edge| {
                                let p =
                                    lmb.edge_signatures[*edge].compute_momentum(&native_loops, &spatial);
                                (p.norm_squared().0 + masses[*edge].0.powi(2)).sqrt()
                            })
                            .sum::<f64>()
                    };
                    let boundary = energy_sum(surface) - energy_sum(host)
                        + surface.compute_shift_part_from_momenta(&externals, lmb).0
                        - host.compute_shift_part_from_momenta(&externals, lmb).0;
                    let index = parent.iter().position(|edge| *edge == active_edge).unwrap();
                    let p = &native_loops.0[index];
                    let q_minus_p = &externals[crate::momentum::sample::ExternalIndex(0)].spatial - p;
                    let analytic =
                        (p.norm_squared().0 + 1.0).sqrt() + (q_minus_p.norm_squared().0 + 1.0).sqrt() - 5.0;
                    assert!((global - analytic).abs() < 1.0e-9);
                    assert!((global - native_global).abs() < 1.0e-9);
                    assert!((global - subspace_value).abs() < 1.0e-9);
                    assert!((global - boundary).abs() < 1.0e-9);
                }
                let mut changed = coordinates;
                changed[3] += 0.07;
                changed[6] -= 0.09;
                let other = bridge.forward(channel, &changed).unwrap();
                assert!((lu(&other.raw_coordinates) - tau).abs() < 1.0e-10);
                assert_eq!(&mapped.raw_coordinates[3..6], &other.raw_coordinates[3..6]);

                // Full 9D finite differences include cut-to-side off-diagonal
                // derivatives; multiplying only independent radial Jacobians
                // would miss the two physical-to-raw t^(-3) factors.
                let check_jacobian = |bridge: &SamplingChannelBridge| {
                    let mapped = bridge.forward(channel, &coordinates).unwrap();
                    let step = 2.0e-6;
                    let mut jacobian = vec![vec![0.0; 9]; 9];
                    for axis in 0..9 {
                        let mut plus = coordinates;
                        let mut minus = coordinates;
                        plus[axis] += step;
                        minus[axis] -= step;
                        let plus = bridge.forward(channel, &plus).unwrap();
                        let minus = bridge.forward(channel, &minus).unwrap();
                        for (row, values) in jacobian.iter_mut().enumerate() {
                            values[axis] =
                                (plus.raw_coordinates[row] - minus.raw_coordinates[row]) / (2.0 * step);
                        }
                    }
                    assert!(jacobian[0][0].abs() > 1.0e-3);
                    let finite_difference = SamplingMapAffine::new(jacobian, vec![0.0; 9])
                        .unwrap()
                        .determinant();
                    assert!(
                        (finite_difference / mapped.map.jacobian - 1.0).abs() < 3.0e-5,
                        "{finite_difference} vs {}",
                        mapped.map.jacobian
                    );
                };
                check_jacobian(&bridge);
                // Attach the actual runtime LU h profile to the reduced cut
                // block, then retain both prepared side maps. Profile-free and
                // matched variants coexist without overwriting shared geometry.
                let mut matched = parameterization.clone();
                matched
                    .sampling_channels
                    .channel_definitions
                    .get_mut(&graph_name)
                    .unwrap()
                    .get_mut("both")
                    .unwrap()
                    .radial_profile = Some(SamplingRadialProfile::default());
                let matched_bridge = term
                    .compile_sampling_bridge(&matched, &runtime, &external_arrays, None)
                    .unwrap();
                let matched_point = matched_bridge.forward(channel, &coordinates).unwrap();
                assert!(
                    matched_point
                        .map
                        .diagnostics
                        .iter()
                        .any(|value| value.contains("max_occurrence=2")),
                    "{:?}",
                    matched_point.map.diagnostics
                );
                assert_ne!(matched_point.raw_coordinates, mapped.raw_coordinates);
                check_jacobian(&matched_bridge);
                for id in [SamplingChannelId(0), SamplingChannelId(1)] {
                    let inverse = matched_bridge
                        .inverse(id, &matched_point.raw_coordinates)
                        .unwrap().expect("full-support map must contain the supplied point");
                    let recovered = matched_bridge
                        .forward(id, &inverse.map.coordinates)
                        .unwrap();
                    assert!(
                        recovered
                            .raw_coordinates
                            .iter()
                            .zip(&matched_point.raw_coordinates)
                            .all(|(a, b)| (a - b).abs() < 1.0e-7)
                    );
                }
                for id in [SamplingChannelId(0), SamplingChannelId(1)] {
                    let inverse = bridge.inverse(id, &mapped.raw_coordinates).unwrap().expect("full-support map must contain the supplied point");
                    let recovered = bridge.forward(id, &inverse.map.coordinates).unwrap();
                    assert!(
                        recovered
                            .raw_coordinates
                            .iter()
                            .zip(&mapped.raw_coordinates)
                            .all(|(a, b)| (a - b).abs() < 1.0e-7)
                    );
                }
                let report = SamplingChannelBridgeAcceptanceReport::normalized_gaussian(
                    &bridge,
                    8192,
                    1.5,
                    &[0.2, -0.1, 0.3, 0.1, 0.2, -0.1, -0.2, 0.1, 0.2],
                )
                .unwrap();
                assert!((report.normalization - 1.0).abs() < 0.06, "{report:?}");
                assert!(
                    (report.second_moment / report.expected_second_moment - 1.0).abs() < 0.08,
                    "{report:?}"
                );

                if boost == 0.0 {
                    // Finite I/t entries can still have an unrepresentable
                    // determinant. This must request precision rescue, never
                    // declare an absent surface or accept zero proposal support.
                    let tiny = [0.3, 0.2, 0.1, 1.0, -0.5, 0.2, 0.7, 0.1, -0.3].map(|v| v * 1.0e-110);
                    let error = bridge.inverse(channel, &tiny).unwrap_err();
                    assert!(
                        matches!(
                            error.downcast_ref::<SamplingEvaluationError>(),
                            Some(SamplingEvaluationError::Unrepresentable { .. })
                        ),
                        "{error:?}"
                    );
                    let native_externals = external_arrays
                        .iter()
                        .map(|p| p.map(|v| F::<ArbPrec>::from_f64(v).0))
                        .collect_vec();
                    let native_bridge = term
                        .compile_sampling_bridge(&parameterization, &runtime, &native_externals, None)
                        .unwrap();
                    let native_tiny = tiny.map(|v| F::<ArbPrec>::from_f64(v).0);
                    let recovered = native_bridge.inverse(channel, &native_tiny).unwrap().expect("full-support map must contain the supplied point");
                    assert!(
                        recovered.map.jacobian.is_finite() && recovered.map.inverse_jacobian.is_finite()
                    );
                    assert!(recovered.map.jacobian > F::<ArbPrec>::default().zero().0);
                }

                // The same complete physical raised-cut/CT estimator is used
                // by summed channels, explicit channel MC, and direct momenta.
                runtime.sampling = SamplingSettings::MultiChanneling(MultiChannelingSettings {
                    parameterization_settings: parameterization.clone(),
                    ..Default::default()
                });
                runtime.stability.rotation_axis.clear();
                runtime.stability.levels = vec![StabilityLevelSetting::default_double()];
                *integrand.get_mut_settings() = runtime.clone();
                integrand.warm_up(&model).unwrap();
                let sample = Sample::Continuous(F(1.0), coordinates.iter().copied().map(F).collect());
                let reference = GaussianReferenceFunction::new(1.5, vec![0.2; 9]).unwrap();
                let summed_reference = integrand
                    .evaluate_reference_sample_detailed(&sample, &reference)
                    .unwrap();
                let summed = integrand
                    .evaluate_samples_raw(
                        EvaluationTarget::Physical(&model),
                        std::slice::from_ref(&sample),
                        0,
                        false,
                        false,
                        Complex::new(F(0.0), F(0.0)),
                    )
                    .unwrap()
                    .samples
                    .remove(0);
                assert!(
                    !summed.evaluation_metadata.is_nan,
                    "{}",
                    summed.evaluation_metadata
                );
                assert!(summed.integrand_result.re.0.abs() + summed.integrand_result.im.0.abs() > 0.0);
                assert!(!summed.event_groups.is_empty());
                let mut direct_sum = Complex::new(F(0.0), F(0.0));
                let mut direct_results = Vec::new();
                for id in [SamplingChannelId(0), SamplingChannelId(1)] {
                    let point = bridge.forward(id, &coordinates).unwrap();
                    let direct = integrand
                        .evaluate_momentum_configuration(
                            &model,
                            &MomentumSpaceEvaluationInput {
                                loop_momenta: to_loops(&point.raw_coordinates).0,
                                integrator_weight: F(1.0),
                                graph_id: Some(0),
                                group_id: None,
                                orientation: None,
                                channel_id: None,
                            },
                            false,
                        )
                        .unwrap();
                    assert!(
                        !direct.evaluation_metadata.is_nan,
                        "{}",
                        direct.evaluation_metadata
                    );
                    assert!(!direct.event_groups.is_empty());
                    let factor = point.selected_factor().unwrap();
                    direct_sum += direct.integrand_result * F(factor);
                    direct_results.push((direct, factor));
                }
                for (a, b) in [
                    (summed.integrand_result.re.0, direct_sum.re.0),
                    (summed.integrand_result.im.0, direct_sum.im.0),
                ] {
                    assert!(
                        (a - b).abs() <= 1.0e-8 * a.abs().max(b.abs()).max(1.0e-25),
                        "{a} != {b}"
                    );
                }
                let SamplingSettings::MultiChanneling(multichanneling) = runtime.sampling.clone() else {
                    unreachable!()
                };
                integrand.get_mut_settings().sampling =
                    SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
                        sample_orientations: false,
                        sampling_type: DiscreteGraphSamplingType::SamplingMultiChanneling(multichanneling),
                        ..Default::default()
                    });
                integrand.warm_up(&model).unwrap();
                let mut explicit_sum = Complex::new(F(0.0), F(0.0));
                let mut reference_sum = 0.0;
                let mut moment_sum = 0.0;
                for (id, (direct, factor)) in direct_results.iter().enumerate() {
                    let discrete = Sample::Discrete(
                        F(1.0),
                        0,
                        Some(Box::new(Sample::Discrete(
                            F(1.0),
                            id,
                            Some(Box::new(sample.clone())),
                        ))),
                    );
                    if id == 1 {
                        let ProcessIntegrand::CrossSection(inner) = &mut integrand else { unreachable!() };
                        crate::integrands::process::tests::check_host_source_transport(
                            inner, &model, &discrete, 0, SamplingChannelId(id),
                        ).unwrap();
                    }
                    let reference = integrand
                        .evaluate_reference_sample_detailed(&discrete, &reference)
                        .unwrap();
                    reference_sum += reference.evaluation.integrand_result.re.0;
                    moment_sum += reference.moments.second_moment.0;
                    let selected = integrand
                        .evaluate_samples_raw(
                            EvaluationTarget::Physical(&model),
                            &[discrete],
                            0,
                            false,
                            false,
                            Complex::new(F(0.0), F(0.0)),
                        )
                        .unwrap()
                        .samples
                        .remove(0);
                    assert!(
                        !selected.evaluation_metadata.is_nan,
                        "{}",
                        selected.evaluation_metadata
                    );
                    let direct_events = direct
                        .event_groups
                        .iter()
                        .flat_map(|group| group.iter())
                        .collect_vec();
                    let selected_events = selected
                        .event_groups
                        .iter()
                        .flat_map(|group| group.iter())
                        .collect_vec();
                    assert!(!selected_events.is_empty());
                    assert_eq!(selected_events.len(), direct_events.len());
                    assert_eq!(selected.event_groups.len(), direct.event_groups.len());
                    for (selected, direct) in selected_events.iter().zip(direct_events) {
                        assert_eq!(selected.cut_info.cut_id, direct.cut_info.cut_id);
                        for (a, b) in [
                            (selected.weight.re.0, direct.weight.re.0 * factor),
                            (selected.weight.im.0, direct.weight.im.0 * factor),
                        ] {
                            assert!(
                                (a - b).abs() < 1.0e-8 * a.abs().max(b.abs()).max(1.0e-25),
                                "event {a} != {b}"
                            );
                        }
                    }
                    explicit_sum += selected.integrand_result;
                }
                assert!(
                    (reference_sum - summed_reference.evaluation.integrand_result.re.0).abs() < 1.0e-11
                );
                assert!((moment_sum - summed_reference.moments.second_moment.0).abs() < 1.0e-10);
                assert!(
                    (explicit_sum.re.0 - summed.integrand_result.re.0).abs()
                        < 1.0e-8 * summed.integrand_result.re.0.abs().max(1.0e-25)
                );
                assert!(
                    (explicit_sum.im.0 - summed.integrand_result.im.0).abs()
                        < 1.0e-8 * summed.integrand_result.im.0.abs().max(1.0e-25)
                );
            }
        }).unwrap().join().unwrap();
}
