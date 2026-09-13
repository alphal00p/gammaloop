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
                        &model,
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
                    &model,
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
                    .compile_sampling_bridge(
                        &parameterization,
                        runtime.kinematics.e_cm,
                        &externals,
                        None,
                    )
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

                // The same actual cut chart is used by the complete physical
                // source/stability route. Keep the draw binary64 and reconstruct
                // it natively after the focused radius collapses in Double.
                {
                    use crate::{
                        integrands::evaluation::{PreciseEvaluationResult, StabilityStatus},
                        settings::runtime::{
                            MultiChannelingSettings, SamplingSettings, StabilityLevelSetting,
                        },
                    };
                    use symbolica::prelude::SingleFloat;
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
                    let rescued = runtime
                        .evaluate_sample_precise(
                            &source,
                            &model,
                            F(1.0),
                            false,
                            Complex::new_zero(),
                        )
                        .unwrap();
                    let PreciseEvaluationResult::Quad(rescued) = rescued else {
                        panic!("the actual cut chart must reconstruct at Quad precision")
                    };
                    assert_eq!(rescued.evaluation_metadata.stability_results.len(), 2);
                    assert!(matches!(
                        rescued.evaluation_metadata.stability_results[0].status,
                        StabilityStatus::Unstable(0)
                    ));
                    assert!(
                        rescued.integrand_result.re.0.is_finite()
                            && rescued.integrand_result.im.0.is_finite()
                    );
                    assert_ne!(
                        rescued.integrand_result,
                        Complex::new_re(rescued.integrand_result.re.zero())
                    );
                    assert!(
                        rescued
                            .parameterization_jacobian
                            .as_ref()
                            .is_some_and(|jacobian| jacobian == &jacobian.one())
                    );
                    runtime.get_mut_settings().stability.levels =
                        vec![StabilityLevelSetting::default_quad()];
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
                    let PreciseEvaluationResult::Quad(direct) = direct else {
                        unreachable!()
                    };
                    assert_eq!(rescued.integrand_result, direct.integrand_result);
                    let rescued_events = rescued
                        .event_groups
                        .iter()
                        .flat_map(|group| group.iter())
                        .collect::<Vec<_>>();
                    let direct_events = direct
                        .event_groups
                        .iter()
                        .flat_map(|group| group.iter())
                        .collect::<Vec<_>>();
                    assert!(!rescued_events.is_empty());
                    assert_eq!(rescued_events.len(), direct_events.len());
                    for (rescued, direct) in rescued_events.iter().zip(direct_events) {
                        assert_eq!(rescued.cut_info.sampling_channel_id, Some(0));
                        assert_eq!(rescued.cut_info.sampling_channel_edge_ids, None);
                        assert_eq!(rescued.weight, direct.weight);
                        assert_eq!(
                            rescued.additional_weights.weights,
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
                    term.compile_sampling_bridge(
                        &parameterization,
                        runtime.kinematics.e_cm,
                        &externals,
                        None,
                    )
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
