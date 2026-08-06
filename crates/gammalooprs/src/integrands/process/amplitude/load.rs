//#!/usr/bin/env -S rust-script
//! ```cargo
//! [dependencies]
//! bincode = "2"
//! bincode-trait-derive = "0.1.1"
//! eyre = "0.6"
//! rand = "0.9"
//! serde_json = "1"
//! serde = { version = "1.0", features = ["derive"] }
//! symbolica = { git = "https://github.com/symbolica-dev/symbolica", rev = "0441bd7a511209dce2ca99925fe87f8b18e4bf03", default-features = false, features = ["bincode", "gmp", "native_code_generation", "serde"] }
//! [patch.crates-io]
//! numerica = { git = "https://github.com/symbolica-dev/symbolica", rev = "0441bd7a511209dce2ca99925fe87f8b18e4bf03" }
//! graphica = { git = "https://github.com/symbolica-dev/symbolica", rev = "0441bd7a511209dce2ca99925fe87f8b18e4bf03" }
//! ```

#![allow(dead_code)]
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::Cursor,
    ops::Neg,
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use bincode_trait_derive::{Decode, Encode};
use eyre::{Context, Result, eyre};
use rand::Rng;
use serde::{Deserialize, Serialize};
use symbolica::{
    domains::rational::Fraction, evaluate::JITCompiledEvaluator, prelude::*, state::StateMap,
};

use crate::integrands::process::cross_section::load::{
    LoadedStandaloneThresholdMultiplierCollection, StandaloneThresholdMultiplierCollectionArchive,
    build_threshold_multiplier_collection, validate_threshold_multiplier_archive,
};
use crate::processes::{
    StandaloneNumericTarget, ThresholdCountertermComponentKind,
    ThresholdCountertermMetadataRegistry, ThresholdCountertermOrigin, ThresholdCountertermSide,
};

pub const STANDALONE_EVALUATORS_VERSION: u32 = 8;
pub const STANDALONE_MODE_RUST: u8 = 0;

#[derive(
    Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Encode, Decode, Serialize, Deserialize,
)]
pub struct StandaloneCutCFFIndex {
    pub(crate) left_threshold_order: Option<usize>,
    pub(crate) right_threshold_order: Option<usize>,
    pub(crate) lu_cut_order: Option<usize>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneEvaluatorArchive<S = Vec<u8>, T = Vec<u8>> {
    pub(crate) version: u32,
    pub(crate) numeric_target: StandaloneNumericTarget,
    pub(crate) symbolica_state: S,
    pub(crate) graph_terms: Vec<StandaloneGraphTermArchive<T>>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneComplexInput {
    pub(crate) re: String,
    pub(crate) im: String,
}

impl StandaloneComplexInput {
    pub(crate) fn to_f64(&self) -> Result<Complex<f64>> {
        Ok(Complex::new(self.re.parse()?, self.im.parse()?))
    }
}
impl StandaloneEvaluatorArchive<(), String> {
    pub fn load(self) -> Result<LoadedStandaloneEvaluators> {
        if self.version != STANDALONE_EVALUATORS_VERSION {
            return Err(eyre!(
                "Unsupported version {} (expected {})",
                self.version,
                STANDALONE_EVALUATORS_VERSION
            ));
        }

        let mut symbolica_state = Vec::new();
        State::export(&mut symbolica_state)
            .with_context(|| "Failed to export Symbolica state for standalone evaluators")?;

        let mut state_cursor = Cursor::new(&symbolica_state);
        let state_map = State::import(&mut state_cursor, None)?;

        self.load_impl(&state_map)
    }
}

impl<S, A> StandaloneEvaluatorArchive<S, A>
where
    A: ImportWithMap
        + crate::integrands::process::cross_section::load::ImportWithMap
        + Clone
        + PartialEq,
{
    #[allow(clippy::type_complexity)]
    pub fn load_impl(self, state_map: &StateMap) -> Result<LoadedStandaloneEvaluators> {
        let mut graph_terms = Vec::new();

        for graph in self.graph_terms {
            let graph_name = graph.graph_name.as_str();
            let params = graph
                .param_builder_params
                .iter()
                .map(|b| ImportWithMap::import_with_map(b, state_map))
                .collect::<Result<Vec<_>>>()?;

            for p in params.iter() {
                println!("Loaded param builder param: {}", p);
            }
            let replacements = parse_fn_map_entries(&graph.fn_map_entries, state_map)?;
            let timed_build = |label: &str,
                               payload: StandaloneGenericEvaluatorArchive<A>,
                               iterate: bool|
             -> Result<LoadedGenericEvaluator> {
                let started = Instant::now();
                let evaluator =
                    build_evaluator(payload, &params, replacements.clone(), state_map, iterate)?;
                println!(
                    "[timing] build_evaluator {}::{} took {:?}",
                    graph_name,
                    label,
                    started.elapsed()
                );
                Ok(evaluator)
            };

            let (exprs, all_reps, single, result) = timed_build(
                "original.parametric",
                graph.original_integrand.single_parametric,
                false,
            )?;

            let iterative = graph
                .original_integrand
                .iterative
                .map(|payload| timed_build("original.iterative", payload, true))
                .transpose()?;

            let summed = graph
                .original_integrand
                .summed
                .map(|payload| timed_build("original.summed", payload, false))
                .transpose()?;

            let mut fnmap_integrand = None;

            let summed_fnmap = graph
                .original_integrand
                .summed_function_map
                .map(|payload| {
                    let (_, rhs, _, _) =
                        &parse_fn_map_entries(&payload.additional_fn_map_entries, state_map)?[0];
                    fnmap_integrand = Some(rhs.clone());
                    // parse_lit!(gammaloop::integrand(1,1,1,1,1,-1,1,-1,-1,-1,1,-1,-1))
                    timed_build("original.summed_fnmap", payload, false)
                })
                .transpose()?;

            if let Some(a) = fnmap_integrand {
                println!("Comparing fnmap summed fn and parametric epression");

                if a != exprs[0] {
                    println!("They are the different:\n {}!", (&a - &exprs[0]).expand());
                } else {
                    println!("They are the same!")
                }
            }

            let original_integrand = LoadedStandaloneEvaluatorStack {
                parametric: (exprs, all_reps, single, result),
                orientation_start: graph.original_integrand.start,
                mult_offset: graph.original_integrand.mult_offset,
                representative_input: graph
                    .original_integrand
                    .representative_input
                    .iter()
                    .map(StandaloneComplexInput::to_f64)
                    .collect::<Result<Vec<_>>>()?,
                iterative,
                summed,
                summed_fnmap,
            };
            let threshold_counterterms = graph
                .threshold_counterterms
                .into_iter()
                .enumerate()
                .map(|(ct_idx, ct_orders)| {
                    ct_orders
                        .into_iter()
                        .enumerate()
                        .map(|(slot_idx, ct)| {
                            let StandaloneIndexedEvaluatorStackArchive {
                                cut_cff_index,
                                evaluator_stack: ct,
                            } = ct;

                            let ct_parametric_label =
                                format!("ct[{ct_idx}].slot[{slot_idx}].parametric");
                            let parametric =
                                timed_build(&ct_parametric_label, ct.single_parametric, false)?;
                            let iterative = ct
                                .iterative
                                .map(|payload| {
                                    let label = format!("ct[{ct_idx}].slot[{slot_idx}].iterative");
                                    timed_build(&label, payload, true)
                                })
                                .transpose()?;
                            let summed = ct
                                .summed
                                .map(|payload| {
                                    let label = format!("ct[{ct_idx}].slot[{slot_idx}].summed");
                                    timed_build(&label, payload, false)
                                })
                                .transpose()?;
                            let summed_fnmap = ct
                                .summed_function_map
                                .map(|payload| {
                                    let label =
                                        format!("ct[{ct_idx}].slot[{slot_idx}].summed_fnmap");
                                    timed_build(&label, payload, false)
                                })
                                .transpose()?;

                            Ok((
                                cut_cff_index,
                                LoadedStandaloneEvaluatorStack {
                                    orientation_start: ct.start,
                                    mult_offset: ct.mult_offset,
                                    representative_input: ct
                                        .representative_input
                                        .iter()
                                        .map(StandaloneComplexInput::to_f64)
                                        .collect::<Result<Vec<_>>>()?,
                                    parametric,
                                    iterative,
                                    summed,
                                    summed_fnmap,
                                },
                            ))
                        })
                        .collect::<Result<BTreeMap<_, _>>>()
                })
                .collect::<Result<Vec<_>>>()?;

            validate_amplitude_threshold_archive(
                &graph.graph_name,
                graph.threshold_counterterms_are_variants,
                threshold_counterterms.len(),
                &graph.threshold_variants,
                graph.threshold_multipliers.as_ref(),
            )?;
            if let Some(registry) = &graph.metadata_registry {
                validate_amplitude_threshold_metadata_archive(
                    registry,
                    &graph.graph_name,
                    graph.threshold_counterterms_are_variants,
                    threshold_counterterms.len(),
                    &graph.threshold_variants,
                    graph.threshold_multipliers.as_ref(),
                )?;
            } else if graph.threshold_counterterms_are_variants {
                return Err(eyre!(
                    "generalized amplitude graph '{}' is missing its threshold metadata registry",
                    graph.graph_name,
                ));
            }
            let threshold_multipliers = graph
                .threshold_multipliers
                .map(|collection| {
                    build_threshold_multiplier_collection(
                        collection,
                        graph.threshold_variants.len(),
                        0,
                        state_map,
                    )
                })
                .transpose()?;
            if let Some(registry) = &graph.metadata_registry {
                validate_loaded_amplitude_threshold_evaluators(
                    registry,
                    threshold_multipliers.as_ref(),
                    &graph.graph_name,
                )?;
            }

            println!("Loaded evaluators for graph {}", graph.graph_name);
            graph_terms.push(LoadedStandaloneGraphTerm {
                orientations: graph.orientations,
                graph_name: graph.graph_name,
                param_builder_params: params,
                original_integrand,
                threshold_counterterms,
                threshold_counterterms_are_variants: graph.threshold_counterterms_are_variants,
                threshold_variants: graph.threshold_variants,
                threshold_multipliers,
                metadata_registry: graph.metadata_registry,
            });
        }

        Ok(LoadedStandaloneEvaluators { graph_terms })
    }
}

impl StandaloneEvaluatorArchive {
    pub fn load(self) -> Result<LoadedStandaloneEvaluators> {
        if self.version != STANDALONE_EVALUATORS_VERSION {
            return Err(eyre!(
                "Unsupported version {} (expected {})",
                self.version,
                STANDALONE_EVALUATORS_VERSION
            ));
        }

        let mut state_cursor = Cursor::new(&self.symbolica_state);
        let state_map = State::import(&mut state_cursor, None)?;

        self.load_impl(&state_map)
    }
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneGraphTermArchive<A = Vec<u8>> {
    pub(crate) graph_name: String,
    pub(crate) orientations: Vec<Vec<i8>>,
    pub(crate) param_builder_params: Vec<A>,
    pub(crate) fn_map_entries: Vec<SerializedFnMapEntry<A>>,
    pub(crate) original_integrand: StandaloneEvaluatorStackArchive<A>,
    pub(crate) threshold_counterterms: Vec<Vec<StandaloneIndexedEvaluatorStackArchive<A>>>,
    pub(crate) threshold_counterterms_are_variants: bool,
    pub(crate) threshold_variants: Vec<StandaloneAmplitudeThresholdVariant>,
    pub(crate) threshold_multipliers: Option<StandaloneThresholdMultiplierCollectionArchive<A>>,
    pub(crate) metadata_registry: Option<ThresholdCountertermMetadataRegistry>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneAmplitudeThresholdVariant {
    pub(crate) variant_id: usize,
    pub(crate) name: String,
    pub(crate) raised_esurface_id: usize,
    pub(crate) generated: bool,
    pub(crate) active: bool,
    pub(crate) requested_subspace: Option<Vec<usize>>,
    pub(crate) requested_parent_lmb: Option<Vec<usize>>,
    pub(crate) resolved_parent_lmb: Vec<usize>,
    pub(crate) resolved_subspace: Vec<usize>,
    pub(crate) subspace_loop_count: usize,
    pub(crate) multiplier_expression: Option<String>,
    pub(crate) multiplier_symmetrize: bool,
    pub(crate) multiplier_opaque_derivatives: bool,
    pub(crate) threshold_edge_sets: Vec<Vec<usize>>,
    pub(crate) explicit_associations: Vec<bool>,
}

fn validate_amplitude_threshold_archive<A: PartialEq>(
    graph_name: &str,
    generalized: bool,
    evaluator_slots: usize,
    variants: &[StandaloneAmplitudeThresholdVariant],
    multipliers: Option<&StandaloneThresholdMultiplierCollectionArchive<A>>,
) -> Result<()> {
    if !generalized {
        if !variants.is_empty() || multipliers.is_some() {
            return Err(eyre!(
                "legacy amplitude graph '{graph_name}' cannot contain generalized threshold variant metadata or multipliers",
            ));
        }
        return Ok(());
    }
    if variants.len() != evaluator_slots {
        return Err(eyre!(
            "amplitude graph '{graph_name}' standalone archive has {evaluator_slots} variant evaluator slots but {} variant metadata entries",
            variants.len(),
        ));
    }

    for (variant_id, variant) in variants.iter().enumerate() {
        if variant.variant_id != variant_id {
            return Err(eyre!(
                "amplitude graph '{graph_name}' standalone threshold variant IDs are not contiguous and ordered",
            ));
        }
        if variant.name.trim().is_empty() {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold variant {variant_id} has a blank name",
            ));
        }
        if variant.subspace_loop_count == 0
            || variant.resolved_subspace.len() != variant.subspace_loop_count
            || variant.resolved_parent_lmb.is_empty()
            || variant
                .resolved_parent_lmb
                .iter()
                .collect::<BTreeSet<_>>()
                .len()
                != variant.resolved_parent_lmb.len()
            || variant
                .resolved_subspace
                .iter()
                .collect::<BTreeSet<_>>()
                .len()
                != variant.resolved_subspace.len()
            || variant
                .resolved_subspace
                .iter()
                .any(|edge| !variant.resolved_parent_lmb.contains(edge))
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold variant {variant_id} has inconsistent resolved parent/subspace metadata",
            ));
        }
        for (label, edges) in [
            ("requested subspace", variant.requested_subspace.as_deref()),
            (
                "requested parent LMB",
                variant.requested_parent_lmb.as_deref(),
            ),
        ] {
            if edges.is_some_and(|edges| {
                edges.is_empty() || edges.iter().collect::<BTreeSet<_>>().len() != edges.len()
            }) {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' threshold variant {variant_id} has an invalid {label}",
                ));
            }
        }
        if variant.threshold_edge_sets.len() != variant.explicit_associations.len()
            || variant.threshold_edge_sets.is_empty()
            || variant
                .threshold_edge_sets
                .iter()
                .any(|edges| edges.is_empty() || edges.windows(2).any(|pair| pair[0] >= pair[1]))
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold variant {variant_id} has invalid association metadata",
            ));
        }
        if variant.multiplier_symmetrize
            || !variant.multiplier_opaque_derivatives
            || variant
                .multiplier_expression
                .as_ref()
                .is_some_and(|expression| expression.trim().is_empty())
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold variant {variant_id} has unsupported multiplier metadata",
            ));
        }
    }

    let has_multiplier = variants
        .iter()
        .any(|variant| variant.multiplier_expression.is_some());
    if has_multiplier != multipliers.is_some() {
        return Err(eyre!(
            "amplitude graph '{graph_name}' threshold multiplier metadata and evaluator payload disagree",
        ));
    }
    if let Some(multipliers) = multipliers
        && (multipliers.left_variants.len() != variants.len()
            || !multipliers.right_variants.is_empty()
            || multipliers
                .left_variants
                .iter()
                .enumerate()
                .any(|(variant_id, reference)| {
                    reference.variant_id != variant_id
                        || reference.evaluator_id.is_some()
                            != variants[variant_id].multiplier_expression.is_some()
                }))
    {
        return Err(eyre!(
            "amplitude graph '{graph_name}' standalone multiplier references disagree with its variant registry",
        ));
    }
    if let Some(multipliers) = multipliers {
        validate_threshold_multiplier_archive(multipliers, variants.len(), 0).with_context(
            || format!("amplitude graph '{graph_name}' has invalid threshold multipliers"),
        )?;
    }
    Ok(())
}

fn validate_amplitude_threshold_metadata_archive<A: PartialEq>(
    registry: &ThresholdCountertermMetadataRegistry,
    graph_name: &str,
    generalized: bool,
    evaluator_slots: usize,
    variants: &[StandaloneAmplitudeThresholdVariant],
    multipliers: Option<&StandaloneThresholdMultiplierCollectionArchive<A>>,
) -> Result<()> {
    registry.validate().with_context(|| {
        format!("amplitude graph '{graph_name}' has an invalid threshold metadata registry")
    })?;
    if registry.graph_name != graph_name {
        return Err(eyre!(
            "amplitude graph '{graph_name}' threshold metadata registry belongs to '{}'",
            registry.graph_name,
        ));
    }

    let expected_variant_count = if generalized {
        variants.len()
    } else {
        evaluator_slots
    };
    if registry.variants.len() != expected_variant_count {
        return Err(eyre!(
            "amplitude graph '{graph_name}' threshold metadata registry has {} variants but its runtime payload has {expected_variant_count}",
            registry.variants.len(),
        ));
    }

    let mut raised_to_esurfaces = BTreeMap::<usize, Vec<usize>>::new();
    let mut esurfaces_to_raised = BTreeMap::<Vec<usize>, usize>::new();
    for (variant_id, metadata) in registry.variants.iter().enumerate() {
        if metadata.cut_group_id.is_some() || metadata.side != ThresholdCountertermSide::Amplitude {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata variant {variant_id} has a cut group or non-amplitude side",
            ));
        }
        if metadata.name.trim().is_empty()
            || metadata.subspace_loop_count == 0
            || metadata.resolved_subspace.len() != metadata.subspace_loop_count
            || metadata.resolved_parent_lmb.is_empty()
            || metadata
                .resolved_parent_lmb
                .iter()
                .collect::<BTreeSet<_>>()
                .len()
                != metadata.resolved_parent_lmb.len()
            || metadata
                .resolved_subspace
                .iter()
                .collect::<BTreeSet<_>>()
                .len()
                != metadata.resolved_subspace.len()
            || metadata
                .resolved_subspace
                .iter()
                .any(|edge| !metadata.resolved_parent_lmb.contains(edge))
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata variant {variant_id} has invalid name or resolved subspace metadata",
            ));
        }
        for (label, edges) in [
            ("requested subspace", metadata.requested_subspace.as_deref()),
            (
                "requested parent LMB",
                metadata.requested_parent_lmb.as_deref(),
            ),
        ] {
            if edges.is_some_and(|edges| {
                edges.is_empty() || edges.iter().collect::<BTreeSet<_>>().len() != edges.len()
            }) {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' threshold metadata variant {variant_id} has an invalid {label}",
                ));
            }
        }
        if metadata.associations.is_empty()
            || metadata.threshold_esurface_ids
                != metadata
                    .associations
                    .iter()
                    .map(|association| association.esurface_id)
                    .collect::<Vec<_>>()
            || metadata
                .threshold_esurface_ids
                .iter()
                .collect::<BTreeSet<_>>()
                .len()
                != metadata.threshold_esurface_ids.len()
            || metadata.associations.iter().any(|association| {
                association.cut_id.is_some()
                    || !association.cut_edges.is_empty()
                    || association.threshold_edges.is_empty()
                    || association
                        .threshold_edges
                        .windows(2)
                        .any(|pair| pair[0] >= pair[1])
                    || !association.eligible
            })
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata variant {variant_id} has invalid amplitude associations",
            ));
        }
        if metadata.multiplier.as_ref().is_some_and(|multiplier| {
            multiplier.expression.trim().is_empty()
                || multiplier.symmetrize
                || !multiplier.opaque_derivatives
        }) {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata variant {variant_id} has unsupported multiplier metadata",
            ));
        }

        if !generalized {
            if metadata.multiplier.is_some() {
                return Err(eyre!(
                    "legacy amplitude graph '{graph_name}' threshold metadata variant {variant_id} cannot have a multiplier",
                ));
            }
            continue;
        }

        let archived = &variants[variant_id];
        if metadata.name != archived.name
            || metadata.requested_subspace != archived.requested_subspace
            || metadata.requested_parent_lmb != archived.requested_parent_lmb
            || metadata.resolved_parent_lmb != archived.resolved_parent_lmb
            || metadata.resolved_subspace != archived.resolved_subspace
            || metadata.subspace_loop_count != archived.subspace_loop_count
            || metadata.generated != archived.generated
            || metadata.active != archived.active
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata variant {variant_id} disagrees with its standalone runtime record",
            ));
        }
        if metadata.associations.len() != archived.threshold_edge_sets.len()
            || metadata
                .associations
                .iter()
                .zip(&archived.threshold_edge_sets)
                .zip(&archived.explicit_associations)
                .any(|((association, threshold_edges), explicit)| {
                    association.threshold_edges != *threshold_edges
                        || (association.origin == ThresholdCountertermOrigin::Explicit) != *explicit
                })
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata variant {variant_id} association provenance disagrees with its standalone runtime record",
            ));
        }
        match (&metadata.multiplier, &archived.multiplier_expression) {
            (None, None) => {}
            (Some(multiplier), Some(expression))
                if multiplier.expression == *expression
                    && multiplier.symmetrize == archived.multiplier_symmetrize
                    && multiplier.opaque_derivatives == archived.multiplier_opaque_derivatives => {}
            _ => {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' threshold metadata variant {variant_id} multiplier disagrees with its standalone runtime record",
                ));
            }
        }

        if let Some(previous) = raised_to_esurfaces.insert(
            archived.raised_esurface_id,
            metadata.threshold_esurface_ids.clone(),
        ) && previous != metadata.threshold_esurface_ids
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' raised E-surface {} is associated with inconsistent threshold E-surface groups",
                archived.raised_esurface_id,
            ));
        }
        if let Some(previous) = esurfaces_to_raised.insert(
            metadata.threshold_esurface_ids.clone(),
            archived.raised_esurface_id,
        ) && previous != archived.raised_esurface_id
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold E-surface group {:?} refers to inconsistent raised E-surfaces",
                metadata.threshold_esurface_ids,
            ));
        }
    }

    let collection_references = multipliers.map(|collection| &collection.left_variants);
    let mut variant_evaluators = vec![None; registry.variants.len()];
    let mut collection_evaluators = BTreeSet::new();
    for (evaluator_id, evaluator) in registry.evaluators.iter().enumerate() {
        if evaluator.cut_group_id.is_some() || evaluator.expression.trim().is_empty() {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata evaluator {evaluator_id} has a cut group or blank expression",
            ));
        }
        if !collection_evaluators.insert(evaluator.collection_evaluator_id) {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata refers to cut-local evaluator {} more than once",
                evaluator.collection_evaluator_id,
            ));
        }
        let collection = multipliers.ok_or_else(|| {
            eyre!(
                "amplitude graph '{graph_name}' threshold metadata evaluator {evaluator_id} has no multiplier archive",
            )
        })?;
        if evaluator.collection_evaluator_id >= collection.evaluators.len() {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata evaluator {evaluator_id} refers to missing cut-local evaluator {}",
                evaluator.collection_evaluator_id,
            ));
        }
        let archived_variant_ids = collection
            .left_variants
            .iter()
            .filter_map(|reference| {
                (reference.evaluator_id == Some(evaluator.collection_evaluator_id))
                    .then_some(reference.variant_id)
            })
            .collect::<Vec<_>>();
        if evaluator.variant_ids != archived_variant_ids {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata evaluator {evaluator_id} variant links disagree with the multiplier archive",
            ));
        }
        for &variant_id in &evaluator.variant_ids {
            let variant = registry.variants.get(variant_id).ok_or_else(|| {
                eyre!(
                    "amplitude graph '{graph_name}' threshold metadata evaluator {evaluator_id} refers to missing variant {variant_id}",
                )
            })?;
            if variant.multiplier.is_none()
                || variant_evaluators[variant_id]
                    .replace(evaluator_id)
                    .is_some()
            {
                return Err(eyre!(
                    "amplitude graph '{graph_name}' threshold metadata variant {variant_id} has an invalid evaluator assignment",
                ));
            }
        }
    }
    if collection_evaluators.len()
        != multipliers.map_or(0, |collection| collection.evaluators.len())
    {
        return Err(eyre!(
            "amplitude graph '{graph_name}' threshold metadata does not cover every multiplier evaluator",
        ));
    }
    for (variant_id, variant) in registry.variants.iter().enumerate() {
        let archived_evaluator = collection_references
            .and_then(|references| references.get(variant_id))
            .and_then(|reference| reference.evaluator_id);
        let metadata_evaluator = variant_evaluators[variant_id]
            .map(|evaluator_id| registry.evaluators[evaluator_id].collection_evaluator_id);
        if archived_evaluator != metadata_evaluator
            || variant.multiplier.is_some() != metadata_evaluator.is_some()
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata variant {variant_id} multiplier/evaluator links disagree",
            ));
        }
    }

    let mut actual_components = BTreeSet::new();
    for (component_id, component) in registry.components.iter().enumerate() {
        if component.cut_group_id.is_some()
            || component.variant_ids.len() != 1
            || component.evaluator_ids.len() != 1
        {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata component {component_id} is not a single-variant amplitude component",
            ));
        }
        let variant_id = component.variant_ids[0];
        if component.evaluator_ids[0] != variant_evaluators.get(variant_id).copied().flatten() {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata component {component_id} has an inconsistent evaluator link",
            ));
        }
        if !actual_components.insert((component.kind, variant_id)) {
            return Err(eyre!(
                "amplitude graph '{graph_name}' has duplicate threshold component metadata for variant {variant_id}",
            ));
        }
    }
    let expected_components = (0..registry.variants.len())
        .flat_map(|variant_id| {
            [
                (ThresholdCountertermComponentKind::Local, variant_id),
                (ThresholdCountertermComponentKind::Integrated, variant_id),
            ]
        })
        .collect::<BTreeSet<_>>();
    if actual_components != expected_components {
        return Err(eyre!(
            "amplitude graph '{graph_name}' threshold component registry does not contain exactly one local and integrated component per variant",
        ));
    }
    Ok(())
}

fn validate_loaded_amplitude_threshold_evaluators(
    registry: &ThresholdCountertermMetadataRegistry,
    multipliers: Option<&LoadedStandaloneThresholdMultiplierCollection>,
    graph_name: &str,
) -> Result<()> {
    for evaluator in &registry.evaluators {
        let archived_expression = multipliers
            .and_then(|collection| collection.evaluators.get(evaluator.collection_evaluator_id))
            .and_then(|loaded| loaded.0.first())
            .ok_or_else(|| {
                eyre!(
                    "amplitude graph '{graph_name}' threshold metadata evaluator {} has no loaded expression",
                    evaluator.evaluator_id,
                )
            })?;
        let metadata_expression = try_parse!(evaluator.expression.as_str()).map_err(|error| {
            eyre!(
                "amplitude graph '{graph_name}' threshold metadata evaluator {} has an invalid expression: {error}",
                evaluator.evaluator_id,
            )
        })?;
        if archived_expression != &metadata_expression {
            return Err(eyre!(
                "amplitude graph '{graph_name}' threshold metadata evaluator {} expression disagrees with the multiplier archive",
                evaluator.evaluator_id,
            ));
        }
    }
    Ok(())
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneIndexedEvaluatorStackArchive<A = Vec<u8>> {
    pub(crate) cut_cff_index: StandaloneCutCFFIndex,
    pub(crate) evaluator_stack: StandaloneEvaluatorStackArchive<A>,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneEvaluatorStackArchive<A = Vec<u8>> {
    pub(crate) single_parametric: StandaloneGenericEvaluatorArchive<A>,
    pub(crate) iterative: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) summed_function_map: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) summed: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) representative_input: Vec<StandaloneComplexInput>,
    pub(crate) start: usize,
    pub(crate) mult_offset: usize,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneGenericEvaluatorArchive<A = Vec<u8>> {
    pub(crate) exprs: Vec<A>,
    pub(crate) additional_fn_map_entries: Vec<SerializedFnMapEntry<A>>,
    pub(crate) dual_shape: Option<Vec<Vec<usize>>>,
}

type SerializedFnMapEntry<A> = (A, A, Vec<A>, Vec<A>);
type ParsedFnMapEntry = (Atom, Atom, Vec<Atom>, Vec<Indeterminate>);
type LoadedGenericEvaluator = (
    Vec<Atom>,
    Vec<Replacement>,
    ExpressionEvaluator<Complex<f64>>,
    Vec<Complex<f64>>,
);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum StandaloneBackend {
    Eager,
    Cpp,
    Assembly,
    Symjit,
}

impl StandaloneBackend {
    fn parse(value: &str) -> Result<Self> {
        match value {
            "eager" => Ok(Self::Eager),
            "c++" | "cpp" => Ok(Self::Cpp),
            "assembly" => Ok(Self::Assembly),
            "symjit" => Ok(Self::Symjit),
            _ => Err(eyre!(
                "Unsupported backend '{value}', expected eager, c++, assembly, or symjit"
            )),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Eager => "eager",
            Self::Cpp => "c++",
            Self::Assembly => "assembly",
            Self::Symjit => "symjit",
        }
    }

    fn inline_asm(self) -> InlineASM {
        match self {
            Self::Assembly => InlineASM::default(),
            Self::Eager | Self::Cpp | Self::Symjit => InlineASM::None,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum StandaloneMethod {
    SingleParametric,
    Iterative,
    SummedFunctionMap,
    Summed,
}

impl StandaloneMethod {
    fn parse(value: &str) -> Result<Self> {
        match value {
            "single_parametric" | "parametric" => Ok(Self::SingleParametric),
            "iterative" => Ok(Self::Iterative),
            "summed_function_map" | "summed_fnmap" => Ok(Self::SummedFunctionMap),
            "summed" => Ok(Self::Summed),
            _ => Err(eyre!(
                "Unsupported method '{value}', expected single_parametric, iterative, summed_function_map, or summed"
            )),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::SingleParametric => "single_parametric",
            Self::Iterative => "iterative",
            Self::SummedFunctionMap => "summed_function_map",
            Self::Summed => "summed",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum StandaloneStackSelection {
    Original,
    ThresholdCounterterm((usize, usize)),
}

impl StandaloneStackSelection {
    fn parse(value: &str) -> Result<Self> {
        if value == "original" {
            return Ok(Self::Original);
        }

        if let Some(rest) = value.strip_prefix("ct:") {
            let mut parts = rest.split(',');
            let first = parts
                .next()
                .ok_or_else(|| eyre!("Invalid threshold counterterm format"))?
                .parse::<usize>()?;
            let second = parts
                .next()
                .ok_or_else(|| eyre!("Invalid threshold counterterm format"))?
                .parse::<usize>()?;
            return Ok(Self::ThresholdCounterterm((first, second)));
        }

        Err(eyre!(
            "Unsupported stack '{value}', expected original or ct:<first>,<second>"
        ))
    }

    fn label(&self) -> String {
        match self {
            Self::Original => "original".to_string(),
            Self::ThresholdCounterterm((first, second)) => format!("ct_{first}_{second}"),
        }
    }
}

#[derive(Debug)]
struct StandaloneCliOptions {
    input: PathBuf,
    input_json: Option<PathBuf>,
    graph_index: usize,
    graph_name: Option<String>,
    stack: StandaloneStackSelection,
    method: StandaloneMethod,
    orientation_index: Option<usize>,
    backend: StandaloneBackend,
    compare_backends: Vec<StandaloneBackend>,
    artifact_dir: Option<PathBuf>,
    print_input: bool,
}

impl Default for StandaloneCliOptions {
    fn default() -> Self {
        Self {
            input: PathBuf::from("standalone_evaluators.bin"),
            input_json: None,
            graph_index: 0,
            graph_name: None,
            stack: StandaloneStackSelection::Original,
            method: StandaloneMethod::SingleParametric,
            orientation_index: None,
            backend: StandaloneBackend::Eager,
            compare_backends: Vec::new(),
            artifact_dir: None,
            print_input: false,
        }
    }
}

enum StandaloneRuntimeEvaluator<'a> {
    Eager(&'a mut ExpressionEvaluator<Complex<f64>>),
    Compiled(CompiledComplexEvaluator),
    Symjit(Box<JITCompiledEvaluator<Complex<f64>>>),
}

impl<'a> StandaloneRuntimeEvaluator<'a> {
    fn build(
        evaluator: &'a mut ExpressionEvaluator<Complex<f64>>,
        backend: StandaloneBackend,
        artifact_root: &Path,
        label: &str,
    ) -> Result<Self> {
        match backend {
            StandaloneBackend::Eager => Ok(Self::Eager(evaluator)),
            StandaloneBackend::Cpp | StandaloneBackend::Assembly => {
                fs::create_dir_all(artifact_root)?;
                let function_name = format!(
                    "standalone_{}_{}",
                    sanitize_label(label),
                    sanitize_label(backend.as_str())
                );
                let source_path = artifact_root.join(format!("{function_name}.cpp"));
                let library_path = artifact_root.join(format!("{function_name}.so"));
                let compiled = evaluator
                    .export_cpp::<Complex<f64>>(
                        &source_path,
                        &function_name,
                        ExportSettings::new()
                            .include_header(true)
                            .inline_asm(backend.inline_asm())
                            .custom_header(None),
                    )
                    .map_err(|err| eyre!(err))?
                    .compile(&library_path, CompileOptions::default())
                    .map_err(|err| eyre!(err))?
                    .load()
                    .map_err(|err| eyre!(err))?;
                Ok(Self::Compiled(compiled))
            }
            StandaloneBackend::Symjit => Ok(Self::Symjit(Box::new(
                evaluator
                    // SymJIT 2.21 supports optimization levels up to O2 and cannot compact some
                    // complex temporary layouts.
                    .jit_compile(
                        JITCompilationSettings::new()
                            .optimization_level(2)
                            .with_option("compact", "false"),
                    )
                    .map_err(|err| eyre!(err))?,
            ))),
        }
    }

    fn evaluate(&mut self, args: &[Complex<f64>], out: &mut [Complex<f64>]) {
        match self {
            Self::Eager(evaluator) => evaluator.evaluate(args, out),
            Self::Compiled(evaluator) => evaluator.evaluate(args, out),
            Self::Symjit(evaluator) => evaluator.evaluate(args, out),
        }
    }
}

pub trait ImportWithMap {
    fn import_with_map(&self, state_map: &StateMap) -> Result<Atom>;
}

impl ImportWithMap for Vec<u8> {
    fn import_with_map(&self, state_map: &StateMap) -> Result<Atom> {
        let mut cursor = Cursor::new(self);
        Atom::import_with_map(&mut cursor, state_map).map_err(|e| eyre!(e))
    }
}

impl ImportWithMap for String {
    fn import_with_map(&self, _: &StateMap) -> Result<Atom> {
        try_parse!(self).map_err(|e| eyre!(e))
    }
}

fn parse_fn_map_entries<A: ImportWithMap>(
    entries: &[SerializedFnMapEntry<A>],
    state_map: &StateMap,
) -> Result<Vec<ParsedFnMapEntry>> {
    entries
        .iter()
        .map(|(lhs, rhs, tags, args)| {
            let lhs_atom = lhs.import_with_map(state_map)?;
            let rhs_atom = rhs.import_with_map(state_map)?;
            let tags = tags
                .iter()
                .map(|t| t.import_with_map(state_map))
                .collect::<Result<Vec<_>>>()?;
            let args = args
                .iter()
                .map(|a| {
                    let a = a.import_with_map(state_map)?;

                    if let Ok(s) = a.clone().try_into() {
                        Ok(s)
                    } else {
                        Err(eyre!(
                            "Expected indeterminate in function argument, got {}",
                            a
                        ))
                    }
                })
                .collect::<Result<Vec<_>>>()?;

            Ok((lhs_atom, rhs_atom, tags, args))
        })
        .collect()
}

fn apply_fn_map_entries(
    parsed_entries: Vec<ParsedFnMapEntry>,
) -> Result<(Vec<Replacement>, Vec<Replacement>, FunctionMap)> {
    let mut all_replacements: Vec<Replacement> = vec![];
    let mut fn_map: FunctionMap = FunctionMap::new();
    let mut replacements: Vec<Replacement> = vec![];
    fn_map
        .add_aliases([(parse_lit!(gammalooprs::x), Atom::Zero)])
        .map_err(|e| eyre!(e))?;
    for (lhs, rhs, tags, args) in parsed_entries {
        if let AtomView::Var(_) = lhs.as_view() {
            if let Ok(t) = Complex::<Rational>::try_from(rhs.as_view()) {
                fn_map
                    .add_aliases([(lhs.clone(), Atom::num(t))])
                    .map_err(|e| eyre!(e))?;

                all_replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
            } else {
                replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
            }
        } else if let AtomView::Fun(f) = lhs.as_view() {
            if tags.is_empty() {
                let mut wildcards = Vec::new();
                for (i, a) in args.iter().enumerate() {
                    let atom: Atom = a.clone().into();
                    wildcards.push(
                        Replacement::new(atom.to_pattern(), Atom::var(symbol!(format!("x{i}_"))))
                            .allow_new_wildcards_on_rhs(true),
                    )
                }

                fn_map
                    .add_function(f.get_symbol(), args, rhs.clone())
                    .map_err(|e| eyre!(e))?;

                all_replacements.push(Replacement::new(
                    lhs.replace_multiple(&wildcards).to_pattern(),
                    rhs.replace_multiple(&wildcards),
                ));
            } else {
                fn_map
                    .add_tagged_function(f.get_symbol(), tags, args, rhs.clone())
                    .map_err(|e| eyre!(e))?;
                all_replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
            }
        } else {
            all_replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
            replacements.push(Replacement::new(lhs.to_pattern(), rhs.clone()));
        }
    }

    Ok((replacements, all_replacements, fn_map))
}

#[allow(clippy::type_complexity)]
fn build_evaluator<A: ImportWithMap>(
    payload: StandaloneGenericEvaluatorArchive<A>,
    params: &[Atom],
    mut fn_map_entries: Vec<ParsedFnMapEntry>,
    state_map: &StateMap,
    iterate: bool,
) -> Result<LoadedGenericEvaluator> {
    let optimization_settings = OptimizationSettings::new()
        .horner_iterations(10)
        .cores(10)
        .abort_check(Some(Box::new(
            crate::is_interrupt_requested as fn() -> bool,
        )));
    let cpe_iterations = None;
    let exprs = payload
        .exprs
        .iter()
        .map(|b| b.import_with_map(state_map))
        .collect::<Result<Vec<_>>>()?;

    let additional_reps = parse_fn_map_entries(&payload.additional_fn_map_entries, state_map)?;
    fn_map_entries.extend(additional_reps);

    // for (a, _, _, _) in &fn_map_entries {
    //     exprs.push(a.clone());
    // }

    let result = vec![Complex::new(0.0, 0.); exprs.len()];

    let (replacements, all_replacements, fn_map) = apply_fn_map_entries(fn_map_entries)?;
    // for e in exprs.iter() {
    //     println!("Loaded expression: {}", e);
    // }

    if iterate {
        let mut tree: Option<(
            Vec<Atom>,
            ExpressionEvaluator<Complex<Fraction<IntegerRing>>>,
        )> = None;

        for expr in &exprs {
            let eval = expr
                .replace_multiple(&replacements)
                .evaluator(params)
                .function_map(fn_map.clone())
                .optimization_settings(optimization_settings.clone())
                .build()
                .map_err(|e| eyre!("{e} for {expr}:{}", expr.replace_multiple(&replacements)))?;

            tree = Some(if let Some((e, mut t)) = tree {
                t.merge(eval, cpe_iterations).map_err(|e| eyre!(e))?;
                (e, t)
            } else {
                (exprs.clone(), eval)
            });
        }
        tree.map(|(a, eval)| {
            (
                a,
                all_replacements,
                eval.map_coeff(&|r| Complex {
                    re: r.re.to_f64(),
                    im: r.im.to_f64(),
                }),
                result,
            )
        })
        .ok_or_else(|| eyre!("No expressions in evaluator payload"))
    } else {
        let mut replaced_exprs = vec![];
        for expr in &exprs {
            replaced_exprs.push(expr.replace_multiple(&replacements))
        }

        Atom::evaluator_multiple(&replaced_exprs, params)
            .function_map(fn_map)
            .optimization_settings(optimization_settings)
            .build()
            .map(|eval| {
                (
                    exprs,
                    all_replacements,
                    eval.map_coeff(&|r| Complex {
                        re: r.re.to_f64(),
                        im: r.im.to_f64(),
                    }),
                    result,
                )
            })
            .map_err(|e| eyre!("{e}"))
    }
}

pub struct LoadedStandaloneEvaluators {
    pub graph_terms: Vec<LoadedStandaloneGraphTerm>,
}

pub struct LoadedStandaloneGraphTerm {
    pub graph_name: String,
    pub orientations: Vec<Vec<i8>>,
    pub param_builder_params: Vec<Atom>,
    pub original_integrand: LoadedStandaloneEvaluatorStack,
    pub threshold_counterterms:
        Vec<BTreeMap<StandaloneCutCFFIndex, LoadedStandaloneEvaluatorStack>>,
    pub threshold_counterterms_are_variants: bool,
    pub threshold_variants: Vec<StandaloneAmplitudeThresholdVariant>,
    pub threshold_multipliers: Option<LoadedStandaloneThresholdMultiplierCollection>,
    pub metadata_registry: Option<ThresholdCountertermMetadataRegistry>,
}

#[allow(clippy::type_complexity)]
pub struct LoadedStandaloneEvaluatorStack {
    pub(crate) representative_input: Vec<Complex<f64>>,
    pub(crate) orientation_start: usize,
    pub(crate) mult_offset: usize,
    pub parametric: LoadedGenericEvaluator,
    pub iterative: Option<LoadedGenericEvaluator>,
    pub summed: Option<LoadedGenericEvaluator>,
    pub summed_fnmap: Option<LoadedGenericEvaluator>,
}

pub(crate) fn set_orientation_values_impl<A: Clone + Neg<Output = A>>(
    values: &mut [A],
    one: A,
    zero: A,
    mult_offset: usize,
    start: usize,
    orientation: &[i8],
) {
    let minusone = -(one.clone());
    let mut o_start = start * mult_offset;

    for i in orientation {
        match i {
            1 => {
                values[o_start] = one.clone();
                o_start += mult_offset;
            }
            -1 => {
                values[o_start] = minusone.clone();
                o_start += mult_offset;
            }
            0 => {
                values[o_start] = zero.clone();
                o_start += mult_offset;
            }
            _ => panic!("Should be -1,0,1"),
        }
    }
}

impl LoadedStandaloneEvaluators {
    fn select_graph_term_mut(
        &mut self,
        graph_index: usize,
        graph_name: Option<&str>,
    ) -> Result<&mut LoadedStandaloneGraphTerm> {
        if let Some(graph_name) = graph_name {
            return self
                .graph_terms
                .iter_mut()
                .find(|term| term.graph_name == graph_name)
                .ok_or_else(|| eyre!("Unknown graph '{graph_name}'"));
        }

        let graph_count = self.graph_terms.len();
        self.graph_terms.get_mut(graph_index).ok_or_else(|| {
            eyre!(
                "Graph index {} is out of range for {} graph terms",
                graph_index,
                graph_count
            )
        })
    }
}

impl LoadedStandaloneGraphTerm {
    fn stack_mut(
        &mut self,
        selection: &StandaloneStackSelection,
    ) -> Result<&mut LoadedStandaloneEvaluatorStack> {
        match selection {
            StandaloneStackSelection::Original => Ok(&mut self.original_integrand),
            StandaloneStackSelection::ThresholdCounterterm((ct_idx, order_idx)) => {
                let graph_name = self.graph_name.clone();
                let Some(orders) = self.threshold_counterterms.get_mut(*ct_idx) else {
                    return Err(eyre!(
                        "Threshold counterterm index {},{} is out of range for graph {}",
                        ct_idx,
                        order_idx,
                        graph_name
                    ));
                };

                let order_count = orders.len();
                orders.values_mut().nth(*order_idx).ok_or_else(|| {
                    eyre!(
                        "Threshold counterterm index {},{} is out of range for graph {} ({} keyed stacks)",
                        ct_idx,
                        order_idx,
                        graph_name,
                        order_count
                    )
                })
            }
        }
    }
}

fn parse_backend_list(value: &str) -> Result<Vec<StandaloneBackend>> {
    value
        .split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(StandaloneBackend::parse)
        .collect()
}

fn sanitize_label(value: &str) -> String {
    value
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() {
                ch.to_ascii_lowercase()
            } else {
                '_'
            }
        })
        .collect()
}

fn print_usage(program: &str) {
    eprintln!(
        "Usage: {program} [standalone_evaluators.bin|json] [options]\n\
         \n\
         Options:\n\
           --backend <eager|c++|assembly|symjit>\n\
           --compare-backends <backend[,backend...]>\n\
           --input-json <path>\n\
           --graph-index <usize>\n\
           --graph-name <name>\n\
           --stack <original|ct:N,M>\n\
           --method <single_parametric|iterative|summed_function_map|summed>\n\
           --orientation-index <usize> (single_parametric only)\n\
           --artifact-dir <path>\n\
           --print-input\n\
           --help"
    );
}

fn parse_cli_options() -> Result<StandaloneCliOptions> {
    let mut options = StandaloneCliOptions::default();
    let mut args = std::env::args().skip(1);

    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--help" | "-h" => {
                let program = std::env::args()
                    .next()
                    .unwrap_or_else(|| "standalone_evaluators_rust.rs".to_string());
                print_usage(&program);
                std::process::exit(0);
            }
            "--backend" => {
                options.backend = StandaloneBackend::parse(
                    &args
                        .next()
                        .ok_or_else(|| eyre!("Missing value for --backend"))?,
                )?;
            }
            "--compare-backends" => {
                options.compare_backends = parse_backend_list(
                    &args
                        .next()
                        .ok_or_else(|| eyre!("Missing value for --compare-backends"))?,
                )?;
            }
            "--input-json" => {
                options.input_json = Some(PathBuf::from(
                    args.next()
                        .ok_or_else(|| eyre!("Missing value for --input-json"))?,
                ));
            }
            "--graph-index" => {
                options.graph_index = args
                    .next()
                    .ok_or_else(|| eyre!("Missing value for --graph-index"))?
                    .parse()?;
            }
            "--graph-name" => {
                options.graph_name = Some(
                    args.next()
                        .ok_or_else(|| eyre!("Missing value for --graph-name"))?,
                );
            }
            "--stack" => {
                options.stack = StandaloneStackSelection::parse(
                    &args
                        .next()
                        .ok_or_else(|| eyre!("Missing value for --stack"))?,
                )?;
            }
            "--method" => {
                options.method = StandaloneMethod::parse(
                    &args
                        .next()
                        .ok_or_else(|| eyre!("Missing value for --method"))?,
                )?;
            }
            "--orientation-index" => {
                options.orientation_index = Some(
                    args.next()
                        .ok_or_else(|| eyre!("Missing value for --orientation-index"))?
                        .parse()?,
                );
            }
            "--artifact-dir" => {
                options.artifact_dir = Some(PathBuf::from(
                    args.next()
                        .ok_or_else(|| eyre!("Missing value for --artifact-dir"))?,
                ));
            }
            "--print-input" => {
                options.print_input = true;
            }
            _ if arg.starts_with("--") => {
                return Err(eyre!("Unsupported option '{arg}'"));
            }
            _ => {
                options.input = PathBuf::from(arg);
            }
        }
    }

    if options.orientation_index.is_some() && options.method != StandaloneMethod::SingleParametric {
        return Err(eyre!(
            "--orientation-index can only be used with --method single_parametric"
        ));
    }

    Ok(options)
}

fn diff_ratio(lhs: Complex<f64>, rhs: Complex<f64>) -> Option<Complex<f64>> {
    if rhs.re == 0.0 && rhs.im == 0.0 {
        None
    } else {
        Some(lhs / rhs)
    }
}

struct StandaloneEvaluationRequest<'a> {
    backend: StandaloneBackend,
    method: StandaloneMethod,
    orientations: &'a [Vec<i8>],
    orientation_index: Option<usize>,
    custom_input: Option<&'a [Complex<f64>]>,
    artifact_root: &'a Path,
    label: &'a str,
}

impl LoadedStandaloneEvaluatorStack {
    fn selected_evaluator_mut(
        &mut self,
        method: StandaloneMethod,
    ) -> Result<&mut LoadedGenericEvaluator> {
        match method {
            StandaloneMethod::SingleParametric => Ok(&mut self.parametric),
            StandaloneMethod::Iterative => self
                .iterative
                .as_mut()
                .ok_or_else(|| eyre!("Missing iterative evaluator in standalone archive")),
            StandaloneMethod::SummedFunctionMap => self.summed_fnmap.as_mut().ok_or_else(|| {
                eyre!("Missing summed_function_map evaluator in standalone archive")
            }),
            StandaloneMethod::Summed => self
                .summed
                .as_mut()
                .ok_or_else(|| eyre!("Missing summed evaluator in standalone archive")),
        }
    }

    fn evaluate_with_backend(
        &mut self,
        request: StandaloneEvaluationRequest<'_>,
    ) -> Result<Vec<Complex<f64>>> {
        let inputs = if let Some(custom_input) = request.custom_input {
            if custom_input.len() != self.representative_input.len() {
                return Err(eyre!(
                    "Custom input length {} does not match evaluator input length {}",
                    custom_input.len(),
                    self.representative_input.len()
                ));
            }
            vec![custom_input.to_vec()]
        } else {
            match request.method {
                StandaloneMethod::SingleParametric => {
                    if let Some(index) = request.orientation_index {
                        let orientation = request.orientations.get(index).ok_or_else(|| {
                            eyre!(
                                "Orientation index {} is out of range for {} orientations",
                                index,
                                request.orientations.len()
                            )
                        })?;
                        vec![self.set_orientation(orientation)]
                    } else {
                        request
                            .orientations
                            .iter()
                            .map(|orientation| self.set_orientation(orientation))
                            .collect()
                    }
                }
                StandaloneMethod::Iterative
                | StandaloneMethod::SummedFunctionMap
                | StandaloneMethod::Summed => vec![self.representative_input.clone()],
            }
        };

        let evaluator = self.selected_evaluator_mut(request.method)?;
        let (_, _, eval, result_template) = evaluator;
        let mut runtime = StandaloneRuntimeEvaluator::build(
            eval,
            request.backend,
            request.artifact_root,
            &format!("{}_{}", request.label, request.method.as_str()),
        )?;
        let mut accumulated = vec![Complex::new(0.0, 0.0); result_template.len()];

        for input in inputs {
            let mut current = vec![Complex::new(0.0, 0.0); result_template.len()];
            runtime.evaluate(&input, &mut current);
            for (accumulated_value, current_value) in accumulated.iter_mut().zip(current) {
                *accumulated_value += current_value;
            }
        }

        Ok(accumulated)
    }

    // fn benchmark_parametric(&self){
    //     self.parametric.evaluate_single(params)
    // }
    //
    fn benchmark_summed<R: Rng + ?Sized>(
        &mut self,
        rng: &mut R,
        n_samples: usize,
        compile: bool,
    ) -> Option<(Duration, Duration)> {
        let samples: Vec<_> = (0..n_samples).map(|_| self.scramble_input(rng)).collect();
        let Some((_e, _r, eval, result)) = &mut self.summed else {
            return None;
        };

        let mut sum = Duration::ZERO;
        let mut max = Duration::ZERO;
        if compile {
            let compile_started = Instant::now();
            let e = eval
                .export_cpp::<Complex<f64>>(
                    "bench_summed.cpp",
                    "bench_summed",
                    ExportSettings::new()
                        .include_header(true)
                        .inline_asm(symbolica::evaluate::InlineASM::AArch64)
                        .custom_header(None),
                )
                .unwrap()
                .compile("bench_summed.so", CompileOptions::default())
                .unwrap();
            println!(
                "[timing] benchmark_summed compile took {:?}",
                compile_started.elapsed()
            );

            let mut eval = e.load().unwrap();

            for s in &samples {
                let instant = Instant::now();
                eval.evaluate(s, result);
                let duration = instant.elapsed();
                if max < duration {
                    max = duration;
                }
                sum += duration;
            }
        } else {
            for s in &samples {
                let instant = Instant::now();
                eval.evaluate(s, result);
                let duration = instant.elapsed();
                if max < duration {
                    max = duration;
                }
                sum += duration;
            }
        }
        Some((sum / (n_samples as u32), max))
    }

    fn benchmark_iterative<R: Rng + ?Sized>(
        &mut self,
        rng: &mut R,
        n_samples: usize,
        compile: bool,
    ) -> Option<(Duration, Duration, Complex<f64>)> {
        let samples: Vec<_> = (0..n_samples).map(|_| self.scramble_input(rng)).collect();
        let Some((_e, _r, eval, result)) = &mut self.iterative else {
            return None;
        };

        let mut sum = Duration::ZERO;
        let mut max = Duration::ZERO;
        let mut orientation_sum = Complex::new(0.0, 0.0);
        if compile {
            let compile_started = Instant::now();
            let e = eval
                .export_cpp::<Complex<f64>>(
                    "bench_summed.cpp",
                    "bench_summed",
                    ExportSettings::new()
                        .include_header(true)
                        .inline_asm(symbolica::evaluate::InlineASM::AArch64)
                        .custom_header(None),
                )
                .unwrap()
                .compile("bench_summed.so", CompileOptions::default())
                .unwrap();
            println!(
                "[timing] benchmark_iterative compile took {:?}",
                compile_started.elapsed()
            );

            let mut eval = e.load().unwrap();

            for s in &samples {
                let instant = Instant::now();
                eval.evaluate(s, result);
                orientation_sum = Complex::new(0.0, 0.0);
                for _r in result.iter() {
                    orientation_sum += Complex::new(0.0, 0.0);
                }
                let duration = instant.elapsed();
                if max < duration {
                    max = duration;
                }
                sum += duration;
            }
        } else {
            let mut orientation_sum;
            for s in &samples {
                let instant = Instant::now();
                eval.evaluate(s, result);
                orientation_sum = Complex::new(0.0, 0.0);
                for _r in result.iter() {
                    orientation_sum += Complex::new(0.0, 0.0);
                }
                let duration = instant.elapsed();
                if max < duration {
                    max = duration;
                }
                sum += duration;
            }
        }
        Some((sum / (n_samples as u32), max, orientation_sum))
    }

    #[allow(clippy::type_complexity)]
    fn benchmark_parametric<R: Rng + ?Sized>(
        &mut self,
        orientations: &[Vec<i8>],
        rng: &mut R,
        n_samples: usize,
    ) -> (
        Vec<Vec<Complex<f64>>>,
        Vec<Vec<Complex<f64>>>,
        Duration,
        Duration,
    ) {
        let samples: Vec<_> = (0..n_samples)
            .map(|_| {
                let mut samples = vec![];
                for o in orientations {
                    samples.push(self.scramble_input_with_orientation(o, rng))
                }
                samples
            })
            .collect();
        let (_, _, eval, result) = &mut self.parametric;
        let mut result_per_orientation =
            vec![vec![Complex::new_zero(); result.len()]; orientations.len()];

        let mut sum = Duration::ZERO;
        let mut max = Duration::ZERO;
        for s in &samples {
            for r in result.iter_mut() {
                *r = Complex::new(0.0, 0.0);
            }
            let instant = Instant::now();
            for (i, o) in s.iter().enumerate() {
                eval.evaluate(o, &mut result_per_orientation[i]);

                for (r, a) in result.iter_mut().zip(&result_per_orientation[i]) {
                    *r += a
                }
            }
            let duration = instant.elapsed();
            if max < duration {
                max = duration;
            }
            sum += duration;
        }

        (
            result_per_orientation,
            samples.last().unwrap().clone(),
            sum / (n_samples as u32),
            max,
        )
    }

    fn benchmark_summed_fnmap<R: Rng + ?Sized>(
        &mut self,

        rng: &mut R,
        n_samples: usize,
    ) -> Option<(Vec<Complex<f64>>, Duration, Duration)> {
        let samples: Vec<_> = (0..n_samples).map(|_| self.scramble_input(rng)).collect();
        let Some((_e, _r, eval, result)) = &mut self.summed_fnmap else {
            return None;
        };

        // for e in e {
        //     println!("{:120}", e);
        // }

        let mut sum = Duration::ZERO;
        let mut max = Duration::ZERO;

        for s in &samples {
            let instant = Instant::now();
            eval.evaluate(s, result);
            let duration = instant.elapsed();
            if max < duration {
                max = duration;
            }
            sum += duration;
        }

        Some((
            samples.last().unwrap().clone(),
            sum / (n_samples as u32),
            max,
        ))
    }

    fn scramble_input_with_orientation<R: Rng + ?Sized>(
        &self,
        orientation: &[i8],
        _rng: &mut R,
    ) -> Vec<Complex<f64>> {
        let mut new_input = self.representative_input.clone();
        // for n in &mut new_input {
        //     *n = *n * rng.random_range(0.8..1.2);
        // }
        set_orientation_values_impl(
            &mut new_input,
            Complex::new(1., 0.),
            Complex::new(0., 0.),
            self.mult_offset,
            self.orientation_start,
            orientation,
        );
        new_input
    }

    fn scramble_input<R: Rng + ?Sized>(&self, _rng: &mut R) -> Vec<Complex<f64>> {
        self.representative_input.clone()
    }

    fn set_orientation(&self, orientation: &[i8]) -> Vec<Complex<f64>> {
        let mut sample = self.representative_input.clone();
        set_orientation_values_impl(
            &mut sample,
            Complex::new(1., 0.),
            Complex::new(0., 0.),
            self.mult_offset,
            self.orientation_start,
            orientation,
        );
        sample
    }
}

fn load_bin(path: impl AsRef<Path>) -> Result<LoadedStandaloneEvaluators> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let (archive, _): (StandaloneEvaluatorArchive, _) =
        bincode::decode_from_slice(&binary, bincode::config::standard())?;

    archive.load()
}

fn load_json(path: impl AsRef<Path>) -> Result<LoadedStandaloneEvaluators> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let archive: StandaloneEvaluatorArchive<(), String> = serde_json::from_slice(&binary)?;

    archive.load()
}

fn load_custom_input(path: impl AsRef<Path>) -> Result<Vec<Complex<f64>>> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let raw: Vec<[f64; 2]> = serde_json::from_slice(&binary)
        .with_context(|| format!("Failed to parse {}", path.as_ref().display()))?;
    Ok(raw
        .into_iter()
        .map(|[re, im]| Complex::new(re, im))
        .collect())
}

fn print_backend_result(backend: StandaloneBackend, values: &[Complex<f64>]) {
    println!("backend={}", backend.as_str());
    for (index, value) in values.iter().enumerate() {
        println!("  result[{index}] = {value}");
    }
}

fn main() -> Result<()> {
    let options = parse_cli_options()?;
    let input = options.input.clone();

    let Some(ext) = input.extension() else {
        return Err(eyre!("No extension, expected .bin or .json"));
    };
    let mut loaded = match ext.to_string_lossy().as_ref() {
        "bin" => load_bin(&input)?,
        "json" => load_json(&input)?,
        _ => {
            return Err(eyre!(
                "Unsupported file extension {}, expected .bin or .json",
                ext.to_string_lossy()
            ));
        }
    };

    let compare_backends = if options.compare_backends.is_empty() {
        vec![options.backend]
    } else {
        options.compare_backends.clone()
    };
    let custom_input = options
        .input_json
        .as_ref()
        .map(load_custom_input)
        .transpose()?;
    let artifact_root = options.artifact_dir.clone().unwrap_or_else(|| {
        input
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join("standalone_backend_artifacts")
    });

    let graph = loaded.select_graph_term_mut(options.graph_index, options.graph_name.as_deref())?;
    let graph_name = graph.graph_name.clone();
    let orientations = graph.orientations.clone();
    let stack_label = options.stack.label();
    let stack = graph.stack_mut(&options.stack)?;

    println!(
        "graph={} stack={} method={} orientation={} artifact_dir={}",
        graph_name,
        stack_label,
        options.method.as_str(),
        options
            .orientation_index
            .map(|index| index.to_string())
            .unwrap_or_else(|| "all".to_string()),
        artifact_root.display()
    );

    if options.print_input {
        match options.method {
            StandaloneMethod::SingleParametric => {
                if let Some(custom_input) = custom_input.as_ref() {
                    println!("input={custom_input:?}");
                } else if let Some(index) = options.orientation_index {
                    let orientation = orientations.get(index).ok_or_else(|| {
                        eyre!(
                            "Orientation index {} is out of range for {} orientations",
                            index,
                            orientations.len()
                        )
                    })?;
                    println!("input={:?}", stack.set_orientation(orientation));
                } else {
                    for (index, orientation) in orientations.iter().enumerate() {
                        println!("orientation[{index}]={orientation:?}");
                        println!("input[{index}]={:?}", stack.set_orientation(orientation));
                    }
                }
            }
            StandaloneMethod::Iterative
            | StandaloneMethod::SummedFunctionMap
            | StandaloneMethod::Summed => {
                println!(
                    "input={:?}",
                    custom_input
                        .as_ref()
                        .cloned()
                        .unwrap_or_else(|| stack.representative_input.clone())
                );
            }
        }
    }

    let mut results = Vec::new();
    for backend in compare_backends {
        let values = stack.evaluate_with_backend(StandaloneEvaluationRequest {
            backend,
            method: options.method,
            orientations: &orientations,
            orientation_index: options.orientation_index,
            custom_input: custom_input.as_deref(),
            artifact_root: &artifact_root.join(sanitize_label(&graph_name)),
            label: &format!("{}_{}", graph_name, stack_label),
        })?;
        print_backend_result(backend, &values);
        results.push((backend, values));
    }

    if results.len() == 2 {
        let (lhs_backend, lhs_values) = &results[0];
        let (rhs_backend, rhs_values) = &results[1];
        println!(
            "comparison {} -> {}",
            lhs_backend.as_str(),
            rhs_backend.as_str()
        );
        for (index, (lhs_value, rhs_value)) in lhs_values.iter().zip(rhs_values).enumerate() {
            println!("  diff[{index}] = {}", *rhs_value - *lhs_value);
            match diff_ratio(*rhs_value, *lhs_value) {
                Some(ratio) => println!("  ratio[{index}] = {ratio}"),
                None => println!("  ratio[{index}] = undefined"),
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod threshold_variant_archive_tests {
    use super::*;
    use crate::{
        initialisation::test_initialise,
        integrands::process::{
            cross_section::{
                export::export_threshold_multiplier_collection,
                load::{
                    StandaloneGenericEvaluatorArchive as CrossSectionGenericEvaluatorArchive,
                    StandaloneThresholdMultiplierLayoutArchive,
                    StandaloneThresholdMultiplierVariantReference,
                },
            },
            threshold_multiplier::{
                ThresholdMultiplierEvaluatorCollection, ThresholdMultiplierLayout,
            },
        },
        processes::{
            EvaluatorSettings, ThresholdCountertermAssociationMetadata,
            ThresholdCountertermComponentMetadata, ThresholdCountertermEvaluatorMetadata,
            ThresholdCountertermMultiplierMetadata, ThresholdCountertermVariantId,
            ThresholdCountertermVariantMetadata,
        },
    };
    use symbolica::state::State;

    fn variant() -> StandaloneAmplitudeThresholdVariant {
        StandaloneAmplitudeThresholdVariant {
            variant_id: 0,
            name: "default".to_string(),
            raised_esurface_id: 0,
            generated: true,
            active: true,
            requested_subspace: None,
            requested_parent_lmb: None,
            resolved_parent_lmb: vec![1, 2],
            resolved_subspace: vec![1],
            subspace_loop_count: 1,
            multiplier_expression: None,
            multiplier_symmetrize: false,
            multiplier_opaque_derivatives: true,
            threshold_edge_sets: vec![vec![3, 4]],
            explicit_associations: vec![true],
        }
    }

    fn identity_registry() -> ThresholdCountertermMetadataRegistry {
        ThresholdCountertermMetadataRegistry {
            graph_name: "graph".to_string(),
            variants: vec![ThresholdCountertermVariantMetadata {
                variant_id: 0,
                name: "default".to_string(),
                cut_group_id: None,
                associations: vec![ThresholdCountertermAssociationMetadata {
                    cut_id: None,
                    cut_edges: Vec::new(),
                    threshold_edges: vec![3, 4],
                    esurface_id: 7,
                    eligible: true,
                    origin: ThresholdCountertermOrigin::Explicit,
                }],
                side: ThresholdCountertermSide::Amplitude,
                threshold_esurface_ids: vec![7],
                requested_subspace: None,
                resolved_subspace: vec![1],
                requested_parent_lmb: None,
                resolved_parent_lmb: vec![1, 2],
                subspace_loop_count: 1,
                multiplier: None,
                generated: true,
                active: true,
            }],
            evaluators: Vec::new(),
            components: vec![
                ThresholdCountertermComponentMetadata {
                    component_id: 0,
                    cut_group_id: None,
                    kind: ThresholdCountertermComponentKind::Local,
                    variant_ids: vec![0],
                    evaluator_ids: vec![None],
                    sign: -1,
                },
                ThresholdCountertermComponentMetadata {
                    component_id: 1,
                    cut_group_id: None,
                    kind: ThresholdCountertermComponentKind::Integrated,
                    variant_ids: vec![0],
                    evaluator_ids: vec![None],
                    sign: -1,
                },
            ],
        }
    }

    fn multiplier_registry(expression: &str) -> ThresholdCountertermMetadataRegistry {
        let mut registry = identity_registry();
        registry.variants[0].multiplier = Some(ThresholdCountertermMultiplierMetadata {
            expression: expression.to_string(),
            symmetrize: false,
            opaque_derivatives: true,
        });
        registry.evaluators = vec![ThresholdCountertermEvaluatorMetadata {
            evaluator_id: 0,
            cut_group_id: None,
            collection_evaluator_id: 0,
            expression: expression.to_string(),
            variant_ids: vec![0],
        }];
        for component in &mut registry.components {
            component.evaluator_ids = vec![Some(0)];
        }
        registry
    }

    fn assert_metadata_rejected(
        mutate: impl FnOnce(&mut ThresholdCountertermMetadataRegistry),
        expected: &str,
    ) {
        let mut registry = identity_registry();
        mutate(&mut registry);
        let error = validate_amplitude_threshold_metadata_archive::<String>(
            &registry,
            "graph",
            true,
            1,
            &[variant()],
            None,
        )
        .unwrap_err();
        assert!(
            format!("{error:?}").contains(expected),
            "unexpected validation error: {error:?}",
        );
    }

    #[test]
    fn identity_variant_archive_needs_no_multiplier_payload() {
        validate_amplitude_threshold_archive::<String>("graph", true, 1, &[variant()], None)
            .unwrap();
        validate_amplitude_threshold_metadata_archive::<String>(
            &identity_registry(),
            "graph",
            true,
            1,
            &[variant()],
            None,
        )
        .unwrap();
    }

    #[test]
    fn normalized_legacy_metadata_uses_the_runtime_evaluator_dimensions() {
        validate_amplitude_threshold_metadata_archive::<String>(
            &identity_registry(),
            "graph",
            false,
            1,
            &[],
            None,
        )
        .unwrap();

        let error = validate_amplitude_threshold_metadata_archive::<String>(
            &identity_registry(),
            "graph",
            false,
            2,
            &[],
            None,
        )
        .unwrap_err();
        assert!(format!("{error:?}").contains("runtime payload has 2"));
    }

    #[test]
    fn amplitude_metadata_rejects_tampered_variant_fields() {
        assert_metadata_rejected(
            |registry| registry.variants[0].variant_id = 1,
            "stored at index",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].name = "renamed".to_string(),
            "runtime record",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].cut_group_id = Some(0),
            "cut group",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].side = ThresholdCountertermSide::Left,
            "non-amplitude side",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].associations[0].cut_id = Some(0),
            "invalid amplitude associations",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].associations[0].cut_edges = vec![1],
            "invalid amplitude associations",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].associations[0].threshold_edges = vec![3, 5],
            "association provenance",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].associations[0].eligible = false,
            "invalid amplitude associations",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].associations[0].esurface_id = 8,
            "invalid amplitude associations",
        );
        assert_metadata_rejected(
            |registry| {
                registry.variants[0].associations[0].origin =
                    ThresholdCountertermOrigin::Autogenerated
            },
            "association provenance",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].requested_subspace = Some(vec![2]),
            "runtime record",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].requested_parent_lmb = Some(vec![1, 2]),
            "runtime record",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].resolved_subspace = vec![2],
            "runtime record",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].resolved_parent_lmb = vec![1],
            "runtime record",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].subspace_loop_count = 2,
            "invalid name or resolved subspace",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].generated = false,
            "runtime record",
        );
        assert_metadata_rejected(
            |registry| registry.variants[0].active = false,
            "runtime record",
        );
    }

    #[test]
    fn amplitude_metadata_rejects_tampered_components_and_evaluator_links() {
        assert_metadata_rejected(
            |registry| {
                registry.components.pop();
            },
            "exactly one local and integrated",
        );
        assert_metadata_rejected(
            |registry| registry.components[1].kind = ThresholdCountertermComponentKind::Local,
            "duplicate threshold component",
        );
        assert_metadata_rejected(
            |registry| registry.components[1].component_id = 0,
            "stored at index",
        );
        assert_metadata_rejected(
            |registry| registry.components[0].cut_group_id = Some(0),
            "different cut groups",
        );
        assert_metadata_rejected(
            |registry| registry.components[0].variant_ids = vec![1],
            "missing variant",
        );

        let mut standalone_variant = variant();
        standalone_variant.multiplier_expression = Some("1".to_string());
        let multiplier_archive = StandaloneThresholdMultiplierCollectionArchive {
            layout: StandaloneThresholdMultiplierLayoutArchive {
                model_parameter_count: 0,
                additional_parameters: Vec::new(),
                external_count: 0,
                edges: Vec::new(),
                esurfaces: Vec::new(),
                inputs: Vec::new(),
                parameters: Vec::<String>::new(),
            },
            evaluators: vec![CrossSectionGenericEvaluatorArchive {
                exprs: vec!["1".to_string()],
                additional_fn_map_entries: Vec::new(),
                dual_shape: None,
            }],
            left_variants: vec![StandaloneThresholdMultiplierVariantReference {
                variant_id: 0,
                evaluator_id: Some(0),
            }],
            right_variants: Vec::new(),
        };
        let mut registry = multiplier_registry("1");
        validate_amplitude_threshold_metadata_archive(
            &registry,
            "graph",
            true,
            1,
            &[standalone_variant.clone()],
            Some(&multiplier_archive),
        )
        .unwrap();

        let mut duplicate_id = multiplier_registry("1");
        duplicate_id.evaluators[0].evaluator_id = 1;
        let error = validate_amplitude_threshold_metadata_archive(
            &duplicate_id,
            "graph",
            true,
            1,
            &[standalone_variant.clone()],
            Some(&multiplier_archive),
        )
        .unwrap_err();
        assert!(format!("{error:?}").contains("stored at index"));

        let mut changed_expression = multiplier_registry("2");
        changed_expression.evaluators[0].expression = "1".to_string();
        let error = validate_amplitude_threshold_metadata_archive(
            &changed_expression,
            "graph",
            true,
            1,
            &[standalone_variant.clone()],
            Some(&multiplier_archive),
        )
        .unwrap_err();
        assert!(format!("{error:?}").contains("multiplier disagrees"));

        let mut unsupported_flags = multiplier_registry("1");
        unsupported_flags.variants[0]
            .multiplier
            .as_mut()
            .unwrap()
            .opaque_derivatives = false;
        let error = validate_amplitude_threshold_metadata_archive(
            &unsupported_flags,
            "graph",
            true,
            1,
            &[standalone_variant.clone()],
            Some(&multiplier_archive),
        )
        .unwrap_err();
        assert!(format!("{error:?}").contains("unsupported multiplier"));

        registry.evaluators[0].variant_ids.clear();
        registry.components.clear();
        let error = validate_amplitude_threshold_metadata_archive(
            &registry,
            "graph",
            true,
            1,
            &[standalone_variant.clone()],
            Some(&multiplier_archive),
        )
        .unwrap_err();
        assert!(format!("{error:?}").contains("variant links"));

        let mut registry = multiplier_registry("1");
        registry.components[0].evaluator_ids[0] = None;
        let error = validate_amplitude_threshold_metadata_archive(
            &registry,
            "graph",
            true,
            1,
            &[standalone_variant],
            Some(&multiplier_archive),
        )
        .unwrap_err();
        assert!(format!("{error:?}").contains("omits the evaluator"));
    }

    #[test]
    fn default_name_can_repeat_for_distinct_thresholds() {
        let first = variant();
        let mut second = variant();
        second.variant_id = 1;
        second.threshold_edge_sets = vec![vec![5, 6]];

        validate_amplitude_threshold_archive::<String>("graph", true, 2, &[first, second], None)
            .unwrap();
    }

    #[test]
    fn multiplier_references_must_match_variant_metadata() {
        let mut variant = variant();
        variant.multiplier_expression = Some("1".to_string());
        let collection = StandaloneThresholdMultiplierCollectionArchive {
            layout: StandaloneThresholdMultiplierLayoutArchive {
                model_parameter_count: 0,
                additional_parameters: Vec::new(),
                external_count: 0,
                edges: Vec::new(),
                esurfaces: Vec::new(),
                inputs: Vec::new(),
                parameters: Vec::<String>::new(),
            },
            evaluators: vec![CrossSectionGenericEvaluatorArchive {
                exprs: vec!["1".to_string()],
                additional_fn_map_entries: Vec::new(),
                dual_shape: None,
            }],
            left_variants: vec![StandaloneThresholdMultiplierVariantReference {
                variant_id: 0,
                evaluator_id: None,
            }],
            right_variants: Vec::new(),
        };

        let error =
            validate_amplitude_threshold_archive("graph", true, 1, &[variant], Some(&collection))
                .unwrap_err();
        assert!(error.to_string().contains("disagree"));
    }

    #[test]
    fn standalone_roundtrip_preserves_multiplier_evaluator_and_layout() {
        test_initialise().unwrap();
        let additional_parameters = [
            Atom::var(symbol!("archive_weight_b")),
            Atom::var(symbol!("archive_weight_a")),
        ];
        let layout = ThresholdMultiplierLayout::new(
            Vec::new(),
            additional_parameters.to_vec(),
            0,
            vec![0],
            Vec::new(),
        )
        .unwrap();
        let expression = layout
            .parse_expression("2 * archive_weight_b + archive_weight_a")
            .unwrap();
        let collection = ThresholdMultiplierEvaluatorCollection::build(
            layout,
            vec![(ThresholdCountertermVariantId(0), Some(expression))],
            Vec::new(),
            &EvaluatorSettings::default(),
        )
        .unwrap()
        .unwrap();
        let archive = export_threshold_multiplier_collection::<String>(&collection).unwrap();
        let mut multiplier_variant = variant();
        multiplier_variant.multiplier_expression =
            Some(collection.evaluators()[0].expression().to_string());
        let registry =
            multiplier_registry(multiplier_variant.multiplier_expression.as_deref().unwrap());
        validate_amplitude_threshold_archive(
            "archive_graph",
            true,
            1,
            &[multiplier_variant.clone()],
            Some(&archive),
        )
        .unwrap();
        let mut archive_registry = registry.clone();
        archive_registry.graph_name = "archive_graph".to_string();
        validate_amplitude_threshold_metadata_archive(
            &archive_registry,
            "archive_graph",
            true,
            1,
            &[multiplier_variant.clone()],
            Some(&archive),
        )
        .unwrap();
        let mut renamed = archive.clone();
        renamed.layout.parameters[0] = "renamed_weight".to_owned();
        let error = validate_amplitude_threshold_archive(
            "archive_graph",
            true,
            1,
            &[multiplier_variant],
            Some(&renamed),
        )
        .unwrap_err();
        let error = format!("{error:?}");
        assert!(error.contains("amplitude graph 'archive_graph'"));
        assert!(error.contains("additional-parameter atoms and order"));

        let expected_inputs = archive.layout.inputs.clone();
        let expected_additional_parameters = archive.layout.additional_parameters.clone();
        let encoded = bincode::encode_to_vec(&archive, bincode::config::standard()).unwrap();
        let (decoded, _): (StandaloneThresholdMultiplierCollectionArchive<String>, _) =
            bincode::decode_from_slice(&encoded, bincode::config::standard()).unwrap();
        assert_eq!(
            decoded.layout.additional_parameters,
            expected_additional_parameters,
        );
        let mut state_bytes = Vec::new();
        State::export(&mut state_bytes).unwrap();
        let state_map = State::import(&mut Cursor::new(state_bytes), None).unwrap();
        let mut loaded = build_threshold_multiplier_collection(decoded, 1, 0, &state_map).unwrap();
        validate_loaded_amplitude_threshold_evaluators(
            &archive_registry,
            Some(&loaded),
            "archive_graph",
        )
        .unwrap();
        let mut tampered_registry = archive_registry;
        tampered_registry.evaluators[0].expression = "0".to_string();
        let error = validate_loaded_amplitude_threshold_evaluators(
            &tampered_registry,
            Some(&loaded),
            "archive_graph",
        )
        .unwrap_err();
        assert!(format!("{error:?}").contains("expression disagrees"));

        assert_eq!(loaded.layout.edges, vec![0]);
        assert_eq!(loaded.layout.inputs, expected_inputs);
        assert_eq!(loaded.layout.additional_parameters, additional_parameters);
        assert_eq!(loaded.evaluators.len(), 1);
        assert_eq!(loaded.left_variants.len(), 1);
        assert_eq!(loaded.left_variants[0].variant_id, 0);
        assert_eq!(loaded.left_variants[0].evaluator_id, Some(0));
        assert!(loaded.right_variants.is_empty());

        let mut values = vec![Complex::new(0.0, 0.0); loaded.layout.parameters.len()];
        values[0] = Complex::new(3.5, 0.0);
        values[1] = Complex::new(1.25, 0.0);
        let (_, _, evaluator, result) = &mut loaded.evaluators[0];
        evaluator.evaluate(&values, result);
        assert_eq!(result.as_slice(), &[Complex::new(8.25, 0.0)]);
    }
}
