//#!/usr/bin/env -S rust-script
//! ```cargo
//! [dependencies]
//! bincode = "2"
//! bincode-trait-derive = "0.1.1"
//! eyre = "0.6"
//! rand = "0.9"
//! serde_json = "1"
//! serde = { version = "1.0", features = ["derive"] }
//! symbolica = { version = "=3.0.1", default-features = false, features = ["bincode", "float-mpfr", "integer-gmp", "native_code_generation", "serde"] }
//! # Use the same stock SymJIT version as GammaLoop.
//! symjit = "=2.27.0"
//! ```

#![allow(dead_code)]
use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    fs,
    io::Cursor,
    ops::Neg,
    path::{Path, PathBuf},
    time::{Duration, Instant},
};

use crate::integrands::process::retained_dual::{RetainedFunctionDefinition, build_dual_evaluator};
use bincode_trait_derive::{Decode, Encode};
use eyre::{Context, Result, eyre};
use rand::Rng;
use serde::{Deserialize, Serialize};
use symbolica::{
    domains::rational::Fraction,
    evaluate::{FunctionRegistrationOptions, InliningPolicy, JITCompiledEvaluator},
    prelude::*,
    state::StateMap,
};

use crate::integrands::process::cross_section::load::{
    LoadedStandaloneThresholdMultiplierCollection, StandaloneThresholdMultiplierCollectionArchive,
    build_threshold_multiplier_collection, validate_threshold_multiplier_archive,
};
use crate::integrands::process::standalone::{self, StandaloneParametricResidueRows};
use crate::processes::{
    StandaloneNumericTarget, ThresholdCountertermComponentKind,
    ThresholdCountertermMetadataRegistry, ThresholdCountertermOrigin, ThresholdCountertermSide,
};

// Includes representation-indexed rows and physical auxiliary programs; reject
// archives predating this combined layout and the OnShellEnergy/Esurface split.
pub const STANDALONE_EVALUATORS_VERSION: u32 = 14;
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
    pub(crate) three_d_representations: Vec<String>,
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
        let representations = self.three_d_representations.iter().collect::<BTreeSet<_>>();
        if representations.is_empty()
            || representations.len() != self.three_d_representations.len()
            || representations
                .iter()
                .any(|representation| !matches!(representation.as_str(), "cff" | "ltd"))
        {
            return Err(eyre!(
                "Standalone archive has an invalid generated representation list"
            ));
        }
        let mut graph_terms = Vec::new();

        for graph in self.graph_terms {
            if graph.original_integrand.keys().collect::<BTreeSet<_>>() != representations
                || graph
                    .threshold_counterterms
                    .iter()
                    .any(|slot| slot.keys().collect::<BTreeSet<_>>() != representations)
            {
                return Err(eyre!(
                    "Standalone graph '{}' has evaluator representations inconsistent with its generated list",
                    graph.graph_name
                ));
            }
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
                               iterate: bool,
                               parameters: &[Atom]|
             -> Result<LoadedGenericEvaluator> {
                let started = Instant::now();
                let evaluator = build_evaluator(
                    payload,
                    parameters,
                    replacements.clone(),
                    state_map,
                    iterate,
                )?;
                println!(
                    "[timing] build_evaluator {}::{} took {:?}",
                    graph_name,
                    label,
                    started.elapsed()
                );
                Ok(evaluator)
            };

            let load_stack = |label: &str,
                              stack: StandaloneEvaluatorStackArchive<A>|
             -> Result<LoadedStandaloneEvaluatorStack> {
                standalone::validate_standalone_input_layout(
                    params.len(),
                    Some(stack.representative_input.len()),
                    stack.mult_offset,
                    [stack.start, stack.residue_map_id_start],
                    std::iter::once(&stack.single_parametric)
                        .chain(stack.iterative.as_ref())
                        .chain(stack.summed_function_map.as_ref())
                        .chain(stack.summed.as_ref())
                        .chain(
                            stack
                                .single_parametric_auxiliary
                                .as_ref()
                                .map(|(_, payload)| payload),
                        )
                        .map(|payload| payload.dual_shape.as_deref()),
                )?;
                let parametric_rows = stack
                    .parametric_rows
                    .map(|catalog| -> Result<_> {
                        catalog.validate(
                            stack
                                .single_parametric_auxiliary
                                .as_ref()
                                .map(|(ports, payload)| (ports.len(), payload.exprs.len())),
                        )?;
                        Ok(StandaloneParametricResidueRows {
                            parameters: catalog
                                .parameters
                                .iter()
                                .map(|parameter| {
                                    ImportWithMap::import_with_map(parameter, state_map)
                                })
                                .collect::<Result<Vec<_>>>()?,
                            rows: catalog.rows,
                            input_columns: catalog.input_columns,
                            auxiliary_parameters: catalog
                                .auxiliary_parameters
                                .iter()
                                .map(|parameter| {
                                    ImportWithMap::import_with_map(parameter, state_map)
                                })
                                .collect::<Result<Vec<_>>>()?,
                            auxiliary_values: catalog.auxiliary_values,
                            auxiliary_shared_outputs: catalog.auxiliary_shared_outputs,
                            auxiliary_gathers: catalog.auxiliary_gathers,
                        })
                    })
                    .transpose()?;
                let mut parametric_parameters = params.clone();
                let single_parametric_auxiliary = stack
                    .single_parametric_auxiliary
                    .map(|(parameters, payload)| -> Result<_> {
                        if parametric_rows.is_none()
                            || parameters.is_empty()
                            || payload.dual_shape != stack.single_parametric.dual_shape
                            || payload.dual_shape.as_ref().map_or(1, Vec::len) != stack.mult_offset
                        {
                            return Err(eyre!(
                                "Standalone auxiliary payload has inconsistent ports, rows, or dual shape"
                            ));
                        }
                        let parameters = parameters
                            .iter()
                            .map(|parameter| ImportWithMap::import_with_map(parameter, state_map))
                            .collect::<Result<Vec<_>>>()?;
                        let mut auxiliary_inputs = params.clone();
                        auxiliary_inputs.extend(parametric_rows.as_ref().unwrap().auxiliary_parameters.iter().cloned());
                        let output_count = payload.exprs.len();
                        let evaluator = timed_build(
                            &format!("{label}.single_parametric_auxiliary"),
                            payload,
                            false,
                            &auxiliary_inputs,
                        )?;
                        if evaluator.2.get_input_len() != auxiliary_inputs.len() * stack.mult_offset
                            || evaluator.2.get_output_len() != output_count * stack.mult_offset
                        {
                            return Err(eyre!(
                                "Standalone auxiliary evaluator has inconsistent input or output dimensions"
                            ));
                        }
                        parametric_parameters.extend(parameters.iter().cloned());
                        Ok((parameters, evaluator))
                    })
                    .transpose()?;
                if let Some(catalog) = &parametric_rows {
                    parametric_parameters.extend(
                        catalog
                            .input_columns
                            .iter()
                            .map(|&column| catalog.parameters[column].clone()),
                    );
                }
                Ok(LoadedStandaloneEvaluatorStack {
                    single_parametric_auxiliary,
                    parametric_rows,
                    explicit_orientation_sum_only: stack.explicit_orientation_sum_only,
                    production_orientation_ids: stack.production_orientation_ids,
                    parametric: timed_build(
                        &format!("{label}.parametric"),
                        stack.single_parametric,
                        false,
                        &parametric_parameters,
                    )?,
                    orientation_start: stack.start,
                    residue_map_id_start: stack.residue_map_id_start,
                    mult_offset: stack.mult_offset,
                    representative_input: stack
                        .representative_input
                        .iter()
                        .map(StandaloneComplexInput::to_f64)
                        .collect::<Result<Vec<_>>>()?,
                    iterative: stack
                        .iterative
                        .map(|payload| {
                            timed_build(&format!("{label}.iterative"), payload, true, &params)
                        })
                        .transpose()?,
                    summed: stack
                        .summed
                        .map(|payload| {
                            timed_build(&format!("{label}.summed"), payload, false, &params)
                        })
                        .transpose()?,
                    summed_fnmap: stack
                        .summed_function_map
                        .map(|payload| {
                            timed_build(&format!("{label}.summed_fnmap"), payload, false, &params)
                        })
                        .transpose()?,
                })
            };
            let original_integrand = graph
                .original_integrand
                .into_iter()
                .map(|(representation, stack)| {
                    let loaded = load_stack(&format!("{representation}.original"), stack)?;
                    Ok((representation, loaded))
                })
                .collect::<Result<BTreeMap<_, _>>>()?;
            let threshold_counterterms = graph
                .threshold_counterterms
                .into_iter()
                .enumerate()
                .map(|(ct_idx, by_rep)| {
                    by_rep
                        .into_iter()
                        .map(|(representation, stacks)| {
                            let loaded = stacks
                                .into_iter()
                                .enumerate()
                                .map(|(slot_idx, stack)| {
                                    Ok((
                                        stack.cut_cff_index,
                                        load_stack(
                                            &format!(
                                                "{representation}.ct[{ct_idx}].slot[{slot_idx}]"
                                            ),
                                            stack.evaluator_stack,
                                        )?,
                                    ))
                                })
                                .collect::<Result<BTreeMap<_, _>>>()?;
                            Ok((representation, loaded))
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

        Ok(LoadedStandaloneEvaluators {
            three_d_representations: self.three_d_representations,
            graph_terms,
        })
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
    pub(crate) original_integrand: BTreeMap<String, StandaloneEvaluatorStackArchive<A>>,
    pub(crate) threshold_counterterms:
        Vec<BTreeMap<String, Vec<StandaloneIndexedEvaluatorStackArchive<A>>>>,
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
    pub(crate) explicit_orientation_sum_only: bool,
    pub(crate) production_orientation_ids: Vec<usize>,
    pub(crate) parametric_rows: Option<StandaloneParametricResidueRows<A>>,
    pub(crate) single_parametric_auxiliary: Option<(Vec<A>, StandaloneGenericEvaluatorArchive<A>)>,
    pub(crate) single_parametric: StandaloneGenericEvaluatorArchive<A>,
    pub(crate) iterative: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) summed_function_map: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) summed: Option<StandaloneGenericEvaluatorArchive<A>>,
    pub(crate) representative_input: Vec<StandaloneComplexInput>,
    pub(crate) start: usize,
    pub(crate) residue_map_id_start: usize,
    pub(crate) mult_offset: usize,
}

#[derive(Clone, Encode, Decode, Serialize, Deserialize)]
pub struct StandaloneGenericEvaluatorArchive<A = Vec<u8>> {
    pub(crate) exprs: Vec<A>,
    pub(crate) additional_fn_map_entries: Vec<SerializedFnMapEntry<A>>,
    pub(crate) dual_shape: Option<Vec<Vec<usize>>>,
    #[serde(default)]
    pub(crate) zero_components: Vec<(usize, usize)>,
}

type SerializedFnMapEntry<A> = (A, A, Vec<A>, Vec<A>, InliningPolicy, bool);
type ParsedFnMapEntry = (
    Atom,
    Atom,
    Vec<Atom>,
    Vec<Indeterminate>,
    InliningPolicy,
    bool,
);
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
    representation: Option<String>,
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
            representation: None,
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
                let exported = evaluator
                    .export_cpp::<Complex<f64>>(
                        &source_path,
                        &function_name,
                        ExportSettings::new()
                            .include_header(true)
                            .inline_asm(backend.inline_asm())
                            .custom_header(None),
                    )
                    .map_err(|err| eyre!(err))?;

                if matches!(backend, StandaloneBackend::Cpp) {
                    // Symbolica 3.0.1 wraps complex f64 components in T twice,
                    // producing invalid T(T(re), T(im)) for std::complex<double>.
                    // Correct only exact constants from this exported program,
                    // including retained bodies; keep the stock source beside it
                    // so this upstream export defect remains reproducible.
                    let instructions = evaluator.export_instructions();
                    let mut pending = vec![&instructions];
                    let mut visited = HashSet::new();
                    let mut replacements = BTreeMap::new();
                    while let Some(program) = pending.pop() {
                        for coefficient in &program.constants {
                            if coefficient.im != 0.0 {
                                replacements.insert(
                                    coefficient.export_wrapped(),
                                    format!("T({})", coefficient.export()),
                                );
                            }
                        }
                        for callee in &program.sub_evaluators {
                            if visited.insert(std::sync::Arc::as_ptr(&callee.instructions)) {
                                pending.push(callee.instructions.as_ref());
                            }
                        }
                    }
                    let original = fs::read_to_string(&source_path)?;
                    fs::write(source_path.with_extension("symbolica.cpp"), &original)?;
                    if !replacements.is_empty() {
                        let mut corrected = original.clone();
                        for (malformed, valid) in replacements {
                            corrected = corrected.replace(&malformed, &valid);
                        }
                        if corrected != original {
                            fs::write(&source_path, corrected)?;
                        }
                    }
                }
                let compiled = exported
                    .compile(&library_path, CompileOptions::default())
                    .map_err(|err| eyre!(err))?
                    .load()
                    .map_err(|err| eyre!(err))?;
                Ok(Self::Compiled(compiled))
            }
            StandaloneBackend::Symjit => Ok(Self::Symjit(Box::new(
                evaluator
                    // Match the runtime's conservative optimization settings; these
                    // retain the policy introduced for older SymJIT versions.
                    .jit_compile(
                        JITCompilationSettings::new()
                            .direct_translation(true)
                            .optimization_level(2)
                            .with_option("compact", "false")
                            .with_option("cse", "false"),
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
        .map(|(lhs, rhs, tags, args, inlining, is_alias)| {
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

            Ok((lhs_atom, rhs_atom, tags, args, *inlining, *is_alias))
        })
        .collect()
}

fn apply_fn_map_entries(
    parsed_entries: Vec<ParsedFnMapEntry>,
) -> Result<(Vec<Replacement>, Vec<Replacement>, FunctionMap)> {
    let mut fn_map = FunctionMap::new();
    let mut replacements = Vec::new();
    let mut all_replacements = Vec::new();
    fn_map
        .add_aliases([(parse_lit!(gammalooprs::x), Atom::Zero)])
        .map_err(|error| eyre!(error))?;
    // Graph and evaluator archives can share definitions. Policy and alias scope
    // are part of their identity; conflicting definitions remain errors.
    let mut seen = HashSet::new();
    for (lhs, rhs, tags, args, inlining, is_alias) in parsed_entries
        .into_iter()
        .filter(|entry| seen.insert(entry.clone()))
    {
        if is_alias {
            if !args.is_empty() || inlining != InliningPolicy::Always {
                return Err(eyre!(
                    "Caller-scope aliases require no formal arguments and Always inlining"
                ));
            }
            fn_map
                .add_aliases([(lhs.clone(), rhs.clone())])
                .map_err(|error| eyre!(error))?;
            all_replacements.push(Replacement::new(lhs.to_pattern(), rhs));
            continue;
        }
        let options = FunctionRegistrationOptions::new().inlining(inlining);
        let name = match lhs.as_view() {
            AtomView::Fun(function) => Some(function.get_symbol()),
            AtomView::Var(variable) if tags.is_empty() => Some(variable.get_symbol()),
            _ => None,
        };
        if let Some(name) = name {
            fn_map
                .add_tagged_function_with_options(
                    name,
                    tags.clone(),
                    args.clone(),
                    rhs.clone(),
                    options,
                )
                .map_err(|error| eyre!(error))?;
            if inlining == InliningPolicy::Never {
                continue;
            }
            let wildcards = args
                .iter()
                .enumerate()
                .map(|(index, arg)| {
                    let atom = Atom::from(arg.clone());
                    Replacement::new(atom.to_pattern(), Atom::var(symbol!(format!("x{index}_"))))
                        .allow_new_wildcards_on_rhs(true)
                })
                .collect::<Vec<_>>();
            let replacement = Replacement::new(
                lhs.replace_multiple(&wildcards).to_pattern(),
                rhs.replace_multiple(&wildcards),
            );
            all_replacements.push(replacement);
        } else {
            if inlining == InliningPolicy::Never {
                return Err(eyre!(
                    "Non-inlined definitions require a symbol or function call"
                ));
            }
            let replacement = Replacement::new(lhs.to_pattern(), rhs);
            all_replacements.push(replacement.clone());
            replacements.push(replacement);
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
    if params.iter().collect::<HashSet<_>>().len() != params.len() {
        return Err(eyre!("Standalone evaluator parameters must be distinct"));
    }
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

    let result = vec![Complex::new(0.0, 0.); exprs.len()];

    let definitions = fn_map_entries
        .iter()
        .map(
            |(lhs, rhs, tags, args, inlining, is_alias)| RetainedFunctionDefinition {
                lhs: lhs.clone(),
                rhs: rhs.clone(),
                tags: tags.clone(),
                args: args.clone(),
                retained: *inlining == InliningPolicy::Never,
                is_alias: *is_alias,
            },
        )
        .collect::<Vec<_>>();
    let (replacements, all_replacements, fn_map) = apply_fn_map_entries(fn_map_entries)?;
    if let Some(shape) = payload.dual_shape {
        let replaced = exprs
            .iter()
            .map(|expr| expr.replace_multiple(&replacements))
            .collect::<Vec<_>>();
        let evaluator = build_dual_evaluator(
            &replaced,
            params,
            &fn_map,
            &definitions,
            shape,
            payload.zero_components,
            optimization_settings,
        )
        .map_err(|error| eyre!(error))?;
        let result = vec![Complex::new(0.0, 0.0); evaluator.get_output_len()];
        return Ok((
            exprs,
            all_replacements,
            evaluator.map_coeff(&|r| Complex::new(r.re.to_f64(), r.im.to_f64())),
            result,
        ));
    }
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
    pub three_d_representations: Vec<String>,
    pub graph_terms: Vec<LoadedStandaloneGraphTerm>,
}

pub struct LoadedStandaloneGraphTerm {
    pub graph_name: String,
    pub orientations: Vec<Vec<i8>>,
    pub param_builder_params: Vec<Atom>,
    pub original_integrand: BTreeMap<String, LoadedStandaloneEvaluatorStack>,
    pub threshold_counterterms:
        Vec<BTreeMap<String, BTreeMap<StandaloneCutCFFIndex, LoadedStandaloneEvaluatorStack>>>,
    pub threshold_counterterms_are_variants: bool,
    pub threshold_variants: Vec<StandaloneAmplitudeThresholdVariant>,
    pub threshold_multipliers: Option<LoadedStandaloneThresholdMultiplierCollection>,
    pub metadata_registry: Option<ThresholdCountertermMetadataRegistry>,
}

#[allow(clippy::type_complexity)]
pub struct LoadedStandaloneEvaluatorStack {
    pub(crate) parametric_rows: Option<StandaloneParametricResidueRows<Atom>>,
    pub single_parametric_auxiliary: Option<(Vec<Atom>, LoadedGenericEvaluator)>,
    pub(crate) explicit_orientation_sum_only: bool,
    pub(crate) production_orientation_ids: Vec<usize>,
    pub(crate) representative_input: Vec<Complex<f64>>,
    pub(crate) orientation_start: usize,
    pub(crate) residue_map_id_start: usize,
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
        representation: &str,
        selection: &StandaloneStackSelection,
    ) -> Result<&mut LoadedStandaloneEvaluatorStack> {
        match selection {
            StandaloneStackSelection::Original => self
                .original_integrand
                .get_mut(representation)
                .ok_or_else(|| eyre!("Representation {representation} was not generated")),
            StandaloneStackSelection::ThresholdCounterterm((ct_idx, order_idx)) => {
                let graph_name = self.graph_name.clone();
                let Some(orders) = self
                    .threshold_counterterms
                    .get_mut(*ct_idx)
                    .and_then(|reps| reps.get_mut(representation))
                else {
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
           --three-dimensional-representation <cff|ltd>\n\
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
            "--three-dimensional-representation" => {
                options.representation = Some(args.next().ok_or_else(|| {
                    eyre!("Missing value for --three-dimensional-representation")
                })?);
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
        if self.explicit_orientation_sum_only && request.orientation_index.is_some() {
            return Err(eyre!(
                "`--orientation-index` is invalid for an explicit orientation-sum evaluator because its single-parametric expression already contains the complete orientation sum"
            ));
        }

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
                StandaloneMethod::SingleParametric if self.explicit_orientation_sum_only => {
                    vec![self.representative_input.clone()]
                }
                StandaloneMethod::SingleParametric => {
                    if self.production_orientation_ids.len() != request.orientations.len() {
                        return Err(eyre!(
                            "Standalone evaluator has {} production residue-map IDs for {} physical orientations",
                            self.production_orientation_ids.len(),
                            request.orientations.len()
                        ));
                    }
                    if let Some(index) = request.orientation_index {
                        let orientation = request.orientations.get(index).ok_or_else(|| {
                            eyre!(
                                "Orientation index {} is out of range for {} orientations",
                                index,
                                request.orientations.len()
                            )
                        })?;
                        vec![self.set_orientation(index, orientation)?]
                    } else {
                        request
                            .orientations
                            .iter()
                            .enumerate()
                            .map(|(index, orientation)| self.set_orientation(index, orientation))
                            .collect::<Result<Vec<_>>>()?
                    }
                }
                StandaloneMethod::Iterative
                | StandaloneMethod::SummedFunctionMap
                | StandaloneMethod::Summed => vec![self.representative_input.clone()],
            }
        };

        let rows = if request.method == StandaloneMethod::SingleParametric {
            self.parametric_rows.clone()
        } else {
            None
        };
        let multiplicity = self.mult_offset;
        let auxiliary_width = if rows.is_some() {
            self.single_parametric_auxiliary
                .as_ref()
                .map_or(0, |(ports, _)| ports.len() * multiplicity)
        } else {
            0
        };
        let inputs = if request.method == StandaloneMethod::SingleParametric
            && let Some((_, (_, _, evaluator, result_template))) =
                &mut self.single_parametric_auxiliary
        {
            let mut runtime = StandaloneRuntimeEvaluator::build(
                evaluator,
                request.backend,
                request.artifact_root,
                &format!("{}_single_parametric_auxiliary", request.label),
            )?;
            inputs
                .into_iter()
                .map(|physical_input| {
                    let mut input = physical_input.clone();
                    for value in &rows.as_ref().unwrap().auxiliary_values {
                        input.push(Complex::new(value.to_f64(), 0.0));
                        input.extend((1..multiplicity).map(|_| Complex::new(0.0, 0.0)));
                    }
                    let mut values = vec![Complex::new(0.0, 0.0); result_template.len()];
                    runtime.evaluate(&input, &mut values);
                    (physical_input, values)
                })
                .collect::<Vec<_>>()
        } else {
            inputs
                .into_iter()
                .map(|input| (input, Vec::new()))
                .collect()
        };
        let (_, _, eval, result_template) = self.selected_evaluator_mut(request.method)?;
        let mut runtime = StandaloneRuntimeEvaluator::build(
            eval,
            request.backend,
            request.artifact_root,
            &format!("{}_{}", request.label, request.method.as_str()),
        )?;
        let mut accumulated = vec![Complex::new(0.0, 0.0); result_template.len()];
        let mut current = vec![Complex::new(0.0, 0.0); result_template.len()];
        for (physical_input, auxiliary_values) in inputs {
            let physical_len = physical_input.len();
            let mut input = physical_input;
            input.resize(
                physical_len
                    + auxiliary_width
                    + rows
                        .as_ref()
                        .map_or(0, |catalog| catalog.input_columns.len() * multiplicity),
                Complex::new(0.0, 0.0),
            );
            let shared = rows
                .as_ref()
                .map_or(0, |catalog| catalog.auxiliary_shared_outputs * multiplicity);
            input[physical_len..physical_len + shared].copy_from_slice(&auxiliary_values[..shared]);
            for row_index in 0..rows.as_ref().map_or(1, |catalog| catalog.rows.len()) {
                if let Some(catalog) = &rows {
                    for (port, output) in catalog.auxiliary_gathers[row_index].iter().enumerate() {
                        if row_index > 0
                            && *output == catalog.auxiliary_gathers[row_index - 1][port]
                        {
                            continue;
                        }
                        let offset = physical_len + shared + port * multiplicity;
                        let target = &mut input[offset..offset + multiplicity];
                        if let Some(output) = output {
                            target.copy_from_slice(
                                &auxiliary_values
                                    [output * multiplicity..(output + 1) * multiplicity],
                            );
                        } else {
                            target.fill(Complex::new(0.0, 0.0));
                        }
                    }
                    for (port, &column) in catalog.input_columns.iter().enumerate() {
                        input[physical_len + auxiliary_width + port * multiplicity] =
                            Complex::new(catalog.rows[row_index][column].to_f64(), 0.0);
                    }
                }
                runtime.evaluate(&input, &mut current);
                for (sum, value) in accumulated.iter_mut().zip(&current) {
                    *sum += value;
                }
            }
        }

        Ok(accumulated)
    }

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
        orientation_index: usize,
        orientation: &[i8],
        _rng: &mut R,
    ) -> Vec<Complex<f64>> {
        let mut new_input = self.representative_input.clone();
        // for n in &mut new_input {
        //     *n = *n * rng.random_range(0.8..1.2);
        // }
        new_input[self.residue_map_id_start * self.mult_offset] = Complex::new(
            self.production_orientation_ids[orientation_index] as f64,
            0.0,
        );
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

    fn set_orientation(
        &self,
        orientation_index: usize,
        orientation: &[i8],
    ) -> Result<Vec<Complex<f64>>> {
        let mut sample = self.representative_input.clone();
        let production_orientation_id = self
            .production_orientation_ids
            .get(orientation_index)
            .ok_or_else(|| {
                eyre!(
                    "Missing production residue-map ID for runtime orientation {}",
                    orientation_index
                )
            })?;
        sample[self.residue_map_id_start * self.mult_offset] =
            Complex::new(*production_orientation_id as f64, 0.0);
        set_orientation_values_impl(
            &mut sample,
            Complex::new(1., 0.),
            Complex::new(0., 0.),
            self.mult_offset,
            self.orientation_start,
            orientation,
        );
        Ok(sample)
    }
}

fn load_bin(path: impl AsRef<Path>) -> Result<LoadedStandaloneEvaluators> {
    let binary =
        fs::read(&path).with_context(|| format!("Cannot read {}", path.as_ref().display()))?;
    let (version, _): (u32, _) = bincode::decode_from_slice(&binary, bincode::config::standard())?;
    if version != STANDALONE_EVALUATORS_VERSION {
        return Err(eyre!(
            "Unsupported version {version} (expected {}); regenerate the standalone archive",
            STANDALONE_EVALUATORS_VERSION
        ));
    }
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

    let representation = options
        .representation
        .as_ref()
        .or_else(|| loaded.three_d_representations.first())
        .ok_or_else(|| eyre!("No generated 3D representation"))?
        .clone();
    if !loaded.three_d_representations.contains(&representation) {
        return Err(eyre!("Representation {representation} was not generated"));
    }
    let graph = loaded.select_graph_term_mut(options.graph_index, options.graph_name.as_deref())?;
    let graph_name = graph.graph_name.clone();
    let orientations = graph.orientations.clone();
    let stack_label = options.stack.label();
    let stack = graph.stack_mut(&representation, &options.stack)?;

    if stack.explicit_orientation_sum_only && options.orientation_index.is_some() {
        return Err(eyre!(
            "`--orientation-index` is invalid for an explicit orientation-sum evaluator because its single-parametric expression already contains the complete orientation sum"
        ));
    }

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
            StandaloneMethod::SingleParametric if stack.explicit_orientation_sum_only => {
                println!(
                    "input={:?}",
                    custom_input
                        .as_ref()
                        .cloned()
                        .unwrap_or_else(|| stack.representative_input.clone())
                );
            }
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
                    println!("input={:?}", stack.set_orientation(index, orientation)?);
                } else {
                    for (index, orientation) in orientations.iter().enumerate() {
                        println!("orientation[{index}]={orientation:?}");
                        println!(
                            "input[{index}]={:?}",
                            stack.set_orientation(index, orientation)?
                        );
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
            label: &format!("{graph_name}_{representation}_{stack_label}"),
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
                group_id: None,
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
            function_map: Default::default(),
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
                zero_components: Vec::new(),
                parameter_override: None,
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
                zero_components: Vec::new(),
                parameter_override: None,
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
            .parse_expression(
                "scaled(archive_weight_b)",
                &BTreeMap::from([
                    (
                        "scaled(value)".into(),
                        "double(value)+archive_weight_a".into(),
                    ),
                    ("double(value)".into(), "2*value".into()),
                ]),
            )
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

#[cfg(test)]
mod tests {
    use linnet::half_edge::involution::EdgeIndex;
    use symbolica::printer::PrintOptions;

    use crate::utils::GS;

    use super::*;

    #[test]
    fn standalone_loader_cpp_complex_constants_match_eager_and_assembly() -> Result<()> {
        let root = std::env::temp_dir().join(format!(
            "gammaloop_standalone_loader_complex_export_{}",
            std::process::id()
        ));
        let x = symbolica::parse!("standalone_loader_complex_export::x");
        let q = symbolica::parse!("standalone_loader_complex_export::q");
        for imaginary in [true, false] {
            for retained in [false, true] {
                let coefficient = Atom::num(Rational::new(1, 8))
                    * if imaginary {
                        Atom::one() + Atom::i()
                    } else {
                        Atom::one()
                    };
                let offset = Atom::num(Rational::new(1, 16))
                    * if imaginary { Atom::i() } else { Atom::one() };
                let mut functions = FunctionMap::new();
                let product = if retained {
                    let inner = symbol!("standalone_loader_complex_export::inner");
                    let outer = symbol!("standalone_loader_complex_export::outer");
                    for (function, body) in [
                        (inner, &coefficient * &q),
                        (outer, inner.call_args([q.clone()])),
                    ] {
                        functions.add_function_with_options(
                            function,
                            vec![Indeterminate::try_from(q.clone()).map_err(|error| eyre!(error))?],
                            body,
                            FunctionRegistrationOptions::new().inlining(InliningPolicy::Never),
                        )?;
                    }
                    outer.call_args([x.clone()])
                } else {
                    &coefficient * &x
                };
                let expressions = [
                    product + &x * Atom::num(Rational::new(1, 2)) + offset.clone(),
                    offset,
                    Atom::num(Rational::new(3, 8)),
                ];
                let mut trees = expressions.iter().map(|expression| {
                    expression
                        .evaluator(std::slice::from_ref(&x))
                        .function_map(functions.clone())
                        .optimization_settings(OptimizationSettings::new().cores(1))
                        .build()
                        .map_err(|error| eyre!(error))
                });
                let mut tree = trees.next().unwrap()?;
                for next in trees {
                    tree.merge(next?, Some(1)).map_err(|error| eyre!(error))?;
                }
                let mut evaluator =
                    tree.map_coeff(&|value| Complex::new(value.re.to_f64(), value.im.to_f64()));
                let label = format!("complex_{imaginary}_retained_{retained}");
                for backend in [
                    StandaloneBackend::Eager,
                    StandaloneBackend::Cpp,
                    StandaloneBackend::Assembly,
                ] {
                    if !imaginary && matches!(backend, StandaloneBackend::Cpp) {
                        // A real-only regeneration must refresh an existing sidecar.
                        fs::create_dir_all(&root)?;
                        fs::write(
                            root.join(format!("standalone_{label}_c__.symbolica.cpp")),
                            "stale export from an earlier complex program",
                        )?;
                    }
                    let mut runtime =
                        StandaloneRuntimeEvaluator::build(&mut evaluator, backend, &root, &label)?;
                    for input in [2.0, -2.0] {
                        let mut output = [Complex::new(0.0, 0.0); 3];
                        runtime.evaluate(&[Complex::new(input, 0.0)], &mut output);
                        assert_eq!(
                            output,
                            [
                                Complex::new(
                                    input * 0.625 + if imaginary { 0.0 } else { 0.0625 },
                                    if imaginary {
                                        input * 0.125 + 0.0625
                                    } else {
                                        0.0
                                    },
                                ),
                                if imaginary {
                                    Complex::new(0.0, 0.0625)
                                } else {
                                    Complex::new(0.0625, 0.0)
                                },
                                Complex::new(0.375, 0.0),
                            ]
                        );
                    }
                    if matches!(backend, StandaloneBackend::Cpp) {
                        let path = root.join(format!("standalone_{label}_c__.cpp"));
                        let corrected = fs::read_to_string(&path)?;
                        // The unrelated real constructor must survive byte for byte.
                        assert!(corrected.contains("T(3.75e-1)"));
                        let original_path = path.with_extension("symbolica.cpp");
                        if imaginary {
                            let original = fs::read_to_string(original_path)?;
                            assert!(original.contains("T(T("));
                            // Exact spellings expected from these three coefficients;
                            // no other source content may change.
                            let expected = original
                                .replace("T(T(1.25e-1), T(1.25e-1))", "T(1.25e-1, 1.25e-1)")
                                .replace("T(T(6.25e-1), T(1.25e-1))", "T(6.25e-1, 1.25e-1)")
                                .replace("T(T(0e0), T(6.25e-2))", "T(0e0, 6.25e-2)");
                            assert_eq!(corrected, expected);
                        } else {
                            assert_eq!(corrected, fs::read_to_string(original_path)?);
                        }
                    }
                }
            }
        }
        fs::remove_dir_all(root)?;
        Ok(())
    }

    #[test]
    fn standalone_function_policy_roundtrip_preserves_shared_bodies_and_alias_scope() -> Result<()>
    {
        let x = symbolica::parse!("standalone_policy::x");
        let global = symbolica::parse!("standalone_policy::global");
        let alias = symbolica::parse!("standalone_policy::alias(0)");
        let function = symbol!("standalone_policy::function");
        let print = |atom: &Atom| {
            atom.printer(symbolica::printer::PrintOptions::file())
                .to_string()
        };
        for tags in [Vec::new(), vec![Atom::num(7)]] {
            let call = |argument: Atom| {
                FunctionBuilder::new(function)
                    .add_args(&tags)
                    .add_arg(argument)
                    .finish()
            };
            for inlining in [InliningPolicy::Always, InliningPolicy::Never] {
                let entries: Vec<SerializedFnMapEntry<String>> = vec![
                    (
                        print(&alias),
                        print(&x.pow(2)),
                        vec!["0".into()],
                        vec![],
                        InliningPolicy::Always,
                        true,
                    ),
                    (
                        print(&call(x.clone())),
                        print(&(&alias + &global)),
                        tags.iter().map(print).collect(),
                        vec![print(&x)],
                        inlining,
                        false,
                    ),
                ];
                let json = serde_json::to_vec(&entries)?;
                let parsed: Vec<SerializedFnMapEntry<String>> = serde_json::from_slice(&json)?;
                let bytes = bincode::encode_to_vec(&parsed, bincode::config::standard())?;
                let (restored, consumed): (Vec<SerializedFnMapEntry<String>>, _) =
                    bincode::decode_from_slice(&bytes, bincode::config::standard())?;
                assert_eq!(consumed, bytes.len());
                assert_eq!(restored, entries);
                let mut parsed = parse_fn_map_entries(&restored, &StateMap::default())?;
                parsed.extend(parsed.clone()); // Graph and evaluator metadata overlap.
                let (replacements, all_replacements, fn_map) = apply_fn_map_entries(parsed)?;
                let root = call(Atom::num(2)) + call(Atom::num(3));
                assert_eq!(root.replace_multiple(&replacements), root);
                if inlining == InliningPolicy::Never {
                    assert_eq!(root.replace_multiple(&all_replacements), root);
                }
                let mut evaluator = root
                    .replace_multiple(&replacements)
                    .evaluator(&[global.clone(), x.clone()])
                    .function_map(fn_map)
                    .horner_iterations(1)
                    .build()
                    .map_err(|error| eyre!(error))?
                    .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
                assert_eq!(
                    evaluator.export_instructions().sub_evaluators.len(),
                    usize::from(inlining == InliningPolicy::Never)
                );
                assert_eq!(
                    evaluator.evaluate_single(&[Complex::new(10.0, 0.0), Complex::new(11.0, 0.0)]),
                    Complex::new(33.0, 0.0)
                );
            }
        }
        Ok(())
    }

    #[test]
    fn standalone_roundtrip_restores_non_dense_residue_map_keys() -> Result<()> {
        let residue_map_id = Atom::var(GS.residue_map_id);
        let selector = |id: i64, coefficient: i64| {
            Symbol::IF.call_args([
                residue_map_id.clone() - Atom::num(id),
                Atom::Zero,
                Atom::num(coefficient),
            ])
        };
        let expression = selector(4, 2) + selector(9, 3);
        let print = |atom: &Atom| atom.printer(PrintOptions::file()).to_string();
        let generic = StandaloneGenericEvaluatorArchive {
            exprs: vec![print(&expression)],
            additional_fn_map_entries: Vec::new(),
            dual_shape: None,
            zero_components: Vec::new(),
        };
        let archive = StandaloneEvaluatorArchive {
            version: STANDALONE_EVALUATORS_VERSION,
            three_d_representations: vec!["cff".to_owned()],
            numeric_target: StandaloneNumericTarget::Double,
            symbolica_state: (),
            graph_terms: vec![StandaloneGraphTermArchive {
                graph_name: "duplicate_orientations".to_string(),
                orientations: vec![vec![1], vec![1]],
                param_builder_params: vec![print(&residue_map_id), print(&GS.sign(EdgeIndex(0)))],
                fn_map_entries: Vec::new(),
                original_integrand: BTreeMap::from([(
                    "cff".to_owned(),
                    StandaloneEvaluatorStackArchive {
                        explicit_orientation_sum_only: false,
                        production_orientation_ids: vec![4, 9],
                        parametric_rows: None,
                        single_parametric_auxiliary: None,
                        single_parametric: generic,
                        iterative: None,
                        summed_function_map: None,
                        summed: None,
                        representative_input: vec![
                            StandaloneComplexInput {
                                re: "0".to_string(),
                                im: "0".to_string(),
                            },
                            StandaloneComplexInput {
                                re: "0".to_string(),
                                im: "0".to_string(),
                            },
                        ],
                        start: 1,
                        residue_map_id_start: 0,
                        mult_offset: 1,
                    },
                )]),
                threshold_counterterms: Vec::new(),
                threshold_counterterms_are_variants: false,
                threshold_variants: Vec::new(),
                threshold_multipliers: None,
                metadata_registry: None,
            }],
        };

        let serialized = serde_json::to_vec(&archive)?;
        let roundtripped: StandaloneEvaluatorArchive<(), String> =
            serde_json::from_slice(&serialized)?;
        let mut loaded = roundtripped.load()?;
        let orientations = loaded.graph_terms[0].orientations.clone();
        let stack = loaded.graph_terms[0]
            .original_integrand
            .get_mut("cff")
            .unwrap();

        for (orientation_index, expected) in [(Some(0), 2.0), (Some(1), 3.0), (None, 5.0)] {
            let result = stack.evaluate_with_backend(StandaloneEvaluationRequest {
                backend: StandaloneBackend::Eager,
                method: StandaloneMethod::SingleParametric,
                orientations: &orientations,
                orientation_index,
                custom_input: None,
                artifact_root: Path::new("."),
                label: "standalone_residue_map_roundtrip",
            })?;
            assert_eq!(result, vec![Complex::new(expected, 0.0)]);
        }

        Ok(())
    }

    #[test]
    fn exported_threshold_inputs_preserve_native_rows_and_nonzero_jets() -> Result<()> {
        use crate::{
            cff::expression::OrientationID,
            initialisation::test_initialise,
            integrands::{
                evaluation::EvaluationMetaData,
                process::{
                    ParamBuilder,
                    amplitude::export::{ExportAtomTo, export_evaluator_stack},
                    cross_section::{export as cross_export, load as cross},
                    evaluators::{
                        EvaluatorMethod, EvaluatorStack, InputParams, SingleOrAllOrientations,
                        SliceMut,
                    },
                    param_builder::FnMapEntry,
                },
            },
            processes::EvaluatorSettings,
            settings::RuntimeSettings,
            utils::{F, hyperdual_utils::DualOrNot},
            uv::{ParametricIntegrandRow, ParametricIntegrandTerm},
        };
        use linnet::half_edge::{
            involution::{EdgeVec, Orientation},
            subgraph::{SubSetLike, subset::SubSet},
        };
        use spenso::algebra::complex::Complex as NativeComplex;
        use std::sync::Arc;
        use typed_index_collections::TiVec;

        test_initialise()?;
        let x = Atom::var(symbol!("standalone_export_test::x"));
        let y = Atom::var(symbol!("standalone_export_test::y"));
        let q = Atom::var(symbol!("standalone_export_test::q"));
        let mut builder = ParamBuilder::new_empty();
        builder.pairs.additional_params = [x.clone(), y.clone()].into_iter().collect();
        assert_eq!(builder.pairs.update_ranges(), 2);
        let definition = Arc::new(FnMapEntry {
            lhs: function!(symbol!("standalone_export_test::numerator"), 0, q.clone()),
            rhs: &q * &x,
            tags: vec![Atom::num(0)],
            args: vec![Indeterminate::try_from(q).unwrap()],
            inlining: InliningPolicy::Always,
            is_alias: false,
        });
        let settings = EvaluatorSettings {
            store_atom: true,
            summed_function_map: true,
            iterative_orientation_optimization: false,
            ..Default::default()
        };
        let orientations = TiVec::<OrientationID, EdgeVec<Orientation>>::new();
        let filter = SubSet::empty(0);
        // Sum q*x*(x^2+y) for q=1,3. Probe x=2+t,y=3+2t and
        // x=2+s,y=3+t; these are full Taylor jets, not derivative values.
        for (shape, input, expected) in [
            (vec![vec![0], vec![1]], vec![2., 1., 3., 2.], vec![56., 76.]),
            (
                vec![vec![0], vec![1], vec![2]],
                vec![2., 1., 0., 3., 2., 0.],
                vec![56., 76., 32.],
            ),
            (
                vec![vec![0, 0], vec![1, 0], vec![0, 1], vec![1, 1]],
                vec![2., 1., 0., 0., 3., 0., 1., 0.],
                vec![56., 60., 8., 4.],
            ),
        ] {
            let width = shape.len();
            let term = ParametricIntegrandTerm {
                numerator: Arc::clone(&definition),
                rows: [1, 3]
                    .into_iter()
                    .map(|q| ParametricIntegrandRow {
                        coefficients: vec![Rational::from(q)],
                        carrier: x.pow(2) + &y,
                        source_keys: Vec::new(),
                    })
                    .collect(),
            };
            let (mut native, _) = EvaluatorStack::from_integrand_with_timings(
                &Atom::zero(),
                &builder,
                &[],
                &[term],
                None,
                Some(shape.clone()),
                &settings,
            )?;
            assert_eq!(native.parametric_rows.as_ref().unwrap().rows.len(), 2);
            assert!(native.single_parametric_auxiliary.is_some());
            let exported = export_evaluator_stack::<String>(
                &native,
                0,
                0,
                [2, 3]
                    .into_iter()
                    .map(|value| StandaloneComplexInput {
                        re: value.to_string(),
                        im: "0".into(),
                    })
                    .collect(),
            )?;
            assert_eq!(exported.mult_offset, width);
            assert_eq!(exported.representative_input.len(), 2 * width);
            let defaults = exported
                .representative_input
                .iter()
                .map(StandaloneComplexInput::to_f64)
                .collect::<Result<Vec<_>>>()?;
            for (jet, base) in defaults.chunks_exact(width).zip([2., 3.]) {
                assert_eq!(jet[0], Complex::new(base, 0.));
                assert!(jet[1..].iter().all(|v| *v == Complex::new(0., 0.)));
            }
            assert_eq!(
                exported
                    .single_parametric_auxiliary
                    .as_ref()
                    .unwrap()
                    .1
                    .dual_shape,
                Some(shape.clone())
            );
            let cross_stack = cross_export::export_evaluator_stack::<String>(&native, 0, 0)?;
            assert_eq!(cross_stack.mult_offset, width);
            assert!(cross_stack.representative_input.is_empty());
            let indexed = cross::StandaloneIndexedEvaluatorStackArchive {
                cut_cff_index: cross::StandaloneCutCFFIndex {
                    left_threshold_order: None,
                    right_threshold_order: None,
                    lu_cut_order: Some(width),
                },
                evaluator_stack: cross_stack,
            };
            let cross_archive = cross::StandaloneCrossSectionArchive {
                version: cross::STANDALONE_EVALUATORS_VERSION,
                three_d_representations: vec!["cff".into(), "ltd".into()],
                symbolica_state: (),
                graph_terms: vec![cross::StandaloneCrossSectionGraphTermArchive {
                    graph_name: "threshold_jets".into(),
                    orientations: vec![],
                    param_builder_params: vec![
                        String::export_atom_to(&x)?,
                        String::export_atom_to(&y)?,
                    ],
                    fn_map_entries: vec![],
                    cut_group_integrands: BTreeMap::from(
                        ["cff", "ltd"].map(|rep| (rep.into(), vec![vec![indexed.clone()]])),
                    ),
                    counterterms: vec![cross::StandaloneCountertermArchive {
                        integrands: BTreeMap::from(["cff", "ltd"].map(|rep| {
                            (
                                rep.into(),
                                cross::StandaloneCountertermIntegrandsArchive {
                                    left_thresholds_evaluator: vec![vec![indexed.clone()]],
                                    right_thresholds_evaluator: vec![vec![indexed.clone()]],
                                    iterated_evaluator:
                                        cross::StandaloneIteratedCollectionArchive {
                                            data: vec![vec![indexed.clone()]],
                                            num_right_thresholds: 1,
                                        },
                                },
                            )
                        })),
                        left_threshold_helpers: vec![vec![]],
                        right_threshold_helpers: vec![vec![]],
                        iterated_helpers: cross::StandaloneIteratedCollectionArchive {
                            data: vec![vec![]],
                            num_right_thresholds: 1,
                        },
                        threshold_multipliers: None,
                        pass_two_evaluator: vec![],
                    }],
                    threshold_counterterm_metadata: None,
                }],
            };
            let cross_json: cross::StandaloneCrossSectionArchive<(), String> =
                serde_json::from_slice(&serde_json::to_vec(&cross_archive)?)?;
            let (cross_binary, _): (cross::StandaloneCrossSectionArchive<(), String>, _) =
                bincode::decode_from_slice(
                    &bincode::encode_to_vec(&cross_archive, bincode::config::standard())?,
                    bincode::config::standard(),
                )?;
            let mut cross_loaded = Vec::new();
            for restored in [cross_json, cross_binary] {
                // Every physical/left/right/iterated entry reaches the shared validator,
                // even though cross sections intentionally have no default input sample.
                for location in 0..4 {
                    for corruption in 0..6 {
                        let mut invalid = restored.clone();
                        let graph = &mut invalid.graph_terms[0];
                        let stack = match location {
                            0 => {
                                &mut graph.cut_group_integrands.get_mut("cff").unwrap()[0][0]
                                    .evaluator_stack
                            }
                            1 => {
                                &mut graph.counterterms[0]
                                    .integrands
                                    .get_mut("cff")
                                    .unwrap()
                                    .left_thresholds_evaluator[0][0]
                                    .evaluator_stack
                            }
                            2 => {
                                &mut graph.counterterms[0]
                                    .integrands
                                    .get_mut("cff")
                                    .unwrap()
                                    .right_thresholds_evaluator[0][0]
                                    .evaluator_stack
                            }
                            _ => {
                                &mut graph.counterterms[0]
                                    .integrands
                                    .get_mut("cff")
                                    .unwrap()
                                    .iterated_evaluator
                                    .data[0][0]
                                    .evaluator_stack
                            }
                        };
                        match corruption {
                            0 => stack.mult_offset = 1,
                            1 => stack.representative_input = vec![Complex::new(0., 0.)],
                            2 => stack.summed_function_map.as_mut().unwrap().dual_shape = None,
                            3 => stack.start = 3,
                            4 => {
                                stack
                                    .single_parametric_auxiliary
                                    .as_mut()
                                    .unwrap()
                                    .1
                                    .dual_shape = None
                            }
                            _ => stack.single_parametric.dual_shape = Some(vec![]),
                        }
                        assert!(
                            invalid
                                .load()
                                .err()
                                .expect("invalid cross jet layout")
                                .to_string()
                                .contains("physical input or dual layout")
                        );
                    }
                }
                cross_loaded.push(restored.load()?);
            }
            let archive = StandaloneEvaluatorArchive {
                version: STANDALONE_EVALUATORS_VERSION,
                three_d_representations: vec!["cff".into(), "ltd".into()],
                numeric_target: StandaloneNumericTarget::Double,
                symbolica_state: (),
                graph_terms: vec![StandaloneGraphTermArchive {
                    graph_name: "threshold_jets".into(),
                    orientations: vec![],
                    param_builder_params: vec![
                        String::export_atom_to(&x)?,
                        String::export_atom_to(&y)?,
                    ],
                    fn_map_entries: vec![],
                    original_integrand: BTreeMap::from([
                        ("cff".into(), exported.clone()),
                        ("ltd".into(), exported.clone()),
                    ]),
                    threshold_counterterms: vec![BTreeMap::from(["cff", "ltd"].map(|rep| {
                        (
                            rep.into(),
                            vec![StandaloneIndexedEvaluatorStackArchive {
                                cut_cff_index: StandaloneCutCFFIndex {
                                    left_threshold_order: Some(if shape[0].len() == 2 {
                                        2
                                    } else {
                                        width
                                    }),
                                    right_threshold_order: (shape[0].len() == 2).then_some(2),
                                    lu_cut_order: None,
                                },
                                evaluator_stack: exported.clone(),
                            }],
                        )
                    }))],
                    threshold_counterterms_are_variants: false,
                    threshold_variants: vec![],
                    threshold_multipliers: None,
                    metadata_registry: None,
                }],
            };
            let binary = bincode::encode_to_vec(&archive, bincode::config::standard())?;
            let (binary, _): (StandaloneEvaluatorArchive<(), String>, _) =
                bincode::decode_from_slice(&binary, bincode::config::standard())?;
            let json = serde_json::from_slice::<StandaloneEvaluatorArchive<(), String>>(
                &serde_json::to_vec(&archive)?,
            )?;
            for restored in [binary, json] {
                for corruption in 0..6 {
                    let mut malformed = restored.clone();
                    let stack = &mut malformed.graph_terms[0].threshold_counterterms[0]
                        .get_mut("cff")
                        .unwrap()[0]
                        .evaluator_stack;
                    match corruption {
                        0 => stack.mult_offset = 1,
                        1 => {
                            stack.representative_input.pop();
                        }
                        2 => stack.summed_function_map.as_mut().unwrap().dual_shape = None,
                        3 => stack.start = 3,
                        4 => {
                            stack
                                .single_parametric_auxiliary
                                .as_mut()
                                .unwrap()
                                .1
                                .dual_shape = None
                        }
                        5 => stack.single_parametric.dual_shape = Some(Vec::new()),
                        _ => unreachable!(),
                    }
                    let error = malformed.load().err().expect("invalid physical jet layout");
                    assert!(error.to_string().contains("physical input or dual layout"));
                }
                let mut loaded = restored.load()?;
                for representation in ["cff", "ltd"] {
                    let loaded = loaded.graph_terms[0].stack_mut(
                        representation,
                        &StandaloneStackSelection::ThresholdCounterterm((0, 0)),
                    )?;
                    for (method, native_method) in [
                        (
                            StandaloneMethod::SingleParametric,
                            EvaluatorMethod::SingleParametric,
                        ),
                        (
                            StandaloneMethod::SummedFunctionMap,
                            EvaluatorMethod::SummedFunctionMap,
                        ),
                    ] {
                        for custom in [false, true] {
                            let values = if custom {
                                input
                                    .iter()
                                    .map(|&v| Complex::new(v, 0.))
                                    .collect::<Vec<_>>()
                            } else {
                                defaults.clone()
                            };
                            let mut runtime = RuntimeSettings::default();
                            runtime.general.evaluator_method = native_method.clone();
                            let original = native.evaluate(
                                InputParams {
                                    values: SliceMut::Owned(
                                        values
                                            .iter()
                                            .map(|v| NativeComplex::new(F(v.re), F(v.im)))
                                            .collect(),
                                    ),
                                    residue_map_id_start: 0,
                                    orientations_start: 0,
                                    multiplicative_offset: width,
                                },
                                SingleOrAllOrientations::All {
                                    all: &orientations,
                                    filter: &filter,
                                },
                                &runtime,
                                &mut EvaluationMetaData::new_empty(),
                            )?;
                            let [DualOrNot::Dual(original)] = original.as_slice() else {
                                panic!("expected one full jet")
                            };
                            let original = original
                                .values
                                .iter()
                                .map(|v| Complex::new(v.re.0, v.im.0))
                                .collect::<Vec<_>>();
                            let actual =
                                loaded.evaluate_with_backend(StandaloneEvaluationRequest {
                                    backend: StandaloneBackend::Eager,
                                    method,
                                    orientations: &[],
                                    orientation_index: None,
                                    custom_input: custom.then_some(values.as_slice()),
                                    artifact_root: Path::new("."),
                                    label: "exported_threshold_jets",
                                })?;
                            assert_eq!(actual, original);
                            for cross in &mut cross_loaded {
                                let graph = &mut cross.graph_terms[0];
                                let counterterm = graph.counterterms[0]
                                    .integrands
                                    .get_mut(representation)
                                    .unwrap();
                                let key = indexed.cut_cff_index;
                                for stack in [
                                    graph.cut_group_integrands.get_mut(representation).unwrap()[0]
                                        .get_mut(&key)
                                        .unwrap(),
                                    counterterm.left_thresholds_evaluator[0]
                                        .get_mut(&key)
                                        .unwrap(),
                                    counterterm.right_thresholds_evaluator[0]
                                        .get_mut(&key)
                                        .unwrap(),
                                    counterterm.iterated_evaluator.data[0]
                                        .get_mut(&key)
                                        .unwrap(),
                                ] {
                                    assert!(stack.representative_input.is_empty());
                                    let observed = match method {
                                        StandaloneMethod::SingleParametric => {
                                            stack.evaluate_single_parametric(&values)?
                                        }
                                        StandaloneMethod::SummedFunctionMap => {
                                            let (_, _, evaluator, output) =
                                                stack.summed_function_map.as_mut().unwrap();
                                            evaluator.evaluate(&values, output);
                                            output.clone()
                                        }
                                        _ => unreachable!(),
                                    };
                                    assert_eq!(observed, original);
                                }
                            }

                            if custom {
                                assert_eq!(
                                    actual,
                                    expected
                                        .iter()
                                        .map(|&v| Complex::new(v, 0.))
                                        .collect::<Vec<_>>()
                                );
                            } else {
                                assert_eq!(actual[0], Complex::new(56., 0.));
                                assert!(actual[1..].iter().all(|v| *v == Complex::new(0., 0.)));
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn standalone_roundtrip_keeps_representation_order_and_selects_each_program() -> Result<()> {
        // Noncontiguous port names and output gathers preserve the same x+c
        // argument: the two dynamic ports cancel their runtime q binding.
        // The final inactive row must clear both complete gathered jets.
        let stack = |value: i64| {
            let parametric = StandaloneGenericEvaluatorArchive {
                exprs: vec![format!(
                    "{value}*if(active,standalone_joint_numerator(standalone_auxiliary_port(2)-2*standalone_auxiliary_port(7)+standalone_auxiliary_port(9)+c),0)"
                )],
                additional_fn_map_entries: vec![(
                    "standalone_joint_numerator(u)".to_owned(),
                    "u^2".to_owned(),
                    Vec::new(),
                    vec!["u".to_owned()],
                    InliningPolicy::Never,
                    false,
                )],
                dual_shape: Some(vec![vec![0], vec![1]]),
                zero_components: vec![(4, 1), (5, 1)],
            };
            let mut summed = parametric.clone();
            summed.exprs = vec![format!(
                "{value}*(standalone_joint_numerator(x+1)+standalone_joint_numerator(x+3))"
            )];
            summed.zero_components.clear();
            StandaloneEvaluatorStackArchive {
                explicit_orientation_sum_only: true,
                production_orientation_ids: Vec::new(),
                parametric_rows: Some(StandaloneParametricResidueRows {
                    parameters: vec!["unused".to_owned(), "c".to_owned(), "active".to_owned()],
                    rows: vec![
                        vec![Rational::from(8), Rational::from(1), Rational::from(1)],
                        vec![Rational::from(9), Rational::from(3), Rational::from(1)],
                        vec![Rational::from(10), Rational::from(999), Rational::from(0)],
                    ],
                    input_columns: vec![1, 2],
                    auxiliary_parameters: vec!["q1".to_owned(), "q3".to_owned()],
                    auxiliary_values: vec![Rational::from(1), Rational::from(3)],
                    auxiliary_shared_outputs: 1,
                    auxiliary_gathers: vec![
                        vec![Some(2), Some(4)],
                        vec![Some(1), Some(3)],
                        vec![None, None],
                    ],
                }),
                single_parametric_auxiliary: Some((
                    vec![
                        "standalone_auxiliary_port(7)".to_owned(),
                        "standalone_auxiliary_port(2)".to_owned(),
                        "standalone_auxiliary_port(9)".to_owned(),
                    ],
                    StandaloneGenericEvaluatorArchive {
                        exprs: vec![
                            "x+1".to_owned(),
                            "3*x+2+q3".to_owned(),
                            "3*x+2+q1".to_owned(),
                            "-q3".to_owned(),
                            "-q1".to_owned(),
                        ],
                        additional_fn_map_entries: Vec::new(),
                        dual_shape: Some(vec![vec![0], vec![1]]),
                        zero_components: vec![(1, 1), (2, 1)],
                    },
                )),
                single_parametric: parametric,
                iterative: None,
                summed_function_map: Some(summed),
                summed: None,
                representative_input: vec![
                    StandaloneComplexInput {
                        re: "2".to_owned(),
                        im: "0".to_owned(),
                    },
                    StandaloneComplexInput {
                        re: "1".to_owned(),
                        im: "0".to_owned(),
                    },
                ],
                start: 0,
                residue_map_id_start: 0,
                mult_offset: 2,
            }
        };
        let archive = StandaloneEvaluatorArchive {
            version: STANDALONE_EVALUATORS_VERSION,
            three_d_representations: vec!["ltd".to_owned(), "cff".to_owned()],
            numeric_target: StandaloneNumericTarget::Double,
            symbolica_state: (),
            graph_terms: vec![StandaloneGraphTermArchive {
                graph_name: "representation_dispatch".to_owned(),
                orientations: Vec::new(),
                param_builder_params: vec!["x".to_owned()],
                fn_map_entries: Vec::new(),
                original_integrand: BTreeMap::from([
                    ("cff".to_owned(), stack(3)),
                    ("ltd".to_owned(), stack(7)),
                ]),
                threshold_counterterms: Vec::new(),
                threshold_counterterms_are_variants: false,
                threshold_variants: Vec::new(),
                threshold_multipliers: None,
                metadata_registry: None,
            }],
        };
        let encoded = bincode::encode_to_vec(&archive, bincode::config::standard())?;
        let (binary, _): (StandaloneEvaluatorArchive<(), String>, _) =
            bincode::decode_from_slice(&encoded, bincode::config::standard())?;
        let json = serde_json::from_slice::<StandaloneEvaluatorArchive<(), String>>(
            &serde_json::to_vec(&archive)?,
        )?;
        for restored in [binary, json] {
            let mut loaded = restored.load()?;
            assert_eq!(loaded.three_d_representations, ["ltd", "cff"]);
            for (representation, coefficient) in [("cff", 3.0), ("ltd", 7.0)] {
                let stack = loaded.graph_terms[0]
                    .stack_mut(representation, &StandaloneStackSelection::Original)?;
                assert_eq!(
                    stack.single_parametric_auxiliary.as_ref().unwrap().0,
                    [
                        parse_lit!(standalone_auxiliary_port(7)),
                        parse_lit!(standalone_auxiliary_port(2)),
                        parse_lit!(standalone_auxiliary_port(9))
                    ]
                );
                for method in [
                    StandaloneMethod::SingleParametric,
                    StandaloneMethod::SummedFunctionMap,
                ] {
                    for (custom_input, expected) in [
                        (None, [34.0, 16.0]),
                        (
                            Some([Complex::new(3.0, 0.0), Complex::new(1.0, 0.0)]),
                            [52.0, 20.0],
                        ),
                    ] {
                        let result = stack.evaluate_with_backend(StandaloneEvaluationRequest {
                            backend: StandaloneBackend::Eager,
                            method,
                            orientations: &[],
                            orientation_index: None,
                            custom_input: custom_input.as_ref().map(|input| input.as_slice()),
                            artifact_root: Path::new("."),
                            label: "representation_dispatch",
                        })?;
                        assert_eq!(
                            result,
                            expected.map(|value| Complex::new(coefficient * value, 0.0))
                        );
                    }
                }
            }
            assert!(
                loaded.graph_terms[0]
                    .stack_mut("unavailable", &StandaloneStackSelection::Original)
                    .is_err()
            );
        }
        for corrupt_shape in [false, true] {
            let mut malformed = archive.clone();
            let auxiliary = malformed.graph_terms[0]
                .original_integrand
                .get_mut("cff")
                .unwrap()
                .single_parametric_auxiliary
                .as_mut()
                .unwrap();
            if corrupt_shape {
                auxiliary.1.dual_shape = None;
            } else {
                auxiliary.0.clear();
            }
            assert!(malformed.load().is_err());
        }
        for conflicting_port in ["x", "c"] {
            let mut malformed = archive.clone();
            malformed.graph_terms[0]
                .original_integrand
                .get_mut("cff")
                .unwrap()
                .single_parametric_auxiliary
                .as_mut()
                .unwrap()
                .0[0] = conflicting_port.to_owned();
            let error = match malformed.load() {
                Ok(_) => panic!("auxiliary port collides with an existing input"),
                Err(error) => error,
            };
            assert!(error.to_string().contains("parameters must be distinct"));
        }
        for corruption in 0..8 {
            let mut malformed = archive.clone();
            let catalog = malformed.graph_terms[0]
                .original_integrand
                .get_mut("cff")
                .unwrap()
                .parametric_rows
                .as_mut()
                .unwrap();
            match corruption {
                0 => catalog.input_columns = vec![2, 1],
                1 => catalog.input_columns.push(3),
                2 => catalog.auxiliary_values.clear(),
                3 => catalog.auxiliary_shared_outputs = 4,
                4 => catalog.auxiliary_gathers.truncate(1),
                5 => catalog.auxiliary_gathers[0].push(Some(0)),
                6 => catalog.auxiliary_gathers[0][0] = Some(5),
                _ => catalog.input_columns = vec![1, 1],
            }
            assert!(
                malformed.load().is_err(),
                "accepted malformed auxiliary layout {corruption}"
            );
        }
        let mut inconsistent = archive;
        inconsistent.three_d_representations = vec!["cff".to_owned()];
        assert!(inconsistent.load().is_err());
        Ok(())
    }
}
