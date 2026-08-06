use std::{collections::BTreeMap, fs, path::Path};

use color_eyre::Result;
use eyre::{Context, eyre};
use linnet::half_edge::involution::Orientation;
use symbolica::state::State;

use super::{
    CrossSectionIntegrand,
    load::{
        STANDALONE_EVALUATORS_VERSION, StandaloneCountertermArchive, StandaloneCrossSectionArchive,
        StandaloneCrossSectionGraphTermArchive, StandaloneCutCFFIndex,
        StandaloneEvaluatorStackArchive, StandaloneGenericEvaluatorArchive,
        StandaloneIndexedEvaluatorStackArchive, StandaloneIndexedGenericEvaluatorArchive,
        StandaloneIteratedCollectionArchive, StandaloneThresholdCountertermMetadataRegistry,
        StandaloneThresholdMultiplierCollectionArchive, StandaloneThresholdMultiplierEsurface,
        StandaloneThresholdMultiplierInput, StandaloneThresholdMultiplierLayoutArchive,
        StandaloneThresholdMultiplierPoint, StandaloneThresholdMultiplierVariantReference,
    },
};
use crate::{
    cff::CutCFFIndex,
    integrands::process::{
        GenericEvaluator,
        amplitude::export::ExportAtomTo,
        threshold_multiplier::{
            ThresholdMultiplierEvaluatorCollection, ThresholdMultiplierInput,
            ThresholdMultiplierPoint,
        },
    },
    processes::{
        IteratedCtCollection, StandaloneDataFormat, StandaloneExportMode, StandaloneExportSettings,
        StandaloneNumericTarget, ThresholdCountertermMetadataRegistry,
    },
    subtraction::lu_counterterm::LUCounterTermEvaluators,
};

const STANDALONE_DATA_FILE: &str = "standalone_cross_section";
const STANDALONE_RUST_SCRIPT_FILE: &str = "standalone_cross_section_rust.rs";

#[cfg(unix)]
fn make_script_executable(path: &Path) -> Result<()> {
    use std::os::unix::fs::PermissionsExt;

    let mut permissions = fs::metadata(path)?.permissions();
    permissions.set_mode(0o755);
    fs::set_permissions(path, permissions)?;
    Ok(())
}

#[cfg(not(unix))]
fn make_script_executable(_: &Path) -> Result<()> {
    Ok(())
}

fn export_generic_evaluator<T: ExportAtomTo>(
    evaluator: &GenericEvaluator,
) -> Result<StandaloneGenericEvaluatorArchive<T>> {
    let exprs = evaluator
        .exprs
        .as_ref()
        .ok_or_else(|| {
            eyre!("Standalone export requires stored evaluator atoms (`store_atom=true`)")
        })?
        .iter()
        .map(T::export_atom_to)
        .collect::<Result<Vec<_>>>()?;

    let additional_fn_map_entries = evaluator
        .fn_map_entries
        .iter()
        .map(|entry| entry.archive())
        .collect::<Result<Vec<_>>>()?;

    Ok(StandaloneGenericEvaluatorArchive {
        exprs,
        additional_fn_map_entries,
        dual_shape: evaluator.dual_shape.clone(),
    })
}

fn export_evaluator_stack<T: ExportAtomTo>(
    stack: &crate::integrands::process::evaluators::EvaluatorStack,
) -> Result<StandaloneEvaluatorStackArchive<T>> {
    Ok(StandaloneEvaluatorStackArchive {
        single_parametric: export_generic_evaluator(&stack.single_parametric)?,
        iterative: stack
            .iterative
            .as_ref()
            .map(|(eval, _)| export_generic_evaluator(eval))
            .transpose()?,
        summed_function_map: stack
            .summed_function_map
            .as_ref()
            .map(export_generic_evaluator)
            .transpose()?,
        summed: stack
            .summed
            .as_ref()
            .map(export_generic_evaluator)
            .transpose()?,
        representative_input: Vec::new(),
        start: 0,
        mult_offset: 0,
    })
}

fn export_cut_cff_index(index: &CutCFFIndex) -> StandaloneCutCFFIndex {
    StandaloneCutCFFIndex {
        left_threshold_order: index.left_threshold_order,
        right_threshold_order: index.right_threshold_order,
        lu_cut_order: index.lu_cut_order,
    }
}

fn export_evaluator_map<T: ExportAtomTo>(
    stacks: &BTreeMap<CutCFFIndex, crate::integrands::process::evaluators::EvaluatorStack>,
) -> Result<Vec<StandaloneIndexedEvaluatorStackArchive<T>>> {
    stacks
        .iter()
        .map(|(cut_cff_index, stack)| {
            Ok(StandaloneIndexedEvaluatorStackArchive {
                cut_cff_index: export_cut_cff_index(cut_cff_index),
                evaluator_stack: export_evaluator_stack(stack)?,
            })
        })
        .collect()
}

fn export_generic_evaluator_map<T: ExportAtomTo>(
    evaluators: &BTreeMap<CutCFFIndex, GenericEvaluator>,
) -> Result<Vec<StandaloneIndexedGenericEvaluatorArchive<T>>> {
    evaluators
        .iter()
        .map(|(cut_cff_index, evaluator)| {
            Ok(StandaloneIndexedGenericEvaluatorArchive {
                cut_cff_index: export_cut_cff_index(cut_cff_index),
                evaluator: export_generic_evaluator(evaluator)?,
            })
        })
        .collect()
}

fn export_iterated_collection<U, V, F>(
    collection: &IteratedCtCollection<U>,
    export_item: F,
) -> Result<StandaloneIteratedCollectionArchive<V>>
where
    F: Fn(&U) -> Result<V>,
{
    Ok(StandaloneIteratedCollectionArchive {
        data: collection
            .iter()
            .map(export_item)
            .collect::<Result<Vec<_>>>()?,
        num_right_thresholds: collection.num_right_thresholds(),
    })
}

fn export_threshold_multiplier_point(
    point: ThresholdMultiplierPoint,
) -> StandaloneThresholdMultiplierPoint {
    match point {
        ThresholdMultiplierPoint::Effective => StandaloneThresholdMultiplierPoint::Effective,
        ThresholdMultiplierPoint::Star => StandaloneThresholdMultiplierPoint::Star,
    }
}

fn export_threshold_multiplier_input(
    input: &ThresholdMultiplierInput,
) -> StandaloneThresholdMultiplierInput {
    match *input {
        ThresholdMultiplierInput::ModelParameter { index } => {
            StandaloneThresholdMultiplierInput::ModelParameter { index }
        }
        ThresholdMultiplierInput::AdditionalParameter { index } => {
            StandaloneThresholdMultiplierInput::AdditionalParameter { index }
        }
        ThresholdMultiplierInput::ExternalMomentum {
            position,
            component,
        } => StandaloneThresholdMultiplierInput::ExternalMomentum {
            position,
            component,
        },
        ThresholdMultiplierInput::EdgeMomentum {
            point,
            edge,
            component,
        } => StandaloneThresholdMultiplierInput::EdgeMomentum {
            point: export_threshold_multiplier_point(point),
            edge,
            component,
        },
        ThresholdMultiplierInput::Esurface { point, esurface } => {
            StandaloneThresholdMultiplierInput::Esurface {
                point: export_threshold_multiplier_point(point),
                esurface,
            }
        }
    }
}

pub(crate) fn export_threshold_multiplier_collection<T: ExportAtomTo>(
    collection: &ThresholdMultiplierEvaluatorCollection,
) -> Result<StandaloneThresholdMultiplierCollectionArchive<T>> {
    collection.validate(
        collection.left_variants().len(),
        collection.right_variants().len(),
    )?;
    let layout = collection.layout();
    Ok(StandaloneThresholdMultiplierCollectionArchive {
        layout: StandaloneThresholdMultiplierLayoutArchive {
            model_parameter_count: layout.model_parameters().len(),
            additional_parameters: layout
                .additional_parameters()
                .iter()
                .map(T::export_atom_to)
                .collect::<Result<Vec<_>>>()?,
            external_count: layout.external_count(),
            edges: layout.edges().to_vec(),
            esurfaces: layout
                .esurfaces()
                .iter()
                .map(|esurface| StandaloneThresholdMultiplierEsurface {
                    edges: esurface.edges.clone(),
                    external_shift: esurface.external_shift.clone(),
                })
                .collect(),
            inputs: layout
                .inputs()
                .iter()
                .map(export_threshold_multiplier_input)
                .collect(),
            parameters: layout
                .parameters()
                .iter()
                .map(T::export_atom_to)
                .collect::<Result<Vec<_>>>()?,
        },
        evaluators: collection
            .evaluators()
            .iter()
            .map(|evaluator| {
                Ok(StandaloneGenericEvaluatorArchive {
                    exprs: vec![T::export_atom_to(evaluator.expression())?],
                    additional_fn_map_entries: Vec::new(),
                    dual_shape: None,
                })
            })
            .collect::<Result<Vec<_>>>()?,
        left_variants: collection
            .left_variants()
            .iter()
            .map(|reference| StandaloneThresholdMultiplierVariantReference {
                variant_id: reference.variant_id.0,
                evaluator_id: reference.evaluator_id.map(|id| id.0),
            })
            .collect(),
        right_variants: collection
            .right_variants()
            .iter()
            .map(|reference| StandaloneThresholdMultiplierVariantReference {
                variant_id: reference.variant_id.0,
                evaluator_id: reference.evaluator_id.map(|id| id.0),
            })
            .collect(),
    })
}

fn export_counterterm<T: ExportAtomTo>(
    evaluators: &LUCounterTermEvaluators,
) -> Result<StandaloneCountertermArchive<T>> {
    if let Some(multipliers) = &evaluators.threshold_multipliers {
        multipliers.validate(
            evaluators.left_thresholds_evaluator.len(),
            evaluators.right_thresholds_evaluator.len(),
        )?;
    }
    Ok(StandaloneCountertermArchive {
        left_thresholds_evaluator: evaluators
            .left_thresholds_evaluator
            .iter()
            .map(export_evaluator_map)
            .collect::<Result<Vec<_>>>()?,
        right_thresholds_evaluator: evaluators
            .right_thresholds_evaluator
            .iter()
            .map(export_evaluator_map)
            .collect::<Result<Vec<_>>>()?,
        iterated_evaluator: export_iterated_collection(
            &evaluators.iterated_evaluator,
            export_evaluator_map,
        )?,
        left_threshold_helpers: evaluators
            .threshold_helpers
            .left_thresholds
            .iter()
            .map(export_generic_evaluator_map)
            .collect::<Result<Vec<_>>>()?,
        right_threshold_helpers: evaluators
            .threshold_helpers
            .right_thresholds
            .iter()
            .map(export_generic_evaluator_map)
            .collect::<Result<Vec<_>>>()?,
        iterated_helpers: export_iterated_collection(
            &evaluators.threshold_helpers.iterated,
            export_generic_evaluator_map,
        )?,
        threshold_multipliers: evaluators
            .threshold_multipliers
            .as_ref()
            .map(export_threshold_multiplier_collection)
            .transpose()?,
        pass_two_evaluator: evaluators
            .residue_from_e_surface_evaluators
            .iter()
            .map(export_generic_evaluator)
            .collect::<Result<Vec<_>>>()?,
    })
}

fn export_threshold_counterterm_metadata(
    registry: &ThresholdCountertermMetadataRegistry,
) -> Result<StandaloneThresholdCountertermMetadataRegistry> {
    serde_json::from_value(serde_json::to_value(registry)?).with_context(|| {
        format!(
            "Failed to convert threshold-counterterm metadata for graph '{}' to the standalone schema",
            registry.graph_name,
        )
    })
}

fn standalone_rust_script() -> String {
    let mut script = include_str!("load.rs").to_string();
    if let Some(rest) = script.strip_prefix("//#!/usr/bin/env -S rust-script\n") {
        script = format!("#!/usr/bin/env -S rust-script\n{rest}");
    }

    script
}

impl CrossSectionIntegrand {
    pub(crate) fn export_standalone(
        &self,
        path: impl AsRef<Path>,
        settings: &StandaloneExportSettings,
    ) -> Result<()> {
        if settings.precision != StandaloneNumericTarget::Double {
            return Err(eyre!(
                "Cross-section standalone export currently only supports double precision",
            ));
        }

        let mut symbolica_state = Vec::new();
        State::export(&mut symbolica_state)
            .with_context(|| "Failed to export Symbolica state for standalone cross section")?;

        match settings.format {
            StandaloneDataFormat::Binary => {
                let standalone: StandaloneCrossSectionArchive<Vec<u8>, Vec<u8>> =
                    self.create_archive(symbolica_state)?;
                let binary = bincode::encode_to_vec(&standalone, bincode::config::standard())?;
                let standalone_path = path.as_ref().join(format!("{STANDALONE_DATA_FILE}.bin"));
                fs::write(&standalone_path, binary).with_context(|| {
                    format!(
                        "Failed writing standalone cross section binary to {}",
                        standalone_path.display()
                    )
                })?;
            }
            StandaloneDataFormat::Json => {
                let standalone: StandaloneCrossSectionArchive<(), String> =
                    self.create_archive(())?;
                let json = serde_json::to_vec_pretty(&standalone)?;
                let standalone_path = path.as_ref().join(format!("{STANDALONE_DATA_FILE}.json"));
                fs::write(&standalone_path, json).with_context(|| {
                    format!(
                        "Failed writing standalone cross section JSON to {}",
                        standalone_path.display()
                    )
                })?;
            }
        }

        match settings.mode {
            StandaloneExportMode::Python => {
                return Err(eyre!("Python Export mode not implemented"));
            }
            StandaloneExportMode::Rust => {
                let script = standalone_rust_script();
                let script_path = path.as_ref().join(STANDALONE_RUST_SCRIPT_FILE);
                fs::write(&script_path, script).with_context(|| {
                    format!(
                        "Failed writing standalone cross section rust-script loader to {}",
                        script_path.display()
                    )
                })?;
                make_script_executable(&script_path)?;
            }
        }

        Ok(())
    }

    fn create_archive<T: ExportAtomTo, S>(
        &self,
        symbolica_state: S,
    ) -> Result<StandaloneCrossSectionArchive<S, T>> {
        let graph_terms = self
            .data
            .graph_terms
            .iter()
            .map(|term| {
                let param_builder_params = (&term.param_builder.pairs)
                    .into_iter()
                    .flat_map(|pair| pair.params.iter())
                    .map(T::export_atom_to)
                    .collect::<Result<Vec<_>>>()?;

                let orientations = term
                    .orientations
                    .iter()
                    .map(|orientation| {
                        orientation
                            .iter()
                            .map(|(_, value)| match value {
                                Orientation::Default => 1i8,
                                Orientation::Reversed => -1i8,
                                Orientation::Undirected => 0,
                            })
                            .collect()
                    })
                    .collect();

                let fn_map_entries = term
                    .param_builder
                    .reps
                    .iter()
                    .map(|entry| entry.archive())
                    .collect::<Result<Vec<_>>>()?;

                let cut_group_integrands = term
                    .integrand
                    .iter()
                    .map(export_evaluator_map)
                    .collect::<Result<Vec<_>>>()?;

                let counterterms = term
                    .counterterm
                    .evaluators
                    .iter()
                    .map(export_counterterm)
                    .collect::<Result<Vec<_>>>()?;
                let threshold_counterterm_metadata = term
                    .counterterm
                    .metadata_registry
                    .as_ref()
                    .map(export_threshold_counterterm_metadata)
                    .transpose()?;

                Ok(StandaloneCrossSectionGraphTermArchive {
                    graph_name: term.graph.name.clone(),
                    orientations,
                    param_builder_params,
                    fn_map_entries,
                    cut_group_integrands,
                    counterterms,
                    threshold_counterterm_metadata,
                })
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(StandaloneCrossSectionArchive {
            version: STANDALONE_EVALUATORS_VERSION,
            symbolica_state,
            graph_terms,
        })
    }
}
