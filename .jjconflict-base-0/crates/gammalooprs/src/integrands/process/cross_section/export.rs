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
        StandaloneIteratedCollectionArchive,
    },
};
use crate::{
    cff::CutCFFIndex,
    integrands::process::{GenericEvaluator, amplitude::export::ExportAtomTo},
    processes::{
        IteratedCtCollection, StandaloneDataFormat, StandaloneExportMode, StandaloneExportSettings,
        StandaloneNumericTarget,
        cross_section::{CrossSectionGraph, params_for_derivative_order},
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
    parameter_override: Option<&[symbolica::atom::Atom]>,
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
        parameter_override: parameter_override
            .map(|params| {
                params
                    .iter()
                    .map(T::export_atom_to)
                    .collect::<Result<Vec<_>>>()
            })
            .transpose()?,
        additional_fn_map_entries,
        dual_shape: evaluator.dual_shape.clone(),
    })
}

fn export_evaluator_stack<T: ExportAtomTo>(
    stack: &crate::integrands::process::evaluators::EvaluatorStack,
    orientation_start: usize,
    residue_map_id_start: usize,
    mult_offset: usize,
) -> Result<StandaloneEvaluatorStackArchive<T>> {
    Ok(StandaloneEvaluatorStackArchive {
        explicit_orientation_sum_only: stack.explicit_orientation_sum_only,
        production_orientation_ids: stack
            .production_orientation_ids()
            .iter()
            .map(|id| id.0)
            .collect(),
        single_parametric: export_generic_evaluator(&stack.single_parametric, None)?,
        iterative: stack
            .iterative
            .as_ref()
            .map(|(eval, _)| export_generic_evaluator(eval, None))
            .transpose()?,
        summed_function_map: stack
            .summed_function_map
            .as_ref()
            .map(|eval| export_generic_evaluator(eval, None))
            .transpose()?,
        summed: stack
            .summed
            .as_ref()
            .map(|eval| export_generic_evaluator(eval, None))
            .transpose()?,
        representative_input: Vec::new(),
        start: orientation_start,
        residue_map_id_start,
        mult_offset,
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
    orientation_start: usize,
    residue_map_id_start: usize,
    mult_offset: usize,
) -> Result<Vec<StandaloneIndexedEvaluatorStackArchive<T>>> {
    stacks
        .iter()
        .map(|(cut_cff_index, stack)| {
            Ok(StandaloneIndexedEvaluatorStackArchive {
                cut_cff_index: export_cut_cff_index(cut_cff_index),
                evaluator_stack: export_evaluator_stack(
                    stack,
                    orientation_start,
                    residue_map_id_start,
                    mult_offset,
                )?,
            })
        })
        .collect()
}

fn export_generic_evaluator_map<T: ExportAtomTo>(
    evaluators: &BTreeMap<CutCFFIndex, GenericEvaluator>,
    parameter_override: impl Fn(&CutCFFIndex) -> Result<Vec<symbolica::atom::Atom>>,
) -> Result<Vec<StandaloneIndexedGenericEvaluatorArchive<T>>> {
    evaluators
        .iter()
        .map(|(cut_cff_index, evaluator)| {
            let params = parameter_override(cut_cff_index)?;
            Ok(StandaloneIndexedGenericEvaluatorArchive {
                cut_cff_index: export_cut_cff_index(cut_cff_index),
                evaluator: export_generic_evaluator(evaluator, Some(&params))?,
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

fn export_counterterm<T: ExportAtomTo>(
    evaluators: &LUCounterTermEvaluators,
    orientation_start: usize,
    residue_map_id_start: usize,
    mult_offset: usize,
) -> Result<StandaloneCountertermArchive<T>> {
    Ok(StandaloneCountertermArchive {
        left_thresholds_evaluator: evaluators
            .left_thresholds_evaluator
            .iter()
            .map(|stacks| {
                export_evaluator_map(stacks, orientation_start, residue_map_id_start, mult_offset)
            })
            .collect::<Result<Vec<_>>>()?,
        right_thresholds_evaluator: evaluators
            .right_thresholds_evaluator
            .iter()
            .map(|stacks| {
                export_evaluator_map(stacks, orientation_start, residue_map_id_start, mult_offset)
            })
            .collect::<Result<Vec<_>>>()?,
        iterated_evaluator: export_iterated_collection(&evaluators.iterated_evaluator, |stacks| {
            export_evaluator_map(stacks, orientation_start, residue_map_id_start, mult_offset)
        })?,
        left_threshold_helpers: evaluators
            .threshold_helpers
            .left_thresholds
            .iter()
            .map(|helpers| {
                export_generic_evaluator_map(helpers, |index| {
                    let order = index.left_threshold_order.ok_or_else(|| {
                        eyre!("Left threshold helper is missing its threshold order")
                    })?;
                    Ok(CrossSectionGraph::single_th_prefactor_helper_params(
                        order as u8,
                        false,
                    ))
                })
            })
            .collect::<Result<Vec<_>>>()?,
        right_threshold_helpers: evaluators
            .threshold_helpers
            .right_thresholds
            .iter()
            .map(|helpers| {
                export_generic_evaluator_map(helpers, |index| {
                    let order = index.right_threshold_order.ok_or_else(|| {
                        eyre!("Right threshold helper is missing its threshold order")
                    })?;
                    Ok(CrossSectionGraph::single_th_prefactor_helper_params(
                        order as u8,
                        true,
                    ))
                })
            })
            .collect::<Result<Vec<_>>>()?,
        iterated_helpers: export_iterated_collection(
            &evaluators.threshold_helpers.iterated,
            |helpers| {
                export_generic_evaluator_map(helpers, |index| {
                    let left_order = index.left_threshold_order.ok_or_else(|| {
                        eyre!("Iterated threshold helper is missing its left threshold order")
                    })?;
                    let right_order = index.right_threshold_order.ok_or_else(|| {
                        eyre!("Iterated threshold helper is missing its right threshold order")
                    })?;
                    Ok(CrossSectionGraph::iterated_th_prefactor_helper_params(
                        left_order as u8,
                        right_order as u8,
                    ))
                })
            },
        )?,
        pass_two_evaluator: evaluators
            .residue_from_e_surface_evaluators
            .iter()
            .enumerate()
            .map(|(order, evaluator)| {
                let params = params_for_derivative_order((order + 1) as u8);
                export_generic_evaluator(evaluator, Some(&params))
            })
            .collect::<Result<Vec<_>>>()?,
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
                let orientation_start = term.param_builder.pairs.orientations.value_range.start;
                let residue_map_id_start =
                    term.param_builder.pairs.residue_map_id.value_range.start;
                let multiplicative_offset = 1;
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
                    .map(|stacks| {
                        export_evaluator_map(
                            stacks,
                            orientation_start,
                            residue_map_id_start,
                            multiplicative_offset,
                        )
                    })
                    .collect::<Result<Vec<_>>>()?;

                let counterterms = term
                    .counterterm
                    .evaluators
                    .iter()
                    .map(|evaluators| {
                        export_counterterm(
                            evaluators,
                            orientation_start,
                            residue_map_id_start,
                            multiplicative_offset,
                        )
                    })
                    .collect::<Result<Vec<_>>>()?;

                Ok(StandaloneCrossSectionGraphTermArchive {
                    graph_name: term.graph.name.clone(),
                    orientations,
                    param_builder_params,
                    fn_map_entries,
                    cut_group_integrands,
                    counterterms,
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
