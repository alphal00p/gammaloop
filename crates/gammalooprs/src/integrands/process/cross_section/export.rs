use std::{fs, path::Path};

use color_eyre::Result;
use eyre::{Context, eyre};
use linnet::half_edge::involution::Orientation;
use symbolica::state::State;

use super::{
    CrossSectionIntegrand,
    load::{
        STANDALONE_EVALUATORS_VERSION, StandaloneCountertermArchive, StandaloneCrossSectionArchive,
        StandaloneCrossSectionGraphTermArchive, StandaloneEvaluatorStackArchive,
        StandaloneGenericEvaluatorArchive, StandaloneIteratedCollectionArchive,
    },
};
use crate::{
    integrands::process::{GenericEvaluator, amplitude::export::ExportAtomTo},
    processes::{
        IteratedCtCollection, StandaloneDataFormat, StandaloneExportMode, StandaloneExportSettings,
        StandaloneNumericTarget,
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
        override_pos: 0,
        mult_offset: 0,
    })
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
        num_left_thresholds: collection.num_left_thresholds(),
    })
}

fn export_counterterm<T: ExportAtomTo>(
    evaluators: &LUCounterTermEvaluators,
) -> Result<StandaloneCountertermArchive<T>> {
    Ok(StandaloneCountertermArchive {
        left_thresholds_evaluator: evaluators
            .left_thresholds_evaluator
            .iter()
            .map(|stacks| {
                stacks
                    .iter()
                    .map(export_evaluator_stack)
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?,
        right_thresholds_evaluator: evaluators
            .right_thresholds_evaluator
            .iter()
            .map(|stacks| {
                stacks
                    .iter()
                    .map(export_evaluator_stack)
                    .collect::<Result<Vec<_>>>()
            })
            .collect::<Result<Vec<_>>>()?,
        iterated_evaluator: export_iterated_collection(&evaluators.iterated_evaluator, |stacks| {
            stacks
                .iter()
                .map(export_evaluator_stack)
                .collect::<Result<Vec<_>>>()
        })?,
        pass_two_evaluator: evaluators
            .residue_from_e_surface_evaluators
            .iter()
            .map(export_generic_evaluator)
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

                let raised_cut_integrands = term
                    .integrand
                    .iter()
                    .map(|integrands| {
                        integrands
                            .iter()
                            .map(export_evaluator_stack)
                            .collect::<Result<Vec<_>>>()
                    })
                    .collect::<Result<Vec<_>>>()?;

                let counterterms = term
                    .counterterm
                    .evaluators
                    .iter()
                    .map(export_counterterm)
                    .collect::<Result<Vec<_>>>()?;

                Ok(StandaloneCrossSectionGraphTermArchive {
                    graph_name: term.graph.name.clone(),
                    orientations,
                    param_builder_params,
                    fn_map_entries,
                    raised_cut_integrands,
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
