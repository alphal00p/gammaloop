use std::{
    fs::{self},
    path::Path,
};

use color_eyre::Result;
use eyre::{Context, eyre};
use linnet::half_edge::involution::Orientation;
use rand::Rng;
use symbolica::{
    atom::{Atom, AtomCore},
    printer::PrintOptions,
    state::State,
};

use crate::{
    integrands::process::{
        GenericEvaluator, GenericEvaluatorFloat, ProcessIntegrandImpl,
        amplitude::{
            AmplitudeIntegrand,
            load::{
                STANDALONE_EVALUATORS_VERSION, StandaloneEvaluatorArchive,
                StandaloneEvaluatorStackArchive, StandaloneGenericEvaluatorArchive,
                StandaloneGraphTermArchive,
            },
        },
    },
    momentum::ThreeMomentum,
    momentum::sample::{LoopMomenta, MomentumSample},
    processes::{StandaloneDataFormat, StandaloneExportMode, StandaloneExportSettings},
    utils::F,
};

const STANDALONE_DATA_FILE: &str = "standalone_evaluators";
const STANDALONE_RUST_SCRIPT_FILE: &str = "standalone_evaluators_rust.rs";

fn atom_to_bytes(atom: &Atom) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    atom.as_view()
        .write(&mut out)
        .with_context(|| "Failed to write atom to standalone binary stream")?;
    Ok(out)
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
        .map(|a| T::export_atom_to(a))
        .collect::<Result<Vec<_>>>()?;

    let additional_fn_map_entries = evaluator
        .fn_map_entries
        .iter()
        .map(|e| e.archive())
        .collect::<Result<_>>()?;

    Ok(StandaloneGenericEvaluatorArchive {
        exprs,
        additional_fn_map_entries,
        dual_shape: evaluator.dual_shape.clone(),
    })
}

fn standalone_rust_script() -> String {
    let mut script = include_str!("load.rs").to_string();
    if let Some(rest) = script.strip_prefix("//#!/usr/bin/env -S rust-script\n") {
        script = format!("#!/usr/bin/env -S rust-script\n{rest}");
    }

    script
}

pub trait CreateArchive<T, S> {
    fn create_archive(&self, state: S) -> Result<StandaloneEvaluatorArchive<S, T>>;
}

pub trait ExportAtomTo: Sized {
    fn export_atom_to(atom: &Atom) -> Result<Self>;
}

impl ExportAtomTo for Vec<u8> {
    fn export_atom_to(atom: &Atom) -> Result<Self> {
        atom_to_bytes(atom)
    }
}

impl ExportAtomTo for String {
    fn export_atom_to(atom: &Atom) -> Result<Self> {
        Ok(atom.printer(PrintOptions::file()).to_string())
    }
}

impl<T: ExportAtomTo, S> CreateArchive<T, S> for AmplitudeIntegrand {
    fn create_archive(&self, symbolica_state: S) -> Result<StandaloneEvaluatorArchive<S, T>> {
        let sample_inputs = self.representative_input()?;

        let graph_terms = self
            .data
            .graph_terms
            .iter()
            .zip(sample_inputs)
            .map(
                |(term, (representative_input, mult_offset, start, override_pos))| {
                    let param_builder_params = (&term.param_builder.pairs)
                        .into_iter()
                        .flat_map(|pair| pair.params.iter())
                        .map(|a| T::export_atom_to(a))
                        .collect::<Result<Vec<_>>>()?;

                    let orientations: Vec<Vec<_>> = term
                        .orientations
                        .iter()
                        .map(|a| {
                            a.iter()
                                .map(|(_, o)| match o {
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

                    let original_integrand = StandaloneEvaluatorStackArchive {
                        start,
                        mult_offset,
                        override_pos,
                        representative_input: representative_input.clone(),
                        single_parametric: export_generic_evaluator(
                            &term.original_integrand.single_parametric,
                        )?,
                        iterative: term
                            .original_integrand
                            .iterative
                            .as_ref()
                            .map(|(eval, _)| export_generic_evaluator(eval))
                            .transpose()?,
                        summed_function_map: term
                            .original_integrand
                            .summed_function_map
                            .as_ref()
                            .map(export_generic_evaluator)
                            .transpose()?,
                        summed: term
                            .original_integrand
                            .summed
                            .as_ref()
                            .map(export_generic_evaluator)
                            .transpose()?,
                    };

                    let threshold_counterterms = term
                        .threshold_counterterm
                        .evaluators
                        .iter()
                        .map(|counterterm| {
                            let evaluator_stack = counterterm.evaluator_stack.borrow();
                            Ok(StandaloneEvaluatorStackArchive {
                                start,
                                override_pos,
                                mult_offset,
                                representative_input: representative_input.clone(),
                                single_parametric: export_generic_evaluator(
                                    &evaluator_stack.single_parametric,
                                )?,
                                iterative: evaluator_stack
                                    .iterative
                                    .as_ref()
                                    .map(|(eval, _)| export_generic_evaluator(eval))
                                    .transpose()?,
                                summed_function_map: evaluator_stack
                                    .summed_function_map
                                    .as_ref()
                                    .map(export_generic_evaluator)
                                    .transpose()?,
                                summed: evaluator_stack
                                    .summed
                                    .as_ref()
                                    .map(export_generic_evaluator)
                                    .transpose()?,
                            })
                        })
                        .collect::<Result<Vec<_>>>()?;

                    Ok(StandaloneGraphTermArchive {
                        graph_name: term.graph.name.clone(),
                        orientations,
                        param_builder_params,
                        fn_map_entries,
                        original_integrand,
                        threshold_counterterms,
                    })
                },
            )
            .collect::<Result<Vec<_>>>()?;
        Ok(StandaloneEvaluatorArchive {
            version: STANDALONE_EVALUATORS_VERSION,
            symbolica_state,
            graph_terms,
        })
    }
}

impl AmplitudeIntegrand {
    #[allow(clippy::type_complexity)]
    pub(crate) fn representative_input(
        &self,
    ) -> Result<
        Vec<(
            Vec<symbolica::domains::float::Complex<f64>>,
            usize,
            usize,
            usize,
        )>,
    > {
        let hel = self.settings.kinematics.externals.get_helicities();
        let mut rng = rand::rng();

        let mut mom_samples = Vec::new();
        {
            let dependent_momenta_constructor = self.get_dependent_momenta_constructor();
            for t in self.data.graph_terms.iter() {
                let loop_moms = LoopMomenta::from_iter(
                    (0..t.graph.loop_momentum_basis.loop_edges.len()).map(|_| ThreeMomentum {
                        px: F(rng.random_range(-1.0..1.0)),
                        py: F(rng.random_range(-1.0..1.0)),
                        pz: F(rng.random_range(-1.0..1.0)),
                    }),
                );
                let momentum_sample = MomentumSample::new(
                    loop_moms,
                    0,
                    &self.settings.kinematics.externals,
                    0,
                    F(1.0),
                    dependent_momenta_constructor,
                    None,
                )?;
                mom_samples.push(momentum_sample);
            }
        }
        let mut inputs = Vec::new();

        for (t, momentum_sample) in self.data.graph_terms.iter().zip(&mom_samples) {
            let mut p_build = t.param_builder.clone();
            let input = f64::get_parameters(
                &mut p_build,
                (false, false),
                &t.graph,
                momentum_sample,
                hel,
                &self.settings.additional_params(),
                None,
                None,
                None,
            );

            inputs.push((
                input
                    .as_slice()
                    .iter()
                    .map(|c| symbolica::domains::float::Complex::new(c.re.0, c.im.0))
                    .collect(),
                input.multiplicative_offset,
                input.orientations_start,
                input.override_pos,
            ));
        }

        Ok(inputs)
    }

    // pub(crate) fn create_archive<

    pub(crate) fn export_standalone(
        &self,
        path: impl AsRef<Path>,
        settings: &StandaloneExportSettings,
    ) -> Result<()> {
        let mut symbolica_state = Vec::new();
        State::export(&mut symbolica_state)
            .with_context(|| "Failed to export Symbolica state for standalone evaluators")?;

        let mut standalone_path = path.as_ref().join(STANDALONE_DATA_FILE);

        match settings.format {
            StandaloneDataFormat::Binary => {
                let standalone: StandaloneEvaluatorArchive<Vec<u8>, Vec<u8>> =
                    self.create_archive(symbolica_state)?;
                let binary = bincode::encode_to_vec(&standalone, bincode::config::standard())?;
                standalone_path.add_extension("bin");
                fs::write(&standalone_path, binary).with_context(|| {
                    format!(
                        "Failed writing standalone evaluator binary to {}",
                        standalone_path.display()
                    )
                })?;
            }
            StandaloneDataFormat::Json => {
                let standalone: StandaloneEvaluatorArchive<(), String> = self.create_archive(())?;
                let json = serde_json::to_vec_pretty(&standalone)?;
                standalone_path.add_extension("json");
                fs::write(&standalone_path, json).with_context(|| {
                    format!(
                        "Failed writing standalone evaluator JSON to {}",
                        standalone_path.display()
                    )
                })?
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
                        "Failed writing standalone rust-script loader to {}",
                        script_path.display()
                    )
                })?;
            }
        }

        Ok(())
    }
}
