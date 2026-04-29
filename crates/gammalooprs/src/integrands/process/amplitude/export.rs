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
                STANDALONE_EVALUATORS_VERSION, StandaloneComplexInput, StandaloneEvaluatorArchive,
                StandaloneEvaluatorStackArchive, StandaloneGenericEvaluatorArchive,
                StandaloneGraphTermArchive,
            },
        },
    },
    momentum::ThreeMomentum,
    momentum::sample::{LoopMomenta, MomentumSample},
    processes::{
        StandaloneDataFormat, StandaloneExportMode, StandaloneExportSettings,
        StandaloneNumericTarget,
    },
    utils::{ArbPrec, F, FloatLike, f128},
};

const STANDALONE_DATA_FILE: &str = "standalone_evaluators";
const STANDALONE_RUST_SCRIPT_FILE: &str = "standalone_evaluators_rust.rs";

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

fn export_evaluator_stack<T: ExportAtomTo>(
    evaluator_stack: &crate::integrands::process::evaluators::EvaluatorStack,
    start: usize,
    override_pos: usize,
    mult_offset: usize,
    representative_input: Vec<StandaloneComplexInput>,
) -> Result<StandaloneEvaluatorStackArchive<T>> {
    Ok(StandaloneEvaluatorStackArchive {
        start,
        override_pos,
        mult_offset,
        representative_input,
        single_parametric: export_generic_evaluator(&evaluator_stack.single_parametric)?,
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
}

fn standalone_rust_script() -> String {
    include_str!("standalone_template.rs").to_string()
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

trait StandaloneRepresentativeInputPrecision: FloatLike + GenericEvaluatorFloat + Default {
    fn sample_from_base(sample: &MomentumSample<f64>) -> MomentumSample<Self>;
}

impl StandaloneRepresentativeInputPrecision for f64 {
    fn sample_from_base(sample: &MomentumSample<f64>) -> MomentumSample<Self> {
        sample.clone()
    }
}

impl StandaloneRepresentativeInputPrecision for f128 {
    fn sample_from_base(sample: &MomentumSample<f64>) -> MomentumSample<Self> {
        sample.higher_precision()
    }
}

impl StandaloneRepresentativeInputPrecision for ArbPrec {
    fn sample_from_base(sample: &MomentumSample<f64>) -> MomentumSample<Self> {
        sample.higher_precision().higher_precision()
    }
}

fn serialize_representative_input<T: FloatLike>(
    values: &[spenso::algebra::complex::Complex<F<T>>],
) -> Vec<StandaloneComplexInput> {
    values
        .iter()
        .map(|value| StandaloneComplexInput {
            re: value.re.to_string(),
            im: value.im.to_string(),
        })
        .collect()
}

impl AmplitudeIntegrand {
    fn create_archive<T: ExportAtomTo, S>(
        &self,
        symbolica_state: S,
        numeric_target: StandaloneNumericTarget,
    ) -> Result<StandaloneEvaluatorArchive<S, T>> {
        let sample_inputs = match numeric_target {
            StandaloneNumericTarget::Double => self.representative_input_for::<f64>()?,
            StandaloneNumericTarget::Quad => self.representative_input_for::<f128>()?,
            StandaloneNumericTarget::Arb => self.representative_input_for::<ArbPrec>()?,
        };

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
                            counterterm
                                .evaluator_stacks
                                .iter()
                                .map(|evaluator_stack| {
                                    export_evaluator_stack(
                                        evaluator_stack,
                                        start,
                                        override_pos,
                                        mult_offset,
                                        representative_input.clone(),
                                    )
                                })
                                .collect::<Result<Vec<_>>>()
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
            numeric_target,
            symbolica_state,
            graph_terms,
        })
    }

    #[allow(clippy::type_complexity)]
    fn representative_input_for<T>(
        &self,
    ) -> Result<Vec<(Vec<StandaloneComplexInput>, usize, usize, usize)>>
    where
        T: StandaloneRepresentativeInputPrecision,
    {
        let hel = self.settings.kinematics.externals.get_helicities();
        let mut rng = rand::rng();

        let mut base_mom_samples = Vec::new();
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
                base_mom_samples.push(momentum_sample);
            }
        }
        let mut inputs = Vec::new();

        for (t, base_momentum_sample) in self.data.graph_terms.iter().zip(&base_mom_samples) {
            let mut p_build = t.param_builder.clone();
            let momentum_sample = T::sample_from_base(base_momentum_sample);
            let input = T::get_parameters(
                &mut p_build,
                (false, false),
                &t.graph,
                &momentum_sample,
                hel,
                &self.settings.additional_params(),
                None,
                None,
                None,
            );

            inputs.push((
                serialize_representative_input(input.as_slice()),
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
                    self.create_archive(symbolica_state, settings.precision)?;
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
                let standalone: StandaloneEvaluatorArchive<(), String> =
                    self.create_archive((), settings.precision)?;
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
                make_script_executable(&script_path)?;
            }
        }

        Ok(())
    }
}
