use color_eyre::Result;
use eyre::eyre;
use gammalooprs::{
    integrands::{
        evaluation::{
            BatchSampleEvaluationResult, PreciseBatchSampleEvaluationResult,
            PreciseSampleEvaluationResult, PreciseSingleSampleEvaluationResult,
            SampleEvaluationResult, SingleSampleEvaluationResult,
        },
        process::MomentumSpaceEvaluationInput,
    },
    momentum::ThreeMomentum,
    observables::ObservableSnapshotBundle,
    utils::F,
};
use ndarray::ArrayView2;
use spenso::algebra::complex::Complex;
use symbolica::numerical_integration::Sample;

use crate::state::{ProcessRef, State};

#[derive(Debug, Clone)]
pub struct EvaluateSamples<'a> {
    pub process_id: Option<usize>,
    pub integrand_name: Option<String>,
    pub use_arb_prec: bool,
    pub minimal_output: bool,
    pub momentum_space: bool,
    pub points: ArrayView2<'a, f64>,
    pub discrete_dims: Option<ArrayView2<'a, usize>>,
    pub graph_names: Option<Vec<Option<String>>>,
    pub orientations: Option<Vec<Option<usize>>>,
}

impl<'a> EvaluateSamples<'a> {
    pub fn run(&self, state: &mut State) -> Result<BatchSampleEvaluationResult> {
        let process_ref = self.process_id.map(ProcessRef::Id);
        let (process_id, integrand_name) =
            state.find_integrand_ref(process_ref.as_ref(), self.integrand_name.as_ref())?;
        let model = state.resolve_model_for_integrand(process_id, &integrand_name)?;

        let integrand = state
            .process_list
            .get_integrand_mut(process_id, &integrand_name)?;
        integrand.warm_up(&model)?;

        let batch_len = self.points.nrows();
        let graph_names =
            normalize_optional_per_sample(self.graph_names.as_ref(), batch_len, "graph_names")?;
        let orientations =
            normalize_optional_per_sample(self.orientations.as_ref(), batch_len, "orientations")?;

        if !self.momentum_space
            && (graph_names.iter().any(Option::is_some) || orientations.iter().any(Option::is_some))
        {
            return Err(eyre!(
                "Graph and orientation selection are only supported in momentum-space evaluation."
            ));
        }

        if let Some(discrete_dims) = &self.discrete_dims {
            if discrete_dims.nrows() != batch_len {
                return Err(eyre!(
                    "Expected {} rows in discrete_dims, got {}.",
                    batch_len,
                    discrete_dims.nrows()
                ));
            }
        }

        let samples = if self.momentum_space {
            let inputs = build_momentum_inputs(
                integrand,
                &self.points,
                self.discrete_dims.as_ref(),
                &graph_names,
                &orientations,
            )?;
            integrand
                .evaluate_momentum_configurations_raw(&model, &inputs, self.use_arb_prec)?
                .samples
        } else {
            let samples =
                build_havana_samples(integrand, &self.points, self.discrete_dims.as_ref())?;
            integrand
                .evaluate_samples_raw(&model, &samples, 1, self.use_arb_prec, Complex::new_zero())?
                .samples
        };

        let observables = integrand
            .observable_snapshot_bundle()
            .unwrap_or_else(ObservableSnapshotBundle::default);

        Ok(BatchSampleEvaluationResult {
            samples: samples
                .into_iter()
                .map(|evaluation| SampleEvaluationResult {
                    evaluation: evaluation.into_output(self.minimal_output),
                })
                .collect(),
            observables,
        })
    }
}

#[derive(Debug, Clone)]
pub struct EvaluateSamplesPrecise<'a> {
    pub process_id: Option<usize>,
    pub integrand_name: Option<String>,
    pub use_arb_prec: bool,
    pub minimal_output: bool,
    pub momentum_space: bool,
    pub points: ArrayView2<'a, f64>,
    pub discrete_dims: Option<ArrayView2<'a, usize>>,
    pub graph_names: Option<Vec<Option<String>>>,
    pub orientations: Option<Vec<Option<usize>>>,
}

impl<'a> EvaluateSamplesPrecise<'a> {
    pub fn run(&self, state: &mut State) -> Result<PreciseBatchSampleEvaluationResult> {
        let process_ref = self.process_id.map(ProcessRef::Id);
        let (process_id, integrand_name) =
            state.find_integrand_ref(process_ref.as_ref(), self.integrand_name.as_ref())?;
        let model = state.resolve_model_for_integrand(process_id, &integrand_name)?;

        let integrand = state
            .process_list
            .get_integrand_mut(process_id, &integrand_name)?;
        integrand.warm_up(&model)?;

        let batch_len = self.points.nrows();
        let graph_names =
            normalize_optional_per_sample(self.graph_names.as_ref(), batch_len, "graph_names")?;
        let orientations =
            normalize_optional_per_sample(self.orientations.as_ref(), batch_len, "orientations")?;

        if !self.momentum_space
            && (graph_names.iter().any(Option::is_some) || orientations.iter().any(Option::is_some))
        {
            return Err(eyre!(
                "Graph and orientation selection are only supported in momentum-space evaluation."
            ));
        }

        if let Some(discrete_dims) = &self.discrete_dims {
            if discrete_dims.nrows() != batch_len {
                return Err(eyre!(
                    "Expected {} rows in discrete_dims, got {}.",
                    batch_len,
                    discrete_dims.nrows()
                ));
            }
        }

        let mut raw_samples = if self.momentum_space {
            let inputs = build_momentum_inputs(
                integrand,
                &self.points,
                self.discrete_dims.as_ref(),
                &graph_names,
                &orientations,
            )?;
            integrand
                .evaluate_momentum_configurations_precise_raw(&model, &inputs, self.use_arb_prec)?
                .samples
        } else {
            let samples =
                build_havana_samples(integrand, &self.points, self.discrete_dims.as_ref())?;
            integrand
                .evaluate_samples_precise_raw(
                    &model,
                    &samples,
                    self.use_arb_prec,
                    Complex::new_zero(),
                )?
                .samples
        };

        let observables = merge_precise_observables(integrand, &mut raw_samples)?;

        Ok(PreciseBatchSampleEvaluationResult {
            samples: raw_samples
                .into_iter()
                .map(|evaluation| PreciseSampleEvaluationResult {
                    evaluation: evaluation.into_output(self.minimal_output),
                })
                .collect(),
            observables,
        })
    }
}

pub fn evaluate_samples_precise<'a>(
    state: &mut State,
    request: &EvaluateSamplesPrecise<'a>,
) -> Result<PreciseBatchSampleEvaluationResult> {
    request.run(state)
}

pub fn evaluate_sample_precise<'a>(
    state: &mut State,
    request: &EvaluateSamplesPrecise<'a>,
) -> Result<PreciseSingleSampleEvaluationResult> {
    let mut results = request.run(state)?;
    if results.samples.len() != 1 {
        return Err(eyre!(
            "evaluate_sample_precise expects exactly one sample, got {}.",
            results.samples.len()
        ));
    }
    Ok(PreciseSingleSampleEvaluationResult {
        sample: results.samples.remove(0),
        observables: results.observables,
    })
}

pub fn evaluate_samples<'a>(
    state: &mut State,
    request: &EvaluateSamples<'a>,
) -> Result<BatchSampleEvaluationResult> {
    request.run(state)
}

pub fn evaluate_sample<'a>(
    state: &mut State,
    request: &EvaluateSamples<'a>,
) -> Result<SingleSampleEvaluationResult> {
    let mut results = request.run(state)?;
    if results.samples.len() != 1 {
        return Err(eyre!(
            "evaluate_sample expects exactly one sample, got {}.",
            results.samples.len()
        ));
    }
    Ok(SingleSampleEvaluationResult {
        sample: results.samples.remove(0),
        observables: results.observables,
    })
}

fn normalize_optional_per_sample<T: Clone>(
    values: Option<&Vec<Option<T>>>,
    batch_len: usize,
    label: &str,
) -> Result<Vec<Option<T>>> {
    match values {
        Some(values) => {
            if values.len() != batch_len {
                return Err(eyre!(
                    "Expected {} entries in {}, got {}.",
                    batch_len,
                    label,
                    values.len()
                ));
            }
            Ok(values.clone())
        }
        None => Ok(vec![None; batch_len]),
    }
}

fn build_havana_sample(
    integrand: &gammalooprs::integrands::process::ProcessIntegrand,
    point: &[f64],
    discrete_dim: &[usize],
) -> Result<Sample<F<f64>>> {
    let expected_dimension = integrand.expected_x_space_dimension(discrete_dim)?;
    if point.len() != expected_dimension {
        return Err(eyre!(
            "Expected {} x-space coordinates for this integrand selection, got {}.",
            expected_dimension,
            point.len()
        ));
    }

    Ok(havana_sample(
        point.iter().map(|&x| F(x)).collect::<Vec<F<f64>>>(),
        discrete_dim,
    ))
}

fn build_havana_samples(
    integrand: &gammalooprs::integrands::process::ProcessIntegrand,
    points: &ArrayView2<'_, f64>,
    discrete_dims: Option<&ArrayView2<'_, usize>>,
) -> Result<Vec<Sample<F<f64>>>> {
    (0..points.nrows())
        .map(|sample_index| {
            let discrete_dim = discrete_dims
                .as_ref()
                .map(|dims| dims.row(sample_index).iter().copied().collect())
                .unwrap_or_else(Vec::new);
            build_havana_sample(
                integrand,
                points
                    .row(sample_index)
                    .as_slice()
                    .expect("point rows are contiguous"),
                &discrete_dim,
            )
        })
        .collect()
}

fn build_momentum_input(
    integrand: &gammalooprs::integrands::process::ProcessIntegrand,
    point: &[f64],
    discrete_dim: &[usize],
    graph_name: Option<&str>,
    orientation: Option<usize>,
) -> Result<MomentumSpaceEvaluationInput> {
    if point.len() % 3 != 0 {
        return Err(eyre!(
            "Momentum-space evaluation expects a multiple of 3 coordinates, got {}.",
            point.len()
        ));
    }

    let loop_momenta = point
        .chunks_exact(3)
        .map(|coords| ThreeMomentum {
            px: F(coords[0]),
            py: F(coords[1]),
            pz: F(coords[2]),
        })
        .collect();

    match &integrand.get_settings().sampling {
        gammalooprs::settings::runtime::SamplingSettings::Default(_)
        | gammalooprs::settings::runtime::SamplingSettings::MultiChanneling(_) => {
            if !discrete_dim.is_empty() {
                return Err(eyre!(
                    "This integrand does not use discrete graph sampling; expected no discrete dimensions, got {:?}.",
                    discrete_dim
                ));
            }

            if orientation.is_some() && graph_name.is_none() {
                return Err(eyre!(
                    "A specific orientation requires a specific graph in momentum-space evaluation."
                ));
            }

            let graph_id = graph_name
                .map(|graph_name| {
                    integrand.find_graph_id_by_name(graph_name).ok_or_else(|| {
                        eyre!(
                            "Unknown graph '{}' in momentum-space evaluation.",
                            graph_name
                        )
                    })
                })
                .transpose()?;

            if let (Some(graph_id), Some(orientation)) = (graph_id, orientation) {
                let orientation_count =
                    integrand.graph_orientation_count(graph_id).ok_or_else(|| {
                        eyre!(
                            "Graph id {} is invalid for momentum-space evaluation.",
                            graph_id
                        )
                    })?;
                if orientation >= orientation_count {
                    return Err(eyre!(
                        "Orientation {} is out of range for graph '{}'; the graph has {} orientations.",
                        orientation,
                        graph_name.expect("graph name should exist when graph_id exists"),
                        orientation_count
                    ));
                }
            }

            Ok(MomentumSpaceEvaluationInput {
                loop_momenta,
                graph_id,
                group_id: None,
                orientation,
                channel_id: None,
            })
        }
        gammalooprs::settings::runtime::SamplingSettings::DiscreteGraphs(settings) => {
            let (group_from_dims, orientation_from_dims, channel_id) = if discrete_dim.is_empty() {
                (None, None, None)
            } else {
                integrand.resolve_discrete_selection(discrete_dim)?
            };

            let group_from_name = graph_name
                .map(|graph_name| integrand.resolve_group_id_by_master_name(graph_name))
                .transpose()?;

            let group_id = match (group_from_dims, group_from_name) {
                (Some(group_id), Some(group_from_name)) if group_id != group_from_name => {
                    return Err(eyre!(
                        "Discrete dimensions select graph group {}, but graph '{}' resolves to group {}.",
                        group_id.0,
                        graph_name.expect("graph name must exist"),
                        group_from_name.0
                    ));
                }
                (Some(group_id), _) => Some(group_id),
                (None, Some(group_id)) => Some(group_id),
                (None, None) => None,
            };

            let orientation = match (orientation_from_dims, orientation) {
                (Some(from_dims), Some(from_arg)) if from_dims != from_arg => {
                    return Err(eyre!(
                        "Discrete dimensions select orientation {}, but the explicit momentum-space orientation is {}.",
                        from_dims,
                        from_arg
                    ));
                }
                (Some(from_dims), _) => Some(from_dims),
                (None, from_arg) => from_arg,
            };

            if settings.sample_orientations && group_id.is_some() && orientation.is_none() {
                return Err(eyre!(
                    "Momentum-space evaluation for this integrand requires an orientation selection."
                ));
            }

            Ok(MomentumSpaceEvaluationInput {
                loop_momenta,
                graph_id: None,
                group_id,
                orientation,
                channel_id,
            })
        }
    }
}

fn build_momentum_inputs(
    integrand: &gammalooprs::integrands::process::ProcessIntegrand,
    points: &ArrayView2<'_, f64>,
    discrete_dims: Option<&ArrayView2<'_, usize>>,
    graph_names: &[Option<String>],
    orientations: &[Option<usize>],
) -> Result<Vec<MomentumSpaceEvaluationInput>> {
    (0..points.nrows())
        .map(|sample_index| {
            let discrete_dim = discrete_dims
                .as_ref()
                .map(|dims| dims.row(sample_index).iter().copied().collect())
                .unwrap_or_else(Vec::new);
            build_momentum_input(
                integrand,
                points
                    .row(sample_index)
                    .as_slice()
                    .expect("point rows are contiguous"),
                &discrete_dim,
                graph_names[sample_index].as_deref(),
                orientations[sample_index],
            )
        })
        .collect()
}

fn clear_precise_event_groups(
    result: &mut gammalooprs::integrands::evaluation::PreciseEvaluationResult,
) {
    match result {
        gammalooprs::integrands::evaluation::PreciseEvaluationResult::Double(result) => {
            result.event_groups.clear()
        }
        gammalooprs::integrands::evaluation::PreciseEvaluationResult::Quad(result) => {
            result.event_groups.clear()
        }
        gammalooprs::integrands::evaluation::PreciseEvaluationResult::Arb(result) => {
            result.event_groups.clear()
        }
    }
}

fn merge_precise_observables(
    integrand: &gammalooprs::integrands::process::ProcessIntegrand,
    results: &mut [gammalooprs::integrands::evaluation::PreciseEvaluationResult],
) -> Result<ObservableSnapshotBundle> {
    let mut observables = ObservableSnapshotBundle::default();
    let mut have_observables = false;

    for evaluation in results.iter_mut() {
        let sample_observables = integrand
            .build_observable_snapshots_for_precise_result(evaluation)
            .unwrap_or_else(ObservableSnapshotBundle::default);
        if !sample_observables.histograms.is_empty() {
            if have_observables {
                observables.merge_in_place(&sample_observables)?;
            } else {
                observables = sample_observables;
                have_observables = true;
            }
        }
        if !integrand.get_settings().general.generate_events {
            clear_precise_event_groups(evaluation);
        }
    }

    Ok(observables)
}

fn havana_sample(cont: Vec<F<f64>>, discrete_dimensions: &[usize]) -> Sample<F<f64>> {
    let mut sample = Sample::Continuous(F(1.0), cont);

    for &discrete_dimension in discrete_dimensions.iter().rev() {
        sample = Sample::Discrete(F(1.0), discrete_dimension, Some(Box::new(sample)));
    }

    sample
}
