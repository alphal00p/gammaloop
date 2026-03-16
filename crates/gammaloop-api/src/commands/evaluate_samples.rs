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
        HasIntegrand,
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
    pub force_radius: bool,
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

        let integrand = state
            .process_list
            .get_integrand_mut(process_id, &integrand_name)?;
        integrand.warm_up(&state.model)?;

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

        let mut samples = Vec::with_capacity(batch_len);
        let mut observables = ObservableSnapshotBundle::default();
        let mut have_observables = false;
        for sample_index in 0..batch_len {
            let point = self.points.row(sample_index);
            let graph_name = graph_names[sample_index].as_deref();
            let orientation = orientations[sample_index];

            let mut evaluation = if self.momentum_space {
                evaluate_momentum_point(
                    integrand,
                    &state.model,
                    point.as_slice().expect("point rows are contiguous"),
                    graph_name,
                    orientation,
                    self.use_arb_prec,
                )?
            } else {
                let discrete_dim = self
                    .discrete_dims
                    .as_ref()
                    .map(|dims| dims.row(sample_index).iter().copied().collect())
                    .unwrap_or_else(Vec::new);
                evaluate_x_point(
                    integrand,
                    &state.model,
                    point.as_slice().expect("point rows are contiguous"),
                    &discrete_dim,
                    self.force_radius,
                    self.use_arb_prec,
                )?
            };

            let sample_observables = integrand
                .build_observable_snapshots_for_result(&evaluation)
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
                evaluation.event_groups.clear();
            }
            samples.push(SampleEvaluationResult {
                evaluation: evaluation.into_output(self.minimal_output),
            });
        }

        Ok(BatchSampleEvaluationResult {
            samples,
            observables,
        })
    }
}

#[derive(Debug, Clone)]
pub struct EvaluateSamplesPrecise<'a> {
    pub process_id: Option<usize>,
    pub integrand_name: Option<String>,
    pub use_arb_prec: bool,
    pub force_radius: bool,
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

        let integrand = state
            .process_list
            .get_integrand_mut(process_id, &integrand_name)?;
        integrand.warm_up(&state.model)?;

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

        let mut samples = Vec::with_capacity(batch_len);
        let mut observables = ObservableSnapshotBundle::default();
        let mut have_observables = false;
        for sample_index in 0..batch_len {
            let point = self.points.row(sample_index);
            let graph_name = graph_names[sample_index].as_deref();
            let orientation = orientations[sample_index];

            let mut evaluation = if self.momentum_space {
                evaluate_momentum_point_precise(
                    integrand,
                    &state.model,
                    point.as_slice().expect("point rows are contiguous"),
                    graph_name,
                    orientation,
                    self.use_arb_prec,
                )?
            } else {
                let discrete_dim = self
                    .discrete_dims
                    .as_ref()
                    .map(|dims| dims.row(sample_index).iter().copied().collect())
                    .unwrap_or_else(Vec::new);
                evaluate_x_point_precise(
                    integrand,
                    &state.model,
                    point.as_slice().expect("point rows are contiguous"),
                    &discrete_dim,
                    self.force_radius,
                    self.use_arb_prec,
                )?
            };

            let sample_observables = integrand
                .build_observable_snapshots_for_precise_result(&evaluation)
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
                clear_precise_event_groups(&mut evaluation);
            }
            samples.push(PreciseSampleEvaluationResult {
                evaluation: evaluation.into_output(self.minimal_output),
            });
        }

        Ok(PreciseBatchSampleEvaluationResult {
            samples,
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

fn evaluate_x_point(
    integrand: &mut gammalooprs::integrands::process::ProcessIntegrand,
    model: &gammalooprs::model::Model,
    point: &[f64],
    discrete_dim: &[usize],
    mut force_radius: bool,
    use_arb_prec: bool,
) -> Result<gammalooprs::integrands::evaluation::EvaluationResult> {
    if integrand.get_n_dim() as isize == point.len() as isize - 1 {
        force_radius = true;
    }

    let point = point.iter().map(|&x| F(x)).collect::<Vec<F<f64>>>();
    let continuous_sample = if force_radius {
        point
            .get(1..)
            .ok_or_else(|| {
                eyre!("Force-radius x-space evaluation requires at least one radial coordinate.")
            })?
            .to_vec()
    } else {
        point
    };
    let sample = havana_sample(continuous_sample, discrete_dim);

    integrand.evaluate_sample(&sample, model, F(1.0), 1, use_arb_prec, Complex::new_zero())
}

fn evaluate_momentum_point(
    integrand: &mut gammalooprs::integrands::process::ProcessIntegrand,
    model: &gammalooprs::model::Model,
    point: &[f64],
    graph_name: Option<&str>,
    orientation: Option<usize>,
    use_arb_prec: bool,
) -> Result<gammalooprs::integrands::evaluation::EvaluationResult> {
    if orientation.is_some() && graph_name.is_none() {
        return Err(eyre!(
            "A specific orientation requires a specific graph in momentum-space evaluation."
        ));
    }

    if point.len() % 3 != 0 {
        return Err(eyre!(
            "Momentum-space evaluation expects a multiple of 3 coordinates, got {}.",
            point.len()
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
        let orientation_count = integrand.graph_orientation_count(graph_id).ok_or_else(|| {
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

    let loop_momenta = point
        .chunks_exact(3)
        .map(|coords| ThreeMomentum {
            px: F(coords[0]),
            py: F(coords[1]),
            pz: F(coords[2]),
        })
        .collect();
    integrand.evaluate_momentum_configuration(
        model,
        &MomentumSpaceEvaluationInput {
            loop_momenta,
            graph_id,
            orientation,
        },
        use_arb_prec,
    )
}

fn evaluate_x_point_precise(
    integrand: &mut gammalooprs::integrands::process::ProcessIntegrand,
    model: &gammalooprs::model::Model,
    point: &[f64],
    discrete_dim: &[usize],
    mut force_radius: bool,
    use_arb_prec: bool,
) -> Result<gammalooprs::integrands::evaluation::PreciseEvaluationResult> {
    if integrand.get_n_dim() as isize == point.len() as isize - 1 {
        force_radius = true;
    }

    let point = point.iter().map(|&x| F(x)).collect::<Vec<F<f64>>>();
    let continuous_sample = if force_radius {
        point
            .get(1..)
            .ok_or_else(|| {
                eyre!("Force-radius x-space evaluation requires at least one radial coordinate.")
            })?
            .to_vec()
    } else {
        point
    };
    let sample = havana_sample(continuous_sample, discrete_dim);

    integrand.evaluate_sample_precise(&sample, model, F(1.0), use_arb_prec, Complex::new_zero())
}

fn evaluate_momentum_point_precise(
    integrand: &mut gammalooprs::integrands::process::ProcessIntegrand,
    model: &gammalooprs::model::Model,
    point: &[f64],
    graph_name: Option<&str>,
    orientation: Option<usize>,
    use_arb_prec: bool,
) -> Result<gammalooprs::integrands::evaluation::PreciseEvaluationResult> {
    if orientation.is_some() && graph_name.is_none() {
        return Err(eyre!(
            "A specific orientation requires a specific graph in momentum-space evaluation."
        ));
    }

    if point.len() % 3 != 0 {
        return Err(eyre!(
            "Momentum-space evaluation expects a multiple of 3 coordinates, got {}.",
            point.len()
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
        let orientation_count = integrand.graph_orientation_count(graph_id).ok_or_else(|| {
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

    let loop_momenta = point
        .chunks_exact(3)
        .map(|coords| ThreeMomentum {
            px: F(coords[0]),
            py: F(coords[1]),
            pz: F(coords[2]),
        })
        .collect();
    integrand.evaluate_momentum_configuration_precise(
        model,
        &MomentumSpaceEvaluationInput {
            loop_momenta,
            graph_id,
            orientation,
        },
        use_arb_prec,
    )
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

fn havana_sample(cont: Vec<F<f64>>, discrete_dimensions: &[usize]) -> Sample<F<f64>> {
    let mut sample = Sample::Continuous(F(1.0), cont);

    for &discrete_dimension in discrete_dimensions.iter().rev() {
        sample = Sample::Discrete(F(1.0), discrete_dimension, Some(Box::new(sample)));
    }

    sample
}
