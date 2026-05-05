use bincode_trait_derive::{Decode, Encode};
use eyre::Result;
use itertools::Itertools;
use linnet::half_edge::involution::EdgeVec;
use spenso::algebra::complex::Complex;
use symbolica::{
    atom::{Atom, AtomCore},
    domains::{dual::HyperDual, float::Real},
    evaluate::{FunctionMap, OptimizationSettings},
    function, parse, symbol,
};
use tracing::debug;
use typed_index_collections::TiVec;

use crate::{
    GammaLoopContext,
    cff::esurface::Esurface,
    graph::{LmbIndex, LoopMomentumBasis},
    integrands::process::GenericEvaluator,
    momentum::{
        Energy, FourMomentum,
        sample::{MomentumSample, SubspaceData},
    },
    processes::EvaluatorSettings,
    settings::runtime::{
        IntegratedCounterTermRange, IntegratedCounterTermSettings, UVLocalisationSettings,
    },
    utils::{
        self, F, FloatLike,
        hyperdual_utils::{DualOrNot, new_constant},
    },
};

pub mod amplitude_counterterm;
pub mod lu_counterterm;
pub mod overlap;
pub mod overlap_subspace;

fn evaluate_uv_damper<T: FloatLike>(
    radius: &F<T>,
    radius_star: &F<T>,
    e_cm: &F<T>,
    settings: &UVLocalisationSettings,
) -> F<T> {
    let normalizing_scale = match settings.dynamic_width {
        true => radius_star,
        false => e_cm,
    };

    let delta_r = radius - radius_star;

    if delta_r.abs() > F::from_f64(settings.sliver_width) * normalizing_scale {
        return radius.zero();
    }

    let delta_r_sq = &delta_r * &delta_r;
    let width = F::from_f64(settings.gaussian_width) * normalizing_scale;
    let width_sq = &width * &width;

    (-delta_r_sq / width_sq).exp()
}

fn evaluate_uv_damper_dual<T: FloatLike>(
    radius: &HyperDual<F<T>>,
    radius_star: &HyperDual<F<T>>,
    e_cm: &F<T>,
    settings: &UVLocalisationSettings,
) -> HyperDual<F<T>> {
    let normalizing_scale = if settings.dynamic_width {
        radius_star.clone()
    } else {
        new_constant(radius, e_cm)
    };

    let delta_r = radius.clone() - radius_star.clone();

    if delta_r.values[0].abs() > F::from_f64(settings.sliver_width) * &normalizing_scale.values[0] {
        return new_constant(radius, &radius.values[0].zero());
    }

    let delta_r_sq = delta_r.clone() * delta_r;
    let width = new_constant(radius, &F::from_f64(settings.gaussian_width)) * normalizing_scale;
    let width_sq = width.clone() * width;

    (-delta_r_sq / width_sq).exp()
}

fn evaluate_integrated_ct_normalisation<T: FloatLike>(
    radius: &F<T>,
    radius_star: &F<T>,
    _e_cm: &F<T>,
    settings: &IntegratedCounterTermSettings,
) -> F<T> {
    match &settings.range {
        IntegratedCounterTermRange::Infinite {
            h_function_settings,
        } => {
            let h = utils::h(&(radius_star / radius), None, None, h_function_settings);
            h * (radius_star).inv()
        }
        IntegratedCounterTermRange::Compact {} => {
            todo!();
        }
    }
}

fn evaluate_integrated_ct_normalisation_dual<T: FloatLike>(
    radius: &HyperDual<F<T>>,
    radius_star: &HyperDual<F<T>>,
    _e_cm: &F<T>,
    settings: &IntegratedCounterTermSettings,
) -> HyperDual<F<T>> {
    match &settings.range {
        IntegratedCounterTermRange::Infinite {
            h_function_settings,
        } => {
            let h = utils::h_dual(
                &(radius_star.clone() / radius.clone()),
                None,
                None,
                h_function_settings,
            );
            h / radius_star.clone()
        }
        IntegratedCounterTermRange::Compact {} => {
            todo!();
        }
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct RstarTDependenceEvaluator {
    dual_shape_for_esurface_evaluation: Vec<Vec<usize>>,
    implicit_function_theorem: Option<GenericEvaluator>,
}

impl RstarTDependenceEvaluator {
    pub(crate) fn supports_t_derivatives(&self) -> bool {
        self.implicit_function_theorem.is_some()
    }

    fn evaluate<T: FloatLike>(&mut self, input: RstarTDependenceInput<'_, T>) -> HyperDual<F<T>> {
        let RstarTDependenceInput {
            t_star,
            radius_star,
            overlap_center,
            subspace,
            unrescaled_momentum_sample,
            masses,
            threshold_esurface,
            lmb,
            all_lmbs,
        } = input;
        debug!("t-star: {}", t_star);
        debug!("r-star: {}", radius_star);

        let dual = HyperDual::new(self.dual_shape_for_esurface_evaluation.clone());
        let dual_rstar = dual.variable(1, radius_star.clone());
        let dual_t = dual.variable(0, t_star.clone());

        let rescale_tstar = unrescaled_momentum_sample
            .loop_moms()
            .rescale_with_hyper_dual(&dual_t, None);

        let dualized_externals_three_momenta = unrescaled_momentum_sample
            .external_moms()
            .iter()
            .map(|fm: &FourMomentum<F<T>>| fm.spatial.map_ref(&|x| new_constant(&dual, x)))
            .collect();

        let lmb_transform = rescale_tstar.lmb_transform(
            lmb,
            subspace.get_lmb(all_lmbs),
            &dualized_externals_three_momenta,
        );

        let dualized_overlap_center = overlap_center
            .iter()
            .map(|momentum| momentum.map_ref(&|x| new_constant(&dual, x)))
            .collect::<crate::momentum::sample::LoopMomenta<_>>();

        let shifted_t_dependent_momenta = lmb_transform
            .iter()
            .zip(dualized_overlap_center.iter())
            .map(|(momentum, center)| momentum.clone() - center.clone())
            .collect::<crate::momentum::sample::LoopMomenta<_>>();

        let zero: HyperDual<F<T>> = new_constant(&dual, &radius_star.zero());

        let shifted_radius_squared: HyperDual<F<T>> = match subspace.as_subspace_simple() {
            None => shifted_t_dependent_momenta
                .iter()
                .map(|momentum| momentum.norm_squared())
                .fold(zero.clone(), |acc, norm_squared| acc + norm_squared),
            Some(indices) => shifted_t_dependent_momenta
                .iter_enumerated()
                .filter(|(loop_index, _)| indices.contains(loop_index))
                .map(|(_, momentum)| momentum.norm_squared())
                .fold(zero, |acc, norm_squared| acc + norm_squared),
        };

        let shifted_radius = shifted_radius_squared.sqrt();
        let inverse_shifted_radius =
            new_constant(&shifted_radius, &radius_star.one()) / shifted_radius.clone();

        let unit_shifted_t_dependent_momenta = shifted_t_dependent_momenta
            .rescale(&inverse_shifted_radius, subspace.as_subspace_simple());

        let dualized_external_fourmomenta = dualized_externals_three_momenta
            .into_iter()
            .zip(unrescaled_momentum_sample.external_moms())
            .map(|(spatial, fm)| FourMomentum {
                spatial,
                temporal: Energy {
                    value: new_constant(&dual, &fm.temporal.value),
                },
            })
            .collect();

        let rescale_rstar = unit_shifted_t_dependent_momenta
            .rescale(&dual_rstar, subspace.as_subspace_simple())
            .iter()
            .zip(dualized_overlap_center.iter())
            .map(|(momentum, center)| momentum.clone() + center.clone())
            .collect::<crate::momentum::sample::LoopMomenta<_>>();

        let dual_esurface = threshold_esurface.compute_from_dual_momenta(
            subspace.get_lmb(all_lmbs),
            masses,
            &rescale_rstar,
            &dualized_external_fourmomenta,
        );

        debug!("Dual e-surface: {}", dual_esurface);

        let params = dual_esurface.values[1..]
            .iter()
            .map(|x| Complex::new_re(x.clone()))
            .collect_vec();

        debug!("Parameters for implicit function theorem: {:#?}", params);

        let result = T::get_evaluator(
            self.implicit_function_theorem
                .as_mut()
                .expect("r_star(t) evaluator requested without t-derivative support"),
        )(&params)
        .into_iter()
        .map(DualOrNot::unwrap_real)
        .collect_vec();

        debug!("Result from implicit function theorem: {:#?}", result);

        let mut dual_values = vec![radius_star.clone()];
        let mut n_factorial = 1;

        for (i, result) in result.into_iter().enumerate() {
            if i > 0 {
                n_factorial *= i;
                dual_values.push(result.re / F::from_f64(n_factorial as f64));
            } else {
                dual_values.push(result.re);
            }
        }

        todo!("construct hyperdual with t derivatives and r_star componenents")
        //       HyperDual::from_values(shape_for_t_derivatives(dual_values.len() - 1), dual_values)
    }
}

pub(crate) struct RstarTDependenceInput<'a, T: FloatLike> {
    pub t_star: &'a F<T>,
    pub radius_star: &'a F<T>,
    pub overlap_center: &'a crate::momentum::sample::LoopMomenta<F<T>>,
    pub subspace: &'a SubspaceData,
    pub unrescaled_momentum_sample: &'a MomentumSample<T>,
    pub masses: &'a EdgeVec<F<T>>,
    pub threshold_esurface: &'a Esurface,
    pub lmb: &'a LoopMomentumBasis,
    pub all_lmbs: &'a TiVec<LmbIndex, LoopMomentumBasis>,
}

// use the chain rule to express the t-derivatives of r_star in terms of the t and r derivatives of η(r_star(t), t)
pub(crate) fn generate_rstar_t_dependence_evaluator(
    num_t_derivatives: usize,
) -> Result<RstarTDependenceEvaluator> {
    if num_t_derivatives == 0 {
        return Ok(RstarTDependenceEvaluator {
            dual_shape_for_esurface_evaluation: Vec::new(),
            implicit_function_theorem: None,
        });
    }

    let t = symbol!("t");

    let e_surface = parse!("η(r_star(t), t)");
    let rstar = parse!("r_star(t)");

    let mut rstar_derivatives = vec![];
    for i in 0..num_t_derivatives {
        if i == 0 {
            rstar_derivatives.push(rstar.derivative(t));
        } else {
            rstar_derivatives.push(rstar_derivatives.last().unwrap().derivative(t));
        }
    }

    let mut equations = vec![];
    for i in 0..num_t_derivatives {
        if i == 0 {
            equations.push(e_surface.derivative(t));
        } else {
            equations.push(equations.last().unwrap().derivative(t));
        }
    }

    let mut solutions = equations
        .iter()
        .zip(&rstar_derivatives)
        .map(|(eq, variable)| {
            Atom::solve_linear_system::<u8, _, _>(&[eq], &[variable])
                .unwrap()
                .pop()
                .unwrap()
        })
        .collect_vec();

    for i in 1..solutions.len() {
        for j in 0..i {
            solutions[i] = solutions[i]
                .replace(rstar_derivatives[j].clone())
                .with(solutions[j].clone());
        }
    }

    // dual shape is for e-surface derivatives, implict function theorem should NOT be dualized with this
    let mut dual_shape = vec![vec![0, 0]];
    let mut params = vec![];
    for i in 1..=num_t_derivatives {
        let mut eta_derivatives_at_this_order = vec![];

        let mut current_r_derivative_counter = 0;
        let mut current_t_derivative_counter = i;

        loop {
            dual_shape.push(vec![
                current_r_derivative_counter,
                current_t_derivative_counter,
            ]);
            let eta_derivative = function!(
                symbolica::atom::Symbol::DERIVATIVE,
                current_r_derivative_counter,
                current_t_derivative_counter,
                e_surface.clone()
            );

            eta_derivatives_at_this_order.push(eta_derivative);

            if current_t_derivative_counter == 0 {
                break;
            }
            current_r_derivative_counter += 1;
            current_t_derivative_counter -= 1;
        }

        for param in &eta_derivatives_at_this_order {
            println!("  {}", param);
        }

        params.extend(eta_derivatives_at_this_order);
    }

    let fn_map = FunctionMap::new();
    let fn_map_entries = vec![];
    let implict_function_theorem = GenericEvaluator::new_from_raw_params(
        solutions,
        &params,
        &fn_map,
        fn_map_entries,
        OptimizationSettings::default(),
        None,
        &EvaluatorSettings::default(),
    )?;

    Ok(RstarTDependenceEvaluator {
        dual_shape_for_esurface_evaluation: dual_shape,
        implicit_function_theorem: Some(implict_function_theorem),
    })
}

#[test]
fn test_rstar_t_dependence_evaluator() {
    generate_rstar_t_dependence_evaluator(3).unwrap();
}

#[test]
fn test_rstar_t_dependence_evaluator_zero_derivatives() {
    let evaluator = generate_rstar_t_dependence_evaluator(0).unwrap();
    assert!(!evaluator.supports_t_derivatives());
}
