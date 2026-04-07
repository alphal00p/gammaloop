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
        hyperdual_utils::{self, DualOrNot, new_constant},
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

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct RstarTDependenceEvaluator {
    dual_shape: Vec<Vec<usize>>,
    implicit_function_theorem: GenericEvaluator,
}

impl RstarTDependenceEvaluator {
    fn evaluate<T: FloatLike>(
        &mut self,
        t_star: &F<T>,
        radius_star: &F<T>,
        subspace: &SubspaceData,
        unrescaled_momentum_sample: &MomentumSample<T>,
        masses: &EdgeVec<F<T>>,
        e_surface: &Esurface,
        lmb: &LoopMomentumBasis,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    ) -> Vec<Complex<F<T>>> {
        let dual = HyperDual::new(self.dual_shape.clone());
        let dual_rstar = dual.variable(0, radius_star.clone());
        let dual_t = dual.variable(1, t_star.clone());

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

        let rescale_rstar = lmb_transform.rescale(&dual_rstar, subspace.as_subspace_simple());

        let dual_esurface = e_surface.compute_from_dual_momenta(
            lmb,
            masses,
            &rescale_rstar,
            &dualized_external_fourmomenta,
        );

        let params = dual_esurface.values[1..]
            .iter()
            .map(|x| Complex::new_re(x.clone()))
            .collect_vec();

        let result = T::get_evaluator(&mut self.implicit_function_theorem)(&params)
            .into_iter()
            .map(DualOrNot::unwrap_real)
            .collect_vec();

        result
    }
}

// use the chain rule to express the t-derivatives of r_star in terms of the t and r derivatives of η(r_star(t), t)
pub(crate) fn generate_rstar_t_dependence_evaluator(
    num_t_derivatives: usize,
) -> Result<RstarTDependenceEvaluator> {
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

    for (solution, rstar_derivative) in solutions.iter().zip(&rstar_derivatives) {
        println!("{} = {}", rstar_derivative, solution);
    }

    // dual shape is for e-surface derivatives, implict function theorem should NOT be dualized with this
    let mut dual_shape = vec![];
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

        println!("eta derivatives for order {}", i);
        for param in &eta_derivatives_at_this_order {
            println!("  {}", param);
        }

        params.extend(eta_derivatives_at_this_order);
    }

    println!("Dual shape: {:#?}", dual_shape);

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
        dual_shape,
        implicit_function_theorem: implict_function_theorem,
    })
}

#[test]
fn test_rstar_t_dependence_evaluator() {
    generate_rstar_t_dependence_evaluator(3).unwrap();
}
