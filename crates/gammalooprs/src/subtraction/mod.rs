use bincode_trait_derive::{Decode, Encode};
use eyre::Result;
use itertools::Itertools;
use linnet::half_edge::involution::EdgeVec;
use spenso::algebra::complex::Complex;
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
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
        self, F, FloatLike, GS,
        hyperdual_utils::{DualOrNot, new_constant, simple_n_deriv_shape},
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
    if settings.force_uv_dampers_to_one {
        return radius.one();
    }

    let normalizing_scale = match settings.dynamic_width {
        true => radius_star,
        false => e_cm,
    };

    let delta_r = radius - radius_star;
    let sliver_width = F::from_f64(settings.sliver_width) * normalizing_scale;

    if delta_r.abs() > sliver_width || (settings.smooth_sliver && delta_r.abs() == sliver_width) {
        return radius.zero();
    }

    let delta_r_sq = &delta_r * &delta_r;
    let width = F::from_f64(settings.gaussian_width) * normalizing_scale;
    let width_sq = &width * &width;

    let mut exponent = -&delta_r_sq / width_sq;
    if settings.smooth_sliver {
        // This even bump is one at the pole and flat to all orders at the edge.
        let gap = sliver_width.square() - &delta_r_sq;
        exponent -= delta_r_sq / gap;
    }
    exponent.exp()
}

fn evaluate_uv_damper_dual<T: FloatLike>(
    radius: &HyperDual<F<T>>,
    radius_star: &HyperDual<F<T>>,
    e_cm: &F<T>,
    settings: &UVLocalisationSettings,
) -> HyperDual<F<T>> {
    if settings.force_uv_dampers_to_one {
        return new_constant(radius, &radius.values[0].one());
    }

    let normalizing_scale = if settings.dynamic_width {
        radius_star.clone()
    } else {
        new_constant(radius, e_cm)
    };

    let delta_r = radius.clone() - radius_star.clone();
    let sliver_width =
        new_constant(radius, &F::from_f64(settings.sliver_width)) * normalizing_scale.clone();

    if delta_r.values[0].abs() > sliver_width.values[0]
        || (settings.smooth_sliver && delta_r.values[0].abs() == sliver_width.values[0])
    {
        return new_constant(radius, &radius.values[0].zero());
    }

    let delta_r_sq = delta_r.clone() * delta_r;
    let width = new_constant(radius, &F::from_f64(settings.gaussian_width)) * normalizing_scale;
    let width_sq = width.clone() * width;

    let mut exponent = -delta_r_sq.clone() / width_sq;
    if settings.smooth_sliver {
        let gap = sliver_width.clone() * sliver_width - delta_r_sq.clone();
        exponent -= delta_r_sq / gap;
    }
    exponent.exp()
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
            // The helper's inverse radial measure leaves dr, so this density must
            // integrate to one in r: dr / r_star = d(r / r_star).
            let h = utils::h(&(radius / radius_star), None, None, h_function_settings);
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
                &(radius.clone() / radius_star.clone()),
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

        HyperDual::from_values(simple_n_deriv_shape(dual_values.len() - 1), dual_values)
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

    let rstar = parse!("r_star(t)");
    let e_surface = function!(GS.eta, rstar.clone(), Atom::var(t));

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
                Symbol::DERIVATIVE,
                current_r_derivative_counter,
                current_t_derivative_counter,
                GS.eta,
                rstar.clone(),
                Atom::var(t)
            );

            eta_derivatives_at_this_order.push(eta_derivative);

            if current_t_derivative_counter == 0 {
                break;
            }
            current_r_derivative_counter += 1;
            current_t_derivative_counter -= 1;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        integrands::process::GenericEvaluatorFloat,
        processes::cross_section::CrossSectionGraph,
        settings::runtime::{HFunction, HFunctionSettings},
        utils::hyperdual_utils::new_from_values,
    };

    #[test]
    fn integrated_ct_profiles_preserve_the_radial_residue() {
        for (function, power) in [
            (HFunction::Exponential, None),
            (HFunction::PolyExponential, None),
            (HFunction::PolyExponential, Some(4)),
            (HFunction::PolyExponential, Some(16)),
        ] {
            for sigma in [0.5, 2.0] {
                let settings = IntegratedCounterTermSettings {
                    range: IntegratedCounterTermRange::Infinite {
                        h_function_settings: HFunctionSettings {
                            function: function.clone(),
                            sigma,
                            power,
                            ..Default::default()
                        },
                    },
                };
                for radius_star in [0.3, 3.0] {
                    // The generated helper cancels r^(3L-1) from the measure.
                    // Both normalized profiles have negligible tails after 8 sigma.
                    let step = 8.0 * sigma * radius_star / 2048.0;
                    let integral: f64 = (0..2048)
                        .map(|i| {
                            evaluate_integrated_ct_normalisation(
                                &F((i as f64 + 0.5) * step),
                                &F(radius_star),
                                &F(1.0),
                                &settings,
                            )
                            .0
                        })
                        .sum::<f64>()
                        * step;
                    assert!(
                        (integral - 1.0).abs() < 2.0e-12,
                        "{function:?}, power={power:?}, sigma={sigma}, r_star={radius_star}: {integral}"
                    );
                }
            }
        }
    }

    #[test]
    fn integrated_ct_profile_dual_matches_normalized_density_derivatives() {
        let sigma = 0.6_f64;
        let radius_star = 1.7_f64;
        let settings = IntegratedCounterTermSettings {
            range: IntegratedCounterTermRange::Infinite {
                h_function_settings: HFunctionSettings {
                    function: HFunction::Exponential,
                    sigma,
                    ..Default::default()
                },
            },
        };
        let shape = HyperDual::new(simple_n_deriv_shape(2));
        let dual_star = new_from_values(&shape, &[F(radius_star), F(1.0), F(0.0)]);
        for radius in [0.2_f64, 1.4, 4.0] {
            let dual_radius = new_from_values(&shape, &[F(radius), F(0.0), F(0.0)]);
            let actual = evaluate_integrated_ct_normalisation_dual(
                &dual_radius,
                &dual_star,
                &F(1.0),
                &settings,
            );
            // Analytic half-Gaussian density and its first two r_star derivatives.
            let x = radius / (sigma * radius_star);
            let density =
                2.0 * (-x * x).exp() / (std::f64::consts::PI.sqrt() * sigma * radius_star);
            let log_derivative = (2.0 * x * x - 1.0) / radius_star;
            let second_log_derivative = (1.0 - 6.0 * x * x) / radius_star.powi(2);
            let expected = [
                density,
                density * log_derivative,
                density * (log_derivative.powi(2) + second_log_derivative) / 2.0,
            ];
            for (actual, expected) in actual.values.iter().zip(expected) {
                assert!((actual.0 - expected).abs() < 2.0e-13 * expected.abs().max(1.0));
            }
        }
    }

    #[test]
    fn uv_smooth_sliver_preserves_pole_symmetry_and_hard_mode() {
        for dynamic_width in [false, true] {
            let settings = UVLocalisationSettings {
                smooth_sliver: true,
                sliver_width: 0.5,
                dynamic_width,
                ..Default::default()
            };
            let radius_star = F(2.0_f64);
            let e_cm = F(4.0_f64);
            let half_width = if dynamic_width { 1.0 } else { 2.0 };
            assert_eq!(
                evaluate_uv_damper(&radius_star, &radius_star, &e_cm, &settings),
                F(1.0)
            );
            for fraction in [0.25, 0.75] {
                let delta = F(fraction * half_width);
                assert_eq!(
                    evaluate_uv_damper(&(radius_star + delta), &radius_star, &e_cm, &settings),
                    evaluate_uv_damper(&(radius_star - delta), &radius_star, &e_cm, &settings)
                );
            }
            let hard = UVLocalisationSettings {
                smooth_sliver: false,
                ..settings.clone()
            };
            for sign in [-1.0, 1.0] {
                let edge = F(radius_star.0 + sign * half_width);
                assert_eq!(
                    evaluate_uv_damper(&edge, &radius_star, &e_cm, &settings),
                    F(0.0)
                );
                // The existing hard sliver includes its endpoints.
                assert_eq!(
                    evaluate_uv_damper(&edge, &radius_star, &e_cm, &hard),
                    F((-0.25_f64).exp())
                );
            }
            let forced = UVLocalisationSettings {
                force_uv_dampers_to_one: true,
                ..settings
            };
            assert_eq!(
                evaluate_uv_damper(&F(100.0), &radius_star, &e_cm, &forced),
                F(1.0)
            );
        }
    }

    #[test]
    fn uv_smooth_sliver_dual_matches_derivatives_and_flat_boundary() {
        let shape = HyperDual::new(simple_n_deriv_shape(2));
        let radius_star = new_from_values(&shape, &[F(2.0_f64), F(0.7), F(0.0)]);
        for dynamic_width in [false, true] {
            let settings = UVLocalisationSettings {
                smooth_sliver: true,
                sliver_width: 0.5,
                dynamic_width,
                ..Default::default()
            };
            let scale = if dynamic_width { 2.0 } else { 4.0 };
            let scale_derivative = if dynamic_width { 0.7 } else { 0.0 };
            for u in [-0.4_f64, 0.0, 0.125] {
                let radius = new_from_values(&shape, &[F(2.0 + u * scale), F(0.3), F(0.0)]);
                let actual = evaluate_uv_damper_dual(&radius, &radius_star, &F(4.0), &settings);
                // Differentiate the dimensionless even profile, including S(t).
                let gap = 0.25 - u * u;
                let value = (-u * u - u * u / gap).exp();
                let log_first = -2.0 * u - 0.5 * u / gap.powi(2);
                let log_second = -2.0 - 0.5 / gap.powi(2) - 2.0 * u * u / gap.powi(3);
                let u_first = (-0.4 - u * scale_derivative) / scale;
                let u_second = -2.0 * scale_derivative * u_first / scale;
                let expected = [
                    value,
                    value * log_first * u_first,
                    value
                        * ((log_first.powi(2) + log_second) * u_first.powi(2)
                            + log_first * u_second)
                        / 2.0,
                ];
                let scalar = evaluate_uv_damper(
                    &radius.values[0],
                    &radius_star.values[0],
                    &F(4.0),
                    &settings,
                );
                assert!((actual.values[0].0 - scalar.0).abs() < 2.0e-14);
                for (actual, expected) in actual.values.iter().zip(expected) {
                    assert!((actual.0 - expected).abs() < 2.0e-13 * expected.abs().max(1.0));
                }
            }
            for u in [-0.6_f64, -0.5, -0.4995, 0.4995, 0.5, 0.6] {
                let radius = new_from_values(&shape, &[F(2.0 + u * scale), F(0.3), F(0.0)]);
                let actual = evaluate_uv_damper_dual(&radius, &radius_star, &F(4.0), &settings);
                for coefficient in actual.values {
                    assert!(coefficient.0.abs() < 1.0e-190);
                    if u.abs() >= 0.5 {
                        assert_eq!(coefficient, F(0.0));
                    }
                }
            }
            let forced = UVLocalisationSettings {
                force_uv_dampers_to_one: true,
                ..settings
            };
            let radius = new_from_values(&shape, &[F(100.0), F(0.3), F(0.0)]);
            assert_eq!(
                evaluate_uv_damper_dual(&radius, &radius_star, &F(4.0), &forced).values,
                vec![F(1.0), F(0.0), F(0.0)]
            );
        }
    }

    #[test]
    fn uv_smooth_sliver_helper_preserves_signed_radial_principal_value() {
        crate::initialisation::test_initialise().unwrap();
        let settings = UVLocalisationSettings {
            smooth_sliver: true,
            sliver_width: 1.5,
            ..Default::default()
        };
        for loop_count in [1, 2] {
            let (pieces, _) =
                CrossSectionGraph::single_th_prefactor_helper_atoms(1, loop_count, false, false);
            let mut helper = GenericEvaluator::new_from_raw_params(
                [pieces.local],
                &CrossSectionGraph::single_th_prefactor_helper_params(1, false),
                &FunctionMap::new(),
                vec![],
                OptimizationSettings::default(),
                None,
                &EvaluatorSettings::default(),
            )
            .unwrap()
            .into_eager_only();
            for radius_star in [0.5, 2.0] {
                // Equal midpoint bins pair around the pole. W > r_star in the
                // first case requires the helper's negative-radius mirror term.
                let step = 1.0 / 2048.0;
                let steps = ((radius_star + settings.sliver_width) / step) as usize;
                let mut integral = 0.0;
                for i in 0..steps {
                    let radius = F((i as f64 + 0.5) * step);
                    let plus = evaluate_uv_damper(&radius, &F(radius_star), &F(1.0), &settings);
                    let minus = evaluate_uv_damper(&(-radius), &F(radius_star), &F(1.0), &settings);
                    let params = [F(1.0), F(1.0), radius, F(radius_star), plus, minus, F(0.0)]
                        .map(Complex::new_re);
                    let value = f64::get_evaluator_single(&mut helper)(&params);
                    integral += value.re.0 * radius.0.powi(3 * loop_count as i32 - 1) * step;
                }
                assert!(
                    integral.abs() < 1.0e-10,
                    "{loop_count} loops, r_star={radius_star}: PV={integral}"
                );
            }
        }
    }

    #[test]
    fn uv_damper_can_be_forced_to_one() {
        let settings = UVLocalisationSettings {
            force_uv_dampers_to_one: true,
            ..Default::default()
        };

        assert_eq!(
            evaluate_uv_damper(&F(100.0_f64), &F(1.0_f64), &F(1.0_f64), &settings),
            F(1.0_f64)
        );
    }

    #[test]
    fn uv_damper_dual_can_be_forced_to_constant_one() {
        let settings = UVLocalisationSettings {
            force_uv_dampers_to_one: true,
            ..Default::default()
        };
        let shape = HyperDual::new(simple_n_deriv_shape(1));
        let radius = new_from_values(&shape, &[F(100.0_f64), F(5.0_f64)]);
        let radius_star = new_from_values(&shape, &[F(1.0_f64), F(3.0_f64)]);

        let damper = evaluate_uv_damper_dual(&radius, &radius_star, &F(1.0_f64), &settings);

        assert_eq!(damper.values, vec![F(1.0_f64), F(0.0_f64)]);
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
}
