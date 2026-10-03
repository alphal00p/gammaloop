use color_eyre::Result;
use linnet::half_edge::involution::EdgeIndex;
use spenso::algebra::complex::Complex;
use symbolica::{
    atom::AtomCore,
    domains::dual::{DualNumberStructure, HyperDual},
    prelude::{FunctionMap, Real},
};
use three_dimensional_reps::ThermalDistributionFactor;

use super::{FermiSurfaceProduct, routing::test_graph};
use crate::{
    integrands::{
        evaluation::EvaluationMetaData,
        process::{GenericEvaluator, evaluators::evaluate_evaluator_single},
    },
    momentum::{
        ThreeMomentum,
        sample::{LoopIndex, LoopMomenta},
    },
    processes::cross_section::{
        build_derivative_structure, build_derivative_structure_atom, params_for_derivative_order,
    },
    settings::GlobalSettings,
    utils::{
        F,
        hyperdual_utils::{
            extract_t_derivatives, extract_t_derivatives_complex, new_constant,
            simple_n_deriv_shape,
        },
    },
};

#[test]
fn fermi_localization_matches_raised_cut_residues_on_nonlinear_shells() -> Result<()> {
    let graph = test_graph()?;
    let settings = GlobalSettings::default().generation.evaluator;
    let mut metadata = EvaluationMetaData::new_empty();
    let coefficient = |momenta: &LoopMomenta<HyperDual<F<f64>>>| {
        let shell = &momenta[LoopIndex(0)];
        let spectator = &momenta[LoopIndex(1)];
        let real = shell.px.clone() * &shell.px + shell.py.clone() * &shell.pz + &spectator.px;
        let imaginary = shell.px.clone() * &shell.pz + shell.py.clone() * &spectator.pz;
        HyperDual::from_values(
            shell
                .px
                .get_shape()
                .into_iter()
                .map(<[usize]>::to_vec)
                .collect(),
            real.values
                .into_iter()
                .zip(imaginary.values)
                .map(|(real, imaginary)| Complex::new(real, imaginary))
                .collect(),
        )
    };

    for order in 1..=5 {
        let product = FermiSurfaceProduct::new(
            &graph,
            &[ThermalDistributionFactor {
                edge_id: EdgeIndex(1),
                sign: 1,
                derivative_order: order,
            }],
        )?;
        let mut residue = if order <= 3 {
            build_derivative_structure(order as u8, -1, &settings)
        } else {
            // The shared symbolic builder leaves cancelling negative delta_t
            // powers unexpanded at these orders. Normalize the exact expression
            // for this reference; dropping or substituting those terms would
            // conceal a genuine residue error.
            let expression = build_derivative_structure_atom(order as u8, -1).expand();
            GenericEvaluator::new_from_raw_params(
                [expression],
                &params_for_derivative_order(order as u8),
                &FunctionMap::default(),
                vec![],
                settings.optimization_settings(),
                None,
                &settings,
            )?
            .into_eager_only()
        };
        for (mass, nu, shell) in [
            (1.25_f64, 2.5_f64, [0.4_f64, -0.8_f64, 1.3_f64]),
            (0.6_f64, 1.7_f64, [-1.1_f64, 0.6_f64, 0.5_f64]),
        ] {
            let mass = F(mass);
            let nu = F(nu);
            let momenta = LoopMomenta::from_iter([
                ThreeMomentum::new(F(shell[0]), F(shell[1]), F(shell[2])),
                ThreeMomentum::new(F(0.3_f64), F(-0.7_f64), F(0.9_f64)),
            ]);
            let localized = product.localize(
                &momenta,
                std::slice::from_ref(&mass),
                std::slice::from_ref(&nu),
                |_, scale| Ok((-scale.clone()).exp()),
                |momenta| Ok(coefficient(momenta)),
            )?;

            // The established cut path differentiates in t, including the
            // nonlinear g(t), whereas Fermi localization differentiates in E.
            let norm_squared = momenta[LoopIndex(0)].norm_squared();
            let root = ((nu * nu - mass * mass) / norm_squared).sqrt();
            let mut scale_values = vec![F(0.0_f64); order];
            scale_values[0] = root;
            if order > 1 {
                scale_values[1] = F(1.0_f64);
            }
            let scale = HyperDual::from_values(simple_n_deriv_shape(order - 1), scale_values);
            let mut rescaled = LoopMomenta::from_iter(
                momenta
                    .iter()
                    .map(|momentum| momentum.map_ref(&|value| new_constant(&scale, value))),
            );
            rescaled[LoopIndex(0)] =
                momenta[LoopIndex(0)].map_ref(&|value| new_constant(&scale, value) * &scale);
            let profile_and_volume = (-scale.clone()).exp() * &scale * &scale * &scale;
            let profile_and_volume = HyperDual::from_values(
                simple_n_deriv_shape(order - 1),
                profile_and_volume
                    .values
                    .into_iter()
                    .map(Complex::new_re)
                    .collect(),
            );
            let weighted = coefficient(&rescaled) * profile_and_volume;
            let mut params = extract_t_derivatives_complex(weighted);

            let surface_scale = HyperDual::new(simple_n_deriv_shape(order)).variable(0, root);
            let surface =
                (new_constant(&surface_scale, &norm_squared) * &surface_scale * &surface_scale
                    + new_constant(&surface_scale, &(mass * mass)))
                .sqrt()
                    - new_constant(&surface_scale, &nu);
            params.extend(
                extract_t_derivatives(surface)
                    .into_iter()
                    .skip(1)
                    .map(Complex::new_re),
            );
            let raw_residue = evaluate_evaluator_single(&mut residue, &params, &mut metadata);
            // Res A/g^n = (1/(n-1)!) d_E^(n-1)(A/g'). Here g' > 0,
            // so delta^(n-1)(g) adds precisely (-1)^(n-1) (n-1)!.
            let normalization = (1..order).fold(1.0_f64, |value, factor| -value * factor as f64);
            let expected = raw_residue * Complex::new_re(F(normalization));
            for (actual, expected) in [(localized.re, expected.re), (localized.im, expected.im)] {
                let error = (actual.0 - expected.0).abs();
                assert!(
                    error < 1.0e-10 * expected.0.abs().max(1.0),
                    "distribution order {order}, mass {mass}, nu {nu}: {actual} != {expected}"
                );
            }
        }
    }
    Ok(())
}
