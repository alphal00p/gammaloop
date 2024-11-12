use colored::Colorize;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use spenso::complex::Complex;
use symbolica::domains::float::NumericalFloatLike;
use symbolica::numerical_integration::Sample;

use crate::debug_info::DEBUG_LOGGER;
use crate::integrands::HasIntegrand;
use crate::momentum::ThreeMomentum;
use crate::utils;
use crate::utils::f128;
use crate::utils::F;
use crate::Integrand;
use crate::SamplingSettings;
use crate::Settings;

pub fn inspect(
    settings: &Settings,
    integrand: &mut Integrand,
    mut pt: Vec<F<f64>>,
    _term: &[usize],
    mut force_radius: bool,
    is_momentum_space: bool,
    use_f128: bool,
) -> Complex<F<f64>> {
    DEBUG_LOGGER.reset();

    if integrand.get_n_dim() == pt.len() - 1 {
        force_radius = true;
    }

    let xs_f128 = if is_momentum_space {
        let (xs, inv_jac) = utils::global_inv_parameterize::<f128>(
            &pt.chunks_exact_mut(3)
                .map(|x| ThreeMomentum::new(x[0], x[1], x[2]).higher())
                .collect::<Vec<ThreeMomentum<F<f128>>>>(),
            settings.kinematics.e_cm.square().higher(),
            settings,
            force_radius,
        );
        if settings.general.debug > 1 {
            info!(
                "f128 sampling jacobian for this point = {:+.32e}",
                inv_jac.inv()
            );
        };
        xs
    } else {
        pt.iter()
            .map(|x| F::<f128>::from_ff64(*x))
            .collect::<Vec<_>>()
    };
    let xs_f64 = xs_f128.iter().map(|x| F(x.into_f64())).collect::<Vec<_>>();

    let sample = {
        let cont_sample = Sample::Continuous(
            F(1.),
            if force_radius {
                xs_f64.clone()[1..].to_vec()
            } else {
                xs_f64.clone()
            },
        );
        match &settings.sampling {
            SamplingSettings::DiscreteGraphs(_) => {
                let graph_id = _term[0];
                Sample::Discrete(F(1.), graph_id, Some(Box::new(cont_sample)))
            }
            _ => cont_sample,
        }
    };

    let eval_result = integrand.evaluate_sample(&sample, F(0.), 1, use_f128, F(0.0));
    let eval = eval_result.integrand_result;

    info!(
        "\nFor input point xs: \n\n{}\n\nThe evaluation of integrand '{}' is:\n\n{}\n",
        format!(
            "( {} )",
            xs_f64
                .iter()
                .map(|&x| format!("{:.16}", x))
                .collect::<Vec<_>>()
                .join(", ")
        )
        .blue(),
        format!("{}", settings.hard_coded_integrand).green(),
        format!("( {:+.16e}, {:+.16e} i)", eval.re, eval.im).blue(),
    );

    eval
}
