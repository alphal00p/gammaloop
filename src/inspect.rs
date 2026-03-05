use color_eyre::Result;
use colored::Colorize;
use momtrop::float::MomTropFloat;
use spenso::algebra::complex::Complex;
use symbolica::domains::float::FloatLike;
use symbolica::numerical_integration::Sample;
#[allow(unused_imports)]
use tracing::{debug, error, info, instrument, trace, warn};

use crate::integrands::HasIntegrand;
use crate::model::Model;
use crate::momentum::ThreeMomentum;
use crate::settings::RuntimeSettings;
use crate::utils;
use crate::utils::F;
use crate::utils::f128;

#[allow(clippy::too_many_arguments)]
#[instrument(skip_all)]
pub fn inspect<I: HasIntegrand>(
    settings: &RuntimeSettings,
    integrand: &mut I,
    model: &Model,
    mut pt: Vec<F<f64>>,
    discrete_dimensions: &[usize],
    mut force_radius: bool,
    is_momentum_space: bool,
    use_arb_prec: bool,
) -> Result<(Option<f64>, Complex<F<f64>>)> {
    if integrand.get_n_dim() as isize == pt.len() as isize - 1 {
        force_radius = true;
    }

    let (jac, xs_f128) = if is_momentum_space {
        let (xs, inv_jac) = utils::global_inv_parameterize::<f128>(
            &pt.chunks_exact_mut(3)
                .map(|x| ThreeMomentum::new(x[0], x[1], x[2]).higher())
                .collect::<Vec<ThreeMomentum<F<f128>>>>(),
            F(settings.kinematics.e_cm).higher(),
            // TODO! return an error instead of panic
            &settings
                .sampling
                .get_parameterization_settings()
                .unwrap_or_else(|| panic!("Invalid sampling method for momentum-space inspect.")),
            force_radius,
        );

        info!(
            "f128 sampling jacobian for this point = {:+.32e}",
            FloatLike::inv(&inv_jac)
        );

        (Some(FloatLike::inv(&inv_jac).to_f64()), xs)
    } else {
        (
            None,
            pt.iter()
                .map(|x| F::<f128>::from_ff64(*x))
                .collect::<Vec<_>>(),
        )
    };
    let xs_f64 = xs_f128.iter().map(|x| F(x.into_f64())).collect::<Vec<_>>();

    let sample = {
        let cont_sample = if force_radius {
            xs_f64.clone()[1..].to_vec()
        } else {
            xs_f64.clone()
        };
        havana_sample(cont_sample, discrete_dimensions)
    };

    let eval_result =
        integrand.evaluate_sample(&sample, model, F(0.), 1, use_arb_prec, Complex::new_zero())?;
    let eval = eval_result.integrand_result;

    info!(
        "\nInput point in unit hypercube xs: \n\n{}\n\nThe evaluation of integrand '{}' is:\n\n{}\n",
        format!(
            "( {} )",
            xs_f64
                .iter()
                .map(|&x| format!("{:.16}", x))
                .collect::<Vec<_>>()
                .join(", ")
        )
        .blue(),
        integrand.name().to_string().green(),
        format!("( {:+.16e}, {:+.16e} i)", eval.re, eval.im).blue(),
    );

    Ok((jac, eval))
}

pub(crate) fn havana_sample(cont: Vec<F<f64>>, discrete_dimensions: &[usize]) -> Sample<F<f64>> {
    let mut sample = Sample::Continuous(F(1.), cont);

    for &d in discrete_dimensions.iter().rev() {
        sample = Sample::Discrete(F(1.), d, Some(Box::new(sample)));
    }

    sample
}
