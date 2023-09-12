use colored::Colorize;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use symbolica::numerical_integration::Sample;

use crate::integrands::HasIntegrand;
use crate::utils;
use crate::Integrand;
use crate::Settings;
use lorentz_vector::LorentzVector;
use num::Complex;
use num_traits::ToPrimitive;

pub fn inspect(
    settings: &Settings,
    integrand: &mut Integrand,
    mut pt: Vec<f64>,
    mut force_radius: bool,
    is_momentum_space: bool,
    use_f128: bool,
) -> Complex<f64> {
    if integrand.get_n_dim() == pt.len() - 1 {
        force_radius = true;
    }

    let xs_f128 = if is_momentum_space {
        let (xs, inv_jac) = utils::global_inv_parameterize::<f128::f128>(
            &pt.chunks_exact_mut(3)
                .map(|x| LorentzVector::from_args(0., x[0], x[1], x[2]).cast())
                .collect::<Vec<LorentzVector<f128::f128>>>(),
            (settings.kinematics.e_cm * settings.kinematics.e_cm).into(),
            settings,
            force_radius,
        );
        if settings.general.debug > 1 {
            info!(
                "f128 sampling jacobian for this point = {:+.32e}",
                f128::f128::ONE / inv_jac
            );
        };
        xs
    } else {
        pt.iter().map(|x| f128::f128::from(*x)).collect::<Vec<_>>()
    };
    let xs_f64 = xs_f128
        .iter()
        .map(|x| f128::f128::to_f64(x).unwrap())
        .collect::<Vec<_>>();

    let eval = integrand.evaluate_sample(
        &Sample::Continuous(
            1.,
            if force_radius {
                xs_f64.clone()[1..].iter().map(|x| *x).collect::<Vec<_>>()
            } else {
                xs_f64.clone()
            },
        ),
        1.,
        1,
        use_f128,
    );
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
