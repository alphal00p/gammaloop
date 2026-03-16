use symbolica::domains::float::Real;

use crate::{
    settings::runtime::{
        IntegratedCounterTermRange, IntegratedCounterTermSettings, UVLocalisationSettings,
    },
    utils::{self, F, FloatLike},
};

pub mod overlap;
//pub mod static_counterterm;
pub mod amplitude_counterterm;
pub mod lu_counterterm;
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
