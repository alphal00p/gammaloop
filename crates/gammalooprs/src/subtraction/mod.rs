use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::Real,
    function, parse, symbol,
};

use crate::{
    settings::runtime::{
        IntegratedCounterTermRange, IntegratedCounterTermSettings, UVLocalisationSettings,
    },
    utils::{self, F, FloatLike},
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

#[allow(dead_code)]
fn generate_rstar_t_dependence(num_t_derivatives: usize) -> Vec<Atom> {
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

    let solutions = equations
        .iter()
        .zip(&rstar_derivatives)
        .map(|(eq, variable)| {
            Atom::solve_linear_system::<u8, _, _>(&[eq], &[variable])
                .unwrap()
                .pop()
                .unwrap()
        })
        .collect();

    for i in 1..=num_t_derivatives {
        let mut params = vec![];

        let mut current_r_derivative_counter = 0;
        let mut current_t_derivative_counter = i;

        loop {
            let eta_derivative = function!(
                symbolica::atom::Symbol::DERIVATIVE,
                current_r_derivative_counter,
                current_t_derivative_counter,
                e_surface.clone()
            );

            params.push(eta_derivative);

            if current_t_derivative_counter == 0 {
                break;
            }
            current_r_derivative_counter += 1;
            current_t_derivative_counter -= 1;
        }

        println!("params for derivative order {}", i);
        for param in params {
            println!("  {}", param);
        }
    }

    solutions
}
