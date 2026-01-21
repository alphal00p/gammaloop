use std::iter;

use crate::{
    momentum::{Dep, Energy, ExternalMomenta, FourMomentum, SignOrZero, ThreeMomentum},
    momentum_sample::ExternalIndex,
    settings::runtime::kinematic::Externals,
    signature::SignatureLike,
    utils::{F, FloatLike, newton_solver::newton_iteration_and_derivative},
};

use crate::utils::serde_utils::{IsDefault, is_false};
use bincode_trait_derive::{Decode, Encode};
use eyre::Result;
use itertools::Itertools;
use rand::Rng;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use symbolica::{
    domains::float::{FloatLike as SymFloatLike, Real},
    numerical_integration::MonteCarloRng,
};
use tracing::{debug, warn};
use typed_index_collections::TiVec;

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Deserialize, Serialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct PhaseSpaceImprovementSettings {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub mode: ImprovementMode,
    #[serde(skip_serializing_if = "_default_large_def_check")]
    pub large_deformation_check: Option<F<f64>>,
    #[serde(skip_serializing_if = "is_false")]
    pub only_warn_on_large_deformation: bool,
}

fn _default_large_def_check(val: &Option<F<f64>>) -> bool {
    *val == Some(F(1e-12))
}

impl Default for PhaseSpaceImprovementSettings {
    fn default() -> Self {
        Self {
            mode: ImprovementMode::None,
            large_deformation_check: Some(F(1e-12)),
            only_warn_on_large_deformation: false,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Deserialize, Serialize, Encode, Decode, PartialEq, JsonSchema, Default)]
#[serde(deny_unknown_fields)]
pub enum ImprovementMode {
    #[serde(rename = "vh")]
    Vh,
    #[serde(rename = "mf")]
    Mf,
    #[default]
    #[serde(rename = "none")]
    None,
}

fn dimensionless_metric<T: FloatLike>(
    externals_1: &TiVec<ExternalIndex, FourMomentum<F<T>>>,
    externals_2: &TiVec<ExternalIndex, FourMomentum<F<T>>>,
    e_cm: &F<T>,
) -> F<T> {
    externals_1
        .iter()
        .zip(externals_2)
        .fold(e_cm.zero(), |acc, (p1, p2)| acc + (p1 - p2).square())
        / e_cm.square()
        / F::from_f64(externals_1.len() as f64)
}

pub(crate) fn improve_ps<T: FloatLike>(
    dependent_momenta: &TiVec<ExternalIndex, FourMomentum<F<T>>>,
    external_masses: &TiVec<ExternalIndex, F<T>>,
    external_signature: &SignatureLike<ExternalIndex>,
    e_cm: &F<T>,
    settings: &PhaseSpaceImprovementSettings,
) -> Result<TiVec<ExternalIndex, FourMomentum<F<T>>>> {
    let result = match settings.mode {
        ImprovementMode::Vh => {
            improve_ps_vh(dependent_momenta, external_masses, external_signature, e_cm)
        }
        ImprovementMode::Mf => {
            improve_ps_mf(dependent_momenta, external_masses, external_signature, e_cm)
        }
        ImprovementMode::None => Ok(dependent_momenta.clone()),
    }?;

    let deformation_size = dimensionless_metric(dependent_momenta, &result, e_cm);

    if let Some(threshold) = settings.large_deformation_check
        && deformation_size > F::from_ff64(threshold)
    {
        if settings.only_warn_on_large_deformation {
            warn!(
                "Phase space improvement resulted in a large deformation: {:+e} > {:+e}",
                deformation_size, threshold
            );
        } else {
            return Err(eyre::eyre!(
                "Phase space improvement resulted in a large deformation: {:+e} > {:+e}",
                deformation_size,
                threshold
            ));
        }
    }

    Ok(result)
}

pub(crate) fn generate_default_momenta(
    external_masses: &TiVec<ExternalIndex, F<f64>>,
    external_signature: &SignatureLike<ExternalIndex>,
    e_cm: &F<f64>,
) -> Result<Externals> {
    let initial_states = external_signature
        .iter_enumerated()
        .filter_map(|(id, s)| {
            if *s == SignOrZero::Plus {
                Some(id)
            } else {
                None
            }
        })
        .collect_vec();

    let final_states = external_signature
        .iter_enumerated()
        .filter_map(|(id, s)| {
            if *s == SignOrZero::Minus {
                Some(id)
            } else {
                None
            }
        })
        .collect_vec();

    let num_initial = initial_states.len();
    let num_final = final_states.len();

    if num_initial == 0 || num_final == 0 {
        return Err(eyre::eyre!(
            "Only processes with at least one initial and one final state are supported"
        ));
    }

    let final_state_mass_sum = final_states
        .iter()
        .map(|i| &external_masses[*i])
        .fold(F(0.0), |a, b| a + b);

    let initial_state_mass_sum = initial_states
        .iter()
        .map(|i| &external_masses[*i])
        .fold(F(0.0), |a, b| a + b);

    // maybe add a fudge factor?
    let incoming_energy = e_cm.max(final_state_mass_sum).max(initial_state_mass_sum);

    let capital_p = FourMomentum {
        temporal: Energy {
            value: incoming_energy,
        },
        spatial: ThreeMomentum {
            px: F(0.0),
            py: F(0.0),
            pz: F(0.0),
        },
    };

    if num_initial == 1 && num_final == 1 {
        return Ok(Externals::Constant {
            momenta: vec![
                ExternalMomenta::Independent([incoming_energy, F(0.0), F(0.0), incoming_energy]),
                ExternalMomenta::Dependent(Dep::Dep),
            ],
            improvement_settings: PhaseSpaceImprovementSettings::default(),
            helicities: vec![SignOrZero::Plus, SignOrZero::Plus],
            f_64_cache: None,
            f_128_cache: None,
        });
    }

    // start a rng with a fixed seed
    let mut rng = MonteCarloRng::new(1234567, 0);

    let mut final_state_momenta = (0..num_final - 1)
        .map(|i| {
            let sign = if i % 2 == 0 { F(1.0) } else { F(-1.0) };
            let spatial = ThreeMomentum {
                px: F(rng.random::<f64>() + 1.0) * incoming_energy / F(num_final as f64 * 10.0),
                py: F(rng.random::<f64>() + 1.0) * incoming_energy / F(num_final as f64 * 10.0),
                pz: F(rng.random::<f64>() + 1.0) * incoming_energy / F(num_final as f64 * 10.0),
            };
            let signed_spatial = spatial * sign;
            let ose = spatial.on_shell_energy(Some(external_masses[final_states[i]]));
            FourMomentum {
                temporal: ose,
                spatial: signed_spatial,
            }
        })
        .collect_vec();

    let last_final_state = final_state_momenta
        .iter()
        .fold(capital_p, |acc, mom| acc - *mom);
    final_state_momenta.push(last_final_state);

    let mut initial_state_momenta = (0..num_initial - 1)
        .map(|i| {
            let (sign_x, sign_y, sign_z) = (0..3)
                .map(|_| {
                    if rng.random_bool(0.5) {
                        F(1.0)
                    } else {
                        F(-1.0)
                    }
                })
                .collect_tuple()
                .unwrap();

            let spatial = ThreeMomentum {
                px: F(rng.random::<f64>() + 1.0) * incoming_energy / F(num_final as f64 * 10.0)
                    * sign_x,
                py: F(rng.random::<f64>() + 1.0) * incoming_energy / F(num_final as f64 * 10.0)
                    * sign_y,
                pz: F(rng.random::<f64>() + 1.0) * incoming_energy / F(num_final as f64 * 10.0)
                    * sign_z,
            };
            let ose = spatial.on_shell_energy(Some(external_masses[final_states[i]]));
            FourMomentum {
                temporal: ose,
                spatial,
            }
        })
        .collect_vec();

    let last_initial_state = initial_state_momenta
        .iter()
        .fold(capital_p, |acc, mom| acc - *mom);
    initial_state_momenta.push(last_initial_state);

    let mut initial_state_iter = initial_state_momenta.into_iter();
    let mut final_state_iter = final_state_momenta.into_iter();

    let momenta_to_rescale = external_signature
        .iter()
        .map(|s| match s {
            SignOrZero::Plus => initial_state_iter.next().unwrap(),
            SignOrZero::Minus => final_state_iter.next().unwrap(),
            SignOrZero::Zero => unreachable!(),
        })
        .collect::<TiVec<ExternalIndex, FourMomentum<F<f64>>>>();

    let initial_states_rescaled = find_rescaling(
        &momenta_to_rescale,
        external_masses,
        external_signature,
        e_cm,
        &incoming_energy,
        SignOrZero::Plus,
    );

    let all_states_rescaled = find_rescaling(
        &initial_states_rescaled,
        external_masses,
        external_signature,
        e_cm,
        &incoming_energy,
        SignOrZero::Minus,
    );

    let mut new_externals_momenta = vec![ExternalMomenta::Dependent(Dep::Dep)];
    for mom in all_states_rescaled.iter().skip(1) {
        new_externals_momenta.push(ExternalMomenta::Independent([
            mom.temporal.value,
            mom.spatial.px,
            mom.spatial.py,
            mom.spatial.pz,
        ]));
    }

    let helicities = iter::from_fn(|| Some(SignOrZero::Plus))
        .take(num_final + num_initial)
        .collect();

    Ok(Externals::Constant {
        momenta: new_externals_momenta,
        improvement_settings: PhaseSpaceImprovementSettings::default(),
        helicities,
        f_64_cache: None,
        f_128_cache: None,
    })
}

fn find_rescaling<T: FloatLike>(
    momenta: &TiVec<ExternalIndex, FourMomentum<F<T>>>,
    external_masses: &TiVec<ExternalIndex, F<T>>,
    external_signature: &SignatureLike<ExternalIndex>,
    e_cm: &F<T>,
    total_energy: &F<T>,
    selected_sign: SignOrZero,
) -> TiVec<ExternalIndex, FourMomentum<F<T>>> {
    let function = |x: &F<T>| {
        let (energy_sum, der) = momenta
            .iter()
            .enumerate()
            .zip(external_signature)
            .zip(external_masses)
            .map(|(((_id, mom), sign), mass)| {
                if selected_sign == sign {
                    let rescaled_mom = mom.rescale_spatial(x);
                    let ose = rescaled_mom
                        .spatial
                        .on_shell_energy(Some(mass.clone()))
                        .value;
                    let der = x * mom.spatial.norm_squared() / &ose;
                    // println!(
                    //     "Rescaling x: {:?}, mom: {:?}, ose: {:?}, der: {:?}",
                    //     x, mom, ose, der
                    // );

                    (ose, der)
                } else {
                    (F::from_f64(0.0), F::from_f64(0.0))
                }
            })
            .reduce(|(ose_a, der_a), (ose_b, der_b)| ((ose_a + ose_b), (der_a + der_b)))
            .unwrap_or((F::from_f64(0.0), F::from_f64(0.0)));

        (energy_sum - total_energy, der)
    };

    let solution = newton_iteration_and_derivative(
        &F::from_f64(1.0),
        function,
        &F::<T>::from_f64(1.000),
        400,
        e_cm,
    );

    debug!("solution: {:?}", solution);

    momenta
        .iter()
        .enumerate()
        .zip(external_signature)
        .zip(external_masses)
        .map(|(((_id, mom), sign), mass)| {
            if selected_sign == sign {
                let mut rescaled_mom = mom.rescale_spatial(&solution.solution);
                let ose = rescaled_mom.spatial.on_shell_energy(Some(mass.clone()));
                rescaled_mom.temporal = ose;
                rescaled_mom
            } else {
                mom.clone()
            }
        })
        .collect::<TiVec<ExternalIndex, FourMomentum<F<T>>>>()
}

// Here we use Lorentz invariance to align sum of incoming momenta along the z-axis, and boost to the CM frame.
// We can then do a psmc rescaling in just the final states without the need to add an x^2 term
// In principle, it would be possible to do an inverse transformation to go back to the original frame
fn improve_ps_mf<T: FloatLike>(
    dependent_momenta: &TiVec<ExternalIndex, FourMomentum<F<T>>>,
    external_masses: &TiVec<ExternalIndex, F<T>>,
    external_signature: &SignatureLike<ExternalIndex>,
    e_cm: &F<T>,
) -> Result<TiVec<ExternalIndex, FourMomentum<F<T>>>> {
    let initial_states = external_signature
        .iter_enumerated()
        .filter_map(|(id, s)| {
            if *s == SignOrZero::Plus {
                Some(id)
            } else {
                None
            }
        })
        .collect_vec();

    let initial_state_sum = initial_states
        .iter()
        .map(|i| &dependent_momenta[*i])
        .fold(dependent_momenta[ExternalIndex::from(0)].zero(), |a, b| {
            a.clone() + b.clone()
        });

    let z_axis = ThreeMomentum {
        px: e_cm.zero(),
        py: e_cm.zero(),
        pz: e_cm.one(),
    };

    let initial_state_norm = initial_state_sum.spatial.norm();
    let initial_state_unit_vector = if initial_state_norm > e_cm.epsilon() * e_cm {
        &initial_state_sum.spatial * initial_state_sum.spatial.norm().inv()
    } else {
        z_axis.clone()
    };

    debug!("Initial state unit vector: {:?}", initial_state_unit_vector);

    let cos_theta_to_z_axis = initial_state_sum.spatial.get_cos_theta_with(&z_axis);

    debug!("cos theta to z axis: {:?}", cos_theta_to_z_axis);

    let rotated_dependent_momenta = dependent_momenta
        .iter()
        .map(|m| FourMomentum {
            temporal: m.temporal.clone(),
            spatial: m.spatial.axis_angle_rotation(&cos_theta_to_z_axis, &z_axis),
        })
        .collect::<TiVec<ExternalIndex, FourMomentum<F<T>>>>();

    debug!(
        "Rotated dependent momenta: {:#?}",
        rotated_dependent_momenta
    );

    let rotated_initial_state_sum = initial_states
        .iter()
        .map(|i| &rotated_dependent_momenta[*i])
        .fold(
            rotated_dependent_momenta[ExternalIndex::from(0)].zero(),
            |a, b| a.clone() + b.clone(),
        );

    let x = rotated_initial_state_sum.spatial.pz / &initial_state_sum.temporal.value;

    if x.abs() > e_cm.one() {
        return Err(eyre::eyre!(
            "Initial states with total spacelike momentum can not be boosted to rest frame"
        ));
    }
    let two = e_cm.one() + e_cm.one();
    let rapiditiy = ((e_cm.one() + &x) / (e_cm.one() - &x)).ln() / &two;
    let cosh_eta = rapiditiy.cosh();
    let sinh_eta = rapiditiy.sinh();

    let mut boosted_momenta = rotated_dependent_momenta
        .iter()
        .map(|m| FourMomentum {
            temporal: Energy {
                value: &m.temporal.value * &cosh_eta - &m.spatial.pz * &sinh_eta,
            },
            spatial: ThreeMomentum {
                px: m.spatial.px.clone(),
                py: m.spatial.py.clone(),
                pz: -&m.temporal.value * &sinh_eta + &m.spatial.pz * &cosh_eta,
            },
        })
        .collect::<TiVec<ExternalIndex, FourMomentum<F<T>>>>();

    debug!("Boosted momenta: {:#?}", boosted_momenta);

    let mut total_energy = e_cm.zero();
    for ((boosted_momenta, sign), mass) in boosted_momenta
        .iter_mut()
        .zip(external_signature)
        .zip(external_masses)
    {
        // ensure the initial states are on-shell
        if sign == SignOrZero::Plus {
            let ose = boosted_momenta.spatial.on_shell_energy(Some(mass.clone()));
            boosted_momenta.temporal = ose;
            total_energy += &boosted_momenta.temporal.value;
        }
    }

    let rotated_boosted_rescaled_momenta = find_rescaling(
        &boosted_momenta,
        external_masses,
        external_signature,
        e_cm,
        &total_energy,
        SignOrZero::Minus,
    );

    debug!(
        "Rotated, boosted and rescaled momenta: {:#?}",
        rotated_boosted_rescaled_momenta
    );

    // transform back to original frame
    let rescaled_momenta = rotated_boosted_rescaled_momenta
        .into_iter()
        .map(|m| FourMomentum {
            temporal: Energy {
                value: &m.temporal.value * &cosh_eta + &m.spatial.pz * &sinh_eta,
            },
            spatial: ThreeMomentum {
                px: m.spatial.px.clone(),
                py: m.spatial.py.clone(),
                pz: &m.temporal.value * &sinh_eta + &m.spatial.pz * &cosh_eta,
            },
        })
        .map(|m| FourMomentum {
            temporal: m.temporal,
            spatial: m
                .spatial
                .axis_angle_rotation(&cos_theta_to_z_axis, &initial_state_unit_vector),
        })
        .collect::<TiVec<ExternalIndex, FourMomentum<F<T>>>>();

    debug!("Rescaled momenta: {:#?}", rescaled_momenta);
    Ok(rescaled_momenta)
}

fn improve_ps_vh<T: FloatLike>(
    dependent_momenta: &TiVec<ExternalIndex, FourMomentum<F<T>>>,
    external_masses: &TiVec<ExternalIndex, F<T>>,
    external_signature: &SignatureLike<ExternalIndex>,
    e_cm: &F<T>,
) -> Result<TiVec<ExternalIndex, FourMomentum<F<T>>>> {
    let initial_states = external_signature
        .iter_enumerated()
        .filter_map(|(id, s)| {
            if *s == SignOrZero::Plus {
                Some(id)
            } else {
                None
            }
        })
        .collect_vec();

    let final_states = external_signature
        .iter_enumerated()
        .filter_map(|(id, s)| {
            if *s == SignOrZero::Minus {
                Some(id)
            } else {
                None
            }
        })
        .collect_vec();

    if initial_states.len() != 2 {
        return Err(eyre::eyre!(
            "Only initial states with two incoming particles are supported for vh improvement"
        ));
    }

    if initial_states.iter().any(|i| {
        dependent_momenta[*i].spatial.px > e_cm * e_cm.epsilon()
            || dependent_momenta[*i].spatial.py > e_cm * e_cm.epsilon()
    }) {
        return Err(eyre::eyre!(
            "Initial states with transverse momentum are not supported for vh improvement"
        ));
    }

    if initial_states
        .iter()
        .any(|i| external_masses[*i] != F::from_f64(0.0))
    {
        return Err(eyre::eyre!(
            "Initial states with mass are not supported for vh improvement"
        ));
    }

    if final_states
        .iter()
        .any(|i| dependent_momenta[*i].temporal.value < F::from_f64(0.0))
    {
        return Err(eyre::eyre!(
            "Final states with negative energy are not supported for vh improvement"
        ));
    }

    let mut new_momenta = dependent_momenta
        .iter()
        .cloned()
        .collect::<TiVec<ExternalIndex, _>>();

    let mut pt = dependent_momenta[ExternalIndex::from(0)].zero();
    let two = e_cm.one() + e_cm.one();

    for final_state in final_states.iter() {
        let mom = &dependent_momenta[*final_state];
        let sign = if mom.spatial.pz > F::from_f64(0.0) {
            e_cm.one()
        } else {
            -e_cm.one()
        };
        let mass = &external_masses[*final_state];
        let new_z_value = (mom.temporal.value.square()
            - mom.spatial.px.square()
            - mom.spatial.py.square()
            - mass.square())
        .abs()
        .sqrt()
            * sign;

        new_momenta[*final_state].spatial.pz = new_z_value;
        pt += &new_momenta[*final_state];
    }

    let (p1, p2) = if dependent_momenta[initial_states[0]].spatial.pz > F::from_f64(0.0) {
        (initial_states[0], initial_states[1])
    } else if dependent_momenta[initial_states[1]].spatial.pz > F::from_f64(0.0) {
        (initial_states[1], initial_states[0])
    } else {
        return Err(eyre::eyre!(
            "At least one initial state should have positive pz for vh improvement"
        ));
    };

    let p1_vec = &dependent_momenta[p1];
    let p2_vec = &dependent_momenta[p2];

    let discr = pt.square().abs();

    let shift_e_1 = (&pt.temporal.value
        * (-&two * &p1_vec.temporal.value * &pt.temporal.value
            + pt.temporal.value.square()
            + pt.spatial.px.square()
            + pt.spatial.py.square())
        + (&two * &p1_vec.temporal.value - &pt.temporal.value) * pt.spatial.pz.square()
        + &pt.spatial.pz * &discr)
        / (&two * (&pt.temporal.value - &pt.spatial.pz) * (&pt.temporal.value + &pt.spatial.pz));

    let shift_e_2 = -(&pt.temporal.value
        * (&two * &p2_vec.temporal.value * &pt.temporal.value - pt.temporal.value.square()
            + pt.spatial.px.square()
            + pt.spatial.py.square())
        + (-&two * &p2_vec.temporal.value + &pt.temporal.value) * pt.spatial.pz.square()
        + &pt.spatial.pz * &discr)
        / (&two * (&pt.temporal.value - &pt.spatial.pz) * (&pt.temporal.value + &pt.spatial.pz));

    let shift_z_1 =
        (-&two * &p1_vec.spatial.pz * (pt.temporal.value.square() - pt.spatial.pz.square())
            + &pt.spatial.pz
                * (pt.temporal.value.square() + pt.spatial.px.square() + pt.spatial.py.square()
                    - pt.spatial.pz.square())
            + &pt.temporal.value * &discr)
            / (&two * (pt.temporal.value.square() - pt.spatial.pz.square()));

    let shift_z_2 =
        -(&two * &p2_vec.spatial.pz * (pt.temporal.value.square() - pt.spatial.pz.square())
            + &pt.spatial.pz
                * (-pt.temporal.value.square()
                    + pt.spatial.px.square()
                    + pt.spatial.py.square()
                    + pt.spatial.pz.square())
            + &pt.temporal.value * &discr)
            / (&two * (pt.temporal.value.square() - pt.spatial.pz.square()));

    new_momenta[p1].temporal.value = &p1_vec.temporal.value + &shift_e_1;
    new_momenta[p1].spatial.pz = &p1_vec.spatial.pz + &shift_z_1;
    new_momenta[p2].temporal.value = &p2_vec.temporal.value + &shift_e_2;
    new_momenta[p2].spatial.pz = &p2_vec.spatial.pz + &shift_z_2;
    new_momenta[p2].spatial.px = p2_vec.spatial.px.clone();
    new_momenta[p2].spatial.py = p2_vec.spatial.py.clone();

    let final_state_x_component = final_states
        .iter()
        .map(|s| &dependent_momenta[*s].spatial.px)
        .fold(e_cm.zero(), |a, b| a + b);

    let final_state_y_component = final_states
        .iter()
        .map(|s| &dependent_momenta[*s].spatial.py)
        .fold(e_cm.zero(), |a, b| a + b);

    let ref_x = final_state_x_component - &p2_vec.spatial.px;
    let ref_y = final_state_y_component - &p2_vec.spatial.py;

    new_momenta[p1].spatial.px = ref_x;
    new_momenta[p1].spatial.py = ref_y;

    Ok(new_momenta)
}

#[cfg(test)]
mod tests {
    use eyre::Result;
    use typed_index_collections::TiVec;

    use crate::{
        DependentMomentaConstructor, dot,
        graph::{FeynmanGraph, Graph, parse::IntoGraph},
        initialisation::test_initialise,
        model::Model,
        momentum::{FourMomentum, SignOrZero},
        momentum_sample::ExternalIndex,
        signature::SignatureLike,
        utils::{F, FloatLike, f128, test_utils::load_generic_model},
    };

    fn test_default_momenta_graph(graph: &Graph, model: &Model, e_cm: &F<f64>) -> Result<()> {
        let external_signature = graph.get_external_signature();
        let external_masses = graph.get_external_masses::<f64>(model);

        let default_momenta =
            super::generate_default_momenta(&external_masses, &external_signature, e_cm).unwrap();

        let constructor = DependentMomentaConstructor::Amplitude(&external_signature);
        let dependent_momenta = default_momenta
            .get_dependent_externals::<f64>(constructor)
            .unwrap();

        tracing::debug!("Default momenta: {:#?}", dependent_momenta);

        test_kinematic_validity(
            &dependent_momenta,
            &external_signature,
            &external_masses,
            e_cm,
        )?;

        let dependent_momenta_f128 = default_momenta
            .get_dependent_externals::<f128>(constructor)
            .unwrap();

        tracing::debug!("upcasted momenta f128: {:#?}", dependent_momenta_f128);
        let masses_f128 = graph.get_external_masses::<f128>(model);
        let improve_ps_settings = super::PhaseSpaceImprovementSettings::default();

        let improved_point = super::improve_ps(
            &dependent_momenta_f128,
            &masses_f128,
            &external_signature,
            &F::from_ff64(e_cm.clone()),
            &improve_ps_settings,
        )?;
        tracing::debug!("Improved momenta f128: {:#?}", improved_point);

        test_kinematic_validity(
            &improved_point,
            &external_signature,
            &masses_f128,
            &F::from_ff64(e_cm.clone()),
        )?;

        Ok(())
    }

    fn test_kinematic_validity<T: FloatLike>(
        externals: &TiVec<ExternalIndex, FourMomentum<F<T>>>,
        external_signature: &SignatureLike<ExternalIndex>,
        masses: &TiVec<ExternalIndex, F<T>>,
        e_cm: &F<T>,
    ) -> Result<()> {
        let mom_sum = externals.iter().zip(external_signature).fold(
            externals[ExternalIndex::from(0)].zero(),
            |acc, (mom, sign)| match sign {
                SignOrZero::Plus => acc + mom.clone(),
                SignOrZero::Minus => acc - mom.clone(),
                SignOrZero::Zero => unreachable!(),
            },
        );

        //println!("Total momentum after improvement: {:?}", mom_sum);
        if mom_sum.temporal.value.abs() > e_cm.epsilon() * e_cm {
            return Err(eyre::eyre!(
                "Energy is not conserved: {:+e} > {:+e}",
                mom_sum.temporal.value,
                e_cm.epsilon() * e_cm
            ));
        }

        if mom_sum.spatial.px.abs() > e_cm.epsilon() * e_cm {
            return Err(eyre::eyre!(
                "Momentum in x direction is not conserved: {:+e} > {:+e}",
                mom_sum.spatial.px,
                e_cm.epsilon() * e_cm
            ));
        }

        if mom_sum.spatial.py.abs() > e_cm.epsilon() * e_cm {
            return Err(eyre::eyre!(
                "Momentum in y direction is not conserved: {:+e} > {:+e}",
                mom_sum.spatial.py,
                e_cm.epsilon() * e_cm
            ));
        }

        if mom_sum.spatial.pz.abs() > e_cm.epsilon() * e_cm {
            return Err(eyre::eyre!(
                "Momentum in z direction is not conserved: {:+e} > {:+e}",
                mom_sum.spatial.pz,
                e_cm.epsilon() * e_cm
            ));
        }

        for (i, (mom, mass)) in externals.iter().zip(masses).enumerate() {
            let diff = (mom.square() - mass.square()).abs();
            let threshold = F::from_f64(2.0) * e_cm.epsilon() * e_cm * e_cm;
            if diff > threshold {
                return Err(eyre::eyre!(
                    "External {} is not on-shell: {:+e} > {:+e}",
                    i,
                    diff,
                    threshold
                ));
            }
        }

        Ok(())
    }

    #[test]
    fn test_photon_box() {
        test_initialise().unwrap();

        let photon_box: Graph = dot!(
            digraph photon_box {
                ext [style=invis];
                ext -> v1:0 [particle = "a", id=0];
                ext -> v2:1 [particle = "a", id=1];
                v3:2 -> ext [particle = "a", id=2];
                v4:3 -> ext [particle = "a", id=3];
                v1 -> v2 [particle = "t", id=4];
                v2 -> v3 [particle = "t", id=5];
                v3 -> v4 [particle = "t", id=6];
                v4 -> v1 [particle = "t", id=7];
            },
            "sm"
        )
        .unwrap();

        let sm = load_generic_model("sm");
        let e_cm = F(500.0);

        let photon_box_res = test_default_momenta_graph(&photon_box, &sm, &e_cm);
        if let Some(err) = photon_box_res.err() {
            panic!("Error in photon box test: {:?}", err);
        }
    }

    #[test]
    fn test_aa_tt() {
        test_initialise().unwrap();
        let sm = load_generic_model("sm");
        let e_cm = F(600.0);

        let aa_tt: Graph = dot!(
            digraph aa_tt {
                ext [style=invis];
                ext -> v1:0 [particle = "a", id=0];
                ext -> v2:1 [particle = "a", id=1];
                v1:2 -> ext [particle = "t", id=2];
                v2:3 -> ext [particle = "t~", id=3];
                v2 -> v1 [particle = "t", id=4];
            },
            "sm"
        )
        .unwrap();

        let aa_tt_res = test_default_momenta_graph(&aa_tt, &sm, &e_cm);
        if let Some(err) = aa_tt_res.err() {
            panic!("Error in aa_tt test: {:?}", err);
        }
    }

    #[test]
    fn test_gt_gt() {
        test_initialise().unwrap();
        let sm = load_generic_model("sm");
        let e_cm = F(700.0);

        let gt_gt: Graph = dot!(
            digraph gt_gt {
                ext [style=invis];
                ext -> v1:0 [particle = "g", id=0];
                ext -> v1:1 [particle = "t", id=1];
                v2:2 -> ext [particle = "g", id=2];
                v2:3 -> ext [particle = "t", id=3];
                v1 -> v2 [particle = "t", id=4];
            }, "sm"
        )
        .unwrap();

        let gt_gt_res = test_default_momenta_graph(&gt_gt, &sm, &e_cm);
        if let Some(err) = gt_gt_res.err() {
            panic!("Error in gt_gt test: {:?}", err);
        }
    }
}
