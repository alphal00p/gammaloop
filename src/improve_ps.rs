use crate::{
    graph::{FeynmanGraph, Graph},
    model::Model,
    momentum::{Energy, FourMomentum, Sign, SignOrZero, ThreeMomentum},
    momentum_sample::ExternalIndex,
    settings::runtime::kinematic::Externals,
    signature::SignatureLike,
    status_debug,
    utils::{newton_solver::newton_iteration_and_derivative, FloatLike, F},
    DependentMomentaConstructor,
};

use eyre::{Ok, Result};
use itertools::Itertools;
use symbolica::domains::float::{NumericalFloatLike, Real};
use typed_index_collections::TiVec;

pub struct PhaseSpaceImprovementSettings {
    pub mode: ImprovementMode,
    pub large_deformation_failure: Option<F<f64>>,
}

pub enum ImprovementMode {
    Psmc,
    Vh,
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
        ImprovementMode::Psmc => {
            improve_ps_psmc(dependent_momenta, external_masses, external_signature, e_cm)
        }
        ImprovementMode::Vh => {
            improve_ps_vh(dependent_momenta, external_masses, external_signature, e_cm)
        }
    }?;

    let deformation_size = dimensionless_metric(dependent_momenta, &result, e_cm);

    if let Some(threshold) = settings.large_deformation_failure {
        if deformation_size > F::from_ff64(threshold) {
            return Err(eyre::eyre!(
                "Phase space improvement resulted in a large deformation: {:+e} > {:+e}",
                deformation_size,
                threshold
            ));
        }
    }

    Ok(result)
}

fn improve_ps_psmc<T: FloatLike>(
    dependent_momenta: &TiVec<ExternalIndex, FourMomentum<F<T>>>,
    external_masses: &TiVec<ExternalIndex, F<T>>,
    external_signature: &SignatureLike<ExternalIndex>,
    e_cm: &F<T>,
) -> Result<TiVec<ExternalIndex, FourMomentum<F<T>>>> {
    let function = |x: &F<T>| {
        dependent_momenta
            .iter()
            .enumerate()
            .zip(external_signature)
            .zip(external_masses)
            .map(|(((id, mom), sign), mass)| {
                let rescaled_mom = if id == 0 {
                    mom.rescale_spatial(x)
                } else {
                    mom.rescale_spatial(x)
                };

                let ose = rescaled_mom.spatial.on_shell_energy(Some(mass.clone()));
                let der = if id == 0 {
                    x * mom.spatial.norm_squared() / &ose.value
                } else {
                    x * mom.spatial.norm_squared() / &ose.value
                };

                let signed_ose = match sign {
                    SignOrZero::Plus => ose.value.clone(),
                    SignOrZero::Minus => -&ose.value,
                    SignOrZero::Zero => unreachable!(),
                };

                let signed_der = match sign {
                    SignOrZero::Plus => der,
                    SignOrZero::Minus => -&der,
                    SignOrZero::Zero => unreachable!(),
                };

                (signed_ose, signed_der)
            })
            .reduce(|(ose_a, der_a), (ose_b, der_b)| ((ose_a + ose_b), (der_a + der_b)))
            .unwrap_or((F::from_f64(0.0), F::from_f64(0.0)))
    };

    let solution = newton_iteration_and_derivative(
        &F::from_f64(1.0),
        function,
        &F::<T>::from_f64(1.000),
        400,
        &e_cm,
    );

    status_debug!("solution: {:?}", solution);

    Ok(dependent_momenta
        .iter()
        .enumerate()
        .zip(external_signature)
        .zip(external_masses)
        .map(|(((id, mom), sign), mass)| {
            let mut rescaled_mom = if id == 0 {
                mom.rescale_spatial(&solution.solution)
            } else {
                mom.rescale_spatial(&solution.solution)
            };

            let ose = rescaled_mom.spatial.on_shell_energy(Some(mass.clone()));
            rescaled_mom.temporal = ose;
            rescaled_mom
        })
        .collect::<TiVec<ExternalIndex, FourMomentum<F<T>>>>())
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

    let cos_theta_to_z_axis = initial_state_sum.spatial.get_cos_theta_with(&z_axis);

    let rotated_dependent_momenta = dependent_momenta
        .iter()
        .map(|m| {
            let mut new_momenta = FourMomentum {
                temporal: m.temporal.clone(),
                spatial: m.spatial.axis_angle_rotation(&cos_theta_to_z_axis, &z_axis),
            };
            // force px and py to be zero:
            new_momenta.spatial.px = e_cm.zero();
            new_momenta.spatial.py = e_cm.zero();
            new_momenta
        })
        .collect::<TiVec<ExternalIndex, FourMomentum<F<T>>>>();

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

    let boosted_momenta = rotated_dependent_momenta
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

    let function = |x: &F<T>| {
        boosted_momenta
            .iter()
            .enumerate()
            .zip(external_signature)
            .zip(external_masses)
            .map(|(((_id, mom), sign), mass)| {
                let rescaled_mom = match sign {
                    SignOrZero::Plus => mom.clone(),
                    SignOrZero::Minus => mom.rescale_spatial(x),
                    SignOrZero::Zero => unreachable!(),
                };

                let ose = rescaled_mom.spatial.on_shell_energy(Some(mass.clone()));
                let der = match sign {
                    SignOrZero::Plus => e_cm.zero(),
                    SignOrZero::Minus => -x * mom.spatial.norm_squared() / &ose.value,
                    SignOrZero::Zero => unreachable!(),
                };

                let signed_ose = match sign {
                    SignOrZero::Plus => ose.value.clone(),
                    SignOrZero::Minus => -&ose.value,
                    SignOrZero::Zero => unreachable!(),
                };

                (signed_ose, der)
            })
            .reduce(|(ose_a, der_a), (ose_b, der_b)| ((ose_a + ose_b), (der_a + der_b)))
            .unwrap_or((F::from_f64(0.0), F::from_f64(0.0)))
    };

    let solution = newton_iteration_and_derivative(
        &F::from_f64(1.0),
        function,
        &F::<T>::from_f64(1.000),
        400,
        &e_cm,
    );

    status_debug!("solution: {:?}", solution);

    Ok(boosted_momenta
        .iter()
        .enumerate()
        .zip(external_signature)
        .zip(external_masses)
        .map(|(((id, mom), sign), mass)| {
            let mut rescaled_mom = match sign {
                SignOrZero::Plus => mom.clone(),
                SignOrZero::Minus => mom.rescale_spatial(&solution.solution),
                SignOrZero::Zero => unreachable!(),
            };

            let ose = rescaled_mom.spatial.on_shell_energy(Some(mass.clone()));
            rescaled_mom.temporal = ose;
            rescaled_mom
        })
        .collect::<TiVec<ExternalIndex, FourMomentum<F<T>>>>())
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
        .map(|m| m.clone())
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

    use symbolica::create_hyperdual_single_derivative;
    use symbolica::domains::float::NumericalFloatLike;
    use symbolica::domains::rational::Rational;
    use typed_index_collections::TiVec;

    use crate::{
        dot,
        graph::{ext, parse::IntoGraph, FeynmanGraph, Graph},
        initialisation::test_initialise,
        momentum::{Dep, ExternalMomenta, FourMomentum, SignOrZero},
        momentum_sample::ExternalIndex,
        settings::runtime::kinematic::Externals,
        signature::SignatureLike,
        utils::{f128, test_utils::load_generic_model, FloatLike, F},
        DependentMomentaConstructor,
    };

    fn test_kinematic_validity<T: FloatLike>(
        externals: &TiVec<ExternalIndex, FourMomentum<F<T>>>,
        external_signature: &SignatureLike<ExternalIndex>,
        masses: &TiVec<ExternalIndex, F<T>>,
        e_cm: &F<T>,
    ) {
        let mom_sum = externals.iter().zip(external_signature).fold(
            externals[ExternalIndex::from(0)].zero(),
            |acc, (mom, sign)| match sign {
                SignOrZero::Plus => acc + mom.clone(),
                SignOrZero::Minus => acc - mom.clone(),
                SignOrZero::Zero => unreachable!(),
            },
        );

        println!("Total momentum after improvement: {:?}", mom_sum);
        assert!(mom_sum.temporal.value.abs() < e_cm.epsilon() * e_cm);
        assert!(mom_sum.spatial.px.abs() < e_cm.epsilon() * e_cm);
        assert!(mom_sum.spatial.py.abs() < e_cm.epsilon() * e_cm);
        assert!(mom_sum.spatial.pz.abs() < e_cm.epsilon() * e_cm);

        for (mom, mass) in externals.iter().zip(masses) {
            let diff = (mom.square() - mass.square()).abs();
            assert!(diff < e_cm.epsilon() * e_cm.square());
        }
    }

    #[test]
    fn test_massless_photon_box() {
        test_initialise().unwrap();

        let graph: Graph = dot!(
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

        let external = Externals::Constant {
            momenta: vec![
                ExternalMomenta::Independent([F(10.0), F(0.0), F(0.00), F(10.0)].into()),
                ExternalMomenta::Independent([F(10.0), F(0.0), F(0.0), F(-10.0)].into()),
                ExternalMomenta::Independent([F(5.0), F(0.0), F(4.0), F(3.0)].into()),
                ExternalMomenta::Dependent(Dep::Dep),
            ],
            helicities: vec![
                SignOrZero::Plus,
                SignOrZero::Plus,
                SignOrZero::Minus,
                SignOrZero::Minus,
            ],
        };

        let external_signature = graph.get_external_signature();
        let external_masses = graph.get_external_masses::<f64>(&sm);

        let constructor = DependentMomentaConstructor::Amplitude(&external_signature);

        let dependent_momenta = external
            .get_dependent_externals::<f64>(constructor)
            .unwrap();

        let e_cm = F(20.0);
        let improved_momenta = super::improve_ps_psmc(
            &dependent_momenta,
            &external_masses,
            &external_signature,
            &e_cm,
        )
        .unwrap();

        println!("After improvement:");
        test_kinematic_validity(
            &improved_momenta,
            &external_signature,
            &external_masses,
            &e_cm,
        );

        let vh_improved_momenta = super::improve_ps_vh(
            &dependent_momenta,
            &external_masses,
            &external_signature,
            &e_cm,
        )
        .unwrap();

        test_kinematic_validity(
            &vh_improved_momenta,
            &external_signature,
            &external_masses,
            &e_cm,
        );
    }

    #[test]
    fn test_aa_tt() {
        test_initialise().unwrap();
        let graph: Graph = dot!(
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

        let sm = load_generic_model("sm");
        let externals = Externals::Constant {
            momenta: vec![
                ExternalMomenta::Independent([F(170.0), F(0.0), F(0.0), F(-170.0)].into()),
                ExternalMomenta::Independent([F(180.0), F(0.0), F(0.0), F(180.0)].into()),
                ExternalMomenta::Independent(
                    [
                        F((173.0 * 173.0 + 134.0 * 134.0f64).sqrt()),
                        F(0.0),
                        F(0.0),
                        F(134.0),
                    ]
                    .into(),
                ),
                ExternalMomenta::Dependent(Dep::Dep),
            ],
            helicities: vec![
                SignOrZero::Plus,
                SignOrZero::Plus,
                SignOrZero::Minus,
                SignOrZero::Minus,
            ],
        };

        let external_signature = graph.get_external_signature();
        let external_masses = graph.get_external_masses::<f64>(&sm);
        let constructor = DependentMomentaConstructor::Amplitude(&external_signature);

        let dependent_momenta = externals
            .get_dependent_externals::<f64>(constructor)
            .unwrap();

        let e_cm = F(440.0);
        let improved_momenta = super::improve_ps_mf(
            &dependent_momenta,
            &external_masses,
            &external_signature,
            &e_cm,
        )
        .unwrap();

        test_kinematic_validity(
            &improved_momenta,
            &external_signature,
            &external_masses,
            &e_cm,
        );
    }

    #[test]
    fn buh() {
        let test: F<f128> = F(0.0324).higher();

        println!("test: {:?}", test);
    }
}
