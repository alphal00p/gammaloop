use crate::{
    graph::{FeynmanGraph, Graph},
    model::Model,
    momentum::{FourMomentum, Sign, SignOrZero, ThreeMomentum},
    momentum_sample::ExternalIndex,
    settings::runtime::kinematic::Externals,
    utils::{newton_solver::newton_iteration_and_derivative, FloatLike, F},
    DependentMomentaConstructor,
};

use eyre::{Ok, Result};
use itertools::Itertools;
use symbolica::domains::float::{NumericalFloatLike, Real};
use typed_index_collections::TiVec;

pub fn improve_ps_psmc<T: FloatLike>(
    externals: Externals,
    e_cm: &F<T>,
    graph: &Graph,
    model: &Model,
) -> Result<TiVec<ExternalIndex, FourMomentum<F<T>>>> {
    let external_signature = graph.get_external_signature();
    let mom_constructor = DependentMomentaConstructor::Amplitude(&external_signature);
    let external_masses = graph.get_external_masses::<T>(model);

    let dependent_momenta = externals.get_dependent_externals::<T>(mom_constructor)?;

    let function = |x: &F<T>| {
        dependent_momenta
            .iter()
            .zip(&external_signature)
            .zip(&external_masses)
            .map(|((mom, sign), mass)| {
                let rescaled_mom = match sign {
                    SignOrZero::Minus => mom.rescale_spatial(x),
                    _ => mom.clone(),
                };
                let ose = rescaled_mom.spatial.on_shell_energy(Some(mass.clone()));

                let signed_ose = match sign {
                    SignOrZero::Plus => ose.value.clone(),
                    SignOrZero::Minus => -&ose.value,
                    SignOrZero::Zero => unreachable!(),
                };

                let der = match sign {
                    SignOrZero::Minus => -x * mom.spatial.norm_squared() / &ose.value,
                    SignOrZero::Plus => F::<T>::from_f64(0.0),
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
        &F::<T>::from_f64(1.0),
        20,
        &e_cm,
    );

    Ok(dependent_momenta
        .iter()
        .zip(&external_signature)
        .zip(&external_masses)
        .map(|((mom, sign), mass)| {
            let mut rescaled_mom = match sign {
                SignOrZero::Minus => mom.rescale_spatial(&solution.solution),
                _ => mom.clone(),
            };

            let ose = rescaled_mom.spatial.on_shell_energy(Some(mass.clone()));
            rescaled_mom.temporal = ose;
            rescaled_mom
        })
        .collect::<TiVec<ExternalIndex, FourMomentum<F<T>>>>())
}

fn improve_ps_vh<T: FloatLike>(
    externals: Externals,
    e_cm: &F<T>,
    graph: &Graph,
    model: &Model,
) -> Result<TiVec<ExternalIndex, FourMomentum<F<T>>>> {
    let external_signature = graph.get_external_signature();
    let mom_constructor = DependentMomentaConstructor::Amplitude(&external_signature);
    let external_masses = graph.get_external_masses::<T>(model);
    let dependent_momenta = externals.get_dependent_externals::<T>(mom_constructor)?;

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
    use typed_index_collections::TiVec;

    use crate::{
        dot,
        graph::{parse::IntoGraph, FeynmanGraph, Graph},
        initialisation::test_initialise,
        momentum::{Dep, ExternalMomenta, FourMomentum, Sign, SignOrZero},
        momentum_sample::ExternalIndex,
        settings::runtime::kinematic::Externals,
        signature::SignatureLike,
        utils::{f128, test_utils::load_generic_model, FloatLike, F},
        DependentMomentaConstructor,
    };

    fn test_kinematic_validity<T: FloatLike>(
        externals: TiVec<ExternalIndex, FourMomentum<F<T>>>,
        external_signature: SignatureLike<ExternalIndex>,
        masses: TiVec<ExternalIndex, F<T>>,
    ) {
        let mom_sum = externals.iter().zip(&external_signature).fold(
            externals[ExternalIndex::from(0)].zero(),
            |acc, (mom, sign)| match sign {
                SignOrZero::Plus => acc + mom.clone(),
                SignOrZero::Minus => acc - mom.clone(),
                SignOrZero::Zero => unreachable!(),
            },
        );

        println!("Sum of momenta: {:+16e}", mom_sum);

        for (mom, mass) in externals.iter().zip(masses) {
            let p2 = mom.square();
            let mass2 = &mass * &mass;
            let diff = (p2 - mass2).abs();
            println!("on-shell ness: {:+16e}", diff);
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
                ExternalMomenta::Independent([F(10.0), F(0.0), F(0.0), F(10.0)].into()),
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

        println!("Before improvement:");
        test_kinematic_validity(
            dependent_momenta.clone(),
            external_signature.clone(),
            external_masses.clone(),
        );

        let e_cm = F(20.0);
        let improved_momenta =
            super::improve_ps_psmc(external.clone(), &e_cm, &graph, &sm).unwrap();

        println!("After improvement:");
        test_kinematic_validity(
            improved_momenta,
            external_signature.clone(),
            external_masses.clone(),
        );

        let vh_improved_momenta =
            super::improve_ps_vh(external.clone(), &e_cm, &graph, &sm).unwrap();

        println!("vh improved momenta: {:?}", vh_improved_momenta);
        println!("After vh improvement:");
        test_kinematic_validity(vh_improved_momenta, external_signature, external_masses);
        println!("eps: {:+16e}", e_cm.epsilon());
    }

    #[test]
    fn buh() {
        let test: F<f128> = F(0.0324).higher();

        println!("test: {:?}", test);
    }
}
