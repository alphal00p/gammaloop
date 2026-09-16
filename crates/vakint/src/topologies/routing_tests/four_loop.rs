use super::*;
use crate::{RustRedEvaluationError, RustRedEvaluationOptions, VakintSettings};

// Physical-slot fixtures corresponding to the four external K=10 family
// inputs. Names are deliberately absent: neither validation nor parent
// selection may depend on a topology label. Auxiliary slots are not physical
// matcher propagators and do not participate in this routing-only contract.
fn physical_parent_momenta() -> [Vec<Atom>; 4] {
    [
        vec![
            "k(1)",
            "k(2)",
            "k(3)",
            "k(4)",
            "k(1)-k(3)",
            "k(2)-k(3)",
            "k(3)-k(1)+k(4)",
            "k(3)-k(2)+k(4)",
            "k(3)+k(4)",
        ],
        vec![
            "k(1)",
            "k(2)",
            "k(3)",
            "k(4)",
            "k(1)-k(3)",
            "k(2)-k(3)",
            "k(3)-k(1)+k(4)",
            "k(3)-k(2)+k(4)",
            "k(3)-k(1)-k(2)+k(4)",
        ],
        vec![
            "k(1)",
            "k(2)",
            "k(3)",
            "k(4)",
            "k(1)-k(2)",
            "k(3)-k(4)",
            "k(2)+k(3)-k(1)",
            "k(3)-k(4)-k(1)",
        ],
        vec![
            "k(1)",
            "k(2)",
            "k(3)",
            "k(1)-k(3)",
            "k(4)",
            "k(2)-k(3)",
            "k(1)-k(3)+k(4)",
            "k(1)-k(2)",
        ],
    ]
    .map(|momenta| {
        momenta
            .into_iter()
            .map(|momentum| vk_parse!(momentum).unwrap())
            .collect()
    })
}

#[test]
fn registered_four_loop_census_validates_typed_parent_witnesses_without_admission() {
    let vakint = Vakint::new().unwrap();
    let parents = physical_parent_momenta();
    let identity = (1..=4)
        .map(|axis| function!(S.k, Atom::num(axis)))
        .collect::<Vec<_>>();
    let settings = VakintSettings {
        form_exe_path: "/this/path/must/not/be/invoked/by-routing-validation".to_owned(),
        ..VakintSettings::default()
    };
    let mut classes = 0;
    let mut slots = 0;
    let mut parent_counts = [0; 4];
    let mut nonidentity = 0;
    for registered in &vakint.topologies.0 {
        if !matches!(registered, Topology::FourLoop(_)) {
            continue;
        }
        let mut topology = registered.clone();
        let integral = topology.get_integral_mut();
        integral.name = "not-a-dispatch-key".into();
        let matches = parents
            .iter()
            .enumerate()
            .filter(|(_, expected)| integral.validate_parent_routing(expected).is_ok())
            .map(|(parent, _)| parent)
            .collect::<Vec<_>>();
        assert_eq!(matches.len(), 1);
        parent_counts[matches[0]] += 1;
        let expression = integral.canonical_expression.as_ref().unwrap();
        slots += (1..=integral.n_props)
            .filter(|&slot| get_prop_with_id(expression.as_view(), slot).is_some())
            .count();
        let (_, coordinates) = integral.parent_routing.as_ref().unwrap();
        assert_eq!(coordinates.len(), 4);
        nonidentity += usize::from(coordinates.as_ref() != identity.as_slice());
        classes += 1;

        // A validated witness is not a closing artifact. Even with master
        // substitution disabled, no four-loop RustRed backend is available.
        assert!(!crate::rustred_evaluation::supports(
            &settings,
            &topology,
            &RustRedEvaluationOptions {
                substitute_masters: false,
            },
        ));
    }
    assert_eq!(classes, 19);
    assert_eq!(slots, 123);
    assert_eq!(parent_counts, [1, 1, 1, 16]);
    assert_eq!(nonidentity, 6);
}

#[test]
fn typed_parent_validation_rejects_corrupted_four_loop_evidence() {
    let vakint = Vakint::new().unwrap();
    let parents = physical_parent_momenta();
    let original = vakint
        .topologies
        .0
        .iter()
        .filter(|topology| matches!(topology, Topology::FourLoop(_)))
        .map(Topology::get_integral)
        .find(|integral| {
            integral.validate_parent_routing(&parents[3]).is_ok()
                && Topologies::count_propagators_in_integral(
                    integral.canonical_expression.as_ref().unwrap().as_view(),
                ) < integral.n_props
        })
        .unwrap();
    let expected = &parents[3];
    let absent_slot = (1..=original.n_props)
        .find(|&slot| {
            get_prop_with_id(
                original.canonical_expression.as_ref().unwrap().as_view(),
                slot,
            )
            .is_none()
        })
        .unwrap();

    let mut missing = original.clone();
    missing.parent_routing = None;
    let mut short = original.clone();
    short.parent_routing.as_mut().unwrap().1 = Arc::from([vk_parse!("k(1)").unwrap()]);
    let mut wrong_coordinates = original.clone();
    let coordinates = &mut wrong_coordinates.parent_routing.as_mut().unwrap().1;
    let mut changed = coordinates.to_vec();
    changed[0] = &changed[0] + vk_parse!("k(1)").unwrap();
    *coordinates = changed.into();
    let mut wrong_parent = original.clone();
    let parent = &mut wrong_parent.parent_routing.as_mut().unwrap().0;
    *parent = parent
        .replace(
            vk_parse!(&format!(
                "prop({absent_slot},edge(left_,right_),q_,mass_,power_)"
            ))
            .unwrap(),
        )
        .with(
            vk_parse!(&format!(
                "prop({absent_slot},edge(left_,right_),2*q_,mass_,power_)"
            ))
            .unwrap(),
        );
    assert_ne!(parent, &original.parent_routing.as_ref().unwrap().0);
    let mut wrong_actual = original.clone();
    let actual = wrong_actual.canonical_expression.as_mut().unwrap();
    *actual = actual
        .replace(vk_parse!("k(1)").unwrap())
        .with(vk_parse!("2*k(1)").unwrap());
    for invalid in [
        missing,
        short,
        wrong_coordinates,
        wrong_parent,
        wrong_actual,
    ] {
        assert!(matches!(
            invalid.validate_parent_routing(expected),
            Err(RustRedEvaluationError::InvalidMatchedFamily { .. })
        ));
    }
    for invalid_descriptor in [Vec::new(), expected[..expected.len() - 1].to_vec()] {
        assert!(matches!(
            original.validate_parent_routing(&invalid_descriptor),
            Err(RustRedEvaluationError::InvalidMatchedFamily { .. })
        ));
    }
    // Squared denominators alone would accept this sign flip, but the shared
    // numerator-routing witness binds the signed physical momentum itself.
    let mut wrong_orientation = expected.clone();
    wrong_orientation[0] = -&wrong_orientation[0];
    assert!(matches!(
        original.validate_parent_routing(&wrong_orientation),
        Err(RustRedEvaluationError::InvalidMatchedFamily { .. })
    ));
}
