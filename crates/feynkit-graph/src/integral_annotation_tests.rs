use std::collections::BTreeMap;

use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use symbolica::atom::{Atom, AtomCore};

use crate::{DiagramError, EdgeId, FeynmanDiagram, NumeratorInterface};

fn master_diagram() -> FeynmanDiagram {
    FeynmanDiagram::from_dot(
        Model::phi3(),
        r#"digraph master {
            graph [num="(x+y)^2*(x-y)^3"];
            a -> b [particle="phi", power=2, lmb_id=0];
            a -> b [particle="phi", power=3, lmb_id=1];
            b -> a [particle="phi", num="7"];
        }"#,
    )
    .unwrap()
}

#[test]
fn powers_and_factored_numerators_survive_semantic_round_trips() {
    let diagram = master_diagram()
        .with_loop_momentum_edges(&[EdgeId(1), EdgeId(0)])
        .unwrap();
    let mut serialized = diagram.serde_view();
    serialized.half_edge_order.rotate_left(2);
    let diagram = FeynmanDiagram::from_serde(diagram.model_arc(), serialized).unwrap();
    let json = diagram.to_json().unwrap();
    let dot = diagram.to_dot().unwrap();
    assert!(dot.contains("power=2"));
    assert!(dot.contains("lmb_id=1"));
    for restored in [
        FeynmanDiagram::from_json(diagram.model_arc(), &json).unwrap(),
        FeynmanDiagram::from_dot(diagram.model_arc(), &dot).unwrap(),
    ] {
        assert_eq!(restored.to_json().unwrap(), json);
        assert!(
            restored
                .numerator_prefactor()
                .same_value(diagram.numerator_prefactor())
        );
        assert_eq!(restored.numerator().expression(), &Atom::num(7));
        assert_eq!(
            restored.loop_momentum_basis(),
            diagram.loop_momentum_basis()
        );
        assert_eq!(
            restored.half_edges().collect::<Vec<_>>(),
            diagram.half_edges().collect::<Vec<_>>()
        );
        assert_eq!(
            restored
                .edges()
                .map(|(_, _, edge)| edge.propagator_power)
                .collect::<Vec<_>>(),
            [2, 3, 1]
        );
    }
}

#[test]
fn physical_power_is_used_by_denominators_but_not_family_forms() {
    let dotted = master_diagram();
    let unit = dotted
        .map_data(
            |_, vertex| vertex.clone(),
            |_, _, edge| {
                let mut edge = edge.clone();
                edge.propagator_power = 1;
                edge
            },
        )
        .unwrap();
    assert_eq!(dotted.id(), unit.id());
    assert_eq!(
        dotted.canonical_key().unwrap(),
        unit.canonical_key().unwrap()
    );
    let explicit = BTreeMap::from([(EdgeId(0), 2), (EdgeId(1), 3)]);
    assert_eq!(
        dotted.denominator_expression().unwrap(),
        unit.denominator_of(&unit.internal_subgraph(), &explicit)
            .unwrap()
    );
    assert_eq!(
        dotted.superficial_degree_of_divergence(4).unwrap(),
        unit.superficial_degree_of_divergence(4).unwrap() - 6
    );
    let kinematics = Kinematics::new();
    assert_eq!(
        dotted.ordered_propagators(&kinematics).unwrap(),
        unit.ordered_propagators(&kinematics).unwrap()
    );
}

#[test]
fn nonpositive_powers_and_inconsistent_loop_labels_are_rejected() {
    for power in [0, -1] {
        let dot = format!(
            "digraph {{ a -> b [particle=phi, power={power}]; a -> b [particle=phi]; b -> a [particle=phi]; }}"
        );
        assert!(matches!(
            FeynmanDiagram::from_dot(Model::phi3(), &dot),
            Err(DiagramError::InvalidPropagatorPower { .. })
        ));
    }
    let diagram = master_diagram();
    let invalid_id = diagram
        .to_dot()
        .unwrap()
        .replace("id=0, particle_id", "id=9, particle_id");
    assert!(FeynmanDiagram::from_dot(diagram.model_arc(), &invalid_id).is_err());
    let dot = diagram.to_dot().unwrap().replace("lmb_id=1", "lmb_id=0");
    assert!(matches!(
        FeynmanDiagram::from_dot(diagram.model_arc(), &dot),
        Err(DiagramError::InvalidLoopMomentumBasis(_))
    ));
}

#[test]
fn stored_power_defaults_to_one_for_existing_json() {
    let diagram = master_diagram();
    // Retain the integer-valued diagram ID exactly; generic JSON Value would
    // pass its u128 through floating point when arbitrary_precision is disabled.
    let json = diagram
        .to_json()
        .unwrap()
        .lines()
        .filter(|line| !line.contains("\"propagator_power\""))
        .collect::<Vec<_>>()
        .join("\n");
    let restored = FeynmanDiagram::from_json(diagram.model_arc(), &json).unwrap();
    assert!(
        restored
            .edges()
            .all(|(_, _, edge)| edge.propagator_power == 1)
    );
    assert!(
        NumeratorInterface::export(restored.numerator_prefactor()).0
            == NumeratorInterface::export(diagram.numerator_prefactor()).0
    );
}

#[test]
fn signed_overrides_replace_only_the_requested_stored_powers() {
    let diagram = master_diagram();
    let unit_denominators = (0..3)
        .map(|edge| {
            let powers = (0..3)
                .map(|other| (EdgeId(other), if other == edge { 1 } else { 0 }))
                .collect();
            diagram
                .denominator_of(&diagram.internal_subgraph(), &powers)
                .unwrap()
        })
        .collect::<Vec<_>>();
    for power in [-3, 0, 1, 4] {
        let actual = diagram
            .denominator_of(
                &diagram.internal_subgraph(),
                &BTreeMap::from([(EdgeId(0), power)]),
            )
            .unwrap();
        assert_eq!(
            actual,
            unit_denominators[0].pow(power as i64)
                * unit_denominators[1].pow(3)
                * &unit_denominators[2],
        );
    }
}

#[test]
fn uv_expansion_uses_stored_powers_and_signed_overrides() {
    let diagram = FeynmanDiagram::from_dot(
        Model::phi3(),
        "digraph { ext [style=invis]; a -> b [particle=phi, power=2]; b -> a [particle=phi]; ext -> a [particle=phi]; b -> ext [particle=phi]; }",
    )
    .unwrap();
    let region = diagram.internal_subgraph();
    let mass = Atom::num(3);
    let stored = diagram
        .uv_expansion_of(&region, &mass, 4, None, &BTreeMap::new())
        .unwrap();
    assert_eq!(stored, Atom::Zero);
    assert_eq!(
        stored,
        diagram
            .uv_expansion_of(&region, &mass, 4, None, &BTreeMap::from([(EdgeId(0), 2)]),)
            .unwrap()
    );
    let unit = diagram
        .map_data(
            |_, vertex| vertex.clone(),
            |_, _, edge| {
                let mut edge = edge.clone();
                edge.propagator_power = 1;
                edge
            },
        )
        .unwrap();
    for power in [-1, 0, 1] {
        let overrides = BTreeMap::from([(EdgeId(0), power)]);
        let actual = diagram
            .uv_expansion_of(&region, &mass, 4, None, &overrides)
            .unwrap();
        assert_ne!(actual, Atom::Zero);
        assert_eq!(
            actual,
            unit.uv_expansion_of(&region, &mass, 4, None, &overrides)
                .unwrap()
        );
    }
}

#[test]
fn invalid_powers_are_rejected_at_every_import_boundary() {
    let diagram = master_diagram();
    for power in ["0", "-1", "2147483648", "4294967297", "1.5", "x"] {
        let compact = format!(
            "digraph {{ a -> b [particle=phi, power=\"{power}\"]; a -> b [particle=phi]; b -> a [particle=phi]; }}"
        );
        assert!(
            FeynmanDiagram::from_dot(Model::phi3(), &compact).is_err(),
            "{power}"
        );
        let dot = diagram
            .to_dot()
            .unwrap()
            .replace("power=2", &format!("power=\"{power}\""));
        assert!(
            FeynmanDiagram::from_dot(diagram.model_arc(), &dot).is_err(),
            "{power}"
        );
        let json = diagram.to_json().unwrap().replace(
            "\"propagator_power\": 2",
            &format!("\"propagator_power\": {power}"),
        );
        assert!(
            FeynmanDiagram::from_json(diagram.model_arc(), &json).is_err(),
            "{power}"
        );
    }
    for power in [0, -1] {
        let invalid = diagram
            .map_data(
                |_, vertex| vertex.clone(),
                |_, _, edge| {
                    let mut edge = edge.clone();
                    edge.propagator_power = power;
                    edge
                },
            )
            .unwrap();
        assert!(matches!(
            invalid.validate(),
            Err(DiagramError::InvalidPropagatorPower { .. })
        ));
        assert!(matches!(
            invalid.denominator_expression(),
            Err(DiagramError::InvalidPropagatorPower { .. })
        ));
        for degree in [
            invalid.superficial_degree_of_divergence(4),
            invalid.superficial_degree_of_divergence_of(&invalid.internal_subgraph(), 4),
        ] {
            assert!(matches!(
                degree,
                Err(DiagramError::InvalidPropagatorPower { power: actual, .. }) if actual == power
            ));
        }
        assert!(matches!(
            invalid.uv_expansion_of(
                &invalid.internal_subgraph(),
                &Atom::num(3),
                4,
                None,
                &BTreeMap::new()
            ),
            Err(DiagramError::InvalidPropagatorPower { .. })
        ));
    }
    let cloned = diagram.clone();
    assert_eq!(cloned.to_json().unwrap(), diagram.to_json().unwrap());
    let json = diagram.to_json().unwrap().replace(
        "\"propagator_power\": 2",
        "\"propagator_power\": 2147483647",
    );
    let largest = FeynmanDiagram::from_json(diagram.model_arc(), &json).unwrap();
    assert_eq!(largest.edges().next().unwrap().2.propagator_power, i32::MAX);
    assert!(largest.superficial_degree_of_divergence(4).is_err());
}

#[test]
fn sewing_validates_both_source_powers_before_discarding_a_leg() {
    for outgoing in [0, -1, 3] {
        let dot = format!(
            "digraph {{ graph [final_state=\"phi,phi\"]; ext [style=invis]; ext -> a [particle=phi, is_cut=9, power=2]; b -> ext [particle=phi, is_cut=9, power={outgoing}]; a -> b [particle=phi]; b -> a [particle=phi]; }}"
        );
        let error = FeynmanDiagram::from_dot(Model::phi3(), &dot).unwrap_err();
        if outgoing <= 0 {
            assert!(
                matches!(error, DiagramError::InvalidPropagatorPower { power, .. } if power == outgoing)
            );
        } else {
            assert!(matches!(
                error,
                DiagramError::InconsistentPropagatorPowers {
                    incoming_power: 2,
                    outgoing_power: 3,
                    ..
                }
            ));
        }
    }
    let diagram = FeynmanDiagram::from_dot(Model::phi3(), "digraph { graph [final_state=\"phi,phi\"]; ext [style=invis]; ext -> a [particle=phi, is_cut=9, power=2]; b -> ext [particle=phi, is_cut=9, power=2]; a -> b [particle=phi]; b -> a [particle=phi]; }").unwrap();
    assert_eq!(
        diagram
            .edges()
            .find(|(_, _, edge)| edge.external.is_some())
            .unwrap()
            .2
            .propagator_power,
        2
    );
    assert_eq!(
        FeynmanDiagram::from_json(diagram.model_arc(), &diagram.to_json().unwrap())
            .unwrap()
            .to_json()
            .unwrap(),
        diagram.to_json().unwrap()
    );
}
