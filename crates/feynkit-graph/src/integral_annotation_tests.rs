use std::collections::BTreeMap;

use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
#[cfg(feature = "svg")]
use serde_json::Value;
use symbolica::atom::Atom;

#[cfg(feature = "svg")]
use crate::SceneOptions;
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

#[cfg(feature = "svg")]
#[test]
fn native_scene_displays_dots_without_an_automatic_title() {
    let diagram = master_diagram();
    let scene = diagram
        .to_scene(None, &Default::default(), None, &SceneOptions::default())
        .unwrap();
    assert!(scene.pages.iter().any(|label| {
        label.mathml.contains("<mn>2</mn>")
            && label
                .typst
                .as_ref()
                .is_some_and(|source| source.contains("^(2)"))
    }));
    assert!(scene.title.is_none());
    assert!(
        scene
            .edges
            .iter()
            .any(|edge| edge.details.get("propagator-power") == Some(&Value::from(3)))
    );
    assert!(
        scene
            .edges
            .iter()
            .any(|edge| edge.details.get("lmb-id") == Some(&Value::from(1)))
    );
}
