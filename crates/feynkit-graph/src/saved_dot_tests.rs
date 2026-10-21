//! Saved integral diagrams retain their numerical meaning, not just topology.

use symbolica::{
    atom::{Atom, AtomCore},
    symbol,
};

use super::*;

fn four_loop_integral() -> FeynmanDiagram {
    let model = crate::tests::scalar_model();
    let particle = model.particle_id("phi").unwrap();
    let numerator = Atom::var(symbol!("saved_dot_test::n")).pow(2);
    let mut builder =
        FeynmanDiagram::builder(model, "dotted four-loop master").numerator(numerator.clone());
    let left = builder.add_vertex(DiagramVertex {
        name: "left".into(),
        interaction: None,
        external: None,
        numerator,
    });
    let right = builder.add_vertex(DiagramVertex {
        name: "right".into(),
        interaction: None,
        external: None,
        numerator: Atom::one(),
    });
    for index in 0..5 {
        let mut edge = DiagramEdge::new(particle, false);
        if index == 0 {
            edge.propagator_power = 2;
        }
        builder.add_edge(left, right, edge).unwrap();
    }
    let mut diagram = builder
        .build()
        .unwrap()
        .with_loop_momentum_edges(&[EdgeId(2), EdgeId(3), EdgeId(0), EdgeId(1)])
        .unwrap();
    let mut order = diagram.half_edge_order();
    order.rotate_left(4);
    order.swap(1, 7);
    diagram.restore_half_edge_order(&order).unwrap();
    diagram.validate().unwrap();
    diagram
}

#[test]
fn saved_integral_round_trip_retains_powers_routing_numerator_and_half_edges() {
    let diagram = four_loop_integral();
    let dot = diagram.to_dot().unwrap();
    assert!(dot.contains("power=2"));
    assert!(dot.contains("lmb_id=3"));
    let restored = FeynmanDiagram::from_dot(diagram.model_arc(), &dot).unwrap();
    assert_eq!(restored.to_json().unwrap(), diagram.to_json().unwrap());
    assert_eq!(restored.half_edge_order(), diagram.half_edge_order());
    assert_eq!(
        restored.loop_momentum_basis(),
        diagram.loop_momentum_basis()
    );
    assert_eq!(restored.numerator(), diagram.numerator());
    assert_eq!(restored.edges().next().unwrap().2.propagator_power, 2);
    let from_json =
        FeynmanDiagram::from_json(diagram.model_arc(), &diagram.to_json().unwrap()).unwrap();
    assert_eq!(from_json.to_dot().unwrap(), dot);
    let mapped = diagram
        .map_data(|_, vertex| vertex.clone(), |_, _, edge| edge.clone())
        .unwrap();
    assert_eq!(mapped.to_json().unwrap(), diagram.to_json().unwrap());
}

#[test]
fn rejects_corrupt_stored_routing_and_inconsistent_loop_labels() {
    let diagram = four_loop_integral();
    let dot = diagram.to_dot().unwrap();
    let mut invalid_basis = serde_json::to_value(diagram.loop_momentum_basis()).unwrap();
    invalid_basis["edge_signatures"]["4"]["loops"][0] = 0.into();
    let corrupt = dot.replace(
        &FeynmanDiagram::dot_string(&serde_json::to_string(diagram.loop_momentum_basis()).unwrap())
            .unwrap(),
        &FeynmanDiagram::dot_string(&serde_json::to_string(&invalid_basis).unwrap()).unwrap(),
    );
    assert!(matches!(
        FeynmanDiagram::from_dot(diagram.model_arc(), &corrupt),
        Err(DiagramError::InvalidLoopMomentumBasis(_))
    ));
    for corrupt in [
        dot.replace("lmb_id=0", "lmb_id=4"),
        dot.replace("power=2", "power=0"),
        dot.replace("id=0, particle_id", "id=9, particle_id"),
    ] {
        assert!(FeynmanDiagram::from_dot(diagram.model_arc(), &corrupt).is_err());
    }
    let mut order = diagram.half_edge_order();
    order[1] = order[0];
    let corrupt = dot.replace(
        &FeynmanDiagram::dot_string(&serde_json::to_string(&diagram.half_edge_order()).unwrap())
            .unwrap(),
        &FeynmanDiagram::dot_string(&serde_json::to_string(&order).unwrap()).unwrap(),
    );
    assert!(FeynmanDiagram::from_dot(diagram.model_arc(), &corrupt).is_err());
}

#[test]
fn metadata_mapping_preserves_nonascending_tree_and_loop_orders() {
    let model = crate::tests::scalar_model();
    let particle = model.particle_id("phi").unwrap();
    let mut builder = FeynmanDiagram::builder(model, "ordered forest");
    let vertices = (0..3)
        .map(|i| {
            builder.add_vertex(DiagramVertex {
                name: format!("v{i}"),
                interaction: None,
                external: None,
                numerator: Atom::one(),
            })
        })
        .collect::<Vec<_>>();
    for (source, target) in [(0, 1), (0, 2), (1, 2), (0, 1), (1, 2), (0, 2)] {
        builder
            .add_edge(
                vertices[source],
                vertices[target],
                DiagramEdge::new(particle, false),
            )
            .unwrap();
    }
    let mut diagram = builder
        .build()
        .unwrap()
        .with_loop_momentum_tree_edges(&[EdgeId(4), EdgeId(1)])
        .unwrap();
    diagram.loop_momentum_basis = diagram
        .loop_momentum_basis
        .with_loop_edge_order(&[EdgeId(5), EdgeId(0), EdgeId(3), EdgeId(2)])
        .unwrap();
    diagram.validate().unwrap();
    let mapped = diagram
        .map_data(|_, vertex| vertex.clone(), |_, _, edge| edge.clone())
        .unwrap();
    assert_eq!(mapped.to_json().unwrap(), diagram.to_json().unwrap());
    let restored =
        FeynmanDiagram::from_dot(diagram.model_arc(), &diagram.to_dot().unwrap()).unwrap();
    assert_eq!(restored.to_json().unwrap(), diagram.to_json().unwrap());
}

#[test]
fn explicitly_missing_routing_signatures_are_rejected_instead_of_recomputed() {
    let diagram = four_loop_integral()
        .with_loop_momentum_edges(&[EdgeId(1), EdgeId(2), EdgeId(3), EdgeId(4)])
        .unwrap();
    let mut invalid_basis = serde_json::to_value(diagram.loop_momentum_basis()).unwrap();
    invalid_basis["edge_signatures"] = serde_json::json!({});
    let corrupt = diagram.to_dot().unwrap().replace(
        &FeynmanDiagram::dot_string(&serde_json::to_string(diagram.loop_momentum_basis()).unwrap())
            .unwrap(),
        &FeynmanDiagram::dot_string(&serde_json::to_string(&invalid_basis).unwrap()).unwrap(),
    );
    assert!(matches!(
        FeynmanDiagram::from_dot(diagram.model_arc(), &corrupt),
        Err(DiagramError::InvalidLoopMomentumBasis(_))
    ));
}
