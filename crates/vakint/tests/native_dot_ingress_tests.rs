#![cfg(feature = "feynkit-ingress")]

use feynkit_graph::FeynmanDiagram;
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use symbolica::atom::{Atom, AtomCore};
use vakint::diagram_integral::DiagramIntegralOptions;
use vakint::{
    RustRedMasterMode, TensorReductionMethod, Vakint, VakintExpression, VakintSettings,
    vakint_parse,
};

#[path = "support/rustred.rs"]
#[allow(dead_code)]
mod rustred_support;

/// Symbolic RustRed after the native FeynKit tensor reduction, without FORM.
fn symbolic_rustred() -> VakintSettings {
    VakintSettings {
        tensor_reduction_method: TensorReductionMethod::FeynKit,
        ..rustred_support::rustred_only(
            RustRedMasterMode::Symbolic,
            "/this/executable/must/never/run",
        )
    }
}

#[test]
fn native_dot_diagram_reaches_symbolic_rustred() {
    let mut definition: serde_json::Value =
        serde_json::from_str(&Model::phi3().to_json().unwrap()).unwrap();
    definition["name"] = "two_point_insertion".into();
    definition["vertex_rules"][0]["particles"]
        .as_array_mut()
        .unwrap()
        .truncate(2);
    definition["lorentz_structures"][0]["spins"]
        .as_array_mut()
        .unwrap()
        .truncate(2);
    let diagram = FeynmanDiagram::from_dot(
        Model::from_json(&definition.to_string()).unwrap(),
        r#"digraph tadpole { a -> a [particle="phi"]; }"#,
    )
    .unwrap();
    let family = diagram.propagator_family(&Kinematics::new()).unwrap();
    let settings = symbolic_rustred();
    let momentum = &family.loop_momenta()[0];
    let mass_squared = (family
        .kinematics()
        .scalar_product(momentum, momentum)
        .unwrap()
        - &family.denominators()[0])
        .expand();
    let engine = Vakint::new().unwrap();
    for power in [1, 2] {
        let input = VakintExpression::from_diagram(
            &diagram,
            &family,
            &Atom::one(),
            &DiagramIntegralOptions {
                powers: Some(vec![power.into()]),
                ..Default::default()
            },
        )
        .unwrap();
        let result = engine
            .evaluate(&settings, Atom::from(input).as_view())
            .unwrap();
        let direct = vakint_parse!(&format!("topo(I1L(muvsq,{power}))"))
            .unwrap()
            .replace(vakint_parse!("muvsq").unwrap())
            .with(mass_squared.to_pattern());
        let expected = engine.evaluate(&settings, direct.as_view()).unwrap();
        assert!(!result.is_zero());
        assert!((result - expected).together().expand().is_zero());
    }
}

#[test]
fn symbolic_one_loop_powers_obey_mass_derivative_identity() {
    let settings = symbolic_rustred();
    let engine = Vakint::new().unwrap();
    let result = engine
        .evaluate_integral(
            &settings,
            vakint_parse!(
                "dot(k(1),k(1))*topo(I1L(muvsq,2))-topo(I1L(muvsq,1))-muvsq*topo(I1L(muvsq,2))"
            )
            .unwrap()
            .as_view(),
        )
        .unwrap();
    assert!(result.together().expand().is_zero());
}
