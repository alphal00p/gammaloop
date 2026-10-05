use super::*;
use feynkit_model::Model;
use symbolica::symbol;

fn sunset() -> (FeynmanDiagram, IntegralFamily) {
    let diagram = FeynmanDiagram::from_dot(
        Model::phi_3_4(),
        "digraph sunset { a -> b [particle=\"phi\"]; a -> b [particle=\"phi\"]; a -> b [particle=\"phi\"]; }",
    ).unwrap();
    let family = diagram.propagator_family(&Kinematics::new()).unwrap();
    (diagram, family)
}

#[test]
fn graph_ingress_preserves_scalar_product_and_auxiliary_power() {
    let (diagram, family) = sunset();
    let product = family
        .kinematics()
        .scalar_product(&family.loop_momenta()[0], &family.loop_momenta()[1])
        .unwrap();
    let extended = IntegralFamily::new(
        family.loop_momenta().to_vec(),
        vec![],
        family
            .denominators()
            .iter()
            .cloned()
            .chain([product.clone()])
            .collect(),
        family.kinematics(),
    )
    .unwrap();
    let ordinary = VakintExpression::from_diagram(
        &diagram,
        &family,
        &product,
        &DiagramIntegralOptions::default(),
    )
    .unwrap();
    let auxiliary = VakintExpression::from_diagram(
        &diagram,
        &extended,
        &Atom::one(),
        &DiagramIntegralOptions {
            powers: Some(vec![1.into(), 1.into(), 1.into(), (-1).into()]),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(Atom::from(ordinary), Atom::from(auxiliary));
}

#[test]
fn graph_ingress_substitutions_are_simultaneous() {
    let (diagram, family) = sunset();
    let x = Atom::var(symbol!("vakint_ingress_test::x"));
    let y = Atom::var(symbol!("vakint_ingress_test::y"));
    let exchanged = VakintExpression::from_diagram(
        &diagram,
        &family,
        &(&x + Atom::num(2) * &y),
        &DiagramIntegralOptions {
            parameter_substitutions: vec![(x.clone(), y.clone()), (y.clone(), x.clone())],
            ..Default::default()
        },
    )
    .unwrap();
    let expected = VakintExpression::from_diagram(
        &diagram,
        &family,
        &(&y + Atom::num(2) * &x),
        &DiagramIntegralOptions::default(),
    )
    .unwrap();
    assert_eq!(Atom::from(exchanged), Atom::from(expected));
    assert!(
        VakintExpression::from_diagram(
            &diagram,
            &family,
            &Atom::one(),
            &DiagramIntegralOptions {
                parameter_substitutions: vec![(x.clone(), y.clone()), (x, y)],
                ..Default::default()
            }
        )
        .unwrap_err()
        .to_string()
        .contains("sources must be distinct")
    );
}

#[test]
fn graph_ingress_refuses_family_sign_and_spectator_ambiguity() {
    let (diagram, family) = sunset();
    let mut denominators = family.denominators().to_vec();
    denominators[0] = -denominators[0].clone();
    let wrong = IntegralFamily::new(
        family.loop_momenta().to_vec(),
        vec![],
        denominators,
        family.kinematics(),
    )
    .unwrap();
    assert!(
        VakintExpression::from_diagram(
            &diagram,
            &wrong,
            &Atom::one(),
            &DiagramIntegralOptions::default()
        )
        .unwrap_err()
        .to_string()
        .contains("order/sign/masses")
    );
    assert!(
        VakintExpression::from_diagram(
            &diagram,
            &family,
            &Atom::one(),
            &DiagramIntegralOptions {
                external_momenta: vec![family.loop_momenta()[0].clone()],
                ..Default::default()
            }
        )
        .unwrap_err()
        .to_string()
        .contains("distinct")
    );
    assert!(
        VakintExpression::from_diagram(
            &diagram,
            &family,
            &Atom::one(),
            &DiagramIntegralOptions {
                external_momenta: vec![
                    family.loop_momenta()[0].clone() + &family.loop_momenta()[1]
                ],
                ..Default::default()
            }
        )
        .unwrap_err()
        .to_string()
        .contains("momentum name")
    );
}

#[test]
fn graph_ingress_refuses_undeclared_spectator_products() {
    let (diagram, family) = sunset();
    let spectator = Atom::var(symbol!("vakint_ingress_test::spectator"));
    let free = Kinematics::new()
        .with_momenta(
            family
                .loop_momenta()
                .iter()
                .cloned()
                .chain([spectator.clone()]),
        )
        .unwrap();
    let numerator = free
        .scalar_product(&family.loop_momenta()[0], &spectator)
        .unwrap();
    assert!(
        VakintExpression::from_diagram(
            &diagram,
            &family,
            &numerator,
            &DiagramIntegralOptions::default()
        )
        .unwrap_err()
        .to_string()
        .contains("unconverted")
    );
    assert!(
        VakintExpression::from_diagram(
            &diagram,
            &family,
            &numerator,
            &DiagramIntegralOptions {
                external_momenta: vec![spectator],
                ..Default::default()
            }
        )
        .is_ok()
    );
}
