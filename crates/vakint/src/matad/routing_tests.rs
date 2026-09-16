use super::*;

#[test]
fn signed_three_loop_edges_preserve_matad_momentum_identities() {
    Vakint::initialize_vakint_symbols();
    let [a, b, c] = ["routing_a", "routing_b", "routing_c"].map(|name| vk_parse!(name).unwrap());
    // Independent MATAD convention: p2=p1+p3, p4=p1+p6, p5=p6-p3.
    // Fix p4=a, p5=b, p6=c, just as the uncontracted Vakint parent does.
    let matad = [&a - &c, &a - &b, &c - &b, a.clone(), b.clone(), c.clone()];
    let vakint = [a.clone(), b.clone(), c.clone(), &c - &a, &a - &b, &b - &c];
    for ((edge, orientation), expected) in MATAD::edge_momenta("I3L_pinch_1_6")
        .unwrap()
        .iter()
        .zip(&vakint)
    {
        let mapped = Atom::num(*orientation) * &matad[edge - 1];
        assert!((mapped - expected).expand().is_zero());
    }
    assert_eq!(
        MATAD::edge_momenta("I3L")
            .unwrap()
            .iter()
            .map(|(edge, _)| *edge)
            .collect::<Vec<_>>(),
        [4, 5, 6, 1, 2, 3]
    );
}

#[test]
fn one_and_two_loop_matad_routes_are_unchanged() {
    assert_eq!(MATAD::edge_momenta("I1L").unwrap(), &[(1, 1)]);
    assert_eq!(
        MATAD::edge_momenta("I2L").unwrap(),
        &[(2, 1), (3, 1), (1, 1)]
    );
    assert_eq!(
        MATAD::edge_momenta("I2L_pinch_3").unwrap(),
        MATAD::edge_momenta("I2L").unwrap()
    );
    assert!(MATAD::edge_momenta("I4L_H").is_err());
}

#[test]
fn contracted_basketball_basis_retains_the_negative_matad_edge() {
    Vakint::initialize_vakint_symbols();
    let topologies = crate::Topologies::generate_topologies().unwrap();
    let matched = topologies
        .match_topologies_to_user_input(
            vk_parse!("topo(I3L_pinch_1_6(muvsq,0,1,1,1,1,0))")
                .unwrap()
                .as_view(),
            false,
        )
        .unwrap()
        .unwrap();
    let integral = matched.canonical_topology.get_integral();
    let (_, coordinates) = integral.parent_routing.as_ref().unwrap();
    assert_eq!(
        coordinates.as_ref(),
        &[
            vk_parse!("k(2)").unwrap(),
            vk_parse!("k(3)").unwrap(),
            vk_parse!("k(3)-k(1)").unwrap(),
        ]
    );
    let expression = integral.canonical_expression.as_ref().unwrap();
    let edge4 = crate::get_prop_with_id(expression.as_view(), 4).unwrap();
    assert_eq!(edge4[&vk_symbol!("q_")], vk_parse!("k(3)").unwrap());
    let edge5 = crate::get_prop_with_id(expression.as_view(), 5).unwrap();
    assert!(
        (&edge5[&vk_symbol!("q_")] - vk_parse!("-k(1)+k(2)-k(3)").unwrap())
            .expand()
            .is_zero()
    );
    assert_eq!(MATAD::edge_momenta("I3L_pinch_1_6").unwrap()[3], (1, -1));
}
