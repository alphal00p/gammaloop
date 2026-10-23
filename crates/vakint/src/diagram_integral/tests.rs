use super::*;
use feynkit_model::Model;
use symbolica::symbol;

fn scalar_model() -> Model {
    // Native model parameter printers have process-wide identities. Initialize
    // this shared fixture once before parallel diagram tests clone it.
    static MODEL: std::sync::LazyLock<Model> = std::sync::LazyLock::new(|| {
        Model::from_json(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/phi_3_4.json"
        )))
        .expect("valid phi_3_4 fixture")
    });
    MODEL.clone()
}

fn sunset() -> (FeynmanDiagram, IntegralFamily) {
    let diagram = FeynmanDiagram::from_dot(
        scalar_model(),
        "digraph sunset { a -> b [particle=\"phi\"]; a -> b [particle=\"phi\"]; a -> b [particle=\"phi\"]; }",
    ).unwrap();
    let family = diagram.propagator_family(&Kinematics::new()).unwrap();
    (diagram, family)
}

#[test]
fn native_graph_ingress_preserves_high_rank_factored_numerator() {
    let (diagram, family) = sunset();
    let first = &family.denominators()[0] + &family.denominators()[1];
    let second = &family.denominators()[1] - &family.denominators()[2];
    let convert = |numerator: &Atom| {
        VakintExpression::from_diagram(
            &diagram,
            &family,
            numerator,
            &DiagramIntegralOptions::default(),
        )
        .unwrap()
    };
    let first_mapped = convert(&first);
    let second_mapped = convert(&second);
    let complete = convert(&(first.pow(8) * second.pow(9)));
    assert_eq!(complete.0.len(), 1);
    assert_eq!(complete.0[0].integral, first_mapped.0[0].integral);
    assert_eq!(complete.0[0].integral, second_mapped.0[0].integral);
    assert_eq!(
        complete.0[0].numerator,
        first_mapped.0[0].numerator.pow(8) * second_mapped.0[0].numerator.pow(9),
        "native scalar-product conversion must preserve exact numerator factors and powers",
    );
}

#[test]
fn graph_propagator_powers_supply_defaults_and_explicit_powers_override_them() {
    let (diagram, family) = sunset();
    let dotted = diagram
        .map_data(
            |_, vertex| vertex.clone(),
            |id, _, edge| {
                let mut edge = edge.clone();
                if id.0 == 0 {
                    edge.propagator_power = 2;
                }
                edge
            },
        )
        .unwrap();
    let convert = |graph: &FeynmanDiagram, powers: Option<Vec<Integer>>| {
        Atom::from(
            VakintExpression::from_diagram(
                graph,
                &family,
                &Atom::one(),
                &DiagramIntegralOptions {
                    powers,
                    ..Default::default()
                },
            )
            .unwrap(),
        )
    };
    assert_eq!(
        convert(&dotted, None),
        convert(&diagram, Some(vec![2.into(), 1.into(), 1.into()]))
    );
    assert_eq!(
        convert(&dotted, Some(vec![1.into(); 3])),
        convert(&diagram, None)
    );
}

#[test]
fn native_graph_ingress_preserves_high_rank_factored_numerator() {
    let (diagram, family) = sunset();
    let first = &family.denominators()[0] + &family.denominators()[1];
    let second = &family.denominators()[1] - &family.denominators()[2];
    let convert = |numerator: &Atom| {
        VakintExpression::from_diagram(
            &diagram,
            &family,
            numerator,
            &DiagramIntegralOptions::default(),
        )
        .unwrap()
    };
    let first_mapped = convert(&first);
    let second_mapped = convert(&second);
    let complete = convert(&(first.pow(8) * second.pow(9)));
    assert_eq!(complete.0.len(), 1);
    assert_eq!(complete.0[0].integral, first_mapped.0[0].integral);
    assert_eq!(complete.0[0].integral, second_mapped.0[0].integral);
    assert_eq!(
        complete.0[0].numerator,
        first_mapped.0[0].numerator.pow(8) * second_mapped.0[0].numerator.pow(9),
        "native scalar-product conversion must preserve exact numerator factors and powers",
    );
}

#[test]
fn repeated_native_propagators_coalesce_exact_powers_before_matching() {
    let diagram = FeynmanDiagram::from_dot(scalar_model(),
        "digraph repeated { a -> b [particle=\"phi\"]; a -> b [particle=\"phi\"]; a -> a [particle=\"phi\"]; b -> b [particle=\"phi\"]; }").unwrap();
    let family = diagram.propagator_family(&Kinematics::new()).unwrap();
    assert_eq!(diagram.loop_count(), 3);
    let expression = Atom::from(
        VakintExpression::from_diagram(
            &diagram,
            &family,
            &Atom::one(),
            &DiagramIntegralOptions::default(),
        )
        .unwrap(),
    );
    let terms = VakintExpression::split_integrals(expression.as_view()).unwrap();
    let vakint = crate::Vakint::new().unwrap();
    let matched = vakint
        .topologies
        .match_topologies_to_user_input(terms[0].integral.as_view(), false)
        .unwrap()
        .unwrap();
    assert_eq!(matched.canonical_topology.get_integral().n_loops, 3);
    let pattern =
        crate::utils::vakint_macros::vk_parse!("prop(id_,edge(nl_,nr_),q_,m_,power_)").unwrap();
    let mut powers = terms[0]
        .integral
        .pattern_match(&pattern.to_pattern(), None, None)
        .map(|matched| {
            crate::get_integer_from_atom(
                matched[&crate::utils::vakint_macros::vk_symbol!("power_")].as_view(),
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    powers.sort();
    assert_eq!(powers, [1, 1, 2]);
}

#[test]
fn massive_zero_flow_bridge_is_a_scalar_factor_without_losing_loops() {
    let diagram = FeynmanDiagram::from_dot(scalar_model(),
        "digraph two_tadpoles { a -> a [particle=\"phi\"]; a -> b [particle=\"phi\"]; b -> b [particle=\"phi\"]; }").unwrap();
    let original = diagram.propagator_family(&Kinematics::new()).unwrap();
    assert_eq!(diagram.loop_count(), 2);
    let mut edges = diagram.edges().collect::<Vec<_>>();
    edges.sort_by_key(|(id, _, _)| *id);
    let bridge = edges
        .iter()
        .position(|(id, _, _)| {
            diagram.loop_momentum_basis().edge_signatures[id]
                .loops
                .integer_coefficients()
                .iter()
                .all(|value| *value == 0)
        })
        .unwrap();
    let symbols = original.denominators()[bridge].get_all_symbols(true);
    assert_eq!(
        symbols.len(),
        1,
        "bridge denominator has exactly the native scalar mass symbol"
    );
    let mass = Atom::var(*symbols.iter().next().unwrap());
    let vakint = crate::Vakint::new().unwrap();
    for mass_value in [None, Some(3)] {
        let replacements = mass_value
            .map(|value| Replacement::new(mass.clone(), Atom::num(value)))
            .into_iter()
            .collect::<Vec<_>>();
        let family = IntegralFamily::new(
            original.loop_momenta().to_vec(),
            vec![],
            original
                .denominators()
                .iter()
                .map(|denominator| denominator.replace_multiple(&replacements))
                .collect(),
            original.kinematics(),
        )
        .unwrap();
        let options = DiagramIntegralOptions {
            parameter_substitutions: mass_value
                .map(|value| (mass.clone(), Atom::num(value)))
                .into_iter()
                .collect(),
            ..Default::default()
        };
        let convert = |bridge_power: i64| {
            let mut powers = vec![Integer::from(1); edges.len()];
            powers[bridge] = bridge_power.into();
            Atom::from(
                VakintExpression::from_diagram(
                    &diagram,
                    &family,
                    &Atom::one(),
                    &DiagramIntegralOptions {
                        powers: Some(powers),
                        ..options.clone()
                    },
                )
                .unwrap(),
            )
        };
        let baseline = convert(0);
        for power in [1, 2, 0, -1] {
            let actual = convert(power);
            let expected = family.denominators()[bridge].clone().pow(-power) * &baseline;
            assert!((actual.clone() - expected).expand().together().is_zero());
            let terms = VakintExpression::split_integrals(actual.as_view()).unwrap();
            assert_eq!(terms.len(), 1);
            let matched = vakint
                .topologies
                .match_topologies_to_user_input(terms[0].integral.as_view(), false)
                .unwrap()
                .unwrap();
            assert_eq!(matched.canonical_topology.get_integral().n_loops, 2);
        }
    }
}

#[test]
fn five_loop_self_edge_pinch_never_downgrades_the_integration_dimension() {
    let diagram = FeynmanDiagram::from_dot(
        scalar_model(),
        "digraph five_tadpoles {\
        a -> a [particle=\"phi\"]; a -> b [particle=\"phi\"];\
        b -> b [particle=\"phi\"]; b -> c [particle=\"phi\"];\
        c -> c [particle=\"phi\"]; c -> d [particle=\"phi\"];\
        d -> d [particle=\"phi\"]; d -> e [particle=\"phi\"];\
        e -> e [particle=\"phi\"]; }",
    )
    .unwrap();
    assert_eq!(diagram.loop_count(), 5);
    let family = diagram.propagator_family(&Kinematics::new()).unwrap();
    let mut edges = diagram.edges().collect::<Vec<_>>();
    edges.sort_by_key(|(id, _, _)| *id);
    let self_edge = edges
        .iter()
        .position(|(_, endpoints, _)| endpoints.source == endpoints.target)
        .unwrap();
    let mut powers = vec![1.into(); family.denominators().len()];
    powers[self_edge] = 0.into();
    let error = VakintExpression::from_diagram(
        &diagram,
        &family,
        &Atom::one(),
        &DiagramIntegralOptions {
            powers: Some(powers),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("pinch removes an independent vacuum loop")
    );
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

#[test]
fn graph_ingress_explicit_slots_transport_pinches_and_negative_powers() {
    let (diagram, family) = sunset();
    let reordered = IntegralFamily::new(
        family.loop_momenta().to_vec(),
        vec![],
        [2, 0, 1]
            .map(|slot| family.denominators()[slot].clone())
            .to_vec(),
        family.kinematics(),
    )
    .unwrap();
    let options = DiagramIntegralOptions {
        powers: Some(vec![0.into(), 2.into(), 3.into()]),
        physical_slots: Some(vec![1, 2, 0]),
        ..Default::default()
    };
    let pinched =
        VakintExpression::from_diagram(&diagram, &reordered, &family.denominators()[2], &options)
            .unwrap();
    let negative = VakintExpression::from_diagram(
        &diagram,
        &family,
        &Atom::one(),
        &DiagramIntegralOptions {
            powers: Some(vec![2.into(), 3.into(), (-1).into()]),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(Atom::from(pinched), Atom::from(negative));

    let original = VakintExpression::from_diagram(
        &diagram,
        &family,
        &Atom::one(),
        &DiagramIntegralOptions::default(),
    )
    .unwrap();
    let permuted = VakintExpression::from_diagram(
        &diagram,
        &reordered,
        &Atom::one(),
        &DiagramIntegralOptions {
            physical_slots: Some(vec![1, 2, 0]),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(Atom::from(original), Atom::from(permuted));
}

#[test]
fn graph_ingress_refuses_ambiguous_slots_and_loop_removing_pinches() {
    let (diagram, family) = sunset();
    for slots in [vec![0, 0, 2], vec![0, 1, 3], vec![0, 1]] {
        let error = VakintExpression::from_diagram(
            &diagram,
            &family,
            &Atom::one(),
            &DiagramIntegralOptions {
                physical_slots: Some(slots),
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(error.to_string().contains("physical slots"));
    }
    let error = VakintExpression::from_diagram(
        &diagram,
        &family,
        &Atom::one(),
        &DiagramIntegralOptions {
            powers: Some(vec![0.into(), 0.into(), 1.into()]),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(error.to_string().contains("independent vacuum loop"));
}
