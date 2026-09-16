use super::*;
use crate::Vakint;
use symbolica::id::Replacement;

#[path = "routing_tests/four_loop.rs"]
mod four_loop;

#[test]
fn forced_basis_preserves_mass_and_power_under_a_negative_orientation() {
    Vakint::initialize_vakint_symbols();
    let input = vk_parse!("topo(prop(1,edge(1,1),-k(1),mass_squared,-2))").unwrap();
    let result = Topology::force_an_lmb(input.as_view(), 1).unwrap();
    assert_eq!(
        result,
        vk_parse!("topo(prop(1,edge(1,1),k(1),mass_squared,-2))").unwrap()
    );
}

#[test]
fn retained_basis_coordinates_have_the_canonical_to_parent_direction() {
    Vakint::initialize_vakint_symbols();
    let input = vk_parse!("topo(prop(1,edge(1,1),-k(1),mass_squared,-2))").unwrap();
    let (canonical, coordinates) = Topology::canonicalize_lmb(input.as_view(), 1).unwrap();
    assert_eq!(coordinates.as_ref(), &[vk_parse!("-k(1)").unwrap()]);
    assert_eq!(
        canonical,
        Topology::force_an_lmb(input.as_view(), 1).unwrap()
    );
    let (unchanged, identity) = Topology::canonicalize_lmb(canonical.as_view(), 1).unwrap();
    assert_eq!(unchanged, canonical);
    assert_eq!(identity.as_ref(), &[vk_parse!("k(1)").unwrap()]);

    let defining_input = vk_parse!("topo(prop(1,edge(1,1),-k(1),msq(1),pow(1)))").unwrap();
    let topology = Topology::generate_topology_with_contraction(
        1,
        defining_input.as_view(),
        vk_parse!("I1L(msq(1),pow(1))").unwrap().as_view(),
        Vec::new(),
        EvaluationOrder::empty(),
    )
    .unwrap();
    let (parent, retained) = topology.get_integral().parent_routing.as_ref().unwrap();
    assert_eq!(parent, &defining_input);
    assert_eq!(retained, &coordinates);
}

#[test]
fn forced_basis_rejects_a_parameter_conditional_witness() {
    Vakint::initialize_vakint_symbols();
    let input = vk_parse!("topo(prop(1,edge(1,1),a*k(1),mass_squared,1))").unwrap();
    assert!(Topology::force_an_lmb(input.as_view(), 1).is_err());
}

#[test]
fn forced_basis_rejects_dependent_momentum_routes() {
    Vakint::initialize_vakint_symbols();
    let input = vk_parse!(
        "topo(prop(1,edge(1,2),k(1)+k(2),mass_squared,1)\
         *prop(2,edge(1,2),k(1)+k(2),mass_squared,1)\
         *prop(3,edge(2,1),2*k(1)+2*k(2),mass_squared,1))"
    )
    .unwrap();
    assert!(Topology::force_an_lmb(input.as_view(), 3).is_err());
}

#[test]
fn forced_basis_rejects_a_nonlinear_triangular_route() {
    Vakint::initialize_vakint_symbols();
    let input = vk_parse!(
        "topo(prop(1,edge(1,2),-k(1),mass_squared,1)\
         *prop(2,edge(1,2),k(2)+k(1)^2,mass_squared,1)\
         *prop(3,edge(2,1),-k(1)+k(2)+k(1)^2,mass_squared,1))"
    )
    .unwrap();
    assert!(Topology::force_an_lmb(input.as_view(), 3).is_err());
}

#[test]
fn all_three_loop_matcher_classes_have_exact_parent_basis_transport() {
    Vakint::initialize_vakint_symbols();
    let topologies = Topologies::generate_topologies().unwrap();
    let classes = topologies
        .0
        .iter()
        .filter(|topology| matches!(topology, Topology::ThreeLoop(_)))
        .collect::<Vec<_>>();
    assert_eq!(classes.len(), 5);
    let parent = classes
        .iter()
        .map(|topology| topology.get_integral())
        .find(|integral| {
            Topologies::count_propagators_in_integral(
                integral.canonical_expression.as_ref().unwrap().as_view(),
            ) == 6
        })
        .unwrap();
    let current = (1..=3)
        .map(|index| function!(S.k, Atom::num(index)))
        .collect::<Vec<_>>();
    let parent_coordinates = (1..=3)
        .map(|index| Atom::var(vk_symbol!(format!("parent_loop_{index}"))))
        .collect::<Vec<_>>();
    let to_parent_coordinates = current
        .iter()
        .zip(&parent_coordinates)
        .map(|(source, target)| Replacement::new(source.to_pattern(), target.to_pattern()))
        .collect::<Vec<_>>();
    let from_parent_coordinates = parent_coordinates
        .iter()
        .zip(&current)
        .map(|(source, target)| Replacement::new(source.to_pattern(), target.to_pattern()))
        .collect::<Vec<_>>();
    let powers = [2, -1, 3, 0, 4, -2];
    let mass = vk_parse!("13/10").unwrap();
    let mut nonidentity_maps = 0;
    for topology in classes {
        let integral = topology.get_integral();
        assert_eq!(integral.n_props, 6);
        let mut expression = integral.canonical_expression.as_ref().unwrap().clone();
        expression = expression
            .replace(vk_parse!("msq(1)").unwrap())
            .with(mass.clone());
        for (position, power) in powers.iter().enumerate() {
            expression = expression
                .replace(vk_parse!(&format!("pow({})", position + 1)).unwrap())
                .with(Atom::num(*power));
        }
        let mut source_momenta = Vec::new();
        let mut parent_momenta = Vec::new();
        let mut equations = Vec::new();
        for prop_id in 1..=6 {
            let Some(properties) = get_prop_with_id(expression.as_view(), prop_id) else {
                continue;
            };
            let source = properties[&vk_symbol!("q_")].clone();
            let expected = get_prop_with_id(
                parent.canonical_expression.as_ref().unwrap().as_view(),
                prop_id,
            )
            .unwrap()[&vk_symbol!("q_")]
                .clone();
            equations.push(&source - expected.replace_multiple(&to_parent_coordinates));
            source_momenta.push(source);
            parent_momenta.push(expected);
        }
        let solutions = Atom::solve(&equations).wrt(&current).unwrap();
        assert_eq!(solutions.coverage(), SolveCoverage::Complete);
        assert!(solutions.coverage_guard().is_empty());
        assert_eq!(
            solutions.len(),
            1,
            "{} has no unique parent routing",
            integral.name
        );
        let solution = &solutions[0];
        assert!(solution.is_point());
        let values = solution
            .coordinates()
            .iter()
            .map(|(_, value)| value.clone())
            .collect::<Vec<_>>();
        let (matrix, rhs) =
            Atom::system_to_matrix::<u8, _, _>(&values, &parent_coordinates).unwrap();
        assert!(
            rhs.into_vec()
                .iter()
                .all(|coefficient| coefficient.is_zero())
        );
        let determinant = matrix.det().unwrap().to_expression();
        assert!(determinant == Atom::num(1) || determinant == Atom::num(-1));
        let values = values
            .iter()
            .map(|value| value.replace_multiple(&from_parent_coordinates))
            .collect::<Vec<_>>();
        let (retained_parent, retained_coordinates) = integral.parent_routing.as_ref().unwrap();
        assert_eq!(
            retained_parent,
            parent.canonical_expression.as_ref().unwrap()
        );
        assert_eq!(retained_coordinates.len(), current.len());
        for (retained, independently_solved) in retained_coordinates.iter().zip(&values) {
            assert!((retained - independently_solved).expand().is_zero());
        }
        // Exercise the retained production witness below. The independent
        // solve above is a test oracle, not another runtime routing pass.
        let values = retained_coordinates.as_ref();
        nonidentity_maps += usize::from(values != current.as_slice());
        let witness = current
            .iter()
            .zip(values.iter())
            .map(|(source, target)| Replacement::new(source.to_pattern(), target.to_pattern()))
            .collect::<Vec<_>>();

        // The same simultaneous witness acts on all denominators and a
        // nontrivial scalar numerator; repeated substitutions would be wrong.
        for (source, expected) in source_momenta.iter().zip(&parent_momenta) {
            assert!(
                (source.replace_multiple(&witness) - expected)
                    .expand()
                    .is_zero()
            );
        }
        let numerators = [&source_momenta, &parent_momenta].map(|momenta| {
            (&mass + vk_parse!("2*mursq").unwrap())
                * (function!(S.dot, &momenta[0], &momenta[1])
                    * function!(S.dot, &momenta[1], &momenta[2])
                    + function!(S.dot, &momenta[0], &momenta[0]))
        });
        assert!(
            (numerators[0].replace_multiple(&witness) - &numerators[1])
                .expand()
                .is_zero()
        );

        // Scalar lowering accepts atomic vector arguments, not dot(sum, sum).
        // Exercise the existing component/dot conversion boundary with the
        // same simultaneous witness, including external-vector spectators.
        let component_index = Atom::var(vk_symbol!("routing_component_"));
        let to_components = current
            .iter()
            .enumerate()
            .map(|(axis, source)| {
                Replacement::new(
                    source.to_pattern(),
                    function!(S.k, Atom::num(axis + 1), &component_index).to_pattern(),
                )
                .allow_new_wildcards_on_rhs(true)
            })
            .collect::<Vec<_>>();
        let component_witness = values
            .iter()
            .enumerate()
            .map(|(axis, value)| {
                Replacement::new(
                    function!(S.k, Atom::num(axis + 1), &component_index).to_pattern(),
                    value.replace_multiple(&to_components).to_pattern(),
                )
            })
            .collect::<Vec<_>>();
        let scalar = vk_parse!(
            "(13/10+2*mursq)*(\
             dot(k(1),k(2))*dot(k(2),k(3))\
             +dot(k(1),k(1))*dot(p(1),p(2)))"
        )
        .unwrap();
        let routed_components = Vakint::convert_from_dot_notation(scalar.as_view())
            .replace_multiple(&component_witness)
            .expand();
        let routed_scalar = Vakint::convert_to_dot_notation(routed_components.as_view());

        // Independent component construction uses fixed dummy indices and
        // does not rely on the wildcard component transport being checked.
        let component_values = [101, 102].map(|index| {
            let replacements = current
                .iter()
                .enumerate()
                .map(|(axis, source)| {
                    Replacement::new(
                        source.to_pattern(),
                        function!(S.k, Atom::num(axis + 1), Atom::num(index)).to_pattern(),
                    )
                })
                .collect::<Vec<_>>();
            values
                .iter()
                .map(|value| value.replace_multiple(&replacements))
                .collect::<Vec<_>>()
        });
        let expected_components = (&mass + vk_parse!("2*mursq").unwrap())
            * (&component_values[0][0]
                * &component_values[0][1]
                * &component_values[1][1]
                * &component_values[1][2]
                + component_values[0][0].clone().pow(Atom::num(2))
                    * vk_parse!("p(1,102)*p(2,102)").unwrap());
        let expected_scalar =
            Vakint::convert_to_dot_notation(expected_components.expand().as_view());
        assert!(
            (routed_scalar.clone() - expected_scalar).expand().is_zero(),
            "component transport differs for {}",
            integral.name
        );
        assert!(!routed_scalar.is_zero());
        assert!(routed_scalar.contains_symbol(S.dot));
        for matched in routed_scalar.pattern_match(
            &vk_parse!("dot(left_,right_)").unwrap().to_pattern(),
            None,
            None,
        ) {
            for argument in ["left_", "right_"] {
                let argument = matched.get(&vk_symbol!(argument)).unwrap().clone();
                assert!(
                    current.contains(&argument)
                        || argument == vk_parse!("p(1)").unwrap()
                        || argument == vk_parse!("p(2)").unwrap(),
                    "a non-atomic scalar-product argument escaped: {argument}"
                );
            }
        }
        let transported = expression.replace_multiple(&witness);
        for prop_id in 1..=6 {
            let before = get_prop_with_id(expression.as_view(), prop_id);
            let after = get_prop_with_id(transported.as_view(), prop_id);
            assert_eq!(before.is_some(), after.is_some());
            if let Some(after) = after {
                assert_eq!(after[&vk_symbol!("mUVsq_")], mass);
                assert_eq!(after[&vk_symbol!("pow_")], Atom::num(powers[prop_id - 1]));
            }
        }
    }
    assert!(
        nonidentity_maps > 0,
        "the test must exercise contracted reroutings"
    );
}
