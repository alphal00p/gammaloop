use super::*;
use crate::{Vakint, VakintSettings, VakintTerm};

fn root(vakint: &Vakint) -> &Topology {
    vakint
        .topologies
        .0
        .iter()
        .find(|topology| {
            matches!(topology, Topology::TwoLoop(_))
                && topology.get_integral().graph.edges.len() == 3
        })
        .unwrap()
}
fn powers(candidate: &ReplacementRules) -> Vec<i64> {
    (1..=candidate.canonical_topology.get_integral().n_props)
        .map(|slot| {
            get_integer_from_atom(
                candidate.canonical_expression_substitutions[&S.pow.call(slot)].as_view(),
            )
            .unwrap()
        })
        .collect()
}
fn mapped_numerator(candidate: &ReplacementRules, numerator: Atom) -> Atom {
    let mut term = VakintTerm {
        numerator,
        integral: Atom::one(),
        vectors: vec![],
    };
    term.apply_numerator_replacement_rules(candidate, &VakintSettings::default())
        .unwrap();
    term.numerator.expand()
}
fn base() -> Atom {
    vk_parse!("topo(prop(1,edge(1,2),k(1),m2,1)*prop(2,edge(1,2),k(2),m2,2)*prop(3,edge(2,1),k(1)+k(2),m2,3))").unwrap()
}

#[test]
fn scope_rejection_of_first_witness_selects_verified_alternative() {
    let vakint = Vakint::new().unwrap();
    let topology = root(&vakint);
    let expression = base();
    let source = MomentumInput::from_expression(expression.as_view(), 2)
        .unwrap()
        .unwrap();
    let original = source.match_topologies([topology]).unwrap().unwrap();
    let mut desired = powers(&original);
    desired.reverse();
    assert_ne!(powers(&original), desired);
    let mut first = None;
    let accepted = source
        .match_topologies_with_acceptance([topology], Default::default(), |candidate| {
            let candidate = powers(candidate);
            first.get_or_insert(candidate.clone());
            Ok((candidate == desired).then_some(candidate))
        })
        .unwrap()
        .unwrap();
    assert_eq!(first.unwrap(), powers(&original));
    assert_eq!(powers(&accepted), desired);
    let refused = source
        .match_topologies_with_acceptance([topology], Default::default(), |_| Ok(None::<()>))
        .unwrap()
        .unwrap();
    assert_eq!(powers(&refused), powers(&original));
}

#[test]
fn accepted_witness_output_is_unimodular_orientation_and_label_invariant() {
    let vakint = Vakint::new().unwrap();
    let topology = root(&vakint);
    let fixtures=[
        (base(),vk_parse!("k(1)^2").unwrap()),
        (vk_parse!("topo(prop(1,edge(1,2),k(1)+k(2),m2,1)*prop(2,edge(1,2),-k(2),m2,2)*prop(3,edge(2,1),k(1),m2,3))").unwrap(),vk_parse!("(k(1)+k(2))^2").unwrap()),
        (vk_parse!("topo(prop(30,edge(2,7),-k(1),m2,1)*prop(10,edge(7,2),k(2),m2,2)*prop(20,edge(2,7),k(1)+k(2),m2,3))").unwrap(),vk_parse!("k(1)^2").unwrap()),
    ];
    let mut expected = None;
    for (input, numerator) in fixtures {
        let source = MomentumInput::from_expression(input.as_view(), 2)
            .unwrap()
            .unwrap();
        let accepted = source
            .match_topologies_with_acceptance([topology], Default::default(), |candidate| {
                Ok(Some((
                    powers(candidate),
                    mapped_numerator(candidate, numerator.clone()),
                )))
            })
            .unwrap()
            .unwrap();
        let result = (powers(&accepted), mapped_numerator(&accepted, numerator));
        if let Some(expected) = &expected {
            assert_eq!(&result, expected);
        } else {
            expected = Some(result);
        }
    }
}

#[test]
fn alternative_witness_limits_and_acceptance_errors_are_typed_and_propagated() {
    let vakint = Vakint::new().unwrap();
    let topology = root(&vakint);
    let expression = base();
    let source = MomentumInput::from_expression(expression.as_view(), 2)
        .unwrap()
        .unwrap();
    for limits in [
        MomentumMatchLimits {
            max_graph_automorphisms: 0,
            ..Default::default()
        },
        MomentumMatchLimits {
            max_orientation_candidates: 0,
            ..Default::default()
        },
        MomentumMatchLimits {
            max_verified_witnesses: 0,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            source.match_topologies_with_acceptance(
                [topology],
                limits,
                |_| -> Result<Option<()>, VakintError> { panic!("exhaustion precedes acceptance") }
            ),
            Err(VakintError::MomentumMatchResourceLimit { .. })
        ));
    }
    let error = source
        .match_topologies_with_acceptance(
            [topology],
            Default::default(),
            |_| -> Result<Option<()>, VakintError> {
                Err(VakintError::InvalidNumerator("predicate-sentinel".into()))
            },
        )
        .unwrap_err();
    assert!(
        matches!(error,VakintError::InvalidNumerator(message) if message=="predicate-sentinel")
    );
}

#[test]
fn native_edge_marking_includes_parallel_and_self_loop_permutations() {
    for self_loops in [false, true] {
        let mut graph = NativeGraph::new();
        let a = graph.add_node(Empty);
        let b = if self_loops { a } else { graph.add_node(Empty) };
        for _ in 0..3 {
            graph.add_edge(a, b, false, Atom::one()).unwrap();
        }
        let marked = edge_orbits::EdgeGraph::new(&graph).unwrap();
        let mut maps = BTreeSet::new();
        let mut count = 0;
        edge_orbits::visit(&marked, &marked, 6, &mut count, |map| {
            maps.insert(map.to_vec());
            Ok(())
        })
        .unwrap();
        assert_eq!(maps.len(), 6);
        assert_eq!(count, 6);
        assert!(matches!(
            edge_orbits::visit(&marked, &marked, 5, &mut 0, |_| Ok(())),
            Err(VakintError::MomentumMatchResourceLimit { .. })
        ));
    }
}
