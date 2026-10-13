//! Generated contraction coverage is distinct from the registered TIDE controls.

use std::collections::{BTreeMap, BTreeSet};

use feynkit_generator::{GenerationOptions, Process};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use symbolica::atom::{Atom, AtomCore};

use crate::graph::Graph;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{Vakint, VakintExpression, diagram_integral::DiagramIntegralOptions, symbols::S};

fn dense_unit_mass_integral(input: &Atom) -> (Atom, usize) {
    let pattern = vk_parse!("prop(id_,edge(a_,b_),q_,m_,p_)").unwrap();
    let mut records = input
        .pattern_match(&pattern.to_pattern(), None, None)
        .collect::<Vec<_>>();
    records.sort_by_key(|record| {
        crate::get_integer_from_atom(record[&vk_symbol!("id_")].as_view()).unwrap()
    });
    let count = records.len();
    let product = records
        .iter()
        .enumerate()
        .map(|(index, record)| {
            S.prop.call((
                index + 1,
                S.edge
                    .call((&record[&vk_symbol!("a_")], &record[&vk_symbol!("b_")])),
                &record[&vk_symbol!("q_")],
                Atom::one(),
                &record[&vk_symbol!("p_")],
            ))
        })
        .product::<Atom>();
    (S.topo.call(product), count)
}

#[test_log::test]
fn generated_contraction_closure_has_native_exact_registered_witnesses() {
    let options = GenerationOptions::default()
        .with_loop_count(5, 5)
        .unwrap()
        .threads(1)
        .max_vertices(8)
        .allow_self_loops(true)
        .allow_zero_flow_edges(false);
    let generated = Process::new(Vec::<&str>::new(), Vec::<&str>::new())
        .generate_diagrams(Model::phi_3_4(), &options)
        .unwrap();
    let vakint = Vakint::new().unwrap();
    vakint
        .topologies
        .prepare_five_loop_catalogue_for_tests()
        .unwrap();
    let mut graphs = BTreeMap::new();
    let mut source_representatives = 0usize;
    let loop_labels = (1..=5).map(|axis| S.k.call(axis)).collect::<Vec<_>>();
    for (graph_index, diagram) in generated.diagrams.iter().enumerate() {
        let family = diagram.propagator_family(&Kinematics::new()).unwrap();
        let basis = diagram.loop_momentum_basis();
        let mut edges = diagram.edges().collect::<Vec<_>>();
        edges.sort_by_key(|(id, _, _)| *id);
        // The raw generated incidence is retained for contraction enumeration.
        // Momentum rows come from FeynKit's stored routing, never a new solver.
        let raw = S.topo.call(
            edges
                .iter()
                .enumerate()
                .map(|(index, (id, endpoints, _))| {
                    let momentum = basis.edge_signatures[id]
                        .loops
                        .apply(&loop_labels)
                        .unwrap()
                        .unwrap_or_default();
                    assert!(
                        !momentum.is_zero(),
                        "zero-flow source edge escaped configured generation filter"
                    );
                    S.prop.call((
                        index + 1,
                        S.edge
                            .call((endpoints.source.unwrap().0, endpoints.target.unwrap().0)),
                        momentum,
                        Atom::one(),
                        1,
                    ))
                })
                .product::<Atom>(),
        );
        let raw_graph = Graph::new_from_atom(raw.as_view(), edges.len()).unwrap();
        assert!(raw_graph.to_symbolica_graph(false).is_connected());
        assert_eq!(raw_graph.to_symbolica_graph(false).num_loops(), 5);
        let contractions = std::iter::once(Vec::new()).chain(raw_graph.find_unique_contractions());
        for contracted in contractions {
            source_representatives += 1;
            let mut powers = vec![1.into(); edges.len()];
            for edge in &contracted {
                powers[*edge - 1] = 0.into();
            }
            // Native ingress applies the actual pinches and combines repeated
            // denominators, using the same original diagram/family routing.
            let input = VakintExpression::from_diagram(
                diagram,
                &family,
                &Atom::one(),
                &DiagramIntegralOptions {
                    powers: Some(powers),
                    ..Default::default()
                },
            )
            .unwrap_or_else(|error| {
                panic!("generated graph {graph_index}, contractions {contracted:?}: {error}")
            });
            let (input, count) = dense_unit_mass_integral(&Atom::from(input));
            let graph = Graph::new_from_atom(input.as_view(), count).unwrap();
            assert_eq!(graph.to_symbolica_graph(false).num_loops(), 5);
            let canonical = graph.to_symbolica_graph(false).canonize().graph;
            graphs
                .entry(canonical)
                .or_insert((graph_index, contracted, input));
            assert!(
                graphs.len() <= 10_000,
                "native contraction census exceeded its graph allowance"
            );
        }
    }
    let mut classes = BTreeMap::<String, usize>::new();
    let mut unknown = Vec::new();
    for (graph_index, contracted, input) in graphs.values() {
        match vakint
            .topologies
            .match_topologies_to_user_input(input.as_view(), false)
            .unwrap()
        {
            Some(matched) => {
                assert_eq!(matched.canonical_topology.get_integral().n_loops, 5);
                *classes
                    .entry(matched.canonical_topology.get_integral().name.clone())
                    .or_default() += 1;
            }
            None => unknown.push((*graph_index, contracted.clone())),
        }
    }
    let registry = super::generate().unwrap();
    let registered = registry
        .iter()
        .map(|topology| topology.get_integral().name.clone())
        .collect::<BTreeSet<_>>();
    let factorized_classes = registry
        .iter()
        .zip(&super::data::CENSUS)
        .filter(|(topology, witness)| {
            witness.factorized && classes.contains_key(&topology.get_integral().name)
        })
        .count();
    println!(
        "generated_graphs={} source_contraction_representatives={source_representatives} native_incidence_classes={} registered_integral_classes={} factorized_classes={factorized_classes} unknown={unknown:?} scaleless_reduction_claim=false filters=phi3+phi4,connected,max_vertices8,self_loops_allowed,zero_flow_edges_forbidden,loop_preserving_contractions_only",
        generated.diagrams.len(),
        graphs.len(),
        classes.len()
    );
    println!("contracted_classes={classes:?}");
    assert!(
        unknown.is_empty(),
        "unregistered native generated contractions"
    );
    assert_eq!(classes.keys().cloned().collect::<BTreeSet<_>>(), registered);
    assert_eq!(factorized_classes, 19);
}
