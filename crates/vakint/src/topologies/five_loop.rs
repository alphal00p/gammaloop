//! Five-loop topology data uses the same contraction and momentum matcher as
//! the lower-loop registry. A census representative retains its defining
//! parent's simultaneous routing and explicit common-basis physical slots.

#[path = "five_loop_data.rs"]
mod data;

#[cfg(all(test, feature = "feynkit-ingress"))]
#[path = "five_loop_census_tests.rs"]
mod census_tests;

use std::{collections::BTreeSet, sync::LazyLock};
use symbolica::atom::{Atom, AtomView, Symbol};

use super::Topology;
use crate::{EvaluationOrder, VakintError, symbols::S, utils::vakint_macros::vk_symbol};

pub(crate) struct ParentDescriptor {
    pub momenta: Vec<Atom>,
    pub physical_slots: &'static [usize],
}

pub(crate) fn common_momentum_rows() -> &'static [[i64; 5]; 15] {
    &data::MOMENTA
}

pub(crate) fn descriptors() -> Vec<ParentDescriptor> {
    data::PARENTS
        .iter()
        .map(|parent| ParentDescriptor {
            momenta: parent
                .physical_slots
                .iter()
                .map(|slot| {
                    data::MOMENTA[*slot]
                        .iter()
                        .enumerate()
                        .map(|(axis, coefficient)| Atom::num(*coefficient) * S.k.call(axis + 1))
                        .sum()
                })
                .collect(),
            physical_slots: &parent.physical_slots,
        })
        .collect()
}

/// Symbol identities come from the same compiled census and contraction
/// convention as generation, without constructing any Integral or graph.
pub(super) fn is_known_short_form(input: AtomView) -> bool {
    static HEADS: LazyLock<BTreeSet<Symbol>> = LazyLock::new(|| {
        data::CENSUS
            .iter()
            .map(|witness| {
                let parent = &data::PARENTS[witness.parent];
                let contraction = witness
                    .deleted_slots
                    .iter()
                    .map(|slot| {
                        parent
                            .physical_slots
                            .iter()
                            .position(|physical| physical == slot)
                            .unwrap()
                            + 1
                    })
                    .collect::<Vec<_>>();
                Topology::contracted_short_symbol(
                    vk_symbol!(format!("I5L_{}", parent.sector)),
                    &contraction,
                )
            })
            .collect()
    });
    let AtomView::Fun(topo) = input else {
        return false;
    };
    if topo.get_symbol() != S.topo || topo.get_nargs() != 1 {
        return false;
    }
    matches!(topo.iter().next(), Some(AtomView::Fun(head)) if HEADS.contains(&head.get_symbol()))
}

pub(super) fn generate() -> Result<Vec<Topology>, VakintError> {
    let descriptors = descriptors();
    let parents = data::PARENTS
        .iter()
        .zip(&descriptors)
        .map(|(parent, descriptor)| {
            let mut expression = Atom::one();
            let mut arguments = vec![vk_symbol!("msq").call(1)];
            for (physical, slot) in parent.physical_slots.iter().enumerate() {
                let signed_slot = (*slot + 1) as i64;
                let source = parent
                    .vertices
                    .iter()
                    .position(|v| v.contains(&signed_slot))
                    .unwrap();
                let target = parent
                    .vertices
                    .iter()
                    .position(|v| v.contains(&-signed_slot))
                    .unwrap();
                expression *= S.prop.call((
                    physical + 1,
                    S.edge.call((source + 1, target + 1)),
                    &descriptor.momenta[physical],
                    vk_symbol!("msq").call(1),
                    S.pow.call(physical + 1),
                ));
                arguments.push(S.pow.call(physical + 1));
            }
            let short = vk_symbol!(format!("I5L_{}", parent.sector)).call(&arguments);
            (S.topo.call(expression), short)
        })
        .collect::<Vec<_>>();
    data::CENSUS
        .iter()
        .map(|witness| {
            let descriptor = &descriptors[witness.parent];
            let (expression, short) = &parents[witness.parent];
            let contractions = witness
                .deleted_slots
                .iter()
                .map(|slot| {
                    descriptor
                        .physical_slots
                        .iter()
                        .position(|p| p == slot)
                        .unwrap()
                        + 1
                })
                .collect();
            Topology::generate_topology_with_contraction(
                descriptor.physical_slots.len(),
                expression.as_view(),
                short.as_view(),
                contractions,
                EvaluationOrder(vec![
                    crate::EvaluationMethod::RustRed(Default::default()),
                    crate::EvaluationMethod::PySecDec(Default::default()),
                ]),
            )
            .map_err(|error| {
                VakintError::InvalidIntegralFormat(format!(
                    "five-loop census sector {}: {error}",
                    witness.sector
                ))
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::atom::AtomCore;

    #[test]
    fn registry_keeps_default_numerical_fallback_eligible() {
        let vakint = crate::Vakint::new().unwrap();
        vakint
            .topologies
            .prepare_five_loop_catalogue_for_tests()
            .unwrap();
        let settings = crate::VakintSettings::default();
        for topology in vakint
            .topologies
            .iter()
            .filter(|topology| matches!(topology, Topology::FiveLoop(_)))
        {
            assert!(settings.evaluation_order.0.iter().any(|method| matches!(
                method,
                crate::EvaluationMethod::PySecDec(_)
            )
                && method.supports(&settings, topology)));
        }
    }

    #[test]
    fn census_routes_preserve_all_physical_momenta_and_factorized_classes() {
        for parent in &data::PARENTS {
            for vertex in &parent.vertices {
                for axis in 0..5 {
                    let balance: i64 = vertex
                        .iter()
                        .map(|slot| {
                            slot.signum() * data::MOMENTA[slot.unsigned_abs() as usize - 1][axis]
                        })
                        .sum();
                    assert_eq!(balance, 0, "parent {}", parent.sector);
                }
            }
        }
        assert_eq!(data::CENSUS.iter().filter(|w| w.factorized).count(), 19);
        for witness in &data::CENSUS {
            let parent = &data::PARENTS[witness.parent];
            assert_eq!(witness.source_slots.len(), witness.target_slots.len());
            assert_eq!(witness.source_slots.len() + witness.deleted_slots.len(), 12);
            for ((source, target), sign) in witness
                .source_slots
                .iter()
                .zip(witness.target_slots)
                .zip(witness.signs)
            {
                assert!(parent.physical_slots.contains(target));
                assert!(!witness.deleted_slots.contains(target));
                for axis in 0..5 {
                    let transformed: i64 = (0..5)
                        .map(|i| data::MOMENTA[*source][i] * witness.routing[i][axis])
                        .sum();
                    assert_eq!(
                        transformed,
                        sign * data::MOMENTA[*target][axis],
                        "sector {}",
                        witness.sector
                    );
                }
            }
        }
    }

    #[test]
    fn census_registration_retains_all_parent_routings() {
        crate::Vakint::initialize_vakint_symbols();
        let descriptors = descriptors();
        let topologies = generate().unwrap();
        assert_eq!(topologies.len(), 67);
        for (topology, witness) in topologies.iter().zip(&data::CENSUS) {
            assert!(
                matches!(topology, Topology::FiveLoop(_)),
                "sector {}",
                witness.sector
            );
            topology
                .get_integral()
                .validate_parent_routing(&descriptors[witness.parent].momenta)
                .unwrap();
            let graph = &topology.get_integral().graph;
            assert_eq!(
                graph.get_one_lmb().unwrap().len(),
                5,
                "sector {}",
                witness.sector
            );
            for node in graph.nodes.values() {
                let balance: Atom = node
                    .edges
                    .iter()
                    .map(|(edge, direction)| {
                        let momentum = graph.edges[edge].momentum.clone();
                        if direction.is_incoming() {
                            momentum
                        } else {
                            -momentum
                        }
                    })
                    .sum();
                assert!(
                    balance.expand().is_zero(),
                    "sector {} vertex {} violates momentum conservation",
                    witness.sector,
                    node.id
                );
            }
        }
    }

    #[test]
    fn five_loop_basis_changes_retain_simultaneous_exact_witnesses() {
        use crate::get_prop_with_id;
        use symbolica::id::Replacement;
        crate::Vakint::initialize_vakint_symbols();
        let topologies = generate().unwrap();
        let routing = [
            S.k.call(2) + S.k.call(3),
            -S.k.call(1),
            S.k.call(3),
            S.k.call(5),
            S.k.call(4),
        ];
        let replacements = routing
            .iter()
            .enumerate()
            .map(|(axis, target)| Replacement::new(S.k.call(axis + 1), target.clone()))
            .collect::<Vec<_>>();
        for (topology, witness) in topologies
            .iter()
            .zip(&data::CENSUS)
            .filter(|(_, witness)| witness.deleted_slots.is_empty())
        {
            let original = topology
                .get_integral()
                .parent_routing
                .as_ref()
                .unwrap()
                .0
                .clone();
            let transformed = original.replace_multiple(&replacements).expand();
            let (canonical, coordinates) =
                Topology::canonicalize_lmb(transformed.as_view(), 12).unwrap();
            let retained = coordinates
                .iter()
                .enumerate()
                .map(|(axis, target)| Replacement::new(S.k.call(axis + 1), target.clone()))
                .collect::<Vec<_>>();
            for physical in 1..=12 {
                let actual = get_prop_with_id(canonical.as_view(), physical).unwrap()
                    [&vk_symbol!("q_")]
                    .clone();
                let expected = get_prop_with_id(transformed.as_view(), physical).unwrap()
                    [&vk_symbol!("q_")]
                    .clone();
                assert!(
                    (actual.replace_multiple(&retained) - expected)
                        .expand()
                        .is_zero(),
                    "sector {} physical slot {physical}",
                    witness.sector
                );
            }
        }
    }

    #[test]
    fn factorized_graph_attachments_share_an_exact_momentum_match() {
        use crate::utils::vakint_macros::vk_parse;
        use crate::{Vakint, VakintSettings, VakintTerm};
        let vakint = Vakint::new().unwrap();
        vakint
            .topologies
            .prepare_five_loop_catalogue_for_tests()
            .unwrap();
        for last_attachment in [1, 2] {
            let expression = vk_parse!(format!(
                "topo(\
                prop(1,edge(1,2),k(1),m2,2)*\
                prop(2,edge(1,2),k(2),m2,3)*\
                prop(3,edge(2,1),k(1)+k(2),m2,4)*\
                prop(4,edge(1,1),k(3),m2,5)*\
                prop(5,edge(1,1),k(4),m2,6)*\
                prop(6,edge({last_attachment},{last_attachment}),k(5),m2,7))"
            ))
            .unwrap();
            let mut matched = vakint
                .topologies
                .match_topologies_to_user_input(expression.as_view(), false)
                .unwrap()
                .expect("nonisomorphic vacuum-factor attachment must retain its integral class");
            matched.apply_replacement_rules().unwrap();
            assert_eq!(matched.canonical_topology.get_integral().n_loops, 5);
            let mut term = VakintTerm {
                integral: expression,
                numerator: vk_parse!("dot(k(1),k(5))^2").unwrap(),
                vectors: Vec::new(),
            };
            term.apply_numerator_replacement_rules(
                &matched,
                &VakintSettings {
                    verify_numerator_identification: true,
                    ..Default::default()
                },
            )
            .unwrap();
            assert!(!term.numerator.is_zero());
        }
    }

    #[test]
    fn disconnected_vacuum_products_match_all_independent_loops() {
        use crate::Vakint;
        use crate::utils::vakint_macros::vk_parse;
        let vakint = Vakint::new().unwrap();
        vakint
            .topologies
            .prepare_five_loop_catalogue_for_tests()
            .unwrap();
        let mut expressions = Vec::new();
        for nodes in [[1, 2, 3, 4, 5], [1, 1, 1, 2, 2]] {
            let mut expression = Atom::one();
            for (axis, node) in nodes.into_iter().enumerate() {
                expression *= S.prop.call((
                    axis + 1,
                    S.edge.call((node, node)),
                    S.k.call(axis + 1),
                    Atom::var(vk_symbol!("m2")),
                    1,
                ));
            }
            expressions.push(S.topo.call(expression));
        }
        expressions.push(
            vk_parse!(
                "topo(\
            prop(1,edge(1,2),k(1),m2,1)*prop(2,edge(1,2),k(2),m2,1)*\
            prop(3,edge(1,2),k(3),m2,1)*prop(4,edge(2,1),k(1)+k(2)+k(3),m2,1)*\
            prop(5,edge(3,4),k(4),m2,1)*prop(6,edge(3,4),k(5),m2,1)*\
            prop(7,edge(4,3),k(4)+k(5),m2,1))"
            )
            .unwrap(),
        );
        for expression in expressions {
            let matched = vakint
                .topologies
                .match_topologies_to_user_input(expression.as_view(), false)
                .unwrap()
                .expect("a disconnected product retains every independent integration loop");
            assert_eq!(matched.canonical_topology.get_integral().n_loops, 5);
        }
    }

    #[test]
    fn native_five_loop_matching_refuses_nonconservation_and_preserves_unknown_fallback() {
        use crate::utils::vakint_macros::vk_parse;
        let vakint = crate::Vakint::new().unwrap();
        vakint
            .topologies
            .prepare_five_loop_catalogue_for_tests()
            .unwrap();
        let malformed = vk_parse!(
            "topo(prop(1,edge(1,2),k(1),m2,1)*\
            prop(2,edge(1,2),k(2),m2,1)*\
            prop(3,edge(2,1),k(1)+2*k(2),m2,1)*\
            prop(4,edge(1,1),k(3),m2,1)*\
            prop(5,edge(1,1),k(4),m2,1)*\
            prop(6,edge(1,1),k(5),m2,1))"
        )
        .unwrap();
        let error = vakint
            .topologies
            .match_topologies_to_user_input(malformed.as_view(), true)
            .unwrap_err();
        assert!(error.to_string().contains("momentum conservation"));

        // This graph has five independent loops but its physical basis carries
        // determinant 2. A polynomial permutation cannot prove unit measure.
        let nonunit = vk_parse!(
            "topo(prop(1,edge(1,1),2*k(1),m2,1)*\
            prop(2,edge(1,1),k(2),m2,1)*\
            prop(3,edge(1,1),k(3),m2,1)*\
            prop(4,edge(1,1),k(4),m2,1)*\
            prop(5,edge(1,1),k(5),m2,1))"
        )
        .unwrap();
        assert!(
            vakint
                .topologies
                .match_topologies_to_user_input(nonunit.as_view(), false)
                .unwrap()
                .is_none()
        );
        let fallback = vakint
            .topologies
            .match_topologies_to_user_input(nonunit.as_view(), true)
            .unwrap()
            .expect("unknown-integral fallback remains available after failed exact proof");
        assert!(matches!(fallback.canonical_topology, Topology::Unknown(_)));
    }

    #[cfg(feature = "feynkit-ingress")]
    #[test_log::test]
    fn native_dot_roundtrips_all_census_pinches() {
        native_dot_roundtrips(data::CENSUS.iter());
    }

    #[cfg(feature = "feynkit-ingress")]
    #[test_log::test]
    fn native_dot_four_plus_one_factor_routing() {
        native_dot_roundtrips(
            data::CENSUS
                .iter()
                .filter(|witness| witness.sector == 32652),
        );
    }

    #[cfg(feature = "feynkit-ingress")]
    #[test_log::test]
    fn generated_scalar_vacuum_graphs_match_registered_integral_classes() {
        use crate::{Vakint, VakintExpression, diagram_integral::DiagramIntegralOptions};
        use feynkit_generator::{GenerationOptions, Process};
        use feynkit_graph::FeynmanDiagram;
        use feynkit_kinematics::Kinematics;
        use feynkit_model::Model;

        let options = GenerationOptions::default()
            .with_loop_count(5, 5)
            .unwrap()
            .threads(1)
            .allow_self_loops(true)
            .max_vertices(8);
        let generated = Process::new(Vec::<&str>::new(), Vec::<&str>::new())
            .generate_diagrams(Model::phi_3_4(), &options)
            .unwrap();
        assert!(!generated.diagrams.is_empty());
        let vakint = Vakint::new().unwrap();
        vakint
            .topologies
            .prepare_five_loop_catalogue_for_tests()
            .unwrap();
        let mut classes = std::collections::BTreeMap::<String, usize>::new();
        for (index, graph) in generated.diagrams.iter().enumerate() {
            let started = std::time::Instant::now();
            log::debug!(target: "vakint::five_loop_match", "generated_graph={index} phase=begin");
            let dot = graph.to_dot().unwrap();
            let decoded = FeynmanDiagram::from_dot(Model::phi_3_4(), &dot).unwrap();
            let family = decoded.propagator_family(&Kinematics::new()).unwrap();
            let input = Atom::from(
                VakintExpression::from_diagram(
                    &decoded,
                    &family,
                    &Atom::one(),
                    &DiagramIntegralOptions::default(),
                )
                .unwrap(),
            );
            let terms = VakintExpression::split_integrals(input.as_view()).unwrap();
            assert_eq!(terms.len(), 1);
            let matched = vakint
                .topologies
                .match_topologies_to_user_input(terms[0].integral.as_view(), false)
                .unwrap_or_else(|error| panic!("generated graph {index}: {error}\n{dot}"))
                .unwrap_or_else(|| {
                    panic!("generated graph {index} has no registered class\n{dot}")
                });
            assert_eq!(matched.canonical_topology.get_integral().n_loops, 5);
            *classes
                .entry(matched.canonical_topology.get_integral().name.clone())
                .or_default() += 1;
            log::debug!(target: "vakint::five_loop_match", "generated_graph={index} elapsed_seconds={:.6}", started.elapsed().as_secs_f64());
        }
        println!(
            "generated_graphs={} matched_classes={} registry_classes={} classes={classes:?}",
            generated.diagrams.len(),
            classes.len(),
            data::CENSUS.len(),
        );
        // FeynGen's connected scalar graph census and the integral-class census
        // need not have identical cardinalities: bridge denominators can coincide,
        // and disconnected factors may attach to different vertices.
    }

    #[cfg(feature = "feynkit-ingress")]
    fn native_dot_roundtrips<'a>(witnesses: impl IntoIterator<Item = &'a data::Witness>) {
        use crate::{Vakint, VakintExpression, diagram_integral::DiagramIntegralOptions};
        use feynkit_graph::FeynmanDiagram;
        use feynkit_kinematics::Kinematics;
        use feynkit_model::Model;

        let vakint = Vakint::new().unwrap();
        vakint
            .topologies
            .prepare_five_loop_catalogue_for_tests()
            .unwrap();
        let parents = data::PARENTS
            .iter()
            .map(|parent| {
                let edges = parent
                    .physical_slots
                    .iter()
                    .map(|slot| {
                        let signed = (*slot + 1) as i64;
                        let source = parent
                            .vertices
                            .iter()
                            .position(|vertex| vertex.contains(&signed))
                            .unwrap();
                        let target = parent
                            .vertices
                            .iter()
                            .position(|vertex| vertex.contains(&-signed))
                            .unwrap();
                        format!("v{source} -> v{target} [particle=\"phi\"];")
                    })
                    .collect::<String>();
                FeynmanDiagram::from_dot(Model::phi_3_4(), &format!("digraph {{ {edges} }}"))
                    .unwrap()
            })
            .collect::<Vec<_>>();

        for witness in witnesses {
            let started = std::time::Instant::now();
            log::debug!(target: "vakint::five_loop_match", "census_sector={} phase=native_dot_roundtrip", witness.sector);
            let diagram = &parents[witness.parent];
            let roundtrip =
                FeynmanDiagram::from_dot(Model::phi_3_4(), &diagram.to_dot().unwrap()).unwrap();
            let mut expressions = Vec::new();
            for graph in [diagram, &roundtrip] {
                // DOT parsers may reorder edges. Resolve physical slots through
                // the retained vertex names, never through insertion positions.
                let parent = &data::PARENTS[witness.parent];
                let mut edges = graph.edges().collect::<Vec<_>>();
                edges.sort_by_key(|(id, _, _)| *id);
                let powers = edges
                    .iter()
                    .map(|(_, endpoints, _)| {
                        let source = &graph.vertex(endpoints.source).unwrap().name;
                        let target = &graph.vertex(endpoints.target).unwrap().name;
                        let slot = parent
                            .physical_slots
                            .iter()
                            .find(|slot| {
                                let signed = (**slot + 1) as i64;
                                let left = format!(
                                    "v{}",
                                    parent
                                        .vertices
                                        .iter()
                                        .position(|vertex| vertex.contains(&signed))
                                        .unwrap()
                                );
                                let right = format!(
                                    "v{}",
                                    parent
                                        .vertices
                                        .iter()
                                        .position(|vertex| vertex.contains(&-signed))
                                        .unwrap()
                                );
                                (source == &left && target == &right)
                                    || (source == &right && target == &left)
                            })
                            .expect("native edge has an explicit parent physical slot");
                        i64::from(!witness.deleted_slots.contains(slot)).into()
                    })
                    .collect();
                let options = DiagramIntegralOptions {
                    powers: Some(powers),
                    ..Default::default()
                };
                let family = graph.propagator_family(&Kinematics::new()).unwrap();
                log::debug!(target: "vakint::five_loop_match", "census_sector={} phase=diagram_family elapsed_seconds={:.6}", witness.sector, started.elapsed().as_secs_f64());
                let expression =
                    VakintExpression::from_diagram(graph, &family, &Atom::one(), &options)
                        .unwrap_or_else(|error| panic!("sector {}: {error}", witness.sector));
                log::debug!(target: "vakint::five_loop_match", "census_sector={} phase=diagram_ingress elapsed_seconds={:.6}", witness.sector, started.elapsed().as_secs_f64());
                let expression = Atom::from(expression);
                let terms = VakintExpression::split_integrals(expression.as_view()).unwrap();
                assert_eq!(terms.len(), 1);
                let matched = vakint
                    .topologies
                    .match_topologies_to_user_input(terms[0].integral.as_view(), false)
                    .unwrap()
                    .unwrap_or_else(|| panic!("missing census sector {}", witness.sector));
                assert_eq!(matched.canonical_topology.get_integral().n_loops, 5);
                expressions.push(expression);
            }
            assert_eq!(
                expressions[0], expressions[1],
                "DOT roundtrip sector {}",
                witness.sector
            );
            log::debug!(target: "vakint::five_loop_match", "census_sector={} elapsed_seconds={:.6}", witness.sector, started.elapsed().as_secs_f64());
        }
    }
}
