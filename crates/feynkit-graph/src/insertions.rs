//! Proper two-point regions, with attachment identities in the parent diagram.

use std::collections::BTreeSet;

use linnet::half_edge::{
    involution::{EdgeIndex, Flow},
    subgraph::{ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};

use crate::{DiagramEndpoint, DiagramError, DiagramHalfEdge, EdgeId, FeynmanDiagram};

/// An amputated, connected, one-particle-irreducible propagator insertion.
///
/// The support and both boundary endpoints refer to the original diagram.
/// Boundary endpoints lie inside the region; their opposite endpoints belong
/// to the surrounding graph. Their order follows [`DiagramHalfEdge`], not a
/// particle-flow convention. Species, flow and vertex slots remain available
/// through the parent, including off-diagonal two-point functions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PropagatorInsertion {
    /// Both half-edges of each internal propagator in this region.
    pub support: SuBitGraph,
    /// The two unselected half-edges incident to the region's vertices.
    pub boundary: [DiagramHalfEdge; 2],
    pub loops: usize,
}

impl FeynmanDiagram {
    /// Find all proper, loop-bearing two-point regions without changing the graph.
    ///
    /// Original external attachments and physical cut carriers cannot be
    /// swallowed or become insertion boundaries. Dummy attachments do not carry
    /// momentum and are ignored. Topology-only threshold candidates impose no
    /// restriction. Boundary species need not agree: eligibility for a particular
    /// propagator evaluator is separate from this structural classification.
    ///
    /// Nested regions are all returned. A caller extracting several regions at
    /// once must select vertex-disjoint supports; edge disjointness alone is not
    /// sufficient. Serial insertions are returned separately because their union
    /// contains a bridge. Each support can be passed to [`Self::excise`].
    ///
    /// A proper two-point region is separated by removing its two attachment
    /// edges. Enumerating these pairs avoids a powerset of internal edges.
    pub fn propagator_insertions(&self) -> Result<Vec<PropagatorInsertion>, DiagramError> {
        self.validate()?;
        let momentum = self.momentum_subgraph();
        let mut protected: SuBitGraph = self.graph.empty_subgraph();
        for (pair, _, edge) in self.graph.iter_edges() {
            if edge.data.external.is_some() && !edge.data.is_dummy {
                protected.add(pair);
            }
        }
        for cut in &self.cuts {
            for half_edge in &cut.cut {
                protected.add(self.graph[&EdgeIndex(half_edge.edge.0)].1);
            }
        }
        let protected_vertices = self
            .graph
            .nodes(&protected)
            .into_iter()
            .collect::<BTreeSet<_>>();
        let internal = self.internal_subgraph().subtract(&protected);
        let edges = self
            .graph
            .iter_edges_of(&internal)
            .map(|(pair, _, _)| pair)
            .collect::<Vec<_>>();
        let mut seen = BTreeSet::new();
        let mut insertions = Vec::new();
        for (index, first) in edges.iter().enumerate() {
            for second in &edges[index + 1..] {
                let mut remainder = internal.clone();
                remainder.sub(*first);
                remainder.sub(*second);
                for support in self.graph.connected_components(&remainder) {
                    if !seen.insert(support.clone()) {
                        continue;
                    }
                    let crown = self.graph.full_crown(&support).intersection(&momentum);
                    if crown.n_included() != 2
                        // Opening one omitted edge inside the region is not
                        // an insertion attached to the surrounding graph.
                        || self.graph.iter_edges_of(&crown).count() != 2
                        || self
                            .graph
                            .iter_nodes_of(&support)
                            .any(|(vertex, _, _)| protected_vertices.contains(&vertex))
                    {
                        continue;
                    }
                    let loops = self.graph.cyclotomatic_number(&support);
                    if loops == 0 || !self.graph.bridges_of(&support).is_empty() {
                        continue;
                    }
                    let mut boundary = crown
                        .included_iter()
                        .map(|hedge| DiagramHalfEdge {
                            edge: EdgeId(self.graph[&hedge].0),
                            endpoint: match self.graph.flow(hedge) {
                                Flow::Source => DiagramEndpoint::Source,
                                Flow::Sink => DiagramEndpoint::Target,
                            },
                        })
                        .collect::<Vec<_>>();
                    boundary.sort();
                    insertions.push(PropagatorInsertion {
                        support,
                        boundary: [boundary[0], boundary[1]],
                        loops,
                    });
                }
            }
        }
        insertions.sort_by(|left, right| left.support.cmp(&right.support));
        Ok(insertions)
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use feynkit_model::Model;
    use idenso::tensor::SymbolicTensor;
    use linnet::half_edge::subgraph::Inclusion;
    use symbolica::atom::Atom;

    use super::*;
    use crate::{DiagramEdge, DiagramVertex, ExternalLeg, ExternalState, VertexId};

    fn model() -> Arc<Model> {
        Arc::new(Model::from_json(r#"{
            "name":"propagator_insertions","restriction":null,"orders":[],
            "parameters":[{"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null}],
            "particles":[
                {"pdg_code":25,"name":"phi","antiname":"phi","spin":1,"color":1,"mass":"ZERO","width":"ZERO","texname":"phi","antitexname":"phi","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                {"pdg_code":1,"name":"f","antiname":"f~","spin":2,"color":1,"mass":"ZERO","width":"ZERO","texname":"f","antitexname":"fbar","charge":0.0,"ghost_number":0,"lepton_number":1,"y_charge":0},
                {"pdg_code":-1,"name":"f~","antiname":"f","spin":2,"color":1,"mass":"ZERO","width":"ZERO","texname":"fbar","antitexname":"f","charge":0.0,"ghost_number":0,"lepton_number":-1,"y_charge":0}
            ],
            "propagators":[],"lorentz_structures":[],"couplings":[],"vertex_rules":[]
        }"#).unwrap())
    }

    fn diagram(edges: &[(usize, usize)], external_vertices: &[usize]) -> FeynmanDiagram {
        let model = model();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "insertion_test");
        let nodes = edges
            .iter()
            .flat_map(|&(a, b)| [a, b])
            .chain(external_vertices.iter().copied())
            .max()
            .unwrap()
            + 1;
        for vertex in 0..nodes {
            builder.add_vertex(DiagramVertex {
                name: format!("v{vertex}"),
                interaction: None,
                numerator: SymbolicTensor::infer(Atom::one()).unwrap(),
            });
        }
        for &(source, target) in edges {
            builder
                .add_edge(
                    VertexId(source),
                    VertexId(target),
                    DiagramEdge::new(particle, false),
                )
                .unwrap();
        }
        for (index, &vertex) in external_vertices.iter().enumerate() {
            let mut edge = DiagramEdge::new(particle, false);
            edge.external = Some(ExternalLeg {
                name: format!("p{index}"),
                index,
                state: ExternalState::Outgoing,
                connection: index,
            });
            builder
                .add_edge(Some(VertexId(vertex)), None, edge)
                .unwrap();
        }
        builder.build().unwrap()
    }

    fn edge_ids(diagram: &FeynmanDiagram, insertion: &PropagatorInsertion) -> Vec<usize> {
        diagram
            .graph
            .iter_edges_of(&insertion.support)
            .map(|(_, edge, _)| edge.0)
            .collect()
    }

    fn single_bubble() -> FeynmanDiagram {
        diagram(&[(0, 1), (0, 2), (2, 3), (2, 3), (3, 1)], &[0, 1])
    }

    #[test]
    fn single_insertion_preserves_parent_ports_and_can_be_excised() {
        let parent = single_bubble();
        let before = parent.to_json().unwrap();
        let insertions = parent.propagator_insertions().unwrap();
        assert_eq!(insertions.len(), 1);
        let insertion = &insertions[0];
        assert_eq!(edge_ids(&parent, insertion), [2, 3]);
        assert_eq!(insertion.loops, 1);
        assert_eq!(
            insertion.boundary,
            [
                DiagramHalfEdge {
                    edge: EdgeId(1),
                    endpoint: DiagramEndpoint::Target
                },
                DiagramHalfEdge {
                    edge: EdgeId(4),
                    endpoint: DiagramEndpoint::Source
                },
            ]
        );
        for port in insertion.boundary {
            let hedge = parent.half_edge_id(port).unwrap();
            assert!(!insertion.support.includes(&hedge));
            assert!(
                parent
                    .graph
                    .nodes(&insertion.support)
                    .contains(&parent.graph.node_id(hedge))
            );
            assert!(
                !parent
                    .graph
                    .nodes(&insertion.support)
                    .contains(&parent.graph.node_id(parent.graph.inv(hedge)))
            );
        }
        let child = parent.excise(&insertion.support, &[]).unwrap();
        assert_eq!(child.loop_count(), 1);
        assert_eq!(
            child
                .edges()
                .filter(|(_, _, edge)| edge.external.is_some())
                .count(),
            2
        );
        assert!(child.propagator_insertions().unwrap().is_empty());
        assert_eq!(parent.to_json().unwrap(), before);
    }

    #[test]
    fn nested_insertions_are_both_reported() {
        let parent = diagram(
            &[
                (0, 1),
                (0, 2),
                (2, 3),
                (2, 4),
                (4, 5),
                (4, 5),
                (5, 3),
                (3, 1),
            ],
            &[0, 1],
        );
        let insertions = parent.propagator_insertions().unwrap();
        assert_eq!(insertions.len(), 2);
        let inner = insertions
            .iter()
            .find(|insertion| insertion.loops == 1)
            .unwrap();
        let outer = insertions
            .iter()
            .find(|insertion| insertion.loops == 2)
            .unwrap();
        assert_eq!(edge_ids(&parent, inner), [4, 5]);
        assert_eq!(edge_ids(&parent, outer), [2, 3, 4, 5, 6]);
        assert!(outer.support.includes(&inner.support));
        assert_eq!(
            parent
                .excise(&outer.support, &[])
                .unwrap()
                .propagator_insertions()
                .unwrap()
                .len(),
            1
        );
    }

    #[test]
    fn disjoint_insertions_remain_independent() {
        let parent = diagram(
            &[
                (0, 2),
                (2, 3),
                (2, 3),
                (3, 1),
                (0, 4),
                (4, 5),
                (4, 5),
                (5, 1),
            ],
            &[0, 1],
        );
        let insertions = parent.propagator_insertions().unwrap();
        assert_eq!(insertions.len(), 2);
        assert!(insertions.iter().all(|insertion| insertion.loops == 1));
        let left = parent.graph.nodes(&insertions[0].support);
        let right = parent.graph.nodes(&insertions[1].support);
        assert!(left.iter().all(|vertex| !right.contains(vertex)));
    }

    #[test]
    fn serial_insertions_do_not_include_the_reducible_union() {
        let parent = diagram(
            &[
                (0, 1),
                (0, 2),
                (2, 3),
                (2, 3),
                (3, 4),
                (4, 5),
                (4, 5),
                (5, 1),
            ],
            &[0, 1],
        );
        let insertions = parent.propagator_insertions().unwrap();
        assert_eq!(insertions.len(), 2);
        assert!(insertions.iter().all(|insertion| insertion.loops == 1));
        assert!(
            insertions
                .iter()
                .all(|insertion| !edge_ids(&parent, insertion).contains(&4))
        );
    }

    #[test]
    fn external_regions_trees_and_three_point_regions_are_excluded() {
        for parent in [
            diagram(&[(0, 1), (0, 1)], &[0, 1]),
            diagram(&[(0, 2), (2, 3), (3, 1)], &[0, 1]),
            diagram(
                &[(0, 2), (1, 3), (5, 4), (2, 3), (3, 4), (4, 2)],
                &[0, 1, 5],
            ),
        ] {
            assert!(parent.propagator_insertions().unwrap().is_empty());
        }
    }

    #[test]
    fn opening_an_internal_edge_of_a_vacuum_component_is_not_an_insertion() {
        let parent = diagram(&[(0, 1), (0, 1), (0, 1), (2, 3)], &[]);
        assert!(parent.propagator_insertions().unwrap().is_empty());
    }

    #[test]
    fn two_ports_at_one_vertex_can_bound_a_tadpole_insertion() {
        let parent = diagram(&[(0, 1), (0, 2), (2, 2), (2, 1)], &[0, 1]);
        let insertions = parent.propagator_insertions().unwrap();
        assert_eq!(insertions.len(), 1);
        assert_eq!(edge_ids(&parent, &insertions[0]), [2]);
        assert_eq!(insertions[0].loops, 1);
        assert_eq!(
            parent
                .excise(&insertions[0].support, &[])
                .unwrap()
                .loop_count(),
            1
        );
    }

    #[test]
    fn physical_cut_carriers_cannot_be_swallowed() {
        let parent = diagram(&[(0, 1), (0, 2), (2, 3), (2, 3), (3, 1)], &[]);
        let parent = parent
            .map_data(
                |_, vertex| vertex.clone(),
                |id, _, edge| {
                    let mut edge = edge.clone();
                    if id == EdgeId(0) {
                        edge.external = Some(ExternalLeg {
                            name: "sewn".into(),
                            index: 0,
                            state: ExternalState::Incoming,
                            connection: 0,
                        });
                    }
                    edge
                },
            )
            .unwrap()
            .with_canonical_loop_momentum_basis()
            .unwrap();
        let (left, right) = parent.half_edges().partition(|port| {
            let vertex = parent.graph.node_id(parent.half_edge_id(*port).unwrap());
            vertex.0 == 0 || vertex.0 == 2
        });
        let parent = parent.with_cut_partitions(vec![(left, right)]).unwrap();
        assert_eq!(parent.cuts()[0].cut.len(), 2);
        assert!(parent.propagator_insertions().unwrap().is_empty());
    }

    #[test]
    fn mixed_species_and_reversed_flow_preserve_attachment_slots() {
        let parent = single_bubble();
        let fermion = parent.model().particle_id("f").unwrap();
        let parent = parent
            .map_data(
                |_, vertex| vertex.clone(),
                |id, _, edge| {
                    let mut edge = edge.clone();
                    if id == EdgeId(1) {
                        edge.particle = fermion;
                        edge.directed = true;
                    }
                    edge
                },
            )
            .unwrap();
        let old_slot = parent.graph[EdgeIndex(1)].target_slot();
        let parent = parent.reverse_edge(EdgeId(1)).unwrap();
        let insertions = parent.propagator_insertions().unwrap();
        assert_eq!(insertions.len(), 1);
        assert_eq!(
            insertions[0].boundary[0],
            DiagramHalfEdge {
                edge: EdgeId(1),
                endpoint: DiagramEndpoint::Source
            }
        );
        assert_eq!(parent.graph[EdgeIndex(1)].source_slot(), old_slot);
        let child = parent.excise(&insertions[0].support, &[]).unwrap();
        let species = child
            .edges()
            .filter(|(_, _, edge)| edge.external.is_some())
            .map(|(_, _, edge)| edge.particle)
            .collect::<BTreeSet<_>>();
        assert_eq!(species.len(), 2);
    }

    #[test]
    fn dummy_attachments_do_not_change_physical_boundaries() {
        let parent = diagram(&[(0, 1), (0, 2), (2, 3), (2, 3), (3, 1), (2, 4)], &[0, 1]);
        let parent = parent
            .map_data(
                |_, vertex| vertex.clone(),
                |id, _, edge| {
                    let mut edge = edge.clone();
                    edge.is_dummy = id == EdgeId(5);
                    edge
                },
            )
            .unwrap()
            .with_canonical_loop_momentum_basis()
            .unwrap();
        let insertions = parent.propagator_insertions().unwrap();
        assert_eq!(insertions.len(), 1);
        assert_eq!(edge_ids(&parent, &insertions[0]), [2, 3]);
        assert!(
            insertions[0]
                .boundary
                .iter()
                .all(|port| port.edge != EdgeId(5))
        );
    }
}
