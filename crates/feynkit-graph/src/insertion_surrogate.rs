//! Scalar representatives for already integrated two-point insertions.

use std::collections::{BTreeMap, BTreeSet};

use feynkit_model::ParticleId;
use idenso::tensor::SymbolicTensor;
use linnet::half_edge::{
    NodeIndex,
    involution::{Flow, HedgePair, Orientation},
    subgraph::Inclusion,
};
use spenso::structure::partial::PartialStructure;
use symbolica::atom::Atom;

use crate::{
    DiagramEdge, DiagramEndpoint, DiagramError, DiagramHalfEdge, DiagramRelabeling, DiagramVertex,
    FeynmanDiagram, PropagatorInsertion, VertexId, VertexSlot,
};

/// An L-loop insertion represented by L scalar one-loop bubbles in series.
///
/// The scalar chain has unit local numerators; its particle is only massless
/// denominator metadata. No Feynman rule or closed-fermion sign is inferred.
/// `with_coefficient` attaches the original insertion's open tensor divided by
/// the scalar chain's bare integral, retaining its full regulator dependence.
/// This is a bare-integral identity, not permission to substitute inside a
/// local counterterm operation or to discard the original subdivergences.
#[derive(Clone, Debug)]
pub struct ScalarBubbleSurrogate {
    pub diagram: FeynmanDiagram,
    pub relabeling: DiagramRelabeling,
    pub coefficient_vertex: VertexId,
}

impl ScalarBubbleSurrogate {
    /// Attach an insertion coefficient expressed in the original diagram's owners.
    /// The map remains available to restore the final integrated open tensor.
    pub fn with_coefficient(
        &self,
        coefficient: &SymbolicTensor<PartialStructure>,
    ) -> Result<FeynmanDiagram, DiagramError> {
        let coefficient = self.relabeling.translate_tensor(coefficient, None)?;
        self.diagram
            .clone()
            .with_numerator_fragments([(self.coefficient_vertex, coefficient)], [])
    }
}

impl FeynmanDiagram {
    /// Replace a proper two-point region by an equally many-loop scalar chain.
    ///
    /// Both original boundary propagators, their directions and outside vertex
    /// slots survive. The removed region has no symbolic owner image; the two
    /// inner boundary half-edges do. Overall factors, numerator prefactors and
    /// projectors are transported without inventing symmetry or species factors.
    /// `carrier` may have any spin or color, but must be massless in the model.
    /// Its artificial edges have explicit unit numerators and power one, and
    /// artificial vertices have no interaction rule. Do not regenerate Feynman
    /// rules on the resulting auxiliary diagram.
    ///
    /// Physical cuts are unsupported. Topology-only threshold candidates are
    /// cleared because partitions through the removed interior have no image.
    pub fn scalar_bubble_surrogate(
        &self,
        insertion: &PropagatorInsertion,
        carrier: ParticleId,
    ) -> Result<ScalarBubbleSurrogate, DiagramError> {
        let error = |message: &str| DiagramError::Invariant {
            operation: "constructing a scalar insertion representative",
            message: message.into(),
        };
        self.validate()?;
        if !self.cuts.is_empty() {
            return Err(error(
                "physical cuts are unsupported for scalar insertion representatives",
            ));
        }
        let particle = self.model.particle_by_id(carrier)?;
        if !particle.is_massless(&self.model) {
            return Err(error("the scalar denominator carrier must be massless"));
        }
        // Unit numerators do not remove the carrier's model flow metadata.
        let directed = !self.model.particle_is_self_conjugate(carrier);
        let orientation = if !directed {
            Orientation::Undirected
        } else if particle.is_antiparticle() {
            Orientation::Reversed
        } else {
            Orientation::Default
        };
        if !self.propagator_insertions()?.contains(insertion) {
            return Err(error(
                "the selected region is not a proper protected two-point insertion",
            ));
        }
        let removed = self
            .graph
            .nodes(&insertion.support)
            .into_iter()
            .collect::<BTreeSet<_>>();
        let one = SymbolicTensor::infer(Atom::one())?;
        let mut builder =
            Self::builder(self.model_arc(), format!("{}_scalar_insertion", self.name))
                .symmetry_factor(self.symmetry_factor);
        let mut vertices = BTreeMap::new();
        for (id, vertex) in self.vertices() {
            if !removed.contains(&NodeIndex(id.0)) {
                vertices.insert(id.0, builder.add_vertex(vertex.clone()).0);
            }
        }
        let chain = (0..=insertion.loops)
            .map(|i| {
                builder.add_vertex(DiagramVertex {
                    name: format!("scalar_insertion_{i}"),
                    interaction: None,
                    numerator: one.clone(),
                })
            })
            .collect::<Vec<_>>();
        let mut chain_slots = vec![0; chain.len()];
        let mut edges = BTreeMap::new();
        let mut half_edges = BTreeMap::new();
        let mut orientations = Vec::new();
        for (pair, id, data) in self.graph.iter_edges() {
            if insertion.support.includes(&pair) {
                continue;
            }
            let old_ends = match pair {
                HedgePair::Paired { source, sink } => vec![source, sink],
                HedgePair::Unpaired { hedge, .. } => vec![hedge],
                HedgePair::Split { .. } => unreachable!("a full graph has no split edges"),
            };
            if old_ends
                .iter()
                .all(|h| removed.contains(&self.graph.node_id(*h)))
            {
                // Dummy attachments carry no momentum and may lie entirely
                // inside the selected region without belonging to its support.
                if data.data.is_dummy {
                    continue;
                }
                return Err(error(
                    "an unselected physical edge lies inside the insertion",
                ));
            }
            let mut source = None;
            let mut target = None;
            let mut source_slot = data.data.source_slot();
            let mut target_slot = data.data.target_slot();
            for &hedge in &old_ends {
                let endpoint = match self.graph.flow(hedge) {
                    Flow::Source => DiagramEndpoint::Source,
                    Flow::Sink => DiagramEndpoint::Target,
                };
                let old_vertex = self.graph.node_id(hedge);
                let mut slot = None;
                let vertex = if let Some(&vertex) = vertices.get(&old_vertex.0) {
                    VertexId(vertex)
                } else {
                    let boundary = DiagramHalfEdge {
                        edge: crate::EdgeId(id.0),
                        endpoint,
                    };
                    let index = if boundary == insertion.boundary[0] {
                        0
                    } else if boundary == insertion.boundary[1] {
                        insertion.loops
                    } else if data.data.is_dummy {
                        0
                    } else {
                        return Err(error("an insertion has more than two momentum attachments"));
                    };
                    slot = Some(VertexSlot(chain_slots[index]));
                    chain_slots[index] += 1;
                    chain[index]
                };
                match endpoint {
                    DiagramEndpoint::Source => {
                        source = Some(vertex);
                        if let Some(slot) = slot {
                            source_slot = slot;
                        }
                    }
                    DiagramEndpoint::Target => {
                        target = Some(vertex);
                        if let Some(slot) = slot {
                            target_slot = slot;
                        }
                    }
                }
            }
            let new = builder.add_edge_with_slots(
                source,
                target,
                data.data.clone(),
                source_slot,
                target_slot,
            )?;
            edges.insert(id.0, vec![new.0]);
            orientations.push(data.orientation);
            for hedge in old_ends {
                half_edges.insert(
                    hedge,
                    DiagramHalfEdge {
                        edge: new,
                        endpoint: match self.graph.flow(hedge) {
                            Flow::Source => DiagramEndpoint::Source,
                            Flow::Sink => DiagramEndpoint::Target,
                        },
                    },
                );
            }
        }
        for index in 0..insertion.loops {
            for _ in 0..2 {
                let source_slot = VertexSlot(chain_slots[index]);
                let target_slot = VertexSlot(chain_slots[index + 1]);
                chain_slots[index] += 1;
                chain_slots[index + 1] += 1;
                builder.add_edge_with_slots(
                    chain[index],
                    chain[index + 1],
                    DiagramEdge::new(carrier, directed),
                    source_slot,
                    target_slot,
                )?;
                orientations.push(orientation);
            }
        }
        builder.edge_orientations = Some(orientations);
        let mut diagram = builder.build()?;
        let hedges = half_edges
            .into_iter()
            .map(|(old, new)| {
                (
                    old.0,
                    diagram.half_edge_id(new).expect("retained endpoint").0,
                )
            })
            .collect();
        let relabeling = DiagramRelabeling::new(
            self,
            &diagram,
            vertices.clone(),
            edges,
            hedges,
            BTreeMap::new(),
        );
        let originals = vertices
            .into_iter()
            .map(|(old, new)| (new, VertexId(old)))
            .collect::<BTreeMap<_, _>>();
        diagram.graph = diagram.graph.map_result(
            |_, id, mut vertex| {
                if let Some(original) = originals.get(&id.0) {
                    vertex.numerator =
                        relabeling.translate_tensor(&vertex.numerator, Some(*original))?;
                }
                Ok::<_, DiagramError>(vertex)
            },
            |_, _, _, _, mut data| {
                data.data.numerator = relabeling.translate_tensor(&data.data.numerator, None)?;
                Ok(data)
            },
            |_, data| Ok(data),
        )?;
        diagram.overall_factor = relabeling.translate_atom(&self.overall_factor, None)?;
        diagram.numerator_prefactor =
            relabeling.translate_tensor(&self.numerator_prefactor, None)?;
        diagram.projector = relabeling.translate_tensor(&self.projector, None)?;
        diagram.numerator = Self::fragment_numerator(&diagram.graph)?;
        if diagram.loop_count() != self.loop_count() {
            return Err(error("the scalar representative changed the loop count"));
        }
        diagram.validate()?;
        Ok(ScalarBubbleSurrogate {
            diagram,
            relabeling,
            coefficient_vertex: chain[0],
        })
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use feynkit_model::Model;
    use linnet::half_edge::subgraph::{ModifySubSet, SuBitGraph, SubSetLike};
    use symbolica::function;

    use super::*;
    use crate::{EdgeId, ExternalLeg, ExternalState, NumeratorInterface, symbols};

    fn model() -> Arc<Model> {
        Arc::new(Model::from_json(r#"{
            "name":"scalar_insertion_surrogate","restriction":null,"orders":[],
            "parameters":[
                {"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null},
                {"name":"M","lhablock":"MASS","lhacode":[25],"nature":"external","parameter_type":"real","value":[1.0,0.0],"expression":null}
            ],
            "particles":[
                {"pdg_code":21,"name":"g","antiname":"g","spin":3,"color":8,"mass":"ZERO","width":"ZERO","texname":"g","antitexname":"g","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                {"pdg_code":25,"name":"h","antiname":"h","spin":1,"color":1,"mass":"M","width":"ZERO","texname":"h","antitexname":"h","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                {"pdg_code":9000001,"name":"c","antiname":"c~","spin":-1,"color":8,"mass":"ZERO","width":"ZERO","texname":"c","antitexname":"cbar","charge":0.0,"ghost_number":1,"lepton_number":0,"y_charge":0},
                {"pdg_code":-9000001,"name":"c~","antiname":"c","spin":-1,"color":8,"mass":"ZERO","width":"ZERO","texname":"cbar","antitexname":"c","charge":0.0,"ghost_number":-1,"lepton_number":0,"y_charge":0},
                {"pdg_code":1,"name":"q","antiname":"q~","spin":2,"color":3,"mass":"ZERO","width":"ZERO","texname":"q","antitexname":"qbar","charge":0.0,"ghost_number":0,"lepton_number":1,"y_charge":0},
                {"pdg_code":-1,"name":"q~","antiname":"q","spin":2,"color":-3,"mass":"ZERO","width":"ZERO","texname":"qbar","antitexname":"q","charge":0.0,"ghost_number":0,"lepton_number":-1,"y_charge":0}
            ],
            "propagators":[],"lorentz_structures":[],"couplings":[],"vertex_rules":[]
        }"#).unwrap())
    }

    fn diagram(edges: &[(usize, usize)]) -> FeynmanDiagram {
        let model = model();
        let particle = model.particle_id("g").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "insertion_parent")
            .overall_factor(Atom::num(7))
            .numerator_prefactor(SymbolicTensor::infer(Atom::num(13)).unwrap())
            .projector(SymbolicTensor::infer(Atom::num(17)).unwrap())
            .symmetry_factor(2);
        for vertex in 0..=edges.iter().flat_map(|&(a, b)| [a, b]).max().unwrap() {
            builder.add_vertex(DiagramVertex {
                name: format!("v{vertex}"),
                interaction: None,
                numerator: SymbolicTensor::infer(Atom::one()).unwrap(),
            });
        }
        for &(a, b) in edges {
            builder
                .add_edge(VertexId(a), VertexId(b), DiagramEdge::new(particle, false))
                .unwrap();
        }
        for vertex in 0..2 {
            let mut edge = DiagramEdge::new(particle, false);
            edge.external = Some(ExternalLeg {
                name: format!("p{vertex}"),
                index: vertex,
                state: ExternalState::Outgoing,
                connection: vertex,
            });
            builder.add_edge(VertexId(vertex), None, edge).unwrap();
        }
        builder.build().unwrap()
    }

    fn banana(loops: usize) -> FeynmanDiagram {
        let mut edges = vec![(0, 1), (0, 2)];
        edges.extend(std::iter::repeat_n((2, 3), loops + 1));
        edges.push((3, 1));
        diagram(&edges)
    }

    #[test]
    fn scalar_chain_preserves_loops_carriers_and_global_factors() {
        for loops in 1..=4 {
            let parent = banana(loops)
                .with_numerator_fragments(
                    [
                        (VertexId(0), SymbolicTensor::infer(Atom::num(2)).unwrap()),
                        (VertexId(2), SymbolicTensor::infer(Atom::num(7)).unwrap()),
                    ],
                    [
                        (EdgeId(0), SymbolicTensor::infer(Atom::num(3)).unwrap()),
                        (EdgeId(2), SymbolicTensor::infer(Atom::num(5)).unwrap()),
                    ],
                )
                .unwrap();
            let insertion = parent
                .propagator_insertions()
                .unwrap()
                .into_iter()
                .find(|i| i.loops == loops)
                .unwrap();
            let result = parent
                .scalar_bubble_surrogate(&insertion, parent.model.particle_id("g").unwrap())
                .unwrap();
            assert_eq!(result.diagram.loop_count(), parent.loop_count());
            assert_eq!(result.diagram.graph.n_nodes(), 2 + loops + 1);
            assert_eq!(result.diagram.graph.n_edges(), 5 + 2 * loops);
            assert_eq!(result.diagram.numerator().expression(), &Atom::num(6));
            assert_eq!(result.diagram.overall_factor(), parent.overall_factor());
            assert_eq!(
                result.diagram.numerator_prefactor(),
                parent.numerator_prefactor()
            );
            assert_eq!(result.diagram.projector(), parent.projector());
            assert_eq!(result.diagram.symmetry_factor(), parent.symmetry_factor());
            for boundary in insertion.boundary {
                let mapped = result.relabeling.translate_half_edge(boundary).unwrap();
                assert_eq!(mapped.endpoint, boundary.endpoint);
                let old = &parent.graph[linnet::half_edge::involution::EdgeIndex(boundary.edge.0)];
                let new =
                    &result.diagram.graph[linnet::half_edge::involution::EdgeIndex(mapped.edge.0)];
                assert_eq!(new.particle, old.particle);
                assert_eq!(new.propagator_power, old.propagator_power);
                assert_eq!(new.directed, old.directed);
            }
            let coefficient = SymbolicTensor::infer(Atom::num(19)).unwrap();
            let attached = result.with_coefficient(&coefficient).unwrap();
            assert_eq!(attached.numerator().expression(), &Atom::num(114));
            attached.validate().unwrap();
            let roundtrip =
                FeynmanDiagram::from_json(parent.model_arc(), &attached.to_json().unwrap())
                    .unwrap();
            assert_eq!(roundtrip.to_json().unwrap(), attached.to_json().unwrap());
        }
    }

    #[test]
    fn ghost_and_fermion_carriers_preserve_orientation_with_scalar_numerators() {
        use linnet::half_edge::{EdgeAccessors, involution::EdgeIndex};

        let parent = banana(2);
        let insertion = parent.propagator_insertions().unwrap().remove(0);
        let scalar = parent
            .scalar_bubble_surrogate(&insertion, parent.model.particle_id("g").unwrap())
            .unwrap();
        let denominator = scalar.diagram.denominator_expression().unwrap();
        for name in ["c", "c~", "q", "q~"] {
            let carrier = parent.model.particle_id(name).unwrap();
            let result = parent.scalar_bubble_surrogate(&insertion, carrier).unwrap();
            let expected = if parent
                .model
                .particle_by_id(carrier)
                .unwrap()
                .is_antiparticle()
            {
                Orientation::Reversed
            } else {
                Orientation::Default
            };
            let artificial = result
                .diagram
                .edges()
                .filter(|(id, _, _)| {
                    !parent
                        .edges()
                        .any(|(old, _, _)| result.relabeling.edge(old) == Some(*id))
                })
                .collect::<Vec<_>>();
            assert_eq!(artificial.len(), 2 * insertion.loops);
            for (id, _, edge) in artificial {
                assert_eq!(edge.particle, carrier);
                assert!(edge.directed, "{name} needs directed model metadata");
                assert_eq!(result.diagram.graph.orientation(EdgeIndex(id.0)), expected);
                assert_eq!(edge.propagator_power, 1);
                assert_eq!(edge.numerator.expression(), &Atom::one());
            }
            assert_eq!(result.diagram.numerator(), scalar.diagram.numerator());
            assert_eq!(result.diagram.overall_factor(), parent.overall_factor());
            assert_eq!(
                result.diagram.denominator_expression().unwrap(),
                denominator
            );
            result.diagram.validate().unwrap();
            let restored =
                FeynmanDiagram::from_dot(parent.model_arc(), &result.diagram.to_dot().unwrap())
                    .unwrap();
            assert_eq!(
                restored.to_json().unwrap(),
                result.diagram.to_json().unwrap()
            );
        }
    }

    #[test]
    fn open_tensor_owners_survive_excision_surrogate_and_restoration() {
        let mut parent = banana(2);
        // Reverse a boundary carrier: no orientation-dependent minus sign may
        // be invented by opening, substituting or restoring its momentum.
        parent = parent.reverse_edge(EdgeId(1)).unwrap();
        let insertion = parent
            .propagator_insertions()
            .unwrap()
            .into_iter()
            .find(|i| i.loops == 2)
            .unwrap();
        let [left, right] = insertion.boundary;
        let h1 = parent.half_edge_id(left).unwrap().0;
        let h2 = parent.half_edge_id(right).unwrap().0;
        let coefficient = NumeratorInterface::admit(symbolica::parse!(&format!(
            "gammalooprs::Q({},spenso::mink(4,gammalooprs::hedge({h1},1)))*gammalooprs::Q({},spenso::mink(4,gammalooprs::hedge({h2},1)))",left.edge.0,right.edge.0)),None).unwrap();
        let excised = parent.excise_with_mapping(&insertion.support, &[]).unwrap();
        let child = excised
            .relabeling
            .translate_tensor(&coefficient, None)
            .unwrap();
        let restored = excised.relabeling.restore_tensor(&child).unwrap();
        assert_eq!(
            NumeratorInterface::export(&restored).0,
            NumeratorInterface::export(&coefficient).0
        );
        assert_eq!(
            NumeratorInterface::export(&restored).1.unwrap().to_json(),
            NumeratorInterface::export(&coefficient)
                .1
                .unwrap()
                .to_json()
        );
        let result = parent
            .scalar_bubble_surrogate(&insertion, parent.model.particle_id("g").unwrap())
            .unwrap();
        let attached = result.with_coefficient(&restored).unwrap();
        let returned = result
            .relabeling
            .restore_tensor(attached.numerator())
            .unwrap();
        assert_eq!(
            NumeratorInterface::export(&returned).0,
            NumeratorInterface::export(&coefficient).0
        );
        assert_eq!(
            NumeratorInterface::export(&returned).1.unwrap().to_json(),
            NumeratorInterface::export(&coefficient)
                .1
                .unwrap()
                .to_json()
        );
        let zero = coefficient.try_sub(&coefficient).unwrap();
        let attached_zero = result.with_coefficient(&zero).unwrap();
        let returned_zero = result
            .relabeling
            .restore_tensor(attached_zero.numerator())
            .unwrap();
        assert!(returned_zero.expression().is_zero());
        assert_eq!(
            NumeratorInterface::export(&returned_zero)
                .1
                .unwrap()
                .to_json(),
            NumeratorInterface::export(&zero).1.unwrap().to_json()
        );
        for (index, edge) in excised
            .diagram
            .loop_momentum_basis
            .external_edges
            .iter()
            .enumerate()
        {
            let mapped = function!(symbols::external_momentum(), index);
            let restored = excised.relabeling.restore_atom(&mapped).unwrap();
            let expected = excised
                .relabeling
                .restore_atom(&function!(symbols::momentum(), edge.0))
                .unwrap();
            assert_eq!(restored, expected);
        }
        let artificial = result
            .diagram
            .edges()
            .find(|(id, _, _)| {
                !parent
                    .edges()
                    .any(|(old, _, _)| result.relabeling.edge(old) == Some(*id))
            })
            .unwrap()
            .0;
        assert!(
            result
                .relabeling
                .restore_atom(&function!(symbols::momentum(), artificial.0))
                .is_err()
        );
    }

    #[test]
    fn outside_fragments_and_global_tensor_payloads_are_relabelled() {
        let mut parent = banana(2);
        let boundary = DiagramHalfEdge {
            edge: EdgeId(5),
            endpoint: DiagramEndpoint::Source,
        };
        let outside = DiagramHalfEdge {
            edge: EdgeId(7),
            endpoint: DiagramEndpoint::Source,
        };
        let inner_hedge = parent.half_edge_id(boundary).unwrap().0;
        let outer_hedge = parent.half_edge_id(outside).unwrap().0;
        let tensor = |edge: usize, hedge: usize, slot: usize| {
            NumeratorInterface::admit(
                symbolica::parse!(&format!(
                    "gammalooprs::Q({edge},spenso::mink(4,gammalooprs::hedge({hedge},{slot})))"
                )),
                None,
            )
            .unwrap()
        };
        let fragment = tensor(5, outer_hedge, 1);
        let prefactor = tensor(7, outer_hedge, 2);
        let projector = tensor(5, inner_hedge, 3);
        parent = parent
            .with_numerator_fragments([(VertexId(1), fragment.clone())], [])
            .unwrap()
            .with_numerator_prefactor(prefactor.clone())
            .with_projector(projector.clone());
        let insertion = parent
            .propagator_insertions()
            .unwrap()
            .into_iter()
            .find(|i| i.loops == 2)
            .unwrap();
        let result = parent
            .scalar_bubble_surrogate(&insertion, parent.model.particle_id("g").unwrap())
            .unwrap();
        for (original, derived) in [
            (&fragment, result.diagram.numerator()),
            (&prefactor, result.diagram.numerator_prefactor()),
            (&projector, result.diagram.projector()),
        ] {
            let returned = result.relabeling.restore_tensor(derived).unwrap();
            assert_eq!(
                NumeratorInterface::export(&returned).0,
                NumeratorInterface::export(original).0
            );
            assert_eq!(
                NumeratorInterface::export(&returned).1.unwrap().to_json(),
                NumeratorInterface::export(original).1.unwrap().to_json()
            );
        }
    }

    #[test]
    fn disjoint_insertions_transport_for_sequential_replacement() {
        let parent = diagram(&[
            (0, 1),
            (0, 2),
            (2, 3),
            (2, 3),
            (3, 4),
            (4, 5),
            (4, 5),
            (5, 1),
        ]);
        let insertions = parent.propagator_insertions().unwrap();
        assert_eq!(insertions.len(), 2);
        let carrier = parent.model.particle_id("g").unwrap();
        let first = parent
            .scalar_bubble_surrogate(&insertions[0], carrier)
            .unwrap();
        let original = &insertions[1];
        let translated = PropagatorInsertion {
            support: first
                .relabeling
                .translate_subgraph(&original.support)
                .unwrap(),
            boundary: original
                .boundary
                .map(|h| first.relabeling.translate_half_edge(h).unwrap()),
            loops: original.loops,
        };
        let second = first
            .diagram
            .scalar_bubble_surrogate(&translated, carrier)
            .unwrap();
        assert_eq!(second.diagram.loop_count(), parent.loop_count());
        assert!(
            first
                .relabeling
                .translate_subgraph(&insertions[0].support)
                .is_err()
        );
        assert!(
            first
                .relabeling
                .translate_subgraph(&SuBitGraph::empty(1))
                .is_err()
        );
    }

    #[test]
    fn invalid_region_and_massive_carrier_are_rejected() {
        let parent = banana(1);
        let insertion = parent.propagator_insertions().unwrap().remove(0);
        let massive = parent.model.particle_id("h").unwrap();
        assert!(parent.scalar_bubble_surrogate(&insertion, massive).is_err());
        let mut invalid = insertion.clone();
        invalid.loops += 1;
        let carrier = parent.model.particle_id("g").unwrap();
        assert!(parent.scalar_bubble_surrogate(&invalid, carrier).is_err());
        let mut partial = insertion;
        let hedge = partial.support.included_iter().next().unwrap();
        partial.support.sub(hedge);
        assert!(parent.scalar_bubble_surrogate(&partial, carrier).is_err());
    }
}
