//! Independent diagrams obtained by opening a selected region's boundary.

use std::collections::{BTreeMap, BTreeSet};

use linnet::half_edge::{
    NodeIndex,
    involution::{Flow, HedgePair},
    subgraph::{Inclusion, SuBitGraph, SubSetLike},
};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};

use crate::{
    DiagramEndpoint, DiagramError, DiagramHalfEdge, ExternalLeg, ExternalState, FeynmanDiagram,
    VertexId, symbols,
};

impl FeynmanDiagram {
    /// Excise an independent diagram, preserving the valence of retained vertices.
    ///
    /// Vertices touched by the selection retain their whole incident crown. Only
    /// internal edges with both halves selected remain paired; every other
    /// incident half becomes a dangling boundary leg. In particular, an omitted
    /// edge joining two retained vertices becomes two independent external legs.
    /// Sewn initial-state carriers are opened using the amplitude particle-flow
    /// convention. No external wavefunctions are inserted at new boundaries.
    ///
    /// Local numerator fragments follow [`Self::numerator_of`]. Half-edge, edge
    /// and vertex indices are relabeled, and an independent momentum basis is
    /// computed. Ambiguous references to an edge split into two legs at the same
    /// vertex are rejected rather than assigned to an arbitrary leg.
    ///
    /// A proper region has unit overall factor, numerator prefactor, projector
    /// and symmetry factor, and does not inherit whole-diagram cuts. Selecting
    /// every half-edge and isolated vertex returns an exact independent clone,
    /// including the original routing and global metadata.
    pub fn excise(
        &self,
        subgraph: &SuBitGraph,
        isolated_vertices: &[VertexId],
    ) -> Result<Self, DiagramError> {
        let error = |message: String| DiagramError::Invariant {
            operation: "excising a diagram region",
            message,
        };
        if subgraph.size() != self.graph.n_hedges() {
            return Err(error("selection has a different half-edge universe".into()));
        }
        let mut retained = self
            .graph
            .nodes(subgraph)
            .into_iter()
            .collect::<BTreeSet<_>>();
        for &VertexId(vertex) in isolated_vertices {
            if vertex >= self.graph.n_nodes() {
                return Err(DiagramError::UnknownVertex {
                    vertex,
                    vertices: self.graph.n_nodes(),
                });
            }
            let node = NodeIndex(vertex);
            if self.graph.iter_crown(node).next().is_some() {
                return Err(error(format!("vertex {vertex} is not isolated")));
            }
            retained.insert(node);
        }
        if subgraph.n_included() == self.graph.n_hedges() && retained.len() == self.graph.n_nodes()
        {
            return Ok(self.clone());
        }

        let mut contributing = BTreeSet::new();
        for (pair, _, _) in self.graph.iter_edges_of(subgraph) {
            if let HedgePair::Paired { source, sink } = pair {
                contributing.extend([self.graph.node_id(source), self.graph.node_id(sink)]);
            }
        }
        let mut builder = Self::builder(self.model_arc(), format!("{}_subgraph", self.name));
        let mut vertices = BTreeMap::new();
        let mut original_vertices = Vec::new();
        for node in retained.iter().copied() {
            let mut vertex = self.graph[node].clone();
            if !contributing.contains(&node)
                && !self.graph.iter_crown(node).all(|hedge| {
                    subgraph.includes(&hedge) || self.graph[self.graph[&hedge]].is_dummy
                })
            {
                vertex.numerator = Atom::one();
            }
            vertices.insert(node.0, builder.add_vertex(vertex).0);
            original_vertices.push(VertexId(node.0));
        }

        let mut edges: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        let mut vertex_edges: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
        let mut half_edges = BTreeMap::new();
        let mut orientations = Vec::new();
        let mut external_ids = self
            .edges()
            .filter_map(|(_, _, edge)| edge.external.as_ref())
            .flat_map(|leg| [leg.index, leg.connection])
            .collect::<BTreeSet<_>>();
        let mut next_external = 0_usize;
        for (pair, id, data) in self.graph.iter_edges() {
            let keep_paired =
                pair.is_paired() && subgraph.includes(&pair) && data.data.external.is_none();
            let attached = match pair {
                HedgePair::Paired { source, sink } => vec![source, sink],
                HedgePair::Unpaired { hedge, .. } => vec![hedge],
                HedgePair::Split { .. } => unreachable!("a full graph has no split edges"),
            };
            let attached = attached
                .into_iter()
                .filter(|hedge| retained.contains(&self.graph.node_id(*hedge)))
                .collect::<Vec<_>>();
            let groups = if keep_paired {
                vec![attached]
            } else {
                attached.into_iter().map(|hedge| vec![hedge]).collect()
            };
            for group in groups {
                let mut edge = data.data.clone();
                if !keep_paired {
                    let hedge = group[0];
                    edge.numerator = Atom::one();
                    if pair.is_paired() && edge.external.is_some() {
                        // A sewn carrier's particle assignment uses the opposite
                        // convention from each of its amplitude half-edges.
                        edge.particle = self.model.particle_by_id(edge.particle)?.antiparticle;
                    }
                    if !edge.is_dummy && pair.is_paired() {
                        while !external_ids.insert(next_external) {
                            next_external = next_external
                                .checked_add(1)
                                .ok_or_else(|| error("external index overflow".into()))?;
                        }
                        edge.external = Some(ExternalLeg {
                            name: format!("boundary_{}", hedge.0),
                            index: next_external,
                            state: match self.graph.flow(hedge) {
                                Flow::Source => ExternalState::Outgoing,
                                Flow::Sink => ExternalState::Incoming,
                            },
                            connection: next_external,
                        });
                    }
                }
                let mut source = None;
                let mut target = None;
                for &hedge in &group {
                    let vertex = Some(VertexId(vertices[&self.graph.node_id(hedge).0]));
                    match self.graph.flow(hedge) {
                        Flow::Source => source = vertex,
                        Flow::Sink => target = vertex,
                    }
                }
                let orientation = if pair.is_paired() && data.data.external.is_some() {
                    data.orientation.reverse()
                } else {
                    data.orientation
                };
                let (source_slot, target_slot) = (edge.source_slot(), edge.target_slot());
                let new =
                    builder.add_edge_with_slots(source, target, edge, source_slot, target_slot)?;
                orientations.push(orientation);
                edges.entry(id.0).or_default().push(new.0);
                for hedge in group {
                    vertex_edges
                        .entry((self.graph.node_id(hedge).0, id.0))
                        .or_default()
                        .push(new.0);
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
        }
        for candidates in vertex_edges.values_mut() {
            candidates.dedup();
        }
        builder.half_edge_order = Some(half_edges.values().copied().collect());
        builder.edge_orientations = Some(orientations);
        let mut diagram = builder.build()?;
        let hedges = half_edges
            .into_iter()
            .map(|(old, new)| (old.0, diagram.half_edge_id(new).unwrap().0))
            .collect::<BTreeMap<_, _>>();
        let translate = |atom: &Atom, vertex: Option<VertexId>| -> Result<Atom, DiagramError> {
            let mut failure = None;
            let result = atom.replace_map_bottom_up(|term, _, out| {
                let AtomView::Fun(function) = term else {
                    return;
                };
                let mut head = function.get_symbol();
                let args = function.iter().collect::<Vec<_>>();
                let Some(mut owner) = args.first().and_then(|arg| usize::try_from(*arg).ok())
                else {
                    return;
                };
                let mapped = if head == symbols::hedge_index() {
                    hedges.get(&owner).copied()
                } else if head == symbols::vertex_index() {
                    vertices.get(&owner).copied()
                } else {
                    if head == symbols::loop_momentum() || head == symbols::external_momentum() {
                        let coordinates = if head == symbols::loop_momentum() {
                            &self.loop_momentum_basis.loop_edges
                        } else {
                            &self.loop_momentum_basis.external_edges
                        };
                        let Some(edge) = coordinates.get(owner) else {
                            failure = Some(format!("unknown parent momentum coordinate {owner}"));
                            return;
                        };
                        owner = edge.0;
                        head = symbols::momentum();
                    }
                    if ![
                        symbols::edge_index(),
                        symbols::momentum(),
                        symbols::denominator(),
                        symbols::u(),
                        symbols::ubar(),
                        symbols::v(),
                        symbols::vbar(),
                        symbols::epsilon(),
                        symbols::epsilonbar(),
                    ]
                    .contains(&head)
                    {
                        return;
                    }
                    let candidates = vertex
                        .and_then(|v| vertex_edges.get(&(v.0, owner)))
                        .or_else(|| edges.get(&owner));
                    match candidates.map(Vec::as_slice) {
                        Some([edge]) => Some(*edge),
                        Some(_) => {
                            failure = Some(format!(
                                "symbolic edge {owner} has multiple boundary copies"
                            ));
                            return;
                        }
                        None => None,
                    }
                };
                let Some(mapped) = mapped else {
                    failure = Some(format!(
                        "symbolic index {owner} lies outside the excised region"
                    ));
                    return;
                };
                **out = FunctionBuilder::new(head)
                    .add_arg(mapped)
                    .add_args(args[1..].iter().copied())
                    .finish();
            });
            failure.map_or(Ok(result), |message| Err(error(message)))
        };
        diagram.graph = diagram.graph.map_result(
            |_, id, mut vertex| {
                vertex.numerator = translate(&vertex.numerator, Some(original_vertices[id.0]))?;
                Ok::<_, DiagramError>(vertex)
            },
            |_, _, _, _, mut data| {
                data.data.numerator = translate(&data.data.numerator, None)?;
                Ok(data)
            },
            |_, data| Ok(data),
        )?;
        diagram.numerator = diagram
            .vertices()
            .map(|(_, vertex)| vertex.numerator.clone())
            .chain(diagram.edges().map(|(_, _, edge)| edge.numerator.clone()))
            .product();
        diagram.validate()?;
        Ok(diagram)
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use feynkit_model::Model;
    use linnet::half_edge::{
        involution::{EdgeIndex, Hedge},
        subgraph::{BaseSubgraph, ModifySubSet},
    };
    use symbolica::function;

    use super::*;
    use crate::{DiagramEdge, DiagramVertex, EdgeId};

    fn model() -> Arc<Model> {
        Arc::new(Model::from_json(r#"{
            "name":"excision_phi3","restriction":null,"orders":[],
            "parameters":[{"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null}],
            "particles":[{"pdg_code":25,"name":"phi","antiname":"phi","spin":1,"color":1,"mass":"ZERO","width":"ZERO","texname":"phi","antitexname":"phi","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0}],
            "propagators":[],"lorentz_structures":[{"name":"L1","spins":[1,1,1],"structure":"1"}],"couplings":[],
            "vertex_rules":[{"name":"V_1","particles":["phi","phi","phi"],"color_structures":["1"],"lorentz_structures":["L1"],"couplings":[[null]]}]
        }"#).unwrap())
    }

    fn bubble() -> FeynmanDiagram {
        let model = model();
        let particle = model.particle_id("phi").unwrap();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "bubble")
            .overall_factor(Atom::num(7))
            .numerator_prefactor(Atom::num(11))
            .projector(Atom::num(13))
            .symmetry_factor(2)
            .numerator(Atom::num(6));
        let mut left = DiagramVertex::interaction("left", rule);
        left.numerator = Atom::num(2);
        let left = builder.add_vertex(left);
        let mut right = DiagramVertex::interaction("right", rule);
        right.numerator = Atom::num(3);
        let right = builder.add_vertex(right);
        for (index, state, source, target) in [
            (0, ExternalState::Incoming, None, Some(left)),
            (1, ExternalState::Outgoing, Some(right), None),
        ] {
            let mut edge = DiagramEdge::new(particle, false);
            edge.external = Some(ExternalLeg {
                name: format!("p{index}"),
                index,
                state,
                connection: index,
            });
            builder.add_edge(source, target, edge).unwrap();
        }
        for _ in 0..2 {
            builder
                .add_edge(left, right, DiagramEdge::new(particle, false))
                .unwrap();
        }
        let diagram = builder.build().unwrap();
        diagram.validate().unwrap();
        diagram
    }

    #[test]
    fn full_excision_preserves_routing_and_global_metadata() {
        let parent = bubble();
        let result = parent.excise(&parent.graph.full_filter(), &[]).unwrap();
        assert_eq!(result.to_json().unwrap(), parent.to_json().unwrap());
        assert_eq!(result.loop_momentum_basis(), parent.loop_momentum_basis());
    }

    #[test]
    fn excision_opens_omitted_parallel_edge_without_losing_vertex_slots() {
        let parent = bubble();
        let mut selected: SuBitGraph = parent.graph.empty_subgraph();
        selected.add(parent.graph[&EdgeIndex(2)].1);
        let result = parent.excise(&selected, &[]).unwrap();
        assert_eq!(result.graph.n_nodes(), 2);
        assert_eq!(result.graph.n_edges(), 5);
        assert_eq!(result.graph.n_internals(), 1);
        assert_eq!(result.graph.n_hedges(), parent.graph.n_hedges());
        assert!(result.loop_momentum_basis().loop_edges.is_empty());
        assert_eq!(
            result.numerator(),
            &parent.numerator_of(&selected, &parent.graph.empty_subgraph::<SuBitGraph>())
        );
        assert_eq!(result.overall_factor(), &Atom::one());
        assert_eq!(result.numerator_prefactor(), &Atom::one());
        assert_eq!(result.projector(), &Atom::one());
        assert_eq!(result.symmetry_factor(), 1);
        assert!(result.cuts().is_empty());
        result.validate().unwrap();
        let restored =
            FeynmanDiagram::from_json(parent.model_arc(), &result.to_json().unwrap()).unwrap();
        assert_eq!(restored.to_json().unwrap(), result.to_json().unwrap());
    }

    #[test]
    fn single_half_edge_excision_has_complete_interaction_and_local_numerator() {
        let parent = bubble();
        let mut selected: SuBitGraph = parent.graph.empty_subgraph();
        selected.add(
            parent
                .half_edge_id(DiagramHalfEdge {
                    edge: EdgeId(2),
                    endpoint: DiagramEndpoint::Target,
                })
                .unwrap(),
        );
        let result = parent.excise(&selected, &[]).unwrap();
        assert_eq!(result.graph.n_nodes(), 1);
        assert_eq!(result.graph.n_edges(), 3);
        assert_eq!(result.graph.n_internals(), 0);
        assert_eq!(result.numerator(), &Atom::one());
        assert!(result.edges().all(|(_, _, edge)| edge.external.is_some()));
        result.validate().unwrap();
    }

    #[test]
    fn disconnected_excision_routes_each_boundary_component() {
        let parent = bubble();
        let mut selected: SuBitGraph = parent.graph.empty_subgraph();
        selected.add(parent.graph[&EdgeIndex(0)].1);
        selected.add(parent.graph[&EdgeIndex(1)].1);
        let result = parent.excise(&selected, &[]).unwrap();
        assert_eq!(
            result
                .graph
                .count_connected_components(&result.graph.full_filter()),
            2
        );
        assert_eq!(result.graph.n_edges(), 6);
        assert_eq!(result.loop_momentum_basis().dependent_externals.len(), 2);
        result.validate().unwrap();
    }

    #[test]
    fn excision_remaps_tensor_owners_and_parent_routing_coordinates() {
        let mut parent = bubble();
        let selected = SuBitGraph::from_hedge_iter(
            parent.graph.iter_crown(NodeIndex(1)),
            parent.graph.n_hedges(),
        );
        let old_hedge = parent.graph.iter_crown(NodeIndex(1)).last().unwrap();
        let old_edge = parent.graph[&old_hedge].0;
        let old_loop = parent.loop_momentum_basis.loop_edges[0].0;
        parent.graph[NodeIndex(1)].numerator = function!(symbols::hedge_index(), old_hedge.0, 1)
            * function!(symbols::vertex_index(), 1, 2)
            * function!(symbols::edge_index(), old_edge, 3)
            * function!(symbols::loop_momentum(), 0);
        parent.numerator = Atom::num(2) * &parent.graph[NodeIndex(1)].numerator;
        let result = parent.excise(&selected, &[]).unwrap();
        let new_hedge = result.graph.iter_crown(NodeIndex(0)).last().unwrap();
        let new_edge = result.graph[&new_hedge].0;
        let new_loop_edge = result
            .edges()
            .find_map(|(edge, _, data)| {
                (data.external.as_ref().unwrap().name
                    == format!(
                        "boundary_{}",
                        parent
                            .half_edge_id(DiagramHalfEdge {
                                edge: EdgeId(old_loop),
                                endpoint: DiagramEndpoint::Target
                            })
                            .unwrap()
                            .0
                    ))
                .then_some(edge.0)
            })
            .unwrap();
        let expected = function!(symbols::hedge_index(), new_hedge.0, 1)
            * function!(symbols::vertex_index(), 0, 2)
            * function!(symbols::edge_index(), new_edge, 3)
            * function!(symbols::momentum(), new_loop_edge);
        assert_eq!(result.numerator(), &expected);
        result.validate().unwrap();
    }

    #[test]
    fn isolated_vertices_and_empty_excision_are_independent() {
        let mut builder = FeynmanDiagram::builder(model(), "isolated");
        for name in ["first", "second"] {
            builder.add_vertex(DiagramVertex {
                name: name.into(),
                interaction: None,
                numerator: Atom::one(),
            });
        }
        let parent = builder.build().unwrap();
        let empty: SuBitGraph = parent.graph.empty_subgraph();
        let none = parent.excise(&empty, &[]).unwrap();
        let one = parent.excise(&empty, &[VertexId(1)]).unwrap();
        assert_eq!(none.graph.n_nodes(), 0);
        assert_eq!(one.graph.n_nodes(), 1);
        assert_eq!(one.vertex(VertexId(0)).unwrap().name, "second");
        none.validate().unwrap();
        one.validate().unwrap();
        assert!(parent.excise(&empty, &[VertexId(2)]).is_err());
        assert!(parent.excise(&SuBitGraph::empty(1), &[]).is_err());
        let parent = bubble();
        assert!(
            parent
                .excise(&parent.graph.empty_subgraph::<SuBitGraph>(), &[VertexId(0)])
                .is_err()
        );
    }

    #[test]
    fn excision_splits_self_loop_into_separate_boundary_legs() {
        let model = model();
        let particle = model.particle_id("phi").unwrap();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "tadpole");
        let vertex = builder.add_vertex(DiagramVertex::interaction("v", rule));
        builder
            .add_edge(vertex, vertex, DiagramEdge::new(particle, false))
            .unwrap();
        let mut external = DiagramEdge::new(particle, false);
        external.external = Some(ExternalLeg {
            name: "p".into(),
            index: 0,
            state: ExternalState::Incoming,
            connection: 0,
        });
        builder.add_edge(None, vertex, external).unwrap();
        let parent = builder.build().unwrap();
        let mut selected = SuBitGraph::empty(parent.graph.n_hedges());
        selected.add(Hedge(2));
        let result = parent.excise(&selected, &[]).unwrap();
        assert_eq!(result.graph.n_hedges(), 3);
        assert_eq!(result.graph.n_edges(), 3);
        assert_eq!(result.graph.n_internals(), 0);
        assert_eq!(result.loop_momentum_basis.external_edges.len(), 3);
        result.validate().unwrap();
    }

    #[test]
    fn excision_opens_sewn_fermions_with_amplitude_particle_flow() {
        let model = Arc::new(Model::from_json(r#"{
            "name":"excision_fermions","restriction":null,"orders":[],
            "parameters":[{"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null}],
            "particles":[
                {"pdg_code":1,"name":"f","antiname":"f~","spin":2,"color":1,"mass":"ZERO","width":"ZERO","texname":"f","antitexname":"fbar","charge":0.0,"ghost_number":0,"lepton_number":1,"y_charge":0},
                {"pdg_code":-1,"name":"f~","antiname":"f","spin":2,"color":1,"mass":"ZERO","width":"ZERO","texname":"fbar","antitexname":"f","charge":0.0,"ghost_number":0,"lepton_number":-1,"y_charge":0}
            ],
            "propagators":[],"lorentz_structures":[{"name":"L","spins":[2,2],"structure":"1"}],"couplings":[],
            "vertex_rules":[
                {"name":"V","particles":["f","f~"],"color_structures":["1"],"lorentz_structures":["L"],"couplings":[[null]]},
                {"name":"Vbar","particles":["f~","f"],"color_structures":["1"],"lorentz_structures":["L"],"couplings":[[null]]}
            ]
        }"#).unwrap());
        let fermion = model.particle_id("f").unwrap();
        let antifermion = model.particle_id("f~").unwrap();
        let mut builder = FeynmanDiagram::builder(model.clone(), "sewn");
        let left = builder.add_vertex(DiagramVertex::interaction(
            "left",
            model.vertex_rule_id("V").unwrap(),
        ));
        let right = builder.add_vertex(DiagramVertex::interaction(
            "right",
            model.vertex_rule_id("Vbar").unwrap(),
        ));
        let mut carrier = DiagramEdge::new(fermion, true);
        carrier.external = Some(ExternalLeg {
            name: "p".into(),
            index: 0,
            state: ExternalState::Incoming,
            connection: 0,
        });
        builder.add_edge(left, right, carrier).unwrap();
        builder
            .add_edge(left, right, DiagramEdge::new(fermion, true))
            .unwrap();
        let source = |edge| DiagramHalfEdge {
            edge: EdgeId(edge),
            endpoint: DiagramEndpoint::Source,
        };
        let target = |edge| DiagramHalfEdge {
            edge: EdgeId(edge),
            endpoint: DiagramEndpoint::Target,
        };
        let parent = builder
            .build()
            .unwrap()
            .with_cut_partitions(vec![(
                vec![target(0), target(1)],
                vec![source(0), source(1)],
            )])
            .unwrap();
        parent.validate().unwrap();
        let mut selected = parent.graph.empty_subgraph::<SuBitGraph>();
        selected.add(parent.graph[&EdgeIndex(1)].1);
        let result = parent.excise(&selected, &[]).unwrap();
        assert_eq!(result.graph.n_internals(), 1);
        assert_eq!(result.graph.n_edges(), 3);
        assert_eq!(
            result
                .edges()
                .filter(|(_, _, edge)| edge.external.is_some() && edge.particle == antifermion)
                .count(),
            2
        );
        assert!(result.cuts().is_empty());
        result.validate().unwrap();
    }

    #[test]
    fn excision_keeps_dummy_boundaries_out_of_momentum_routing() {
        let model = model();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "dummy_boundary");
        let left = builder.add_vertex(DiagramVertex {
            name: "left".into(),
            interaction: None,
            numerator: Atom::one(),
        });
        let right = builder.add_vertex(DiagramVertex {
            name: "right".into(),
            interaction: None,
            numerator: Atom::one(),
        });
        builder
            .add_edge(left, right, DiagramEdge::new(particle, false))
            .unwrap();
        let mut dummy = DiagramEdge::new(particle, false);
        dummy.is_dummy = true;
        builder.add_edge(left, None, dummy).unwrap();
        let mut external = DiagramEdge::new(particle, false);
        external.external = Some(ExternalLeg {
            name: "p".into(),
            index: 0,
            state: ExternalState::Incoming,
            connection: 0,
        });
        builder.add_edge(None, right, external).unwrap();
        let parent = builder.build().unwrap();
        parent.validate().unwrap();
        let mut selected = parent.graph.empty_subgraph::<SuBitGraph>();
        selected.add(parent.graph[&EdgeIndex(0)].1);
        let result = parent.excise(&selected, &[]).unwrap();
        assert_eq!(result.graph.n_edges(), 3);
        let (id, _, dummy) = result.edges().find(|(_, _, edge)| edge.is_dummy).unwrap();
        assert!(dummy.external.is_none());
        assert!(!result.loop_momentum_basis.external_edges.contains(&id));
        let signature = &result.loop_momentum_basis.edge_signatures[&id];
        assert!(
            signature
                .loops
                .iter()
                .chain(signature.external.iter())
                .all(|sign| sign.is_zero())
        );
        result.validate().unwrap();
    }
}
