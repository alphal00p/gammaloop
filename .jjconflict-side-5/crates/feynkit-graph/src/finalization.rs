//! Finalize generated topology once, in GammaLoop's native half-edge convention.
use super::*;
use linnet::{
    half_edge::{
        EdgeAccessors,
        involution::Hedge,
        subgraph::{ModifySubSet, OrientedCut},
        swap::Swap,
    },
    permutation::Permutation,
};
use symbolica::atom::{AtomView, FunctionBuilder};

#[derive(Clone)]
struct PendingEdge {
    original: EdgeId,
    data: DiagramEdge,
}

impl FeynmanDiagramBuilder {
    pub(super) fn finalize(self) -> Result<FeynmanDiagram, DiagramError> {
        let generated = !self.generation_externals.is_empty();
        let mut external_indices = BTreeSet::new();
        for (vertex, external) in &self.generation_externals {
            if !external_indices.insert(external.index) {
                return Err(DiagramError::DuplicateExternalIndex(external.index));
            }
            let degree = self
                .edges
                .iter()
                .map(|(ends, _)| {
                    usize::from(ends.source == Some(*vertex))
                        + usize::from(ends.target == Some(*vertex))
                })
                .sum();
            if degree != 1 {
                return Err(DiagramError::InvalidExternalDegree {
                    vertex: vertex.0,
                    degree,
                });
            }
            if self.vertices[vertex.0].numerator != Atom::one() {
                return Err(DiagramError::ExternalVertexNumerator { vertex: vertex.0 });
            }
        }
        let mut builder = HedgeGraphBuilder::new();
        let mut vertices = BTreeMap::new();
        for (id, vertex) in self.vertices.iter().enumerate() {
            if !self.generation_externals.contains_key(&VertexId(id)) {
                vertices.insert(VertexId(id), builder.add_node(vertex.clone()));
            }
        }
        let mut connections = BTreeMap::<usize, (Option<usize>, Option<usize>)>::new();
        let mut externals = BTreeMap::new();
        let mut internal = Vec::new();
        for (id, (endpoints, edge)) in self.edges.iter().enumerate() {
            let external = endpoints
                .source
                .and_then(|v| self.generation_externals.get(&v))
                .map(|leg| (leg.clone(), endpoints.target, Flow::Sink, edge.target_slot))
                .or_else(|| {
                    endpoints
                        .target
                        .and_then(|v| self.generation_externals.get(&v))
                        .map(|leg| {
                            (
                                leg.clone(),
                                endpoints.source,
                                Flow::Source,
                                edge.source_slot,
                            )
                        })
                });
            if let Some((leg, attachment, old_flow, slot)) = external {
                if edge.numerator != Atom::one() {
                    return Err(DiagramError::ExternalEdgeNumerator { edge: id });
                }
                let attachment = attachment
                    .filter(|v| vertices.contains_key(v))
                    .ok_or(DiagramError::ExternalToExternalEdge { edge: id })?;
                let entry = connections.entry(leg.connection).or_default();
                let target = match leg.state {
                    ExternalState::Incoming => &mut entry.0,
                    ExternalState::Outgoing => &mut entry.1,
                };
                if target.replace(id).is_some() {
                    return Err(DiagramError::InvalidExternalConnectionStates {
                        connection: leg.connection,
                        state: leg.state.as_str(),
                    });
                }
                externals.insert(id, (leg, attachment, old_flow, slot));
            } else {
                internal.push(id);
            }
        }
        let orientation = |edge: &DiagramEdge| -> Result<Orientation, DiagramError> {
            Ok(if !edge.directed {
                Orientation::Undirected
            } else if self.model.particle_by_id(edge.particle)?.is_antiparticle() {
                Orientation::Reversed
            } else {
                Orientation::Default
            })
        };
        let mut logical = BTreeMap::new();
        let mut next_edge = 0;
        let add_external = |id: usize,
                            flow: Flow,
                            builder: &mut HedgeGraphBuilder<PendingEdge, DiagramVertex>|
         -> Result<(), DiagramError> {
            let (leg, vertex, old_flow, slot) = &externals[&id];
            let mut edge = self.edges[id].1.clone();
            edge.external = Some(leg.clone());
            if old_flow != &flow {
                edge.particle = self.model.particle_by_id(edge.particle)?.antiparticle;
            }
            match flow {
                Flow::Source => edge.source_slot = *slot,
                Flow::Sink => edge.target_slot = *slot,
            }
            let direction = orientation(&edge)?;
            builder.add_external_edge(
                vertices[vertex],
                PendingEdge {
                    original: EdgeId(id),
                    data: edge,
                },
                direction,
                flow,
            );
            Ok(())
        };
        // Half-edge insertion order controls channel coordinates. Logical edge
        // identity is assigned separately, as in GammaLoop's finalization.
        if generated {
            for (incoming, _) in connections.values() {
                if let Some(id) = incoming {
                    logical.insert(EdgeId(*id), next_edge);
                    next_edge += 1;
                    add_external(*id, Flow::Sink, &mut builder)?;
                }
            }
            for (incoming, outgoing) in connections.values() {
                if incoming.is_none()
                    && let Some(id) = outgoing
                {
                    logical.insert(EdgeId(*id), next_edge);
                    next_edge += 1;
                    add_external(*id, Flow::Source, &mut builder)?;
                }
            }
            for id in &internal {
                logical.insert(EdgeId(*id), next_edge);
                next_edge += 1;
            }
            internal.sort_by_key(|id| {
                let (endpoints, edge) = &self.edges[*id];
                let (source, target, particle) = if endpoints.source <= endpoints.target {
                    (endpoints.source, endpoints.target, edge.particle)
                } else {
                    (
                        endpoints.target,
                        endpoints.source,
                        self.model
                            .particle_by_id(edge.particle)
                            .expect("resolved particle")
                            .antiparticle,
                    )
                };
                (
                    source,
                    target,
                    self.model
                        .particle_by_id(particle)
                        .expect("resolved particle")
                        .pdg_code,
                    *id,
                )
            });
        }
        for id in internal {
            let (endpoints, edge) = &self.edges[id];
            let data = PendingEdge {
                original: EdgeId(id),
                data: edge.clone(),
            };
            let direction = orientation(edge)?;
            match (endpoints.source, endpoints.target) {
                (Some(source), Some(target)) => {
                    builder.add_edge(vertices[&source], vertices[&target], data, direction);
                }
                (Some(vertex), None) => {
                    builder.add_external_edge(vertices[&vertex], data, direction, Flow::Source);
                }
                (None, Some(vertex)) => {
                    builder.add_external_edge(vertices[&vertex], data, direction, Flow::Sink);
                }
                (None, None) => unreachable!("builder validates attached endpoints"),
            }
            if !generated {
                logical.insert(EdgeId(id), id);
            }
        }
        for (incoming, outgoing) in connections.values() {
            if let (Some(incoming), Some(outgoing)) = (incoming, outgoing) {
                logical.insert(EdgeId(*outgoing), logical[&EdgeId(*incoming)]);
                add_external(*outgoing, Flow::Source, &mut builder)?;
            }
        }
        let mut pending: HedgeGraph<PendingEdge, DiagramVertex> = builder.into();
        if generated {
            pending
                .sew(
                    |_, left, _, right| match (&left.data.data.external, &right.data.data.external)
                    {
                        (Some(left), Some(right)) => {
                            left.connection == right.connection && left.state != right.state
                        }
                        _ => false,
                    },
                    |left_flow, mut left, right_flow, right| match (left_flow, right_flow) {
                        (Flow::Sink, Flow::Source) => {
                            left.data.data.source_slot = left.data.data.target_slot;
                            left.data.data.target_slot = right.data.data.source_slot;
                            (Flow::Source, left)
                        }
                        (Flow::Source, Flow::Sink) => {
                            let mut incoming = right;
                            incoming.data.data.source_slot = incoming.data.data.target_slot;
                            incoming.data.data.target_slot = left.data.data.source_slot;
                            (Flow::Sink, incoming)
                        }
                        _ => unreachable!("external connections have opposite process flows"),
                    },
                )
                .map_err(|error| DiagramError::Invariant {
                    operation: "sewing initial-state edges",
                    message: format!("{error:?}"),
                })?;
            let cut_hedges = pending
                .iter_edges()
                .filter_map(|(pair, _, edge)| {
                    if let HedgePair::Paired { sink, .. } = pair {
                        edge.data
                            .data
                            .external
                            .as_ref()
                            .map(|external| (external.connection, sink))
                    } else {
                        None
                    }
                })
                .collect::<BTreeMap<_, _>>();
            let permutation = Permutation::from_mappings(
                cut_hedges
                    .values()
                    .enumerate()
                    .map(|(target, hedge)| (hedge.0, target)),
                pending.n_hedges(),
            )
            .map_err(|error| DiagramError::Invariant {
                operation: "ordering initial-state half-edges",
                message: error.to_string(),
            })?;
            <HedgeGraph<_, _> as Swap<Hedge>>::permute(&mut pending, &permutation);
        }
        let permutation = Permutation::from_mappings(
            pending
                .iter_edges()
                .map(|(_, id, edge)| (id.0, logical[&edge.data.original])),
            pending.n_edges(),
        )
        .map_err(|error| DiagramError::Invariant {
            operation: "ordering logical edges",
            message: error.to_string(),
        })?;
        <HedgeGraph<_, _> as Swap<EdgeIndex>>::permute(&mut pending, &permutation);
        let mut graph = pending.map_data_ref(
            |_, _, vertex| vertex.clone(),
            |_, _, _, edge| edge.map(|edge| edge.data.clone()),
            |_, data| *data,
        );
        if let Some(orientations) = &self.edge_orientations {
            if orientations.len() != graph.n_edges() {
                return Err(DiagramError::Invariant {
                    operation: "restoring edge orientations",
                    message: "orientation count differs from edge count".into(),
                });
            }
            for (edge, orientation) in orientations.iter().enumerate() {
                graph.set_orientation(EdgeIndex(edge), *orientation);
            }
        }
        if let Some(order) = &self.half_edge_order {
            let mut mappings = Vec::new();
            for (target, half_edge) in order.iter().enumerate() {
                let Some(pair) = (half_edge.edge.0 < graph.n_edges())
                    .then(|| graph[&EdgeIndex(half_edge.edge.0)].1)
                else {
                    return Err(DiagramError::UnknownEdge {
                        edge: half_edge.edge.0,
                        edges: graph.n_edges(),
                    });
                };
                let hedge = match (pair, half_edge.endpoint) {
                    (HedgePair::Paired { source, .. }, DiagramEndpoint::Source)
                    | (
                        HedgePair::Unpaired {
                            hedge: source,
                            flow: Flow::Source,
                        },
                        DiagramEndpoint::Source,
                    ) => source,
                    (HedgePair::Paired { sink, .. }, DiagramEndpoint::Target)
                    | (
                        HedgePair::Unpaired {
                            hedge: sink,
                            flow: Flow::Sink,
                        },
                        DiagramEndpoint::Target,
                    ) => sink,
                    _ => {
                        return Err(DiagramError::Invariant {
                            operation: "restoring half-edge order",
                            message: "serialized endpoint does not exist".into(),
                        });
                    }
                };
                mappings.push((hedge.0, target));
            }
            if order.len() != graph.n_hedges()
                || order.iter().collect::<BTreeSet<_>>().len() != order.len()
            {
                return Err(DiagramError::Invariant {
                    operation: "restoring half-edge order",
                    message: "serialized order must contain every half-edge exactly once".into(),
                });
            }
            let permutation =
                Permutation::from_mappings(mappings, graph.n_hedges()).map_err(|error| {
                    DiagramError::Invariant {
                        operation: "restoring half-edge order",
                        message: error.to_string(),
                    }
                })?;
            <HedgeGraph<_, _> as Swap<Hedge>>::permute(&mut graph, &permutation);
        }
        let mut half_edges = BTreeMap::new();
        let mut signs = BTreeMap::new();
        for (old, (endpoints, _)) in self.edges.iter().enumerate() {
            let edge = EdgeId(old);
            let pair = graph[&EdgeIndex(logical[&edge])].1;
            let candidates = match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    vec![source, sink]
                }
                HedgePair::Unpaired { hedge, .. } => vec![hedge],
            };
            for (vertex, endpoint, flow) in [
                (endpoints.source, DiagramEndpoint::Source, Flow::Source),
                (endpoints.target, DiagramEndpoint::Target, Flow::Sink),
            ] {
                let attached = candidates
                    .iter()
                    .copied()
                    .filter(|hedge| {
                        vertex
                            .and_then(|v| vertices.get(&v))
                            .is_none_or(|v| graph.node_id(*hedge) == *v)
                    })
                    .collect::<Vec<_>>();
                let hedge = attached
                    .iter()
                    .copied()
                    .find(|h| graph.flow(*h) == flow)
                    .or_else(|| attached.first().copied())
                    .or_else(|| candidates.iter().copied().find(|h| graph.flow(*h) == flow))
                    .unwrap_or(candidates[0]);
                half_edges.insert(DiagramHalfEdge { edge, endpoint }, hedge);
            }
            let original_flow = if endpoints.source.is_some_and(|v| vertices.contains_key(&v)) {
                Flow::Source
            } else {
                Flow::Sink
            };
            let endpoint = if original_flow == Flow::Source {
                DiagramEndpoint::Source
            } else {
                DiagramEndpoint::Target
            };
            signs.insert(
                edge,
                if graph.flow(half_edges[&DiagramHalfEdge { edge, endpoint }]) == original_flow {
                    1
                } else {
                    -1
                },
            );
        }
        let translate = |atom: &Atom| -> Result<Atom, DiagramError> {
            if !generated {
                return Ok(atom.clone());
            }
            let mut failure = None;
            let translated = atom.replace_map(|term, _, out| {
                let AtomView::Fun(function) = term else {
                    return;
                };
                let symbol = function.get_symbol();
                let arguments = function.iter().collect::<Vec<_>>();
                let kind = [
                    (symbols::hedge_index(), 0),
                    (symbols::edge_index(), 2),
                    (symbols::vertex_index(), 3),
                    (momentum_symbol(), 4),
                    (symbols::u(), 5),
                    (symbols::ubar(), 5),
                    (symbols::v(), 5),
                    (symbols::vbar(), 5),
                    (symbols::epsilon(), 5),
                    (symbols::epsilonbar(), 5),
                ]
                .into_iter()
                .find_map(|(head, kind)| (head == symbol).then_some(kind));
                let Some(kind) = kind else {
                    return;
                };
                let Some(raw_owner) = arguments
                    .first()
                    .and_then(|owner| usize::try_from(*owner).ok())
                else {
                    return;
                };
                let owner = if kind == 0 { raw_owner / 2 } else { raw_owner };
                let mapped = if kind == 3 {
                    vertices.get(&VertexId(owner)).map(|v| v.0)
                } else {
                    logical.get(&EdgeId(owner)).copied()
                };
                let Some(mapped) = mapped else {
                    failure = Some(format!("unknown expression index owner {owner}"));
                    return;
                };
                let (head, index) = match kind {
                    0 => (
                        symbols::hedge_index(),
                        half_edges[&DiagramHalfEdge {
                            edge: EdgeId(owner),
                            endpoint: if raw_owner % 2 == 0 {
                                DiagramEndpoint::Source
                            } else {
                                DiagramEndpoint::Target
                            },
                        }]
                            .0,
                    ),
                    2 => (symbols::edge_index(), mapped),
                    3 => (symbols::vertex_index(), mapped),
                    5 => (symbol, mapped),
                    _ => (momentum_symbol(), mapped),
                };
                let value = FunctionBuilder::new(head)
                    .add_arg(index)
                    .add_args(arguments[1..].iter().copied())
                    .finish();
                **out = if kind == 4 && signs[&EdgeId(owner)] < 0 {
                    -value
                } else {
                    value
                };
            });
            if let Some(message) = failure {
                return Err(DiagramError::Invariant {
                    operation: "finalizing symbolic indices",
                    message,
                });
            }
            Ok(translated)
        };
        graph = graph.map_data_ref_result(
            |_, _, vertex| {
                let mut vertex = vertex.clone();
                vertex.numerator = translate(&vertex.numerator)?;
                Ok::<_, DiagramError>(vertex)
            },
            |_, _, _, edge| {
                let numerator = translate(&edge.data.numerator)?;
                Ok(edge.map(|edge| {
                    let mut edge = edge.clone();
                    edge.numerator = numerator;
                    edge
                }))
            },
            |(_, data)| Ok(*data),
        )?;
        let mut diagram = FeynmanDiagram {
            model: self.model.clone(),
            id: DiagramId(0),
            name: self.name,
            graph,
            symmetry_factor: self.symmetry_factor,
            overall_factor: translate(&self.overall_factor)?,
            numerator: translate(&self.numerator)?,
            numerator_prefactor: translate(&self.numerator_prefactor)?,
            projector: translate(&self.projector)?,
            loop_momentum_basis: LoopMomentumBasis {
                tree_edges: vec![],
                loop_edges: vec![],
                external_edges: vec![],
                dependent_externals: vec![],
                edge_signatures: BTreeMap::new(),
            },
            cuts: self.cuts.clone(),
            topology_threshold_candidates: self.topology_threshold_candidates.clone(),
        };
        if generated {
            diagram.finalize_cut_indices(&half_edges)?;
        }
        diagram = if let Some(basis) = self.loop_momentum_basis {
            if generated {
                let requested = basis
                    .loop_edges
                    .iter()
                    .map(|id| {
                        logical
                            .get(id)
                            .copied()
                            .map(EdgeId)
                            .ok_or(DiagramError::UnknownEdge {
                                edge: id.0,
                                edges: self.edges.len(),
                            })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                diagram.with_loop_momentum_edges(&requested)?
            } else {
                basis.validate(&diagram)?;
                diagram.loop_momentum_basis = basis;
                diagram
            }
        } else {
            let basis = diagram
                .loop_momentum_bases_with_limit(1)?
                .into_iter()
                .next()
                .ok_or(DiagramError::MissingLoopMomentumBasis)?;
            diagram.loop_momentum_basis = basis;
            diagram
        };
        if generated {
            let ignored = diagram.initial_state_tree().0;
            let excluded = diagram
                .half_edges()
                .filter(|half_edge| {
                    ignored.includes(
                        &diagram
                            .half_edge_id(*half_edge)
                            .expect("existing half-edge"),
                    )
                })
                .collect::<BTreeSet<_>>();
            for candidate in &mut diagram.topology_threshold_candidates {
                // The loop-independent initial-state attachment tree belongs
                // to neither threshold side, as in GammaLoop finalization.
                candidate
                    .left
                    .retain(|half_edge| !excluded.contains(half_edge));
                candidate
                    .right
                    .retain(|half_edge| !excluded.contains(half_edge));
            }
        }
        let mut cuts = std::mem::take(&mut diagram.cuts);
        cuts.sort_by_cached_key(|cut| diagram.cut_order_key(&cut.cut));
        diagram.cuts = cuts;
        let mut thresholds = std::mem::take(&mut diagram.topology_threshold_candidates);
        thresholds.sort_by_cached_key(|candidate| diagram.cut_order_key(&candidate.cut));
        diagram.topology_threshold_candidates = thresholds;
        diagram.id = DiagramId::from_key(diagram.model.fingerprint(), &diagram.structural_key()?)?;
        Ok(diagram)
    }
}

impl FeynmanDiagram {
    pub(super) fn cut_order_key(&self, half_edges: &[DiagramHalfEdge]) -> OrientedCut {
        let mut left = self.graph.empty_subgraph::<SuBitGraph>();
        let mut right = self.graph.empty_subgraph::<SuBitGraph>();
        for half_edge in half_edges {
            if let Some(hedge) = self.half_edge_id(*half_edge) {
                left.add(hedge);
                right.add(self.graph.inv(hedge));
            }
        }
        OrientedCut { left, right }
    }

    fn finalize_cut_indices(
        &mut self,
        mapping: &BTreeMap<DiagramHalfEdge, Hedge>,
    ) -> Result<(), DiagramError> {
        let native = |half_edge: &DiagramHalfEdge| -> DiagramHalfEdge {
            let hedge = mapping[half_edge];
            DiagramHalfEdge {
                edge: EdgeId(self.graph[&hedge].0),
                endpoint: if self.graph.flow(hedge) == Flow::Source {
                    DiagramEndpoint::Source
                } else {
                    DiagramEndpoint::Target
                },
            }
        };
        let initial_left = self
            .graph
            .iter_edges()
            .filter_map(|(pair, id, edge)| {
                (edge.data.external.is_some() && pair.is_paired()).then_some(DiagramHalfEdge {
                    edge: EdgeId(id.0),
                    endpoint: DiagramEndpoint::Target,
                })
            })
            .collect::<BTreeSet<_>>();
        let initial_right = initial_left
            .iter()
            .map(|h| DiagramHalfEdge {
                edge: h.edge,
                endpoint: DiagramEndpoint::Source,
            })
            .collect::<BTreeSet<_>>();
        for cut in &mut self.cuts {
            // GammaLoop's physical left side contains the sewn outgoing endpoint.
            let mut left = cut
                .left
                .half_edges
                .iter()
                .map(native)
                .collect::<BTreeSet<_>>();
            let mut right = cut
                .right
                .half_edges
                .iter()
                .map(native)
                .collect::<BTreeSet<_>>();
            let reversed = !initial_left.is_subset(&left) && initial_left.is_subset(&right);
            if reversed {
                std::mem::swap(&mut left, &mut right);
                std::mem::swap(
                    &mut cut.left.coupling_orders,
                    &mut cut.right.coupling_orders,
                );
                std::mem::swap(&mut cut.left.loop_count, &mut cut.right.loop_count);
            }
            cut.cut = cut
                .cut
                .iter()
                .map(native)
                .map(|mut h| {
                    if reversed {
                        h.endpoint = match h.endpoint {
                            DiagramEndpoint::Source => DiagramEndpoint::Target,
                            DiagramEndpoint::Target => DiagramEndpoint::Source,
                        };
                    }
                    h
                })
                .collect();
            cut.left.half_edges = left.difference(&initial_right).copied().collect();
            cut.right.half_edges = right.difference(&initial_left).copied().collect();
            cut.cut.sort();
            cut.cut.dedup();
        }
        for candidate in &mut self.topology_threshold_candidates {
            let mut left = candidate.left.iter().map(native).collect::<BTreeSet<_>>();
            let mut right = candidate.right.iter().map(native).collect::<BTreeSet<_>>();
            let reversed = !initial_left.is_subset(&left) && initial_left.is_subset(&right);
            if reversed {
                std::mem::swap(&mut left, &mut right);
            }
            candidate.cut = candidate
                .cut
                .iter()
                .map(native)
                .map(|mut h| {
                    if reversed {
                        h.endpoint = match h.endpoint {
                            DiagramEndpoint::Source => DiagramEndpoint::Target,
                            DiagramEndpoint::Target => DiagramEndpoint::Source,
                        };
                    }
                    h
                })
                .filter(|h| !initial_right.contains(h) && !initial_left.contains(h))
                .collect();
            candidate.left = left.difference(&initial_right).copied().collect();
            candidate.right = right.difference(&initial_left).copied().collect();
            candidate.cut.sort();
            candidate.cut.dedup();
        }
        self.cuts.sort();
        self.cuts.dedup();
        self.topology_threshold_candidates.sort();
        self.topology_threshold_candidates.dedup();
        Ok(())
    }
}
