//! Finalize generated topology once, in GammaLoop's native half-edge convention.
use super::*;
use linnet::{
    half_edge::{
        EdgeAccessors,
        involution::Hedge,
        remap::HedgeRemap,
        subgraph::{ModifySubSet, OrientedCut},
        swap::Swap,
    },
    permutation::Permutation,
};
use std::cell::LazyCell;
use symbolica::atom::{AtomView, FunctionBuilder, Symbol};

/// Where finalization put each provisional edge of a generated diagram, indexed
/// by provisional edge id. Vertices keep their ids.
struct FinalFrame {
    /// The provisional endpoints.
    endpoints: Vec<EdgeEndpoints>,
    /// The final (logical) edge index.
    edges: Vec<usize>,
    /// The final half-edges at the provisional source and target; a dangling
    /// edge has one half-edge for both.
    half_edges: Vec<[Hedge; 2]>,
    /// Whether the final orientation reverses the provisional momentum.
    reversed: Vec<bool>,
    /// Internal non-dummy edges in canonical insertion order, only when that
    /// routing was requested. Logical edge ids need not be in this order.
    canonical_order: Option<Vec<EdgeId>>,
    vertices: usize,
}

impl FinalFrame {
    fn half_edge(&self, half_edge: DiagramHalfEdge) -> Option<Hedge> {
        let end = match half_edge.endpoint {
            DiagramEndpoint::Source => 0,
            DiagramEndpoint::Target => 1,
        };
        self.half_edges.get(half_edge.edge.0).map(|ends| ends[end])
    }

    /// Resolve provisional edge `id` on the finalized graph. Each provisional
    /// end takes a half-edge of the final edge at its vertex, preferring the
    /// end's flow (source or sink). The momentum is reversed when the
    /// half-edge at the provisional source, or at the target of a leg without
    /// a source, no longer has its provisional flow.
    fn resolve(&mut self, graph: &HedgeGraph<DiagramEdge, DiagramVertex>, id: usize) {
        let candidates = match graph[&EdgeIndex(self.edges[id])].1 {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                [source, sink]
            }
            HedgePair::Unpaired { hedge, .. } => [hedge; 2],
        };
        let endpoints = self.endpoints[id];
        let at = |vertex: Option<VertexId>, flow| {
            let attached = || {
                candidates.into_iter().filter(move |&hedge| {
                    vertex.is_none_or(|vertex| graph.node_id(hedge) == NodeIndex(vertex.0))
                })
            };
            attached()
                .find(|&hedge| graph.flow(hedge) == flow)
                .or_else(|| attached().next())
                .or_else(|| {
                    candidates
                        .into_iter()
                        .find(|&hedge| graph.flow(hedge) == flow)
                })
                .unwrap_or(candidates[0])
        };
        self.half_edges[id] = [
            at(endpoints.source, Flow::Source),
            at(endpoints.target, Flow::Sink),
        ];
        let (end, flow) = if endpoints.source.is_some() {
            (0, Flow::Source)
        } else {
            (1, Flow::Sink)
        };
        self.reversed[id] = graph.flow(self.half_edges[id][end]) != flow;
    }
}

/// How a provisional edge of a generated diagram enters its final graph.
#[derive(Clone, Copy)]
enum Insertion {
    /// A leg, dangling with its final process flow at its vertex.
    Leg(Flow),
    /// An internal edge, attached at both ends.
    Internal,
}

/// What the first argument of an index head, or of a momentum or external-state
/// carrier, names.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum IndexOwner {
    HalfEdge,
    Edge,
    Vertex,
    /// The momentum `Q` of an edge, whose sign follows the edge's orientation.
    Momentum,
    /// A spinor or polarization vector of an external edge.
    ExternalState,
}

impl IndexOwner {
    /// Carriers hold tensor indices in their other arguments.
    fn is_carrier(self) -> bool {
        matches!(self, Self::Momentum | Self::ExternalState)
    }
}

/// The index heads and carriers whose first argument names a half-edge, an edge or a
/// vertex.
pub(super) struct IndexOwners([(Symbol, IndexOwner); 10]);

impl IndexOwners {
    pub(super) fn new() -> Self {
        Self([
            (symbols::hedge_index(), IndexOwner::HalfEdge),
            (symbols::edge_index(), IndexOwner::Edge),
            (symbols::vertex_index(), IndexOwner::Vertex),
            (momentum_symbol(), IndexOwner::Momentum),
            (symbols::u(), IndexOwner::ExternalState),
            (symbols::ubar(), IndexOwner::ExternalState),
            (symbols::v(), IndexOwner::ExternalState),
            (symbols::vbar(), IndexOwner::ExternalState),
            (symbols::epsilon(), IndexOwner::ExternalState),
            (symbols::epsilonbar(), IndexOwner::ExternalState),
        ])
    }

    /// Whether an expression contains an index head or carrier to transport.
    fn contains(&self, atom: &Atom) -> bool {
        let mut found = false;
        atom.visitor(&mut |term| {
            found = found
                || matches!(term, AtomView::Fun(function)
                    if self.0.iter().any(|(head, _)| *head == function.get_symbol()));
            !found
        });
        found
    }

    /// Relabel the owner of every index head and carrier in `atom`. `relabel` returns
    /// the new owner and whether the carrier changes sign, or `None` for an unknown
    /// owner, which is an error. A zero local label of an index head is dropped:
    /// `hedge(h, 0)` becomes `hedge(h')`.
    pub(super) fn relabel(
        &self,
        atom: &Atom,
        mut relabel: impl FnMut(IndexOwner, usize) -> Option<(usize, bool)>,
    ) -> Result<Atom, String> {
        let mut failure = None;
        let mut relabelled = atom.clone();
        // Relabel the indices inside carrier arguments before their carriers. A
        // bottom-up signed replacement can collapse -Q to Q and then visit that new
        // carrier again, applying the relabelling twice.
        for carriers in [false, true] {
            relabelled = relabelled.replace_map(|term, _, out| {
                let AtomView::Fun(function) = term else {
                    return;
                };
                let symbol = function.get_symbol();
                let Some(owner) = self.0.iter().find_map(|&(head, owner)| {
                    (head == symbol && owner.is_carrier() == carriers).then_some(owner)
                }) else {
                    return;
                };
                let mut arguments = function.iter();
                let Some(index) = arguments
                    .next()
                    .and_then(|index| usize::try_from(index).ok())
                else {
                    return;
                };
                let Some((index, negated)) = relabel(owner, index) else {
                    failure = Some(format!("{symbol} refers to unknown index {index}"));
                    return;
                };
                let local = arguments.collect::<Vec<_>>();
                let local = match local.as_slice() {
                    [label] if !carriers && label.is_zero() => &[],
                    local => local,
                };
                let value = FunctionBuilder::new(symbol)
                    .add_arg(index)
                    .add_args(local.iter().copied())
                    .finish();
                **out = if negated { -value } else { value };
            });
        }
        failure.map_or(Ok(relabelled), Err)
    }
}

impl DiagramHalfEdge {
    /// The index of a builder half-edge in the provisional expressions of a generated
    /// diagram: twice its edge, plus one at its target.
    pub(super) fn provisional_index(self) -> usize {
        2 * self.edge.0 + usize::from(self.endpoint == DiagramEndpoint::Target)
    }

    fn from_provisional_index(index: usize) -> Self {
        Self {
            edge: EdgeId(index / 2),
            endpoint: if index.is_multiple_of(2) {
                DiagramEndpoint::Source
            } else {
                DiagramEndpoint::Target
            },
        }
    }
}

impl FeynmanDiagramBuilder {
    /// Finalize the diagram and return, for every builder edge, its finalized edge.
    /// Both legs of a sewn initial state map to its one edge.
    pub(super) fn finalize(self) -> Result<(FeynmanDiagram, Vec<EdgeId>), DiagramError> {
        // The vertices and edges move into the graph. Expressions, cuts and a
        // requested basis of a generated diagram refer to its provisional
        // frame; an imported diagram is already final.
        let (mut graph, frame) = if self.generation_externals.is_empty() {
            let graph = imported_graph(
                &self.model,
                self.vertices,
                self.edges,
                self.edge_orientations.as_deref(),
                self.half_edge_order.as_deref(),
            )?;
            (graph, None)
        } else if self.edge_orientations.is_some() || self.half_edge_order.is_some() {
            return Err(DiagramError::Invariant {
                operation: "finalizing generated legs",
                message: "only imported diagrams restore a serialized frame".into(),
            });
        } else {
            let (graph, frame) = generated_graph(
                &self.model,
                &self.generation_externals,
                self.vertices,
                self.edges,
                self.canonical_loop_momentum_basis && self.loop_momentum_basis.is_none(),
            )?;
            (graph, Some(frame))
        };
        // Resolve heads once per generated diagram. Imported diagrams bypass
        // transport and should not initialize its symbols.
        let owners = LazyCell::new(IndexOwners::new);
        let translate = |atom: &mut Atom| -> Result<(), DiagramError> {
            // Imported diagrams keep their frame, and atoms without an index
            // head, such as unset fragments and the generated overall factor,
            // have nothing to remap.
            let Some(frame) = &frame else {
                return Ok(());
            };
            if !owners.contains(atom) {
                return Ok(());
            }
            *atom = owners
                .relabel(atom, |owner, index| match owner {
                    IndexOwner::HalfEdge => frame
                        .half_edge(DiagramHalfEdge::from_provisional_index(index))
                        .map(|half_edge| (half_edge.0, false)),
                    IndexOwner::Vertex => (index < frame.vertices).then_some((index, false)),
                    IndexOwner::Momentum => frame
                        .edges
                        .get(index)
                        .map(|edge| (*edge, frame.reversed[index])),
                    IndexOwner::Edge | IndexOwner::ExternalState => {
                        frame.edges.get(index).map(|edge| (*edge, false))
                    }
                })
                .map_err(|message| DiagramError::Invariant {
                    operation: "finalizing symbolic indices",
                    message,
                })?;
            Ok(())
        };
        for (_, _, vertex) in graph.iter_nodes_mut() {
            translate(&mut vertex.numerator)?;
        }
        for edge in 0..graph.n_edges() {
            translate(&mut graph[EdgeIndex(edge)].numerator)?;
        }
        let numerator = match self.numerator {
            Some(mut numerator) => {
                translate(&mut numerator)?;
                DiagramNumerator::Explicit(numerator)
            }
            None => DiagramNumerator::default(),
        };
        let mut diagram = FeynmanDiagram {
            model: self.model,
            id: OnceCell::new(),
            name: self.name,
            graph,
            symmetry_factor: self.symmetry_factor,
            overall_factor: self.overall_factor,
            numerator,
            numerator_prefactor: self.numerator_prefactor,
            projector: self.projector,
            loop_momentum_basis: LoopMomentumBasis {
                tree_edges: vec![],
                loop_edges: vec![],
                external_edges: vec![],
                dependent_externals: vec![],
                edge_signatures: BTreeMap::new(),
            },
            cuts: self.cuts,
            topology_threshold_candidates: self.topology_threshold_candidates,
            validation: OnceLock::new(),
            routing_validation: OnceLock::new(),
        };
        for atom in [
            &mut diagram.overall_factor,
            &mut diagram.numerator_prefactor,
            &mut diagram.projector,
        ] {
            translate(atom)?;
        }
        // Only cross sections have cuts and threshold candidates.
        if let Some(frame) = &frame
            && !(diagram.cuts.is_empty() && diagram.topology_threshold_candidates.is_empty())
        {
            let mut incoming = BTreeSet::new();
            let mut outgoing = BTreeSet::new();
            for (edge, external) in &self.generation_externals {
                // The provisional end of the leg away from its attachment vertex.
                let half_edge = DiagramHalfEdge {
                    edge: *edge,
                    endpoint: if frame.endpoints[edge.0].source.is_none() {
                        DiagramEndpoint::Source
                    } else {
                        DiagramEndpoint::Target
                    },
                };
                match external.state {
                    ExternalState::Incoming => incoming.insert(half_edge),
                    ExternalState::Outgoing => outgoing.insert(half_edge),
                };
            }
            diagram.finalize_cut_indices(frame, &incoming, &outgoing)?;
        }
        diagram = if let Some(basis) = self.loop_momentum_basis {
            if let Some(frame) = &frame {
                let requested = basis
                    .loop_edges
                    .iter()
                    .map(|id| {
                        frame.edges.get(id.0).copied().map(EdgeId).ok_or(
                            DiagramError::UnknownEdge {
                                edge: id.0,
                                edges: frame.edges.len(),
                            },
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                diagram.with_loop_momentum_edges(&requested)?
            } else {
                basis.validate(&diagram)?;
                diagram.loop_momentum_basis = basis;
                diagram
            }
        } else if self.canonical_loop_momentum_basis {
            if let Some(frame) = &frame {
                let ordered = frame
                    .canonical_order
                    .as_deref()
                    .expect("canonical generated routing collected its order");
                diagram.with_canonical_loop_momentum_basis_finalized(ordered)?
            } else {
                diagram.with_canonical_loop_momentum_basis()?
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
        // Each branch above either constructs a routing from a spanning forest
        // or validates an imported routing against that same construction.
        let _ = diagram.routing_validation.set(());
        // Project threshold candidates outside the initial-state attachment tree.
        if frame.is_some() && !diagram.topology_threshold_candidates.is_empty() {
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
            if !excluded.is_empty() {
                let vertex_half_edges = diagram.vertex_half_edges();
                let mut candidates = std::mem::take(&mut diagram.topology_threshold_candidates);
                for (index, candidate) in candidates.iter_mut().enumerate() {
                    // The loop-independent initial-state attachment tree belongs
                    // to neither threshold side, as in GammaLoop finalization.
                    candidate
                        .left
                        .retain(|half_edge| !excluded.contains(half_edge));
                    candidate
                        .right
                        .retain(|half_edge| !excluded.contains(half_edge));
                    // A fixed-momentum crossing into an omitted crown disappears,
                    // while loop-carrying crossings retain their oriented boundary.
                    candidate.cut = diagram
                        .partition_crossing(
                            index,
                            &candidate.left,
                            &candidate.right,
                            &vertex_half_edges,
                            Some(&ignored),
                        )
                        .map_err(|error| match error {
                            DiagramError::InvalidCut { message, .. } => {
                                DiagramError::InvalidThresholdCandidate {
                                    candidate: index,
                                    message,
                                }
                            }
                            error => error,
                        })?
                        .crossing
                        .into_iter()
                        .collect();
                }
                diagram.topology_threshold_candidates = candidates;
            }
        }
        let mut cuts = std::mem::take(&mut diagram.cuts);
        cuts.sort_by_cached_key(|cut| diagram.cut_order_key(&cut.cut));
        diagram.cuts = cuts;
        let mut thresholds = std::mem::take(&mut diagram.topology_threshold_candidates);
        thresholds.sort_by_cached_key(|candidate| diagram.cut_order_key(&candidate.cut));
        diagram.topology_threshold_candidates = thresholds;
        let edges = match frame {
            Some(frame) => frame.edges.into_iter().map(EdgeId).collect(),
            None => (0..diagram.graph.n_edges()).map(EdgeId).collect(),
        };
        Ok((diagram, edges))
    }
}

/// The orientation of `edge`: undirected, or against an antiparticle's edge so
/// that directed edges follow the particle flow.
fn orientation(model: &Model, edge: &DiagramEdge) -> Result<Orientation, DiagramError> {
    let particle = model.particle_by_id(edge.particle)?;
    Ok(if !edge.directed {
        Orientation::Undirected
    } else if particle.is_antiparticle() {
        Orientation::Reversed
    } else {
        Orientation::Default
    })
}

/// Build the graph of a generated diagram once, in its final frame, with the
/// payloads moved in, and record where each provisional edge went.
///
/// Half-edge insertion order controls channel coordinates: the legs by
/// connection (the incoming legs, then the outgoing legs of connections
/// without one), the internal edges in canonical order, and last the outgoing
/// copies of cross-section initial states, which sewing merges into their
/// incoming copies. Logical edge identity is assigned separately, as in
/// GammaLoop's finalization: the legs as inserted, then the internal edges in
/// provisional order. The frame follows from the insertion order; only sewn
/// legs are resolved on the graph, and every edge when sewing moves
/// half-edges.
fn generated_graph(
    model: &Model,
    legs: &BTreeMap<EdgeId, ExternalLeg>,
    vertices: Vec<DiagramVertex>,
    edges: Vec<(EdgeEndpoints, DiagramEdge)>,
    canonical_loop_momentum_basis: bool,
) -> Result<(HedgeGraph<DiagramEdge, DiagramVertex>, FinalFrame), DiagramError> {
    let mut external_indices = BTreeSet::new();
    for leg in legs.values() {
        if !external_indices.insert(leg.index) {
            return Err(DiagramError::DuplicateExternalIndex(leg.index));
        }
    }
    // The incoming and outgoing leg of each connection.
    let mut connections = BTreeMap::<usize, (Option<usize>, Option<usize>)>::new();
    let mut internal = Vec::new();
    for (id, (endpoints, edge)) in edges.iter().enumerate() {
        model.particle_by_id(edge.particle)?;
        let Some(leg) = legs.get(&EdgeId(id)) else {
            if endpoints.source.is_none() || endpoints.target.is_none() || edge.external.is_some() {
                return Err(DiagramError::Invariant {
                    operation: "finalizing generated legs",
                    message: format!(
                        "edge {id} is dangling or carries a leg but is not a generated leg"
                    ),
                });
            }
            internal.push(id);
            continue;
        };
        if edge.numerator != *ONE {
            return Err(DiagramError::ExternalEdgeNumerator { edge: id });
        }
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
    }
    // Every provisional edge in insertion order.
    let mut order = connections
        .values()
        .filter_map(|&(incoming, _)| Some((incoming?, Insertion::Leg(Flow::Sink))))
        .collect::<Vec<_>>();
    order.extend(connections.values().filter_map(|&(incoming, outgoing)| {
        Some((
            outgoing.filter(|_| incoming.is_none())?,
            Insertion::Leg(Flow::Source),
        ))
    }));
    let n_edges = edges.len();
    let mut final_edges = vec![0; n_edges];
    for (index, &id) in order.iter().map(|(id, _)| id).chain(&internal).enumerate() {
        final_edges[id] = index;
    }
    internal.sort_by_cached_key(|&id| {
        let (endpoints, edge) = &edges[id];
        let (source, target, particle) = if endpoints.source <= endpoints.target {
            (endpoints.source, endpoints.target, edge.particle)
        } else {
            (
                endpoints.target,
                endpoints.source,
                model
                    .particle_by_id(edge.particle)
                    .expect("resolved particle")
                    .antiparticle,
            )
        };
        (
            source,
            target,
            model
                .particle_by_id(particle)
                .expect("resolved particle")
                .pdg_code,
            id,
        )
    });
    // Edge indices follow the insertion order unless the provisional order of
    // the internal edges is not canonical.
    let relabel = !internal.is_sorted();
    let canonical_order = canonical_loop_momentum_basis.then(|| Vec::with_capacity(internal.len()));
    order.extend(internal.into_iter().map(|id| (id, Insertion::Internal)));
    let logical_edges = order.len();
    let sewn = connections
        .into_values()
        .filter_map(|(incoming, outgoing)| Some((incoming?, outgoing?)))
        .collect::<Vec<_>>();
    for &(incoming, outgoing) in &sewn {
        final_edges[outgoing] = final_edges[incoming];
        order.push((outgoing, Insertion::Leg(Flow::Source)));
    }
    let mut builder = HedgeGraphBuilder::with_capacity(vertices.len(), 2 * n_edges);
    let mut frame = FinalFrame {
        endpoints: Vec::with_capacity(n_edges),
        edges: final_edges,
        half_edges: vec![[Hedge(0); 2]; n_edges],
        reversed: vec![false; n_edges],
        canonical_order,
        vertices: vertices.len(),
    };
    for vertex in vertices {
        if let Some(interaction) = vertex.interaction {
            model.vertex_rule_by_id(interaction)?;
        }
        builder.add_node(vertex);
    }
    let mut payloads = Vec::with_capacity(n_edges);
    for (endpoints, edge) in edges {
        frame.endpoints.push(endpoints);
        payloads.push(Some(edge));
    }
    let mut next_hedge = 0;
    for &(id, insertion) in &order {
        let mut edge = payloads[id].take().expect("edges are inserted once");
        let endpoints = frame.endpoints[id];
        let Insertion::Leg(flow) = insertion else {
            let (Some(source), Some(target)) = (endpoints.source, endpoints.target) else {
                unreachable!("internal edges of generated diagrams are attached at both ends")
            };
            let direction = orientation(model, &edge)?;
            if !edge.is_dummy
                && let Some(ordered) = &mut frame.canonical_order
            {
                // Internal insertion order is canonical, while logical ids
                // retain provisional order even after the graph is relabelled.
                ordered.push(EdgeId(frame.edges[id]));
            }
            builder.add_edge(NodeIndex(source.0), NodeIndex(target.0), edge, direction);
            // The source, then the sink.
            frame.half_edges[id] = [Hedge(next_hedge), Hedge(next_hedge + 1)];
            next_hedge += 2;
            continue;
        };
        // The flow and slot of the leg at its attachment vertex.
        let (vertex, old_flow, slot) = match (endpoints.source, endpoints.target) {
            (Some(vertex), None) => (vertex, Flow::Source, edge.source_slot),
            (None, Some(vertex)) => (vertex, Flow::Sink, edge.target_slot),
            _ => unreachable!("generated legs are dangling edges"),
        };
        edge.external = Some(legs[&EdgeId(id)].clone());
        if old_flow != flow {
            edge.particle = model.particle_by_id(edge.particle)?.antiparticle;
        }
        match flow {
            Flow::Source => edge.source_slot = slot,
            Flow::Sink => edge.target_slot = slot,
        }
        let direction = orientation(model, &edge)?;
        builder.add_external_edge(NodeIndex(vertex.0), edge, direction, flow);
        frame.half_edges[id] = [Hedge(next_hedge); 2];
        frame.reversed[id] = flow != old_flow;
        next_hedge += 1;
    }
    let mut graph: HedgeGraph<DiagramEdge, DiagramVertex> = builder.into();
    // Without a sewn connection, no two legs match.
    let mut moved_hedges = false;
    if !sewn.is_empty() {
        graph
            .sew(
                |_, left, _, right| match (&left.data.external, &right.data.external) {
                    (Some(left), Some(right)) => {
                        left.connection == right.connection && left.state != right.state
                    }
                    _ => false,
                },
                |left_flow, mut left, right_flow, right| match (left_flow, right_flow) {
                    (Flow::Sink, Flow::Source) => {
                        left.data.source_slot = right.data.source_slot;
                        (Flow::Sink, left)
                    }
                    (Flow::Source, Flow::Sink) => {
                        let mut incoming = right;
                        incoming.data.source_slot = left.data.source_slot;
                        (Flow::Source, incoming)
                    }
                    _ => unreachable!("external connections have opposite process flows"),
                },
            )
            .map_err(|error| DiagramError::Invariant {
                operation: "sewing initial-state edges",
                message: format!("{error:?}"),
            })?;
        // The sinks of the sewn edges come first, by connection. They do
        // already unless an incoming leg without an outgoing copy precedes a
        // sewn one.
        let cut_hedges = graph
            .iter_edges()
            .filter_map(|(pair, _, edge)| {
                if let HedgePair::Paired { sink, .. } = pair {
                    edge.data
                        .external
                        .as_ref()
                        .map(|external| (external.connection, sink))
                } else {
                    None
                }
            })
            .collect::<BTreeMap<_, _>>();
        moved_hedges = cut_hedges
            .values()
            .enumerate()
            .any(|(target, hedge)| hedge.0 != target);
        if moved_hedges {
            let permutation = Permutation::from_mappings(
                cut_hedges
                    .values()
                    .enumerate()
                    .map(|(target, hedge)| (hedge.0, target)),
                graph.n_hedges(),
            )
            .map_err(|error| DiagramError::Invariant {
                operation: "ordering initial-state half-edges",
                message: error.to_string(),
            })?;
            graph
                .permute_hedges(&HedgeRemap::from_permutation(&permutation))
                .map_err(|error| DiagramError::Invariant {
                    operation: "ordering initial-state half-edges",
                    message: error.to_string(),
                })?;
        }
    }
    // Sewing keeps the edge index of the incoming copy and drops the outgoing
    // copies, which come last.
    if relabel {
        let permutation = Permutation::from_mappings(
            order[..logical_edges]
                .iter()
                .enumerate()
                .map(|(index, &(id, _))| (index, frame.edges[id])),
            graph.n_edges(),
        )
        .map_err(|error| DiagramError::Invariant {
            operation: "ordering logical edges",
            message: error.to_string(),
        })?;
        <HedgeGraph<_, _> as Swap<EdgeIndex>>::permute(&mut graph, &permutation);
    }
    if moved_hedges {
        (0..n_edges).for_each(|id| frame.resolve(&graph, id));
    } else {
        for (incoming, outgoing) in sewn {
            frame.resolve(&graph, incoming);
            frame.resolve(&graph, outgoing);
        }
    }
    Ok((graph, frame))
}

/// Build the graph of an imported diagram in its serialized edge order and
/// restore its serialized orientations and half-edge order.
fn imported_graph(
    model: &Model,
    vertices: Vec<DiagramVertex>,
    edges: Vec<(EdgeEndpoints, DiagramEdge)>,
    orientations: Option<&[Orientation]>,
    half_edge_order: Option<&[DiagramHalfEdge]>,
) -> Result<HedgeGraph<DiagramEdge, DiagramVertex>, DiagramError> {
    let mut builder = HedgeGraphBuilder::with_capacity(vertices.len(), 2 * edges.len());
    for vertex in vertices {
        if let Some(interaction) = vertex.interaction {
            model.vertex_rule_by_id(interaction)?;
        }
        builder.add_node(vertex);
    }
    for (endpoints, edge) in edges {
        let direction = orientation(model, &edge)?;
        match (endpoints.source, endpoints.target) {
            (Some(source), Some(target)) => {
                builder.add_edge(NodeIndex(source.0), NodeIndex(target.0), edge, direction);
            }
            (Some(vertex), None) => {
                builder.add_external_edge(NodeIndex(vertex.0), edge, direction, Flow::Source);
            }
            (None, Some(vertex)) => {
                builder.add_external_edge(NodeIndex(vertex.0), edge, direction, Flow::Sink);
            }
            (None, None) => unreachable!("builder validates attached endpoints"),
        }
    }
    let mut graph: HedgeGraph<DiagramEdge, DiagramVertex> = builder.into();
    if let Some(orientations) = orientations {
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
    if let Some(order) = half_edge_order {
        let mut mappings = Vec::new();
        for (target, half_edge) in order.iter().enumerate() {
            let Some(pair) =
                (half_edge.edge.0 < graph.n_edges()).then(|| graph[&EdgeIndex(half_edge.edge.0)].1)
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
        graph
            .permute_hedges(&HedgeRemap::from_permutation(&permutation))
            .map_err(|error| DiagramError::Invariant {
                operation: "restoring half-edge order",
                message: error.to_string(),
            })?;
    }
    Ok(graph)
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
        frame: &FinalFrame,
        incoming: &BTreeSet<DiagramHalfEdge>,
        outgoing: &BTreeSet<DiagramHalfEdge>,
    ) -> Result<(), DiagramError> {
        // Decide the process side before sewing: the two provisional external
        // edges become aliases, so their original membership is not recoverable
        // from the sewn half-edge sets. GammaLoop retains the incoming sink on
        // the physical left side and the positive external momentum frame.
        let requires_swap = |left: &[DiagramHalfEdge],
                             right: &[DiagramHalfEdge]|
         -> Result<bool, DiagramError> {
            let left = left.iter().copied().collect::<BTreeSet<_>>();
            let right = right.iter().copied().collect::<BTreeSet<_>>();
            match (incoming.is_subset(&left) && outgoing.is_subset(&right), incoming.is_subset(&right) && outgoing.is_subset(&left)) {
                (true, false) => Ok(false),
                (false, true) => Ok(true),
                _ => Err(DiagramError::Invariant { operation: "normalizing cut sides", message: "the provisional partition does not have an unambiguous incoming/outgoing side".into() }),
            }
        };
        let native = |half_edge: &DiagramHalfEdge| -> Result<DiagramHalfEdge, DiagramError> {
            let hedge = frame
                .half_edge(*half_edge)
                .ok_or_else(|| DiagramError::Invariant {
                    operation: "finalizing cut endpoints",
                    message: format!("unknown provisional endpoint {half_edge:?}"),
                })?;
            Ok(DiagramHalfEdge {
                edge: EdgeId(self.graph[&hedge].0),
                endpoint: if self.graph.flow(hedge) == Flow::Source {
                    DiagramEndpoint::Source
                } else {
                    DiagramEndpoint::Target
                },
            })
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
            // GammaLoop's physical left side contains the sewn incoming endpoint.
            let reversed = requires_swap(&cut.left.half_edges, &cut.right.half_edges)?;
            let mut left = cut
                .left
                .half_edges
                .iter()
                .map(native)
                .collect::<Result<BTreeSet<_>, _>>()?;
            let mut right = cut
                .right
                .half_edges
                .iter()
                .map(native)
                .collect::<Result<BTreeSet<_>, _>>()?;
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
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
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
            let reversed = requires_swap(&candidate.left, &candidate.right)?;
            let mut left = candidate
                .left
                .iter()
                .map(native)
                .collect::<Result<BTreeSet<_>, _>>()?;
            let mut right = candidate
                .right
                .iter()
                .map(native)
                .collect::<Result<BTreeSet<_>, _>>()?;
            if reversed {
                std::mem::swap(&mut left, &mut right);
            }
            candidate.cut = candidate
                .cut
                .iter()
                .map(native)
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tests::scalar_model;

    /// An edge of a generated diagram, in provisional order.
    enum Provisional {
        /// A leg at a vertex where its end has the flow, with its particle,
        /// process state and connection.
        Leg(usize, Flow, &'static str, ExternalState, usize),
        /// An internal edge from a vertex to a vertex, with its particle.
        Internal(usize, usize, &'static str),
    }
    use ExternalState::{Incoming as In, Outgoing as Out};
    use Flow::{Sink, Source};
    use Provisional::{Internal, Leg};

    fn generated(
        model: Arc<Model>,
        vertices: usize,
        edges: &[Provisional],
    ) -> FeynmanDiagramBuilder {
        let rule = model.vertex_rule_id_at(0).unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "final-frame");
        let vertices = (0..vertices)
            .map(|vertex| {
                builder.add_vertex(DiagramVertex::interaction(format!("v{vertex}"), rule))
            })
            .collect::<Vec<_>>();
        let edge = |particle| {
            let particle = model.particle_id(particle).unwrap();
            DiagramEdge::new(particle, model.particle_by_id(particle).unwrap().spin == 2)
        };
        for (index, provisional) in edges.iter().enumerate() {
            match *provisional {
                Leg(vertex, flow, particle, state, connection) => {
                    let leg = ExternalLeg {
                        name: format!("p{index}"),
                        index,
                        state,
                        connection,
                    };
                    builder
                        .add_generation_external(
                            vertices[vertex],
                            flow,
                            edge(particle),
                            VertexSlot(0),
                            leg,
                        )
                        .unwrap();
                }
                Internal(source, target, particle) => {
                    builder
                        .add_edge(vertices[source], vertices[target], edge(particle))
                        .unwrap();
                }
            }
        }
        builder
    }

    #[test]
    fn generated_build_rejects_unknown_model_references_without_demanding_id() {
        for directed in [false, true] {
            for edge in [0, 1] {
                let mut builder = case("relabelled").canonical_loop_momentum_basis();
                builder.edges[edge].1.particle = ParticleId::from_index(usize::MAX);
                builder.edges[edge].1.directed = directed;
                assert!(matches!(
                    builder.build(),
                    Err(DiagramError::Model(ModelError::NotFound { .. }))
                ));
            }
        }
        let mut builder = case("relabelled");
        builder.vertices[0].interaction = Some(VertexRuleId::from_index(usize::MAX));
        assert!(matches!(
            builder.build(),
            Err(DiagramError::Model(ModelError::NotFound { .. }))
        ));
    }

    /// Generated diagrams that exercise every path of the final-frame build.
    fn cases() -> Vec<(&'static str, FeynmanDiagramBuilder)> {
        let qed = || Arc::new(Model::qed());
        vec![
            // Internal edges out of canonical order: edge indices are relabelled.
            (
                "relabelled",
                generated(
                    scalar_model(),
                    4,
                    &[
                        Leg(0, Sink, "phi", In, 0),
                        Internal(0, 1, "phi"),
                        Internal(1, 2, "phi"),
                        Internal(1, 2, "phi"),
                        Internal(2, 3, "phi"),
                        Internal(0, 3, "phi"),
                        Leg(3, Source, "phi", Out, 1),
                    ],
                ),
            ),
            // Two components whose external-root priority opposes canonical
            // internal order, with interleaved and reversed provisional edges.
            (
                "disconnected",
                generated(
                    scalar_model(),
                    4,
                    &[
                        Leg(2, Sink, "phi", In, 0),
                        Internal(2, 3, "phi"),
                        Leg(0, Sink, "phi", In, 1),
                        Internal(0, 1, "phi"),
                        Internal(3, 2, "phi"),
                        Leg(3, Source, "phi", Out, 3),
                        Internal(1, 0, "phi"),
                        Leg(1, Source, "phi", Out, 2),
                    ],
                ),
            ),
            // A self-loop before its bridge.
            (
                "tadpole",
                generated(
                    scalar_model(),
                    2,
                    &[
                        Leg(0, Sink, "phi", In, 0),
                        Internal(1, 1, "phi"),
                        Internal(1, 0, "phi"),
                        Leg(0, Source, "phi", Out, 1),
                    ],
                ),
            ),
            // Legs out of connection order, the incoming one against the process flow.
            (
                "reversed_legs",
                generated(
                    scalar_model(),
                    1,
                    &[
                        Leg(0, Source, "phi", Out, 2),
                        Leg(0, Source, "phi", In, 0),
                        Leg(0, Source, "phi", Out, 1),
                    ],
                ),
            ),
            // A cross section whose initial states are all sewn, one against the process flow.
            (
                "sewn",
                generated(
                    scalar_model(),
                    4,
                    &[
                        Leg(0, Sink, "phi", In, 0),
                        Leg(1, Source, "phi", In, 1),
                        Leg(2, Source, "phi", Out, 0),
                        Leg(3, Source, "phi", Out, 1),
                        Internal(0, 1, "phi"),
                        Internal(1, 2, "phi"),
                        Internal(2, 3, "phi"),
                        Internal(3, 0, "phi"),
                    ],
                ),
            ),
            // An incoming leg without an outgoing copy, against the process
            // flow, before a sewn one: the sewn sink moves to the front, and the
            // internal edges are relabelled.
            (
                "partially_sewn",
                generated(
                    scalar_model(),
                    4,
                    &[
                        Leg(0, Source, "phi", In, 0),
                        Leg(1, Sink, "phi", In, 1),
                        Leg(2, Source, "phi", Out, 1),
                        Leg(3, Source, "phi", Out, 2),
                        Internal(2, 3, "phi"),
                        Internal(0, 1, "phi"),
                        Internal(1, 2, "phi"),
                        Internal(3, 0, "phi"),
                        Internal(1, 3, "phi"),
                    ],
                ),
            ),
            // A cross section with both initial states at one vertex, and a
            // bubble between two loop-independent bridges.
            (
                "s_channel",
                generated(
                    scalar_model(),
                    4,
                    &[
                        Leg(0, Sink, "phi", In, 0),
                        Leg(0, Sink, "phi", In, 1),
                        Leg(3, Source, "phi", Out, 0),
                        Leg(3, Source, "phi", Out, 1),
                        Internal(0, 1, "phi"),
                        Internal(1, 2, "phi"),
                        Internal(1, 2, "phi"),
                        Internal(2, 3, "phi"),
                    ],
                ),
            ),
            // Fermion legs against the process flow and an s-channel photon.
            (
                "fermion_legs",
                generated(
                    qed(),
                    2,
                    &[
                        Leg(1, Sink, "mu+", Out, 2),
                        Leg(0, Source, "e+", In, 0),
                        Internal(1, 0, "a"),
                        Leg(0, Sink, "e+", In, 1),
                        Leg(1, Source, "mu+", Out, 3),
                    ],
                ),
            ),
            // Antiparticle and anti-parallel fermion edges.
            (
                "fermion_loop",
                generated(
                    qed(),
                    2,
                    &[
                        Leg(0, Sink, "a", In, 0),
                        Internal(1, 0, "e+"),
                        Internal(0, 1, "e+"),
                        Leg(1, Source, "a", Out, 1),
                    ],
                ),
            ),
        ]
    }

    fn case(name: &str) -> FeynmanDiagramBuilder {
        cases()
            .into_iter()
            .find_map(|(case, builder)| (case == name).then_some(builder))
            .unwrap()
    }

    /// The final graph and frame of a generated builder.
    fn final_graph(
        builder: &FeynmanDiagramBuilder,
    ) -> (HedgeGraph<DiagramEdge, DiagramVertex>, FinalFrame) {
        generated_graph(
            &builder.model,
            &builder.generation_externals,
            builder.vertices.clone(),
            builder.edges.clone(),
            builder.canonical_loop_momentum_basis && builder.loop_momentum_basis.is_none(),
        )
        .unwrap()
    }

    /// The graph of `diagram` rebuilt by the import path from its edges,
    /// orientations and half-edge order, as JSON import restores it.
    fn reimported(diagram: &FeynmanDiagram) -> HedgeGraph<DiagramEdge, DiagramVertex> {
        let orientations = (0..diagram.graph.n_edges())
            .map(|edge| diagram.graph.orientation(EdgeIndex(edge)))
            .collect::<Vec<_>>();
        let order = diagram.half_edges().collect::<Vec<_>>();
        imported_graph(
            &diagram.model,
            diagram
                .vertices()
                .map(|(_, vertex)| vertex.clone())
                .collect(),
            diagram
                .edges()
                .map(|(_, endpoints, edge)| (endpoints, edge.clone()))
                .collect(),
            Some(&orientations),
            Some(&order),
        )
        .unwrap()
    }

    #[test]
    fn final_frame_matches_the_frame_resolved_on_the_graph() {
        for (name, builder) in cases() {
            let (graph, frame) = final_graph(&builder);
            let edges = frame.edges.len();
            let mut resolved = FinalFrame {
                endpoints: frame.endpoints.clone(),
                edges: frame.edges.clone(),
                half_edges: vec![[Hedge(0); 2]; edges],
                reversed: vec![false; edges],
                canonical_order: None,
                vertices: frame.vertices,
            };
            for id in 0..edges {
                resolved.resolve(&graph, id);
                // Each provisional end, as cuts name it, lies on its final edge
                // at its vertex.
                let endpoints = frame.endpoints[id];
                for (endpoint, vertex) in [
                    (DiagramEndpoint::Source, endpoints.source),
                    (DiagramEndpoint::Target, endpoints.target),
                ] {
                    let hedge = frame
                        .half_edge(DiagramHalfEdge {
                            edge: EdgeId(id),
                            endpoint,
                        })
                        .unwrap();
                    assert_eq!(graph[&hedge], EdgeIndex(frame.edges[id]), "{name}");
                    if let Some(vertex) = vertex {
                        assert_eq!(graph.node_id(hedge), NodeIndex(vertex.0), "{name}");
                    }
                }
            }
            assert_eq!(frame.half_edges, resolved.half_edges, "{name}");
            assert_eq!(frame.reversed, resolved.reversed, "{name}");
        }
    }

    #[test]
    fn final_frame_relabels_edges_and_moves_cut_half_edges_only_when_needed() {
        // The second internal edge inserted is the last in provisional order.
        let (graph, frame) = final_graph(&case("relabelled"));
        assert_eq!(frame.half_edges[5], [Hedge(4), Hedge(5)]);
        assert_eq!(graph[&Hedge(4)], EdgeIndex(6));
        // The sewn incoming leg was inserted after the unsewn one.
        let (graph, frame) = final_graph(&case("partially_sewn"));
        let sinks = graph
            .iter_edges()
            .filter_map(|(pair, _, edge)| match pair {
                HedgePair::Paired { sink, .. } if edge.data.external.is_some() => Some(sink),
                _ => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(sinks, [Hedge(0)]);
        assert_eq!(frame.half_edges[1][1], Hedge(0));
        assert_ne!(frame.half_edges[0], [Hedge(0); 2]);
        // Fully sewn initial states keep their half-edges.
        let (_, frame) = final_graph(&case("sewn"));
        assert_eq!(frame.half_edges[0][1], Hedge(0));
        assert_eq!(frame.half_edges[1][1], Hedge(1));
    }

    #[test]
    fn generated_graphs_equal_their_import_reconstruction() {
        for canonical in [false, true] {
            for (name, mut builder) in cases() {
                if canonical {
                    builder = builder.canonical_loop_momentum_basis();
                }
                let (_, frame) = final_graph(&builder);
                assert_eq!(frame.canonical_order.is_some(), canonical, "{name}");
                let diagram = builder.build().unwrap();
                assert!(diagram.id.get().is_none(), "{name}");
                assert!(
                    matches!(&diagram.numerator, DiagramNumerator::Derived(value) if value.get().is_none()),
                    "{name}"
                );
                if let Some(ordered) = frame.canonical_order {
                    assert_eq!(
                        ordered,
                        diagram.canonical_internal_edge_order().unwrap(),
                        "{name}"
                    );
                    let inserted = diagram
                        .graph
                        .iter_edges_of(&diagram.graph.full_filter())
                        .filter_map(|(_, id, edge)| {
                            (edge.data.external.is_none() && !edge.data.is_dummy)
                                .then_some(EdgeId(id.0))
                        })
                        .collect::<Vec<_>>();
                    assert_eq!(ordered, inserted, "{name}: canonical half-edge subsequence");
                    // The final graph's crown order determines the DFS stack
                    // order, including after partial sewing and relabelling.
                    for (_, crown, _) in diagram.graph.iter_nodes() {
                        assert!(crown.is_sorted(), "{name}: node crowns stay ascending");
                    }
                    let reference = diagram
                        .clone()
                        .with_canonical_loop_momentum_basis()
                        .unwrap();
                    assert_eq!(
                        diagram.loop_momentum_basis(),
                        reference.loop_momentum_basis(),
                        "{name}"
                    );
                }
                assert_eq!(
                    format!("{:?}", diagram.graph),
                    format!("{:?}", reimported(&diagram)),
                    "{name}"
                );
            }
        }
    }

    #[test]
    fn canonical_order_excludes_dummy_edges() {
        let mut builder = case("relabelled").canonical_loop_momentum_basis();
        builder.edges[1].1.is_dummy = true;
        let (_, frame) = final_graph(&builder);
        // Provisional edge 1 becomes logical edge 2. The external edges are
        // logical 0 and 1, and internal insertion still differs from edge ids.
        assert_eq!(
            frame.canonical_order.unwrap(),
            [EdgeId(6), EdgeId(3), EdgeId(4), EdgeId(5)]
        );
    }

    #[test]
    fn explicit_loop_basis_skips_canonical_order_collection() {
        let builder = case("relabelled").canonical_loop_momentum_basis();
        let (_, frame) = final_graph(&builder);
        let diagram = builder.build().unwrap();
        let expected = diagram
            .loop_momentum_bases_with_limit(2)
            .unwrap()
            .into_iter()
            .find(|basis| basis.loop_edges != diagram.loop_momentum_basis().loop_edges)
            .expect("the two-loop fixture has a noncanonical basis");
        let mut provisional = expected.clone();
        for edge in &mut provisional.loop_edges {
            edge.0 = frame.edges.iter().position(|&id| id == edge.0).unwrap();
        }
        let builder = case("relabelled")
            .canonical_loop_momentum_basis()
            .loop_momentum_basis(provisional);
        assert!(final_graph(&builder).1.canonical_order.is_none());
        assert_eq!(builder.build().unwrap().loop_momentum_basis(), &expected);
    }

    #[test]
    fn finalization_translates_fragments_in_the_final_frame() {
        // Both half-edges, the momentum and the index of every provisional edge.
        let fragment = |frame: Option<&FinalFrame>| {
            let terms = (0..9)
                .map(|id| {
                    let (sign, edge, [source, target]) = frame.map_or((1, id, [2 * id, 2 * id + 1]), |frame| {
                        (
                            if frame.reversed[id] { -1 } else { 1 },
                            frame.edges[id],
                            frame.half_edges[id].map(|hedge| hedge.0),
                        )
                    });
                    format!(
                        "{sign}*gammalooprs::Q({edge},spenso::mink(4,gammalooprs::hedge({source},1)))\
                         +{sign}*gammalooprs::Q({edge},spenso::mink(4,gammalooprs::hedge({target},2)))\
                         +gammalooprs::edge({edge},3)"
                    )
                })
                .collect::<Vec<_>>();
            Atom::parse(
                format!("({})*gammalooprs::vertex(3,4)", terms.join("+")),
                "feynkit_graph",
                ParseSettings::default(),
            )
            .unwrap()
        };
        let (_, frame) = final_graph(&case("partially_sewn"));
        let mut builder = case("partially_sewn");
        builder.vertices[0].numerator = fragment(None);
        let diagram = builder.build().unwrap();
        assert_eq!(diagram.numerator(), &fragment(Some(&frame)));
    }

    #[test]
    fn threshold_candidates_without_cuts_reach_the_final_frame() {
        // The bubble of the s-channel cross section cut, vertices 0 and 1
        // against 2 and 3, in provisional half-edges as the generator
        // partitions it: each side holds the half-edges at its vertices and
        // both ends of its legs, the cut the left ends of the bubble.
        let (source, target) = (DiagramEndpoint::Source, DiagramEndpoint::Target);
        let half_edges = |ends: &[(usize, DiagramEndpoint)]| {
            ends.iter()
                .map(|&(edge, endpoint)| DiagramHalfEdge {
                    edge: EdgeId(edge),
                    endpoint,
                })
                .collect::<Vec<_>>()
        };
        let candidate = DiagramThresholdCandidate {
            cut: half_edges(&[(5, source), (6, source)]),
            left: half_edges(&[
                (0, source),
                (0, target),
                (1, source),
                (1, target),
                (4, source),
                (4, target),
                (5, source),
                (6, source),
            ]),
            right: half_edges(&[
                (2, source),
                (2, target),
                (3, source),
                (3, target),
                (5, target),
                (6, target),
                (7, source),
                (7, target),
            ]),
        };
        let diagram = case("s_channel")
            .topology_threshold_candidates(vec![candidate.clone()])
            .build()
            .unwrap();
        let [finalized] = diagram.topology_threshold_candidates() else {
            panic!("one candidate");
        };
        let vertex = |half_edge| {
            diagram
                .graph
                .node_id(diagram.half_edge_id(half_edge).unwrap())
                .0
        };
        // The crowns of the initial-state vertices 0 and 3, attached through
        // the loop-independent bridges, belong to neither side.
        assert_eq!(finalized.cut.len(), 2);
        assert!(
            finalized
                .cut
                .iter()
                .chain(&finalized.left)
                .all(|&half_edge| vertex(half_edge) == 1)
        );
        assert!(
            finalized
                .right
                .iter()
                .all(|&half_edge| vertex(half_edge) == 2)
        );
        assert_eq!((finalized.left.len(), finalized.right.len()), (3, 3));
        // Physical cuts do not change how candidates are finalized.
        let side = |half_edges: &[DiagramHalfEdge]| DiagramCutSide {
            half_edges: half_edges.to_vec(),
            coupling_orders: BTreeMap::new(),
            loop_count: 0,
        };
        let cut = DiagramCut {
            cut: candidate.cut.clone(),
            left: side(&candidate.left),
            right: side(&candidate.right),
        };
        let with_cut = case("s_channel")
            .topology_threshold_candidates(vec![candidate])
            .cuts(vec![cut])
            .build()
            .unwrap();
        assert_eq!(
            with_cut.topology_threshold_candidates(),
            diagram.topology_threshold_candidates()
        );
    }

    #[test]
    fn generated_thresholds_project_fixed_and_loop_carrying_crossings() {
        let (source, target) = (DiagramEndpoint::Source, DiagramEndpoint::Target);
        // A bubble between attachment bridges, with a second bubble on the
        // incoming attachment crown. The second bubble has no crossing before
        // projection, but both of its edges cross the omitted crown afterward.
        let mut bridged = generated(
            scalar_model(),
            5,
            &[
                Leg(1, Sink, "phi", In, 0),
                Leg(0, Source, "phi", Out, 0),
                Internal(1, 2, "phi"),
                Internal(2, 3, "phi"),
                Internal(2, 3, "phi"),
                Internal(3, 0, "phi"),
                Internal(1, 4, "phi"),
                Internal(4, 1, "phi"),
            ],
        );
        for vertex in [0, 1, 4] {
            bridged.vertices[vertex].interaction = None;
        }
        for (mut builder, left_vertices, physical_left, expected, before, ignored_crossings) in [
            (
                // The captured phi3 two-loop failure: the 2 -> 4 bridge used
                // to remain in the boundary after crown 2 had been removed.
                generated(
                    scalar_model(),
                    6,
                    &[
                        Leg(2, Sink, "phi", In, 0),
                        Leg(1, Sink, "phi", In, 1),
                        Leg(5, Source, "phi", Out, 0),
                        Leg(5, Source, "phi", Out, 1),
                        Internal(0, 0, "phi"),
                        Internal(0, 2, "phi"),
                        Internal(1, 3, "phi"),
                        Internal(1, 4, "phi"),
                        Internal(2, 4, "phi"),
                        Internal(3, 4, "phi"),
                        Internal(3, 5, "phi"),
                    ],
                ),
                vec![0, 1, 2],
                vec![0, 1, 2, 4],
                vec![(6, source), (7, source)],
                3,
                0,
            ),
            (
                bridged,
                vec![1, 2, 4],
                vec![1, 2, 4],
                vec![(3, source), (4, source), (6, target), (7, source)],
                2,
                2,
            ),
        ] {
            let mut slots = vec![0; builder.vertices.len()];
            let mut candidate = DiagramThresholdCandidate {
                cut: Vec::new(),
                left: Vec::new(),
                right: Vec::new(),
            };
            for (edge, (ends, data)) in builder.edges.iter_mut().enumerate() {
                for (vertex, slot) in [
                    (ends.source, &mut data.source_slot),
                    (ends.target, &mut data.target_slot),
                ] {
                    if let Some(vertex) = vertex {
                        *slot = VertexSlot(slots[vertex.0]);
                        slots[vertex.0] += 1;
                    }
                }
                let mut left = Vec::new();
                for (endpoint, vertex) in [
                    (source, ends.source.or(ends.target).unwrap()),
                    (target, ends.target.or(ends.source).unwrap()),
                ] {
                    let half_edge = DiagramHalfEdge {
                        edge: EdgeId(edge),
                        endpoint,
                    };
                    if left_vertices.contains(&vertex.0) {
                        candidate.left.push(half_edge);
                        left.push(half_edge);
                    } else {
                        candidate.right.push(half_edge);
                    }
                }
                if left.len() == 1 {
                    candidate.cut.extend(left);
                }
            }
            assert_eq!(candidate.cut.len(), before);
            let builder = builder
                .canonical_loop_momentum_basis()
                .topology_threshold_candidates(vec![candidate]);
            let (graph, frame) = final_graph(&builder);
            let expected = expected
                .into_iter()
                .map(|(edge, endpoint)| {
                    let hedge = frame
                        .half_edge(DiagramHalfEdge {
                            edge: EdgeId(edge),
                            endpoint,
                        })
                        .unwrap();
                    DiagramHalfEdge {
                        edge: EdgeId(graph[&hedge].0),
                        endpoint: if graph.flow(hedge) == Flow::Source {
                            source
                        } else {
                            target
                        },
                    }
                })
                .collect::<BTreeSet<_>>();
            let diagram = builder.build().unwrap();
            // Full cross-section validation also requires its physical cut.
            // These are the captured cut and the existing bridged-bubble cut.
            let crowns = diagram.vertex_half_edges();
            let (left, right): (Vec<_>, Vec<_>) =
                (0..crowns.len()).partition(|vertex| physical_left.contains(vertex));
            let [left, right] = [left, right].map(|vertices| {
                vertices
                    .into_iter()
                    .flat_map(|vertex| crowns[vertex].clone())
                    .collect()
            });
            let diagram = diagram.with_cut_partitions(vec![(left, right)]).unwrap();
            diagram.validate().unwrap();
            let [candidate] = diagram.topology_threshold_candidates() else {
                panic!("one projected threshold");
            };
            assert_eq!(
                candidate.cut.iter().copied().collect::<BTreeSet<_>>(),
                expected
            );
            let ignored = diagram.initial_state_tree().0;
            assert_eq!(
                candidate
                    .cut
                    .iter()
                    .filter(|&&half_edge| {
                        let hedge = diagram.half_edge_id(half_edge).unwrap();
                        ignored.includes(&hedge) || ignored.includes(&diagram.graph.inv(hedge))
                    })
                    .count(),
                ignored_crossings,
            );
        }
    }

    #[test]
    fn generated_diagrams_reject_other_dangling_edges_legs_and_restored_frames() {
        let rejected = |builder: FeynmanDiagramBuilder| {
            matches!(
                builder.build(),
                Err(DiagramError::Invariant {
                    operation: "finalizing generated legs",
                    ..
                })
            )
        };
        let mut builder = case("relabelled");
        let particle = builder.edges[1].1.particle;
        builder
            .add_edge(VertexId(1), None, DiagramEdge::new(particle, false))
            .unwrap();
        assert!(rejected(builder));
        let mut builder = case("relabelled");
        builder.edges[1].1.external = builder.generation_externals.get(&EdgeId(0)).cloned();
        assert!(rejected(builder));
        let mut builder = case("relabelled");
        builder.half_edge_order = Some(Vec::new());
        assert!(rejected(builder));
    }
}
