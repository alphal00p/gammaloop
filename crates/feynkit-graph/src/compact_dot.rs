//! Resolve compact physics DOT against a model before ordinary diagram finalization.

use super::finalization::{IndexOwner, IndexOwners};
use super::*;
use std::sync::LazyLock;
use symbolica::atom::{AtomView, FunctionBuilder, Symbol};
use symbolica::symbol;

/// Edge-numerator shorthand for the edge index and its endpoint half-edge indices.
static LOCAL_EDGE: LazyLock<Symbol> = LazyLock::new(|| symbol!("gammalooprs::eid"));
static LOCAL_SOURCE: LazyLock<Symbol> = LazyLock::new(|| symbol!("gammalooprs::source"));
static LOCAL_SINK: LazyLock<Symbol> = LazyLock::new(|| symbol!("gammalooprs::sink"));

impl FeynmanDiagram {
    /// Read compact DOT: particles by name or PDG code, optional numerators, cut tags
    /// and loop slots, with interactions and slots matched against the model.
    ///
    /// Matching `is_cut` tags on an incoming and an outgoing dangling leg, or on one
    /// internal edge, mark a cross-section initial state. Such a diagram keeps every
    /// separating cut as a physical cut, or only those whose particles are the
    /// graph's `final_state`. Expressions of a cross section refer to the DOT edge
    /// and half-edge indices; sewing renumbers both and transports the expressions.
    pub(super) fn from_compact_dot(
        model: Arc<Model>,
        parsed: DotGraph,
    ) -> Result<DotDiagram, DiagramError> {
        let attributes = &parsed.global_data.statements;
        let mut builder = Self::builder(Arc::clone(&model), parsed.global_data.name.clone());
        for (attribute, destination) in [
            ("num", &mut builder.numerator_prefactor),
            ("overall_factor", &mut builder.overall_factor),
            ("projector", &mut builder.projector),
        ] {
            if let Some(value) = attributes.get(attribute) {
                *destination = Self::parse_expression(attribute, value.clone())?;
            }
        }
        let mut explicit = ExplicitFactors {
            projector: attributes.contains_key("projector"),
            ..ExplicitFactors::default()
        };
        let mut vertex_attributes = Vec::new();
        let mut nodes = BTreeMap::new();
        for (node, _, data) in parsed.iter_nodes() {
            let numerator = data
                .statements
                .get("num")
                .map(|num| Self::parse_expression("vertex numerator", num.clone()))
                .transpose()?;
            let vertex = builder.add_vertex(DiagramVertex {
                name: data.name.clone().unwrap_or_else(|| format!("v{}", node.0)),
                interaction: None,
                numerator: numerator.clone().unwrap_or_else(Atom::one),
            });
            if numerator.is_some() {
                explicit.vertices.insert(vertex);
            }
            nodes.insert(node, vertex);
            vertex_attributes.push(data.statements.clone());
        }
        let mut cut_tags = BTreeMap::new();
        for (_, edge_id, data) in parsed.iter_edges() {
            if let Some(value) = data.data.statements.get("is_cut") {
                let tag =
                    value
                        .parse::<usize>()
                        .map_err(|_| DiagramError::InvalidDotAttribute {
                            target: format!("edge {}", edge_id.0),
                            attribute: "is_cut",
                            value: value.clone(),
                        })?;
                cut_tags.insert(tag, 0);
            }
        }
        for (connection, value) in cut_tags.values_mut().enumerate() {
            *value = connection;
        }
        let sewn = !cut_tags.is_empty();
        let mut next_connection = cut_tags.len();
        let mut next_external = 0;
        let mut loop_edges = BTreeMap::new();
        let mut internal = Vec::new();
        let mut explicit_edges = BTreeSet::new();
        let mut tagged_legs = BTreeMap::<usize, Vec<(Flow, ParticleId)>>::new();
        // Builder edges and half-edges of every DOT edge and half-edge, to transport
        // cross-section expressions and to index the edge attributes.
        let mut builder_edges = BTreeMap::new();
        let mut builder_half_edges = BTreeMap::new();
        let half_edge = |edge, endpoint| DiagramHalfEdge { edge, endpoint };
        let mut edge_attributes = Vec::new();
        for (pair, edge_id, data) in parsed.iter_edges() {
            let attributes = &data.data.statements;
            let target = format!("edge {}", edge_id.0);
            let by_name = attributes
                .get("particle")
                .map(|name| model.particle_id(name))
                .transpose()?;
            let by_pdg = attributes
                .get("pdg")
                .map(|value| {
                    let pdg =
                        value
                            .parse::<i64>()
                            .map_err(|_| DiagramError::InvalidDotAttribute {
                                target: target.clone(),
                                attribute: "pdg",
                                value: value.clone(),
                            })?;
                    model.particle_id_by_pdg(pdg).map_err(DiagramError::from)
                })
                .transpose()?;
            if by_name.is_some() && by_pdg.is_some() && by_name != by_pdg {
                return Err(DiagramError::InvalidDotAttribute {
                    target,
                    attribute: "pdg",
                    value: "particle and pdg refer to different model particles".into(),
                });
            }
            let particle = by_name
                .or(by_pdg)
                .ok_or_else(|| {
                    if attributes.contains_key("mass") {
                        DiagramError::InvalidDotAttribute {
                            target: target.clone(),
                            attribute: "mass",
                            value: "mass-only edges are unsupported; specify a model particle with `particle` or `pdg`".into(),
                        }
                    } else {
                        DiagramError::MissingDotAttribute {
                            target: target.clone(),
                            attribute: "particle or pdg",
                        }
                    }
                })?;
            let expected_orientation = if model.particle_is_self_conjugate(particle) {
                Orientation::Undirected
            } else if model.particle_by_id(particle)?.is_antiparticle() {
                Orientation::Reversed
            } else {
                Orientation::Default
            };
            // Only an edge's own `dir` is checked; a default `edge [dir=...]` is drawing style.
            if data.data.local_statements.contains_key("dir")
                && data.orientation != expected_orientation
            {
                return Err(DiagramError::InvalidDotAttribute {
                    target,
                    attribute: "dir",
                    value: "explicit orientation disagrees with the model particle".into(),
                });
            }
            let mut edge = DiagramEdge::new(particle, !model.particle_is_self_conjugate(particle));
            if let Some(value) = attributes.get("is_dummy") {
                edge.is_dummy =
                    value
                        .parse::<bool>()
                        .map_err(|_| DiagramError::InvalidDotAttribute {
                            target: target.clone(),
                            attribute: "is_dummy",
                            value: value.clone(),
                        })?;
            }
            let tag = attributes
                .get("is_cut")
                .map(|value| value.parse::<usize>().expect("validated cut tag"));
            let leg = |state, next_external: &mut usize, next_connection: &mut usize| {
                let index = *next_external;
                *next_external += 1;
                let connection = tag.map(|tag| cut_tags[&tag]).unwrap_or_else(|| {
                    *next_connection += 1;
                    *next_connection - 1
                });
                ExternalLeg {
                    name: attributes
                        .get("name")
                        .cloned()
                        .unwrap_or_else(|| format!("ext{index}")),
                    index,
                    state,
                    connection,
                }
            };
            let added = match (pair, tag) {
                (
                    HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. },
                    Some(tag),
                ) => {
                    // A cut internal edge is a sewn initial state: its target receives the
                    // incoming copy and its source emits the outgoing one.
                    let incoming = builder.add_generation_external(
                        nodes[&parsed.node_id(sink)],
                        Flow::Sink,
                        edge.clone(),
                        VertexSlot(0),
                        leg(
                            ExternalState::Incoming,
                            &mut next_external,
                            &mut next_connection,
                        ),
                    )?;
                    let outgoing = builder.add_generation_external(
                        nodes[&parsed.node_id(source)],
                        Flow::Source,
                        edge,
                        VertexSlot(0),
                        leg(
                            ExternalState::Outgoing,
                            &mut next_external,
                            &mut next_connection,
                        ),
                    )?;
                    tagged_legs
                        .entry(tag)
                        .or_default()
                        .extend([(Flow::Sink, particle), (Flow::Source, particle)]);
                    builder_half_edges.insert(sink.0, half_edge(incoming, DiagramEndpoint::Target));
                    builder_half_edges
                        .insert(source.0, half_edge(outgoing, DiagramEndpoint::Source));
                    vec![incoming, outgoing]
                }
                (
                    HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. },
                    None,
                ) => {
                    // Only propagators carry a numerator fragment: the `num` of an
                    // external leg or a sewn initial state is ignored like any other
                    // attribute FeynKit does not use.
                    let numerator = attributes
                        .get("num")
                        .map(|num| Self::parse_expression("edge numerator", num.clone()))
                        .transpose()?;
                    let mut edge = edge;
                    let numerator_given = numerator.is_some();
                    if let Some(numerator) = numerator {
                        edge.numerator = numerator;
                    }
                    let added = builder.add_edge(
                        nodes[&parsed.node_id(source)],
                        nodes[&parsed.node_id(sink)],
                        edge,
                    )?;
                    if numerator_given {
                        explicit_edges.insert(added);
                    }
                    internal.push(added);
                    builder_half_edges.insert(source.0, half_edge(added, DiagramEndpoint::Source));
                    builder_half_edges.insert(sink.0, half_edge(added, DiagramEndpoint::Target));
                    vec![added]
                }
                (HedgePair::Unpaired { hedge, flow }, _) => {
                    if let Some(tag) = tag {
                        tagged_legs.entry(tag).or_default().push((flow, particle));
                    }
                    let state = match flow {
                        Flow::Sink => ExternalState::Incoming,
                        Flow::Source => ExternalState::Outgoing,
                    };
                    let internal_vertex = nodes[&parsed.node_id(hedge)];
                    let added = if sewn {
                        builder.add_generation_external(
                            internal_vertex,
                            flow,
                            edge,
                            VertexSlot(0),
                            leg(state, &mut next_external, &mut next_connection),
                        )?
                    } else {
                        if !edge.is_dummy {
                            edge.external =
                                Some(leg(state, &mut next_external, &mut next_connection));
                        }
                        let (source, target) = match flow {
                            Flow::Source => (Some(internal_vertex), None),
                            Flow::Sink => (None, Some(internal_vertex)),
                        };
                        builder.add_edge(source, target, edge)?
                    };
                    let endpoint = match flow {
                        Flow::Source => DiagramEndpoint::Source,
                        Flow::Sink => DiagramEndpoint::Target,
                    };
                    builder_half_edges.insert(hedge.0, half_edge(added, endpoint));
                    vec![added]
                }
            };
            builder_edges.insert(edge_id.0, added[0]);
            for added in &added {
                edge_attributes.push((*added, attributes.clone()));
            }
            if let Some(value) = attributes.get("lmb_id") {
                let slot =
                    value
                        .parse::<usize>()
                        .map_err(|_| DiagramError::InvalidDotAttribute {
                            target: target.clone(),
                            attribute: "lmb_id",
                            value: value.clone(),
                        })?;
                if !internal.contains(&added[0]) || loop_edges.insert(slot, added[0]).is_some() {
                    return Err(DiagramError::InvalidDotAttribute {
                        target,
                        attribute: "lmb_id",
                        value: "loop slots must be unique and belong to uncut internal edges"
                            .into(),
                    });
                }
            }
        }
        for (tag, legs) in &tagged_legs {
            if legs.len() != 2 || legs[0].0 == legs[1].0 || legs[0].1 != legs[1].1 {
                return Err(DiagramError::InvalidDotAttribute {
                    target: "external legs".into(),
                    attribute: "is_cut",
                    value: format!(
                        "tag {tag} must pair one incoming and one outgoing leg of the same particle"
                    ),
                });
            }
        }
        Self::assign_compact_interactions(&model, &parsed, &nodes, &mut builder, &explicit)?;
        if sewn {
            // Finalization relabels the expressions of a sewn diagram from their
            // provisional indices, so move the DOT indices to those first.
            let owners = IndexOwners::new();
            let transport = |atom: &Atom| -> Result<Atom, DiagramError> {
                owners
                    .relabel(atom, |owner, index| {
                        match owner {
                            IndexOwner::HalfEdge => builder_half_edges
                                .get(&index)
                                .map(|half_edge| half_edge.provisional_index()),
                            IndexOwner::Vertex => {
                                nodes.get(&NodeIndex(index)).map(|vertex| vertex.0)
                            }
                            IndexOwner::Edge | IndexOwner::Momentum | IndexOwner::ExternalState => {
                                builder_edges.get(&index).map(|edge| edge.0)
                            }
                        }
                        .map(|index| (index, false))
                    })
                    .map_err(|message| DiagramError::Invariant {
                        operation: "transporting DOT expression indices",
                        message,
                    })
            };
            for vertex in &mut builder.vertices {
                vertex.numerator = transport(&vertex.numerator)?;
            }
            for (_, edge) in &mut builder.edges {
                edge.numerator = transport(&edge.numerator)?;
            }
            builder.numerator_prefactor = transport(&builder.numerator_prefactor)?;
            builder.projector = transport(&builder.projector)?;
            builder.overall_factor = transport(&builder.overall_factor)?;
        } else {
            builder.half_edge_order = Some(
                parsed
                    .iter_hedges()
                    .map(|(hedge, _)| DiagramHalfEdge {
                        edge: EdgeId(parsed[&hedge].0),
                        endpoint: match parsed.flow(hedge) {
                            Flow::Source => DiagramEndpoint::Source,
                            Flow::Sink => DiagramEndpoint::Target,
                        },
                    })
                    .collect(),
            );
        }
        // The two legs of a sewn initial state merge into one edge, whose attributes
        // are those of the incoming leg, completed by the outgoing one's.
        edge_attributes.sort_by_key(|(edge, _)| {
            builder
                .generation_externals
                .get(edge)
                .is_some_and(|leg| leg.state == ExternalState::Outgoing)
        });
        let (mut diagram, finalized) = builder.finalize()?;
        let mut edges = vec![BTreeMap::new(); diagram.graph.n_edges()];
        for (edge, attributes) in edge_attributes {
            let merged: &mut BTreeMap<String, String> = &mut edges[finalized[edge.0].0];
            for (key, value) in attributes {
                merged.entry(key).or_insert(value);
            }
        }
        explicit.edges = explicit_edges
            .into_iter()
            .map(|edge| finalized[edge.0])
            .collect();
        let requested = loop_edges
            .values()
            .map(|edge| finalized[edge.0])
            .collect::<Vec<_>>();
        if !loop_edges.is_empty() {
            if loop_edges.keys().copied().ne(0..loop_edges.len()) {
                return Err(DiagramError::InvalidLoopMomentumBasis(
                    "lmb_id slots must be contiguous from zero".into(),
                ));
            }
            diagram = diagram.with_loop_momentum_edges(&requested)?;
        }
        if sewn {
            diagram = diagram.with_compact_cuts(attributes.get("final_state"))?;
        } else if attributes.contains_key("final_state") {
            return Err(DiagramError::InvalidDotAttribute {
                target: "graph".into(),
                attribute: "final_state",
                value: "final_state requires paired is_cut external legs".into(),
            });
        }
        // Edge shorthand refers to the finalized frame, like transported expressions.
        let local = diagram
            .edges()
            .filter_map(|(edge, _, data)| {
                let resolved = diagram.resolve_local_edge_indices(&data.numerator, edge);
                (resolved != data.numerator).then_some((edge, resolved))
            })
            .collect::<Vec<_>>();
        if !local.is_empty() {
            diagram = diagram.with_numerator_fragments([], local)?;
        }
        if loop_edges.is_empty() {
            // Without explicit slots, route with generation's canonical basis for the
            // cuts kept here, which keeps beam momenta off cut edges where possible.
            diagram = diagram.with_canonical_loop_momentum_basis()?;
        }
        diagram.validate()?;
        Ok(DotDiagram {
            diagram,
            attributes: DotAttributes {
                graph: attributes.clone(),
                vertices: vertex_attributes,
                edges,
            },
            explicit,
        })
    }

    /// Match every vertex to a model interaction and assign its endpoint slots in the
    /// interaction's leg order. Dummy edges take the slots after the legs. A vertex
    /// with an explicit numerator may match no interaction; its slots then follow
    /// the incident half-edges.
    fn assign_compact_interactions(
        model: &Model,
        parsed: &DotGraph,
        nodes: &BTreeMap<NodeIndex, VertexId>,
        builder: &mut FeynmanDiagramBuilder,
        explicit: &ExplicitFactors,
    ) -> Result<(), DiagramError> {
        for (node, _, data) in parsed.iter_nodes() {
            let vertex = nodes[&node];
            let mut incident = Vec::new();
            let mut dummies = Vec::new();
            for (edge, (endpoints, data)) in builder.edges.iter().enumerate() {
                for (attached, source, particle) in [
                    (
                        endpoints.source,
                        true,
                        model.particle_by_id(data.particle)?.antiparticle,
                    ),
                    (endpoints.target, false, data.particle),
                ] {
                    if attached == Some(vertex) {
                        if data.is_dummy {
                            dummies.push((edge, source));
                        } else {
                            incident.push((edge, source, particle));
                        }
                    }
                }
            }
            let mut signature = incident
                .iter()
                .map(|(_, _, particle)| *particle)
                .collect::<Vec<_>>();
            signature.sort();
            let requested = data
                .statements
                .get("int_id")
                .map(|name| model.vertex_rule_id(name))
                .transpose()?;
            let matches = model
                .vertex_rules()
                .iter()
                .enumerate()
                .filter_map(|(id, rule)| {
                    let mut particles = rule.particles.clone();
                    particles.sort();
                    (particles == signature
                        && requested.is_none_or(|requested| requested.index() == id))
                    .then_some(id)
                })
                .collect::<Vec<_>>();
            let mut slots = Vec::new();
            match matches.as_slice() {
                [rule] => {
                    builder.vertices[vertex.0].interaction = Some(model.vertex_rule_id_at(*rule)?);
                    for particle in &model.vertex_rules()[*rule].particles {
                        let index = incident
                            .iter()
                            .position(|(_, _, candidate)| candidate == particle)
                            .expect("matched particle signature");
                        let (edge, source, _) = incident.remove(index);
                        slots.push((edge, source));
                    }
                }
                [] if requested.is_none() && explicit.vertices.contains(&vertex) => {
                    slots.extend(incident.into_iter().map(|(edge, source, _)| (edge, source)));
                }
                _ => {
                    return Err(DiagramError::InvalidDotAttribute {
                        target: format!("vertex {}", builder.vertices[vertex.0].name),
                        attribute: "int_id",
                        value: format!(
                            "incident particles match {} model interactions; specify a matching int_id{}",
                            matches.len(),
                            if matches.is_empty() {
                                " or the vertex numerator"
                            } else {
                                ""
                            }
                        ),
                    });
                }
            }
            for (slot, (edge, source)) in slots.into_iter().chain(dummies).enumerate() {
                let data = &mut builder.edges[edge].1;
                if source {
                    data.source_slot = VertexSlot(slot);
                } else {
                    data.target_slot = VertexSlot(slot);
                }
            }
        }
        Ok(())
    }

    /// Resolve an edge numerator's shorthand: `eid`, `source` and `sink` are the index
    /// of `edge` and of its endpoint half-edges in this diagram, and `eid(i)`,
    /// `source(i)` and `sink(i)` the corresponding tensor-index identities.
    fn resolve_local_edge_indices(&self, numerator: &Atom, edge: EdgeId) -> Atom {
        let endpoint = |endpoint| {
            self.half_edge_id(DiagramHalfEdge { edge, endpoint })
                .map(|hedge| hedge.0)
        };
        let shorthand = [
            (*LOCAL_EDGE, Some(edge.0), symbols::edge_index()),
            (
                *LOCAL_SOURCE,
                endpoint(DiagramEndpoint::Source),
                symbols::hedge_index(),
            ),
            (
                *LOCAL_SINK,
                endpoint(DiagramEndpoint::Target),
                symbols::hedge_index(),
            ),
        ];
        numerator.replace_map(|term, _, out| {
            let (symbol, local) = match term {
                AtomView::Var(var) => (var.get_symbol(), None),
                AtomView::Fun(function) if function.get_nargs() == 1 => (
                    function.get_symbol(),
                    Some(function.iter().next().expect("one argument")),
                ),
                _ => return,
            };
            let Some((_, Some(index), head)) = shorthand.iter().find(|(name, ..)| *name == symbol)
            else {
                return;
            };
            **out = match local {
                None => Atom::num(*index as i64),
                Some(local) if local.is_zero() => {
                    FunctionBuilder::new(*head).add_arg(*index).finish()
                }
                Some(local) => FunctionBuilder::new(*head)
                    .add_arg(*index)
                    .add_arg(local)
                    .finish(),
            };
        })
    }

    /// Install the separating cuts of a sewn compact diagram: those whose particles
    /// are `final_state`, if given, or else all of them. Every cut through at least
    /// two edges is also a topology threshold candidate.
    fn with_compact_cuts(self, final_state: Option<&String>) -> Result<Self, DiagramError> {
        let mut particles = final_state
            .map(|final_state| {
                final_state
                    .split(',')
                    .map(str::trim)
                    .map(|name| self.model.particle_id(name))
                    .collect::<Result<Vec<_>, _>>()
            })
            .transpose()?;
        if let Some(particles) = &mut particles {
            particles.sort();
        }
        let mut left_nodes = BTreeSet::new();
        let mut right_nodes = BTreeSet::new();
        for (_, endpoints, edge) in self.edges() {
            if edge.external.is_some() {
                match (endpoints.source, endpoints.target) {
                    (Some(source), Some(target)) => {
                        left_nodes.insert(NodeIndex(target.0));
                        right_nodes.insert(NodeIndex(source.0));
                    }
                    _ => {
                        return Err(DiagramError::InvalidDotAttribute {
                            target: "external legs".into(),
                            attribute: "is_cut",
                            value: "all cross-section external legs must have a matching cut tag"
                                .into(),
                        });
                    }
                }
            }
        }
        if !left_nodes.is_disjoint(&right_nodes) {
            return Err(DiagramError::InvalidCut {
                cut: 0,
                message: "incoming and outgoing attachments must lie on separate vertices".into(),
            });
        }
        let vertex_half_edges = self.vertex_half_edges();
        let mut cuts = Vec::new();
        let mut thresholds = Vec::new();
        for (left, _, right) in self.graph.all_cuts_from_ids(
            &left_nodes.into_iter().collect::<Vec<_>>(),
            &right_nodes.into_iter().collect::<Vec<_>>(),
        ) {
            let side = |selected: &SuBitGraph| {
                self.half_edges()
                    .filter(|half| {
                        selected.includes(&self.half_edge_id(*half).expect("existing half-edge"))
                    })
                    .collect::<Vec<_>>()
            };
            let (left, right) = (side(&left), side(&right));
            let cut = self.cut_from_partitions(cuts.len(), &left, &right, &vertex_half_edges)?;
            let selected = match &particles {
                Some(particles) => {
                    let mut content = self.cut_particles(&cut)?;
                    content.sort();
                    &content == particles
                }
                None => true,
            };
            if cut.cut.len() >= 2 {
                thresholds.push((left, right));
            }
            if selected {
                cuts.push(cut);
            }
        }
        if cuts.is_empty() {
            return Err(DiagramError::InvalidDotAttribute {
                target: "graph".into(),
                attribute: "final_state",
                value: format!(
                    "no separating cut has final state {}",
                    final_state.map_or("of any particles", String::as_str)
                ),
            });
        }
        self.with_cuts(cuts)?
            .with_topology_threshold_partitions(thresholds)
    }
}
