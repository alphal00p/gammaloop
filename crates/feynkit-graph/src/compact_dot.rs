//! Resolve compact physics DOT against a model before ordinary diagram finalization.

use super::*;

impl FeynmanDiagram {
    pub(super) fn from_compact_dot(
        model: Arc<Model>,
        parsed: DotGraph,
    ) -> Result<Self, DiagramError> {
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
        let mut nodes = BTreeMap::new();
        for (node, _, data) in parsed.iter_nodes() {
            let numerator = Self::parse_expression(
                "vertex numerator",
                data.statements
                    .get("num")
                    .cloned()
                    .unwrap_or_else(|| "1".into()),
            )?;
            builder.numerator *= &numerator;
            nodes.insert(
                node,
                builder.add_vertex(DiagramVertex {
                    name: data.name.clone().unwrap_or_else(|| format!("v{}", node.0)),
                    interaction: None,
                    numerator,
                }),
            );
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
        let mut internal_order = Vec::new();
        let mut tagged_legs = BTreeMap::<usize, Vec<(Flow, ParticleId)>>::new();
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
                .ok_or_else(|| DiagramError::MissingDotAttribute {
                    target: target.clone(),
                    attribute: "particle or pdg",
                })?;
            let expected_orientation = if model.particle_is_self_conjugate(particle) {
                Orientation::Undirected
            } else if model.particle_by_id(particle)?.is_antiparticle() {
                Orientation::Reversed
            } else {
                Orientation::Default
            };
            if attributes.contains_key("dir") && data.orientation != expected_orientation {
                return Err(DiagramError::InvalidDotAttribute {
                    target,
                    attribute: "dir",
                    value: "explicit orientation disagrees with the model particle".into(),
                });
            }
            let mut edge = DiagramEdge::new(particle, !model.particle_is_self_conjugate(particle));
            edge.numerator = Self::parse_expression(
                "edge numerator",
                attributes.get("num").cloned().unwrap_or_else(|| "1".into()),
            )?;
            builder.numerator *= &edge.numerator;
            let tag = attributes
                .get("is_cut")
                .map(|value| value.parse::<usize>().expect("validated cut tag"));
            let (source, sink) = match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    if tag.is_some() {
                        return Err(DiagramError::InvalidDotAttribute {
                            target, attribute: "is_cut", value: "compact cut tags belong on matching incoming/outgoing dangling legs".into(),
                        });
                    }
                    internal_order.push(edge_id.0);
                    (
                        Some(nodes[&parsed.node_id(source)]),
                        Some(nodes[&parsed.node_id(sink)]),
                    )
                }
                HedgePair::Unpaired { hedge, flow } => {
                    if let Some(tag) = tag {
                        tagged_legs.entry(tag).or_default().push((flow, particle));
                    }
                    let index = next_external;
                    next_external += 1;
                    let connection = tag.map(|tag| cut_tags[&tag]).unwrap_or_else(|| {
                        let connection = next_connection;
                        next_connection += 1;
                        connection
                    });
                    let state = match flow {
                        Flow::Sink => ExternalState::Incoming,
                        Flow::Source => ExternalState::Outgoing,
                    };
                    let name = attributes
                        .get("name")
                        .cloned()
                        .unwrap_or_else(|| format!("ext{index}"));
                    let external = if sewn {
                        Some(builder.add_generation_external(name, index, state, connection))
                    } else {
                        edge.external = Some(ExternalLeg {
                            name,
                            index,
                            state,
                            connection,
                        });
                        None
                    };
                    let internal = Some(nodes[&parsed.node_id(hedge)]);
                    match flow {
                        Flow::Sink => (external, internal),
                        Flow::Source => (internal, external),
                    }
                }
            };
            let added = builder.add_edge(source, sink, edge)?;
            if let Some(value) = attributes.get("lmb_id") {
                let slot =
                    value
                        .parse::<usize>()
                        .map_err(|_| DiagramError::InvalidDotAttribute {
                            target: target.clone(),
                            attribute: "lmb_id",
                            value: value.clone(),
                        })?;
                if !internal_order.contains(&edge_id.0) || loop_edges.insert(slot, added).is_some()
                {
                    return Err(DiagramError::InvalidDotAttribute {
                        target,
                        attribute: "lmb_id",
                        value: "loop slots must be unique and belong to internal edges".into(),
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
        // Compact DOT omits vertex slots. Match the incident particles to a UFO
        // rule and assign slots in that rule's order before sewing external legs.
        for (node, _, data) in parsed.iter_nodes() {
            let vertex = nodes[&node];
            let mut incident = Vec::new();
            for (edge, (endpoints, data)) in builder.edges.iter().enumerate() {
                if endpoints.source == Some(vertex) {
                    incident.push((
                        edge,
                        true,
                        model.particle_by_id(data.particle)?.antiparticle,
                    ));
                }
                if endpoints.target == Some(vertex) {
                    incident.push((edge, false, data.particle));
                }
            }
            let mut signature = incident
                .iter()
                .map(|(_, _, particle)| *particle)
                .collect::<Vec<_>>();
            signature.sort();
            let explicit = data
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
                        && explicit.is_none_or(|explicit| explicit.index() == id))
                    .then_some(id)
                })
                .collect::<Vec<_>>();
            let [rule] = matches.as_slice() else {
                return Err(DiagramError::InvalidDotAttribute {
                    target: format!("vertex {}", builder.vertices[vertex.0].name),
                    attribute: "int_id",
                    value: format!(
                        "incident particles match {} model interactions; specify a matching int_id",
                        matches.len()
                    ),
                });
            };
            builder.vertices[vertex.0].interaction = Some(model.vertex_rule_id_at(*rule)?);
            for (slot, particle) in model.vertex_rules()[*rule].particles.iter().enumerate() {
                let index = incident
                    .iter()
                    .position(|(_, _, candidate)| candidate == particle)
                    .expect("matched particle signature");
                let (edge, source, _) = incident.remove(index);
                if source {
                    builder.edges[edge].1.source_slot = VertexSlot(slot);
                } else {
                    builder.edges[edge].1.target_slot = VertexSlot(slot);
                }
            }
        }
        if !sewn {
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
        let mut diagram = builder.build()?;
        if !loop_edges.is_empty() {
            if loop_edges.keys().copied().ne(0..loop_edges.len()) {
                return Err(DiagramError::InvalidLoopMomentumBasis(
                    "lmb_id slots must be contiguous from zero".into(),
                ));
            }
            let requested = loop_edges
                .into_values()
                .map(|edge| {
                    if sewn {
                        EdgeId(
                            diagram.graph.n_edges() - internal_order.len()
                                + internal_order
                                    .iter()
                                    .position(|id| *id == edge.0)
                                    .expect("internal loop edge"),
                        )
                    } else {
                        edge
                    }
                })
                .collect::<Vec<_>>();
            diagram = diagram.with_loop_momentum_edges(&requested)?;
        }
        if sewn {
            let final_state =
                attributes
                    .get("final_state")
                    .ok_or_else(|| DiagramError::MissingDotAttribute {
                        target: "compact cross-section graph".into(),
                        attribute: "final_state",
                    })?;
            let mut particles = final_state
                .split(',')
                .map(str::trim)
                .map(|name| model.particle_id(name))
                .collect::<Result<Vec<_>, _>>()?;
            particles.sort();
            let mut left_nodes = BTreeSet::new();
            let mut right_nodes = BTreeSet::new();
            for (_, endpoints, edge) in diagram.edges() {
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
                                value:
                                    "all cross-section external legs must have a matching cut tag"
                                        .into(),
                            });
                        }
                    }
                }
            }
            if !left_nodes.is_disjoint(&right_nodes) {
                return Err(DiagramError::InvalidCut {
                    cut: 0,
                    message: "incoming and outgoing attachments must lie on separate vertices"
                        .into(),
                });
            }
            let vertex_half_edges = diagram.vertex_half_edges();
            let mut cuts = Vec::new();
            let mut thresholds = Vec::new();
            for (left, _, right) in diagram.graph.all_cuts_from_ids(
                &left_nodes.into_iter().collect::<Vec<_>>(),
                &right_nodes.into_iter().collect::<Vec<_>>(),
            ) {
                let left = diagram
                    .half_edges()
                    .filter(|half| {
                        left.includes(&diagram.half_edge_id(*half).expect("existing half-edge"))
                    })
                    .collect::<Vec<_>>();
                let right = diagram
                    .half_edges()
                    .filter(|half| {
                        right.includes(&diagram.half_edge_id(*half).expect("existing half-edge"))
                    })
                    .collect::<Vec<_>>();
                let cut = diagram.cut_from_partitions(
                    cuts.len(),
                    &left,
                    &right,
                    &vertex_half_edges,
                    None,
                )?;
                let mut content = diagram.cut_particles(&cut)?;
                content.sort();
                if content == particles {
                    cuts.push(cut.clone());
                }
                if cut.cut.len() >= 2 {
                    thresholds.push((left, right));
                }
            }
            if cuts.is_empty() {
                return Err(DiagramError::InvalidDotAttribute {
                    target: "graph".into(),
                    attribute: "final_state",
                    value: format!("no separating cut has final state {final_state}"),
                });
            }
            diagram = diagram
                .with_cuts(cuts)?
                .with_topology_threshold_partitions(thresholds)?;
        } else if attributes.contains_key("final_state") {
            return Err(DiagramError::InvalidDotAttribute {
                target: "graph".into(),
                attribute: "final_state",
                value: "final_state requires paired is_cut external legs".into(),
            });
        }
        diagram.validate()?;
        Ok(diagram)
    }
}
