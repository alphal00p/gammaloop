//! Resolve compact physics DOT against a model before ordinary diagram finalization.

use super::*;
use linnet::half_edge::involution::Hedge;
use linnet::half_edge::subgraph::ModifySubSet;
use std::ops::RangeInclusive;

impl FeynmanDiagram {
    /// Resolve compact physics DOT and retain the exact source-to-final map
    /// used by finalization. Consumers can transport graph-local metadata
    /// without reconstructing ownership from momentum or particle labels.
    /// `final_state` lists comma-separated particles; `|` joins allowed alternatives.
    pub fn from_compact_dot(
        model: Arc<Model>,
        parsed: DotGraph,
    ) -> Result<(Self, DiagramSourceMap), DiagramError> {
        let attributes = &parsed.global_data.statements;
        let mut builder = Self::builder(Arc::clone(&model), parsed.global_data.name.clone());
        if let Some(value) = attributes.get("overall_factor") {
            builder.overall_factor = Self::parse_expression("overall_factor", value.clone())?;
        }
        for (attribute, metadata, destination) in [
            (
                "num",
                "numerator_prefactor_interface",
                &mut builder.numerator_prefactor,
            ),
            ("projector", "projector_interface", &mut builder.projector),
        ] {
            if let Some(value) = attributes.get(attribute) {
                *destination = NumeratorInterface::admit(
                    Self::parse_expression(attribute, value.clone())?,
                    attributes
                        .get(metadata)
                        .map(|value| NumeratorInterface::parse(value))
                        .transpose()?,
                )?;
            }
        }
        let mut nodes = BTreeMap::new();
        for (node, _, data) in parsed.iter_nodes() {
            let numerator = NumeratorInterface::admit(
                Self::parse_expression(
                    "vertex numerator",
                    data.statements
                        .get("num")
                        .cloned()
                        .unwrap_or_else(|| "1".into()),
                )?,
                data.statements
                    .get("num_interface")
                    .map(|value| NumeratorInterface::parse(value))
                    .transpose()?,
            )?;
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
            edge.numerator = NumeratorInterface::admit(
                Self::parse_expression(
                    "edge numerator",
                    attributes.get("num").cloned().unwrap_or_else(|| "1".into()),
                )?,
                data.data
                    .statements
                    .get("num_interface")
                    .map(|value| NumeratorInterface::parse(value))
                    .transpose()?,
            )?;
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
        } else {
            // Supplied expressions use the parsed DOT half-edge chart. The
            // generated-vertex finalizer consumes the builder's edge/end chart.
            // Transport this interface from recorded incidence before sewing;
            // the finalizer then applies its existing signed ownership map.
            let source_hedges = parsed
                .iter_hedges()
                .map(|(hedge, _)| {
                    let edge = parsed[&hedge];
                    let endpoint = usize::from(parsed.flow(hedge) == Flow::Sink);
                    (hedge.0, 2 * edge.0 + endpoint)
                })
                .collect::<BTreeMap<_, _>>();
            let translate = |atom: &Atom| -> Result<Atom, DiagramError> {
                let mut invalid = None;
                let translated = atom.replace_map(|part, _, output| {
                    let symbolica::atom::AtomView::Fun(index) = part else {
                        return;
                    };
                    if index.get_symbol() != symbols::hedge_index() {
                        return;
                    }
                    let Some(source) = index
                        .iter()
                        .next()
                        .and_then(|arg| usize::try_from(arg).ok())
                    else {
                        return;
                    };
                    let Some(&target) = source_hedges.get(&source) else {
                        invalid = Some(source);
                        return;
                    };
                    **output = symbolica::atom::FunctionBuilder::new(symbols::hedge_index())
                        .add_arg(target)
                        .add_args(index.iter().skip(1))
                        .finish();
                });
                if let Some(source) = invalid {
                    return Err(DiagramError::InvalidDotAttribute {
                        target: "compact numerator interface".into(),
                        attribute: "num",
                        value: format!("half-edge {source} is outside the parsed source graph"),
                    });
                }
                Ok(translated)
            };
            builder.overall_factor = translate(&builder.overall_factor)?;
            builder.numerator_prefactor =
                NumeratorInterface::translate(&builder.numerator_prefactor, translate)?;
            builder.projector = NumeratorInterface::translate(&builder.projector, translate)?;
            for vertex in &mut builder.vertices {
                vertex.numerator = NumeratorInterface::translate(&vertex.numerator, translate)?;
            }
            for (_, edge) in &mut builder.edges {
                edge.numerator = NumeratorInterface::translate(&edge.numerator, translate)?;
            }
        }
        let (mut diagram, source_map) = builder.finalize()?;
        if attributes.contains_key("projector") && !attributes.contains_key("projector_interface") {
            // No logical order was declared. Use named algebra order in the
            // final chart after checked transport has retained index history.
            diagram.projector = diagram.projector.canonicalize_named_interface();
        }
        if !loop_edges.is_empty() {
            if loop_edges.keys().copied().ne(0..loop_edges.len()) {
                return Err(DiagramError::InvalidLoopMomentumBasis(
                    "lmb_id slots must be contiguous from zero".into(),
                ));
            }
            let requested = loop_edges
                .into_values()
                .map(|edge| source_map.edges[edge.0])
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
            let final_states = final_state
                .split('|')
                .map(|alternative| {
                    let mut particles = alternative
                        .split(',')
                        .map(str::trim)
                        .map(|name| model.particle_id(name))
                        .collect::<Result<Vec<_>, _>>()?;
                    particles.sort();
                    Ok::<_, ModelError>(particles)
                })
                .collect::<Result<BTreeSet<_>, _>>()?;
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
            let blob_range = match attributes.get("n_cut_blobs") {
                None => 1..=1,
                Some(value) => {
                    parse_blob_range(value).ok_or_else(|| DiagramError::InvalidDotAttribute {
                        target: "graph".into(),
                        attribute: "n_cut_blobs",
                        value: format!(
                            "{value:?} is not a blob count or an ordered range \"min,max\""
                        ),
                    })?
                }
            };
            // The unsewn structure: every sewn external edge r -> l becomes
            // the legs leg_in -> l and r -> leg_out, after the diagram
            // vertices. Internal edges keep their diagram edge.
            let n_nodes = diagram.graph.n_nodes();
            let mut unsewn = HedgeGraphBuilder::new();
            for _ in 0..n_nodes {
                unsewn.add_node(());
            }
            let mut incoming = Vec::new();
            let mut outgoing = Vec::new();
            for (edge_id, endpoints, edge) in diagram.edges() {
                let (Some(source), Some(target)) = (endpoints.source, endpoints.target) else {
                    continue;
                };
                let (source, target) = (NodeIndex(source.0), NodeIndex(target.0));
                if edge.external.is_some() {
                    let leg = unsewn.add_node(());
                    unsewn.add_edge(leg, target, None, Orientation::Default);
                    incoming.push(leg);
                    let leg = unsewn.add_node(());
                    unsewn.add_edge(source, leg, None, Orientation::Default);
                    outgoing.push(leg);
                } else {
                    unsewn.add_edge(source, target, Some(edge_id), Orientation::Default);
                }
            }
            let unsewn: HedgeGraph<Option<EdgeId>, ()> = unsewn.build();
            let forward = ForwardGraph::new(&unsewn, &incoming, &outgoing).map_err(|error| {
                DiagramError::InvalidCut {
                    cut: 0,
                    message: error.to_string(),
                }
            })?;
            let unsewn_half_edges: BTreeMap<EdgeId, Hedge> = unsewn
                .iter_edges()
                .filter_map(|(pair, _, data)| match (pair, *data.data) {
                    (HedgePair::Paired { source, .. }, Some(edge)) => Some((edge, source)),
                    _ => None,
                })
                .collect();
            // (partition, found only by the search with the legs fixed)
            let mut partitions = Vec::new();
            let mut enumerated = BTreeSet::new();
            for (_, cut, _) in diagram.graph.all_cuts_from_ids(
                &left_nodes.iter().copied().collect::<Vec<_>>(),
                &right_nodes.iter().copied().collect::<Vec<_>>(),
            ) {
                // The sewn search keeps the far half of a bridge it skips on
                // the left. Rebuild whole vertex sides on the unsewn structure
                // from the internal edges the cut crosses.
                let mut boundary: SuBitGraph = unsewn.empty_subgraph();
                for hedge in cut.left.included_iter() {
                    if let Some(&half) = unsewn_half_edges.get(&EdgeId(diagram.graph[&hedge].0)) {
                        boundary.add(half);
                    }
                }
                if let Some(partition) = forward.partition(&boundary) {
                    enumerated.insert(partition.cut().clone());
                    partitions.push((partition, false));
                }
            }
            // As in generation, also enumerate with one fixed degree-one vertex
            // per leg. The search above grows the left side one vertex at a
            // time and drops every step whose right side has more components
            // than the fixed right vertices. With adjacent outgoing
            // attachments, a subgraph attached only to incoming attachments
            // (such as a self-energy on the internal line between them) is a
            // second right component from the start, so it never reaches a
            // partition that places it on the left.
            for found in forward.leg_seeded_cuts(|boundary| enumerated.contains(boundary)) {
                if let LegSeededCut::Partition(partition) = found {
                    partitions.push((partition, true));
                }
            }
            let vertex_half_edges = diagram.vertex_half_edges();
            let mut cuts = Vec::new();
            let mut thresholds = Vec::new();
            for (partition, recovered) in partitions {
                let left_vertices = (0..n_nodes)
                    .map(NodeIndex)
                    .filter(|&node| {
                        unsewn
                            .iter_crown(node)
                            .any(|hedge| partition.left().includes(&hedge))
                    })
                    .collect::<BTreeSet<_>>();
                let (left, right): (Vec<_>, Vec<_>) = diagram.half_edges().partition(|half| {
                    let hedge = diagram.half_edge_id(*half).expect("existing half-edge");
                    left_vertices.contains(&diagram.graph.node_id(hedge))
                });
                // Blobs of each side, counted with its legs as generation
                // counts them.
                let left_blobs = unsewn.count_connected_components(partition.left());
                let right_blobs = unsewn.count_connected_components(partition.right());
                let cut = diagram.cut_from_partitions(
                    cuts.len(),
                    &left,
                    &right,
                    &vertex_half_edges,
                    None,
                )?;
                let mut content = diagram.cut_particles(&cut)?;
                content.sort();
                // Physical cuts follow generation: both sides within the blob
                // range, with the requested final state.
                if blob_range.contains(&left_blobs)
                    && blob_range.contains(&right_blobs)
                    && final_states.contains(&content)
                {
                    cuts.push(cut.clone());
                }
                // A recovered threshold candidate needs connected sides,
                // independently of the blob range.
                if cut.cut.len() >= 2 && (!recovered || (left_blobs == 1 && right_blobs == 1)) {
                    thresholds.push((left, right));
                }
            }
            if cuts.is_empty() {
                return Err(DiagramError::InvalidDotAttribute {
                    target: "graph".into(),
                    attribute: "final_state",
                    value: format!(
                        "no separating cut has final state {final_state} with {}..={} blobs per side",
                        blob_range.start(),
                        blob_range.end()
                    ),
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
        Ok((diagram, source_map))
    }
}

/// A compact `n_cut_blobs` value: one blob count, or the ordered bounds
/// `min,max`, as `--n-cut-blobs` accepts them in generation.
fn parse_blob_range(value: &str) -> Option<RangeInclusive<usize>> {
    let mut bounds = value
        .split(',')
        .map(|bound| bound.trim().parse::<usize>().ok());
    let minimum = bounds.next()??;
    let maximum = bounds.next().unwrap_or(Some(minimum))?;
    (bounds.next().is_none() && minimum <= maximum).then_some(minimum..=maximum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::structure::abstract_index::AbstractIndex;
    use spenso::structure::partial::PartialIndex;
    use spenso::structure::representation::LibraryRep;
    use spenso::structure::slot::IsAbstractSlot;
    use spenso::tensor_symbol;
    use symbolica::function;

    #[test]
    fn compact_projector_layout_distinguishes_inference_from_declaration() {
        let model = Arc::new(
            Model::from_json(include_str!(
                "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
            ))
            .unwrap(),
        );
        let input = r#"digraph projector_bubble {
            final_state="scalar_1,scalar_1"
            ext [style=invis]
            ext -> a:5 [id=0 particle=scalar_1 is_cut=7]
            a:1 -> b:2 [id=1 particle=scalar_1 lmb_id=0]
            b:0 -> ext [id=2 particle=scalar_1 is_cut=7]
            a:3 -> b:4 [id=3 particle=scalar_1]
        }"#;
        let rep = Minkowski {}.new_rep(4).cast::<LibraryRep>();
        let index = |owner| AbstractIndex::Named(symbols::hedge_index().into(), owner, 1);
        let slot = |owner| rep.slot(PartialIndex::Explicit(index(owner)));
        let slot_atom = |owner| rep.slot::<AbstractIndex, _>(index(owner)).to_atom();
        let head = tensor_symbol!("compact_projector_order");
        let source_owners = [0, 5];
        let payload = function!(head, slot_atom(0), slot_atom(5));
        for declaration in [None, Some([0, 1]), Some([1, 0])] {
            let mut parsed = DotGraph::from_string(input).unwrap();
            parsed
                .global_data
                .statements
                .insert("projector".into(), payload.to_canonical_string());
            if let Some(order) = declaration {
                let interface = NumeratorInterface(
                    order
                        .into_iter()
                        .map(|position| slot_atom(source_owners[position]))
                        .collect(),
                );
                parsed
                    .global_data
                    .statements
                    .insert("projector_interface".into(), interface.to_json());
            }
            let (diagram, source_map) =
                FeynmanDiagram::from_compact_dot(Arc::clone(&model), parsed.clone()).unwrap();
            // Source h0 belongs to the outgoing leg, h5 to the incoming leg.
            // Resolve their finalized endpoints independently of tensor transport.
            let target_owners = [
                DiagramHalfEdge {
                    edge: source_map.edges()[2],
                    endpoint: DiagramEndpoint::Source,
                },
                DiagramHalfEdge {
                    edge: source_map.edges()[0],
                    endpoint: DiagramEndpoint::Target,
                },
            ]
            .map(|half| diagram.half_edge_id(half).unwrap().0);
            assert_eq!(target_owners, [5, 0]);
            let target_ports = target_owners.map(slot);
            let expected_payload = function!(
                head,
                slot_atom(target_owners[0]),
                slot_atom(target_owners[1])
            );
            let expected_order = declaration.unwrap_or([1, 0]);
            let expected_ports = expected_order.map(|position| target_ports[position]);
            assert_eq!(diagram.projector().expression(), &expected_payload);
            assert_eq!(
                diagram.projector().structure().logical_slots(),
                expected_ports
            );
            if declaration.is_none() {
                assert_eq!(diagram.projector().presentation_atom(), expected_payload);
                let (exported, interface) = NumeratorInterface::export(diagram.projector());
                assert!(
                    interface.is_some(),
                    "leaf syntax alone has a different logical order"
                );
                let restored = NumeratorInterface::admit(exported, interface).unwrap();
                assert!(restored.same_value(diagram.projector()));
            } else {
                // A declared ranked zero has no payload from which to recover order.
                parsed
                    .global_data
                    .statements
                    .insert("projector".into(), "0".into());
                let (zero, _) =
                    FeynmanDiagram::from_compact_dot(Arc::clone(&model), parsed).unwrap();
                assert!(zero.projector().is_zero());
                assert_eq!(zero.projector().structure().logical_slots(), expected_ports);
            }
        }

        // An undeclared positional interface must bypass named-axis sorting.
        let mut parsed = DotGraph::from_string(input).unwrap();
        let payload = function!(head, rep.to_symbolic([]), slot_atom(0));
        let source = NumeratorInterface::admit(payload.clone(), None).unwrap();
        parsed
            .global_data
            .statements
            .insert("projector".into(), payload.to_canonical_string());
        let (diagram, source_map) = FeynmanDiagram::from_compact_dot(model, parsed).unwrap();
        let target_owner = diagram
            .half_edge_id(DiagramHalfEdge {
                edge: source_map.edges()[2],
                endpoint: DiagramEndpoint::Source,
            })
            .unwrap()
            .0;
        let mut expected_ports = source.structure().logical_slots();
        expected_ports[1] = slot(target_owner);
        let expected_payload = source
            .expression()
            .replace(slot_atom(0).to_pattern())
            .with(slot_atom(target_owner));
        assert_eq!(
            diagram.projector().structure().logical_slots(),
            expected_ports
        );
        assert_eq!(diagram.projector().structure().open_positions(), [0]);
        assert_eq!(diagram.projector().expression(), &expected_payload);
    }

    #[test]
    fn compact_final_state_alternatives_preserve_the_union_of_physical_cuts() {
        let model = Arc::new(
            Model::from_json(include_str!(
                "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
            ))
            .unwrap(),
        );
        let input = r#"digraph alternatives {
            final_state="FINAL_STATE";
            ext_in [style=invis]; ext_out [style=invis];
            ext_in -> a [id=0 particle=scalar_1 is_cut=7];
            a -> b [id=1 particle=scalar_1];
            a -> d [id=2 particle=scalar_1];
            b -> c [id=3 particle=scalar_1];
            b -> c [id=4 particle=scalar_2];
            c -> d [id=5 particle=scalar_1];
            d -> ext_out [id=6 particle=scalar_1 is_cut=7];
        }"#;
        let load = |states: &str| {
            FeynmanDiagram::from_compact_dot(
                Arc::clone(&model),
                DotGraph::from_string(input.replace("FINAL_STATE", states)).unwrap(),
            )
            .map(|(diagram, _)| diagram)
        };
        let two = load("scalar_1,scalar_1").unwrap();
        let three = load("scalar_1,scalar_1,scalar_2").unwrap();
        let both =
            load("scalar_1,scalar_1 | scalar_1,scalar_2,scalar_1 | scalar_1,scalar_1").unwrap();
        assert_eq!(two.cuts().len(), 2);
        assert_eq!(three.cuts().len(), 1);
        assert_eq!(
            both.cuts().iter().cloned().collect::<BTreeSet<_>>(),
            two.cuts()
                .iter()
                .chain(three.cuts())
                .cloned()
                .collect::<BTreeSet<_>>(),
        );
        assert_eq!(both.loop_momentum_basis(), two.loop_momentum_basis());
        assert_eq!(both.loop_momentum_basis(), three.loop_momentum_basis());
        assert_eq!(
            both.topology_threshold_candidates(),
            two.topology_threshold_candidates()
        );
        assert_eq!(
            both.topology_threshold_candidates(),
            three.topology_threshold_candidates()
        );
        assert_eq!(both.numerator(), two.numerator());
        assert!(load("scalar_1,scalar_1|").is_err());
    }

    #[test]
    fn compact_self_energies_between_attachments_keep_their_cut_on_either_side() {
        let model = Arc::new(
            Model::from_json(include_str!(
                "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
            ))
            .unwrap(),
        );
        // Attachments a, b (incoming) and c, d (outgoing) on one loop
        // a - c - d - b, closed by a self-energy e = f on the internal line
        // between either the incoming or the outgoing pair.
        let incoming_side = r#"digraph incoming_side {
            final_state="scalar_1,scalar_1";
            in0 [style=invis]; in1 [style=invis]; out0 [style=invis]; out1 [style=invis];
            in0 -> a [id=0 particle=scalar_1 is_cut=7];
            in1 -> b [id=1 particle=scalar_1 is_cut=8];
            c -> out0 [id=2 particle=scalar_1 is_cut=7];
            d -> out1 [id=3 particle=scalar_1 is_cut=8];
            a -> c [id=4 particle=scalar_1];
            c -> d [id=5 particle=scalar_1];
            d -> b [id=6 particle=scalar_1];
            a -> e [id=7 particle=scalar_1];
            e -> f [id=8 particle=scalar_1];
            e -> f [id=9 particle=scalar_1];
            f -> b [id=10 particle=scalar_1];
        }"#;
        let outgoing_side = r#"digraph outgoing_side {
            final_state="scalar_1,scalar_1";
            in0 [style=invis]; in1 [style=invis]; out0 [style=invis]; out1 [style=invis];
            in0 -> a [id=0 particle=scalar_1 is_cut=7];
            in1 -> b [id=1 particle=scalar_1 is_cut=8];
            c -> out0 [id=2 particle=scalar_1 is_cut=7];
            d -> out1 [id=3 particle=scalar_1 is_cut=8];
            a -> c [id=4 particle=scalar_1];
            a -> b [id=5 particle=scalar_1];
            b -> d [id=6 particle=scalar_1];
            c -> e [id=7 particle=scalar_1];
            e -> f [id=8 particle=scalar_1];
            e -> f [id=9 particle=scalar_1];
            f -> d [id=10 particle=scalar_1];
        }"#;
        for input in [incoming_side, outgoing_side] {
            let (diagram, _) = FeynmanDiagram::from_compact_dot(
                Arc::clone(&model),
                DotGraph::from_string(input).unwrap(),
            )
            .unwrap();
            assert_eq!(diagram.cuts().len(), 1, "{input}");
        }
    }

    #[test]
    fn compact_multi_blob_cuts_follow_the_declared_blob_range() {
        let model = Arc::new(
            Model::from_json(include_str!(
                "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
            ))
            .unwrap(),
        );
        // The incoming-side self-energy graph. With only a and b on the left,
        // the partition crosses four lines and has two blobs on each side.
        let input = r#"digraph incoming_side {
            final_state="STATE"; BLOBS
            in0 [style=invis]; in1 [style=invis]; out0 [style=invis]; out1 [style=invis];
            in0 -> a [id=0 particle=scalar_1 is_cut=7];
            in1 -> b [id=1 particle=scalar_1 is_cut=8];
            c -> out0 [id=2 particle=scalar_1 is_cut=7];
            d -> out1 [id=3 particle=scalar_1 is_cut=8];
            a -> c [id=4 particle=scalar_1];
            c -> d [id=5 particle=scalar_1];
            d -> b [id=6 particle=scalar_1];
            a -> e [id=7 particle=scalar_1];
            e -> f [id=8 particle=scalar_1];
            e -> f [id=9 particle=scalar_1];
            f -> b [id=10 particle=scalar_1];
        }"#;
        let load = |state: &str, blobs: &str| {
            FeynmanDiagram::from_compact_dot(
                Arc::clone(&model),
                DotGraph::from_string(input.replace("STATE", state).replace("BLOBS", blobs))
                    .unwrap(),
            )
            .map(|(diagram, _)| diagram)
        };
        let four = "scalar_1,scalar_1,scalar_1,scalar_1";
        assert!(matches!(
            load(four, ""),
            Err(DiagramError::InvalidDotAttribute {
                attribute: "final_state",
                ..
            })
        ));
        let two_blobs = load(four, r#"n_cut_blobs="1,2";"#).unwrap();
        assert_eq!(two_blobs.cuts().len(), 1);
        assert_eq!(two_blobs.cuts()[0].cut.len(), 4);
        // Threshold candidates need connected sides whatever the blob range:
        // only the two-body boundary qualifies.
        let two_body = load("scalar_1,scalar_1", "").unwrap();
        assert_eq!(two_blobs.topology_threshold_candidates().len(), 1);
        assert_eq!(
            two_blobs.topology_threshold_candidates(),
            two_body.topology_threshold_candidates()
        );
        // The same ordered bounds as generation's --n-cut-blobs.
        for blobs in ["1,2", "0,1", "1"] {
            assert_eq!(
                load("scalar_1,scalar_1", &format!(r#"n_cut_blobs="{blobs}";"#))
                    .unwrap()
                    .cuts(),
                two_body.cuts(),
                "{blobs}"
            );
        }
        for invalid in ["2,1", "1,2,3", "one", "-1"] {
            assert!(matches!(
                load(four, &format!(r#"n_cut_blobs="{invalid}";"#)),
                Err(DiagramError::InvalidDotAttribute {
                    attribute: "n_cut_blobs",
                    ..
                })
            ));
        }
        // Canonical DOT keeps the resolved inventories, not the input filter,
        // and reproduces itself after import.
        let canonical = two_blobs.to_dot().unwrap();
        assert!(!canonical.contains("n_cut_blobs"));
        let reimported = FeynmanDiagram::from_dot(Arc::clone(&model), &canonical).unwrap();
        assert_eq!(reimported.cuts(), two_blobs.cuts());
        assert_eq!(
            reimported.topology_threshold_candidates(),
            two_blobs.topology_threshold_candidates()
        );
        assert_eq!(reimported.to_dot().unwrap(), canonical);
    }

    #[test]
    fn compact_cuts_keep_whole_vertices_behind_bridges() {
        let model = Arc::new(
            Model::from_json(include_str!(
                "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
            ))
            .unwrap(),
        );
        // The line a - x - c, the line b - d, a self-energy y = w between c
        // and d, and a self-loop vertex z attached to x by the bridge x - z.
        // The attachment search finds both two-body cuts by itself; behind the
        // cut {x-c, b-d}, z and its self-loop stay with x on the left.
        let input = r#"digraph attachment_bridge {
            final_state="scalar_1,scalar_1";
            in0 [style=invis]; in1 [style=invis]; out0 [style=invis]; out1 [style=invis];
            in0 -> a [id=0 particle=scalar_1 is_cut=7];
            in1 -> b [id=1 particle=scalar_1 is_cut=8];
            c -> out0 [id=2 particle=scalar_1 is_cut=7];
            d -> out1 [id=3 particle=scalar_1 is_cut=8];
            a -> b [id=4 particle=scalar_1];
            a -> x [id=5 particle=scalar_1];
            x -> c [id=6 particle=scalar_1];
            b -> d [id=7 particle=scalar_1];
            c -> y [id=8 particle=scalar_1];
            y -> w [id=9 particle=scalar_1];
            y -> w [id=10 particle=scalar_1];
            w -> d [id=11 particle=scalar_1];
            x -> z [id=12 particle=scalar_1];
            z -> z [id=13 particle=scalar_1];
        }"#;
        let (diagram, _) = FeynmanDiagram::from_compact_dot(
            Arc::clone(&model),
            DotGraph::from_string(input).unwrap(),
        )
        .unwrap();
        let loops = diagram
            .cuts()
            .iter()
            .map(|cut| [cut.left.loop_count, cut.right.loop_count])
            .collect::<BTreeSet<_>>();
        assert_eq!(loops, BTreeSet::from([[0, 2], [1, 1]]));
    }

    #[test]
    fn compact_dot_source_map_records_sewing_and_loop_ownership() {
        let model = Arc::new(
            Model::from_json(include_str!(
                "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
            ))
            .unwrap(),
        );
        let input = r#"digraph mapped_bubble {
            final_state="scalar_1,scalar_1"
            ext [style=invis]
            ext -> a [id=0 particle=scalar_1 is_cut=7]
            a -> b [id=1 particle=scalar_1 lmb_id=0]
            b -> ext [id=2 particle=scalar_1 is_cut=7]
            a -> b [id=3 particle=scalar_1]
        }"#;
        let parsed = DotGraph::from_string(input).unwrap();
        let (diagram, source) = FeynmanDiagram::from_compact_dot(model, parsed).unwrap();
        assert_eq!(
            source.edges(),
            &[EdgeId(0), EdgeId(1), EdgeId(0), EdgeId(2)]
        );
        assert_eq!(
            &source.vertices()[..2],
            &[Some(VertexId(0)), Some(VertexId(1))]
        );
        assert_eq!(
            diagram.loop_momentum_basis().loop_edges,
            [source.edges()[1]]
        );
        assert_eq!(diagram.edges().count(), 3);
        assert!(!diagram.cuts().is_empty());
        for cut in diagram.cuts() {
            let edges = cut
                .cut
                .iter()
                .map(|half| half.edge)
                .collect::<BTreeSet<_>>();
            assert_eq!(
                edges,
                BTreeSet::from([source.edges()[1], source.edges()[3]])
            );
        }
    }
}
