//! Construct an uncut model-aware topology without accepting saved-diagram metadata.

use super::*;

impl FeynmanDiagram {
    /// Construct an uncut topology whose edges name model particles or PDG codes.
    ///
    /// Invisible DOT nodes denote external half-edges; each becomes a distinct
    /// external vertex. `lmb_id` selects an ordered loop basis. Graph `num`,
    /// `overall_factor` and `projector`, and node/edge `num` remain explicit
    /// factors; interaction numerators are never inserted. Unique matching UFO
    /// interactions supply vertex-slot order, while unmatched structural vertices
    /// remain untyped. Ambiguous or explicit interactions must match `int_id`.
    /// Saved diagrams must use [`Self::from_dot`], which verifies their model
    /// fingerprint and serialized state. Cut sewing is not inferred here.
    pub fn from_topology_dot(
        model: impl Into<Arc<Model>>,
        input: &str,
    ) -> Result<Self, DiagramError> {
        let model = model.into();
        let parsed: DotGraph = DotGraph::from_string(input)
            .map_err(|error| DiagramError::DotParse(error.to_string()))?;
        let attributes = &parsed.global_data.statements;
        let invalid = |target: String, attribute, value: &str| DiagramError::InvalidDotAttribute {
            target,
            attribute,
            value: value.into(),
        };
        for key in [
            "model_fingerprint",
            "feynkit_name",
            "half_edge_order",
            "loop_momentum_basis",
            "loop_momentum_edges",
            "cuts",
            "topology_threshold_candidates",
            "numerator_prefactor",
            "numerator",
            "symmetry_factor",
        ] {
            if attributes.contains_key(key) {
                return Err(invalid(
                    "graph".into(),
                    key,
                    "saved diagram metadata requires from_dot",
                ));
            }
        }
        if attributes.contains_key("final_state") {
            return Err(invalid(
                "graph".into(),
                "final_state",
                "topology construction does not infer cuts",
            ));
        }
        let mut builder = Self::builder(Arc::clone(&model), parsed.global_data.name.clone());
        for (name, destination) in [
            ("num", &mut builder.numerator_prefactor),
            ("overall_factor", &mut builder.overall_factor),
            ("projector", &mut builder.projector),
        ] {
            if let Some(value) = attributes.get(name) {
                *destination = Self::parse_expression(name, value.clone())?;
            }
        }
        let mut nodes = BTreeMap::new();
        for (node, _, data) in parsed.iter_nodes() {
            for key in [
                "vertex_id",
                "external_index",
                "external_state",
                "external_connection",
                "interaction",
                "interaction_id",
                "feynkit_name",
                "numerator",
            ] {
                if data.statements.contains_key(key) {
                    return Err(invalid(
                        format!("vertex {}", node.0),
                        key,
                        "saved diagram metadata requires from_dot",
                    ));
                }
            }
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
                    external: None,
                    numerator,
                }),
            );
        }
        let mut loop_edges = BTreeMap::new();
        let mut next_external = 0usize;
        for (pair, edge_id, data) in parsed.iter_edges() {
            let attributes = &data.data.statements;
            let target = format!("edge {}", edge_id.0);
            for key in [
                "particle_id",
                "source_slot",
                "target_slot",
                "directed",
                "numerator",
            ] {
                if attributes.contains_key(key) {
                    return Err(invalid(
                        target,
                        key,
                        "saved diagram metadata requires from_dot",
                    ));
                }
            }
            if attributes.contains_key("is_cut") {
                return Err(invalid(
                    target,
                    "is_cut",
                    "topology construction does not infer cuts",
                ));
            }
            let by_name = attributes
                .get("particle")
                .map(|name| model.particle_id(name))
                .transpose()?;
            let by_pdg = attributes
                .get("pdg")
                .map(|value| {
                    let pdg = value
                        .parse::<i64>()
                        .map_err(|_| invalid(target.clone(), "pdg", value))?;
                    model.particle_id_by_pdg(pdg).map_err(DiagramError::from)
                })
                .transpose()?;
            if by_name.is_some() && by_pdg.is_some() && by_name != by_pdg {
                return Err(invalid(
                    target,
                    "pdg",
                    "particle and pdg refer to different model particles",
                ));
            }
            let particle = by_name
                .or(by_pdg)
                .ok_or_else(|| DiagramError::MissingDotAttribute {
                    target: target.clone(),
                    attribute: "particle or pdg",
                })?;
            let directed = !model.particle_is_self_conjugate(particle);
            if attributes.contains_key("dir")
                && (data.orientation != Orientation::Undirected) != directed
            {
                return Err(invalid(
                    target,
                    "dir",
                    "orientation disagrees with the model particle",
                ));
            }
            let mut edge = DiagramEdge::new(particle, directed);
            edge.numerator = Self::parse_expression(
                "edge numerator",
                attributes.get("num").cloned().unwrap_or_else(|| "1".into()),
            )?;
            builder.numerator *= &edge.numerator;
            let internal = matches!(pair, HedgePair::Paired { .. } | HedgePair::Split { .. });
            let (source, sink) = match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    let mut source = nodes[&parsed.node_id(source)];
                    let mut sink = nodes[&parsed.node_id(sink)];
                    if data.orientation == Orientation::Reversed {
                        std::mem::swap(&mut source, &mut sink);
                    }
                    (source, sink)
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let index = next_external;
                    next_external = next_external
                        .checked_add(1)
                        .ok_or(DiagramError::DotExternalIndexOverflow(index))?;
                    let state = match flow {
                        Flow::Sink => ExternalState::Incoming,
                        Flow::Source => ExternalState::Outgoing,
                    };
                    let name = attributes
                        .get("name")
                        .cloned()
                        .unwrap_or_else(|| format!("ext{index}"));
                    let external = builder.add_vertex(DiagramVertex::external(name, index, state));
                    let node = nodes[&parsed.node_id(hedge)];
                    match flow {
                        Flow::Sink => (external, node),
                        Flow::Source => (node, external),
                    }
                }
            };
            let added = builder.add_edge(source, sink, edge)?;
            if let Some(value) = attributes.get("lmb_id") {
                let slot = value
                    .parse::<usize>()
                    .map_err(|_| invalid(target.clone(), "lmb_id", value))?;
                if !internal || loop_edges.insert(slot, added).is_some() {
                    return Err(invalid(
                        target,
                        "lmb_id",
                        "loop slots must be unique and belong to internal edges",
                    ));
                }
            }
        }
        for (node, _, data) in parsed.iter_nodes() {
            let vertex = nodes[&node];
            let mut incident = Vec::new();
            for (edge, (endpoints, data)) in builder.edges.iter().enumerate() {
                if endpoints.source == vertex {
                    incident.push((
                        edge,
                        true,
                        model.particle_by_id(data.particle)?.antiparticle,
                    ));
                }
                if endpoints.target == vertex {
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
                    (particles == signature && explicit.is_none_or(|chosen| chosen.index() == id))
                        .then_some(id)
                })
                .collect::<Vec<_>>();
            if matches.is_empty() && explicit.is_none() {
                // A propagator topology need not carry a model interaction at
                // every vertex (for example a one-line vacuum tadpole).
                continue;
            }
            let [rule] = matches.as_slice() else {
                return Err(invalid(
                    format!("vertex {}", builder.vertices[vertex.0].name),
                    "int_id",
                    "incident particles need a unique matching interaction",
                ));
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
        let mut diagram = builder.build()?;
        if !loop_edges.is_empty() {
            if loop_edges.keys().copied().ne(0..loop_edges.len()) {
                return Err(DiagramError::InvalidLoopMomentumBasis(
                    "lmb_id slots must be contiguous from zero".into(),
                ));
            }
            diagram =
                diagram.with_loop_momentum_edges(&loop_edges.into_values().collect::<Vec<_>>())?;
        }
        diagram.validate()?;
        Ok(diagram)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compact_vacuum_and_external_topologies_preserve_factors_and_basis() {
        let model = crate::tests::scalar_model();
        let tadpole = FeynmanDiagram::from_topology_dot(
            model.clone(),
            r#"digraph { a -> a [particle="phi",lmb_id=0]; }"#,
        )
        .unwrap();
        assert_eq!(tadpole.loop_count(), 1);
        assert_eq!(tadpole.edges().count(), 1);
        let graph = FeynmanDiagram::from_topology_dot(model.clone(), r#"digraph triangle {
            graph [num="5"]; ext [style=invis]; a [num="2"];
            ext -> a [id=0,pdg=25]; b -> ext [id=1,particle="phi"]; c -> ext [id=2,particle="phi"];
            a -> b [id=3,particle="phi",num="3",lmb_id=0]; b -> c [id=4,particle="phi"]; c -> a [id=5,particle="phi"];
        }"#).unwrap();
        assert_eq!(graph.loop_count(), 1);
        assert_eq!(graph.vertices().filter(|(_, v)| v.is_external()).count(), 3);
        assert_eq!(graph.numerator(), &Atom::num(6));
        assert_eq!(graph.numerator_prefactor(), &Atom::num(5));
        assert_eq!(graph.loop_momentum_basis().loop_edges, vec![EdgeId(3)]);
        let saved = graph.to_dot().unwrap();
        let restored = FeynmanDiagram::from_dot(model.clone(), &saved).unwrap();
        assert_eq!(restored.to_json().unwrap(), graph.to_json().unwrap());
        assert!(FeynmanDiagram::from_topology_dot(model, &saved).is_err());
    }

    #[test]
    fn compact_import_rejects_inconsistent_or_serialized_state() {
        let model = crate::tests::scalar_model();
        for source in [
            r#"digraph { graph [cuts="[]"]; a -> a [particle="phi"]; }"#,
            r#"digraph { a -> a [particle="phi",particle_id=0]; }"#,
            r#"digraph { a -> a [particle="phi",pdg=999]; }"#,
            r#"digraph { a -> a [particle="phi",lmb_id=1]; }"#,
            r#"digraph { a -> a [particle="phi",is_cut=0]; }"#,
            r#"digraph { a [interaction_id=0]; a -> a [particle="phi"]; }"#,
            r#"digraph { a [external_connection=0]; a -> a [particle="phi"]; }"#,
            r#"digraph { a [feynkit_name="saved"]; a -> a [particle="phi"]; }"#,
            r#"digraph { ext [style=invis]; ext -> a [particle="phi",lmb_id=0]; a -> a [particle="phi"]; }"#,
        ] {
            assert!(
                FeynmanDiagram::from_topology_dot(model.clone(), source).is_err(),
                "{source}"
            );
        }
        assert!(
            FeynmanDiagram::from_dot(model, r#"digraph { a -> a [particle="phi"]; }"#).is_err()
        );
    }
}
