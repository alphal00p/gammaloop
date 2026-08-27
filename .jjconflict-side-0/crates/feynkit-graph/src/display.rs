use std::{collections::BTreeMap, fmt::Write};

use crate::{DiagramVertex, ExternalState, FeynmanDiagram, VertexId};

const EXTERNAL_Y_SCALE: f64 = 10.0;

fn typst_string(value: &str) -> String {
    let mut output = String::with_capacity(value.len() + 2);
    output.push('"');
    for character in value.chars() {
        match character {
            '\\' => output.push_str("\\\\"),
            '"' => output.push_str("\\\""),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            character if character.is_control() => {
                write!(output, "\\u{{{:x}}}", character as u32)
                    .expect("writing to a string cannot fail");
            }
            character => output.push(character),
        }
    }
    output.push('"');
    output
}

fn typst_optional_string(value: Option<&str>) -> String {
    value.map_or_else(|| "none".to_owned(), typst_string)
}

fn centered_external_y(count: usize, rank: usize) -> f64 {
    ((count as f64 - 1.0) / 2.0 - rank as f64) * EXTERNAL_Y_SCALE
}

fn external_vertices<'a>(
    vertices: &'a BTreeMap<VertexId, &'a DiagramVertex>,
    state: ExternalState,
) -> Vec<(VertexId, &'a DiagramVertex)> {
    vertices
        .iter()
        .filter_map(|(id, vertex)| {
            vertex
                .external
                .as_ref()
                .is_some_and(|external| external.state == state)
                .then_some((*id, *vertex))
        })
        .collect()
}

impl FeynmanDiagram {
    /// Emit a complete Typst document that renders this diagram with Linnest.
    ///
    /// The document imports the canonical Linnest package tree at
    /// `crates/linnest/typst`, which must be available below the Typst project
    /// root together with its sibling Kurvst package. Interaction vertices are
    /// densely remapped to valid Linnest node identifiers. FeynKit external
    /// vertices become Linnest dangling half-edges and retain their names,
    /// indices, and incoming/outgoing states as edge data.
    ///
    /// Layout uses the deterministic seed and force-layout settings from
    /// GammaLoop's amplitude renderer. Incoming and outgoing legs are constrained
    /// to the left and right respectively, while force layout separates loops and
    /// parallel propagators.
    pub fn to_linnest(&self) -> String {
        let vertices: BTreeMap<_, _> = self.vertices().collect();
        let internal_vertices: Vec<_> = vertices
            .iter()
            .filter_map(|(id, vertex)| (!vertex.is_external()).then_some((*id, *vertex)))
            .collect();
        let internal_ids: BTreeMap<_, _> = internal_vertices
            .iter()
            .enumerate()
            .map(|(dense_id, (id, _))| (*id, dense_id))
            .collect();

        let mut incoming = external_vertices(&vertices, ExternalState::Incoming);
        let mut outgoing = external_vertices(&vertices, ExternalState::Outgoing);
        incoming.sort_by_key(|(id, vertex)| {
            (
                vertex
                    .external
                    .as_ref()
                    .expect("external vertex has external metadata")
                    .index,
                id.0,
            )
        });
        outgoing.sort_by_key(|(id, vertex)| {
            (
                vertex
                    .external
                    .as_ref()
                    .expect("external vertex has external metadata")
                    .index,
                id.0,
            )
        });
        let incoming_ranks: BTreeMap<_, _> = incoming
            .iter()
            .enumerate()
            .map(|(rank, (id, _))| (*id, rank))
            .collect();
        let outgoing_ranks: BTreeMap<_, _> = outgoing
            .iter()
            .enumerate()
            .map(|(rank, (id, _))| (*id, rank))
            .collect();

        let mut output = String::from(
            r##"#set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm))
#set text(size: 9pt)
#import "crates/linnest/typst/src/draw.typ": draw
#import "crates/linnest/typst/src/graph.typ" as graph
#import "crates/linnest/typst/src/layout.typ" as layout
#import graph: build, edge, node, sink, source

#context {
  let raw = build({
"##,
        );

        for (id, vertex) in &internal_vertices {
            let dense_id = internal_ids[id];
            writeln!(
                output,
                "    node(<v{dense_id}>, id: {dense_id}, label: none, feynkit-id: {}, feynkit-name: {}, interaction: {}, numerator: {})",
                id.0,
                typst_string(&vertex.name),
                typst_optional_string(vertex.interaction.as_deref()),
                typst_optional_string(vertex.numerator.as_deref()),
            )
            .expect("writing to a string cannot fail");
        }

        for (id, endpoints, edge) in self.edges() {
            let source_internal = internal_ids.get(&endpoints.source).copied();
            let target_internal = internal_ids.get(&endpoints.target).copied();
            match (source_internal, target_internal) {
                (Some(source), Some(target)) => {
                    let orientation = if edge.directed {
                        "default"
                    } else {
                        "undirected"
                    };
                    writeln!(
                        output,
                        "    edge(source(<v{source}>), <e{}>, sink(<v{target}>), id: {}, orientation: {orientation:?}, label: text({}), particle: {}, pdg: {}, directed: {}, numerator: {}, feynkit-source: {}, feynkit-target: {})",
                        id.0,
                        id.0,
                        typst_string(&edge.particle.name),
                        typst_string(&edge.particle.name),
                        edge.particle.pdg,
                        edge.directed,
                        typst_optional_string(edge.numerator.as_deref()),
                        endpoints.source.0,
                        endpoints.target.0,
                    )
                    .expect("writing to a string cannot fail");
                }
                (Some(internal), None) | (None, Some(internal)) => {
                    let external_id = if source_internal.is_none() {
                        endpoints.source
                    } else {
                        endpoints.target
                    };
                    let external_vertex = vertices[&external_id];
                    let external = external_vertex
                        .external
                        .as_ref()
                        .expect("an edge with one internal endpoint has one external endpoint");
                    let original_external_is_source = endpoints.source == external_id;
                    let represented_external_is_source = external.state == ExternalState::Incoming;
                    let orientation = if !edge.directed {
                        "undirected"
                    } else if original_external_is_source == represented_external_is_source {
                        "default"
                    } else {
                        "reversed"
                    };
                    let (rank, count, side, side_constraint, label_anchor) = match external.state {
                        ExternalState::Incoming => (
                            incoming_ranks[&external_id],
                            incoming.len(),
                            "left",
                            "-",
                            "east",
                        ),
                        ExternalState::Outgoing => (
                            outgoing_ranks[&external_id],
                            outgoing.len(),
                            "right",
                            "+",
                            "west",
                        ),
                    };
                    let y = centered_external_y(count, rank);
                    let display_label =
                        format!("{} ({})", edge.particle.name, external_vertex.name);
                    let endpoint_spec = match external.state {
                        ExternalState::Incoming => format!("<e{}>, sink(<v{internal}>)", id.0),
                        ExternalState::Outgoing => format!("source(<v{internal}>), <e{}>", id.0),
                    };
                    writeln!(
                        output,
                        "    edge({endpoint_spec}, id: {}, orientation: {orientation:?}, label: text({}), label-anchor: {label_anchor:?}, particle: {}, pdg: {}, directed: {}, numerator: {}, external-state: {:?}, external-index: {}, external-name: {}, feynkit-source: {}, feynkit-target: {}, pos: graph.pos(x: graph.group({side:?}, side: {side_constraint:?}), y: graph.start({y:.1})))",
                        id.0,
                        typst_string(&display_label),
                        typst_string(&edge.particle.name),
                        edge.particle.pdg,
                        edge.directed,
                        typst_optional_string(edge.numerator.as_deref()),
                        external.state.as_str(),
                        external.index,
                        typst_string(&external_vertex.name),
                        endpoints.source.0,
                        endpoints.target.0,
                    )
                    .expect("writing to a string cannot fail");
                }
                (None, None) => {
                    unreachable!("validated diagrams cannot connect two external vertices")
                }
            }
        }

        let roots = if internal_vertices.is_empty() {
            "()"
        } else {
            "(0,)"
        };
        writeln!(
            output,
            "  }}, name: {}, data: (symmetry-factor: {}, overall-factor: {}, numerator: {}, loop-count: {}))",
            typst_string(self.name()),
            self.symmetry_factor(),
            typst_string(self.overall_factor()),
            typst_optional_string(self.numerator()),
            self.loop_count(),
        )
        .expect("writing to a string cannot fail");
        write!(
            output,
            r##"
  let edge-label(edge) = edge.data.at("label", default: none)
  let edge-label-style(edge) = (
    anchor: edge.data.at("label-anchor", default: "south"),
    padding: 0.08,
  )
  let styled = graph.style(
    raw,
    unit: 1.5,
    node-label: none,
    edge-label: edge-label,
    edge-label-style: edge-label-style,
  )
  let positioned = layout.layout(
    styled,
    seed: 2,
    steps: 30,
    epochs: 30,
    step: 0.6,
    k-spring: 4.5,
    eps: 1e-7,
    gamma-dangling: 2.3,
    directional-force: 4.5,
    label-length-scale: 1.2,
    label-steps: 100,
    label-layout: "dangling-tangent",
    layout-algo: "force",
    layout-direction: "right",
    layout-roots: {roots},
  )
  let directed-edge-style(edge) = (
    stroke: (paint: rgb("#315b8a"), thickness: 0.75pt, cap: "round"),
    mark: (
      end: (
        symbol: ">",
        fill: rgb("#315b8a"),
        anchor: "center",
        shorten-to: auto,
      ),
      scale: 0.7,
    ),
    mark-position: "center-if-dangling",
    mark-orientation: "edge",
  )
  draw(
    positioned,
    unit: 1.5,
    title: auto,
    node-label: none,
    node-radius: 0.10,
    node-fill: black,
    node-stroke: black,
    source-style: directed-edge-style,
    sink-style: directed-edge-style,
    edge-label: edge-label,
    edge-label-style: edge-label-style,
    padding: 0.55,
  )
}}
"##,
        )
        .expect("writing to a string cannot fail");
        output
    }
}

#[cfg(test)]
mod tests {
    use crate::{DiagramEdge, DiagramVertex, ExternalState, FeynmanDiagram, ParticleReference};

    fn one_loop() -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder("bubble");
        let incoming =
            builder.add_vertex(DiagramVertex::external("p1", 0, ExternalState::Incoming));
        let outgoing =
            builder.add_vertex(DiagramVertex::external("p2", 1, ExternalState::Outgoing));
        let left = builder.add_vertex(DiagramVertex::interaction("left", "V_1"));
        let right = builder.add_vertex(DiagramVertex::interaction("right", "V_1"));
        let scalar = || DiagramEdge::new(ParticleReference::new("phi", 25), false);
        builder.add_edge(incoming, left, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(right, outgoing, scalar()).unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn emits_deterministic_complete_linnest_source_for_a_loop() {
        let diagram = one_loop();
        let source = diagram.to_linnest();

        assert_eq!(source, diagram.to_linnest());
        assert!(source.starts_with("#set page(width: auto"));
        assert!(source.contains("#import \"crates/linnest/typst/src/draw.typ\": draw"));
        assert!(source.contains("#import \"crates/linnest/typst/src/graph.typ\" as graph"));
        assert!(source.contains("#import \"crates/linnest/typst/src/layout.typ\" as layout"));
        assert_eq!(source.matches("    node(").count(), 2);
        assert!(source.contains("node(<v0>, id: 0"));
        assert!(source.contains("node(<v1>, id: 1"));
        assert!(source.contains("edge(<e0>, sink(<v0>), id: 0"));
        assert!(source.contains("edge(source(<v1>), <e3>, id: 3"));
        assert_eq!(source.matches("edge(source(<v0>), <e").count(), 2);
        assert!(source.contains("graph.group(\"left\", side: \"-\")"));
        assert!(source.contains("graph.group(\"right\", side: \"+\")"));
        assert!(source.contains("layout-algo: \"force\""));
        assert!(source.contains("layout.layout("));
        assert!(source.contains("k-spring: 4.5"));
        assert!(source.contains("label-layout: \"dangling-tangent\""));
        assert!(source.contains("mark-orientation: \"edge\""));
        assert!(source.ends_with("}\n"));
    }

    #[test]
    fn sorts_external_legs_and_centers_each_side() {
        let mut builder = FeynmanDiagram::builder("one-to-two");
        let incoming =
            builder.add_vertex(DiagramVertex::external("in", 0, ExternalState::Incoming));
        let outgoing_high = builder.add_vertex(DiagramVertex::external(
            "out-high",
            2,
            ExternalState::Outgoing,
        ));
        let outgoing_low = builder.add_vertex(DiagramVertex::external(
            "out-low",
            1,
            ExternalState::Outgoing,
        ));
        let interaction = builder.add_vertex(DiagramVertex::interaction("v", "V_3"));
        let scalar = || DiagramEdge::new(ParticleReference::new("phi", 25), false);
        builder.add_edge(incoming, interaction, scalar()).unwrap();
        builder
            .add_edge(interaction, outgoing_high, scalar())
            .unwrap();
        builder
            .add_edge(interaction, outgoing_low, scalar())
            .unwrap();
        let source = builder.build().unwrap().to_linnest();

        let low = source.find("external-name: \"out-low\"").unwrap();
        let high = source.find("external-name: \"out-high\"").unwrap();
        assert!(low > high, "edge emission remains stable by edge id");
        assert!(source.contains("external-name: \"out-high\", feynkit-source: 3, feynkit-target: 1, pos: graph.pos(x: graph.group(\"right\", side: \"+\"), y: graph.start(-5.0))"));
        assert!(source.contains("external-name: \"out-low\", feynkit-source: 3, feynkit-target: 2, pos: graph.pos(x: graph.group(\"right\", side: \"+\"), y: graph.start(5.0))"));
    }

    #[test]
    fn escapes_typst_strings_and_preserves_directed_orientation() {
        let mut builder = FeynmanDiagram::builder("quote \" and \\ slash");
        let incoming =
            builder.add_vertex(DiagramVertex::external("p\n1", 0, ExternalState::Incoming));
        let interaction = builder.add_vertex(DiagramVertex::interaction("v", "V\"1"));
        builder
            .add_edge(
                interaction,
                incoming,
                DiagramEdge::new(ParticleReference::new("f\\\"", 1), true),
            )
            .unwrap();
        let source = builder.build().unwrap().to_linnest();

        assert!(source.contains("name: \"quote \\\" and \\\\ slash\""));
        assert!(source.contains("feynkit-name: \"v\""));
        assert!(source.contains("interaction: \"V\\\"1\""));
        assert!(source.contains("external-name: \"p\\n1\""));
        assert!(source.contains("orientation: \"reversed\""));
        assert!(!source.contains("p\n1"));
    }
}
