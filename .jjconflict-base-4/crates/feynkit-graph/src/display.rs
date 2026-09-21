use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Write,
};

use crate::FeynmanDiagram;
use linnet::half_edge::{
    EdgeAccessors,
    involution::{EdgeIndex, Orientation},
};
use symbolica::atom::AtomCore;

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

impl FeynmanDiagram {
    /// Emit a complete Typst document that renders this diagram with Linnest.
    ///
    /// The document imports the canonical Linnest package tree at
    /// `crates/linnest/typst`, which must be available below the Typst project
    /// root together with its sibling Kurvst package and the shared physics
    /// styles in `assets/embedded/drawing/templates`. Particle spin, color,
    /// charge, mass, and TeX names select line patterns, arrows, and labels.
    /// Interaction vertices retain their native identifiers. Amplitude external
    /// legs are dangling half-edges and retain their names, indices, and
    /// incoming/outgoing states as edge data.
    ///
    /// GammaLoop's shared physics layout owns particle styling, label measurement,
    /// force settings, and left/right amplitude placement. Finalized cross sections
    /// already contain paired initial-state edges with their sewing connection IDs.
    pub fn to_linnest(&self) -> String {
        let internal_vertices: Vec<_> = self.vertices().collect();
        let internal_ids: BTreeMap<_, _> = internal_vertices
            .iter()
            .map(|(id, _)| (*id, id.0))
            .collect();

        let mut output = String::from(
            r##"#set page(width: auto, height: auto, margin: (x: 2mm, y: 2mm), fill: none)
#set text(size: 9pt)
#import "crates/linnest/typst/src/graph.typ" as graph
#import "crates/linnest/typst/src/render/layout.typ" as renderer
#import "assets/embedded/drawing/templates/layout-core.typ" as physics-layout
#import "assets/embedded/drawing/templates/physics-edge-style.typ" as physics
#import physics: mi, palette, massive, massless, dashed, dotted, source-stroke, sink-stroke, fermion-flow, wave, coil, zigzag
#import graph: build, edge, node, sink, source
#set text(fill: palette.ink)

#let particle-map = (
"##,
        );

        let particles = self
            .edges()
            .map(|(_, _, edge)| edge.particle)
            .collect::<BTreeSet<_>>();
        if particles.is_empty() {
            output.push(':');
        }
        for particle in particles {
            let particle = self
                .model()
                .particle_by_id(particle)
                .expect("validated particle ID");
            writeln!(
                output,
                "  {}: {},",
                typst_string(&particle.name),
                particle.generate_edge_typst_dict(self.model())
            )
            .expect("writing to a string cannot fail");
        }
        output.push_str(")\n\n#context {\n  let raw = build({\n");

        for (id, vertex) in &internal_vertices {
            let dense_id = internal_ids[id];
            let interaction = vertex
                .interaction
                .and_then(|rule| self.model().vertex_rule_by_id(rule).ok())
                .map(|rule| typst_string(&rule.name))
                .unwrap_or_else(|| "none".to_owned());
            writeln!(
                output,
                "    node(<v{dense_id}>, id: {dense_id}, label: none, feynkit-id: {}, feynkit-name: {}, interaction: {}, numerator: {})",
                id.0,
                typst_string(&vertex.name),
                interaction,
                typst_string(&vertex.numerator.to_canonical_string()),
            )
            .expect("writing to a string cannot fail");
        }

        for (id, endpoints, edge) in self.edges() {
            let particle = self
                .model()
                .particle_by_id(edge.particle)
                .expect("validated diagram particle IDs resolve in the owned model");
            let source_internal = endpoints.source.map(|vertex| internal_ids[&vertex]);
            let target_internal = endpoints.target.map(|vertex| internal_ids[&vertex]);
            let orientation = match self.underlying().orientation(EdgeIndex(id.0)) {
                Orientation::Default => "default",
                Orientation::Reversed => "reversed",
                Orientation::Undirected => "undirected",
            };
            let endpoint_spec = match (source_internal, target_internal) {
                (Some(source), Some(target)) => {
                    format!("source(<v{source}>), <e{}>, sink(<v{target}>)", id.0)
                }
                (Some(source), None) => format!("source(<v{source}>), <e{}>", id.0),
                (None, Some(target)) => format!("<e{}>, sink(<v{target}>)", id.0),
                (None, None) => unreachable!("every native edge has an incident interaction"),
            };
            let external_metadata = edge.external.as_ref().map(|external| format!(
                ", external-state: {:?}, external-index: {}, external-name: {}, external-connection: {}",
                external.state.as_str(), external.index, typst_string(&external.name), external.connection,
            )).unwrap_or_default();
            let cut_metadata = edge
                .external
                .as_ref()
                .filter(|_| endpoints.source.is_some() && endpoints.target.is_some())
                .map(|external| format!(", is_cut: {}", external.connection))
                .unwrap_or_default();
            writeln!(
                output,
                "    edge({endpoint_spec}, id: {}, orientation: {orientation:?}, particle: {}, pdg: {}, directed: {}, numerator: {}, feynkit-source: {}, feynkit-target: {}{external_metadata}{cut_metadata})",
                id.0,
                typst_string(&particle.name),
                particle.pdg_code,
                edge.directed,
                typst_string(&edge.numerator.to_canonical_string()),
                endpoints.source.map(|vertex| vertex.0.to_string()).unwrap_or_else(|| "none".to_owned()),
                endpoints.target.map(|vertex| vertex.0.to_string()).unwrap_or_else(|| "none".to_owned()),
            )
            .expect("writing to a string cannot fail");
        }

        writeln!(
            output,
            "  }}, name: {}, data: (symmetry-factor: {}, overall-factor: {}, numerator: {}, loop-count: {}))",
            typst_string(self.name()),
            self.symmetry_factor(),
            typst_string(&self.overall_factor().to_canonical_string()),
            typst_string(&self.numerator().to_canonical_string()),
            self.loop_count(),
        )
        .expect("writing to a string cannot fail");
        write!(
            output,
            r##"
  physics-layout.layout(
    raw,
    graph: graph,
    renderer: renderer,
    physics: physics,
    edge-style: (map: particle-map, default-edge: physics.default-edge),
    unit: 1.5,
    amplitude-mode: {},
    cross-section-mode: {},
    style-options: (node-label: none),
    edge-style-options: (label-fill: palette.ink),
  )
}}
"##,
            self.cuts().is_empty(),
            !self.cuts().is_empty(),
        )
        .expect("writing to a string cannot fail");
        output
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::{DiagramEdge, DiagramVertex, ExternalLeg, ExternalState, FeynmanDiagram};
    use feynkit_model::Model;

    fn display_model() -> Arc<Model> {
        Arc::new(Model::from_json(
            r#"{
                "name":"display","restriction":null,"orders":[],
                "parameters":[
                    {"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null},
                    {"name":"M","lhablock":"MASS","lhacode":[25],"nature":"external","parameter_type":"real","value":[1.0,0.0],"expression":null}
                ],
                "particles":[{"pdg_code":25,"name":"phi","antiname":"phi","spin":1,"color":1,"mass":"M","width":"ZERO","texname":"phi","antitexname":"phi","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0}],
                "propagators":[{"name":"phi_prop","particle":"phi","numerator":"1","denominator":"P^2-M^2"}],
                "lorentz_structures":[
                    {"name":"L3","spins":[1,1,1],"structure":"1"},
                    {"name":"L1","spins":[1],"structure":"1"}
                ],
                "couplings":[],
                "vertex_rules":[
                    {"name":"V_1","particles":["phi","phi","phi"],"color_structures":["1"],"lorentz_structures":["L3"],"couplings":[[null]]},
                    {"name":"V_3","particles":["phi","phi","phi"],"color_structures":["1"],"lorentz_structures":["L3"],"couplings":[[null]]},
                    {"name":"V\"1","particles":["phi"],"color_structures":["1"],"lorentz_structures":["L1"],"couplings":[[null]]}
                ]
            }"#,
        ).unwrap())
    }

    fn one_loop() -> FeynmanDiagram {
        let model = display_model();
        let rule = model.vertex_rule_id("V_1").unwrap();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "bubble");
        let mut incoming = DiagramEdge::new(particle, false);
        incoming.external = Some(ExternalLeg {
            name: "p1".into(),
            index: 0,
            state: ExternalState::Incoming,
            connection: 0,
        });
        let mut outgoing = DiagramEdge::new(particle, false);
        outgoing.external = Some(ExternalLeg {
            name: "p2".into(),
            index: 1,
            state: ExternalState::Outgoing,
            connection: 1,
        });
        let left = builder.add_vertex(DiagramVertex::interaction("left", rule));
        let right = builder.add_vertex(DiagramVertex::interaction("right", rule));
        let scalar = || DiagramEdge::new(particle, false);
        builder.add_edge(None, left, incoming).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(right, None, outgoing).unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn emits_deterministic_complete_linnest_source_for_a_loop() {
        let diagram = one_loop();
        let source = diagram.to_linnest();

        assert_eq!(source, diagram.to_linnest());
        assert!(source.starts_with("#set page(width: auto"));
        assert!(source.contains("#import \"crates/linnest/typst/src/graph.typ\" as graph"));
        assert_eq!(source.matches("    node(").count(), 2);
        assert!(source.contains("node(<v0>, id: 0"));
        assert!(source.contains("node(<v1>, id: 1"));
        assert!(source.contains("edge(<e0>, sink(<v0>), id: 0"));
        assert!(source.contains("edge(source(<v1>), <e3>, id: 3"));
        assert_eq!(source.matches("edge(source(<v0>), <e").count(), 2);
        assert!(source.contains("physics-layout.layout("));
        assert!(source.contains("amplitude-mode: true"));
        assert!(source.contains("cross-section-mode: false"));
        assert!(!source.contains("is_cut:"));
        assert!(!source.contains("pos: graph.pos"));
        assert!(source.contains("dash: dashed"));
        assert!(source.ends_with("}\n"));
    }

    #[test]
    fn delegates_external_placement_without_reordering_edges() {
        let model = display_model();
        let rule = model.vertex_rule_id("V_3").unwrap();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "one-to-two");
        let interaction = builder.add_vertex(DiagramVertex::interaction("v", rule));
        for (name, index, state) in [
            ("in", 0, ExternalState::Incoming),
            ("out-high", 2, ExternalState::Outgoing),
            ("out-low", 1, ExternalState::Outgoing),
        ] {
            let mut edge = DiagramEdge::new(particle, false);
            edge.external = Some(ExternalLeg {
                name: name.into(),
                index,
                state,
                connection: index,
            });
            match state {
                ExternalState::Incoming => builder.add_edge(None, interaction, edge),
                ExternalState::Outgoing => builder.add_edge(interaction, None, edge),
            }
            .unwrap();
        }
        let source = builder.build().unwrap().to_linnest();

        let low = source.find("external-name: \"out-low\"").unwrap();
        let high = source.find("external-name: \"out-high\"").unwrap();
        assert!(low > high, "edge emission remains stable by edge id");
        assert!(source.contains("amplitude-mode: true"));
        assert!(!source.contains("pos: graph.pos"));
    }

    #[test]
    fn escapes_typst_strings_and_preserves_directed_orientation() {
        let model = display_model();
        let rule = model.vertex_rule_id("V\"1").unwrap();
        let particle = model.particle_id("phi").unwrap();
        let mut builder = FeynmanDiagram::builder(model, "quote \" and \\ slash");
        let interaction = builder.add_vertex(DiagramVertex::interaction("v", rule));
        let mut incoming = DiagramEdge::new(particle, true);
        incoming.external = Some(ExternalLeg {
            name: "p\n1".into(),
            index: 0,
            state: ExternalState::Incoming,
            connection: 0,
        });
        builder.add_edge(None, interaction, incoming).unwrap();
        builder.edge_orientations =
            Some(vec![linnet::half_edge::involution::Orientation::Reversed]);
        let source = builder.build().unwrap().to_linnest();

        assert!(source.contains("name: \"quote \\\" and \\\\ slash\""));
        assert!(source.contains("feynkit-name: \"v\""));
        assert!(source.contains("interaction: \"V\\\"1\""));
        assert!(source.contains("external-name: \"p\\n1\""));
        assert!(source.contains("orientation: \"reversed\""));
        assert!(!source.contains("p\n1"));
    }
}
