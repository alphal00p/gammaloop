//! Native drawing of a diagram: the physics scene `to_linnest` describes in
//! Typst, laid out and painted by `linnest::svg`. Typst typesets only labels.
use std::collections::{BTreeMap, BTreeSet};

use feynkit_model::{LineKind, Model, Particle, ParticleId};
use linnest::{
    TypstEdgeSpec, TypstEndpointSpec, TypstGraphSpec, TypstNodeSpec,
    svg::{Config, Dash, Details, EdgeDrawing, NodeDrawing, Pattern, Scene, Stroke},
};
use linnet::half_edge::{
    EdgeAccessors,
    involution::{EdgeIndex, Orientation},
    subgraph::{SuBitGraph, SubSetLike},
};
use serde_json::{Map, Value};

use crate::{DiagramError, FeynmanDiagram, LoopMomentumBasis};

/// The physics palette of `physics-edge-style.typ` and the highlight paints.
const INK: &str = "#3d2645";
const ACCENT: &str = "#6f4d85";
const HIGHLIGHT: &str = "#ffd166";
const OUTSIDE: &str = "#77777773";
/// External legs start one layout group either side, ten units apart.
const EXTERNAL_SPACING: f64 = 10.0;

/// Physics drawing options the native renderer honours.
#[derive(Clone, Debug, PartialEq)]
pub struct SceneOptions {
    /// Momentum arrows beside every edge.
    pub momentum_arrows: bool,
    /// Momenta in labels; `None` follows the arrows.
    pub show_momentum: Option<bool>,
    pub show_particle: bool,
    pub show_edge_index: bool,
    pub show_node_index: bool,
    pub split_initial_state: bool,
    /// Refine the layout around drawn labels in short warm passes.
    pub label_feedback: bool,
    /// Explicit `impred-*` layout options.
    pub layout: Map<String, Value>,
}

impl Default for SceneOptions {
    fn default() -> Self {
        Self {
            momentum_arrows: false,
            show_momentum: None,
            show_particle: true,
            show_edge_index: false,
            show_node_index: false,
            split_initial_state: true,
            label_feedback: false,
            layout: Map::new(),
        }
    }
}

impl FeynmanDiagram {
    /// The native drawing of this diagram, as `to_linnest` renders it with
    /// Typst: particle lines and flow arrows, momentum arrows, highlighting and
    /// inspection details. Cross sections open their initial states before
    /// layout while retaining the original inspection identities.
    pub fn to_scene(
        &self,
        highlight: Option<&SuBitGraph>,
        isolated: &BTreeSet<usize>,
        lmb: Option<&LoopMomentumBasis>,
        options: &SceneOptions,
    ) -> Result<Scene, DiagramError> {
        let basis = lmb.unwrap_or_else(|| self.loop_momentum_basis());
        if lmb.is_some() {
            basis.validate(self)?;
        }
        let show_momentum = options.show_momentum.unwrap_or(options.momentum_arrows);
        let dense: BTreeMap<_, _> = self
            .vertices()
            .enumerate()
            .map(|(index, (id, _))| (id, index))
            .collect();

        // Highlighted edges carry any selected half; their vertices stay selected.
        let mut selected_edges = BTreeSet::new();
        let mut selected_nodes = isolated.clone();
        if let Some(highlight) = highlight {
            for hedge in highlight.included_iter() {
                selected_nodes.insert(self.graph.node_id(hedge).0);
                selected_edges.insert(self.graph[&hedge].0);
            }
        }

        let nodes = self
            .vertices()
            .map(|(id, vertex)| {
                let incident: Vec<Value> = self
                    .edges()
                    .filter(|(_, ends, _)| ends.source == Some(id) || ends.target == Some(id))
                    .map(|(edge, _, _)| edge.0.into())
                    .collect();
                let stroke = match highlight {
                    None => Stroke {
                        paint: INK.to_owned(),
                        width: 1.45,
                        dash: Dash::Solid,
                        round_cap: false,
                    },
                    Some(_) if selected_nodes.contains(&id.0) => Stroke {
                        paint: HIGHLIGHT.to_owned(),
                        width: 1.2,
                        dash: Dash::Solid,
                        round_cap: false,
                    },
                    Some(_) => Self::outside_stroke(),
                };
                let details = [
                    ("node", Value::from(id.0)),
                    ("edges", incident.into()),
                    ("name", vertex.name.clone().into()),
                ]
                .into_iter()
                .collect();
                (
                    TypstNodeSpec {
                        name: Some(format!("v{}", dense[&id])),
                        index: Some(dense[&id]),
                        data: None,
                        pos: None,
                        statements: [("layout-width", "0.16"), ("layout-height", "0.16")]
                            .map(|(key, value)| (key.to_owned(), value.to_owned()))
                            .into(),
                    },
                    NodeDrawing {
                        radius: 0.16,
                        label: None,
                        rectangular: false,
                        fill: "none".into(),
                        stroke,
                        details,
                    },
                )
            })
            .collect::<Vec<_>>();

        // GammaLoop's amplitude convention: incoming legs on the left and outgoing
        // legs on the right, each side in half-edge order with one free x group.
        let (incoming, outgoing): (Vec<_>, Vec<_>) = self
            .edges()
            .filter(|(_, ends, _)| ends.source.is_none() != ends.target.is_none())
            .partition(|(_, ends, _)| ends.source.is_none());
        let place = !incoming.is_empty() && !outgoing.is_empty();

        let mut pages: Vec<String> = Vec::new();
        let mut hedge = 0;
        let mut edges = Vec::new();
        for (id, ends, edge) in self.edges() {
            let particle = self
                .model()
                .particle_by_id(edge.particle)
                .expect("validated diagram particle IDs resolve in the owned model");
            let orientation = match self.underlying().orientation(EdgeIndex(id.0)) {
                Orientation::Default => "default",
                Orientation::Reversed => "reversed",
                Orientation::Undirected => "undirected",
            };
            // Builder half-edges follow source/sink input order.
            let mut endpoint = |vertex: Option<_>| {
                vertex.map(|vertex| {
                    hedge += 1;
                    (
                        hedge - 1,
                        TypstEndpointSpec {
                            node: dense[&vertex],
                            statement: None,
                            id: None,
                            data: None,
                            port_label: None,
                            compass: None,
                            in_subgraph: false,
                            route_points: Vec::new(),
                        },
                    )
                })
            };
            let (source, sink) = (endpoint(ends.source), endpoint(ends.target));

            let mut statements = BTreeMap::new();
            if source.is_none() || sink.is_none() {
                statements.insert("pos-z".to_owned(), "0".to_owned());
                statements.insert("pos-z-mode".to_owned(), "pin".to_owned());
            }
            if place {
                let (side, legs, dx) = match (&source, &sink) {
                    (None, Some(_)) => ("-left", &incoming, -EXTERNAL_SPACING),
                    (Some(_), None) => ("+right", &outgoing, EXTERNAL_SPACING),
                    _ => ("", &incoming, 0.0),
                };
                if let Some(rank) = legs
                    .iter()
                    .position(|(leg, _, _)| *leg == id)
                    .filter(|_| dx != 0.0)
                {
                    let y = ((legs.len() as f64 - 1.0) / 2.0 - rank as f64) * EXTERNAL_SPACING;
                    for (key, value) in [
                        ("group-start-x", "true".to_owned()),
                        ("pin", format!("x:@{side}")),
                        ("pos", format!("{dx},{y}")),
                        ("pos-mode", "pin".to_owned()),
                        ("pos-x-set", "true".to_owned()),
                        ("pos-y-set", "true".to_owned()),
                    ] {
                        statements.insert(key.to_owned(), value);
                    }
                }
            }

            let drawing = options.particle_drawing(particle, self.model(), orientation);
            let stroke = match highlight {
                None => drawing.stroke,
                Some(_) if selected_edges.contains(&id.0) => Stroke {
                    paint: HIGHLIGHT.into(),
                    ..drawing.stroke
                },
                Some(_) => Self::outside_stroke(),
            };
            let momentum = &basis.edge_signatures[&id];
            let label =
                (options.show_particle || show_momentum || options.show_edge_index).then(|| {
                    let mut components = Vec::new();
                    if options.show_particle {
                        components.push(particle.typst_label());
                    }
                    if show_momentum {
                        let (loops, external) = momentum.integer_coefficients();
                        let mut terms = Vec::new();
                        for (head, coefficients) in [("k", loops), ("p", external)] {
                            for (index, coefficient) in coefficients.into_iter().enumerate() {
                                if coefficient == 0 {
                                    continue;
                                }
                                let sign = if coefficient < 0 {
                                    "-"
                                } else if terms.is_empty() {
                                    ""
                                } else {
                                    "+"
                                };
                                let factor = if coefficient.abs() == 1 {
                                    String::new()
                                } else {
                                    coefficient.abs().to_string()
                                };
                                terms.push(format!("{sign}{factor} {head}_{index}"));
                            }
                        }
                        components.push(if terms.is_empty() {
                            "0".into()
                        } else {
                            terms.join(" ")
                        });
                    }
                    if options.show_edge_index {
                        components.push(format!("upright(e)_{}", id.0));
                    }
                    let content = format!("[$ {} $]", components.join(" quad "));
                    pages
                        .iter()
                        .position(|page| *page == content)
                        .unwrap_or_else(|| {
                            pages.push(content);
                            pages.len() - 1
                        })
                });

            let vertex =
                |vertex: Option<crate::VertexId>| vertex.map_or(Value::Null, |v| v.0.into());
            let mut details: Details = [
                ("edge", Value::from(id.0)),
                ("name", Value::from(format!("e{}", id.0))),
                ("source", vertex(ends.source)),
                ("sink", vertex(ends.target)),
                (
                    "source-hedge",
                    source.as_ref().map_or(Value::Null, |(h, _)| (*h).into()),
                ),
                (
                    "sink-hedge",
                    sink.as_ref().map_or(Value::Null, |(h, _)| (*h).into()),
                ),
                ("particle", particle.name.clone().into()),
                ("pdg", particle.pdg_code.into()),
            ]
            .into_iter()
            .collect();
            if let Some(external) = &edge.external {
                details.insert("external-state", external.state.as_str());
                details.insert("external-index", external.index);
                details.insert("external-name", external.name.clone());
                if ends.source.is_some() && ends.target.is_some() {
                    details.insert("is_cut", external.connection);
                }
            }
            details.insert("momentum", momentum.format_momentum());

            edges.push((
                TypstEdgeSpec {
                    name: None,
                    source: source.map(|(_, end)| end),
                    sink: sink.map(|(_, end)| end),
                    data: None,
                    orientation: Some(orientation.to_owned()),
                    flow: None,
                    id: Some(id.0),
                    pos: None,
                    statements,
                },
                EdgeDrawing {
                    stroke,
                    pattern: drawing.pattern,
                    flow: drawing.flow,
                    momentum: options.momentum_arrows,
                    label,
                    details,
                },
            ));
        }

        // Open sewn initial-state connections into their two original halves.
        // The inspection IDs remain those of the unsplit physics graph.
        if options.split_initial_state {
            let initial_count = edges
                .iter()
                .filter(|(spec, drawing)| {
                    spec.source.is_some()
                        && spec.sink.is_some()
                        && drawing.details.get("external-state").is_some()
                })
                .count();
            let mut initial_index = 0;
            let mut opened = Vec::new();
            for (spec, drawing) in edges {
                let initial = drawing
                    .details
                    .get("external-state")
                    .and_then(Value::as_str)
                    .is_some();
                if initial && spec.source.is_some() && spec.sink.is_some() {
                    let mut left = spec.clone();
                    left.sink = None;
                    let mut right = spec;
                    right.source = None;
                    let y = (initial_index as f64 - (initial_count - 1) as f64 / 2.0) * 4.0;
                    initial_index += 1;
                    for (edge, side) in [(&mut left, -1.0), (&mut right, 1.0)] {
                        edge.statements
                            .insert("pos".into(), format!("{},{y}", side * EXTERNAL_SPACING));
                        edge.statements.insert("pos-mode".into(), "pin".into());
                        edge.statements.insert("pos-z".into(), "0".into());
                        edge.statements.insert("pos-z-mode".into(), "pin".into());
                    }
                    let mut left_drawing = drawing.clone();
                    let mut right_drawing = drawing;
                    let name = left_drawing
                        .details
                        .get("name")
                        .and_then(Value::as_str)
                        .unwrap_or("edge")
                        .to_owned();
                    left_drawing
                        .details
                        .insert("name", format!("{name}-source"));
                    right_drawing.details.insert("name", format!("{name}-sink"));
                    opened.push((left, left_drawing));
                    opened.push((right, right_drawing));
                } else {
                    opened.push((spec, drawing));
                }
            }
            edges = opened;
        }
        for (index, (edge, _)) in edges.iter_mut().enumerate() {
            edge.id = Some(index);
        }
        let (node_specs, mut node_drawings): (Vec<_>, Vec<_>) = nodes.into_iter().unzip();
        if options.show_node_index {
            for (index, node) in node_drawings.iter_mut().enumerate() {
                node.label = Some(pages.len());
                pages.push(format!("[$ upright(v)_{index} $]"));
            }
        }
        let (edge_specs, edge_drawings) = edges.into_iter().unzip();
        Ok(Scene {
            graph: TypstGraphSpec {
                name: Some(self.name().to_owned()),
                data: None,
                statements: BTreeMap::new(),
                default_edge_statements: BTreeMap::new(),
                default_node_statements: BTreeMap::new(),
                nodes: node_specs,
                edges: edge_specs,
            },
            nodes: node_drawings,
            edges: edge_drawings,
            preamble: format!("#set text(size: 9pt, fill: rgb({INK:?}))"),
            // Notebook captions carry the name; callers can request an SVG title.
            title: None,
            pages,
            layout: options.layout.clone(),
            layout_edges: None,
            label_feedback: options.label_feedback,
        })
    }

    /// The faded, dotted stroke of everything outside a highlighted region.
    fn outside_stroke() -> Stroke {
        Stroke {
            paint: OUTSIDE.to_owned(),
            width: 0.6,
            dash: Dash::Dotted,
            round_cap: false,
        }
    }
}

impl SceneOptions {
    fn particle_drawing(
        &self,
        particle: &Particle,
        model: &Model,
        orientation: &str,
    ) -> EdgeDrawing {
        let style = particle.line_style(model);
        let base = Stroke {
            paint: if style.charged { ACCENT } else { INK }.to_owned(),
            width: if style.massive { 1.55 } else { 1.0 },
            dash: match style.kind {
                // `dashed` is (0.1em, 0.45em) at the 9pt drawing text size.
                LineKind::Dashed => Dash::Dashed(0.9, 4.05),
                LineKind::Dotted => Dash::Dotted,
                _ => Dash::Solid,
            },
            round_cap: true,
        };
        let pattern = match style.kind {
            LineKind::Wave => Some(Pattern::Wave {
                amplitude: 0.14,
                wavelength: 0.55,
            }),
            LineKind::Coil => Some(Pattern::Coil {
                amplitude: 0.15,
                wavelength: 0.45,
                longitudinal_scale: 1.4,
            }),
            LineKind::Zigzag => Some(Pattern::Zigzag {
                amplitude: 0.14,
                wavelength: 0.55,
            }),
            _ => None,
        };

        EdgeDrawing {
            stroke: base,
            pattern,
            flow: if style.fermion_flow {
                match orientation {
                    "default" => Some(true),
                    "reversed" => Some(false),
                    _ => None,
                }
            } else {
                None
            },
            momentum: self.momentum_arrows,
            label: None,
            details: Details::default(),
        }
    }

    /// A process is a star graph with a finite interaction region.
    pub fn process_scene(
        &self,
        model: &Model,
        incoming: &[ParticleId],
        outgoing: &[ParticleId],
        config: &Config,
    ) -> Result<Scene, String> {
        let node = TypstNodeSpec {
            name: Some("process".into()),
            index: Some(0),
            data: None,
            pos: None,
            statements: BTreeMap::from([
                ("pos".into(), "0,0".into()),
                ("pin".into(), "x:0,y:0".into()),
            ]),
        };
        let mut scene = Scene {
            graph: TypstGraphSpec {
                name: None,
                data: None,
                statements: BTreeMap::new(),
                default_edge_statements: BTreeMap::new(),
                default_node_statements: BTreeMap::new(),
                nodes: vec![node],
                edges: vec![],
            },
            nodes: vec![NodeDrawing {
                // Preserve the original 3 × 1.5em blob at 10pt in the native
                // renderer's 13.5pt drawing units.
                radius: 10.0 / 3.0,
                label: None,
                rectangular: false,
                fill: "url(#linnest-process-hatch)".into(),
                stroke: Stroke {
                    paint: INK.into(),
                    width: 0.7,
                    dash: Dash::Dashed(3.0, 3.0),
                    round_cap: false,
                },
                details: [("title", "Process")].into_iter().collect(),
            }],
            edges: vec![],
            preamble: format!("#set text(size: 10pt, fill: rgb({INK:?}))"),
            title: None,
            pages: vec![],
            layout: self.layout.clone(),
            label_feedback: self.label_feedback,
            layout_edges: None,
        };
        for (is_incoming, particles) in [(true, incoming), (false, outgoing)] {
            for (rank, id) in particles.iter().enumerate() {
                let index = scene.edges.len();
                let particle = model
                    .particle_by_id(*id)
                    .expect("validated process particle");
                let orientation = if particle.antiparticle == *id {
                    "undirected"
                } else if particle.is_antiparticle() {
                    "reversed"
                } else {
                    "default"
                };
                let end = TypstEndpointSpec {
                    node: 0,
                    statement: None,
                    id: Some(index),
                    data: None,
                    port_label: None,
                    compass: None,
                    in_subgraph: false,
                    route_points: vec![],
                };
                let (source, sink) = if is_incoming {
                    (None, Some(end))
                } else {
                    (Some(end), None)
                };
                scene.graph.edges.push(TypstEdgeSpec {
                    name: None,
                    source,
                    sink,
                    data: None,
                    orientation: Some(orientation.into()),
                    flow: None,
                    id: Some(index),
                    pos: None,
                    statements: BTreeMap::new(),
                });
                let mut drawing = self.particle_drawing(particle, model, orientation);
                if self.show_particle {
                    drawing.label = Some(scene.pages.len());
                    scene
                        .pages
                        .push(format!("[$ {} $]", particle.typst_label()));
                }
                drawing.details = [
                    ("particle", Value::from(particle.name.clone())),
                    ("pdg", particle.pdg_code.into()),
                    (
                        "external-state",
                        if is_incoming { "incoming" } else { "outgoing" }.into(),
                    ),
                    ("external-index", rank.into()),
                ]
                .into_iter()
                .collect();
                scene.edges.push(drawing);
            }
        }
        config.apply(&mut scene)?;
        // A process has a finite interaction region, not a point vertex. Place
        // straight legs around its final styled radius, leaving visible line
        // outside even when the caller enlarges the blob. Mixed states occupy
        // opposite semicircles in input order; one-sided states use the full ring.
        let reach = scene.nodes[0].radius.max(1.0) * 1.65;
        let mixed = !incoming.is_empty() && !outgoing.is_empty();
        for (index, edge) in scene.graph.edges.iter_mut().enumerate() {
            let is_incoming = index < incoming.len();
            let (rank, count) = if is_incoming {
                (index, incoming.len())
            } else {
                (index - incoming.len(), outgoing.len())
            };
            let angle = if mixed {
                std::f64::consts::PI * (rank + 1) as f64 / (count + 1) as f64
            } else {
                std::f64::consts::TAU * rank as f64 / count as f64
            };
            let x = reach * angle.sin() * if is_incoming { -1.0 } else { 1.0 };
            let y = reach * angle.cos();
            edge.statements = [
                ("pos", format!("{x},{y}")),
                ("pin", format!("x:{x},y:{y}")),
                ("pos-z", "0".into()),
                ("pos-z-mode", "pin".into()),
            ]
            .into_iter()
            .map(|(k, v)| (k.into(), v))
            .collect();
        }
        Ok(scene)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use linnest::svg::Dash;
    use linnet::half_edge::{
        involution::Hedge,
        subgraph::{ModifySubSet, SuBitGraph},
    };

    use super::{HIGHLIGHT, INK, OUTSIDE};
    use crate::{SceneOptions, display::tests::one_loop};

    #[test]
    fn describes_amplitudes_as_the_typst_layout_prepares_them() {
        let diagram = one_loop();
        let scene = diagram
            .to_scene(None, &BTreeSet::new(), None, &SceneOptions::default())
            .unwrap();
        // Incoming legs start in the left group and outgoing legs in the right
        // one, every external endpoint at raw depth zero.
        assert_eq!(scene.graph.nodes.len(), 2);
        let statements = |edge: usize| &scene.graph.edges[edge].statements;
        assert_eq!(statements(0)["pin"], "x:@-left");
        assert_eq!(statements(0)["pos"], "-10,0");
        assert_eq!(statements(3)["pin"], "x:@+right");
        assert_eq!(statements(3)["pos-z-mode"], "pin");
        assert!(statements(1).is_empty());
        // Massive scalars are dashed 1.55pt ink lines without arrowheads.
        let scalar = &scene.edges[1];
        assert_eq!(scalar.stroke.paint, INK);
        assert_eq!(scalar.stroke.width, 1.55);
        assert_eq!(scalar.stroke.dash, Dash::Dashed(0.9, 4.05));
        assert_eq!(
            (scalar.pattern, scalar.flow, scalar.momentum),
            (None, None, false)
        );
        // Every edge shares one particle label; the name stays outside the SVG.
        assert_eq!(scene.pages, ["[$ phi $]"]);
        assert!(scene.edges.iter().all(|edge| edge.label == Some(0)));
        assert!(scene.title.is_none());
        // Half-edges are numbered in builder order, as Typst's `build` does.
        assert_eq!(scalar.details.get("source-hedge"), Some(&1.into()));
        assert_eq!(scalar.details.get("sink-hedge"), Some(&2.into()));
        let document = scene.label_document();
        assert!(!document.contains("#import"));
        assert!(document.contains("[$ phi $]"));
    }

    #[test]
    fn labels_momenta_and_highlights_regions() {
        let diagram = one_loop();
        let options = SceneOptions {
            momentum_arrows: true,
            ..SceneOptions::default()
        };
        let mut selected = diagram.graph.empty_subgraph::<SuBitGraph>();
        selected.add(Hedge(1));
        let scene = diagram
            .to_scene(Some(&selected), &BTreeSet::new(), None, &options)
            .unwrap();
        // Momenta follow the arrows into the labels; both legs carry the same one.
        assert!(scene.edges.iter().all(|edge| edge.momentum));
        assert_eq!(scene.pages.len(), 3);
        assert_eq!(scene.edges[0].label, scene.edges[3].label);
        let internal = &scene.pages[scene.edges[1].label.unwrap()];
        assert!(internal.contains("k_0"));
        // A selected half highlights its edge and vertex; the rest fades.
        assert_eq!(scene.edges[1].stroke.paint, HIGHLIGHT);
        assert_eq!(scene.edges[1].stroke.dash, Dash::Dashed(0.9, 4.05));
        assert_eq!(scene.edges[2].stroke.paint, OUTSIDE);
        assert_eq!(scene.edges[2].stroke.dash, Dash::Dotted);
        assert_eq!(scene.nodes[0].stroke.paint, HIGHLIGHT);
        assert_eq!(scene.nodes[1].stroke.paint, OUTSIDE);
    }
}
