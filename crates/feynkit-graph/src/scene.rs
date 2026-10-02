//! Native drawing of a diagram: the physics scene `to_linnest` describes in
//! Typst, laid out and painted by `linnest::svg`. Typst typesets only labels.
use std::collections::{BTreeMap, BTreeSet};

use feynkit_model::LineKind;
use linnest::{
    TypstEdgeSpec, TypstEndpointSpec, TypstGraphSpec, TypstNodeSpec,
    svg::{Dash, Details, EdgeDrawing, NodeDrawing, Pattern, Scene, Stroke},
};
use linnet::half_edge::{
    EdgeAccessors,
    involution::{EdgeIndex, Orientation},
    subgraph::{SuBitGraph, SubSetLike},
};
use serde_json::{Map, Value};

use crate::{DiagramError, FeynmanDiagram, LoopMomentumBasis, display::typst_string};

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
            label_feedback: false,
            layout: Map::new(),
        }
    }
}

impl FeynmanDiagram {
    /// The native drawing of this diagram, as `to_linnest` renders it with
    /// Typst: particle lines and flow arrows, momentum arrows, highlighting and
    /// inspection details. Cross sections open their initial states before
    /// layout, which only the Typst renderer does; they return `None`.
    pub fn to_scene(
        &self,
        highlight: Option<&SuBitGraph>,
        isolated: &BTreeSet<usize>,
        lmb: Option<&LoopMomentumBasis>,
        options: &SceneOptions,
    ) -> Result<Option<Scene>, DiagramError> {
        if !self.cuts().is_empty() {
            return Ok(None);
        }
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
                    NodeDrawing { stroke, details },
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

            let style = particle.line_style(self.model());
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
            let stroke = match highlight {
                None => base,
                Some(_) if selected_edges.contains(&id.0) => Stroke {
                    paint: HIGHLIGHT.to_owned(),
                    ..base
                },
                Some(_) => Self::outside_stroke(),
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

            let momentum = &basis.edge_signatures[&id];
            let label =
                (options.show_particle || show_momentum || options.show_edge_index).then(|| {
                    let mut record = String::from("(");
                    if options.show_edge_index {
                        record.push_str(&format!("eid: {}, ", id.0));
                    }
                    record.push_str(&format!("particle: {}, ", typst_string(&particle.name)));
                    if show_momentum {
                        let (loops, external) = momentum.integer_coefficients();
                        let list = |values: Vec<isize>| {
                            values
                                .iter()
                                .map(|value| format!("{value},"))
                                .collect::<String>()
                        };
                        record.push_str(&format!(
                            "momentum-signature: (loops: ({}), external: ({})), momentum: {}, ",
                            list(loops),
                            list(external),
                            typst_string(&momentum.format_momentum())
                        ));
                    }
                    let content = format!("label({record}))");
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
                    pattern,
                    flow: match orientation {
                        _ if !style.fermion_flow => None,
                        "default" => Some(true),
                        "reversed" => Some(false),
                        _ => None,
                    },
                    momentum: options.momentum_arrows,
                    label,
                    details,
                },
            ));
        }

        let (node_specs, node_drawings) = nodes.into_iter().unzip();
        let (edge_specs, edge_drawings) = edges.into_iter().unzip();
        let preamble = format!(
            "#import \"assets/embedded/drawing/templates/physics-edge-style.typ\" as physics\n\
             #import physics: mi, palette, massive, massless, dashed, dotted, source-stroke, sink-stroke, fermion-flow, wave, coil, zigzag\n\
             #set text(size: 9pt, fill: palette.ink)\n\
             {}\
             #let label(edge) = physics.edge-label(edge, map: particle-map, show-momentum: {show_momentum}, show-particle: {}, show-edge-index: {}, label-fill: palette.ink)\n",
            self.typst_particle_map(),
            if options.show_particle {
                "auto"
            } else {
                "false"
            },
            options.show_edge_index,
        );
        Ok(Some(Scene {
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
            preamble,
            // The title row is always shown; an empty name still takes a line.
            title: Some(if self.name().is_empty() {
                "#hide[X]".to_owned()
            } else {
                format!("#{}", typst_string(self.name()))
            }),
            pages,
            layout: options.layout.clone(),
            label_feedback: options.label_feedback,
        }))
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
            .unwrap()
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
        // Every edge shares one particle label below the title row.
        assert_eq!(scene.pages, ["label((particle: \"phi\", ))"]);
        assert!(scene.edges.iter().all(|edge| edge.label == Some(0)));
        assert_eq!(scene.title.as_deref(), Some("#\"bubble\""));
        // Half-edges are numbered in builder order, as Typst's `build` does.
        assert_eq!(scalar.details.get("source-hedge"), Some(&1.into()));
        assert_eq!(scalar.details.get("sink-hedge"), Some(&2.into()));
        let document = scene.label_document();
        assert!(document.contains("#let particle-map = (\n  \"phi\": "));
        assert!(document.contains("show-momentum: false, show-particle: auto"));
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
            .unwrap()
            .unwrap();
        // Momenta follow the arrows into the labels; both legs carry the same one.
        assert!(scene.edges.iter().all(|edge| edge.momentum));
        assert_eq!(scene.pages.len(), 3);
        assert_eq!(scene.edges[0].label, scene.edges[3].label);
        let internal = &scene.pages[scene.edges[1].label.unwrap()];
        assert!(internal.contains("momentum-signature: (loops: (1,), external: ("));
        // A selected half highlights its edge and vertex; the rest fades.
        assert_eq!(scene.edges[1].stroke.paint, HIGHLIGHT);
        assert_eq!(scene.edges[1].stroke.dash, Dash::Dashed(0.9, 4.05));
        assert_eq!(scene.edges[2].stroke.paint, OUTSIDE);
        assert_eq!(scene.edges[2].stroke.dash, Dash::Dotted);
        assert_eq!(scene.nodes[0].stroke.paint, HIGHLIGHT);
        assert_eq!(scene.nodes[1].stroke.paint, OUTSIDE);
    }
}
