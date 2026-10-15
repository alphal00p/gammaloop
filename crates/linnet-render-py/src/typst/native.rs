//! Convert typed graph drawings to native SVG scenes.
use super::{NativeValue, RenderSettingsTransport};
use crate::drawing::{DrawingKind, EDGE_FIELDS, NODE_FIELDS};
use linnest::svg::{Config, Dash, EdgeDrawing, NodeDrawing, Scene, Stroke};
use linnest::TypstGraphSpec;
use std::collections::BTreeMap;

impl NativeValue {
    fn fields(&self) -> Option<BTreeMap<String, Self>> {
        match self {
            Self::None | Self::Inherit => Some(BTreeMap::new()),
            Self::Dict(values) | Self::Stroke(values) => Some(
                values
                    .iter()
                    .filter(|(_, v)| !matches!(v, Self::Inherit))
                    .map(|(k, v)| (k.clone(), v.clone()))
                    .collect(),
            ),
            _ => None,
        }
    }

    fn number(&self) -> Option<f64> {
        self.json().ok()?.as_f64()
    }

    fn label(&self, pages: &mut Vec<String>, default: String) -> Option<Option<usize>> {
        if matches!(self, Self::None) {
            return Some(None);
        }
        let source = match self {
            Self::Auto | Self::Inherit => default,
            Self::String(_) | Self::Text(_) | Self::Math(_) => {
                self.source(&BTreeMap::new()).ok()?
            }
            _ => return None,
        };
        let index = pages.len();
        pages.push(source);
        Some(Some(index))
    }
}

fn ink(width: f64) -> Stroke {
    Stroke {
        paint: "#000000".into(),
        width,
        dash: Dash::Solid,
        round_cap: false,
    }
}

// Reuse the domain renderer's interpretation of colors, lengths and strokes.
fn apply_style(
    scene: &mut Scene,
    node: bool,
    values: &BTreeMap<String, NativeValue>,
) -> Option<()> {
    let group = if node { "node-style" } else { "edge-style" };
    let mut config = Config::default();
    config
        .style
        .insert(group.into(), NativeValue::Dict(values.clone()).json().ok()?);
    config.apply(scene).ok()
}

impl RenderSettingsTransport {
    pub fn native_scene(
        &self,
        graph: TypstGraphSpec,
        selection: Option<&(Vec<bool>, Vec<usize>)>,
        defaults: Option<Scene>,
    ) -> pyo3::PyResult<Scene> {
        self.scene(graph, selection, defaults).ok_or_else(|| pyo3::exceptions::PyValueError::new_err(
            "drawing settings are not supported by the native SVG renderer; Typst templates, imports, and executable styles are not supported"
        ))
    }

    fn scene(
        &self,
        graph: TypstGraphSpec,
        selection: Option<&(Vec<bool>, Vec<usize>)>,
        defaults: Option<Scene>,
    ) -> Option<Scene> {
        if !self.imports.is_empty() {
            return None;
        }
        let mut config = self.native.fields()?;
        config.remove("schema");
        config.remove("version");
        let elements = config.remove("elements")?.fields()?;
        if !elements.get("graph")?.fields()?.is_empty() {
            return None;
        }
        let mut draw = config
            .remove("draw")
            .unwrap_or(NativeValue::None)
            .fields()?;
        let mut style = config
            .remove("style")
            .unwrap_or(NativeValue::None)
            .fields()?;
        if !config
            .remove("options")
            .unwrap_or(NativeValue::None)
            .fields()?
            .is_empty()
        {
            return None;
        }
        let title = draw.remove("title").or_else(|| config.remove("title"));
        config.remove("title");
        let layouts = config.remove("layouts");
        if !config.is_empty() {
            return None;
        }
        let mut scene = defaults.unwrap_or_else(|| Scene {
            graph,
            nodes: Vec::new(),
            edges: Vec::new(),
            preamble: "#set text(size: 9pt, fill: black)".into(),
            title: None,
            pages: Vec::new(),
            layout: Default::default(),
            layout_passes: Vec::new(),
            layout_edges: None,
            label_feedback: false,
        });
        let default_nodes = std::mem::take(&mut scene.nodes);
        let default_edges = std::mem::take(&mut scene.edges);
        if let Some(title) = title {
            if !matches!(
                title,
                NativeValue::None | NativeValue::Auto | NativeValue::Inherit
            ) {
                let index = title.label(&mut scene.pages, String::new())??;
                scene.title = Some(format!("#{}", scene.pages.remove(index)));
            }
        }
        if let Some(layouts) = layouts {
            let passes = match layouts {
                NativeValue::None | NativeValue::Inherit => Vec::new(),
                NativeValue::Array(passes) => passes
                    .iter()
                    .map(|pass| Some(pass.json().ok()?.as_object()?.clone()))
                    .collect::<Option<Vec<_>>>()?,
                _ => return None,
            };
            Config {
                layouts: passes,
                ..Default::default()
            }
            .apply(&mut scene)
            .ok()?;
        }
        // The default drawing unit is shared with native domain scenes. Other
        // units and executable style/layout functions are rejected here.
        if let Some(unit) = draw.remove("unit") {
            if !matches!(unit, NativeValue::Auto) {
                return None;
            }
        }
        let node_label = draw
            .remove("node-label")
            .or_else(|| style.remove("node-label"))
            .unwrap_or(NativeValue::Auto);
        let edge_label = draw
            .remove("edge-label")
            .or_else(|| style.remove("edge-label"))
            .unwrap_or(NativeValue::Auto);
        let mut node_style = style
            .remove("node-style")
            .unwrap_or(NativeValue::None)
            .fields()?;
        node_style.extend(
            draw.remove("node-style")
                .unwrap_or(NativeValue::None)
                .fields()?,
        );
        let mut edge_style = style
            .remove("edge-style")
            .unwrap_or(NativeValue::None)
            .fields()?;
        edge_style.extend(
            draw.remove("edge-style")
                .unwrap_or(NativeValue::None)
                .fields()?,
        );
        if !style.is_empty() {
            return None;
        }
        let radius = draw.remove("node-radius").unwrap_or(NativeValue::Auto);
        let minimum = draw
            .remove("node-min-radius")
            .unwrap_or(NativeValue::Float(0.16))
            .number()?;
        let padding = draw
            .remove("node-label-padding")
            .unwrap_or(NativeValue::Float(0.08))
            .number()?;
        let automatic = matches!(radius, NativeValue::Auto);
        let radius = if automatic { minimum } else { radius.number()? };
        let override_fill = draw.contains_key("node-fill");
        let fill = draw
            .remove("node-fill")
            .unwrap_or(NativeValue::String("#ffffff".into()));
        let node_stroke = draw.remove("node-stroke");
        let edge_stroke = draw.remove("edge-stroke");
        if !draw.is_empty() {
            return None;
        }
        let NativeValue::Array(nodes) = elements.get("nodes")? else {
            return None;
        };
        let NativeValue::Array(edges) = elements.get("edges")? else {
            return None;
        };
        let NativeValue::Array(hedges) = elements.get("hedges")? else {
            return None;
        };
        if hedges
            .iter()
            .any(|h| h.fields().is_none_or(|f| !f.is_empty()))
        {
            return None;
        }
        for (index, record) in nodes.iter().enumerate() {
            let mut record = record.fields()?;
            record.remove("inspection");
            let label = record.remove("label").unwrap_or_else(|| node_label.clone());
            let mut local_style = node_style.clone();
            local_style.extend(
                record
                    .remove("node-style")
                    .unwrap_or(NativeValue::None)
                    .fields()?,
            );
            if let Some(position) = record.remove("pos") {
                scene.graph.nodes[index].pos =
                    serde_json::from_value(position.json().ok()?).ok()?;
            }
            if record.keys().any(|key| {
                NODE_FIELDS
                    .iter()
                    .any(|field| DrawingKind::Node.transport_key(field) == key)
            }) {
                return None;
            }
            let override_radius = !automatic || local_style.contains_key("radius");
            let local_radius = local_style.remove("radius");
            let (radius, label_padding) = match local_radius {
                Some(NativeValue::Auto) => (minimum, Some(padding)),
                Some(value) => (value.number()?, None),
                None => (radius, automatic.then_some(padding)),
            };
            let label = if matches!(label, NativeValue::Auto | NativeValue::Inherit)
                && !default_nodes.is_empty()
            {
                default_nodes.get(index)?.label
            } else {
                label.label(&mut scene.pages, format!("$n_({index})$"))?
            };
            let mut node = default_nodes.get(index).cloned().unwrap_or(NodeDrawing {
                radius,
                label_padding,
                label,
                rectangular: false,
                fill: fill.json().ok()?.as_str()?.into(),
                stroke: ink(0.5),
                details: scene.graph.nodes[index]
                    .name
                    .iter()
                    .map(|name| ("name", name.clone()))
                    .collect(),
            });
            node.label = label;
            if override_fill {
                node.fill = fill.json().ok()?.as_str()?.into();
            }
            if override_radius {
                node.radius = radius;
                node.label_padding = label_padding;
            }
            for (key, value) in record {
                node.details.insert(&key, value.json().ok()?);
            }
            scene.nodes.push(node);
            // Config applies to every node, so isolate this record while merging
            // its explicit style; earlier records must remain untouched.
            let previous = std::mem::take(&mut scene.nodes);
            scene.nodes = vec![previous.last()?.clone()];
            if let Some(stroke) = &node_stroke {
                apply_style(
                    &mut scene,
                    true,
                    &BTreeMap::from([("stroke".into(), stroke.clone())]),
                )?;
            }
            apply_style(&mut scene, true, &local_style)?;
            let node = scene.nodes.pop()?;
            scene.nodes = previous;
            *scene.nodes.last_mut()? = node;
        }
        for (index, record) in edges.iter().enumerate() {
            let mut record = record.fields()?;
            record.remove("inspection");
            let label = record.remove("label").unwrap_or_else(|| edge_label.clone());
            let mut local_style = edge_style.clone();
            local_style.extend(
                record
                    .remove("edge-style")
                    .unwrap_or(NativeValue::None)
                    .fields()?,
            );
            if let Some(position) = record.remove("pos") {
                scene.graph.edges[index].pos =
                    serde_json::from_value(position.json().ok()?).ok()?;
            }
            if record.keys().any(|key| {
                EDGE_FIELDS
                    .iter()
                    .any(|field| DrawingKind::Edge.transport_key(field) == key)
            }) {
                return None;
            }
            let label = if matches!(label, NativeValue::Auto | NativeValue::Inherit)
                && !default_edges.is_empty()
            {
                default_edges.get(index)?.label
            } else {
                label.label(&mut scene.pages, format!("$e_({index})$"))?
            };
            let previous = std::mem::take(&mut scene.edges);
            let mut edge = default_edges.get(index).cloned().unwrap_or(EdgeDrawing {
                stroke: ink(0.9),
                pattern: None,
                flow: None,
                momentum: false,
                momentum_offset: EdgeDrawing::DEFAULT_MOMENTUM_OFFSET,
                label,
                details: scene.graph.edges[index]
                    .name
                    .iter()
                    .map(|name| ("name", name.clone()))
                    .collect(),
            });
            edge.label = label;
            for (key, value) in record {
                edge.details.insert(&key, value.json().ok()?);
            }
            scene.edges.push(edge);
            if let Some(stroke) = &edge_stroke {
                apply_style(
                    &mut scene,
                    false,
                    &BTreeMap::from([("stroke".into(), stroke.clone())]),
                )?;
            }
            apply_style(&mut scene, false, &local_style)?;
            let edge = scene.edges.pop()?;
            scene.edges = previous;
            scene.edges.push(edge);
        }
        if let Some((hedges, nodes)) = selection {
            let selected = Stroke {
                paint: "#c58b13".into(),
                width: 1.2,
                ..ink(1.2)
            };
            let outside = Stroke {
                paint: "#77777773".into(),
                width: 0.6,
                dash: Dash::Dotted,
                round_cap: false,
            };
            for (index, drawing) in scene.nodes.iter_mut().enumerate() {
                let inside = nodes.contains(&index);
                drawing.stroke = if inside {
                    selected.clone()
                } else {
                    outside.clone()
                };
                drawing.fill = if inside { "#ffd16659" } else { "#77777726" }.into();
            }
            for (edge, drawing) in scene.graph.edges.iter_mut().zip(&mut scene.edges) {
                let mut included = None;
                for endpoint in [&mut edge.source, &mut edge.sink].into_iter().flatten() {
                    let inside = *hedges.get(endpoint.id?)?;
                    // Different half-edge paints require the full drawing API.
                    if included.is_some_and(|previous| previous != inside) {
                        return None;
                    }
                    included = Some(inside);
                    endpoint.in_subgraph = inside;
                }
                drawing.stroke = if included? {
                    selected.clone()
                } else {
                    outside.clone()
                };
            }
        }
        Some(scene)
    }
}
