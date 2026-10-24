//! Convert typed graph drawings to native SVG scenes.
use super::{
    LengthUnit, MathScript, MathValue, NativeEnum, NativeValue, RenderSettingsTransport, TextValue,
};
use crate::drawing::{DrawingKind, EDGE_FIELDS, NODE_FIELDS};
use linnest::svg::{Config, Dash, EdgeDrawing, LabelSource, NodeDrawing, Scene, Stroke};
use linnest::TypstGraphSpec;
use pyo3::exceptions::PyValueError;
use pyo3::{PyErr, PyResult, Python};
use serde_json::{Error as JsonError, Value};
use symbolica_notation::IndexDisplay;

use std::collections::BTreeMap;
use std::mem;
use std::num::TryFromIntError;
use symbolica::printer::AnsiHtmlFormatter;
use thiserror::Error;

/// Private conversion failures, retained until the Python presentation boundary.
#[derive(Debug, Error)]
enum NativeSceneError {
    #[error("drawing settings are not supported by the native SVG renderer: {0}")]
    Unsupported(&'static str),
    #[error("invalid native drawing record: {0}")]
    InvalidShape(&'static str),
    #[error("invalid label style property {property:?}")]
    InvalidStyle { property: String },
    #[error("label size must be finite and nonnegative: {0}")]
    InvalidSize(f64),
    #[error("expected a native number")]
    InvalidNumber,
    #[error("label has no supported MathML representation")]
    UnsupportedMathml,
    #[error("cannot encode native drawing value: {0}")]
    JsonValue(#[source] PyErr),
    #[error("cannot encode Typst label source: {0}")]
    TypstSource(#[source] PyErr),
    #[error("invalid drawing position: {0}")]
    Position(#[from] JsonError),
    #[error("label index exceeds the supported integer range: {0}")]
    LabelIndex(#[from] TryFromIntError),
    // Config's existing API originates strings, so there is no structured cause.
    #[error("invalid native drawing configuration: {0}")]
    Configuration(String),
    #[error("internal bug while preparing a drawing: {0}")]
    InternalInvariant(&'static str),
}

impl NativeSceneError {
    fn into_py_err(self) -> PyErr {
        let error = PyValueError::new_err(self.to_string());
        if let Self::JsonValue(cause) | Self::TypstSource(cause) = self {
            Python::attach(|py| error.set_cause(py, Some(cause)));
        }
        error
    }
}

impl MathScript {
    fn display(&self) -> IndexDisplay {
        match self {
            Self::Symbol(symbol) => IndexDisplay::Symbol(symbol.clone()),
            Self::Integer(value) => IndexDisplay::Number(*value),
        }
    }
}

impl MathValue {
    fn display(&self) -> IndexDisplay {
        let mut display = IndexDisplay::Symbol(self.symbol.clone());
        if let Some(script) = &self.subscript {
            display = display.with_bottom(script.display());
        }
        if let Some(script) = &self.superscript {
            display = display.with_top(script.display());
        }
        display
    }
}

impl TextValue {
    fn style_label(
        mut source: LabelSource,
        styles: &[&BTreeMap<String, NativeValue>],
    ) -> Result<LabelSource, NativeSceneError> {
        let mut css = BTreeMap::new();
        let mut arguments = Vec::new();
        let mut em = 1.0;
        // Resolve outer-to-inner typed styles before emitting MathML. Absolute
        // sizes replace the inherited size; em sizes multiply it. A single
        // wrapper keeps absolute sizes from being scaled by an outer wrapper.
        for options in styles {
            let mut group = Vec::new();
            for (key, value) in options
                .iter()
                .filter(|(_, value)| !matches!(value, NativeValue::Inherit))
            {
                let property = match key.as_str() {
                    "fill" => {
                        let color = value.json().map_err(NativeSceneError::JsonValue)?;
                        let color =
                            color
                                .as_str()
                                .ok_or_else(|| NativeSceneError::InvalidStyle {
                                    property: key.clone(),
                                })?;
                        format!("color:{color}")
                    }
                    "size" => {
                        // The probe and SVG label both use the fixed 9-unit font size.
                        let size = match value {
                            NativeValue::Length(value, LengthUnit::Em) => *value,
                            _ => value.number()? / 9.0,
                        };
                        if !size.is_finite() || size < 0.0 {
                            return Err(NativeSceneError::InvalidSize(size));
                        }
                        em = if matches!(value, NativeValue::Length(_, LengthUnit::Em)) {
                            em * size
                        } else {
                            size
                        };
                        if !em.is_finite() {
                            return Err(NativeSceneError::InvalidSize(em));
                        }
                        format!("font-size:{em}em")
                    }
                    "weight" => {
                        let weight = match value {
                            NativeValue::Int(value) if (1..=1000).contains(value) => *value,
                            NativeValue::String(value) => match value.as_str() {
                                "thin" => 100,
                                "extralight" => 200,
                                "light" => 300,
                                "regular" | "normal" => 400,
                                "medium" => 500,
                                "semibold" => 600,
                                "bold" => 700,
                                "extrabold" => 800,
                                "black" | "heavy" => 900,
                                _ => {
                                    return Err(NativeSceneError::InvalidStyle {
                                        property: key.clone(),
                                    })
                                }
                            },
                            _ => {
                                return Err(NativeSceneError::InvalidStyle {
                                    property: key.clone(),
                                })
                            }
                        };
                        format!("font-weight:{weight}")
                    }
                    "style" => {
                        let NativeValue::Enum(NativeEnum::TextStyle(style)) = value else {
                            return Err(NativeSceneError::InvalidStyle {
                                property: key.clone(),
                            });
                        };
                        format!("font-style:{}", style.typst_name())
                    }
                    _ => {
                        return Err(NativeSceneError::InvalidStyle {
                            property: key.clone(),
                        })
                    }
                };
                css.insert(key, property);
                let typst = value
                    .source(&BTreeMap::new())
                    .map_err(NativeSceneError::TypstSource)?;
                group.push(format!("{key}: {typst}"));
            }
            if !group.is_empty() {
                arguments.push(group);
            }
        }
        if !css.is_empty() {
            source.mathml = format!(
                "<mstyle style=\"{}\">{}</mstyle>",
                AnsiHtmlFormatter::escape_html(&css.into_values().collect::<Vec<_>>().join(";")),
                source.mathml
            );
            for group in arguments.into_iter().rev() {
                if let Some(typst) = &mut source.typst {
                    *typst = format!("text({typst}, {})", group.join(", "));
                }
            }
        }
        Ok(source)
    }
}

impl NativeValue {
    fn fields(&self) -> Result<BTreeMap<String, Self>, NativeSceneError> {
        match self {
            Self::None | Self::Inherit => Ok(BTreeMap::new()),
            // Scene preparation consumes and overrides its own settings map;
            // the reusable transport retains the caller's original values.
            Self::Dict(values) | Self::Stroke(values) => Ok(values
                .iter()
                .filter(|(_, v)| !matches!(v, Self::Inherit))
                .map(|(k, v)| (k.clone(), v.clone()))
                .collect()),
            _ => Err(NativeSceneError::InvalidShape("expected an option group")),
        }
    }

    fn number(&self) -> Result<f64, NativeSceneError> {
        self.json()
            .map_err(NativeSceneError::JsonValue)?
            .as_f64()
            .ok_or(NativeSceneError::InvalidNumber)
    }

    fn label_source(
        &self,
        default: LabelSource,
        style: &BTreeMap<String, NativeValue>,
    ) -> Result<LabelSource, NativeSceneError> {
        let source = match self {
            Self::Auto | Self::Inherit => default,
            Self::String(text) => LabelSource::text(text),
            Self::Text(value) => {
                return TextValue::style_label(
                    LabelSource::text(&value.text),
                    &[style, &value.options],
                );
            }
            Self::Math(value) => LabelSource {
                mathml: value
                    .display()
                    .to_mathml()
                    .ok_or(NativeSceneError::UnsupportedMathml)?,
                typst: Some(
                    self.source(&BTreeMap::new())
                        .map_err(NativeSceneError::TypstSource)?,
                ),
            },
            _ => {
                return Err(NativeSceneError::Unsupported(
                    "expected a text or math label",
                ))
            }
        };
        TextValue::style_label(source, &[style])
    }

    fn label(
        &self,
        pages: &mut Vec<LabelSource>,
        default: LabelSource,
        style: &BTreeMap<String, NativeValue>,
    ) -> Result<Option<usize>, NativeSceneError> {
        if matches!(self, Self::None) {
            return Ok(None);
        }
        let source = self.label_source(default, style)?;
        let index = pages.len();
        pages.push(source);
        Ok(Some(index))
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

impl RenderSettingsTransport {
    // Reuse the domain renderer's interpretation of colors, lengths and strokes.
    fn apply_style(
        scene: &mut Scene,
        group: &str,
        values: &BTreeMap<String, NativeValue>,
    ) -> Result<(), NativeSceneError> {
        let mut config = Config::default();
        let value = NativeValue::Dict(values.clone())
            .json()
            .map_err(NativeSceneError::JsonValue)?;
        config.style.insert(group.into(), value);
        config.apply(scene).map_err(NativeSceneError::Configuration)
    }

    /// Resolve native drawing options and labels; see [`crate::PyDiagramView`]
    /// for the resulting scene's browser/vector measurement lifecycle.
    ///
    /// # Errors
    /// Returns `ValueError` for invalid values or unsupported executable options,
    /// styles, label sources, positions, and topology records. Errors from Python
    /// source/JSON conversion retain their original exception as the cause.
    pub fn native_scene(
        &self,
        graph: TypstGraphSpec,
        selection: Option<&(Vec<bool>, Vec<usize>)>,
        defaults: Option<Scene>,
    ) -> PyResult<Scene> {
        self.scene(graph, selection, defaults)
            .map_err(NativeSceneError::into_py_err)
    }

    fn scene(
        &self,
        graph: TypstGraphSpec,
        selection: Option<&(Vec<bool>, Vec<usize>)>,
        defaults: Option<Scene>,
    ) -> Result<Scene, NativeSceneError> {
        if !self.imports.is_empty() {
            return Err(NativeSceneError::Unsupported(
                "Typst templates, imports, and executable styles are not supported",
            ));
        }
        let mut config = self.native.fields()?;
        config.remove("schema");
        config.remove("version");
        let elements = config
            .remove("elements")
            .ok_or(NativeSceneError::InvalidShape("missing elements"))?
            .fields()?;
        if !elements
            .get("graph")
            .ok_or(NativeSceneError::InvalidShape("missing graph"))?
            .fields()?
            .is_empty()
        {
            return Err(NativeSceneError::Unsupported(
                "unrecognized native drawing option",
            ));
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
            return Err(NativeSceneError::Unsupported(
                "unrecognized native drawing option",
            ));
        }
        let title = draw.remove("title").or_else(|| config.remove("title"));
        config.remove("title");
        let layouts = config.remove("layouts");
        if !config.is_empty() {
            return Err(NativeSceneError::Unsupported(
                "unrecognized native drawing option",
            ));
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
        let default_nodes = mem::take(&mut scene.nodes);
        let default_edges = mem::take(&mut scene.edges);
        let has_default_nodes = !default_nodes.is_empty();
        let has_default_edges = !default_edges.is_empty();
        let mut default_nodes = default_nodes.into_iter();
        let mut default_edges = default_edges.into_iter();
        if let Some(title) = title {
            if !matches!(
                title,
                NativeValue::None | NativeValue::Auto | NativeValue::Inherit
            ) {
                scene.title = Some(title.label_source(LabelSource::text(""), &BTreeMap::new())?);
            }
        }
        if let Some(layouts) = layouts {
            let passes = match layouts {
                NativeValue::None | NativeValue::Inherit => Vec::new(),
                NativeValue::Array(passes) => passes
                    .iter()
                    .map(
                        |pass| match pass.json().map_err(NativeSceneError::JsonValue)? {
                            Value::Object(options) => Ok(options),
                            _ => Err(NativeSceneError::InvalidShape(
                                "layout pass must be an object",
                            )),
                        },
                    )
                    .collect::<Result<Vec<_>, _>>()?,
                _ => return Err(NativeSceneError::InvalidShape("layouts must be a sequence")),
            };
            Config {
                layouts: passes,
                ..Default::default()
            }
            .apply(&mut scene)
            .map_err(NativeSceneError::Configuration)?;
        }
        // The default drawing unit is shared with native domain scenes. Other
        // units and executable style/layout functions are rejected here.
        if let Some(unit) = draw.remove("unit") {
            if !matches!(unit, NativeValue::Auto) {
                return Err(NativeSceneError::Unsupported(
                    "unrecognized native drawing option",
                ));
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
        let mut node_label_style = style
            .remove("node-label-style")
            .unwrap_or(NativeValue::None)
            .fields()?;
        node_label_style.extend(
            draw.remove("node-label-style")
                .unwrap_or(NativeValue::None)
                .fields()?,
        );
        let mut edge_label_style = style
            .remove("edge-label-style")
            .unwrap_or(NativeValue::None)
            .fields()?;
        edge_label_style.extend(
            draw.remove("edge-label-style")
                .unwrap_or(NativeValue::None)
                .fields()?,
        );
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
            return Err(NativeSceneError::Unsupported(
                "unrecognized native drawing option",
            ));
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
            return Err(NativeSceneError::Unsupported(
                "unrecognized native drawing option",
            ));
        }
        let NativeValue::Array(nodes) = elements
            .get("nodes")
            .ok_or(NativeSceneError::InvalidShape("missing nodes"))?
        else {
            return Err(NativeSceneError::Unsupported(
                "unrecognized native drawing option",
            ));
        };
        let NativeValue::Array(edges) = elements
            .get("edges")
            .ok_or(NativeSceneError::InvalidShape("missing edges"))?
        else {
            return Err(NativeSceneError::Unsupported(
                "unrecognized native drawing option",
            ));
        };
        let NativeValue::Array(hedges) = elements
            .get("hedges")
            .ok_or(NativeSceneError::InvalidShape("missing hedges"))?
        else {
            return Err(NativeSceneError::Unsupported(
                "unrecognized native drawing option",
            ));
        };
        for hedge in hedges {
            if !hedge.fields()?.is_empty() {
                return Err(NativeSceneError::Unsupported("half-edge drawing styles"));
            }
        }
        if nodes.len() != scene.graph.nodes.len() || edges.len() != scene.graph.edges.len() {
            return Err(NativeSceneError::InvalidShape(
                "drawing records must match the graph's nodes and edges",
            ));
        }
        for (index, record) in nodes.iter().enumerate() {
            let default_node = default_nodes.next();
            let mut record = record.fields()?;
            record.remove("inspection");
            let label = record.remove("label").unwrap_or_else(|| node_label.clone());
            let mut label_style = node_label_style.clone();
            label_style.extend(
                record
                    .remove("node-label-style")
                    .unwrap_or(NativeValue::None)
                    .fields()?,
            );
            let mut local_style = node_style.clone();
            local_style.extend(
                record
                    .remove("node-style")
                    .unwrap_or(NativeValue::None)
                    .fields()?,
            );
            if let Some(position) = record.remove("pos") {
                scene.graph.nodes[index].pos =
                    serde_json::from_value(position.json().map_err(NativeSceneError::JsonValue)?)?;
            }
            if record.keys().any(|key| {
                NODE_FIELDS
                    .iter()
                    .any(|field| DrawingKind::Node.transport_key(field) == key)
            }) {
                return Err(NativeSceneError::Unsupported(
                    "unrecognized native drawing option",
                ));
            }
            let override_radius = !automatic || local_style.contains_key("radius");
            let local_radius = local_style.remove("radius");
            let (radius, label_padding) = match local_radius {
                Some(NativeValue::Auto) => (minimum, Some(padding)),
                Some(value) => (value.number()?, None),
                None => (radius, automatic.then_some(padding)),
            };
            let label =
                if matches!(label, NativeValue::Auto | NativeValue::Inherit) && has_default_nodes {
                    match default_node
                        .as_ref()
                        .ok_or(NativeSceneError::InvalidShape(
                            "missing default node drawing",
                        ))?
                        .label
                    {
                        Some(page) if !label_style.is_empty() => {
                            let source = TextValue::style_label(
                                scene
                                    .pages
                                    .get(page)
                                    .ok_or(NativeSceneError::InvalidShape(
                                        "default label page is out of bounds",
                                    ))?
                                    .clone(),
                                &[&label_style],
                            )?;
                            let page = scene.pages.len();
                            scene.pages.push(source);
                            Some(page)
                        }
                        label => label,
                    }
                } else {
                    let default = NativeValue::Math(MathValue {
                        symbol: "n".into(),
                        subscript: Some(MathScript::Integer(index.try_into()?)),
                        superscript: None,
                    });
                    label.label(
                        &mut scene.pages,
                        default.label_source(LabelSource::text(""), &BTreeMap::new())?,
                        &label_style,
                    )?
                };
            let mut node = default_node.unwrap_or(NodeDrawing {
                radius,
                label_padding,
                label,
                rectangular: false,
                fill: fill
                    .json()
                    .map_err(NativeSceneError::JsonValue)?
                    .as_str()
                    .ok_or(NativeSceneError::InvalidShape(
                        "node fill must be a color string",
                    ))?
                    .into(),
                stroke: ink(0.5),
                details: scene.graph.nodes[index]
                    .name
                    .iter()
                    .map(|name| ("name", name.clone()))
                    .collect(),
            });
            node.label = label;
            if override_fill {
                node.fill = fill
                    .json()
                    .map_err(NativeSceneError::JsonValue)?
                    .as_str()
                    .ok_or(NativeSceneError::InvalidShape(
                        "node fill must be a color string",
                    ))?
                    .into();
            }
            if override_radius {
                node.radius = radius;
                node.label_padding = label_padding;
            }
            for (key, value) in record {
                node.details
                    .insert(&key, value.json().map_err(NativeSceneError::JsonValue)?);
            }
            // Config applies to every node, so isolate this record while merging
            // its explicit style; earlier records must remain untouched.
            let previous = mem::take(&mut scene.nodes);
            scene.nodes.push(node);
            if let Some(stroke) = &node_stroke {
                Self::apply_style(
                    &mut scene,
                    "node-style",
                    &BTreeMap::from([("stroke".into(), stroke.clone())]),
                )?;
            }
            Self::apply_style(&mut scene, "node-style", &local_style)?;
            let node = scene
                .nodes
                .pop()
                .ok_or(NativeSceneError::InternalInvariant(
                    "style application preserves the isolated node",
                ))?;
            scene.nodes = previous;
            scene.nodes.push(node);
        }
        for (index, record) in edges.iter().enumerate() {
            let default_edge = default_edges.next();
            let mut record = record.fields()?;
            record.remove("inspection");
            let label = record.remove("label").unwrap_or_else(|| edge_label.clone());
            let mut label_style = edge_label_style.clone();
            label_style.extend(
                record
                    .remove("edge-label-style")
                    .unwrap_or(NativeValue::None)
                    .fields()?,
            );
            let mut local_style = edge_style.clone();
            local_style.extend(
                record
                    .remove("edge-style")
                    .unwrap_or(NativeValue::None)
                    .fields()?,
            );
            if let Some(position) = record.remove("pos") {
                scene.graph.edges[index].pos =
                    serde_json::from_value(position.json().map_err(NativeSceneError::JsonValue)?)?;
            }
            if record.keys().any(|key| {
                EDGE_FIELDS
                    .iter()
                    .any(|field| DrawingKind::Edge.transport_key(field) == key)
            }) {
                return Err(NativeSceneError::Unsupported(
                    "unrecognized native drawing option",
                ));
            }
            let label =
                if matches!(label, NativeValue::Auto | NativeValue::Inherit) && has_default_edges {
                    match default_edge
                        .as_ref()
                        .ok_or(NativeSceneError::InvalidShape(
                            "missing default edge drawing",
                        ))?
                        .label
                    {
                        Some(page) if !label_style.is_empty() => {
                            let source = TextValue::style_label(
                                scene
                                    .pages
                                    .get(page)
                                    .ok_or(NativeSceneError::InvalidShape(
                                        "default label page is out of bounds",
                                    ))?
                                    .clone(),
                                &[&label_style],
                            )?;
                            let page = scene.pages.len();
                            scene.pages.push(source);
                            Some(page)
                        }
                        label => label,
                    }
                } else {
                    let default = NativeValue::Math(MathValue {
                        symbol: "e".into(),
                        subscript: Some(MathScript::Integer(index.try_into()?)),
                        superscript: None,
                    });
                    label.label(
                        &mut scene.pages,
                        default.label_source(LabelSource::text(""), &BTreeMap::new())?,
                        &label_style,
                    )?
                };
            let previous = mem::take(&mut scene.edges);
            let mut edge = default_edge.unwrap_or(EdgeDrawing {
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
                edge.details
                    .insert(&key, value.json().map_err(NativeSceneError::JsonValue)?);
            }
            scene.edges.push(edge);
            if let Some(stroke) = &edge_stroke {
                Self::apply_style(
                    &mut scene,
                    "edge-style",
                    &BTreeMap::from([("stroke".into(), stroke.clone())]),
                )?;
            }
            Self::apply_style(&mut scene, "edge-style", &local_style)?;
            let edge = scene
                .edges
                .pop()
                .ok_or(NativeSceneError::InternalInvariant(
                    "style application preserves the isolated edge",
                ))?;
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
                    let id = endpoint.id.ok_or(NativeSceneError::InvalidShape(
                        "selection endpoint has no half-edge id",
                    ))?;
                    let inside = *hedges.get(id).ok_or(NativeSceneError::InvalidShape(
                        "selection half-edge id is out of bounds",
                    ))?;
                    // Different half-edge paints require the full drawing API.
                    if included.is_some_and(|previous| previous != inside) {
                        return Err(NativeSceneError::Unsupported(
                            "unrecognized native drawing option",
                        ));
                    }
                    included = Some(inside);
                    endpoint.in_subgraph = inside;
                }
                drawing.stroke = if included.ok_or(NativeSceneError::InvalidShape(
                    "selected edge has no endpoints",
                ))? {
                    selected.clone()
                } else {
                    outside.clone()
                };
            }
        }
        Ok(scene)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::typst::{ColorValue, PyTextStyle, SelectorCallbacks, TypstModuleSource};
    use roxmltree::Document;
    use serde_json::json;
    use std::error::Error as _;

    #[test]
    fn label_absence_is_distinct_from_invalid_content_and_styles() {
        let mut pages = Vec::new();
        assert_eq!(
            NativeValue::None
                .label(&mut pages, LabelSource::text("default"), &BTreeMap::new())
                .unwrap(),
            None
        );
        assert!(pages.is_empty());
        assert!(matches!(
            NativeValue::Bool(true).label(
                &mut pages,
                LabelSource::text("default"),
                &BTreeMap::new()
            ),
            Err(NativeSceneError::Unsupported(_))
        ));
        for (property, value) in [
            ("weight", NativeValue::Int(1001)),
            ("style", NativeValue::String("italic".into())),
            ("font", NativeValue::String("arbitrary".into())),
        ] {
            let error = TextValue::style_label(
                LabelSource::text("x"),
                &[&BTreeMap::from([(property.into(), value)])],
            )
            .unwrap_err();
            assert!(matches!(
                error,
                NativeSceneError::InvalidStyle { property: actual } if actual == property
            ));
        }
        assert!(matches!(
            NativeValue::Math(MathValue {
                symbol: "\0".into(),
                subscript: None,
                superscript: None,
            })
            .label_source(LabelSource::text("default"), &BTreeMap::new()),
            Err(NativeSceneError::UnsupportedMathml)
        ));
        assert!(pages.is_empty(), "rejected labels must not append a page");
    }

    #[test]
    fn native_label_json_errors_retain_the_python_cause() {
        Python::initialize();
        Python::attach(|py| {
            // Relative lengths cannot be emitted as a native JSON paint. The
            // originating NativeValue error must survive style propagation.
            let error = TextValue::style_label(
                LabelSource::text("x"),
                &[&BTreeMap::from([(
                    "fill".into(),
                    NativeValue::Length(1.0, LengthUnit::Em),
                )])],
            )
            .unwrap_err();
            assert!(matches!(error, NativeSceneError::JsonValue(_)));
            assert!(error.source().unwrap().is::<PyErr>());
            let error = error.into_py_err();
            assert!(error.is_instance_of::<PyValueError>(py));
            let cause = error
                .cause(py)
                .expect("native conversion cause is retained");
            assert!(cause.is_instance_of::<PyValueError>(py));
            assert!(cause.to_string().contains("this option is not supported"));
        });
    }

    #[test]
    fn invalid_positions_preserve_the_json_parser_error() {
        let graph: TypstGraphSpec = serde_json::from_value(json!({
            "nodes": [{}], "edges": [],
        }))
        .unwrap();
        let transport = RenderSettingsTransport {
            native: NativeValue::Dict(BTreeMap::from([(
                "elements".into(),
                NativeValue::Dict(BTreeMap::from([
                    ("graph".into(), NativeValue::None),
                    (
                        "nodes".into(),
                        NativeValue::Array(vec![NativeValue::Dict(BTreeMap::from([(
                            "pos".into(),
                            NativeValue::String("invalid".into()),
                        )]))]),
                    ),
                    ("edges".into(), NativeValue::Array(Vec::new())),
                    ("hedges".into(), NativeValue::Array(Vec::new())),
                ])),
            )])),
            imports: Default::default(),
            selectors: SelectorCallbacks {
                node: None,
                edge: None,
                source: None,
                sink: None,
            },
        };
        let Err(error) = transport.scene(graph, None, None) else {
            panic!("invalid position must fail scene preparation");
        };
        assert!(matches!(error, NativeSceneError::Position(_)));
        assert!(error.source().unwrap().is::<JsonError>());
    }

    #[test]
    fn literal_labels_escape_mathml_and_keep_literal_typst() {
        let text = "<&>\"' #import \"evil.typ\"\n\\";
        let source = NativeValue::String(text.into())
            .label_source(LabelSource::text("unused"), &BTreeMap::new())
            .unwrap();
        let wrapped = format!("<math>{}</math>", source.mathml);
        let document = Document::parse(&wrapped).unwrap();
        assert_eq!(
            document
                .root_element()
                .first_element_child()
                .unwrap()
                .text(),
            Some(text)
        );
        assert_eq!(source.typst, Some(LabelSource::typst_text(text)));
        assert!(!source.mathml.contains("<import"));
    }

    #[test]
    fn typed_math_uses_shared_greek_and_script_rendering() {
        for (superscript, glyph) in [("Gamma", "Gamma"), ("Γ", "Γ"), ("mu", "μ")] {
            let value = NativeValue::Math(MathValue {
                symbol: "alpha".into(),
                subscript: Some(MathScript::Integer(-12)),
                superscript: Some(MathScript::Symbol(superscript.into())),
            });
            let source = value
                .label_source(LabelSource::text(""), &BTreeMap::new())
                .unwrap();
            assert_eq!(
                source.mathml,
                format!("<msubsup><mi>α</mi><mn>-12</mn><mi>{glyph}</mi></msubsup>")
            );
            assert_eq!(source.typst, Some(format!("$alpha_(-12)^({superscript})$")));
        }
    }

    #[test]
    fn typed_text_styles_preserve_color_size_weight_and_slant() {
        let source = NativeValue::Text(TextValue {
            text: "styled <&>".into(),
            options: BTreeMap::from([
                (
                    "fill".into(),
                    NativeValue::Color(ColorValue::Rgba(20, 40, 60, Some(128))),
                ),
                ("size".into(), NativeValue::Length(2.0, LengthUnit::Em)),
                ("weight".into(), NativeValue::String("semibold".into())),
                (
                    "style".into(),
                    NativeValue::Enum(NativeEnum::TextStyle(PyTextStyle::Oblique)),
                ),
            ]),
        })
        .label_source(LabelSource::text(""), &BTreeMap::new())
        .unwrap();
        let wrapped = format!("<math>{}</math>", source.mathml);
        let document = Document::parse(&wrapped).unwrap();
        let style = document.root_element().first_element_child().unwrap();
        let css = style.attribute("style").unwrap();
        for property in [
            "color:#14283c80",
            "font-size:2em",
            "font-weight:600",
            "font-style:oblique",
        ] {
            assert!(css.contains(property), "{css}");
        }
        assert_eq!(
            style.first_element_child().unwrap().text(),
            Some("styled <&>")
        );
        assert_eq!(
            source.typst.as_deref(),
            Some(
                r#"text("styled <&>", fill: rgb(20, 40, 60, 128), size: 2.0em, style: "oblique", weight: "semibold")"#
            )
        );
    }

    #[test]
    fn absolute_text_sizes_share_probe_and_svg_scale_without_changing_typst() {
        for (value, points, typst_size) in [
            (NativeValue::Int(18), 18.0, "18"),
            (NativeValue::Float(18.0), 18.0, "18.0"),
            (NativeValue::Length(18.0, LengthUnit::Pt), 18.0, "18.0pt"),
            (NativeValue::Length(25.4, LengthUnit::Mm), 72.0, "25.4mm"),
            (NativeValue::Length(2.54, LengthUnit::Cm), 72.0, "2.54cm"),
            (NativeValue::Length(1.0, LengthUnit::In), 72.0, "1.0in"),
        ] {
            assert!((value.number().unwrap() - points).abs() < 1e-12);
            let source = TextValue::style_label(
                LabelSource::text("x"),
                &[&BTreeMap::from([("size".into(), value)])],
            )
            .unwrap();
            let wrapped = format!("<math>{}</math>", source.mathml);
            let document = Document::parse(&wrapped).unwrap();
            let css = document
                .root_element()
                .first_element_child()
                .unwrap()
                .attribute("style")
                .unwrap();
            let em: f64 = css
                .strip_prefix("font-size:")
                .unwrap()
                .strip_suffix("em")
                .unwrap()
                .parse()
                .unwrap();
            assert!((12.0 * em * 0.75 - points).abs() < 1e-12, "{css}");
            assert!((9.0 * em - points).abs() < 1e-12, "{css}");
            assert_eq!(
                source.typst,
                Some(format!("text(\"x\", size: {typst_size})"))
            );
        }
    }

    #[test]
    fn relative_text_sizes_stay_relative_and_invalid_sizes_are_rejected() {
        for (size, typst_size) in [(0.0, "0.0em"), (0.5, "0.5em"), (2.0, "2.0em")] {
            let source = TextValue::style_label(
                LabelSource::text("x"),
                &[&BTreeMap::from([(
                    "size".into(),
                    NativeValue::Length(size, LengthUnit::Em),
                )])],
            )
            .unwrap();
            assert!(source.mathml.contains(&format!("font-size:{size}em")));
            assert_eq!(
                source.typst,
                Some(format!("text(\"x\", size: {typst_size})"))
            );
        }
        for size in [-1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for unit in [
                LengthUnit::Pt,
                LengthUnit::Mm,
                LengthUnit::Cm,
                LengthUnit::In,
                LengthUnit::Em,
            ] {
                assert!(TextValue::style_label(
                    LabelSource::text("x"),
                    &[&BTreeMap::from([(
                        "size".into(),
                        NativeValue::Length(size, unit)
                    )])],
                )
                .is_err());
            }
            assert!(TextValue::style_label(
                LabelSource::text("x"),
                &[&BTreeMap::from([("size".into(), NativeValue::Float(size))])],
            )
            .is_err());
        }
        assert!(TextValue::style_label(
            LabelSource::text("x"),
            &[&BTreeMap::from([("size".into(), NativeValue::Int(-1))])],
        )
        .is_err());
        assert!(TextValue::style_label(
            LabelSource::text("x"),
            &[&BTreeMap::from([(
                "size".into(),
                NativeValue::Length(f64::MAX, LengthUnit::In)
            )])],
        )
        .is_err());
        assert!(TextValue::style_label(
            LabelSource::text("x"),
            &[
                &BTreeMap::from([("size".into(), NativeValue::Length(f64::MAX, LengthUnit::Em))]),
                &BTreeMap::from([("size".into(), NativeValue::Length(2.0, LengthUnit::Em))]),
            ],
        )
        .is_err());
    }

    #[test]
    fn nested_typed_sizes_resolve_absolutes_and_compose_relative_em() {
        for (outer, inner, em, typst) in [
            (
                NativeValue::Length(18.0, LengthUnit::Pt),
                NativeValue::Length(18.0, LengthUnit::Pt),
                2.0,
                "text(text(\"X\", size: 18.0pt), size: 18.0pt)",
            ),
            (
                NativeValue::Length(27.0, LengthUnit::Pt),
                NativeValue::Length(18.0, LengthUnit::Pt),
                2.0,
                "text(text(\"X\", size: 18.0pt), size: 27.0pt)",
            ),
            (
                NativeValue::Length(2.0, LengthUnit::Em),
                NativeValue::Length(18.0, LengthUnit::Pt),
                2.0,
                "text(text(\"X\", size: 18.0pt), size: 2.0em)",
            ),
            (
                NativeValue::Length(18.0, LengthUnit::Pt),
                NativeValue::Length(2.0, LengthUnit::Em),
                4.0,
                "text(text(\"X\", size: 2.0em), size: 18.0pt)",
            ),
            (
                NativeValue::Length(2.0, LengthUnit::Em),
                NativeValue::Length(2.0, LengthUnit::Em),
                4.0,
                "text(text(\"X\", size: 2.0em), size: 2.0em)",
            ),
            (
                NativeValue::Length(18.0, LengthUnit::Pt),
                NativeValue::Length(0.5, LengthUnit::Em),
                1.0,
                "text(text(\"X\", size: 0.5em), size: 18.0pt)",
            ),
            (
                NativeValue::Length(0.0, LengthUnit::Em),
                NativeValue::Length(18.0, LengthUnit::Pt),
                2.0,
                "text(text(\"X\", size: 18.0pt), size: 0.0em)",
            ),
            (
                NativeValue::Length(0.0, LengthUnit::Pt),
                NativeValue::Length(2.0, LengthUnit::Em),
                0.0,
                "text(text(\"X\", size: 2.0em), size: 0.0pt)",
            ),
            (
                NativeValue::Length(2.0, LengthUnit::Em),
                NativeValue::Inherit,
                2.0,
                "text(\"X\", size: 2.0em)",
            ),
        ] {
            let label = NativeValue::Text(TextValue {
                text: "X".into(),
                options: BTreeMap::from([("size".into(), inner)]),
            });
            let mut pages = Vec::new();
            assert_eq!(
                label
                    .label(
                        &mut pages,
                        LabelSource::text("unused"),
                        &BTreeMap::from([("size".into(), outer)]),
                    )
                    .unwrap(),
                Some(0)
            );
            assert_eq!(
                pages[0].mathml,
                format!("<mstyle style=\"font-size:{em}em\"><mtext>X</mtext></mstyle>")
            );
            assert_eq!(pages[0].typst.as_deref(), Some(typst));
        }
    }

    #[test]
    fn global_and_local_label_styles_respect_inline_absolute_sizes() {
        let graph: TypstGraphSpec = serde_json::from_value(json!({
            "nodes": [{}, {}, {}],
            "edges": [{"source": {"node": 0, "id": 0}}],
        }))
        .unwrap();
        let label = NativeValue::Text(TextValue {
            text: "X".into(),
            options: BTreeMap::from([("size".into(), NativeValue::Length(18.0, LengthUnit::Pt))]),
        });
        let local = NativeValue::Dict(BTreeMap::from([(
            "size".into(),
            NativeValue::Length(27.0, LengthUnit::Pt),
        )]));
        let global = NativeValue::Dict(BTreeMap::from([(
            "size".into(),
            NativeValue::Length(18.0, LengthUnit::Pt),
        )]));
        let transport = RenderSettingsTransport {
            native: NativeValue::Dict(BTreeMap::from([
                (
                    "elements".into(),
                    NativeValue::Dict(BTreeMap::from([
                        ("graph".into(), NativeValue::None),
                        (
                            "nodes".into(),
                            NativeValue::Array(vec![
                                NativeValue::Dict(BTreeMap::from([(
                                    "label".into(),
                                    label.clone(),
                                )])),
                                NativeValue::Dict(BTreeMap::from([
                                    ("label".into(), label.clone()),
                                    ("node-label-style".into(), local.clone()),
                                ])),
                                NativeValue::Dict(BTreeMap::from([
                                    ("label".into(), NativeValue::String("X".into())),
                                    ("node-label-style".into(), local.clone()),
                                ])),
                            ]),
                        ),
                        (
                            "edges".into(),
                            NativeValue::Array(vec![NativeValue::Dict(BTreeMap::from([
                                ("label".into(), label),
                                ("edge-label-style".into(), local),
                            ]))]),
                        ),
                        ("hedges".into(), NativeValue::Array(Vec::new())),
                    ])),
                ),
                (
                    "style".into(),
                    NativeValue::Dict(BTreeMap::from([
                        ("node-label-style".into(), global.clone()),
                        ("edge-label-style".into(), global),
                    ])),
                ),
            ])),
            imports: Default::default(),
            selectors: SelectorCallbacks {
                node: None,
                edge: None,
                source: None,
                sink: None,
            },
        };
        let scene = transport.scene(graph, None, None).unwrap();
        for (page, em, typst) in [
            (
                scene.nodes[0].label.unwrap(),
                2.0,
                "text(text(\"X\", size: 18.0pt), size: 18.0pt)",
            ),
            (
                scene.nodes[1].label.unwrap(),
                2.0,
                "text(text(\"X\", size: 18.0pt), size: 27.0pt)",
            ),
            (
                scene.nodes[2].label.unwrap(),
                3.0,
                "text(\"X\", size: 27.0pt)",
            ),
            (
                scene.edges[0].label.unwrap(),
                2.0,
                "text(text(\"X\", size: 18.0pt), size: 27.0pt)",
            ),
        ] {
            assert_eq!(
                scene.pages[page].mathml,
                format!("<mstyle style=\"font-size:{em}em\"><mtext>X</mtext></mstyle>")
            );
            assert_eq!(scene.pages[page].typst.as_deref(), Some(typst));
        }
    }

    #[test]
    fn automatic_scene_labels_are_mathml_and_imports_are_rejected() {
        let graph: TypstGraphSpec = serde_json::from_value(json!({
            "nodes": [{}],
            "edges": [{"source": {"node": 0, "id": 0}}],
        }))
        .unwrap();
        let mut transport = RenderSettingsTransport {
            native: NativeValue::Dict(BTreeMap::from([
                (
                    "elements".into(),
                    NativeValue::Dict(BTreeMap::from([
                        ("graph".into(), NativeValue::None),
                        ("nodes".into(), NativeValue::Array(vec![NativeValue::None])),
                        ("edges".into(), NativeValue::Array(vec![NativeValue::None])),
                        ("hedges".into(), NativeValue::Array(Vec::new())),
                    ])),
                ),
                ("title".into(), NativeValue::String("A <&> title".into())),
            ])),
            imports: Default::default(),
            selectors: SelectorCallbacks {
                node: None,
                edge: None,
                source: None,
                sink: None,
            },
        };
        let scene = transport.scene(graph.clone(), None, None).unwrap();
        for (index, letter) in ["n", "e"].into_iter().enumerate() {
            assert_eq!(
                scene.pages[index].mathml,
                format!("<msub><mi>{letter}</mi><mn>0</mn></msub>")
            );
            assert_eq!(
                scene.pages[index].typst.as_deref(),
                Some(format!("${letter}_(0)$").as_str())
            );
        }
        assert_eq!(
            scene.title.unwrap().mathml,
            "<mtext>A &lt;&amp;&gt; title</mtext>"
        );
        transport
            .imports
            .insert(TypstModuleSource::File("evil.typ".into()));
        assert!(transport.scene(graph, None, None).is_err());
    }
}
