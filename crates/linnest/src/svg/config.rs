use super::{Dash, Scene, Stroke};
use serde::Deserialize;
use serde_json::{Map, Value};

/// Serializable drawing options for native notebook graphs. No executable
/// Typst templates or callbacks cross this boundary.
#[derive(Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Config {
    pub title: Option<String>,
    pub layouts: Map<String, Value>,
    pub drawing: Map<String, Value>,
    pub style: Map<String, Value>,
    pub template_options: Map<String, Value>,
}

impl Config {
    pub fn from_json(source: &str) -> Result<Self, String> {
        serde_json::from_str(source).map_err(|e| format!("invalid SVG configuration: {e}"))
    }

    pub fn apply(&self, scene: &mut Scene) -> Result<(), String> {
        if let Some(title) = &self.title {
            scene.title = Some(format!("#{title:?}"));
        }
        for (key, value) in &self.layouts {
            let key = key.replace('_', "-");
            let key = if key == "steps" {
                "impred-steps".to_owned()
            } else {
                key
            };
            if ![
                "layout-algo",
                "tree-dx",
                "tree-dy",
                "impred-steps",
                "impred-step-scale",
                "impred-spacing",
                "impred-repulsion",
                "impred-attraction",
                "impred-parallel-balance",
                "impred-pull",
                "impred-pull-balance",
                "impred-external-max-points",
                "impred-split-length-ratio",
                "impred-contract-chord-ratio",
                "impred-edge-clearance",
                "impred-node-edge-strength",
                "impred-labels",
                "impred-level",
                "impred-label-room",
                "impred-label-room-margin",
                "impred-label-room-max-scale",
            ]
            .contains(&key.as_str())
            {
                return Err(format!(
                    "layout option {key:?} is not supported by the native SVG renderer"
                ));
            }
            if key == "layout-algo" && !matches!(value.as_str(), Some("dot" | "impred")) {
                return Err("layout-algo must be dot or impred".into());
            }
            if matches!(key.as_str(), "tree-dx" | "tree-dy")
                && !value.as_f64().is_some_and(|v| v.is_finite() && v > 0.0)
            {
                return Err(format!("{key} must be positive"));
            }
            // Native label-aware layouts refine around the drawn labels.
            if key == "impred-labels" {
                scene.label_feedback = value.as_bool().ok_or("impred-labels must be boolean")?;
                continue;
            }
            scene.layout.insert(key, value.clone());
        }
        for (key, value) in &self.drawing {
            match key.replace('_', "-").as_str() {
                "label-leaders" => {
                    scene.label_leaders = value.as_bool().ok_or("label-leaders must be boolean")?;
                }
                "node-radius" => {
                    let radius = value
                        .as_f64()
                        .filter(|r| r.is_finite() && *r >= 0.0)
                        .ok_or("node-radius must be nonnegative")?;
                    for node in &mut scene.nodes {
                        node.radius = radius;
                    }
                }
                _ => {
                    return Err(format!(
                        "drawing option {key:?} is not supported by the native SVG renderer"
                    ));
                }
            }
        }
        for (key, value) in &self.style {
            let fields = value
                .as_object()
                .ok_or("style groups must be dictionaries")?;
            match key.replace('_', "-").as_str() {
                "node-style" => {
                    for (key, value) in fields {
                        match key.as_str() {
                            "stroke" => {
                                let stroke = StrokePatch::parse(value)?;
                                for node in &mut scene.nodes {
                                    stroke.apply(&mut node.stroke);
                                }
                            }
                            "fill" => {
                                let fill = value.as_str().ok_or("fill must be a CSS color")?;
                                for node in &mut scene.nodes {
                                    fill.clone_into(&mut node.fill);
                                }
                            }
                            "radius" => {
                                let radius = value
                                    .as_f64()
                                    .filter(|r| r.is_finite() && *r >= 0.0)
                                    .ok_or("radius must be nonnegative")?;
                                for node in &mut scene.nodes {
                                    node.radius = radius;
                                }
                            }
                            _ => return Err(format!("unsupported node style {key:?}")),
                        }
                    }
                }
                "edge-style" => {
                    for (key, value) in fields {
                        match key.as_str() {
                            "stroke" => {
                                let stroke = StrokePatch::parse(value)?;
                                for edge in &mut scene.edges {
                                    stroke.apply(&mut edge.stroke);
                                }
                            }
                            "flow-arrow" | "momentum-arrow" => {
                                let mark = Self::mark(value)?;
                                for edge in &mut scene.edges {
                                    if key == "flow-arrow" {
                                        edge.flow_arrow = mark.clone();
                                    } else {
                                        edge.momentum_arrow = mark.clone();
                                    }
                                }
                            }
                            "flow-arrow-paints" | "momentum-arrow-paints" => {
                                let paints =
                                    serde_json::from_value::<Vec<super::MarkPaint>>(value.clone())
                                        .map_err(|e| format!("invalid mark paints: {e}"))?;
                                for edge in &mut scene.edges {
                                    if key == "flow-arrow-paints" {
                                        edge.flow_arrow_paints = paints.clone();
                                    } else {
                                        edge.momentum_arrow_paints = paints.clone();
                                    }
                                }
                            }
                            _ => return Err(format!("unsupported edge style {key:?}")),
                        }
                    }
                    for edge in &scene.edges {
                        for (mark, paints, stroke) in [
                            (&edge.flow_arrow, &edge.flow_arrow_paints, &edge.stroke),
                            (
                                &edge.momentum_arrow,
                                &edge.momentum_arrow_paints,
                                &super::arrow_stroke(),
                            ),
                        ] {
                            let prepared = mark.prepare(super::marks::Batch::context(stroke))?;
                            if !paints.is_empty() && paints.len() != prepared.paths().len() {
                                return Err(
                                    "mark paints must match the flattened mark path count".into()
                                );
                            }
                        }
                    }
                }
                _ => return Err(format!("unsupported graph style {key:?}")),
            }
        }
        Ok(())
    }

    fn mark(value: &Value) -> Result<kurvst::marks::MarkSpec, String> {
        serde_json::from_value(value.clone()).map_err(|error| {
            format!("invalid mark: {error}; replace CeTZ mark dictionaries with Mark(name, ...) normalized shape data")
        })
    }
}

/// A stroke option: a paint alone, or the fields a stroke dictionary sets.
/// It is parsed once, so an invalid option fails even without elements.
#[derive(Default)]
struct StrokePatch {
    paint: Option<String>,
    width: Option<f64>,
    dash: Option<Dash>,
}

impl StrokePatch {
    fn parse(value: &Value) -> Result<Self, String> {
        if let Some(paint) = value.as_str() {
            return Ok(Self {
                paint: Some(paint.into()),
                ..Self::default()
            });
        }
        let mut patch = Self::default();
        for (key, value) in value
            .as_object()
            .ok_or("stroke must be a color or dictionary")?
        {
            match key.as_str() {
                "paint" => {
                    patch.paint = Some(value.as_str().ok_or("paint must be a CSS color")?.into())
                }
                "thickness" => {
                    patch.width = Some(
                        value
                            .as_f64()
                            .filter(|v| v.is_finite() && *v >= 0.0)
                            .ok_or("thickness must be nonnegative points")?,
                    )
                }
                "dash" => {
                    patch.dash = Some(match value.as_str() {
                        Some("solid") => Dash::Solid,
                        Some("dotted") => Dash::Dotted,
                        Some("dashed") => Dash::Dashed(3.0, 3.0),
                        _ => return Err("dash must be solid, dotted, or dashed".into()),
                    })
                }
                _ => return Err(format!("unsupported stroke field {key:?}")),
            }
        }
        Ok(patch)
    }

    fn apply(&self, stroke: &mut Stroke) {
        if let Some(paint) = &self.paint {
            stroke.paint.clone_from(paint);
        }
        if let Some(width) = self.width {
            stroke.width = width;
        }
        if let Some(dash) = self.dash {
            stroke.dash = dash;
        }
    }
}
