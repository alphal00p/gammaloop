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
            ]
            .contains(&key.as_str())
            {
                return Err(format!(
                    "layout option {key:?} is not supported by the native SVG renderer"
                ));
            }
            if key == "impred-labels" {
                scene.label_feedback = value.as_bool().ok_or("impred-labels must be boolean")?;
            }
            scene.layout.insert(key, value.clone());
        }
        for (key, value) in &self.drawing {
            match key.replace('_', "-").as_str() {
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
                    for node in &mut scene.nodes {
                        for (key, value) in fields {
                            match key.as_str() {
                                "stroke" => Self::stroke(&mut node.stroke, value)?,
                                "fill" => {
                                    node.fill =
                                        value.as_str().ok_or("fill must be a CSS color")?.to_owned()
                                }
                                "radius" => {
                                    node.radius = value
                                        .as_f64()
                                        .filter(|r| r.is_finite() && *r >= 0.0)
                                        .ok_or("radius must be nonnegative")?
                                }
                                _ => return Err(format!("unsupported node style {key:?}")),
                            }
                        }
                    }
                }
                "edge-style" => {
                    for edge in &mut scene.edges {
                        for (key, value) in fields {
                            match key.as_str() {
                                "stroke" => Self::stroke(&mut edge.stroke, value)?,
                                _ => return Err(format!("unsupported edge style {key:?}")),
                            }
                        }
                    }
                }
                _ => return Err(format!("unsupported graph style {key:?}")),
            }
        }
        Ok(())
    }

    fn stroke(stroke: &mut Stroke, value: &Value) -> Result<(), String> {
        if let Some(paint) = value.as_str() {
            stroke.paint = paint.into();
            return Ok(());
        }
        for (key, value) in value
            .as_object()
            .ok_or("stroke must be a color or dictionary")?
        {
            match key.as_str() {
                "paint" => stroke.paint = value.as_str().ok_or("paint must be a CSS color")?.into(),
                "thickness" => {
                    stroke.width = value
                        .as_f64()
                        .filter(|v| v.is_finite() && *v >= 0.0)
                        .ok_or("thickness must be nonnegative points")?
                }
                "dash" => {
                    stroke.dash = match value.as_str() {
                        Some("solid") => Dash::Solid,
                        Some("dotted") => Dash::Dotted,
                        Some("dashed") => Dash::Dashed(3.0, 3.0),
                        _ => return Err("dash must be solid, dotted, or dashed".into()),
                    }
                }
                _ => return Err(format!("unsupported stroke field {key:?}")),
            }
        }
        Ok(())
    }
}
