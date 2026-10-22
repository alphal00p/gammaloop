//! Convert typed graph drawings to native SVG scenes.
use super::{
    LengthUnit, MathScript, MathValue, NativeEnum, NativeValue, RenderSettingsTransport, TextValue,
};
use crate::drawing::{DrawingKind, EDGE_FIELDS, NODE_FIELDS};
use linnest::svg::{Config, Dash, EdgeDrawing, LabelSource, NodeDrawing, Scene, Stroke};
use linnest::TypstGraphSpec;
use spenso::structure::representation::IndexDisplay;
use std::collections::BTreeMap;
use symbolica::printer::AnsiHtmlFormatter;

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
    ) -> Option<LabelSource> {
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
                    "fill" => format!("color:{}", value.json().ok()?.as_str()?),
                    "size" => {
                        // The 12px probe and 9-unit SVG label share a 9pt base.
                        let size = match value {
                            NativeValue::Length(value, LengthUnit::Em) => *value,
                            _ => value.number()? / 9.0,
                        };
                        if !size.is_finite() || size < 0.0 {
                            return None;
                        }
                        em = if matches!(value, NativeValue::Length(_, LengthUnit::Em)) {
                            em * size
                        } else {
                            size
                        };
                        if !em.is_finite() {
                            return None;
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
                                _ => return None,
                            },
                            _ => return None,
                        };
                        format!("font-weight:{weight}")
                    }
                    "style" => {
                        let NativeValue::Enum(NativeEnum::TextStyle(style)) = value else {
                            return None;
                        };
                        format!("font-style:{}", style.typst_name())
                    }
                    _ => return None,
                };
                css.insert(key, property);
                group.push(format!("{key}: {}", value.source(&BTreeMap::new()).ok()?));
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
        Some(source)
    }
}

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

    fn label_source(
        &self,
        default: LabelSource,
        style: &BTreeMap<String, NativeValue>,
    ) -> Option<LabelSource> {
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
                mathml: value.display().to_mathml()?,
                typst: Some(self.source(&BTreeMap::new()).ok()?),
            },
            _ => return None,
        };
        TextValue::style_label(source, &[style])
    }

    fn label(
        &self,
        pages: &mut Vec<LabelSource>,
        default: LabelSource,
        style: &BTreeMap<String, NativeValue>,
    ) -> Option<Option<usize>> {
        if matches!(self, Self::None) {
            return Some(None);
        }
        let source = self.label_source(default, style)?;
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
                scene.title = Some(title.label_source(LabelSource::text(""), &BTreeMap::new())?);
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
                match default_nodes.get(index)?.label {
                    Some(page) if !label_style.is_empty() => {
                        let source = TextValue::style_label(
                            scene.pages.get(page)?.clone(),
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
                    subscript: Some(MathScript::Integer(index.try_into().ok()?)),
                    superscript: None,
                });
                label.label(
                    &mut scene.pages,
                    default.label_source(LabelSource::text(""), &BTreeMap::new())?,
                    &label_style,
                )?
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
                match default_edges.get(index)?.label {
                    Some(page) if !label_style.is_empty() => {
                        let source = TextValue::style_label(
                            scene.pages.get(page)?.clone(),
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
                    subscript: Some(MathScript::Integer(index.try_into().ok()?)),
                    superscript: None,
                });
                label.label(
                    &mut scene.pages,
                    default.label_source(LabelSource::text(""), &BTreeMap::new())?,
                    &label_style,
                )?
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::typst::{ColorValue, PyTextStyle, SelectorCallbacks, TypstModuleSource};

    #[test]
    fn literal_labels_escape_mathml_and_keep_literal_typst() {
        let text = "<&>\"' #import \"evil.typ\"\n\\";
        let source = NativeValue::String(text.into())
            .label_source(LabelSource::text("unused"), &BTreeMap::new())
            .unwrap();
        let wrapped = format!("<math>{}</math>", source.mathml);
        let document = roxmltree::Document::parse(&wrapped).unwrap();
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
        let document = roxmltree::Document::parse(&wrapped).unwrap();
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
            let document = roxmltree::Document::parse(&wrapped).unwrap();
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
                .is_none());
            }
            assert!(TextValue::style_label(
                LabelSource::text("x"),
                &[&BTreeMap::from([("size".into(), NativeValue::Float(size))])],
            )
            .is_none());
        }
        assert!(TextValue::style_label(
            LabelSource::text("x"),
            &[&BTreeMap::from([("size".into(), NativeValue::Int(-1))])],
        )
        .is_none());
        assert!(TextValue::style_label(
            LabelSource::text("x"),
            &[&BTreeMap::from([(
                "size".into(),
                NativeValue::Length(f64::MAX, LengthUnit::In)
            )])],
        )
        .is_none());
        assert!(TextValue::style_label(
            LabelSource::text("x"),
            &[
                &BTreeMap::from([("size".into(), NativeValue::Length(f64::MAX, LengthUnit::Em))]),
                &BTreeMap::from([("size".into(), NativeValue::Length(2.0, LengthUnit::Em))]),
            ],
        )
        .is_none());
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
                label.label(
                    &mut pages,
                    LabelSource::text("unused"),
                    &BTreeMap::from([("size".into(), outer)]),
                ),
                Some(Some(0))
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
        let graph: TypstGraphSpec = serde_json::from_value(serde_json::json!({
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
        let graph: TypstGraphSpec = serde_json::from_value(serde_json::json!({
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
        assert!(transport.scene(graph, None, None).is_none());
    }
}
