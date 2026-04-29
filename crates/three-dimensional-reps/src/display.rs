use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Orientation};
use nu_ansi_term::{Color, Style as AnsiStyle};
use tabled::{builder::Builder, settings::Style};

use crate::{
    GraphInfo, HybridSurfaceID, LinearSurfaceKind, NumeratorSamplingScaleMode,
    OrientationExpression, OrientationID, RepresentationMode, ThreeDExpression,
};

#[derive(Debug, Clone, Default)]
pub struct DisplayOptions {
    pub use_color: bool,
    pub details_for_orientation: Option<String>,
}

pub fn render_expression_summary(
    expression: &ThreeDExpression<OrientationID>,
    family: RepresentationMode,
    graph: &GraphInfo,
    energy_degree_bounds: &[(usize, usize)],
    numerator_expr: Option<&str>,
    sampling_scale: NumeratorSamplingScaleMode,
    options: &DisplayOptions,
) -> String {
    let mut out = Vec::new();
    out.push(title(
        &format!("{} structure", representation_label(family).to_uppercase()),
        options.use_color,
    ));
    out.push(summary_table(
        expression,
        family,
        graph,
        energy_degree_bounds,
        numerator_expr,
        sampling_scale,
        options.use_color,
    ));
    out.push(surface_table(expression, options.use_color));
    out.push(orientation_table(expression, options.use_color));

    if let Some(selector) = options.details_for_orientation.as_deref() {
        let details = orientation_details(expression, selector, options.use_color);
        if !details.is_empty() {
            out.push(details);
        }
    }

    out.join("\n\n")
}

fn summary_table(
    expression: &ThreeDExpression<OrientationID>,
    family: RepresentationMode,
    graph: &GraphInfo,
    energy_degree_bounds: &[(usize, usize)],
    numerator_expr: Option<&str>,
    sampling_scale: NumeratorSamplingScaleMode,
    use_color: bool,
) -> String {
    let mut table = Builder::new();
    table.push_record(vec![h("field", use_color), h("value", use_color)]);
    table.push_record(vec![
        "family".to_string(),
        c(representation_label(family), Color::Green, use_color),
    ]);
    table.push_record(vec![
        "internal edges".to_string(),
        graph.n_internal_edges.to_string(),
    ]);
    table.push_record(vec![
        "external edges".to_string(),
        graph.n_external_edges.to_string(),
    ]);
    table.push_record(vec![
        "external symbols".to_string(),
        graph.ext_names.len().to_string(),
    ]);
    table.push_record(vec![
        "linear surfaces".to_string(),
        expression.surfaces.linear_surface_cache.len().to_string(),
    ]);
    table.push_record(vec![
        "orientations".to_string(),
        expression.orientations.len().to_string(),
    ]);
    table.push_record(vec![
        "variants".to_string(),
        expression
            .orientations
            .iter()
            .map(|orientation| orientation.variants.len())
            .sum::<usize>()
            .to_string(),
    ]);
    table.push_record(vec![
        "repeated groups".to_string(),
        graph
            .repeated_groups
            .iter()
            .map(|group| format!("[{}]", group.iter().join(",")))
            .join(" "),
    ]);
    table.push_record(vec![
        "energy degree bounds".to_string(),
        format_energy_degree_bounds(energy_degree_bounds),
    ]);
    table.push_record(vec![
        "numerator sampling scale".to_string(),
        sampling_scale_label(sampling_scale).to_string(),
    ]);
    table.push_record(vec![
        "numerator".to_string(),
        numerator_expr.unwrap_or("-").to_string(),
    ]);
    table.build().with(Style::rounded()).to_string()
}

fn surface_table(expression: &ThreeDExpression<OrientationID>, use_color: bool) -> String {
    let mut table = Builder::new();
    table.push_record(vec![
        h("sid", use_color),
        h("class", use_color),
        h("origin", use_color),
        h("num only", use_color),
        h("expr", use_color),
    ]);

    if expression.surfaces.linear_surface_cache.is_empty() {
        table.push_record(vec![
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
        ]);
    } else {
        for (id, surface) in expression.surfaces.linear_surface_cache.iter().enumerate() {
            table.push_record(vec![
                c(&format!("L{id}"), Color::Green, use_color),
                c(
                    match surface.kind {
                        LinearSurfaceKind::Esurface => "e",
                        LinearSurfaceKind::Hsurface => "h",
                    },
                    Color::Purple,
                    use_color,
                ),
                format!("{:?}", surface.origin).to_lowercase(),
                surface.numerator_only.to_string(),
                surface.expression.to_string(),
            ]);
        }
    }

    format!(
        "{}\n{}",
        title("Surfaces", use_color),
        table.build().with(Style::rounded())
    )
}

fn orientation_table(expression: &ThreeDExpression<OrientationID>, use_color: bool) -> String {
    let mut table = Builder::new();
    table.push_record(vec![
        h("id", use_color),
        h("orient", use_color),
        h("variant", use_color),
        h("origin", use_color),
        h("pref", use_color),
        h("half edges", use_color),
        h("M power", use_color),
        h("num surfaces", use_color),
        h("den tree", use_color),
    ]);

    if expression.orientations.is_empty() {
        table.push_record(vec![
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
            "-".to_string(),
        ]);
    }

    for (orientation_id, orientation) in expression.orientations.iter().enumerate() {
        let label = orientation_label(orientation);
        for (variant_id, variant) in orientation.variants.iter().enumerate() {
            table.push_record(vec![
                if variant_id == 0 {
                    c(&orientation_id.to_string(), Color::Green, use_color)
                } else {
                    String::new()
                },
                if variant_id == 0 {
                    color_orientation_label(&label, use_color)
                } else {
                    String::new()
                },
                variant_id.to_string(),
                variant.origin.as_deref().unwrap_or("term").to_string(),
                c(&variant.prefactor.to_string(), Color::Yellow, use_color),
                format_edge_list(&variant.half_edges),
                inverse_scale_text(variant.uniform_scale_power),
                format_surface_list(&variant.numerator_surfaces),
                truncate_for_table(&tree_shape(&variant.denominator), 140),
            ]);
        }
    }

    format!(
        "{}\n{}",
        title("Orientations", use_color),
        table.build().with(Style::rounded())
    )
}

fn orientation_details(
    expression: &ThreeDExpression<OrientationID>,
    selector: &str,
    use_color: bool,
) -> String {
    let matching = expression
        .orientations
        .iter()
        .enumerate()
        .filter(|(id, orientation)| orientation_matches(*id, orientation, selector))
        .collect_vec();
    if matching.is_empty() {
        return format!(
            "{}\nNo orientation matched `{selector}`.",
            title("Details", use_color)
        );
    }

    let mut sections = Vec::new();
    sections.push(title("Details", use_color));
    for (id, orientation) in matching {
        sections.push(format!(
            "{} {}",
            c("Orientation", Color::Blue, use_color),
            c(
                &format!("{id}:{}", orientation_label(orientation)),
                Color::Green,
                use_color
            )
        ));
        sections.push(energy_map_table(
            "loop q0 map",
            orientation
                .loop_energy_map
                .iter()
                .enumerate()
                .map(|(i, expr)| (format!("k{i}"), expr.to_string())),
            use_color,
        ));
        sections.push(energy_map_table(
            "edge q0 map",
            orientation
                .edge_energy_map
                .iter()
                .enumerate()
                .map(|(i, expr)| (format!("q{i}"), expr.to_string())),
            use_color,
        ));
        for (variant_id, variant) in orientation.variants.iter().enumerate() {
            sections.push(format!(
                "{} {}  origin={}  pref={}  half_edges={}{}",
                c("variant", Color::Blue, use_color),
                variant_id,
                variant.origin.as_deref().unwrap_or("term"),
                variant.prefactor,
                format_edge_list(&variant.half_edges),
                if variant.uniform_scale_power == 0 {
                    String::new()
                } else {
                    format!("  {}", inverse_scale_text(variant.uniform_scale_power))
                }
            ));
            sections.push(tree_node_table(&variant.denominator, use_color));
        }
    }

    sections.join("\n")
}

fn energy_map_table(
    title_text: &str,
    rows: impl IntoIterator<Item = (String, String)>,
    use_color: bool,
) -> String {
    let mut table = Builder::new();
    table.push_record(vec![h("id", use_color), h("expr", use_color)]);
    let mut count = 0usize;
    for (id, expr) in rows {
        table.push_record(vec![id, expr]);
        count += 1;
    }
    if count == 0 {
        table.push_record(vec!["-".to_string(), "-".to_string()]);
    }
    format!(
        "{}\n{}",
        c(title_text, Color::Blue, use_color),
        table.build().with(Style::rounded())
    )
}

fn tree_node_table(tree: &crate::tree::Tree<HybridSurfaceID>, use_color: bool) -> String {
    let mut table = Builder::new();
    table.push_record(vec![
        h("node", use_color),
        h("surface", use_color),
        h("children", use_color),
    ]);
    for node in tree.iter_nodes() {
        table.push_record(vec![
            node.node_id.0.to_string(),
            format_surface_id(node.data),
            format!(
                "[{}]",
                node.children
                    .iter()
                    .map(|child| child.0.to_string())
                    .join(",")
            ),
        ]);
    }
    table.build().with(Style::rounded()).to_string()
}

fn title(text: &str, use_color: bool) -> String {
    c(text, Color::Cyan, use_color)
}

fn h(text: &str, use_color: bool) -> String {
    if use_color {
        AnsiStyle::new()
            .bold()
            .paint(Color::Cyan.paint(text).to_string())
            .to_string()
    } else {
        text.to_string()
    }
}

fn c(text: &str, color: Color, use_color: bool) -> String {
    if use_color {
        color.paint(text).to_string()
    } else {
        text.to_string()
    }
}

fn representation_label(family: RepresentationMode) -> &'static str {
    match family {
        RepresentationMode::Ltd => "ltd",
        RepresentationMode::Cff => "cff",
        #[cfg(any(feature = "diagnostics", feature = "test-support"))]
        RepresentationMode::PureLtd => "pure_ltd",
        #[cfg(feature = "old-cff")]
        RepresentationMode::OldCff => "old_cff",
    }
}

fn sampling_scale_label(mode: NumeratorSamplingScaleMode) -> &'static str {
    match mode {
        NumeratorSamplingScaleMode::None => "none",
        NumeratorSamplingScaleMode::BeyondQuadratic => "beyond-quadratic",
        NumeratorSamplingScaleMode::All => "all",
    }
}

fn format_energy_degree_bounds(bounds: &[(usize, usize)]) -> String {
    if bounds.is_empty() {
        "-".to_string()
    } else {
        bounds
            .iter()
            .map(|(edge, degree)| format!("{edge}:{degree}"))
            .join(",")
    }
}

fn orientation_label(orientation: &OrientationExpression) -> String {
    orientation.data.label.clone().unwrap_or_else(|| {
        orientation
            .data
            .orientation
            .iter()
            .map(|(_, orientation)| match orientation {
                Orientation::Default => '+',
                Orientation::Reversed => '-',
                Orientation::Undirected => 'x',
            })
            .collect()
    })
}

fn color_orientation_label(label: &str, use_color: bool) -> String {
    if !use_color {
        return label.to_string();
    }
    label
        .chars()
        .map(|ch| match ch {
            '+' => Color::Green.paint(ch.to_string()).to_string(),
            '-' => Color::Red.paint(ch.to_string()).to_string(),
            '0' => Color::Yellow.paint(ch.to_string()).to_string(),
            'x' => Color::Blue.paint(ch.to_string()).to_string(),
            _ => ch.to_string(),
        })
        .join("")
}

fn format_edge_list(edges: &[EdgeIndex]) -> String {
    if edges.is_empty() {
        "[]".to_string()
    } else {
        format!(
            "[{}]",
            edges.iter().map(|edge| edge.0.to_string()).join(",")
        )
    }
}

fn format_surface_list(surfaces: &[HybridSurfaceID]) -> String {
    if surfaces.is_empty() {
        "[]".to_string()
    } else {
        format!(
            "[{}]",
            surfaces.iter().copied().map(format_surface_id).join(",")
        )
    }
}

fn format_surface_id(surface_id: HybridSurfaceID) -> String {
    match surface_id {
        HybridSurfaceID::Esurface(id) => format!("E{}", id.0),
        HybridSurfaceID::Hsurface(id) => format!("H{}", id.0),
        HybridSurfaceID::Linear(id) => format!("L{}", id.0),
        HybridSurfaceID::Unit => "1".to_string(),
        HybridSurfaceID::Infinite => "inf".to_string(),
    }
}

fn inverse_scale_text(power: usize) -> String {
    if power == 0 {
        "-".to_string()
    } else {
        format!("M^-{power}")
    }
}

fn tree_shape(tree: &crate::tree::Tree<HybridSurfaceID>) -> String {
    tree.iter_nodes()
        .map(|node| {
            format!(
                "{}:{}->[{}]",
                node.node_id.0,
                format_surface_id(node.data),
                node.children
                    .iter()
                    .map(|child| child.0.to_string())
                    .join(",")
            )
        })
        .join("; ")
}

fn truncate_for_table(text: &str, max_chars: usize) -> String {
    if text.chars().count() <= max_chars {
        return text.to_string();
    }
    let mut truncated = text
        .chars()
        .take(max_chars.saturating_sub(3))
        .collect::<String>();
    truncated.push_str("...");
    truncated
}

fn orientation_matches(id: usize, orientation: &OrientationExpression, selector: &str) -> bool {
    if selector == id.to_string() {
        return true;
    }
    let label = orientation_label(orientation);
    if label == selector {
        return true;
    }
    label
        .split_once('|')
        .is_some_and(|(base, _)| base == selector)
}

#[cfg(test)]
mod tests {
    use crate::{
        Generate3DExpressionOptions, RepresentationMode,
        generation::generate_3d_expression_from_parsed, graph_io::graph_info,
    };

    use super::*;

    #[test]
    fn display_renders_summary_and_orientation_details() {
        let parsed = crate::graph_io::test_graphs::box_graph();
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                ..Default::default()
            },
        )
        .unwrap();

        let rendered = render_expression_summary(
            &expression,
            RepresentationMode::Cff,
            &graph_info(&parsed),
            &[],
            None,
            NumeratorSamplingScaleMode::None,
            &DisplayOptions {
                use_color: false,
                details_for_orientation: Some("0".to_string()),
            },
        );

        assert!(rendered.contains("CFF structure"));
        assert!(rendered.contains("Surfaces"));
        assert!(rendered.contains("Orientations"));
        assert!(rendered.contains("loop q0 map"));
        assert!(rendered.contains("edge q0 map"));
    }
}
