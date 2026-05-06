use std::collections::{BTreeMap, BTreeSet};

use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Orientation};
use nu_ansi_term::{Color, Style as AnsiStyle};
use symbolica::atom::AtomCore;
use tabled::{builder::Builder, settings::Style};

use crate::{
    GraphInfo, HybridSurfaceID, LinearEnergyExpr, LinearSurfaceKind, NumeratorSamplingScaleMode,
    OrientationExpression, OrientationID, RepresentationMode, ThreeDExpression,
    surface::{RationalAtomExt, SurfaceOrigin},
};

#[derive(Debug, Clone, Default)]
pub struct DisplayOptions {
    pub use_color: bool,
    pub details_for_orientation: Option<String>,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct NumeratorDisplay<'a> {
    pub original: Option<&'a str>,
    pub simplified: Option<&'a str>,
}

pub fn render_expression_summary(
    expression: &ThreeDExpression<OrientationID>,
    family: RepresentationMode,
    graph: &GraphInfo,
    energy_degree_bounds: &[(usize, usize)],
    numerator: NumeratorDisplay<'_>,
    sampling_scale: NumeratorSamplingScaleMode,
    options: &DisplayOptions,
) -> String {
    if let Some(selector) = options.details_for_orientation.as_deref() {
        return orientation_details(expression, selector, graph, options.use_color);
    }

    let mut out = Vec::new();
    out.push(format!(
        "{}\n{}",
        title(
            &format!("{} structure", representation_label(family).to_uppercase()),
            options.use_color,
        ),
        summary_table(
            expression,
            family,
            graph,
            energy_degree_bounds,
            numerator,
            sampling_scale,
            options.use_color,
        )
    ));
    out.push(surface_table(expression, options.use_color));
    out.push(orientation_table(expression, graph, options.use_color));

    out.join("\n\n")
}

fn summary_table(
    expression: &ThreeDExpression<OrientationID>,
    family: RepresentationMode,
    graph: &GraphInfo,
    energy_degree_bounds: &[(usize, usize)],
    numerator: NumeratorDisplay<'_>,
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
        "# nodes".to_string(),
        expression
            .orientations
            .iter()
            .flat_map(|orientation| &orientation.variants)
            .map(|variant| variant.denominator.get_num_nodes())
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
        format_energy_degree_bounds(energy_degree_bounds, use_color),
    ]);
    table.push_record(vec![
        "numerator sampling scale".to_string(),
        sampling_scale_label(sampling_scale).to_string(),
    ]);
    table.push_record(vec![
        "numerator".to_string(),
        numerator.original.unwrap_or("-").to_string(),
    ]);
    table.push_record(vec![
        "simplified_numerator".to_string(),
        numerator.simplified.unwrap_or("-").to_string(),
    ]);
    table.build().with(Style::rounded()).to_string()
}

fn surface_table(expression: &ThreeDExpression<OrientationID>, use_color: bool) -> String {
    let mut table = Builder::new();
    table.push_record(vec![
        h("sid", use_color),
        h("type", use_color),
        h("expr", use_color),
    ]);

    if expression.surfaces.linear_surface_cache.is_empty() {
        table.push_record(vec!["-".to_string(), "-".to_string(), "-".to_string()]);
    } else {
        let groups_by_surface = expression
            .surfaces
            .linear_surface_cache
            .iter()
            .map(|surface| format_linear_energy_groups(&surface.expression, use_color))
            .collect::<Vec<_>>();
        let internal_width = groups_by_surface
            .iter()
            .map(|groups| groups.internal_plain.chars().count())
            .max()
            .unwrap_or(0);
        for ((id, surface), groups) in expression
            .surfaces
            .linear_surface_cache
            .iter()
            .enumerate()
            .zip(groups_by_surface)
        {
            table.push_record(vec![
                surface_id_label(
                    HybridSurfaceID::Linear(crate::surface::LinearSurfaceID(id)),
                    use_color,
                ),
                surface_type_label(
                    surface.kind,
                    surface.origin,
                    surface.numerator_only,
                    use_color,
                ),
                format!(
                    "{}{}{}",
                    groups.internal_colored,
                    " ".repeat(
                        internal_width.saturating_sub(groups.internal_plain.chars().count()) + 3
                    ),
                    groups.external_colored
                )
                .trim_end()
                .to_string(),
            ]);
        }
    }

    format!(
        "{}\n{}",
        title("Surfaces", use_color),
        table.build().with(Style::rounded())
    )
}

fn orientation_table(
    expression: &ThreeDExpression<OrientationID>,
    graph: &GraphInfo,
    use_color: bool,
) -> String {
    let mut table = Builder::new();
    table.push_record(vec![
        h("id", use_color),
        h("orientation", use_color),
        h("variant", use_color),
        h("origin", use_color),
        h("pref", use_color),
        h("factors", use_color),
        h("surfaces", use_color),
        h("# nodes", use_color),
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

    for entry in orientation_display_entries(expression) {
        for (variant_id, variant) in entry.orientation.variants.iter().enumerate() {
            table.push_record(vec![
                if variant_id == 0 {
                    c(&entry.id.to_string(), Color::Green, use_color)
                } else {
                    String::new()
                },
                if variant_id == 0 {
                    color_orientation_label(&entry.label, graph.n_external_edges, use_color)
                } else {
                    String::new()
                },
                variant_id.to_string(),
                origin_label(variant.origin.as_deref().unwrap_or("term"), use_color),
                coefficient_label(&variant.prefactor, use_color),
                factor_list(variant.uniform_scale_power, &variant.half_edges, use_color),
                variant_surface_list(variant, use_color),
                variant.denominator.get_num_nodes().to_string(),
                truncate_for_table(&tree_shape(&variant.denominator, use_color), 30),
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
    graph: &GraphInfo,
    use_color: bool,
) -> String {
    let parsed_selector = DetailSelector::parse(selector);
    let matching = orientation_display_entries(expression)
        .into_iter()
        .filter(|entry| {
            orientation_matches(entry, &parsed_selector)
                && parsed_selector
                    .variant_id
                    .is_none_or(|variant_id| variant_id < entry.orientation.variants.len())
        })
        .collect_vec();
    if matching.is_empty() {
        return format!(
            "{}\nNo orientation matched `{selector}`.",
            title("Details", use_color)
        );
    }

    let mut sections = Vec::new();
    sections.push(title("Details", use_color));
    for entry in matching {
        sections.push(format!(
            "{} {}:{}",
            c("Orientation", Color::Blue, use_color),
            c(&entry.id.to_string(), Color::Green, use_color),
            color_orientation_label(&entry.label, graph.n_external_edges, use_color)
        ));
        sections.push(energy_map_table(
            "loop q0 map",
            entry
                .orientation
                .loop_energy_map
                .iter()
                .enumerate()
                .map(|(i, expr)| (format!("k{i}"), format_linear_energy_expr(expr, use_color))),
            use_color,
        ));
        sections.push(energy_map_table(
            "edge q0 map",
            entry
                .orientation
                .edge_energy_map
                .iter()
                .enumerate()
                .map(|(i, expr)| (format!("q{i}"), format_linear_energy_expr(expr, use_color))),
            use_color,
        ));
        for (variant_id, variant) in
            entry
                .orientation
                .variants
                .iter()
                .enumerate()
                .filter(|(variant_id, _)| {
                    parsed_selector
                        .variant_id
                        .is_none_or(|requested| requested == *variant_id)
                })
        {
            sections.push(format!(
                "{} {}  origin={}  pref={}  factors={}  surfaces={}  #nodes={}",
                c("variant", Color::Blue, use_color),
                variant_id,
                origin_label(variant.origin.as_deref().unwrap_or("term"), use_color),
                coefficient_label(&variant.prefactor, use_color),
                factor_list(variant.uniform_scale_power, &variant.half_edges, use_color),
                variant_surface_list(variant, use_color),
                variant.denominator.get_num_nodes(),
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
            surface_id_label(node.data, use_color),
            format!(
                "[{}]",
                node.children
                    .iter()
                    .map(|child| node_id_label(child.0, use_color))
                    .join(" ")
            ),
        ]);
    }
    table.build().with(Style::rounded()).to_string()
}

struct DetailSelector<'a> {
    orientation_selector: &'a str,
    map_index: Option<usize>,
    variant_id: Option<usize>,
}

impl<'a> DetailSelector<'a> {
    fn parse(selector: &'a str) -> Self {
        let parts = selector.split('|').collect::<Vec<_>>();
        if let [orientation_selector, map_selector, variant_selector] = parts.as_slice()
            && let Some(map_index) = parse_map_selector(map_selector)
            && let Ok(variant_id) = variant_selector.parse::<usize>()
        {
            return Self {
                orientation_selector,
                map_index: Some(map_index),
                variant_id: Some(variant_id),
            };
        }
        if let [orientation_selector, map_selector] = parts.as_slice()
            && let Some(map_index) = parse_map_selector(map_selector)
        {
            return Self {
                orientation_selector,
                map_index: Some(map_index),
                variant_id: None,
            };
        }
        if let [orientation_selector, variant_selector] = parts.as_slice()
            && let Ok(variant_id) = variant_selector.parse::<usize>()
        {
            return Self {
                orientation_selector,
                map_index: None,
                variant_id: Some(variant_id),
            };
        }
        Self {
            orientation_selector: selector,
            map_index: None,
            variant_id: None,
        }
    }
}

fn parse_map_selector(selector: &str) -> Option<usize> {
    selector
        .strip_prefix('N')
        .and_then(|value| value.parse::<usize>().ok())
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
        #[cfg(feature = "old_cff")]
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

fn format_energy_degree_bounds(bounds: &[(usize, usize)], use_color: bool) -> String {
    if bounds.is_empty() {
        "-".to_string()
    } else {
        bounds
            .iter()
            .map(|(edge, degree)| {
                format!(
                    "{}:{}",
                    c(&edge.to_string(), Color::Green, use_color),
                    c(
                        &degree.to_string(),
                        if *degree > 1 { Color::Red } else { Color::Blue },
                        use_color
                    )
                )
            })
            .join(" ")
    }
}

struct OrientationDisplayEntry<'a> {
    id: usize,
    orientation: &'a OrientationExpression,
    base_label: String,
    map_index: Option<usize>,
    label: String,
}

fn orientation_display_entries(
    expression: &ThreeDExpression<OrientationID>,
) -> Vec<OrientationDisplayEntry<'_>> {
    let mut entries = expression
        .orientations
        .iter()
        .enumerate()
        .map(|(id, orientation)| {
            let label = orientation_label(orientation);
            let (base_label, parsed_map_index) = orientation_label_parts(&label);
            let map_index = orientation.data.numerator_map_index.or(parsed_map_index);
            let label = if let Some(map_index) = map_index {
                format!("{base_label}|N{map_index}")
            } else {
                base_label.clone()
            };
            OrientationDisplayEntry {
                id,
                orientation,
                base_label,
                map_index,
                label,
            }
        })
        .collect_vec();
    entries.sort_by(|lhs, rhs| {
        lhs.base_label
            .cmp(&rhs.base_label)
            .then(lhs.map_index.cmp(&rhs.map_index))
            .then(lhs.id.cmp(&rhs.id))
    });
    entries
}

fn orientation_label_parts(label: &str) -> (String, Option<usize>) {
    let Some((base, suffix)) = label.split_once('|') else {
        return (label.to_string(), None);
    };
    (base.to_string(), parse_map_selector(suffix))
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

fn color_orientation_label(
    label: &str,
    external_half_edge_count: usize,
    use_color: bool,
) -> String {
    if !use_color {
        return label.to_string();
    }
    let (base, suffix) = label
        .split_once('|')
        .map_or((label, None), |(base, suffix)| (base, Some(suffix)));
    let colored_base = base
        .chars()
        .enumerate()
        .map(|(position, ch)| {
            if position < external_half_edge_count {
                return Color::Fixed(8).paint(ch.to_string()).to_string();
            }
            match ch {
                '+' => Color::Green.paint(ch.to_string()).to_string(),
                '-' => Color::Red.paint(ch.to_string()).to_string(),
                '0' => Color::Yellow.paint(ch.to_string()).to_string(),
                'x' => Color::Blue.paint(ch.to_string()).to_string(),
                _ => ch.to_string(),
            }
        })
        .join("");
    if let Some(suffix) = suffix {
        format!("{}|{}", colored_base, c(suffix, Color::Purple, use_color))
    } else {
        colored_base
    }
}

fn origin_label(origin: &str, use_color: bool) -> String {
    let color = match origin {
        "mixed" => Color::Yellow,
        "cff" | "pure_cff" => Color::Green,
        "ltd" | "pure_ltd" => Color::Blue,
        _ => Color::Purple,
    };
    c(&short_origin_label(origin), color, use_color)
}

fn short_origin_label(origin: &str) -> String {
    origin
        .replace("ltd_confluent", "ltd_cflt")
        .replace("bounded_degree", "bd")
        .replace("known_factor", "known")
        .replace("e_surface", "e_surf")
        .replace("quadratic_recursive", "quad_rec")
        .replace("lower_sector", "low_sec")
        .replace("component_product", "comp_prod")
        .replace(":beta=", ":β=")
        .replace(":gamma=", ":γ=")
}

fn coefficient_label(coeff: &symbolica::atom::Atom, use_color: bool) -> String {
    c(&coeff.to_canonical_string(), Color::Yellow, use_color)
}

fn factor_list(power: usize, half_edges: &[EdgeIndex], use_color: bool) -> String {
    let mut factors = Vec::new();
    if power > 0 {
        let label = if power == 1 {
            "M".to_string()
        } else {
            format!("M^{power}")
        };
        factors.push(c(&label, Color::Purple, use_color));
    }
    let mut counts = BTreeMap::<usize, usize>::new();
    for edge in half_edges {
        *counts.entry(edge.0).or_default() += 1;
    }
    factors.extend(counts.into_iter().map(|(edge, count)| {
        let label = if count == 1 {
            edge.to_string()
        } else {
            format!("{edge}^{count}")
        };
        c(&label, Color::Blue, use_color)
    }));
    if factors.is_empty() {
        "[]".to_string()
    } else {
        format!("[{}]", factors.join(","))
    }
}

fn variant_surface_list(variant: &crate::expression::CFFVariant, use_color: bool) -> String {
    let mut all_surface_ids = BTreeSet::<HybridSurfaceID>::new();
    let denominator_powers = denominator_surface_powers(&variant.denominator);
    all_surface_ids.extend(denominator_powers.keys().copied());
    for surface_id in &variant.numerator_surfaces {
        all_surface_ids.insert(*surface_id);
    }
    if all_surface_ids.is_empty() {
        return "[]".to_string();
    }

    let mut numerator_counts = BTreeMap::<HybridSurfaceID, usize>::new();
    for surface_id in &variant.numerator_surfaces {
        *numerator_counts.entry(*surface_id).or_default() += 1;
    }

    let entries = all_surface_ids
        .into_iter()
        .map(|surface_id| {
            let denominator_labels = denominator_powers
                .get(&surface_id)
                .into_iter()
                .flat_map(|powers| powers.iter())
                .map(|power| {
                    let mut label = surface_list_id(surface_id, use_color);
                    if *power > 1 {
                        label.push_str(&format!("^{}", power));
                    }
                    label
                })
                .join(",");
            if let Some(count) = numerator_counts.get(&surface_id) {
                let mut numerator_label = surface_list_id(surface_id, use_color);
                if *count > 1 {
                    numerator_label.push_str(&format!("^{}", count));
                }
                if denominator_labels.is_empty() {
                    format!("({numerator_label})")
                } else {
                    format!("{denominator_labels} ({numerator_label})")
                }
            } else {
                denominator_labels
            }
        })
        .join(",");
    format!("[{entries}]")
}

fn denominator_surface_powers(
    tree: &crate::tree::Tree<HybridSurfaceID>,
) -> BTreeMap<HybridSurfaceID, BTreeSet<usize>> {
    let mut powers = BTreeMap::<HybridSurfaceID, BTreeSet<usize>>::new();
    if tree.get_num_nodes() == 0 {
        return powers;
    }
    let mut counts = BTreeMap::<HybridSurfaceID, usize>::new();
    collect_denominator_surface_powers(tree, crate::tree::NodeId::root(), &mut counts, &mut powers);
    powers
}

fn collect_denominator_surface_powers(
    tree: &crate::tree::Tree<HybridSurfaceID>,
    node_id: crate::tree::NodeId,
    counts: &mut BTreeMap<HybridSurfaceID, usize>,
    powers: &mut BTreeMap<HybridSurfaceID, BTreeSet<usize>>,
) {
    let node = tree.get_node(node_id);
    let counted = !matches!(node.data, HybridSurfaceID::Unit | HybridSurfaceID::Infinite);
    if counted {
        *counts.entry(node.data).or_default() += 1;
    }
    if node.children.is_empty() {
        for (surface_id, count) in counts.iter() {
            if *count > 0 {
                powers.entry(*surface_id).or_default().insert(*count);
            }
        }
    } else {
        for child in &node.children {
            collect_denominator_surface_powers(tree, *child, counts, powers);
        }
    }
    if counted && let Some(count) = counts.get_mut(&node.data) {
        *count -= 1;
        if *count == 0 {
            counts.remove(&node.data);
        }
    }
}

fn surface_list_id(surface_id: HybridSurfaceID, use_color: bool) -> String {
    match surface_id {
        HybridSurfaceID::Linear(id) => c(&id.0.to_string(), Color::Blue, use_color),
        _ => c(&surface_id_text(surface_id), Color::Blue, use_color),
    }
}

fn surface_type_label(
    kind: LinearSurfaceKind,
    origin: SurfaceOrigin,
    numerator_only: bool,
    use_color: bool,
) -> String {
    let mut label = match kind {
        LinearSurfaceKind::Esurface => "e".to_string(),
        LinearSurfaceKind::Hsurface => "h".to_string(),
    };
    if matches!(origin, SurfaceOrigin::Helper) {
        label.push_str("_sp");
    }
    if numerator_only {
        label = format!("({label})");
    }
    let color = match kind {
        LinearSurfaceKind::Esurface => Color::Green,
        LinearSurfaceKind::Hsurface => Color::Purple,
    };
    c(&label, color, use_color)
}

struct FormattedEnergyGroups {
    internal_plain: String,
    internal_colored: String,
    external_colored: String,
}

fn format_linear_energy_expr(expr: &LinearEnergyExpr, use_color: bool) -> String {
    let groups = format_linear_energy_groups(expr, use_color);
    match (
        groups.internal_colored.is_empty(),
        groups.external_colored.is_empty(),
    ) {
        (true, true) => "0".to_string(),
        (false, true) => groups.internal_colored,
        (true, false) => groups.external_colored,
        (false, false) => format!("{} {}", groups.internal_colored, groups.external_colored),
    }
}

fn format_linear_energy_groups(expr: &LinearEnergyExpr, use_color: bool) -> FormattedEnergyGroups {
    let internal = expr
        .internal_terms
        .iter()
        .sorted_by_key(|(edge_id, _)| edge_id.0)
        .map(|(edge_id, coeff)| signed_energy_term("OSE", edge_id.0, coeff, use_color))
        .collect::<Vec<_>>();
    let mut external = expr
        .external_terms
        .iter()
        .sorted_by_key(|(edge_id, _)| edge_id.0)
        .map(|(edge_id, coeff)| signed_energy_term("E", edge_id.0, coeff, use_color))
        .collect::<Vec<_>>();
    if !expr.uniform_scale_coeff.is_zero_coeff() {
        external.push(signed_named_term("M", &expr.uniform_scale_coeff, use_color));
    }
    if !expr.constant.is_zero_coeff() {
        external.insert(0, signed_named_term("1", &expr.constant, use_color));
    }

    FormattedEnergyGroups {
        internal_plain: internal.iter().map(|term| term.plain.as_str()).join(" "),
        internal_colored: internal.iter().map(|term| term.colored.as_str()).join(" "),
        external_colored: external.iter().map(|term| term.colored.as_str()).join(" "),
    }
}

struct SignedTerm {
    plain: String,
    colored: String,
}

fn signed_energy_term(
    function_name: &str,
    id: usize,
    coeff: &symbolica::atom::Atom,
    use_color: bool,
) -> SignedTerm {
    signed_named_term(&format!("{function_name}[{id}]"), coeff, use_color)
}

fn signed_named_term(name: &str, coeff: &symbolica::atom::Atom, use_color: bool) -> SignedTerm {
    let negative = coeff.is_negative_coeff();
    let sign = if negative { "-" } else { "+" };
    let magnitude = coeff.abs_coeff();
    let coefficient = if magnitude.is_one_coeff() {
        String::new()
    } else {
        format!("{}*", magnitude.to_canonical_string())
    };
    let plain = format!("{sign} {coefficient}{name}");
    let colored = format!(
        "{} {}{}",
        c(
            sign,
            if negative { Color::Red } else { Color::Green },
            use_color
        ),
        if coefficient.is_empty() {
            String::new()
        } else {
            c(&coefficient, Color::Yellow, use_color)
        },
        color_function_call(name, use_color)
    );
    SignedTerm { plain, colored }
}

fn color_function_call(name: &str, use_color: bool) -> String {
    let Some((head, rest)) = name.split_once('[') else {
        return c(name, Color::Purple, use_color);
    };
    let argument = rest.trim_end_matches(']');
    format!("{}[{}]", head, c(argument, Color::Purple, use_color))
}

fn surface_id_text(surface_id: HybridSurfaceID) -> String {
    match surface_id {
        HybridSurfaceID::Esurface(id) => format!("E{}", id.0),
        HybridSurfaceID::Hsurface(id) => format!("H{}", id.0),
        HybridSurfaceID::Linear(id) => format!("S{}", id.0),
        HybridSurfaceID::Unit => "1".to_string(),
        HybridSurfaceID::Infinite => "inf".to_string(),
    }
}

fn surface_id_label(surface_id: HybridSurfaceID, use_color: bool) -> String {
    match surface_id {
        HybridSurfaceID::Linear(id) => {
            format!(
                "{}{}",
                c("S", Color::Green, use_color),
                c(&id.0.to_string(), Color::Blue, use_color)
            )
        }
        HybridSurfaceID::Esurface(id) => {
            format!(
                "{}{}",
                c("E", Color::Green, use_color),
                c(&id.0.to_string(), Color::Blue, use_color)
            )
        }
        HybridSurfaceID::Hsurface(id) => {
            format!(
                "{}{}",
                c("H", Color::Green, use_color),
                c(&id.0.to_string(), Color::Blue, use_color)
            )
        }
        HybridSurfaceID::Unit => "1".to_string(),
        HybridSurfaceID::Infinite => "inf".to_string(),
    }
}

fn node_id_label(node_id: usize, use_color: bool) -> String {
    c(&node_id.to_string(), Color::Yellow, use_color)
}

fn tree_shape(tree: &crate::tree::Tree<HybridSurfaceID>, use_color: bool) -> String {
    tree.iter_nodes()
        .map(|node| {
            format!(
                "{}: {} → [{}]",
                node_id_label(node.node_id.0, use_color),
                surface_id_label(node.data, use_color),
                node.children
                    .iter()
                    .map(|child| node_id_label(child.0, use_color))
                    .join(" ")
            )
        })
        .join(" ")
}

fn truncate_for_table(text: &str, max_chars: usize) -> String {
    let mut visible = 0usize;
    let mut out = String::new();
    let mut chars = text.chars().peekable();
    let visible_limit = max_chars.saturating_sub(3);
    let mut truncated = false;

    while let Some(ch) = chars.next() {
        if ch == '\u{1b}' {
            out.push(ch);
            for seq_ch in chars.by_ref() {
                out.push(seq_ch);
                if seq_ch == 'm' {
                    break;
                }
            }
            continue;
        }
        if visible >= visible_limit {
            truncated = true;
            break;
        }
        out.push(ch);
        visible += 1;
    }

    if !truncated && chars.next().is_none() {
        return text.to_string();
    }
    if text.contains('\u{1b}') {
        out.push_str("\u{1b}[0m");
    }
    out.push_str("...");
    if text.contains('\u{1b}') {
        out.push_str("\u{1b}[0m");
    }
    out
}

fn orientation_matches(entry: &OrientationDisplayEntry<'_>, selector: &DetailSelector<'_>) -> bool {
    if selector.orientation_selector == entry.id.to_string() {
        return selector
            .map_index
            .is_none_or(|map_index| entry.map_index == Some(map_index));
    }
    if entry.base_label == selector.orientation_selector
        || entry.label == selector.orientation_selector
    {
        return selector
            .map_index
            .is_none_or(|map_index| entry.map_index == Some(map_index));
    }
    if let Some((base, suffix)) = selector.orientation_selector.split_once('|')
        && base == entry.base_label
        && let Some(map_index) = parse_map_selector(suffix)
    {
        return entry.map_index == Some(map_index);
    }
    if selector.orientation_selector == entry.label {
        return true;
    }
    false
}

#[cfg(test)]
mod tests {
    use crate::{
        Generate3DExpressionOptions, HybridSurfaceID, LinearSurfaceID, RepresentationMode,
        expression::CFFVariant,
        generation::generate_3d_expression_from_parsed,
        graph_io::graph_info,
        surface::rational_coeff_one,
        tree::{NodeId, Tree},
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
            NumeratorDisplay::default(),
            NumeratorSamplingScaleMode::None,
            &DisplayOptions {
                use_color: false,
                details_for_orientation: None,
            },
        );

        assert!(rendered.contains("CFF structure"));
        assert!(rendered.contains("Surfaces"));
        assert!(rendered.contains("Orientations"));
        assert!(rendered.contains("0: 1 → ["));
        assert!(!rendered.contains("1 → [S"));

        let rendered = render_expression_summary(
            &expression,
            RepresentationMode::Cff,
            &graph_info(&parsed),
            &[],
            NumeratorDisplay::default(),
            NumeratorSamplingScaleMode::None,
            &DisplayOptions {
                use_color: false,
                details_for_orientation: Some("0".to_string()),
            },
        );

        assert!(!rendered.contains("CFF structure"));
        assert!(rendered.contains("loop q0 map"));
        assert!(rendered.contains("edge q0 map"));
    }

    #[test]
    fn display_accepts_map_and_variant_detail_selectors() {
        let parsed = crate::graph_io::test_graphs::box_pow3_graph();
        let expression = generate_3d_expression_from_parsed(
            &parsed,
            &Generate3DExpressionOptions {
                representation: RepresentationMode::Cff,
                energy_degree_bounds: vec![(0, 1), (1, 1), (3, 4)],
                ..Default::default()
            },
        )
        .unwrap();
        let label = orientation_label(expression.orientations.iter().next().unwrap());
        let (base, map_index) = orientation_label_parts(&label);
        let selector = format!("{}|N{}|0", base, map_index.unwrap_or_default());

        let rendered = render_expression_summary(
            &expression,
            RepresentationMode::Cff,
            &graph_info(&parsed),
            &[],
            NumeratorDisplay::default(),
            NumeratorSamplingScaleMode::None,
            &DisplayOptions {
                use_color: false,
                details_for_orientation: Some(selector),
            },
        );

        assert!(rendered.contains("Orientation"));
        assert!(rendered.contains("variant 0"));
        assert!(rendered.contains("children"));
        assert_eq!(rendered.matches("Orientation ").count(), 1);
    }

    #[test]
    fn surface_list_reports_all_branch_powers() {
        let surface = HybridSurfaceID::Linear(LinearSurfaceID(5));
        let mut denominator = Tree::from_root(HybridSurfaceID::Unit);
        denominator.insert_node(NodeId::root(), surface);
        denominator.insert_node(NodeId::root(), surface);
        denominator.insert_node(NodeId(2), surface);
        let variant = CFFVariant {
            origin: Some("test".to_string()),
            prefactor: rational_coeff_one(),
            half_edges: Vec::new(),
            denominator_edges: Vec::new(),
            uniform_scale_power: 0,
            numerator_surfaces: Vec::new(),
            denominator,
        };

        assert_eq!(variant_surface_list(&variant, false), "[5,5^2]");
    }
}
