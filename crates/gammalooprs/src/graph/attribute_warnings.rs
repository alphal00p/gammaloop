use ahash::AHashSet;
use linnet::parser::{DotGraph, DotVertexData};
use tracing::warn;

struct AttributeSpec {
    entity: &'static str,
    known: &'static [&'static str],
    ignored: &'static [&'static str],
    aliases: &'static [(&'static str, &'static str)],
}

const GRAPH_SPEC: AttributeSpec = AttributeSpec {
    entity: "graph",
    known: &[
        "group_id",
        "is_group_master",
        "num",
        "overall_factor",
        "params",
        "projector",
    ],
    ignored: &[
        "bgcolor",
        "color",
        "label",
        "layout",
        "nodesep",
        "overlap",
        "overall_factor_evaluated",
        "pad",
        "rankdir",
        "ranksep",
        "ratio",
        "size",
        "splines",
        "start",
    ],
    aliases: &[
        ("group_master", "is_group_master"),
        ("numprojector", "projector"),
    ],
};

const EDGE_SPEC: AttributeSpec = AttributeSpec {
    entity: "edge",
    known: &[
        "dod",
        "is_cut",
        "is_dummy",
        "lmb_id",
        "mass",
        "momtrop_edge_power",
        "name",
        "num",
        "particle",
        "pdg",
        "vakint_edge_power",
    ],
    ignored: &[
        "color",
        "constraint",
        "decorate",
        "dir",
        "flow",
        "fontcolor",
        "fontsize",
        "id",
        "label",
        "lmb_rep",
        "penwidth",
        "pin",
        "sink",
        "source",
        "style",
    ],
    aliases: &[
        ("dummy", "is_dummy"),
        ("iscut", "is_cut"),
        ("lmb_index", "lmb_id"),
        ("momtrop_power", "momtrop_edge_power"),
        ("numerator", "num"),
        ("pdg_code", "pdg"),
        ("vakint_power", "vakint_edge_power"),
    ],
};

const NODE_SPEC: AttributeSpec = AttributeSpec {
    entity: "node",
    known: &["dod", "int_id", "name", "num"],
    ignored: &[
        "color",
        "fillcolor",
        "fontcolor",
        "fontname",
        "fontsize",
        "height",
        "label",
        "penwidth",
        "pin",
        "pos",
        "shape",
        "style",
        "width",
        "xlabel",
    ],
    aliases: &[
        ("inter_id", "int_id"),
        ("numerator", "num"),
        ("vertex_rule", "int_id"),
    ],
};

pub(crate) fn warn_about_unknown_attributes(graph: &DotGraph) {
    let mut emitted = AHashSet::new();

    warn_attributes(
        graph.global_data.statements.keys().map(String::as_str),
        &GRAPH_SPEC,
        &graph.global_data.name,
        "graph",
        None,
        &mut emitted,
    );
    warn_attributes(
        graph.global_data.edge_statements.keys().map(String::as_str),
        &EDGE_SPEC,
        &graph.global_data.name,
        "default edge",
        None,
        &mut emitted,
    );
    warn_attributes(
        graph.global_data.node_statements.keys().map(String::as_str),
        &NODE_SPEC,
        &graph.global_data.name,
        "default node",
        None,
        &mut emitted,
    );

    for (_, _, edge) in graph.graph.iter_edges() {
        warn_attributes(
            edge.data.local_statements.keys().map(String::as_str),
            &EDGE_SPEC,
            &graph.global_data.name,
            "edge",
            edge_identifier(edge.data),
            &mut emitted,
        );
    }

    for (node_id, _, vertex) in graph.graph.iter_nodes() {
        warn_attributes(
            vertex
                .statements
                .keys()
                .filter(|key| !graph.global_data.node_statements.contains_key(*key))
                .map(String::as_str),
            &NODE_SPEC,
            &graph.global_data.name,
            "node",
            Some(node_identifier(node_id, vertex)),
            &mut emitted,
        );
    }
}

fn warn_attributes<'a>(
    attrs: impl IntoIterator<Item = &'a str>,
    spec: &AttributeSpec,
    graph_name: &str,
    owner: &str,
    identifier: Option<String>,
    emitted: &mut AHashSet<String>,
) {
    for attr in attrs {
        if spec
            .known
            .iter()
            .chain(spec.ignored.iter())
            .any(|known| *known == attr)
        {
            continue;
        }

        let mut message = format!("Unknown {owner} attribute '{attr}' on graph '{graph_name}'");
        if let Some(identifier) = &identifier {
            message.push_str(&format!(", {} '{identifier}'", spec.entity));
        }
        if let Some(suggestion) = best_match(attr, spec) {
            message.push_str(&format!("; did you mean '{suggestion}'?"));
        } else {
            message.push('.');
        }

        if emitted.insert(message.clone()) {
            warn!("{message}");
        }
    }
}

fn edge_identifier(edge: &linnet::parser::DotEdgeData) -> Option<String> {
    edge.statements
        .get("name")
        .cloned()
        .or_else(|| edge.edge_id.map(|edge_id| edge_id.0.to_string()))
}

fn node_identifier(node_id: linnet::half_edge::NodeIndex, vertex: &DotVertexData) -> String {
    vertex
        .name()
        .map(ToOwned::to_owned)
        .or_else(|| vertex.index.map(|index| index.0.to_string()))
        .unwrap_or_else(|| node_id.0.to_string())
}

fn best_match<'a>(attr: &str, spec: &'a AttributeSpec) -> Option<&'a str> {
    if let Some((_, target)) = spec.aliases.iter().find(|(alias, _)| *alias == attr) {
        return Some(*target);
    }

    let normalized_attr = normalize(attr);
    spec.known
        .iter()
        .copied()
        .filter_map(|candidate| {
            let normalized_candidate = normalize(candidate);
            let distance = edit_distance(&normalized_attr, &normalized_candidate);
            let threshold = normalized_attr
                .len()
                .max(normalized_candidate.len())
                .clamp(1, 12)
                / 3
                + 1;
            let looks_related = normalized_candidate.contains(&normalized_attr)
                || normalized_attr.contains(&normalized_candidate)
                || distance <= threshold;
            looks_related.then_some((distance, candidate))
        })
        .min_by_key(|(distance, candidate)| (*distance, candidate.len()))
        .map(|(_, candidate)| candidate)
}

fn normalize(value: &str) -> String {
    value
        .chars()
        .filter(|ch| ch.is_ascii_alphanumeric())
        .flat_map(char::to_lowercase)
        .collect()
}

fn edit_distance(lhs: &str, rhs: &str) -> usize {
    if lhs.is_empty() {
        return rhs.len();
    }
    if rhs.is_empty() {
        return lhs.len();
    }

    let rhs_chars = rhs.chars().collect::<Vec<_>>();
    let mut prev = (0..=rhs_chars.len()).collect::<Vec<_>>();
    let mut curr = vec![0; rhs_chars.len() + 1];

    for (i, lhs_char) in lhs.chars().enumerate() {
        curr[0] = i + 1;
        for (j, rhs_char) in rhs_chars.iter().enumerate() {
            let substitution_cost = usize::from(lhs_char != *rhs_char);
            curr[j + 1] = (prev[j + 1] + 1)
                .min(curr[j] + 1)
                .min(prev[j] + substitution_cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    prev[rhs_chars.len()]
}

#[cfg(test)]
mod tests {
    use super::{EDGE_SPEC, GRAPH_SPEC, NODE_SPEC, best_match};

    #[test]
    fn matches_common_graph_typos() {
        assert_eq!(best_match("numprojector", &GRAPH_SPEC), Some("projector"));
        assert_eq!(
            best_match("group_master", &GRAPH_SPEC),
            Some("is_group_master")
        );
    }

    #[test]
    fn matches_common_edge_typos() {
        assert_eq!(best_match("lmb_index", &EDGE_SPEC), Some("lmb_id"));
        assert_eq!(
            best_match("vakint_power", &EDGE_SPEC),
            Some("vakint_edge_power")
        );
    }

    #[test]
    fn matches_common_node_typos() {
        assert_eq!(best_match("inter_id", &NODE_SPEC), Some("int_id"));
        assert_eq!(best_match("numerator", &NODE_SPEC), Some("num"));
    }
}
