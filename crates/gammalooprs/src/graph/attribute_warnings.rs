//! Warn about DOT attributes that neither FeynKit nor GammaLoop interprets.

use ahash::AHashSet;
use feynkit_graph::{DotDiagram, VertexId};
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
        "final_state",
        "group_id",
        "is_group_master",
        "num",
        "overall_factor",
        "pole_part",
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
        "pos",
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

/// Warn once per unknown attribute of a compact DOT graph, naming its first owner by
/// its DOT name. FeynKit's annotated export, which carries a model fingerprint, is not
/// checked.
pub(crate) fn warn_about_unknown_attributes(dot: &DotDiagram) {
    for message in unknown_attribute_warnings(dot) {
        warn!("{message}");
    }
}

fn unknown_attribute_warnings(dot: &DotDiagram) -> Vec<String> {
    let DotDiagram {
        diagram,
        attributes,
        ..
    } = dot;
    if attributes.graph.contains_key("model_fingerprint") {
        return Vec::new();
    }
    let vertex_name = |vertex: Option<VertexId>| {
        vertex
            .and_then(|vertex| diagram.vertices().find(|(id, _)| *id == vertex))
            .map_or("ext", |(_, vertex)| vertex.name.as_str())
    };
    // Default statements are merged into every owner: report them on the first.
    let owners =
        std::iter::once((&GRAPH_SPEC, None, &attributes.graph))
            .chain(
                diagram
                    .vertices()
                    .zip(&attributes.vertices)
                    .map(|((_, vertex), attributes)| {
                        (&NODE_SPEC, Some(vertex.name.clone()), attributes)
                    }),
            )
            .chain(diagram.edges().zip(&attributes.edges).map(
                |((_, endpoints, _), attributes)| {
                    // An edge without a name is identified by its endpoints.
                    let name = attributes.get("name").cloned().unwrap_or_else(|| {
                        format!(
                            "{} -> {}",
                            vertex_name(endpoints.source),
                            vertex_name(endpoints.target)
                        )
                    });
                    (&EDGE_SPEC, Some(name), attributes)
                },
            ));
    let mut emitted = AHashSet::new();
    let mut warnings = Vec::new();
    for (spec, owner, attributes) in owners {
        for attr in attributes.keys() {
            if spec
                .known
                .iter()
                .chain(spec.ignored)
                .any(|known| known == attr)
                || !emitted.insert((spec.entity, attr))
            {
                continue;
            }
            let mut message = format!(
                "Unknown {} attribute '{attr}' on graph '{}'",
                spec.entity,
                diagram.name().unwrap_or("<unnamed>")
            );
            if let Some(owner) = &owner {
                message.push_str(&format!(", {} '{owner}'", spec.entity));
            }
            if let Some(suggestion) = best_match(attr, spec) {
                message.push_str(&format!("; did you mean '{suggestion}'?"));
            } else {
                message.push('.');
            }
            warnings.push(message);
        }
    }
    warnings
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
    use feynkit_graph::DotDiagram;

    use super::{EDGE_SPEC, GRAPH_SPEC, NODE_SPEC, best_match, unknown_attribute_warnings};
    use crate::{initialisation::test_initialise, utils::load_generic_model};

    #[test]
    fn warns_once_per_unknown_attribute_naming_its_dot_owner() {
        test_initialise().unwrap();
        let model = load_generic_model("scalars");
        let dot = |statements: &str| {
            DotDiagram::from_dot(
                model.clone(),
                &format!(
                    r#"digraph forward {{
                        {statements}
                        ext [style=invis];
                        ext -> a [particle="scalar_0", is_cut=0, colour=red];
                        a -> b [particle="scalar_0", lmb_index=0];
                        a -> b [particle="scalar_0", name=lower];
                        b -> ext [particle="scalar_0", is_cut=0];
                    }}"#
                ),
            )
            .unwrap()
        };
        let warnings = unknown_attribute_warnings(&dot(
            r#"node [inter_id=V_3]; edge [weight=2]; graph [numprojector=1];"#,
        ));
        assert_eq!(
            warnings,
            [
                "Unknown graph attribute 'numprojector' on graph 'forward'; did you mean 'projector'?",
                "Unknown node attribute 'inter_id' on graph 'forward', node 'a'; did you mean 'int_id'?",
                "Unknown edge attribute 'colour' on graph 'forward', edge 'b -> a'.",
                "Unknown edge attribute 'weight' on graph 'forward', edge 'b -> a'.",
                "Unknown edge attribute 'lmb_index' on graph 'forward', edge 'a -> b'; did you mean 'lmb_id'?",
            ]
        );
        // FeynKit's annotated export is not checked.
        let annotated =
            DotDiagram::from_dot(model.clone(), &dot("").diagram.to_dot().unwrap()).unwrap();
        assert!(unknown_attribute_warnings(&annotated).is_empty());
    }

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
