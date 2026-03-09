use std::{borrow::Borrow, collections::BTreeMap, fmt::Display, str::FromStr};

use itertools::Either;

use crate::half_edge::{
    builder::HedgeData,
    involution::{EdgeIndex, Flow, Orientation},
};

use super::{strip_quotes, subgraph_free::Edge, DotHedgeData, GlobalData, NodeIdOrDangling};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DotEdgeData {
    pub statements: BTreeMap<String, String>,
    pub local_statements: BTreeMap<String, String>,
    pub edge_id: Option<EdgeIndex>,
}

impl FromIterator<(String, String)> for DotEdgeData {
    fn from_iter<T: IntoIterator<Item = (String, String)>>(iter: T) -> Self {
        let mut edge_id = None;
        let statements: BTreeMap<String, String> = iter
            .into_iter()
            .filter_map(|(k, v)| match k.as_str() {
                "id" => {
                    edge_id = Some(EdgeIndex::from(strip_quotes(&v).parse::<usize>().unwrap()));
                    None
                }
                _ => Some((k, strip_quotes(&v).to_string())),
            })
            .collect();

        DotEdgeData {
            local_statements: BTreeMap::new(),
            statements,

            edge_id,
        }
    }
}

impl Display for DotEdgeData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut first = true;
        for (key, value) in &self.statements {
            if !first {
                write!(f, " ")?;
            }
            write!(f, "{key}=\"{value}\"")?;
            first = false;
        }

        // if let Some(id) = self.edge_id {
        //     write!(f, " id={}", id.0)?;
        // }

        Ok(())
    }
}

impl DotEdgeData {
    pub fn remove_common(&mut self, common: &GlobalData) {
        self.statements.retain(|k, v| {
            if let Some(common_value) = common.edge_statements.get(k) {
                v != common_value
            } else {
                true
            }
        });
    }
    pub fn empty() -> Self {
        DotEdgeData {
            statements: BTreeMap::new(),
            local_statements: BTreeMap::new(),
            edge_id: None,
        }
    }

    pub fn extend(&mut self, other: Self) {
        self.statements.extend(other.statements);
    }

    pub fn format(&self, template: impl AsRef<str>) -> String {
        let mut result = template.as_ref().to_owned();

        // Find all occurrences of {key} in the template
        for (key, value) in &self.statements {
            let placeholder = format!("{{{key}}} ");
            result = result.replace(&placeholder, value);
        }

        // Remove any remaining {whatever} patterns
        while let Some(start) = result.find('{') {
            if let Some(end) = result[start..].find('}') {
                result.replace_range(start..=start + end, "");
            } else {
                break;
            }
        }

        result
    }

    pub fn add_statement(&mut self, key: impl ToString, value: impl ToString) {
        self.statements.insert(key.to_string(), value.to_string());
    }

    pub fn get<Q: Ord + ?Sized, F: FromStr>(&self, key: &Q) -> Option<Result<F, F::Err>>
    where
        String: Borrow<Q>,
    {
        self.statements.get(key).map(|s| s.parse())
    }
}

impl Default for DotEdgeData {
    fn default() -> Self {
        Self::empty()
    }
}

impl DotEdgeData {
    pub fn from_parser(
        edge: Edge,
        map: &BTreeMap<String, NodeIdOrDangling>,
        is_digraph: impl Into<Orientation>,
        global_data: &GlobalData,
    ) -> (
        Self,
        Orientation,
        HedgeData<DotHedgeData>,
        Either<HedgeData<DotHedgeData>, Flow>,
    ) {
        let mut orientation = is_digraph.into();
        let mut source_data = DotHedgeData::from(edge.source_port());
        let mut sink_data = DotHedgeData::from(edge.sink_port());
        let local_statements = edge.attr.clone().into_iter().collect();
        let mut statements = global_data.edge_statements.clone();
        statements.extend(edge.attr.into_iter().filter_map(|(key, value)| {
            let stripped_value = strip_quotes(&value);
            match key.as_str() {
                "dir" => match stripped_value {
                    "forward" => orientation = Orientation::Default,
                    "back" => orientation = Orientation::Reversed,
                    "none" => orientation = Orientation::Undirected,
                    _ => panic!("Invalid edge direction"),
                },
                "source" => {
                    source_data.statement = Some(stripped_value.to_string());
                }
                "sink" => {
                    sink_data.statement = Some(stripped_value.to_string());
                }
                _ => {
                    return Some((key, stripped_value.to_string()));
                }
            }
            None
        }));

        let source = map[&edge.from.id].clone();

        let target = map[&edge.to.id].clone();

        let (edge, source, target) = match (source, target) {
            (NodeIdOrDangling::Id(source), NodeIdOrDangling::Id(target)) => {
                //Full edge
                let mut dot_edge: DotEdgeData = statements.into_iter().collect();
                dot_edge.local_statements = local_statements;
                (
                    dot_edge,
                    source.add_data(source_data),
                    Either::Left(target.add_data(sink_data)),
                )
            }
            (NodeIdOrDangling::Id(source), NodeIdOrDangling::Dangling { statements: states }) => {
                statements.extend(
                    states
                        .into_iter()
                        .filter(|(a, _)| !(a.as_str() == "shape" || a.as_str() == "label")),
                );
                let mut dot_edge: DotEdgeData = statements.into_iter().collect();
                dot_edge.local_statements = local_statements;
                if !sink_data.is_none() {
                    panic!("Sink edge cannot have data:{sink_data}");
                }
                (
                    dot_edge,
                    source.add_data(source_data),
                    Either::Right(Flow::Source),
                )
            }
            (NodeIdOrDangling::Dangling { statements: states }, NodeIdOrDangling::Id(sink)) => {
                statements.extend(
                    states
                        .into_iter()
                        .filter(|(a, _)| !(a.as_str() == "shape" || a.as_str() == "label")),
                );
                let mut dot_edge: DotEdgeData = statements.into_iter().collect();
                dot_edge.local_statements = local_statements;
                if !source_data.is_none() {
                    panic!("Source edge cannot have data:{source_data}");
                }

                (
                    dot_edge,
                    sink.add_data(sink_data),
                    Either::Right(Flow::Sink),
                )
            }
            _ => panic!("Cannot connect an edge to two external nodes"),
        };

        (edge, orientation, source, target)
    }
}
