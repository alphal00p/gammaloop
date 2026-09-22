use std::{borrow::Borrow, collections::BTreeMap, fmt::Display, str::FromStr};

use itertools::Either;

use crate::half_edge::{
    builder::HedgeData,
    involution::{EdgeIndex, Flow, Orientation},
};

use super::{
    dot_id, dot_value, strip_quotes, subgraph_free::Edge, DotHedgeData, DotParseError,
    ExplicitIdKind, GlobalData, NodeIdOrDangling,
};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
#[cfg_attr(feature = "rkyv", archive(check_bytes))]
pub struct DotEdgeData {
    pub payload: Option<Vec<u8>>,
    pub statements: BTreeMap<String, String>,
    pub local_statements: BTreeMap<String, String>,
    pub edge_id: Option<EdgeIndex>,
}

impl DotEdgeData {
    fn from_statements(
        iter: impl IntoIterator<Item = (String, String)>,
    ) -> Result<Self, DotParseError> {
        let mut edge_id = None;
        let mut statements = BTreeMap::new();
        for (key, value) in iter {
            match key.as_str() {
                "id" => {
                    edge_id = Some(EdgeIndex::from(value.parse::<usize>().map_err(
                        |source| DotParseError::InvalidExplicitId {
                            kind: ExplicitIdKind::Edge,
                            value: value.to_string(),
                            source,
                        },
                    )?));
                }
                _ => {
                    statements.insert(key, value);
                }
            }
        }

        Ok(DotEdgeData {
            payload: None,
            local_statements: BTreeMap::new(),
            statements,
            edge_id,
        })
    }
}

impl Display for DotEdgeData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut first = true;
        for (key, value) in &self.statements {
            if key.starts_with("__linnest-") {
                continue;
            }
            if !first {
                write!(f, " ")?;
            }
            write!(f, "{}={}", dot_id(key), dot_value(value))?;
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
            payload: None,
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
    #[allow(clippy::type_complexity)]
    pub fn from_parser(
        edge: Edge,
        map: &BTreeMap<String, NodeIdOrDangling>,
        is_digraph: impl Into<Orientation>,
        global_data: &GlobalData,
    ) -> Result<
        (
            Self,
            Orientation,
            HedgeData<DotHedgeData>,
            Either<HedgeData<DotHedgeData>, Flow>,
        ),
        DotParseError,
    > {
        let mut orientation = is_digraph.into();
        let mut source_data = DotHedgeData::from_parser(edge.source_port())?;
        let mut sink_data = DotHedgeData::from_parser(edge.sink_port())?;
        let local_statements = edge
            .attr
            .iter()
            .map(|(key, value)| (strip_quotes(key), strip_quotes(value)))
            .collect();
        let mut statements = global_data.edge_statements.clone();
        if let Some(value) = statements.remove("dir") {
            orientation = match value.as_str() {
                "forward" => Orientation::Default,
                "back" => Orientation::Reversed,
                "none" => Orientation::Undirected,
                _ => return Err(DotParseError::InvalidEdgeDirection { value }),
            };
        }
        for (key, value) in edge.attr {
            let key = strip_quotes(&key);
            let stripped_value = strip_quotes(&value);
            match key.as_str() {
                "dir" => match stripped_value.as_str() {
                    "forward" => orientation = Orientation::Default,
                    "back" => orientation = Orientation::Reversed,
                    "none" => orientation = Orientation::Undirected,
                    _ => {
                        return Err(DotParseError::InvalidEdgeDirection {
                            value: stripped_value,
                        });
                    }
                },
                "source" => {
                    source_data.statement = Some(stripped_value);
                }
                "sink" => {
                    sink_data.statement = Some(stripped_value);
                }
                _ => {
                    statements.insert(key, stripped_value);
                }
            }
        }

        let source = map[&edge.from.id].clone();

        let target = map[&edge.to.id].clone();

        let (edge, source, target) = match (source, target) {
            (NodeIdOrDangling::Id(source), NodeIdOrDangling::Id(target)) => {
                //Full edge
                let mut dot_edge = DotEdgeData::from_statements(statements)?;
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
                let mut dot_edge = DotEdgeData::from_statements(statements)?;
                dot_edge.local_statements = local_statements;
                if !sink_data.is_none() {
                    return Err(DotParseError::ExternalSinkEndpointData);
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
                let mut dot_edge = DotEdgeData::from_statements(statements)?;
                dot_edge.local_statements = local_statements;
                if !source_data.is_none() {
                    return Err(DotParseError::ExternalSourceEndpointData);
                }

                (
                    dot_edge,
                    sink.add_data(sink_data),
                    Either::Right(Flow::Sink),
                )
            }
            (NodeIdOrDangling::Dangling { .. }, NodeIdOrDangling::Dangling { .. }) => {
                return Err(DotParseError::EdgeBetweenExternalNodes)
            }
        };

        Ok((edge, orientation, source, target))
    }
}
