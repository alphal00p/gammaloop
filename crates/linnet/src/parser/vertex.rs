use std::{borrow::Borrow, collections::BTreeMap, fmt::Display, str::FromStr};

use itertools::Either;

use crate::half_edge::NodeIndex;

use super::{dot_id, escape_dot_string, strip_quotes, DotParseError, ExplicitIdKind, GlobalData};

#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
#[cfg_attr(feature = "rkyv", archive(check_bytes))]
pub struct DotVertexData {
    pub name: Option<String>,
    pub index: Option<NodeIndex>,
    pub payload: Option<Vec<u8>>,
    pub statements: BTreeMap<String, String>,
}
impl DotVertexData {
    pub fn remove_common(&mut self, common: &GlobalData) {
        self.statements.retain(|k, v| {
            if let Some(common_value) = common.node_statements.get(k) {
                v != common_value
            } else {
                true
            }
        });
    }
    pub fn name(&self) -> Option<&str> {
        if let Some(d) = self.statements.get("name") {
            Some(d.as_str())
        } else {
            self.name.as_deref()
        }
    }

    pub fn format(&self, template: impl AsRef<str>) -> String {
        let mut result = template.as_ref().to_owned();

        // Find all occurrences of {key} in the template
        for (key, value) in &self.statements {
            let placeholder = format!("{{{key}}}");
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

    pub fn empty() -> Self {
        DotVertexData {
            index: None,
            name: None,
            payload: None,
            statements: BTreeMap::new(),
        }
    }

    pub fn extend(&mut self, other: Self) {
        self.statements.extend(other.statements);
    }

    pub fn get<Q: Ord + ?Sized, F: FromStr>(&self, key: &Q) -> Option<Result<F, F::Err>>
    where
        String: Borrow<Q>,
    {
        self.statements.get(key).map(|s| s.parse())
    }

    pub fn add_statement(&mut self, key: impl ToString, value: impl ToString) {
        self.statements.insert(key.to_string(), value.to_string());
    }

    pub fn try_id_from_name(&mut self) {
        if self.index.is_none() {
            if let Some(name) = self.name.clone() {
                if let Ok(id) = name.parse::<usize>() {
                    self.index = Some(NodeIndex(id));
                    self.name = None;
                }
            }
        }
    }

    pub fn from_parser(
        value: dot_parser::canonical::Node<(String, String)>,
        global: &GlobalData,
    ) -> Result<Either<Self, BTreeMap<String, String>>, DotParseError> {
        let mut is_dangling = false;
        let mut index = global
            .node_statements
            .get("id")
            .map(|value| {
                value.parse::<usize>().map(NodeIndex).map_err(|source| {
                    DotParseError::InvalidExplicitId {
                        kind: ExplicitIdKind::Node,
                        value: value.clone(),
                        source,
                    }
                })
            })
            .transpose()?;
        let mut node_statements = BTreeMap::new();
        for (key, value) in value.attr {
            let key = strip_quotes(&key);
            match key.as_str() {
                "style" => {
                    if value.as_str() == "invis" {
                        is_dangling = true
                    }
                }
                "id" => {
                    let value = strip_quotes(&value);
                    index = Some(NodeIndex(value.parse::<usize>().map_err(|source| {
                        DotParseError::InvalidExplicitId {
                            kind: ExplicitIdKind::Node,
                            value: value.to_string(),
                            source,
                        }
                    })?));
                }
                _ => {
                    node_statements.insert(key, strip_quotes(&value));
                }
            }
        }

        if is_dangling {
            Ok(Either::Right(node_statements))
        } else {
            let mut statements = global.node_statements.clone();
            statements.remove("id");
            statements.extend(node_statements);

            let mut node = DotVertexData {
                name: Some(strip_quotes(&value.id)),
                index,
                payload: None,
                statements,
            };

            node.try_id_from_name();
            Ok(Either::Left(node))
        }
    }
}

impl Display for DotVertexData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut first = true;
        // if let Some(id) = self.name() {
        //     first = false;
        //     write!(f, "name={id}")?;
        // }
        for (key, value) in &self.statements {
            if key.starts_with("__linnest-") {
                continue;
            }
            if !first {
                write!(f, " ")?;
            }
            write!(f, "{}=\"{}\"", dot_id(key), escape_dot_string(value))?;
            first = false;
        }
        Ok(())
    }
}
