use std::{borrow::Borrow, collections::BTreeMap, fmt::Display, str::FromStr};

use itertools::Either;

use crate::half_edge::NodeIndex;

use super::{strip_quotes, GlobalData};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DotVertexData {
    pub name: Option<String>,
    pub index: Option<NodeIndex>,
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
    ) -> Either<Self, BTreeMap<String, String>> {
        let mut is_dangling = false;
        let mut index = None;
        let node_statements: BTreeMap<String, String> = value
            .attr
            .into_iter()
            .filter_map(|(key, value)| {
                match key.as_str() {
                    "style" => {
                        if value.as_str() == "invis" {
                            is_dangling = true
                        }
                    }
                    "id" => index = Some(NodeIndex(strip_quotes(&value).parse::<usize>().unwrap())),
                    _ => return Some((key, strip_quotes(&value).to_string())),
                }
                None
            })
            .collect();

        if is_dangling {
            Either::Right(node_statements)
        } else {
            let mut statements = global.node_statements.clone();
            statements.extend(node_statements);

            let mut node = DotVertexData {
                name: Some(value.id),
                index,
                statements,
            };

            node.try_id_from_name();
            Either::Left(node)
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
            // if key == "name" {
            //     continue;
            // }
            if !first {
                write!(f, " ")?;
            }
            write!(f, "{key}=\"{value}\"")?;
            first = false;
        }
        Ok(())
    }
}
