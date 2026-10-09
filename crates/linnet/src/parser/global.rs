use std::{collections::BTreeMap, fmt::Display};

use figment::Figment;

use super::{dot_id, dot_value, syntax::Attr};

#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
#[cfg_attr(feature = "rkyv", archive(check_bytes))]
pub struct GlobalData {
    pub name: String,
    pub payload: Option<Vec<u8>>,
    pub statements: BTreeMap<String, String>,
    pub edge_statements: BTreeMap<String, String>,
    pub node_statements: BTreeMap<String, String>,
}

impl GlobalData {
    pub fn add_name(&mut self, name: String) {
        self.name = name;
    }

    /// Graph metadata from the flattened attribute statements of a DOT graph.
    /// Later statements win; a `name` statement becomes the graph name.
    pub(crate) fn from_statements<'a>(
        statements: impl IntoIterator<Item = Attr<'a>>,
        node_defaults: impl IntoIterator<Item = Attr<'a>>,
        edge_defaults: impl IntoIterator<Item = Attr<'a>>,
    ) -> Self {
        fn owned<'a>(attrs: impl IntoIterator<Item = Attr<'a>>) -> BTreeMap<String, String> {
            (attrs.into_iter())
                .map(|(key, value)| (key.into_owned(), value.into_owned()))
                .collect()
        }
        let mut statements = owned(statements);
        GlobalData {
            name: statements.remove("name").unwrap_or_default(),
            payload: None,
            statements,
            edge_statements: owned(edge_defaults),
            node_statements: owned(node_defaults),
        }
    }

    #[cfg(feature = "serde")]
    pub fn with_figment(mut self, figment: Figment) -> Self {
        self.extract_figment_data(&figment);
        self
    }

    #[cfg(feature = "serde")]
    pub fn set_figment(&mut self, figment: Figment) {
        self.extract_figment_data(&figment);
    }

    #[cfg(not(feature = "serde"))]
    pub fn set_figment(&mut self, _figment: Figment) {
        // Incremental CLI always passes a Figment instance, but without the
        // serde feature we cannot deserialize it, so we silently ignore it.
    }

    #[cfg(feature = "serde")]
    fn extract_figment_data(&mut self, figment: &Figment) {
        // Extract graph.* keys using focus
        if let Ok(graph_data) = figment.focus("graph").extract::<BTreeMap<String, String>>() {
            self.statements.extend(graph_data);
        }

        // Extract edge.* keys using focus
        if let Ok(edge_data) = figment.focus("edge").extract::<BTreeMap<String, String>>() {
            self.edge_statements.extend(edge_data);
        }

        // Extract node.* keys using focus
        if let Ok(node_data) = figment.focus("node").extract::<BTreeMap<String, String>>() {
            self.node_statements.extend(node_data);
        }

        // Extract top-level keys (no nested structure). Runtime parameters
        // override existing global statements when callers explicitly merge a
        // Figment into graph metadata.
        if let Ok(all_data) = figment.extract::<BTreeMap<String, figment::value::Value>>() {
            for (key, value) in all_data {
                if !key.contains('.') {
                    if let Some(value_str) = value.into_string() {
                        self.statements.insert(key, value_str);
                    }
                }
            }
        }
    }
}

impl From<()> for GlobalData {
    fn from(_: ()) -> Self {
        GlobalData {
            name: String::new(),
            payload: None,
            statements: BTreeMap::new(),
            edge_statements: BTreeMap::new(),
            node_statements: BTreeMap::new(),
        }
    }
}

impl Display for GlobalData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if !self.statements.is_empty() {
            for (key, value) in &self.statements {
                if let Some(indent) = f.width() {
                    writeln!(
                        f,
                        "{}{} = {};",
                        vec![" "; indent].join(""),
                        dot_id(key),
                        dot_value(value)
                    )?;
                } else {
                    write!(f, "\n{} = {};", dot_id(key), dot_value(value))?;
                }
            }
        }

        if !self.edge_statements.is_empty() {
            if let Some(indent) = f.width() {
                writeln!(f, "{}edge\t [", vec![" "; indent].join(""))?;
            } else {
                write!(f, "edge\t [")?;
            }

            let mut first = true;
            for (key, value) in &self.edge_statements {
                if !first {
                    write!(f, ", ")?;
                }
                write!(f, "{} = {}", dot_id(key), dot_value(value))?;
                first = false;
            }
            writeln!(f, "]")?;
        }

        if !self.node_statements.is_empty() {
            if let Some(indent) = f.width() {
                writeln!(f, "{}node\t [", vec![" "; indent].join(""))?;
            } else {
                write!(f, "node\t [")?;
            }

            let mut first = true;
            for (key, value) in &self.node_statements {
                if !first {
                    write!(f, ", ")?;
                }
                write!(f, "{} = {}", dot_id(key), dot_value(value))?;
                first = false;
            }
            writeln!(f, "]")?;
        }

        Ok(())
    }
}
