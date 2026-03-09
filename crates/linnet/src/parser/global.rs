use std::{collections::BTreeMap, fmt::Display};

use dot_parser::{ast::AttrStmt, canonical::IDEq};
use figment::Figment;

use super::strip_quotes;

#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GlobalData {
    pub name: String,
    pub statements: BTreeMap<String, String>,
    pub edge_statements: BTreeMap<String, String>,
    pub node_statements: BTreeMap<String, String>,
}

impl GlobalData {
    pub fn add_name(&mut self, name: String) {
        self.name = name;
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

        // Extract top-level keys (no nested structure)
        // Only add CLI parameters if they don't already exist in DOT statements
        // This ensures DOT parameters take precedence over CLI parameters
        if let Ok(all_data) = figment.extract::<BTreeMap<String, figment::value::Value>>() {
            for (key, value) in all_data {
                if !key.contains('.') {
                    if let Some(value_str) = value.into_string() {
                        // Only insert if the key doesn't already exist (preserves local parameters)
                        self.statements.entry(key).or_insert(value_str);
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
                    writeln!(f, "{}{key} = \"{value}\";", vec![" "; indent].join(""))?;
                } else {
                    write!(f, "\n{key} = \"{value}\";")?;
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
                write!(f, "{key} = \"{value}\"")?;
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
                write!(f, "{key} = \"{value}\"")?;
                first = false;
            }
            writeln!(f, "]")?;
        }

        Ok(())
    }
}

impl TryFrom<(Vec<AttrStmt<(String, String)>>, Vec<IDEq>)> for GlobalData {
    type Error = ();

    fn try_from(value: (Vec<AttrStmt<(String, String)>>, Vec<IDEq>)) -> Result<Self, Self::Error> {
        let mut statements = BTreeMap::new();
        let mut edge_statements = BTreeMap::new();
        let mut node_statements = BTreeMap::new();

        for attr_stmt in value.0 {
            match attr_stmt {
                AttrStmt::Graph(l) => {
                    for l in l.elems {
                        statements.extend(
                            l.into_iter()
                                .map(|(k, v)| (k, strip_quotes(&v).to_string())),
                        );
                    }
                }
                AttrStmt::Node(l) => {
                    for l in l.elems {
                        node_statements.extend(
                            l.into_iter()
                                .map(|(k, v)| (k, strip_quotes(&v).to_string())),
                        );
                    }
                }
                AttrStmt::Edge(l) => {
                    for l in l.elems {
                        edge_statements.extend(
                            l.into_iter()
                                .map(|(k, v)| (k, strip_quotes(&v).to_string())),
                        );
                    }
                }
            }
        }

        for e in value.1 {
            statements.insert(e.lhs, strip_quotes(&e.rhs).to_string());
        }

        let mut name = String::new();
        statements.retain(|k, v| {
            if k.as_str() == "name" {
                name = strip_quotes(v).to_string();
                false
            } else {
                true
            }
        });

        Ok(GlobalData {
            name,
            statements,
            edge_statements,
            node_statements,
        })
    }
}
