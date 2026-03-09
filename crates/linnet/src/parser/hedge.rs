use std::fmt::Display;

use dot_parser::ast::{CompassPt, Port};

use crate::half_edge::involution::{Flow, Hedge};

use super::{strip_quotes, subgraph_free::PortExt};

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct DotHedgeData {
    pub statement: Option<String>,
    pub id: Option<Hedge>,
    pub port_label: Option<String>,
    pub compasspt: Option<CompassPt>,
}

impl DotHedgeData {
    pub fn is_none(&self) -> bool {
        self.statement.is_none()
            && self.id.is_none()
            && self.port_label.is_none()
            && self.compasspt.is_none()
    }

    #[allow(clippy::useless_asref)]
    pub fn dot_serialize(&self) -> Option<String> {
        let mut out = String::new();
        let mut info = false;
        if let Some(statement) = &self.statement {
            info = true;
            out.push_str(statement);
        }

        if let Some(id) = &self.id {
            info = true;
            out.push_str(&format!(" [id={id}]"));
        }
        if info {
            Some(out)
        } else {
            None
        }
    }

    pub fn with_statement(mut self, statement: String) -> Self {
        self.statement = Some(statement);
        self
    }

    pub fn with_port(mut self, port: &Port) -> Self {
        self.port_label = port.id().map(|a| a.to_string());
        self.compasspt = port.compass().cloned();
        self.id = port.hedge();
        self
    }

    pub fn with_id(mut self, id: Hedge) -> Self {
        self.id = Some(id);
        self
    }
}

impl Display for DotHedgeData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(statement) = &self.statement {
            write!(f, "{statement}")?;
        }
        Ok(())
    }
}

impl From<Option<String>> for DotHedgeData {
    fn from(statement: Option<String>) -> Self {
        DotHedgeData {
            statement: statement.map(|s| strip_quotes(&s).to_string()),
            ..Default::default()
        }
    }
}

impl From<Option<&Port>> for DotHedgeData {
    fn from(port: Option<&Port>) -> Self {
        if let Some(port) = port {
            DotHedgeData::default().with_port(port)
        } else {
            DotHedgeData::default()
        }
    }
}

pub enum ParsingHedgePair {
    Unpaired {
        hedge: Option<Hedge>,
        flow: Flow,
        data: DotHedgeData,
    },
    Paired {
        source: Option<Hedge>,
        source_data: DotHedgeData,
        sink: Option<Hedge>,
        sink_data: DotHedgeData,
    },
}
