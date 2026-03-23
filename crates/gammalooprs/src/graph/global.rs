use eyre::Context;
use linnet::parser::GlobalData;
use spenso::network::library::TensorLibraryData;
use symbolica::atom::Atom;

use crate::{feyngen::diagram_generator::evaluate_overall_factor, graph::GroupId};

use super::{
    Graph,
    parse::{ParseGraph, StripParse, ToQuoted},
};

#[derive(Clone, Debug)]
pub struct ParseData {
    pub name: String,
    pub overall_factor: Atom,
    pub projectors: Option<Atom>,
    pub num: Atom,
    pub parameters: Vec<Atom>,
    pub group_id: Option<GroupId>,
    pub is_group_master: bool,
}

impl Default for ParseData {
    fn default() -> Self {
        ParseData {
            name: String::new(),
            overall_factor: Atom::one(),
            projectors: None,
            parameters: Vec::new(),
            num: Atom::one(),
            group_id: None,
            is_group_master: false,
        }
    }
}

impl ParseData {
    pub(crate) fn with_overall_factor(self, overall_factor: Atom) -> Self {
        ParseData {
            name: self.name,
            overall_factor,
            projectors: self.projectors,
            num: self.num,
            parameters: self.parameters,
            group_id: self.group_id,
            is_group_master: self.is_group_master,
        }
    }

    pub(crate) fn with_projectors(self, polarizations: Atom) -> Self {
        ParseData {
            name: self.name,
            overall_factor: self.overall_factor,
            projectors: Some(polarizations),
            num: self.num,
            parameters: self.parameters,
            group_id: self.group_id,
            is_group_master: self.is_group_master,
        }
    }

    pub(crate) fn with_num(self, num: Atom) -> Self {
        ParseData {
            name: self.name,
            overall_factor: self.overall_factor,
            projectors: self.projectors,
            num,
            parameters: self.parameters,
            group_id: self.group_id,
            is_group_master: self.is_group_master,
        }
    }
}

impl From<linnet::parser::GlobalData> for ParseData {
    fn from(value: linnet::parser::GlobalData) -> Self {
        let mut parse_data = ParseData {
            name: value.name,
            ..Default::default()
        };

        if let Some(factor) = value.statements.get("overall_factor") {
            parse_data = parse_data
                .with_overall_factor(factor.strip_parse().context("overall_factor").unwrap());
        }

        if let Some(polarizations) = value.statements.get("projector") {
            parse_data = parse_data
                .with_projectors(polarizations.strip_parse().context("projector").unwrap());
        }

        if let Some(factor) = value.statements.get("num") {
            parse_data = parse_data.with_num(factor.strip_parse().context("num").unwrap());
        }

        if let Some(is_group_master) = value.statements.get("is_group_master") {
            parse_data.is_group_master = is_group_master
                .strip_parse()
                .context("is_group_master")
                .unwrap();
        }

        if let Some(group_id) = value.statements.get("group_id") {
            parse_data.group_id =
                Some(GroupId(group_id.strip_parse().context("group_id").unwrap()));
        }

        if let Some(params) = value.statements.get("params") {
            let params: String = params.strip_parse().context("params").unwrap();
            parse_data.parameters = params
                .split(';')
                .map(str::trim)
                .filter(|param| !param.is_empty())
                .map(|param| {
                    param
                        .strip_parse()
                        .with_context(|| format!("params entry {param}"))
                        .unwrap()
                })
                .collect();
        }

        parse_data
    }
}

impl Graph {
    pub(crate) fn global_data(&self) -> GlobalData {
        let mut g = GlobalData::from(());

        // println!("Name: {}", self.name);
        g.add_name(self.name.clone());

        g.statements
            .insert("num".to_string(), self.global_prefactor.num.to_quoted());

        g.statements.insert(
            "projector".to_string(),
            self.global_prefactor.projector.to_quoted(),
        );
        // g.statements.insert(
        //     "overall_factor".to_string(),
        //     self.global_prefactor.color.to_canonical_string(),
        // );
        g.statements.insert(
            "overall_factor".to_string(),
            self.overall_factor.to_quoted(),
        );

        g.statements.insert(
            "overall_factor_evaluated".to_string(),
            evaluate_overall_factor(self.overall_factor.as_view()).to_quoted(),
        );

        if !self.param_builder.pairs.additional_params.params.is_empty() {
            let params = self
                .param_builder
                .pairs
                .additional_params
                .params
                .iter()
                .map(ToQuoted::to_quoted)
                .collect::<Vec<_>>()
                .join(";");
            g.statements.insert("params".to_string(), params);
        }

        g
    }
}

impl ParseGraph {
    pub(crate) fn global_data(&self) -> GlobalData {
        let mut g = GlobalData::from(());

        // println!("Name: {}", self.name);
        g.add_name(self.global_data.name.clone());

        g.statements
            .insert("num".to_string(), self.global_data.num.to_quoted());
        if let Some(proj) = &self.global_data.projectors {
            g.statements
                .insert("projector".to_string(), proj.to_quoted());
        }

        // g.statements.insert(
        //     "overall_factor".to_string(),
        //     self.global_prefactor.color.to_canonical_string(),
        // );
        g.statements.insert(
            "overall_factor".to_string(),
            self.global_data.overall_factor.to_quoted(),
        );

        if !self.global_data.parameters.is_empty() {
            let params = self
                .global_data
                .parameters
                .iter()
                .map(ToQuoted::to_quoted)
                .collect::<Vec<_>>()
                .join(";");
            g.statements.insert("params".to_string(), params);
        }

        g
    }
}

#[cfg(test)]
mod tests {
    use linnet::{
        half_edge::nodestore::NodeStorageVec,
        parser::{DotGraph, DotVertexData},
    };

    use crate::{
        dot,
        graph::{Graph, parse::IntoGraph},
        initialisation::test_initialise,
        processes::DotExportSettings,
    };

    #[test]
    fn params_roundtrip_in_global_data() {
        test_initialise().unwrap();
        match dot!(digraph params_roundtrip {
            graph [
                overall_factor = 1;
                multiplicity_factor = 1;
                params = "a;b;c";
            ]
            edge [pdg=1000]
            ext [style=invis]
            ext -> v4
            ext -> v5
            v6 -> ext
            v5 -> v4 [lmb_index=0];
            v6 -> v5;
            v4 -> v6;
        },"scalars")
        {
            Ok(g) => {
                let g: Graph = g;
                // g.to_dot_graph_with_settings().dot()
                let serialized = g.dot_serialize(&DotExportSettings::default());
                let parsed: DotGraph<NodeStorageVec<DotVertexData>> =
                    DotGraph::from_string(serialized).unwrap();

                assert_eq!(
                    parsed.global_data.statements.get("params"),
                    Some(&"a;b;c".to_string())
                );
            }
            Err(e) => {
                eprintln!("Graph parsing failed: {:?}", e);
            }
        }
    }
}
