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
    pub group_id: Option<GroupId>,
    pub is_group_master: bool,
}

impl Default for ParseData {
    fn default() -> Self {
        ParseData {
            name: String::new(),
            overall_factor: Atom::one(),
            projectors: None,
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
            group_id: self.group_id,
            is_group_master: self.is_group_master,
        }
    }

    pub(crate) fn with_polarizations(self, polarizations: Atom) -> Self {
        ParseData {
            name: self.name,
            overall_factor: self.overall_factor,
            projectors: Some(polarizations),
            num: self.num,
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
            group_id: self.group_id,
            is_group_master: self.is_group_master,
        }
    }
}

impl From<linnet::parser::GlobalData> for ParseData {
    fn from(value: linnet::parser::GlobalData) -> Self {
        let mut parse_data = ParseData::default();

        parse_data.name = value.name;

        if let Some(factor) = value.statements.get("overall_factor") {
            parse_data = parse_data
                .with_overall_factor(factor.strip_parse().context("overall_factor").unwrap());
        }

        if let Some(polarizations) = value.statements.get("projector") {
            parse_data = parse_data
                .with_polarizations(polarizations.strip_parse().context("projector").unwrap());
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

        g
    }
}
