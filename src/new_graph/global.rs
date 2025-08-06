use linnet::parser::GlobalData;
use log::info;
use spenso::network::library::TensorLibraryData;
use symbolica::atom::{Atom, AtomCore};

use super::{
    parse::{StripParse, ToQuoted},
    Graph,
};

#[derive(Clone, Debug)]
pub struct ParseData {
    pub name: String,
    pub overall_factor: Atom,
    pub projectors: Option<Atom>,
    pub num: Atom,
}

impl Default for ParseData {
    fn default() -> Self {
        ParseData {
            name: String::new(),
            overall_factor: Atom::one(),
            projectors: None,
            num: Atom::one(),
        }
    }
}

impl ParseData {
    pub fn with_overall_factor(self, overall_factor: Atom) -> Self {
        ParseData {
            name: self.name,
            overall_factor,
            projectors: self.projectors,
            num: self.num,
        }
    }

    pub fn with_polarizations(self, polarizations: Atom) -> Self {
        ParseData {
            name: self.name,
            overall_factor: self.overall_factor,
            projectors: Some(polarizations),
            num: self.num,
        }
    }

    pub fn with_num(self, num: Atom) -> Self {
        ParseData {
            name: self.name,
            overall_factor: self.overall_factor,
            projectors: self.projectors,
            num,
        }
    }
}

impl From<linnet::parser::GlobalData> for ParseData {
    fn from(value: linnet::parser::GlobalData) -> Self {
        let mut parse_data = ParseData::default();

        parse_data.name = value.name;

        if let Some(factor) = value.statements.get("overall_factor") {
            parse_data = parse_data.with_overall_factor(factor.strip_parse());
        }

        if let Some(polarizations) = value.statements.get("projector") {
            parse_data = parse_data.with_polarizations(polarizations.strip_parse());
        }

        if let Some(factor) = value.statements.get("num") {
            parse_data = parse_data.with_num(factor.strip_parse());
        }

        parse_data
    }
}

impl Graph {
    pub fn global_data(&self) -> GlobalData {
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

        g
    }
}
