use linnet::parser::GlobalData;
use log::info;
use spenso::network::library::TensorLibraryData;
use symbolica::atom::{Atom, AtomCore};

use super::{parse::StripParse, Graph};

#[derive(Clone, Debug)]
pub struct ParseData {
    pub name: String,
    pub overall_factor: Atom,
    pub multiplicity_factor: Atom,
    pub color: Atom,
    pub colorless: Atom,
}

impl Default for ParseData {
    fn default() -> Self {
        ParseData {
            name: String::new(),
            overall_factor: Atom::one(),
            multiplicity_factor: Atom::one(),
            color: Atom::one(),
            colorless: Atom::one(),
        }
    }
}

impl ParseData {
    pub fn with_overall_factor(self, overall_factor: Atom) -> Self {
        ParseData {
            name: self.name,
            overall_factor,
            multiplicity_factor: self.multiplicity_factor,
            color: self.color,
            colorless: self.colorless,
        }
    }

    pub fn with_multiplicity_factor(self, multiplicity_factor: Atom) -> Self {
        ParseData {
            name: self.name,
            overall_factor: self.overall_factor,
            multiplicity_factor,
            color: self.color,
            colorless: self.colorless,
        }
    }

    pub fn with_color(self, color: Atom) -> Self {
        ParseData {
            name: self.name,
            overall_factor: self.overall_factor,
            multiplicity_factor: self.multiplicity_factor,
            color,
            colorless: self.colorless,
        }
    }

    pub fn with_colorless(self, colorless: Atom) -> Self {
        ParseData {
            name: self.name,
            overall_factor: self.overall_factor,
            multiplicity_factor: self.multiplicity_factor,
            color: self.color,
            colorless,
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

        if let Some(factor) = value.statements.get("multiplicity_factor") {
            parse_data = parse_data.with_multiplicity_factor(factor.strip_parse());
        }

        if let Some(color) = value.statements.get("num") {
            parse_data = parse_data.with_color(color.strip_parse());
        }

        if let Some(colorless) = value.statements.get("color_num") {
            parse_data = parse_data.with_colorless(colorless.strip_parse());
        }

        parse_data
    }
}

impl Graph {
    pub fn global_data(&self) -> GlobalData {
        let mut g = GlobalData::from(());

        // info!("Name: {}", self.name);
        g.add_name(self.name.clone());

        g.statements.insert(
            "num".to_string(),
            self.global_prefactor.color.to_canonical_string(),
        );
        g.statements.insert(
            "color_num".to_string(),
            self.global_prefactor.colorless.to_canonical_string(),
        );
        // g.statements.insert(
        //     "overall_factor".to_string(),
        //     self.global_prefactor.color.to_canonical_string(),
        // );
        g.statements.insert(
            "multiplicity_factor".to_string(),
            self.multiplicity.to_canonical_string(),
        );

        g
    }
}
