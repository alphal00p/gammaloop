use eyre::Context;
use linnet::parser::GlobalData;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function, symbol,
};

use crate::graph::GroupId;

use super::{
    Graph,
    parse::{ParseGraph, StripParse, ToQuoted},
};

/// Evaluate the bookkeeping heads in a finalized diagram's overall factor.
pub fn evaluate_overall_factor(factor: AtomView<'_>) -> Atom {
    let mut result = factor.to_owned();
    for head in [
        "AutG",
        "CouplingsMultiplicity",
        "InternalFermionLoopSign",
        "ExternalFermionOrderingSign",
        "AntiFermionSpinSumSign",
        "NumeratorIndependentSymmetryGrouping",
    ] {
        for symbol in [
            symbol!(head),
            symbol!(&format!("feynkit_generator_factor::{head}")),
        ] {
            result = result
                .replace(function!(symbol, Atom::var(symbol!("x_"))).to_pattern())
                .with(Atom::var(symbol!("x_")).to_pattern());
        }
    }
    for head in [
        symbol!("NumeratorDependentGrouping"),
        symbol!("feynkit_generator::NumeratorDependentGrouping"),
        symbol!("feynkit_generator_factor::NumeratorDependentGrouping"),
    ] {
        result = result
            .replace(
                function!(
                    head,
                    Atom::var(symbol!("GraphId_")),
                    Atom::var(symbol!("ratio_")),
                    Atom::var(symbol!("GraphSymmetryFactor_"))
                )
                .to_pattern(),
            )
            .with(
                (Atom::var(symbol!("ratio_")) * Atom::var(symbol!("GraphSymmetryFactor_")))
                    .to_pattern(),
            );
    }
    result.expand()
}

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

        if !self.finalized_cuts.is_empty() {
            // Runtime DOT intentionally does not serialize canonical physical
            // cuts. Cross-section artifacts must round-trip through
            // FeynmanDiagram DOT, which owns the typed cut partitions.
            g.statements.insert(
                "canonical_cuts_required".to_owned(),
                "feynkit_diagram_dot".to_owned(),
            );
        }

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
    use symbolica::{atom::Atom, parser::ParseSettings};

    use super::evaluate_overall_factor;

    use crate::{
        finalized_runtime_dot,
        graph::{Graph, parse::IntoFinalizedRuntimeGraph},
        initialisation::test_initialise,
        processes::DotExportSettings,
    };

    #[test]
    fn evaluates_canonical_feynkit_factor_namespace() {
        test_initialise().unwrap();
        let factor = Atom::parse(
            "InternalFermionLoopSign(-1)*CouplingsMultiplicity(2)/AutG(4)\
             +NumeratorDependentGrouping(7,2,3)",
            "feynkit_generator_factor",
            ParseSettings::default(),
        )
        .unwrap();
        assert_eq!(
            evaluate_overall_factor(factor.as_view()),
            Atom::parse("11/2", "factor_expected", ParseSettings::default()).unwrap()
        );
    }

    #[test]
    fn params_roundtrip_in_global_data() {
        test_initialise().unwrap();
        match finalized_runtime_dot!(digraph params_roundtrip {
            graph [
                overall_factor = 1;
                multiplicity_factor = 1;
                params = "a;b;c";
                projector = 1;
            ]
            edge [pdg=1000 num=1]
            node [num=1]
            ext [style=invis]
            ext -> v4 [sink="{ufo_order:0}"]
            ext -> v5 [sink="{ufo_order:0}"]
            v6 -> ext [source="{ufo_order:0}"]
            v5 -> v4 [lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
            v6 -> v5 [source="{ufo_order:1}" sink="{ufo_order:2}"];
            v4 -> v6 [source="{ufo_order:2}" sink="{ufo_order:2}"];
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
