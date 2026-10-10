use feynkit_graph::expressions::evaluate_overall_factor;
use linnet::parser::GlobalData;
use symbolica::atom::Atom;

use crate::graph::GroupId;

use super::{Graph, parse::ToQuoted};

/// Graph-level payloads of a diagram converted from FeynKit.
#[derive(Clone, Debug)]
pub struct ParseData {
    pub name: String,
    pub overall_factor: Atom,
    pub projector: Atom,
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
            projector: Atom::one(),
            parameters: Vec::new(),
            num: Atom::one(),
            group_id: None,
            is_group_master: false,
        }
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

#[cfg(test)]
mod tests {
    use linnet::{half_edge::nodestore::NodeForest, parser::DotGraph};
    use symbolica::{atom::Atom, parser::ParseSettings};

    use super::evaluate_overall_factor;

    use crate::{
        dot,
        graph::{Graph, parse::IntoGraph},
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
        match dot!(digraph params_roundtrip {
            graph [
                overall_factor = 1;
                multiplicity_factor = 1;
                params = "a;b;c";
                projector = 1;
            ]
            edge [pdg=1000 num=1]
            node [num=1]
            ext [style=invis]
            ext -> v4 
            ext -> v5 
            v6 -> ext 
            v5 -> v4 [lmb_id=0];
            v6 -> v5 ;
            v4 -> v6 ;
        },"scalars")
        {
            Ok(g) => {
                let g: Graph = g;
                // g.to_dot_graph_with_settings().dot()
                let serialized = g.dot_serialize(&DotExportSettings::default());
                let parsed: DotGraph<NodeForest> = DotGraph::from_string(serialized).unwrap();

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
