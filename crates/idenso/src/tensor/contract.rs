#[cfg(test)]
pub mod test {
    use super::super::{AbstractIndex, Contract, ShadowedStructure, SymbolicTensor};
    use crate::representations::initialize;
    use crate::schoonschip::Schoonschip;
    use spenso::network::library::symbolic::ETS;
    use spenso::network::parsing::StructureFromAtom;
    use symbolica::{parse, symbol};

    #[test]
    fn parse() {
        let _ = ETS.metric;
        let expr = parse!("g(mink(4,6),mink(4,7))");

        let parsed = ShadowedStructure::<AbstractIndex>::parse(expr.as_view()).unwrap();
        let structure = SymbolicTensor::from_permuted(&parsed).unwrap();

        Contract::<SymbolicTensor, ()>::contract(&structure, &structure).unwrap();
    }

    #[test]
    fn schoonschip_simplifies_metrics() {
        initialize();
        symbol!("metric_simplify_mu", "metric_simplify_nu"; tags = ["spenso::index"]);
        let expr = parse!(
            "spenso::g(spenso::mink(4,metric_simplify_mu),spenso::mink(4,metric_simplify_nu))*p(spenso::mink(4,metric_simplify_nu))"
        );

        assert_eq!(
            expr.schoonschip_net::<AbstractIndex>(),
            parse!("p(spenso::mink(4,metric_simplify_mu))")
        );
    }
}
