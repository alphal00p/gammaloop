#[cfg(test)]
pub mod test {
    use super::super::AbstractIndex;
    use crate::representations::initialize;
    use crate::schoonschip::Schoonschip;
    use symbolica::{parse, symbol};

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
