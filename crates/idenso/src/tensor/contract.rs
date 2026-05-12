#[cfg(test)]
pub mod test {
    use super::super::AbstractIndex;
    use crate::schoonschip::Schoonschip;
    use crate::test_support::test_initialize;
    use symbolica::{parse, symbol};

    #[test]
    fn schoonschip_simplifies_metrics() {
        test_initialize();
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
