#[cfg(test)]
pub mod test {
    use super::super::AbstractIndex;
    use crate::schoonschip::Schoonschip;
    use crate::test_support::test_initialize;
    use spenso::{g, mink, p};
    use symbolica::symbol;

    #[test]
    fn schoonschip_simplifies_metrics() {
        test_initialize();
        let (mu, nu) =
            symbol!("metric_simplify_mu", "metric_simplify_nu"; tags = ["spenso::index"]);
        let expr = g!(mink!(4, mu.clone()), mink!(4, nu.clone())) * p!(mink!(4, nu.clone()));

        assert_eq!(
            expr.schoonschip_net::<AbstractIndex>(),
            p!(mink!(4, mu.clone()))
        );
    }
}
