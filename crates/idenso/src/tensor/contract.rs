#[cfg(test)]
pub mod test {
    use super::super::AbstractIndex;
    use crate::shorthands::schoonschip::Schoonschip;
    use crate::test_support::test_initialize;
    use spenso::{g, mink, p};

    #[test]
    fn schoonschip_simplifies_metrics() {
        test_initialize();
        let expr = g!(mink!(4, mu), mink!(4, nu)) * p!(mink!(4, nu));

        assert_eq!(expr.schoonschip_net::<AbstractIndex>(), p!(mink!(4, mu)));
    }
}
