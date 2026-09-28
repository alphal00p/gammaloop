use idenso::{CookMode, CookSettings, dirac::GammaSimplifySettings, tensor::SymbolicTensor};
use std::sync::Arc;
use symbolica::parse_lit;

use crate::initialisation::test_initialise;
use spenso::shadowing::symbolica_utils::LogPrint;

#[test]
fn algebra() {
    test_initialise().unwrap();

    let expr = parse_lit!(
        UFO::GC_11
            ^ 2 * Q(7, spenso::mink(dim, edge(7, 1)))
                * spenso::g(spenso::mink(dim, hedge(14)), spenso::mink(dim, hedge(15)))
                * spenso::g(spenso::coad(8, hedge(14)), spenso::coad(8, hedge(15)))
                * spenso::g(
                    spenso::dind(spenso::cof(3, hedge(13))),
                    spenso::cof(3, hedge(12))
                )
                * spenso::gamma(
                    spenso::bis(4, hedge(4)),
                    spenso::bis(4, hedge(13)),
                    spenso::mink(dim, hedge(15))
                )
                * spenso::gamma(
                    spenso::bis(4, hedge(12)),
                    spenso::bis(4, hedge(9)),
                    spenso::mink(dim, hedge(14))
                )
                * spenso::gamma(
                    spenso::bis(4, hedge(13)),
                    spenso::bis(4, hedge(12)),
                    spenso::mink(dim, edge(7, 1))
                )
                * spenso::t(
                    spenso::coad(8, hedge(14)),
                    spenso::cof(3, hedge(9)),
                    spenso::dind(spenso::cof(3, hedge(12)))
                )
                * spenso::t(
                    spenso::coad(8, hedge(15)),
                    spenso::cof(3, hedge(13)),
                    spenso::dind(spenso::cof(3, hedge(4)))
                )
    );

    println!("{}", expr.log_print(Some(120)));

    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let source = SymbolicTensor::infer(cooking.try_cook(expr.as_view()).unwrap()).unwrap();
    let contracted = Arc::new(source.contract(Default::default()).unwrap());
    println!(
        "{}",
        contracted
            .resolved()
            .unwrap()
            .expression()
            .log_print(Some(120))
    );
    let simplified = contracted
        .simplify_gamma(GammaSimplifySettings::default())
        .unwrap();
    println!(
        "{}",
        simplified
            .resolved()
            .unwrap()
            .expression()
            .log_print(Some(120))
    );
}
