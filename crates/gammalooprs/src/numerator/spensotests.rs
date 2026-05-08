use idenso::schoonschip::{Schoonschip, SchoonschipSettings};
use symbolica::parse_lit;

use crate::{
    initialisation::test_initialise, numerator::aind::Aind, utils::symbolica_ext::LogPrint,
};

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

    println!("{}", expr.schoonschip().log_print(Some(120)));
    println!("{}", expr.schoonschip().schoonschip().log_print(Some(120)));
    println!(
        "{}",
        expr.schoonschip_with_net::<true, false, Aind>(&SchoonschipSettings::breadth_first(Some(
            1
        )))
        .log_print(Some(120))
    );
}
