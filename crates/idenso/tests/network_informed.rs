use idenso::{metric::MetricSimplifier, representations::initialize, schoonschip::Schoonschip};
use spenso::structure::abstract_index::AbstractIndex;
use symbolica::{atom::AtomCore, id::Pattern, parse, symbol, transformer::Transformer};

fn run_network_informed() {
    initialize();
    symbol!("mu1", "mu2", "mu3", "mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "mu10", "mu11";
        tags=["spenso::index"]);

    let mut r = parse!(
        "vx(1,-k0, k0-k1, k1, k10, mu1, mu8)
    vx(2,-k1, k2, k1-k2, mu1, mu2, mu9)
    vx(3,-k2, k3, k2-k3, mu2, mu3, mu10)
    vx(4,-k3, k4, k3-k4, mu3, mu4, mu11)
    vx(5,-k4, k0, k4-k0, mu4, k20, mu5)
    vx(6,-k4+k0, -k3+k4, k3-k0, mu5, mu11, mu6)
    vx(7,-k3+k0, -k2+k3, k2-k0, mu6, mu10, mu7)
    vx(8,-k2+k0, -k1+k2, k1-k0, mu7, mu9, mu8)"
    );

    let gluon_rule = parse!(
        "(- spenso::mink(4,mu1_, mu3_) * spenso::mink(4,k1_, mu2_)
                + spenso::mink(4,mu1_, mu2_) * spenso::mink(4,k1_, mu3_)
                + spenso::mink(4,mu2_, mu3_) * spenso::mink(4,k2_, mu1_)
                - spenso::mink(4,mu1_, mu2_) * spenso::mink(4,k2_, mu3_)
                - spenso::mink(4,mu2_, mu3_) * spenso::mink(4,k3_, mu1_)
                + spenso::mink(4,mu1_, mu3_) * spenso::mink(4,k3_, mu2_)
                )"
    )
    .to_pattern();

    let rhs_subs = Pattern::Transformer(Box::new((
        Some(gluon_rule),
        vec![Transformer::Map(Box::new(move |input, _state, out| {
            *out = input.expand().simplify_metrics(); //.to_dots();
            Ok(())
        }))],
    )));

    for i in 0..8 {
        r = r
            .replace(parse!(format!(
                "vx({}, k1_, k2_, k3_, mu1_, mu2_, mu3_)",
                i + 1
            )))
            .level_range((0, Some(0)))
            .rhs_cache_size(1000)
            .with(&rhs_subs);
    }

    let _ = r.schoonschip_net::<AbstractIndex>();
}

#[test]
fn network_informed() {
    run_network_informed();
}
