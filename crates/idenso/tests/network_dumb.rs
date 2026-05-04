use idenso::{
    representations::initialize,
    schoonschip::{Schoonschip, SchoonschipSettings},
};
use spenso::{
    shadowing::symbolica_utils::{AtomCoreExt, SpensoPrintSettings},
    structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
    },
};
use symbolica::{atom::AtomCore, parse, symbol};

fn contains_index(result: &str, index: &str) -> bool {
    result.contains(&format!("{index})")) || result.contains(&format!("{index},"))
}

fn run_network_informed() {
    initialize();
    let _mink = Minkowski {}.new_rep(4);

    symbol!("mu1", "mu2", "mu3", "mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "mu10", "mu11";
        tags=["spenso::index"]);

    let settings = SchoonschipSettings::full();

    let metric_vector = parse!("spenso::mink(4, mu1, mu2)*k1(spenso::mink(4, mu1))")
        .schoonschip_with_net::<false, true, AbstractIndex>(&settings)
        .to_bare_ordered_string();
    assert_eq!(metric_vector, "k1(mink(4,mu2))");

    let metric_metric = parse!("spenso::mink(4, mu1, mu2)*spenso::mink(4, mu2, mu3)")
        .schoonschip_with_net::<false, true, AbstractIndex>(&settings)
        .to_bare_ordered_string();
    assert_eq!(metric_metric, "g(mink(4,mu1),mink(4,mu3))");

    let metric_chain = parse!(
        "spenso::mink(4, mu1, mu2)
            * spenso::mink(4, mu2, mu3)
            * k1(spenso::mink(4, mu3))"
    )
    .schoonschip_with_net::<false, true, AbstractIndex>(&settings)
    .to_bare_ordered_string();
    assert_eq!(metric_chain, "k1(mink(4,mu1))");
    assert!(!contains_index(&metric_chain, "mu2"));
    assert!(!contains_index(&metric_chain, "mu3"));

    let metric_into_vertex = parse!(
        "spenso::mink(4, mu1, mu2)
            * vx(3, -k2, k3, k2-k3, spenso::mink(4, mu1), spenso::mink(4, mu3), spenso::mink(4, mu10))"
    )
    .schoonschip_with_net::<false, true, AbstractIndex>(&settings)
    .to_bare_ordered_string();

    assert!(!contains_index(&metric_into_vertex, "mu1"));
    assert!(contains_index(&metric_into_vertex, "mu2"));
}

#[test]
fn network_informed() {
    run_network_informed();
}

#[test]
fn spenso_bare_symb_vertex_substitution() {
    initialize();
    let _mink = Minkowski {}.new_rep(4);

    let (mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10, mu11) = symbol!("mu1", "mu2", "mu3", "mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "mu10", "mu11";
        tags=["spenso::index"]);

    symbol!("k";
        tags=["spenso::tensor","spenso::rank1"]);

    let v1 =
        parse!("vx(1,-k(0), k(0)-k(1), k(1), k(10), spenso::mink(4,mu1), spenso::mink(4,mu8))");
    let v2 = parse!(
        "vx(2,-k(1), k(2), k(1)-k(2), spenso::mink(4,mu1), spenso::mink(4,mu2), spenso::mink(4,mu9))"
    );
    let v3 = parse!(
        "vx(3,-k(2), k(3), k(2)-k(3), spenso::mink(4,mu2), spenso::mink(4,mu3), spenso::mink(4,mu10))"
    );
    let v4 = parse!(
        "vx(4,-k(3), k(4), k(3)-k(4), spenso::mink(4,mu3), spenso::mink(4,mu4), spenso::mink(4,mu11))"
    );
    let v5 =
        parse!("vx(5,-k(4), k(0), k(4)-k(0), spenso::mink(4,mu4), k(20), spenso::mink(4,mu5))");
    let v6 = parse!(
        "vx(6,-k(4)+k(0), -k(3)+k(4), k(3)-k(0), spenso::mink(4,mu5), spenso::mink(4,mu11), spenso::mink(4,mu6))"
    );
    let v7 = parse!(
        "vx(7,-k(3)+k(0), -k(2)+k(3), k(2)-k(0), spenso::mink(4,mu6), spenso::mink(4,mu10), spenso::mink(4,mu7))"
    );
    let v8 = parse!(
        "vx(8,-k(2)+k(0), -k(1)+k(2), k(1)-k(0), spenso::mink(4,mu7), spenso::mink(4,mu9), spenso::mink(4,mu8))"
    );

    let gluon_rule = parse!(
        "(- spenso::g(mu1_,mu3_) * spenso::g(k1_,mu2_)
                + spenso::g(mu1_,mu2_) * spenso::g(k1_,mu3_)
                + spenso::g(mu2_,mu3_) * spenso::g(k2_,mu1_)
                - spenso::g(mu1_,mu2_) * spenso::g(k2_,mu3_)
                - spenso::g(mu2_,mu3_) * spenso::g(k3_,mu1_)
                + spenso::g(mu1_,mu3_) * spenso::g(k3_,mu2_)
                )"
    )
    .to_pattern();

    let mut r = v1 * v2 * v3 * v4 * v5; // * v6 * v7 * v8;

    for i in 0..8 {
        let gluon_rule = gluon_rule.clone();
        r = r
            .replace(parse!(format!(
                "vx({}, k1_, k2_, k3_, mu1_, mu2_, mu3_)",
                i + 1
            )))
            .level_range((0, Some(0)))
            .rhs_cache_size(1000)
            .with_map(move |matches| {
                gluon_rule
                    .replace_wildcards_with_matches(matches)
                    .normalize_dots()
            });
    }

    let result = r.to_string();
    assert!(!result.contains("vx("), "{result}");
    assert!(result.contains("mink"), "{result}");
    let mut settings = SpensoPrintSettings::compact().nice_symbolica();
    settings.max_line_length = Some(80);
    println!("in:{}", r.printer(settings));

    let out = r.schoonschip_with_net::<false, false, AbstractIndex>(
        &SchoonschipSettings::partial().with_expanded_contracted_sums(),
    );

    println!("out:{}", out.printer(settings));

    for mu in [mu1, mu2, mu3, mu4] {
        assert!(
            out.replace(mu).match_iter().next().is_none(),
            "{}",
            out.printer(settings)
        );
    }
}
