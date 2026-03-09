use symbolica::{

    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    id::{Match, MatchSettings, MatchStack, Pattern, Replacement},
    parse, symbol,
    transformer::Transformer,
};

fn run_network_informed() {
    let d = symbol!("d"; Symmetric, Linear); // metric

    symbol!("mu1", "mu2", "mu3", "mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "mu10", "mu11";
        tags=["spenso::index"]);

    let k1_ = symbol!("k1_");
    let muw1_ = symbol!("muw1_", tags = ["spenso::index"]);
    let x___ = symbol!("x___");

    let mut r = parse!(
        "vx(1,-k(0), k(0)-k(1), k(1), k(10), mu1, mu8)
    vx(2,-k(1), k(2), k(1)-k(2), mu1, mu2, mu9)
    vx(3,-k(2), k(3), k(2)-k(3), mu2, mu3, mu10)
    vx(4,-k(3), k(4), k(3)-k(4), mu3, mu4, mu11)
    vx(5,-k(4), k(0), k(4)-k(0), mu4, k(20), mu5)
    vx(6,-k(4)+k(0), -k(3)+k(4), k(3)-k(0), mu5, mu11, mu6)
    vx(7,-k(3)+k(0), -k(2)+k(3), k(2)-k(0), mu6, mu10, mu7)
    vx(8,-k(2)+k(0), -k(1)+k(2), k(1)-k(0), mu7, mu9, mu8)"
    );

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
        "(- d(mu1_, mu3_) * d(k1_, mu2_)
                + d(mu1_, mu2_) * d(k1_, mu3_)
                + d(mu2_, mu3_) * d(k2_, mu1_)
                - d(mu1_, mu2_) * d(k2_, mu3_)
                - d(mu2_, mu3_) * d(k3_, mu1_)
                + d(mu1_, mu3_) * d(k3_, mu2_)
                )"
    )
    .to_pattern();

    let rhs_subs = Pattern::Transformer(Box::new((
        Some(gluon_rule),
        vec![Transformer::Map(Box::new(move |input, _state, out| {
            input
                .expand()
                .replace(parse!("d(muw1_,k1_)*d(muw1_,k2_)"))
                .level_range((0, Some(0)))
                .repeat()
                .with(parse!("d(k1_,k2_)"))
                .replace(parse!("d(muw1_,k1_)^2"))
                .level_range((0, Some(0)))
                .with(parse!("d(k1_,k1_)"))
                .replace(parse!("d(muw1_,muw1_)"))
                .level_range((0, Some(0)))
                .with_into(4, out);

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
        r = r.expand();

        // 479.0 ms, 454.6ms when momenta are symbols, 5% faster if not using functions
        r = r
            .replace(parse!("d(muw1_, k1_)*vx(x___,muw1_,y___)"))
            .level_range((0, Some(0)))
            .repeat()
            .with(parse!("vx(x___,k1_,y___)"));

        // 548.9 ms 14% slower than above, 492.8 ms when momenta are symbols, 8% slower than above
        // 516ms when momenta symbols and no replace_map
        // r = r
        //     .replace(parse!("d(muw1_, k1_)*x___"))
        //     .level_range((0, Some(0)))
        //     .repeat()
        //     .with_map(move |m| {
        //         let a1 = m.get(muw1_).unwrap().to_atom();
        //         let a2 = m.get(k1_).unwrap().to_atom();
        //         let dest = m.get(x___).unwrap().to_atom(); // PREVENT!

        //         dest.replace(a1)
        //             .level_range((1, Some(1)))
        //             .rhs_cache_size(0)
        //             .with(a2)

        //         // if let Some(Match::Single(a1)) = m.get(muw1_) {
        //         //     if let Some(Match::Single(a2)) = m.get(k1_) {
        //         //         l

        //         //         return dest.replace_map(|input, _, out| {
        //         //             if input == *a1 {
        //         //                 out.set_from_view(a2);
        //         //             }
        //         //         });
        //         //     }
        //         // }

        //         // Atom::Zero
        //     });

        println!("{}", r.nterms());
        //println!("{:>}", r);
    }
}

#[test]
fn network_informed() {
    run_network_informed();
}
