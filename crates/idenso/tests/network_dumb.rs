use idenso::{
    representations::initialize,
    schoonschip::{Schoonschip, SchoonschipContractionOrder, SchoonschipSettings},
};
use spenso::{
    shadowing::symbolica_utils::{AtomCoreExt, SpensoPrintSettings},
    structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
    },
};
use symbolica::{
    atom::{Atom, AtomCore},
    parse, symbol,
};

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

fn substituted_three_vertex_reproducer() -> (Atom, [(&'static str, Atom); 3]) {
    initialize();
    let _mink = Minkowski {}.new_rep(4);

    let (mu1, _mu2, _mu7, mu8, mu9) =
        symbol!("mu1", "mu2", "mu7", "mu8", "mu9"; tags=["spenso::index"]);

    symbol!("k";
        tags=["spenso::tensor","spenso::rank1"]);

    let v1 =
        parse!("vx(1,-k(0), k(0)-k(1), k(1), k(10), spenso::mink(4,mu1), spenso::mink(4,mu8))");
    let v2 = parse!(
        "vx(2,-k(1), k(2), k(1)-k(2), spenso::mink(4,mu1), spenso::mink(4,mu2), spenso::mink(4,mu9))"
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

    let mut r = v1 * v2 * v8;

    for i in [1, 2, 8] {
        let gluon_rule = gluon_rule.clone();
        r = r
            .replace(parse!(format!("vx({i}, k1_, k2_, k3_, mu1_, mu2_, mu3_)",)))
            .level_range((0, Some(0)))
            .rhs_cache_size(1000)
            .with_map(move |matches| {
                gluon_rule
                    .replace_wildcards_with_matches(matches)
                    .normalize_dots()
            });
    }

    (
        r,
        [
            ("mu1", mu1.into()),
            ("mu8", mu8.into()),
            ("mu9", mu9.into()),
        ],
    )
}

fn residual_dummy_names<'a>(out: &Atom, dummies: &'a [(&'a str, Atom)]) -> Vec<&'a str> {
    dummies
        .iter()
        .filter_map(|(name, dummy)| {
            out.replace(dummy.clone())
                .match_iter()
                .next()
                .map(|_| *name)
        })
        .collect()
}

fn print_three_vertex_method(name: &str, out: Atom, dummies: &[(&str, Atom)]) {
    let residuals = residual_dummy_names(&out, dummies);
    println!(
        "{name:<56} ok={:<5} residual={residuals:?} terms={} bytes={}",
        residuals.is_empty(),
        out.nterms(),
        out.as_view().get_byte_size()
    );
}

fn cleanup_with_smallest_degree(mut result: Atom) -> Atom {
    let cleanup_settings = SchoonschipSettings::partial()
        .with_expanded_contracted_sums()
        .with_contraction_order(SchoonschipContractionOrder::SmallestDegree);
    for _ in 0..4 {
        let next = result.schoonschip_with_net::<false, false, AbstractIndex>(&cleanup_settings);
        if next == result {
            break;
        }
        result = next;
    }
    result
}

fn print_two_dummy_method(name: &str, out: Atom, mu1: &Atom, mu9: &Atom) {
    let dummies = [("mu1", mu1.clone()), ("mu9", mu9.clone())];
    print_three_vertex_method(name, out, &dummies);
}

#[test]
fn min_product_terms_three_vertex_still_has_residual_after_boundary_cleanup() {
    let (r, dummies) = substituted_three_vertex_reproducer();
    let out = r.schoonschip_with_net::<false, false, AbstractIndex>(
        &SchoonschipSettings::partial()
            .with_expanded_contracted_sums()
            .with_contraction_order(SchoonschipContractionOrder::MinProductTerms),
    );

    assert_eq!(residual_dummy_names(&out, &dummies), ["mu9"]);
}

#[test]
fn compare_two_slot_boundary_shapes() {
    initialize();
    let _mink = Minkowski {}.new_rep(4);

    let (mu1, mu9) = symbol!("mu1", "mu9"; tags=["spenso::index"]);
    symbol!("k"; tags=["spenso::tensor","spenso::rank1"]);

    let mu1: Atom = mu1.into();
    let mu9: Atom = mu9.into();
    let settings = SchoonschipSettings::partial()
        .with_expanded_contracted_sums()
        .with_contraction_order(SchoonschipContractionOrder::MinProductTerms);

    let cases = [
        (
            "sum side scalar metric terms times simple target sum",
            parse!(
                "(a * spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9))
                   + b * spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9)))
                 * (spenso::g(k(2), spenso::mink(4,mu1))
                   * spenso::g(k(3), spenso::mink(4,mu9))
                   + spenso::g(k(4), spenso::mink(4,mu1))
                   * spenso::g(k(5), spenso::mink(4,mu9)))"
            ),
        ),
        (
            "metric times simple vector product",
            parse!(
                "spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9))
                 * spenso::g(k(0), spenso::mink(4,mu1))
                 * spenso::g(k(1), spenso::mink(4,mu9))"
            ),
        ),
        (
            "metric times vector product with summed momenta",
            parse!(
                "spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9))
                 * spenso::g(k(0)-k(1), spenso::mink(4,mu1))
                 * spenso::g(k(1)-k(0), spenso::mink(4,mu9))"
            ),
        ),
        (
            "sum side metric term times simple target sum",
            parse!(
                "(spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9))
                   + spenso::g(k(0), spenso::mink(4,mu1))
                     * spenso::g(k(1), spenso::mink(4,mu9)))
                 * (spenso::g(k(2), spenso::mink(4,mu1))
                   * spenso::g(k(3), spenso::mink(4,mu9))
                   + spenso::g(k(4), spenso::mink(4,mu1))
                   * spenso::g(k(5), spenso::mink(4,mu9)))"
            ),
        ),
        (
            "sum side metric term times summed-momentum target sum",
            parse!(
                "(spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9))
                   + spenso::g(k(0), spenso::mink(4,mu1))
                     * spenso::g(k(1), spenso::mink(4,mu9)))
                 * (spenso::g(k(2)-k(3), spenso::mink(4,mu1))
                   * spenso::g(k(3)-k(2), spenso::mink(4,mu9))
                   + spenso::g(k(4)-k(5), spenso::mink(4,mu1))
                   * spenso::g(k(5)-k(4), spenso::mink(4,mu9)))"
            ),
        ),
    ];

    println!("\ntwo-slot contraction boundary shape comparison");
    for (name, expr) in cases {
        print_two_dummy_method(
            &format!("{name} / normalize_dots"),
            expr.normalize_dots(),
            &mu1,
            &mu9,
        );
        print_two_dummy_method(
            &format!("{name} / network"),
            expr.schoonschip_with_net::<false, false, AbstractIndex>(&settings),
            &mu1,
            &mu9,
        );
        print_two_dummy_method(
            &format!("{name} / expanded-input network"),
            expr.expand()
                .schoonschip_with_net::<false, false, AbstractIndex>(&settings),
            &mu1,
            &mu9,
        );
    }
}

#[test]
fn metric_sum_boundary_uses_pattern_schoonschip_cleanup() {
    initialize();
    let _mink = Minkowski {}.new_rep(4);

    let (mu1, mu9) = symbol!("mu1", "mu9"; tags=["spenso::index"]);
    symbol!("k"; tags=["spenso::tensor","spenso::rank1"]);

    let mu1: Atom = mu1.into();
    let mu9: Atom = mu9.into();
    let dummies = [("mu1", mu1.clone()), ("mu9", mu9.clone())];
    let settings = SchoonschipSettings::partial()
        .with_expanded_contracted_sums()
        .with_contraction_order(SchoonschipContractionOrder::MinProductTerms);

    let target_after_metric_identification = parse!(
        "spenso::g(k(2), spenso::mink(4,mu9))
         * spenso::g(k(3), spenso::mink(4,mu9))
         + spenso::g(k(4), spenso::mink(4,mu9))
         * spenso::g(k(5), spenso::mink(4,mu9))"
    );
    assert_eq!(
        residual_dummy_names(
            &target_after_metric_identification.normalize_dots(),
            &dummies
        ),
        ["mu9"]
    );
    assert!(
        residual_dummy_names(
            &target_after_metric_identification
                .schoonschip_with_net::<false, false, AbstractIndex>(&settings),
            &dummies
        )
        .is_empty()
    );

    let boundary_expression = parse!(
        "(a * spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9))
           + b * spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9)))
         * (spenso::g(k(2), spenso::mink(4,mu1))
           * spenso::g(k(3), spenso::mink(4,mu9))
           + spenso::g(k(4), spenso::mink(4,mu1))
           * spenso::g(k(5), spenso::mink(4,mu9)))"
    );
    assert!(
        residual_dummy_names(
            &boundary_expression.schoonschip_with_net::<false, false, AbstractIndex>(&settings),
            &dummies
        )
        .is_empty()
    );
    assert!(
        residual_dummy_names(
            &boundary_expression
                .expand()
                .schoonschip_with_net::<false, false, AbstractIndex>(&settings),
            &dummies
        )
        .is_empty()
    );
}

#[test]
fn non_linear_metric_simplifies_summed_momentum_boundary_without_expansion() {
    initialize();
    let _mink = Minkowski {}.new_rep(4);

    let (mu1, mu9) = symbol!("mu1", "mu9"; tags=["spenso::index"]);
    symbol!("k"; tags=["spenso::tensor","spenso::rank1"]);

    let mu1: Atom = mu1.into();
    let mu9: Atom = mu9.into();
    let dummies = [("mu1", mu1.clone()), ("mu9", mu9.clone())];
    let settings = SchoonschipSettings::partial()
        .with_expanded_contracted_sums()
        .with_contraction_order(SchoonschipContractionOrder::MinProductTerms);

    let metric_identified_target = parse!(
        "spenso::g(k(2)-k(3), spenso::mink(4,mu9))
         * spenso::g(k(3)-k(2), spenso::mink(4,mu9))
         + spenso::g(k(4)-k(5), spenso::mink(4,mu9))
         * spenso::g(k(5)-k(4), spenso::mink(4,mu9))"
    );
    assert!(residual_dummy_names(&metric_identified_target.schoonschip(), &dummies).is_empty());
    assert!(
        residual_dummy_names(&metric_identified_target.expand().schoonschip(), &dummies).is_empty()
    );

    let boundary_expression = parse!(
        "(spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9))
           + spenso::g(k(0), spenso::mink(4,mu1))
             * spenso::g(k(1), spenso::mink(4,mu9)))
         * (spenso::g(k(2)-k(3), spenso::mink(4,mu1))
           * spenso::g(k(3)-k(2), spenso::mink(4,mu9))
           + spenso::g(k(4)-k(5), spenso::mink(4,mu1))
           * spenso::g(k(5)-k(4), spenso::mink(4,mu9)))"
    );

    assert!(
        residual_dummy_names(
            &boundary_expression.schoonschip_with_net::<false, false, AbstractIndex>(&settings),
            &dummies
        )
        .is_empty()
    );
    assert!(
        residual_dummy_names(
            &boundary_expression
                .expand()
                .schoonschip_with_net::<false, false, AbstractIndex>(&settings),
            &dummies
        )
        .is_empty()
    );
}

#[test]
fn metric_vector_product_with_free_metric_slot_simplifies_in_bare_cleanup() {
    initialize();
    let _mink = Minkowski {}.new_rep(4);

    let (_mu7, mu9) = symbol!("mu7", "mu9"; tags=["spenso::index"]);
    symbol!("k"; tags=["spenso::tensor","spenso::rank1"]);

    let mu9: Atom = mu9.into();
    let dummies = [("mu9", mu9.clone())];
    let settings = SchoonschipSettings::partial().with_expanded_contracted_sums();
    let expr = parse!(
        "spenso::g(spenso::mink(4,mu7), spenso::mink(4,mu9))
         * spenso::g(k(0)-k(1), spenso::mink(4,mu9))"
    );

    assert!(residual_dummy_names(&expr.schoonschip(), &dummies).is_empty());
    assert!(
        residual_dummy_names(
            &expr.schoonschip_with_net::<false, false, AbstractIndex>(&settings),
            &dummies
        )
        .is_empty()
    );
}

#[test]
fn compare_three_vertex_residual_methods() {
    let (r, dummies) = substituted_three_vertex_reproducer();
    let orders = [
        (
            "smallest_degree",
            SchoonschipContractionOrder::SmallestDegree,
        ),
        ("largest_degree", SchoonschipContractionOrder::LargestDegree),
        (
            "min_largest_operand_bytes",
            SchoonschipContractionOrder::MinLargestOperandBytes,
        ),
        (
            "min_product_terms",
            SchoonschipContractionOrder::MinProductTerms,
        ),
        (
            "min_product_bytes",
            SchoonschipContractionOrder::MinProductBytes,
        ),
        (
            "smallest_degree_min_largest_operand_bytes",
            SchoonschipContractionOrder::SmallestDegreeMinLargestOperandBytes,
        ),
        (
            "smallest_degree_min_product_terms",
            SchoonschipContractionOrder::SmallestDegreeMinProductTerms,
        ),
        (
            "smallest_degree_min_product_bytes",
            SchoonschipContractionOrder::SmallestDegreeMinProductBytes,
        ),
    ];

    println!("\nthree-vertex residual contraction comparison");
    print_three_vertex_method("normalize_dots", r.normalize_dots(), &dummies);
    print_three_vertex_method("bare schoonschip", r.schoonschip(), &dummies);
    print_three_vertex_method(
        "expanded bare schoonschip",
        r.expand().schoonschip(),
        &dummies,
    );

    for (name, order) in orders {
        let one_pass = r.schoonschip_with_net::<false, false, AbstractIndex>(
            &SchoonschipSettings::partial()
                .with_expanded_contracted_sums()
                .with_contraction_order(order),
        );
        print_three_vertex_method(
            &format!("net one-pass partial expanded {name}"),
            one_pass.clone(),
            &dummies,
        );
        print_three_vertex_method(
            &format!("net one-pass + smallest cleanup {name}"),
            cleanup_with_smallest_degree(one_pass),
            &dummies,
        );

        let expanded_input = r
            .expand()
            .schoonschip_with_net::<false, false, AbstractIndex>(
                &SchoonschipSettings::partial()
                    .with_expanded_contracted_sums()
                    .with_contraction_order(order),
            );
        print_three_vertex_method(
            &format!("expanded-input net one-pass partial {name}"),
            expanded_input,
            &dummies,
        );

        let expanded_input_full = r
            .expand()
            .schoonschip_with_net::<false, true, AbstractIndex>(
                &SchoonschipSettings::full()
                    .with_expanded_contracted_sums()
                    .with_contraction_order(order),
            );
        print_three_vertex_method(
            &format!("expanded-input net full {name}"),
            expanded_input_full,
            &dummies,
        );

        let full = r.schoonschip_with_net::<false, true, AbstractIndex>(
            &SchoonschipSettings::full()
                .with_expanded_contracted_sums()
                .with_contraction_order(order),
        );
        print_three_vertex_method(&format!("net full expanded {name}"), full.clone(), &dummies);
        print_three_vertex_method(
            &format!("net full + smallest cleanup {name}"),
            cleanup_with_smallest_degree(full),
            &dummies,
        );
    }
}
