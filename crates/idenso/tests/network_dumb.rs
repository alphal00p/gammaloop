use idenso::{
    representations::initialize,
    shorthands::schoonschip::{Schoonschip, SchoonschipContractionOrder, SchoonschipSettings},
};
use spenso::{
    shadowing::{TensorCollectExt, symbolica_utils::SpensoPrintSettings},
    structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
    },
    symbol_set,
};
use symbolica::{
    atom::{Atom, AtomCore},
    parse, symbol,
};

// Generate TestSymbols with all alphabet characters and some multi-character symbols
symbol_set!(TestSymbols, TS;
    mu1 mu2 mu3 mu4 mu5 mu6 mu7 mu8 mu9 mu10 mu11
);

#[test]
fn spenso_bare_symb_vertex_substitution() {
    initialize();
    let _mu1 = TS.mu1;
    let _mink = Minkowski {}.new_rep(4);

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
    let _v6 = parse!(
        "vx(6,-k(4)+k(0), -k(3)+k(4), k(3)-k(0), spenso::mink(4,mu5), spenso::mink(4,mu11), spenso::mink(4,mu6))"
    );
    let _v7 = parse!(
        "vx(7,-k(3)+k(0), -k(2)+k(3), k(2)-k(0), spenso::mink(4,mu6), spenso::mink(4,mu10), spenso::mink(4,mu7))"
    );
    let _v8 = parse!(
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

    let out = r.schoonschip_with_net::<false, AbstractIndex>(
        &SchoonschipSettings::partial()
            .into_single_pass()
            .with_expanded_contracted_sums(),
    );

    println!("out:{}", out.printer(settings));

    for mu in [TS.mu1, TS.mu2, TS.mu3, TS.mu4] {
        assert!(
            out.replace(mu).match_iter().next().is_none(),
            "{}",
            out.printer(settings)
        );
    }
}

fn substituted_three_vertex_reproducer() -> (Atom, [(&'static str, Atom); 3]) {
    initialize();
    let mu1 = TS.mu1;
    let _mink = Minkowski {}.new_rep(4);

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
            ("mu8", TS.mu8.into()),
            ("mu9", TS.mu9.into()),
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
        .into_single_pass()
        .with_expanded_contracted_sums()
        .with_contraction_order(SchoonschipContractionOrder::SmallestDegree);
    for _ in 0..4 {
        let next = result.schoonschip_with_net::<false, AbstractIndex>(&cleanup_settings);
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
fn min_product_terms_three_vertex_simplifies_after_boundary_cleanup() {
    initialize();
    let _ = TS.mu1;
    let (r, dummies) = substituted_three_vertex_reproducer();
    let out = r.schoonschip_with_net::<false, AbstractIndex>(
        &SchoonschipSettings::partial()
            .into_single_pass()
            .with_expanded_contracted_sums()
            .with_contraction_order(SchoonschipContractionOrder::MinProductTerms),
    );

    assert!(residual_dummy_names(&out, &dummies).is_empty());
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
        .into_single_pass()
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
            expr.schoonschip_with_net::<false, AbstractIndex>(&settings),
            &mu1,
            &mu9,
        );
        print_two_dummy_method(
            &format!("{name} / expanded-input network"),
            expr.expand()
                .schoonschip_with_net::<false, AbstractIndex>(&settings),
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
        .into_single_pass()
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
                .schoonschip_with_net::<false, AbstractIndex>(&settings),
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
            &boundary_expression.schoonschip_with_net::<false, AbstractIndex>(&settings),
            &dummies
        )
        .is_empty()
    );
    assert!(
        residual_dummy_names(
            &boundary_expression
                .expand()
                .schoonschip_with_net::<false, AbstractIndex>(&settings),
            &dummies
        )
        .is_empty()
    );
}

#[test]
#[ignore = "diagnostic repro for direct-boundary sum handling"]
fn direct_boundary_cleanup_rewrites_target_sum_before_expanding() {
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

    let boundary_expression = parse!(
        "(a * spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9))
          + b * spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9))
          + c * spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9)))
         * (spenso::g(k(2), spenso::mink(4,mu1))
          * spenso::g(k(3), spenso::mink(4,mu9))
          + spenso::g(k(4), spenso::mink(4,mu1))
          * spenso::g(k(5), spenso::mink(4,mu9))
          + spenso::g(k(6), spenso::mink(4,mu1))
          * spenso::g(k(7), spenso::mink(4,mu9)))"
    );
    let compact_expected = parse!(
        "(a + b + c)
         * (spenso::g(k(2, spenso::mink(4)), k(3, spenso::mink(4)))
          + spenso::g(k(4, spenso::mink(4)), k(5, spenso::mink(4)))
          + spenso::g(k(6, spenso::mink(4)), k(7, spenso::mink(4))))"
    );

    let out = boundary_expression.schoonschip_with_net::<false, AbstractIndex>(&settings);
    let expanded_expected = compact_expected.expand().schoonschip();

    println!(
        "direct-boundary term handling: compact_terms={} expanded_terms={} network_terms={} network_bytes={}",
        compact_expected.nterms(),
        expanded_expected.nterms(),
        out.nterms(),
        out.as_view().get_byte_size()
    );

    assert!(residual_dummy_names(&out, &dummies).is_empty());
    assert_eq!(compact_expected.nterms(), 1);
    assert_eq!(expanded_expected.nterms(), 9);
    assert_eq!(out.nterms(), 3);
    assert!(out.nterms() < expanded_expected.nterms());
}

#[test]
#[ignore = "minimal repro for multi-slot direct-boundary cleanup"]
fn multi_slot_boundary_mwe_needs_target_boundary_expansion() {
    initialize();
    let _mink = Minkowski {}.new_rep(4);

    let (mu1, mu2, mu3) = symbol!("mu1", "mu2", "mu3"; tags=["spenso::index"]);
    symbol!("k"; tags=["spenso::tensor","spenso::rank1"]);

    let dummies = [
        ("mu1", Atom::from(mu1)),
        ("mu2", Atom::from(mu2)),
        ("mu3", Atom::from(mu3)),
    ];
    let settings = SchoonschipSettings::partial()
        .with_expanded_contracted_sums()
        .with_contraction_order(SchoonschipContractionOrder::SmallestDegree);

    let target_after_direct_replacements = parse!(
        "(spenso::g(k(10), spenso::mink(4,mu2))
          + spenso::g(k(11), spenso::mink(4,mu2)))
         * (spenso::g(k(20), spenso::mink(4,mu2))
          + spenso::g(k(21), spenso::mink(4,mu2)))
         * (spenso::g(k(30), k(0, spenso::mink(4)))
          + spenso::g(k(31), k(0, spenso::mink(4))))"
    );
    let compact_target_cleanup = target_after_direct_replacements.schoonschip();
    let expanded_target_cleanup = target_after_direct_replacements.expand().schoonschip();

    print_three_vertex_method(
        "mwe target after direct replacements / compact cleanup",
        compact_target_cleanup.clone(),
        &dummies,
    );
    print_three_vertex_method(
        "mwe target after direct replacements / expanded cleanup",
        expanded_target_cleanup.clone(),
        &dummies,
    );

    let boundary_expression = parse!(
        "spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu2))
         * spenso::g(k(0), spenso::mink(4,mu3))
         * (spenso::g(k(10), spenso::mink(4,mu1))
          + spenso::g(k(11), spenso::mink(4,mu1)))
         * (spenso::g(k(20), spenso::mink(4,mu2))
          + spenso::g(k(21), spenso::mink(4,mu2)))
         * (spenso::g(k(30), spenso::mink(4,mu3))
          + spenso::g(k(31), spenso::mink(4,mu3)))"
    );
    let network_cleanup =
        boundary_expression.schoonschip_with_net::<false, AbstractIndex>(&settings);

    print_three_vertex_method(
        "mwe full boundary / network cleanup",
        network_cleanup.clone(),
        &dummies,
    );

    assert_eq!(
        residual_dummy_names(&compact_target_cleanup, &dummies),
        ["mu2"]
    );
    assert!(residual_dummy_names(&expanded_target_cleanup, &dummies).is_empty());
    assert!(residual_dummy_names(&network_cleanup, &dummies).is_empty());
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
        .into_single_pass()
        .with_expanded_contracted_sums()
        .with_contraction_order(SchoonschipContractionOrder::MinProductTerms);

    let metric_identified_target = parse!(
        "spenso::g(k(2,spenso::mink(4))-k(3,spenso::mink(4)), spenso::mink(4,mu9))
         * spenso::g(k(3,spenso::mink(4))-k(2,spenso::mink(4)), spenso::mink(4,mu9))
         + spenso::g(k(4,spenso::mink(4))-k(5,spenso::mink(4)), spenso::mink(4,mu9))
         * spenso::g(k(5,spenso::mink(4))-k(4,spenso::mink(4)), spenso::mink(4,mu9))"
    );

    let simplified = metric_identified_target
        .normalize_dots()
        .collect_tensors()
        .schoonschip();

    let res = residual_dummy_names(&simplified, &dummies);
    assert!(
        res.is_empty(),
        "residual dummy names: {}, for {}",
        res.join(","),
        simplified
    );

    let boundary_expression = parse!(
        "(spenso::g(spenso::mink(4,mu1), spenso::mink(4,mu9))
           + spenso::g(k(0,spenso::mink(4)), spenso::mink(4,mu1))
             * spenso::g(k(1,spenso::mink(4)), spenso::mink(4,mu9)))
         * (spenso::g(k(2,spenso::mink(4))-k(3,spenso::mink(4)), spenso::mink(4,mu1))
           * spenso::g(k(3,spenso::mink(4))-k(2,spenso::mink(4)), spenso::mink(4,mu9))
           + spenso::g(k(4,spenso::mink(4))-k(5,spenso::mink(4)), spenso::mink(4,mu1))
           * spenso::g(k(5,spenso::mink(4))-k(4,spenso::mink(4)), spenso::mink(4,mu9)))"
    );

    assert!(
        residual_dummy_names(
            &boundary_expression.schoonschip_with_net::<false, AbstractIndex>(&settings),
            &dummies
        )
        .is_empty()
    );
    assert!(
        residual_dummy_names(
            &boundary_expression
                .expand()
                .schoonschip_with_net::<false, AbstractIndex>(&settings),
            &dummies
        )
        .is_empty()
    );
}

#[test]
fn metric_vector_product_with_free_metric_slot_simplifies_in_bare_cleanup() {
    initialize();
    let _mu1 = TS.mu1;
    let _mink = Minkowski {}.new_rep(4);

    symbol!("k"; tags=["spenso::tensor","spenso::rank1"]);

    let mu9: Atom = TS.mu9.into();
    let dummies = [("mu9", mu9.clone())];
    let settings = SchoonschipSettings::partial()
        .into_single_pass()
        .with_expanded_contracted_sums();
    let expr = parse!(
        "spenso::g(spenso::mink(4,mu7), spenso::mink(4,mu9))
         * spenso::g(k(0,spenso::mink(4))-k(1,spenso::mink(4)), spenso::mink(4,mu9))"
    );

    assert!(residual_dummy_names(&expr.collect_metrics().schoonschip(), &dummies).is_empty());
    assert!(
        residual_dummy_names(
            &expr.schoonschip_with_net::<false, AbstractIndex>(&settings),
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
        let one_pass = r.schoonschip_with_net::<false, AbstractIndex>(
            &SchoonschipSettings::partial()
                .into_single_pass()
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

        let expanded_input = r.expand().schoonschip_with_net::<false, AbstractIndex>(
            &SchoonschipSettings::partial()
                .into_single_pass()
                .with_expanded_contracted_sums()
                .with_contraction_order(order),
        );
        print_three_vertex_method(
            &format!("expanded-input net one-pass partial {name}"),
            expanded_input,
            &dummies,
        );

        let expanded_input_full = r.expand().schoonschip_with_net::<false, AbstractIndex>(
            &SchoonschipSettings::full()
                .with_expanded_contracted_sums()
                .with_contraction_order(order),
        );
        print_three_vertex_method(
            &format!("expanded-input net full {name}"),
            expanded_input_full,
            &dummies,
        );

        let full = r.schoonschip_with_net::<false, AbstractIndex>(
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
