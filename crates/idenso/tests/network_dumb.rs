use idenso::{
    IndexTooling,
    representations::initialize,
    shorthands::schoonschip::{Schoonschip, SchoonschipContractionOrder, SchoonschipSettings},
    tensor::SymbolicNetParse,
};
use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
use spenso::{
    network::{
        ExecutionResult, Network, Sequential, SmallestDegree,
        library::{
            panicing::ErroringLibrary,
            symbolic::{ExplicitKey, TensorLibrary},
        },
        parsing::{ParseSettings, ShadowedStructure, StrictTensorFilter},
        store::NetworkStore,
    },
    shadowing::TensorCollectExt,
    structure::{
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        representation::{LibraryRep, Lorentz, Minkowski, RepName},
        slot::{DualSlotTo, DummyAind, IsAbstractSlot},
    },
    symbol_set,
    tensors::data::DataTensor,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function, parse, symbol,
};

// Generate TestSymbols with all alphabet characters and some multi-character symbols
symbol_set!(TestSymbols, TS;
    mu1 mu2 mu3 mu4 mu5 mu6 mu7 mu8 mu9 mu10 mu11
);

fn assert_factorized_contraction(input: &Atom, output: &Atom, settings: &SchoonschipSettings) {
    // Sum-by-sum boundaries deliberately remain factorized in this mode. A
    // metric/rank-one probe under the same settings still requires contraction
    // progress, so returning every input unchanged cannot satisfy this oracle.
    let probe =
        parse!("spenso::g(spenso::mink(4,mu6),spenso::mink(4,mu7))*k(99,spenso::mink(4,mu6))");
    assert_eq!(
        probe.schoonschip_with_net::<false, AbstractIndex>(settings),
        parse!("k(99,spenso::mink(4,mu7))"),
    );

    // Reuse the component-network boundary, independently of the symbolic
    // Schoonschip rewrite. Close external slots with independent probe vectors
    // and compare exact integer assignments after finite component contraction.
    let mut external = input
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap()
        .graph
        .dangling_indices();
    external.sort();
    let probes = external
        .iter()
        .enumerate()
        .fold(Atom::one(), |product, (index, slot)| {
            product * function!(symbol!("k"), 100 + index, slot.to_atom())
        });
    type ComponentLibrary =
        TensorLibrary<DataTensor<Atom, ExplicitKey<AbstractIndex>>, AbstractIndex>;
    let mut library = ComponentLibrary::new();
    library.update_ids();
    let components = [input, output].map(|expression| {
        // Bare k(i) in an existing scalar g(k(i),k(j)) has an implicit
        // Minkowski representation. Make it explicit for the component parser,
        // whose existing shorthand materializer then contracts that dot with
        // the same metric as the indexed vectors, preserving surrounding factors.
        let expression = (expression * &probes).replace_map(|view, _, out| {
            if let AtomView::Fun(vector) = view
                && vector.get_symbol() == symbol!("k")
                && vector.get_nargs() == 1
            {
                **out = function!(
                    symbol!("k"),
                    vector.get(0),
                    function!(symbol!("spenso::mink"), 4)
                );
            }
        });
        let mut network = Network::<
            NetworkStore<DataTensor<Atom, ShadowedStructure<AbstractIndex>>, Atom>,
            ExplicitKey<AbstractIndex>,
            Symbol,
            AbstractIndex,
        >::try_from_view::<ShadowedStructure<AbstractIndex>, ComponentLibrary>(
            expression.as_view(),
            &library,
            &ParseSettings::default().with_strict_tensor_filter(StrictTensorFilter::ContainsReps),
        )
        .unwrap();
        network
            .execute::<
                Sequential,
                SmallestDegree,
                DataTensor<Atom, ExplicitKey<AbstractIndex>>,
                ComponentLibrary,
                ErroringLibrary<Symbol>,
            >(&library, &ErroringLibrary::new())
            .unwrap();
        match network.result_scalar().unwrap() {
            ExecutionResult::One => Atom::one(),
            ExecutionResult::Zero => Atom::zero(),
            ExecutionResult::Val(value) => value.into_owned(),
        }
    });
    let mut nonzero = false;
    for seed in 1..=3 {
        let [before, after] = components.each_ref().map(|expression| {
            expression.replace_map(|view, _, out| {
                let AtomView::Fun(vector) = view else {
                    return;
                };
                if vector.get_symbol() != symbol!("k") {
                    return;
                }
                let Some(AtomView::Fun(component)) = vector.iter().last() else {
                    return;
                };
                if component.get_symbol() != AIND_SYMBOLS.cind {
                    return;
                }
                let id = i64::try_from(vector.get(0)).unwrap();
                let component = i64::try_from(component.get(0)).unwrap();
                **out = Atom::num(
                    ((id + 1) * (component + 2) * (seed + 3) + id * id + component * component) % 7
                        - 3,
                );
            })
        });
        assert!(
            matches!(before.as_view(), AtomView::Num(_)),
            "unassigned tensor components: {before}"
        );
        nonzero |= !before.is_zero();
        assert_eq!(before, after, "component assignment {seed}");
    }
    assert!(nonzero, "the component comparisons must not all vanish");
}

#[test]
fn spenso_bare_symb_vertex_substitution() {
    initialize();
    let _mu1 = TS.mu1;
    let mink = Minkowski {}.new_rep(4);

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
    println!("in:{}", r.printer(settings.clone()));

    let contraction_settings = SchoonschipSettings::partial().into_single_pass();
    let out = r.schoonschip_with_net::<false, AbstractIndex>(&contraction_settings);

    println!("out:{}", out.printer(settings.clone()));

    // The four internal indices remain bound even when a contraction keeps
    // products of tensor sums. Check physical external slots, not whether the
    // bound names disappear from a distributed expression.
    let mut expected = [TS.mu5, TS.mu8, TS.mu9, TS.mu10, TS.mu11]
        .map(|index| mink.slot::<AbstractIndex, _>(index).cast::<LibraryRep>());
    expected.sort();
    for expression in [&r, &out] {
        let net = expression
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                take_first_term_from_sum: false,
                ..Default::default()
            })
            .unwrap();
        let mut external = net.graph.dangling_indices();
        external.sort();
        assert_eq!(external, expected, "{expression}");
    }
    assert_factorized_contraction(&r, &out, &contraction_settings);
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
    let (r, _) = substituted_three_vertex_reproducer();
    let contraction_settings = SchoonschipSettings::partial()
        .into_single_pass()
        .with_contraction_order(SchoonschipContractionOrder::MinProductTerms);
    let out = r.schoonschip_with_net::<false, AbstractIndex>(&contraction_settings);

    // Only mu2 and mu7 are physical external slots. Parsing every sum branch
    // validates the internal contractions without distributing the numerator.
    let mink = Minkowski {}.new_rep(4);
    let mut expected =
        [TS.mu2, TS.mu7].map(|index| mink.slot::<AbstractIndex, _>(index).cast::<LibraryRep>());
    expected.sort();
    for expression in [&r, &out] {
        let net = expression
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                take_first_term_from_sum: false,
                ..Default::default()
            })
            .unwrap();
        let mut external = net.graph.dangling_indices();
        external.sort();
        assert_eq!(external, expected, "{expression}");
    }
    assert_factorized_contraction(&r, &out, &contraction_settings);
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

    for (name, order) in orders {
        let one_pass = r.schoonschip_with_net::<false, AbstractIndex>(
            &SchoonschipSettings::partial()
                .into_single_pass()
                .with_contraction_order(order),
        );
        print_three_vertex_method(
            &format!("net one-pass partial factorized {name}"),
            one_pass.clone(),
            &dummies,
        );
        print_three_vertex_method(
            &format!("net one-pass + smallest cleanup {name}"),
            cleanup_with_smallest_degree(one_pass),
            &dummies,
        );

        let full = r.schoonschip_with_net::<false, AbstractIndex>(
            &SchoonschipSettings::full().with_contraction_order(order),
        );
        print_three_vertex_method(
            &format!("net full factorized {name}"),
            full.clone(),
            &dummies,
        );
        print_three_vertex_method(
            &format!("net full + smallest cleanup {name}"),
            cleanup_with_smallest_degree(full),
            &dummies,
        );
    }
}

#[test]
fn canonicalization_keeps_explicit_dummies_distinct_from_compact_dot_indices() {
    initialize();
    let mink = Minkowski {}.new_rep(4);
    let vector = symbol!("capture_probe");
    let spectator = parse!("(capture_a+capture_b)*(capture_c+capture_d)");
    for dots in 1..=2 {
        let compact = (0..dots).fold(Atom::one(), |product, index| {
            product
                * function!(
                    spenso::network::tags::SPENSO_TAG.dot,
                    function!(vector, 3 + 2 * index, mink.to_symbolic([])),
                    function!(vector, 4 + 2 * index, mink.to_symbolic([]))
                )
        });
        // The control differs only by the explicit dummy's name. Both forms
        // contain three or fewer independent closed Lorentz contractions.
        let mut expected = None;
        for reserved in [17, 1_000_000, 1_000_001] {
            let slot = mink.slot::<AbstractIndex, _>(AbstractIndex::new_dummy_at(reserved));
            let input = &spectator
                * &compact
                * function!(vector, 1, slot.to_atom())
                * function!(vector, 2, slot.to_atom());
            let actual = input.canonize(AbstractIndex::Dummy);
            if let Some(expected) = &expected {
                assert_eq!(
                    &actual, expected,
                    "reserved dummy {reserved}, {dots} compact dots"
                );
            } else {
                expected = Some(actual.clone());
            }
            for factor in [parse!("capture_a+capture_b"), parse!("capture_c+capture_d")] {
                assert!(
                    actual
                        .pattern_match(&factor.to_pattern(), None, None)
                        .next()
                        .is_some()
                );
            }
        }
    }
}

#[test]
fn canonicalization_keeps_external_dummy_names_distinct_from_canonical_dummies() {
    initialize();
    let mink = Minkowski {}.new_rep(4);
    for reserved in [0, 1] {
        let external = mink
            .slot::<AbstractIndex, _>(AbstractIndex::new_dummy_at(reserved))
            .to_atom();
        let expression = |dummy| {
            let slot = mink
                .slot::<AbstractIndex, _>(AbstractIndex::new_dummy_at(dummy))
                .to_atom();
            (function!(symbol!("capture_open_a"), slot.clone())
                + function!(symbol!("capture_open_b"), slot.clone()))
                * (function!(symbol!("capture_open_c"), slot.clone(), external.clone())
                    + function!(symbol!("capture_open_d"), slot, external.clone()))
        };
        let canonical = expression(17).canonize(AbstractIndex::Dummy);
        assert_eq!(canonical, expression(18).canonize(AbstractIndex::Dummy));
        assert_eq!(canonical, canonical.canonize(AbstractIndex::Dummy));
        assert_eq!(canonical.list_dangling::<AbstractIndex>(), vec![external]);
        assert!(matches!(canonical.as_view(), AtomView::Mul(product)
            if product.iter().filter(|factor| matches!(factor, AtomView::Add(_))).count() == 2));
    }
}

#[test]
fn canonicalization_keeps_dual_external_slots_distinct_from_canonical_dummies() {
    initialize();
    let rep = Lorentz {}.new_rep(4);
    let external = rep
        .slot::<AbstractIndex, _>(AbstractIndex::Dummy(0))
        .dual()
        .to_atom();
    let spectator = parse!("(capture_a+capture_b)*(capture_c+capture_d)");
    let expression = |dummy| {
        let slot = rep.slot::<AbstractIndex, _>(AbstractIndex::Dummy(dummy));
        &spectator
            * function!(symbol!("capture_dual_t"), slot.to_atom(), external.clone())
            * function!(symbol!("capture_dual_u"), slot.dual().to_atom())
    };
    let canonical = expression(17).canonize(AbstractIndex::Dummy);
    assert_eq!(canonical, expression(18).canonize(AbstractIndex::Dummy));
    assert_eq!(canonical, canonical.canonize(AbstractIndex::Dummy));
    assert_eq!(canonical.list_dangling::<AbstractIndex>(), vec![external]);
    assert!(
        canonical
            .pattern_match(&spectator.to_pattern(), None, None)
            .next()
            .is_some()
    );
}

#[test]
fn canonicalization_ignores_indices_from_canceled_representation_groups() {
    initialize();
    let expressions = [17, 18].map(|dummy| {
        let index = AbstractIndex::Dummy(dummy);
        let mink = Minkowski {}.new_rep(4).slot::<AbstractIndex, _>(index);
        let lorentz = Lorentz {}.new_rep(4).slot::<AbstractIndex, _>(index);
        [
            (mink.to_atom(), mink.to_atom()),
            (lorentz.to_atom(), lorentz.dual().to_atom()),
        ]
        .map(|(slot, dual)| {
            function!(symbol!("capture_cancel_a"), slot)
                * function!(symbol!("capture_cancel_b"), dual)
        })
    });
    // Exercise both representation orderings: canceled terms must not consume
    // a canonical dummy number ahead of the surviving contraction.
    for canceled in 0..2 {
        let survivor = &expressions[0][1 - canceled];
        let input = &expressions[0][canceled] - &expressions[1][canceled] + survivor;
        let canonical = input.canonize(AbstractIndex::Dummy);
        assert_eq!(canonical, survivor.canonize(AbstractIndex::Dummy));
        assert_eq!(canonical, canonical.canonize(AbstractIndex::Dummy));
        assert!(canonical.list_dangling::<AbstractIndex>().is_empty());
    }
}

#[test]
fn tensor_canonicalization_closes_contractions_inside_nested_sums() {
    initialize();
    let input = parse!(
        "nested_spectator*((nested_a(nested_mu)+nested_b(nested_mu))
            *(nested_c(nested_mu)+nested_d(nested_mu))
          +(nested_e(nested_mu)+nested_f(nested_mu))
            *(nested_g(nested_mu)+nested_h(nested_mu)))"
    );
    // Each alternative closes the same index inside two vector sums. Keep
    // both products of sums intact while validating the canonicalizer's scope.
    let canonical = input.canonize_tensors([(parse!("nested_mu"), 0)]).unwrap();
    assert!(canonical.external_indices.is_empty());
    assert_eq!(canonical.canonical_form, input);
}
