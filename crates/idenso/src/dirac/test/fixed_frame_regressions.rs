//! The OS-branch fixed-frame algebra inputs, using Feynkit's typed entry points.
//! Scalar spectators deliberately remain factored in inputs and certificates.
use super::*;

#[track_caller]
fn simplify(expression: &Atom, settings: GammaSimplifySettings) -> Atom {
    SymbolicTensor::infer(expression.clone())
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            gamma: Some(settings),
            epsilon: settings.output == crate::dirac::GammaOutput::Reduced,
            ..Default::default()
        })
        .unwrap()
        .into_expression()
}

fn collect(expression: &Atom) -> Atom {
    SymbolicTensor::infer(expression.clone())
        .unwrap()
        .contract(crate::tensor::ContractSettings::default())
        .unwrap()
        .into_expression()
}

#[test]
fn mixed_temporal_traces_reduce_without_fixing_formal_dimensions() {
    use spenso::{
        network::library::symbolic::ETS, shadowing::IntoAtom,
        structure::abstract_index::AIND_SYMBOLS,
    };
    use symbolica::function;

    let r = test_initialize();
    let mu = slot!(r.mink4, mu).into_atom();
    let nu = slot!(r.mink4, nu).into_atom();
    let time_mu = function!(ETS.delta, function!(AIND_SYMBOLS.cind, 0), mu.clone());
    let time_nu = function!(ETS.delta, function!(AIND_SYMBOLS.cind, 0), nu.clone());
    let expression = trace!(
        r.bis4.to_symbolic([]),
        gamma0!(),
        gamma!(mu.clone()),
        gamma0!(),
        gamma!(nu.clone()),
    );
    let expected = Atom::num(8) * time_mu * time_nu - Atom::num(4) * g!(mu, nu);
    let reduced = simplify(&expression, GammaSimplifySettings::default());
    assert!(!reduced.contains_symbol(spenso::network::tags::SPENSO_TAG.trace));
    assert_eq!(reduced, expected);
    assert_eq!(
        simplify(&reduced, GammaSimplifySettings::default()),
        reduced
    );

    let odd = trace!(
        r.bis4.to_symbolic([]),
        gamma0!(),
        gamma!(slot!(r.mink4, mu)),
        gamma!(slot!(r.mink4, nu)),
    );
    assert!(simplify(&odd, GammaSimplifySettings::default()).is_zero());

    for (spin, vector) in [(r.bis_d, r.mink4), (r.bis4, r.mink_d)] {
        let formal = trace!(
            spin.to_symbolic([]),
            gamma0!(),
            gamma!(slot!(vector, mu)),
            gamma0!(),
            gamma!(slot!(vector, nu)),
        );
        assert!(
            simplify(&formal, GammaSimplifySettings::default())
                .contains_symbol(spenso::network::tags::SPENSO_TAG.trace)
        );
    }
}

#[test]
fn fixed_temporal_slash_uses_gamma_zero_without_fixing_lorentz_dimension() {
    use spenso::{
        network::library::symbolic::ETS, shadowing::IntoAtom,
        structure::abstract_index::AIND_SYMBOLS,
    };
    use symbolica::{function, symbol};

    let r = test_initialize();
    let dimension = Atom::var(symbol!("temporal_slash_dimension"));
    let vector_slot = Minkowski {}.to_symbolic([&dimension]);
    let temporal = function!(ETS.delta, function!(AIND_SYMBOLS.cind, 0), &vector_slot);
    let slash = gamma!(temporal);
    let start = slot!(r.bis4, a).into_atom();
    let end = slot!(r.bis4, b).into_atom();
    let middle = slot!(r.bis4, c).into_atom();
    let evanescent = &dimension - Atom::num(4);
    for settings in [
        GammaSimplifySettings::canonical(),
        GammaSimplifySettings::repeated_pairs(),
    ] {
        for factors in [
            vec![gamma0!(), slash.clone()],
            vec![slash.clone(), gamma0!()],
            vec![slash.clone(), slash.clone()],
        ] {
            let expression = &evanescent * chain!(&start, &end; factors);
            assert_eq!(
                simplify(&expression, settings),
                &evanescent * g!(&start, &end),
            );
        }
        let projected =
            &evanescent * chain!(&middle, &start, gamma0!()) * chain!(&middle, &end, slash.clone());
        assert_eq!(
            simplify(&projected, settings),
            &evanescent * g!(&start, &end),
        );
    }
    let spatial = gamma!(function!(
        ETS.delta,
        function!(AIND_SYMBOLS.cind, 1),
        &vector_slot
    ));
    assert!(
        simplify(
            &chain!(&start, &end, spatial),
            GammaSimplifySettings::default()
        )
        .contains_symbol(AGS.gamma)
    );
}

#[test]
fn fixed_vector_sums_remain_admissible_after_formal_dirac_reduction() {
    use crate::tensor::{SymbolicNetExt, SymbolicNetParse};
    use spenso::{
        network::{
            library::symbolic::ETS,
            parsing::{ParseSettings, SchoonschipExpansionMode, ShorthandParsing},
            tags::SPENSO_TAG,
        },
        structure::abstract_index::{AIND_SYMBOLS, AbstractIndex},
    };
    use symbolica::function;

    let r = test_initialize();
    let a = r.bis4.to_symbolic([Atom::num(98201)]);
    let b = r.bis4.to_symbolic([Atom::num(98202)]);
    let c = r.bis4.to_symbolic([Atom::num(98203)]);
    let e = r.bis4.to_symbolic([Atom::num(98204)]);
    let q = SPENSO_TAG.rank_one_tensor_symbol("fixed_vector_sum::q");
    let mass = Atom::var(symbol!("fixed_vector_sum::M"));
    let spectator = function!(
        symbol!("fixed_vector_sum::opaque"; Scalar),
        (Atom::var(symbol!("fixed_vector_sum::x")) + Atom::var(symbol!("fixed_vector_sum::y")))
            .pow(11)
    );
    let settings = crate::tensor::AlgebraSettings {
        gamma: Some(GammaSimplifySettings::canonical()),
        ..Default::default()
    };
    let expose = |value: &Atom| {
        value
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Expand {
                    schoonschip: SchoonschipExpansionMode {
                        inner_products: true,
                        expand_inside_chains: true,
                        expand_schoonship: true,
                    },
                    trace: false,
                    chain: false,
                },
                ..Default::default()
            })
            .unwrap()
            .simple_execute::<()>()
            .unwrap()
    };
    for dimension in [Atom::var(symbol!("fixed_vector_sum::D")), Atom::num(4)] {
        let compact = Minkowski {}.to_symbolic([&dimension]);
        let mu = Minkowski {}.to_symbolic([&dimension, &Atom::num(98205)]);
        let nu = Minkowski {}.to_symbolic([&dimension, &Atom::num(98206)]);
        let fixed = |i, port: &Atom| function!(ETS.delta, function!(AIND_SYMBOLS.cind, i), port);
        let momentum = function!(q, &compact);
        let spatial = (1..4).map(|i| fixed(i, &nu)).sum::<Atom>();
        let velocity = (0..4).map(|i| fixed(i, &compact)).sum::<Atom>();
        for component in 0..4 {
            let selected = function!(
                ETS.delta,
                function!(AIND_SYMBOLS.cind, component),
                &mass * &velocity
            );
            let selected = SymbolicTensor::infer(selected).unwrap().to_dots().unwrap();
            let sign = Atom::num(if component == 0 { 1 } else { -1 });
            assert_eq!(expose(selected.expression()), &sign * &mass);
            let weighted = function!(
                ETS.delta,
                function!(AIND_SYMBOLS.cind, component),
                &spectator * &mass * &velocity
            );
            let weighted = SymbolicTensor::infer(weighted).unwrap().to_dots().unwrap();
            assert!(weighted.expression().contains(spectator.as_view()));
            assert_eq!(expose(weighted.expression()), &sign * &mass * &spectator);
        }
        let prefix = spatial * gamma!(&a, &c, &nu) + function!(AGS.gamma0, &c, &a);
        let word = chain!(
            &c,
            &e,
            gamma!(&momentum),
            gamma!(&mu),
            gamma!(&momentum),
            gamma!(&(&mass * &velocity))
        );
        let source =
            SymbolicTensor::infer(&spectator * prefix * gamma!(&e, &b, &mu) * word).unwrap();
        assert_eq!(source.rank(), 2);
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(result.structure(), source.structure());
        assert!(result.expression().contains(spectator.as_view()));
        assert_eq!(
            SymbolicTensor::infer(result.expression().clone())
                .unwrap()
                .structure(),
            source.structure()
        );

        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
        if dimension != Atom::num(4) {
            let evanescent = (&dimension - Atom::num(4)) * g!(&a, &b);
            assert_eq!(
                SymbolicTensor::infer(evanescent.clone())
                    .unwrap()
                    .simplify_algebra(&settings)
                    .unwrap()
                    .expression(),
                &evanescent
            );
        }
        let other = Minkowski {}.to_symbolic([Atom::num(7)]);
        let mismatched = function!(
            ETS.delta,
            function!(AIND_SYMBOLS.cind, 0),
            fixed(0, &compact) + fixed(0, &other)
        );
        assert!(SymbolicTensor::infer(mismatched).is_err());
        let invalid = function!(ETS.delta, function!(AIND_SYMBOLS.cind, 0), g!(&mu, &nu));
        assert!(
            SymbolicTensor::infer(invalid).is_err(),
            "a rank-two operand is not a compact vector"
        );
    }
}

#[test]
fn fixed_spatial_slashes_anticommute_with_gamma_zero_in_symbolic_dimension() {
    use spenso::{
        network::library::symbolic::ETS, shadowing::IntoAtom,
        structure::abstract_index::AIND_SYMBOLS,
    };
    use symbolica::{function, symbol};
    let r = test_initialize();
    let dimension = Atom::var(symbol!("fixed_spatial_slash::d"));
    let vector_slot = Minkowski {}.to_symbolic([&dimension]);
    let start = slot!(r.bis4, a).into_atom();
    let end = slot!(r.bis4, b).into_atom();
    let spectator = (&dimension - Atom::num(2))
        * (Atom::var(symbol!("fixed_spatial_slash::x"))
            + Atom::var(symbol!("fixed_spatial_slash::y")))
        .pow(7);
    for settings in [
        GammaSimplifySettings::canonical(),
        GammaSimplifySettings::repeated_pairs(),
    ] {
        for component in 1..4 {
            let spatial = gamma!(function!(
                ETS.delta,
                function!(AIND_SYMBOLS.cind, component),
                &vector_slot
            ));
            let temporal_first = chain!(&start, &end, gamma0!(), spatial.clone());
            let temporal_last = chain!(&start, &end, spatial.clone(), gamma0!());
            let canonical = simplify(&(&spectator * &temporal_first), settings);
            assert_eq!(
                canonical,
                -simplify(&(&spectator * &temporal_last), settings)
            );
            assert_eq!(simplify(&canonical, settings), canonical);
            assert!(simplify(&(temporal_first + temporal_last), settings).is_zero());
            assert_eq!(
                simplify(
                    &(&spectator * chain!(&start, &end, spatial.clone(), spatial.clone())),
                    settings
                ),
                -&spectator * g!(&start, &end)
            );
            let conjugated = simplify(
                &(&spectator * chain!(&start, &end, gamma0!(), spatial.clone(), gamma0!())),
                settings,
            );
            assert_eq!(
                conjugated,
                simplify(&(-&spectator * chain!(&start, &end, spatial)), settings)
            );
        }
        let momentum = function!(
            spenso::network::tags::SPENSO_TAG.rank_one_tensor_symbol("fixed_spatial_slash::p"),
            &vector_slot
        );
        let generic = chain!(&start, &end, gamma0!(), gamma!(momentum.clone()))
            + chain!(&start, &end, gamma!(momentum), gamma0!());
        assert!(
            !simplify(&generic, settings).is_zero(),
            "a general momentum slash has a nonzero temporal anticommutator"
        );
        let mu =
            Minkowski {}.to_symbolic([&dimension, &Atom::var(symbol!("fixed_spatial_slash::mu"))]);
        let unrestricted_trace = trace!(r.bis4.to_symbolic([]), gamma!(mu.clone()), gamma!(mu));
        assert_eq!(
            simplify(&unrestricted_trace, settings),
            Atom::num(4) * &dimension
        );
    }
}

#[test]
fn gamma_zero_chain_composition_preserves_order_for_small_and_large_dummy_labels() {
    use spenso::shadowing::IntoAtom;
    use symbolica::{function, symbol};
    let r = test_initialize();
    let mu = slot!(r.mink_d, mu).into_atom();
    let nu = slot!(r.mink_d, nu).into_atom();
    let tags = &*spenso::network::tags::SPENSO_TAG;
    let transposed_mu = function!(AGS.gamma, tags.chain_out, tags.chain_in, &mu);
    let transposed_nu = function!(AGS.gamma, tags.chain_out, tags.chain_in, &nu);
    let spectator = (Atom::var(symbol!("gamma_zero_chain_composition::x"))
        + Atom::var(symbol!("gamma_zero_chain_composition::y")))
    .pow(11);
    for dummy_label in [0, 1000000] {
        // Typed inference admits atomic/scoped indices, not the old arbitrary
        // index(label) metadata. Keep the same label ordering with typed indices.
        let [left, right, middle] = [1, 2, dummy_label].map(|label| {
            r.bis4
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(label))
                .into_atom()
        });
        for settings in [
            GammaSimplifySettings::canonical(),
            GammaSimplifySettings::repeated_pairs(),
        ] {
            for (source, expected, reversed) in [
                (
                    gamma0!(&left, &middle) * gamma!(&middle, &right, &mu),
                    chain!(&left, &right, gamma0!(), gamma!(&mu)),
                    chain!(&right, &left, &transposed_mu, gamma0!()),
                ),
                (
                    gamma!(&left, &middle, &mu) * gamma0!(&middle, &right),
                    chain!(&left, &right, gamma!(&mu), gamma0!()),
                    chain!(&right, &left, gamma0!(), &transposed_mu),
                ),
            ] {
                let actual = simplify(&(&spectator * source), settings);
                // The self-dual collector may choose either traversal. These
                // are exactly equivalent by (AB)^T = B^T A^T and gamma0^T = gamma0,
                // with the external endpoints exchanged in the reversed form.
                assert!(
                    actual == simplify(&(&spectator * expected), settings)
                        || actual == simplify(&(&spectator * reversed), settings),
                    "unexpected chain orientation: {actual}"
                );
                assert_eq!(simplify(&actual, settings), actual);
                let collected = collect(&actual);
                assert_eq!(collect(&collected), collected);
            }
            let arbitrary_same_input = gamma!(&middle, &left, &mu) * gamma!(&middle, &right, &nu);
            let actual = simplify(&arbitrary_same_input, settings);
            for ordinary_product in [
                chain!(&left, &right, gamma!(&mu), gamma!(&nu)),
                chain!(&right, &left, &transposed_nu, &transposed_mu),
            ] {
                assert_ne!(
                    actual,
                    simplify(&ordinary_product, settings),
                    "ordinary gamma matrices must not acquire a transpose symmetry"
                );
            }
            // Distinct chain lengths keep this closed contraction explicit,
            // rather than encoding a tensor square as an opaque scalar power.
            let closed = chain!(&left, &middle, gamma0!())
                * chain!(&left, &middle, gamma0!(), gamma0!(), gamma0!());
            assert_eq!(simplify(&closed, settings), Atom::num(4));
        }
    }
}

#[test]
fn gamma_zero_symbolic_spin_does_not_use_four_dimensional_identity() {
    use spenso::shadowing::IntoAtom;
    let r = test_initialize();
    let mu = slot!(r.mink_d, mu).into_atom();
    // The special symmetry join is confined to the existing four-dimensional
    // spinor convention, even when the Lorentz dimension remains symbolic.
    let a = slot!(r.bis_d, a).into_atom();
    let b = slot!(r.bis_d, b).into_atom();
    let c = slot!(r.bis_d, c).into_atom();
    let generic_spin = chain!(&a, &b, gamma0!()) * chain!(&a, &c, gamma!(&mu));
    // The new collector represents common-end joins by an explicit transpose.
    // Keep that normal-form change, but require the symbolic-spin word to stay
    // unreduced: the four-dimensional gamma0 symmetry cannot be used here.
    let collected = collect(&generic_spin);
    assert_eq!(
        simplify(&generic_spin, GammaSimplifySettings::canonical()),
        collected
    );
}

#[test]
fn gamma_zero_common_end_joins_preserve_physical_transpose_symmetry() {
    use spenso::shadowing::IntoAtom;
    let r = test_initialize();
    let left = slot!(r.bis4, gamma_zero_common_left).into_atom();
    let right = slot!(r.bis4, gamma_zero_common_right).into_atom();
    let middle = slot!(r.bis4, gamma_zero_common_middle).into_atom();
    let mu = slot!(r.mink_d, gamma_zero_common_mu).into_atom();
    let tags = &*spenso::network::tags::SPENSO_TAG;
    let transposed_mu = symbolica::function!(AGS.gamma, tags.chain_out, tags.chain_in, &mu);
    for settings in [
        GammaSimplifySettings::canonical(),
        GammaSimplifySettings::repeated_pairs(),
    ] {
        for (source, expected, reversed) in [
            (
                gamma0!(&middle, &left) * gamma!(&middle, &right, &mu),
                chain!(&left, &right, gamma0!(), gamma!(&mu)),
                chain!(&right, &left, &transposed_mu, gamma0!()),
            ),
            (
                gamma!(&left, &middle, &mu) * gamma0!(&right, &middle),
                chain!(&left, &right, gamma!(&mu), gamma0!()),
                chain!(&right, &left, gamma0!(), &transposed_mu),
            ),
        ] {
            let actual = simplify(&source, settings);
            assert!(
                actual == simplify(&expected, settings) || actual == simplify(&reversed, settings),
                "unexpected chain orientation: {actual}"
            );
        }
    }
}

#[test]
fn fixed_minkowski_basis_products_contract_through_typed_pipeline() {
    use spenso::{dot, network::library::symbolic::ETS, structure::abstract_index::AIND_SYMBOLS};
    use symbolica::{function, symbol};
    test_initialize();
    let dimension = Atom::var(symbol!("fixed_basis_pipeline::d"));
    let vector_slot = Minkowski {}.to_symbolic([&dimension]);
    let basis = |component| {
        function!(
            ETS.delta,
            function!(AIND_SYMBOLS.cind, component),
            &vector_slot
        )
    };
    let spectator = (&dimension - Atom::num(2))
        * (Atom::var(symbol!("fixed_basis_pipeline::x"))
            + Atom::var(symbol!("fixed_basis_pipeline::y")))
        .pow(9);
    for left in 0..4 {
        for right in 0..4 {
            let sign = if left != right {
                0
            } else if left == 0 {
                1
            } else {
                -1
            };
            for product in [
                g!(basis(left), basis(right)),
                dot!(basis(left), basis(right)),
            ] {
                let input = &spectator * product;
                let actual = collect(&input);
                assert_eq!(actual, Atom::num(sign) * &spectator);
                assert_eq!(collect(&actual), actual);
            }
        }
    }
}
