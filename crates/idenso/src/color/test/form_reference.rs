use super::*;
use crate::{
    coad,
    shorthands::{UndoShorthands, chain::Chain},
};
use insta::assert_snapshot;
use spenso::{antisym, g, mink, shadowing::IntoAtom, structure::representation::RepName};
use symbolica::{function, symbol};

#[test]
fn closed_quartic_invariants_with_two_structure_constants() {
    test_initialize();
    let source = parse_lit!(
        d(
            cof(Nc),
            coad(Nc ^ 2 - 1, a),
            coad(Nc ^ 2 - 1, b),
            coad(Nc ^ 2 - 1, c),
            coad(Nc ^ 2 - 1, d)
        ) * d(
            coad(Nc ^ 2 - 1),
            coad(Nc ^ 2 - 1, e),
            coad(Nc ^ 2 - 1, f),
            coad(Nc ^ 2 - 1, c),
            coad(Nc ^ 2 - 1, d)
        ) * f(
            coad(Nc ^ 2 - 1, x),
            coad(Nc ^ 2 - 1, a),
            coad(Nc ^ 2 - 1, e)
        ) * f(
            coad(Nc ^ 2 - 1, x),
            coad(Nc ^ 2 - 1, b),
            coad(Nc ^ 2 - 1, f)
        ),
        default_namespace = "spenso"
    );
    // Invariance on either symmetric rank-four tensor gives CA/3 times
    // their Gram contraction. Its independent exact SU(3) matrix value is 15/2.
    let n = Atom::var(CS.nc);
    let expected =
        n.clone().pow(2) * (n.clone().pow(2) - Atom::num(1)) * (n.clone().pow(2) + Atom::num(6))
            / Atom::num(144);
    let settings = ColorSimplifySettings::default().with_cof_dimension_invariants();
    let result = simplify_formal_color(&source, settings);
    assert_eq!(result.expand(), expected.expand());
    assert_eq!(simplify_formal_color(&result, settings), result);
    assert_eq!(
        result
            .replace(n.to_pattern())
            .with(Atom::num(3).to_pattern()),
        Atom::num(15) / Atom::num(2)
    );
}

#[test]
fn seven_generator_trace_with_three_structure_constants() {
    test_initialize();
    let n = Atom::var(CS.nc);
    let adjoint_dimension = n.clone().pow(2) - Atom::num(1);
    let [a, b, c, d, e, f, g, h] = ["a", "b", "c", "d", "e", "f", "g", "h"].map(|label| {
        ColorAdjoint {}.to_symbolic([
            adjoint_dimension.clone(),
            Atom::var(symbol!(&format!("quartic_capture::{label}"))),
        ])
    });
    let fundamental = ColorFundamental {}.to_symbolic([n.clone()]);
    let trace = trace!(&fundamental; [&a, &b, &d, &g, &a, &e, &f].map(|a| color_t!(a)));
    let source =
        Atom::i() * trace * color_f!(&b, &e, &h) * color_f!(&f, &c, &d) * color_f!(&h, &c, &g);
    // Expand compact projector words before signed canonicalization. This
    // removes the zero odd-invariant remainder using the existing tensor
    // canonicalizer, without a separate color identity for that remainder.
    let cooking = crate::CookSettings::reversible().with_representation_payloads(true, true);
    let settings = crate::tensor::AlgebraSettings {
        color: Some(
            ColorSimplifySettings::default()
                .without_one_shot_traces()
                .with_cof_dimension_invariants(),
        ),
        ..Default::default()
    };
    let result = SymbolicTensor::infer(cooking.try_cook(source.as_view()).unwrap())
        .unwrap()
        .simplify_algebra(&settings)
        .unwrap()
        .contract(Default::default())
        .unwrap()
        .expand_projectors()
        .unwrap()
        .undo_trace()
        .unwrap()
        .undo_chain()
        .unwrap()
        .expanded(None, false)
        .unwrap()
        .canonize()
        .unwrap()
        .simplify_algebra(&settings)
        .unwrap()
        .into_expression();
    let result = cooking.uncook(result.as_view());
    let expected =
        -(n.clone().pow(2) - Atom::num(1)) * (n.clone().pow(2) - Atom::num(2)) / Atom::num(32);
    assert_eq!(result.expand(), expected.expand());
    let fixed = result
        .replace(n.to_pattern())
        .with(Atom::num(3).to_pattern());
    assert_eq!(fixed, -Atom::num(7) / Atom::num(4), "{result}");
}

macro_rules! fco {
    ($r:ident, $a:tt, $b:tt, $c:tt) => {
        color_f!(
            slot!($r.coad_da, $a),
            slot!($r.coad_da, $b),
            slot!($r.coad_da, $c),
        )
    };
}

macro_rules! tco {
    ($r:ident, $a:tt) => {
        color_t!(slot!($r.coad_da, $a))
    };
}

#[test]
fn one_generator_trace_vanishes() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, i))),
        default_namespace = "spenso"
    );

    // Strict representation dimensions are atomic. Restore the FORM dimension
    // after the typed operation, before comparing with the original reference.
    let dimension = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let atomic_dimension = TestReps::new().coad_da.dim.to_symbolic();
    let admitted = expr
        .replace(dimension.to_pattern())
        .with(atomic_dimension.to_pattern());
    assert_eq!(
        admitted
            .replace(atomic_dimension.to_pattern())
            .with(dimension.to_pattern()),
        expr
    );
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression()
        .replace(atomic_dimension.to_pattern())
        .with(dimension.to_pattern());

    assert!(result.is_zero());
}

#[test]
fn two_generator_trace_normalizes_to_tr_metric() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, i))),
        default_namespace = "spenso"
    );

    // Strict representation dimensions are atomic. Restore the FORM dimension
    // after the typed operation, before comparing with the original reference.
    let dimension = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let atomic_dimension = TestReps::new().coad_da.dim.to_symbolic();
    let admitted = expr
        .replace(dimension.to_pattern())
        .with(atomic_dimension.to_pattern());
    assert_eq!(
        admitted
            .replace(atomic_dimension.to_pattern())
            .with(dimension.to_pattern()),
        expr
    );
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression()
        .replace(atomic_dimension.to_pattern())
        .with(dimension.to_pattern());

    assert_snapshot!(result.to_bare_ordered_string(), @"g(coad(-1+Nc^2,a),coad(-1+Nc^2,b))*idx(2,cof(Nc))");
}

#[test]
fn invariant_open_fundamental_chain_reduces_with_symbolic_nc() {
    test_initialize();
    let n = Atom::var(CS.nc);
    let adjoints = [
        parse_lit!(coad(Nc ^ 2 - 1, a), default_namespace = "spenso"),
        parse_lit!(coad(Nc ^ 2 - 1, b), default_namespace = "spenso"),
        parse_lit!(coad(Nc ^ 2 - 1, c), default_namespace = "spenso"),
    ];
    let factors = adjoints
        .into_iter()
        .map(|adjoint| color_t!(adjoint))
        .collect::<Vec<_>>();
    let closed = trace!(ColorFundamental {}.to_symbolic([n.clone()]); factors.clone());
    let start = parse_lit!(cof(Nc, i), default_namespace = "spenso");
    let end = parse_lit!(dind(cof(Nc, j)), default_namespace = "spenso");
    let metric = function!(spenso::network::library::symbolic::ETS.metric, &start, &end);
    let spectator = parse_lit!((opaque(x) + opaque(y)) ^ 5 * (opaque(z) + opaque(w)));
    let settings = ColorSimplifySettings::default().with_cof_dimension_invariants();
    for reverse in [false, true] {
        let mut word = factors.clone();
        if reverse {
            word.reverse();
        }
        let open = chain!(&start, &end; word);
        // Independent Fierz contraction gives the scalar operator; the two
        // original free fundamental slots and spectator remain factorized.
        let coefficient = if reverse {
            (n.clone().pow(2) - Atom::num(1)) * (n.clone().pow(2) - Atom::num(2))
                / (Atom::num(8) * n.clone().pow(2))
        } else {
            -(n.clone().pow(2) - Atom::num(1)) / (Atom::num(4) * n.clone().pow(2))
        };
        let simplified = simplify_formal_color(&(&closed * &open * &spectator), settings);
        let scalar = simplified.coefficient((&metric * &spectator).as_view());
        assert_eq!(scalar.expand(), coefficient.expand(), "{simplified}");
        assert_eq!(simplified, scalar * &metric * &spectator);
        assert_eq!(simplify_formal_color(&simplified, settings), simplified);
        assert_eq!(
            crate::CookSettings::reversible()
                .with_representation_payloads(true, true)
                .cook(simplified.as_view())
                .list_dangling::<AbstractIndex>()
                .expect("formal reference slots"),
            crate::CookSettings::reversible()
                .with_representation_payloads(true, true)
                .cook(metric.as_view())
                .list_dangling::<AbstractIndex>()
                .expect("formal reference slots"),
        );
    }
}

#[test]
fn repeated_generators_in_long_fundamental_trace_use_fierz() {
    test_initialize();
    let n = Atom::var(CS.nc);
    let factors = [
        parse_lit!(coad(Nc ^ 2 - 1, a), default_namespace = "spenso"),
        parse_lit!(coad(Nc ^ 2 - 1, b), default_namespace = "spenso"),
        parse_lit!(coad(Nc ^ 2 - 1, c), default_namespace = "spenso"),
    ]
    .into_iter()
    .map(|adjoint| color_t!(adjoint))
    .collect::<Vec<_>>();
    let word = factors.iter().chain(&factors).cloned().collect::<Vec<_>>();
    let trace = trace!(ColorFundamental {}.to_symbolic([n.clone()]); word);
    // Sum_abc Tr(Ta Tb Tc Ta Tb Tc) = (Nc^2 - Nc^-2)/8.
    let expected = (n.clone().pow(2) - n.pow(-2)) / Atom::num(8);
    let settings = ColorSimplifySettings::default().with_cof_dimension_invariants();
    assert_eq!(
        simplify_formal_color(&trace, settings).expand(),
        expected.expand()
    );
    assert_eq!(
        simplify_formal_color(&trace, settings.without_trace_evaluation()),
        trace,
    );
}

#[test]
fn equivariant_open_fundamental_chain_keeps_its_adjoint_generator() {
    test_initialize();
    let n = Atom::var(CS.nc);
    let adjoints = [
        parse_lit!(coad(Nc ^ 2 - 1, a), default_namespace = "spenso"),
        parse_lit!(coad(Nc ^ 2 - 1, b), default_namespace = "spenso"),
        parse_lit!(coad(Nc ^ 2 - 1, c), default_namespace = "spenso"),
    ];
    let factors = adjoints
        .iter()
        .map(|adjoint| color_t!(adjoint.clone()))
        .collect::<Vec<_>>();
    let start = parse_lit!(cof(Nc, i), default_namespace = "spenso");
    let end = parse_lit!(dind(cof(Nc, j)), default_namespace = "spenso");
    let generator = color_t!(adjoints[0].clone(), &start, &end);
    let open = chain!(&start, &end; factors[1..].iter().cloned());
    let spectator = parse_lit!((opaque(x) + opaque(y)) ^ 5 * (opaque(z) + opaque(w)));
    let settings = ColorSimplifySettings::default().with_cof_dimension_invariants();
    for reverse in [false, true] {
        let mut word = factors.clone();
        if reverse {
            word.reverse();
        }
        let closed = trace!(ColorFundamental {}.to_symbolic([n.clone()]); word);
        // The two-loop DIS networks contain Tr(Ta Tb Tc) Tb Tc and its
        // reversed quark trace. Independent Fierz contraction gives these
        // coefficients of the same free generator, with ratio -(Nc^2-2)/2.
        let coefficient = if reverse {
            (n.clone().pow(2) - Atom::num(2)) / (Atom::num(4) * n.clone())
        } else {
            -Atom::one() / (Atom::num(2) * n.clone())
        };
        let input = &closed * &open * &spectator;
        let simplified = simplify_formal_color(&input, settings);
        let scalar = simplified.coefficient((&generator * &spectator).as_view());
        assert_eq!(scalar.expand(), coefficient.expand(), "{simplified}");
        assert_eq!(simplified, scalar * &generator * &spectator);
        assert_eq!(simplify_formal_color(&simplified, settings), simplified);
        assert_eq!(
            crate::CookSettings::reversible()
                .with_representation_payloads(true, true)
                .cook(simplified.as_view())
                .list_dangling::<AbstractIndex>()
                .expect("formal reference slots"),
            crate::CookSettings::reversible()
                .with_representation_payloads(true, true)
                .cook(generator.as_view())
                .list_dangling::<AbstractIndex>()
                .expect("formal reference slots"),
        );
        // Trace evaluation and cross-line Fierz are independent policies.
        // Fierz may shorten a retained trace, but must not evaluate it.
        let unevaluated = simplify_formal_color(&input, settings.without_trace_evaluation());
        assert!(unevaluated.contains_symbol(spenso::network::tags::SPENSO_TAG.trace));
        assert_eq!(simplify_formal_color(&unevaluated, settings), simplified);
        assert_eq!(
            simplify_formal_color(
                &input,
                settings
                    .without_trace_evaluation()
                    .without_cross_chain_fierz_expansion()
            ),
            input,
        );
        // The preserved free index must not be captured by a trace dummy,
        // including when the input already reserves the first dummy labels.
        for index in [
            Atom::var(symbol!("idenso::x")),
            Atom::from(AbstractIndex::Dummy(0)),
            // This is a literal tensor label, although Symbolica patterns
            // interpret its trailing underscore as a wildcard.
            Atom::var(symbol!("free_adjoint_")),
        ] {
            let adjoint = ColorAdjoint {}.to_symbolic([n.clone().pow(2) - 1, index]);
            let internal = ColorAdjoint {}
                .to_symbolic([n.clone().pow(2) - 1, Atom::from(AbstractIndex::Dummy(1))]);
            let replacements = [
                symbolica::id::Replacement::new(
                    symbolica::id::Pattern::Literal(adjoints[0].clone()),
                    symbolica::id::Pattern::Literal(adjoint),
                ),
                symbolica::id::Replacement::new(
                    symbolica::id::Pattern::Literal(adjoints[1].clone()),
                    symbolica::id::Pattern::Literal(internal),
                ),
            ];
            let renamed = input.replace_multiple(&replacements);
            let expected = simplified.replace_multiple(&replacements);
            let actual = simplify_formal_color(&renamed, settings);
            assert_eq!(actual, expected);
            assert_eq!(
                crate::CookSettings::reversible()
                    .with_representation_payloads(true, true)
                    .cook(actual.as_view())
                    .list_dangling::<AbstractIndex>()
                    .expect("formal reference slots"),
                crate::CookSettings::reversible()
                    .with_representation_payloads(true, true)
                    .cook(expected.as_view())
                    .list_dangling::<AbstractIndex>()
                    .expect("formal reference slots"),
            );
        }
    }
}

#[test]
fn equivariant_projection_rejects_additional_adjoints_and_unknown_operators() {
    test_initialize();
    let _ = spenso::tensor_symbol!("spenso::unknown_color_operator");
    let expression = parse!(
        "chain(cof(Nc,i), dind(cof(Nc,j)),
            t(coad(Nc^2-1,a), spenso::in, spenso::out),
            t(coad(Nc^2-1,b), spenso::in, spenso::out))",
        default_namespace = "spenso"
    );
    assert_eq!(
        simplify_formal_color(&expression, ColorSimplifySettings::default()),
        expression
    );
    let unknown = parse_lit!(
        unknown_color_operator(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j))),
        default_namespace = "spenso"
    );
    assert_eq!(
        simplify_formal_color(&unknown, ColorSimplifySettings::default())
            .undo_chain::<AbstractIndex>()
            .unwrap(),
        unknown,
    );
}

#[test]
fn cubic_trace_contraction_preserves_two_open_color_lines() {
    test_initialize();
    let n = Atom::var(CS.nc);
    let adjoints = [
        parse_lit!(coad(Nc ^ 2 - 1, a), default_namespace = "spenso"),
        parse_lit!(coad(Nc ^ 2 - 1, b), default_namespace = "spenso"),
        parse_lit!(coad(Nc ^ 2 - 1, c), default_namespace = "spenso"),
    ];
    let factors = adjoints
        .iter()
        .map(|adjoint| color_t!(adjoint.clone()))
        .collect::<Vec<_>>();
    let start = parse_lit!(cof(Nc, i), default_namespace = "spenso");
    let end = parse_lit!(dind(cof(Nc, j)), default_namespace = "spenso");
    let other = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );
    let generator = color_t!(adjoints[0].clone(), &start, &end);
    let spectator = parse_lit!((opaque(x) + opaque(y)) ^ 5 * (opaque(z) + opaque(w)));
    // The adjoint joining the two open fundamental lines is summed. Compare
    // that internal label in one canonical chart; i,j,k,l remain free and ordered.
    let canonical = |expression: &Atom| {
        let cooking = crate::CookSettings::reversible().with_representation_payloads(true, true);
        let value = SymbolicTensor::infer(cooking.try_cook(expression.as_view()).unwrap())
            .unwrap()
            .canonize()
            .unwrap()
            .into_expression();
        cooking.uncook(value.as_view())
    };
    let tensor = canonical(&(&generator * &other * &spectator));
    let default_settings = ColorSimplifySettings::default().with_cof_dimension_invariants();
    let settings = default_settings.without_cross_chain_fierz_expansion();
    let fierz_tensor = simplify_formal_color(&(&generator * &other), default_settings);
    for reverse_trace in [false, true] {
        for reverse_chain in [false, true] {
            let mut word = factors.clone();
            let mut pair = factors[1..].to_vec();
            if reverse_trace {
                word.reverse();
            }
            if reverse_chain {
                pair.reverse();
            }
            let closed = trace!(ColorFundamental {}.to_symbolic([n.clone()]); word);
            let open = chain!(&start, &end; pair);
            let input = &closed * &open * &other * &spectator;
            // Independent Fierz contraction gives -1/(2 Nc) and
            // (Nc^2-2)/(4 Nc), hence the four-incoming DIS ratio -7/2
            // at Nc=3. The second open quark line is a spectator.
            let coefficient = if reverse_trace != reverse_chain {
                (n.clone().pow(2) - Atom::num(2)) / (Atom::num(4) * n.clone())
            } else {
                -Atom::one() / (Atom::num(2) * n.clone())
            };
            let simplified = canonical(&simplify_formal_color(&input, settings));
            let scalar = simplified.coefficient(tensor.as_view());
            assert_eq!(
                scalar.expand(),
                coefficient.expand(),
                "reverse_trace={reverse_trace}, reverse_chain={reverse_chain}; simplified={simplified}; tensor={tensor}"
            );
            let default_result = simplify_formal_color(&input, default_settings);
            let default_color = default_result.coefficient(spectator.as_view());
            assert_eq!(default_result, &default_color * &spectator);
            // Compare the independent coefficient in the existing Fierz basis.
            // Only the isolated color expression is expanded; the momentum
            // spectator remains factorized in both simplification modes.
            assert_eq!(
                default_color.expand(),
                (&coefficient * &fierz_tensor).expand()
            );
            assert_eq!(simplified, scalar * &tensor);
            assert_eq!(simplify_formal_color(&simplified, settings), simplified);
            assert_eq!(
                crate::CookSettings::reversible()
                    .with_representation_payloads(true, true)
                    .cook(simplified.as_view())
                    .list_dangling::<AbstractIndex>()
                    .expect("formal reference slots"),
                crate::CookSettings::reversible()
                    .with_representation_payloads(true, true)
                    .cook(tensor.as_view())
                    .list_dangling::<AbstractIndex>()
                    .expect("formal reference slots"),
            );
            assert_eq!(
                simplify_formal_color(&input, settings.without_trace_evaluation()),
                input,
            );
            for index in [
                Atom::var(symbol!("idenso::x")),
                Atom::from(AbstractIndex::Dummy(0)),
                Atom::var(symbol!("boundary_adjoint_")),
            ] {
                let adjoint = ColorAdjoint {}.to_symbolic([n.clone().pow(2) - 1, index]);
                let internal = ColorAdjoint {}
                    .to_symbolic([n.clone().pow(2) - 1, Atom::from(AbstractIndex::Dummy(1))]);
                let replacements = [
                    symbolica::id::Replacement::new(
                        symbolica::id::Pattern::Literal(adjoints[0].clone()),
                        symbolica::id::Pattern::Literal(adjoint),
                    ),
                    symbolica::id::Replacement::new(
                        symbolica::id::Pattern::Literal(adjoints[1].clone()),
                        symbolica::id::Pattern::Literal(internal),
                    ),
                ];
                assert_eq!(
                    canonical(&simplify_formal_color(
                        &input.replace_multiple(&replacements),
                        settings
                    )),
                    canonical(&simplified.replace_multiple(&replacements)),
                );
            }
        }
    }

    let before = parse_lit!(coad(Nc ^ 2 - 1, d), default_namespace = "spenso");
    let after = parse_lit!(coad(Nc ^ 2 - 1, e), default_namespace = "spenso");
    let open = chain!(&start, &end;
        [color_t!(before.clone()), factors[1].clone(), factors[2].clone(), color_t!(after.clone())]
    );
    let preserved_color = chain!(&start, &end;
        [color_t!(before), factors[0].clone(), color_t!(after)]
    ) * &other;
    let preserved = canonical(&(&preserved_color * &spectator));
    let symmetric =
        CS.symmetric_generator_trace(ColorFundamental {}.to_symbolic([n.clone()]), adjoints);
    let input = &symmetric * open * &other * &spectator;
    let simplified_color = simplify_formal_color(&input, settings);
    assert_eq!(
        simplify_formal_color(&simplified_color, settings),
        simplified_color
    );
    let simplified = canonical(&simplified_color);
    let scalar = simplified.coefficient(preserved.as_view());
    let coefficient = (n.clone().pow(2) - Atom::num(4)) / (Atom::num(8) * n);
    assert_eq!(scalar.expand(), coefficient.expand(), "{simplified}");
    let default_result = simplify_formal_color(&input, default_settings);
    let default_color = default_result.coefficient(spectator.as_view());
    assert_eq!(default_result, &default_color * &spectator);
    let fierz_tensor = simplify_formal_color(&preserved_color, default_settings);
    assert_eq!(
        default_color.expand(),
        (&coefficient * fierz_tensor).expand()
    );
    assert_eq!(simplified, scalar * preserved);
    assert_eq!(
        canonical(&simplify_formal_color(&simplified, settings)),
        simplified
    );
}

#[test]
fn adjoint_chain_endpoint_contraction_vanishes_before_materialization() {
    test_initialize();
    let settings = ColorSimplifySettings::default().without_trace_evaluation();
    let spectator = parse_lit!((x + y) ^ 7);
    for word in [
        "chain(coad(dA,a),coad(dA,b),f(in,out,coad(dA,a)),f(in,out,coad(dA,c)))",
        "chain(coad(dA,a),coad(dA,b),f(in,out,coad(dA,c)),f(in,out,coad(dA,b)))",
        "chain(coad(dA,a),coad(dA,b),f(out,in,coad(dA,c)),f(out,in,coad(dA,b)))",
    ] {
        let word = Atom::parse(word, "spenso", Default::default()).unwrap();
        assert_eq!(
            simplify_formal_color(&(&word * &spectator), settings),
            Atom::Zero
        );
    }
    let zero = parse!(
        "chain(coad(dA,b),coad(dA,d_0),f(in,out,coad(dA,a)),f(in,out,coad(dA,d_0)))",
        default_namespace = "spenso"
    );
    let expected = parse_lit!(
        g(coad(dA, a), coad(dA, b)) * cas(2, coad(dA)) ^ 3 / 2,
        default_namespace = "spenso"
    );
    let source = (&zero * parse_lit!(cas(2, coad(dA)) ^ 2 / 12, default_namespace = "spenso")
        + &expected)
        * &spectator;
    assert_eq!(
        simplify_formal_color(&source, settings),
        &expected * &spectator
    );
    let open = parse!(
        "chain(coad(dA,a),coad(dA,b),f(in,out,coad(dA,c)))",
        default_namespace = "spenso"
    );
    assert_eq!(simplify_formal_color(&open, settings), open);
}

#[test]
fn adjoint_chain_repeated_pair_agrees_with_explicit_structure_contraction() {
    test_initialize();
    let settings = ColorSimplifySettings::default().without_trace_evaluation();
    let spectator = parse_lit!((x + y) ^ 7);
    let expected = -parse_lit!(
        cas(2, coad(dA)) * g(coad(dA, i), coad(dA, j)),
        default_namespace = "spenso"
    ) * &spectator;
    for (first, second, explicit, sign) in [
        (
            "f(in,out,coad(dA,a))",
            "f(in,out,coad(dA,a))",
            "f(coad(dA,i),coad(dA,k),coad(dA,a))*f(coad(dA,k),coad(dA,j),coad(dA,a))",
            1,
        ),
        (
            "f(out,in,coad(dA,a))",
            "f(in,out,coad(dA,a))",
            "f(coad(dA,k),coad(dA,i),coad(dA,a))*f(coad(dA,k),coad(dA,j),coad(dA,a))",
            -1,
        ),
        (
            "f(out,in,coad(dA,a))",
            "f(out,in,coad(dA,a))",
            "f(coad(dA,k),coad(dA,i),coad(dA,a))*f(coad(dA,j),coad(dA,k),coad(dA,a))",
            1,
        ),
    ] {
        let input = Atom::parse(
            format!("chain(coad(dA,i),coad(dA,j),{first},{second})"),
            "spenso",
            Default::default(),
        )
        .unwrap()
            * &spectator;
        let explicit = Atom::parse(explicit, "spenso", Default::default()).unwrap() * &spectator;
        let result = simplify_formal_color(&input, settings);
        assert_eq!(result, &expected * sign);
        assert_eq!(result, simplify_formal_color(&explicit, settings));
        assert_eq!(simplify_formal_color(&result, settings), result);
    }
}

#[test]
fn cubic_trace_contraction_rejects_unowned_or_incompatible_slots() {
    test_initialize();
    let _ = spenso::tensor_symbol!("spenso::unknown_color_operator");
    let settings = ColorSimplifySettings::default().without_cross_chain_fierz_expansion();
    let symmetric = parse!(
        "trace(cof(Nc), sym(t(coad(Nc^2-1,a), spenso::in, spenso::out),
            t(coad(Nc^2-1,b), spenso::in, spenso::out),
            t(coad(Nc^2-1,c), spenso::in, spenso::out)))",
        default_namespace = "spenso"
    );
    let spectator = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );
    for open in [
        // The contracted pair must be adjacent in the ordered chain.
        parse!(
            "chain(cof(Nc,i), dind(cof(Nc,j)),
            t(coad(Nc^2-1,b), spenso::in, spenso::out),
            t(coad(Nc^2-1,d), spenso::in, spenso::out),
            t(coad(Nc^2-1,c), spenso::in, spenso::out))",
            default_namespace = "spenso"
        ),
        // An unknown matrix cannot certify a known invariant color network.
        parse!(
            "chain(cof(Nc,i), dind(cof(Nc,j)),
            unknown_color_operator(spenso::in, spenso::out),
            t(coad(Nc^2-1,b), spenso::in, spenso::out),
            t(coad(Nc^2-1,c), spenso::in, spenso::out))",
            default_namespace = "spenso"
        ),
        parse!(
            "chain(cof(Nc+1,i), dind(cof(Nc+1,j)),
            t(coad(Nc^2-1,b), spenso::in, spenso::out),
            t(coad(Nc^2-1,c), spenso::in, spenso::out))",
            default_namespace = "spenso"
        ),
    ] {
        let input = &symmetric * open * &spectator;
        assert_eq!(simplify_formal_color(&input, settings), input);
    }
    let open = parse!(
        "chain(cof(Nc,i), dind(cof(Nc,j)),
        t(coad(Nc^2-1,b), spenso::in, spenso::out),
        t(coad(Nc^2-1,c), spenso::in, spenso::out))",
        default_namespace = "spenso"
    );
    let third_occurrence = parse_lit!(
        t(coad(Nc ^ 2 - 1, b), cof(Nc, m), dind(cof(Nc, n))),
        default_namespace = "spenso"
    );
    let input = &symmetric * &open * &spectator * third_occurrence;
    let cooking = crate::CookSettings::reversible().with_representation_payloads(true, true);
    let admitted = cooking.try_cook(input.as_view()).unwrap();
    assert!(matches!(
        SymbolicTensor::infer(admitted),
        Err(crate::tensor::inference::TensorInferenceError::Composition(
            crate::tensor::composition::TensorCompositionError::InvalidExplicitMultiplicity {
                occurrences: 3,
                ..
            }
        ))
    ));
    for incompatible_spectator in [
        parse_lit!(
            t(coad(Nc ^ 2 - 1, a), cof(Nc + 1, k), dind(cof(Nc + 1, l))),
            default_namespace = "spenso"
        ),
        parse_lit!(
            t(coad(Nc ^ 2, a), cof(Nc, k), dind(cof(Nc, l))),
            default_namespace = "spenso"
        ),
    ] {
        let input = &symmetric * &open * incompatible_spectator;
        // Unsupported dimensions may retain a one-generator chain. Normalize
        // only that notation while requiring the contraction to stay untouched.
        assert_eq!(
            simplify_formal_color(&input, settings).undo_single_length(),
            input.undo_single_length(),
        );
    }
    let wrong_dimension = symmetric
        .replace(symbolica::id::Pattern::Literal(parse_lit!(
            coad(Nc ^ 2 - 1, b),
            default_namespace = "spenso"
        )))
        .with(symbolica::id::Pattern::Literal(parse_lit!(
            coad(dA, b),
            default_namespace = "spenso"
        )));
    let input = wrong_dimension * open * spectator;
    assert_eq!(simplify_formal_color(&input, settings), input);
}

#[test]
fn schur_reduction_keeps_free_adjoint_slots_and_unknown_tensors() {
    test_initialize();
    let _ = spenso::tensor_symbol!("spenso::unknown_color_operator");
    let generator = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j))),
        default_namespace = "spenso"
    );
    // A free adjoint label makes this a covariant operator, whose trace is
    // zero although the original operator is nonzero.
    assert_eq!(
        simplify_formal_color(&generator, ColorSimplifySettings::default()),
        generator
    );
    let unknown = parse_lit!(
        unknown_color_operator(cof(Nc, i), dind(cof(Nc, j))),
        default_namespace = "spenso"
    );
    // Simplification may represent the same operator as a chain shorthand.
    assert_eq!(
        simplify_formal_color(&unknown, ColorSimplifySettings::default())
            .undo_chain::<AbstractIndex>()
            .unwrap(),
        unknown,
    );
}

#[test]
fn fierz_generator_contraction() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );

    // Strict representation dimensions are atomic. Restore the FORM dimension
    // after the typed operation, before comparing with the original reference.
    let dimension = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let atomic_dimension = TestReps::new().coad_da.dim.to_symbolic();
    let admitted = expr
        .replace(dimension.to_pattern())
        .with(atomic_dimension.to_pattern());
    assert_eq!(
        admitted
            .replace(atomic_dimension.to_pattern())
            .with(dimension.to_pattern()),
        expr
    );
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression()
        .replace(atomic_dimension.to_pattern())
        .with(dimension.to_pattern());

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((result.expand()).as_atom_view().to_owned()).unwrap().contract(crate::tensor::ContractSettings::default().without_rank_one_tensors()).unwrap().into_expression().to_bare_ordered_string(), @"-1*Nc^(-1)*g(cof(Nc,i),dind(cof(Nc,j)))*g(cof(Nc,k),dind(cof(Nc,l)))*idx(2,cof(Nc))+g(cof(Nc,i),dind(cof(Nc,l)))*g(cof(Nc,k),dind(cof(Nc,j)))*idx(2,cof(Nc))");
}

#[test]
fn separated_generator_casimir_shortcut() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, k)))
            * t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );

    // Strict representation dimensions are atomic. Restore the FORM dimension
    // after the typed operation, before comparing with the original reference.
    let dimension = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");
    let atomic_dimension = TestReps::new().coad_da.dim.to_symbolic();
    let admitted = expr
        .replace(dimension.to_pattern())
        .with(atomic_dimension.to_pattern());
    assert_eq!(
        admitted
            .replace(atomic_dimension.to_pattern())
            .with(dimension.to_pattern()),
        expr
    );
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression()
        .replace(atomic_dimension.to_pattern())
        .with(dimension.to_pattern());

    use crate::shorthands::UndoShorthands;

    // Keep the original FORM identity, allowing the equivalent Fierz result
    // -T_R/N. Factor out exactly the same generator before reducing only its
    // scalar coefficient; the tensor numerator is never distributed.
    let expected = parse!(
        "-1/2*cas(2,coad(-1+Nc^2))*t(coad(-1+Nc^2,b),cof(Nc,i),dind(cof(Nc,l)))+cas(2,cof(Nc))*t(coad(-1+Nc^2,b),cof(Nc,i),dind(cof(Nc,l)))",
        default_namespace = "spenso"
    );
    let generator = parse_lit!(
        t(coad(Nc ^ 2 - 1, b), cof(Nc, i), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );
    let coefficient = |expression: Atom| {
        // The graph's slot grammar still uses atomic dimensions when undoing
        // the chain. Restore the original dimension before scalar comparison.
        let factored = expression
            .replace(dimension.to_pattern())
            .with(atomic_dimension.to_pattern())
            .undo_chain::<AbstractIndex>()
            .unwrap()
            .replace(atomic_dimension.to_pattern())
            .with(dimension.to_pattern())
            .collect_factors();
        let scalar = factored
            .replace(generator.to_pattern())
            .with(Atom::num(1).to_pattern());
        assert_eq!(factored, &scalar * &generator);
        scalar.to_cof_dimension_invariants().together().cancel()
    };
    assert_eq!(coefficient(result), coefficient(expected));
}

#[test]
fn chain_one_generator_trace_vanishes() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_da, a)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(crate::color::ColorSimplifySettings::default()), ..Default::default() }).unwrap().contract(crate::tensor::ContractSettings { collect_chains: false, collect_traces: false, ..Default::default() }).unwrap().into_expression().to_bare_ordered_string(), @"0");
}

#[test]
fn separated_fundamental_generators_contract_against_structure_constant() {
    test_initialize();
    let expression = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, k)))
            * t(coad(Nc ^ 2 - 1, c), cof(Nc, k), dind(cof(Nc, l)))
            * f(
                coad(Nc ^ 2 - 1, c),
                coad(Nc ^ 2 - 1, e),
                coad(Nc ^ 2 - 1, a)
            ),
        default_namespace = "spenso"
    );
    // Keep TR symbolic, as in the other trace identities. For Pauli/2
    // generators TR=1/2 and the coefficient is independently -i/4.
    let expected = -Atom::i()
        * parse_lit!(
            idx(2, cof(Nc))
                ^ 2 * g(coad(Nc ^ 2 - 1, b), coad(Nc ^ 2 - 1, e)) * g(cof(Nc, i), dind(cof(Nc, l))),
            default_namespace = "spenso"
        );
    assert_eq!(
        simplify_formal_color(&expression, ColorSimplifySettings::default()),
        expected
    );

    let reversed = parse_lit!(
        t(coad(Nc ^ 2 - 1, c), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, k)))
            * t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l)))
            * f(
                coad(Nc ^ 2 - 1, c),
                coad(Nc ^ 2 - 1, e),
                coad(Nc ^ 2 - 1, a)
            ),
        default_namespace = "spenso"
    );
    assert_eq!(
        simplify_formal_color(&reversed, ColorSimplifySettings::default()),
        -expected
    );

    // The three missed DIS color zeros close the same subchain through a
    // fourth generator and f(b,d,e); the metric sets b=e and kills f.
    let closing_factor = parse_lit!(
        t(coad(Nc ^ 2 - 1, d), cof(Nc, l), dind(cof(Nc, m)))
            * f(
                coad(Nc ^ 2 - 1, b),
                coad(Nc ^ 2 - 1, d),
                coad(Nc ^ 2 - 1, e)
            ),
        default_namespace = "spenso"
    );
    assert!(
        simplify_formal_color(
            &(expression * closing_factor),
            ColorSimplifySettings::default()
        )
        .is_zero()
    );
}

#[test]
fn chain_two_generator_trace_normalizes() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(crate::color::ColorSimplifySettings::default()), ..Default::default() }).unwrap().contract(crate::tensor::ContractSettings { collect_chains: false, collect_traces: false, ..Default::default() }).unwrap().into_expression().to_bare_ordered_string(), @"g(coad(dA,a),coad(dA,b))*idx(2,cof(Nc))");
}

#[test]
fn three_generator_trace_terminal() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, c)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(crate::color::ColorSimplifySettings::default()), ..Default::default() }).unwrap().contract(crate::tensor::ContractSettings { collect_chains: false, collect_traces: false, ..Default::default() }).unwrap().into_expression().to_bare_ordered_string(), @"(𝑖/2)*f(coad(dA,a),coad(dA,b),coad(dA,c))*idx(2,cof(Nc))+trace(cof(Nc),sym(t(coad(dA,a),in,out),t(coad(dA,b),in,out),t(coad(dA,c),in,out)))");
}

#[test]
fn adjacent_generator_casimir_chain() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), k),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, a)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(crate::color::ColorSimplifySettings::default()), ..Default::default() }).unwrap().contract(crate::tensor::ContractSettings { collect_chains: false, collect_traces: false, ..Default::default() }).unwrap().into_expression().to_bare_ordered_string(), @"cas(2,cof(Nc))*g(cof(Nc,i),dind(cof(Nc,k)))");
}

#[test]
fn separated_generator_casimir_trace() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, c)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(crate::color::ColorSimplifySettings::default()), ..Default::default() }).unwrap().contract(crate::tensor::ContractSettings { collect_chains: false, collect_traces: false, ..Default::default() }).unwrap().into_expression().to_bare_ordered_string(), @"-1/2*cas(2,coad(dA))*g(coad(dA,b),coad(dA,c))*idx(2,cof(Nc))+cas(2,cof(Nc))*g(coad(dA,b),coad(dA,c))*idx(2,cof(Nc))");
}

#[test]
fn two_f_loop_contracts_to_ca_metric() {
    test_initialize();
    let expr = parse_lit!(
        f(coad(dA, a), coad(dA, c), coad(dA, x)) * f(coad(dA, b), coad(dA, c), coad(dA, x)),
        default_namespace = "spenso"
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(crate::color::ColorSimplifySettings::default()), ..Default::default() }).unwrap().contract(crate::tensor::ContractSettings { collect_chains: false, collect_traces: false, ..Default::default() }).unwrap().into_expression().to_bare_ordered_string(), @"cas(2,coad(dA))*g(coad(dA,a),coad(dA,b))");
}

#[test]
fn metric_contraction_respects_antisymmetric_normalization() {
    test_initialize();
    let coad8 = ColorAdjoint {}.new_rep(8);
    // Seed the symbol order from the three-loop graph that exposed the lost
    // permutation parity in the old sequence-wildcard contraction.
    let a = slot!(coad8, standalone_metric_a);
    let x = slot!(coad8, standalone_metric_x);
    let u = slot!(coad8, standalone_metric_u);
    let v = slot!(coad8, standalone_metric_v);
    let _ = crate::tensor::SymbolicTensor::infer(
        (color_f!(a, x, u) * color_f!(a, x, v))
            .as_atom_view()
            .to_owned(),
    )
    .unwrap()
    .contract(crate::tensor::ContractSettings::default().without_rank_one_tensors())
    .unwrap()
    .into_expression();
    let r = slot!(coad8, standalone_metric_r);
    let s = slot!(coad8, standalone_metric_s);
    let j = slot!(coad8, standalone_metric_j);

    let contracted = crate::tensor::SymbolicTensor::infer(
        (g!(r, s) * color_f!(u, j, r)).as_atom_view().to_owned(),
    )
    .unwrap()
    .contract(crate::tensor::ContractSettings::default().without_rank_one_tensors())
    .unwrap()
    .into_expression();
    // Re-normalizing f after direct slot substitution may change its displayed
    // argument order and coefficient, so compare with a freshly built target.
    assert_eq!(
        contracted.to_bare_ordered_string(),
        color_f!(u, j, s).to_bare_ordered_string(),
    );

    let closed = g!(u, v) * g!(r, s) * color_f!(u, j, r) * color_f!(v, j, s);
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((closed).as_atom_view().to_owned())
            .unwrap()
            .simplify_algebra(&crate::tensor::AlgebraSettings {
                color: Some(crate::color::ColorSimplifySettings::default()),
                ..Default::default()
            })
            .unwrap()
            .contract(crate::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression(),
        Atom::num(8) * color_cas!(2, &coad8)
    );
}

#[test]
fn metric_contraction_only_replaces_the_immediate_slot() {
    test_initialize();
    let coad8 = ColorAdjoint {}.new_rep(8);
    let u = slot!(coad8, standalone_nested_metric_u);
    let j = slot!(coad8, standalone_nested_metric_j);
    let r = slot!(coad8, standalone_nested_metric_r);
    let s = slot!(coad8, standalone_nested_metric_s);
    let probe = symbol!("spenso::standalone_metric_probe");
    let structure = color_f!(u, j, r);
    let r_atom = r.into_atom();
    let s_atom = s.into_atom();

    // This raw owner regression intentionally leaves the outer wrapper untagged:
    // its nested tensor-looking metadata must not be treated as a second port.
    let source = g!(r, s) * function!(probe, &structure, r_atom);
    let contracted =
        crate::shorthands::schoonschip::SlotContraction::new().run(source.as_view(), true, false);
    assert_eq!(contracted, function!(probe, structure, s_atom));
}

#[test]
fn three_f_loop_contracts_to_ca_f() {
    test_initialize();
    let expr = parse_lit!(
        f(coad(dA, a), coad(dA, b), coad(dA, e))
            * f(coad(dA, b), coad(dA, c), coad(dA, f_))
            * f(coad(dA, c), coad(dA, a), coad(dA, g_)),
        default_namespace = "spenso"
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(crate::color::ColorSimplifySettings::default()), ..Default::default() }).unwrap().contract(crate::tensor::ContractSettings { collect_chains: false, collect_traces: false, ..Default::default() }).unwrap().into_expression().to_bare_ordered_string(), @"1/2*cas(2,coad(dA))*f(coad(dA,e),coad(dA,f_),coad(dA,g_))");
}

#[test]
fn three_f_loop_preserves_antisymmetric_orientation() {
    test_initialize();
    let coad8 = ColorAdjoint {}.new_rep(8);
    // Seed the symbol order that made the unoriented triangle reduction pick
    // the wrong sign in the three-loop graph.
    let r = slot!(coad8, standalone_triangle_r);
    let s = slot!(coad8, standalone_triangle_s);
    let x = slot!(coad8, standalone_triangle_x);
    let z = slot!(coad8, standalone_triangle_z);
    let _ = g!(r, s) * color_f!(x, z, r);
    let u = slot!(coad8, standalone_triangle_u);
    let y = slot!(coad8, standalone_triangle_y);
    let triangle = color_f!(u, x, r) * color_f!(y, z, r) * color_f!(u, z, s);
    let closing_structure = color_f!(s, x, y);
    let expected = Atom::num(-1) * color_cas!(2, &coad8) / Atom::num(2) * color_f!(s, x, y);

    assert_eq!(
        crate::tensor::SymbolicTensor::infer((triangle).as_atom_view().to_owned())
            .unwrap()
            .simplify_algebra(&crate::tensor::AlgebraSettings {
                color: Some(crate::color::ColorSimplifySettings::default()),
                ..Default::default()
            })
            .unwrap()
            .contract(crate::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression(),
        expected
    );
    // Closing the loop before or after triangle reduction must keep the same
    // antisymmetric orientation.
    assert_eq!(
        crate::tensor::SymbolicTensor::infer(
            (&triangle * &closing_structure).as_atom_view().to_owned()
        )
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression(),
        crate::tensor::SymbolicTensor::infer(
            (crate::tensor::SymbolicTensor::infer((triangle).as_atom_view().to_owned())
                .unwrap()
                .simplify_algebra(&crate::tensor::AlgebraSettings {
                    color: Some(crate::color::ColorSimplifySettings::default()),
                    ..Default::default()
                })
                .unwrap()
                .contract(crate::tensor::ContractSettings {
                    collect_chains: false,
                    collect_traces: false,
                    ..Default::default()
                })
                .unwrap()
                .into_expression()
                * closing_structure)
                .as_atom_view()
                .to_owned()
        )
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression(),
    );
}

#[test]
fn six_f_k33_odd_automorphism_simplifies_to_zero() {
    test_initialize();
    let r = TestReps::new();
    let contraction = fco!(r, a, b, c)
        * fco!(r, d, e, f)
        * fco!(r, h, i, j)
        * fco!(r, a, d, h)
        * fco!(r, b, e, i)
        * fco!(r, c, f, j);
    // Exchanging the first two K3,3 vertices is a dummy relabeling, but it
    // reverses three f tensors. Thus C=-C and the contraction is zero.
    let odd_relabeling = fco!(r, d, e, f)
        * fco!(r, a, b, c)
        * fco!(r, h, i, j)
        * fco!(r, d, a, h)
        * fco!(r, e, b, i)
        * fco!(r, f, c, j);

    assert_eq!(odd_relabeling, -&contraction);
    // Color simplification applies signed tensor canonicalization only to its
    // extracted color factor.
    assert!(
        crate::tensor::SymbolicTensor::infer((contraction).as_atom_view().to_owned())
            .unwrap()
            .simplify_algebra(&crate::tensor::AlgebraSettings {
                color: Some(crate::color::ColorSimplifySettings::default()),
                ..Default::default()
            })
            .unwrap()
            .contract(crate::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression()
            .is_zero()
    );
}

#[test]
fn six_f_k33_with_structured_indices_simplifies_to_zero() {
    test_initialize();
    let cooking = crate::CookSettings::reversible()
        .with_index_payload_filter(None)
        .with_output_tags(["idenso::canonical_color_index"]);
    let hedge = symbol!("spenso::hedge");
    let structured_slot = |index| coad!(8, function!(hedge, Atom::num(index as i64)));
    let a = structured_slot(1);
    let b = structured_slot(2);
    let c = structured_slot(3);
    let d = structured_slot(4);
    let e = structured_slot(5);
    let f = structured_slot(6);
    let h = structured_slot(7);
    let i = structured_slot(8);
    let j = structured_slot(9);
    let contraction = color_f!(&a, &b, &c)
        * color_f!(&d, &e, &f)
        * color_f!(&h, &i, &j)
        * color_f!(&a, &d, &h)
        * color_f!(&b, &e, &i)
        * color_f!(&c, &f, &j);

    assert!(
        crate::tensor::SymbolicTensor::infer(cooking.cook_indices(contraction.as_view()))
            .unwrap()
            .simplify_algebra(&crate::tensor::AlgebraSettings {
                color: Some(crate::color::ColorSimplifySettings::default()),
                ..Default::default()
            })
            .unwrap()
            .contract(crate::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression()
            .is_zero()
    );
    // Pruning a closed zero factor must restore every surviving structured
    // label, including free ports on the remaining color tensor.
    let surviving = color_f!(
        structured_slot(21),
        structured_slot(22),
        structured_slot(23)
    );
    let input = (contraction + Atom::one()) * &surviving;
    let simplified = crate::tensor::SymbolicTensor::infer(cooking.cook_indices(input.as_view()))
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression();
    let simplified = cooking.uncook(simplified.as_view());
    assert_eq!(simplified, surviving);
    assert!(simplified.contains_symbol(hedge));
}

#[test]
fn signed_pruning_reenters_color_simplification() {
    use spenso::structure::slot::IsAbstractSlot;

    test_initialize();
    let r = TestReps::new();
    let odd = fco!(r, a, b, c)
        * fco!(r, d, e, f)
        * fco!(r, h, i, j)
        * fco!(r, a, d, h)
        * fco!(r, b, e, i)
        * fco!(r, c, f, j);
    let vanishing_factor = fco!(r, u, v, w);
    let surviving_factor = fco!(r, x, y, z);
    // An algebraically vanishing coefficient does not excuse incompatible
    // source ports: validate before signed-zero pruning.
    let incompatible = (&odd * &vanishing_factor + &surviving_factor) * &surviving_factor;
    assert!(crate::tensor::SymbolicTensor::infer(incompatible).is_err());
    let vanishing_factor = vanishing_factor
        .replace(slot!(r.coad_da, u).to_atom().to_pattern())
        .with(slot!(r.coad_da, x).to_atom().to_pattern())
        .replace(slot!(r.coad_da, v).to_atom().to_pattern())
        .with(slot!(r.coad_da, y).to_atom().to_pattern())
        .replace(slot!(r.coad_da, w).to_atom().to_pattern())
        .with(slot!(r.coad_da, z).to_atom().to_pattern());
    assert_eq!(vanishing_factor, surviving_factor);
    let expected = crate::tensor::SymbolicTensor::infer(
        (surviving_factor.clone().pow(2)).as_atom_view().to_owned(),
    )
    .unwrap()
    .simplify_algebra(&crate::tensor::AlgebraSettings {
        color: Some(crate::color::ColorSimplifySettings::default()),
        ..Default::default()
    })
    .unwrap()
    .contract(crate::tensor::ContractSettings {
        collect_chains: false,
        collect_traces: false,
        ..Default::default()
    })
    .unwrap()
    .into_expression();

    assert_eq!(
        crate::tensor::SymbolicTensor::infer(
            ((odd * vanishing_factor + &surviving_factor) * surviving_factor)
                .as_atom_view()
                .to_owned()
        )
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression(),
        expected
    );
}

#[test]
fn color_simplification_does_not_canonicalize_lorentz_tensors() {
    // This test selects family prerequisites only, preserving its explicit output scope.

    test_initialize();
    let r = TestReps::new();
    let a = mink!(4, color_only_canonicalization_a);
    let b = mink!(4, color_only_canonicalization_b);
    let c = mink!(4, color_only_canonicalization_c);
    let odd_lorentz =
        antisym!(a.clone(), b.clone()) * antisym!(b.clone(), c.clone()) * antisym!(c, a);
    let color = fco!(r, x, y, z);
    let expected = crate::tensor::SymbolicTensor::infer((color).as_atom_view().to_owned())
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            contract: crate::tensor::AlgebraContraction::None,
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression()
        * &odd_lorentz;

    assert!(!expected.is_zero());
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((color * &odd_lorentz).as_atom_view().to_owned())
            .unwrap()
            .simplify_algebra(&crate::tensor::AlgebraSettings {
                contract: crate::tensor::AlgebraContraction::None,
                color: Some(crate::color::ColorSimplifySettings::default()),
                ..Default::default()
            })
            .unwrap()
            .contract(crate::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap()
            .into_expression(),
        expected
    );
    assert!(
        odd_lorentz
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .expect("test expression should canonicalize")
            .is_zero()
    );
}

#[test]
fn two_f_even_automorphism_survives_canonicalization() {
    test_initialize();
    let r = TestReps::new();
    let contraction = fco!(r, a, b, c) * fco!(r, a, b, c);

    // Exchanging two contracted edges reverses both f tensors, so the total
    // graph-automorphism sign is positive.
    assert!(
        !contraction
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .expect("test expression should canonicalize")
            .is_zero()
    );
}

#[test]
fn four_f_closed_ghost_topology_preserves_antisymmetric_orientation() {
    test_initialize();
    let expr = parse_lit!(
        -1 / 8
            * g(coad(8, hedge(0)), coad(8, hedge(1)))
            * f(coad(8, hedge(3)), coad(8, hedge(0)), coad(8, hedge(5)))
            * f(coad(8, hedge(3)), coad(8, hedge(7)), coad(8, hedge(9)))
            * f(coad(8, hedge(7)), coad(8, hedge(11)), coad(8, hedge(1)))
            * f(coad(8, hedge(9)), coad(8, hedge(5)), coad(8, hedge(11))),
        default_namespace = "spenso"
    );

    // The FORM fixture labels indices with hedge(n); typed admission uses the
    // existing reversible index encoding and restores it at the Atom boundary.
    let cooking = crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
    let admitted = cooking.try_cook_indices(expr.as_view()).unwrap();
    assert_eq!(cooking.uncook(admitted.as_view()), expr);
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression();
    let result = cooking.uncook(result.as_view());

    assert_snapshot!(result.to_bare_ordered_string(), @"(cas(2,coad(8)))^2*1/2");
}

#[test]
fn factored_four_f_closed_ghost_topology_preserves_antisymmetric_orientation() {
    test_initialize();
    let expr = parse_lit!(
        -1 / 8
            * g(coad(8, hedge(0)), coad(8, hedge(1)))
            * (5 / 16
                * x
                * f(coad(8, hedge(3)), coad(8, hedge(0)), coad(8, hedge(5)))
                * f(coad(8, hedge(3)), coad(8, hedge(7)), coad(8, hedge(9)))
                * f(coad(8, hedge(7)), coad(8, hedge(11)), coad(8, hedge(1)))
                * f(coad(8, hedge(9)), coad(8, hedge(5)), coad(8, hedge(11)))
                + 3 / 8
                    * y
                    * f(coad(8, hedge(3)), coad(8, hedge(0)), coad(8, hedge(5)))
                    * f(coad(8, hedge(3)), coad(8, hedge(7)), coad(8, hedge(9)))
                    * f(coad(8, hedge(7)), coad(8, hedge(11)), coad(8, hedge(1)))
                    * f(coad(8, hedge(9)), coad(8, hedge(5)), coad(8, hedge(11)))),
        default_namespace = "spenso"
    );

    // The FORM fixture labels indices with hedge(n); typed admission uses the
    // existing reversible index encoding and restores it at the Atom boundary.
    let cooking = crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
    let admitted = cooking.try_cook_indices(expr.as_view()).unwrap();
    assert_eq!(cooking.uncook(admitted.as_view()), expr);
    let result = crate::tensor::SymbolicTensor::infer(admitted)
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression();
    let result = cooking.uncook(result.as_view());

    assert_snapshot!(result.expand().to_bare_ordered_string(), @"(cas(2,coad(8)))^2*3/16*y+(cas(2,coad(8)))^2*5/32*x");
}

#[test]
fn mixed_trace_structure_contraction() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, d)),
    ) * color_f!(
        slot!(r.coad_da, a),
        slot!(r.coad_da, b),
        slot!(r.coad_da, c),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(crate::color::ColorSimplifySettings::default()), ..Default::default() }).unwrap().contract(crate::tensor::ContractSettings { collect_chains: false, collect_traces: false, ..Default::default() }).unwrap().into_expression().to_bare_ordered_string(), @"(𝑖/2)*cas(2,coad(dA))*g(coad(dA,c),coad(dA,d))*idx(2,cof(Nc))");
}

#[test]
fn symmetric_invariant_d33_partial_contraction() {
    test_initialize();
    let expr = parse_lit!(
        d(sym_x, coad(dA, a), coad(dA, b), coad(dA, c))
            * d(sym_y, coad(dA, a), coad(dA, b), coad(dA, d_)),
        default_namespace = "spenso"
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(crate::color::ColorSimplifySettings::default()), ..Default::default() }).unwrap().contract(crate::tensor::ContractSettings { collect_chains: false, collect_traces: false, ..Default::default() }).unwrap().into_expression().to_bare_ordered_string(), @"dA^(-1)*g(coad(dA,c),coad(dA,d_))*gram(3,sym_x,sym_y)");
}

#[test]
fn four_generator_trace_terminal() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_da, a)),
        color_t!(slot!(r.coad_da, b)),
        color_t!(slot!(r.coad_da, c)),
        color_t!(slot!(r.coad_da, d)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(crate::color::ColorSimplifySettings::default()), ..Default::default() }).unwrap().contract(crate::tensor::ContractSettings { collect_chains: false, collect_traces: false, ..Default::default() }).unwrap().into_expression().to_bare_ordered_string(), @"(𝑖/2)*f(coad(dA,a),coad(dA,b),coad(dA,d_0))*trace(cof(Nc),sym(t(coad(dA,c),in,out),t(coad(dA,d),in,out),t(coad(dA,d_0),in,out)))+(𝑖/2)*f(coad(dA,c),coad(dA,d),coad(dA,d_0))*trace(cof(Nc),sym(t(coad(dA,a),in,out),t(coad(dA,b),in,out),t(coad(dA,d_0),in,out)))+-1/3*chain(coad(dA,c),coad(dA,d),f(coad(dA,b),in,out),f(coad(dA,a),in,out))*idx(2,cof(Nc))+1/6*chain(coad(dA,c),coad(dA,d),f(coad(dA,a),in,out),f(coad(dA,b),in,out))*idx(2,cof(Nc))+trace(cof(Nc),sym(t(coad(dA,a),in,out),t(coad(dA,b),in,out),t(coad(dA,c),in,out),t(coad(dA,d),in,out)))");
}

// FORM color.h reference families. The expressions follow FORM's color.frm
// size-5 port of color.h's tloop.frm
// (https://github.com/form-dev/form/blob/fdf6af0ce520f7da2dfe0b0b61bf4cef396770be/check/extra/color.frm),
// and tloop.frm and su.frm of color.tar.gz
// (https://www.nikhef.nl/~form/maindir/packages/color/). Every reference is
// FORM 5.0.0 running color.h (docolor) or SUn.prc on the same input; the
// closed networks reduce at once with the closed-network reducer.

fn form_color_tloop_q10(r: &TestReps) -> Atom {
    trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    )
}

fn form_color_tloop_g10(r: &TestReps) -> Atom {
    fco!(r, i1, i2, j1)
        * fco!(r, i2, i3, j2)
        * fco!(r, i3, i4, j3)
        * fco!(r, i4, i5, j4)
        * fco!(r, i5, i6, j5)
        * fco!(r, i6, i7, j1)
        * fco!(r, i7, i8, j2)
        * fco!(r, i8, i9, j3)
        * fco!(r, i9, i10, j4)
        * fco!(r, i10, i1, j5)
}

fn form_color_tloop_qq5(r: &TestReps) -> Atom {
    trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    ) * trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    )
}

fn form_color_tloop_qg5(r: &TestReps) -> Atom {
    trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    ) * fco!(r, k1, k2, j1)
        * fco!(r, k2, k3, j2)
        * fco!(r, k3, k4, j3)
        * fco!(r, k4, k5, j4)
        * fco!(r, k5, k1, j5)
}

fn form_color_tloop_gg5(r: &TestReps) -> Atom {
    fco!(r, i1, i2, j1)
        * fco!(r, i2, i3, j2)
        * fco!(r, i3, i4, j3)
        * fco!(r, i4, i5, j4)
        * fco!(r, i5, i1, j5)
        * fco!(r, k1, k2, j1)
        * fco!(r, k2, k3, j2)
        * fco!(r, k3, k4, j3)
        * fco!(r, k4, k5, j4)
        * fco!(r, k5, k1, j5)
}

fn simplify_closed_reference(expr: Atom) -> Atom {
    crate::tensor::SymbolicTensor::infer(expr)
        .unwrap()
        .simplify_algebra(&crate::tensor::AlgebraSettings {
            color: Some(crate::color::ColorSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap()
        .contract(crate::tensor::ContractSettings {
            collect_chains: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap()
        .into_expression()
        .expand()
}

/// color.h's scalar notation in idenso's invariants: NA, I2R, cR, cA, and the
/// named invariants by color.h name and argument representations (`A` adjoint).
fn color_h_reference(r: &TestReps, form: &str) -> Atom {
    let fundamental = r.cof_nc.to_symbolic([]);
    let adjoint = r.coad_da.to_symbolic([]);
    let mut text = form.to_string();
    // Named invariants first: d44(R1,A1) -> named_invariant("d44", [R, A]).
    let mut out = String::new();
    while let Some(start) = text.find(|c: char| c == 'd') {
        let (head, rest) = text.split_at(start);
        out.push_str(head);
        let name_end = rest.find('(').unwrap();
        let close = rest.find(')').unwrap();
        let name = &rest[..name_end];
        let adjoint_args: Vec<bool> = rest[name_end + 1..close].split(',').map(|arg| arg.starts_with('A')).collect();
        let value = crate::color::network::named_invariant(name, &adjoint_args, fundamental.clone(), adjoint.clone())
            .unwrap_or_else(|| panic!("unknown invariant {name}"));
        out.push_str(&format!("({})", value.to_canonical_string()));
        text = rest[close + 1..].to_string();
    }
    out.push_str(&text);
    let replacements = [
        ("NA", r.coad_da.dim.to_symbolic().to_canonical_string()),
        ("I2R", color_idx!(2, &fundamental).to_canonical_string()),
        ("cR", color_cas!(2, &fundamental).to_canonical_string()),
        ("cA", color_cas!(2, &adjoint).to_canonical_string()),
    ];
    for (token, value) in replacements {
        out = out.replace(token, &format!("({value})"));
    }
    symbolica::try_parse!(&out.replace("i_", "𝑖")).unwrap()
}

fn assert_matches_color_h(r: &TestReps, actual: &Atom, form: &str) {
    let reference = color_h_reference(r, form);
    TestReps::assert_color_invariant_equivalence(
        actual,
        &reference,
        &r.cof_nc.to_symbolic([]),
        &r.coad_da.to_symbolic([]),
    );
}

#[test]
fn form_github_color_tloop_q10_size_5() {
    test_initialize();
    let r = TestReps::new();
    let actual = simplify_closed_reference(form_color_tloop_q10(&r));
    assert_matches_color_h(
        &r,
        &actual,
        "NA*I2R*cR^4-5*NA*I2R*cA*cR^3+35/4*NA*I2R*cA^2*cR^2-155/24*NA*I2R*cA^3*cR+125/72*NA*I2R*cA^4+5*d44(R1,A1)*cR-6*d44(R1,A1)*cA+1/3*d44(A1,A2)*I2R",
    );
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^2*(cas(2,cof(Nc)))^3*35/4*Nc+(cas(2,coad(dA)))^3*(cas(2,cof(Nc)))^2*-155/24*Nc+(cas(2,coad(dA)))^4*125/72*Nc*cas(2,cof(Nc))+(cas(2,cof(Nc)))^4*-5*Nc*cas(2,coad(dA))+(cas(2,cof(Nc)))^5*Nc+-6*cas(2,coad(dA))*gram(4,coad(dA),cof(Nc))+1/3*gram(4,coad(dA),coad(dA))*idx(2,cof(Nc))+5*cas(2,cof(Nc))*gram(4,coad(dA),cof(Nc))");
}

#[test]
fn form_github_color_tloop_g10_size_5() {
    test_initialize();
    let r = TestReps::new();
    let actual = simplify_closed_reference(form_color_tloop_g10(&r));
    assert_matches_color_h(&r, &actual, "-1/36*NA*cA^5+2/3*d44(A1,A2)*cA");
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^5*-1/36*dA+2/3*cas(2,coad(dA))*gram(4,coad(dA),coad(dA))");
}

#[test]
fn form_github_color_tloop_qq5_size_5() {
    test_initialize();
    let r = TestReps::new();
    let actual = simplify_closed_reference(form_color_tloop_qq5(&r));
    assert_matches_color_h(
        &r,
        &actual,
        "-5/144*NA*I2R^2*cA^3+11/48*d33(R1,R2)*cA^2-5/6*d44(R1,R2)*cA+1/6*d44(R1,A1)*I2R+d55(R1,R2)",
    );
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^2*11/48*gram(3,cof(Nc),cof(Nc))+(cas(2,coad(dA)))^3*-5/144*Nc*cas(2,cof(Nc))*idx(2,cof(Nc))+-5/6*cas(2,coad(dA))*gram(4,cof(Nc),cof(Nc))+1/6*gram(4,coad(dA),cof(Nc))*idx(2,cof(Nc))+gram(5,cof(Nc),cof(Nc))");
}

#[test]
fn form_github_color_tloop_qg5_size_5() {
    test_initialize();
    let r = TestReps::new();
    let actual = simplify_closed_reference(form_color_tloop_qg5(&r));
    assert_matches_color_h(
        &r,
        &actual,
        "-5/144*i_*NA*I2R*cA^4-3/4*d44(R1,A1)*i_*cA+1/12*d44(A1,A2)*i_*I2R",
    );
    assert_snapshot!(actual.to_bare_ordered_string(), @"(-3𝑖/4)*cas(2,coad(dA))*gram(4,coad(dA),cof(Nc))+(-5𝑖/144)*(cas(2,coad(dA)))^4*Nc*cas(2,cof(Nc))+(𝑖/12)*gram(4,coad(dA),coad(dA))*idx(2,cof(Nc))");
}

#[test]
fn form_github_color_tloop_gg5_size_5() {
    test_initialize();
    let r = TestReps::new();
    let actual = simplify_closed_reference(form_color_tloop_gg5(&r));
    assert_matches_color_h(&r, &actual, "5/144*NA*cA^5+2/3*d44(A1,A2)*cA");
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^5*5/144*dA+2/3*cas(2,coad(dA))*gram(4,coad(dA),coad(dA))");
}

/// color.h's simpli reduces f's in an environment of invariants; here the
/// environment is a foreign function of the two open adjoint labels.
#[test]
fn simpli_contracts_projected_f_pair() {
    test_initialize();
    let expr = parse_lit!(
        f(coad(dA, x), coad(dA, a), coad(dA, c))
            * f(coad(dA, x), coad(dA, b), coad(dA, c))
            * invariant_environment(a, b),
        default_namespace = "spenso"
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(crate::color::ColorSimplifySettings::default()), ..Default::default() }).unwrap().contract(crate::tensor::ContractSettings { collect_chains: false, collect_traces: false, ..Default::default() }).unwrap().into_expression().to_bare_ordered_string(), @"cas(2,coad(dA))*g(coad(dA,a),coad(dA,b))*invariant_environment(a,b)");
}

// color.tar.gz tloop.frm with SIZE = 3.

#[test]
fn tloop_qloop_size_3() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(&r.cof_nc, tco!(r, j1), tco!(r, j2), tco!(r, j3), tco!(r, j1), tco!(r, j2), tco!(r, j3));
    let actual = simplify_closed_reference(expr);
    assert_matches_color_h(&r, &actual, "NA*I2R*cR^2-3/2*NA*I2R*cA*cR+1/2*NA*I2R*cA^2");
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^2*1/2*Nc*cas(2,cof(Nc))+(cas(2,cof(Nc)))^2*-3/2*Nc*cas(2,coad(dA))+(cas(2,cof(Nc)))^3*Nc");
}

#[test]
fn tloop_gloop_size_3() {
    test_initialize();
    let r = TestReps::new();
    let expr = fco!(r, i1, i2, j1)
        * fco!(r, i2, i3, j2)
        * fco!(r, i3, i4, j3)
        * fco!(r, i4, i5, j1)
        * fco!(r, i5, i6, j2)
        * fco!(r, i6, i1, j3);
    assert_snapshot!(simplify_closed_reference(expr).to_bare_ordered_string(), @"0");
}

#[test]
fn tloop_qqloop_size_3() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(&r.cof_nc, tco!(r, j1), tco!(r, j2), tco!(r, j3))
        * trace!(&r.cof_nc, tco!(r, j1), tco!(r, j2), tco!(r, j3));
    let actual = simplify_closed_reference(expr);
    assert_matches_color_h(&r, &actual, "-1/4*NA*I2R^2*cA+d33(R1,R2)");
    assert_snapshot!(actual.to_bare_ordered_string(), @"-1/4*Nc*cas(2,coad(dA))*cas(2,cof(Nc))*idx(2,cof(Nc))+gram(3,cof(Nc),cof(Nc))");
}

#[test]
fn tloop_qgloop_size_3() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(&r.cof_nc, tco!(r, j1), tco!(r, j2), tco!(r, j3))
        * fco!(r, k1, k2, j1)
        * fco!(r, k2, k3, j2)
        * fco!(r, k3, k1, j3);
    let actual = simplify_closed_reference(expr);
    assert_matches_color_h(&r, &actual, "1/4*i_*NA*I2R*cA^2");
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^2*(𝑖/4)*Nc*cas(2,cof(Nc))");
}

#[test]
fn tloop_ggloop_size_3() {
    test_initialize();
    let r = TestReps::new();
    let expr = fco!(r, i1, i2, j1)
        * fco!(r, i2, i3, j2)
        * fco!(r, i3, i1, j3)
        * fco!(r, k1, k2, j1)
        * fco!(r, k2, k3, j2)
        * fco!(r, k3, k1, j3);
    let actual = simplify_closed_reference(expr);
    assert_matches_color_h(&r, &actual, "1/4*NA*cA^3");
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^3*1/4*dA");
}

/// The Coxeter (Heawood) graph: 14 structure constants, girth six.
#[test]
fn tloop_g14() {
    test_initialize();
    let r = TestReps::new();
    let expr = fco!(r, i1, i2, i3)
        * fco!(r, i1, i4, i5)
        * fco!(r, i2, i6, i7)
        * fco!(r, i3, i8, i9)
        * fco!(r, i4, i10, i11)
        * fco!(r, i5, i12, i13)
        * fco!(r, i6, i14, i15)
        * fco!(r, i7, i16, i17)
        * fco!(r, i8, i18, i19)
        * fco!(r, i9, i20, i21)
        * fco!(r, i10, i21, i15)
        * fco!(r, i13, i19, i14)
        * fco!(r, i17, i11, i18)
        * fco!(r, i12, i16, i20);
    let actual = simplify_closed_reference(expr);
    assert_matches_color_h(&r, &actual, "1/648*NA*cA^7-8/15*d444(A1,A2,A3)*cA+16/9*d644(A1,A2,A3)");
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^7*1/648*dA+-8/15*cas(2,coad(dA))*d444(coad(dA),coad(dA),coad(dA))+16/9*d644(coad(dA),coad(dA),coad(dA))");
}

/// Five quark loops and eight gluons (16 vertices).
#[test]
fn tloop_fiveq() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(&r.cof_nc, tco!(r, j1), tco!(r, j2), tco!(r, j3))
        * trace!(&r.cof_nc, tco!(r, j2), tco!(r, j4), tco!(r, j5))
        * trace!(&r.cof_nc, tco!(r, j4), tco!(r, j6), tco!(r, j7))
        * trace!(&r.cof_nc, tco!(r, j6), tco!(r, j1), tco!(r, j8))
        * trace!(&r.cof_nc, tco!(r, j3), tco!(r, j5), tco!(r, j7), tco!(r, j8));
    let actual = simplify_closed_reference(expr);
    assert_matches_color_h(
        &r,
        &actual,
        "1/192*NA*I2R^5*cA^3+1/4*d33(R1,R2)*I2R^3*cA^2+5/48*d33(R1,R2)*d33(R3,R4)*NA^-1*I2R*cA+5/48*d33(R1,R3)*d33(R2,R4)*NA^-1*I2R*cA+1/8*d33(R1,R4)*d33(R2,R3)*NA^-1*I2R*cA+1/16*d44(R1,A1)*I2R^4+3/8*d433(R3,R1,R2)*I2R^2*cA+1/2*d3333(R1,R2,R3,R4)*I2R*cA+d43333a(R5,R2,R1,R4,R3)",
    );
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^2*(idx(2,cof(Nc)))^3*1/4*gram(3,cof(Nc),cof(Nc))+(cas(2,coad(dA)))^3*(idx(2,cof(Nc)))^4*1/192*Nc*cas(2,cof(Nc))+(gram(3,cof(Nc),cof(Nc)))^2*1/3*cas(2,coad(dA))*dA^(-1)*idx(2,cof(Nc))+(idx(2,cof(Nc)))^2*3/8*cas(2,coad(dA))*d433(cof(Nc),cof(Nc),cof(Nc))+(idx(2,cof(Nc)))^4*1/16*gram(4,coad(dA),cof(Nc))+1/2*cas(2,coad(dA))*d3333(cof(Nc),cof(Nc),cof(Nc),cof(Nc))*idx(2,cof(Nc))+d43333a(cof(Nc),cof(Nc),cof(Nc),cof(Nc),cof(Nc))");
}

// color.tar.gz su.frm (SIMPLE, n = 3): SU(NF) with Tr(T^a T^b) = a delta^{ab}
// and nf per closed fundamental trace, from SUn.prc.

/// The generic result in SU(NF) for Tr(T^a T^b) = a: the T_F = 1/2 values,
/// times (2a)^(V/2) for V generators and structure constants, and nf per trace.
fn su_n_reference_value(actual: &Atom, r: &TestReps, vertices: i64, traces: i64) -> Atom {
    use crate::color::ColorSimplifier;
    let nf_dim = Atom::var(symbol!("color_su_reference::NF"));
    let a = Atom::var(symbol!("color_su_reference::a"));
    let nf = Atom::var(symbol!("color_su_reference::nf"));
    let fundamental = r.cof_nc.to_symbolic([]);
    let adjoint = r.coad_da.to_symbolic([]);
    let replacements = [
        symbolica::id::Replacement::new(fundamental.to_pattern(), ColorFundamental {}.to_symbolic([nf_dim.clone()]).to_pattern()),
        symbolica::id::Replacement::new(adjoint.to_pattern(), ColorAdjoint {}.to_symbolic([nf_dim.clone().pow(Atom::num(2)) - Atom::num(1)]).to_pattern()),
        symbolica::id::Replacement::new(r.cof_nc.dim.to_symbolic().to_pattern(), nf_dim.to_pattern()),
        symbolica::id::Replacement::new(r.coad_da.dim.to_symbolic().to_pattern(), (nf_dim.clone().pow(Atom::num(2)) - Atom::num(1)).to_pattern()),
    ];
    let value = actual.replace_multiple(&replacements).to_cof_dimension_invariants();
    let scale = (Atom::num(2) * a).pow(Atom::num(vertices) / Atom::num(2)) * nf.pow(Atom::num(traces));
    (value * scale).expand()
}

fn assert_su_reference(actual: &Atom, r: &TestReps, vertices: i64, traces: i64, form: &str) {
    let expected = symbolica::try_parse!(&form.replace("i_", "𝑖"), default_namespace = "color_su_reference").unwrap();
    let value = su_n_reference_value(actual, r, vertices, traces);
    assert!((value.clone() - &expected).expand().is_zero(), "{value} != {expected}");
}

fn adjoint_matrix(r: &TestReps, index: Symbol) -> Atom {
    CS.structure_f(
        r.coad_da.to_symbolic([Atom::var(index)]),
        Atom::var(spenso::network::tags::SPENSO_TAG.chain_in),
        Atom::var(spenso::network::tags::SPENSO_TAG.chain_out),
    )
}

fn sym_fundamental_trace(r: &TestReps, indices: &[Symbol]) -> Atom {
    CS.symmetric_generator_trace(r.cof_nc.to_symbolic([]), indices.iter().map(|i| r.coad_da.to_symbolic([Atom::var(*i)])))
}

fn sym_adjoint_trace(r: &TestReps, indices: &[Symbol]) -> Atom {
    spenso::shadowing::trace_sym(r.coad_da.to_symbolic([]), indices.iter().map(|i| adjoint_matrix(r, *i)))
}

fn adjoint_trace(r: &TestReps, indices: &[Symbol]) -> Atom {
    trace!(r.coad_da.to_symbolic([]); indices.iter().map(|i| adjoint_matrix(r, *i)))
}

fn fundamental_trace(r: &TestReps, indices: &[Symbol]) -> Atom {
    trace!(r.cof_nc.to_symbolic([]); indices.iter().map(|i| color_t!(r.coad_da.to_symbolic([Atom::var(*i)]))))
}

fn su_indices<const N: usize>() -> [Symbol; N] {
    std::array::from_fn(|k| symbol!(&format!("color_su_reference::j{}", k + 1)))
}

#[test]
fn su_f3f3() {
    test_initialize();
    let r = TestReps::new();
    let j = su_indices::<3>();
    let actual = simplify_closed_reference(sym_fundamental_trace(&r, &j) * fundamental_trace(&r, &j));
    // (Tr(j1,j2,j3)+Tr(j1,j3,j2))/2*Tr(j1,j2,j3)/nf^2
    assert_su_reference(&actual, &r, 6, 0, "2*a^3*NF^-1-5/2*a^3*NF+1/2*a^3*NF^3");
    assert_snapshot!(actual.to_bare_ordered_string(), @"gram(3,cof(Nc),cof(Nc))");
}

#[test]
fn su_f4f4() {
    test_initialize();
    let r = TestReps::new();
    let j = su_indices::<4>();
    let actual = simplify_closed_reference(sym_fundamental_trace(&r, &j) * fundamental_trace(&r, &j));
    assert_su_reference(&actual, &r, 8, 2, "-3*a^4*nf^2*NF^-2+4*a^4*nf^2-7/6*a^4*nf^2*NF^2+1/6*a^4*nf^2*NF^4");
    assert_snapshot!(actual.to_bare_ordered_string(), @"gram(4,cof(Nc),cof(Nc))");
}

#[test]
fn su_f4a4() {
    test_initialize();
    let r = TestReps::new();
    let j = su_indices::<4>();
    let actual = simplify_closed_reference(sym_fundamental_trace(&r, &j) * adjoint_trace(&r, &j));
    assert_su_reference(&actual, &r, 8, 1, "-2*a^4*nf*NF+5/3*a^4*nf*NF^3+1/3*a^4*nf*NF^5");
    assert_snapshot!(actual.to_bare_ordered_string(), @"gram(4,coad(dA),cof(Nc))");
}

#[test]
fn su_a4a4() {
    test_initialize();
    let r = TestReps::new();
    let j = su_indices::<4>();
    let actual = simplify_closed_reference(sym_adjoint_trace(&r, &j) * adjoint_trace(&r, &j));
    assert_su_reference(&actual, &r, 8, 0, "-24*a^4*NF^2+70/3*a^4*NF^4+2/3*a^4*NF^6");
    assert_snapshot!(actual.to_bare_ordered_string(), @"gram(4,coad(dA),coad(dA))");
}

#[test]
fn su_fnfn_n3() {
    test_initialize();
    let r = TestReps::new();
    let j = su_indices::<3>();
    let actual = simplify_closed_reference(fundamental_trace(&r, &j) * fundamental_trace(&r, &j));
    assert_su_reference(&actual, &r, 6, 2, "2*a^3*nf^2*NF^-1-2*a^3*nf^2*NF");
    assert_snapshot!(actual.to_bare_ordered_string(), @"-1/4*Nc*cas(2,coad(dA))*cas(2,cof(Nc))*idx(2,cof(Nc))+gram(3,cof(Nc),cof(Nc))");
}

#[test]
fn su_fnan_n3() {
    test_initialize();
    let r = TestReps::new();
    let j = su_indices::<3>();
    let actual = simplify_closed_reference(fundamental_trace(&r, &j) * adjoint_trace(&r, &j));
    assert_su_reference(&actual, &r, 6, 1, "-i_*a^3*nf*NF^2+i_*a^3*nf*NF^4");
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^2*(𝑖/4)*Nc*cas(2,cof(Nc))");
}

#[test]
fn su_anan_n3() {
    test_initialize();
    let r = TestReps::new();
    let j = su_indices::<3>();
    let actual = simplify_closed_reference(adjoint_trace(&r, &j) * adjoint_trace(&r, &j));
    assert_su_reference(&actual, &r, 6, 0, "-2*a^3*NF^3+2*a^3*NF^5");
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^3*1/4*dA");
}

#[test]
fn su_qloop_n3() {
    test_initialize();
    let r = TestReps::new();
    let [j1, j2, j3] = su_indices::<3>();
    let actual = simplify_closed_reference(fundamental_trace(&r, &[j1, j2, j3, j1, j2, j3]));
    assert_su_reference(&actual, &r, 6, 1, "-a^3*nf*NF^-2+a^3*nf*NF^2");
    assert_snapshot!(actual.to_bare_ordered_string(), @"(cas(2,coad(dA)))^2*1/2*Nc*cas(2,cof(Nc))+(cas(2,cof(Nc)))^2*-3/2*Nc*cas(2,coad(dA))+(cas(2,cof(Nc)))^3*Nc");
}

#[test]
fn su_gloop_n3() {
    test_initialize();
    let r = TestReps::new();
    let [j1, j2, j3] = su_indices::<3>();
    let actual = simplify_closed_reference(adjoint_trace(&r, &[j1, j2, j3, j1, j2, j3]));
    assert_su_reference(&actual, &r, 6, 0, "0");
    assert_snapshot!(actual.to_bare_ordered_string(), @"0");
}

/// The value of a generic result in SU(n) with T_F = 1/2.
fn su_n_number(actual: &Atom, r: &TestReps, n: i64) -> Atom {
    use crate::color::ColorSimplifier;
    let replacements = [
        symbolica::id::Replacement::new(r.cof_nc.to_symbolic([]).to_pattern(), ColorFundamental {}.to_symbolic([Atom::num(n)]).to_pattern()),
        symbolica::id::Replacement::new(r.coad_da.to_symbolic([]).to_pattern(), ColorAdjoint {}.to_symbolic([Atom::num(n * n - 1)]).to_pattern()),
        symbolica::id::Replacement::new(r.cof_nc.dim.to_symbolic().to_pattern(), Atom::num(n).to_pattern()),
        symbolica::id::Replacement::new(r.coad_da.dim.to_symbolic().to_pattern(), Atom::num(n * n - 1).to_pattern()),
    ];
    actual.replace_multiple(&replacements).to_cof_dimension_invariants().expand()
}

/// d6433a(R,R,R,R) = d_R^{a1..a6} d_R^{a1a2a7a8} d_R^{a3a4a5} d_R^{a6a7a8}
/// (16 vertices) reaches color.h's d6433a relation, whose d433 term has
/// -1/40 (color.h's makeds writes -1/4). Explicit SU(3) and SU(4)
/// generators give 103/20736 and 915/8192.
#[test]
fn d6433a_network_has_its_su_n_value() {
    test_initialize();
    let r = TestReps::new();
    let [j1, j2, j3, j4, j5, j6, j7, j8] = su_indices::<8>();
    let expr = sym_fundamental_trace(&r, &[j1, j2, j3, j4, j5, j6])
        * sym_fundamental_trace(&r, &[j1, j2, j7, j8])
        * sym_fundamental_trace(&r, &[j3, j4, j5])
        * sym_fundamental_trace(&r, &[j6, j7, j8]);
    let actual = simplify_closed_reference(expr);
    assert_eq!(su_n_number(&actual, &r, 3), Atom::num(103) / Atom::num(20736), "{actual}");
    assert_eq!(su_n_number(&actual, &r, 4), Atom::num(915) / Atom::num(8192), "{actual}");
}

/// d_R^{a a b1..b5} d_R^{b1..b5}: color.h's contract of a rank-7 d_R with a
/// summed pair gives (cR - 5/12 cA) d_R^(5) plus pure-number multiples of
/// d_R^(3) d_A^(4) (color.h writes I2R/120). Explicit SU(3) and SU(4)
/// generators give 35/324 and 3115/4096.
#[test]
fn rank_seven_contracted_pair_has_its_su_n_value() {
    test_initialize();
    let r = TestReps::new();
    let [a, b1, b2, b3, b4, b5] = su_indices::<6>();
    let expr = sym_fundamental_trace(&r, &[a, a, b1, b2, b3, b4, b5]) * sym_fundamental_trace(&r, &[b1, b2, b3, b4, b5]);
    let actual = simplify_closed_reference(expr);
    assert_eq!(su_n_number(&actual, &r, 3), Atom::num(35) / Atom::num(324), "{actual}");
    assert_eq!(su_n_number(&actual, &r, 4), Atom::num(3115) / Atom::num(4096), "{actual}");
}

/// A colour factor the closed-network reducer does not read (here explicit
/// symmetric tensors beside a closed f loop) keeps the term with the local
/// rules instead of being passed through as a spectator.
#[test]
fn colour_spectators_of_a_closed_network_are_reduced() {
    test_initialize();
    let expr = parse_lit!(
        f(coad(dA, a), coad(dA, b), coad(dA, c))
            ^ 2 * d(cof(Nc), coad(dA, e), coad(dA, g), coad(dA, h))
                * d(cof(Nc), coad(dA, e), coad(dA, g), coad(dA, k))
                * g(coad(dA, h), coad(dA, k)),
        default_namespace = "spenso"
    );
    let simplify = |settings: crate::color::ColorSimplifySettings| {
        crate::tensor::SymbolicTensor::infer(expr.clone())
            .unwrap()
            .simplify_algebra(&crate::tensor::AlgebraSettings { color: Some(settings), ..Default::default() })
            .unwrap()
            .into_expression()
            .expand()
    };
    let with = simplify(crate::color::ColorSimplifySettings::default());
    let without = simplify(crate::color::ColorSimplifySettings::default().without_closed_networks());
    assert_eq!(with, without);
    let printed = with.to_bare_ordered_string();
    assert!(printed.contains("gram(3,cof(Nc),cof(Nc))") && !printed.contains("d(cof"), "{printed}");
}

/// Five-loop background-field QCD self-energies leave this closed network
/// (gammaloop issue #255): d_F^{d0 d1 d2} d_A^{d2 d3 d4 d5} (F^d0 F^d5 F^d1)_{d3 d4}.
/// Transposing the odd chain between the symmetric legs d3, d4 reverses it
/// with a sign; relabelling the symmetric legs d0, d1 restores it, so it is 0.
/// Explicit SU(3) generators give 0 as well.
#[test]
fn odd_adjoint_chain_between_symmetric_legs_vanishes() {
    test_initialize();
    let expr = parse!(
        "chain(coad(dA,d3),coad(dA,d4),f(coad(dA,d0),in,out),f(coad(dA,d5),in,out),f(coad(dA,d1),in,out))\
         *trace(cof(Nc),sym(t(coad(dA,d0),in,out),t(coad(dA,d1),in,out),t(coad(dA,d2),in,out)))\
         *trace(coad(dA),sym(f(coad(dA,d2),in,out),f(coad(dA,d3),in,out),f(coad(dA,d4),in,out),f(coad(dA,d5),in,out)))",
        default_namespace = "spenso"
    );
    let settings = crate::tensor::AlgebraSettings {
        color: Some(crate::color::ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: crate::tensor::AlgebraContraction::Dots,
        ..Default::default()
    };
    let result = crate::tensor::SymbolicTensor::infer_indexed(expr)
        .unwrap()
        .simplify_algebra(&settings)
        .unwrap();
    assert!(result.expression().is_zero(), "{}", result.expression());
}
