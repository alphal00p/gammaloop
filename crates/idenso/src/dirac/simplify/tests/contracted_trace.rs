use super::*;

#[test]
fn contracted_sixteen_trace_matches_clifford_components_and_dimension_specialization() {
    use super::super::trace_kernel::tests::clifford_trace;
    use spenso::network::library::symbolic::ETS;

    let r = test_initialize();
    let settings = GammaSimplifySettings::default();
    let mut outputs = Vec::new();
    for mink in [&r.mink4, &r.mink_d] {
        let p = momenta(&mink.to_symbolic([]));
        let slots: [_; 4] = std::array::from_fn(|i| mink.pattern(Atom::num(100 + i)));
        let [a, b, c, d] = &slots;
        let input = trace!(r.bis4.to_symbolic([]); [
            a, &p[0], b, &p[1], c, &p[2], d, &p[3],
            a, &p[4], d, &p[5], c, &p[6], b, &p[7],
        ].map(|argument| gamma!(argument)));
        let result = crate::tensor::SymbolicTensor::infer((input).as_atom_view().to_owned())
            .unwrap()
            .simplify_algebra(&crate::tensor::AlgebraSettings {
                gamma: Some(settings),
                epsilon: settings.output == crate::dirac::GammaOutput::Reduced,
                ..Default::default()
            })
            .unwrap()
            .into_expression();
        assert_ne!(result, input);
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((result).as_atom_view().to_owned())
                .unwrap()
                .simplify_algebra(&crate::tensor::AlgebraSettings {
                    gamma: Some(settings),
                    epsilon: settings.output == crate::dirac::GammaOutput::Reduced,
                    ..Default::default()
                })
                .unwrap()
                .into_expression(),
            result
        );
        result.visitor(&mut |node| {
            assert!(!slots.iter().any(|slot| slot.as_view() == node));
            true
        });
        outputs.push(result);
    }
    let dimension = minkowski_dimension(r.mink_d.to_symbolic([]).as_view())
        .unwrap()
        .to_owned();
    let specialized = outputs[1]
        .replace(dimension.to_pattern())
        .with(Atom::num(4).to_pattern())
        .expand();
    assert_eq!(specialized, outputs[0].expand());

    // Reuse the independent blade-product oracle, specializing every metric
    // to a Euclidean orthonormal 4D basis. No trace recurrence is used here.
    let p = momenta(&r.mink4.to_symbolic([]));
    let mut seed = 71u64;
    let mut nonzero = false;
    for _ in 0..3 {
        let vectors: [[i64; 4]; 8] = std::array::from_fn(|_| {
            std::array::from_fn(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                ((seed >> 32) % 3) as i64 - 1
            })
        });
        let mut expected = 0;
        for assignment in 0usize..256 {
            let basis: [[i64; 4]; 4] = std::array::from_fn(|slot| {
                std::array::from_fn(|axis| i64::from(axis == ((assignment >> (2 * slot)) & 3)))
            });
            expected += clifford_trace(
                &[
                    basis[0], vectors[0], basis[1], vectors[1], basis[2], vectors[2], basis[3],
                    vectors[3], basis[0], vectors[4], basis[3], vectors[5], basis[2], vectors[6],
                    basis[1], vectors[7],
                ],
                false,
            );
        }
        nonzero |= expected != 0;
        let actual = outputs[0].replace_map(|node, _, out| {
            if let AtomView::Fun(metric) = node
                && metric.get_symbol() == ETS.metric
            {
                let indices: Vec<_> = metric
                    .iter()
                    .map(|argument| p.iter().position(|p| p.as_view() == argument).unwrap())
                    .collect();
                assert_eq!(indices.len(), 2);
                **out = Atom::num(
                    vectors[indices[0]]
                        .iter()
                        .zip(vectors[indices[1]])
                        .map(|(a, b)| a * b)
                        .sum::<i64>(),
                );
            }
        });
        assert_eq!(actual, Atom::num(expected));
    }
    assert!(
        nonzero,
        "the component oracle must not only test zero traces"
    );
}

#[test]
fn contracted_sixteen_mixed_components_keep_normalizer_semantics() {
    fn component(value: AtomView<'_>, out: &mut Settable<Atom>) {
        if let AtomView::Fun(vector) = value
            && let Some(AtomView::Fun(slot)) = vector.iter().last()
            && slot.get_symbol() == *MINKOWSKI_SYMBOL
            && slot.get_nargs() == 2
        {
            **out = Atom::num(1);
        }
    }
    let r = test_initialize();
    let rep = r.mink_d.to_symbolic([]);
    let p = symbolica::function!(
        spenso::vector_symbol!("idenso::contracted_sixteen_callback_p", norm = component),
        &rep
    );
    let q = symbolica::function!(
        spenso::vector_symbol!("idenso::contracted_sixteen_callback_q", norm = component),
        &rep
    );
    let a = slot!(r.mink_d, contracted_sixteen_a).into_atom();
    let b = slot!(r.mink_d, contracted_sixteen_b).into_atom();
    let c = slot!(r.mink_d, contracted_sixteen_c).into_atom();
    let input = trace!(r.bis4.to_symbolic([]);
        [&a, &a, &b, &p, &c, &q].into_iter()
            .chain(std::iter::repeat_n(&p, 10))
            .map(|argument| gamma!(argument)));
    let dimension = mink_slot_dimension(a.as_view()).unwrap();
    // The surviving components p(b), p(c), q(b), q(c) all normalize to one.
    // Treating their pre-normalization metrics as polynomial variables fails.
    let expected = (dimension
        * g!(&p, &p).pow(5).as_view()
        * (Atom::num(8) - Atom::num(4) * g!(&b, &c) * g!(&p, &q)))
    .expand();
    let settings = GammaSimplifySettings::default();
    let source = SymbolicTensor::<PartialStructure>::infer(input.clone()).unwrap();
    let result = DiracSimplifier::new(&settings)
        .evaluate_terminal_trace::<false>(input.as_view())
        .map(|(expression, _)| expression)
        .unwrap();
    assert_eq!(result.expand(), expected.expand());
    // The callback deliberately mixes a scalar with surviving b,c ports.
    // Keep its raw kernel oracle while rejecting the changed typed interface.
    assert!(
        source
            .simplify_algebra(&crate::tensor::AlgebraSettings {
                gamma: Some(settings),
                epsilon: settings.output == crate::dirac::GammaOutput::Reduced,
                ..Default::default()
            })
            .is_err()
    );
    assert!(source.with_rewritten_expression(result).is_err());
}

#[test]
fn contracted_sixteen_unsupported_coefficients_keep_atom_fallback() {
    test_initialize();
    let settings = GammaSimplifySettings::default();
    for (dimension, unit) in [
        (
            Atom::num(4),
            Atom::var(symbolica::symbol!("contracted_sixteen_spin")),
        ),
        (Atom::num(4), Atom::num(1) / Atom::num(2)),
        (Atom::num(6), Atom::num(4)),
    ] {
        let rep = symbolica::function!(*MINKOWSKI_SYMBOL, &dimension);
        let a = symbolica::function!(
            *MINKOWSKI_SYMBOL,
            &dimension,
            Atom::var(symbolica::symbol!("contracted_sixteen_fallback_a"))
        );
        let p = &momenta(&rep)[0];
        let spin = symbolica::function!(*BISPINOR_SYMBOL, &unit);
        let input = trace!(&spin;
            [&a, &a].into_iter().chain(std::iter::repeat_n(p, 14))
                .map(|argument| gamma!(argument)));
        let expected = &unit * &dimension * g!(p, p).pow(7);
        // These raw trace units intentionally exceed typed gamma admission.
        // Keep the explicit sparse oracle and its Atom-emission fallback covered
        // through the same local trace evaluator, without a public raw route.
        let result = DiracSimplifier {
            settings: &settings,
            output: trace_kernel::TraceOutput::Expanded,
        }
        .evaluate_terminal_trace::<false>(input.as_view())
        .map(|(expression, _)| expression)
        .expect("the local trace owner accepts these coefficient domains");
        assert_eq!(result.expand(), expected.expand());
        assert_eq!(
            settings.rewrite_expression(result.clone(), trace_kernel::TraceOutput::Expanded),
            result
        );
    }
}

#[test]
fn sixteen_trace_preserves_open_words_and_metadata() {
    let r = test_initialize();
    let rep = r.mink_d.to_symbolic([]);
    let metadata = symbolica::parse_lit!((contracted_sixteen_x + contracted_sixteen_y) ^ 3);
    let vectors: Vec<_> = (0..16)
        .map(|i| {
            FunctionBuilder::new(T.rank_one_tensor_symbol("idenso::contracted_sixteen_metadata"))
                .add_arg(&metadata)
                .add_arg(Atom::num(i))
                .add_arg(&rep)
                .finish()
        })
        .collect();
    let start = slot!(r.bis4, contracted_sixteen_start).into_atom();
    let end = slot!(r.bis4, contracted_sixteen_end).into_atom();
    let open = chain!(&start, &end; vectors.iter().map(|p| gamma!(p)));
    let closed = trace!(r.bis4.to_symbolic([]); vectors.iter().map(|p| gamma!(p)));
    let settings = GammaSimplifySettings::default();
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((open).as_atom_view().to_owned())
            .unwrap()
            .simplify_algebra(&crate::tensor::AlgebraSettings {
                gamma: Some(settings),
                epsilon: settings.output == crate::dirac::GammaOutput::Reduced,
                ..Default::default()
            })
            .unwrap()
            .into_expression(),
        open
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer((closed).as_atom_view().to_owned())
            .unwrap()
            .simplify_algebra(&crate::tensor::AlgebraSettings {
                gamma: Some(settings.without_trace_evaluation()),
                epsilon: settings.output == crate::dirac::GammaOutput::Reduced,
                ..Default::default()
            })
            .unwrap()
            .into_expression(),
        closed
    );
}
// Reuses the existing independent Clifford-blade oracle; no production numerator
// expansion and no duplicated gamma trace recurrence in the expected result.
#[test]
fn massive_routed_ladder_preserves_open_compound_ports_and_component_values() {
    use super::super::trace_kernel::tests::clifford_trace;
    use crate::tensor::{AlgebraSettings, ReductionStatus};
    use spenso::network::library::symbolic::ETS;

    let r = test_initialize();
    let label = spenso::index_symbol!("idenso::massive_ladder_test::edge");
    let mass = Atom::var(symbolica::symbol!("idenso::massive_ladder_test::mass"));
    let spin: [Atom; 8] = std::array::from_fn(|i| {
        r.bis4
            .to_symbolic([symbolica::function!(label, i as i64, 2)])
    });
    let coordinates = [[2, 1, -1, 0], [1, -2, 0, 1], [-1, 1, 2, 0], [2, 0, 1, -1]];
    let routed: [[i64; 4]; 4] = std::array::from_fn(|i| {
        std::array::from_fn(|j| match i {
            0 => coordinates[0][j],
            1 => coordinates[1][j] + coordinates[2][j],
            2 => coordinates[2][j] - coordinates[3][j],
            _ => coordinates[3][j],
        })
    });

    for mink in [&r.mink4, &r.mink_d] {
        let momenta = momenta(&mink.to_symbolic([]));
        let [mu, nu, alpha] =
            [40, 10, 30].map(|i| mink.to_symbolic([symbolica::function!(label, i, 2)]));
        let external = [&mu, &nu];
        let vertices = [&mu, &alpha, &nu, &alpha];
        let indexed = |i: usize, port: &Atom| {
            let AtomView::Fun(momentum) = momenta[i].as_view() else {
                unreachable!()
            };
            FunctionBuilder::new(momentum.get_symbol())
                .add_arg(port)
                .finish()
        };
        let mut factors = Vec::new();
        for i in 0..4 {
            let rho = mink.to_symbolic([symbolica::function!(label, 100 + i as i64, 2)]);
            let momentum = match i {
                0 => indexed(0, &rho),
                1 => indexed(1, &rho) + indexed(2, &rho),
                2 => indexed(2, &rho) - indexed(3, &rho),
                _ => indexed(3, &rho),
            };
            factors.push(
                &mass * g!(&spin[2 * i], &spin[2 * i + 1])
                    + momentum * gamma!(&spin[2 * i], &spin[2 * i + 1], &rho),
            );
            factors.push(gamma!(
                &spin[2 * i + 1],
                &spin[(2 * i + 2) % 8],
                vertices[i]
            ));
        }
        let source = SymbolicTensor::infer(Atom::mul_many(factors)).unwrap();
        let settings = AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            ..Default::default()
        };
        let reduced = source.simplify_algebra(&settings).unwrap();
        assert_eq!(reduced.structure, source.structure);
        assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
        let again = reduced.simplify_algebra(&settings).unwrap();
        assert_eq!(again.expression, reduced.expression);
        assert_eq!(again.structure, source.structure);

        let dimension = minkowski_dimension(r.mink_d.to_symbolic([]).as_view())
            .unwrap()
            .to_owned();
        let specialize = |value: &Atom| {
            value
                .replace(dimension.to_pattern())
                .with(Atom::num(4).to_pattern())
        };
        let component_expression = specialize(&reduced.expression);
        let component_momenta = momenta.iter().map(specialize).collect::<Vec<_>>();
        let component_external = external.map(specialize);
        for m in [0_i64, 2] {
            for first in 0..4 {
                for second in 0..4 {
                    let basis = |axis| std::array::from_fn(|i| i64::from(i == axis));
                    let endpoints = [basis(first), basis(second)];
                    let mut expected = 0_i64;
                    // Expand only the independent finite oracle's 4 scalar/vector
                    // matrix choices. The production input remains factorized.
                    for mask in 0_u32..16 {
                        for axis in 0..4 {
                            let mut word = Vec::new();
                            for (position, momentum) in routed.iter().enumerate() {
                                if mask & (1 << position) != 0 {
                                    word.push(*momentum);
                                }
                                word.push(match position {
                                    0 => endpoints[0],
                                    2 => endpoints[1],
                                    _ => basis(axis),
                                });
                            }
                            expected += m.pow(4 - mask.count_ones()) * clifford_trace(&word, false);
                        }
                    }
                    let vector = |argument: AtomView<'_>| -> [i64; 4] {
                        if let Some(i) = component_momenta
                            .iter()
                            .position(|p| p.as_view() == argument)
                        {
                            coordinates[i]
                        } else {
                            let i = component_external
                                .iter()
                                .position(|p| p.as_view() == argument)
                                .unwrap();
                            endpoints[i]
                        }
                    };
                    let actual = component_expression.replace_map(|node, _, out| {
                        if node == mass.as_view() {
                            **out = Atom::num(m);
                        } else if let AtomView::Fun(function) = node {
                            if function.get_symbol() == ETS.metric {
                                let arguments = function.iter().collect::<Vec<_>>();
                                let a = vector(arguments[0]);
                                let b = vector(arguments[1]);
                                **out = Atom::num(a.into_iter().zip(b).map(|(a,b)| a*b).sum::<i64>());
                            } else if let Some(i) = component_momenta.iter().position(|p| {
                                matches!(p.as_view(), AtomView::Fun(p) if p.get_symbol() == function.get_symbol())
                            })
                                && let Some(argument) = function.iter().next()
                                && let Some(axis) = component_external.iter().position(|p| p.as_view() == argument) {
                                **out = Atom::num(coordinates[i][[first, second][axis]]);
                            }
                        }
                    });
                    assert_eq!(
                        actual,
                        Atom::num(expected),
                        "mass={m}, indices=({first},{second})"
                    );
                }
            }
        }
    }
}

#[test]
fn gamma_collection_consumes_incident_coefficients_once_and_keeps_spectators() {
    use crate::tensor::{AlgebraSettings, ReductionStatus};

    let r = test_initialize();
    let label = spenso::index_symbol!("idenso::gamma_coefficient_consumption::edge");
    let [left, right] = [1, 2].map(|i| r.bis4.to_symbolic([symbolica::function!(label, i, 7)]));
    let [mu, alpha] = [3, 4].map(|i| r.mink_d.to_symbolic([symbolica::function!(label, i, 7)]));
    let p = T.rank_one_tensor_symbol("idenso::gamma_coefficient_consumption::p");
    let q = T.rank_one_tensor_symbol("idenso::gamma_coefficient_consumption::q");
    let vector = |head, port: &Atom| FunctionBuilder::new(head).add_arg(port).finish();
    let foreign =
        FunctionBuilder::new(T.tensor_symbol("idenso::gamma_coefficient_consumption::foreign"))
            .add_arg(r.coad_da.to_symbolic([symbolica::function!(label, 5, 7)]))
            .finish();
    let scalar =
        symbolica::parse_lit!((coefficient_consumption_x + coefficient_consumption_y) ^ 12);
    let spectator = &scalar * &foreign;
    let source = SymbolicTensor::infer(
        &spectator
            * (Atom::num(2) * vector(p, &alpha) - Atom::num(3) * vector(q, &alpha))
            * gamma!(&left, &right, &alpha)
            * gamma!(&right, &left, &mu),
    )
    .unwrap();
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        ..Default::default()
    };
    let reduced = source.simplify_algebra(&settings).unwrap();
    assert_eq!(reduced.structure, source.structure);
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    let again = reduced.simplify_algebra(&settings).unwrap();
    assert_eq!(again.expression, reduced.expression);
    assert_eq!(again.structure, reduced.structure);

    // Tr(gamma_alpha gamma_mu) = 4 g_alpha_mu. Compare its tiny selected
    // coefficient polynomial exactly while keeping both spectators opaque.
    assert!(
        reduced
            .expression
            .replace(scalar.to_pattern())
            .with(Atom::Zero)
            .is_zero()
    );
    assert!(
        reduced
            .expression
            .replace(foreign.to_pattern())
            .with(Atom::Zero)
            .is_zero()
    );
    let actual = reduced
        .expression
        .replace(scalar.to_pattern())
        .with(Atom::num(1).to_pattern())
        .replace(foreign.to_pattern())
        .with(Atom::num(1).to_pattern());
    let expected = Atom::num(4) * (Atom::num(2) * vector(p, &mu) - Atom::num(3) * vector(q, &mu));
    crate::test_support::assert_factored_snapshot_eq(&actual.to_string(), &expected.to_string());
}

#[test]
fn repeated_weighted_routes_with_distinct_dummy_ports_preserve_components() {
    use super::super::trace_kernel::tests::clifford_trace;
    use crate::tensor::{AlgebraSettings, ReductionStatus};

    let r = test_initialize();
    let edge = spenso::index_symbol!("idenso::repeated_weighted_route::edge");
    let spin: [Atom; 4] = std::array::from_fn(|i| {
        r.bis4
            .to_symbolic([symbolica::function!(edge, i as i64, 11)])
    });
    let p = T.rank_one_tensor_symbol("idenso::repeated_weighted_route::p");
    let q = T.rank_one_tensor_symbol("idenso::repeated_weighted_route::q");
    let component = |head, port: &Atom| FunctionBuilder::new(head).add_arg(port).finish();
    let route = |port: &Atom| Atom::num(2) * component(p, port) - component(q, port);
    let scalar = symbolica::parse_lit!((weighted_route_x + weighted_route_y) ^ 12);
    let coordinates = [[2_i64, 1, -1, 0], [1_i64, -2, 0, 1]];
    let routed: [i64; 4] = std::array::from_fn(|i| 2 * coordinates[0][i] - coordinates[1][i]);
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        ..Default::default()
    };

    for mink in [&r.mink4, &r.mink_d] {
        let [mu, nu, rho, sigma] =
            [20, 10, 30, 40].map(|i| mink.to_symbolic([symbolica::function!(edge, i, 11)]));
        // The two copies have identical vector bodies but distinct dummy labels.
        // Neither vector sum nor the unrelated scalar sum is distributed here.
        let source = SymbolicTensor::infer(
            &scalar
                * route(&rho)
                * route(&sigma)
                * gamma!(&spin[0], &spin[1], &rho)
                * gamma!(&spin[1], &spin[2], &mu)
                * gamma!(&spin[2], &spin[3], &sigma)
                * gamma!(&spin[3], &spin[0], &nu),
        )
        .unwrap();
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
        assert!(
            result
                .expression
                .replace(scalar.to_pattern())
                .with(Atom::Zero)
                .is_zero()
        );

        let dot = |a: usize, b: usize| {
            coordinates[a]
                .into_iter()
                .zip(coordinates[b])
                .map(|(a, b)| a * b)
                .sum::<i64>()
        };
        let compact = mink.to_symbolic([]);
        let dimension = minkowski_dimension(compact.as_view()).unwrap();
        let value = result
            .expression
            .replace(scalar.to_pattern())
            .with(Atom::num(1).to_pattern());
        for first in 0..4 {
            for second in 0..4 {
                let basis = |axis| std::array::from_fn(|i| i64::from(i == axis));
                let expected =
                    clifford_trace(&[routed, basis(first), routed, basis(second)], false);
                let mut actual = value.clone();
                for (a, head) in [p, q].into_iter().enumerate() {
                    for (port, coordinate) in [(&mu, first), (&nu, second)] {
                        actual = actual
                            .replace(component(head, port).to_pattern())
                            .with(Atom::num(coordinates[a][coordinate]).to_pattern());
                    }
                    for (b, other) in [p, q].into_iter().enumerate().take(a + 1) {
                        actual = actual
                            .replace(
                                g!(component(head, &compact), component(other, &compact))
                                    .to_pattern(),
                            )
                            .with(Atom::num(dot(a, b)).to_pattern());
                    }
                }
                actual = actual
                    .replace(g!(&mu, &nu).to_pattern())
                    .with(Atom::num(i64::from(first == second)).to_pattern())
                    .replace(dimension.to_pattern())
                    .with(Atom::num(4).to_pattern());
                assert_eq!(actual, Atom::num(expected), "indices=({first},{second})");
            }
        }
    }
}
