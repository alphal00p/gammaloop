use super::*;
use crate::tensor::{AlgebraSettings, ReductionStatus};

#[test]
fn massive_routed_trace_collects_coefficients_without_opening_spectators() {
    use super::super::trace_kernel::tests::clifford_trace;

    let reps = test_initialize();
    let label = spenso::index_symbol!("idenso::collected_trace_test::edge");
    let mass = Atom::var(symbolica::symbol!("idenso::collected_trace_test::mass"));
    let p = T.rank_one_tensor_symbol("idenso::collected_trace_test::p");
    let q = T.rank_one_tensor_symbol("idenso::collected_trace_test::q");
    let component = |head, port: &Atom| FunctionBuilder::new(head).add_arg(port).finish();
    let scalar = symbolica::parse_lit!((collected_trace_x + collected_trace_y) ^ 12);
    let foreign = FunctionBuilder::new(T.tensor_symbol("idenso::collected_trace_test::foreign"))
        .add_arg(symbolica::parse_lit!(
            (collected_metadata_x + collected_metadata_y) ^ 5
        ))
        .add_arg(
            reps.coad_da
                .to_symbolic([symbolica::function!(label, 17, 3)]),
        )
        .finish();
    let spin: [Atom; 12] = std::array::from_fn(|i| {
        reps.bis4
            .to_symbolic([symbolica::function!(label, i as i64, 2)])
    });
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        ..Default::default()
    };
    assert!(settings.collect_coefficients);
    let uncollected_settings = AlgebraSettings {
        collect_coefficients: false,
        ..settings.clone()
    };
    let coordinates = [[2_i64, 1, -1, 0], [1_i64, -2, 0, 1]];
    let routed: [[i64; 4]; 2] = [
        std::array::from_fn(|i| 2 * coordinates[0][i] - coordinates[1][i]),
        std::array::from_fn(|i| coordinates[0][i] + 3 * coordinates[1][i]),
    ];

    for mink in [&reps.mink4, &reps.mink_d] {
        let compact = mink.to_symbolic([]);
        let dimension = minkowski_dimension(compact.as_view()).unwrap();
        let [mu, nu, alpha, beta] =
            [20, 10, 30, 40].map(|i| mink.to_symbolic([symbolica::function!(label, i, 7)]));
        let vertices = [&mu, &alpha, &beta, &nu, &alpha, &beta];
        let mut factors = vec![scalar.clone(), foreign.clone()];
        for (i, vertex) in vertices.into_iter().enumerate() {
            // Equal routed vectors belong to separate propagators and use
            // distinct compound dummy labels at every occurrence.
            let rho = mink.to_symbolic([symbolica::function!(label, 100 + i as i64, 7)]);
            let route = if i % 2 == 0 {
                Atom::num(2) * component(p, &rho) - component(q, &rho)
            } else {
                component(p, &rho) + Atom::num(3) * component(q, &rho)
            };
            factors.push(
                &mass * g!(&spin[2 * i], &spin[2 * i + 1])
                    + route * gamma!(&spin[2 * i], &spin[2 * i + 1], &rho),
            );
            factors.push(gamma!(&spin[2 * i + 1], &spin[(2 * i + 2) % 12], (vertex)));
        }
        let source = SymbolicTensor::infer(Atom::mul_many(factors)).unwrap();
        let result = source.simplify_algebra(&settings).unwrap();
        let uncollected = source.simplify_algebra(&uncollected_settings).unwrap();
        for (value, policy) in [(&result, &settings), (&uncollected, &uncollected_settings)] {
            assert_eq!(value.structure, source.structure);
            assert_eq!(value.reduction_status(), ReductionStatus::Complete);
            assert_eq!(value.simplify_algebra(policy).unwrap(), *value);
            // Match whole spectator subtrees: neither their arithmetic nor the
            // tensor's scalar metadata may be distributed by trace collection.
            for spectator in [&scalar, &foreign] {
                assert!(
                    value
                        .expression
                        .replace(spectator.to_pattern())
                        .with(Atom::Zero)
                        .is_zero()
                );
            }
        }
        let compact_bytes = result.expression.as_view().get_byte_size();
        let nested_bytes = uncollected.expression.as_view().get_byte_size();
        assert!(
            compact_bytes * 4 < nested_bytes,
            "generated coefficients were not collected: {compact_bytes} versus {nested_bytes} bytes"
        );

        let values = [&result, &uncollected].map(|value| {
            value
                .expression
                .replace(scalar.to_pattern())
                .with(Atom::num(1).to_pattern())
                .replace(foreign.to_pattern())
                .with(Atom::num(1).to_pattern())
        });
        let dot = |a: usize, b: usize| {
            coordinates[a]
                .into_iter()
                .zip(coordinates[b])
                .map(|(a, b)| a * b)
                .sum::<i64>()
        };
        for m in [0_i64, 2] {
            for first in 0..4 {
                for second in 0..4 {
                    let basis = |axis| std::array::from_fn(|i| i64::from(i == axis));
                    let mut expected = 0_i64;
                    // Only the finite independent oracle enumerates scalar vs
                    // vector matrix entries; the symbolic input stays factored.
                    for mask in 0_u32..64 {
                        for assignment in 0..16 {
                            let vertices = [
                                basis(first),
                                basis(assignment % 4),
                                basis(assignment / 4),
                                basis(second),
                                basis(assignment % 4),
                                basis(assignment / 4),
                            ];
                            let mut word = Vec::new();
                            for (i, vertex) in vertices.into_iter().enumerate() {
                                if mask & (1 << i) != 0 {
                                    word.push(routed[i % 2]);
                                }
                                word.push(vertex);
                            }
                            expected += m.pow(6 - mask.count_ones()) * clifford_trace(&word, false);
                        }
                    }
                    for value in &values {
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
                            .replace(mass.to_pattern())
                            .with(Atom::num(m).to_pattern())
                            .replace(dimension.to_pattern())
                            .with(Atom::num(4).to_pattern());
                        assert_eq!(
                            actual,
                            Atom::num(expected),
                            "mass={m}, ports=({first},{second})"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn cancellation_between_routed_trace_sectors_preserves_typed_zero_ports() {
    let reps = test_initialize();
    let label = spenso::index_symbol!("idenso::collected_trace_zero::edge");
    let p = T.rank_one_tensor_symbol("idenso::collected_trace_zero::p");
    let q = T.rank_one_tensor_symbol("idenso::collected_trace_zero::q");
    let component = |head, port: &Atom| FunctionBuilder::new(head).add_arg(port).finish();
    let route = |port: &Atom| Atom::num(2) * component(p, port) - component(q, port);
    let scalar = symbolica::parse_lit!((collected_zero_x + collected_zero_y) ^ 12);
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        ..Default::default()
    };

    for mink in [&reps.mink4, &reps.mink_d] {
        let [mu, nu, alpha, rho, sigma, tau, eta, zeta] = [20, 10, 30, 40, 50, 60, 70, 80]
            .map(|i| mink.to_symbolic([symbolica::function!(label, i, 5)]));
        let dimension = minkowski_dimension(mink.to_symbolic([]).as_view())
            .unwrap()
            .to_owned();
        let word = |ports: &[&Atom], copy: i64| {
            Atom::mul_many(ports.iter().enumerate().map(|(i, port)| {
                let left = reps
                    .bis4
                    .to_symbolic([symbolica::function!(label, i as i64, copy)]);
                let right = reps.bis4.to_symbolic([symbolica::function!(
                    label,
                    ((i + 1) % ports.len()) as i64,
                    copy
                )]);
                gamma!(&left, &right, *port)
            }))
        };
        // gamma_alpha r_slash gamma_mu gamma^alpha
        // = 4 r_mu + (D - 4) r_slash gamma_mu. Each term has its own
        // spinor and vector dummy scope, while mu and nu retain their identities.
        let six = route(&rho) * route(&sigma) * word(&[&alpha, &rho, &mu, &alpha, &sigma, &nu], 11);
        let four = route(&tau) * route(&eta) * word(&[&tau, &mu, &eta, &nu], 12);
        let two = route(&zeta) * word(&[&zeta, &nu], 13);
        let source = SymbolicTensor::infer(
            &scalar * (six - (dimension - Atom::num(4)) * four - Atom::num(4) * route(&mu) * two),
        )
        .unwrap();
        assert!(!source.expression.is_zero());
        let result = source.simplify_algebra(&settings).unwrap();
        assert!(
            result.expression.is_zero(),
            "uncancelled result: {}",
            result.expression
        );
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
    }
}

#[test]
fn collected_scalar_weights_keep_opaque_powers_and_rounded_arithmetic() {
    use super::super::trace_kernel::tests::clifford_trace;

    let reps = test_initialize();
    let label = spenso::index_symbol!("idenso::collected_trace_weight::edge");
    let [mu, nu, rho, sigma] =
        [20, 10, 30, 40].map(|i| reps.mink_d.to_symbolic([symbolica::function!(label, i, 7)]));
    let spin: [Atom; 4] = std::array::from_fn(|i| {
        reps.bis4
            .to_symbolic([symbolica::function!(label, i as i64, 2)])
    });
    let p = T.rank_one_tensor_symbol("idenso::collected_trace_weight::p");
    let q = T.rank_one_tensor_symbol("idenso::collected_trace_weight::q");
    let component = |head, port: &Atom| FunctionBuilder::new(head).add_arg(port).finish();
    let base = symbolica::parse_lit!(collected_weight_x + collected_weight_y);
    let squared_weight = base.clone().pow(20);
    let rounded = Atom::num(symbolica::domains::float::Float::parse("0.1", Some(11)).unwrap());
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        ..Default::default()
    };
    let uncollected_settings = AlgebraSettings {
        collect_coefficients: false,
        ..settings.clone()
    };

    for weight in [base.clone().pow(10), rounded.clone()] {
        let route = |port: &Atom| &weight * component(p, port) + component(q, port);
        let source = SymbolicTensor::infer(
            route(&rho)
                * route(&sigma)
                * gamma!(&spin[0], &spin[1], &rho)
                * gamma!(&spin[1], &spin[2], &mu)
                * gamma!(&spin[2], &spin[3], &sigma)
                * gamma!(&spin[3], &spin[0], &nu),
        )
        .unwrap();
        let mut source_has_squared_weight = false;
        source.expression.visitor(&mut |node| {
            source_has_squared_weight |= node == squared_weight.as_view();
            true
        });
        assert!(!source_has_squared_weight);
        let result = source.simplify_algebra(&settings).unwrap();
        let uncollected = source.simplify_algebra(&uncollected_settings).unwrap();
        assert_eq!(result.structure, source.structure);
        if weight == rounded {
            // Polynomial reassociation must not turn rounded arithmetic into
            // exact rational arithmetic or change its operation ordering.
            // The existing nonrational kernel may defer structural work; the
            // collection option must preserve that exact result and status.
            assert_eq!(uncollected.structure, source.structure);
            assert_eq!(result.expression, uncollected.expression);
            assert_eq!(result.reduction_status(), uncollected.reduction_status());
            continue;
        }
        assert_eq!(
            result.reduction_status(),
            ReductionStatus::Complete,
            "weight={weight}, result={}",
            result.expression
        );
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);

        let mut has_base = false;
        let mut has_squared_weight = false;
        result.expression.visitor(&mut |node| {
            has_base |= node == base.as_view();
            has_squared_weight |= node == squared_weight.as_view();
            true
        });
        assert!(has_base, "the scalar weight's sum must remain opaque");
        assert!(
            has_squared_weight,
            "the newly generated power must retain its base"
        );

        let compact = reps.mink_d.to_symbolic([]);
        let coordinates = [[2_i64, 1, -1, 0], [1_i64, -2, 0, 1]];
        for scalar in [0_i64, 1, 2] {
            let routed: [i64; 4] =
                std::array::from_fn(|i| scalar.pow(10) * coordinates[0][i] + coordinates[1][i]);
            for first in 0..4 {
                for second in 0..4 {
                    let basis = |axis| std::array::from_fn(|i| i64::from(i == axis));
                    let expected =
                        clifford_trace(&[routed, basis(first), routed, basis(second)], false);
                    let mut actual = result
                        .expression
                        .replace(base.to_pattern())
                        .with(Atom::num(scalar).to_pattern());
                    for (a, head) in [p, q].into_iter().enumerate() {
                        for (port, coordinate) in [(&mu, first), (&nu, second)] {
                            actual = actual
                                .replace(component(head, port).to_pattern())
                                .with(Atom::num(coordinates[a][coordinate]).to_pattern());
                        }
                        for (b, other) in [p, q].into_iter().enumerate().take(a + 1) {
                            let dot = coordinates[a]
                                .into_iter()
                                .zip(coordinates[b])
                                .map(|(a, b)| a * b)
                                .sum::<i64>();
                            actual = actual
                                .replace(
                                    g!(component(head, &compact), component(other, &compact))
                                        .to_pattern(),
                                )
                                .with(Atom::num(dot).to_pattern());
                        }
                    }
                    actual = actual
                        .replace(g!(&mu, &nu).to_pattern())
                        .with(Atom::num(i64::from(first == second)).to_pattern());
                    assert_eq!(
                        actual,
                        Atom::num(expected),
                        "weight={scalar}^10, ports=({first},{second})"
                    );
                }
            }
        }
    }
}

#[test]
fn collected_trace_coefficients_distribute_numerical_factors() {
    let reps = test_initialize();
    let label = spenso::index_symbol!("idenso::collected_trace_expand_num::edge");
    let [mu, rho, sigma, tau] =
        [1, 10, 11, 12].map(|i| reps.mink_d.to_symbolic([symbolica::function!(label, i, 7)]));
    let p = T.rank_one_tensor_symbol("idenso::collected_trace_expand_num::p");
    let component = |port: &Atom| FunctionBuilder::new(p).add_arg(port).finish();
    let word = |port: &Atom, copy: i64| {
        let left = reps
            .bis4
            .to_symbolic([symbolica::function!(label, 0, copy)]);
        let right = reps
            .bis4
            .to_symbolic([symbolica::function!(label, 1, copy)]);
        component(port) * gamma!(&left, &right, (port)) * gamma!(&right, &left, &mu)
    };
    let x = Atom::var(symbolica::symbol!("idenso::collected_trace_expand_num::x"));
    let y = Atom::var(symbolica::symbol!("idenso::collected_trace_expand_num::y"));
    let sum = &x + &y;
    let spectator = symbolica::parse_lit!((collected_expand_num_a + collected_expand_num_b) ^ 12);
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        ..Default::default()
    };
    let uncollected_settings = AlgebraSettings {
        collect_coefficients: false,
        ..settings.clone()
    };

    for remainder in [0, 1] {
        let source = SymbolicTensor::infer(
            &spectator
                * ((Atom::num(2) * &sum + Atom::num(remainder)) * word(&rho, 11)
                    - Atom::num(2) * &x * word(&sigma, 12)
                    - Atom::num(2) * &y * word(&tau, 13)),
        )
        .unwrap();
        let result = source.simplify_algebra(&settings).unwrap();
        let uncollected = source.simplify_algebra(&uncollected_settings).unwrap();
        // The independent two-gamma identity gives 4 p_mu in each sector.
        // Only distributing the numerical coefficient cancels
        // 2*(x+y)-2*x-2*y; products of sums and spectator powers stay intact.
        let expected = Atom::num(4 * remainder) * &spectator * component(&mu);
        assert_eq!(result.expression, expected);
        assert_ne!(uncollected.expression, expected);
        assert_eq!(uncollected.expression.expand_num(), expected);
        for value in [&result, &uncollected] {
            assert_eq!(value.structure, source.structure);
            assert_eq!(value.reduction_status(), ReductionStatus::Complete);
        }
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
        if remainder != 0 {
            assert!(
                result
                    .expression
                    .replace(spectator.to_pattern())
                    .with(Atom::Zero)
                    .is_zero()
            );
        }
    }
}
