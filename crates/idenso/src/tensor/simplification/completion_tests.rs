use super::*;
use crate::color::CS;
use crate::representations::{ColorAdjoint, ColorFundamental};
use crate::test_support;
use crate::{color_f, color_t, gamma};
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::partial::PartialStructureExt;
use spenso::structure::representation::{Minkowski, RepName};
use spenso::{chain, dind, g, mink, p, q, tensor_symbol, trace};
use std::array;
use std::sync::atomic::{AtomicUsize, Ordering};
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder};
use symbolica::{function, parse_lit, symbol};

#[test]
fn alternating_completed_policies_reuse_the_shared_carrier_and_contraction_proof() {
    crate::test_support::test_initialize();
    let value =
        SymbolicTensor::infer((Atom::one() + p!(mink!(4, 99501)) * q!(mink!(4, 99501))).pow(3))
            .unwrap();
    let contracted = value.contract(Default::default()).unwrap();
    let settings = AlgebraSettings::hep();
    let completed = contracted.simplify_algebra(&settings).unwrap();
    assert_eq!(completed, contracted);
    assert!(completed.contraction_complete());

    for _ in 0..3 {
        let structural = completed.contract(Default::default()).unwrap();
        assert!(Arc::ptr_eq(
            completed.reduction_observations(),
            structural.reduction_observations()
        ));
        let algebra = structural.simplify_algebra(&settings).unwrap();
        assert!(Arc::ptr_eq(
            completed.reduction_observations(),
            algebra.reduction_observations()
        ));
        assert_eq!(algebra.reduction_status(), ReductionStatus::Complete);
    }

    // A newly admitted expression does not inherit completion facts.
    let mapped = SymbolicTensor::infer(completed.expression.clone()).unwrap();
    assert!(mapped.proofs.reduction.is_none());
    assert!(!mapped.contraction_complete());
}

#[test]
fn expression_changes_invalidate_previously_completed_policies() {
    let reps = crate::test_support::test_initialize();
    let mu = mink!(4, 99511);
    let nu = mink!(4, 99512);
    let expression = trace!(reps.bis4.to_symbolic([]); [gamma!(&mu), gamma!(&nu)]);
    let value = SymbolicTensor::infer(expression).unwrap();
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let structural_policy = ReductionRequest::Contract(Default::default()).policy();
    let contracted = value.contract(Default::default()).unwrap();
    assert!(contracted.contraction_complete());
    let reduced = contracted.simplify_algebra(&settings).unwrap();
    assert_eq!(reduced.expression.expand(), Atom::num(4) * g!(&mu, &nu));
    assert!(!reduced.contraction_complete());
    assert!(
        !reduced
            .proofs
            .reduction
            .as_ref()
            .unwrap()
            .completed
            .contains(&structural_policy)
    );
    let completed = reduced.contract(Default::default()).unwrap();
    assert_eq!(reduced, completed);
    assert!(completed.contraction_complete());
}

#[test]
fn capped_policy_is_revisited_without_losing_other_completion_facts() {
    let reps = crate::test_support::test_initialize();
    let body = SymbolicTensor::infer(
        trace!(reps.bis4.to_symbolic([]); [gamma!(mink!(4, 99521)), gamma!(mink!(4, 99522))]),
    )
    .unwrap();
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        max_passes: Some(0),
        ..Default::default()
    };
    let disabled = AlgebraSettings {
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let completed = body.simplify_algebra(&disabled).unwrap();
    let capped = completed.simplify_algebra(&settings).unwrap();
    let repeated = capped.simplify_algebra(&settings).unwrap();
    assert_eq!(capped.reduction_status(), ReductionStatus::Capped);
    assert_eq!(repeated.reduction_status(), ReductionStatus::Capped);
    assert_eq!(repeated, completed);
    assert_eq!(
        repeated.proofs.reduction.as_ref().unwrap().completed.len(),
        1
    );
    let restored = repeated.simplify_algebra(&disabled).unwrap();
    assert_eq!(restored.reduction_status(), ReductionStatus::Complete);
    assert!(Arc::ptr_eq(
        restored.reduction_observations(),
        restored
            .simplify_algebra(&disabled)
            .unwrap()
            .reduction_observations()
    ));
}

#[test]
fn deferred_policy_is_revisited_without_losing_other_completion_facts() {
    crate::test_support::test_initialize();
    let port = mink!(4, 99531);
    let r = spenso::vector_symbol!("completion_cache_inexact::r");
    let t = spenso::tensor_symbol!("completion_cache_inexact::t");
    let u = spenso::tensor_symbol!("completion_cache_inexact::u");
    let rounded = Atom::num(symbolica::domains::float::Float::parse("0.1", Some(11)).unwrap());
    // Neither sum offers a source in every term: contraction stays deferred.
    let expression = (rounded * p!(&port) + symbolica::function!(t, &port))
        * (symbolica::function!(r, &port) + symbolica::function!(u, &port));
    let source = SymbolicTensor::infer(expression.clone()).unwrap();
    let disabled = AlgebraSettings {
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let settings = ContractSettings {
        collect_chains: false,
        collect_traces: false,
        ..Default::default()
    };
    let completed = source.simplify_algebra(&disabled).unwrap();
    let deferred = completed.contract(settings).unwrap();
    let repeated = deferred.contract(settings).unwrap();
    assert_eq!(deferred.reduction_status(), ReductionStatus::Deferred);
    assert_eq!(repeated.reduction_status(), ReductionStatus::Deferred);
    assert_eq!(repeated.expression, expression);
    assert_eq!(
        repeated.proofs.reduction.as_ref().unwrap().completed.len(),
        1
    );
    let restored = repeated.simplify_algebra(&disabled).unwrap();
    assert_eq!(restored.reduction_status(), ReductionStatus::Complete);
    assert!(Arc::ptr_eq(
        restored.reduction_observations(),
        restored
            .simplify_algebra(&disabled)
            .unwrap()
            .reduction_observations()
    ));
}

#[test]
fn higher_symmetric_color_invariants_are_stable_completed_outputs() {
    let reps = crate::test_support::test_initialize();
    let representation = reps.cof_nc.to_symbolic([]);
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    for indices in [
        vec![99601, 99602, 99603, 99604, 99605, 99606],
        // A higher contracted symmetric invariant is an exact terminal too;
        // completion does not claim a formula solely in quadratic Casimirs.
        vec![99601, 99601, 99603, 99604, 99605, 99606],
    ] {
        let expression = crate::color::CS.symmetric_generator_trace(
            &representation,
            indices
                .into_iter()
                .map(|index| reps.coad_da.to_symbolic([Atom::num(index)])),
        );
        let source = SymbolicTensor::infer(expression.clone()).unwrap();
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(result.expression, expression);
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        let repeated = result.simplify_algebra(&settings).unwrap();
        assert_eq!(repeated, result);
        assert_eq!(repeated.reduction_status(), ReductionStatus::Complete);
        assert!(Arc::ptr_eq(
            result.reduction_observations(),
            repeated.reduction_observations()
        ));
    }
}

#[test]
fn a_higher_color_gram_invariant_remains_symbolic_and_complete() {
    let reps = crate::test_support::test_initialize();
    let representation = reps.cof_nc.to_symbolic([]);
    let symmetric = crate::color::CS.symmetric_generator_trace(
        &representation,
        (99611..99616).map(|index| reps.coad_da.to_symbolic([Atom::num(index)])),
    );
    let source = SymbolicTensor::infer(symmetric.pow(2)).unwrap();
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let result = source.simplify_algebra(&settings).unwrap();
    assert_eq!(
        result.expression,
        crate::color::CS.gram(Atom::num(5), &representation, &representation)
    );
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
}

#[test]
fn squared_symmetric_traces_keep_internal_dummy_pairs_in_each_copy() {
    let reps = crate::test_support::test_initialize();
    let representation = reps.cof_nc.to_symbolic([]);
    let symmetric = crate::color::CS.symmetric_generator_trace(
        &representation,
        [99631, 99631, 99632, 99633, 99634, 99635]
            .into_iter()
            .map(|index| reps.coad_da.to_symbolic([Atom::num(index)])),
    );
    let one_copy = SymbolicTensor::infer(symmetric.clone()).unwrap();
    assert_eq!(one_copy.structure.slots().unwrap().len(), 4);
    let source = SymbolicTensor::infer(symmetric.pow(2)).unwrap();
    assert!(source.structure.slots().unwrap().is_empty());
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let result = source.simplify_algebra(&settings).unwrap();
    // Only the four external slots connect the two copies. Treating all six
    // slots as shared would incorrectly replace two internal traces by gram(6).
    assert_eq!(result.expression, source.expression);
    assert_ne!(
        result.expression,
        crate::color::CS.gram(Atom::num(6), &representation, &representation)
    );
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
}

#[test]
fn mixed_rank_six_gram_contraction_preserves_the_raw_adjoint_phase() {
    let reps = crate::test_support::test_initialize();
    let fundamental = reps.cof_nc.to_symbolic([]);
    let adjoint = reps.coad_da.to_symbolic([]);
    let ports = (99641..99647)
        .map(|index| reps.coad_da.to_symbolic([Atom::num(index)]))
        .collect::<Vec<_>>();
    let fundamental_trace = crate::color::CS.symmetric_generator_trace(&fundamental, &ports);
    let raw_adjoint_trace = spenso::trace_sym!(&adjoint; ports.iter().map(|port| {
        crate::color_f!(
            Atom::var(spenso::network::tags::SPENSO_TAG.chain_in),
            Atom::var(spenso::network::tags::SPENSO_TAG.chain_out),
            port
        )
    }));
    let source = SymbolicTensor::infer(fundamental_trace * raw_adjoint_trace).unwrap();
    assert!(source.structure.slots().unwrap().is_empty());
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let result = source.simplify_algebra(&settings).unwrap();
    // Raw adjoint matrices are i times the Hermitian generators, so their
    // symmetric word contributes i^6 = -1 to the normalized invariant.
    assert_eq!(
        result.expression,
        -crate::color::CS.gram(Atom::num(6), &adjoint, &fundamental)
    );
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
}

#[test]
fn color_only_reduction_is_independent_of_cooked_labels_and_spectators() {
    crate::test_support::test_initialize();
    let settings = [
        ColorSimplifySettings::default(),
        ColorSimplifySettings::default().without_trace_evaluation(),
        ColorSimplifySettings::default().with_cof_dimension_invariants(),
        ColorSimplifySettings::default()
            .with_cof_dimension_invariants()
            .without_trace_evaluation(),
    ]
    .into_iter()
    .flat_map(|color| {
        [false, true].map(|collect_coefficients| AlgebraSettings {
            color: Some(color),
            contract: AlgebraContraction::None,
            collect_coefficients,
            ..Default::default()
        })
    })
    .collect::<Vec<_>>();
    let cooking = crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
    let scalar = symbolica::parse_lit!((color_label_x + color_label_y) ^ 30);
    let momentum = spenso::network::tags::SPENSO_TAG.rank_one_tensor_symbol("Q");
    let mut failures = Vec::new();
    for shift in 0..6 {
        for reversed in [false, true] {
            // Fresh function-valued labels are cooked in each term's own order,
            // just like the half-edge labels admitted by the generator.
            let offset = 1000 + 10 * (2 * shift + i64::from(reversed));
            let labels = [0, 1, 2, 3, 4, 5].map(|position| {
                let position = (position + shift) % 6;
                function!(
                    symbol!("completion_color_labels::label"),
                    offset + if reversed { 5 - position } else { position }
                )
            });
            let [a, b, c, d, e] = [0, 1, 2, 3, 4].map(|position| {
                ColorAdjoint {}.to_symbolic([Atom::num(8), labels[position].clone()])
            });
            let mu = Minkowski {}.to_symbolic([Atom::num(4), labels[5].clone()]);
            let q = |edge: i64| {
                FunctionBuilder::new(momentum)
                    .add_arg(edge)
                    .add_arg(&mu)
                    .finish()
            };
            let bubble = crate::color_f!(&a, &b, &c) * crate::color_f!(&a, &b, &d);
            for (spectator_kind, spectator) in
                [("scalar", scalar.clone()), ("indexed", q(0) * q(2))]
            {
                for (case, color, expected) in [
                    ("zero", &bubble * crate::color_f!(&c, &d, &e), Atom::Zero),
                    ("nonzero", bubble.clone(), g!(&c, &d) * &spectator),
                ] {
                    let expression = cooking.try_cook((color * &spectator).as_view()).unwrap();
                    let source = SymbolicTensor::infer(expression).unwrap();
                    for settings in &settings {
                        let color_settings = settings.color.unwrap();
                        let casimir = if color_settings.substitute_cof_dimension_invariants {
                            Atom::num(3)
                        } else {
                            crate::color::CS
                                .cas(Atom::num(2), ColorAdjoint {}.to_symbolic([Atom::num(8)]))
                        };
                        let expected = cooking.try_cook((&expected * casimir).as_view()).unwrap();
                        let result = source.simplify_algebra(settings).unwrap();
                        if result.expression != expected
                            || result.structure != source.structure
                            || result.reduction_status() != ReductionStatus::Complete
                        {
                            failures.push(format!(
                                "shift={shift}, reversed={reversed}, spectator={spectator_kind}, \
                                 case={case}, evaluate_traces={}, cof_invariants={}, \
                                 collect_coefficients={}: source={}, result={}, \
                                 expected={}, status={:?}, interface_preserved={}",
                                color_settings.evaluate_traces,
                                color_settings.substitute_cof_dimension_invariants,
                                settings.collect_coefficients,
                                source.expression,
                                result.expression,
                                expected,
                                result.reduction_status(),
                                result.structure == source.structure,
                            ));
                        }
                    }
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// With trace evaluation disabled a symmetric trace is not evaluated, but its
/// row may absorb the metric between its own ports: sym(x,y,z) g(x,y) becomes
/// sym(y,y,z). That form stays exact. Evaluating its traces afterwards gives
/// the result of evaluating them directly, the zero of an invariant with one
/// open adjoint index. Active identities in the same product or sum reduce.
#[test]
fn disabled_trace_evaluation_stays_exact_beside_active_identities() {
    let reps = test_support::test_initialize();
    let [a, b, c, d, x, y, z] =
        array::from_fn(|index| reps.coad_da.to_symbolic([Atom::num(99801 + index as i64)]));
    let symmetric = CS.symmetric_generator_trace(reps.cof_nc.to_symbolic([]), [&x, &y, &z]);
    let retained = &symmetric * g!(&x, &y);
    let ordinary = trace!(reps.cof_nc.to_symbolic([]); [color_t!(&z)]);
    // Each active product, a factor its identity keeps and one it removes.
    let active_cases = [
        (
            g!(&a, &b) * color_f!(&b, &c, &d),
            color_f!(&a, &c, &d),
            g!(&a, &b),
        ),
        (
            color_f!(&a, &b, &c) * color_f!(&a, &b, &d),
            g!(&c, &d),
            color_f!(&a, &b, &c),
        ),
    ];
    let disabled = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().without_trace_evaluation()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let enabled = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        ..disabled.clone()
    };
    let colourless = AlgebraSettings {
        color: None,
        ..disabled.clone()
    };
    for spectator in [
        parse_lit!((local_color_x + local_color_y) ^ 30),
        p!(mink!(4, 99811)) * q!(mink!(4, 99811)),
    ] {
        for (active, kept, removed) in &active_cases {
            for expression in [
                active * &retained * &spectator,
                active * (&retained + &ordinary) * &spectator,
            ] {
                let source = SymbolicTensor::infer(expression.clone()).unwrap();
                assert_eq!(
                    source.simplify_algebra(&colourless).unwrap().expression,
                    expression
                );
                let result = source.simplify_algebra(&disabled).unwrap();
                assert!(
                    result.expression.contains(kept.as_view())
                        && !result.expression.contains(removed.as_view()),
                    "{}",
                    result.expression
                );
                assert_eq!(result.structure, source.structure);
                assert_eq!(result.reduction_status(), ReductionStatus::Complete);
                assert_eq!(result.simplify_algebra(&disabled).unwrap(), result);
                let direct = source.simplify_algebra(&enabled).unwrap();
                assert_eq!(direct.expression, Atom::Zero);
                assert_eq!(
                    result.simplify_algebra(&enabled).unwrap().expression,
                    direct.expression,
                    "{}",
                    result.expression
                );
            }
        }
    }
}

#[test]
fn color_component_selection_keeps_unrelated_sums_factored() {
    let reps = crate::test_support::test_initialize();
    let [a, b, c, x, y, z] =
        std::array::from_fn(|index| reps.coad_da.to_symbolic([Atom::num(99831 + index as i64)]));
    let symmetric =
        crate::color::CS.symmetric_generator_trace(reps.cof_nc.to_symbolic([]), [&x, &y, &z]);
    let retained = symmetric * g!(&x, &y);
    let ordinary = trace!(reps.cof_nc.to_symbolic([]); [crate::color_t!(&z)]);
    let unrelated = crate::color_f!(&a, &b, &c)
        + trace!(reps.cof_nc.to_symbolic([]);
            [crate::color_t!(&a), crate::color_t!(&b), crate::color_t!(&c)]);
    let expression = (retained + ordinary)
        * unrelated
        * symbolica::parse_lit!((unrelated_color_x + unrelated_color_y) ^ 30);
    let source = SymbolicTensor::infer(expression.clone()).unwrap();
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().without_trace_evaluation()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let result = source.simplify_algebra(&settings).unwrap();
    assert_eq!(result.expression, expression);
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
}

/// An identity in one summand reduces that summand in place: the unrelated
/// colour sum and the scalar power stay whole factors, and the disabled
/// symmetric trace in the other summand stays exact.
#[test]
fn disabled_trace_rows_keep_unrelated_sums_factored() {
    let reps = test_support::test_initialize();
    let [a, b, c, x, y, z, p, q, r] =
        array::from_fn(|index| reps.coad_da.to_symbolic([Atom::num(99841 + index as i64)]));
    let symmetric = CS.symmetric_generator_trace(reps.cof_nc.to_symbolic([]), [&x, &y, &z]);
    let retained = symmetric * g!(&x, &y);
    // f(z,p,q) f(r,p,q) = CA g(z,r), so this summand becomes CA Tr(T^z).
    let active = color_f!(&z, &p, &q)
        * color_f!(&r, &p, &q)
        * trace!(reps.cof_nc.to_symbolic([]); [color_t!(&r)]);
    let reduced = CS.cas(Atom::num(2), reps.coad_da.to_symbolic([]))
        * trace!(reps.cof_nc.to_symbolic([]); [color_t!(&z)]);
    let unrelated = color_f!(&a, &b, &c)
        + trace!(reps.cof_nc.to_symbolic([]); [color_t!(&a), color_t!(&b), color_t!(&c)]);
    let scalar = parse_lit!((unrelated_color_x + unrelated_color_y) ^ 30);
    let source = SymbolicTensor::infer((&retained + active) * &unrelated * &scalar).unwrap();
    let disabled = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().without_trace_evaluation()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let result = source.simplify_algebra(&disabled).unwrap();
    let AtomView::Mul(product) = result.expression.as_view() else {
        panic!("{}", result.expression)
    };
    let factors = product.iter().collect::<Vec<_>>();
    assert!(
        factors.contains(&unrelated.as_view()) && factors.contains(&scalar.as_view()),
        "{}",
        result.expression
    );
    assert!(
        result.expression.contains(reduced.as_view()),
        "{}",
        result.expression
    );
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&disabled).unwrap(), result);
    let enabled = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        ..disabled
    };
    assert_eq!(
        result.simplify_algebra(&enabled).unwrap().expression,
        source.simplify_algebra(&enabled).unwrap().expression
    );
}

#[test]
fn four_loop_color_zero_is_independent_of_spectators() {
    crate::test_support::test_initialize();
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let cooking = crate::CookSettings::indices().with_mode(crate::CookMode::ReversibleEncoding);
    // Color graph of D0108 in four_loop_numerator. The same exact filter
    // removes this graph under one index registration order and retains it
    // under another. Both forms below denote the same color product.
    let vertices = [
        [0, 7, 8],
        [1, 2, 5],
        [10, 17, 18],
        [11, 12, 14],
        [13, 21, 6],
        [15, 22, 3],
        [16, 20, 23],
        [19, 4, 9],
    ];
    let connections = [
        [10, 11],
        [12, 13],
        [14, 15],
        [16, 17],
        [18, 19],
        [2, 3],
        [20, 21],
        [22, 23],
        [4, 5],
        [6, 7],
        [8, 9],
    ];
    let mut failures = Vec::new();
    for shift in [0, 1, 7, 13] {
        let slots: [Atom; 24] = std::array::from_fn(|index| {
            ColorAdjoint {}.to_symbolic([
                Atom::num(8),
                function!(
                    symbol!("four_loop_color_test::label"),
                    10_000 + 100 * shift + (index + shift) % 24
                ),
            ])
        });
        for contracted in [false, true] {
            let slot = |index: usize| {
                let index = if contracted {
                    connections
                        .iter()
                        .find_map(|pair| pair.contains(&index).then_some(pair[0].min(pair[1])))
                        .unwrap_or(index)
                } else {
                    index
                };
                &slots[index]
            };
            let mut color = Atom::mul_many(
                vertices.map(|[a, b, c]| crate::color_f!(slot(a), slot(b), slot(c))),
            );
            if !contracted {
                color *= Atom::mul_many(connections.map(|[a, b]| g!(slot(a), slot(b))));
            }
            for (spectator_kind, spectator) in [
                ("none", Atom::one()),
                (
                    "scalar",
                    symbolica::parse_lit!((four_loop_x + four_loop_y) ^ 30),
                ),
                ("indexed", p!(mink!(4, 99_899)) * q!(mink!(4, 99_899))),
            ] {
                let expression = cooking
                    .try_cook((color.clone() * spectator).as_view())
                    .unwrap();
                let source = SymbolicTensor::infer(expression).unwrap();
                let result = source.simplify_algebra(&settings).unwrap();
                if !result.expression.is_zero()
                    || result.structure != source.structure
                    || result.reduction_status() != ReductionStatus::Complete
                {
                    failures.push(format!(
                        "shift={shift}, contracted={contracted}, spectator={spectator_kind}, \
                         status={:?}, result={}",
                        result.reduction_status(),
                        result.expression,
                    ));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// A summand may reuse an outer dummy pair, since a sum is its own scope:
/// f(p,q,r) f(p,q,s) [g(z,x) g(x,y) + f(p,q,z) f(p,q,y)] = CA (1 + CA) g(r,s) g(z,y).
/// The pairs stay apart in every contraction mode, also for labels that
/// interleave them, so contracting with f(r,s,w) gives an exact zero.
#[test]
fn reused_summand_dummies_stay_apart_in_every_contraction_mode() {
    test_support::test_initialize();
    let adjoint = |index: i64| ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(index)]);
    // The labels interleave the outer and the inner pair: r < z < s < y.
    let [p, q, r, z, s, y, x, w] =
        [98101, 98102, 98103, 98104, 98105, 98106, 98107, 98108].map(adjoint);
    let product = color_f!(&p, &q, &r)
        * color_f!(&p, &q, &s)
        * (g!(&z, &x) * g!(&x, &y) + color_f!(&p, &q, &z) * color_f!(&p, &q, &y));
    for contract in [
        AlgebraContraction::None,
        AlgebraContraction::Minimal,
        AlgebraContraction::Fully,
    ] {
        let settings = AlgebraSettings {
            color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
            contract,
            ..Default::default()
        };
        let source = SymbolicTensor::infer(product.clone()).unwrap();
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(
            result.expression,
            Atom::num(12) * g!(&r, &s) * g!(&z, &y),
            "{:?}",
            settings.contract
        );
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        let source = SymbolicTensor::infer(&product * color_f!(&r, &s, &w)).unwrap();
        assert_eq!(
            source.simplify_algebra(&settings).unwrap().expression,
            Atom::Zero,
            "{:?}",
            settings.contract
        );
    }
}

/// Colour deltas have exact rules of their own: a closed quark loop is N, a
/// chain contracts to one delta, and foreign tensors joined by a delta are
/// canonicalized with it. They stay colour work under AlgebraContraction::None
/// even where no generator or structure constant touches them.
#[test]
fn colour_deltas_without_generators_stay_colour_work() {
    test_support::test_initialize();
    let fundamental =
        |index: i64| ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(index)]);
    let adjoint = |index: i64| ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(index)]);
    let numeric = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let symbolic = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        ..numeric.clone()
    };
    let simplify = |expression: Atom, settings: &AlgebraSettings| {
        let result = SymbolicTensor::infer(expression)
            .unwrap()
            .simplify_algebra(settings)
            .unwrap();
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        result.expression
    };
    // δ^i_j δ^j_i = N = 3 and Tr(T^a T^a) = T_F (N² - 1) = 4.
    let [i, j, k, l] = [98201, 98202, 98203, 98204].map(fundamental);
    let a = adjoint(98205);
    let quark_loop = g!(dind!(&i), &j) * g!(dind!(&j), &i);
    let pair = color_t!(&a, &k, dind!(&l)) * color_t!(&a, &l, dind!(&k));
    assert_eq!(simplify(&quark_loop * &pair, &numeric), Atom::num(12));
    assert_eq!(
        simplify(&quark_loop * &pair - Atom::num(12), &numeric),
        Atom::Zero
    );
    // f(b,c,d) f(b,c,e) g(x,y) g(y,z) = CA g(d,e) g(x,z).
    let [b, c, d, e, x, y, z] = [98211, 98212, 98213, 98214, 98215, 98216, 98217].map(adjoint);
    let ca = CS.cas(Atom::num(2), ColorAdjoint {}.to_symbolic([Atom::num(8)]));
    let chain = color_f!(&b, &c, &d) * color_f!(&b, &c, &e) * g!(&x, &y) * g!(&y, &z)
        - ca * g!(&d, &e) * g!(&x, &z);
    assert_eq!(simplify(chain, &symbolic), Atom::Zero);
    // X V(m) g(m,n) W(n) - X' V(u) g(u,v) W(v) = 0 for X = f(p,q,r)^2 = X'.
    let v_head = tensor_symbol!("completion_colour_deltas::V");
    let w_head = tensor_symbol!("completion_colour_deltas::W");
    let [m, n, u, v, p, q, r, s, t, o] = [
        98221, 98222, 98223, 98224, 98225, 98226, 98227, 98228, 98229, 98230,
    ]
    .map(adjoint);
    let joined = |loop_: Atom, left: &Atom, right: &Atom| {
        loop_ * function!(v_head, left) * g!(left, right) * function!(w_head, right)
    };
    let first = joined(color_f!(&p, &q, &r) * color_f!(&p, &q, &r), &m, &n);
    let second = joined(color_f!(&s, &t, &o) * color_f!(&s, &t, &o), &u, &v);
    assert_eq!(simplify(first - second, &numeric), Atom::Zero);
}

/// Invariant identities need no trace evaluation: a closed network of
/// symmetric traces with one open adjoint index is an invariant vector of the
/// adjoint representation, which has none but zero.
#[test]
fn disabled_trace_evaluation_keeps_invariant_identities() {
    let reps = test_support::test_initialize();
    let [a, b, c, e] =
        array::from_fn(|index| reps.coad_da.to_symbolic([Atom::num(98301 + index as i64)]));
    let fundamental = reps.cof_nc.to_symbolic([]);
    let network = CS.symmetric_generator_trace(fundamental.clone(), [&a, &b, &c])
        * CS.symmetric_generator_trace(fundamental, [&a, &b, &c, &e]);
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().without_trace_evaluation()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let result = SymbolicTensor::infer(network)
        .unwrap()
        .simplify_algebra(&settings)
        .unwrap();
    assert_eq!(result.expression, Atom::Zero);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
}

/// An unrelated foreign tensor beside an isolated trace keeps the trace's
/// one-shot certificate: a single colour pass completes the sum.
#[test]
fn foreign_summands_keep_isolated_trace_certificates() {
    test_support::test_initialize();
    let adjoint = |index: i64| ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(index)]);
    let indices = 98801..98808;
    let word = trace!(ColorFundamental {}.to_symbolic([Atom::num(3)]);
        indices.clone().map(|index| color_t!(adjoint(index))));
    let foreign = indices
        .map(adjoint)
        .fold(
            FunctionBuilder::new(tensor_symbol!("completion_isolated_trace::U")),
            |builder, slot| builder.add_arg(&slot),
        )
        .finish();
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::None,
        max_passes: Some(1),
        ..Default::default()
    };
    let adjoint_word = trace!(ColorAdjoint {}.to_symbolic([Atom::num(8)]);
    (98801..98808).map(|index| color_f!(
        Atom::var(SPENSO_TAG.chain_in),
        Atom::var(SPENSO_TAG.chain_out),
        adjoint(index)
    )));
    for word in [word, adjoint_word] {
        let source = SymbolicTensor::infer(word + &foreign).unwrap();
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.structure, source.structure);
        let fresh = SymbolicTensor::infer(result.expression.clone()).unwrap();
        assert_eq!(fresh.simplify_algebra(&settings).unwrap(), result);
    }
}

/// Completion does not depend on the order in which symbols were created.
/// Occurrences competing for a changed coefficient open together, and
/// contexts competing for a term stay as written. Like terms merge whatever
/// their signs, as do equal canonical states published by several rows.
#[test]
fn completed_coefficients_are_independent_of_symbol_registration_order() {
    test_support::test_initialize();
    let fundamental =
        |index: i64| ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(index)]);
    let adjoint = |index: i64| ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(index)]);
    // The crossed deltas give sum_a Tr(T^a T^a) = 4 in SU(3).
    let crossed = |base: i64| {
        let [i, j, k, l] = [0, 1, 2, 3].map(|offset| fundamental(base + offset));
        let a = adjoint(base + 4);
        color_t!(&a, &i, dind!(&j))
            * color_t!(&a, &k, dind!(&l))
            * g!(dind!(&i), &l)
            * g!(dind!(&k), &j)
    };
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let simplify = |expression: Atom, settings: &AlgebraSettings| {
        SymbolicTensor::infer(expression)
            .unwrap()
            .simplify_algebra(settings)
            .unwrap()
            .expression
    };
    let name = |label: &str| symbol!(format!("completion_registration_order::{label}"));
    // b is created before c for the first copy, and c before b for the second.
    let [a1, b1, c1] = ["a1", "b1", "c1"].map(name);
    let [a2, c2, b2] = ["a2", "c2", "b2"].map(name);
    let shared = |[a, b, c]: [Symbol; 3], base: i64| {
        let [a, b, c] = [a, b, c].map(Atom::var);
        -crossed(base) * &a / Atom::num(2) + Atom::num(2) * (&a + &b) + Atom::num(2) * (&a + &c)
    };
    // -2a + 2(a+b) + 2(a+c): both occurrences hold a, so neither owns it.
    let first = simplify(shared([a1, b1, c1], 98401), &settings);
    let second = simplify(shared([a2, b2, c2], 98411), &settings);
    let [a, b, c] = [a1, b1, c1].map(Atom::var);
    assert_eq!(first, (Atom::num(2) * (&a + &b + &c)).expand());
    assert_eq!(
        second
            .replace(a2)
            .with(a.clone())
            .replace(b2)
            .with(b.clone())
            .replace(c2)
            .with(c.clone()),
        first
    );
    // x and y compete for x*y*(..): neither regroups, whichever came first.
    let [x1, y1] = ["x1", "y1"].map(name);
    let [y2, x2] = ["y2", "x2"].map(name);
    let [u, v, t] = ["u", "v", "t"].map(|label| Atom::var(name(label)));
    let contexts = |x: Symbol, y: Symbol, base: i64| {
        let [x, y] = [x, y].map(Atom::var);
        &x * &y * (crossed(base) + &u)
            + &x * (crossed(base + 10) + &v)
            + &y * (crossed(base + 20) + &t)
    };
    let first = simplify(contexts(x1, y1, 98501), &settings);
    let second = simplify(contexts(x2, y2, 98531), &settings);
    let [x, y] = [x1, y1].map(Atom::var);
    let four = Atom::num(4);
    assert_eq!(
        first,
        &x * &y * (&four + &u) + &x * (&four + &v) + &y * (&four + &t)
    );
    assert_eq!(
        second
            .replace(x2)
            .with(x.clone())
            .replace(y2)
            .with(y.clone()),
        first
    );
    // A changed term merges with an untouched term of equal sign.
    assert_eq!(
        simplify(
            crossed(98601) * (&x + &y) + Atom::num(2) * (&x + &y),
            &settings
        ),
        Atom::num(6) * &x + Atom::num(6) * &y
    );
    // Equal canonical states of several rows merge with equal-sign weights
    // and cancel exactly against one equal row. A row may publish its states
    // factored or distributed; the comparison sees through either.
    let symbolic = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        ..settings.clone()
    };
    let representation = ColorFundamental {}.to_symbolic([Atom::num(3)]);
    let word = |dummy: i64| {
        trace!(representation.clone(); [98701, 98702, 98703, 98704, dummy, dummy]
            .map(|index| color_t!(adjoint(index))))
    };
    let rows = word(98711) + word(98712) + Atom::num(2) * word(98713);
    assert_eq!(
        simplify(rows.clone(), &symbolic).expand(),
        simplify(Atom::num(4) * word(98714), &symbolic).expand()
    );
    assert_eq!(
        simplify(rows - Atom::num(4) * word(98714), &symbolic),
        Atom::Zero
    );
}

#[test]
fn fully_contracted_structure_constants_preserve_literal_port_orientation() {
    crate::test_support::test_initialize();
    let adjoint = crate::representations::ColorAdjoint {};
    let representation = adjoint.to_symbolic([Atom::num(8)]);
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    for labels in [
        [
            symbolica::symbol!("completion_color_orientation::a"),
            symbolica::symbol!("completion_color_orientation::b"),
            symbolica::symbol!("completion_color_orientation::c"),
        ],
        [
            symbolica::symbol!("completion_color_orientation::a_"),
            symbolica::symbol!("completion_color_orientation::b_"),
            symbolica::symbol!("completion_color_orientation::c_"),
        ],
    ] {
        let [a, b, c] = labels.map(|label| adjoint.to_symbolic([Atom::num(8), Atom::var(label)]));
        let source =
            SymbolicTensor::infer(crate::color_f!(&a, &b, &c) * crate::color_f!(&b, &a, &c))
                .unwrap();
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(
            result.expression,
            -Atom::num(8) * crate::color::CS.cas(Atom::num(2), &representation)
        );
        assert!(result.structure.slots().unwrap().is_empty());
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
    }
}

#[test]
fn rank_five_adjoint_traces_normalize_transposed_and_signed_factors() {
    let reps = crate::test_support::test_initialize();
    let representation = reps.coad_da.to_symbolic([]);
    let labels = [
        symbolica::symbol!("completion_adjoint_trace::a_"),
        symbolica::symbol!("completion_adjoint_trace::b_"),
        symbolica::symbol!("completion_adjoint_trace::c_"),
        symbolica::symbol!("completion_adjoint_trace::d_"),
        symbolica::symbol!("completion_adjoint_trace::e_"),
    ];
    let ports = labels.map(|label| reps.coad_da.to_symbolic([Atom::var(label)]));
    let incoming = Atom::var(spenso::network::tags::SPENSO_TAG.chain_in);
    let outgoing = Atom::var(spenso::network::tags::SPENSO_TAG.chain_out);
    let factors = ports
        .iter()
        .map(|port| crate::color_f!(&incoming, &outgoing, port))
        .collect::<Vec<_>>();
    let source = SymbolicTensor::infer(trace!(&representation; &factors)).unwrap();
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let ordinary = source.simplify_algebra(&settings).unwrap();
    assert_ne!(ordinary.expression, source.expression);
    assert!(!ordinary.expression.is_zero());
    assert_eq!(ordinary.structure, source.structure);
    assert_eq!(ordinary.reduction_status(), ReductionStatus::Complete);
    let expected = (-ordinary.clone()).canonize().unwrap();

    for first in [
        crate::color_f!(&outgoing, &incoming, &ports[0]),
        -&factors[0],
    ] {
        let mut altered = factors.clone();
        altered[0] = first;
        let altered = SymbolicTensor::infer(trace!(&representation; altered)).unwrap();
        let result = altered.simplify_algebra(&settings).unwrap();
        assert_eq!(result.structure, altered.structure);
        // Independent admission can rotate a cyclic word's logical interface.
        // Align every explicit port position before comparing the tensors.
        let actual_slots = result.structure.logical_slots();
        let axes = expected
            .structure
            .logical_slots()
            .iter()
            .map(|slot| {
                actual_slots
                    .iter()
                    .position(|actual| actual == slot)
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let aligned = result.permuted(&axes).unwrap().canonize().unwrap();
        // F^T=-F. Canonicalization identifies generated dummy names without
        // evaluating trace tensors. The colour kernel returns canonical
        // states as a flat sum, which the negated reference keeps factored.
        assert_eq!(aligned.structure, expected.structure);
        assert!(
            (&aligned.expression - &expected.expression)
                .expand()
                .is_zero(),
            "{}\n{}",
            aligned.expression,
            expected.expression
        );
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
    }
}

#[test]
fn arbitrary_color_trace_completion_respects_disabled_work_and_zero_budget() {
    let reps = crate::test_support::test_initialize();
    let word = trace!(reps.cof_nc.to_symbolic([]);
        (99621..99626).map(|index| crate::color_t!(reps.coad_da.to_symbolic([Atom::num(index)]))));
    let source = SymbolicTensor::infer(word.clone()).unwrap();
    let disabled = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().without_trace_evaluation()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let unchanged = source.simplify_algebra(&disabled).unwrap();
    assert_eq!(unchanged.expression, word);
    assert_eq!(unchanged.reduction_status(), ReductionStatus::Complete);
    assert_eq!(unchanged.simplify_algebra(&disabled).unwrap(), unchanged);

    let enabled = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        ..disabled
    };
    let capped_settings = AlgebraSettings {
        max_passes: Some(0),
        ..enabled.clone()
    };
    let capped = unchanged.simplify_algebra(&capped_settings).unwrap();
    assert_eq!(capped.expression, word);
    assert_eq!(capped.structure, source.structure);
    assert_eq!(capped.reduction_status(), ReductionStatus::Capped);
    let still_capped = capped.simplify_algebra(&capped_settings).unwrap();
    assert_eq!(still_capped.expression, word);
    assert_eq!(still_capped.reduction_status(), ReductionStatus::Capped);

    let reduced = capped.simplify_algebra(&enabled).unwrap();
    assert_ne!(reduced.expression, word);
    assert_eq!(reduced.structure, source.structure);
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    let repeated = reduced.simplify_algebra(&enabled).unwrap();
    assert_eq!(repeated, reduced);
    assert_eq!(repeated.reduction_status(), ReductionStatus::Complete);
    assert!(Arc::ptr_eq(
        reduced.reduction_observations(),
        repeated.reduction_observations()
    ));
}

#[test]
fn single_bridge_color_sums_and_scalar_spectators_remain_factored() {
    let reps = crate::test_support::test_initialize();
    let [a, b, x, c, d, e] =
        std::array::from_fn(|i| reps.coad_da.to_symbolic([Atom::num(99701 + i as i64)]));
    let ports = [&x, &c, &d, &e];
    let fundamental =
        crate::color::CS.symmetric_generator_trace(reps.cof_nc.to_symbolic([]), ports);
    let adjoint = spenso::trace_sym!(reps.coad_da.to_symbolic([]); ports.map(|port| {
        crate::color_f!(
            Atom::var(spenso::network::tags::SPENSO_TAG.chain_in),
            Atom::var(spenso::network::tags::SPENSO_TAG.chain_out),
            port
        )
    }));
    let sum = fundamental + adjoint;
    let scalar = symbolica::parse_lit!((single_bridge_x + single_bridge_y) ^ 12);
    let bridge = crate::color_f!(&a, &b, &x);
    for outside in [&scalar, &bridge, &(&scalar * &bridge)] {
        let source = SymbolicTensor::infer(outside * &sum).unwrap();
        let kernel = crate::color::simplify::ColorAlgebraSimplifier::new(
            ColorSimplifySettings::default(),
            SymbolicTensor::<PartialStructure>::reserved_dummies([&source]),
        );
        assert_eq!(
            kernel.step(source.expression.as_view(), true),
            source.expression
        );
        let settings = AlgebraSettings {
            color: Some(ColorSimplifySettings::default()),
            contract: AlgebraContraction::None,
            ..Default::default()
        };
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(result.expression, source.expression);
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
    }
}

#[test]
fn two_color_connections_still_distribute_a_selected_sum() {
    let reps = crate::test_support::test_initialize();
    let [a, b, c, x] =
        std::array::from_fn(|i| reps.coad_da.to_symbolic([Atom::num(99711 + i as i64)]));
    let symmetric =
        crate::color::CS.symmetric_generator_trace(reps.cof_nc.to_symbolic([]), [&a, &b, &c]);
    let source = SymbolicTensor::infer(
        crate::color_f!(&a, &b, &x) * (crate::color_f!(&a, &b, &c) + symmetric),
    )
    .unwrap();
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let result = source.simplify_algebra(&settings).unwrap();
    // The f*f branch gives C_A delta; f contracted into two symmetric legs vanishes.
    assert_eq!(
        result.expression,
        crate::color::CS.cas(Atom::num(2), reps.coad_da.to_symbolic([])) * g!(&x, &c)
    );
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
}

#[test]
fn a_single_fundamental_connection_still_joins_a_sum_of_chains() {
    let reps = crate::test_support::test_initialize();
    let [i, j, k] = std::array::from_fn(|n| reps.cof_nc.to_symbolic([Atom::num(99721 + n as i64)]));
    let [a, b, c] = std::array::from_fn(|n| {
        crate::color_t!(reps.coad_da.to_symbolic([Atom::num(99731 + n as i64)]))
    });
    let left = spenso::chain!(&i, spenso::dind!(&j); [&a, &b])
        + spenso::chain!(&i, spenso::dind!(&j); [&b, &a]);
    let right = spenso::chain!(&j, spenso::dind!(&k); [&c]);
    let source = SymbolicTensor::infer(left * right).unwrap();
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().without_trace_evaluation()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let result = source.simplify_algebra(&settings).unwrap();
    let expected = spenso::chain!(&i, spenso::dind!(&k); [&a, &b, &c])
        + spenso::chain!(&i, spenso::dind!(&k); [&b, &a, &c]);
    assert_eq!(result.expression, expected);
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
}

#[test]
fn a_rank_twelve_symmetric_color_trace_completes_without_materializing_permutations() {
    let reps = crate::test_support::test_initialize();
    let expression = crate::color::CS.symmetric_generator_trace(
        reps.cof_nc.to_symbolic([]),
        (99741..99753).map(|index| reps.coad_da.to_symbolic([Atom::num(index)])),
    );
    let source = SymbolicTensor::infer(expression.clone()).unwrap();
    assert_eq!(source.structure.slots().unwrap().len(), 12);
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    // The already symmetric output is terminal. In particular, the signed-zero
    // prepass must not open its 12! permutations merely to return it unchanged.
    IDENTITY_KERNEL_CALLS.with(|calls| calls.set([0; 3]));
    let result = source.simplify_algebra(&settings).unwrap();
    assert_eq!(result.expression, expression);
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(IDENTITY_KERNEL_CALLS.with(|calls| calls.get()[1]), 1);
    let repeated = result.simplify_algebra(&settings).unwrap();
    assert_eq!(repeated, result);
    assert_eq!(IDENTITY_KERNEL_CALLS.with(|calls| calls.get()[1]), 1);
}

#[test]
fn an_isolated_generator_trace_needs_one_color_kernel_call() {
    let reps = crate::test_support::test_initialize();
    let rank = 7;
    let word = trace!(reps.cof_nc.to_symbolic([]);
        (99_881..99_881 + rank).map(|i| crate::color_t!(reps.coad_da.to_symbolic([Atom::num(i)]))));
    let source = SymbolicTensor::infer(word).unwrap();
    // Incrementally, each of the rank-2 insertions is one round, followed by
    // a round confirming that nothing changes.
    for (color, expected) in [
        (ColorSimplifySettings::default(), 1),
        (
            ColorSimplifySettings::default().without_one_shot_traces(),
            rank as usize - 1,
        ),
    ] {
        let settings = AlgebraSettings {
            color: Some(color),
            contract: AlgebraContraction::None,
            ..Default::default()
        };
        IDENTITY_KERNEL_CALLS.with(|calls| calls.set([0; 3]));
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(IDENTITY_KERNEL_CALLS.with(|calls| calls.get()[1]), expected);
        // The certificate must not skip work: a fresh admission is unchanged.
        let fresh = SymbolicTensor::infer(result.expression.clone()).unwrap();
        assert_eq!(
            fresh.simplify_algebra(&settings).unwrap().expression,
            result.expression
        );
    }
}

#[test]
fn a_single_port_adjoint_chain_leaves_a_color_sum_factored() {
    use crate::color::simplify::ColorAlgebraSimplifier;
    use spenso::network::tags::SPENSO_TAG;
    let reps = crate::test_support::test_initialize();
    let [a, b, x, c, d, e, y] =
        std::array::from_fn(|i| reps.coad_da.to_symbolic([Atom::num(99_891 + i as i64)]));
    let matrix = |port: &Atom| {
        crate::color_f!(
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            port
        )
    };
    let ports = [&x, &c, &d, &e];
    let sum = crate::color::CS.symmetric_generator_trace(reps.cof_nc.to_symbolic([]), ports)
        + spenso::trace_sym!(reps.coad_da.to_symbolic([]); ports.map(matrix));
    // The shared contractor writes f(a,b,z) f(z,y,x) as this adjoint chain.
    // It reaches the sum through x alone, so no colour rule can span both.
    let chain = spenso::chain!(&a, &x; [matrix(&b), matrix(&y)]);
    let source = SymbolicTensor::infer(&chain * &sum).unwrap();
    let step = |settings: ColorSimplifySettings| {
        ColorAlgebraSimplifier::new(
            settings,
            SymbolicTensor::<PartialStructure>::reserved_dummies([&source]),
        )
        .step(source.expression.as_view(), true)
    };
    assert_eq!(step(ColorSimplifySettings::default()), source.expression);
    // Incremental rounds spread the chain over the summands, whose
    // states take canonical names for the contracted index x.
    use crate::IndexTooling;
    let canonical = |expression: Atom| {
        expression
            .expand()
            .canonize(spenso::structure::abstract_index::AbstractIndex::Dummy)
            .unwrap()
    };
    let stepped = step(ColorSimplifySettings::default().without_one_shot_traces());
    assert_ne!(stepped, source.expression);
    assert_eq!(canonical(stepped), canonical(source.expression.clone()));
}

#[test]
fn a_contracted_isolated_adjoint_trace_needs_one_color_kernel_call() {
    use spenso::network::tags::SPENSO_TAG;
    let reps = crate::test_support::test_initialize();
    let word = trace!(reps.coad_da.to_symbolic([]); (99_901..99_908).map(|i| {
        crate::color_f!(
            Atom::var(SPENSO_TAG.chain_in),
            Atom::var(SPENSO_TAG.chain_out),
            reps.coad_da.to_symbolic([Atom::num(i)])
        )
    }));
    let source = SymbolicTensor::infer(word).unwrap();
    for color in [
        ColorSimplifySettings::default(),
        ColorSimplifySettings::default().without_one_shot_traces(),
    ] {
        let settings = AlgebraSettings {
            color: Some(color),
            ..Default::default()
        };
        IDENTITY_KERNEL_CALLS.with(|calls| calls.set([0; 3]));
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        if color.one_shot_traces {
            // Collecting the decomposition's f pairs into adjoint chains neither
            // revokes the certificate nor makes the chains distribute the sums.
            assert_eq!(IDENTITY_KERNEL_CALLS.with(|calls| calls.get()[1]), 1);
        }
        let fresh = SymbolicTensor::infer(result.expression.clone()).unwrap();
        assert_eq!(
            fresh.simplify_algebra(&settings).unwrap().expression,
            result.expression
        );
    }
}

#[test]
fn a_cof_dimension_substitution_keeps_an_isolated_trace_certified() {
    use crate::representations::{ColorAdjoint, ColorFundamental};
    crate::test_support::test_initialize();
    let word = trace!(ColorFundamental {}.to_symbolic([Atom::num(3)]); (99_911..99_918).map(|i| {
        crate::color_t!(ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(i)]))
    }));
    let source = SymbolicTensor::infer(word).unwrap();
    for color in [
        ColorSimplifySettings::default(),
        ColorSimplifySettings::default().without_one_shot_traces(),
    ] {
        let settings = AlgebraSettings {
            color: Some(color.with_cof_dimension_invariants()),
            contract: AlgebraContraction::None,
            ..Default::default()
        };
        IDENTITY_KERNEL_CALLS.with(|calls| calls.set([0; 3]));
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert!(!result.expression.contains_symbol(crate::color::CS.cas));
        if color.one_shot_traces {
            // Writing the invariants as SU(3) numbers keeps every sum around the
            // colour factors, so no round confirms the decomposition.
            assert_eq!(IDENTITY_KERNEL_CALLS.with(|calls| calls.get()[1]), 1);
        }
        let fresh = SymbolicTensor::infer(result.expression.clone()).unwrap();
        assert_eq!(
            fresh.simplify_algebra(&settings).unwrap().expression,
            result.expression
        );
    }
}

#[test]
fn nested_vertex_sums_reach_their_loops_in_one_color_pass() {
    let reps = crate::test_support::test_initialize();
    let slots: [Atom; 12] =
        std::array::from_fn(|i| reps.coad_da.to_symbolic([Atom::num(99_921 + i as i64)]));
    let f =
        |a: usize, b: usize, c: usize| crate::color_f!(&slots[a - 1], &slots[b - 1], &slots[c - 1]);
    let v: [[Atom; 3]; 2] = std::array::from_fn(|vertex| {
        std::array::from_fn(|s| {
            Atom::var(symbolica::symbol!(format!(
                "nested_vertex_v{}s{s}",
                vertex + 1
            )))
        })
    });
    // Four-loop graph FK2338 with two four-gluon vertex sums, one inside the
    // other's distributed terms.
    let first = &v[0][0] * f(11, 1, 8) * f(11, 6, 2)
        + &v[0][1] * f(11, 1, 2) * f(11, 6, 8)
        + &v[0][2] * f(11, 1, 6) * f(11, 2, 8);
    let second = &v[1][0] * f(12, 7, 10) * f(12, 3, 9)
        + &v[1][1] * f(12, 7, 9) * f(12, 3, 10)
        + &v[1][2] * f(12, 7, 3) * f(12, 9, 10);
    let source =
        SymbolicTensor::infer(f(1, 2, 3) * f(4, 5, 6) * f(4, 5, 7) * f(8, 9, 10) * first * second)
            .unwrap();
    // FORM color.h: NA cA^4 times this polynomial in the vertex markers.
    let markers = (&v[0][0] * &v[1][0] + &v[0][2] * &v[1][0]) / Atom::num(4)
        + &v[0][1] * &v[1][0] / Atom::num(2)
        - (&v[0][0] * &v[1][1] + &v[0][2] * &v[1][1]) / Atom::num(4)
        - &v[0][1] * &v[1][1] / Atom::num(2)
        - (&v[0][0] * &v[1][2] + &v[0][2] * &v[1][2]) / Atom::num(2)
        - &v[0][1] * &v[1][2];
    let expected = Atom::var(spenso::s!(dA))
        * crate::color::CS
            .cas(Atom::num(2), reps.coad_da.to_symbolic([]))
            .pow(4)
        * markers;
    for color in [
        ColorSimplifySettings::default(),
        ColorSimplifySettings::default().without_one_shot_traces(),
    ] {
        let settings = AlgebraSettings {
            color: Some(color),
            contract: AlgebraContraction::None,
            ..Default::default()
        };
        IDENTITY_KERNEL_CALLS.with(|calls| calls.set([0; 3]));
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(result.expression.expand(), expected.expand());
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        if color.one_shot_traces {
            // Incremental rounds need six calls: each nested sum waits one round.
            assert!(IDENTITY_KERNEL_CALLS.with(|calls| calls.get()[1]) <= 4);
        }
    }
}

#[test]
fn an_internal_trace_pair_stays_independent_of_a_single_external_color_connection() {
    crate::test_support::test_initialize();
    let slots = (99801..99806)
        .map(|i| crate::representations::ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(i)]))
        .collect::<Vec<_>>();
    let word = trace!(crate::representations::ColorFundamental {}.to_symbolic([Atom::num(3)]);
        [0, 1, 0, 2].map(|i| crate::color_t!(&slots[i])));
    let source =
        SymbolicTensor::infer(word * crate::color_f!(&slots[3], &slots[1], &slots[4])).unwrap();
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let result = source.simplify_algebra(&settings).unwrap();
    assert_eq!(
        result.expression,
        -crate::color_f!(&slots[3], &slots[2], &slots[4]) / Atom::num(12)
    );
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
}

#[test]
fn factored_color_sum_contracts_only_its_color_connections() {
    crate::test_support::test_initialize();
    let [a, b, c, d] = [99821, 99822, 99823, 99824].map(|index| {
        crate::representations::ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(index)])
    });
    let mu = mink!(4, 99825);
    let nu = mink!(4, 99826);
    let left = g!(&mu, &nu) * p!(&mu);
    let right = g!(&mu, &nu) * q!(&mu);
    let source = SymbolicTensor::infer(
        crate::color_f!(&a, &b, &c) * g!(&c, &d) * &left + crate::color_f!(&a, &b, &d) * &right,
    )
    .unwrap();
    let settings = AlgebraSettings {
        color: Some(Default::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let reduced = source.simplify_algebra(&settings).unwrap();
    let expected = crate::color_f!(&a, &b, &d) * left + crate::color_f!(&a, &b, &d) * right;
    assert_eq!(reduced.expression, expected);
    assert_eq!(reduced.structure, source.structure);
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    assert_eq!(reduced.simplify_algebra(&settings).unwrap(), reduced);
}

#[test]
fn symmetric_gram_contractions_require_a_uniform_adjoint_index_space() {
    crate::test_support::test_initialize();
    let fundamental = crate::representations::ColorFundamental {}.to_symbolic([Atom::num(3)]);
    let settings = AlgebraSettings {
        color: Some(Default::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    for degree in [4, 6] {
        for uniform in [false, true] {
            let slots = (0..degree)
                .map(|index| {
                    crate::representations::ColorAdjoint {}.to_symbolic([
                        Atom::num(if uniform || index != 0 { 8 } else { 3 }),
                        Atom::num(99831 + index),
                    ])
                })
                .collect::<Vec<_>>();
            let invariant = crate::color::CS.symmetric_generator_trace(&fundamental, &slots);
            let source = SymbolicTensor::infer(invariant.pow(2)).unwrap();
            assert!(source.structure.slots().unwrap().is_empty());
            let reduced = source.simplify_algebra(&settings).unwrap();
            let gram = crate::color::CS.gram(Atom::num(degree), &fundamental, &fundamental);
            if uniform {
                assert_eq!(reduced.expression, gram);
            } else {
                assert_eq!(reduced.expression, source.expression);
                assert_ne!(reduced.expression, gram);
            }
            assert_eq!(reduced.structure, source.structure);
            assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
            assert_eq!(reduced.simplify_algebra(&settings).unwrap(), reduced);
        }
    }
}

#[test]
fn color_metric_coefficient_sums_complete_to_exact_zero() {
    crate::test_support::test_initialize();
    let [i, j, k, l] = [99851, 99852, 99853, 99854].map(|index| {
        crate::representations::ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(index)])
    });
    let a = crate::representations::ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(99855)]);
    let mu = mink!(4, 99856);
    let nu = mink!(4, 99857);
    let g0 = g!(spenso::dind!(&i), &j);
    let g1 = g!(spenso::dind!(&k), &l);
    let generators =
        crate::color_t!(&a, &i, spenso::dind!(&j)) * crate::color_t!(&a, &k, spenso::dind!(&l));
    let [x, y, z, w] = ["x", "y", "z", "w"].map(|name| {
        Atom::var(symbolica::symbol!(format!(
            "color_row_coefficients::{name}"
        )))
    });
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    for (coefficients, spectator) in [
        ([x.clone(), y.clone(), z.clone(), w.clone()], Atom::one()),
        ([p!(&mu), q!(&mu), p!(&nu), q!(&nu)], p!(&mu) * p!(&nu)),
        ([p!(&mu), q!(&mu), z.clone(), w.clone()], p!(&mu)),
        // The same cancellation must retain an open logical interface too.
        ([p!(&mu), q!(&mu), p!(&nu), q!(&nu)], Atom::one()),
        // A generated 3*z+3*w sum overlaps an original coefficient scope.
        (
            [
                Atom::num(3) * &z + Atom::num(3) * &w,
                y.clone(),
                z.clone(),
                w.clone(),
            ],
            Atom::one(),
        ),
        // SM-like coefficients: the original Lorentz sums remain opaque,
        // while the colour-generated common factor must be collected.
        (
            [
                &x * (p!(&mu) + q!(&mu)),
                &y * (p!(&mu) - q!(&mu)),
                &z * (p!(&nu) + q!(&nu)),
                &w * (p!(&nu) - q!(&nu)),
            ],
            p!(&mu) * p!(&nu),
        ),
    ] {
        let source = SymbolicTensor::infer(
            (&coefficients[0] * &g0 + &coefficients[1] * &g0)
                * (&coefficients[2] * &g1 + &coefficients[3] * &g1)
                * &generators
                * spectator,
        )
        .unwrap();
        // Sum_a tr(T^a) tr(T^a) = 0. No foreign algebra or test-side
        // expansion/factorization is needed to obtain the literal zero.
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(result.expression, Atom::Zero, "{}", source.expression);
        assert_eq!(
            result.structure.logical_slots(),
            source.structure.logical_slots()
        );
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
    }
}

#[test]
fn completed_color_coefficients_preserve_factored_spectators_and_unchanged_rows() {
    crate::test_support::test_initialize();
    let [i, j, k, l] = [99861, 99862, 99863, 99864].map(|index| {
        crate::representations::ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(index)])
    });
    let a = crate::representations::ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(99865)]);
    let mu = mink!(4, 99866);
    let nu = mink!(4, 99867);
    let generators =
        crate::color_t!(&a, &i, spenso::dind!(&j)) * crate::color_t!(&a, &k, spenso::dind!(&l));
    let g0 = g!(spenso::dind!(&i), &j);
    let g1 = g!(spenso::dind!(&k), &l);
    let zero = (p!(&mu) * &g0 + q!(&mu) * &g0)
        * (p!(&nu) * &g1 + q!(&nu) * &g1)
        * &generators
        * p!(&mu)
        * p!(&nu);
    // The crossed metrics give sum_a tr(T^a T^a) = 4 in SU(3).
    let nonzero = &generators * g!(spenso::dind!(&i), &l) * g!(spenso::dind!(&k), &j);
    let rho = mink!(4, 99868);
    let sigma = mink!(4, 99869);
    let scalar = symbolica::parse_lit!(
        (2 * (color_spectator_x + color_spectator_y) + color_spectator_z)
            ^ 12 * color_spectator_function(2 * (color_spectator_x + color_spectator_y))
    );
    let foreign = g!(&rho, &sigma) * (p!(&rho) + q!(&rho));
    let spectator = &scalar * &foreign;
    // A source coefficient with cancellable linear arithmetic must remain
    // unchanged when no colour identity owns its sum.
    let unchanged = symbolica::parse_lit!(
        1 + 2 * (color_spectator_x + color_spectator_y)
            - 2 * color_spectator_x
            - 2 * color_spectator_y
    ) * &foreign;
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    for remainder in [0, 1] {
        let source = SymbolicTensor::infer(
            (&zero + Atom::num(remainder) * &nonzero) * &spectator + &unchanged,
        )
        .unwrap();
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(
            result.expression,
            Atom::add_many([p!(&rho), q!(&rho)].map(|vector| {
                Atom::num(4 * remainder) * scalar.expand_num() * g!(&rho, &sigma) * vector
            })) + &unchanged
        );
        assert_eq!(
            result.structure.logical_slots(),
            source.structure.logical_slots()
        );
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
    }
    let [x, y, z] = [
        symbolica::parse_lit!(color_spectator_x),
        symbolica::parse_lit!(color_spectator_y),
        symbolica::parse_lit!(color_spectator_z),
    ];
    let numeric_sum = Atom::num(2) * (&x + &y);
    for untouched in [
        numeric_sum.clone(),
        Atom::num(3) * (&z + &numeric_sum),
        // Numerical coefficient collection also removes an internally zero
        // child of a changed, factored sum.
        &numeric_sum - Atom::num(2) * &x - Atom::num(2) * &y,
    ] {
        for common in [Atom::one(), foreign.clone()] {
            for (body, expected) in [
                // A sibling's numerical weight takes the collector's canonical
                // form when the completed rows are merged.
                (&nonzero + &untouched, Atom::num(4) + &untouched),
                // The same boundary holds inside another selected row, and
                // before a common indexed spectator is reattached.
                (
                    &nonzero + &z * (&nonzero + &untouched),
                    Atom::num(4) + &z * (Atom::num(4) + &untouched),
                ),
                // A vanished row no longer carries its coefficient scope.
                (
                    &zero * &untouched + &nonzero + &untouched,
                    Atom::num(4) + &untouched,
                ),
            ] {
                let source = SymbolicTensor::infer(body * &common).unwrap();
                let result = source.simplify_algebra(&settings).unwrap();
                assert_eq!(
                    result.expression,
                    expected.expand_num() * &common,
                    "{}",
                    source.expression
                );
                assert_eq!(
                    result.structure.logical_slots(),
                    source.structure.logical_slots()
                );
                assert_eq!(result.reduction_status(), ReductionStatus::Complete);
                assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
            }
        }
    }
}

#[test]
fn color_coefficient_temporaries_do_not_enter_function_normalizers() {
    crate::test_support::test_initialize();
    let exposed = Arc::new(AtomicUsize::new(0));
    let seen = Arc::clone(&exposed);
    let head = symbolica::symbol!("color_row_function_opacity"; Scalar;
    norm = move |value, output| {
        if let symbolica::atom::AtomView::Fun(function) = value
            && matches!(function.iter().next(), Some(symbolica::atom::AtomView::Var(_)))
        {
            seen.fetch_add(1, Ordering::Relaxed);
            **output = Atom::Zero;
        }
    });
    let [x, y, z] = ["x", "y", "z"]
        .map(|name| Atom::var(symbolica::symbol!(format!("color_row_function::{name}"))));
    let sum = &x + &y;
    let function = symbolica::function!(head, &sum);
    let [i, j, k, l] = [99911, 99912, 99913, 99914].map(|index| {
        crate::representations::ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(index)])
    });
    let a = crate::representations::ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(99915)]);
    let color = crate::color_t!(&a, &i, spenso::dind!(&j))
        * crate::color_t!(&a, &k, spenso::dind!(&l))
        * g!(spenso::dind!(&i), &l)
        * g!(spenso::dind!(&k), &j);
    let source = SymbolicTensor::infer(color + &z * &sum + &function).unwrap();
    let result = source
        .simplify_algebra(&AlgebraSettings {
            color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
            contract: AlgebraContraction::None,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(result.expression, Atom::num(4) + z * sum + function);
    assert_eq!(exposed.load(Ordering::Relaxed), 0);
    assert_eq!(
        result.structure.logical_slots(),
        source.structure.logical_slots()
    );
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
}

#[test]
fn completed_color_coefficients_keep_overlapping_powers_opaque() {
    crate::test_support::test_initialize();
    let color = |base: i64| {
        let [i, j, k, l] = [0, 1, 2, 3].map(|index| {
            crate::representations::ColorFundamental {}
                .to_symbolic([Atom::num(3), Atom::num(base + index)])
        });
        let a = crate::representations::ColorAdjoint {}
            .to_symbolic([Atom::num(8), Atom::num(base + 4)]);
        crate::color_t!(&a, &i, spenso::dind!(&j))
            * crate::color_t!(&a, &k, spenso::dind!(&l))
            * g!(spenso::dind!(&i), &l)
            * g!(spenso::dind!(&k), &j)
    };
    let [a, b, c, x] = ["a", "b", "c", "x"]
        .map(|name| Atom::var(symbolica::symbol!(format!("color_row_power::{name}"))));
    let power = (Atom::num(2) * (&a + &b) + &c).pow(2);
    let sum = &x + &power;
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let mu = mink!(4, 99961);
    let nu = mink!(4, 99962);
    let common = g!(&mu, &nu) * (p!(&mu) + q!(&mu));
    for phase in [Atom::i(), Atom::one() + Atom::i()] {
        let sum = &a + &b;
        for spectator in [Atom::one(), common.clone()] {
            for (body, expected) in [
                (
                    -&phase * color(99921) * &sum / Atom::num(4)
                        + &phase * color(99931) * &a / Atom::num(4)
                        + &phase * color(99941) * &b / Atom::num(4),
                    Atom::Zero,
                ),
                (color(99971) + &phase * &sum, Atom::num(4) + &phase * &sum),
                (
                    &phase * color(99921) * &sum / Atom::num(4) - Atom::num(2) * &phase * &sum,
                    -&phase * &sum,
                ),
            ] {
                let source = SymbolicTensor::infer(body * &spectator).unwrap();
                let result = source.simplify_algebra(&settings).unwrap();
                assert_eq!(
                    result.expression,
                    expected.expand_num() * &spectator,
                    "{}",
                    source.expression
                );
                assert_eq!(
                    result.structure.logical_slots(),
                    source.structure.logical_slots()
                );
                assert_eq!(result.reduction_status(), ReductionStatus::Complete);
                assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
            }
        }
    }
    for remainder in [0, 1] {
        let source = SymbolicTensor::infer(
            -color(99921) * &sum / Atom::num(4)
                + color(99931) * &power / Atom::num(4)
                + color(99941) * &x / Atom::num(4)
                + Atom::num(remainder) * &power,
        )
        .unwrap();
        let result = source.simplify_algebra(&settings).unwrap();
        // Each independently indexed color factor is sum_a tr(T^a T^a)=4.
        // Cancellation must compare the original power before opening its base.
        assert_eq!(
            result.expression,
            (Atom::num(remainder) * &power).expand_num()
        );
        assert_eq!(
            result.structure.logical_slots(),
            source.structure.logical_slots()
        );
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
    }
    let [i, j, k, l] = [99951, 99952, 99953, 99954].map(|index| {
        crate::representations::ColorFundamental {}.to_symbolic([Atom::num(3), Atom::num(index)])
    });
    let adjoint =
        crate::representations::ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(99955)]);
    // The uncrossed contraction is sum_a tr(T^a) tr(T^a) = 0.
    let zero = crate::color_t!(&adjoint, &i, spenso::dind!(&j))
        * crate::color_t!(&adjoint, &k, spenso::dind!(&l))
        * g!(spenso::dind!(&i), &j)
        * g!(spenso::dind!(&k), &l);
    let sum = &a + &b;
    let residual = Atom::num(4) * &sum - Atom::num(4) * &a - Atom::num(4) * &b;
    for spectator in [Atom::one(), common.clone()] {
        let cancelled =
            color(99921) * &sum - color(99931) * &a - color(99941) * &b + &zero * &residual;
        for (body, expected) in [
            (cancelled.clone(), Atom::Zero),
            (color(99971) + &x * &cancelled, Atom::num(4)),
            (&x * (color(99971) + &cancelled), Atom::num(4) * &x),
            (color(99971) + &c + &x * &cancelled, Atom::num(4) + &c),
            (
                &x * (color(99971) + &c + &a * &cancelled),
                &x * (Atom::num(4) + &c),
            ),
        ] {
            let source = SymbolicTensor::infer(body * &spectator).unwrap();
            // The newly generated root equals a coefficient of the vanished
            // row, also inside a factored row or beside another changed row.
            let result = source.simplify_algebra(&settings).unwrap();
            assert_eq!(
                result.expression,
                expected * &spectator,
                "{}",
                source.expression
            );
            assert_eq!(
                result.structure.logical_slots(),
                source.structure.logical_slots()
            );
            assert_eq!(result.reduction_status(), ReductionStatus::Complete);
            assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
        }
    }
    for (coefficient, expanded) in [
        (
            Atom::num(2) * (&a + &b),
            Atom::num(2) * &a + Atom::num(2) * &b,
        ),
        (
            Atom::num(2) * (&c + Atom::num(3) * (&a + &b)),
            Atom::num(2) * &c + Atom::num(6) * &a + Atom::num(6) * &b,
        ),
    ] {
        let source = SymbolicTensor::infer(
            -color(99921) * (&x + &coefficient) / Atom::num(4)
                + color(99931) * &expanded / Atom::num(4)
                + color(99941) * &x / Atom::num(4),
        )
        .unwrap();
        // Direct numerical multiples of retained sums must also reopen,
        // including through several nested coefficient scopes.
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(result.expression, Atom::Zero);
        assert_eq!(
            result.structure.logical_slots(),
            source.structure.logical_slots()
        );
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
        for remainder in [&coefficient, &expanded] {
            let source = SymbolicTensor::infer(
                -color(99921) * (&x + &coefficient) / Atom::num(4)
                    + color(99931) * &x / Atom::num(4)
                    + remainder,
            )
            .unwrap();
            // An untouched occurrence must still cancel against a changed
            // row that carries the same scope, also in its expanded form.
            let result = source.simplify_algebra(&settings).unwrap();
            assert_eq!(result.expression, Atom::Zero, "{}", source.expression);
            assert_eq!(
                result.structure.logical_slots(),
                source.structure.logical_slots()
            );
            assert_eq!(result.reduction_status(), ReductionStatus::Complete);
            assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
        }
        for spectator in [Atom::one(), common.clone()] {
            for (body, expected) in [
                // The first two changed contributions vanish collectively.
                // The surviving coefficient still has its exact numerical value.
                (
                    -color(99921) * &coefficient / Atom::num(4)
                        + color(99931) * &coefficient / Atom::num(4)
                        + color(99941)
                        + &coefficient,
                    Atom::num(4) + &coefficient,
                ),
                // The same boundary is required within a factored selected
                // row, not only when the collector groups top-level rows.
                (
                    color(99941)
                        + &x * (-color(99921) * &coefficient / Atom::num(4)
                            + color(99931) * &coefficient / Atom::num(4)
                            + color(99941)
                            + &coefficient),
                    Atom::num(4) + &x * (Atom::num(4) + &coefficient),
                ),
                // Opening another coefficient can expose the cancellation
                // later in completion, with an unchanged remainder surviving.
                (
                    -color(99921) * (&x + &coefficient) / Atom::num(4)
                        + color(99931) * &coefficient / Atom::num(4)
                        + color(99941) * &x / Atom::num(4)
                        + color(99971)
                        + &coefficient,
                    Atom::num(4) + &coefficient,
                ),
                // Equal coefficients can cancel only after numerical
                // distribution; the surviving occurrence is normalized numerically.
                (
                    -color(99921) * &expanded / Atom::num(4)
                        + color(99931) * &coefficient / Atom::num(4)
                        + color(99941)
                        + &coefficient,
                    Atom::num(4) + &coefficient,
                ),
                // The same spelling has distinct owners in these contexts:
                // the outer occurrence must cancel, and the inner one survives.
                (
                    color(99921) * &coefficient / Atom::num(4)
                        - color(99931) * &expanded / Atom::num(4)
                        + &x * (color(99941) + &coefficient),
                    &x * (Atom::num(4) + &coefficient),
                ),
                // The cancellation must precede opening the retained child,
                // including when an unrelated sibling prevents global factoring.
                (
                    color(99971)
                        + &x * (-color(99921) * &expanded / Atom::num(4)
                            + color(99931) * &coefficient / Atom::num(4)
                            + color(99941)
                            + &coefficient),
                    Atom::num(4) + &x * (Atom::num(4) + &coefficient),
                ),
                // The same contexts can multiply into a power after a child
                // changes; that power must still share its cancellation key.
                (
                    color(99971)
                        + &x * (color(99921) * &coefficient / Atom::num(4)
                            - color(99931) * &expanded / Atom::num(4)
                            + &x * (color(99941) + &coefficient)),
                    Atom::num(4) + x.pow(2) * (Atom::num(4) + &coefficient),
                ),
                // Consuming only one of two copies must preserve the exact
                // surviving numerical multiple.
                (
                    color(99971) - color(99921) * &expanded / Atom::num(4)
                        + Atom::num(2) * &coefficient,
                    Atom::num(4) + &coefficient,
                ),
            ] {
                let source = SymbolicTensor::infer(body * &spectator).unwrap();
                let result = source.simplify_algebra(&settings).unwrap();
                assert_eq!(
                    result.expression,
                    expected.expand_num() * &spectator,
                    "{}",
                    source.expression
                );
                assert_eq!(
                    result.structure.logical_slots(),
                    source.structure.logical_slots()
                );
                assert_eq!(result.reduction_status(), ReductionStatus::Complete);
                assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
            }
        }
    }
}

#[test]
fn canonical_color_rows_cancel_against_their_original_expression() {
    let reps = crate::test_support::test_initialize();
    let word = trace!(reps.cof_nc.to_symbolic([]);
        (99871..99875).map(|index| crate::color_t!(reps.coad_da.to_symbolic([Atom::num(index)]))));
    let source = SymbolicTensor::infer(word.clone()).unwrap();
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let canonical = source.simplify_algebra(&settings).unwrap();
    assert_ne!(canonical.expression, word);
    assert!(!canonical.expression.is_zero());
    let difference = SymbolicTensor::infer(Atom::num(2) * (&word - &canonical.expression)).unwrap();
    let result = difference.simplify_algebra(&settings).unwrap();
    assert_eq!(result.expression, Atom::Zero);
    assert_eq!(
        result.structure.logical_slots(),
        difference.structure.logical_slots()
    );
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
}

/// Weighted rows publish equal canonical states once: the two disabled
/// components of 2*(sym(x,y,z) g(x,y) + sym(u,v,z) g(u,v)) relabel to one
/// state that merges, beside the reduced active summand. The result stays
/// exact and keeps its scalar spectator as a factor.
#[test]
fn weighted_disabled_trace_rows_merge_and_stay_exact() {
    let reps = test_support::test_initialize();
    let [x, y, u, v, z, a, b, c] =
        array::from_fn(|index| reps.coad_da.to_symbolic([Atom::num(99881 + index as i64)]));
    let first =
        CS.symmetric_generator_trace(reps.cof_nc.to_symbolic([]), [&x, &y, &z]) * g!(&x, &y);
    let second =
        CS.symmetric_generator_trace(reps.cof_nc.to_symbolic([]), [&u, &v, &z]) * g!(&u, &v);
    let retained = Atom::num(2) * (first + second);
    let active = color_f!(&z, &a, &b)
        * color_f!(&c, &a, &b)
        * trace!(reps.cof_nc.to_symbolic([]); [color_t!(&c)]);
    let reduced = CS.cas(Atom::num(2), reps.coad_da.to_symbolic([]))
        * trace!(reps.cof_nc.to_symbolic([]); [color_t!(&z)]);
    let spectator = parse_lit!((disabled_row_x + disabled_row_y) ^ 12);
    let source = SymbolicTensor::infer((active + &retained) * &spectator).unwrap();
    let disabled = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().without_trace_evaluation()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let result = source.simplify_algebra(&disabled).unwrap();
    let AtomView::Mul(product) = result.expression.as_view() else {
        panic!("{}", result.expression)
    };
    let factors = product.iter().collect::<Vec<_>>();
    assert!(
        factors.contains(&spectator.as_view()),
        "{}",
        result.expression
    );
    // The reduced summand and a single merged disabled state remain.
    assert!(
        factors
            .iter()
            .any(|factor| matches!(factor, AtomView::Add(sum)
            if sum.get_nargs() == 2 && sum.iter().any(|term| term == reduced.as_view()))),
        "{}",
        result.expression
    );
    assert_eq!(
        result.structure.logical_slots(),
        source.structure.logical_slots()
    );
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&disabled).unwrap(), result);
    let enabled = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        ..disabled
    };
    assert_eq!(
        result.simplify_algebra(&enabled).unwrap().expression,
        source.simplify_algebra(&enabled).unwrap().expression
    );
}

#[test]
fn short_generator_traces_require_compatible_adjoint_spaces() {
    crate::test_support::test_initialize();
    let representation = crate::representations::ColorFundamental {}.to_symbolic([Atom::num(3)]);
    let settings = AlgebraSettings {
        color: Some(Default::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    for dimensions in [
        vec![3, 8],
        vec![3, 8, 8],
        vec![3, 8, 8, 8],
        vec![3, 3, 8, 8],
    ] {
        let slots = dimensions
            .iter()
            .enumerate()
            .map(|(index, dimension)| {
                let index = if dimensions == [3, 3, 8, 8] && index == 1 {
                    0
                } else {
                    index
                };
                crate::representations::ColorAdjoint {}
                    .to_symbolic([Atom::num(*dimension), Atom::num(99841 + index)])
            })
            .collect::<Vec<_>>();
        let source = SymbolicTensor::infer(
            trace!(&representation; slots.iter().map(|slot| crate::color_t!(slot))),
        )
        .unwrap();
        let reduced = source.simplify_algebra(&settings).unwrap();
        assert_eq!(reduced.expression, source.expression);
        assert_eq!(reduced.structure, source.structure);
        assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
        assert_eq!(reduced.simplify_algebra(&settings).unwrap(), reduced);
    }
}

#[test]
fn completed_color_coefficients_do_not_oscillate_with_contraction() {
    let reps = test_support::test_initialize();
    let [a, b] = [99991, 99992].map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
    let [i, j] = [99993, 99994].map(|index| reps.cof_nc.to_symbolic([Atom::num(index)]));
    let source = SymbolicTensor::infer(chain!(&i, dind!(&j);
        [color_t!(&a), color_t!(&b), color_t!(&a)]))
    .unwrap();
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default()),
        max_passes: Some(8),
        ..Default::default()
    };
    let result = source.simplify_algebra(&settings).unwrap();
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.structure, source.structure);
    let coefficient = CS.cas(Atom::num(2), reps.cof_nc.to_symbolic([]))
        - CS.cas(Atom::num(2), reps.coad_da.to_symbolic([])) / Atom::num(2);
    let expected = coefficient * chain!(&i, dind!(&j); [color_t!(&b)]);
    assert_eq!(result.expression.expand(), expected.expand());
    let fresh = SymbolicTensor::infer(result.expression.clone()).unwrap();
    assert_eq!(fresh.simplify_algebra(&settings).unwrap(), result);
}

#[test]
fn color_metric_sum_factors_collect_their_indexed_coefficients() {
    test_support::test_initialize();
    let [a, b, c, d, r, s] = array::from_fn(|index| {
        ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::num(100101 + index as i64)])
    });
    let [mu, nu, rho] = [100111, 100112, 100113]
        .map(|index| Minkowski {}.to_symbolic([Atom::num(4), Atom::num(index)]));
    let x = p!(&mu) * g!(&nu, &rho);
    let y = p!(&nu) * g!(&mu, &rho);
    let z = p!(&rho) * g!(&mu, &nu);
    let metric = g!(&c, &d);
    let bubble = color_f!(&a, &b, &r) * color_f!(&a, &b, &s);
    // A metric-only color sum is coefficient arithmetic, even when nested
    // inside the complete row of another identity. This is the quartic-vertex
    // cancellation in four-loop diagrams D3971/D3991/D4013/D4033.
    let zero = -&x * &metric + Atom::num(3) * &y * &metric
        - Atom::num(2) * &z * &metric
        - Atom::num(2) * (&y - &z) * &metric
        - (-&x + &y) * &metric;
    assert!(!zero.is_zero());
    let spectator = parse_lit!((metric_sum_x + metric_sum_y) ^ 12);
    let settings = AlgebraSettings {
        color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    for remainder in [0, 1] {
        let source = SymbolicTensor::infer(
            &bubble * (&zero + Atom::num(remainder) * &x * &metric) * &spectator,
        )
        .unwrap();
        let result = source.simplify_algebra(&settings).unwrap();
        assert_eq!(
            result.expression,
            Atom::num(3 * remainder) * &x * &metric * g!(&r, &s) * &spectator
        );
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
    }
}
