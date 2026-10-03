use super::*;
use crate::gamma;
use spenso::{g, mink, p, q, structure::partial::PartialStructureExt, trace};
use symbolica::atom::{Atom, AtomCore};

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
    let s = spenso::vector_symbol!("completion_cache_inexact::s");
    let rounded = Atom::num(symbolica::domains::float::Float::parse("0.1", Some(11)).unwrap());
    let expression = (rounded * p!(&port) + q!(&port))
        * (symbolica::function!(r, &port) + symbolica::function!(s, &port));
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
            SymbolicTensor::reserved_dummies([&source]),
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
        ColorAlgebraSimplifier::new(settings, SymbolicTensor::reserved_dummies([&source]))
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
