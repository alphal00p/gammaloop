use super::*;
use crate::{gamma, gamma5};
use spenso::structure::partial::PartialStructureExt;
use spenso::{g, mink, p, q, trace};
use std::collections::HashSet;
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    symbol,
};

#[test]
fn spinor_vectors_contract_into_the_matrix_ports_of_a_collected_gamma_chain() {
    let reps = crate::test_support::test_initialize();
    let spin = reps.bis4.to_symbolic([]);
    let left = reps.bis4.to_symbolic([Atom::num(99841)]);
    let right = reps.bis4.to_symbolic([Atom::num(99843)]);
    let incoming = spenso::vector_symbol!("chain_endpoint_incoming");
    let outgoing = spenso::vector_symbol!("chain_endpoint_outgoing");
    let spectator = symbolica::parse_lit!((chain_endpoint_x + chain_endpoint_y) ^ 12);
    let word = spenso::chain!(&left, &right, gamma!(p!(mink!(4))));
    let source = SymbolicTensor::infer(
        &spectator
            * symbolica::function!(incoming, &left)
            * word
            * symbolica::function!(outgoing, &right),
    )
    .unwrap();
    let expected = &spectator
        * symbolica::function!(
            crate::dirac::AGS.gamma,
            symbolica::function!(incoming, &spin),
            symbolica::function!(outgoing, &spin),
            p!(mink!(4))
        );
    for result in [
        source.contract(ContractSettings::default()).unwrap(),
        source
            .contract(ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap(),
        source
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::default()),
                ..Default::default()
            })
            .unwrap(),
    ] {
        assert_eq!(result.expression, expected);
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        let admitted = SymbolicTensor::infer(result.expression.clone()).unwrap();
        assert_eq!(admitted.structure, result.structure);
    }
    assert_eq!(
        source
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::default()),
                contract: AlgebraContraction::None,
                ..Default::default()
            })
            .unwrap()
            .expression,
        source.expression
    );
}

#[test]
fn gamma_identity_precedes_attached_spinor_endpoint_substitutions() {
    let reps = crate::test_support::test_initialize();
    let left = reps.bis4.to_symbolic([Atom::num(99861)]);
    let right = reps.bis4.to_symbolic([Atom::num(99863)]);
    let incoming = spenso::vector_symbol!("gamma_pair_incoming");
    let outgoing = spenso::vector_symbol!("gamma_pair_outgoing");
    let word = spenso::chain!(
        &left,
        &right,
        gamma!(mink!(4, 99865)),
        gamma!(mink!(4, 99865))
    );
    let source = SymbolicTensor::infer(
        symbolica::function!(incoming, &left) * word * symbolica::function!(outgoing, &right),
    )
    .unwrap();
    for contract in [AlgebraContraction::Fully, AlgebraContraction::None] {
        let result = source
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::default()),
                contract,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(
            result.expression,
            Atom::num(4)
                * g!(
                    symbolica::function!(incoming, reps.bis4.to_symbolic([])),
                    symbolica::function!(outgoing, reps.bis4.to_symbolic([]))
                )
        );
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    }
}

#[test]
fn gamma_prerequisites_still_substitute_lorentz_vectors() {
    let reps = crate::test_support::test_initialize();
    let left = reps.bis4.to_symbolic([Atom::num(99871)]);
    let right = reps.bis4.to_symbolic([Atom::num(99873)]);
    let mu = mink!(4, 99875);
    let nu = mink!(4, 99877);
    let source = SymbolicTensor::infer(
        p!(&mu) * p!(&nu) * spenso::chain!(&left, &right, gamma!(&mu), gamma!(&nu)),
    )
    .unwrap();
    let result = source
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            contract: AlgebraContraction::None,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        result.expression,
        g!(&left, &right) * g!(p!(mink!(4)), p!(mink!(4)))
    );
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
}

#[test]
fn contracted_chain_copies_keep_order_and_reserve_unrelated_dummy_names() {
    use spenso::structure::{
        abstract_index::AbstractIndex,
        slot::{DummyAind, ParseableAind},
    };
    use symbolica::atom::FunctionBuilder;

    let reps = crate::test_support::test_initialize();
    let AbstractIndex::Dummy(seed) = AbstractIndex::new_dummy() else {
        unreachable!()
    };
    let reserved = reps
        .bis4
        .to_symbolic([AbstractIndex::new_dummy_at(seed + 1).to_atom()]);
    let left = reps.bis4.to_symbolic([Atom::num(99851)]);
    let right = reps.bis4.to_symbolic([Atom::num(99853)]);
    let incoming = spenso::vector_symbol!("chain_scope_incoming");
    let outgoing = spenso::vector_symbol!("chain_scope_outgoing");
    let spectator = FunctionBuilder::new(spenso::tensor_symbol!("chain_scope_spectator"))
        .add_arg(&reserved)
        .finish();
    let word = spenso::chain!(
        &left,
        &right,
        spenso::chain_factor!(chain_scope_first, in, out),
        spenso::chain_factor!(chain_scope_second, in, out)
    );
    let source = SymbolicTensor::infer(
        &spectator
            * symbolica::function!(incoming, &left)
            * word
            * symbolica::function!(outgoing, &right),
    )
    .unwrap();
    let settings = ContractSettings {
        collect_chains: false,
        collect_traces: false,
        ..Default::default()
    };
    let first = source.contract(settings).unwrap();
    let second = source.contract(settings).unwrap();
    let mut internals = Vec::new();
    for result in [&first, &second] {
        assert_eq!(result.structure, source.structure);
        SymbolicTensor::validate_atom(&result.expression).unwrap();
        let mut ports = std::collections::HashMap::new();
        result.expression.visitor(&mut |value| {
            if let AtomView::Fun(function) = value {
                let symbol = function.get_symbol();
                if symbol == spenso::tensor_symbol!(chain_scope_first)
                    || symbol == spenso::tensor_symbol!(chain_scope_second)
                {
                    ports.insert(
                        symbol,
                        function
                            .iter()
                            .map(|argument| argument.to_owned())
                            .collect::<Vec<_>>(),
                    );
                }
            }
            true
        });
        let a = &ports[&spenso::tensor_symbol!(chain_scope_first)];
        let b = &ports[&spenso::tensor_symbol!(chain_scope_second)];
        assert_eq!(
            a[0],
            symbolica::function!(incoming, reps.bis4.to_symbolic([]))
        );
        assert_eq!(
            b[1],
            symbolica::function!(outgoing, reps.bis4.to_symbolic([]))
        );
        assert_eq!(a[1], b[0]);
        assert_ne!(a[1], reserved);
        internals.push(a[1].clone());
    }
    assert_ne!(internals[0], internals[1]);
    SymbolicTensor::validate_atom(&(&first.expression * &second.expression)).unwrap();
}

#[test]
fn a_zero_colour_identity_discharges_a_removed_deferred_gamma_scope() {
    let reps = crate::test_support::test_initialize();
    let hook = symbol!("deferred_gamma_scalar_hook"; Scalar; norm = |_, _| {});
    let gamma_trace = trace!(
        reps.bis4.to_symbolic([]),
        gamma!(p!(mink!(4))),
        gamma!(q!(mink!(4)))
    );
    let coefficient = (symbolica::function!(hook, gamma_trace)
        + Atom::var(symbol!("deferred_gamma_scalar_weight")))
    .pow(2);
    let colour_trace = trace!(
        reps.cof_nc.to_symbolic([]),
        crate::color_t!(reps.coad_da.to_symbolic([Atom::num(99831)]))
    );
    let source = SymbolicTensor::infer(coefficient * colour_trace).unwrap();
    let gamma_only = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let deferred = source.simplify_algebra(&gamma_only).unwrap();
    assert_eq!(deferred.reduction_status(), ReductionStatus::Deferred);
    assert_eq!(deferred.expression, source.expression);

    // Tr(T^a)=0 removes the entire callback-sensitive coefficient exactly.
    // No deferred gamma work survives in the resulting typed zero.
    let result = source
        .simplify_algebra(&AlgebraSettings {
            color: Some(ColorSimplifySettings::default()),
            ..gamma_only
        })
        .unwrap();
    assert!(result.expression.is_zero());
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
}

#[test]
fn free_trace_results_discharge_cleanup_without_entering_contraction() {
    use crate::tensor::contract::CONTRACT_DOMAIN_CALLS;
    let reps = crate::test_support::test_initialize();
    let slots = (0..8)
        .map(|index| mink!(4, { Atom::num(98600 + index) }))
        .collect::<Vec<_>>();
    let source = SymbolicTensor::infer(
        trace!(reps.bis4.to_symbolic([]); slots.iter().map(|slot| gamma!(slot))),
    )
    .unwrap();
    CONTRACT_DOMAIN_CALLS.with(|count| count.set(0));
    let reduced = source
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            // The trace identity uses this budget. Its free-port result
            // requires no additional contraction.
            max_passes: Some(1),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    assert_eq!(reduced.structure, source.structure);
    assert_eq!(
        CONTRACT_DOMAIN_CALLS.with(|count| count.get()),
        0,
        "free-port results must use their observations before contraction setup"
    );
}

#[test]
fn contracted_scalar_results_reuse_producer_observations() {
    use observation::{INITIAL_SCANS, REPLACEMENT_SCANS};
    use spenso::network::tags::SPENSO_TAG;

    let reps = crate::test_support::test_initialize();
    let vector = spenso::vector_symbol!("contracted_scalar_observation_vector");
    let spectator = symbolica::parse_lit!((contracted_scalar_x + contracted_scalar_y) ^ 20);
    for representation in [&reps.mink4, &reps.mink_d] {
        let compact = representation.to_symbolic([]);
        let indexed = representation.to_symbolic([Atom::num(99731)]);
        let component = |index, port: &Atom| symbolica::function!(vector, Atom::num(index), port);
        let coefficient = &spectator
            * symbolica::function!(
                SPENSO_TAG.dot,
                component(4, &compact),
                component(5, &compact)
            );
        let source = SymbolicTensor::infer(
            &coefficient
                * (component(0, &indexed) + component(1, &indexed))
                * (component(2, &indexed) + component(3, &indexed)),
        )
        .unwrap();
        source.reduction_observations();
        let scans = INITIAL_SCANS.with(|count| count.get());
        REPLACEMENT_SCANS.with(|count| count.set(0));

        let contracted = source.contract(Default::default()).unwrap();
        assert_eq!(contracted.reduction_status(), ReductionStatus::Complete);
        assert!(contracted.reduction_observations().is_closed_scalar());
        let dotted = contracted.to_dots().unwrap();
        assert!(dotted.reduction_observations().is_closed_scalar());
        assert_eq!(REPLACEMENT_SCANS.with(|count| count.get()), 0);
        assert_eq!(INITIAL_SCANS.with(|count| count.get()), scans);
        assert!(Arc::ptr_eq(
            dotted.reduction_observations(),
            dotted
                .contract(Default::default())
                .unwrap()
                .reduction_observations()
        ));

        // Bilinearity supplies four pairings; the unrelated scalar power is
        // cancelled before comparing the small tensor-dependent polynomial.
        let mut pairings = Vec::with_capacity(4);
        for left in 0..2 {
            for right in 2..4 {
                pairings.push(symbolica::function!(
                    SPENSO_TAG.dot,
                    component(left, &compact),
                    component(right, &compact)
                ));
            }
        }
        let expected = Atom::add_many(pairings);
        assert_eq!((&dotted.expression / &coefficient).expand(), expected);
    }
}

#[test]
fn contraction_certificates_preserve_unselected_trace_candidates() {
    let reps = crate::test_support::test_initialize();
    let word = trace!(reps.bis4.to_symbolic([]);
        [gamma!(p!(mink!(4))), gamma!(q!(mink!(4)))]);
    let mu = mink!(4, 99732);
    let source = SymbolicTensor::infer(&word * p!(&mu) * q!(&mu)).unwrap();
    let contracted = source.contract(Default::default()).unwrap();
    assert_eq!(
        contracted.expression,
        &word * g!(p!(mink!(4)), q!(mink!(4)))
    );
    assert!(!contracted.reduction_observations().is_closed_scalar());
    let reduced = contracted
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    assert_eq!(
        reduced.expression,
        Atom::num(4) * g!(p!(mink!(4)), q!(mink!(4))).pow(2)
    );
}

#[test]
fn closed_trace_results_reuse_observations_through_dot_conversion() {
    use observation::{INITIAL_DOMAINS, REPLACEMENT_SCANS};
    use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};

    let reps = crate::test_support::test_initialize();
    let vector = spenso::vector_symbol!("closed_trace_observation_vector");
    let spectator = symbolica::parse_lit!((closed_trace_x + closed_trace_y) ^ 8);
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        ..Default::default()
    };
    for representation in [&reps.mink4, &reps.mink_d] {
        let vectors = (0..4)
            .map(|index| {
                FunctionBuilder::new(vector)
                    .add_arg(Atom::num(index))
                    .add_arg(representation.to_symbolic([]))
                    .finish()
            })
            .collect::<Vec<_>>();
        for length in [2, 4] {
            let source = SymbolicTensor::infer(
                &spectator
                    * trace!(reps.bis4.to_symbolic([]);
                        vectors[..length].iter().map(|vector| gamma!(vector))),
            )
            .unwrap();
            source.reduction_observations();
            INITIAL_DOMAINS.with(|domains| domains.borrow_mut().clear());
            REPLACEMENT_SCANS.with(|count| count.set(0));

            let reduced = source.simplify_algebra(&settings).unwrap();
            assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
            assert!(Arc::ptr_eq(
                reduced.reduction_observations(),
                reduced
                    .simplify_algebra(&settings)
                    .unwrap()
                    .reduction_observations()
            ));
            let dotted = reduced.to_dots().unwrap();
            assert!(Arc::ptr_eq(
                dotted.reduction_observations(),
                dotted.to_dots().unwrap().reduction_observations()
            ));
            assert_eq!(REPLACEMENT_SCANS.with(|count| count.get()), 0);
            INITIAL_DOMAINS.with(|domains| {
                assert!(
                    domains.borrow().iter().all(|domain| {
                        !domain.contains_symbol(ETS.metric)
                            && !domain.contains_symbol(SPENSO_TAG.dot)
                    }),
                    "generated pairing results must retain their producer observations: {:?}",
                    *domains.borrow()
                );
            });

            // These two Clifford identities are independent of the production
            // recurrence. Only the test oracle distributes the small pairing sum.
            for (result, head) in [(reduced, ETS.metric), (dotted, SPENSO_TAG.dot)] {
                let pair = |left: usize, right: usize| {
                    FunctionBuilder::new(head)
                        .add_arg(&vectors[left])
                        .add_arg(&vectors[right])
                        .finish()
                };
                let pairing = if length == 2 {
                    pair(0, 1)
                } else {
                    pair(0, 1) * pair(2, 3) - pair(0, 2) * pair(1, 3) + pair(0, 3) * pair(1, 2)
                };
                let expected = Atom::num(4) * &spectator * pairing;
                let resolved = result;
                assert_eq!(resolved.structure, source.structure);
                assert!((&resolved.expression - expected).expand().is_zero());
                assert!(
                    matches!(resolved.expression.as_view(), AtomView::Mul(product)
                    if product.iter().any(|factor| factor == spectator.as_view()))
                );
            }
        }
    }
}

#[test]
fn closed_trace_observations_preserve_a_surviving_epsilon_identity() {
    let reps = crate::test_support::test_initialize();
    let compact = reps.mink4.to_symbolic([]);
    let epsilon = crate::epsilon!(mink!(2, 99501), mink!(2, 99502));
    let spectator = symbolica::parse_lit!((closed_epsilon_x + closed_epsilon_y) ^ 8);
    let source = SymbolicTensor::infer(
        &spectator
            * epsilon.pow(2)
            * trace!(
                reps.bis4.to_symbolic([]),
                gamma!(p!(&compact)),
                gamma!(q!(&compact))
            ),
    )
    .unwrap();
    let gamma_only = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        ..Default::default()
    };
    let reduced = source.simplify_algebra(&gamma_only).unwrap();
    assert!(
        reduced
            .expression
            .contains_symbol(*crate::epsilon::EPSILON_SYMBOL)
    );
    let settings = AlgebraSettings {
        epsilon: true,
        ..gamma_only
    };
    let completed = reduced.simplify_algebra(&settings).unwrap();
    assert_eq!(completed.reduction_status(), ReductionStatus::Complete);
    assert_eq!(
        completed.expression,
        Atom::num(8) * spectator * g!(p!(&compact), q!(&compact))
    );
    assert!(Arc::ptr_eq(
        completed.reduction_observations(),
        completed
            .simplify_algebra(&settings)
            .unwrap()
            .reduction_observations()
    ));
}

#[test]
fn closed_trace_observations_keep_generated_leaves_available_to_collection() {
    use spenso::{
        network::{library::symbolic::ETS, tags::SPENSO_TAG},
        shadowing::TensorCollectFilter,
    };

    let reps = crate::test_support::test_initialize();
    let vector = spenso::vector_symbol!("closed_trace_collection_vector");
    let vectors = (0..4)
        .map(|index| {
            FunctionBuilder::new(vector)
                .add_arg(Atom::num(index))
                .add_arg(reps.mink4.to_symbolic([]))
                .finish()
        })
        .collect::<Vec<_>>();
    for head in [ETS.metric, SPENSO_TAG.dot] {
        let source = SymbolicTensor::infer(
            trace!(reps.bis4.to_symbolic([]); vectors.iter().map(|vector| gamma!(vector))),
        )
        .unwrap();
        let reduced = source
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::default()),
                ..Default::default()
            })
            .unwrap();
        let reduced = if head == SPENSO_TAG.dot {
            reduced.to_dots().unwrap()
        } else {
            reduced
        };
        let target = FunctionBuilder::new(head)
            .add_arg(&vectors[0])
            .add_arg(&vectors[1])
            .finish();
        let body = reduced;
        let mut selected_target = false;
        let mut mapped_target = false;
        body.collect_with_map(
            crate::tensor::CollectionMode::Monomials,
            None,
            |value| {
                let selected = value == target.as_view();
                selected_target |= selected;
                selected
            },
            |value, _, _| {
                mapped_target |= value.expression == target;
                Ok(value)
            },
        )
        .unwrap();
        assert!(selected_target, "the selector must see generated leaves");
        assert!(
            mapped_target,
            "the collector must process the selected leaf"
        );

        if head == ETS.metric {
            let rows = body
                .coefficient_list(TensorCollectFilter::<0>::Metrics)
                .unwrap();
            assert_eq!(rows.len(), 3);
            assert!(rows.iter().all(|(selected, coefficient)| {
                selected.expression.contains_symbol(ETS.metric)
                    && !coefficient.expression.contains_symbol(ETS.metric)
            }));
        }
    }
}

#[test]
fn unchanged_scalar_scope_is_observed_once_and_completed_reruns_reuse_it() {
    use observation::{INITIAL_DOMAINS, INITIAL_SCANS};
    crate::test_support::test_initialize();
    let value =
        SymbolicTensor::infer((Atom::one() + p!(mink!(4, 98701)) * q!(mink!(4, 98701))).pow(3))
            .unwrap();
    INITIAL_SCANS.with(|count| count.set(0));
    INITIAL_DOMAINS.with(|domains| domains.borrow_mut().clear());
    value.reduction_observations();
    assert_eq!(INITIAL_SCANS.with(|count| count.get()), 1);
    let reduced = value.contract(Default::default()).unwrap();
    let domains = 1;
    assert_eq!(
        INITIAL_SCANS.with(|count| count.get()),
        domains,
        "the original domain is observed once: {:?}",
        INITIAL_DOMAINS.with(|domains| domains.borrow().clone())
    );
    assert_eq!(
        reduced.expression,
        (Atom::one() + g!(p!(mink!(4)), q!(mink!(4)))).pow(3)
    );
    let repeated = reduced.contract(Default::default()).unwrap();
    assert!(Arc::ptr_eq(
        reduced.reduction_observations(),
        repeated.reduction_observations()
    ));
    assert_eq!(INITIAL_SCANS.with(|count| count.get()), domains);
}

#[test]
fn regional_inventory_preserves_a_surviving_head_and_reuses_scalar_spectators() {
    use observation::{DomainObservations, REUSED_REGIONS};
    crate::test_support::test_initialize();
    let epsilon = *crate::epsilon::EPSILON_SYMBOL;
    let left = FunctionBuilder::new(epsilon)
        .add_arg(Atom::var(symbol!("inventory_left")))
        .finish();
    let right = FunctionBuilder::new(epsilon)
        .add_arg(Atom::var(symbol!("inventory_right")))
        .finish();
    let spectator =
        (Atom::var(symbol!("inventory_x")) + Atom::var(symbol!("inventory_y"))).pow(100);
    let original = &left * &right * &spectator;
    let tensor =
        SymbolicTensor::from_normalized_parts(original, PartialStructure::from_logical_slots([]));
    let observed = tensor.reduction_observations();
    assert_eq!(observed.candidates.counts[4], 2);
    assert_eq!(observed.epsilon_degree, 2);
    REUSED_REGIONS.with(|count| count.set(0));
    let next: DomainObservations = observed.updated((&right * &spectator).as_view());
    assert_eq!(next.candidates.counts[4], 1);
    assert_eq!(next.epsilon_degree, 1);
    assert!(next.candidates.symbols[4]);
    assert_eq!(REUSED_REGIONS.with(|count| count.get()), 2);
    let sum = observed.updated((&left + &right).as_view());
    assert_eq!(sum.epsilon_degree, 1);
    let squared = sum.updated((&left + &right).pow(2).as_view());
    assert_eq!(squared.epsilon_degree, 2);
}

#[test]
fn gamma_prerequisites_do_not_enable_epsilon_identities() {
    crate::test_support::test_initialize();
    let spin = crate::representations::Bispinor {};
    use spenso::structure::representation::RepName;
    let slots = (0..4)
        .map(|i| mink!(4, { Atom::num(98800 + i) }))
        .collect::<Vec<_>>();
    let trace = trace!(spin.to_symbolic([4]); std::iter::once(gamma5!()).chain(slots.iter().map(|slot| gamma!(slot))));
    let source = SymbolicTensor::infer(trace).unwrap();
    let gamma_only = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        ..Default::default()
    };
    let result = source.simplify_algebra(&gamma_only).unwrap();
    assert!(
        result
            .expression
            .contains_symbol(*crate::epsilon::EPSILON_SYMBOL)
    );
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert!(Arc::ptr_eq(
        result.reduction_observations(),
        result
            .simplify_algebra(&gamma_only)
            .unwrap()
            .reduction_observations()
    ));
}

#[test]
fn a_capped_algebra_result_retains_its_unprocessed_trace() {
    let reps = crate::test_support::test_initialize();
    let expression =
        trace!(reps.bis4.to_symbolic([]); [gamma!(mink!(4, 98901)), gamma!(mink!(4, 98902))]);
    let source = SymbolicTensor::infer(expression).unwrap();
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        max_passes: Some(0),
        ..Default::default()
    };
    let capped = source.simplify_algebra(&settings).unwrap();
    assert_eq!(capped.reduction_status(), ReductionStatus::Capped);
    assert_eq!(source, capped);
    let finished = capped
        .simplify_algebra(&AlgebraSettings {
            max_passes: Some(16),
            ..settings
        })
        .unwrap();
    assert_eq!(
        finished.expression.expand(),
        Atom::num(4) * g!(mink!(4, 98901), mink!(4, 98902))
    );
}

#[test]
fn gamma_identity_contracts_its_prerequisite_without_opening_scalar_spectators() {
    let reps = crate::test_support::test_initialize();
    let left = reps.bis4.to_symbolic([Atom::num(99001)]);
    let middle = reps.bis4.to_symbolic([Atom::num(99002)]);
    let right = reps.bis4.to_symbolic([Atom::num(99003)]);
    let mu = mink!(4, 99004);
    let nu = mink!(4, 99005);
    let spectator = (Atom::one() + p!(mink!(4, 99006)) * q!(mink!(4, 99006))).pow(20);
    let source = SymbolicTensor::infer(
        &spectator * g!(&mu, &nu) * gamma!(&left, &middle, &mu) * gamma!(&middle, &right, &nu),
    )
    .unwrap();
    let result = source
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            contract: AlgebraContraction::None,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.structure, source.structure);
    // Clifford identity gamma^mu gamma_mu = 4 I; the independent scalar
    // coefficient retains even its internal contractions and power notation.
    assert_eq!(
        result.expression,
        Atom::num(4) * g!(&left, &right) * spectator
    );
}

#[test]
fn zero_budget_is_complete_when_no_permitted_work_exists() {
    crate::test_support::test_initialize();
    let constant = SymbolicTensor::infer(Atom::num(7)).unwrap();
    let settings = ContractSettings {
        max_passes: Some(0),
        ..Default::default()
    };
    assert_eq!(
        constant.contract(settings).unwrap().reduction_status(),
        ReductionStatus::Complete
    );
    let indexed = SymbolicTensor::infer(p!(mink!(4, 99101)) * q!(mink!(4, 99101))).unwrap();
    let result = indexed
        .contract(ContractSettings {
            representations: Some(&[]),
            ..settings
        })
        .unwrap();
    assert_eq!(result, indexed);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
}

#[test]
fn dummy_reservation_reuses_observed_domains_including_opaque_metadata() {
    use super::super::composition::DUMMY_RESERVATION_WALKS;
    use spenso::structure::abstract_index::AbstractIndex;
    crate::test_support::test_initialize();
    let mu = mink!(4, 99201);
    let hidden = mink!(4, 99202);
    let metadata = FunctionBuilder::new(symbol!("reservation_metadata"; Scalar))
        .add_arg(&hidden)
        .finish();
    let source = SymbolicTensor::infer(p!(&mu) * q!(&mu) * &metadata).unwrap();
    let observed = source.reduction_observations();
    assert_eq!(
        observed.reserved_indices,
        HashSet::from([AbstractIndex::Normal(99201), AbstractIndex::Normal(99202),])
    );
    DUMMY_RESERVATION_WALKS.with(|count| count.set(0));
    let _ = SymbolicTensor::reserved_dummies([&source, &source]);
    assert_eq!(DUMMY_RESERVATION_WALKS.with(|count| count.get()), 0);
    let replaced = source.with_identity_result(metadata, None).unwrap();
    assert_eq!(
        replaced.reduction_observations().reserved_indices,
        HashSet::from([AbstractIndex::Normal(99202)])
    );
}

#[test]
fn unchanged_kernel_regions_are_reused_after_an_independent_family_changes() {
    use observation::REUSED_KERNEL_REGIONS;
    crate::test_support::test_initialize();
    let left = crate::bis!(4, 99301);
    let right = crate::bis!(4, 99302);
    let matrix = gamma!(&left, &right, mink!(4, 99303));
    let epsilon = crate::epsilon!(mink!(2, 99304), mink!(2, 99305));
    let source = SymbolicTensor::infer(&matrix * epsilon.pow(2)).unwrap();
    REUSED_KERNEL_REGIONS.with(|count| count.set(0));
    let result = source
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            epsilon: true,
            ..Default::default()
        })
        .unwrap();
    // The two-dimensional determinant contracts to 2; its independent scalar
    // reduction does not change the open Clifford factor's inputs or boundary.
    let explicit = result.undo_chain().unwrap().expression.expand();
    assert!(
        (&explicit - Atom::num(2) * matrix).expand().is_zero(),
        "the explicit result is {explicit}"
    );
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert!(REUSED_KERNEL_REGIONS.with(|count| count.get()) > 0);
}

#[test]
fn changed_domain_regions_keep_independent_colour_work_settled() {
    use observation::{INITIAL_SCANS, REPLACEMENT_SCANS};

    crate::test_support::test_initialize();
    let label = spenso::index_symbol!("dirty_family::edge");
    let port = |index| mink!(4, { symbolica::function!(label, 99500, Atom::num(index)) });
    let color = crate::color_t!([8, 99501], [3, 99502], [3, 99503]);
    let gamma = gamma!(crate::bis!(4, 99504), crate::bis!(4, 99505), port(1));
    let source = SymbolicTensor::infer(&color * &gamma).unwrap();
    let replacement = SymbolicTensor::infer(
        &color * gamma!(crate::bis!(4, 99504), crate::bis!(4, 99505), port(2)),
    )
    .unwrap();
    let before = source.reduction_observations();
    replacement.observe_replacement(&source);
    let after = replacement.reduction_observations();
    assert!(
        before.candidates.complete,
        "compound labels do not hide tensor candidates"
    );
    assert!(before.candidates.traversal_complete);
    let scans = (
        INITIAL_SCANS.with(|count| count.get()),
        REPLACEMENT_SCANS.with(|count| count.get()),
    );
    assert_eq!(after.changed_identity_families(before), Some(GAMMA));
    assert_eq!(
        scans,
        (
            INITIAL_SCANS.with(|count| count.get()),
            REPLACEMENT_SCANS.with(|count| count.get()),
        ),
        "scheduling must consume observations without another candidate scan"
    );
    assert!(source.normalization_is_intrinsic());
    assert!(replacement.normalization_is_intrinsic());
}

#[test]
fn changed_connections_revisit_unchanged_identity_bodies() {
    crate::test_support::test_initialize();
    let a = mink!(4, 99511);
    let b = mink!(4, 99512);
    let c = mink!(4, 99513);
    let gamma = gamma!(crate::bis!(4, 99514), crate::bis!(4, 99515), &a);
    let source = SymbolicTensor::infer(&gamma * g!(&b, &c)).unwrap();
    let replacement = SymbolicTensor::infer(&gamma * g!(&a, &c)).unwrap();
    replacement.observe_replacement(&source);
    assert_eq!(
        replacement
            .reduction_observations()
            .changed_identity_families(source.reduction_observations()),
        Some(GAMMA),
        "a new Lorentz connection can expose work in the unchanged gamma"
    );
}

#[test]
fn newly_emitted_and_removed_identity_occurrences_invalidate_their_families() {
    let reps = crate::test_support::test_initialize();
    let color = crate::color_t!([8, 99521], [3, 99522], [3, 99523]);
    let first = trace!(reps.bis4.to_symbolic([]), gamma!(mink!(4, 99524)));
    let second = trace!(reps.bis4.to_symbolic([]), gamma!(mink!(4, 99525)));
    let source = SymbolicTensor::infer(&color * &first * &second).unwrap();
    let survivor = SymbolicTensor::infer(&color * &second).unwrap();
    survivor.observe_replacement(&source);
    assert_eq!(
        survivor
            .reduction_observations()
            .changed_identity_families(source.reduction_observations()),
        Some(GAMMA)
    );
    assert_ne!(
        survivor.reduction_observations().identity_candidates() & GAMMA,
        0
    );
    let epsilon = crate::epsilon!(
        mink!(4, 99526),
        mink!(4, 99527),
        mink!(4, 99528),
        mink!(4, 99529)
    );
    let emitted = SymbolicTensor::infer(&color * &second * epsilon).unwrap();
    emitted.observe_replacement(&survivor);
    let dirty = emitted
        .reduction_observations()
        .changed_identity_families(survivor.reduction_observations())
        .unwrap();
    assert_eq!(dirty & (GAMMA | EPSILON | COLOR), GAMMA | EPSILON);
}

#[test]
fn a_gamma_rewrite_does_not_recollect_an_unchanged_colour_sector() {
    let reps = crate::test_support::test_initialize();
    let color = SymbolicTensor::infer(crate::color_t!([8, 99531], [3, 99532], [3, 99533]))
        .unwrap()
        .simplify_algebra(&AlgebraSettings {
            color: Some(ColorSimplifySettings::default()),
            contract: AlgebraContraction::None,
            ..Default::default()
        })
        .unwrap()
        .expression;
    let word = trace!(reps.bis4.to_symbolic([]);
        (0..8).map(|index| gamma!(mink!(4, { Atom::num(99540 + index) }))));
    let source = SymbolicTensor::infer(&color * &word).unwrap();
    IDENTITY_KERNEL_CALLS.with(|calls| calls.set([0; 3]));
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        color: Some(ColorSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let result = source.simplify_algebra(&settings).unwrap();
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(IDENTITY_KERNEL_CALLS.with(|calls| calls.get()[1]), 1);
    let gamma = SymbolicTensor::infer(word)
        .unwrap()
        .simplify_algebra(&AlgebraSettings {
            color: None,
            ..settings
        })
        .unwrap();
    assert!(
        (&result.expression - color * gamma.expression)
            .expand()
            .is_zero()
    );
}

#[test]
fn trusted_replacements_reuse_unchanged_opaque_boundary_interfaces() {
    crate::test_support::test_initialize();
    let mu = mink!(4, 99401);
    let nu = mink!(4, 99402);
    let spectator = (p!(mink!(4, 99403)) * q!(mink!(4, 99403))
        + p!(mink!(4, 99404)) * q!(mink!(4, 99404)) * Atom::num(3))
    .pow(12);
    let source = SymbolicTensor::infer(g!(&mu, &nu) * p!(&nu) * &spectator).unwrap();
    source.reduction_observations();
    source.shallow_graph().unwrap();
    let interfaces = source.proofs.leaf_interfaces.get().unwrap();
    let prior = interfaces.lock().unwrap().clone();
    assert!(prior.contains_key(&spectator));

    // This is the metric substitution's exact result. Its spectator keeps its
    // admitted scalar interface even though its interior still has dummy pairs.
    let rewritten = source
        .with_identity_result(p!(&mu) * &spectator, None)
        .unwrap();
    assert_ne!(rewritten.expression, source.expression);
    assert!(Arc::ptr_eq(
        interfaces,
        rewritten.proofs.leaf_interfaces.get().unwrap()
    ));
    rewritten.shallow_graph().unwrap();
    let current = interfaces.lock().unwrap();
    assert_eq!(current.get(&spectator), prior.get(&spectator));
    assert_eq!(current.len(), prior.len() + 1);
    assert!(current.contains_key(&p!(&mu)));
}
