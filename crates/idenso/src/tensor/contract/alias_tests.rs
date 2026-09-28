use super::*;
use spenso::structure::{abstract_index::AbstractIndex, slot::IsAbstractSlot};
use spenso::{g, mink};
use symbolica::atom::FunctionBuilder;

fn vector(head: &str, slot: &Atom) -> Atom {
    FunctionBuilder::new(spenso::network::tags::SPENSO_TAG.rank_one_tensor_symbol(head))
        .add_arg(slot)
        .finish()
}

#[test]
fn ordered_fallback_registers_alias_ports_after_intrinsic_refusal() {
    crate::test_support::test_initialize();
    let a = mink!(4, 97111);
    let b = mink!(4, 97113);
    let body = SymbolicTensor::infer(vector("ordered_alias_fallback_p", &a)).unwrap();
    let handle = body.alias_handle().unwrap();
    // The callback deliberately declines intrinsic collection. Its no-op
    // behavior does not authorize dropping the checked ordered fallback.
    let callback = spenso::tensor_symbol!("ordered_alias_fallback_callback", norm = |_, _| {});
    let spectator = FunctionBuilder::new(callback).finish();
    let expression = g!(a, b.clone()) * &handle.expression * &spectator;
    assert!(
        SlotContraction::new()
            .contract_factorized(expression.as_view(), None, true)
            .is_none()
    );
    let root = SymbolicTensor::infer(expression).unwrap();
    let value = Arc::new(root.with_aliases([(handle, body)]).unwrap());
    let result = value.contract(Default::default()).unwrap();
    assert_eq!(
        result.resolved().unwrap().expression,
        vector("ordered_alias_fallback_p", &b) * spectator
    );
    assert_eq!(result.root().structure, value.root().structure);
    assert_eq!(
        result
            .contract(Default::default())
            .unwrap()
            .resolved()
            .unwrap(),
        result.resolved().unwrap()
    );
}

#[test]
fn alias_contraction_registers_nested_literal_labels_and_reuses_completion() {
    crate::test_support::test_initialize();
    let a = mink!(4, 91301);
    let b = mink!(4, 91303);
    let body =
        SymbolicTensor::infer(vector("alias_contract_p", &a) + vector("alias_contract_q", &a))
            .unwrap();
    let inner = body.alias_handle().unwrap();
    let outer_body = inner.clone();
    let outer = outer_body.alias_handle().unwrap();
    let root = SymbolicTensor::infer(g!(a.clone(), b.clone()) * outer.expression.clone()).unwrap();
    let value = Arc::new(
        root.with_aliases([(inner, body), (outer, outer_body)])
            .unwrap(),
    );
    let expected = value
        .resolved()
        .unwrap()
        .contract(Default::default())
        .unwrap()
        .expanded()
        .unwrap();
    let result = value.contract(Default::default()).unwrap();
    assert_eq!(result.expanded().unwrap(), expected);
    assert!(result.aliases().unwrap().len() >= 4);
    assert!(Arc::ptr_eq(
        &result,
        &result.contract(Default::default()).unwrap()
    ));
}

#[test]
fn literal_registration_uses_the_exact_source_not_the_alias_owner() {
    crate::test_support::test_initialize();
    let a = mink!(4, 91401);
    let b = mink!(4, 91403);
    let c = mink!(4, 91405);
    let body = SymbolicTensor::infer(vector("alias_exact_p", &a)).unwrap();
    let handle = body.alias_handle().unwrap();
    let other = handle
        .reindex_interface_ports(&HashMap::from([(0, AbstractIndex::from(91405))]))
        .unwrap();
    let foreign = SymbolicTensor::infer(vector("alias_exact_q", &c)).unwrap();
    let root = SymbolicTensor::infer(g!(a, b.clone()) * handle.expression.clone()).unwrap();
    let value = Arc::new(
        root.with_aliases([(handle, body), (other.clone(), foreign.clone())])
            .unwrap(),
    );
    let result = value.contract(Default::default()).unwrap();
    assert_eq!(
        result.resolved().unwrap().expression,
        vector("alias_exact_p", &b)
    );
    let definitions = result.aliases().unwrap();
    assert!(
        definitions
            .iter()
            .any(|(handle, body)| handle == &other && body == &foreign)
    );
}

#[test]
fn alias_contraction_binds_compact_vectors_into_literal_definitions() {
    crate::test_support::test_initialize();
    let a = mink!(4, 91411);
    let b = mink!(4, 91413);
    let body = SymbolicTensor::infer(
        FunctionBuilder::new(spenso::tensor_symbol!("alias_bound_T"))
            .add_arg(&a)
            .add_arg(&b)
            .finish(),
    )
    .unwrap();
    let handle = body.alias_handle().unwrap();
    let momentum = vector("alias_bound_p", &a);
    let source = SymbolicTensor::infer(&momentum * &handle.expression).unwrap();
    let value = Arc::new(source.with_aliases([(handle, body.clone())]).unwrap());
    let result = value.contract(Default::default()).unwrap();
    let compact = vector("alias_bound_p", &mink!(4));
    let expected = body
        .expression
        .replace(a.to_pattern())
        .with(compact.to_pattern());
    assert_eq!(result.resolved().unwrap().expression, expected);
    assert_eq!(
        result.root().structure.logical_slots(),
        body.structure.logical_slots()[1..]
    );
    assert!(Arc::ptr_eq(
        &result,
        &result.contract(Default::default()).unwrap()
    ));
}

#[test]
fn literal_registration_uses_separate_reversed_open_layouts() {
    crate::test_support::test_initialize();
    let body = SymbolicTensor::infer(
        FunctionBuilder::new(spenso::tensor_symbol!("alias_open_layout_T"))
            .add_arg(mink!(4))
            .add_arg(spenso::euc!(6))
            .finish(),
    )
    .unwrap();
    let mut handle = body.alias_handle().unwrap();
    handle.structure =
        PartialStructure::from_logical_slots(handle.structure.logical_slots().into_iter().rev());
    let AtomView::Fun(function) = handle.expression.as_view() else {
        panic!("alias handle");
    };
    let arguments = function
        .iter()
        .map(|argument| argument.to_owned())
        .collect::<Vec<_>>();
    let compact = vector("alias_open_layout_p", &mink!(4));
    let target = FunctionBuilder::new(function.get_symbol())
        .add_arg(&arguments[0])
        .add_arg(&compact)
        .add_arg(&arguments[2])
        .finish();
    let mut definitions = HashMap::from([(handle.expression.clone(), (handle.clone(), body))]);
    SymbolicTensor::register_literal_use(&handle.expression, &target, &mut definitions).unwrap();
    let (next_handle, next_body) = &definitions[&target];
    assert_eq!(
        next_body.expression,
        FunctionBuilder::new(spenso::tensor_symbol!("alias_open_layout_T"))
            .add_arg(&compact)
            .add_arg(spenso::euc!(6))
            .finish()
    );
    assert_eq!(
        next_handle.structure.logical_slots(),
        next_body.structure.logical_slots()
    );
    assert_eq!(next_body.structure.logical_slots().len(), 1);
    assert_eq!(
        next_body.structure.logical_slots()[0].rep().to_symbolic([]),
        spenso::euc!(6)
    );
}

#[test]
fn alias_port_rewrite_checks_callback_rank_loss_before_registry_publication() {
    crate::test_support::test_initialize();
    let a = mink!(4, 91501);
    let b = mink!(4, 91503);
    let target = b.clone();
    let head = spenso::tensor_symbol!(
        "alias_rankloss",
        norm = move |value, out| {
            if let AtomView::Fun(function) = value
                && function.iter().any(|argument| argument == target.as_view())
            {
                **out = Atom::one();
            }
        }
    );
    let body = SymbolicTensor::infer(FunctionBuilder::new(head).add_arg(&a).finish()).unwrap();
    let handle = body.alias_handle().unwrap();
    let root = SymbolicTensor::infer(g!(a, b) * handle.expression.clone()).unwrap();
    let value = Arc::new(root.with_aliases([(handle, body)]).unwrap());
    let original = value.expression.clone();
    assert!(value.contract(Default::default()).is_err());
    assert_eq!(value.expression.get_root(), original.get_root());
    assert_eq!(value.expression.get_aliases(), original.get_aliases());
}

#[test]
fn alias_contraction_exposes_open_metric_products_and_powers() {
    use spenso::structure::slot::ParseableAind;
    crate::test_support::test_initialize();
    let head = symbolica::symbol!(
        "alias_frontier::port",
        tags = [spenso::network::tags::SPENSO_TAG.index.clone()]
    );
    let scope = symbolica::symbol!("alias_frontier::bra");
    for named in [false, true] {
        for scoped in [false, true] {
            let indices = (0..4)
                .map(|i| {
                    let index = if named {
                        AbstractIndex::Named(head.into(), 94000 + i, 0)
                    } else {
                        AbstractIndex::from(94000 + i)
                    };
                    let index = if scoped { index.scoped(scope) } else { index };
                    mink!(4, index.to_atom())
                })
                .collect::<Vec<_>>();
            let [a, b, c, d] = indices.as_slice() else {
                unreachable!()
            };
            let body = SymbolicTensor::infer(
                Atom::num(4) * (g!(a, b) * g!(c, d) - g!(a, c) * g!(b, d) + g!(a, d) * g!(b, c)),
            )
            .unwrap();
            assert_eq!(body.structure.logical_slots().len(), 4);
            let first = body.alias_handle().unwrap();
            let second = body.alias_handle().unwrap();
            for root in [
                &first.expression * &first.expression,
                &first.expression * &second.expression,
            ] {
                let root = SymbolicTensor::infer(root).unwrap();
                let value = Arc::new(
                    root.with_aliases([
                        (first.clone(), body.clone()),
                        (second.clone(), body.clone()),
                    ])
                    .unwrap(),
                );
                let result = value.contract(Default::default()).unwrap();
                assert!(result.contraction_complete());
                assert_eq!(result.resolved().unwrap().expression, Atom::num(640));
                assert!(Arc::ptr_eq(
                    &result,
                    &result.contract(Default::default()).unwrap()
                ));
            }
            assert_eq!(
                SymbolicTensor::infer(body.expression.pow(2))
                    .unwrap()
                    .contract(Default::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .expression,
                Atom::num(640)
            );
        }
    }
}

#[test]
fn alias_contraction_exposes_two_independent_dirac_traces() {
    crate::test_support::test_initialize();
    let head = symbolica::symbol!(
        "alias_trace_frontier::port",
        tags = [spenso::network::tags::SPENSO_TAG.index.clone()]
    );
    let scope = symbolica::symbol!("alias_trace_frontier::bra");
    use spenso::structure::slot::ParseableAind;
    for mode in 0..3 {
        let slots = (0..4)
            .map(|i| {
                let index = match mode {
                    0 => AbstractIndex::from(95100 + i),
                    1 => AbstractIndex::Named(head.into(), 95100 + i, 0),
                    _ => AbstractIndex::Named(head.into(), 95100 + i, 0).scoped(scope),
                };
                mink!(4, index.to_atom())
            })
            .collect::<Vec<_>>();
        let trace = |offset: usize| {
            Atom::mul_many((0..4).map(|i| {
                crate::gamma!(
                    [4usize, offset + i],
                    [4usize, offset + (i + 1) % 4],
                    &slots[i]
                )
            }))
        };
        let source = SymbolicTensor::infer(trace(95200) * trace(95300)).unwrap();
        let result = source.simplify_gamma(Default::default()).unwrap();
        assert!(result.contraction_complete());
        assert_eq!(result.resolved().unwrap().expression, Atom::num(640));
    }
}

#[test]
fn alias_contraction_exposure_preserves_disconnected_callbacks_and_spectators() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    crate::test_support::test_initialize();
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&calls);
    let callback = spenso::vector_symbol!(
        "alias_frontier_unused_callback",
        norm = move |_, _| {
            observed.fetch_add(1, Ordering::Relaxed);
        }
    );
    let unused = SymbolicTensor::infer(
        FunctionBuilder::new(callback)
            .add_arg(mink!(4, 95401))
            .finish(),
    )
    .unwrap();
    let unused_handle = unused.alias_handle().unwrap();
    let body = SymbolicTensor::infer(g!(mink!(4, 95403), mink!(4, 95405))).unwrap();
    let first = body.alias_handle().unwrap();
    let second = body.alias_handle().unwrap();
    let spectator = symbolica::parse_lit!((alias_frontier_x + alias_frontier_y) ^ 20);
    let root = SymbolicTensor::infer(&spectator * &first.expression * &second.expression).unwrap();
    let value = Arc::new(
        root.with_aliases([
            (first, body.clone()),
            (second, body),
            (unused_handle.clone(), unused.clone()),
        ])
        .unwrap(),
    );
    calls.store(0, Ordering::Relaxed);
    let result = value.contract(Default::default()).unwrap();
    assert!(result.contraction_complete());
    assert_eq!(calls.load(Ordering::Relaxed), 0);
    assert_eq!(
        result.resolved().unwrap().expression,
        Atom::num(4) * spectator
    );
    assert!(result.aliases().unwrap().contains(&(unused_handle, unused)));
}

#[test]
fn alias_contraction_exposure_does_not_restart_an_incomplete_frontier() {
    crate::test_support::test_initialize();
    let p = spenso::vector_symbol!("alias_frontier_budget_p");
    let q = spenso::vector_symbol!("alias_frontier_budget_q");
    let t = spenso::tensor_symbol!("alias_frontier_budget_t");
    let slots = (0..9)
        .map(|i| mink!(4, Atom::num(95500 + i)))
        .collect::<Vec<_>>();
    let weight = Atom::var(symbolica::symbol!("alias_frontier_budget_weight"));
    let expression = FunctionBuilder::new(t).add_args(&slots).finish()
        * Atom::mul_many(slots.iter().map(|slot| {
            FunctionBuilder::new(p).add_arg(slot).finish()
                + &weight * FunctionBuilder::new(q).add_arg(slot).finish()
        }));
    let body = SymbolicTensor::infer(expression).unwrap();
    let direct = body.contract(Default::default()).unwrap();
    assert!(!direct.contraction_complete());
    let handle = body.alias_handle().unwrap();
    let value = Arc::new(handle.clone().with_aliases([(handle, body)]).unwrap());
    let result = value.contract(Default::default()).unwrap();
    assert!(!result.contraction_complete());
    assert_eq!(
        result.aliases().unwrap().len(),
        direct.aliases().unwrap().len() + 1
    );
    assert_eq!(
        result.resolved().unwrap().expression,
        direct.resolved().unwrap().expression
    );
}

#[test]
fn alias_contraction_exposure_reports_a_capped_cross_definition_frontier() {
    crate::test_support::test_initialize();
    let p = spenso::vector_symbol!("alias_exposure_budget_p");
    let q = spenso::vector_symbol!("alias_exposure_budget_q");
    let weight = Atom::var(symbolica::symbol!("alias_exposure_budget_weight"));
    let body = SymbolicTensor::infer(Atom::mul_many((0..9).map(|i| {
        let slot = mink!(4, Atom::num(95600 + i));
        FunctionBuilder::new(p).add_arg(&slot).finish()
            + &weight * FunctionBuilder::new(q).add_arg(&slot).finish()
    })))
    .unwrap();
    assert!(
        body.contract(Default::default())
            .unwrap()
            .contraction_complete()
    );
    let handle = body.alias_handle().unwrap();
    let root = SymbolicTensor::infer(handle.expression.pow(2)).unwrap();
    let interface = root.structure.clone();
    let value = Arc::new(root.with_aliases([(handle, body)]).unwrap());
    let result = value.contract(Default::default()).unwrap();
    assert!(!result.contraction_complete());
    assert_eq!(result.root().structure, interface);
    assert!(!result.aliases().unwrap().is_empty());
}

#[test]
fn alias_contraction_keeps_disjoint_open_definitions_factored() {
    crate::test_support::test_initialize();
    let a = mink!(4, 94801);
    let b = mink!(4, 94803);
    let c = mink!(4, 94805);
    let d = mink!(4, 94807);
    let e = mink!(4, 94809);
    let f = mink!(4, 94811);
    let h = mink!(4, 94813);
    let metric_body =
        SymbolicTensor::infer(g!(&a, &b) * (g!(&c, &d) * g!(&e, &f) - g!(&c, &e) * g!(&d, &f)))
            .unwrap();
    let vector_body =
        SymbolicTensor::infer(vector("alias_disjoint_p", &h) + vector("alias_disjoint_q", &h))
            .unwrap();
    let metric = metric_body.alias_handle().unwrap();
    let vector = vector_body.alias_handle().unwrap();
    let root = SymbolicTensor::infer(&metric.expression * &vector.expression).unwrap();
    let value = Arc::new(
        root.with_aliases([(metric, metric_body), (vector, vector_body)])
            .unwrap(),
    );
    let result = value.contract(Default::default()).unwrap();
    assert!(result.contraction_complete());
    assert_eq!(result.root(), value.root());
    let definitions = |value: &Arc<SymbolicTensor<AliasInterfaces, AliasedAtom>>| {
        value
            .aliases()
            .unwrap()
            .into_iter()
            .map(|(handle, body)| (handle.expression, body))
            .collect::<HashMap<_, _>>()
    };
    assert_eq!(definitions(&result), definitions(&value));
    assert!(Arc::ptr_eq(
        &result,
        &result.contract(Default::default()).unwrap()
    ));
}

#[test]
fn alias_contraction_registers_same_stage_templates_before_publication() {
    crate::test_support::test_initialize();
    let a = mink!(4, 97201);
    let b = mink!(4, 97203);
    let c = mink!(4, 97205);
    let compact = mink!(4);
    let p = |slot: &Atom| vector("alias_stage_p", slot);
    let q = |slot: &Atom| vector("alias_stage_q", slot);
    let r = |slot: &Atom| vector("alias_stage_r", slot);
    let t = |slot: &Atom| vector("alias_stage_t", slot);
    let literal = |id: i64, port: Atom| {
        SymbolicTensor::infer(
            FunctionBuilder::new(*crate::tensor::aliases::TENSOR_ALIAS_SYMBOL)
                .add_arg(Atom::num(id))
                .add_arg(port)
                .finish(),
        )
        .unwrap()
    };
    for (source_id, outer_id) in [(97210, 97211), (97211, 97210)] {
        for conflict in [false, true] {
            let body = SymbolicTensor::infer(p(&a) * q(&a) * p(&b) + p(&a).pow(2) * q(&b)).unwrap();
            let source = literal(source_id, b.clone());
            let target = literal(source_id, r(&compact));
            let mut definitions =
                HashMap::from([(source.expression.clone(), (source.clone(), body))]);
            SymbolicTensor::register_literal_use(
                &source.expression,
                &target.expression,
                &mut definitions,
            )
            .unwrap();
            if conflict {
                let body = &mut definitions.get_mut(&target.expression).unwrap().1;
                *body = SymbolicTensor::infer(&body.expression + Atom::one()).unwrap();
            }
            let outer = literal(outer_id, c.clone());
            let outer_body = SymbolicTensor::infer(r(&b) * &source.expression * t(&c)).unwrap();
            let root =
                SymbolicTensor::infer(&target.expression * t(&c) + &outer.expression).unwrap();
            definitions.insert(outer.expression.clone(), (outer, outer_body));
            let value = Arc::new(root.with_aliases(definitions.into_values()).unwrap());
            let result = value.contract(Default::default());
            if conflict {
                assert!(
                    result
                        .unwrap_err()
                        .to_string()
                        .contains("conflicting relabelled alias")
                );
                continue;
            }
            let result = result.unwrap();
            assert!(result.contraction_complete());
            assert_eq!(result.root().structure, value.root().structure);
            assert_eq!(
                result.resolved().unwrap().expression,
                Atom::num(2)
                    * (g!(p(&compact), q(&compact)) * g!(p(&compact), r(&compact))
                        + g!(p(&compact), p(&compact)) * g!(q(&compact), r(&compact)))
                    * t(&c)
            );
            assert!(Arc::ptr_eq(
                &result,
                &result.contract(Default::default()).unwrap()
            ));
        }
    }
}

#[test]
fn alias_contraction_exposes_pairs_with_exact_imaginary_weights() {
    crate::test_support::test_initialize();
    let a = mink!(4, 97301);
    let b = mink!(4, 97303);
    let c = mink!(4, 97305);
    let d = mink!(4, 97307);
    let first_body = SymbolicTensor::infer(g!(&a, &b)).unwrap();
    let second_body = SymbolicTensor::infer(g!(&c, &d)).unwrap();
    let handles = [
        first_body.alias_handle().unwrap(),
        first_body.alias_handle().unwrap(),
        second_body.alias_handle().unwrap(),
        second_body.alias_handle().unwrap(),
    ];
    let denominator = Atom::var(symbolica::symbol!("alias_imaginary_denominator"));
    let root = SymbolicTensor::infer(
        Atom::i() * &handles[0].expression * &handles[1].expression / &denominator
            + Atom::num(2) * Atom::i() * &handles[2].expression * &handles[3].expression
                / denominator.pow(2),
    )
    .unwrap();
    let value = Arc::new(
        root.with_aliases(handles.into_iter().zip([
            first_body.clone(),
            first_body,
            second_body.clone(),
            second_body,
        ]))
        .unwrap(),
    );
    let result = value.contract(Default::default()).unwrap();
    assert!(result.contraction_complete());
    let expected = SymbolicTensor::infer(
        Atom::i() * Atom::num(4) / &denominator + Atom::i() * Atom::num(8) / denominator.pow(2),
    )
    .unwrap();
    // The finite metric traces leave only this scalar denominator algebra.
    assert!(
        (result.resolved().unwrap().expression - expected.expression)
            .together()
            .is_zero()
    );
    assert!(Arc::ptr_eq(
        &result,
        &result.contract(Default::default()).unwrap()
    ));
}

#[test]
fn deferred_epsilon_alias_frontier_finishes_in_the_shared_fixed_point() {
    crate::test_support::test_initialize();
    let a = mink!(4, 97401);
    let b = mink!(4, 97403);
    let compact = mink!(4);
    let body = SymbolicTensor::infer(g!(&a, &b)).unwrap();
    let first = body.alias_handle().unwrap();
    let second = body.alias_handle().unwrap();
    let third = body.alias_handle().unwrap();
    let epsilon = crate::epsilon!(
        &a,
        &b,
        vector("alias_deferred_p", &compact),
        vector("alias_deferred_q", &compact)
    );
    let root =
        SymbolicTensor::infer(&first.expression * &second.expression + epsilon * &third.expression)
            .unwrap();
    let parts = root.contract_parts(Default::default()).unwrap();
    assert!(parts.status == ContractionStatus::Deferred);
    let value = Arc::new(
        root.with_aliases([(first, body.clone()), (second, body.clone()), (third, body)])
            .unwrap(),
    );
    let result = value.simplify_gamma(Default::default()).unwrap();
    assert!(result.contraction_complete());
    assert_eq!(result.resolved().unwrap().expression, Atom::num(4));
    assert!(Arc::ptr_eq(
        &result,
        &result.contract(Default::default()).unwrap()
    ));
}

#[test]
fn capped_alias_frontier_before_first_emission_is_not_reentered() {
    crate::test_support::test_initialize();
    let p = spenso::vector_symbol!("alias_initial_cap_p");
    let q = spenso::vector_symbol!("alias_initial_cap_q");
    let slot = mink!(4, 97501);
    let expression = Atom::add_many((0..256).map(|position| {
        FunctionBuilder::new(p)
            .add_arg(Atom::num(position))
            .add_arg(&slot)
            .finish()
            * FunctionBuilder::new(q)
                .add_arg(Atom::num(position))
                .add_arg(&slot)
                .finish()
    }));
    let body = SymbolicTensor::infer(expression).unwrap();
    let direct = body.contract_parts(Default::default()).unwrap();
    assert!(direct.status == ContractionStatus::Capped);
    assert!(direct.aliases.is_empty());
    assert_eq!(direct.root, body);
    let handle = body.alias_handle().unwrap();
    let value = Arc::new(handle.clone().with_aliases([(handle, body)]).unwrap());
    let result = value.contract(Default::default()).unwrap();
    assert!(!result.contraction_complete());
    assert_eq!(result.root(), value.root());
    assert_eq!(result.aliases().unwrap(), value.aliases().unwrap());
}

#[test]
fn epsilon_context_exposes_symmetric_and_epsilon_aliases() {
    crate::test_support::test_initialize();
    let a = mink!(4, 97601);
    let b = mink!(4, 97603);
    let compact = mink!(4);
    let p = |slot: &Atom| vector("epsilon_context_p", slot);
    let q = |slot: &Atom| vector("epsilon_context_q", slot);
    let r = vector("epsilon_context_r", &compact);
    let s = vector("epsilon_context_s", &compact);
    let symmetric = SymbolicTensor::infer(p(&a) * q(&b) + p(&b) * q(&a)).unwrap();
    let handle = symmetric.alias_handle().unwrap();
    let epsilon = crate::epsilon!(&a, &b, &r, &s);
    let root = SymbolicTensor::infer(&epsilon * &handle.expression).unwrap();
    let value = Arc::new(root.with_aliases([(handle, symmetric)]).unwrap());
    let result = value.simplify_epsilon().unwrap();
    assert!(result.contraction_complete());
    assert!(result.resolved().unwrap().expression.is_zero());

    let c = mink!(4, 97605);
    let d = mink!(4, 97607);
    let body = SymbolicTensor::infer(crate::epsilon!(&a, &b, &c, &d)).unwrap();
    let first = body.alias_handle().unwrap();
    let second = body.alias_handle().unwrap();
    let root = SymbolicTensor::infer(&first.expression * &second.expression).unwrap();
    let value = Arc::new(
        root.with_aliases([(first, body.clone()), (second, body)])
            .unwrap(),
    );
    let result = value.simplify_epsilon().unwrap();
    assert!(result.contraction_complete());
    assert_eq!(result.resolved().unwrap().expression, Atom::num(24));
}

#[test]
fn epsilon_stage_does_not_collect_metrics_from_opaque_index_labels() {
    use spenso::structure::slot::ParseableAind;
    crate::test_support::test_initialize();
    let label = symbolica::symbol!(
        "epsilon_scope_index",
        tags = [spenso::network::tags::SPENSO_TAG.index.clone()]
    );
    let a = mink!(4, AbstractIndex::Named(label.into(), 1, 0).to_atom());
    let b = mink!(4, AbstractIndex::Named(label.into(), 2, 0).to_atom());
    let body = SymbolicTensor::infer(g!(&a, &b)).unwrap();
    let first = body.alias_handle().unwrap();
    let second = body.alias_handle().unwrap();
    let root = SymbolicTensor::infer(&first.expression * &second.expression).unwrap();
    let value = Arc::new(
        root.with_aliases([(first, body.clone()), (second, body)])
            .unwrap(),
    );
    let result = value
        .simplify(&crate::tensor::simplification::SimplifySettings {
            metrics: false,
            epsilon: true,
            ..Default::default()
        })
        .unwrap();
    assert!(Arc::ptr_eq(&value, &result));
}
