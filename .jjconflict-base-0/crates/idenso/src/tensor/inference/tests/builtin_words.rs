use super::*;

#[test]
fn admitted_slot_metadata_does_not_become_a_standalone_tensor() {
    crate::test_support::test_initialize();
    let rep = LibraryRep::from(ColorAdjoint {}).new_rep(8);
    let scope = symbolica::symbol!("slot_metadata_scope");
    let index = AbstractIndex::Symbol(CS.d.into());
    for index in [index, index.scoped(scope)] {
        let slot = rep.slot::<AbstractIndex, _>(index);
        let expression = FunctionBuilder::new(CS.f)
            .add_arg(slot.to_atom())
            .add_arg(
                rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99501))
                    .to_atom(),
            )
            .add_arg(
                rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99502))
                    .to_atom(),
            )
            .finish();
        let value = SymbolicTensor::infer(expression.clone()).unwrap();
        assert_eq!(value.expression, expression);
        assert_eq!(value.structure.logical_slots().len(), 3);
        assert!(
            value
                .structure
                .logical_slots()
                .contains(&rep.slot(PartialIndex::Explicit(index)))
        );
    }
    assert!(SymbolicTensor::infer(Atom::var(CS.d)).is_err());

    let nested = FunctionBuilder::new(symbolica::symbol!("slot_metadata_nested"))
        .add_arg(Atom::var(CS.d))
        .finish();
    let malformed_scope = FunctionBuilder::new(AbstractIndex::scope_symbol())
        .add_arg(Atom::var(scope))
        .add_arg(&nested)
        .finish();
    for index in [nested, malformed_scope] {
        let invalid_slot = FunctionBuilder::new(rep.rep.symbol())
            .add_arg(Atom::num(8))
            .add_arg(index)
            .finish();
        let expression = FunctionBuilder::new(spenso::tensor_symbol!("slot_metadata_tensor"))
            .add_arg(invalid_slot)
            .finish();
        assert!(SymbolicTensor::infer(expression).is_err());
    }
}

#[test]
fn placeholder_gamma_words_use_the_factory_logical_order() {
    crate::test_support::test_initialize();
    for dim in [
        Dimension::Concrete(4),
        Dimension::from(symbolica::symbol!("builtin_word_D")),
    ] {
        let mink = Minkowski {}.new_rep(dim);
        let mu = mink
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(99101))
            .to_atom();
        let word = FunctionBuilder::new(AGS.gamma)
            .add_arg(Atom::var(SPENSO_TAG.chain_in))
            .add_arg(Atom::var(SPENSO_TAG.chain_out))
            .add_arg(&mu)
            .finish();
        let trace = FunctionBuilder::new(SPENSO_TAG.trace)
            .add_arg(Bispinor {}.new_rep(4).to_symbolic([]))
            .add_arg(&word)
            .add_arg(&word)
            .finish();
        let value = SymbolicTensor::infer(trace.clone()).unwrap();
        assert!(value.is_scalar());
        assert_eq!(value.expression, trace);
        let chain = FunctionBuilder::new(SPENSO_TAG.chain)
            .add_arg(
                Bispinor {}
                    .new_rep(4)
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(99102))
                    .to_atom(),
            )
            .add_arg(
                Bispinor {}
                    .new_rep(4)
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(99103))
                    .to_atom(),
            )
            .add_arg(&word)
            .finish();
        let value = SymbolicTensor::infer(chain).unwrap();
        assert_eq!(value.structure.logical_slots().len(), 3);
        assert!(
            SymbolicTensor::infer(word).is_err(),
            "unbound placeholders remain invalid"
        );
    }
    let invalid_word = FunctionBuilder::new(AGS.gamma)
        .add_arg(Atom::var(SPENSO_TAG.chain_in))
        .add_arg(Atom::var(SPENSO_TAG.chain_out))
        .add_arg(
            ExtendibleReps::EUCLIDEAN
                .new_rep(4)
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(99104))
                .to_atom(),
        )
        .finish();
    let invalid_trace = FunctionBuilder::new(SPENSO_TAG.trace)
        .add_arg(Bispinor {}.new_rep(4).to_symbolic([]))
        .add_arg(invalid_word)
        .finish();
    assert!(matches!(
        SymbolicTensor::infer(invalid_trace),
        Err(TensorInferenceError::InvalidBuiltinSignature {
            factory: "gamma",
            ..
        })
    ));
    let standalone =
        SymbolicTensor::from_signature(&AGS.gamma_strct(Dimension::Concrete(4))).unwrap();
    assert_eq!(
        SymbolicTensor::infer(standalone.expression.clone())
            .unwrap()
            .structure,
        standalone.structure
    );
    assert!(
        SymbolicTensor::infer(
            FunctionBuilder::new(AGS.gamma)
                .add_arg(Atom::one())
                .finish()
        )
        .is_err()
    );
}

#[test]
fn placeholder_color_words_validate_the_existing_channel_signature() {
    crate::test_support::test_initialize();
    let adjoint = ColorAdjoint {}.new_rep(8);
    let a = adjoint
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99201))
        .to_atom();
    let word = FunctionBuilder::new(CS.t)
        .add_arg(&a)
        .add_arg(Atom::var(SPENSO_TAG.chain_in))
        .add_arg(Atom::var(SPENSO_TAG.chain_out))
        .finish();
    let trace = FunctionBuilder::new(SPENSO_TAG.trace)
        .add_arg(ColorFundamental {}.new_rep(3).to_symbolic([]))
        .add_arg(&word)
        .add_arg(&word)
        .finish();
    let value = SymbolicTensor::infer(trace.clone()).unwrap();
    assert!(value.is_scalar());
    assert_eq!(value.expression, trace);
    let wrong_channel = FunctionBuilder::new(SPENSO_TAG.trace)
        .add_arg(Bispinor {}.new_rep(4).to_symbolic([]))
        .add_arg(&word)
        .finish();
    assert!(matches!(
        SymbolicTensor::infer(wrong_channel),
        Err(TensorInferenceError::InvalidBuiltinSignature { factory: "t", .. })
    ));
    assert!(SymbolicTensor::infer(word).is_err());
}

#[test]
fn scalar_arithmetic_and_color_invariants_have_an_empty_interface() {
    crate::test_support::test_initialize();
    let invariants = FunctionBuilder::new(CS.cas)
        .add_arg(Atom::num(2))
        .add_arg(ColorAdjoint {}.new_rep(8).to_symbolic([]))
        .finish()
        * FunctionBuilder::new(CS.idx)
            .add_arg(Atom::num(2))
            .add_arg(ColorFundamental {}.new_rep(3).to_symbolic([]))
            .finish();
    let x = Atom::var(symbolica::symbol!("builtin_scalar_x"));
    for expression in [
        Atom::Zero,
        Atom::one(),
        invariants.clone(),
        (&x + invariants).pow(3),
    ] {
        let value = SymbolicTensor::infer(expression.clone()).unwrap();
        assert!(value.is_scalar());
        assert_eq!(value.expression, expression);
        assert_eq!(
            SymbolicTensor::infer(value.expression.clone()).unwrap(),
            value
        );
    }
    let bare_tensor = Atom::var(spenso::tensor_symbol!("builtin_scalar_invalid_tensor"));
    assert!(SymbolicTensor::infer(bare_tensor).is_err());
    assert!(SymbolicTensor::infer(Atom::var(SPENSO_TAG.chain_in)).is_err());
}

#[test]
fn formal_gamma_and_sigma_words_inherit_atomic_spin_channels() {
    crate::test_support::test_initialize();
    let mink = Minkowski {}.new_rep(4);
    let [mu, nu] = [99301, 99302].map(|index| {
        mink.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom()
    });
    for dimension in [
        Dimension::Concrete(6),
        Dimension::from(symbolica::symbol!("formal_spin_Ns")),
    ] {
        let spin = Bispinor {}.new_rep(dimension);
        for (head, lorentz) in [
            (AGS.gamma, vec![mu.clone()]),
            (AGS.sigma, vec![mu.clone(), nu.clone()]),
        ] {
            let spin_ports = [
                Atom::var(SPENSO_TAG.chain_in),
                Atom::var(SPENSO_TAG.chain_out),
            ];
            // The factories retain different logical argument orders:
            // gamma(in, out, mu), sigma(mu, nu, in, out).
            let word = if head == AGS.gamma {
                FunctionBuilder::new(head)
                    .add_args(&spin_ports)
                    .add_args(&lorentz)
                    .finish()
            } else {
                FunctionBuilder::new(head)
                    .add_args(&lorentz)
                    .add_args(&spin_ports)
                    .finish()
            };
            let trace = FunctionBuilder::new(SPENSO_TAG.trace)
                .add_arg(spin.to_symbolic([]))
                .add_arg(&word)
                .finish();
            let value = SymbolicTensor::infer(trace.clone()).unwrap();
            assert_eq!(value.expression, trace);
            assert_eq!(value.structure.logical_slots().len(), lorentz.len());
            let chain = FunctionBuilder::new(SPENSO_TAG.chain)
                .add_arg(
                    spin.slot::<AbstractIndex, _>(AbstractIndex::Normal(99303))
                        .to_atom(),
                )
                .add_arg(
                    spin.slot::<AbstractIndex, _>(AbstractIndex::Normal(99304))
                        .to_atom(),
                )
                .add_arg(&word)
                .finish();
            let value = SymbolicTensor::infer(chain.clone()).unwrap();
            assert_eq!(value.expression, chain);
            assert_eq!(value.structure.logical_slots().len(), lorentz.len() + 2);
            let unequal = FunctionBuilder::new(SPENSO_TAG.chain)
                .add_arg(
                    spin.slot::<AbstractIndex, _>(AbstractIndex::Normal(99303))
                        .to_atom(),
                )
                .add_arg(
                    Bispinor {}
                        .new_rep(4)
                        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99304))
                        .to_atom(),
                )
                .add_arg(&word)
                .finish();
            assert!(SymbolicTensor::infer(unequal).is_err());
            for invalid in [
                ColorFundamental {}.new_rep(3).to_symbolic([]),
                FunctionBuilder::new(LibraryRep::from(Bispinor {}).symbol())
                    .add_arg(Atom::num(1) / Atom::num(2))
                    .finish(),
            ] {
                let trace = FunctionBuilder::new(SPENSO_TAG.trace)
                    .add_arg(invalid)
                    .add_arg(&word)
                    .finish();
                assert!(SymbolicTensor::infer(trace).is_err());
            }
        }
        let signature = AGS.gamma_strct::<AbstractIndex>(4);
        let arguments = signature
            .canonical()
            .external_reps_iter()
            .map(|mut port| {
                if port.rep == LibraryRep::from(Bispinor {}) {
                    port.dim = dimension;
                }
                port.to_symbolic([])
            })
            .collect::<Vec<_>>();
        let standalone = FunctionBuilder::new(AGS.gamma).add_args(arguments).finish();
        assert!(
            SymbolicTensor::infer(standalone).is_err(),
            "standalone factory signature remains unchanged"
        );
    }
}

#[test]
fn inferred_closed_chain_interface_excludes_its_contracted_endpoints() {
    crate::test_support::test_initialize();
    let spin = ColorFundamental {}.new_rep(3);
    let slot = LibraryRep::from(ColorAdjoint {})
        .new_rep(8)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99401));
    let word = FunctionBuilder::new(CS.t)
        .add_arg(slot.to_atom())
        .add_arg(Atom::var(SPENSO_TAG.chain_in))
        .add_arg(Atom::var(SPENSO_TAG.chain_out))
        .finish();
    let expression = FunctionBuilder::new(SPENSO_TAG.chain)
        .add_arg(
            spin.slot::<AbstractIndex, _>(AbstractIndex::Normal(99402))
                .to_atom(),
        )
        .add_arg(
            spin.dual()
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(99402))
                .to_atom(),
        )
        .add_arg(word)
        .finish();
    let value = SymbolicTensor::infer(expression.clone()).unwrap();
    assert_eq!(value.expression, expression);
    assert_eq!(
        value.structure.logical_slots(),
        vec![slot.rep().slot(PartialIndex::Explicit(slot.aind()))]
    );
}

#[test]
fn builtin_compact_ports_are_observed_without_constructing_components() {
    use std::sync::{Arc, Mutex};

    crate::test_support::test_initialize();
    let calls = Arc::new(Mutex::new(Vec::new()));
    let recorded = Arc::clone(&calls);
    let vector = spenso::vector_symbol!(
        "builtin_observed_component",
        norm = move |node, out| {
            recorded.lock().unwrap().push(node.to_owned());
            if let AtomView::Fun(function) = node
                && let Some(AtomView::Fun(slot)) = function.iter().last()
                && slot.get_nargs() == 2
            {
                **out = Atom::one();
            }
        }
    );
    let rep = LibraryRep::from(Minkowski {}).new_rep(4);
    let slot = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99701));
    let compact = FunctionBuilder::new(vector)
        .add_arg(rep.to_symbolic([]))
        .finish();
    let expression = shadowing::trace(
        Bispinor {}.new_rep(4).to_symbolic([]),
        [crate::gamma!(&compact), crate::gamma!(slot.to_atom())],
    );
    calls.lock().unwrap().clear();
    let observed = InterfaceInference::default()
        .infer_validated(expression.as_view())
        .unwrap();
    assert!(calls.lock().unwrap().is_empty());
    let value = SymbolicTensor::infer(expression.clone()).unwrap();
    assert_eq!(value.structure, observed);
    // Power lowering retains its established one compact normalizer call;
    // signature observation never constructs an explicit vector component.
    assert_eq!(*calls.lock().unwrap(), vec![compact]);
    assert_eq!(value.expression, expression);
    assert_eq!(
        value.structure.logical_slots(),
        vec![slot.rep().slot(PartialIndex::Explicit(slot.aind()))]
    );
    assert!(
        !InterfaceInference::default()
            .terminal_trace_preserves_interface(&expression, &value.structure)
    );
    assert_eq!(calls.lock().unwrap().len(), 1);

    // Observing the consumed port still checks earlier scalar metadata and
    // rejects a vector from the wrong representation.
    let extra = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99702));
    let tensor = FunctionBuilder::new(spenso::tensor_symbol!("builtin_compact_metadata"))
        .add_arg(extra.to_atom())
        .finish();
    let invalid_metadata = FunctionBuilder::new(vector)
        .add_arg(tensor)
        .add_arg(rep.to_symbolic([]))
        .finish();
    let wrong_rep = FunctionBuilder::new(vector)
        .add_arg(ExtendibleReps::EUCLIDEAN.new_rep(4).to_symbolic([]))
        .finish();
    for invalid in [invalid_metadata, wrong_rep] {
        let expression = shadowing::trace(
            Bispinor {}.new_rep(4).to_symbolic([]),
            [crate::gamma!(&invalid), crate::gamma!(slot.to_atom())],
        );
        calls.lock().unwrap().clear();
        assert!(
            InterfaceInference::default()
                .infer_validated(expression.as_view())
                .is_err()
        );
        assert!(calls.lock().unwrap().is_empty());
    }
}

#[test]
fn typed_tensor_ports_require_tagged_heads() {
    use spenso::network::{library::DummyLibrary, parsing::ParseSettings};
    use symbolica::atom::AtomOrView;

    let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
    let slot = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99701));
    let untagged = symbolica::symbol!("strict_admission_untagged_vector");
    let tagged = spenso::vector_symbol!("strict_admission_tagged_vector");
    for port in [slot.to_atom(), rep.to_symbolic([])] {
        let raw = FunctionBuilder::new(untagged)
            .add_arg(Atom::num(2))
            .add_arg(&port)
            .finish();
        let error = SymbolicTensor::<PartialStructure>::infer(raw.clone()).unwrap_err();
        assert!(matches!(error, TensorInferenceError::Invalid(ref message)
            if message.contains("not tagged as a tensor")));
        assert!(InterfaceInference::replacement_interface(raw.as_view()).is_err());

        let registered = FunctionBuilder::new(tagged)
            .add_arg(Atom::num(2))
            .add_arg(&port)
            .finish();
        let tensor = SymbolicTensor::<PartialStructure>::infer(registered.clone()).unwrap();
        assert_eq!(tensor.expression(), &registered);
        assert_eq!(tensor.structure().canonical().order(), 1);

        // The explicit raw-network opt-in retains its broader grammar.
        let settings = ParseSettings {
            strict_tensor_filter: StrictTensorFilter::ContainsReps,
            ..ParseSettings::default()
        };
        let network = crate::tensor::SymbolicNet::<AbstractIndex, AtomOrView<'_>>::try_from_view::<
            OrderedStructure,
            _,
        >(
            raw.as_view(),
            &DummyLibrary::<SymbolicTensor<OrderedStructure, AtomOrView<'_>>>::new(),
            &settings,
        );
        assert!(network.is_ok());
    }
}

#[test]
fn scalar_function_metadata_remains_opaque_at_typed_admission() {
    let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
    let slot = rep
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99702))
        .to_atom();
    let vector = FunctionBuilder::new(spenso::vector_symbol!("strict_admission_opaque_vector"))
        .add_arg(&slot)
        .finish();
    let opaque = symbolica::symbol!("strict_admission_opaque_function");
    for argument in [Atom::num(2), vector] {
        let source = FunctionBuilder::new(opaque).add_arg(argument).finish();
        let tensor = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
        assert!(tensor.is_scalar());
        assert_eq!(tensor.expression(), &source);
    }

    let untagged = FunctionBuilder::new(symbolica::symbol!("strict_admission_scalar_metadata"))
        .add_arg(&slot)
        .finish();
    for head in [
        symbolica::symbol!("strict_admission_declared_scalar"; Scalar),
        SPENSO_TAG.scalar,
        SPENSO_TAG.pure_scalar,
    ] {
        for argument in [&slot, &untagged] {
            let source = FunctionBuilder::new(head).add_arg(argument).finish();
            let tensor = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
            assert!(tensor.is_scalar());
            assert_eq!(tensor.expression(), &source);
        }
    }
}

#[test]
fn compact_inner_products_observe_ports_without_component_callbacks() {
    use std::sync::{Arc, Mutex};

    crate::test_support::test_initialize();
    let calls = Arc::new(Mutex::new(Vec::new()));
    let recorded = Arc::clone(&calls);
    let vector = spenso::vector_symbol!(
        "inner_product_observed_component",
        norm = move |node, out| {
            recorded.lock().unwrap().push(node.to_owned());
            if let AtomView::Fun(function) = node
                && let Some(AtomView::Fun(slot)) = function.iter().last()
                && slot.get_nargs() == 2
            {
                **out = Atom::one();
            }
        }
    );
    let rep = LibraryRep::from(Minkowski {}).new_rep(4);
    let compact = FunctionBuilder::new(vector)
        .add_arg(rep.to_symbolic([]))
        .finish();
    let expression = spenso::g!(&compact, &compact);
    calls.lock().unwrap().clear();
    let mut inference = InterfaceInference::default();
    let interface = inference.infer_validated(expression.as_view()).unwrap();
    assert!(interface.logical_slots().is_empty());
    assert!(inference.reusable_interfaces.is_empty());
    assert!(calls.lock().unwrap().is_empty());
    SymbolicTensor::<PartialStructure>::validate_interface(&expression, &interface).unwrap();
    assert!(calls.lock().unwrap().is_empty());

    let value = SymbolicTensor::infer(expression.clone()).unwrap();
    assert_eq!(value.expression, expression);
    assert!(value.structure.logical_slots().is_empty());
    // Public power lowering may reconstruct a compact input as before. It must
    // not manufacture the explicit component that would normalize to scalar 1.
    assert!(calls.lock().unwrap().iter().all(|call| call == &compact));

    let slot = rep
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99731))
        .to_atom();
    let tensor = FunctionBuilder::new(spenso::tensor_symbol!("inner_product_metadata"))
        .add_arg(&slot)
        .finish();
    let invalid_metadata = FunctionBuilder::new(vector)
        .add_arg(tensor)
        .add_arg(rep.to_symbolic([]))
        .finish();
    let wrong_rep = FunctionBuilder::new(vector)
        .add_arg(ExtendibleReps::EUCLIDEAN.new_rep(4).to_symbolic([]))
        .finish();
    let rank_two = FunctionBuilder::new(spenso::tensor_symbol!("inner_product_rank_two"))
        .add_arg(rep.to_symbolic([]))
        .add_arg(rep.to_symbolic([]))
        .finish();
    for invalid in [invalid_metadata, wrong_rep, rank_two] {
        let expression = spenso::g!(&compact, invalid);
        assert!(
            InterfaceInference::default()
                .infer_validated(expression.as_view())
                .is_err()
        );
    }
}

#[test]
fn scalar_finishers_share_strict_admission_and_preserve_identity_and_zero() {
    let scalar = SymbolicTensor::<PartialStructure>::infer(Atom::one()).unwrap();
    let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
    let slot = rep
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99703))
        .to_atom();
    let invalid = FunctionBuilder::new(symbolica::symbol!("strict_finisher_untagged"))
        .add_arg(&slot)
        .finish();
    assert!(scalar.with_checked_expression(invalid.clone()).is_err());
    assert!(scalar.with_transformed_expression(invalid.clone()).is_err());
    assert!(SymbolicTensor::validate_interface(&invalid, scalar.structure()).is_err());

    let vector = FunctionBuilder::new(spenso::vector_symbol!("strict_finisher_vector"))
        .add_arg(slot)
        .finish();
    let opaque = FunctionBuilder::new(symbolica::symbol!("strict_finisher_opaque"))
        .add_arg(&vector)
        .finish();
    assert!(
        scalar
            .with_checked_expression(opaque.clone())
            .unwrap()
            .is_scalar()
    );
    assert!(
        scalar
            .with_transformed_expression(opaque)
            .unwrap()
            .0
            .is_scalar()
    );

    let tensor = SymbolicTensor::<PartialStructure>::infer(vector).unwrap();
    INFERENCE_CALLS.with(|count| count.set(0));
    for value in [&scalar, &tensor] {
        assert_eq!(
            value
                .with_checked_expression(value.expression().clone())
                .unwrap(),
            *value
        );
        assert_eq!(
            value
                .with_transformed_expression(value.expression().clone())
                .unwrap(),
            (value.clone(), false)
        );
        let zero = value.with_checked_expression(Atom::Zero).unwrap();
        assert_eq!(zero.expression(), &Atom::Zero);
        assert_eq!(zero.structure(), value.structure());
        let (zero, contracted) = value.with_transformed_expression(Atom::Zero).unwrap();
        assert_eq!(zero.expression(), &Atom::Zero);
        assert_eq!(zero.structure(), value.structure());
        assert!(!contracted);
    }
    INFERENCE_CALLS.with(|count| assert_eq!(count.get(), 0));
}

#[test]
fn checked_error_precedence_keeps_strict_structured_and_scalar_failures() {
    crate::test_support::test_initialize();
    let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
    let slot = rep
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99803))
        .to_atom();
    let head = spenso::tensor_symbol!("strict_precedence_tensor");
    let tensor = SymbolicTensor::infer(FunctionBuilder::new(head).add_arg(&slot).finish()).unwrap();
    let scalar = SymbolicTensor::infer(Atom::one()).unwrap();
    let malformed_gamma = FunctionBuilder::new(AGS.gamma).add_arg(&slot).finish();
    for source in [&tensor, &scalar] {
        assert!(matches!(
            source.with_checked_expression(malformed_gamma.clone()),
            Err(TensorInferenceError::InvalidBuiltinSignature {
                factory: "gamma",
                ..
            })
        ));
    }
    for (invalid, message) in [
        (Atom::var(SPENSO_TAG.chain_in), "chain placeholders"),
        (
            FunctionBuilder::new(symbolica::symbol!("strict_precedence_untagged"))
                .add_arg(&slot)
                .finish(),
            "not tagged as a tensor",
        ),
        (Atom::var(head), "must be called"),
    ] {
        for error in [
            scalar.with_checked_expression(invalid.clone()).unwrap_err(),
            scalar
                .with_transformed_expression(invalid.clone())
                .unwrap_err(),
        ] {
            let TensorInferenceError::Invalid(reason) = error else {
                panic!("unexpected strict admission error: {error}");
            };
            assert!(reason.contains(message), "{reason}");
        }
    }
}

#[test]
fn custom_dual_metrics_keep_observed_port_order() {
    crate::test_support::test_initialize();
    let representations = [
        ExtendibleReps::new_dual("builtin_metric_order::R").unwrap(),
        ExtendibleReps::new_dual("builtin_metric_order::S").unwrap(),
    ];
    let vector = spenso::vector_symbol!("builtin_metric_order::q");
    for representation in representations {
        for representation in [representation, representation.dual()] {
            let rep = representation.new_rep(4);
            let i = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99901));
            let j = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99902));
            let dual_j = rep
                .dual()
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(99902));
            let i_atom = i.to_atom();
            let dual_j_atom = dual_j.to_atom();
            let target = FunctionBuilder::new(vector).add_arg(&i_atom).finish();
            for arguments in [[&i_atom, &dual_j_atom], [&dual_j_atom, &i_atom]] {
                // Check the original orders directly as well as the normalized
                // public input: metric symmetry may make those Atoms identical.
                let arguments = arguments.map(Atom::as_view);
                let mut inference = InterfaceInference::default();
                let signature = inference
                    .builtin_tensor_structure(ETS.metric, &arguments, false)
                    .unwrap()
                    .unwrap();
                assert_eq!(
                    signature.layout().canonical_to_logical(
                        &signature
                            .canonical()
                            .external_reps_iter()
                            .collect::<Vec<_>>()
                    ),
                    arguments.map(|argument| {
                        inference
                            .slots
                            .parse::<LibraryRep, AbstractIndex>(argument)
                            .unwrap()
                            .rep()
                    })
                );
                let metric = FunctionBuilder::new(ETS.metric)
                    .add_arg(arguments[0])
                    .add_arg(arguments[1])
                    .finish();
                let source = metric * FunctionBuilder::new(vector).add_arg(j.to_atom()).finish();
                let tensor = SymbolicTensor::infer(source).unwrap();
                assert_eq!(
                    tensor.structure().logical_slots(),
                    vec![rep.slot(PartialIndex::Explicit(i.aind()))]
                );
                // The ordered fallback finishes this dual space. The
                // coefficient collector deliberately refuses its unsupported
                // leaf, and the general rewrite proof requires self-dual ports;
                // neither refusal implies remaining contraction work.
                let target_tensor = SymbolicTensor::infer(target.clone()).unwrap();
                assert!(target_tensor.established_interface_is_valid());
                assert!(!target_tensor.rewrites_preserve_interface());
                assert!(
                    crate::shorthands::schoonschip::SlotContraction::new()
                        .contract_factorized(target.as_view(), None, true)
                        .is_none()
                );
                let contracted = tensor.contract(Default::default()).unwrap();
                assert!(contracted.contraction_complete());
                let result = contracted.resolved().unwrap();
                assert_eq!(result.expression(), &target);
                assert_eq!(result.structure(), tensor.structure());
                assert!(
                    result
                        .contract(Default::default())
                        .unwrap()
                        .contraction_complete()
                );
            }
        }
    }
}

#[test]
fn custom_dual_metrics_still_reject_incompatible_ports() {
    crate::test_support::test_initialize();
    let first = ExtendibleReps::new_dual("builtin_metric_invalid::R").unwrap();
    let other = ExtendibleReps::new_dual("builtin_metric_invalid::S").unwrap();
    for rep in [first, first.dual()] {
        let left = rep
            .new_rep(4)
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(99911))
            .to_atom();
        for wrong in [
            rep.dual().new_rep(5),
            other.new_rep(4),
            other.dual().new_rep(4),
        ] {
            let right = wrong
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(99912))
                .to_atom();
            for arguments in [[&left, &right], [&right, &left]] {
                let expression = FunctionBuilder::new(ETS.metric)
                    .add_arg(arguments[0])
                    .add_arg(arguments[1])
                    .finish();
                assert!(matches!(
                    SymbolicTensor::infer(expression),
                    Err(TensorInferenceError::InvalidBuiltinSignature { factory: "g", .. })
                ));
            }
        }
    }
}

#[test]
fn compact_dots_keep_explicit_spectators_and_reject_ambiguous_axes() {
    crate::test_support::test_initialize();
    let rep = LibraryRep::from(Minkowski {}).new_rep(4);
    let spin = LibraryRep::from(Bispinor {}).new_rep(4);
    let slots = [99801, 99802, 99803, 99804]
        .map(|index| spin.slot::<AbstractIndex, _>(AbstractIndex::Normal(index)));
    let head = spenso::tensor_symbol!("compact_dot_spectators");
    let operand = |ports: &[Atom], compact: Atom| {
        FunctionBuilder::new(head)
            .add_args(ports)
            .add_arg(compact)
            .finish()
    };
    for dot_head in [SPENSO_TAG.dot, ETS.metric] {
        for right_start in [1, 2] {
            let left = operand(
                &[slots[0].to_atom(), slots[1].to_atom()],
                rep.to_symbolic([]),
            );
            let right = operand(
                &[slots[right_start].to_atom(), slots[3].to_atom()],
                rep.to_symbolic([]),
            );
            let expression = FunctionBuilder::new(dot_head)
                .add_arg(&left)
                .add_arg(&right)
                .finish();
            let value = SymbolicTensor::infer(expression.clone()).unwrap();
            let expected = if right_start == 1 {
                vec![slots[0], slots[3]]
            } else {
                slots.to_vec()
            };
            let expected = PartialStructure::from_logical_slots(
                expected
                    .into_iter()
                    .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
            );
            assert_eq!(value.expression, expression);
            assert_eq!(value.structure, expected);
            assert_eq!(
                value.with_checked_expression(expression.clone()).unwrap(),
                value
            );
            let opened = value.undo_dots().unwrap();
            assert_eq!(opened.structure, expected);
            assert_eq!(
                SymbolicTensor::infer(opened.expression.clone())
                    .unwrap()
                    .structure,
                expected
            );
        }
        let compact = operand(
            &[slots[0].to_atom(), slots[1].to_atom()],
            rep.to_symbolic([]),
        );
        let multiple = operand(&[rep.to_symbolic([])], rep.to_symbolic([]));
        let wrong = operand(
            &[slots[2].to_atom(), slots[3].to_atom()],
            ExtendibleReps::EUCLIDEAN.new_rep(4).to_symbolic([]),
        );
        for invalid in [multiple, wrong] {
            assert!(
                SymbolicTensor::infer(
                    FunctionBuilder::new(dot_head)
                        .add_arg(&compact)
                        .add_arg(invalid)
                        .finish()
                )
                .is_err()
            );
        }
    }
    let dimension = symbolica::symbol!("compact_dot_gamma_D");
    let lorentz = LibraryRep::from(Minkowski {}).new_rep(Dimension::from(dimension));
    let gamma = |left: Atom, right: Atom| {
        FunctionBuilder::new(AGS.gamma)
            .add_arg(left)
            .add_arg(right)
            .add_arg(lorentz.to_symbolic([]))
            .finish()
    };
    let expression = FunctionBuilder::new(SPENSO_TAG.dot)
        .add_arg(gamma(slots[0].to_atom(), slots[1].to_atom()))
        .add_arg(gamma(slots[1].to_atom(), slots[3].to_atom()))
        .finish();
    let value = SymbolicTensor::infer(expression.clone()).unwrap();
    assert_eq!(value.expression, expression);
    assert_eq!(
        value.structure.logical_slots(),
        vec![
            spin.slot(PartialIndex::Explicit(slots[0].aind())),
            spin.slot(PartialIndex::Explicit(slots[3].aind())),
        ]
    );
}

#[test]
fn compact_dot_tensor_spectators_preserve_materialization_and_observation_callbacks() {
    use std::sync::{Arc, Mutex};
    crate::test_support::test_initialize();
    let calls = Arc::new(Mutex::new(Vec::new()));
    let observed = Arc::clone(&calls);
    let head = spenso::tensor_symbol!(
        "compact_dot_tensor_callback",
        norm = move |node, _out| {
            observed.lock().unwrap().push(node.to_owned());
        }
    );
    let rep = LibraryRep::from(Minkowski {}).new_rep(4);
    let spectator = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99811));
    let compact = FunctionBuilder::new(head)
        .add_arg(spectator.to_atom())
        .add_arg(rep.to_symbolic([]))
        .finish();
    let other = FunctionBuilder::new(spenso::vector_symbol!("compact_dot_plain_vector"))
        .add_arg(rep.to_symbolic([]))
        .finish();
    let dot = FunctionBuilder::new(SPENSO_TAG.dot)
        .add_arg(&compact)
        .add_arg(other)
        .finish();
    calls.lock().unwrap().clear();
    let mut observer = InterfaceInference {
        leaf_inference: LeafInference::Observe,
        ..Default::default()
    };
    let interface = observer.infer_validated(dot.as_view()).unwrap();
    assert_eq!(
        interface.logical_slots(),
        vec![rep.slot(PartialIndex::Explicit(spectator.aind()))]
    );
    assert_eq!(observer.infer_validated(dot.as_view()).unwrap(), interface);
    assert!(calls.lock().unwrap().is_empty());
    // Ordinary tensor operands retain their existing constructor schedule;
    // only tagged rank-one consumed vectors use observation-only admission.
    let mut constructor = InterfaceInference::default();
    assert_eq!(
        constructor.infer_validated(dot.as_view()).unwrap(),
        interface
    );
    let first = calls.lock().unwrap().len();
    assert!(first > 0);
    assert_eq!(
        constructor.infer_validated(dot.as_view()).unwrap(),
        interface
    );
    assert!(calls.lock().unwrap().len() > first);
}
