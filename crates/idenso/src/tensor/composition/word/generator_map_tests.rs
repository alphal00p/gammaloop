use super::*;
use crate::{bis, dirac::AGS};
use spenso::structure::representation::Minkowski;
use spenso::structure::slot::DummyAind;
use symbolica::{function, symbol};

type Tensor = SymbolicTensor<PartialStructure>;

fn head() -> Symbol {
    crate::test_support::test_initialize();
    symbol!("word_generator_map::exterior"; Flat, Linear, Antisymmetric;
        tag = &SPENSO_TAG.tensor_word)
}

fn slot(index: usize) -> Atom {
    Minkowski {}
        .new_rep(4)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
        .to_atom()
}

fn leaf(name: &str, slots: &[Atom]) -> Tensor {
    Tensor::infer(
        FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
            .add_args(slots)
            .finish(),
    )
    .unwrap()
}

/// This fixture's explicit family map preserves the original native arguments.
fn renamed(generator: &Tensor, name: &str) -> Tensor {
    let AtomView::Fun(function) = generator.expression.as_view() else {
        panic!("fixture generator must be a tensor function")
    };
    Tensor::from_declared_parts(
        FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
            .add_args(function.iter())
            .finish(),
        generator.structure.clone(),
    )
    .unwrap()
}

fn marked(head: Symbol, value: &Tensor) -> Tensor {
    Tensor::word(head, std::slice::from_ref(value)).unwrap()
}

fn canonical(value: &Tensor) -> Atom {
    value.canonize().unwrap().expression
}

#[test]
fn generator_map_preserves_identity_and_coefficients_without_callbacks() {
    let head = head();
    let a = leaf("word_generator_map::identity_a", &[slot(11)]);
    let b = leaf("word_generator_map::identity_b", &[slot(11)]);
    let coefficient = Tensor::infer(Atom::num((3, 7)) + Atom::i()).unwrap();
    let source = coefficient
        .multiply_word(head, &Tensor::word(head, &[a, b]).unwrap())
        .unwrap()
        .wrap_indices(symbol!("word_generator_map::original_scope"), true)
        .unwrap();
    let mut calls = 0;
    let actual = source
        .map_word_generators(head, |generator| {
            calls += 1;
            Ok(marked(head, generator))
        })
        .unwrap();
    assert_eq!(calls, 2);
    assert!(actual.same_value(&source));
    let actual = coefficient
        .map_word_generators(head, |_| {
            panic!("a word-free coefficient has no generators")
        })
        .unwrap();
    assert!(actual.same_value(&coefficient));
}

#[test]
fn generator_map_extracts_rich_dirac_tensor_coefficients() {
    let head = head();
    let psi = leaf("word_generator_map::psi", &[bis!(4, 21)]);
    let chi = leaf("word_generator_map::chi", &[bis!(4, 22)]);
    let current = leaf("word_generator_map::current", &[slot(23)]);
    let gamma = Tensor::infer(function!(AGS.gamma, bis!(4, 21), bis!(4, 22), slot(23))).unwrap();
    let rhs = current
        .multiply(&gamma)
        .unwrap()
        .multiply_word(head, &marked(head, &chi))
        .unwrap();
    assert_eq!(rhs.structure.logical_slots(), psi.structure.logical_slots());
    let source = marked(head, &psi);
    let actual = source
        .map_word_generators(head, |_| Ok(rhs.clone()))
        .unwrap();
    let expected = rhs
        .wrap_indices(symbol!("word_generator_map::expected_dirac"), true)
        .unwrap();
    assert_eq!(
        actual.structure.logical_slots(),
        source.structure.logical_slots()
    );
    assert_eq!(canonical(&actual), canonical(&expected));
}

#[test]
fn generator_map_reserves_source_and_all_raw_values_before_freshening() {
    let head = head();
    let a = leaf("word_generator_map::fresh_a", &[slot(31)]);
    let b = leaf("word_generator_map::fresh_b", &[slot(31)]);
    // The original later generator uses the written internal RHS label.
    let spectator = leaf("word_generator_map::fresh_spectator", &[slot(99)]);
    let source = Tensor::word(head, &[a.clone(), b.clone(), spectator.clone()]).unwrap();
    let pair = leaf("word_generator_map::J", &[slot(99)])
        .multiply(&leaf("word_generator_map::K", &[slot(99)]))
        .unwrap();
    let left = pair
        .multiply_word(
            head,
            &marked(head, &renamed(&a, "word_generator_map::fresh_chi")),
        )
        .unwrap();
    let right = pair
        .multiply_word(
            head,
            &marked(head, &renamed(&b, "word_generator_map::fresh_eta")),
        )
        .unwrap();
    let actual = source
        .map_word_generators(head, |generator| {
            Ok(
                if generator.expression.get_symbol() == a.expression.get_symbol() {
                    left.clone()
                } else if generator.expression.get_symbol() == b.expression.get_symbol() {
                    right.clone()
                } else {
                    marked(head, generator)
                },
            )
        })
        .unwrap();
    let expected = left
        .wrap_indices(symbol!("word_generator_map::expected_left"), true)
        .unwrap()
        .multiply_word(
            head,
            &right
                .wrap_indices(symbol!("word_generator_map::expected_right"), true)
                .unwrap(),
        )
        .unwrap()
        .multiply_word(head, &marked(head, &spectator))
        .unwrap();
    assert_eq!(
        actual.structure.logical_slots(),
        source.structure.logical_slots()
    );
    assert!(!actual.is_zero());
    assert_eq!(canonical(&actual), canonical(&expected));
}

#[test]
fn generator_map_is_simultaneous_and_does_not_revisit_inserted_generators() {
    let head = head();
    let a = leaf("word_generator_map::sim_a", &[]);
    let b = leaf("word_generator_map::sim_b", &[]);
    let source = Tensor::word(head, &[a.clone(), b.clone()]).unwrap();
    let mut calls = 0;
    let actual = source
        .map_word_generators(head, |generator| {
            calls += 1;
            Ok(marked(head, if generator.same_value(&a) { &b } else { &a }))
        })
        .unwrap();
    assert_eq!(calls, 2);
    assert_eq!(actual.expression, -source.expression.clone());

    let c = leaf("word_generator_map::sim_c", &[]);
    let triple = Tensor::word(head, &[a.clone(), b, c]).unwrap();
    calls = 0;
    let actual = marked(head, &a)
        .map_word_generators(head, |_| {
            calls += 1;
            Ok(triple.clone())
        })
        .unwrap();
    assert_eq!(calls, 1);
    assert!(actual.same_value(&triple));
}

#[test]
fn generator_map_distributes_sums_and_repeated_closed_word_powers() {
    let head = head();
    let a = leaf("word_generator_map::sum_a", &[]);
    let b = leaf("word_generator_map::sum_b", &[]);
    let c = leaf("word_generator_map::sum_c", &[]);
    let d = leaf("word_generator_map::sum_d", &[]);
    let source = Tensor::word(head, &[a.clone(), b.clone()]).unwrap();
    let rhs = marked(head, &c).try_add(&marked(head, &d)).unwrap();
    let actual = source
        .map_word_generators(head, |generator| {
            Ok(if generator.same_value(&a) {
                rhs.clone()
            } else {
                marked(head, generator)
            })
        })
        .unwrap();
    let expected = Tensor::word(head, &[c, b.clone()])
        .unwrap()
        .try_add(&Tensor::word(head, &[d, b]).unwrap())
        .unwrap();
    assert!(actual.same_value(&expected));

    let a = leaf("word_generator_map::power_a", &[slot(41)]);
    let b = leaf("word_generator_map::power_b", &[slot(41)]);
    let source = Tensor::word(head, &[a.clone(), b.clone()])
        .unwrap()
        .pow_word(head, 2)
        .unwrap();
    let mut calls = 0;
    let actual = source
        .map_word_generators(head, |generator| {
            calls += 1;
            let name = if generator.expression.get_symbol() == a.expression.get_symbol() {
                "word_generator_map::power_c"
            } else {
                "word_generator_map::power_d"
            };
            Ok(marked(head, &renamed(generator, name)))
        })
        .unwrap();
    let expected = Tensor::word(
        head,
        &[
            renamed(&a, "word_generator_map::power_c"),
            renamed(&b, "word_generator_map::power_d"),
        ],
    )
    .unwrap()
    .pow_word(head, 2)
    .unwrap();
    assert_eq!(calls, 4);
    assert_eq!(
        actual.structure.logical_slots(),
        source.structure.logical_slots()
    );
    assert!(!actual.is_zero());
    assert_eq!(canonical(&actual), canonical(&expected));
}

#[test]
fn generator_map_preserves_typed_zero_and_refuses_unresolved_nonzero() {
    let head = head();
    let a = leaf("word_generator_map::zero_a", &[slot(51)]);
    let source = marked(head, &a);
    let zero = Tensor::from_declared_parts(Atom::Zero, a.structure.clone()).unwrap();
    let actual = source
        .map_word_generators(head, |_| Ok(zero.clone()))
        .unwrap();
    assert!(actual.is_zero());
    assert_eq!(
        actual.structure.logical_slots(),
        source.structure.logical_slots()
    );

    let representation = Minkowski {}.new_rep(4).cast::<LibraryRep>();
    let unresolved = Tensor::from_declared_parts(
        function!(
            SPENSO_TAG.tensor_symbol("word_generator_map::unresolved"),
            representation.to_symbolic([])
        ),
        PartialStructure::from_logical_slots([representation.slot(PartialIndex::open(0))]),
    )
    .unwrap();
    let zero = Tensor::from_declared_parts(Atom::Zero, unresolved.structure.clone()).unwrap();
    let actual = zero
        .map_word_generators(head, |_| panic!("zero has no generators"))
        .unwrap();
    assert!(actual.same_value(&zero));
    let source = unresolved
        .multiply_word(
            head,
            &marked(head, &leaf("word_generator_map::open_odd", &[])),
        )
        .unwrap();
    assert!(
        source
            .map_word_generators(head, |_| panic!("refuse before callbacks"))
            .is_err()
    );
}

#[test]
fn generator_map_checks_complete_slot_identity_variance_and_odd_markers() {
    let head = head();
    let rep = LibraryRep::new_dual("word_generator_map::flavor")
        .unwrap()
        .new_rep(3);
    let named = |representation: Representation<LibraryRep>, index| {
        representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
            .to_atom()
    };
    let a = leaf("word_generator_map::charged", &[named(rep, 61)]);
    let source = marked(head, &a);
    for wrong in [
        marked(
            head,
            &leaf("word_generator_map::wrong_label", &[named(rep, 62)]),
        ),
        marked(
            head,
            &leaf("word_generator_map::wrong_dual", &[named(rep.dual(), 61)]),
        ),
        a.clone(), // A raw tensor is an even coefficient, never an implicit odd marker.
    ] {
        assert!(
            source
                .map_word_generators(head, |_| Ok(wrong.clone()))
                .is_err()
        );
    }
    let scalar = leaf("word_generator_map::odd_scalar", &[]);
    let even = Tensor::word(head, &[a.clone(), scalar]).unwrap();
    assert!(
        source
            .map_word_generators(head, |_| Ok(even.clone()))
            .is_err()
    );
    let even_term = a.clone();
    let mixed = marked(head, &a).try_add(&even_term).unwrap();
    assert!(
        source
            .map_word_generators(head, |_| Ok(mixed.clone()))
            .is_err()
    );
}

#[test]
fn generator_map_refuses_foreign_hidden_and_powered_word_grammars() {
    let head = head();
    let a = leaf("word_generator_map::grammar_a", &[]);
    let b = leaf("word_generator_map::grammar_b", &[]);
    let foreign = symbol!("word_generator_map::foreign"; Flat, Linear, Antisymmetric;
        tag = &SPENSO_TAG.tensor_word);
    let source = marked(head, &a);
    let foreign = Tensor::word(foreign, std::slice::from_ref(&b)).unwrap();
    assert!(
        source
            .map_word_generators(head, |_| Ok(foreign.clone()))
            .is_err()
    );
    for expression in [
        function!(
            symbol!("word_generator_map::opaque"; Scalar),
            source.expression()
        ),
        source.expression.clone().pow(Atom::num(-1)),
    ] {
        let hidden = Tensor::from_declared_parts(expression, source.structure.clone()).unwrap();
        assert!(
            hidden
                .map_word_generators(head, |_| panic!("invalid source grammar"))
                .is_err()
        );
        assert!(
            source
                .map_word_generators(head, |_| Ok(hidden.clone()))
                .is_err()
        );
    }
    let ordered = symbol!("word_generator_map::ordered"; Flat, Linear;
        tag = &SPENSO_TAG.tensor_word);
    assert!(
        Tensor::word(ordered, &[a])
            .unwrap()
            .map_word_generators(ordered, |_| { panic!("non-exterior head") })
            .is_err()
    );
}

#[test]
fn generator_map_preserves_opaque_metadata_and_unchanged_normalizers() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    let head = head();
    let metadata = symbol!("word_generator_map::metadata"; Scalar);
    let normalizer = symbol!(
        "word_generator_map::normalized_generator",
        tag = &SPENSO_TAG.tensor,
        norm = |_, _| {
            CALLS.fetch_add(1, Ordering::Relaxed);
        }
    );
    let index = slot(71);
    let opaque = function!(metadata, index.clone());
    let a = Tensor::infer(function!(normalizer, opaque.clone(), index)).unwrap();
    let source = marked(head, &a);
    let before = CALLS.load(Ordering::Relaxed);
    let actual = source
        .map_word_generators(head, |generator| Ok(marked(head, generator)))
        .unwrap();
    assert_eq!(CALLS.load(Ordering::Relaxed), before);
    assert!(actual.same_value(&source));
    assert!(actual.expression.contains(opaque.as_view()));

    let opaque_bound = function!(metadata, slot(72));
    let mapped = leaf(
        "word_generator_map::metadata_response",
        &[opaque_bound.clone(), slot(71)],
    );
    let pair = leaf("word_generator_map::metadata_J", &[slot(72)])
        .multiply(&leaf("word_generator_map::metadata_K", &[slot(72)]))
        .unwrap();
    let rhs = pair.multiply_word(head, &marked(head, &mapped)).unwrap();
    let actual = source
        .map_word_generators(head, |_| Ok(rhs.clone()))
        .unwrap();
    let expected = pair
        .wrap_indices(symbol!("word_generator_map::metadata_expected"), true)
        .unwrap()
        .multiply_word(head, &marked(head, &mapped))
        .unwrap();
    assert!(actual.expression.contains(opaque_bound.as_view()));
    assert_eq!(
        actual.structure.logical_slots(),
        source.structure.logical_slots()
    );
    assert_eq!(canonical(&actual), canonical(&expected));
}

#[test]
fn generator_map_reserves_existing_and_nested_scope_symbols() {
    let head = head();
    let a = leaf("word_generator_map::scoped_a", &[]);
    let current = leaf("word_generator_map::scoped_J", &[slot(81)])
        .multiply(&leaf("word_generator_map::scoped_K", &[slot(81)]))
        .unwrap();
    let dummy_scope = |index| AbstractIndex::Dummy(index).to_atom().get_symbol().unwrap();
    for ambient in [
        current.wrap_indices(dummy_scope(0), true).unwrap(),
        current
            .wrap_indices(dummy_scope(0), true)
            .unwrap()
            .wrap_indices(dummy_scope(1), true)
            .unwrap(),
    ] {
        let source = ambient.multiply_word(head, &marked(head, &a)).unwrap();
        let rhs = current
            .multiply_word(
                head,
                &marked(head, &leaf("word_generator_map::scoped_b", &[])),
            )
            .unwrap();
        let actual = source
            .map_word_generators(head, |_| Ok(rhs.clone()))
            .unwrap();
        let expected = ambient
            .multiply_word(
                head,
                &rhs.wrap_indices(symbol!("word_generator_map::fresh_expected_scope"), true)
                    .unwrap(),
            )
            .unwrap();
        assert!(!actual.is_zero());
        assert_eq!(canonical(&actual), canonical(&expected));
    }
}

#[test]
fn generator_map_refuses_rank_changes_from_scoping_normalizers() {
    let head = head();
    let normalizer = symbol!(
        "word_generator_map::scope_sensitive_coefficient",
        tag = &SPENSO_TAG.tensor,
        norm = |node, output| {
            if let AtomView::Fun(function) = node
                && function.iter().any(|argument| {
                    matches!(Slot::<LibraryRep, AbstractIndex>::try_from(argument),
                    Ok(slot) if matches!(slot.aind, AbstractIndex::Scoped(_)))
                })
            {
                **output = Atom::one();
            }
        }
    );
    let coefficient = Tensor::infer(function!(normalizer, slot(91)))
        .unwrap()
        .multiply(&leaf(
            "word_generator_map::scope_sensitive_partner",
            &[slot(91)],
        ))
        .unwrap();
    let source = marked(head, &leaf("word_generator_map::scope_sensitive_a", &[]));
    let rhs = coefficient
        .multiply_word(
            head,
            &marked(head, &leaf("word_generator_map::scope_sensitive_b", &[])),
        )
        .unwrap();
    // Scoping must re-enter the checked owner: the callback erases only one
    // member of a contracted pair, leaving an unplanned free coefficient port.
    assert!(
        source
            .map_word_generators(head, |_| Ok(rhs.clone()))
            .is_err()
    );
}

#[test]
fn generator_map_scope_collision_at_the_actual_allocator_frontier() {
    const CHILD: &str = "IDENSO_WORD_MAP_SCOPE_FRONTIER_CHILD";
    if std::env::var_os(CHILD).is_none() {
        let output = std::process::Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "tensor::composition::word::generator_map_tests::generator_map_scope_collision_at_the_actual_allocator_frontier",
                "--test-threads=1",
                "--nocapture",
            ])
            .env(CHILD, "1")
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(String::from_utf8_lossy(&output.stdout).contains("1 passed"));
        return;
    }
    let head = head();
    let pair = leaf("word_generator_map::frontier_J", &[slot(101)])
        .multiply(&leaf("word_generator_map::frontier_K", &[slot(101)]))
        .unwrap();
    let a = marked(head, &leaf("word_generator_map::frontier_a", &[]));
    let rhs = pair
        .multiply_word(
            head,
            &marked(head, &leaf("word_generator_map::frontier_b", &[])),
        )
        .unwrap();
    // Public native allocation, not a counter reset or a parsed dummy name.
    // Everything constructed after this point has explicit or no ports, so
    // no unresolved-port materialization consumes the next native dummy.
    let AbstractIndex::Dummy(previous) = AbstractIndex::new_dummy() else {
        unreachable!()
    };
    let occupied_scope = AbstractIndex::Dummy(previous.checked_add(1).unwrap())
        .to_atom()
        .get_symbol()
        .unwrap();
    let ambient = pair.wrap_indices(occupied_scope, true).unwrap();
    let source = ambient.multiply_word(head, &a).unwrap();
    let actual = source
        .map_word_generators(head, |_| Ok(rhs.clone()))
        .unwrap();
    let expected = ambient
        .multiply_word(
            head,
            &rhs.wrap_indices(symbol!("word_generator_map::frontier_independent"), true)
                .unwrap(),
        )
        .unwrap();
    assert_eq!(
        actual.structure.logical_slots(),
        source.structure.logical_slots()
    );
    assert_eq!(canonical(&actual), canonical(&expected));
}

#[test]
fn word_constructor_reuses_explicit_operands_but_materializes_unresolved_ones() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static CALLS: AtomicUsize = AtomicUsize::new(0);
    let head = head();
    let normalizer = symbol!(
        "word_generator_map::constructor_normalizer",
        tag = &SPENSO_TAG.tensor,
        norm = |_, _| {
            CALLS.fetch_add(1, Ordering::Relaxed);
        }
    );
    let explicit = Tensor::infer(function!(normalizer, slot(111))).unwrap();
    let before = CALLS.load(Ordering::Relaxed);
    let word = marked(head, &explicit);
    assert_eq!(CALLS.load(Ordering::Relaxed), before);
    assert_eq!(
        word.structure().logical_slots(),
        explicit.structure().logical_slots()
    );

    let representation = Minkowski {}.new_rep(4).cast::<LibraryRep>();
    let unresolved = Tensor::from_declared_parts(
        function!(normalizer, representation.to_symbolic([])),
        PartialStructure::from_logical_slots([representation.slot(PartialIndex::open(0))]),
    )
    .unwrap();
    let before = CALLS.load(Ordering::Relaxed);
    let materialized = marked(head, &unresolved);
    assert_eq!(CALLS.load(Ordering::Relaxed), before + 1);
    assert_eq!(materialized.rank(), 1);
    assert!(materialized.structure().open_positions().is_empty());

    let invalid = Tensor::new(
        function!(
            SPENSO_TAG.tensor_symbol("word_generator_map::unchecked"),
            SPENSO_TAG.chain_in
        ),
        PartialStructure::from_logical_slots([]),
    );
    assert!(Tensor::word(head, &[invalid]).is_err());
    let erased = symbol!(
        "word_generator_map::constructor_erased",
        tag = &SPENSO_TAG.tensor_word,
        norm = |_, output| {
            **output = Atom::one();
        }
    );
    assert!(Tensor::word(erased, std::slice::from_ref(&explicit)).is_err());
    let annihilated = symbol!(
        "word_generator_map::constructor_zero",
        tag = &SPENSO_TAG.tensor_word,
        norm = |_, output| {
            **output = Atom::Zero;
        }
    );
    let zero = Tensor::word(annihilated, &[explicit]).unwrap();
    assert!(zero.is_zero());
    assert_eq!(zero.rank(), 1);
}
