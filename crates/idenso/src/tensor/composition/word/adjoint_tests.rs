use super::*;
use crate::dirac::AGS;
use crate::tensor::AlgebraSettings;
use crate::test_support;
use crate::{bis, gamma};
use spenso::mink;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::PartialIndex;
use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
use spenso::structure::slot::{IsAbstractSlot, SlotMatcher};
use std::collections::HashMap;
use std::slice;
use symbolica::{function, symbol};

fn head() -> Symbol {
    test_support::test_initialize();
    symbol!("word_adjoint::exterior"; Flat, Linear, Antisymmetric;
        tag = &SPENSO_TAG.tensor_word)
}

fn leaf(name: &str, slots: &[Atom]) -> SymbolicTensor<PartialStructure> {
    SymbolicTensor::infer(
        FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
            .add_args(slots)
            .finish(),
    )
    .unwrap()
}

/// A test domain's explicit immutable field/bar pairing, including jet metadata.
/// This supplies no coefficient conjugation, word reversal, or sign rules.
fn bar(
    source: &SymbolicTensor<PartialStructure>,
) -> Result<SymbolicTensor<PartialStructure>, AdjointError> {
    let pairs = [
        ("word_adjoint::a", "word_adjoint::abar"),
        ("word_adjoint::b", "word_adjoint::bbar"),
        ("word_adjoint::jet", "word_adjoint::barjet"),
    ];
    let AtomView::Fun(function) = source.expression.as_view() else {
        return Err(TensorInferenceError::invalid("unsupported field policy").into());
    };
    let target = pairs.into_iter().find_map(|(left, right)| {
        let left = SPENSO_TAG.tensor_symbol(left);
        let right = SPENSO_TAG.tensor_symbol(right);
        if function.get_symbol() == left {
            Some(right)
        } else if function.get_symbol() == right {
            Some(left)
        } else {
            None
        }
    });
    let Some(target) = target else {
        return Err(TensorInferenceError::invalid("unsupported field policy").into());
    };
    let mut slots = SlotMatcher::default();
    let expression = FunctionBuilder::new(target)
        .add_args(function.iter().map(|argument| {
            slots
                .parse::<LibraryRep, AbstractIndex>(argument)
                .map_or_else(|_| argument.to_owned(), |slot| slot.dual().to_atom())
        }))
        .finish();
    Ok(SymbolicTensor::from_declared_parts(
        expression,
        PartialStructure::from_logical_slots(
            source
                .structure
                .logical_slots()
                .into_iter()
                .map(|slot| slot.dual()),
        ),
    )?)
}

#[test]
fn word_adjoint_conjugates_exact_scalar_weights_and_reverses_generators() {
    let head = head();
    let a = leaf("word_adjoint::a", &[mink!(4, 1)]);
    let b = leaf("word_adjoint::b", &[mink!(4, 1)]);
    let word = SymbolicTensor::word(head, &[a.clone(), b.clone()]).unwrap();
    let reversed = SymbolicTensor::word(head, &[bar(&b).unwrap(), bar(&a).unwrap()]).unwrap();
    let unreversed = SymbolicTensor::word(head, &[bar(&a).unwrap(), bar(&b).unwrap()]).unwrap();
    assert_eq!(reversed.expression, -unreversed.expression);
    for weight in [
        Atom::num((2, 3)) + Atom::i(),
        Atom::var(symbol!("word_adjoint::real"; Real, Scalar)),
        Atom::var(symbol!("word_adjoint::complex"; Scalar)),
    ] {
        let scalar = SymbolicTensor::infer(weight.clone()).unwrap();
        let source = scalar.multiply_word(head, &word).unwrap();
        let actual = source.dirac_adjoint_word(head, bar).unwrap();
        let expected = SymbolicTensor::infer(weight.conj())
            .unwrap()
            .multiply_word(head, &reversed)
            .unwrap();
        assert!(actual.same_value(&expected));
        assert!(
            actual
                .dirac_adjoint_word(head, bar)
                .unwrap()
                .same_value(&source)
        );
    }
}

#[test]
fn word_adjoint_preserves_free_slot_order_duality_labels_and_jet_metadata() {
    let head = head();
    let representation = LibraryRep::new_dual("word_adjoint::flavor")
        .unwrap()
        .new_rep(3);
    let a = leaf(
        "word_adjoint::a",
        &[representation
            .slot::<AbstractIndex, _>(AbstractIndex::Normal(31))
            .to_atom()],
    );
    let jet = leaf(
        "word_adjoint::jet",
        &[
            Atom::num(2),
            mink!(4, 32),
            mink!(4, 33),
            representation
                .dual()
                .slot::<AbstractIndex, _>(AbstractIndex::Normal(34))
                .to_atom(),
        ],
    );
    let source = SymbolicTensor::word(head, &[a, jet])
        .unwrap()
        .wrap_indices(symbol!("word_adjoint::owned_free"), false)
        .unwrap()
        .permuted(&[3, 1, 0, 2])
        .unwrap();
    let adjoint = source.dirac_adjoint_word(head, bar).unwrap();
    assert_eq!(
        adjoint.structure.logical_slots(),
        source
            .structure
            .logical_slots()
            .into_iter()
            .map(|slot| slot.dual())
            .collect::<Vec<_>>()
    );
    assert!(
        adjoint
            .expression
            .contains_symbol(SPENSO_TAG.tensor_symbol("word_adjoint::barjet"))
    );
    assert!(
        adjoint
            .expression
            .contains_symbol(symbol!("word_adjoint::owned_free"))
    );
    assert!(
        adjoint
            .dirac_adjoint_word(head, bar)
            .unwrap()
            .same_value(&source)
    );
    assert!(
        source
            .dirac_adjoint_word(head, |leaf| Ok(leaf.clone()))
            .is_err()
    );
    let closed = SymbolicTensor::word(
        head,
        &[
            leaf(
                "word_adjoint::a",
                &[representation
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(35))
                    .to_atom()],
            ),
            leaf(
                "word_adjoint::b",
                &[representation
                    .dual()
                    .slot::<AbstractIndex, _>(AbstractIndex::Normal(35))
                    .to_atom()],
            ),
        ],
    )
    .unwrap()
    .wrap_indices(symbol!("word_adjoint::owned_bound"), true)
    .unwrap();
    let adjoint = closed.dirac_adjoint_word(head, bar).unwrap();
    assert!(adjoint.is_scalar());
    assert!(
        adjoint
            .expression
            .contains_symbol(symbol!("word_adjoint::owned_bound"))
    );
    assert!(
        adjoint
            .dirac_adjoint_word(head, bar)
            .unwrap()
            .same_value(&closed)
    );
}

#[test]
fn word_adjoint_uses_native_gamma_and_projector_coefficient_adjoint() {
    let head = head();
    let canonical = |value: &SymbolicTensor<PartialStructure>| {
        value
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(Default::default()),
                ..Default::default()
            })
            .unwrap()
            .canonize()
            .unwrap()
            .expression
    };
    let a = leaf("word_adjoint::abar", &[bis!(4, 41)]);
    let b = leaf("word_adjoint::b", &[bis!(4, 42)]);
    let word = SymbolicTensor::word(head, &[a.clone(), b.clone()]).unwrap();
    let reversed = SymbolicTensor::word(head, &[bar(&b).unwrap(), bar(&a).unwrap()]).unwrap();
    for (coefficient, adjoint_coefficient) in [
        (gamma!(41, 42, 43), gamma!(42, 41, 43)),
        (
            function!(AGS.projp, bis!(4, 41), bis!(4, 42)),
            function!(AGS.projm, bis!(4, 42), bis!(4, 41)),
        ),
        (
            gamma!(41, 44, 43) * function!(AGS.projp, bis!(4, 44), bis!(4, 42)),
            function!(AGS.projm, bis!(4, 42), bis!(4, 44)) * gamma!(44, 41, 43),
        ),
    ] {
        let coefficient = SymbolicTensor::infer(coefficient).unwrap();
        let expected = SymbolicTensor::infer(adjoint_coefficient)
            .unwrap()
            .multiply_word(head, &reversed)
            .unwrap();
        let source = coefficient.multiply_word(head, &word).unwrap();
        let actual = source.dirac_adjoint_word(head, bar).unwrap();
        assert_eq!(canonical(&actual), canonical(&expected));
        let restored = actual.dirac_adjoint_word(head, bar).unwrap();
        assert_eq!(
            restored.structure.logical_slots(),
            source.structure.logical_slots()
        );
        assert_eq!(canonical(&restored), canonical(&source));
    }
}

#[test]
fn word_adjoint_distributes_sums_with_even_and_odd_terms() {
    let head = head();
    let a = leaf("word_adjoint::a", &[mink!(4, 51)]);
    let b = leaf("word_adjoint::b", &[mink!(4, 51)]);
    let odd = SymbolicTensor::word(head, slice::from_ref(&a))
        .unwrap()
        .multiply(&leaf("word_adjoint::even", &[mink!(4, 51)]))
        .unwrap();
    let word = SymbolicTensor::word(head, &[a, b]).unwrap();
    let scalar = SymbolicTensor::infer(Atom::num((3, 5)) + Atom::i()).unwrap();
    let source = scalar
        .try_add(&scalar.multiply_word(head, &word).unwrap())
        .unwrap()
        .try_add(&odd)
        .unwrap();
    let expected = scalar
        .dirac_adjoint(true)
        .unwrap()
        .try_add(
            &scalar
                .dirac_adjoint(true)
                .unwrap()
                .multiply_word(head, &word.dirac_adjoint_word(head, bar).unwrap())
                .unwrap(),
        )
        .unwrap()
        .try_add(&odd.dirac_adjoint_word(head, bar).unwrap())
        .unwrap();
    let actual = source.dirac_adjoint_word(head, bar).unwrap();
    assert!(actual.same_value(&expected));
    assert!(
        actual
            .dirac_adjoint_word(head, bar)
            .unwrap()
            .same_value(&source)
    );
}

#[test]
fn word_adjoint_materializes_open_coefficients_and_retains_ranked_zero() {
    let head = head();
    let representation = Minkowski {}.new_rep(4).cast::<LibraryRep>();
    let coefficient = SymbolicTensor::from_declared_parts(
        function!(
            symbol!("word_adjoint::real_even"; Real; tag = &SPENSO_TAG.tensor),
            representation.to_symbolic([])
        ),
        PartialStructure::from_logical_slots([representation.slot(PartialIndex::open(0))]),
    )
    .unwrap();
    let word = SymbolicTensor::word(
        head,
        &[
            leaf("word_adjoint::a", &[mink!(4, 61)]),
            leaf("word_adjoint::b", &[mink!(4, 61)]),
        ],
    )
    .unwrap();
    let source = coefficient.multiply_word(head, &word).unwrap();
    let actual = source.dirac_adjoint_word(head, bar).unwrap();
    assert_eq!(actual.rank(), 1);
    assert!(actual.structure.open_positions().is_empty());
    let PartialIndex::Explicit(index) = actual.structure.logical_slots()[0].aind else {
        unreachable!()
    };
    let named_source = source
        .reindex_interface_ports(&HashMap::from([(0, index)]))
        .unwrap();
    assert!(
        actual
            .dirac_adjoint_word(head, bar)
            .unwrap()
            .same_value(&named_source)
    );
    let zero = source.with_algebra_result(Atom::Zero).unwrap();
    let actual = zero
        .dirac_adjoint_word(head, |_| panic!("zero has no generators"))
        .unwrap();
    assert!(actual.is_zero());
    assert_eq!(actual.rank(), 1);
    assert!(actual.structure.open_positions().is_empty());
}

#[test]
fn word_adjoint_rejects_wrong_callback_interfaces_words_and_unsupported_fields() {
    let head = head();
    let a = leaf("word_adjoint::a", &[mink!(4, 71)]);
    let word = SymbolicTensor::word(head, slice::from_ref(&a)).unwrap();
    assert!(
        word.dirac_adjoint_word(head, |_| Ok(leaf("word_adjoint::abar", &[mink!(4, 72)])))
            .is_err()
    );
    assert!(
        word.dirac_adjoint_word(head, |generator| Ok(SymbolicTensor::word(
            head,
            &[bar(generator)?],
        )?))
        .is_err()
    );
    let even_composite = leaf("word_adjoint::even_c", &[mink!(4, 81)])
        .multiply(&leaf("word_adjoint::even_d", &[mink!(4, 81)]))
        .unwrap();
    assert!(
        word.dirac_adjoint_word(head, |generator| {
            Ok(bar(generator)?.multiply(&even_composite)?)
        })
        .is_err()
    );
    for unsupported in ["word_adjoint::majorana", "word_adjoint::ghost"] {
        let word = SymbolicTensor::word(head, &[leaf(unsupported, &[mink!(4, 71)])]).unwrap();
        assert!(matches!(
            word.dirac_adjoint_word(head, bar),
            Err(AdjointError::Tensor(TensorInferenceError::Invalid(reason)))
                if reason == "unsupported field policy"
        ));
    }
    let b = SymbolicTensor::word(head, &[leaf("word_adjoint::b", &[mink!(4, 72)])]).unwrap();
    assert!(
        word.multiply(&b)
            .unwrap()
            .dirac_adjoint_word(head, bar)
            .is_err()
    );
    let foreign_head = symbol!("word_adjoint::foreign"; Flat, Linear, Antisymmetric;
        tag = &SPENSO_TAG.tensor_word);
    let foreign = SymbolicTensor::word(foreign_head, &[a]).unwrap();
    assert!(foreign.dirac_adjoint_word(head, bar).is_err());
    assert!(matches!(
        word.dirac_adjoint(true),
        Err(AdjointError::ImplicitTensorWordAdjoint)
    ));
}

#[test]
fn word_adjoint_retains_hidden_occurrences_for_zero_and_nonzero_results() {
    let head = head();
    let a = leaf("word_adjoint_history::a", &[mink!(4, 991)]);
    let zero = SymbolicTensor::word(head, &[a.clone(), a]).unwrap();
    // Additive terms bind their own internal contractions. Use the native
    // algebra-result owner to retain hidden history in a nonzero scalar.
    let scalar = zero.with_algebra_result(Atom::one()).unwrap();
    let odd =
        SymbolicTensor::word(head, &[leaf("word_adjoint_history::b", &[mink!(4, 992)])]).unwrap();
    let nonzero = scalar.multiply_word(head, &odd).unwrap();
    let probe = leaf("word_adjoint_history::probe", &[mink!(4, 991)]);
    for source in [zero, nonzero] {
        assert!(
            source
                .explicit_indices()
                .unwrap()
                .any(|index| index == AbstractIndex::Normal(991))
        );
        let adjoint = source
            .dirac_adjoint_word(head, |generator| Ok(generator.clone()))
            .unwrap();
        assert!(matches!(
            adjoint.multiply_word(head, &probe),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
        ));
    }
}
