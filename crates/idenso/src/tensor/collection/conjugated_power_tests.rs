use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

use super::{CollectionMode, SELECTED_OPENS, SymbolicTensor};
use spenso::{
    algebra::ScalarMul,
    network::{library::function_lib::INBUILTS, tags::SPENSO_TAG},
    structure::{
        abstract_index::AbstractIndex,
        partial::{PartialStructure, PartialStructureExt},
        representation::{LibraryRep, Minkowski, RepName},
        slot::IsAbstractSlot,
    },
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Symbol},
    function, symbol,
};

type Tensor = SymbolicTensor<PartialStructure>;

fn closed_conjugation(name: &str, index: usize, dual: bool, coefficient: &Atom) -> Tensor {
    let rep = LibraryRep::new_dual(&format!("conjugated_power::{name}_rep"))
        .unwrap()
        .new_rep(3);
    let rep = if dual { rep.dual() } else { rep };
    let index = AbstractIndex::Normal(index);
    let head = SPENSO_TAG.tensor_symbol(&format!("conjugated_power::{name}"));
    let direct = function!(head, rep.slot::<AbstractIndex, _>(index).to_atom());
    let opposite = function!(head, rep.dual().slot::<AbstractIndex, _>(index).to_atom());
    Tensor::infer(INBUILTS.conj(coefficient * INBUILTS.conj(direct) * opposite)).unwrap()
}

fn independent_square(value: &Tensor, left: Symbol, right: Symbol) -> Tensor {
    value
        .wrap_indices(left, true)
        .unwrap()
        .multiply(&value.wrap_indices(right, true).unwrap())
        .unwrap()
}

fn same(actual: &Tensor, expected: &Tensor) {
    assert_eq!(
        actual.structure().logical_slots(),
        expected.structure().logical_slots()
    );
    let actual = actual.canonize().unwrap();
    let expected = expected.canonize().unwrap();
    assert!(actual.same_value(&expected));
    assert!(actual.canonize().unwrap().same_value(&actual));
    assert_eq!(
        Tensor::infer(actual.expression().clone())
            .unwrap()
            .structure()
            .logical_slots(),
        actual.structure().logical_slots()
    );
}

#[test]
fn closed_conjugation_powers_match_independently_scoped_copies() {
    crate::test_support::test_initialize();
    let two = Tensor::infer(Atom::num(2)).unwrap();
    for dual in [false, true] {
        let value = closed_conjugation("norm", 99131, dual, &Atom::one());
        let source = value.pow(&two).unwrap();
        let before = source.expression().clone();
        let expected = independent_square(
            &value,
            symbol!("conjugated_power::left"),
            symbol!("conjugated_power::right"),
        );
        assert!(source.is_scalar());
        same(&source, &expected);
        assert_eq!(source.expression(), &before);
    }
}

#[test]
fn closed_conjugation_power_keeps_complex_coefficients_and_external_ports() {
    crate::test_support::test_initialize();
    let complex = Atom::var(symbol!("conjugated_power::complex"; Scalar));
    let value = closed_conjugation("weighted", 99137, false, &complex);
    let port = Minkowski {}
        .new_rep(4)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99139));
    let spectator = Tensor::infer(function!(
        SPENSO_TAG.tensor_symbol("conjugated_power::spectator"),
        port.to_atom()
    ))
    .unwrap();
    let source = value
        .pow(&Tensor::infer(Atom::num(2)).unwrap())
        .unwrap()
        .multiply(&spectator)
        .unwrap()
        .scalar_mul(&complex.conj())
        .unwrap();
    let expected = independent_square(
        &value,
        symbol!("conjugated_power::weighted_left"),
        symbol!("conjugated_power::weighted_right"),
    )
    .multiply(&spectator)
    .unwrap()
    .scalar_mul(&complex.conj())
    .unwrap();
    same(&source, &expected);
    assert_eq!(
        source.structure().logical_slots(),
        spectator.structure().logical_slots()
    );
    let symbols = source.expression().get_all_symbols(true);
    assert!(symbols.contains(&Symbol::CONJ));
    assert!(symbols.contains(&complex.get_symbol().unwrap()));
}

#[test]
fn conjugated_power_opening_keeps_unknown_functions_metadata_and_words_opaque() {
    crate::test_support::test_initialize();
    let value = closed_conjugation("opaque", 99143, false, &Atom::one());
    let metadata = symbol!("conjugated_power::metadata"; Scalar);
    let broadcast = spenso::broadcast_symbol!("conjugated_power::unknown_broadcast");
    let word_head = symbol!("conjugated_power::word"; Flat, Linear, Antisymmetric;
        tag = &SPENSO_TAG.tensor_word);
    let operands = ["conjugated_power::a", "conjugated_power::b"].map(|name| {
        Tensor::infer(FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name)).finish()).unwrap()
    });
    let word = Tensor::word(word_head, &operands).unwrap();
    for atom in [
        function!(metadata, value.expression()),
        function!(broadcast, value.expression()),
        INBUILTS.conj(word.expression()),
    ] {
        let source = Tensor::infer(atom)
            .unwrap()
            .pow(&Tensor::infer(Atom::num(2)).unwrap())
            .unwrap();
        SELECTED_OPENS.with(|count| count.set(0));
        let result = source
            .collect_with_map(
                CollectionMode::Factored,
                None,
                |_| true,
                |value, _, _| Ok(value),
            )
            .unwrap();
        assert!(result.same_value(&source));
        assert_eq!(SELECTED_OPENS.with(|count| count.get()), 0);
    }
}

#[test]
fn closed_conjugation_power_does_not_replay_sensitive_normalizers() {
    crate::test_support::test_initialize();
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&calls);
    let head = spenso::tensor_symbol!(
        "conjugated_power::callback",
        norm = move |_, _| {
            observed.fetch_add(1, Ordering::Relaxed);
        }
    );
    let rep = LibraryRep::new_dual("conjugated_power::callback_rep")
        .unwrap()
        .new_rep(3);
    let index = AbstractIndex::Normal(99149);
    let body = function!(head, rep.slot::<AbstractIndex, _>(index).to_atom())
        * function!(head, rep.dual().slot::<AbstractIndex, _>(index).to_atom());
    let source = Tensor::infer(INBUILTS.conj(body))
        .unwrap()
        .pow(&Tensor::infer(Atom::num(2)).unwrap())
        .unwrap();
    calls.store(0, Ordering::Relaxed);
    SELECTED_OPENS.with(|count| count.set(0));
    let result = source
        .collect_with_map(
            CollectionMode::Factored,
            None,
            |_| true,
            |value, _, _| Ok(value),
        )
        .unwrap();
    assert!(result.same_value(&source));
    assert_eq!(SELECTED_OPENS.with(|count| count.get()), 0);
    assert_eq!(calls.load(Ordering::Relaxed), 0);
}
