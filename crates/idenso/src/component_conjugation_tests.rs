use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::IndexTooling;
use crate::tensor::SymbolicTensor;
use crate::test_support::test_initialize;
use spenso::broadcast_symbol;
use spenso::network::library::function_lib::INBUILTS;
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialStructure, PartialStructureExt};
use spenso::structure::representation::{LibraryRep, RepName};
use spenso::structure::slot::IsAbstractSlot;
use symbolica::atom::{Atom, AtomCore, FunctionBuilder, SymbolAttribute, SymbolBuilder};
use symbolica::{function, symbol, wrap_symbol};

type Tensor = SymbolicTensor<PartialStructure>;

fn components(name: &str, dual: bool, index: usize) -> (Atom, Atom) {
    let rep = LibraryRep::new_dual(&format!("component_conjugation::{name}_rep"))
        .unwrap()
        .new_rep(3);
    let rep = if dual { rep.dual() } else { rep };
    let index = AbstractIndex::Normal(index);
    let head = SPENSO_TAG.tensor_symbol(&format!("component_conjugation::{name}"));
    (
        function!(head, rep.slot::<AbstractIndex, _>(index).to_atom()),
        function!(head, rep.dual().slot::<AbstractIndex, _>(index).to_atom()),
    )
}

fn same(actual: &Tensor, expected: &Tensor) {
    assert_eq!(
        actual.structure().logical_slots(),
        expected.structure().logical_slots()
    );
    let canonical = actual.canonize().unwrap();
    let expected_canonical = expected.canonize().unwrap();
    assert!(
        canonical.same_value(&expected_canonical),
        "input={}\nactual={}\nexpected={}",
        actual.expression(),
        canonical.expression(),
        expected_canonical.expression()
    );
    assert!(canonical.canonize().unwrap().same_value(&canonical));
    assert_eq!(
        Tensor::infer(canonical.expression().clone())
            .unwrap()
            .structure()
            .logical_slots(),
        canonical.structure().logical_slots()
    );
}

#[test]
fn closed_broadcast_conjugation_preserves_component_orientation_and_complex_weights() {
    test_initialize();
    let z = Atom::var(symbol!("component_conjugation::z"; Scalar));
    for dual in [false, true] {
        let (direct, opposite) = components("weighted", dual, 99203);
        for (weight, conjugated_weight) in [
            (Atom::one(), Atom::one()),
            (z.clone(), z.conj()),
            (Atom::i(), -Atom::i()),
            (&z + Atom::i(), z.conj() - Atom::i()),
        ] {
            let argument = weight * INBUILTS.conj(&direct) * &opposite;
            let literal = INBUILTS.conj(&argument);
            // Raw construction retains tensor broadcast semantics. Only the
            // explicit canonicalizer normalizes this closed mathematical value.
            assert_ne!(literal, argument);
            let source = Tensor::infer(literal.clone()).unwrap();
            let expected =
                Tensor::infer(conjugated_weight * &direct * INBUILTS.conj(&opposite)).unwrap();
            same(&source, &expected);
            assert_eq!(source.expression(), &literal);
        }
    }
}

#[test]
fn conjugated_norm_powers_match_independently_fresh_norms_and_sums() {
    test_initialize();
    let (h, h_dual) = components("h", false, 99209);
    let (k, k_dual) = components("k", false, 99211);
    let p = Tensor::infer(&h * INBUILTS.conj(&h_dual)).unwrap();
    let q = Tensor::infer(&k * INBUILTS.conj(&k_dual)).unwrap();
    for (argument, expected_base) in [
        (INBUILTS.conj(&h) * &h_dual, p.clone()),
        (
            INBUILTS.conj(&h) * &h_dual + INBUILTS.conj(&k) * &k_dual,
            p.try_add(&q).unwrap(),
        ),
    ] {
        let source = Tensor::infer(INBUILTS.conj(argument))
            .unwrap()
            .pow(&Tensor::infer(Atom::num(2)).unwrap())
            .unwrap();
        let expected = expected_base
            .wrap_indices(symbol!("component_conjugation::left_copy"), true)
            .unwrap()
            .multiply(
                &expected_base
                    .wrap_indices(symbol!("component_conjugation::right_copy"), true)
                    .unwrap(),
            )
            .unwrap();
        same(&source, &expected);
    }
}

#[test]
fn real_scalar_metadata_with_hidden_universal_conjugation_remains_opaque() {
    test_initialize();
    let (h, h_dual) = components("metadata_norm", false, 99217);
    let (private, _) = components("private_payload", false, 99221);
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&calls);
    let metadata_hook = symbol!("component_conjugation::real_metadata_hook"; Scalar, Real;
        norm = move |_, _| { observed.fetch_add(1, Ordering::Relaxed); });
    let hidden = private.conj();
    for metadata_head in [
        symbol!("component_conjugation::real_metadata"; Scalar, Real),
        metadata_hook,
        symbol!("component_conjugation::untagged_real_metadata"; Real),
    ] {
        let metadata = function!(metadata_head, &hidden);
        // Constructing this comparison fixture can itself invoke the hook.
        let forbidden_metadata = function!(metadata_head, INBUILTS.conj(&private));
        // The first case is intrinsic, so callback refusal alone cannot protect
        // this Real+Scalar container from an unguarded global CONJ conversion.
        assert_eq!(metadata.conj(), metadata);
        let literal = INBUILTS.conj(&metadata * INBUILTS.conj(&h) * &h_dual);
        assert_eq!(literal.get_symbol(), Some(INBUILTS.conj));
        let source = Tensor::infer(literal)
            .unwrap()
            .pow(&Tensor::infer(Atom::num(2)).unwrap())
            .unwrap();
        calls.store(0, Ordering::Relaxed);
        let canonical = source.canonize().unwrap();
        assert!(canonical.expression().contains(&metadata));
        assert!(canonical.expression().contains(&hidden));
        assert!(!canonical.expression().contains(&forbidden_metadata));
        assert_eq!(
            canonical.structure().logical_slots(),
            source.structure().logical_slots()
        );
        assert_eq!(calls.load(Ordering::Relaxed), 0);
    }
}

#[test]
fn real_tensor_metadata_keeps_private_universal_conjugation() {
    test_initialize();
    let rep = LibraryRep::new_dual("component_conjugation::real_tensor_rep")
        .unwrap()
        .new_rep(3);
    let index = AbstractIndex::Normal(99225);
    let head = SymbolBuilder::new(wrap_symbol!("component_conjugation::real_tensor"))
        .with_attributes(vec![SymbolAttribute::Real])
        .with_tags([&SPENSO_TAG.tensor])
        .build()
        .unwrap();
    let scalar = Atom::var(symbol!("component_conjugation::private_scalar"));
    let hidden = scalar.conj();
    let port = rep.slot::<AbstractIndex, _>(index).to_atom();
    let real = function!(head, &hidden, &port);
    let forbidden = function!(head, INBUILTS.conj(&scalar), &port);
    let opposite = function!(
        SPENSO_TAG.tensor_symbol("component_conjugation::real_tensor_partner"),
        rep.dual().slot::<AbstractIndex, _>(index).to_atom()
    );
    let source = Tensor::infer(INBUILTS.conj(&real * &opposite)).unwrap();
    let expected = Tensor::infer(&real * INBUILTS.conj(&opposite)).unwrap();
    same(&source, &expected);
    let canonical = source.canonize().unwrap();
    assert!(canonical.expression().contains(&hidden));
    assert!(!canonical.expression().contains(&forbidden));
}

#[test]
fn component_normalization_keeps_open_broadcasts_words_and_unknown_functions_opaque() {
    test_initialize();
    let (direct, opposite) = components("open", false, 99227);
    let source = Tensor::infer(INBUILTS.conj(&direct)).unwrap();
    assert_eq!(source.rank(), 1);
    assert!(source.canonize().unwrap().same_value(&source));
    // Full physical conjugation still dualizes the representation exactly once.
    let physical = Tensor::infer(direct.spenso_conj()).unwrap();
    same(&physical, &Tensor::infer(INBUILTS.conj(&opposite)).unwrap());

    let word_head = symbol!("component_conjugation::word"; Flat, Linear, Antisymmetric;
        tag = &SPENSO_TAG.tensor_word);
    let operands = ["component_conjugation::a", "component_conjugation::b"].map(|name| {
        Tensor::infer(FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name)).finish()).unwrap()
    });
    let word = Tensor::word(word_head, &operands).unwrap();
    let wrapped = Tensor::infer(INBUILTS.conj(word.expression())).unwrap();
    assert!(wrapped.canonize().unwrap().same_value(&wrapped));

    let broadcast = broadcast_symbol!("component_conjugation::unknown");
    let atom = function!(
        broadcast,
        Atom::var(symbol!("component_conjugation::x"; Scalar))
    );
    let unknown = Tensor::infer(INBUILTS.conj(&atom)).unwrap();
    assert!(unknown.canonize().unwrap().expression().contains(&atom));
}
