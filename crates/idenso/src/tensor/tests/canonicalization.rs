use crate::{IndexTooling, test_support::test_initialize};
use spenso::{antisym, bracket, broadcast_symbol, cyclic, dind, euc, lor, mink, sym, tensor};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
};

use spenso::structure::abstract_index::AbstractIndex;

fn canonicalize(expression: &Atom) -> Atom {
    expression.canonize::<AbstractIndex>(AbstractIndex::Dummy)
}

fn odd_cycle() -> Atom {
    let a = mink!(4, signed_cycle_a);
    let b = mink!(4, signed_cycle_b);
    let c = mink!(4, signed_cycle_c);
    antisym!(a.clone(), b.clone()) * antisym!(b.clone(), c.clone()) * antisym!(c, a)
}

#[test]
fn odd_antisymmetric_cycle_canonicalizes_to_zero() {
    test_initialize();
    assert!(canonicalize(&odd_cycle()).is_zero());
}

#[test]
fn nested_odd_summand_is_removed() {
    test_initialize();
    let scalar = tensor!(signed_nested_scalar);
    let expression = &scalar * (Atom::num(1) + odd_cycle());

    assert_eq!(canonicalize(&expression), scalar);
}

#[test]
fn factored_expression_is_preserved_when_no_term_is_removed() {
    test_initialize();
    let a = mink!(4, signed_factored_a);
    let b = mink!(4, signed_factored_b);
    let scalar = tensor!(signed_factored_scalar);
    let expression = scalar * (Atom::num(1) + antisym!(a, b).pow(2));
    let canonical = canonicalize(&expression);

    assert!(!canonical.is_zero());
    assert!(matches!(canonical.as_view(), AtomView::Mul(product)
        if product.iter().any(|factor| matches!(factor, AtomView::Add(_)))));
}

#[test]
fn sum_inside_function_remains_opaque() {
    test_initialize();
    let expression = tensor!(signed_opaque, Atom::num(1) + odd_cycle());
    let canonical = canonicalize(&expression);

    assert!(!canonical.is_zero());
    assert!(matches!(canonical.as_view(), AtomView::Fun(function)
        if function.iter().any(|argument| matches!(argument, AtomView::Add(_)))));
}

#[test]
fn broadcast_boundary_does_not_hide_independent_zero() {
    test_initialize();
    let a = mink!(4, signed_opaque_a);
    let b = mink!(4, signed_opaque_b);
    let c = mink!(4, signed_opaque_c);
    let inner = antisym!(a.clone(), b.clone()) * antisym!(b.clone(), c.clone()) * antisym!(c, a);
    let opaque = function!(
        broadcast_symbol!(signed_opaque_broadcast),
        Atom::num(1) + inner
    );

    // A nonlinear function must retain the vanishing summand rather than
    // applying its zero before evaluating the enclosing function.
    let canonical = canonicalize(&opaque);
    assert!(!canonical.is_zero());
    assert!(matches!(canonical.as_view(), AtomView::Fun(function)
        if function.iter().any(|argument| matches!(argument, AtomView::Add(_)))));

    // The opaque function boundary must not suppress an independent zero.
    assert!(canonicalize(&(odd_cycle() * opaque)).is_zero());
}

#[test]
fn even_antisymmetric_automorphism_survives_canonicalization() {
    test_initialize();
    let a = mink!(4, signed_square_a);
    let b = mink!(4, signed_square_b);

    assert!(!canonicalize(&antisym!(a, b).pow(2)).is_zero());
}

#[test]
fn power_copies_closed_components_as_a_whole() {
    test_initialize();

    // Keep the closed contraction as one Power child until Spenso parses it.
    assert!(canonicalize(&bracket!(odd_cycle()).pow(2)).is_zero());
}

#[test]
fn external_indices_prevent_slot_automorphism() {
    test_initialize();
    let a = mink!(4, signed_external_a);
    let b = mink!(4, signed_external_b);

    assert!(!canonicalize(&antisym!(a, b)).is_zero());
}

#[test]
fn ordered_slots_prevent_antisymmetric_automorphism() {
    test_initialize();
    let a = mink!(4, signed_ordered_a);
    let b = mink!(4, signed_ordered_b);
    let contraction = tensor!(signed_ordered, a.clone(), b.clone()) * antisym!(a, b);

    assert!(!canonicalize(&contraction).is_zero());
}

#[test]
fn symmetric_slots_allow_antisymmetric_automorphism() {
    test_initialize();
    let a = mink!(4, signed_symmetric_a);
    let b = mink!(4, signed_symmetric_b);
    let contraction = sym!(a.clone(), b.clone()) * antisym!(a, b);

    assert!(canonicalize(&contraction).is_zero());
}

#[test]
fn different_representation_groups_break_odd_automorphism() {
    test_initialize();
    let a = mink!(4, signed_group_a);
    let b = euc!(4, signed_group_b);
    let c = lor!(4, signed_group_c);
    let contraction =
        antisym!(a.clone(), b.clone()) * antisym!(b, c.clone()) * antisym!(dind!(c), a);

    assert!(!canonicalize(&contraction).is_zero());
}

#[test]
fn cyclic_symmetry_preserves_only_rotations() {
    test_initialize();
    let a = mink!(4, signed_cyclic_a);
    let b = mink!(4, signed_cyclic_b);
    let c = mink!(4, signed_cyclic_c);
    let rank_three =
        cyclic!(a.clone(), b.clone(), c.clone()) * antisym!(a.clone(), b.clone(), c.clone());

    // A three-cycle is even, so the permitted cyclic rotation has no sign.
    assert!(!canonicalize(&rank_three).is_zero());

    let d = mink!(4, signed_cyclic_d);
    let rank_four = cyclic!(a.clone(), b.clone(), c.clone(), d.clone()) * antisym!(a, b, c, d);

    // A four-cycle is odd. Reflections are not part of cyclic symmetry.
    assert!(canonicalize(&rank_four).is_zero());
}
