use crate::{IndexTooling, test_support::test_initialize};
use spenso::{antisym, bracket, broadcast_symbol, cyclic, dind, euc, lor, mink, sym, tensor};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function,
};

use spenso::structure::abstract_index::AbstractIndex;

fn canonicalize(expression: &Atom) -> Atom {
    expression
        .canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .expect("test expression should canonicalize")
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
fn canonicalization_preserves_spectators_after_pruning() {
    test_initialize();
    assert!(canonicalize(&odd_cycle()).is_zero());
    let spectator = tensor!(signed_metric_spectator_a) + tensor!(signed_metric_spectator_b);
    let expression = &spectator * (Atom::one() + odd_cycle());
    let canonical = canonicalize(&expression);
    assert_eq!(canonical, spectator);
    assert_eq!(canonicalize(&canonical), canonical);
}

#[test]
fn nested_odd_summand_is_removed() {
    test_initialize();
    let scalar = tensor!(signed_nested_scalar);
    let expression = &scalar * (Atom::num(1) + odd_cycle());

    assert_eq!(canonicalize(&expression), scalar);
}

#[test]
fn pruning_nested_zero_preserves_unrelated_sum_factors() {
    test_initialize();
    let spectator = tensor!(signed_spectator_a) + tensor!(signed_spectator_b);
    let expected = tensor!(signed_nested_scalar) * spectator;
    let expression = &expected * (Atom::one() + odd_cycle());

    assert_eq!(canonicalize(&expression), expected);
}

#[test]
fn products_of_tensor_sums_keep_closed_indices_and_factorization() {
    test_initialize();
    let expression = |i: Atom, j: Atom| {
        (tensor!(factorized_a, i.clone()) + tensor!(factorized_b, i.clone()))
            * (tensor!(factorized_c, i.clone()) + tensor!(factorized_d, i))
            + tensor!(factorized_e, j.clone()) * tensor!(factorized_f, j)
    };
    let original = expression(mink!(4, factorized_i), mink!(4, factorized_j));
    let renamed = expression(mink!(4, factorized_l), mink!(4, factorized_k));
    let canonical = canonicalize(&original);

    assert_eq!(canonical, canonicalize(&renamed));
    assert_eq!(canonical, canonicalize(&canonical));
    assert!(
        canonical
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .is_empty()
    );
    let AtomView::Add(sum) = canonical.as_view() else {
        panic!("the two original scalar summands must remain present");
    };
    assert!(sum.iter().any(|term| matches!(term, AtomView::Mul(product)
        if product.iter().filter(|factor| matches!(factor, AtomView::Add(_))).count() == 2)));
    assert_ne!(canonical, canonicalize(&(original + Atom::one())));
}

#[test]
fn products_of_tensor_sums_preserve_genuine_external_indices() {
    test_initialize();
    // Native allocation must reserve true external names even when they sort
    // before all dummy candidates, as well as when they sort after them.
    for external in [
        mink!(4, a_factorized_external),
        mink!(4, z_factorized_external),
    ] {
        let expression = |i: Atom| {
            (tensor!(factorized_open_a, i.clone()) + tensor!(factorized_open_b, i.clone()))
                * (tensor!(factorized_open_c, i.clone(), external.clone())
                    + tensor!(factorized_open_d, i, external.clone()))
        };
        let original = expression(mink!(4, factorized_open_i));
        let renamed = expression(mink!(4, factorized_open_j));
        let canonical = canonicalize(&original);

        assert_eq!(canonical, canonicalize(&renamed));
        assert_eq!(canonical, canonicalize(&canonical));
        assert_eq!(
            canonical.list_dangling::<AbstractIndex>().unwrap(),
            vec![external]
        );
        assert!(matches!(canonical.as_view(), AtomView::Mul(product)
            if product.iter().filter(|factor| matches!(factor, AtomView::Add(_))).count() == 2));
    }
}

#[test]
fn factorized_canonicalization_rejects_mismatched_external_indices() {
    test_initialize();
    let expression = tensor!(factorized_bad_a, mink!(4, factorized_bad_i))
        + tensor!(factorized_bad_b, mink!(4, factorized_bad_j));

    assert!(std::panic::catch_unwind(|| canonicalize(&expression)).is_err());
}

#[test]
fn factorized_canonicalization_rejects_overcontracted_indices() {
    test_initialize();
    let index = mink!(4, factorized_overcontracted);
    let expression = tensor!(factorized_over_a, index.clone())
        * tensor!(factorized_over_b, index.clone())
        * tensor!(factorized_over_c, index);

    assert!(std::panic::catch_unwind(|| canonicalize(&expression)).is_err());
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

#[test]
fn inverse_scalar_contractions_preserve_their_boundary_during_canonicalization() {
    test_initialize();
    let p = spenso::vector_symbol!("inverse_scalar_canon::p");
    let q = spenso::vector_symbol!("inverse_scalar_canon::q");
    let minkowski = spenso::structure::representation::Minkowski {};
    use spenso::structure::representation::RepName;
    let denominator = minkowski.new_rep(4).inner_product(p, q) - Atom::one();
    let numerator = tensor!(inverse_scalar_a, mink!(4, inverse_scalar_i))
        * tensor!(inverse_scalar_b, mink!(4, inverse_scalar_i));
    let expression = numerator / denominator;
    let canonical = canonicalize(&expression);
    assert_eq!(canonical, canonicalize(&canonical));
    assert!(
        canonical
            .list_dangling::<AbstractIndex>()
            .unwrap()
            .is_empty()
    );
    assert!(!canonical.to_plain_string().contains("scalar_store"));
}

#[test]
fn pruning_zero_arm_exposes_outer_signed_contraction() {
    test_initialize();
    let a = mink!(4, signed_exposed_a);
    let b = mink!(4, signed_exposed_b);
    let sum = sym!(a.clone(), b.clone())
        + odd_cycle() * tensor!(signed_exposed_other, a.clone(), b.clone());
    let expression = antisym!(a, b) * sum;
    assert!(canonicalize(&expression).is_zero());
}

#[test]
fn pruning_open_zero_arm_preserves_free_and_dual_ports() {
    test_initialize();
    let a = lor!(4, signed_pruned_dual_a);
    let b = mink!(4, signed_pruned_free_b);
    let survivor = tensor!(signed_pruned_survivor, a.clone(), b.clone());
    let sum = &survivor + odd_cycle() * tensor!(signed_pruned_zero, a.clone(), b.clone());
    let closing = tensor!(signed_pruned_close, dind!(a));
    let expected = canonicalize(&(&closing * survivor));
    let result = canonicalize(&(closing * sum));
    assert_eq!(result, expected);
    assert_eq!(canonicalize(&result), result);
    assert_eq!(result.list_dangling::<AbstractIndex>().unwrap(), vec![b]);
}

#[test]
fn pruning_signed_term_does_not_mutate_a_shared_tensor_occurrence() {
    test_initialize();
    let a = mink!(4, signed_cycle_a);
    let b = mink!(4, signed_cycle_b);
    let even = antisym!(a, b).pow(2);
    let expected = canonicalize(&even);
    assert!(!expected.is_zero());
    assert_eq!(canonicalize(&(odd_cycle() + even)), expected);
}

#[test]
fn unparseable_indices_are_rejected_before_signed_zero_pruning() {
    test_initialize();
    use crate::{CanonicalizationError, Cookable, NetworkToolingError};
    use symbolica::symbol;
    let slot = |n| mink!(4, function!(symbol!("signed_raw_index"), Atom::num(n)));
    let expression =
        antisym!(slot(1), slot(2)) * antisym!(slot(2), slot(3)) * antisym!(slot(3), slot(1));
    assert!(!expression.is_zero());
    let raw = expression.canonize::<AbstractIndex>(AbstractIndex::Dummy);
    assert!(
        matches!(
            raw,
            Err(CanonicalizationError::Network(
                NetworkToolingError::Parse { .. }
            ))
        ),
        "structured index admission result: {raw:?}"
    );
    // A caller with structured labels explicitly admits them once. Zero and
    // nonzero terms obey the same abstract-index boundary.
    assert!(canonicalize(&expression.cook_indices()).is_zero());
}
