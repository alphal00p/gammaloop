use spenso::{
    antisym, bracket, broadcast_symbol, cyclic, dind, euc, lor, mink,
    network::library::symbolic::ETS,
    slot,
    structure::{abstract_index::AbstractIndex, representation::LibrarySlot, slot::IsAbstractSlot},
    sym, tensor, tensor_symbol, vector,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    function, symbol,
};

use crate::{
    IndexTooling,
    tensor::{CanonicalizationError, SymbolicNet, SymbolicNetExt},
    test_support::test_initialize,
};

fn canonicalize(expression: &Atom) -> Result<Atom, CanonicalizationError> {
    expression.try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
}

fn unsigned_reference(
    expression: &Atom,
    indices: impl IntoIterator<Item = LibrarySlot<AbstractIndex>>,
) -> Atom {
    expression
        .canonize_tensors(
            indices
                .into_iter()
                .map(|slot| (slot.to_atom(), slot.rep().base())),
        )
        .unwrap()
        .canonical_form
}

fn odd_cycle() -> Atom {
    let rep = test_initialize().mink4;
    let a = slot!(rep, signed_cycle_a);
    let b = slot!(rep, signed_cycle_b);
    let c = slot!(rep, signed_cycle_c);
    let antisymmetric = tensor_symbol!(signed_cycle_tensor; Antisymmetric);

    function!(antisymmetric, a.to_atom(), b.to_atom())
        * function!(antisymmetric, b.to_atom(), c.to_atom())
        * function!(antisymmetric, c.to_atom(), a.to_atom())
}

fn structural_odd_cycle() -> Atom {
    let a = mink!(4, signed_structural_cycle_a);
    let b = mink!(4, signed_structural_cycle_b);
    let c = mink!(4, signed_structural_cycle_c);
    antisym!(a.clone(), b.clone()) * antisym!(b.clone(), c.clone()) * antisym!(c, a)
}

#[test]
fn odd_antisymmetric_stabilizer_zeroes_the_term() {
    assert!(canonicalize(&odd_cycle()).unwrap().is_zero());
}

#[test]
fn even_antisymmetric_stabilizer_survives() {
    let rep = test_initialize().mink4;
    let a = slot!(rep, signed_even_a);
    let b = slot!(rep, signed_even_b);
    let antisymmetric = tensor_symbol!(signed_even_tensor; Antisymmetric);
    let expression = function!(antisymmetric, a.to_atom(), b.to_atom()).pow(2);

    assert!(!canonicalize(&expression).unwrap().is_zero());
}

#[test]
fn higher_positive_powers_are_canonical_and_idempotent() {
    let rep = test_initialize().mink4;
    let a = slot!(rep, signed_power_a);
    let b = slot!(rep, signed_power_b);
    let symmetric = tensor_symbol!(signed_power_tensor; Symmetric);
    for power in [3, 4] {
        let expression = function!(symmetric, a.to_atom(), b.to_atom()).pow(power);
        let once = canonicalize(&expression).unwrap();
        assert_eq!(canonicalize(&once).unwrap(), once);
    }

    let leaf = tensor_symbol!(signed_power_leaf);
    let expression = function!(symmetric, a.to_atom(), b.to_atom()).pow(3)
        * function!(leaf, a.to_atom())
        * function!(leaf, b.to_atom());
    let once = canonicalize(&expression).unwrap();
    assert_eq!(canonicalize(&once).unwrap(), once);
}

#[test]
fn contractions_across_sums_ignore_original_dummy_names() {
    let rep = test_initialize().mink4;
    let left = tensor_symbol!(signed_boundary_left);
    let first = tensor_symbol!(signed_boundary_first);
    let second = tensor_symbol!(signed_boundary_second);
    let build = |slot: LibrarySlot<AbstractIndex>| {
        function!(left, slot.to_atom())
            * (function!(first, slot.to_atom()) + function!(second, slot.to_atom()))
    };

    assert_eq!(
        canonicalize(&build(slot!(rep, signed_boundary_a).cast())).unwrap(),
        canonicalize(&build(slot!(rep, signed_boundary_b).cast())).unwrap()
    );
}

#[test]
fn independent_summands_reuse_canonical_dummy_scope() {
    let rep = test_initialize().mink4;
    let left = tensor_symbol!(signed_sum_scope_left);
    let right = tensor_symbol!(signed_sum_scope_right);
    let term = |slot: LibrarySlot<AbstractIndex>| {
        function!(left, slot.to_atom()) * function!(right, slot.to_atom())
    };
    let first: LibrarySlot<_> = slot!(rep, signed_sum_scope_a).cast();
    let second: LibrarySlot<_> = slot!(rep, signed_sum_scope_b).cast();

    assert_eq!(
        canonicalize(&(term(first) + term(second))).unwrap(),
        canonicalize(&(term(first) + term(first))).unwrap()
    );
}

#[test]
fn antisymmetric_boundary_slots_do_not_create_local_stabilizers() {
    let rep = test_initialize().mink4;
    let i: LibrarySlot<_> = slot!(rep, signed_barrier_i).cast();
    let j: LibrarySlot<_> = slot!(rep, signed_barrier_j).cast();
    let antisymmetric = tensor_symbol!(signed_barrier_antisymmetric; Antisymmetric);
    let first_i = tensor_symbol!(signed_barrier_first_i);
    let first_j = tensor_symbol!(signed_barrier_first_j);
    let second_i = tensor_symbol!(signed_barrier_second_i);
    let second_j = tensor_symbol!(signed_barrier_second_j);
    let expression = function!(antisymmetric, i.to_atom(), j.to_atom())
        * (function!(first_i, i.to_atom()) * function!(first_j, j.to_atom())
            + function!(second_i, i.to_atom()) * function!(second_j, j.to_atom()));

    let canonical = canonicalize(&expression).unwrap();
    assert!(!canonical.is_zero());
    assert_eq!(canonicalize(&canonical).unwrap(), canonical);
}

#[test]
fn unsigned_ordered_and_partial_cyclic_match_symbolica() {
    let rep = test_initialize().mink4;
    let i: LibrarySlot<_> = slot!(rep, AbstractIndex::Dummy(0)).cast();
    let j: LibrarySlot<_> = slot!(rep, AbstractIndex::Dummy(1)).cast();
    let k: LibrarySlot<_> = slot!(rep, AbstractIndex::Dummy(2)).cast();
    let outer = tensor_symbol!(unsigned_partial_outer);
    let leaf = tensor_symbol!(unsigned_partial_leaf);
    let parameter = Atom::var(symbol!("unsigned_parameter"));
    let expression = FunctionBuilder::new(outer)
        .add_arg(parameter)
        .add_arg(cyclic!(i, j, k))
        .finish()
        * function!(leaf, i.to_atom())
        * function!(leaf, j.to_atom())
        * function!(leaf, k.to_atom());

    assert_eq!(
        canonicalize(&expression).unwrap(),
        unsigned_reference(&expression, [i, j, k])
    );

    let symmetric = tensor_symbol!(unsigned_intrinsic_symmetric; Symmetric);
    let expression = function!(symmetric, j.to_atom(), i.to_atom())
        * function!(leaf, i.to_atom())
        * function!(leaf, j.to_atom());
    assert_eq!(
        canonicalize(&expression).unwrap(),
        unsigned_reference(&expression, [i, j])
    );
}

#[test]
fn reflected_cycles_remain_distinct() {
    let rep = test_initialize().mink4;
    let i: LibrarySlot<_> = slot!(rep, AbstractIndex::Dummy(0)).cast();
    let j: LibrarySlot<_> = slot!(rep, AbstractIndex::Dummy(1)).cast();
    let k: LibrarySlot<_> = slot!(rep, AbstractIndex::Dummy(2)).cast();
    let outer = tensor_symbol!(unsigned_reflected_cycle_outer);
    let first = tensor_symbol!(unsigned_reflected_cycle_first);
    let second = tensor_symbol!(unsigned_reflected_cycle_second);
    let third = tensor_symbol!(unsigned_reflected_cycle_third);
    let build = |group| {
        FunctionBuilder::new(outer).add_arg(group).finish()
            * function!(first, i.to_atom())
            * function!(second, j.to_atom())
            * function!(third, k.to_atom())
    };
    let forward = build(cyclic!(i, j, k));
    let reflected = build(cyclic!(i, k, j));
    let canonical_forward = canonicalize(&forward).unwrap();
    let canonical_reflected = canonicalize(&reflected).unwrap();

    assert_eq!(canonical_forward, unsigned_reference(&forward, [i, j, k]));
    assert_eq!(
        canonical_reflected,
        unsigned_reference(&reflected, [i, j, k])
    );
    assert_ne!(canonical_forward, canonical_reflected);

    let intrinsic = tensor_symbol!(unsigned_reflected_cycle_intrinsic; Cyclesymmetric);
    let intrinsic_forward = function!(intrinsic, i.to_atom(), j.to_atom(), k.to_atom())
        * function!(first, i.to_atom())
        * function!(second, j.to_atom())
        * function!(third, k.to_atom());
    let intrinsic_reflected = function!(intrinsic, i.to_atom(), k.to_atom(), j.to_atom())
        * function!(first, i.to_atom())
        * function!(second, j.to_atom())
        * function!(third, k.to_atom());
    let canonical_intrinsic_forward = canonicalize(&intrinsic_forward).unwrap();
    let canonical_intrinsic_reflected = canonicalize(&intrinsic_reflected).unwrap();

    assert_eq!(
        canonical_intrinsic_forward,
        unsigned_reference(&intrinsic_forward, [i, j, k])
    );
    assert_eq!(
        canonical_intrinsic_reflected,
        unsigned_reference(&intrinsic_reflected, [i, j, k])
    );
    assert_ne!(canonical_intrinsic_forward, canonical_intrinsic_reflected);
}

#[test]
fn rank_one_cycles_are_trivial() {
    let rep = test_initialize().mink4;
    let index: LibrarySlot<_> = slot!(rep, signed_rank_one_cycle).cast();
    let outer = tensor_symbol!(signed_rank_one_cycle_outer);
    let leaf = tensor_symbol!(signed_rank_one_cycle_leaf);
    let partial = FunctionBuilder::new(outer).add_arg(cyclic!(index)).finish()
        * function!(leaf, index.to_atom());
    let intrinsic = tensor_symbol!(signed_rank_one_cycle_intrinsic; Cyclesymmetric);
    let intrinsic = function!(intrinsic, index.to_atom()) * function!(leaf, index.to_atom());

    for expression in [partial, intrinsic] {
        let canonical = canonicalize(&expression).unwrap();
        assert_eq!(canonicalize(&canonical).unwrap(), canonical);
    }
}

#[test]
fn zero_summands_are_removed_without_distributing_products() {
    let scalar = Atom::var(symbol!("signed_scalar"));
    let expression = &scalar * (Atom::num(1) + odd_cycle());

    assert_eq!(canonicalize(&expression).unwrap(), scalar);
}

#[test]
fn partial_antisymmetric_sign_stays_inside_nonlinear_tensor() {
    let rep = test_initialize().mink4;
    let a = slot!(rep, signed_nested_a);
    let b = slot!(rep, signed_nested_b);
    let outer = tensor_symbol!(signed_nonlinear_outer);
    let expression = FunctionBuilder::new(outer).add_arg(antisym!(b, a)).finish();

    let canonical = canonicalize(&expression).unwrap();
    assert!(
        matches!(canonical.as_view(), AtomView::Fun(function) if function.get_symbol() == outer)
    );
    assert_eq!(
        canonical,
        FunctionBuilder::new(outer)
            .add_arg(-antisym!(a, b))
            .finish()
    );
}

#[test]
fn bracketed_parameters_remain_opaque() {
    let rep = test_initialize().mink4;
    let a = slot!(rep, signed_bracket_a);
    let b = slot!(rep, signed_bracket_b);
    let outer = tensor_symbol!(signed_linear_outer; Linear);
    let parameter = bracket!(Atom::var(symbol!("p")) + Atom::var(symbol!("q")));
    let expression = FunctionBuilder::new(outer)
        .add_arg(parameter)
        .add_arg(antisym!(a, b))
        .finish();

    assert_eq!(canonicalize(&expression).unwrap(), expression);
}

#[test]
fn mixed_intrinsic_symmetry_returns_a_typed_error() {
    let rep = test_initialize().mink4;
    let slot = slot!(rep, signed_invalid_slot);
    let symmetric = tensor_symbol!(signed_invalid_tensor; Symmetric);
    let expression = function!(symmetric, Atom::num(1), slot.to_atom());

    assert!(matches!(
        canonicalize(&expression),
        Err(CanonicalizationError::InvalidIntrinsicArgument {
            head,
            argument: 0,
            ..
        }) if head == symmetric
    ));
}

#[test]
fn scalar_network_leaves_use_intrinsic_validation() {
    let rep = test_initialize().mink4;
    let slot = slot!(rep, signed_invalid_scalar_slot);
    let symmetric = tensor_symbol!(signed_invalid_scalar_tensor; Symmetric);
    let expression = function!(symmetric, Atom::num(1), slot.to_atom());
    let network = SymbolicNet::<AbstractIndex>::from_scalar(expression);

    assert!(matches!(
        network.canonize(AbstractIndex::Dummy),
        Err(CanonicalizationError::InvalidIntrinsicArgument {
            head,
            argument: 0,
            ..
        }) if head == symmetric
    ));

    let network = SymbolicNet::<AbstractIndex>::from_scalar(Atom::num(1)).fun(symmetric);
    assert!(matches!(
        network.canonize(AbstractIndex::Dummy),
        Err(CanonicalizationError::InvalidIntrinsicArgument {
            head,
            argument: 0,
            ..
        }) if head == symmetric
    ));
}

#[test]
fn compact_metric_shorthand_is_validated_after_normalization() {
    let rep = test_initialize().mink4;
    let expression = function!(
        ETS.metric,
        vector!(signed_metric_left, rep.to_symbolic([])),
        vector!(signed_metric_right, rep.to_symbolic([]))
    );

    let canonical = canonicalize(&expression).unwrap();
    assert_eq!(canonicalize(&canonical).unwrap(), canonical);
}

#[test]
fn canonicalization_is_idempotent() {
    let once = canonicalize(&odd_cycle()).unwrap();
    let twice = canonicalize(&once).unwrap();
    assert_eq!(once, twice);
}

#[test]
fn odd_antisymmetric_cycle_canonicalizes_to_zero() {
    test_initialize();
    assert!(canonicalize(&structural_odd_cycle()).unwrap().is_zero());
}

#[test]
fn nested_odd_summand_is_removed() {
    test_initialize();
    let scalar = tensor!(signed_nested_scalar);
    let expression = &scalar * (Atom::num(1) + structural_odd_cycle());

    assert_eq!(canonicalize(&expression).unwrap(), scalar);
}

#[test]
fn factored_expression_is_preserved_when_no_term_is_removed() {
    test_initialize();
    let a = mink!(4, signed_factored_a);
    let b = mink!(4, signed_factored_b);
    let scalar = tensor!(signed_factored_scalar);
    let expression = scalar * (Atom::num(1) + antisym!(a, b).pow(2));
    let canonical = canonicalize(&expression).unwrap();

    assert!(!canonical.is_zero());
    assert!(matches!(canonical.as_view(), AtomView::Mul(product)
        if product.iter().any(|factor| matches!(factor, AtomView::Add(_)))));
}

#[test]
fn zero_summand_inside_function_is_removed() {
    test_initialize();
    let expression = tensor!(signed_opaque, Atom::num(1) + structural_odd_cycle());
    let expected = tensor!(signed_opaque, Atom::num(1));

    assert_eq!(
        canonicalize(&expression).unwrap(),
        canonicalize(&expected).unwrap()
    );
}

#[test]
fn nonlinear_function_retains_a_zero_argument() {
    test_initialize();
    let expression = tensor!(signed_zero_argument, structural_odd_cycle());
    let expected = tensor!(signed_zero_argument, Atom::Zero);
    let canonical = canonicalize(&expression).unwrap();

    assert!(!canonical.is_zero());
    assert_eq!(canonical, canonicalize(&expected).unwrap());
}

#[test]
fn broadcast_boundary_does_not_hide_independent_zero() {
    test_initialize();
    let a = mink!(4, signed_opaque_a);
    let b = mink!(4, signed_opaque_b);
    let c = mink!(4, signed_opaque_c);
    let inner = antisym!(a.clone(), b.clone()) * antisym!(b.clone(), c.clone()) * antisym!(c, a);
    let broadcast = broadcast_symbol!(signed_opaque_broadcast);
    let opaque = function!(broadcast, Atom::num(1) + &inner);

    // A nonlinear function retains its boundary, while ordinary Sum execution
    // removes a vanishing addend inside its argument.
    assert_eq!(
        canonicalize(&opaque).unwrap(),
        canonicalize(&function!(broadcast, Atom::num(1))).unwrap()
    );

    // A zero whole argument does not erase the nonlinear function itself.
    let zero_argument = canonicalize(&function!(broadcast, inner)).unwrap();
    assert!(!zero_argument.is_zero());
    assert_eq!(
        zero_argument,
        canonicalize(&function!(broadcast, Atom::Zero)).unwrap()
    );

    // The opaque function boundary must not suppress an independent zero.
    assert!(
        canonicalize(&(structural_odd_cycle() * opaque))
            .unwrap()
            .is_zero()
    );
}

#[test]
fn even_antisymmetric_automorphism_survives_canonicalization() {
    test_initialize();
    let a = mink!(4, signed_square_a);
    let b = mink!(4, signed_square_b);

    assert!(!canonicalize(&antisym!(a, b).pow(2)).unwrap().is_zero());
}

#[test]
fn power_copies_closed_components_as_a_whole() {
    test_initialize();

    // Keep the closed contraction as one Power child until Spenso parses it.
    assert!(
        canonicalize(&bracket!(structural_odd_cycle()).pow(2))
            .unwrap()
            .is_zero()
    );
}

#[test]
fn external_indices_prevent_slot_automorphism() {
    test_initialize();
    let a = mink!(4, signed_external_a);
    let b = mink!(4, signed_external_b);

    assert!(!canonicalize(&antisym!(a, b)).unwrap().is_zero());
}

#[test]
fn ordered_slots_prevent_antisymmetric_automorphism() {
    test_initialize();
    let a = mink!(4, signed_ordered_a);
    let b = mink!(4, signed_ordered_b);
    let contraction = tensor!(signed_ordered, a.clone(), b.clone()) * antisym!(a, b);

    assert!(!canonicalize(&contraction).unwrap().is_zero());
}

#[test]
fn symmetric_slots_allow_antisymmetric_automorphism() {
    test_initialize();
    let a = mink!(4, signed_symmetric_a);
    let b = mink!(4, signed_symmetric_b);
    let contraction = sym!(a.clone(), b.clone()) * antisym!(a, b);

    assert!(canonicalize(&contraction).unwrap().is_zero());
}

#[test]
fn different_representation_groups_break_odd_automorphism() {
    test_initialize();
    let a = mink!(4, signed_group_a);
    let b = euc!(4, signed_group_b);
    let c = lor!(4, signed_group_c);
    let contraction =
        antisym!(a.clone(), b.clone()) * antisym!(b, c.clone()) * antisym!(dind!(c), a);

    assert!(!canonicalize(&contraction).unwrap().is_zero());
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
    assert!(!canonicalize(&rank_three).unwrap().is_zero());

    let d = mink!(4, signed_cyclic_d);
    let rank_four = cyclic!(a.clone(), b.clone(), c.clone(), d.clone()) * antisym!(a, b, c, d);

    // A four-cycle is odd. Reflections are not part of cyclic symmetry.
    assert!(canonicalize(&rank_four).unwrap().is_zero());
}
