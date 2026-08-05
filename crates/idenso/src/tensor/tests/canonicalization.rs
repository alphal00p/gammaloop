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

use crate::{IndexTooling, tensor::CanonicalizationError, test_support::test_initialize};

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
fn signed_negative_powers_are_canonical_and_idempotent() {
    let rep = test_initialize().mink4;
    let a = slot!(rep, signed_negative_power_a);
    let b = slot!(rep, signed_negative_power_b);
    let symmetric = tensor_symbol!(signed_negative_power_tensor; Symmetric);
    for power in [-4, -2] {
        let expression = function!(symmetric, a.to_atom(), b.to_atom()).pow(power);
        let once = canonicalize(&expression).unwrap();
        assert_eq!(canonicalize(&once).unwrap(), once);
    }

    let reciprocal_zero = canonicalize(&bracket!(odd_cycle()).pow(-1)).unwrap();
    assert!(!reciprocal_zero.is_zero());
    assert_eq!(canonicalize(&reciprocal_zero).unwrap(), reciprocal_zero);
}

#[test]
fn product_powers_follow_execution_without_losing_sign_or_division_polarity() {
    let rep = test_initialize().mink4;
    let index: LibrarySlot<_> = slot!(rep, signed_product_power_index).cast();
    let left = tensor_symbol!(signed_product_power_left);
    let right = tensor_symbol!(signed_product_power_right);
    let base = function!(left, index.to_atom()) * function!(right, index.to_atom());

    for exponent in [2, -2] {
        let expression = bracket!(base.clone()).pow(exponent);
        let canonical = canonicalize(&expression).unwrap();
        assert_eq!(canonicalize(&canonical).unwrap(), canonical);
    }

    let squared = canonicalize(&bracket!(base.clone()).pow(2)).unwrap();
    let squared_negative_base = canonicalize(&bracket!(-base.clone()).pow(2)).unwrap();
    let negative_squared = canonicalize(&(-bracket!(base).pow(2))).unwrap();
    assert_eq!(squared_negative_base, squared);
    assert_ne!(negative_squared, squared);
}

#[test]
fn signed_execution_matches_symbolica_for_powers_sums_and_opaque_controls() {
    let rep = test_initialize().mink4;
    let first: LibrarySlot<_> = slot!(rep, signed_execution_first).cast();
    let second: LibrarySlot<_> = slot!(rep, signed_execution_second).cast();
    let left = tensor_symbol!(signed_execution_left);
    let right = tensor_symbol!(signed_execution_right);
    let branch = |index: LibrarySlot<AbstractIndex>| {
        function!(left, index.to_atom()) * function!(right, index.to_atom())
    };
    let first_branch = branch(first);
    let second_branch = branch(second);
    let canonical_branch = canonicalize(&first_branch).unwrap();
    assert_eq!(canonicalize(&second_branch).unwrap(), canonical_branch);

    // The bracket is the parser's grouping syntax. Positive powers cover odd and
    // even sign parity; the negative even power covers reciprocal polarity.
    // Negative odd powers of open tensor expressions are deliberately rejected.
    for exponent in [2, 3, -2] {
        let positive = canonicalize(&bracket!(first_branch.clone()).pow(exponent)).unwrap();
        let negative = canonicalize(&bracket!(-first_branch.clone()).pow(exponent)).unwrap();
        let expected = if exponent % 2 == 0 {
            positive.clone()
        } else {
            -positive.clone()
        };
        assert_eq!(negative, expected);
    }

    // Reusing one exact Atom would make Symbolica reduce 2*A + 3*A and A - A
    // before the canonicalizer is called. Alpha-equivalent dummy spellings keep
    // those sums intact until canonical reconstruction aligns both A branches.
    let numeric_sum = Atom::num(2) * first_branch.clone() + Atom::num(3) * second_branch.clone();
    let symbolica_numeric =
        Atom::num(2) * canonical_branch.clone() + Atom::num(3) * canonical_branch.clone();
    assert_ne!(numeric_sum, symbolica_numeric);
    assert_eq!(canonicalize(&numeric_sum).unwrap(), symbolica_numeric);

    let cancellation = first_branch.clone() - second_branch.clone();
    assert!(!cancellation.is_zero());
    assert!(canonicalize(&cancellation).unwrap().is_zero());

    let x = Atom::var(symbol!("signed_execution_x"));
    let y = Atom::var(symbol!("signed_execution_y"));
    let symbolic_sum = &x * &first_branch + &y * &second_branch;
    let symbolica_symbolic = &x * &canonical_branch + &y * &canonical_branch;
    assert_ne!(symbolic_sum, symbolica_symbolic);
    assert_eq!(canonicalize(&symbolic_sum).unwrap(), symbolica_symbolic);

    let symbolic_cancellation = &x * &first_branch - &x * &second_branch;
    assert!(!symbolic_cancellation.is_zero());
    assert!(canonicalize(&symbolic_cancellation).unwrap().is_zero());

    let opaque_index = slot!(rep, signed_execution_opaque_index);
    let parameter = Atom::var(symbol!("signed_execution_opaque_parameter"));
    let nested = function!(
        symbol!("signed_execution_opaque_nested"),
        Atom::var(symbol!("signed_execution_opaque_nested_argument"))
    );
    let untagged = symbol!("signed_execution_untagged_antisymmetric"; Antisymmetric);
    let opaque = function!(
        untagged,
        parameter.clone(),
        opaque_index.to_atom(),
        nested.clone()
    );
    let odd_permutation = function!(untagged, opaque_index.to_atom(), parameter, nested);
    assert_eq!(opaque, -odd_permutation.clone());
    assert_eq!(canonicalize(&opaque).unwrap(), opaque);
    assert_eq!(canonicalize(&odd_permutation).unwrap(), odd_permutation);
}

#[test]
fn outer_zero_and_zero_denominator_are_normalized_by_execution() {
    let expression = odd_cycle() * bracket!(structural_odd_cycle()).pow(-1);
    let canonical = canonicalize(&expression).unwrap();

    assert!(!canonical.is_zero());
    assert_eq!(canonicalize(&canonical).unwrap(), canonical);
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
fn signed_branch_exchange_zeroes_the_enclosing_product() {
    let rep = test_initialize().mink4;
    let i: LibrarySlot<_> = slot!(rep, signed_exchange_i).cast();
    let j: LibrarySlot<_> = slot!(rep, signed_exchange_j).cast();
    let antisymmetric = tensor_symbol!(signed_exchange_antisymmetric; Antisymmetric);
    let left = tensor_symbol!(signed_exchange_left);
    let right = tensor_symbol!(signed_exchange_right);
    let branch = |first, second| function!(left, first) * function!(right, second);
    let expression = function!(antisymmetric, i.to_atom(), j.to_atom())
        * (branch(i.to_atom(), j.to_atom()) + branch(j.to_atom(), i.to_atom()));

    assert!(canonicalize(&expression).unwrap().is_zero());
}

#[test]
fn visible_scalars_receive_canonical_tensor_signs_through_execution() {
    let rep = test_initialize().mink4;
    let i: LibrarySlot<_> = slot!(rep, signed_scalar_sink_i).cast();
    let j: LibrarySlot<_> = slot!(rep, signed_scalar_sink_j).cast();
    let antisymmetric = tensor_symbol!(signed_scalar_sink_tensor; Antisymmetric);
    let ordered = function!(antisymmetric, i.to_atom(), j.to_atom());
    let reversed = Atom::num(2) * function!(antisymmetric, j.to_atom(), i.to_atom());
    let expected = Atom::num(-2) * ordered;

    assert_eq!(
        canonicalize(&reversed).unwrap(),
        canonicalize(&expected).unwrap()
    );
}

#[test]
fn signed_sum_cancellation_is_owned_by_network_execution() {
    let rep = test_initialize().mink4;
    let i: LibrarySlot<_> = slot!(rep, signed_cancellation_i).cast();
    let j: LibrarySlot<_> = slot!(rep, signed_cancellation_j).cast();
    let antisymmetric = tensor_symbol!(signed_cancellation_tensor; Antisymmetric);
    let forward = Atom::num(2) * function!(antisymmetric, i.to_atom(), j.to_atom());
    let reverse = Atom::num(2) * function!(antisymmetric, j.to_atom(), i.to_atom());

    assert!(canonicalize(&(forward + reverse)).unwrap().is_zero());
}

#[test]
fn common_signed_sum_branches_combine_after_leaf_local_sign_reconstruction() {
    let rep = test_initialize().mink4;
    let i: LibrarySlot<_> = slot!(rep, AbstractIndex::Dummy(0)).cast();
    let j: LibrarySlot<_> = slot!(rep, AbstractIndex::Dummy(1)).cast();
    let k: LibrarySlot<_> = slot!(rep, AbstractIndex::Dummy(2)).cast();
    let l: LibrarySlot<_> = slot!(rep, AbstractIndex::Dummy(3)).cast();
    let antisymmetric = tensor_symbol!(signed_common_sum_tensor; Antisymmetric);
    let left = tensor_symbol!(signed_common_sum_left);
    let right = tensor_symbol!(signed_common_sum_right);
    let scalar = Atom::var(symbol!("signed_common_sum_scalar"));
    let branch = |first: LibrarySlot<_>, second: LibrarySlot<_>| {
        function!(antisymmetric, first.to_atom(), second.to_atom())
            * function!(left, second.to_atom())
            * function!(right, first.to_atom())
    };
    let expression = scalar.clone() * (branch(i, j) + branch(k, l));
    let expected = Atom::num(-2)
        * scalar
        * function!(antisymmetric, i.to_atom(), j.to_atom())
        * function!(left, i.to_atom())
        * function!(right, j.to_atom());

    let canonical = canonicalize(&expression).unwrap();
    assert_eq!(canonical, expected);
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
fn partial_antisymmetric_sign_lifts_through_an_explicitly_linear_tensor() {
    let rep = test_initialize().mink4;
    let a = slot!(rep, signed_linear_nested_a);
    let b = slot!(rep, signed_linear_nested_b);
    let linear = tensor_symbol!(signed_linear_nested_outer; Linear);
    let expression = FunctionBuilder::new(linear)
        .add_arg(antisym!(b, a))
        .finish();
    let expected = -FunctionBuilder::new(linear)
        .add_arg(antisym!(a, b))
        .finish();

    let canonical = canonicalize(&expression).unwrap();
    assert_eq!(canonical, expected);
    assert_eq!(canonicalize(&canonical).unwrap(), canonical);
}

#[test]
fn unbracketed_parameter_sum_uses_declared_linear_distribution() {
    let rep = test_initialize().mink4;
    let a = slot!(rep, signed_linear_sum_a);
    let b = slot!(rep, signed_linear_sum_b);
    let linear = tensor_symbol!(signed_linear_sum_outer; Linear);
    let p = Atom::var(symbol!("signed_linear_sum_p"));
    let q = Atom::var(symbol!("signed_linear_sum_q"));
    let group = antisym!(a, b);
    let expression = FunctionBuilder::new(linear)
        .add_arg(&p + &q)
        .add_arg(group.clone())
        .finish();
    let distributed = FunctionBuilder::new(linear)
        .add_arg(p)
        .add_arg(group.clone())
        .finish()
        + FunctionBuilder::new(linear)
            .add_arg(q)
            .add_arg(group)
            .finish();

    assert_eq!(
        canonicalize(&expression).unwrap(),
        canonicalize(&distributed).unwrap()
    );
}

#[test]
fn nested_signed_zero_respects_declared_function_linearity() {
    let nonlinear = tensor_symbol!(signed_nested_zero_nonlinear);
    let linear = tensor_symbol!(signed_nested_zero_linear; Linear);
    let zero_scope = odd_cycle();
    let nonlinear_expression = function!(nonlinear, zero_scope.clone());
    let linear_expression = function!(linear, zero_scope);
    let expected_nonlinear = function!(nonlinear, Atom::Zero);

    let nonlinear_canonical = canonicalize(&nonlinear_expression).unwrap();
    assert_eq!(
        nonlinear_canonical,
        canonicalize(&expected_nonlinear).unwrap()
    );
    assert_eq!(
        canonicalize(&nonlinear_canonical).unwrap(),
        nonlinear_canonical
    );
    assert!(canonicalize(&linear_expression).unwrap().is_zero());
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
