use spenso::{
    slot,
    structure::{abstract_index::AbstractIndex, representation::LibrarySlot, slot::IsAbstractSlot},
    sym, tensor_symbol,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    function, symbol,
};

use crate::{IndexTooling, test_support::test_initialize};

#[test]
fn compound_structured_indices_round_trip_through_canonicalization() {
    let rep = test_initialize().mink4;
    let tensor = tensor_symbol!(unified_structured_index_tensor);

    for payload in [
        Atom::var(symbol!("unified_structured_i")) + Atom::var(symbol!("unified_structured_j")),
        Atom::var(symbol!("unified_structured_i")) * Atom::var(symbol!("unified_structured_j")),
        Atom::var(symbol!("unified_structured_i")).pow(Atom::num(2)),
    ] {
        let expression = function!(tensor, rep.to_symbolic([payload]));
        assert_eq!(
            expression
                .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .unwrap(),
            expression,
        );
    }
}

#[test]
fn ordered_scalar_arguments_keep_their_positions_around_slots() {
    let rep = test_initialize().mink4;
    let first: LibrarySlot<_> = slot!(rep, unified_ordered_argument_first).cast();
    let second: LibrarySlot<_> = slot!(rep, unified_ordered_argument_second).cast();
    let outer = tensor_symbol!(unified_ordered_argument_outer);
    let left = tensor_symbol!(unified_ordered_argument_left);
    let right = tensor_symbol!(unified_ordered_argument_right);
    let before = Atom::var(symbol!("unified_ordered_argument_before"));
    let between = Atom::var(symbol!("unified_ordered_argument_between"));
    let after = Atom::var(symbol!("unified_ordered_argument_after"));
    let expression = FunctionBuilder::new(outer)
        .add_arg(before.clone())
        .add_arg(first.to_atom())
        .add_arg(between.clone())
        .add_arg(second.to_atom())
        .add_arg(after.clone())
        .finish()
        * function!(left, first.to_atom())
        * function!(right, second.to_atom());

    let canonical = expression
        .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .unwrap();
    let AtomView::Mul(product) = canonical.as_view() else {
        panic!("the closed contraction must remain a Product")
    };
    let tensor = product
        .iter()
        .find_map(|factor| match factor {
            AtomView::Fun(function) if function.get_symbol() == outer => Some(function),
            _ => None,
        })
        .expect("the ordered tensor must remain visible");
    let arguments = tensor.iter().collect::<Vec<_>>();
    assert_eq!(arguments[0], before.as_view());
    assert_eq!(arguments[2], between.as_view());
    assert_eq!(arguments[4], after.as_view());
}

#[test]
fn free_index_renaming_is_equivariant_while_dummies_are_canonical() {
    let rep = test_initialize().mink4;
    let first_external: LibrarySlot<_> = slot!(rep, unified_free_external_first).cast();
    let second_external: LibrarySlot<_> = slot!(rep, unified_free_external_second).cast();
    let dummy: LibrarySlot<_> = slot!(rep, unified_free_dummy).cast();
    let outer = tensor_symbol!(unified_free_outer);
    let leaf = tensor_symbol!(unified_free_leaf);
    let build = |external: LibrarySlot<_>| {
        function!(outer, external.to_atom(), dummy.to_atom()) * function!(leaf, dummy.to_atom())
    };

    let first = build(first_external)
        .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .unwrap();
    let second = build(second_external)
        .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .unwrap();
    assert_eq!(
        first
            .replace(first_external.to_atom())
            .with(second_external.to_atom()),
        second
    );
}

#[test]
fn association_and_commutative_order_preserve_the_canonical_result() {
    let rep = test_initialize().mink4;
    let path_left = tensor_symbol!(unified_order_path_left);
    let path_middle = tensor_symbol!(unified_order_path_middle);
    let path_right = tensor_symbol!(unified_order_path_right);
    let edge_left = tensor_symbol!(unified_order_edge_left);
    let edge_right = tensor_symbol!(unified_order_edge_right);

    let i: LibrarySlot<_> = slot!(rep, unified_order_i).cast();
    let j: LibrarySlot<_> = slot!(rep, unified_order_j).cast();
    let k: LibrarySlot<_> = slot!(rep, unified_order_k).cast();
    let first_path = (function!(path_left, i.to_atom())
        * function!(path_middle, i.to_atom(), j.to_atom()))
        * function!(path_right, j.to_atom());
    let first_edge = function!(edge_left, k.to_atom()) * function!(edge_right, k.to_atom());
    let first = (first_path + first_edge) + Atom::num(7);

    let x: LibrarySlot<_> = slot!(rep, unified_order_x).cast();
    let y: LibrarySlot<_> = slot!(rep, unified_order_y).cast();
    let z: LibrarySlot<_> = slot!(rep, unified_order_z).cast();
    let reordered_path = function!(path_right, y.to_atom())
        * (function!(path_middle, x.to_atom(), y.to_atom()) * function!(path_left, x.to_atom()));
    let reordered_edge = function!(edge_right, z.to_atom()) * function!(edge_left, z.to_atom());
    let reordered = Atom::num(7) + (reordered_edge + reordered_path);

    assert_eq!(
        first
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap(),
        reordered
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap()
    );
}

#[test]
fn shared_sum_interfaces_use_canonical_line_order_for_dummies() {
    let rep = test_initialize().mink4;
    let outer = tensor_symbol!(unified_sum_interface_outer);
    let first_left = tensor_symbol!(unified_sum_interface_first_left);
    let first_right = tensor_symbol!(unified_sum_interface_first_right);
    let second_left = tensor_symbol!(unified_sum_interface_second_left);
    let second_right = tensor_symbol!(unified_sum_interface_second_right);
    let build = |left, right| {
        let left = rep.slot::<AbstractIndex, _>(left).to_atom();
        let right = rep.slot::<AbstractIndex, _>(right).to_atom();
        function!(outer, left.clone())
            * function!(outer, right.clone())
            * (function!(first_left, left.clone()) * function!(first_right, right.clone())
                + function!(second_left, left) * function!(second_right, right))
    };

    // Identical outer heads reverse source traversal by dummy name, while the
    // distinct branch heads give the two lines fixed canonical roles.
    let forward = build(AbstractIndex::Dummy(0), AbstractIndex::Dummy(1));
    let reversed = build(AbstractIndex::Dummy(1), AbstractIndex::Dummy(0));
    assert_eq!(
        forward
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap(),
        reversed
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap(),
    );
}

#[test]
fn identical_local_fragments_retain_different_global_wiring() {
    let rep = test_initialize().mink4;
    let tensor = tensor_symbol!(unified_wiring_tensor);
    let a: LibrarySlot<_> = slot!(rep, unified_wiring_a).cast();
    let b: LibrarySlot<_> = slot!(rep, unified_wiring_b).cast();
    let c: LibrarySlot<_> = slot!(rep, unified_wiring_c).cast();
    let d: LibrarySlot<_> = slot!(rep, unified_wiring_d).cast();

    let cycle = function!(tensor, a.to_atom(), b.to_atom())
        * function!(tensor, b.to_atom(), c.to_atom())
        * function!(tensor, c.to_atom(), d.to_atom())
        * function!(tensor, d.to_atom(), a.to_atom());
    let double_edges = function!(tensor, a.to_atom(), b.to_atom())
        * function!(tensor, b.to_atom(), a.to_atom())
        * function!(tensor, c.to_atom(), d.to_atom())
        * function!(tensor, d.to_atom(), c.to_atom());

    let canonical_cycle = cycle
        .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .unwrap();
    let canonical_double_edges = double_edges
        .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .unwrap();

    assert_ne!(canonical_cycle, canonical_double_edges);
    assert_eq!(
        canonical_cycle
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap(),
        canonical_cycle
    );
    assert_eq!(
        canonical_double_edges
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap(),
        canonical_double_edges
    );
}

#[test]
fn partial_groups_keep_distinct_stable_layout_paths() {
    let rep = test_initialize().mink4;
    let outer = tensor_symbol!(unified_layout_outer);
    let first_leaf = tensor_symbol!(unified_layout_first_leaf);
    let second_leaf = tensor_symbol!(unified_layout_second_leaf);
    let third_leaf = tensor_symbol!(unified_layout_third_leaf);
    let fourth_leaf = tensor_symbol!(unified_layout_fourth_leaf);
    let parameter = Atom::var(symbol!("unified_layout_parameter"));

    let build = |i: LibrarySlot<_>,
                 j: LibrarySlot<_>,
                 k: LibrarySlot<_>,
                 l: LibrarySlot<_>,
                 reverse_members: bool,
                 swap_groups: bool| {
        let first_group = if reverse_members {
            sym!(j, i)
        } else {
            sym!(i, j)
        };
        let second_group = if reverse_members {
            sym!(l, k)
        } else {
            sym!(k, l)
        };
        let (first_group, second_group) = if swap_groups {
            (second_group, first_group)
        } else {
            (first_group, second_group)
        };
        FunctionBuilder::new(outer)
            .add_arg(first_group)
            .add_arg(parameter.clone())
            .add_arg(second_group)
            .finish()
            * function!(first_leaf, i.to_atom())
            * function!(second_leaf, j.to_atom())
            * function!(third_leaf, k.to_atom())
            * function!(fourth_leaf, l.to_atom())
    };

    let original = build(
        slot!(rep, unified_layout_a).cast(),
        slot!(rep, unified_layout_b).cast(),
        slot!(rep, unified_layout_c).cast(),
        slot!(rep, unified_layout_d).cast(),
        false,
        false,
    );
    let member_permuted = build(
        slot!(rep, unified_layout_w).cast(),
        slot!(rep, unified_layout_x).cast(),
        slot!(rep, unified_layout_y).cast(),
        slot!(rep, unified_layout_z).cast(),
        true,
        false,
    );
    let group_swapped = build(
        slot!(rep, unified_layout_e).cast(),
        slot!(rep, unified_layout_f).cast(),
        slot!(rep, unified_layout_g).cast(),
        slot!(rep, unified_layout_h).cast(),
        false,
        true,
    );

    let canonical = original
        .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .unwrap();
    assert_eq!(
        member_permuted
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap(),
        canonical
    );
    assert_ne!(
        group_swapped
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap(),
        canonical
    );
}

#[test]
fn repeated_partial_group_slots_remain_distinct_occurrences() {
    let rep = test_initialize().mink4;
    let outer = tensor_symbol!(unified_repeated_group_outer);
    let first: LibrarySlot<_> = slot!(rep, unified_repeated_group_first).cast();
    let second: LibrarySlot<_> = slot!(rep, unified_repeated_group_second).cast();
    let expression = FunctionBuilder::new(outer)
        .add_arg(sym!(first, first))
        .finish();
    let renamed = FunctionBuilder::new(outer)
        .add_arg(sym!(second, second))
        .finish();

    let canonical = expression
        .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
        .unwrap();
    assert_eq!(
        renamed
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap(),
        canonical
    );

    let AtomView::Fun(tensor) = canonical.as_view() else {
        panic!("the repeated partial group must remain inside its tensor")
    };
    let Some(AtomView::Fun(group)) = tensor.iter().next() else {
        panic!("the partial group must be reconstructed as a function")
    };
    let members = group.iter().collect::<Vec<_>>();
    assert_eq!(members.len(), 2);
    assert_eq!(members[0], members[1]);
}

#[test]
fn untagged_symbolica_symmetries_remain_opaque_controls() {
    let rep = test_initialize().mink4;
    let index: LibrarySlot<_> = slot!(rep, unified_untagged_index).cast();
    let parameter = Atom::var(symbol!("unified_untagged_parameter"));
    let nested = function!(
        symbol!("unified_untagged_nested"),
        Atom::var(symbol!("unified_untagged_nested_argument"))
    );
    let symmetric = symbol!("unified_untagged_symmetric"; Symmetric);
    let cyclesymmetric = symbol!("unified_untagged_cyclic"; Cyclesymmetric);

    let symmetric_value = function!(
        symmetric,
        parameter.clone(),
        index.to_atom(),
        nested.clone()
    );
    let symmetric_permuted = function!(
        symmetric,
        nested.clone(),
        parameter.clone(),
        index.to_atom()
    );
    assert_eq!(symmetric_value, symmetric_permuted);
    assert_eq!(
        symmetric_value
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap(),
        symmetric_value
    );

    let cyclic_value = function!(
        cyclesymmetric,
        parameter.clone(),
        index.to_atom(),
        nested.clone()
    );
    let cyclic_rotation = function!(
        cyclesymmetric,
        nested.clone(),
        parameter.clone(),
        index.to_atom()
    );
    let cyclic_reflection = function!(cyclesymmetric, parameter, nested, index.to_atom());
    assert_eq!(cyclic_value, cyclic_rotation);
    assert_ne!(cyclic_value, cyclic_reflection);
    assert_eq!(
        cyclic_value
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap(),
        cyclic_value
    );
    assert_eq!(
        cyclic_reflection
            .try_canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .unwrap(),
        cyclic_reflection
    );
}
