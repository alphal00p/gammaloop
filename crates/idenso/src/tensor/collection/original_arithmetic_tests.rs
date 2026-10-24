//! Source-only controls for original arithmetic cuts. This slice is UNRUN.

use super::{OriginalArithmeticAction as Action, OriginalArithmeticOperation as Operation};
use crate::tensor::SymbolicTensor;
use spenso::{
    network::{
        library::{function_lib::INBUILTS, symbolic::ExplicitKey},
        tags::SPENSO_TAG,
    },
    structure::{
        abstract_index::AbstractIndex,
        partial::{PartialStructure, PartialStructureExt},
        representation::{LibraryRep, Minkowski, RepName, Representation},
        slot::IsAbstractSlot,
    },
};
use std::{
    collections::HashMap,
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    function, symbol,
};

type Tensor = SymbolicTensor<PartialStructure>;

fn initialize() {
    crate::test_support::test_initialize();
}

fn scalar(atom: Atom) -> Tensor {
    Tensor::from_declared_parts(atom, PartialStructure::from_logical_slots([])).unwrap()
}

fn leaf(name: &str, arguments: &[Atom]) -> Tensor {
    Tensor::infer(
        FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
            .add_args(arguments)
            .finish(),
    )
    .unwrap()
}

fn oriented_slot(rep: Representation<LibraryRep>, index: usize) -> Atom {
    rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(index))
        .to_atom()
}

fn slot(index: usize) -> Atom {
    oriented_slot(Minkowski {}.new_rep(4).cast::<LibraryRep>(), index)
}

/// Test callback replay uses only the existing typed arithmetic owners.
fn arithmetic(_: &Tensor, op: Operation, children: &[Tensor]) -> super::Result<Tensor> {
    let mut values = children.iter();
    let first = values
        .next()
        .expect("an original Add/Mul is nonempty")
        .clone();
    values.try_fold(first, |left, right| match op {
        Operation::Sum => left.try_add(right),
        Operation::Product => Ok(left.multiply(right)?),
    })
}

#[test]
fn complete_root_cuts_run_once_and_do_not_visit_replacement_payloads() {
    initialize();
    let a = Atom::var(symbol!("original_arithmetic::root_a"));
    let b = Atom::var(symbol!("original_arithmetic::root_b"));
    let source = scalar((&a + &b) * (&a - &b));
    let replacement = scalar(function!(
        symbol!("original_arithmetic::generated_root"; Scalar),
        &source.expression,
        a.pow(2)
    ));
    for replace in [false, true] {
        let mut calls = 0;
        let actual = source
            .map_original_arithmetic(
                &mut |original| {
                    calls += 1;
                    assert!(original.same_value(&source));
                    Ok(if replace {
                        Action::Replace(replacement.clone())
                    } else {
                        Action::Retain
                    })
                },
                &mut |_, _, _| panic!("a complete source cut has no emitted operation"),
            )
            .unwrap();
        assert_eq!(calls, 1);
        assert!(actual.same_value(if replace { &replacement } else { &source }));
    }
}

#[test]
fn nested_depth_one_opening_uses_exact_stored_originals_and_respects_child_cuts() {
    initialize();
    let a = Atom::var(symbol!("original_arithmetic::nested_a"));
    let b = Atom::var(symbol!("original_arithmetic::nested_b"));
    let c = Atom::var(symbol!("original_arithmetic::nested_c"));
    let d = Atom::var(symbol!("original_arithmetic::nested_d"));
    let left = &a + &b;
    let right = &c + &d;
    let original = &left * &right;
    let source = scalar(original.clone());
    let replacement = scalar(a.pow(2));
    let mut seen = Vec::new();
    let mut operations = Vec::new();
    let actual = source
        .map_original_arithmetic(
            &mut |value| {
                seen.push(value.expression.clone());
                Ok(if value.expression == left {
                    Action::Retain
                } else if value.expression == c {
                    Action::Replace(replacement.clone())
                } else if matches!(
                    value.expression.as_view(),
                    AtomView::Add(_) | AtomView::Mul(_)
                ) {
                    Action::Descend
                } else {
                    Action::Retain
                })
            },
            &mut |value, op, children| {
                operations.push(value.expression.clone());
                if value.expression == right {
                    assert_eq!(op, Operation::Sum);
                    assert_eq!(children.len(), 2);
                } else {
                    assert_eq!(value.expression, original);
                    assert_eq!(op, Operation::Product);
                    assert!(children.iter().any(|child| child.expression == left));
                }
                arithmetic(value, op, children)
            },
        )
        .unwrap();
    assert_eq!(seen.len(), 5); // root, both original sums, original c and d.
    for expected in [&original, &left, &right, &c, &d] {
        assert_eq!(seen.iter().filter(|actual| *actual == expected).count(), 1);
    }
    assert!(!seen.contains(&a)); // Inside the retained original left sum.
    assert!(!seen.contains(&replacement.expression)); // Newly emitted syntax.
    assert_eq!(operations.len(), 2);
    assert!(operations.contains(&right));
    assert!(operations.contains(&original));
    assert!(actual.same_value(&scalar(left * (&replacement.expression + d))));
}

#[test]
fn unchanged_nested_arithmetic_never_calls_the_operation_hook() {
    initialize();
    let a = Atom::var(symbol!("original_arithmetic::unchanged_a"));
    let b = Atom::var(symbol!("original_arithmetic::unchanged_b"));
    let c = Atom::var(symbol!("original_arithmetic::unchanged_c"));
    let source = scalar((&a + &b) * c);
    let mut nodes = 0;
    let actual = source
        .map_original_arithmetic(
            &mut |value| {
                nodes += 1;
                Ok(
                    if matches!(
                        value.expression.as_view(),
                        AtomView::Add(_) | AtomView::Mul(_)
                    ) {
                        Action::Descend
                    } else {
                        Action::Retain
                    },
                )
            },
            &mut |_, _, _| panic!("unchanged operations retain their original carrier"),
        )
        .unwrap();
    assert_eq!(nodes, 5);
    assert!(actual.same_value(&source));
}

#[test]
fn changing_a_sibling_does_not_replay_an_opaque_normalizer_or_visit_its_body() {
    initialize();
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&calls);
    let hooked = spenso::tensor_symbol!(
        "original_arithmetic::opaque_callback",
        norm = move |_, _| {
            observed.fetch_add(1, Ordering::Relaxed);
        }
    );
    let metadata = symbol!("original_arithmetic::opaque_metadata"; Scalar);
    let h = Atom::var(symbol!("original_arithmetic::opaque_H"));
    let q = Atom::var(symbol!("original_arithmetic::opaque_q"));
    let private = function!(metadata, &h, slot(99201));
    let opaque = leaf("original_arithmetic::opaque_carrier", &[private.clone()]);
    let hooked_atom = function!(hooked, private.clone());
    let source = scalar(&h * &hooked_atom * &opaque.expression);
    let replacement = scalar(q.clone());
    calls.store(0, Ordering::Relaxed);
    let mut h_calls = 0;
    let mut opaque_calls = 0;
    let actual = source
        .map_original_arithmetic(
            &mut |value| {
                Ok(if value.expression == source.expression {
                    Action::Descend
                } else if value.expression == h {
                    h_calls += 1;
                    Action::Replace(replacement.clone())
                } else {
                    opaque_calls += 1;
                    assert!(
                        value.expression == hooked_atom || value.expression == opaque.expression
                    );
                    Action::Retain
                })
            },
            &mut arithmetic,
        )
        .unwrap();
    assert_eq!(h_calls, 1);
    assert_eq!(opaque_calls, 2);
    assert_eq!(calls.load(Ordering::Relaxed), 0);
    assert!(actual.expression.contains(private.as_view()));
    assert_eq!(actual.expression, q * hooked_atom * opaque.expression);
    assert_eq!(calls.load(Ordering::Relaxed), 0);
}

#[test]
fn changed_sum_hook_receives_original_root_and_uncombined_typed_operands() {
    initialize();
    let a = Atom::var(symbol!("original_arithmetic::sum_a"));
    let b = Atom::var(symbol!("original_arithmetic::sum_b"));
    let source = scalar(&a + &b);
    let generated = symbol!("original_arithmetic::sum_generated"; Scalar);
    let mut visited = Vec::new();
    let mut operations = 0;
    let actual = source
        .map_original_arithmetic(
            &mut |value| {
                visited.push(value.expression.clone());
                Ok(if value.expression == source.expression {
                    Action::Descend
                } else if value.expression == a {
                    Action::Replace(scalar(b.clone()))
                } else {
                    assert_eq!(value.expression, b);
                    Action::Retain
                })
            },
            &mut |original, op, children| {
                operations += 1;
                assert!(original.same_value(&source));
                assert_eq!(op, Operation::Sum);
                assert_eq!(children.len(), 2);
                assert!(children.iter().all(|child| child.expression == b));
                // The ordinary Add constructor would have combined these to
                // 2*b. The checked hook owns this explicit output boundary.
                Ok(scalar(function!(
                    generated,
                    &children[0].expression,
                    &children[1].expression
                )))
            },
        )
        .unwrap();
    assert_eq!(operations, 1);
    assert_eq!(visited.len(), 3);
    assert_eq!(visited.iter().filter(|value| **value == b).count(), 1);
    assert_eq!(actual.expression, function!(generated, &b, &b));
}

#[test]
fn a_depth_one_power_child_stays_the_complete_original_power() {
    initialize();
    let a = Atom::var(symbol!("original_arithmetic::power_a"));
    let b = Atom::var(symbol!("original_arithmetic::power_b"));
    let c = Atom::var(symbol!("original_arithmetic::power_c"));
    let z = Atom::var(symbol!("original_arithmetic::power_z"));
    let base = &a + &b;
    let power = base.pow(2);
    let source = scalar(&power * &c);
    let mut seen = Vec::new();
    let actual = source
        .map_original_arithmetic(
            &mut |original| {
                seen.push(original.expression.clone());
                Ok(if original.expression == source.expression {
                    Action::Descend
                } else if original.expression == c {
                    Action::Replace(scalar(z.clone()))
                } else {
                    assert_eq!(original.expression, power);
                    Action::Retain
                })
            },
            &mut arithmetic,
        )
        .unwrap();
    assert_eq!(seen.len(), 3);
    assert!(!seen.contains(&base));
    assert_eq!(actual.expression, power * z);
}

#[test]
fn scalar_times_held_word_can_emit_the_native_current_square_without_an_atom_power() {
    initialize();
    let head = symbol!("original_arithmetic::exterior"; Flat, Linear, Antisymmetric;
        tag = &SPENSO_TAG.tensor_word);
    let h = scalar(Atom::var(symbol!("original_arithmetic::word_H")));
    let k = scalar(Atom::var(symbol!("original_arithmetic::word_k")));
    let bar = leaf("original_arithmetic::word_barpsi", &[crate::bis!(4, 99211)]);
    let psi = leaf("original_arithmetic::word_psi", &[crate::bis!(4, 99211)]);
    let current = Tensor::word(head, &[bar, psi]).unwrap();
    // Scope every original/RHS closed contraction before emission. This API
    // deliberately does not invent an operation-wide simultaneous RHS scope.
    let old = current
        .wrap_indices(symbol!("original_arithmetic::old_current"), true)
        .unwrap();
    let image = current
        .wrap_indices(symbol!("original_arithmetic::new_current"), true)
        .unwrap();
    let replacement = image.divide(&k).unwrap();
    let source = h.multiply_word(head, &old).unwrap();
    let expected = image.multiply_word(head, &old).unwrap().divide(&k).unwrap();
    assert!(!expected.is_zero());
    let mut visits = 0;
    let mut operations = 0;
    let actual = source
        .map_original_arithmetic(
            &mut |original| {
                visits += 1;
                Ok(if original.expression == source.expression {
                    Action::Descend
                } else if original.same_value(&h) {
                    Action::Replace(replacement.clone())
                } else {
                    assert!(original.same_value(&old));
                    Action::Retain
                })
            },
            &mut |original, op, children| {
                operations += 1;
                assert!(original.same_value(&source));
                assert_eq!(op, Operation::Product);
                assert_eq!(children.len(), 2);
                Ok(children[0].multiply_word(head, &children[1])?)
            },
        )
        .unwrap();
    assert_eq!(visits, 3); // The two original generators are never traversed.
    assert_eq!(operations, 1);
    assert!(!actual.is_zero());
    assert!(
        actual
            .canonize()
            .unwrap()
            .same_value(&expected.canonize().unwrap())
    );
}

#[test]
fn complete_nonselfdual_interfaces_reject_wrong_duality_dimension_label_and_order() {
    initialize();
    let family = LibraryRep::new_dual("original_arithmetic::interface_flavor").unwrap();
    let rep = family.new_rep(3);
    let other = Minkowski {}.new_rep(4).cast::<LibraryRep>();
    let source = leaf(
        "original_arithmetic::interface_T",
        &[oriented_slot(rep, 99221), oriented_slot(other, 99223)],
    );
    let original = Tensor::from_declared_parts(
        source.expression.clone(),
        PartialStructure::from_logical_slots(source.structure.logical_slots().into_iter().rev()),
    )
    .unwrap();
    let image = Tensor::from_declared_parts(
        leaf(
            "original_arithmetic::interface_U",
            &[oriented_slot(rep, 99221), oriented_slot(other, 99223)],
        )
        .expression,
        original.structure.clone(),
    )
    .unwrap();
    let actual = original
        .map_original_arithmetic(
            &mut |_| Ok(Action::Replace(image.clone())),
            &mut |_, _, _| panic!("a complete nonselfdual replacement is a cut"),
        )
        .unwrap();
    assert!(actual.same_value(&image));
    let reversed = |value: Tensor| {
        Tensor::from_declared_parts(
            value.expression,
            PartialStructure::from_logical_slots(value.structure.logical_slots().into_iter().rev()),
        )
        .unwrap()
    };
    let wrongs = [
        source.clone(), // Same expression, different public logical order.
        reversed(leaf(
            "original_arithmetic::interface_dual",
            &[
                oriented_slot(rep.dual(), 99221),
                oriented_slot(other, 99223),
            ],
        )),
        reversed(leaf(
            "original_arithmetic::interface_dimension",
            &[
                oriented_slot(family.new_rep(5), 99221),
                oriented_slot(other, 99223),
            ],
        )),
        reversed(leaf(
            "original_arithmetic::interface_label",
            &[oriented_slot(rep, 99225), oriented_slot(other, 99223)],
        )),
    ];
    for wrong in wrongs {
        let error = original
            .map_original_arithmetic(
                &mut |_| Ok(Action::Replace(wrong.clone())),
                &mut |_, _, _| panic!("invalid image cannot trigger descent"),
            )
            .unwrap_err();
        assert!(error.to_string().contains("ordered tensor interface"));
    }
}

#[test]
fn nested_nonselfdual_children_keep_their_actual_ports_and_root_publication_is_checked() {
    initialize();
    let rep = LibraryRep::new_dual("original_arithmetic::nested_flavor")
        .unwrap()
        .new_rep(3);
    let value = leaf(
        "original_arithmetic::nested_T",
        &[oriented_slot(rep, 99231), slot(99233)],
    );
    let coefficient = Atom::var(symbol!("original_arithmetic::nested_coefficient"));
    let expression = &coefficient * &value.expression;
    let original = Tensor::from_declared_parts(
        expression,
        PartialStructure::from_logical_slots(value.structure.logical_slots().into_iter().rev()),
    )
    .unwrap();
    let mut visit = |node: &Tensor| {
        Ok(if node.expression == original.expression {
            Action::Descend
        } else if node.expression == coefficient {
            Action::Replace(scalar(Atom::num(2)))
        } else {
            assert!(node.same_value(&value));
            Action::Retain
        })
    };
    // A mathematically valid product with the wrong public port ordering is
    // still an invalid callback result for this exact original boundary.
    assert!(
        original
            .map_original_arithmetic(&mut visit, &mut arithmetic)
            .is_err()
    );
    let actual = original
        .map_original_arithmetic(&mut visit, &mut |source, operation, children| {
            let output = arithmetic(source, operation, children)?;
            Tensor::from_declared_parts(output.expression, source.structure.clone())
        })
        .unwrap();
    let expected =
        Tensor::from_declared_parts(2 * &value.expression, original.structure.clone()).unwrap();
    assert!(actual.same_value(&expected));
    assert_eq!(
        actual.structure.logical_slots(),
        original.structure.logical_slots()
    );
}

#[test]
fn ranked_zero_cuts_keep_exact_logical_ports_and_reject_scalar_or_reordered_zero() {
    initialize();
    let rep = LibraryRep::new_dual("original_arithmetic::zero_flavor")
        .unwrap()
        .new_rep(3);
    let source = leaf(
        "original_arithmetic::zero_T",
        &[oriented_slot(rep, 99241), slot(99243)],
    );
    let source = Tensor::from_declared_parts(
        Atom::Zero,
        PartialStructure::from_logical_slots(source.structure.logical_slots().into_iter().rev()),
    )
    .unwrap();
    let mut calls = 0;
    let actual = source
        .map_original_arithmetic(
            &mut |original| {
                calls += 1;
                assert!(original.same_value(&source));
                Ok(Action::Retain)
            },
            &mut |_, _, _| panic!("a retained zero has no arithmetic operation"),
        )
        .unwrap();
    assert_eq!(calls, 1);
    assert_eq!(actual.rank(), 2);
    assert!(actual.same_value(&source));
    for wrong in [
        scalar(Atom::Zero),
        Tensor::from_declared_parts(
            Atom::Zero,
            PartialStructure::from_logical_slots(
                source.structure.logical_slots().into_iter().rev(),
            ),
        )
        .unwrap(),
    ] {
        assert!(
            source
                .map_original_arithmetic(
                    &mut |_| Ok(Action::Replace(wrong.clone())),
                    &mut |_, _, _| panic!("a zero replacement is a complete cut"),
                )
                .is_err()
        );
    }
    assert!(
        source
            .map_original_arithmetic(
                &mut |_| Ok(Action::Replace(source.clone())),
                &mut arithmetic,
            )
            .unwrap()
            .same_value(&source)
    );
}

#[test]
fn powers_functions_and_metadata_require_explicit_cuts_instead_of_descent() {
    initialize();
    let a = Atom::var(symbol!("original_arithmetic::refused_a"));
    let b = Atom::var(symbol!("original_arithmetic::refused_b"));
    let metadata = symbol!("original_arithmetic::refused_metadata"; Scalar);
    let body = &a + &b;
    let component = leaf("original_arithmetic::refused_conjugate", &[slot(99251)]);
    let conjugated = Tensor::infer(INBUILTS.conj(&component.expression)).unwrap();
    assert!(matches!(conjugated.expression.as_view(), AtomView::Fun(_)));
    for source in [
        scalar(body.pow(2)),
        scalar(function!(metadata, &body)),
        scalar(function!(SPENSO_TAG.pure_scalar, &body)),
        scalar(function!(SPENSO_TAG.bracket, &body)),
        conjugated,
    ] {
        let mut calls = 0;
        let error = source
            .map_original_arithmetic(
                &mut |original| {
                    calls += 1;
                    assert!(original.same_value(&source));
                    Ok(Action::Descend)
                },
                &mut |_, _, _| panic!("unsupported original headers cannot emit an operation"),
            )
            .unwrap_err();
        assert_eq!(calls, 1);
        assert!(error.to_string().contains("explicit semantic cut"));
        assert!(
            source
                .map_original_arithmetic(&mut |_| Ok(Action::Retain), &mut arithmetic,)
                .unwrap()
                .same_value(&source)
        );
    }
}

#[test]
fn changed_native_sum_keeps_a_ranked_zero_and_rejects_a_scalar_zero_result() {
    initialize();
    let rep = LibraryRep::new_dual("original_arithmetic::sum_zero_flavor")
        .unwrap()
        .new_rep(3);
    let ports = [oriented_slot(rep, 99261), slot(99263)];
    let a = leaf("original_arithmetic::sum_zero_A", &ports);
    let b = leaf("original_arithmetic::sum_zero_B", &ports);
    let source = a.try_add(&b).unwrap();
    let negative = Tensor::from_declared_parts(-&b.expression, b.structure.clone()).unwrap();
    let mut visit = |value: &Tensor| {
        Ok(if value.expression == source.expression {
            Action::Descend
        } else if value.same_value(&a) {
            Action::Replace(negative.clone())
        } else {
            assert!(value.same_value(&b));
            Action::Retain
        })
    };
    let actual = source
        .map_original_arithmetic(&mut visit, &mut arithmetic)
        .unwrap();
    let expected = Tensor::from_declared_parts(Atom::Zero, source.structure.clone()).unwrap();
    assert!(actual.same_value(&expected));
    assert_eq!(actual.rank(), 2);
    let error = source
        .map_original_arithmetic(&mut visit, &mut |original, op, children| {
            assert!(original.same_value(&source));
            assert_eq!(op, Operation::Sum);
            assert_eq!(children.len(), 2);
            Ok(scalar(Atom::Zero))
        })
        .unwrap_err();
    assert!(error.to_string().contains("ordered tensor interface"));
}

#[test]
fn closed_bracket_identity_preserves_literal_factor_order_without_product_callback() {
    initialize();
    let a = scalar(Atom::var(symbol!(
        "original_arithmetic::bracket_identity_a"
    )));
    let b = scalar(Atom::var(symbol!(
        "original_arithmetic::bracket_identity_b"
    )));
    // The intentionally reversed order is mathematical input, not an order
    // reconstructed from a parser's arithmetic graph.
    let source = scalar(function!(SPENSO_TAG.bracket, &b.expression, &a.expression));
    let mut seen = Vec::new();
    let actual = source
        .map_closed_bracket_factors(
            |original| {
                seen.push(original.expression.clone());
                Ok(original.clone())
            },
            |_, _| panic!("an unchanged original bracket never constructs a product"),
        )
        .unwrap();
    assert_eq!(seen, vec![b.expression, a.expression]);
    assert!(actual.same_value(&source));
}

#[test]
fn closed_bracket_odd_word_products_follow_literal_argument_order_and_its_sign() {
    initialize();
    let head = symbol!("original_arithmetic::bracket_exterior"; Flat, Linear, Antisymmetric;
        tag = &SPENSO_TAG.tensor_word);
    let a = leaf("original_arithmetic::bracket_odd_a", &[]);
    let b = leaf("original_arithmetic::bracket_odd_b", &[]);
    let c = leaf("original_arithmetic::bracket_odd_c", &[]);
    let marked = |value: &Tensor| Tensor::word(head, std::slice::from_ref(value)).unwrap();
    let old = marked(&a);
    let spectator = marked(&b);
    let replacement = marked(&c);
    let forward = scalar(function!(
        SPENSO_TAG.bracket,
        &old.expression,
        &spectator.expression
    ));
    let reverse = scalar(function!(
        SPENSO_TAG.bracket,
        &spectator.expression,
        &old.expression
    ));
    let mut output_bodies = Vec::new();
    for (source, originals, expected_factors) in [
        (&forward, [&old, &spectator], [&c, &b]),
        (&reverse, [&spectator, &old], [&b, &c]),
    ] {
        let mut seen = Vec::new();
        let mut products = 0;
        let actual = source
            .map_closed_bracket_factors(
                |original| {
                    seen.push(original.expression.clone());
                    Ok(if original.same_value(&old) {
                        replacement.clone()
                    } else {
                        assert!(original.same_value(&spectator));
                        original.clone()
                    })
                },
                |original, images| {
                    products += 1;
                    assert!(original.same_value(source));
                    assert_eq!(images.len(), 2);
                    for (image, original) in images.iter().zip(originals) {
                        assert!(image.same_value(if original.same_value(&old) {
                            &replacement
                        } else {
                            &spectator
                        }));
                    }
                    Ok(images[0].multiply_word(head, &images[1])?)
                },
            )
            .unwrap();
        assert_eq!(seen, originals.map(|value| value.expression.clone()));
        assert_eq!(products, 1);
        let expected_word = Tensor::word(head, &expected_factors.map(Clone::clone)).unwrap();
        assert!(!expected_word.is_zero());
        let expected = scalar(function!(SPENSO_TAG.bracket, &expected_word.expression));
        assert!(actual.same_value(&expected));
        let AtomView::Fun(bracket) = actual.expression.as_view() else {
            panic!("the owner must retain the closed bracket")
        };
        assert_eq!(bracket.get_symbol(), SPENSO_TAG.bracket);
        assert_eq!(bracket.get_nargs(), 1);
        output_bodies.push(bracket.get(0).to_owned());
    }
    // The sign comes entirely from the native exterior owner.
    assert!((&output_bodies[0] + &output_bodies[1]).is_zero());
}

#[test]
fn closed_bracket_images_keep_separate_bound_scopes_and_opaque_metadata_without_hooks() {
    initialize();
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&calls);
    let hook = spenso::tensor_symbol!(
        "original_arithmetic::bracket_callback",
        norm = move |_, _| {
            observed.fetch_add(1, Ordering::Relaxed);
        }
    );
    let rep = LibraryRep::new_dual("original_arithmetic::bracket_bound_flavor")
        .unwrap()
        .new_rep(3);
    let pair = leaf(
        "original_arithmetic::bracket_B",
        &[oriented_slot(rep, 99271)],
    )
    .multiply(&leaf(
        "original_arithmetic::bracket_C",
        &[oriented_slot(rep.dual(), 99271)],
    ))
    .unwrap();
    let first_scope = symbol!("original_arithmetic::bracket_first_scope");
    let second_scope = symbol!("original_arithmetic::bracket_second_scope");
    let new_scope = symbol!("original_arithmetic::bracket_new_scope");
    let first = pair.wrap_indices(first_scope, true).unwrap();
    let second = pair.wrap_indices(second_scope, true).unwrap();
    let image = pair.wrap_indices(new_scope, true).unwrap();
    let metadata = function!(
        symbol!("original_arithmetic::bracket_metadata"; Scalar),
        oriented_slot(rep, 99271)
    );
    let opaque = scalar(function!(hook, &metadata));
    let source = scalar(function!(
        SPENSO_TAG.bracket,
        &first.expression,
        &opaque.expression,
        &second.expression
    ));
    calls.store(0, Ordering::Relaxed);
    let mut seen = Vec::new();
    let actual = source
        .map_closed_bracket_factors(
            |original| {
                seen.push(original.expression.clone());
                Ok(if original.same_value(&first) {
                    image.clone()
                } else {
                    assert!(original.same_value(&second) || original.same_value(&opaque));
                    original.clone()
                })
            },
            |original, images| {
                assert!(original.same_value(&source));
                assert_eq!(images.len(), 3);
                images[1..]
                    .iter()
                    .try_fold(images[0].clone(), |left, right| Ok(left.multiply(right)?))
            },
        )
        .unwrap();
    assert_eq!(
        seen,
        vec![
            first.expression,
            opaque.expression.clone(),
            second.expression.clone()
        ]
    );
    assert_eq!(calls.load(Ordering::Relaxed), 0);
    let expected_body = image.multiply(&opaque).unwrap().multiply(&second).unwrap();
    let expected = scalar(function!(SPENSO_TAG.bracket, &expected_body.expression));
    assert!(actual.same_value(&expected));
    assert!(actual.expression.contains(metadata.as_view()));
    assert!(actual.expression.contains_symbol(second_scope));
    assert!(actual.expression.contains_symbol(new_scope));
    assert!(!actual.expression.contains_symbol(first_scope));
    let AtomView::Fun(bracket) = actual.expression.as_view() else {
        panic!("a changed product remains inside one original closed scope")
    };
    assert_eq!(bracket.get_nargs(), 1);
    assert_eq!(calls.load(Ordering::Relaxed), 0);
}

#[test]
fn closed_bracket_mapping_refuses_open_nonbracket_and_wrong_factor_interfaces() {
    initialize();
    let a = scalar(Atom::var(symbol!("original_arithmetic::not_bracket_a")));
    let b = scalar(Atom::var(symbol!("original_arithmetic::not_bracket_b")));
    let component = leaf("original_arithmetic::open_bracket", &[slot(99281)]);
    let open = Tensor::from_declared_parts(
        function!(SPENSO_TAG.bracket, &component.expression),
        component.structure.clone(),
    )
    .unwrap();
    for original in [a.multiply(&b).unwrap(), open] {
        let error = original
            .map_closed_bracket_factors(
                |_| panic!("invalid original bracket must fail before factor callbacks"),
                |_, _| panic!("invalid original bracket cannot construct a product"),
            )
            .unwrap_err();
        assert!(
            error.to_string().contains("closed scalar interface")
                || error.to_string().contains("original native bracket")
        );
    }
    let partner = leaf(
        "original_arithmetic::closed_bracket_partner",
        &[slot(99281)],
    );
    let source = scalar(function!(
        SPENSO_TAG.bracket,
        &component.expression,
        &partner.expression
    ));
    let error = source
        .map_closed_bracket_factors(
            |_| Ok(scalar(Atom::one())),
            |_, _| panic!("each factor interface must be admitted before product emission"),
        )
        .unwrap_err();
    assert!(error.to_string().contains("ordered tensor interface"));
}

#[test]
fn two_same_representation_auto_operands_keep_their_actual_global_occurrences() {
    initialize();
    let rep = Minkowski {}.new_rep(4);
    let first_head = SPENSO_TAG.tensor_symbol("original_arithmetic::AUTO_first");
    let second_head = SPENSO_TAG.tensor_symbol("original_arithmetic::AUTO_second");
    let changed_head = SPENSO_TAG.tensor_symbol("original_arithmetic::AUTO_changed");
    let first = Tensor::from_signature(&ExplicitKey::<AbstractIndex>::from_iter(
        [rep],
        first_head,
        None,
    ))
    .unwrap();
    let second = Tensor::from_signature(&ExplicitKey::<AbstractIndex>::from_iter(
        [rep],
        second_head,
        None,
    ))
    .unwrap();
    // The native signature owner creates independent original AUTO identities.
    // Infer their literal Product as one original source, retaining both axes.
    let source = Tensor::infer(&first.expression * &second.expression).unwrap();
    assert!(matches!(source.expression.as_view(), AtomView::Mul(_)));
    assert_eq!(source.structure.open_positions(), vec![0, 1]);
    let renaming = HashMap::from([(first_head, changed_head)]);
    let first_image = first.with_renamed_symbols(&renaming).unwrap();
    let expected = source.with_renamed_symbols(&renaming).unwrap();
    let mut first_seen = 0;
    let mut second_seen = 0;
    let result = source.map_original_arithmetic(
        &mut |original| {
            Ok(if original.expression == source.expression {
                Action::Descend
            } else if original.expression == first.expression {
                first_seen += 1;
                assert!(original.same_value(&first));
                Action::Replace(first_image.clone())
            } else {
                second_seen += 1;
                assert_eq!(original.expression, second.expression);
                assert!(original.same_value(&second));
                Action::Retain
            })
        },
        &mut |original, op, images| {
            assert!(original.same_value(&source));
            assert_eq!(op, Operation::Product);
            assert_eq!(images.len(), 2);
            assert!(images.iter().any(|image| image.same_value(&first_image)));
            assert!(images.iter().any(|image| image.same_value(&second)));
            // This original literal Product has two external AUTO axes. Use
            // the native explicit outer owner rather than infer a contraction.
            Ok(images[0].outer(&images[1])?)
        },
    );
    let actual = result.unwrap_or_else(|error| {
        panic!("same-representation AUTO source opening remains an unresolved native capability: {error}")
    });
    assert_eq!(first_seen, 1);
    assert_eq!(second_seen, 1);
    assert_eq!(
        actual.structure.logical_slots(),
        source.structure.logical_slots()
    );
    // Materialize both full global interfaces independently. This detects a
    // detached second occurrence silently becoming another copy of axis zero.
    let ports = HashMap::from([
        (0, AbstractIndex::Normal(99291)),
        (1, AbstractIndex::Normal(99293)),
    ]);
    let actual = actual.reindex_interface_ports(&ports).unwrap();
    let expected = expected.reindex_interface_ports(&ports).unwrap();
    assert!(
        actual
            .canonize()
            .unwrap()
            .same_value(&expected.canonize().unwrap())
    );
}
