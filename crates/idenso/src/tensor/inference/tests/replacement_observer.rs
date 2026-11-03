//! Public replacement observation is syntax-preserving, not materialization.

use crate::tensor::inference::InterfaceInference;
use crate::test_support;
use spenso::structure::TensorStructure;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::{PartialIndex, PartialStructureExt};
use spenso::structure::representation::{LibraryRep, Minkowski, RepName};
use spenso::structure::slot::{DualSlotTo, IsAbstractSlot};
use spenso::{tensor_symbol, vector_symbol};
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use symbolica::atom::{Atom, AtomCore, FunctionBuilder};
use symbolica::{function, symbol};

#[test]
fn replacement_observer_does_not_replay_scalar_metadata_or_tensor_callbacks() {
    test_support::test_initialize();
    let metadata_calls = Arc::new(AtomicUsize::new(0));
    let tensor_calls = Arc::new(AtomicUsize::new(0));
    let seen = Arc::clone(&metadata_calls);
    let metadata = symbol!(
        "replacement_observer::private_metadata"; Scalar;
        norm = move |_, _| {
            seen.fetch_add(1, Ordering::Relaxed);
        }
    );
    let seen = Arc::clone(&tensor_calls);
    let head = tensor_symbol!(
        "replacement_observer::observed_tensor",
        norm = move |_, _| {
            seen.fetch_add(1, Ordering::Relaxed);
        }
    );
    let rep = LibraryRep::from(Minkowski {}).new_rep(4);
    let external = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(991001));
    let hidden = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(991002));
    let private_body = function!(
        tensor_symbol!("replacement_observer::private_body"),
        hidden.to_atom()
    )
    .pow(2);
    let expression = function!(head, function!(metadata, &private_body), external.to_atom());
    let original = expression.clone();
    assert!(metadata_calls.load(Ordering::Relaxed) > 0);
    assert!(tensor_calls.load(Ordering::Relaxed) > 0);
    metadata_calls.store(0, Ordering::Relaxed);
    tensor_calls.store(0, Ordering::Relaxed);

    let observed = InterfaceInference::replacement_interface(expression.as_view()).unwrap();
    assert_eq!(
        observed.logical_slots(),
        vec![external.rep().slot(PartialIndex::Explicit(external.aind()))]
    );
    assert_eq!(expression, original);
    assert_eq!(metadata_calls.load(Ordering::Relaxed), 0);
    assert_eq!(tensor_calls.load(Ordering::Relaxed), 0);
}

#[test]
fn replacement_observer_retains_encoded_open_owners_and_axis_order() {
    test_support::test_initialize();
    let rep = LibraryRep::from(Minkowski {}).new_rep(4);
    let first = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
        owner: 991011,
        axis: 7,
    });
    let second = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
        owner: 991012,
        axis: 3,
    });
    let expression = function!(
        tensor_symbol!("replacement_observer::encoded_axes"),
        first.to_atom(),
        second.to_atom()
    );
    let observed = InterfaceInference::replacement_interface(expression.as_view()).unwrap();
    assert_eq!(
        observed.logical_slots(),
        vec![
            first.rep().slot(PartialIndex::Explicit(first.aind())),
            second.rep().slot(PartialIndex::Explicit(second.aind())),
        ]
    );
    assert!(observed.open_positions().is_empty());
}

#[test]
fn replacement_observer_distinguishes_direct_dual_and_closed_ports() {
    test_support::test_initialize();
    let rep = LibraryRep::new_dual("replacement_observer::R")
        .unwrap()
        .new_rep(3);
    let slot = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(991021));
    let dual = slot.dual();
    let left = function!(tensor_symbol!("replacement_observer::left"), slot.to_atom());
    let right = function!(
        tensor_symbol!("replacement_observer::right"),
        dual.to_atom()
    );
    assert_eq!(
        InterfaceInference::replacement_interface(left.as_view())
            .unwrap()
            .logical_slots(),
        vec![slot.rep().slot(PartialIndex::Explicit(slot.aind()))]
    );
    assert_eq!(
        InterfaceInference::replacement_interface(right.as_view())
            .unwrap()
            .logical_slots(),
        vec![dual.rep().slot(PartialIndex::Explicit(dual.aind()))]
    );
    let closed = &left * &right;
    assert!(
        InterfaceInference::replacement_interface(closed.as_view())
            .unwrap()
            .canonical()
            .is_scalar()
    );
}

#[test]
fn replacement_observer_observes_compact_ports_without_materialization() {
    test_support::test_initialize();
    let calls = Arc::new(AtomicUsize::new(0));
    let seen = Arc::clone(&calls);
    let head = vector_symbol!(
        "replacement_observer::compact_vector",
        norm = move |_, _| {
            seen.fetch_add(1, Ordering::Relaxed);
        }
    );
    let rep = LibraryRep::from(Minkowski {}).new_rep(4);
    let expression = FunctionBuilder::new(head)
        .add_arg(Atom::num(7))
        .add_arg(rep.to_symbolic([]))
        .finish();
    let original = expression.clone();
    assert!(calls.load(Ordering::Relaxed) > 0);
    calls.store(0, Ordering::Relaxed);
    let observed = InterfaceInference::replacement_interface(expression.as_view()).unwrap();
    assert_eq!(observed.canonical().order(), 1);
    assert_eq!(observed.open_positions(), vec![0]);
    assert_eq!(expression, original);
    assert_eq!(calls.load(Ordering::Relaxed), 0);
}

#[test]
fn replacement_observer_accepts_ordinary_scalars_without_a_carried_zero_rank() {
    test_support::test_initialize();
    for expression in [
        Atom::Zero,
        Atom::one(),
        Atom::var(symbol!("replacement_observer::a")),
    ] {
        assert!(
            InterfaceInference::replacement_interface(expression.as_view())
                .unwrap()
                .canonical()
                .is_scalar()
        );
    }
}

#[test]
fn replacement_observer_reports_malformed_syntax_separately_from_valid_open_ports() {
    test_support::test_initialize();
    let rep = LibraryRep::from(Minkowski {}).new_rep(4);
    let slot = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(991031));
    let head = tensor_symbol!("replacement_observer::malformed");
    let open = function!(head, slot.to_atom());
    assert_eq!(
        InterfaceInference::replacement_interface(open.as_view())
            .unwrap()
            .canonical()
            .order(),
        1
    );
    let other = function!(
        tensor_symbol!("replacement_observer::other"),
        slot.to_atom()
    );
    let third = function!(
        tensor_symbol!("replacement_observer::third"),
        slot.to_atom()
    );
    // Generic tensor functions may have rank zero. Unlike a bare tensor
    // symbol, the zero-argument function is already a complete closed leaf.
    let closed = InterfaceInference::replacement_interface(function!(head).as_view()).unwrap();
    assert!(closed.canonical().is_scalar());
    assert!(closed.logical_slots().is_empty());

    // Rank-one heads have the stronger native one-final-port contract.
    let missing_rank_one_port = function!(vector_symbol!(
        "replacement_observer::rank_one_missing_port"
    ));
    for (label, malformed) in [
        ("bare_tensor_variable", Atom::var(head)),
        ("rank_one_missing_port", missing_rank_one_port),
        ("triple_repeated_slot", open * other * third),
    ] {
        let result = InterfaceInference::replacement_interface(malformed.as_view());
        assert!(result.is_err(), "malformed case {label} was accepted");
        if let Err(error) = result {
            eprintln!("{label}: {error}");
        }
    }
}
