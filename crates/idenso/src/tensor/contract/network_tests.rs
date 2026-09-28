//! Regressions carried from the retired symbolic-network strategies.
use super::*;
use crate::shorthands::schoonschip::Schoonschip;
use spenso::{
    network::{library::symbolic::ETS, parsing::StructureFromAtom},
    structure::{
        OrderedStructure, StructureContract, TensorStructure,
        partial::PartialIndex,
        representation::LibrarySlot,
        slot::{AbsInd, IsAbstractSlot, ParseableAind},
    },
};
use std::sync::{Arc, Mutex};
use symbolica::atom::{Atom, AtomCore, FunctionBuilder};

fn callback_contraction(
    metric: &SymbolicTensor,
    tensor: &SymbolicTensor,
) -> Result<SymbolicTensor, TensorInferenceError> {
    let structure = metric
        .structure
        .merge(&tensor.structure)
        .map_err(|error| TensorInferenceError::invalid(error.to_string()))?
        .0;
    let interface = PartialStructure::from_logical_slots(
        structure
            .external_structure()
            .into_iter()
            .map(|slot| slot.rep().slot(PartialIndex::Explicit(slot.aind()))),
    );
    let source = SymbolicTensor::from_normalized_parts(
        Atom::mul_many([&metric.expression, &tensor.expression]),
        interface,
    );
    let expression = source
        .contract(Default::default())?
        .resolved()?
        .into_expression();
    let value = SymbolicTensor {
        expression,
        structure,
        is_metric: false,
        is_composite: true,
        proofs: Default::default(),
    };
    value.validate_rewritten_interface(&value.expression)?;
    Ok(value)
}

#[test]
fn network_rejects_callback_rank_loss_before_reusing_merged_ports() {
    crate::test_support::test_initialize();
    let parse = |source| {
        Atom::parse(
            source,
            "network_callback_rank_loss",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap()
    };
    let a = parse("spenso::mink(4,a)");
    let b = parse("spenso::mink(4,b)");
    let events = Arc::new(Mutex::new(Vec::new()));
    let seen = Arc::clone(&events);
    let target = b.clone();
    let head = spenso::tensor_symbol!(
        "network_callback_rank_loss_T",
        norm = move |node, out| {
            if let AtomView::Fun(function) = node
                && function.iter().last() == Some(target.as_view())
            {
                seen.lock().unwrap().push(node.to_owned());
                **out = Atom::one();
            }
        }
    );
    let metric_atom = FunctionBuilder::new(ETS.metric)
        .add_arg(&a)
        .add_arg(&b)
        .finish();
    let tensor_atom = FunctionBuilder::new(head).add_arg(&a).finish();
    let metric = SymbolicTensor::parse(metric_atom.as_view())
        .unwrap()
        .into_canonical();
    let tensor = SymbolicTensor::parse(tensor_atom.as_view())
        .unwrap()
        .into_canonical();

    events.lock().unwrap().clear();
    let error = callback_contraction(&metric, &tensor).unwrap_err();
    assert!(error.to_string().contains("compatible tensor interface"));
    assert_eq!(events.lock().unwrap().len(), 1);
}

#[test]
fn network_callback_validation_preserves_zero_ports_and_does_not_replay() {
    crate::test_support::test_initialize();
    let parse = |source| {
        Atom::parse(
            source,
            "network_callback_zero",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap()
    };
    let a = parse("spenso::mink(4,a)");
    let b = parse("spenso::mink(4,b)");
    let events = Arc::new(Mutex::new(Vec::new()));
    let seen = Arc::clone(&events);
    let target = b.clone();
    let head = spenso::tensor_symbol!(
        "network_callback_zero_T",
        norm = move |node, out| {
            if let AtomView::Fun(function) = node
                && function.iter().last() == Some(target.as_view())
            {
                seen.lock().unwrap().push(node.to_owned());
                **out = Atom::Zero;
            }
        }
    );
    let metric_atom = FunctionBuilder::new(ETS.metric)
        .add_arg(&a)
        .add_arg(&b)
        .finish();
    let tensor_atom = FunctionBuilder::new(head).add_arg(&a).finish();
    let vector_atom = FunctionBuilder::new(spenso::vector_symbol!("network_callback_zero_Q"))
        .add_arg(&b)
        .finish();
    let metric = SymbolicTensor::parse(metric_atom.as_view())
        .unwrap()
        .into_canonical();
    let tensor = SymbolicTensor::parse(tensor_atom.as_view())
        .unwrap()
        .into_canonical();
    let vector = SymbolicTensor::parse(vector_atom.as_view())
        .unwrap()
        .into_canonical();
    let planned = metric.structure.merge(&tensor.structure).unwrap().0;

    events.lock().unwrap().clear();
    let result = callback_contraction(&metric, &tensor).unwrap();
    assert_eq!(result.expression, Atom::Zero);
    assert_eq!(result.structure, planned);
    assert_eq!(events.lock().unwrap().len(), 1);
    let closed = callback_contraction(&result, &vector).unwrap();
    assert_eq!(closed.expression, Atom::Zero);
    assert!(closed.structure.is_scalar());
    assert_eq!(events.lock().unwrap().len(), 1);
}

#[test]
fn network_callback_validation_observes_indices_and_retains_planned_order() {
    crate::test_support::test_initialize();
    let parse = |source| {
        Atom::parse(
            source,
            "network_callback_order",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap()
    };
    let a = parse("spenso::mink(4,a)");
    let b = parse("spenso::mink(4,b)");
    let c = parse("spenso::mink(4,c)");
    let events = Arc::new(Mutex::new(Vec::new()));
    let seen = Arc::clone(&events);
    let head = spenso::tensor_symbol!(
        "network_callback_order_T",
        norm = move |node, _out| {
            seen.lock().unwrap().push(node.to_owned());
        }
    );
    let metric_atom = FunctionBuilder::new(ETS.metric)
        .add_arg(&a)
        .add_arg(&b)
        .finish();
    let tensor_atom = FunctionBuilder::new(head).add_arg(&c).add_arg(&a).finish();
    let metric = SymbolicTensor::parse(metric_atom.as_view())
        .unwrap()
        .into_canonical();
    let tensor = SymbolicTensor::parse(tensor_atom.as_view())
        .unwrap()
        .into_canonical();
    let planned = metric.structure.merge(&tensor.structure).unwrap().0;
    events.lock().unwrap().clear();
    let expected = tensor_atom
        .replace(a.to_pattern())
        .with(b.to_pattern())
        .normalize_dots();
    let expected_events = std::mem::take(&mut *events.lock().unwrap());
    let result = callback_contraction(&metric, &tensor).unwrap();
    assert_eq!(result.expression, expected);
    assert_eq!(result.structure, planned);
    assert_eq!(*events.lock().unwrap(), expected_events);
    result
        .validate_rewritten_interface(&result.expression)
        .unwrap();
    assert_eq!(*events.lock().unwrap(), expected_events);
    let wrong_index = result
        .expression
        .replace(b.to_pattern())
        .with(a.to_pattern());
    assert!(result.validate_rewritten_interface(&wrong_index).is_err());
}

#[test]
fn observed_network_interface_accepts_generic_index_storage() {
    use spenso::structure::abstract_index::{AbstractIndex, AbstractIndexError};

    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
    struct WrappedIndex(AbstractIndex);
    impl std::fmt::Display for WrappedIndex {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            std::fmt::Display::fmt(&self.0, f)
        }
    }
    impl AbsInd for WrappedIndex {}
    impl ParseableAind for WrappedIndex {
        type Error = AbstractIndexError;
        fn from_view(view: AtomView<'_>) -> Result<Self, Self::Error> {
            AbstractIndex::from_view(view).map(Self)
        }
        fn to_atom(&self) -> Atom {
            self.0.to_atom()
        }
    }

    crate::test_support::test_initialize();
    let head = spenso::tensor_symbol!("network_generic_index_T");
    let slot = Atom::parse(
        "mink(4,network_generic_index_a)",
        "spenso",
        symbolica::parser::ParseSettings::symbolica(),
    )
    .unwrap();
    let expression = FunctionBuilder::new(head).add_arg(&slot).finish();
    let port = LibrarySlot::<WrappedIndex>::try_from(slot.as_view()).unwrap();
    let value = SymbolicTensor {
        proofs: Default::default(),
        structure: OrderedStructure::new(vec![port]).into_canonical(),
        expression,
        is_composite: false,
        is_metric: false,
    };
    value
        .validate_rewritten_interface(&value.expression)
        .unwrap();
    assert!(value.validate_rewritten_interface(&Atom::one()).is_err());
    value.validate_rewritten_interface(&Atom::Zero).unwrap();
    assert_eq!(value.structure.external_structure(), vec![port]);
}

#[test]
fn observed_network_interface_keeps_distinct_open_owners() {
    use spenso::structure::abstract_index::AbstractIndex;

    crate::test_support::test_initialize();
    let slot = Atom::parse(
        "mink(4,network_open_owner_a)",
        "spenso",
        symbolica::parser::ParseSettings::symbolica(),
    )
    .unwrap();
    let rep = LibrarySlot::<AbstractIndex>::try_from(slot.as_view())
        .unwrap()
        .rep();
    let first = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
        owner: 101,
        axis: 0,
    });
    let second = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
        owner: 102,
        axis: 0,
    });
    let events = Arc::new(Mutex::new(Vec::new()));
    let seen = Arc::clone(&events);
    let head = spenso::tensor_symbol!(
        "network_distinct_open_owners_T",
        norm = move |node, _out| {
            seen.lock().unwrap().push(node.to_owned());
        }
    );
    let expression = FunctionBuilder::new(head)
        .add_arg(first.to_atom())
        .add_arg(second.to_atom())
        .finish();
    let value = SymbolicTensor {
        proofs: Default::default(),
        structure: OrderedStructure::new(vec![first, second]).into_canonical(),
        expression,
        is_composite: false,
        is_metric: false,
    };
    let construction_events = events.lock().unwrap().clone();
    value
        .validate_rewritten_interface(&value.expression)
        .unwrap();
    assert_eq!(*events.lock().unwrap(), construction_events);
    assert_eq!(value.structure.order(), 2);
    assert_ne!(first.aind(), second.aind());
    let planned = value.structure.clone();
    for port in [first, second] {
        let dropped = FunctionBuilder::new(head).add_arg(port.to_atom()).finish();
        let before_validation = events.lock().unwrap().clone();
        assert!(value.validate_rewritten_interface(&dropped).is_err());
        assert_eq!(*events.lock().unwrap(), before_validation);
        assert_eq!(value.structure, planned);
    }
}

#[test]
fn network_callback_rejects_changed_encoded_open_owner() {
    use spenso::structure::{
        abstract_index::AbstractIndex,
        dimension::Dimension,
        representation::{ExtendibleReps, RepName},
    };

    crate::test_support::test_initialize();
    let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
    let first = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
        owner: 301,
        axis: 0,
    });
    let second = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
        owner: 302,
        axis: 0,
    });
    let changed = rep.slot::<AbstractIndex, _>(AbstractIndex::Open {
        owner: 303,
        axis: 0,
    });
    let contracted = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(41));
    let target = first.to_atom();
    let replacement = FunctionBuilder::new(spenso::tensor_symbol!("network_changed_owner_U"))
        .add_arg(second.to_atom())
        .add_arg(changed.to_atom())
        .finish();
    let events = Arc::new(Mutex::new(Vec::new()));
    let seen = Arc::clone(&events);
    let head = spenso::tensor_symbol!(
        "network_changed_owner_T",
        norm = move |node, out| {
            if let AtomView::Fun(function) = node
                && function.iter().last() == Some(target.as_view())
            {
                seen.lock().unwrap().push(node.to_owned());
                **out = replacement.clone();
            }
        }
    );
    let metric_atom = FunctionBuilder::new(ETS.metric)
        .add_arg(contracted.to_atom())
        .add_arg(first.to_atom())
        .finish();
    let tensor_atom = FunctionBuilder::new(head)
        .add_arg(second.to_atom())
        .add_arg(contracted.to_atom())
        .finish();
    let metric = SymbolicTensor::parse(metric_atom.as_view())
        .unwrap()
        .into_canonical();
    let tensor = SymbolicTensor::parse(tensor_atom.as_view())
        .unwrap()
        .into_canonical();

    events.lock().unwrap().clear();
    let error = callback_contraction(&metric, &tensor).unwrap_err();
    assert!(error.to_string().contains("compatible tensor interface"));
    assert_eq!(events.lock().unwrap().len(), 1);
}

#[test]
fn network_callback_preserves_both_mixed_representation_argument_orders() {
    use spenso::structure::{
        ScalarTensor,
        abstract_index::AbstractIndex,
        dimension::Dimension,
        representation::{ExtendibleReps, RepName},
    };

    crate::test_support::test_initialize();
    let mink = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
    let euc = ExtendibleReps::EUCLIDEAN.new_rep(Dimension::Concrete(3));
    let first = mink.slot::<AbstractIndex, _>(AbstractIndex::Open {
        owner: 401,
        axis: 0,
    });
    let second = euc.slot::<AbstractIndex, _>(AbstractIndex::Open {
        owner: 402,
        axis: 0,
    });
    let events = Arc::new(Mutex::new(Vec::new()));
    let seen = Arc::clone(&events);
    let head = spenso::tensor_symbol!(
        "network_mixed_owner_order_T",
        norm = move |node, _out| {
            seen.lock().unwrap().push(node.to_owned());
        }
    );
    let scalar = SymbolicTensor::new_scalar(Atom::num(2));
    for ports in [[first, second], [second, first]] {
        let expression = FunctionBuilder::new(head)
            .add_args(ports.map(|slot| slot.to_atom()))
            .finish();
        let value = SymbolicTensor::parse(expression.as_view())
            .unwrap()
            .into_canonical();
        let expected = (&scalar.expression * &expression).normalize_dots();
        let before = events.lock().unwrap().clone();
        let result = callback_contraction(&scalar, &value).unwrap();
        assert_eq!(result.expression, expected);
        assert_eq!(result.structure, value.structure);
        assert_eq!(*events.lock().unwrap(), before);
    }
}
