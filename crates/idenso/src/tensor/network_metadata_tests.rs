//! Network execution must preserve the admission owner's scalar metadata.

use super::*;
use crate::representations::ColorFundamental;
use crate::tensor::inference::tests::{
    INFERENCE_CALLS, INFERENCE_VISITS, OCCURRENCE_FUNCTION_VISITS,
};
use crate::test_support;
use spenso::network::parsing::{LeafInterfaceCache, ShorthandParsing};
use spenso::structure::representation::RepName;
use spenso::structure::slot::DualSlotTo;
use spenso::{g, mink, tensor_symbol, vector_symbol};
use symbolica::symbol;

thread_local! {
    pub(super) static RAW_PARSER_ADMISSIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

fn contracted_pair(index: usize) -> Atom {
    let port = mink!(4, Atom::num(index));
    function!(tensor_symbol!("network_metadata_left"), &port)
        * function!(tensor_symbol!("network_metadata_right"), &port)
}

#[test]
fn planning_leaves_share_source_facts_until_selected() {
    test_support::test_initialize();
    let first = function!(tensor_symbol!("planning_source_first"), mink!(4, 99857));
    let second = function!(tensor_symbol!("planning_source_second"), mink!(4, 99859));
    let tensor = SymbolicTensor::infer(&first * &second).unwrap();
    OCCURRENCE_FUNCTION_VISITS.with(|visits| visits.set(0));
    let plan = tensor.shallow_graph().unwrap();
    assert!(OCCURRENCE_FUNCTION_VISITS.with(|visits| visits.get()) > 0);
    let leaves = &plan.store.tensors;
    assert_eq!(leaves.len(), 2);
    assert!(leaves.iter().all(|leaf| leaf.local.is_none()));
    let shared = Arc::downgrade(plan.admission.as_ref().unwrap());
    let selected = leaves.iter().find(|leaf| leaf.expression == first).unwrap();
    // A selected occurrence can retain a dummy which cancelled out of its
    // source bytes. Reopening checks it only at this same root.
    let mut scoped = selected.clone();
    let hidden = selected
        .structure
        .get_slot(0)
        .unwrap()
        .rep()
        .slot(AbstractIndex::Normal(99861));
    Arc::make_mut(&mut scoped.history).add(hidden, 2);
    OCCURRENCE_FUNCTION_VISITS.with(|visits| visits.set(0));
    let reopened = planning::parse(
        first.as_view(),
        None,
        Some(scoped.source(first.as_view(), plan.admission.as_ref())),
        &ParseSettings::default(),
    )
    .unwrap();
    assert!(Arc::ptr_eq(
        reopened.admission.as_ref().unwrap(),
        plan.admission.as_ref().unwrap()
    ));
    assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|visits| visits.get()), 0);
    assert!(Arc::ptr_eq(
        &reopened.store.tensors[0].history,
        &scoped.history
    ));
    assert_eq!(reopened.store.tensors[0].history.count(hidden), 2);
    drop(reopened);
    let mut materialized = SymbolicTensor::from_validated_parts(
        selected.expression.clone(),
        TensorProofs::parser_boundary(&selected.structure).unwrap(),
        Arc::clone(&selected.history),
    );
    materialized
        .proofs
        .inherit_interfaces([plan.evidence(selected)]);
    let proofs = materialized.proofs;
    assert!(proofs.owns_expression(first.as_view()));
    let written = proofs.parser_leaf_interfaces();
    let written = written.lock().unwrap();
    assert_eq!(written[&first].1, [mink!(4, 99857)]);
    assert!(!written.contains_key(&second));
    assert!(proofs.retained_interface(second.as_view()).is_none());
    drop(plan);
    drop(tensor);
    assert!(shared.upgrade().is_none());
}

#[test]
fn planning_checks_materialized_leaf_history_before_reuse() {
    test_support::test_initialize();
    let expression = function!(
        tensor_symbol!("planning_materialized_source"),
        mink!(4, 99869)
    );
    let settings = ParseSettings::default();
    let network = expression
        .parse_to_symbolic_net::<AbstractIndex>(&settings)
        .unwrap();
    let plan = planning::PlanningNetwork::from_materialized(network);
    let leaf = &plan.store.tensors[0];
    OCCURRENCE_FUNCTION_VISITS.with(|visits| visits.set(0));
    let reopened = planning::parse(
        expression.as_view(),
        None,
        Some(leaf.source(expression.as_view(), plan.admission.as_ref())),
        &settings,
    )
    .unwrap();
    assert!(OCCURRENCE_FUNCTION_VISITS.with(|visits| visits.get()) > 0);
    assert_eq!(reopened.store.tensors[0].structure, leaf.structure);
}

#[test]
fn planning_keeps_graph_storage_order_and_written_constructor_ports() {
    test_support::test_initialize();
    let first = mink!(4, 99871);
    let second = mink!(4, 99873);
    let expression = function!(tensor_symbol!("planning_reverse_ports"), &second, &first);
    let declared = SymbolicTensor::infer(expression.clone()).unwrap().structure;
    let contracted = contracted_pair(99875);
    let input = SymbolicTensor::from_declared_parts(&expression * contracted, declared).unwrap();
    let tensor = input
        .with_materialized_expression("remove a contracted scalar", expression.clone())
        .unwrap();
    let plan = tensor.shallow_graph().unwrap();
    let leaf = &plan.store.tensors[0];
    let storage = TensorProofs::parser_boundary(&leaf.structure).unwrap();
    assert_eq!(tensor.structure.canonical(), storage.canonical());
    assert_ne!(tensor.structure.layout(), storage.layout());
    let hidden = leaf
        .structure
        .get_slot(0)
        .unwrap()
        .rep()
        .slot(AbstractIndex::Normal(99875));
    assert_eq!(leaf.history.count(hidden), 2);
    assert!(leaf.local.is_some());
    let facts = plan.evidence(leaf);
    assert_eq!(
        facts.retained_interface(expression.as_view()),
        Some(storage.clone())
    );
    assert_eq!(
        facts.parser_leaf_interfaces().lock().unwrap()[&expression].1,
        [second, first]
    );
    let reopened = planning::parse(
        expression.as_view(),
        None,
        Some(leaf.source(expression.as_view(), plan.admission.as_ref())),
        &ParseSettings::default(),
    )
    .unwrap();
    let reopened_leaf = &reopened.store.tensors[0];
    assert!(reopened_leaf.local.is_none());
    assert_eq!(reopened_leaf.history.count(hidden), 2);
    assert!(Arc::ptr_eq(&leaf.history, &reopened_leaf.history));
    assert_eq!(
        reopened
            .evidence(reopened_leaf)
            .retained_interface(expression.as_view()),
        Some(storage)
    );
}

#[test]
fn collection_rechecks_product_readiness_for_changed_expression() {
    test_support::test_initialize();
    let port = mink!(4, 99887);
    let selected = function!(tensor_symbol!("collection_original_leaf"), &port);
    let replacement = function!(tensor_symbol!("collection_rewritten_leaf"), &port);
    let weight = Atom::var(symbol!("collection_memo_weight"));
    let source = SymbolicTensor::infer(&weight * &selected).unwrap();
    assert!(source.borrowed_parts().product_is_normalized());
    let output = source
        .collect_with_map(
            CollectionMode::Monomials,
            None,
            |leaf| leaf == selected.as_view(),
            |value, _, _| value.with_checked_expression(replacement.clone()),
        )
        .unwrap();
    assert_eq!(output.expression, &weight * &replacement);
    assert_eq!(source.proofs.normalized_product.get(), Some(&true));
    assert!(output.proofs.normalized_product.get().is_none());
    assert!(output.borrowed_parts().product_is_normalized());
}

#[test]
fn collection_keeps_the_enclosing_declaration_in_both_leaf_orders() {
    test_support::test_initialize();
    let reversed = function!(
        tensor_symbol!("collection_reverse_ports"),
        mink!(4, 99883),
        mink!(4, 99881)
    );
    let sibling = function!(tensor_symbol!("collection_single_port"), mink!(4, 99885));
    let tensor = SymbolicTensor::infer(&reversed * &sibling).unwrap();
    let plan = tensor.shallow_graph().unwrap();
    let leaf = plan
        .store
        .tensors
        .iter()
        .find(|leaf| leaf.expression == reversed)
        .unwrap();
    let sibling = plan
        .store
        .tensors
        .iter()
        .find(|leaf| leaf.expression == sibling)
        .unwrap();
    let declared = tensor
        .proofs
        .retained_interface(reversed.as_view())
        .unwrap();
    assert_ne!(
        declared,
        plan.evidence(leaf)
            .retained_interface(reversed.as_view())
            .unwrap()
    );
    for inputs in [[leaf, sibling], [sibling, leaf]] {
        let mut output = SymbolicTensor::from_validated_parts(
            tensor.expression.clone(),
            tensor.structure.clone(),
            Arc::clone(tensor.proofs.occurrences.get().unwrap()),
        );
        plan.inherit_leaf_interfaces(&mut output.proofs, inputs);
        assert_eq!(
            output
                .proofs
                .retained_interface(reversed.as_view())
                .as_ref(),
            Some(&declared)
        );
    }
}

#[test]
fn planning_keeps_local_facts_for_a_derived_only_leaf() {
    test_support::test_initialize();
    let slot = mink!(4, 99863);
    let expression = function!(tensor_symbol!("planning_derived_leaf"), &slot);
    let producer = TensorProofs::default();
    producer
        .retain_leaf_ports(expression.as_view(), vec![slot.clone()])
        .unwrap();
    let plan = planning::parse(
        expression.as_view(),
        None,
        Some(planning::PlanningSource::root(&producer)),
        &ParseSettings::default(),
    )
    .unwrap();
    let leaf = &plan.store.tensors[0];
    assert!(leaf.local.is_some());
    producer
        .retain_leaf_ports(expression.as_view(), vec![mink!(4, 99865)])
        .unwrap();
    assert_eq!(
        plan.evidence(leaf).parser_leaf_interfaces().lock().unwrap()[&expression].1,
        [slot],
    );
}

#[test]
fn parsed_history_keeps_written_orientation_and_independent_storage_slots() {
    test_support::test_initialize();
    let group = LibraryRep::from(ColorFundamental {}).new_rep(3);
    let written = group.slot::<AbstractIndex, _>(AbstractIndex::Normal(99851));
    let storage = group.slot::<AbstractIndex, _>(AbstractIndex::Normal(99853));
    let head = tensor_symbol!("network_history_storage");
    let expression = function!(head, written.dual().to_atom(), 7, written.to_atom());
    let structure = OrderedStructure::new(vec![storage.dual(), storage]);
    let proofs = TensorProofs::parsed_leaf(
        expression.as_view(),
        &structure,
        Some(&TensorProofs::default()),
    )
    .unwrap();
    let occurrences = proofs.occurrences.get().unwrap();
    assert_eq!(
        occurrences.oriented_entries().next(),
        Some((written.dual(), 2, 1, 1)),
    );
    assert_eq!(occurrences.count(storage), 2);
    let reassigned = group.slot::<AbstractIndex, _>(AbstractIndex::Normal(99854));
    let changed = OrderedStructure::new(vec![reassigned.dual(), reassigned]);
    let reopened =
        TensorProofs::with_parser_ports(Arc::clone(occurrences), changed.canonical()).unwrap();
    assert_eq!(reopened.count(reassigned), 2);
    assert_eq!(reopened.count(storage), 2);
    assert_eq!(reopened.count(written), 2);
    assert!(
        proofs
            .retained_interface(expression.as_view())
            .unwrap()
            .canonical()
            .is_scalar()
    );

    let invalid = function!(
        head,
        written.to_atom(),
        written.dual().to_atom(),
        written.to_atom()
    );
    assert!(TensorProofs::parsed_leaf(invalid.as_view(), &structure, Some(&proofs)).is_err());
}

#[test]
fn admitted_compact_constructors_keep_raw_ports_for_typed_parsing() {
    test_support::test_initialize();
    let compact = mink!(4);
    let metadata = function!(symbol!("network_compact_metadata"; Scalar), 7);
    let generic = tensor_symbol!("network_compact_generic");
    let inputs = [
        (function!(generic, &compact), vec![compact.clone()]),
        (
            function!(generic, metadata, &compact),
            vec![compact.clone()],
        ),
        (
            function!(vector_symbol!("network_compact_vector"), 0, &compact),
            vec![compact.clone()],
        ),
        (
            g!(&compact, &compact),
            vec![compact.clone(), compact.clone()],
        ),
    ];
    for (expression, ports) in inputs {
        let source = SymbolicTensor::infer(expression).unwrap();
        assert_eq!(
            source
                .proofs
                .leaf_interfaces
                .get()
                .unwrap()
                .lock()
                .unwrap()
                .get(source.expression())
                .unwrap()
                .1,
            ports
        );
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        for _ in 0..2 {
            let network = source
                .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
                .unwrap();
            assert_eq!(network.graph.dangling_indices().len(), ports.len());
        }
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }
}

#[test]
fn named_products_publish_their_original_and_compact_vector_ports() {
    test_support::test_initialize();
    let port = mink!(4, 99845);
    let compact = mink!(4);
    let generic = function!(tensor_symbol!("network_batched_generic"), &port);
    let left = vector_symbol!("network_batched_left");
    let right = vector_symbol!("network_batched_right");
    let left_vector = SymbolicTensor::infer(function!(left, &port)).unwrap();
    let right_vector = SymbolicTensor::infer(function!(right, &port)).unwrap();
    let scalar = SymbolicTensor::infer(Atom::num(3)).unwrap();
    let generic_value = SymbolicTensor::infer(generic.clone()).unwrap();
    INFERENCE_CALLS.with(|count| count.set(0));
    RAW_PARSER_ADMISSIONS.with(|count| count.set(0));

    let batched = SymbolicTensor::multiply_many([&scalar, &generic_value]).unwrap();
    let ports = batched.proofs.parser_leaf_interfaces();
    assert_eq!(ports.lock().unwrap()[&generic].1, [port]);
    for product in [
        left_vector.dot(&right_vector).unwrap(),
        left_vector.multiply(&right_vector).unwrap(),
        SymbolicTensor::multiply_many([&left_vector, &right_vector]).unwrap(),
    ] {
        let ports = product.proofs.parser_leaf_interfaces();
        for head in [left, right] {
            assert_eq!(
                ports.lock().unwrap()[&function!(head, &compact)]
                    .1
                    .as_slice(),
                std::slice::from_ref(&compact)
            );
        }
    }
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(RAW_PARSER_ADMISSIONS.with(|count| count.get()), 0);
}

#[test]
fn fresh_known_producers_publish_raw_ports_without_inference() {
    test_support::test_initialize();
    let left = mink!(4, 99841);
    let right = mink!(4, 99843);
    let vector = vector_symbol!("network_produced_vector");
    for expression in [
        g!(&left, &right),
        function!(vector, 0, &left),
        function!(vector, 0, &left) * function!(vector, 1, &left),
    ] {
        let declared = SymbolicTensor::infer(expression.clone()).unwrap().structure;
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let value = SymbolicTensor::from_produced_parts(expression, declared).unwrap();
        let network = value
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        assert_eq!(network.graph.dangling_indices().len(), value.rank());
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }
}

fn execute(network: &mut SymbolicNet<AbstractIndex>) -> SymbolicTensor {
    let library = DummyLibrary::<SymbolicTensor>::new();
    network
        .execute::<Sequential, SmallestDegree, _, _, _>(&library, &Wrap {})
        .unwrap();
    match network.result_tensor(&library).unwrap() {
        ExecutionResult::Val(value) => value.into_owned(),
        _ => panic!("the result must retain its scalar tensor carrier"),
    }
}

#[test]
fn admitted_network_scalar_sums_retain_facts_without_inference() {
    test_support::test_initialize();
    let source = SymbolicTensor::infer(contracted_pair(99781) + 3).unwrap();
    for depth_limit in [None, Some(0)] {
        let settings = ParseSettings {
            depth_limit,
            precontract_scalars: false,
            shorthand_parsing: ShorthandParsing::Opaque,
            ..ParseSettings::default()
        };
        RAW_PARSER_ADMISSIONS.with(|count| count.set(0));
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let mut network = source
            .parse_to_symbolic_net::<AbstractIndex>(&settings)
            .unwrap();
        let result = execute(&mut network);
        assert_eq!(result.expression, source.expression);
        assert!(result.structure.is_scalar());
        assert!(result.proofs.occurrences.get().unwrap().is_empty());
        assert!(result.proofs.owns_expression(result.expression.as_view()));
        assert_eq!(RAW_PARSER_ADMISSIONS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }
}

#[test]
fn executed_scalar_power_keeps_its_contraction_scope_without_inference() {
    test_support::test_initialize();
    let source = SymbolicTensor::infer(contracted_pair(99783)).unwrap();
    let settings = ParseSettings {
        precontract_scalars: false,
        shorthand_parsing: ShorthandParsing::Opaque,
        ..ParseSettings::default()
    };
    RAW_PARSER_ADMISSIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));
    let mut network = source
        .parse_to_symbolic_net::<AbstractIndex>(&settings)
        .unwrap()
        .pow(2);
    let result = execute(&mut network);
    assert_eq!(result.expression, bracket!(&source.expression).pow(2));
    assert!(result.structure.is_scalar());
    assert!(result.proofs.occurrences.get().unwrap().is_empty());
    assert!(result.proofs.owns_expression(result.expression.as_view()));
    assert_eq!(RAW_PARSER_ADMISSIONS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn a_raw_scalar_power_retains_its_first_admission_for_later_powers() {
    test_support::test_initialize();
    let expression = contracted_pair(99784);
    let mut network = SymbolicNet::<AbstractIndex>::from_scalar(expression.clone()).pow(2);
    execute(&mut network);
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));
    let mut network = network.pow(2);
    let result = execute(&mut network);
    assert_eq!(result.expression, bracket!(expression).pow(4));
    assert!(result.proofs.occurrences.get().unwrap().is_empty());
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn opaque_public_parser_observes_missing_leaf_ports_and_enforces_its_settings() {
    test_support::test_initialize();
    let port = mink!(4, 99783);
    let expression = function!(tensor_symbol!("public_opaque_observation"), &port);
    let library = DummyLibrary::<SymbolicTensor>::new();
    let cache = LeafInterfaceCache::default();
    let settings = ParseSettings {
        shorthand_parsing: ShorthandParsing::Opaque,
        depth_limit: Some(1),
        ..ParseSettings::default()
    };
    let unbounded = ParseSettings {
        depth_limit: None,
        ..settings.clone()
    };
    let expanding = ParseSettings {
        shorthand_parsing: ShorthandParsing::expand_all(),
        ..settings.clone()
    };
    let parse = |settings| {
        SymbolicNet::<AbstractIndex>::try_from_admitted_view::<OrderedStructure, _>(
            expression.as_view(),
            &library,
            settings,
            Arc::clone(&cache),
        )
    };
    let network = parse(&settings).unwrap();
    assert_eq!(network.store.tensors[0].structure.order(), 1);
    assert_eq!(cache.lock().unwrap()[&expression].1, [port]);
    assert!(parse(&unbounded).is_err());
    assert!(parse(&expanding).is_err());
}

#[test]
fn admitted_parser_rejects_missing_facts_without_inference() {
    test_support::test_initialize();
    let expression = function!(tensor_symbol!("network_metadata_missing"), mink!(4, 99785));
    let library = DummyLibrary::<SymbolicTensor>::new();
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));
    let error =
        SymbolicNet::<AbstractIndex>::try_from_admitted_view_with_context::<OrderedStructure, _>(
            expression.as_view(),
            &library,
            &ParseSettings::default(),
            LeafInterfaceCache::default(),
            Some(&TensorProofs::default()),
        )
        .unwrap_err();
    assert!(error.to_string().contains("missing its raw parser ports"));
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn a_scalar_root_cannot_supply_missing_constructor_self_loop_ports() {
    test_support::test_initialize();
    let port = mink!(4, 99787);
    let source = SymbolicTensor::from_validated_parts(
        function!(
            tensor_symbol!("network_metadata_missing_loop"),
            &port,
            &port
        ),
        PartialStructure::from_logical_slots([]),
        Arc::default(),
    );
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));
    let error = source
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap_err();
    assert!(error.to_string().contains("missing its raw parser ports"));
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn repeated_admitted_parsing_keeps_the_written_generic_port_order() {
    test_support::test_initialize();
    let logical = vec![mink!(4, 99799), mink!(4, 99791)];
    let expression = FunctionBuilder::new(tensor_symbol!("network_metadata_port_order"))
        .add_args(&logical)
        .finish();
    let source = SymbolicTensor::infer(expression.clone()).unwrap();
    let cache = source.proofs.leaf_interfaces.get().unwrap();
    assert_eq!(cache.lock().unwrap().get(&expression).unwrap().1, logical);
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));
    for _ in 0..2 {
        let mut network = source
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        let result = execute(&mut network);
        assert_eq!(result.expression, expression);
        assert_eq!(cache.lock().unwrap().get(&expression).unwrap().1, logical);
        assert_eq!(
            result
                .proofs
                .leaf_interfaces
                .get()
                .unwrap()
                .lock()
                .unwrap()
                .get(&expression)
                .unwrap()
                .1,
            logical
        );
    }
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn storage_reindexing_transports_permutations_and_repeated_axis_occurrences() {
    test_support::test_initialize();
    let head = tensor_symbol!("network_metadata_reindexed");
    for source_ports in [[99811, 99813], [99811, 99811]] {
        let expression = function!(head, mink!(4, source_ports[0]), mink!(4, source_ports[1]));
        let source =
            SymbolicTensor::structure_from_atom(expression.as_view(), &mut SlotMatcher::default())
                .unwrap()
                .into_canonical();
        INFERENCE_CALLS.with(|count| count.set(0));
        INFERENCE_VISITS.with(|count| count.set(0));
        let result = source
            .reindex_storage(&[AbstractIndex::Normal(99819), AbstractIndex::Normal(99817)])
            .unwrap()
            .apply();
        assert_eq!(
            result.expression,
            function!(head, mink!(4, 99819), mink!(4, 99817))
        );
        assert_eq!(
            result.structure.external_indices_iter().collect::<Vec<_>>(),
            vec![AbstractIndex::Normal(99817), AbstractIndex::Normal(99819)]
        );
        assert_eq!(
            result
                .proofs
                .leaf_interfaces
                .get()
                .unwrap()
                .lock()
                .unwrap()
                .get(&result.expression)
                .unwrap()
                .1,
            vec![mink!(4, 99819), mink!(4, 99817)]
        );
        let closed = result
            .reindex_storage(&[AbstractIndex::Normal(99821), AbstractIndex::Normal(99821)])
            .unwrap()
            .apply()
            .internal_contract();
        assert!(closed.structure.is_scalar());
        assert_eq!(
            closed.expression,
            function!(head, mink!(4, 99821), mink!(4, 99821))
        );
        assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
        assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    }
}

#[test]
fn storage_dual_transports_oriented_ports_and_restores_them_on_second_dual() {
    test_support::test_initialize();
    let rep = LibraryRep::from(ColorFundamental {}).new_rep(3);
    let left = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99831));
    let right = rep
        .dual()
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99833));
    let head = tensor_symbol!("network_metadata_dual");
    let expression = function!(head, left.to_atom(), right.to_atom());
    let source: SymbolicTensor =
        SymbolicTensor::structure_from_atom(expression.as_view(), &mut SlotMatcher::default())
            .unwrap()
            .into_canonical();
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));
    let dual = source.dual();
    assert_eq!(
        dual.expression,
        function!(
            head,
            left.rep()
                .dual()
                .slot::<AbstractIndex, _>(left.aind())
                .to_atom(),
            right
                .rep()
                .dual()
                .slot::<AbstractIndex, _>(right.aind())
                .to_atom()
        )
    );
    assert_eq!(
        dual.proofs
            .occurrences
            .get()
            .unwrap()
            .count(left.rep().dual().slot::<AbstractIndex, _>(left.aind())),
        1
    );
    assert_eq!(dual.dual().expression, expression);
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}
