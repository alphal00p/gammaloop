//! Network execution must preserve the admission owner's scalar metadata.

use super::*;
use crate::representations::ColorFundamental;
use crate::tensor::inference::tests::{INFERENCE_CALLS, INFERENCE_VISITS};
use crate::test_support;
use spenso::network::parsing::{LeafInterfaceCache, ShorthandParsing};
use spenso::structure::representation::RepName;
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
fn admitted_parser_rejects_missing_facts_without_inference() {
    test_support::test_initialize();
    let expression = function!(tensor_symbol!("network_metadata_missing"), mink!(4, 99785));
    let library = DummyLibrary::<SymbolicTensor>::new();
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));
    let error = SymbolicNet::<AbstractIndex>::try_from_admitted_view::<OrderedStructure, _>(
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
