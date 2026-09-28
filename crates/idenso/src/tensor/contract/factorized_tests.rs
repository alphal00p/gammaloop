use crate::{
    shorthands::schoonschip::{ContractionStatus, Schoonschip},
    tensor::SymbolicTensor,
};
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::partial::{PartialStructure, PartialStructureExt},
};
use std::sync::{Arc, Mutex};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    parser::ParseSettings,
};

fn setup() {
    crate::test_support::test_initialize();
    for head in ["p", "q"] {
        SPENSO_TAG.rank_one_tensor_symbol(&format!("fused_raw_routing::{head}"));
    }
    for head in ["t", "u"] {
        SPENSO_TAG.tensor_symbol(&format!("fused_raw_routing::{head}"));
    }
}

fn input(source: &str) -> Atom {
    Atom::parse(source, "fused_raw_routing", ParseSettings::symbolica()).unwrap()
}

fn contract(source: &Atom) -> SymbolicTensor<PartialStructure> {
    SymbolicTensor::<PartialStructure>::infer(source.clone())
        .unwrap()
        .contract(Default::default())
        .unwrap()
        .resolved()
        .unwrap()
}

#[test]
fn factored_contracts_opaque_sum_products() {
    setup();
    let source = input(
        "(t(spenso::mink(4,a),spenso::mink(4,b))+u(spenso::mink(4,a),spenso::mink(4,b)))*(p(spenso::mink(4,a))+q(spenso::mink(4,a)))",
    );
    let expected = input(
        "t(p(spenso::mink(4)),spenso::mink(4,b))+t(q(spenso::mink(4)),spenso::mink(4,b))+u(p(spenso::mink(4)),spenso::mink(4,b))+u(q(spenso::mink(4)),spenso::mink(4,b))",
    );
    assert_eq!(contract(&source).expression.expand(), expected.expand());
    assert_ne!(
        source, expected,
        "success must perform contraction, not only agree after expansion"
    );
    assert_eq!(contract(&expected).expression.expand(), expected.expand());
    assert_eq!(source.expand().schoonschip(), expected);
}

#[test]
fn factored_contract_preserves_common_scalar_spectator_factorization() {
    setup();
    let spectator = input(
        "(x+y)^3*(spenso::g(p(spenso::mink(4)),p(spenso::mink(4)))+spenso::g(q(spenso::mink(4)),q(spenso::mink(4))))",
    );
    let tensor = input(
        "(t(spenso::mink(4,a),spenso::mink(4,b))+u(spenso::mink(4,a),spenso::mink(4,b)))*(p(spenso::mink(4,a))+q(spenso::mink(4,a)))",
    );
    let contracted = input(
        "t(p(spenso::mink(4)),spenso::mink(4,b))+t(q(spenso::mink(4)),spenso::mink(4,b))+u(p(spenso::mink(4)),spenso::mink(4,b))+u(q(spenso::mink(4)),spenso::mink(4,b))",
    );
    let source = &spectator * tensor;
    let expected = &spectator * contracted;
    let actual = contract(&source).expression;
    assert_eq!(
        actual.expand(),
        expected.expand(),
        "resolved definitions retain the exact contracted polynomial"
    );
    assert_ne!(
        actual,
        actual.expand(),
        "this fixture distinguishes factor protection from global expansion"
    );
    assert_eq!(actual.expand(), source.expand().schoonschip().expand());
    assert_eq!(contract(&actual).expression.expand(), actual.expand());
}

#[test]
fn factored_contract_preserves_typed_logical_order_and_zero() {
    setup();
    let product = input(
        "(t(spenso::mink(4,a),spenso::mink(4,b),spenso::mink(6,c))+u(spenso::mink(4,a),spenso::mink(4,b),spenso::mink(6,c)))*(p(spenso::mink(4,a))+q(spenso::mink(4,a)))",
    );
    let expected = input(
        "t(p(spenso::mink(4)),spenso::mink(4,b),spenso::mink(6,c))+t(q(spenso::mink(4)),spenso::mink(4,b),spenso::mink(6,c))+u(p(spenso::mink(4)),spenso::mink(4,b),spenso::mink(6,c))+u(q(spenso::mink(4)),spenso::mink(4,b),spenso::mink(6,c))",
    );
    for (source, expected) in [
        (product.clone(), expected.clone()),
        (product - expected, Atom::Zero),
    ] {
        let inferred = SymbolicTensor::<PartialStructure>::infer(source).unwrap();
        let mut tensor =
            SymbolicTensor::checked_parts(inferred.expression, inferred.structure).unwrap();
        tensor.structure = PartialStructure::from_logical_slots(
            tensor.structure.logical_slots().into_iter().rev(),
        );
        let original_ports = tensor.structure.logical_slots();
        assert_eq!(original_ports.len(), 2);
        let result = tensor
            .contract(Default::default())
            .unwrap()
            .resolved()
            .unwrap();
        assert_eq!(result.expression.expand(), expected.expand());
        assert_eq!(result.structure.logical_slots(), original_ports);
        assert_eq!(
            result
                .contract(Default::default())
                .unwrap()
                .resolved()
                .unwrap()
                .expression
                .expand(),
            result.expression.expand()
        );
    }
}

#[test]
fn factored_contract_coalesces_vector_powers_before_counting_ports() {
    setup();
    let source = input("(p(spenso::mink(4,a))+q(spenso::mink(4,a)))^2");
    let expected = input("spenso::g(p(spenso::mink(4)),p(spenso::mink(4)))+2*spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+spenso::g(q(spenso::mink(4)),q(spenso::mink(4)))").expand();
    let actual = contract(&source).expression;
    assert_eq!(actual.expand(), expected, "{source}");
    assert_eq!(
        actual.expand(),
        source.expand().schoonschip().expand(),
        "normalized duplicate powers retain their local-pair rule"
    );
    assert_eq!(contract(&actual).expression.expand(), actual.expand());
    // Three source occurrences are ambiguous before distribution, even when
    // cancellation would expose an admissible polynomial. Typed admission must
    // reject that source rather than infer its interface after cancellation.
    assert!(
        SymbolicTensor::<PartialStructure>::infer(input(
            "p(spenso::mink(4,a))*(p(spenso::mink(4,a))+q(spenso::mink(4,a)))*(p(spenso::mink(4,a))-q(spenso::mink(4,a)))"
        ))
        .is_err()
    );
}

#[test]
fn factored_contract_does_not_relax_missing_sum_branch_interfaces() {
    setup();
    let source = input("p(spenso::mink(4,a))*(t(spenso::mink(4,a))+u(spenso::mink(4,b)))");
    assert!(SymbolicTensor::<PartialStructure>::infer(source.clone()).is_err());
    assert!(SymbolicTensor::<PartialStructure>::infer(source.expand()).is_err());
    // Unlike the retired raw expansion flag, typed contraction must not
    // distribute first and invent a missing branch interface.
}

#[test]
fn factored_contract_preserves_callback_fallback_and_hidden_metadata_schedule() {
    setup();
    let calls = Arc::new(Mutex::new(Vec::new()));
    let observed = Arc::clone(&calls);
    let head = spenso::tensor_symbol!(
        "fused_raw_callback_leaf",
        norm = move |value, _| observed.lock().unwrap().push(value.to_owned())
    );
    let a = spenso::mink!(4, 78311);
    let b = spenso::mink!(4, 78317);
    let leaf = FunctionBuilder::new(head).add_arg(&a).finish();
    let product = ETS.metric(&a, &b) * leaf;
    let scalar = symbolica::symbol!("fused_raw_routing::opaque_metadata"; Scalar);
    let metadata = FunctionBuilder::new(scalar).add_arg(&product).finish();
    let outside = input("(p(spenso::mink(4,c))+q(spenso::mink(4,c)))*t(spenso::mink(4,c))");
    for source in [product, metadata * outside] {
        // Compare operations on an established carrier. Construction is its
        // own checked boundary and may normalize callback-bearing leaves.
        let typed = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
        calls.lock().unwrap().clear();
        let expected = source.schoonschip();
        let transcript = calls.lock().unwrap().clone();
        assert!(!transcript.is_empty());
        calls.lock().unwrap().clear();
        assert_eq!(
            typed
                .contract(Default::default())
                .unwrap()
                .expanded()
                .unwrap()
                .expression,
            expected.expand()
        );
        assert_eq!(
            *calls.lock().unwrap(),
            transcript,
            "declined planning must not replay normalization"
        );
        let typed = SymbolicTensor::<PartialStructure>::infer(expected.clone()).unwrap();
        calls.lock().unwrap().clear();
        assert_eq!(
            typed
                .contract(Default::default())
                .unwrap()
                .expanded()
                .unwrap()
                .expression,
            expected.expand()
        );
        assert!(calls.lock().unwrap().is_empty());
    }
}

#[test]
fn factored_contract_keeps_checked_metric_callback_rank_loss_rejection() {
    setup();
    let calls = Arc::new(Mutex::new(Vec::new()));
    let observed = Arc::clone(&calls);
    let a = spenso::mink!(4, 78401);
    let b = spenso::mink!(4, 78403);
    let target = b.clone();
    let head = spenso::tensor_symbol!(
        "fused_raw_metric_rank_loss",
        norm = move |value, output| {
            observed.lock().unwrap().push(value.to_owned());
            if let AtomView::Fun(function) = value
                && function.iter().next() == Some(target.as_view())
            {
                **output = Atom::num(7);
            }
        }
    );
    let source = ETS.metric(&a, &b) * FunctionBuilder::new(head).add_arg(&a).finish();
    let inferred = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
    let tensor = SymbolicTensor::checked_parts(inferred.expression, inferred.structure).unwrap();
    calls.lock().unwrap().clear();
    let expected = source.schoonschip();
    let transcript = calls.lock().unwrap().clone();
    assert_eq!(expected, Atom::num(7));
    assert!(!transcript.is_empty());
    calls.lock().unwrap().clear();
    assert!(tensor.contract(Default::default()).is_err());
    assert_eq!(
        *calls.lock().unwrap(),
        transcript,
        "checked contraction must not replay callbacks while rejecting rank loss"
    );
}

#[test]
fn factored_contract_preserves_vector_metadata_in_merged_states() {
    setup();
    let source = input(
        "(p(1,spenso::mink(4,a))+p(2,spenso::mink(4,a)))*(q(3,label,spenso::mink(4,a))+q(4,label,spenso::mink(4,a)))",
    );
    let expected = source.expand().schoonschip().expand();
    let typed = SymbolicTensor::<PartialStructure>::infer(source).unwrap();
    let result = typed.contract(Default::default()).unwrap();
    assert_eq!(result.expanded().unwrap().expression, expected);
    assert!(!result.expression.get_aliases().is_empty());
    assert_eq!(expected.nterms(), 4);
    assert_eq!(
        result
            .resolved()
            .unwrap()
            .contract(Default::default())
            .unwrap()
            .expanded()
            .unwrap()
            .expression,
        expected
    );
}

#[test]
fn factor_occurrences_follow_transparent_graph_products() {
    setup();
    let body =
        input("(p(spenso::mink(4,a))+q(spenso::mink(4,a)))*t(spenso::mink(4,a),spenso::mink(4,b))");
    let source = input(
        "spenso::bracket(p(spenso::mink(4,a))+q(spenso::mink(4,a)),t(spenso::mink(4,a),spenso::mink(4,b)))",
    );
    let expected = body.expand().schoonschip().expand();
    let typed = SymbolicTensor::<PartialStructure>::infer(source).unwrap();
    assert!(
        crate::shorthands::schoonschip::SlotContraction::new()
            .contract_factorized(typed.expression.as_view(), None, true)
            .is_some(),
        "the parsed transparent product must use graph intake"
    );
    let result = typed.contract(Default::default()).unwrap();
    assert_eq!(result.expanded().unwrap().expression, expected);
    assert_eq!(result.root().structure, typed.structure);
}

#[test]
fn powered_scalar_vector_sum_preserves_dummy_scope() {
    setup();
    let base = input("p(spenso::mink(4,a))*q(spenso::mink(4,a))+x");
    let expected_base = input("spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+x");
    for power in [2, 3] {
        let source = base.pow(power);
        let expected = expected_base.pow(power);
        let typed = SymbolicTensor::infer(source.clone()).unwrap();
        let collected = typed.contract(Default::default()).unwrap();
        assert!(
            collected.contraction_complete(),
            "power scope must complete"
        );
        let actual = collected.resolved().unwrap().expression;
        let raw = source.schoonschip();
        assert_eq!(raw, expected, "legacy raw control");
        assert_eq!(
            actual, expected,
            "each copy owns its internal dummy; source={source}; legacy raw={raw}"
        );
        assert_eq!(contract(&actual).expression, actual);
    }
}

#[test]
fn powered_scalar_metric_path_sum_preserves_dummy_scope() {
    setup();
    let base = input(
        "spenso::g(spenso::mink(4,a),spenso::mink(4,b))*p(spenso::mink(4,a))*q(spenso::mink(4,b))+x",
    );
    let expected_base = input("spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+x");
    for power in [2, 3] {
        let source = base.pow(power);
        let expected = expected_base.pow(power);
        let typed = SymbolicTensor::infer(source.clone()).unwrap();
        let collected = typed.contract(Default::default()).unwrap();
        assert!(
            collected.contraction_complete(),
            "power scope must complete"
        );
        let actual = collected.resolved().unwrap().expression;
        let raw = source.schoonschip();
        assert_eq!(raw, expected, "legacy raw control");
        assert_eq!(
            actual, expected,
            "each metric path copy owns its internal dummies; source={source}; legacy raw={raw}"
        );
        assert_eq!(contract(&actual).expression, actual);
    }
}

#[test]
fn powered_scalar_metric_cycle_sum_preserves_dummy_scope() {
    setup();
    let base = input(
        "spenso::g(spenso::mink(4,a),spenso::mink(4,b))*spenso::g(spenso::mink(4,b),spenso::mink(4,c))*spenso::g(spenso::mink(4,c),spenso::mink(4,a))+x",
    );
    let expected_base = input("4+x");
    for power in [2, 3] {
        let source = base.pow(power);
        let expected = expected_base.pow(power);
        let typed = SymbolicTensor::infer(source.clone()).unwrap();
        let collected = typed.contract(Default::default()).unwrap();
        assert!(
            collected.contraction_complete(),
            "power scope must complete"
        );
        let actual = collected.resolved().unwrap().expression;
        let raw = source.schoonschip();
        assert_eq!(raw, expected, "legacy raw control");
        assert_eq!(
            actual, expected,
            "each metric cycle copy owns its internal dummies; source={source}; legacy raw={raw}"
        );
        let rerun =
            SymbolicTensor::checked_parts(actual.clone(), PartialStructure::from_logical_slots([]))
                .unwrap()
                .contract(Default::default())
                .unwrap()
                .resolved()
                .unwrap();
        assert_eq!(rerun.expression, actual);
    }
}

#[test]
fn powered_scalar_sum_keeps_disjoint_spectators_factored() {
    setup();
    let source =
        input("(p(spenso::mink(4,a))*q(spenso::mink(4,a))+x)^2*(y+z)^3*t(spenso::mink(6,c))");
    let expected = input(
        "(spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+x)^2*(y+z)^3*t(spenso::mink(6,c))",
    );
    let actual = contract(&source);
    assert_eq!(actual.expression, expected);
    assert_eq!(actual.structure.logical_slots().len(), 1);
    assert_eq!(contract(&actual.expression).expression, actual.expression);
}

#[test]
fn powered_open_metric_path_keeps_internal_and_external_pairs_distinct() {
    setup();
    let source = input(
        "(spenso::g(spenso::mink(4,a),spenso::mink(4,b))*p(spenso::mink(4,a))+x*q(spenso::mink(4,b)))^2",
    );
    let expected = input(
        "spenso::g(p(spenso::mink(4)),p(spenso::mink(4)))+2*x*spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+x^2*spenso::g(q(spenso::mink(4)),q(spenso::mink(4)))",
    );
    let typed = SymbolicTensor::infer(source.clone()).unwrap();
    assert_eq!(typed.expression, source, "constructor control");
    let direct = crate::shorthands::schoonschip::SlotContraction::new()
        .contract_factorized(source.as_view(), None, true)
        .expect("open power must be admitted directly, before fallback rewriting");
    assert!(
        direct.status == ContractionStatus::Complete,
        "direct graph contraction must complete"
    );
    let mut direct_expression = symbolica::atom::AliasedAtom::from(direct.root);
    for (handle, body) in direct.aliases {
        direct_expression.register_alias(handle, body);
    }
    assert_eq!(direct_expression.into_inner().expand(), expected.expand());
    let collected = typed.contract(Default::default()).unwrap();
    assert!(
        collected.contraction_complete(),
        "open power must complete through the graph owner, not silently refuse"
    );
    let actual = collected.resolved().unwrap();
    let raw = source.schoonschip();
    assert_eq!(
        actual.expression.expand(),
        expected.expand(),
        "internal a is local to each copy; only b contracts across copies; legacy raw={raw}"
    );
    assert_eq!(contract(&actual.expression).expression, actual.expression);
}

#[test]
fn powered_open_vector_base_preserves_internal_scalar_pair() {
    setup();
    let source = input(
        "(p(spenso::mink(4,a))*q(spenso::mink(4,a))*p(spenso::mink(4,b))+x*q(spenso::mink(4,b)))^2",
    );
    let expected = input(
        "spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))^2*spenso::g(p(spenso::mink(4)),p(spenso::mink(4)))+2*x*spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))^2+x^2*spenso::g(q(spenso::mink(4)),q(spenso::mink(4)))",
    );
    let typed = SymbolicTensor::infer(source.clone()).unwrap();
    assert_eq!(typed.expression, source, "constructor control");
    let direct = crate::shorthands::schoonschip::SlotContraction::new()
        .contract_factorized(source.as_view(), None, true)
        .expect("open internal pair must use the graph owner");
    assert!(direct.status == ContractionStatus::Complete);
    let mut direct_expression = symbolica::atom::AliasedAtom::from(direct.root);
    for (handle, body) in direct.aliases {
        direct_expression.register_alias(handle, body);
    }
    assert_eq!(direct_expression.into_inner().expand(), expected.expand());
    let actual = typed
        .contract(Default::default())
        .unwrap()
        .resolved()
        .unwrap();
    assert_eq!(actual.expression.expand(), expected.expand());
    assert_eq!(contract(&actual.expression).expression, actual.expression);
}

#[test]
fn budget_limited_alias_contraction_reports_incomplete_without_reentering_its_frontier() {
    setup();
    let slots = (0..9)
        .map(|i| format!("spenso::mink(4,budget_port_{i})"))
        .collect::<Vec<_>>();
    let expression = input(&format!(
        "t({})*{}",
        slots.join(","),
        slots
            .iter()
            .map(|slot| format!("(x*p({slot})+y*q({slot}))"))
            .collect::<Vec<_>>()
            .join("*")
    ));
    let body = SymbolicTensor::infer(expression).unwrap();
    let direct = body.contract(Default::default()).unwrap();
    assert!(
        !direct.contraction_complete(),
        "the fixed test frontier budget must expose its incomplete status"
    );
    let handle = body.alias_handle().unwrap();
    let source = Arc::new(handle.clone().with_aliases([(handle, body)]).unwrap());
    assert!(
        !source.contraction_complete(),
        "attaching a registry alone cannot certify contraction"
    );
    let result = source.contract(Default::default()).unwrap();
    assert!(!result.contraction_complete());
    assert_eq!(
        result.expression.get_aliases().len(),
        direct.expression.get_aliases().len() + 1,
        "a single call must retain the bounded frontier once, not recursively restart generated templates"
    );
    assert_eq!(
        result.resolved().unwrap().expression.expand(),
        direct.resolved().unwrap().expression.expand()
    );
}
