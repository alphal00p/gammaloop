use crate::{
    shorthands::schoonschip::{ReductionStatus, Schoonschip},
    tensor::SymbolicTensor,
};
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        TensorStructure,
        partial::{PartialStructure, PartialStructureExt},
    },
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
}

#[test]
fn disabled_metric_elimination_preserves_supplied_vectors() {
    use crate::tensor::ContractSettings;
    setup();
    let settings = ContractSettings {
        metrics: false,
        ..Default::default()
    };
    for (source, expected) in [
        (
            "spenso::g(spenso::mink(4,mu),spenso::mink(4,nu))*p(spenso::mink(4,mu))",
            "spenso::g(p(spenso::mink(4)),spenso::mink(4,nu))",
        ),
        (
            "spenso::g(spenso::mink(4,mu),spenso::mink(4,nu))*p(spenso::mink(4,mu))*q(spenso::mink(4,nu))",
            "spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))",
        ),
        (
            "(x+y)^7*spenso::g(spenso::mink(4,mu),spenso::mink(4,nu))*p(spenso::mink(4,mu))*q(spenso::mink(4,nu))",
            "(x+y)^7*spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))",
        ),
    ] {
        let source = SymbolicTensor::infer(input(source)).unwrap();
        let result = source.contract(settings).unwrap();
        assert_eq!(result.expression, input(expected));
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.contract(settings).unwrap(), result);
    }
}

#[test]
fn contraction_forwards_observations_after_normalization_and_ordered_rewriting() {
    use crate::tensor::{ContractSettings, SELECTED_OPENS};
    setup();
    let settings = ContractSettings {
        collect_chains: false,
        collect_traces: false,
        ..Default::default()
    };
    for (expression, normalized_first) in [
        // The empty chain becomes a metric before planning, changing the first
        // branch key while retaining its scalar coefficient and tensor leaves.
        (
            "(x+y)*spenso::chain(spenso::mink(4,a),spenso::mink(4,b))*t(spenso::mink(4,a))*u(spenso::mink(4,b))-(x+y)*t(spenso::mink(4,c))*u(spenso::mink(4,c))",
            true,
        ),
        // A tensor argument with pending contractions initially declines the
        // polynomial intake. The ordered rewrite finishes its metadata; the
        // retry must still collect the two alpha-equivalent dummy networks.
        (
            "t(p(spenso::mink(4,x))*q(spenso::mink(4,x)),spenso::mink(4,a))*u(spenso::mink(4,a))-t(spenso::g(p(spenso::mink(4)),q(spenso::mink(4))),spenso::mink(4,b))*u(spenso::mink(4,b))",
            false,
        ),
    ] {
        let source = SymbolicTensor::infer(input(expression)).unwrap();
        assert!(!source.expression.is_zero());
        assert_eq!(
            source.normalized_contract_notation(settings) != source.expression,
            normalized_first,
        );
        source.reduction_observations();
        SELECTED_OPENS.with(|count| count.set(0));
        let result = source.contract(settings).unwrap();
        assert_eq!(result.expression, Atom::Zero);
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(SELECTED_OPENS.with(|count| count.get()), 0);
        assert_eq!(result.contract(settings).unwrap(), result);
    }
}

#[test]
fn chain_admission_skips_finished_open_tensor_coefficients() {
    use crate::Cookable;
    use crate::shorthands::chain::CHAINIFY_CALLS;
    use crate::tensor::simplification::observation::INITIAL_SCANS;
    setup();
    let reps = crate::test_support::test_initialize();
    let label = spenso::index_symbol!("fused_raw_routing::chain_channel_edge");
    let [a, b] = [1, 2].map(|i| reps.mink4.to_symbolic([symbolica::function!(label, i)]));
    let coefficient = Atom::add_many((0..128).map(|i| {
        let x = Atom::var(symbolica::symbol!(&format!(
            "fused_raw_routing::channel_coefficient_{i}"
        )));
        x * input("spenso::g(p(17,spenso::mink(4)),q(23,spenso::mink(4)))").pow(i % 4 + 1)
    }));
    let vector = |name, slot: &Atom| {
        FunctionBuilder::new(SPENSO_TAG.rank_one_tensor_symbol(name))
            .add_arg(slot)
            .finish()
    };
    let source = SymbolicTensor::infer(
        (&coefficient * spenso::g!(&a, &b)
            + coefficient.pow(2)
                * vector("fused_raw_routing::p", &a)
                * vector("fused_raw_routing::q", &b))
        .cook_indices(),
    )
    .unwrap();
    source.reduction_observations();
    let scans = INITIAL_SCANS.with(|count| count.get());
    CHAINIFY_CALLS.with(|count| count.set(0));
    let result = source.contract(Default::default()).unwrap();
    assert_eq!(result.expression, source.expression);
    assert_eq!(result.structure, source.structure);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(CHAINIFY_CALLS.with(|count| count.get()), 0);
    assert_eq!(INITIAL_SCANS.with(|count| count.get()), scans);
}

#[test]
fn chain_admission_updates_channels_after_vector_substitution() {
    use crate::shorthands::chain::CHAINIFY_CALLS;
    setup();
    let source = input(
        "(t(spenso::mink(4,a),spenso::mink(4,i))+u(spenso::mink(4,a),spenso::mink(4,i)))*p(spenso::mink(4,i))*q(spenso::mink(4,b))",
    );
    let expected = input(
        "(t(spenso::mink(4,a),p(spenso::mink(4)))+u(spenso::mink(4,a),p(spenso::mink(4))))*q(spenso::mink(4,b))",
    );
    CHAINIFY_CALLS.with(|count| count.set(0));
    let result = contract(&source);
    assert_eq!(result.expression, expected);
    assert_eq!(CHAINIFY_CALLS.with(|count| count.get()), 0);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
}

#[test]
fn chain_admission_retains_generic_channels_inside_scalar_scopes() {
    use crate::{shorthands::chain::CHAINIFY_CALLS, tensor::ContractSettings};
    setup();
    let source = SymbolicTensor::infer(input("(x+y)^7*(t(spenso::mink(4,i),spenso::mink(4,j))*u(spenso::mink(4,j),spenso::mink(4,i))+z)^2")).unwrap();
    CHAINIFY_CALLS.with(|count| count.set(0));
    let untraced = source
        .contract(ContractSettings {
            collect_traces: false,
            ..Default::default()
        })
        .unwrap();
    assert!(!untraced.expression.contains_symbol(SPENSO_TAG.trace));
    assert!(CHAINIFY_CALLS.with(|count| count.get()) > 0);
    let traced = untraced.contract(Default::default()).unwrap();
    assert!(traced.expression.contains_symbol(SPENSO_TAG.trace));
    assert_eq!(traced.reduction_status(), ReductionStatus::Complete);
    assert_eq!(traced.structure, source.structure);
    let expected = input(
        "(x+y)^7*(spenso::trace(spenso::mink(4),spenso::cyclic(t(spenso::in,spenso::out),u(spenso::in,spenso::out)))+z)^2",
    );
    assert_eq!(traced.expression, expected);
}

#[test]
fn contraction_collects_matrix_words_with_supplied_endpoints() {
    use crate::tensor::ContractSettings;
    setup();
    for (source, expected, reversed, rank) in [
        (
            "t(p(spenso::mink(4)),spenso::mink(4,b))*u(spenso::mink(4,b),spenso::mink(4,c))",
            "spenso::chain(p(spenso::mink(4)),spenso::mink(4,c),t(spenso::in,spenso::out),u(spenso::in,spenso::out))",
            None,
            1,
        ),
        (
            "t(spenso::mink(4,a),spenso::mink(4,b))*u(spenso::mink(4,b),q(spenso::mink(4)))",
            "spenso::chain(spenso::mink(4,a),q(spenso::mink(4)),t(spenso::in,spenso::out),u(spenso::in,spenso::out))",
            None,
            1,
        ),
        (
            "t(p(spenso::mink(4)),spenso::mink(4,b))*u(spenso::mink(4,b),q(spenso::mink(4)))",
            "spenso::chain(p(spenso::mink(4)),q(spenso::mink(4)),t(spenso::in,spenso::out),u(spenso::in,spenso::out))",
            None,
            0,
        ),
        (
            "t(p(spenso::mink(4)),spenso::mink(4,b))*u(spenso::mink(4,c),spenso::mink(4,b))",
            "spenso::chain(p(spenso::mink(4)),spenso::mink(4,c),t(spenso::in,spenso::out),u(spenso::out,spenso::in))",
            Some(
                "spenso::chain(spenso::mink(4,c),p(spenso::mink(4)),u(spenso::in,spenso::out),t(spenso::out,spenso::in))",
            ),
            1,
        ),
        (
            "t(spenso::mink(4,b),p(spenso::mink(4)))*u(spenso::mink(4,b),spenso::mink(4,c))",
            "spenso::chain(p(spenso::mink(4)),spenso::mink(4,c),t(spenso::out,spenso::in),u(spenso::in,spenso::out))",
            Some(
                "spenso::chain(spenso::mink(4,c),p(spenso::mink(4)),u(spenso::out,spenso::in),t(spenso::in,spenso::out))",
            ),
            1,
        ),
    ] {
        let source = SymbolicTensor::infer(input(source)).unwrap();
        let result = source.contract(Default::default()).unwrap();
        // Self-dual common-start/end contractions may traverse the word in
        // either direction. Reversing the word must also transpose its factors.
        assert!(
            result.expression == input(expected)
                || reversed.is_some_and(|expected| result.expression == input(expected)),
            "unexpected chain orientation: {}",
            result.expression,
        );
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.structure.canonical().order(), rank);
        assert_eq!(result.contract(Default::default()).unwrap(), result);
        let indexed = source
            .contract(ContractSettings {
                collect_chains: false,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(indexed.expression, source.expression);
        assert_eq!(indexed.structure, source.structure);
    }
}

#[test]
fn supplied_chain_endpoints_are_validated_and_do_not_expose_ports() {
    setup();
    for (expression, rank) in [
        (
            "spenso::chain(p(spenso::mink(4)),spenso::mink(4,c),t(spenso::in,spenso::out))",
            1,
        ),
        (
            "spenso::chain(p(spenso::mink(4)),q(spenso::mink(4)),t(spenso::in,spenso::out))",
            0,
        ),
    ] {
        let tensor = SymbolicTensor::infer(input(expression)).unwrap();
        assert_eq!(tensor.structure.canonical().order(), rank);
    }
    for expression in [
        "spenso::chain(p(spenso::mink(3)),spenso::mink(4,c),t(spenso::in,spenso::out))",
        "spenso::chain(p(spenso::mink(4,a)),spenso::mink(4,c),t(spenso::in,spenso::out))",
        "spenso::chain(t(spenso::mink(4)),spenso::mink(4,c),t(spenso::in,spenso::out))",
        "spenso::chain(p(t(spenso::mink(4,a)),spenso::mink(4)),spenso::mink(4,c),t(spenso::in,spenso::out))",
    ] {
        assert!(
            SymbolicTensor::infer(input(expression)).is_err(),
            "{expression}"
        );
    }
}

#[test]
fn contraction_reuses_source_free_color_invariants_without_opening_scopes() {
    use crate::tensor::{ContractSettings, SELECTED_OPENS, collection::SHALLOW_BUILDS};
    let reps = crate::test_support::test_initialize();
    let [a, b, c, d, e, f] = [97601, 97602, 97603, 97604, 97605, 97606]
        .map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
    let rep = reps.cof_nc.to_symbolic([] as [Atom; 0]);
    let symmetric = |x: &Atom| {
        spenso::trace_sym!(
            &rep,
            crate::color_t!(x),
            crate::color_t!(&c),
            crate::color_t!(&d)
        )
    };
    let invariant =
        crate::color_f!(&a, &b, &e) * symmetric(&e) * crate::color_f!(&a, &b, &f) * symmetric(&f);
    let weight = symbolica::parse_lit!((source_free_weight_x + source_free_weight_y) ^ 12);
    let source = SymbolicTensor::infer(weight * (Atom::one() + invariant).pow(2)).unwrap();
    let settings = ContractSettings {
        collect_chains: false,
        collect_traces: false,
        expand: false,
        ..Default::default()
    };
    let observations = source.reduction_observations();
    assert!(observations.candidates.repeated_indices);
    assert!(observations.excludes_contraction_sources(settings));
    SHALLOW_BUILDS.with(|count| count.set(0));
    SELECTED_OPENS.with(|count| count.set(0));
    // Minimal contraction keeps structurally finished invariants opaque.
    // Expanding contraction also collects alpha-equivalent terms and therefore
    // cannot rule out its work solely from the absence of metric/vector sources.
    let result = source.contract_parts(settings).unwrap();
    assert_eq!(result.root, source);
    assert_eq!(result.status, ReductionStatus::Complete);
    assert!(Arc::ptr_eq(
        observations,
        result.root.reduction_observations()
    ));
    assert_eq!(SHALLOW_BUILDS.with(|count| count.get()), 0);
    assert_eq!(SELECTED_OPENS.with(|count| count.get()), 0);
    assert!(
        source
            .contract_parts(ContractSettings {
                order: Some(&[]),
                ..settings
            })
            .is_err()
    );
}

#[test]
fn contraction_source_proof_keeps_new_metric_connections_eligible() {
    use crate::tensor::ContractSettings;
    let reps = crate::test_support::test_initialize();
    let [a, b, c, d, e, z] = [97611, 97612, 97613, 97614, 97615, 97616]
        .map(|index| reps.coad_da.to_symbolic([Atom::num(index)]));
    let rep = reps.cof_nc.to_symbolic([] as [Atom; 0]);
    let symmetric = spenso::trace_sym!(
        &rep,
        crate::color_t!(&e),
        crate::color_t!(&c),
        crate::color_t!(&d)
    );
    let source =
        SymbolicTensor::infer(spenso::g!(&a, &z) * crate::color_f!(&a, &b, &e) * &symmetric)
            .unwrap();
    let settings = ContractSettings {
        collect_chains: false,
        collect_traces: false,
        ..Default::default()
    };
    assert!(
        !source
            .reduction_observations()
            .excludes_contraction_sources(settings)
    );
    let result = source.contract_parts(settings).unwrap();
    assert_eq!(
        result.root.expression,
        crate::color_f!(&z, &b, &e) * symmetric
    );
    assert_eq!(result.root.structure, source.structure);
    assert_eq!(result.status, ReductionStatus::Complete);
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
        "the factorized result retains the exact contracted polynomial"
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
fn disconnected_antisymmetric_spectators_preserve_sum_contractions() {
    setup();
    let pair = input(
        "(p(spenso::mink(4,mu))+q(spenso::mink(4,mu)))*(p(17,spenso::mink(4,mu))-q(23,spenso::mink(4,mu)))",
    );
    let contracted = input(
        "spenso::g(p(spenso::mink(4)),p(17,spenso::mink(4)))-spenso::g(p(spenso::mink(4)),q(23,spenso::mink(4)))+spenso::g(q(spenso::mink(4)),p(17,spenso::mink(4)))-spenso::g(q(spenso::mink(4)),q(23,spenso::mink(4)))",
    );
    let color = input("spenso::f(spenso::coad(8,a),spenso::coad(8,b),spenso::coad(8,c))");
    let coefficient = input("(x+y)^3");
    let settings = crate::tensor::ContractSettings {
        collect_chains: false,
        collect_traces: false,
        ..Default::default()
    };
    for power in [1, 2] {
        let spectator = &coefficient * color.pow(power);
        let source = SymbolicTensor::infer(&pair * &spectator).unwrap();
        let expected = &contracted * spectator;
        let direct = crate::shorthands::schoonschip::SlotContraction::new()
            .contract_factorized(source.expression.as_view(), None, true)
            .expect("an unsupported component must not discard independent work");
        assert_eq!(direct.root, expected);
        assert!(direct.observations.is_none(), "colour indices still remain");
        let result = source.contract(settings).unwrap();
        assert_eq!(result.expression, expected);
        assert_eq!(result.structure, source.structure);
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.contract(settings).unwrap().expression, expected);
    }
}

#[test]
fn disconnected_refusal_retains_local_metric_fallback_and_finished_components() {
    setup();
    let source = input(
        "(x+y)^3*(3*p(spenso::mink(4,mu))+2*q(spenso::mink(4,mu)))*q(spenso::mink(4,mu))*p(11,spenso::mink(4,nu))*q(13,spenso::mink(4,nu))*spenso::g(spenso::coad(8,a),spenso::coad(8,d))*spenso::f(spenso::coad(8,a),spenso::coad(8,b),spenso::coad(8,c))",
    );
    let expected = input(
        "(x+y)^3*(3*spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+2*spenso::g(q(spenso::mink(4)),q(spenso::mink(4))))*spenso::g(p(11,spenso::mink(4)),q(13,spenso::mink(4)))*spenso::f(spenso::coad(8,d),spenso::coad(8,b),spenso::coad(8,c))",
    );
    let typed = SymbolicTensor::infer(source.clone()).unwrap();
    let AtomView::Mul(product) = source.as_view() else {
        unreachable!();
    };
    let count = product.iter().len();
    for order in [(0..count).collect::<Vec<_>>(), (0..count).rev().collect()] {
        let direct = crate::shorthands::schoonschip::SlotContraction::new()
            .contract_factorized(source.as_view(), Some(&order), true)
            .expect("component refusal retains an exact partial result");
        assert_eq!(direct.root, expected);
        assert!(direct.observations.is_none());
        let result = typed
            .contract(crate::tensor::ContractSettings {
                order: Some(&order),
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(result.expression, expected);
        assert_eq!(result.structure, typed.structure);
    }
}

#[test]
fn disconnected_components_do_not_bypass_callback_interface_checks() {
    setup();
    let a = spenso::mink!(4, 78901);
    let b = spenso::mink!(4, 78903);
    let target = b.clone();
    let calls = Arc::new(Mutex::new(Vec::new()));
    let observed = Arc::clone(&calls);
    let head = spenso::tensor_symbol!(
        "fused_raw_disconnected_callback_rank_loss",
        norm = move |value, output| {
            observed.lock().unwrap().push(value.to_owned());
            if let AtomView::Fun(function) = value
                && function.iter().next() == Some(target.as_view())
            {
                **output = Atom::num(7);
            }
        }
    );
    let source = ETS.metric(&a, &b)
        * FunctionBuilder::new(head).add_arg(&a).finish()
        * input(
            "(p(spenso::mink(4,mu))+q(spenso::mink(4,mu)))*q(spenso::mink(4,mu))*spenso::f(spenso::coad(8,a),spenso::coad(8,b),spenso::coad(8,c))^2",
        );
    let typed = SymbolicTensor::infer(source.clone()).unwrap();
    calls.lock().unwrap().clear();
    let _ = source.schoonschip();
    let expected_calls = calls.lock().unwrap().clone();
    assert!(!expected_calls.is_empty());
    calls.lock().unwrap().clear();
    assert!(
        typed
            .contract(crate::tensor::ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            })
            .is_err()
    );
    assert_eq!(*calls.lock().unwrap(), expected_calls);
}

#[test]
fn factored_contract_keeps_completed_dot_coefficients_opaque() {
    setup();
    let coefficient = input(
        &(0..256)
            .map(|i| format!("x{i}*spenso::g(p({i},spenso::mink(4)),q({i},spenso::mink(4)))"))
            .collect::<Vec<_>>()
            .join("+"),
    );
    let tensor = input("(p(spenso::mink(4,a))+q(spenso::mink(4,a)))*q(spenso::mink(4,a))");
    crate::tensor::SELECTED_OPENS.with(|count| count.set(0));
    let baseline = contract(&tensor);
    let baseline_opens = crate::tensor::SELECTED_OPENS.with(|count| count.get());
    let source = &coefficient * tensor;
    let expected = &coefficient
        * input(
            "spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+spenso::g(q(spenso::mink(4)),q(spenso::mink(4)))",
        );
    crate::tensor::SELECTED_OPENS.with(|count| count.set(0));
    let actual = contract(&source);
    assert_eq!(
        crate::tensor::SELECTED_OPENS.with(|count| count.get()),
        baseline_opens,
        "the completed scalar coefficient must remain unopened"
    );
    assert_eq!(actual.expression, &coefficient * baseline.expression);
    assert_eq!(actual.reduction_status(), ReductionStatus::Complete);
    assert_eq!(actual.expression, expected);
    assert_eq!(contract(&actual.expression).expression, actual.expression);

    // Scalar interfaces can still hide genuine index contractions. Keep the
    // large finished coefficient closed while opening the internal dummy pair.
    let internal = input("p(spenso::mink(4,b))*q(spenso::mink(4,b))+z");
    let source = &coefficient * internal;
    let expected = coefficient * input("spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+z");
    let actual = contract(&source);
    assert_eq!(actual.reduction_status(), ReductionStatus::Complete);
    assert_eq!(actual.expression, expected);
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
        let result = tensor.contract(Default::default()).unwrap();
        assert_eq!(result.expression.expand(), expected.expand());
        assert_eq!(result.structure.logical_slots(), original_ports);
        assert_eq!(
            result
                .contract(Default::default())
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
                .expression
                .expand(),
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
                .expression
                .expand(),
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
    assert_eq!(result.expression.expand(), expected);
    assert_eq!(result.structure, typed.structure);
    assert_eq!(expected.nterms(), 4);
    assert_eq!(
        result
            .contract(Default::default())
            .unwrap()
            .expression
            .expand(),
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
    assert_eq!(result.expression.expand(), expected);
    assert_eq!(result.structure, typed.structure);
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
        let actual = collected.expression;
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
        let actual = collected.expression;
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
        let actual = collected.expression;
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
        direct.status == ReductionStatus::Complete,
        "direct graph contraction must complete"
    );
    assert_eq!(direct.root.expand(), expected.expand());
    let collected = typed.contract(Default::default()).unwrap();
    assert!(
        collected.contraction_complete(),
        "open power must complete through the graph owner, not silently refuse"
    );
    let actual = collected;
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
    assert!(direct.status == ReductionStatus::Complete);
    assert_eq!(direct.root.expand(), expected.expand());
    let actual = typed.contract(Default::default()).unwrap();
    assert_eq!(actual.expression.expand(), expected.expand());
    assert_eq!(contract(&actual.expression).expression, actual.expression);
}

#[test]
fn default_contraction_completes_the_selected_product() {
    setup();
    let slots = (0..9)
        .map(|i| format!("spenso::mink(4,product_port_{i})"))
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
    assert!(direct.contraction_complete());
    assert_eq!(direct.reduction_status(), ReductionStatus::Complete);
    assert_eq!(direct.structure, body.structure);
    assert_eq!(
        direct.expression.expand(),
        body.expression.expand().schoonschip().expand()
    );
    assert_eq!(direct.contract(Default::default()).unwrap(), direct);
}
