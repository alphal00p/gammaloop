use crate::tensor::{
    AlgebraContraction, AlgebraSettings, ContractSettings, ReductionStatus, SymbolicTensor,
};
use spenso::{network::tags::SPENSO_TAG, structure::partial::PartialStructure};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    parser::ParseSettings,
};

fn input(value: &str) -> Atom {
    crate::test_support::test_initialize();
    for name in ["p", "q", "r", "s"] {
        SPENSO_TAG.rank_one_tensor_symbol(&format!("minimal_contract::{name}"));
    }
    for name in ["A", "B", "C"] {
        SPENSO_TAG.tensor_symbol(&format!("minimal_contract::{name}"));
    }
    Atom::parse(value, "minimal_contract", ParseSettings::symbolica()).unwrap()
}

fn tensor(value: &str) -> SymbolicTensor<PartialStructure> {
    SymbolicTensor::infer(input(value)).unwrap()
}

#[test]
fn minimal_contract_substitutes_across_existing_sum_branches() {
    for (source, expected) in [
        (
            "spenso::g(spenso::mink(4,a),spenso::mink(4,b))*(p(spenso::mink(4,a))+q(spenso::mink(4,a)))",
            "p(spenso::mink(4,b))+q(spenso::mink(4,b))",
        ),
        (
            "p(spenso::mink(4,a))*(q(spenso::mink(4,a))+r(spenso::mink(4,a)))",
            "spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+spenso::g(p(spenso::mink(4)),r(spenso::mink(4)))",
        ),
        (
            "(x*spenso::g(spenso::mink(4,a),spenso::mink(4,b))+p(spenso::mink(4,a))*q(spenso::mink(4,b))+B(spenso::mink(4,a),spenso::mink(4,b)))*A(spenso::mink(4,a))",
            "x*A(spenso::mink(4,b))+A(p(spenso::mink(4)))*q(spenso::mink(4,b))+B(spenso::mink(4,a),spenso::mink(4,b))*A(spenso::mink(4,a))",
        ),
        (
            "(x+y)^7*(spenso::g(spenso::mink(4,a),spenso::mink(4,b))*p(spenso::mink(4,a))+q(spenso::mink(4,b)))",
            "(x+y)^7*(p(spenso::mink(4,b))+q(spenso::mink(4,b)))",
        ),
        (
            "(p(spenso::mink(4,a))*q(spenso::mink(4,a))+x)^3*(A(spenso::mink(4,b))+B(spenso::mink(4,b)))",
            "(spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+x)^3*(A(spenso::mink(4,b))+B(spenso::mink(4,b)))",
        ),
    ] {
        let source = tensor(source);
        let result = source
            .contract(ContractSettings {
                expand: false,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(result.expression, input(expected));
        assert_eq!(result.reduction_status(), ReductionStatus::Complete);
        assert_eq!(result.structure, source.structure);
        assert_eq!(
            result
                .contract(ContractSettings {
                    expand: false,
                    ..Default::default()
                })
                .unwrap(),
            result
        );
    }
}

#[test]
fn minimal_contract_preserves_independent_alternatives_and_full_policy_remains_eligible() {
    let source = tensor(
        "(p(spenso::mink(4,a))+q(spenso::mink(4,a)))*(r(spenso::mink(4,a))+s(spenso::mink(4,a)))",
    );
    let expected = input(
        "spenso::g(p(spenso::mink(4)),r(spenso::mink(4)))+spenso::g(p(spenso::mink(4)),s(spenso::mink(4)))+spenso::g(q(spenso::mink(4)),r(spenso::mink(4)))+spenso::g(q(spenso::mink(4)),s(spenso::mink(4)))",
    );
    let result = source
        .contract(ContractSettings {
            expand: false,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(result.expression, source.expression);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert!(!result.contraction_complete());
    assert_eq!(
        result.contract(Default::default()).unwrap().expression,
        expected
    );
    let algebra = source
        .simplify_algebra(&AlgebraSettings {
            contract: AlgebraContraction::Minimal,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(algebra.expression, source.expression);
    assert_eq!(algebra.reduction_status(), ReductionStatus::Complete);
    assert_eq!(
        algebra.contract(Default::default()).unwrap().expression,
        expected
    );
}

#[test]
fn minimal_contract_keeps_tensor_powers_without_multiplying_alternatives() {
    let atomic = tensor("p(spenso::mink(4,a))^2")
        .contract(ContractSettings {
            expand: false,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(
        atomic.expression,
        input("spenso::g(p(spenso::mink(4)),p(spenso::mink(4)))")
    );
    assert_eq!(atomic.reduction_status(), ReductionStatus::Complete);
    let nested = tensor("(p(spenso::mink(4,a))+spenso::g(spenso::mink(4,a),spenso::mink(4,b))*q(spenso::mink(4,b)))^2")
        .contract(ContractSettings { expand: false, ..Default::default() }).unwrap();
    assert_eq!(
        nested.expression,
        input("(p(spenso::mink(4,a))+q(spenso::mink(4,a)))^2")
    );
    assert_eq!(nested.reduction_status(), ReductionStatus::Complete);
    let readmitted = SymbolicTensor::infer(nested.expression.clone()).unwrap();
    assert_eq!(readmitted.structure, nested.structure);
    let source = tensor("(p(spenso::mink(4,a))+q(spenso::mink(4,a)))^2");
    let result = source
        .contract(ContractSettings {
            expand: false,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(result.expression, source.expression);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    let algebra = source
        .simplify_algebra(&AlgebraSettings {
            contract: AlgebraContraction::Minimal,
            max_passes: Some(0),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(algebra.expression, source.expression);
    assert_eq!(algebra.reduction_status(), ReductionStatus::Complete);
    assert_eq!(
        result.contract(Default::default()).unwrap().expression,
        input(
            "spenso::g(p(spenso::mink(4)),p(spenso::mink(4)))+2*spenso::g(p(spenso::mink(4)),q(spenso::mink(4)))+spenso::g(q(spenso::mink(4)),q(spenso::mink(4)))"
        )
    );
}

#[test]
fn minimal_contract_zero_budget_completes_without_permitted_work() {
    let source = tensor(
        "(p(spenso::mink(4,a))+q(spenso::mink(4,a)))*(r(spenso::mink(4,a))+s(spenso::mink(4,a)))",
    );
    let result = source
        .contract(ContractSettings {
            expand: false,
            max_passes: Some(0),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(result.expression, source.expression);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    let active = tensor("p(spenso::mink(4,a))*q(spenso::mink(4,a))");
    assert_eq!(
        active
            .contract(ContractSettings {
                expand: false,
                max_passes: Some(0),
                ..Default::default()
            })
            .unwrap()
            .reduction_status(),
        ReductionStatus::Capped
    );
    let disabled = tensor("p(spenso::mink(4,a))^2")
        .contract(ContractSettings {
            expand: false,
            rank_one: false,
            max_passes: Some(0),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(disabled.reduction_status(), ReductionStatus::Complete);
}

#[test]
fn minimal_contract_keeps_identity_prerequisites_and_scalar_spectators() {
    use crate::{dirac::GammaSimplifySettings, gamma};
    use spenso::g;
    let reps = crate::test_support::test_initialize();
    let [left, middle, right] =
        [99501, 99502, 99503].map(|i| reps.bis4.to_symbolic([Atom::num(i)]));
    let [mu, nu] = [99504, 99505].map(|i| reps.mink4.to_symbolic([Atom::num(i)]));
    let spectator = input(
        "(p(spenso::mink(4,a))+q(spenso::mink(4,a)))*(r(spenso::mink(4,a))+s(spenso::mink(4,a)))",
    );
    let source = SymbolicTensor::infer(
        &spectator * g!(&mu, &nu) * gamma!(&left, &middle, &mu) * gamma!(&middle, &right, &nu),
    )
    .unwrap();
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        contract: AlgebraContraction::Minimal,
        ..Default::default()
    };
    let result = source.simplify_algebra(&settings).unwrap();
    assert_eq!(
        result.expression,
        Atom::num(4) * g!(&left, &right) * spectator
    );
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.simplify_algebra(&settings).unwrap(), result);
}

#[test]
fn minimal_contract_respects_connection_filters_and_trace_collection() {
    use spenso::structure::{
        abstract_index::AbstractIndex,
        representation::{LibraryRep, RepName},
        slot::IsAbstractSlot,
    };
    let reps = crate::test_support::test_initialize();
    let lorentz = LibraryRep::from(reps.mink4.rep);
    let [i, j] = [99511, 99512].map(|i| reps.bis4.to_symbolic([Atom::num(i)]));
    let head = SPENSO_TAG.tensor_symbol("minimal_contract::mixed");
    let matrix = |slot: &Atom| {
        FunctionBuilder::new(head)
            .add_arg(&i)
            .add_arg(&j)
            .add_arg(slot)
            .finish()
    };
    let a = input("spenso::mink(4,a)");
    let b = input("spenso::mink(4,b)");
    let source = SymbolicTensor::infer(spenso::g!(&a, &b) * matrix(&a)).unwrap();
    let empty = source
        .contract(ContractSettings {
            expand: false,
            representations: Some(&[]),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(empty.expression, source.expression);
    assert_eq!(empty.reduction_status(), ReductionStatus::Complete);
    let filtered = source
        .contract(ContractSettings {
            expand: false,
            representations: Some(&[lorentz]),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(filtered.expression, matrix(&b));
    assert_eq!(filtered.structure, source.structure);

    let rep = LibraryRep::new_dual("minimal_contract::R")
        .unwrap()
        .new_rep(3);
    let [i, j] = [99521, 99522].map(|i| rep.slot::<AbstractIndex, _>(i).to_atom());
    let [di, dj] = [99521, 99522].map(|i| rep.dual().slot::<AbstractIndex, _>(i).to_atom());
    let [a, b] = ["minimal_contract::matrixA", "minimal_contract::matrixB"]
        .map(|name| SPENSO_TAG.tensor_symbol(name));
    let source =
        SymbolicTensor::infer(symbolica::function!(a, &i, &dj) * symbolica::function!(b, &j, &di))
            .unwrap();
    let untraced = source
        .contract(ContractSettings {
            expand: false,
            collect_traces: false,
            ..Default::default()
        })
        .unwrap();
    assert!(!untraced.expression.contains_symbol(SPENSO_TAG.trace));
    let untouched = source
        .contract(ContractSettings {
            expand: false,
            collect_traces: false,
            max_passes: Some(0),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(untouched.reduction_status(), ReductionStatus::Complete);
    let traced = untraced
        .contract(ContractSettings {
            expand: false,
            ..Default::default()
        })
        .unwrap();
    assert!(traced.expression.contains_symbol(SPENSO_TAG.trace));
    assert!(!matches!(traced.expression.as_view(), AtomView::Num(_)));
    assert_eq!(traced.reduction_status(), ReductionStatus::Complete);
}

#[test]
fn minimal_contract_checks_callback_interface_changes() {
    crate::test_support::test_initialize();
    let a = input("spenso::mink(4,a)");
    let b = input("spenso::mink(4,b)");
    let target = b.clone();
    let head = spenso::tensor_symbol!(
        "minimal_contract::callback_rank_loss",
        norm = move |value, output| {
            if let AtomView::Fun(function) = value
                && function.iter().next() == Some(target.as_view())
            {
                **output = Atom::num(7);
            }
        }
    );
    let source =
        SymbolicTensor::infer(spenso::g!(&a, &b) * symbolica::function!(head, &a)).unwrap();
    assert!(
        source
            .contract(ContractSettings {
                expand: false,
                ..Default::default()
            })
            .is_err()
    );
}

#[test]
fn minimal_contract_leaves_unrelated_scalar_arithmetic_opaque() {
    let source = tensor("p(spenso::mink(4,a))*(q(spenso::mink(4,a))+r(spenso::mink(4,a)))");
    let settings = ContractSettings {
        expand: false,
        ..Default::default()
    };
    crate::tensor::SELECTED_OPENS.with(|count| count.set(0));
    let bare = source.contract(settings).unwrap();
    let baseline = crate::tensor::SELECTED_OPENS.with(|count| count.get());
    let spectator = input(&format!(
        "({})^7",
        (0..256)
            .map(|i| format!("x{i}"))
            .collect::<Vec<_>>()
            .join("+")
    ));
    let source = SymbolicTensor::infer(&spectator * source.expression).unwrap();
    crate::tensor::SELECTED_OPENS.with(|count| count.set(0));
    let result = source.contract(settings).unwrap();
    assert_eq!(result.expression, spectator * bare.expression);
    assert_eq!(
        crate::tensor::SELECTED_OPENS.with(|count| count.get()),
        baseline
    );
}

#[test]
fn metrics_between_external_ports_complete_beside_internal_dummies() {
    use crate::color::ColorSimplifySettings;
    // The internal f-f dummy does not make a metric on two open ports a
    // contraction source; the control has no internal dummy at all.
    for source in [
        "spenso::f(spenso::coad(8,a1),spenso::coad(8,a4),spenso::coad(8,d))*spenso::f(spenso::coad(8,a2),spenso::coad(8,a5),spenso::coad(8,d))*spenso::g(spenso::coad(8,b0),spenso::coad(8,b1))",
        "spenso::g(spenso::coad(8,b0),spenso::coad(8,b1))*spenso::f(spenso::coad(8,a1),spenso::coad(8,a2),spenso::coad(8,a3))",
    ] {
        let source = tensor(source);
        for settings in [
            ContractSettings::default(),
            ContractSettings {
                collect_chains: false,
                collect_traces: false,
                ..Default::default()
            },
        ] {
            let contracted = source.contract(settings).unwrap();
            assert_eq!(
                contracted.reduction_status(),
                ReductionStatus::Complete,
                "{source:?}"
            );
            assert_eq!(contracted.expression, source.expression);
        }
        for (color, contract) in [
            (ColorSimplifySettings::default(), AlgebraContraction::None),
            (
                ColorSimplifySettings::default().with_cof_dimension_invariants(),
                AlgebraContraction::Fully,
            ),
        ] {
            let reduced = source
                .simplify_algebra(&AlgebraSettings {
                    color: Some(color),
                    contract,
                    ..Default::default()
                })
                .unwrap();
            assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
            assert_eq!(reduced.structure, source.structure);
        }
    }
}
