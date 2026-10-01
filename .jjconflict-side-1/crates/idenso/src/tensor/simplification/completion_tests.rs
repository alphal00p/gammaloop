use super::*;
use crate::gamma;
use spenso::{g, mink, p, q, trace};
use symbolica::atom::{Atom, AtomCore};

#[test]
fn alternating_completed_policies_reuse_the_shared_carrier_and_contraction_proof() {
    crate::test_support::test_initialize();
    let value =
        SymbolicTensor::infer((Atom::one() + p!(mink!(4, 99501)) * q!(mink!(4, 99501))).pow(3))
            .unwrap();
    let contracted = value.contract(Default::default()).unwrap();
    let settings = AlgebraSettings::hep();
    let completed = contracted.simplify_algebra(&settings).unwrap();
    assert_eq!(completed, contracted);
    assert!(completed.contraction_complete());

    for _ in 0..3 {
        let structural = completed.contract(Default::default()).unwrap();
        assert!(Arc::ptr_eq(
            completed.reduction_observations(),
            structural.reduction_observations()
        ));
        let algebra = structural.simplify_algebra(&settings).unwrap();
        assert!(Arc::ptr_eq(
            completed.reduction_observations(),
            algebra.reduction_observations()
        ));
        assert_eq!(algebra.reduction_status(), ReductionStatus::Complete);
    }

    // A newly admitted expression does not inherit completion facts.
    let mapped = SymbolicTensor::infer(completed.expression.clone()).unwrap();
    assert!(mapped.proofs.reduction.is_none());
    assert!(!mapped.contraction_complete());
}

#[test]
fn expression_changes_invalidate_previously_completed_policies() {
    let reps = crate::test_support::test_initialize();
    let mu = mink!(4, 99511);
    let nu = mink!(4, 99512);
    let expression = trace!(reps.bis4.to_symbolic([]); [gamma!(&mu), gamma!(&nu)]);
    let value = SymbolicTensor::infer(expression).unwrap();
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let structural_policy = ReductionRequest::Contract(Default::default()).policy();
    let contracted = value.contract(Default::default()).unwrap();
    assert!(contracted.contraction_complete());
    let reduced = contracted.simplify_algebra(&settings).unwrap();
    assert_eq!(reduced.expression.expand(), Atom::num(4) * g!(&mu, &nu));
    assert!(!reduced.contraction_complete());
    assert!(
        !reduced
            .proofs
            .reduction
            .as_ref()
            .unwrap()
            .completed
            .contains(&structural_policy)
    );
    let completed = reduced.contract(Default::default()).unwrap();
    assert_eq!(reduced, completed);
    assert!(completed.contraction_complete());
}

#[test]
fn capped_policy_is_revisited_without_losing_other_completion_facts() {
    let reps = crate::test_support::test_initialize();
    let body = SymbolicTensor::infer(
        trace!(reps.bis4.to_symbolic([]); [gamma!(mink!(4, 99521)), gamma!(mink!(4, 99522))]),
    )
    .unwrap();
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        max_passes: Some(0),
        ..Default::default()
    };
    let disabled = AlgebraSettings {
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let completed = body.simplify_algebra(&disabled).unwrap();
    let capped = completed.simplify_algebra(&settings).unwrap();
    let repeated = capped.simplify_algebra(&settings).unwrap();
    assert_eq!(capped.reduction_status(), ReductionStatus::Capped);
    assert_eq!(repeated.reduction_status(), ReductionStatus::Capped);
    assert_eq!(repeated, completed);
    assert_eq!(
        repeated.proofs.reduction.as_ref().unwrap().completed.len(),
        1
    );
    let restored = repeated.simplify_algebra(&disabled).unwrap();
    assert_eq!(restored.reduction_status(), ReductionStatus::Complete);
    assert!(Arc::ptr_eq(
        restored.reduction_observations(),
        restored
            .simplify_algebra(&disabled)
            .unwrap()
            .reduction_observations()
    ));
}

#[test]
fn deferred_policy_is_revisited_without_losing_other_completion_facts() {
    crate::test_support::test_initialize();
    let port = mink!(4, 99531);
    let r = spenso::vector_symbol!("completion_cache_inexact::r");
    let s = spenso::vector_symbol!("completion_cache_inexact::s");
    let rounded = Atom::num(symbolica::domains::float::Float::parse("0.1", Some(11)).unwrap());
    let expression = (rounded * p!(&port) + q!(&port))
        * (symbolica::function!(r, &port) + symbolica::function!(s, &port));
    let source = SymbolicTensor::infer(expression.clone()).unwrap();
    let disabled = AlgebraSettings {
        contract: AlgebraContraction::None,
        ..Default::default()
    };
    let settings = ContractSettings {
        collect_chains: false,
        collect_traces: false,
        ..Default::default()
    };
    let completed = source.simplify_algebra(&disabled).unwrap();
    let deferred = completed.contract(settings).unwrap();
    let repeated = deferred.contract(settings).unwrap();
    assert_eq!(deferred.reduction_status(), ReductionStatus::Deferred);
    assert_eq!(repeated.reduction_status(), ReductionStatus::Deferred);
    assert_eq!(repeated.expression, expression);
    assert_eq!(
        repeated.proofs.reduction.as_ref().unwrap().completed.len(),
        1
    );
    let restored = repeated.simplify_algebra(&disabled).unwrap();
    assert_eq!(restored.reduction_status(), ReductionStatus::Complete);
    assert!(Arc::ptr_eq(
        restored.reduction_observations(),
        restored
            .simplify_algebra(&disabled)
            .unwrap()
            .reduction_observations()
    ));
}
