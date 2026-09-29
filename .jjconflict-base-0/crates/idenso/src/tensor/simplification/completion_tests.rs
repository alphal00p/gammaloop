use super::*;
use crate::gamma;
use spenso::{g, mink, p, q, trace};
use symbolica::atom::AtomCore;

#[test]
fn alternating_completed_policies_reuse_the_shared_carrier_and_contraction_proof() {
    crate::test_support::test_initialize();
    let body = SymbolicTensor::infer(p!(mink!(4, 99501)) * q!(mink!(4, 99501))).unwrap();
    let handle = body.alias_handle().unwrap();
    let root = SymbolicTensor::infer(handle.expression.pow(3)).unwrap();
    let value = Arc::new(root.with_aliases([(handle, body)]).unwrap());
    let contracted = value.contract(Default::default()).unwrap();
    let settings = AlgebraSettings::hep();
    let completed = contracted.reduce_algebra(&settings).unwrap();
    assert_eq!(completed.root(), contracted.root());
    assert_eq!(completed.aliases().unwrap(), contracted.aliases().unwrap());
    assert!(completed.contraction_complete());

    for _ in 0..3 {
        let structural = completed.contract(Default::default()).unwrap();
        assert!(Arc::ptr_eq(&completed, &structural));
        let algebra = structural.reduce_algebra(&settings).unwrap();
        assert!(Arc::ptr_eq(&completed, &algebra));
        assert_eq!(algebra.reduction_status(), ReductionStatus::Complete);
    }

    // An explicit callback boundary invalidates facts even if the callback
    // happens to preserve every expression in this invocation.
    let mut calls = 0;
    let mapped = completed
        .map_aliases(|_, body| {
            calls += 1;
            Ok(body)
        })
        .unwrap();
    assert!(calls > 0);
    assert!(mapped.proofs.reduction.is_none());
    assert!(!mapped.contraction_complete());
}

#[test]
fn root_or_definition_changes_invalidate_previously_completed_policies() {
    let reps = crate::test_support::test_initialize();
    let mu = mink!(4, 99511);
    let nu = mink!(4, 99512);
    let expression = trace!(reps.bis4.to_symbolic([]); [gamma!(&mu), gamma!(&nu)]);
    let body = SymbolicTensor::infer(expression).unwrap();
    let handle = body.alias_handle().unwrap();
    let values = [
        Arc::new(body.clone().with_aliases([]).unwrap()),
        Arc::new(handle.clone().with_aliases([(handle, body)]).unwrap()),
    ];
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        ..Default::default()
    };
    let structural_policy = ReductionRequest::Contract(Default::default()).policy();
    for value in values {
        let contracted = value.contract(Default::default()).unwrap();
        assert!(contracted.contraction_complete());
        let reduced = contracted.reduce_algebra(&settings).unwrap();
        assert_eq!(
            reduced.expanded().unwrap().expression,
            Atom::num(4) * g!(&mu, &nu)
        );
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
        assert!(!Arc::ptr_eq(&reduced, &completed));
        assert!(completed.contraction_complete());
    }
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
        max_passes: 0,
        ..Default::default()
    };
    let disabled = AlgebraSettings::default();
    let completed = body.reduce_algebra(&disabled).unwrap();
    let capped = completed.reduce_algebra(&settings).unwrap();
    let repeated = capped.reduce_algebra(&settings).unwrap();
    assert_eq!(capped.reduction_status(), ReductionStatus::Capped);
    assert_eq!(repeated.reduction_status(), ReductionStatus::Capped);
    assert!(!Arc::ptr_eq(&capped, &repeated));
    assert_eq!(repeated.root(), completed.root());
    assert_eq!(repeated.aliases().unwrap(), completed.aliases().unwrap());
    assert_eq!(
        repeated.proofs.reduction.as_ref().unwrap().completed.len(),
        1
    );
    let restored = repeated.reduce_algebra(&disabled).unwrap();
    assert_eq!(restored.reduction_status(), ReductionStatus::Complete);
    assert!(Arc::ptr_eq(
        &restored,
        &restored.reduce_algebra(&disabled).unwrap()
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
    let disabled = AlgebraSettings::default();
    let settings = ContractSettings {
        collect_chains: false,
        collect_traces: false,
        ..Default::default()
    };
    let completed = source.reduce_algebra(&disabled).unwrap();
    let deferred = completed.contract(settings).unwrap();
    let repeated = deferred.contract(settings).unwrap();
    assert_eq!(deferred.reduction_status(), ReductionStatus::Deferred);
    assert_eq!(repeated.reduction_status(), ReductionStatus::Deferred);
    assert!(!Arc::ptr_eq(&deferred, &repeated));
    assert_eq!(repeated.resolved().unwrap().expression, expression);
    assert_eq!(
        repeated.proofs.reduction.as_ref().unwrap().completed.len(),
        1
    );
    let restored = repeated.reduce_algebra(&disabled).unwrap();
    assert_eq!(restored.reduction_status(), ReductionStatus::Complete);
    assert!(Arc::ptr_eq(
        &restored,
        &restored.reduce_algebra(&disabled).unwrap()
    ));
}
