use super::*;
use crate::{gamma, gamma5};
use spenso::structure::partial::PartialStructureExt;
use spenso::{g, mink, p, q, trace};
use symbolica::{
    atom::{AtomCore, FunctionBuilder},
    symbol,
};

#[test]
fn free_trace_definitions_discharge_cleanup_without_entering_contraction() {
    use crate::tensor::contract::CONTRACT_DOMAIN_CALLS;
    let reps = crate::test_support::test_initialize();
    let slots = (0..8)
        .map(|index| mink!(4, { Atom::num(98600 + index) }))
        .collect::<Vec<_>>();
    let source = SymbolicTensor::infer(
        trace!(reps.bis4.to_symbolic([]); slots.iter().map(|slot| gamma!(slot))),
    )
    .unwrap();
    CONTRACT_DOMAIN_CALLS.with(|count| count.set(0));
    let reduced = source
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            // The trace identity uses this budget. Its free-port definitions
            // require no additional contraction, even at their own boundaries.
            max_passes: 1,
            ..Default::default()
        })
        .unwrap();
    assert_eq!(reduced.reduction_status(), ReductionStatus::Complete);
    assert!(!reduced.aliases().unwrap().is_empty());
    assert_eq!(reduced.root().structure, source.structure);
    assert_eq!(
        CONTRACT_DOMAIN_CALLS.with(|count| count.get()),
        0,
        "free-port alias definitions must use their observations before contraction setup"
    );
}

#[test]
fn unchanged_shared_definition_is_observed_once_and_completed_reruns_reuse_it() {
    use observation::{INITIAL_DOMAINS, INITIAL_SCANS};
    crate::test_support::test_initialize();
    let body = SymbolicTensor::infer(p!(mink!(4, 98701)) * q!(mink!(4, 98701))).unwrap();
    let handle = body.alias_handle().unwrap();
    let root = SymbolicTensor::infer(handle.expression.pow(3)).unwrap();
    INITIAL_SCANS.with(|count| count.set(0));
    INITIAL_DOMAINS.with(|domains| domains.borrow_mut().clear());
    let value = Arc::new(root.with_aliases([(handle, body)]).unwrap());
    assert_eq!(INITIAL_SCANS.with(|count| count.get()), 2);
    let reduced = value.contract(Default::default()).unwrap();
    let domains = 1 + reduced.aliases().unwrap().len();
    assert_eq!(
        INITIAL_SCANS.with(|count| count.get()),
        domains,
        "each original or newly introduced domain is observed once: {:?}",
        (
            INITIAL_DOMAINS.with(|domains| domains.borrow().clone()),
            reduced.aliases().unwrap()
        )
    );
    assert_eq!(
        reduced.resolved().unwrap().expression,
        g!(p!(mink!(4)), q!(mink!(4))).pow(3)
    );
    let repeated = reduced.contract(Default::default()).unwrap();
    assert!(Arc::ptr_eq(&reduced, &repeated));
    assert_eq!(INITIAL_SCANS.with(|count| count.get()), domains);
}

#[test]
fn regional_inventory_preserves_a_surviving_head_and_reuses_scalar_spectators() {
    use observation::{DomainObservations, REUSED_REGIONS};
    crate::test_support::test_initialize();
    let epsilon = *crate::epsilon::EPSILON_SYMBOL;
    let left = FunctionBuilder::new(epsilon)
        .add_arg(Atom::var(symbol!("inventory_left")))
        .finish();
    let right = FunctionBuilder::new(epsilon)
        .add_arg(Atom::var(symbol!("inventory_right")))
        .finish();
    let spectator =
        (Atom::var(symbol!("inventory_x")) + Atom::var(symbol!("inventory_y"))).pow(100);
    let original = &left * &right * &spectator;
    let tensor =
        SymbolicTensor::from_normalized_parts(original, PartialStructure::from_logical_slots([]));
    let observed = tensor.reduction_observations();
    assert_eq!(observed.candidates.counts[4], 2);
    assert_eq!(observed.epsilon_degree, 2);
    REUSED_REGIONS.with(|count| count.set(0));
    let next: DomainObservations = observed.updated((&right * &spectator).as_view());
    assert_eq!(next.candidates.counts[4], 1);
    assert_eq!(next.epsilon_degree, 1);
    assert!(next.candidates.symbols[4]);
    assert_eq!(REUSED_REGIONS.with(|count| count.get()), 2);
    let sum = observed.updated((&left + &right).as_view());
    assert_eq!(sum.epsilon_degree, 1);
    let squared = sum.updated((&left + &right).pow(2).as_view());
    assert_eq!(squared.epsilon_degree, 2);
}

#[test]
fn gamma_prerequisites_do_not_enable_epsilon_identities() {
    crate::test_support::test_initialize();
    let spin = crate::representations::Bispinor {};
    use spenso::structure::representation::RepName;
    let slots = (0..4)
        .map(|i| mink!(4, { Atom::num(98800 + i) }))
        .collect::<Vec<_>>();
    let trace = trace!(spin.to_symbolic([4]); std::iter::once(gamma5!()).chain(slots.iter().map(|slot| gamma!(slot))));
    let source = SymbolicTensor::infer(trace).unwrap();
    let gamma_only = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        ..Default::default()
    };
    let result = source.simplify_algebra(&gamma_only).unwrap();
    assert!(
        result
            .resolved()
            .unwrap()
            .expression
            .contains_symbol(*crate::epsilon::EPSILON_SYMBOL)
    );
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert!(Arc::ptr_eq(
        &result,
        &result.simplify_algebra(&gamma_only).unwrap()
    ));
}

#[test]
fn a_capped_algebra_result_retains_its_unprocessed_trace_and_aliases() {
    let reps = crate::test_support::test_initialize();
    let expression =
        trace!(reps.bis4.to_symbolic([]); [gamma!(mink!(4, 98901)), gamma!(mink!(4, 98902))]);
    let body = SymbolicTensor::infer(expression).unwrap();
    let handle = body.alias_handle().unwrap();
    let source = Arc::new(handle.clone().with_aliases([(handle, body)]).unwrap());
    let settings = AlgebraSettings {
        gamma: Some(GammaSimplifySettings::default()),
        max_passes: 0,
        ..Default::default()
    };
    let capped = source.simplify_algebra(&settings).unwrap();
    assert_eq!(capped.reduction_status(), ReductionStatus::Capped);
    assert_eq!(source.root(), capped.root());
    assert_eq!(source.aliases().unwrap(), capped.aliases().unwrap());
    let finished = capped
        .simplify_algebra(&AlgebraSettings {
            max_passes: 16,
            ..settings
        })
        .unwrap();
    assert_eq!(
        finished.expanded().unwrap().expression,
        Atom::num(4) * g!(mink!(4, 98901), mink!(4, 98902))
    );
}

#[test]
fn gamma_identity_contracts_its_prerequisite_without_opening_scalar_spectators() {
    let reps = crate::test_support::test_initialize();
    let left = reps.bis4.to_symbolic([Atom::num(99001)]);
    let middle = reps.bis4.to_symbolic([Atom::num(99002)]);
    let right = reps.bis4.to_symbolic([Atom::num(99003)]);
    let mu = mink!(4, 99004);
    let nu = mink!(4, 99005);
    let spectator = (Atom::one() + p!(mink!(4, 99006)) * q!(mink!(4, 99006))).pow(20);
    let source = SymbolicTensor::infer(
        &spectator * g!(&mu, &nu) * gamma!(&left, &middle, &mu) * gamma!(&middle, &right, &nu),
    )
    .unwrap();
    let result = source
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            ..Default::default()
        })
        .unwrap();
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert_eq!(result.root().structure, source.structure);
    // Clifford identity gamma^mu gamma_mu = 4 I; the independent scalar
    // coefficient retains even its internal contractions and power notation.
    assert_eq!(
        result.resolved().unwrap().expression,
        Atom::num(4) * g!(&left, &right) * spectator
    );
}

#[test]
fn zero_budget_is_complete_when_no_permitted_work_exists() {
    crate::test_support::test_initialize();
    let constant = SymbolicTensor::infer(Atom::num(7)).unwrap();
    let settings = ContractSettings {
        max_passes: 0,
        ..Default::default()
    };
    assert_eq!(
        constant.contract(settings).unwrap().reduction_status(),
        ReductionStatus::Complete
    );
    let indexed = SymbolicTensor::infer(p!(mink!(4, 99101)) * q!(mink!(4, 99101))).unwrap();
    let result = indexed
        .contract(ContractSettings {
            representations: Some(&[]),
            ..settings
        })
        .unwrap();
    assert_eq!(result.root(), indexed);
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
}

#[test]
fn dummy_reservation_reuses_observed_domains_including_opaque_metadata() {
    use super::super::composition::DUMMY_RESERVATION_WALKS;
    use spenso::structure::abstract_index::AbstractIndex;
    crate::test_support::test_initialize();
    let mu = mink!(4, 99201);
    let hidden = mink!(4, 99202);
    let metadata = FunctionBuilder::new(symbol!("reservation_metadata"; Scalar))
        .add_arg(&hidden)
        .finish();
    let source = SymbolicTensor::infer(p!(&mu) * q!(&mu) * &metadata).unwrap();
    let observed = source.reduction_observations();
    assert_eq!(
        observed.reserved_indices,
        HashSet::from([AbstractIndex::Normal(99201), AbstractIndex::Normal(99202),])
    );
    DUMMY_RESERVATION_WALKS.with(|count| count.set(0));
    let _ = SymbolicTensor::reserved_dummies([&source, &source]);
    assert_eq!(DUMMY_RESERVATION_WALKS.with(|count| count.get()), 0);
    let replaced = source.with_identity_result(metadata).unwrap();
    assert_eq!(
        replaced.reduction_observations().reserved_indices,
        HashSet::from([AbstractIndex::Normal(99202)])
    );
}

#[test]
fn unchanged_kernel_regions_are_reused_after_an_independent_family_changes() {
    use observation::REUSED_KERNEL_REGIONS;
    crate::test_support::test_initialize();
    let left = crate::bis!(4, 99301);
    let right = crate::bis!(4, 99302);
    let matrix = gamma!(&left, &right, mink!(4, 99303));
    let epsilon = crate::epsilon!(mink!(2, 99304), mink!(2, 99305));
    let source = SymbolicTensor::infer(&matrix * epsilon.pow(2)).unwrap();
    REUSED_KERNEL_REGIONS.with(|count| count.set(0));
    let result = source
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            epsilon: true,
            ..Default::default()
        })
        .unwrap();
    // The two-dimensional determinant contracts to 2; its independent scalar
    // reduction does not change the open Clifford factor's inputs or boundary.
    let explicit = result.undo_chain().unwrap().expanded().unwrap().expression;
    assert!(
        (&explicit - Atom::num(2) * matrix).expand().is_zero(),
        "the explicit result is {explicit}"
    );
    assert_eq!(result.reduction_status(), ReductionStatus::Complete);
    assert!(REUSED_KERNEL_REGIONS.with(|count| count.get()) > 0);
}

#[test]
fn trusted_replacements_reuse_unchanged_opaque_boundary_interfaces() {
    crate::test_support::test_initialize();
    let mu = mink!(4, 99401);
    let nu = mink!(4, 99402);
    let spectator = (p!(mink!(4, 99403)) * q!(mink!(4, 99403))
        + p!(mink!(4, 99404)) * q!(mink!(4, 99404)) * Atom::num(3))
    .pow(12);
    let source = SymbolicTensor::infer(g!(&mu, &nu) * p!(&nu) * &spectator).unwrap();
    source.reduction_observations();
    source.shallow_graph().unwrap();
    let interfaces = source.proofs.leaf_interfaces.get().unwrap();
    let prior = interfaces.lock().unwrap().clone();
    assert!(prior.contains_key(&spectator));

    // This is the metric substitution's exact result. Its spectator keeps its
    // admitted scalar interface even though its interior still has dummy pairs.
    let rewritten = source.with_identity_result(p!(&mu) * &spectator).unwrap();
    assert_ne!(rewritten.expression, source.expression);
    assert!(Arc::ptr_eq(
        interfaces,
        rewritten.proofs.leaf_interfaces.get().unwrap()
    ));
    rewritten.shallow_graph().unwrap();
    let current = interfaces.lock().unwrap();
    assert_eq!(current.get(&spectator), prior.get(&spectator));
    assert_eq!(current.len(), prior.len() + 1);
    assert!(current.contains_key(&p!(&mu)));
}
