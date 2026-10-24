//! Independent typed comparisons against shipped lower-loop programs.
//! No five-loop package, rule generation, or numerical terminal values enter.
use super::*;
use rustred::algebra::Coefficient;
use rustred::family::AffineDenominator;
use rustred::reduction::terminal_normalization::{
    VacuumFactor, VacuumProductReducer, VacuumProductReduction,
};
use std::collections::HashSet;

type QualifiedProduct = Vec<(String, Vec<i64>)>;

fn key(powers: impl IntoIterator<Item = i64>) -> IntegralKey {
    IntegralKey::try_new(powers).unwrap()
}

fn assert_native_result(
    actual: &VacuumFactorReduction,
    expected_family: &IntegralFamily,
    expected_terms: &BTreeMap<IntegralKey, Coefficient>,
    expected_conditions: &[Coefficient],
    expected_work: ReductionStatistics,
) {
    let terms = actual
        .terms
        .iter()
        .map(|term| {
            assert_eq!(term.family.fingerprint(), expected_family.fingerprint());
            assert!(
                term.family
                    .coefficient_context()
                    .has_same_variable_map(expected_family.coefficient_context())
            );
            assert_eq!(
                term.integral.powers().len(),
                expected_family.denominator_count()
            );
            (term.integral.clone(), term.coefficient.clone())
        })
        .collect::<BTreeMap<_, _>>();
    assert_eq!(
        terms.len(),
        actual.terms.len(),
        "provider outputs must not repeat a terminal"
    );
    assert_eq!(&terms, expected_terms);
    assert_eq!(actual.nonzero_conditions, expected_conditions);
    assert_eq!(
        actual.work,
        VacuumFactorWork {
            rule_applications: expected_work.rule_applications(),
            coalescing_additions: expected_work.coalescing_additions(),
            cache_hits: expected_work.cache_hits(),
        }
    );
}

fn renamed_family(
    source: &IntegralFamily,
    name: &str,
    denominators: Vec<AffineDenominator>,
) -> Arc<IntegralFamily> {
    Arc::new(
        IntegralFamily::new_with_limits(
            name.to_owned(),
            source.loop_momenta().to_vec(),
            source.external_momenta().to_vec(),
            source.coefficient_context().clone(),
            source.dimension().clone(),
            denominators,
            source.external_gram().to_vec(),
            source.power_shifts().to_vec(),
            IntegralFamilyLimits::default(),
        )
        .unwrap(),
    )
}

#[test]
fn three_loop_provider_matches_typed_corner_and_dotted_reductions() {
    crate::Vakint::initialize_vakint_symbols();
    let artifact = shipped_family_for_loop_count(3)
        .unwrap()
        .artifact()
        .unwrap();
    let family = share_family(artifact.family()).unwrap();
    let renamed = renamed_family(
        &family,
        "provider-test-name-is-not-a-selector",
        family.denominators().to_vec(),
    );
    assert_ne!(renamed.fingerprint(), family.fingerprint());
    let corner = key([1; 6]);
    let dotted = key([2, 1, 1, 1, 1, 1]);
    assert!(artifact.masters().contains(&corner));
    assert!(!artifact.masters().contains(&dotted));
    let mut provider = ShippedFactorProvider::default();
    for target in [&corner, &dotted] {
        let direct =
            Reducer::reduce_artifact_unit_mass_guarded(artifact, target, Default::default())
                .unwrap();
        if target == &dotted {
            assert!(direct.2.rule_applications() > 0);
        }
        for source in [&family, &renamed] {
            let actual = provider
                .reduce(Arc::clone(source), target, Default::default())
                .unwrap();
            assert_native_result(&actual, &family, direct.0.terms(), &direct.1, direct.2);
            assert!(
                actual
                    .terms
                    .iter()
                    .all(|term| artifact.masters().contains(&term.integral))
            );
            let warm = provider
                .reduce(Arc::clone(source), target, Default::default())
                .unwrap();
            assert_native_result(&warm, &family, direct.0.terms(), &direct.1, direct.2);
            assert!(
                Arc::ptr_eq(&actual.terms[0].family, &warm.terms[0].family),
                "only immutable family preparation should survive calls"
            );
        }
    }
}

#[test]
fn three_loop_provider_refuses_missing_mass_map_and_recovers_after_limit_error() {
    crate::Vakint::initialize_vakint_symbols();
    let artifact = shipped_family_for_loop_count(3)
        .unwrap()
        .artifact()
        .unwrap();
    let family = share_family(artifact.family()).unwrap();
    let mut different_mass = family.denominators().to_vec();
    different_mass[0] = AffineDenominator::new(
        family.coefficient_context().integer(-2),
        different_mass[0].coefficients().to_vec(),
    );
    let foreign = renamed_family(&family, "unavailable-positive-mass-pattern", different_mass);
    let mut provider = ShippedFactorProvider::default();
    let missing = provider
        .reduce(foreign, &key([1; 6]), Default::default())
        .unwrap_err();
    assert!(
        matches!(&missing, FactorProviderError::Invalid(detail) if detail.contains("no exact momentum map")),
        "{missing}"
    );
    assert!(
        provider
            .reduce(family.clone(), &key([1]), Default::default())
            .is_err()
    );
    let target = key([2, 1, 1, 1, 1, 1]);
    assert!(
        provider
            .reduce(
                family.clone(),
                &target,
                ReductionLimits {
                    max_rule_applications: 0,
                    ..Default::default()
                }
                .with_diagnostics()
            )
            .is_err()
    );
    let direct =
        Reducer::reduce_artifact_unit_mass_guarded(artifact, &target, Default::default()).unwrap();
    let recovered = provider
        .reduce(family, &target, Default::default())
        .unwrap();
    assert_native_result(
        &recovered,
        artifact.family(),
        direct.0.terms(),
        &direct.1,
        direct.2,
    );
}

#[test]
fn signed_lower_factor_matches_the_independent_native_program() {
    crate::Vakint::initialize_vakint_symbols();
    let artifact = shipped_family_for_loop_count(2)
        .unwrap()
        .artifact()
        .unwrap();
    let family = share_family(artifact.family()).unwrap();
    let target = key([1, 1, -1]);
    let direct =
        Reducer::reduce_artifact_unit_mass_guarded(artifact, &target, Default::default()).unwrap();
    let mut provider = ShippedFactorProvider::default();
    let actual = provider
        .reduce(family.clone(), &target, Default::default())
        .unwrap();
    let terms = actual
        .terms
        .iter()
        .map(|term| {
            assert_eq!(term.family.fingerprint(), family.fingerprint());
            (term.integral.clone(), term.coefficient.clone())
        })
        .collect::<BTreeMap<_, _>>();
    assert_eq!(&terms, direct.0.terms());
    let parameter_guards = |conditions: &[Coefficient]| {
        conditions
            .iter()
            .filter_map(|condition| {
                family
                    .coefficient_context()
                    .validate_with_limits(condition, Default::default())
                    .unwrap();
                assert!(!condition.is_zero());
                (!condition.numerator.is_constant() || !condition.denominator.is_constant())
                    .then_some(condition.clone())
            })
            .collect::<HashSet<_>>()
    };
    assert_eq!(
        parameter_guards(&actual.nonzero_conditions),
        parameter_guards(&direct.1)
    );
    assert_eq!(actual.work.rule_applications, direct.2.rule_applications());
}

#[test]
fn four_loop_provider_matches_h_positive_corner_and_dotted_reductions() {
    crate::Vakint::initialize_vakint_symbols();
    // H is deliberately the first installed descriptor. Full physical support
    // exercises the exact mapper without an unrelated failed parent search.
    let engine = crate::rustred_evaluation::cache_tests::external_engine();
    let data = engine.four_loop_data().unwrap();
    let (native, physical) = data.four_loop.as_ref().unwrap().native_program(0).unwrap();
    let family = native.family_shared();
    assert_eq!(physical, 9);
    let corner = key((0..10).map(|axis| i64::from(axis < physical)));
    assert!(
        native.output_terminals().contains(&corner),
        "fixture must be a declared positive terminal"
    );
    let dotted = (0..physical)
        .map(|axis| {
            let mut powers = corner.powers().to_vec();
            powers[axis] += 1;
            key(powers)
        })
        .find(|target| !native.terminals().contains(target))
        .expect("an admitted one-dot fixture must apply a saved formula");
    let mut provider = ShippedFactorProvider::with_data(data);
    for target in [&corner, &dotted] {
        let direct = native
            .reduce_unit_mass_guarded(target, Default::default())
            .unwrap();
        if target == &dotted {
            assert!(direct.2.rule_applications() > 0);
        }
        for _ in 0..2 {
            let actual = provider
                .reduce(family.clone(), target, Default::default())
                .unwrap();
            assert_native_result(&actual, &family, direct.0.terms(), &direct.1, direct.2);
            assert!(
                actual
                    .terms
                    .iter()
                    .all(|term| native.output_terminals().contains(&term.integral))
            );
            assert!(
                actual
                    .terms
                    .iter()
                    .all(|term| Arc::ptr_eq(&term.family, &family))
            );
        }
    }
    assert!(
        provider
            .reduce(
                family.clone(),
                &dotted,
                ReductionLimits {
                    max_rule_applications: 0,
                    ..Default::default()
                }
                .with_diagnostics()
            )
            .is_err()
    );
    let direct = native
        .reduce_unit_mass_guarded(&dotted, Default::default())
        .unwrap();
    let recovered = provider
        .reduce(family, &dotted, Default::default())
        .unwrap();
    assert_native_result(
        &recovered,
        &native.family_shared(),
        direct.0.terms(),
        &direct.1,
        direct.2,
    );
}

fn product_map(result: &VacuumProductReduction) -> BTreeMap<QualifiedProduct, Coefficient> {
    result
        .terms()
        .iter()
        .map(|(factors, coefficient)| {
            let mut product = factors
                .iter()
                .map(|factor| {
                    (
                        factor.family_fingerprint().to_owned(),
                        factor.integral().powers().to_vec(),
                    )
                })
                .collect::<Vec<_>>();
            product.sort();
            (product, coefficient.clone())
        })
        .collect()
}

#[test]
fn five_loop_three_times_two_product_matches_independent_typed_lower_expansion() {
    crate::Vakint::initialize_vakint_symbols();
    let three = shipped_family_for_loop_count(3)
        .unwrap()
        .artifact()
        .unwrap();
    let two = shipped_family_for_loop_count(2)
        .unwrap()
        .artifact()
        .unwrap();
    let one = shipped_family_for_loop_count(1)
        .unwrap()
        .artifact()
        .unwrap();
    let three_family = share_family(three.family()).unwrap();
    let two_family = share_family(two.family()).unwrap();
    let three_target = key([1; 6]);
    let two_target = key([2, 2, 1]);
    let direct_three =
        Reducer::reduce_artifact_unit_mass_guarded(three, &three_target, Default::default())
            .unwrap();
    let direct_two =
        Reducer::reduce_artifact_unit_mass_guarded(two, &two_target, Default::default()).unwrap();
    let direct_one =
        Reducer::reduce_artifact_unit_mass_guarded(one, &key([1]), Default::default()).unwrap();
    let c = three.family().coefficient_context();
    assert_eq!(
        direct_three.0.terms(),
        &BTreeMap::from([(three_target.clone(), c.one())])
    );
    assert_eq!(direct_one.0.terms(), &BTreeMap::from([(key([1]), c.one())]));
    assert_eq!(
        direct_two.0.terms().len(),
        2,
        "the dotted two-loop fixture must exercise connected and product terms"
    );

    // Independently replay the only product conversion needed by this
    // reference: I(0,1,1) = I_1(1)^2 under q1 -> k2-k1, q2 -> k1,
    // whose native verified determinant is -1. No product reducer is used
    // to build the expected sum, and no terminal values are substituted.
    let source = momentum_family(two.family(), &[1, 2]).unwrap();
    let destination = physical_family(&[S.k.call(1), S.k.call(2)], 2).unwrap();
    assert!(
        source
            .mapping_to(&destination, &[S.k.call(2) - S.k.call(1), S.k.call(1)])
            .unwrap()
            .is_some()
    );
    let three_key = (
        three.family().fingerprint().to_owned(),
        three_target.powers().to_vec(),
    );
    let one_key = (one.family().fingerprint().to_owned(), vec![1]);
    let mut expected = BTreeMap::new();
    for (terminal, coefficient) in direct_two.0.terms() {
        let mut factors = vec![three_key.clone()];
        match terminal.powers() {
            [1, 1, 1] => factors.push((two.family().fingerprint().to_owned(), vec![1, 1, 1])),
            [0, 1, 1] => factors.extend([one_key.clone(), one_key.clone()]),
            powers => panic!("unexpected shipped two-loop terminal {powers:?}"),
        }
        factors.sort();
        expected.insert(
            factors,
            c.try_mul(&c.integer(7), coefficient, Default::default())
                .unwrap(),
        );
    }
    // Native family-coordinate proofs may retain nonzero rational constants
    // (for example basis determinants), whereas closed guarded application
    // omits these vacuous conditions. Validate every guard before comparing
    // the remaining parameter-dependent obligations; never discard zero.
    let nonconstant_guard = |condition: &&Coefficient| {
        c.validate_with_limits(condition, Default::default())
            .unwrap();
        assert!(
            !condition.is_zero(),
            "a retained nonzero guard cannot vanish"
        );
        !(condition.numerator.is_constant() && condition.denominator.is_constant())
    };
    let expected_guards = direct_three
        .1
        .iter()
        .chain(&direct_two.1)
        .chain(&direct_one.1)
        .filter(nonconstant_guard)
        .collect::<HashSet<_>>();
    let input = vec![(
        vec![
            VacuumFactor {
                family: three_family.clone(),
                integral: three_target,
            },
            VacuumFactor {
                family: two_family,
                integral: two_target,
            },
        ],
        c.integer(7),
    )];
    let mut reducer = VacuumProductReducer::try_new(
        ShippedFactorProvider::default(),
        three_family,
        Default::default(),
    )
    .unwrap();
    let actual = reducer.reduce_products(input.clone(), []).unwrap();
    assert_eq!(product_map(&actual), expected);
    assert_eq!(
        actual
            .nonzero_conditions()
            .iter()
            .filter(nonconstant_guard)
            .collect::<HashSet<_>>(),
        expected_guards
    );
    assert_eq!(
        actual.statistics().provider_calls,
        3,
        "one cached provider call for each distinct 3L, 2L and 1L factor"
    );
    assert!(actual.terms().keys().all(|factors| {
        factors
            .iter()
            .map(|factor| actual.families()[factor.family_fingerprint()].loop_count())
            .sum::<usize>()
            == 5
    }));
    let warm = reducer.reduce_products(input.clone(), []).unwrap();
    assert_eq!(product_map(&warm), expected);
    assert_eq!(warm.nonzero_conditions(), actual.nonzero_conditions());
    assert_eq!(warm.statistics().provider_calls, 0);
    reducer.clear_cache().unwrap();
    assert_eq!(reducer.cache_statistics().cached_factors, 0);
    let cold_again = reducer.reduce_products(input, []).unwrap();
    assert_eq!(product_map(&cold_again), expected);
    assert_eq!(cold_again.statistics().provider_calls, 3);
}
