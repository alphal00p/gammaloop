use super::*;
use crate::{Topologies, Vakint};

fn registry() -> Topologies {
    Vakint::initialize_vakint_symbols();
    Topologies::generate_topologies().unwrap()
}

fn accepted() -> AcceptedKey {
    vec![(IntegralKey::try_new([1]).unwrap(), Atom::one())]
}

fn insert(cache: &mut WitnessCache, registry: &Topologies, value: i64) {
    cache.insert(
        registry.registry_identity(),
        Policy::from(&VakintSettings::default()),
        Atom::num(value).as_view(),
        Atom::one().as_view(),
        &ReplacementRules::default(),
        accepted(),
    );
}

#[test]
fn witness_memo_lru_eviction_and_clear_do_not_change_selected_payload() {
    let registry = registry();
    let policy = Policy::from(&VakintSettings::default());
    let mut cache = WitnessCache::new(2, 1024 * 1024);
    insert(&mut cache, &registry, 1);
    insert(&mut cache, &registry, 2);
    let first = cache
        .get(
            &registry.registry_identity(),
            &policy,
            Atom::num(1).as_view(),
            Atom::one().as_view(),
        )
        .unwrap();
    assert_eq!(first.accepted, accepted());
    insert(&mut cache, &registry, 3);
    assert!(
        cache
            .get(
                &registry.registry_identity(),
                &policy,
                Atom::num(2).as_view(),
                Atom::one().as_view()
            )
            .is_none()
    );
    assert!(
        cache
            .get(
                &registry.registry_identity(),
                &policy,
                Atom::num(1).as_view(),
                Atom::one().as_view()
            )
            .is_some()
    );
    assert_eq!(cache.statistics().retained_entries, 2);
    assert_eq!(cache.statistics().evictions, 1);
    cache.clear();
    assert_eq!(cache.statistics().retained_entries, 0);
    assert_eq!(cache.statistics().tracked_atom_bytes, 0);
    // A reader retaining an Arc remains valid after eviction/clear.
    assert_eq!(first.accepted, accepted());
}

#[test]
fn witness_memo_binds_registry_incarnation_and_installed_catalogue() {
    let original = registry();
    let cloned_before_load = original.clone();
    assert_eq!(
        original.registry_identity(),
        cloned_before_load.registry_identity()
    );
    assert_ne!(original.registry_identity(), registry().registry_identity());
    let before = original.registry_identity();
    let catalogue = crate::topologies::FiveLoopCatalogue::prepare().unwrap();
    original.install_five_loop_catalogue(Arc::clone(&catalogue));
    assert_ne!(before, original.registry_identity());
    assert_ne!(
        original.registry_identity(),
        cloned_before_load.registry_identity()
    );
    cloned_before_load.install_five_loop_catalogue(catalogue);
    assert_eq!(
        original.registry_identity(),
        cloned_before_load.registry_identity()
    );
}

#[test]
fn witness_memo_oversize_and_disabled_capacity_skip_without_retaining_atoms() {
    let registry = registry();
    for mut cache in [WitnessCache::new(0, 1024), WitnessCache::new(32, 0)] {
        insert(&mut cache, &registry, 1);
        assert_eq!(cache.statistics().oversized_skips, 1);
        assert_eq!(cache.statistics().retained_entries, 0);
        assert_eq!(cache.statistics().tracked_atom_bytes, 0);
    }
}

#[test]
fn witness_memo_encoded_byte_budget_evicts_without_growing_past_limit() {
    let registry = registry();
    let mut measured = WitnessCache::default();
    // Make Atom payload larger than the separately bounded epsilon name, so
    // this control exercises aggregate payload eviction, not name oversize.
    insert(&mut measured, &registry, i64::MAX);
    let first_bytes = measured.statistics().tracked_atom_bytes;
    measured.clear();
    insert(&mut measured, &registry, i64::MAX - 1);
    let one_entry_bytes = first_bytes.max(measured.statistics().tracked_atom_bytes);
    assert!(one_entry_bytes > 0);
    assert!(one_entry_bytes >= Policy::from(&VakintSettings::default()).epsilon.len());
    let mut cache = WitnessCache::new(32, one_entry_bytes);
    insert(&mut cache, &registry, i64::MAX);
    insert(&mut cache, &registry, i64::MAX - 1);
    assert_eq!(cache.statistics().retained_entries, 1);
    assert_eq!(cache.statistics().evictions, 1);
    assert!(cache.statistics().tracked_atom_bytes <= one_entry_bytes);
}

#[test]
fn witness_memo_exact_inputs_and_policy_mutations_never_alias() {
    let registry = registry();
    let defaults = VakintSettings::default();
    let mut cache = WitnessCache::default();
    insert(&mut cache, &registry, 1);
    let mut variants = Vec::new();
    let mut changed = defaults.clone();
    changed.rustred_momentum_match_limits.max_verified_witnesses = 0;
    variants.push(changed);
    let mut changed = defaults.clone();
    changed.rustred_scalar_numerator_limits.max_input_nodes = 0;
    variants.push(changed);
    let mut changed = defaults.clone();
    changed.epsilon_symbol.push('x');
    variants.push(changed);
    let mut changed = defaults.clone();
    changed.verify_numerator_identification = !changed.verify_numerator_identification;
    variants.push(changed);
    let mut changed = defaults.clone();
    changed.allow_unknown_integrals = !changed.allow_unknown_integrals;
    variants.push(changed);
    let mut changed = defaults.clone();
    changed.use_dot_product_notation = !changed.use_dot_product_notation;
    variants.push(changed);
    let mut changed = defaults.clone();
    changed.tensor_reduction_method = TensorReductionMethod::AlphaLoop;
    variants.push(changed);
    for settings in variants {
        assert!(
            cache
                .get(
                    &registry.registry_identity(),
                    &Policy::from(&settings),
                    Atom::one().as_view(),
                    Atom::one().as_view()
                )
                .is_none()
        );
    }
    assert!(
        cache
            .get(
                &registry.registry_identity(),
                &Policy::from(&defaults),
                Atom::num(2).as_view(),
                Atom::one().as_view()
            )
            .is_none()
    );
    assert!(
        cache
            .get(
                &registry.registry_identity(),
                &Policy::from(&defaults),
                Atom::one().as_view(),
                Atom::num(2).as_view()
            )
            .is_none()
    );
    let valid = cache
        .get(
            &registry.registry_identity(),
            &Policy::from(&defaults),
            Atom::one().as_view(),
            Atom::one().as_view(),
        )
        .unwrap();
    cache.invalidate(&valid);
    assert_eq!(cache.statistics().invalidations, 1);
    assert_eq!(cache.statistics().retained_entries, 0);
}

fn same_witness(left: &ReplacementRules, right: &ReplacementRules) {
    assert_eq!(
        left.edge_ids_canonical_to_input_map,
        right.edge_ids_canonical_to_input_map
    );
    assert_eq!(
        left.canonical_expression_substitutions,
        right.canonical_expression_substitutions
    );
    assert_eq!(left.numerator_substitutions, right.numerator_substitutions);
    let left = left.canonical_topology.get_integral();
    let right = right.canonical_topology.get_integral();
    assert_eq!(left.canonical_expression, right.canonical_expression);
    assert_eq!(left.parent_routing, right.parent_routing);
    assert_eq!(left.short_expression, right.short_expression);
}

#[test]
#[ignore = "requires the R0/D10 saved package admitting the scalar banana; no IBP generation"]
fn external_package_witness_memo_preserves_native_minimum_and_policy() {
    use crate::{EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions, vakint_parse};
    let path = std::env::var_os("VAKINT_SAVED_APPLICATION").expect("external application required");
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/witness-memo-test-must-not-use-form".into(),
        ..Default::default()
    };
    let engine = Vakint::new().unwrap();
    let mut early_clone = engine.clone();
    assert!(
        early_clone
            .topologies
            .five_loop_catalogue_statistics()
            .is_none()
    );
    let engine = engine
        .with_rustred_application(&mut settings, path)
        .unwrap();
    let input = vakint_parse!("topo(prop(1,edge(1,2),k(1),1,1)*prop(2,edge(1,2),k(2),1,1)*prop(3,edge(1,2),k(3),1,1)*prop(4,edge(1,2),k(4),1,1)*prop(5,edge(1,2),k(5),1,1)*prop(6,edge(1,2),-k(1)-k(2)-k(3)-k(4)-k(5),1,4))").unwrap();
    let numerator = Atom::one();
    let uncached = engine
        .match_rustred_input_uncached(&settings, input.as_view(), numerator.as_view())
        .unwrap()
        .unwrap();
    let expected = engine
        .rustred_matched_input_key(&settings, numerator.as_view(), &uncached)
        .unwrap()
        .expect("actual native admission");
    let first = engine
        .match_rustred_input(&settings, input.as_view(), numerator.as_view())
        .unwrap()
        .unwrap();
    same_witness(&uncached, &first);
    assert_eq!(
        engine
            .rustred_matched_input_key(&settings, numerator.as_view(), &first)
            .unwrap(),
        Some(expected.clone())
    );
    assert_eq!(
        engine
            .rustred_application_witness_cache_statistics()
            .unwrap()
            .retained_entries,
        1
    );

    // A clone created before loading must install the same immutable catalogue
    // before looking up the shared package's memo.
    let cloned = early_clone
        .match_rustred_input(&settings, input.as_view(), numerator.as_view())
        .unwrap()
        .unwrap();
    same_witness(&first, &cloned);
    assert!(
        early_clone
            .topologies
            .five_loop_catalogue_statistics()
            .is_some()
    );
    assert_eq!(
        engine
            .rustred_application_witness_cache_statistics()
            .unwrap()
            .hits,
        1
    );

    let mut strict = settings.clone();
    strict
        .rustred_momentum_match_limits
        .max_orientation_candidates = 0;
    assert!(matches!(
        engine.match_rustred_input(&strict, input.as_view(), numerator.as_view()),
        Err(crate::VakintError::MomentumMatchResourceLimit { .. })
    ));
    let mut strict = settings.clone();
    strict.rustred_scalar_numerator_limits.max_input_nodes = 0;
    assert!(
        engine
            .match_rustred_input(&strict, input.as_view(), numerator.as_view())
            .is_err()
    );
    let mut disabled = settings.clone();
    disabled.rustred_artifact_path = None;
    assert!(
        engine
            .evaluate_rustred_symbolic(&disabled, input.as_view())
            .is_err()
    );
    assert_eq!(
        engine
            .rustred_application_witness_cache_statistics()
            .unwrap()
            .retained_entries,
        1
    );

    // The no-admission short-form fallback is useful for the final error but
    // must never become a successful witness memo entry.
    let unsupported = vakint_parse!("topo(I5L_32745(1,1,1,1,1,1,1,1,1,1,1,1,1))").unwrap();
    let fallback = engine
        .match_rustred_input(&settings, unsupported.as_view(), numerator.as_view())
        .unwrap()
        .unwrap();
    assert!(
        engine
            .rustred_matched_input_key(&settings, numerator.as_view(), &fallback)
            .unwrap()
            .is_none()
    );
    assert_eq!(
        engine
            .rustred_application_witness_cache_statistics()
            .unwrap()
            .retained_entries,
        1
    );

    std::thread::scope(|scope| {
        let jobs = (0..2)
            .map(|_| {
                scope.spawn(|| {
                    engine
                        .match_rustred_input(&settings, input.as_view(), numerator.as_view())
                        .unwrap()
                        .unwrap()
                })
            })
            .collect::<Vec<_>>();
        for job in jobs {
            same_witness(&first, &job.join().unwrap());
        }
    });
    assert_eq!(
        engine
            .rustred_application_witness_cache_statistics()
            .unwrap()
            .hits,
        3
    );
    // Public registry replacement must miss even when its contents are equal.
    early_clone.topologies = registry();
    let replaced = early_clone
        .match_rustred_input(&settings, input.as_view(), numerator.as_view())
        .unwrap()
        .unwrap();
    same_witness(&first, &replaced);
    assert_eq!(
        engine
            .rustred_application_witness_cache_statistics()
            .unwrap()
            .retained_entries,
        2
    );
    assert_eq!(
        engine
            .rustred_application_witness_cache_statistics()
            .unwrap()
            .hits,
        3
    );
    // Deliberately corrupt only the remembered comparison key, not the native
    // artifact. Hit revalidation must invalidate and recompute the same exact
    // minimum through the ordinary matcher instead of trusting stale metadata.
    let application = engine.saved_application(&settings).unwrap();
    {
        let registry = engine.topologies.registry_identity();
        let mut cache = application.witnesses.lock().unwrap();
        let entry = cache
            .entries
            .iter_mut()
            .find(|entry| entry.registry == registry)
            .unwrap();
        Arc::get_mut(entry)
            .expect("no memo readers remain")
            .accepted[0]
            .1 = Atom::num(999);
    }
    let repaired = engine
        .match_rustred_input(&settings, input.as_view(), numerator.as_view())
        .unwrap()
        .unwrap();
    same_witness(&first, &repaired);
    assert_eq!(
        engine
            .rustred_application_witness_cache_statistics()
            .unwrap()
            .invalidations,
        1
    );
    println!("native_minimum_witness_exact_parity=passed");
    println!(
        "witness_memo={:?}",
        engine
            .rustred_application_witness_cache_statistics()
            .unwrap()
    );
    engine.clear_rustred_application_cache().unwrap();
    let cleared = engine
        .rustred_application_witness_cache_statistics()
        .unwrap();
    assert_eq!(cleared.retained_entries, 0);
    assert_eq!(cleared.tracked_atom_bytes, 0);

    // A panic in optional observational memo state must not prevent the
    // unchanged exact matching path from succeeding afterward.
    let application = engine.saved_application(&settings).unwrap();
    let poisoned = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _guard = application.witnesses.lock().unwrap();
        panic!("intentional witness-memo poisoning control");
    }));
    assert!(poisoned.is_err());
    let without_memo = engine
        .match_rustred_input(&settings, input.as_view(), numerator.as_view())
        .unwrap()
        .unwrap();
    same_witness(&first, &without_memo);
    println!("poisoned_optional_witness_memo_bypassed=passed");
}
