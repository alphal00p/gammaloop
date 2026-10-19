//! Explicit first-stage four-parent cache acceptance, not a large primer or
//! a matched performance benchmark. No FORM/FMFT or new master computation.
use super::{cache::*, four_loop::cache as shipped};
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions, Vakint,
    VakintSettings,
};
use rustred::{
    family::IntegralKey, persistence::application_cache::ApplicationCacheKey,
    reduction::CompletedCachePolicy,
};
use rustred_app::{CachePreparationStatus, SavedApplicationCacheBundle};
use std::{
    path::{Path, PathBuf},
    sync::atomic::{AtomicBool, AtomicUsize, Ordering},
    time::Instant,
};
use symbolica::atom::{Atom, AtomCore};

const CHILD: &str = "VAKINT_FOUR_LOOP_CACHE_TEST_CHILD";
const OUTPUT: &str = "VAKINT_FOUR_LOOP_CACHE_TEST_OUTPUT";
const FILTER: &str = "rustred_evaluation::cache_four_loop_tests::four_loop_first38_fresh_cache_preserves_public_values_guards_and_limits";

fn settings(substitute_masters: bool) -> VakintSettings {
    VakintSettings {
        form_exe_path: "/four-cache-acceptance-must-not-invoke-form".into(),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: if substitute_masters {
                    crate::RustRedMasterMode::Numerical
                } else {
                    crate::RustRedMasterMode::Retain
                },
            },
        )]),
        number_of_terms_in_epsilon_expansion: 5,
        run_time_decimal_precision: 32,
        use_dot_product_notation: true,
        ..Default::default()
    }
}

fn parent_powers() -> Vec<[i64; 10]> {
    let masks = shipped::parent_masks(&super::cache_tests::external_engine()).unwrap();
    assert_eq!(masks.len(), 4);
    assert_eq!(
        masks
            .iter()
            .map(|(_, mask)| mask.iter().filter(|&&p| p).count())
            .collect::<Vec<_>>(),
        [9, 8, 8, 9]
    );
    masks
        .into_iter()
        .enumerate()
        .map(|(namespace, (actual, mask))| {
            assert_eq!(actual, namespace);
            assert_eq!(mask.len(), 10);
            std::array::from_fn(|axis| i64::from(mask[axis]))
        })
        .collect()
}

fn first38(parents: &[[i64; 10]]) -> Vec<RustRedCacheRoot> {
    let mut roots = Vec::new();
    for (namespace, parent) in parents.iter().enumerate() {
        roots.push(RustRedCacheRoot {
            namespace,
            integral: IntegralKey::try_new(*parent).unwrap(),
        });
        for axis in 0..10 {
            if parent[axis] == 0 {
                continue;
            }
            let mut powers = *parent;
            powers[axis] += 1;
            roots.push(RustRedCacheRoot {
                namespace,
                integral: IntegralKey::try_new(powers).unwrap(),
            });
        }
    }
    assert_eq!(roots.len(), 38);
    roots
}

fn exact(left: Atom, right: Atom, label: &str) {
    assert!(
        (left - right).expand().together().is_zero(),
        "four-cache exact mismatch: {label}"
    );
}

fn ordinary_statistics() -> Vec<rustred::solver::CandidateStatistics> {
    (0..4).map(shipped::ordinary_statistics).collect()
}

/// A fresh Vakint object alone does not isolate the ordinary static reducers.
/// Check their full native work+retained census around each cached call, before
/// running any ordinary reference which may legitimately warm that static cache.
fn cached_public(
    engine: &Vakint,
    settings: &VakintSettings,
    input: &Atom,
    label: &str,
) -> (Atom, f64) {
    let ordinary_before = ordinary_statistics();
    let started = Instant::now();
    let value = engine
        .evaluate_integral(settings, input.as_view())
        .unwrap_or_else(|error| panic!("{label}: {error}"));
    let seconds = started.elapsed().as_secs_f64();
    assert_eq!(
        ordinary_statistics(),
        ordinary_before,
        "{label}: cache path mutated ordinary shared reducers"
    );
    (value, seconds)
}

fn numerical_parity(
    engine: &Vakint,
    reference: &Vakint,
    namespace: usize,
    powers: &[i64; 10],
    mass: i64,
) {
    let settings = settings(true);
    let input = shipped::test_input(namespace, powers, mass);
    let (cached, seconds) = cached_public(
        engine,
        &settings,
        &input,
        "existing master-table substitution",
    );
    let ordinary = reference
        .evaluate_integral(&settings, input.as_view())
        .unwrap();
    exact(
        cached.clone(),
        ordinary.clone(),
        "master-substituted public result",
    );
    let params = crate::params_from_f64(
        &[("mursq".into(), 1.0)].into_iter().collect(),
        settings.run_time_decimal_precision,
    );
    let cached_numeric = engine
        .numerical_evaluation(
            &settings,
            cached.as_view(),
            &params,
            &Default::default(),
            None,
        )
        .unwrap();
    let ordinary_numeric = reference
        .numerical_evaluation(
            &settings,
            ordinary.as_view(),
            &params,
            &Default::default(),
            None,
        )
        .unwrap();
    assert!(!cached_numeric.0.get_epsilon_coefficients().is_empty());
    let errors = match (&cached_numeric.1, &ordinary_numeric.1) {
        (Some(left), Some(right)) => Some(left.aggregate_errors(right)),
        (Some(value), None) | (None, Some(value)) => Some(value.clone()),
        (None, None) => None,
    };
    let (matches, detail) =
        cached_numeric
            .0
            .does_approx_match(&ordinary_numeric.0, errors.as_ref(), 1e-24, 1.0);
    assert!(matches, "existing-table numerical parity: {detail}");
    println!(
        "CACHE_4L stage=master_substitution namespace={namespace} mass_squared={mass} public_seconds={seconds:.9} table_source=existing_shipped no_new_master_computation=true"
    );
}

fn fresh_reader(path: &Path) {
    let started = Instant::now();
    let engine = super::cache_tests::external_engine();
    println!(
        "CACHE_4L stage=fresh_engine_new seconds={:.9}",
        started.elapsed().as_secs_f64()
    );
    assert!(!engine.cache_statistics(4).unwrap().attached);
    let started = Instant::now();
    let programs = shipped::programs(&engine).unwrap();
    println!(
        "CACHE_4L stage=fresh_program_load seconds={:.9}",
        started.elapsed().as_secs_f64()
    );
    let parents = parent_powers();
    let roots = first38(&parents);
    let ordinary_before = ordinary_statistics();
    let started = Instant::now();
    engine.load_cache(4, path, &Default::default()).unwrap();
    println!(
        "CACHE_4L stage=attach_index seconds={:.9} program_load_excluded=true",
        started.elapsed().as_secs_f64()
    );
    assert_eq!(ordinary_statistics(), ordinary_before);
    let attached = engine.cache_statistics(4).unwrap();
    assert!(attached.attached);
    assert_eq!(attached.namespaces.len(), 4);
    for (namespace, program) in programs.iter().enumerate() {
        assert_eq!(&attached.namespaces[namespace].binding, program.binding());
    }
    assert!(
        attached
            .namespaces
            .iter()
            .all(|ns| ns.disk.as_ref().is_some_and(|io| io.chunks_read == 0))
    );
    let reference = super::cache_tests::external_engine();
    assert!(!reference.cache_statistics(4).unwrap().attached);
    let raw = settings(false);
    for root in &roots {
        let powers: [i64; 10] = root.integral.powers().try_into().unwrap();
        for mass in [1, 4] {
            let label = format!("namespace{} {:?} mass{mass}", root.namespace, powers);
            let input = shipped::test_input(root.namespace, &powers, mass);
            let before = engine.cache_statistics(4).unwrap();
            let (cached, first_seconds) = cached_public(&engine, &raw, &input, &label);
            assert!(
                !cached.is_zero(),
                "prepared nonzero parent/dot was zero: {label}"
            );
            let first = engine.cache_statistics(4).unwrap();
            assert_eq!(
                first.native.completed_reductions, 0,
                "prepared root missed: {label}"
            );
            assert_eq!(first.native.work.rule_applications(), 0);
            assert!(
                first.native.backing_hits + first.native.ram_hits
                    > before.native.backing_hits + before.native.ram_hits
            );
            let started = Instant::now();
            let ordinary = reference.evaluate_integral(&raw, input.as_view()).unwrap();
            let reference_seconds = started.elapsed().as_secs_f64();
            exact(cached.clone(), ordinary, &label);
            let (warm, warm_seconds) = cached_public(&engine, &raw, &input, &label);
            exact(cached, warm, "RAM repeat");
            let after = engine.cache_statistics(4).unwrap();
            assert!(after.native.ram_hits > first.native.ram_hits);
            println!(
                "CACHE_4L stage=prepared_public namespace={} powers={powers:?} mass_squared={mass} first_loaded_seconds={first_seconds:.9} ram_repeat_seconds={warm_seconds:.9} ordinary_reference_seconds={reference_seconds:.9} native={:?} timing_scope=shared_host_unpaired_reference_may_reuse_ordinary_cache",
                root.namespace, after.native
            );
        }
    }
    let hits = engine.cache_statistics(4).unwrap();
    assert!(hits.namespaces.iter().all(|ns| ns.native.backing_hits > 0));
    println!("CACHE_4L prepared_hit_counters={hits:?}");

    // Scalar numerator cancellation and its physical pinch are separate from
    // the 38 prepared parent/dot roots: any native misses are reported honestly.
    for (namespace, parent) in parents.iter().enumerate() {
        for mass in [1, 4] {
            let input = (shipped::test_denominator(namespace, 6, mass)
                * shipped::test_input(namespace, parent, mass))
            .expand();
            let mut pinched = *parent;
            pinched[6] = 0;
            let pinch = shipped::test_input(namespace, &pinched, mass);
            let (numerator, seconds) =
                cached_public(&engine, &raw, &input, "expanded denominator cancellation");
            let (pinch_value, _) = cached_public(&engine, &raw, &pinch, "physical pinch");
            assert!(
                !pinch_value.is_zero(),
                "vacuous cancellation fixture namespace={namespace}"
            );
            exact(
                numerator.clone(),
                pinch_value,
                "denominator cancellation equals pinch",
            );
            exact(
                numerator,
                reference.evaluate_integral(&raw, input.as_view()).unwrap(),
                "cancellation ordinary parity",
            );
            println!(
                "CACHE_4L stage=cancellation namespace={namespace} mass_squared={mass} seconds={seconds:.9} counters={:?}",
                engine.cache_statistics(4).unwrap()
            );
        }
    }

    for (namespace, parent) in parents.iter().enumerate() {
        let mut powers = *parent;
        powers[0] = 3;
        let target = IntegralKey::try_new(powers).unwrap();
        assert!(
            !roots
                .iter()
                .any(|root| root.namespace == namespace && root.integral == target)
        );
        let (native, _) = engine
            .four_loop_data()
            .unwrap()
            .four_loop
            .as_ref()
            .unwrap()
            .native_program(namespace)
            .unwrap();
        assert!(!native.terminals().contains(&target));
        let input = shipped::test_input(namespace, &powers, 4);
        let before = engine.cache_statistics(4).unwrap();
        let (value, seconds) = cached_public(&engine, &raw, &input, "held-out cubed line");
        assert!(!value.is_zero());
        let after = engine.cache_statistics(4).unwrap();
        assert!(after.native.completed_reductions > before.native.completed_reductions);
        assert!(after.native.work.rule_applications() > before.native.work.rule_applications());
        exact(
            value.clone(),
            reference.evaluate_integral(&raw, input.as_view()).unwrap(),
            "held-out ordinary parity",
        );
        let (warm, warm_seconds) = cached_public(&engine, &raw, &input, "held-out RAM repeat");
        exact(value, warm, "held-out warm parity");
        println!(
            "CACHE_4L stage=held_out_cube namespace={namespace} first_seconds={seconds:.9} ram_repeat_seconds={warm_seconds:.9} counters={after:?} candidate_descendant_reuse_claim=false"
        );
    }
    for namespace in [0, 1] {
        for mass in [1, 4] {
            numerical_parity(&engine, &reference, namespace, &parents[namespace], mass);
        }
    }
    let before_clear = engine.cache_statistics(4).unwrap();
    engine.clear_cache(4).unwrap();
    assert_eq!(
        engine.cache_statistics(4).unwrap().native.cached_integrals,
        0
    );
    let _ = cached_public(
        &engine,
        &raw,
        &shipped::test_input(0, &parents[0], 1),
        "clear then backing reload",
    );
    assert!(
        engine.cache_statistics(4).unwrap().native.backing_hits > before_clear.native.backing_hits
    );
    println!(
        "CACHE_4L final_counters={:?}",
        engine.cache_statistics(4).unwrap()
    );
}

struct Scratch(PathBuf);
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
#[ignore = "explicit 38-root four-family native preparation and fresh-process acceptance"]
fn four_loop_first38_fresh_cache_preserves_public_values_guards_and_limits() {
    if let Some(path) = std::env::var_os(CHILD) {
        fresh_reader(Path::new(&path));
        return;
    }
    let started = Instant::now();
    let engine = super::cache_tests::external_engine();
    println!(
        "CACHE_4L stage=engine_new seconds={:.9}",
        started.elapsed().as_secs_f64()
    );
    let started = Instant::now();
    let programs = shipped::programs(&engine).unwrap();
    println!(
        "CACHE_4L stage=program_load seconds={:.9}",
        started.elapsed().as_secs_f64()
    );
    let parents = parent_powers();
    let roots = first38(&parents);
    let (output, _scratch) = if let Some(output) = std::env::var_os(OUTPUT) {
        (PathBuf::from(output), None)
    } else {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "vakint-four-cache-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&path).unwrap();
        (path.join("prepared"), Some(Scratch(path)))
    };
    assert!(
        !output.exists(),
        "first-stage output must be new; preserve any previous receipt"
    );
    let options = RustRedCachePreparationOptions {
        workload: RustRedCacheWorkload::Explicit(roots.clone()),
        ..Default::default()
    };
    let ordinary_before = ordinary_statistics();
    let started = Instant::now();
    let report = engine
        .prepare_cache(4, &output, &options, &AtomicBool::new(false), |_| Ok(()))
        .unwrap();
    println!(
        "CACHE_4L stage=prepare seconds={:.9} program_load_excluded=true report={}",
        started.elapsed().as_secs_f64(),
        report.to_document()
    );
    assert_eq!(
        report.status,
        CachePreparationStatus::Complete,
        "{}",
        report.to_document()
    );
    assert_eq!(report.total_roots, 38);
    assert_eq!(report.completed_roots, 38);
    assert_eq!(report.newly_stored_roots, 38);
    assert_eq!(report.reused_roots, 0);
    assert_eq!(
        ordinary_statistics(),
        ordinary_before,
        "preparation changed ordinary native caches"
    );
    assert!(!engine.cache_statistics(4).unwrap().attached);

    // Strong native record controls outside all timings: at least one dotted
    // nonterminal per namespace performs actual cold native work and has a
    // nonzero result equal to the prepared exact terms and ordered guards.
    let bindings = programs
        .iter()
        .map(|program| program.binding().clone())
        .collect::<Vec<_>>();
    let resolved = options.storage.resolve(4, Default::default()).unwrap();
    let bundle =
        SavedApplicationCacheBundle::open_generated(&output, &bindings, resolved.storage_limits)
            .unwrap();
    for (namespace, program) in programs.iter().enumerate() {
        let (native, _) = engine
            .four_loop_data()
            .unwrap()
            .four_loop
            .as_ref()
            .unwrap()
            .native_program(namespace)
            .unwrap();
        let selected = roots
            .iter()
            .find(|root| {
                root.namespace == namespace
                    && root.integral.powers().contains(&2)
                    && !native.terminals().contains(&root.integral)
            })
            .expect("each family needs a nonterminal dotted control");
        let limits = rustred::reduction::ReductionLimits {
            max_cached_coefficient_bytes: resolved.retained_coefficient_bytes / programs.len(),
            ..options.native
        };
        let mut cache = program.new_cache(limits, CompletedCachePolicy::EvictLeastRecentlyUsed);
        let reference = program
            .reduce_guarded_root(&selected.integral, &mut cache, limits, None)
            .unwrap();
        assert!(
            reference.work.rule_applications() > 0,
            "namespace{namespace} control is not a native application"
        );
        assert!(!reference.record.terms.is_empty());
        let stored = bundle
            .open_reader(program.binding())
            .unwrap()
            .lookup(
                &ApplicationCacheKey::guarded_root(&selected.integral),
                program.context(),
            )
            .unwrap()
            .unwrap();
        program.validate_completed_record(&stored, limits).unwrap();
        assert_eq!(stored, reference.record);
        println!(
            "CACHE_4L native_record_control namespace={namespace} target={:?} rule_applications={} guards={} binding={:?}",
            selected.integral.powers(),
            reference.work.rule_applications(),
            stored.nonzero_conditions.len(),
            program.binding()
        );
    }
    assert_eq!(ordinary_statistics(), ordinary_before);
    let result = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            FILTER,
            "--ignored",
            "--nocapture",
            "--test-threads=1",
        ])
        .env(CHILD, &output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "fresh reader failed:\n{}\n{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        String::from_utf8_lossy(&result.stdout).contains("1 passed"),
        "exact fresh-child filter must execute one test"
    );
    print!("{}", String::from_utf8_lossy(&result.stdout));
}
