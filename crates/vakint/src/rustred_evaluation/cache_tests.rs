//! Real shipped-program smoke: no rule generation or numerical master values.
use super::cache::*;
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions, Vakint,
    VakintSettings, vakint_parse,
};
use rustred::family::IntegralKey;
use rustred_app::CachePreparationStatus;
use std::{
    path::{Path, PathBuf},
    sync::atomic::{AtomicBool, AtomicUsize, Ordering},
};
use symbolica::atom::AtomCore;

fn settings() -> VakintSettings {
    VakintSettings {
        form_exe_path: "/native-cache-test-must-not-use-form".into(),
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                substitute_masters: false,
            },
        )]),
        use_dot_product_notation: true,
        ..Default::default()
    }
}
fn input(power: i64) -> symbolica::atom::Atom {
    vakint_parse!(&format!("topo(I1L(4,{power}))")).unwrap()
}
fn assert_same(left: symbolica::atom::Atom, right: symbolica::atom::Atom) {
    assert!((left - right).expand().together().is_zero());
}
fn fresh_reader(path: &Path) {
    let engine = Vakint::new().unwrap();
    let reference = Vakint::new().unwrap();
    let settings = settings();
    assert!(!engine.cache_statistics(1).unwrap().attached);
    engine.load_cache(1, path, &Default::default()).unwrap();
    assert_eq!(
        engine
            .cache_statistics(1)
            .unwrap()
            .disk
            .unwrap()
            .chunks_read,
        0
    );
    let cached = engine
        .evaluate_integral(&settings, input(4).as_view())
        .unwrap();
    assert_same(
        cached.clone(),
        reference
            .evaluate_integral(&settings, input(4).as_view())
            .unwrap(),
    );
    let first = engine.cache_statistics(1).unwrap();
    assert_eq!(first.native.backing_hits, 1);
    assert_eq!(first.native.completed_reductions, 0);
    assert_eq!(first.native.work.rule_applications(), 0);
    assert_same(
        cached,
        engine
            .evaluate_integral(&settings, input(4).as_view())
            .unwrap(),
    );
    assert_eq!(engine.cache_statistics(1).unwrap().native.ram_hits, 1);
    let held_out = engine
        .evaluate_integral(&settings, input(5).as_view())
        .unwrap();
    assert_same(
        held_out,
        reference
            .evaluate_integral(&settings, input(5).as_view())
            .unwrap(),
    );
    let missed = engine.cache_statistics(1).unwrap();
    assert_eq!(missed.native.completed_reductions, 1);
    assert!(missed.native.work.rule_applications() > 0);
    engine.clear_cache(1).unwrap();
    assert_eq!(
        engine.cache_statistics(1).unwrap().native.cached_integrals,
        0
    );
    let _ = engine
        .evaluate_integral(&settings, input(4).as_view())
        .unwrap();
    assert_eq!(engine.cache_statistics(1).unwrap().native.backing_hits, 2);
    // Loading without an explicit override retains the active native RAM cap.
    let mut strict = RustRedCacheLoadOptions::default();
    let mut native = rustred::reduction::ReductionLimits::default();
    native.max_cached_coefficient_bytes = 1;
    strict.native = Some(native);
    engine.load_cache(1, path, &strict).unwrap();
    engine.load_cache(1, path, &Default::default()).unwrap();
    assert!(matches!(
        engine.evaluate_integral(&settings, input(4).as_view()),
        Err(crate::VakintError::RustRedEvaluation(super::RustRedEvaluationError::CacheOperation { source }))
            if source.kind() == rustred_app::AppErrorKind::Limit
    ));
}

struct Scratch(PathBuf);
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn one_loop_preparation_and_fresh_process_hits_preserve_public_evaluation() {
    const CHILD: &str = "VAKINT_ONE_LOOP_CACHE_TEST_CHILD";
    if let Some(path) = std::env::var_os(CHILD) {
        fresh_reader(Path::new(&path));
        return;
    }
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let scratch_path = std::env::temp_dir().join(format!(
        "vakint-cache-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&scratch_path).unwrap();
    let scratch = Scratch(scratch_path);
    let output = scratch.0.join("prepared");
    let engine = Vakint::new().unwrap();
    let options = RustRedCachePreparationOptions {
        workload: RustRedCacheWorkload::Explicit(
            (0..=4)
                .map(|power| RustRedCacheRoot {
                    namespace: 0,
                    integral: IntegralKey::try_new([power]).unwrap(),
                })
                .collect(),
        ),
        ..Default::default()
    };
    let report = engine
        .prepare_cache(1, &output, &options, &AtomicBool::new(false), |_| Ok(()))
        .unwrap();
    assert_eq!(
        report.status,
        CachePreparationStatus::Complete,
        "{:?}",
        report.to_document()
    );
    assert_eq!(report.completed_roots, 5);
    assert!(
        !engine.cache_statistics(1).unwrap().attached,
        "preparation must not implicitly attach a cache"
    );
    let result = std::process::Command::new(std::env::current_exe().unwrap())
        .args(["--exact", "rustred_evaluation::cache_tests::one_loop_preparation_and_fresh_process_hits_preserve_public_evaluation",
            "--nocapture", "--test-threads=1"])
        .env(CHILD, &output).output().unwrap();
    assert!(
        result.status.success(),
        "fresh cache reader failed: {}\n{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stdout).contains("1 passed"));
}
