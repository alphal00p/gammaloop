//! Representative 2/3-loop native cache parity, not numerical master work.
use super::{artifact::shipped_family_for_loop_count, cache::*};
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions, Vakint,
    VakintSettings, vakint_parse,
};
use rustred::{
    family::IntegralKey,
    persistence::application_cache::ApplicationCacheKey,
    reduction::{ClosedUnitMassProgram, Reducer},
};
use rustred_app::{CachePreparationStatus, SavedApplicationCacheBundle};
use std::{
    collections::HashSet,
    path::{Path, PathBuf},
    sync::atomic::{AtomicBool, AtomicUsize, Ordering},
    time::Instant,
};
use symbolica::atom::{Atom, AtomCore};

struct Scratch(PathBuf);
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn settings() -> VakintSettings {
    VakintSettings {
        form_exe_path: "/lower-cache-must-not-use-form".into(),
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
fn input(loops: usize, mass: i64, power: i64, pinch: bool) -> Atom {
    let expression = match (loops, pinch) {
        (2, false) => format!("topo(I2L({mass},{power},1,1))"),
        (2, true) => format!("topo(I2L_pinch_3({mass},{power},1,0))"),
        (3, false) => format!("topo(I3L({mass},{power},1,1,1,1,1))"),
        (3, true) => format!("topo(I3L_pinch_1_6({mass},0,1,1,1,1,0))"),
        _ => panic!("fixture loop count"),
    };
    vakint_parse!(&expression).unwrap()
}
fn parity(left: Atom, right: Atom, label: &str) {
    assert!(
        (left - right).expand().together().is_zero(),
        "{label}: cached/uncached mismatch"
    );
}
fn fresh_reader(path: &Path, loops: usize) {
    let engine = Vakint::new().unwrap();
    let reference = Vakint::new().unwrap();
    let settings = settings();
    let started = Instant::now();
    engine.load_cache(loops, path, &Default::default()).unwrap();
    println!(
        "LOWER_CACHE loops={loops} stage=load seconds={:.9}",
        started.elapsed().as_secs_f64()
    );
    assert_eq!(
        engine
            .cache_statistics(loops)
            .unwrap()
            .disk
            .unwrap()
            .chunks_read,
        0
    );
    for pinch in [false, true] {
        for mass in [1, 4] {
            let integral = input(loops, mass, 2, pinch);
            let started = Instant::now();
            let cached = engine
                .evaluate_integral(&settings, integral.as_view())
                .unwrap();
            let elapsed = started.elapsed();
            let uncached = reference
                .evaluate_integral(&settings, integral.as_view())
                .unwrap();
            parity(cached.clone(), uncached, "prepared root");
            let warm = engine
                .evaluate_integral(&settings, integral.as_view())
                .unwrap();
            parity(cached, warm, "warm root");
            let stats = engine.cache_statistics(loops).unwrap();
            assert_eq!(
                stats.native.completed_reductions, 0,
                "selected prepared roots must not miss"
            );
            assert_eq!(stats.native.work.rule_applications(), 0);
            println!(
                "LOWER_CACHE loops={loops} stage=prepared_public mass_squared={mass} pinch={pinch} seconds={:.9} counters={:?}",
                elapsed.as_secs_f64(),
                stats.native
            );
        }
    }
    let integral = input(loops, 4, 3, false);
    let before = engine.cache_statistics(loops).unwrap();
    let started = Instant::now();
    let held_out = engine
        .evaluate_integral(&settings, integral.as_view())
        .unwrap();
    let elapsed = started.elapsed();
    parity(
        held_out,
        reference
            .evaluate_integral(&settings, integral.as_view())
            .unwrap(),
        "held-out root",
    );
    let after = engine.cache_statistics(loops).unwrap();
    assert_eq!(
        after.native.completed_reductions,
        before.native.completed_reductions + 1
    );
    assert!(after.native.work.rule_applications() > before.native.work.rule_applications());
    println!(
        "LOWER_CACHE loops={loops} stage=held_out_public seconds={:.9} counters={:?}",
        elapsed.as_secs_f64(),
        after.native
    );

    // Inspect an exact prepared root separately. Stored source guards survive
    // independently of the public value/materialization comparison above.
    let snapshot = shipped_family_for_loop_count(loops)
        .unwrap()
        .snapshot()
        .unwrap();
    let program = ClosedUnitMassProgram::new(snapshot.artifact(), snapshot.source_digest());
    let bundle = SavedApplicationCacheBundle::open_generated(
        path,
        std::slice::from_ref(program.binding()),
        Default::default(),
    )
    .unwrap();
    let mut reader = bundle.open_reader(program.binding()).unwrap();
    let mut powers = vec![1; program.artifact().arity()];
    powers[0] = 2;
    let target = IntegralKey::try_new(powers).unwrap();
    let stored = reader
        .lookup(
            &ApplicationCacheKey::guarded_root(&target),
            program.context(),
        )
        .unwrap()
        .unwrap();
    let (reference, guards, _) =
        Reducer::reduce_artifact_unit_mass_guarded(program.artifact(), &target, Default::default())
            .unwrap();
    assert!(
        stored.terms == *reference.terms(),
        "native prepared coefficients differ"
    );
    assert!(
        stored.nonzero_conditions.iter().collect::<HashSet<_>>()
            == guards.iter().collect::<HashSet<_>>(),
        "native prepared guard obligations differ"
    );
}

fn run(loops: usize) {
    const CHILD: &str = "VAKINT_LOWER_CACHE_CHILD";
    if let Some(path) = std::env::var_os(CHILD) {
        fresh_reader(Path::new(&path), loops);
        return;
    }
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let path = std::env::temp_dir().join(format!(
        "vakint-lower-cache-{}-{loops}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&path).unwrap();
    let scratch = Scratch(path);
    let output = scratch.0.join("prepared");
    let engine = Vakint::new().unwrap();
    let options = RustRedCachePreparationOptions {
        workload: RustRedCacheWorkload::RelativePrimer {
            max_extra_dots: 1,
            max_numerator_rank: 0,
        },
        ..Default::default()
    };
    let started = Instant::now();
    let report = engine
        .prepare_cache(
            loops,
            &output,
            &options,
            &AtomicBool::new(false),
            |_| Ok(()),
        )
        .unwrap();
    assert_eq!(
        report.status,
        CachePreparationStatus::Complete,
        "{}",
        report.to_document()
    );
    assert_eq!(report.total_roots, if loops == 2 { 20 } else { 256 });
    println!(
        "LOWER_CACHE loops={loops} stage=prepare total_seconds={:.9} report={}",
        started.elapsed().as_secs_f64(),
        report.to_document()
    );
    let filter = match loops {
        2 => {
            "rustred_evaluation::cache_lower_tests::two_loop_prepared_fresh_cache_preserves_values_and_guards"
        }
        3 => {
            "rustred_evaluation::cache_lower_tests::three_loop_prepared_fresh_cache_preserves_values_and_guards"
        }
        _ => unreachable!(),
    };
    let result = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            filter,
            "--ignored",
            "--nocapture",
            "--test-threads=1",
        ])
        .env(CHILD, &output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "fresh child failed: {}\n{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stdout).contains("1 passed"));
    print!("{}", String::from_utf8_lossy(&result.stdout));
}

#[test]
#[ignore = "explicit native 20-root cache preparation and fresh-process gate"]
fn two_loop_prepared_fresh_cache_preserves_values_and_guards() {
    run(2);
}

#[test]
#[ignore = "explicit native 256-root cache preparation and fresh-process gate"]
fn three_loop_prepared_fresh_cache_preserves_values_and_guards() {
    run(3);
}
