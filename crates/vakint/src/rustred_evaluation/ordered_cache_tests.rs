//! Ordered step cache: cold and warm evaluations are identical, persisted
//! caches replay every step, foreign caches are rejected.
use std::path::Path;

use symbolica::atom::Atom;

use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, OrderedCacheAutosave,
    OrderedCacheMode, OrderedCacheSettings, RustRedEvaluationOptions, RustRedMasterMode, Vakint,
    VakintSettings, vakint_parse,
};

const THREE_LOOP: [&str; 3] = [
    "topo(I3L(muvsq,2,2,2,1,1,1))",
    "dot(k(1),k(2))*dot(k(2),k(3))*topo(I3L(muvsq,2,1,1,1,1,1))",
    "topo(I3L_pinch_1_6(muvsq,0,1,1,1,2,0))",
];

fn settings(cache: OrderedCacheSettings, arithmetic: &str) -> VakintSettings {
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        use_dot_product_notation: true,
        form_exe_path: "/ordered-cache-tests-must-not-call-form".into(),
        ordered_cache: cache,
        ..Default::default()
    };
    settings
        .rustred_application_options
        .ordered_weight_arithmetic =
        crate::parse_rustred_ordered_weight_arithmetic(arithmetic).unwrap();
    settings
}

fn evaluate_all(engine: &Vakint, settings: &VakintSettings, inputs: &[String]) -> Vec<Atom> {
    inputs
        .iter()
        .map(|input| {
            engine
                .evaluate(settings, vakint_parse!(input).unwrap().as_view())
                .unwrap()
        })
        .collect()
}

fn four_loop_benchmark_inputs(family: &str) -> Vec<String> {
    let cases: serde_json::Value = serde_json::from_str(include_str!(
        "../../tests/fixtures/benchmark_reduction.json"
    ))
    .unwrap();
    cases
        .as_array()
        .unwrap()
        .iter()
        .filter(|case| case["family"] == family)
        .map(|case| case["input"].as_str().unwrap().to_owned())
        .collect()
}

fn opaque_hits(statistics: &crate::OrderedCacheStatistics) -> u64 {
    statistics
        .programs
        .iter()
        .flat_map(|program| &program.native.shards)
        .map(|shard| shard.opaque_hits)
        .sum()
}

fn program_directories(root: &Path) -> Vec<String> {
    let mut names: Vec<_> = std::fs::read_dir(root)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// Uncached, cold (populating), in-process warm and reloaded warm runs return
/// identical retained expressions in both weight arithmetics.
fn assert_cold_and_warm_runs_agree(engine: impl Fn() -> Vakint, inputs: &[String]) {
    let directory = tempfile::tempdir().unwrap();
    let uncached = evaluate_all(
        &engine(),
        &settings(OrderedCacheSettings::off(), "exact"),
        inputs,
    );
    let populating = engine();
    let read_write = settings(
        OrderedCacheSettings::read_write(directory.path()),
        "modular:threads=1",
    );
    assert_eq!(evaluate_all(&populating, &read_write, inputs), uncached);
    let cold = populating.ordered_cache_statistics();
    assert!(cold.inserts > 0 && cold.records > 0, "{cold:?}");
    // Every miss is recorded (opaque factorization steps included).
    assert_eq!(cold.inserts, cold.misses);
    // Autosave wrote every record after each evaluation.
    assert!(
        cold.programs
            .iter()
            .all(|program| program.native.shards.iter().all(|shard| !shard.dirty))
    );
    // Warm in the same engine, in the other arithmetic: every step replays.
    let read_write_exact = settings(OrderedCacheSettings::read_write(directory.path()), "exact");
    assert_eq!(
        evaluate_all(&populating, &read_write_exact, inputs),
        uncached
    );
    let warm = populating.ordered_cache_statistics();
    assert_eq!(warm.misses, cold.misses);
    assert_eq!(warm.inserts, cold.inserts);
    assert!(warm.hits > cold.hits);
    drop(populating);
    // A fresh engine reads the saved shards without recording.
    for arithmetic in ["exact", "modular:threads=1"] {
        let reloaded = engine();
        let read_only = settings(
            OrderedCacheSettings::read_only(directory.path()),
            arithmetic,
        );
        assert_eq!(evaluate_all(&reloaded, &read_only, inputs), uncached);
        let statistics = reloaded.ordered_cache_statistics();
        assert_eq!(statistics.misses, 0, "{statistics:?}");
        assert_eq!(statistics.inserts, 0);
        assert_eq!(statistics.hits, warm.hits - cold.hits);
        assert!(
            statistics
                .programs
                .iter()
                .all(|program| program.mode == OrderedCacheMode::ReadOnly
                    && program.native.loaded_shards > 0)
        );
    }
}

#[test]
fn three_loop_cold_and_warm_ordered_cache_runs_are_identical() {
    let inputs: Vec<String> = THREE_LOOP.iter().map(|input| input.to_string()).collect();
    assert_cold_and_warm_runs_agree(|| Vakint::new().unwrap(), &inputs);
}

#[test]
fn four_loop_cold_and_warm_ordered_cache_runs_are_identical() {
    let inputs = four_loop_benchmark_inputs("FG");
    assert_eq!(inputs.len(), 2);
    assert_cold_and_warm_runs_agree(super::super::test_data::external_engine, &inputs);
}

#[test]
fn warmed_cache_replays_a_batch_and_held_out_inputs_reuse_shared_steps() {
    let directory = tempfile::tempdir().unwrap();
    let inputs: Vec<String> = THREE_LOOP.iter().map(|input| input.to_string()).collect();
    let engine = Vakint::new().unwrap();
    // Warming ignores a disabled mode and saves into the configured directory.
    let mut warming = settings(OrderedCacheSettings::read_write(directory.path()), "exact");
    warming.ordered_cache.mode = OrderedCacheMode::Off;
    warming.ordered_cache.autosave = OrderedCacheAutosave::Never;
    let parsed: Vec<_> = inputs
        .iter()
        .map(|input| vakint_parse!(input).unwrap())
        .collect();
    let report = engine
        .warm_ordered_cache(&warming, parsed.iter().map(Atom::as_view))
        .unwrap();
    assert_eq!(report.integrals, 3);
    assert!(report.inserts > 0);
    let save = report.save.as_ref().unwrap();
    assert_eq!(save.records_written as u64, report.inserts);
    assert_eq!(
        program_directories(directory.path())
            .iter()
            .filter(|name| name.starts_with("3L-unit-mass-vacuum-k6."))
            .count(),
        1
    );
    // The warmed batch replays completely in a fresh engine.
    let fresh = Vakint::new().unwrap();
    let read_only = settings(OrderedCacheSettings::read_only(directory.path()), "exact");
    let warm = evaluate_all(&fresh, &read_only, &inputs);
    let statistics = fresh.ordered_cache_statistics();
    assert_eq!(statistics.misses, 0);
    // Opaque (recomputed) factorization records are counted apart.
    assert_eq!(
        statistics.hits + opaque_hits(&statistics),
        report.inserts + report.hits + opaque_hits(&report.statistics)
    );
    assert_eq!(
        warm,
        evaluate_all(
            &Vakint::new().unwrap(),
            &settings(OrderedCacheSettings::off(), "exact"),
            &inputs
        )
    );
    // A held-out integral of the same family replays the steps it shares.
    let held_out = ["topo(I3L(muvsq,3,2,2,2,1,2))".to_owned()];
    let before = fresh.ordered_cache_statistics();
    let value = evaluate_all(&fresh, &read_only, &held_out);
    let after = fresh.ordered_cache_statistics();
    assert!(after.hits > before.hits && after.misses > before.misses);
    assert_eq!(after.inserts, 0);
    assert_eq!(
        value,
        evaluate_all(
            &Vakint::new().unwrap(),
            &settings(OrderedCacheSettings::off(), "exact"),
            &held_out
        )
    );
}

#[test]
fn a_cache_directory_of_another_program_is_rejected() {
    let directory = tempfile::tempdir().unwrap();
    let engine = Vakint::new().unwrap();
    let read_write = settings(OrderedCacheSettings::read_write(directory.path()), "exact");
    evaluate_all(
        &engine,
        &read_write,
        &[
            "topo(I2L(muvsq,2,1,1))".to_owned(),
            "topo(I3L(muvsq,2,1,1,1,1,1))".to_owned(),
        ],
    );
    let names = program_directories(directory.path());
    let two_loop = names.iter().find(|name| name.starts_with("2L-")).unwrap();
    let three_loop = names.iter().find(|name| name.starts_with("3L-")).unwrap();
    // Replace the 3L manifest by the 2L one: its stamp names another program.
    std::fs::copy(
        directory
            .path()
            .join(two_loop)
            .join("ordered-step-cache.txt"),
        directory
            .path()
            .join(three_loop)
            .join("ordered-step-cache.txt"),
    )
    .unwrap();
    let fresh = Vakint::new().unwrap();
    let read_only = settings(OrderedCacheSettings::read_only(directory.path()), "exact");
    let error = fresh
        .evaluate(
            &read_only,
            vakint_parse!("topo(I3L(muvsq,2,1,1,1,1,1))")
                .unwrap()
                .as_view(),
        )
        .unwrap_err()
        .to_string();
    assert!(
        error.contains("Ordered step cache") && error.contains(three_loop.as_str()),
        "{error}"
    );
    assert!(fresh.load_ordered_cache(&read_only).is_err());
    // The other programs are unaffected.
    let two_loop_value = fresh
        .evaluate(
            &read_only,
            vakint_parse!("topo(I2L(muvsq,2,1,1))").unwrap().as_view(),
        )
        .unwrap();
    assert!(!two_loop_value.is_zero());
    assert_eq!(fresh.ordered_cache_statistics().misses, 0);
}

#[test]
fn ordered_cache_settings_select_the_attached_caches() {
    let directory = tempfile::tempdir().unwrap();
    let input = "topo(I2L(muvsq,2,1,1))".to_owned();
    let engine = Vakint::new().unwrap();
    // Off attaches nothing; read-only without a directory has nothing to read.
    evaluate_all(
        &engine,
        &settings(OrderedCacheSettings::off(), "exact"),
        std::slice::from_ref(&input),
    );
    let mut read_only = settings(OrderedCacheSettings::default(), "exact");
    read_only.ordered_cache.mode = OrderedCacheMode::ReadOnly;
    evaluate_all(&engine, &read_only, std::slice::from_ref(&input));
    assert!(engine.ordered_cache_statistics().programs.is_empty());
    // The default records in memory only.
    let in_memory = settings(OrderedCacheSettings::default(), "exact");
    evaluate_all(&engine, &in_memory, std::slice::from_ref(&input));
    let statistics = engine.ordered_cache_statistics();
    assert_eq!(statistics.programs.len(), 1);
    assert_eq!(statistics.programs[0].program, "2L-unit-mass-vacuum-k3");
    assert!(statistics.programs[0].directory.is_none());
    assert!(engine.save_ordered_cache().unwrap().programs.is_empty());
    // Without autosave nothing is written until an explicit save; clones
    // share the caches.
    let mut manual = settings(OrderedCacheSettings::read_write(directory.path()), "exact");
    manual.ordered_cache.autosave = OrderedCacheAutosave::Never;
    let clone = engine.clone();
    evaluate_all(&clone, &manual, std::slice::from_ref(&input));
    assert!(program_directories(directory.path()).is_empty());
    let saved = engine.save_ordered_cache().unwrap();
    assert_eq!(saved.programs.len(), 1);
    assert!(saved.records_written > 0);
    assert_eq!(program_directories(directory.path()).len(), 1);
    // Loading attaches every loaded program's cache.
    let loaded = engine.load_ordered_cache(&manual).unwrap();
    assert_eq!(
        loaded
            .programs
            .iter()
            .map(|program| program.program.as_str())
            .collect::<Vec<_>>(),
        [
            "1L-unit-mass-vacuum-k1",
            "2L-unit-mass-vacuum-k3",
            "3L-unit-mass-vacuum-k6"
        ]
    );
    clone.clear_ordered_cache();
    assert!(engine.ordered_cache_statistics().programs.is_empty());
    // On-drop autosave writes when the last clone is released.
    let on_drop_directory = tempfile::tempdir().unwrap();
    let mut on_drop = settings(
        OrderedCacheSettings::read_write(on_drop_directory.path()),
        "exact",
    );
    on_drop.ordered_cache.autosave = OrderedCacheAutosave::OnDrop;
    let engine = Vakint::new().unwrap();
    let clone = engine.clone();
    evaluate_all(&engine, &on_drop, std::slice::from_ref(&input));
    drop(engine);
    assert!(program_directories(on_drop_directory.path()).is_empty());
    drop(clone);
    assert_eq!(program_directories(on_drop_directory.path()).len(), 1);
}

#[test]
fn ordered_cache_modes_and_autosave_policies_round_trip() {
    for mode in [
        OrderedCacheMode::Off,
        OrderedCacheMode::ReadOnly,
        OrderedCacheMode::ReadWrite,
    ] {
        assert_eq!(mode.to_string().parse::<OrderedCacheMode>().unwrap(), mode);
    }
    for autosave in [
        OrderedCacheAutosave::Never,
        OrderedCacheAutosave::AfterEvaluation,
        OrderedCacheAutosave::OnDrop,
    ] {
        assert_eq!(
            autosave
                .to_string()
                .parse::<OrderedCacheAutosave>()
                .unwrap(),
            autosave
        );
    }
    assert_eq!(
        "read_only".parse::<OrderedCacheMode>().unwrap(),
        OrderedCacheMode::ReadOnly
    );
    assert!("cached".parse::<OrderedCacheMode>().is_err());
    assert!("sometimes".parse::<OrderedCacheAutosave>().is_err());
    let defaults = OrderedCacheSettings::default();
    assert_eq!(defaults.mode, OrderedCacheMode::ReadWrite);
    assert!(defaults.directory.is_none());
    assert_eq!(VakintSettings::default().ordered_cache, defaults);
}
