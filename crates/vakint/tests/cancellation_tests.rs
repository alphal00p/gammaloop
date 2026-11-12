#![cfg(unix)]
//! A cancelled lower-loop evaluation stops promptly and leaves its engine
//! (and the engine's Ordered step cache) exactly as usable as a fresh one.

#[path = "support/rustred.rs"]
#[allow(dead_code)]
mod rustred_support;

use std::{
    sync::{
        Arc, OnceLock,
        atomic::{AtomicBool, Ordering},
    },
    time::{Duration, Instant},
};

use symbolica::atom::{Atom, AtomCore};
use vakint::{
    Cancellation, LoopNormalizationFactor, OrderedCacheSettings, RustRedMasterMode,
    RustRedMonitorOptions, RustRedProgressSink, Vakint, VakintSettings, vakint_parse,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

/// The notebook gallery's three-loop quadratic ladder coefficient: ten scalar
/// terms whose Ordered reduction runs long enough to be interrupted inside.
fn three_loop_quadratic_ladder() -> Result<Atom> {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/examples/notebooks/data/ladder/3L-quadratic.vakint-input.symbolica"
    );
    let source = vakint_parse!(&std::fs::read_to_string(path)?)?;
    let dimension = vakint_parse!("rustred::d")?;
    let binding = vakint_parse!("4-2*ε")?;
    Ok(source
        .replace(dimension.to_pattern())
        .with(binding.to_pattern()))
}

fn settings() -> VakintSettings {
    VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        use_dot_product_notation: true,
        n_cores: 1,
        // The slowest arithmetic: the most native work after the first census.
        ordered_weight_arithmetic: vakint::OrderedWeightArithmetic::exact(),
        // In-memory read-write cache, never seeded by shipped records.
        ordered_cache: OrderedCacheSettings {
            shipped: false,
            ..Default::default()
        },
        ..rustred_support::rustred_only(
            RustRedMasterMode::Retain,
            "/cancellation-test-must-not-use-form",
        )
    }
}

#[test]
fn cancelled_lower_loop_evaluation_stops_and_keeps_the_engine_reusable() -> Result<()> {
    let input = three_loop_quadratic_ladder()?;
    let settings = settings();
    let engine = Vakint::new()?;
    let sink = RustRedProgressSink::default();
    let options = RustRedMonitorOptions {
        progress_sink: Some(sink.clone()),
        interval: Duration::ZERO,
        ..Default::default()
    };
    let cancellation = Cancellation::default();
    let cancelled_at = OnceLock::new();
    let finished = Arc::new(AtomicBool::new(false));
    let (result, returned_at) = std::thread::scope(|scope| {
        let watcher = {
            let (sink, cancellation, finished, cancelled_at) =
                (&sink, &cancellation, &finished, &cancelled_at);
            scope.spawn(move || {
                let mut seen = 0;
                while !finished.load(Ordering::Acquire) {
                    if let Some((sequence, model)) = sink.latest_since(seen) {
                        seen = sequence;
                        // Inside the native Ordered reduction.
                        if !model["ordered"].is_null() || model["stage"]["kind"] == "native" {
                            cancelled_at.set(Instant::now()).unwrap();
                            cancellation.cancel();
                            return Some(model);
                        }
                    }
                    std::thread::sleep(Duration::from_millis(2));
                }
                None
            })
        };
        let result =
            cancellation.run(|| engine.evaluate_with_monitor(&settings, input.as_view(), &options));
        let returned_at = Instant::now();
        finished.store(true, Ordering::Release);
        let trigger = watcher.join().unwrap();
        assert!(
            trigger.is_some(),
            "the evaluation finished before its native reduction was observed"
        );
        (result, returned_at)
    });
    assert!(cancellation.is_cancelled());
    let error = result.expect_err("a cancelled request must not return a result");
    let latency = returned_at.duration_since(*cancelled_at.get().unwrap());
    assert!(latency < Duration::from_secs(2), "{latency:?}: {error}");
    // The monitor's final model records the failed request.
    let (_, last) = sink.latest_since(0).unwrap();
    assert_eq!(last["status"], "Failed");

    // The same engine (and its partially populated step cache) without
    // cancellation agrees with a fresh engine.
    let again = engine.evaluate(&settings, input.as_view())?;
    let fresh = Vakint::new()?.evaluate(&settings, input.as_view())?;
    assert_eq!(again, fresh);
    assert!(!fresh.is_zero());
    // A pending cancellation fails the next request before it starts work.
    let started = Instant::now();
    assert!(
        cancellation
            .run(|| engine.evaluate(&settings, input.as_view()))
            .is_err()
    );
    assert!(started.elapsed() < Duration::from_secs(1));
    Ok(())
}
