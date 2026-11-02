//! A saved Ordered step cache is replayed by a fresh process.

#[path = "support/rustred.rs"]
#[allow(dead_code)]
mod rustred_support;

use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore};
use vakint::{
    LoopNormalizationFactor, OrderedCacheSettings, RustRedMasterMode, Vakint, VakintSettings,
    vakint_parse,
};

/// Set by the parent test for its child process.
const CHILD_DIRECTORY: &str = "VAKINT_ORDERED_CACHE_TEST_DIRECTORY";
const CHILD_DIGEST: &str = "VAKINT_ORDERED_CACHE_TEST_DIGEST";

const INTEGRALS: [&str; 3] = [
    "topo(I2L(muvsq,3,2,1))",
    "topo(I3L(muvsq,2,2,2,1,1,1))",
    "dot(k(1),k(2))*dot(k(2),k(3))*topo(I3L(muvsq,2,1,1,1,1,1))",
];

fn settings(cache: OrderedCacheSettings) -> VakintSettings {
    VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        use_dot_product_notation: true,
        ordered_cache: cache,
        ..rustred_support::rustred_only(
            RustRedMasterMode::Retain,
            "/ordered-step-cache-test-must-not-call-form",
        )
    }
}

/// SHA-256 of the canonical retained results, in input order.
fn evaluate_digest(engine: &Vakint, settings: &VakintSettings) -> String {
    let mut digest = Sha256::new();
    for integral in INTEGRALS {
        let value: Atom = engine
            .evaluate(settings, vakint_parse!(integral).unwrap().as_view())
            .unwrap();
        digest.update(value.to_canonical_string().as_bytes());
        digest.update(b"\n");
    }
    format!("{:x}", digest.finalize())
}

#[test]
fn persisted_ordered_cache_is_replayed_by_a_fresh_process() {
    let directory = tempfile::tempdir().unwrap();
    let uncached = evaluate_digest(
        &Vakint::new().unwrap(),
        &settings(OrderedCacheSettings::off()),
    );
    let engine = Vakint::new().unwrap();
    let cold = evaluate_digest(
        &engine,
        &settings(OrderedCacheSettings::read_write(directory.path())),
    );
    assert_eq!(cold, uncached);
    let statistics = engine.ordered_cache_statistics();
    assert!(statistics.inserts > 0);
    drop(engine);
    let output = std::process::Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "replay_saved_ordered_cache_in_a_child_process",
            "--ignored",
            "--test-threads=1",
        ])
        .env(CHILD_DIRECTORY, directory.path())
        .env(CHILD_DIGEST, &uncached)
        .output()
        .unwrap();
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        output.status.success() && stdout.contains("1 passed"),
        "child failed:\n{stdout}\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
#[ignore = "run by persisted_ordered_cache_is_replayed_by_a_fresh_process in a child process"]
fn replay_saved_ordered_cache_in_a_child_process() {
    let directory = std::env::var_os(CHILD_DIRECTORY).expect("run by the parent test");
    let expected = std::env::var(CHILD_DIGEST).expect("run by the parent test");
    let engine = Vakint::new().unwrap();
    let read_only = settings(OrderedCacheSettings::read_only(directory));
    assert_eq!(evaluate_digest(&engine, &read_only), expected);
    let statistics = engine.ordered_cache_statistics();
    assert!(statistics.hits > 0, "{statistics:?}");
    assert_eq!(statistics.misses, 0, "{statistics:?}");
    assert_eq!(statistics.inserts, 0);
    assert!(
        statistics
            .programs
            .iter()
            .all(|program| program.native.loaded_shards > 0 && program.native.load_seconds >= 0.0)
    );
}
