//! `GAMMALOOP_DELTA_KEY`: the routing-free graph key of R* child Δs
//! (delta-store lane, step 1: key, transport, in-process cache, verify).
//!
//! `legacy` (the default): today's key, bit for bit. `graph`: the key of
//! [`super::counterterm_identity::CountertermIdentity::for_graph_key`], which
//! sees the child graph with its Feynman-rule content only (no parent routing,
//! numbering, boson orientation, rule signs, closed quark-loop direction, and
//! open quark-line direction with the transposition transport).
//!
//! `GAMMALOOP_DELTA_VERIFY=off|on|<N>`: `on` recomputes every graph-key hit
//! (in the process cache) and compares it with the transported stored value;
//! `N >= 2` recomputes a salted 1-in-N sample of the hits in a private cold
//! cache. A mismatch fails the lookup (`DELTA_VERIFY_MISMATCH`).
//!
//! `GAMMALOOP_DELTA_KEY_DUMP=stderr`: one `DELTA_KEY` line per graph-key lookup.
//! `GAMMALOOP_DELTA_STATS=1`: the counters below also in legacy mode
//! (`run.json` `delta_cache_stats.delta_key`).

use std::{
    collections::BTreeMap,
    sync::{
        LazyLock, Mutex,
        atomic::{AtomicU64, Ordering},
    },
};

use eyre::{Result, eyre};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DeltaKeyMode {
    Legacy,
    Graph,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DeltaVerify {
    Off,
    On,
    Sample(u64),
}

fn parse_mode() -> Result<DeltaKeyMode> {
    match std::env::var("GAMMALOOP_DELTA_KEY").as_deref().map(str::trim) {
        Err(_) | Ok("") | Ok("legacy") => Ok(DeltaKeyMode::Legacy),
        Ok("graph") => Ok(DeltaKeyMode::Graph),
        Ok(other) => Err(eyre!(
            "GAMMALOOP_DELTA_KEY={other}: expected legacy or graph"
        )),
    }
}

fn parse_verify() -> Result<DeltaVerify> {
    match std::env::var("GAMMALOOP_DELTA_VERIFY")
        .as_deref()
        .map(str::trim)
    {
        Err(_) | Ok("") | Ok("off") | Ok("0") => Ok(DeltaVerify::Off),
        Ok("on") | Ok("1") => Ok(DeltaVerify::On),
        Ok(other) => match other.parse::<u64>() {
            Ok(n) if n >= 2 => Ok(DeltaVerify::Sample(n)),
            _ => Err(eyre!(
                "GAMMALOOP_DELTA_VERIFY={other}: expected off, on or an integer N >= 2"
            )),
        },
    }
}

/// Startup check of the switches (beta calls it before any work).
pub fn check_env() -> Result<()> {
    parse_mode()?;
    parse_verify()?;
    match std::env::var("GAMMALOOP_DELTA_KEY_DUMP").as_deref() {
        Err(_) | Ok("") | Ok("stderr") | Ok("0") | Ok("off") => {}
        Ok(other) => return Err(eyre!("GAMMALOOP_DELTA_KEY_DUMP={other}: expected stderr")),
    }
    Ok(())
}

pub fn mode() -> DeltaKeyMode {
    static MODE: LazyLock<DeltaKeyMode> =
        LazyLock::new(|| parse_mode().unwrap_or_else(|error| panic!("{error}")));
    *MODE
}

pub fn verify() -> DeltaVerify {
    static VERIFY: LazyLock<DeltaVerify> =
        LazyLock::new(|| parse_verify().unwrap_or_else(|error| panic!("{error}")));
    *VERIFY
}

pub(crate) fn dump() -> bool {
    static DUMP: LazyLock<bool> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_DELTA_KEY_DUMP").as_deref() == Ok("stderr")
    });
    *DUMP
}

/// Whether the counters are kept: graph mode, or `GAMMALOOP_DELTA_STATS=1`.
pub(crate) fn stats_enabled() -> bool {
    static ON: LazyLock<bool> = LazyLock::new(|| {
        mode() == DeltaKeyMode::Graph
            || matches!(
                std::env::var("GAMMALOOP_DELTA_STATS").as_deref(),
                Ok("1") | Ok("on")
            )
    });
    *ON
}

/// 64 random-enough bits per process for the sampled verify (time and pid).
pub(crate) fn run_salt() -> u64 {
    static SALT: LazyLock<u64> = LazyLock::new(|| {
        use std::hash::{BuildHasher, Hasher};
        let mut hasher = std::collections::hash_map::RandomState::new().build_hasher();
        hasher.write_u128(
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_nanos())
                .unwrap_or_default(),
        );
        hasher.write_u32(std::process::id());
        hasher.finish()
    });
    *SALT
}

/// Whether the hit of `key_digest` is in the salted 1-in-`n` sample.
pub(crate) fn sampled(key_digest: u64, n: u64) -> bool {
    use std::hash::{Hash, Hasher};
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    (run_salt(), key_digest).hash(&mut hasher);
    hasher.finish() < u64::MAX / n
}

macro_rules! counters {
    ($($name:ident),* $(,)?) => {
        #[allow(non_snake_case)]
        pub(crate) struct Counters { $(pub(crate) $name: AtomicU64,)* }
        pub(crate) static C: Counters = Counters { $($name: AtomicU64::new(0),)* };
        fn counters_json() -> serde_json::Map<String, serde_json::Value> {
            let mut map = serde_json::Map::new();
            $(map.insert(stringify!($name).to_owned(), C.$name.load(Ordering::Relaxed).into());)*
            map
        }
    };
}

counters!(
    lookups,
    graph_lookups,
    graph_hits,
    graph_misses,
    legacy_lookups,
    legacy_hits,
    legacy_misses,
    key_ns,
    key_fallback_nonincident,
    key_fallback_coef,
    key_fallback_error,
    scoped_flag_lookup,
    lines_open_quotiented,
    lines_closed_quotiented,
    line_not_quotiented,
    transport_unsupported,
    fill_ns_outer,
    verify_checked,
    verify_identical,
    verify_equal_nf,
    verify_mismatch,
    verify_errors,
    verify_ns,
    nf_conservation_used,
    nf_conservation_skipped,
    hit_legacy_would_miss,
);

pub(crate) fn add(counter: &AtomicU64, n: u64) {
    counter.fetch_add(n, Ordering::Relaxed);
}

/// The symmetry features a hit shows against the occurrence that filled its
/// entry (design §3.4); every feature is counted independently.
pub(crate) const FEATURES: [&str; 9] = [
    "routing",
    "relabel",
    "leg_permutation",
    "boson_orientation",
    "closed_loop_reversed",
    "open_line_reversed",
    "sigma_negative",
    "crown_sign",
    "legacy_key_differs",
];

#[derive(Default)]
struct FeatureCounts {
    hits: [u64; FEATURES.len()],
    verified: [u64; FEATURES.len()],
    /// hits that show no feature at all (an exact repeat)
    plain_hits: u64,
    plain_verified: u64,
}

/// Per child loop order: lookups, hits, fills and the exclusive fill seconds.
#[derive(Default, Clone, Copy)]
struct LoopCounts {
    lookups: u64,
    hits: u64,
    fills: u64,
    fill_ns_exclusive: u64,
    fill_ns_inclusive: u64,
}

static FEATURE_COUNTS: LazyLock<Mutex<FeatureCounts>> =
    LazyLock::new(|| Mutex::new(FeatureCounts::default()));
static LOOP_COUNTS: LazyLock<Mutex<BTreeMap<usize, LoopCounts>>> =
    LazyLock::new(|| Mutex::new(BTreeMap::new()));

pub(crate) fn count_features(shown: &[bool; FEATURES.len()], verified: bool) {
    let mut counts = FEATURE_COUNTS.lock().expect("delta key counters poisoned");
    let any = shown.iter().any(|shown| *shown);
    if !any {
        counts.plain_hits += 1;
        if verified {
            counts.plain_verified += 1;
        }
    }
    for (k, shown) in shown.iter().enumerate() {
        if *shown {
            counts.hits[k] += 1;
            if verified {
                counts.verified[k] += 1;
            }
        }
    }
}

pub(crate) fn count_loop(loops: usize, hit: bool) {
    let mut counts = LOOP_COUNTS.lock().expect("delta key counters poisoned");
    let entry = counts.entry(loops).or_default();
    entry.lookups += 1;
    if hit {
        entry.hits += 1;
    }
}

thread_local! {
    /// The inclusive nanoseconds of the fills nested in the fill running on
    /// this thread (one slot per open fill).
    static FILL_STACK: std::cell::RefCell<Vec<u64>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// A running fill: its inclusive time, minus that of its nested fills, is
/// charged to its loop order on drop.
pub(crate) struct FillTimer {
    loops: usize,
    started: std::time::Instant,
}

impl FillTimer {
    pub(crate) fn start(loops: usize) -> Self {
        FILL_STACK.with(|stack| stack.borrow_mut().push(0));
        Self {
            loops,
            started: std::time::Instant::now(),
        }
    }
}

impl Drop for FillTimer {
    fn drop(&mut self) {
        let inclusive = self.started.elapsed().as_nanos() as u64;
        let nested = FILL_STACK.with(|stack| {
            let mut stack = stack.borrow_mut();
            let nested = stack.pop().unwrap_or(0);
            if let Some(parent) = stack.last_mut() {
                *parent += inclusive;
            }
            if stack.is_empty() {
                add(&C.fill_ns_outer, inclusive);
            }
            nested
        });
        let mut counts = LOOP_COUNTS.lock().expect("delta key counters poisoned");
        let entry = counts.entry(self.loops).or_default();
        entry.fills += 1;
        entry.fill_ns_inclusive += inclusive;
        entry.fill_ns_exclusive += inclusive.saturating_sub(nested);
    }
}

/// The counters as JSON, or `None` when they are not kept.
pub fn stats_json() -> Option<serde_json::Value> {
    if !stats_enabled() {
        return None;
    }
    let mut map = counters_json();
    map.insert(
        "mode".into(),
        format!("{:?}", mode()).to_lowercase().into(),
    );
    map.insert("verify".into(), format!("{:?}", verify()).into());
    let features = FEATURE_COUNTS.lock().expect("delta key counters poisoned");
    let mut by_feature = serde_json::Map::new();
    for (k, name) in FEATURES.iter().enumerate() {
        by_feature.insert(
            (*name).into(),
            serde_json::json!({"hits": features.hits[k], "verified": features.verified[k]}),
        );
    }
    by_feature.insert(
        "none".into(),
        serde_json::json!({"hits": features.plain_hits, "verified": features.plain_verified}),
    );
    map.insert("features".into(), by_feature.into());
    let loops = LOOP_COUNTS.lock().expect("delta key counters poisoned");
    let mut by_loop = serde_json::Map::new();
    for (loops, counts) in loops.iter() {
        by_loop.insert(
            loops.to_string(),
            serde_json::json!({
                "lookups": counts.lookups, "hits": counts.hits, "fills": counts.fills,
                "fill_seconds_exclusive": counts.fill_ns_exclusive as f64 * 1e-9,
                "fill_seconds_inclusive": counts.fill_ns_inclusive as f64 * 1e-9,
            }),
        );
    }
    map.insert("by_child_loops".into(), by_loop.into());
    Some(map.into())
}
