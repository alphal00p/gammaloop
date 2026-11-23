//! `GAMMALOOP_DELTA_STORE`: the persistent store of R* child Δs under the
//! graph key (delta-store lane, step 2).
//!
//! gammalooprs owns the key bytes, the context, the Symbolica export and
//! import, the transport and the checks; the store itself is a
//! [`DeltaBackingStore`] that sees bytes only (beta: SQLite, `src/delta_store.rs`).
//! The process installs one with [`set_backing_store`]; without one (the
//! default, mode `off`) every lookup is step 1's, bit for bit.
//!
//! Exactness never depends on the store: a row that fails any acceptance
//! check is a miss, and in mode `verify` every store hit is recomputed in a
//! private cold cache and compared (`DELTA_STORE_VERIFY_MISMATCH` fails the
//! lookup). A value another writer stored first that differs from this
//! process's value is compared the same way; a difference is a dispute and
//! fails the lookup.

use std::{
    cell::RefCell,
    collections::{BTreeMap, BTreeSet, HashMap},
    sync::{
        Arc, LazyLock, Mutex, OnceLock,
        atomic::{AtomicU64, Ordering},
    },
};

use eyre::{Result, eyre};
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore};

/// The store mode a process runs under.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StoreMode {
    Off,
    /// read only: hits are served, nothing is written
    Ro,
    /// read and write
    Rw,
    /// rw, and every store hit is recomputed and compared
    Verify,
}

impl StoreMode {
    pub fn parse(text: &str) -> Result<Self> {
        match text.trim() {
            "" | "off" | "0" => Ok(Self::Off),
            "ro" => Ok(Self::Ro),
            "rw" => Ok(Self::Rw),
            "verify" => Ok(Self::Verify),
            other => Err(eyre!(
                "GAMMALOOP_DELTA_STORE={other}: expected off, ro, rw or verify"
            )),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::Ro => "ro",
            Self::Rw => "rw",
            Self::Verify => "verify",
        }
    }
}

/// The key format and Δ semantics of this build (design §1.4, §1.5).
pub const KEY_FORMAT: u16 = 1;
pub const DELTA_SEMANTICS: u16 = 1;
pub const KIND_RSTAR_CHILD: &str = "rstar-child-delta";

/// A stored row as the backend returns it.
#[derive(Clone, Debug)]
pub struct StoredRow {
    pub tier: u8,
    pub value_format: String,
    pub value_export: Vec<u8>,
    pub blob_sha: [u8; 32],
    pub value_sha: [u8; 32],
    pub note: Option<String>,
    pub writer_transport: String,
    pub writer_build: String,
    pub writer_run_id: String,
    pub dep_builds: Vec<String>,
    /// the writer's fill seconds (inclusive of its nested fills)
    pub compute_s: f64,
}

pub enum StoreGet {
    Hit(StoredRow),
    Miss,
    /// The store could not be read (busy, corrupt, ...): a miss with a reason.
    Error(String),
}

/// What a write did. The write is synchronous (one small transaction on the
/// computing thread), so a row another process stored first comes back to
/// the thread that holds the occurrence and can compare it exactly.
pub enum PutResult {
    Inserted,
    /// Another writer stored the key first: its row.
    Existing(StoredRow),
    /// Not written (read-only, full, busy, an error): the reason.
    NotWritten(String),
}

/// One value for the store, from the computing thread.
#[derive(Clone, Debug)]
pub struct DeltaRecord {
    pub key_sha: [u8; 32],
    pub context_sha: [u8; 32],
    pub context_json: Arc<str>,
    pub graph_sha: [u8; 32],
    pub model_fp: [u8; 32],
    pub child_loops: usize,
    pub is_zero: bool,
    pub value_format: String,
    pub value_export: Vec<u8>,
    pub blob_sha: [u8; 32],
    pub value_sha: [u8; 32],
    pub note: Option<String>,
    pub writer_transport: String,
    pub writer_graph: String,
    pub writer_support: String,
    pub compute_s: f64,
    pub deps: Vec<[u8; 32]>,
    pub dep_builds: Vec<String>,
    pub key_bytes: Option<Arc<[u8]>>,
    pub text: Option<String>,
}

/// What the computing thread asks the backend to record besides inserts.
pub enum StoreEvent {
    /// Another writer stored the key first; this process's value agrees
    /// (`identical`: equal value_sha; otherwise equal after the normal form).
    Confirm { key_sha: [u8; 32], tier: u8, identical: bool },
    /// A disagreement (`kind`: conflict | verify): the stored row is
    /// quarantined when `quarantine` is set.
    Dispute {
        kind: &'static str,
        key_sha: [u8; 32],
        tier: u8,
        stored_value_sha: [u8; 32],
        record: Box<DeltaRecord>,
        quarantine: bool,
    },
}

/// The persistent store behind the Δ cache.
pub trait DeltaBackingStore: Send + Sync {
    fn get(&self, key_sha: &[u8; 32]) -> StoreGet;
    fn put(&self, record: &DeltaRecord) -> PutResult;
    fn event(&self, event: StoreEvent);
    /// Whether this process writes (rw / verify and not degraded or full).
    fn writable(&self) -> bool;
    fn stats(&self) -> serde_json::Value;
}

struct Installed {
    mode: StoreMode,
    store: Box<dyn DeltaBackingStore>,
}

static INSTALLED: OnceLock<Installed> = OnceLock::new();

/// Install the process's store (once, before any lookup).
pub fn set_backing_store(mode: StoreMode, store: Box<dyn DeltaBackingStore>) -> Result<()> {
    if mode == StoreMode::Off {
        return Ok(());
    }
    INSTALLED
        .set(Installed { mode, store })
        .map_err(|_| eyre!("the Δ store is installed twice"))
}

/// The installed store, if any.
pub(crate) fn installed() -> Option<(StoreMode, &'static dyn DeltaBackingStore)> {
    INSTALLED
        .get()
        .map(|installed| (installed.mode, installed.store.as_ref()))
}

pub fn mode() -> StoreMode {
    INSTALLED.get().map_or(StoreMode::Off, |installed| installed.mode)
}

/// Children with fewer loops are not persisted (`GAMMALOOP_DELTA_STORE_MIN_LOOPS`, default 2).
pub(crate) fn min_loops() -> usize {
    static MIN: LazyLock<usize> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_DELTA_STORE_MIN_LOOPS")
            .ok()
            .and_then(|value| value.trim().parse().ok())
            .unwrap_or(2)
    });
    *MIN
}

/// `GAMMALOOP_DELTA_STORE_AUDIT=N` (default 256; 0 off): a store hit in the
/// salted 1-in-N sample is recomputed in a private cold cache.
pub(crate) fn audit_every() -> u64 {
    static N: LazyLock<u64> = LazyLock::new(|| {
        std::env::var("GAMMALOOP_DELTA_STORE_AUDIT")
            .ok()
            .and_then(|value| value.trim().parse().ok())
            .unwrap_or(256)
    });
    *N
}

fn flag(name: &str) -> bool {
    matches!(std::env::var(name).as_deref(), Ok("1") | Ok("on"))
}

pub(crate) fn store_keys() -> bool {
    static ON: LazyLock<bool> = LazyLock::new(|| flag("GAMMALOOP_DELTA_STORE_KEYS"));
    *ON
}

pub(crate) fn store_text() -> bool {
    static ON: LazyLock<bool> = LazyLock::new(|| flag("GAMMALOOP_DELTA_STORE_TEXT"));
    *ON
}

/// The Symbolica export format of the values (design §4.4 `value_format`).
pub fn value_format() -> &'static str {
    "symbolica-export/1/394f6e16d97823ac3920e95a5537741a26622dcd"
}

// ---------------------------------------------------------------------------
// Context (design §1.5)
// ---------------------------------------------------------------------------

static CONTEXT_PARTS: LazyLock<Mutex<BTreeMap<String, serde_json::Value>>> =
    LazyLock::new(|| Mutex::new(BTreeMap::new()));
static CONTEXT_FROZEN: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);

/// A part of the context (`data` of design §1.5: master data, the effective
/// vakint and UV settings of the engine). The same part set twice must be
/// equal; a part after the first key is an error.
pub fn set_context_part(name: &str, value: serde_json::Value) -> Result<()> {
    let mut parts = CONTEXT_PARTS.lock().expect("delta store context poisoned");
    if let Some(previous) = parts.get(name) {
        if *previous == value {
            return Ok(());
        }
        return Err(eyre!(
            "the Δ store context part {name} is set twice with different values"
        ));
    }
    if CONTEXT_FROZEN.load(Ordering::SeqCst) {
        return Err(eyre!(
            "the Δ store context part {name} is set after the first key"
        ));
    }
    parts.insert(name.to_owned(), value);
    Ok(())
}

/// The value-path environment (design §1.5 `env`).
pub fn context_env() -> Vec<String> {
    let mut lines = std::env::vars()
        .filter(|(name, _)| {
            ["GAMMALOOP_", "IDENSO_", "VAKINT_", "BETA_P_", "RUSTRED_"]
                .iter()
                .any(|prefix| name.starts_with(prefix))
                && name != "SYMBOLICA_LICENSE"
                && !name.starts_with("GAMMALOOP_DELTA_")
                && !name.ends_with("_LOG")
                && !name.ends_with("_GIB")
                && !name.ends_with("_BYTES")
                && name != "BETA_P_LIBRARY"
        })
        .map(|(name, value)| format!("{name}={value}"))
        .collect::<Vec<_>>();
    lines.sort();
    lines
}

/// The names of the symbols the stored values are written in (design §1.5
/// `vocab`; a rename changes the key).
fn vocab() -> Vec<String> {
    use crate::utils::GS;
    let mut names = vec![
        GS.dim_epsilon.get_name().to_owned(),
        GS.mu_r_sq.get_name().to_owned(),
        GS.emr_mom.get_name().to_owned(),
        GS.dim.get_name().to_owned(),
        "uv_delta_cache_private".to_owned(),
        "uv_delta_cache_port".to_owned(),
        "uv_delta_cache_momentum".to_owned(),
        "uv_delta_cache_hedge_momentum".to_owned(),
    ];
    names.sort();
    names
}

pub(crate) struct Context {
    pub sha: [u8; 32],
    pub json: Arc<str>,
}

static CONTEXTS: LazyLock<Mutex<HashMap<String, Arc<Context>>>> =
    LazyLock::new(|| Mutex::new(HashMap::new()));

/// The context of a lookup with today's settings string (memoised).
pub(crate) fn context(settings_key: &str) -> Result<Arc<Context>> {
    if let Some(context) = CONTEXTS
        .lock()
        .expect("delta store context poisoned")
        .get(settings_key)
    {
        return Ok(context.clone());
    }
    CONTEXT_FROZEN.store(true, Ordering::SeqCst);
    let parts = CONTEXT_PARTS
        .lock()
        .expect("delta store context poisoned")
        .clone();
    if parts.is_empty() {
        return Err(eyre!(
            "the Δ store context has no data part (set_context_part was never called)"
        ));
    }
    // serde_json's Map is a BTreeMap here: keys sorted, no whitespace.
    let json = serde_json::to_string(&serde_json::json!({
        "data": parts,
        "env": context_env(),
        "key_format": KEY_FORMAT,
        "semantics": DELTA_SEMANTICS,
        "settings": settings_key,
        "value_format": value_format(),
        "vocab": vocab(),
    }))?;
    let context = Arc::new(Context {
        sha: Sha256::digest(json.as_bytes()).into(),
        json: json.into(),
    });
    CONTEXTS
        .lock()
        .expect("delta store context poisoned")
        .insert(settings_key.to_owned(), context.clone());
    Ok(context)
}

/// GLDK/1 key bytes and their hashes: `(key_sha, graph_sha)`.
pub(crate) fn key_sha(
    model_fp: &[u8; 32],
    context_sha: &[u8; 32],
    graph_bytes: &[u8],
) -> ([u8; 32], [u8; 32], Vec<u8>) {
    let mut bytes = Vec::with_capacity(graph_bytes.len() + 128);
    bytes.extend_from_slice(b"GLDK");
    bytes.extend_from_slice(&KEY_FORMAT.to_le_bytes());
    bytes.extend_from_slice(&DELTA_SEMANTICS.to_le_bytes());
    bytes.extend_from_slice(&(KIND_RSTAR_CHILD.len() as u32).to_le_bytes());
    bytes.extend_from_slice(KIND_RSTAR_CHILD.as_bytes());
    bytes.extend_from_slice(model_fp);
    bytes.extend_from_slice(context_sha);
    bytes.extend_from_slice(graph_bytes);
    let graph_sha: [u8; 32] = Sha256::digest(graph_bytes).into();
    let key_sha: [u8; 32] = Sha256::digest(&bytes).into();
    (key_sha, graph_sha, bytes)
}

pub fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

// ---------------------------------------------------------------------------
// Values: export, import, acceptance (design §4.5)
// ---------------------------------------------------------------------------

pub(crate) struct Exported {
    pub export: Vec<u8>,
    pub blob_sha: [u8; 32],
    pub value_sha: [u8; 32],
    pub text: String,
}

/// The export of a canonical value, or why it is not stored.
pub(crate) fn export(value: &Atom) -> std::result::Result<Exported, &'static str> {
    let text = value.to_canonical_string();
    if text.contains("UNKNOWN") {
        return Err("unknown");
    }
    let mut export = Vec::new();
    value
        .as_view()
        .export(&mut export)
        .map_err(|_| "export_error")?;
    Ok(Exported {
        blob_sha: Sha256::digest(&export).into(),
        value_sha: Sha256::digest(text.as_bytes()).into(),
        export,
        text,
    })
}

/// The canonical value of a row, or the reason it is refused (a miss).
pub(crate) fn import(row: &StoredRow) -> std::result::Result<Atom, &'static str> {
    if row.value_format != value_format() {
        return Err("value_format");
    }
    let blob_sha: [u8; 32] = Sha256::digest(&row.value_export).into();
    if blob_sha != row.blob_sha {
        return Err("blob_sha");
    }
    let value = Atom::import(&mut &row.value_export[..], None).map_err(|_| "import")?;
    let text = value.to_canonical_string();
    let value_sha: [u8; 32] = Sha256::digest(text.as_bytes()).into();
    if value_sha != row.value_sha {
        return Err("value_sha");
    }
    Ok(value)
}

// ---------------------------------------------------------------------------
// Dependencies (design §3.2): the store keys a value was derived from
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Default)]
pub(crate) struct Deps {
    pub keys: BTreeSet<[u8; 32]>,
    pub builds: BTreeSet<String>,
}

impl Deps {
    pub(crate) fn merge(&mut self, other: &Deps) {
        self.keys.extend(other.keys.iter().copied());
        self.builds.extend(other.builds.iter().cloned());
    }
}

thread_local! {
    static COLLECTORS: RefCell<Vec<Deps>> = const { RefCell::new(Vec::new()) };
}

/// A fill's dependency collector: nested lookups on this thread add to it.
pub(crate) struct Collector {
    depth: usize,
}

impl Collector {
    pub(crate) fn push() -> Self {
        let depth = COLLECTORS.with(|stack| {
            let mut stack = stack.borrow_mut();
            stack.push(Deps::default());
            stack.len()
        });
        Self { depth }
    }

    pub(crate) fn finish(self) -> Deps {
        let deps = COLLECTORS.with(|stack| {
            let mut stack = stack.borrow_mut();
            debug_assert_eq!(stack.len(), self.depth, "delta store collectors out of order");
            stack.pop().unwrap_or_default()
        });
        std::mem::forget(self);
        deps
    }
}

impl Drop for Collector {
    fn drop(&mut self) {
        COLLECTORS.with(|stack| {
            let mut stack = stack.borrow_mut();
            if stack.len() == self.depth {
                stack.pop();
            }
        });
    }
}

/// A finished lookup adds its entry's dependencies to the enclosing fill.
pub(crate) fn add_to_enclosing(deps: &Deps) {
    if deps.keys.is_empty() && deps.builds.is_empty() {
        return;
    }
    COLLECTORS.with(|stack| {
        if let Some(top) = stack.borrow_mut().last_mut() {
            top.merge(deps);
        }
    });
}

// ---------------------------------------------------------------------------
// Counters and the hit list (design §7)
// ---------------------------------------------------------------------------

macro_rules! store_counters {
    ($($name:ident),* $(,)?) => {
        pub(crate) struct Counters { $(pub(crate) $name: AtomicU64,)* }
        pub(crate) static S: Counters = Counters { $($name: AtomicU64::new(0),)* };
        fn counters_json() -> serde_json::Map<String, serde_json::Value> {
            let mut map = serde_json::Map::new();
            $(map.insert(stringify!($name).to_owned(), S.$name.load(Ordering::Relaxed).into());)*
            map
        }
    };
}

store_counters!(
    store_lookups,
    store_hits,
    store_misses,
    store_below_min_loops,
    store_not_portable,
    store_context_error,
    read_error,
    reject_value_format,
    reject_blob_sha,
    reject_import,
    reject_value_sha,
    reject_note,
    reject_restore,
    stored,
    not_stored_unknown,
    not_stored_export,
    not_stored_ro,
    not_written,
    existing_identical,
    existing_agree,
    disputes,
    verify_checked,
    verify_identical,
    verify_equal_nf,
    verify_mismatch,
    audit_checked,
    key_ns,
    read_ns,
    import_ns,
    export_ns,
    verify_ns,
    hit_compute_s_saved_ms,
);

pub(crate) fn add(counter: &AtomicU64, n: u64) {
    counter.fetch_add(n, Ordering::Relaxed);
}

#[derive(Clone, Debug)]
struct HitLine {
    tier: u8,
    value_sha: [u8; 32],
    writer_build: String,
    writer_run_id: String,
    uses: u64,
    audited: bool,
}

static HITS: LazyLock<Mutex<BTreeMap<[u8; 32], HitLine>>> =
    LazyLock::new(|| Mutex::new(BTreeMap::new()));

pub(crate) fn record_hit(key_sha: &[u8; 32], row: &StoredRow, audited: bool) {
    let mut hits = HITS.lock().expect("delta store hits poisoned");
    let line = hits.entry(*key_sha).or_insert_with(|| HitLine {
        tier: row.tier,
        value_sha: row.value_sha,
        writer_build: row.writer_build.clone(),
        writer_run_id: row.writer_run_id.clone(),
        uses: 0,
        audited: false,
    });
    line.uses += 1;
    line.audited |= audited;
}

/// The hit list (`delta_store_hits.tsv`): one sorted line per store row used.
pub fn hits_tsv() -> String {
    let hits = HITS.lock().expect("delta store hits poisoned");
    let mut text = String::from("key_sha\ttier\tvalue_sha\twriter_build\twriter_run_id\tuses\taudited\n");
    for (key, line) in hits.iter() {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            hex(key),
            line.tier,
            hex(&line.value_sha),
            line.writer_build,
            line.writer_run_id,
            line.uses,
            line.audited as u8
        ));
    }
    text
}

/// Per child loop order: store hits, store misses, values stored.
static BY_LOOPS: LazyLock<Mutex<BTreeMap<usize, (u64, u64, u64)>>> =
    LazyLock::new(|| Mutex::new(BTreeMap::new()));

pub(crate) fn count_loops(loops: usize, hit: bool, stored: bool) {
    let mut by = BY_LOOPS.lock().expect("delta store counters poisoned");
    let entry = by.entry(loops).or_default();
    if hit {
        entry.0 += 1;
    } else if stored {
        entry.2 += 1;
    } else {
        entry.1 += 1;
    }
}

/// The process's store block for run.json (`None` with the store off).
pub fn stats_json() -> Option<serde_json::Value> {
    let (mode, store) = installed()?;
    let mut map = counters_json();
    map.insert("mode".into(), mode.as_str().into());
    let by = BY_LOOPS.lock().expect("delta store counters poisoned");
    let mut by_loops = serde_json::Map::new();
    for (loops, (hits, misses, stored)) in by.iter() {
        by_loops.insert(
            loops.to_string(),
            serde_json::json!({"store_hits": hits, "store_misses": misses, "stored": stored}),
        );
    }
    map.insert("by_child_loops".into(), by_loops.into());
    let contexts = CONTEXTS.lock().expect("delta store context poisoned");
    let mut shas = contexts
        .values()
        .map(|context| hex(&context.sha))
        .collect::<Vec<_>>();
    shas.sort();
    shas.dedup();
    map.insert("context_shas".into(), shas.into());
    map.insert("backend".into(), store.stats());
    Some(map.into())
}
