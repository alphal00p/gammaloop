//! Reuse complete local UV counterterms, including their subdivergences.

use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    sync::{
        Arc, Condvar, Mutex, OnceLock,
        atomic::{AtomicU64, Ordering},
    },
    time::Instant,
};

use eyre::{Result, ensure, eyre};
use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};
use symbolica::atom::Atom;
use vakint::{Vakint, stage_timing::StageSpan};

use crate::{
    debug_tags,
    graph::Graph,
    model::Model,
    uv::{
        ApproximationType, Spinney, UVgenerationSettings, UltravioletGraph,
        approx::{
            UVCtx,
            integrated::Integrated,
            local_4d::{self, Full4dCts},
        },
        counterterm_identity::{
            CountertermIdentity, CountertermKey,
            graph_key::{self, OccurrenceFeatures},
        },
        delta_key::{self, DeltaKeyMode, DeltaVerify},
        delta_store::{self, StoreGet, StoreMode},
        hedge_poset::ForestNode,
        literal_zero_if_vanishing,
    },
};

/// A complete pole-only UV map with its exact dressed graph and settings.
/// Constructed by R* planning; callers cannot change the identity after planning.
#[derive(Clone)]
pub struct CountertermRequest<'a> {
    graph: Arc<Graph>,
    support: SuBitGraph,
    model: &'a Model,
    settings: Arc<UVgenerationSettings>,
}

impl<'a> CountertermRequest<'a> {
    pub(crate) fn new(
        graph: Arc<Graph>,
        support: SuBitGraph,
        model: &'a Model,
        settings: Arc<UVgenerationSettings>,
    ) -> Self {
        Self {
            graph,
            support,
            model,
            settings,
        }
    }

    pub fn loops(&self) -> usize {
        self.graph.n_loops(&self.support)
    }

    /// Publish through the ordinary cache, preserving its boundary validation.
    pub fn evaluate(&self, cache: &CountertermCache, vakint: &Vakint) -> Result<Atom> {
        cache.counterterm(
            &self.graph,
            &self.support,
            self.model,
            vakint,
            &self.settings,
        )
    }
}

/// An exact static-UV closure shared by parallel discovery workers.
/// Only unique non-ready requests retain their dressed graph families.
/// The mutex protects insertion only; identity and classification run outside it.
#[derive(Default)]
pub struct CountertermPrewarmPlan<'a> {
    state: Mutex<CountertermPrewarmState<'a>>,
}

#[derive(Default)]
struct CountertermPrewarmState<'a> {
    seen: HashSet<Key>,
    levels: BTreeMap<usize, Vec<CountertermRequest<'a>>>,
}

impl<'a> CountertermPrewarmPlan<'a> {
    pub fn unique_requests_by_loop(&self) -> BTreeMap<usize, usize> {
        self.state
            .lock()
            .expect("counterterm discovery poisoned")
            .levels
            .iter()
            .map(|(loops, requests)| (*loops, requests.len()))
            .collect()
    }

    /// Allocation count, not a byte estimate: graph size depends on numerators.
    pub fn retained_graph_families(&self) -> usize {
        self.state
            .lock()
            .expect("counterterm discovery poisoned")
            .levels
            .values()
            .flatten()
            .map(|request| Arc::as_ptr(&request.graph))
            .collect::<HashSet<_>>()
            .len()
    }

    pub fn into_levels(self) -> BTreeMap<usize, Vec<CountertermRequest<'a>>> {
        self.state
            .into_inner()
            .expect("counterterm discovery poisoned")
            .levels
    }
}

/// A failed computation is stored with its whole error, and every lookup
/// replays it as a [`CachedError`]: the classes of its causes survive the
/// cache (a fatal IRH-DEPTH alarm, `rstar::depth::DepthAlarm`, stays one on
/// every later hit; a failed comparison,
/// (`rstar::ir_dirac::diagnostic_failure`), which no fallback may decline
/// (a decorated quotient's R* runs inside the child guard of an IR hard
/// part), stays one too: `rstar::is_diagnostic_failure` reads through the
/// replay), its text unchanged. A transient failure
/// ([`CachedError::is_transient`], an exhausted resource budget) is not
/// stored: the next lookup computes again.
type Entry = Arc<OnceLock<std::result::Result<Atom, CachedError>>>;

/// An error of a cached computation, replayed on a lookup. Its message is the
/// stored error's whole chain (`{:#}`, the text the cache replayed before);
/// [`CachedError::cause`] is the stored error itself, whose chain callers
/// inspect for typed causes.
#[derive(Clone, Debug)]
pub(crate) struct CachedError(Arc<eyre::Report>);

impl CachedError {
    /// The stored error.
    pub(crate) fn cause(&self) -> &eyre::Report {
        &self.0
    }

    /// The replay on a lookup.
    fn replay(&self) -> eyre::Report {
        eyre::Report::new(self.clone())
    }

    /// A transient failure, never kept by the cache: an exhausted resource
    /// budget (a `ResourceLimitError` in the chain, also one replayed).
    pub(crate) fn is_transient(&self) -> bool {
        is_resource_limit(&self.0)
    }
}

/// Tests: `error` as the cache replays a stored error.
#[cfg(test)]
pub(crate) fn replay_for_tests(error: eyre::Report) -> eyre::Report {
    CachedError(Arc::new(error)).replay()
}

/// Whether `error` carries a symbolica `ResourceLimitError` anywhere in its
/// chain, also one inside a replayed [`CachedError`].
pub(crate) fn is_resource_limit(error: &eyre::Report) -> bool {
    use symbolica::utils::resources::ResourceLimitError;
    error.chain().any(|cause| {
        cause.downcast_ref::<ResourceLimitError>().is_some()
            || cause
                .downcast_ref::<CachedError>()
                .is_some_and(|cached| is_resource_limit(cached.cause()))
    })
}

impl std::fmt::Display for CachedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:#}", self.0)
    }
}

impl std::error::Error for CachedError {}

/// What a fill noted about itself, kept with its entry
/// ([`CountertermCache::counterterm_noted`]; the uncontracted ledger's
/// dispositions of the records the fill evaluated).
pub(crate) type CacheNote = Arc<dyn std::any::Any + Send + Sync>;

/// How [`CountertermCache::counterterm_noted`] found its entry.
pub(crate) enum CacheLookup {
    /// This call computed the entry.
    Filled,
    /// The entry was there (or being computed by another call): the note of
    /// its fill, if it left one.
    Hit(Option<CacheNote>),
}

/// The single-flight rounds a lookup that only WAITED on another lookup's
/// transient failure may retry before it is refused.
const WAIT_ROUNDS: usize = 3;

/// IRH-DEPTH round 8 (r47, the three-waiter retry limit): a lookup that
/// waited [`WAIT_ROUNDS`] times on an entry another lookup failed with a
/// transient error (an exhausted resource budget), without ever running its
/// own computation, is REFUSED with this typed error (logged
/// `DELTA_CACHE_RETRIES_EXHAUSTED`). Before round 8 it returned the other
/// lookup's error as its own. The value: never a success; the cause (the
/// other lookup's transient error) is in the chain, so a guarded child
/// attempt declines the whole correction to today's path (closed::guarded,
/// a resource-limit decline) and any other caller fails the diagram.
#[derive(Debug)]
pub(crate) struct RetriesExhausted {
    rounds: usize,
    cause: CachedError,
}

impl std::fmt::Display for RetriesExhausted {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "complete UV counterterm cache: refused after {} single-flight rounds that each ended on another lookup's transient failure, without computing this lookup's own entry",
            self.rounds
        )
    }
}

impl std::error::Error for RetriesExhausted {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.cause)
    }
}

/// A cache key: today's certificate (`Legacy`), or the routing-free graph
/// certificate of an R* child Δ (`Graph`, `GAMMALOOP_DELTA_KEY=graph`).
#[derive(Clone, PartialEq, Eq, Hash)]
enum Key {
    Legacy {
        graph: CountertermKey,
        settings: String,
    },
    Graph {
        graph: CountertermKey,
        settings: String,
    },
}

/// GAMMALOOP_DELTA_STORE: the store key of one graph-key lookup.
struct StoreKey {
    kind: &'static str,
    key_sha: [u8; 32],
    graph_sha: [u8; 32],
    model_fp: [u8; 32],
    context: Arc<delta_store::Context>,
    loops: usize,
    key_bytes: Option<Arc<[u8]>>,
}

/// An accepted store row, ready for this occurrence.
struct StoreHit {
    value: Atom,
    note: Option<CacheNote>,
    features: Option<OccurrenceFeatures>,
    writer_build: String,
    dep_builds: Vec<String>,
}

/// What a lookup is: the graph certificate applies to R* child Δs only.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum KeyKind {
    Standard,
    RstarChildDelta,
    BareInsertion,
    RenormalizedInsertion,
}

impl KeyKind {
    fn insertion(self) -> bool {
        matches!(self, Self::BareInsertion | Self::RenormalizedInsertion)
    }

    fn store_kind(self) -> &'static str {
        match self {
            Self::BareInsertion => "bare-sigma",
            Self::RenormalizedInsertion => "renormalized-sigma",
            _ => delta_store::KIND_RSTAR_CHILD,
        }
    }
}

thread_local! {
    /// Sampled verify (`GAMMALOOP_DELTA_VERIFY=N`): the private cold cache
    /// every lookup of this thread goes to while a fresh value is computed.
    static COLD_CACHE: std::cell::RefCell<Option<Arc<CountertermCache>>> =
        const { std::cell::RefCell::new(None) };
}

/// Optional elapsed-time diagnostics. Initializer time includes recursive child
/// calls; the phase sums therefore are not disjoint process wall time. Reuse
/// measures a non-initializing `get_or_init`, including waiting and handoff, not
/// CPU time. Pending reuse is its subset whose entry was initially incomplete.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize)]
pub struct CacheTimingStats {
    pub calls: u64,
    pub initializers: u64,
    pub reuses: u64,
    pub pending_reuses: u64,
    pub key_ns: u64,
    pub map_lock_wait_ns: u64,
    pub map_lock_hold_ns: u64,
    pub initializer_ns: u64,
    pub store_ns: u64,
    pub reuse_ns: u64,
    pub pending_reuse_ns: u64,
    pub restore_ns: u64,
}

#[derive(Clone, Copy)]
pub(crate) enum CacheTimingPhase {
    Calls,
    Initializers,
    Reuses,
    PendingReuses,
    Key,
    MapLockWait,
    MapLockHold,
    Initializer,
    Store,
    Reuse,
    PendingReuse,
    Restore,
}

const TIMING_FIELDS: usize = CacheTimingPhase::Restore as usize + 1;

#[derive(Default)]
pub(crate) struct CacheTimings {
    values: Option<Box<[AtomicU64; TIMING_FIELDS]>>,
}

impl CacheTimings {
    pub(crate) fn new(enabled: bool) -> Self {
        Self {
            values: enabled.then(|| Box::new(std::array::from_fn(|_| AtomicU64::new(0)))),
        }
    }

    pub(crate) fn stats(&self) -> Option<CacheTimingStats> {
        let values = self.values.as_ref()?;
        let get = |field: CacheTimingPhase| values[field as usize].load(Ordering::Relaxed);
        Some(CacheTimingStats {
            calls: get(CacheTimingPhase::Calls),
            initializers: get(CacheTimingPhase::Initializers),
            reuses: get(CacheTimingPhase::Reuses),
            pending_reuses: get(CacheTimingPhase::PendingReuses),
            key_ns: get(CacheTimingPhase::Key),
            map_lock_wait_ns: get(CacheTimingPhase::MapLockWait),
            map_lock_hold_ns: get(CacheTimingPhase::MapLockHold),
            initializer_ns: get(CacheTimingPhase::Initializer),
            store_ns: get(CacheTimingPhase::Store),
            reuse_ns: get(CacheTimingPhase::Reuse),
            pending_reuse_ns: get(CacheTimingPhase::PendingReuse),
            restore_ns: get(CacheTimingPhase::Restore),
        })
    }

    pub(crate) fn start(&self) -> CacheTimingSample<'_> {
        CacheTimingSample {
            totals: self,
            active: self.values.as_ref().map(|_| {
                let mut values = [0; TIMING_FIELDS];
                values[CacheTimingPhase::Calls as usize] = 1;
                (Instant::now(), values)
            }),
        }
    }
}

pub(crate) struct CacheTimingSample<'a> {
    totals: &'a CacheTimings,
    active: Option<(Instant, [u64; TIMING_FIELDS])>,
}

impl CacheTimingSample<'_> {
    pub(crate) fn lap(&mut self, phase: CacheTimingPhase) {
        if let Some((started, values)) = &mut self.active {
            let now = Instant::now();
            values[phase as usize] += now.duration_since(*started).as_nanos() as u64;
            *started = now;
        }
    }

    pub(crate) fn restart(&mut self) {
        if let Some((started, _)) = &mut self.active {
            *started = Instant::now();
        }
    }

    pub(crate) fn resolve<'a, T>(
        &mut self,
        entry: &'a OnceLock<T>,
        compute: impl FnOnce(&mut Self) -> T,
    ) -> &'a T {
        if self.active.is_none() {
            return entry.get_or_init(|| compute(self));
        }
        let pending = entry.get().is_none();
        let started = Instant::now();
        let mut initialized = false;
        let result = entry.get_or_init(|| {
            initialized = true;
            let (mark, values) = self.active.as_mut().unwrap();
            *mark = Instant::now();
            values[CacheTimingPhase::Initializers as usize] += 1;
            compute(self)
        });
        if !initialized {
            let elapsed = started.elapsed().as_nanos() as u64;
            let (_, values) = self.active.as_mut().unwrap();
            values[CacheTimingPhase::Reuses as usize] += 1;
            values[CacheTimingPhase::Reuse as usize] += elapsed;
            if pending {
                values[CacheTimingPhase::PendingReuses as usize] += 1;
                values[CacheTimingPhase::PendingReuse as usize] += elapsed;
            }
        }
        result
    }
}

impl Drop for CacheTimingSample<'_> {
    fn drop(&mut self) {
        if let Some((_, values)) = &self.active {
            // Publication occurs after every map/entry critical section ends.
            for (total, value) in self.totals.values.as_ref().unwrap().iter().zip(values) {
                if *value != 0 {
                    total.fetch_add(*value, Ordering::Relaxed);
                }
            }
        }
    }
}

/// A shared cache of unprojected local counterterms in canonical graph slots.
///
/// Graph prefactors are excluded. Each entry retains the complete nested pole
/// projections and both auxiliary-mass parameters. Settings and the dressed
/// particle graph form the key; a hit transports its open indices and momenta
/// to the requesting graph. Concurrent requests for one entry share one fill.
#[derive(Default)]
pub struct CountertermCache {
    entries: Mutex<HashMap<Key, Entry>>,
    /// The notes of the fills made by [`Self::counterterm_noted`].
    notes: Mutex<HashMap<Key, CacheNote>>,
    hits: AtomicU64,
    misses: AtomicU64,
    stores: AtomicU64,
    waits: AtomicU64,
    timings: CacheTimings,
    compute_limit: Option<ComputeLimit>,
    /// `GAMMALOOP_DELTA_KEY=graph`: the symmetry features of the occurrence
    /// that filled each graph-key entry.
    fills: Mutex<HashMap<Key, Arc<OccurrenceFeatures>>>,
    /// Verify or dump: the legacy digests seen, for `hit_legacy_would_miss`.
    legacy_seen: Mutex<HashSet<u64>>,
    /// `GAMMALOOP_DELTA_STORE`: the store keys each entry's value depends on.
    store_deps: Mutex<HashMap<Key, Arc<delta_store::Deps>>>,
}

/// At most `permits` concurrent computations of counterterms with at least
/// `min_loops` loops. Their transient algebra expands the complete local
/// numerator, which at five loops reaches gigabytes per computation, while the
/// stored results stay small.
struct ComputeLimit {
    min_loops: usize,
    available: Mutex<usize>,
    released: Condvar,
}

thread_local! {
    /// Whether this thread already holds a compute permit. Child counterterms
    /// of a computation reuse it, so a permit holder never waits for another.
    static HOLDS_COMPUTE_PERMIT: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// A held compute permit, returned on drop.
struct ComputePermit<'a>(&'a ComputeLimit);

impl ComputeLimit {
    /// A permit for a computation of `loops` loops, unless the limit does not
    /// apply or this thread already holds one. The permit is taken before the
    /// cache entry is claimed: a thread that waits here owns no unfinished
    /// entry, and every claimed large entry belongs to a permit holder.
    fn acquire(&self, loops: usize) -> Option<ComputePermit<'_>> {
        if loops < self.min_loops || HOLDS_COMPUTE_PERMIT.get() {
            return None;
        }
        let mut available = self.available.lock().expect("counterterm permits poisoned");
        while *available == 0 {
            available = self
                .released
                .wait(available)
                .expect("counterterm permits poisoned");
        }
        *available -= 1;
        HOLDS_COMPUTE_PERMIT.set(true);
        Some(ComputePermit(self))
    }
}

impl Drop for ComputePermit<'_> {
    fn drop(&mut self) {
        HOLDS_COMPUTE_PERMIT.set(false);
        *self
            .0
            .available
            .lock()
            .expect("counterterm permits poisoned") += 1;
        self.0.released.notify_one();
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CountertermCacheStats {
    pub hits: u64,
    pub misses: u64,
    pub stores: u64,
    pub waits: u64,
    pub entries: usize,
    pub timing: Option<CacheTimingStats>,
}

impl CountertermCache {
    /// Allow at most `permits` concurrent computations of counterterms with at
    /// least `min_loops` loops. Results are unchanged; only scheduling is.
    pub fn with_compute_limit(mut self, permits: usize, min_loops: usize) -> Self {
        self.compute_limit = Some(ComputeLimit {
            min_loops,
            available: Mutex::new(permits.max(1)),
            released: Condvar::new(),
        });
        self
    }

    /// Enable optional cache-phase accounting; disabled caches take no clocks.
    pub fn with_timing(mut self, enabled: bool) -> Self {
        self.timings = CacheTimings::new(enabled);
        self
    }

    pub fn stats(&self) -> CountertermCacheStats {
        CountertermCacheStats {
            hits: self.hits.load(Ordering::Relaxed),
            misses: self.misses.load(Ordering::Relaxed),
            stores: self.stores.load(Ordering::Relaxed),
            waits: self.waits.load(Ordering::Relaxed),
            entries: self
                .entries
                .lock()
                .expect("counterterm cache poisoned")
                .len(),
            timing: self.timings.stats(),
        }
    }

    /// Number of entries and bytes of their stored expressions (diagnostics).
    pub fn footprint(&self) -> (usize, usize) {
        let entries = self.entries.lock().expect("counterterm cache poisoned");
        let bytes = entries
            .values()
            .filter_map(|entry| entry.get())
            .map(|value| match value {
                Ok(atom) => atom.as_view().get_byte_size(),
                Err(error) => error.to_string().len(),
            })
            .sum();
        (entries.len(), bytes)
    }

    /// Extend one global static closure from an existing discovery worker.
    /// Concurrent workers claim each exact identity once, without waiting for its
    /// classification. All workers must finish successfully before consuming the
    /// plan. Only key and result insertion hold the plan mutex; graph algebra and
    /// recursion stay on the worker. This creates no threads or integrals.
    /// Completing each returned loop level before the next prevents recursive
    /// misses for this static set. Decorated R* quotients remain outside it.
    pub fn extend_prewarm_plan<'a>(
        &self,
        plan: &CountertermPrewarmPlan<'a>,
        requests: impl IntoIterator<Item = CountertermRequest<'a>>,
    ) -> Result<()> {
        let _stage = StageSpan::new("gamma.delta.prewarm_plan");
        let mut pending = requests.into_iter().collect::<Vec<_>>();
        while let Some(request) = pending.pop() {
            let identity =
                CountertermIdentity::new(&request.graph, &request.support, request.model)?;
            let key = Key::Legacy {
                graph: identity.key,
                settings: serde_json::to_string(&*request.settings)?,
            };
            // Clone before acquiring the mutex to avoid copying the graph
            // identity while holding the shared state lock.
            let discovery_key = key.clone();
            let claimed = plan
                .state
                .lock()
                .expect("counterterm discovery poisoned")
                .seen
                .insert(discovery_key);
            if !claimed {
                continue;
            }
            let entry = self
                .entries
                .lock()
                .expect("counterterm cache poisoned")
                .get(&key)
                .cloned();
            if let Some(result) = entry.as_ref().and_then(|entry| entry.get()) {
                result.as_ref().map_err(CachedError::replay)?;
                continue;
            }
            let loops = request.loops();
            let graph = &request.graph;
            let catalogue = graph.classified_spinneys(
                &request.support,
                request.model,
                &request.settings,
                &graph.loop_momentum_basis,
            );
            ensure!(
                catalogue
                    .iter()
                    .all(|spinney| spinney.renormalization_scheme == ApproximationType::PolePart),
                "complete UV-counterterm prewarming requires pole-part subtraction"
            );
            // The real miss path returns zero or an error before recursing
            // when the catalogue has no overall support. Keep that request so
            // compute() still validates it; do not discover unused descendants.
            let components = if catalogue
                .iter()
                .any(|spinney| spinney.filter() == &request.support)
            {
                catalogue
                    .iter()
                    .filter(|spinney| {
                        spinney.filter() != &request.support && !spinney.filter().is_empty()
                    })
                    .flat_map(|spinney| graph.connected_components(spinney.filter()))
                    .collect::<BTreeSet<_>>()
            } else {
                BTreeSet::new()
            };
            for support in components {
                ensure!(
                    graph.n_loops(&support) < loops,
                    "a proper UV child must consume fewer loops than its parent"
                );
                pending.push(CountertermRequest::new(
                    Arc::clone(graph),
                    support,
                    request.model,
                    Arc::clone(&request.settings),
                ));
            }
            plan.state
                .lock()
                .expect("counterterm discovery poisoned")
                .levels
                .entry(loops)
                .or_default()
                .push(request);
        }
        Ok(())
    }

    pub fn prewarm_levels<'a>(
        &self,
        requests: impl IntoIterator<Item = CountertermRequest<'a>>,
    ) -> Result<BTreeMap<usize, Vec<CountertermRequest<'a>>>> {
        let plan = CountertermPrewarmPlan::default();
        self.extend_prewarm_plan(&plan, requests)?;
        Ok(plan.into_levels())
    }

    pub(crate) fn counterterm(
        &self,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<Atom> {
        // GAMMALOOP_MASSIVE_ROUTES=p-only: every massive R request (family
        // graph with a reference mass, decorated quotient, prewarm) passes here.
        crate::uv::rstar::MassiveRoutes::from_environment().catch_all(
            "massive_counterterm_cache",
            concat!(file!(), ":", line!()),
            || {
                format!(
                    "graph={} child_edges={} loops={} why=a massive R request reached the counterterm cache",
                    graph.name,
                    graph
                        .iter_edges_of(support)
                        .map(|(_, edge, _)| usize::from(edge).to_string())
                        .collect::<Vec<_>>()
                        .join(","),
                    graph.n_loops(support)
                )
            },
        )?;
        self.lookup(
            graph,
            support,
            model,
            serde_json::to_string(settings)?,
            KeyKind::Standard,
            || self.compute(graph, support, model, vakint, settings),
        )
    }

    /// A complete counterterm of `support` computed by another route.
    ///
    /// `settings_key` must identify that route and its settings; it must differ
    /// from every key of [`Self::counterterm`], whose raw entries keep both
    /// auxiliary masses apart for massive parents. The entry shares this
    /// cache's single-flight fill, compute permits, canonical boundary
    /// transport and statistics. `compute` runs only on a miss.
    pub(crate) fn counterterm_with(
        &self,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        settings_key: String,
        compute: impl FnOnce() -> Result<Atom>,
    ) -> Result<Atom> {
        self.lookup(
            graph,
            support,
            model,
            settings_key,
            KeyKind::Standard,
            compute,
        )
    }

    /// [`Self::counterterm_with`] for a lookup of `kind`: an R* child Δ
    /// takes the graph certificate under `GAMMALOOP_DELTA_KEY=graph`.
    pub(crate) fn counterterm_with_kind(
        &self,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        settings_key: String,
        kind: KeyKind,
        compute: impl FnOnce() -> Result<Atom>,
    ) -> Result<Atom> {
        self.lookup(graph, support, model, settings_key, kind, compute)
    }

    /// [`Self::counterterm_with`] whose fill leaves a note with its entry
    /// (stored before the entry is published, so a waiting hit finds it):
    /// the value and whether this call filled the entry or hit it, with the
    /// note of the fill (`None` for a fill that left none, e.g. one made by
    /// [`Self::counterterm_with`]). The uncontracted ledger's Δ-cache hits
    /// (`GAMMALOOP_UNCONTRACTED_NUMERATOR`, U7a delta).
    pub(crate) fn counterterm_noted(
        &self,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        settings_key: String,
        kind: KeyKind,
        compute: impl FnOnce() -> Result<(Atom, CacheNote)>,
    ) -> Result<(Atom, CacheLookup)> {
        self.lookup_noted(
            graph,
            support,
            model,
            settings_key,
            kind,
            || compute().map(|(value, note)| (value, Some(note))),
            true,
        )
    }

    fn lookup(
        &self,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        settings_key: String,
        kind: KeyKind,
        compute: impl FnOnce() -> Result<Atom>,
    ) -> Result<Atom> {
        self.lookup_noted(
            graph,
            support,
            model,
            settings_key,
            kind,
            || compute().map(|value| (value, None)),
            false,
        )
        .map(|(value, _)| value)
    }

    /// The identity and key of a lookup: the graph certificate for an R*
    /// child Δ under `GAMMALOOP_DELTA_KEY=graph` (unless it falls back; also
    /// a noted lookup of the uncontracted ledger: the fill's note is about the
    /// child's computation, which the certificate identifies), the
    /// legacy one otherwise. The bool: the graph certificate is used.
    fn identity_and_key(
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        settings_key: String,
        kind: KeyKind,
    ) -> Result<(CountertermIdentity, Key, bool)> {
        if kind.insertion() {
            let identity = CountertermIdentity::for_insertion(graph, support, model)?;
            let key = Key::Graph {
                graph: identity.key.clone(),
                settings: settings_key,
            };
            return Ok((identity, key, false));
        }
        if kind == KeyKind::RstarChildDelta && delta_key::mode() == DeltaKeyMode::Graph {
            if crate::uv::rstar::ir_hard_fixed_mass_scope() {
                delta_key::add(&delta_key::C.scoped_flag_lookup, 1);
                debug_assert!(
                    false,
                    "an R* child Δ lookup inside the IR-hard fixed-mass scope"
                );
            } else {
                match CountertermIdentity::for_graph_key(graph, support, model) {
                    Ok(Ok(identity)) => {
                        let key = Key::Graph {
                            graph: identity.key.clone(),
                            settings: settings_key,
                        };
                        return Ok((identity, key, true));
                    }
                    Ok(Err(graph_key::Fallback::NonIncident)) => {
                        delta_key::add(&delta_key::C.key_fallback_nonincident, 1)
                    }
                    Ok(Err(graph_key::Fallback::Coefficient)) => {
                        delta_key::add(&delta_key::C.key_fallback_coef, 1)
                    }
                    Err(error) => {
                        if delta_key::C
                            .key_fallback_error
                            .fetch_add(1, Ordering::Relaxed)
                            < 5
                        {
                            eprintln!(
                                "DELTA_KEY_FALLBACK_ERROR graph={} support={} error={error:#}",
                                graph.name,
                                support.string_label()
                            );
                        }
                    }
                }
            }
        }
        let identity = CountertermIdentity::new(graph, support, model)?;
        let key = Key::Legacy {
            graph: identity.key.clone(),
            settings: settings_key,
        };
        Ok((identity, key, false))
    }

    /// The lookup; `noted`: the fill's note is stored with its entry and a
    /// hit returns it (no note is kept or read otherwise). `kind`: an R*
    /// child Δ takes the graph certificate under `GAMMALOOP_DELTA_KEY=graph`.
    #[allow(clippy::too_many_arguments)]
    fn lookup_noted(
        &self,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        settings_key: String,
        kind: KeyKind,
        compute: impl FnOnce() -> Result<(Atom, Option<CacheNote>)>,
        noted: bool,
    ) -> Result<(Atom, CacheLookup)> {
        // Sampled verify: a fresh value is computed in a private cold cache.
        if let Some(cold) = COLD_CACHE.with(|cold| cold.borrow().clone())
            && !std::ptr::eq(Arc::as_ptr(&cold), self)
        {
            return cold.lookup_noted(graph, support, model, settings_key, kind, compute, noted);
        }
        let _stage = StageSpan::new("gamma.delta.lookup");
        let stats = delta_key::stats_enabled();
        let key_started = stats.then(Instant::now);
        let mut timing = self.timings.start();
        let identity_stage = StageSpan::new("gamma.delta.identity");
        let (mut identity, key, graph_key) =
            Self::identity_and_key(graph, support, model, settings_key, kind)?;
        let noted_key = noted.then(|| key.clone());
        timing.lap(CacheTimingPhase::Key);
        drop(identity_stage);
        let loops = stats.then(|| graph.n_loops(support));
        let verify_mode = if graph_key {
            delta_key::verify()
        } else {
            DeltaVerify::Off
        };
        let in_cold = COLD_CACHE.with(|cold| cold.borrow().is_some());
        let key_digest = (graph_key && (delta_key::dump() || verify_mode != DeltaVerify::Off))
            .then(|| identity.key_digest());
        // Verify or dump: the legacy certificate beside the graph one.
        let legacy_new = if graph_key && (delta_key::dump() || verify_mode != DeltaVerify::Off) {
            use std::hash::{Hash, Hasher};
            let legacy = CountertermIdentity::new(graph, support, model)?;
            let mut hasher = std::collections::hash_map::DefaultHasher::new();
            legacy.key.hash(&mut hasher);
            if let Key::Graph { settings, .. } = &key {
                settings.hash(&mut hasher);
            }
            let digest = hasher.finish();
            identity.set_legacy_digest(digest);
            Some(
                self.legacy_seen
                    .lock()
                    .expect("counterterm cache poisoned")
                    .insert(digest),
            )
        } else {
            None
        };
        if let Some(started) = key_started {
            delta_key::add(&delta_key::C.key_ns, started.elapsed().as_nanos() as u64);
            delta_key::add(&delta_key::C.lookups, 1);
        }
        // A transient error (an exhausted resource budget, e.g. the scoped
        // budget of a guarded child attempt) is never kept: the entry is
        // removed, the computing lookup returns the error unstored, and a
        // lookup that only waited on that entry computes again under its own
        // budget (MATRIX-R45 review, nested-and-cache finding 1).
        let mut compute = Some(compute);
        let mut rounds = 0_usize;
        // The filler's own value when the open-line transport cannot store it.
        let mut unsupported_value: Option<Atom> = None;
        // GAMMALOOP_DELTA_STORE (graph key, not in a private cold cache): the
        // store; `from_store`: the filler took the value from it.
        let store = if (graph_key || kind.insertion()) && !in_cold {
            delta_store::installed()
        } else {
            None
        };
        let mut from_store = false;
        let (entry, existing, filled) = loop {
            rounds += 1;
            let (entry, existing) = {
                let _stage = StageSpan::new("gamma.delta.cache_map");
                let mut entries = self.entries.lock().expect("counterterm cache poisoned");
                timing.lap(CacheTimingPhase::MapLockWait);
                let existing = entries.contains_key(&key);
                let entry = entries.entry(key.clone()).or_default().clone();
                drop(entries);
                timing.lap(CacheTimingPhase::MapLockHold);
                (entry, existing)
            };
            if existing {
                self.hits.fetch_add(1, Ordering::Relaxed);
                if entry.get().is_none() {
                    self.waits.fetch_add(1, Ordering::Relaxed);
                }
            }
            let _permit = if entry.get().is_none() {
                let _stage = StageSpan::new("gamma.delta.compute_permit_wait");
                self.compute_limit
                    .as_ref()
                    .and_then(|limit| limit.acquire(graph.n_loops(support)))
            } else {
                None
            };
            let resolve_stage = StageSpan::new("gamma.delta.cache_resolve");
            let mut computed = false;
            let result = timing.resolve(&entry, |timing| {
                computed = true;
                // GAMMALOOP_DELTA_STORE: the single filler reads the store
                // first; the fill's nested lookups add their store keys.
                let store_key = store.and_then(|_| self.store_key(kind, &key, graph, support));
                if let (Some((mode, backing)), Some(store_key)) = (store, &store_key) {
                    match self.store_read(
                        mode,
                        backing,
                        store_key,
                        &identity,
                        graph,
                        support,
                        model,
                        noted,
                        &mut compute,
                    ) {
                        Ok(Some(hit)) => {
                            if let (Some(key), Some(note)) = (&noted_key, hit.note) {
                                self.notes
                                    .lock()
                                    .expect("counterterm cache notes poisoned")
                                    .insert(key.clone(), note);
                            }
                            if let Some(features) = hit.features {
                                self.fills
                                    .lock()
                                    .expect("counterterm cache poisoned")
                                    .insert(key.clone(), Arc::new(features));
                            }
                            let mut deps = delta_store::Deps::default();
                            deps.keys.insert(store_key.key_sha);
                            deps.builds.insert(hit.writer_build);
                            deps.builds.extend(hit.dep_builds);
                            self.store_deps
                                .lock()
                                .expect("counterterm cache poisoned")
                                .insert(key.clone(), Arc::new(deps));
                            from_store = true;
                            self.stores.fetch_add(1, Ordering::Relaxed);
                            timing.lap(CacheTimingPhase::Store);
                            return Ok(hit.value);
                        }
                        Ok(None) => {}
                        Err(error) => return Err(CachedError(Arc::new(error))),
                    }
                }
                let collector = store.map(|_| delta_store::Collector::push());
                self.misses.fetch_add(1, Ordering::Relaxed);
                let fill_timer = loops.map(delta_key::FillTimer::start);
                let fill_started = Instant::now();
                let result = compute.take().expect("one computation per lookup")();
                let compute_s = fill_started.elapsed().as_secs_f64();
                drop(fill_timer);
                let deps = collector.map(delta_store::Collector::finish);
                timing.lap(CacheTimingPhase::Initializer);
                let _stage = StageSpan::new("gamma.delta.cache_store");
                let result = result
                    .and_then(|(value, note)| {
                        let note_json = store_key
                            .as_ref()
                            .and(note.as_ref())
                            .and_then(|note| crate::uv::rstar::note_to_json(note.as_ref()));
                        if let (Some(key), Some(note)) = (&noted_key, note) {
                            self.notes
                                .lock()
                                .expect("counterterm cache notes poisoned")
                                .insert(key.clone(), note);
                        }
                        if let Some(mut inherited) = deps.clone() {
                            // A parent depends on this child even on its first (cold) fill.
                            if let Some(child) = &store_key {
                                inherited.keys.insert(child.key_sha);
                            }
                            self.store_deps
                                .lock()
                                .expect("counterterm cache poisoned")
                                .insert(key.clone(), Arc::new(inherited));
                        }
                        let canonical = identity.canonicalize_delta(&value);
                        if let (Some((_, backing)), Some(store_key), Ok(canonical)) =
                            (store, &store_key, &canonical)
                        {
                            self.store_write(
                                backing,
                                store_key,
                                &identity,
                                graph,
                                support,
                                model,
                                canonical,
                                &value,
                                note_json.filter(|_| noted),
                                compute_s,
                                deps.as_ref(),
                            )?;
                        }
                        if graph_key {
                            match &canonical {
                                Err(error) if graph_key::is_transport_unsupported(error) => {
                                    unsupported_value = Some(value);
                                }
                                Ok(_) => {
                                    if let Some(transport) = identity.graph_transport() {
                                        self.fills
                                            .lock()
                                            .expect("counterterm cache poisoned")
                                            .insert(
                                                key.clone(),
                                                Arc::new(transport.features.clone()),
                                            );
                                    }
                                }
                                Err(_) => {}
                            }
                        }
                        canonical
                    })
                    .inspect(|_| {
                        self.stores.fetch_add(1, Ordering::Relaxed);
                    })
                    .map_err(|error| CachedError(Arc::new(error)));
                timing.lap(CacheTimingPhase::Store);
                result
            });
            drop(resolve_stage);
            if let Err(error) = result
                && error.is_transient()
            {
                let mut entries = self.entries.lock().expect("counterterm cache poisoned");
                if entries
                    .get(&key)
                    .is_some_and(|stored| Arc::ptr_eq(stored, &entry))
                {
                    entries.remove(&key);
                }
                drop(entries);
                if computed || compute.is_none() {
                    return Err(error.replay());
                }
                if rounds >= WAIT_ROUNDS {
                    eprintln!(
                        "DELTA_CACHE_RETRIES_EXHAUSTED rounds={rounds} loops={} support={} refused=true cause={}",
                        graph.n_loops(support),
                        support.string_label(),
                        format!("{error}").chars().take(200).collect::<String>()
                    );
                    return Err(eyre::Report::new(RetriesExhausted {
                        rounds,
                        cause: error.clone(),
                    }));
                }
                continue;
            }
            break (entry, existing, computed);
        };
        // A store hit is a hit for the ledger and the counters.
        let store_hit = filled && from_store;
        if store.is_some()
            && let Some(deps) = self
                .store_deps
                .lock()
                .expect("counterterm cache poisoned")
                .get(&key)
                .cloned()
        {
            delta_store::add_to_enclosing(&deps);
        }
        if let Some(loops) = loops {
            delta_key::count_loop(loops, !filled || store_hit);
            let (lookups, hits) = if graph_key {
                (&delta_key::C.graph_lookups, &delta_key::C.graph_hits)
            } else {
                (&delta_key::C.legacy_lookups, &delta_key::C.legacy_hits)
            };
            delta_key::add(lookups, 1);
            if !filled || store_hit {
                delta_key::add(hits, 1);
            } else if graph_key {
                delta_key::add(&delta_key::C.graph_misses, 1);
            } else {
                delta_key::add(&delta_key::C.legacy_misses, 1);
            }
        }
        if !filled && legacy_new == Some(true) {
            delta_key::add(&delta_key::C.hit_legacy_would_miss, 1);
        }
        let result = entry.get().expect("resolved above");
        // GAMMALOOP_DELTA_KEY=graph: an entry whose open-line transport the
        // value does not support serves no one; every lookup computes.
        if graph_key
            && let Err(error) = result
            && graph_key::is_transport_unsupported(error.cause())
        {
            delta_key::add(&delta_key::C.transport_unsupported, 1);
            if let Some(value) = unsupported_value.take() {
                return Ok((value, CacheLookup::Filled));
            }
            let (value, _) = compute.take().expect("a hit keeps its computation")()?;
            return Ok((value, CacheLookup::Filled));
        }
        let canonical = result.as_ref().map_err(CachedError::replay)?;
        debug_tags!(#uv, #integrated, #cache;
            stage = "uv_counterterm_cache_lookup",
            hit = existing,
            loops = graph.n_loops(support),
            support = %support.string_label(),
            "Complete UV counterterm cache lookup"
        );
        // GAMMALOOP_DELTA_KEY=graph: the features of a hit, and its verification.
        if graph_key && !filled {
            let verified = match verify_mode {
                DeltaVerify::On => !in_cold,
                DeltaVerify::Sample(n) => !in_cold && delta_key::sampled(identity.key_digest(), n),
                DeltaVerify::Off => false,
            };
            if verified {
                let started = Instant::now();
                let fresh = match verify_mode {
                    DeltaVerify::Sample(_) => {
                        let cold = Arc::new(CountertermCache::default());
                        let previous = COLD_CACHE.with(|slot| slot.replace(Some(cold)));
                        let fresh = compute.take().expect("a hit keeps its computation")();
                        COLD_CACHE.with(|slot| slot.replace(previous));
                        fresh
                    }
                    _ => compute.take().expect("a hit keeps its computation")(),
                };
                let fresh = match fresh {
                    Ok((value, _)) => value,
                    Err(error) => {
                        delta_key::add(&delta_key::C.verify_errors, 1);
                        return Err(error.wrap_err("DELTA_VERIFY: the fresh computation failed"));
                    }
                };
                match identity.verify_agree(graph, support, model, canonical, &fresh) {
                    Ok((identical, true)) => {
                        delta_key::add(&delta_key::C.verify_checked, 1);
                        delta_key::add(
                            if identical {
                                &delta_key::C.verify_identical
                            } else {
                                &delta_key::C.verify_equal_nf
                            },
                            1,
                        );
                    }
                    Ok((_, false)) => {
                        delta_key::add(&delta_key::C.verify_mismatch, 1);
                        let fill = self
                            .fills
                            .lock()
                            .expect("counterterm cache poisoned")
                            .get(&key)
                            .cloned();
                        let message = format!(
                            "DELTA_VERIFY_MISMATCH key={:016x} graph={} support={} loops={} hit=[{}] fill=[{}]\nstored(transported): {}\nfresh: {}",
                            identity.key_digest(),
                            graph.name,
                            support.string_label(),
                            graph.n_loops(support),
                            identity
                                .graph_transport()
                                .map(|transport| transport.features.summary())
                                .unwrap_or_default(),
                            fill.map(|fill| fill.summary()).unwrap_or_default(),
                            identity
                                .restore_delta(canonical)
                                .map(|value| value.to_string())
                                .unwrap_or_else(|error| format!("{error:#}")),
                            fresh
                        );
                        eprintln!("{message}");
                        return Err(eyre!(message));
                    }
                    Err(error) => {
                        delta_key::add(&delta_key::C.verify_errors, 1);
                        return Err(error.wrap_err("DELTA_VERIFY: the comparison failed"));
                    }
                }
                delta_key::add(&delta_key::C.verify_ns, started.elapsed().as_nanos() as u64);
            }
            if let Some(transport) = identity.graph_transport() {
                let fill = self
                    .fills
                    .lock()
                    .expect("counterterm cache poisoned")
                    .get(&key)
                    .cloned();
                if let Some(fill) = fill {
                    delta_key::count_features(&transport.features.shown(&fill), verified);
                }
            }
        }
        if let (true, Some(digest)) = (delta_key::dump(), key_digest) {
            eprintln!(
                "DELTA_KEY key={digest:016x} graph={} support={} loops={} {} hit={}",
                graph.name,
                support.string_label(),
                graph.n_loops(support),
                identity
                    .graph_transport()
                    .map(|transport| transport.features.summary())
                    .unwrap_or_default(),
                if store_hit {
                    "store"
                } else if filled {
                    "miss"
                } else if existing {
                    "mem"
                } else {
                    "wait"
                }
            );
        }
        timing.restart();
        let _stage = StageSpan::new("gamma.delta.restore_boundary");
        let result = match identity.restore_delta(canonical) {
            Err(error) if graph_key && graph_key::is_transport_unsupported(&error) => {
                delta_key::add(&delta_key::C.transport_unsupported, 1);
                match compute.take() {
                    Some(compute) => compute().map(|(value, _)| value),
                    None => Err(error),
                }
            }
            other => other,
        };
        timing.lap(CacheTimingPhase::Restore);
        let lookup = if filled && !store_hit {
            CacheLookup::Filled
        } else {
            CacheLookup::Hit(noted_key.and_then(|key| {
                self.notes
                    .lock()
                    .expect("counterterm cache notes poisoned")
                    .get(&key)
                    .cloned()
            }))
        };
        result.map(|value| (value, lookup))
    }

    /// GAMMALOOP_DELTA_STORE: the store key of a graph-key lookup (`None`:
    /// below the store's minimum loop count, a certificate that is not
    /// portable, or no context: never persisted, counted).
    fn store_key(
        &self,
        kind: KeyKind,
        key: &Key,
        graph: &Graph,
        support: &SuBitGraph,
    ) -> Option<StoreKey> {
        let Key::Graph {
            graph: certificate,
            settings,
        } = key
        else {
            return None;
        };
        let loops = graph.n_loops(support);
        if loops < delta_store::min_loops(if kind.insertion() { 1 } else { 2 }) {
            delta_store::add(&delta_store::S.store_below_min_loops, 1);
            return None;
        }
        let started = Instant::now();
        let Some((model_fp, graph_bytes)) = certificate.portable_bytes() else {
            delta_store::add(&delta_store::S.store_not_portable, 1);
            return None;
        };
        let context = match delta_store::context(settings) {
            Ok(context) => context,
            Err(error) => {
                if delta_store::S
                    .store_context_error
                    .fetch_add(1, Ordering::Relaxed)
                    == 0
                {
                    eprintln!("DELTA_STORE_CONTEXT_ERROR {error:#}");
                }
                return None;
            }
        };
        let (key_sha, graph_sha, key_bytes) =
            delta_store::key_sha(kind.store_kind(), &model_fp, &context.sha, &graph_bytes);
        delta_store::add(&delta_store::S.key_ns, started.elapsed().as_nanos() as u64);
        Some(StoreKey {
            kind: kind.store_kind(),
            key_sha,
            graph_sha,
            model_fp,
            context,
            loops,
            key_bytes: delta_store::store_keys().then(|| Arc::from(key_bytes)),
        })
    }

    /// The filler's store read (design §4.5, §5.1): the transported-ready
    /// canonical value of an accepted row, or `None` (a miss). In mode
    /// `verify`, and for the salted audit sample, the hit is recomputed in a
    /// private cold cache and compared; a disagreement is an error.
    #[allow(clippy::too_many_arguments)]
    fn store_read<F>(
        &self,
        mode: StoreMode,
        backing: &dyn delta_store::DeltaBackingStore,
        store_key: &StoreKey,
        identity: &CountertermIdentity,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        noted: bool,
        compute: &mut Option<F>,
    ) -> Result<Option<StoreHit>>
    where
        F: FnOnce() -> Result<(Atom, Option<CacheNote>)>,
    {
        use delta_store::S;
        delta_store::add(&S.store_lookups, 1);
        let started = Instant::now();
        let got = backing.get(&store_key.key_sha);
        delta_store::add(&S.read_ns, started.elapsed().as_nanos() as u64);
        let row = match got {
            StoreGet::Hit(row) => row,
            StoreGet::Miss => {
                delta_store::add(&S.store_misses, 1);
                delta_store::count_loops(store_key.loops, false, false);
                return Ok(None);
            }
            StoreGet::Error(reason) => {
                if S.read_error.fetch_add(1, Ordering::Relaxed) < 5 {
                    eprintln!("DELTA_STORE_READ_ERROR {reason}");
                }
                if mode == StoreMode::Verify {
                    return Err(eyre!("DELTA_STORE_READ_ERROR (mode verify) {reason}"));
                }
                delta_store::add(&S.store_misses, 1);
                delta_store::count_loops(store_key.loops, false, false);
                return Ok(None);
            }
        };
        let miss = |counter: &AtomicU64| {
            delta_store::add(counter, 1);
            delta_store::add(&S.store_misses, 1);
            delta_store::count_loops(store_key.loops, false, false);
            Ok(None)
        };
        let started = Instant::now();
        let value = match delta_store::import(&row) {
            Ok(value) => value,
            Err("value_format") => return miss(&S.reject_value_format),
            Err("blob_sha") => return miss(&S.reject_blob_sha),
            Err("value_sha") => return miss(&S.reject_value_sha),
            Err(_) => return miss(&S.reject_import),
        };
        delta_store::add(&S.import_ns, started.elapsed().as_nanos() as u64);
        let note = match (&row.note, noted) {
            (_, false) => None,
            (Some(text), true) => match crate::uv::rstar::note_from_json(text) {
                Some(note) => Some(note),
                None => return miss(&S.reject_note),
            },
            (None, true) => return miss(&S.reject_note),
        };
        // the value must restore to this occurrence (open-line transport)
        if identity.restore_delta(&value).is_err() {
            return miss(&S.reject_restore);
        }
        let verify_all = mode == StoreMode::Verify || delta_key::verify() == DeltaVerify::On;
        let audited = verify_all
            || (delta_store::audit_every() > 0
                && delta_key::sampled(
                    u64::from_le_bytes(store_key.key_sha[..8].try_into().expect("8 bytes")),
                    delta_store::audit_every(),
                ));
        if audited && let Some(compute) = compute.take() {
            let started = Instant::now();
            let cold = Arc::new(CountertermCache::default());
            let previous = COLD_CACHE.with(|slot| slot.replace(Some(cold)));
            let fresh = compute();
            COLD_CACHE.with(|slot| slot.replace(previous));
            let (fresh, _) = fresh.map_err(|error| {
                error.wrap_err("DELTA_STORE_VERIFY: the fresh computation failed")
            })?;
            let agree = identity.verify_agree(graph, support, model, &value, &fresh)?;
            delta_store::add(&S.verify_ns, started.elapsed().as_nanos() as u64);
            if verify_all {
                delta_store::add(&S.verify_checked, 1);
            } else {
                delta_store::add(&S.audit_checked, 1);
            }
            match agree {
                (true, true) => delta_store::add(&S.verify_identical, 1),
                (false, true) => delta_store::add(&S.verify_equal_nf, 1),
                (_, false) => {
                    delta_store::add(&S.verify_mismatch, 1);
                    delta_store::add(&S.disputes, 1);
                    let canonical = identity.canonicalize_delta(&fresh)?;
                    let record = self.store_record(
                        store_key, identity, graph, support, &canonical, None, 0.0, None,
                    );
                    let message = format!(
                        "DELTA_STORE_VERIFY_MISMATCH key_sha={} graph_sha={} graph={} support={} loops={} writer_build={} writer_run_id={}\nstored(transported): {}\nfresh: {}",
                        delta_store::hex(&store_key.key_sha),
                        delta_store::hex(&store_key.graph_sha),
                        graph.name,
                        support.string_label(),
                        store_key.loops,
                        row.writer_build,
                        row.writer_run_id,
                        identity
                            .restore_delta(&value)
                            .map(|value| value.to_string())
                            .unwrap_or_else(|error| format!("{error:#}")),
                        fresh
                    );
                    eprintln!("{message}");
                    if let Ok(record) = record {
                        backing.event(delta_store::StoreEvent::Dispute {
                            kind: "verify",
                            key_sha: store_key.key_sha,
                            tier: row.tier,
                            stored_value_sha: row.value_sha,
                            record: Box::new(record),
                            quarantine: true,
                        });
                    }
                    return Err(eyre!(message));
                }
            }
        }
        delta_store::add(&S.store_hits, 1);
        delta_store::add(&S.hit_compute_s_saved_ms, (row.compute_s * 1e3) as u64);
        delta_store::count_loops(store_key.loops, true, false);
        delta_store::record_hit(backing, &store_key.key_sha, &row, audited)?;
        let features = serde_json::from_str::<OccurrenceFeatures>(&row.writer_transport).ok();
        // the symmetry features this hit shows against the writer's occurrence
        if let (Some(writer), Some(transport)) = (&features, identity.graph_transport()) {
            delta_key::count_features(&transport.features.shown(writer), audited);
        }
        Ok(Some(StoreHit {
            value,
            note,
            features,
            writer_build: row.writer_build,
            dep_builds: row.dep_builds,
        }))
    }

    /// The record of a canonical value for the store.
    #[allow(clippy::too_many_arguments)]
    fn store_record(
        &self,
        store_key: &StoreKey,
        identity: &CountertermIdentity,
        graph: &Graph,
        support: &SuBitGraph,
        canonical: &Atom,
        note: Option<String>,
        compute_s: f64,
        deps: Option<&delta_store::Deps>,
    ) -> std::result::Result<delta_store::DeltaRecord, &'static str> {
        let started = Instant::now();
        let exported = delta_store::export(canonical)?;
        delta_store::add(
            &delta_store::S.export_ns,
            started.elapsed().as_nanos() as u64,
        );
        Ok(delta_store::DeltaRecord {
            kind: store_key.kind,
            key_sha: store_key.key_sha,
            context_sha: store_key.context.sha,
            context_json: store_key.context.json.clone(),
            graph_sha: store_key.graph_sha,
            model_fp: store_key.model_fp,
            child_loops: store_key.loops,
            is_zero: canonical.is_zero(),
            value_format: delta_store::value_format().to_owned(),
            blob_sha: exported.blob_sha,
            value_sha: exported.value_sha,
            value_export: exported.export,
            note,
            writer_transport: identity
                .graph_transport()
                .and_then(|transport| serde_json::to_string(&transport.features).ok())
                .unwrap_or_default(),
            writer_graph: graph.name.to_string(),
            writer_support: support.string_label(),
            compute_s,
            deps: deps
                .map(|deps| deps.keys.iter().copied().collect())
                .unwrap_or_default(),
            dep_builds: deps
                .map(|deps| deps.builds.iter().cloned().collect())
                .unwrap_or_default(),
            key_bytes: store_key.key_bytes.clone(),
            text: delta_store::store_text().then_some(exported.text),
        })
    }

    /// The filler's store write (design §5.1-5.2): one synchronous insert.
    /// A row another writer stored first is compared with this value in
    /// this occurrence; a disagreement is a dispute and an error.
    #[allow(clippy::too_many_arguments)]
    fn store_write(
        &self,
        backing: &dyn delta_store::DeltaBackingStore,
        store_key: &StoreKey,
        identity: &CountertermIdentity,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        canonical: &Atom,
        value: &Atom,
        note: Option<String>,
        compute_s: f64,
        deps: Option<&delta_store::Deps>,
    ) -> Result<()> {
        use delta_store::S;
        if !backing.writable() {
            delta_store::add(&S.not_stored_ro, 1);
            return Ok(());
        }
        let record = match self.store_record(
            store_key, identity, graph, support, canonical, note, compute_s, deps,
        ) {
            Ok(record) => record,
            Err("unknown") => {
                delta_store::add(&S.not_stored_unknown, 1);
                return Ok(());
            }
            Err(_) => {
                delta_store::add(&S.not_stored_export, 1);
                return Ok(());
            }
        };
        match backing.put(&record) {
            delta_store::PutResult::Inserted => {
                delta_store::add(&S.stored, 1);
                delta_store::count_loops(store_key.loops, false, true);
                Ok(())
            }
            delta_store::PutResult::NotWritten(reason) => {
                if S.not_written.fetch_add(1, Ordering::Relaxed) < 5 {
                    eprintln!("DELTA_STORE_NOT_WRITTEN {reason}");
                }
                Ok(())
            }
            delta_store::PutResult::Existing(row) => {
                if row.value_sha == record.value_sha {
                    delta_store::add(&S.existing_identical, 1);
                    backing.event(delta_store::StoreEvent::Confirm {
                        key_sha: store_key.key_sha,
                        tier: row.tier,
                        identical: true,
                    });
                    return Ok(());
                }
                let Ok(stored) = delta_store::import(&row) else {
                    // an unreadable row is never served (a miss for readers)
                    delta_store::add(&S.reject_import, 1);
                    return Ok(());
                };
                match identity.verify_agree(graph, support, model, &stored, value)? {
                    (_, true) => {
                        delta_store::add(&S.existing_agree, 1);
                        backing.event(delta_store::StoreEvent::Confirm {
                            key_sha: store_key.key_sha,
                            tier: row.tier,
                            identical: false,
                        });
                        Ok(())
                    }
                    (_, false) => {
                        delta_store::add(&S.disputes, 1);
                        let message = format!(
                            "DELTA_STORE_DISPUTE key_sha={} graph={} support={} loops={} stored_build={} stored_run_id={}",
                            delta_store::hex(&store_key.key_sha),
                            graph.name,
                            support.string_label(),
                            store_key.loops,
                            row.writer_build,
                            row.writer_run_id
                        );
                        eprintln!("{message}");
                        backing.event(delta_store::StoreEvent::Dispute {
                            kind: "conflict",
                            key_sha: store_key.key_sha,
                            tier: row.tier,
                            stored_value_sha: row.value_sha,
                            record: Box::new(record),
                            quarantine: true,
                        });
                        Err(eyre!(message))
                    }
                }
            }
        }
    }

    fn compute(
        &self,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<Atom> {
        let _stage = StageSpan::new("gamma.delta.compute");
        let classification_stage = StageSpan::new("gamma.delta.classify_spinneys");
        let catalogue =
            graph.classified_spinneys(support, model, settings, &graph.loop_momentum_basis);
        ensure!(
            catalogue
                .iter()
                .all(|spinney| spinney.renormalization_scheme == ApproximationType::PolePart),
            "complete UV-counterterm reuse requires pole-part subtraction on every support"
        );
        let overall = match catalogue.iter().find(|spinney| spinney.filter() == support) {
            Some(overall) => overall,
            None => {
                // The complete MS map of an overall convergent argument is zero.
                // Cache that answer too: checking the degree itself builds and
                // rescales its dressed integrand, so it belongs on the miss path.
                // A divergent miss already has its degree in the catalogue.
                let degree = graph.compute_dod(support);
                ensure!(
                    graph.approximation_scheme(support, model, settings, degree)
                        == ApproximationType::PolePart,
                    "complete UV-counterterm reuse requires pole-part subtraction on every support"
                );
                if degree < 0 {
                    return Ok(Atom::Zero);
                }
                return Err(eyre!(
                    "the UV argument has no overall spinney with a compatible momentum basis"
                ));
            }
        };
        drop(classification_stage);
        let loops = graph.n_loops(support);
        let mut integration_settings = settings.vakint.true_settings();
        let ctx = UVCtx::new(graph, model, settings);
        let current = ForestNode {
            spinney: overall,
            topo_order: 0,
        };
        let empty = Spinney::empty(graph);
        let mut result = Atom::Zero;
        // Bogoliubov's preparation sums products of *complete* connected Δs.
        // Replacing those products by their unfolded forests gives the existing
        // Wood evaluator exactly; keeping them aggregated permits reuse before
        // either the child Taylor expansion or its vacuum reduction is built.
        for spinney in std::iter::once(&empty).chain(
            catalogue
                .iter()
                .filter(|spinney| spinney.filter() != support && !spinney.filter().is_empty()),
        ) {
            let children_stage = StageSpan::new("gamma.delta.child_products");
            let mut poles = Atom::num(1);
            for component in graph.connected_components(spinney.filter()) {
                ensure!(
                    graph.n_loops(&component) < loops,
                    "a proper UV child must consume fewer loops than its parent"
                );
                poles *= self.counterterm(graph, &component, model, vakint, settings)?;
                if poles.is_zero() {
                    break;
                }
            }
            drop(children_stage);
            if poles.is_zero() {
                continue;
            }
            let consumed_loops = graph.n_loops(spinney.filter());
            let input = Full4dCts::from_counterterm_product(poles, consumed_loops);
            let given = ForestNode {
                spinney,
                topo_order: 0,
            };
            let local = local_4d::uv_limit(&input, &ctx, &current, &given, &current, &given)?;
            let remaining_loops = loops - consumed_loops;
            // Vakint truncates the complete quotient, including child poles,
            // through N - remaining_loops - 1. Request only its poles; the
            // master valuation planner supplies any finite or positive master
            // orders promoted by insertion or reduction-coefficient poles.
            // Preserve the existing depth for a zero-loop boundary case.
            integration_settings.number_of_terms_in_epsilon_expansion = if remaining_loops > 0 {
                remaining_loops as i64
            } else {
                loops as i64 + 1
            };
            let integrated = Integrated::new(vakint, &integration_settings);
            result += integrated
                .run(&local, &ctx, &current, &given, &current, &given)?
                .physical_pole_atom();
        }
        // Store an exact zero literally: parents then skip every spinney that
        // contains this support, and R* drops the preparation terms using it.
        Ok(literal_zero_if_vanishing(result))
    }
}

#[cfg(test)]
mod timing_tests {
    use std::{sync::mpsc, thread, time::Duration};

    use super::*;

    #[test]
    fn cache_timing_counts_single_flight_ownership_and_reuse() {
        assert!(CountertermCache::default().stats().timing.is_none());
        let timings = CacheTimings::new(true);
        let entry = OnceLock::new();
        let computations = AtomicU64::new(0);
        thread::scope(|scope| {
            let timings = &timings;
            let entry = &entry;
            let computations = &computations;
            let (started_tx, started_rx) = mpsc::channel();
            let (release_tx, release_rx) = mpsc::channel();
            let owner = scope.spawn(move || {
                let mut timing = timings.start();
                *timing.resolve(entry, |timing| {
                    computations.fetch_add(1, Ordering::Relaxed);
                    started_tx.send(()).unwrap();
                    release_rx.recv().unwrap();
                    timing.lap(CacheTimingPhase::Initializer);
                    42
                })
            });
            started_rx.recv().unwrap();
            let (waiting_tx, waiting_rx) = mpsc::channel();
            let (result_tx, result_rx) = mpsc::channel();
            let waiter = scope.spawn(move || {
                waiting_tx.send(()).unwrap();
                let result = {
                    let mut timing = timings.start();
                    *timing.resolve(entry, |_| panic!("one key must have one initializer"))
                };
                result_tx.send(result).unwrap();
            });
            waiting_rx.recv().unwrap();
            // A competing request cannot return before the owner publishes.
            assert!(matches!(
                result_rx.recv_timeout(Duration::from_millis(20)),
                Err(mpsc::RecvTimeoutError::Timeout)
            ));
            release_tx.send(()).unwrap();
            assert_eq!(owner.join().unwrap(), 42);
            assert_eq!(result_rx.recv().unwrap(), 42);
            waiter.join().unwrap();
        });
        {
            let mut timing = timings.start();
            assert_eq!(
                *timing.resolve(&entry, |_| panic!("ready hit recomputed")),
                42
            );
        }
        let stats = timings.stats().unwrap();
        assert_eq!(computations.load(Ordering::Relaxed), 1);
        assert_eq!((stats.calls, stats.initializers, stats.reuses), (3, 1, 2));
        // Incompleteness is a snapshot, not a scheduler guarantee: a delayed
        // contender can observe a ready entry after its owner is released.
        assert!(stats.pending_reuses <= 1);
        assert!(stats.pending_reuse_ns <= stats.reuse_ns);
        assert!(stats.initializer_ns > 0);
        assert!(stats.reuse_ns > 0);
    }
}
