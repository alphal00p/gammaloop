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
        counterterm_identity::{CountertermIdentity, CountertermKey},
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
/// [`crate::uv::rstar::diagnostic_failure`], which no fallback may decline
/// (a decorated quotient's R* runs inside the child guard of an IR hard
/// part), stays one too: `rstar::is_diagnostic_failure` reads through the
/// replay), its text unchanged.
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

#[derive(Clone, PartialEq, Eq, Hash)]
struct Key {
    graph: CountertermKey,
    settings: String,
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
            let key = Key {
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
        self.lookup(
            graph,
            support,
            model,
            serde_json::to_string(settings)?,
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
        self.lookup(graph, support, model, settings_key, compute)
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
        compute: impl FnOnce() -> Result<(Atom, CacheNote)>,
    ) -> Result<(Atom, CacheLookup)> {
        self.lookup_noted(
            graph,
            support,
            model,
            settings_key,
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
        compute: impl FnOnce() -> Result<Atom>,
    ) -> Result<Atom> {
        self.lookup_noted(
            graph,
            support,
            model,
            settings_key,
            || compute().map(|value| (value, None)),
            false,
        )
        .map(|(value, _)| value)
    }

    /// The lookup; `noted`: the fill's note is stored with its entry and a
    /// hit returns it (no note is kept or read otherwise).
    fn lookup_noted(
        &self,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        settings_key: String,
        compute: impl FnOnce() -> Result<(Atom, Option<CacheNote>)>,
        noted: bool,
    ) -> Result<(Atom, CacheLookup)> {
        let _stage = StageSpan::new("gamma.delta.lookup");
        let mut timing = self.timings.start();
        let identity_stage = StageSpan::new("gamma.delta.identity");
        let identity = CountertermIdentity::new(graph, support, model)?;
        let key = Key {
            graph: identity.key.clone(),
            settings: settings_key,
        };
        let noted_key = noted.then(|| key.clone());
        timing.lap(CacheTimingPhase::Key);
        drop(identity_stage);
        let (entry, existing) = {
            let _stage = StageSpan::new("gamma.delta.cache_map");
            let mut entries = self.entries.lock().expect("counterterm cache poisoned");
            timing.lap(CacheTimingPhase::MapLockWait);
            let existing = entries.contains_key(&key);
            let entry = entries.entry(key).or_default().clone();
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
        let mut filled = false;
        let result = timing.resolve(&entry, |timing| {
            filled = true;
            self.misses.fetch_add(1, Ordering::Relaxed);
            let result = compute();
            timing.lap(CacheTimingPhase::Initializer);
            let _stage = StageSpan::new("gamma.delta.cache_store");
            let result = result
                .and_then(|(value, note)| {
                    if let (Some(key), Some(note)) = (&noted_key, note) {
                        self.notes
                            .lock()
                            .expect("counterterm cache notes poisoned")
                            .insert(key.clone(), note);
                    }
                    identity.canonicalize_delta(&value)
                })
                .inspect(|_| {
                    self.stores.fetch_add(1, Ordering::Relaxed);
                })
                .map_err(|error| CachedError(Arc::new(error)));
            timing.lap(CacheTimingPhase::Store);
            result
        });
        drop(resolve_stage);
        let canonical = result.as_ref().map_err(CachedError::replay)?;
        debug_tags!(#uv, #integrated, #cache;
            stage = "uv_counterterm_cache_lookup",
            hit = existing,
            loops = graph.n_loops(support),
            support = %support.string_label(),
            "Complete UV counterterm cache lookup"
        );
        timing.restart();
        let _stage = StageSpan::new("gamma.delta.restore_boundary");
        let result = identity.restore_delta(canonical);
        timing.lap(CacheTimingPhase::Restore);
        let lookup = if filled {
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
