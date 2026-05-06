use std::{
    env,
    sync::{
        Mutex, OnceLock,
        atomic::{AtomicU64, Ordering},
    },
    time::Instant,
};

const REPORT_INTERVAL_MS: u64 = 5_000;

#[derive(Clone, Copy, Debug)]
#[repr(usize)]
pub(crate) enum Counter {
    ParseView,
    ParseMul,
    ParseAdd,
    ParsePow,
    ParseFun,
    MulFactor,
    AddTerm,
    ScalarMulAccum,
    ScalarAddAccum,
    ParseStructureAttempt,
    ParseStructureOk,
    ParseStructureErr,
    NetworkNMul,
    NetworkNAdd,
    GraphNMul,
    GraphNAdd,
    GraphJoin,
    ShiftScalars,
    ShiftTensors,
    DanglingScan,
    StoreExtend,
    FromScalar,
    FromTensor,
    LibraryTensor,
    ExecuteIteration,
    ExecuteFound,
    ExecuteNeg,
    ExecuteProduct,
    ExecuteSum,
    ExecuteFunction,
    ExecutePower,
    Splice,
    MergeOps,
}

impl Counter {
    const COUNT: usize = Counter::MergeOps as usize + 1;

    const fn label(self) -> &'static str {
        match self {
            Counter::ParseView => "parse.view",
            Counter::ParseMul => "parse.mul",
            Counter::ParseAdd => "parse.add",
            Counter::ParsePow => "parse.pow",
            Counter::ParseFun => "parse.fun",
            Counter::MulFactor => "parse.mul.factor",
            Counter::AddTerm => "parse.add.term",
            Counter::ScalarMulAccum => "parse.scalar_mul_accum",
            Counter::ScalarAddAccum => "parse.scalar_add_accum",
            Counter::ParseStructureAttempt => "parse.structure.attempt",
            Counter::ParseStructureOk => "parse.structure.ok",
            Counter::ParseStructureErr => "parse.structure.err",
            Counter::NetworkNMul => "network.n_mul",
            Counter::NetworkNAdd => "network.n_add",
            Counter::GraphNMul => "graph.n_mul",
            Counter::GraphNAdd => "graph.n_add",
            Counter::GraphJoin => "graph.join",
            Counter::ShiftScalars => "graph.shift_scalars",
            Counter::ShiftTensors => "graph.shift_tensors",
            Counter::DanglingScan => "graph.dangling_scan",
            Counter::StoreExtend => "store.extend",
            Counter::FromScalar => "network.from_scalar",
            Counter::FromTensor => "network.from_tensor",
            Counter::LibraryTensor => "network.library_tensor",
            Counter::ExecuteIteration => "execute.iteration",
            Counter::ExecuteFound => "execute.found",
            Counter::ExecuteNeg => "execute.op.neg",
            Counter::ExecuteProduct => "execute.op.product",
            Counter::ExecuteSum => "execute.op.sum",
            Counter::ExecuteFunction => "execute.op.function",
            Counter::ExecutePower => "execute.op.power",
            Counter::Splice => "execute.splice",
            Counter::MergeOps => "execute.merge_ops",
        }
    }
}

#[derive(Clone, Copy, Debug)]
#[repr(usize)]
pub(crate) enum Timer {
    ParseMul,
    ParseAdd,
    ParsePow,
    ParseFun,
    ParseStructure,
    ScalarMulAccum,
    ScalarAddAccum,
    NetworkNMul,
    NetworkNMulStorePrep,
    NetworkNMulGraph,
    NetworkNMulDangling,
    NetworkNAdd,
    NetworkNAddStorePrep,
    NetworkNAddGraph,
    GraphNMul,
    GraphNAdd,
    GraphJoin,
    ShiftScalars,
    ShiftTensors,
    DanglingScan,
    SyncOrder,
    BuildGraph,
    StoreExtend,
    FromScalar,
    FromTensor,
    LibraryTensor,
    ExecuteFindReady,
    ExecuteCacheRoots,
    ExecuteReadyBatch,
    ExecuteReadyPreorder,
    ExecuteReadyParts,
    ExecuteReadySubgraph,
    ExecuteReadyIntersection,
    ExecuteReadyUnion,
    ExecutePlanOps,
    ExecuteOp,
    ExecuteNeg,
    ExecuteProduct,
    ExecuteSum,
    ExecuteFunction,
    ExecutePower,
    Splice,
    GraphExtract,
    IdentifyNodes,
    MergeOps,
    ExprTree,
}

impl Timer {
    const COUNT: usize = Timer::ExprTree as usize + 1;

    const fn label(self) -> &'static str {
        match self {
            Timer::ParseMul => "parse.mul",
            Timer::ParseAdd => "parse.add",
            Timer::ParsePow => "parse.pow",
            Timer::ParseFun => "parse.fun",
            Timer::ParseStructure => "parse.structure",
            Timer::ScalarMulAccum => "parse.scalar_mul_accum",
            Timer::ScalarAddAccum => "parse.scalar_add_accum",
            Timer::NetworkNMul => "network.n_mul",
            Timer::NetworkNMulStorePrep => "network.n_mul.store_prep",
            Timer::NetworkNMulGraph => "network.n_mul.graph",
            Timer::NetworkNMulDangling => "network.n_mul.dangling",
            Timer::NetworkNAdd => "network.n_add",
            Timer::NetworkNAddStorePrep => "network.n_add.store_prep",
            Timer::NetworkNAddGraph => "network.n_add.graph",
            Timer::GraphNMul => "graph.n_mul",
            Timer::GraphNAdd => "graph.n_add",
            Timer::GraphJoin => "graph.join",
            Timer::ShiftScalars => "graph.shift_scalars",
            Timer::ShiftTensors => "graph.shift_tensors",
            Timer::DanglingScan => "graph.dangling_scan",
            Timer::SyncOrder => "graph.sync_order",
            Timer::BuildGraph => "graph.build",
            Timer::StoreExtend => "store.extend",
            Timer::FromScalar => "network.from_scalar",
            Timer::FromTensor => "network.from_tensor",
            Timer::LibraryTensor => "network.library_tensor",
            Timer::ExecuteFindReady => "execute.find_ready",
            Timer::ExecuteCacheRoots => "execute.cache_roots",
            Timer::ExecuteReadyBatch => "execute.ready_batch",
            Timer::ExecuteReadyPreorder => "execute.ready.preorder",
            Timer::ExecuteReadyParts => "execute.ready.parts",
            Timer::ExecuteReadySubgraph => "execute.ready.subgraph",
            Timer::ExecuteReadyIntersection => "execute.ready.intersection",
            Timer::ExecuteReadyUnion => "execute.ready.union",
            Timer::ExecutePlanOps => "execute.plan_ops",
            Timer::ExecuteOp => "execute.op",
            Timer::ExecuteNeg => "execute.op.neg",
            Timer::ExecuteProduct => "execute.op.product",
            Timer::ExecuteSum => "execute.op.sum",
            Timer::ExecuteFunction => "execute.op.function",
            Timer::ExecutePower => "execute.op.power",
            Timer::Splice => "execute.splice",
            Timer::GraphExtract => "graph.extract",
            Timer::IdentifyNodes => "graph.identify_nodes",
            Timer::MergeOps => "execute.merge_ops",
            Timer::ExprTree => "graph.expr_tree",
        }
    }
}

struct Profile {
    counters: Vec<AtomicU64>,
    timer_counts: Vec<AtomicU64>,
    timer_nanos: Vec<AtomicU64>,
    start: Mutex<Instant>,
    last_report_ms: AtomicU64,
}

impl Profile {
    fn new() -> Self {
        Self {
            counters: (0..Counter::COUNT).map(|_| AtomicU64::new(0)).collect(),
            timer_counts: (0..Timer::COUNT).map(|_| AtomicU64::new(0)).collect(),
            timer_nanos: (0..Timer::COUNT).map(|_| AtomicU64::new(0)).collect(),
            start: Mutex::new(Instant::now()),
            last_report_ms: AtomicU64::new(0),
        }
    }

    fn elapsed_ms(&self) -> u64 {
        self.start.lock().unwrap().elapsed().as_millis() as u64
    }
}

pub(crate) struct Span {
    timer: Timer,
    start: Instant,
}

impl Drop for Span {
    fn drop(&mut self) {
        let profile = profile();
        let idx = self.timer as usize;
        profile.timer_counts[idx].fetch_add(1, Ordering::Relaxed);
        profile.timer_nanos[idx].fetch_add(
            self.start.elapsed().as_nanos().min(u128::from(u64::MAX)) as u64,
            Ordering::Relaxed,
        );
        maybe_report();
    }
}

static ENABLED: OnceLock<bool> = OnceLock::new();
static VERBOSE: OnceLock<bool> = OnceLock::new();
static ATOM_BYTES: OnceLock<bool> = OnceLock::new();
static PROFILE: OnceLock<Profile> = OnceLock::new();

fn profile() -> &'static Profile {
    PROFILE.get_or_init(Profile::new)
}

pub fn enabled() -> bool {
    *ENABLED.get_or_init(|| env::var_os("SPENSO_NETWORK_PROFILE").is_some())
}

pub fn verbose() -> bool {
    *VERBOSE.get_or_init(|| env::var_os("SPENSO_NETWORK_PROFILE_VERBOSE").is_some())
}

pub fn atom_bytes() -> bool {
    *ATOM_BYTES.get_or_init(|| env::var_os("SPENSO_NETWORK_PROFILE_ATOM_BYTES").is_some())
}

pub fn reset() {
    if !enabled() {
        return;
    }

    let profile = profile();
    for counter in &profile.counters {
        counter.store(0, Ordering::Relaxed);
    }
    for timer_count in &profile.timer_counts {
        timer_count.store(0, Ordering::Relaxed);
    }
    for timer_nanos in &profile.timer_nanos {
        timer_nanos.store(0, Ordering::Relaxed);
    }
    *profile.start.lock().unwrap() = Instant::now();
    profile.last_report_ms.store(0, Ordering::Relaxed);
}

pub(crate) fn bump(counter: Counter, amount: u64) {
    if !enabled() {
        return;
    }

    profile().counters[counter as usize].fetch_add(amount, Ordering::Relaxed);
    maybe_report();
}

pub(crate) fn span(timer: Timer) -> Option<Span> {
    enabled().then(|| Span {
        timer,
        start: Instant::now(),
    })
}

pub fn report(label: &str) {
    if enabled() {
        report_impl(label);
    }
}

fn maybe_report() {
    let profile = profile();
    let elapsed_ms = profile.elapsed_ms();
    let last = profile.last_report_ms.load(Ordering::Relaxed);

    if elapsed_ms < last + REPORT_INTERVAL_MS {
        return;
    }

    if profile
        .last_report_ms
        .compare_exchange(last, elapsed_ms, Ordering::Relaxed, Ordering::Relaxed)
        .is_ok()
    {
        report_impl("progress");
    }
}

fn report_impl(label: &str) {
    let profile = profile();
    let elapsed = profile.start.lock().unwrap().elapsed();
    eprintln!("spenso_profile {label}: elapsed={elapsed:.3?}");

    let mut timers = (0..Timer::COUNT)
        .filter_map(|idx| {
            let count = profile.timer_counts[idx].load(Ordering::Relaxed);
            let nanos = profile.timer_nanos[idx].load(Ordering::Relaxed);
            (count > 0).then_some((idx, count, nanos))
        })
        .collect::<Vec<_>>();
    timers.sort_by_key(|(_, _, nanos)| std::cmp::Reverse(*nanos));

    for (idx, count, nanos) in timers.iter().take(16) {
        let total_ms = *nanos as f64 / 1_000_000.0;
        let avg_us = *nanos as f64 / *count as f64 / 1_000.0;
        eprintln!(
            "spenso_profile timer {:<28} count={:<12} total_ms={:<12.3} avg_us={:.3}",
            timer_label(*idx),
            count,
            total_ms,
            avg_us,
        );
    }

    let counters = (0..Counter::COUNT)
        .filter_map(|idx| {
            let count = profile.counters[idx].load(Ordering::Relaxed);
            (count > 0).then_some((counter_label(idx), count))
        })
        .collect::<Vec<_>>();

    for (label, count) in counters {
        eprintln!("spenso_profile counter {label:<28} count={count}");
    }
}

fn timer_label(idx: usize) -> &'static str {
    match idx {
        idx if idx == Timer::ParseMul as usize => Timer::ParseMul.label(),
        idx if idx == Timer::ParseAdd as usize => Timer::ParseAdd.label(),
        idx if idx == Timer::ParsePow as usize => Timer::ParsePow.label(),
        idx if idx == Timer::ParseFun as usize => Timer::ParseFun.label(),
        idx if idx == Timer::ParseStructure as usize => Timer::ParseStructure.label(),
        idx if idx == Timer::ScalarMulAccum as usize => Timer::ScalarMulAccum.label(),
        idx if idx == Timer::ScalarAddAccum as usize => Timer::ScalarAddAccum.label(),
        idx if idx == Timer::NetworkNMul as usize => Timer::NetworkNMul.label(),
        idx if idx == Timer::NetworkNMulStorePrep as usize => Timer::NetworkNMulStorePrep.label(),
        idx if idx == Timer::NetworkNMulGraph as usize => Timer::NetworkNMulGraph.label(),
        idx if idx == Timer::NetworkNMulDangling as usize => Timer::NetworkNMulDangling.label(),
        idx if idx == Timer::NetworkNAdd as usize => Timer::NetworkNAdd.label(),
        idx if idx == Timer::NetworkNAddStorePrep as usize => Timer::NetworkNAddStorePrep.label(),
        idx if idx == Timer::NetworkNAddGraph as usize => Timer::NetworkNAddGraph.label(),
        idx if idx == Timer::GraphNMul as usize => Timer::GraphNMul.label(),
        idx if idx == Timer::GraphNAdd as usize => Timer::GraphNAdd.label(),
        idx if idx == Timer::GraphJoin as usize => Timer::GraphJoin.label(),
        idx if idx == Timer::ShiftScalars as usize => Timer::ShiftScalars.label(),
        idx if idx == Timer::ShiftTensors as usize => Timer::ShiftTensors.label(),
        idx if idx == Timer::DanglingScan as usize => Timer::DanglingScan.label(),
        idx if idx == Timer::SyncOrder as usize => Timer::SyncOrder.label(),
        idx if idx == Timer::BuildGraph as usize => Timer::BuildGraph.label(),
        idx if idx == Timer::StoreExtend as usize => Timer::StoreExtend.label(),
        idx if idx == Timer::FromScalar as usize => Timer::FromScalar.label(),
        idx if idx == Timer::FromTensor as usize => Timer::FromTensor.label(),
        idx if idx == Timer::LibraryTensor as usize => Timer::LibraryTensor.label(),
        idx if idx == Timer::ExecuteFindReady as usize => Timer::ExecuteFindReady.label(),
        idx if idx == Timer::ExecuteCacheRoots as usize => Timer::ExecuteCacheRoots.label(),
        idx if idx == Timer::ExecuteReadyBatch as usize => Timer::ExecuteReadyBatch.label(),
        idx if idx == Timer::ExecuteReadyPreorder as usize => Timer::ExecuteReadyPreorder.label(),
        idx if idx == Timer::ExecuteReadyParts as usize => Timer::ExecuteReadyParts.label(),
        idx if idx == Timer::ExecuteReadySubgraph as usize => Timer::ExecuteReadySubgraph.label(),
        idx if idx == Timer::ExecuteReadyIntersection as usize => {
            Timer::ExecuteReadyIntersection.label()
        }
        idx if idx == Timer::ExecuteReadyUnion as usize => Timer::ExecuteReadyUnion.label(),
        idx if idx == Timer::ExecutePlanOps as usize => Timer::ExecutePlanOps.label(),
        idx if idx == Timer::ExecuteOp as usize => Timer::ExecuteOp.label(),
        idx if idx == Timer::ExecuteNeg as usize => Timer::ExecuteNeg.label(),
        idx if idx == Timer::ExecuteProduct as usize => Timer::ExecuteProduct.label(),
        idx if idx == Timer::ExecuteSum as usize => Timer::ExecuteSum.label(),
        idx if idx == Timer::ExecuteFunction as usize => Timer::ExecuteFunction.label(),
        idx if idx == Timer::ExecutePower as usize => Timer::ExecutePower.label(),
        idx if idx == Timer::Splice as usize => Timer::Splice.label(),
        idx if idx == Timer::GraphExtract as usize => Timer::GraphExtract.label(),
        idx if idx == Timer::IdentifyNodes as usize => Timer::IdentifyNodes.label(),
        idx if idx == Timer::MergeOps as usize => Timer::MergeOps.label(),
        idx if idx == Timer::ExprTree as usize => Timer::ExprTree.label(),
        _ => "unknown",
    }
}

fn counter_label(idx: usize) -> &'static str {
    match idx {
        idx if idx == Counter::ParseView as usize => Counter::ParseView.label(),
        idx if idx == Counter::ParseMul as usize => Counter::ParseMul.label(),
        idx if idx == Counter::ParseAdd as usize => Counter::ParseAdd.label(),
        idx if idx == Counter::ParsePow as usize => Counter::ParsePow.label(),
        idx if idx == Counter::ParseFun as usize => Counter::ParseFun.label(),
        idx if idx == Counter::MulFactor as usize => Counter::MulFactor.label(),
        idx if idx == Counter::AddTerm as usize => Counter::AddTerm.label(),
        idx if idx == Counter::ScalarMulAccum as usize => Counter::ScalarMulAccum.label(),
        idx if idx == Counter::ScalarAddAccum as usize => Counter::ScalarAddAccum.label(),
        idx if idx == Counter::ParseStructureAttempt as usize => {
            Counter::ParseStructureAttempt.label()
        }
        idx if idx == Counter::ParseStructureOk as usize => Counter::ParseStructureOk.label(),
        idx if idx == Counter::ParseStructureErr as usize => Counter::ParseStructureErr.label(),
        idx if idx == Counter::NetworkNMul as usize => Counter::NetworkNMul.label(),
        idx if idx == Counter::NetworkNAdd as usize => Counter::NetworkNAdd.label(),
        idx if idx == Counter::GraphNMul as usize => Counter::GraphNMul.label(),
        idx if idx == Counter::GraphNAdd as usize => Counter::GraphNAdd.label(),
        idx if idx == Counter::GraphJoin as usize => Counter::GraphJoin.label(),
        idx if idx == Counter::ShiftScalars as usize => Counter::ShiftScalars.label(),
        idx if idx == Counter::ShiftTensors as usize => Counter::ShiftTensors.label(),
        idx if idx == Counter::DanglingScan as usize => Counter::DanglingScan.label(),
        idx if idx == Counter::StoreExtend as usize => Counter::StoreExtend.label(),
        idx if idx == Counter::FromScalar as usize => Counter::FromScalar.label(),
        idx if idx == Counter::FromTensor as usize => Counter::FromTensor.label(),
        idx if idx == Counter::LibraryTensor as usize => Counter::LibraryTensor.label(),
        idx if idx == Counter::ExecuteIteration as usize => Counter::ExecuteIteration.label(),
        idx if idx == Counter::ExecuteFound as usize => Counter::ExecuteFound.label(),
        idx if idx == Counter::ExecuteNeg as usize => Counter::ExecuteNeg.label(),
        idx if idx == Counter::ExecuteProduct as usize => Counter::ExecuteProduct.label(),
        idx if idx == Counter::ExecuteSum as usize => Counter::ExecuteSum.label(),
        idx if idx == Counter::ExecuteFunction as usize => Counter::ExecuteFunction.label(),
        idx if idx == Counter::ExecutePower as usize => Counter::ExecutePower.label(),
        idx if idx == Counter::Splice as usize => Counter::Splice.label(),
        idx if idx == Counter::MergeOps as usize => Counter::MergeOps.label(),
        _ => "unknown",
    }
}
