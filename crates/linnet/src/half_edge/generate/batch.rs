//! Binding transport for typed predicates. Logical searches suspend before a verdict is
//! needed; the calling thread dispatches a whole wave after native workers have returned.
//! This is deliberately separate from the public scalar predicate signatures.

use std::{
    collections::VecDeque,
    sync::{
        atomic::{AtomicBool, AtomicU64, Ordering},
        Arc, OnceLock,
    },
};

use super::{
    decorate::{Decoration, Decorator, Edge},
    model::VertexState,
    output::SkeletonData,
    suspend::{Step, Work},
    topology::{Search, Topology},
    view::{State, View},
    Compiled, Counts, Diagram, DiagramView, DressingView, GenerationError, HedgeGraph, Limits,
    OutputTables, Plan, Progress, Run, RunContext, Scratch, SkeletonSnapshot, SkeletonView,
};

/// Which scalar predicates are serviced by the calling-thread transport.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct HookMask {
    pub skeletons: bool,
    pub dressings: bool,
    pub select: bool,
}

impl HookMask {
    pub(super) fn any(self) -> bool {
        self.skeletons || self.dressings || self.select
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Hook {
    Skeleton,
    Dressing,
    Select,
}

/// An execution-scoped continuation. No pointer or recycled slot address crosses the binding.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Ticket {
    run: u64,
    serial: u64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Reply {
    ticket: Ticket,
    hook: Hook,
    keep: bool,
}

/// Compact, owned callback state. It is released when its scoped binding view expires.
pub struct Request<F, L> {
    ticket: Ticket,
    hook: Hook,
    snapshot: Snapshot<F, L>,
}

impl<F, L> Request<F, L> {
    /// Sorted original keys used by the compact state's field and label ranks.
    pub fn fields(&self) -> &[F] {
        &self.snapshot.compiled.tables.fields
    }
    pub fn labels(&self) -> &[L] {
        &self.snapshot.compiled.tables.labels
    }

    /// Translate keys to their existing dense ranks while preserving the continuation.
    /// The model and its automata stay shared; only small output tables are cached once.
    /// Graph caches have key-dependent payloads: ranked requests cache their own repeated reads,
    /// but cannot transfer that graph back to an output with the original key types.
    pub fn into_ranked(self) -> Request<u32, u32> {
        let compiled = self.snapshot.compiled.ranked.get_or_init(|| {
            let tables = &self.snapshot.compiled.tables;
            Arc::new(Compiled {
                model: Arc::clone(&self.snapshot.compiled.model),
                tables: Arc::new(OutputTables {
                    fields: (0..tables.fields.len() as u32).collect(),
                    labels: (0..tables.labels.len() as u32).collect(),
                    packing: tables.packing.clone(),
                    signatures: Arc::clone(&tables.signatures),
                    weights: Arc::clone(&tables.weights),
                    legs: tables.legs.clone(),
                }),
                ranked: OnceLock::new(),
            })
        });
        Request {
            ticket: self.ticket,
            hook: self.hook,
            snapshot: Snapshot {
                compiled: Arc::clone(compiled),
                loops: self.snapshot.loops,
                topology: self.snapshot.topology,
                decoration: self.snapshot.decoration,
                complete: self
                    .snapshot
                    .complete
                    .as_ref()
                    .map(|_| Arc::new(OnceLock::new())),
            },
        }
    }

    pub fn ticket(&self) -> Ticket {
        self.ticket
    }
    pub fn hook(&self) -> Hook {
        self.hook
    }
    pub fn reply(&self, keep: bool) -> Reply {
        Reply {
            ticket: self.ticket,
            hook: self.hook,
            keep,
        }
    }
    pub fn skeleton(&self) -> Option<SkeletonView<'_, F, L>> {
        (self.hook == Hook::Skeleton).then(|| SkeletonView(self.snapshot.view()))
    }
    pub fn dressing(&self) -> Option<DressingView<'_, F, L>> {
        (self.hook == Hook::Dressing).then(|| DressingView(self.snapshot.view()))
    }
    pub fn diagram(&self) -> Option<DiagramView<'_, F, L>> {
        (self.hook == Hook::Select).then(|| DiagramView(self.snapshot.view()))
    }
}

struct Snapshot<F, L> {
    compiled: Arc<Compiled<F, L>>,
    loops: usize,
    topology: Topology,
    decoration: Option<Assignment>,
    complete: Option<Arc<OnceLock<HedgeGraph<F, usize, Option<L>>>>>,
}

struct Assignment {
    edges: Box<[Edge]>,
    labels: Box<[u8]>,
    states: Box<[VertexState]>,
    rules: Box<[u32]>,
    automorphisms: Option<u128>,
    vetoed: bool,
}

impl<F, L> Snapshot<F, L> {
    fn new(compiled: &Arc<Compiled<F, L>>, loops: usize, state: State<'_>) -> Self {
        let (topology, decoration) = match state {
            State::Topology(t) => (t, None),
            State::Decoration(d) => (
                d.topology,
                Some(Assignment {
                    edges: d.edges.into(),
                    labels: d.labels.into(),
                    states: d.states.into(),
                    rules: d.rules.into(),
                    automorphisms: d.automorphism_count,
                    vetoed: d.vetoed,
                }),
            ),
        };
        Self {
            compiled: Arc::clone(compiled),
            loops,
            topology: topology.snapshot(),
            decoration,
            complete: None,
        }
    }

    fn view(&self) -> View<'_, F, L> {
        let state = match &self.decoration {
            None => State::Topology(&self.topology),
            Some(d) => State::Decoration(Decoration {
                topology: &self.topology,
                edges: &d.edges,
                labels: &d.labels,
                states: &d.states,
                rules: &d.rules,
                automorphism_count: d.automorphisms,
                vetoed: d.vetoed,
            }),
        };
        View::new(&self.compiled, self.loops, state).with_complete_cache(self.complete.as_deref())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Output {
    Diagrams,
    Skeletons,
    Count,
}

pub enum Event<F, L> {
    Batch(Vec<Request<F, L>>),
    Diagram(Diagram<F, L>),
    Skeleton(SkeletonSnapshot<F, L>),
    Progress(Progress),
    Yield,
    Done,
}

/// Fixed logical capacity and work quantum keep dispatch order independent of worker count.
const WINDOW: usize = 64;
const QUANTUM: usize = 128;
static NEXT_RUN: AtomicU64 = AtomicU64::new(1);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Phase {
    Search,
    StartDressing,
    Decorate,
    Complete,
    Select,
    Emit,
    AwaitSkeleton,
    AwaitDressing,
    AwaitComplete,
    AwaitSelect,
    Done,
}

struct Cursor<F, L> {
    search: Option<Search>,
    decorator: Decorator,
    skeleton: Option<Arc<SkeletonData<F, L>>>,
    phase: Phase,
    /// Only the current complete assignment may share a demanded graph across hooks.
    complete: Option<Arc<OnceLock<HedgeGraph<F, usize, Option<L>>>>>,
}

enum Job<F, L> {
    Source(Search),
    Candidate(Topology),
    Cursor(Cursor<F, L>),
}

enum Item<F, L> {
    Diagram(Diagram<F, L>),
    Skeleton(SkeletonSnapshot<F, L>),
    Count(Option<u128>),
}

struct Slot<F, L> {
    id: u64,
    job: Job<F, L>,
    state: SlotState<F, L>,
}

struct SlotState<F, L> {
    pending: Option<(Hook, Snapshot<F, L>)>,
    waiting: bool,
    outputs: VecDeque<Item<F, L>>,
    topologies: usize,
    done: bool,
}

impl<F, L> Slot<F, L> {
    fn new(id: u64, job: Job<F, L>) -> Self {
        Self {
            id,
            job,
            state: SlotState {
                pending: None,
                waiting: false,
                outputs: VecDeque::new(),
                topologies: 0,
                done: false,
            },
        }
    }
}

/// A pull/resume execution. Each `poll` joins its native wave before returning; dropping this
/// owner never leaves workers waiting for Python or snapshots borrowing mutable scratch.
pub struct Driver<'a, F, L> {
    plan: &'a Plan<F, L>,
    run: &'a Run<'a>,
    output: Output,
    orders: VecDeque<(usize, Arc<Limits>)>,
    current: Option<(usize, Arc<Limits>)>,
    slots: VecDeque<Slot<F, L>>,
    pending: Vec<(Ticket, Hook, u64)>,
    generation: u64,
    serial: u64,
    next_slot: u64,
    stop: AtomicBool,
    counts: Counts,
    count_overflow: bool,
    progress: Progress,
    ready: Option<Item<F, L>>,
    alone: bool,
    finished: bool,
}

impl<F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> Plan<F, L> {
    pub fn batched<'a>(
        &'a self,
        run: &'a Run<'a>,
        output: Output,
    ) -> Result<Driver<'a, F, L>, GenerationError> {
        Ok(Driver::new(self, run, output))
    }

    pub fn batched_dressings<'a>(
        &'a self,
        skeleton: &SkeletonSnapshot<F, L>,
        run: &'a Run<'a>,
    ) -> Result<Driver<'a, F, L>, GenerationError> {
        if !Arc::ptr_eq(&self.compiled.tables, &skeleton.data.tables) {
            return Err(GenerationError::ForeignSkeleton);
        }
        let mut driver = Driver::new(self, run, Output::Diagrams);
        driver.orders.clear();
        driver.current = self
            .orders
            .iter()
            .find(|(l, _)| *l == skeleton.data.loops)
            .cloned();
        driver.progress.loops = skeleton.data.loops;
        let mut decorator = Decorator::default();
        decorator.start(
            &self.compiled.model,
            &skeleton.data.topology,
            &skeleton.data.generators,
            skeleton.data.node_group_order,
        );
        let mut slot = Slot::new(
            0,
            Job::Cursor(Cursor {
                search: None,
                decorator,
                skeleton: Some(Arc::clone(&skeleton.data)),
                phase: Phase::Decorate,
                complete: None,
            }),
        );
        slot.state.topologies = 1;
        driver.slots.push_back(slot);
        driver.next_slot = 1;
        Ok(driver)
    }
}

impl<'a, F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> Driver<'a, F, L> {
    fn new(plan: &'a Plan<F, L>, run: &'a Run<'a>, output: Output) -> Self {
        Self {
            plan,
            run,
            output,
            orders: plan.orders.iter().cloned().collect(),
            current: None,
            slots: VecDeque::new(),
            pending: Vec::new(),
            generation: NEXT_RUN.fetch_add(1, Ordering::Relaxed),
            serial: 0,
            next_slot: 0,
            stop: AtomicBool::new(false),
            counts: Counts::default(),
            count_overflow: false,
            progress: Progress::default(),
            ready: None,
            alone: true,
            finished: false,
        }
    }

    /// Exact totals. Diagram iteration remains available when an automorphism order is too
    /// large for the native certificate; asking for its reciprocal totals then reports that.
    pub fn counts(&self) -> Result<&Counts, GenerationError> {
        if self.count_overflow {
            return Err(GenerationError::AutomorphismOverflow);
        }
        Ok(&self.counts)
    }
    pub fn progress(&self) -> &Progress {
        &self.progress
    }
    pub fn aborted(&self) -> bool {
        self.stop.load(Ordering::Relaxed)
    }
    pub fn abort(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        self.pending.clear();
        self.slots.clear();
        self.ready = None;
        self.finished = true;
        self.counts.complete = false;
    }

    /// A progress veto stops future work but delivers the accepted item already counted at
    /// that notification, matching scalar iteration. Exceptions use `abort` instead.
    pub fn stop_after_progress(&mut self) {
        let ready = self.ready.take();
        self.abort();
        self.ready = ready;
    }

    /// Validate the entire reply before consuming any continuation. Stale, reordered, duplicate,
    /// wrong-hook or truncated batches cannot partially resume a run.
    pub fn resolve(&mut self, replies: &[Reply]) -> Result<(), GenerationError> {
        if self.aborted()
            || self.pending.is_empty()
            || replies.len() != self.pending.len()
            || replies
                .iter()
                .zip(&self.pending)
                .any(|(r, (ticket, hook, _))| r.ticket != *ticket || r.hook != *hook)
        {
            self.abort();
            return Err(GenerationError::InvalidReply);
        }
        let cx = RunContext::new(self.plan, self.run, &self.current, &self.stop);
        for (reply, (_, _, id)) in replies.iter().zip(self.pending.drain(..)) {
            let slot = self
                .slots
                .iter_mut()
                .find(|s| s.id == id)
                .expect("pending slot remains owned");
            match &mut slot.job {
                Job::Candidate(topology) => {
                    if topology.resolve_prune(reply.keep, self.run.check_heredity) {
                        let scratch = Scratch::take();
                        slot.job = Job::Cursor(Cursor::new(
                            &self.plan.compiled,
                            topology.clone(),
                            scratch,
                        ));
                    } else {
                        slot.state.done = true;
                    }
                }
                Job::Cursor(cursor) => cursor.resume(&cx, reply.keep, self.run.check_heredity),
                Job::Source(_) => unreachable!("sources transfer pending children before dispatch"),
            }
            slot.state.waiting = false;
        }
        Ok(())
    }

    pub fn poll(&mut self) -> Result<Event<F, L>, GenerationError> {
        if self.run.abort_check.as_ref().is_some_and(|check| check()) {
            self.abort();
        }
        if let Some(item) = self.ready.take() {
            return Ok(Self::event(item));
        }
        if self.aborted() || self.finished {
            return Ok(Event::Done);
        }
        if !self.pending.is_empty() {
            return Err(GenerationError::InvalidReply);
        }
        if let Some(event) = self.notify_topologies() {
            return Ok(event);
        }
        self.discard_finished();
        if self.slots.is_empty() {
            let Some((loops, limits)) = self.orders.pop_front() else {
                self.finished = true;
                self.counts.complete = true;
                return Ok(Event::Done);
            };
            self.progress.loops = loops;
            let root = Topology::empty(&self.plan.compiled.model, &limits);
            self.current = Some((loops, limits));
            let scratch = Scratch::take();
            let source = Search::child_source(&self.plan.compiled.model, root, scratch.search);
            self.slots
                .push_back(Slot::new(self.next_slot, Job::Source(source)));
            self.next_slot += 1;
        }
        self.split_one();
        self.fill_sources();
        let cx = RunContext::new(self.plan, self.run, &self.current, &self.stop);
        let compiled = &self.plan.compiled;
        let output = self.output;
        let ahead = if self.alone { 1 } else { super::AHEAD };
        // The first output is a notification boundary even for exact counting. Later slots may
        // supply hook requests, but do not emit before the entire earlier subtree.
        let first_cursor = self
            .slots
            .iter()
            .position(|s| !s.state.done && matches!(s.job, Job::Cursor(_)));
        let alone = self.alone;
        let step = |(index, slot): (usize, &mut Slot<F, L>)| {
            if slot.state.done
                || slot.state.waiting
                || slot.state.pending.is_some()
                || slot.state.outputs.len() >= ahead
                || (alone && first_cursor != Some(index))
            {
                return;
            }
            if let Job::Cursor(cursor) = &mut slot.job {
                let mut work = Work::bounded(QUANTUM);
                cursor.advance(&cx, compiled, output, &mut work, &mut slot.state);
            }
        };
        #[cfg(all(feature = "rayon", not(target_arch = "wasm32")))]
        if self.run.parallel {
            use rayon::prelude::*;
            self.slots
                .make_contiguous()
                .par_iter_mut()
                .enumerate()
                .for_each(step);
        } else {
            self.slots.iter_mut().enumerate().for_each(step);
        }
        #[cfg(not(all(feature = "rayon", not(target_arch = "wasm32"))))]
        self.slots.iter_mut().enumerate().for_each(step);
        if let Some(event) = self.notify_topologies() {
            return Ok(event);
        }
        let loops = self.progress.loops;
        if self.aborted() {
            self.abort();
            return Ok(Event::Done);
        }
        let mut requests = Vec::new();
        for slot in &mut self.slots {
            if let Some((hook, snapshot)) = slot.state.pending.take() {
                let ticket = Ticket {
                    run: self.generation,
                    serial: self.serial,
                };
                self.serial = self
                    .serial
                    .checked_add(1)
                    .expect("callback ticket space exhausted");
                self.pending.push((ticket, hook, slot.id));
                slot.state.waiting = true;
                requests.push(Request {
                    ticket,
                    hook,
                    snapshot,
                });
            }
        }
        if !requests.is_empty() {
            return Ok(Event::Batch(requests));
        }
        self.discard_finished();
        if let Some(item) = self
            .slots
            .front_mut()
            .and_then(|slot| slot.state.outputs.pop_front())
        {
            self.alone = false;
            let aut = match &item {
                Item::Diagram(d) => Some(d.automorphism_count),
                Item::Count(n) => Some(*n),
                Item::Skeleton(_) => None,
            };
            if let Some(aut) = aut {
                if let Err(error) = self.counts.diagram(loops, aut) {
                    self.count_overflow = true;
                    if self.output == Output::Count {
                        self.abort();
                        return Err(error);
                    }
                    self.counts.diagrams += 1u8;
                    self.counts.per_loop.entry(loops).or_default().diagrams += 1u8;
                }
                self.progress.diagrams += 1u8;
            }
            super::report(self.run.progress.as_deref(), &self.progress, &self.stop);
            if self.aborted() {
                self.abort();
                return Ok(Self::event(item));
            }
            if self.run.deferred_progress {
                self.ready = Some(item);
                return Ok(Event::Progress(self.progress.clone()));
            }
            return Ok(Self::event(item));
        }
        Ok(Event::Yield)
    }

    fn event(item: Item<F, L>) -> Event<F, L> {
        match item {
            Item::Diagram(d) => Event::Diagram(d),
            Item::Skeleton(s) => Event::Skeleton(s),
            Item::Count(_) => Event::Yield,
        }
    }

    fn notify_topologies(&mut self) -> Option<Event<F, L>> {
        let mut changed = false;
        for slot in &mut self.slots {
            let count = std::mem::take(&mut slot.state.topologies);
            if count == 0 {
                continue;
            }
            changed = true;
            self.progress.topologies += count;
            self.counts.skeletons += count;
            self.counts
                .per_loop
                .entry(self.progress.loops)
                .or_default()
                .skeletons += count;
        }
        if !changed {
            return None;
        }
        super::report(self.run.progress.as_deref(), &self.progress, &self.stop);
        if self.aborted() {
            self.abort();
            return Some(Event::Done);
        }
        self.run
            .deferred_progress
            .then(|| Event::Progress(self.progress.clone()))
    }

    fn discard_finished(&mut self) {
        self.slots
            .retain(|s| !s.state.done || !s.state.outputs.is_empty());
    }

    fn split_one(&mut self) {
        // A running subtree can donate its residual sibling enumerator without replaying the
        // accepted prefix. Sources and pending candidates consume the same capacity as cursors.
        if self.slots.len() < WINDOW {
            let cx = RunContext::new(self.plan, self.run, &self.current, &self.stop);
            for i in 0..self.slots.len() {
                let slot = &mut self.slots[i];
                if slot.state.waiting || slot.state.pending.is_some() || slot.state.done {
                    continue;
                }
                if let Job::Cursor(cursor) = &mut slot.job {
                    if let Some(search) = &mut cursor.search {
                        let mut scratch = Scratch::take();
                        if let Some(source) = search.split_source(&cx.search(), &mut scratch.search)
                        {
                            self.slots
                                .insert(i + 1, Slot::new(self.next_slot, Job::Source(source)));
                            self.next_slot += 1;
                            scratch.give();
                            break;
                        }
                        scratch.give();
                    }
                }
            }
        }
    }

    fn fill_sources(&mut self) {
        let cx = RunContext::new(self.plan, self.run, &self.current, &self.stop);
        let mut work = Work::bounded(QUANTUM);
        let mut index = 0;
        while index < self.slots.len() && self.slots.len() < WINDOW {
            let slot = &mut self.slots[index];
            if slot.state.done {
                index += 1;
                continue;
            }
            let Job::Source(source) = &mut slot.job else {
                index += 1;
                continue;
            };
            let deferred = self.plan.hooks.deferred.skeletons;
            let (topology, pending) = match source.poll(&cx.search(), &mut work, deferred) {
                Step::Pending => (source.take_pending_child(), true),
                Step::Ready => (source.topology().clone(), false),
                Step::Done => {
                    slot.state.done = true;
                    index += 1;
                    continue;
                }
                Step::Yield => break,
            };
            let mut topology = topology;
            let snapshot = pending
                .then(|| Snapshot::new(&self.plan.compiled, cx.loops, State::Topology(&topology)));
            let native_keeps = snapshot.as_ref().is_none_or(|s| {
                self.plan
                    .hooks
                    .skeletons
                    .as_ref()
                    .is_none_or(|hook| hook(&SkeletonView(s.view())))
            });
            if !native_keeps && !topology.resolve_prune(false, self.run.check_heredity) {
                continue;
            }
            let (job, pending) = if native_keeps && pending {
                (
                    Job::Candidate(topology),
                    snapshot.map(|s| (Hook::Skeleton, s)),
                )
            } else {
                (
                    Job::Cursor(Cursor::new(&self.plan.compiled, topology, Scratch::take())),
                    None,
                )
            };
            let mut child = Slot::new(self.next_slot, job);
            child.state.pending = pending;
            self.next_slot += 1;
            self.slots.insert(index, child);
            index += 1;
        }
    }
}

impl<F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> Cursor<F, L> {
    fn new(compiled: &Compiled<F, L>, topology: Topology, scratch: Scratch) -> Self {
        Self {
            search: Some(Search::new(&compiled.model, topology, scratch.search)),
            decorator: scratch.decorator,
            skeleton: None,
            phase: Phase::Search,
            complete: None,
        }
    }

    fn topology(&self) -> &Topology {
        self.search
            .as_ref()
            .map(Search::topology)
            .unwrap_or_else(|| {
                &self
                    .skeleton
                    .as_ref()
                    .expect("retained dressing skeleton")
                    .topology
            })
    }

    fn resume(&mut self, cx: &RunContext<'_, F, L>, keep: bool, check: bool) {
        match self.phase {
            Phase::AwaitSkeleton => {
                self.search.as_mut().unwrap().resume(&cx.search(), keep);
                self.phase = Phase::Search;
            }
            Phase::AwaitDressing => {
                self.decorator.resume(&cx.compiled.model, keep, check);
                self.phase = Phase::Decorate;
            }
            Phase::AwaitComplete => {
                let d = self.decorator.decoration(self.topology());
                assert!(
                    !(keep && d.vetoed),
                    "prune_dressings is not hereditary: it accepts a diagram but rejected an ancestor assignment"
                );
                self.phase = if keep { Phase::Select } else { Phase::Decorate };
            }
            Phase::AwaitSelect => self.phase = if keep { Phase::Emit } else { Phase::Decorate },
            _ => unreachable!("only a parked continuation accepts a verdict"),
        }
        if self.phase == Phase::Decorate {
            self.complete = None;
        }
    }

    fn advance(
        &mut self,
        cx: &RunContext<'_, F, L>,
        compiled: &Arc<Compiled<F, L>>,
        output: Output,
        work: &mut Work,
        state: &mut SlotState<F, L>,
    ) {
        loop {
            if cx.stopped() {
                self.phase = Phase::Done;
            }
            match self.phase {
                Phase::Done => {
                    self.complete = None;
                    state.done = true;
                    return;
                }
                Phase::Search => {
                    let search = self.search.as_mut().expect("search phase has a search");
                    match search.poll(&cx.search(), work, cx.hooks.deferred.skeletons) {
                        Step::Pending => {
                            let snapshot = Snapshot::new(
                                compiled,
                                cx.loops,
                                State::Topology(search.topology()),
                            );
                            self.phase = Phase::AwaitSkeleton;
                            if cx
                                .hooks
                                .skeletons
                                .as_ref()
                                .is_some_and(|hook| !hook(&SkeletonView(snapshot.view())))
                            {
                                self.resume(
                                    cx,
                                    false,
                                    cx.filter.as_ref().is_some_and(|p| p.check_heredity),
                                );
                                continue;
                            }
                            state.pending = Some((Hook::Skeleton, snapshot));
                            return;
                        }
                        Step::Yield => return,
                        Step::Done => {
                            self.phase = Phase::Done;
                            continue;
                        }
                        Step::Ready => {
                            state.topologies += 1;
                            if output != Output::Count {
                                self.skeleton =
                                    Some(Arc::new(SkeletonData::new(compiled, search, cx.loops)));
                            }
                            if output == Output::Skeletons {
                                state.outputs.push_back(Item::Skeleton(SkeletonSnapshot {
                                    data: Arc::clone(self.skeleton.as_ref().unwrap()),
                                }));
                                return;
                            }
                            self.phase = Phase::StartDressing;
                            return;
                        }
                    }
                }
                Phase::StartDressing => {
                    let search = self.search.as_ref().expect("a discovered topology");
                    self.decorator.start(
                        &compiled.model,
                        search.topology(),
                        search.generators(),
                        search.node_group_order(),
                    );
                    self.phase = Phase::Decorate;
                }
                Phase::Decorate => {
                    self.complete = None;
                    let topology = self
                        .search
                        .as_ref()
                        .map(Search::topology)
                        .unwrap_or_else(|| &self.skeleton.as_ref().unwrap().topology);
                    match self.decorator.poll(
                        &compiled.model,
                        topology,
                        cx.prune(),
                        cx.hooks.deferred.dressings,
                        work,
                    ) {
                        Step::Pending => {
                            let snapshot = Snapshot::new(
                                compiled,
                                cx.loops,
                                State::Decoration(self.decorator.pending(topology)),
                            );
                            self.phase = Phase::AwaitDressing;
                            if cx
                                .hooks
                                .dressings
                                .as_ref()
                                .is_some_and(|hook| !hook(&DressingView(snapshot.view())))
                            {
                                self.resume(
                                    cx,
                                    false,
                                    cx.filter.as_ref().is_some_and(|p| p.check_heredity),
                                );
                                continue;
                            }
                            state.pending = Some((Hook::Dressing, snapshot));
                            return;
                        }
                        Step::Yield => return,
                        Step::Done => {
                            self.phase = if self.search.is_some() {
                                Phase::Search
                            } else {
                                Phase::Done
                            }
                        }
                        Step::Ready => self.phase = Phase::Complete,
                    }
                }
                Phase::Complete | Phase::Select => {
                    if !work.spend() {
                        return;
                    }
                    let hook = if self.phase == Phase::Complete {
                        Hook::Dressing
                    } else {
                        Hook::Select
                    };
                    let (present, deferred) = match hook {
                        Hook::Dressing => {
                            (cx.hooks.dressings.is_some(), cx.hooks.deferred.dressings)
                        }
                        Hook::Select => (cx.hooks.select.is_some(), cx.hooks.deferred.select),
                        Hook::Skeleton => unreachable!(),
                    };
                    if present || deferred {
                        self.complete
                            .get_or_insert_with(|| Arc::new(OnceLock::new()));
                    }
                    let d = self.decorator.decoration(self.topology());
                    let view = View::new(compiled, cx.loops, State::Decoration(d))
                        .with_complete_cache(self.complete.as_deref());
                    let keep = match hook {
                        Hook::Dressing => cx
                            .hooks
                            .dressings
                            .as_ref()
                            .is_none_or(|h| h(&DressingView(view))),
                        Hook::Select => cx
                            .hooks
                            .select
                            .as_ref()
                            .is_none_or(|h| h(&DiagramView(view))),
                        Hook::Skeleton => unreachable!(),
                    };
                    self.phase = if hook == Hook::Dressing {
                        Phase::AwaitComplete
                    } else {
                        Phase::AwaitSelect
                    };
                    if keep && deferred {
                        let mut snapshot = Snapshot::new(
                            compiled,
                            cx.loops,
                            State::Decoration(self.decorator.decoration(self.topology())),
                        );
                        snapshot.complete = self.complete.clone();
                        state.pending = Some((hook, snapshot));
                        return;
                    }
                    self.resume(
                        cx,
                        keep,
                        cx.filter.as_ref().is_some_and(|p| p.check_heredity),
                    );
                }
                Phase::Emit => {
                    let graph = self.complete.take().and_then(|cache| {
                        if output == Output::Count {
                            return None;
                        }
                        match Arc::try_unwrap(cache) {
                            Ok(cache) => cache.into_inner(),
                            // Native consumers may retain a request after resolving it. Preserve
                            // that snapshot; clone only its already-demanded graph in this case.
                            Err(cache) => cache.get().cloned(),
                        }
                    });
                    let d = self.decorator.decoration(self.topology());
                    let item = if output == Output::Count {
                        Item::Count(d.automorphism_count)
                    } else {
                        Item::Diagram(
                            compiled.diagram(
                                d,
                                cx.loops,
                                View::new(compiled, cx.loops, State::Decoration(d))
                                    .with_complete_graph(graph),
                                Arc::clone(self.skeleton.as_ref().unwrap()),
                            ),
                        )
                    };
                    state.outputs.push_back(item);
                    self.phase = Phase::Decorate;
                    return;
                }
                Phase::AwaitSkeleton
                | Phase::AwaitDressing
                | Phase::AwaitComplete
                | Phase::AwaitSelect => return,
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::half_edge::generate::{Legs, Query, Stub, Theory};
    use crate::half_edge::NodeIndex;
    use std::collections::BTreeMap;
    use std::sync::Mutex;

    fn input(rules: usize) -> (Theory<u8>, Query<u8, u8>) {
        let mut theory = Theory::builder();
        for _ in 0..rules {
            theory = theory.vertex(vec![Stub::undirected(0u8); 3], []);
        }
        let legs = (0u8..3).fold(Legs::new(), |legs, label| {
            legs.leg(label, Stub::undirected(0))
        });
        (theory.build().unwrap(), Query::new(&legs).loops(0..=0))
    }

    fn node_pointer(graph: &HedgeGraph<u8, usize, Option<u8>>) -> usize {
        &graph[NodeIndex(0)] as *const usize as usize
    }

    #[test]
    fn complete_graph_cache_survives_repeated_views_selection_and_emission() {
        let (theory, query) = input(2);
        let native = Arc::new(Mutex::new(Vec::new()));
        let dressing = Arc::clone(&native);
        let selection = Arc::clone(&native);
        let query = query
            .prune_dressings(move |view| {
                if let Some(graph) = view.complete_graph() {
                    dressing.lock().unwrap().push((
                        Hook::Dressing,
                        graph[NodeIndex(0)],
                        node_pointer(graph),
                    ));
                }
                true
            })
            .select(move |view| {
                let graph = view.complete_graph().unwrap();
                selection.lock().unwrap().push((
                    Hook::Select,
                    graph[NodeIndex(0)],
                    node_pointer(graph),
                ));
                true
            })
            .deferred_hooks(HookMask {
                dressings: true,
                select: true,
                ..HookMask::default()
            });
        let plan = Plan::new(&theory, &query).unwrap();
        let run = Run::default();
        let mut driver = plan.batched(&run, Output::Diagrams).unwrap();
        let mut warm = BTreeMap::new();
        let mut emitted = Vec::new();
        loop {
            match driver.poll().unwrap() {
                Event::Batch(requests) => {
                    let mut replies = Vec::new();
                    for request in &requests {
                        let pointer = || match request.hook() {
                            Hook::Dressing => request
                                .dressing()
                                .unwrap()
                                .complete_graph()
                                .map(|g| g as *const _ as usize),
                            Hook::Select => request
                                .diagram()
                                .unwrap()
                                .complete_graph()
                                .map(|g| g as *const _ as usize),
                            Hook::Skeleton => unreachable!(),
                        };
                        if let Some(first) = pointer() {
                            assert_eq!(Some(first), pointer(), "fresh views share one graph");
                            let cache = request.snapshot.complete.as_ref().unwrap();
                            let graph = cache.get().unwrap();
                            let rule = graph[NodeIndex(0)];
                            let address = node_pointer(graph);
                            if request.hook() == Hook::Dressing {
                                assert!(warm
                                    .insert(rule, (Arc::downgrade(cache), address))
                                    .is_none());
                            } else {
                                let (previous, expected) = &warm[&rule];
                                assert!(std::ptr::eq(previous.as_ptr(), Arc::as_ptr(cache)));
                                assert_eq!(address, *expected);
                            }
                        } else {
                            assert!(
                                request.snapshot.complete.is_none(),
                                "partial hooks stay cold"
                            );
                        }
                        replies.push(request.reply(true));
                    }
                    drop(requests);
                    driver.resolve(&replies).unwrap();
                }
                Event::Diagram(diagram) => {
                    let graph = diagram
                        .graph
                        .get()
                        .expect("warm graph transferred before any accessor");
                    let rule = graph[NodeIndex(0)];
                    assert_eq!(
                        node_pointer(graph),
                        warm[&rule].1,
                        "move preserves graph storage"
                    );
                    assert!(
                        warm[&rule].0.upgrade().is_none(),
                        "unretained request cache was moved"
                    );
                    emitted.push(rule);
                }
                Event::Done => break,
                Event::Yield | Event::Progress(_) => {}
                Event::Skeleton(_) => unreachable!(),
            }
        }
        assert_eq!(emitted, vec![0, 1]);
        let native = native.lock().unwrap();
        assert_eq!(native.len(), 4);
        for &(_, rule, pointer) in native.iter() {
            assert_eq!(
                pointer, warm[&rule].1,
                "native and deferred hooks share the cache"
            );
        }
    }

    #[test]
    fn rejected_assignments_drop_their_graph_and_do_not_warm_the_next_rule() {
        let (theory, query) = input(3);
        let query = query.deferred_hooks(HookMask {
            dressings: true,
            select: true,
            ..HookMask::default()
        });
        let plan = Plan::new(&theory, &query).unwrap();
        let run = Run::default();
        let mut driver = plan.batched(&run, Output::Diagrams).unwrap();
        let mut seen = Vec::new();
        let mut emitted = Vec::new();
        loop {
            match driver.poll().unwrap() {
                Event::Batch(requests) => {
                    let mut replies = Vec::new();
                    for request in &requests {
                        let view = request.snapshot.view();
                        let keep = if view.is_complete() {
                            let cache = request.snapshot.complete.as_ref().unwrap();
                            if request.hook() == Hook::Dressing {
                                assert!(cache.get().is_none(), "new rule starts cold");
                            }
                            let rule =
                                request.snapshot.decoration.as_ref().unwrap().rules[0] as usize;
                            let graph = view.complete_graph().unwrap();
                            assert_eq!(graph[NodeIndex(0)], rule, "no stale rule payload");
                            seen.push((request.hook(), rule));
                            match request.hook() {
                                Hook::Dressing => rule != 0,
                                Hook::Select => rule != 1,
                                Hook::Skeleton => unreachable!(),
                            }
                        } else {
                            true
                        };
                        replies.push(request.reply(keep));
                    }
                    drop(requests);
                    driver.resolve(&replies).unwrap();
                    for slot in &driver.slots {
                        if let Job::Cursor(cursor) = &slot.job {
                            if cursor.phase == Phase::Decorate {
                                assert!(
                                    cursor.complete.is_none(),
                                    "rejection releases the current cache"
                                );
                            }
                        }
                    }
                }
                Event::Diagram(diagram) => emitted.push(diagram.rules()[0]),
                Event::Done => break,
                Event::Yield | Event::Progress(_) => {}
                Event::Skeleton(_) => unreachable!(),
            }
        }
        assert_eq!(
            seen,
            vec![
                (Hook::Dressing, 0),
                (Hook::Dressing, 1),
                (Hook::Select, 1),
                (Hook::Dressing, 2),
                (Hook::Select, 2)
            ]
        );
        assert_eq!(emitted, vec![2]);

        let (theory, query) = input(3);
        let plan = Plan::new(&theory, &query).unwrap();
        let mut count = plan.batched(&run, Output::Count).unwrap();
        loop {
            let event = count.poll().unwrap();
            for slot in &count.slots {
                if let Job::Cursor(cursor) = &slot.job {
                    assert!(
                        cursor.complete.is_none(),
                        "count without hooks allocates no graph cache"
                    );
                }
            }
            if matches!(event, Event::Done) {
                break;
            }
            assert!(matches!(event, Event::Yield | Event::Progress(_)));
        }
        assert_eq!(count.counts().unwrap().diagrams, 3u8.into());
    }

    #[test]
    fn ranked_requests_share_automata_and_tables_and_keep_original_tickets() {
        let theory = Theory::builder()
            .vertex(
                [
                    Stub::undirected('a'),
                    Stub::undirected('a'),
                    Stub::undirected('z'),
                ],
                [],
            )
            .build()
            .unwrap();
        let legs = Legs::new()
            .leg('k', Stub::undirected('a'))
            .leg('b', Stub::undirected('a'))
            .leg('k', Stub::undirected('z'));
        let query = Query::new(&legs).loops(0..=0).deferred_hooks(HookMask {
            skeletons: true,
            dressings: true,
            select: true,
        });
        let plan = Plan::new(&theory, &query).unwrap();
        let run = Run::default();
        let mut driver = plan.batched(&run, Output::Diagrams).unwrap();
        let mut ranked_owner: Option<Arc<Compiled<u32, u32>>> = None;
        let mut requests_seen = 0;
        let mut diagrams = 0;
        loop {
            match driver.poll().unwrap() {
                Event::Batch(requests) => {
                    let mut replies = Vec::new();
                    for request in requests {
                        assert_eq!(request.fields(), &['a', 'z']);
                        assert_eq!(request.labels(), &['b', 'k']);
                        let old = Arc::clone(&request.snapshot.compiled);
                        let ticket = request.ticket();
                        let hook = request.hook();
                        let expected_legs = {
                            let view = request.snapshot.view();
                            (0..view.n_vertices())
                                .flat_map(|v| {
                                    view.legs(v).map(move |(&l, s)| {
                                        (v, u32::from(l == 'k'), u32::from(*s.field == 'z'), s.flow)
                                    })
                                })
                                .collect::<Vec<_>>()
                        };
                        // A ranked graph cannot be reused as an original-key output graph.
                        if request.snapshot.complete.is_some() {
                            request.snapshot.view().complete_graph().unwrap();
                        }
                        let ranked = request.into_ranked();
                        assert_eq!(ranked.ticket(), ticket);
                        assert_eq!(ranked.hook(), hook);
                        assert_eq!(ranked.fields(), &[0, 1]);
                        assert_eq!(ranked.labels(), &[0, 1]);
                        let compiled = &ranked.snapshot.compiled;
                        assert!(Arc::ptr_eq(&old.model, &compiled.model));
                        assert!(Arc::ptr_eq(
                            &old.tables.signatures,
                            &compiled.tables.signatures
                        ));
                        assert!(Arc::ptr_eq(&old.tables.weights, &compiled.tables.weights));
                        if let Some(previous) = &ranked_owner {
                            assert!(Arc::ptr_eq(previous, compiled));
                            assert!(Arc::ptr_eq(&previous.tables, &compiled.tables));
                        } else {
                            ranked_owner = Some(Arc::clone(compiled));
                        }
                        let view = ranked.snapshot.view();
                        let actual_legs = (0..view.n_vertices())
                            .flat_map(|v| view.legs(v).map(move |(&l, s)| (v, l, *s.field, s.flow)))
                            .collect::<Vec<_>>();
                        assert_eq!(actual_legs, expected_legs);
                        if let Some(cache) = &ranked.snapshot.complete {
                            assert!(
                                cache.get().is_none(),
                                "key translation starts its own cache"
                            );
                            let first = view.complete_graph().unwrap() as *const _;
                            let second =
                                ranked.snapshot.view().complete_graph().unwrap() as *const _;
                            assert_eq!(first, second, "ranked repeated reads stay cached");
                        }
                        replies.push(ranked.reply(true));
                        requests_seen += 1;
                    }
                    driver.resolve(&replies).unwrap();
                }
                Event::Diagram(_) => diagrams += 1,
                Event::Done => break,
                Event::Yield | Event::Progress(_) => {}
                Event::Skeleton(_) => unreachable!(),
            }
        }
        assert!(requests_seen >= 3);
        assert_eq!(diagrams, 1);
    }
}
