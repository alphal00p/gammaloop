//! Compact outputs of the common skeleton/decoration traversal.

use std::{
    collections::BTreeMap,
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc, OnceLock,
    },
};

use num_bigint::{BigInt, BigUint};
use num_rational::BigRational;
use num_traits::One;

use super::super::canon::{factorial, scaled};

use super::{
    decorate::{Decoration, Decorator},
    topology::{Search, Topology},
    view::View,
    CanonicalForm, Compiled, Diagram, GenerationError, HedgeGraph, OutputTables, Plan, Progress,
    Run, RunContext, SkeletonView, Traversal,
};

/// Exact totals for one loop order. The inverse automorphism sum has no fermion signs or
/// process-specific orbit multiplicities.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct LoopCounts {
    pub skeletons: BigUint,
    pub diagrams: BigUint,
    pub sum_inv_aut: BigRational,
}

/// Exact totals, including whether execution exhausted the requested search.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub skeletons: BigUint,
    pub diagrams: BigUint,
    pub sum_inv_aut: BigRational,
    pub per_loop: BTreeMap<usize, LoopCounts>,
    pub complete: bool,
}

impl Counts {
    pub(super) fn diagram(
        &mut self,
        loops: usize,
        automorphisms: Option<u128>,
    ) -> Result<(), GenerationError> {
        let aut = automorphisms
            .filter(|&n| n != 0)
            .ok_or(GenerationError::AutomorphismOverflow)?;
        let weight = BigRational::new(BigInt::one(), BigInt::from(aut));
        let row = self.per_loop.entry(loops).or_default();
        row.diagrams += 1u8;
        row.sum_inv_aut += &weight;
        self.diagrams += 1u8;
        self.sum_inv_aut += weight;
        Ok(())
    }
}

/// A topology snapshot has no theory automata, callbacks, running search or decorator buffers.
/// Only the enumerating handle borrows a plan; diagrams keep this lightweight shared data.
pub(super) struct SkeletonData<F, L> {
    pub(super) topology: Topology,
    pub(super) tables: Arc<OutputTables<F, L>>,
    pub(super) loops: usize,
    pub(super) generators: Box<[u32]>,
    pub(super) node_group_order: Option<u128>,
    automorphisms: Option<u128>,
    graph: OnceLock<HedgeGraph<Option<F>, Option<usize>, Option<L>>>,
    canonical_form: OnceLock<CanonicalForm<Option<usize>, Option<F>, Option<L>>>,
}

impl<F, L> SkeletonData<F, L> {
    pub(super) fn new(compiled: &Compiled<F, L>, search: &Search, loops: usize) -> Self {
        let topology = search.topology();
        // The search group permutes vertices. Its kernel fixes vertices and permutes identical
        // legs/parallel edges, and may swap the two ends of each untyped self-loop.
        let mut aut = search.node_group_order();
        for u in 0..topology.n {
            let vertex = &topology.vertices[u];
            for class in 0..compiled.model.legs.n_classes() {
                aut = scaled(
                    aut,
                    factorial(compiled.model.legs.count(vertex.labelled, class) as usize),
                );
            }
            let loops = vertex.loops as usize;
            aut = scaled(scaled(aut, factorial(loops)), Some(1u128 << loops));
            for v in u + 1..topology.n {
                aut = scaled(aut, factorial(topology.multiplicity(u, v) as usize));
            }
        }
        Self {
            topology: topology.snapshot(),
            tables: Arc::clone(&compiled.tables),
            loops,
            generators: search.generators().into(),
            node_group_order: search.node_group_order(),
            automorphisms: aut,
            graph: OnceLock::new(),
            canonical_form: OnceLock::new(),
        }
    }

    pub(super) fn view(&self) -> SkeletonView<'_, F, L> {
        SkeletonView(View::stored(self))
    }
}

impl<F: Clone + Ord, L: Clone + Ord> SkeletonData<F, L> {
    fn graph(&self) -> &HedgeGraph<Option<F>, Option<usize>, Option<L>> {
        self.graph.get_or_init(|| self.view().graph())
    }

    /// Canonical incidence, not the traversal ordinal, supplies stable skeleton identity.
    pub(super) fn id(&self) -> &CanonicalForm<Option<usize>, Option<F>, Option<L>> {
        self.canonical_form.get_or_init(|| {
            self.graph()
                .canonize_by(Clone::clone, Clone::clone, Clone::clone)
                .canonical_form()
                .clone()
        })
    }
}

/// One accepted complete skeleton and the compiled request used to dress it.
/// The handle is immutable as its originating iterator advances.
pub struct Skeleton<'a, F, L> {
    snapshot: SkeletonSnapshot<F, L>,
    plan: &'a Plan<F, L>,
    run: &'a Run<'a>,
}

/// An owned skeleton without a live execution or plan. Retain its originating plan to start
/// another dressing traversal with [`Plan::dressings`]. Cloning shares the compact snapshot.
pub struct SkeletonSnapshot<F, L> {
    pub(super) data: Arc<SkeletonData<F, L>>,
}

impl<F, L> SkeletonSnapshot<F, L> {
    /// Identity of this shared allocation, valid only while an owning snapshot is retained.
    /// This is an execution-local storage key, never a canonical identifier or persisted ID.
    #[doc(hidden)]
    pub fn storage_id(&self) -> usize {
        Arc::as_ptr(&self.data) as usize
    }
}

impl<F, L> Clone for SkeletonSnapshot<F, L> {
    fn clone(&self) -> Self {
        Self {
            data: Arc::clone(&self.data),
        }
    }
}

impl<F, L> std::ops::Deref for Skeleton<'_, F, L> {
    type Target = SkeletonSnapshot<F, L>;
    fn deref(&self) -> &Self::Target {
        &self.snapshot
    }
}

impl<F: Clone + Ord, L: Clone + Ord> SkeletonSnapshot<F, L> {
    pub fn id(&self) -> &CanonicalForm<Option<usize>, Option<F>, Option<L>> {
        self.data.id()
    }
    pub fn graph(&self) -> &HedgeGraph<Option<F>, Option<usize>, Option<L>> {
        self.data.graph()
    }
    pub fn automorphisms(&self) -> Option<u128> {
        self.data.automorphisms
    }
    pub fn view(&self) -> SkeletonView<'_, F, L> {
        self.data.view()
    }
}

impl<F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> Skeleton<'_, F, L> {
    pub fn snapshot(&self) -> SkeletonSnapshot<F, L> {
        self.snapshot.clone()
    }

    /// Start a fresh dressing traversal of this snapshot; no topology search is repeated.
    pub fn diagrams(&self) -> Dressings<'_, F, L> {
        Dressings::new(self.plan, Arc::clone(&self.snapshot.data), self.run)
    }
}

impl<'a, F, L> Dressings<'a, F, L> {
    fn new(plan: &'a Plan<F, L>, data: Arc<SkeletonData<F, L>>, run: &'a Run<'a>) -> Self {
        let mut decorator = Decorator::default();
        decorator.start(
            &plan.compiled.model,
            &data.topology,
            &data.generators,
            data.node_group_order,
        );
        let loops = data.loops;
        Self {
            data,
            plan,
            run,
            decorator,
            started: false,
            stop: AtomicBool::new(false),
            progress: Progress {
                loops,
                ..Progress::default()
            },
        }
    }
}

/// Accepted diagrams from one complete skeleton, sharing final acceptance with `Plan::diagrams`.
pub struct Dressings<'a, F, L> {
    data: Arc<SkeletonData<F, L>>,
    plan: &'a Plan<F, L>,
    run: &'a Run<'a>,
    decorator: Decorator,
    started: bool,
    stop: AtomicBool,
    progress: Progress,
}

impl<F, L> Dressings<'_, F, L> {
    pub fn aborted(&self) -> bool {
        self.stop.load(Ordering::Relaxed)
    }
}

impl<F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> Iterator for Dressings<'_, F, L> {
    type Item = Diagram<F, L>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.aborted() {
            return None;
        }
        if self.run.abort_check.as_ref().is_some_and(|check| check()) {
            self.stop.store(true, Ordering::Relaxed);
            return None;
        }
        if !self.started {
            self.started = true;
            self.progress.topologies = 1u8.into();
            super::report(self.run.progress.as_deref(), &self.progress, &self.stop);
            if self.aborted() {
                return None;
            }
        }
        let current = self
            .plan
            .orders
            .iter()
            .find(|(loops, _)| *loops == self.data.loops)
            .cloned();
        let cx = RunContext::new(self.plan, self.run, &current, &self.stop);
        loop {
            if cx.stopped() {
                return None;
            }
            let Some(d) = self.decorator.next(
                &self.plan.compiled.model,
                &self.data.topology,
                cx.prune(),
                || cx.stopped(),
            ) else {
                cx.stopped();
                return None;
            };
            if let Some(view) = cx.accept(d) {
                let diagram =
                    self.plan
                        .compiled
                        .diagram(d, self.data.loops, view, Arc::clone(&self.data));
                self.progress.diagrams += 1u8;
                super::report(self.run.progress.as_deref(), &self.progress, &self.stop);
                return Some(diagram);
            }
        }
    }
}

pub struct Skeletons<'a, F: Send + Sync, L: Send + Sync> {
    pub(super) traversal: Traversal<'a, F, L, SkeletonOutput>,
}

impl<F: Send + Sync, L: Send + Sync> Skeletons<'_, F, L> {
    pub fn aborted(&self) -> bool {
        self.traversal.stop.load(Ordering::Relaxed)
    }
}

impl<'a, F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> Iterator
    for Skeletons<'a, F, L>
{
    type Item = Skeleton<'a, F, L>;
    fn next(&mut self) -> Option<Self::Item> {
        self.traversal.next().map(|data| Skeleton {
            snapshot: SkeletonSnapshot { data },
            plan: self.traversal.plan,
            run: self.traversal.controls,
        })
    }
}

// Private output policies let the same search and acceptance boundary count without freezing a
// Diagram, or yield a skeleton before starting its decorator. They do not alter search decisions.
pub(super) trait Emit<F, L>: Send {
    type Item: Send;
    const DECORATE: bool;
    const SNAPSHOT: bool;
    fn skeleton(data: Arc<SkeletonData<F, L>>) -> Option<Self::Item>;
    fn diagram(
        compiled: &Compiled<F, L>,
        d: Decoration<'_>,
        loops: usize,
        view: View<'_, F, L>,
        skeleton: Option<&Arc<SkeletonData<F, L>>>,
    ) -> Self::Item;
}

pub(super) struct DiagramOutput;
pub(super) struct SkeletonOutput;
pub(super) struct CountOutput;

impl<F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> Emit<F, L> for DiagramOutput {
    type Item = Diagram<F, L>;
    const DECORATE: bool = true;
    const SNAPSHOT: bool = true;
    fn skeleton(_: Arc<SkeletonData<F, L>>) -> Option<Self::Item> {
        None
    }
    fn diagram(
        compiled: &Compiled<F, L>,
        d: Decoration<'_>,
        loops: usize,
        view: View<'_, F, L>,
        skeleton: Option<&Arc<SkeletonData<F, L>>>,
    ) -> Self::Item {
        compiled.diagram(
            d,
            loops,
            view,
            Arc::clone(skeleton.expect("diagram mode snapshots its skeleton")),
        )
    }
}

impl<F: Send + Sync, L: Send + Sync> Emit<F, L> for SkeletonOutput {
    type Item = Arc<SkeletonData<F, L>>;
    const DECORATE: bool = false;
    const SNAPSHOT: bool = true;
    fn skeleton(data: Arc<SkeletonData<F, L>>) -> Option<Self::Item> {
        Some(data)
    }
    fn diagram(
        _: &Compiled<F, L>,
        _: Decoration<'_>,
        _: usize,
        _: View<'_, F, L>,
        _: Option<&Arc<SkeletonData<F, L>>>,
    ) -> Self::Item {
        unreachable!("skeleton mode never starts a decorator")
    }
}

impl<F, L> Emit<F, L> for CountOutput {
    type Item = (usize, Option<u128>);
    const DECORATE: bool = true;
    const SNAPSHOT: bool = false;
    fn skeleton(_: Arc<SkeletonData<F, L>>) -> Option<Self::Item> {
        None
    }
    fn diagram(
        _: &Compiled<F, L>,
        d: Decoration<'_>,
        loops: usize,
        _: View<'_, F, L>,
        _: Option<&Arc<SkeletonData<F, L>>>,
    ) -> Self::Item {
        (loops, d.automorphism_count)
    }
}

impl<F: Clone + Ord + Send + Sync, L: Clone + Ord + Send + Sync> Plan<F, L> {
    /// Dress a retained snapshot from this plan with fresh execution controls. The snapshot
    /// does not keep a theory automaton alive and no topology search is repeated.
    pub fn dressings<'a>(
        &'a self,
        skeleton: &SkeletonSnapshot<F, L>,
        run: &'a Run<'a>,
    ) -> Result<Dressings<'a, F, L>, GenerationError> {
        self.direct_hooks()?;
        if !Arc::ptr_eq(&self.compiled.tables, &skeleton.data.tables) {
            return Err(GenerationError::ForeignSkeleton);
        }
        Ok(Dressings::new(self, Arc::clone(&skeleton.data), run))
    }

    pub fn skeletons<'a>(
        &'a self,
        run: &'a Run<'a>,
    ) -> Result<Skeletons<'a, F, L>, GenerationError> {
        self.direct_hooks()?;
        Ok(Skeletons {
            traversal: Traversal::new(self, run),
        })
    }

    /// Count accepted assignments without constructing diagrams, graphs or canonical forms.
    /// A hook may explicitly request a graph; that demand retains its usual cost.
    pub fn count(&self, run: &Run<'_>) -> Result<Counts, GenerationError> {
        self.direct_hooks()?;
        let mut traversal = Traversal::<_, _, CountOutput>::new(self, run);
        let mut counts = Counts::default();
        for (loops, aut) in traversal.by_ref() {
            counts.diagram(loops, aut)?;
        }
        for (loops, skeletons) in &traversal.topologies {
            counts.per_loop.entry(*loops).or_default().skeletons = skeletons.clone();
            counts.skeletons += skeletons;
        }
        counts.complete = !traversal.aborted();
        Ok(counts)
    }
}
