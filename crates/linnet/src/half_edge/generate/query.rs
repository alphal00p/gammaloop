//! Reusable theory input and the structural part of a generation request.

use std::{
    ops::{Bound, ControlFlow, RangeBounds},
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc,
    },
};

use super::{
    budget::inclusive, model, Countable, DiagramView, DressingView, GenerationError, Progress,
    SkeletonView, Stub,
};

type SkeletonPredicate<F, L> = dyn Fn(&SkeletonView<'_, F, L>) -> bool + Send + Sync;
type DressingPredicate<F, L> = dyn Fn(&DressingView<'_, F, L>) -> bool + Send + Sync;
type SelectionPredicate<F, L> = dyn Fn(&DiagramView<'_, F, L>) -> bool + Send + Sync;

#[derive(Clone)]
pub(super) struct Hooks<F, L> {
    pub(super) skeletons: Option<Arc<SkeletonPredicate<F, L>>>,
    pub(super) dressings: Option<Arc<DressingPredicate<F, L>>>,
    pub(super) select: Option<Arc<SelectionPredicate<F, L>>>,
    pub(super) deferred: super::batch::HookMask,
}

impl<F, L> Default for Hooks<F, L> {
    fn default() -> Self {
        Self {
            skeletons: None,
            dressings: None,
            select: None,
            deferred: super::batch::HookMask::default(),
        }
    }
}

/// Ordered vertex rules and their reusable structural automata.
///
/// Rules with equal stubs remain distinct: their insertion indices label generated vertices.
/// A theory has no external legs, selected order bounds, execution callbacks or search scratch.
#[derive(Clone)]
pub struct Theory<F> {
    pub(super) fields: Arc<[F]>,
    pub(super) dense: Arc<model::Theory>,
    /// Named columns in declaration order, each with one value per rule.
    pub(super) weights: Arc<[(String, Vec<i32>)]>,
    pub(super) n_rules: usize,
}

impl<F> Theory<F> {
    pub fn builder() -> TheoryBuilder<F> {
        TheoryBuilder {
            names: Vec::new(),
            signatures: Vec::new(),
            weights: Vec::new(),
        }
    }

    pub fn n_rules(&self) -> usize {
        self.n_rules
    }

    pub fn fields(&self) -> &[F] {
        &self.fields
    }

    pub fn weight_names(&self) -> impl ExactSizeIterator<Item = &str> {
        self.weights.iter().map(|(name, _)| name.as_str())
    }
}

/// Collect rules before compiling one immutable [`Theory`].
pub struct TheoryBuilder<F> {
    names: Vec<String>,
    signatures: Vec<Vec<Stub<F>>>,
    weights: Vec<Vec<i32>>,
}

impl<F> TheoryBuilder<F> {
    pub fn weights(mut self, names: impl IntoIterator<Item = impl Into<String>>) -> Self {
        self.names = names.into_iter().map(Into::into).collect();
        self
    }

    pub fn vertex(
        mut self,
        stubs: impl IntoIterator<Item = Stub<F>>,
        weights: impl IntoIterator<Item = i32>,
    ) -> Self {
        self.signatures.push(stubs.into_iter().collect());
        self.weights.push(weights.into_iter().collect());
        self
    }
}

impl<F: Clone + Ord> TheoryBuilder<F> {
    pub fn build(self) -> Result<Theory<F>, GenerationError> {
        let mut sorted_names: Vec<_> = self.names.iter().collect();
        sorted_names.sort_unstable();
        if let Some(duplicate) = sorted_names.windows(2).find(|pair| pair[0] == pair[1]) {
            return Err(GenerationError::DuplicateWeight(duplicate[0].clone()));
        }
        for (rule, weights) in self.weights.iter().enumerate() {
            if weights.len() != self.names.len() {
                return Err(GenerationError::VertexWeights {
                    rule,
                    columns: self.names.len(),
                    weights: weights.len(),
                });
            }
        }
        let mut fields: Vec<_> = self
            .signatures
            .iter()
            .flatten()
            .map(|stub| stub.field.clone())
            .collect();
        fields.sort();
        fields.dedup();
        if fields.len() > model::MAX_FIELDS {
            return Err(GenerationError::TooManyFields);
        }
        let signatures: Vec<Vec<_>> = self
            .signatures
            .iter()
            .map(|signature| {
                signature
                    .iter()
                    .map(|stub| {
                        let field = fields.binary_search(&stub.field).expect("interned field");
                        3 * field as u8 + stub.direction()
                    })
                    .collect()
            })
            .collect();
        let dense = Arc::new(model::Theory::new(fields.len(), &signatures)?);
        let weights = self
            .names
            .into_iter()
            .enumerate()
            .map(|(column, name)| (name, self.weights.iter().map(|row| row[column]).collect()))
            .collect();
        Ok(Theory {
            fields: fields.into(),
            dense,
            weights,
            n_rules: signatures.len(),
        })
    }
}

/// External slots. Equal labels and stubs form one interchangeable class.
#[derive(Clone, Debug)]
pub struct Legs<F, L> {
    pub(super) slots: Vec<(L, Stub<F>)>,
}

impl<F, L> Default for Legs<F, L> {
    fn default() -> Self {
        Self { slots: Vec::new() }
    }
}

impl<F, L> Legs<F, L> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn leg(mut self, label: L, stub: Stub<F>) -> Self {
        self.slots.push((label, stub));
        self
    }

    /// Repeating the same slot uses the same normalization as repeated [`Self::leg`] calls.
    pub fn identical(mut self, label: L, stub: Stub<F>, multiplicity: usize) -> Self
    where
        F: Clone,
        L: Clone,
    {
        self.slots
            .extend(std::iter::repeat_n((label, stub), multiplicity));
        self
    }

    pub fn len(&self) -> usize {
        self.slots.len()
    }

    pub fn is_empty(&self) -> bool {
        self.slots.is_empty()
    }
}

/// One home for bridge restrictions; dangling external legs are never bridges.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Bridges {
    #[default]
    Any,
    AtMost(usize),
    OnePI,
    /// Keep trees unrestricted; diagrams of positive loop order must be one-particle irreducible.
    OnePIInLoops,
}

impl Bridges {
    pub(super) fn limit(self, loops: usize) -> Option<usize> {
        match self {
            Self::Any => None,
            Self::AtMost(n) => Some(n),
            Self::OnePI => Some(0),
            Self::OnePIInLoops => (loops > 0).then_some(0),
        }
    }
}

/// Topology restrictions, independent of fields assigned during dressing.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Structure {
    pub bridges: Bridges,
    pub self_loops: bool,
    /// Allow bridges with no external leg on one side, including vacuum bridges.
    pub zero_flow_edges: bool,
    /// Required when a rule has degree at most two; otherwise a finite loop bound also suffices.
    pub max_vertices: Option<usize>,
}

impl Structure {
    pub fn one_pi() -> Self {
        Self {
            bridges: Bridges::OnePI,
            ..Self::default()
        }
    }
}

/// What to enumerate. Compiling a [`super::Plan`] validates this request once against a theory.
///
/// Pruning hooks must reject every refinement of a rejected state. Final selection may use any
/// predicate and does not affect which complete skeletons are exposed by skeleton streaming.
#[derive(Clone)]
pub struct Query<F, L> {
    pub(super) legs: Legs<F, L>,
    pub(super) loops: (usize, Option<usize>),
    pub(super) structure: Structure,
    pub(super) weights: Vec<(String, [i64; 2])>,
    pub(super) counts: Vec<(Countable<F>, [i64; 2])>,
    pub(super) hooks: Hooks<F, L>,
}

impl<F: Clone, L: Clone> Query<F, L> {
    /// Binding installation boundary; public predicates still receive one typed view.
    #[doc(hidden)]
    pub fn deferred_hooks(mut self, hooks: super::batch::HookMask) -> Self {
        self.hooks.deferred = hooks;
        self
    }

    pub fn new(legs: &Legs<F, L>) -> Self {
        Self {
            legs: legs.clone(),
            loops: (0, None),
            structure: Structure::default(),
            weights: Vec::new(),
            counts: Vec::new(),
            hooks: Hooks::default(),
        }
    }

    pub fn loops(mut self, range: impl RangeBounds<usize>) -> Self {
        let low = match range.start_bound() {
            Bound::Included(&n) => n,
            Bound::Excluded(&n) => match n.checked_add(1) {
                Some(next) => next,
                None => {
                    self.loops = (1, Some(0));
                    return self;
                }
            },
            Bound::Unbounded => 0,
        };
        let high = match range.end_bound() {
            Bound::Included(&n) => Some(n),
            Bound::Excluded(&n) => {
                // `..0` is empty, including when the lower bound is unbounded.
                if n == 0 {
                    self.loops = (1, Some(0));
                    return self;
                }
                Some(n - 1)
            }
            Bound::Unbounded => None,
        };
        self.loops = (low, high);
        self
    }

    pub fn structure(mut self, structure: Structure) -> Self {
        self.structure = structure;
        self
    }

    /// Select a declared weight column. Repeated constraints on a column are all enforced.
    pub fn weight(mut self, name: impl Into<String>, range: impl RangeBounds<i32>) -> Self {
        self.weights
            .push((name.into(), inclusive(&range, i64::from)));
        self
    }

    /// Counts share the same budget propagation as weight constraints.
    pub fn count(mut self, count: Countable<F>, range: impl RangeBounds<usize>) -> Self {
        self.counts.push((
            count,
            inclusive(&range, |n| i64::try_from(n).unwrap_or(i64::MAX)),
        ));
        self
    }

    pub fn prune_skeletons(
        mut self,
        predicate: impl Fn(&SkeletonView<'_, F, L>) -> bool + Send + Sync + 'static,
    ) -> Self {
        self.hooks.skeletons = Some(Arc::new(predicate));
        self
    }

    pub fn prune_dressings(
        mut self,
        predicate: impl Fn(&DressingView<'_, F, L>) -> bool + Send + Sync + 'static,
    ) -> Self {
        self.hooks.dressings = Some(Arc::new(predicate));
        self
    }

    pub fn select(
        mut self,
        predicate: impl Fn(&DiagramView<'_, F, L>) -> bool + Send + Sync + 'static,
    ) -> Self {
        self.hooks.select = Some(Arc::new(predicate));
        self
    }
}

/// Transient execution controls. Plans and owned diagrams do not retain a `Run`.
/// Live iterators and enumerating skeleton handles borrow it for their executions.
pub struct Run<'a> {
    pub(super) parallel: bool,
    pub(super) deferred_progress: bool,
    pub(super) check_heredity: bool,
    pub(super) abort_check: Option<Box<dyn Fn() -> bool + Send + Sync + 'a>>,
    pub(super) progress: Option<Box<dyn Fn(&Progress) -> ControlFlow<()> + Send + Sync + 'a>>,
}

impl Default for Run<'_> {
    fn default() -> Self {
        Self {
            parallel: false,
            deferred_progress: false,
            check_heredity: false,
            abort_check: None,
            progress: None,
        }
    }
}

impl<'a> Run<'a> {
    /// Ask the binding driver to yield progress notifications to the calling thread.
    #[doc(hidden)]
    pub fn deferred_progress(mut self, enabled: bool) -> Self {
        self.deferred_progress = enabled;
        self
    }

    /// On wasm and builds without Rayon, execution remains sequential with the same callbacks.
    pub fn parallel(mut self, parallel: bool) -> Self {
        self.parallel = parallel;
        self
    }

    pub fn abort(self, flag: &'a AtomicBool) -> Self {
        self.abort_check(move || flag.load(Ordering::Relaxed))
    }

    pub fn abort_check(mut self, check: impl Fn() -> bool + Send + Sync + 'a) -> Self {
        self.abort_check = Some(Box::new(check));
        self
    }

    /// Search past user vetoes and diagnose a later acceptance by the same pruning hook.
    /// Built-in pruning stays active; arbitrary final selection is checked separately.
    pub fn check_heredity(mut self, check: bool) -> Self {
        self.check_heredity = check;
        self
    }

    pub fn progress(
        mut self,
        callback: impl Fn(&Progress) -> ControlFlow<()> + Send + Sync + 'a,
    ) -> Self {
        self.progress = Some(Box::new(callback));
        self
    }
}

/// Where an effective constraint executes in a compiled plan.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlanStage {
    Planning,
    Skeletons,
    Dressing,
    Rules,
    Selection,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlanRole {
    Restriction,
    Propagator,
    FinalCheck,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PlanCriterion {
    pub criterion: String,
    pub effective: String,
    pub stage: PlanStage,
    pub hereditary: bool,
    pub role: PlanRole,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PlanReport {
    pub criteria: Vec<PlanCriterion>,
    pub empty_partitions: Vec<(usize, String)>,
    pub infeasible: Option<String>,
}

impl PlanReport {
    pub(super) fn criterion(
        &mut self,
        criterion: impl Into<String>,
        effective: impl Into<String>,
        stage: PlanStage,
        hereditary: bool,
        role: PlanRole,
    ) {
        self.criteria.push(PlanCriterion {
            criterion: criterion.into(),
            effective: effective.into(),
            stage,
            hereditary,
            role,
        });
    }

    pub(super) fn infeasible(&mut self, reason: impl Into<String>) -> GenerationError {
        self.infeasible = Some(reason.into());
        GenerationError::Infeasible(self.clone())
    }
}

impl std::fmt::Display for PlanReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for entry in &self.criteria {
            writeln!(
                f,
                "{} = {} -> {:?}; hereditary={}; {:?}",
                entry.criterion, entry.effective, entry.stage, entry.hereditary, entry.role
            )?;
        }
        for (loops, reason) in &self.empty_partitions {
            writeln!(f, "loops {loops}: empty ({reason})")?;
        }
        if let Some(reason) = &self.infeasible {
            writeln!(f, "infeasible: {reason}")?;
        }
        Ok(())
    }
}
