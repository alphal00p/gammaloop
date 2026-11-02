//! Owned physical preparation, explicit numerator construction, selection and grouping.

use super::{
    Arc, BTreeMap, BTreeSet, BigInt, BigRational, BigUint, CompactDiagram, ControlFlow, DiagramCut,
    DiagramGroup, DiagramInput, DiagramPair, DiagramThresholdCandidate, EdgeId, EngineRun,
    ExternalLabel, ExternalNode, FermionSigns, FeynmanDiagram, GeneratedEdge, GenerationError,
    GenerationReport, GenerationType, Generator, Grouping, Hash, HashMap, HedgeGraph, Model, Mutex,
    Naming, NumeratorOptions, ParticleId, PreparedEdge, PreparedSearch, Request, Run, Selection,
    SkeletonSnapshot, SourceIdentity, Stage, ThreadPoolBuilder, VertexRuleId,
    unresolved_cut_content,
};
#[cfg(test)]
use crate::generation::signed_routing::SignedRoutedDiagram;
use crate::grouping;
use crate::grouping::DiagramLabels;
use crate::{CutAmplitude, Range, Restrict, Symmetrize};
use feynkit_amplitude::{Amplitude, AmplitudeOptions, SquaredPairing};
use linnet::half_edge::generate::PlanReport;
use linnet::half_edge::generate::batch::{Event, Output};
use linnet::half_edge::generate::{PlanCriterion, PlanRole, PlanStage};
#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use std::collections::HashSet;
use std::fmt::{Debug, Formatter, Result as FmtResult};
use std::slice::Iter;

#[cfg(test)]
use super::cp_merge_diagnostic;
#[cfg(test)]
use serde_json::json;

/// Exact physical totals at one loop order, before numerator construction.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct LoopCounts {
    pub skeletons: BigUint,
    pub dressings: BigUint,
    pub diagrams: BigUint,
    pub sum_inv_aut: BigRational,
    pub physical_weight: BigRational,
}

/// Exact totals for the observed search. An interrupted prefix is never a complete result.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub skeletons: BigUint,
    pub dressings: BigUint,
    pub diagrams: BigUint,
    pub sum_inv_aut: BigRational,
    pub physical_weight: BigRational,
    pub per_loop: BTreeMap<usize, LoopCounts>,
    pub complete: bool,
}

/// A physically prepared assignment, with cuts and signs but no diagram numerator or ID.
#[derive(Clone)]
pub struct PreparedDiagram {
    input: Arc<DiagramInput>,
    compact: Arc<CompactDiagram<ParticleId, ExternalLabel>>,
    skeleton: SkeletonSnapshot<ParticleId, ExternalLabel>,
    orders: BTreeMap<String, usize>,
}

impl Debug for PreparedDiagram {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.debug_struct("PreparedDiagram")
            .field("loops", &self.loops())
            .field("aut", &self.aut())
            .field("orders", &self.orders)
            .finish_non_exhaustive()
    }
}

impl PreparedDiagram {
    pub fn rules(&self) -> impl ExactSizeIterator<Item = VertexRuleId> + '_ {
        self.input.representative.graph.vertices()
    }
    pub fn edges(&self) -> impl ExactSizeIterator<Item = PreparedEdge> + '_ {
        self.input.representative.graph.edges()
    }
    pub fn legs(&self) -> impl Iterator<Item = (usize, ExternalNode, usize)> + '_ {
        self.input.representative.graph.externals()
    }
    pub fn loops(&self) -> usize {
        self.compact.loops
    }
    pub fn aut(&self) -> u128 {
        self.input.representative.symmetry
    }
    pub fn skeleton(&self) -> &SkeletonSnapshot<ParticleId, ExternalLabel> {
        &self.skeleton
    }
    /// An owned graph of the physical representative, built on demand: prepared diagrams
    /// share their skeleton's incidence and keep no graph of their own.
    pub fn to_graph(&self) -> HedgeGraph<ParticleId, VertexRuleId, Option<ExternalNode>> {
        self.input.representative.graph.to_hedge_graph()
    }
    pub fn cuts(&self) -> &[DiagramCut] {
        &self.input.representative.cut_partitions
    }
    pub fn orders(&self) -> &BTreeMap<String, usize> {
        &self.orders
    }
}

/// The one resolved relation survives selection, color filtering and grouping, including empty
/// partner classes. It never changes an absent partner into unrestricted interference.
#[derive(Clone, Debug)]
struct OrderClasses {
    orders: Vec<BTreeMap<String, usize>>,
    allowed: BTreeSet<(usize, usize)>,
}

impl OrderClasses {
    fn resolve(
        request: &Request,
        orders: impl IntoIterator<Item = BTreeMap<String, usize>>,
    ) -> Option<Self> {
        if request.kind != GenerationType::Amplitude || request.orders.squared.is_empty() {
            return None;
        }
        let orders: Vec<_> = orders
            .into_iter()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let mut allowed = BTreeSet::new();
        for (a, left) in orders.iter().enumerate() {
            for (b, right) in orders.iter().enumerate() {
                // A single partner must satisfy every order constraint simultaneously.
                if request.orders.squared.iter().all(|(name, range)| {
                    left.get(name)
                        .copied()
                        .unwrap_or(0)
                        .checked_add(right.get(name).copied().unwrap_or(0))
                        .is_some_and(|sum| range.contains(sum))
                }) {
                    allowed.insert((a, b));
                }
            }
        }
        Some(Self { orders, allowed })
    }
    fn class(&self, orders: &BTreeMap<String, usize>) -> usize {
        self.orders
            .binary_search(orders)
            .expect("observed order class")
    }
    fn retained(&self, class: usize) -> bool {
        self.allowed.iter().any(|&(a, _)| a == class)
    }
    fn pairing(
        &self,
        classes: impl IntoIterator<Item = usize>,
    ) -> Result<SquaredPairing, GenerationError> {
        SquaredPairing::new(self.orders.len(), classes, self.allowed.iter().copied())
            .map_err(|error| GenerationError::Amplitude(Box::new(error)))
    }
}

struct SetContext {
    generator: Generator,
    search: PreparedSearch,
}

/// Compact physical outputs. Cloning a set shares its immutable prepared assignments.
#[derive(Clone)]
pub struct DiagramSet {
    context: Arc<SetContext>,
    diagrams: Vec<PreparedDiagram>,
    skeletons: Vec<SkeletonSnapshot<ParticleId, ExternalLabel>>,
    report: GenerationReport,
    classes: Option<OrderClasses>,
}
impl Debug for DiagramSet {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.debug_struct("DiagramSet")
            .field("diagrams", &self.diagrams)
            .field("report", &self.report)
            .finish_non_exhaustive()
    }
}
impl DiagramSet {
    pub fn len(&self) -> usize {
        self.diagrams.len()
    }
    pub fn is_empty(&self) -> bool {
        self.diagrams.is_empty()
    }
    pub fn get(&self, index: usize) -> Option<&PreparedDiagram> {
        self.diagrams.get(index)
    }
    pub fn iter(&self) -> Iter<'_, PreparedDiagram> {
        self.diagrams.iter()
    }
    pub fn report(&self) -> &GenerationReport {
        &self.report
    }
    pub fn skeletons(&self) -> &[SkeletonSnapshot<ParticleId, ExternalLabel>] {
        &self.skeletons
    }

    pub fn build(
        &self,
        options: &NumeratorOptions,
        run: &Run,
    ) -> Result<BuiltDiagramSet, GenerationError> {
        if matches!(options.naming, Naming::None) && !options.named_loop_momentum_bases.is_empty() {
            return Err(GenerationError::NamesRequireNamedOutput);
        }
        run.report_progress(Stage::Numerators, 0usize, Some(self.len()));
        let built = Mutex::new(0usize);
        let build = |(index, input): (usize, &PreparedDiagram)| -> Result<_, GenerationError> {
            if run.cancellation_requested() {
                return Ok(None);
            }
            let name = options.naming.name(index, self.len());
            // Equal options make every pair share one finalized physical owner.
            let mut diagram = self
                .context
                .generator
                .diagram_pair(
                    name,
                    self.context.search.request.kind,
                    (*input.input).clone(),
                    (options, options),
                )?
                .into_representative()?;
            let basis = if options.loop_momentum_bases.is_empty() {
                None
            } else {
                options.loop_momentum_bases.get(&diagram.id()?)
            }
            .or_else(|| {
                diagram
                    .name()
                    .and_then(|name| options.named_loop_momentum_bases.get(name))
            });
            if let Some(basis) = basis {
                diagram = diagram.with_loop_momentum_edges(basis)?;
            }
            #[cfg(test)]
            cp_merge_diagnostic::record(|| {
                json!({
                    "stage": "built_source", "prepared_source": index, "diagram": cp_merge_diagnostic::operator(&diagram),
                })
            });
            let mut built = built.lock().unwrap();
            *built += 1;
            run.report_progress(Stage::Numerators, *built, Some(self.len()));
            Ok(Some((diagram, index)))
        };
        // A cancelled build keeps the diagrams it completed, in their prepared order.
        let (diagrams, sources): (Vec<_>, Vec<_>) = run
            .in_pool(|| run.map_ordered(self.diagrams.iter().enumerate().collect(), build))?
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
            .unzip();
        let mut report = self.report.clone();
        report.completed &= !run.cancellation_requested();
        let mut result = BuiltDiagramSet {
            diagrams,
            groups: Vec::new(),
            report,
            prepared: Some(Arc::new(self.clone())),
            relation: self.classes.clone(),
            options: options.clone(),
            sources,
            classes: Vec::new(),
        };
        result.classes = result
            .sources
            .iter()
            .map(|&source| {
                self.classes
                    .as_ref()
                    .map_or(0, |classes| classes.class(&self.diagrams[source].orders))
            })
            .collect();
        if options.filter_zero_color {
            let filtered = run.in_pool(|| {
                grouping::group_diagrams(
                    std::mem::take(&mut result.diagrams)
                        .into_iter()
                        .map(DiagramPair::Shared)
                        .collect(),
                    &self.context.generator.model,
                    &Grouping::None,
                    false,
                    run,
                    true,
                    DiagramLabels::default(),
                )
            })??;
            let indices = filtered
                .groups
                .iter()
                .map(|group| group.members[0].source_diagram)
                .collect::<Vec<_>>();
            result.sources = indices.iter().map(|&i| result.sources[i]).collect();
            result.classes = indices.iter().map(|&i| result.classes[i]).collect();
            result.diagrams = filtered.diagrams;
            result.report.zero_numerator_count = filtered.zero_numerator_count.into();
        }
        result.singleton_groups()?;
        result.report.retained_count = result.len().into();
        Ok(result.done(run))
    }
}

impl DiagramSet {
    /// The finalized signed geometries, named as [`Self::build`] names its diagrams, without
    /// local Feynman rules or projectors. The routed stage is crate-internal until it gets its
    /// public `route()` boundary.
    #[cfg(test)]
    pub(super) fn routed(
        &self,
        naming: &Naming,
    ) -> Result<Vec<SignedRoutedDiagram>, GenerationError> {
        self.diagrams
            .iter()
            .enumerate()
            .map(|(index, diagram)| {
                self.context.generator.routed_diagram(
                    naming.name(index, self.len()),
                    self.context.search.request.kind,
                    (*diagram.input).clone(),
                )
            })
            .collect()
    }
}

/// Built physical diagrams and their source provenance. Selection and grouping keep the
/// resolved squared-order relation; raw lists passed to `Amplitude::new` remain unrestricted.
#[derive(Clone)]
pub struct BuiltDiagramSet {
    diagrams: Vec<FeynmanDiagram>,
    groups: Vec<DiagramGroup>,
    report: GenerationReport,
    /// The physical preparations that grouping rebuilds its comparison frames from; released
    /// once the set is grouped.
    prepared: Option<Arc<DiagramSet>>,
    relation: Option<OrderClasses>,
    options: NumeratorOptions,
    sources: Vec<usize>,
    classes: Vec<usize>,
}
impl Debug for BuiltDiagramSet {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        f.debug_struct("BuiltDiagramSet")
            .field("diagrams", &self.diagrams)
            .field("groups", &self.groups)
            .field("report", &self.report)
            .finish_non_exhaustive()
    }
}
impl BuiltDiagramSet {
    pub fn len(&self) -> usize {
        self.diagrams.len()
    }
    pub fn is_empty(&self) -> bool {
        self.diagrams.is_empty()
    }
    pub fn get(&self, index: usize) -> Option<&FeynmanDiagram> {
        self.diagrams.get(index)
    }
    pub fn iter(&self) -> Iter<'_, FeynmanDiagram> {
        self.diagrams.iter()
    }
    pub fn diagrams(&self) -> &[FeynmanDiagram] {
        &self.diagrams
    }
    pub fn groups(&self) -> &[DiagramGroup] {
        &self.groups
    }
    pub fn report(&self) -> &GenerationReport {
        &self.report
    }
    pub fn into_amplitude(self, options: AmplitudeOptions) -> Result<Amplitude, GenerationError> {
        if !self.report.completed {
            return Err(GenerationError::IncompleteAmplitude);
        }
        let pairing = self
            .relation
            .as_ref()
            .map(|relation| relation.pairing(self.classes.iter().copied()))
            .transpose()?;
        let amplitude =
            Amplitude::new(self.diagrams.into_iter().map(Arc::new), options).map_err(Box::new)?;
        Ok(match pairing {
            Some(pairing) => amplitude.with_pairing(pairing).map_err(Box::new)?,
            None => amplitude,
        })
    }
    fn singleton_groups(&mut self) -> Result<(), GenerationError> {
        self.groups = self
            .diagrams
            .iter()
            .zip(&self.sources)
            .enumerate()
            .map(|(i, (diagram, &source))| {
                Ok(DiagramGroup::singleton(
                    i,
                    source,
                    source_identity(diagram)?,
                    diagram.overall_factor().clone(),
                ))
            })
            .collect::<Result<_, GenerationError>>()?;
        self.validate_groups()
    }
    /// Whether grouping collapsed this set, which then no longer retains its preparations.
    fn grouped(&self) -> bool {
        self.prepared.is_none()
    }
    /// Finish a stage: report it done and keep its interruption state.
    fn done(mut self, run: &Run) -> Self {
        run.report_progress_with_zero_count(
            Stage::Done,
            self.len(),
            Some(self.len()),
            self.report.zero_numerator_count.clone(),
        );
        self.report.completed &= !run.cancellation_requested();
        self
    }
    /// Select a subset while reporting scanned inputs and retaining interruption state.
    pub fn select(&self, selection: &Selection, run: &Run) -> Result<Self, GenerationError> {
        if matches!(self.options.naming, Naming::None)
            && (selection
                .selected_diagram_names
                .as_ref()
                .is_some_and(|names| !names.is_empty())
                || !selection.vetoed_diagram_names.is_empty())
        {
            return Err(GenerationError::NamesRequireNamedOutput);
        }
        if self.grouped() && !selection.forced_cuts.is_empty() {
            return Err(GenerationError::InvalidRequest(
                "forced cuts must be selected before grouping".into(),
            ));
        }
        let mut result = Self {
            diagrams: Vec::new(),
            groups: Vec::new(),
            report: self.report.clone(),
            prepared: self.prepared.clone(),
            relation: self.relation.clone(),
            options: self.options.clone(),
            sources: Vec::new(),
            classes: Vec::new(),
        };
        let mut matched = vec![false; selection.forced_cuts.len()];
        run.report_progress_with_zero_count(
            Stage::Selection,
            0usize,
            Some(self.len()),
            result.report.zero_numerator_count.clone(),
        );
        run.in_pool(|| -> Result<(), GenerationError> {
            for (index, diagram) in self.diagrams.iter().enumerate() {
                if run.cancellation_requested() {
                    break;
                }
                'candidate: {
                    let selected_id = selection
                        .selected_diagram_ids
                        .as_ref()
                        .filter(|ids| !ids.is_empty())
                        .map(|ids| diagram.id().map(|id| ids.contains(&id)))
                        .transpose()?
                        .unwrap_or(false);
                    let selected_name = diagram.name().is_some_and(|name| {
                        selection
                            .selected_diagram_names
                            .as_ref()
                            .is_some_and(|names| names.contains(name))
                    });
                    if (selection.selected_diagram_ids.is_some()
                        || selection.selected_diagram_names.is_some())
                        && !selected_id
                        && !selected_name
                    {
                        break 'candidate;
                    }
                    if !selection.vetoed_diagram_ids.is_empty()
                        && selection.vetoed_diagram_ids.contains(&diagram.id()?)
                    {
                        break 'candidate;
                    }
                    if diagram
                        .name()
                        .is_some_and(|name| selection.vetoed_diagram_names.contains(name))
                    {
                        break 'candidate;
                    }
                    let mut diagram = diagram.clone();
                    if !selection.forced_cuts.is_empty() {
                        let cuts = diagram
                            .cuts()
                            .iter()
                            .filter(|cut| {
                                let edges = cut.cut.iter().map(|h| h.edge).collect::<BTreeSet<_>>();
                                let mut keep = false;
                                for (i, requested) in selection.forced_cuts.iter().enumerate() {
                                    if *requested == edges {
                                        matched[i] = true;
                                        keep = true
                                    }
                                }
                                keep
                            })
                            .cloned()
                            .collect::<Vec<_>>();
                        if cuts.is_empty() {
                            break 'candidate;
                        }
                        diagram = diagram.with_cuts(cuts)?;
                    }
                    let output = result.len();
                    let mut group = self.groups[index].clone();
                    group.master = output;
                    for member in &mut group.members {
                        member.diagram = output
                    }
                    result.diagrams.push(diagram);
                    result.sources.push(self.sources[index]);
                    result.classes.push(self.classes[index]);
                    result.groups.push(group);
                }
                run.report_progress_with_zero_count(
                    Stage::Selection,
                    index + 1,
                    Some(self.len()),
                    result.report.zero_numerator_count.clone(),
                );
            }
            Ok(())
        })??;
        result.report.completed &= !run.cancellation_requested();
        if result.report.completed
            && let Some((_, cut)) = selection
                .forced_cuts
                .iter()
                .enumerate()
                .find(|(i, _)| !matched[*i])
        {
            return Err(GenerationError::UnknownForcedCut {
                edges: cut.iter().map(|e| e.0).collect(),
            });
        }
        result.report.retained_count = result.len().into();
        result.validate_groups()?;
        Ok(result.done(run))
    }
}

/// The identity a group member records for its named source diagram.
fn source_identity(diagram: &FeynmanDiagram) -> Result<Option<SourceIdentity>, GenerationError> {
    diagram
        .name()
        .map(|name| {
            Ok(SourceIdentity {
                id: diagram.id()?,
                name: name.to_owned(),
            })
        })
        .transpose()
}

/// Enumerated assignments awaiting physical preparation. A pool with several threads prepares
/// them in bounded batches; one thread prepares each assignment as it is found.
struct Batches<A> {
    pending: Vec<CompactDiagram<ParticleId, ExternalLabel>>,
    size: usize,
    accept: A,
}

impl<A> Batches<A>
where
    A: FnMut(Vec<CompactDiagram<ParticleId, ExternalLabel>>) -> Result<(), GenerationError>,
{
    /// Batches for the run's pool.
    fn new(run: &Run, accept: A) -> Self {
        let threads = if run.on_caller() || cfg!(target_arch = "wasm32") {
            1
        } else {
            rayon::current_num_threads()
        };
        let size = if threads == 1 { 1 } else { 8 * threads };
        Self {
            pending: Vec::with_capacity(size),
            size,
            accept,
        }
    }

    fn push(
        &mut self,
        diagram: CompactDiagram<ParticleId, ExternalLabel>,
    ) -> Result<(), GenerationError> {
        self.pending.push(diagram);
        if self.pending.len() < self.size {
            return Ok(());
        }
        self.flush()
    }

    fn flush(&mut self) -> Result<(), GenerationError> {
        if self.pending.is_empty() {
            return Ok(());
        }
        (self.accept)(std::mem::replace(
            &mut self.pending,
            Vec::with_capacity(self.size),
        ))
    }
}

impl Run {
    /// Whether an explicit single thread runs every stage serially on the caller, without
    /// creating a Rayon pool.
    fn on_caller(&self) -> bool {
        self.threads == Some(1)
    }

    /// Run `work` in a pool of the requested threads, or on the caller for a single thread.
    pub(crate) fn in_pool<T: Send>(
        &self,
        work: impl FnOnce() -> T + Send,
    ) -> Result<T, GenerationError> {
        if self.threads == Some(0) {
            return Err(GenerationError::InvalidThreadCount(0));
        }
        #[cfg(not(target_arch = "wasm32"))]
        if !self.on_caller() {
            let mut builder = ThreadPoolBuilder::new();
            if let Some(threads) = self.threads {
                builder = builder.num_threads(threads)
            }
            return Ok(builder.build()?.install(work));
        }
        Ok(work())
    }

    /// Map `items` in their order: in parallel on the current Rayon pool, serially on the caller
    /// for a single thread and on WebAssembly. The output order never depends on the thread
    /// count. Every item is mapped, so callers read failures in item order and report the first
    /// one, as a serial loop would.
    pub(crate) fn map_ordered<T: Send, U: Send>(
        &self,
        items: Vec<T>,
        map: impl Fn(T) -> U + Send + Sync,
    ) -> Vec<U> {
        #[cfg(not(target_arch = "wasm32"))]
        if !self.on_caller() {
            return items.into_par_iter().map(map).collect();
        }
        items.into_iter().map(map).collect()
    }
}

impl PreparedSearch {
    /// Prepare one enumerated assignment physically, or reject it.
    fn prepare(
        &self,
        generator: &Generator,
        compact: CompactDiagram<ParticleId, ExternalLabel>,
    ) -> Result<Option<PreparedDiagram>, GenerationError> {
        let Some(multiplicity) = self.externals.sewn_pair_orbit(&compact) else {
            return Ok(None);
        };
        let aut = compact.symmetry_factor()?;
        let graph = self.externals.graph(&compact, &self.signatures);
        let Some(input) = generator.diagram_input(
            &self.resolved,
            graph,
            aut,
            multiplicity,
            &self.request,
            &Grouping::None,
        )?
        else {
            return Ok(None);
        };
        let mut orders = BTreeMap::<String, usize>::new();
        for rule in input.representative.graph.vertices() {
            for (name, power) in generator
                .model
                .vertex_rule_by_id(rule)?
                .coupling_orders(&generator.model)
            {
                let entry = orders.entry(name).or_default();
                *entry = entry
                    .checked_add(power)
                    .ok_or(GenerationError::ArithmeticOverflow("coupling order"))?
            }
        }
        orders.retain(|_, power| *power != 0);
        let skeleton = compact.skeleton_snapshot();
        Ok(Some(PreparedDiagram {
            input: Arc::new(input),
            compact: Arc::new(compact),
            skeleton,
            orders,
        }))
    }

    /// Visit the prepared assignments in enumeration order; batches are prepared in parallel.
    /// Count uses this same physical preparation boundary and retains only exact orbit
    /// keys/class accumulators, never all prepared graphs or numerators.
    fn visit(
        &self,
        generator: &Generator,
        run: &Run,
        mut accept: impl FnMut(PreparedDiagram) -> Result<(), GenerationError> + Send,
    ) -> Result<(GenerationReport, BTreeMap<usize, LoopCounts>), GenerationError> {
        let mut prepared = BigUint::default();
        self.enumerate(run, |batch| {
            if run.cancellation_requested() {
                return Ok(());
            }
            // Accept in order, stopping where a serial run would: at the first failure, or
            // once cancellation follows a report.
            for diagram in run.map_ordered(batch, |compact| self.prepare(generator, compact)) {
                if let Some(diagram) = diagram? {
                    accept(diagram)?;
                }
                prepared += 1u8;
                run.report_progress(Stage::Prepare, prepared.clone(), None);
                if run.cancellation_requested() {
                    break;
                }
            }
            Ok(())
        })
    }

    /// Enumerate assignments in their deterministic order, in batches for the run's pool.
    pub(super) fn enumerate(
        &self,
        run: &Run,
        accept: impl FnMut(
            Vec<CompactDiagram<ParticleId, ExternalLabel>>,
        ) -> Result<(), GenerationError>
        + Send,
    ) -> Result<(GenerationReport, BTreeMap<usize, LoopCounts>), GenerationError> {
        let progress = Mutex::new((BigUint::default(), BTreeMap::<usize, LoopCounts>::new()));
        let engine = EngineRun::default()
            .parallel(run.parallel && !run.on_caller())
            .check_heredity(run.check_heredity)
            .progress(|update| {
                let mut progress = progress.lock().unwrap();
                let delta = &update.topologies - &progress.0;
                progress.0 = update.topologies.clone();
                progress.1.entry(update.loops).or_default().skeletons += delta;
                run.report_progress(Stage::Skeletons, update.topologies.clone(), None);
                if run.cancellation_requested() {
                    ControlFlow::Break(())
                } else {
                    ControlFlow::Continue(())
                }
            })
            .abort_check(|| run.cancellation_requested());
        let mut dressings = BigUint::default();
        let mut visit = |batches: &mut Batches<_>,
                         compact: CompactDiagram<ParticleId, ExternalLabel>|
         -> Result<(), GenerationError> {
            dressings += 1u8;
            progress
                .lock()
                .unwrap()
                .1
                .entry(compact.loops)
                .or_default()
                .dressings += 1u8;
            if run.cancellation_requested() {
                return Ok(());
            }
            run.report_progress(Stage::Dressing, dressings.clone(), None);
            if run.cancellation_requested() {
                return Ok(());
            }
            batches.push(compact)
        };
        let enumerate = || -> Result<(), GenerationError> {
            if self.externals.vetoed(&self.restrictions) {
                return Ok(());
            }
            let Some(plan) = &self.plan else {
                return Ok(());
            };
            let mut batches = Batches::new(run, accept);
            if let Some((_, dispatch)) = &self.request.deferred_hooks {
                let mut driver = plan.batched(&engine, Output::Diagrams)?;
                loop {
                    match driver.poll()? {
                        Event::Batch(requests) => match dispatch(requests) {
                            ControlFlow::Continue(replies) => driver.resolve(&replies)?,
                            ControlFlow::Break(()) => {
                                run.cancellation.cancel();
                                driver.abort();
                            }
                        },
                        Event::Diagram(diagram) => visit(&mut batches, *diagram)?,
                        Event::Yield | Event::Progress(_) => {
                            if run.cancellation_requested() {
                                driver.abort()
                            }
                        }
                        Event::Done => break,
                        Event::Skeleton(_) => unreachable!("diagram output requested"),
                    }
                }
            } else {
                for diagram in plan.diagrams(&engine)? {
                    visit(&mut batches, diagram)?;
                }
            }
            batches.flush()
        };
        run.in_pool(enumerate)??;
        drop(engine);
        let (topology_count, per_loop) = progress.into_inner().unwrap();
        Ok((
            GenerationReport {
                topology_count,
                interaction_assignment_count: dressings,
                completed: !run.cancellation_requested(),
                ..Default::default()
            },
            per_loop,
        ))
    }
}

/// A collision-free compact equality key used only when the requested mirror quotient requires
/// merging physical representatives. It owns no graph, rule table, numerator or canonical ID.
#[derive(PartialEq, Eq, Hash)]
struct OrbitKey {
    rules: Vec<VertexRuleId>,
    edges: Vec<GeneratedEdge>,
    cuts: Vec<DiagramCut>,
    thresholds: Vec<DiagramThresholdCandidate>,
    aut: u128,
    signs: FermionSigns,
}
impl OrbitKey {
    fn of(input: &DiagramInput) -> Self {
        let graph = &input.representative;
        Self {
            rules: graph.graph.vertices().collect(),
            edges: graph.graph.edges().collect(),
            cuts: graph.cut_partitions.clone(),
            thresholds: graph.topology_threshold_candidates.clone(),
            aut: graph.symmetry,
            signs: input.signs,
        }
    }
}

impl Request {
    pub fn explain(&self, model: &Model) -> Result<PlanReport, GenerationError> {
        let search = Generator::new(Arc::new(model.clone())).search(self, |_, _| {}, |_, _| {})?;
        let mut report = search.report;
        let mut criterion = |name: &str, effective: String, stage, role| {
            report.criteria.push(PlanCriterion {
                criterion: name.into(),
                effective,
                stage,
                // These entries describe completed-state checks or compiled restrictions,
                // never an additional hereditary veto on a partial native search state.
                hereditary: false,
                role,
            });
        };
        if self.restrictions != Restrict::default() {
            criterion(
                "physical restrictions",
                format!(
                    "resolved before theory compilation: {:?}; external legs checked before enumeration",
                    search.restrictions
                ),
                PlanStage::Planning,
                PlanRole::Restriction,
            );
            if search.restrictions.vertices.is_some() || !search.restrictions.vertex_veto.is_empty()
            {
                criterion(
                    "interaction restriction after CP repair",
                    format!(
                        "allow={:?}, exclude={:?}",
                        search.restrictions.vertices, search.restrictions.vertex_veto
                    ),
                    PlanStage::Preparation,
                    PlanRole::FinalCheck,
                );
            }
        }
        let topology = &self.topology;
        if topology.self_energies.enabled()
            || topology.tadpoles.enabled()
            || topology.snails.enabled()
            || topology.factorized_loops != Range::Unrestricted
        {
            for stage in [PlanStage::Skeletons, PlanStage::Dressing] {
                criterion(
                    "physical topology vetoes",
                    format!(
                        "complete states only; self_energies={:?}, tadpoles={:?}, snails={:?}, factorized_loops={:?}; unknown mass vetoed only when both masses are excluded",
                        topology.self_energies,
                        topology.tadpoles,
                        topology.snails,
                        topology.factorized_loops
                    ),
                    stage,
                    PlanRole::FinalCheck,
                );
            }
        }
        if self.symmetrization != Symmetrize::default() {
            criterion(
                "external symmetrization",
                format!(
                    "{:?}; interchangeable leg classes in the native query, sewn-pair orbit representatives and multiplicities during preparation",
                    self.symmetrization
                ),
                PlanStage::Planning,
                PlanRole::Restriction,
            );
            if self.kind == GenerationType::CrossSection && self.symmetrization.initial {
                criterion(
                    "sewn-pair orbit",
                    "canonical orbit representative with exact orbit multiplicity".into(),
                    PlanStage::Preparation,
                    PlanRole::FinalCheck,
                );
            }
        }
        if self.kind == GenerationType::CrossSection {
            criterion(
                "physical cuts",
                format!(
                    "resolved final states={:?}; covariant partners included; blobs={:?}, spectators={:?}; retain a candidate only with a matching physical cut",
                    search
                        .resolved
                        .outgoing
                        .iter()
                        .map(|state| state
                            .iter()
                            .map(|particle| particle.pdg_code)
                            .collect::<Vec<_>>())
                        .collect::<Vec<_>>(),
                    self.cuts.blobs,
                    self.cuts.spectators
                ),
                PlanStage::Preparation,
                PlanRole::FinalCheck,
            );
            if self.cuts.amplitude != CutAmplitude::default() {
                criterion(
                    "cut amplitude",
                    format!(
                        "sum(left loops, right loops)={:?}; each side's coupling orders={:?}",
                        self.cuts.amplitude.loops, self.cuts.amplitude.orders.plain
                    ),
                    PlanStage::Preparation,
                    PlanRole::FinalCheck,
                );
            }
            if let Some(unresolved) = unresolved_cut_content(model, self)? {
                criterion(
                    "perturbative unresolved cut content",
                    format!(
                        "orders={:?}; additional massless particles={:?}; maximum multiplicity={}",
                        self.orders.perturbative,
                        unresolved.particles,
                        unresolved.maximum_multiplicity
                    ),
                    PlanStage::Preparation,
                    PlanRole::FinalCheck,
                );
            }
            if self.cuts.sewn_tadpoles {
                criterion(
                    "sewn tadpoles",
                    "exclude tadpole creation under initial-state sewing".into(),
                    PlanStage::Preparation,
                    PlanRole::FinalCheck,
                );
            }
            if self.cuts.symmetrize_left_right {
                criterion(
                    "left/right quotient",
                    "one global exchange of every sewn initial pair; canonical representatives merge with summed multiplicities".into(),
                    PlanStage::Preparation,
                    PlanRole::FinalCheck,
                );
            }
        } else if !self.orders.squared.is_empty() {
            criterion(
                "squared-order partner relation",
                format!(
                    "{:?}; whole-class final check after physical preparation: one common observed partner must satisfy every summed-order bound; retain the symmetric ordered relation through selection and grouping",
                    self.orders.squared
                ),
                PlanStage::Preparation,
                PlanRole::FinalCheck,
            );
        }
        criterion(
            "physical preparation",
            "native assignments still require canonicalization, CP-rule validation and fermion-flow/sign preparation; no numerator or final diagram is constructed by explain".into(),
            PlanStage::Preparation,
            PlanRole::FinalCheck,
        );
        Ok(report)
    }
    pub fn generate(&self, model: Arc<Model>, run: &Run) -> Result<DiagramSet, GenerationError> {
        let generator = Generator::new(model);
        let mut search = generator.search(self, |_, _| {}, |_, _| {})?;
        let mut diagrams: Vec<PreparedDiagram> = Vec::new();
        let mut firsts = HashMap::new();
        let (mut report, _) = search.visit(&generator, run, |diagram| {
            if search.resolved.identifies_mirror_images() {
                let index = diagrams.len();
                let first = *firsts.entry(OrbitKey::of(&diagram.input)).or_insert(index);
                if first != index {
                    let input = Arc::make_mut(&mut diagrams[first].input);
                    input.multiplicity = input
                        .multiplicity
                        .checked_add(diagram.input.multiplicity)
                        .ok_or(GenerationError::ArithmeticOverflow("orbit multiplicity"))?;
                    return Ok(());
                }
            }
            diagrams.push(diagram);
            Ok(())
        })?;
        let classes = OrderClasses::resolve(self, diagrams.iter().map(|d| d.orders.clone()));
        if let Some(classes) = &classes {
            diagrams.retain(|diagram| classes.retained(classes.class(&diagram.orders)));
        }
        let mut skeletons = Vec::new();
        let mut seen = HashSet::new();
        for diagram in &diagrams {
            if seen.insert(diagram.skeleton.storage_id()) {
                skeletons.push(diagram.skeleton.clone());
            }
        }
        report.retained_count = diagrams.len().into();
        run.report_progress(Stage::Done, diagrams.len(), Some(diagrams.len()));
        report.completed &= !run.cancellation_requested();
        // Results retain physical tables and model ownership, never an executable plan or
        // callback closure (which may itself own Python objects and execution-local state).
        search.plan = None;
        search.request.prune_skeletons = None;
        search.request.prune_dressings = None;
        search.request.select = None;
        search.request.deferred_hooks = None;
        Ok(DiagramSet {
            context: Arc::new(SetContext { generator, search }),
            diagrams,
            skeletons,
            report,
            classes,
        })
    }
    pub fn count(&self, model: Arc<Model>, run: &Run) -> Result<Counts, GenerationError> {
        let generator = Generator::new(model);
        let search = generator.search(self, |_, _| {}, |_, _| {})?;
        let mut by_class = BTreeMap::<BTreeMap<String, usize>, BTreeMap<usize, LoopCounts>>::new();
        let mut seen = HashSet::new();
        let (report, mut per_loop) = search.visit(&generator, run, |diagram| {
            let input = &diagram.input;
            let first =
                !search.resolved.identifies_mirror_images() || seen.insert(OrbitKey::of(input));
            let row = by_class
                .entry(diagram.orders.clone())
                .or_default()
                .entry(diagram.loops())
                .or_default();
            if first {
                row.diagrams += 1u8;
                row.sum_inv_aut += BigRational::new(BigInt::from(1), BigInt::from(diagram.aut()));
            }
            let negative = input.signs.internal_loop_negative
                ^ (input.signs.include_external_ordering && input.signs.external_ordering_negative)
                ^ (input.signs.include_antifermion_spin_sum
                    && input.signs.antifermion_spin_sum_negative);
            let numerator = BigInt::from(input.multiplicity) * if negative { -1 } else { 1 };
            row.physical_weight += BigRational::new(numerator, BigInt::from(diagram.aut()));
            Ok(())
        })?;
        let classes = OrderClasses::resolve(self, by_class.keys().cloned());
        for (orders, counts) in by_class {
            if classes
                .as_ref()
                .is_some_and(|classes| !classes.retained(classes.class(&orders)))
            {
                continue;
            }
            for (loops, row) in counts {
                let total = per_loop.entry(loops).or_default();
                total.diagrams += row.diagrams;
                total.sum_inv_aut += row.sum_inv_aut;
                total.physical_weight += row.physical_weight;
            }
        }
        let mut result = Counts {
            skeletons: report.topology_count,
            dressings: report.interaction_assignment_count,
            complete: report.completed,
            per_loop,
            ..Default::default()
        };
        for row in result.per_loop.values() {
            result.diagrams += &row.diagrams;
            result.sum_inv_aut += &row.sum_inv_aut;
            result.physical_weight += &row.physical_weight;
        }
        run.report_progress(Stage::Done, result.diagrams.clone(), None);
        result.complete &= !run.cancellation_requested();
        Ok(result)
    }
}

impl BuiltDiagramSet {
    /// Validate collapsed group provenance against the finalized master list.
    pub fn validate_groups(&self) -> Result<(), GenerationError> {
        if self.groups.len() != self.diagrams.len() {
            return Err(GenerationError::InvalidGroups(format!(
                "{} groups do not cover {} collapsed diagrams",
                self.groups.len(),
                self.diagrams.len()
            )));
        }
        let mut masters = BTreeSet::new();
        for (group_index, group) in self.groups.iter().enumerate() {
            if group.master >= self.diagrams.len() || !masters.insert(group.master) {
                return Err(GenerationError::InvalidGroups(format!(
                    "group {group_index} has invalid or duplicate master {}",
                    group.master
                )));
            }
            if group.members.is_empty()
                || group
                    .members
                    .iter()
                    .any(|member| member.diagram != group.master)
            {
                return Err(GenerationError::InvalidGroups(format!(
                    "group {group_index} has no members or a member assigned to another master"
                )));
            }
            self.diagrams[group.master].validate()?;
        }
        if masters != (0..self.diagrams.len()).collect::<BTreeSet<_>>() {
            return Err(GenerationError::InvalidGroups(
                "collapsed diagram masters are not covered exactly once".to_owned(),
            ));
        }
        Ok(())
    }
}

impl BuiltDiagramSet {
    /// Compare numerators only within a resolved squared-order class. Rebuild the requested
    /// comparison frame through the original canonical owner, then transport selected cuts and
    /// explicitly chosen loop edges through its actual edge permutation.
    ///
    /// Grouping consumes the set: [`Grouping::None`] returns it as it is, and a collapsed set
    /// moves the masters' owners and releases the physical preparations it was grouped from.
    pub fn group(mut self, grouping: &Grouping, run: &Run) -> Result<Self, GenerationError> {
        if matches!(grouping, Grouping::None) {
            return Ok(self.done(run));
        }
        let Some(prepared) = self.prepared.clone() else {
            return Err(GenerationError::InvalidRequest("group an ungrouped built set; grouping a collapsed set would lose its source terms".into()));
        };
        // Canonical comparison alone is an exact grouping mode without samples.
        if matches!(
            grouping,
            Grouping::Identical(options) | Grouping::UpToSign(options) | Grouping::UpToScalar(options)
                if options.number_of_numerical_samples == 0 && !options.check_canonical_numerator
        ) {
            return Err(GenerationError::InvalidNumericalSampleCount);
        }
        let context = &prepared.context;
        let mut order = (0..self.len()).collect::<Vec<_>>();
        // This is the former pre-preparation ordering boundary, demanded only by grouping.
        order.sort_by(|&a, &b| {
            prepared.diagrams[self.sources[a]]
                .compact
                .canonical_form()
                .cmp(prepared.diagrams[self.sources[b]].compact.canonical_form())
        });
        let prepare = |&index: &usize| -> Result<Option<DiagramPair>, GenerationError> {
            if run.cancellation_requested() {
                return Ok(None);
            }
            let source = &prepared.diagrams[self.sources[index]];
            let raw = context
                .search
                .externals
                .graph(&source.compact, &context.search.signatures);
            let input = context
                .generator
                .diagram_input(
                    &context.search.resolved,
                    raw,
                    source.aut(),
                    source.input.multiplicity,
                    &context.search.request,
                    grouping,
                )?
                .ok_or_else(|| {
                    GenerationError::InvalidGroups(
                        "grouping changed physical preparation acceptance".into(),
                    )
                })?;
            // Both preparations reframe the same assignment: move the requested factors from
            // the ungrouped physical frame into the comparison and representative frames.
            let ungrouped = &source.input.representative.graph;
            let comparison_options = ungrouped
                .frame_map(&input.comparison.graph)
                .options(&self.options)?;
            let representative_frame = ungrouped.frame_map(&input.representative.graph);
            let representative_options = representative_frame.options(&self.options)?;
            let new_internal = input
                .representative
                .graph
                .edges()
                .enumerate()
                .filter_map(|(edge, data)| (!data.is_external()).then_some(edge))
                .collect::<Vec<_>>();
            let mut pair = context.generator.diagram_pair(
                self.diagrams[index].name().map(str::to_owned),
                context.search.request.kind,
                input,
                (&comparison_options, &representative_options),
            )?;
            let old_internal = ungrouped
                .edges()
                .enumerate()
                .filter_map(|(edge, data)| (!data.is_external()).then_some(edge));
            let old_final = self.diagrams[index]
                .edges()
                .filter_map(|(edge, _, data)| data.external.is_none().then_some(edge));
            let new_final = pair
                .representative()
                .edges()
                .filter_map(|(edge, _, data)| data.external.is_none().then_some(edge))
                .collect::<Vec<_>>();
            let new_ids = new_internal
                .into_iter()
                .zip(new_final)
                .collect::<BTreeMap<_, _>>();
            let map = old_internal
                .zip(old_final)
                .map(|(edge, final_edge)| {
                    (final_edge, new_ids[&representative_frame.edges[edge].0])
                })
                .collect::<BTreeMap<_, _>>();
            let map_edge = |edge: &EdgeId| {
                map.get(edge).copied().ok_or_else(|| {
                    GenerationError::InvalidGroups(format!(
                        "cut or loop edge {} is not internal",
                        edge.0
                    ))
                })
            };
            let selected = self.diagrams[index]
                .cuts()
                .iter()
                .map(|cut| {
                    cut.cut
                        .iter()
                        .map(|half_edge| map_edge(&half_edge.edge))
                        .collect::<Result<BTreeSet<_>, _>>()
                })
                .collect::<Result<BTreeSet<_>, _>>()?;
            let cuts = pair
                .representative()
                .cuts()
                .iter()
                .filter(|cut| {
                    selected.contains(&cut.cut.iter().map(|h| h.edge).collect::<BTreeSet<_>>())
                })
                .cloned()
                .collect::<Vec<_>>();
            let found = cuts
                .iter()
                .map(|cut| cut.cut.iter().map(|h| h.edge).collect::<BTreeSet<_>>())
                .collect::<BTreeSet<_>>();
            if found != selected {
                return Err(GenerationError::InvalidGroups(
                    "canonical frame lost a selected physical cut".into(),
                ));
            }
            if cuts != pair.representative().cuts() {
                pair = pair.with_cuts(cuts)?;
            }
            let named = self.diagrams[index]
                .name()
                .is_some_and(|name| self.options.named_loop_momentum_bases.contains_key(name));
            let by_id = !self.options.loop_momentum_bases.is_empty()
                && self
                    .options
                    .loop_momentum_bases
                    .contains_key(&self.diagrams[index].id()?);
            if named || by_id {
                let basis = self.diagrams[index]
                    .loop_momentum_basis()
                    .loop_edges
                    .iter()
                    .map(map_edge)
                    .collect::<Result<Vec<_>, _>>()?;
                pair = pair.with_loop_momentum_edges(&basis)?;
            }
            Ok(Some(pair))
        };
        let pairs = run
            .in_pool(|| run.map_ordered(order.iter().collect(), prepare))?
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?;
        if run.cancellation_requested() {
            self.report.completed = false;
            return Ok(self.done(run));
        }
        let classes = order.iter().map(|&i| self.classes[i]).collect::<Vec<_>>();
        let source_indices = order
            .iter()
            .map(|&index| self.sources[index])
            .collect::<Vec<_>>();
        let source_lookup = source_indices
            .iter()
            .copied()
            .enumerate()
            .map(|(local, source)| (source, local))
            .collect::<BTreeMap<_, _>>();
        let pairs = pairs
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .expect("an uncancelled preparation builds every pair");
        // The pairs replace the ungrouped diagrams. Members keep only their sources' identities,
        // whose failures surface for a member, as when they were read from the source itself.
        let mut identities = std::mem::take(&mut self.diagrams)
            .iter()
            .map(source_identity)
            .collect::<Vec<_>>();
        let grouped = run.in_pool(|| {
            grouping::group_diagrams(
                pairs,
                &context.generator.model,
                grouping,
                context.search.resolved.identifies_mirror_images(),
                run,
                false,
                DiagramLabels {
                    classes: self.relation.as_ref().map(|_| classes.as_slice()),
                    source_indices: Some(&source_indices),
                },
            )
        })??;
        let (ungrouped_sources, ungrouped_classes) = (
            std::mem::take(&mut self.sources),
            std::mem::take(&mut self.classes),
        );
        self.groups = grouped.groups;
        self.prepared = None;
        self.report.zero_numerator_count += BigUint::from(grouped.zero_numerator_count);
        // Grouping moved each master's physical owner to its output, carrying the merged
        // partitions, routing and factor of its comparison frame.
        for (output, (representative, group)) in grouped
            .diagrams
            .into_iter()
            .zip(&mut self.groups)
            .enumerate()
        {
            let first = group
                .members
                .first()
                .ok_or_else(|| {
                    GenerationError::InvalidGroups("group has no representative".into())
                })?
                .source_diagram;
            let original = order[source_lookup[&first]];
            self.diagrams.push(representative);
            self.sources.push(ungrouped_sources[original]);
            self.classes.push(ungrouped_classes[original]);
            group.master = output;
            for member in &mut group.members {
                let original = order[source_lookup[&member.source_diagram]];
                member.source_identity = std::mem::replace(&mut identities[original], Ok(None))?;
                member.source_diagram = ungrouped_sources[original];
                member.diagram = output;
            }
        }
        self.report.retained_count = self.len().into();
        let result = self.done(run);
        result.validate_groups()?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::super::cases::Case;
    use super::super::cp_merge_diagnostic;
    use super::super::{InteractionStars, ResolvedProcess};
    use super::*;
    use crate::Classes;
    use crate::generation::graph::{Endpoint, GeneratedGraph};
    use crate::{Cuts, GenerationControl, GraphGroupingOptions};
    use crate::{Loops, Orders, Process, Restrict, TopologyRules};
    use feynkit_graph as graph_api;
    use feynkit_graph::symbols;
    use feynkit_graph::{DiagramEndpoint, DiagramId, ExternalState};
    use serde_json::Value;
    use serde_json::json;
    use std::env;
    use std::fs;
    use std::mem;
    use std::path::PathBuf;
    use symbolica::atom::AtomCore;

    /// Four independently sewn tree/one-loop QED sources (domain 0, 78, 132, 155).
    /// The external photons attach to vertices 0..4 in exact label order. Neither
    /// these arrows nor the B/M actions below use the preparation canonicalizer.
    fn beam_cp_sources(model: &Model) -> [GeneratedGraph; 4] {
        let photon = model.particle_id_by_pdg(22).unwrap();
        let electron = model.particle_id_by_pdg(11).unwrap();
        let rule = model.vertex_rule_id("V_98").unwrap();
        [
            [(0, 2), (1, 0), (2, 4), (3, 1), (4, 5), (5, 3)],
            [(0, 1), (1, 3), (2, 0), (3, 4), (4, 5), (5, 2)],
            [(0, 2), (1, 4), (2, 3), (3, 1), (4, 5), (5, 0)],
            [(0, 4), (1, 3), (2, 0), (3, 2), (4, 5), (5, 1)],
        ]
        .map(|arrows| {
            GeneratedGraph::new(
                [rule; 6],
                [GeneratedEdge {
                    source: Endpoint::Vertex(4),
                    sink: Endpoint::Vertex(5),
                    particle: photon,
                    directed: false,
                }]
                .into_iter()
                .chain(arrows.map(|(source, sink)| GeneratedEdge {
                    source: Endpoint::Vertex(source),
                    sink: Endpoint::Vertex(sink),
                    particle: electron,
                    directed: true,
                }))
                .chain((0..4).map(|index| GeneratedEdge {
                    source: Endpoint::External(ExternalNode {
                        index,
                        state: if index < 2 {
                            ExternalState::Incoming
                        } else {
                            ExternalState::Outgoing
                        },
                        particle: photon,
                    }),
                    sink: Endpoint::Vertex(index),
                    particle: photon,
                    directed: false,
                })),
            )
        })
    }

    #[test]
    fn joint_beam_cp_orbit_preserves_complete_preparation_and_source_coverage() {
        let model = Arc::new(Model::qed());
        let generator = Generator::new(Arc::clone(&model));
        let graphs = beam_cp_sources(&model);
        let physical = |graph: &GeneratedGraph| {
            // A storage reversal with charge conjugation preserves a directed line.
            // This comparison retains physical arrows, exact external states and rules.
            graph
                .reframed(
                    &(0..graph.n_vertices()).collect::<Vec<_>>(),
                    graph.edges().enumerate().map(|(edge, data)| {
                        let particle = model.particle_by_id(data.particle).unwrap();
                        let reverse = data.directed && particle.is_antiparticle();
                        (
                            edge,
                            reverse,
                            if reverse {
                                particle.antiparticle
                            } else {
                                data.particle
                            },
                        )
                    }),
                    |external| external,
                )
                .to_hedge_graph()
                .canonical_form_by(|&rule| rule, |&particle| particle, |&external| external)
        };
        // Explicit independent B=(01)(23), M=CP, BM maps. CP conjugates a
        // charged line; storing it again as an electron reverses its endpoints.
        for (images, vertices, beam, mirror) in [
            ([1, 0, 3, 2], [1, 0, 3, 2, 4, 5], true, false),
            ([2, 3, 0, 1], [2, 3, 0, 1, 5, 4], false, true),
            ([3, 2, 1, 0], [3, 2, 1, 0, 5, 4], true, true),
        ] {
            for (source, &image) in images.iter().enumerate() {
                let transformed = graphs[source].reframed(
                    &vertices,
                    graphs[source]
                        .edges()
                        .enumerate()
                        .map(|(edge, data)| (edge, mirror && data.directed, data.particle)),
                    |mut external| {
                        if beam {
                            external.index ^= 1;
                        }
                        if mirror {
                            external.index = (external.index + 2) % 4;
                            external.state = match external.state {
                                ExternalState::Incoming => ExternalState::Outgoing,
                                ExternalState::Outgoing => ExternalState::Incoming,
                            };
                        }
                        external
                    },
                );
                assert_eq!(physical(&transformed), physical(&graphs[image]));
            }
        }
        let process = Process::new(["a", "a"], ["e+", "e-"]);
        let run = Run::default().threads(1);
        for classes in [Classes::None, Classes::Auto] {
            for (beam, mirror, partition) in [
                (false, false, [0, 1, 2, 3]),
                (true, false, [0, 0, 2, 2]),
                (false, true, [0, 1, 0, 1]),
                (true, true, [0, 0, 0, 0]),
            ] {
                let request = process
                    .clone()
                    .cross_section(Loops::forward(2..=2))
                    .with_classes(classes)
                    .with_topology(TopologyRules::default().with_self_loops(false))
                    .with_symmetrize(Symmetrize {
                        initial: beam,
                        final_state: false,
                        external_fermions: false,
                    })
                    .with_cuts(
                        Cuts::new()
                            .with_blobs(1..=1)
                            .with_spectators(0..=0)
                            .with_symmetrize_left_right(mirror),
                    );
                let search = generator.search(&request, |_, _| {}, |_, _| {}).unwrap();
                let mut members = Vec::new();
                search
                    .enumerate(&run, |batch| {
                        for compact in batch {
                            let raw = search.externals.graph(&compact, &search.signatures);
                            if let Some(source) = graphs
                                .iter()
                                .position(|graph| physical(graph) == physical(&raw))
                            {
                                members.push((source, compact));
                            }
                        }
                        Ok(())
                    })
                    .unwrap();
                members.sort_by_key(|(source, _)| *source);
                assert_eq!(
                    members
                        .iter()
                        .map(|(source, _)| *source)
                        .collect::<Vec<_>>(),
                    [0, 1, 2, 3]
                );
                let selected = members
                    .iter()
                    .map(|(_, compact)| compact.canonical_form().clone())
                    .collect::<Vec<_>>();
                let request = request.select(move |view| {
                    selected.contains(&view.complete_graph().unwrap().canonical_form_by(
                        |&rule| rule,
                        |&particle| particle,
                        |&label| label,
                    ))
                });
                let mut inputs = Vec::new();
                let mut events = Vec::new();
                for (source, compact) in members {
                    cp_merge_diagnostic::start();
                    let prepared = search.prepare(&generator, compact).unwrap();
                    events.push(json!({"source": source, "events": cp_merge_diagnostic::finish()}));
                    if let Some(prepared) = prepared {
                        let input = prepared.input;
                        assert_eq!(input.representative.symmetry, 1);
                        assert_eq!(input.representative.cut_partitions.len(), 1);
                        assert_eq!(input.representative.topology_threshold_candidates.len(), 1);
                        assert!(input.signs.internal_loop_negative);
                        assert!(input.signs.include_external_ordering);
                        assert!(!input.signs.external_ordering_negative);
                        assert!(input.signs.include_antifermion_spin_sum);
                        assert!(!input.signs.antifermion_spin_sum_negative);
                        assert_eq!(input.comparison, input.representative);
                        assert!(input.reversed_edges.is_empty());
                        let graph = &input.representative.graph;
                        let canonical = search
                            .resolved
                            .canonicalize_numerator_graph(
                                graph,
                                &model,
                                &generator.cache.stars,
                                &Grouping::None,
                            )
                            .unwrap();
                        assert_eq!(
                            canonical,
                            search
                                .resolved
                                .canonicalize_numerator_graph(
                                    &canonical,
                                    &model,
                                    &generator.cache.stars,
                                    &Grouping::None
                                )
                                .unwrap()
                        );
                        let repeated = generator
                            .diagram_input(
                                &search.resolved,
                                graph.clone(),
                                1,
                                input.multiplicity,
                                &request,
                                &Grouping::None,
                            )
                            .unwrap()
                            .unwrap();
                        assert_eq!(input.representative, repeated.representative);
                        assert_eq!(input.signs, repeated.signs);
                        inputs.push((source, input));
                    }
                }
                assert_eq!(inputs.len(), if beam { 2 } else { 4 });
                assert_eq!(
                    inputs
                        .iter()
                        .map(|(_, input)| input.multiplicity)
                        .sum::<u64>(),
                    4
                );
                // Admission must cover each B pair once, with its complete multiplicity.
                for pair in [[0, 1], [2, 3]] {
                    assert_eq!(
                        inputs
                            .iter()
                            .filter(|(source, _)| pair.contains(source))
                            .count(),
                        if beam { 1 } else { 2 }
                    );
                }
                eprintln!(
                    "joint preparation beam={beam} cp={mirror} classes={classes:?}: {}",
                    json!(events)
                );
                for (a, left) in &inputs {
                    for (b, right) in &inputs {
                        assert_eq!(
                            OrbitKey::of(left) == OrbitKey::of(right),
                            partition[*a] == partition[*b],
                            "complete joint key: beam={beam} cp={mirror} classes={classes:?}, sources {a}/{b}"
                        );
                    }
                }
                let set = request.generate(Arc::clone(&model), &run).unwrap();
                let counts = request.count(Arc::clone(&model), &run).unwrap();
                let expected = if beam && mirror {
                    1
                } else if beam || mirror {
                    2
                } else {
                    4
                };
                assert_eq!(set.len(), expected);
                assert_eq!(counts.diagrams, BigUint::from(expected));
                assert!(set.report.completed && counts.complete);
                assert_eq!(
                    counts.sum_inv_aut,
                    BigRational::from_integer(expected.into())
                );
                assert_eq!(
                    counts.physical_weight,
                    BigRational::from_integer((-4).into())
                );
                assert_eq!(
                    set.diagrams
                        .iter()
                        .map(|diagram| diagram.input.multiplicity)
                        .sum::<u64>(),
                    4
                );
                for diagram in &set.diagrams {
                    assert!(
                        inputs
                            .iter()
                            .any(|(_, input)| OrbitKey::of(input) == OrbitKey::of(&diagram.input))
                    );
                }
            }
        }
    }

    #[test]
    fn cp_merge_preserves_complete_preparation_and_counts_each_graph_once() {
        let model = Arc::new(Model::qed());
        let photon = model.particle_id_by_pdg(22).unwrap();
        let electron = model.particle_id_by_pdg(11).unwrap();
        let rule = model.vertex_rule_id("V_98").unwrap();
        // Two charged triangles joined by a photon. The four labelled photon legs
        // are the unsewn form of the witnessed aaee forward graph.
        let graph = GeneratedGraph::new(
            [rule; 6],
            [
                (0, 1, photon),
                (3, 0, electron),
                (0, 5, electron),
                (1, 2, electron),
                (4, 1, electron),
                (2, 4, electron),
                (5, 3, electron),
            ]
            .into_iter()
            .map(|(source, sink, particle)| GeneratedEdge {
                source: Endpoint::Vertex(source),
                sink: Endpoint::Vertex(sink),
                particle,
                directed: particle == electron,
            })
            .chain(
                [(0, 5), (1, 3), (2, 2), (3, 4)].map(|(index, vertex)| GeneratedEdge {
                    source: Endpoint::External(ExternalNode {
                        index,
                        state: if index < 2 {
                            ExternalState::Incoming
                        } else {
                            ExternalState::Outgoing
                        },
                        particle: photon,
                    }),
                    sink: Endpoint::Vertex(vertex),
                    particle: photon,
                    directed: false,
                }),
            ),
        );
        // Explicit CP, independent of the canonicalizer: exchange all sewn legs,
        // conjugate charged fields, and apply the certified vertex permutation.
        let cp = |graph: &GeneratedGraph| {
            let endpoint = |endpoint| match endpoint {
                Endpoint::Vertex(vertex) => Endpoint::Vertex([1, 0, 5, 4, 3, 2][vertex]),
                Endpoint::External(external) => Endpoint::External(ExternalNode {
                    index: (external.index + 2) % 4,
                    state: if external.index < 2 {
                        ExternalState::Outgoing
                    } else {
                        ExternalState::Incoming
                    },
                    ..external
                }),
            };
            GeneratedGraph::new(
                graph.vertices(),
                graph.edges().map(|edge| GeneratedEdge {
                    source: endpoint(edge.source),
                    sink: endpoint(edge.sink),
                    particle: if edge.directed {
                        model.particle_by_id(edge.particle).unwrap().antiparticle
                    } else {
                        edge.particle
                    },
                    ..edge
                }),
            )
        };
        let mirror = cp(&graph);
        assert_eq!(cp(&mirror), graph);
        let process = Process::new(["a", "a"], ["e+", "e-"]);
        let generator = Generator::new(Arc::clone(&model));
        let request = |enabled| {
            process
                .clone()
                .cross_section(Loops::forward(2..=2))
                .with_topology(TopologyRules::new().with_self_loops(false))
                .with_cuts(
                    Cuts::new()
                        .with_blobs(1..=1)
                        .with_spectators(0..=0)
                        .with_symmetrize_left_right(enabled),
                )
        };
        let off = request(false);
        let resolved =
            ResolvedProcess::new(&model, &process, GenerationType::CrossSection, &off).unwrap();
        let mut inventory = resolved.cut_partitions(&graph, &model, &off).unwrap();
        let mut other = resolved.cut_partitions(&mirror, &model, &off).unwrap();
        assert_eq!(inventory.physical.len(), 2);
        assert_eq!(inventory.topology_threshold_candidates.len(), 2);
        // CP interchanges the complete sides, including orders and loop counts.
        // The selected cut end becomes its complement; storage itself is unchanged.
        let complement = |endpoint| match endpoint {
            DiagramEndpoint::Source => DiagramEndpoint::Target,
            DiagramEndpoint::Target => DiagramEndpoint::Source,
        };
        for cut in &mut inventory.physical {
            mem::swap(&mut cut.left, &mut cut.right);
            for hedge in &mut cut.cut {
                hedge.endpoint = complement(hedge.endpoint);
            }
            cut.cut.sort();
        }
        for threshold in &mut inventory.topology_threshold_candidates {
            mem::swap(&mut threshold.left, &mut threshold.right);
            for hedge in &mut threshold.cut {
                hedge.endpoint = complement(hedge.endpoint);
            }
            threshold.cut.sort();
        }
        inventory.physical.sort();
        inventory.topology_threshold_candidates.sort();
        other.physical.sort();
        other.topology_threshold_candidates.sort();
        assert_eq!(inventory.physical, other.physical);
        assert_eq!(
            inventory.topology_threshold_candidates,
            other.topology_threshold_candidates
        );

        for classes in [Classes::None, Classes::Auto] {
            for enabled in [false, true] {
                let request = request(enabled).with_classes(classes);
                let resolved =
                    ResolvedProcess::new(&model, &process, GenerationType::CrossSection, &request)
                        .unwrap();
                let prepared = [&graph, &mirror].map(|graph| {
                    let canonical = resolved
                        .canonicalize_numerator_graph(
                            graph,
                            &model,
                            &InteractionStars::new(&model),
                            &Grouping::None,
                        )
                        .unwrap();
                    assert_eq!(
                        canonical,
                        resolved
                            .canonicalize_numerator_graph(
                                &canonical,
                                &model,
                                &InteractionStars::new(&model),
                                &Grouping::None
                            )
                            .unwrap()
                    );
                    let input = generator
                        .diagram_input(&resolved, graph.clone(), 1, 1, &request, &Grouping::None)
                        .unwrap()
                        .unwrap();
                    assert_eq!(input.representative.cut_partitions.len(), 2);
                    assert_eq!(input.representative.topology_threshold_candidates.len(), 2);
                    assert_eq!(input.representative.symmetry, 1);
                    assert_eq!(
                        input.signs,
                        FermionSigns {
                            internal_loop_negative: false,
                            include_external_ordering: true,
                            external_ordering_negative: false,
                            include_antifermion_spin_sum: true,
                            antifermion_spin_sum_negative: false,
                        }
                    );
                    let repeated = generator
                        .diagram_input(
                            &resolved,
                            input.representative.graph.clone(),
                            1,
                            1,
                            &request,
                            &Grouping::None,
                        )
                        .unwrap()
                        .unwrap();
                    assert_eq!(input.representative, repeated.representative);
                    assert_eq!(input.signs, repeated.signs);
                    assert!(input.reversed_edges.is_empty());
                    assert_eq!(input.comparison, input.representative);
                    input
                });
                assert_eq!(
                    OrbitKey::of(&prepared[0]) == OrbitKey::of(&prepared[1]),
                    enabled,
                    "CP={enabled}, classes={classes:?}: compare every ordered OrbitKey field"
                );
                if enabled {
                    assert_eq!(prepared[0].representative, prepared[1].representative);
                }
                let set = request
                    .generate(Arc::clone(&model), &Run::default().threads(1))
                    .unwrap();
                let witnesses = set
                    .diagrams
                    .iter()
                    .filter(|diagram| {
                        prepared
                            .iter()
                            .any(|input| OrbitKey::of(&diagram.input) == OrbitKey::of(input))
                    })
                    .collect::<Vec<_>>();
                assert_eq!(witnesses.len(), if enabled { 1 } else { 2 });
                assert_eq!(
                    witnesses
                        .iter()
                        .map(|diagram| diagram.input.multiplicity)
                        .sum::<u64>(),
                    2,
                    "two physical cuts must not multiply the number of source graphs"
                );
            }
        }
    }

    #[test]
    fn cp_rule_colors_are_repaired_before_frame_selection() {
        let model = Arc::new(Model::standard_model());
        let rules = ["V_82", "V_117", "V_133", "V_133", "V_71"]
            .map(|name| model.vertex_rule_id(name).unwrap());
        // Two flavour-changing quark pentagons. The two labelled charm-photon
        // attachments in each pentagon distinguish its two physical orientations.
        // Goldstones join the conjugate weak vertices; the remaining photons join
        // the down-quark vertices. All stars are actual Standard Model rules.
        let graph = GeneratedGraph::new(
            rules.into_iter().chain(rules),
            [
                (0, 2, 4),
                (2, 3, 4),
                (3, 1, 4),
                (1, 4, 1),
                (4, 0, 1),
                (5, 7, 4),
                (7, 8, 4),
                (8, 6, 4),
                (6, 9, 1),
                (9, 5, 1),
                (6, 0, 251),
                (1, 5, 251),
                (4, 9, 22),
            ]
            .map(|(source, sink, pdg)| GeneratedEdge {
                source: Endpoint::Vertex(source),
                sink: Endpoint::Vertex(sink),
                particle: model.particle_id_by_pdg(pdg).unwrap(),
                directed: pdg != 22,
            })
            .into_iter()
            .chain(
                [(0, 2), (1, 3), (2, 7), (3, 8)].map(|(index, vertex)| GeneratedEdge {
                    source: Endpoint::External(ExternalNode {
                        index,
                        state: if index < 2 {
                            ExternalState::Incoming
                        } else {
                            ExternalState::Outgoing
                        },
                        particle: model.particle_id_by_pdg(22).unwrap(),
                    }),
                    sink: Endpoint::Vertex(vertex),
                    particle: model.particle_id_by_pdg(22).unwrap(),
                    directed: false,
                }),
            ),
        );
        // An explicit full CP map: exchange pentagons and sewn partners,
        // conjugate every charged field AND its interaction rule.
        let cp = |graph: &GeneratedGraph| {
            let mut cp_rules = vec![rules[0]; 10];
            for (vertex, rule) in graph.vertices().enumerate() {
                cp_rules[(vertex + 5) % 10] = if rule == rules[0] {
                    rules[1]
                } else if rule == rules[1] {
                    rules[0]
                } else {
                    rule
                };
            }
            let endpoint = |endpoint| match endpoint {
                Endpoint::Vertex(vertex) => Endpoint::Vertex((vertex + 5) % 10),
                Endpoint::External(external) => Endpoint::External(ExternalNode {
                    index: (external.index + 2) % 4,
                    state: if external.index < 2 {
                        ExternalState::Outgoing
                    } else {
                        ExternalState::Incoming
                    },
                    ..external
                }),
            };
            GeneratedGraph::new(
                cp_rules,
                graph.edges().map(|edge| GeneratedEdge {
                    source: endpoint(edge.source),
                    sink: endpoint(edge.sink),
                    particle: model.particle_by_id(edge.particle).unwrap().antiparticle,
                    ..edge
                }),
            )
        };
        let mirror = cp(&graph);
        assert_eq!(cp(&mirror), graph);
        let generator = Generator::new(Arc::clone(&model));
        for graph in [&graph, &mirror] {
            assert_eq!(
                *graph,
                generator
                    .repair_cp_interaction_rules(graph.clone())
                    .unwrap()
            );
            for (vertex, rule) in graph.vertices().enumerate() {
                generator
                    .interaction_half_edges(graph, vertex, model.vertex_rule_by_id(rule).unwrap())
                    .unwrap();
            }
        }
        let process = Process::new(["a", "a"], ["c", "c~"]);
        let request = process.clone().cross_section(Loops::forward(0..=6));
        let resolved =
            ResolvedProcess::new(&model, &process, GenerationType::CrossSection, &request).unwrap();
        let mut inventory = resolved.cut_partitions(&graph, &model, &request).unwrap();
        let mut other = resolved.cut_partitions(&mirror, &model, &request).unwrap();
        let complement = |endpoint| match endpoint {
            DiagramEndpoint::Source => DiagramEndpoint::Target,
            DiagramEndpoint::Target => DiagramEndpoint::Source,
        };
        for cut in &mut inventory.physical {
            mem::swap(&mut cut.left, &mut cut.right);
            for hedge in &mut cut.cut {
                hedge.endpoint = complement(hedge.endpoint);
            }
            cut.cut.sort();
        }
        for threshold in &mut inventory.topology_threshold_candidates {
            mem::swap(&mut threshold.left, &mut threshold.right);
            for hedge in &mut threshold.cut {
                hedge.endpoint = complement(hedge.endpoint);
            }
            threshold.cut.sort();
        }
        inventory.physical.sort();
        inventory.topology_threshold_candidates.sort();
        other.physical.sort();
        other.topology_threshold_candidates.sort();
        assert_eq!(inventory.physical, other.physical);
        assert_eq!(
            inventory.topology_threshold_candidates,
            other.topology_threshold_candidates
        );
        for enabled in [false, true] {
            if enabled {
                cp_merge_diagnostic::start();
            }
            let request = process
                .clone()
                .cross_section(Loops::forward(0..=6))
                .with_cuts(Cuts::new().with_symmetrize_left_right(enabled));
            let resolved =
                ResolvedProcess::new(&model, &process, GenerationType::CrossSection, &request)
                    .unwrap();
            let prepared = [&graph, &mirror].map(|graph| {
                let canonical = resolved
                    .canonicalize_numerator_graph(
                        graph,
                        &model,
                        &InteractionStars::new(&model),
                        &Grouping::None,
                    )
                    .unwrap();
                assert_eq!(
                    canonical,
                    generator
                        .repair_cp_interaction_rules(canonical.clone())
                        .unwrap(),
                    "CP frame selection must expose actual repaired rule colors"
                );
                assert_eq!(
                    canonical,
                    resolved
                        .canonicalize_numerator_graph(
                            &canonical,
                            &model,
                            &InteractionStars::new(&model),
                            &Grouping::None
                        )
                        .unwrap()
                );
                let input = generator
                    .diagram_input(&resolved, graph.clone(), 1, 1, &request, &Grouping::None)
                    .unwrap()
                    .unwrap();
                assert!(!input.representative.cut_partitions.is_empty());
                assert!(
                    !input
                        .representative
                        .topology_threshold_candidates
                        .is_empty()
                );
                assert_eq!(
                    input.representative.graph,
                    generator
                        .repair_cp_interaction_rules(input.representative.graph.clone())
                        .unwrap()
                );
                let repeated = generator
                    .diagram_input(
                        &resolved,
                        input.representative.graph.clone(),
                        1,
                        1,
                        &request,
                        &Grouping::None,
                    )
                    .unwrap()
                    .unwrap();
                assert_eq!(input.representative, repeated.representative);
                assert_eq!(input.signs, repeated.signs);
                input
            });
            assert_eq!(
                OrbitKey::of(&prepared[0]) == OrbitKey::of(&prepared[1]),
                enabled
            );
            if enabled {
                let events = cp_merge_diagnostic::finish();
                let choices = events
                    .iter()
                    .filter(|event| event["stage"] == "frame_choice" && event["remapping"] == true)
                    .collect::<Vec<_>>();
                assert!(!choices.is_empty());
                assert!(choices.iter().all(|event| event["tie"] == true));
            }
        }
    }

    #[test]
    #[ignore = "one focused scalar observation; requires a new FEYNKIT_ORACLE_DUMP directory"]
    fn scalar_preparation_and_group_components() {
        let directory = PathBuf::from(env::var("FEYNKIT_ORACLE_DUMP").unwrap());
        assert!(!directory.exists(), "diagnostic output must be new");
        let case = Case::named("xs_qcdqed_addx_canonical_scalar_lr_L3");
        let model = Arc::new(case.model());
        let request = case.request(&model).with_classes(Classes::None);
        let run = Run::default().threads(1);
        cp_merge_diagnostic::start();
        let prepared = request.generate(Arc::clone(&model), &run).unwrap();
        let prepared_sources = prepared.diagrams.iter().enumerate().map(|(index, diagram)| {
            let context = &prepared.context;
            let raw = context.search.externals.graph(&diagram.compact, &context.search.signatures);
            json!({"index": index, "raw": cp_merge_diagnostic::graph(&raw), "input": cp_merge_diagnostic::input(&diagram.input)})
        }).collect::<Vec<_>>();
        let options = case.numerators().with_naming(Naming::prefix("D"));
        let built = prepared.build(&options, &run).unwrap();
        // Operator observations occur at each actual owner before zero admission;
        // the final grouping records carry indices, factors and ratios only.
        let sources = built.sources.clone();
        let grouped = built.group(&case.grouping(), &run).unwrap();
        let groups = grouped.groups().iter().map(|group| json!({
            "master": group.master,
            "members": group.members.iter().map(|member| json!({
                "source": member.source_diagram, "ratio": member.ratio.to_canonical_string(),
                "factor": member.overall_factor.to_canonical_string(),
            })).collect::<Vec<_>>(),
            "factor": grouped.diagrams()[group.master].overall_factor().to_canonical_string(),
        })).collect::<Vec<_>>();
        let events = cp_merge_diagnostic::finish();
        fs::create_dir_all(&directory).unwrap();
        fs::write(
            directory.join("diagnostic.json"),
            serde_json::to_vec(&json!({
                "case": case.name, "classes": "None", "pool": 1,
                "events": events, "prepared_sources": prepared_sources,
                "sources": sources, "groups": groups,
            }))
            .unwrap(),
        )
        .unwrap();
    }

    fn model() -> Arc<Model> {
        Arc::new(
            Model::from_json(include_str!(
                "../../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
            ))
            .unwrap(),
        )
    }
    fn request() -> Request {
        Process::new(["scalar_0", "scalar_0"], ["scalar_0", "scalar_0"])
            .amplitude(Loops::total(0..=0))
            .with_restrict(Restrict::new().only_vertices(["V_3_SCALAR_000"]))
            .with_topology(TopologyRules::default().with_self_loops(false))
    }

    #[test]
    fn compact_generation_and_count_share_physical_preparation_without_numerators() {
        let model = model();
        let run = Run::default().threads(1);
        let request = request();
        let counts = request.count(model.clone(), &run).unwrap();
        let set = request.generate(model, &run).unwrap();
        assert!(counts.complete && set.report().completed);
        assert_eq!(counts.diagrams, BigUint::from(set.len()));
        assert_eq!(counts.dressings, set.report().interaction_assignment_count);
        assert_eq!(counts.skeletons, set.report().topology_count);
        assert_eq!(set.len(), 3);
        assert_eq!(counts.sum_inv_aut, BigRational::from_integer(3.into()));
        assert_eq!(counts.physical_weight, counts.sum_inv_aut);
        let cache = &set.context.generator.cache;
        let numerators = &set.context.generator.numerators;
        assert!(numerators.interactions.read().unwrap().is_empty());
        assert!(numerators.propagators.read().unwrap().is_empty());
        assert!(numerators.interaction_templates.read().unwrap().is_empty());
        assert!(numerators.propagator_templates.read().unwrap().is_empty());
        assert!(cache.projectors.read().unwrap().is_empty());
        let cloned = set.clone();
        assert!(Arc::ptr_eq(
            &set.diagrams[0].input,
            &cloned.diagrams[0].input
        ));
        for diagram in &set.diagrams {
            assert!(
                set.skeletons()
                    .iter()
                    .any(|s| s.storage_id() == diagram.skeleton().storage_id())
            );
        }
        let built = set.build(&NumeratorOptions::default(), &run).unwrap();
        assert_eq!(built.len(), set.len());
        assert!(built.iter().all(|d| d.name().is_none()));
        built.validate_groups().unwrap();
    }

    #[test]
    fn empty_and_cancelled_results_report_completion_honestly() {
        let model = model();
        let empty =
            request().with_restrict(Restrict::new().only_particles(std::iter::empty::<&str>()));
        let stages = Arc::new(Mutex::new(Vec::new()));
        let captured = stages.clone();
        let run = Run::default().threads(1).progress(move |p| {
            captured.lock().unwrap().push(p.stage);
            GenerationControl::Continue
        });
        let count = empty.count(model.clone(), &run).unwrap();
        assert!(count.complete);
        assert_eq!(count.diagrams, BigUint::default());
        assert_eq!(stages.lock().unwrap().last(), Some(&Stage::Done));
        let run = Run::default().threads(1);
        run.cancellation().cancel();
        let cancelled = request().generate(model, &run).unwrap();
        assert!(!cancelled.report().completed);
        assert!(cancelled.is_empty());
        let built = cancelled.build(&NumeratorOptions::default(), &run).unwrap();
        assert!(built.into_amplitude(Default::default()).is_err());
    }

    #[test]
    fn dressing_progress_precedes_preparation_and_can_cancel_it() {
        for parallel in [false, true] {
            for count_only in [false, true] {
                for cancel in [false, true] {
                    let updates = Arc::new(Mutex::new(Vec::new()));
                    let captured = updates.clone();
                    let run = Run::default()
                        .threads(2)
                        .parallel(parallel)
                        .progress(move |p| {
                            let stop = cancel && p.stage == Stage::Dressing;
                            captured.lock().unwrap().push(p);
                            if stop {
                                GenerationControl::Cancel
                            } else {
                                GenerationControl::Continue
                            }
                        });
                    let (diagrams, completed) = if count_only {
                        let result = request().count(model(), &run).unwrap();
                        (result.diagrams, result.complete)
                    } else {
                        let result = request().generate(model(), &run).unwrap();
                        (BigUint::from(result.len()), result.report.completed)
                    };
                    assert_eq!(completed, !cancel);
                    assert_eq!(diagrams, BigUint::from(if cancel { 0usize } else { 3 }));
                    let updates = updates.lock().unwrap();
                    let visits = updates
                        .iter()
                        .filter(|p| matches!(p.stage, Stage::Dressing | Stage::Prepare))
                        .map(|p| (p.stage, p.completed.clone()))
                        .collect::<Vec<_>>();
                    if cancel {
                        assert_eq!(visits, [(Stage::Dressing, BigUint::from(1usize))]);
                    } else {
                        // Two workers prepare the dressed assignments in batches: each stage
                        // counts in order, and an assignment is dressed before it is prepared.
                        for stage in [Stage::Dressing, Stage::Prepare] {
                            assert_eq!(
                                visits
                                    .iter()
                                    .filter(|(s, _)| *s == stage)
                                    .map(|(_, n)| n.clone())
                                    .collect::<Vec<_>>(),
                                (1usize..=3).map(BigUint::from).collect::<Vec<_>>()
                            );
                        }
                        let position = |visit: (Stage, usize)| {
                            visits
                                .iter()
                                .position(|v| *v == (visit.0, BigUint::from(visit.1)))
                                .unwrap()
                        };
                        for n in 1..=3 {
                            assert!(position((Stage::Dressing, n)) < position((Stage::Prepare, n)));
                        }
                    }
                    assert_eq!(updates.last().unwrap().stage, Stage::Done);
                    assert_eq!(updates.last().unwrap().completed, diagrams);
                }
            }
        }
    }

    #[test]
    fn preparation_progress_cancels_the_rest_of_its_batch() {
        for threads in [1, 8] {
            for count_only in [false, true] {
                let prepared = Arc::new(Mutex::new(Vec::new()));
                let captured = prepared.clone();
                let run = Run::default().threads(threads).progress(move |p| {
                    if p.stage != Stage::Prepare {
                        return GenerationControl::Continue;
                    }
                    captured.lock().unwrap().push(p.completed);
                    GenerationControl::Cancel
                });
                // All three assignments share one batch with eight threads.
                let (diagrams, completed) = if count_only {
                    let result = request().count(model(), &run).unwrap();
                    (result.diagrams, result.complete)
                } else {
                    let result = request().generate(model(), &run).unwrap();
                    (BigUint::from(result.len()), result.report.completed)
                };
                assert_eq!((diagrams, completed), (BigUint::from(1usize), false));
                assert_eq!(*prepared.lock().unwrap(), [BigUint::from(1usize)]);
            }
        }
    }

    #[test]
    fn ordered_maps_report_the_first_failure_in_item_order() {
        let run = Run::default().threads(8);
        for _ in 0..20 {
            let failure = run
                .in_pool(|| {
                    run.map_ordered((0..1000).collect(), |item: usize| {
                        if item % 100 == 37 {
                            Err(GenerationError::InvalidRequest(item.to_string()))
                        } else {
                            Ok(item)
                        }
                    })
                })
                .unwrap()
                .into_iter()
                .collect::<Result<Vec<_>, _>>();
            assert!(matches!(failure, Err(GenerationError::InvalidRequest(item)) if item == "37"));
        }
    }

    #[test]
    fn selection_reports_scanned_inputs_and_preserves_cancelled_prefixes() {
        let normal = Run::default().threads(1);
        let built = request()
            .generate(model(), &normal)
            .unwrap()
            .build(&NumeratorOptions::default(), &normal)
            .unwrap();
        assert_eq!(built.len(), 3);
        for stop in [None, Some(0usize), Some(1)] {
            let updates = Arc::new(Mutex::new(Vec::new()));
            let captured = updates.clone();
            let run = Run::default().threads(2).progress(move |p| {
                let cancel = p.stage == Stage::Selection
                    && stop.is_some_and(|n| p.completed == BigUint::from(n));
                captured.lock().unwrap().push(p);
                if cancel {
                    GenerationControl::Cancel
                } else {
                    GenerationControl::Continue
                }
            });
            let selected = built.select(&Selection::default(), &run).unwrap();
            assert_eq!(selected.len(), stop.unwrap_or(3));
            assert_eq!(selected.report.completed, stop.is_none());
            assert_eq!(selected.sources, built.sources[..selected.len()]);
            assert_eq!(selected.classes, built.classes[..selected.len()]);
            selected.validate_groups().unwrap();
            let updates = updates.lock().unwrap();
            assert_eq!(updates.first().unwrap().stage, Stage::Selection);
            assert_eq!(updates.first().unwrap().completed, BigUint::default());
            assert_eq!(updates.first().unwrap().total, Some(BigUint::from(3usize)));
            assert_eq!(updates.last().unwrap().stage, Stage::Done);
            assert_eq!(
                updates.last().unwrap().completed,
                BigUint::from(selected.len())
            );
            if stop.is_some() {
                let again = selected.select(&Selection::default(), &normal).unwrap();
                assert!(!again.report.completed);
                assert!(matches!(
                    again.into_amplitude(Default::default()),
                    Err(GenerationError::IncompleteAmplitude)
                ));
            }
        }
        let updates = Arc::new(Mutex::new(Vec::new()));
        let captured = updates.clone();
        let run = Run::default().threads(1).progress(move |p| {
            captured.lock().unwrap().push(p);
            GenerationControl::Continue
        });
        let empty = built
            .select(&Selection::default().select_ids([]), &run)
            .unwrap();
        assert!(empty.is_empty() && empty.report.completed);
        assert_eq!(
            updates
                .lock()
                .unwrap()
                .iter()
                .filter(|p| p.stage == Stage::Selection)
                .map(|p| p.completed.clone())
                .collect::<Vec<_>>(),
            (0usize..=3).map(BigUint::from).collect::<Vec<_>>()
        );
        updates.lock().unwrap().clear();
        assert!(
            empty
                .select(&Selection::default(), &run)
                .unwrap()
                .report
                .completed
        );
        assert_eq!(updates.lock().unwrap().len(), 2);
        assert!(matches!(
            built.select(&Selection::default(), &Run::default().threads(0)),
            Err(GenerationError::InvalidThreadCount(0))
        ));
        let run = Run::default().threads(1).progress(|p| {
            if p.stage == Stage::Done {
                GenerationControl::Cancel
            } else {
                GenerationControl::Continue
            }
        });
        assert!(
            !built
                .select(&Selection::default(), &run)
                .unwrap()
                .report
                .completed
        );
        let forced = Selection::default().with_forced_cuts([[EdgeId(usize::MAX)]]);
        assert!(matches!(
            built.select(&forced, &normal),
            Err(GenerationError::UnknownForcedCut { .. })
        ));
        let run = Run::default()
            .threads(1)
            .progress(|_| GenerationControl::Cancel);
        let cancelled = built.select(&forced, &run).unwrap();
        assert!(cancelled.is_empty() && !cancelled.report.completed);
        // An incomplete input cannot prove that a requested cut is globally absent.
        assert!(!cancelled.select(&forced, &normal).unwrap().report.completed);
    }

    #[test]
    fn explain_separates_native_search_from_physical_and_partner_checks() {
        let model = model();
        let amplitude = request().with_orders(Orders::new().squared_exact("QCD", 4));
        let report = amplitude.explain(&model).unwrap();
        let partner = report
            .criteria
            .iter()
            .find(|c| c.criterion == "squared-order partner relation")
            .unwrap();
        assert_eq!(partner.stage, PlanStage::Preparation);
        assert_eq!(partner.role, PlanRole::FinalCheck);
        assert!(!partner.hereditary);
        assert!(partner.effective.contains("one common observed partner"));
        let cross = Process::new(["scalar_0", "scalar_0"], ["scalar_0", "scalar_0"])
            .cross_section(Loops::forward(1..=1))
            // This fixture includes two-point vertices, so its search needs a vertex bound.
            .with_topology(TopologyRules::default().with_max_vertices(4))
            .with_cuts(Cuts::default().with_cut_amplitude(CutAmplitude {
                loops: Range::exact(0),
                orders: Orders::new().exact("QCD", 2),
            }));
        let report = cross.explain(&model).unwrap();
        let cut = report
            .criteria
            .iter()
            .find(|c| c.criterion == "cut amplitude")
            .unwrap();
        assert_eq!(cut.stage, PlanStage::Preparation);
        assert_eq!(cut.role, PlanRole::FinalCheck);
        assert!(!cut.hereditary);
        assert_eq!(
            cut.effective,
            "sum(left loops, right loops)=Inclusive { min: 0, max: Some(0) }; each side's coupling orders={\"QCD\": Inclusive { min: 2, max: Some(2) }}"
        );
        assert!(
            !report
                .criteria
                .iter()
                .any(|c| c.criterion == "squared-order partner relation")
        );
    }

    #[test]
    fn squared_constraints_require_one_common_partner() {
        let request = request().with_orders(
            Orders::new()
                .squared_exact("QCD", 2)
                .squared_exact("QED", 2),
        );
        let a = BTreeMap::from([("QCD".into(), 1), ("QED".into(), 1)]);
        let b = BTreeMap::from([("QCD".into(), 1), ("QED".into(), 3)]);
        let c = BTreeMap::from([("QCD".into(), 3), ("QED".into(), 1)]);
        // a would incorrectly pass by taking QCD from b and QED from c; remove its self-pair.
        let isolated = BTreeMap::from([("QCD".into(), 0), ("QED".into(), 0)]);
        let classes = OrderClasses::resolve(
            &request,
            [
                isolated.clone(),
                BTreeMap::from([("QCD".into(), 2), ("QED".into(), 4)]),
                BTreeMap::from([("QCD".into(), 4), ("QED".into(), 2)]),
            ],
        )
        .unwrap();
        assert!(!classes.retained(classes.class(&isolated)));
        let classes = OrderClasses::resolve(&request, [a.clone(), b, c]).unwrap();
        assert_eq!(
            classes.allowed,
            BTreeSet::from([(classes.class(&a), classes.class(&a))])
        );
        let pairing = classes.pairing([]).unwrap();
        assert!(pairing.term_classes().is_empty());
        assert_eq!(pairing.allowed_pairs(), &classes.allowed);
    }

    #[test]
    fn selection_and_grouping_retain_squared_relation_and_source_terms() {
        let model = model();
        let run = Run::default().threads(1);
        let plain = request().generate(model.clone(), &run).unwrap();
        let orders = plain.get(0).unwrap().orders();
        let (name, power) = orders.iter().next().unwrap();
        let constrained = request().with_orders(Orders::new().squared_exact(name, 2 * power));
        let built = constrained
            .generate(model, &run)
            .unwrap()
            .build(
                &NumeratorOptions::default().with_naming(Naming::named()),
                &run,
            )
            .unwrap();
        let selected = built
            .select(
                &Selection::default().select_names([built.get(1).unwrap().name().unwrap()]),
                &run,
            )
            .unwrap();
        assert_eq!(selected.len(), 1);
        let grouping = Grouping::Identical(GraphGroupingOptions {
            check_canonical_numerator: true,
            ..Default::default()
        });
        let grouped = selected.group(&grouping, &run).unwrap();
        grouped.validate_groups().unwrap();
        assert_eq!(grouped.groups[0].members[0].source_diagram, 1);
        let amplitude = grouped.into_amplitude(Default::default()).unwrap();
        assert!(amplitude.pairing().is_some());
        assert_eq!(amplitude.pairing().unwrap().term_classes().len(), 1);
        let empty = built
            .select(
                &Selection::default().select_ids(std::iter::empty::<DiagramId>()),
                &run,
            )
            .unwrap();
        assert!(empty.is_empty());
        assert!(empty.into_amplitude(Default::default()).is_err());
    }

    #[test]
    fn distinct_order_classes_keep_only_oriented_interference_after_selection_and_grouping() {
        let mut definition: Value = serde_json::from_str(include_str!(
            "../../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
        ))
        .unwrap();
        let low = definition["vertex_rules"]
            .as_array()
            .unwrap()
            .iter()
            .find(|rule| rule["name"] == "V_3_SCALAR_000")
            .unwrap()
            .clone();
        let mut high = low.clone();
        high["name"] = "V_3_HIGH".into();
        high["couplings"] = json!([["HIGH"]]);
        definition["vertex_rules"] = json!([low, high]);
        let mut coupling = definition["couplings"][0].clone();
        coupling["name"] = "HIGH".into();
        coupling["orders"] = json!([["QCD", 2]]);
        definition["couplings"]
            .as_array_mut()
            .unwrap()
            .push(coupling);
        let model = Arc::new(Model::from_json(&definition.to_string()).unwrap());
        let run = Run::default().threads(1);
        let request = Process::new(["scalar_0"], ["scalar_0", "scalar_0"])
            .amplitude(Loops::total(0..=0))
            .with_orders(Orders::new().squared_exact("QCD", 3));
        let prepared = request.generate(model, &run).unwrap();
        assert_eq!(prepared.len(), 2);
        let built = prepared
            .build(
                &NumeratorOptions::default().with_naming(Naming::named()),
                &run,
            )
            .unwrap();
        let stop = Run::default().threads(1).progress(|p| {
            if p.stage == Stage::Selection && p.completed == BigUint::from(1usize) {
                GenerationControl::Cancel
            } else {
                GenerationControl::Continue
            }
        });
        let prefix = built.select(&Selection::default(), &stop).unwrap();
        assert_eq!(prefix.len(), 1);
        assert!(!prefix.report.completed);
        assert!(Arc::ptr_eq(
            prefix.prepared.as_ref().unwrap(),
            built.prepared.as_ref().unwrap()
        ));
        assert_eq!(prefix.classes, built.classes[..1]);
        let relation = prefix.relation.as_ref().unwrap();
        assert_eq!(relation.allowed.len(), 2);
        assert_eq!(
            relation
                .pairing(prefix.classes.iter().copied())
                .unwrap()
                .term_classes()
                .len(),
            1
        );
        prefix.validate_groups().unwrap();
        let grouping = Grouping::Identical(GraphGroupingOptions {
            check_canonical_numerator: true,
            ..Default::default()
        });
        let grouped = built.clone().group(&grouping, &run).unwrap();
        assert_eq!(
            grouped.len(),
            2,
            "equal numerators of different classes must not collapse"
        );
        assert!(
            grouped.prepared.is_none(),
            "grouping releases the preparations"
        );
        let amplitude = grouped.into_amplitude(Default::default()).unwrap();
        let pairing = amplitude.pairing().unwrap();
        assert_eq!(pairing.allowed_pairs(), &BTreeSet::from([(0, 1), (1, 0)]));
        let unrestricted =
            Amplitude::new(amplitude.diagrams().iter().cloned(), Default::default()).unwrap();
        assert_eq!(
            amplitude.squared().unwrap().expression().expression() * 2,
            unrestricted
                .squared()
                .unwrap()
                .expression()
                .expression()
                .clone()
        );
        let selected = built
            .select(
                &Selection::default().select_names([built.get(0).unwrap().name().unwrap()]),
                &run,
            )
            .unwrap();
        let selected = selected
            .group(&grouping, &run)
            .unwrap()
            .into_amplitude(Default::default())
            .unwrap();
        assert_eq!(
            selected.pairing().unwrap().allowed_pairs(),
            pairing.allowed_pairs()
        );
        assert!(
            selected.squared().unwrap().expression().is_zero(),
            "the missing partner remains zero, not unrestricted"
        );
    }

    #[test]
    fn signed_frame_transport_composes_without_rewriting_twice() {
        let particle = ParticleId::from_index(0);
        let edge = |source, sink| GeneratedEdge {
            source: Endpoint::Vertex(source),
            sink: Endpoint::Vertex(sink),
            particle,
            directed: true,
        };
        let raw = GeneratedGraph::new(
            [0, 1].map(VertexRuleId::from_index),
            [edge(0, 1), edge(1, 0)],
        );
        let reframe = |graph: &GeneratedGraph, edges: [(usize, bool); 2]| {
            graph.reframed(
                &[1, 0],
                edges.map(|(edge, reversed)| (edge, reversed, particle)),
                |leg| leg,
            )
        };
        let first_graph = reframe(&raw, [(1, false), (0, true)]);
        let second_graph = reframe(&first_graph, [(1, true), (0, true)]);
        let first = raw.frame_map(&first_graph);
        assert_eq!(first.vertices, vec![1, 0]);
        assert_eq!(first.edges, vec![(1, true), (0, false)]);
        let second = first_graph.frame_map(&second_graph);
        assert_eq!(second.vertices, vec![1, 0]);
        assert_eq!(second.edges, vec![(1, true), (0, true)]);
        // Each frame is kept relative to the base assignment, so maps compose without
        // rewriting an intermediate frame.
        let composed = raw.frame_map(&second_graph);
        assert_eq!(composed.vertices, vec![0, 1]);
        assert_eq!(composed.edges, vec![(0, false), (1, true)]);
        let q = graph_api::momentum_symbol().call(0);
        let h = symbols::hedge_index().call(0);
        let v = symbols::vertex_index().call(0);
        let options = NumeratorOptions::default()
            .with_prefactor(q.clone() * h.clone())
            .with_projector(q * v);
        let moved = first.options(&options).unwrap();
        assert_eq!(
            moved.prefactor,
            -graph_api::momentum_symbol().call(1) * symbols::hedge_index().call(3)
        );
        assert_eq!(
            moved.projector.unwrap(),
            -graph_api::momentum_symbol().call(1) * symbols::vertex_index().call(1)
        );
        assert_eq!(
            composed.options(&options).unwrap(),
            second.options(&first.options(&options).unwrap()).unwrap()
        );
    }
}
