//! Owned physical preparation, explicit numerator construction, selection and grouping.

use super::*;
use std::collections::HashSet;

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

impl fmt::Debug for PreparedDiagram {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
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
    pub fn graph(&self) -> &HedgeGraph<ParticleId, VertexRuleId, Option<ExternalNode>> {
        self.input.representative.graph.graph()
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
    ) -> Result<feynkit_amplitude::SquaredPairing, GenerationError> {
        feynkit_amplitude::SquaredPairing::new(
            self.orders.len(),
            classes,
            self.allowed.iter().copied(),
        )
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
impl fmt::Debug for DiagramSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
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
    pub fn iter(&self) -> std::slice::Iter<'_, PreparedDiagram> {
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
        let width = self.len().saturating_sub(1).to_string().len();
        let mut diagrams = Vec::with_capacity(self.len());
        let mut sources = Vec::with_capacity(self.len());
        run.report_progress(Stage::Numerators, 0usize, Some(self.len()));
        let build = || -> Result<_, GenerationError> {
            for (index, input) in self.diagrams.iter().enumerate() {
                if run.cancellation_requested() {
                    break;
                }
                let name = match &options.naming {
                    Naming::None => None,
                    Naming::Named { prefix } => Some(format!("{prefix}{index:0width$}")),
                };
                let pair = self.context.generator.diagram_pair(
                    name,
                    self.context.search.request.kind,
                    (*input.input).clone(),
                    (options, options),
                )?;
                let mut diagram = pair.representative;
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
                diagrams.push(diagram);
                sources.push(index);
                run.report_progress(Stage::Numerators, index + 1, Some(self.len()));
            }
            Ok(())
        };
        run.in_pool(build)??;
        let mut report = self.report.clone();
        report.completed &= !run.cancellation_requested();
        let mut result = BuiltDiagramSet {
            diagrams,
            groups: Vec::new(),
            report,
            prepared: Arc::new(self.clone()),
            options: options.clone(),
            sources,
            classes: Vec::new(),
            grouped: false,
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
                    result.diagrams.clone(),
                    &self.context.generator.model,
                    &Grouping::None,
                    false,
                    run,
                    true,
                    None,
                    None,
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
        run.report_progress_with_zero_count(
            Stage::Done,
            result.len(),
            Some(result.len()),
            result.report.zero_numerator_count.clone(),
        );
        result.report.completed &= !run.cancellation_requested();
        Ok(result)
    }
}

/// Built physical diagrams and their source provenance. Selection and grouping keep the
/// resolved squared-order relation; raw lists passed to `Amplitude::new` remain unrestricted.
#[derive(Clone)]
pub struct BuiltDiagramSet {
    diagrams: Vec<FeynmanDiagram>,
    groups: Vec<DiagramGroup>,
    report: GenerationReport,
    prepared: Arc<DiagramSet>,
    options: NumeratorOptions,
    sources: Vec<usize>,
    classes: Vec<usize>,
    grouped: bool,
}
impl fmt::Debug for BuiltDiagramSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
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
    pub fn iter(&self) -> std::slice::Iter<'_, FeynmanDiagram> {
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
    pub fn into_amplitude(
        self,
        options: feynkit_amplitude::AmplitudeOptions,
    ) -> Result<feynkit_amplitude::Amplitude, GenerationError> {
        if !self.report.completed {
            return Err(GenerationError::IncompleteAmplitude);
        }
        let pairing = self
            .prepared
            .classes
            .as_ref()
            .map(|classes| classes.pairing(self.classes.iter().copied()))
            .transpose()?;
        let amplitude =
            feynkit_amplitude::Amplitude::new(self.diagrams.into_iter().map(Arc::new), options)
                .map_err(Box::new)?;
        Ok(match pairing {
            Some(pairing) => amplitude.with_pairing(pairing).map_err(Box::new)?,
            None => amplitude,
        })
    }
    fn singleton_groups(&mut self) -> Result<(), GenerationError> {
        self.groups = self
            .diagrams
            .iter()
            .enumerate()
            .map(|(i, diagram)| {
                Ok(DiagramGroup {
                    master: i,
                    members: vec![GroupMember {
                        source_diagram: self.sources[i],
                        source_identity: diagram
                            .name()
                            .map(|name| {
                                Ok::<_, GenerationError>(SourceIdentity {
                                    id: diagram.id()?,
                                    name: name.to_owned(),
                                })
                            })
                            .transpose()?,
                        diagram: i,
                        ratio: Atom::one(),
                        overall_factor: diagram.overall_factor().clone(),
                    }],
                })
            })
            .collect::<Result<_, GenerationError>>()?;
        self.validate_groups()
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
        if self.grouped && !selection.forced_cuts.is_empty() {
            return Err(GenerationError::InvalidRequest(
                "forced cuts must be selected before grouping".into(),
            ));
        }
        let mut result = self.clone();
        result.diagrams.clear();
        result.sources.clear();
        result.classes.clear();
        result.groups.clear();
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
        run.report_progress_with_zero_count(
            Stage::Done,
            result.len(),
            Some(result.len()),
            result.report.zero_numerator_count.clone(),
        );
        result.report.completed &= !run.cancellation_requested();
        Ok(result)
    }
}

impl Run {
    fn in_pool<T: Send>(&self, work: impl FnOnce() -> T + Send) -> Result<T, GenerationError> {
        if self.threads == Some(0) {
            return Err(GenerationError::InvalidThreadCount(0));
        }
        #[cfg(target_arch = "wasm32")]
        {
            Ok(work())
        }
        #[cfg(not(target_arch = "wasm32"))]
        {
            let mut builder = ThreadPoolBuilder::new();
            if let Some(threads) = self.threads {
                builder = builder.num_threads(threads)
            }
            Ok(builder.build()?.install(work))
        }
    }
}

impl PreparedSearch {
    /// Visit one assignment at a time. Count uses this same physical preparation boundary and
    /// retains only exact orbit keys/class accumulators, never all prepared graphs or numerators.
    fn visit(
        &self,
        generator: &Generator,
        run: &Run,
        mut accept: impl FnMut(PreparedDiagram) -> Result<(), GenerationError> + Send,
    ) -> Result<(GenerationReport, BTreeMap<usize, LoopCounts>), GenerationError> {
        let mut prepared = BigUint::default();
        self.enumerate(run, |compact| {
            if let Some(multiplicity) = self.externals.sewn_pair_orbit(&compact) {
                let aut = compact.symmetry_factor()?;
                let graph = self.externals.graph(&compact, &self.signatures);
                if let Some(input) = generator.diagram_input(
                    &self.resolved,
                    graph,
                    aut,
                    multiplicity,
                    &self.request,
                    &Grouping::None,
                )? {
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
                    accept(PreparedDiagram {
                        input: Arc::new(input),
                        compact: Arc::new(compact),
                        skeleton,
                        orders,
                    })?;
                }
            }
            prepared += 1u8;
            run.report_progress(Stage::Prepare, prepared.clone(), None);
            Ok(())
        })
    }

    pub(super) fn enumerate(
        &self,
        run: &Run,
        mut accept: impl FnMut(CompactDiagram<ParticleId, ExternalLabel>) -> Result<(), GenerationError>
        + Send,
    ) -> Result<(GenerationReport, BTreeMap<usize, LoopCounts>), GenerationError> {
        let progress = Mutex::new((BigUint::default(), BTreeMap::<usize, LoopCounts>::new()));
        let engine = EngineRun::default()
            .parallel(run.parallel)
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
        let mut visit =
            |compact: CompactDiagram<ParticleId, ExternalLabel>| -> Result<(), GenerationError> {
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
                accept(compact)
            };
        let enumerate = || -> Result<(), GenerationError> {
            if self.externals.vetoed(&self.restrictions) {
                return Ok(());
            }
            let Some(plan) = &self.plan else {
                return Ok(());
            };
            if let Some((_, dispatch)) = &self.request.deferred_hooks {
                let mut driver = plan.batched(&engine, batch::Output::Diagrams)?;
                loop {
                    match driver.poll()? {
                        batch::Event::Batch(requests) => match dispatch(requests) {
                            ControlFlow::Continue(replies) => driver.resolve(&replies)?,
                            ControlFlow::Break(()) => {
                                run.cancellation.cancel();
                                driver.abort();
                            }
                        },
                        batch::Event::Diagram(diagram) => visit(diagram)?,
                        batch::Event::Yield | batch::Event::Progress(_) => {
                            if run.cancellation_requested() {
                                driver.abort()
                            }
                        }
                        batch::Event::Done => break,
                        batch::Event::Skeleton(_) => unreachable!("diagram output requested"),
                    }
                }
            } else {
                for diagram in plan.diagrams(&engine)? {
                    visit(diagram)?;
                }
            }
            Ok(())
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
    pub fn explain(
        &self,
        model: &Model,
    ) -> Result<linnet::half_edge::generate::PlanReport, GenerationError> {
        use linnet::half_edge::generate::{PlanCriterion, PlanRole, PlanStage};

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
        if self.restrictions != crate::Restrict::default() {
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
            || topology.factorized_loops != crate::Range::Unrestricted
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
        if self.symmetrization != crate::Symmetrize::default() {
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
            if self.cuts.amplitude != crate::CutAmplitude::default() {
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
    pub fn group(&self, grouping: &Grouping, run: &Run) -> Result<Self, GenerationError> {
        if matches!(grouping, Grouping::None) {
            let mut result = self.clone();
            run.report_progress_with_zero_count(
                Stage::Done,
                result.len(),
                Some(result.len()),
                result.report.zero_numerator_count.clone(),
            );
            result.report.completed &= !run.cancellation_requested();
            return Ok(result);
        }
        if self.grouped {
            return Err(GenerationError::InvalidRequest("group an ungrouped built set; grouping a collapsed set would lose its source terms".into()));
        }
        if matches!(grouping,Grouping::Identical(options)|Grouping::UpToSign(options)|Grouping::UpToScalar(options) if options.number_of_numerical_samples==0)
        {
            return Err(GenerationError::InvalidNumericalSampleCount);
        }
        let context = &self.prepared.context;
        let mut order = (0..self.len()).collect::<Vec<_>>();
        // This is the former pre-preparation ordering boundary, demanded only by grouping.
        order.sort_by(|&a, &b| {
            self.prepared.diagrams[self.sources[a]]
                .compact
                .canonical_form()
                .cmp(
                    self.prepared.diagrams[self.sources[b]]
                        .compact
                        .canonical_form(),
                )
        });
        let prepare = || -> Result<Vec<DiagramPair>, GenerationError> {
            let mut pairs = Vec::with_capacity(order.len());
            for &index in &order {
                if run.cancellation_requested() {
                    break;
                }
                let source = &self.prepared.diagrams[self.sources[index]];
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
                let source_edges = input
                    .representative_frame
                    .edges
                    .iter()
                    .map(|&(edge, _)| edge)
                    .collect::<Vec<_>>();
                let comparison_options = source
                    .input
                    .representative_frame
                    .relative_to(&input.source_frame)
                    .options(&self.options)?;
                let representative_options = source
                    .input
                    .representative_frame
                    .relative_to(&input.representative_frame)
                    .options(&self.options)?;
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
                let old_internal = source
                    .input
                    .representative
                    .graph
                    .edges()
                    .enumerate()
                    .filter_map(|(edge, data)| (!data.is_external()).then_some(edge));
                let old_final = self.diagrams[index]
                    .edges()
                    .filter_map(|(edge, _, data)| data.external.is_none().then_some(edge));
                let new_final = pair
                    .representative
                    .edges()
                    .filter_map(|(edge, _, data)| data.external.is_none().then_some(edge))
                    .collect::<Vec<_>>();
                let new_ids = new_internal
                    .into_iter()
                    .zip(new_final)
                    .collect::<BTreeMap<_, _>>();
                let inverse = source
                    .input
                    .representative_frame
                    .edges
                    .iter()
                    .enumerate()
                    .map(|(old, &(new, _))| (new, old))
                    .collect::<BTreeMap<_, _>>();
                let map = old_internal
                    .zip(old_final)
                    .map(|(edge, final_edge)| (final_edge, new_ids[&source_edges[inverse[&edge]]]))
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
                    .representative
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
                if cuts != pair.representative.cuts() {
                    pair.representative = pair.representative.with_cuts(cuts)?;
                    let cuts =
                        transport_cut_metadata(pair.representative.cuts(), &pair.reversed_edges);
                    pair.comparison = pair.comparison.with_cut_partitions(
                        cuts.into_iter()
                            .map(|cut| (cut.left.half_edges, cut.right.half_edges))
                            .collect(),
                    )?;
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
                    pair.comparison = pair.comparison.with_loop_momentum_edges(&basis)?;
                    pair.representative = pair.representative.with_loop_momentum_edges(&basis)?;
                }
                pairs.push(pair);
            }
            Ok(pairs)
        };
        let pairs = run.in_pool(prepare)??;
        if run.cancellation_requested() {
            let mut result = self.clone();
            result.report.completed = false;
            run.report_progress_with_zero_count(
                Stage::Done,
                result.len(),
                Some(result.len()),
                result.report.zero_numerator_count.clone(),
            );
            return Ok(result);
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
        let grouped = run.in_pool(|| {
            grouping::group_diagrams(
                pairs.iter().map(|pair| pair.comparison.clone()).collect(),
                &context.generator.model,
                grouping,
                context.search.resolved.identifies_mirror_images(),
                run,
                false,
                self.prepared.classes.as_ref().map(|_| classes.as_slice()),
                Some(&source_indices),
            )
        })??;
        let mut result = self.clone();
        result.diagrams.clear();
        result.sources.clear();
        result.classes.clear();
        result.groups = grouped.groups;
        result.grouped = true;
        result.report.zero_numerator_count += BigUint::from(grouped.zero_numerator_count);
        for (output, (comparison, group)) in grouped
            .diagrams
            .into_iter()
            .zip(&mut result.groups)
            .enumerate()
        {
            let first = group
                .members
                .first()
                .ok_or_else(|| {
                    GenerationError::InvalidGroups("group has no representative".into())
                })?
                .source_diagram;
            let first = source_lookup[&first];
            let original = order[first];
            let source = &pairs[first];
            let mut representative = source.representative.clone();
            let cuts = transport_cut_metadata(comparison.cuts(), &source.reversed_edges);
            if cuts != representative.cuts() {
                representative = representative.with_cut_partitions(
                    cuts.into_iter()
                        .map(|cut| (cut.left.half_edges, cut.right.half_edges))
                        .collect(),
                )?;
            }
            let thresholds = transport_threshold_candidate_metadata(
                comparison.topology_threshold_candidates(),
                &source.reversed_edges,
            );
            if thresholds != representative.topology_threshold_candidates() {
                representative = representative.with_topology_threshold_partitions(
                    thresholds
                        .into_iter()
                        .map(|cut| (cut.left, cut.right))
                        .collect(),
                )?;
            }
            if representative.loop_momentum_basis().loop_edges
                != comparison.loop_momentum_basis().loop_edges
            {
                representative = representative
                    .with_loop_momentum_edges(&comparison.loop_momentum_basis().loop_edges)?;
            }
            representative =
                representative.with_overall_factor(comparison.overall_factor().clone());
            result.diagrams.push(representative);
            result.sources.push(self.sources[original]);
            result.classes.push(self.classes[original]);
            group.master = output;
            for member in &mut group.members {
                let original = order[source_lookup[&member.source_diagram]];
                let source = &self.diagrams[original];
                member.source_identity = source
                    .name()
                    .map(|name| {
                        Ok::<_, GenerationError>(SourceIdentity {
                            id: source.id()?,
                            name: name.to_owned(),
                        })
                    })
                    .transpose()?;
                member.source_diagram = self.sources[original];
                member.diagram = output;
            }
        }
        result.report.retained_count = result.len().into();
        run.report_progress_with_zero_count(
            Stage::Done,
            result.len(),
            Some(result.len()),
            result.report.zero_numerator_count.clone(),
        );
        result.report.completed &= !run.cancellation_requested();
        result.validate_groups()?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Loops, Orders, Restrict, TopologyRules};

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
        let cache = &set.context.generator.numerators;
        assert!(cache.interactions.read().unwrap().is_empty());
        assert!(cache.propagators.read().unwrap().is_empty());
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
            crate::GenerationControl::Continue
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
                                crate::GenerationControl::Cancel
                            } else {
                                crate::GenerationControl::Continue
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
                    let expected = if cancel {
                        vec![(Stage::Dressing, BigUint::from(1usize))]
                    } else {
                        (1usize..=3)
                            .flat_map(|n| {
                                [
                                    (Stage::Dressing, BigUint::from(n)),
                                    (Stage::Prepare, BigUint::from(n)),
                                ]
                            })
                            .collect()
                    };
                    assert_eq!(visits, expected);
                    assert_eq!(updates.last().unwrap().stage, Stage::Done);
                    assert_eq!(updates.last().unwrap().completed, diagrams);
                }
            }
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
                    crate::GenerationControl::Cancel
                } else {
                    crate::GenerationControl::Continue
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
            crate::GenerationControl::Continue
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
                crate::GenerationControl::Cancel
            } else {
                crate::GenerationControl::Continue
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
            .progress(|_| crate::GenerationControl::Cancel);
        let cancelled = built.select(&forced, &run).unwrap();
        assert!(cancelled.is_empty() && !cancelled.report.completed);
        // An incomplete input cannot prove that a requested cut is globally absent.
        assert!(!cancelled.select(&forced, &normal).unwrap().report.completed);
    }

    #[test]
    fn explain_separates_native_search_from_physical_and_partner_checks() {
        use linnet::half_edge::generate::{PlanRole, PlanStage};
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
            .with_cuts(
                crate::Cuts::default().with_cut_amplitude(crate::CutAmplitude {
                    loops: crate::Range::exact(0),
                    orders: Orders::new().exact("QCD", 2),
                }),
            );
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
        let grouping = Grouping::Identical(crate::GraphGroupingOptions {
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
        let mut definition: serde_json::Value = serde_json::from_str(include_str!(
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
        high["couplings"] = serde_json::json!([["HIGH"]]);
        definition["vertex_rules"] = serde_json::json!([low, high]);
        let mut coupling = definition["couplings"][0].clone();
        coupling["name"] = "HIGH".into();
        coupling["orders"] = serde_json::json!([["QCD", 2]]);
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
                crate::GenerationControl::Cancel
            } else {
                crate::GenerationControl::Continue
            }
        });
        let prefix = built.select(&Selection::default(), &stop).unwrap();
        assert_eq!(prefix.len(), 1);
        assert!(!prefix.report.completed);
        assert!(Arc::ptr_eq(&prefix.prepared, &built.prepared));
        assert_eq!(prefix.classes, built.classes[..1]);
        let relation = prefix.prepared.classes.as_ref().unwrap();
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
        let grouping = Grouping::Identical(crate::GraphGroupingOptions {
            check_canonical_numerator: true,
            ..Default::default()
        });
        let grouped = built.group(&grouping, &run).unwrap();
        assert_eq!(
            grouped.len(),
            2,
            "equal numerators of different classes must not collapse"
        );
        let amplitude = grouped.into_amplitude(Default::default()).unwrap();
        let pairing = amplitude.pairing().unwrap();
        assert_eq!(pairing.allowed_pairs(), &BTreeSet::from([(0, 1), (1, 0)]));
        let unrestricted = feynkit_amplitude::Amplitude::new(
            amplitude.diagrams().iter().cloned(),
            Default::default(),
        )
        .unwrap();
        assert_eq!(
            amplitude.squared().unwrap().expression() * 2,
            unrestricted.squared().unwrap().expression().clone()
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
        assert_eq!(
            selected.squared().unwrap().expression(),
            &Atom::Zero,
            "the missing partner remains zero, not unrestricted"
        );
    }

    #[test]
    fn signed_frame_transport_composes_without_rewriting_twice() {
        let first = GeneratedFrame {
            vertices: vec![1, 0],
            edges: vec![(1, true), (0, false)],
        };
        let second = GeneratedFrame {
            vertices: vec![1, 0],
            edges: vec![(1, true), (0, true)],
        };
        let composed = first.then(&second);
        assert_eq!(composed.vertices, vec![0, 1]);
        assert_eq!(composed.edges, vec![(0, false), (1, true)]);
        let relative = first.relative_to(&composed);
        assert_eq!(relative.vertices, second.vertices);
        assert_eq!(relative.edges, second.edges);
        let q = momentum_symbol().call(0);
        let h = feynkit_graph::symbols::hedge_index().call(0);
        let v = feynkit_graph::symbols::vertex_index().call(0);
        let options = NumeratorOptions::default()
            .with_prefactor(q.clone() * h.clone())
            .with_projector(q * v);
        let moved = first.options(&options).unwrap();
        assert_eq!(
            moved.prefactor,
            -momentum_symbol().call(1) * feynkit_graph::symbols::hedge_index().call(3)
        );
        assert_eq!(
            moved.projector.unwrap(),
            -momentum_symbol().call(1) * feynkit_graph::symbols::vertex_index().call(1)
        );
        assert_eq!(
            first.then(&second).options(&options).unwrap(),
            second.options(&first.options(&options).unwrap()).unwrap()
        );
    }
}
