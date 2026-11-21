//! Bare propagator-insertion identities with the original UV forest kept separate.

use std::{
    collections::{BTreeSet, HashMap},
    sync::{Arc, LazyLock, Mutex, OnceLock},
    time::{Duration, Instant},
};

use eyre::{Result, WrapErr, ensure, eyre};
use feynkit_graph::{
    DiagramEdge, DiagramRelabeling, DiagramVertex, ExternalLeg, ExternalState, FeynmanDiagram,
    NumeratorInterface, PropagatorInsertion, VertexId,
};
use idenso::{CookMode, CookSettings, tensor::SymbolicTensor};
use spenso::{network::parsing::ParseState, structure::abstract_index::AbstractIndex};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function, symbol,
};
use vakint::Vakint;

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::FeynmanGraph,
    numerator::aind::Aind,
    utils::{GS, W_},
    uv::{
        UVgenerationSettings,
        counterterm_cache::{CacheTimingPhase, CacheTimingStats, CacheTimings},
        counterterm_identity::{CountertermIdentity, CountertermKey},
    },
};

use linnet::half_edge::{
    NodeIndex,
    subgraph::{SuBitGraph, SubSetOps},
};

use super::{
    CountertermCache, InsertionSkeletonIdentity, RStarPlan, expression::RStarIntegrand,
    integrate_propagator,
};

type CachedIntegral = Arc<OnceLock<Result<Atom, String>>>;

/// Every slashed momentum `gamma(i, j, Q(e, R(D)))` of an insertion
/// coefficient (idenso's canonical form, e.g. a quark self-energy `c slash(p)`)
/// written as `gamma(i, j, R(D, nu_k)) Q(e, R(D, nu_k))` with a contracted symbol
/// index `nu_k`, unique in the process. Applied after
/// [`freshen_insertion_indices`], whose `dummy(N)` indices Spenso does not parse. The momentum then
/// reaches the surrogate's coefficient vertex only through `Q`, as in the
/// Feynman rules, so the UV-child degree truncation can grade it
/// (`polynomial_in_roots`); a function of the momentum declines it. Chains
/// are undone first, so no product lands inside a chain. Exact: the same tensor.
pub(super) fn open_slashes(coefficient: &Atom) -> Atom {
    let is_slash = |view: AtomView<'_>| {
        matches!(view, AtomView::Fun(gamma) if gamma.get_symbol() == idenso::dirac::AGS.gamma
            && gamma.get_nargs() == 3
            && matches!(gamma.get(2), AtomView::Fun(momentum) if momentum.get_symbol() == GS.emr_mom
                && momentum.get_nargs() == 2
                && matches!(momentum.get(1), AtomView::Fun(_))))
    };
    let mut found = false;
    coefficient.visitor(&mut |view| {
        found |= is_slash(view);
        !found
    });
    if !found {
        return coefficient.clone();
    }
    let unchained = if coefficient.contains_symbol(spenso::network::tags::SPENSO_TAG.chain) {
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        match cooking
            .try_cook(coefficient.as_view())
            .map_err(|error| eyre!("{error}"))
            .and_then(|cooked| Ok(SymbolicTensor::infer_indexed(cooked)?.undo_chain()?))
        {
            Ok(tensor) => cooking.uncook(tensor.expression().as_view()),
            // Not representable as a tensor: leave it, the truncation declines.
            Err(_) => return coefficient.clone(),
        }
    } else {
        coefficient.clone()
    };
    // Process-unique: a surrogate may carry several coefficients.
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    unchained.replace_map(|view, _, out| {
        if !is_slash(view) {
            return;
        }
        let AtomView::Fun(gamma) = view else { return };
        let AtomView::Fun(momentum) = gamma.get(2) else {
            return;
        };
        let AtomView::Fun(bare) = momentum.get(1) else {
            return;
        };
        let next = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let index = Atom::from(Aind::Symbol(symbol!(&format!(
            "gammalooprs::insertion_slash_{next}"
        ))));
        let mut slot = symbolica::atom::FunctionBuilder::new(bare.get_symbol());
        for argument in bare.iter() {
            slot = slot.add_arg(argument);
        }
        let slot = slot.add_arg(&index).finish();
        **out = function!(idenso::dirac::AGS.gamma, gamma.get(0), gamma.get(1), &slot)
            * function!(GS.emr_mom, momentum.get(0), &slot);
    })
}

/// Cached open tensors share integration-local names. Freshen only contracted
/// indices on each use, preserving their established public boundary ports.
fn freshen_insertion_indices(coefficient: &Atom, reservation_context: &[Atom]) -> Result<Atom> {
    let state = ParseState::<Aind>::default();
    for expression in reservation_context {
        state.reserve_indices(expression.as_view());
    }
    state.reserve_indices(coefficient.as_view());
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    static MARKER: LazyLock<Symbol> = LazyLock::new(|| symbol!("rstar_insertion_internal_scope"));
    let marker = *MARKER;
    let scoped = SymbolicTensor::infer(cooking.try_cook(coefficient.as_view())?)?
        .wrap_indices(marker, true)?;
    let mut names = HashMap::new();
    let expression = scoped.expression().replace_map(|term, _, out| {
        if let AtomView::Fun(fun) = term
            && fun.get_symbol() == AbstractIndex::scope_symbol()
            && matches!(fun.iter().next(), Some(AtomView::Var(scope)) if scope.get_symbol() == marker)
        {
            // Lower the temporary native scope to the existing runtime Aind
            // representation, which the R* tensor/parser bridge accepts.
            **out = names
                .entry(term.to_owned())
                .or_insert_with(|| Atom::from(state.fresh_index()))
                .clone();
        }
    });
    Ok(cooking.uncook(expression.as_view()))
}

/// `GAMMALOOP_SURROGATE_CLOSED_ROOT=1` (default on: production since r5r-nf1,
/// results/surrogate-closed-root): above four loops, a caller that closes the
/// root metric ([`PropagatorInsertionEvaluator::pole_part_closing_root`] with
/// `close_root`) receives the surrogate's bare poles already root-closed
/// ([`RStarPlan::closed_pole_part`]), so the surrogate's root-scope
/// remainders take the remainder-boundary path. Any other value (`0`,
/// `off`): off.
pub fn surrogate_closed_root() -> bool {
    static MODE: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| {
        matches!(
            std::env::var("GAMMALOOP_SURROGATE_CLOSED_ROOT")
                .as_deref()
                .map(str::trim),
            Err(std::env::VarError::NotPresent) | Ok("1")
        )
    });
    *MODE
}

/// `GAMMALOOP_INSERTION_CLOSURE=off|on|check` (default on: production since r6,
/// results/ins4-family, check runs insclose-a..f): how each excised
/// self-energy insertion is integrated in [`PropagatorInsertionEvaluator::prepare`].
/// `off` (an explicit `off`, `0` or empty value) integrates the open tensor
/// (the path before the production default). `on` (also unset)
/// integrates its closed form factors and rebuilds the open tensor exactly
/// ([`super::propagator::integrate_propagator_closed`]: quark `A slash(p)`,
/// gluon `a g + b p p` as one integral with channel markers); a ghost or
/// scalar insertion is already closed and an unrecognized shape declines to
/// the open path. `check` computes both and requires equal normalised
/// coefficients (`INSERTION_CLOSURE_CHECK` lines, an error otherwise),
/// returning the closed one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InsertionClosure {
    Off,
    On,
    Check,
}

impl InsertionClosure {
    /// The mode of a `GAMMALOOP_INSERTION_CLOSURE` value (`None`: unset).
    pub fn parse(value: Option<&str>) -> Result<Self> {
        match value.map(str::trim) {
            Some("" | "off" | "0") => Ok(Self::Off),
            None | Some("on" | "1") => Ok(Self::On),
            Some("check") => Ok(Self::Check),
            Some(other) => Err(eyre!(
                "GAMMALOOP_INSERTION_CLOSURE must be off, on or check, not `{other}`"
            )),
        }
    }

    pub fn from_environment() -> Self {
        static MODE: LazyLock<InsertionClosure> = LazyLock::new(|| {
            match std::env::var("GAMMALOOP_INSERTION_CLOSURE") {
                Ok(value) => InsertionClosure::parse(Some(&value)),
                Err(std::env::VarError::NotPresent) => InsertionClosure::parse(None),
                Err(std::env::VarError::NotUnicode(value)) => Err(eyre!(
                    "GAMMALOOP_INSERTION_CLOSURE={value:?} is not unicode"
                )),
            }
            .unwrap_or_else(|error| panic!("{error:#}"))
        });
        *MODE
    }
}

#[cfg(test)]
#[test]
fn insertion_closure_values_parse() {
    assert_eq!(InsertionClosure::parse(None).unwrap(), InsertionClosure::On);
    for off in ["", "0", "off", " off "] {
        assert_eq!(
            InsertionClosure::parse(Some(off)).unwrap(),
            InsertionClosure::Off
        );
    }
    for on in ["1", "on"] {
        assert_eq!(
            InsertionClosure::parse(Some(on)).unwrap(),
            InsertionClosure::On
        );
    }
    assert_eq!(
        InsertionClosure::parse(Some("check")).unwrap(),
        InsertionClosure::Check
    );
    assert!(InsertionClosure::parse(Some("yes")).is_err());
}

/// The bare poles of the original graph, without its overall factor or projector.
pub struct InsertionPoleResult {
    pub poles: Atom,
    pub insertions: usize,
    pub insertion_loops: usize,
    /// Whether `poles` are root-closed (the pole part of the root metric
    /// times the graph, a closed scalar): feed them to
    /// [`RStarPlan::projected_counterterm_from_closed_pole_part_closing_root`].
    /// Otherwise they keep the open root slots.
    pub closed: bool,
}

/// One selected insertion, expressed in the final surrogate's owner frame.
pub struct PreparedInsertionHole {
    pub vertex: VertexId,
    pub loops: usize,
    pub coefficient: Atom,
    /// The selected support and boundary in the original input graph.
    pub original: PropagatorInsertion,
}

/// A bare insertion identity before the outer propagator is integrated.
///
/// Local coefficient ownership is retained for the ordinary R* fallback above
/// four loops. Group integration instead transports the complete numerator,
/// including the correlated product of all holes, into one physical integrand.
pub struct PreparedInsertion {
    surrogate: FeynmanDiagram,
    holes: Vec<PreparedInsertionHole>,
    maps: Vec<DiagramRelabeling>,
    identity: OnceLock<Result<InsertionSkeletonIdentity, String>>,
    context: Arc<()>,
    extraction_elapsed: Duration,
    /// U7c: what the surrogate's records need from the original diagram
    /// (`None` unless the uncontracted kind `surrogate` is enabled).
    k5: Option<Arc<super::uncontracted::k5::Provenance>>,
}

impl PreparedInsertion {
    pub fn surrogate(&self) -> &FeynmanDiagram {
        &self.surrogate
    }

    pub fn holes(&self) -> &[PreparedInsertionHole] {
        &self.holes
    }

    pub fn loop_count(&self) -> usize {
        self.surrogate.loop_count()
    }

    /// The complete local numerator, excluding global factors and projectors.
    pub fn complete_numerator(&self) -> Atom {
        NumeratorInterface::export(self.surrogate.numerator()).0
    }

    /// Remove only the selected coefficients; all outside local rules survive.
    pub fn skeleton(&self) -> Result<FeynmanDiagram> {
        let one = SymbolicTensor::infer(Atom::one())?;
        Ok(self.surrogate.clone().with_numerator_fragments(
            self.holes.iter().map(|hole| (hole.vertex, one.clone())),
            [],
        )?)
    }

    /// Exact dressed skeleton identity, computed only when grouping is requested.
    pub fn skeleton_identity(&self) -> Result<&InsertionSkeletonIdentity> {
        self.identity
            .get_or_init(|| {
                let build = || -> Result<InsertionSkeletonIdentity> {
                    let skeleton = self.skeleton()?;
                    let graph = skeleton.to_gamma_loop_graph(None, true)?;
                    // The native bridge maps graph data in place, preserving node IDs.
                    let holes = self
                        .holes
                        .iter()
                        .map(|hole| (NodeIndex(hole.vertex.0), hole.loops))
                        .collect::<Vec<_>>();
                    InsertionSkeletonIdentity::new(&graph, skeleton.model(), &holes)
                };
                build().map_err(|error| format!("{error:#}"))
            })
            .as_ref()
            .map_err(|error| eyre!(error.clone()))
    }

    /// Restore an integrated tensor in this surrogate's external frame.
    pub fn restore_atom(&self, expression: &Atom) -> Result<Atom> {
        let mut restored = expression.clone();
        for map in self.maps.iter().rev() {
            restored = map.restore_atom(&restored)?;
        }
        Ok(restored)
    }
}

/// A weighted aggregate, in the first supplied member's surrogate frame.
///
/// This cannot be recorded as an individual member's integral. Every original
/// graph still needs its own proper-UV forest contribution, even when `poles`
/// is zero. No global graph factor or projector is added implicitly.
pub struct GroupedInsertionPoleResult {
    pub poles: Atom,
    pub members: usize,
    pub insertions: usize,
    pub insertion_loops: usize,
}

/// Evaluate proper two-point regions once and replace them by scalar bubbles.
///
/// This implements the *bare-amplitude* identity in arXiv:1801.06084, Eq. (3.2).
/// It must not be used to replace the local counterterm by that of the surrogate:
/// epsilon-dependent insertion coefficients do not commute with the child pole
/// maps. Feed the returned poles to the original plan's
/// [`RStarPlan::counterterm_from_pole_part`] instead.
///
/// All open tensors are retained, including longitudinal terms and fermion
/// chains. No particle-specific self-energy ansatz or transversality is assumed.
/// Settings are fixed for the evaluator's lifetime. Its cache matches dressed
/// particle graphs up to isomorphism and contracted-index renaming, transporting
/// the stored open amplitude into each caller's boundary slots and momenta.
pub struct PropagatorInsertionEvaluator {
    settings: UVgenerationSettings,
    cache: Mutex<HashMap<(String, usize, CountertermKey), CachedIntegral>>,
    counterterm_cache: Option<Arc<CountertermCache>>,
    timings: CacheTimings,
    context: Arc<()>,
    closure: InsertionClosure,
}

impl PropagatorInsertionEvaluator {
    pub fn new(settings: UVgenerationSettings) -> Self {
        Self {
            settings,
            cache: Mutex::new(HashMap::new()),
            counterterm_cache: None,
            timings: CacheTimings::default(),
            context: Arc::new(()),
            closure: InsertionClosure::from_environment(),
        }
    }

    /// Override the insertion closure ([`InsertionClosure`]) taken from
    /// `GAMMALOOP_INSERTION_CLOSURE`.
    pub fn with_insertion_closure(mut self, closure: InsertionClosure) -> Self {
        self.closure = closure;
        self
    }

    /// Enable optional cache-phase accounting without changing cached values.
    pub fn with_timing(mut self, enabled: bool) -> Self {
        self.timings = CacheTimings::new(enabled);
        self
    }

    pub fn cache_timing(&self) -> Option<CacheTimingStats> {
        self.timings.stats()
    }

    /// Share complete UV counterterms with surrounding graph evaluations.
    pub fn with_counterterm_cache(mut self, cache: Arc<CountertermCache>) -> Self {
        self.counterterm_cache = Some(cache);
        self
    }

    /// Extract a maximal vertex-disjoint set, preferring larger regions.
    ///
    /// Nested regions are included in their largest selected ancestor's direct
    /// integral. Each selected region must have at most four loops. A graph
    /// without eligible insertions returns `None`, allowing the full-graph path.
    pub fn pole_part(
        &self,
        diagram: &FeynmanDiagram,
        vakint: &Vakint,
    ) -> Result<Option<InsertionPoleResult>> {
        self.prepare(diagram, vakint)?
            .map(|prepared| self.pole_part_of_prepared(&prepared, vakint))
            .transpose()
    }

    /// [`Self::pole_part`] for a caller that closes the root with the metric
    /// when `close_root` is true: with [`surrogate_closed_root`] on, a
    /// surrogate above four loops whose root closes then returns closed poles
    /// (`closed` in the result). Otherwise, as [`Self::pole_part`].
    pub fn pole_part_closing_root(
        &self,
        diagram: &FeynmanDiagram,
        vakint: &Vakint,
        close_root: bool,
    ) -> Result<Option<InsertionPoleResult>> {
        let close = close_root && surrogate_closed_root();
        self.prepare(diagram, vakint)?
            .map(|prepared| self.pole_part_of_prepared_with(&prepared, vakint, close))
            .transpose()
    }

    /// Prepare the same bare identity as [`Self::pole_part`] without integrating
    /// its outer skeleton. A preparation belongs to this evaluator's settings.
    pub fn prepare(
        &self,
        diagram: &FeynmanDiagram,
        vakint: &Vakint,
    ) -> Result<Option<PreparedInsertion>> {
        ensure!(
            diagram.cuts().is_empty(),
            "insertion evaluation requires an uncut amplitude"
        );
        let mut candidates = diagram.propagator_insertions()?;
        candidates.retain(|candidate| candidate.loops <= 4);
        candidates.sort_by_key(|candidate| std::cmp::Reverse(candidate.loops));
        let mut occupied = BTreeSet::new();
        let selected = candidates
            .into_iter()
            .filter(|candidate| {
                let vertices = diagram.underlying().nodes(&candidate.support);
                if vertices.iter().any(|vertex| occupied.contains(vertex)) {
                    false
                } else {
                    occupied.extend(vertices);
                    true
                }
            })
            .collect::<Vec<_>>();
        if selected.is_empty() {
            return Ok(None);
        }
        let loops = diagram.loop_count();
        let extraction_started = Instant::now();
        let counted = vakint::timing::Snapshot::now();
        let mut surrogate = diagram.clone();
        let mut maps: Vec<DiagramRelabeling> = Vec::new();
        let mut holes: Vec<PreparedInsertionHole> = Vec::new();
        for insertion in &selected {
            let excised = diagram.excise_with_mapping(&insertion.support, &[])?;
            let coefficient =
                self.normalized_insertion(&excised.diagram, vakint, loops, diagram.name())?;
            let vanishes = coefficient.is_zero();
            // Slashes are opened after the freshening, with process-unique
            // symbol indices: freshened `dummy(N)` indices are not Spenso indices.
            let mut coefficient = open_slashes(&freshen_insertion_indices(
                &excised.relabeling.restore_atom(&coefficient)?,
                &[
                    NumeratorInterface::export(surrogate.numerator()).0,
                    NumeratorInterface::export(surrogate.projector()).0,
                    NumeratorInterface::export(surrogate.numerator_prefactor()).0,
                ],
            )?);
            let mut current = insertion.clone();
            for map in &maps {
                coefficient = map.translate_atom(&coefficient, None)?;
                current = PropagatorInsertion {
                    support: map.translate_subgraph(&current.support)?,
                    boundary: [
                        map.translate_half_edge(current.boundary[0])?,
                        map.translate_half_edge(current.boundary[1])?,
                    ],
                    loops: current.loops,
                };
            }
            let particle = surrogate
                .edges()
                .find(|(id, _, _)| *id == current.boundary[0].edge)
                .ok_or_else(|| eyre!("missing insertion boundary edge"))?
                .2
                .particle;
            let replacement = surrogate.scalar_bubble_surrogate(&current, particle)?;
            surrogate = replacement.with_coefficient(&SymbolicTensor::infer(coefficient)?)?;
            for hole in &mut holes {
                hole.vertex = replacement
                    .relabeling
                    .vertex(hole.vertex)
                    .ok_or_else(|| eyre!("a disjoint insertion lost its coefficient vertex"))?;
                hole.coefficient = replacement
                    .relabeling
                    .translate_atom(&hole.coefficient, None)?;
            }
            let coefficient = surrogate
                .vertices()
                .find(|(id, _)| *id == replacement.coefficient_vertex)
                .ok_or_else(|| eyre!("missing insertion coefficient vertex"))?
                .1;
            holes.push(PreparedInsertionHole {
                vertex: replacement.coefficient_vertex,
                loops: insertion.loops,
                coefficient: NumeratorInterface::export(&coefficient.numerator).0,
                original: insertion.clone(),
            });
            maps.push(replacement.relabeling);
            // The complete bare product already vanishes. Preserve the old
            // early exit without evaluating further insertion coefficients.
            if vanishes {
                break;
            }
        }
        if super::plan::uv_children_log() {
            // The insertions' own integrals (normalisation), outside every R*
            // evaluation of the diagram.
            eprintln!(
                "INSERTION_PREPARE graph={} insertions={} seconds={:.3}{}",
                diagram.name(),
                holes.len(),
                extraction_started.elapsed().as_secs_f64(),
                counted.elapsed().fields()
            );
        }
        let k5 = super::uncontracted::k5::provenance(
            diagram,
            holes
                .iter()
                .map(|hole| (hole.vertex.0, hole.loops, &hole.original.support)),
        )?;
        Ok(Some(PreparedInsertion {
            surrogate,
            holes,
            maps,
            identity: OnceLock::new(),
            context: self.context.clone(),
            extraction_elapsed: extraction_started.elapsed(),
            k5,
        }))
    }

    /// Integrate one prepared identity, retaining local ownership above four
    /// loops and returning the original graph's external tensor frame.
    pub fn pole_part_of_prepared(
        &self,
        prepared: &PreparedInsertion,
        vakint: &Vakint,
    ) -> Result<InsertionPoleResult> {
        self.pole_part_of_prepared_with(prepared, vakint, false)
    }

    /// [`Self::pole_part_of_prepared`]; with `close_root`, a surrogate above
    /// four loops whose root closes is evaluated root-closed
    /// ([`RStarPlan::closed_pole_part`]) and the result is flagged `closed`.
    /// Through four loops, and for a vanishing numerator, the poles stay open
    /// (`closed` false): the caller closes them as before.
    fn pole_part_of_prepared_with(
        &self,
        prepared: &PreparedInsertion,
        vakint: &Vakint,
        close_root: bool,
    ) -> Result<InsertionPoleResult> {
        self.check_preparation(prepared)?;
        let surrogate = &prepared.surrogate;
        let loops = prepared.loop_count();
        let surrogate_started = Instant::now();
        let mut closed = false;
        // Through four loops the physical propagator backend gives KG directly.
        // Reconstructing it through R* would compute a second, unnecessary UV
        // forest: the caller still needs the original graph's forest for Delta.
        let poles = if prepared.complete_numerator().is_zero() {
            Ok(Atom::Zero)
        } else if loops <= 4 {
            let graph = surrogate.to_gamma_loop_graph(None, true)?;
            // U7c: the K5 direct record (kind surrogate; off: just the call)
            super::uncontracted::k5::direct(
                prepared.k5.as_ref(),
                &graph,
                surrogate.model(),
                vakint,
                &self.settings,
                loops,
                || {
                    integrate_propagator(
                        &graph,
                        surrogate.model(),
                        vakint,
                        &self.settings,
                        loops as i64,
                    )
                },
            )
        } else {
            let graph = surrogate.to_gamma_loop_graph(None, true)?;
            let mut plan =
                RStarPlan::with_automatic_edge(&graph, surrogate.model(), &self.settings)?;
            if let Some(cache) = &self.counterterm_cache {
                plan = plan.with_counterterm_cache(cache.clone());
            }
            // GAMMALOOP_SURROGATE_CLOSED_ROOT: the closed surrogate root, when
            // its numerator has the two open external slots of the closure.
            closed = close_root && plan.closes_root(&graph);
            if close_root && !closed {
                eprintln!(
                    "SURROGATE_CLOSED_ROOT graph={} closed=false reason=root_closure",
                    graph.name
                );
            }
            // U7c: the surrogate's records (kind surrogate; off: just the call)
            super::uncontracted::k5::scoped(prepared.k5.as_ref(), &graph.name, closed, &self.settings, || {
                super::expression::in_dressed_scope(|| {
                    if closed {
                        plan.closed_pole_part(&graph, surrogate.model(), vakint, &self.settings)
                    } else {
                        plan.pole_part(&graph, surrogate.model(), vakint, &self.settings)
                    }
                })
            })
        }
        .wrap_err("while evaluating the scalar-bubble surrogate's bare poles")?;
        // Closed poles are a scalar in the surrogate's external momentum: the
        // owner map still carries that momentum label to the original's.
        let poles = prepared.restore_atom(&poles)?;
        if closed {
            eprintln!(
                "SURROGATE_CLOSED_ROOT graph={} closed=true loops={loops} seconds={:.3}",
                surrogate.name(),
                surrogate_started.elapsed().as_secs_f64()
            );
        }
        tracing::info!(
            uv = true,
            rstar = true,
            timing = true,
            insertions = prepared.holes.len(),
            insertion_loops = prepared.holes.iter().map(|hole| hole.loops).sum::<usize>(),
            extraction_seconds = prepared.extraction_elapsed.as_secs_f64(),
            surrogate_seconds = surrogate_started.elapsed().as_secs_f64(),
            "Evaluated bare poles with propagator insertions extracted"
        );
        Ok(InsertionPoleResult {
            poles,
            insertions: prepared.holes.len(),
            insertion_loops: prepared.holes.iter().map(|hole| hole.loops).sum(),
            closed,
        })
    }

    /// Integrate one verified skeleton with the sum of the supplied weighted
    /// complete numerators. Products at different sites remain correlated:
    /// this sums actual members and never forms a Cartesian product of holes.
    ///
    /// Weights must be momentum-independent scalars regular at epsilon=0 and
    /// include any desired original graph factors exactly once. At most four
    /// loops are admitted: placing a multisite sum on one local vertex is only
    /// valid for direct physical integration, not for a surrogate UV forest.
    /// Call [`Self::pole_part_of_prepared`] separately above four loops.
    pub fn weighted_pole_part(
        &self,
        members: &[(&PreparedInsertion, Atom)],
        vakint: &Vakint,
    ) -> Result<GroupedInsertionPoleResult> {
        let (representative, _) = members
            .first()
            .ok_or_else(|| eyre!("an insertion group needs at least one member"))?;
        // U7c: a grouped skeleton carries the weights inside its integrand
        // (U7-PLAN requirement 12): it never builds K5 records
        ensure!(
            !super::uncontracted::switch().enabled(super::uncontracted::Kind::Surrogate),
            "{}: the grouped insertion skeleton is not supported with the kind surrogate (its weights are inside the integrand: no record); run without skeleton grouping",
            super::uncontracted::VARIABLE
        );
        self.check_preparation(representative)?;
        let loops = representative.loop_count();
        ensure!(
            (1..=4).contains(&loops),
            "grouped insertion integration supports at most four loops; retain per-member local ownership above four loops"
        );
        let identity = representative.skeleton_identity()?;
        let mut numerator = Atom::Zero;
        let mut insertions = 0;
        let mut insertion_loops = 0;
        for (member, weight) in members {
            self.check_preparation(member)?;
            ensure!(
                member.loop_count() == loops,
                "insertion group has different parent loop orders"
            );
            Self::validate_weight(weight)?;
            let source_identity = member.skeleton_identity()?;
            ensure!(
                identity.key() == source_identity.key(),
                "insertion group contains different dressed skeletons"
            );
            numerator +=
                weight * identity.transport_from(source_identity, &member.complete_numerator())?;
            insertions += member.holes.len();
            insertion_loops += member.holes.iter().map(|hole| hole.loops).sum::<usize>();
        }
        // `with_numerator` keeps the complete tensor as one expression. This
        // graph is used solely by the physical integrator, never by R* forests.
        let poles = if numerator.is_zero() {
            Atom::Zero
        } else {
            let diagram = representative
                .surrogate
                .clone()
                // Identity transport restores source-owned indices, including
                // functional dummy(n) labels. Re-admit that notation through
                // the same reversible boundary as a retained graph numerator.
                .with_numerator(NumeratorInterface::admit(numerator, None)?)?;
            let graph = diagram.to_gamma_loop_graph(None, true)?;
            integrate_propagator(
                &graph,
                diagram.model(),
                vakint,
                &self.settings,
                loops as i64,
            )
            .wrap_err("while integrating a weighted insertion skeleton")?
        };
        Ok(GroupedInsertionPoleResult {
            poles,
            members: members.len(),
            insertions,
            insertion_loops,
        })
    }

    fn check_preparation(&self, prepared: &PreparedInsertion) -> Result<()> {
        ensure!(
            Arc::ptr_eq(&self.context, &prepared.context),
            "insertion preparation belongs to a different evaluator/settings context"
        );
        Ok(())
    }

    fn validate_weight(weight: &Atom) -> Result<()> {
        use spenso::structure::partial::PartialStructureExt;
        ensure!(
            ![GS.emr_mom, GS.loop_mom, GS.external_mom]
                .into_iter()
                .any(|momentum| weight.contains_symbol(momentum)),
            "insertion group weights must not depend on momenta"
        );
        ensure!(
            SymbolicTensor::infer(weight.clone())?
                .structure()
                .logical_slots()
                .is_empty(),
            "insertion group weights must be scalars"
        );
        // MATRIX-R45: a closed trace (a scalar with summed dummy indices,
        // e.g. g_ab g^ab = d or tr(γ_a γ^a) = 4d) in a denominator hides a
        // pole from the substitution of d below. A weight that holds a tensor
        // structure is normalised by the IRH-DEPTH certificate first (closed
        // traces contracted, d = 4−2ε; an unresolved one is refused). Weights
        // without a tensor structure (all production weights: colour, nf and
        // the overall factor) keep the plain substitution.
        let normalised;
        let weight = if super::regularity::has_tensor_structure(weight) {
            normalised = super::regularity::normalise(
                super::regularity::Factor::InsertionCoefficient,
                weight,
            )
            .map_err(|refused| {
                eyre!("insertion group weights must be provably regular in epsilon: {refused}")
            })?;
            &normalised
        } else {
            weight
        };
        let epsilon_weight = weight
            .replace(Atom::var(GS.dim))
            .with(Atom::num(4) - Atom::num(2) * Atom::var(GS.dim_epsilon));
        ensure!(
            epsilon_weight
                .series(GS.dim_epsilon, 0, -1)?
                .to_atom()
                .expand()
                .together()
                .is_zero(),
            "insertion group weights must be regular in epsilon; singular weights require additional bare finite terms"
        );
        Ok(())
    }

    fn normalized_insertion(
        &self,
        insertion: &FeynmanDiagram,
        vakint: &Vakint,
        parent_loops: usize,
        parent: &str,
    ) -> Result<Atom> {
        let loops = insertion.loop_count();
        // GAMMALOOP_INSERTION_CLOSURE: the closed form factors, when the
        // insertion has an open Lorentz or spinor structure to close.
        let shape = match self.closure {
            InsertionClosure::Off => None,
            InsertionClosure::On | InsertionClosure::Check => {
                let graph = insertion.to_gamma_loop_graph(None, true)?;
                match super::propagator::insertion_shape(&graph)? {
                    Ok(super::propagator::InsertionShape::Scalar) => {
                        eprintln!(
                            "INSERTION_CLOSURE graph={parent} kind=scalar loops={loops} closure=none"
                        );
                        None
                    }
                    Ok(shape) => Some(shape),
                    Err(reason) => {
                        eprintln!(
                            "INSERTION_CLOSURE graph={parent} kind=open loops={loops} closure=declined reason={reason}"
                        );
                        None
                    }
                }
            }
        };
        let Some(shape) = shape else {
            return self.normalized_insertion_with(insertion, vakint, parent_loops, false, parent);
        };
        let started = Instant::now();
        let closed =
            self.normalized_insertion_with(insertion, vakint, parent_loops, true, parent)?;
        if self.closure != InsertionClosure::Check {
            return Ok(closed);
        }
        let closed_seconds = started.elapsed().as_secs_f64();
        let started = Instant::now();
        let open =
            self.normalized_insertion_with(insertion, vakint, parent_loops, false, parent)?;
        let open_seconds = started.elapsed().as_secs_f64();
        let (identical, equal) = insertion_coefficients_agree(&closed, &open)?;
        eprintln!(
            "INSERTION_CLOSURE_CHECK graph={parent} kind={} loops={loops} terms={parent_loops} identical={identical} equal={equal} zero={} closed_seconds={closed_seconds:.3} open_seconds={open_seconds:.3}",
            shape.kind(),
            open.is_zero(),
        );
        if !equal {
            // A failed comparison, which no fallback declines.
            return Err(super::ir_dirac::diagnostic_failure(eyre!(
                "INSERTION_CLOSURE_CHECK: the closed {loops}-loop {} insertion of {parent} differs from the open tensor\nclosed: {closed}\nopen: {open}",
                shape.kind()
            )));
        }
        Ok(closed)
    }

    /// The normalised coefficient of one insertion, its integral open
    /// (`closed` false) or through its closed form factors.
    fn normalized_insertion_with(
        &self,
        insertion: &FeynmanDiagram,
        vakint: &Vakint,
        parent_loops: usize,
        closed: bool,
        parent: &str,
    ) -> Result<Atom> {
        let loops = insertion.loop_count();
        // An N-loop surrogate has at most N poles. Its regular coefficient is
        // therefore needed through eps^(N-1), requiring the L-loop insertion
        // through eps^(N-L-1): precisely N terms in Vakint's convention.
        let integral =
            self.integrate_cached_with(insertion, vakint, parent_loops, closed, parent)?;
        if integral.is_zero() {
            return Ok(Atom::Zero);
        }
        let bubble = scalar_bubble(insertion)?;
        let bubble_integral = self.integrate_cached(&bubble, vakint, parent_loops)?;
        let external = insertion
            .edges()
            .find_map(|(id, _, data)| data.external.as_ref().map(|_| id))
            .ok_or_else(|| eyre!("excised insertion has no external momentum"))?;
        let bubble_integral = bubble_integral
            .replace(function!(GS.emr_mom, 0, W_.x___))
            .with(function!(GS.emr_mom, external.0, W_.x___));
        // IRH-DEPTH round 8b (C1): d = 4 − 2ε before the truncation at
        // ε^(N−1) (the closed form factors carry 1/(d−1), an insertion vertex
        // may carry any scalar d): the coefficient is a pure ε polynomial.
        let ratio = integral / bubble_integral.pow(loops);
        let coefficient = crate::uv::approx::integrated::Integrated::regulated_for_expansion(
            &ratio,
            "insertion_coefficient_truncation",
        )?
        .series(GS.dim_epsilon, 0, parent_loops as i64 - 1)?
        .to_atom()
        .expand()
        .together();
        ensure!(
            coefficient
                .series(GS.dim_epsilon, 0, -1)?
                .to_atom()
                .expand()
                .is_zero(),
            "normalized insertion has epsilon poles; scalar-bubble bare poles alone are insufficient"
        );
        ensure!(
            !coefficient.contains_symbol(Symbol::LOG) && !coefficient.contains_symbol(GS.mu_r_sq),
            "insertion scale did not cancel against the scalar bubbles: {coefficient}"
        );
        // The surrogate R* forest requires local polynomial vertices. A more
        // general insertion can produce inverse external invariants, which
        // must be encoded as line powers before that representation is valid.
        let mut polynomial = true;
        let _ = coefficient.replace_map_bottom_up(|term, _, _| {
            if let AtomView::Pow(power) = term
                && power.get_base().contains_symbol(GS.emr_mom)
                && i64::try_from(power.get_exp())
                    .ok()
                    .is_none_or(|exponent| exponent < 0)
            {
                polynomial = false;
            }
        });
        ensure!(
            polynomial,
            "normalized insertion is not polynomial in its external momentum"
        );
        Ok(coefficient)
    }

    fn integrate_cached(
        &self,
        diagram: &FeynmanDiagram,
        vakint: &Vakint,
        terms: usize,
    ) -> Result<Atom> {
        self.integrate_cached_with(diagram, vakint, terms, false, "")
    }

    /// [`Self::integrate_cached`]; with `closed`, the open tensor is rebuilt
    /// from its closed form factors
    /// ([`super::propagator::integrate_propagator_closed`]), in its own cache
    /// entries. An unrecognized shape falls back to the open integral.
    fn integrate_cached_with(
        &self,
        diagram: &FeynmanDiagram,
        vakint: &Vakint,
        terms: usize,
        closed: bool,
        parent: &str,
    ) -> Result<Atom> {
        let mut timing = self.timings.start();
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let boundary = graph.external_filter::<SuBitGraph>();
        let support = graph.no_dummy().subtract(&boundary);
        let identity = CountertermIdentity::new(&graph, &support, diagram.model())?;
        // The complete bare tensor is stored, not just its poles. Expansion
        // depth belongs to the key; lower-loop finite terms enter parent poles.
        let fingerprint = format!("{:?}", diagram.model().fingerprint());
        let key = (
            if closed {
                format!("{fingerprint}#insertion-closure")
            } else {
                fingerprint
            },
            terms,
            identity.key.clone(),
        );
        timing.lap(CacheTimingPhase::Key);
        let entry = {
            let mut cache = self
                .cache
                .lock()
                .map_err(|_| eyre!("insertion cache lock poisoned"))?;
            timing.lap(CacheTimingPhase::MapLockWait);
            let entry = cache.entry(key).or_default().clone();
            drop(cache);
            timing.lap(CacheTimingPhase::MapLockHold);
            entry
        };
        let canonical = timing
            .resolve(&entry, |timing| {
                let started = Instant::now();
                let integral = if closed {
                    let counted = vakint::timing::Snapshot::now();
                    super::propagator::integrate_propagator_closed(
                        &graph,
                        diagram.model(),
                        vakint,
                        &self.settings,
                        terms as i64,
                    )
                    .and_then(|closed| match closed {
                        Ok((integral, shape)) => {
                            eprintln!(
                                "INSERTION_CLOSURE graph={parent} kind={} loops={} terms={terms} closure=closed seconds={:.3}{}",
                                shape.kind(),
                                diagram.loop_count(),
                                started.elapsed().as_secs_f64(),
                                counted.elapsed().fields()
                            );
                            Ok(integral)
                        }
                        Err(reason) => {
                            eprintln!(
                                "INSERTION_CLOSURE graph={parent} kind=open loops={} terms={terms} closure=declined reason={reason}",
                                diagram.loop_count()
                            );
                            integrate_propagator(
                                &graph,
                                diagram.model(),
                                vakint,
                                &self.settings,
                                terms as i64,
                            )
                        }
                    })
                } else {
                    integrate_propagator(
                        &graph,
                        diagram.model(),
                        vakint,
                        &self.settings,
                        terms as i64,
                    )
                };
                timing.lap(CacheTimingPhase::Initializer);
                let result = integral
                    .and_then(|integral| identity.canonicalize_delta(&integral))
                    .map_err(|error| format!("{error:#}"));
                timing.lap(CacheTimingPhase::Store);
                tracing::info!(
                    uv = true,
                    rstar = true,
                    timing = true,
                    loops = diagram.loop_count(),
                    terms = terms,
                    seconds = started.elapsed().as_secs_f64(),
                    "Computed an uncached bare propagator insertion"
                );
                result
            })
            .clone()
            .map_err(|error| eyre!(error))?;
        // An isomorphism can send the stored momentum to either external leg.
        // Restore the direct integrator's one-momentum convention before the
        // amplitude is divided by its scalar-bubble normalization. This also
        // identifies logarithms and finite terms, not just polynomial poles.
        timing.restart();
        let result = RStarIntegrand::physical_external_routing(&graph, &boundary).and_then(
            |(_, routing)| {
                Ok(super::linearize::linearize(
                    identity
                        .restore_delta(&canonical)?
                        .replace_multiple(&routing),
                ))
            },
        );
        timing.lap(CacheTimingPhase::Restore);
        result
    }
}

/// Whether two normalised insertion coefficients in the same caller frame
/// agree: `identical` as expressions, `equal` after gamma and colour
/// simplification and canonical dummy relabelling of their difference.
pub(super) fn insertion_coefficients_agree(closed: &Atom, open: &Atom) -> Result<(bool, bool)> {
    use idenso::{
        color::ColorSimplifySettings, dirac::GammaSimplifySettings, tensor::AlgebraSettings,
    };
    let difference = (closed - open).expand();
    if difference.is_zero() {
        return Ok((true, true));
    }
    if difference.together().is_zero() {
        return Ok((false, true));
    }
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(difference.as_view())?)?
        .undo_trace()?
        .undo_chain()?
        .simplify_algebra(&AlgebraSettings {
            gamma: Some(GammaSimplifySettings::default()),
            color: Some(ColorSimplifySettings::default()),
            ..Default::default()
        })?
        .to_dots()?;
    let reduced = cooking.uncook(tensor.expression().as_view()).expand();
    if reduced.is_zero() {
        return Ok((false, true));
    }
    let canonical = SymbolicTensor::infer(cooking.try_cook(reduced.as_view())?)?.canonize()?;
    let canonical = cooking.uncook(canonical.expression().as_view()).expand();
    Ok((false, canonical.together().is_zero()))
}

fn scalar_bubble(source: &FeynmanDiagram) -> Result<FeynmanDiagram> {
    let particle = source
        .edges()
        .find(|(_, _, edge)| !edge.is_dummy)
        .ok_or_else(|| eyre!("insertion has no propagators"))?
        .2
        .particle;
    let directed = !source.model().particle_is_self_conjugate(particle);
    let mut builder = FeynmanDiagram::builder(source.model_arc(), "insertion_normalization");
    let mut vertex = |name: &str| {
        builder.add_vertex(DiagramVertex {
            name: name.into(),
            interaction: None,
            numerator: SymbolicTensor::infer(Atom::one()).unwrap(),
        })
    };
    let left = vertex("left");
    let right = vertex("right");
    let mut incoming = DiagramEdge::new(particle, directed);
    incoming.external = Some(ExternalLeg {
        name: "in".into(),
        index: 0,
        state: ExternalState::Incoming,
        connection: 0,
    });
    builder.add_edge(None, left, incoming)?;
    builder.add_edge(left, right, DiagramEdge::new(particle, directed))?;
    builder.add_edge(left, right, DiagramEdge::new(particle, directed))?;
    let mut outgoing = DiagramEdge::new(particle, directed);
    outgoing.external = Some(ExternalLeg {
        name: "out".into(),
        index: 1,
        state: ExternalState::Outgoing,
        connection: 1,
    });
    builder.add_edge(right, None, outgoing)?;
    Ok(builder.build()?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{initialisation::test_initialise, model::Model};
    use feynkit_graph::VertexId;
    use idenso::representations::{Bispinor, ColorAdjoint};
    use spenso::{
        network::{library::symbolic::ETS, tags::SPENSO_TAG},
        structure::{
            partial::PartialStructureExt,
            representation::{Minkowski, RepName},
        },
    };
    use std::collections::HashSet;

    /// A quark self-energy coefficient `c slash(Q_7)` (idenso's slashed form)
    /// declines the UV-child truncation; [`open_slashes`] gives the same
    /// tensor with the momentum only in `Q`, which the truncation grades.
    #[test]
    fn open_slashes_make_a_slashed_insertion_coefficient_polynomial() -> Result<()> {
        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        // Symbol indices parse both as Spenso indices and as `Aind`.
        let i = spin.to_symbolic([Atom::from(Aind::Symbol(symbol!("open_slashes_i")))]);
        let j = spin.to_symbolic([Atom::from(Aind::Symbol(symbol!("open_slashes_j")))]);
        let slash = |edge: i64| {
            function!(
                idenso::dirac::AGS.gamma,
                &i,
                &j,
                function!(GS.emr_mom, edge, mink.to_symbolic([]))
            )
        };
        let c = Atom::var(symbol!("c_open_slashes"));
        let coefficient = (&c * slash(7) + Atom::num(3) * slash(7)).expand();
        let roots = BTreeSet::from([7_usize]);
        assert!(!super::super::expression::polynomial_in_roots(
            coefficient.as_view(),
            &roots
        ));
        // The production order: freshen, then open, then the surrogate's
        // Spenso parse (`with_coefficient`), which rejected `dummy(N)` indices.
        let freshened = freshen_insertion_indices(&coefficient, &[])?;
        let opened = open_slashes(&freshened);
        SymbolicTensor::infer(coefficient.clone())
            .wrap_err_with(|| format!("slashed input: {coefficient}"))?;
        SymbolicTensor::infer(opened.clone())
            .wrap_err_with(|| format!("freshened: {freshened}\nopened: {opened}"))?;
        assert!(super::super::expression::polynomial_in_roots(
            opened.as_view(),
            &roots
        ));
        let mut slashed = false;
        opened.visitor(&mut |view| {
            if let AtomView::Fun(gamma) = view
                && gamma.get_symbol() == idenso::dirac::AGS.gamma
                && matches!(gamma.get(2), AtomView::Fun(momentum) if momentum.get_symbol() == GS.emr_mom)
            {
                slashed = true;
            }
            !slashed
        });
        assert!(!slashed, "a slash survived: {opened}");
        // Contracting the opened index restores the slashed tensor.
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let contract = |atom: &Atom| -> Result<Atom> {
            let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(atom.as_view())?)?
                .simplify_algebra(&idenso::tensor::AlgebraSettings {
                    gamma: Some(idenso::dirac::GammaSimplifySettings::default()),
                    ..Default::default()
                })?
                .to_dots()?;
            Ok(cooking.uncook(tensor.expression().as_view()).expand())
        };
        assert_eq!(contract(&opened)?, contract(&coefficient)?);
        // Without a slash the coefficient is unchanged.
        assert_eq!(open_slashes(&c), c);
        Ok(())
    }

    #[test]
    fn cached_open_dirac_and_color_tensors_get_independent_dummy_indices() -> Result<()> {
        test_initialise()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let color = ColorAdjoint {}.new_rep(8);
        let mu = mink.to_symbolic([Atom::from(Aind::UVTerm(0, 1))]);
        let nu = mink.to_symbolic([Atom::from(Aind::UVTerm(0, 2))]);
        let internal_spin = spin.to_symbolic([Atom::from(Aind::UVTerm(0, 3))]);
        let internal_color = color.to_symbolic([Atom::from(Aind::UVTerm(0, 4))]);
        let cached = function!(
            idenso::dirac::AGS.gamma,
            spin.to_symbolic([Atom::from(Aind::Hedge(8, 0))]),
            &internal_spin,
            &mu
        ) * function!(
            idenso::dirac::AGS.gamma,
            &internal_spin,
            spin.to_symbolic([Atom::from(Aind::Hedge(9, 0))]),
            &nu
        ) * function!(GS.emr_mom, 0, &mu)
            * function!(GS.emr_mom, 0, &nu)
            * function!(
                ETS.metric,
                color.to_symbolic([Atom::from(Aind::Hedge(8, 1))]),
                &internal_color
            )
            * function!(
                ETS.metric,
                &internal_color,
                color.to_symbolic([Atom::from(Aind::Hedge(9, 1))])
            );
        // Simulate restoring one cache entry onto two disjoint graph regions.
        let second = cached.replace_map(|term, _, out| {
            if let AtomView::Fun(fun) = term
                && fun.get_symbol() == GS.hedgeaind
                && let Ok(Aind::Hedge(owner, local)) = Aind::try_from(term)
            {
                // Local slot zero is serialized as hedge(owner), whereas the
                // color slot uses hedge(owner, 1); both must move together.
                **out = Atom::from(Aind::Hedge(owner + 10, local));
            }
        });
        for local in [0, 1] {
            for owner in [8, 9] {
                assert!(!second.contains(Atom::from(Aind::Hedge(owner, local))));
                assert!(second.contains(Atom::from(Aind::Hedge(owner + 10, local))));
            }
        }
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let mut tensors = Vec::new();
        let mut private_names = Vec::new();
        let mut public_ports = Vec::new();
        // Imported dummy labels need not have passed through the runtime Aind
        // parser. Reserve all surrounding fragments before allocating names.
        let imported_max = 1_000_002;
        let mut context = (1_000_000..=imported_max)
            .map(|index| {
                function!(
                    GS.emr_mom,
                    31,
                    mink.to_symbolic([Atom::from(Aind::Dummy(index))])
                )
            })
            .collect::<Vec<_>>();
        for source in [&cached, &second] {
            let fresh = freshen_insertion_indices(source, &context)?;
            assert!(!fresh.contains_symbol(GS.uvaind));
            let mut names = HashSet::new();
            fresh.as_view().visitor(&mut |term| {
                if let AtomView::Fun(fun) = term
                    && fun.get_symbol() == GS.dummyaind
                {
                    names.insert(Aind::try_from(term).expect("runtime-compatible dummy"));
                    return false;
                }
                true
            });
            assert_eq!(names.len(), 4);
            assert!(
                names
                    .iter()
                    .all(|index| { matches!(index, Aind::Dummy(index) if *index > imported_max) })
            );
            let original = SymbolicTensor::infer(cooking.try_cook(source.as_view())?)?;
            context.push(fresh.clone());
            let fresh = SymbolicTensor::infer(cooking.try_cook(fresh.as_view())?)?;
            let ports = fresh
                .structure()
                .logical_slots()
                .into_iter()
                .collect::<HashSet<_>>();
            assert_eq!(ports.len(), 4);
            assert_eq!(
                original
                    .structure()
                    .logical_slots()
                    .into_iter()
                    .collect::<HashSet<_>>(),
                ports
            );
            public_ports.push(ports);
            private_names.push(names);
            tensors.push(fresh);
        }
        assert!(private_names[0].is_disjoint(&private_names[1]));
        assert!(public_ports[0].is_disjoint(&public_ports[1]));
        let product = tensors[0].multiply(&tensors[1])?;
        assert_eq!(product.structure().logical_slots().len(), 8);
        Ok(())
    }

    fn serial_bubbles(count: usize) -> Result<FeynmanDiagram> {
        serial_bubbles_with(count, Atom::one() + Atom::var(GS.dim_epsilon))
    }

    /// [`serial_bubbles`] with each insertion's numerator `scalar · (Q·k)`.
    fn serial_bubbles_with(count: usize, scalar: Atom) -> Result<FeynmanDiagram> {
        let model = Arc::new(Model::from_json(
            r#"{
            "name":"insertion_identity","restriction":null,"orders":[],
            "parameters":[{"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null}],
            "particles":[{"pdg_code":25,"name":"phi","antiname":"phi","spin":1,"color":1,"mass":"ZERO","width":"ZERO","texname":"phi","antitexname":"phi","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0}],
            "propagators":[],"lorentz_structures":[],"couplings":[],"vertex_rules":[]
        }"#,
        )?);
        let particle = model.particle_id("phi")?;
        let mut builder = FeynmanDiagram::builder(model, "serial_dressed_bubbles");
        for vertex in 0..2 + 2 * count {
            builder.add_vertex(DiagramVertex {
                name: format!("v{vertex}"),
                interaction: None,
                numerator: SymbolicTensor::infer(Atom::one())?,
            });
        }
        for (index, (source, target, state)) in [
            (None, Some(VertexId(0)), ExternalState::Incoming),
            (Some(VertexId(1)), None, ExternalState::Outgoing),
        ]
        .into_iter()
        .enumerate()
        {
            let mut edge = DiagramEdge::new(particle, false);
            edge.external = Some(ExternalLeg {
                name: format!("external_{index}"),
                index,
                state,
                connection: index,
            });
            builder.add_edge(source, target, edge)?;
        }
        builder.add_edge(VertexId(0), VertexId(1), DiagramEdge::new(particle, false))?;
        let mut previous = VertexId(0);
        let mut numerators = Vec::new();
        let vector = |edge: usize| {
            function!(
                GS.emr_mom,
                edge,
                Minkowski {}.new_rep(GS.dim).to_symbolic([])
            )
        };
        for index in 0..count {
            let left = VertexId(2 + 2 * index);
            let right = VertexId(3 + 2 * index);
            let boundary = builder.add_edge(previous, left, DiagramEdge::new(particle, false))?;
            let internal = builder.add_edge(left, right, DiagramEdge::new(particle, false))?;
            builder.add_edge(left, right, DiagramEdge::new(particle, false))?;
            // A real loop-momentum numerator, localized within each insertion.
            // Its integral is q²/2 times the bubble. The evanescent factor must
            // survive both normalization and the original child pole maps.
            numerators.push((
                left,
                SymbolicTensor::infer(
                    &scalar * function!(SPENSO_TAG.dot, vector(boundary.0), vector(internal.0)),
                )?,
            ));
            previous = right;
        }
        builder.add_edge(previous, VertexId(1), DiagramEdge::new(particle, false))?;
        Ok(builder.build()?.with_numerator_fragments(numerators, [])?)
    }

    fn check_serial(count: usize) -> Result<()> {
        test_initialise()?;
        let diagram = serial_bubbles(count)?;
        let mut settings = super::super::tests::settings();
        settings.vakint.normalization = "GScheme".into();
        let evaluator = PropagatorInsertionEvaluator::new(settings.clone());
        let vakint = Vakint::new()?;
        let extracted = evaluator
            .pole_part(&diagram, &vakint)?
            .expect("proper insertions");
        assert_eq!(extracted.insertions, count);
        assert_eq!(extracted.insertion_loops, count);
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let direct = integrate_propagator(
            &graph,
            diagram.model(),
            &vakint,
            &settings,
            diagram.loop_count() as i64,
        )?;
        super::super::tests::assert_same_counterterm(&extracted.poles, &direct)?;
        let plan = RStarPlan::with_automatic_edge(&graph, diagram.model(), &settings)?;
        let recovered = plan.counterterm_from_pole_part(
            &graph,
            diagram.model(),
            &vakint,
            &settings,
            &extracted.poles,
        )?;
        let full = plan.counterterm(&graph, diagram.model(), &vakint, &settings)?;
        super::super::tests::assert_same_counterterm(&recovered, &full)?;
        // Cache reuse must not change tensor owners or the regulator depth.
        let repeated = evaluator.pole_part(&diagram, &vakint)?.unwrap();
        super::super::tests::assert_same_counterterm(&repeated.poles, &direct)
    }

    /// IRH-DEPTH round 8b (C1): a factor equal to −1 at d = 4 − 2ε inside a
    /// propagator insertion, spelled `1/(1+(d−4)/ε)` or `ε/(ε+d−4)`, reaches
    /// the surrogate's normalised insertion coefficient (its truncation at
    /// ε^(N−1)) and the original counterterm: both equal the literal −1's.
    #[test]
    fn round8b_an_insertion_coefficient_spelling_equals_the_literal() -> Result<()> {
        test_initialise()?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let dimension = Atom::var(GS.dim);
        let mut settings = super::super::tests::settings();
        settings.vakint.normalization = "GScheme".into();
        let vakint = Vakint::new()?;
        let run = |scalar: Atom| -> Result<(Atom, Atom)> {
            let diagram = serial_bubbles_with(2, scalar)?;
            let evaluator = PropagatorInsertionEvaluator::new(settings.clone());
            let extracted = evaluator
                .pole_part(&diagram, &vakint)?
                .expect("proper insertions");
            let graph = diagram.to_gamma_loop_graph(None, true)?;
            let plan = RStarPlan::with_automatic_edge(&graph, diagram.model(), &settings)?;
            let recovered = plan.counterterm_from_pole_part(
                &graph,
                diagram.model(),
                &vakint,
                &settings,
                &extracted.poles,
            )?;
            Ok((extracted.poles, recovered))
        };
        let (literal_poles, literal_counterterm) = run(Atom::num(-1))?;
        let cross = Atom::one() / (Atom::one() + (&dimension - Atom::num(4)) / &epsilon);
        let form = &epsilon / (&epsilon + &dimension - Atom::num(4));
        let mut failures = Vec::new();
        for (name, scalar) in [("cross", cross), ("form", form)] {
            let (poles, counterterm) = run(scalar)?;
            let difference = |a: &Atom, b: &Atom| -> Result<bool> {
                let difference = super::super::tests::canonical_tensor(
                    &crate::uv::approx::integrated::Integrated::dimensionally_regularized(&(a - b)),
                )?;
                Ok(difference
                    .series(GS.dim_epsilon, 0, -1)?
                    .to_atom()
                    .expand()
                    .together()
                    .is_zero())
            };
            if !difference(&poles, &literal_poles)? {
                failures.push(format!(
                    "{name}: bare poles {poles} against {literal_poles}"
                ));
            }
            if !difference(&counterterm, &literal_counterterm)? {
                failures.push(format!(
                    "{name}: counterterm {counterterm} against {literal_counterterm}"
                ));
            }
        }
        assert!(failures.is_empty(), "{failures:#?}");
        Ok(())
    }

    #[test]
    fn extracted_dressed_bubble_preserves_bare_poles_and_original_counterterm() -> Result<()> {
        check_serial(1)
    }

    #[test]
    fn two_disjoint_dressed_insertions_preserve_three_loop_poles() -> Result<()> {
        check_serial(2)
    }

    #[test]
    fn grouped_multisite_sum_preserves_correlations_and_epsilon_depth() -> Result<()> {
        test_initialise()?;
        let template = serial_bubbles(2)?;
        let epsilon = Atom::var(GS.dim_epsilon);
        let factors = [
            [Atom::one() + 2 * &epsilon, Atom::num(3) - &epsilon],
            [Atom::num(4) - &epsilon, Atom::num(2) + 3 * &epsilon],
        ];
        let weights = [Atom::num(2) + &epsilon, Atom::num(-3) + 2 * &epsilon];
        let mut settings = super::super::tests::settings();
        settings.vakint.normalization = "GScheme".into();
        let evaluator = PropagatorInsertionEvaluator::new(settings.clone());
        let vakint = Vakint::new()?;
        let mut prepared = Vec::new();
        let mut direct = Vec::new();
        for factors in &factors {
            let fragments = [VertexId(2), VertexId(4)]
                .into_iter()
                .zip(factors)
                .map(|(id, factor)| {
                    let vertex = template
                        .vertices()
                        .find(|(vertex, _)| *vertex == id)
                        .unwrap()
                        .1;
                    Ok((
                        id,
                        SymbolicTensor::infer(
                            factor * NumeratorInterface::export(&vertex.numerator).0,
                        )?,
                    ))
                })
                .collect::<Result<Vec<_>>>()?;
            let diagram = template.clone().with_numerator_fragments(fragments, [])?;
            let graph = diagram.to_gamma_loop_graph(None, true)?;
            direct.push(integrate_propagator(
                &graph,
                diagram.model(),
                &vakint,
                &settings,
                3,
            )?);
            let term = evaluator
                .prepare(&diagram, &vakint)?
                .expect("two proper insertions");
            assert_eq!(term.holes().len(), 2);
            assert_eq!(term.loop_count(), 3);
            // Both coefficients have been carried through the second replacement.
            for hole in term.holes() {
                let fragment = term
                    .surrogate()
                    .vertices()
                    .find(|(id, _)| *id == hole.vertex)
                    .unwrap()
                    .1;
                assert_eq!(
                    hole.coefficient,
                    NumeratorInterface::export(&fragment.numerator).0
                );
            }
            let individual = evaluator.pole_part_of_prepared(&term, &vakint)?;
            super::super::tests::assert_same_counterterm(
                &individual.poles,
                direct.last().unwrap(),
            )?;
            prepared.push(term);
        }
        let group = evaluator.weighted_pole_part(
            &[
                (&prepared[0], weights[0].clone()),
                (&prepared[1], weights[1].clone()),
            ],
            &vakint,
        )?;
        assert_eq!(
            (group.members, group.insertions, group.insertion_loops),
            (2, 4, 4)
        );
        let actual = prepared[0].restore_atom(&group.poles)?;
        let poles = |value: Atom| -> Result<Atom> {
            Ok(value
                .series(GS.dim_epsilon, 0, -1)?
                .to_atom()
                .expand()
                .together())
        };
        let expected = poles(&weights[0] * &direct[0] + &weights[1] * &direct[1])?;
        assert!(!expected.is_zero());
        super::super::tests::assert_same_counterterm(&actual, &expected)?;
        // Independent family sums would invent the absent (A1,B2)/(A2,B1)
        // tuples. Their nonzero discrepancy makes this a correlation test.
        let cartesian = poles(
            &direct[0]
                * (&weights[0] * &factors[0][0] + &weights[1] * &factors[1][0])
                * (&factors[0][1] + &factors[1][1])
                / (&factors[0][0] * &factors[0][1]),
        )?;
        assert!(!(cartesian - &expected).expand().together().is_zero());
        let constant_weights = poles(2 * &direct[0] - 3 * &direct[1])?;
        assert!(
            !(constant_weights - expected).expand().together().is_zero(),
            "positive epsilon weights must contribute to lower poles"
        );
        Ok(())
    }

    #[test]
    fn grouped_weights_reject_singular_dimension_and_open_tensors() -> Result<()> {
        test_initialise()?;
        PropagatorInsertionEvaluator::validate_weight(&(Atom::num(2) + Atom::var(GS.dim_epsilon)))?;
        for weight in [
            Atom::one() / Atom::var(GS.dim_epsilon),
            Atom::one() / (Atom::var(GS.dim) - 4),
            function!(GS.emr_mom, 0, Minkowski {}.new_rep(GS.dim).to_symbolic([])),
            ETS.metric(
                Minkowski {}
                    .new_rep(GS.dim)
                    .to_symbolic([Atom::from(Aind::Hedge(0, 0))]),
                Minkowski {}
                    .new_rep(GS.dim)
                    .to_symbolic([Atom::from(Aind::Hedge(1, 0))]),
            ),
        ] {
            assert!(PropagatorInsertionEvaluator::validate_weight(&weight).is_err());
        }
        Ok(())
    }

    /// MATRIX-R45, entry insertion_coefficient: an insertion group weight
    /// carrying ε⁻¹ is rejected in every spelling (explicit, `1/(d−4)`, a
    /// closed metric trace and a closed Dirac trace with fresh slots, which
    /// equal `d` and `4d`); their regular neighbours are admitted.
    #[test]
    fn matrix_r45_insertion_weights_reject_every_spelling_of_a_pole() -> Result<()> {
        test_initialise()?;
        let eps = "gammalooprs::{}::ε";
        let dim = "gammalooprs::{}::dim";
        let hedge = "gammalooprs::{spenso::index}::hedge";
        let mink = |n: usize| {
            format!("spenso::{{spenso::representation,spenso::self_dual}}::mink(4,{hedge}({n},1))")
        };
        let bis = |n: usize| {
            format!("spenso::{{spenso::representation,spenso::self_dual}}::bis(4,{hedge}({n},1))")
        };
        let metric = format!(
            "(spenso::{{symmetric,linear,real,spenso::tensor}}::g({},{}))^2",
            mink(81),
            mink(82)
        );
        let gamma = "spenso::{linear,spenso::tensor}::gamma";
        let dirac = format!(
            "{gamma}({},{},{})*{gamma}({},{},{})",
            bis(83),
            bis(84),
            mink(85),
            bis(84),
            bis(83),
            mink(85)
        );
        let parse = |text: &str| -> Result<Atom> {
            symbolica::try_parse!(text).map_err(|error| eyre!("{error}"))
        };
        for (label, pole, regular) in [
            ("eps", format!("({eps})^(-1)"), format!("1+{eps}")),
            ("dm4", format!("(-2)/({dim}-4)"), format!("(-2)/({dim}-3)")),
            (
                "gtr",
                format!("(-2)/({metric}-4)"),
                format!("(-2)/({metric}-3)"),
            ),
            (
                "dtr",
                format!("(-8)/({dirac}-16)"),
                format!("(-8)/({dirac}-12)"),
            ),
        ] {
            let pole = parse(&pole)?;
            let outcome = PropagatorInsertionEvaluator::validate_weight(&pole);
            eprintln!(
                "IRH_R45_WEIGHT {label} {:?}",
                outcome.as_ref().map_err(|e| format!("{e:#}"))
            );
            assert!(
                outcome.is_err(),
                "{label}: a pole in an insertion weight is accepted"
            );
            PropagatorInsertionEvaluator::validate_weight(&parse(&regular)?)
                .map_err(|error| eyre!("{label} regular neighbour: {error:#}"))?;
        }
        Ok(())
    }

    #[test]
    fn grouped_tensor_numerator_readmits_restored_dummy_indices() -> Result<()> {
        test_initialise()?;
        let minkowski = Minkowski {}.new_rep(GS.dim);
        let mu = minkowski.to_symbolic([Atom::from(Aind::Dummy(1_000_001))]);
        let nu = minkowski.to_symbolic([Atom::from(Aind::Dummy(1_000_002))]);
        let trace = ETS.metric(&mu, &nu).pow(2);
        let diagram = serial_bubbles(1)?.with_numerator_fragments(
            [(VertexId(0), NumeratorInterface::admit(trace, None)?)],
            [],
        )?;
        let mut settings = super::super::tests::settings();
        settings.vakint.normalization = "GScheme".into();
        let evaluator = PropagatorInsertionEvaluator::new(settings.clone());
        let vakint = Vakint::new()?;
        let prepared = evaluator.prepare(&diagram, &vakint)?.unwrap();
        let identity = prepared.skeleton_identity()?;
        let transported = identity.transport_from(identity, &prepared.complete_numerator())?;
        assert!(
            transported.contains_symbol(GS.dummyaind),
            "the fixture must exercise restored functional dummy indices"
        );
        let group = evaluator.weighted_pole_part(
            &[(&prepared, Atom::one()), (&prepared, Atom::num(2))],
            &vakint,
        )?;
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let direct = integrate_propagator(&graph, diagram.model(), &vakint, &settings, 2)?;
        assert!(!direct.is_zero());
        super::super::tests::assert_same_counterterm(
            &prepared.restore_atom(&group.poles)?,
            &(3 * direct),
        )
    }

    #[test]
    fn zero_grouped_bare_sum_does_not_replace_original_forests() -> Result<()> {
        test_initialise()?;
        let diagram = serial_bubbles(1)?;
        let mut settings = super::super::tests::settings();
        settings.vakint.normalization = "GScheme".into();
        let evaluator = PropagatorInsertionEvaluator::new(settings.clone());
        let vakint = Vakint::new()?;
        let prepared = evaluator.prepare(&diagram, &vakint)?.unwrap();
        let grouped = evaluator.weighted_pole_part(
            &[(&prepared, Atom::one()), (&prepared, Atom::num(-1))],
            &vakint,
        )?;
        assert!(grouped.poles.is_zero());
        assert_eq!(grouped.members, 2);
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let plan = RStarPlan::with_automatic_edge(&graph, diagram.model(), &settings)?;
        let forest = plan.counterterm_from_pole_part(
            &graph,
            diagram.model(),
            &vakint,
            &settings,
            &Atom::Zero,
        )?;
        assert!(
            !forest.is_zero(),
            "a zero bare group is not a zero-forest certificate"
        );
        let other_evaluator = PropagatorInsertionEvaluator::new(settings);
        assert!(
            other_evaluator
                .pole_part_of_prepared(&prepared, &vakint)
                .is_err()
        );
        Ok(())
    }

    #[test]
    fn isomorphic_bare_insertions_share_complete_amplitudes_at_fixed_depth() -> Result<()> {
        test_initialise()?;
        let diagram = serial_bubbles(2)?;
        let insertions = diagram
            .propagator_insertions()?
            .into_iter()
            .filter(|insertion| insertion.loops == 1)
            .map(|insertion| diagram.excise_with_mapping(&insertion.support, &[]))
            .collect::<Result<Vec<_>, _>>()?;
        assert_eq!(insertions.len(), 2);
        let mut settings = super::super::tests::settings();
        settings.vakint.normalization = "GScheme".into();
        let evaluator = PropagatorInsertionEvaluator::new(settings.clone());
        assert!(evaluator.cache_timing().is_none());
        let evaluator = evaluator.with_timing(true);
        let vakint = Vakint::new()?;
        for insertion in &insertions {
            let cached = evaluator.integrate_cached(&insertion.diagram, &vakint, 3)?;
            let graph = insertion.diagram.to_gamma_loop_graph(None, true)?;
            let direct = integrate_propagator(&graph, diagram.model(), &vakint, &settings, 3)?;
            super::super::tests::assert_same_counterterm(&cached, &direct)?;
            assert_eq!(evaluator.cache.lock().unwrap().len(), 1);
        }
        let deeper = evaluator.integrate_cached(&insertions[0].diagram, &vakint, 4)?;
        assert_eq!(evaluator.cache.lock().unwrap().len(), 2);
        let graph = insertions[0].diagram.to_gamma_loop_graph(None, true)?;
        let direct = integrate_propagator(&graph, diagram.model(), &vakint, &settings, 4)?;
        super::super::tests::assert_same_counterterm(&deeper, &direct)?;
        let timing = evaluator.cache_timing().unwrap();
        assert_eq!(
            (timing.calls, timing.initializers, timing.reuses),
            (3, 2, 1)
        );
        assert_eq!(timing.pending_reuses, 0);
        assert!(timing.key_ns > 0 && timing.restore_ns > 0);
        Ok(())
    }

    #[test]
    fn cached_bare_amplitude_restores_the_callers_single_external_momentum() -> Result<()> {
        test_initialise()?;
        let template = scalar_bubble(&serial_bubbles(1)?)?.with_numerator_fragments(
            [(
                VertexId(0),
                SymbolicTensor::infer(Atom::one() + Atom::var(GS.dim_epsilon))?,
            )],
            [],
        )?;
        let mut settings = super::super::tests::settings();
        settings.vakint.normalization = "GScheme".into();
        let evaluator = PropagatorInsertionEvaluator::new(settings.clone());
        let vakint = Vakint::new()?;
        let _ = evaluator.integrate_cached(&template, &vakint, 3)?;
        for reverse_outgoing in [false, true] {
            // Keep the decorated vertex fixed while making the other boundary
            // edge the first external momentum chosen by the direct integrator.
            // Equal-flow boundary legs additionally require Q_other = -Q_first.
            let mut builder =
                FeynmanDiagram::builder(template.model_arc(), "reordered_bare_bubble");
            for (_, vertex) in template.vertices() {
                builder.add_vertex(vertex.clone());
            }
            let edges = template.edges().collect::<Vec<_>>();
            for (_, endpoints, edge) in edges.into_iter().rev() {
                if reverse_outgoing && endpoints.target.is_none() {
                    let mut incoming = edge.clone();
                    incoming.external.as_mut().unwrap().state = ExternalState::Incoming;
                    builder.add_edge_with_slots(
                        endpoints.target,
                        endpoints.source,
                        incoming,
                        edge.target_slot(),
                        edge.source_slot(),
                    )?;
                } else {
                    builder.add_edge_with_slots(
                        endpoints.source,
                        endpoints.target,
                        edge.clone(),
                        edge.source_slot(),
                        edge.target_slot(),
                    )?;
                }
            }
            let reordered = builder.build()?;
            let actual = evaluator.integrate_cached(&reordered, &vakint, 3)?;
            assert_eq!(evaluator.cache.lock().unwrap().len(), 1);
            let graph = reordered.to_gamma_loop_graph(None, true)?;
            let expected = integrate_propagator(&graph, reordered.model(), &vakint, &settings, 3)?;
            let canonical = evaluator
                .cache
                .lock()
                .unwrap()
                .values()
                .next()
                .unwrap()
                .get()
                .unwrap()
                .as_ref()
                .unwrap()
                .clone();
            let support = graph
                .no_dummy()
                .subtract(&graph.external_filter::<SuBitGraph>());
            let identity = CountertermIdentity::new(&graph, &support, reordered.model())?;
            let before_routing = identity.restore_delta(&canonical)?;
            assert!(
                !(&before_routing - &expected).expand().together().is_zero(),
                "the fixture must expose the different external-momentum convention before routing"
            );
            super::super::tests::assert_same_counterterm(&actual, &expected)?;
        }
        Ok(())
    }
}
