//! The numerator engine of the helper's numerator assignment (lane engine,
//! Phase 2 ideas (a), (b), (c)): a [`FactorizedNumerator`] (one colour branch
//! of one graph, its numerator still a product of spenso factors) goes
//! straight to the folded channel inputs RustRed receives, without ever
//! building the contracted scalar numerator as an atom and without monomials
//! of canonical scalar products.
//!
//! 1. MATCH FIRST, on the denominators only ([`match_denominators`]): the
//!    graph's propagators become the `topo(prop(..))` that vakint's ordinary
//!    p path (R2, `VakintExpression::from_diagram`) builds; it is normalised
//!    and matched once (`p_lowering::match_cut`, numerator 1). The matched
//!    images of `k(i)` and `p(1)` turn every edge momentum of the input into a
//!    vector over the canonical `k(1..L)` and the head: one linear map applied
//!    to the factors.
//! 2. LOWER WHILE CONTRACTING ([`kernel`]): the running object is an open
//!    tensor whose coefficients are sparse polynomials in the presentation's
//!    SLOT variables with exact coefficients in `Q[D]`. A scalar product of
//!    two momenta becomes its slot form the moment it is formed (forms of the
//!    canonical products from RustRed's `express_scalar_product`, combined
//!    bilinearly; at a three-point vertex this is the propagator-only form).
//!    Metric renames (Feynman-gauge propagators) are applied first; the other
//!    factors are contracted one by one in a parameterised order
//!    ([`ContractionOrder`]).
//! 3. PRUNE EARLY: a monomial whose cancelled lines already make its sector
//!    zero (RustRed's census, lane prune's oracle `p_zero::ZeroSectors`) is
//!    dropped when formed, with everything it would generate.
//! 4. DIRAC: each closed fermion loop is a trace word handed to a
//!    [`trace::TraceEngine`] (idenso's production gamma algebra by default)
//!    and enters the contraction as one more factor.
//! 5. FOLD exactly as the reference does (`p_lowering::lower_closure_source`):
//!    key = line powers minus slot exponents, `scale_power = 2L − Σ powers +
//!    r` (`r` = loop-dependent canonical products), channel coefficient =
//!    spectator monomial × `dot(head, head)^j`, ε basis × value; the result
//!    is handed to the same reduction entry (`reduce_folded_channels`).
//!
//! Nothing existing changes behaviour: this is a new entry point.

mod kernel;
mod parse;
pub mod trace;

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use ahash::AHashMap as HashMap;
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::domains::rational::Rational;
use symbolica::function;

use self::kernel::{KernelCounters, Kinematics, Poly, SlotForm, Slots, Tensor, contract};
use self::parse::{Parser, split_fermion_factor, trace_word};
use self::trace::{TraceBackend, TraceEngine};
use super::p_data::PLibraryHandle;
use super::p_gate::{GateView, ZeroSectorCount};
use super::p_integral::{Channel, FoldedInputs, spectator_terms};
use crate::p_integral::factorized::{Endpoint, FactorizedNumerator, LineKind, names};
use crate::p_integral::{P_FAMILIES, PropagatorRequest};
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{RustRedMasterMode, Vakint, VakintError, VakintSettings, VakintTerm};

/// Why the engine could not lower a value.
#[derive(Debug, thiserror::Error)]
pub enum ContractError {
    #[error(transparent)]
    Vakint(#[from] VakintError),
    /// The value is outside what the engine handles (stated reason).
    #[error("p_contract: {0}")]
    Unsupported(String),
}

impl From<String> for ContractError {
    fn from(detail: String) -> Self {
        Self::Unsupported(detail)
    }
}

impl From<&str> for ContractError {
    fn from(detail: &str) -> Self {
        Self::Unsupported(detail.to_string())
    }
}

type Result<T> = std::result::Result<T, ContractError>;

/// The order in which the factors are multiplied into the running tensor.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum ContractionOrder {
    /// Next = the factor that leaves the fewest open indices (ties: most
    /// shared indices, then fewest terms, then the natural position).
    #[default]
    Greedy,
    /// The natural order: projector, vertices by id, edges by id, traces.
    Natural,
    /// Positions in the natural order (a permutation of all factors).
    Explicit(Vec<usize>),
}

/// Options of [`contract_and_lower`].
#[derive(Clone, Debug)]
pub struct ContractOptions {
    pub order: ContractionOrder,
    /// Early zero-sector pruning (needs the presentation's census table).
    pub prune: bool,
    pub trace: TraceBackend,
}

impl Default for ContractOptions {
    fn default() -> Self {
        Self {
            order: ContractionOrder::Greedy,
            prune: true,
            trace: TraceBackend::Idenso,
        }
    }
}

/// One contraction step (statistics).
#[derive(Clone, Debug, Default)]
pub struct StepStats {
    /// What was multiplied in (`vertex 3`, `edge 7`, `trace 0`, `projector`).
    pub factor: String,
    /// Terms (structure, monomial) of the factor and of the running tensor
    /// after the step, and its open indices.
    pub factor_terms: usize,
    pub terms: usize,
    pub open: usize,
    /// Candidate monomials formed and dropped (zero sector) in this step.
    pub formed: u64,
    pub dropped: u64,
}

/// Statistics of one lowering.
#[derive(Clone, Debug, Default)]
pub struct ContractStats {
    /// Seconds: match on denominators (+ momentum images), preparation
    /// (forms, zero table, renames, factor parsing), traces (the trace engine
    /// alone), contraction, fold.
    pub match_seconds: f64,
    pub prepare_seconds: f64,
    pub trace_seconds: f64,
    pub contract_seconds: f64,
    pub fold_seconds: f64,
    /// Factors multiplied, metric renames applied, traces and their lengths.
    pub factors: usize,
    pub renames: usize,
    pub trace_lengths: Vec<usize>,
    /// Largest running tensor (terms) and the final polynomial's terms.
    pub peak_terms: usize,
    pub final_terms: usize,
    /// Candidate monomials formed / dropped by the zero-sector test, in all.
    pub formed: u64,
    pub dropped: u64,
    /// Lowered (key, channel, basis) entries emitted and distinct keys.
    pub emitted: usize,
    pub keys: usize,
    /// Why there is no zero-sector oracle (empty when there is one).
    pub oracle: String,
    pub steps: Vec<StepStats>,
}

impl ContractStats {
    fn absorb(&mut self, other: &ContractStats) {
        self.match_seconds += other.match_seconds;
        self.prepare_seconds += other.prepare_seconds;
        self.trace_seconds += other.trace_seconds;
        self.contract_seconds += other.contract_seconds;
        self.fold_seconds += other.fold_seconds;
        self.factors += other.factors;
        self.renames += other.renames;
        self.trace_lengths.extend(&other.trace_lengths);
        self.peak_terms = self.peak_terms.max(other.peak_terms);
        self.final_terms += other.final_terms;
        self.formed += other.formed;
        self.dropped += other.dropped;
        self.emitted += other.emitted;
        if self.oracle.is_empty() {
            self.oracle.clone_from(&other.oracle);
        }
    }

    /// Seconds of contraction + lowering (everything but the match).
    pub fn lower_seconds(&self) -> f64 {
        self.prepare_seconds + self.trace_seconds + self.contract_seconds + self.fold_seconds
    }
}

/// The folded channel inputs of one or several values (branches of a graph):
/// what [`reduce_folded_channels`](super::p_integral::reduce_folded_channels)
/// receives.
#[derive(Clone, Debug, Default)]
pub struct ContractedInputs {
    pub(crate) folded: BTreeMap<Channel, FoldedInputs>,
    names: BTreeMap<(usize, Atom), BTreeSet<String>>,
    pub stats: ContractStats,
    /// Values absorbed (branches).
    pub values: usize,
    /// Values that were exact zeros (scaleless graph or zero numerator).
    pub zeros: usize,
}

impl ContractedInputs {
    /// Add another value's inputs (another colour branch of the same graph):
    /// the folded inputs are summed as `fold_requests` sums requests.
    pub fn absorb(&mut self, other: ContractedInputs) {
        self.stats.absorb(&other.stats);
        self.values += other.values;
        self.zeros += other.zeros;
        for (key, names) in other.names {
            self.names.entry(key).or_default().extend(names);
        }
        for (channel, inputs) in other.folded {
            let target = self.folded.entry(channel).or_default();
            target.integrals.extend(inputs.integrals);
            for (basis, keys) in inputs.by_basis {
                let sum = target.by_basis.entry(basis).or_default();
                for (key, value) in keys {
                    *sum.entry(key).or_insert_with(Rational::zero) += &value;
                }
            }
        }
        for inputs in self.folded.values_mut() {
            for keys in inputs.by_basis.values_mut() {
                keys.retain(|_, value| !value.is_zero());
            }
            inputs.by_basis.retain(|_, keys| !keys.is_empty());
        }
        self.folded.retain(|_, inputs| !inputs.by_basis.is_empty());
    }

    /// Folded channels.
    pub fn channels(&self) -> usize {
        self.folded.len()
    }

    /// (key, channel, basis) entries.
    pub fn entries(&self) -> usize {
        self.folded
            .values()
            .flat_map(|inputs| inputs.by_basis.values())
            .map(BTreeMap::len)
            .sum()
    }

    /// The gate view of the inputs (`p_gate::gate_view`), with the
    /// zero-sector keys dropped by the library's census.
    pub fn gate_view(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<(GateView, ZeroSectorCount)> {
        let handle = vakint.p_library().map_err(VakintError::from)?;
        let view = super::p_gate::gate_view(&self.folded, vk_symbol!(&settings.epsilon_symbol));
        Ok(super::p_gate::drop_zero_sectors(&handle, view))
    }

    /// Reduce the inputs through the same entry as the reference
    /// (`reduce_folded_channels`, the settings' fold mode, one RustRed
    /// reduction per folded channel) and finalise them (`finalize_reduced`
    /// with the master mode of the settings' RustRed method): the value
    /// vakint's evaluation returns for the same integrand. Also the seconds
    /// of the reduction (vakint `prepare` + library) and of the finalisation.
    pub fn reduce(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
    ) -> Result<(Atom, ReduceSeconds)> {
        use super::p_integral::{
            FoldMode, finalize_reduced, fold_mode, last_stats, reduce_folded_channels,
        };
        let handle = vakint.p_library().map_err(VakintError::from)?;
        let mode = settings
            .evaluation_order
            .0
            .iter()
            .find_map(|method| match method {
                crate::EvaluationMethod::RustRed(options) => Some(options.master_mode),
                _ => None,
            })
            .unwrap_or(RustRedMasterMode::Symbolic);
        let fold = match fold_mode() {
            FoldMode::Off => FoldMode::Coefficient,
            other => other,
        };
        let before = last_stats();
        let started = Instant::now();
        let reduced = reduce_folded_channels(
            vakint,
            settings,
            &handle,
            self.folded.clone(),
            None,
            fold,
            false,
        )?;
        let reduce_wall = started.elapsed().as_secs_f64();
        let after = last_stats();
        let finalizing = Instant::now();
        let value = finalize_reduced(vakint, settings, reduced, &self.names, mode)?;
        Ok((
            value,
            ReduceSeconds {
                reduce: reduce_wall,
                prepare: after.prepare_seconds - before.prepare_seconds,
                library: after.library_seconds - before.library_seconds,
                render: after.render_seconds - before.render_seconds,
                finalise: finalizing.elapsed().as_secs_f64(),
                reductions: after.reductions - before.reductions,
                pairs: after.pairs - before.pairs,
            },
        ))
    }
}

/// Seconds of [`ContractedInputs::reduce`].
#[derive(Clone, Debug, Default)]
pub struct ReduceSeconds {
    /// Wall time of `reduce_folded_channels` (prepare + library + render).
    pub reduce: f64,
    pub prepare: f64,
    pub library: f64,
    pub render: f64,
    pub finalise: f64,
    pub reductions: usize,
    pub pairs: usize,
}

/// `p(1)` of vakint's p ingress.
fn head_atom() -> Atom {
    function!(S.p, 1)
}

/// The topology R2 builds for `value` (`VakintExpression::from_diagram` of
/// a propagator graph): `topo(Π prop(id, edge(source, target), momentum, 0,
/// power))` over the internal edges sorted by id, momenta in `k(i)` and
/// `orientation · p(1)` with `p(1)` the momentum flowing INTO the graph
/// through the first external leg (smallest edge id). Also the orientation.
pub fn r2_topology(value: &FactorizedNumerator) -> Result<(Atom, i64)> {
    let loops = value.loops();
    if value.externals() != 1 {
        return Err(format!(
            "a two-point graph with one independent external momentum, not {}",
            value.externals()
        )
        .into());
    }
    let first = value
        .edges
        .iter()
        .filter(|edge| edge.is_external())
        .min_by_key(|edge| edge.id)
        .ok_or("a graph without external legs")?;
    if first.momentum[..loops].iter().any(|c| *c != 0) {
        return Err("an external leg carries a loop momentum".into());
    }
    let along = first.momentum[loops];
    let into = if first.source == Endpoint::External {
        along
    } else {
        -along
    };
    let orientation = match into {
        1 => 1,
        -1 => -1,
        _ => return Err("the first external leg does not carry ±p(1)".into()),
    };
    let head = head_atom();
    let mut topology = Atom::num(1);
    for edge in value.edges.iter().filter(|edge| !edge.is_external()) {
        let (Endpoint::Vertex(source), Endpoint::Vertex(target)) = (edge.source, edge.target)
        else {
            unreachable!("internal edges have two vertices");
        };
        if !edge.mass_squared.is_zero() {
            return Err("a massive line".into());
        }
        let mut momentum = Atom::Zero;
        for (i, c) in edge.momentum[..loops].iter().enumerate() {
            if *c != 0 {
                momentum += Atom::num(*c) * function!(S.k, i as i64 + 1);
            }
        }
        let c = edge.momentum[loops] * orientation;
        if c != 0 {
            momentum += Atom::num(c) * &head;
        }
        topology *= function!(
            S.prop,
            edge.id as i64,
            function!(S.edge, source as i64, target as i64),
            momentum,
            Atom::Zero,
            Atom::num(edge.power)
        );
    }
    Ok((function!(S.topo, topology), orientation))
}

/// The integer (rational) coefficients of a linear momentum over the
/// canonical `k(1..L)` and the head.
fn canonical_vector(momentum: &Atom, loops: usize, head: &Atom) -> Result<Vec<Rational>> {
    let expanded = momentum.expand();
    let mut vector = vec![Rational::zero(); loops + 1];
    let mut add = |term: AtomView| -> Result<()> {
        let (scale, vector_atom) = match term {
            AtomView::Mul(product) => {
                let factors = product.iter().collect::<Vec<_>>();
                match factors.as_slice() {
                    [n @ AtomView::Num(_), v] | [v, n @ AtomView::Num(_)] => (
                        parse::rational(*n).ok_or("a non-rational routing coefficient")?,
                        *v,
                    ),
                    _ => return Err("a nonlinear routing".into()),
                }
            }
            other => (Rational::one(), other),
        };
        let slot = if vector_atom == head.as_view() {
            loops
        } else {
            let AtomView::Fun(function) = vector_atom else {
                return Err(format!("a routing term {vector_atom} is not a momentum").into());
            };
            if function.get_symbol() != S.k || function.get_nargs() != 1 {
                return Err(format!("a routing term {vector_atom} is not a canonical k").into());
            }
            let index = i64::try_from(function.iter().next().unwrap())
                .map_err(|_| "a symbolic loop momentum index")?;
            if index < 1 || index as usize > loops {
                return Err("a loop momentum outside the canonical range".into());
            }
            index as usize - 1
        };
        vector[slot] += &scale;
        Ok(())
    };
    match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().try_for_each(&mut add)?,
        v if v.is_zero() => {}
        other => add(other)?,
    }
    Ok(vector)
}

/// The match of a graph on its denominators (idea (a)).
pub(crate) struct DenominatorMatch {
    pub(crate) family: usize,
    pub(crate) request_family: &'static crate::p_integral::PFamilyData,
    pub(crate) powers: Vec<i64>,
    pub(crate) external: Atom,
    pub(crate) head: Atom,
    pub(crate) sigma: i64,
    /// Canonical vector of every edge (internal and external) by edge id.
    pub(crate) edge_vectors: BTreeMap<usize, Vec<Rational>>,
}

/// Match `value`'s denominators exactly as R2 matches its graph (module
/// documentation, step 1). `Ok(None)`: the graph is scaleless (an exact zero).
pub(crate) fn match_denominators(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
) -> Result<Option<DenominatorMatch>> {
    let (topology, orientation) = r2_topology(value)?;
    // The steps of `p_lowering::match_cut` (normalise, match with numerator 1,
    // propagator request, head, scaleless split), without its per-call
    // `validate_settings` (a property of the settings, not of the graph: it
    // expands the G-scheme normalisation series every call; the caller's
    // evaluation validates them once).
    let mut term = VakintTerm {
        integral: topology,
        numerator: Atom::num(1),
        vectors: Vec::new(),
    };
    if !term.normalize_propagator()? {
        return Ok(None);
    }
    let factor = term.numerator.clone();
    if factor != Atom::num(1) {
        return Err(format!(
            "a normalisation factor {factor} of the graph (tree parts, negative powers)"
        )
        .into());
    }
    let _external = term
        .propagator_external()?
        .ok_or("the normalised graph is not a p-integral")?;
    let rules = term.prepare_integral_specs(vakint, settings, None)?;
    if !matches!(rules.canonical_topology, crate::Topology::Propagator(_)) {
        return Err("not a p-integral topology".into());
    }
    match vakint.evaluation_method(settings, &rules.canonical_topology) {
        Some(crate::EvaluationMethod::RustRed(_))
            if crate::p_integral::massive::rustred_max_rank().is_none() => {}
        _ => return Err("not reduced by RustRed alone".into()),
    }
    let one = Atom::num(1);
    let request = rules
        .propagator_request(one.as_view())
        .ok_or("a p-integral match without a propagator request")?;
    let (head, sigma) = crate::p_integral::external_head(&request.external)?;
    let scaleless = crate::p_integral::parts::split(
        request.family,
        request.powers,
        Atom::num(1).as_view(),
        &request.external,
    )?
    .is_none();
    let cut = super::p_lowering::MatchedCut {
        request_family: request.family,
        powers: request.powers.to_vec(),
        external: request.external.clone(),
        head,
        sigma,
        factor,
        rules: rules.clone(),
        scaleless,
    };
    if cut.scaleless {
        return Ok(None);
    }
    if cut.powers.iter().any(|power| *power < 0) {
        return Err("a negative line power".into());
    }
    let loops = cut.request_family.loops;
    if loops != value.loops() {
        return Err("the match changed the loop count".into());
    }
    // The canonical images of k(i) and p(1) through the matched rules.
    let image = |momentum: Atom| -> Result<Vec<Rational>> {
        let mut probe = VakintTerm {
            integral: Atom::num(1),
            numerator: momentum,
            vectors: Vec::new(),
        };
        probe.apply_numerator_replacement_rules(&cut.rules, settings)?;
        canonical_vector(&probe.numerator, loops, &cut.head)
    };
    let mut basis = Vec::with_capacity(loops + 1);
    for i in 1..=loops as i64 {
        basis.push(image(function!(S.k, i))?);
    }
    let mut p1 = image(head_atom())?;
    for entry in &mut p1 {
        *entry *= &Rational::from(orientation);
    }
    basis.push(p1);
    let mut edge_vectors = BTreeMap::new();
    for edge in &value.edges {
        let mut vector = vec![Rational::zero(); loops + 1];
        for (coefficient, image) in edge.momentum.iter().zip(&basis) {
            if *coefficient == 0 {
                continue;
            }
            let c = Rational::from(*coefficient);
            for (target, entry) in vector.iter_mut().zip(image) {
                *target += &(&c * entry);
            }
        }
        edge_vectors.insert(edge.id, vector);
    }
    let family = P_FAMILIES
        .iter()
        .position(|family| family.name == cut.request_family.name)
        .ok_or("an unknown family")?;
    Ok(Some(DenominatorMatch {
        family,
        request_family: cut.request_family,
        powers: cut.powers,
        external: cut.external,
        head: cut.head,
        sigma: cut.sigma,
        edge_vectors,
    }))
}

/// The slot forms of the canonical products `k_i·k_j` (i <= j, rows) then
/// `k_i·Q`, as the reference builds them (`p_lowering`, `canonical_rows`).
fn canonical_forms(
    presentation: &rustred::propagator::presentation::Presentation,
    loops: usize,
) -> Result<Vec<SlotForm>> {
    let unit = |position: usize| {
        let mut row = vec![0i64; loops + 1];
        row[position] = 1;
        row
    };
    let mut rows = Vec::new();
    for i in 0..loops {
        for j in i..loops {
            rows.push((unit(i), unit(j)));
        }
    }
    for i in 0..loops {
        rows.push((unit(i), unit(loops)));
    }
    rows.into_iter()
        .map(|(a, b)| {
            presentation
                .express_scalar_product(&a, &b)
                .map(|form| SlotForm {
                    constant: form.constant,
                    terms: form.terms,
                })
                .map_err(|_| "a scalar product the presentation cannot express".into())
        })
        .collect()
}

/// One factor of the contraction plan.
struct Factor {
    label: String,
    tensor: Tensor,
}

/// Apply the metric renames: every edge factor that is exactly one metric
/// `g(mink(D,a), mink(D,b))` (a Feynman-gauge propagator once its scalars are
/// pulled) identifies `a` with `b` in all other factors and disappears.
fn metric_renames(value: &FactorizedNumerator) -> (Vec<(usize, Atom)>, Vec<(String, Atom)>, usize) {
    let is_metric = |factor: &Atom| -> Option<(Atom, Atom)> {
        let AtomView::Fun(function) = factor.as_view() else {
            return None;
        };
        if function.get_symbol().get_name() != names::METRIC || function.get_nargs() != 2 {
            return None;
        }
        let args = function.iter().collect::<Vec<_>>();
        let (a, b) = (
            crate::p_integral::factorized::as_slot(args[0])?,
            crate::p_integral::factorized::as_slot(args[1])?,
        );
        (a.rep == names::MINK && b.rep == names::MINK).then(|| (a.atom, b.atom))
    };
    // Union-find over slot atoms.
    let mut parent = HashMap::<Atom, Atom>::default();
    fn find(parent: &HashMap<Atom, Atom>, slot: &Atom) -> Atom {
        let mut current = slot.clone();
        while let Some(next) = parent.get(&current) {
            current = next.clone();
        }
        current
    }
    let mut renames = 0;
    let mut kept_edges = Vec::new();
    for edge in &value.edges {
        if edge.factor == Atom::num(1) {
            continue;
        }
        if edge.line != LineKind::Fermion
            && let Some((a, b)) = is_metric(&edge.factor)
        {
            let (ra, rb) = (find(&parent, &a), find(&parent, &b));
            if ra != rb {
                parent.insert(ra, rb);
            }
            renames += 1;
            continue;
        }
        kept_edges.push((edge.id, edge.factor.clone()));
    }
    let apply = |atom: &Atom| -> Atom {
        if parent.is_empty() {
            return atom.clone();
        }
        atom.replace_map(|part, _, out| {
            if let Some(slot) = crate::p_integral::factorized::as_slot(part)
                && slot.rep == names::MINK
                && parent.contains_key(&slot.atom)
            {
                **out = find(&parent, &slot.atom);
            }
        })
    };
    let kept_edges = kept_edges
        .into_iter()
        .map(|(id, factor)| (id, apply(&factor)))
        .collect();
    let mut others = Vec::new();
    for vertex in &value.vertices {
        if vertex.factor != Atom::num(1) {
            others.push((format!("vertex {}", vertex.id), apply(&vertex.factor)));
        }
    }
    if let Some(projector) = &value.projector {
        others.push(("projector".to_string(), apply(&projector.tensor)));
    }
    (kept_edges, others, renames)
}

/// Whether an atom contains a spinor slot.
fn has_spinor(atom: &Atom) -> bool {
    let mut found = false;
    atom.visitor(&mut |part| {
        if found {
            return false;
        }
        if crate::p_integral::factorized::as_slot(part).is_some_and(|slot| slot.rep == names::BIS) {
            found = true;
            return false;
        }
        true
    });
    found
}

/// The half-edge ids of the spinor slots of an atom.
fn spinor_hedges(atom: &Atom) -> BTreeSet<usize> {
    let mut hedges = BTreeSet::new();
    atom.visitor(&mut |part| {
        if let Some(slot) = crate::p_integral::factorized::as_slot(part)
            && slot.rep == names::BIS
        {
            if let Some(h) = slot.hedge() {
                hedges.insert(h);
            }
            return false;
        }
        true
    });
    hedges
}

/// Lower one value (one colour branch of one graph) to its folded channel
/// inputs (module documentation). `prefactor` is the scalar that multiplies
/// the contracted value (in the value's dimension symbol; e.g.
/// `value.prefactor.total()` times a projection normalisation);
/// `dot(head, head)^p2_power` of the projector is multiplied in here.
pub fn contract_and_lower(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
    options: &ContractOptions,
) -> Result<ContractedInputs> {
    contract_and_lower_with(
        vakint,
        settings,
        value,
        prefactor,
        options,
        &mut GraphMemo::default(),
    )
}

/// What the values of one graph share (its colour branches have the same
/// denominators): the match on the denominators, keyed by the R2 topology.
#[derive(Default)]
pub struct GraphMemo {
    matches: HashMap<String, Option<std::sync::Arc<DenominatorMatch>>>,
    /// Per matched graph: the canonical kinematics (slot forms of the
    /// presentation, interned vectors, pair-form memo) and the vector of every
    /// edge id.
    setups: HashMap<String, (Kinematics, HashMap<i64, Option<u16>>)>,
    /// Matches served from the memo / computed.
    pub hits: usize,
    pub misses: usize,
}

impl GraphMemo {
    fn matched(
        &mut self,
        vakint: &Vakint,
        settings: &VakintSettings,
        value: &FactorizedNumerator,
    ) -> Result<(String, Option<std::sync::Arc<DenominatorMatch>>)> {
        let (topology, _) = r2_topology(value)?;
        let key = topology.to_canonical_string();
        if let Some(found) = self.matches.get(&key) {
            self.hits += 1;
            return Ok((key, found.clone()));
        }
        self.misses += 1;
        let found = match_denominators(vakint, settings, value)?.map(std::sync::Arc::new);
        self.matches.insert(key.clone(), found.clone());
        Ok((key, found))
    }
}

/// [`contract_and_lower`] with a memo shared by the values of one graph.
pub fn contract_and_lower_with(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
    options: &ContractOptions,
    memo: &mut GraphMemo,
) -> Result<ContractedInputs> {
    let mut out = ContractedInputs {
        values: 1,
        ..Default::default()
    };
    let mut stats = ContractStats::default();
    let started = Instant::now();
    let (key, matched) = memo.matched(vakint, settings, value)?;
    stats.match_seconds = started.elapsed().as_secs_f64();
    let Some(matched) = matched else {
        out.zeros = 1;
        out.stats = stats;
        return Ok(out);
    };
    if prefactor.is_zero() {
        out.zeros = 1;
        out.stats = stats;
        return Ok(out);
    }

    // Preparation: presentation, forms, zero table, renames, factors.
    let started = Instant::now();
    let handle: std::sync::Arc<PLibraryHandle> = vakint.p_library().map_err(VakintError::from)?;
    let top = handle.top(matched.family).map_err(VakintError::from)?;
    let presentation = handle.library.presentation(top.presentation);
    let loops = matched.request_family.loops;
    if usize::from(presentation.loops) != loops {
        return Err("the presentation's loop count differs".into());
    }
    let base = (0..top.arity)
        .map(|slot| {
            if slot < top.line_count {
                matched.powers[slot].max(0)
            } else {
                0
            }
        })
        .collect::<Vec<_>>();
    let zero = match handle.zero_sectors(top) {
        Ok(table) => Some(table),
        Err(reason) => {
            stats.oracle = reason.to_string();
            None
        }
    };
    let slots = Slots::new(top.arity, top.line_count, base, zero, options.prune)?;
    if !memo.setups.contains_key(&key) {
        let forms = canonical_forms(presentation, loops)?;
        let mut kinematics = Kinematics::new(loops, matched.sigma, forms);
        let mut edge_vectors = HashMap::default();
        for (id, vector) in &matched.edge_vectors {
            let interned = kinematics.intern(vector.clone())?;
            edge_vectors.insert(*id as i64, interned);
        }
        memo.setups.insert(key.clone(), (kinematics, edge_vectors));
    }
    let (kinematics, edge_vectors) = memo.setups.get_mut(&key).expect("inserted above");
    let kinematics = &mut *kinematics;
    let mut parser = Parser {
        dimension: value.dimension,
        edge_vectors: edge_vectors.clone(),
        indices: HashMap::default(),
    };
    let mut counters = KernelCounters::default();
    let (edges, others, renames) = metric_renames(value);
    stats.renames = renames;
    // Fermion loops: their edges' half-edges; every factor with spinor slots
    // goes to the loop whose half-edges it touches.
    let loops_of_graph = value
        .closed_loops()
        .map_err(|error| error.to_string())?
        .into_iter()
        .filter(|l| l.line == LineKind::Fermion)
        .collect::<Vec<_>>();
    let mut loop_of_hedge = BTreeMap::new();
    for (index, fermion_loop) in loops_of_graph.iter().enumerate() {
        for edge in &fermion_loop.edges {
            let edge = value.edge(*edge).ok_or("an unknown loop edge")?;
            for h in edge.source_hedge.into_iter().chain(edge.target_hedge) {
                loop_of_hedge.insert(h, index);
            }
        }
    }
    let mut fermion_factors = vec![Vec::new(); loops_of_graph.len()];
    let mut factors = Vec::<Factor>::new();
    let mut lorentz = Vec::<(String, Atom)>::new();
    for (label, atom) in edges
        .into_iter()
        .map(|(id, atom)| (format!("edge {id}"), atom))
        .chain(others)
    {
        if has_spinor(&atom) {
            let hedges = spinor_hedges(&atom);
            let owners = hedges
                .iter()
                .filter_map(|h| loop_of_hedge.get(h))
                .collect::<BTreeSet<_>>();
            let [owner] = owners.into_iter().collect::<Vec<_>>()[..] else {
                return Err(format!("{label}: spinor slots of several or no fermion loops").into());
            };
            fermion_factors[*owner].push((label, atom));
        } else {
            lorentz.push((label, atom));
        }
    }
    for (label, atom) in lorentz {
        let tensor = parser.parse(atom.as_view(), kinematics, &slots, &mut counters)?;
        factors.push(Factor { label, tensor });
    }
    stats.prepare_seconds = started.elapsed().as_secs_f64();

    // Traces.
    let mut engine: Box<dyn TraceEngine> = options.trace.engine();
    for (index, members) in fermion_factors.into_iter().enumerate() {
        let started = Instant::now();
        let split = members
            .iter()
            .map(|(label, atom)| {
                split_fermion_factor(atom.as_view()).map_err(|error| format!("{label}: {error}"))
            })
            .collect::<std::result::Result<Vec<_>, String>>()?;
        let (word, rest) = trace_word(split)?;
        stats.prepare_seconds += started.elapsed().as_secs_f64();
        let started = Instant::now();
        let traced = engine.trace(&word, value.dimension)?;
        stats.trace_seconds += started.elapsed().as_secs_f64();
        stats.trace_lengths.push(word.len());
        let started = Instant::now();
        let mut tensor = parser
            .parse(traced.as_view(), kinematics, &slots, &mut counters)
            .map_err(|error| format!("trace {index}: {error}"))?;
        for atom in rest {
            let part = parser.parse(atom.as_view(), kinematics, &slots, &mut counters)?;
            tensor = contract(&tensor, &part, kinematics, &slots, &mut counters)?;
        }
        stats.prepare_seconds += started.elapsed().as_secs_f64();
        factors.push(Factor {
            label: format!("trace {index} (length {})", word.len()),
            tensor,
        });
    }
    stats.factors = factors.len();

    // Contraction.
    let started = Instant::now();
    let order = plan(&factors, &options.order)?;
    let mut running = Tensor::one();
    let mut open = BTreeSet::<u32>::new();
    let mut slots_factors = factors.into_iter().map(Some).collect::<Vec<_>>();
    for position in order {
        let factor = slots_factors[position]
            .take()
            .ok_or("a factor used twice")?;
        let (formed, dropped) = (counters.formed, counters.dropped);
        running = contract(&running, &factor.tensor, kinematics, &slots, &mut counters)?;
        for index in factor.tensor.open_indices() {
            if !open.remove(&index) {
                open.insert(index);
            }
        }
        let size = running.size();
        stats.peak_terms = stats.peak_terms.max(size);
        stats.steps.push(StepStats {
            factor: factor.label,
            factor_terms: factor.tensor.size(),
            terms: size,
            open: open.len(),
            formed: counters.formed - formed,
            dropped: counters.dropped - dropped,
        });
        if running.terms.is_empty() {
            break;
        }
    }
    stats.contract_seconds = started.elapsed().as_secs_f64();
    stats.formed = counters.formed;
    stats.dropped = counters.dropped;
    if running.terms.keys().any(|structure| !structure.is_empty()) {
        return Err("the contraction left open indices (ports without a projector?)".into());
    }
    let polynomial = running.terms.remove(&Vec::new()).unwrap_or_default();
    stats.final_terms = polynomial.len();

    // Fold.
    let started = Instant::now();
    let projector_power = value.projector.as_ref().map_or(0, |p| p.p2_power);
    fold(
        settings,
        value,
        &matched,
        top.presentation,
        &slots,
        &polynomial,
        prefactor,
        projector_power,
        &mut out,
        &mut stats,
    )?;
    stats.fold_seconds = started.elapsed().as_secs_f64();
    stats.keys = out
        .folded
        .values()
        .flat_map(|inputs| inputs.by_basis.values())
        .flat_map(BTreeMap::keys)
        .collect::<BTreeSet<_>>()
        .len();
    out.stats = stats;
    Ok(out)
}

/// The order of the factors (positions in `factors`).
fn plan(factors: &[Factor], order: &ContractionOrder) -> Result<Vec<usize>> {
    let n = factors.len();
    match order {
        ContractionOrder::Natural => Ok((0..n).collect()),
        ContractionOrder::Explicit(positions) => {
            let mut sorted = positions.clone();
            sorted.sort_unstable();
            if sorted != (0..n).collect::<Vec<_>>() {
                return Err(format!("an explicit order must be a permutation of 0..{n}").into());
            }
            Ok(positions.clone())
        }
        ContractionOrder::Greedy => {
            let indices = factors
                .iter()
                .map(|factor| factor.tensor.open_indices())
                .collect::<Vec<_>>();
            let mut open = BTreeSet::<u32>::new();
            let mut left = (0..n).collect::<Vec<_>>();
            let mut out = Vec::with_capacity(n);
            while !left.is_empty() {
                let best = left
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, position)| {
                        let own = &indices[**position];
                        let shared = own.iter().filter(|i| open.contains(i)).count();
                        let after = open.len() + own.len() - 2 * shared;
                        (
                            after,
                            std::cmp::Reverse(shared),
                            factors[**position].tensor.size(),
                            **position,
                        )
                    })
                    .map(|(place, _)| place)
                    .unwrap();
                let position = left.remove(best);
                for index in &indices[position] {
                    if !open.remove(index) {
                        open.insert(*index);
                    }
                }
                out.push(position);
            }
            Ok(out)
        }
    }
}

/// Fold the final slot polynomial into the folded channels exactly as the
/// reference does (module documentation, step 5).
#[allow(clippy::too_many_arguments)]
fn fold(
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    matched: &DenominatorMatch,
    presentation: rustred::propagator::PresentationId,
    slots: &Slots,
    polynomial: &Poly,
    prefactor: &Atom,
    projector_power: i64,
    out: &mut ContractedInputs,
    stats: &mut ContractStats,
) -> Result<()> {
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    let in_epsilon = Atom::num(4) - Atom::num(2) * Atom::var(epsilon);
    let head_square = function!(S.dot, &matched.head, &matched.head);
    let total = prefactor * head_square.pow(Atom::num(projector_power));
    let total = total
        .replace(Atom::var(value.dimension).to_pattern())
        .with(in_epsilon.to_pattern());
    let pre_terms = spectator_terms(&total, epsilon);
    let loops = matched.request_family.loops;
    let total_power: i64 = matched.powers.iter().filter(|power| **power > 0).sum();
    let integral = super::p_integral::integral_name(&PropagatorRequest {
        family: matched.request_family,
        powers: &matched.powers,
        canonical_numerator: Atom::num(1).as_view(),
        external: matched.external.clone(),
    });
    let mut channels = HashMap::<(usize, u16), Atom>::default();
    let mut bases = HashMap::<(usize, usize), Atom>::default();
    let mut keys = HashMap::<u128, Vec<i64>>::default();
    for (mono, coefficient) in polynomial {
        let key = keys
            .entry(mono.slots)
            .or_insert_with(|| slots.key(mono.slots))
            .clone();
        let scale_power = 2 * loops as i64 - total_power + i64::from(mono.r);
        for (power, c) in coefficient.0.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            for (term, (monomial, basis, v)) in pre_terms.iter().enumerate() {
                let channel_coefficient = channels
                    .entry((term, mono.j))
                    .or_insert_with(|| {
                        if mono.j == 0 {
                            monomial.clone()
                        } else {
                            monomial * head_square.pow(Atom::num(i64::from(mono.j)))
                        }
                    })
                    .clone();
                let basis_atom = bases
                    .entry((term, power))
                    .or_insert_with(|| basis * in_epsilon.pow(Atom::num(power as i64)))
                    .clone();
                let target = out
                    .folded
                    .entry(Channel {
                        loops,
                        head: matched.head.clone(),
                        scale_power,
                        coefficient: channel_coefficient,
                    })
                    .or_default();
                target.integrals.insert(integral.clone());
                *target
                    .by_basis
                    .entry(basis_atom)
                    .or_default()
                    .entry((presentation, key.clone()))
                    .or_insert_with(Rational::zero) += &(v * c);
                stats.emitted += 1;
            }
        }
    }
    for inputs in out.folded.values_mut() {
        for keys in inputs.by_basis.values_mut() {
            keys.retain(|_, value| !value.is_zero());
        }
        inputs.by_basis.retain(|_, keys| !keys.is_empty());
    }
    out.folded.retain(|_, inputs| !inputs.by_basis.is_empty());
    let powers = matched
        .powers
        .iter()
        .map(i64::to_string)
        .collect::<Vec<_>>()
        .join(",");
    out.names
        .entry((loops, matched.head.clone()))
        .or_default()
        .insert(format!("{}({powers})", matched.request_family.name));
    Ok(())
}

#[cfg(test)]
mod tests;
