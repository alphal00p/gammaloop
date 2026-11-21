//! The R* CLOSURE kind of the uncontracted-numerator input mode: one R*
//! remainder (root remainder K1, UV child K2, IR hard part K3) as a
//! [`FactorizedNumerator`] of kind [`FactorizedKind::Closure`]: a vacuum graph
//! with ONE massive (closing) line among massless ones, external legs carrying
//! the independent external momenta `p(1..X)`, its numerator factorised and
//! uncontracted, and the R* Taylor step UNAPPLIED ([`TaylorSpec`]), or one
//! pre-expanded Taylor sector (`taylor = None`).
//!
//! What the record denotes (N/notes/wh/WH-RECORD-FORMAT.md section 4.1): with
//! `E = projector · PROD factors · PROD_e den_e^(-power_e)`, `den_e = q_e² −
//! m_e²`, every scaled `p(j) → λ p(j)` (in every edge momentum and every
//! `Q(x_j, ..)`), the reference mass `R → (1 − λ²) M_aux²` (compensate) or
//! `R → M_aux²` (fixed) everywhere, `E` expanded to total degree `w` in `λ`,
//! the orders of [`TaylorOrders`] kept, `λ = 1`, integrated over the loop
//! momenta; times the caller's prefactor and `dot(p(1), p(1))^p2_power`.
//!
//! How it is evaluated (the composition of the R* boundary path, by existing
//! code only; no tensor algebra here):
//! 1. idenso contracts `projector · PROD factors` ONCE per record
//!    ([`super::contract_with_idenso`], streamed when large), `R` replaced by
//!    `(1 − λ²) M_aux²` / `M_aux²` first (a scalar substitution; for `w = 0`
//!    the external momenta are removed from the factors first: `Q(x, ..) → 0`
//!    on the legs, `Q(e, ..) → Q(e, ..) − r_e` on the scaled lines, an exact
//!    root-grade-0 truncation); the closed scalar is mapped to vakint's
//!    notation (`dot` of the edges' momenta, open ports as raw indices);
//! 2. GRADED SPLIT: every summand is split into its q-structure `T` (the
//!    factors carrying an external momentum, an open port or `λ`) and a closed
//!    loop-dot scalar; the scalars of one `T` form one BLOCK `S_T`;
//! 3. the DENOMINATORS alone are Taylor expanded (scalar binomial series of
//!    `(A + B λ + C λ²)^(−n)`, `B = 2 l_e·r_e`, `C = r_e·r_e` (+ `M_aux²` on the
//!    compensated closing line)): every product term is a power PATTERN (the
//!    lines' raised powers, a key change) with a coefficient of known λ-degree;
//! 4. per pattern and block, the small coefficient `T × Σ_d d_{π,d}` (degrees
//!    kept by the Taylor degree and `orders`) is tensor-reduced by vakint's
//!    vacuum tensor reduction (memoised per q-monomial: the "q averaging" and
//!    the open-port projection onto external tensors), external tensors
//!    carried as scalar placeholders;
//! 5. the cut of the closing line is derived from the lines
//!    ([`super::vacuum::vacuum_backend`]); the blocks and patterns form ONE
//!    [`RemainderBoundary`] per record, evaluated by the existing boundary path
//!    (block reuse across patterns, one key plan per closing key, one
//!    reduction): the lowering efficiency of the R* boundary path is kept.
//!
//! Refusals: [`FactorizedError::Closure`] with the reason (see
//! `validate_closure`). Insertions (child counterterms, surrogate coefficient
//! vertices) are factors of the routed network; their summed pole order
//! raises the integral's ε depth.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::time::Instant;

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::function;

use super::vacuum::{
    self, VacuumBackend, VacuumLine, closure_inputs, contract_first_bound_bytes,
    distributed_lorentz_monomials, map_to_vakint_vacuum, stream_pieces_for, vacuum_backend,
    vacuum_lines, vacuum_topology,
};
use super::{
    ClosureSpec, FactorizedError, FactorizedKind, FactorizedNumerator, Result, TaylorSpec, names,
};
use crate::p_integral::{BoundaryOptions, BoundaryPattern, RemainderBoundary};
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{Vakint, VakintSettings};

/// Whether `settings` keep the masters as symbols (master mode Retain), the
/// mode the massive closure backend applies (the first RustRed or Forcer
/// method's master mode, Symbolic by default): its value keeps the exact ε
/// dependence, which no truncation may cut (format section 7).
fn retains_masters(settings: &VakintSettings) -> bool {
    settings
        .evaluation_order
        .0
        .iter()
        .find_map(|method| match method {
            crate::EvaluationMethod::RustRed(options) => Some(options.master_mode),
            crate::EvaluationMethod::Forcer(options) => Some(options.master_mode),
            _ => None,
        })
        == Some(crate::RustRedMasterMode::Retain)
}

/// Highest Taylor degree accepted (WH-RECORD-FORMAT OQ-5: QCD R* needs <= 2).
pub const MAX_CLOSURE_TAYLOR_DEGREE: u32 = 2;

fn refuse(detail: impl Into<String>) -> FactorizedError {
    FactorizedError::Closure {
        detail: detail.into(),
    }
}

/// An INTERNAL failure of the evaluation of an accepted record (lowering,
/// reduction, the ordinary closure, the truncation): `Evaluation`, never a
/// refusal variant, so that `Closure` means only "refused by the format".
fn internal(stage: &str, error: impl std::fmt::Display) -> FactorizedError {
    FactorizedError::Evaluation {
        detail: format!("closure {stage}: {error}"),
    }
}

#[cfg(test)]
thread_local! {
    /// TEST-ONLY fault injection: the internal-failure stage to fail (one of
    /// the four `internal(..)` call sites of the evaluation), per thread.
    pub(super) static INJECT_INTERNAL_FAILURE: std::cell::Cell<Option<&'static str>> =
        const { std::cell::Cell::new(None) };
}

/// TEST-ONLY fault injection INSIDE a wrapped call, a no-op in a library
/// build: when the tests asked for a failure at `stage`
/// (`INJECT_INTERNAL_FAILURE`), the wrapped operation's `result` is replaced
/// by its own NATIVE error `native()` (a `VakintError`, a `SeriesError`), so
/// that the call site's own `map_err` produces the reported error.
fn injected<T, E>(
    stage: &'static str,
    native: impl FnOnce() -> E,
    result: std::result::Result<T, E>,
) -> std::result::Result<T, E> {
    #[cfg(test)]
    if INJECT_INTERNAL_FAILURE.with(|c| c.get()) == Some(stage) {
        return Err(native());
    }
    let _ = (stage, native);
    result
}

/// The native error the injection gives a vakint call.
fn injected_vakint_error() -> crate::VakintError {
    crate::p_integral::invalid("injected failure")
}

fn lambda_symbol() -> Symbol {
    symbolica::symbol!("vakint::closure_lambda")
}

fn placeholder(n: usize) -> Atom {
    Atom::var(symbolica::symbol!(format!("vakint::closure_tensor_{n}")))
}

/// Whether `view` mentions a symbol whose name satisfies `test`.
fn mentions(view: AtomView<'_>, test: &dyn Fn(&str) -> bool) -> bool {
    let mut found = false;
    view.visitor(&mut |part| {
        match part {
            AtomView::Var(v) => found |= test(v.get_symbol().get_name()),
            AtomView::Fun(f) => found |= test(f.get_symbol().get_name()),
            _ => {}
        }
        !found
    });
    found
}

/// Whether `needle` occurs as a subexpression of `hay`.
pub(super) fn contains_atom(hay: AtomView<'_>, needle: AtomView<'_>) -> bool {
    let mut found = false;
    hay.visitor(&mut |part| {
        found |= part == needle;
        !found
    });
    found
}

/// `(fermion, ghost)` closed cycles of the carried graph; a line through
/// external legs (an open quark line of a hard part) is not a cycle.
fn carried_cycles(value: &FactorizedNumerator) -> (u32, u32) {
    let mut out = [0u32; 2];
    for (slot, kind) in [super::LineKind::Fermion, super::LineKind::Ghost]
        .into_iter()
        .enumerate()
    {
        let lines = value
            .edges
            .iter()
            .filter(|e| e.line == kind)
            .collect::<Vec<_>>();
        let mut visited = BTreeSet::new();
        for start in &lines {
            if visited.contains(&start.id) {
                continue;
            }
            let mut current = *start;
            visited.insert(current.id);
            loop {
                let next = match current.target {
                    super::Endpoint::Vertex(v) => lines
                        .iter()
                        .find(|e| e.source == super::Endpoint::Vertex(v)),
                    super::Endpoint::External => None,
                };
                match next {
                    Some(next) if next.id == start.id => {
                        out[slot] += 1;
                        break;
                    }
                    Some(next) if !visited.contains(&next.id) => {
                        visited.insert(next.id);
                        current = next;
                    }
                    _ => break,
                }
            }
        }
    }
    (out[0], out[1])
}

/// The record with its open slots ([`FactorizedNumerator::open_slots`]) as
/// its ports and no projector (the vakint notation's raw port indices `PORT +
/// i` and the caller notation follow them; the projector is contracted
/// already when this view is used).
fn opened(value: &FactorizedNumerator) -> Result<FactorizedNumerator> {
    let mut opened = value.clone();
    opened.ports = value.open_slots()?;
    opened.projector = None;
    Ok(opened)
}

/// The edge ids `Q(e, ..)` in `view`.
fn momentum_ids(view: AtomView<'_>) -> BTreeSet<usize> {
    let mut ids = BTreeSet::new();
    view.visitor(&mut |part| {
        if let AtomView::Fun(f) = part
            && f.get_symbol().get_name() == names::MOMENTUM
            && f.get_nargs() == 2
            && let Some(e) = f
                .iter()
                .next()
                .and_then(|a| i64::try_from(a).ok())
                .and_then(|a| usize::try_from(a).ok())
        {
            ids.insert(e);
        }
        true
    });
    ids
}

impl FactorizedNumerator {
    /// The closure rules (the record of N/notes/wh/WH-RECORD-FORMAT.md):
    /// CONTRACT-v1's graph, factor, index, locality and fermion rules
    /// (momentum conservation on all components for a conserved record with a
    /// Taylor step, on the loop components otherwise), then the closure kind's.
    pub(super) fn validate_closure(&self, spec: &ClosureSpec) -> Result<()> {
        let loops = self.loops();
        let externals = self.externals();
        let conserved = if spec.taylor.is_some() && !spec.independent_externals {
            loops + externals
        } else {
            loops
        };
        if loops == 0 {
            return Err(refuse(
                "a closure record has at least one loop (a zero-loop term, e.g. an inserted counterterm alone, stays with the caller)",
            ));
        }
        // MASSES ARE PLAIN SYMBOLS: the reference mass, the auxiliary mass
        // (unless zero: refused below) and every nonzero edge mass
        let plain = |a: &Atom| matches!(a.as_view(), AtomView::Var(_));
        if !plain(&spec.mass_squared) {
            return Err(refuse(format!(
                "the reference mass {} is not a symbol (it is substituted everywhere)",
                spec.mass_squared
            )));
        }
        if !spec.auxiliary_mass_squared.is_zero() && !plain(&spec.auxiliary_mass_squared) {
            return Err(refuse(format!(
                "the auxiliary mass {} is not a plain symbol",
                spec.auxiliary_mass_squared
            )));
        }
        if let Some(edge) = self
            .edges
            .iter()
            .find(|e| !e.mass_squared.is_zero() && !plain(&e.mass_squared))
        {
            return Err(refuse(format!(
                "edge {}: mass {} is neither 0 nor a plain symbol",
                edge.id, edge.mass_squared
            )));
        }
        // ids are written as Q(id, ..) (an i64), and a loop without a loop
        // edge gets the fresh id max + 1 + i ([`loop_vector_id`]): every id,
        // fresh ones included, must be an i64 (no wrap onto an existing id)
        // loop_vector_id(i) = max + 1 + i: bounded by the LARGEST loop index
        // without a loop edge, plus one
        let fresh = self
            .loop_edges
            .iter()
            .rposition(|e| e.is_none())
            .map_or(0, |i| i + 1);
        let top = self
            .edges
            .iter()
            .map(|e| e.id)
            .chain(self.constant_vectors())
            .max()
            .unwrap_or(0);
        if top
            .checked_add(fresh)
            .is_none_or(|last| last > i64::MAX as usize)
        {
            return Err(refuse(format!(
                "edge id {top}: the ids (with fresh loop-vector ids up to {top} + {fresh} above them) must not exceed i64::MAX"
            )));
        }
        self.validate_contract(Some(conserved))?;
        // weights, exactly once (WH-RECORD-FORMAT section 3): with
        // `original_cycles` the original graph's weights stay with the caller;
        // without, the carried graph's closed loops are applied here
        let p = &self.prefactor;
        // cycles inside insertions are the original graph's too
        let (inner_fermion, inner_ghost) = insertion_cycles(&self.insertions)
            .ok_or_else(|| refuse("the insertions' cycles overflow u32"))?;
        // a Taylor step raises a power by at most its degree w, with the
        // coefficient C(n + w - 1, w) w!: both must be i64 numbers
        if let Some(taylor) = &spec.taylor {
            let w = i64::from(taylor.degree);
            for edge in self.edges.iter().filter(|e| !e.is_external()) {
                if taylor_coefficient_bound(edge.power, w).is_none() {
                    return Err(refuse(format!(
                        "edge {}: power {} is too large for a Taylor step of degree {w} (its coefficients overflow i64)",
                        edge.id, edge.power
                    )));
                }
            }
        }
        match spec.original_cycles {
            Some((fermion, ghost)) => {
                if p.loop_sign != 1 || p.nf_power != 0 || p.overall != Atom::num(1) {
                    return Err(refuse(format!(
                        "with original_cycles the prefactor has loop_sign=1, nf_power=0, overall 1 (the original graph's weights stay with the caller); found loop_sign={} nf_power={} overall {}",
                        p.loop_sign, p.nf_power, p.overall
                    )));
                }
                let (f, g) = carried_cycles(self);
                let (Some(f), Some(g)) = (f.checked_add(inner_fermion), g.checked_add(inner_ghost))
                else {
                    return Err(refuse(
                        "the carried graph's and the insertions' cycles overflow u32",
                    ));
                };
                if f > fermion || g > ghost {
                    return Err(refuse(format!(
                        "the carried graph and its insertions have {f} fermion and {g} ghost cycles, more than original_cycles fermion={fermion} ghost={ghost}"
                    )));
                }
            }
            None => {
                let cycles = self.closed_loops()?;
                let f = cycles
                    .iter()
                    .filter(|c| c.line == super::LineKind::Fermion)
                    .count()
                    + inner_fermion as usize;
                let n = cycles.len() + inner_fermion as usize + inner_ghost as usize;
                let sign = if n.is_multiple_of(2) { 1 } else { -1 };
                if p.loop_sign != sign {
                    return Err(FactorizedError::LoopSign {
                        found: p.loop_sign,
                        loops: n,
                    });
                }
                if p.nf_power as usize != f {
                    return Err(FactorizedError::NfPower {
                        found: p.nf_power,
                        loops: f,
                    });
                }
            }
        }
        for (name, by_vakint, atom) in &spec.weights {
            if *by_vakint && *atom != p.scalars {
                return Err(refuse(format!(
                    "weight {name} applied_by=vakint is {atom}, prefactor_scalars is {}",
                    p.scalars
                )));
            }
        }
        let all = self
            .edges
            .iter()
            .map(|e| e.factor.clone())
            .chain(self.vertices.iter().map(|v| v.factor.clone()))
            .chain(self.projector.iter().map(|p| p.tensor.clone()))
            .chain(
                self.insertions
                    .iter()
                    .flat_map(|i| [i.coefficient.clone(), i.factor.clone()]),
            )
            .chain([p.scalars.clone(), p.colour.clone()])
            .collect::<Vec<_>>();
        for atom in &all {
            if mentions(atom.as_view(), &|n| n == names::NF || n.ends_with("::nf")) {
                return Err(refuse(
                    "a flavour symbol in the record (n_f stays with the caller)",
                ));
            }
            if mentions(atom.as_view(), &|n| n.ends_with("mUVexp")) {
                return Err(refuse("mUVexp in the record (the producer maps it to mUV)"));
            }
            if mentions(atom.as_view(), &|n| n == "vakint::closure_lambda") {
                return Err(refuse(
                    "vakint's internal Taylor parameter occurs in the record",
                ));
            }
        }
        // the closing line: the one massive internal edge
        let Some(closing) = self.edge(spec.edge) else {
            return Err(refuse(format!("closing edge {} does not exist", spec.edge)));
        };
        if closing.is_external() {
            return Err(refuse(format!(
                "closing edge {} is an external leg",
                spec.edge
            )));
        }
        if closing.mass_squared.is_zero() || closing.mass_squared != spec.mass_squared {
            return Err(refuse(format!(
                "closing edge {}: mass {} is zero or differs from the closure line's {}",
                spec.edge, closing.mass_squared, spec.mass_squared
            )));
        }
        if closing.power != spec.power || spec.power < 1 {
            return Err(refuse(format!(
                "closing edge {}: power {} differs from the closure line's {} (or < 1)",
                spec.edge, closing.power, spec.power
            )));
        }
        let massive = self
            .edges
            .iter()
            .filter(|e| !e.is_external() && !e.mass_squared.is_zero())
            .map(|e| e.id)
            .collect::<Vec<_>>();
        if massive.len() != 1 {
            return Err(refuse(format!(
                "several massive lines ({massive:?}): a closure record has exactly one, the closing edge"
            )));
        }
        if spec.auxiliary_mass_squared.is_zero() {
            return Err(refuse("the auxiliary mass is zero"));
        }
        if contains_atom(p.scalars.as_view(), spec.mass_squared.as_view())
            || contains_atom(p.colour.as_view(), spec.mass_squared.as_view())
        {
            return Err(refuse(format!(
                "the reference mass {} in the prefactor (it is substituted in the integrand only)",
                spec.mass_squared
            )));
        }
        // projector: no loop momentum
        if let Some(projector) = &self.projector {
            for e in momentum_ids(projector.tensor.as_view()) {
                if spec.taylor.is_some() {
                    return Err(refuse(format!(
                        "the projector contains Q({e}, ..): with a Taylor step it acts after the step and must be momentum-free"
                    )));
                }
                if self.edge(e).is_some_and(|edge| !edge.is_external()) {
                    return Err(refuse(format!(
                        "the projector contains the momentum of internal edge {e}"
                    )));
                }
            }
            if projector.p2_power != 0 && externals == 0 {
                return Err(refuse("p2_power without an external momentum"));
            }
        }
        // open slots (ports the projector does not close, the projector's own
        // free slots): Lorentz only
        for slot in self.open_slots()? {
            let lorentz = matches!(slot.as_view(), AtomView::Fun(f) if f.get_symbol().get_name() == names::MINK);
            if !lorentz {
                return Err(refuse(format!(
                    "open slot {slot} is not a Lorentz slot (open bispinor ports are not supported; close the line with the projector)"
                )));
            }
        }
        for soft in &spec.softs {
            if soft.external == 0
                || soft.external > externals
                || self.external_edges[soft.external - 1] != Some(soft.edge)
            {
                return Err(refuse(format!(
                    "soft line {} edge {}: not the leg of p({})",
                    soft.external, soft.edge, soft.external
                )));
            }
        }
        match &spec.taylor {
            None => {
                if spec.auxiliary_mass_squared != spec.mass_squared {
                    return Err(refuse(
                        "a pre-expanded sector (taylor none) has no auxiliary mass: it must equal the closing mass",
                    ));
                }
                if let Some(edge) = self
                    .edges
                    .iter()
                    .find(|e| !e.is_external() && e.momentum[loops..].iter().any(|c| *c != 0))
                {
                    return Err(refuse(format!(
                        "taylor none (lambda^0 denominators): internal edge {} carries an external momentum",
                        edge.id
                    )));
                }
            }
            Some(taylor) => {
                if taylor.degree > MAX_CLOSURE_TAYLOR_DEGREE {
                    return Err(refuse(format!(
                        "Taylor degree {} above {MAX_CLOSURE_TAYLOR_DEGREE}",
                        taylor.degree
                    )));
                }
                if taylor.orders != super::TaylorOrders::All {
                    return Err(refuse(format!(
                        "orders={} is not supported: only orders=all is evaluated (the top-mass drop modes have no verified record)",
                        taylor.orders.as_str()
                    )));
                }
                if !orders_keep_any(taylor) {
                    return Err(refuse(format!(
                        "orders={} keeps no order up to degree {}",
                        taylor.orders.as_str(),
                        taylor.degree
                    )));
                }
                let mut scaled = BTreeSet::new();
                let mut seen = BTreeSet::new();
                for id in &taylor.scaled_edges {
                    let Some(leg) = self.edge(*id).filter(|e| e.is_external()) else {
                        return Err(refuse(format!("scaled_edges: {id} is not an external leg")));
                    };
                    if !seen.insert(*id) {
                        return Err(refuse(format!("scaled_edges: {id} twice")));
                    }
                    for (j, c) in leg.momentum[loops..].iter().enumerate() {
                        if *c != 0 {
                            scaled.insert(j + 1);
                        }
                    }
                }
                let expected = (1..=externals).collect::<BTreeSet<_>>();
                if scaled != expected {
                    return Err(refuse(format!(
                        "scaled_edges={:?} scale p{:?}: every independent external momentum p(1..{externals}) must be scaled",
                        taylor.scaled_edges, scaled
                    )));
                }
                if let Some(j) = self.external_edges.iter().position(Option::is_none) {
                    return Err(refuse(format!(
                        "p({}) has no external leg: the Taylor step needs Q(x_j) for every scaled momentum",
                        j + 1
                    )));
                }
                if !self.constant_vectors().is_empty() {
                    return Err(refuse(format!(
                        "constant external vectors {:?} (edges outside the graph) with a Taylor step: which momenta scale is not defined",
                        self.constant_vectors()
                    )));
                }
                if !matches!(spec.mass_squared.as_view(), AtomView::Var(_)) {
                    return Err(refuse(format!(
                        "the reference mass {} is not a symbol (it is substituted everywhere)",
                        spec.mass_squared
                    )));
                }
            }
        }
        // the vacuum lines the evaluation integrates (external momenta set to
        // zero, the closing line at the auxiliary mass, series merged): one
        // closing line among massless ones, or scaleless; a Taylor step only
        // raises powers, so every pattern has the same back end and the
        // evaluation never refuses for the lines' shape
        let mut sector = self.clone();
        sector.kind = FactorizedKind::Direct;
        for edge in sector.edges.iter_mut().filter(|e| !e.is_external()) {
            for c in edge.momentum[loops..].iter_mut() {
                *c = 0;
            }
            if edge.id == spec.edge {
                edge.mass_squared = spec.auxiliary_mass_squared.clone();
            }
        }
        // the back end's errors here are FORMAT defects of the record (several
        // massive lines after merging, a closing line without loop momentum or
        // without a unit coefficient): a Closure refusal, not "failed inside"
        let backend =
            vacuum_backend(&vacuum_lines(&sector), loops).map_err(|error| match error {
                FactorizedError::Evaluation { detail } => refuse(format!(
                    "the closure's vacuum lines are not one closing line among massless ones: {}",
                    detail.trim_start_matches("vacuum: ")
                )),
                other => other,
            })?;
        if let VacuumBackend::Family = backend {
            return Err(refuse(
                "every line is massive: not a closure (after series merging the closing line is the only line)",
            ));
        }
        // MASSES BY POSITION: the reference mass only as the closure mass and
        // the closing edge's mass, the auxiliary mass only as the taylor
        // mass (the prefactor's and the insertions' own rules ran above)
        for (field, atom) in self.record_atoms() {
            if field.ends_with(" mass") || field.starts_with("insertion ") {
                continue;
            }
            let prefactor = field.starts_with("prefactor_");
            if !prefactor && contains_atom(atom.as_view(), spec.mass_squared.as_view()) {
                return Err(refuse(format!(
                    "the reference mass {} in {field}: it may occur only as the closure mass and as the closing edge's mass",
                    spec.mass_squared
                )));
            }
            if spec.auxiliary_mass_squared != spec.mass_squared
                && contains_atom(atom.as_view(), spec.auxiliary_mass_squared.as_view())
            {
                return Err(refuse(format!(
                    "the auxiliary mass {} in {field}: it may occur only as the taylor line's mass",
                    spec.auxiliary_mass_squared
                )));
            }
        }
        Ok(())
    }
}

/// The value of one record of [`evaluate_factorized_closure`].
#[derive(Clone, Debug)]
pub struct ClosureRecordValue {
    /// In the caller's notation, times the prefactor and `(p²)^p2_power`.
    pub value: Atom,
    /// The same in vakint notation.
    pub vakint_value: Atom,
    /// Power patterns (nonzero, non-scaleless) and blocks of the boundary.
    pub patterns: usize,
    pub blocks: usize,
    /// Distinct q-monomials tensor-reduced.
    pub reduced_monomials: usize,
    /// Seconds: idenso contraction; graded split, den expansion and the
    /// coefficients' tensor reduction; the boundary evaluation (lowering,
    /// reduction, finalisation).
    pub contract_seconds: f64,
    pub split_seconds: f64,
    pub evaluate_seconds: f64,
    /// The boundary path's own stage seconds (lowering vs reduction).
    pub lower_seconds: f64,
    pub reduce_seconds: f64,
    /// What the boundary hands to RustRed: key plans (one reduction each),
    /// lowered key entries, distinct canonical loop monomials.
    pub key_plans: usize,
    pub lowered_entries: usize,
    pub distinct_monomials: usize,
}

/// The value of [`evaluate_factorized_closure`].
#[derive(Clone, Debug)]
pub struct ClosureValue {
    /// The sum over the records, caller's notation.
    pub value: Atom,
    /// The sum over the records, vakint notation.
    pub vakint_value: Atom,
    pub records: Vec<ClosureRecordValue>,
}

/// The vakint vector `Σ c_i v(i)` of `coefficients` over `head`.
fn vector(head: Symbol, coefficients: &[i64]) -> Atom {
    let mut out = Atom::Zero;
    for (i, c) in coefficients.iter().enumerate() {
        if *c != 0 {
            out += Atom::num(*c) * function!(head, Atom::num(i as i64 + 1));
        }
    }
    out
}

/// Whether the order filter of `taylor` keeps at least one order `0..=degree`.
pub(super) fn orders_keep_any(taylor: &TaylorSpec) -> bool {
    (0..=taylor.degree).any(|n| taylor.orders.keeps(n))
}

/// The summed `(fermion, ghost)` cycles of the insertions, `None` on a u32
/// overflow.
pub(super) fn insertion_cycles(insertions: &[super::Insertion]) -> Option<(u32, u32)> {
    insertions.iter().try_fold((0u32, 0u32), |(f, g), i| {
        Some((f.checked_add(i.cycles.0)?, g.checked_add(i.cycles.1)?))
    })
}

/// The largest Taylor coefficient factor of a line of power `n` raised by at
/// most `w`: `n + w` and `C(n + w - 1, w) w!` as i64, `None` on an overflow
/// (`C(n + m - 1, m) m!` grows with m for n >= 1).
fn taylor_coefficient_bound(n: i64, w: i64) -> Option<i64> {
    n.checked_add(w)?;
    binomial(n.checked_add(w - 1)?, w)?.checked_mul(factorial(w))
}

/// `C(n, k)` exactly (in i128), `None` when it is not an i64.
fn binomial(n: i64, k: i64) -> Option<i64> {
    let mut out = 1i128;
    for i in 0..k {
        out = out.checked_mul(i128::from(n) - i128::from(i))? / (i128::from(i) + 1);
    }
    i64::try_from(out).ok()
}

/// `n!` (`k <= MAX_CLOSURE_TAYLOR_DEGREE`: never overflows).
fn factorial(n: i64) -> i64 {
    (1..=n).product::<i64>().max(1)
}

/// One term of the denominator expansion: raised powers by edge id, λ-degree,
/// coefficient (vakint notation).
#[derive(Clone)]
struct DenTerm {
    raise: BTreeMap<usize, i64>,
    degree: u32,
    coefficient: Atom,
}

/// Step 3 of the module documentation: the Taylor terms of the denominator
/// product up to degree `w`.
fn denominator_terms(
    value: &FactorizedNumerator,
    spec: &ClosureSpec,
    taylor: &TaylorSpec,
) -> Result<Vec<DenTerm>> {
    let loops = value.loops();
    let w = taylor.degree as i64;
    let mut terms = vec![DenTerm {
        raise: BTreeMap::new(),
        degree: 0,
        coefficient: Atom::num(1),
    }];
    for edge in value.edges.iter().filter(|e| !e.is_external()) {
        let l = vector(S.k, &edge.momentum[..loops]);
        let r = vector(S.p, &edge.momentum[loops..]);
        let closing = edge.id == spec.edge;
        let compensated = closing && taylor.compensate_mass;
        if r.is_zero() && !compensated {
            continue;
        }
        let b = if r.is_zero() {
            Atom::Zero
        } else {
            Atom::num(2) * S.dot.call((&l, &r))
        };
        let mut c = if r.is_zero() {
            Atom::Zero
        } else {
            S.dot.call((&r, &r))
        };
        if compensated {
            c += &spec.auxiliary_mass_squared;
        }
        let n = edge.power;
        // (A + B λ + C λ²)^(−n) = Σ (−1)^m C(n+m−1, m) m!/(a! b!) B^a C^b λ^(a+2b) / A^(n+m)
        let mut own = Vec::new();
        for b_power in 0..=w / 2 {
            for a_power in 0..=(w - 2 * b_power) {
                if (a_power > 0 && b.is_zero()) || (b_power > 0 && c.is_zero()) {
                    continue;
                }
                let m = a_power + b_power;
                let number = n
                    .checked_add(m - 1)
                    .and_then(|top| binomial(top, m))
                    .and_then(|c| c.checked_mul(factorial(m)))
                    .ok_or_else(|| {
                        internal(
                            "denominator series",
                            format!("edge {}: power {n} overflows i64 (validated)", edge.id),
                        )
                    })?
                    / (factorial(a_power) * factorial(b_power))
                    * if m % 2 == 0 { 1 } else { -1 };
                // B^0 = C^0 = 1 (also where B or C vanishes: no 0^0)
                let mut coefficient = Atom::num(number);
                if a_power > 0 {
                    coefficient *= b.pow(Atom::num(a_power));
                }
                if b_power > 0 {
                    coefficient *= c.pow(Atom::num(b_power));
                }
                own.push((m, (a_power + 2 * b_power) as u32, coefficient));
            }
        }
        let mut next = Vec::new();
        for term in &terms {
            for (m, degree, coefficient) in &own {
                if term.degree + degree > taylor.degree {
                    continue;
                }
                let mut raise = term.raise.clone();
                if *m != 0 {
                    raise.insert(edge.id, *m);
                }
                next.push(DenTerm {
                    raise,
                    degree: term.degree + degree,
                    coefficient: (&term.coefficient * coefficient).expand(),
                });
            }
        }
        terms = next;
    }
    Ok(terms)
}

/// The q-structure factors of a summand (step 2): anything carrying an
/// external momentum `p`, an open index (two-argument `k`, `g`) or `λ`.
fn is_structure(view: AtomView<'_>, lambda: Symbol) -> bool {
    let base = match view {
        AtomView::Pow(p) => p.get_base_exp().0,
        other => other,
    };
    match base {
        AtomView::Var(v) => v.get_symbol() == lambda,
        AtomView::Fun(f) => {
            let symbol = f.get_symbol();
            symbol == S.g
                || (symbol == S.p)
                || (symbol == S.k && f.get_nargs() == 2)
                || (symbol == S.dot && base.contains_symbol(S.p))
        }
        _ => false,
    }
}

/// Whether `view` contains a q-structure atom anywhere ([`is_structure`]).
fn contains_structure(view: AtomView<'_>, lambda: Symbol) -> bool {
    let mut found = false;
    view.visitor(&mut |part| {
        found |= is_structure(part, lambda);
        !found
    });
    found
}

/// `(structure, its λ-degree, block scalar)` of an expanded scalar.
fn graded_split(expanded: &Atom, lambda: Symbol) -> Vec<(Atom, u32, Atom)> {
    let summands = match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().map(|t| t.to_owned()).collect::<Vec<_>>(),
        AtomView::Num(_) if expanded.is_zero() => Vec::new(),
        single => vec![single.to_owned()],
    };
    // a summand whose non-structure factors still hide a structure inside a
    // sum is expanded on its own (the contracted scalar is a flat sum of
    // monomials in practice: no full expansion of the whole scalar)
    let mut flat = Vec::with_capacity(summands.len());
    for summand in summands {
        let hidden = match summand.as_view() {
            AtomView::Mul(product) => product
                .iter()
                .any(|f| !is_structure(f, lambda) && contains_structure(f, lambda)),
            AtomView::Add(_) => true,
            single => !is_structure(single, lambda) && contains_structure(single, lambda),
        };
        if hidden {
            let expanded = summand.expand();
            match expanded.as_view() {
                AtomView::Add(sum) => flat.extend(sum.iter().map(|t| t.to_owned())),
                AtomView::Num(_) if expanded.is_zero() => {}
                _ => flat.push(expanded),
            }
        } else {
            flat.push(summand);
        }
    }
    let mut groups: Vec<(Atom, Vec<Atom>)> = Vec::new();
    let mut index = HashMap::<Atom, usize>::new();
    for summand in &flat {
        let factors = match summand.as_view() {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            single => vec![single],
        };
        let (structure, scalar): (Vec<_>, Vec<_>) =
            factors.into_iter().partition(|f| is_structure(*f, lambda));
        let structure = Atom::mul_many(structure.iter().copied());
        let scalar = Atom::mul_many(scalar.iter().copied());
        match index.get(&structure) {
            Some(i) => groups[*i].1.push(scalar),
            None => {
                index.insert(structure.clone(), groups.len());
                groups.push((structure, vec![scalar]));
            }
        }
    }
    groups
        .into_iter()
        .map(|(structure, scalars)| {
            let degree = lambda_degree(structure.as_view(), lambda);
            (
                structure.replace(Atom::var(lambda)).with(Atom::num(1)),
                degree,
                Atom::add_many(scalars.iter()),
            )
        })
        .collect()
}

/// The λ-degree of a q-structure monomial: its number of external momenta
/// (the routing's lambda sat on each of them) plus its power of λ (the mass
/// compensation's λ² on M_aux²).
fn lambda_degree(view: AtomView<'_>, lambda: Symbol) -> u32 {
    // external momenta: p(j) (an open index or inside a dot) count one each
    let one = |view: AtomView<'_>| -> u32 {
        match view {
            AtomView::Fun(f) if f.get_symbol() == S.p => 1,
            AtomView::Fun(f) if f.get_symbol() == S.dot => f
                .iter()
                .filter(|a| matches!(a, AtomView::Fun(g) if g.get_symbol() == S.p))
                .count() as u32,
            _ => 0,
        }
    };
    let factor = |view: AtomView<'_>| -> u32 {
        match view {
            AtomView::Pow(p) => {
                let (base, exponent) = p.get_base_exp();
                one(base).saturating_mul(degree_exponent(exponent))
            }
            other => one(other),
        }
    };
    let momenta = match view {
        AtomView::Mul(m) => m.iter().map(factor).fold(0, u32::saturating_add),
        other => factor(other),
    };
    // plus the lambda^2 of a compensated reference mass inside the factors
    momenta.saturating_add(lambda_power(view, lambda))
}

/// The per-record preparation: the boundary (or `None`: exactly zero) and
/// its statistics.
struct Prepared {
    boundary: Option<RemainderBoundary>,
    tensors: Vec<(Atom, Atom)>,
    blocks: usize,
    reduced_monomials: usize,
    contract_seconds: f64,
    split_seconds: f64,
}

/// The edge id whose `Q(id, ..)` stands for the loop momentum `k(i+1)` in the
/// routed factors: the record's loop edge `b_i` (its momentum is exactly
/// `k(i+1)`, as the producer routes onto its basis edges), else a fresh id
/// above every edge id of the record.
fn loop_vector_id(value: &FactorizedNumerator, i: usize) -> usize {
    match value.loop_edges.get(i).copied().flatten() {
        Some(edge) => edge,
        // above every edge id AND every constant external vector Q(x, ..) of
        // an absent edge the record names (a pre-expanded sector may carry
        // one): a fresh id no caller atom can contain
        // (validate_closure bounds max + fresh count by i64::MAX: no wrap)
        None => value
            .edges
            .iter()
            .map(|e| e.id)
            .chain(value.constant_vectors())
            .max()
            .unwrap_or(0)
            .saturating_add(1 + i),
    }
}

/// The ROUTING of a factor before the contraction (a linear substitution of
/// its momenta, no tensor algebra): every `Q(e, s)` of an internal edge becomes
/// `SUM_i eta_ei Q(b_i, s) + lam SUM_j xi_ej Q(x_j, s)` (x_j = the
/// leg of p(j)), every `Q(x, s)` of an external leg `lam Q(x, s)`; `lam` is
/// the Taylor parameter (lambda), 0 for degree 0, 1 for a pre-expanded sector.
/// Then, for a Taylor step of degree `w`, the factor's terms of lambda-degree
/// above `w` are dropped (they cannot survive the expansion: every other
/// factor is regular at lambda = 0).
fn routed_factor(
    value: &FactorizedNumerator,
    factor: &Atom,
    lam: &Atom,
    lambda: Symbol,
    degree: Option<u32>,
) -> Atom {
    let loops = value.loops();
    let routed = factor.replace_map(|part, _, out| {
        let AtomView::Fun(f) = part else { return };
        if f.get_symbol().get_name() != names::MOMENTUM || f.get_nargs() != 2 {
            return;
        }
        let mut args = f.iter();
        let (Some(e), Some(slot)) = (args.next(), args.next()) else {
            return;
        };
        let Some(edge) = i64::try_from(e)
            .ok()
            .and_then(|e| usize::try_from(e).ok())
            .and_then(|e| value.edge(e))
        else {
            return;
        };
        let q = |id: usize| function!(f.get_symbol(), Atom::num(id as i64), slot.to_owned());
        if edge.is_external() {
            // a leg through the legs of the independent momenta it carries
            // (the dependent leg of a two-point root becomes -Q(x0): one
            // external vector per p(j), as the R* sigma does)
            let mut image = Atom::Zero;
            let mut own = true;
            for (j, c) in edge.momentum[loops..].iter().enumerate() {
                if *c != 0 {
                    match value.external_edges[j] {
                        Some(leg) if leg != edge.id || *c != 1 => {
                            own = false;
                            image += Atom::num(*c) * q(leg);
                        }
                        Some(_) => image += q(edge.id),
                        None => own = true,
                    }
                }
            }
            **out = if own {
                lam * part.to_owned()
            } else {
                lam * image
            };
            return;
        }
        let mut image = Atom::Zero;
        for (i, c) in edge.momentum[..loops].iter().enumerate() {
            if *c != 0 {
                image += Atom::num(*c) * q(loop_vector_id(value, i));
            }
        }
        for (j, c) in edge.momentum[loops..].iter().enumerate() {
            if *c != 0 {
                let leg = value.external_edges[j].expect("validated: every p(j) has a leg");
                image += Atom::num(*c) * lam * q(leg);
            }
        }
        **out = image;
    });
    // flatten the substituted sums (`c (sum)` inside a sum), as the R* producer's
    // linearisation does before its contraction: the factor keeps its products
    let routed = linearized(routed);
    match degree {
        // truncate only a factor whose lambda-degree can exceed w (a bound read
        // off the tree without expanding): the usual vertex and propagator
        // factors are linear in the momenta and keep their factored form
        Some(w) if w > 0 && lambda_degree_bound(routed.as_view(), lambda) > w => {
            let expanded = routed.expand();
            let summands = match expanded.as_view() {
                AtomView::Add(sum) => sum.iter().map(|t| t.to_owned()).collect::<Vec<_>>(),
                _ => vec![expanded.clone()],
            };
            Atom::add_many(
                summands
                    .iter()
                    .filter(|t| lambda_power(t.as_view(), lambda) <= w),
            )
        }
        _ => routed,
    }
}

/// An upper bound of the lambda-degree of an expression, without expanding:
/// products add, sums take the maximum, `x^n` multiplies by n.
fn lambda_degree_bound(view: AtomView<'_>, lambda: Symbol) -> u32 {
    match view {
        AtomView::Var(v) => u32::from(v.get_symbol() == lambda),
        AtomView::Num(_) => 0,
        AtomView::Add(a) => a
            .iter()
            .map(|t| lambda_degree_bound(t, lambda))
            .max()
            .unwrap_or(0),
        AtomView::Mul(m) => m
            .iter()
            .map(|t| lambda_degree_bound(t, lambda))
            .fold(0, u32::saturating_add),
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            lambda_degree_bound(base, lambda).saturating_mul(degree_exponent(exponent))
        }
        AtomView::Fun(f) => f
            .iter()
            .map(|t| lambda_degree_bound(t, lambda))
            .fold(0, u32::saturating_add),
    }
}

/// `c (sum)` with a number `c`.
fn is_scaled_sum(view: AtomView<'_>) -> bool {
    match view {
        AtomView::Mul(product) => {
            product.get_nargs() == 2
                && product.iter().any(|f| matches!(f, AtomView::Num(_)))
                && product.iter().any(|f| matches!(f, AtomView::Add(_)))
        }
        _ => false,
    }
}

/// A sum with a nested (scaled) sum among its terms.
fn is_nested_sum(view: AtomView<'_>) -> bool {
    matches!(view, AtomView::Add(sum) if sum
        .iter()
        .any(|t| is_scaled_sum(t) || matches!(t, AtomView::Add(_))))
}

/// The expression with every nested sum flattened by distributing its numbers
/// (`q5 - q8 + q11 - (-q9 - q11)` -> `q5 + 2 q11 - q8 + q9`); products of
/// sums stay products. The same normalisation as the R* producer's
/// momentum linearisation; the value is unchanged.
fn linearized(atom: Atom) -> Atom {
    let mut nested = false;
    atom.as_view().visitor(&mut |part| {
        nested |= is_nested_sum(part);
        !nested
    });
    if !nested {
        return atom;
    }
    atom.replace_map(|term, _, out| {
        if is_nested_sum(term) {
            **out = term.expand_num();
        }
    })
}

/// The power of `lambda` in a product.
fn lambda_power(view: AtomView<'_>, lambda: Symbol) -> u32 {
    let one = |view: AtomView<'_>| -> u32 {
        match view {
            AtomView::Var(v) if v.get_symbol() == lambda => 1,
            AtomView::Pow(p) => match p.get_base_exp() {
                (AtomView::Var(v), exp) if v.get_symbol() == lambda => degree_exponent(exp),
                _ => 0,
            },
            _ => 0,
        }
    };
    match view {
        AtomView::Mul(m) => m.iter().map(one).fold(0, u32::saturating_add),
        other => one(other),
    }
}

/// The exponent of a power as a DEGREE multiplier: a nonnegative integer,
/// saturated at `u32::MAX` (a degree is only compared with the Taylor degree
/// `w <= 2`: a saturated degree is dropped as the true one would be, where a
/// truncation to the low 32 bits read `x^(2^32)` as order 0); a negative or
/// non-integer exponent counts 0 (as before).
fn degree_exponent(exponent: AtomView<'_>) -> u32 {
    match i64::try_from(exponent) {
        Ok(n) if n >= 0 => u32::try_from(n).unwrap_or(u32::MAX),
        Ok(_) => 0,
        Err(_) => match exponent {
            AtomView::Num(n)
                if matches!(
                    n.get_coeff_view(),
                    symbolica::coefficient::CoefficientView::Large(..)
                ) && exponent.to_owned().is_positive().is_true() =>
            {
                u32::MAX
            }
            _ => 0,
        },
    }
}

/// The record for the mapping to vakint notation after a routed contraction:
/// its open slots as ports and the loop-momentum vectors of loops without a loop edge as
/// edges carrying `k(i+1)` (external, so that they are never lines).
fn mapping_view(value: &FactorizedNumerator) -> Result<FactorizedNumerator> {
    let mut view = opened(value)?;
    let n = value.loops() + value.externals();
    for i in (0..value.loops()).filter(|i| value.loop_edges[*i].is_none()) {
        let mut momentum = vec![0; n];
        momentum[i] = 1;
        view.edges.push(super::FactorizedEdge {
            id: loop_vector_id(value, i),
            source: super::Endpoint::External,
            target: super::Endpoint::External,
            source_hedge: None,
            target_hedge: None,
            particle: "loop".into(),
            line: super::LineKind::Scalar,
            power: 0,
            mass_squared: Atom::Zero,
            momentum,
            factor: Atom::num(1),
        });
    }
    Ok(view)
}

fn prepare(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    spec: &ClosureSpec,
    prefactor: &Atom,
) -> Result<Prepared> {
    let _step = crate::timing::contract();
    let lambda = lambda_symbol();
    let loops = value.loops();
    let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
    let in_epsilon = Atom::num(4) - Atom::num(2) * &epsilon;
    let d = Atom::var(value.dimension);
    // the caller's prefactor rides in every pattern coefficient (inside the
    // Laurent expansion, as on the boundary path); (p²)^p2_power with the
    // external tensors
    let prefactor = prefactor.replace(d.clone()).with(in_epsilon.clone());
    let p2 = match &value.projector {
        Some(projector) if projector.p2_power != 0 => S
            .dot
            .call((S.p.call(1), S.p.call(1)))
            .pow(Atom::num(projector.p2_power)),
        _ => Atom::num(1),
    };
    // step 1: the network, R substituted, w = 0 truncated; idenso
    let started = Instant::now();
    let degree_zero = spec.taylor.as_ref().is_some_and(|t| t.degree == 0);
    let lam = match &spec.taylor {
        Some(_) if degree_zero => Atom::Zero,
        Some(_) => Atom::var(lambda),
        None => Atom::num(1),
    };
    let degree = spec.taylor.as_ref().map(|t| t.degree);
    let mut network = Atom::num(1);
    for factor in value
        .edges
        .iter()
        .map(|e| &e.factor)
        .chain(value.vertices.iter().map(|v| &v.factor))
        .chain(value.projector.iter().map(|p| &p.tensor))
    {
        network *= routed_factor(value, factor, &lam, lambda, degree);
    }
    // insertions (validated by `check_insertions`): coefficient x factor, a
    // factor like any other: routed and lambda-scaled (its boundary momenta),
    // truncated, R-substituted, contracted with the rest
    for insertion in &value.insertions {
        let product = &insertion.coefficient * &insertion.factor;
        network *= routed_factor(value, &product, &lam, lambda, degree);
    }
    // lambda has done its work (the per-factor truncation): every lambda
    // accompanies exactly one external leg vector, so the grade of a term is
    // read after the contraction as its number of external momenta; idenso
    // gets the producer's product without a symbolic coefficient on the
    // external momenta (with Dirac traces that coefficient costs ~25 %)
    network = network.replace(Atom::var(lambda)).with(Atom::num(1));
    if let Some(taylor) = &spec.taylor {
        let image = if taylor.compensate_mass && !degree_zero {
            (Atom::num(1) - Atom::var(lambda).pow(Atom::num(2))) * &spec.auxiliary_mass_squared
        } else {
            spec.auxiliary_mass_squared.clone()
        };
        network = network.replace(spec.mass_squared.clone()).with(image);
    }
    let routing_seconds = started.elapsed().as_secs_f64();
    if let Ok(path) = std::env::var("VAKINT_CLOSURE_DUMP_NETWORK") {
        let _ = std::fs::write(path, network.to_canonical_string());
    }
    // The routed product is contracted WHOLE, as the R* boundary path does
    // (one idenso call); the streamed contraction (pieces contracted one after
    // another) is only the memory fallback when the caller sets a process
    // budget that the contract-first bound of the size predictor exceeds.
    let predicted = distributed_lorentz_monomials(network.as_view());
    let pieces = match settings.rustred_resources.max_process_rss_bytes {
        Some(limit) if contract_first_bound_bytes(predicted) > limit as f64 => {
            stream_pieces_for(predicted)
        }
        _ => 1,
    };
    let contracted = match pieces {
        1 => super::oracle::contract_with_idenso_unexpanded(&network, false)?,
        n => vacuum::contract_streamed(settings, &network, n)?.0,
    };
    let idenso_seconds = started.elapsed().as_secs_f64();
    // the dimension is substituted per block after the split (no expansion of
    // the D-dependent coefficients into epsilon terms)
    let scalar = map_to_vakint_vacuum(&mapping_view(value)?, &contracted)?;
    let contract_seconds = started.elapsed().as_secs_f64();
    if std::env::var("VAKINT_CLOSURE_PROFILE").is_ok() {
        eprintln!(
            "CLOSURE_PROFILE label=\"{}\" loop_terms={} root_terms={} predicted_monomials={predicted:.3e} pieces={} network_bytes={} routing_s={routing_seconds:.4} idenso_s={idenso_seconds:.4} map_expand_s={:.4} contracted_terms={} contracted_bytes={} scalar_terms={}",
            value.label,
            value
                .edges
                .iter()
                .filter(|e| !e.is_external())
                .map(|e| e.momentum[..value.loops()]
                    .iter()
                    .filter(|c| **c != 0)
                    .count())
                .sum::<usize>(),
            value
                .edges
                .iter()
                .filter(|e| !e.is_external() && e.momentum[value.loops()..].iter().any(|c| *c != 0))
                .count(),
            pieces,
            network.as_view().get_byte_size(),
            contract_seconds - idenso_seconds,
            crate::census::summands(contracted.as_view()),
            contracted.as_view().get_byte_size(),
            crate::census::summands(scalar.as_view())
        );
    }

    // step 2: blocks
    let started = Instant::now();
    let split = graded_split(&scalar, lambda);
    if split.is_empty() {
        return Ok(Prepared {
            boundary: None,
            tensors: Vec::new(),
            blocks: 0,
            reduced_monomials: 0,
            contract_seconds,
            split_seconds: started.elapsed().as_secs_f64(),
        });
    }
    // step 3: the denominator terms (or the sector as it is)
    let (den, w, orders) = match &spec.taylor {
        Some(taylor) => (
            denominator_terms(value, spec, taylor)?,
            taylor.degree,
            taylor.orders,
        ),
        None => (
            vec![DenTerm {
                raise: BTreeMap::new(),
                degree: 0,
                coefficient: Atom::num(1),
            }],
            // a pre-expanded sector: every grade is kept as written
            u32::MAX,
            super::TaylorOrders::All,
        ),
    };
    // patterns: the raised powers; per pattern Σ over den terms by degree
    let mut patterns = BTreeMap::<BTreeMap<usize, i64>, Vec<&DenTerm>>::new();
    for term in &den {
        patterns.entry(term.raise.clone()).or_default().push(term);
    }
    // the lines of every pattern, its cut
    let base_powers = value
        .edges
        .iter()
        .filter(|e| !e.is_external())
        .map(|e| (e.id, e.power))
        .collect::<BTreeMap<_, _>>();
    let closing_mass = spec.auxiliary_mass_squared.clone();
    let head = S.p.call(value.externals() as i64 + 1);
    let lines_of = |raise: &BTreeMap<usize, i64>| -> Result<Vec<VacuumLine>> {
        let mut sector = value.clone();
        sector.kind = FactorizedKind::Direct;
        for edge in sector.edges.iter_mut().filter(|e| !e.is_external()) {
            for c in edge.momentum[loops..].iter_mut() {
                *c = 0;
            }
            edge.power = base_powers[&edge.id]
                .checked_add(raise.get(&edge.id).copied().unwrap_or(0))
                .ok_or_else(|| {
                    internal(
                        "raised powers",
                        format!(
                            "edge {}: power {} overflows i64 (validated)",
                            edge.id, base_powers[&edge.id]
                        ),
                    )
                })?;
            if edge.id == spec.edge {
                edge.mass_squared = closing_mass.clone();
            }
        }
        Ok(vacuum_lines(&sector))
    };
    // step 4: per pattern and block the reduced coefficient
    let base_lines = lines_of(&BTreeMap::new())?;
    let base_topology = S.topo.call(vacuum_topology(&base_lines));
    let mut memo = HashMap::<Atom, Vec<(Atom, Atom)>>::new();
    let mut tensor_index = BTreeMap::<Atom, usize>::new();
    let mut tensors = Vec::<(Atom, Atom)>::new();
    let base_cut = match vacuum_backend(&base_lines, loops)? {
        VacuumBackend::Closure(cut) => Some(cut),
        VacuumBackend::Scaleless => None,
        VacuumBackend::Family => {
            return Err(internal(
                "vacuum lines",
                "every line is massive (validated)",
            ));
        }
    };
    let mut reduce = |monomial: &Atom| -> Result<Vec<(Atom, Atom)>> {
        if let Some(found) = memo.get(monomial) {
            return Ok(found.clone());
        }
        let input = monomial * &base_topology;
        let Some(cut) = &base_cut else {
            return Err(internal(
                "vacuum lines",
                "a scaleless base with a pattern to reduce (validated)",
            ));
        };
        let (external, inputs) = closure_inputs(vakint, settings, &input, cut)?;
        let reduced = external
            .into_iter()
            .zip(inputs)
            .map(|(tensor, (_, input))| (tensor, input.numerator))
            .collect::<Vec<_>>();
        memo.insert(monomial.clone(), reduced.clone());
        Ok(reduced)
    };
    let mut boundary_patterns = Vec::new();
    let mut routing = None::<Vec<(Atom, Atom)>>;
    let mut blocks_used = BTreeSet::new();
    for (raise, terms) in &patterns {
        let lines = lines_of(raise)?;
        let cut = match vacuum_backend(&lines, loops)? {
            VacuumBackend::Closure(cut) => cut,
            VacuumBackend::Scaleless => continue,
            VacuumBackend::Family => {
                return Err(internal(
                    "vacuum lines",
                    "every line is massive (validated)",
                ));
            }
        };
        if cut.line.mass_squared != closing_mass {
            return Err(internal(
                "vacuum lines",
                format!(
                    "the derived closing line {} is not the closure's",
                    cut.line.id
                ),
            ));
        }
        let to_head = |atom: &Atom| atom.replace(S.p.call(1)).with(head.clone());
        let images = cut
            .routing
            .iter()
            .map(|(k, image)| (k.clone(), to_head(image)))
            .collect::<Vec<_>>();
        match &routing {
            None => routing = Some(images),
            Some(first) if *first == images => {}
            Some(_) => {
                return Err(internal(
                    "vacuum lines",
                    "the cut routing differs between patterns",
                ));
            }
        }
        let mut pattern_terms = Vec::new();
        for (block, (structure, degree, _)) in split.iter().enumerate() {
            let mut den_sum = Atom::Zero;
            for term in terms {
                let total = degree.saturating_add(term.degree);
                if spec.taylor.is_none() || (total <= w && orders.keeps(total)) {
                    den_sum += &term.coefficient;
                }
            }
            if den_sum.is_zero() {
                continue;
            }
            let small = (structure * den_sum).expand();
            let mut coefficient = Atom::Zero;
            for (monomial, _, number) in graded_split(&small, lambda) {
                for (tensor, reduced) in reduce(&monomial)? {
                    let slot = match tensor_index.get(&tensor) {
                        Some(i) => *i,
                        None => {
                            let i = tensors.len();
                            tensor_index.insert(tensor.clone(), i);
                            tensors.push((placeholder(i), &tensor * &p2));
                            i
                        }
                    };
                    coefficient += &number * reduced * &tensors[slot].0 * &prefactor;
                }
            }
            let coefficient = coefficient.expand();
            if !coefficient.is_zero() {
                blocks_used.insert(block);
                pattern_terms.push((block, coefficient));
            }
        }
        if pattern_terms.is_empty() {
            continue;
        }
        boundary_patterns.push(BoundaryPattern {
            topology: to_head(&cut.topology),
            closure_power: cut.line.power,
            terms: pattern_terms,
        });
    }
    let split_seconds = started.elapsed().as_secs_f64();
    let blocks = split
        .into_iter()
        .map(|(_, _, block)| block.replace(d.clone()).with(in_epsilon.clone()))
        .collect::<Vec<_>>();
    if boundary_patterns.is_empty() {
        return Ok(Prepared {
            boundary: None,
            tensors,
            blocks: 0,
            reduced_monomials: memo.len(),
            contract_seconds,
            split_seconds,
        });
    }
    let boundary = injected(
        "boundary",
        injected_vakint_error,
        RemainderBoundary::new(
            value.label.clone(),
            head,
            closing_mass,
            routing.expect("one pattern"),
            blocks,
            boundary_patterns,
        ),
    )
    .map_err(|e| internal("boundary", e))?;
    Ok(Prepared {
        boundary: Some(boundary),
        tensors,
        blocks: blocks_used.len(),
        reduced_monomials: memo.len(),
        contract_seconds,
        split_seconds,
    })
}

/// Evaluate R* closure records (module documentation): every record of kind
/// [`FactorizedKind::Closure`] (validated), times `prefactors[r]` (normally
/// `records[r].prefactor.total()`), `D -> 4 - 2ε`. One boundary per record,
/// all its power patterns lowered and reduced together.
pub fn evaluate_factorized_closure(
    vakint: &Vakint,
    settings: &VakintSettings,
    records: &[FactorizedNumerator],
    prefactors: &[Atom],
) -> Result<ClosureValue> {
    // Lane unknown-masters: strict refusal once, on this call's summed value.
    if let Some(_boundary) = crate::public_result_boundary() {
        let inner = crate::unknown_master_coefficients_carried(settings);
        let result = evaluate_factorized_closure(vakint, &inner, records, prefactors)?;
        vakint
            .refuse_unknown_master_coefficients(settings, &result.vakint_value)
            .map_err(|error| FactorizedError::Evaluation {
                detail: error.to_string(),
            })?;
        return Ok(result);
    }
    if records.len() != prefactors.len() {
        return Err(refuse(format!(
            "{} records and {} prefactors",
            records.len(),
            prefactors.len()
        )));
    }
    let mut out = ClosureValue {
        value: Atom::Zero,
        vakint_value: Atom::Zero,
        records: Vec::with_capacity(records.len()),
    };
    for (value, prefactor) in records.iter().zip(prefactors) {
        let FactorizedKind::Closure(spec) = &value.kind else {
            return Err(refuse(format!(
                "{}: not a closure record (kind direct: evaluate_factorized / evaluate_factorized_vacuum)",
                value.label
            )));
        };
        value.validate()?;
        // a marked value: the separately supplied prefactor gets the checks of
        // the record's prefactor fields, the settings' epsilon the caller-symbol
        // rule (no-ops for an unmarked value: the base behaviour)
        super::check_marked_entry(value, prefactor, vk_symbol!(&settings.epsilon_symbol), None)?;
        super::check_entry_names(value, prefactor, settings)?;
        let record = evaluate_validated(vakint, settings, value, spec, prefactor)?;
        super::backstop::check_result(value, prefactor, settings, &record.value, false)?;
        out.value += &record.value;
        out.vakint_value += &record.vakint_value;
        out.records.push(record);
    }
    Ok(out)
}

/// ONE closure record, VALIDATED by the caller (its `validate` and the entry
/// check passed), times `prefactor`: the evaluation below the validation
/// (the entry's per-record body; also the tests' handle on the evaluation of
/// record shapes the validation refuses, e.g. the order filters).
pub(super) fn evaluate_validated(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    spec: &ClosureSpec,
    prefactor: &Atom,
) -> Result<ClosureRecordValue> {
    // the value must not depend on the caller's tensor-reduction method or
    // notation (WH-RECORD-FORMAT v1.4 section 7): pinned to what the R*
    // boundary path uses (FeynKit, dot notation)
    let mut local = settings.clone();
    local.tensor_reduction_method = crate::TensorReductionMethod::FeynKit;
    local.use_dot_product_notation = true;
    // the record's depth, else the settings' number of terms, by the depth
    // rule of a marked value (no cast: a depth beyond i64 is refused)
    let requested = match spec.depth {
        Some(depth) => i128::try_from(depth).unwrap_or(i128::MAX),
        None => i128::from(settings.number_of_terms_in_epsilon_expansion),
    };
    local.number_of_terms_in_epsilon_expansion = super::marked_depth(requested).map_err(refuse)?;
    // `depth` counts the terms of the RETURNED product: insertions with
    // poles of total order n promote n more orders of the integral, which
    // is evaluated n terms deeper and the product truncated back
    let epsilon_symbol = vk_symbol!(&settings.epsilon_symbol);
    value.check_insertion_epsilon(epsilon_symbol.get_name())?;
    let pole_order = value.insertion_depth_raise(epsilon_symbol)?;
    let returned_terms = local.number_of_terms_in_epsilon_expansion;
    local.number_of_terms_in_epsilon_expansion = i64::try_from(pole_order)
        .ok()
        .and_then(|n| returned_terms.checked_add(n))
        .filter(|terms| *terms <= super::MARKED_MAX_TERMS)
        .ok_or_else(|| {
            refuse(format!(
                "the depth {returned_terms} is out of range: with the insertions' raise {pole_order} a marked value asks for at most {} ε terms",
                super::MARKED_MAX_TERMS
            ))
        })?;
    let prepared = prepare(vakint, &local, value, spec, prefactor)?;
    let started = Instant::now();
    let mut counts = (0, 0, 0);
    let (mut vakint_value, lower_seconds, reduce_seconds, patterns) = match &prepared.boundary {
        None => (Atom::Zero, 0.0, 0.0, 0),
        Some(boundary) => {
            let direct = boundary.patterns.iter().all(
                |p| matches!(p.topology.as_view(), AtomView::Fun(f) if f.get_symbol() == S.topo),
            );
            if direct {
                let (atom, stats) = injected(
                    "lowering and reduction",
                    injected_vakint_error,
                    boundary.evaluate_with_stats(vakint, &local, BoundaryOptions::default()),
                )
                .map_err(|e| internal("lowering and reduction", e))?;
                let lower = stats.match_seconds
                    + stats.read_seconds
                    + stats.classify_seconds
                    + stats.slot_seconds
                    + stats.finish_seconds;
                counts = (
                    stats.key_plans,
                    stats.lowered_entries,
                    stats.distinct_monomials,
                );
                (
                    atom,
                    lower,
                    stats.reduce_finalize_seconds,
                    boundary.patterns.len(),
                )
            } else {
                // a loop-free cut (one loop): the ordinary closure per pattern
                let mut atom = Atom::Zero;
                for (closure, input) in boundary.ordinary() {
                    atom += injected(
                        "ordinary closure",
                        injected_vakint_error,
                        closure.evaluate(vakint, &local, input.as_view()),
                    )
                    .map_err(|e| internal("ordinary closure", e))?;
                }
                (atom, 0.0, 0.0, boundary.patterns.len())
            }
        }
    };
    for (symbol, tensor) in &prepared.tensors {
        vakint_value = vakint_value.replace(symbol.clone()).with(tensor.clone());
    }
    // the product truncated back to the requested depth; in Retain the masters
    // are symbols and the value is the exact product: no series
    let vakint_value = if pole_order > 0 && !retains_masters(settings) {
        let target = i64::try_from(value.loops())
            .ok()
            .and_then(|loops| returned_terms.checked_sub(loops))
            .and_then(|t| t.checked_sub(1))
            .ok_or_else(|| {
                refuse(format!(
                    "the depth {returned_terms} is out of range: no truncation order for {} loops",
                    value.loops()
                ))
            })?;
        injected(
            "truncation of the insertion product",
            || symbolica::poly::series::SeriesError::NonPositiveRelativeDepth {
                depth: symbolica::poly::series::SeriesDepth::Relative(
                    symbolica::domains::rational::Rational::from(0),
                ),
            },
            vakint_value.series(
                epsilon_symbol,
                Atom::Zero.as_view(),
                symbolica::domains::rational::Rational::from(target),
            ),
        )
        .map_err(|e| internal("truncation of the insertion product", e))?
        .to_atom()
        .expand()
    } else {
        vakint_value.expand()
    };
    let evaluate_seconds = started.elapsed().as_secs_f64();
    // the caller's notation; with no open slot, scalar products of the
    // external momenta still come back as spenso::dot(Q(x, mink(D)), ..):
    // the representation is taken from a placeholder slot
    let mut caller_view = opened(value)?;
    if caller_view.ports.is_empty() {
        caller_view.ports.push(function!(
            symbolica::symbol!(names::MINK),
            Atom::var(value.dimension),
            Atom::num(0)
        ));
    }
    let caller = vacuum::to_caller_notation(&caller_view, &vakint_value);
    Ok(ClosureRecordValue {
        value: caller,
        vakint_value,
        patterns,
        blocks: prepared.blocks,
        reduced_monomials: prepared.reduced_monomials,
        contract_seconds: prepared.contract_seconds,
        split_seconds: prepared.split_seconds,
        evaluate_seconds,
        lower_seconds,
        reduce_seconds,
        key_plans: counts.0,
        lowered_entries: counts.1,
        distinct_monomials: counts.2,
    })
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod internal_failure_tests {
    use super::*;

    /// An internal failure is `Evaluation` with the stage named; a designed
    /// refusal stays `Closure`.
    #[test]
    fn internal_failures_are_evaluation_errors_and_refusals_stay_closure() {
        for stage in [
            "boundary",
            "lowering and reduction",
            "ordinary closure",
            "truncation of the insertion product",
        ] {
            match internal(stage, "inner failure") {
                FactorizedError::Evaluation { detail } => {
                    assert!(
                        detail.starts_with(&format!("closure {stage}: ")),
                        "{detail}"
                    )
                }
                other => panic!("{other:?}"),
            }
        }
        assert!(matches!(refuse("x"), FactorizedError::Closure { .. }));
    }

    /// The four internal-failure CALL SITES use `internal(..)`, not
    /// `refuse(..)` (none can be triggered from a valid record by a cheap
    /// input: the boundary, the lowering / reduction, the ordinary closure and
    /// the truncation fail only on an inconsistent internal state): a site
    /// reverted to a refusal fails this test.
    #[test]
    fn the_internal_failure_sites_return_evaluation() {
        // the whole file: the patterns below appear in this test only escaped
        let code = include_str!("closure.rs");
        for stage in [
            "internal(\"boundary\", e)",
            "internal(\"lowering and reduction\", e)",
            "internal(\"ordinary closure\", e)",
            "internal(\"truncation of the insertion product\", e)",
        ] {
            assert_eq!(code.matches(stage).count(), 1, "{stage}");
        }
        let refusal = ["map_err(|e| ", "refuse("].concat();
        assert!(
            !code.contains(&refusal),
            "an evaluation failure mapped to a refusal"
        );
    }
}
