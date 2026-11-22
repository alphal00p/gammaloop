//! U7c: K5 insertion-surrogate records (kind `surrogate` of
//! `GAMMALOOP_UNCONTRACTED_NUMERATOR`; WH-RECORD-FORMAT v1.4 section 2 row K5,
//! sections 5 and 7; the beta's fixtures k5-direct, k5-surrogate, k5-stress).
//!
//! An insertion surrogate ([`super::super::insertions`]) replaces each
//! selected two-point region of a diagram by an L-bubble scalar chain
//! (`scalar_insertion_0..L`) whose first vertex carries the region's
//! normalised coefficient. Its records:
//!
//! * <= 4 loops: ONE `kind direct` record, the surrogate's direct propagator
//!   integral (`integrate_propagator(surrogate)`, insertions.rs): every
//!   internal line massless, no Taylor step, no projector (both root slots
//!   stay open), depth = loops; comparison point = the return of
//!   `integrate_propagator` (before `restore_atom`), ours = vakint's
//!   `evaluate_factorized_ports` then `finish_tensor(value, 0)` ([`direct`]).
//! * above four loops: the remainders of the surrogate's own R*
//!   (`closed_pole_part` / `pole_part` inside the dressed scope) are built by
//!   the K1 / K2 emitter ([`super::build`] with [`frame`]); K3-shaped IR hard
//!   terms inside a surrogate decline by design (`surrogate_ir_hard`) before
//!   any record is built; K7 and K4 stay contracted (no kind).
//!
//! In every record ([`finish`]): a chain whose coefficient vertex chain[0] is
//! a record vertex becomes `vertex v rule=ins:<rule>` factor 1 with
//! `insertion v hedges=<inner half-edge at chain[0]>,<inner half-edge at
//! chain[L]> cycles=<the excised region's closed cycles> 1` and the vertex's
//! factor verbatim as `insertion_factor`; bubble edges are `line=scalar`
//! (factor 1, checked); `original_cycles` is the ORIGINAL diagram's
//! ([`provenance`], computed at `prepare`); a `ct:` child gets the closed
//! cycles of its own support (bubble lines excluded, each chain one node) plus
//! those of the holes whose chain[0] it contracts, and a hole is never both.
//!
//! A surrogate evaluated under the enabled kind is a SITE of its diagram's
//! ledger ([`super::super::uncontracted_ledger::k5_site`]): outside every
//! diagram ledger it is an error (no implicit ledger), and a diagram with a
//! site whose kind expected no record fails. Errors (never designed):
//! `surrogate_context_missing`, `surrogate_chain`, `surrogate_coefficient_mass`,
//! `surrogate_hole_twice`, and every U7a error code.

use std::{
    cell::RefCell,
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
    time::Instant,
};

use eyre::{Result, eyre};
use linnet::half_edge::{
    NodeIndex,
    involution::{EdgeIndex, Flow, HedgePair},
    subgraph::{SuBitGraph, SubSetOps},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
};
use vakint::{
    Vakint,
    p_integral::factorized::{
        DirectRecord, Endpoint, FactorizedEdge, FactorizedKind, FactorizedNumerator,
        FactorizedVertex, Insertion, LineKind, Prefactor, evaluate_factorized_ports,
    },
};

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::{FeynmanGraph, Graph},
    model::{Model, VertexRuleIdGammaLoopExt},
    numerator::symbolica_ext::NumeratorAtomExt,
    utils::GS,
    uv::{UVgenerationSettings, UltravioletGraph},
};

use super::super::{
    GraphRegion, expression::RStarIntegrand, top_mass::TopMassTerms, uncontracted_ledger,
};
use super::{
    Attempt, Decline, EdgeInfo, Frame, Kind, Mode, Outcome, Site, check_failure, child_frame,
    closed_cycles, compare, edge_table, error_failure, insertion_masses, is_budget_text,
    is_unknown_master_text, mentions_mass, q_coefficients, report, slots_of, switch,
    tally_unknown_masters, unknown_masters_failure, valid_order,
};

// ---------------------------------------------------------------------------
// Provenance (computed at prepare) and the surrogate context
// ---------------------------------------------------------------------------

/// What a surrogate's records need from the ORIGINAL diagram: its closed
/// cycles and, per hole, the coefficient vertex and the closed cycles inside
/// the excised region.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Provenance {
    /// The original diagram's name (diagnostics).
    pub original: String,
    /// `(fermion, ghost)` closed cycles of the original diagram.
    pub original_cycles: (u32, u32),
    pub holes: Vec<Hole>,
}

/// One hole of a surrogate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Hole {
    /// chain[0], the coefficient vertex (a vertex id of the surrogate graph).
    pub vertex: usize,
    /// L, the bubbles of the chain (= the region's loops).
    pub loops: usize,
    /// `(fermion, ghost)` closed cycles inside the excised original region.
    pub cycles: (u32, u32),
}

/// The provenance of a surrogate of `diagram` whose holes are `holes`
/// (`(coefficient vertex, loops, the region in the original graph)`), when
/// the switch enables `surrogate`; `None` otherwise (nothing is computed).
pub(crate) fn provenance<'h>(
    diagram: &feynkit_graph::FeynmanDiagram,
    holes: impl IntoIterator<Item = (usize, usize, &'h SuBitGraph)>,
) -> Result<Option<Arc<Provenance>>> {
    if !switch().enabled(Kind::Surrogate) {
        return Ok(None);
    }
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let table = edge_table(&graph, diagram.model());
    let internal = table
        .iter()
        .filter(|(_, info)| info.source.is_some() && info.target.is_some())
        .map(|(edge, _)| *edge)
        .collect::<BTreeSet<_>>();
    let cycles_of = |edges: &BTreeSet<usize>| {
        (
            closed_cycles(&table, edges, LineKind::Fermion),
            closed_cycles(&table, edges, LineKind::Ghost),
        )
    };
    let mut out = Vec::new();
    for (vertex, loops, support) in holes {
        let region = graph
            .iter_edges_of(support)
            .filter(|(pair, _, _)| matches!(pair, HedgePair::Paired { .. }))
            .map(|(_, edge, _)| usize::from(edge))
            .collect::<BTreeSet<_>>();
        if region.is_empty() || !region.is_subset(&internal) {
            return Err(eyre!(
                "U7c: the excised region {region:?} of {} (hole at surrogate vertex {vertex}) is not a set of internal edges of the original graph",
                graph.name
            ));
        }
        out.push(Hole {
            vertex,
            loops,
            cycles: cycles_of(&region),
        });
    }
    Ok(Some(Arc::new(Provenance {
        original: graph.name.clone(),
        original_cycles: cycles_of(&internal),
        holes: out,
    })))
}

thread_local! {
    static CONTEXT: RefCell<Option<Arc<Provenance>>> = const { RefCell::new(None) };
}

/// The surrogate context of this thread (set by [`scoped`] / [`direct`]).
fn context() -> Option<Arc<Provenance>> {
    CONTEXT.with(|slot| slot.borrow().clone())
}

/// Run `f` with `provenance` as this thread's surrogate context (restored
/// on unwinding too).
fn with_context<T>(provenance: Arc<Provenance>, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<Arc<Provenance>>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0.take();
            let _ = CONTEXT.try_with(|slot| *slot.borrow_mut() = previous);
        }
    }
    let _restore = Restore(CONTEXT.with(|slot| slot.replace(Some(provenance))));
    f()
}

/// The R* of a surrogate above four loops (`graph` its name; `closed`:
/// its root is closed, [`super::super::RStarPlan::closed_pole_part`]): with
/// the kind enabled, a site of the diagram's ledger and the surrogate context
/// for its records; otherwise just `f`. A setting under which every
/// surrogate root record would decline ([`settings_refusal`]) is an error
/// here, and so is a site that adds no expected record of the kind (its
/// records missed the hook: [`uncontracted_ledger::expected_now`]).
pub(crate) fn scoped<T>(
    provenance: Option<&Arc<Provenance>>,
    graph: &str,
    closed: bool,
    uv: &UVgenerationSettings,
    f: impl FnOnce() -> Result<T>,
) -> Result<T> {
    if !switch().enabled(Kind::Surrogate) {
        return f();
    }
    let provenance = provenance.ok_or_else(|| missing_context(graph))?;
    settings_refusal(graph, Some(closed), uv)?;
    uncontracted_ledger::k5_site(Kind::Surrogate, graph)?;
    let before = uncontracted_ledger::expected_now(Kind::Surrogate);
    let out = with_context(provenance.clone(), f)?;
    let after = uncontracted_ledger::expected_now(Kind::Surrogate);
    if after <= before {
        return Err(eyre!(
            "UNCONTRACTED kind=surrogate graph={graph} outcome=error:surrogate_site_without_records: the surrogate's R* reached no record of the kind (a missed hook or a misclassified record)"
        ));
    }
    Ok(out)
}

/// The settings under which EVERY record of a surrogate would decline by
/// design (pre-review of U7c, finding 3: a setting-wide decline must not
/// hide behind the children's comparisons): an open surrogate root
/// (`GAMMALOOP_SURROGATE_CLOSED_ROOT=0`, a non-metric projection or a root
/// that does not close; `closed` is `None` for a direct record), top-mass
/// terms other than keep, an `additional_normalization` other than 1. An
/// error with the kind enabled, never a fallback.
fn settings_refusal(graph: &str, closed: Option<bool>, uv: &UVgenerationSettings) -> Result<()> {
    #[cfg(test)]
    if !SETTINGS_CHECK.with(std::cell::Cell::get) {
        return Ok(());
    }
    let refuse = |why: String| {
        Err(eyre!(
            "UNCONTRACTED kind=surrogate graph={graph} outcome=error:surrogate_settings: {why}: every record of the surrogate would decline by design; run the kind with the production settings or drop it from the list"
        ))
    };
    if closed == Some(false) {
        return refuse(
            "the surrogate root is not closed (GAMMALOOP_SURROGATE_CLOSED_ROOT=0, a non-metric projection or a root that does not close)"
                .into(),
        );
    }
    if uv.vakint.additional_normalization.trim() != "1" {
        return refuse(format!(
            "additional_normalization is {}, not 1",
            uv.vakint.additional_normalization
        ));
    }
    if closed.is_some() {
        let mode = super::super::top_mass::top_mass_terms()?;
        if mode != TopMassTerms::Keep {
            return refuse(format!("top-mass terms {mode:?}, not keep"));
        }
    }
    Ok(())
}

fn missing_context(graph: &str) -> eyre::Report {
    eyre!(
        "UNCONTRACTED kind=surrogate graph={graph} outcome=error:surrogate_context_missing: the surrogate was prepared without its provenance (prepare ran with the kind not enabled)"
    )
}

// ---------------------------------------------------------------------------
// GAMMALOOP_UNCONTRACTED_K5_HOLE_CUT: how a kept and cut hole is written
// ---------------------------------------------------------------------------

/// The switch of a kept and cut hole (a record that holds a hole's anchor
/// but not its far end): `form` writes the single-hedge `ins:` with the far
/// end's slot as a spectator port (WH-RECORD-FORMAT v1.5 proposal (a), beta's
/// ruling of 21 Nov; for the vakint step after G-D), `decline` (unset, the
/// DEFAULT: the v1.5 rules stay out of the G-D gate) declines such a record
/// by design (`surrogate_hole_cut`, today's path), `auto` probes the reader
/// once: `form` when vakint's reader accepts a single-hedge `ins:`, else
/// `decline` (the gated reader ccba5e0c8a requires exactly two hedges).
pub(crate) const HOLE_CUT_VARIABLE: &str = "GAMMALOOP_UNCONTRACTED_K5_HOLE_CUT";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum HoleCut {
    Form,
    Decline,
}

/// The hole-cut mode of the value `value` of [`HOLE_CUT_VARIABLE`] (`None`:
/// unset). `probe`: whether `auto` may probe the reader (the kind
/// `surrogate` enabled); without it `auto` stands for `decline` (the mode is
/// then never used). Production calls it only through [`resolve_hole_cut`].
pub(crate) fn parse_hole_cut(
    value: Option<&str>,
    probe: bool,
) -> std::result::Result<HoleCut, String> {
    match value.map(str::trim) {
        None | Some("" | "decline") => Ok(HoleCut::Decline),
        Some("auto") => Ok(if probe && reader_accepts_single_hedge() {
            HoleCut::Form
        } else {
            HoleCut::Decline
        }),
        Some("form") => Ok(HoleCut::Form),
        Some(other) => Err(format!(
            "{HOLE_CUT_VARIABLE} must be decline (the default), form or auto, not `{other}`"
        )),
    }
}

/// The name of a hole-cut mode (identity, cache key, lines).
pub(crate) fn hole_cut_name(mode: HoleCut) -> &'static str {
    match mode {
        HoleCut::Form => "form",
        HoleCut::Decline => "decline",
    }
}

/// The process's hole-cut mode: written ONCE, by [`resolve_hole_cut`] at the
/// switch's initialisation (uncontracted.rs `ENVIRONMENT`), never again
/// (a `OnceLock`: no setter, no second write).
static RESOLVED_HOLE_CUT: std::sync::OnceLock<HoleCut> = std::sync::OnceLock::new();

/// Resolve the hole-cut mode of the value `value` of [`HOLE_CUT_VARIABLE`]
/// and store it for the process: called exactly once, at the switch's
/// initialisation (a bad value refuses the run at start-up; U7c delta
/// pre-review N3); a second resolution is an error.
pub(crate) fn resolve_hole_cut(
    value: Option<&str>,
    probe: bool,
) -> std::result::Result<HoleCut, String> {
    let mode = parse_hole_cut(value, probe)?;
    RESOLVED_HOLE_CUT.set(mode).map_err(|_| {
        format!(
            "the hole-cut mode ({HOLE_CUT_VARIABLE}) is resolved once, at the switch's initialisation: a second resolution is refused"
        )
    })?;
    Ok(mode)
}

/// The hole-cut mode at the switch's start-up (uncontracted.rs
/// `ENVIRONMENT`): read (`read`) and resolved ONLY when the kind `surrogate`
/// is enabled; otherwise the variable is not read at all and nothing is
/// resolved (a K5 variable never affects a run without the kind, as with
/// the switch off). A bad value is then an error only for a K5 run.
pub(crate) fn start_up_hole_cut(
    switch: super::Switch,
    read: impl FnOnce() -> std::result::Result<Option<String>, String>,
) -> std::result::Result<Option<HoleCut>, String> {
    if switch.mode == Mode::Off || !switch.enabled(Kind::Surrogate) {
        return Ok(None);
    }
    resolve_hole_cut(read()?.as_deref(), true).map(Some)
}

/// The mode resolved at start-up (`None`: not resolved: the switch off, or
/// the kind `surrogate` not enabled).
pub(crate) fn resolved_hole_cut() -> Option<HoleCut> {
    RESOLVED_HOLE_CUT.get().copied()
}

#[cfg(test)]
thread_local! {
    static TEST_HOLE_CUT: std::cell::Cell<Option<HoleCut>> = const { std::cell::Cell::new(None) };
}

/// The mode in force: the one resolved at the switch's initialisation
/// (constant for the process; tests may override it on a thread). It is
/// part of the resume identity; it is NOT part of the Delta-cache key: the
/// cache is in-process only (CountertermCache, uv/counterterm_cache.rs, one
/// per engine run), so a per-process constant cannot split its fills.
pub(crate) fn hole_cut() -> HoleCut {
    #[cfg(test)]
    if let Some(mode) = TEST_HOLE_CUT.with(std::cell::Cell::get) {
        return mode;
    }
    super::environment_hole_cut()
}

/// Tests: run `f` with the hole-cut mode `mode` on this thread.
#[cfg(test)]
pub(crate) fn with_hole_cut<T>(mode: HoleCut, f: impl FnOnce() -> T) -> T {
    struct Restore(Option<HoleCut>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            TEST_HOLE_CUT.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(TEST_HOLE_CUT.with(|slot| slot.replace(Some(mode))));
    f()
}

/// Whether vakint's reader accepts an `ins:` insertion with ONE hedge (a
/// minimal one-loop direct record: two scalar lines between an `ins:`
/// vertex and a plain one, a leg at each).
pub(crate) fn reader_accepts_single_hedge() -> bool {
    static ACCEPTS: std::sync::LazyLock<bool> = std::sync::LazyLock::new(|| {
        let leg = |id, source, target, source_hedge, target_hedge| FactorizedEdge {
            id,
            source,
            target,
            source_hedge,
            target_hedge,
            particle: "c".into(),
            line: LineKind::Scalar,
            power: 1,
            mass_squared: Atom::Zero,
            momentum: vec![0, 1],
            factor: Atom::num(1),
        };
        let line = |id, source_hedge, target_hedge, momentum: Vec<i64>| FactorizedEdge {
            id,
            source: Endpoint::Vertex(0),
            target: Endpoint::Vertex(1),
            source_hedge: Some(source_hedge),
            target_hedge: Some(target_hedge),
            particle: "c".into(),
            line: LineKind::Scalar,
            power: 1,
            mass_squared: Atom::Zero,
            momentum,
            factor: Atom::num(1),
        };
        let record = FactorizedNumerator {
            label: "u7c reader probe: single-hedge ins:".into(),
            kind: FactorizedKind::Direct,
            dimension: GS.dim,
            loop_edges: vec![Some(2)],
            external_edges: vec![Some(0)],
            edges: vec![
                leg(0, Endpoint::External, Endpoint::Vertex(0), None, Some(0)),
                leg(1, Endpoint::Vertex(1), Endpoint::External, Some(1), None),
                line(2, 2, 3, vec![1, 0]),
                line(3, 4, 5, vec![-1, 1]),
            ],
            vertices: vec![
                FactorizedVertex {
                    id: 0,
                    rule: "ins:scalar_insertion_0".into(),
                    factor: Atom::num(1),
                },
                FactorizedVertex {
                    id: 1,
                    rule: "scalar_insertion_1".into(),
                    factor: Atom::num(1),
                },
            ],
            ports: Vec::new(),
            projector: None,
            prefactor: Prefactor {
                overall: Atom::num(1),
                loop_sign: 1,
                nf_power: 0,
                scalars: Atom::num(1),
                colour: Atom::num(1),
            },
            colour_branch: None,
            insertions: vec![Insertion {
                vertex: 0,
                hedges: vec![0],
                cycles: (0, 0),
                coefficient: Atom::num(1),
                factor: Atom::num(1),
            }],
            record: Some(DirectRecord {
                original_cycles: Some((0, 0)),
                depth: Some(1),
                weights: Vec::new(),
            }),
        };
        record
            .to_text()
            .ok()
            .is_some_and(|text| FactorizedNumerator::from_text(&text).is_ok())
    });
    *ACCEPTS
}

// ---------------------------------------------------------------------------
// The frame and the K5 rewrite of a K1 / K2 record
// ---------------------------------------------------------------------------

/// The frame of a record of kind `surrogate` at the ordinary hook: a
/// surrogate root remainder is K1-shaped, a remainder of a UV child of the
/// surrogate K2-shaped; an IR hard term declines by design.
pub(super) fn frame<'s>(
    site: &Site<'s>,
    remainder: &RStarIntegrand,
) -> std::result::Result<Frame<'s>, Decline> {
    if site.ir_hard.is_some() {
        return Err(Decline::new(
            "surrogate_ir_hard",
            "a K3-shaped IR hard term inside an insertion surrogate's R* (no verified fixture; today's path)",
        ));
    }
    if context().is_none() {
        return Err(Decline::new(
            "surrogate_context_missing",
            "a surrogate record outside the surrogate context",
        ));
    }
    match site.plan.root_boundary {
        None => {
            let (Some(closure), true) = (site.root_closure, site.root_projector) else {
                return Err(Decline::new(
                    "open_root",
                    "a surrogate root without its root closure or projector (not closed)",
                ));
            };
            let Some(closed_root) = site.closed_root else {
                return Err(Decline::new(
                    "root_normalisation_opaque",
                    "a surrogate root projector without the closed-root normalisation as data",
                ));
            };
            if site.top_mass != Some(TopMassTerms::Keep) {
                return Err(Decline::new(
                    "top_mass_mode",
                    format!("top-mass terms {:?}", site.top_mass),
                ));
            }
            Ok(Frame::Root {
                closure,
                closed_root,
            })
        }
        Some(_) => child_frame(site, remainder),
    }
}

/// One scalar chain of the surrogate graph.
#[derive(Clone, Debug)]
pub(super) struct Chain {
    /// chain[0..=L] (vertex ids).
    pub vertices: Vec<usize>,
    /// The bubble edges.
    pub bubbles: BTreeSet<usize>,
    /// The inner boundary half-edges at chain[0] and chain[L].
    pub inner: [usize; 2],
    /// Their edges (e_in at chain[0], e_out at chain[L]) and whether each
    /// half-edge is its edge's sink (the edge enters the chain there).
    pub inner_edges: [(usize, bool); 2],
    /// The closed cycles inside the excised region.
    pub cycles: (u32, u32),
    /// Whether the coefficient C_H (chain[0]'s factor in the surrogate) has
    /// a slot of the far end's inner half-edge `inner[1]` (a gluon or quark
    /// carrier; a ghost's has none).
    pub far_slot: bool,
}

fn chain_error(detail: impl Into<String>) -> Decline {
    Decline::new("surrogate_chain", detail)
}

/// The chains of `provenance`'s holes in the surrogate `graph`: from chain[0]
/// (`scalar_insertion_0`) through `scalar_insertion_i`, each link two
/// parallel edges with factor 1; the inner half-edges are the one non-bubble,
/// non-dummy half-edge at chain[0] and at chain[L].
pub(super) fn chains(
    graph: &Graph,
    provenance: &Provenance,
) -> std::result::Result<Vec<Chain>, Decline> {
    let mut at = BTreeMap::<usize, Vec<(usize, usize)>>::new();
    let mut ends = BTreeMap::<usize, (usize, usize)>::new();
    let mut nodes = BTreeMap::<usize, NodeIndex>::new();
    let mut dummy = BTreeSet::new();
    let mut unit = BTreeSet::new();
    let mut hedge_edge = BTreeMap::<usize, (usize, bool)>::new();
    for (pair, edge, data) in graph.iter_edges() {
        let e = usize::from(edge);
        if data.data.is_dummy {
            dummy.insert(e);
        }
        if data.data.num.value == Atom::num(1) {
            unit.insert(e);
        }
        let hedges = match pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                hedge_edge.insert(usize::from(source), (e, false));
                hedge_edge.insert(usize::from(sink), (e, true));
                let (s, t) = (graph.node_id(source), graph.node_id(sink));
                nodes.insert(usize::from(s), s);
                nodes.insert(usize::from(t), t);
                ends.insert(e, (usize::from(s), usize::from(t)));
                vec![source, sink]
            }
            HedgePair::Unpaired { hedge, flow } => {
                hedge_edge.insert(usize::from(hedge), (e, flow == Flow::Sink));
                vec![hedge]
            }
        };
        for hedge in hedges {
            let node = graph.node_id(hedge);
            nodes.insert(usize::from(node), node);
            at.entry(usize::from(node))
                .or_default()
                .push((usize::from(hedge), e));
        }
    }
    let name = |vertex: usize| {
        nodes
            .get(&vertex)
            .map(|node| graph[*node].name.value.clone())
            .unwrap_or_default()
    };
    let mut out = Vec::new();
    for hole in &provenance.holes {
        if name(hole.vertex) != "scalar_insertion_0" {
            return Err(chain_error(format!(
                "the coefficient vertex {} is `{}`, not scalar_insertion_0",
                hole.vertex,
                name(hole.vertex)
            )));
        }
        let mut vertices = vec![hole.vertex];
        let mut bubbles = BTreeSet::new();
        for i in 0..hole.loops {
            let current = vertices[i];
            let next_name = format!("scalar_insertion_{}", i + 1);
            let mut links = BTreeMap::<usize, Vec<usize>>::new();
            for (_, e) in at.get(&current).into_iter().flatten() {
                if let Some((a, b)) = ends.get(e) {
                    let other = if *a == current { *b } else { *a };
                    if other != current && !vertices.contains(&other) && name(other) == next_name {
                        links.entry(other).or_default().push(*e);
                    }
                }
            }
            let mut links = links.into_iter().filter(|(_, edges)| edges.len() == 2);
            let (Some((next, edges)), None) = (links.next(), links.next()) else {
                return Err(chain_error(format!(
                    "no unique {next_name} with two bubble edges after chain vertex {current} (hole at {})",
                    hole.vertex
                )));
            };
            for e in &edges {
                if !unit.contains(e) {
                    return Err(chain_error(format!(
                        "bubble edge {e} of the chain at {} has a factor other than 1",
                        hole.vertex
                    )));
                }
            }
            bubbles.extend(edges);
            vertices.push(next);
        }
        let inner_at = |vertex: usize| -> std::result::Result<usize, Decline> {
            let found = at
                .get(&vertex)
                .into_iter()
                .flatten()
                .filter(|(_, e)| !bubbles.contains(e) && !dummy.contains(e))
                .map(|(h, _)| *h)
                .collect::<Vec<_>>();
            match found[..] {
                [h] => Ok(h),
                _ => Err(chain_error(format!(
                    "chain vertex {vertex} (hole at {}) has the inner half-edges {found:?}, not one",
                    hole.vertex
                ))),
            }
        };
        let inner = [inner_at(vertices[0])?, inner_at(*vertices.last().unwrap())?];
        let edge_of = |h: usize| {
            hedge_edge
                .get(&h)
                .copied()
                .ok_or_else(|| chain_error(format!("half-edge {h} has no edge")))
        };
        let inner_edges = [edge_of(inner[0])?, edge_of(inner[1])?];
        let far_slot = nodes
            .get(&vertices[0])
            .is_some_and(|node| hedge_slots(&graph[*node].get_num()).contains_key(&inner[1]));
        out.push(Chain {
            vertices,
            bubbles,
            inner,
            inner_edges,
            cycles: hole.cycles,
            far_slot,
        });
    }
    Ok(out)
}

/// The K5 rewrite of a record of kind `surrogate` (built by the K1 / K2
/// emitter, or by [`build_direct`]) in the surrogate context: see the module
/// documentation. Returns the record's `original_cycles`.
#[allow(clippy::too_many_arguments)]
pub(super) fn finish(
    graph: &Graph,
    table: &BTreeMap<usize, EdgeInfo>,
    region: &GraphRegion,
    merge: &BTreeMap<usize, usize>,
    nodes: &[(usize, usize)],
    edges: &mut [FactorizedEdge],
    vertices: &mut [FactorizedVertex],
    insertions: &mut Vec<Insertion>,
) -> std::result::Result<Finished, Decline> {
    let provenance = context().ok_or_else(|| {
        Decline::new(
            "surrogate_context_missing",
            "a surrogate record outside the surrogate context",
        )
    })?;
    finish_with(
        &provenance,
        graph,
        table,
        region,
        merge,
        nodes,
        edges,
        vertices,
        insertions,
        false,
    )
}

/// How `atom` depends on the symbol `eps` (WH-RECORD-FORMAT v1.4
/// correction 4: only through sums, products and integer powers; Minkowski
/// and bispinor slot arguments aside): `Ok(true)` polynomially, `Ok(false)`
/// through a negative integer power, `Err` otherwise (inside a function, a
/// non-integer power).
fn eps_dependence(atom: AtomView<'_>, eps: Symbol) -> std::result::Result<bool, String> {
    let contains = |view: AtomView<'_>| view.to_owned().contains_symbol(eps);
    match atom {
        AtomView::Num(_) | AtomView::Var(_) => Ok(true),
        AtomView::Add(sum) => sum.iter().try_fold(true, |polynomial, term| {
            Ok(eps_dependence(term, eps)? && polynomial)
        }),
        AtomView::Mul(product) => product.iter().try_fold(true, |polynomial, factor| {
            Ok(eps_dependence(factor, eps)? && polynomial)
        }),
        AtomView::Pow(power) => {
            let (base, exponent) = power.get_base_exp();
            if !contains(base) {
                return if contains(exponent) {
                    Err(format!("eps in the exponent of {}", atom))
                } else {
                    Ok(true)
                };
            }
            let Ok(n) = i64::try_from(exponent) else {
                return Err(format!("a non-integer power {}", atom));
            };
            Ok(eps_dependence(base, eps)? && n >= 0)
        }
        AtomView::Fun(function) => {
            let symbol = function.get_symbol();
            let name = symbol.get_name();
            if name.ends_with("::mink") || name.ends_with("::bis") || !contains(atom) {
                Ok(true)
            } else {
                Err(format!("eps inside the function {}", atom))
            }
        }
    }
}

/// The record's id of the edge that owns the hole anchor's inner half-edge
/// `hedge` (of the graph edge `edge`, e_in): `edge` itself, or the incoming
/// leg `PAIRED_LEG_OFFSET + edge` of a paired crown edge whose sink is the
/// anchor (uncontracted.rs `paired_leg_id`); any other owner, or none, is
/// the error `surrogate_chain`.
pub(crate) fn anchor_leg(
    edges: &[FactorizedEdge],
    hedge: usize,
    edge: usize,
) -> std::result::Result<usize, Decline> {
    match edges
        .iter()
        .find(|e| e.source_hedge == Some(hedge) || e.target_hedge == Some(hedge))
    {
        Some(owner) if owner.id == edge || owner.id == super::PAIRED_LEG_OFFSET + edge => {
            Ok(owner.id)
        }
        Some(owner) => Err(chain_error(format!(
            "the anchor's inner half-edge {hedge} belongs to the record edge {}, not to {edge} or its paired leg",
            owner.id
        ))),
        None => Err(chain_error(format!(
            "the anchor's inner half-edge {hedge} is not a half-edge of the record"
        ))),
    }
}

/// C_H written in the anchor's record edge `anchor_edge` ([`anchor_leg`]):
/// `Q(e_out)` -> `s Q(anchor_edge)` (the chain conserves the momentum, `s = +1`
/// when one end enters the chain and the other leaves it) and, when the
/// anchor's leg is a paired crown edge's incoming leg, `Q(e_in)` ->
/// `Q(anchor_edge)` (both legs carry e_in's momentum in its orientation).
pub(crate) fn in_anchor_edge(
    factor: Atom,
    inner_edges: [(usize, bool); 2],
    anchor_edge: usize,
) -> Atom {
    let [(e_in, in_sink), (e_out, out_sink)] = inner_edges;
    if !factor.contains_symbol(GS.emr_mom) {
        return factor;
    }
    let sign = if in_sink != out_sink { 1 } else { -1 };
    let mut replacements = Vec::new();
    if e_out != e_in {
        replacements.push(symbolica::id::Replacement::new(
            function!(GS.emr_mom, e_out, crate::utils::W_.x___).to_pattern(),
            (Atom::num(sign) * function!(GS.emr_mom, anchor_edge, crate::utils::W_.x___))
                .to_pattern(),
        ));
    }
    if anchor_edge != e_in {
        replacements.push(symbolica::id::Replacement::new(
            function!(GS.emr_mom, e_in, crate::utils::W_.x___).to_pattern(),
            function!(GS.emr_mom, anchor_edge, crate::utils::W_.x___).to_pattern(),
        ));
    }
    if replacements.is_empty() {
        factor
    } else {
        factor.replace_multiple(&replacements)
    }
}

#[allow(clippy::too_many_arguments)]
fn finish_with(
    provenance: &Provenance,
    graph: &Graph,
    table: &BTreeMap<usize, EdgeInfo>,
    region: &GraphRegion,
    merge: &BTreeMap<usize, usize>,
    nodes: &[(usize, usize)],
    edges: &mut [FactorizedEdge],
    vertices: &mut [FactorizedVertex],
    insertions: &mut Vec<Insertion>,
    direct: bool,
) -> std::result::Result<Finished, Decline> {
    let chains = chains(graph, provenance)?;
    let bubbles = chains
        .iter()
        .flat_map(|chain| chain.bubbles.iter().copied())
        .collect::<BTreeSet<_>>();
    // bubble lines: scalar, factor 1 (INSERT FR-1)
    for edge in edges.iter_mut() {
        if bubbles.contains(&edge.id) {
            if edge.factor != Atom::num(1) {
                return Err(chain_error(format!(
                    "bubble edge {} carries the factor {}",
                    edge.id, edge.factor
                )));
            }
            edge.line = LineKind::Scalar;
        }
    }
    // the cycle table: bubble lines excluded, each chain one node
    let mut node_of = BTreeMap::new();
    for chain in &chains {
        for vertex in &chain.vertices {
            node_of.insert(*vertex, chain.vertices[0]);
        }
    }
    let map = |vertex: Option<usize>| vertex.map(|v| *node_of.get(&v).unwrap_or(&v));
    let merged = table
        .iter()
        .map(|(edge, info)| {
            (
                *edge,
                EdgeInfo {
                    line: if bubbles.contains(edge) {
                        LineKind::Scalar
                    } else {
                        info.line
                    },
                    particle: String::new(),
                    source: map(info.source),
                    target: map(info.target),
                },
            )
        })
        .collect::<BTreeMap<_, _>>();
    let record_hedges = edges
        .iter()
        .flat_map(|edge| [edge.source_hedge, edge.target_hedge])
        .flatten()
        .collect::<BTreeSet<_>>();
    let masses = insertion_masses(edges);
    let mut free = BTreeSet::new();
    for chain in &chains {
        let start = chain.vertices[0];
        if merge.contains_key(&start) {
            continue;
        }
        let Some(vertex) = vertices.iter_mut().find(|vertex| vertex.id == start) else {
            // the hole lies outside this record (a child that does not
            // contain it)
            continue;
        };
        // the anchor's own inner half-edge is always in the record (its edge
        // is a remaining edge or a crown leg there); the far end's is when
        // the record holds the chain's far end (else its slot is a
        // spectator port, below)
        if !record_hedges.contains(&chain.inner[0]) {
            return Err(chain_error(format!(
                "the anchor's inner half-edge {} of the hole at {start} is not a half-edge of the record",
                chain.inner[0]
            )));
        }
        if mentions_mass(&vertex.factor, &masses) {
            return Err(Decline::new(
                "surrogate_coefficient_mass",
                format!(
                    "the insertion-surrogate coefficient at vertex {start} contains a mass symbol (the normalised insertion is a massless integral)"
                ),
            ));
        }
        let factor = std::mem::replace(&mut vertex.factor, Atom::num(1));
        // C_H carries the momentum of the excised diagram's FIRST external
        // edge (insertions.rs normalized_insertion), e_in or e_out: write it in
        // e_in, the anchor's own edge (the chain conserves the momentum: Q(e_out)
        // = s Q(e_in), s = +1 when one end enters the chain and the other
        // leaves it), so that an absent e_out never stands as a constant vector
        let (e_in, _) = chain.inner_edges[0];
        // ... and in the RECORD's name of that edge at the anchor: `e_in`, or
        // the incoming leg `PAIRED_LEG_OFFSET + e_in` when e_in is a paired
        // crown edge of the record whose sink is the anchor (beta's pick 8;
        // the same rename as `crown_legs_at_node` makes in a ct: Delta)
        let anchor_edge = anchor_leg(edges, chain.inner[0], e_in)?;
        let factor = in_anchor_edge(factor, chain.inner_edges, anchor_edge);
        vertex.rule = format!("ins:{}", vertex.rule);
        let mut hedges = chain
            .inner
            .iter()
            .copied()
            .filter(|h| record_hedges.contains(h))
            .collect::<Vec<_>>();
        hedges.sort_unstable();
        hedges.dedup();
        // a KEPT and CUT hole (the anchor in the record, the far end's inner
        // half-edge not): the single-hedge ins: with the far-end slot as a
        // spectator port (WH-RECORD-FORMAT v1.5 proposal (a)), or, on a reader
        // that refuses that form, the designed decline `surrogate_hole_cut`
        // ([`hole_cut`])
        if !record_hedges.contains(&chain.inner[1]) && hole_cut() == HoleCut::Decline {
            return Err(Decline::new(
                "surrogate_hole_cut",
                format!(
                    "the record holds the anchor {start} of a {}-loop hole but not its far end (inner half-edge {}): the single-hedge ins: form is off ({})",
                    chain.vertices.len() - 1,
                    chain.inner[1],
                    HOLE_CUT_VARIABLE
                ),
            ));
        }
        insertions.push(Insertion {
            vertex: start,
            hedges,
            cycles: chain.cycles,
            coefficient: Atom::num(1),
            factor,
        });
        free.insert(start);
    }
    // ct: children: their own support's cycles plus the holes they contract
    for (node, index) in nodes {
        let support =
            super::super::plan::uv_child_edge_ids(graph, &region.insertions[*index].support)
                .into_iter()
                .collect::<BTreeSet<_>>();
        let mut cycles = (
            closed_cycles(&merged, &support, LineKind::Fermion),
            closed_cycles(&merged, &support, LineKind::Ghost),
        );
        for chain in &chains {
            if merge.get(&chain.vertices[0]) == Some(node) {
                if free.contains(&chain.vertices[0]) {
                    return Err(Decline::new(
                        "surrogate_hole_twice",
                        format!(
                            "the hole at {} is both a free ins: vertex and inside the ct: child at {node}",
                            chain.vertices[0]
                        ),
                    ));
                }
                cycles.0 += chain.cycles.0;
                cycles.1 += chain.cycles.1;
            }
        }
        if let Some(insertion) = insertions
            .iter_mut()
            .find(|insertion| insertion.vertex == *node)
        {
            insertion.cycles = cycles;
        }
    }
    // correction 4: every insertion's coefficient and factor depends on eps
    // only through sums, products and integer powers; a direct record's
    // coefficient vertex is regular (its entry raises no depth)
    for insertion in insertions.iter() {
        for (what, atom) in [
            ("coefficient", &insertion.coefficient),
            ("factor", &insertion.factor),
        ] {
            match eps_dependence(atom.as_view(), GS.dim_epsilon) {
                Err(why) => {
                    return Err(Decline::new(
                        "surrogate_insertion_eps",
                        format!("insertion at {} {what}: {why}", insertion.vertex),
                    ));
                }
                Ok(false) if direct && free.contains(&insertion.vertex) => {
                    return Err(Decline::new(
                        "surrogate_insertion_eps",
                        format!(
                            "the direct record's coefficient at {} has an eps pole (its entry raises no depth)",
                            insertion.vertex
                        ),
                    ));
                }
                Ok(_) => {}
            }
        }
    }
    // A hole cut by a child boundary (a chain of two or more loops): a slot
    // of the coefficient, or of a chain end's edge, whose partner lies outside
    // the record is a SPECTATOR: free in the record's product and of a
    // half-edge the record does not have (a port, as vakint admits it), or
    // of a half-edge the record has without a factor that contracts it
    // (a port too). A free slot of an INSERTION that is a half-edge of the
    // record but not one of its own (`hedges`) cannot be stated (a ct:
    // child that contracts a hole's anchor but not its far end): declined
    // by design, `surrogate_split_chain`, before the record exists.
    let legs = edges
        .iter()
        .filter(|edge| edge.is_external())
        .flat_map(|edge| [edge.source_hedge, edge.target_hedge])
        .flatten()
        .collect::<BTreeSet<_>>();
    let mut atoms = Vec::new();
    for edge in edges.iter() {
        atoms.push((None, &edge.factor));
    }
    for vertex in vertices.iter() {
        atoms.push((None, &vertex.factor));
    }
    for (index, insertion) in insertions.iter().enumerate() {
        atoms.push((Some(index), &insertion.factor));
    }
    let mut owners = BTreeMap::<usize, BTreeSet<usize>>::new();
    let per_atom = atoms
        .iter()
        .map(|(_, atom)| hedge_slots(atom))
        .collect::<Vec<_>>();
    let mut bispinor = BTreeSet::new();
    for (n, slots) in per_atom.iter().enumerate() {
        for (h, bis) in slots {
            owners.entry(*h).or_default().insert(n);
            if *bis {
                bispinor.insert(*h);
            }
        }
    }
    // The reader's locality rule for an insertion (factorized.rs
    // check_insertions, b6bb69835 and ccba5e0c8a alike): every free slot of
    // its factor is a slot of one of its `hedges`, or an open port of a
    // half-edge the record does not have. Checked here for every insertion
    // (crown legs included), so that a record the reader would refuse is a
    // designed decline before it exists, never an IndexNotLocal error.
    let inner = chains
        .iter()
        .flat_map(|chain| chain.inner)
        .collect::<BTreeSet<_>>();
    let mut spectators = BTreeSet::new();
    for insertion in insertions.iter() {
        let contracted = vertices
            .iter()
            .any(|v| v.id == insertion.vertex && v.rule.starts_with("ct:"));
        for (h, bis) in hedge_slots(&insertion.factor) {
            if insertion.hedges.contains(&h) {
                continue;
            }
            if record_hedges.contains(&h) {
                // a ct: that contracts a hole's anchor while the record keeps
                // its far end: C_H's far slot is a half-edge of the record
                // that is not the ct: node's (beta's checked rule, point 4)
                return Err(Decline::new(
                    "surrogate_split_chain",
                    format!(
                        "the insertion at {} carries the slot of half-edge {h}, a half-edge of the record that is not one of its own (a child that contracts a hole's anchor but not its far end)",
                        insertion.vertex
                    ),
                ));
            }
            let holders = owners.get(&h).map_or(0, BTreeSet::len);
            if holders > 1 {
                return Err(Decline::new(
                    "surrogate_split_chain",
                    format!(
                        "the slot of half-edge {h} (not in the record) is shared by the insertion at {} and another factor (the two ends of one hole in different contracted children)",
                        insertion.vertex
                    ),
                ));
            }
            if contracted && bis && inner.contains(&h) {
                // a contracted cut QUARK hole: an open spinor port no
                // projector of the record closes (point 4)
                return Err(Decline::new(
                    "surrogate_split_chain",
                    format!(
                        "the ct: child at {} carries the far-end bispinor slot of half-edge {h} of a hole it cuts (a contracted cut quark hole)",
                        insertion.vertex
                    ),
                ));
            }
            spectators.insert(h);
        }
    }
    // a free slot of an edge or vertex factor whose partner lies outside the
    // record (a chain end's edge without the anchor): a port
    for (h, holders) in &owners {
        if legs.contains(h) || holders.len() != 1 {
            continue;
        }
        let holder = *holders.iter().next().unwrap();
        if atoms[holder].0.is_none() {
            spectators.insert(*h);
        }
    }
    // U7c delta pre-review N4: a spectator port is only ever a hole's inner
    // half-edge (any other open slot is a builder error, never a port), and
    // the ports are exactly those the record's STRUCTURE implies. Per hole
    // whose coefficient C_H has a slot of the far end's inner half-edge
    // h_out, that slot pairs C_H with the far edge e_out's factor; it is a
    // port when exactly one side is in the record:
    // - C_H's side (owner in the record: the kept anchor, or the ct: node
    //   holding it) with h_out absent from the record, e_out neither inside
    //   the owner's child (the pair is then closed inside that Delta) nor
    //   inside another child;
    // - e_out's side with C_H's owner absent: e_out a remaining (non-leg)
    //   edge of the record, or inside a contracted child of it.
    if let Some(h) = spectators.iter().find(|h| !inner.contains(h)) {
        return Err(Decline::new(
            "surrogate_stray_slot",
            format!(
                "the slot of half-edge {h} is open in the record and is not a hole's inner half-edge"
            ),
        ));
    }
    let supports = nodes
        .iter()
        .map(|(node, index)| {
            (
                *node,
                super::super::plan::uv_child_edge_ids(graph, &region.insertions[*index].support)
                    .into_iter()
                    .collect::<BTreeSet<_>>(),
            )
        })
        .collect::<Vec<_>>();
    let expected = chains
        .iter()
        .filter(|chain| {
            let (anchor, h_out, e_out) =
                (chain.vertices[0], chain.inner[1], chain.inner_edges[1].0);
            let owner_child = merge.get(&anchor).copied();
            let owner = free.contains(&anchor) || owner_child.is_some();
            let in_child = |node: Option<usize>| {
                supports.iter().any(|(n, support)| {
                    node.is_none_or(|node| *n == node) && support.contains(&e_out)
                })
            };
            let in_owner_child = owner_child.is_some() && in_child(owner_child);
            let in_other_child = supports
                .iter()
                .any(|(n, support)| Some(*n) != owner_child && support.contains(&e_out));
            let remaining = edges
                .iter()
                .any(|edge| edge.id == e_out && !edge.is_external());
            chain.far_slot
                && if owner {
                    !record_hedges.contains(&h_out) && !in_owner_child && !in_other_child
                } else {
                    (record_hedges.contains(&h_out) && remaining) || in_child(None)
                }
        })
        .map(|chain| chain.inner[1])
        .collect::<BTreeSet<_>>();
    if spectators != expected {
        return Err(Decline::new(
            "surrogate_ports",
            format!(
                "the record's spectator ports {spectators:?} are not those its holes imply {expected:?}"
            ),
        ));
    }
    vertices.sort_by_key(|vertex| vertex.id);
    insertions.sort_by_key(|insertion| insertion.vertex);
    Ok(Finished {
        original_cycles: provenance.original_cycles,
        spectators,
    })
}

/// What [`finish`] gives the record builder: the record's
/// `original_cycles` and its spectator ports (half-edges).
pub(super) struct Finished {
    pub original_cycles: (u32, u32),
    pub spectators: BTreeSet<usize>,
}

/// The half-edges `h` of the slots `mink(.., hedge(h, k))` / `bis(..,
/// hedge(h, k))` in `atom`.
fn hedge_slots(atom: &Atom) -> BTreeMap<usize, bool> {
    let mink = symbolica::symbol!("spenso::mink");
    let bis = symbolica::symbol!("spenso::bis");
    let mut found = BTreeMap::new();
    atom.as_view().visitor(&mut |view| {
        if let AtomView::Fun(function) = view
            && (function.get_symbol() == mink || function.get_symbol() == bis)
            && function.get_nargs() == 2
            && let Some(AtomView::Fun(index)) = function.iter().nth(1)
            && index.get_symbol() == GS.hedgeaind
            && let Some(h) = index.iter().next()
            && let Ok(h) = usize::try_from(h)
        {
            *found.entry(h).or_insert(false) |= function.get_symbol() == bis;
            return false;
        }
        true
    });
    found
}

// ---------------------------------------------------------------------------
// The direct record (<= 4 loops)
// ---------------------------------------------------------------------------

/// What [`build_direct`] built.
struct BuiltDirect {
    record: FactorizedNumerator,
    text: String,
    insertions: usize,
}

/// The `kind direct` record of a surrogate at <= 4 loops: the graph of
/// `integrate_propagator_with`, routed as there (its loop basis, p(1) the
/// first external leg, the other leg ±p(1)).
fn build_direct(
    provenance: &Provenance,
    graph: &Graph,
    model: &Model,
    uv: &UVgenerationSettings,
    loops: usize,
) -> std::result::Result<BuiltDirect, Decline> {
    if uv.vakint.additional_normalization.trim() != "1" {
        return Err(Decline::new(
            "additional_normalization",
            uv.vakint.additional_normalization.clone(),
        ));
    }
    let routing_error = |detail: String| Decline::new("routing", detail);
    let boundary = graph.external_filter::<SuBitGraph>();
    let (external, routing) = RStarIntegrand::physical_external_routing(graph, &boundary)
        .map_err(|error| routing_error(format!("{error:#}")))?;
    let support = graph.no_dummy().subtract(&boundary);
    if graph.n_loops(&support) != loops {
        return Err(routing_error(format!(
            "{} loops for a {loops}-loop surrogate",
            graph.n_loops(&support)
        )));
    }
    let marked = graph
        .iter_edges_of(&support)
        .next()
        .map(|(_, edge, _)| edge)
        .ok_or_else(|| routing_error("no internal edge".into()))?;
    let region = GraphRegion {
        remaining: support.clone(),
        support,
        insertions: Vec::new(),
        boundary,
        loops,
    };
    let integrand =
        RStarIntegrand::new_with_closure(graph, model, &region, &[], marked, Atom::Zero, None)
            .map_err(|error| routing_error(format!("{error:#}")))?;
    let basis = integrand.basis();
    let loop_edges = basis
        .loop_edges
        .iter()
        .map(|edge| usize::from(*edge))
        .collect::<Vec<_>>();
    if loop_edges.len() != loops {
        return Err(routing_error(format!(
            "{} loop edges for {loops} loops",
            loop_edges.len()
        )));
    }
    let loop_index = loop_edges
        .iter()
        .enumerate()
        .map(|(i, edge)| (*edge, i))
        .collect::<BTreeMap<_, _>>();
    let x0 = usize::from(external);
    let p = function!(GS.emr_mom, x0);
    let sign_of = |edge: usize| -> std::result::Result<i64, Decline> {
        let image = function!(GS.emr_mom, edge).replace_multiple(&routing);
        if image == p {
            Ok(1)
        } else if image == -p.clone() {
            Ok(-1)
        } else {
            Err(routing_error(format!(
                "the external routing maps Q({edge}) to {image}, not ±Q({x0})"
            )))
        }
    };
    let width = loops + 1;
    let vector = |edge: EdgeIndex| -> std::result::Result<Vec<i64>, Decline> {
        let mut out = vec![0i64; width];
        let parse = |part: &Atom| {
            q_coefficients(part).ok_or_else(|| {
                routing_error(format!("edge {}: momentum {part}", usize::from(edge)))
            })
        };
        for (e, c) in parse(&basis.loop_atom::<Atom>(edge, GS.emr_mom, &[], true))? {
            let Some(i) = loop_index.get(&e) else {
                return Err(routing_error(format!(
                    "edge {}: Q({e}) is not a loop edge",
                    usize::from(edge)
                )));
            };
            out[*i] += c;
        }
        for (e, c) in parse(&basis.ext_atom::<Atom>(edge, GS.emr_mom, &[], true))? {
            out[loops] += c * sign_of(e)?;
        }
        Ok(out)
    };
    let table = edge_table(graph, model);
    let mut edges = Vec::new();
    let mut legs = BTreeSet::new();
    let mut node_index = BTreeMap::new();
    for (pair, edge, _) in graph.iter_edges_of(&region.boundary) {
        let HedgePair::Unpaired { hedge, flow } = pair else {
            return Err(routing_error(format!(
                "external edge {} is not unpaired",
                usize::from(edge)
            )));
        };
        let id = usize::from(edge);
        let info = &table[&id];
        let node = graph.node_id(hedge);
        node_index.insert(usize::from(node), node);
        let v = usize::from(node);
        let h = usize::from(hedge);
        let mut momentum = vec![0i64; width];
        momentum[loops] = sign_of(id)?;
        let (source, target, source_hedge, target_hedge) = if flow == Flow::Sink {
            (Endpoint::External, Endpoint::Vertex(v), None, Some(h))
        } else {
            (Endpoint::Vertex(v), Endpoint::External, Some(h), None)
        };
        edges.push(FactorizedEdge {
            id,
            source,
            target,
            source_hedge,
            target_hedge,
            particle: info.particle.clone(),
            line: info.line,
            power: 1,
            mass_squared: Atom::Zero,
            momentum,
            factor: Atom::num(1),
        });
        legs.insert(h);
    }
    for (pair, edge, data) in graph.iter_edges_of(&region.remaining) {
        let HedgePair::Paired { source, sink } = pair else {
            continue;
        };
        let id = usize::from(edge);
        let info = &table[&id];
        let (s, t) = (graph.node_id(source), graph.node_id(sink));
        node_index.insert(usize::from(s), s);
        node_index.insert(usize::from(t), t);
        edges.push(FactorizedEdge {
            id,
            source: Endpoint::Vertex(usize::from(s)),
            target: Endpoint::Vertex(usize::from(t)),
            source_hedge: Some(usize::from(source)),
            target_hedge: Some(usize::from(sink)),
            particle: info.particle.clone(),
            line: info.line,
            power: 1,
            mass_squared: Atom::Zero,
            momentum: vector(edge)?,
            factor: data.data.num.value.map_mink_dim(GS.dim),
        });
    }
    edges.sort_by_key(|edge| edge.id);
    let mut vertices = node_index
        .iter()
        .map(|(raw, node)| {
            let vertex = &graph[*node];
            FactorizedVertex {
                id: *raw,
                rule: vertex.vertex_rule.map_or_else(
                    || vertex.name.value.clone(),
                    |rule| rule.resolve(model).name.clone(),
                ),
                factor: vertex.get_num().map_mink_dim(GS.dim),
            }
        })
        .collect::<Vec<_>>();
    let mut insertions = Vec::new();
    let finished = finish_with(
        provenance,
        graph,
        &table,
        &region,
        &BTreeMap::new(),
        &[],
        &mut edges,
        &mut vertices,
        &mut insertions,
        true,
    )?;
    let factor_atoms = edges
        .iter()
        .map(|edge| &edge.factor)
        .chain(vertices.iter().map(|vertex| &vertex.factor))
        .chain(insertions.iter().map(|insertion| &insertion.factor))
        .collect::<Vec<_>>();
    let mut port_hedges = legs.clone();
    port_hedges.extend(finished.spectators.iter().copied());
    let ports = slots_of(&factor_atoms, &port_hedges);
    let original_cycles = finished.original_cycles;
    let n_insertions = insertions.len();
    let record = FactorizedNumerator {
        label: format!(
            "u7c kind=surrogate graph={} original={} scope=direct loops={loops} holes={n_insertions} id=0",
            graph.name, provenance.original
        )
        .replace(['\t', '\n', '\r'], " "),
        kind: FactorizedKind::Direct,
        dimension: GS.dim,
        loop_edges: loop_edges.iter().map(|edge| Some(*edge)).collect(),
        external_edges: vec![Some(x0)],
        edges,
        vertices,
        ports,
        projector: None,
        prefactor: Prefactor {
            overall: Atom::num(1),
            loop_sign: 1,
            nf_power: 0,
            scalars: Atom::num(1),
            colour: Atom::num(1),
        },
        colour_branch: None,
        insertions,
        record: Some(DirectRecord {
            original_cycles: Some(original_cycles),
            depth: Some(loops),
            weights: Vec::new(),
        }),
    };
    let text = record
        .to_text()
        .map_err(|error| Decline::new("text", error.to_string()))?;
    let record = FactorizedNumerator::from_text(&text).map_err(|error| Decline::vakint(&error))?;
    // beta's guard (pick 8), as for every K1 / K2 record: a direct record has
    // no Taylor step, where vakint reads an undeclared Q(x) as a constant
    // vector, so `undeclared_momentum` is a real check here
    super::declared_momenta(&record)?;
    Ok(BuiltDirect {
        record,
        text,
        insertions: n_insertions,
    })
}

/// Our attempt at the direct record (no state changes). `Err`: a hard error
/// (unknown masters, an exhausted budget), which fails the diagram.
fn attempt_direct(
    provenance: &Provenance,
    graph: &Graph,
    model: &Model,
    vakint: &Vakint,
    uv: &UVgenerationSettings,
    loops: usize,
) -> Result<Attempt> {
    let started = Instant::now();
    let kind = Kind::Surrogate;
    let mut attempt = Attempt {
        kind,
        graph: graph.name.clone(),
        id: 0,
        loops,
        depth: loops,
        order: valid_order(loops, loops),
        closure: "none",
        // the direct record's own record class for the ledger's setting-wide
        // rule (U7a follow-ups, item B)
        shape: "direct",
        text: None,
        outcome: Outcome::Error(Decline::new("unbuilt", "")),
        fields: String::new(),
        build_seconds: 0.0,
        evaluate_seconds: 0.0,
        finish_seconds: 0.0,
        delta: None,
        conserved: 0,
    };
    let built = if super::super::pint6_dump::settings().is_some() {
        Err(Decline::new(
            "pint6_dump",
            "PINT6 dumps keep today's records",
        ))
    } else if loops == 0 {
        Err(Decline::new("zero_loops", "a zero-loop surrogate"))
    } else {
        build_direct(provenance, graph, model, uv, loops)
    };
    let built = match built {
        Ok(built) => built,
        Err(decline) => {
            attempt.outcome = Outcome::refused(decline);
            attempt.build_seconds = started.elapsed().as_secs_f64();
            return Ok(attempt);
        }
    };
    let build_seconds = started.elapsed().as_secs_f64();
    attempt.build_seconds = build_seconds;
    attempt.text = Some(built.text.clone());
    let mut settings = uv.vakint.true_settings();
    settings.number_of_terms_in_epsilon_expansion = i64::try_from(loops)?;
    let prefactor = built.record.prefactor.total();
    let evaluate_started = Instant::now();
    let result = evaluate_factorized_ports(vakint, &settings, &built.record, &prefactor);
    attempt.evaluate_seconds = evaluate_started.elapsed().as_secs_f64();
    let value = match result {
        Ok(value) => value,
        Err(error) => {
            let text = error.to_string();
            if is_budget_text(&text) {
                return Err(eyre!(
                    "resource budget exhausted in the uncontracted record kind=surrogate graph={} id=0: {text}",
                    attempt.graph
                ));
            }
            if is_unknown_master_text(&text) {
                tally_unknown_masters(kind);
                eprintln!(
                    "UNCONTRACTED kind=surrogate graph={} id=0 outcome=error:unknown_masters detail=\"{}\"",
                    attempt.graph,
                    super::escape(&text)
                );
                return Err(unknown_masters_failure(kind, &attempt.graph, 0, &text));
            }
            attempt.outcome = Outcome::refused(Decline::vakint(&error));
            return Ok(attempt);
        }
    };
    attempt.delta = Some(value.timing);
    let finish_started = Instant::now();
    let finished = match RStarIntegrand::finish_tensor(value.value.clone(), 0) {
        Ok(finished) => finished,
        Err(error) => {
            attempt.outcome = Outcome::refused(Decline::new("finish_tensor", format!("{error:#}")));
            return Ok(attempt);
        }
    };
    attempt.finish_seconds = finish_started.elapsed().as_secs_f64();
    let rss = vakint::timing::Stage::ALL
        .iter()
        .filter_map(|stage| value.timing.stage(*stage).peak_rss_bytes)
        .max()
        .map_or_else(
            || "-".to_owned(),
            |bytes| format!("{:.2}", bytes as f64 / (1u64 << 30) as f64),
        );
    attempt.fields = format!(
        " depth={loops} degree=0 externals=1 insertions={} conserved=0 today=direct patterns=- blocks=- key_plans=- lowered_entries=- build_s={build_seconds:.3} finish_s={:.3}{}{} wall_s={:.3} rss_gib={rss} record_bytes={}",
        built.insertions,
        attempt.finish_seconds,
        value.timing.fields(),
        value.timing.stage_fields(),
        started.elapsed().as_secs_f64(),
        built.text.len()
    );
    attempt.outcome = Outcome::Value(finished);
    Ok(attempt)
}

/// The direct integral of a surrogate at <= 4 loops (`graph`, `loops`):
/// `today` is `integrate_propagator(surrogate)`. Off: `today`, nothing else
/// (with `GAMMALOOP_UNCONTRACTED_TIMING=1` also its occurrence's timing
/// line). On / check with `surrogate` enabled: one record (a site of the
/// diagram's ledger; outside one an error); on: vakint's value, or today's
/// for a designed decline; check: ours first, then today's, compared
/// through `ε^-1`, today's returned. An error outcome fails the diagram.
#[allow(clippy::too_many_arguments)]
pub(crate) fn direct(
    provenance: Option<&Arc<Provenance>>,
    graph: &Graph,
    model: &Model,
    vakint: &Vakint,
    uv: &UVgenerationSettings,
    loops: usize,
    today: impl FnOnce() -> Result<Atom>,
) -> Result<Atom> {
    let kind = Kind::Surrogate;
    let switch = switch();
    let term = || format!("direct{{loops={loops}}}");
    let finish_hook = |hook: Option<uncontracted_ledger::Hook>, route: &str| {
        if let Some(hook) = hook {
            hook.finish(route);
        }
    };
    if switch.mode == Mode::Off {
        if !uncontracted_ledger::timing_enabled() {
            return today();
        }
        let hook = uncontracted_ledger::hook(graph, kind, 0, false, term);
        let value = today();
        finish_hook(hook, if value.is_ok() { "today" } else { "failed" });
        return value;
    }
    let enabled = switch.enabled(kind);
    if !enabled {
        let hook = uncontracted_ledger::hook(graph, kind, 0, false, term);
        uncontracted_ledger::inert(kind);
        let value = today();
        finish_hook(hook, if value.is_ok() { "today" } else { "failed" });
        return value;
    }
    let provenance = provenance.ok_or_else(|| missing_context(&graph.name))?;
    settings_refusal(&graph.name, None, uv)?;
    uncontracted_ledger::k5_site(kind, &graph.name)?;
    let hook = uncontracted_ledger::hook(graph, kind, 0, true, term);
    uncontracted_ledger::expect(kind);
    let started = Instant::now();
    let attempt = with_context(provenance.clone(), || {
        attempt_direct(provenance, graph, model, vakint, uv, loops)
    })?;
    #[cfg(test)]
    let attempt = uncontracted_ledger::inject(attempt);
    let uncontracted_elapsed = started.elapsed();
    if switch.mode == Mode::Check {
        if attempt.is_error() {
            report(&attempt, None);
            finish_hook(hook, &attempt.outcome_name());
            return Err(super::super::ir_dirac::diagnostic_failure(error_failure(
                &attempt,
                switch.mode,
            )));
        }
        let today_started = Instant::now();
        let today = match today() {
            Ok(today) => today,
            Err(error) => {
                report(&attempt, None);
                finish_hook(hook, "failed");
                return Err(error);
            }
        };
        let check = attempt.value().map(|ours| {
            compare(
                &today,
                &ours,
                GS.dim_epsilon,
                attempt.order,
                today_started.elapsed().as_secs_f64(),
                uncontracted_elapsed.as_secs_f64(),
            )
        });
        report(&attempt, check.as_ref());
        finish_hook(
            hook,
            &if check.is_some() {
                "check".to_owned()
            } else {
                attempt.outcome_name()
            },
        );
        if let Some(check) = &check
            && !check.equal
        {
            return Err(super::super::ir_dirac::diagnostic_failure(check_failure(
                &attempt,
                check,
                GS.dim_epsilon,
            )));
        }
        return Ok(today);
    }
    report(&attempt, None);
    if attempt.is_error() {
        finish_hook(hook, &attempt.outcome_name());
        return Err(super::super::ir_dirac::diagnostic_failure(error_failure(
            &attempt,
            switch.mode,
        )));
    }
    match attempt.value() {
        Some(value) => {
            finish_hook(hook, &attempt.outcome_name());
            Ok(value)
        }
        None => {
            let value = today();
            finish_hook(hook, &attempt.outcome_name());
            value
        }
    }
}

#[cfg(test)]
thread_local! {
    static SETTINGS_CHECK: std::cell::Cell<bool> = const { std::cell::Cell::new(true) };
}

/// Tests: run `f` WITHOUT the K5 site's settings check (to show what the
/// ledger's generic rules alone catch).
#[cfg(test)]
pub(crate) fn without_settings_check<T>(f: impl FnOnce() -> T) -> T {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            let previous = self.0;
            SETTINGS_CHECK.with(|slot| slot.set(previous));
        }
    }
    let _restore = Restore(SETTINGS_CHECK.with(|slot| slot.replace(false)));
    f()
}

/// Tests: run `f` with `provenance` as this thread's surrogate context.
#[cfg(test)]
pub(crate) fn with_test_context<T>(provenance: Arc<Provenance>, f: impl FnOnce() -> T) -> T {
    with_context(provenance, f)
}
