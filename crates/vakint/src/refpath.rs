//! The reference path R2 of the numerator assignment (helper VH, lane
//! `refpath`): one FeynKit DOT graph of the background-field `B -> B`
//! self-energy, from its per-edge / per-vertex spenso numerators to the
//! Laurent expansion of its contribution, through today's production pieces
//! only, with the stage times and the gate view of what RustRed receives.
//!
//! 1. parse: `Model::from_json`, `FeynmanDiagram::from_dot`;
//! 2. colour: every local tensor is split into its colour rows
//!    (`coefficient_list` over the colour representations, as the beta
//!    thread's `color.rs` `split_local`); the Cartesian product of the rows
//!    gives the colour branches (a quartic vertex has three); each branch's
//!    closed colour network (rows × δ of the two external adjoint ports) is
//!    reduced to SU(N) invariants with idenso's
//!    `try_reduce_closed_color_network` and divided by `N_A`; branches with a
//!    zero coefficient are dropped;
//! 3. contract: per branch, the product of the Lorentz/Dirac residuals is
//!    routed onto FeynKit's loop-momentum basis, continued to `D` dimensions
//!    (`mink(4)` → `mink(D)`), cooked as gammaloop's R* simplifier cooks
//!    (`CookSettings::indices()` with `ReversibleEncoding` and representation
//!    payloads), closed with the projector on its two dangling Lorentz ports,
//!    and simplified by idenso's production chain (`infer_indexed` →
//!    `simplify_algebra{gamma: canonical(), epsilon: true}` → `to_dots`);
//!    couplings are expanded (`gs = 1`), and the beta thread's normalisation
//!    `1/(i q²)` (metric and transverse also `1/(D - 1)`) is applied;
//! 4. lower, reduce, finalise: the colour-weighted sum of the branches, times
//!    `overall_factor` (symmetry factor, fermion and ghost loop signs) and
//!    `n_f` per closed quark loop, enters vakint through its FeynKit ingress
//!    (`VakintExpression::from_diagram`: edge momenta → `k(i)`, `±p(1)`) and
//!    is evaluated by `Vakint::evaluate` with RustRed alone (the ordinary p
//!    path R2: direct lowering of plain requests, `fold_requests`, one
//!    RustRed reduction per folded channel, the p finalizer), observed by
//!    [`crate::rustred_evaluation::p_gate`].
//!
//! The normalisation dictionary to the renormalisation ledger's `Π_l` is in
//! `N/notes/REFPATH.md` (helper VH).

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use std::time::Instant;

use eyre::{Result, ensure, eyre};
use feynkit_graph::{FeynmanDiagram, NumeratorInterface};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use idenso::color::CS;
use idenso::dirac::GammaSimplifySettings;
use idenso::representations::{ColorAdjoint, ColorFundamental, ColorSextet};
use idenso::tensor::{AlgebraSettings, SymbolicTensor};
use idenso::{CookMode, CookSettings, IndexTooling, W_};
use spenso::network::library::symbolic::ETS;
use spenso::shadowing::TensorCollectFilter;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::partial::PartialStructure;
use spenso::structure::representation::{LibrarySlot, Minkowski, RepName};
use spenso::structure::slot::IsAbstractSlot;
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use symbolica::symbol;

use crate::diagram_integral::DiagramIntegralOptions;
use crate::rustred_evaluation::p_gate::{self, GateObservation, GateView, ZeroSectorCount};
use crate::utils::vakint_macros::vk_symbol;
use crate::{Vakint, VakintExpression, VakintSettings};

type Tensor = SymbolicTensor<PartialStructure>;

/// The symbols of the reference path's output.
pub struct RefpathSymbols {
    /// The space-time dimension of the contraction (also the gate view's `D`).
    pub d: Symbol,
    /// SU(N) invariants: `C_A`, `C_F`, `T_F` (= I2R), `N_A`, `N_R`.
    pub ca: Symbol,
    pub cf: Symbol,
    pub tf: Symbol,
    pub na: Symbol,
    pub nr: Symbol,
    /// Number of quark flavours (one power per closed quark loop).
    pub nf: Symbol,
    /// Projection tags of the combined contraction.
    pub metric_tag: Symbol,
    pub qq_tag: Symbol,
}

pub static SYMBOLS: std::sync::LazyLock<RefpathSymbols> =
    std::sync::LazyLock::new(|| RefpathSymbols {
        d: p_gate::gate_dimension(),
        ca: symbol!("refpath::CA"),
        cf: symbol!("refpath::CF"),
        tf: symbol!("refpath::TF"),
        na: symbol!("refpath::NA"),
        nr: symbol!("refpath::NR"),
        nf: symbol!("refpath::nf"),
        metric_tag: symbol!("refpath::tag_metric"),
        qq_tag: symbol!("refpath::tag_qq"),
    });

/// Which projection closes the two external Lorentz ports.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Projection {
    /// `g^{μν}/(D-1)`: one contraction; equals the transverse projection on
    /// the sum of the graphs of one loop order (transversality).
    Metric,
    /// `(g^{μν} - q^μ q^ν/q²)/(D-1)`.
    Transverse,
    /// `q^μ q^ν/q²` (must vanish on the sum of the graphs).
    Longitudinal,
    /// Transverse and longitudinal from one contraction (two evaluations).
    Both,
}

impl std::str::FromStr for Projection {
    type Err = String;

    fn from_str(text: &str) -> std::result::Result<Self, String> {
        match text {
            "metric" => Ok(Self::Metric),
            "transverse" => Ok(Self::Transverse),
            "longitudinal" => Ok(Self::Longitudinal),
            "both" => Ok(Self::Both),
            _ => Err(format!("unknown projection {text}")),
        }
    }
}

impl Projection {
    pub fn name(self) -> &'static str {
        match self {
            Self::Metric => "metric",
            Self::Transverse => "transverse",
            Self::Longitudinal => "longitudinal",
            Self::Both => "both",
        }
    }
}

/// One colour branch: its colour coefficient (SU(N) invariants, already
/// divided by `N_A`) and the Lorentz/Dirac residuals of every local, as
/// authored atoms.
pub struct ColourBranch {
    pub coefficient: Atom,
    pub residuals: Vec<Atom>,
}

/// Seconds of the stages of one graph (wall time on this thread).
#[derive(Clone, Debug, Default)]
pub struct StageSeconds {
    pub parse: f64,
    pub colour: f64,
    pub contract: f64,
    /// vakint: FeynKit ingress (`from_diagram`) + front end (canonicalise,
    /// match) + direct reading and lowering + fold: the evaluation's wall
    /// time minus `reduce`, `finalise` and the observer's copies.
    pub lower: f64,
    /// vakint's per-key Q(d) coefficients (`prepare`) + the RustRed library.
    pub reduce: f64,
    /// Rendering of the reductions + the p finalizer.
    pub finalise: f64,
    /// Of `lower`: the FeynKit ingress alone.
    pub ingress: f64,
    /// Of `lower`: vakint's split/read and lowering + fold (its statistics).
    pub split: f64,
    pub lowering: f64,
    /// Of `reduce`.
    pub prepare: f64,
    pub library: f64,
    /// Of `finalise`.
    pub render: f64,
    pub finalize: f64,
    /// Term count of the expanded vakint numerator (statistics only, not
    /// part of any stage above).
    pub count: f64,
}

/// The evaluated contribution of one projection.
#[derive(Clone, Debug)]
pub struct ProjectedResult {
    pub projection: &'static str,
    /// The Laurent expansion (ε symbol of the settings), colour invariants
    /// and `n_f` symbolic, `gs = 1`.
    pub value: Atom,
    /// What RustRed received, every key.
    pub observation: GateObservation,
}

/// Everything recorded for one graph.
pub struct GraphRun {
    pub name: String,
    pub loops: usize,
    pub quark_loops: usize,
    pub overall_factor: Atom,
    /// Colour branches before and after dropping zero colour coefficients.
    pub branches_total: usize,
    pub branches: usize,
    /// Why the graph is an exact zero without evaluation, if it is.
    pub zero_reason: Option<String>,
    pub seconds: StageSeconds,
    /// Terms of the expanded numerator vakint reads (per projection, summed).
    pub contracted_terms: usize,
    /// The numerator handed to vakint (per projection), unexpanded.
    pub numerators: Vec<(&'static str, Atom)>,
    pub results: Vec<ProjectedResult>,
}

/// The reference path's settings for vakint: RustRed alone, `Symbolic` master
/// mode, the G-scheme, `terms` Laurent terms, one core.
pub fn vakint_settings(terms: i64) -> VakintSettings {
    VakintSettings {
        evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
            crate::RustRedEvaluationOptions {
                master_mode: crate::RustRedMasterMode::Symbolic,
            },
        )]),
        form_exe_path: "/refpath-must-not-call-form".into(),
        integral_normalization_factor: crate::LoopNormalizationFactor::GScheme,
        number_of_terms_in_epsilon_expansion: terms,
        run_time_decimal_precision: 40,
        use_dot_product_notation: true,
        n_cores: 1,
        ..VakintSettings::default()
    }
}

/// A new engine that shares `template`'s loaded p library (the immutable
/// `PLibrary`, loaded once) through a fresh handle: a cold request memo and
/// step cache per engine, the package's master catalogue, directory and
/// relations root kept. Repetitions on fresh engines then measure cold
/// reductions without reloading the package.
pub fn fresh_engine(template: &Vakint) -> Result<Vakint> {
    use crate::rustred_evaluation::p_data::PLibraryHandle;
    let mut engine = Vakint::new()?;
    let handle = template.p_library()?;
    let mut fresh = PLibraryHandle::new(Arc::clone(&handle.library), handle.origin)?;
    fresh.master_catalogue = handle.master_catalogue.clone();
    fresh.package_directory = handle.package_directory.clone();
    fresh.relations_root = handle.relations_root.clone();
    engine.p_library = Some(Arc::new(fresh));
    Ok(engine)
}

/// Load the model of a prepared directory (`model.json`).
pub fn load_model(directory: &std::path::Path) -> Result<Arc<Model>> {
    let text = std::fs::read_to_string(directory.join("model.json"))?;
    Ok(Arc::new(Model::from_json(&text)?))
}

fn colour_filter() -> TensorCollectFilter<3> {
    TensorCollectFilter::Reps([
        ColorAdjoint {}.into(),
        ColorFundamental {}.into(),
        ColorSextet {}.into(),
    ])
}

/// The colour rows `(colour, residual)` of one local tensor; a local without
/// colour gives one row `(1, local)`.
fn split_local(tensor: &Tensor) -> Result<Vec<(Atom, Atom)>> {
    let rows = tensor.coefficient_list(colour_filter())?;
    if rows.is_empty() {
        return Ok(vec![(Atom::num(1), NumeratorInterface::export(tensor).0)]);
    }
    let mut out = Vec::with_capacity(rows.len());
    for (colour, residual) in rows {
        let residual_rows = residual.coefficient_list(colour_filter())?;
        ensure!(
            residual_rows
                .iter()
                .all(|(selected, _)| selected.expression().is_one()),
            "a colour row kept a colour tensor in its residual"
        );
        out.push((
            NumeratorInterface::export(&colour).0,
            NumeratorInterface::export(&residual).0,
        ));
    }
    Ok(out)
}

/// `coad(8, i)` → `coad(N_A, i)`, `cof(3, i)` → `cof(N_R, i)`: symbolic SU(N).
fn symbolic_colour_dimensions(colour: &Atom) -> Atom {
    let s = &*SYMBOLS;
    colour
        .replace(ColorAdjoint {}.to_symbolic([Atom::num(8), Atom::var(W_.a___)]))
        .with(ColorAdjoint {}.to_symbolic([Atom::var(s.na), Atom::var(W_.a___)]))
        .replace(ColorFundamental {}.to_symbolic([Atom::num(3), Atom::var(W_.a___)]))
        .with(ColorFundamental {}.to_symbolic([Atom::var(s.nr), Atom::var(W_.a___)]))
}

/// The adjoint slots `coad(N, i)` of a product of colour factors that occur
/// once (the external ports).
fn open_adjoint_slots(product: &Atom) -> Vec<Atom> {
    let adjoint = ColorAdjoint {}.to_symbolic([Atom::var(W_.a_), Atom::var(W_.b_)]);
    let symbol = match adjoint.as_view() {
        AtomView::Fun(f) => f.get_symbol(),
        _ => unreachable!(),
    };
    let mut counts = BTreeMap::<Atom, usize>::new();
    product.visitor(&mut |part| {
        if let AtomView::Fun(f) = part
            && f.get_symbol() == symbol
            && f.get_nargs() == 2
        {
            *counts.entry(part.to_owned()).or_default() += 1;
            return false;
        }
        true
    });
    counts
        .into_iter()
        .filter(|(_, count)| *count == 1)
        .map(|(slot, _)| slot)
        .collect()
}

/// The name of a higher invariant `name(args)` (idenso's color.h names, `gram`)
/// as a plain symbol `refpath::<name>_<args>`, an argument `coad(N_A)` written
/// `A`, `cof(N_R)` written `R`, a number as it is: e.g. `gram(3, cof, cof)` →
/// `refpath::gram_3RR`. vakint's ingress refuses `spenso::` symbols in a
/// numerator, and the gate view needs plain spectators.
fn invariant_symbol(function: AtomView<'_>) -> Option<Atom> {
    let AtomView::Fun(f) = function else {
        return None;
    };
    let symbol = f.get_symbol();
    let full = symbol.get_name();
    let name = full.rsplit("::").next().unwrap_or(full);
    let mut arguments = String::new();
    for argument in f.iter() {
        match argument {
            AtomView::Fun(rep) if rep.get_symbol() == CS.adjoint_rep => arguments.push('A'),
            AtomView::Fun(rep) if rep.get_symbol() == CS.fundamental_rep => arguments.push('R'),
            AtomView::Num(_) => arguments.push_str(&argument.to_string()),
            _ => return None,
        }
    }
    Some(Atom::var(symbol!(&format!("refpath::{name}_{arguments}"))))
}

/// idenso's invariants → `CA`, `CF`, `TF`, the higher ones → plain symbols
/// ([`invariant_symbol`]), then `N_R = T_F N_A / C_F` (the beta thread's
/// colour basis v1). `adjoint` and `fundamental` are the dimension atoms of the
/// reduced network's representations (`N_A`, `N_R` here; `spenso::dA`,
/// `spenso::Nc` in the input type's reduction).
fn map_invariants(value: &Atom, adjoint_dimension: &Atom, fundamental_dimension: &Atom) -> Atom {
    let s = &*SYMBOLS;
    let adjoint = ColorAdjoint {}.to_symbolic([adjoint_dimension.clone()]);
    let fundamental = ColorFundamental {}.to_symbolic([fundamental_dimension.clone()]);
    let named = value
        .replace(CS.cas(Atom::num(2), adjoint))
        .with(Atom::var(s.ca))
        .replace(CS.cas(Atom::num(2), fundamental.clone()))
        .with(Atom::var(s.cf))
        .replace(CS.idx(Atom::num(2), fundamental))
        .with(Atom::var(s.tf));
    let named = named.replace_map(|view, _, out| {
        if let AtomView::Fun(_) = view
            && (view.contains_symbol(CS.adjoint_rep) || view.contains_symbol(CS.fundamental_rep))
            && let Some(symbol) = invariant_symbol(view)
        {
            **out = symbol;
        }
    });
    named
        .replace(adjoint_dimension.to_pattern())
        .with(Atom::var(s.na))
        .replace(fundamental_dimension.to_pattern())
        .with(Atom::var(s.tf) * Atom::var(s.na) / Atom::var(s.cf))
}

fn named_invariants(value: &Atom) -> Atom {
    let s = &*SYMBOLS;
    map_invariants(value, &Atom::var(s.na), &Atom::var(s.nr))
}

/// The colour coefficient of a branch (lane's split): its closed network with
/// the external δ, reduced to invariants, over `N_A`.
fn colour_coefficient_lane(colour: &Atom) -> Result<Atom> {
    let colour = symbolic_colour_dimensions(colour);
    if colour.is_zero() {
        return Ok(Atom::Zero);
    }
    let ports = open_adjoint_slots(&colour);
    ensure!(
        ports.len() == 2,
        "a colour branch with {} open adjoint slots: {colour}",
        ports.len()
    );
    let closed = &colour * ETS.metric(&ports[0], &ports[1]);
    let reduced = idenso::color::network::try_reduce_closed_color_network(closed.as_view())
        .map_err(|error| eyre!("colour network: {error}"))?
        .ok_or_else(|| eyre!("not a closed colour network: {closed}"))?;
    let coefficient = (named_invariants(&reduced) / Atom::var(SYMBOLS.na))
        .expand()
        .together();
    // Higher invariants (four loops on) stay as idenso's named functions of
    // the representations; an open colour tensor cannot survive a complete
    // closed-network reduction (`Some` above).
    Ok(coefficient)
}

/// The colour branches of a diagram by the input type's colour split
/// (`vakint::p_integral::factorized::{split_local_colour, colour_coefficient}`
/// in `ColourMode::Invariants`, the helper's one implementation), the
/// coefficients written in this chain's symbols ([`map_invariants`]); zero
/// coefficients dropped; with the number of branches before the drop.
pub fn colour_branches(diagram: &FeynmanDiagram) -> Result<(Vec<ColourBranch>, usize)> {
    use crate::p_integral::factorized::{ColourMode, colour_coefficient, split_local_colour};
    let mut locals = Vec::new();
    for (_, vertex) in diagram.vertices() {
        locals.push(split_local_colour(&NumeratorInterface::export(&vertex.numerator).0)?);
    }
    for (_, _, edge) in diagram.edges() {
        locals.push(split_local_colour(&NumeratorInterface::export(&edge.numerator).0)?);
    }
    let mut partial: Vec<(Atom, Vec<Atom>)> = vec![(Atom::num(1), Vec::new())];
    for rows in &locals {
        let mut next = Vec::with_capacity(partial.len() * rows.len());
        for (colour, residuals) in &partial {
            for row in rows {
                let mut residuals = residuals.clone();
                residuals.push(row.residual.clone());
                next.push((colour * &row.colour, residuals));
            }
        }
        partial = next;
    }
    let total = partial.len();
    let adjoint = Atom::var(symbol!("spenso::dA"));
    let fundamental = Atom::var(CS.nc);
    let mut branches = Vec::new();
    for (colour, residuals) in partial {
        let coefficient = colour_coefficient(&colour, ColourMode::Invariants)?;
        let coefficient = map_invariants(&coefficient, &adjoint, &fundamental)
            .expand()
            .together();
        if std::env::var_os("REFPATH_DEBUG_COLOUR").is_some() {
            eprintln!(
                "[refpath] colour {} -> {}",
                colour.to_canonical_string(),
                coefficient.to_canonical_string()
            );
        }
        if !coefficient.is_zero() {
            branches.push(ColourBranch {
                coefficient,
                residuals,
            });
        }
    }
    Ok((branches, total))
}

/// The lane's own colour split (idenso `coefficient_list` on the cooked
/// tensors, 2026-11-18 until the switch to the input type's split): kept only
/// for [`compare_colour_splits`].
pub fn colour_branches_lane(diagram: &FeynmanDiagram) -> Result<(Vec<ColourBranch>, usize)> {
    let mut locals = Vec::new();
    for (_, vertex) in diagram.vertices() {
        locals.push(split_local(&vertex.numerator)?);
    }
    for (_, _, edge) in diagram.edges() {
        locals.push(split_local(&edge.numerator)?);
    }
    // Cartesian product of the rows.
    let mut partial: Vec<(Atom, Vec<Atom>)> = vec![(Atom::num(1), Vec::new())];
    for rows in &locals {
        let mut next = Vec::with_capacity(partial.len() * rows.len());
        for (colour, residuals) in &partial {
            for (row_colour, row_residual) in rows {
                let mut residuals = residuals.clone();
                residuals.push(row_residual.clone());
                next.push((colour * row_colour, residuals));
            }
        }
        partial = next;
    }
    let total = partial.len();
    let mut branches = Vec::new();
    for (colour, residuals) in partial {
        let coefficient = colour_coefficient_lane(&colour)?;
        if !coefficient.is_zero() {
            branches.push(ColourBranch {
                coefficient,
                residuals,
            });
        }
    }
    Ok((branches, total))
}

/// The numeric factor of a product and the product without it.
fn split_number(atom: &Atom) -> (Atom, Atom) {
    match atom.as_view() {
        AtomView::Num(_) => (atom.clone(), Atom::num(1)),
        AtomView::Mul(m) => {
            let (mut number, mut rest) = (Atom::num(1), Atom::num(1));
            for part in m.iter() {
                if let AtomView::Num(_) = part {
                    number *= part.to_owned();
                } else {
                    rest *= part.to_owned();
                }
            }
            (number, rest)
        }
        _ => (Atom::num(1), atom.clone()),
    }
}

/// The branch-by-branch comparison of the lane's colour split with the input
/// type's on one diagram.
#[derive(Clone, Debug, Default)]
pub struct ColourSplitComparison {
    pub lane_total: usize,
    pub lane_nonzero: usize,
    pub input_total: usize,
    pub input_nonzero: usize,
    /// Non-zero branches present in one split only (by the colour key: the
    /// canonical strings of the colour rows without their numeric factors).
    pub only_lane: usize,
    pub only_input: usize,
    /// Branches in both whose coefficient per unit colour product differs.
    pub coefficient_differs: Vec<(String, String, String)>,
    /// Branches in both whose full content `coefficient × Π residuals`
    /// differs (expanded difference non-zero).
    pub content_differs: usize,
}

/// Compare the lane's colour split ([`colour_branches_lane`]) with the input
/// type's ([`colour_branches`]) on `diagram`, branch by branch: branches are
/// matched by their colour rows (canonical strings without numeric factors);
/// per matched branch the coefficient per unit colour product
/// (`coefficient / numeric factor of the colour product`) and the full
/// content `coefficient × Π residuals` (expanded) are compared exactly.
pub fn compare_colour_splits(diagram: &FeynmanDiagram) -> Result<ColourSplitComparison> {
    use crate::p_integral::factorized::{ColourMode, colour_coefficient, split_local_colour};
    let s = &*SYMBOLS;
    let mut comparison = ColourSplitComparison::default();
    // Rows of both splits, per local, with their keys.
    let locals = diagram
        .vertices()
        .map(|(_, v)| &v.numerator)
        .chain(diagram.edges().map(|(_, _, e)| &e.numerator))
        .collect::<Vec<_>>();
    type Rows = Vec<(String, Atom, Atom, Atom)>; // key, number, unit colour, residual
    let mut lane_rows = Vec::<Rows>::new();
    let mut input_rows = Vec::<Rows>::new();
    for local in &locals {
        let mut rows = Rows::new();
        for (colour, residual) in split_local(local)? {
            let (number, unit) = split_number(&colour);
            rows.push((unit.to_canonical_string(), number, unit, residual));
        }
        lane_rows.push(rows);
        let mut rows = Rows::new();
        for row in split_local_colour(&NumeratorInterface::export(*local).0)? {
            let (number, unit) = split_number(&row.colour);
            rows.push((unit.to_canonical_string(), number, unit, row.residual));
        }
        input_rows.push(rows);
    }
    let adjoint = Atom::var(symbol!("spenso::dA"));
    let fundamental = Atom::var(CS.nc);
    // key -> (coefficient per unit colour, content)
    let branches = |rows: &[Rows], lane: bool| -> Result<(BTreeMap<Vec<String>, (Atom, Atom)>, usize)> {
        let mut partial: Vec<(Vec<String>, Atom, Atom, Atom)> =
            vec![(Vec::new(), Atom::num(1), Atom::num(1), Atom::num(1))];
        for local in rows {
            let mut next = Vec::new();
            for (keys, number, unit, residual) in &partial {
                for (key, n, u, r) in local {
                    let mut keys = keys.clone();
                    keys.push(key.clone());
                    next.push((keys, number * n, unit * u, residual * r));
                }
            }
            partial = next;
        }
        let total = partial.len();
        let mut out = BTreeMap::new();
        for (keys, number, unit, residual) in partial {
            let colour = &number * &unit;
            let coefficient = if lane {
                colour_coefficient_lane(&colour)?
            } else {
                map_invariants(
                    &colour_coefficient(&colour, ColourMode::Invariants)?,
                    &adjoint,
                    &fundamental,
                )
                .expand()
                .together()
            };
            if coefficient.is_zero() {
                continue;
            }
            let per_unit = (&coefficient / &number).expand().together();
            let content = &coefficient * residual;
            out.insert(keys, (per_unit, content));
        }
        Ok((out, total))
    };
    let (lane, lane_total) = branches(&lane_rows, true)?;
    let (input, input_total) = branches(&input_rows, false)?;
    comparison.lane_total = lane_total;
    comparison.input_total = input_total;
    comparison.lane_nonzero = lane.len();
    comparison.input_nonzero = input.len();
    comparison.only_lane = lane.keys().filter(|k| !input.contains_key(*k)).count();
    comparison.only_input = input.keys().filter(|k| !lane.contains_key(*k)).count();
    let continued = |atom: &Atom| {
        atom.replace(Minkowski {}.to_symbolic([Atom::var(W_.d_), Atom::var(W_.a___)]))
            .with(Minkowski {}.to_symbolic([Atom::var(s.d), Atom::var(W_.a___)]))
    };
    for (keys, (lane_unit, lane_content)) in &lane {
        let Some((input_unit, input_content)) = input.get(keys) else {
            continue;
        };
        if !(lane_unit - input_unit).expand().together().is_zero() {
            comparison.coefficient_differs.push((
                keys.join(" * "),
                lane_unit.to_canonical_string(),
                input_unit.to_canonical_string(),
            ));
        }
        if !(continued(lane_content) - continued(input_content)).expand().is_zero() {
            comparison.content_differs += 1;
        }
    }
    Ok(comparison)
}

/// Closed quark loops of a diagram (components of the quark subgraph).
pub fn closed_quark_loops(diagram: &FeynmanDiagram) -> Result<usize> {
    let mut adjacent = BTreeMap::<usize, Vec<usize>>::new();
    for (_, endpoints, edge) in diagram.edges() {
        let particle = diagram.model().particle_by_id(edge.particle)?;
        if particle.spin != 2 {
            continue;
        }
        ensure!(edge.external.is_none(), "external quarks");
        let (Some(left), Some(right)) = (endpoints.source, endpoints.target) else {
            return Err(eyre!("an internal quark edge without two endpoints"));
        };
        adjacent.entry(left.0).or_default().push(right.0);
        adjacent.entry(right.0).or_default().push(left.0);
    }
    ensure!(
        adjacent.values().all(|neighbors| neighbors.len() == 2),
        "the quark subgraph is not a union of closed loops"
    );
    let mut seen = BTreeSet::new();
    let mut components = 0;
    for &start in adjacent.keys() {
        if seen.contains(&start) {
            continue;
        }
        components += 1;
        let mut stack = vec![start];
        while let Some(vertex) = stack.pop() {
            if seen.insert(vertex) {
                stack.extend(&adjacent[&vertex]);
            }
        }
    }
    Ok(components)
}

/// Whether no internal edge carries the external momentum (a contact graph
/// whose two background legs meet at one vertex: a scaleless vacuum graph,
/// an exact zero in dimensional regularisation).
pub fn scaleless_contact(diagram: &FeynmanDiagram) -> bool {
    let basis = diagram.loop_momentum_basis();
    diagram
        .edges()
        .filter(|(_, _, edge)| edge.external.is_none())
        .all(|(id, _, _)| {
            basis.edge_signatures.get(&id).is_none_or(|signature| {
                signature
                    .external
                    .integer_coefficients()
                    .iter()
                    .all(|&value| value == 0)
            })
        })
}

fn indexed(momentum: &Atom, slot: &Atom) -> Result<Atom> {
    match momentum.as_view() {
        AtomView::Var(variable) => Ok(variable.get_symbol().call(slot)),
        AtomView::Fun(fun) => Ok(FunctionBuilder::new(fun.get_symbol())
            .add_args(fun.iter())
            .add_arg(slot)
            .finish()),
        _ => Err(eyre!("expected a momentum head, got {momentum}")),
    }
}

/// The two contracted scalars of one branch: `S_g` (metric between the two
/// external Lorentz ports) and, unless `metric_only`, `S_qq` (`q^μ q^ν`),
/// from one simplification of `T·g + L·q q` (tags `T`, `L`).
fn contract_branch(
    product: &Atom,
    momentum: &Atom,
    metric_only: bool,
) -> Result<(Atom, Atom)> {
    let s = &*SYMBOLS;
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let cooked = cooking.try_cook(product.as_view())?;
    let ports = cooked
        .list_dangling::<AbstractIndex>()?
        .iter()
        .map(|port| LibrarySlot::<AbstractIndex>::try_from(port.as_view()))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    ensure!(
        ports.len() == 2
            && ports
                .iter()
                .all(|slot| slot.rep_name() == (Minkowski {}).into()),
        "a branch must have two dangling Lorentz ports, found {ports:?}"
    );
    let (left, right) = (ports[0].to_atom(), ports[1].to_atom());
    let metric = ETS.metric(&left, &right);
    let projector = if metric_only {
        metric
    } else {
        Atom::var(s.metric_tag) * metric
            + Atom::var(s.qq_tag) * indexed(momentum, &left)? * indexed(momentum, &right)?
    };
    let simplified = SymbolicTensor::infer_indexed(&cooked * &projector)?.simplify_algebra(
        &AlgebraSettings {
            gamma: Some(GammaSimplifySettings::canonical()),
            epsilon: true,
            ..Default::default()
        },
    )?;
    let status = simplified.reduction_status();
    let reduced = simplified.to_dots()?;
    // gammaloop's R* simplifier (`simplify_view`) takes the `to_dots` result
    // as it is; its scalar-ness is checked below and by vakint's ingress
    // (which refuses unconverted tensor notation).
    if std::env::var_os("REFPATH_DEBUG").is_some() {
        eprintln!(
            "[refpath] simplify status {status:?}, after to_dots {:?}: {}",
            reduced.reduction_status(),
            cooking.uncook(reduced.expression().as_view())
        );
    }
    let scalar = cooking.uncook(reduced.expression().as_view());
    ensure!(
        scalar.list_dangling::<AbstractIndex>().is_ok_and(|ports| ports.is_empty()),
        "the contraction kept open indices"
    );
    if metric_only {
        return Ok((scalar, Atom::Zero));
    }
    let g = scalar
        .replace(Atom::var(s.qq_tag))
        .with(Atom::Zero)
        .replace(Atom::var(s.metric_tag))
        .with(Atom::num(1));
    let qq = scalar
        .replace(Atom::var(s.metric_tag))
        .with(Atom::Zero)
        .replace(Atom::var(s.qq_tag))
        .with(Atom::num(1));
    Ok((g, qq))
}

/// The scalar numerators (before vakint) of the requested projections, the
/// sum over the colour branches weighted by their coefficients, with the beta
/// normalisation `1/(i q²)` (metric and transverse also `1/(D-1)`), couplings
/// expanded and `gs = 1`.
fn contract(
    diagram: &FeynmanDiagram,
    branches: &[ColourBranch],
    momentum: &Atom,
    q_squared: &Atom,
    projection: Projection,
) -> Result<Vec<(&'static str, Atom)>> {
    let s = &*SYMBOLS;
    let d = Atom::var(s.d);
    let metric_only = projection == Projection::Metric;
    let mut g_sum = Atom::Zero;
    let mut qq_sum = Atom::Zero;
    let basis = diagram.loop_momentum_basis();
    for branch in branches {
        let product = branch
            .residuals
            .iter()
            .fold(Atom::num(1), |product, residual| product * residual);
        let routed = basis.route_expression(&product);
        let continued = routed
            .replace(Minkowski {}.to_symbolic([Atom::var(W_.d_), Atom::var(W_.a___)]))
            .with(Minkowski {}.to_symbolic([d.clone(), Atom::var(W_.a___)]));
        let (g, qq) = contract_branch(&continued, momentum, metric_only)?;
        g_sum += &branch.coefficient * g;
        qq_sum += &branch.coefficient * qq;
    }
    let couplings = |value: Atom| {
        diagram
            .model()
            .expand_couplings(&value)
            .replace(Atom::var(symbol!("UFO::gs")))
            .with(Atom::num(1))
    };
    let i_q2 = Atom::i() * q_squared;
    let mut out = Vec::new();
    let transverse = || (&g_sum - &qq_sum / q_squared) / (&d - Atom::num(1)) / &i_q2;
    let longitudinal = || &qq_sum / q_squared / &i_q2;
    match projection {
        Projection::Metric => out.push(("metric", couplings(&g_sum / (&d - Atom::num(1)) / &i_q2))),
        Projection::Transverse => out.push(("transverse", couplings(transverse()))),
        Projection::Longitudinal => out.push(("longitudinal", couplings(longitudinal()))),
        Projection::Both => {
            out.push(("transverse", couplings(transverse())));
            out.push(("longitudinal", couplings(longitudinal())));
        }
    }
    Ok(out)
}

fn reset_peak_rss() -> bool {
    std::fs::write("/proc/self/clear_refs", "5").is_ok()
}

/// This process's peak RSS since the last reset, in bytes.
pub fn peak_rss_bytes() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find_map(|line| line.strip_prefix("VmHWM:"))?;
    let kilobytes: u64 = line.trim().trim_end_matches("kB").trim().parse().ok()?;
    Some(kilobytes * 1024)
}

/// Reset the peak RSS (Linux `clear_refs`); false where unsupported.
pub fn start_peak_rss() -> bool {
    reset_peak_rss()
}

/// Options of one graph run.
#[derive(Clone, Debug)]
pub struct RunOptions {
    pub projection: Projection,
    /// Laurent terms (`L + 1`: poles and the finite part).
    pub extra_terms: i64,
    /// Count the terms of the expanded numerator (statistics; costs an
    /// expansion outside every stage time).
    pub count_terms: bool,
    /// Stop after the contraction (no vakint evaluation).
    pub contract_only: bool,
}

impl Default for RunOptions {
    fn default() -> Self {
        Self {
            projection: Projection::Metric,
            extra_terms: 1,
            count_terms: true,
            contract_only: false,
        }
    }
}

/// Run the reference path on one DOT graph with `vakint` (an engine whose p
/// library covers the graph's loop count).
pub fn run_graph(
    vakint: &Vakint,
    model: &Arc<Model>,
    dot: &str,
    options: &RunOptions,
) -> Result<GraphRun> {
    let s = &*SYMBOLS;
    let mut seconds = StageSeconds::default();

    let started = Instant::now();
    let diagram = FeynmanDiagram::from_dot(Arc::clone(model), dot)?;
    seconds.parse = started.elapsed().as_secs_f64();
    ensure!(
        diagram.projector().expression().is_one()
            && diagram.numerator_prefactor().expression().is_one(),
        "the reference path expects projector 1 and numerator_prefactor 1"
    );
    let loops = diagram.loop_count();
    let quark_loops = closed_quark_loops(&diagram)?;
    let overall_factor =
        feynkit_graph::expressions::evaluate_overall_factor(diagram.overall_factor().as_view());
    let mut run = GraphRun {
        name: diagram.name().to_string(),
        loops,
        quark_loops,
        overall_factor: overall_factor.clone(),
        branches_total: 0,
        branches: 0,
        zero_reason: None,
        seconds,
        contracted_terms: 0,
        numerators: Vec::new(),
        results: Vec::new(),
    };
    if scaleless_contact(&diagram) {
        run.zero_reason = Some("no internal edge carries the external momentum: scaleless".into());
        return Ok(run);
    }

    let started = Instant::now();
    let (branches, total) = colour_branches(&diagram)?;
    run.seconds.colour = started.elapsed().as_secs_f64();
    run.branches_total = total;
    run.branches = branches.len();
    if branches.is_empty() {
        run.zero_reason = Some("every colour branch vanishes".into());
        return Ok(run);
    }

    let dimension = Atom::var(s.d);
    let kinematics = Kinematics::in_dimension(&dimension)?;
    let family = diagram.propagator_family(&kinematics)?;
    let [momentum] = family.external_momenta() else {
        return Err(eyre!("a two-point graph needs one independent external momentum"));
    };
    let momentum = momentum.clone();
    let q_squared = Kinematics::in_dimension(&dimension)?
        .with_momenta([momentum.clone()])?
        .scalar_product(&momentum, &momentum)?;

    let started = Instant::now();
    let numerators = contract(&diagram, &branches, &momentum, &q_squared, options.projection)?;
    run.seconds.contract = started.elapsed().as_secs_f64();
    let weight = overall_factor * Atom::var(s.nf).pow(Atom::num(quark_loops as i64));
    run.numerators = numerators
        .into_iter()
        .map(|(name, numerator)| (name, numerator * &weight))
        .collect();
    if options.count_terms {
        let started = Instant::now();
        run.contracted_terms = run
            .numerators
            .iter()
            .map(|(_, numerator)| match numerator.expand().as_view() {
                AtomView::Add(add) => add.get_nargs(),
                view if view.is_zero() => 0,
                _ => 1,
            })
            .sum();
        run.seconds.count = started.elapsed().as_secs_f64();
    }
    if options.contract_only {
        return Ok(run);
    }

    let settings = vakint_settings(loops as i64 + options.extra_terms);
    let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
    let in_epsilon = Atom::num(4) - Atom::num(2) * &epsilon;
    for (name, numerator) in &run.numerators {
        if numerator.is_zero() {
            continue;
        }
        let started = Instant::now();
        let input = match VakintExpression::from_diagram(
            &diagram,
            &family,
            numerator,
            &DiagramIntegralOptions::default(),
        ) {
            Ok(input) => Atom::from(input),
            Err(error) => {
                if let Some(directory) = std::env::var_os("REFPATH_DEBUG_DIR") {
                    let path = std::path::Path::new(&directory)
                        .join(format!("{}.{name}.numerator", run.name));
                    let _ = std::fs::write(path, numerator.to_canonical_string());
                }
                return Err(error.into());
            }
        };
        let input = input.replace(Atom::var(s.d)).with(in_epsilon.clone());
        let ingress = started.elapsed().as_secs_f64();
        let (value, observation) =
            p_gate::observe(|| vakint.evaluate(&settings, input.as_view()));
        let total = started.elapsed().as_secs_f64();
        let value = value?;
        let seconds = &mut run.seconds;
        let reduce = observation.prepare_seconds + observation.library_seconds;
        let finalise = observation.render_seconds + observation.finalize_seconds;
        seconds.ingress += ingress;
        seconds.split += observation.split_seconds;
        seconds.lowering += observation.lower_seconds;
        seconds.prepare += observation.prepare_seconds;
        seconds.library += observation.library_seconds;
        seconds.render += observation.render_seconds;
        seconds.finalize += observation.finalize_seconds;
        seconds.reduce += reduce;
        seconds.finalise += finalise;
        seconds.lower += total - reduce - finalise - observation.observer_seconds;
        run.results.push(ProjectedResult {
            projection: name,
            value,
            observation,
        });
    }
    Ok(run)
}

/// The gate view of a projected result without its zero-sector keys, and
/// the count of dropped keys.
pub fn gate_view(result: &ProjectedResult) -> (GateView, ZeroSectorCount) {
    result.observation.view_without_zero_sectors()
}

/// The coefficients of a Laurent expansion in `epsilon` from `ε^low` to
/// `ε^high`, each expanded (colour invariants and `n_f` symbolic).
pub fn laurent_coefficients(value: &Atom, epsilon: Symbol, low: i64, high: i64) -> Vec<(i64, Atom)> {
    let expanded = value.expand();
    (low..=high)
        .map(|order| {
            let coefficient = expanded
                .coefficient(Atom::var(epsilon).pow(Atom::num(order)).as_view())
                .expand();
            (order, coefficient)
        })
        .collect()
}
