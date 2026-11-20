//! The reference path R2 of the numerator assignment (helper VH, lane
//! `refpath`): one FeynKit DOT graph of the background-field `B -> B`
//! self-energy, from its per-edge / per-vertex spenso numerators to the
//! Laurent expansion of its contribution, through today's production pieces
//! only, with the stage times and the gate view of what RustRed receives.
//!
//! 1. parse: `Model::from_json`, `FeynmanDiagram::from_dot`;
//! 2. colour: the input type's colour split
//!    (`vakint::p_integral::factorized::{split_local_colour,
//!    colour_coefficient}`, the helper's one implementation: every local in
//!    source notation split into colour rows as the beta thread's `color.rs`
//!    `split_local`; the Cartesian product of the rows gives the branches, a
//!    quartic vertex has three; each branch's closed colour network, with the δ
//!    of the two external adjoint ports, reduced by idenso's
//!    `try_reduce_closed_color_network`, over `N_A`), written in this chain's
//!    symbols; branches with a zero coefficient are dropped;
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
use idenso::representations::{ColorAdjoint, ColorFundamental};
use idenso::tensor::{AlgebraSettings, SymbolicTensor};
use idenso::{CookMode, CookSettings, IndexTooling, W_};
use spenso::network::library::symbolic::ETS;
use spenso::structure::abstract_index::AbstractIndex;
use spenso::structure::representation::{LibrarySlot, Minkowski, RepName};
use spenso::structure::slot::IsAbstractSlot;
use symbolica::atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol};
use symbolica::symbol;

use crate::diagram_integral::DiagramIntegralOptions;
use crate::p_integral::FactorizedNumerator;
use crate::rustred_evaluation::p_gate::{self, GateObservation, GateView, ZeroSectorCount};
use crate::utils::vakint_macros::vk_symbol;
use crate::{Vakint, VakintExpression, VakintSettings};

/// R3: idenso's contraction in the edge basis handed directly to the lowering
/// (helper VH lane idenso-audit).
pub mod direct;

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

/// A new `Vakint` instance that shares `template`'s loaded p library (the immutable
/// `PLibrary`, loaded once) through a fresh handle: a cold request memo and
/// step cache per instance, the package's master catalogue, directory and
/// relations root kept. Repetitions on fresh instances then measure cold
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
pub fn map_invariants(
    value: &Atom,
    adjoint_dimension: &Atom,
    fundamental_dimension: &Atom,
) -> Atom {
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

/// The colour branches of a diagram by the input type's colour split
/// (`vakint::p_integral::factorized::{split_local_colour, colour_coefficient}`
/// in `ColourMode::Invariants`, the helper's one implementation), the
/// coefficients written in this chain's symbols ([`map_invariants`]); zero
/// coefficients dropped; with the number of branches before the drop.
pub fn colour_branches(diagram: &FeynmanDiagram) -> Result<(Vec<ColourBranch>, usize)> {
    use crate::p_integral::factorized::{ColourMode, colour_coefficient, split_local_colour};
    let mut locals = Vec::new();
    for (_, vertex) in diagram.vertices() {
        locals.push(split_local_colour(
            &NumeratorInterface::export(&vertex.numerator).0,
        )?);
    }
    for (_, _, edge) in diagram.edges() {
        locals.push(split_local_colour(
            &NumeratorInterface::export(&edge.numerator).0,
        )?);
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

/// The colour branches of the input type's values of one graph
/// ([`FactorizedNumerator::from_feynkit`] in `ColourMode::Invariants`, with
/// `dimension` = this chain's `D` and no projector: one value per non-zero
/// colour branch), for [`contract`]: the coefficient is the value's colour
/// coefficient in this chain's symbols ([`map_invariants`], as in
/// [`colour_branches`]) times its pulled scalars (couplings, `𝑖`, numbers,
/// the graph's `numerator_prefactor`); the one residual is the product of its
/// factors. The weight is [`factorized_weight`].
pub fn factorized_branches(values: &[FactorizedNumerator]) -> Result<Vec<ColourBranch>> {
    values
        .iter()
        .map(|value| {
            ensure!(
                value.projector.is_none() && value.dimension == SYMBOLS.d,
                "{}: the chain closes the ports itself (no projector) in its own D",
                value.label
            );
            Ok(ColourBranch {
                coefficient: factorized_colour(value) * &value.prefactor.scalars,
                residuals: vec![value.factor_product()],
            })
        })
        .collect()
}

/// The colour coefficient of an input-type value (`ColourMode::Invariants`)
/// in this chain's symbols, as [`colour_branches`] writes it: `C_A`, `C_F`,
/// `T_F`, `N_A`, higher invariants named (`gram(4, coad, cof)` ->
/// `refpath::gram_4AR`, ...).
pub fn factorized_colour(value: &FactorizedNumerator) -> Atom {
    let adjoint = Atom::var(symbol!("spenso::dA"));
    let fundamental = Atom::var(CS.nc);
    map_invariants(&value.prefactor.colour, &adjoint, &fundamental)
        .expand()
        .together()
}

/// The weight of the input type's values of one graph in this chain's
/// symbols: `overall × loop_sign × n_f^nf_power` ([`Prefactor`] without
/// scalars and colour), the same for every value (checked).
///
/// [`Prefactor`]: crate::p_integral::factorized::Prefactor
pub fn factorized_weight(values: &[FactorizedNumerator]) -> Result<Atom> {
    let mut weights = BTreeMap::new();
    for value in values {
        let p = &value.prefactor;
        let weight = &p.overall
            * Atom::num(p.loop_sign)
            * Atom::var(SYMBOLS.nf).pow(Atom::num(p.nf_power as i64));
        weights.insert(weight.to_canonical_string(), weight);
    }
    match weights.len() {
        1 => Ok(weights.into_values().next().unwrap()),
        n => Err(eyre!("{n} different weights over the colour branches")),
    }
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
fn contract_branch(product: &Atom, momentum: &Atom, metric_only: bool) -> Result<(Atom, Atom)> {
    let _step = crate::timing::contract();
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
        scalar
            .list_dangling::<AbstractIndex>()
            .is_ok_and(|ports| ports.is_empty()),
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
pub fn contract(
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

/// D6 (lane d6-open-indices): the colour-weighted sum of the branches
/// contracted by [`contract`]'s idenso chain with the external Lorentz ports
/// left OPEN (the indices that occur once, found as [`contract`] finds them),
/// couplings expanded (`gs = 1`), no projector and no normalisation. Returns
/// the open numerator and its ports (`mink(D, ..)` slots, canonical order).
pub fn contract_open(
    diagram: &FeynmanDiagram,
    branches: &[ColourBranch],
) -> Result<(Atom, Vec<Atom>)> {
    let s = &*SYMBOLS;
    let d = Atom::var(s.d);
    let basis = diagram.loop_momentum_basis();
    let mut sum = Atom::Zero;
    let mut ports: Option<Vec<Atom>> = None;
    for branch in branches {
        let product = branch
            .residuals
            .iter()
            .fold(Atom::num(1), |product, residual| product * residual);
        let routed = basis.route_expression(&product);
        let continued = routed
            .replace(Minkowski {}.to_symbolic([Atom::var(W_.d_), Atom::var(W_.a___)]))
            .with(Minkowski {}.to_symbolic([d.clone(), Atom::var(W_.a___)]));
        let (open, branch_ports) = contract_branch_open(&continued)?;
        match &ports {
            None => ports = Some(branch_ports),
            Some(known) => ensure!(
                *known == branch_ports,
                "the colour branches have different open ports: {known:?} and {branch_ports:?}"
            ),
        }
        sum += &branch.coefficient * open;
    }
    let value = diagram
        .model()
        .expand_couplings(&sum)
        .replace(Atom::var(symbol!("UFO::gs")))
        .with(Atom::num(1));
    Ok((value, ports.unwrap_or_default()))
}

/// [`contract_branch`] without a projector: the open Lorentz ports stay.
fn contract_branch_open(product: &Atom) -> Result<(Atom, Vec<Atom>)> {
    let _step = crate::timing::contract();
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let cooked = cooking.try_cook(product.as_view())?;
    let mut ports = Vec::new();
    for port in cooked.list_dangling::<AbstractIndex>()?.iter() {
        let slot = LibrarySlot::<AbstractIndex>::try_from(port.as_view())?;
        ensure!(
            slot.rep_name() == (Minkowski {}).into(),
            "an open index outside the Lorentz representation: {}",
            cooking.uncook(slot.to_atom().as_view())
        );
        ports.push(cooking.uncook(slot.to_atom().as_view()));
    }
    ports.sort_by_key(|port| port.to_canonical_string());
    let simplified = SymbolicTensor::infer_indexed(cooked)?.simplify_algebra(&AlgebraSettings {
        gamma: Some(GammaSimplifySettings::canonical()),
        epsilon: true,
        ..Default::default()
    })?;
    let reduced = simplified.to_dots()?;
    Ok((cooking.uncook(reduced.expression().as_view()), ports))
}

/// D6 (lane d6-open-indices): one graph through the OPEN-INDEX route, the
/// counterpart of [`run_graph`] without a projector: the colour branches
/// contracted with the external ports left open ([`contract_open`]), the
/// colour-weighted sum times `overall_factor` and `n_f` per closed quark loop
/// through vakint's FeynKit ingress (open indices mapped to indexed momenta;
/// ports renamed to the indices `1`, `2`), evaluated by `Vakint::evaluate`
/// (denominator match, tensor reduction with the external momentum declared,
/// lowering, RustRed, finaliser). The open value is `Π^{μν}` with `μ = 1`,
/// `ν = 2`; the caller closes it. Stage seconds as in [`run_graph`] (`lower`
/// = evaluation wall minus reduction, finalisation and observer: it contains
/// the denominator match and the tensor reduction); `tensor_probe` = a separate
/// timing of the denominator match plus the tensor reduction of the same
/// input on the same engine (`VakintExpression::canonicalize` +
/// `tensor_reduce`), reported apart.
pub struct OpenRun {
    pub name: String,
    pub loops: usize,
    pub zero_reason: Option<String>,
    pub seconds: StageSeconds,
    pub canonicalize_probe: f64,
    pub tensor_probe: f64,
    pub value: Atom,
    pub observation: Option<GateObservation>,
}

pub fn run_graph_open(
    vakint: &Vakint,
    model: &Arc<Model>,
    dot: &str,
    extra_terms: i64,
) -> Result<OpenRun> {
    let s = &*SYMBOLS;
    let started = Instant::now();
    let diagram = FeynmanDiagram::from_dot(Arc::clone(model), dot)?;
    let mut run = OpenRun {
        name: diagram.name().to_string(),
        loops: diagram.loop_count(),
        zero_reason: None,
        seconds: StageSeconds::default(),
        canonicalize_probe: 0.0,
        tensor_probe: 0.0,
        value: Atom::Zero,
        observation: None,
    };
    run.seconds.parse = started.elapsed().as_secs_f64();
    ensure!(
        diagram.projector().expression().is_one()
            && diagram.numerator_prefactor().expression().is_one(),
        "the open-index route expects projector 1 and numerator_prefactor 1"
    );
    if scaleless_contact(&diagram) {
        run.zero_reason = Some("no internal edge carries the external momentum: scaleless".into());
        return Ok(run);
    }
    let started = Instant::now();
    let (branches, _) = colour_branches(&diagram)?;
    run.seconds.colour = started.elapsed().as_secs_f64();
    if branches.is_empty() {
        run.zero_reason = Some("every colour branch vanishes".into());
        return Ok(run);
    }
    let started = Instant::now();
    let (open, ports) = contract_open(&diagram, &branches)?;
    run.seconds.contract = started.elapsed().as_secs_f64();
    ensure!(
        ports.len() == 2,
        "a two-point graph has two open ports, found {ports:?}"
    );
    let mut open = open;
    for (port, label) in ports.iter().zip([1, 2]) {
        let AtomView::Fun(slot) = port.as_view() else {
            return Err(eyre!("port {port} is not a slot"));
        };
        let index = slot
            .iter()
            .nth(1)
            .ok_or_else(|| eyre!("port {port} has no index"))?
            .to_owned();
        open = open
            .replace(index.to_pattern())
            .with(Atom::num(label).to_pattern());
    }
    let quark_loops = closed_quark_loops(&diagram)?;
    let weight =
        feynkit_graph::expressions::evaluate_overall_factor(diagram.overall_factor().as_view())
            * Atom::var(s.nf).pow(Atom::num(quark_loops as i64));
    let dimension = Atom::var(s.d);
    let kinematics = Kinematics::in_dimension(&dimension)?;
    let family = diagram.propagator_family(&kinematics)?;
    let settings = vakint_settings(run.loops as i64 + extra_terms);
    let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
    let in_epsilon = Atom::num(4) - Atom::num(2) * &epsilon;
    let started = Instant::now();
    let input = Atom::from(VakintExpression::from_diagram(
        &diagram,
        &family,
        &(open * &weight),
        &DiagramIntegralOptions::default(),
    )?)
    .replace(Atom::var(s.d))
    .with(in_epsilon);
    run.seconds.ingress = started.elapsed().as_secs_f64();
    // The probe: denominator match, then the tensor reduction, timed apart.
    let started = Instant::now();
    let mut probe = VakintExpression::try_from(input.as_view())?;
    probe.canonicalize(&settings, &vakint.topologies, false)?;
    run.canonicalize_probe = started.elapsed().as_secs_f64();
    let started = Instant::now();
    probe.tensor_reduce(vakint, &settings)?;
    run.tensor_probe = started.elapsed().as_secs_f64();
    // The evaluation proper.
    let started = Instant::now();
    let (value, observation) = p_gate::observe(|| vakint.evaluate(&settings, input.as_view()));
    let total = started.elapsed().as_secs_f64();
    run.value = value?;
    let seconds = &mut run.seconds;
    let reduce = observation.prepare_seconds + observation.library_seconds;
    let finalise = observation.render_seconds + observation.finalize_seconds;
    seconds.split = observation.split_seconds;
    seconds.lowering = observation.lower_seconds;
    seconds.prepare = observation.prepare_seconds;
    seconds.library = observation.library_seconds;
    seconds.render = observation.render_seconds;
    seconds.finalize = observation.finalize_seconds;
    seconds.reduce = reduce;
    seconds.finalise = finalise;
    seconds.lower = total - reduce - finalise - observation.observer_seconds;
    run.observation = Some(observation);
    Ok(run)
}

fn reset_peak_rss() -> bool {
    std::fs::write("/proc/self/clear_refs", "5").is_ok()
}

/// This process's peak RSS since the last reset, in bytes.
pub fn peak_rss_bytes() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let line = status
        .lines()
        .find_map(|line| line.strip_prefix("VmHWM:"))?;
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

/// Run the reference path on one DOT graph with `vakint` (a `Vakint` instance whose p
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
    // FeynKit reads a DOT without `overall_factor` as weight 1 (no 1/AutG, no
    // loop or ordering sign): the beta thread's DOT input must carry it.
    ensure!(
        crate::p_integral::factorized::dot_sets_graph_attribute(dot, "overall_factor"),
        "graph {}: the DOT graph attribute overall_factor is missing (required: 1/AutG and the \
         closed-loop and external ordering signs)",
        diagram.name()
    );
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
        return Err(eyre!(
            "a two-point graph needs one independent external momentum"
        ));
    };
    let momentum = momentum.clone();
    let q_squared = Kinematics::in_dimension(&dimension)?
        .with_momenta([momentum.clone()])?
        .scalar_product(&momentum, &momentum)?;

    let started = Instant::now();
    let numerators = contract(
        &diagram,
        &branches,
        &momentum,
        &q_squared,
        options.projection,
    )?;
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
        let (value, observation) = p_gate::observe(|| vakint.evaluate(&settings, input.as_view()));
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
pub fn laurent_coefficients(
    value: &Atom,
    epsilon: Symbol,
    low: i64,
    high: i64,
) -> Vec<(i64, Atom)> {
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
