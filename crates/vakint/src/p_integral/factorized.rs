//! A graph whose numerator is still a PRODUCT of spenso tensor factors.
//!
//! [`FactorizedNumerator`] is the input type of the numerator-in-vakint work
//! (lane `input-type`, deliverable D1): one value per graph and colour branch,
//! holding
//!
//! * the graph: edges with endpoints, half-edges, direction, line kind,
//!   power, mass and momentum as an INTEGER vector over the caller's loop
//!   momenta `k(1..L)` and independent external momenta `p(1..X)`;
//! * one spenso tensor expression per edge and per vertex (its *factor*) with
//!   explicit indices in the notation of FeynGen graphs: Lorentz slots
//!   `mink(D, hedge(h,1))`, spinor slots `bis(4, hedge(h,1))`, local dummies
//!   `mink(D, edge(e,k))` / `vertex(v,k)`, momenta `Q(e, mink(D, ..))` along
//!   the direction of edge `e`;
//! * the open external Lorentz indices (`ports`) and an optional projector
//!   over them;
//! * a scalar [`Prefactor`] (symmetry factor, closed-loop sign, couplings,
//!   colour coefficient, n_f power): free of loop momenta and of indices.
//!
//! Colour is never in a factor: one value per colour branch, the branch's
//! colour coefficient in the prefactor. The binding description of what a
//! value may contain is the written contract `CONTRACT-v1.md` of the helper
//! (N/notes/CONTRACT-v1.md); [`FactorizedNumerator::validate`] checks it and
//! returns a precise [`FactorizedError`] for the first violation found.
//!
//! The versioned text form ([`FACTORIZED_SCHEMA`]) is written by
//! [`FactorizedNumerator::to_text`] and read by
//! [`FactorizedNumerator::from_text`]; the round trip is exact.
//!
//! Reserved for the R* path (Phase 4) and REJECTED by v1: the closure kind
//! ([`FactorizedKind::Closure`], with an unapplied [`TaylorSpec`]) and
//! counterterm insertion vertices ([`Insertion`]).
//!
//! Nothing here changes an existing entry point: the existing
//! `DirectClosureInput` / `RemainderBoundary` paths are untouched.
//!
//! Numerator algebra (Lorentz and Dirac contraction, traces, colour) is done
//! by idenso / spenso only, never by vakint itself: a value is contracted with
//! idenso's production chain ([`contract_with_idenso`],
//! [`FactorizedNumerator::contract_naively`]) and the contracted numerator is
//! handed to vakint's existing lowering: [`evaluate_factorized`] does both
//! (value -> idenso contraction x prefactor -> `Vakint::evaluate`, the folded
//! inputs behind the gate observers). An earlier experimental contraction
//! engine for these values is kept for reference on the branch
//! `parked/vakint-numerator-engine` and is deprecated.

mod colour;
mod evaluate;
#[cfg(feature = "feynkit-ingress")]
mod feynkit;
#[cfg(feature = "feynkit-ingress")]
pub mod graph;
mod oracle;
mod text;
pub mod vacuum;

use std::collections::{BTreeMap, BTreeSet};

use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::coefficient::CoefficientView;

pub use colour::{
    ColourMode, ColourRow, MergedColourTerm, colour_coefficient, merge_colour_channels,
    split_local_colour,
};
pub use evaluate::{FactorizedValue, contract_to_vakint_dots, evaluate_factorized, vakint_input};
#[cfg(feature = "feynkit-ingress")]
pub use feynkit::{
    FromFeynkitOptions, OverallFactorSource, ProjectorChoice, dot_sets_graph_attribute,
};
#[cfg(feature = "feynkit-ingress")]
pub use graph::{GraphIntegral, GraphValue, evaluate_graph};
pub use oracle::contract_with_idenso;
pub use text::{parse_factorized_records, read_factorized_records};
pub use vacuum::{VacuumValue, evaluate_factorized_vacuum};

/// Schema line of the text form (first line of every record).
pub const FACTORIZED_SCHEMA: &str = "pint6.factorized-numerator.v1";

/// Fully qualified names of the symbols the contract refers to (DOT
/// notation of FeynGen graphs). Classification is by name, so values built in
/// another process (or read from text) are recognised.
pub(crate) mod names {
    pub const MINK: &str = "spenso::mink";
    pub const BIS: &str = "spenso::bis";
    pub const COAD: &str = "spenso::coad";
    pub const COF: &str = "spenso::cof";
    pub const COS: &str = "spenso::cos";
    pub const DIND: &str = "spenso::dind";
    pub const METRIC: &str = "spenso::g";
    pub const GAMMA: &str = "spenso::gamma";
    pub const GAMMA5: &str = "spenso::gamma5";
    pub const PROJP: &str = "spenso::projp";
    pub const PROJM: &str = "spenso::projm";
    pub const EPSILON: &str = "spenso::epsilon";
    pub const STRUCTURE_F: &str = "spenso::f";
    pub const GENERATOR_T: &str = "spenso::t";
    pub const MOMENTUM: &str = "gammalooprs::Q";
    pub const HEDGE: &str = "gammalooprs::hedge";
    /// The default `D` of [`super::FactorizedNumerator::dimension`]
    /// (gammaloop's `GS.dim`, what `Numerator::to_d_dim` maps `mink(4)` to).
    pub const DIMENSION: &str = "gammalooprs::dim";
    /// The symbol of n_f in [`super::Prefactor::total`].
    pub const NF: &str = "vakint::nf";
}

/// What a value denotes.
#[derive(Clone, Debug, PartialEq)]
pub enum FactorizedKind {
    /// A whole bare graph (the direct path of the assignment): the value is
    /// `prefactor × projector · Π factors / Π_e D_e^{power_e}` integrated over
    /// the loop momenta. Accepted by v1.
    Direct,
    /// RESERVED (Phase 4, R* remainders): a graph with a massive closing line
    /// and an unapplied Taylor expansion. [`FactorizedNumerator::validate`]
    /// and the text reader reject it in v1.
    Closure(ClosureSpec),
}

/// RESERVED (Phase 4): the closing line of an R* remainder (the v1 keys
/// `closure_momentum` / `closure_mass_squared` of `pint6.remainder-boundary.v1`)
/// and its unapplied Taylor step. Not accepted by v1.
#[derive(Clone, Debug, PartialEq)]
pub struct ClosureSpec {
    /// Index of the closing edge in [`FactorizedNumerator::edges`] order.
    pub edge: usize,
    pub mass_squared: Atom,
    pub power: i64,
    pub taylor: Option<TaylorSpec>,
}

/// RESERVED (Phase 4): an unapplied Taylor expansion in `lambda` of the
/// momenta of `scaled_edges` to `degree`, with mass compensation `M² → (1 −
/// λ²) M²` when `compensate_mass`. Not accepted by v1.
#[derive(Clone, Debug, PartialEq)]
pub struct TaylorSpec {
    pub lambda: Symbol,
    pub degree: u32,
    pub scaled_edges: Vec<usize>,
    pub compensate_mass: bool,
}

/// RESERVED (Phase 4): a contracted child counterterm acting as an extra
/// local vertex with open boundary slots, times its coefficient. Not accepted
/// by v1 (a non-empty [`FactorizedNumerator::insertions`] is rejected).
#[derive(Clone, Debug, PartialEq)]
pub struct Insertion {
    pub vertex: usize,
    pub coefficient: Atom,
    pub factor: Atom,
}

/// One end of an edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Endpoint {
    /// A vertex id of [`FactorizedNumerator::vertices`].
    Vertex(usize),
    /// The edge is an external leg on this side.
    External,
}

/// The spinor / flow type of a line; decides the fermion-line checks.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum LineKind {
    Scalar,
    /// A spin-1/2 line: directed, two spinor indices, closed loops traced.
    Fermion,
    Vector,
    /// A Faddeev-Popov ghost: directed, no spinor index; closed loops carry
    /// the same (−1) as fermion loops.
    Ghost,
}

impl LineKind {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Scalar => "scalar",
            Self::Fermion => "fermion",
            Self::Vector => "vector",
            Self::Ghost => "ghost",
        }
    }

    pub(crate) fn parse(text: &str) -> Option<Self> {
        Some(match text {
            "scalar" => Self::Scalar,
            "fermion" => Self::Fermion,
            "vector" => Self::Vector,
            "ghost" => Self::Ghost,
            _ => return None,
        })
    }

    fn is_directed_loop(self) -> bool {
        matches!(self, Self::Fermion | Self::Ghost)
    }
}

/// One edge (internal propagator or external leg).
#[derive(Clone, Debug, PartialEq)]
pub struct FactorizedEdge {
    /// The edge id `e` of `Q(e, ..)` and `edge(e, k)` (the DOT edge id).
    pub id: usize,
    pub source: Endpoint,
    pub target: Endpoint,
    /// Half-edge ids `h` (slots `hedge(h, k)`) at the source / target vertex;
    /// `None` exactly on an external side.
    pub source_hedge: Option<usize>,
    pub target_hedge: Option<usize>,
    /// Particle name (informational).
    pub particle: String,
    pub line: LineKind,
    /// Propagator power (internal edges: >= 1 in v1; ignored on external legs).
    pub power: i64,
    /// Mass squared (loop-free scalar; `0` for massless lines).
    pub mass_squared: Atom,
    /// `Q(e) = Σ_i momentum[i] k(i+1) + Σ_j momentum[L+j] p(j+1)` along
    /// source → target.
    pub momentum: Vec<i64>,
    /// The edge's tensor factor (colour-free; `1` when it has none).
    pub factor: Atom,
}

impl FactorizedEdge {
    pub fn is_external(&self) -> bool {
        self.source == Endpoint::External || self.target == Endpoint::External
    }
}

/// One vertex and its tensor factor.
#[derive(Clone, Debug, PartialEq)]
pub struct FactorizedVertex {
    pub id: usize,
    /// Vertex rule name (informational), e.g. `V_BQQ`.
    pub rule: String,
    /// The vertex's tensor factor (colour-free, one colour channel).
    pub factor: Atom,
}

/// A projector over the open external indices: the value is contracted with
/// `tensor × (p²)^p2_power`, `p = p(1)` (the first independent external
/// momentum; only meaningful when there is one, e.g. a two-point function).
#[derive(Clone, Debug, PartialEq)]
pub struct Projector {
    /// Free indices: exactly the ports. Allowed: metrics, external momenta
    /// `Q(x, mink(D, ..))` of EXTERNAL edges `x`, scalar coefficients in `D`.
    pub tensor: Atom,
    pub p2_power: i64,
}

/// The scalar that multiplies the contracted value, kept in named parts so
/// that the checkable ones can be checked. A consumer multiplies
/// [`Prefactor::total`] exactly once.
#[derive(Clone, Debug, PartialEq)]
pub struct Prefactor {
    /// Graph weight without the closed-loop sign (from a FeynGen graph:
    /// `overall_factor` with `InternalFermionLoopSign` removed, evaluated:
    /// `1/AutG × ExternalFermionOrderingSign`).
    pub overall: Atom,
    /// `(−1)^(closed fermion loops + closed ghost loops)`. Checked by
    /// `validate` against the loops of the graph: a consumer never adds it.
    pub loop_sign: i64,
    /// Power of n_f (= closed fermion loops). Checked by `validate`.
    pub nf_power: u32,
    /// Couplings, `𝑖`, rationals pulled out of the factors, the graph's
    /// `numerator_prefactor`, normalisations (may contain `D`).
    pub scalars: Atom,
    /// Colour coefficient of this branch (`δ_ab C_ab / N_A` for two adjoint
    /// external legs).
    pub colour: Atom,
}

impl Default for Prefactor {
    fn default() -> Self {
        Self {
            overall: Atom::num(1),
            loop_sign: 1,
            nf_power: 0,
            scalars: Atom::num(1),
            colour: Atom::num(1),
        }
    }
}

impl Prefactor {
    /// [`Self::total`] with the caller's flavour symbol in place of
    /// `vakint::nf` (e.g. the beta thread's `beta::nf`): the flavour factor is
    /// written once, in the symbol of whoever extracts its coefficient.
    pub fn total_with_flavour(&self, nf: &Atom) -> Atom {
        &self.overall
            * Atom::num(self.loop_sign)
            * nf.pow(Atom::num(self.nf_power as i64))
            * &self.scalars
            * &self.colour
    }

    /// `overall × loop_sign × nf^nf_power × scalars × colour`.
    pub fn total(&self) -> Atom {
        let nf = Atom::var(symbolica::symbol!(names::NF));
        &self.overall
            * Atom::num(self.loop_sign)
            * nf.pow(Atom::num(self.nf_power as i64))
            * &self.scalars
            * &self.colour
    }
}

/// Provenance of a colour branch (informational; not used in values).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ColourBranchInfo {
    /// Index of this branch among the non-zero branches of the graph.
    pub index: usize,
    /// Number of non-zero branches of the graph.
    pub count: usize,
    /// The colour row chosen per local factor, vertices (by id) then edges
    /// (by id); rows are ordered by the canonical string of their colour part.
    pub rows: Vec<usize>,
}

/// A closed fermion or ghost loop of the graph, edges in flow order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ClosedLoop {
    pub line: LineKind,
    pub edges: Vec<usize>,
}

/// The input type: one graph (one colour branch) with a factorised numerator.
/// See the module documentation and N/notes/CONTRACT-v1.md.
#[derive(Clone, Debug, PartialEq)]
pub struct FactorizedNumerator {
    /// Free-form provenance (graph name, colour branch): no TAB, no newline.
    pub label: String,
    pub kind: FactorizedKind,
    /// The symbol `D` of every Minkowski slot `mink(D, ..)`.
    pub dimension: Symbol,
    /// Loop momenta `k(1..L)`: the edge whose momentum is exactly `k(i)`, if
    /// any (checked to be the unit vector).
    pub loop_edges: Vec<Option<usize>>,
    /// Independent external momenta `p(1..X)`: the edge whose momentum is
    /// exactly `p(j)`, if any.
    pub external_edges: Vec<Option<usize>>,
    /// Sorted by id.
    pub edges: Vec<FactorizedEdge>,
    /// Sorted by id.
    pub vertices: Vec<FactorizedVertex>,
    /// Open external Lorentz indices of `Π factors`, sorted by canonical string.
    pub ports: Vec<Atom>,
    pub projector: Option<Projector>,
    pub prefactor: Prefactor,
    /// Informational.
    pub colour_branch: Option<ColourBranchInfo>,
    /// RESERVED: must be empty in v1.
    pub insertions: Vec<Insertion>,
}

/// Plain data for [`FactorizedNumerator::from_parts`].
#[derive(Clone, Debug)]
pub struct FactorizedParts {
    pub label: String,
    pub dimension: Symbol,
    pub loop_edges: Vec<Option<usize>>,
    pub external_edges: Vec<Option<usize>>,
    pub edges: Vec<FactorizedEdge>,
    pub vertices: Vec<FactorizedVertex>,
    pub ports: Vec<Atom>,
    pub projector: Option<Projector>,
    pub prefactor: Prefactor,
}

/// A violation of the contract, a text-form error or a constructor error.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum FactorizedError {
    #[error("{location}: index {slot} occurs {count} times (contracted more than once)")]
    IndexOverContracted {
        slot: String,
        count: usize,
        location: String,
    },
    #[error("index {slot} is not contracted (it occurs once, in {location}) and is not a port")]
    IndexUncontracted { slot: String, location: String },
    #[error("{location}: free index {slot} is not a half-edge index of this location")]
    IndexNotLocal { slot: String, location: String },
    #[error("dummy index {slot} of {location} also occurs in {other}")]
    DummyIndexShared {
        slot: String,
        location: String,
        other: String,
    },
    #[error("{location}: the terms of a sum have different free indices ({left} vs {right})")]
    InconsistentSum {
        location: String,
        left: String,
        right: String,
    },
    #[error("ports: {detail}")]
    Ports { detail: String },
    #[error("vertex {vertex}: momentum is not conserved (incoming − outgoing = {residual:?})")]
    MomentumNotConserved { vertex: usize, residual: Vec<i64> },
    #[error("edge {edge}: momentum vector of length {found}, expected {expected}")]
    MomentumLength {
        edge: usize,
        expected: usize,
        found: usize,
    },
    #[error("{detail}")]
    MomentumBasis { detail: String },
    #[error("{location}: Q({edge}, ..) refers to an unknown edge")]
    UnknownEdge { location: String, edge: String },
    #[error("graph: {detail}")]
    Graph { detail: String },
    #[error("{location}: colour tensor {tensor} in a factor (colour must be stripped)")]
    ColourInFactor { location: String, tensor: String },
    #[error("{location}: {tensor} (gamma5 / chiral projector) is not allowed in v1")]
    Gamma5 { location: String, tensor: String },
    #[error("{location}: Levi-Civita tensor {tensor} is not allowed in v1")]
    LeviCivita { location: String, tensor: String },
    #[error("{location}: {tensor} is not an allowed factor kind")]
    UnsupportedTensor { location: String, tensor: String },
    #[error("{location}: slot {slot} has dimension {found}, expected {expected}")]
    WrongDimension {
        location: String,
        slot: String,
        found: String,
        expected: String,
    },
    #[error("{location}: inexact or non-finite number {number}")]
    InexactNumber { location: String, number: String },
    #[error("{location}: {detail}")]
    NotScalar { location: String, detail: String },
    #[error("{location}: fermion line orientation: {detail}")]
    FermionOrientation { location: String, detail: String },
    #[error("vertex {vertex}: {line} flow not conserved ({incoming} in, {outgoing} out)")]
    FlowNotConserved {
        vertex: usize,
        line: String,
        incoming: usize,
        outgoing: usize,
    },
    #[error("edge {edge}: open {line} line (external fermion / ghost legs are not in v1)")]
    OpenLine { edge: usize, line: String },
    #[error("prefactor: loop_sign {found} but the graph has {loops} closed fermion/ghost loops")]
    LoopSign { found: i64, loops: usize },
    #[error("prefactor: nf_power {found} but the graph has {loops} closed fermion loops")]
    NfPower { found: u32, loops: usize },
    #[error("edge {edge}: power {power} (v1 needs >= 1 on internal edges)")]
    Power { edge: usize, power: i64 },
    #[error("{feature} is reserved for the R* path and not accepted by v1")]
    Reserved { feature: String },
    #[error("text line {line}: unknown schema {found:?} (this reader knows {FACTORIZED_SCHEMA})")]
    TextSchema { line: usize, found: String },
    #[error("text line {line}: `lines` asserts {asserted} lines, {found} were read")]
    TextLineCount {
        line: usize,
        asserted: usize,
        found: usize,
    },
    #[error("text line {line}: section {section} declares {declared} entries, found {found}")]
    TextSectionCount {
        line: usize,
        section: String,
        declared: usize,
        found: usize,
    },
    #[error("text line {line}: {detail}")]
    TextSyntax { line: usize, detail: String },
    #[error("constructor: {detail}")]
    Constructor { detail: String },
    #[error("graph {graph}: overall_factor: {detail}")]
    OverallFactor { graph: String, detail: String },
    #[error("colour: {detail}")]
    Colour { detail: String },
    #[error("contraction oracle: {detail}")]
    Oracle { detail: String },
    #[error("evaluation: {detail}")]
    Evaluation { detail: String },
}

type Result<T> = std::result::Result<T, FactorizedError>;

/// Where a factor sits (for error messages).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum Location {
    Edge(usize),
    Vertex(usize),
    Projector,
}

impl std::fmt::Display for Location {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Edge(id) => write!(f, "edge {id}"),
            Self::Vertex(id) => write!(f, "vertex {id}"),
            Self::Projector => write!(f, "projector"),
        }
    }
}

/// A representation slot `rep(dim, index)` found in a factor.
#[derive(Clone, Debug)]
pub(crate) struct Slot {
    pub atom: Atom,
    pub rep: String,
    pub dim: Atom,
    pub index: Atom,
}

impl Slot {
    pub(crate) fn key(&self) -> String {
        self.atom.to_canonical_string()
    }

    /// The half-edge `h` of an index `hedge(h, k)`.
    pub(crate) fn hedge(&self) -> Option<usize> {
        let AtomView::Fun(f) = self.index.as_view() else {
            return None;
        };
        if f.get_symbol().get_name() != names::HEDGE || f.get_nargs() != 2 {
            return None;
        }
        let first = f.iter().next()?;
        i64::try_from(first)
            .ok()
            .and_then(|h| usize::try_from(h).ok())
    }
}

pub(crate) fn is_rep_name(name: &str) -> bool {
    matches!(
        name,
        names::MINK | names::BIS | names::COAD | names::COF | names::COS | names::DIND
    )
}

pub(crate) fn is_colour_rep(name: &str) -> bool {
    matches!(name, names::COAD | names::COF | names::COS)
}

/// Read a slot `rep(dim, index)` or `dind(rep(dim, index))`.
pub(crate) fn as_slot(view: AtomView<'_>) -> Option<Slot> {
    let AtomView::Fun(f) = view else {
        return None;
    };
    let symbol = f.get_symbol();
    let name = symbol.get_name();
    if name == names::DIND {
        let inner = f.iter().next()?;
        let mut slot = as_slot(inner)?;
        slot.atom = view.to_owned();
        slot.rep = format!("dind({})", slot.rep);
        return Some(slot);
    }
    if !is_rep_name(name) || f.get_nargs() != 2 {
        return None;
    }
    let mut args = f.iter();
    let dim = args.next()?.to_owned();
    let index = args.next()?.to_owned();
    Some(Slot {
        atom: view.to_owned(),
        rep: name.to_string(),
        dim,
        index,
    })
}

/// Whether `view` contains a representation slot `rep(dim, index)` or a
/// momentum `Q(..)` anywhere.
pub(crate) fn contains_tensor_content(view: AtomView<'_>) -> bool {
    let mut found = false;
    view.visitor(&mut |part| {
        if found {
            return false;
        }
        if let AtomView::Fun(f) = part {
            let symbol = f.get_symbol();
            let name = symbol.get_name();
            // A slot is `rep(dim, index)`; `rep(dim)` alone (as in a colour
            // invariant `cas(2, coad(dA))`) is not an index.
            if (is_rep_name(name) && name != names::DIND && f.get_nargs() == 2)
                || name == names::MOMENTUM
            {
                found = true;
                return false;
            }
        }
        true
    });
    found
}

/// Whether `view` contains a colour slot `rep(dim, index)` anywhere.
pub(crate) fn contains_colour_slot(view: AtomView<'_>) -> bool {
    let mut found = false;
    view.visitor(&mut |part| {
        if found {
            return false;
        }
        if let AtomView::Fun(f) = part
            && is_colour_rep(f.get_symbol().get_name())
            && f.get_nargs() == 2
        {
            found = true;
            return false;
        }
        true
    });
    found
}

fn check_exact_numbers(view: AtomView<'_>, location: &str) -> Result<()> {
    let mut bad = None;
    view.visitor(&mut |part| {
        if bad.is_some() {
            return false;
        }
        if let AtomView::Num(n) = part
            && !matches!(
                n.get_coeff_view(),
                CoefficientView::Natural(..) | CoefficientView::Large(..)
            )
        {
            bad = Some(part.to_owned().to_canonical_string());
            return false;
        }
        true
    });
    match bad {
        Some(number) => Err(FactorizedError::InexactNumber {
            location: location.to_string(),
            number,
        }),
        None => Ok(()),
    }
}

/// A scalar part of the value: exact numbers, no slot, no momentum.
fn check_scalar(view: AtomView<'_>, location: &str) -> Result<()> {
    check_exact_numbers(view, location)?;
    if contains_tensor_content(view) {
        return Err(FactorizedError::NotScalar {
            location: location.to_string(),
            detail: format!(
                "must be free of indices and momenta: {}",
                view.to_owned().to_canonical_string()
            ),
        });
    }
    Ok(())
}

/// Free and contracted (dummy) indices of one factor, by slot key.
#[derive(Clone, Debug, Default)]
pub(crate) struct IndexSets {
    pub free: BTreeMap<String, Slot>,
    pub dummies: BTreeMap<String, Slot>,
}

/// Elementary tensors seen in a factor (for the fermion-orientation check):
/// `(row slot, column slot)` of every `gamma`.
#[derive(Clone, Debug, Default)]
pub(crate) struct Elements {
    pub gammas: Vec<(Slot, Slot)>,
}

/// What a factor at `location` may contain.
struct FactorRules<'a> {
    location: Location,
    dimension: Symbol,
    edges: &'a BTreeMap<usize, &'a FactorizedEdge>,
}

impl FactorRules<'_> {
    fn loc(&self) -> String {
        self.location.to_string()
    }

    fn check_slot(&self, slot: &Slot) -> Result<()> {
        let rep = slot.rep.as_str();
        if rep.starts_with("dind(") || is_colour_rep(rep) {
            return Err(FactorizedError::ColourInFactor {
                location: self.loc(),
                tensor: slot.atom.to_canonical_string(),
            });
        }
        let expected = match rep {
            names::MINK => Atom::var(self.dimension),
            names::BIS => Atom::num(4),
            _ => {
                return Err(FactorizedError::UnsupportedTensor {
                    location: self.loc(),
                    tensor: slot.atom.to_canonical_string(),
                });
            }
        };
        if slot.dim != expected {
            return Err(FactorizedError::WrongDimension {
                location: self.loc(),
                slot: slot.key(),
                found: slot.dim.to_canonical_string(),
                expected: expected.to_canonical_string(),
            });
        }
        Ok(())
    }

    /// The slots of one elementary tensor `f`, after checking its kind.
    fn element(&self, view: AtomView<'_>, elements: &mut Elements) -> Result<Vec<Slot>> {
        let AtomView::Fun(f) = view else {
            unreachable!()
        };
        let symbol = f.get_symbol();
        let name = symbol.get_name();
        let text = || view.to_owned().to_canonical_string();
        let args = f.iter().collect::<Vec<_>>();
        let slots = args.iter().filter_map(|a| as_slot(*a)).collect::<Vec<_>>();
        // Colour first: a colour tensor is reported as such, whatever its head.
        if slots
            .iter()
            .any(|s| s.rep.starts_with("dind(") || is_colour_rep(&s.rep))
        {
            return Err(FactorizedError::ColourInFactor {
                location: self.loc(),
                tensor: text(),
            });
        }
        match name {
            names::GAMMA5 | names::PROJP | names::PROJM => {
                return Err(FactorizedError::Gamma5 {
                    location: self.loc(),
                    tensor: text(),
                });
            }
            names::EPSILON => {
                return Err(FactorizedError::LeviCivita {
                    location: self.loc(),
                    tensor: text(),
                });
            }
            names::STRUCTURE_F | names::GENERATOR_T => {
                return Err(FactorizedError::ColourInFactor {
                    location: self.loc(),
                    tensor: text(),
                });
            }
            _ => {}
        }
        for slot in &slots {
            self.check_slot(slot)?;
        }
        let unsupported = || FactorizedError::UnsupportedTensor {
            location: self.loc(),
            tensor: text(),
        };
        match name {
            names::METRIC => {
                if args.len() != 2 || slots.len() != 2 || slots[0].rep != slots[1].rep {
                    return Err(unsupported());
                }
                if self.location == Location::Projector && slots[0].rep != names::MINK {
                    return Err(unsupported());
                }
                Ok(slots)
            }
            names::GAMMA => {
                if self.location == Location::Projector
                    || args.len() != 3
                    || slots.len() != 3
                    || slots[0].rep != names::BIS
                    || slots[1].rep != names::BIS
                    || slots[2].rep != names::MINK
                {
                    return Err(unsupported());
                }
                elements.gammas.push((slots[0].clone(), slots[1].clone()));
                Ok(slots)
            }
            names::MOMENTUM => {
                if args.len() != 2 || slots.len() != 1 || slots[0].rep != names::MINK {
                    return Err(unsupported());
                }
                let edge = i64::try_from(args[0])
                    .ok()
                    .and_then(|e| usize::try_from(e).ok());
                let Some(edge) = edge else {
                    return Err(FactorizedError::UnknownEdge {
                        location: self.loc(),
                        edge: args[0].to_owned().to_canonical_string(),
                    });
                };
                // `Q(x, ..)` of an edge that is not in the graph is a constant
                // external vector (reader rule 2: derived, not declared).
                let Some(known) = self.edges.get(&edge) else {
                    return Ok(slots);
                };
                if self.location == Location::Projector && !known.is_external() {
                    return Err(FactorizedError::UnsupportedTensor {
                        location: self.loc(),
                        tensor: format!("{} (loop momentum in the projector)", text()),
                    });
                }
                Ok(slots)
            }
            _ => {
                if contains_tensor_content(view) {
                    Err(unsupported())
                } else {
                    // A scalar function (e.g. a coupling function): allowed.
                    check_exact_numbers(view, &self.loc())?;
                    Ok(Vec::new())
                }
            }
        }
    }

    /// The slots of `view` that are not contracted inside a product or a sum
    /// term below it (a multiset: a product folds its own slots, so an
    /// enclosing product sees each slot as often as it occurs), collecting
    /// the contracted ones into `sets.dummies`.
    fn free(
        &self,
        view: AtomView<'_>,
        sets: &mut IndexSets,
        elements: &mut Elements,
    ) -> Result<Vec<Slot>> {
        match view {
            AtomView::Num(_) => {
                check_exact_numbers(view, &self.loc())?;
                Ok(Vec::new())
            }
            AtomView::Var(_) => Ok(Vec::new()),
            AtomView::Fun(f) => {
                if is_rep_name(f.get_symbol().get_name()) {
                    return Err(FactorizedError::UnsupportedTensor {
                        location: self.loc(),
                        tensor: format!("{} (a bare slot)", view.to_owned().to_canonical_string()),
                    });
                }
                self.element(view, elements)
            }
            AtomView::Mul(m) => {
                let mut all = Vec::new();
                for part in m.iter() {
                    all.extend(self.free(part, sets, elements)?);
                }
                self.fold(all, sets)
            }
            AtomView::Add(a) => {
                let mut first: Option<(BTreeSet<String>, Vec<Slot>)> = None;
                for term in a.iter() {
                    let free = self.free(term, sets, elements)?;
                    let free = self.fold(free, sets)?;
                    let keys = free.iter().map(Slot::key).collect::<BTreeSet<_>>();
                    match &first {
                        None => first = Some((keys, free)),
                        Some((left, _)) if *left != keys => {
                            return Err(FactorizedError::InconsistentSum {
                                location: self.loc(),
                                left: left.iter().cloned().collect::<Vec<_>>().join(","),
                                right: keys.into_iter().collect::<Vec<_>>().join(","),
                            });
                        }
                        Some(_) => {}
                    }
                }
                Ok(first.map(|(_, free)| free).unwrap_or_default())
            }
            AtomView::Pow(p) => {
                let (base, exponent) = p.get_base_exp();
                let base_free = self.free(base, sets, elements)?;
                if base_free.is_empty() && !contains_tensor_content(base) {
                    check_scalar(exponent, &self.loc())?;
                    return Ok(Vec::new());
                }
                let n = i64::try_from(exponent).ok().filter(|n| *n >= 1);
                let Some(n) = n else {
                    return Err(FactorizedError::UnsupportedTensor {
                        location: self.loc(),
                        tensor: format!(
                            "{} (a tensor to a non-positive or symbolic power)",
                            view.to_owned().to_canonical_string()
                        ),
                    });
                };
                let mut all = Vec::new();
                for _ in 0..n {
                    all.extend(base_free.iter().cloned());
                }
                Ok(all)
            }
        }
    }

    /// Pair up repeated slots of a product: twice = contracted, more = error.
    fn fold(&self, all: Vec<Slot>, sets: &mut IndexSets) -> Result<Vec<Slot>> {
        let mut counts = BTreeMap::<String, (usize, Slot)>::new();
        let mut order = Vec::new();
        for slot in all {
            let key = slot.key();
            let entry = counts.entry(key.clone()).or_insert_with(|| {
                order.push(key.clone());
                (0, slot.clone())
            });
            entry.0 += 1;
        }
        let mut free = Vec::new();
        for key in order {
            let (count, slot) = counts.remove(&key).unwrap();
            match count {
                1 => free.push(slot),
                2 => {
                    sets.dummies.insert(key, slot);
                }
                _ => {
                    return Err(FactorizedError::IndexOverContracted {
                        slot: key,
                        count,
                        location: self.loc(),
                    });
                }
            }
        }
        Ok(free)
    }

    fn analyse(&self, factor: &Atom) -> Result<(IndexSets, Elements)> {
        let mut sets = IndexSets::default();
        let mut elements = Elements::default();
        let free = self.free(factor.as_view(), &mut sets, &mut elements)?;
        let free = self.fold(free, &mut sets)?;
        for slot in free {
            let key = slot.key();
            if sets.free.insert(key.clone(), slot).is_some() {
                return Err(FactorizedError::IndexOverContracted {
                    slot: key,
                    count: 2,
                    location: self.loc(),
                });
            }
        }
        Ok((sets, elements))
    }
}

impl FactorizedNumerator {
    /// A value from plain data, checked by [`Self::validate`].
    pub fn from_parts(parts: FactorizedParts) -> Result<Self> {
        let FactorizedParts {
            label,
            dimension,
            loop_edges,
            external_edges,
            mut edges,
            mut vertices,
            mut ports,
            projector,
            prefactor,
        } = parts;
        edges.sort_by_key(|edge| edge.id);
        vertices.sort_by_key(|vertex| vertex.id);
        ports.sort_by_key(|port| port.to_canonical_string());
        let value = Self {
            label,
            kind: FactorizedKind::Direct,
            dimension,
            loop_edges,
            external_edges,
            edges,
            vertices,
            ports,
            projector,
            prefactor,
            colour_branch: None,
            insertions: Vec::new(),
        };
        value.validate()?;
        Ok(value)
    }

    /// The default `D` of the Minkowski slots, `gammalooprs::dim`
    /// (gammaloop's `GS.dim`, what its `Numerator::to_d_dim` uses).
    pub fn default_dimension() -> Symbol {
        symbolica::symbol!(names::DIMENSION)
    }

    /// Number of loop momenta `L`.
    pub fn loops(&self) -> usize {
        self.loop_edges.len()
    }

    /// Number of independent external momenta `X`.
    pub fn externals(&self) -> usize {
        self.external_edges.len()
    }

    pub fn edge(&self, id: usize) -> Option<&FactorizedEdge> {
        self.edges.iter().find(|edge| edge.id == id)
    }

    pub fn vertex(&self, id: usize) -> Option<&FactorizedVertex> {
        self.vertices.iter().find(|vertex| vertex.id == id)
    }

    /// The product of all edge and vertex factors (colour-free), as one atom.
    pub fn factor_product(&self) -> Atom {
        let mut product = Atom::num(1);
        for edge in &self.edges {
            product *= &edge.factor;
        }
        for vertex in &self.vertices {
            product *= &vertex.factor;
        }
        product
    }

    /// Check the contract (CONTRACT-v1.md). Returns the first violation.
    pub fn validate(&self) -> Result<()> {
        if let FactorizedKind::Closure(_) = self.kind {
            return Err(FactorizedError::Reserved {
                feature: "the closure kind (R* remainder with a closing line)".into(),
            });
        }
        if !self.insertions.is_empty() {
            return Err(FactorizedError::Reserved {
                feature: "insertion vertices (child counterterms)".into(),
            });
        }
        if self.label.contains(['\t', '\n', '\r']) {
            return Err(FactorizedError::Graph {
                detail: "the label contains a TAB or a newline".into(),
            });
        }
        let edges = self.check_graph()?;
        self.check_momenta()?;
        // Factors: kinds, dimensions, indices.
        let mut analysed = Vec::<(Location, IndexSets, Elements)>::new();
        for edge in &self.edges {
            let rules = FactorRules {
                location: Location::Edge(edge.id),
                dimension: self.dimension,
                edges: &edges,
            };
            let (sets, elements) = rules.analyse(&edge.factor)?;
            analysed.push((rules.location, sets, elements));
        }
        for vertex in &self.vertices {
            let rules = FactorRules {
                location: Location::Vertex(vertex.id),
                dimension: self.dimension,
                edges: &edges,
            };
            let (sets, elements) = rules.analyse(&vertex.factor)?;
            analysed.push((rules.location, sets, elements));
        }
        if let Some(projector) = &self.projector {
            let rules = FactorRules {
                location: Location::Projector,
                dimension: self.dimension,
                edges: &edges,
            };
            let (sets, elements) = rules.analyse(&projector.tensor)?;
            analysed.push((rules.location, sets, elements));
        }
        self.check_locality(&analysed)?;
        self.check_contraction(&analysed)?;
        self.check_fermion_orientation(&analysed)?;
        // Prefactor and masses.
        let p = &self.prefactor;
        check_scalar(p.overall.as_view(), "prefactor overall")?;
        check_scalar(p.scalars.as_view(), "prefactor scalars")?;
        check_scalar(p.colour.as_view(), "prefactor colour")?;
        if contains_colour_slot(p.colour.as_view()) {
            return Err(FactorizedError::NotScalar {
                location: "prefactor colour".into(),
                detail: "an open colour index".into(),
            });
        }
        for edge in &self.edges {
            check_scalar(
                edge.mass_squared.as_view(),
                &format!("edge {} mass", edge.id),
            )?;
        }
        let loops = self.closed_loops()?;
        let n_loops = loops.len();
        let n_fermion = loops.iter().filter(|l| l.line == LineKind::Fermion).count();
        let expected_sign = if n_loops % 2 == 0 { 1 } else { -1 };
        if p.loop_sign != expected_sign {
            return Err(FactorizedError::LoopSign {
                found: p.loop_sign,
                loops: n_loops,
            });
        }
        if p.nf_power as usize != n_fermion {
            return Err(FactorizedError::NfPower {
                found: p.nf_power,
                loops: n_fermion,
            });
        }
        Ok(())
    }

    /// Ids unique, endpoints and half-edges consistent, powers; returns the
    /// edges by id.
    fn check_graph(&self) -> Result<BTreeMap<usize, &FactorizedEdge>> {
        let graph = |detail: String| FactorizedError::Graph { detail };
        let mut edges = BTreeMap::new();
        for edge in &self.edges {
            if edges.insert(edge.id, edge).is_some() {
                return Err(graph(format!("edge id {} occurs twice", edge.id)));
            }
        }
        let mut vertices = BTreeSet::new();
        for vertex in &self.vertices {
            if !vertices.insert(vertex.id) {
                return Err(graph(format!("vertex id {} occurs twice", vertex.id)));
            }
            if vertex.rule.contains(['\t', '\n', '=', ',']) {
                return Err(graph(format!("vertex {}: bad rule name", vertex.id)));
            }
        }
        let mut hedges = BTreeSet::new();
        for edge in &self.edges {
            if edge.particle.contains(['\t', '\n', '=', ',']) || edge.particle.is_empty() {
                return Err(graph(format!("edge {}: bad particle name", edge.id)));
            }
            for (end, hedge, side) in [
                (edge.source, edge.source_hedge, "source"),
                (edge.target, edge.target_hedge, "target"),
            ] {
                match (end, hedge) {
                    (Endpoint::Vertex(v), Some(h)) => {
                        if !vertices.contains(&v) {
                            return Err(graph(format!(
                                "edge {}: {side} vertex {v} does not exist",
                                edge.id
                            )));
                        }
                        if !hedges.insert(h) {
                            return Err(graph(format!("half-edge {h} occurs twice")));
                        }
                    }
                    (Endpoint::External, None) => {}
                    _ => {
                        return Err(graph(format!(
                            "edge {}: the {side} half-edge must be given exactly on a vertex side",
                            edge.id
                        )));
                    }
                }
            }
            if edge.source == Endpoint::External && edge.target == Endpoint::External {
                return Err(graph(format!("edge {}: both ends external", edge.id)));
            }
            if !edge.is_external() && edge.power < 1 {
                return Err(FactorizedError::Power {
                    edge: edge.id,
                    power: edge.power,
                });
            }
        }
        for window in self.edges.windows(2) {
            if window[0].id >= window[1].id {
                return Err(graph("edges are not sorted by id".into()));
            }
        }
        for window in self.vertices.windows(2) {
            if window[0].id >= window[1].id {
                return Err(graph("vertices are not sorted by id".into()));
            }
        }
        Ok(edges)
    }

    fn check_momenta(&self) -> Result<()> {
        let n = self.loops() + self.externals();
        for edge in &self.edges {
            if edge.momentum.len() != n {
                return Err(FactorizedError::MomentumLength {
                    edge: edge.id,
                    expected: n,
                    found: edge.momentum.len(),
                });
            }
        }
        let basis = self
            .loop_edges
            .iter()
            .enumerate()
            .map(|(i, e)| (format!("k({})", i + 1), i, *e))
            .chain(
                self.external_edges
                    .iter()
                    .enumerate()
                    .map(|(j, e)| (format!("p({})", j + 1), self.loops() + j, *e)),
            );
        for (name, position, edge) in basis {
            let Some(id) = edge else { continue };
            let Some(edge) = self.edge(id) else {
                return Err(FactorizedError::MomentumBasis {
                    detail: format!("{name} names edge {id}, which does not exist"),
                });
            };
            let unit = (0..n).map(|i| i64::from(i == position)).collect::<Vec<_>>();
            if edge.momentum != unit {
                return Err(FactorizedError::MomentumBasis {
                    detail: format!(
                        "{name} names edge {id}, whose momentum is {:?}",
                        edge.momentum
                    ),
                });
            }
        }
        for vertex in &self.vertices {
            let mut residual = vec![0i64; n];
            for edge in &self.edges {
                if edge.target == Endpoint::Vertex(vertex.id) {
                    for (r, c) in residual.iter_mut().zip(&edge.momentum) {
                        *r += c;
                    }
                }
                if edge.source == Endpoint::Vertex(vertex.id) {
                    for (r, c) in residual.iter_mut().zip(&edge.momentum) {
                        *r -= c;
                    }
                }
            }
            if residual.iter().any(|c| *c != 0) {
                return Err(FactorizedError::MomentumNotConserved {
                    vertex: vertex.id,
                    residual,
                });
            }
        }
        Ok(())
    }

    /// Every half-edge of the graph (both ends of every edge that has them).
    pub(crate) fn hedges(&self) -> BTreeSet<usize> {
        self.edges
            .iter()
            .flat_map(|e| e.source_hedge.into_iter().chain(e.target_hedge))
            .collect()
    }

    /// The CONSTANT EXTERNAL VECTORS: the ids `x` of `Q(x, ..)` in the factors
    /// and the projector with `x` not an edge of the graph (reader rule 2),
    /// sorted. A vacuum value takes them as the external vectors after the
    /// graph's own external momenta.
    pub fn constant_vectors(&self) -> Vec<usize> {
        let mut ids = BTreeSet::new();
        let factors = self
            .edges
            .iter()
            .map(|e| &e.factor)
            .chain(self.vertices.iter().map(|v| &v.factor))
            .chain(self.projector.iter().map(|p| &p.tensor));
        for factor in factors {
            factor.visitor(&mut |part| {
                if let AtomView::Fun(f) = part
                    && f.get_symbol().get_name() == names::MOMENTUM
                    && f.get_nargs() == 2
                    && let Some(x) = f
                        .iter()
                        .next()
                        .and_then(|a| i64::try_from(a).ok())
                        .and_then(|a| usize::try_from(a).ok())
                    && self.edge(x).is_none()
                {
                    ids.insert(x);
                }
                true
            });
        }
        ids.into_iter().collect()
    }

    /// Half-edges at a location: an edge's own two, a vertex's incident ones.
    fn local_hedges(&self, location: Location) -> BTreeSet<usize> {
        match location {
            Location::Edge(id) => self
                .edge(id)
                .map(|e| e.source_hedge.into_iter().chain(e.target_hedge).collect())
                .unwrap_or_default(),
            Location::Vertex(v) => self
                .edges
                .iter()
                .flat_map(|e| {
                    let s = (e.source == Endpoint::Vertex(v)).then_some(e.source_hedge);
                    let t = (e.target == Endpoint::Vertex(v)).then_some(e.target_hedge);
                    s.into_iter().chain(t).flatten()
                })
                .collect(),
            Location::Projector => BTreeSet::new(),
        }
    }

    /// Every free index of an edge / vertex factor is a half-edge index of
    /// that location.
    fn check_locality(&self, analysed: &[(Location, IndexSets, Elements)]) -> Result<()> {
        let present = self.hedges();
        let ports = self
            .ports
            .iter()
            .map(|p| p.to_canonical_string())
            .collect::<BTreeSet<_>>();
        for (location, sets, _) in analysed {
            if *location == Location::Projector {
                continue;
            }
            let local = self.local_hedges(*location);
            for (key, slot) in &sets.free {
                // A slot of a half-edge the graph does not have is an open port
                // (reader rule 1), accepted when it is one of the ports.
                let absent_port =
                    slot.hedge().is_some_and(|h| !present.contains(&h)) && ports.contains(key);
                if !absent_port && !slot.hedge().is_some_and(|h| local.contains(&h)) {
                    return Err(FactorizedError::IndexNotLocal {
                        slot: key.clone(),
                        location: location.to_string(),
                    });
                }
            }
        }
        Ok(())
    }

    /// Every index is contracted exactly once across the factors, or is a
    /// port (once in the factors, once in the projector when there is one).
    fn check_contraction(&self, analysed: &[(Location, IndexSets, Elements)]) -> Result<()> {
        let mut free = BTreeMap::<String, Vec<Location>>::new();
        let mut dummies = BTreeMap::<String, Location>::new();
        for (location, sets, _) in analysed {
            for key in sets.free.keys() {
                free.entry(key.clone()).or_default().push(*location);
            }
            for key in sets.dummies.keys() {
                if let Some(other) = dummies.insert(key.clone(), *location) {
                    return Err(FactorizedError::DummyIndexShared {
                        slot: key.clone(),
                        location: other.to_string(),
                        other: location.to_string(),
                    });
                }
            }
        }
        for (key, location) in &dummies {
            if let Some(others) = free.get(key) {
                return Err(FactorizedError::DummyIndexShared {
                    slot: key.clone(),
                    location: location.to_string(),
                    other: others[0].to_string(),
                });
            }
        }
        let ports = self
            .ports
            .iter()
            .map(|p| p.to_canonical_string())
            .collect::<BTreeSet<_>>();
        if ports.len() != self.ports.len() {
            return Err(FactorizedError::Ports {
                detail: "a port is listed twice".into(),
            });
        }
        for port in &self.ports {
            match as_slot(port.as_view()) {
                Some(slot) if slot.rep == names::MINK || slot.rep == names::BIS => {}
                _ => {
                    return Err(FactorizedError::Ports {
                        detail: format!(
                            "{} is not a Lorentz or spinor slot",
                            port.to_canonical_string()
                        ),
                    });
                }
            }
        }
        let projector_free = analysed
            .iter()
            .find(|(location, _, _)| *location == Location::Projector)
            .map(|(_, sets, _)| sets.free.keys().cloned().collect::<BTreeSet<_>>());
        if let Some(projector_free) = &projector_free
            && *projector_free != ports
        {
            return Err(FactorizedError::Ports {
                detail: format!(
                    "the projector's free indices {projector_free:?} are not the ports {ports:?}"
                ),
            });
        }
        for (key, locations) in &free {
            let in_factors = locations
                .iter()
                .filter(|l| **l != Location::Projector)
                .count();
            let is_port = ports.contains(key);
            let expected = if is_port { 1 } else { 2 };
            if in_factors > expected {
                return Err(FactorizedError::IndexOverContracted {
                    slot: key.clone(),
                    count: in_factors,
                    location: locations
                        .iter()
                        .map(|l| l.to_string())
                        .collect::<Vec<_>>()
                        .join(", "),
                });
            }
            if in_factors < expected {
                if is_port {
                    return Err(FactorizedError::Ports {
                        detail: format!("port {key} does not occur in the factors"),
                    });
                }
                return Err(FactorizedError::IndexUncontracted {
                    slot: key.clone(),
                    location: locations[0].to_string(),
                });
            }
        }
        for port in &ports {
            if !free.contains_key(port) {
                return Err(FactorizedError::Ports {
                    detail: format!("port {port} does not occur in the factors"),
                });
            }
        }
        Ok(())
    }

    /// `gamma(row, column, μ)`: at an edge the row is its target half-edge
    /// and the column its source half-edge; at a vertex the row is the source
    /// half-edge of the outgoing fermion line and the column the target
    /// half-edge of the incoming one (spinor matrices multiply against the
    /// fermion flow). Free spinor indices of a fermion edge are its own two.
    fn check_fermion_orientation(
        &self,
        analysed: &[(Location, IndexSets, Elements)],
    ) -> Result<()> {
        for (location, sets, elements) in analysed {
            let (rows, columns) = match location {
                Location::Edge(id) => {
                    let edge = self.edge(*id).unwrap();
                    let spinor_free = sets
                        .free
                        .values()
                        .filter(|s| s.rep == names::BIS)
                        .filter_map(Slot::hedge)
                        .collect::<BTreeSet<_>>();
                    let own = edge
                        .source_hedge
                        .into_iter()
                        .chain(edge.target_hedge)
                        .collect::<BTreeSet<_>>();
                    if edge.line == LineKind::Fermion {
                        if spinor_free != own {
                            return Err(FactorizedError::FermionOrientation {
                                location: location.to_string(),
                                detail: format!(
                                    "free spinor half-edges {spinor_free:?}, expected {own:?}"
                                ),
                            });
                        }
                    } else if !spinor_free.is_empty() {
                        return Err(FactorizedError::FermionOrientation {
                            location: location.to_string(),
                            detail: format!("a {} line carries spinor indices", edge.line.as_str()),
                        });
                    }
                    (
                        edge.target_hedge.into_iter().collect::<BTreeSet<_>>(),
                        edge.source_hedge.into_iter().collect::<BTreeSet<_>>(),
                    )
                }
                Location::Vertex(v) => {
                    let fermions = self.edges.iter().filter(|e| e.line == LineKind::Fermion);
                    let mut rows = BTreeSet::new();
                    let mut columns = BTreeSet::new();
                    for e in fermions {
                        if e.source == Endpoint::Vertex(*v) {
                            rows.extend(e.source_hedge);
                        }
                        if e.target == Endpoint::Vertex(*v) {
                            columns.extend(e.target_hedge);
                        }
                    }
                    (rows, columns)
                }
                Location::Projector => continue,
            };
            let local = self.local_hedges(*location);
            for (row, column) in &elements.gammas {
                if let Some(h) = row.hedge().filter(|h| local.contains(h))
                    && !rows.contains(&h)
                {
                    return Err(FactorizedError::FermionOrientation {
                        location: location.to_string(),
                        detail: format!(
                            "gamma row index {} is not the outgoing side of a fermion line",
                            row.key()
                        ),
                    });
                }
                if let Some(h) = column.hedge().filter(|h| local.contains(h))
                    && !columns.contains(&h)
                {
                    return Err(FactorizedError::FermionOrientation {
                        location: location.to_string(),
                        detail: format!(
                            "gamma column index {} is not the incoming side of a fermion line",
                            column.key()
                        ),
                    });
                }
            }
        }
        Ok(())
    }

    /// The closed fermion and ghost loops, each in flow order, sorted by
    /// their smallest edge id. Errors on open lines and on vertices where a
    /// flow is not conserved or branches.
    pub fn closed_loops(&self) -> Result<Vec<ClosedLoop>> {
        let mut loops = Vec::new();
        for line in [LineKind::Fermion, LineKind::Ghost] {
            let lines = self
                .edges
                .iter()
                .filter(|e| e.line == line)
                .collect::<Vec<_>>();
            if !line.is_directed_loop() {
                continue;
            }
            for edge in &lines {
                if edge.is_external() {
                    return Err(FactorizedError::OpenLine {
                        edge: edge.id,
                        line: line.as_str().into(),
                    });
                }
            }
            for vertex in &self.vertices {
                let incoming = lines
                    .iter()
                    .filter(|e| e.target == Endpoint::Vertex(vertex.id))
                    .count();
                let outgoing = lines
                    .iter()
                    .filter(|e| e.source == Endpoint::Vertex(vertex.id))
                    .count();
                if incoming != outgoing || incoming > 1 {
                    return Err(FactorizedError::FlowNotConserved {
                        vertex: vertex.id,
                        line: line.as_str().into(),
                        incoming,
                        outgoing,
                    });
                }
            }
            let mut unused = lines.iter().map(|e| e.id).collect::<BTreeSet<_>>();
            while let Some(start) = unused.pop_first() {
                let mut cycle = vec![start];
                let mut current = lines.iter().find(|e| e.id == start).unwrap();
                loop {
                    let next = lines
                        .iter()
                        .find(|e| e.source == current.target)
                        .expect("conserved flow");
                    if next.id == start {
                        break;
                    }
                    unused.remove(&next.id);
                    cycle.push(next.id);
                    current = next;
                }
                loops.push(ClosedLoop { line, edges: cycle });
            }
        }
        loops.sort_by_key(|l| l.edges.iter().min().copied());
        Ok(loops)
    }

    /// The versioned text form ([`FACTORIZED_SCHEMA`]).
    pub fn to_text(&self) -> Result<String> {
        text::write_record(self)
    }

    /// Read one record of the text form and [`Self::validate`] it.
    pub fn from_text(input: &str) -> Result<Self> {
        let value = text::parse_record(input)?;
        value.validate()?;
        Ok(value)
    }

    /// Read one record of the text form without validating it.
    pub fn parse_text(input: &str) -> Result<Self> {
        text::parse_record(input)
    }
}

#[cfg(test)]
mod tests;
