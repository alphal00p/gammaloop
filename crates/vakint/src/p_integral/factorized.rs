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
//! The R* closure kind ([`FactorizedKind::Closure`]: a massive closing line
//! and an unapplied [`TaylorSpec`]) is evaluated by
//! [`evaluate_factorized_closure`] (module `closure`, text schema
//! [`FACTORIZED_SCHEMA_V3`]). Counterterm insertion vertices ([`Insertion`])
//! exist in schema v3 only: a v1 text with `insertions` other than 0 is
//! `Reserved`, and a value with insertions is written as v3.
//!
//! A value WITHOUT a new feature (kind closure, insertions, a v3 direct
//! record: [`FactorizedNumerator::has_new_features`]) is validated and
//! evaluated by exactly the rules it had before the R* record kinds: every
//! pre-existing entry runs the earlier body for it. The widened grammar of
//! the record kinds applies to marked values only. The existing
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

mod closure;
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

pub use closure::{
    ClosureRecordValue, ClosureValue, MAX_CLOSURE_TAYLOR_DEGREE, evaluate_factorized_closure,
};
pub use colour::{
    ColourMode, ColourRow, MergedColourTerm, colour_coefficient, merge_colour_channels,
    split_local_colour,
};
pub use evaluate::{
    FactorizedValue, MAX_OPEN_PORTS, PortsValue, contract_to_vakint_dots, evaluate_factorized,
    evaluate_factorized_ports, to_caller_notation_p, vakint_input,
};
#[cfg(feature = "feynkit-ingress")]
pub use feynkit::{
    FromFeynkitOptions, OverallFactorSource, ProjectorChoice, dot_sets_graph_attribute,
};
#[cfg(feature = "feynkit-ingress")]
pub use graph::{GraphIntegral, GraphRoute, GraphValue, evaluate_graph, evaluate_graph_adapter};
pub use oracle::contract_with_idenso;
pub use text::{parse_factorized_records, read_factorized_records};
pub use vacuum::{VacuumValue, evaluate_factorized_vacuum};

/// Schema line of the text form (first line of every record).
pub const FACTORIZED_SCHEMA: &str = "pint6.factorized-numerator.v1";

/// Schema line of the v3 text form: v1 plus the R* lines (closure section,
/// `open`, `softs`, `original_cycles`, `weights`, `depth`). Records of kind
/// closure are written and read in it; direct records stay v1.
pub const FACTORIZED_SCHEMA_V3: &str = "pint6.factorized-numerator.v3";

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
    /// idenso's scalar product of two compact momenta `dot(Q(e, mink(D)), Q(f, mink(D)))`.
    pub const DOT: &str = "spenso::dot";
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
    /// An R* remainder (or one pre-expanded Taylor sector of one) closed by a
    /// massive line among massless ones, evaluated by
    /// [`evaluate_factorized_closure`]; the record form of
    /// N/notes/wh/WH-RECORD-FORMAT.md (kinds K1-K3). The adapter and the vacuum
    /// entry refuse it.
    Closure(ClosureSpec),
}

/// The closing line of an R* remainder and its Taylor step.
#[derive(Clone, Debug, PartialEq)]
pub struct ClosureSpec {
    /// Id of the closing (marked) edge.
    pub edge: usize,
    /// The closing edge's mass atom as the record writes it: the REFERENCE
    /// symbol `R` (equal to that edge's `mass_squared`).
    pub mass_squared: Atom,
    /// `M_aux²`: with a Taylor step, `R -> (1 - λ²) M_aux²` (compensate) or
    /// `R -> M_aux²` (fixed) everywhere in the integrand; without one (a
    /// pre-expanded sector) it equals `mass_squared`.
    pub auxiliary_mass_squared: Atom,
    /// Power of the closing line (equal to that edge's `power`).
    pub power: i64,
    /// `externals=independent`: every crown leg carries its own `p(j)`;
    /// momentum conservation is checked on the loop momenta only.
    pub independent_externals: bool,
    /// `None` = a pre-expanded Taylor sector (`denominators=lambda0`): every
    /// internal edge carries its `λ = 0` momentum.
    pub taylor: Option<TaylorSpec>,
    /// Informational provenance of the crown legs (`soft` lines).
    pub softs: Vec<SoftLeg>,
    /// `(fermion, ghost)` closed cycles of the ORIGINAL graph (provenance).
    pub original_cycles: Option<(u32, u32)>,
    /// Laurent terms wanted (vakint `number_of_terms_in_epsilon_expansion`);
    /// `None`: the caller's settings.
    pub depth: Option<usize>,
    /// Optional weight provenance (`weight<TAB><name><TAB>applied_by=..`):
    /// `(name, applied by vakint, atom)`; an `applied_by=vakint` weight must
    /// equal `prefactor_scalars`.
    pub weights: Vec<(String, bool, Atom)>,
}

/// The unapplied Taylor step of an R* remainder: the momentum `p(j)` of every
/// external leg in `scaled_edges` becomes `λ p(j)` (in every edge momentum and
/// every `Q(x, ..)`), the reference mass `R` becomes `(1 − λ²) M_aux²`
/// (`compensate_mass`) or `M_aux²`, the product of factors and denominators
/// (not the projector) is expanded to total degree `degree` in `λ`, the
/// `orders` are kept, then `λ = 1`; the projector and `(p²)^p2_power` act
/// after the step.
#[derive(Clone, Debug, PartialEq)]
pub struct TaylorSpec {
    pub degree: u32,
    /// Edge ids of the external legs whose momenta are scaled (the R*
    /// `root_edges`); together they carry every independent external momentum.
    pub scaled_edges: Vec<usize>,
    pub compensate_mass: bool,
    pub orders: TaylorOrders,
}

/// Which orders in `λ` of the Taylor polynomial are kept.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TaylorOrders {
    /// `0..=degree` (`orders=all`).
    All,
    /// `1..=degree` (`orders=drop0`).
    DropZero,
    /// The even orders `2..=degree` (`orders=drop0-odd`).
    DropZeroOdd,
}

impl TaylorOrders {
    /// Whether order `n` (`<= degree`) is kept.
    pub fn keeps(self, n: u32) -> bool {
        match self {
            Self::All => true,
            Self::DropZero => n >= 1,
            Self::DropZeroOdd => n >= 2 && n.is_multiple_of(2),
        }
    }

    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::DropZero => "drop0",
            Self::DropZeroOdd => "drop0-odd",
        }
    }

    pub(crate) fn parse(text: &str) -> Option<Self> {
        Some(match text {
            "all" => Self::All,
            "drop0" => Self::DropZero,
            "drop0-odd" => Self::DropZeroOdd,
            _ => return None,
        })
    }
}

/// One crown leg of a closure record (`soft<TAB><j><TAB>edge=<x><TAB>role=..`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SoftLeg {
    /// One-based index of its external momentum `p(j)`.
    pub external: usize,
    pub edge: usize,
    /// `true`: `role=root`; `false`: `role=soft`.
    pub root: bool,
}

/// A counterterm INSERTION: an integrated subgraph acting as a local tensor
/// at one spot of the graph, `coefficient × factor`, multiplied into the
/// product of factors like any vertex factor (its boundary slots contracted
/// with the adjacent edges' slots by idenso).
///
/// Two producers have this shape: an R* child counterterm Δ on the vertex its
/// support was contracted to (vertex rule `ct:<label>`; `hedges` = exactly the
/// vertex's incident half-edges), and the coefficient vertex of an
/// insertion surrogate (a proper two-point region replaced by a chain of
/// scalar bubbles; `hedges` = the two inner boundary half-edges, one at each
/// end of the chain, so the factor's slots may sit at two vertices).
///
/// Rules ([`FactorizedNumerator::validate`]): the vertex exists, carries at
/// most one insertion and has vertex factor `1` (the factor is stated once,
/// here); `hedges` are distinct half-edges of the graph, at least one at the
/// vertex; every free slot of `factor` is a slot of one of `hedges`; every
/// momentum `Q(e, ..)` in it belongs to an edge owning one of `hedges` (or is
/// a constant external vector absent from the graph); `coefficient` is a
/// scalar; neither contains colour objects, a flavour symbol (`…::nf`) or
/// `mUVexp`. Both may depend on the dimension, on the evaluation's ε symbol
/// and on any other scalar (couplings, the reference mass squared). For the
/// closed-loop weights the vertices owning `hedges` are one node (a fermion
/// line passes through a surrogate) and `cycles` (closed fermion, closed
/// ghost cycles INSIDE the inserted region) add to those of the carried
/// graph: `loop_sign` / `nf_power` are the whole original diagram's.
#[derive(Clone, Debug, PartialEq)]
pub struct Insertion {
    pub vertex: usize,
    /// Boundary half-edges whose slots the factor carries, sorted.
    pub hedges: Vec<usize>,
    /// (closed fermion cycles, closed ghost cycles) inside the inserted region.
    pub cycles: (u32, u32),
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
    /// A DIRECT record read from (and written as) schema v3: its provenance
    /// lines; `None` for a v1 record and for every closure record (whose
    /// provenance is in its [`ClosureSpec`]).
    pub record: Option<DirectRecord>,
}

/// The v3 provenance of a DIRECT record (`closure none`, `taylor none`,
/// `softs 0` in its closure section): with `original_cycles` the weights
/// (loop sign, n_f, overall factor) are the caller's, applied once after the
/// call (WH-RECORD-FORMAT, the weights rule R2); without, the v1 rule.
#[derive(Clone, Debug, PartialEq)]
pub struct DirectRecord {
    /// `(fermion, ghost)` closed cycles of the ORIGINAL graph.
    pub original_cycles: Option<(u32, u32)>,
    /// Laurent terms wanted; `None`: the caller's settings.
    pub depth: Option<usize>,
    /// Weight provenance `(name, applied by vakint, atom)`.
    pub weights: Vec<(String, bool, Atom)>,
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
    #[error("insertion at vertex {vertex}: {detail}")]
    Insertion { vertex: usize, detail: String },
    #[error("{feature} is reserved for the R* path and not accepted by v1")]
    Reserved { feature: String },
    #[error("closure: {detail}")]
    Closure { detail: String },
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
    /// The factor of the insertion at this vertex.
    Insertion(usize),
    Projector,
}

impl std::fmt::Display for Location {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Edge(id) => write!(f, "edge {id}"),
            Self::Vertex(id) => write!(f, "vertex {id}"),
            Self::Insertion(id) => write!(f, "insertion at vertex {id}"),
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

/// Rule R4 of the record format: the COMPACT slash `gamma(a, b, Q(x,
/// mink(D)))` (a gamma matrix contracted with a momentum, as the R* sites
/// write the quark closure) written out as `gamma(a, b, mink(D, s)) Q(x,
/// mink(D, s))` with a fresh dummy index `s = vakint::slash(scope, n)`, `n =
/// 0, 1, ..` in the order of the expression tree. `scope` keeps the dummies
/// of different factors apart (the contract wants a dummy local to one
/// factor): the projector uses -1, an insertion factor its vertex id.
/// Everything else is returned unchanged. A substitution, no algebra.
pub fn expand_compact_slash(expression: &Atom, scope: i64) -> Atom {
    let counter = std::cell::Cell::new(0_i64);
    let slash = symbolica::symbol!("vakint::slash");
    expression.replace_map(|part, _, out| {
        let AtomView::Fun(gamma) = part else {
            return;
        };
        if gamma.get_symbol().get_name() != names::GAMMA || gamma.get_nargs() != 3 {
            return;
        }
        let args = gamma.iter().collect::<Vec<_>>();
        let AtomView::Fun(momentum) = args[2] else {
            return;
        };
        if momentum.get_symbol().get_name() != names::MOMENTUM || momentum.get_nargs() != 2 {
            return;
        }
        let margs = momentum.iter().collect::<Vec<_>>();
        let AtomView::Fun(rep) = margs[1] else {
            return;
        };
        if rep.get_symbol().get_name() != names::MINK || rep.get_nargs() != 1 {
            return;
        }
        let n = counter.get();
        counter.set(n + 1);
        let dimension = rep.iter().next().unwrap().to_owned();
        let index = symbolica::function!(slash, Atom::num(scope), Atom::num(n));
        let slot = symbolica::function!(rep.get_symbol(), dimension, index);
        **out = symbolica::function!(
            gamma.get_symbol(),
            args[0].to_owned(),
            args[1].to_owned(),
            &slot
        ) * symbolica::function!(momentum.get_symbol(), margs[0].to_owned(), &slot);
    })
}

/// The one symbol of vakint's namespace a record of the R* interface may
/// carry: vakint's default epsilon (`VakintSettings::epsilon_symbol`), used by
/// internal code only as the evaluation's epsilon.
pub const ALLOWED_VAKINT_SYMBOL: &str = "vakint::ε";

/// What [`first_forbidden_symbol`] found.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ForbiddenSymbol {
    /// A symbol of the reserved `vakint::` namespace.
    Reserved(String),
    /// A symbol whose name ends in `_`: a symbolica pattern WILDCARD wherever
    /// the entries use a caller atom as a replacement pattern (the dimension,
    /// the reference mass, ...), where it would match any subexpression.
    Wildcard(String),
}

/// The first FORBIDDEN symbol (variable or function head, at any nesting
/// level) of `view`: one of vakint's namespace other than the names in
/// `exempt` (callers pass [`ALLOWED_VAKINT_SYMBOL`] when vakint's epsilon is
/// the evaluation's), or one whose name ends in `_`. The dummy index vakint's reader writes when it opens a compact
/// slash, `mink(D, vakint::slash(scope, n))`, is accepted in exactly that slot
/// form (a forged one can only clash with a reader dummy, which the
/// contraction rule refuses as an index used more than twice or shared
/// between factors; no code substitutes it).
pub(crate) fn first_forbidden_symbol(
    view: AtomView<'_>,
    exempt: &[&str],
) -> Option<ForbiddenSymbol> {
    let classify = |name: &str| -> Option<ForbiddenSymbol> {
        if name.ends_with('_') {
            Some(ForbiddenSymbol::Wildcard(name.to_string()))
        } else if name.starts_with("vakint::") && !exempt.contains(&name) {
            Some(ForbiddenSymbol::Reserved(name.to_string()))
        } else {
            None
        }
    };
    let mut found = None;
    view.visitor(&mut |part| {
        if found.is_some() {
            return false;
        }
        match part {
            AtomView::Var(v) => found = classify(v.get_symbol().get_name()),
            AtomView::Fun(f) => {
                let symbol = f.get_symbol();
                let name = symbol.get_name();
                if name == names::MINK && f.get_nargs() == 2 {
                    let args = f.iter().collect::<Vec<_>>();
                    if let AtomView::Fun(index) = args[1]
                        && index.get_symbol().get_name() == "vakint::slash"
                        && index.iter().all(|a| matches!(a, AtomView::Num(_)))
                    {
                        found = first_forbidden_symbol(args[0], exempt);
                        return false;
                    }
                }
                found = classify(name);
            }
            _ => {}
        }
        found.is_none()
    });
    found
}

/// Where a checked atom comes from (it decides the refusal's advice).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum AtomSource {
    /// A field of the record (written by the reader or the caller): a bare
    /// name in a record text is read into `vakint::`.
    Record,
    /// An argument of the entry (the prefactor argument, a setting).
    Argument,
}

/// The refusal of a [`ForbiddenSymbol`] in `field` (variant `Graph`): in a
/// record field `reserved symbol <name> in <field>: the vakint:: namespace is
/// reserved for vakint's internal symbols; a name written without a
/// namespace is read into vakint::` (true whether the caller wrote
/// `vakint::x` or a bare `x`: the atom cannot tell them apart), in an
/// argument `reserved symbol <name> in <field>: the vakint:: namespace is
/// reserved for vakint's internal symbols`; or `a symbol ending in _ cannot
/// be used as <field>: <name> is a pattern wildcard`.
pub(crate) fn forbidden_symbol_error(
    found: ForbiddenSymbol,
    field: &str,
    source: AtomSource,
) -> FactorizedError {
    let detail = match (found, source) {
        (ForbiddenSymbol::Reserved(name), AtomSource::Record) => format!(
            "reserved symbol {name} in {field}: the vakint:: namespace is reserved for vakint's internal symbols; a name written without a namespace is read into vakint::"
        ),
        (ForbiddenSymbol::Reserved(name), AtomSource::Argument) => format!(
            "reserved symbol {name} in {field}: the vakint:: namespace is reserved for vakint's internal symbols"
        ),
        (ForbiddenSymbol::Wildcard(name), _) => {
            format!("a symbol ending in _ cannot be used as {field}: {name} is a pattern wildcard")
        }
    };
    FactorizedError::Graph { detail }
}

/// THE caller-symbol refusal of a marked value's atoms (one function for the
/// record walk [`FactorizedNumerator::check_reserved_symbols`] and the
/// entries' [`check_prefactor_argument`]): no symbol of the reserved
/// `vakint::` namespace (but `exempt`), no wildcard name (ending in `_`).
pub(crate) fn check_caller_atom(
    atom: &Atom,
    field: &str,
    exempt: &[&str],
    source: AtomSource,
) -> Result<()> {
    match first_forbidden_symbol(atom.as_view(), exempt) {
        Some(found) => Err(forbidden_symbol_error(found, field, source)),
        None => Ok(()),
    }
}

/// THE PREFACTOR ARGUMENT of a marked value (the atom an entry multiplies in
/// besides the record, normally `value.prefactor.total()`): the same checks as
/// the record's prefactor fields, called by every entry that takes one
/// (`evaluate_factorized`, hence `evaluate_factorized_ports`,
/// `evaluate_factorized_vacuum[_streamed]`, `evaluate_factorized_closure`),
/// in this order:
/// - caller symbols ([`check_caller_atom`], variant `Graph`): no reserved
///   `vakint::` symbol (vakint's internal placeholders, the Taylor parameter
///   `vakint::closure_lambda`, the notation heads `p`, `k`, `g`, `dot`,
///   `prop`, `topo`, at any nesting, function heads too) and no wildcard name;
///   `vakint::nf` is exempt here only, being the flavour symbol
///   [`Prefactor::total`] writes (no internal code substitutes it), and
///   `vakint::ε` only when it is the settings' epsilon `epsilon`;
/// - a scalar (`InexactNumber` / `NotScalar`, location `prefactor argument`):
///   exact numbers, no index (no open colour index either), no momentum (no
///   loop momentum); colour invariants as in `prefactor_colour` are allowed;
/// - (`NotScalar`) a flavour symbol (`nf`, `…::nf`) only when the record's
///   own `nf_power` is nonzero (the factor `total` writes; otherwise n_f
///   stays with the caller), no `mUVexp`;
/// - (`NotScalar`, closure kind) no reference mass `R` (substituted in the
///   integrand only, as in `prefactor_scalars` / `prefactor_colour`).
///
/// Unmarked values are not checked (the base behaviour).
pub(crate) fn check_prefactor_argument(
    value: &FactorizedNumerator,
    prefactor: &Atom,
    epsilon: Symbol,
) -> Result<()> {
    if !is_marked(value) {
        return Ok(());
    }
    const FIELD: &str = "prefactor argument";
    let exempt: &[&str] = if epsilon.get_name() == ALLOWED_VAKINT_SYMBOL {
        &[names::NF, ALLOWED_VAKINT_SYMBOL]
    } else {
        &[names::NF]
    };
    check_caller_atom(prefactor, FIELD, exempt, AtomSource::Argument)?;
    check_scalar(prefactor.as_view(), FIELD)?;
    let refuse = |detail: String| FactorizedError::NotScalar {
        location: FIELD.to_string(),
        detail,
    };
    let mut flavour = None;
    let mut exp = None;
    prefactor.as_view().visitor(&mut |part| {
        let symbol = match part {
            AtomView::Var(v) => v.get_symbol(),
            AtomView::Fun(f) => f.get_symbol(),
            _ => return true,
        };
        let name = symbol.get_name();
        if name == "nf" || name.ends_with("::nf") {
            flavour.get_or_insert_with(|| name.to_string());
        }
        if name.ends_with("mUVexp") {
            exp.get_or_insert_with(|| name.to_string());
        }
        true
    });
    if let Some(name) = flavour
        && value.prefactor.nf_power == 0
    {
        return Err(refuse(format!(
            "the flavour symbol {name} with the record's nf_power=0 (n_f stays with the caller)"
        )));
    }
    if let Some(name) = exp {
        return Err(refuse(format!("{name} (the producer maps it to mUV)")));
    }
    if let FactorizedKind::Closure(spec) = &value.kind
        && closure::contains_atom(prefactor.as_view(), spec.mass_squared.as_view())
    {
        return Err(refuse(format!(
            "the reference mass {} (it is substituted in the integrand only)",
            spec.mass_squared
        )));
    }
    Ok(())
}

/// The settings' epsilon of a marked evaluation: not a wildcard, not a
/// reserved `vakint::` name but [`ALLOWED_VAKINT_SYMBOL`] (it is the
/// replacement of `D -> 4 - 2 epsilon` and the expansion variable). When it
/// is NOT vakint's epsilon, the record may not carry `vakint::ε` either
/// (validate exempts it, not knowing the settings: there it would be a
/// silent constant): the record is walked once more without the exemption.
pub(crate) fn check_epsilon_setting(value: &FactorizedNumerator, epsilon: Symbol) -> Result<()> {
    if !is_marked(value) {
        return Ok(());
    }
    check_caller_atom(
        &Atom::var(epsilon),
        "the epsilon_symbol setting",
        &[ALLOWED_VAKINT_SYMBOL],
        AtomSource::Argument,
    )?;
    if epsilon.get_name() != ALLOWED_VAKINT_SYMBOL {
        value
            .check_reserved_symbols_except(&[])
            .map_err(|error| match error {
                FactorizedError::Graph { detail }
                    if detail.starts_with(&format!("reserved symbol {ALLOWED_VAKINT_SYMBOL} in ")) =>
                {
                    let field = detail
                        .trim_start_matches(&format!("reserved symbol {ALLOWED_VAKINT_SYMBOL} in "))
                        .split(": ")
                        .next()
                        .unwrap_or_default()
                        .to_string();
                    FactorizedError::Graph {
                        detail: format!(
                            "reserved symbol {ALLOWED_VAKINT_SYMBOL} in {field}: the settings' epsilon is {}, so vakint's epsilon would be a constant here (write the settings' epsilon)",
                            epsilon.get_name()
                        ),
                    }
                }
                other => other,
            })?;
    }
    Ok(())
}

/// THE ENTRY CHECK of a marked value ([`check_prefactor_argument`] then
/// [`check_epsilon_setting`]), called by every public function that turns a
/// marked value and a prefactor into vakint's input: the closure entry,
/// `vakint_input` (hence `evaluate_factorized` and the ports entry),
/// `vakint_vacuum_input_streamed` (hence `vakint_vacuum_input` and the
/// vacuum entries), and at the top of `evaluate_factorized` and the vacuum
/// entry. A no-op for an unmarked value.
pub(crate) fn check_marked_entry(
    value: &FactorizedNumerator,
    prefactor: &Atom,
    epsilon: Symbol,
    terms: Option<i64>,
) -> Result<()> {
    check_prefactor_argument(value, prefactor, epsilon)?;
    check_epsilon_setting(value, epsilon)?;
    // the settings' number of terms (an entry that does not take it from the
    // record), by the depth rule of [`MARKED_MAX_TERMS`]; marked values only
    // (the vacuum entries call this for every value)
    if let Some(terms) = terms.filter(|_| is_marked(value)) {
        marked_depth(i128::from(terms)).map_err(|detail| FactorizedError::Graph { detail })?;
    }
    Ok(())
}

/// Whether `value` is a MARKED value (a record kind of the R* interface:
/// kind closure, insertions, or a v3 direct record), to which the reserved
/// namespace applies.
pub(crate) fn is_marked(value: &FactorizedNumerator) -> bool {
    matches!(value.kind, FactorizedKind::Closure(_)) || evaluate::is_new_kind(value)
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

/// Limits of an insertion factor (R* child counterterms of QCD: at most four
/// boundary slots, degree at most two in the boundary momenta).
pub const INSERTION_MAX_SLOTS: usize = 4;
pub const INSERTION_MAX_DEGREE: u64 = 2;

/// An upper bound of the degree of `view` in the momenta `Q(..)`: the largest
/// number of momentum factors of a term (a sum takes the maximum, a product
/// the sum, an integer power multiplies). Not expanded: cancellations are not
/// seen, so the bound may only be too large for a polynomial written with
/// cancelling terms.
fn momentum_degree(view: AtomView<'_>) -> u64 {
    match view {
        AtomView::Fun(f) if f.get_symbol().get_name() == names::MOMENTUM => 1,
        AtomView::Fun(f) if f.get_symbol().get_name() == names::DOT => {
            f.iter().map(momentum_degree).fold(0, u64::saturating_add)
        }
        AtomView::Add(a) => a.iter().map(momentum_degree).max().unwrap_or(0),
        AtomView::Mul(m) => m.iter().map(momentum_degree).fold(0, u64::saturating_add),
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            let degree = momentum_degree(base);
            match i64::try_from(exponent) {
                Ok(n) if n >= 0 => degree.saturating_mul(n as u64),
                _ if degree == 0 => 0,
                _ => u64::MAX,
            }
        }
        _ => 0,
    }
}

/// The pole order's rules for one insertion, `depends` selecting ε and the
/// dimension: the dependence on them is a Laurent polynomial
/// ([`non_laurent_dependence`]) and its bound is at most
/// [`INSERTION_MAX_EPSILON_POWER`] (`an ε power out of range`).
fn insertion_epsilon_rules(insertion: &Insertion, depends: &dyn Fn(Symbol) -> bool) -> Result<()> {
    let v = insertion.vertex;
    let error = |detail: String| FactorizedError::Insertion { vertex: v, detail };
    // The pole order (`insertion_pole_order`) is that of coefficient x
    // factor as a Laurent expression in ε, the dimension set to
    // 4 - 2ε: ε and the dimension may occur only through sums,
    // products and integer powers (a representation slot's dimension
    // argument, mink(D, ..), is not a dependence). Inside a function's
    // argument, in an exponent or under a non-integer power the
    // expansion has no pole order (log ε would silently count as 0).
    for (what, atom) in [
        ("coefficient", &insertion.coefficient),
        ("factor", &insertion.factor),
    ] {
        if let Some(found) = non_laurent_dependence(atom.as_view(), depends) {
            return Err(error(format!(
                "the dependence on ε is not a Laurent polynomial: {found} in the {what}"
            )));
        }
    }
    // bounded exponents: the degree count of coefficient x factor (not the
    // pole order, which insertion_pole_order computes exactly)
    let bound = epsilon_power_bound(insertion.coefficient.as_view(), depends)
        .and_then(|c| epsilon_power_bound(insertion.factor.as_view(), depends)?.checked_add(c));
    match bound {
        Some(b) if b <= INSERTION_MAX_EPSILON_POWER => {}
        _ => {
            return Err(error(format!(
                "an ε power out of range: the powers of ε and of the dimension in coefficient x factor reach {}, at most {INSERTION_MAX_EPSILON_POWER} are supported",
                bound.map_or_else(|| "beyond 2^64".to_string(), |b| b.to_string())
            )));
        }
    }
    Ok(())
}

/// The cap on the DEGREE COUNT of an insertion's coefficient x factor in ε
/// (the dimension outside representation slots counted as 4 - 2ε; see
/// [`epsilon_power_bound`], the coefficient's and the factor's counts added):
/// 64. The count is NOT the pole order: it bounds the size of the exponents
/// the series code is handed (no exponent beyond i64, no overflow in the
/// count). It can refuse an insertion without a high pole ((1+ε)^-100: count
/// 100, pole order 0) and admit one whose pole order exceeds it when a
/// negative power acts on a sum whose terms cancel at ε = 0
/// (((1-ε)^-1 - 1 - ε)^-40: count 40, pole order 80). The exact pole order is
/// computed by [`FactorizedNumerator::insertion_pole_order`] and used, and
/// the depth raise it causes is capped by [`INSERTION_MAX_TOTAL_POLE_ORDER`];
/// for such an unsimplified cancelling power symbolica's series is slow
/// (about 1 s at ^-10, 50 s at ^-20, not finished in minutes at ^-40).
pub(crate) const INSERTION_MAX_EPSILON_POWER: u64 = 64;

/// The cap on the SUM of the exact pole orders of a record's insertions, the
/// closure entry's depth raise (the only entry that raises a depth by it: the
/// direct, ports and vacuum entries hand the insertions to vakint inside the
/// integrand): 64. Real records stay far below: the 83 insertions of the R*
/// producer's fixture sets have totals of at most 3 per record (the five-loop
/// k1-l5 records), with at most 4 insertions per record (k5-stress). An R*
/// record's children are disjoint subgraphs and an l-loop child's
/// counterterm has poles up to ε^-l, so the total is at most L - 1, 5 at six
/// loops; even the producer's largest six-loop pole part, ε^-6, on each of 4
/// insertions gives 24. Above it the record is refused (`an ε power out of
/// range`) instead of being evaluated some 64 or more terms deeper.
pub(crate) const INSERTION_MAX_TOTAL_POLE_ORDER: u64 = 64;

/// An upper bound on |the power of a symbol `depends` selects| in any term of
/// `view` expanded (sums: the largest; products: the sum; an integer power:
/// |n| times its base's), with a representation slot's dimension argument
/// and, after the Laurent rule, any function counting 0; None beyond u64 or
/// for an exponent beyond i64. Not a bound on the pole order (see
/// [`INSERTION_MAX_EPSILON_POWER`]): a negative power of a sum whose terms
/// cancel at ε = 0 has a higher one.
fn epsilon_power_bound(view: AtomView<'_>, depends: &dyn Fn(Symbol) -> bool) -> Option<u64> {
    match view {
        AtomView::Var(v) => Some(u64::from(depends(v.get_symbol()))),
        AtomView::Add(a) => a
            .iter()
            .try_fold(0u64, |m, t| Some(m.max(epsilon_power_bound(t, depends)?))),
        AtomView::Mul(m) => m
            .iter()
            .try_fold(0u64, |s, t| s.checked_add(epsilon_power_bound(t, depends)?)),
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            let b = epsilon_power_bound(base, depends)?;
            if b == 0 {
                return Some(0);
            }
            b.checked_mul(i64::try_from(exponent).ok()?.unsigned_abs())
        }
        _ => Some(0),
    }
}

/// The most ε terms a marked value may ask for: its depth (the record's
/// `depth`, or the settings' number of terms) is at least 1 and, with the
/// closure entry's insertion raise added, at most 128 (a depth of 64 plus the
/// raise cap [`INSERTION_MAX_TOTAL_POLE_ORDER`]). The R* sites ask for at
/// most 12 (six-loop K1 / K2 at depth 6, K3 with a quotient's retries). Below
/// 1 the entries returned an empty (zero) value; far above, an evaluation that
/// does not end; beyond i64 the depth wrapped.
pub(crate) const MARKED_MAX_TERMS: i64 = 128;

/// The depth rule of a marked value ([`MARKED_MAX_TERMS`]): `depth` (any
/// integer type, widened without loss) as the number of terms, or the detail
/// of the refusal.
pub(crate) fn marked_depth(depth: i128) -> std::result::Result<i64, String> {
    i64::try_from(depth)
        .ok()
        .filter(|d| (1..=MARKED_MAX_TERMS).contains(d))
        .ok_or_else(|| {
            format!(
                "the depth {depth} is out of range: a marked value asks for at least 1 and at most {MARKED_MAX_TERMS} ε terms"
            )
        })
}

/// A symbol named `ε` in any namespace (the evaluation's epsilon; which one is
/// checked against the settings when the input is built).
fn is_epsilon_symbol(symbol: Symbol) -> bool {
    let name = symbol.get_name();
    name == "ε" || name.ends_with("::ε")
}

/// Whether `view` contains a symbol `depends` selects anywhere except as the
/// dimension argument (the first) of a representation slot, `mink(D, ..)`,
/// `bis(D, ..)`, a colour representation (`dind(..)` is looked through).
fn depends_outside_slots(view: AtomView<'_>, depends: &dyn Fn(Symbol) -> bool) -> bool {
    match view {
        AtomView::Var(v) => depends(v.get_symbol()),
        AtomView::Fun(f) => {
            let symbol = f.get_symbol();
            let name = symbol.get_name();
            let slot = is_rep_name(name) && name != names::DIND;
            f.iter()
                .enumerate()
                .any(|(i, a)| !(slot && i == 0) && depends_outside_slots(a, depends))
        }
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            depends_outside_slots(base, depends) || depends_outside_slots(exponent, depends)
        }
        AtomView::Mul(m) => m.iter().any(|t| depends_outside_slots(t, depends)),
        AtomView::Add(a) => a.iter().any(|t| depends_outside_slots(t, depends)),
        _ => false,
    }
}

/// The first sub-expression of `view` that depends on a symbol `depends`
/// selects (outside representation slots' dimension arguments) other than
/// through sums, products and integer powers (of any sign and size, of
/// expressions that are themselves such): a function whose argument depends
/// on it, a power whose exponent does, or a non-integer or symbolic power of a
/// base that does; with the reason. None when the dependence is a Laurent
/// polynomial (or a ratio of polynomials) in it.
fn non_laurent_dependence(view: AtomView<'_>, depends: &dyn Fn(Symbol) -> bool) -> Option<String> {
    let contains = |part: AtomView<'_>| depends_outside_slots(part, depends);
    let shown = |why: &str| Some(format!("{} ({why})", view.to_owned().to_canonical_string()));
    match view {
        AtomView::Add(a) => a.iter().find_map(|t| non_laurent_dependence(t, depends)),
        AtomView::Mul(m) => m.iter().find_map(|t| non_laurent_dependence(t, depends)),
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            if contains(exponent) {
                shown("in an exponent")
            } else if !contains(base) {
                None
            } else if symbolica::domains::integer::Integer::try_from(exponent).is_ok() {
                non_laurent_dependence(base, depends)
            } else {
                shown("under a non-integer or symbolic power")
            }
        }
        AtomView::Fun(_) if contains(view) => shown("inside a function argument"),
        _ => None,
    }
}

/// What an insertion's coefficient or factor may not contain: a colour
/// object, a flavour symbol (`…::nf`) or the auxiliary `mUVexp`; the first
/// one found, by name.
fn forbidden_in_insertion(view: AtomView<'_>) -> Option<String> {
    let mut found = None;
    view.visitor(&mut |part| {
        if found.is_some() {
            return false;
        }
        let symbol = match part {
            AtomView::Fun(f) => f.get_symbol(),
            AtomView::Var(v) => v.get_symbol(),
            _ => return true,
        };
        let name = symbol.get_name();
        if matches!(
            name,
            names::COAD
                | names::COF
                | names::COS
                | names::DIND
                | names::STRUCTURE_F
                | names::GENERATOR_T
        ) {
            found = Some(format!(
                "the colour object {}",
                part.to_owned().to_canonical_string()
            ));
        } else if name == "nf" || name.ends_with("::nf") {
            found = Some(format!("the flavour symbol {name}"));
        } else if name.contains("mUVexp") {
            found = Some(format!("the auxiliary mass symbol {name}"));
        }
        found.is_none()
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
    /// For an insertion factor: the edges owning one of its boundary
    /// half-edges (its momenta must be theirs).
    boundary: Option<BTreeSet<usize>>,
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
                // In a projector: the caller's Dirac closure between a
                // bispinor port pair (its spinor slots must be ports,
                // `check_contraction`; its Lorentz slot a port or open).
                if args.len() != 3
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
                if let Some(boundary) = &self.boundary
                    && !boundary.contains(&edge)
                {
                    return Err(FactorizedError::UnsupportedTensor {
                        location: self.loc(),
                        tensor: format!("{} (not a boundary momentum of the insertion)", text()),
                    });
                }
                if self.location == Location::Projector && !known.is_external() {
                    return Err(FactorizedError::UnsupportedTensor {
                        location: self.loc(),
                        tensor: format!("{} (loop momentum in the projector)", text()),
                    });
                }
                Ok(slots)
            }
            names::DOT => {
                // A scalar product of two compact momenta `Q(e, mink(D))`
                // (the beta's insertion factors write p² so); no slot.
                if args.len() != 2 {
                    return Err(unsupported());
                }
                for arg in &args {
                    let AtomView::Fun(q) = arg else {
                        return Err(unsupported());
                    };
                    let parts = q.iter().collect::<Vec<_>>();
                    let compact = q.get_symbol().get_name() == names::MOMENTUM
                        && parts.len() == 2
                        && matches!(parts[1], AtomView::Fun(rep)
                            if rep.get_symbol().get_name() == names::MINK
                                && rep.get_nargs() == 1
                                && rep.iter().next() == Some(Atom::var(self.dimension).as_view()));
                    if !compact {
                        return Err(unsupported());
                    }
                    let Some(edge) = i64::try_from(parts[0])
                        .ok()
                        .and_then(|e| usize::try_from(e).ok())
                    else {
                        return Err(FactorizedError::UnknownEdge {
                            location: self.loc(),
                            edge: parts[0].to_owned().to_canonical_string(),
                        });
                    };
                    let Some(known) = self.edges.get(&edge) else {
                        continue;
                    };
                    if let Some(boundary) = &self.boundary
                        && !boundary.contains(&edge)
                    {
                        return Err(FactorizedError::UnsupportedTensor {
                            location: self.loc(),
                            tensor: format!(
                                "{} (not a boundary momentum of the insertion)",
                                text()
                            ),
                        });
                    }
                    if self.location == Location::Projector && !known.is_external() {
                        return Err(FactorizedError::UnsupportedTensor {
                            location: self.loc(),
                            tensor: format!("{} (loop momentum in the projector)", text()),
                        });
                    }
                }
                Ok(Vec::new())
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
                // a power above 2 of a tensor with free slots repeats each
                // slot n times: over-contracted, refused before the n copies
                if base_free.is_empty() {
                    return Ok(Vec::new());
                }
                if n > 2 {
                    return Err(FactorizedError::IndexOverContracted {
                        slot: base_free[0].key(),
                        count: usize::try_from(n).unwrap_or(usize::MAX),
                        location: self.loc(),
                    });
                }
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

/// The base's factor rules, verbatim (values without new features:
/// [`FactorizedNumerator::has_new_features`]).
struct FactorRulesV1<'a> {
    location: Location,
    dimension: Symbol,
    edges: &'a BTreeMap<usize, &'a FactorizedEdge>,
}

impl FactorRulesV1<'_> {
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
            record: None,
        };
        // from_parts builds a value without new features: the base's rules
        // (no compact-slash expansion, no widened grammar)
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

    /// The slots that stay OPEN in the value: without a projector every port
    /// (in the order of `ports`); with one (sorted by canonical string), the
    /// ports it does not close
    /// and its own free slots that are not ports (e.g. a fresh Lorentz index
    /// of a Dirac projector `gamma(bis_c, bis_r, mink(D, nu))`). Empty for a
    /// closed scalar value. Assumes a valid value ([`Self::validate`]).
    pub fn open_slots(&self) -> Result<Vec<Atom>> {
        let Some(projector) = &self.projector else {
            // every port, in the value's own order (the vacuum entry numbers
            // the open slots by position, as it numbered the ports)
            return Ok(self.ports.clone());
        };
        let edges = self.check_graph()?;
        let rules = FactorRules {
            location: Location::Projector,
            dimension: self.dimension,
            edges: &edges,
            boundary: None,
        };
        let (sets, _) = rules.analyse(&projector.tensor)?;
        let ports = self
            .ports
            .iter()
            .map(|port| (port.to_canonical_string(), port.clone()))
            .collect::<BTreeMap<_, _>>();
        let mut open = ports
            .iter()
            .filter(|(key, _)| !sets.free.contains_key(*key))
            .map(|(key, port)| (key.clone(), port.clone()))
            .collect::<BTreeMap<_, _>>();
        for (key, slot) in &sets.free {
            if !ports.contains_key(key) {
                open.insert(key.clone(), slot.atom.clone());
            }
        }
        Ok(open.into_values().collect())
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
        for insertion in &self.insertions {
            product *= &insertion.coefficient;
            product *= &insertion.factor;
        }
        product
    }

    /// Check the contract (CONTRACT-v1.md). Returns the first violation. A
    /// closure record is checked by the closure rules
    /// (`validate_closure`: the contract's graph, factor and index rules plus
    /// those of the closure kind).
    pub fn validate(&self) -> Result<()> {
        if !self.has_new_features() {
            return self.validate_v1();
        }
        // the record kinds of the R* interface may not name vakint's own
        // symbols (the entries' internal placeholders live there)
        if is_marked(self) {
            self.check_reserved_symbols()?;
        }
        let closure = matches!(self.kind, FactorizedKind::Closure(_));
        let in_kind = |detail: String| {
            if closure {
                FactorizedError::Closure { detail }
            } else {
                FactorizedError::Graph { detail }
            }
        };
        // the record's depth (text or typed), checked without a cast
        let depth = match &self.kind {
            FactorizedKind::Closure(spec) => spec.depth,
            _ => self.record.as_ref().and_then(|r| r.depth),
        };
        if let Some(depth) = depth {
            marked_depth(i128::try_from(depth).unwrap_or(i128::MAX)).map_err(in_kind)?;
        }
        // ids are written as Q(id, ..) and edge(source, target): i64 (the
        // closure kind has its own bound, with its own text, in
        // validate_closure)
        if !closure && let Some(detail) = self.marked_id_overflow() {
            return Err(in_kind(detail));
        }
        if let FactorizedKind::Closure(spec) = &self.kind {
            return self.validate_closure(spec);
        }
        self.validate_contract(None)
    }

    /// A marked value's ids must be i64 (they are written as Q(id, ..),
    /// edge(source, target), with fresh loop-vector ids max + 1 + i above the
    /// largest): the detail of the refusal when one is not.
    fn marked_id_overflow(&self) -> Option<String> {
        let fresh = self
            .loop_edges
            .iter()
            .rposition(|e| e.is_none())
            .map_or(0, |i| i + 1);
        let top = self
            .edges
            .iter()
            .flat_map(|e| {
                let vertex = |end: &Endpoint| match end {
                    Endpoint::Vertex(v) => Some(*v),
                    Endpoint::External => None,
                };
                [Some(e.id), vertex(&e.source), vertex(&e.target)]
                    .into_iter()
                    .flatten()
            })
            .chain(self.vertices.iter().map(|v| v.id))
            .chain(self.constant_vectors())
            .max()
            .unwrap_or(0);
        let limit = usize::try_from(i64::MAX).unwrap_or(usize::MAX);
        top.checked_add(fresh)
            .is_none_or(|last| last > limit)
            .then(|| {
                format!(
                    "id {top}: the edge and vertex ids (with fresh loop-vector ids up to {top} + {fresh} above them) must not exceed i64::MAX"
                )
            })
    }

    /// RESERVED SYMBOLS: no atom of the record (factors, masses, ports,
    /// projector, prefactor, insertions, the closure line, the weights, the
    /// dimension) may contain a symbol of vakint's namespace `vakint::`, where
    /// the entries create their internal placeholders and the vakint notation
    /// they map to (`closure_tensor_<n>`, `closure_lambda`, `slash`, `nf`,
    /// `p`, `k`, `dot`, `g`, `topo`, `prop`, `edge`, ...): a caller's symbol
    /// of that name would be substituted or read as one of them. The single
    /// exception is [`ALLOWED_VAKINT_SYMBOL`], vakint's default epsilon,
    /// which internal code only ever uses as the evaluation's epsilon. Nor a
    /// symbol whose name ends in `_` (a symbolica pattern wildcard: the
    /// dimension and the reference mass are replacement patterns). The first
    /// offending atom is reported ([`forbidden_symbol_error`]).
    pub(crate) fn check_reserved_symbols(&self) -> Result<()> {
        self.check_reserved_symbols_except(&[ALLOWED_VAKINT_SYMBOL])
    }

    /// [`Self::check_reserved_symbols`] with the exempt `vakint::` names
    /// `exempt` (none: also `vakint::ε`, for settings with another epsilon).
    pub(crate) fn check_reserved_symbols_except(&self, exempt: &[&str]) -> Result<()> {
        let check =
            |atom: &Atom, field: &str| check_caller_atom(atom, field, exempt, AtomSource::Record);
        check(&Atom::var(self.dimension), "dimension")?;
        for edge in &self.edges {
            check(&edge.factor, &format!("edge {} factor", edge.id))?;
            check(&edge.mass_squared, &format!("edge {} mass", edge.id))?;
        }
        for vertex in &self.vertices {
            check(&vertex.factor, &format!("vertex {} factor", vertex.id))?;
        }
        for port in &self.ports {
            check(port, "ports")?;
        }
        if let Some(projector) = &self.projector {
            check(&projector.tensor, "projector")?;
        }
        check(&self.prefactor.overall, "prefactor_overall")?;
        check(&self.prefactor.scalars, "prefactor_scalars")?;
        check(&self.prefactor.colour, "prefactor_colour")?;
        for insertion in &self.insertions {
            check(
                &insertion.coefficient,
                &format!("insertion {} coefficient", insertion.vertex),
            )?;
            check(
                &insertion.factor,
                &format!("insertion {} factor", insertion.vertex),
            )?;
        }
        let weights = match (&self.kind, &self.record) {
            (FactorizedKind::Closure(spec), _) => {
                check(&spec.mass_squared, "closure mass")?;
                check(&spec.auxiliary_mass_squared, "taylor auxiliary mass")?;
                spec.weights.as_slice()
            }
            (_, Some(record)) => record.weights.as_slice(),
            _ => &[],
        };
        for (name, _, atom) in weights {
            check(atom, &format!("weight {name}"))?;
        }
        Ok(())
    }

    /// The contract's checks. `closure = Some(conserved)`: a closure record,
    /// momentum conservation on the first `conserved` momentum components
    /// only, and no graph-derived loop-sign / n_f rule (its own rule).
    pub(crate) fn validate_contract(&self, closure: Option<usize>) -> Result<()> {
        if self.label.contains(['\t', '\n', '\r']) {
            return Err(FactorizedError::Graph {
                detail: "the label contains a TAB or a newline".into(),
            });
        }
        let edges = self.check_graph()?;
        self.check_momenta(closure.unwrap_or(self.loops() + self.externals()))?;
        // Factors: kinds, dimensions, indices.
        let mut analysed = Vec::<(Location, IndexSets, Elements)>::new();
        for edge in &self.edges {
            let rules = FactorRules {
                location: Location::Edge(edge.id),
                dimension: self.dimension,
                edges: &edges,
                boundary: None,
            };
            let (sets, elements) = rules.analyse(&edge.factor)?;
            analysed.push((rules.location, sets, elements));
        }
        for vertex in &self.vertices {
            let rules = FactorRules {
                location: Location::Vertex(vertex.id),
                dimension: self.dimension,
                edges: &edges,
                boundary: None,
            };
            let (sets, elements) = rules.analyse(&vertex.factor)?;
            analysed.push((rules.location, sets, elements));
        }
        analysed.extend(self.check_insertions(&edges)?);
        if let Some(projector) = &self.projector {
            let rules = FactorRules {
                location: Location::Projector,
                dimension: self.dimension,
                edges: &edges,
                boundary: None,
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
        if closure.is_some() {
            return Ok(());
        }
        let loops = self.closed_loops()?;
        // Cycles inside insertions are the original diagram's too.
        let (inner_fermion, inner_ghost) =
            self.insertions.iter().fold((0usize, 0usize), |(f, g), i| {
                (f + i.cycles.0 as usize, g + i.cycles.1 as usize)
            });
        let n_loops = loops.len() + inner_fermion + inner_ghost;
        let n_fermion =
            loops.iter().filter(|l| l.line == LineKind::Fermion).count() + inner_fermion;
        // a v3 direct record with `original_cycles`: the weights are the
        // caller's (R2): loop_sign 1, nf_power 0, overall 1, and the carried
        // closed cycles do not exceed the original graph's
        if let Some(DirectRecord {
            original_cycles: Some((fermion, ghost)),
            ..
        }) = &self.record
        {
            let n_ghost = n_loops - n_fermion;
            if p.loop_sign != 1 || p.nf_power != 0 || p.overall != Atom::num(1) {
                return Err(FactorizedError::NotScalar {
                    location: "prefactor".into(),
                    detail: "with original_cycles the weights are the caller's: loop_sign=1, \
                             nf_power=0, prefactor_overall 1"
                        .into(),
                });
            }
            if n_fermion > *fermion as usize || n_ghost > *ghost as usize {
                return Err(FactorizedError::Graph {
                    detail: format!(
                        "the record and its insertions carry {n_fermion} fermion / {n_ghost} ghost cycles, more than \
                         original_cycles fermion={fermion} ghost={ghost}"
                    ),
                });
            }
            return Ok(());
        }
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

    /// The insertion at vertex `v`, if any.
    pub fn insertion(&self, vertex: usize) -> Option<&Insertion> {
        self.insertions.iter().find(|i| i.vertex == vertex)
    }

    /// The edges owning one of `hedges` (at either end).
    fn edges_owning(&self, hedges: &[usize]) -> BTreeSet<usize> {
        self.edges
            .iter()
            .filter(|e| {
                e.source_hedge.is_some_and(|h| hedges.contains(&h))
                    || e.target_hedge.is_some_and(|h| hedges.contains(&h))
            })
            .map(|e| e.id)
            .collect()
    }

    /// The vertex at which half-edge `h` sits.
    fn hedge_vertex(&self, h: usize) -> Option<usize> {
        self.edges.iter().find_map(|e| match (e.source, e.target) {
            (Endpoint::Vertex(v), _) if e.source_hedge == Some(h) => Some(v),
            (_, Endpoint::Vertex(v)) if e.target_hedge == Some(h) => Some(v),
            _ => None,
        })
    }

    /// The structural rules of the insertions ([`Insertion`]) and the index
    /// analysis of their factors (for the locality, contraction and
    /// orientation checks of the whole value).
    pub(crate) fn check_insertions(
        &self,
        edges: &BTreeMap<usize, &FactorizedEdge>,
    ) -> Result<Vec<(Location, IndexSets, Elements)>> {
        let mut analysed = Vec::new();
        let present = self.hedges();
        for window in self.insertions.windows(2) {
            if window[0].vertex >= window[1].vertex {
                return Err(FactorizedError::Insertion {
                    vertex: window[1].vertex,
                    detail: "insertions are not sorted by vertex, or two share a vertex".into(),
                });
            }
        }
        for insertion in &self.insertions {
            let v = insertion.vertex;
            let error = |detail: String| FactorizedError::Insertion { vertex: v, detail };
            let Some(vertex) = self.vertex(v) else {
                return Err(error("the vertex does not exist".into()));
            };
            if vertex.factor != Atom::num(1) {
                return Err(error(format!(
                    "the vertex factor must be 1 (the insertion's factor is stated once), found {}",
                    vertex.factor.to_canonical_string()
                )));
            }
            if insertion.hedges.is_empty() || insertion.hedges.windows(2).any(|w| w[0] >= w[1]) {
                return Err(error(
                    "hedges must be non-empty, sorted and distinct".into(),
                ));
            }
            if let Some(h) = insertion.hedges.iter().find(|h| !present.contains(h)) {
                return Err(error(format!("half-edge {h} is not in the graph")));
            }
            let at_vertex = self.local_hedges(Location::Vertex(v));
            if !insertion.hedges.iter().any(|h| at_vertex.contains(h)) {
                return Err(error("no boundary half-edge sits at the vertex".into()));
            }
            if vertex.rule.starts_with("ct:") {
                if insertion.hedges.iter().copied().collect::<BTreeSet<_>>() != at_vertex {
                    return Err(error(format!(
                        "a contracted vertex (rule {}) has hedges {:?}, its incident half-edges are {:?}",
                        vertex.rule, insertion.hedges, at_vertex
                    )));
                }
            } else if vertex.rule.starts_with("ins:") {
                if insertion.hedges.len() != 2 {
                    return Err(error(format!(
                        "a surrogate coefficient vertex (rule {}) has {} hedges, not the two inner boundary half-edges",
                        vertex.rule,
                        insertion.hedges.len()
                    )));
                }
            } else {
                return Err(error(format!(
                    "the vertex rule {} is neither `ct:<label>` (a contracted child) nor `ins:<label>` (a surrogate's coefficient vertex)",
                    vertex.rule
                )));
            }
            for (what, atom) in [
                ("coefficient", &insertion.coefficient),
                ("factor", &insertion.factor),
            ] {
                // the format's variants (section 8): colour -> ColourInFactor,
                // a flavour symbol or mUVexp -> Graph
                if let Some(found) = forbidden_in_insertion(atom.as_view()) {
                    let location = format!("insertion at vertex {v} {what}");
                    return Err(if found.starts_with("the colour object") {
                        FactorizedError::ColourInFactor {
                            location,
                            tensor: found,
                        }
                    } else {
                        FactorizedError::Graph {
                            detail: format!("{location} contains {found}"),
                        }
                    });
                }
            }
            if let Some(name) = self
                .mass_symbol_in(&insertion.coefficient)
                .or_else(|| self.mass_symbol_in(&insertion.factor))
            {
                return Err(error(format!("contains a mass symbol {name}")));
            }
            // the pole order's rules (Laurent in ε, bounded), ε recognised by
            // its name here (the settings are not known)
            let dimension = self.dimension;
            let depends = |s: Symbol| is_epsilon_symbol(s) || s == dimension;
            insertion_epsilon_rules(insertion, &depends)?;
            check_scalar(
                insertion.coefficient.as_view(),
                &format!("insertion at vertex {v} coefficient"),
            )?;
            let rules = FactorRules {
                location: Location::Insertion(v),
                dimension: self.dimension,
                edges,
                boundary: Some(self.edges_owning(&insertion.hedges)),
            };
            let (sets, elements) = rules.analyse(&insertion.factor)?;
            // Locality: every free slot is one of the boundary half-edges'
            // (or an open port of a half-edge the graph does not have).
            for (key, slot) in &sets.free {
                let local = slot.hedge().is_some_and(|h| insertion.hedges.contains(&h));
                let absent_port = slot.hedge().is_some_and(|h| !present.contains(&h))
                    && self.ports.iter().any(|p| p.to_canonical_string() == *key);
                if !local && !absent_port {
                    return Err(FactorizedError::IndexNotLocal {
                        slot: key.clone(),
                        location: Location::Insertion(v).to_string(),
                    });
                }
            }
            if sets.free.len() > INSERTION_MAX_SLOTS {
                return Err(error(format!(
                    "the factor has {} open slots, at most {INSERTION_MAX_SLOTS} are supported",
                    sets.free.len()
                )));
            }
            let degree = momentum_degree(insertion.factor.as_view());
            if degree > INSERTION_MAX_DEGREE {
                return Err(error(format!(
                    "the factor has degree {degree} in the boundary momenta, at most {INSERTION_MAX_DEGREE} is supported"
                )));
            }
            analysed.push((rules.location, sets, elements));
        }
        Ok(analysed)
    }

    /// An insertion may depend on the evaluation's ε (`epsilon`, the full
    /// name of the settings' symbol) but on no other symbol named `ε` /
    /// `…::ε` (a second ε would be silently constant).
    pub(crate) fn check_insertion_epsilon(&self, epsilon: &str) -> Result<()> {
        for insertion in &self.insertions {
            let mut other = None;
            for atom in [&insertion.coefficient, &insertion.factor] {
                atom.visitor(&mut |part| {
                    if let AtomView::Var(v) = part {
                        let symbol = v.get_symbol();
                        let name = symbol.get_name();
                        if (name == "ε" || name.ends_with("::ε")) && name != epsilon {
                            other = Some(name.to_string());
                        }
                    }
                    other.is_none()
                });
            }
            if let Some(name) = other {
                return Err(FactorizedError::Insertion {
                    vertex: insertion.vertex,
                    detail: format!("the symbol {name} is not the evaluation's epsilon {epsilon}"),
                });
            }
        }
        // validation recognises ε by its name; an epsilon setting of another
        // name gets the pole order's rules here, where it is known
        if !(epsilon == "ε" || epsilon.ends_with("::ε")) {
            let dimension = self.dimension;
            let depends = |s: Symbol| s.get_name() == epsilon || s == dimension;
            for insertion in &self.insertions {
                insertion_epsilon_rules(insertion, &depends)?;
            }
        }
        Ok(())
    }

    /// The first mass symbol occurring in `atom`: a symbol of an edge's mass,
    /// or R*'s reference / auxiliary mass symbols by name.
    fn mass_symbol_in(&self, atom: &Atom) -> Option<String> {
        let mut masses = BTreeSet::new();
        for edge in &self.edges {
            edge.mass_squared.visitor(&mut |part| {
                if let AtomView::Var(v) = part {
                    masses.insert(v.get_symbol().get_name().to_string());
                }
                true
            });
        }
        let mut found = None;
        atom.visitor(&mut |part| {
            if let AtomView::Var(v) = part {
                let symbol = v.get_symbol();
                let name = symbol.get_name();
                if masses.contains(name)
                    || name.ends_with("::rstar_reference_mass_squared")
                    || name.ends_with("::rstar_mass_squared")
                {
                    found = Some(name.to_string());
                }
            }
            found.is_none()
        });
        found
    }

    /// The POLE ORDER of the insertions: the SUM over the insertions of
    /// `max(0, −v)`, `v` the trailing power of `epsilon` of `coefficient ×
    /// factor` (the dimension `D → 4 − 2 ε` first; tensors are constants in
    /// ε), i.e. an upper bound of the pole order of their PRODUCT (exact
    /// unless leading terms cancel or another insertion has a zero), which
    /// multiplies the rest of the record; 0 when no insertion has a pole. An entry that integrates the rest
    /// of the record separately and multiplies the insertions afterwards must
    /// raise the integral's ε depth by it (a pole promotes higher orders of
    /// the integral); one that hands the insertions to vakint inside the
    /// integrand (the direct and vacuum entries) needs nothing.
    pub fn insertion_pole_order(&self, epsilon: Symbol) -> Result<u64> {
        use symbolica::poly::series::SeriesDepth;
        let in_epsilon = Atom::num(4) - Atom::num(2) * Atom::var(epsilon);
        let mut order = 0u64;
        for insertion in &self.insertions {
            let product = (&insertion.coefficient * &insertion.factor)
                .replace(Atom::var(self.dimension))
                .with(in_epsilon.clone())
                .expand();
            if product.is_zero() {
                continue;
            }
            let failure = |detail: String| FactorizedError::Insertion {
                vertex: insertion.vertex,
                detail,
            };
            // a failing series expansion is internal (`Evaluation`), not a
            // refusal of the record
            let series = product
                .series(epsilon, Atom::Zero.as_view(), SeriesDepth::relative(1))
                .map_err(|e| FactorizedError::Evaluation {
                    detail: format!(
                        "insertion at vertex {}: the Laurent series in ε failed: {e}",
                        insertion.vertex
                    ),
                })?;
            if series.is_zero() {
                continue;
            }
            let exponent = series.get_trailing_exponent();
            if !exponent.is_integer() {
                return Err(failure("a fractional power of ε".into()));
            }
            let v = exponent
                .numerator()
                .to_i64()
                .ok_or_else(|| failure("an ε power out of range".into()))?;
            if v < 0 {
                order = order.checked_add(v.unsigned_abs()).ok_or_else(|| {
                    failure("an ε power out of range: the insertions' pole order overflows".into())
                })?;
            }
        }
        Ok(order)
    }

    /// The closure entry's depth raise: the insertions' exact pole order
    /// ([`Self::insertion_pole_order`]), refused above
    /// [`INSERTION_MAX_TOTAL_POLE_ORDER`] (`Insertion` at the first
    /// insertion's vertex, `an ε power out of range: the insertions' pole
    /// orders add up to <n>, above 64`).
    pub(crate) fn insertion_depth_raise(&self, epsilon: Symbol) -> Result<u64> {
        let order = self.insertion_pole_order(epsilon)?;
        if order > INSERTION_MAX_TOTAL_POLE_ORDER {
            return Err(FactorizedError::Insertion {
                vertex: self.insertions.first().map_or(0, |i| i.vertex),
                detail: format!(
                    "an ε power out of range: the insertions' pole orders add up to {order}, above {INSERTION_MAX_TOTAL_POLE_ORDER}"
                ),
            });
        }
        Ok(order)
    }

    /// The node of every vertex for the closed-loop flow: the vertices that
    /// own the boundary half-edges of one insertion (and its vertex) are one
    /// node; every other vertex is its own.
    fn flow_nodes(&self) -> BTreeMap<usize, usize> {
        let mut node = self
            .vertices
            .iter()
            .map(|v| (v.id, v.id))
            .collect::<BTreeMap<_, _>>();
        fn find(node: &BTreeMap<usize, usize>, mut v: usize) -> usize {
            while let Some(&parent) = node.get(&v) {
                if parent == v {
                    break;
                }
                v = parent;
            }
            v
        }
        for insertion in &self.insertions {
            let members = std::iter::once(insertion.vertex)
                .chain(
                    insertion
                        .hedges
                        .iter()
                        .filter_map(|h| self.hedge_vertex(*h)),
                )
                .collect::<Vec<_>>();
            for member in members {
                let (a, b) = (find(&node, insertion.vertex), find(&node, member));
                if a != b {
                    let (low, high) = (a.min(b), a.max(b));
                    node.insert(high, low);
                }
            }
        }
        let ids = node.keys().copied().collect::<Vec<_>>();
        ids.into_iter().map(|v| (v, find(&node, v))).collect()
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

    /// Momentum vectors, basis edges and conservation of the first
    /// `conserved` components at every vertex.
    fn check_momenta(&self, conserved: usize) -> Result<()> {
        let n = self.loops() + self.externals();
        for edge in &self.edges {
            if edge.momentum.len() != n {
                return Err(FactorizedError::MomentumLength {
                    edge: edge.id,
                    expected: n,
                    found: edge.momentum.len(),
                });
            }
            // coefficients are multiplied by ±1 and summed per vertex: bounded
            // so that no product or sum of fewer than 2^32 edges can wrap
            if let Some(c) = edge.momentum.iter().find(|c| c.unsigned_abs() > 1 << 31) {
                return Err(FactorizedError::MomentumBasis {
                    detail: format!(
                        "edge {}: momentum coefficient {c} out of range (at most 2^31 in absolute value)",
                        edge.id
                    ),
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
            if residual.iter().take(conserved).any(|c| *c != 0) {
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
            .chain(self.insertions.iter().map(|i| &i.factor))
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
            Location::Insertion(v) => self
                .insertion(v)
                .map(|i| i.hedges.iter().copied().collect())
                .unwrap_or_default(),
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
        // The projector closes some or all ports; a free slot of the projector
        // that is not a port must not occur in the factors: it stays OPEN
        // (`open_slots`). Ports it does not close stay open too.
        if let Some(projector_free) = &projector_free {
            for key in projector_free.difference(&ports) {
                if free
                    .get(key)
                    .is_some_and(|locations| locations.iter().any(|l| *l != Location::Projector))
                {
                    return Err(FactorizedError::Ports {
                        detail: format!(
                            "the projector's free index {key} occurs in the factors but is not a \
                             port (the projector's free indices {projector_free:?}, the ports \
                             {ports:?})"
                        ),
                    });
                }
            }
        }
        for (key, locations) in &free {
            if locations.iter().all(|l| *l == Location::Projector) && !ports.contains(key) {
                // An open slot of the projector that is not a port (checked
                // above); a port must still occur once in the factors.
                continue;
            }
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
                        // An amputated external fermion leg (factor without
                        // spinors): its slot at the vertex is a port.
                        let amputated_leg = edge.is_external() && spinor_free.is_empty();
                        if spinor_free != own && !amputated_leg {
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
                Location::Insertion(v) => {
                    // The boundary half-edges of the insertion: rows on the
                    // outgoing, columns on the incoming fermion lines.
                    let hedges = self.local_hedges(Location::Insertion(*v));
                    let mut rows = BTreeSet::new();
                    let mut columns = BTreeSet::new();
                    for e in self.edges.iter().filter(|e| e.line == LineKind::Fermion) {
                        rows.extend(e.source_hedge.filter(|h| hedges.contains(h)));
                        columns.extend(e.target_hedge.filter(|h| hedges.contains(h)));
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
        if !self.has_new_features() {
            return self.closed_loops_v1();
        }
        let mut loops = Vec::new();
        let open = if self.projector.is_some() {
            self.open_slots()?
        } else {
            self.ports.clone()
        };
        let nodes = self.flow_nodes();
        for line in [LineKind::Fermion, LineKind::Ghost] {
            let lines = self
                .edges
                .iter()
                .filter(|e| e.line == line)
                .collect::<Vec<_>>();
            if !line.is_directed_loop() {
                continue;
            }
            // An open FERMION line whose two legs' spinor slots are ports the
            // projector closes (the caller's Dirac closure, e.g. the quark
            // form factor Tr[slash(p) X]) is traced by the contraction: it is
            // neither an open line nor a closed loop (no loop sign, no n_f).
            let closed_by_projector = |edge: &FactorizedEdge| {
                line == LineKind::Fermion
                    && self.projector.is_some()
                    && edge.source_hedge.or(edge.target_hedge).is_some_and(|h| {
                        self.ports.iter().any(|port| {
                            as_slot(port.as_view()).is_some_and(|slot| {
                                slot.rep == names::BIS
                                    && slot.hedge() == Some(h)
                                    && !open.contains(port)
                            })
                        })
                    })
            };
            let legs = lines.iter().filter(|e| e.is_external()).collect::<Vec<_>>();
            let incoming = legs
                .iter()
                .filter(|e| e.source == Endpoint::External)
                .count();
            for edge in &legs {
                if !closed_by_projector(edge) || legs.len() != 2 || incoming != 1 {
                    return Err(FactorizedError::OpenLine {
                        edge: edge.id,
                        line: line.as_str().into(),
                    });
                }
            }
            let node_of = |end: Endpoint| match end {
                Endpoint::Vertex(v) => Endpoint::Vertex(nodes.get(&v).copied().unwrap_or(v)),
                Endpoint::External => Endpoint::External,
            };
            let node_ids = nodes.values().copied().collect::<BTreeSet<_>>();
            for node in node_ids {
                let incoming = lines
                    .iter()
                    .filter(|e| node_of(e.target) == Endpoint::Vertex(node))
                    .count();
                let outgoing = lines
                    .iter()
                    .filter(|e| node_of(e.source) == Endpoint::Vertex(node))
                    .count();
                if incoming != outgoing || incoming > 1 {
                    return Err(FactorizedError::FlowNotConserved {
                        vertex: node,
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
                        .find(|e| node_of(e.source) == node_of(current.target))
                        .expect("conserved flow");
                    if next.id == start {
                        break;
                    }
                    unused.remove(&next.id);
                    cycle.push(next.id);
                    current = next;
                }
                // the projector-closed line runs through the external legs
                if cycle
                    .iter()
                    .any(|id| lines.iter().any(|e| e.id == *id && e.is_external()))
                {
                    continue;
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
        let value = Self::parse_text(input)?;
        value.validate()?;
        Ok(value)
    }

    /// Read one record of the text form without validating it (for a value
    /// with new features, compact slashes in the projector are written out on
    /// read: [`Self::expand_compact_slashes`]).
    pub fn parse_text(input: &str) -> Result<Self> {
        let mut value = text::parse_record(input)?;
        if value.has_new_features() {
            value.expand_compact_slashes()?;
        }
        Ok(value)
    }

    /// Rule R4 of the record format for the projector: every COMPACT slash
    /// `gamma(a, b, Q(x, mink(D)))` is written out by [`expand_compact_slash`]
    /// (scope -1). Refused under a Taylor step (a compact slash of an external
    /// momentum in the projector of a closure record with an unapplied Taylor
    /// specification). Idempotent; a value without compact slash is unchanged.
    pub fn expand_compact_slashes(&mut self) -> Result<()> {
        let Some(projector) = &self.projector else {
            return Ok(());
        };
        let expanded = expand_compact_slash(&projector.tensor, -1);
        if expanded == projector.tensor {
            return Ok(());
        }
        if let FactorizedKind::Closure(spec) = &self.kind
            && spec.taylor.is_some()
        {
            return Err(FactorizedError::Ports {
                detail: "a compact slash gamma(a, b, Q(x, mink(D))) in the projector of a record \
                         with a Taylor step (write it out, or close after the Taylor step)"
                    .into(),
            });
        }
        if let Some(projector) = &mut self.projector {
            projector.tensor = expanded;
        }
        Ok(())
    }
}

/// The base's validation, for values without new features.
impl FactorizedNumerator {
    /// Whether the value uses a feature the base did not have: the closure
    /// kind, insertion vertices or a schema-v3 direct record. A value without
    /// any is validated and evaluated by exactly the base's rules (every
    /// pre-existing entry branches on this once); the widened grammar of the
    /// new record kinds applies to marked values only.
    pub fn has_new_features(&self) -> bool {
        matches!(self.kind, FactorizedKind::Closure(_))
            || !self.insertions.is_empty()
            || self.record.is_some()
    }

    /// The base's `validate`, verbatim (the rules of a value without new
    /// features: [`FactorizedNumerator::has_new_features`]).
    pub(crate) fn validate_v1(&self) -> Result<()> {
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
        self.check_momenta_v1()?;
        // Factors: kinds, dimensions, indices.
        let mut analysed = Vec::<(Location, IndexSets, Elements)>::new();
        for edge in &self.edges {
            let rules = FactorRulesV1 {
                location: Location::Edge(edge.id),
                dimension: self.dimension,
                edges: &edges,
            };
            let (sets, elements) = rules.analyse(&edge.factor)?;
            analysed.push((rules.location, sets, elements));
        }
        for vertex in &self.vertices {
            let rules = FactorRulesV1 {
                location: Location::Vertex(vertex.id),
                dimension: self.dimension,
                edges: &edges,
            };
            let (sets, elements) = rules.analyse(&vertex.factor)?;
            analysed.push((rules.location, sets, elements));
        }
        if let Some(projector) = &self.projector {
            let rules = FactorRulesV1 {
                location: Location::Projector,
                dimension: self.dimension,
                edges: &edges,
            };
            let (sets, elements) = rules.analyse(&projector.tensor)?;
            analysed.push((rules.location, sets, elements));
        }
        self.check_locality(&analysed)?;
        self.check_contraction_v1(&analysed)?;
        self.check_fermion_orientation_v1(&analysed)?;
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
        let loops = self.closed_loops_v1()?;
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

    /// The base's `check_momenta`, verbatim (the rules of a value without new
    /// features: [`FactorizedNumerator::has_new_features`]).
    fn check_momenta_v1(&self) -> Result<()> {
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

    /// The base's `check_contraction`, verbatim (the rules of a value without new
    /// features: [`FactorizedNumerator::has_new_features`]).
    fn check_contraction_v1(&self, analysed: &[(Location, IndexSets, Elements)]) -> Result<()> {
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

    /// The base's `check_fermion_orientation`, verbatim (the rules of a value without new
    /// features: [`FactorizedNumerator::has_new_features`]).
    fn check_fermion_orientation_v1(
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
                // no insertion in a value without new features
                Location::Insertion(_) | Location::Projector => continue,
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

    /// The base's `closed_loops`, verbatim (the rules of a value without new
    /// features: [`FactorizedNumerator::has_new_features`]).
    pub(crate) fn closed_loops_v1(&self) -> Result<Vec<ClosedLoop>> {
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
}

#[cfg(test)]
mod insertion_records_tests;
#[cfg(test)]
mod insertion_tests;
#[cfg(test)]
mod tests;
