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

mod backstop;
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

/// The schemas a TextSchema refusal names: the base's text (only v1) for any
/// line but the v3 name with surrounding whitespace (e.g. a trailing space);
/// for that line both, and how they must be written.
fn known_schemas(found: &str) -> String {
    if found != FACTORIZED_SCHEMA_V3 && found.trim() == FACTORIZED_SCHEMA_V3 {
        format!(
            "{FACTORIZED_SCHEMA} and {FACTORIZED_SCHEMA_V3}, written exactly, without trailing whitespace"
        )
    } else {
        FACTORIZED_SCHEMA.to_string()
    }
}

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
    /// the record format (kinds K1-K3). The adapter and the vacuum
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
/// call (the record format, the weights rule R2); without, the v1 rule.
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
    #[error(
        "text line {line}: unknown schema {found:?} (this reader knows {})",
        known_schemas(found)
    )]
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
/// internal code only as the evaluation's epsilon. The one ATOM of vakint's
/// namespace it may carry besides: an expected but unknown master
/// coefficient exactly as vakint's results write it,
/// `UNKNOWN(master, k)` ([`crate::truncation_guard::unknown_coefficient_symbol`]
/// with a master symbol and an integer: [`is_unknown_coefficient`]), a plain
/// scalar carried linearly.
pub const ALLOWED_VAKINT_SYMBOL: &str = "vakint::ε";

/// The largest |order| of an accepted `UNKNOWN(master, k)`.
pub const UNKNOWN_MAX_ORDER: i64 = 64;

/// THE VAKINT DATA ATOMS a marked record may carry (the one exception to the
/// reserved `vakint::` namespace): the constants vakint's own p-integral
/// values and master catalogues write, each with its fixed argument shape,
/// accepted only in an insertion's coefficient and factor, never inside a
/// function argument or an exponent ([`is_vakint_data_atom`]):
/// - `vakint::UNKNOWN(vakint::<master>, k)`: an expected but unknown master
///   coefficient (written by `truncation_guard::unknown_coefficient`, from a
///   version-2 catalogue's `UNKNOWN(symbol, k)` entries); `<master>` a
///   catalogue master name ([`is_master_name`]; at an entry, the loaded
///   catalogue's), `k` an integer of modulus at most [`UNKNOWN_MAX_ORDER`];
///   carried linearly;
/// - `vakint::MZV(n_1, .., n_r)`: a multiple zeta value (a catalogue
///   coefficient's constant, kept by `master_catalogue::parse_constant`,
///   evaluated by `mzv`), 1 to 16 positive integers of at most 64, the first
///   at least 2: exactly the convergent ones vakint's evaluator accepts
///   (`mzv::evaluable_indices`; MZV(1, ..) diverges).
///
/// ζ(n) and π are Symbolica's own symbols (`symbolica::zeta`, `𝜋`), not
/// vakint's.
pub const VAKINT_DATA_ATOMS: [(&str, &str); 2] = [
    (
        "vakint::UNKNOWN",
        "(vakint::<catalogue master>, <integer k, |k| <= 64>)",
    ),
    (
        "vakint::MZV",
        "(<1 to 16 positive integers <= 64, the first >= 2>)",
    ),
];

/// Whether `name` is a master name of vakint's catalogues by their naming
/// pattern (`vakint::ggppss_int_<n>`, `vakint::x_t<n>(_<n>)*`, the names of
/// the built-in and the shipped catalogues); at an entry with a catalogue
/// loaded the name must also be one of its masters.
fn is_master_name(name: &str) -> bool {
    let Some(rest) = name.strip_prefix("vakint::") else {
        return false;
    };
    let digits = |t: &str| !t.is_empty() && t.bytes().all(|b| b.is_ascii_digit());
    if let Some(n) = rest.strip_prefix("ggppss_int_") {
        return digits(n);
    }
    if let Some(groups) = rest.strip_prefix("x_t") {
        return groups.split('_').all(digits);
    }
    false
}

/// Whether `view` is vakint's `MZV(n_1, .., n_r)` in its data shape.
fn is_mzv(view: AtomView<'_>) -> bool {
    let AtomView::Fun(f) = view else {
        return false;
    };
    if f.get_symbol().get_name() != "vakint::MZV" || !(1..=16).contains(&f.get_nargs()) {
        return false;
    }
    let Some(indices) = f
        .iter()
        .map(|a| i64::try_from(a).ok().filter(|n| (1..=64).contains(n)))
        .collect::<Option<Vec<i64>>>()
    else {
        return false;
    };
    // the evaluator's own rule (a leading 1 diverges)
    crate::p_integral::mzv::evaluable_indices(&indices)
}

/// Whether `view` is one of [`VAKINT_DATA_ATOMS`] in its exact shape.
pub(crate) fn is_vakint_data_atom(view: AtomView<'_>) -> bool {
    is_unknown_coefficient(view) || is_mzv(view)
}

/// The first data atom of `view` inside a function argument or an exponent
/// (where no data atom may stand), shown.
fn misplaced_data_atom(view: AtomView<'_>, inside: bool) -> Option<String> {
    if is_vakint_data_atom(view) {
        return inside.then(|| view.to_owned().to_canonical_string());
    }
    match view {
        AtomView::Fun(f) => f.iter().find_map(|a| misplaced_data_atom(a, true)),
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            misplaced_data_atom(base, inside).or_else(|| misplaced_data_atom(exponent, true))
        }
        AtomView::Add(a) => a.iter().find_map(|t| misplaced_data_atom(t, inside)),
        AtomView::Mul(m) => m.iter().find_map(|t| misplaced_data_atom(t, inside)),
        _ => None,
    }
}

/// Whether `view` is vakint's `UNKNOWN(master, k)` atom in exactly the shape
/// its results write it: the head of
/// [`crate::truncation_guard::unknown_coefficient_symbol`], a master name
/// ([`is_master_name`]) and an integer order of at most
/// [`UNKNOWN_MAX_ORDER`] in modulus.
pub(crate) fn is_unknown_coefficient(view: AtomView<'_>) -> bool {
    let AtomView::Fun(f) = view else {
        return false;
    };
    let args = f.iter().collect::<Vec<_>>();
    f.get_symbol() == crate::truncation_guard::unknown_coefficient_symbol()
        && args.len() == 2
        && matches!(args[0], AtomView::Var(v) if is_master_name(v.get_symbol().get_name()))
        && i64::try_from(args[1])
            .is_ok_and(|k| k.unsigned_abs() <= UNKNOWN_MAX_ORDER.unsigned_abs())
}

/// The degree of `view` in vakint's DATA ATOMS ([`is_vakint_data_atom`]:
/// `UNKNOWN(master, k)` and `MZV(..)`; a sum takes the maximum, a product the
/// sum); `None` when a data atom stands inside the base of a power (whatever
/// the exponent: a power, a denominator). A data atom is carried as a LINEAR
/// FACTOR only: degree at most 1 and never under a power (every monomial of
/// both catalogues' expanded coefficients has degree 0 or 1).
fn data_atom_degree(view: AtomView<'_>) -> Option<u64> {
    match view {
        AtomView::Fun(_) if is_vakint_data_atom(view) => Some(1),
        AtomView::Add(a) => a
            .iter()
            .map(data_atom_degree)
            .try_fold(0, |m, d| Some(m.max(d?))),
        AtomView::Mul(m) => m
            .iter()
            .map(data_atom_degree)
            .try_fold(0u64, |s, d| s.checked_add(d?)),
        AtomView::Pow(p) => match data_atom_degree(p.get_base_exp().0)? {
            0 => Some(0),
            _ => None,
        },
        _ => Some(0),
    }
}

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
    first_forbidden_symbol_with(view, exempt, false)
}

/// [`first_forbidden_symbol`]; `data_atoms`: [`VAKINT_DATA_ATOMS`] are
/// accepted (an insertion's coefficient and factor).
pub(crate) fn first_forbidden_symbol_with(
    view: AtomView<'_>,
    exempt: &[&str],
    data_atoms: bool,
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
                // vakint's data atoms, where the field admits them
                if data_atoms && is_vakint_data_atom(part) {
                    return false;
                }
                if name == names::MINK && f.get_nargs() == 2 {
                    let args = f.iter().collect::<Vec<_>>();
                    if let AtomView::Fun(index) = args[1]
                        && index.get_symbol().get_name() == "vakint::slash"
                        && index.iter().all(|a| matches!(a, AtomView::Num(_)))
                    {
                        found = first_forbidden_symbol_with(args[0], exempt, data_atoms);
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

/// [`check_caller_atom`] for an insertion's coefficient or factor, where
/// [`VAKINT_DATA_ATOMS`] are accepted.
pub(crate) fn check_caller_atom_with_data_atoms(
    atom: &Atom,
    field: &str,
    exempt: &[&str],
    source: AtomSource,
) -> Result<()> {
    match first_forbidden_symbol_with(atom.as_view(), exempt, true) {
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
    value.check_opaque_dependence(prefactor, FIELD)?;
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
    if let FactorizedKind::Closure(spec) = &value.kind
        && spec.auxiliary_mass_squared != spec.mass_squared
        && closure::contains_atom(prefactor.as_view(), spec.auxiliary_mass_squared.as_view())
    {
        return Err(refuse(format!(
            "the auxiliary mass {} (it may occur only as the taylor line's mass)",
            spec.auxiliary_mass_squared
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
///
/// INVARIANT: no rule that COMPUTES (a zero test, `together`, `expand`, a
/// series) runs on an atom that has not passed the SHAPE rules (names, exact
/// numbers, data atoms, exponents of at most 64, heads, structural shapes):
/// every entry calls [`check_master_modes`] (the settings) before this
/// function and [`check_entry_names`] (the non-ASCII names: a walk, no
/// computation) after it; here the record's atoms get
/// the reserved-name walk and [`FactorizedNumerator::check_marked_atoms`]
/// first (the entries may run this before `validate`), the prefactor argument
/// its own in [`check_prefactor_argument`], before the defined-and-finite
/// checks.
pub(crate) fn check_marked_entry(
    value: &FactorizedNumerator,
    prefactor: &Atom,
    epsilon: Symbol,
    terms: Option<i64>,
) -> Result<()> {
    if is_marked(value) {
        // the order of `validate`: the reserved names, then the atom rules
        value.check_reserved_symbols()?;
        value.check_marked_atoms()?;
    }
    check_prefactor_argument(value, prefactor, epsilon)?;
    check_epsilon_setting(value, epsilon)?;
    if is_marked(value) {
        // the entry's prefactor argument: defined and finite at ε = 0
        check_finite_expression(
            prefactor,
            "prefactor argument",
            value.dimension,
            epsilon,
            true,
        )?;
        // the record fields with the settings' ε (validation took the
        // insertions' ε, or vakint's own without insertions)
        value.check_marked_fields(Some(epsilon))?;
    }
    // the settings' number of terms (an entry that does not take it from the
    // record), by the depth rule of [`MARKED_MAX_TERMS`]; marked values only
    // (the vacuum entries call this for every value)
    if let Some(terms) = terms.filter(|_| is_marked(value)) {
        marked_depth(i128::from(terms)).map_err(|detail| FactorizedError::Graph { detail })?;
        // a v3 direct record's depth sets the number of terms (format
        // section 7); these entries take the settings' number: they must agree
        if let Some(depth) = value.record.as_ref().and_then(|r| r.depth)
            && i128::try_from(depth).ok() != Some(i128::from(terms))
        {
            return Err(FactorizedError::Graph {
                detail: format!(
                    "the depth {depth} of the record differs from the settings' number of terms {terms} (they must be equal)"
                ),
            });
        }
    }
    Ok(())
}

/// THE MASTER MODES of the settings' evaluation order, for a marked value:
/// every method that has one (RustRed, Forcer) must have the SAME mode. The
/// entries decide from the first such method whether masters are retained
/// (no truncation of the closure product, Retain's output vocabulary), while
/// the method that actually runs depends on the record (a two-loop closure
/// ran RustRed in Symbolic under [Forcer(Retain), RustRed(Symbolic)] and came
/// back untruncated): a mixed order is refused instead of guessing (`Graph`,
/// prefix `the evaluation methods mix master modes: `). Methods without a
/// master mode (AlphaLoop, MATAD, FMFT, pySecDec) are not compared.
pub(crate) fn check_master_modes(settings: &crate::VakintSettings) -> Result<()> {
    let modes = settings
        .evaluation_order
        .0
        .iter()
        .filter_map(|method| match method {
            crate::EvaluationMethod::RustRed(options) => Some(options.master_mode),
            crate::EvaluationMethod::Forcer(options) => Some(options.master_mode),
            _ => None,
        })
        .collect::<Vec<_>>();
    if modes.windows(2).any(|w| w[0] != w[1]) {
        return Err(FactorizedError::Graph {
            detail: format!(
                "the evaluation methods mix master modes: {} (a marked value needs one mode for every RustRed and Forcer method)",
                settings
                    .evaluation_order
                    .0
                    .iter()
                    .map(|m| format!("{m:?}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        });
    }
    Ok(())
}

/// IDENTIFIERS of a marked value at an entry (where the settings are known):
/// a symbol or function head whose name has a character outside ASCII is
/// refused (`Graph`, detail prefix `a non-ASCII name `) unless it is the
/// settings' epsilon, the settings' scale symbol or one of Symbolica's own
/// symbols (`symbolica::…`: 𝜋, 𝑒); the imaginary unit 𝑖 is a number. A
/// lookalike of the epsilon (ϵ, 𝜀, 𝜖, or ε in another namespace) is thereby
/// refused instead of being read as a constant. A no-op for an unmarked value.
pub(crate) fn check_entry_names(
    value: &FactorizedNumerator,
    prefactor: &Atom,
    settings: &crate::VakintSettings,
) -> Result<()> {
    if !is_marked(value) {
        return Ok(());
    }
    let epsilon = crate::utils::vakint_macros::vk_symbol!(&settings.epsilon_symbol);
    let scale = crate::utils::vakint_macros::vk_symbol!(&settings.mu_r_sq_symbol);
    // Symbolica's non-ASCII built-ins by name (correction 3): the constants
    // 𝜋 and 𝑒; its built-in function heads are ASCII
    let allowed = |name: &str| {
        name.is_ascii()
            || name == epsilon.get_name()
            || name == scale.get_name()
            || name == "symbolica::𝜋"
            || name == "symbolica::𝑒"
    };
    let dimension = Atom::var(value.dimension);
    let mut atoms = value.record_atoms();
    atoms.push(("dimension".to_string(), &dimension));
    atoms.push(("prefactor argument".to_string(), prefactor));
    for (field, atom) in atoms {
        // an insertion: an epsilon spelled ε of another name is the
        // insertion rule's (section 8, `Insertion`); a lookalike is refused
        // here with the insertion's variant
        let insertion = field
            .strip_prefix("insertion ")
            .and_then(|rest| rest.split(' ').next())
            .and_then(|v| v.parse::<usize>().ok());
        let mut found = None;
        atom.as_view().visitor(&mut |part| {
            let name = match part {
                AtomView::Var(v) => v.get_symbol().get_name().to_string(),
                AtomView::Fun(f) => f.get_symbol().get_name().to_string(),
                _ => return true,
            };
            let insertion_epsilon = insertion.is_some() && (name == "ε" || name.ends_with("::ε"));
            if !allowed(&name) && !insertion_epsilon {
                found = Some(name);
            }
            found.is_none()
        });
        if let Some(name) = found {
            let detail = format!(
                "a non-ASCII name {name} in {field}: outside ASCII only the settings' epsilon {} and scale {} and Symbolica's built-ins (𝜋, 𝑒) are accepted",
                epsilon.get_name(),
                scale.get_name()
            );
            return Err(match insertion {
                Some(vertex) => FactorizedError::Insertion { vertex, detail },
                None => FactorizedError::Graph { detail },
            });
        }
    }
    Ok(())
}

/// The UNKNOWN atoms of a marked value's insertions name masters of the
/// engine's LOADED p master catalogue, when one is loaded (validation
/// checks the catalogues' naming pattern only): `Insertion`, detail prefix
/// `the UNKNOWN master `. A no-op for an unmarked value.
pub(crate) fn check_data_atoms_catalogue(
    value: &FactorizedNumerator,
    vakint: &crate::Vakint,
) -> Result<()> {
    if !is_marked(value) {
        return Ok(());
    }
    let Some(catalogue) = vakint.p_master_catalogue() else {
        return Ok(());
    };
    for insertion in &value.insertions {
        let mut missing = None;
        for atom in [&insertion.coefficient, &insertion.factor] {
            atom.as_view().visitor(&mut |part| {
                if is_unknown_coefficient(part)
                    && let AtomView::Fun(f) = part
                    && let Some(AtomView::Var(master)) = f.iter().next()
                    && catalogue
                        .id_of_symbol(master.get_symbol().get_stripped_name())
                        .is_none()
                {
                    missing.get_or_insert_with(|| master.get_symbol().get_name().to_string());
                }
                missing.is_none()
            });
        }
        if let Some(master) = missing {
            return Err(FactorizedError::Insertion {
                vertex: insertion.vertex,
                detail: format!(
                    "the UNKNOWN master {master} is not a master of the loaded catalogue {}",
                    catalogue.data_version()
                ),
            });
        }
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

/// Whether `base` is IDENTICALLY ZERO as a rational function of its symbols
/// (D, ε, a caller's scalars, function atoms alike) after expansion. First a
/// cheap exact pre-filter: fixed distinct integers for its variables,
/// symbolica's own normalisation (`replace` + exact rational arithmetic). Only
/// a finite NONZERO NUMBER is conclusive (a rational function that is nonzero
/// at one point is not identically zero, so one assignment suffices); every
/// other outcome (zero, a non-finite number because the point is a pole of an
/// inner power, a non-number such as a remaining function atom) falls through
/// to the exact normal forms: `together` (a common denominator) and `expand`. A
/// zero that holds only by an identity of special functions
/// (sin(x)^2 + cos(x)^2 - 1) is not seen: such a base is not zero here.
fn identically_zero(base: AtomView<'_>) -> bool {
    let mut symbols = Vec::new();
    base.visitor(&mut |part| {
        if let AtomView::Var(v) = part
            && !symbols.contains(&v.get_symbol())
        {
            symbols.push(v.get_symbol());
        }
        true
    });
    let mut point = base.to_owned();
    for (i, symbol) in symbols.iter().enumerate() {
        let value = 10007 + 101 * i64::try_from(i).unwrap_or(i64::MAX / 202);
        point = point.replace(Atom::var(*symbol)).with(Atom::num(value));
    }
    if let AtomView::Num(n) = point.as_view()
        && matches!(
            n.get_coeff_view(),
            CoefficientView::Natural(..) | CoefficientView::Large(..)
        )
        && !point.is_zero()
    {
        return false;
    }
    base.to_owned().together().is_zero() || base.to_owned().expand().is_zero()
}

/// Whether `exponent` is a positive rational number (a literal): the only
/// exponents under which a vanishing base is harmless (0^(1/2) = 0^2 = 0).
fn positive_rational(exponent: AtomView<'_>) -> bool {
    symbolica::domains::rational::Rational::try_from(exponent)
        .is_ok_and(|r| !r.is_negative() && !r.is_zero())
}

/// The first power in `view` whose exponent is NOT a positive rational
/// number (a negative rational, or no rational at all: sqrt(2), π, a symbol,
/// an expression such as P(y) - 1) and whose base is identically zero
/// ([`identically_zero`]) as written, or once the dimension is set to 4 - 2ε
/// (`epsilon`): the expression is undefined there, or of an undecidable kind
/// (0^π is 0, refused too: harmless, since no record needs it), and not a
/// pole of any order (a base that vanishes only at ε = 0, such as D - 4, is
/// an ordinary pole; one that vanishes only for some value of a caller's
/// scalar, x - 3, is defined). Every base is tested, whatever symbols it
/// contains; inner powers first, in bases, exponents and function arguments.
fn vanishing_denominator(
    view: AtomView<'_>,
    dimension: Symbol,
    epsilon: Symbol,
) -> Option<(bool, String)> {
    let inner = |v: AtomView<'_>| vanishing_denominator(v, dimension, epsilon);
    match view {
        AtomView::Add(a) => a.iter().find_map(inner),
        AtomView::Mul(m) => m.iter().find_map(inner),
        AtomView::Fun(f) => f.iter().find_map(inner),
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            if let Some(found) = inner(base).or_else(|| inner(exponent)) {
                return Some(found);
            }
            if positive_rational(exponent) {
                return None;
            }
            let shown = view.to_owned().to_canonical_string();
            if identically_zero(base) {
                return Some((
                    false,
                    format!("{shown} (its base is identically zero as written)"),
                ));
            }
            let mut dimension_inside = false;
            base.visitor(&mut |part| {
                if let AtomView::Var(v) = part
                    && v.get_symbol() == dimension
                {
                    dimension_inside = true;
                }
                !dimension_inside
            });
            if dimension_inside {
                let substituted = base
                    .to_owned()
                    .replace(Atom::var(dimension))
                    .with(Atom::num(4) - Atom::num(2) * Atom::var(epsilon));
                if identically_zero(substituted.as_view()) {
                    return Some((
                        true,
                        format!("{shown} (its base is identically zero after D -> 4-2ε)"),
                    ));
                }
            }
            None
        }
        _ => None,
    }
}

/// The first power with a negative exponent in `view` whose base, with the
/// dimension set to 4 - 2ε, vanishes at ε = 0: a CANDIDATE pole (1/(D-4),
/// ε^-1, 1/(ε + ε^2)); whether the expression has one is decided on the whole
/// expression ([`check_finite_expression`]). Inner powers first.
fn pole_denominator(view: AtomView<'_>, dimension: Symbol, epsilon: Symbol) -> Option<String> {
    let inner = |v: AtomView<'_>| pole_denominator(v, dimension, epsilon);
    match view {
        AtomView::Add(a) => a.iter().find_map(inner),
        AtomView::Mul(m) => m.iter().find_map(inner),
        AtomView::Fun(f) => f.iter().find_map(inner),
        AtomView::Pow(p) => {
            let (base, exponent) = p.get_base_exp();
            if let Some(found) = inner(base).or_else(|| inner(exponent)) {
                return Some(found);
            }
            let negative = symbolica::domains::rational::Rational::try_from(exponent)
                .is_ok_and(|r| r.is_negative());
            if !negative {
                return None;
            }
            let at_zero = base
                .to_owned()
                .replace(Atom::var(dimension))
                .with(Atom::num(4) - Atom::num(2) * Atom::var(epsilon))
                .replace(Atom::var(epsilon))
                .with(Atom::Zero);
            identically_zero(at_zero.as_view()).then(|| view.to_owned().to_canonical_string())
        }
        _ => None,
    }
}

/// The first symbol named ε (in any namespace) in `view`.
fn epsilon_symbol_in(view: AtomView<'_>) -> Option<Symbol> {
    let mut found = None;
    view.visitor(&mut |part| {
        if found.is_none()
            && let AtomView::Var(v) = part
            && is_epsilon_symbol(v.get_symbol())
        {
            found = Some(v.get_symbol());
        }
        found.is_none()
    });
    found
}

/// The pole order's rules for one insertion, `depends` selecting ε and the
/// dimension: the dependence on them is a Laurent polynomial
/// ([`non_laurent_dependence`]) and its bound is at most
/// [`INSERTION_MAX_EPSILON_POWER`] (`an ε power out of range`).
fn insertion_epsilon_rules(
    insertion: &Insertion,
    depends: &dyn Fn(Symbol) -> bool,
    dimension: Symbol,
    epsilon: Option<Symbol>,
) -> Result<()> {
    let v = insertion.vertex;
    let error = |detail: String| FactorizedError::Insertion { vertex: v, detail };
    // The pole order (`insertion_pole_order`) is that of coefficient x
    // factor as a Laurent expression in ε, the dimension set to
    // 4 - 2ε: ε and the dimension may occur only through sums,
    // products and integer powers (a representation slot's dimension
    // argument, mink(D, ..), is not a dependence). Inside a function's
    // argument, in an exponent or under a non-integer power the
    // expansion has no pole order (log ε would silently count as 0).
    // the exponent bound of the atom rules FIRST: the rules below compute
    // (an integer point raised to each power, exactly), so a power beyond
    // the bound is refused before any of them sees it
    for (what, atom) in [
        ("coefficient", &insertion.coefficient),
        ("factor", &insertion.factor),
    ] {
        if let Some(found) = first_exponent_beyond_bound(atom.as_view()) {
            return Err(error(format!(
                "a non-integer exponent in insertion {v} {what}: {found}"
            )));
        }
        if let Some(found) = first_oversized_number(atom.as_view()) {
            return Err(error(format!(
                "a number beyond {MAX_NUMBER_BITS} bits in insertion {v} {what}: {found}"
            )));
        }
    }
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
    // defined: no denominator identically zero, as written or after
    // D -> 4-2ε, whatever it depends on (ε given by the settings at the
    // entries, by its name at validation, else vakint's own: the test does not
    // depend on which, ε being a free variable there)
    let epsilon = epsilon
        .or_else(|| {
            epsilon_symbol_in(insertion.coefficient.as_view())
                .or_else(|| epsilon_symbol_in(insertion.factor.as_view()))
        })
        .unwrap_or_else(|| crate::utils::vakint_macros::vk_symbol!(ALLOWED_VAKINT_SYMBOL));
    for (what, atom) in [
        ("coefficient", &insertion.coefficient),
        ("factor", &insertion.factor),
    ] {
        if let Some((after, found)) = vanishing_denominator(atom.as_view(), dimension, epsilon) {
            let when = if after { " after D -> 4-2ε" } else { "" };
            return Err(error(format!(
                "a denominator vanishes identically{when} (undefined, not a pole of any order): {found} in the {what}"
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

/// The cap on the SUM of the exact pole orders of a record's insertions (the
/// closure entry raises the record's depth by it) (the only entry that raises a depth by it: the
/// direct, ports and vacuum entries hand the insertions to vakint inside the
/// integrand): 64. UNREACHABLE behind the depth cap of a marked value
/// ([`MARKED_MAX_TERMS`], 32: depth plus this raise), which refuses first;
/// kept as defence. Real records stay far below: the totals are at most 4 in
/// the producer's records (L - 1 at L loops: the five-loop L5_000121), with
/// at most 4 insertions per record (k5-stress). An R*
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

/// One caller expression outside the insertions (a record field or the
/// entry's prefactor argument) is defined and finite at ε = 0 after
/// D -> 4-2ε: `NotScalar` at `location`, "a denominator vanishes identically
/// .." or "a pole in ε outside an insertion ..".
pub(crate) fn check_finite_expression(
    atom: &Atom,
    location: &str,
    dimension: Symbol,
    epsilon: Symbol,
    poles: bool,
) -> Result<()> {
    let refuse = |detail: String| FactorizedError::NotScalar {
        location: location.to_string(),
        detail,
    };
    if let Some((after, found)) = vanishing_denominator(atom.as_view(), dimension, epsilon) {
        let when = if after { " after D -> 4-2ε" } else { "" };
        return Err(refuse(format!(
            "a denominator vanishes identically{when} (undefined, not a pole of any order): {found}"
        )));
    }
    // a pole: a candidate base vanishing at ε = 0 first (none in a real
    // record: no series then); for a candidate the WHOLE expression decides,
    // in lowest terms (a zero of a numerator cancels one of a denominator:
    // ε/(D-4) is finite), by the trailing power of its Laurent series in ε
    if !poles {
        return Ok(());
    }
    if let Some(candidate) = pole_denominator(atom.as_view(), dimension, epsilon) {
        let reduced = atom
            .replace(Atom::var(dimension))
            .with(Atom::num(4) - Atom::num(2) * Atom::var(epsilon))
            .together();
        let pole = match reduced.series(
            epsilon,
            Atom::Zero.as_view(),
            symbolica::poly::series::SeriesDepth::relative(1),
        ) {
            Ok(series) => !series.is_zero() && series.get_trailing_exponent().is_negative(),
            Err(e) => {
                return Err(refuse(format!(
                    "a pole in ε outside an insertion (no Laurent series at ε = 0: {e}): {candidate}"
                )));
            }
        };
        if pole {
            return Err(refuse(format!(
                "a pole in ε outside an insertion (only insertion coefficients may carry poles): {candidate}"
            )));
        }
    }
    Ok(())
}

/// The most ε terms a marked value may ask for: its depth (the record's
/// `depth`, or the settings' number of terms) is at least 1 and, with the
/// closure entry's insertion raise added, at most 32. A SANITY BOUND far above
/// the producer's maximum (12: six-loop K1 / K2 at depth 6, K3 with a
/// quotient's retries; the largest real depth + raise is 5 + 4), NOT a promise
/// that an evaluation at the cap is fast: the one-loop K2 takes 2.2 s at 16
/// terms (a 2.96 M character value) and did not finish in 200 s at 32. Below 1
/// the entries returned an empty (zero) value; far above, an evaluation that
/// does not end; beyond i64 the depth wrapped.
pub(crate) const MARKED_MAX_TERMS: i64 = 32;

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

/// The largest |exponent| a caller may write in a marked value (the real
/// records reach -8..9).
pub const MAX_CALLER_EXPONENT: i64 = 64;

/// The first power of `view` (at any depth, inside function arguments too)
/// whose exponent is not WRITTEN as an integer: after Symbolica's automatic
/// normalisation the exponent atom is not an integer number (`x^(2-3)` is
/// read as `x^-1` and accepted; an expression equal to an integer, `π`,
/// `sqrt(2)`, a symbol, a rational `1/2` are refused; nothing is expanded
/// to decide); with the exponent named.
/// The prefactor field a weight is checked against ([`WEIGHT_RULES`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum WeightField {
    Scalars,
    Overall,
    Colour,
    /// the prefactor's loop_sign (neutral 1)
    LoopSign,
    /// the prefactor's nf_power (neutral 0)
    NfPower,
}

/// THE WEIGHTS TABLE of the record format (section 9: `weight <name>
/// applied_by=beta|vakint <atom>`, the kinds colour, flavour, loop_sign,
/// overall and scalars; OQ-17: a weight the caller applies ("beta") stays out
/// of the prefactor, one vakint applies is in it): each weight name with ITS
/// OWN prefactor field. Applied by vakint, the field equals the weight;
/// applied by the caller, the field is neutral (1; nf_power 0). A name not in
/// the table, and a name listed twice, are refused.
pub(crate) const WEIGHT_RULES: [(&str, WeightField); 5] = [
    ("scalars", WeightField::Scalars),
    ("overall", WeightField::Overall),
    ("colour", WeightField::Colour),
    ("loop_sign", WeightField::LoopSign),
    ("flavour", WeightField::NfPower),
];

/// The largest numerator or denominator a marked value's numbers may have,
/// in bits: far above the real records (whose largest number has 35 bits);
/// it keeps a number that a text's literal power made huge at parse time
/// (2^1000000000 is a 10^9-bit integer once parsed) out of the evaluation.
pub const MAX_NUMBER_BITS: usize = 4096;

/// The first number of `view` whose numerator or denominator has more than
/// [`MAX_NUMBER_BITS`] bits, described by its decimal size.
fn first_oversized_number(view: AtomView<'_>) -> Option<String> {
    let mut found = None;
    view.visitor(&mut |part| {
        if found.is_some() {
            return false;
        }
        if let AtomView::Num(n) = part
            && let CoefficientView::Large(re, im) = n.get_coeff_view()
        {
            // the exact bit lengths of |numerator| and denominator, real and
            // imaginary parts
            for (what, r) in [("real", re.to_rat()), ("imaginary", im.to_rat())] {
                let bits = [
                    ("numerator", r.numerator_ref().significant_bits()),
                    ("denominator", r.denominator_ref().significant_bits()),
                ];
                if let Some((which, b)) = bits.iter().find(|(_, b)| *b > MAX_NUMBER_BITS as u64) {
                    found = Some(format!("a {what} part whose {which} has {b} bits"));
                    break;
                }
            }
        }
        found.is_none()
    });
    found
}

/// The first integer exponent of `view` beyond [`MAX_CALLER_EXPONENT`] in
/// modulus (the cheap part of [`first_non_integer_exponent`], for the rules
/// that compute and must not see such a power), shown with its power.
fn first_exponent_beyond_bound(view: AtomView<'_>) -> Option<String> {
    let mut found = None;
    view.visitor(&mut |part| {
        if found.is_some() {
            return false;
        }
        if let AtomView::Pow(p) = part {
            let (_, exponent) = p.get_base_exp();
            let beyond = match i64::try_from(exponent) {
                Ok(n) => n.unsigned_abs() > MAX_CALLER_EXPONENT.unsigned_abs(),
                Err(_) => symbolica::domains::integer::Integer::try_from(exponent).is_ok(),
            };
            if beyond {
                found = Some(format!(
                    "the exponent {} of {} exceeds {MAX_CALLER_EXPONENT} in modulus",
                    exponent.to_owned().to_canonical_string(),
                    part.to_owned().to_canonical_string()
                ));
            }
        }
        found.is_none()
    });
    found
}

fn first_non_integer_exponent(view: AtomView<'_>) -> Option<String> {
    let mut found = None;
    view.visitor(&mut |part| {
        if found.is_some() {
            return false;
        }
        if let AtomView::Pow(p) = part {
            let (_, exponent) = p.get_base_exp();
            match i64::try_from(exponent) {
                Ok(n) if n.unsigned_abs() <= MAX_CALLER_EXPONENT.unsigned_abs() => {}
                Ok(_) => {
                    found = Some(format!(
                        "the exponent {} of {} exceeds {MAX_CALLER_EXPONENT} in modulus",
                        exponent.to_owned().to_canonical_string(),
                        part.to_owned().to_canonical_string()
                    ))
                }
                Err(_) if symbolica::domains::integer::Integer::try_from(exponent).is_ok() => {
                    found = Some(format!(
                        "the exponent {} of {} exceeds {MAX_CALLER_EXPONENT} in modulus",
                        exponent.to_owned().to_canonical_string(),
                        part.to_owned().to_canonical_string()
                    ))
                }
                Err(_) => {
                    found = Some(format!(
                        "the exponent {} of {} is not an integer literal",
                        exponent.to_owned().to_canonical_string(),
                        part.to_owned().to_canonical_string()
                    ))
                }
            }
        }
        found.is_none()
    });
    found
}

/// The shape of one argument of a STRUCTURAL head (the table
/// [`STRUCTURAL_HEADS`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ArgShape {
    /// A non-negative integer literal of at most [`MAX_INDEX_ID`] (an id).
    Id,
    /// An integer literal of at most [`MAX_INDEX_ID`] in modulus (the compact
    /// slash's scope, `-1` for the projector).
    SignedId,
    /// A positive integer literal (an invariant's order).
    Positive,
    /// The record's dimension symbol (`mink`'s dimension).
    Dimension,
    /// The literal `4` (`bis`'s dimension).
    Four,
    /// A colour dimension: a positive integer literal, `spenso::dA` or
    /// `spenso::Nc`.
    ColourDimension,
    /// An index: an index label (`hedge`, `edge`, `vertex`, `dummy`, the
    /// compact slash), a plain caller symbol (a free index label such as
    /// `ir_dirac_closure_nu`) that is not ε, the dimension or a mass, or an
    /// integer literal id (the base grammar's integer indices).
    Index,
    /// A Lorentz or spinor slot `mink(D, index)` / `bis(4, index)`.
    Slot,
    /// `mink(D, index)` or the index-free `mink(D)` (a momentum's slot).
    MinkSlot,
    /// A momentum `Q(id, mink slot)`.
    Momentum,
    /// A Lorentz slot `mink(D, index)` or a momentum (the compact slash
    /// `gamma(a, b, Q(..))` before it is opened).
    SlotOrMomentum,
    /// An index-free colour representation `coad(dim)` / `cof(dim)`.
    ColourRep,
}

/// What a structural head is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum HeadRole {
    Slot,
    IndexLabel,
    Momentum,
    Metric,
    Gamma,
    Dot,
    ColourAtom,
    CompactSlash,
}

/// THE STRUCTURAL HEADS of a marked value: the only function heads besides
/// vakint's data atoms and `zeta(n)` ([`first_disallowed_head`]), each with
/// its role and its argument shapes by position (one list per allowed
/// arity). Index labels and the compact slash stand only as the index of a
/// slot ([`first_out_of_shape`]).
const STRUCTURAL_HEADS: &[(&str, HeadRole, &[&[ArgShape]])] = {
    use ArgShape::*;
    &[
        (
            names::MINK,
            HeadRole::Slot,
            &[&[Dimension, Index], &[Dimension]],
        ),
        (names::BIS, HeadRole::Slot, &[&[Four, Index]]),
        (names::HEDGE, HeadRole::IndexLabel, &[&[Id, Id]]),
        ("gammalooprs::edge", HeadRole::IndexLabel, &[&[Id, Id]]),
        ("gammalooprs::vertex", HeadRole::IndexLabel, &[&[Id, Id]]),
        ("gammalooprs::dummy", HeadRole::IndexLabel, &[&[Id]]),
        ("vakint::slash", HeadRole::CompactSlash, &[&[SignedId, Id]]),
        (names::MOMENTUM, HeadRole::Momentum, &[&[Id, MinkSlot]]),
        (names::METRIC, HeadRole::Metric, &[&[Slot, Slot]]),
        (
            names::GAMMA,
            HeadRole::Gamma,
            &[&[Slot, Slot, SlotOrMomentum]],
        ),
        (names::DOT, HeadRole::Dot, &[&[Momentum, Momentum]]),
        (
            "spenso::cas",
            HeadRole::ColourAtom,
            &[&[Positive, ColourRep]],
        ),
        (names::COAD, HeadRole::ColourAtom, &[&[ColourDimension]]),
        (names::COF, HeadRole::ColourAtom, &[&[ColourDimension]]),
    ]
};

/// The largest id an index label or a momentum may carry.
const MAX_INDEX_ID: i64 = u32::MAX as i64;

fn structural_head(name: &str) -> Option<(HeadRole, &'static [&'static [ArgShape]])> {
    STRUCTURAL_HEADS
        .iter()
        .find(|(n, _, _)| *n == name)
        .map(|(_, role, shapes)| (*role, *shapes))
}

/// The one special function a marked value may write: Symbolica's `zeta`
/// with ONE integer literal argument `n >= 2` (a number-like constant; the
/// producer writes `zeta(3)`).
fn is_zeta_constant(view: AtomView<'_>) -> bool {
    let AtomView::Fun(f) = view else {
        return false;
    };
    f.get_symbol().get_name() == "symbolica::zeta"
        && f.get_nargs() == 1
        && f.iter().all(|a| i64::try_from(a).is_ok_and(|n| n >= 2))
}

/// THE HEAD ALLOW-LIST of a marked value's caller atoms: every function head
/// (at any depth) is a structural head ([`STRUCTURAL_HEADS`]), one of
/// [`VAKINT_DATA_ATOMS`] in its exact shape (where it may stand is
/// [`misplaced_data_atom`]'s rule), or [`is_zeta_constant`]. Every other head
/// (Symbolica's built-ins `log`, `exp`, `gamma`, `sin`, …, `zeta` of another
/// argument, and every caller-namespaced function) is refused BY NAME,
/// whatever its argument: the first one found, named, with its atom.
fn first_disallowed_head(view: AtomView<'_>) -> Option<String> {
    let mut found = None;
    view.visitor(&mut |part| {
        if found.is_some() {
            return false;
        }
        if let AtomView::Fun(f) = part {
            if is_vakint_data_atom(part) || is_zeta_constant(part) {
                return false;
            }
            let head = f.get_symbol();
            let name = head.get_name();
            if structural_head(name).is_none() {
                found = Some(format!(
                    "{name} (in {}; allowed: the structural heads, vakint's data atoms, zeta(n) with an integer literal n >= 2)",
                    part.to_owned().to_canonical_string()
                ));
            }
        }
        found.is_none()
    });
    found
}

/// Whether the structural head `view` has one of its argument lists of
/// [`STRUCTURAL_HEADS`]; the first argument out of shape otherwise.
fn head_out_of_shape(
    view: AtomView<'_>,
    dimension: Symbol,
    depends: &dyn Fn(Symbol) -> bool,
) -> Option<String> {
    let AtomView::Fun(f) = view else {
        return None;
    };
    let head = f.get_symbol();
    let (_, shapes) = structural_head(head.get_name())?;
    let args = f.iter().collect::<Vec<_>>();
    let Some(shape) = shapes.iter().find(|s| s.len() == args.len()) else {
        return Some(format!(
            "{} ({} arguments; expected {})",
            view.to_owned().to_canonical_string(),
            args.len(),
            shapes
                .iter()
                .map(|s| s.len().to_string())
                .collect::<Vec<_>>()
                .join(" or ")
        ));
    };
    for (i, (a, s)) in args.iter().zip(shape.iter()).enumerate() {
        if let Some(why) = argument_out_of_shape(*a, *s, dimension, depends) {
            return Some(format!(
                "{} (argument {}: {why})",
                view.to_owned().to_canonical_string(),
                i + 1
            ));
        }
    }
    None
}

/// Whether `a` has the shape `shape` (a head inside it checked against its
/// own shapes); why not otherwise.
fn argument_out_of_shape(
    a: AtomView<'_>,
    shape: ArgShape,
    dimension: Symbol,
    depends: &dyn Fn(Symbol) -> bool,
) -> Option<String> {
    use ArgShape::*;
    let shown = || a.to_owned().to_canonical_string();
    let head_role = |a: AtomView<'_>| match a {
        AtomView::Fun(f) => structural_head(f.get_symbol().get_name())
            .map(|(r, _)| (r, f.get_symbol().get_name() == names::MINK, f.get_nargs())),
        _ => None,
    };
    let fits = match shape {
        Id => i64::try_from(a).is_ok_and(|n| (0..=MAX_INDEX_ID).contains(&n)),
        SignedId => i64::try_from(a).is_ok_and(|n| n.unsigned_abs() <= MAX_INDEX_ID.unsigned_abs()),
        Positive => i64::try_from(a).is_ok_and(|n| n >= 1),
        Dimension => matches!(a, AtomView::Var(v) if v.get_symbol() == dimension),
        Four => i64::try_from(a) == Ok(4),
        ColourDimension => {
            i64::try_from(a).is_ok_and(|n| n >= 1)
                || matches!(a, AtomView::Var(v) if matches!(v.get_symbol().get_name(), "spenso::dA" | "spenso::Nc"))
        }
        Index => match a {
            AtomView::Var(v) => !depends(v.get_symbol()) && v.get_symbol() != dimension,
            AtomView::Num(_) => i64::try_from(a).is_ok_and(|n| (0..=MAX_INDEX_ID).contains(&n)),
            _ => matches!(
                head_role(a),
                Some((HeadRole::IndexLabel | HeadRole::CompactSlash, _, _))
            ),
        },
        Slot => matches!(head_role(a), Some((HeadRole::Slot, _, 2))),
        MinkSlot => matches!(head_role(a), Some((HeadRole::Slot, true, _))),
        Momentum => matches!(head_role(a), Some((HeadRole::Momentum, _, _))),
        SlotOrMomentum => matches!(
            head_role(a),
            Some((HeadRole::Slot, true, 2) | (HeadRole::Momentum, _, _))
        ),
        ColourRep => matches!(head_role(a), Some((HeadRole::ColourAtom, _, 1))),
    };
    if !fits {
        return Some(format!("{} is not {shape:?}", shown()));
    }
    head_out_of_shape(a, dimension, depends)
}

/// The first STRUCTURAL head of `view` out of its shape
/// ([`STRUCTURAL_HEADS`]: an argument that is not of the stated shape, e.g.
/// `hedge(ε, 1)`, `mink(D, ε)`, `vertex(x, 1)`, `dummy(1/2)`, `Q(1, x)`), or
/// an index label or the compact slash outside a slot (a free factor): with
/// the reason. Heads that are not structural are [`first_disallowed_head`]'s.
fn first_out_of_shape(
    view: AtomView<'_>,
    dimension: Symbol,
    depends: &dyn Fn(Symbol) -> bool,
) -> Option<String> {
    match view {
        AtomView::Add(a) => a
            .iter()
            .find_map(|t| first_out_of_shape(t, dimension, depends)),
        AtomView::Mul(m) => m
            .iter()
            .find_map(|t| first_out_of_shape(t, dimension, depends)),
        AtomView::Pow(p) => first_out_of_shape(p.get_base_exp().0, dimension, depends),
        AtomView::Fun(f) => match structural_head(f.get_symbol().get_name()) {
            Some((HeadRole::IndexLabel | HeadRole::CompactSlash, _)) => Some(format!(
                "{} (an index label outside a slot)",
                view.to_owned().to_canonical_string()
            )),
            Some(_) => head_out_of_shape(view, dimension, depends),
            None => None,
        },
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
        // a closure value's provenance is its ClosureSpec: a direct record
        // next to it (typed only: the text form cannot write both) is an
        // inconsistent pair
        if matches!(self.kind, FactorizedKind::Closure(_)) && self.record.is_some() {
            return Err(FactorizedError::Closure {
                detail: "a closure value carries a direct record (record must be None: a closure's provenance is its ClosureSpec)".into(),
            });
        }
        // the record kinds of the R* interface may not name vakint's own
        // symbols (the entries' internal placeholders live there)
        if is_marked(self) {
            self.check_reserved_symbols()?;
        }
        self.check_marked_bounds(false)?;
        if is_marked(self) {
            self.check_marked_prefactor()?;
            self.check_paired_legs()?;
            // a projector power of p^2 needs an external momentum (the
            // closure kind's rule, for the direct kind too: with no external
            // p, (p^2)^n has no value)
            if matches!(self.kind, FactorizedKind::Direct)
                && self.projector.as_ref().is_some_and(|p| p.p2_power != 0)
                && self.externals() == 0
            {
                return Err(FactorizedError::Graph {
                    detail: "p2_power without an external momentum".into(),
                });
            }
        }
        // the contract / closure rules; their insertion rules compute (zero
        // tests, together / expand, series) only after the cheap bounds they
        // depend on (exponents of at most 64, numbers of at most 4096 bits:
        // insertion_epsilon_rules); the remaining atom rules after them
        if let FactorizedKind::Closure(spec) = &self.kind {
            self.validate_closure(spec)?;
        } else {
            self.validate_contract(None)?;
        }
        if is_marked(self) {
            self.check_marked_atoms()?;
            // every caller expression outside the insertions: defined, finite
            self.check_marked_fields(None)?;
        }
        Ok(())
    }

    /// The mutable twin of [`Self::record_atoms`] (tests: put a marker into
    /// every caller-written atom), destructured the same way.
    #[cfg(test)]
    pub(crate) fn record_atoms_mut(&mut self) -> Vec<(String, &mut Atom)> {
        let FactorizedNumerator {
            label: _,
            kind,
            dimension: _,
            loop_edges: _,
            external_edges: _,
            edges,
            vertices,
            ports,
            projector,
            prefactor,
            colour_branch: _,
            insertions,
            record,
        } = self;
        let mut out = Vec::new();
        for edge in edges {
            let FactorizedEdge {
                id,
                source: _,
                target: _,
                source_hedge: _,
                target_hedge: _,
                particle: _,
                line: _,
                power: _,
                mass_squared,
                momentum: _,
                factor,
            } = edge;
            out.push((format!("edge {id} factor"), factor));
            out.push((format!("edge {id} mass"), mass_squared));
        }
        for vertex in vertices {
            let FactorizedVertex {
                id,
                rule: _,
                factor,
            } = vertex;
            out.push((format!("vertex {id} factor"), factor));
        }
        for port in ports {
            out.push(("ports".to_string(), port));
        }
        if let Some(Projector {
            tensor,
            p2_power: _,
        }) = projector
        {
            out.push(("projector".to_string(), tensor));
        }
        let Prefactor {
            overall,
            loop_sign: _,
            nf_power: _,
            scalars,
            colour,
        } = prefactor;
        out.push(("prefactor_overall".to_string(), overall));
        out.push(("prefactor_scalars".to_string(), scalars));
        out.push(("prefactor_colour".to_string(), colour));
        for insertion in insertions {
            let Insertion {
                vertex,
                hedges: _,
                cycles: _,
                coefficient,
                factor,
            } = insertion;
            out.push((format!("insertion {vertex} coefficient"), coefficient));
            out.push((format!("insertion {vertex} factor"), factor));
        }
        let weights = match kind {
            FactorizedKind::Closure(spec) => {
                let ClosureSpec {
                    edge: _,
                    mass_squared,
                    auxiliary_mass_squared,
                    power: _,
                    independent_externals: _,
                    taylor: _,
                    softs: _,
                    original_cycles: _,
                    depth: _,
                    weights,
                } = spec;
                out.push(("closure mass".to_string(), mass_squared));
                out.push(("taylor auxiliary mass".to_string(), auxiliary_mass_squared));
                Some(weights)
            }
            FactorizedKind::Direct => None,
        };
        for (name, _, atom) in weights.into_iter().flatten() {
            out.push((format!("weight {name}"), atom));
        }
        if let Some(record) = record {
            for (name, _, atom) in record.weights.iter_mut() {
                out.push((format!("weight {name}"), atom));
            }
        }
        out
    }

    /// THE WEIGHTS of a marked value, by [`WEIGHT_RULES`] (the record format's
    /// weights rule, section
    /// 3), for every kind and both weight sources (a closure's ClosureSpec,
    /// a direct record): a weight `applied_by=vakint` must equal
    /// prefactor_scalars (atom equality after normalisation); a weight applied
    /// by the caller is 1 in the prefactor: `overall` -> prefactor_overall 1,
    /// `colour` -> prefactor_colour 1, `loop_sign` -> loop_sign 1, `flavour`
    /// -> nf_power 0 (a caller weight of another name has no counterpart).
    /// And prefactor_scalars carries no flavour symbol (`nf`, `…::nf`: n_f
    /// stays with the caller; section 3). `Closure` for the closure kind, else
    /// `Graph`. Shared by `validate` and the vacuum mode's validation.
    pub(crate) fn check_marked_prefactor(&self) -> Result<()> {
        let closure = matches!(self.kind, FactorizedKind::Closure(_));
        let refuse = |detail: String| {
            if closure {
                FactorizedError::Closure { detail }
            } else {
                FactorizedError::Graph { detail }
            }
        };
        let p = &self.prefactor;
        let mut flavour = false;
        p.scalars.as_view().visitor(&mut |part| {
            if let AtomView::Var(v) = part {
                let symbol = v.get_symbol();
                let name = symbol.get_name();
                flavour |= name == names::NF || name.ends_with("::nf");
            }
            !flavour
        });
        if flavour {
            return Err(refuse(format!(
                "a flavour symbol in the record (n_f stays with the caller): in prefactor_scalars {}",
                p.scalars
            )));
        }
        let from_spec = match &self.kind {
            FactorizedKind::Closure(spec) => spec.weights.as_slice(),
            FactorizedKind::Direct => &[],
        };
        let from_record = self
            .record
            .as_ref()
            .map_or(&[][..], |r| r.weights.as_slice());
        let mut seen = std::collections::BTreeSet::new();
        for (name, by_vakint, atom) in from_spec.iter().chain(from_record) {
            let applied = if *by_vakint { "vakint" } else { "beta" };
            if !seen.insert(name.as_str()) {
                return Err(refuse(format!(
                    "weight {name} is listed twice (each weight is applied once)"
                )));
            }
            let Some((_, field)) = WEIGHT_RULES.iter().find(|(n, _)| *n == name.as_str()) else {
                return Err(refuse(format!(
                    "weight {name} applied_by={applied}: not a weight of the record format (scalars, overall, colour, loop_sign, flavour)"
                )));
            };
            // the field's value as an atom, and its neutral value
            let (label, value, neutral) = match field {
                WeightField::Scalars => ("prefactor_scalars", p.scalars.clone(), Atom::num(1)),
                WeightField::Overall => ("prefactor_overall", p.overall.clone(), Atom::num(1)),
                WeightField::Colour => ("prefactor_colour", p.colour.clone(), Atom::num(1)),
                WeightField::LoopSign => ("loop_sign", Atom::num(p.loop_sign), Atom::num(1)),
                WeightField::NfPower => {
                    ("nf_power", Atom::num(i64::from(p.nf_power)), Atom::num(0))
                }
            };
            let (expected, why) = if *by_vakint {
                (
                    atom.clone(),
                    "a weight vakint applies is its prefactor field",
                )
            } else {
                (
                    neutral,
                    "a weight the caller applies leaves its prefactor field neutral",
                )
            };
            if value != expected {
                return Err(refuse(format!(
                    "weight {name} applied_by={applied} {atom}: {label} is {value}, not {expected} ({why})"
                )));
            }
        }
        Ok(())
    }

    /// ORIGINAL CYCLES of a marked DIRECT record (the record format, section
    /// 3, "every R* record"), for the vacuum mode, which does not run the
    /// contract rules (`validate_contract` checks the same for the direct and
    /// ports entries): with `original_cycles` the prefactor has loop_sign 1,
    /// nf_power 0, prefactor_overall 1 (the caller applies the original
    /// graph's weights), and the carried graph's closed cycles plus the
    /// insertions' do not exceed the original's; without `original_cycles`
    /// the loop sign and the n_f power are those of the carried graph's closed
    /// loops (with the insertions' cycles). The text form of the same value
    /// runs these rules through `validate`; the typed vacuum entry runs them
    /// here, so the two cannot disagree.
    pub(crate) fn check_vacuum_original_cycles(&self) -> Result<()> {
        let p = &self.prefactor;
        let Some(DirectRecord {
            original_cycles: Some((fermion, ghost)),
            ..
        }) = &self.record
        else {
            // without original_cycles: the carried graph's closed loops (the
            // contract's rule, which the text form of the same value runs)
            let loops = self.closed_loops()?;
            let (inner_fermion, inner_ghost) =
                self.insertions.iter().fold((0usize, 0usize), |(f, g), i| {
                    (f + i.cycles.0 as usize, g + i.cycles.1 as usize)
                });
            let n_loops = loops.len() + inner_fermion + inner_ghost;
            let n_fermion =
                loops.iter().filter(|l| l.line == LineKind::Fermion).count() + inner_fermion;
            let sign = if n_loops % 2 == 0 { 1 } else { -1 };
            if p.loop_sign != sign {
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
            return Ok(());
        };
        if p.loop_sign != 1 || p.nf_power != 0 || p.overall != Atom::num(1) {
            return Err(FactorizedError::NotScalar {
                location: "prefactor".into(),
                detail: "with original_cycles the weights are the caller's: loop_sign=1, \
                         nf_power=0, prefactor_overall 1"
                    .into(),
            });
        }
        let loops = self.closed_loops()?;
        let (inner_fermion, inner_ghost) =
            self.insertions.iter().fold((0usize, 0usize), |(f, g), i| {
                (f + i.cycles.0 as usize, g + i.cycles.1 as usize)
            });
        let n_fermion =
            loops.iter().filter(|l| l.line == LineKind::Fermion).count() + inner_fermion;
        let n_ghost = loops.len() + inner_fermion + inner_ghost - n_fermion;
        if n_fermion > *fermion as usize || n_ghost > *ghost as usize {
            return Err(FactorizedError::Graph {
                detail: format!(
                    "the record and its insertions carry {n_fermion} fermion / {n_ghost} ghost cycles, more than \
                     original_cycles fermion={fermion} ghost={ghost}"
                ),
            });
        }
        Ok(())
    }

    /// PAIRED LEGS (the record format, section 3: a crown line with both ends in
    /// the scope is two legs, `e` and `100000 + e`): when a marked record has
    /// an external edge `100000 + e` and an external edge `e`, the two carry
    /// exactly the same momentum vector (each leg's Q is routed through the
    /// momentum it declares: unequal ones evaluate to another value). `Closure`
    /// for the closure kind, else `Graph`.
    pub(crate) fn check_paired_legs(&self) -> Result<()> {
        let external = |id: usize| self.edges.iter().find(|e| e.id == id && e.is_external());
        for leg in self
            .edges
            .iter()
            .filter(|e| e.is_external() && e.id >= 100_000)
        {
            if let Some(partner) = external(leg.id - 100_000)
                && partner.momentum != leg.momentum
            {
                let (id, a, other, b) = (leg.id, &leg.momentum, partner.id, &partner.momentum);
                return Err(if matches!(self.kind, FactorizedKind::Closure(_)) {
                    FactorizedError::Closure {
                        detail: format!(
                            "the two legs of a paired line carry one momentum: edge {id} carries {a:?}, its partner {other} {b:?}"
                        ),
                    }
                } else {
                    FactorizedError::Graph {
                        detail: format!(
                            "the two legs of a paired line carry one momentum: edge {id} carries {a:?}, its partner {other} {b:?}"
                        ),
                    }
                });
            }
        }
        Ok(())
    }

    /// THE BOUNDS of a value with new features, shared by `validate` and the
    /// vacuum mode's validation (which never calls `validate`): the record's
    /// depth ([`marked_depth`]), the ids (i64: written as Q(id, ..) and
    /// edge(source, target); the closure kind bounds its edge ids itself, in
    /// `validate_closure`, here its vertex ids), the edge powers (summed when
    /// lines are merged in series: at most 2^31 in absolute value) and, with
    /// `momenta` (the vacuum mode: `validate` checks them in the contract's
    /// momentum rules), the momentum coefficients (at most 2^31). `Closure`
    /// for the closure kind, else `Graph` (`MomentumBasis` for coefficients).
    pub(crate) fn check_marked_bounds(&self, momenta: bool) -> Result<()> {
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
        if let Some(detail) = self.marked_id_overflow(!closure) {
            return Err(in_kind(detail));
        }
        if let Some(edge) = self.edges.iter().find(|e| e.power.unsigned_abs() > 1 << 31) {
            return Err(in_kind(format!(
                "edge {}: power {} out of range (at most 2^31 in absolute value)",
                edge.id, edge.power
            )));
        }
        if momenta {
            for edge in &self.edges {
                if let Some(c) = edge.momentum.iter().find(|c| c.unsigned_abs() > 1 << 31) {
                    return Err(FactorizedError::MomentumBasis {
                        detail: format!(
                            "edge {}: momentum coefficient {c} out of range (at most 2^31 in absolute value)",
                            edge.id
                        ),
                    });
                }
            }
        }
        Ok(())
    }

    /// The atom rules of a MARKED value beyond the namespace walk, for
    /// `validate` and the vacuum mode alike: every non-mass record atom
    /// through [`Self::check_opaque_dependence`] (exponents, opaque function
    /// arguments, data atoms, data-atom linearity), and every nonzero edge mass
    /// a plain symbol (`Graph`, prefix `edge <id>: mass `).
    pub(crate) fn check_marked_atoms(&self) -> Result<()> {
        for (field, atom) in self.record_atoms() {
            if !field.ends_with(" mass") {
                self.check_opaque_dependence(atom, &field)?;
            }
        }
        if let Some(edge) = self.edges.iter().find(|e| {
            !e.mass_squared.is_zero() && !matches!(e.mass_squared.as_view(), AtomView::Var(_))
        }) {
            return Err(FactorizedError::Graph {
                detail: format!(
                    "edge {}: mass {} is neither 0 nor a plain symbol",
                    edge.id, edge.mass_squared
                ),
            });
        }
        Ok(())
    }

    /// EXACT NUMBERS, EXPONENTS, DATA ATOMS AND FUNCTION HEADS (marked
    /// values): vakint's data atoms are linear factors ([`data_atom_degree`],
    /// prefix `a non-linear vakint data atom in <field>: `); every
    /// number is exact and finite ([`check_exact_numbers`], `InexactNumber`
    /// at `<field>`: Symbolica writes a literal singular form, `zeta(1)`,
    /// `log(0)`, `gamma(0)`, as a non-finite number); every exponent is
    /// written as an integer ([`first_non_integer_exponent`]: detail prefix
    /// `a non-integer exponent in <field>: `); every function head is a
    /// structural head, a vakint data atom or `zeta(n)` with an integer
    /// literal `n >= 2` ([`first_disallowed_head`]: detail prefix `a function
    /// head that is not allowed in <field>: <head>`); every structural head
    /// has its arguments of the shapes of [`STRUCTURAL_HEADS`], an index label
    /// only inside a slot ([`first_out_of_shape`]: detail prefix `a structural
    /// head out of its shape in <field>: `). `Graph`, `Insertion` in an
    /// insertion. Caller-namespaced
    /// SCALARS (plain symbols) stay accepted as opaque (the format has no
    /// declared symbol table); caller-namespaced FUNCTIONS are refused.
    pub(crate) fn check_opaque_dependence(&self, atom: &Atom, field: &str) -> Result<()> {
        // exact numbers in EVERY field this walk visits (the weights too: a
        // literal singular form such as zeta(1) or log(0) is a non-finite
        // number after Symbolica's normalisation)
        check_exact_numbers(atom.as_view(), field)?;
        let masses = match &self.kind {
            FactorizedKind::Closure(spec) => [&spec.mass_squared, &spec.auxiliary_mass_squared]
                .into_iter()
                .filter_map(|m| match m.as_view() {
                    AtomView::Var(v) => Some(v.get_symbol()),
                    _ => None,
                })
                .collect::<Vec<_>>(),
            _ => Vec::new(),
        };
        let dimension = self.dimension;
        let depends = |symbol: Symbol| {
            let name = symbol.get_name();
            name == "ε" || name.ends_with("::ε") || symbol == dimension || masses.contains(&symbol)
        };
        let found = misplaced_data_atom(atom.as_view(), false)
            .map(|found| {
                format!("a misplaced vakint data atom in {field}: {found} (only as a factor of a coefficient, not inside a function or an exponent)")
            })
            .or_else(|| {
                first_non_integer_exponent(atom.as_view())
                    .map(|found| format!("a non-integer exponent in {field}: {found}"))
            })
            .or_else(|| match data_atom_degree(atom.as_view()) {
                Some(d) if d <= 1 => None,
                _ => Some(format!(
                    "a non-linear vakint data atom in {field}: vakint's data atoms (UNKNOWN(master, k), MZV(..)) are carried as linear factors (a product of two, a power, a denominator, or one inside the base of a power is refused): {}",
                    atom.to_canonical_string()
                )),
            })
            .or_else(|| {
                first_disallowed_head(atom.as_view())
                    .map(|found| format!("a function head that is not allowed in {field}: {found}"))
            })
            .or_else(|| {
                first_oversized_number(atom.as_view()).map(|found| {
                    format!("a number beyond {MAX_NUMBER_BITS} bits in {field}: {found}")
                })
            })
            .or_else(|| {
                first_out_of_shape(atom.as_view(), dimension, &depends).map(|found| {
                    format!("a structural head out of its shape in {field}: {found}")
                })
            });
        match found {
            None => Ok(()),
            Some(detail) => {
                let insertion = field
                    .strip_prefix("insertion ")
                    .and_then(|rest| rest.split(' ').next())
                    .and_then(|v| v.parse::<usize>().ok());
                Err(match insertion {
                    Some(vertex) => FactorizedError::Insertion { vertex, detail },
                    None => FactorizedError::Graph { detail },
                })
            }
        }
    }

    /// Every caller expression of a marked value outside its insertions (every
    /// atom of [`Self::record_atoms`] but the insertions' coefficients and
    /// factors: edge factors and masses, vertex factors, ports, the projector,
    /// the prefactor fields, the closure line's masses, the WEIGHTS) is defined
    /// and finite at ε = 0 once the dimension is set to 4 - 2ε: no denominator
    /// vanishes identically, as written or after the substitution, and none
    /// makes a pole (poles belong in insertion coefficients, whose pole order
    /// raises the depth; R1 counts nothing else): `NotScalar` at the field. ε
    /// is `epsilon` when given (every entry, with the settings' ε), else (at
    /// validation, without settings) the first symbol named ε in the
    /// insertions (the record's own ε); without insertions validation does
    /// not know the record's ε: a field naming a symbol ε gets only the zero
    /// test there, its pole test at the entry (with the settings' ε: a
    /// lookalike, such as symbolica::ε, is then a plain symbol, refused by the
    /// entry's name rule); a field written in the dimension alone gets both.
    pub(crate) fn check_marked_fields(&self, epsilon: Option<Symbol>) -> Result<()> {
        let fields = self
            .record_atoms()
            .into_iter()
            .filter(|(field, _)| !field.starts_with("insertion "))
            .collect::<Vec<_>>();
        let known = epsilon.or_else(|| {
            self.insertions.iter().find_map(|i| {
                epsilon_symbol_in(i.coefficient.as_view())
                    .or_else(|| epsilon_symbol_in(i.factor.as_view()))
            })
        });
        let epsilon =
            known.unwrap_or_else(|| crate::utils::vakint_macros::vk_symbol!(ALLOWED_VAKINT_SYMBOL));
        for (location, atom) in fields {
            let poles = known.is_some() || epsilon_symbol_in(atom.as_view()).is_none();
            check_finite_expression(atom, &location, self.dimension, epsilon, poles)?;
        }
        Ok(())
    }

    /// A marked value's ids must be i64 (they are written as Q(id, ..),
    /// edge(source, target), with fresh loop-vector ids max + 1 + i above the
    /// largest): the detail of the refusal when one is not.
    fn marked_id_overflow(&self, edges: bool) -> Option<String> {
        // the closure kind bounds its edge ids itself (validate_closure, with
        // the fresh loop-vector ids): here only its vertex ids
        let fresh = if edges {
            self.loop_edges
                .iter()
                .rposition(|e| e.is_none())
                .map_or(0, |i| i + 1)
        } else {
            0
        };
        let top = self
            .edges
            .iter()
            .flat_map(|e| {
                let vertex = |end: &Endpoint| match end {
                    Endpoint::Vertex(v) => Some(*v),
                    Endpoint::External => None,
                };
                [edges.then_some(e.id), vertex(&e.source), vertex(&e.target)]
                    .into_iter()
                    .flatten()
            })
            .chain(self.vertices.iter().map(|v| v.id))
            .chain(self.constant_vectors().into_iter().filter(|_| edges))
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

    /// EVERY caller-written atom of the record with its field name: edge
    /// factors and masses, vertex factors, ports, the projector, the three
    /// prefactor fields, insertion coefficients and factors, the closure
    /// line's masses, the weights (a closure's or a v3 direct record's). The
    /// structs are destructured WITHOUT `..`: a new field does not compile
    /// until it is listed here (an atom) or named as no atom, so that every
    /// rule that walks this list (the atom rules, the namespace walk's
    /// counterpart, the defined-and-finite check) reaches it.
    pub(crate) fn record_atoms(&self) -> Vec<(String, &Atom)> {
        let FactorizedNumerator {
            label: _,
            kind,
            dimension: _,
            loop_edges: _,
            external_edges: _,
            edges,
            vertices,
            ports,
            projector,
            prefactor,
            colour_branch: _,
            insertions,
            record,
        } = self;
        let mut out = Vec::new();
        for edge in edges {
            let FactorizedEdge {
                id,
                source: _,
                target: _,
                source_hedge: _,
                target_hedge: _,
                particle: _,
                line: _,
                power: _,
                mass_squared,
                momentum: _,
                factor,
            } = edge;
            out.push((format!("edge {id} factor"), factor));
            out.push((format!("edge {id} mass"), mass_squared));
        }
        for vertex in vertices {
            let FactorizedVertex {
                id,
                rule: _,
                factor,
            } = vertex;
            out.push((format!("vertex {id} factor"), factor));
        }
        for port in ports {
            out.push(("ports".to_string(), port));
        }
        if let Some(Projector {
            tensor,
            p2_power: _,
        }) = projector
        {
            out.push(("projector".to_string(), tensor));
        }
        let Prefactor {
            overall,
            loop_sign: _,
            nf_power: _,
            scalars,
            colour,
        } = prefactor;
        out.push(("prefactor_overall".to_string(), overall));
        out.push(("prefactor_scalars".to_string(), scalars));
        out.push(("prefactor_colour".to_string(), colour));
        for insertion in insertions {
            let Insertion {
                vertex,
                hedges: _,
                cycles: _,
                coefficient,
                factor,
            } = insertion;
            out.push((format!("insertion {vertex} coefficient"), coefficient));
            out.push((format!("insertion {vertex} factor"), factor));
        }
        let weights = match (kind, record) {
            (FactorizedKind::Closure(spec), _) => {
                let ClosureSpec {
                    edge: _,
                    mass_squared,
                    auxiliary_mass_squared,
                    power: _,
                    independent_externals: _,
                    taylor: _,
                    softs: _,
                    original_cycles: _,
                    depth: _,
                    weights,
                } = spec;
                out.push(("closure mass".to_string(), mass_squared));
                out.push(("taylor auxiliary mass".to_string(), auxiliary_mass_squared));
                weights.as_slice()
            }
            (FactorizedKind::Direct, _) => &[],
        };
        for (name, _, atom) in weights {
            out.push((format!("weight {name}"), atom));
        }
        // the direct record's weights whatever the kind (a closure value with
        // a record is refused by validate; no atom is skipped because of it)
        if let Some(DirectRecord {
            original_cycles: _,
            depth: _,
            weights,
        }) = record
        {
            for (name, _, atom) in weights {
                out.push((format!("weight {name}"), atom));
            }
        }
        out
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
            for (atom, part) in [
                (&insertion.coefficient, "coefficient"),
                (&insertion.factor, "factor"),
            ] {
                check_caller_atom_with_data_atoms(
                    atom,
                    &format!("insertion {} {part}", insertion.vertex),
                    exempt,
                    AtomSource::Record,
                )?;
            }
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
            insertion_epsilon_rules(insertion, &depends, dimension, None)?;
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
            let settings_epsilon = crate::utils::vakint_macros::vk_symbol!(epsilon);
            for insertion in &self.insertions {
                insertion_epsilon_rules(insertion, &depends, dimension, Some(settings_epsilon))?;
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

    /// The versioned text form ([`FACTORIZED_SCHEMA`]). A MARKED value is
    /// validated first: a value validation refuses is not written, and the
    /// error is EXACTLY validation's (the same variant and detail), so that a
    /// caller can route it like the reader's.
    pub fn to_text(&self) -> Result<String> {
        text::write_record(self)
    }

    /// The text form WITHOUT validation, to keep the text of a refused value
    /// for diagnosis: NOT evidence of validity (the reader may refuse it, and
    /// for a value validation refuses it may not be that value's text: the
    /// writer sorts the insertions and drops what the format cannot express);
    /// for a valid value the same text as [`Self::to_text`].
    pub fn to_text_unvalidated_for_diagnostics(&self) -> Result<String> {
        text::write_record_unvalidated(self)
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
