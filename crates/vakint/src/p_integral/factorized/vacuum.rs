//! The uncontracted-numerator input mode for SINGLE-SCALE VACUUM integrals
//! (helper VH, numerator assignment, lane rstar-kinds): the vacuum sibling of
//! [`super::evaluate_factorized`] (which takes massless propagator graphs).
//!
//! A [`FactorizedNumerator`] whose internal edges carry loop momenta only
//! (no external momentum flows through a denominator) and massive lines of
//! one mass scale; its numerator is the product of its spenso factors (gamma
//! chains, metrics, momenta), optionally closed by the caller's projector.
//! OPEN INDICES are the declared `ports` (`mink(D, ..)` slots of the factor
//! product); EXTERNAL VECTORS are the independent external momenta `p(j)`,
//! written in the factors as `Q(x, ..)` of the external leg `x` with momentum
//! `p(j)`, and the CONSTANT external vectors `Q(x, ..)` of edges `x` that are
//! not in the graph (reader rule 2; after the legs' `p(j)`, in id order). The
//! same declaration serves the propagator case (lane d6) and the R* hard parts
//! (D5 section 14), in the leg form and in the flag-free form (the hard part as
//! a vacuum graph of its own lines: slots of its absent boundary half-edges
//! are ports, reader rule 1 of the FeynKit ingress and of `validate`).
//!
//! 1. CONTRACT with idenso: `projector · PROD factors` through
//!    [`super::contract_with_idenso`] (cook, `infer_indexed`,
//!    `simplify_algebra{gamma canonical, epsilon}`, `to_dots`, uncook): every
//!    trace, chain and contraction among the factors; what is left is the
//!    scope vakint's vacuum evaluation always supported: scalar products,
//!    momenta with an open index and metrics between open indices.
//! 2. MAP to vakint's vacuum notation (linear substitution): an edge dot
//!    becomes `dot` of the edges' momenta in `k(1..L)` and `p(j)`, an edge
//!    momentum with open port i becomes the same combination with the raw
//!    index `PORT + i` (`k(l, PORT + i)`, `p(j, PORT + i)`), a metric between
//!    ports `g(PORT + i, PORT + j)`; anything else is an error.
//! 3. EVALUATE, by the back end the graph's lines select ([`vacuum_backend`],
//!    after series merging, [`vacuum_lines`]):
//!    - EVERY line massive (one mass): vakint's EXISTING vacuum path
//!      (`Vakint::evaluate`): its vacuum tensor reduction (method of the
//!      settings: FeynKit by default, feynkit `TensorReducer`) projects the
//!      loop-momentum tensors onto tensors of purely external objects, which
//!      ride along as constant coefficients; matching, the backend of the
//!      evaluation order (RustRed first by default: one to three loops from
//!      its vendored artifacts, four with its external data, five with an
//!      application package; FORM-based backends only where the order names
//!      them) and the finaliser.
//!    - ONE massive line among massless ones (the R* hard part's shape): the
//!      same vacuum tensor reduction (`Vakint::tensor_reduce`, which takes
//!      mixed-mass vacua), then per external tensor its scalar coefficient
//!      through vakint's MASSIVE CLOSURE ([`MassiveClosure`]): the massive line
//!      is cut, its momentum becomes the head `p(1)` of the massless
//!      p-integral that remains ([`VacuumCut`]: routing and cut graph derived
//!      here from the lines), which the closure reduces (its p backend) and
//!      closes with the massive line.
//!    - one massive line that no massless line depends on: a massive tadpole
//!      times a scaleless massless vacuum, exactly zero.
//! 4. The value keeps its OPEN indices; it is returned in vakint notation and
//!    mapped back to the caller's notation (`g(slot, slot)`, `Q(x, slot)`,
//!    `spenso::dot(Q(x, mink(D)), Q(y, mink(D)))`).

use std::time::Instant;

use spenso::network::{library::symbolic::ETS, tags::SPENSO_TAG};
use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;

use super::{Endpoint, FactorizedError, FactorizedKind, FactorizedNumerator, Result, names};
use crate::p_integral::{DirectClosureInput, MassiveClosure};
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{Vakint, VakintExpression, VakintSettings, VakintTerm};

/// First raw index of the open ports in vakint notation (`PORT + i`).
pub const PORT: i64 = 9001;

pub(super) fn error(detail: impl Into<String>) -> FactorizedError {
    FactorizedError::Evaluation {
        detail: format!("vacuum: {}", detail.into()),
    }
}

/// A value of [`evaluate_factorized_vacuum`] and what it was made of.
#[derive(Clone, Debug)]
pub struct VacuumValue {
    /// The Laurent expansion with the open ports and external vectors in the
    /// caller's notation.
    pub value: Atom,
    /// The same in vakint notation (raw port indices `PORT + i`, `p(j)`).
    pub vakint_value: Atom,
    /// The integrand (in epsilon): handed to `Vakint::evaluate` (family) or
    /// to `Vakint::tensor_reduce` (closure).
    pub input: Atom,
    /// The back end the lines selected.
    pub backend: VacuumBackend,
    /// The size predictor of the uncontracted network ([`distributed_lorentz_monomials`]).
    pub predicted_monomials: f64,
    /// Pieces of the streamed contraction (0: unstreamed).
    pub pieces: usize,
    /// Wall seconds of the idenso contraction and of the vakint evaluation.
    pub contract_seconds: f64,
    pub evaluate_seconds: f64,
}

/// The SIZE PREDICTOR of an uncontracted network (D5-HARD-KIND 8.2): the
/// number of Lorentz monomials a full distribution of it would create, read off
/// the expression tree in one pass without expanding: a sum that carries a
/// Minkowski slot adds its terms' counts, a sum free of one counts 1 (a scalar
/// coefficient), a product multiplies, `x^n` (n > 0) raises to n, negative or
/// symbolic powers and `denom(..)` count 1, any other function multiplies its
/// arguments' counts. Equal to the number of expanded monomials when no two of
/// them merge, an upper bound otherwise; saturates at `f64::MAX`.
pub fn distributed_lorentz_monomials(expression: AtomView<'_>) -> f64 {
    fn walk(view: AtomView<'_>) -> (f64, bool) {
        match view {
            AtomView::Num(_) | AtomView::Var(_) => (1.0, false),
            AtomView::Fun(f) => {
                let symbol = f.get_symbol();
                let name = symbol.get_name();
                if name == names::MINK {
                    return (1.0, true);
                }
                if name.ends_with("::denom") {
                    return (1.0, false);
                }
                f.iter()
                    .map(walk)
                    .fold((1.0, false), |(d, l), (x, y)| (d * x, l || y))
            }
            AtomView::Add(a) => {
                let children = a.iter().map(walk).collect::<Vec<_>>();
                if children.iter().any(|(_, l)| *l) {
                    (children.iter().map(|(d, _)| d).sum(), true)
                } else {
                    (1.0, false)
                }
            }
            AtomView::Mul(m) => m
                .iter()
                .map(walk)
                .fold((1.0, false), |(d, l), (x, y)| (d * x, l || y)),
            AtomView::Pow(p) => {
                let (base, exponent) = p.get_base_exp();
                match i64::try_from(exponent) {
                    Ok(n) if n > 0 => {
                        let (d, l) = walk(base);
                        (d.powi(n.min(i64::from(i32::MAX)) as i32), l)
                    }
                    _ => (1.0, false),
                }
            }
        }
    }
    walk(expression).0.min(f64::MAX)
}

/// The resource guard's memory bound for the contract-first order: `60 MB +
/// 350 bytes x D` (D = [`distributed_lorentz_monomials`]), the envelope of
/// the four measured R* hard Taylor polynomials of nf0 #32 (peaks 58 MB,
/// 166 MB, 10.2 GB, 20.4 GB at D = 3.4e4, 5.3e5, 3.2e7, 1.4e8; D5-HARD-KIND 8.1).
pub fn contract_first_bound_bytes(monomials: f64) -> f64 {
    60.0e6 + 350.0 * monomials
}

/// The vakint momentum `SUM c_i k(i) + SUM d_j p(j)` of every edge, with an
/// optional raw index on every vector; for an id that is not an edge of the
/// graph, its CONSTANT EXTERNAL VECTOR (reader rule 2): `p(X + 1 + position)`
/// after the graph's `X` external momenta, in the order of
/// [`FactorizedNumerator::constant_vectors`].
fn momentum(value: &FactorizedNumerator, id: i64, index: Option<i64>) -> Option<Atom> {
    let Some(edge) = value.edges.iter().find(|edge| edge.id as i64 == id) else {
        let position = value
            .constant_vectors()
            .iter()
            .position(|x| *x as i64 == id)?;
        let n = Atom::num((value.externals() + 1 + position) as i64);
        return Some(match index {
            Some(index) => function!(S.p, n, Atom::num(index)),
            None => function!(S.p, n),
        });
    };
    let loops = value.loops();
    let mut out = Atom::Zero;
    for (i, c) in edge.momentum.iter().enumerate() {
        if *c == 0 {
            continue;
        }
        let (head, n) = if i < loops {
            (S.k, i as i64 + 1)
        } else {
            (S.p, (i - loops) as i64 + 1)
        };
        let vector = match index {
            Some(index) => function!(head, Atom::num(n), Atom::num(index)),
            None => function!(head, Atom::num(n)),
        };
        out += Atom::num(*c) * vector;
    }
    Some(out)
}

/// `(edge id, its mink argument)` of `Q(e, mink(D, ..))` / `Q(e, mink(D))`.
fn momentum_parts(view: AtomView<'_>) -> Option<(i64, AtomView<'_>)> {
    let AtomView::Fun(function) = view else {
        return None;
    };
    if function.get_symbol().get_name() != names::MOMENTUM || function.get_nargs() != 2 {
        return None;
    }
    let mut args = function.iter();
    let edge = i64::try_from(args.next()?).ok()?;
    Some((edge, args.next()?))
}

fn is_compact(rep: AtomView<'_>) -> bool {
    matches!(rep, AtomView::Fun(f) if f.get_symbol().get_name() == names::MINK && f.get_nargs() == 1)
}

/// The checks of the vacuum mode (the contract's graph checks apply where they
/// do not assume a propagator graph).
fn validate_vacuum(value: &FactorizedNumerator) -> Result<()> {
    // a value with new features: its shape refusals are the caller's
    // mistakes (`Graph`, the same "vacuum: " text), not internal failures
    let refuse = |detail: &str| FactorizedError::Graph {
        detail: format!("vacuum: {detail}"),
    };
    if !value.has_new_features() {
        return validate_vacuum_v1(value);
    }
    if value.kind != FactorizedKind::Direct {
        return Err(FactorizedError::Reserved {
            feature: "the closure kind in the vacuum mode".into(),
        });
    }
    // MARKED SIDE: a marked vacuum value's own atoms get the reserved
    // namespace and wildcard walk of `validate` (the vacuum mode never calls
    // `validate`); an unmarked value is not walked (the base behaviour)
    if super::is_marked(value) {
        value.check_reserved_symbols()?;
    }
    // the bounds of `validate` (depth, ids, edge powers) and the momentum
    // coefficients, which `validate` checks in the contract's momentum rules
    value.check_marked_bounds(true)?;
    if super::is_marked(value) {
        value.check_marked_prefactor()?;
        value.check_paired_legs()?;
        value.check_vacuum_original_cycles()?;
        value.check_marked_atoms()?;
        // every caller field outside the insertions: defined, finite (the
        // entry repeats it with the settings' ε)
        value.check_marked_fields(None)?;
    }
    if !value.insertions.is_empty() {
        // The insertion rules (`Insertion`); the vacuum mode does not run the
        // whole contract (hard parts keep half-edges of a larger graph).
        let edges = value.edges.iter().map(|e| (e.id, e)).collect();
        value.check_insertions(&edges)?;
    }
    let loops = value.loops();
    if loops == 0 {
        return Err(refuse("a vacuum value has at least one loop"));
    }
    let mut mass = None::<&Atom>;
    for edge in &value.edges {
        if edge.momentum.len() != loops + value.externals() {
            return Err(refuse(&format!(
                "edge {}: momentum of the wrong length",
                edge.id
            )));
        }
        if edge.is_external() {
            continue;
        }
        if edge.momentum[loops..].iter().any(|c| *c != 0) {
            return Err(refuse(&format!(
                "edge {} carries an external momentum (not a vacuum integral)",
                edge.id
            )));
        }
        if !edge.mass_squared.is_zero() {
            match mass {
                None => mass = Some(&edge.mass_squared),
                Some(m) if *m == edge.mass_squared => {}
                Some(m) => {
                    return Err(refuse(&format!(
                        "two mass scales ({m} and {}): single-scale vacua only",
                        edge.mass_squared
                    )));
                }
            }
        }
    }
    if mass.is_none() {
        return Err(refuse(
            "no massive line: a massless vacuum integral is scaleless",
        ));
    }
    // The OPEN slots (ports no projector closes, fresh projector slots):
    // Lorentz only. Bispinor ports closed by the caller's gamma projector
    // (a Dirac line closed into a trace) are fine.
    for port in &value.open_slots()? {
        if !matches!(port.as_view(), AtomView::Fun(f) if f.get_symbol().get_name() == names::MINK && f.get_nargs() == 2)
        {
            return Err(refuse(&format!(
                "open slot {port}: only Lorentz slots mink(D, index) stay open in the vacuum scope"
            )));
        }
    }
    Ok(())
}

/// The closed-except-ports scalar of `projector · PROD factors` contracted by
/// idenso and mapped to vakint's vacuum notation (steps 1 and 2), in the
/// value's dimension symbol; and the seconds of the contraction.
pub fn contract_to_vakint_vacuum(value: &FactorizedNumerator) -> Result<(Atom, f64)> {
    validate_vacuum(value)?;
    let started = Instant::now();
    let contracted = value.contract_naively()?;
    let seconds = started.elapsed().as_secs_f64();
    Ok((map_to_vakint_vacuum(value, &contracted)?, seconds))
}

/// [`contract_to_vakint_vacuum`] with the contraction STREAMED in `pieces`
/// pieces ([`contract_streamed`]): the same scalar, and the seconds and the
/// number of pieces actually contracted.
pub fn contract_to_vakint_vacuum_streamed(
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    pieces: usize,
) -> Result<(Atom, f64, usize)> {
    validate_vacuum(value)?;
    let started = Instant::now();
    let mut network = value.factor_product();
    if let Some(projector) = &value.projector {
        network *= &projector.tensor;
    }
    let (contracted, contracted_pieces) = contract_streamed(settings, &network, pieces)?;
    let seconds = started.elapsed().as_secs_f64();
    Ok((
        map_to_vakint_vacuum(value, &contracted)?,
        seconds,
        contracted_pieces,
    ))
}

/// Step 2 of the module documentation: an idenso-contracted network in vakint's
/// vacuum notation.
pub(super) fn map_to_vakint_vacuum(value: &FactorizedNumerator, contracted: &Atom) -> Result<Atom> {
    if !value.has_new_features() {
        return map_to_vakint_vacuum_v1(value, contracted);
    }
    // the open slots: the ports when there is no projector (raw index PORT + i
    // by position, as always), what a projector leaves open otherwise
    let open = value.open_slots()?;
    let port_index = |slot: AtomView<'_>| {
        open.iter()
            .position(|port| port.as_view() == slot)
            .map(|i| PORT + i as i64)
    };
    let failure = std::cell::Cell::new(None::<String>);
    let mapped = contracted.replace_map(|part, _, out| {
        let AtomView::Fun(function) = part else {
            return;
        };
        let symbol = function.get_symbol();
        let name = symbol.get_name();
        // a scalar product of two edge momenta
        if (name == names::METRIC || name == "spenso::dot") && function.get_nargs() == 2 {
            let args = function.iter().collect::<Vec<_>>();
            if let (Some((a, ra)), Some((b, rb))) =
                (momentum_parts(args[0]), momentum_parts(args[1]))
                && is_compact(ra)
                && is_compact(rb)
            {
                match (momentum(value, a, None), momentum(value, b, None)) {
                    (Some(x), Some(y)) => **out = S.dot.call((&x, &y)),
                    _ => failure.set(Some(format!("a dot of unknown edges {a}, {b}"))),
                }
                return;
            }
            // a metric between two open ports
            if name == names::METRIC
                && let (Some(i), Some(j)) = (port_index(args[0]), port_index(args[1]))
            {
                **out = S.g.call((Atom::num(i), Atom::num(j)));
                return;
            }
            return;
        }
        // an edge momentum carrying an open port
        if let Some((edge, slot)) = momentum_parts(part)
            && let Some(index) = port_index(slot)
        {
            match momentum(value, edge, Some(index)) {
                Some(vector) => **out = vector,
                None => failure.set(Some(format!("unknown edge {edge}"))),
            }
        }
    });
    if let Some(detail) = failure.take() {
        return Err(error(detail));
    }
    let mut leftover = None;
    mapped.as_view().visitor(&mut |part| {
        if let AtomView::Fun(function) = part {
            let symbol = function.get_symbol();
            let name = symbol.get_name();
            if [
                names::MOMENTUM,
                "spenso::dot",
                names::MINK,
                names::BIS,
                names::METRIC,
                names::GAMMA,
            ]
            .contains(&name)
            {
                leftover = Some(part.to_owned());
            }
        }
        leftover.is_none()
    });
    if let Some(part) = leftover {
        return Err(error(format!(
            "outside the vacuum scope (metrics and momenta with open indices) after idenso: {part}"
        )));
    }
    Ok(mapped)
}

/// How the vacuum entry contracts its network (D5-HARD-KIND 9.4 items 1-2).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Streaming {
    /// One idenso contraction of the whole network (the unstreamed path).
    Off,
    /// The network split into (at most) this many pieces ([`split_network`]),
    /// contracted one after another, added into one sum ([`contract_streamed`]).
    Pieces(usize),
    /// [`stream_pieces_for`] of the size predictor: `Off` when it gives one
    /// piece, else that many pieces.
    Auto,
}

/// Lorentz monomials per piece of the streaming rule.
pub const STREAM_MONOMIALS_PER_PIECE: f64 = 1.0e5;
/// Largest piece count of the streaming rule.
pub const STREAM_MAX_PIECES: usize = 32;

/// THE STREAMING RULE: the piece count for a network that distributes to
/// `monomials` Lorentz monomials ([`distributed_lorentz_monomials`]): 1 up to
/// [`STREAM_MONOMIALS_PER_PIECE`], else the smallest power of two not below
/// `monomials / STREAM_MONOMIALS_PER_PIECE`, at most [`STREAM_MAX_PIECES`].
/// (Measured, D5-HARD-KIND 9.3: nf0 #32's hard terms at D = 3.4e4 / 5.3e5 /
/// 3.2e7 / 1.4e8 get 1 / 8 / 32 / 32 pieces.)
pub fn stream_pieces_for(monomials: f64) -> usize {
    if monomials.is_nan() || monomials <= STREAM_MONOMIALS_PER_PIECE {
        return 1;
    }
    let wanted = (monomials / STREAM_MONOMIALS_PER_PIECE).ceil();
    let mut pieces = 1;
    while (pieces as f64) < wanted && pieces < STREAM_MAX_PIECES {
        pieces *= 2;
    }
    pieces
}

/// One split: a sum into the halves of its terms; a product at its largest
/// sum factor (by byte size, the first of equals) into the halves of that
/// sum's terms, each times the other factors; anything else does not split.
fn split_once(piece: &Atom) -> Option<[Atom; 2]> {
    let halves = |terms: Vec<Atom>| {
        let half = terms.len() / 2;
        (half > 0).then(|| {
            [
                Atom::add_many(terms[..half].iter()),
                Atom::add_many(terms[half..].iter()),
            ]
        })
    };
    match piece.as_view() {
        AtomView::Add(sum) => halves(sum.iter().map(|t| t.to_owned()).collect()),
        AtomView::Mul(product) => {
            let factors = product.iter().collect::<Vec<_>>();
            let mut largest = None::<(usize, usize)>;
            for (i, factor) in factors.iter().enumerate() {
                if let AtomView::Add(_) = factor {
                    let size = factor.get_byte_size();
                    if largest.is_none_or(|(_, best)| size > best) {
                        largest = Some((i, size));
                    }
                }
            }
            let (at, _) = largest?;
            let AtomView::Add(sum) = factors[at] else {
                return None;
            };
            let [a, b] = halves(sum.iter().map(|t| t.to_owned()).collect())?;
            let mut rest = Atom::num(1);
            for (i, factor) in factors.iter().enumerate() {
                if i != at {
                    rest *= factor.to_owned();
                }
            }
            Some([a * &rest, b * &rest])
        }
        _ => None,
    }
}

/// The network as a list of pieces whose sum it is: the largest piece (by byte
/// size, the first of equals) is split ([`split_once`]) until there are
/// `pieces` pieces or none splits. Deterministic: Symbolica's normal form
/// fixes the order of every sum and product.
pub fn split_network(network: &Atom, pieces: usize) -> Vec<Atom> {
    let mut list = vec![network.clone()];
    while list.len() < pieces {
        let mut largest = 0;
        for (i, piece) in list.iter().enumerate() {
            if piece.as_view().get_byte_size() > list[largest].as_view().get_byte_size() {
                largest = i;
            }
        }
        let Some([a, b]) = split_once(&list[largest]) else {
            break;
        };
        list[largest] = a;
        list.insert(largest + 1, b);
    }
    list
}

/// The STREAMED contraction (D5-HARD-KIND 9.4 items 1-2): the network split
/// into `pieces` pieces ([`split_network`]), each contracted by idenso with the
/// open slots kept ([`super::contract_with_idenso`]), one after another, in
/// order, with a resource checkpoint before every piece (the caller's budget
/// can stop the job between pieces); the contracted pieces are added into ONE
/// expanded sum (like terms combine), so that ONE tensor reduction follows.
/// The sum does not depend on the piece count. Returns the sum and the number
/// of pieces contracted.
pub fn contract_streamed(
    settings: &VakintSettings,
    network: &Atom,
    pieces: usize,
) -> Result<(Atom, usize)> {
    let list = split_network(network, pieces.max(1));
    let count = list.len();
    let streamed = settings
        .with_resource_budget(|| {
            let mut sum = Atom::Zero;
            for (i, piece) in list.iter().enumerate() {
                if let Err(limit) = crate::execution::checkpoint() {
                    return Ok(Err(error(format!(
                        "streamed contraction, before piece {} of {count}: {limit}",
                        i + 1
                    ))));
                }
                match super::contract_with_idenso(piece, false) {
                    Ok(closed) => sum += closed,
                    Err(e) => return Ok(Err(e)),
                }
            }
            Ok(Ok(sum.expand()))
        })
        .map_err(|e| error(format!("streamed contraction: {e}")))?;
    Ok((streamed?, count))
}

/// One propagator line of the vacuum topology handed to vakint.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumLine {
    pub id: usize,
    pub source: usize,
    pub target: usize,
    /// Over the loop momenta only.
    pub momentum: Vec<i64>,
    pub mass_squared: Atom,
    pub power: i64,
}

/// The lines of the vacuum graph: its internal edges with every chain of
/// lines in SERIES merged (a vertex with exactly two internal line ends whose
/// lines carry the same momentum through it and the same mass is removed, the
/// two lines become one with the summed power). This is graph normalisation
/// (vakint's matching takes one line per distinct propagator), as the merged
/// lines of D5 section 1.3; external legs carry no momentum through a vacuum
/// and do not count.
pub fn vacuum_lines(value: &FactorizedNumerator) -> Vec<VacuumLine> {
    let loops = value.loops();
    let mut lines = value
        .edges
        .iter()
        .filter(|edge| !edge.is_external())
        .map(|edge| {
            let (Endpoint::Vertex(source), Endpoint::Vertex(target)) = (edge.source, edge.target)
            else {
                unreachable!("internal edges have two vertices")
            };
            VacuumLine {
                id: edge.id,
                source,
                target,
                momentum: edge.momentum[..loops].to_vec(),
                mass_squared: edge.mass_squared.clone(),
                power: edge.power,
            }
        })
        .collect::<Vec<_>>();
    loop {
        let mut merged = false;
        let mut vertices = lines
            .iter()
            .flat_map(|line| [line.source, line.target])
            .collect::<Vec<_>>();
        vertices.sort_unstable();
        vertices.dedup();
        for v in vertices {
            let ends = lines
                .iter()
                .enumerate()
                .filter(|(_, l)| l.source == v || l.target == v)
                .collect::<Vec<_>>();
            if ends.len() != 2 || ends.iter().any(|(_, l)| l.source == l.target) {
                continue;
            }
            let (ia, a) = (ends[0].0, ends[0].1.clone());
            let (ib, b) = (ends[1].0, ends[1].1.clone());
            if a.mass_squared != b.mass_squared {
                continue;
            }
            // a into v, b out of v (reversing a line flips its momentum)
            let (a_from, a_in) = if a.target == v {
                (a.source, a.momentum.clone())
            } else {
                (a.target, a.momentum.iter().map(|c| -c).collect())
            };
            let (b_to, b_out) = if b.source == v {
                (b.target, b.momentum.clone())
            } else {
                (b.source, b.momentum.iter().map(|c| -c).collect())
            };
            if a_in != b_out {
                continue;
            }
            let line = VacuumLine {
                id: a.id.min(b.id),
                source: a_from,
                target: b_to,
                momentum: a_in,
                mass_squared: a.mass_squared.clone(),
                power: a.power + b.power,
            };
            let (first, second) = (ia.min(ib), ia.max(ib));
            lines.remove(second);
            lines[first] = line;
            merged = true;
            break;
        }
        if !merged {
            break;
        }
    }
    lines.sort_by_key(|line| line.id);
    lines
}

/// The massive closure of a vacuum graph with ONE massive line among massless
/// ones ([`vacuum_backend`]): the massive (closing) line is cut and its
/// momentum becomes the head `p(1)` of the massless p-integral of the other
/// lines.
#[derive(Clone, Debug, PartialEq)]
pub struct VacuumCut {
    /// The closing line (merged: its power is the chain's).
    pub line: VacuumLine,
    /// Every vacuum loop momentum `k(i)` (i = 1..L) and its image over the cut
    /// loop momenta `k(1..L-1)` and the head `p(1)`: a unimodular change of
    /// basis in which the closing line's momentum is `p(1)`.
    pub routing: Vec<(Atom, Atom)>,
    /// The cut graph, the massless lines with their cut momenta: when a line
    /// still carries a cut loop momentum, `topo(PROD prop(id, edge(s,t),
    /// momentum, 0, power))` of ALL of them, the lines that carry the head
    /// alone included (the cut stays one connected two-point graph, as
    /// gammalooprs's `VacuumCut` writes it); when none does (one loop),
    /// `PROD dot(p(1),p(1))^-power` (the loop-free closure input).
    pub topology: Atom,
}

/// How [`evaluate_factorized_vacuum`] integrates a vacuum graph, derived from
/// its lines (module documentation, step 3).
#[derive(Clone, Debug, PartialEq)]
pub enum VacuumBackend {
    /// Every line massive with the one mass: vakint's single-scale families.
    Family,
    /// One massive line among massless ones: the massive closure.
    Closure(VacuumCut),
    /// One massive line that no massless line depends on: exactly zero.
    Scaleless,
}

/// The back end of the lines of a vacuum graph with `loops` loops (module
/// documentation, step 3): every line massive -> family; exactly one massive
/// line -> its cut ([`VacuumCut`]); several massive lines among massless ones
/// -> an error (neither a single-scale family nor one closure).
///
/// The routing: with `m` the closing line's momentum over `k(1..L)` and `j`
/// its last loop momentum with a nonzero coefficient `s = m_j = +-1`,
/// `k(j) -> s (p(1) - SUM_{i != j} m_i k(i))` and every other `k(i)` keeps
/// its order among `k(1..L-1)`; a massless line `SUM_i c_i k(i)` becomes
/// `SUM_{i != j} (c_i - c_j s m_i) k(i) + c_j s p(1)`. A line on `p(1)` alone
/// stays a propagator of the cut graph unless the cut has no loop left.
pub fn vacuum_backend(lines: &[VacuumLine], loops: usize) -> Result<VacuumBackend> {
    let massive = lines
        .iter()
        .filter(|line| !line.mass_squared.is_zero())
        .collect::<Vec<_>>();
    if massive.len() == lines.len() {
        return Ok(VacuumBackend::Family);
    }
    let [closing] = massive.as_slice() else {
        return Err(error(format!(
            "{} massive lines among {} lines (after series merging): neither a single-scale family (every line massive) nor one closing line",
            massive.len(),
            lines.len()
        )));
    };
    let m = &closing.momentum;
    let Some(j) = m.iter().rposition(|c| *c != 0) else {
        return Err(error(format!(
            "line {}: carries no loop momentum",
            closing.id
        )));
    };
    let sign = m[j];
    if sign.abs() != 1 {
        return Err(error(format!(
            "line {}: momentum {m:?} has no unit coefficient to solve for",
            closing.id
        )));
    }
    let position = |i: usize| if i < j { i } else { i - 1 };
    let cut_loop = |i: usize| S.k.call(position(i) as i64 + 1);
    let head = S.p.call(1);
    let routing = (0..loops)
        .map(|i| {
            let image = if i == j {
                let mut image = Atom::num(sign) * &head;
                for (l, c) in m.iter().enumerate() {
                    if l != j && *c != 0 {
                        image += Atom::num(-sign * c) * cut_loop(l);
                    }
                }
                image
            } else {
                cut_loop(i)
            };
            (S.k.call(i as i64 + 1), image)
        })
        .collect::<Vec<_>>();
    let mut propagators = Atom::num(1);
    let mut loop_lines = 0;
    let mut tree = Atom::num(1);
    let propagator = |line: &VacuumLine, momentum: Atom| {
        S.prop.call((
            line.id as i64,
            S.edge.call((line.source as i64, line.target as i64)),
            momentum,
            Atom::Zero,
            Atom::num(line.power),
        ))
    };
    let mut depends = false;
    for line in lines.iter().filter(|line| line.id != closing.id) {
        let c = &line.momentum;
        let mut cut = vec![0i64; loops - 1];
        for i in (0..loops).filter(|i| *i != j) {
            cut[position(i)] = c[i] - c[j] * sign * m[i];
        }
        let on_head = c[j] * sign;
        depends |= on_head != 0;
        if cut.iter().all(|x| *x == 0) {
            if on_head.abs() != 1 {
                return Err(error(format!(
                    "line {}: momentum {c:?} vanishes or is not +-p(1) after the cut",
                    line.id
                )));
            }
            tree *= S.dot.call((&head, &head)).pow(Atom::num(-line.power));
            propagators *= propagator(line, Atom::num(on_head) * &head);
            continue;
        }
        let mut momentum = Atom::Zero;
        for (a, x) in cut.iter().enumerate() {
            if *x != 0 {
                momentum += Atom::num(*x) * S.k.call(a as i64 + 1);
            }
        }
        if on_head != 0 {
            momentum += Atom::num(on_head) * &head;
        }
        propagators *= propagator(line, momentum);
        loop_lines += 1;
    }
    if !depends {
        return Ok(VacuumBackend::Scaleless);
    }
    let topology = if loop_lines > 0 {
        S.topo.call(propagators)
    } else {
        tree
    };
    Ok(VacuumBackend::Closure(VacuumCut {
        line: (*closing).clone(),
        routing,
        topology,
    }))
}

/// `(external tensor, coefficient)` of a reduced vacuum numerator: the
/// monomials in the purely external objects (`g(i, j)`, `p(j, i)`,
/// `dot(p(a), p(b))`) and their scalar coefficients in the loop momenta.
fn external_tensors(numerator: &Atom) -> Vec<(Atom, Atom)> {
    let mut external = std::collections::BTreeSet::new();
    numerator.as_view().visitor(&mut |part| {
        if let AtomView::Fun(function) = part {
            let symbol = function.get_symbol();
            let external_dot = symbol == S.dot
                && function
                    .iter()
                    .all(|v| matches!(v, AtomView::Fun(f) if f.get_symbol() == S.p));
            if symbol == S.g || (symbol == S.p && function.get_nargs() == 2) || external_dot {
                external.insert(part.to_owned());
                return false;
            }
        }
        true
    });
    let external = external.into_iter().collect::<Vec<_>>();
    if external.is_empty() {
        vec![(Atom::num(1), numerator.clone())]
    } else {
        numerator.coefficient_list::<i8>(&external)
    }
}

/// The integrand as vakint's ONE term, read off its top-level product: the
/// entry builds it as `numerator × topo(..)`, so the topology is one factor
/// and the numerator the product of the others (exactly the pair vakint's
/// `split_integrals` returns, without its `coefficient_list` over the whole
/// expanded integrand, which dominated the entry's single-core reduction:
/// D5-HARD-KIND section 12); a bare topology is the term with numerator 1.
/// `None` for any other shape.
fn vacuum_term(input: &Atom) -> Result<Option<VakintTerm>> {
    let is_topo = |view: AtomView<'_>| matches!(view, AtomView::Fun(f) if f.get_symbol() == S.topo);
    if is_topo(input.as_view()) {
        let numerator = Atom::num(1);
        let vectors = VakintTerm::identify_vectors_in_numerator(numerator.as_view())
            .map_err(|e| error(e.to_string()))?;
        return Ok(Some(VakintTerm {
            integral: input.clone(),
            numerator,
            vectors,
        }));
    }
    let AtomView::Mul(product) = input.as_view() else {
        return Ok(None);
    };
    let mentions_topo = |view: AtomView<'_>| {
        let mut found = false;
        view.visitor(&mut |part| {
            found |= is_topo(part);
            !found
        });
        found
    };
    let (topologies, factors): (Vec<_>, Vec<_>) = product.iter().partition(|f| is_topo(*f));
    let [integral] = topologies.as_slice() else {
        return Ok(None);
    };
    if factors.iter().any(|f| mentions_topo(*f)) {
        return Ok(None);
    }
    let numerator = factors
        .iter()
        .fold(Atom::num(1), |product, factor| product * factor.to_owned());
    let vectors = VakintTerm::identify_vectors_in_numerator(numerator.as_view())
        .map_err(|e| error(e.to_string()))?;
    Ok(Some(VakintTerm {
        integral: integral.to_owned(),
        numerator,
        vectors,
    }))
}

/// Vakint's vacuum tensor reduction (dot notation) of the integrand `input`:
/// its one term ([`vacuum_term`]) reduced directly, else `Vakint::tensor_reduce`.
fn reduce_vacuum(vakint: &Vakint, settings: &VakintSettings, input: &Atom) -> Result<Atom> {
    let mut reducing = settings.clone();
    reducing.use_dot_product_notation = true;
    match vacuum_term(input)? {
        Some(term) => {
            let mut expression = VakintExpression(vec![term]);
            expression
                .tensor_reduce(vakint, &reducing)
                .map_err(|e| error(e.to_string()))?;
            Ok(expression.into())
        }
        None => vakint
            .tensor_reduce(&reducing, input.as_view())
            .map_err(|e| error(e.to_string())),
    }
}

/// The closure back end (module documentation, step 3): vakint's vacuum
/// tensor reduction of `input` ([`reduce_vacuum`]), then every external
/// tensor's coefficient through the massive closure of `cut`, all closures of
/// one value in one batch.
fn evaluate_closure(
    vakint: &Vakint,
    settings: &VakintSettings,
    input: &Atom,
    cut: &VacuumCut,
) -> Result<Atom> {
    let (tensors, inputs) = closure_inputs(vakint, settings, input, cut)?;
    let direct = matches!(cut.topology.as_view(), AtomView::Fun(f) if f.get_symbol() == S.topo);
    let values = if direct {
        MassiveClosure::evaluate_direct_batch(vakint, settings, &inputs)
            .map_err(|e| error(e.to_string()))?
    } else {
        inputs
            .iter()
            .map(|(closure, input)| {
                closure
                    .evaluate(vakint, settings, input.ordinary().as_view())
                    .map_err(|e| error(e.to_string()))
            })
            .collect::<Result<Vec<_>>>()?
    };
    let mut value = Atom::Zero;
    for (tensor, coefficient) in tensors.iter().zip(values) {
        value += tensor * coefficient;
    }
    Ok(value)
}

/// The first half of [`evaluate_closure`]: vakint's vacuum tensor reduction of
/// `input` ([`reduce_vacuum`]) and, per external tensor, the tensor and the
/// massive closure of its scalar coefficient through `cut` (not evaluated).
pub(super) type TensorClosures = (Vec<Atom>, Vec<(MassiveClosure, DirectClosureInput)>);

/// See [`TensorClosures`]: the external tensors and their closure inputs.
pub(super) fn closure_inputs(
    vakint: &Vakint,
    settings: &VakintSettings,
    input: &Atom,
    cut: &VacuumCut,
) -> Result<TensorClosures> {
    let reduced = reduce_vacuum(vakint, settings, input)?;
    let expression = VakintExpression::try_from(reduced).map_err(|e| error(e.to_string()))?;
    let mut tensors = Vec::new();
    let mut coefficients = Vec::new();
    for term in &expression.0 {
        for (tensor, coefficient) in external_tensors(&term.numerator.expand()) {
            let mut leftover = None;
            coefficient.as_view().visitor(&mut |part| {
                if let AtomView::Fun(function) = part {
                    let symbol = function.get_symbol();
                    if symbol == S.p
                        || symbol == S.g
                        || (symbol == S.k && function.get_nargs() != 1)
                    {
                        leftover = Some(part.to_owned());
                    }
                }
                leftover.is_none()
            });
            if let Some(part) = leftover {
                return Err(error(format!(
                    "the vacuum tensor reduction left {part} in a scalar coefficient"
                )));
            }
            tensors.push(tensor);
            coefficients.push(coefficient);
        }
    }
    let closure = MassiveClosure {
        momentum: S.p.call(1),
        mass_squared: cut.line.mass_squared.clone(),
        power: cut.line.power,
    };
    let inputs = coefficients
        .into_iter()
        .map(|numerator| {
            (
                closure.clone(),
                DirectClosureInput {
                    topology: cut.topology.clone(),
                    routing: cut.routing.clone(),
                    numerator,
                },
            )
        })
        .collect::<Vec<_>>();
    Ok((tensors, inputs))
}

/// `PROD prop(id, edge(s,t), momentum, mass², power)` of the vacuum lines
/// (the argument of `topo(..)` of the integrand).
pub(super) fn vacuum_topology(lines: &[VacuumLine]) -> Atom {
    let mut topology = Atom::num(1);
    for line in lines {
        let mut momentum = Atom::Zero;
        for (i, c) in line.momentum.iter().enumerate() {
            if *c != 0 {
                momentum += Atom::num(*c) * S.k.call(i as i64 + 1);
            }
        }
        topology *= S.prop.call((
            line.id as i64,
            S.edge.call((line.source as i64, line.target as i64)),
            momentum,
            line.mass_squared.clone(),
            Atom::num(line.power),
        ));
    }
    topology
}

/// The integrand vakint receives (module documentation): `prefactor ×
/// scalar × topo(PROD prop(id, edge(s,t), momentum, mass², power))`, `D → 4 −
/// 2ε`; and the contraction's seconds.
pub fn vakint_vacuum_input(
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
) -> Result<(Atom, f64)> {
    vakint_vacuum_input_streamed(settings, value, prefactor, Streaming::Off)
        .map(|(input, seconds, _)| (input, seconds))
}

/// [`vakint_vacuum_input`] with the contraction streamed as `streaming` says
/// (`Auto` reads as `Off` here: the caller resolves it); also the number of
/// pieces contracted (0: unstreamed).
pub fn vakint_vacuum_input_streamed(
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
    streaming: Streaming,
) -> Result<(Atom, f64, usize)> {
    // a marked value called here directly: the entry check too; the
    // record itself is walked by validate_vacuum
    if !matches!(value.kind, FactorizedKind::Closure(_)) {
        if super::is_marked(value) {
            // the settings' rule before any computing rule
            super::check_master_modes(settings)?;
        }
        super::check_marked_entry(
            value,
            prefactor,
            vk_symbol!(&settings.epsilon_symbol),
            Some(settings.number_of_terms_in_epsilon_expansion),
        )?;
        super::check_entry_names(value, prefactor, settings)?;
    }
    value.check_insertion_epsilon(vk_symbol!(&settings.epsilon_symbol).get_name())?;
    let (scalar, seconds, pieces) = match streaming {
        Streaming::Pieces(n) => contract_to_vakint_vacuum_streamed(settings, value, n)?,
        Streaming::Off | Streaming::Auto => {
            let (scalar, seconds) = contract_to_vakint_vacuum(value)?;
            (scalar, seconds, 0)
        }
    };
    let topology = vacuum_topology(&vacuum_lines(value));
    let epsilon = Atom::var(vk_symbol!(&settings.epsilon_symbol));
    let in_epsilon = Atom::num(4) - Atom::num(2) * epsilon;
    let input = (prefactor * scalar * S.topo.call(topology))
        .replace(Atom::var(value.dimension))
        .with(in_epsilon);
    Ok((input, seconds, pieces))
}

/// The vakint-notation value mapped back to the caller's notation: raw port
/// indices to the port slots, `p(j)` to `Q(x_j, ..)` of the external leg of
/// `p(j)` or of the constant external vector `x_j` (an edge outside the
/// graph); an external vector without either stays `p(j)`.
pub fn to_caller_notation(value: &FactorizedNumerator, vakint_value: &Atom) -> Atom {
    if !value.has_new_features() {
        return to_caller_notation_v1(value, vakint_value);
    }
    let q = symbolica::symbol!(names::MOMENTUM);
    let open = value.open_slots().unwrap_or_else(|_| value.ports.clone());
    // the representation `mink(D)` of the caller's scalar products: from the
    // first open Lorentz slot, else from the first Lorentz port (a projector
    // that closes every port leaves no open slot; the ports still name D)
    let rep = open
        .iter()
        .chain(value.ports.iter())
        .find_map(|port| match port.as_view() {
            AtomView::Fun(f) if f.get_symbol().get_name() == names::MINK && f.get_nargs() == 2 => {
                Some(function!(
                    f.get_symbol(),
                    f.iter().next().unwrap().to_owned()
                ))
            }
            _ => None,
        });
    // a raw index, or a representation slot `mink(D', raw index)`
    let port = |index: AtomView<'_>| {
        let raw = match index {
            AtomView::Fun(f) if f.get_nargs() == 2 => {
                f.iter().nth(1).and_then(|i| i64::try_from(i).ok())
            }
            _ => i64::try_from(index).ok(),
        };
        raw.and_then(|i| usize::try_from(i - PORT).ok())
            .and_then(|i| open.get(i).cloned())
    };
    let constant = value.constant_vectors();
    let leg = |j: AtomView<'_>| {
        let j = i64::try_from(j)
            .ok()
            .and_then(|j| usize::try_from(j - 1).ok())?;
        match value.external_edges.get(j) {
            Some(leg) => *leg,
            None => constant.get(j - value.externals()).copied(),
        }
    };
    vakint_value.replace_map(|part, _, out| {
        let AtomView::Fun(function) = part else {
            return;
        };
        let args = function.iter().collect::<Vec<_>>();
        if function.get_symbol() == S.g && args.len() == 2 {
            if let (Some(a), Some(b)) = (port(args[0]), port(args[1])) {
                **out = function!(ETS.metric, a, b);
            }
        } else if function.get_symbol() == S.p && args.len() == 2 {
            if let (Some(x), Some(slot)) = (leg(args[0]), port(args[1])) {
                **out = function!(q, Atom::num(x as i64), slot);
            }
        } else if function.get_symbol() == S.dot && args.len() == 2 {
            let vector = |v: AtomView<'_>| match v {
                AtomView::Fun(f) if f.get_symbol() == S.p && f.get_nargs() == 1 => {
                    leg(f.iter().next().unwrap())
                }
                _ => None,
            };
            if let (Some(x), Some(y), Some(rep)) = (vector(args[0]), vector(args[1]), rep.as_ref())
            {
                **out = function!(
                    SPENSO_TAG.dot,
                    function!(q, Atom::num(x as i64), rep),
                    function!(q, Atom::num(y as i64), rep)
                );
            }
        }
    })
}

/// Evaluate a single-scale VACUUM value (module documentation) times
/// `prefactor` with vakint's existing vacuum evaluation (a family, or the
/// massive closure when one line is massive); open ports stay open.
pub fn evaluate_factorized_vacuum(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
) -> Result<VacuumValue> {
    evaluate_factorized_vacuum_streamed(vakint, settings, value, prefactor, Streaming::Auto)
}

/// [`evaluate_factorized_vacuum`] with the contraction streamed as `streaming`
/// says ([`Streaming`]; the default entry uses `Auto`). The value does not
/// depend on the choice: every piece count gives the same contracted sum.
pub fn evaluate_factorized_vacuum_streamed(
    vakint: &Vakint,
    settings: &VakintSettings,
    value: &FactorizedNumerator,
    prefactor: &Atom,
    streaming: Streaming,
) -> Result<VacuumValue> {
    // a marked value: the separately supplied prefactor gets the checks of
    // the record's prefactor fields, the settings' epsilon the caller-symbol
    // rule (no-ops for an unmarked value: the base behaviour; the closure
    // kind keeps its earlier refusal of this entry, `Reserved`)
    if !matches!(value.kind, super::FactorizedKind::Closure(_)) {
        if super::is_marked(value) {
            // the settings' rule before any computing rule
            super::check_master_modes(settings)?;
        }
        super::check_marked_entry(
            value,
            prefactor,
            vk_symbol!(&settings.epsilon_symbol),
            Some(settings.number_of_terms_in_epsilon_expansion),
        )?;
        super::check_entry_names(value, prefactor, settings)?;
        super::check_data_atoms_catalogue(value, vakint)?;
    }
    // the new record kinds do not depend on the caller's tensor reduction or
    // notation (evaluate::is_new_kind); values accepted before keep them
    let pinned;
    let settings = if super::evaluate::is_new_kind(value) {
        pinned = super::evaluate::pinned_settings(settings);
        &pinned
    } else {
        settings
    };
    // Lane unknown-masters: strict refusal once, on this call's summed value.
    if let Some(_boundary) = crate::public_result_boundary() {
        let inner = crate::unknown_master_coefficients_carried(settings);
        let result =
            evaluate_factorized_vacuum_streamed(vakint, &inner, value, prefactor, streaming)?;
        vakint
            .refuse_unknown_master_coefficients(settings, &result.vakint_value)
            .map_err(|error| FactorizedError::Evaluation {
                detail: error.to_string(),
            })?;
        return Ok(result);
    }
    // the resource guard (H2): refuse, naming the predicted size, rather than
    // run into the caller's process memory budget
    let mut network = value.factor_product();
    if let Some(projector) = &value.projector {
        network *= &projector.tensor;
    }
    let predicted_monomials = distributed_lorentz_monomials(network.as_view());
    if let Some(limit) = settings.rustred_resources.max_process_rss_bytes {
        let bound = contract_first_bound_bytes(predicted_monomials);
        if bound > limit as f64 {
            return Err(error(format!(
                "the uncontracted network distributes to D = {predicted_monomials:.3e} Lorentz monomials; its contraction is bounded by {bound:.3e} bytes, above the budget of {limit} bytes"
            )));
        }
    }
    let streaming = match streaming {
        Streaming::Auto => match stream_pieces_for(predicted_monomials) {
            1 => Streaming::Off,
            n => Streaming::Pieces(n),
        },
        other => other,
    };
    let (input, contract_seconds, pieces) =
        vakint_vacuum_input_streamed(settings, value, prefactor, streaming)?;
    let backend = vacuum_backend(&vacuum_lines(value), value.loops())?;
    let started = Instant::now();
    let vakint_value = match &backend {
        VacuumBackend::Family => vakint
            .evaluate(settings, input.as_view())
            .map_err(|e| error(e.to_string()))?,
        VacuumBackend::Closure(cut) => evaluate_closure(vakint, settings, &input, cut)?,
        VacuumBackend::Scaleless => Atom::Zero,
    };
    let evaluate_seconds = started.elapsed().as_secs_f64();
    let caller = to_caller_notation(value, &vakint_value);
    // (the vacuum caller notation keeps p(j) of an external vector without a
    // leg, and its scalar products, in vakint notation: allowed here)
    super::backstop::check_result(value, prefactor, settings, &caller, true)?;
    // the value in vakint's notation is returned too (`vakint_value`): the
    // same check, in vakint's notation
    super::backstop::check_result(value, prefactor, settings, &vakint_value, true)?;
    Ok(VacuumValue {
        value: caller,
        vakint_value,
        input,
        backend,
        predicted_monomials,
        pieces,
        contract_seconds,
        evaluate_seconds,
    })
}

// The base's vacuum rules and mappings, verbatim, for values without new
// features.

/// The base's `validate_vacuum`, verbatim (the rules of a value without new
/// features: [`FactorizedNumerator::has_new_features`]).
fn validate_vacuum_v1(value: &FactorizedNumerator) -> Result<()> {
    if value.kind != FactorizedKind::Direct {
        return Err(FactorizedError::Reserved {
            feature: "the closure kind in the vacuum mode".into(),
        });
    }
    if !value.insertions.is_empty() {
        return Err(FactorizedError::Reserved {
            feature: "insertions".into(),
        });
    }
    let loops = value.loops();
    if loops == 0 {
        return Err(error("a vacuum value has at least one loop"));
    }
    let mut mass = None::<&Atom>;
    for edge in &value.edges {
        if edge.momentum.len() != loops + value.externals() {
            return Err(error(format!(
                "edge {}: momentum of the wrong length",
                edge.id
            )));
        }
        if edge.is_external() {
            continue;
        }
        if edge.momentum[loops..].iter().any(|c| *c != 0) {
            return Err(error(format!(
                "edge {} carries an external momentum (not a vacuum integral)",
                edge.id
            )));
        }
        if !edge.mass_squared.is_zero() {
            match mass {
                None => mass = Some(&edge.mass_squared),
                Some(m) if *m == edge.mass_squared => {}
                Some(m) => {
                    return Err(error(format!(
                        "two mass scales ({m} and {}): single-scale vacua only",
                        edge.mass_squared
                    )));
                }
            }
        }
    }
    if mass.is_none() {
        return Err(error(
            "no massive line: a massless vacuum integral is scaleless",
        ));
    }
    for port in &value.ports {
        if !matches!(port.as_view(), AtomView::Fun(f) if f.get_symbol().get_name() == names::MINK && f.get_nargs() == 2)
        {
            return Err(error(format!(
                "port {port}: only Lorentz ports mink(D, index) are in the vacuum scope"
            )));
        }
    }
    Ok(())
}

/// The base's `map_to_vakint_vacuum`, verbatim (the rules of a value without new
/// features: [`FactorizedNumerator::has_new_features`]).
fn map_to_vakint_vacuum_v1(value: &FactorizedNumerator, contracted: &Atom) -> Result<Atom> {
    let port_index = |slot: AtomView<'_>| {
        value
            .ports
            .iter()
            .position(|port| port.as_view() == slot)
            .map(|i| PORT + i as i64)
    };
    let failure = std::cell::Cell::new(None::<String>);
    let mapped = contracted.replace_map(|part, _, out| {
        let AtomView::Fun(function) = part else {
            return;
        };
        let symbol = function.get_symbol();
        let name = symbol.get_name();
        // a scalar product of two edge momenta
        if (name == names::METRIC || name == "spenso::dot") && function.get_nargs() == 2 {
            let args = function.iter().collect::<Vec<_>>();
            if let (Some((a, ra)), Some((b, rb))) =
                (momentum_parts(args[0]), momentum_parts(args[1]))
                && is_compact(ra)
                && is_compact(rb)
            {
                match (momentum(value, a, None), momentum(value, b, None)) {
                    (Some(x), Some(y)) => **out = S.dot.call((&x, &y)),
                    _ => failure.set(Some(format!("a dot of unknown edges {a}, {b}"))),
                }
                return;
            }
            // a metric between two open ports
            if name == names::METRIC
                && let (Some(i), Some(j)) = (port_index(args[0]), port_index(args[1]))
            {
                **out = S.g.call((Atom::num(i), Atom::num(j)));
                return;
            }
            return;
        }
        // an edge momentum carrying an open port
        if let Some((edge, slot)) = momentum_parts(part)
            && let Some(index) = port_index(slot)
        {
            match momentum(value, edge, Some(index)) {
                Some(vector) => **out = vector,
                None => failure.set(Some(format!("unknown edge {edge}"))),
            }
        }
    });
    if let Some(detail) = failure.take() {
        return Err(error(detail));
    }
    let mut leftover = None;
    mapped.as_view().visitor(&mut |part| {
        if let AtomView::Fun(function) = part {
            let symbol = function.get_symbol();
            let name = symbol.get_name();
            if [
                names::MOMENTUM,
                "spenso::dot",
                names::MINK,
                names::BIS,
                names::METRIC,
                names::GAMMA,
            ]
            .contains(&name)
            {
                leftover = Some(part.to_owned());
            }
        }
        leftover.is_none()
    });
    if let Some(part) = leftover {
        return Err(error(format!(
            "outside the vacuum scope (metrics and momenta with open indices) after idenso: {part}"
        )));
    }
    Ok(mapped)
}

/// The base's `to_caller_notation`, verbatim (the rules of a value without new
/// features: [`FactorizedNumerator::has_new_features`]).
fn to_caller_notation_v1(value: &FactorizedNumerator, vakint_value: &Atom) -> Atom {
    let q = symbolica::symbol!(names::MOMENTUM);
    let rep = value.ports.first().and_then(|port| match port.as_view() {
        AtomView::Fun(f) => Some(function!(
            f.get_symbol(),
            f.iter().next().unwrap().to_owned()
        )),
        _ => None,
    });
    // a raw index, or a representation slot `mink(D', raw index)`
    let port = |index: AtomView<'_>| {
        let raw = match index {
            AtomView::Fun(f) if f.get_nargs() == 2 => {
                f.iter().nth(1).and_then(|i| i64::try_from(i).ok())
            }
            _ => i64::try_from(index).ok(),
        };
        raw.and_then(|i| usize::try_from(i - PORT).ok())
            .and_then(|i| value.ports.get(i).cloned())
    };
    let constant = value.constant_vectors();
    let leg = |j: AtomView<'_>| {
        let j = i64::try_from(j)
            .ok()
            .and_then(|j| usize::try_from(j - 1).ok())?;
        match value.external_edges.get(j) {
            Some(leg) => *leg,
            None => constant.get(j - value.externals()).copied(),
        }
    };
    vakint_value.replace_map(|part, _, out| {
        let AtomView::Fun(function) = part else {
            return;
        };
        let args = function.iter().collect::<Vec<_>>();
        if function.get_symbol() == S.g && args.len() == 2 {
            if let (Some(a), Some(b)) = (port(args[0]), port(args[1])) {
                **out = function!(ETS.metric, a, b);
            }
        } else if function.get_symbol() == S.p && args.len() == 2 {
            if let (Some(x), Some(slot)) = (leg(args[0]), port(args[1])) {
                **out = function!(q, Atom::num(x as i64), slot);
            }
        } else if function.get_symbol() == S.dot && args.len() == 2 {
            let vector = |v: AtomView<'_>| match v {
                AtomView::Fun(f) if f.get_symbol() == S.p && f.get_nargs() == 1 => {
                    leg(f.iter().next().unwrap())
                }
                _ => None,
            };
            if let (Some(x), Some(y), Some(rep)) = (vector(args[0]), vector(args[1]), rep.as_ref())
            {
                **out = function!(
                    SPENSO_TAG.dot,
                    function!(q, Atom::num(x as i64), rep),
                    function!(q, Atom::num(y as i64), rep)
                );
            }
        }
    })
}

#[cfg(test)]
mod tests;
