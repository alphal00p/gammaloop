//! Direct hand-over of an idenso-contracted scalar to the p lowering (helper VH,
//! numerator assignment, lane idenso-audit): shared lowering code for the
//! reference chain's direct variant (`refpath::direct`) and the input-type
//! adapter (`FactorizedNumerator` -> idenso -> lowering).
//!
//! The ordinary path re-enters vakint with the contracted scalar written in
//! the routed loop momenta and a `topo(..)`: the front end re-derives the
//! topology, matches it, canonicalises the (expanded) numerator, splits it by
//! terms and only then lowers it. Here the graph is matched ONCE from its
//! lines (the two steps of the front end: match and canonicalise the user's
//! topology, then match the canonical form, so the keys are the ordinary
//! path's), every edge momentum gets its image in the canonical momenta, and
//! the scalar, a polynomial in the scalar products of EDGE momenta, is read by
//! [`super::p_lowering::lower_closure`] without expanding it (products of sums
//! are multiplied as sparse polynomials; zero-sector pruning as configured).
//! The folded channels then take the ordinary reduction and finalisation.
//!
//! Exactness: the folded channels must equal the ordinary path's key by key
//! (`p_gate` views); every input this module does not reproduce exactly
//! declines and the caller keeps the ordinary path.
use std::collections::{BTreeMap, BTreeSet};

use symbolica::atom::{Atom, AtomCore, AtomView};
use symbolica::function;

use super::p_data::PLibraryHandle;
use super::p_gate::{self, GateView, ZeroSectorCount};
use super::p_integral::{Channel, FoldedInputs};
use super::p_lowering::{self, DirectStats, MatchedCut};
use crate::symbols::S;
use crate::utils::vakint_macros::vk_symbol;
use crate::{
    EvaluationMethod, ReplacementRules, Topology, Vakint, VakintError, VakintSettings, VakintTerm,
};

/// First vacuum index of a renamed edge momentum `k(EDGE_OFFSET + e)`.
pub(crate) const EDGE_OFFSET: i64 = 1_000_000;

/// One internal line of a two-point graph: its id, its endpoints, its momentum
/// as integer coefficients over the loop momenta `k(1..L)` and the external
/// head `p(1)` (which flows into the graph through its first leg), its power.
#[derive(Clone, Debug)]
pub struct PLine {
    pub id: i64,
    pub source: i64,
    pub target: i64,
    pub loops: Vec<i64>,
    pub head: i64,
    pub power: i64,
}

/// An edge whose momentum may occur in the scalar (lines and legs), with its
/// momentum in the same coordinates as [`PLine`].
#[derive(Clone, Debug)]
pub struct PEdge {
    pub id: i64,
    pub loops: Vec<i64>,
    pub head: i64,
}

/// A vector of the contracted scalar as the caller identifies it.
#[derive(Clone, Copy, Debug)]
pub enum VectorRef {
    /// The momentum of edge `id`.
    Edge(i64),
    /// `c · p(1)`.
    Head(i64),
}

#[derive(Clone, Debug)]
enum EdgeImage {
    Loop(Atom),
    Head(i64),
}

/// A graph matched once: everything the lowering needs that does not depend
/// on the numerator.
pub struct HandOver {
    cut: MatchedCut,
    routing: Vec<(Atom, Atom)>,
    images: BTreeMap<i64, EdgeImage>,
    names: BTreeMap<(usize, Atom), BTreeSet<String>>,
}

/// The outcome of [`match_lines`].
pub enum Matched {
    /// Scaleless after normalisation or match: an exact zero.
    Scaleless,
    /// Not reproduced exactly by this path (the ordinary path runs).
    Declined(String),
    Ready(Box<HandOver>),
}

/// Match a two-point graph given by its lines as the ordinary front end does
/// (normalise, match the user's topology, canonicalise, normalise again, match
/// the canonical form) and map every edge of `edges` into the canonical
/// momenta. `settings` must have been validated by the caller.
pub fn match_lines(
    vakint: &Vakint,
    settings: &VakintSettings,
    loops: usize,
    lines: &[PLine],
    edges: &[PEdge],
) -> Result<Matched, VakintError> {
    let head = function!(S.p, 1);
    let momentum_of = |l: &[i64], c: i64| -> Atom {
        let mut momentum = Atom::Zero;
        for (i, c) in l.iter().enumerate() {
            if *c != 0 {
                momentum += Atom::num(*c) * function!(S.k, i as i64 + 1);
            }
        }
        if c != 0 {
            momentum += Atom::num(c) * &head;
        }
        momentum
    };
    let mut topology = Atom::num(1);
    for line in lines {
        if line.loops.len() != loops {
            return Ok(Matched::Declined("a line with a wrong loop count".into()));
        }
        topology *= function!(
            S.prop,
            line.id,
            function!(S.edge, line.source, line.target),
            momentum_of(&line.loops, line.head),
            Atom::Zero,
            Atom::num(line.power)
        );
    }
    let topology = function!(S.topo, topology);
    // The two matches of the front end (lane integrate's fix ba2cfea8f of the
    // parked engine): a placeholder numerator k(1)·k(1) keeps the scalar-only
    // master transport off, as it is off for a loop-dependent numerator.
    let mut term = VakintTerm {
        integral: topology,
        numerator: Atom::num(1),
        vectors: Vec::new(),
    };
    if !term.normalize_propagator()? {
        return Ok(Matched::Scaleless);
    }
    if term.numerator != Atom::num(1) {
        return Ok(Matched::Declined(
            "a normalisation factor of the graph".into(),
        ));
    }
    let k1 = function!(S.k, 1);
    term.numerator = function!(S.dot, &k1, &k1);
    let Some(front) = vakint.match_topologies_for_settings(settings, term.integral.as_view())?
    else {
        return Ok(Matched::Declined(
            "the graph's topology is not recognised".into(),
        ));
    };
    term.canonicalize_expanding(settings, &front, false, false)?;
    if term.numerator.is_zero() {
        return Ok(Matched::Scaleless);
    }
    let mut again = VakintTerm {
        integral: term.integral.clone(),
        numerator: Atom::num(1),
        vectors: Vec::new(),
    };
    if !again.normalize_propagator()? {
        return Ok(Matched::Scaleless);
    }
    if again.numerator != Atom::num(1) {
        return Ok(Matched::Declined(
            "a normalisation factor of the canonical form".into(),
        ));
    }
    let rules: ReplacementRules = again.prepare_integral_specs(vakint, settings, None)?;
    if !matches!(rules.canonical_topology, Topology::Propagator(_)) {
        return Ok(Matched::Declined("not a p-integral topology".into()));
    }
    match vakint.evaluation_method(settings, &rules.canonical_topology) {
        Some(EvaluationMethod::RustRed(_))
            if crate::p_integral::massive::rustred_max_rank().is_none() => {}
        _ => return Ok(Matched::Declined("not reduced by RustRed alone".into())),
    }
    let one = Atom::num(1);
    let Some(request) = rules.propagator_request(one.as_view()) else {
        return Ok(Matched::Declined(
            "a p-integral match without a propagator request".into(),
        ));
    };
    if request.family.loops != loops {
        return Ok(Matched::Declined("the match changed the loop count".into()));
    }
    let (cut_head, sigma) = crate::p_integral::external_head(&request.external)?;
    let scaleless = crate::p_integral::parts::split(
        request.family,
        request.powers,
        Atom::num(1).as_view(),
        &request.external,
    )?
    .is_none();
    if scaleless {
        return Ok(Matched::Scaleless);
    }
    let mut names = BTreeMap::<(usize, Atom), BTreeSet<String>>::new();
    names
        .entry((request.family.loops, cut_head.clone()))
        .or_default()
        .insert(format!(
            "{}({})",
            request.family.name,
            request
                .powers
                .iter()
                .map(i64::to_string)
                .collect::<Vec<_>>()
                .join(",")
        ));
    let cut = MatchedCut {
        request_family: request.family,
        powers: request.powers.to_vec(),
        external: request.external.clone(),
        head: cut_head,
        sigma,
        factor: Atom::num(1),
        rules: rules.clone(),
        scaleless: false,
    };
    // Images of k(i), p(1) through the FIRST match (the canonical form's
    // momenta); lower_closure maps those through the second match itself.
    let front_image = |momentum: Atom| -> Result<Atom, VakintError> {
        let mut probe = VakintTerm {
            integral: Atom::num(1),
            numerator: momentum,
            vectors: Vec::new(),
        };
        probe.apply_numerator_replacement_rules(&front, settings)?;
        Ok(probe.numerator.expand())
    };
    let k_images = (1..=loops as i64)
        .map(|i| front_image(function!(S.k, i)))
        .collect::<Result<Vec<_>, _>>()?;
    let p_image = front_image(head.clone())?;
    let mut routing = Vec::new();
    let mut images = BTreeMap::new();
    for edge in edges {
        if edge.loops.len() != loops {
            return Ok(Matched::Declined("an edge with a wrong loop count".into()));
        }
        if edge.loops.iter().all(|c| *c == 0) {
            images.insert(edge.id, EdgeImage::Head(edge.head));
            continue;
        }
        let mut image = Atom::Zero;
        for (c, k) in edge.loops.iter().zip(&k_images) {
            if *c != 0 {
                image += Atom::num(*c) * k;
            }
        }
        if edge.head != 0 {
            image += Atom::num(edge.head) * &p_image;
        }
        let renamed = function!(S.k, EDGE_OFFSET + edge.id);
        routing.push((renamed.clone(), image.expand()));
        images.insert(edge.id, EdgeImage::Loop(renamed));
    }
    routing.push((head, p_image));
    Ok(Matched::Ready(Box::new(HandOver {
        cut,
        routing,
        images,
        names,
    })))
}

impl HandOver {
    fn image(&self, vector: VectorRef) -> Result<EdgeImage, String> {
        match vector {
            VectorRef::Head(c) => Ok(EdgeImage::Head(c)),
            VectorRef::Edge(edge) => self
                .images
                .get(&edge)
                .cloned()
                .ok_or_else(|| format!("an unknown edge {edge}")),
        }
    }

    /// The contracted scalar in the reader's momenta: every scalar product
    /// `dot_symbol(a, b)` whose arguments `vector` identifies becomes
    /// `vakint::dot` of the renamed edge momenta (multiples of `p(1)` for
    /// loop-free edges). `Err`: an argument `vector` does not identify, or a
    /// symbol of `forbidden` is left (unconverted momentum notation).
    pub fn rename(
        &self,
        scalar: &Atom,
        dot_symbol: symbolica::atom::Symbol,
        vector: impl Fn(AtomView<'_>) -> Option<VectorRef>,
        forbidden: &[symbolica::atom::Symbol],
    ) -> Result<Atom, String> {
        let head = function!(S.p, 1);
        let mut failure = None;
        let renamed = scalar.replace_map(|view, _, out| {
            let AtomView::Fun(f) = view else {
                return;
            };
            if f.get_symbol() != dot_symbol || f.get_nargs() != 2 {
                return;
            }
            let mut arguments = f.iter();
            let (a, b) = (arguments.next().unwrap(), arguments.next().unwrap());
            let image = |v: AtomView<'_>| {
                vector(v)
                    .ok_or_else(|| format!("an unknown vector {v}"))
                    .and_then(|r| self.image(r))
            };
            match (image(a), image(b)) {
                (Ok(EdgeImage::Loop(x)), Ok(EdgeImage::Loop(y))) => **out = function!(S.dot, x, y),
                (Ok(EdgeImage::Loop(x)), Ok(EdgeImage::Head(c)))
                | (Ok(EdgeImage::Head(c)), Ok(EdgeImage::Loop(x))) => {
                    **out = Atom::num(c) * function!(S.dot, x, &head)
                }
                (Ok(EdgeImage::Head(c)), Ok(EdgeImage::Head(e))) => {
                    **out = Atom::num(c * e) * function!(S.dot, &head, &head)
                }
                (Err(error), _) | (_, Err(error)) => failure = Some(error),
            }
        });
        if let Some(error) = failure {
            return Err(error);
        }
        if forbidden
            .iter()
            .chain([&dot_symbol])
            .any(|symbol| renamed.contains_symbol(*symbol))
        {
            return Err("momentum notation left after the rename".into());
        }
        // As vakint's FeynKit ingress: no tensor notation may survive (a
        // leftover would otherwise be read as a loop-free spectator).
        if let Some(symbol) = renamed
            .get_all_symbols(true)
            .into_iter()
            .find(|symbol| symbol.get_name().starts_with("spenso::"))
        {
            return Err(format!(
                "unconverted tensor notation {} in the contracted scalar",
                symbol.get_name()
            ));
        }
        Ok(renamed)
    }

    /// The folded channels of a renamed scalar ([`Self::rename`]; ε, not `D`,
    /// as vakint's input): the direct reader, slot substitution and fold of
    /// [`p_lowering::lower_closure`]. `Err(Declined)`: the ordinary path runs.
    pub(crate) fn lower(
        &self,
        settings: &VakintSettings,
        handle: &PLibraryHandle,
        scalar: AtomView<'_>,
        stats: &mut DirectStats,
    ) -> Result<Result<BTreeMap<Channel, FoldedInputs>, p_lowering::Declined>, VakintError> {
        p_lowering::lower_closure(settings, handle, &self.cut, &self.routing, scalar, stats)
    }

    /// The ordinary reduction and finalisation of folded channels (as the
    /// ordinary path's `reduce_requests_folded` + `finalize_reduced`).
    pub(crate) fn reduce(
        &self,
        vakint: &Vakint,
        settings: &VakintSettings,
        handle: &PLibraryHandle,
        folded: BTreeMap<Channel, FoldedInputs>,
    ) -> Result<Atom, VakintError> {
        use super::p_integral::{FoldMode, finalize_reduced, fold_mode, reduce_folded_channels};
        let mode = settings
            .evaluation_order
            .0
            .iter()
            .find_map(|method| match method {
                EvaluationMethod::RustRed(options) => Some(options.master_mode),
                _ => None,
            })
            .unwrap_or(crate::RustRedMasterMode::Symbolic);
        let fold = match fold_mode() {
            FoldMode::Off => FoldMode::Coefficient,
            other => other,
        };
        let reduced = reduce_folded_channels(vakint, settings, handle, folded, None, fold, false)?;
        finalize_reduced(vakint, settings, reduced, &self.names, mode)
    }
}

/// The gate view of folded channels without their zero-sector keys (the form
/// of the ordinary path's `GateObservation::view_without_zero_sectors`).
pub(crate) fn gate_view(
    settings: &VakintSettings,
    handle: &PLibraryHandle,
    folded: &BTreeMap<Channel, FoldedInputs>,
) -> (GateView, ZeroSectorCount) {
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    p_gate::drop_zero_sectors(handle, p_gate::gate_view(folded, epsilon))
}
