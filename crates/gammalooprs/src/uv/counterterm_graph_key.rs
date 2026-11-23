//! `GAMMALOOP_DELTA_KEY=graph`: the graph certificate of an R* child Δ and
//! its exact transport (delta-store design §1.3, §2; step 1, in process).
//!
//! The certificate identifies two children that differ only by
//! - a graph isomorphism (vertices, edges, half-edges, dummy indices, legs);
//! - the parent's loop-momentum routing: every `Q(e)` of a local rule is
//!   rewritten as the momentum `P(h)` flowing into the vertex of the half-edge
//!   `h` of `e` that the rule's vertex owns (`Q(e) = s(h) P(h)`, `s = -1` on
//!   the source half-edge), or `(P(t) - P(s))/2` in `e`'s own rule;
//! - the reversal of a self-conjugate (boson) edge: no flow, no orientation;
//! - the sign of each local rule: its monomials hang under two sign classes
//!   of one rule node (`σ`);
//! - the direction of a quark line whose rules are consistently oriented
//!   gamma chains: the two spinor slots of each of its gammas are unordered.
//!   A closed loop needs no transport (`tr Xᵀ = tr X`); an open line whose
//!   occurrence is crossed against the canonical orientation transposes the
//!   gammas of the stored value's chain between its two ports (`φ`).
//!
//! The stored value is `V = σ Φ(canon(Δ))`; a reader gets
//! `Δ = σ restore(Φ(V))`. A value the transposition cannot handle is never
//! served through a crossed line ([`TransportUnsupported`]).

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use super::*;
use crate::uv::delta_key;
use idenso::{dirac::AGS, representations::Bispinor};
use spenso::network::library::symbolic::ETS;
use symbolica::coefficient::CoefficientView;

static HEDGE_MOMENTUM: LazyLock<Symbol> =
    LazyLock::new(|| symbol!("uv_delta_cache_hedge_momentum"));

/// Why a lookup keeps the legacy certificate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Fallback {
    /// A rule carries the momentum of an edge it does not touch.
    NonIncident,
    /// A rule coefficient outside Q[i].
    Coefficient,
}

/// A crossed open quark line whose value carries a spinor structure other
/// than gamma chains and the spinor metric between its two ports.
#[derive(Debug)]
pub(crate) struct TransportUnsupported;

impl std::fmt::Display for TransportUnsupported {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "DELTA_TRANSPORT_UNSUPPORTED: the open-line transposition does not apply to this value"
        )
    }
}

impl std::error::Error for TransportUnsupported {}

pub(crate) fn is_transport_unsupported(error: &eyre::Report) -> bool {
    error
        .chain()
        .any(|cause| cause.downcast_ref::<TransportUnsupported>().is_some())
}

/// The symmetry data of one occurrence against the canonical frame; a hit
/// compares its own with the filler's (design §3.4).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub(crate) struct OccurrenceFeatures {
    pub(crate) sigma: i64,
    pub(crate) open: Vec<bool>,
    pub(crate) closed: Vec<bool>,
    pub(crate) crown: Vec<i64>,
    pub(crate) relabel: u64,
    pub(crate) legs: Vec<usize>,
    pub(crate) boson: Vec<bool>,
    pub(crate) routing: Vec<usize>,
    pub(crate) legacy: Option<u64>,
}

impl OccurrenceFeatures {
    /// The features of [`delta_key::FEATURES`] this hit shows against `fill`.
    pub(crate) fn shown(&self, fill: &Self) -> [bool; delta_key::FEATURES.len()] {
        [
            self.routing != fill.routing,
            self.relabel != fill.relabel,
            self.legs != fill.legs,
            self.boson != fill.boson,
            self.closed != fill.closed,
            self.open != fill.open,
            self.sigma * fill.sigma == -1,
            self.crown != fill.crown,
            matches!((self.legacy, fill.legacy), (Some(a), Some(b)) if a != b),
        ]
    }

    pub(crate) fn summary(&self) -> String {
        let bits = |bits: &[bool]| {
            bits.iter()
                .map(|bit| if *bit { '1' } else { '0' })
                .collect::<String>()
        };
        format!(
            "sigma={} open={} closed={} crown={} boson={} relabel={:016x} legs={:?} routing={:?}",
            self.sigma,
            bits(&self.open),
            bits(&self.closed),
            self.crown
                .iter()
                .map(|s| if *s > 0 { '+' } else { '-' })
                .collect::<String>(),
            bits(&self.boson),
            self.relabel,
            self.legs,
            self.routing
        )
    }
}

/// The transport of one occurrence.
#[derive(Clone, Debug)]
pub(crate) struct GraphTransport {
    pub(super) sigma: i64,
    /// The canonical port labels of each open line this occurrence crosses.
    crossed_open: Vec<[usize; 2]>,
    pub(crate) features: OccurrenceFeatures,
    /// Single-attached crown edges with the sign of their half-edge in S,
    /// when every half-edge at S is in S or its crown (momentum
    /// conservation `Σ s Q(e) = 0` holds): verify only.
    conservation: Option<Vec<(usize, i64)>>,
}

impl GraphTransport {
    pub(super) fn to_canonical(&self, canonical: Atom) -> Result<Atom> {
        let value = self.phi(&canonical)?;
        Ok(if self.sigma == 1 { value } else { value * -1 })
    }

    pub(super) fn from_canonical(&self, value: &Atom) -> Result<Atom> {
        self.phi(value)
    }

    fn phi(&self, value: &Atom) -> Result<Atom> {
        if self.crossed_open.is_empty() || value.is_zero() {
            return Ok(value.clone());
        }
        transpose_open_lines(value, &self.crossed_open)
            .ok_or_else(|| eyre::Report::new(TransportUnsupported))
    }
}

fn bis_rep() -> LibraryRep {
    Bispinor {}.into()
}

fn bis_head() -> Symbol {
    static HEAD: LazyLock<Symbol> = LazyLock::new(|| {
        match (Bispinor {}).new_rep(4).to_symbolic([Atom::Zero]).as_view() {
            AtomView::Fun(fun) => fun.get_symbol(),
            other => panic!("unexpected bispinor slot form {other}"),
        }
    });
    *HEAD
}

/// A raw (uncooked) bispinor slot `bis(dim, index)`: its index.
fn raw_bis_index(view: AtomView<'_>) -> Option<Atom> {
    match view {
        AtomView::Fun(fun) if fun.get_symbol() == bis_head() && fun.get_nargs() >= 2 => {
            fun.iter().last().map(|index| index.to_owned())
        }
        _ => None,
    }
}

fn raw_has_bis(view: AtomView<'_>) -> bool {
    match view {
        AtomView::Fun(fun) => {
            fun.get_symbol() == bis_head() || fun.iter().any(raw_has_bis)
        }
        AtomView::Add(add) => add.iter().any(raw_has_bis),
        AtomView::Mul(mul) => mul.iter().any(raw_has_bis),
        AtomView::Pow(pow) => {
            let (base, exponent) = pow.get_base_exp();
            raw_has_bis(base) || raw_has_bis(exponent)
        }
        _ => false,
    }
}

fn terms_of(value: &Atom) -> Vec<Atom> {
    match value.as_view() {
        AtomView::Add(add) => add.iter().map(|term| term.to_owned()).collect(),
        _ if value.is_zero() => Vec::new(),
        view => vec![view.to_owned()],
    }
}

fn factors_of(term: AtomView<'_>) -> Vec<Atom> {
    match term {
        AtomView::Mul(mul) => mul.iter().map(|factor| factor.to_owned()).collect(),
        view => vec![view.to_owned()],
    }
}

/// `φ` of every line in `lines` on the canonical value: the gammas of the
/// spinor chain between the line's two ports have their spinor slots swapped.
/// `None` when the value has another spinor structure on such a chain.
pub(crate) fn transpose_open_lines(value: &Atom, lines: &[[usize; 2]]) -> Option<Atom> {
    let expanded = value.expand();
    let mut out = Atom::Zero;
    for term in terms_of(&expanded) {
        let mut factors = factors_of(term.as_view());
        for ports in lines {
            transpose_line(&mut factors, ports)?;
        }
        let mut product = Atom::num(1);
        for factor in factors {
            product = product * factor;
        }
        out = out + product;
    }
    Some(out)
}

fn transpose_line(factors: &mut [Atom], ports: &[usize; 2]) -> Option<()> {
    let gamma = AGS.gamma;
    let metric = ETS.metric;
    let mut occurrences: HashMap<Atom, Vec<(usize, usize)>> = HashMap::new();
    let mut is_gamma = vec![false; factors.len()];
    for (position, factor) in factors.iter().enumerate() {
        if !raw_has_bis(factor.as_view()) {
            continue;
        }
        let AtomView::Fun(fun) = factor.as_view() else {
            return None;
        };
        let args = fun.iter().collect::<Vec<_>>();
        let spinors = if fun.get_symbol() == gamma
            && args.len() == 3
            && !raw_has_bis(args[2])
        {
            is_gamma[position] = true;
            [raw_bis_index(args[0])?, raw_bis_index(args[1])?]
        } else if fun.get_symbol() == metric && args.len() == 2 {
            [raw_bis_index(args[0])?, raw_bis_index(args[1])?]
        } else {
            return None;
        };
        for (slot, index) in spinors.into_iter().enumerate() {
            occurrences.entry(index).or_default().push((position, slot));
        }
    }
    let start = function!(*DELTA_PORT, ports[0] as i64);
    let end = function!(*DELTA_PORT, ports[1] as i64);
    if occurrences.get(&start)?.len() != 1 || occurrences.get(&end)?.len() != 1 {
        return None;
    }
    let mut current = start;
    let mut previous: Option<(usize, usize)> = None;
    let mut swap = BTreeSet::new();
    let mut visited = HashSet::new();
    loop {
        let list = occurrences.get(&current)?;
        let (position, slot) = match previous {
            None => list[0],
            Some(previous) => {
                if list.len() != 2 {
                    return None;
                }
                *list.iter().find(|occurrence| **occurrence != previous)?
            }
        };
        if !visited.insert(position) {
            return None;
        }
        if is_gamma[position] {
            swap.insert(position);
        }
        let AtomView::Fun(fun) = factors[position].as_view() else {
            return None;
        };
        let other = raw_bis_index(fun.iter().nth(1 - slot)?)?;
        previous = Some((position, 1 - slot));
        if other == end {
            break;
        }
        current = other;
    }
    for position in swap {
        let AtomView::Fun(fun) = factors[position].as_view() else {
            return None;
        };
        let args = fun.iter().collect::<Vec<_>>();
        factors[position] = FunctionBuilder::new(fun.get_symbol())
            .add_arg(args[1])
            .add_arg(args[0])
            .add_arg(args[2])
            .finish();
    }
    Some(())
}

/// One edge of `S ∪ crown`.
struct EdgeFacts {
    boundary: bool,
    self_conjugate: bool,
    quark: bool,
    particle: Option<usize>,
    mass: String,
    power: isize,
    dod: i32,
    orientation: i8,
    /// Its half-edges at vertices of S: (half-edge, vertex, s), `s = -1` on
    /// the source half-edge.
    hedges: Vec<(usize, usize, i64)>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Owner {
    Vertex(usize),
    Edge(usize),
}

/// One monomial of a rule: in the `+` class, the coefficient representative,
/// the monomial without its coefficient.
type Monomial = (bool, String, Atom);

/// The coefficient representative: `c` if `Re c > 0`, or `Re c = 0` and
/// `Im c > 0`; `-c` otherwise. `None` outside Q[i].
fn coefficient_class(coefficient: AtomView<'_>) -> Option<(bool, Atom)> {
    let AtomView::Num(number) = coefficient else {
        return None;
    };
    let positive = match number.get_coeff_view() {
        CoefficientView::Natural(re, _, im, _) => re > 0 || (re == 0 && im > 0),
        CoefficientView::Large(re, im) => {
            (!re.is_negative() && !re.is_zero())
                || (re.is_zero() && !im.is_negative() && !im.is_zero())
        }
        _ => return None,
    };
    let owned = coefficient.to_owned();
    Some((positive, if positive { owned } else { owned * -1 }))
}

fn monomials(rule: &Atom) -> Option<Vec<Monomial>> {
    let mut out = Vec::new();
    for term in terms_of(rule) {
        let coefficient = match term.as_view() {
            AtomView::Num(_) => Some(term.clone()),
            AtomView::Mul(mul) => {
                let numbers = mul
                    .iter()
                    .filter(|factor| matches!(factor, AtomView::Num(_)))
                    .collect::<Vec<_>>();
                match numbers.len() {
                    0 => None,
                    1 => Some(numbers[0].to_owned()),
                    _ => return None,
                }
            }
            _ => None,
        };
        let (positive, representative, monomial) = match coefficient {
            None => (true, Atom::num(1), term.clone()),
            Some(coefficient) => {
                let (positive, representative) = coefficient_class(coefficient.as_view())?;
                (positive, representative, &term / &coefficient)
            }
        };
        out.push((positive, representative.to_canonical_string(), monomial));
    }
    Some(out)
}

/// The spinor analysis of a quotiented line candidate (design §1.3 step 6).
/// `Some(ports)`: its rules are consistently oriented gamma chains, connected,
/// with `ports` its public spinor slots.
fn analyse_line(
    rules: &[(Owner, Vec<Monomial>)],
    on_line: &BTreeSet<usize>,
    bispinor: LibraryRep,
) -> Option<Vec<Atom>> {
    let is_bis = |view: AtomView<'_>| {
        LibrarySlot::<Aind>::try_from(view).is_ok_and(|slot| slot.rep().rep == bispinor)
    };
    fn has_bis(view: AtomView<'_>, is_bis: &dyn Fn(AtomView<'_>) -> bool) -> bool {
        if is_bis(view) {
            return true;
        }
        match view {
            AtomView::Fun(fun) => fun.iter().any(|arg| has_bis(arg, is_bis)),
            AtomView::Add(add) => add.iter().any(|arg| has_bis(arg, is_bis)),
            AtomView::Mul(mul) => mul.iter().any(|arg| has_bis(arg, is_bis)),
            AtomView::Pow(pow) => {
                let (base, exponent) = pow.get_base_exp();
                has_bis(base, is_bis) || has_bis(exponent, is_bis)
            }
            _ => false,
        }
    }
    // Union-find over factors and rule-external spinor indices.
    let mut parent: Vec<usize> = Vec::new();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let mut external_nodes: HashMap<Atom, usize> = HashMap::new();
    let mut factor_nodes = Vec::new();
    let mut global: HashMap<Atom, Vec<u8>> = HashMap::new();
    let mut gammas = 0_usize;
    for &rule in on_line {
        let mut external_roles: HashMap<Atom, u8> = HashMap::new();
        let mut internal: HashSet<Atom> = HashSet::new();
        for (_, _, monomial) in &rules[rule].1 {
            let mut local: HashMap<Atom, Vec<(usize, u8)>> = HashMap::new();
            for factor in factors_of(monomial.as_view()) {
                if !has_bis(factor.as_view(), &is_bis) {
                    continue;
                }
                let AtomView::Fun(fun) = factor.as_view() else {
                    return None;
                };
                let args = fun.iter().collect::<Vec<_>>();
                if fun.get_symbol() != AGS.gamma
                    || args.len() != 3
                    || !is_bis(args[0])
                    || !is_bis(args[1])
                    || has_bis(args[2], &is_bis)
                {
                    return None;
                }
                gammas += 1;
                let node = parent.len();
                parent.push(node);
                factor_nodes.push(node);
                for (role, arg) in args[..2].iter().enumerate() {
                    local
                        .entry(arg.to_owned())
                        .or_default()
                        .push((node, role as u8));
                }
            }
            for (index, occurrences) in local {
                match occurrences.as_slice() {
                    [(node, role)] => {
                        if internal.contains(&index) {
                            return None;
                        }
                        if let Some(previous) = external_roles.insert(index.clone(), *role)
                            && previous != *role
                        {
                            return None;
                        }
                        let external = *external_nodes.entry(index).or_insert_with(|| {
                            parent.push(parent.len());
                            parent.len() - 1
                        });
                        let (a, b) = (find(&mut parent, *node), find(&mut parent, external));
                        parent[a] = b;
                    }
                    [(first, first_role), (second, second_role)] => {
                        if first_role == second_role || external_roles.contains_key(&index) {
                            return None;
                        }
                        internal.insert(index);
                        let (a, b) = (find(&mut parent, *first), find(&mut parent, *second));
                        parent[a] = b;
                    }
                    _ => return None,
                }
            }
        }
        for (index, role) in external_roles {
            global.entry(index).or_default().push(role);
        }
    }
    if gammas == 0 {
        return None;
    }
    let root = find(&mut parent, factor_nodes[0]);
    if factor_nodes
        .iter()
        .any(|node| find(&mut parent, *node) != root)
    {
        return None;
    }
    // No other rule may touch the line's spinor indices.
    for (rule, (_, monomials)) in rules.iter().enumerate() {
        if on_line.contains(&rule) {
            continue;
        }
        for (_, _, monomial) in monomials {
            let mut touches = false;
            monomial.visitor(&mut |part| {
                if !touches && is_bis(part) && global.contains_key(&part.to_owned()) {
                    touches = true;
                }
                !touches
            });
            if touches {
                return None;
            }
        }
    }
    let mut ports = Vec::new();
    for (index, roles) in global {
        match roles.as_slice() {
            [_] => ports.push(index),
            [a, b] if a != b => {}
            _ => return None,
        }
    }
    ports.sort();
    Some(ports)
}

impl CountertermIdentity {
    /// The graph certificate of an R* child Δ (`GAMMALOOP_DELTA_KEY=graph`),
    /// or why the lookup keeps the legacy certificate.
    pub(crate) fn for_graph_key(
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
    ) -> Result<std::result::Result<Self, Fallback>> {
        ensure!(
            graph.is_connected(support),
            "counterterm identity needs a connected support"
        );
        let boundary = graph.dummy_less_full_crown(support);
        let conservation_holds = graph.underlying.full_crown(support) == boundary;
        let cooking = Self::cooking();
        let coordinates = [
            (
                feynkit_graph::symbols::loop_momentum(),
                graph
                    .loop_momentum_basis
                    .loop_edges
                    .iter()
                    .copied()
                    .collect::<Vec<_>>(),
            ),
            (
                feynkit_graph::symbols::external_momentum(),
                graph
                    .loop_momentum_basis
                    .ext_edges
                    .iter()
                    .copied()
                    .collect::<Vec<_>>(),
            ),
        ]
        .into_iter()
        .flat_map(|(head, edges)| {
            edges
                .into_iter()
                .enumerate()
                .map(move |(coordinate, edge)| ((head, coordinate), usize::from(edge)))
        })
        .collect::<BTreeMap<_, _>>();
        let numerator = graph
            .numerator(support, &graph.empty_subgraph())
            .to_d_dim(GS.dim)
            .get_single_atom()?;
        let mut public = BTreeSet::new();
        {
            let tensor = SymbolicTensor::infer(cooking.try_cook(numerator.as_view())?)?;
            for slot in tensor.structure().logical_slots() {
                let PartialIndex::Explicit(index) = slot.aind else {
                    return Err(eyre!(
                        "counterterm cache needs explicitly named boundary tensor slots"
                    ));
                };
                let atom = cooking.uncook(index.to_atom().as_view());
                public.insert((slot.rep().rep.base(), Aind::try_from(atom.as_view())?));
            }
        }

        // Vertices, edges and half-edges of S ∪ crown.
        let vertices = graph
            .underlying
            .iter_nodes_of(support)
            .map(|(vertex, _, data)| (usize::from(vertex), data.dod.value))
            .collect::<BTreeMap<_, _>>();
        let boundary_edges = graph
            .iter_edges_of(&boundary)
            .map(|(_, edge, _)| usize::from(edge))
            .collect::<BTreeSet<_>>();
        let all = support.union(&boundary);
        let mut edges: BTreeMap<usize, EdgeFacts> = BTreeMap::new();
        for (pair, edge, data) in graph.iter_edges_of(&all) {
            ensure!(
                !data.data.is_dummy,
                "dummy propagator in counterterm support"
            );
            let edge_id = usize::from(edge);
            let is_boundary = boundary_edges.contains(&edge_id);
            let self_conjugate = data.data.particle.is_self_antiparticle(model);
            let halfedges = match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    vec![source, sink]
                }
                HedgePair::Unpaired { hedge, .. } => vec![hedge],
            };
            let mut hedges = Vec::new();
            for hedge in halfedges {
                if !all.includes(&hedge) {
                    continue;
                }
                let vertex = usize::from(graph.underlying.node_id(hedge));
                if !vertices.contains_key(&vertex) {
                    continue;
                }
                let s = if graph.underlying.flow(hedge) == Flow::Source {
                    -1
                } else {
                    1
                };
                hedges.push((hedge.0, vertex, s));
            }
            edges.insert(
                edge_id,
                EdgeFacts {
                    boundary: is_boundary,
                    self_conjugate,
                    quark: data.data.particle.is_fermion(model) && !self_conjugate,
                    particle: data.data.particle().map(|particle| particle.index()),
                    mass: if is_boundary {
                        Atom::Zero.to_canonical_string()
                    } else {
                        data.data.mass_atom(model).to_canonical_string()
                    },
                    power: if is_boundary {
                        0
                    } else {
                        data.data.extra_data.vakint_edge_power.unwrap_or(1)
                    },
                    dod: if is_boundary { 0 } else { data.data.dod.value },
                    orientation: match data.orientation {
                        _ if self_conjugate => 0,
                        Orientation::Undirected => 0,
                        Orientation::Default => 1,
                        Orientation::Reversed => -1,
                    },
                    hedges,
                },
            );
        }

        // The local rules, with half-edge momenta, as signed monomials.
        let mut rules: Vec<(Owner, Vec<Monomial>)> = Vec::new();
        let sources = graph
            .underlying
            .iter_nodes_of(support)
            .map(|(vertex, _, data)| (Owner::Vertex(usize::from(vertex)), &data.num.value))
            .chain(
                graph
                    .iter_edges_of(support)
                    .map(|(_, edge, data)| (Owner::Edge(usize::from(edge)), &data.data.num.value)),
            )
            .collect::<Vec<_>>();
        for (owner, rule) in sources {
            let rule = Self::edge_momenta(rule, &coordinates)?;
            let rule = cooking.try_cook(rule.as_view())?;
            let mut nonincident = false;
            let rule = rule.replace_map(|term, _, out| {
                if let AtomView::Fun(fun) = term
                    && fun.get_symbol() == GS.emr_mom
                {
                    let mut args = fun.iter();
                    let Some(edge) = args.next().and_then(|arg| usize::try_from(arg).ok())
                    else {
                        nonincident = true;
                        return;
                    };
                    let rest = args.collect::<Vec<_>>();
                    let Some(facts) = edges.get(&edge) else {
                        nonincident = true;
                        return;
                    };
                    let at = facts
                        .hedges
                        .iter()
                        .filter(|(_, vertex, _)| match owner {
                            Owner::Vertex(owner) => *vertex == owner,
                            Owner::Edge(owner) => owner == edge,
                        })
                        .collect::<Vec<_>>();
                    let momentum = |hedge: usize| {
                        let mut function =
                            FunctionBuilder::new(*HEDGE_MOMENTUM).add_arg(hedge as i64);
                        for arg in &rest {
                            function = function.add_arg(*arg);
                        }
                        function.finish()
                    };
                    match at.as_slice() {
                        [(hedge, _, s)] => **out = momentum(*hedge) * *s,
                        [(first, _, s1), (second, _, s2)] => {
                            **out = (momentum(*first) * *s1 + momentum(*second) * *s2)
                                / Atom::num(2)
                        }
                        _ => nonincident = true,
                    }
                }
            });
            if nonincident {
                return Ok(Err(Fallback::NonIncident));
            }
            let Some(monomials) = monomials(&rule.expand()) else {
                return Ok(Err(Fallback::Coefficient));
            };
            rules.push((owner, monomials));
        }

        // Quark lines: components of the quark edges of S at vertices with
        // exactly two quark half-edges in S ∪ crown.
        let bispinor = bis_rep();
        let mut quark_hedges: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
        for (edge, facts) in &edges {
            if facts.quark {
                for (hedge, vertex, _) in &facts.hedges {
                    quark_hedges.entry(*vertex).or_default().push((*hedge, *edge));
                }
            }
        }
        let mut component: BTreeMap<usize, usize> =
            quark_hedges.keys().map(|vertex| (*vertex, *vertex)).collect();
        fn root(component: &mut BTreeMap<usize, usize>, mut x: usize) -> usize {
            while component[&x] != x {
                x = component[&x];
            }
            x
        }
        for facts in edges.values() {
            if facts.quark && !facts.boundary && facts.hedges.len() == 2 {
                let a = root(&mut component, facts.hedges[0].1);
                let b = root(&mut component, facts.hedges[1].1);
                component.insert(a, b);
            }
        }
        let mut lines: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
        for vertex in quark_hedges.keys().copied().collect::<Vec<_>>() {
            let r = root(&mut component, vertex);
            lines.entry(r).or_default().insert(vertex);
        }
        // line id -> (rule indices, open?, ports)
        let mut quotiented: Vec<(BTreeSet<usize>, bool, Vec<(LibraryRep, Aind)>)> = Vec::new();
        let mut rule_line: Vec<Option<usize>> = vec![None; rules.len()];
        let mut dropped_edges: BTreeSet<usize> = BTreeSet::new();
        let mut not_quotiented = 0_u64;
        let mut quotiented_vertices: BTreeSet<usize> = BTreeSet::new();
        for line_vertices in lines.values() {
            let degrees_ok = line_vertices
                .iter()
                .all(|vertex| quark_hedges[vertex].len() == 2);
            let crown = line_vertices
                .iter()
                .flat_map(|vertex| quark_hedges[vertex].iter())
                .filter(|(_, edge)| edges[edge].boundary)
                .count();
            let line_edges = edges
                .iter()
                .filter(|(_, facts)| {
                    facts.quark
                        && !facts.boundary
                        && facts
                            .hedges
                            .iter()
                            .any(|(_, vertex, _)| line_vertices.contains(vertex))
                })
                .map(|(edge, _)| *edge)
                .collect::<BTreeSet<_>>();
            let on_line = rules
                .iter()
                .enumerate()
                .filter(|(_, (owner, _))| match owner {
                    Owner::Vertex(vertex) => line_vertices.contains(vertex),
                    Owner::Edge(edge) => line_edges.contains(edge),
                })
                .map(|(rule, _)| rule)
                .collect::<BTreeSet<_>>();
            let analysed = if degrees_ok && (crown == 0 || crown == 2) {
                analyse_line(&rules, &on_line, bispinor)
            } else {
                None
            };
            let ports = analysed.and_then(|ports| {
                let ports = ports
                    .iter()
                    .map(|port| {
                        let slot = LibrarySlot::<Aind>::try_from(port.as_view()).ok()?;
                        let raw = cooking.uncook(slot.aind.to_atom().as_view());
                        let key = (slot.rep().rep.base(), Aind::try_from(raw.as_view()).ok()?);
                        public.contains(&key).then_some(key)
                    })
                    .collect::<Option<Vec<_>>>()?;
                (ports.len() == crown).then_some(ports)
            });
            let Some(ports) = ports else {
                not_quotiented += 1;
                continue;
            };
            let id = quotiented.len();
            for rule in &on_line {
                rule_line[*rule] = Some(id);
            }
            dropped_edges.extend(line_edges.iter().copied());
            quotiented_vertices.extend(line_vertices.iter().copied());
            quotiented.push((on_line, crown == 2, ports));
        }
        // A crown quark edge loses its direction when every line it ends is quotiented.
        for (edge, facts) in &edges {
            if facts.quark
                && facts.boundary
                && !facts.hedges.is_empty()
                && facts
                    .hedges
                    .iter()
                    .all(|(_, vertex, _)| quotiented_vertices.contains(vertex))
            {
                dropped_edges.insert(*edge);
            }
        }
        delta_key::add(&delta_key::C.line_not_quotiented, not_quotiented);
        let open = quotiented.iter().filter(|(_, open, _)| *open).count() as u64;
        delta_key::add(&delta_key::C.lines_open_quotiented, open);
        delta_key::add(
            &delta_key::C.lines_closed_quotiented,
            quotiented.len() as u64 - open,
        );

        // The certificate.
        let mut builder = IdentityBuilder {
            graph: IncidenceGraph::new(),
            indices: BTreeMap::new(),
            edges: BTreeMap::new(),
            hedges: BTreeMap::new(),
            vertices: BTreeMap::new(),
            public,
            boundary_signs: BTreeMap::new(),
            momentum_signs: BTreeMap::new(),
            momentum_owners: BTreeMap::new(),
            momentum_sources: BTreeMap::new(),
            topology_labels: Vec::new(),
            index_owners: BTreeMap::new(),
            coordinates: coordinates.clone(),
            cooking,
        };
        for (vertex, dod) in &vertices {
            let node = builder.graph.add_node(NodeColor::Vertex(*dod));
            builder.vertices.insert(*vertex, node);
        }
        let mut hedge_sign: BTreeMap<usize, i64> = BTreeMap::new();
        for (edge, facts) in &edges {
            let undirected = facts.self_conjugate || dropped_edges.contains(edge);
            let node = builder.graph.add_node(NodeColor::EdgeP {
                particle: facts.particle,
                mass: facts.mass.clone(),
                power: facts.power,
                dod: facts.dod,
                boundary: facts.boundary,
                orientation: if undirected { 0 } else { facts.orientation },
            });
            builder.edges.insert(*edge, node);
            for (hedge, vertex, s) in &facts.hedges {
                let hnode = builder.graph.add_node(NodeColor::HalfEdge {
                    boundary: facts.boundary,
                    flow: if undirected { None } else { Some(*s == -1) },
                });
                builder.hedges.insert(*hedge, hnode);
                hedge_sign.insert(*hedge, *s);
                builder.link(node, hnode, Link::Incidence);
                builder.link(builder.vertices[vertex], hnode, Link::Incidence);
            }
        }
        let mut sign_classes = Vec::new();
        let mut factors: Vec<(usize, usize, [usize; 2])> = Vec::new();
        for (rule, (owner, monomials)) in rules.iter().enumerate() {
            let owner_node = match owner {
                Owner::Vertex(vertex) => builder.vertices[vertex],
                Owner::Edge(edge) => builder.edges[edge],
            };
            let rule_node = builder.graph.add_node(NodeColor::Rule);
            builder.link(owner_node, rule_node, Link::Numerator);
            let plus = builder.graph.add_node(NodeColor::SignClass);
            let minus = builder.graph.add_node(NodeColor::SignClass);
            builder.link(rule_node, plus, Link::Numerator);
            builder.link(rule_node, minus, Link::Numerator);
            for (positive, coefficient, monomial) in monomials {
                let term = builder
                    .graph
                    .add_node(NodeColor::Term(coefficient.clone()));
                builder.link(if *positive { plus } else { minus }, term, Link::Numerator);
                let root =
                    builder.expression_graph(monomial.as_view(), rule_line[rule], &mut factors)?;
                builder.link(term, root, Link::Numerator);
            }
            sign_classes.push((plus, minus));
        }
        let canonical = builder.graph.canonize();
        let label = |node: usize| canonical.vertex_map[node];
        let sigma = sign_classes
            .iter()
            .map(|(plus, minus)| if label(*plus) < label(*minus) { 1 } else { -1 })
            .product::<i64>();
        // Line orientation bits: the gamma factor with the smallest label.
        let mut line_bits: Vec<(usize, bool, bool, Vec<usize>)> = Vec::new();
        for (id, (_, open, ports)) in quotiented.iter().enumerate() {
            let reference = factors
                .iter()
                .filter(|(line, _, _)| *line == id)
                .min_by_key(|(_, node, _)| label(*node));
            let Some((_, node, slots)) = reference else {
                continue;
            };
            let crossed = label(slots[0]) > label(slots[1]);
            let mut port_labels = ports
                .iter()
                .map(|key| builder.indices.get(key).map(|node| label(*node)))
                .collect::<Option<Vec<_>>>()
                .ok_or_else(|| eyre!("an open-line port has no index node"))?;
            port_labels.sort();
            line_bits.push((label(*node), *open, crossed, port_labels));
        }
        line_bits.sort();
        let crossed_open = line_bits
            .iter()
            .filter(|(_, open, crossed, _)| *open && *crossed)
            .map(|(_, _, _, ports)| [ports[0], ports[1]])
            .collect::<Vec<_>>();
        let ports = builder
            .indices
            .iter()
            .filter(|(index, _)| builder.public.contains(index))
            .map(|(index, node)| (*index, label(*node)))
            .collect();
        // Boundary momenta at the crown half-edge with the smaller label.
        let mut momenta = BTreeMap::new();
        for (edge, facts) in &edges {
            if !facts.boundary {
                continue;
            }
            let Some((hedge, _, s)) = facts
                .hedges
                .iter()
                .min_by_key(|(hedge, _, _)| label(builder.hedges[hedge]))
            else {
                continue;
            };
            momenta.insert(*edge, (label(builder.hedges[hedge]), *s));
        }
        let conservation = conservation_holds.then(|| {
            edges
                .iter()
                .filter(|(_, facts)| facts.boundary && facts.hedges.len() == 1)
                .map(|(edge, facts)| (*edge, facts.hedges[0].2))
                .collect::<Vec<_>>()
        });
        // Features.
        let features = {
            use std::hash::{Hash, Hasher};
            let mut hasher = std::collections::hash_map::DefaultHasher::new();
            for node in builder.vertices.values() {
                label(*node).hash(&mut hasher);
            }
            for node in builder.hedges.values() {
                label(*node).hash(&mut hasher);
            }
            let crown_hedges = edges
                .values()
                .filter(|facts| facts.boundary)
                .flat_map(|facts| facts.hedges.iter())
                .collect::<BTreeSet<_>>();
            let mut crown = crown_hedges
                .iter()
                .map(|(hedge, _, s)| (label(builder.hedges[hedge]), *s))
                .collect::<Vec<_>>();
            crown.sort();
            let mut boson = edges
                .iter()
                .filter(|(_, facts)| {
                    facts.self_conjugate && !facts.boundary && facts.hedges.len() == 2
                })
                .map(|(edge, facts)| {
                    let (source, sink) = if facts.hedges[0].2 == -1 {
                        (facts.hedges[0].0, facts.hedges[1].0)
                    } else {
                        (facts.hedges[1].0, facts.hedges[0].0)
                    };
                    (
                        label(builder.edges[edge]),
                        label(builder.hedges[&source]) < label(builder.hedges[&sink]),
                    )
                })
                .collect::<Vec<_>>();
            boson.sort();
            let mut routing = graph
                .loop_momentum_basis
                .loop_edges
                .iter()
                .map(|edge| usize::from(*edge))
                .filter(|edge| edges.get(edge).is_some_and(|facts| !facts.boundary))
                .map(|edge| label(builder.edges[&edge]))
                .collect::<Vec<_>>();
            routing.sort();
            OccurrenceFeatures {
                sigma,
                open: line_bits
                    .iter()
                    .filter(|(_, open, _, _)| *open)
                    .map(|(_, _, crossed, _)| *crossed)
                    .collect(),
                closed: line_bits
                    .iter()
                    .filter(|(_, open, _, _)| !*open)
                    .map(|(_, _, crossed, _)| *crossed)
                    .collect(),
                crown: crown.into_iter().map(|(_, s)| s).collect(),
                relabel: hasher.finish(),
                legs: crown_hedges
                    .iter()
                    .map(|(hedge, _, _)| label(builder.hedges[hedge]))
                    .collect(),
                boson: boson.into_iter().map(|(_, bit)| bit).collect(),
                routing,
                legacy: None,
            }
        };
        let reservations = graph
            .iter_edges()
            .map(|(_, _, edge)| edge.data.num.value.clone())
            .chain(
                graph
                    .underlying
                    .iter_nodes()
                    .map(|(_, _, vertex)| vertex.num.value.clone()),
            )
            .chain([
                graph.global_prefactor.num.clone(),
                graph.global_prefactor.projector.clone(),
            ])
            .collect();
        Ok(Ok(Self {
            key: CountertermKey {
                model: model.fingerprint(),
                graph: canonical.graph,
            },
            ports,
            momenta,
            coordinates,
            reservations,
            graph_transport: Some(GraphTransport {
                sigma,
                crossed_open,
                features,
                conservation,
            }),
        }))
    }

    pub(crate) fn graph_transport(&self) -> Option<&GraphTransport> {
        self.graph_transport.as_ref()
    }

    pub(crate) fn set_legacy_digest(&mut self, digest: u64) {
        if let Some(transport) = &mut self.graph_transport {
            transport.features.legacy = Some(digest);
        }
    }

    /// Verify: whether the stored canonical value, transported to this
    /// occurrence, agrees with `fresh`: `(identical, equal)` — identical in
    /// the canonical frame, or equal after momentum conservation (when its
    /// precondition holds) and gamma/colour simplification of the difference.
    pub(crate) fn verify_agree(
        &self,
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        stored: &Atom,
        fresh: &Atom,
    ) -> Result<(bool, bool)> {
        use idenso::{
            color::ColorSimplifySettings,
            dirac::GammaSimplifySettings,
            tensor::{AlgebraSettings, SymbolicTensor},
        };
        if let Ok(canonical) = self.canonicalize_delta(fresh)
            && (&canonical - stored).expand().is_zero()
        {
            return Ok((true, true));
        }
        let restored = self.restore_delta(stored)?;
        if (&restored - fresh).expand().is_zero() {
            return Ok((false, true));
        }
        // In this occurrence the boundary momenta are linear combinations of
        // the parent's loop-momentum basis: two writers' representatives of
        // one value agree after that substitution (momentum conservation in
        // the parent, also through a dummy edge at S).
        let (restored, fresh) = (lmb_momenta(graph, &restored), lmb_momenta(graph, fresh));
        let mut difference = (&restored - &fresh).expand();
        delta_key::add(&delta_key::C.nf_conservation_used, 1);
        if difference.is_zero() {
            return Ok((false, true));
        }
        if let Some(legs) = self
            .graph_transport
            .as_ref()
            .and_then(|transport| transport.conservation.as_ref())
            .filter(|legs| legs.len() >= 2)
        {
            difference = conserve(&difference, legs);
            if difference.is_zero() {
                return Ok((false, true));
            }
        } else {
            delta_key::add(&delta_key::C.nf_conservation_skipped, 1);
        }
        let _ = (support, model);
        let cooking = Self::cooking();
        let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(difference.as_view())?)?
            .undo_trace()?
            .undo_chain()?
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(GammaSimplifySettings::default()),
                color: Some(ColorSimplifySettings::default()),
                ..Default::default()
            })?
            .to_dots()?;
        let reduced = cooking.uncook(tensor.expression().as_view()).expand();
        if reduced.is_zero() {
            return Ok((false, true));
        }
        let canonical = SymbolicTensor::infer(cooking.try_cook(reduced.as_view())?)?.canonize()?;
        let canonical = cooking.uncook(canonical.expression().as_view()).expand();
        Ok((false, canonical.together().is_zero()))
    }
}

/// Every edge momentum `Q(e, ..)` as its signature in `graph`'s loop-momentum
/// basis: a combination of the basis edges' `Q(l, ..)` and external `Q(x, ..)`.
fn lmb_momenta(graph: &Graph, value: &Atom) -> Atom {
    let basis = &graph.loop_momentum_basis;
    value
        .replace_map(|term, _, out| {
            if let AtomView::Fun(fun) = term
                && fun.get_symbol() == GS.emr_mom
            {
                let mut args = fun.iter();
                if let Some(edge) = args
                    .next()
                    .and_then(|arg| usize::try_from(arg).ok())
                    .map(linnet::half_edge::involution::EdgeIndex::from)
                    && basis.edge_signatures.get(edge).is_some()
                {
                    let rest = args.map(|arg| arg.to_owned()).collect::<Vec<_>>();
                    **out = basis.loop_atom::<Atom>(edge, GS.emr_mom, &rest, true)
                        + basis.ext_atom::<Atom>(edge, GS.emr_mom, &rest, true);
                }
            }
        })
        .expand()
}

/// `Q(e_last) -> -s_last Σ_{j≠last} s_j Q(e_j)` (momentum conservation of the
/// child; `e_last` the largest edge id).
fn conserve(value: &Atom, legs: &[(usize, i64)]) -> Atom {
    let (last, s_last) = *legs.iter().max().expect("at least two legs");
    value
        .replace_map(|term, _, out| {
            if let AtomView::Fun(fun) = term
                && fun.get_symbol() == GS.emr_mom
            {
                let mut args = fun.iter();
                if args.next().and_then(|arg| usize::try_from(arg).ok()) == Some(last) {
                    let rest = args.collect::<Vec<_>>();
                    let mut sum = Atom::Zero;
                    for (edge, s) in legs {
                        if *edge == last {
                            continue;
                        }
                        let mut function = FunctionBuilder::new(GS.emr_mom).add_arg(*edge as i64);
                        for arg in &rest {
                            function = function.add_arg(*arg);
                        }
                        sum = sum + function.finish() * *s;
                    }
                    **out = sum * -s_last;
                }
            }
        })
        .expand()
}

impl IdentityBuilder {
    /// [`IdentityBuilder::expression`] of the graph certificate: portable
    /// colours, half-edge momenta, and unordered spinor slots of the gammas of
    /// a quotiented quark line (`line`), recorded in `factors`.
    fn expression_graph(
        &mut self,
        expression: AtomView<'_>,
        line: Option<usize>,
        factors: &mut Vec<(usize, usize, [usize; 2])>,
    ) -> Result<usize> {
        if let Ok(slot) = LibrarySlot::<Aind>::try_from(expression) {
            let raw = self.cooking.uncook(slot.aind.to_atom().as_view());
            let original = Aind::try_from(raw.as_view())?;
            let index = self.index(slot.rep().rep, original);
            let node = self.graph.add_node(NodeColor::SlotP(
                slot.rep().to_symbolic([]).to_canonical_string(),
            ));
            self.link(node, index, Link::Index);
            return Ok(node);
        }
        let (color, args, ordered) = match expression {
            AtomView::Fun(fun) if fun.get_symbol() == *HEDGE_MOMENTUM => {
                let mut args = fun.iter();
                let hedge = usize::try_from(
                    args.next()
                        .ok_or_else(|| eyre!("missing half-edge momentum label"))?,
                )
                .map_err(|error| eyre!(error))?;
                let owner = *self
                    .hedges
                    .get(&hedge)
                    .ok_or_else(|| eyre!("half-edge momentum outside S and its crown"))?;
                let node = self.graph.add_node(NodeColor::HedgeMomentum);
                self.link(node, owner, Link::Momentum);
                for (position, arg) in args.enumerate() {
                    let child = self.expression_graph(arg, line, factors)?;
                    self.link(node, child, Link::Argument(position));
                }
                return Ok(node);
            }
            AtomView::Fun(fun) if fun.get_symbol() == GS.emr_mom => {
                return Err(eyre!("an edge momentum survived the half-edge rewrite"));
            }
            AtomView::Fun(fun)
                if line.is_some() && fun.get_symbol() == AGS.gamma && fun.get_nargs() == 3 =>
            {
                let node = self
                    .graph
                    .add_node(NodeColor::FunctionP(fun.get_symbol().get_name().to_owned()));
                let args = fun.iter().collect::<Vec<_>>();
                let row = self.expression_graph(args[0], line, factors)?;
                let column = self.expression_graph(args[1], line, factors)?;
                let lorentz = self.expression_graph(args[2], line, factors)?;
                self.link(node, row, Link::Unordered);
                self.link(node, column, Link::Unordered);
                self.link(node, lorentz, Link::Argument(2));
                factors.push((line.expect("checked"), node, [row, column]));
                return Ok(node);
            }
            AtomView::Fun(fun) => (
                NodeColor::FunctionP(fun.get_symbol().get_name().to_owned()),
                fun.iter().collect::<Vec<_>>(),
                !fun.get_symbol().is_symmetric(),
            ),
            AtomView::Add(sum) => (NodeColor::Sum, sum.iter().collect(), false),
            AtomView::Mul(product) => (NodeColor::Product, product.iter().collect(), false),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                (NodeColor::Power, vec![base, exponent], true)
            }
            _ => {
                return Ok(self
                    .graph
                    .add_node(NodeColor::LeafP(expression.to_canonical_string())));
            }
        };
        let node = self.graph.add_node(color);
        for (position, arg) in args.into_iter().enumerate() {
            let child = self.expression_graph(arg, line, factors)?;
            self.link(
                node,
                child,
                if ordered {
                    Link::Argument(position)
                } else {
                    Link::Unordered
                },
            );
        }
        Ok(node)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        finalized_runtime_dot,
        graph::{Autogen, FeynmanGraph, parse::IntoFinalizedRuntimeGraph},
        initialisation::test_initialise,
        utils::load_generic_model,
    };
    use linnet::half_edge::{NodeIndex, involution::EdgeIndex};
    use spenso::{
        network::library::symbolic::ETS,
        structure::representation::{Minkowski, RepName},
    };

    fn support(graph: &Graph) -> SuBitGraph {
        graph
            .no_dummy()
            .subtract(&graph.external_filter::<SuBitGraph>())
    }

    #[derive(Clone, Copy)]
    enum R {
        Mink,
        Bis,
    }

    fn mink() -> R {
        R::Mink
    }

    fn bis() -> R {
        R::Bis
    }

    fn slot(rep: R, index: Aind) -> Atom {
        match rep {
            R::Mink => Minkowski {}.new_rep(GS.dim).to_symbolic([Atom::from(index)]),
            R::Bis => Bispinor {}.new_rep(4).to_symbolic([Atom::from(index)]),
        }
    }

    fn hedge(h: usize) -> Aind {
        Aind::Hedge(h as u16, 1)
    }

    /// (source half-edge, sink half-edge) of every edge; `None` for the
    /// missing end of an external edge.
    fn ends(graph: &Graph) -> BTreeMap<usize, (Option<usize>, Option<usize>)> {
        graph
            .iter_edges()
            .map(|(pair, edge, _)| {
                let ends = match pair {
                    HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                        (Some(source.0), Some(sink.0))
                    }
                    HedgePair::Unpaired { hedge, .. } => {
                        if graph.underlying.flow(hedge) == Flow::Source {
                            (Some(hedge.0), None)
                        } else {
                            (None, Some(hedge.0))
                        }
                    }
                };
                (usize::from(edge), ends)
            })
            .collect()
    }

    fn vertex_of(graph: &Graph, h: usize) -> NodeIndex {
        graph
            .underlying
            .node_id(linnet::half_edge::involution::Hedge(h))
    }

    /// Feynman-like rules written along each quark edge's own orientation:
    /// vertex `gamma(bis(in), bis(out), mink(gluon))` (times `weights` of the
    /// vertex), quark propagator `gamma(bis(source), bis(sink), mink(e)) Q(e, mink(e))`,
    /// gluon propagator `g(mink(source), mink(sink))`.
    fn dress(graph: &mut Graph, quarks: &[usize], weights: &BTreeMap<usize, i64>) {
        let ends = ends(graph);
        let mut vertex_in: BTreeMap<NodeIndex, usize> = BTreeMap::new();
        let mut vertex_out: BTreeMap<NodeIndex, usize> = BTreeMap::new();
        let mut vertex_gluon: BTreeMap<NodeIndex, usize> = BTreeMap::new();
        for (edge, (source, sink)) in &ends {
            for (h, is_source) in [(source, true), (sink, false)] {
                let Some(h) = h else { continue };
                let v = vertex_of(graph, *h);
                if quarks.contains(edge) {
                    if is_source {
                        vertex_out.insert(v, *h);
                    } else {
                        vertex_in.insert(v, *h);
                    }
                } else {
                    vertex_gluon.insert(v, *h);
                }
            }
        }
        for (v, gluon) in &vertex_gluon {
            let rule = function!(
                AGS.gamma,
                slot(bis(), hedge(vertex_in[v])),
                slot(bis(), hedge(vertex_out[v])),
                slot(mink(), hedge(*gluon))
            ) * weights.get(&usize::from(*v)).copied().unwrap_or(1);
            graph.underlying[*v].num = Autogen::explicit(rule);
        }
        for (edge, (source, sink)) in &ends {
            let (Some(source), Some(sink)) = (source, sink) else {
                continue;
            };
            let index = slot(mink(), Aind::Edge(*edge as u16, 1));
            let rule = if quarks.contains(edge) {
                function!(
                    AGS.gamma,
                    slot(bis(), hedge(*source)),
                    slot(bis(), hedge(*sink)),
                    &index
                ) * function!(GS.emr_mom, *edge as i64, &index)
            } else {
                function!(
                    ETS.metric,
                    slot(mink(), hedge(*source)),
                    slot(mink(), hedge(*sink))
                )
            };
            graph.underlying[EdgeIndex::from(*edge)].num = Autogen::explicit(rule);
        }
    }

    fn graph_identity(graph: &Graph, model: &Model) -> Result<CountertermIdentity> {
        CountertermIdentity::for_graph_key(graph, &support(graph), model)?
            .map_err(|fallback| eyre!("graph key fell back: {fallback:?}"))
    }

    /// Relabelling, a different loop-momentum basis and a negated local rule:
    /// one graph certificate, the legacy ones differ, and the transport
    /// carries the rule sign and the momentum orientation of the legs.
    #[test]
    fn graph_key_ignores_numbering_routing_and_rule_sign_and_transports() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let bubble = |permuted: bool| -> Result<Graph> {
            let mut graph: Graph = if permuted {
                finalized_runtime_dot!(
                    digraph G {
                        graph [projector=19]; node [num=1]; edge [particle="scalar_1" num=1];
                        e [style=invis]; B; A;
                        A:0 -> e [id=0 source="{ufo_order:0}"];
                        e -> B:1 [id=1 sink="{ufo_order:0}"];
                        B -> A [id=3 source="{ufo_order:1}" sink="{ufo_order:1}"];
                        A -> B [id=2 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
                    }, &model
                )?
            } else {
                finalized_runtime_dot!(
                    digraph G {
                        graph [projector=7]; node [num=1]; edge [particle="scalar_1" num=1];
                        e [style=invis]; A; B;
                        e -> A:0 [id=2 sink="{ufo_order:0}"];
                        B:1 -> e [id=3 source="{ufo_order:0}"];
                        A -> B [id=0 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                        A -> B [id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
                    }, &model
                )?
            };
            for edge in graph.iter_edges().map(|(_, edge, _)| edge).collect::<Vec<_>>() {
                let data = &mut graph.underlying[edge];
                data.particle = data.particle.clone().override_mass(Some(Atom::Zero));
            }
            Ok(graph)
        };
        // External edge `ext`, internal edges `first` (rule momentum) and `second`.
        let decorate = |graph: &mut Graph, ext: [usize; 2], first: usize, sign: i64| {
            let ends = ends(graph);
            let mut ports = Vec::new();
            for (position, edge) in ext.into_iter().enumerate() {
                let h = ends[&edge].0.or(ends[&edge].1).unwrap();
                let v = vertex_of(graph, h);
                let port = hedge(h);
                // Each rule carries Q(first) flowing into its vertex (the
                // same physics in every numbering and orientation); B's rule
                // is negated by `sign`.
                let (source, _) = ends[&first];
                let inflow = if source.is_some_and(|s| vertex_of(graph, s) == v) {
                    -1
                } else {
                    1
                };
                let rule = if position == 0 {
                    function!(GS.emr_mom, first as i64, slot(mink(), port)) * (2 * inflow)
                } else {
                    function!(
                        ETS.metric,
                        slot(mink(), port),
                        slot(mink(), Aind::Dummy(77))
                    ) * function!(GS.emr_mom, first as i64, slot(mink(), Aind::Dummy(77)))
                        * (3 * sign * inflow)
                };
                graph.underlying[v].num = Autogen::explicit(rule);
                ports.push(port);
            }
            [ports[0], ports[1]]
        };
        let mut a = bubble(false)?;
        let mut b = bubble(true)?;
        let a_ports = decorate(&mut a, [2, 3], 0, 1);
        let b_ports = decorate(&mut b, [0, 1], 3, -1);
        let (ga, gb) = (graph_identity(&a, &model)?, graph_identity(&b, &model)?);
        assert_eq!(ga.key, gb.key, "one graph certificate");
        let (la, lb) = (
            CountertermIdentity::new(&a, &support(&a), &model)?,
            CountertermIdentity::new(&b, &support(&b), &model)?,
        );
        assert_ne!(la.key, lb.key, "the legacy certificate sees the rule sign");
        let (ta, tb) = (
            ga.graph_transport().unwrap(),
            gb.graph_transport().unwrap(),
        );
        assert_eq!(ta.sigma * tb.sigma, -1);
        assert_ne!(ta.features.routing, tb.features.routing);
        // Δ linear in B's rule: the restored value is the negated image.
        let q = |edge: usize, port: Aind| function!(GS.emr_mom, edge as i64, slot(mink(), port));
        let value = q(2, a_ports[0]) * q(2, a_ports[1]) / Atom::var(GS.dim_epsilon);
        let restored = gb.restore_delta(&ga.canonicalize_delta(&value)?)?;
        // Leg 2 enters A in `a`; leg 0 leaves A in `b`: Q(2) -> -Q(0).
        assert_eq!(
            restored,
            -(q(0, b_ports[0]) * q(0, b_ports[1]) / Atom::var(GS.dim_epsilon))
        );
        let restored = gb.restore_delta(&ga.canonicalize_delta(&q(2, a_ports[0]))?)?;
        assert_eq!(restored, q(0, b_ports[0]));
        Ok(())
    }

    /// A closed quark loop (a triangle, its vertices told apart by weights) and
    /// its reversal, which is no relabelling: one graph certificate (the
    /// legacy ones differ), no transport, opposite loop orientation bits, and
    /// the sign of the three reversed propagators.
    #[test]
    fn graph_key_identifies_reversed_closed_quark_loop() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("sm");
        let forward: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1]; node [num=1]; edge [particle="g" num=1];
                ext [style=invis];
                ext -> A [id=0 sink="{ufo_order:2}"];
                ext -> B [id=1 sink="{ufo_order:2}"];
                C -> ext [id=2 source="{ufo_order:2}"];
                A -> B [id=3 particle="d" lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"];
                B -> C [id=4 particle="d" source="{ufo_order:0}" sink="{ufo_order:1}"];
                C -> A [id=5 particle="d" source="{ufo_order:0}" sink="{ufo_order:1}"];
            }, "sm"
        )?;
        let reversed: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1]; node [num=1]; edge [particle="g" num=1];
                ext [style=invis];
                ext -> A [id=0 sink="{ufo_order:2}"];
                ext -> B [id=1 sink="{ufo_order:2}"];
                C -> ext [id=2 source="{ufo_order:2}"];
                B -> A [id=3 particle="d" lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"];
                C -> B [id=4 particle="d" source="{ufo_order:0}" sink="{ufo_order:1}"];
                A -> C [id=5 particle="d" source="{ufo_order:0}" sink="{ufo_order:1}"];
            }, "sm"
        )?;
        let weights = |graph: &Graph| {
            let ends = ends(graph);
            BTreeMap::from([
                (usize::from(vertex_of(graph, ends[&0].1.unwrap())), 2),
                (usize::from(vertex_of(graph, ends[&1].1.unwrap())), 3),
            ])
        };
        let (mut a, mut b) = (forward, reversed);
        let (wa, wb) = (weights(&a), weights(&b));
        dress(&mut a, &[3, 4, 5], &wa);
        dress(&mut b, &[3, 4, 5], &wb);
        let (ga, gb) = (graph_identity(&a, &model)?, graph_identity(&b, &model)?);
        assert_eq!(ga.key, gb.key, "a reversed closed quark loop has the same certificate");
        assert_ne!(
            CountertermIdentity::new(&a, &support(&a), &model)?.key,
            CountertermIdentity::new(&b, &support(&b), &model)?.key
        );
        let (ta, tb) = (
            ga.graph_transport().unwrap(),
            gb.graph_transport().unwrap(),
        );
        assert_eq!(ta.features.closed.len(), 1);
        assert_ne!(ta.features.closed, tb.features.closed);
        assert!(ta.crossed_open.is_empty() && tb.crossed_open.is_empty());
        // Three propagators change sign: σ_a σ_b = -1 (Furry's sign).
        assert_eq!(ta.sigma * tb.sigma, -1);
        let ea = ends(&a);
        let (p0, p1, p2) = (
            hedge(ea[&0].1.unwrap()),
            hedge(ea[&1].1.unwrap()),
            hedge(ea[&2].0.unwrap()),
        );
        let value = function!(ETS.metric, slot(mink(), p0), slot(mink(), p1))
            * function!(GS.emr_mom, 2, slot(mink(), p2));
        assert_eq!(gb.restore_delta(&ga.canonicalize_delta(&value)?)?, -value);
        Ok(())
    }

    /// An open quark line (self-energy) and its reversal: one certificate,
    /// opposite open-line bits, and the transposition of the chain.
    #[test]
    fn graph_key_transposes_reversed_open_quark_line() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("sm");
        let forward: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1]; node [num=1]; edge [particle="d" num=1];
                ext [style=invis];
                ext -> A [id=0 sink="{ufo_order:1}"];
                B -> ext [id=1 source="{ufo_order:0}"];
                A -> B [id=2 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"];
                A -> B [id=3 particle="g" source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, "sm"
        )?;
        let reversed: Graph = finalized_runtime_dot!(
            digraph G {
                graph [projector=1]; node [num=1]; edge [particle="d" num=1];
                ext [style=invis];
                A -> ext [id=0 source="{ufo_order:0}"];
                ext -> B [id=1 sink="{ufo_order:1}"];
                B -> A [id=2 lmb_id=0 source="{ufo_order:0}" sink="{ufo_order:1}"];
                A -> B [id=3 particle="g" source="{ufo_order:2}" sink="{ufo_order:2}"];
            }, "sm"
        )?;
        let weight_a = |graph: &Graph| {
            let ends = ends(graph);
            let h = ends[&0].0.or(ends[&0].1).unwrap();
            BTreeMap::from([(usize::from(vertex_of(graph, h)), 2)])
        };
        let (mut a, mut b) = (forward, reversed);
        let (wa, wb) = (weight_a(&a), weight_a(&b));
        dress(&mut a, &[0, 1, 2], &wa);
        dress(&mut b, &[0, 1, 2], &wb);
        let (ga, gb) = (graph_identity(&a, &model)?, graph_identity(&b, &model)?);
        assert_eq!(ga.key, gb.key, "a reversed open quark line has the same certificate");
        let (ta, tb) = (
            ga.graph_transport().unwrap(),
            gb.graph_transport().unwrap(),
        );
        assert_eq!(ta.features.open.len(), 1);
        assert_ne!(ta.features.open, tb.features.open);
        assert_eq!(ta.sigma * tb.sigma, -1, "the quark propagator changes sign");
        let (ea, eb) = (ends(&a), ends(&b));
        // ports: the quark legs' half-edges at A and at B
        let (pa, pb) = (hedge(ea[&0].1.unwrap()), hedge(ea[&1].0.unwrap()));
        let (pa2, pb2) = (hedge(eb[&0].0.unwrap()), hedge(eb[&1].1.unwrap()));
        let dummy = slot(mink(), Aind::Dummy(5));
        let value = function!(AGS.gamma, slot(bis(), pa), slot(bis(), pb), &dummy)
            * function!(GS.emr_mom, 0, &dummy);
        let restored = gb.restore_delta(&ga.canonicalize_delta(&value)?)?;
        // the fresh private dummy back to one name
        let restored = restored.replace_map(|term, _, out| {
            if let Ok(found) = LibrarySlot::<Aind>::try_from(term)
                && !matches!(found.aind, Aind::Hedge(..))
            {
                **out = dummy.clone();
            }
        });
        // Δ_b = σ_aσ_b φ(Δ_a∘ψ): the slots swapped, Q(0) reversed at A: (-1)(-1).
        let expected = function!(AGS.gamma, slot(bis(), pb2), slot(bis(), pa2), &dummy)
            * function!(GS.emr_mom, 0, &dummy);
        assert_eq!(restored, expected);
        // φ in the canonical frame: a chain is transposed, an unknown spinor
        // function is refused.
        let port = |n: i64| Bispinor {}.new_rep(4).to_symbolic([function!(*DELTA_PORT, n)]);
        let chain = function!(AGS.gamma, port(0), port(1), &dummy);
        assert_eq!(
            transpose_open_lines(&chain, &[[0, 1]]),
            Some(function!(AGS.gamma, port(1), port(0), &dummy))
        );
        let unknown = function!(symbol!("chain_of_unknown_kind"), port(0), port(1));
        assert_eq!(transpose_open_lines(&unknown, &[[0, 1]]), None);
        Ok(())
    }
}
