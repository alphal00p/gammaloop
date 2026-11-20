//! Isomorphism certificates for dressed UV subgraphs and their open tensors.
//!
//! Graph labels, routing choices and contracted index names are not physics.
//! The certificate retains incidence, particle species, masses, propagator
//! powers and the complete local numerator syntax. Its canonical vertex map
//! also supplies the boundary transport needed to reuse a computed integral.

use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    sync::LazyLock,
};

use eyre::{Result, ensure, eyre};
use idenso::{CookMode, CookSettings, tensor::SymbolicTensor};
use linnet::half_edge::{
    involution::{Flow, HedgePair, Orientation},
    subgraph::{Inclusion, SuBitGraph, SubSetOps},
};
use spenso::{
    network::parsing::ParseState,
    structure::{
        abstract_index::AbstractIndex,
        partial::{PartialIndex, PartialStructureExt},
        representation::{LibraryRep, LibrarySlot, RepName},
        slot::{IsAbstractSlot, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    function,
    graph::Graph as IncidenceGraph,
    symbol,
};

use crate::{
    graph::Graph,
    model::Model,
    numerator::aind::Aind,
    utils::GS,
    uv::{UltravioletGraph, uv_graph::UVE},
};

static DELTA_PRIVATE: LazyLock<Symbol> = LazyLock::new(|| symbol!("uv_delta_cache_private"));
static DELTA_PORT: LazyLock<Symbol> = LazyLock::new(|| symbol!("uv_delta_cache_port"));
static DELTA_MOMENTUM: LazyLock<Symbol> = LazyLock::new(|| symbol!("uv_delta_cache_momentum"));

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
enum NodeColor {
    Vertex(i32),
    Edge {
        particle: Option<usize>,
        mass: Atom,
        power: isize,
        dod: i32,
        boundary: bool,
        orientation: i8,
    },
    HalfEdge {
        boundary: bool,
        flow: Option<bool>,
    },
    Index {
        public: bool,
    },
    Slot(Atom),
    Leaf(Atom),
    Function(Symbol),
    Sum,
    Product,
    Power,
    Momentum,
    GlobalPrefactor,
    GlobalProjector,
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
enum Link {
    Incidence,
    Numerator,
    Index,
    Owner,
    Momentum,
    MomentumSource,
    Argument(usize),
    Unordered,
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub(crate) struct CountertermKey {
    model: feynkit_model::ModelFingerprint,
    graph: IncidenceGraph<NodeColor, Link>,
}

#[cfg(test)]
impl CountertermKey {
    pub(crate) fn incidence_nodes(&self) -> usize {
        self.graph.nodes().len()
    }
}

/// The key and the particular isomorphism from one occurrence to that key.
pub(crate) struct CountertermIdentity {
    pub(crate) key: CountertermKey,
    ports: BTreeMap<(LibraryRep, Aind), usize>,
    momenta: BTreeMap<usize, (usize, i64)>,
    coordinates: BTreeMap<(Symbol, usize), usize>,
    reservations: Vec<Atom>,
}

impl CountertermIdentity {
    pub(crate) fn new(graph: &Graph, support: &SuBitGraph, model: &Model) -> Result<Self> {
        Self::build(graph, support, model, false, false)
    }

    /// Bare skeletons transport integrands, including their internal momenta.
    /// Their global tensors also participate in equivalence; ordinary local
    /// counterterms deliberately retain their existing boundary-only contract.
    pub(crate) fn for_bare_skeleton(
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
    ) -> Result<Self> {
        Self::build(graph, support, model, true, false)
    }

    /// Certify a bare integral in its original external momentum frame.
    /// Fixing incoming/outgoing boundary roles prevents an undirected root
    /// reflection from choosing internal momentum signs before the complete
    /// numerator breaks that symmetry. External-leg exchange is deliberately
    /// not part of this more conservative equivalence.
    pub(crate) fn for_fixed_boundary_skeleton(
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
    ) -> Result<Self> {
        Self::build(graph, support, model, true, true)
    }

    /// `GAMMALOOP_IR_HARD`: the identity of an R* IR hard component. The marked
    /// line carries the reference mass, and each edge in `decorated` (the
    /// contracted spinney of a preparation term) is coloured as well, so the
    /// certificate separates the marked position and the preparation term.
    pub(crate) fn for_rstar_hard(
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        marked: Option<linnet::half_edge::involution::EdgeIndex>,
        decorated: &BTreeSet<usize>,
    ) -> Result<Self> {
        Self::build_decorated(graph, support, model, false, false, marked, decorated)
    }

    /// A 64-bit digest of the key, stable within one process.
    pub(crate) fn key_digest(&self) -> u64 {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        self.key.hash(&mut hasher);
        hasher.finish()
    }

    fn build(
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        bare_skeleton: bool,
        fixed_boundary: bool,
    ) -> Result<Self> {
        Self::build_decorated(
            graph,
            support,
            model,
            bare_skeleton,
            fixed_boundary,
            None,
            &BTreeSet::new(),
        )
    }

    fn build_decorated(
        graph: &Graph,
        support: &SuBitGraph,
        model: &Model,
        bare_skeleton: bool,
        fixed_boundary: bool,
        marked: Option<linnet::half_edge::involution::EdgeIndex>,
        decorated: &BTreeSet<usize>,
    ) -> Result<Self> {
        ensure!(
            graph.is_connected(support),
            "counterterm identity needs a connected support"
        );
        let boundary = graph.dummy_less_full_crown(support);
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
        let mut interfaces = vec![&numerator];
        if bare_skeleton {
            interfaces.extend([
                &graph.global_prefactor.num,
                &graph.global_prefactor.projector,
            ]);
        }
        for interface in interfaces {
            let tensor = SymbolicTensor::infer(cooking.try_cook(interface.as_view())?)?;
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
        for (vertex, _, data) in graph.underlying.iter_nodes_of(support) {
            let node = builder.graph.add_node(NodeColor::Vertex(data.dod.value));
            builder.vertices.insert(usize::from(vertex), node);
        }
        let all = support.union(&boundary);
        for (pair, edge, data) in graph.iter_edges_of(&all) {
            ensure!(
                !data.data.is_dummy,
                "dummy propagator in counterterm support"
            );
            let edge_id = usize::from(edge);
            let is_boundary = graph
                .iter_edges_of(&boundary)
                .any(|(_, candidate, _)| candidate == edge);
            let node = builder.graph.add_node(NodeColor::Edge {
                particle: data.data.particle().map(|particle| particle.index()),
                // The parent carrier is an amputated leg of this integral.
                // Its propagator denominator is not part of the counterterm.
                mass: if is_boundary {
                    Atom::Zero
                } else {
                    let mut mass = data.data.mass_atom(model);
                    if marked == Some(edge) {
                        mass += Atom::var(symbol!("ir_hard_identity_marked"));
                    }
                    if decorated.contains(&edge_id) {
                        mass += Atom::var(symbol!("ir_hard_identity_contracted"));
                    }
                    mass
                },
                power: if is_boundary {
                    0
                } else {
                    data.data.extra_data.vakint_edge_power.unwrap_or(1)
                },
                dod: if is_boundary { 0 } else { data.data.dod.value },
                boundary: is_boundary,
                orientation: match data.orientation {
                    _ if data.data.particle.is_self_antiparticle(model) => 0,
                    Orientation::Undirected => 0,
                    Orientation::Default => 1,
                    Orientation::Reversed => -1,
                },
            });
            builder.edges.insert(edge_id, node);
            let halfedges = match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    vec![source, sink]
                }
                HedgePair::Unpaired { hedge, .. } => vec![hedge],
            };
            for hedge in halfedges {
                if !all.includes(&hedge) {
                    continue;
                }
                let vertex = usize::from(graph.underlying.node_id(hedge));
                let Some(&vertex_node) = builder.vertices.get(&vertex) else {
                    continue;
                };
                let source = graph.underlying.flow(hedge) == Flow::Source;
                // Self-conjugate carrier directions are integration-variable
                // choices. Their signs are fixed after topology canonization.
                let flow = if is_boundary && fixed_boundary {
                    Some(source)
                } else if data.data.particle.is_self_antiparticle(model) {
                    None
                } else {
                    Some(source)
                };
                let hnode = builder.graph.add_node(NodeColor::HalfEdge {
                    boundary: is_boundary,
                    flow,
                });
                builder.hedges.insert(hedge.0, hnode);
                builder.link(node, hnode, Link::Incidence);
                builder.link(vertex_node, hnode, Link::Incidence);
                if is_boundary {
                    builder
                        .boundary_signs
                        .insert(edge_id, if source { -1 } else { 1 });
                }
            }
        }
        builder.topology_labels = builder.graph.canonize().vertex_map;
        builder.momentum_owners = builder
            .edges
            .iter()
            .map(|(edge, node)| (builder.topology_labels[*node], *edge))
            .collect();
        builder.momentum_signs = builder.boundary_signs.clone();
        for (pair, edge, data) in graph.iter_edges_of(support) {
            let sign = if let HedgePair::Paired { source, sink } = pair
                && data.data.particle.is_self_antiparticle(model)
                && builder.topology_labels[builder.hedges[&source.0]]
                    > builder.topology_labels[builder.hedges[&sink.0]]
            {
                -1
            } else {
                1
            };
            builder.momentum_signs.insert(usize::from(edge), sign);
            if let HedgePair::Paired { source, sink } = pair {
                builder.momentum_sources.insert(
                    usize::from(edge),
                    builder.hedges[&if sign == 1 { source.0 } else { sink.0 }],
                );
            }
        }
        // Preserve the location of each local rule. Moving an epsilon-dependent
        // coefficient across a subdivergence would change its local counterterm.
        for (vertex, _, data) in graph.underlying.iter_nodes_of(support) {
            let expression = builder.prepare(&data.num.value)?;
            let root = builder.expression(expression.as_view())?;
            builder.link(
                builder.vertices[&usize::from(vertex)],
                root,
                Link::Numerator,
            );
        }
        for (_, edge, data) in graph.iter_edges_of(support) {
            let expression = builder.prepare(&data.data.num.value)?;
            let root = builder.expression(expression.as_view())?;
            builder.link(builder.edges[&usize::from(edge)], root, Link::Numerator);
        }
        if bare_skeleton {
            for (color, expression) in [
                (NodeColor::GlobalPrefactor, &graph.global_prefactor.num),
                (
                    NodeColor::GlobalProjector,
                    &graph.global_prefactor.projector,
                ),
            ] {
                let expression = builder.prepare(expression)?;
                let root = builder.expression(expression.as_view())?;
                let owner = builder.graph.add_node(color);
                builder.link(owner, root, Link::Numerator);
            }
        }
        let canonical = builder.graph.canonize();
        let ports = builder
            .indices
            .iter()
            .filter(|(index, _)| builder.public.contains(index))
            .map(|(index, node)| (*index, canonical.vertex_map[*node]))
            .collect();
        let mut momenta = builder
            .boundary_signs
            .iter()
            .map(|(edge, sign)| (*edge, (canonical.vertex_map[builder.edges[edge]], *sign)))
            .collect::<BTreeMap<_, _>>();
        if bare_skeleton {
            for (pair, edge, data) in graph.iter_edges_of(support) {
                // The dressed rules can break an automorphism of the initial
                // topology. Orient transported internal momenta using the final
                // certificate rather than that preliminary topology labeling.
                let sign = if let HedgePair::Paired { source, sink } = pair
                    && data.data.particle.is_self_antiparticle(model)
                    && canonical.vertex_map[builder.hedges[&source.0]]
                        > canonical.vertex_map[builder.hedges[&sink.0]]
                {
                    -1
                } else {
                    1
                };
                let edge = usize::from(edge);
                momenta.insert(edge, (canonical.vertex_map[builder.edges[&edge]], sign));
            }
        }
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
        Ok(Self {
            key: CountertermKey {
                model: model.fingerprint(),
                graph: canonical.graph,
            },
            ports,
            momenta,
            coordinates,
            reservations,
        })
    }

    fn cooking() -> CookSettings {
        CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true)
    }

    pub(crate) fn canonicalize_delta(&self, value: &Atom) -> Result<Atom> {
        self.canonicalize_expression(value, *DELTA_MOMENTUM)
    }

    /// Transport a complete integrand into canonical edge labels while
    /// retaining the registered vector head required by typed tensor algebra.
    /// This comparison-only expression is not the opaque Delta-cache format.
    pub(crate) fn canonicalize_integrand(&self, value: &Atom) -> Result<Atom> {
        self.canonicalize_expression(value, GS.emr_mom)
    }

    fn canonicalize_expression(&self, value: &Atom, momentum_marker: Symbol) -> Result<Atom> {
        let cooking = Self::cooking();
        let private_marker = *DELTA_PRIVATE;
        let port_marker = *DELTA_PORT;
        let value = Self::edge_momenta(value, &self.coordinates)?;
        let scoped = SymbolicTensor::infer(cooking.try_cook(value.as_view())?)?
            .wrap_indices(private_marker, true)?;
        let raw = cooking.uncook(scoped.expression().as_view());
        let mut unknown = false;
        let ports = raw.replace_map(|term, _, out| {
            if let Ok(slot) = LibrarySlot::<Aind>::try_from(term) {
                if let Some(node) = self.ports.get(&(slot.rep().rep.base(), slot.aind)) {
                    **out = slot
                        .rep()
                        .to_symbolic([function!(port_marker, *node as i64)]);
                } else {
                    unknown = true;
                }
            }
        });
        ensure!(
            !unknown,
            "computed counterterm exposes an index absent from its input boundary"
        );
        let result = ports.replace_map(|term, _, out| {
            if let AtomView::Fun(fun) = term
                && fun.get_symbol() == GS.emr_mom
            {
                let mut args = fun.iter();
                let edge = args.next().and_then(|a| usize::try_from(a).ok());
                if let Some((canonical, sign)) = edge.and_then(|edge| self.momenta.get(&edge)) {
                    let mut function =
                        FunctionBuilder::new(momentum_marker).add_arg(*canonical as i64);
                    for arg in args {
                        function = function.add_arg(arg);
                    }
                    **out = function.finish() * *sign;
                } else {
                    unknown = true;
                }
            }
        });
        ensure!(
            !unknown,
            "computed counterterm retains a non-boundary momentum"
        );
        Ok(result)
    }

    pub(crate) fn restore_delta(&self, value: &Atom) -> Result<Atom> {
        let private_marker = *DELTA_PRIVATE;
        let port_marker = *DELTA_PORT;
        let momentum_marker = *DELTA_MOMENTUM;
        let scope_marker = AbstractIndex::scope_symbol();
        let ports = self
            .ports
            .iter()
            .map(|((_, index), node)| (*node, *index))
            .collect::<BTreeMap<_, _>>();
        let momenta = self
            .momenta
            .iter()
            .map(|(edge, (node, sign))| (*node, (*edge, *sign)))
            .collect::<BTreeMap<_, _>>();
        let state = ParseState::<Aind>::default();
        for reservation in &self.reservations {
            state.reserve_indices(reservation.as_view());
        }
        let mut unknown = false;
        let restored = value.replace_map(|term, _, out| {
            if let AtomView::Fun(fun) = term
                && fun.get_symbol() == momentum_marker
            {
                let mut args = fun.iter();
                let node = args.next().and_then(|a| usize::try_from(a).ok());
                if let Some((edge, sign)) = node.and_then(|node| momenta.get(&node)) {
                    let mut function = FunctionBuilder::new(GS.emr_mom).add_arg(*edge as i64);
                    for arg in args {
                        function = function.add_arg(arg);
                    }
                    **out = function.finish() * *sign;
                } else {
                    unknown = true;
                }
            }
        });
        let mut private = HashMap::new();
        let restored = restored.replace_map(|term, _, out| {
            if let AtomView::Fun(fun) = term {
                if fun.get_symbol() == port_marker {
                    let node = fun.iter().next().and_then(|a| usize::try_from(a).ok());
                    if let Some(index) = node.and_then(|node| ports.get(&node)) { **out = Atom::from(*index); }
                    else { unknown = true; }
                } else if fun.get_symbol() == scope_marker
                    && matches!(fun.iter().next(), Some(AtomView::Var(marker)) if marker.get_symbol() == private_marker) {
                    **out = private.entry(term.to_owned()).or_insert_with(|| Atom::from(state.fresh_index())).clone();
                }
            }
        });
        ensure!(
            !unknown,
            "counterterm cache boundary does not match its isomorphism certificate"
        );
        Ok(restored)
    }

    fn edge_momenta(value: &Atom, coordinates: &BTreeMap<(Symbol, usize), usize>) -> Result<Atom> {
        let mut unknown = false;
        let result = value.replace_map(|term, _, out| {
            if let AtomView::Fun(fun) = term
                && [
                    feynkit_graph::symbols::loop_momentum(),
                    feynkit_graph::symbols::external_momentum(),
                ]
                .contains(&fun.get_symbol())
            {
                let mut args = fun.iter();
                let coordinate = args.next().and_then(|arg| usize::try_from(arg).ok());
                if let Some(edge) = coordinate
                    .and_then(|coordinate| coordinates.get(&(fun.get_symbol(), coordinate)))
                {
                    let mut function = FunctionBuilder::new(GS.emr_mom).add_arg(*edge as i64);
                    for arg in args {
                        function = function.add_arg(arg);
                    }
                    **out = function.finish();
                } else {
                    unknown = true;
                }
            }
        });
        ensure!(
            !unknown,
            "counterterm expression references an unknown parent momentum coordinate"
        );
        Ok(result)
    }
}

/// Start of the comparison-only label range for indices owned inside a
/// counterterm support. Raw owners of other indices are graph element ids,
/// far below this bound, so the two ranges can never alias. The labels never
/// leave the identity builder: cached values are transported through ports.
const CANONICAL_OWNER_OFFSET: u16 = 1 << 14;

struct IdentityBuilder {
    graph: IncidenceGraph<NodeColor, Link>,
    edges: BTreeMap<usize, usize>,
    hedges: BTreeMap<usize, usize>,
    vertices: BTreeMap<usize, usize>,
    indices: BTreeMap<(LibraryRep, Aind), usize>,
    public: BTreeSet<(LibraryRep, Aind)>,
    boundary_signs: BTreeMap<usize, i64>,
    momentum_signs: BTreeMap<usize, i64>,
    momentum_owners: BTreeMap<usize, usize>,
    momentum_sources: BTreeMap<usize, usize>,
    topology_labels: Vec<usize>,
    index_owners: BTreeMap<(LibraryRep, Aind), Aind>,
    coordinates: BTreeMap<(Symbol, usize), usize>,
    cooking: CookSettings,
}

impl IdentityBuilder {
    fn link(&mut self, from: usize, to: usize, color: Link) {
        self.graph
            .add_edge(from, to, true, color)
            .expect("inserted incidence nodes");
    }

    fn prepare(&mut self, expression: &Atom) -> Result<Atom> {
        // Canonical temporary owner names precede Symbolica's antisymmetric
        // argument ordering. All extracted signs therefore participate in the
        // comparison of complete local rules, independently of input numbering.
        let expression = CountertermIdentity::edge_momenta(expression, &self.coordinates)?;
        let expression = self.cooking.try_cook(expression.as_view())?;
        let mut failure = None;
        let expression = expression.replace_map(|term, _, out| {
            if let Ok(slot) = LibrarySlot::<Aind>::try_from(term) {
                let raw = self.cooking.uncook(slot.aind.to_atom().as_view());
                let Ok(original) = Aind::try_from(raw.as_view()) else {
                    return;
                };
                let owner = match original {
                    Aind::Hedge(owner, _) => self.hedges.get(&usize::from(owner)),
                    Aind::Edge(owner, _) => self.edges.get(&usize::from(owner)),
                    Aind::Vertex(owner, _) => self.vertices.get(&usize::from(owner)),
                    _ => None,
                }
                .copied();
                let normalized = if let Some(owner) = owner {
                    // Canonical owners live in their own label range. Indices
                    // owned outside the support keep their raw graph owner,
                    // which may numerically equal a canonical label.
                    let Some(owner) = u16::try_from(self.topology_labels[owner])
                        .ok()
                        .and_then(|label| label.checked_add(CANONICAL_OWNER_OFFSET))
                    else {
                        failure = Some("canonical graph owner exceeds runtime index capacity");
                        return;
                    };
                    match original {
                        Aind::Hedge(_, local) => Aind::Hedge(owner, local),
                        Aind::Edge(_, local) => Aind::Edge(owner, local),
                        Aind::Vertex(_, local) => Aind::Vertex(owner, local),
                        _ => unreachable!(),
                    }
                } else {
                    if let Aind::Hedge(raw, _) | Aind::Edge(raw, _) | Aind::Vertex(raw, _) = original
                        && raw >= CANONICAL_OWNER_OFFSET
                    {
                        failure = Some("an index outside the support exceeds the raw owner range");
                        return;
                    }
                    original
                };
                if let Some(previous) = self
                    .index_owners
                    .insert((slot.rep().rep.base(), normalized), original)
                    && previous != original
                {
                    failure = Some("an index outside the support aliases a canonical graph owner");
                    return;
                }
                **out = slot.rep().to_symbolic([Atom::from(normalized)]);
            }
        });
        ensure!(failure.is_none(), "{}", failure.unwrap_or_default());
        let expression = expression.replace_map(|term, _, out| {
            if let AtomView::Fun(fun) = term
                && fun.get_symbol() == GS.emr_mom
            {
                let mut args = fun.iter();
                if let Some(edge) = args.next().and_then(|a| usize::try_from(a).ok()) {
                    if let Some(node) = self.edges.get(&edge) {
                        let mut function = FunctionBuilder::new(GS.emr_mom)
                            .add_arg(self.topology_labels[*node] as i64);
                        for arg in args {
                            function = function.add_arg(arg);
                        }
                        **out = function.finish() * self.momentum_signs[&edge];
                    } else {
                        failure = Some(
                            "numerator momentum lies outside counterterm support and boundary",
                        );
                    }
                }
            }
        });
        ensure!(failure.is_none(), "{}", failure.unwrap_or_default());
        // Reorienting Q can move a minus sign through a sum. Expand within
        // this rule only; epsilon-dependent factors never cross UV vertices.
        Ok(expression.expand())
    }

    fn index(&mut self, representation: LibraryRep, index: Aind) -> usize {
        let key = (representation.base(), index);
        if let Some(node) = self.indices.get(&key) {
            return *node;
        }
        let public = self.public.contains(&key);
        let node = self.graph.add_node(NodeColor::Index { public });
        self.indices.insert(key, node);
        if public {
            let owner = match index {
                Aind::Hedge(owner, _) => self.hedges.get(&usize::from(owner)),
                Aind::Edge(owner, _) => self.edges.get(&usize::from(owner)),
                Aind::Vertex(owner, _) => self.vertices.get(&usize::from(owner)),
                _ => None,
            }
            .copied();
            if let Some(owner) = owner {
                self.link(node, owner, Link::Owner);
            }
        }
        node
    }

    fn expression(&mut self, expression: AtomView<'_>) -> Result<usize> {
        if let Ok(slot) = LibrarySlot::<Aind>::try_from(expression) {
            let raw = self.cooking.uncook(slot.aind.to_atom().as_view());
            let normalized = Aind::try_from(raw.as_view())?;
            let original = self
                .index_owners
                .get(&(slot.rep().rep.base(), normalized))
                .copied()
                .unwrap_or(normalized);
            let index = self.index(slot.rep().rep, original);
            let node = self
                .graph
                .add_node(NodeColor::Slot(slot.rep().to_symbolic([])));
            self.link(node, index, Link::Index);
            return Ok(node);
        }
        let (color, args, ordered) = match expression {
            AtomView::Fun(fun) if fun.get_symbol() == GS.emr_mom => {
                let mut args = fun.iter();
                let edge = usize::try_from(args.next().ok_or_else(|| eyre!("missing Q edge"))?)
                    .map_err(|error| eyre!(error))?;
                let original = self
                    .momentum_owners
                    .get(&edge)
                    .copied()
                    .ok_or_else(|| eyre!("unknown canonical momentum owner"))?;
                let owner = *self.edges.get(&original).ok_or_else(|| {
                    eyre!("numerator momentum lies outside counterterm support and boundary")
                })?;
                let node = self.graph.add_node(NodeColor::Momentum);
                self.link(node, owner, Link::Momentum);
                if let Some(source) = self.momentum_sources.get(&original).copied() {
                    self.link(node, source, Link::MomentumSource);
                }
                for (position, arg) in args.enumerate() {
                    let child = self.expression(arg)?;
                    self.link(node, child, Link::Argument(position));
                }
                return Ok(node);
            }
            AtomView::Fun(fun) => (
                NodeColor::Function(fun.get_symbol()),
                fun.iter().collect::<Vec<_>>(),
                !fun.get_symbol().is_symmetric(),
            ),
            AtomView::Add(sum) => (NodeColor::Sum, sum.iter().collect(), false),
            AtomView::Mul(product) => (NodeColor::Product, product.iter().collect(), false),
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                (NodeColor::Power, vec![base, exponent], true)
            }
            _ => return Ok(self.graph.add_node(NodeColor::Leaf(expression.to_owned()))),
        };
        let node = self.graph.add_node(color);
        for (position, arg) in args.into_iter().enumerate() {
            let child = self.expression(arg)?;
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
    use linnet::half_edge::involution::EdgeIndex;
    use spenso::{
        network::{library::symbolic::ETS, tags::SPENSO_TAG},
        structure::representation::{Minkowski, RepName},
    };

    fn bubble(model: &Model, permuted: bool) -> Result<Graph> {
        let mut graph: Graph = if permuted {
            finalized_runtime_dot!(
                digraph G {
                    graph [projector=19]; node [num=1]; edge [particle="scalar_1" num=1];
                    e [style=invis]; B; A;
                    A:0 -> e [id=0 source="{ufo_order:0}"];
                    e -> B:1 [id=1 sink="{ufo_order:0}"];
                    A -> B [id=3 source="{ufo_order:1}" sink="{ufo_order:1}"];
                    A -> B [id=2 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
                }, model
            )?
        } else {
            finalized_runtime_dot!(
                digraph G {
                    graph [projector=7]; node [num=1]; edge [particle="scalar_1" num=1];
                    e [style=invis]; A; B;
                    e -> A:0 [id=2 sink="{ufo_order:0}"];
                    B:1 -> e [id=3 source="{ufo_order:0}"];
                    A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
                    A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
                }, model
            )?
        };
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            let data = &mut graph.underlying[edge];
            data.particle = data.particle.clone().override_mass(Some(Atom::Zero));
        }
        Ok(graph)
    }

    fn support(graph: &Graph) -> SuBitGraph {
        graph
            .no_dummy()
            .subtract(&graph.external_filter::<SuBitGraph>())
    }

    fn decorate(graph: &mut Graph, external: [usize; 2], dummy: usize) -> [Aind; 2] {
        let mink = Minkowski {}.new_rep(GS.dim);
        let mut ports = Vec::new();
        for (position, edge) in external.into_iter().enumerate() {
            let hedge = graph
                .iter_edges()
                .find_map(|(pair, id, _)| {
                    if usize::from(id) != edge {
                        return None;
                    }
                    match pair {
                        HedgePair::Unpaired { hedge, .. } => Some(hedge),
                        _ => None,
                    }
                })
                .unwrap();
            let port = Aind::Hedge(hedge.0 as u16, 1);
            let vertex = graph.underlying.node_id(hedge);
            let metric = function!(
                ETS.metric,
                mink.to_symbolic([Atom::from(port)]),
                mink.to_symbolic([Atom::from(Aind::Dummy(dummy))])
            );
            graph.underlying[vertex].num = Autogen::explicit(metric * (position as i64 + 2));
            ports.push(port);
        }
        [ports[0], ports[1]]
    }

    #[test]
    fn isomorphism_ignores_graph_labels_and_private_indices_and_transports_open_tensor()
    -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut first = bubble(&model, false)?;
        let mut second = bubble(&model, true)?;
        let first_ports = decorate(&mut first, [2, 3], 400);
        let second_ports = decorate(&mut second, [0, 1], 900);
        let first = CountertermIdentity::new(&first, &support(&first), &model)?;
        let second = CountertermIdentity::new(&second, &support(&second), &model)?;
        assert_eq!(first.key, second.key);
        let mink = Minkowski {}.new_rep(GS.dim);
        let tensor = |edge: usize, ports: [Aind; 2]| {
            let q = function!(GS.emr_mom, edge as i64, mink.to_symbolic([]));
            function!(
                GS.emr_mom,
                edge as i64,
                mink.to_symbolic([Atom::from(ports[0])])
            ) * function!(
                GS.emr_mom,
                edge as i64,
                mink.to_symbolic([Atom::from(ports[1])])
            ) * function!(Symbol::LOG, function!(SPENSO_TAG.dot, &q, &q))
                / Atom::var(GS.dim_epsilon)
        };
        let stored = first.canonicalize_delta(&tensor(2, first_ports))?;
        let restored = second.restore_delta(&stored)?;
        assert_eq!(restored, tensor(0, second_ports));
        let vector = function!(
            GS.emr_mom,
            2,
            mink.to_symbolic([Atom::from(first_ports[0])])
        );
        let stored = first.canonicalize_delta(&vector)?;
        assert_eq!(
            second.restore_delta(&stored)?,
            -function!(
                GS.emr_mom,
                0,
                mink.to_symbolic([Atom::from(second_ports[0])])
            )
        );
        Ok(())
    }

    #[test]
    fn isomorphism_and_transport_distinguish_representations_not_dummy_names() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let mut first = bubble(&model, false)?;
        let mut second = bubble(&model, true)?;
        let first_lorentz = decorate(&mut first, [2, 3], 400);
        let second_lorentz = decorate(&mut second, [0, 1], 900);
        let color = idenso::representations::ColorAdjoint {}.new_rep(8);
        let decorate_color = |graph: &mut Graph, ports: [Aind; 2], local: u16, dummy: usize| {
            ports.map(|port| {
                let Aind::Hedge(owner, _) = port else {
                    unreachable!()
                };
                let port = Aind::Hedge(owner, local);
                let vertex = graph
                    .underlying
                    .node_id(linnet::half_edge::involution::Hedge(usize::from(owner)));
                graph.underlying[vertex].num.value *= function!(
                    ETS.metric,
                    color.to_symbolic([Atom::from(port)]),
                    color.to_symbolic([Atom::from(Aind::Dummy(dummy))])
                );
                port
            })
        };
        // The first occurrence reuses the same names in unrelated index
        // representations. The second renames color independently of Lorentz.
        let first_color = decorate_color(&mut first, first_lorentz, 1, 400);
        let second_color = decorate_color(&mut second, second_lorentz, 7, 901);
        let first = CountertermIdentity::new(&first, &support(&first), &model)?;
        let second = CountertermIdentity::new(&second, &support(&second), &model)?;
        assert_eq!(first.key, second.key);
        let mink = Minkowski {}.new_rep(GS.dim);
        let delta = |lorentz: [Aind; 2], colors: [Aind; 2]| {
            function!(
                ETS.metric,
                mink.to_symbolic([Atom::from(lorentz[0])]),
                mink.to_symbolic([Atom::from(lorentz[1])])
            ) * function!(
                ETS.metric,
                color.to_symbolic([Atom::from(colors[0])]),
                color.to_symbolic([Atom::from(colors[1])])
            )
        };
        let stored = first.canonicalize_delta(&delta(first_lorentz, first_color))?;
        assert_eq!(
            second.restore_delta(&stored)?,
            delta(second_lorentz, second_color)
        );
        Ok(())
    }

    #[test]
    fn isomorphism_preserves_particle_mass_power_and_local_numerator() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let graph = bubble(&model, false)?;
        let key = CountertermIdentity::new(&graph, &support(&graph), &model)?.key;
        let mut changed = graph.clone();
        changed.underlying[EdgeIndex::from(0)]
            .extra_data
            .vakint_edge_power = Some(2);
        assert_ne!(
            key,
            CountertermIdentity::new(&changed, &support(&changed), &model)?.key
        );
        let mut changed = graph.clone();
        changed.underlying[EdgeIndex::from(0)].particle = changed.underlying[EdgeIndex::from(0)]
            .particle
            .clone()
            .override_mass(Some(Atom::num(2)));
        assert_ne!(
            key,
            CountertermIdentity::new(&changed, &support(&changed), &model)?.key
        );
        let mut changed = graph.clone();
        changed.underlying[EdgeIndex::from(0)].particle =
            crate::graph::edge::PossibleParticle::JustMass { expr: Atom::Zero };
        assert_ne!(
            key,
            CountertermIdentity::new(&changed, &support(&changed), &model)?.key
        );
        let mut changed = graph.clone();
        changed.underlying[EdgeIndex::from(0)].num =
            Autogen::explicit(Atom::one() + Atom::var(GS.dim_epsilon));
        assert_ne!(
            key,
            CountertermIdentity::new(&changed, &support(&changed), &model)?.key
        );
        Ok(())
    }

    #[test]
    fn cached_private_contractions_are_fresh_on_every_use() -> Result<()> {
        test_initialise()?;
        let model = load_generic_model("scalars");
        let graph = bubble(&model, false)?;
        let identity = CountertermIdentity::new(&graph, &support(&graph), &model)?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let index = mink.to_symbolic([Atom::from(Aind::Dummy(99001))]);
        let contracted = function!(GS.emr_mom, 2, &index).pow(2);
        let stored = identity.canonicalize_delta(&contracted)?;
        let first = identity.restore_delta(&stored)?;
        let second = identity.restore_delta(&stored)?;
        assert_ne!(
            first, second,
            "separate insertions must not contract with each other"
        );
        assert!(!first.contains(Atom::from(Aind::Dummy(99001))));
        assert!(!second.contains(Atom::from(Aind::Dummy(99001))));
        Ok(())
    }
}
