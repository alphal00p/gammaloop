//! Diagram enumeration with Linnet's native generator.
//!
//! Every allowed vertex rule is one generator signature, so rules with identical field content
//! (split UFO vertices, counterterms) give distinct diagrams, each with the automorphism count of
//! its rule-coloured graph. External legs are labelled by [`ExternalLabel`]: legs with equal label
//! and stub are interchangeable, which implements external-state symmetrization.
//!
//! The process filters run inside the search through the generator's filter, at the earliest
//! stage where they are decided: bridges and special topologies on untyped topologies,
//! coupling-order upper bounds on partial field assignments, and the remaining conditions on
//! complete graphs.

use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use feynkit_graph::ExternalState;
use feynkit_model::{Model, Particle, ParticleId, VertexRuleId};
use linnet::half_edge::{
    NodeIndex,
    generate::{Diagram, GraphView, Stage, Stub},
    involution::{Flow, HedgePair, Orientation},
};

use super::{
    GenerationError,
    graph::{Endpoint, ExternalNode, GeneratedEdge, GeneratedGraph},
    special_topologies::{PlainGraph, SpecialTopologyVetoes},
};
use crate::{GenerationFilter, GenerationOptions, ParticleSelector, VertexSelector};

/// The label of an external leg during generation. Legs with equal labels and stubs are
/// interchangeable: symmetrized legs of one state share a label, all others are exact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ExternalLabel {
    Exact { state: ExternalState, index: usize },
    WithinState(ExternalState),
}

/// A search state of diagram generation, as seen by [`GenerationOptions::filter`]: edge fields
/// are the particle of each particle/antiparticle pair, directed along its flow, and node rules
/// index the generation signatures.
pub type GenerationView<'a> = GraphView<'a, ParticleId, ExternalLabel>;

/// A generator stub: the particle of its conjugate pair as field, flowing into the vertex for a
/// particle and out of it for an antiparticle. Self-conjugate particles are undirected.
pub(crate) type ParticleStub = Stub<ParticleId>;

/// The stub of a particle entering a vertex, as the particles of a UFO vertex do.
pub(crate) fn incoming_stub(model: &Model, particle: &Particle) -> ParticleStub {
    let id = model.particle_id(&particle.name).expect("a model particle");
    if particle.antiparticle == id {
        Stub::undirected(id)
    } else if particle.is_antiparticle() {
        Stub::outgoing(particle.antiparticle)
    } else {
        Stub::incoming(id)
    }
}

/// The vertex rules available to a generation request, one generator signature each.
pub(crate) struct VertexSignatures {
    pub(crate) rules: Vec<VertexRuleId>,
    pub(crate) stubs: Vec<Vec<ParticleStub>>,
}

impl VertexSignatures {
    /// The rules allowed by the vertex selections, the particle vetoes and the coupling-order
    /// upper bounds: coupling powers are nonnegative, so a vertex exceeding a graph-wide upper
    /// bound cannot occur in an accepted graph.
    pub(crate) fn new(model: &Model, options: &GenerationOptions) -> Result<Self, GenerationError> {
        let vertex_selection = |veto: bool| {
            {
                options
                    .graph_filters
                    .iter()
                    .find_map(|filter| match filter {
                        GenerationFilter::VertexAllow(selectors) if !veto => Some(selectors),
                        GenerationFilter::VertexVeto(selectors) if veto => Some(selectors),
                        _ => None,
                    })
            }
            .map(|selectors| {
                selectors
                    .iter()
                    .map(resolved_vertex)
                    .collect::<BTreeSet<_>>()
            })
        };
        let (allowed, vetoed) = (vertex_selection(false), vertex_selection(true));
        let particle_vetoes = particle_vetoes(options);
        let bounds = CouplingBounds::new(options);
        let mut signatures = Self {
            rules: vec![],
            stubs: vec![],
        };
        for (rule_index, rule) in model.vertex_rules().iter().enumerate() {
            let rule_id = model.vertex_rule_id_at(rule_index)?;
            if allowed
                .as_ref()
                .is_some_and(|allowed| !allowed.contains(&rule_id))
                || vetoed
                    .as_ref()
                    .is_some_and(|vetoed| vetoed.contains(&rule_id))
            {
                continue;
            }
            if bounds.exceeded_by(&rule.coupling_orders(model)) {
                continue;
            }
            let particles = rule
                .particles
                .iter()
                .map(|particle| model.particle_by_id(*particle))
                .collect::<Result<Vec<_>, _>>()?;
            if particles.iter().any(|particle| {
                particle_vetoes.as_ref().is_some_and(|vetoed| {
                    vetoed.contains(&model.particle_id(&particle.name).expect("a model particle"))
                        || vetoed.contains(&particle.antiparticle)
                })
            }) {
                continue;
            }
            signatures.rules.push(rule_id);
            signatures.stubs.push(
                particles
                    .into_iter()
                    .map(|particle| incoming_stub(model, particle))
                    .collect(),
            );
        }
        Ok(signatures)
    }
}

pub(crate) fn resolved_vertex(selector: &VertexSelector) -> VertexRuleId {
    match selector {
        VertexSelector::Id { vertex, .. } => *vertex,
        VertexSelector::Name(_) => unreachable!("generation selectors are resolved before use"),
    }
}

fn particle_vetoes(options: &GenerationOptions) -> Option<BTreeSet<ParticleId>> {
    options
        .graph_filters
        .iter()
        .find_map(|filter| match filter {
            GenerationFilter::ParticleVeto(selectors) => Some(
                selectors
                    .iter()
                    .map(|selector| match selector {
                        ParticleSelector::Id { particle, .. } => *particle,
                        ParticleSelector::Name(_) | ParticleSelector::Pdg(_) => {
                            unreachable!("generation selectors are resolved before use")
                        }
                    })
                    .collect(),
            ),
            _ => None,
        })
}

/// The graph-wide coupling-order bounds `(name, minimum, maximum)`.
struct CouplingBounds(Vec<(String, usize, Option<usize>)>);

impl CouplingBounds {
    fn new(options: &GenerationOptions) -> Self {
        Self(
            options
                .graph_filters
                .iter()
                .filter_map(|filter| match filter {
                    GenerationFilter::CouplingOrders(required) => Some(required),
                    _ => None,
                })
                .flatten()
                .map(|(name, (minimum, maximum))| (name.clone(), *minimum, *maximum))
                .collect(),
        )
    }

    fn exceeded_by(&self, orders: &BTreeMap<String, usize>) -> bool {
        self.0.iter().any(|(name, _, maximum)| {
            maximum.is_some_and(|maximum| orders.get(name).copied().unwrap_or(0) > maximum)
        })
    }
}

/// The external legs of a process: the generator labels and stubs, and the momentum labels
/// that [`GeneratedGraph`]s assign to them.
pub(crate) struct ProcessExternals {
    legs: Vec<(ExternalLabel, ParticleStub)>,
    /// Interchangeable legs and their momentum labels in increasing index order.
    classes: BTreeMap<(ExternalLabel, ParticleStub), Vec<ExternalNode>>,
    /// The non-identity permutations of interchangeable sewn pairs, as maps of exact leg labels.
    sewn_pair_permutations: Vec<BTreeMap<ExternalLabel, ExternalLabel>>,
}

impl ProcessExternals {
    /// `sewn_pairs` are classes of interchangeable `[incoming, outgoing]` sewn pairs: a
    /// permutation within a class moves both legs of every pair at once.
    pub(crate) fn new(
        externals: impl IntoIterator<Item = (ExternalLabel, ParticleStub, ExternalNode)>,
        sewn_pairs: impl IntoIterator<Item = Vec<[ExternalLabel; 2]>>,
    ) -> Self {
        let mut legs = vec![];
        let mut classes: BTreeMap<_, Vec<_>> = BTreeMap::new();
        for (label, stub, external) in externals {
            legs.push((label, stub.clone()));
            classes.entry((label, stub)).or_default().push(external);
        }
        for externals in classes.values_mut() {
            externals.sort();
        }
        let mut sewn_pair_permutations = vec![BTreeMap::new()];
        for pairs in sewn_pairs {
            let pairs = &pairs;
            sewn_pair_permutations = sewn_pair_permutations
                .into_iter()
                .flat_map(|permutation| {
                    permutations(pairs.len()).map(move |images| {
                        let mut permutation = permutation.clone();
                        for (pair, image) in pairs.iter().zip(images) {
                            permutation.insert(pair[0], pairs[image][0]);
                            permutation.insert(pair[1], pairs[image][1]);
                        }
                        permutation
                    })
                })
                .collect();
        }
        sewn_pair_permutations
            .retain(|permutation| permutation.iter().any(|(from, to)| from != to));
        Self {
            legs,
            classes,
            sewn_pair_permutations,
        }
    }

    /// The number of sewn-pair permutations that fix the isomorphism class of `diagram`, or
    /// `None` when one maps it to a smaller certificate.
    ///
    /// The generator labels sewn pairs exactly, so a class of diagrams related by permuting
    /// interchangeable pairs is generated once per member. The member with the smallest
    /// certificate represents the class, and its automorphism count is multiplied by the
    /// returned stabilizer order to give the automorphisms of the sewn graph.
    pub(crate) fn sewn_pair_stabilizer(
        &self,
        diagram: &Diagram<ParticleId, ExternalLabel>,
    ) -> Option<u128> {
        let mut stabilizer = 1;
        for permutation in &self.sewn_pair_permutations {
            let image = diagram
                .graph
                .canonize_by(
                    |&rule| rule,
                    |field| *field,
                    |label| label.map(|label| permutation.get(&label).copied().unwrap_or(label)),
                )
                .into_certificate();
            match image.cmp(&diagram.certificate) {
                Ordering::Less => return None,
                Ordering::Equal => stabilizer += 1,
                Ordering::Greater => {}
            }
        }
        Some(stabilizer)
    }

    pub(crate) fn legs(&self) -> impl Iterator<Item = (ExternalLabel, ParticleStub)> + '_ {
        self.legs.iter().cloned()
    }

    /// Whether a particle veto excludes an external leg, which excludes every graph.
    pub(crate) fn vetoed(&self, model: &Model, options: &GenerationOptions) -> bool {
        let vetoes = particle_vetoes(options).unwrap_or_default();
        self.legs.iter().any(|(_, stub)| {
            vetoes.contains(&stub.field)
                || model
                    .particle_by_id(stub.field)
                    .is_ok_and(|particle| vetoes.contains(&particle.antiparticle))
        })
    }

    /// The generated diagram with its vertex rules and external momenta. Interchangeable legs
    /// take their momentum labels in canonical half-edge order, so isomorphic diagrams are
    /// labelled alike.
    pub(crate) fn graph(
        &self,
        diagram: &Diagram<ParticleId, ExternalLabel>,
        signatures: &VertexSignatures,
    ) -> GeneratedGraph {
        let graph = &diagram.graph;
        let canonical = self
            .classes
            .values()
            .any(|externals| externals.len() > 1)
            .then(|| graph.canonize_by(|&rule| rule, |field| *field, |label| *label));
        let mut legs = vec![];
        let mut edges = vec![];
        for (pair, _, data) in graph.iter_edges() {
            let directed = data.orientation != Orientation::Undirected;
            match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    edges.push(Some(GeneratedEdge {
                        source: Endpoint::Vertex(graph.node_id(source).0),
                        sink: Endpoint::Vertex(graph.node_id(sink).0),
                        particle: *data.data,
                        directed,
                    }));
                }
                HedgePair::Unpaired { hedge, flow } => {
                    let label = graph[hedge].expect("legs are labelled");
                    let stub = Stub {
                        field: *data.data,
                        flow: directed.then_some(flow),
                    };
                    let order = canonical
                        .as_ref()
                        .map_or(hedge, |canonical| canonical.canonical_hedge(hedge));
                    legs.push(((label, stub), order, edges.len(), hedge, flow));
                    edges.push(None);
                }
            }
        }
        legs.sort_by(|left, right| (&left.0, left.1).cmp(&(&right.0, right.1)));
        let mut next = BTreeMap::<&(ExternalLabel, ParticleStub), usize>::new();
        for (class, _, position, hedge, flow) in &legs {
            let slot = next.entry(class).or_default();
            let external = Endpoint::External(self.classes[class][*slot]);
            *slot += 1;
            let vertex = Endpoint::Vertex(graph.node_id(*hedge).0);
            let (source, sink) = match flow {
                Flow::Sink => (external, vertex),
                Flow::Source => (vertex, external),
            };
            edges[*position] = Some(GeneratedEdge {
                source,
                sink,
                particle: class.1.field,
                directed: class.1.flow.is_some(),
            });
        }
        GeneratedGraph::new(
            (0..graph.n_nodes()).map(|node| signatures.rules[graph[NodeIndex(node)]]),
            edges
                .into_iter()
                .map(|edge| edge.expect("every edge is placed")),
        )
    }
}

type UserFilter = dyn Fn(&GenerationView<'_>) -> bool + Send + Sync;

/// The process filters evaluated during the search; see the module documentation.
pub(crate) struct Pruner {
    one_particle_irreducible_loops: bool,
    vetoes: SpecialTopologyVetoes,
    coupling_bounds: Vec<(usize, Option<usize>)>,
    /// The coupling orders of each signature, in the order of `coupling_bounds`.
    signature_orders: Vec<Vec<usize>>,
    fermion_loop_range: Option<(usize, usize)>,
    /// Per particle index: massless, and a non-ghost fermion.
    massless: Vec<bool>,
    fermion: Vec<bool>,
    user: Option<Arc<UserFilter>>,
}

impl Pruner {
    pub(crate) fn new(
        model: &Model,
        options: &GenerationOptions,
        signatures: &VertexSignatures,
    ) -> Result<Self, GenerationError> {
        let vetoes = SpecialTopologyVetoes {
            self_energy: options
                .graph_filters
                .iter()
                .find_map(|filter| match filter {
                    GenerationFilter::SelfEnergy(options) => Some(*options),
                    _ => None,
                }),
            tadpoles: options
                .graph_filters
                .iter()
                .find_map(|filter| match filter {
                    GenerationFilter::Tadpoles(options) => Some(*options),
                    _ => None,
                }),
            snails: options
                .graph_filters
                .iter()
                .find_map(|filter| match filter {
                    GenerationFilter::ZeroSnails(options) => Some(*options),
                    _ => None,
                }),
            factorized_loop_range: options
                .graph_filters
                .iter()
                .find_map(|filter| match filter {
                    GenerationFilter::FactorizedLoopTopologiesCountRange(range) => Some(*range),
                    _ => None,
                }),
        };
        let bounds = CouplingBounds::new(options);
        let signature_orders = signatures
            .rules
            .iter()
            .map(|rule| {
                let orders = model.vertex_rule_by_id(*rule)?.coupling_orders(model);
                Ok(bounds
                    .0
                    .iter()
                    .map(|(name, _, _)| orders.get(name).copied().unwrap_or(0))
                    .collect())
            })
            .collect::<Result<_, GenerationError>>()?;
        Ok(Self {
            one_particle_irreducible_loops: options
                .graph_filters
                .contains(&GenerationFilter::LoopOneParticleIrreducible),
            vetoes,
            coupling_bounds: bounds
                .0
                .into_iter()
                .map(|(_, minimum, maximum)| (minimum, maximum))
                .collect(),
            signature_orders,
            fermion_loop_range: options
                .graph_filters
                .iter()
                .find_map(|filter| match filter {
                    GenerationFilter::FermionLoopCountRange(range) => Some(*range),
                    _ => None,
                }),
            massless: (0..model.particles().len())
                .map(|index| model.particle_is_massless(ParticleId::from_index(index)))
                .collect(),
            fermion: model
                .particles()
                .iter()
                .map(|particle| particle.is_fermion() && !particle.is_ghost())
                .collect(),
            user: options.filter.clone(),
        })
    }

    pub(crate) fn keeps(&self, view: &GenerationView<'_>) -> bool {
        let keeps = match view.stage() {
            Stage::Topology => {
                !(self.one_particle_irreducible_loops
                    && view.loop_order() > 0
                    && !view.bridges().is_empty())
                    && !(view.is_topology_complete() && self.special_topology_vetoed(view))
            }
            Stage::Fields => self.coupling_orders_allowed(view, false),
            Stage::Complete => {
                self.coupling_orders_allowed(view, true)
                    && !self.special_topology_vetoed(view)
                    && self.fermion_loop_range.is_none_or(|(minimum, maximum)| {
                        self.closed_fermion_loops(view)
                            .is_none_or(|loops| (minimum..=maximum).contains(&loops))
                    })
            }
        };
        keeps && self.user.as_ref().is_none_or(|user| user(view))
    }

    /// Partial field assignments are bounded from below by the cheapest rule of every vertex
    /// whose stubs are typed; complete graphs are checked exactly.
    fn coupling_orders_allowed(&self, view: &GenerationView<'_>, complete: bool) -> bool {
        if self.coupling_bounds.is_empty() {
            return true;
        }
        let mut totals = vec![0; self.coupling_bounds.len()];
        for vertex in 0..view.n_vertices() {
            let Some(rules) = view.rules(vertex) else {
                continue;
            };
            for (order, total) in totals.iter_mut().enumerate() {
                *total += rules
                    .iter()
                    .map(|rule| self.signature_orders[*rule as usize][order])
                    .min()
                    .unwrap_or(0);
            }
        }
        self.coupling_bounds
            .iter()
            .zip(totals)
            .all(|((minimum, maximum), total)| {
                maximum.is_none_or(|maximum| total <= maximum) && (!complete || total >= *minimum)
            })
    }

    fn special_topology_vetoed(&self, view: &GenerationView<'_>) -> bool {
        !self.vetoes.is_empty() && self.vetoes.vetoes(&self.plain(view))
    }

    /// The view as a plain graph with masses where the fields are known.
    fn plain(&self, view: &GenerationView<'_>) -> PlainGraph {
        let mut graph = PlainGraph::default();
        for _ in 0..view.n_vertices() {
            graph.add_vertex();
        }
        for vertex in 0..view.n_vertices() {
            for (_, stub) in view.legs(vertex) {
                graph.add_leg(vertex, self.massless[stub.field.index()]);
            }
        }
        for edge in view.edges() {
            graph.add_edge(
                (edge.vertices[0], edge.vertices[1]),
                edge.stub.map(|stub| self.massless[stub.field.index()]),
            );
        }
        graph
    }

    /// Closed fermion loops of a complete graph; `None` for branching fermion lines, which
    /// diagram construction rejects.
    fn closed_fermion_loops(&self, view: &GenerationView<'_>) -> Option<usize> {
        let mut nodes = view.n_vertices();
        let mut ends = vec![];
        for vertex in 0..view.n_vertices() {
            for (_, stub) in view.legs(vertex) {
                if self.fermion[stub.field.index()] {
                    ends.push((vertex, nodes));
                    nodes += 1;
                }
            }
        }
        for edge in view.edges() {
            let stub = edge.stub.expect("complete graphs are typed");
            if self.fermion[stub.field.index()] {
                ends.push((edge.vertices[0], edge.vertices[1]));
            }
        }
        closed_loops(nodes, &ends).ok()
    }
}

/// The number of closed components of the edges `ends` between `nodes` nodes, or the first node
/// with a degree above two. External legs are degree-one nodes, so components through them are
/// open lines; a self-loop gives its node degree two.
pub(crate) fn closed_loops(nodes: usize, ends: &[(usize, usize)]) -> Result<usize, (usize, usize)> {
    let mut adjacency = vec![vec![]; nodes];
    let mut degrees = vec![0; nodes];
    for &(left, right) in ends {
        degrees[left] += 1;
        degrees[right] += 1;
        adjacency[left].push(right);
        if left != right {
            adjacency[right].push(left);
        }
    }
    if let Some((node, degree)) = degrees
        .iter()
        .copied()
        .enumerate()
        .find(|(_, degree)| *degree > 2)
    {
        return Err((node, degree));
    }
    let mut visited = vec![false; nodes];
    let mut closed_loops = 0;
    for start in 0..nodes {
        if degrees[start] == 0 || visited[start] {
            continue;
        }
        let mut stack = vec![start];
        let mut closed = true;
        while let Some(node) = stack.pop() {
            if visited[node] {
                continue;
            }
            visited[node] = true;
            closed &= degrees[node] == 2;
            stack.extend(
                adjacency[node]
                    .iter()
                    .copied()
                    .filter(|next| !visited[*next]),
            );
        }
        closed_loops += usize::from(closed);
    }
    Ok(closed_loops)
}

/// The permutations of `0..n`, as the image of each element.
fn permutations(n: usize) -> impl Iterator<Item = Vec<usize>> {
    let mut all = vec![vec![]];
    for length in 0..n {
        all = all
            .into_iter()
            .flat_map(|permutation: Vec<usize>| {
                (0..=length).map(move |position| {
                    let mut permutation = permutation.clone();
                    permutation.insert(position, length);
                    permutation
                })
            })
            .collect();
    }
    all.into_iter()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn closed_loops_count_components_without_leg_nodes() {
        // a self-loop, a two-edge cycle, and an open line ending on leg nodes 4 and 5
        let ends = [(0, 0), (1, 2), (2, 1), (3, 4), (3, 5)];
        assert_eq!(closed_loops(6, &ends), Ok(2));
        // two legs at one vertex open its line even without other edges
        assert_eq!(closed_loops(3, &[(0, 1), (0, 2)]), Ok(0));
        assert_eq!(closed_loops(4, &[(0, 1), (0, 2), (0, 3)]), Err((0, 3)));
    }
}
