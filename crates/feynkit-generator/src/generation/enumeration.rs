//! Diagram enumeration with Linnet's native generator.
//!
//! Every allowed vertex rule is one generator signature, so rules with identical field content
//! (split UFO vertices, counterterms) give distinct diagrams, each with the automorphism count of
//! its rule-coloured graph. External legs are labelled by [`ExternalLabel`]: legs with equal label
//! and stub are interchangeable, which implements external-state symmetrization.
//!
//! Coupling orders and the fermion-loop range are budgets of the generator: it bounds them on
//! every partial field assignment, drops the assignments without a completion in range and checks
//! them exactly on the vertex rules (see
//! [`linnet::half_edge::generate::Query::weight`]). The other process filters run inside
//! the search through the generator's filter, after the budgets and at the earliest stage where
//! they are decided: bridges and special topologies on untyped topologies, and the special
//! topologies again with the masses of complete graphs.

use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet},
    ops::Bound,
};

use feynkit_graph::ExternalState;
use feynkit_model::{Model, Particle, ParticleId, VertexRuleId};
use linnet::half_edge::{
    generate::{
        Bridges, Countable, Diagram, DiagramView, DressingView, EdgeView, Legs, Plan, Query,
        SkeletonView, Structure, Stub, Theory,
    },
    involution::{Flow, Hedge},
};

#[cfg(test)]
use linnet::half_edge::{
    NodeIndex,
    involution::{HedgePair, Orientation},
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

/// A partial or complete untyped process topology. Fields identify particle/antiparticle pairs.
pub type GenerationSkeletonView<'a> = SkeletonView<'a, ParticleId, ExternalLabel>;
/// A partial or complete field and vertex-rule assignment on a complete process topology.
pub type GenerationDressingView<'a> = DressingView<'a, ParticleId, ExternalLabel>;
/// A complete process assignment at the arbitrary final-selection boundary.
pub type GenerationDiagramView<'a> = DiagramView<'a, ParticleId, ExternalLabel>;

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

    /// The generator of a request: the legs of `externals`, these signatures, the loop orders and
    /// structural limits of `options`, and its coupling orders and fermion-loop range as budgets;
    /// Typed hooks prune it after the budgets.
    pub(crate) fn plan(
        &self,
        model: &Model,
        externals: &ProcessExternals,
        options: &GenerationOptions,
        skeletons: impl Fn(&GenerationSkeletonView<'_>) -> bool + Send + Sync + 'static,
        dressings: impl Fn(&GenerationDressingView<'_>) -> bool + Send + Sync + 'static,
    ) -> Result<Option<Plan<ParticleId, ExternalLabel>>, GenerationError> {
        let bounds = CouplingBounds::new(options);
        let mut columns = Vec::new();
        for (name, ..) in &bounds.0 {
            if !columns.contains(name) {
                columns.push(name.clone());
            }
        }
        let mut theory = Theory::builder().weights(columns.iter().cloned());
        for (&rule, stubs) in self.rules.iter().zip(&self.stubs) {
            let orders = if columns.is_empty() {
                BTreeMap::new()
            } else {
                model.vertex_rule_by_id(rule)?.coupling_orders(model)
            };
            let weights = columns
                .iter()
                .map(|name| {
                    i32::try_from(orders.get(name).copied().unwrap_or(0))
                        .map_err(|_| GenerationError::ArithmeticOverflow("a coupling order"))
                })
                .collect::<Result<Vec<_>, _>>()?;
            theory = theory.vertex(stubs.iter().cloned(), weights);
        }
        let theory = theory.build()?;
        let legs = externals
            .legs()
            .fold(Legs::new(), |legs, (label, stub)| legs.leg(label, stub));
        let structure = Structure {
            self_loops: options.allow_self_loops,
            zero_flow_edges: options.allow_zero_flow_edges,
            max_vertices: options.max_vertices,
            bridges: options
                .max_bridges(*options.loop_count().start())
                .map_or(Bridges::Any, Bridges::AtMost),
        };
        let mut query = Query::new(&legs)
            .loops(options.loop_count())
            .structure(structure)
            .prune_skeletons(skeletons)
            .prune_dressings(dressings);
        // Repeated constraints on one declared order must all hold.
        let bound = |order: usize| i32::try_from(order).unwrap_or(i32::MAX);
        for (name, minimum, maximum) in &bounds.0 {
            let maximum =
                maximum.map_or(Bound::Unbounded, |maximum| Bound::Included(bound(maximum)));
            query = query.weight(name, (Bound::Included(bound(*minimum)), maximum));
        }
        // Fermion loops are closed lines of fermions; ghost loops do not count.
        let fermions = || {
            (0..model.particles().len())
                .map(ParticleId::from_index)
                .filter(|&id| {
                    model
                        .particle_by_id(id)
                        .is_ok_and(|particle| particle.is_fermion() && !particle.is_ghost())
                })
                .collect()
        };
        for filter in &options.graph_filters {
            if let GenerationFilter::FermionLoopCountRange((minimum, maximum)) = filter {
                query = query.count(Countable::ClosedLines(fermions()), *minimum..=*maximum);
            }
        }
        if let Some(select) = options.select.clone() {
            query = query.select(move |view| select(view));
        }
        if let Some((mask, _)) = &options.deferred_hooks {
            query = query.deferred_hooks(*mask);
        }
        match Plan::new(&theory, &query) {
            Ok(plan) => Ok(Some(plan)),
            Err(linnet::half_edge::generate::GenerationError::Infeasible(_)) => Ok(None),
            Err(error) => Err(error.into()),
        }
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

    /// The number of distinct diagrams that permuting interchangeable sewn pairs makes of
    /// `diagram`, or `None` when one of them has a smaller canonical form.
    ///
    /// The generator labels sewn pairs exactly, so such an orbit of diagrams is generated once
    /// per member. The member with the smallest canonical form represents the orbit, and the orbit
    /// size, the order of the permutation group over the order of the stabilizer of the
    /// diagram's isomorphism class, is the multiplicity the representative stands for.
    pub(crate) fn sewn_pair_orbit(
        &self,
        diagram: &Diagram<ParticleId, ExternalLabel>,
    ) -> Option<u64> {
        let mut stabilizer = 1;
        for permutation in &self.sewn_pair_permutations {
            let image = diagram.graph().canonical_form_by(
                |&rule| rule,
                |field| *field,
                |label| label.map(|label| permutation.get(&label).copied().unwrap_or(label)),
            );
            match image.cmp(diagram.canonical_form()) {
                Ordering::Less => return None,
                Ordering::Equal => stabilizer += 1,
                Ordering::Greater => {}
            }
        }
        // The identity is the one permutation that `sewn_pair_permutations` leaves out.
        Some((self.sewn_pair_permutations.len() as u64 + 1) / stabilizer)
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
        let canonical = self
            .classes
            .values()
            .any(|externals| externals.len() > 1)
            .then(|| {
                diagram
                    .graph()
                    .canonize_by(|&rule| rule, |field| *field, |label| *label)
            });
        // Compact legs use the materializer's native half-edge slot order. Canonical
        // positions are demanded only to label interchangeable external momenta.
        let mut legs = diagram
            .legs()
            .enumerate()
            .map(|(slot, (vertex, label, stub))| {
                let class = (
                    *label,
                    Stub {
                        field: *stub.field,
                        flow: stub.flow,
                    },
                );
                let hedge = Hedge(slot);
                let order = canonical
                    .as_ref()
                    .map_or(hedge, |form| form.canonical_hedge(hedge));
                (
                    class,
                    order,
                    slot,
                    vertex,
                    stub.flow.unwrap_or(Flow::Source),
                )
            })
            .collect::<Vec<_>>();
        let mut edges = vec![None; legs.len()];
        legs.sort_by(|left, right| (&left.0, left.1).cmp(&(&right.0, right.1)));
        let mut next = BTreeMap::<&(ExternalLabel, ParticleStub), usize>::new();
        for (class, _, slot, vertex, flow) in &legs {
            let index = next.entry(class).or_default();
            let external = Endpoint::External(self.classes[class][*index]);
            *index += 1;
            let vertex = Endpoint::Vertex(*vertex);
            let (source, sink) = match flow {
                Flow::Sink => (external, vertex),
                Flow::Source => (vertex, external),
            };
            edges[*slot] = Some(GeneratedEdge {
                source,
                sink,
                particle: class.1.field,
                directed: class.1.flow.is_some(),
            });
        }
        edges.extend(diagram.edges().map(|edge| {
            let stub = edge.stub.expect("a complete diagram assigns every edge");
            let [mut source, mut sink] = edge.vertices.map(Endpoint::Vertex);
            if stub.flow == Some(Flow::Sink) {
                std::mem::swap(&mut source, &mut sink);
            }
            Some(GeneratedEdge {
                source,
                sink,
                particle: *stub.field,
                directed: stub.flow.is_some(),
            })
        }));
        GeneratedGraph::new(
            diagram
                .rules()
                .iter()
                .map(|&rule| signatures.rules[rule as usize]),
            edges
                .into_iter()
                .map(|edge| edge.expect("every compact edge is placed")),
        )
    }

    // Independent old materialized boundary, retained only as a compact-consumer oracle.
    #[cfg(test)]
    fn materialized_reference(
        &self,
        diagram: &Diagram<ParticleId, ExternalLabel>,
        signatures: &VertexSignatures,
    ) -> GeneratedGraph {
        let graph = diagram.graph();
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
            (0..graph.node_bound()).map(|node| signatures.rules[graph[NodeIndex(node)]]),
            edges
                .into_iter()
                .map(|edge| edge.expect("every edge is placed")),
        )
    }
}

/// The process filters evaluated during the search; see the module documentation.
pub(crate) struct Pruner {
    one_particle_irreducible_loops: bool,
    vetoes: SpecialTopologyVetoes,
    /// Per particle index: massless.
    massless: Vec<bool>,
}

impl Pruner {
    pub(crate) fn new(model: &Model, options: &GenerationOptions) -> Self {
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
        Self {
            one_particle_irreducible_loops: options
                .graph_filters
                .contains(&GenerationFilter::LoopOneParticleIrreducible),
            vetoes,
            massless: (0..model.particles().len())
                .map(|index| model.particle_is_massless(ParticleId::from_index(index)))
                .collect(),
        }
    }

    pub(crate) fn keeps_skeleton(&self, view: &GenerationSkeletonView<'_>) -> bool {
        !(self.one_particle_irreducible_loops
            && view.loop_order() > 0
            && !view.bridges().is_empty())
            && !(view.is_topology_complete()
                && !self.vetoes.is_empty()
                && self.vetoes.vetoes(
                    &self.plain(
                        view.n_vertices(),
                        (0..view.n_vertices())
                            .flat_map(|v| view.legs(v).map(move |(_, stub)| (v, *stub.field))),
                        view.edges(),
                    ),
                ))
    }

    pub(crate) fn keeps_dressing(&self, view: &GenerationDressingView<'_>) -> bool {
        !(view.is_complete()
            && !self.vetoes.is_empty()
            && self.vetoes.vetoes(
                &self.plain(
                    view.n_vertices(),
                    (0..view.n_vertices())
                        .flat_map(|v| view.legs(v).map(move |(_, stub)| (v, *stub.field))),
                    view.edges(),
                ),
            ))
    }

    /// Compact topology data as a plain graph, with masses where fields are known.
    fn plain<'a>(
        &self,
        vertices: usize,
        legs: impl Iterator<Item = (usize, ParticleId)>,
        edges: impl Iterator<Item = EdgeView<'a, ParticleId>>,
    ) -> PlainGraph {
        let mut graph = PlainGraph::default();
        for _ in 0..vertices {
            graph.add_vertex();
        }
        for (vertex, field) in legs {
            graph.add_leg(vertex, self.massless[field.index()]);
        }
        for edge in edges {
            graph.add_edge(
                (edge.vertices[0], edge.vertices[1]),
                edge.stub.map(|stub| self.massless[stub.field.index()]),
            );
        }
        graph
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
    fn compact_inputs_preserve_materialized_edges_slots_and_external_order() {
        use crate::generation::{Generator, PreparedSearch, cases::Case};
        for name in [
            "qcd_gggg_L1",
            "qcd_gggg_sym_L1",
            "xs_qed_aaee_sym_lr_L2",
            "qed_eemm_L1",
            "xs_qed_eeaa_L2",
            "qcd_tad_L2",
            "qcd_vac_L2",
        ] {
            let case = Case::named(name);
            let generator = Generator::new(case.model());
            let PreparedSearch {
                signatures,
                externals,
                plan,
                ..
            } = generator
                .search(
                    &case.process(),
                    &case.options(&generator.model, 1),
                    case.kind,
                    |_, _| {},
                    |_, _| {},
                )
                .unwrap();
            let mut count = 0;
            let run = linnet::half_edge::generate::Run::default();
            for diagram in plan.as_ref().unwrap().diagrams(&run).unwrap() {
                assert_eq!(
                    externals.graph(&diagram, &signatures),
                    externals.materialized_reference(&diagram, &signatures),
                    "{name} diagram {count}"
                );
                count += 1;
            }
            assert!(count > 0, "{name}");
        }
    }

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
