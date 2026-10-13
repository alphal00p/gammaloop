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
    collections::{BTreeMap, BTreeSet, HashMap},
    hash::Hash,
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
use crate::options::ResolvedRestrictions;
use crate::{Classes, Request};

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
    pub(crate) fn new(model: &Model, options: &Request) -> Result<Self, GenerationError> {
        let restrictions = options.restrictions.resolve(model)?;
        let bounds = CouplingBounds::new(options);
        let mut signatures = Self {
            rules: vec![],
            stubs: vec![],
        };
        for (rule_index, rule) in model.vertex_rules().iter().enumerate() {
            let rule_id = model.vertex_rule_id_at(rule_index)?;
            if !restrictions.vertex(rule_id) {
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
            if rule
                .particles
                .iter()
                .any(|&particle| !restrictions.particle(particle))
            {
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

    /// Only properties observed by the active physics request refine the core's structural
    /// field orbits. The integer colours are local equality labels, not model identities.
    fn class_distinctions(
        &self,
        model: &Model,
        options: &Request,
    ) -> Result<(Vec<(ParticleId, u32)>, Vec<(usize, u32)>), GenerationError> {
        let cross_section = options.kind == crate::GenerationType::CrossSection;
        let cut_states = if cross_section {
            options.process.covariant_cut_states(model, options.kind)?
        } else {
            Vec::new()
        };
        let unresolved = if cross_section {
            super::unresolved_cut_content(model, options)?
        } else {
            None
        };
        let distinguish_mass = [
            options.topology.self_energies,
            options.topology.tadpoles,
            options.topology.snails,
        ]
        .iter()
        .any(|veto| veto.massive != veto.massless);
        let mut field_keys = HashMap::new();
        let fields = self
            .stubs
            .iter()
            .flatten()
            .map(|stub| stub.field)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .map(|field| {
                let particle = model.particle_by_id(field)?;
                let anti = model.particle_by_id(particle.antiparticle)?;
                // A directed field represents both charges. Preserve charge and multiplicity
                // separately for every allowed final state, including covariant alternatives.
                let cut_content = cut_states
                    .iter()
                    .map(|state| {
                        [particle.pdg_code, anti.pdg_code]
                            .map(|pdg| state.iter().filter(|&&value| value == pdg).count())
                    })
                    .collect::<Vec<_>>();
                let unresolved_content = unresolved.as_ref().map(|content| {
                    [field, particle.antiparticle].map(|id| content.particles.contains(&id))
                });
                let key = (
                    (particle.spin, particle.color, anti.spin, anti.color),
                    distinguish_mass.then(|| {
                        [field, particle.antiparticle].map(|id| model.particle_is_massless(id))
                    }),
                    cut_content,
                    unresolved_content,
                );
                Ok((field, Self::equality_colour(&mut field_keys, key)?))
            })
            .collect::<Result<Vec<_>, GenerationError>>()?;

        let mut rule_keys = HashMap::new();
        let rules = self
            .rules
            .iter()
            .enumerate()
            .map(|(local_rule, &rule)| {
                let rule = model.vertex_rule_by_id(rule)?;
                let lorentz = rule
                    .lorentz_structures
                    .iter()
                    .map(|&id| {
                        let structure = model.lorentz_structure_by_id(id)?;
                        Ok((structure.spins.clone(), structure.structure.clone()))
                    })
                    .collect::<Result<Vec<_>, GenerationError>>()?;
                let occupancy = rule
                    .couplings
                    .iter()
                    .map(|row| row.iter().map(Option::is_some).collect::<Vec<_>>())
                    .collect::<Vec<_>>();
                let orders = if cross_section {
                    rule.coupling_orders(model)
                } else {
                    BTreeMap::new()
                };
                let cut_orders = options
                    .cuts
                    .amplitude
                    .orders
                    .plain
                    .keys()
                    .map(|name| orders.get(name).copied().unwrap_or(0))
                    .collect::<Vec<_>>();
                // Reversing an edge and conjugating its particle preserves both local
                // signatures. Only the optional global amplitude-side exchange conjugates
                // a rule; then CP repair compares its complete coupling-order map.
                let cp_orders =
                    (cross_section && options.cuts.symmetrize_left_right).then_some(orders);
                let key = (
                    rule.color_structures.clone(),
                    lorentz,
                    occupancy,
                    cut_orders,
                    cp_orders,
                );
                // Keys use the surviving Theory position; model rule IDs can have holes.
                Ok((local_rule, Self::equality_colour(&mut rule_keys, key)?))
            })
            .collect::<Result<Vec<_>, GenerationError>>()?;
        Ok((fields, rules))
    }

    /// Intern structured equality without printing symbolic expressions or assigning stable IDs.
    fn equality_colour<K: Eq + Hash>(
        colours: &mut HashMap<K, u32>,
        key: K,
    ) -> Result<u32, GenerationError> {
        if let Some(&colour) = colours.get(&key) {
            return Ok(colour);
        }
        let colour = u32::try_from(colours.len())
            .map_err(|_| GenerationError::ArithmeticOverflow("a field-class distinction"))?;
        colours.insert(key, colour);
        Ok(colour)
    }

    /// The generator of a request: the legs of `externals`, these signatures, the loop orders and
    /// structural limits of `options`, and its coupling orders and fermion-loop range as budgets;
    /// Typed hooks prune it after the budgets.
    pub(crate) fn plan(
        &self,
        model: &Model,
        externals: &ProcessExternals,
        options: &Request,
        skeletons: impl Fn(&GenerationSkeletonView<'_>) -> bool + Send + Sync + 'static,
        dressings: impl Fn(&GenerationDressingView<'_>) -> bool + Send + Sync + 'static,
    ) -> Result<Plan<ParticleId, ExternalLabel>, GenerationError> {
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
            self_loops: options.topology.self_loops,
            zero_flow_edges: options.topology.zero_flow_edges,
            max_vertices: options.topology.max_vertices,
            bridges: options.topology.bridges,
        };
        let mut query = Query::new(&legs)
            .loops(options.loops.range())
            .structure(structure)
            .classes(options.classes)
            .prune_skeletons(skeletons)
            .prune_dressings(dressings);
        if options.classes == Classes::Auto {
            let (fields, rules) = self.class_distinctions(model, options)?;
            query = query.class_distinctions(fields, rules);
        }
        // Repeated constraints on one declared order must all hold.
        let bound = |order: usize| {
            i32::try_from(order)
                .map_err(|_| GenerationError::ArithmeticOverflow("a coupling-order bound"))
        };
        for (name, minimum, maximum) in &bounds.0 {
            let maximum = maximum
                .map(|maximum| bound(maximum).map(Bound::Included))
                .transpose()?
                .unwrap_or(Bound::Unbounded);
            query = query.weight(name, (Bound::Included(bound(*minimum)?), maximum));
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
        let (minimum, maximum) = options.topology.fermion_loops.bounds();
        if minimum != 0 || maximum.is_some() {
            query = query.count(
                Countable::ClosedLines(fermions()),
                (
                    Bound::Included(minimum),
                    maximum.map_or(Bound::Unbounded, Bound::Included),
                ),
            );
        }
        if let Some(select) = options.select.clone() {
            query = query.select(move |view| select(view));
        }
        if let Some((mask, _)) = &options.deferred_hooks {
            query = query.deferred_hooks(*mask);
        }
        Plan::new(&theory, &query).map_err(Into::into)
    }
}

/// The graph-wide coupling-order bounds `(name, minimum, maximum)`.
struct CouplingBounds(Vec<(String, usize, Option<usize>)>);

impl CouplingBounds {
    fn new(options: &Request) -> Self {
        let mut bounds = options.orders.plain.clone();
        // Squared constraints on a forward graph apply directly to its summed order.
        if options.kind == crate::GenerationType::CrossSection {
            for (name, range) in &options.orders.squared {
                bounds
                    .entry(name.clone())
                    .and_modify(|old| *old = old.intersect(*range))
                    .or_insert(*range);
            }
        }
        Self(
            bounds
                .into_iter()
                .map(|(name, range)| {
                    let (min, max) = range.bounds();
                    (name, min, max)
                })
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
    pub(crate) fn vetoed(&self, restrictions: &ResolvedRestrictions) -> bool {
        self.legs
            .iter()
            .any(|(_, stub)| !restrictions.particle(stub.field))
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
    pub(crate) fn new(model: &Model, options: &Request) -> Self {
        let topology = &options.topology;
        let vetoes = SpecialTopologyVetoes {
            self_energy: topology
                .self_energies
                .enabled()
                .then_some(topology.self_energies),
            tadpoles: topology.tadpoles.enabled().then_some(topology.tadpoles),
            snails: topology.snails.enabled().then_some(topology.snails),
            factorized_loop_range: topology.factorized_loops,
        };
        Self {
            one_particle_irreducible_loops: matches!(topology.bridges, Bridges::OnePIInLoops),
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
    fn active_mass_and_cut_membership_split_physics_equality_colours() {
        use crate::{Loops, Process, TopologyRules, Veto};
        let model = Model::qed();
        let up = model.particle_id("u").unwrap();
        let top = model.particle_id("t").unwrap();
        assert!(model.particle_is_massless(up));
        assert!(!model.particle_is_massless(top));
        let request = Process::new(["a"], ["a"]).amplitude(Loops::total(2..=2));
        let colours = |request: &Request| {
            let signatures = VertexSignatures::new(&model, request).unwrap();
            signatures
                .class_distinctions(&model, request)
                .unwrap()
                .0
                .into_iter()
                .collect::<BTreeMap<_, _>>()
        };
        for veto in [
            Veto::default(),
            Veto::all(),
            Veto::massive(),
            Veto::massless(),
        ] {
            for kind in 0..3 {
                let topology = match kind {
                    0 => TopologyRules::new().with_self_energies(veto),
                    1 => TopologyRules::new().with_tadpoles(veto),
                    _ => TopologyRules::new().with_snails(veto),
                };
                let colours = colours(&request.clone().with_topology(topology));
                assert_eq!(
                    colours[&up] == colours[&top],
                    veto.massive == veto.massless,
                    "only an active mass-selective veto may distinguish these field colours"
                );
            }
        }
        let cross_section =
            Process::new(["a", "a"], ["u", "u~"]).cross_section(Loops::forward(1..=1));
        let colours = colours(&cross_section);
        assert_ne!(
            colours[&up], colours[&top],
            "requested cut membership is active"
        );
    }

    #[test]
    fn sparse_original_rules_keep_dense_local_equality_colours() {
        use crate::{Loops, Process, Restrict, VertexSelector};
        let model = Model::qed();
        let surviving = (1..model.vertex_rules().len())
            .step_by(2)
            .map(|index| model.vertex_rule_id_at(index).unwrap())
            .collect::<Vec<_>>();
        assert!(surviving.len() > 1);
        let selectors = surviving
            .iter()
            .map(|&id| VertexSelector::by_id(&model, id).unwrap());
        let request = Process::new(["a"], ["a"])
            .amplitude(Loops::total(1..=1))
            .with_restrict(Restrict::new().only_vertices(selectors));
        let signatures = VertexSignatures::new(&model, &request).unwrap();
        assert_eq!(signatures.rules, surviving);
        let (_, colours) = signatures.class_distinctions(&model, &request).unwrap();
        assert_eq!(
            colours.iter().map(|&(index, _)| index).collect::<Vec<_>>(),
            (0..surviving.len()).collect::<Vec<_>>()
        );
    }

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
                .search(&case.request(&generator.model), |_, _| {}, |_, _| {})
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
