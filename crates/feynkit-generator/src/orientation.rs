//! Closed fermion-loop orientation sums.
//!
//! For a diagram `D` with closed, charged, non-ghost fermion loops
//! `l_1..l_k`, let `D_S` be `D` with the physical flow of the loops in `S`
//! reversed. `D_S` keeps `D`'s storage: every loop edge swaps particle and
//! antiparticle while its endpoints stay put, so edge IDs, half-edge IDs and
//! the momentum symbols `Q(edge)` are shared by all `2^k` variants. The
//! distinct oriented diagrams of one family weigh `1/|Aut(U)|` times the sum
//! over the labelled variants, so the whole family cancels exactly when
//! `sum_S N(D_S)` vanishes. Every member of a family gets the same verdict, so
//! each diagram is tested on its own.
//!
//! Reversing a closed loop transposes every spinor factor of the loop. The
//! loop is closed, so charge conjugation turns the transposed factors back
//! into the original ones exactly: `C gamma^T C^-1 = -gamma` and the spinor
//! metric is symmetric. Spinor factors outside this allow-list fail closed.
//! Color factors are left as built: the color algebra then decides whether
//! the transposed color traces cancel.
//!
//! Every variant, including `D_0`, is rebuilt from the topology. The rebuild
//! fills two legs of one species at a vertex first-fit, so `D_0` may put them
//! in swapped slots relative to `D` (background-field QCD: the two B legs of
//! V_BBQQ or V_ccBB). Such a variant stands for `D` only when every slot
//! exchange is certified to leave the rule's complete numerator invariant
//! (`Generator::leg_exchange_symmetric`, exact, cached per rule and leg
//! pair); its numerator is then used unchanged (`input_frame_mismatch`).
//! Reversed variants need the same certificate for every pair of identical
//! legs at their vertices. Anything uncertified is `Skipped`.
//!
//! What a removed zero means depends on the loop's colour:
//! - colourless loops (QED, Furry's theorem): the orientation sum vanishes
//!   identically, independent of any gauge group;
//! - coloured loops: the zero is certified by
//!   `grouping::exact_color_zero` at the model's numeric representation
//!   dimensions, i.e. at SU(3) for QCD (colour-zero tier Z3). That colour
//!   algebra uses local rules, including Fierz, which hold only for SU(N); a
//!   removed coloured diagram is therefore not claimed to vanish for every
//!   group (Z1) or for every SU(N) (Z2).
//!
//! Sampled values are never accepted as a zero; they can only show a nonzero.

use std::collections::{BTreeMap, BTreeSet};

use feynkit_graph::{
    DiagramEndpoint, DiagramHalfEdge, EdgeId, FeynmanDiagram, NumeratorInterface, VertexId,
};
use feynkit_model::VertexRuleId;
use idenso::{
    CookMode, CookSettings,
    color::ColorSimplifySettings,
    dirac::AGS,
    representations::Bispinor,
    tensor::{AlgebraContraction, AlgebraSettings},
};
use linnet::half_edge::involution::Flow;
use spenso::{network::library::symbolic::ETS, structure::representation::RepName};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol, ZeroTestSettings},
    graph::Graph,
};

use super::{
    ColoredNode, ColoredTopology, EdgeColor, ExternalNode, FermionSigns, GenerationError,
    GenerationOptions, GenerationType, Generator, NumeratorHalfEdge, NumeratorInstantiation,
    NumeratorOwner, ResolvedProcess, grouping,
};

/// Largest number of independently reversed loops; `2^k` variants are built.
const MAX_ORIENTED_LOOPS: usize = 12;

/// Outcome of the closed-loop orientation-sum test for one diagram.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OrientationSumVerdict {
    /// The diagram has no closed loop of a non-ghost, non-Majorana fermion.
    NoClosedLoops,
    /// The numerator summed over all loop orientations vanishes:
    /// group-independently for colourless loops, and at the model's numeric
    /// representation dimensions (SU(3) for QCD, tier Z3) for coloured loops.
    Zero,
    /// The orientation sum was not proven to vanish.
    NonZero,
    /// The test could not be certified for this diagram; it is retained.
    Skipped(String),
}

impl Generator {
    /// Orientation-sum verdict for a finalized diagram of this generator's model.
    pub(crate) fn orientation_sum_verdict_in(
        &self,
        diagram: &FeynmanDiagram,
        generation_type: GenerationType,
        options: &GenerationOptions,
    ) -> Result<OrientationSumVerdict, GenerationError> {
        let skip = |reason: &str| Ok(OrientationSumVerdict::Skipped(reason.to_owned()));
        let topology = match colored_topology(diagram) {
            Ok(topology) => topology,
            Err(reason) => return skip(&reason),
        };
        let mut loops = Vec::new();
        for edges in self.closed_fermion_loops(&topology.graph, false)? {
            let majorana = edges
                .iter()
                .filter(|edge| {
                    self.model
                        .particle_is_self_conjugate(topology.graph.edges()[**edge].data.particle)
                })
                .count();
            // Majorana lines carry no direction; pure Majorana loops are not
            // varied. A loop that mixes Majorana and Dirac lines has no
            // certified reversal (fermion-flow conventions), so it fails closed.
            if majorana == edges.len() {
                continue;
            }
            if majorana > 0 {
                return skip("a closed loop mixes Majorana and Dirac fermions");
            }
            loops.push(edges);
        }
        if loops.is_empty() {
            return Ok(OrientationSumVerdict::NoClosedLoops);
        }
        if generation_type != GenerationType::Amplitude {
            return skip("orientation sums are only certified for amplitudes");
        }
        if loops.len() > MAX_ORIENTED_LOOPS {
            return skip("too many closed fermion loops");
        }
        // Fail closed for a loop species without any propagator. A species
        // that borrows its antiparticle's propagator is fine: the reversed
        // template is read with swapped slots and momentum -Q(edge), which is
        // the antiparticle line exactly.
        let species = loops
            .iter()
            .flatten()
            .map(|edge| topology.graph.edges()[*edge].data.particle)
            .collect::<BTreeSet<_>>();
        for particle in species {
            let particle = self.model.particle_by_id(particle)?;
            if self.model.particle_propagator(particle)?.is_none() {
                return skip(&format!(
                    "{} has no propagator; reversal is not certified",
                    particle.name
                ));
            }
        }

        let hedges = HedgeFrame::of(diagram);
        let loop_hedges = loops
            .iter()
            .map(|edges| {
                edges
                    .iter()
                    .flat_map(|edge| hedges.of_edge(EdgeId(*edge)))
                    .collect::<BTreeSet<_>>()
            })
            .collect::<Vec<_>>();
        let symbols = SpinorSymbols::new();
        // Admitted numerators carry reversibly cooked indices; the transport
        // and transposition below act on the explicit half-edge owners.
        let cooking = CookSettings::indices()
            .with_mode(CookMode::ReversibleEncoding)
            .with_representation_payloads(true, true);
        let mut terms = Vec::with_capacity(1 << loops.len());
        for mask in 0usize..(1 << loops.len()) {
            let variant = match self.oriented_variant(diagram, &topology, &loops, mask, options)? {
                Ok(variant) => variant,
                Err(reason) => return skip(&reason),
            };
            if mask == 0 {
                // N(D_0) stands for N(D): see `input_frame_mismatch`.
                if let Some(reason) = self.input_frame_mismatch(&variant, diagram)? {
                    return skip(&reason);
                }
            } else {
                // Every variant must keep the input's vertices and edge
                // endpoints, so `Q(edge)` and the half-edge transport mean the
                // same thing.
                if !same_endpoints(&variant, diagram) {
                    return skip("a reversed variant moved an edge endpoint");
                }
                // The rebuild fills identical legs first-fit, as generation
                // does; the reversed family member's numerator is independent
                // of that choice only for certified-symmetric rules.
                if let Some(reason) = self.uncertified_identical_legs(&variant)? {
                    return skip(&reason);
                }
            }
            let explicit = cooking.uncook(variant.numerator().expression().as_view());
            let relabeled = match HedgeFrame::of(&variant).relabel_into(&hedges, &explicit) {
                Some(relabeled) => relabeled,
                None => return skip("the reversed diagram has a different half-edge frame"),
            };
            let reversed = loop_hedges
                .iter()
                .enumerate()
                .filter(|(index, _)| mask & (1 << index) != 0)
                .flat_map(|(_, hedges)| hedges.iter().copied())
                .collect::<BTreeSet<_>>();
            match symbols.untranspose(&relabeled, &reversed) {
                Some(term) => terms.push(term),
                None => {
                    return skip(
                        "a reversed loop carries a spinor structure without a known charge-conjugation parity",
                    );
                }
            }
        }
        let sum = Atom::add_many(terms);
        if sum.is_zero() {
            return Ok(OrientationSumVerdict::Zero);
        }
        // Admission cooks the explicit indices again and infers the open
        // interface of the sum, which every variant shares. A zero found
        // below holds at the model's numeric representation dimensions (SU(3),
        // tier Z3), not for every group: see the module documentation.
        let numerator = NumeratorInterface::admit(sum, None)?;
        Ok(
            if grouping::exact_color_zero(diagram, &numerator, &self.model)? {
                OrientationSumVerdict::Zero
            } else {
                OrientationSumVerdict::NonZero
            },
        )
    }

    /// Rebuild `topology` with the loops selected by `mask` reversed by
    /// particle/antiparticle exchange in place. `Err` carries a skip reason.
    fn oriented_variant(
        &self,
        diagram: &FeynmanDiagram,
        topology: &ColoredTopology,
        loops: &[Vec<usize>],
        mask: usize,
        options: &GenerationOptions,
    ) -> Result<Result<FeynmanDiagram, String>, GenerationError> {
        let toggled = loops
            .iter()
            .enumerate()
            .filter(|(index, _)| mask & (1 << index) != 0)
            .flat_map(|(_, edges)| edges.iter().copied())
            .collect::<BTreeSet<_>>();
        let mut graph = Graph::new();
        for node in topology.graph.nodes() {
            graph.add_node(node.data.clone());
        }
        for (index, edge) in topology.graph.edges().iter().enumerate() {
            let mut data = edge.data;
            if toggled.contains(&index) {
                data.particle = self.model.particle_by_id(data.particle)?.antiparticle;
            }
            graph
                .add_edge(edge.vertices.0, edge.vertices.1, edge.directed, data)
                .map_err(|error| GenerationError::GraphConstruction(error.to_string()))?;
        }
        let graph = match ResolvedProcess::repair_cp_interaction_rules(graph, &self.model) {
            Ok(graph) => graph,
            Err(GenerationError::InvalidCpInteractionRule { interaction, .. }) => {
                return Ok(Err(format!(
                    "no unique interaction for the reversed flow at {interaction}"
                )));
            }
            Err(error) => return Err(error),
        };
        let variant = ColoredTopology {
            graph,
            symmetry: topology.symmetry,
            multiplicity: 1,
            cut_partitions: Vec::new(),
            topology_threshold_candidates: Vec::new(),
        };
        let variant = match self.to_diagram(
            diagram.name().to_owned(),
            GenerationType::Amplitude,
            variant,
            FermionSigns::default(),
            options,
        ) {
            Ok(variant) => variant,
            Err(GenerationError::MissingInteractionLeg { interaction, .. }) => {
                return Ok(Err(format!(
                    "interaction {interaction} has no slot for the reversed flow"
                )));
            }
            Err(error) => return Err(error),
        };
        Ok(Ok(variant))
    }
}

/// Topology view of a finalized amplitude diagram in its own storage order:
/// node `i` is vertex `i`, edge `i` is edge `i`, and each external half-edge
/// gets an external node after the vertices. `Err` carries a skip reason.
fn colored_topology(diagram: &FeynmanDiagram) -> Result<ColoredTopology, String> {
    let mut graph = Graph::new();
    for (position, (vertex, data)) in diagram.vertices().enumerate() {
        let Some(rule) = data.interaction else {
            return Err(format!("vertex {} has no interaction", vertex.0));
        };
        if vertex.0 != position {
            return Err("vertex identifiers are not contiguous".to_owned());
        }
        graph.add_node(ColoredNode::Interaction(rule));
    }
    let mut endpoints = Vec::new();
    for (edge, ends, data) in diagram.edges() {
        if edge.0 != endpoints.len() {
            return Err("edge identifiers are not contiguous".to_owned());
        }
        if data.is_dummy || data.propagator_power != 1 {
            return Err(format!("edge {} is not a plain propagator", edge.0));
        }
        let color = EdgeColor {
            particle: data.particle,
            direction: None,
        };
        let (source, target) = match (ends.source, ends.target, &data.external) {
            (Some(source), Some(target), None) => (source.0, target.0),
            (source, target, Some(external)) if source.is_some() != target.is_some() => {
                let node = graph.add_node(ColoredNode::External(ExternalNode {
                    index: external.index,
                    state: external.state,
                    particle: data.particle,
                }));
                match (source, target) {
                    (None, Some(target)) => (node, target.0),
                    (Some(source), None) => (source.0, node),
                    _ => unreachable!("exactly one endpoint is attached"),
                }
            }
            _ => return Err(format!("edge {} has an unsupported attachment", edge.0)),
        };
        endpoints.push((source, target, data.directed, color));
    }
    for (source, target, directed, color) in endpoints {
        graph
            .add_edge(source, target, directed, color)
            .map_err(|error| error.to_string())?;
    }
    Ok(ColoredTopology {
        graph,
        symmetry: diagram.symmetry_factor(),
        multiplicity: 1,
        cut_partitions: Vec::new(),
        topology_threshold_candidates: Vec::new(),
    })
}

impl Generator {
    /// Whether the rebuilt S = 0 variant may stand for the input diagram.
    /// `None` accepts it, `Some` carries a skip reason.
    ///
    /// The variant must have the input's vertices, interactions, edges,
    /// endpoints and species. Its interaction slots may differ from the
    /// input's only by a permutation `pi` of identical legs at a vertex: the
    /// rebuild fills identical legs first-fit in its own edge order, which
    /// need not be the order generation saw.
    ///
    /// Why the variant's numerator is then used unchanged, and not mapped
    /// through `pi`: the interaction numerator at a vertex is
    /// `N_V(l_1, .., l_n)`, with `l_s` the leg data (edge, half-edge, flow)
    /// put in slot `s`. The input has `N_V(l)`, the variant `N_V(l o pi^-1)`.
    /// Every transposition of `pi` exchanges two legs of one species and is
    /// certified (`leg_exchange_symmetric`) to leave the complete rule
    /// numerator invariant for arbitrary leg data, so `N_V(l o pi^-1) =
    /// N_V(l)` exactly, and so are the diagram numerators. Half-edge indices
    /// are transported by (edge, endpoint), not by slot, so no relabelling is
    /// left to do; renaming indices by `pi` would instead give one half-edge
    /// another's indices while `Q(edge)` stayed put.
    fn input_frame_mismatch(
        &self,
        variant: &FeynmanDiagram,
        input: &FeynmanDiagram,
    ) -> Result<Option<String>, GenerationError> {
        if !same_structure_up_to_slots(variant, input) {
            return Ok(Some(
                "the rebuilt diagram does not reproduce the input frame".to_owned(),
            ));
        }
        let (input_slots, variant_slots) = (slot_tables(input), slot_tables(variant));
        for (vertex, data) in input.vertices() {
            let Some(rule_id) = data.interaction else {
                return Ok(Some(format!("vertex {} has no interaction", vertex.0)));
            };
            let rule = self.model.vertex_rule_by_id(rule_id)?;
            let (Some(input_legs), Some(variant_legs)) =
                (input_slots.get(&vertex), variant_slots.get(&vertex))
            else {
                return Ok(Some(format!("vertex {} has no legs", vertex.0)));
            };
            if input_legs == variant_legs {
                continue;
            }
            // pi: input slot -> variant slot of the same half-edge.
            let mut permutation = BTreeMap::new();
            for (hedge, slot) in input_legs {
                let Some(image) = variant_legs.get(hedge) else {
                    return Ok(Some(format!(
                        "vertex {} lost a leg in the rebuild",
                        vertex.0
                    )));
                };
                let species = |slot: usize| rule.particles.get(slot).copied();
                if species(*slot).is_none() || species(*slot) != species(*image) {
                    return Ok(Some(format!(
                        "the rebuild moved a leg of {} to a slot of another species",
                        rule.name
                    )));
                }
                permutation.insert(*slot, *image);
            }
            if permutation
                .keys()
                .ne(permutation.values().collect::<BTreeSet<_>>())
            {
                return Ok(Some(format!(
                    "the rebuilt slots of {} are not a permutation",
                    rule.name
                )));
            }
            for (a, b) in transpositions(&permutation) {
                if !self.leg_exchange_symmetric(rule_id, a, b)? {
                    return Ok(Some(format!(
                        "{} is not certified symmetric under exchanging its identical legs {} and {}",
                        rule.name,
                        a + 1,
                        b + 1
                    )));
                }
            }
        }
        Ok(None)
    }

    /// `Some` skip reason when a vertex of `diagram` has two legs of one
    /// species whose exchange is not certified to leave its rule invariant.
    fn uncertified_identical_legs(
        &self,
        diagram: &FeynmanDiagram,
    ) -> Result<Option<String>, GenerationError> {
        for (vertex, data) in diagram.vertices() {
            let Some(rule_id) = data.interaction else {
                return Ok(Some(format!("vertex {} has no interaction", vertex.0)));
            };
            let rule = self.model.vertex_rule_by_id(rule_id)?;
            for a in 0..rule.particles.len() {
                for b in a + 1..rule.particles.len() {
                    if rule.particles[a] == rule.particles[b]
                        && !self.leg_exchange_symmetric(rule_id, a, b)?
                    {
                        return Ok(Some(format!(
                            "{} is not certified symmetric under exchanging its identical legs {} and {}",
                            rule.name,
                            a + 1,
                            b + 1
                        )));
                    }
                }
            }
        }
        Ok(None)
    }

    /// Whether the complete numerator of interaction `rule_id` (every color x
    /// Lorentz channel times its coupling) is invariant under exchanging legs
    /// `a` and `b` of one species, for arbitrary leg data. Certified once per
    /// (rule, leg pair) for this generator and cached; anything that cannot
    /// be certified exactly is `false`.
    pub(crate) fn leg_exchange_symmetric(
        &self,
        rule_id: VertexRuleId,
        a: usize,
        b: usize,
    ) -> Result<bool, GenerationError> {
        let key = (rule_id, a.min(b), a.max(b));
        if let Some(certified) = self.numerators.leg_exchanges.read().unwrap().get(&key) {
            return Ok(*certified);
        }
        let certified = self.certify_leg_exchange(key.0, key.1, key.2)?;
        self.numerators
            .leg_exchanges
            .write()
            .unwrap()
            .insert(key, certified);
        Ok(certified)
    }

    /// The legs get distinct edges, half-edges and momenta `Q(edge)`, so the
    /// exchanged instantiation swaps every index and momentum of `a` and `b`
    /// at once. All legs are sinks; a source leg only negates its own momentum,
    /// which a generic `Q(edge)` already covers. The difference of the two
    /// numerators must vanish after expansion, or after the color algebra of
    /// `grouping::exact_color_zero` with the exact zero test; samples never
    /// certify.
    fn certify_leg_exchange(
        &self,
        rule_id: VertexRuleId,
        a: usize,
        b: usize,
    ) -> Result<bool, GenerationError> {
        let rule = self.model.vertex_rule_by_id(rule_id)?;
        if a == b {
            return Ok(true);
        }
        match (rule.particles.get(a), rule.particles.get(b)) {
            (Some(left), Some(right)) if left == right => {}
            _ => return Ok(false),
        }
        let legs = rule
            .particles
            .iter()
            .enumerate()
            .map(|(leg, particle)| {
                let particle = self.model.particle_by_id(*particle)?;
                Ok(NumeratorHalfEdge {
                    edge: leg,
                    hedge: 2 * leg + 1,
                    flow: Flow::Sink,
                    spin: particle.spin,
                    color: particle.color,
                })
            })
            .collect::<Result<Vec<_>, GenerationError>>()?;
        let mut exchanged = legs.clone();
        exchanged.swap(a, b);
        let numerator = |legs: &[NumeratorHalfEdge]| {
            self.interaction_expression(
                rule_id,
                &NumeratorInstantiation {
                    owner: NumeratorOwner::Vertex(0),
                    legs,
                    reversed: false,
                },
            )
        };
        let difference = (numerator(&legs)? - numerator(&exchanged)?).expand();
        if difference.is_zero() {
            return Ok(true);
        }
        let Ok(difference) = NumeratorInterface::admit(difference, None) else {
            return Ok(false);
        };
        let Ok(simplified) = difference.simplify_algebra(&AlgebraSettings {
            color: Some(ColorSimplifySettings::default().with_cof_dimension_invariants()),
            contract: AlgebraContraction::None,
            ..Default::default()
        }) else {
            return Ok(false);
        };
        Ok(simplified
            .expression()
            .expand()
            .collect_factors()
            .check_zero_with(&ZeroTestSettings::default().exact_fallback(true))
            .is_true())
    }
}

/// Same vertices, interactions, edges, endpoints and species; slots are
/// compared by `Generator::input_frame_mismatch`.
fn same_structure_up_to_slots(left: &FeynmanDiagram, right: &FeynmanDiagram) -> bool {
    let edges = |diagram: &FeynmanDiagram| {
        diagram
            .edges()
            .map(|(id, ends, edge)| (id, ends.source, ends.target, edge.particle, edge.directed))
            .collect::<Vec<_>>()
    };
    left.vertices()
        .map(|(id, vertex)| (id, vertex.interaction))
        .eq(right
            .vertices()
            .map(|(id, vertex)| (id, vertex.interaction)))
        && edges(left) == edges(right)
}

/// Interaction slot of every attached half-edge, per vertex.
fn slot_tables(
    diagram: &FeynmanDiagram,
) -> BTreeMap<VertexId, BTreeMap<(EdgeId, DiagramEndpoint), usize>> {
    let mut tables = BTreeMap::<_, BTreeMap<_, _>>::new();
    for (edge, ends, data) in diagram.edges() {
        if let Some(vertex) = ends.source {
            tables
                .entry(vertex)
                .or_default()
                .insert((edge, DiagramEndpoint::Source), data.source_slot().0);
        }
        if let Some(vertex) = ends.target {
            tables
                .entry(vertex)
                .or_default()
                .insert((edge, DiagramEndpoint::Target), data.target_slot().0);
        }
    }
    tables
}

/// Transpositions whose product is the permutation `slot -> image`: each
/// cycle `(c_0 c_1 .. c_k)` is `(c_0 c_k) .. (c_0 c_1)`. They stay inside one
/// species because the permutation does.
fn transpositions(permutation: &BTreeMap<usize, usize>) -> Vec<(usize, usize)> {
    let mut seen = BTreeSet::new();
    let mut result = Vec::new();
    for start in permutation.keys() {
        if !seen.insert(*start) {
            continue;
        }
        let mut next = permutation[start];
        while next != *start && seen.insert(next) {
            result.push((*start, next));
            next = permutation[&next];
        }
    }
    result
}

/// Same vertices, edges, edge endpoints and edge directedness. Species,
/// interactions and slots differ between reversed variants and are not compared.
fn same_endpoints(left: &FeynmanDiagram, right: &FeynmanDiagram) -> bool {
    let edges = |diagram: &FeynmanDiagram| {
        diagram
            .edges()
            .map(|(id, ends, edge)| (id, ends.source, ends.target, edge.directed))
            .collect::<Vec<_>>()
    };
    left.vertices()
        .map(|(id, _)| id)
        .eq(right.vertices().map(|(id, _)| id))
        && edges(left) == edges(right)
}

/// Half-edge IDs of a finalized diagram keyed by stable edge endpoints.
struct HedgeFrame(BTreeMap<(EdgeId, DiagramEndpoint), usize>);

impl HedgeFrame {
    fn of(diagram: &FeynmanDiagram) -> Self {
        let mut hedges = BTreeMap::new();
        for (edge, _, _) in diagram.edges() {
            for endpoint in [DiagramEndpoint::Source, DiagramEndpoint::Target] {
                if let Some(hedge) = diagram.half_edge_id(DiagramHalfEdge { edge, endpoint }) {
                    hedges.insert((edge, endpoint), hedge.0);
                }
            }
        }
        Self(hedges)
    }

    fn of_edge(&self, edge: EdgeId) -> impl Iterator<Item = usize> + '_ {
        [DiagramEndpoint::Source, DiagramEndpoint::Target]
            .into_iter()
            .filter_map(move |endpoint| self.0.get(&(edge, endpoint)).copied())
    }

    /// Rename this frame's half-edge indices into `target`'s frame.
    fn relabel_into(&self, target: &Self, expression: &Atom) -> Option<Atom> {
        if self.0.keys().ne(target.0.keys()) {
            return None;
        }
        let map = self
            .0
            .iter()
            .map(|(key, hedge)| (*hedge, target.0[key]))
            .collect::<BTreeMap<_, _>>();
        if map.iter().all(|(from, to)| from == to) {
            return Some(expression.clone());
        }
        let head = feynkit_graph::symbols::hedge_index();
        let mut unknown = false;
        let relabeled = expression.replace_map(|part, _, output| {
            let AtomView::Fun(function) = part else {
                return;
            };
            if function.get_symbol() != head {
                return;
            }
            let Some(hedge) = function
                .iter()
                .next()
                .and_then(|arg| usize::try_from(arg).ok())
            else {
                unknown = true;
                return;
            };
            let Some(mapped) = map.get(&hedge) else {
                unknown = true;
                return;
            };
            **output = FunctionBuilder::new(head)
                .add_arg(*mapped)
                .add_args(function.iter().skip(1))
                .finish();
        });
        (!unknown).then_some(relabeled)
    }
}

struct SpinorSymbols {
    bispinor: Symbol,
    gamma: Symbol,
    metric: Symbol,
    hedge: Symbol,
}

impl SpinorSymbols {
    fn new() -> Self {
        let index = (Bispinor {}).new_rep(4).to_symbolic([Atom::num(0)]);
        let bispinor = match index.as_view() {
            AtomView::Fun(function) => function.get_symbol(),
            _ => unreachable!("a represented index is a function"),
        };
        Self {
            bispinor,
            gamma: AGS.gamma,
            metric: ETS.metric,
            hedge: feynkit_graph::symbols::hedge_index(),
        }
    }

    /// Whether `index` is a bispinor slot owned by one of `hedges`.
    fn owned_spinor(&self, index: AtomView<'_>, hedges: &BTreeSet<usize>) -> bool {
        let AtomView::Fun(rep) = index else {
            return false;
        };
        rep.get_symbol() == self.bispinor
            && rep.iter().any(|arg| {
                matches!(arg, AtomView::Fun(owner)
                    if owner.get_symbol() == self.hedge
                        && owner
                            .iter()
                            .next()
                            .and_then(|hedge| usize::try_from(hedge).ok())
                            .is_some_and(|hedge| hedges.contains(&hedge)))
            })
    }

    /// Undo the transposition of every spinor factor on the reversed loops'
    /// half-edges by charge conjugation. `None` for an unknown structure.
    fn untranspose(&self, expression: &Atom, reversed: &BTreeSet<usize>) -> Option<Atom> {
        if reversed.is_empty() {
            return Some(expression.clone());
        }
        let mut unsupported = false;
        let result = expression.replace_map(|part, _, output| {
            let AtomView::Fun(function) = part else {
                return;
            };
            let head = function.get_symbol();
            if head == self.bispinor || head == self.hedge {
                return;
            }
            let args = function.iter().collect::<Vec<_>>();
            let owned = args
                .iter()
                .map(|arg| self.owned_spinor(*arg, reversed))
                .collect::<Vec<_>>();
            if !owned.iter().any(|owned| *owned) {
                return;
            }
            match (head, owned.as_slice()) {
                (head, [true, true, false]) if head == self.gamma => {
                    **output = -FunctionBuilder::new(head)
                        .add_arg(args[1])
                        .add_arg(args[0])
                        .add_arg(args[2])
                        .finish();
                }
                (head, [true, true]) if head == self.metric => {}
                _ => unsupported = true,
            }
        });
        (!unsupported).then_some(result)
    }
}

#[cfg(test)]
mod leg_exchange_tests {
    use feynkit_model::Model;

    use super::*;

    const QCD: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/orientation/qcd_l5/model.json"
    );

    fn qcd_definition() -> serde_json::Value {
        serde_json::from_str(&std::fs::read_to_string(QCD).unwrap()).unwrap()
    }

    fn identical_pairs(generator: &Generator) -> Vec<(String, usize, usize, bool)> {
        let mut pairs = Vec::new();
        for (index, rule) in generator.model.vertex_rules().iter().enumerate() {
            let id = generator.model.vertex_rule_id_at(index).unwrap();
            for a in 0..rule.particles.len() {
                for b in a + 1..rule.particles.len() {
                    if rule.particles[a] == rule.particles[b] {
                        let certified = generator.leg_exchange_symmetric(id, a, b).unwrap();
                        pairs.push((rule.name.clone(), a, b, certified));
                    }
                }
            }
        }
        pairs
    }

    #[test]
    fn background_qcd_rules_are_certified_symmetric() {
        let generator = Generator::new(Model::from_json(&qcd_definition().to_string()).unwrap());
        let pairs = identical_pairs(&generator);
        // V_QQQ 3, V_BQQ 1, V_QQQQ 6, V_BQQQ 3, V_BBQQ 2, V_ccBB 1.
        assert_eq!(pairs.len(), 16, "{pairs:?}");
        assert!(pairs.iter().all(|pair| pair.3), "{pairs:?}");
        // Cached: the second lookup does not recompute.
        let cached = generator.numerators.leg_exchanges.read().unwrap().len();
        assert_eq!(cached, 16);
        identical_pairs(&generator);
        assert_eq!(generator.numerators.leg_exchanges.read().unwrap().len(), 16);
    }

    #[test]
    fn non_identical_or_out_of_range_legs_are_not_certified() {
        let generator = Generator::new(Model::from_json(&qcd_definition().to_string()).unwrap());
        let rule = generator.model.vertex_rule_id("V_BBQQ").unwrap();
        assert!(!generator.leg_exchange_symmetric(rule, 0, 2).unwrap());
        assert!(!generator.leg_exchange_symmetric(rule, 1, 7).unwrap());
    }

    /// V_BBQQ with its third color x Lorentz channel removed: exchanging the
    /// two B legs maps the second channel onto the missing one.
    fn asymmetric_bbqq() -> Model {
        let mut definition = qcd_definition();
        let rule = definition["vertex_rules"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|rule| rule["name"] == "V_BBQQ")
            .unwrap();
        rule["couplings"][2][2] = serde_json::Value::Null;
        Model::from_json(&definition.to_string()).unwrap()
    }

    #[test]
    fn asymmetric_four_point_rule_is_not_certified() {
        let generator = Generator::new(asymmetric_bbqq());
        let rule = generator.model.vertex_rule_id("V_BBQQ").unwrap();
        assert!(!generator.leg_exchange_symmetric(rule, 0, 1).unwrap());
        assert!(!generator.leg_exchange_symmetric(rule, 2, 3).unwrap());
        // The other rules are untouched.
        let other = generator.model.vertex_rule_id("V_QQQQ").unwrap();
        assert!(generator.leg_exchange_symmetric(other, 0, 3).unwrap());
    }

    #[test]
    fn sign_flipping_exchange_is_not_certified() {
        // f(1,2,a) f(3,4,a) g(1,2) g(3,4) is odd under B <-> B: an exchange
        // that flips the sign must not be accepted.
        let mut definition = qcd_definition();
        definition["lorentz_structures"]
            .as_array_mut()
            .unwrap()
            .push(serde_json::json!({
                "name": "SYM",
                "spins": [3, 3, 3, 3],
                "structure": "Metric(idx(1,1),idx(1,2))*Metric(idx(1,3),idx(1,4))"
            }));
        let rule = definition["vertex_rules"]
            .as_array_mut()
            .unwrap()
            .iter_mut()
            .find(|rule| rule["name"] == "V_BBQQ")
            .unwrap();
        rule["lorentz_structures"] = serde_json::json!(["SYM"]);
        rule["couplings"] = serde_json::json!([["gs2_i"], [null], [null]]);
        let generator = Generator::new(Model::from_json(&definition.to_string()).unwrap());
        let rule = generator.model.vertex_rule_id("V_BBQQ").unwrap();
        assert!(!generator.leg_exchange_symmetric(rule, 0, 1).unwrap());
    }

    #[test]
    fn transpositions_compose_to_the_permutation() {
        // (0 2 3) and the fixed point 1, and a 2-cycle (4 5).
        let permutation = BTreeMap::from([(0, 2), (1, 1), (2, 3), (3, 0), (4, 5), (5, 4)]);
        let pairs = transpositions(&permutation);
        assert_eq!(pairs.len(), 3);
        // Applying the transpositions in order maps every slot to its image.
        for (slot, image) in &permutation {
            let mapped = pairs.iter().fold(*slot, |current, (a, b)| {
                if current == *a {
                    *b
                } else if current == *b {
                    *a
                } else {
                    current
                }
            });
            assert_eq!(mapped, *image, "slot {slot}");
        }
    }
}
