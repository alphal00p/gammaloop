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

use feynkit_graph::{DiagramEndpoint, DiagramHalfEdge, EdgeId, FeynmanDiagram, NumeratorInterface};
use idenso::{CookMode, CookSettings, dirac::AGS, representations::Bispinor};
use spenso::{network::library::symbolic::ETS, structure::representation::RepName};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    graph::Graph,
};

use super::{
    ColoredNode, ColoredTopology, EdgeColor, ExternalNode, FermionSigns, GenerationError,
    GenerationOptions, GenerationType, Generator, ResolvedProcess, grouping,
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
            if mask == 0 && !same_structure(&variant, diagram) {
                return skip("the rebuilt diagram does not reproduce the input frame");
            }
            // Every variant must keep the input's vertices and edge endpoints,
            // so `Q(edge)` and the half-edge transport mean the same thing.
            if mask != 0 && !same_endpoints(&variant, diagram) {
                return skip("a reversed variant moved an edge endpoint");
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

/// Same vertices, edges, endpoints, species and attached interaction slots.
fn same_structure(left: &FeynmanDiagram, right: &FeynmanDiagram) -> bool {
    let edges = |diagram: &FeynmanDiagram| {
        diagram
            .edges()
            .map(|(id, ends, edge)| {
                (
                    id,
                    ends.source,
                    ends.target,
                    edge.particle,
                    ends.source.map(|_| edge.source_slot()),
                    ends.target.map(|_| edge.target_slot()),
                )
            })
            .collect::<Vec<_>>()
    };
    left.vertices()
        .map(|(id, vertex)| (id, vertex.interaction))
        .eq(right
            .vertices()
            .map(|(id, vertex)| (id, vertex.interaction)))
        && edges(left) == edges(right)
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
