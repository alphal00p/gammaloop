use std::collections::{BTreeMap, BTreeSet};

use feynkit_graph::{CanonicalCutSide, ExternalLeg, FeynmanDiagram};
use feynkit_model::{Model, ModelError, ParameterId, ParticleId};
#[cfg(test)]
use symbolica::parser::ParseSettings;
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    graph::Graph,
    id::Replacement,
    symbol,
};
use thiserror::Error;

use crate::{DiagramGroup, GraphGroupingOptions, GroupMember, NumeratorGrouping};

#[derive(Debug, Error)]
pub enum GroupingError {
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error("failed to construct a topology grouping key: {0}")]
    TopologyConstruction(String),
    #[error("source diagram index {0} does not fit in a symbolic integer")]
    SourceDiagramIndexOverflow(usize),
}

pub(crate) struct GroupingOutcome {
    pub diagrams: Vec<FeynmanDiagram>,
    pub groups: Vec<DiagramGroup>,
    pub zero_numerator_count: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum TopologyVertex {
    Internal,
    External(ExternalLeg),
    CutSide {
        side: CanonicalCutSide,
        coupling_orders: BTreeMap<String, usize>,
        loop_count: usize,
    },
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum InternalParticleColor {
    MassAndSpin { mass: ParameterId, spin: i64 },
    Exact(ParticleId),
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum TopologyEdge {
    Internal(InternalParticleColor),
    External {
        particle: ParticleId,
        directed: bool,
    },
    CutPair,
    CutMembership,
}

type TopologyKey = Graph<TopologyVertex, TopologyEdge>;

struct PreparedNumerator {
    exact: Atom,
    samples: Vec<Atom>,
}

struct CanonicalTopology {
    key: TopologyKey,
}

#[derive(Clone, Copy)]
enum ComparisonMode {
    Identical,
    UpToSign,
    UpToScalar,
}

pub(crate) fn group_diagrams(
    diagrams: Vec<FeynmanDiagram>,
    model: &Model,
    grouping: &NumeratorGrouping,
    symmetrize_left_right: bool,
) -> Result<GroupingOutcome, GroupingError> {
    if matches!(grouping, NumeratorGrouping::None) {
        let groups = singleton_groups(&diagrams, &(0..diagrams.len()).collect::<Vec<_>>());
        return Ok(GroupingOutcome {
            diagrams,
            groups,
            zero_numerator_count: 0,
        });
    }

    let scalar_names = scalar_names(model);
    let mut retained = Vec::with_capacity(diagrams.len());
    let mut source_diagrams = Vec::with_capacity(diagrams.len());
    let mut prepared = Vec::with_capacity(diagrams.len());
    let mut zero_numerator_count = 0;
    for (source_diagram, diagram) in diagrams.into_iter().enumerate() {
        // External wavefunctions carry diagram-local edge IDs. The topology
        // key below already compares exact external particles and states, so
        // including those IDs would split physically identical numerators.
        // Keep the projector in the full zero check and on the retained master,
        // but compare only the generated numerator and its common prefactor.
        let exact = diagram.numerator() * diagram.numerator_prefactor();
        if (&exact * diagram.projector()).expand().is_zero() {
            zero_numerator_count += 1;
            continue;
        }
        retained.push(diagram);
        source_diagrams.push(source_diagram);
        prepared.push(PreparedNumerator {
            exact,
            samples: Vec::new(),
        });
    }

    let (options, mode) = match grouping {
        NumeratorGrouping::Identical(options) => (options, ComparisonMode::Identical),
        NumeratorGrouping::UpToSign(options) => (options, ComparisonMode::UpToSign),
        NumeratorGrouping::UpToScalar(options) => (options, ComparisonMode::UpToScalar),
        NumeratorGrouping::None | NumeratorGrouping::OnlyDetectZeroes => {
            return Ok(GroupingOutcome {
                groups: singleton_groups(&retained, &source_diagrams),
                diagrams: retained,
                zero_numerator_count,
            });
        }
    };

    let mut buckets = BTreeMap::<TopologyKey, Vec<usize>>::new();
    for (index, (diagram, numerator)) in retained.iter().zip(&mut prepared).enumerate() {
        let key = topology_key(diagram, model, options, symmetrize_left_right)?;
        numerator.samples = numerical_samples(&numerator.exact, options, &scalar_names);
        buckets.entry(key).or_default().push(index);
    }
    let mut groups = Vec::new();
    for indices in buckets.into_values() {
        let mut bucket_groups: Vec<DiagramGroup> = Vec::new();
        for diagram in indices {
            let matching_group =
                bucket_groups
                    .iter()
                    .enumerate()
                    .find_map(|(group_index, group)| {
                        compare(
                            &prepared[diagram],
                            &prepared[group.master],
                            mode,
                            options,
                            &scalar_names,
                        )
                        .map(|ratio| (group_index, ratio))
                    });
            if let Some((group_index, ratio)) = matching_group {
                bucket_groups[group_index].members.push(GroupMember {
                    source_diagram: source_diagrams[diagram],
                    source_id: retained[diagram].id(),
                    source_name: retained[diagram].name().to_owned(),
                    diagram,
                    ratio,
                    overall_factor: retained[diagram].overall_factor().clone(),
                });
            } else {
                bucket_groups.push(DiagramGroup {
                    master: diagram,
                    members: vec![GroupMember {
                        source_diagram: source_diagrams[diagram],
                        source_id: retained[diagram].id(),
                        source_name: retained[diagram].name().to_owned(),
                        diagram,
                        ratio: Atom::one(),
                        overall_factor: retained[diagram].overall_factor().clone(),
                    }],
                });
            }
        }
        groups.extend(bucket_groups);
    }
    groups.sort_by_key(|group| group.master);

    let (diagrams, groups) = collapse_groups(retained, groups)?;
    Ok(GroupingOutcome {
        diagrams,
        groups,
        zero_numerator_count,
    })
}

fn collapse_groups(
    diagrams: Vec<FeynmanDiagram>,
    groups: Vec<DiagramGroup>,
) -> Result<(Vec<FeynmanDiagram>, Vec<DiagramGroup>), GroupingError> {
    let mut masters = Vec::with_capacity(groups.len());
    let mut collapsed_groups = Vec::with_capacity(groups.len());
    for group in groups {
        let output_index = masters.len();
        let master = diagrams[group.master].clone();
        let factor = if group.members.len() == 1 {
            master.overall_factor().clone()
        } else {
            let mut factor = Atom::Zero;
            for member in &group.members {
                let source_diagram = i64::try_from(member.source_diagram).map_err(|_| {
                    GroupingError::SourceDiagramIndexOverflow(member.source_diagram)
                })?;
                factor += function!(
                    symbol!("NumeratorDependentGrouping"),
                    Atom::num(source_diagram),
                    member.ratio.clone(),
                    member.overall_factor.clone()
                );
            }
            factor
        };
        masters.push(master.with_overall_factor(factor));
        collapsed_groups.push(DiagramGroup {
            master: output_index,
            members: group
                .members
                .into_iter()
                .map(|mut member| {
                    member.diagram = output_index;
                    member
                })
                .collect(),
        });
    }
    Ok((masters, collapsed_groups))
}

fn singleton_groups(diagrams: &[FeynmanDiagram], source_diagrams: &[usize]) -> Vec<DiagramGroup> {
    diagrams
        .iter()
        .zip(source_diagrams.iter().copied())
        .enumerate()
        .map(|(diagram, (source, source_diagram))| DiagramGroup {
            master: diagram,
            members: vec![GroupMember {
                source_diagram,
                source_id: source.id(),
                source_name: source.name().to_owned(),
                diagram,
                ratio: Atom::one(),
                overall_factor: source.overall_factor().clone(),
            }],
        })
        .collect()
}

fn scalar_names(model: &Model) -> BTreeSet<String> {
    model
        .parameters()
        .iter()
        .map(|parameter| parameter.name.clone())
        .chain(
            model
                .couplings()
                .iter()
                .map(|coupling| coupling.name.clone()),
        )
        .chain(
            model
                .functions()
                .iter()
                .map(|function| function.name.clone()),
        )
        .chain(
            model
                .form_factors()
                .iter()
                .map(|form_factor| form_factor.name.clone()),
        )
        .collect()
}

fn numerical_samples(
    numerator: &Atom,
    options: &GraphGroupingOptions,
    scalar_names: &BTreeSet<String>,
) -> Vec<Atom> {
    let mut indeterminates: Vec<_> = numerator
        .get_all_indeterminates(false)
        .into_iter()
        .map(|indeterminate| indeterminate.to_owned())
        .filter(|indeterminate| {
            options.fully_numerical_substitution
                || !expression_is_scalar(indeterminate, scalar_names)
        })
        .collect();
    indeterminates.sort_by_key(|indeterminate| indeterminate.to_canonical_string());

    (0..options.number_of_numerical_samples)
        .map(|sample| {
            let replacements: Vec<_> = indeterminates
                .iter()
                .map(|indeterminate| {
                    Replacement::new(
                        indeterminate.clone(),
                        Atom::num(sample_value(
                            &indeterminate.to_canonical_string(),
                            options.numerical_sample_seed,
                            sample,
                        )),
                    )
                })
                .collect();
            numerator.replace_multiple(&replacements).expand()
        })
        .collect()
}

fn sample_value(expression: &str, seed: u16, sample: usize) -> i64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64
        ^ u64::from(seed)
        ^ (sample as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    for byte in expression.bytes() {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    let magnitude = (hash % 29 + 2) as i64;
    if hash & (1 << 63) == 0 {
        magnitude
    } else {
        -magnitude
    }
}

fn compare(
    candidate: &PreparedNumerator,
    master: &PreparedNumerator,
    mode: ComparisonMode,
    options: &GraphGroupingOptions,
    scalar_names: &BTreeSet<String>,
) -> Option<Atom> {
    if options.check_canonical_numerator
        && let Some(ratio) = compare_atoms(&candidate.exact, &master.exact, mode, scalar_names)
    {
        return Some(ratio);
    }
    let mut ratios = candidate
        .samples
        .iter()
        .zip(&master.samples)
        .map(|(candidate, master)| compare_atoms(candidate, master, mode, scalar_names));
    let first = ratios.next()??;
    if ratios.all(|ratio| ratio.is_some_and(|ratio| expressions_equal(&ratio, &first))) {
        Some(first)
    } else {
        None
    }
}

fn compare_atoms(
    candidate: &Atom,
    master: &Atom,
    mode: ComparisonMode,
    scalar_names: &BTreeSet<String>,
) -> Option<Atom> {
    if expressions_equal(candidate, master) {
        return Some(Atom::num(1));
    }
    if matches!(mode, ComparisonMode::UpToSign) && (candidate + master).expand().is_zero() {
        return Some(Atom::num(-1));
    }
    if !matches!(mode, ComparisonMode::UpToScalar) || master.is_zero() {
        return None;
    }
    let ratio = (candidate / master).cancel();
    let reconstructed = &ratio * master;
    if expression_is_scalar(&ratio, scalar_names) && (candidate - &reconstructed).expand().is_zero()
    {
        Some(ratio)
    } else {
        None
    }
}

fn expressions_equal(left: &Atom, right: &Atom) -> bool {
    (left - right).expand().is_zero()
}

fn expression_is_scalar(expression: &Atom, scalar_names: &BTreeSet<String>) -> bool {
    expression
        .get_all_symbols(true)
        .iter()
        .all(|symbol| symbol.is_scalar() || scalar_names.contains(symbol.get_stripped_name()))
}

fn topology_key(
    diagram: &FeynmanDiagram,
    model: &Model,
    options: &GraphGroupingOptions,
    symmetrize_left_right: bool,
) -> Result<TopologyKey, GroupingError> {
    canonical_topologies(diagram, model, options, symmetrize_left_right)?
        .into_iter()
        .map(|topology| topology.key)
        .min()
        .ok_or_else(|| GroupingError::TopologyConstruction("no canonical topology".to_owned()))
}

fn canonical_topologies(
    diagram: &FeynmanDiagram,
    model: &Model,
    options: &GraphGroupingOptions,
    symmetrize_left_right: bool,
) -> Result<Vec<CanonicalTopology>, GroupingError> {
    let mut topologies = vec![topology_key_variant(diagram, model, options, None)?];
    if symmetrize_left_right && let Some(partners) = left_right_partners(diagram) {
        topologies.push(topology_key_variant(
            diagram,
            model,
            options,
            Some(&partners),
        )?);
    }
    Ok(topologies)
}

fn left_right_partners(diagram: &FeynmanDiagram) -> Option<BTreeMap<ExternalLeg, ExternalLeg>> {
    let legs = diagram
        .vertices()
        .filter_map(|(_, vertex)| vertex.external.clone())
        .collect::<Vec<_>>();
    let by_connection = legs
        .iter()
        .cloned()
        .map(|leg| ((leg.connection, leg.state), leg))
        .collect::<BTreeMap<_, _>>();
    if by_connection.len() != legs.len() {
        return None;
    }
    legs.into_iter()
        .map(|leg| {
            let opposite = match leg.state {
                feynkit_graph::ExternalState::Incoming => feynkit_graph::ExternalState::Outgoing,
                feynkit_graph::ExternalState::Outgoing => feynkit_graph::ExternalState::Incoming,
            };
            by_connection
                .get(&(leg.connection, opposite))
                .cloned()
                .map(|partner| (leg, partner))
        })
        .collect()
}

fn topology_key_variant(
    diagram: &FeynmanDiagram,
    model: &Model,
    options: &GraphGroupingOptions,
    left_right_partners: Option<&BTreeMap<ExternalLeg, ExternalLeg>>,
) -> Result<CanonicalTopology, GroupingError> {
    let mut graph = Graph::new();
    for (_, vertex) in diagram.vertices() {
        let external = vertex.external.clone().map(|leg| {
            left_right_partners
                .and_then(|partners| partners.get(&leg))
                .cloned()
                .unwrap_or(leg)
        });
        graph.add_node(external.map_or(TopologyVertex::Internal, TopologyVertex::External));
    }
    for (_, endpoints, edge) in diagram.edges() {
        let external = diagram
            .vertex(endpoints.source)
            .is_some_and(|vertex| vertex.is_external())
            || diagram
                .vertex(endpoints.target)
                .is_some_and(|vertex| vertex.is_external());
        let (directed, color) = if external {
            // External edges retain both their exact particle species and
            // direction: exchanging different mass-degenerate external states
            // would change the process rather than merely its topology.
            (
                edge.directed,
                TopologyEdge::External {
                    particle: edge.particle,
                    directed: edge.directed,
                },
            )
        } else {
            let particle = model.particle_by_id(edge.particle)?;
            // Spin remains part of the reduced internal color so massless
            // fermions and vectors cannot be interchanged by canonicalization.
            let color = if options.differentiate_particle_masses_only {
                InternalParticleColor::MassAndSpin {
                    mass: particle.mass,
                    spin: particle.spin,
                }
            } else {
                InternalParticleColor::Exact(edge.particle)
            };
            (false, TopologyEdge::Internal(color))
        };
        let (source, target) = if left_right_partners.is_some() && external && directed {
            (endpoints.target.0, endpoints.source.0)
        } else {
            (endpoints.source.0, endpoints.target.0)
        };
        graph
            .add_edge(source, target, directed, color)
            .map_err(|error| GroupingError::TopologyConstruction(error.to_string()))?;
    }
    let endpoints = diagram
        .edges()
        .map(|(edge, endpoints, _)| (edge, endpoints))
        .collect::<BTreeMap<_, _>>();
    for cut in diagram.cuts() {
        let left = graph.add_node(TopologyVertex::CutSide {
            side: if left_right_partners.is_some() {
                CanonicalCutSide::Right
            } else {
                CanonicalCutSide::Left
            },
            coupling_orders: cut.left.coupling_orders.clone(),
            loop_count: cut.left.loop_count,
        });
        let right = graph.add_node(TopologyVertex::CutSide {
            side: if left_right_partners.is_some() {
                CanonicalCutSide::Left
            } else {
                CanonicalCutSide::Right
            },
            coupling_orders: cut.right.coupling_orders.clone(),
            loop_count: cut.right.loop_count,
        });
        graph
            .add_edge(left, right, false, TopologyEdge::CutPair)
            .map_err(|error| GroupingError::TopologyConstruction(error.to_string()))?;
        for (side, half_edges) in [(left, &cut.left.half_edges), (right, &cut.right.half_edges)] {
            let vertices = half_edges
                .iter()
                .filter_map(|half_edge| {
                    endpoints
                        .get(&half_edge.edge)
                        .map(|endpoints| match half_edge.endpoint {
                            feynkit_graph::DiagramEndpoint::Source => endpoints.source,
                            feynkit_graph::DiagramEndpoint::Target => endpoints.target,
                        })
                })
                .collect::<BTreeSet<_>>();
            for vertex in vertices {
                graph
                    .add_edge(side, vertex.0, false, TopologyEdge::CutMembership)
                    .map_err(|error| GroupingError::TopologyConstruction(error.to_string()))?;
            }
        }
    }
    Ok(CanonicalTopology {
        key: graph.canonize().graph,
    })
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use feynkit_graph::{
        DiagramCut, DiagramCutSide, DiagramEdge, DiagramEndpoint, DiagramHalfEdge, DiagramVertex,
        ExternalState, VertexId,
    };

    use super::*;

    fn test_atom(expression: &str) -> Atom {
        Atom::parse(
            expression,
            "feynkit_grouping_test",
            ParseSettings::default(),
        )
        .unwrap()
    }

    fn model() -> Model {
        Model::from_json(
            r#"{
                "name":"grouping",
                "restriction":null,
                "orders":[],
                "parameters":[
                    {"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null},
                    {"name":"M","lhablock":null,"lhacode":null,"nature":"external","parameter_type":"real","value":[1.0,0.0],"expression":null},
                    {"name":"N","lhablock":null,"lhacode":null,"nature":"external","parameter_type":"real","value":[2.0,0.0],"expression":null}
                ],
                "particles":[
                    {"pdg_code":25,"name":"phi","antiname":"phi","spin":1,"color":1,"mass":"M","width":"ZERO","texname":"phi","antitexname":"phi","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":35,"name":"chi","antiname":"chi","spin":1,"color":1,"mass":"M","width":"ZERO","texname":"chi","antitexname":"chi","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0},
                    {"pdg_code":45,"name":"eta","antiname":"eta","spin":1,"color":1,"mass":"N","width":"ZERO","texname":"eta","antitexname":"eta","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0}
                ],
                "propagators":[],
                "lorentz_structures":[
                    {"name":"L","spins":[1,1],"structure":"1"}
                ],
                "couplings":[
                    {"name":"GC","expression":"1","orders":[],"value":null}
                ],
                "vertex_rules":[
                    {"name":"V","particles":["phi","phi"],"color_structures":["1"],"lorentz_structures":["L"],"couplings":[["GC"]]}
                ],
                "functions":[],
                "form_factors":[]
            }"#,
        )
        .unwrap()
    }

    fn particle(model: &Model, pdg: i64) -> ParticleId {
        model.particle_id_by_pdg(pdg).unwrap()
    }

    fn diagram(
        model: &Arc<Model>,
        name: &str,
        numerator: &str,
        internal_pdg: i64,
        external_pdg: i64,
    ) -> FeynmanDiagram {
        diagram_with_projector(model, name, numerator, "1", internal_pdg, external_pdg)
    }

    fn diagram_with_projector(
        model: &Arc<Model>,
        name: &str,
        numerator: &str,
        projector: &str,
        internal_pdg: i64,
        external_pdg: i64,
    ) -> FeynmanDiagram {
        let numerator =
            Atom::parse(numerator, "feynkit_grouping_test", ParseSettings::default()).unwrap();
        let projector =
            Atom::parse(projector, "feynkit_grouping_test", ParseSettings::default()).unwrap();
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), name)
            .numerator(numerator)
            .projector(projector);
        builder.add_vertex(DiagramVertex::external("in", 0, ExternalState::Incoming));
        let interaction = model.vertex_rule_id("V").unwrap();
        builder.add_vertex(DiagramVertex::interaction("left", interaction));
        builder.add_vertex(DiagramVertex::interaction("right", interaction));
        builder.add_vertex(DiagramVertex::external("out", 1, ExternalState::Outgoing));
        builder
            .add_edge(
                VertexId(0),
                VertexId(1),
                DiagramEdge::new(particle(model, external_pdg), false),
            )
            .unwrap();
        builder
            .add_edge(
                VertexId(1),
                VertexId(2),
                DiagramEdge::new(particle(model, internal_pdg), false),
            )
            .unwrap();
        builder
            .add_edge(
                VertexId(2),
                VertexId(3),
                DiagramEdge::new(particle(model, external_pdg), false),
            )
            .unwrap();
        builder.build().unwrap()
    }

    fn cut_diagram(model: &Arc<Model>, mirrored: bool) -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), "cut-line");
        let incoming = builder.add_vertex(DiagramVertex::external_in_connection(
            "in",
            0,
            ExternalState::Incoming,
            0,
        ));
        let interaction = builder.add_vertex(DiagramVertex::interaction(
            "interaction",
            model.vertex_rule_id("V").unwrap(),
        ));
        let outgoing = builder.add_vertex(DiagramVertex::external_in_connection(
            "out",
            1,
            ExternalState::Outgoing,
            0,
        ));
        let scalar = || DiagramEdge::new(particle(model, 25), false);
        let incoming_edge = builder.add_edge(incoming, interaction, scalar()).unwrap();
        let outgoing_edge = builder.add_edge(interaction, outgoing, scalar()).unwrap();
        let incoming_source = DiagramHalfEdge {
            edge: incoming_edge,
            endpoint: DiagramEndpoint::Source,
        };
        let incoming_target = DiagramHalfEdge {
            edge: incoming_edge,
            endpoint: DiagramEndpoint::Target,
        };
        let outgoing_source = DiagramHalfEdge {
            edge: outgoing_edge,
            endpoint: DiagramEndpoint::Source,
        };
        let outgoing_target = DiagramHalfEdge {
            edge: outgoing_edge,
            endpoint: DiagramEndpoint::Target,
        };
        let side = |half_edges| DiagramCutSide {
            half_edges,
            coupling_orders: BTreeMap::new(),
            loop_count: 0,
        };
        let mut cut = DiagramCut {
            cut: vec![incoming_source],
            left: side(vec![incoming_source]),
            right: side(vec![incoming_target, outgoing_source, outgoing_target]),
        };
        if mirrored {
            // The global transformation exchanges both physical sides at
            // once. This is intentionally not just an independent relabeling
            // of the cut nodes.
            cut.cut = vec![outgoing_source];
            cut.left = side(vec![incoming_source, incoming_target, outgoing_source]);
            cut.right = side(vec![outgoing_target]);
        }
        let diagram = builder.cuts(vec![cut]).build().unwrap();
        diagram.validate().unwrap();
        diagram
    }

    #[derive(Clone, Copy)]
    enum LeftRightAttachment {
        Identity,
        GlobalSwap,
        PartialSwap,
    }

    fn left_right_diagram(model: &Arc<Model>, attachment: LeftRightAttachment) -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder(Arc::clone(model), "left-right");
        let external = [
            builder.add_vertex(DiagramVertex::external_in_connection(
                "in0",
                0,
                ExternalState::Incoming,
                0,
            )),
            builder.add_vertex(DiagramVertex::external_in_connection(
                "in1",
                1,
                ExternalState::Incoming,
                1,
            )),
            builder.add_vertex(DiagramVertex::external_in_connection(
                "out0",
                2,
                ExternalState::Outgoing,
                0,
            )),
            builder.add_vertex(DiagramVertex::external_in_connection(
                "out1",
                3,
                ExternalState::Outgoing,
                1,
            )),
        ];
        let interaction = model.vertex_rule_id("V").unwrap();
        let left = builder.add_vertex(DiagramVertex::interaction("left", interaction));
        let right = builder.add_vertex(DiagramVertex::interaction("right", interaction));
        let marker = builder.add_vertex(DiagramVertex::interaction("marker", interaction));
        let targets = match attachment {
            LeftRightAttachment::Identity => [left, left, right, right],
            LeftRightAttachment::GlobalSwap => [right, right, left, left],
            LeftRightAttachment::PartialSwap => [right, left, left, right],
        };
        for (index, (&external, &target)) in external.iter().zip(&targets).enumerate() {
            let (source, target) = if index < 2 {
                (external, target)
            } else {
                (target, external)
            };
            builder
                .add_edge(source, target, DiagramEdge::new(particle(model, 25), false))
                .unwrap();
        }
        builder
            .add_edge(left, right, DiagramEdge::new(particle(model, 25), false))
            .unwrap();
        builder
            .add_edge(left, marker, DiagramEdge::new(particle(model, 25), false))
            .unwrap();
        builder.build().unwrap()
    }

    fn exact_options() -> GraphGroupingOptions {
        GraphGroupingOptions {
            check_canonical_numerator: true,
            ..GraphGroupingOptions::default()
        }
    }

    #[test]
    fn groups_identical_numerators() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x+y", 25, 25),
                diagram(&model, "g1", "y+x", 25, 25),
            ],
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            false,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].master, 0);
        assert_eq!(grouped.groups[0].members[1].ratio, test_atom("1"));
    }

    #[test]
    fn diagram_local_external_projector_ids_do_not_split_numerator_groups() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram_with_projector(&model, "g0", "x+y", "external_state(0)", 25, 25),
                diagram_with_projector(&model, "g1", "y+x", "external_state(4)", 25, 25),
            ],
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            false,
        )
        .unwrap();

        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members.len(), 2);
        let expected = Atom::parse(
            "external_state(0)",
            "feynkit_grouping_test",
            ParseSettings::default(),
        )
        .unwrap();
        assert_eq!(grouped.diagrams[0].projector(), &expected);
    }

    #[test]
    fn groups_numerators_up_to_sign() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x+y", 25, 25),
                diagram(&model, "g1", "-x-y", 25, 25),
            ],
            &model,
            &NumeratorGrouping::UpToSign(exact_options()),
            false,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members[1].ratio, test_atom("-1"));
    }

    #[test]
    fn scalar_ratio_is_member_over_master() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x+y", 25, 25),
                diagram(&model, "g1", "2*x+2*y", 25, 25),
            ],
            &model,
            &NumeratorGrouping::UpToScalar(exact_options()),
            false,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members[1].ratio, test_atom("2"));
    }

    #[test]
    fn topology_uses_internal_mass_and_spin_but_exact_external_species() {
        let model = Arc::new(model());
        let diagrams = vec![
            diagram(&model, "g0", "1", 25, 25),
            diagram(&model, "g1", "1", 35, 25),
            diagram(&model, "g2", "1", 45, 25),
            diagram(&model, "g3", "1", 25, 35),
        ];
        let grouped = group_diagrams(
            diagrams.clone(),
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            false,
        )
        .unwrap();
        assert_eq!(
            grouped
                .groups
                .iter()
                .map(|group| {
                    group
                        .members
                        .iter()
                        .map(|member| member.source_diagram)
                        .collect()
                })
                .collect::<Vec<Vec<_>>>(),
            vec![vec![0, 1], vec![2], vec![3]]
        );

        let exact_species = GraphGroupingOptions {
            differentiate_particle_masses_only: false,
            ..exact_options()
        };
        let grouped = group_diagrams(
            diagrams,
            &model,
            &NumeratorGrouping::Identical(exact_species),
            false,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 4);
    }

    #[test]
    fn mirrored_cut_partitions_group_only_with_left_right_symmetry() {
        let model = Arc::new(model());
        let diagrams = vec![cut_diagram(&model, false), cut_diagram(&model, true)];
        let grouped = group_diagrams(
            diagrams.clone(),
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            false,
        )
        .unwrap();

        assert_eq!(grouped.groups.len(), 2);
        assert!(grouped.groups.iter().all(|group| group.members.len() == 1));

        let grouped = group_diagrams(
            diagrams,
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            true,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members.len(), 2);
    }

    #[test]
    fn only_one_global_left_right_swap_is_in_the_grouping_orbit() {
        let model = Arc::new(model());
        let identity = left_right_diagram(&model, LeftRightAttachment::Identity);
        let global = left_right_diagram(&model, LeftRightAttachment::GlobalSwap);
        let partial = left_right_diagram(&model, LeftRightAttachment::PartialSwap);

        let grouped = group_diagrams(
            vec![identity.clone(), global],
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            true,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members.len(), 2);

        let grouped = group_diagrams(
            vec![identity, partial],
            &model,
            &NumeratorGrouping::Identical(exact_options()),
            true,
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 2);
    }

    #[test]
    fn removes_only_proven_zero_numerators_and_preserves_source_ids() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x-x", 25, 25),
                diagram(&model, "g1", "1", 25, 25),
            ],
            &model,
            &NumeratorGrouping::OnlyDetectZeroes,
            false,
        )
        .unwrap();
        assert_eq!(grouped.zero_numerator_count, 1);
        assert_eq!(grouped.diagrams.len(), 1);
        assert_eq!(grouped.groups[0].members[0].diagram, 0);
        assert_eq!(grouped.groups[0].members[0].source_diagram, 1);
        assert_eq!(grouped.diagrams[0].name(), "g1");
    }

    #[test]
    fn zero_projectors_are_removed_before_grouping() {
        let model = Arc::new(model());
        let grouped = group_diagrams(
            vec![diagram_with_projector(&model, "zero", "1", "0", 25, 25)],
            &model,
            &NumeratorGrouping::OnlyDetectZeroes,
            false,
        )
        .unwrap();

        assert_eq!(grouped.zero_numerator_count, 1);
        assert!(grouped.diagrams.is_empty());
        assert!(grouped.groups.is_empty());
    }
}
