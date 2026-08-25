use std::collections::{BTreeMap, BTreeSet};

use feynkit_graph::{ExternalLeg, FeynmanDiagram, ParticleReference};
use feynkit_model::{Model, ModelError};
use symbolica::{
    atom::{Atom, AtomCore},
    graph::Graph,
    id::Replacement,
    parser::ParseSettings,
};
use thiserror::Error;

use crate::{DiagramGroup, GraphGroupingOptions, GroupMember, NumeratorGrouping};

#[derive(Debug, Error)]
pub enum GroupingError {
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error(
        "failed to parse numerator of diagram '{diagram}': {message}; expression: {expression}"
    )]
    NumeratorParse {
        diagram: String,
        expression: String,
        message: String,
    },
    #[error("failed to construct a topology grouping key: {0}")]
    TopologyConstruction(String),
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
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum InternalParticleColor {
    MassAndSpin { mass: String, spin: i64 },
    Exact(ParticleReference),
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum TopologyEdge {
    Internal(InternalParticleColor),
    External {
        particle: ParticleReference,
        directed: bool,
    },
}

type TopologyKey = Graph<TopologyVertex, TopologyEdge>;

struct PreparedNumerator {
    exact: Atom,
    samples: Vec<Atom>,
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
) -> Result<GroupingOutcome, GroupingError> {
    if matches!(grouping, NumeratorGrouping::None) {
        let groups = singleton_groups(&(0..diagrams.len()).collect::<Vec<_>>());
        return Ok(GroupingOutcome {
            diagrams,
            groups,
            zero_numerator_count: 0,
        });
    }

    let scalar_names = scalar_names(model);
    let options = grouping_options(grouping);
    let mut retained = Vec::with_capacity(diagrams.len());
    let mut source_diagrams = Vec::with_capacity(diagrams.len());
    let mut prepared = Vec::with_capacity(diagrams.len());
    let mut zero_numerator_count = 0;
    for (source_diagram, diagram) in diagrams.into_iter().enumerate() {
        let exact = parse_numerator(&diagram)?;
        if exact.expand().is_zero() {
            zero_numerator_count += 1;
            continue;
        }
        let samples = options.map_or_else(Vec::new, |options| {
            numerical_samples(&exact, options, &scalar_names)
        });
        retained.push(diagram);
        source_diagrams.push(source_diagram);
        prepared.push(PreparedNumerator { exact, samples });
    }

    let (options, mode) = match grouping {
        NumeratorGrouping::Identical(options) => (options, ComparisonMode::Identical),
        NumeratorGrouping::UpToSign(options) => (options, ComparisonMode::UpToSign),
        NumeratorGrouping::UpToScalar(options) => (options, ComparisonMode::UpToScalar),
        NumeratorGrouping::None | NumeratorGrouping::OnlyDetectZeroes => {
            return Ok(GroupingOutcome {
                groups: singleton_groups(&source_diagrams),
                diagrams: retained,
                zero_numerator_count,
            });
        }
    };

    let mut buckets = BTreeMap::<TopologyKey, Vec<usize>>::new();
    for (index, diagram) in retained.iter().enumerate() {
        buckets
            .entry(topology_key(diagram, model, options)?)
            .or_default()
            .push(index);
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
                    diagram,
                    ratio,
                });
            } else {
                bucket_groups.push(DiagramGroup {
                    master: diagram,
                    members: vec![GroupMember {
                        source_diagram: source_diagrams[diagram],
                        diagram,
                        ratio: "1".to_owned(),
                    }],
                });
            }
        }
        groups.extend(bucket_groups);
    }
    groups.sort_by_key(|group| group.master);

    Ok(GroupingOutcome {
        diagrams: retained,
        groups,
        zero_numerator_count,
    })
}

fn grouping_options(grouping: &NumeratorGrouping) -> Option<&GraphGroupingOptions> {
    match grouping {
        NumeratorGrouping::Identical(options)
        | NumeratorGrouping::UpToSign(options)
        | NumeratorGrouping::UpToScalar(options) => Some(options),
        NumeratorGrouping::None | NumeratorGrouping::OnlyDetectZeroes => None,
    }
}

fn singleton_groups(source_diagrams: &[usize]) -> Vec<DiagramGroup> {
    source_diagrams
        .iter()
        .copied()
        .enumerate()
        .map(|(diagram, source_diagram)| DiagramGroup {
            master: diagram,
            members: vec![GroupMember {
                source_diagram,
                diagram,
                ratio: "1".to_owned(),
            }],
        })
        .collect()
}

fn parse_numerator(diagram: &FeynmanDiagram) -> Result<Atom, GroupingError> {
    let expression = diagram.numerator().unwrap_or("1");
    Atom::parse(expression, "feynkit_grouping", ParseSettings::default()).map_err(|message| {
        GroupingError::NumeratorParse {
            diagram: diagram.name().to_owned(),
            expression: expression.to_owned(),
            message,
        }
    })
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
) -> Option<String> {
    if options.check_canonical_numerator
        && let Some(ratio) = compare_atoms(&candidate.exact, &master.exact, mode, scalar_names)
    {
        return Some(ratio.to_canonical_string());
    }
    let mut ratios = candidate
        .samples
        .iter()
        .zip(&master.samples)
        .map(|(candidate, master)| compare_atoms(candidate, master, mode, scalar_names));
    let first = ratios.next()??;
    if ratios.all(|ratio| ratio.is_some_and(|ratio| expressions_equal(&ratio, &first))) {
        Some(first.to_canonical_string())
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
) -> Result<TopologyKey, GroupingError> {
    let mut graph = Graph::new();
    for (_, vertex) in diagram.vertices() {
        graph.add_node(
            vertex
                .external
                .clone()
                .map_or(TopologyVertex::Internal, TopologyVertex::External),
        );
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
                    particle: edge.particle.clone(),
                    directed: edge.directed,
                },
            )
        } else {
            let particle = model.particle_by_pdg(edge.particle.pdg)?;
            // Spin remains part of the reduced internal color so massless
            // fermions and vectors cannot be interchanged by canonicalization.
            let color = if options.differentiate_particle_masses_only {
                InternalParticleColor::MassAndSpin {
                    mass: particle.mass.clone(),
                    spin: particle.spin,
                }
            } else {
                InternalParticleColor::Exact(edge.particle.clone())
            };
            (false, TopologyEdge::Internal(color))
        };
        graph
            .add_edge(endpoints.source.0, endpoints.target.0, directed, color)
            .map_err(|error| GroupingError::TopologyConstruction(error.to_string()))?;
    }
    Ok(graph.canonize().graph)
}

#[cfg(test)]
mod tests {
    use feynkit_graph::{DiagramEdge, DiagramVertex, ExternalState, VertexId};

    use super::*;

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
                "lorentz_structures":[],
                "couplings":[],
                "vertex_rules":[],
                "functions":[],
                "form_factors":[]
            }"#,
        )
        .unwrap()
    }

    fn particle(model: &Model, pdg: i64) -> ParticleReference {
        let particle = model.particle_by_pdg(pdg).unwrap();
        ParticleReference::new(&particle.name, particle.pdg_code)
    }

    fn diagram(
        model: &Model,
        name: &str,
        numerator: &str,
        internal_pdg: i64,
        external_pdg: i64,
    ) -> FeynmanDiagram {
        let mut builder = FeynmanDiagram::builder(name).numerator(numerator);
        builder.add_vertex(DiagramVertex::external("in", 0, ExternalState::Incoming));
        builder.add_vertex(DiagramVertex::interaction("left", "V"));
        builder.add_vertex(DiagramVertex::interaction("right", "V"));
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

    fn exact_options() -> GraphGroupingOptions {
        GraphGroupingOptions {
            check_canonical_numerator: true,
            ..GraphGroupingOptions::default()
        }
    }

    #[test]
    fn groups_identical_numerators() {
        let model = model();
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x+y", 25, 25),
                diagram(&model, "g1", "y+x", 25, 25),
            ],
            &model,
            &NumeratorGrouping::Identical(exact_options()),
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].master, 0);
        assert_eq!(grouped.groups[0].members[1].ratio, "1");
    }

    #[test]
    fn groups_numerators_up_to_sign() {
        let model = model();
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x+y", 25, 25),
                diagram(&model, "g1", "-x-y", 25, 25),
            ],
            &model,
            &NumeratorGrouping::UpToSign(exact_options()),
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members[1].ratio, "-1");
    }

    #[test]
    fn scalar_ratio_is_member_over_master() {
        let model = model();
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x+y", 25, 25),
                diagram(&model, "g1", "2*x+2*y", 25, 25),
            ],
            &model,
            &NumeratorGrouping::UpToScalar(exact_options()),
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 1);
        assert_eq!(grouped.groups[0].members[1].ratio, "2");
    }

    #[test]
    fn topology_uses_internal_mass_and_spin_but_exact_external_species() {
        let model = model();
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
        )
        .unwrap();
        assert_eq!(
            grouped
                .groups
                .iter()
                .map(|group| group.members.iter().map(|member| member.diagram).collect())
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
        )
        .unwrap();
        assert_eq!(grouped.groups.len(), 4);
    }

    #[test]
    fn removes_only_proven_zero_numerators_and_preserves_source_ids() {
        let model = model();
        let grouped = group_diagrams(
            vec![
                diagram(&model, "g0", "x-x", 25, 25),
                diagram(&model, "g1", "1", 25, 25),
            ],
            &model,
            &NumeratorGrouping::OnlyDetectZeroes,
        )
        .unwrap();
        assert_eq!(grouped.zero_numerator_count, 1);
        assert_eq!(grouped.diagrams.len(), 1);
        assert_eq!(grouped.groups[0].members[0].diagram, 0);
        assert_eq!(grouped.groups[0].members[0].source_diagram, 1);
        assert_eq!(grouped.diagrams[0].name(), "g1");
    }

    #[test]
    fn numerator_parse_failures_are_typed() {
        let model = model();
        assert!(matches!(
            group_diagrams(
                vec![diagram(&model, "bad", "(", 25, 25)],
                &model,
                &NumeratorGrouping::OnlyDetectZeroes,
            ),
            Err(GroupingError::NumeratorParse { diagram, .. }) if diagram == "bad"
        ));
    }
}
