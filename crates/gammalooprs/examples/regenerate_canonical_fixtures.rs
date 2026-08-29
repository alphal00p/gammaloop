//! Deterministically regenerate the canonical FeynKit DOT fixtures used by GammaLoop.
//!
//! Run from the repository root:
//!
//! ```text
//! nix develop -c cargo run --locked -p gammalooprs --example regenerate_canonical_fixtures -- --write
//! nix develop -c cargo run --locked -p gammalooprs --example regenerate_canonical_fixtures -- --verify
//! ```

use std::{
    collections::{BTreeMap, BTreeSet},
    env, fs,
    path::{Path, PathBuf},
};

use feynkit_generator::{GenerationFilter, GenerationOptions, Generator, Process, VertexSelector};
use feynkit_graph::{EdgeId, ExternalState, FeynmanDiagram, VertexId};
use feynkit_model::{Model, ParameterCard};
use gammalooprs::{
    model::{InputParamCard, ModelGammaLoopExt},
    utils::F,
};
use symbolica::atom::Atom;

type DynError = Box<dyn std::error::Error>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TargetVertex {
    Interaction(&'static str),
    External(ExternalState),
}

#[derive(Clone, Copy, Debug)]
struct TargetEdge {
    id: usize,
    source: usize,
    target: usize,
    particle: &'static str,
    directed: bool,
}

#[derive(Clone, Debug)]
struct Target {
    name: &'static str,
    path: &'static str,
    vertices: Vec<TargetVertex>,
    edges: Vec<TargetEdge>,
    loop_edges: Option<Vec<usize>>,
    overall_factor: OverallFactor,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OverallFactor {
    Generated,
    Unit,
    NegativeUnit,
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum CandidateVertex {
    Interaction(String),
    External(ExternalState),
}

#[derive(Clone, Debug)]
struct CandidateEdge {
    id: EdgeId,
    source: usize,
    target: usize,
    particle: String,
    directed: bool,
}

#[derive(Clone, Debug)]
struct Candidate {
    vertices: Vec<CandidateVertex>,
    edges: Vec<CandidateEdge>,
}

#[derive(Clone, Debug)]
struct ExactMatch {
    edge_mapping: BTreeMap<usize, EdgeId>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mode {
    Write,
    Verify,
}

fn parse_mode() -> Result<Mode, DynError> {
    let arguments = env::args().skip(1).collect::<Vec<_>>();
    match arguments.as_slice() {
        [argument] if argument == "--write" => Ok(Mode::Write),
        [argument] if argument == "--verify" => Ok(Mode::Verify),
        _ => Err(
            "usage: regenerate_canonical_fixtures (--write | --verify) (run from repository root)"
                .into(),
        ),
    }
}

fn standard_model() -> Result<Model, DynError> {
    let mut model = Model::from_path("assets/models/json/sm/sm.json")?
        .with_restriction(Some("default".to_owned()))?;
    let mut card: InputParamCard<F<f64>> =
        ParameterCard::from_path("assets/models/json/sm/restrict_default.json")?.into();
    model.simplify(&mut card)?;
    Ok(model)
}

fn scalar_model() -> Result<Model, DynError> {
    Ok(Model::from_path(
        "assets/models/json/scalars/scalars_2p_3p.json",
    )?)
}

fn target_edge(
    id: usize,
    source: usize,
    target: usize,
    particle: &'static str,
    directed: bool,
) -> TargetEdge {
    TargetEdge {
        id,
        source,
        target,
        particle,
        directed,
    }
}

fn addbar_target() -> Target {
    Target {
        name: "addbar",
        path: "tests/resources/graphs/addbar.dot",
        vertices: vec![
            TargetVertex::Interaction("V_71"),
            TargetVertex::Interaction("V_71"),
            TargetVertex::Interaction("V_74"),
            TargetVertex::Interaction("V_74"),
            TargetVertex::External(ExternalState::Incoming),
            TargetVertex::External(ExternalState::Outgoing),
        ],
        edges: vec![
            target_edge(0, 0, 2, "d", true),
            target_edge(1, 3, 0, "d", true),
            target_edge(2, 2, 1, "d", true),
            target_edge(3, 1, 3, "d", true),
            target_edge(4, 2, 3, "g", false),
            target_edge(5, 1, 4, "a", false),
            target_edge(6, 0, 5, "a", false),
        ],
        loop_edges: Some(vec![0, 3]),
        overall_factor: OverallFactor::Generated,
    }
}

fn scalar_targets() -> Vec<Target> {
    use ExternalState::{Incoming, Outgoing};
    use TargetVertex::{External, Interaction};

    vec![
        Target {
            name: "bubble_dot_1t1b",
            path: "tests/resources/graphs/bubble_dot_1t1b.dot",
            vertices: vec![
                Interaction("V_3_SCALAR_111"),
                Interaction("V_2_SCALAR_11"),
                Interaction("V_2_SCALAR_11"),
                Interaction("V_3_SCALAR_111"),
                External(Incoming),
                External(Outgoing),
            ],
            edges: vec![
                target_edge(0, 4, 0, "scalar_1", false),
                target_edge(1, 0, 1, "scalar_1", false),
                target_edge(2, 1, 3, "scalar_1", false),
                target_edge(3, 0, 2, "scalar_1", false),
                target_edge(4, 2, 3, "scalar_1", false),
                target_edge(5, 3, 5, "scalar_1", false),
            ],
            loop_edges: None,
            overall_factor: OverallFactor::Unit,
        },
        Target {
            name: "bubble_dot_2t0b",
            path: "tests/resources/graphs/bubble_dot_2t0b.dot",
            vertices: vec![
                Interaction("V_3_SCALAR_111"),
                Interaction("V_2_SCALAR_11"),
                Interaction("V_2_SCALAR_11"),
                Interaction("V_3_SCALAR_111"),
                External(Incoming),
                External(Outgoing),
            ],
            edges: vec![
                target_edge(0, 4, 0, "scalar_1", false),
                target_edge(1, 0, 1, "scalar_1", false),
                target_edge(2, 1, 2, "scalar_1", false),
                target_edge(3, 2, 3, "scalar_1", false),
                target_edge(4, 3, 0, "scalar_1", false),
                target_edge(5, 3, 5, "scalar_1", false),
            ],
            loop_edges: None,
            overall_factor: OverallFactor::Unit,
        },
        Target {
            name: "dotted_bubble",
            path: "tests/resources/graphs/dotted_bubble.dot",
            vertices: vec![
                Interaction("V_3_SCALAR_111"),
                Interaction("V_2_SCALAR_11"),
                Interaction("V_3_SCALAR_111"),
                External(Incoming),
                External(Outgoing),
            ],
            edges: vec![
                target_edge(0, 3, 0, "scalar_1", false),
                target_edge(1, 0, 1, "scalar_1", false),
                target_edge(2, 1, 2, "scalar_1", false),
                target_edge(3, 0, 2, "scalar_1", false),
                target_edge(4, 2, 4, "scalar_1", false),
            ],
            loop_edges: None,
            overall_factor: OverallFactor::Unit,
        },
    ]
}

fn epemttbar_target() -> Target {
    use ExternalState::{Incoming, Outgoing};
    use TargetVertex::{External, Interaction};

    Target {
        name: "GL092",
        path: "tests/resources/graphs/epemttbar.dot",
        vertices: vec![
            Interaction("V_137"),
            Interaction("V_141"),
            Interaction("V_134"),
            Interaction("V_134"),
            Interaction("V_141"),
            Interaction("V_137"),
            Interaction("V_98"),
            Interaction("V_98"),
            Interaction("V_137"),
            Interaction("V_137"),
            External(Outgoing),
            External(Outgoing),
            External(Incoming),
            External(Incoming),
        ],
        edges: vec![
            target_edge(0, 6, 11, "e-", true),
            target_edge(1, 10, 6, "e-", true),
            target_edge(2, 7, 12, "e+", true),
            target_edge(3, 13, 7, "e+", true),
            target_edge(4, 3, 2, "t", true),
            target_edge(5, 2, 0, "t", true),
            target_edge(6, 0, 8, "t", true),
            target_edge(7, 8, 1, "t", true),
            target_edge(8, 1, 4, "t", true),
            target_edge(9, 4, 9, "t", true),
            target_edge(10, 9, 5, "t", true),
            target_edge(11, 5, 3, "t", true),
            target_edge(12, 5, 0, "g", false),
            target_edge(13, 8, 9, "g", false),
            target_edge(14, 4, 1, "H", false),
            target_edge(15, 3, 6, "a", false),
            target_edge(16, 2, 7, "a", false),
        ],
        loop_edges: Some(vec![12, 14, 8, 13]),
        overall_factor: OverallFactor::NegativeUnit,
    }
}

fn candidate_data(diagram: &FeynmanDiagram, model: &Model) -> Result<Candidate, DynError> {
    let vertices = diagram.vertices().collect::<Vec<_>>();
    let positions = vertices
        .iter()
        .enumerate()
        .map(|(position, (id, _))| (*id, position))
        .collect::<BTreeMap<VertexId, usize>>();
    let candidate_vertices = vertices
        .iter()
        .map(|(_, vertex)| {
            if let Some(interaction) = vertex.interaction {
                Ok(CandidateVertex::Interaction(
                    model.vertex_rule_by_id(interaction)?.name.clone(),
                ))
            } else if let Some(external) = &vertex.external {
                Ok(CandidateVertex::External(external.state))
            } else {
                Err("diagram vertex has neither an interaction nor an external state".into())
            }
        })
        .collect::<Result<Vec<_>, DynError>>()?;
    let edges = diagram
        .edges()
        .map(|(id, endpoints, edge)| {
            Ok(CandidateEdge {
                id,
                source: positions[&endpoints.source],
                target: positions[&endpoints.target],
                particle: model.particle_by_id(edge.particle)?.name.clone(),
                directed: edge.directed,
            })
        })
        .collect::<Result<Vec<_>, DynError>>()?;
    Ok(Candidate {
        vertices: candidate_vertices,
        edges,
    })
}

fn vertex_matches(target: TargetVertex, candidate: &CandidateVertex) -> bool {
    match (target, candidate) {
        (TargetVertex::Interaction(target), CandidateVertex::Interaction(candidate)) => {
            target == candidate
        }
        (TargetVertex::External(target), CandidateVertex::External(candidate)) => {
            target == *candidate
        }
        _ => false,
    }
}

fn edge_matches(target: &TargetEdge, candidate: &CandidateEdge, vertex_mapping: &[usize]) -> bool {
    if target.particle != candidate.particle || target.directed != candidate.directed {
        return false;
    }
    let source = vertex_mapping[target.source];
    let target_vertex = vertex_mapping[target.target];
    (candidate.source == source && candidate.target == target_vertex)
        || (!target.directed && candidate.source == target_vertex && candidate.target == source)
}

fn match_edges(
    target: &Target,
    candidate: &Candidate,
    vertex_mapping: &[usize],
) -> Option<BTreeMap<usize, EdgeId>> {
    let mut used = BTreeSet::new();
    let mut mapping = BTreeMap::new();
    for target_edge in &target.edges {
        let edge = candidate
            .edges
            .iter()
            .filter(|edge| !used.contains(&edge.id))
            .filter(|edge| edge_matches(target_edge, edge, vertex_mapping))
            .min_by_key(|edge| edge.id)?;
        used.insert(edge.id);
        mapping.insert(target_edge.id, edge.id);
    }
    (used.len() == candidate.edges.len()).then_some(mapping)
}

fn exact_match(target: &Target, candidate: &Candidate) -> Option<ExactMatch> {
    if target.vertices.len() != candidate.vertices.len()
        || target.edges.len() != candidate.edges.len()
    {
        return None;
    }

    fn recurse(
        target: &Target,
        candidate: &Candidate,
        index: usize,
        used: &mut [bool],
        vertex_mapping: &mut [usize],
    ) -> Option<ExactMatch> {
        if index == target.vertices.len() {
            return match_edges(target, candidate, vertex_mapping)
                .map(|edge_mapping| ExactMatch { edge_mapping });
        }
        for candidate_index in 0..candidate.vertices.len() {
            if used[candidate_index]
                || !vertex_matches(target.vertices[index], &candidate.vertices[candidate_index])
            {
                continue;
            }
            used[candidate_index] = true;
            vertex_mapping[index] = candidate_index;
            if let Some(found) = recurse(target, candidate, index + 1, used, vertex_mapping) {
                return Some(found);
            }
            used[candidate_index] = false;
        }
        None
    }

    recurse(
        target,
        candidate,
        0,
        &mut vec![false; candidate.vertices.len()],
        &mut vec![usize::MAX; target.vertices.len()],
    )
}

fn select_target(
    target: &Target,
    diagrams: &[FeynmanDiagram],
    model: &Model,
) -> Result<FeynmanDiagram, DynError> {
    let mut matches = Vec::new();
    for diagram in diagrams {
        let candidate = candidate_data(diagram, model)?;
        if let Some(exact) = exact_match(target, &candidate) {
            matches.push((diagram, exact));
        }
    }
    matches.sort_by_key(|(diagram, _)| diagram.id());
    if matches.len() != 1 {
        let diagnostics = diagrams
            .iter()
            .filter_map(|diagram| {
                let candidate = candidate_data(diagram, model).ok()?;
                (candidate.vertices.len() == target.vertices.len()
                    && candidate.edges.len() == target.edges.len())
                .then(|| format!("{} ({}): {candidate:?}", diagram.name(), diagram.id()))
            })
            .collect::<Vec<_>>()
            .join("\n");
        return Err(format!(
            "target '{}' matched {} generated diagrams; refine its deterministic selector. Same-size candidates:\n{}",
            target.name,
            matches.len(),
            diagnostics,
        )
        .into());
    }
    let (diagram, exact) = matches.pop().expect("one match was checked");
    if diagram.cuts().is_empty() {
        return Err(format!("target '{}' has no finalized physical cuts", target.name).into());
    }
    if diagram.topology_threshold_candidates().is_empty() {
        return Err(format!(
            "target '{}' has no finalized topology-threshold candidates",
            target.name
        )
        .into());
    }

    let loop_edges = if let Some(loop_edges) = &target.loop_edges {
        loop_edges
            .iter()
            .map(|edge| {
                exact
                    .edge_mapping
                    .get(edge)
                    .copied()
                    .ok_or_else(|| format!("target '{}' has no edge {edge}", target.name))
            })
            .collect::<Result<Vec<_>, _>>()?
    } else {
        diagram.loop_momentum_basis().loop_edges.clone()
    };
    let mut result = diagram
        .clone()
        .with_loop_momentum_edges(&loop_edges)?
        .with_name(target.name);
    match target.overall_factor {
        OverallFactor::Generated => {}
        OverallFactor::Unit => result = result.with_overall_factor(Atom::one()),
        OverallFactor::NegativeUnit => result = result.with_overall_factor(-Atom::one()),
    }
    result.validate()?;
    Ok(result)
}

fn generate_addbar() -> Result<(Model, Vec<(Target, FeynmanDiagram)>), DynError> {
    let model = standard_model()?;
    let process = Process::cross_section(["a"], ["d", "d~"])
        .with_final_state_alternatives([vec!["d", "d~"], vec!["d", "d~", "g"]])?
        .with_loop_count(2, 2)?;
    let options = GenerationOptions::default()
        .threads(1)
        .max_vertices(4)
        .allow_zero_flow_edges(true)
        .graph_prefix("AUD")
        .with_graph_filter(GenerationFilter::VertexAllow(vec![
            VertexSelector::from("V_71"),
            VertexSelector::from("V_74"),
        ]))
        .with_graph_filter(GenerationFilter::BlobRange(1..=1))
        .with_graph_filter(GenerationFilter::SpectatorRange(0..=0));
    let generated = Generator::new(model.clone()).generate(&process, &options)?;
    let target = addbar_target();
    let diagram = select_target(&target, &generated.diagrams, &model)?;
    Ok((model, vec![(target, diagram)]))
}

fn generate_scalars() -> Result<(Model, Vec<(Target, FeynmanDiagram)>), DynError> {
    let model = scalar_model()?;
    let process =
        Process::cross_section(["scalar_1"], ["scalar_1", "scalar_1"]).with_loop_count(1, 1)?;
    let options = GenerationOptions::default()
        .threads(1)
        .max_vertices(4)
        .allow_self_loops(true)
        .with_graph_filter(GenerationFilter::VertexAllow(vec![
            VertexSelector::from("V_3_SCALAR_111"),
            VertexSelector::from("V_2_SCALAR_11"),
        ]))
        .with_graph_filter(GenerationFilter::BlobRange(1..=1))
        .with_graph_filter(GenerationFilter::SpectatorRange(0..=0));
    let generated = Generator::new(model.clone()).generate(&process, &options)?;
    let diagrams = scalar_targets()
        .into_iter()
        .map(|target| {
            let diagram = select_target(&target, &generated.diagrams, &model)?;
            Ok((target, diagram))
        })
        .collect::<Result<Vec<_>, DynError>>()?;
    Ok((model, diagrams))
}

fn generate_epemttbar(model: &Model) -> Result<Vec<(Target, FeynmanDiagram)>, DynError> {
    // GL092 is the tree-level cross-section topology for the double-real
    // e+ e- -> t t~ H g g final state. Its four graph loops are the sewn
    // phase-space loops, not four virtual loops.
    let process = Process::cross_section(["e+", "e-"], ["t", "t~", "H", "g", "g"])
        .with_loop_count(4, 4)?
        .symmetrize_final(true)
        .symmetrize_left_right(true);
    let options = GenerationOptions::default()
        .threads(1)
        .max_vertices(10)
        .with_graph_filter(GenerationFilter::VertexAllow(vec![
            VertexSelector::from("V_98"),
            VertexSelector::from("V_134"),
            VertexSelector::from("V_137"),
            VertexSelector::from("V_141"),
        ]))
        .with_graph_filter(GenerationFilter::CouplingOrders(BTreeMap::from([
            ("QCD".to_owned(), (4, Some(4))),
            ("QED".to_owned(), (6, Some(6))),
        ])))
        .with_graph_filter(GenerationFilter::BlobRange(1..=1))
        .with_graph_filter(GenerationFilter::SpectatorRange(0..=0));
    let generated = Generator::new(model.clone()).generate(&process, &options)?;
    let target = epemttbar_target();
    let diagram = select_target(&target, &generated.diagrams, model)?;
    Ok(vec![(target, diagram)])
}

fn process_fixture(
    mode: Mode,
    model: &Model,
    target: &Target,
    diagram: &FeynmanDiagram,
) -> Result<(), DynError> {
    let path = PathBuf::from(target.path);
    let dot = diagram.to_dot()?;
    let round_trip = FeynmanDiagram::from_dot(diagram.model_arc(), &dot)?;
    round_trip.validate()?;
    if round_trip.cuts().is_empty() || round_trip.topology_threshold_candidates().is_empty() {
        return Err(format!(
            "canonical round trip lost cut metadata for '{}'",
            target.name
        )
        .into());
    }
    match mode {
        Mode::Write => {
            fs::write(&path, &dot)?;
            println!("wrote {}", path.display());
        }
        Mode::Verify => {
            let checked_in = fs::read_to_string(&path)?;
            if checked_in != dot {
                return Err(format!("{} is stale; regenerate with --write", path.display()).into());
            }
            let parsed = FeynmanDiagram::from_dot(model.clone(), &checked_in)?;
            parsed.validate()?;
            println!("verified {}", path.display());
        }
    }
    Ok(())
}

fn process_group(
    mode: Mode,
    model: &Model,
    fixtures: &[(Target, FeynmanDiagram)],
) -> Result<(), DynError> {
    for (target, diagram) in fixtures {
        process_fixture(mode, model, target, diagram)?;
    }
    Ok(())
}

fn ensure_repository_root() -> Result<(), DynError> {
    if !Path::new("Cargo.toml").is_file() || !Path::new("crates/feynkit").is_dir() {
        return Err("run regenerate_canonical_fixtures from the repository root".into());
    }
    Ok(())
}

fn main() -> Result<(), DynError> {
    ensure_repository_root()?;
    let mode = parse_mode()?;
    let (standard_model, addbar) = generate_addbar()?;
    let epemttbar = generate_epemttbar(&standard_model)?;
    let (scalar_model, scalars) = generate_scalars()?;

    process_group(mode, &standard_model, &addbar)?;
    process_group(mode, &standard_model, &epemttbar)?;
    process_group(mode, &scalar_model, &scalars)?;
    Ok(())
}
