//! Prepare the massless QCD numerator from native model rules, or reuse its verified scalar input.
//! The transverse quadratic coefficient is not a complete beta function.

// The fixture is also included by the original 5L collector and library unit
// tests, which do not invoke the generic comparison driver.
#[allow(dead_code)]
pub mod comparison;
mod physical;
mod prepared;
mod projector;
mod vacuum;

use std::{
    collections::{BTreeMap, BTreeSet},
    error::Error,
    path::Path,
    sync::Arc,
    time::Instant,
};

use super::vakint::{
    Vakint, VakintExpression, diagram_integral::DiagramIntegralOptions, symbols::S,
};
use feynkit_generator::{
    EdgeColor, GenerationFilter, GenerationOptions, NodeColor, Process, VertexSelector,
};
use feynkit_graph::IntegralFamily;
use feynkit_graph::{
    DiagramEdge, DiagramVertex, EdgeId, ExternalState, FeynmanDiagram, NumeratorInterface,
    VertexId, symbols,
};
use feynkit_kinematics::Kinematics;
use feynkit_model::Model;
use idenso::{
    IndexTooling,
    color::ColorSimplifySettings,
    tensor::{AlgebraContraction, AlgebraSettings, ReductionStatus, SymbolicTensor},
};
use rustred::algebra::CoefficientContext;
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::slot::IsAbstractSlot,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    domains::integer::Integer,
    graph::Graph,
    symbol,
};

pub use prepared::PreparedInputManifest;

pub const PREPARED_INPUT_MANIFEST: &str = "prepared-input.json";

pub type Result<T> = std::result::Result<T, Box<dyn Error>>;

#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Projection {
    QuadraticTransverse,
    ZeroMomentumMetric,
}

impl Projection {
    pub fn degree(self, loop_count: usize) -> i32 {
        match self {
            Self::QuadraticTransverse => 2 * loop_count as i32,
            Self::ZeroMomentumMetric => 2 * loop_count as i32 - 1,
        }
    }
}

#[derive(Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct PreparationTimings {
    pub graph_generation_seconds: f64,
    pub scalar_contraction_seconds: f64,
    pub taylor_seconds: f64,
    pub angular_seconds: f64,
    pub ingress_seconds: f64,
    pub total_seconds: f64,
}

pub struct LadderFixture {
    pub loop_count: usize,
    pub input: Atom,
    pub preparation_seconds: f64,
    pub scalar_terms: usize,
    pub contributions: usize,
    pub nonzero_conditions: Vec<Atom>,
    pub projection: Projection,
    pub timings: PreparationTimings,
    /// Historical preparation provenance, present only when this input was loaded.
    pub prepared_input: Option<PreparedInputManifest>,
}

type Shape = Graph<u8, i64>;

fn evaluate_overall_factor(factor: AtomView<'_>) -> Atom {
    // The graph retains generation provenance as tagged scalar factors.
    let mut result = factor.to_owned();
    let value = Atom::var(symbol!("qcd_ladder::factor_"));
    for head in [
        "AutG",
        "CouplingsMultiplicity",
        "InternalFermionLoopSign",
        "ExternalFermionOrderingSign",
        "AntiFermionSpinSumSign",
        "NumeratorIndependentSymmetryGrouping",
    ] {
        result = result
            .replace(symbol!(&format!("feynkit_generator_factor::{head}")).call(&value))
            .with(value.clone());
    }
    result.expand()
}

fn is_external_edge(_diagram: &FeynmanDiagram, ends: feynkit_graph::EdgeEndpoints) -> bool {
    ends.source.is_none() || ends.target.is_none()
}

fn external_color(state: ExternalState) -> u8 {
    match state {
        ExternalState::Incoming => 1,
        ExternalState::Outgoing => 2,
    }
}

fn diagram_labeling(diagram: &FeynmanDiagram) -> Result<(Shape, BTreeMap<VertexId, usize>)> {
    let mut graph = Shape::new();
    let vertices = diagram
        .vertices()
        .map(|(id, _)| (id, graph.add_node(0)))
        .collect::<BTreeMap<_, _>>();
    for (_, ends, edge) in diagram.edges() {
        let mut endpoint = |vertex: Option<VertexId>| {
            vertex.map_or_else(
                || {
                    graph.add_node(external_color(
                        edge.external
                            .as_ref()
                            .expect("external edge metadata")
                            .state,
                    ))
                },
                |id| vertices[&id],
            )
        };
        let source = endpoint(ends.source);
        let target = endpoint(ends.target);
        let pdg = diagram.model().particle_by_id(edge.particle)?.pdg_code;
        let (source, target) = if pdg < 0 {
            (target, source)
        } else {
            (source, target)
        };
        graph.add_edge(source, target, edge.directed, pdg.abs())?;
    }
    let canonical = graph.canonize();
    let labels = vertices
        .into_iter()
        .map(|(id, position)| (id, canonical.vertex_map[position]))
        .collect();
    Ok((canonical.graph, labels))
}

fn diagram_shape(diagram: &FeynmanDiagram) -> Result<Shape> {
    Ok(diagram_labeling(diagram)?.0)
}

fn compact_rail_basis(diagram: FeynmanDiagram, target: &FeynmanDiagram) -> Result<FeynmanDiagram> {
    let (target_shape, target_labels) = diagram_labeling(target)?;
    let (actual_shape, actual_labels) = diagram_labeling(&diagram)?;
    if target_shape != actual_shape {
        return Err("rail routing lacks the native exact colored-graph witness".into());
    }
    let named = target
        .vertices()
        .map(|(id, vertex)| (vertex.name.as_str(), target_labels[&id]))
        .collect::<BTreeMap<_, _>>();
    let rail = std::iter::once("left".to_owned())
        .chain((1..target.loop_count()).map(|rung| format!("t{rung}")))
        .chain(["right".to_owned()])
        .collect::<Vec<_>>();
    let chords = rail
        .windows(2)
        .map(|pair| {
            diagram
                .edges()
                .find_map(|(id, ends, edge)| {
                    let pdg = diagram
                        .model()
                        .particle_by_id(edge.particle)
                        .expect("native particle")
                        .pdg_code;
                    if pdg.abs() != 1 || is_external_edge(&diagram, ends) {
                        return None;
                    }
                    let (left, right) = (
                        actual_labels[&ends.source.expect("internal edge")],
                        actual_labels[&ends.target.expect("internal edge")],
                    );
                    let (left, right) = if pdg < 0 {
                        (right, left)
                    } else {
                        (left, right)
                    };
                    (left == named[pair[0].as_str()] && right == named[pair[1].as_str()])
                        .then_some(id)
                })
                .ok_or_else(|| "native target witness did not identify a fermion rail edge".into())
        })
        .collect::<Result<Vec<_>>>()?;
    // Native spanning-tree validation proves these compact p_i are loop
    // coordinates. Dirac algebra then precedes expansion into k-l1-... form.
    Ok(diagram.with_loop_momentum_edges(&chords)?)
}

fn generated_shape(graph: &Graph<NodeColor, EdgeColor>, model: &Model) -> Shape {
    let mut result = Shape::new();
    for node in graph.nodes() {
        result.add_node(
            node.data
                .external
                .as_ref()
                .map_or(0, |external| external_color(external.state)),
        );
    }
    for edge in graph.edges() {
        let pdg = model
            .particle_by_id(edge.data.particle)
            .expect("native generated particle")
            .pdg_code;
        let (source, target) = if pdg < 0 {
            (edge.vertices.1, edge.vertices.0)
        } else {
            edge.vertices
        };
        result
            .add_edge(source, target, edge.directed, pdg.abs())
            .expect("native generated endpoints");
    }
    result.canonize().graph
}

fn ladder_target(loop_count: usize) -> Result<String> {
    if loop_count == 0 {
        return Err("a quark-loop ladder requires at least one loop".into());
    }
    let rungs = loop_count - 1;
    let mut dot = String::from(
        "digraph massless_quark_ladder_self_energy {\n incoming [style=invis];\n outgoing [style=invis];\n incoming -> left [particle=\"g\"];\n right -> outgoing [particle=\"g\"];\n",
    );
    let rail = std::iter::once("left".to_owned())
        .chain((1..=rungs).map(|rung| format!("t{rung}")))
        .chain(["right".to_owned()])
        .chain((1..=rungs).rev().map(|rung| format!("b{rung}")))
        .chain(["left".to_owned()])
        .collect::<Vec<_>>();
    for pair in rail.windows(2) {
        dot.push_str(&format!(" {} -> {} [particle=\"d\"];\n", pair[0], pair[1]));
    }
    for rung in 1..=rungs {
        dot.push_str(&format!(" t{rung} -> b{rung} [particle=\"g\"];\n"));
    }
    dot.push_str("}\n");
    Ok(dot)
}

#[cfg(test)]
fn dressed_diagram() -> Result<feynkit_generator::GenerationResult> {
    dressed_diagram_at_loop_count(5)
}

fn dressed_diagram_at_loop_count(loop_count: usize) -> Result<feynkit_generator::GenerationResult> {
    let model = Arc::new(
        Model::from_json(include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/fixtures/qcd.json"
        )))
        .expect("valid qcd fixture"),
    );
    let target = FeynmanDiagram::from_dot(Arc::clone(&model), &ladder_target(loop_count)?)?;
    if target.loop_count() != loop_count
        || target.edges().count() != 3 * loop_count + 1
        || target.vertices().count() != 2 * loop_count
    {
        return Err("the source QCD ladder has an inconsistent native incidence census".into());
    }
    if loop_count == 5
        && diagram_shape(&target)?
            != diagram_shape(&FeynmanDiagram::from_dot(
                Arc::clone(&model),
                include_str!("ladder.dot"),
            )?)?
    {
        return Err("parameterization changed the original four-rung source".into());
    }
    let shape = diagram_shape(&target)?;
    let mut ddg = vec![
        model.particle_id("d")?,
        model.particle_id("d~")?,
        model.particle_id("g")?,
    ];
    ddg.sort();
    let rules = model
        .vertex_rules()
        .iter()
        .filter(|rule| {
            let mut particles = rule.particles.clone();
            particles.sort();
            particles == ddg
        })
        .map(|rule| VertexSelector::from(rule.name.clone()))
        .collect::<Vec<_>>();
    if rules.len() != 1 {
        return Err("expected exactly one native down-quark QCD interaction".into());
    }
    let filter_model = Arc::clone(&model);
    let options = GenerationOptions::default()
        .with_loop_count(loop_count, loop_count)?
        .max_vertices(2 * loop_count)
        .threads(1)
        .with_graph_filter(GenerationFilter::VertexAllow(rules))
        .with_graph_filter(GenerationFilter::FermionLoopCountRange((1, 1)))
        .filter(move |graph, completed| {
            if graph.nodes().len() > 2 * loop_count + 2 || graph.edges().len() > 3 * loop_count + 1
            {
                return false;
            }
            if completed != graph.nodes().len() {
                return true;
            }
            graph.nodes().len() == 2 * loop_count + 2
                && graph.edges().len() == 3 * loop_count + 1
                && generated_shape(graph, &filter_model) == shape
        });
    let mut generated = Process::new(["g"], ["g"]).generate_diagrams(model, &options)?;
    if !generated.report.completed || generated.diagrams.len() != 1 {
        return Err(format!(
            "native exact ladder selection returned {} diagrams, completed={}",
            generated.diagrams.len(),
            generated.report.completed
        )
        .into());
    }
    if diagram_shape(&generated.diagrams[0])? != diagram_shape(&target)? {
        return Err("native dressing changed the selected colored ladder topology".into());
    }
    let diagram = generated.diagrams.remove(0);
    generated
        .diagrams
        .push(compact_rail_basis(diagram, &target)?);
    Ok(generated)
}

/// Native graph/rule preparation without committing to an external projector.
pub fn dress(evidence: &Path) -> Result<feynkit_generator::GenerationResult> {
    dress_at_loop_count(evidence, 5)
}

pub fn dress_at_loop_count(
    evidence: &Path,
    loop_count: usize,
) -> Result<feynkit_generator::GenerationResult> {
    std::fs::create_dir_all(evidence)?;
    Vakint::initialize_vakint_symbols();
    let generated = dressed_diagram_at_loop_count(loop_count)?;
    let diagram = &generated.diagrams[0];
    if diagram.numerator().expression().is_one()
        || diagram
            .vertices()
            .any(|(_, vertex)| vertex.numerator.expression().is_one())
    {
        return Err("bare structural DOT was mistaken for a dressed QCD numerator".into());
    }
    std::fs::write(evidence.join("dressed-source.dot"), diagram.to_dot()?)?;
    // Cheap native denominator/incidence validation is also exercised by
    // --dress-only; no Dirac algebra or Taylor projection runs in this mode.
    vacuum::Vacuum::new(diagram, &dimension()?)?.validate_initial_ingress(evidence)?;
    Ok(generated)
}

fn dimension() -> Result<Atom> {
    Ok(CoefficientContext::try_new(["d"])?
        .parameter("d")
        .ok_or("missing native dimension")?
        .to_expression())
}

/// Upgrade every Lorentz representation occurrence, including index-free tags.
/// Color stays SU(3), and the native bispinor dimension/trace of unity stays 4.
fn lorentz_dimension(expression: &Atom, dimension: &Atom) -> Atom {
    expression.replace_map(|view, _, output| {
        if let AtomView::Fun(fun) = view
            && fun.get_symbol() == symbol!("spenso::mink")
        {
            let mut args = fun.iter();
            if args.next().is_some() {
                **output = FunctionBuilder::new(fun.get_symbol())
                    .add_arg(dimension)
                    .add_args(args)
                    .finish();
            }
        }
    })
}

fn assert_lorentz_dimension(expression: &Atom, dimension: &Atom) -> Result<()> {
    let mut wrong = false;
    let _ = expression.replace_map(|view, _, _| {
        if let AtomView::Fun(fun) = view
            && fun.get_symbol() == symbol!("spenso::mink")
        {
            wrong |= fun
                .iter()
                .next()
                .is_none_or(|value| value != dimension.as_view());
        }
    });
    if wrong {
        return Err("a source Lorentz representation retained a non-d dimension".into());
    }
    Ok(())
}

fn loop_rows(diagram: &FeynmanDiagram) -> BTreeMap<EdgeId, Atom> {
    let loops = (0..diagram.loop_count())
        .map(|axis| symbols::loop_momentum().call(axis))
        .collect::<Vec<_>>();
    diagram
        .edges()
        .filter(|(_, ends, _)| !is_external_edge(diagram, *ends))
        .map(|(id, _, _)| {
            let row = diagram.loop_momentum_basis().edge_signatures[&id]
                .loops
                .integer_coefficients();
            let momentum = row
                .iter()
                .zip(&loops)
                .map(|(value, k)| Atom::num(*value as i64) * k)
                .sum();
            (id, momentum)
        })
        .collect()
}

pub fn prepare(evidence: &Path, projection: Projection) -> Result<LadderFixture> {
    prepare_at_loop_count(evidence, projection, 5)
}

pub fn prepare_at_loop_count(
    evidence: &Path,
    projection: Projection,
    loop_count: usize,
) -> Result<LadderFixture> {
    let started = Instant::now();
    std::fs::create_dir_all(evidence)?;
    Vakint::initialize_vakint_symbols();
    let generated = dress_at_loop_count(evidence, loop_count)?;
    let mut timings = PreparationTimings {
        graph_generation_seconds: started.elapsed().as_secs_f64(),
        ..Default::default()
    };
    let projected =
        physical::prepare(generated, &dimension()?, evidence, projection, &mut timings)?;
    let ingress_started = Instant::now();
    let mut input = Atom::Zero;
    let mut scalar_terms = 0usize;
    let mut inventory =
        String::from("# unlowered contribution A r_mom D_source=A-r_mom explicit_powers\n");
    for (index, contribution) in projected.contributions.iter().enumerate() {
        // Count the stored representation without an optional expansion pass.
        scalar_terms += if let AtomView::Add(sum) = contribution.numerator.as_view() {
            sum.get_nargs()
        } else {
            1
        };
        let total_power = contribution.powers.iter().sum::<i32>();
        inventory.push_str(&format!(
            "{index} {total_power} {} {} {:?}\n",
            total_power - projection.degree(loop_count),
            projection.degree(loop_count),
            contribution.powers
        ));
        let routed_numerator = projected.vacuum.mapping.apply(&contribution.numerator);
        let diagram = projected
            .vacuum
            .diagram
            .map_data(
                |_, vertex| vertex.clone(),
                |id, _, edge| {
                    let mut edge = edge.clone();
                    edge.propagator_power = contribution.powers[id.0];
                    edge
                },
            )?
            .with_name(format!("qcd_ladder_{projection:?}_{index:03}"))
            .with_numerator(NumeratorInterface::admit(routed_numerator, None)?)?;
        std::fs::write(
            evidence.join(format!("projected-vacuum-{index:03}.dot")),
            diagram.to_dot()?,
        )?;
        input += projected
            .vacuum
            .input(&contribution.numerator, &contribution.powers)?;
    }
    std::fs::write(evidence.join("physical-power-counting.txt"), inventory)?;
    std::fs::write(
        evidence.join("vakint-input.symbolica"),
        input.to_canonical_string(),
    )?;
    std::fs::write(
        evidence.join("projector-conditions.symbolica"),
        projected
            .conditions
            .iter()
            .map(AtomCore::to_canonical_string)
            .collect::<Vec<_>>()
            .join("\n"),
    )?;
    let projection_convention = match projection {
        Projection::QuadraticTransverse => {
            "Complete homogeneous quadratic Taylor coefficient of numerator AND denominators.\nTransverse contraction delta_ab*(g_munu-q_mu*q_nu/q^2)/[8*(d-1)*q^2]; normalized auxiliary angular average at fixed q^2.\n"
        }
        Projection::ZeroMomentumMetric => {
            "Both external current vertices retained; unaveraged g_munu*delta_ab contraction at q=0.\nNo transverse projector, no division by 8(d-1), and no auxiliary angular integration.\n"
        }
    };
    std::fs::write(
        evidence.join("conventions.txt"),
        format!(
            "Native model-rule-generated massless QCD quark loop with {} gluon rungs; fixed SU(3), tr_Dirac(1)=4, generic Lorentz dimension d.\n{projection_convention}L={loop_count}, original A={}, {} vacuum incidence edges, {} distinct quadratic forms.\nUnlowered D_source=A-r_mom={}; omega_UV(d=4)={}.\nCommon UV mass inserted only into denominators; original slash numerators remain massless.\nNot a complete beta function. Verified saved scalar inputs separate preparation from IBP timing.\n",
            loop_count - 1,
            3 * loop_count - 1,
            projected.vacuum.lines.len(),
            2 * loop_count - 1,
            projection.degree(loop_count),
            4 * loop_count as i32 - 2 * projection.degree(loop_count),
        ),
    )?;
    timings.ingress_seconds = ingress_started.elapsed().as_secs_f64();
    timings.total_seconds = started.elapsed().as_secs_f64();
    std::fs::write(
        evidence.join("preparation.json"),
        serde_json::to_vec_pretty(&serde_json::json!({
            "projection": projection, "loop_count": loop_count, "gluon_rungs": loop_count - 1,
            "timings": timings, "scalar_terms": scalar_terms,
            "contributions": projected.contributions.len(), "unlowered_A_minus_r_mom": projection.degree(loop_count),
            "fresh_native_preparation": true,
        }))?,
    )?;
    Ok(LadderFixture {
        loop_count,
        input,
        preparation_seconds: started.elapsed().as_secs_f64(),
        scalar_terms,
        contributions: projected.contributions.len(),
        nonzero_conditions: projected.conditions.into_iter().collect(),
        projection,
        timings,
        prepared_input: None,
    })
}

#[cfg(test)]
#[test]
fn parameterized_ladder_has_one_fermion_cycle_and_expected_incidence() -> Result<()> {
    let model = Arc::new(Model::from_json(include_str!("../qcd.json"))?);
    for loop_count in 1..=5 {
        let target = FeynmanDiagram::from_dot(Arc::clone(&model), &ladder_target(loop_count)?)?;
        assert_eq!(target.loop_count(), loop_count);
        assert_eq!(target.vertices().count(), 2 * loop_count);
        assert_eq!(target.edges().count(), 3 * loop_count + 1);
        assert_eq!(
            target
                .edges()
                .filter(|(_, _, edge)| edge.external.is_none())
                .count(),
            3 * loop_count - 1
        );
    }
    let original = FeynmanDiagram::from_dot(Arc::clone(&model), include_str!("ladder.dot"))?;
    let parameterized = FeynmanDiagram::from_dot(model, &ladder_target(5)?)?;
    assert_eq!(diagram_shape(&original)?, diagram_shape(&parameterized)?);
    assert!(ladder_target(0).is_err());
    Ok(())
}
