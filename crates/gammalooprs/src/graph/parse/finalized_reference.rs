//! Regression snapshot for the finalized FeynKit/GammaLoop graph boundary.
//! It covers scalar and vector loops, fermions and antifermions, amplitude
//! dangling edges and sewn cross-section initial states. The reference includes
//! current main's Standard Model corrections and the native initial-state and
//! ghost-loop conventions. Every coordinate and symbolic fragment is compared
//! exactly, without reimporting or rewriting the stored reference.

use crate::initialisation;
use crate::model::{InputParamCard, ModelGammaLoopExt};
use crate::utils;
use crate::{graph::Graph, uv::UltravioletGraph};
use feynkit_generator::{GenerationOptions, GenerationType, Process};
use linnet::half_edge::{
    involution::Hedge,
    subgraph::{SuBitGraph, SubSetLike},
};
use serde_json::{Value, json};
use symbolica::atom::AtomCore;

fn selection(selected: &SuBitGraph) -> Vec<usize> {
    selected.included_iter().map(|hedge| hedge.0).collect()
}

fn snapshot(graph: &Graph, model: &crate::model::Model) -> Value {
    let edges = graph.underlying.iter_edges().map(|(pair, edge, data)| {
        let signature = &graph.loop_momentum_basis.edge_signatures[edge];
        json!({
            "id":edge.0,"pair":format!("{pair:?}"),"orientation":format!("{:?}",data.orientation),
            "particle":data.data.particle.particle().map(|particle|particle.index()),
            "num":crate::graph::parse::NumeratorInterface::export(&data.data.num.value).0.to_canonical_string(),"dummy":data.data.is_dummy,
            "mass":data.data.particle.mass_atom(model).to_canonical_string(),
            "internal":signature.internal.iter().map(|x|format!("{x:?}")).collect::<Vec<_>>(),
            "external":signature.external.iter().map(|x|format!("{x:?}")).collect::<Vec<_>>()
        })
    }).collect::<Vec<_>>();
    let vertices=graph.underlying.iter_nodes().map(|(node,_,data)| json!({
        "id":node.0,"rule":data.vertex_rule.map(|rule|rule.index()),"num":crate::graph::parse::NumeratorInterface::export(&data.num.value).0.to_canonical_string()
    })).collect::<Vec<_>>();
    let hedges=(0..graph.underlying.n_hedges()).map(|index| {
        let hedge=Hedge(index);
        json!({"id":index,"node":graph.underlying.node_id(hedge).0,"edge":graph.underlying[&hedge].0,
            "flow":format!("{:?}",graph.underlying.flow(hedge)),"ufo_order":graph.underlying[hedge].ufo_order.value})
    }).collect::<Vec<_>>();
    json!({"edges":edges,"vertices":vertices,"hedges":hedges,
        "tree":selection(&graph.loop_momentum_basis.tree),
        "loop_edges":graph.loop_momentum_basis.loop_edges.iter().map(|edge|edge.0).collect::<Vec<_>>(),
        "external_edges":graph.loop_momentum_basis.ext_edges.iter().map(|edge|edge.0).collect::<Vec<_>>(),
        "initial_left":selection(&graph.initial_state_cut.left),"initial_right":selection(&graph.initial_state_cut.right),
        "num":graph.numerator(&graph.full_filter(),&graph.empty_subgraph()).unwrap().get_single_atom().unwrap().to_canonical_string(),
        "overall":graph.overall_factor.to_canonical_string(),"prefactor":crate::graph::parse::NumeratorInterface::export(&graph.global_prefactor.num).0.to_canonical_string(),
        "projector":crate::graph::parse::NumeratorInterface::export(&graph.global_prefactor.projector).0.to_canonical_string(),
        "cuts":graph.finalized_cuts.iter().map(|cut|json!({"cut":selection(&cut.cut.left),"opposite":selection(&cut.cut.right),"left":selection(&cut.left),"right":selection(&cut.right)})).collect::<Vec<_>>(),
        "thresholds":graph.finalized_topology_threshold_candidates.iter().map(|cut|json!({"cut":selection(&cut.cut.left),"opposite":selection(&cut.cut.right),"left":selection(&cut.left),"right":selection(&cut.right)})).collect::<Vec<_>>()
    })
}

#[test]
fn finalized_graph_matches_gammaloop_reference() {
    crate::initialisation::test_initialise().unwrap();
    let cases = [
        (
            "scalar_amplitude",
            "scalars",
            GenerationType::Amplitude,
            Process::new(["scalar_1"], ["scalar_1"]),
            1,
        ),
        (
            "scalar_cross_section",
            "scalars",
            GenerationType::CrossSection,
            Process::new(["scalar_1"], ["scalar_1", "scalar_1"]),
            1,
        ),
        (
            "fermion_amplitude",
            "sm",
            GenerationType::Amplitude,
            Process::new(["e-", "e+"], ["mu-", "mu+"]),
            0,
        ),
        (
            "vector_amplitude",
            "sm",
            GenerationType::Amplitude,
            Process::new(["g"], ["g"]),
            1,
        ),
        (
            "fermion_cross_section",
            "sm",
            GenerationType::CrossSection,
            Process::new(["e-", "e+"], ["a"]),
            0,
        ),
    ];
    let mut result = serde_json::Map::new();
    for (name, model, generation_type, process, loops) in cases {
        let model = std::sync::Arc::new(crate::utils::load_generic_model(model));
        let options = GenerationOptions::default()
            .with_loop_count(loops, loops)
            .unwrap()
            .threads(1)
            .max_vertices(2);
        let generated = match generation_type {
            GenerationType::Amplitude => process.generate_diagrams(model.clone(), &options),
            GenerationType::CrossSection => process.generate_cross_section(model.clone(), &options),
        }
        .unwrap();
        assert!(!generated.diagrams.is_empty(), "{name}");
        let snapshots = generated
            .diagrams
            .iter()
            .map(|diagram| snapshot(&Graph::from_feynkit(diagram, None, true).unwrap(), &model))
            .collect::<Vec<_>>();
        result.insert(name.into(), json!(snapshots));
    }
    let expected: serde_json::Map<String, Value> =
        serde_json::from_str(include_str!("finalized_reference.json")).unwrap();
    assert_eq!(
        result.keys().collect::<Vec<_>>(),
        expected.keys().collect::<Vec<_>>()
    );
    for (name, diagrams) in &result {
        let diagrams = diagrams.as_array().unwrap();
        let reference = expected[name].as_array().unwrap();
        assert_eq!(diagrams.len(), reference.len(), "diagram inventory: {name}");
        for (index, (actual, expected)) in diagrams.iter().zip(reference).enumerate() {
            for (field, expected) in expected.as_object().unwrap() {
                assert_eq!(&actual[field], expected, "{name} diagram {index}: {field}");
            }
        }
    }
}

#[test]
fn frozen_sm_graphs_match_the_run_card_model() {
    initialisation::test_initialise().unwrap();
    let mut model = utils::load_generic_model("sm")
        .with_restriction(Some("default".to_owned()))
        .unwrap();
    let mut card = InputParamCard::from_str(
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../assets/models/json/sm/restrict_default.json"
        ))
        .to_owned(),
        "json",
    )
    .unwrap();
    model.simplify(&mut card).unwrap();
    for (fixture, source) in [
        (
            "addbar",
            include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/resources/graphs/addbar.dot"
            )),
        ),
        (
            "epemttbar",
            include_str!(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/../../tests/resources/graphs/epemttbar.dot"
            )),
        ),
    ] {
        let bindings = source
            .lines()
            .filter_map(|line| {
                let (_, tail) = line.split_once("interaction=\"")?;
                let name = tail.split('"').next()?;
                Some((name, model.vertex_rule_id(name).unwrap().index()))
            })
            .collect::<Vec<_>>();
        let imported = Graph::import_dot(source, &model).unwrap_or_else(|error| {
            panic!(
                "{fixture}: model {}, bindings {bindings:?}: {error}",
                model.fingerprint()
            )
        });
        assert!(!imported.is_empty());
    }
}
