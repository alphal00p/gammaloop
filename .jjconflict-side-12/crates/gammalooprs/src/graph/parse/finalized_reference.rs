//! Regression reference captured from unchanged GammaLoop before shared graph
//! finalization. It covers scalar and vector loops, fermions and antifermions,
//! amplitude dangling edges and sewn cross-section initial states. Every native
//! coordinate and symbolic fragment is compared without importing or rewriting
//! the reference through the new shared implementation.

use crate::{graph::Graph, uv::UltravioletGraph};
use feynkit_generator::{GenerationOptions, Generator, Process};
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
            "num":data.data.num.value.to_canonical_string(),"dummy":data.data.is_dummy,
            "mass":data.data.particle.mass_atom(model).to_canonical_string(),
            "internal":signature.internal.iter().map(|x|format!("{x:?}")).collect::<Vec<_>>(),
            "external":signature.external.iter().map(|x|format!("{x:?}")).collect::<Vec<_>>()
        })
    }).collect::<Vec<_>>();
    let vertices=graph.underlying.iter_nodes().map(|(node,_,data)| json!({
        "id":node.0,"rule":data.vertex_rule.map(|rule|rule.index()),"num":data.num.value.to_canonical_string()
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
        "num":graph.numerator(&graph.full_filter(),&graph.empty_subgraph()).get_single_atom().unwrap().to_canonical_string(),
        "overall":graph.overall_factor.to_canonical_string(),"prefactor":graph.global_prefactor.num.to_canonical_string(),
        "projector":graph.global_prefactor.projector.to_canonical_string(),
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
            Process::amplitude(["scalar_1"], ["scalar_1"])
                .with_loop_count(1, 1)
                .unwrap(),
        ),
        (
            "scalar_cross_section",
            "scalars",
            Process::cross_section(["scalar_1"], ["scalar_1", "scalar_1"])
                .with_loop_count(1, 1)
                .unwrap(),
        ),
        (
            "fermion_amplitude",
            "sm",
            Process::amplitude(["e-", "e+"], ["mu-", "mu+"])
                .with_loop_count(0, 0)
                .unwrap(),
        ),
        (
            "vector_amplitude",
            "sm",
            Process::amplitude(["g"], ["g"])
                .with_loop_count(1, 1)
                .unwrap(),
        ),
        (
            "fermion_cross_section",
            "sm",
            Process::cross_section(["e-", "e+"], ["a"])
                .with_loop_count(0, 0)
                .unwrap(),
        ),
    ];
    let mut result = serde_json::Map::new();
    for (name, model, process) in cases {
        let model = std::sync::Arc::new(crate::utils::load_generic_model(model));
        let generator = Generator::new(std::sync::Arc::clone(&model));
        let generated = generator
            .generate(
                &process,
                &GenerationOptions::default().threads(1).max_vertices(2),
            )
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
