use std::sync::Arc;

use ahash::HashMap;
use ahash::HashMapExt;
use symbolica::graph::Graph;

use crate::model::ColorStructure;
use crate::model::VertexRule;
use crate::tests_from_pytest::load_generic_model;

use super::diagram_generator::NodeColorWithVertexRule;
use super::GenerationType;
use super::{diagram_generator::FeynGen, FeynGenFilter, FeynGenFilters, FeynGenOptions};

#[test]
fn cut_content() {
    let model = load_generic_model("sm");

    let mut coupling = HashMap::new();
    coupling.insert("QED".into(), 6);
    let mut pert = HashMap::new();
    pert.insert("QCD".into(), 1);
    let filters = FeynGen::new(FeynGenOptions {
        generation_type: GenerationType::CrossSection,
        initial_pdgs: vec![-11, 11],
        final_pdgs: vec![5, -5, 25],
        loop_count_range: (3, 3),
        symmetrize_initial_states: true,
        symmetrize_final_states: true,
        symmetrize_left_right_states: true,
        amplitude_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(Default::default()),
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 12, 14, 16, 2, 4, 6, 3, 250, 251, 13,
                15,
            ]),
            FeynGenFilter::TadpolesFilter(Default::default()),
            FeynGenFilter::ZeroSnailsFilter(Default::default()),
        ]),
        cross_section_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(Default::default()),
            FeynGenFilter::ParticleVeto(vec![
                23, 24, 9000001, 9000002, 9000003, 9000004, 12, 14, 16, 2, 4, 6, 3, 250, 251, 13,
                15,
            ]),
            FeynGenFilter::TadpolesFilter(Default::default()),
            FeynGenFilter::ZeroSnailsFilter(Default::default()),
            FeynGenFilter::PerturbativeOrders(pert),
            FeynGenFilter::CouplingOrders(coupling),
            FeynGenFilter::LoopCountRange((3, 3)),
        ]),
    });

    let mut graph = Graph::new();
    #[allow(non_snake_case)]
    let bbH = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule(&"V_78".into()),
    };
    let bba = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule(&"V_73".into()),
    };
    let bbg = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule(&"V_76".into()),
    };
    let epema = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule(&"V_98".into()),
    };

    let dummy_external_vertex_rule = Arc::new(VertexRule {
        name: "external".into(),
        couplings: vec![],
        lorentz_structures: vec![],
        particles: vec![],
        color_structures: ColorStructure::new(vec![]),
    });

    let v0 = graph.add_node(bbH.clone());
    let v1 = graph.add_node(bbH.clone());
    let v2 = graph.add_node(bbg.clone());
    let v3 = graph.add_node(bba.clone());
    let v4 = graph.add_node(bbg.clone());
    let v5 = graph.add_node(bba.clone());
    let v6 = graph.add_node(epema.clone());
    let v7 = graph.add_node(epema.clone());
    let v8 = graph.add_node(NodeColorWithVertexRule {
        external_tag: 1,
        vertex_rule: dummy_external_vertex_rule.clone(),
    });
    let v9 = graph.add_node(NodeColorWithVertexRule {
        external_tag: 2,
        vertex_rule: dummy_external_vertex_rule.clone(),
    });
    let v10 = graph.add_node(NodeColorWithVertexRule {
        external_tag: 3,
        vertex_rule: dummy_external_vertex_rule.clone(),
    });
    let v11 = graph.add_node(NodeColorWithVertexRule {
        external_tag: 4,
        vertex_rule: dummy_external_vertex_rule.clone(),
    });

    graph.add_edge(v0, v1, true, "H").unwrap();
    graph.add_edge(v0, v5, true, "b").unwrap();
    graph.add_edge(v1, v3, true, "b").unwrap();
    graph.add_edge(v2, v4, true, "g").unwrap();
    graph.add_edge(v2, v4, true, "b").unwrap();
    graph.add_edge(v3, v2, true, "b").unwrap();
    graph.add_edge(v3, v7, true, "a").unwrap();
    graph.add_edge(v4, v1, true, "b").unwrap();
    graph.add_edge(v5, v0, true, "b").unwrap();
    graph.add_edge(v5, v6, true, "a").unwrap();
    graph.add_edge(v6, v10, false, "e+").unwrap();
    graph.add_edge(v6, v11, false, "e-").unwrap();
    graph.add_edge(v8, v7, false, "e+").unwrap();
    graph.add_edge(v9, v7, false, "e-").unwrap();

    let (n_unresolved, unresolved_type) = filters.unresolved_cut_content(&model);
    assert!(!filters.contains_cut(&model, &graph, n_unresolved, &unresolved_type));
}
