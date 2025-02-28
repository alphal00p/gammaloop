use std::sync::Arc;

use ahash::HashMap;
use ahash::HashMapExt;
use symbolica::graph::Graph;

use crate::feyngen::diagram_generator::EdgeColor;
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
        allow_symmetrization_of_external_fermions_in_amplitudes: false,
        max_multiplicity_for_fast_cut_filter: 6,
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
    let e1 = NodeColorWithVertexRule {
        external_tag: 1,
        vertex_rule: dummy_external_vertex_rule.clone(),
    };
    let mut e2 = e1.clone();
    e2.external_tag = 2;
    let mut e3 = e1.clone();
    e3.external_tag = 3;
    let mut e4 = e1.clone();
    e4.external_tag = 4;

    let v8 = graph.add_node(e1.clone());
    let v9 = graph.add_node(e2.clone());
    let v10 = graph.add_node(e3.clone());
    let v11 = graph.add_node(e4.clone());

    let h = EdgeColor::from_particle(model.get_particle(&"H".to_string().into()));

    let b = EdgeColor::from_particle(model.get_particle(&"b".to_string().into()));
    let bbar = EdgeColor::from_particle(
        model
            .get_particle(&"b".to_string().into())
            .get_anti_particle(&model),
    );
    let g = EdgeColor::from_particle(model.get_particle(&"g".to_string().into()));
    let a = EdgeColor::from_particle(model.get_particle(&"a".to_string().into()));
    let eplus = EdgeColor::from_particle(model.get_particle(&"e+".to_string().into()));
    let eminus = EdgeColor::from_particle(model.get_particle(&"e-".to_string().into()));
    graph.add_edge(v0, v1, true, h).unwrap();
    graph.add_edge(v0, v5, true, b).unwrap();
    graph.add_edge(v1, v3, true, b).unwrap();
    graph.add_edge(v2, v4, true, g).unwrap();
    graph.add_edge(v2, v4, true, b).unwrap();
    graph.add_edge(v3, v2, true, b).unwrap();
    graph.add_edge(v3, v7, true, a).unwrap();
    graph.add_edge(v4, v1, true, b).unwrap();
    graph.add_edge(v5, v0, true, b).unwrap();
    graph.add_edge(v5, v6, true, a).unwrap();
    graph.add_edge(v6, v10, false, eplus).unwrap();
    graph.add_edge(v6, v11, false, eminus).unwrap();
    graph.add_edge(v8, v7, false, eplus).unwrap();
    graph.add_edge(v9, v7, false, eminus).unwrap();

    let (n_unresolved, unresolved_type) = filters.unresolved_cut_content(&model);
    assert!(!filters.contains_cut(&model, &graph, n_unresolved, &unresolved_type));

    let mut double_double_triangle = Graph::new();
    let v0 = double_double_triangle.add_node(e1.clone());
    let v1 = double_double_triangle.add_node(e2.clone());
    let v2 = double_double_triangle.add_node(e3.clone());
    let v3 = double_double_triangle.add_node(e4.clone());

    let v4 = double_double_triangle.add_node(bba.clone());
    let v5 = double_double_triangle.add_node(bba.clone());
    let v6 = double_double_triangle.add_node(bbH.clone());
    let v7 = double_double_triangle.add_node(bbH.clone());
    let v8 = double_double_triangle.add_node(bbg.clone());
    let v9 = double_double_triangle.add_node(bbg.clone());
    let v10 = double_double_triangle.add_node(bbg.clone());
    let v11 = double_double_triangle.add_node(bbg.clone());
    let v12 = double_double_triangle.add_node(epema.clone());
    let v13 = double_double_triangle.add_node(epema.clone());

    double_double_triangle
        .add_edge(v0, v13, true, eplus)
        .unwrap();
    double_double_triangle
        .add_edge(v1, v13, true, eminus)
        .unwrap();
    double_double_triangle.add_edge(v4, v13, false, a).unwrap();
    double_double_triangle
        .add_edge(v4, v11, true, bbar)
        .unwrap();
    double_double_triangle.add_edge(v4, v10, true, b).unwrap();
    double_double_triangle.add_edge(v10, v11, false, g).unwrap();
    double_double_triangle.add_edge(v6, v11, true, b).unwrap();
    double_double_triangle
        .add_edge(v6, v10, true, bbar)
        .unwrap();
    double_double_triangle.add_edge(v6, v7, true, h).unwrap();
    double_double_triangle.add_edge(v7, v8, true, bbar).unwrap();
    double_double_triangle.add_edge(v7, v9, true, b).unwrap();
    double_double_triangle.add_edge(v5, v8, true, b).unwrap();
    double_double_triangle.add_edge(v5, v9, true, bbar).unwrap();
    double_double_triangle.add_edge(v8, v9, true, g).unwrap();
    double_double_triangle.add_edge(v5, v12, true, a).unwrap();
    double_double_triangle
        .add_edge(v12, v2, true, eminus)
        .unwrap();
    double_double_triangle
        .add_edge(v12, v3, true, eplus)
        .unwrap();

    let mut coupling = HashMap::new();
    coupling.insert("QED".into(), 6);
    let mut pert = HashMap::new();
    pert.insert("QCD".into(), 2);
    let filters = FeynGen::new(FeynGenOptions {
        generation_type: GenerationType::CrossSection,
        initial_pdgs: vec![-11, 11],
        final_pdgs: vec![5, -5, 25],
        loop_count_range: (4, 4),
        symmetrize_initial_states: true,
        symmetrize_final_states: true,
        symmetrize_left_right_states: true,
        allow_symmetrization_of_external_fermions_in_amplitudes: false,
        max_multiplicity_for_fast_cut_filter: 6,
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
            FeynGenFilter::LoopCountRange((4, 4)),
        ]),
    });

    let (n_unresolved, unresolved_type) = filters.unresolved_cut_content(&model);
    assert!(!filters.contains_cut(
        &model,
        &double_double_triangle,
        n_unresolved,
        &unresolved_type
    ));
}
