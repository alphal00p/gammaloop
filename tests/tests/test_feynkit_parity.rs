use feynkit_generator::{GenerationFilter, GenerationOptions, Generator, Process};
use gammalooprs::{
    feyngen::{FeynGenFilter, FeynGenFilters, GenerationType},
    model::Model as GammaLoopModel,
    processes::ProcessDefinition,
    settings::GlobalSettings,
};

#[path = "support/feynkit_fixtures.rs"]
mod fixtures;

const SCALAR_MODEL: &str = include_str!("../../assets/models/json/scalars/scalars_2p_3p.json");
const CUBIC_VERTEX: &str = "V_3_SCALAR_111";

#[test]
fn scalar_tree_generation_matches_fixed_gammaloop_fixture() {
    let gammaloop_model = GammaLoopModel::from_str(SCALAR_MODEL.to_owned(), "json").unwrap();
    let gammaloop_process = ProcessDefinition {
        generation_type: GenerationType::Amplitude,
        initial_pdgs: vec![1001, 1001],
        final_pdgs_lists: vec![vec![1001, 1001]],
        loop_count_range: (0, 0),
        amplitude_filters: FeynGenFilters(vec![FeynGenFilter::VertexAllow(vec![
            CUBIC_VERTEX.to_owned(),
        ])]),
        ..ProcessDefinition::default()
    };
    let gammaloop_graphs = gammaloop_process
        .generate(&gammaloop_model, &GlobalSettings::default())
        .unwrap();
    fixtures::assert_graph_fixture(&gammaloop_graphs, fixtures::SCALAR_TREE);

    let model = feynkit_model::Model::from_json(SCALAR_MODEL).unwrap();
    let process = Process::amplitude([1001_i64, 1001], [1001_i64, 1001])
        .with_loop_count(0, 0)
        .unwrap();
    let options = GenerationOptions::default()
        .max_vertices(2)
        .with_graph_filter(GenerationFilter::VertexAllow(vec![CUBIC_VERTEX.to_owned()]));
    let generated = Generator::new(model).generate(&process, &options).unwrap();

    assert_eq!(generated.diagrams.len(), fixtures::SCALAR_TREE.len());
    assert!(generated.diagrams.iter().all(|diagram| {
        diagram
            .vertices()
            .filter_map(|(_, vertex)| vertex.interaction.as_deref())
            .all(|interaction| interaction == CUBIC_VERTEX)
    }));
}
