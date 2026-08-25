use feynkit_generator::{
    GenerationFilter, GenerationOptions, GenerationResult, Generator, Process,
};
use gammalooprs::{
    feyngen::{FeynGenFilter, FeynGenFilters, GenerationType},
    model::Model as GammaLoopModel,
    processes::ProcessDefinition,
    settings::GlobalSettings,
};
use serial_test::serial;

#[path = "support/feynkit_fixtures.rs"]
mod fixtures;

const STANDARD_MODEL: &str = include_str!("../../assets/models/json/sm/sm.json");
const PHOTON_LEPTON_VERTEX: &str = "V_98";
const PHOTON_DOWN_VERTEX: &str = "V_71";
const QCD_VACUUM_VERTICES: &[&str] = &["V_35", "V_36", "V_74", "V_135", "V_137"];

fn feynkit_fingerprint(result: &GenerationResult) -> u64 {
    result
        .diagrams
        .iter()
        .flat_map(|diagram| diagram.to_json().unwrap().into_bytes())
        .fold(0xcbf29ce484222325, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x100000001b3)
        })
}

#[test]
#[serial]
fn standard_model_tree_amplitude_matches_fixed_gammaloop_fixture() {
    let allowed = vec![
        PHOTON_LEPTON_VERTEX.to_owned(),
        PHOTON_DOWN_VERTEX.to_owned(),
    ];
    let model = feynkit_model::Model::from_json(STANDARD_MODEL).unwrap();
    let process = Process::amplitude([-11_i64, 11], [1_i64, -1])
        .with_loop_count(0, 0)
        .unwrap();
    let generated = Generator::new(model)
        .generate(
            &process,
            &GenerationOptions::default()
                .max_vertices(2)
                .with_graph_filter(GenerationFilter::VertexAllow(allowed.clone())),
        )
        .unwrap();

    assert_eq!(generated.diagrams.len(), 1);
    assert_eq!(feynkit_fingerprint(&generated), 0xbde083d1713f1c25);

    let gammaloop_model = GammaLoopModel::from_str(STANDARD_MODEL.to_owned(), "json").unwrap();
    let gammaloop_graphs = ProcessDefinition {
        generation_type: GenerationType::Amplitude,
        initial_pdgs: vec![-11, 11],
        final_pdgs_lists: vec![vec![1, -1]],
        loop_count_range: (0, 0),
        amplitude_filters: FeynGenFilters(vec![FeynGenFilter::VertexAllow(allowed.clone())]),
        ..ProcessDefinition::default()
    }
    .generate(&gammaloop_model, &GlobalSettings::default())
    .unwrap();
    fixtures::assert_graph_fixture(&gammaloop_graphs, fixtures::STANDARD_MODEL_TREE);
    assert_eq!(
        generated.diagrams.len(),
        fixtures::STANDARD_MODEL_TREE.len()
    );
}

#[test]
#[serial]
fn standard_model_cross_section_matches_fixed_gammaloop_fixture() {
    let allowed = vec![PHOTON_DOWN_VERTEX.to_owned()];
    let model = feynkit_model::Model::from_json(STANDARD_MODEL).unwrap();
    let process = Process::cross_section([22_i64], [1_i64, -1])
        .with_loop_count(1, 1)
        .unwrap()
        .symmetrize_left_right(true);
    let generated = Generator::new(model)
        .generate(
            &process,
            &GenerationOptions::default()
                .max_vertices(2)
                .with_graph_filter(GenerationFilter::VertexAllow(allowed.clone())),
        )
        .unwrap();

    assert_eq!(generated.diagrams.len(), 1);
    assert_eq!(feynkit_fingerprint(&generated), 0xbf7911b97672a1f4);

    let gammaloop_model = GammaLoopModel::from_str(STANDARD_MODEL.to_owned(), "json").unwrap();
    let gammaloop_graphs = ProcessDefinition {
        generation_type: GenerationType::CrossSection,
        initial_pdgs: vec![22],
        final_pdgs_lists: vec![vec![1, -1]],
        loop_count_range: (1, 1),
        symmetrize_left_right_states: true,
        cross_section_filters: FeynGenFilters(vec![FeynGenFilter::VertexAllow(allowed.clone())]),
        ..ProcessDefinition::default()
    }
    .generate(&gammaloop_model, &GlobalSettings::default())
    .unwrap();
    fixtures::assert_graph_fixture(&gammaloop_graphs, fixtures::STANDARD_MODEL_CROSS_SECTION);
    assert_eq!(
        generated.diagrams.len(),
        fixtures::STANDARD_MODEL_CROSS_SECTION.len()
    );
}

#[test]
#[serial]
fn standard_model_two_loop_vacuum_matches_fixed_gammaloop_fixture() {
    let allowed = QCD_VACUUM_VERTICES
        .iter()
        .map(|name| (*name).to_owned())
        .collect::<Vec<_>>();
    let model = feynkit_model::Model::from_json(STANDARD_MODEL).unwrap();
    let process = Process::amplitude(Vec::<i64>::new(), Vec::<i64>::new())
        .with_loop_count(2, 2)
        .unwrap();
    let generated = Generator::new(model)
        .generate(
            &process,
            &GenerationOptions::default()
                .max_vertices(2)
                .with_graph_filter(GenerationFilter::VertexAllow(allowed.clone())),
        )
        .unwrap();

    assert_eq!(generated.diagrams.len(), 5);
    assert_eq!(feynkit_fingerprint(&generated), 0x04b131f9bedb5b6e);

    let gammaloop_model = GammaLoopModel::from_str(STANDARD_MODEL.to_owned(), "json").unwrap();
    let gammaloop_graphs = ProcessDefinition {
        generation_type: GenerationType::Amplitude,
        final_pdgs_lists: vec![vec![]],
        loop_count_range: (2, 2),
        amplitude_filters: FeynGenFilters(vec![FeynGenFilter::VertexAllow(allowed.clone())]),
        ..ProcessDefinition::default()
    }
    .generate(&gammaloop_model, &GlobalSettings::default())
    .unwrap();
    fixtures::assert_graph_fixture(&gammaloop_graphs, fixtures::STANDARD_MODEL_VACUUM);
    assert_eq!(
        generated.diagrams.len(),
        fixtures::STANDARD_MODEL_VACUUM.len()
    );
}
