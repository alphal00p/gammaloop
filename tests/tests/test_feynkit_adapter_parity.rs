use feynkit_generator::GenerationReport;
use gammalooprs::{
    feyngen::{
        FeynGenFilter, FeynGenFilters, GenerationType,
        feynkit::{FeynkitAdapterResult, FeynkitGeneratorAdapter},
    },
    model::Model,
    numerator::GlobalPrefactor,
    processes::ProcessDefinition,
    settings::GlobalSettings,
};
use serial_test::serial;
use symbolica::atom::Atom;

#[path = "support/feynkit_fixtures.rs"]
mod fixtures;

const SCALAR_MODEL: &str = include_str!("../../assets/models/json/scalars/scalars_2p_3p.json");
const STANDARD_MODEL: &str = include_str!("../../assets/models/json/sm/sm.json");

fn adapted(
    model: &Model,
    definition: &ProcessDefinition,
    settings: &GlobalSettings,
) -> FeynkitAdapterResult {
    FeynkitGeneratorAdapter::new(definition, model, settings)
        .generate()
        .unwrap()
}

fn compare_with_fixture(
    model_json: &str,
    definition: ProcessDefinition,
    expected: &[u64],
) -> GenerationReport {
    let model = Model::from_str(model_json.to_owned(), "json").unwrap();
    let settings = GlobalSettings::default();
    let adapted = adapted(&model, &definition, &settings);

    fixtures::assert_graph_fixture(&adapted.graphs, expected);
    assert_eq!(adapted.report.retained_count, expected.len());
    adapted.report
}

#[test]
#[serial]
fn scalar_amplitude_adapter_matches_fixed_fixture() {
    let report = compare_with_fixture(
        SCALAR_MODEL,
        ProcessDefinition {
            generation_type: GenerationType::Amplitude,
            initial_pdgs: vec![1001, 1001],
            final_pdgs_lists: vec![vec![1001, 1001]],
            loop_count_range: (0, 0),
            amplitude_filters: FeynGenFilters(vec![FeynGenFilter::VertexAllow(vec![
                "V_3_SCALAR_111".to_owned(),
            ])]),
            ..ProcessDefinition::default()
        },
        fixtures::SCALAR_TREE,
    );
    assert!(report.completed);
}

#[test]
#[serial]
fn standard_model_cross_section_adapter_matches_fixed_fixture() {
    let report = compare_with_fixture(
        STANDARD_MODEL,
        ProcessDefinition {
            generation_type: GenerationType::CrossSection,
            initial_pdgs: vec![22],
            final_pdgs_lists: vec![vec![1, -1]],
            loop_count_range: (1, 1),
            symmetrize_left_right_states: true,
            cross_section_filters: FeynGenFilters(vec![FeynGenFilter::VertexAllow(vec![
                "V_71".to_owned(),
            ])]),
            ..ProcessDefinition::default()
        },
        fixtures::STANDARD_MODEL_CROSS_SECTION,
    );
    assert!(report.completed);
}

#[test]
#[serial]
fn selected_graph_and_prefactor_are_gamma_loop_policy() {
    let model = Model::from_str(SCALAR_MODEL.to_owned(), "json").unwrap();
    let settings = GlobalSettings::default();
    let definition = ProcessDefinition {
        generation_type: GenerationType::Amplitude,
        initial_pdgs: vec![1001, 1001],
        final_pdgs_lists: vec![vec![1001, 1001]],
        loop_count_range: (0, 0),
        amplitude_filters: FeynGenFilters(vec![FeynGenFilter::VertexAllow(vec![
            "V_3_SCALAR_111".to_owned(),
        ])]),
        graph_prefix: "A".to_owned(),
        selected_graphs: Some(vec!["A0".to_owned()]),
        prefactor: GlobalPrefactor {
            num: Atom::num(7),
            projector: Atom::num(3),
        },
        ..ProcessDefinition::default()
    };
    let result = adapted(&model, &definition, &settings);

    fixtures::assert_graph_fixture(&result.graphs, fixtures::SCALAR_SELECTED_PREFACTOR);
    assert_eq!(result.graphs.len(), 1);
    assert_eq!(result.graphs[0].name, "A0");
    assert_eq!(result.graphs[0].global_prefactor.num, Atom::num(7));
    assert_eq!(result.graphs[0].global_prefactor.projector, Atom::num(3));
}
