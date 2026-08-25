use feynkit_generator::{GenerationFilter, GenerationOptions, Generator, Process};
use gammalooprs::{
    feyngen::{FeynGenFilter, FeynGenFilters, GenerationType},
    model::Model as GammaLoopModel,
    processes::ProcessDefinition,
    settings::GlobalSettings,
};
use serial_test::serial;
use symbolica::atom::AtomCore;

#[path = "support/feynkit_fixtures.rs"]
mod fixtures;

const STANDARD_MODEL: &str = include_str!("../../assets/models/json/sm/sm.json");
const PHOTON_FERMION_VERTICES: &[&str] = &["V_98", "V_71"];

fn feynkit_options(max_vertices: usize) -> GenerationOptions {
    GenerationOptions::default()
        .max_vertices(max_vertices)
        .with_graph_filter(GenerationFilter::VertexAllow(
            PHOTON_FERMION_VERTICES
                .iter()
                .map(|name| (*name).to_owned())
                .collect(),
        ))
}

fn gammaloop_filters() -> FeynGenFilters {
    FeynGenFilters(vec![FeynGenFilter::VertexAllow(
        PHOTON_FERMION_VERTICES
            .iter()
            .map(|name| (*name).to_owned())
            .collect(),
    )])
}

fn has_negative_sign(factor: &str, function: &str) -> bool {
    factor.contains(&format!("{function}(-1)"))
}

#[test]
#[serial]
fn amplitude_external_fermion_permutation_sign_matches_fixed_fixture() {
    let model = feynkit_model::Model::from_json(STANDARD_MODEL).unwrap();
    let gammaloop_model = GammaLoopModel::from_str(STANDARD_MODEL.to_owned(), "json").unwrap();
    let mut saw_positive = false;
    let mut saw_negative = false;

    for (index, (incoming, outgoing)) in [
        (vec![-11, 11], vec![1, -1]),
        (vec![-11, 11], vec![-1, 1]),
        (vec![11, -11], vec![1, -1]),
        (vec![11, -11], vec![-1, 1]),
    ]
    .into_iter()
    .enumerate()
    {
        let generated = Generator::new(model.clone())
            .generate(
                &Process::amplitude(incoming.clone(), outgoing.clone()),
                &feynkit_options(2),
            )
            .unwrap();
        let gammaloop = ProcessDefinition {
            generation_type: GenerationType::Amplitude,
            initial_pdgs: incoming,
            final_pdgs_lists: vec![outgoing],
            loop_count_range: (0, 0),
            amplitude_filters: gammaloop_filters(),
            ..ProcessDefinition::default()
        }
        .generate(&gammaloop_model, &GlobalSettings::default())
        .unwrap();

        assert_eq!(generated.diagrams.len(), 1);
        fixtures::assert_graph_fixture(&gammaloop, fixtures::FERMION_AMPLITUDES[index]);
        let feynkit_negative = has_negative_sign(
            generated.diagrams[0].overall_factor(),
            "ExternalFermionOrderingSign",
        );
        let gammaloop_negative = has_negative_sign(
            &gammaloop[0].overall_factor.to_canonical_string(),
            "ExternalFermionOrderingSign",
        );
        let expected_negative = [false, true, true, false][index];
        assert_eq!(feynkit_negative, expected_negative);
        assert_eq!(gammaloop_negative, expected_negative);
        saw_positive |= !feynkit_negative;
        saw_negative |= feynkit_negative;
    }

    assert!(saw_positive && saw_negative);
}

#[test]
#[serial]
fn cross_section_sewn_loop_and_antifermion_spin_sum_signs_match_fixed_fixture() {
    let model = feynkit_model::Model::from_json(STANDARD_MODEL).unwrap();
    let process = Process::cross_section([-11_i64, 11], [1_i64, -1])
        .with_loop_count(1, 1)
        .unwrap()
        .symmetrize_left_right(true);
    let generated = Generator::new(model)
        .generate(&process, &feynkit_options(4))
        .unwrap();

    let gammaloop_model = GammaLoopModel::from_str(STANDARD_MODEL.to_owned(), "json").unwrap();
    let gammaloop = ProcessDefinition {
        generation_type: GenerationType::CrossSection,
        initial_pdgs: vec![-11, 11],
        final_pdgs_lists: vec![vec![1, -1]],
        loop_count_range: (1, 1),
        symmetrize_left_right_states: true,
        cross_section_filters: gammaloop_filters(),
        ..ProcessDefinition::default()
    }
    .generate(&gammaloop_model, &GlobalSettings::default())
    .unwrap();

    fixtures::assert_graph_fixture(&gammaloop, fixtures::FERMION_CROSS_SECTION);
    assert_eq!(
        generated.diagrams.len(),
        fixtures::FERMION_CROSS_SECTION.len()
    );
    assert!(!generated.diagrams.is_empty());
    for diagram in &generated.diagrams {
        for function in ["ExternalFermionOrderingSign", "AntiFermionSpinSumSign"] {
            assert!(has_negative_sign(diagram.overall_factor(), function));
        }
    }
    let gammaloop_factor = gammaloop[0].overall_factor.to_canonical_string();
    assert!(has_negative_sign(
        &gammaloop_factor,
        "ExternalFermionOrderingSign"
    ));
    assert!(has_negative_sign(
        &gammaloop_factor,
        "AntiFermionSpinSumSign"
    ));
    assert!(generated.diagrams.iter().all(|diagram| {
        has_negative_sign(diagram.overall_factor(), "ExternalFermionOrderingSign")
            && has_negative_sign(diagram.overall_factor(), "AntiFermionSpinSumSign")
    }));
}
