use std::collections::BTreeMap;

use feynkit_generator::{
    GenerationFilter, GenerationOptions, Generator, Process,
    SelfEnergyFilterOptions as FeynkitSelfEnergyOptions, SnailFilterOptions as FeynkitSnailOptions,
    TadpoleFilterOptions as FeynkitTadpoleOptions,
};
use gammalooprs::{
    feyngen::{
        FeynGenFilter, FeynGenFilters, GenerationType,
        SelfEnergyFilterOptions as GammaLoopSelfEnergyOptions,
        SnailFilterOptions as GammaLoopSnailOptions,
        TadpolesFilterOptions as GammaLoopTadpoleOptions,
    },
    graph::Graph,
    model::Model as GammaLoopModel,
    processes::ProcessDefinition,
    settings::GlobalSettings,
};
use serial_test::serial;

#[path = "support/feynkit_fixtures.rs"]
mod fixtures;

const SCALAR_MODEL: &str = include_str!("../../assets/models/json/scalars/scalars_2p_3p.json");
const STANDARD_MODEL: &str = include_str!("../../assets/models/json/sm/sm.json");
const CUBIC_SCALAR_VERTEX: &str = "V_3_SCALAR_111";
const PHOTON_DOWN_VERTEX: &str = "V_71";

fn feynkit_scalar_count(filters: Vec<GenerationFilter>) -> usize {
    let model = feynkit_model::Model::from_json(SCALAR_MODEL).unwrap();
    let process = Process::amplitude([1001_i64, 1001], [1001_i64, 1001])
        .with_loop_count(1, 1)
        .unwrap();
    let options = filters.into_iter().fold(
        GenerationOptions::default()
            .max_vertices(4)
            .with_graph_filter(GenerationFilter::VertexAllow(vec![
                CUBIC_SCALAR_VERTEX.to_owned(),
            ])),
        GenerationOptions::with_graph_filter,
    );
    Generator::new(model)
        .generate(&process, &options)
        .unwrap()
        .diagrams
        .len()
}

fn gammaloop_scalar_graphs(extra_filters: Vec<FeynGenFilter>) -> Vec<Graph> {
    let model = GammaLoopModel::from_str(SCALAR_MODEL.to_owned(), "json").unwrap();
    let mut filters = vec![FeynGenFilter::VertexAllow(vec![
        CUBIC_SCALAR_VERTEX.to_owned(),
    ])];
    filters.extend(extra_filters);
    ProcessDefinition {
        generation_type: GenerationType::Amplitude,
        initial_pdgs: vec![1001, 1001],
        final_pdgs_lists: vec![vec![1001, 1001]],
        loop_count_range: (1, 1),
        amplitude_filters: FeynGenFilters(filters),
        ..ProcessDefinition::default()
    }
    .generate(&model, &GlobalSettings::default())
    .unwrap()
}

#[test]
#[serial]
fn special_topology_filters_match_fixed_gammaloop_fixtures() {
    let cases = [
        (
            GenerationFilter::SelfEnergy(FeynkitSelfEnergyOptions::default()),
            FeynGenFilter::SelfEnergyFilter(GammaLoopSelfEnergyOptions::default()),
            fixtures::SCALAR_ONE_LOOP_FILTERED,
        ),
        (
            GenerationFilter::Tadpoles(FeynkitTadpoleOptions::default()),
            FeynGenFilter::TadpolesFilter(GammaLoopTadpoleOptions::default()),
            fixtures::SCALAR_ONE_LOOP,
        ),
        (
            GenerationFilter::ZeroSnails(FeynkitSnailOptions {
                veto_attached_to_massive: true,
                veto_attached_to_massless: true,
                only_scaleless: false,
            }),
            FeynGenFilter::ZeroSnailsFilter(GammaLoopSnailOptions {
                veto_snails_attached_to_massive_lines: true,
                veto_snails_attached_to_massless_lines: true,
                veto_only_scaleless_snails: false,
            }),
            fixtures::SCALAR_ONE_LOOP,
        ),
        (
            GenerationFilter::FactorizedLoopTopologiesCountRange((0, 0)),
            FeynGenFilter::FactorizedLoopTopologiesCountRange((0, 0)),
            fixtures::SCALAR_ONE_LOOP_FILTERED,
        ),
    ];

    let baseline = feynkit_scalar_count(Vec::new());
    let gammaloop_baseline = gammaloop_scalar_graphs(Vec::new());
    fixtures::assert_graph_fixture(&gammaloop_baseline, fixtures::SCALAR_ONE_LOOP);
    assert_eq!(baseline, fixtures::SCALAR_ONE_LOOP.len());
    let mut changed = false;
    for (feynkit_filter, gammaloop_filter, expected) in cases {
        let feynkit = feynkit_scalar_count(vec![feynkit_filter]);
        let gammaloop = gammaloop_scalar_graphs(vec![gammaloop_filter]);
        fixtures::assert_graph_fixture(&gammaloop, expected);
        assert_eq!(feynkit, expected.len());
        changed |= feynkit != baseline;
    }
    assert!(
        changed,
        "the parity fixture must exercise an effective veto"
    );
}

fn feynkit_cross_section_count(
    coupling_order: usize,
    amplitude_loops: (usize, usize),
    blob_range: std::ops::RangeInclusive<usize>,
    spectator_range: std::ops::RangeInclusive<usize>,
) -> usize {
    let model = feynkit_model::Model::from_json(STANDARD_MODEL).unwrap();
    let process = Process::cross_section([22_i64], [1_i64, -1])
        .with_loop_count(1, 1)
        .unwrap()
        .symmetrize_left_right(true);
    let options = GenerationOptions::default()
        .max_vertices(2)
        .with_graph_filter(GenerationFilter::VertexAllow(vec![
            PHOTON_DOWN_VERTEX.to_owned(),
        ]))
        .with_graph_filter(GenerationFilter::BlobRange(blob_range))
        .with_graph_filter(GenerationFilter::SpectatorRange(spectator_range))
        .with_cut_amplitude_filter(GenerationFilter::CouplingOrders(BTreeMap::from([(
            "QED".to_owned(),
            (coupling_order, Some(coupling_order)),
        )])))
        .with_cut_amplitude_filter(GenerationFilter::LoopCountRange(amplitude_loops));
    Generator::new(model)
        .generate(&process, &options)
        .unwrap()
        .diagrams
        .len()
}

fn gammaloop_cross_section_graphs(
    coupling_order: usize,
    amplitude_loops: (usize, usize),
    blob_range: std::ops::RangeInclusive<usize>,
    spectator_range: std::ops::RangeInclusive<usize>,
) -> Vec<Graph> {
    let model = GammaLoopModel::from_str(STANDARD_MODEL.to_owned(), "json").unwrap();
    let mut coupling_filter = FeynGenFilter::CouplingOrders(Default::default());
    let FeynGenFilter::CouplingOrders(coupling_orders) = &mut coupling_filter else {
        unreachable!()
    };
    coupling_orders.insert("QED".to_owned(), (coupling_order, Some(coupling_order)));
    let amplitude_filters = FeynGenFilters(vec![
        coupling_filter,
        FeynGenFilter::LoopCountRange(amplitude_loops),
    ]);
    let cross_section_filters = FeynGenFilters(vec![
        FeynGenFilter::VertexAllow(vec![PHOTON_DOWN_VERTEX.to_owned()]),
        FeynGenFilter::BlobRange(blob_range),
        FeynGenFilter::SpectatorRange(spectator_range),
    ]);
    ProcessDefinition {
        generation_type: GenerationType::CrossSection,
        initial_pdgs: vec![22],
        final_pdgs_lists: vec![vec![1, -1]],
        loop_count_range: (1, 1),
        symmetrize_left_right_states: true,
        amplitude_filters,
        cross_section_filters,
        ..ProcessDefinition::default()
    }
    .generate(&model, &GlobalSettings::default())
    .unwrap()
}

#[test]
#[serial]
fn retained_cut_partition_filters_match_fixed_gammaloop_fixtures() {
    let cases = [
        (
            1,
            (0, 0),
            1..=1,
            0..=0,
            fixtures::STANDARD_MODEL_CROSS_SECTION,
        ),
        (2, (0, 0), 1..=1, 0..=0, fixtures::EMPTY),
        (1, (1, 1), 1..=1, 0..=0, fixtures::EMPTY),
        (1, (0, 0), 0..=0, 0..=0, fixtures::EMPTY),
        (1, (0, 0), 1..=1, 1..=1, fixtures::EMPTY),
    ];
    let mut saw_retained = false;
    let mut saw_rejected = false;
    for (coupling_order, amplitude_loops, blob_range, spectator_range, expected) in cases {
        let feynkit = feynkit_cross_section_count(
            coupling_order,
            amplitude_loops,
            blob_range.clone(),
            spectator_range.clone(),
        );
        let gammaloop = gammaloop_cross_section_graphs(
            coupling_order,
            amplitude_loops,
            blob_range,
            spectator_range,
        );
        fixtures::assert_graph_fixture(&gammaloop, expected);
        assert_eq!(feynkit, expected.len());
        saw_retained |= feynkit > 0;
        saw_rejected |= feynkit == 0;
    }
    assert!(saw_retained && saw_rejected);
}
