use std::{env, path::PathBuf, time::Duration};

use _gammaloop::{
    graph::Graph,
    numerator::ContractionSettings,
    tests::load_default_settings,
    tests_from_pytest::{kinematics_builder, load_amplitude_output},
    ExportSettings, GammaloopCompileOptions, TropicalSubgraphTableSettings,
};
use criterion::{criterion_group, criterion_main, Criterion};
use pprof::criterion::{Output, PProfProfiler};
const COMPILED_DUMP: &str = "TMP_COMPILED";

fn load_helper(path: &str, use_orientations: bool) -> Graph {
    let (model, mut amplitude, _) = load_amplitude_output(path, true);
    amplitude.amplitude_graphs[0].graph.generate_cff();

    let export_settings = ExportSettings {
        compile_cff: !use_orientations,
        numerator_settings: Default::default(),
        cpe_rounds_cff: Some(1),
        compile_separate_orientations: use_orientations,
        tropical_subgraph_table_settings: TropicalSubgraphTableSettings {
            target_omega: 1.0,
            panic_on_fail: false,
        },
        gammaloop_compile_options: GammaloopCompileOptions {
            inline_asm: env::var("NO_ASM").is_err(),
            optimization_level: 3,
            fast_math: true,
            unsafe_math: true,
            compiler: "g++".to_string(),
            custom: vec![],
        },
    };

    let true_path = PathBuf::from(COMPILED_DUMP).join(path);

    let mut g = amplitude
        .amplitude_graphs
        .remove(0)
        .graph
        .process_numerator(
            &model,
            ContractionSettings::Normal,
            true_path.clone(),
            &export_settings,
        );
    g.build_compiled_expression(true_path, &export_settings)
        .unwrap();
    g
}

fn criterion_benchmark(c: &mut Criterion) {
    // let _ = symbolica::LicenseManager::set_license_key("GAMMALOOP_USER");
    // env_logger::init();

    let mut group = c.benchmark_group("scalar cff benchmarks");
    group.measurement_time(Duration::from_secs(10));

    let default_settings = load_default_settings();

    let settings = &default_settings;
    let mut triangle_graph =
        load_helper("TEST_AMPLITUDE_massless_scalar_triangle/GL_OUTPUT", false);
    let triangle_sample = kinematics_builder(2, 1, &triangle_graph.bare_graph);
    group.bench_function("Triangle", |b| {
        b.iter(|| triangle_graph.evaluate_cff_expression(&triangle_sample, settings))
    });

    let mut box_graph = load_helper("TEST_AMPLITUDE_scalar_massless_box/GL_OUTPUT", false);
    let box_sample = kinematics_builder(3, 1, &box_graph.bare_graph);
    group.bench_function("Box", |b| {
        b.iter(|| box_graph.evaluate_cff_expression(&box_sample, settings))
    });

    let mut double_triangle_graph =
        load_helper("TEST_AMPLITUDE_scalar_double_triangle/GL_OUTPUT", false);
    let double_triangle_sample = kinematics_builder(1, 2, &double_triangle_graph.bare_graph);
    group.bench_function("Double Triangle", |b| {
        b.iter(|| double_triangle_graph.evaluate_cff_expression(&double_triangle_sample, settings))
    });

    let mut isopod_graph = load_helper("TEST_AMPLITUDE_scalar_isopod/GL_OUTPUT", false);
    let isopod_sample = kinematics_builder(2, 3, &isopod_graph.bare_graph);

    group.bench_function("Isopod (Triangle-Box-Box)", move |b| {
        b.iter(|| isopod_graph.evaluate_cff_expression(&isopod_sample, settings))
    });

    let mut fishnet_2x2_graph = load_helper("TEST_AMPLITUDE_scalar_fishnet_2x2/GL_OUTPUT", false);
    let fishnet_2x2_sample = kinematics_builder(3, 4, &fishnet_2x2_graph.bare_graph);
    group.bench_function("Fishnet 2x2", |b| {
        b.iter(|| fishnet_2x2_graph.evaluate_cff_expression(&fishnet_2x2_sample, settings))
    });

    group.bench_function("Fishnet 2x2 cff", |b| {
        b.iter(|| fishnet_2x2_graph.evaluate_cff_all_orientations(&fishnet_2x2_sample, settings))
    });

    group.bench_function("Fishnet 2x2 num", |b| {
        b.iter(|| {
            fishnet_2x2_graph.evaluate_numerator_all_orientations(&fishnet_2x2_sample, settings)
        })
    });

    // // std::fs::remove_dir_all(COMPILED_DUMP).unwrap();
    // let mut fishnet_2x3_graph = load_helper("TEST_AMPLITUDE_scalar_fishnet_2x3/GL_OUTPUT", true);
    // let fishnet_2x3_sample = kinematics_builder(3, 6);

    // group.bench_function("Fishnet 2x3", |b| {
    //     b.iter(|| fishnet_2x3_graph.evaluate_cff_expression(&fishnet_2x3_sample, 0))
    // });
}

criterion_group! {
    name = benches;
    config = Criterion::default().with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)));
    targets = criterion_benchmark
}
criterion_main!(benches);
