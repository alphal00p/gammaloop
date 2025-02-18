use std::{env, path::PathBuf};

use _gammaloop::{
    graph::Graph,
    numerator::{
        ContractionSettings, EvaluatorOptions, GammaAlgebraMode, GlobalPrefactor,
        NumeratorCompileOptions, NumeratorEvaluatorOptions, NumeratorParseMode, NumeratorSettings,
    },
    tests::load_default_settings,
    tests_from_pytest::{kinematics_builder, load_amplitude_output},
    ExportSettings, GammaloopCompileOptions, TropicalSubgraphTableSettings,
};
use criterion::{criterion_group, criterion_main, Criterion};
use pprof::criterion::{Output, PProfProfiler};
const COMPILED_DUMP: &str = "TMP_COMPILED";

fn load_helper(path: &str) -> Graph {
    let (model, mut amplitude, _) = load_amplitude_output(path, true);

    amplitude.amplitude_graphs[0].graph.generate_cff();
    let export_settings = ExportSettings {
        compile_cff: true,
        numerator_settings: NumeratorSettings {
            eval_settings: NumeratorEvaluatorOptions::Single(EvaluatorOptions {
                compile_options: NumeratorCompileOptions::Compiled,
                cpe_rounds: Some(1),
            }),
            dump_expression: None,
            global_prefactor: GlobalPrefactor::default(),
            gamma_algebra: GammaAlgebraMode::Concrete,
            global_numerator: None,
            parse_mode: NumeratorParseMode::Direct,
        },
        cpe_rounds_cff: Some(1),
        compile_separate_orientations: false,
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
    let _ = symbolica::LicenseManager::set_license_key("GAMMALOOP_USER");
    env_logger::init();

    let mut group = c.benchmark_group("3L physical benchmarks");

    let default_settings = load_default_settings();

    let settings = &default_settings;

    let mut one_loop_graph = load_helper("TEST_AMPLITUDE_physical_1L_6photons/GL_OUTPUT");
    let one_loop_sample = kinematics_builder(5, 1, &one_loop_graph.bare_graph);

    let one_loop_quad_sample = one_loop_sample.higher_precision();

    group.bench_function("Inspect 1l", |b| {
        b.iter(|| one_loop_graph.evaluate_cff_expression(&one_loop_sample, settings))
    });

    group.bench_function("Inspect 1l quad", |b| {
        b.iter(|| one_loop_graph.evaluate_cff_expression(&one_loop_quad_sample, settings))
    });

    let mut two_loop_graph = load_helper("TEST_AMPLITUDE_physical_2L_6photons/GL_OUTPUT");
    let two_loop_sample = kinematics_builder(5, 2, &two_loop_graph.bare_graph);

    group.bench_function("Inspect 2l", |b| {
        b.iter(|| two_loop_graph.evaluate_cff_expression(&two_loop_sample, settings))
    });

    // let mut three_loop_graph =
    //     load_helper("TEST_AMPLITUDE_physical_3L_6photons_topology_A/GL_OUTPUT");
    // let three_loop_sample = kinematics_builder(5, 3);

    // group.bench_function("Inspect 3l", |b| {
    //     b.iter(|| three_loop_graph.evaluate_cff_expression(&three_loop_sample, 0))
    // });
}

// criterion_group!(benches, criterion_benchmark);

criterion_group! {
    name = benches;
    config = Criterion::default().with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)));
    targets = criterion_benchmark
}
criterion_main!(benches);
