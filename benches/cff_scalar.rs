use _gammaloop::{
    gammaloop_integrand::DefaultSample,
    graph::Graph,
    momentum::{FourMomentum, ThreeMomentum},
    tests::load_default_settings,
    tests_from_pytest::load_amplitude_output,
    utils::F,
    ExportSettings, GammaloopCompileOptions, TropicalSubgraphTableSettings,
};
use criterion::{criterion_group, criterion_main, Criterion};
use rand::Rng;
use std::{env, path::PathBuf, time::Duration};
const COMPILED_DUMP: &str = "TMP_COMPILED";

fn kinematics_builder(n_indep_externals: usize, n_loops: usize) -> DefaultSample<f64> {
    let mut external_moms = vec![];
    let mut rng = rand::thread_rng();

    for i in 0..n_indep_externals {
        external_moms.push(FourMomentum::from_args(
            F(i as f64),
            F(i as f64 + rng.gen::<f64>() * 0.25),
            F(i as f64 + rng.gen::<f64>() * 0.5),
            F(i as f64 + rng.gen::<f64>() * 0.75),
        ));
    }

    let mut loop_moms = vec![];

    for i in n_indep_externals..n_indep_externals + n_loops {
        loop_moms.push(ThreeMomentum::new(
            F(i as f64 * rng.gen::<f64>()),
            F(i as f64 + 0.33 * rng.gen::<f64>()),
            F(i as f64 + 0.66 * rng.gen::<f64>()),
        ));
    }

    let jacobian = F(1.0);

    DefaultSample {
        loop_moms,
        external_moms,
        jacobian,
    }
}

fn load_helper(path: &str, use_orientations: bool) -> Graph {
    let (_, mut amplitude) = load_amplitude_output(path, true);
    amplitude.amplitude_graphs[0].graph.generate_cff();

    let export_settings = ExportSettings {
        compile_cff: !use_orientations,
        cpe_rounds_cff: 1,
        compile_separate_orientations: use_orientations,
        tropical_subgraph_table_settings: TropicalSubgraphTableSettings {
            target_omega: 1.0,
            panic_on_fail: false,
        },
        gammaloop_compile_options: GammaloopCompileOptions {
            inline_asm: env::var("USE_ASM").is_ok(),
            optimization_level: 3,
            fast_math: true,
            unsafe_math: true,
            compiler: "g++".to_string(),
            custom: vec![],
        },
    };

    let true_path = PathBuf::from(COMPILED_DUMP).join(path);
    amplitude.amplitude_graphs[0]
        .graph
        .build_compiled_expression(true_path, &export_settings)
        .unwrap();

    amplitude.amplitude_graphs.remove(0).graph
}

fn criterion_benchmark(c: &mut Criterion) {
    let _ = symbolica::LicenseManager::set_license_key("GAMMALOOP_USER");
    let settings = load_default_settings();

    let mut group = c.benchmark_group("scalar cff benchmarks");
    group.measurement_time(Duration::from_secs(10));

    let triangle_graph = load_helper("TEST_AMPLITUDE_massless_scalar_triangle/GL_OUTPUT", false);

    group.bench_function("Triangle", |b| {
        b.iter_batched_ref(
            || kinematics_builder(2, 1),
            |sample| triangle_graph.evaluate_cff_expression(sample, &settings),
            criterion::BatchSize::SmallInput,
        )
    });

    let box_graph = load_helper("TEST_AMPLITUDE_scalar_massless_box/GL_OUTPUT", false);

    group.bench_function("Box", |b| {
        b.iter_batched_ref(
            || kinematics_builder(3, 1),
            |sample| box_graph.evaluate_cff_expression(sample, &settings),
            criterion::BatchSize::SmallInput,
        )
    });

    let double_triangle_graph =
        load_helper("TEST_AMPLITUDE_scalar_double_triangle/GL_OUTPUT", false);

    group.bench_function("Double Triangle", |b| {
        b.iter_batched_ref(
            || kinematics_builder(1, 2),
            |sample| double_triangle_graph.evaluate_cff_expression(sample, &settings),
            criterion::BatchSize::SmallInput,
        )
    });

    let isopod_graph = load_helper("TEST_AMPLITUDE_scalar_isopod/GL_OUTPUT", false);

    group.bench_function("Isopod (Triangle-Box-Box)", |b| {
        b.iter_batched_ref(
            || kinematics_builder(2, 3),
            |sample| isopod_graph.evaluate_cff_expression(sample, &settings),
            criterion::BatchSize::SmallInput,
        )
    });

    let fishnet_2x2_graph = load_helper("TEST_AMPLITUDE_scalar_fishnet_2x2/GL_OUTPUT", false);

    group.bench_function("Fishnet 2x2", |b| {
        b.iter_batched_ref(
            || kinematics_builder(3, 4),
            |sample| fishnet_2x2_graph.evaluate_cff_expression(sample, &settings),
            criterion::BatchSize::SmallInput,
        )
    });

    let fishnet_2x3_graph = load_helper("TEST_AMPLITUDE_scalar_fishnet_2x3/GL_OUTPUT", true);

    group.bench_function("Fishnet 2x3", |b| {
        b.iter_batched_ref(
            || (kinematics_builder(3, 6)),
            |sample| fishnet_2x3_graph.evaluate_cff_expression(sample, &settings),
            criterion::BatchSize::SmallInput,
        )
    });

    std::fs::remove_dir_all(COMPILED_DUMP).unwrap();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
