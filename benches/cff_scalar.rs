use std::{env, path::PathBuf, time::Duration};

use _gammaloop::{
    gammaloop_integrand::DefaultSample,
    graph::Graph,
    momentum::{FourMomentum, ThreeMomentum},
    tests_from_pytest::load_amplitude_output,
    utils::F,
    ExportSettings, GammaloopCompileOptions,
};
use criterion::{criterion_group, criterion_main, Criterion};
use pprof::criterion::{Output, PProfProfiler};
use rand::Rng;
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
    let (model, mut amplitude) = load_amplitude_output(path, true);
    amplitude.amplitude_graphs[0].graph.generate_cff();

    let export_settings = ExportSettings {
        compile_cff: !use_orientations,
        cpe_rounds_cff: 1,
        compile_separate_orientations: use_orientations,
        gammaloop_compile_options: GammaloopCompileOptions {
            inline_asm: env::var("NO_ASM").is_err(),
            optimization_level: 3,
            fast_math: true,
            unsafe_math: true,
            compiler: "g++".to_string(),
            custom: vec![],
        },
    };

    amplitude.amplitude_graphs[0].graph.generate_numerator();
    amplitude.amplitude_graphs[0]
        .graph
        .process_numerator(&model);
    let true_path = PathBuf::from(COMPILED_DUMP).join(path);
    amplitude.amplitude_graphs[0]
        .graph
        .build_compiled_expression(true_path, &export_settings)
        .unwrap();

    amplitude.amplitude_graphs.remove(0).graph
}

fn criterion_benchmark(c: &mut Criterion) {
    // let _ = symbolica::LicenseManager::set_license_key("GAMMALOOP_USER");
    // env_logger::init();

    let mut group = c.benchmark_group("scalar cff benchmarks");
    group.measurement_time(Duration::from_secs(10));

    let mut triangle_graph =
        load_helper("TEST_AMPLITUDE_massless_scalar_triangle/GL_OUTPUT", false);
    let triangle_sample = kinematics_builder(2, 1);
    group.bench_function("Triangle", |b| {
        b.iter(|| triangle_graph.evaluate_cff_expression(&triangle_sample, 0))
    });

    let mut box_graph = load_helper("TEST_AMPLITUDE_scalar_massless_box/GL_OUTPUT", false);
    let box_sample = kinematics_builder(3, 1);
    group.bench_function("Box", |b| {
        b.iter(|| box_graph.evaluate_cff_expression(&box_sample, 0))
    });

    let mut double_triangle_graph =
        load_helper("TEST_AMPLITUDE_scalar_double_triangle/GL_OUTPUT", false);
    let double_triangle_sample = kinematics_builder(1, 2);
    group.bench_function("Double Triangle", |b| {
        b.iter(|| double_triangle_graph.evaluate_cff_expression(&double_triangle_sample, 0))
    });

    let mut isopod_graph = load_helper("TEST_AMPLITUDE_scalar_isopod/GL_OUTPUT", false);
    let isopod_sample = kinematics_builder(2, 3);

    group.bench_function("Isopod (Triangle-Box-Box)", move |b| {
        b.iter(|| isopod_graph.evaluate_cff_expression(&isopod_sample, 0))
    });

    let mut fishnet_2x2_graph = load_helper("TEST_AMPLITUDE_scalar_fishnet_2x2/GL_OUTPUT", false);
    let fishnet_2x2_sample = kinematics_builder(3, 4);
    group.bench_function("Fishnet 2x2", |b| {
        b.iter(|| fishnet_2x2_graph.evaluate_cff_expression(&fishnet_2x2_sample, 0))
    });

    group.bench_function("Fishnet 2x2 cff", |b| {
        b.iter(|| fishnet_2x2_graph.evaluate_cff_all_orientations(&fishnet_2x2_sample, 0))
    });

    group.bench_function("Fishnet 2x2 num", |b| {
        b.iter(|| fishnet_2x2_graph.evaluate_numerator_all_orientations(&fishnet_2x2_sample, 0))
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
