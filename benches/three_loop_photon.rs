use std::path::PathBuf;

use _gammaloop::{
    graph::Graph,
    tests_from_pytest::{kinematics_builder, load_amplitude_output},
};
use criterion::{criterion_group, criterion_main, Criterion};
const COMPILED_DUMP: &str = "TMP_COMPILED";

fn load_helper(path: &str) -> Graph {
    let (model, mut amplitude) = load_amplitude_output(path, true);
    amplitude.amplitude_graphs[0].graph.generate_cff();

    amplitude.amplitude_graphs[0].graph.generate_numerator();
    amplitude.amplitude_graphs[0]
        .graph
        .process_numerator(&model);
    let true_path = PathBuf::from(COMPILED_DUMP).join(path);
    amplitude.amplitude_graphs[0]
        .graph
        .build_compiled_expression(true_path, true, false)
        .unwrap();

    amplitude.amplitude_graphs.remove(0).graph
}

fn criterion_benchmark(c: &mut Criterion) {
    // let _ = symbolica::LicenseManager::set_license_key("GAMMALOOP_USER");
    env_logger::init();

    let mut group = c.benchmark_group("3L physical benchmarks");

    let mut three_loop_graph =
        load_helper("TEST_AMPLITUDE_physical_3L_6photons_topology_A/GL_OUTPUT");
    let three_loop_sample = kinematics_builder(5, 3);

    group.bench_function("Inspect", |b| {
        b.iter(|| three_loop_graph.evaluate_cff_expression(&three_loop_sample, 0))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
