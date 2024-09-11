use std::{env, path::PathBuf, time::Duration};

use _gammaloop::{
    gammaloop_integrand::DefaultSample,
    graph::{self, Graph},
    momentum::{FourMomentum, Helicity, ThreeMomentum},
    numerator::ContractionSettings,
    tests::load_default_settings,
    tests_from_pytest::{kinematics_builder, load_amplitude_output},
    utils::F,
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
    env_logger::init();

    let mut group = c.benchmark_group("polarization benchmarks");
    group.measurement_time(Duration::from_secs(10));

    let (model, mut amplitude, _) =
        load_amplitude_output("TEST_AMPLITUDE_physical_1L_6photons/GL_OUTPUT", true);
    let graph = amplitude.amplitude_graphs[0].graph.clone();

    let n_indep_externals = 5;
    let n_loops = 1;
    let mut external_moms = vec![];

    for i in 0..n_indep_externals {
        external_moms.push(FourMomentum::from_args(
            F(i as f64),
            F(i as f64 + 0.25),
            F(i as f64 + 0.5),
            F(i as f64 + 0.75),
        ));
    }

    let mut loop_moms = vec![];

    for i in n_indep_externals..n_indep_externals + n_loops {
        loop_moms.push(ThreeMomentum::new(
            F(i as f64),
            F(i as f64 + 0.33),
            F(i as f64 + 0.66),
        ));
    }

    let jacobian = F(1.0);

    let helicities = vec![Helicity::Plus; n_indep_externals + 1];
    println!("starting benchmark");
    group.bench_function("sample generation", |b| {
        b.iter(|| {
            DefaultSample::new(
                loop_moms.clone(),
                external_moms.clone(),
                jacobian,
                &helicities,
                &graph.bare_graph,
            );
        })
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default().with_profiler(PProfProfiler::new(100, Output::Flamegraph(None)));
    targets = criterion_benchmark
}
criterion_main!(benches);
