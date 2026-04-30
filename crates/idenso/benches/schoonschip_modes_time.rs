mod common;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use idenso::schoonschip::SchoonschipSettings;

fn bench_mode(
    group: &mut criterion::BenchmarkGroup<'_, criterion::measurement::WallTime>,
    name: &str,
    settings: SchoonschipSettings,
) {
    group.bench_function(name, |bench| {
        bench.iter_batched(
            common::nested_dot_expression,
            |expr| common::run_schoonschip(expr, &settings),
            BatchSize::SmallInput,
        );
    });
}

fn schoonschip_modes_time(criterion: &mut Criterion) {
    common::assert_benchmark_outputs_match();

    let mut group = criterion.benchmark_group("schoonschip_modes_time");

    for (name, settings) in common::benchmark_settings() {
        bench_mode(&mut group, name, settings);
    }

    group.finish();
}

criterion_group!(benches, schoonschip_modes_time);
criterion_main!(benches);
