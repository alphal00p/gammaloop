mod common;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};

fn vertex_algebra_time(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("vertex_algebra_time");

    group.bench_function("network_vertex_substitution_5", |bench| {
        bench.iter_batched(
            common::network_vertex_fixture_5,
            common::network_vertex_substitution,
            BatchSize::SmallInput,
        );
    });

    group.bench_function("network_full_algebra_5", |bench| {
        bench.iter_batched(
            common::network_vertex_fixture_5,
            common::network_full_algebra_5,
            BatchSize::SmallInput,
        );
    });

    group.bench_function("network_vertex_substitution_8", |bench| {
        bench.iter_batched(
            common::network_vertex_fixture_8,
            common::network_vertex_substitution,
            BatchSize::SmallInput,
        );
    });

    group.bench_function("network_normalize_substituted_5", |bench| {
        bench.iter_batched(
            common::network_substituted_5,
            common::network_normalize_substituted,
            BatchSize::SmallInput,
        );
    });

    group.bench_function("network_normalize_substituted_8", |bench| {
        bench.iter_batched(
            common::network_substituted_8,
            common::network_normalize_substituted,
            BatchSize::SmallInput,
        );
    });

    group.bench_function("network_schoonschip_substituted_5", |bench| {
        bench.iter_batched(
            common::network_substituted_5,
            common::network_schoonschip_substituted,
            BatchSize::SmallInput,
        );
    });

    group.bench_function("network_parse_normalized_5", |bench| {
        bench.iter_batched(
            || common::network_normalize_substituted(common::network_substituted_5()),
            common::network_parse_normalized,
            BatchSize::SmallInput,
        );
    });

    group.bench_function("network_parse_normalized_8", |bench| {
        bench.iter_batched(
            || common::network_normalize_substituted(common::network_substituted_8()),
            common::network_parse_normalized,
            BatchSize::SmallInput,
        );
    });

    group.bench_function("bare_vertex_substitution_5", |bench| {
        bench.iter_batched(
            common::bare_vertex_fixture_5,
            common::bare_vertex_substitution,
            BatchSize::SmallInput,
        );
    });

    group.bench_function("bare_vertex_substitution_8", |bench| {
        bench.iter_batched(
            common::bare_vertex_fixture_8,
            common::bare_vertex_substitution,
            BatchSize::SmallInput,
        );
    });

    group.finish();
}

criterion_group!(benches, vertex_algebra_time);
criterion_main!(benches);
