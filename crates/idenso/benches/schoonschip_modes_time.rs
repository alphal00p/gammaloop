mod common;
use criterion::{BatchSize, Criterion, criterion_group, criterion_main};

fn contraction_time(criterion: &mut Criterion) {
    common::assert_contraction_invariants();
    criterion.bench_function("contract", |bench| {
        bench.iter_batched(
            common::nested_dot_expression,
            common::run_contraction,
            BatchSize::SmallInput,
        );
    });
}
criterion_group!(benches, contraction_time);
criterion_main!(benches);
