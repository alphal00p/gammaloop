use std::hint::black_box;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use idenso::reference_cases::{ReferenceDomain, reference_cases};

fn gamma_algebra_time(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("gamma_algebra_time");

    for case in reference_cases()
        .iter()
        .filter(|case| case.domain == ReferenceDomain::Dirac4 && case.validation.is_enabled())
    {
        group.bench_function(case.name, |bench| {
            bench.iter_batched(
                || case.expression(),
                |expression| black_box(case.simplify(black_box(&expression))),
                BatchSize::SmallInput,
            );
        });
    }

    group.finish();
}

criterion_group!(benches, gamma_algebra_time);
criterion_main!(benches);
