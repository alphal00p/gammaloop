use std::hint::black_box;

use criterion::{BatchSize, Criterion, criterion_group, criterion_main};
use idenso::reference_cases::young::YoungProjector;
use symbolica::atom::AtomCore;

mod young_common;

fn young_projectors_time(criterion: &mut Criterion) {
    let corpus = young_common::validated_corpus();

    let mut compilation = criterion.benchmark_group("young_projector_compilation");
    for fixture in &corpus.projectors {
        compilation.bench_function(fixture.name, |bench| {
            bench.iter_batched(
                || fixture.projector.tableau().clone(),
                |tableau| black_box(YoungProjector::new(black_box(tableau))),
                BatchSize::SmallInput,
            );
        });
    }
    compilation.finish();

    let mut expansion = criterion.benchmark_group("young_projector_expansion");
    for fixture in &corpus.projectors {
        expansion.bench_function(fixture.name, |bench| {
            bench.iter_batched(
                || fixture.arguments.clone(),
                |arguments| {
                    black_box(
                        fixture
                            .projector
                            .project(black_box(fixture.head), black_box(&arguments)),
                    )
                },
                BatchSize::SmallInput,
            );
        });
    }
    expansion.finish();

    let mut product_expansion = criterion.benchmark_group("young_riemann_product_expansion");
    for case in &corpus.product_expansion {
        product_expansion.bench_function(case.name, |bench| {
            bench.iter_batched(
                || case.expression.clone(),
                |expression| black_box(black_box(expression).expand()),
                BatchSize::SmallInput,
            );
        });
    }
    product_expansion.finish();

    let mut canonicalization = criterion.benchmark_group("young_projector_canonicalization");
    for case in &corpus.canonicalization {
        canonicalization.bench_function(case.name, |bench| {
            bench.iter_batched(
                || case.expression.clone(),
                |expression| black_box(young_common::canonicalize(black_box(expression))),
                BatchSize::SmallInput,
            );
        });
    }
    canonicalization.finish();

    let mut declared = criterion.benchmark_group("young_declared_canonicalization");
    for case in &corpus.declared_canonicalization {
        declared.bench_function(case.name, |bench| {
            bench.iter_batched(
                || case.expression.clone(),
                |expression| black_box(young_common::canonicalize(black_box(expression))),
                BatchSize::SmallInput,
            );
        });
    }
    declared.finish();
}

criterion_group!(benches, young_projectors_time);
criterion_main!(benches);
