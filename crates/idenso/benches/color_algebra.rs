use std::hint::black_box;

use gungraun::{LibraryBenchmarkConfig, prelude::*};
use idenso::reference_cases::{ReferenceCase, ReferenceDomain, reference_cases};
use symbolica::atom::Atom;

fn enabled_color_cases() -> Vec<(&'static ReferenceCase, Atom)> {
    reference_cases()
        .iter()
        .filter(|case| case.domain == ReferenceDomain::ColorSu3 && case.validation.is_enabled())
        .map(|case| (case, case.expression()))
        .collect()
}

#[library_benchmark]
#[bench::reference_suite(setup = enabled_color_cases)]
fn reference_suite(cases: Vec<(&'static ReferenceCase, Atom)>) -> Vec<Atom> {
    cases
        .into_iter()
        .map(|(case, expression)| case.simplify(black_box(&expression)))
        .collect()
}

library_benchmark_group!(
    name = color_algebra;
    benchmarks = reference_suite
);

main!(
    config = LibraryBenchmarkConfig::default().pass_through_envs([
        "SYMBOLICA_LICENSE",
        "SYMBOLICA_HIDE_BANNER",
    ]);
    library_benchmark_groups = color_algebra
);
