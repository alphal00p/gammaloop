use std::hint::black_box;

use gungraun::{LibraryBenchmarkConfig, prelude::*};
use idenso::reference_cases::young::{YoungProjector, YoungProjectorFixture};
use spenso::structure::YoungTableau;
use symbolica::atom::{Atom, AtomCore};

mod young_common;

fn compilation_inputs() -> Vec<YoungTableau> {
    young_common::validated_corpus()
        .projectors
        .into_iter()
        .map(|fixture| fixture.projector.tableau().clone())
        .collect()
}

fn expansion_inputs() -> Vec<YoungProjectorFixture> {
    young_common::validated_corpus().projectors
}

fn product_expansion_inputs() -> Vec<Atom> {
    young_common::validated_corpus()
        .product_expansion
        .into_iter()
        .map(|case| {
            assert!(!case.name.is_empty());
            case.expression
        })
        .collect()
}

fn canonicalization_inputs() -> Vec<Atom> {
    young_common::validated_corpus()
        .canonicalization
        .into_iter()
        .map(|case| {
            assert!(!case.name.is_empty());
            case.expression
        })
        .collect()
}

fn declared_canonicalization_inputs() -> Vec<Atom> {
    young_common::validated_corpus()
        .declared_canonicalization
        .into_iter()
        .map(|case| {
            assert!(!case.name.is_empty());
            case.expression
        })
        .collect()
}

#[library_benchmark]
#[bench::compilation_suite(setup = compilation_inputs)]
fn compilation_suite(tableaux: Vec<YoungTableau>) -> Vec<YoungProjector> {
    tableaux
        .into_iter()
        .map(|tableau| YoungProjector::new(black_box(tableau)))
        .collect()
}

#[library_benchmark]
#[bench::expansion_suite(setup = expansion_inputs)]
fn expansion_suite(fixtures: Vec<YoungProjectorFixture>) -> Vec<Atom> {
    fixtures
        .into_iter()
        .map(|fixture| {
            fixture
                .projector
                .project(black_box(fixture.head), black_box(&fixture.arguments))
        })
        .collect()
}

#[library_benchmark]
#[bench::product_expansion_suite(setup = product_expansion_inputs)]
fn product_expansion_suite(expressions: Vec<Atom>) -> Vec<Atom> {
    expressions
        .into_iter()
        .map(|expression| black_box(expression).expand())
        .collect()
}

#[library_benchmark]
#[bench::canonicalization_suite(setup = canonicalization_inputs)]
fn canonicalization_suite(expressions: Vec<Atom>) -> Vec<Atom> {
    expressions
        .into_iter()
        .map(|expression| young_common::canonicalize(black_box(expression)))
        .collect()
}

#[library_benchmark]
#[bench::declared_canonicalization_suite(setup = declared_canonicalization_inputs)]
fn declared_canonicalization_suite(expressions: Vec<Atom>) -> Vec<Atom> {
    expressions
        .into_iter()
        .map(|expression| young_common::canonicalize(black_box(expression)))
        .collect()
}

library_benchmark_group!(
    name = young_projectors;
    benchmarks =
        compilation_suite, expansion_suite, product_expansion_suite, canonicalization_suite,
        declared_canonicalization_suite
);

main!(
    config = LibraryBenchmarkConfig::default().pass_through_envs([
        "SYMBOLICA_LICENSE",
        "SYMBOLICA_HIDE_BANNER",
    ]);
    library_benchmark_groups = young_projectors
);
