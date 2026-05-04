use std::hint::black_box;

mod common;

use gungraun::{LibraryBenchmarkConfig, prelude::*};
use idenso::tensor::SymbolicNet;
use spenso::structure::abstract_index::AbstractIndex;
use symbolica::atom::Atom;

#[library_benchmark]
#[bench::network_5(setup = common::network_vertex_fixture_5)]
fn network_vertex_substitution_5(fixture: common::NetworkVertexFixture) -> Atom {
    common::network_vertex_substitution(black_box(fixture))
}

#[library_benchmark]
#[bench::network_5(setup = common::network_vertex_fixture_5)]
fn network_full_algebra_5(fixture: common::NetworkVertexFixture) -> Atom {
    common::network_full_algebra_5(black_box(fixture))
}

#[library_benchmark]
#[bench::network_8(setup = common::network_vertex_fixture_8)]
fn network_vertex_substitution_8(fixture: common::NetworkVertexFixture) -> Atom {
    common::network_vertex_substitution(black_box(fixture))
}

#[library_benchmark]
#[bench::network_5(setup = common::network_substituted_5)]
fn network_normalize_substituted_5(expr: Atom) -> Atom {
    common::network_normalize_substituted(black_box(expr))
}

#[library_benchmark]
#[bench::network_8(setup = common::network_substituted_8)]
fn network_normalize_substituted_8(expr: Atom) -> Atom {
    common::network_normalize_substituted(black_box(expr))
}

#[library_benchmark]
#[bench::network_5(setup = common::network_substituted_5)]
fn network_schoonschip_substituted_5(expr: Atom) -> Atom {
    common::network_schoonschip_substituted(black_box(expr))
}

#[library_benchmark]
#[bench::network_5(setup = setup_network_normalized_5)]
fn network_parse_normalized_5(expr: Atom) -> SymbolicNet<AbstractIndex> {
    common::network_parse_normalized(black_box(expr))
}

fn setup_network_normalized_5() -> Atom {
    common::network_normalize_substituted(common::network_substituted_5())
}

#[library_benchmark]
#[bench::network_8(setup = setup_network_normalized_8)]
fn network_parse_normalized_8(expr: Atom) -> SymbolicNet<AbstractIndex> {
    common::network_parse_normalized(black_box(expr))
}

fn setup_network_normalized_8() -> Atom {
    common::network_normalize_substituted(common::network_substituted_8())
}

#[library_benchmark]
#[bench::bare_5(setup = common::bare_vertex_fixture_5)]
fn bare_vertex_substitution_5(fixture: common::BareVertexFixture) -> Atom {
    common::bare_vertex_substitution(black_box(fixture))
}

#[library_benchmark]
#[bench::bare_8(setup = common::bare_vertex_fixture_8)]
fn bare_vertex_substitution_8(fixture: common::BareVertexFixture) -> Atom {
    common::bare_vertex_substitution(black_box(fixture))
}

library_benchmark_group!(
    name = vertex_algebra;
    benchmarks =
        network_vertex_substitution_5,
        network_full_algebra_5,
        network_vertex_substitution_8,
        network_normalize_substituted_5,
        network_normalize_substituted_8,
        network_schoonschip_substituted_5,
        network_parse_normalized_5,
        network_parse_normalized_8,
        bare_vertex_substitution_5,
        bare_vertex_substitution_8
);

main!(
    config = LibraryBenchmarkConfig::default().pass_through_envs([
        "SYMBOLICA_LICENSE",
        "SYMBOLICA_HIDE_BANNER",
    ]);
    library_benchmark_groups = vertex_algebra
);
