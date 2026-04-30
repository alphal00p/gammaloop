use std::hint::black_box;

mod common;

use gungraun::prelude::*;
use idenso::schoonschip::SchoonschipSettings;
use symbolica::atom::Atom;

fn run(expr: Atom, settings: SchoonschipSettings) -> Atom {
    common::run_schoonschip(black_box(expr), black_box(&settings))
}

#[library_benchmark]
#[bench::depth_first_depth_1(setup = common::checked_nested_dot_expression)]
fn depth_first_depth_1(expr: Atom) -> Atom {
    run(expr, SchoonschipSettings::partial())
}

#[library_benchmark]
#[bench::breadth_first_depth_1(setup = common::checked_nested_dot_expression)]
fn breadth_first_depth_1(expr: Atom) -> Atom {
    run(expr, SchoonschipSettings::breadth_first(Some(1)))
}

#[library_benchmark]
#[bench::full(setup = common::checked_nested_dot_expression)]
fn full(expr: Atom) -> Atom {
    run(expr, SchoonschipSettings::full())
}

library_benchmark_group!(
    name = schoonschip_modes;
    benchmarks = depth_first_depth_1, breadth_first_depth_1, full
);

main!(library_benchmark_groups = schoonschip_modes);
