use std::hint::black_box;
mod common;
use gungraun::prelude::*;
use idenso::tensor::SymbolicTensor;
use spenso::structure::partial::PartialStructure;
use symbolica::atom::Atom;

#[library_benchmark]
#[bench::contract(setup = common::checked_nested_dot_expression)]
fn contract(expr: Atom) -> SymbolicTensor<PartialStructure> {
    common::run_contraction(black_box(expr))
}

library_benchmark_group!(name = contraction; benchmarks = contract);
main!(library_benchmark_groups = contraction);
