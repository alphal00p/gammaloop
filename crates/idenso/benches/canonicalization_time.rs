use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use idenso::IndexTooling;
use spenso::{
    bracket, slot,
    structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName},
        slot::IsAbstractSlot,
    },
    tensor, tensor_symbol,
};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};

fn canonicalization_time(criterion: &mut Criterion) {
    idenso::representations::initialize();

    let rep = Minkowski {}.new_rep(4);
    let slots = (0..27)
        .map(|index| slot!(rep, AbstractIndex::Dummy(index)))
        .collect::<Vec<_>>();

    let antisymmetric = tensor_symbol!(canonicalization_benchmark_antisymmetric; Antisymmetric);
    let f = |a: usize, b: usize, c: usize| {
        function!(
            antisymmetric,
            slots[a].to_atom(),
            slots[b].to_atom(),
            slots[c].to_atom()
        )
    };
    let k33_component = |offset: usize| {
        f(offset, offset + 1, offset + 2)
            * f(offset + 3, offset + 4, offset + 5)
            * f(offset + 6, offset + 7, offset + 8)
            * f(offset, offset + 3, offset + 6)
            * f(offset + 1, offset + 4, offset + 7)
            * f(offset + 2, offset + 5, offset + 8)
    };
    let k33 = k33_component(0);
    let three_k33 = k33.clone() * k33_component(9) * k33_component(18);

    let symmetric = tensor_symbol!(canonicalization_benchmark_symmetric; Symmetric);
    let left = tensor_symbol!(canonicalization_benchmark_left);
    let right = tensor_symbol!(canonicalization_benchmark_right);
    let branch = |a: usize, b: usize| {
        function!(symmetric, slots[a].to_atom(), slots[b].to_atom())
            * function!(left, slots[a].to_atom())
            * function!(right, slots[b].to_atom())
    };
    let repeated_sum = branch(0, 1) + branch(2, 3) + branch(4, 5) + branch(6, 7);

    let power_tensor = tensor_symbol!(canonicalization_benchmark_power; Symmetric);
    let power_base = function!(power_tensor, slots[0].to_atom(), slots[1].to_atom());

    let odd_cycle_tensor = tensor_symbol!(canonicalization_benchmark_odd_cycle; Antisymmetric);
    let odd_cycle = function!(odd_cycle_tensor, slots[0].to_atom(), slots[1].to_atom())
        * function!(odd_cycle_tensor, slots[1].to_atom(), slots[2].to_atom())
        * function!(odd_cycle_tensor, slots[2].to_atom(), slots[0].to_atom());
    let nonlinear = tensor_symbol!(canonicalization_benchmark_nonlinear);
    let nonlinear_boundary = function!(nonlinear, odd_cycle);

    let factored_tensor =
        tensor_symbol!(canonicalization_benchmark_factored_antisymmetric; Antisymmetric);
    let factored = tensor!(canonicalization_benchmark_factored_scalar)
        * (Atom::num(1)
            + function!(factored_tensor, slots[0].to_atom(), slots[1].to_atom()).pow(2));

    let product_left = tensor_symbol!(canonicalization_benchmark_product_power_left);
    let product_right = tensor_symbol!(canonicalization_benchmark_product_power_right);
    let product_base =
        function!(product_left, slots[0].to_atom()) * function!(product_right, slots[0].to_atom());

    let visible_i = slot!(rep, canonicalization_benchmark_visible_i);
    let visible_j = slot!(rep, canonicalization_benchmark_visible_j);
    let visible_sign = tensor_symbol!(canonicalization_benchmark_visible_sign; Antisymmetric);
    let visible_signed_cancellation = Atom::num(2)
        * function!(visible_sign, visible_i.to_atom(), visible_j.to_atom())
        + Atom::num(2) * function!(visible_sign, visible_j.to_atom(), visible_i.to_atom());

    let cases = vec![
        ("signed_k33", k33, None),
        ("three_signed_k33_components", three_k33, None),
        ("repeated_symmetric_sum", repeated_sum, None),
        ("positive_power_4", power_base.clone().pow(4), None),
        ("negative_power_4", power_base.clone().pow(-4), None),
        ("positive_power_24", power_base.clone().pow(24), None),
        (
            "nonlinear_function_boundary",
            nonlinear_boundary,
            Some(false),
        ),
        ("factored_product_sum", factored, Some(false)),
        ("odd_power", bracket!(-power_base).pow(3), Some(false)),
        ("product_power", bracket!(product_base).pow(2), Some(false)),
        (
            "visible_signed_cancellation",
            visible_signed_cancellation,
            Some(true),
        ),
    ];

    for (name, expression, expected_zero) in &cases {
        let canonical = expression.canonize::<AbstractIndex>(AbstractIndex::Dummy);
        if let Some(expected_zero) = expected_zero {
            assert_eq!(
                canonical.is_zero(),
                *expected_zero,
                "benchmark case {name} has the wrong zero status"
            );
        }
        assert_eq!(
            canonical.canonize::<AbstractIndex>(AbstractIndex::Dummy),
            canonical,
            "benchmark case {name} must canonicalize to a fixed point"
        );
    }

    let mut group = criterion.benchmark_group("canonicalization_time");
    group.sample_size(20);
    for (name, expression, _) in cases {
        group.bench_function(name, |bench| {
            bench.iter(|| {
                black_box(black_box(&expression).canonize::<AbstractIndex>(AbstractIndex::Dummy))
            });
        });
    }
    group.finish();
}

criterion_group!(benches, canonicalization_time);
criterion_main!(benches);
