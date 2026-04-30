use idenso::schoonschip::{Schoonschip, SchoonschipSettings};
use spenso::{
    network::tags::SPENSO_TAG as T,
    structure::{
        abstract_index::AbstractIndex,
        representation::{Minkowski, RepName, Representation},
        slot::IsAbstractSlot,
    },
};
use symbolica::{atom::Atom, function, symbol};

pub fn nested_dot_expression() -> Atom {
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let q = T.rank_one_tensor_symbol("Q");

    let p1 = function!(p, 1, mink.slot::<AbstractIndex, _>(1).to_atom());
    let p2 = function!(p, 2, mink.slot::<AbstractIndex, _>(1).to_atom());
    let p2_2 = function!(p, 2, mink.slot::<AbstractIndex, _>(2).to_atom());

    let q2 = function!(
        q,
        2,
        symbol!("bla"),
        mink.slot::<AbstractIndex, _>(1).to_atom()
    );
    let q2_2 = function!(
        q,
        2,
        symbol!("bla"),
        mink.slot::<AbstractIndex, _>(2).to_atom()
    );
    let q3_2 = function!(q, 3, mink.slot::<AbstractIndex, _>(2).to_atom());

    &p1 * (&q2 + &p2 * (&q3_2 * &q2_2 + &p2_2 * &q2_2))
}

pub fn run_schoonschip(expr: Atom, settings: &SchoonschipSettings) -> Atom {
    expr.schoonschip_with_net::<false, true, AbstractIndex>(settings)
}

pub fn assert_benchmark_outputs_match() {
    let expr = nested_dot_expression();
    let expected = run_schoonschip(expr.clone(), &SchoonschipSettings::full());

    for (name, settings) in benchmark_settings() {
        let result = run_schoonschip(expr.clone(), &settings);
        assert_eq!(
            result, expected,
            "benchmark mode {name} produced a different output"
        );
    }
}

#[allow(dead_code)]
pub fn checked_nested_dot_expression() -> Atom {
    assert_benchmark_outputs_match();
    nested_dot_expression()
}

pub fn benchmark_settings() -> Vec<(&'static str, SchoonschipSettings)> {
    vec![
        ("depth_first_depth_1", SchoonschipSettings::partial()),
        (
            "breadth_first_depth_1",
            SchoonschipSettings::breadth_first(Some(1)),
        ),
        ("full", SchoonschipSettings::full()),
    ]
}
