use super::*;
use insta::assert_snapshot;
use spenso::{network::library::symbolic::ETS, symbolica_atom::IntoAtom};
use symbolica::function;

// Ported from FeynCalc's framework-independent SUNSimplify tests:
// https://github.com/FeynCalc/feyncalc/blob/027e6a741fe62a21e6d979f9d555bd2736029108/Tests/SUN/SUNSimplify.test
//
// FeynCalc normally substitutes TR=1/2 and can rewrite SUNN in terms of CA/CF.
// These tests keep idenso's symbolic TR, Nc, and NA conventions visible.
// The rest of SUNSimplify.test mostly exercises explicit SUND tensors, full
// Fierz products of two open chains, or SUNNToCACF-specific presentation rules.

fn fundamental_delta(left: impl IntoAtom, right: impl IntoAtom) -> Atom {
    function!(ETS.metric, left.into_atom(), right.into_atom())
}

macro_rules! t {
    ($r:ident, $a:tt) => {
        color_t!(slot!($r.coad_na, $a))
    };
}

macro_rules! f {
    ($r:ident, $a:tt, $b:tt, $c:tt) => {
        color_f!(
            slot!($r.coad_na, $a),
            slot!($r.coad_na, $b),
            slot!($r.coad_na, $c),
        )
    };
}

macro_rules! sun_tf {
    ($r:ident, $start:tt, $end:tt $(, $a:tt)+ $(,)?) => {
        chain!(
            slot!($r.cof_nc, $start),
            slot!($r.cof_nc.dual(), $end),
            $(t!($r, $a)),+
        )
    };
}

macro_rules! sun_trace {
    ($r:ident $(, $a:tt)+ $(,)?) => {
        trace!(& $r.cof_nc, $(t!($r, $a)),+)
    };
}

macro_rules! sdf {
    ($r:ident, $left:tt, $right:tt) => {
        fundamental_delta(slot!($r.cof_nc, $left), slot!($r.cof_nc.dual(), $right))
    };
}

#[test]
fn sun_simplify_sunn_to_cacf_rewrites_sunn_squared_minus_one() {
    initialize();
    let expr = parse_lit!(Nc ^ 2 - 1, default_namespace = "spenso");

    assert_snapshot!(expr.to_color_casimir().to_bare_ordered_string(), @"2*CA*CF");
}

#[test]
fn sun_simplify_sunn_to_cacf_rewrites_structure_square_dimension() {
    initialize();
    let expr = parse_lit!(
        f(
            coad(Nc ^ 2 - 1, 1),
            coad(Nc ^ 2 - 1, 2),
            coad(Nc ^ 2 - 1, 3)
        ) ^ 2,
        default_namespace = "spenso"
    );

    assert_snapshot!(expr
        .simplify_color()
        .to_color_casimir()
        .to_bare_ordered_string(), @"2*CA^2*CF");
}

#[test]
fn sun_simplify_id2_open_chain_separated_casimir() {
    let r = TestReps::new();
    let expr = sun_tf!(r, i, j, a, b, a);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"-1/2*CA*chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(NA,b),in,out))+CF*chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(NA,b),in,out))");
}

#[test]
fn sun_simplify_id3_structure_loop_to_adjoint_delta() {
    let r = TestReps::new();
    let expr = f!(r, a, c, d) * f!(r, b, c, d);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"CA*g(coad(NA,a),coad(NA,b))");
}

#[test]
fn sun_simplify_id4_structure_times_open_chain_contracts() {
    let r = TestReps::new();
    let expr = f!(r, a, b, c) * sun_tf!(r, i, j, b, c);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1𝑖/2*CA*chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(NA,a),in,out))");
}

#[test]
fn sun_simplify_id7_cyclic_structure_times_open_chain_contracts() {
    let r = TestReps::new();
    let expr = f!(r, c, a, b) * sun_tf!(r, i, j, b, c);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1𝑖/2*CA*chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(NA,a),in,out))");
}

#[test]
fn sun_simplify_id13_fundamental_delta_contracts() {
    let r = TestReps::new();
    let expr = sdf!(r, a, b) * sdf!(r, b, d);

    assert_snapshot!(expr.simplify_color().simplify_metrics().to_bare_ordered_string(), @"g(cof(Nc,a),dind(cof(Nc,d)))");
}

#[test]
fn sun_simplify_id17_delta_renames_open_chain_endpoint() {
    let r = TestReps::new();
    let expr = sdf!(r, b, a) * sun_tf!(r, a, d, i) * sun_tf!(r, d, c, j);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"chain(cof(Nc,b),dind(cof(Nc,c)),t(coad(NA,i),in,out),t(coad(NA,j),in,out))");
}

#[test]
fn sun_simplify_id18_delta_closes_doubled_two_generator_chain() {
    let r = TestReps::new();
    let line = sun_tf!(r, a, d, i) * sun_tf!(r, d, b, j);
    let expr = sdf!(r, b, a) * (line.clone() + line);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"2*TR*g(coad(NA,i),coad(NA,j))");
}

#[test]
fn sun_simplify_id21_adjacent_open_chains_join() {
    let r = TestReps::new();
    let expr = sun_tf!(r, a, d, i) * sun_tf!(r, d, c, j);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"chain(cof(Nc,a),dind(cof(Nc,c)),t(coad(NA,i),in,out),t(coad(NA,j),in,out))");
}

#[test]
fn sun_simplify_id27_structure_square_keeps_idenso_dimensions() {
    let r = TestReps::new();
    let expr = f!(r, a, b, c).pow(Atom::num(2));

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"CA*NA");
}

#[test]
fn sun_simplify_id30_three_generator_trace_terminal() {
    let r = TestReps::new();
    let expr = sun_trace!(r, i, j, k);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1𝑖/2*TR*f(coad(NA,i),coad(NA,j),coad(NA,k))+trace(cof(Nc),sym(t(coad(NA,i),in,out),t(coad(NA,j),in,out),t(coad(NA,k),in,out)))");
}

#[test]
fn sun_simplify_id41_nested_adjacent_open_chain_casimirs() {
    let r = TestReps::new();
    let expr = sun_tf!(r, i, j, b, a, a, b, c);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"CF^2*chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(NA,c),in,out))");
}

#[test]
fn sun_simplify_id44_adjacent_open_chain_casimir() {
    let r = TestReps::new();
    let expr = sun_tf!(r, i, j, a, a);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"CF*g(cof(Nc,i),dind(cof(Nc,j)))");
}

#[test]
fn sun_simplify_id45_open_chain_separated_casimir() {
    let r = TestReps::new();
    let expr = sun_tf!(r, i, j, a, b, a);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"-1/2*CA*chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(NA,b),in,out))+CF*chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(NA,b),in,out))");
}

#[test]
fn sun_simplify_id47_two_adjacent_open_chain_casimirs() {
    let r = TestReps::new();
    let expr = sun_tf!(r, i, j, a1, a1) * sun_tf!(r, k, l, a2, a2);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"CF^2*g(cof(Nc,i),dind(cof(Nc,j)))*g(cof(Nc,k),dind(cof(Nc,l)))");
}

#[test]
fn sun_simplify_id48_iterated_adjacent_open_chain_casimirs() {
    let r = TestReps::new();
    let expr = sun_tf!(r, i, j, b, a, a, b);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"CF^2*g(cof(Nc,i),dind(cof(Nc,j)))");
}

#[test]
fn sun_simplify_id52_repeated_trace_pair() {
    let r = TestReps::new();
    let expr = sun_trace!(r, i1, i2, i1, i2);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"-1/2*CA*NA*TR+CF*NA*TR");
}

#[test]
fn sun_simplify_id55_adjacent_casimir_inside_open_chain() {
    let r = TestReps::new();
    let expr = sun_tf!(r, i, j, b, a, a, c);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"CF*chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(NA,b),in,out),t(coad(NA,c),in,out))");
}

#[test]
fn sun_simplify_id56_explicit_i_times_structure_chain() {
    let r = TestReps::new();
    let expr = Atom::i() * f!(r, b, jj, c) * sun_tf!(r, i, j, jj, c);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"-1/2*CA*chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(NA,b),in,out))");
}

#[test]
fn sun_simplify_id66_trace_adjacent_pairs() {
    let r = TestReps::new();
    let expr = sun_trace!(r, a, a, b, b);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"CF*NA*TR");
}

#[test]
fn sun_simplify_id75_structure_times_two_generator_open_chain() {
    let r = TestReps::new();
    let expr = f!(r, a, b, c) * sun_tf!(r, i, j, b, c);

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1𝑖/2*CA*chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(NA,a),in,out))");
}
