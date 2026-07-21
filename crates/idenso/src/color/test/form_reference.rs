use super::*;
use insta::assert_snapshot;

macro_rules! fco {
    ($r:ident, $a:tt, $b:tt, $c:tt) => {
        color_f!(
            slot!($r.coad_na, $a),
            slot!($r.coad_na, $b),
            slot!($r.coad_na, $c),
        )
    };
}

macro_rules! tco {
    ($r:ident, $a:tt) => {
        color_t!(slot!($r.coad_na, $a))
    };
}

#[test]
fn one_generator_trace_vanishes() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, i))),
        default_namespace = "spenso"
    );

    assert_color_zero(expr);
}

#[test]
fn two_generator_trace_normalizes_to_tr_metric() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, i))),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"g(coad(-1+Nc^2,a),coad(-1+Nc^2,b))*idx(2,cof(Nc))");
}

#[test]
fn fierz_generator_contraction() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().expand().simplify_metrics().to_bare_ordered_string(), @"-1*Nc^(-1)*g(cof(Nc,i),dind(cof(Nc,j)))*g(cof(Nc,k),dind(cof(Nc,l)))*idx(2,cof(Nc))+g(cof(Nc,i),dind(cof(Nc,l)))*g(cof(Nc,k),dind(cof(Nc,j)))*idx(2,cof(Nc))");
}

#[test]
fn separated_generator_casimir_shortcut() {
    test_initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, k)))
            * t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().simplify_metrics().to_bare_ordered_string(), @"-1/2*cas(2,coad(-1+Nc^2))*t(coad(-1+Nc^2,b),cof(Nc,i),dind(cof(Nc,l)))+cas(2,cof(Nc))*t(coad(-1+Nc^2,b),cof(Nc,i),dind(cof(Nc,l)))");
}

#[test]
fn chain_one_generator_trace_vanishes() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_na, a)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"0");
}

#[test]
fn chain_two_generator_trace_normalizes() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"g(coad(NA,a),coad(NA,b))*idx(2,cof(Nc))");
}

#[test]
fn three_generator_trace_terminal() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, c)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1𝑖/2*f(coad(NA,a),coad(NA,b),coad(NA,c))*idx(2,cof(Nc))+trace(cof(Nc),sym(t(coad(NA,a),in,out),t(coad(NA,b),in,out),t(coad(NA,c),in,out)))");
}

#[test]
fn adjacent_generator_casimir_chain() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), k),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, a)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"cas(2,cof(Nc))*g(cof(Nc,i),dind(cof(Nc,k)))");
}

#[test]
fn separated_generator_casimir_trace() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, c)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"-1/2*cas(2,coad(NA))*g(coad(NA,b),coad(NA,c))*idx(2,cof(Nc))+cas(2,cof(Nc))*g(coad(NA,b),coad(NA,c))*idx(2,cof(Nc))");
}

#[test]
fn two_f_loop_contracts_to_ca_metric() {
    test_initialize();
    let expr = parse_lit!(
        f(coad(NA, a), coad(NA, c), coad(NA, x)) * f(coad(NA, b), coad(NA, c), coad(NA, x)),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"cas(2,coad(NA))*g(coad(NA,a),coad(NA,b))");
}

#[test]
fn metric_contracted_f_square_is_positive() {
    test_initialize();
    // Seed the symbol order from the three-loop graph that exposed the lost
    // permutation parity in the old sequence-wildcard contraction.
    let _ = parse_lit!(
        f(
            coad(8, metric_orientation_a),
            coad(8, metric_orientation_x),
            coad(8, metric_orientation_u)
        ) * f(
            coad(8, metric_orientation_a),
            coad(8, metric_orientation_x),
            coad(8, metric_orientation_v)
        ),
        default_namespace = "spenso"
    )
    .simplify_metrics();
    let contracted = parse_lit!(
        g(coad(8, metric_orientation_r), coad(8, metric_orientation_s))
            * f(
                coad(8, metric_orientation_u),
                coad(8, metric_orientation_j),
                coad(8, metric_orientation_r)
            ),
        default_namespace = "spenso"
    )
    .simplify_metrics();
    let expected = parse_lit!(
        f(
            coad(8, metric_orientation_u),
            coad(8, metric_orientation_j),
            coad(8, metric_orientation_s)
        ),
        default_namespace = "spenso"
    );
    assert_eq!(
        contracted.to_bare_ordered_string(),
        expected.to_bare_ordered_string()
    );

    let nested = parse_lit!(
        g(coad(8, metric_orientation_r), coad(8, metric_orientation_s))
            * metric_orientation_probe(
                f(
                    coad(8, metric_orientation_u),
                    coad(8, metric_orientation_j),
                    coad(8, metric_orientation_r)
                ),
                coad(8, metric_orientation_r)
            ),
        default_namespace = "spenso"
    )
    .simplify_metrics();
    let nested_expected = parse_lit!(
        metric_orientation_probe(
            f(
                coad(8, metric_orientation_u),
                coad(8, metric_orientation_j),
                coad(8, metric_orientation_r)
            ),
            coad(8, metric_orientation_s)
        ),
        default_namespace = "spenso"
    );
    assert_eq!(
        nested.to_bare_ordered_string(),
        nested_expected.to_bare_ordered_string()
    );

    let expr = parse_lit!(
        g(coad(8, metric_orientation_u), coad(8, metric_orientation_v))
            * g(coad(8, metric_orientation_r), coad(8, metric_orientation_s))
            * f(
                coad(8, metric_orientation_u),
                coad(8, metric_orientation_j),
                coad(8, metric_orientation_r)
            )
            * f(
                coad(8, metric_orientation_v),
                coad(8, metric_orientation_j),
                coad(8, metric_orientation_s)
            ),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_metrics().to_bare_ordered_string(), @"(f(coad(8,metric_orientation_j),coad(8,metric_orientation_s),coad(8,metric_orientation_v)))^2");
    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"8*cas(2,coad(8))");
}

#[test]
fn three_f_loop_contracts_to_ca_f() {
    test_initialize();
    let expr = parse_lit!(
        f(coad(NA, a), coad(NA, b), coad(NA, e))
            * f(coad(NA, b), coad(NA, c), coad(NA, f_))
            * f(coad(NA, c), coad(NA, a), coad(NA, g_)),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1/2*cas(2,coad(NA))*f(coad(NA,e),coad(NA,f_),coad(NA,g_))");
}

#[test]
fn three_f_loop_preserves_antisymmetric_orientation() {
    test_initialize();
    // Seed the symbol order that made the unoriented triangle reduction pick
    // the wrong sign in the three-loop graph.
    let _ = parse_lit!(
        g(
            coad(8, triangle_orientation_r),
            coad(8, triangle_orientation_s)
        ) * f(
            coad(8, triangle_orientation_x),
            coad(8, triangle_orientation_z),
            coad(8, triangle_orientation_r)
        ),
        default_namespace = "spenso"
    );
    let triangle = parse_lit!(
        f(
            coad(8, triangle_orientation_u),
            coad(8, triangle_orientation_x),
            coad(8, triangle_orientation_r)
        ) * f(
            coad(8, triangle_orientation_y),
            coad(8, triangle_orientation_z),
            coad(8, triangle_orientation_r)
        ) * f(
            coad(8, triangle_orientation_u),
            coad(8, triangle_orientation_z),
            coad(8, triangle_orientation_s)
        ),
        default_namespace = "spenso"
    );
    let closing_structure = parse_lit!(
        f(
            coad(8, triangle_orientation_x),
            coad(8, triangle_orientation_y),
            coad(8, triangle_orientation_s)
        ),
        default_namespace = "spenso"
    );

    assert_snapshot!(triangle.simplify_color().to_bare_ordered_string(), @"-1/2*cas(2,coad(8))*f(coad(8,triangle_orientation_s),coad(8,triangle_orientation_x),coad(8,triangle_orientation_y))");
    assert_snapshot!((&triangle * &closing_structure).simplify_color().to_bare_ordered_string(), @"(cas(2,coad(8)))^2*-4");
    assert_snapshot!((triangle.simplify_color() * closing_structure).simplify_color().to_bare_ordered_string(), @"(cas(2,coad(8)))^2*-4");
}

#[test]
fn four_f_closed_ghost_topology_preserves_antisymmetric_orientation() {
    test_initialize();
    let expr = parse_lit!(
        -1 / 8
            * g(coad(8, hedge(0)), coad(8, hedge(1)))
            * f(coad(8, hedge(3)), coad(8, hedge(0)), coad(8, hedge(5)))
            * f(coad(8, hedge(3)), coad(8, hedge(7)), coad(8, hedge(9)))
            * f(coad(8, hedge(7)), coad(8, hedge(11)), coad(8, hedge(1)))
            * f(coad(8, hedge(9)), coad(8, hedge(5)), coad(8, hedge(11))),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"(cas(2,coad(8)))^2*1/2");
}

#[test]
fn factored_four_f_closed_ghost_topology_preserves_antisymmetric_orientation() {
    test_initialize();
    let expr = parse_lit!(
        -1 / 8
            * g(coad(8, hedge(0)), coad(8, hedge(1)))
            * (5 / 16
                * x
                * f(coad(8, hedge(3)), coad(8, hedge(0)), coad(8, hedge(5)))
                * f(coad(8, hedge(3)), coad(8, hedge(7)), coad(8, hedge(9)))
                * f(coad(8, hedge(7)), coad(8, hedge(11)), coad(8, hedge(1)))
                * f(coad(8, hedge(9)), coad(8, hedge(5)), coad(8, hedge(11)))
                + 3 / 8
                    * y
                    * f(coad(8, hedge(3)), coad(8, hedge(0)), coad(8, hedge(5)))
                    * f(coad(8, hedge(3)), coad(8, hedge(7)), coad(8, hedge(9)))
                    * f(coad(8, hedge(7)), coad(8, hedge(11)), coad(8, hedge(1)))
                    * f(coad(8, hedge(9)), coad(8, hedge(5)), coad(8, hedge(11)))),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().expand().to_bare_ordered_string(), @"(cas(2,coad(8)))^2*3/16*y+(cas(2,coad(8)))^2*5/32*x");
}

#[test]
fn mixed_trace_structure_contraction() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, d)),
    ) * color_f!(
        slot!(r.coad_na, a),
        slot!(r.coad_na, b),
        slot!(r.coad_na, c),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1𝑖/2*cas(2,coad(NA))*g(coad(NA,c),coad(NA,d))*idx(2,cof(Nc))");
}

#[test]
fn symmetric_invariant_d33_partial_contraction() {
    test_initialize();
    let expr = parse_lit!(
        d(sym_x, coad(NA, a), coad(NA, b), coad(NA, c))
            * d(sym_y, coad(NA, a), coad(NA, b), coad(NA, d_)),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"NA^(-1)*g(coad(NA,c),coad(NA,d_))*gram(3,sym_x,sym_y)");
}

#[test]
fn four_generator_trace_terminal() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, c)),
        color_t!(slot!(r.coad_na, d)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"-1/6*f(coad(NA,a),coad(NA,c),coad(NA,x))*f(coad(NA,b),coad(NA,d),coad(NA,x))*idx(2,cof(Nc))+1/3*f(coad(NA,a),coad(NA,d),coad(NA,x))*f(coad(NA,b),coad(NA,c),coad(NA,x))*idx(2,cof(Nc))+1𝑖/2*f(coad(NA,a),coad(NA,b),coad(NA,x))*trace(cof(Nc),sym(t(coad(NA,c),in,out),t(coad(NA,d),in,out),t(coad(NA,x),in,out)))+1𝑖/2*f(coad(NA,c),coad(NA,d),coad(NA,x))*trace(cof(Nc),sym(t(coad(NA,a),in,out),t(coad(NA,b),in,out),t(coad(NA,x),in,out)))+trace(cof(Nc),sym(t(coad(NA,a),in,out),t(coad(NA,b),in,out),t(coad(NA,c),in,out),t(coad(NA,d),in,out)))");
}

// FORM's repository includes a valgrind-oriented size-5 port of color.h's
// tloop.frm check:
// https://github.com/form-dev/form/blob/fdf6af0ce520f7da2dfe0b0b61bf4cef396770be/check/extra/color.frm
//
// These ignored tests keep the upstream expressions and expected invariant
// families close to the implementation. They require higher symmetric
// invariant reductions than the current chain/trace simplifier provides.

fn form_color_tloop_q10(r: &TestReps) -> Atom {
    trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    )
}

fn form_color_tloop_g10(r: &TestReps) -> Atom {
    fco!(r, i1, i2, j1)
        * fco!(r, i2, i3, j2)
        * fco!(r, i3, i4, j3)
        * fco!(r, i4, i5, j4)
        * fco!(r, i5, i6, j5)
        * fco!(r, i6, i7, j1)
        * fco!(r, i7, i8, j2)
        * fco!(r, i8, i9, j3)
        * fco!(r, i9, i10, j4)
        * fco!(r, i10, i1, j5)
}

fn form_color_tloop_qq5(r: &TestReps) -> Atom {
    trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    ) * trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    )
}

fn form_color_tloop_qg5(r: &TestReps) -> Atom {
    trace!(
        &r.cof_nc,
        tco!(r, j1),
        tco!(r, j2),
        tco!(r, j3),
        tco!(r, j4),
        tco!(r, j5),
    ) * fco!(r, k1, k2, j1)
        * fco!(r, k2, k3, j2)
        * fco!(r, k3, k4, j3)
        * fco!(r, k4, k5, j4)
        * fco!(r, k5, k1, j5)
}

fn form_color_tloop_gg5(r: &TestReps) -> Atom {
    fco!(r, i1, i2, j1)
        * fco!(r, i2, i3, j2)
        * fco!(r, i3, i4, j3)
        * fco!(r, i4, i5, j4)
        * fco!(r, i5, i1, j5)
        * fco!(r, k1, k2, j1)
        * fco!(r, k2, k3, j2)
        * fco!(r, k3, k4, j3)
        * fco!(r, k4, k5, j4)
        * fco!(r, k5, k1, j5)
}

#[test]
#[ignore = "pending FORM color.frm size-5 Q10 invariant-family reduction"]
fn form_github_color_tloop_q10_size_5() {
    test_initialize();
    let r = TestReps::new();
    let expr = form_color_tloop_q10(&r);

    assert_snapshot!(expr.simplify_color().expand().to_bare_ordered_string(), @"trace(cof(Nc),cyclic(t(coad(NA,j1),in,out),t(coad(NA,j2),in,out),t(coad(NA,j3),in,out),t(coad(NA,j4),in,out),t(coad(NA,j5),in,out),t(coad(NA,j1),in,out),t(coad(NA,j2),in,out),t(coad(NA,j3),in,out),t(coad(NA,j4),in,out),t(coad(NA,j5),in,out)))");
}

#[test]
#[ignore = "pending FORM color.frm size-5 G10 invariant-family reduction"]
fn form_github_color_tloop_g10_size_5() {
    test_initialize();
    let r = TestReps::new();
    let expr = form_color_tloop_g10(&r);

    assert_snapshot!(expr.simplify_color().expand().to_bare_ordered_string(), @"f(coad(NA,i1),coad(NA,i10),coad(NA,j5))*f(coad(NA,i1),coad(NA,i2),coad(NA,j1))*f(coad(NA,i10),coad(NA,i9),coad(NA,j4))*f(coad(NA,i2),coad(NA,i3),coad(NA,j2))*f(coad(NA,i3),coad(NA,i4),coad(NA,j3))*f(coad(NA,i4),coad(NA,i5),coad(NA,j4))*f(coad(NA,i5),coad(NA,i6),coad(NA,j5))*f(coad(NA,i6),coad(NA,i7),coad(NA,j1))*f(coad(NA,i7),coad(NA,i8),coad(NA,j2))*f(coad(NA,i8),coad(NA,i9),coad(NA,j3))");
}

#[test]
#[ignore = "pending FORM color.frm size-5 QQ5 invariant-family reduction"]
fn form_github_color_tloop_qq5_size_5() {
    test_initialize();
    let r = TestReps::new();
    let expr = form_color_tloop_qq5(&r);

    assert_snapshot!(expr.simplify_color().expand().to_bare_ordered_string(), @"(trace(cof(Nc),cyclic(t(coad(NA,j1),in,out),t(coad(NA,j2),in,out),t(coad(NA,j3),in,out),t(coad(NA,j4),in,out),t(coad(NA,j5),in,out))))^2");
}

#[test]
#[ignore = "pending FORM color.frm size-5 QG5 invariant-family reduction"]
fn form_github_color_tloop_qg5_size_5() {
    test_initialize();
    let r = TestReps::new();
    let expr = form_color_tloop_qg5(&r);

    assert_snapshot!(expr.simplify_color().expand().to_bare_ordered_string(), @"-1*f(coad(NA,j1),coad(NA,k1),coad(NA,k2))*f(coad(NA,j2),coad(NA,k2),coad(NA,k3))*f(coad(NA,j3),coad(NA,k3),coad(NA,k4))*f(coad(NA,j4),coad(NA,k4),coad(NA,k5))*f(coad(NA,j5),coad(NA,k1),coad(NA,k5))*trace(cof(Nc),cyclic(t(coad(NA,j1),in,out),t(coad(NA,j2),in,out),t(coad(NA,j3),in,out),t(coad(NA,j4),in,out),t(coad(NA,j5),in,out)))");
}

#[test]
#[ignore = "pending FORM color.frm size-5 GG5 invariant-family reduction"]
fn form_github_color_tloop_gg5_size_5() {
    test_initialize();
    let r = TestReps::new();
    let expr = form_color_tloop_gg5(&r);

    assert_snapshot!(expr.simplify_color().expand().to_bare_ordered_string(), @"f(coad(NA,i1),coad(NA,i2),coad(NA,j1))*f(coad(NA,i1),coad(NA,i5),coad(NA,j5))*f(coad(NA,i2),coad(NA,i3),coad(NA,j2))*f(coad(NA,i3),coad(NA,i4),coad(NA,j3))*f(coad(NA,i4),coad(NA,i5),coad(NA,j4))*f(coad(NA,j1),coad(NA,k1),coad(NA,k2))*f(coad(NA,j2),coad(NA,k2),coad(NA,k3))*f(coad(NA,j3),coad(NA,k3),coad(NA,k4))*f(coad(NA,j4),coad(NA,k4),coad(NA,k5))*f(coad(NA,j5),coad(NA,k1),coad(NA,k5))");
}

#[test]
#[ignore = "pending color.h simpli strategy over invariant environments"]
fn simpli_contracts_projected_f_pair() {
    test_initialize();
    let expr = parse_lit!(
        f(coad(NA, x), coad(NA, a), coad(NA, c))
            * f(coad(NA, x), coad(NA, b), coad(NA, c))
            * invariant_environment(a, b),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"CA*g(coad(NA,a),coad(NA,b))*invariant_environment(a,b)");
}

#[test]
#[ignore = "pending color.tar.gz tloop qloop size 3 family"]
fn tloop_qloop_size_3() {
    test_initialize();
    let result =
        parse_lit!(NA * TR * CF ^ 2 - 3 / 2 * NA * TR * CA * CF + 1 / 2 * NA * TR * CA ^ 2);

    assert_snapshot!(result.to_bare_ordered_string(), @"-3/2*CA*CF*NA*TR+1/2*CA^2*NA*TR+CF^2*NA*TR");
}

#[test]
#[ignore = "pending color.tar.gz tloop gloop size 3 family"]
fn tloop_gloop_size_3() {
    test_initialize();
    let result = parse_lit!(0);

    assert_snapshot!(result.to_bare_ordered_string(), @"0");
}

#[test]
#[ignore = "pending color.tar.gz tloop qqloop size 3 family"]
fn tloop_qqloop_size_3() {
    test_initialize();
    let result = parse_lit!(-1 / 4 * NA * TR ^ 2 * CA + d33(R1, R2));

    assert_snapshot!(result.to_bare_ordered_string(), @"-1/4*CA*NA*TR^2+d33(R1,R2)");
}

#[test]
#[ignore = "pending color.tar.gz tloop qgloop size 3 family"]
fn tloop_qgloop_size_3() {
    test_initialize();
    let result = parse_lit!(1𝑖 / 4 * NA * TR * CA ^ 2);

    assert_snapshot!(result.to_bare_ordered_string(), @"1𝑖/4*CA^2*NA*TR");
}

#[test]
#[ignore = "pending color.tar.gz tloop ggloop size 3 family"]
fn tloop_ggloop_size_3() {
    test_initialize();
    let result = parse_lit!(1 / 4 * NA * CA ^ 3);

    assert_snapshot!(result.to_bare_ordered_string(), @"1/4*CA^3*NA");
}

#[test]
#[ignore = "pending color.tar.gz fixed g14 family"]
fn tloop_g14() {
    test_initialize();
    let result = parse_lit!(
        1 / 648 * NA * CA ^ 7 - 8 / 15 * d444(A1, A2, A3) * CA + 16 / 9 * d644(A1, A2, A3)
    );

    assert_snapshot!(result.to_bare_ordered_string(), @"-8/15*CA*d444(A1,A2,A3)+1/648*CA^7*NA+16/9*d644(A1,A2,A3)");
}

#[test]
#[ignore = "pending color.tar.gz fixed fiveq family"]
fn tloop_fiveq() {
    test_initialize();
    let result = parse_lit!(
        1 / 192 * NA * TR
            ^ 5 * CA
            ^ 3 + 1 / 4 * d33(R1, R2) * TR
            ^ 3 * CA
            ^ 2 + 5 / 48 * d33(R1, R2) * d33(R3, R4) * NA
            ^ -1 * TR * CA + 5 / 48 * d33(R1, R3) * d33(R2, R4) * NA
            ^ -1 * TR * CA + 1 / 8 * d33(R1, R4) * d33(R2, R3) * NA
            ^ -1 * TR * CA + 1 / 16 * d44(R1, A1) * TR
            ^ 4 + 3 / 8 * d433(R3, R1, R2) * TR
            ^ 2 * CA + 1 / 2 * d3333(R1, R2, R3, R4) * TR * CA + d43333a(R5, R2, R1, R4, R3)
    );

    assert_snapshot!(result.to_bare_ordered_string(), @"1/16*TR^4*d44(R1,A1)+1/192*CA^3*NA*TR^5+1/2*CA*TR*d3333(R1,R2,R3,R4)+1/4*CA^2*TR^3*d33(R1,R2)+1/8*CA*NA^(-1)*TR*d33(R1,R4)*d33(R2,R3)+3/8*CA*TR^2*d433(R3,R1,R2)+5/48*CA*NA^(-1)*TR*d33(R1,R2)*d33(R3,R4)+5/48*CA*NA^(-1)*TR*d33(R1,R3)*d33(R2,R4)+d43333a(R5,R2,R1,R4,R3)");
}

#[test]
#[ignore = "pending color.tar.gz su.frm F3F3 case"]
fn su_f3f3() {
    test_initialize();
    let result = parse_lit!(2 * a ^ 3 * NF ^ -1 - 5 / 2 * a ^ 3 * NF + 1 / 2 * a ^ 3 * NF ^ 3);

    assert_snapshot!(result.to_bare_ordered_string(), @"-5/2*NF*a^3+1/2*NF^3*a^3+2*NF^(-1)*a^3");
}

#[test]
#[ignore = "pending color.tar.gz su.frm F4F4 case"]
fn su_f4f4() {
    test_initialize();
    let result = parse_lit!(
        -3 * a
            ^ 4 * nf
            ^ 2 * NF
            ^ -2 + 4 * a
            ^ 4 * nf
            ^ 2 - 7 / 6 * a
            ^ 4 * nf
            ^ 2 * NF
            ^ 2 + 1 / 6 * a
            ^ 4 * nf
            ^ 2 * NF
            ^ 4
    );

    assert_snapshot!(result.to_bare_ordered_string(), @"-3*NF^(-2)*a^4*nf^2+-7/6*NF^2*a^4*nf^2+1/6*NF^4*a^4*nf^2+4*a^4*nf^2");
}

#[test]
#[ignore = "pending color.tar.gz su.frm F4A4 case"]
fn su_f4a4() {
    test_initialize();
    let result = parse_lit!(
        -2 * a ^ 4 * nf * NF + 5 / 3 * a ^ 4 * nf * NF ^ 3 + 1 / 3 * a ^ 4 * nf * NF ^ 5
    );

    assert_snapshot!(result.to_bare_ordered_string(), @"-2*NF*a^4*nf+1/3*NF^5*a^4*nf+5/3*NF^3*a^4*nf");
}

#[test]
#[ignore = "pending color.tar.gz su.frm A4A4 case"]
fn su_a4a4() {
    test_initialize();
    let result =
        parse_lit!(-24 * a ^ 4 * NF ^ 2 + 70 / 3 * a ^ 4 * NF ^ 4 + 2 / 3 * a ^ 4 * NF ^ 6);

    assert_snapshot!(result.to_bare_ordered_string(), @"-24*NF^2*a^4+2/3*NF^6*a^4+70/3*NF^4*a^4");
}

#[test]
#[ignore = "pending color.tar.gz su.frm FnFn n=3 case"]
fn su_fnfn_n3() {
    test_initialize();
    let result = parse_lit!(2 * a ^ 3 * nf ^ 2 * NF ^ -1 - 2 * a ^ 3 * nf ^ 2 * NF);

    assert_snapshot!(result.to_bare_ordered_string(), @"-2*NF*a^3*nf^2+2*NF^(-1)*a^3*nf^2");
}

#[test]
#[ignore = "pending color.tar.gz su.frm FnAn n=3 case"]
fn su_fnan_n3() {
    test_initialize();
    let result = parse_lit!(-1𝑖 * a ^ 3 * nf * NF ^ 2 + 1𝑖 * a ^ 3 * nf * NF ^ 4);

    assert_snapshot!(result.to_bare_ordered_string(), @"-1𝑖*NF^2*a^3*nf+1𝑖*NF^4*a^3*nf");
}

#[test]
#[ignore = "pending color.tar.gz su.frm AnAn n=3 case"]
fn su_anan_n3() {
    test_initialize();
    let result = parse_lit!(-2 * a ^ 3 * NF ^ 3 + 2 * a ^ 3 * NF ^ 5);

    assert_snapshot!(result.to_bare_ordered_string(), @"-2*NF^3*a^3+2*NF^5*a^3");
}

#[test]
#[ignore = "pending color.tar.gz su.frm qloop n=3 case"]
fn su_qloop_n3() {
    test_initialize();
    let result = parse_lit!(-a ^ 3 * nf * NF ^ -2 + a ^ 3 * nf * NF ^ 2);

    assert_snapshot!(result.to_bare_ordered_string(), @"-1*NF^(-2)*a^3*nf+NF^2*a^3*nf");
}

#[test]
#[ignore = "pending color.tar.gz su.frm gloop n=3 case"]
fn su_gloop_n3() {
    test_initialize();
    let result = parse_lit!(0);

    assert_snapshot!(result.to_bare_ordered_string(), @"0");
}
