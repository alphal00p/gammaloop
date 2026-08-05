use super::*;
use crate::coad;
use insta::assert_snapshot;
use spenso::{antisym, g, mink, shadowing::IntoAtom};
use symbolica::{function, symbol};

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
fn metric_contraction_respects_antisymmetric_normalization() {
    test_initialize();
    let coad8 = ColorAdjoint {}.new_rep(8);
    // Seed the symbol order from the three-loop graph that exposed the lost
    // permutation parity in the old sequence-wildcard contraction.
    let a = slot!(coad8, standalone_metric_a);
    let x = slot!(coad8, standalone_metric_x);
    let u = slot!(coad8, standalone_metric_u);
    let v = slot!(coad8, standalone_metric_v);
    let _ = (color_f!(a, x, u) * color_f!(a, x, v)).simplify_metrics();
    let r = slot!(coad8, standalone_metric_r);
    let s = slot!(coad8, standalone_metric_s);
    let j = slot!(coad8, standalone_metric_j);

    let contracted = (g!(r, s) * color_f!(u, j, r)).simplify_metrics();
    // Re-normalizing f after direct slot substitution may change its displayed
    // argument order and coefficient, so compare with a freshly built target.
    assert_eq!(
        contracted.to_bare_ordered_string(),
        color_f!(u, j, s).to_bare_ordered_string(),
    );

    let closed = g!(u, v) * g!(r, s) * color_f!(u, j, r) * color_f!(v, j, s);
    assert_eq!(
        closed.simplify_color(),
        Atom::num(8) * color_cas!(2, &coad8)
    );
}

#[test]
fn metric_contraction_only_replaces_the_immediate_slot() {
    test_initialize();
    let coad8 = ColorAdjoint {}.new_rep(8);
    let u = slot!(coad8, standalone_nested_metric_u);
    let j = slot!(coad8, standalone_nested_metric_j);
    let r = slot!(coad8, standalone_nested_metric_r);
    let s = slot!(coad8, standalone_nested_metric_s);
    let probe = symbol!("spenso::standalone_metric_probe");
    let structure = color_f!(u, j, r);
    let r_atom = r.into_atom();
    let s_atom = s.into_atom();

    let contracted = (g!(r, s) * function!(probe, &structure, r_atom)).simplify_metrics();
    assert_eq!(contracted, function!(probe, structure, s_atom));
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
    let coad8 = ColorAdjoint {}.new_rep(8);
    // Seed the symbol order that made the unoriented triangle reduction pick
    // the wrong sign in the three-loop graph.
    let r = slot!(coad8, standalone_triangle_r);
    let s = slot!(coad8, standalone_triangle_s);
    let x = slot!(coad8, standalone_triangle_x);
    let z = slot!(coad8, standalone_triangle_z);
    let _ = g!(r, s) * color_f!(x, z, r);
    let u = slot!(coad8, standalone_triangle_u);
    let y = slot!(coad8, standalone_triangle_y);
    let triangle = color_f!(u, x, r) * color_f!(y, z, r) * color_f!(u, z, s);
    let closing_structure = color_f!(s, x, y);
    let expected = Atom::num(-1) * color_cas!(2, &coad8) / Atom::num(2) * color_f!(s, x, y);

    assert_eq!(triangle.simplify_color(), expected);
    // Closing the loop before or after triangle reduction must keep the same
    // antisymmetric orientation.
    assert_eq!(
        (&triangle * &closing_structure).simplify_color(),
        (triangle.simplify_color() * closing_structure).simplify_color(),
    );
}

#[test]
fn six_f_k33_odd_automorphism_simplifies_to_zero() {
    test_initialize();
    let r = TestReps::new();
    let contraction = fco!(r, a, b, c)
        * fco!(r, d, e, f)
        * fco!(r, h, i, j)
        * fco!(r, a, d, h)
        * fco!(r, b, e, i)
        * fco!(r, c, f, j);
    // Exchanging the first two K3,3 vertices is a dummy relabeling, but it
    // reverses three f tensors. Thus C=-C and the contraction is zero.
    let odd_relabeling = fco!(r, d, e, f)
        * fco!(r, a, b, c)
        * fco!(r, h, i, j)
        * fco!(r, d, a, h)
        * fco!(r, e, b, i)
        * fco!(r, f, c, j);

    assert_eq!(odd_relabeling, -&contraction);
    // Color simplification applies signed tensor canonicalization only to its
    // extracted color factor.
    assert!(contraction.simplify_color().is_zero());
}

#[test]
fn six_f_k33_with_structured_indices_simplifies_to_zero() {
    test_initialize();
    let hedge = symbol!("spenso::hedge");
    let structured_slot = |index| coad!(8, function!(hedge, Atom::num(index as i64)));
    let a = structured_slot(1);
    let b = structured_slot(2);
    let c = structured_slot(3);
    let d = structured_slot(4);
    let e = structured_slot(5);
    let f = structured_slot(6);
    let h = structured_slot(7);
    let i = structured_slot(8);
    let j = structured_slot(9);
    let contraction = color_f!(&a, &b, &c)
        * color_f!(&d, &e, &f)
        * color_f!(&h, &i, &j)
        * color_f!(&a, &d, &h)
        * color_f!(&b, &e, &i)
        * color_f!(&c, &f, &j);

    assert!(contraction.simplify_color().is_zero());
}

#[test]
fn rqft_ghost_three_loop_d2_six_f_terms_vanish_independently() {
    test_initialize();
    // These are the two six-f contractions emitted by the d2 RQFT ghost
    // three-loop fixture, with their original structured index payloads.
    let contractions = [
        parse_lit!(
            f(coad(8, hedge(1)), coad(8, hedge(10)), coad(8, hedge(7)))
                * f(coad(8, hedge(1)), coad(8, hedge(3)), coad(8, hedge(8)))
                * f(
                    coad(8, hedge(10)),
                    coad(8, hedge(12)),
                    coad(8, vertex(4, 1))
                )
                * f(coad(8, hedge(12)), coad(8, hedge(3)), coad(8, hedge(5)))
                * f(coad(8, hedge(14)), coad(8, hedge(5)), coad(8, hedge(7)))
                * f(coad(8, hedge(14)), coad(8, hedge(8)), coad(8, vertex(4, 1))),
            default_namespace = "spenso"
        ),
        parse_lit!(
            f(coad(8, hedge(1)), coad(8, hedge(10)), coad(8, hedge(6)))
                * f(coad(8, hedge(1)), coad(8, hedge(3)), coad(8, hedge(8)))
                * f(
                    coad(8, hedge(10)),
                    coad(8, hedge(12)),
                    coad(8, vertex(4, 1))
                )
                * f(coad(8, hedge(12)), coad(8, hedge(3)), coad(8, hedge(5)))
                * f(coad(8, hedge(14)), coad(8, hedge(5)), coad(8, hedge(6)))
                * f(coad(8, hedge(14)), coad(8, hedge(8)), coad(8, vertex(4, 1))),
            default_namespace = "spenso"
        ),
    ];

    // They carried coefficients -3/16 and +3/16 in d2, but their vanishing
    // does not depend on cancellation between the two summands.
    for contraction in contractions {
        assert!(
            contraction
                .canonize::<AbstractIndex>(AbstractIndex::Dummy)
                .is_zero()
        );
    }
}

#[test]
fn signed_pruning_reenters_color_simplification() {
    test_initialize();
    let r = TestReps::new();
    let odd = fco!(r, a, b, c)
        * fco!(r, d, e, f)
        * fco!(r, h, i, j)
        * fco!(r, a, d, h)
        * fco!(r, b, e, i)
        * fco!(r, c, f, j);
    let vanishing_factor = fco!(r, u, v, w);
    let surviving_factor = fco!(r, x, y, z);
    let expected = surviving_factor.clone().pow(2).simplify_color();

    assert_eq!(
        ((odd * vanishing_factor + &surviving_factor) * surviving_factor).simplify_color(),
        expected
    );
}

#[test]
fn color_simplification_does_not_canonicalize_lorentz_tensors() {
    test_initialize();
    let r = TestReps::new();
    let a = mink!(4, color_only_canonicalization_a);
    let b = mink!(4, color_only_canonicalization_b);
    let c = mink!(4, color_only_canonicalization_c);
    let odd_lorentz =
        antisym!(a.clone(), b.clone()) * antisym!(b.clone(), c.clone()) * antisym!(c, a);
    let color = fco!(r, x, y, z);
    let expected = color.simplify_color() * &odd_lorentz;

    assert!(!expected.is_zero());
    assert_eq!((color * &odd_lorentz).simplify_color(), expected);
    assert!(
        odd_lorentz
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .is_zero()
    );
}

#[test]
fn two_f_even_automorphism_survives_canonicalization() {
    test_initialize();
    let r = TestReps::new();
    let contraction = fco!(r, a, b, c) * fco!(r, a, b, c);

    // Exchanging two contracted edges reverses both f tensors, so the total
    // graph-automorphism sign is positive.
    assert!(
        !contraction
            .canonize::<AbstractIndex>(AbstractIndex::Dummy)
            .is_zero()
    );
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

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"-1/6*f(coad(NA,a),coad(NA,c),coad(NA,d_0))*f(coad(NA,b),coad(NA,d),coad(NA,d_0))*idx(2,cof(Nc))+1/3*f(coad(NA,a),coad(NA,d),coad(NA,d_0))*f(coad(NA,b),coad(NA,c),coad(NA,d_0))*idx(2,cof(Nc))+1𝑖/2*f(coad(NA,a),coad(NA,b),coad(NA,x))*trace(cof(Nc),sym(t(coad(NA,c),in,out),t(coad(NA,d),in,out),t(coad(NA,x),in,out)))+1𝑖/2*f(coad(NA,c),coad(NA,d),coad(NA,x))*trace(cof(Nc),sym(t(coad(NA,a),in,out),t(coad(NA,b),in,out),t(coad(NA,x),in,out)))+trace(cof(Nc),sym(t(coad(NA,a),in,out),t(coad(NA,b),in,out),t(coad(NA,c),in,out),t(coad(NA,d),in,out)))");
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
