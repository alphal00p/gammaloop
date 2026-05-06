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
    initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, i))),
        default_namespace = "spenso"
    );

    assert_color_zero(expr);
}

#[test]
fn two_generator_trace_normalizes_to_tr_metric() {
    initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, i))),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"TR*g(coad(-1+Nc^2,a),coad(-1+Nc^2,b))");
}

#[test]
fn fierz_generator_contraction() {
    initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().expand().simplify_metrics().to_bare_ordered_string(), @"-1*Nc^(-1)*TR*g(cof(Nc,i),dind(cof(Nc,j)))*g(cof(Nc,k),dind(cof(Nc,l)))+TR*g(cof(Nc,i),dind(cof(Nc,l)))*g(cof(Nc,k),dind(cof(Nc,j)))");
}

#[test]
fn separated_generator_casimir_shortcut() {
    initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, k)))
            * t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().simplify_metrics().to_bare_ordered_string(), @"-1*Nc^(-1)*TR*t(coad(-1+Nc^2,b),cof(Nc,i),dind(cof(Nc,l)))");
}

#[test]
fn chain_one_generator_trace_vanishes() {
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
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"TR*g(coad(NA,a),coad(NA,b))");
}

#[test]
fn three_generator_trace_terminal() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, c)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1𝑖/2*TR*f(coad(NA,a),coad(NA,b),coad(NA,c))+trace(cof(Nc),sym(t(coad(NA,a),in,out),t(coad(NA,b),in,out),t(coad(NA,c),in,out)))");
}

#[test]
fn adjacent_generator_casimir_chain() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), k),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, a)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"CF*g(cof(Nc,i),dind(cof(Nc,k)))");
}

#[test]
fn separated_generator_casimir_trace() {
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, c)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"-1/2*CA*TR*g(coad(NA,b),coad(NA,c))+CF*TR*g(coad(NA,b),coad(NA,c))");
}

#[test]
fn two_f_loop_contracts_to_ca_metric() {
    initialize();
    let expr = parse_lit!(
        f(coad(NA, a), coad(NA, c), coad(NA, d)) * f(coad(NA, b), coad(NA, c), coad(NA, d)),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"CA*g(coad(NA,a),coad(NA,b))");
}

#[test]
fn three_f_loop_contracts_to_ca_f() {
    initialize();
    let expr = parse_lit!(
        f(coad(NA, a), coad(NA, b), coad(NA, e))
            * f(coad(NA, b), coad(NA, c), coad(NA, f_))
            * f(coad(NA, c), coad(NA, a), coad(NA, g_)),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1/2*CA*f(coad(NA,e),coad(NA,f_),coad(NA,g_))");
}

#[test]
fn mixed_trace_structure_contraction() {
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

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1𝑖/2*CA*TR*g(coad(NA,c),coad(NA,d))");
}

#[test]
fn symmetric_invariant_d33_partial_contraction() {
    initialize();
    let expr = parse_lit!(
        d(sym_x, coad(NA, a), coad(NA, b), coad(NA, c))
            * d(sym_y, coad(NA, a), coad(NA, b), coad(NA, d_)),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"NA^(-1)*d33(sym_x,sym_y)*g(coad(NA,c),coad(NA,d_))");
}

#[test]
fn four_generator_trace_terminal() {
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, c)),
        color_t!(slot!(r.coad_na, d)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"-1/6*TR*f(coad(NA,a),coad(NA,c),coad(NA,x))*f(coad(NA,b),coad(NA,d),coad(NA,x))+1/3*TR*f(coad(NA,a),coad(NA,d),coad(NA,x))*f(coad(NA,b),coad(NA,c),coad(NA,x))+1𝑖/2*f(coad(NA,a),coad(NA,b),coad(NA,x))*trace(cof(Nc),sym(t(coad(NA,c),in,out),t(coad(NA,d),in,out),t(coad(NA,x),in,out)))+1𝑖/2*f(coad(NA,c),coad(NA,d),coad(NA,x))*trace(cof(Nc),sym(t(coad(NA,a),in,out),t(coad(NA,b),in,out),t(coad(NA,x),in,out)))+trace(cof(Nc),sym(t(coad(NA,a),in,out),t(coad(NA,b),in,out),t(coad(NA,c),in,out),t(coad(NA,d),in,out)))");
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
    let r = TestReps::new();
    let expr = form_color_tloop_q10(&r);

    assert_snapshot!(expr.simplify_color().expand().to_bare_ordered_string(), @"125/72*CA^4*NA*TR+-155/24*CA^3*CF*NA*TR+35/4*CA^2*CF^2*NA*TR+-5*CA*CF^3*NA*TR+-6*CA*d44(cOlpR1,cOlpA1)+1/3*TR*d44(cOlpA1,cOlpA2)+CF^4*NA*TR+5*CF*d44(cOlpR1,cOlpA1)");
}

#[test]
#[ignore = "pending FORM color.frm size-5 G10 invariant-family reduction"]
fn form_github_color_tloop_g10_size_5() {
    let r = TestReps::new();
    let expr = form_color_tloop_g10(&r);

    assert_snapshot!(expr.simplify_color().expand().to_bare_ordered_string(), @"-1/36*CA^5*NA+2/3*CA*d44(cOlpA1,cOlpA2)");
}

#[test]
#[ignore = "pending FORM color.frm size-5 QQ5 invariant-family reduction"]
fn form_github_color_tloop_qq5_size_5() {
    let r = TestReps::new();
    let expr = form_color_tloop_qq5(&r);

    assert_snapshot!(expr.simplify_color().expand().to_bare_ordered_string(), @"-5/144*CA^3*NA*TR^2+11/48*CA^2*d33(cOlpR1,cOlpR2)+1/6*TR*d44(cOlpR1,cOlpA1)+-5/6*CA*d44(cOlpR1,cOlpR2)+d55(cOlpR1,cOlpR2)");
}

#[test]
#[ignore = "pending FORM color.frm size-5 QG5 invariant-family reduction"]
fn form_github_color_tloop_qg5_size_5() {
    let r = TestReps::new();
    let expr = form_color_tloop_qg5(&r);

    assert_snapshot!(expr.simplify_color().expand().to_bare_ordered_string(), @"-5𝑖/144*CA^4*NA*TR+-3𝑖/4*CA*d44(cOlpR1,cOlpA1)+1𝑖/12*TR*d44(cOlpA1,cOlpA2)");
}

#[test]
#[ignore = "pending FORM color.frm size-5 GG5 invariant-family reduction"]
fn form_github_color_tloop_gg5_size_5() {
    let r = TestReps::new();
    let expr = form_color_tloop_gg5(&r);

    assert_snapshot!(expr.simplify_color().expand().to_bare_ordered_string(), @"5/144*CA^5*NA+2/3*CA*d44(cOlpA1,cOlpA2)");
}

#[test]
#[ignore = "pending color.h simpli strategy over invariant environments"]
fn simpli_contracts_projected_f_pair() {
    initialize();
    let expr = parse_lit!(
        f(coad(NA, x), coad(NA, a), coad(NA, c))
            * f(coad(NA, x), coad(NA, b), coad(NA, c))
            * invariant_environment(a, b),
        default_namespace = "spenso"
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1/2*CA*projected_invariant(a,b,invariant_environment(a,b))");
}

#[test]
#[ignore = "pending color.tar.gz tloop qloop size 3 family"]
fn tloop_qloop_size_3() {
    initialize();
    let result =
        parse_lit!(NA * TR * CF ^ 2 - 3 / 2 * NA * TR * CA * CF + 1 / 2 * NA * TR * CA ^ 2);

    assert_snapshot!(result.to_bare_ordered_string(), @"1/2*CA^2*NA*TR+-3/2*CA*CF*NA*TR+CF^2*NA*TR");
}

#[test]
#[ignore = "pending color.tar.gz tloop gloop size 3 family"]
fn tloop_gloop_size_3() {
    initialize();
    let result = parse_lit!(0);

    assert_snapshot!(result.to_bare_ordered_string(), @"0");
}

#[test]
#[ignore = "pending color.tar.gz tloop qqloop size 3 family"]
fn tloop_qqloop_size_3() {
    initialize();
    let result = parse_lit!(-1 / 4 * NA * TR ^ 2 * CA + d33(R1, R2));

    assert_snapshot!(result.to_bare_ordered_string(), @"-1/4*CA*NA*TR^2+d33(R1,R2)");
}

#[test]
#[ignore = "pending color.tar.gz tloop qgloop size 3 family"]
fn tloop_qgloop_size_3() {
    initialize();
    let result = parse_lit!(1𝑖 / 4 * NA * TR * CA ^ 2);

    assert_snapshot!(result.to_bare_ordered_string(), @"1/4*1𝑖*CA^2*NA*TR");
}

#[test]
#[ignore = "pending color.tar.gz tloop ggloop size 3 family"]
fn tloop_ggloop_size_3() {
    initialize();
    let result = parse_lit!(1 / 4 * NA * CA ^ 3);

    assert_snapshot!(result.to_bare_ordered_string(), @"1/4*CA^3*NA");
}

#[test]
#[ignore = "pending color.tar.gz fixed g14 family"]
fn tloop_g14() {
    initialize();
    let result = parse_lit!(
        1 / 648 * NA * CA ^ 7 - 8 / 15 * d444(A1, A2, A3) * CA + 16 / 9 * d644(A1, A2, A3)
    );

    assert_snapshot!(result.to_bare_ordered_string(), @"1/648*CA^7*NA+-8/15*CA*d444(A1,A2,A3)+16/9*d644(A1,A2,A3)");
}

#[test]
#[ignore = "pending color.tar.gz fixed fiveq family"]
fn tloop_fiveq() {
    initialize();
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

    assert_snapshot!(result.to_bare_ordered_string(), @"1/192*CA^3*NA*TR^5+1/4*CA^2*TR^3*d33(R1,R2)+5/48*CA*NA^(-1)*TR*d33(R1,R2)*d33(R3,R4)+5/48*CA*NA^(-1)*TR*d33(R1,R3)*d33(R2,R4)+1/8*CA*NA^(-1)*TR*d33(R1,R4)*d33(R2,R3)+3/8*CA*TR^2*d433(R3,R1,R2)+1/2*CA*TR*d3333(R1,R2,R3,R4)+1/16*TR^4*d44(R1,A1)+d43333a(R5,R2,R1,R4,R3)");
}

#[test]
#[ignore = "pending color.tar.gz su.frm F3F3 case"]
fn su_f3f3() {
    initialize();
    let result = parse_lit!(2 * a ^ 3 * NF ^ -1 - 5 / 2 * a ^ 3 * NF + 1 / 2 * a ^ 3 * NF ^ 3);

    assert_snapshot!(result.to_bare_ordered_string(), @"2*NF^(-1)*a^3+-5/2*NF*a^3+1/2*NF^3*a^3");
}

#[test]
#[ignore = "pending color.tar.gz su.frm F4F4 case"]
fn su_f4f4() {
    initialize();
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

    assert_snapshot!(result.to_bare_ordered_string(), @"-3*NF^(-2)*a^4*nf^2+4*a^4*nf^2+-7/6*NF^2*a^4*nf^2+1/6*NF^4*a^4*nf^2");
}

#[test]
#[ignore = "pending color.tar.gz su.frm F4A4 case"]
fn su_f4a4() {
    initialize();
    let result = parse_lit!(
        -2 * a ^ 4 * nf * NF + 5 / 3 * a ^ 4 * nf * NF ^ 3 + 1 / 3 * a ^ 4 * nf * NF ^ 5
    );

    assert_snapshot!(result.to_bare_ordered_string(), @"-2*NF*a^4*nf+5/3*NF^3*a^4*nf+1/3*NF^5*a^4*nf");
}

#[test]
#[ignore = "pending color.tar.gz su.frm A4A4 case"]
fn su_a4a4() {
    initialize();
    let result =
        parse_lit!(-24 * a ^ 4 * NF ^ 2 + 70 / 3 * a ^ 4 * NF ^ 4 + 2 / 3 * a ^ 4 * NF ^ 6);

    assert_snapshot!(result.to_bare_ordered_string(), @"-24*NF^2*a^4+70/3*NF^4*a^4+2/3*NF^6*a^4");
}

#[test]
#[ignore = "pending color.tar.gz su.frm FnFn n=3 case"]
fn su_fnfn_n3() {
    initialize();
    let result = parse_lit!(2 * a ^ 3 * nf ^ 2 * NF ^ -1 - 2 * a ^ 3 * nf ^ 2 * NF);

    assert_snapshot!(result.to_bare_ordered_string(), @"2*NF^(-1)*a^3*nf^2+-2*NF*a^3*nf^2");
}

#[test]
#[ignore = "pending color.tar.gz su.frm FnAn n=3 case"]
fn su_fnan_n3() {
    initialize();
    let result = parse_lit!(-1𝑖 * a ^ 3 * nf * NF ^ 2 + 1𝑖 * a ^ 3 * nf * NF ^ 4);

    assert_snapshot!(result.to_bare_ordered_string(), @"-1𝑖*NF^2*a^3*nf+1𝑖*NF^4*a^3*nf");
}

#[test]
#[ignore = "pending color.tar.gz su.frm AnAn n=3 case"]
fn su_anan_n3() {
    initialize();
    let result = parse_lit!(-2 * a ^ 3 * NF ^ 3 + 2 * a ^ 3 * NF ^ 5);

    assert_snapshot!(result.to_bare_ordered_string(), @"-2*NF^3*a^3+2*NF^5*a^3");
}

#[test]
#[ignore = "pending color.tar.gz su.frm qloop n=3 case"]
fn su_qloop_n3() {
    initialize();
    let result = parse_lit!(-a ^ 3 * nf * NF ^ -2 + a ^ 3 * nf * NF ^ 2);

    assert_snapshot!(result.to_bare_ordered_string(), @"-1*NF^(-2)*a^3*nf+NF^2*a^3*nf");
}

#[test]
#[ignore = "pending color.tar.gz su.frm gloop n=3 case"]
fn su_gloop_n3() {
    initialize();
    let result = parse_lit!(0);

    assert_snapshot!(result.to_bare_ordered_string(), @"0");
}
