use insta::assert_snapshot;
use spenso::shadowing::symbolica_utils::AtomCoreExt;
use spenso::structure::IndexlessNamedStructure;
use spenso::structure::PermutedStructure;

static _CF: LazyLock<PermutedStructure<IndexlessNamedStructure<Symbol, ()>>> =
    LazyLock::new(|| {
        IndexlessNamedStructure::from_iter(
            [
                ColorAdjoint {}.new_rep(8),
                ColorAdjoint {}.new_rep(2),
                ColorAdjoint {}.new_rep(2),
                ColorAdjoint {}.new_rep(4),
                ColorAdjoint {}.new_rep(2),
            ],
            CS.f,
            None,
        )
    });

static CT: LazyLock<PermutedStructure<IndexlessNamedStructure<Symbol, ()>>> = LazyLock::new(|| {
    IndexlessNamedStructure::from_iter(
        [
            ColorAntiFundamental {}.new_rep(3).to_lib(),
            ColorFundamental {}.new_rep(3).to_lib(),
            ColorAdjoint {}.new_rep(8).to_lib(),
        ],
        CS.t,
        None,
    )
});
use spenso::{chain, network::parsing::ShadowedStructure, slot, structure::permuted::Perm, trace};
use symbolica::{parse, parse_lit};

use crate::dirac::PS;
use crate::tensor::SymbolicTensor;
use crate::test::test_initialize;
use crate::{
    IndexTooling, dirac::GammaSimplifier, metric::MetricSimplifier, representations::initialize,
};

use super::*;

use crate::test_support::{TestReps, assert_bare_snapshot};

fn assert_color_zero(expr: Atom) {
    assert!(expr.simplify_color().is_zero());
}

#[test]
fn test_color_structures() {
    let f = IndexlessNamedStructure::<Symbol, ()>::from_iter(
        [
            ColorAdjoint {}.new_rep(8),
            ColorAdjoint {}.new_rep(2),
            ColorAdjoint {}.new_rep(2),
            ColorAdjoint {}.new_rep(4),
            ColorAdjoint {}.new_rep(2),
            ColorAdjoint {}.new_rep(7),
        ],
        symbol!("test"),
        None,
    )
    .clone()
    .reindex([5, 4, 2, 3, 1, 0])
    .unwrap()
    .map_structure(|a| SymbolicTensor::from_named(&a).unwrap());

    let f_s = f.structure.structure.clone();

    let f_p = f.clone().permute_inds();

    let f_parsed = PermutedStructure::<ShadowedStructure<AbstractIndex>>::try_from(
        &f_p.expression.simplify_metrics(),
    )
    .unwrap();

    assert_eq!(f.index_permutation, f_parsed.index_permutation);
    assert!(f_parsed.rep_permutation.is_identity());

    let f_p = f.clone().permute_reps_wrapped().permute_inds();

    let f_parsed = PermutedStructure::<ShadowedStructure<AbstractIndex>>::try_from(
        &f_p.expression.simplify_metrics(),
    )
    .unwrap();

    assert_eq!(f.index_permutation, f_parsed.index_permutation);
    assert_eq!(f.rep_permutation, f_parsed.rep_permutation);

    println!("Parsed: {}", f_parsed.index_permutation);
    println!("OG: {}", f.index_permutation);

    println!(
        "Structure:{}\nPermuted:{}\nPermuted Structure{}\nMetric simplified{}",
        f_s,
        f_p,
        f_p.structure,
        f_p.expression.simplify_metrics()
    );

    let t = CT
        .clone()
        .reindex([4, 2, 3])
        .unwrap()
        .map_structure(|a| SymbolicTensor::from_named(&a).unwrap())
        .permute_inds();

    println!("{t}")
}

#[test]
fn test_color_simplification() {
    initialize();
    let atom = parse_lit!(
        f(
            coad(Nc ^ 2 - 1, 1),
            coad(Nc ^ 2 - 1, 2),
            coad(Nc ^ 2 - 1, 3)
        ) ^ 2,
        default_namespace = "spenso"
    );
    let simplified = atom.simplify_color();

    assert_snapshot!(simplified.to_bare_ordered_string(), @"-1*ca+Nc^2*ca");
}

#[test]
fn form_one_generator_trace_vanishes() {
    initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, i))),
        default_namespace = "spenso"
    );

    assert_color_zero(expr);
}

#[test]
fn form_two_generator_trace_normalizes_to_tr_metric() {
    initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, i))),
        default_namespace = "spenso"
    );

    assert_bare_snapshot!(expr.simplify_color(),@"TR*g(coad(-1+Nc^2,a),coad(-1+Nc^2,b))");
}

#[test]
fn form_fierz_generator_contraction() {
    initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );

    assert_bare_snapshot!(expr.simplify_color().expand().simplify_metrics(),@"-1*Nc^(-1)*TR*g(cof(Nc,i),dind(cof(Nc,j)))*g(cof(Nc,k),dind(cof(Nc,l)))+TR*g(cof(Nc,i),dind(cof(Nc,l)))*g(cof(Nc,k),dind(cof(Nc,j)))"
    );
}

#[test]
#[ignore = "pending separated generator Casimir shortcut"]
fn form_separated_generator_casimir_shortcut() {
    initialize();
    let expr = parse_lit!(
        t(coad(Nc ^ 2 - 1, a), cof(Nc, i), dind(cof(Nc, j)))
            * t(coad(Nc ^ 2 - 1, b), cof(Nc, j), dind(cof(Nc, k)))
            * t(coad(Nc ^ 2 - 1, a), cof(Nc, k), dind(cof(Nc, l))),
        default_namespace = "spenso"
    );

    assert_bare_snapshot!(expr.simplify_color().simplify_metrics(),@"-1*Nc^(-1)*TR*t(coad(-1+Nc^2,b),cof(Nc,i),dind(cof(Nc,l)))"
    );
}

#[test]
#[ignore = "pending public chain-based color trace simplification"]
fn form_chain_one_generator_trace_vanishes() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_na, a)),
    );

    assert_bare_snapshot!(expr.simplify_color(),@"0");
}

#[test]
#[ignore = "pending public chain-based color trace simplification"]
fn form_chain_two_generator_trace_normalizes() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
    );

    assert_bare_snapshot!(expr.simplify_color(),@"TR*g(coad(NA,a),coad(NA,b))");
}

#[test]
#[ignore = "pending public chain-based three-generator color trace terminal"]
fn form_three_generator_trace_terminal() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, c)),
    );

    assert_bare_snapshot!(expr.simplify_color(),@"dR(coad(NA,a),coad(NA,b),coad(NA,c))+1/2*1𝑖*TR*f(coad(NA,a),coad(NA,b),coad(NA,c))"
    );
}

#[test]
#[ignore = "pending public chain-based adjacent color Casimir reduction"]
fn form_adjacent_generator_casimir_chain() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), k),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, a)),
    );

    assert_bare_snapshot!(expr.simplify_color(),@"CF*g(cof(Nc,i),dind(cof(Nc,k)))");
}

#[test]
#[ignore = "pending public trace-based separated color Casimir reduction"]
fn form_separated_generator_casimir_trace() {
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, c)),
    );

    assert_bare_snapshot!(expr.simplify_color(),@"(CF+-1/2*CA)*trace(cof(Nc),t(coad(NA,b),in,out),t(coad(NA,c),in,out))"
    );
}

#[test]
#[ignore = "pending structure-constant loop contraction"]
fn form_two_f_loop_contracts_to_ca_metric() {
    initialize();
    let expr = parse_lit!(
        f(coad(NA, a), coad(NA, c), coad(NA, d)) * f(coad(NA, b), coad(NA, c), coad(NA, d)),
        default_namespace = "spenso"
    );

    assert_bare_snapshot!(expr.simplify_color(),@"CA*g(coad(NA,a),coad(NA,b))");
}

#[test]
#[ignore = "pending structure-constant loop contraction"]
fn form_three_f_loop_contracts_to_ca_f() {
    initialize();
    let expr = parse_lit!(
        f(coad(NA, a), coad(NA, b), coad(NA, e))
            * f(coad(NA, b), coad(NA, c), coad(NA, f_))
            * f(coad(NA, c), coad(NA, a), coad(NA, g_)),
        default_namespace = "spenso"
    );

    assert_bare_snapshot!(expr.simplify_color(),@"1/2*CA*f(coad(NA,e),coad(NA,f_),coad(NA,g_))"
    );
}

#[test]
#[ignore = "pending public trace-f shortcut"]
fn form_mixed_trace_structure_contraction() {
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, d)),
    ) * parse_lit!(
        f(coad(NA, a), coad(NA, b), coad(NA, c)),
        default_namespace = "spenso"
    );

    assert_bare_snapshot!(expr.simplify_color(),@"1/2*1𝑖*CA*trace(cof(Nc),t(coad(NA,c),in,out),t(coad(NA,d),in,out))"
    );
}

#[test]
#[ignore = "pending symmetric invariant contraction"]
fn form_symmetric_invariant_d33_partial_contraction() {
    initialize();
    let expr = parse_lit!(
        d(sym_x, coad(NA, a), coad(NA, b), coad(NA, c))
            * d(sym_y, coad(NA, a), coad(NA, b), coad(NA, d_)),
        default_namespace = "spenso"
    );

    assert_bare_snapshot!(expr.simplify_color(),@"NA^(-1)*d33(sym_x,sym_y)*g(coad(NA,c),coad(NA,d_))"
    );
}

#[test]
#[ignore = "pending four-generator color trace terminal"]
fn form_four_generator_trace_terminal() {
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, c)),
        color_t!(slot!(r.coad_na, d)),
    );

    assert_bare_snapshot!(expr.simplify_color(),@"dR(coad(NA,a),coad(NA,b),coad(NA,c),coad(NA,d))+1/2*1𝑖*dR(coad(NA,a),coad(NA,b),coad(NA,x))*f(coad(NA,c),coad(NA,d),coad(NA,x))+1/2*1𝑖*dR(coad(NA,c),coad(NA,d),coad(NA,x))*f(coad(NA,a),coad(NA,b),coad(NA,x))+-1/6*TR*f(coad(NA,a),coad(NA,c),coad(NA,x))*f(coad(NA,b),coad(NA,d),coad(NA,x))+1/3*TR*f(coad(NA,a),coad(NA,d),coad(NA,x))*f(coad(NA,b),coad(NA,c),coad(NA,x))"
    );
}

#[test]
#[ignore = "pending color.h simpli strategy over invariant environments"]
fn form_simpli_contracts_projected_f_pair() {
    initialize();
    let expr = parse_lit!(
        f(coad(NA, x), coad(NA, a), coad(NA, c))
            * f(coad(NA, x), coad(NA, b), coad(NA, c))
            * invariant_environment(a, b),
        default_namespace = "spenso"
    );

    assert_bare_snapshot!(expr.simplify_color(),@"1/2*CA*projected_invariant(a,b,invariant_environment(a,b))"
    );
}

#[test]
#[ignore = "pending color.tar.gz tloop qloop size 3 family"]
fn form_tloop_qloop_size_3() {
    initialize();
    let result =
        parse_lit!(NA * TR * CF ^ 2 - 3 / 2 * NA * TR * CA * CF + 1 / 2 * NA * TR * CA ^ 2);

    assert_bare_snapshot!(
        result,
        @"1/2*CA^2*NA*TR+-3/2*CA*CF*NA*TR+CF^2*NA*TR"
    );
}

#[test]
#[ignore = "pending color.tar.gz tloop gloop size 3 family"]
fn form_tloop_gloop_size_3() {
    initialize();
    let result = parse_lit!(0);

    assert_bare_snapshot!(result, @"0");
}

#[test]
#[ignore = "pending color.tar.gz tloop qqloop size 3 family"]
fn form_tloop_qqloop_size_3() {
    initialize();
    let result = parse_lit!(-1 / 4 * NA * TR ^ 2 * CA + d33(R1, R2));

    assert_bare_snapshot!(
        result,
        @"-1/4*CA*NA*TR^2+d33(R1,R2)"
    );
}

#[test]
#[ignore = "pending color.tar.gz tloop qgloop size 3 family"]
fn form_tloop_qgloop_size_3() {
    initialize();
    let result = parse_lit!(1𝑖 / 4 * NA * TR * CA ^ 2);

    assert_bare_snapshot!(
        result,
        @"1/4*1𝑖*CA^2*NA*TR"
    );
}

#[test]
#[ignore = "pending color.tar.gz tloop ggloop size 3 family"]
fn form_tloop_ggloop_size_3() {
    initialize();
    let result = parse_lit!(1 / 4 * NA * CA ^ 3);

    assert_bare_snapshot!(result, @"1/4*CA^3*NA");
}

#[test]
#[ignore = "pending color.tar.gz fixed g14 family"]
fn form_tloop_g14() {
    initialize();
    let result = parse_lit!(
        1 / 648 * NA * CA ^ 7 - 8 / 15 * d444(A1, A2, A3) * CA + 16 / 9 * d644(A1, A2, A3)
    );

    assert_bare_snapshot!(
        result,
        @"1/648*CA^7*NA+-8/15*CA*d444(A1,A2,A3)+16/9*d644(A1,A2,A3)"
    );
}

#[test]
#[ignore = "pending color.tar.gz fixed fiveq family"]
fn form_tloop_fiveq() {
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

    assert_bare_snapshot!(
        result,
        @"1/192*CA^3*NA*TR^5+1/4*CA^2*TR^3*d33(R1,R2)+5/48*CA*NA^(-1)*TR*d33(R1,R2)*d33(R3,R4)+5/48*CA*NA^(-1)*TR*d33(R1,R3)*d33(R2,R4)+1/8*CA*NA^(-1)*TR*d33(R1,R4)*d33(R2,R3)+3/8*CA*TR^2*d433(R3,R1,R2)+1/2*CA*TR*d3333(R1,R2,R3,R4)+1/16*TR^4*d44(R1,A1)+d43333a(R5,R2,R1,R4,R3)"
    );
}

#[test]
#[ignore = "pending color.tar.gz su.frm F3F3 case"]
fn form_su_f3f3() {
    initialize();
    let result = parse_lit!(2 * a ^ 3 * NF ^ -1 - 5 / 2 * a ^ 3 * NF + 1 / 2 * a ^ 3 * NF ^ 3);

    assert_bare_snapshot!(
        result,
        @"2*NF^(-1)*a^3+-5/2*NF*a^3+1/2*NF^3*a^3"
    );
}

#[test]
#[ignore = "pending color.tar.gz su.frm F4F4 case"]
fn form_su_f4f4() {
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

    assert_bare_snapshot!(
        result,
        @"-3*NF^(-2)*a^4*nf^2+4*a^4*nf^2+-7/6*NF^2*a^4*nf^2+1/6*NF^4*a^4*nf^2"
    );
}

#[test]
#[ignore = "pending color.tar.gz su.frm F4A4 case"]
fn form_su_f4a4() {
    initialize();
    let result = parse_lit!(
        -2 * a ^ 4 * nf * NF + 5 / 3 * a ^ 4 * nf * NF ^ 3 + 1 / 3 * a ^ 4 * nf * NF ^ 5
    );

    assert_bare_snapshot!(
        result,
        @"-2*NF*a^4*nf+5/3*NF^3*a^4*nf+1/3*NF^5*a^4*nf"
    );
}

#[test]
#[ignore = "pending color.tar.gz su.frm A4A4 case"]
fn form_su_a4a4() {
    initialize();
    let result =
        parse_lit!(-24 * a ^ 4 * NF ^ 2 + 70 / 3 * a ^ 4 * NF ^ 4 + 2 / 3 * a ^ 4 * NF ^ 6);

    assert_bare_snapshot!(
        result,
        @"-24*NF^2*a^4+70/3*NF^4*a^4+2/3*NF^6*a^4"
    );
}

#[test]
#[ignore = "pending color.tar.gz su.frm FnFn n=3 case"]
fn form_su_fnfn_n3() {
    initialize();
    let result = parse_lit!(2 * a ^ 3 * nf ^ 2 * NF ^ -1 - 2 * a ^ 3 * nf ^ 2 * NF);

    assert_bare_snapshot!(
        result,
        @"2*NF^(-1)*a^3*nf^2+-2*NF*a^3*nf^2"
    );
}

#[test]
#[ignore = "pending color.tar.gz su.frm FnAn n=3 case"]
fn form_su_fnan_n3() {
    initialize();
    let result = parse_lit!(-1𝑖 * a ^ 3 * nf * NF ^ 2 + 1𝑖 * a ^ 3 * nf * NF ^ 4);

    assert_bare_snapshot!(
        result,
        @"-1𝑖*NF^2*a^3*nf+1𝑖*NF^4*a^3*nf"
    );
}

#[test]
#[ignore = "pending color.tar.gz su.frm AnAn n=3 case"]
fn form_su_anan_n3() {
    initialize();
    let result = parse_lit!(-2 * a ^ 3 * NF ^ 3 + 2 * a ^ 3 * NF ^ 5);

    assert_bare_snapshot!(
        result,
        @"-2*NF^3*a^3+2*NF^5*a^3"
    );
}

#[test]
#[ignore = "pending color.tar.gz su.frm qloop n=3 case"]
fn form_su_qloop_n3() {
    initialize();
    let result = parse_lit!(-a ^ 3 * nf * NF ^ -2 + a ^ 3 * nf * NF ^ 2);

    assert_bare_snapshot!(
        result,
        @"-1*NF^(-2)*a^3*nf+NF^2*a^3*nf"
    );
}

#[test]
#[ignore = "pending color.tar.gz su.frm gloop n=3 case"]
fn form_su_gloop_n3() {
    initialize();
    let result = parse_lit!(0);

    assert_bare_snapshot!(result, @"0");
}

fn colored_matrix_element() -> (Atom, Atom) {
    (
        parse!(
            "-G^2
            *(
                -g(mink(D,5),mink(D,6))*Q(2,mink(D,7))
                +g(mink(D,5),mink(D,6))*Q(3,mink(D,7))
                +g(mink(D,5),mink(D,7))*Q(2,mink(D,6))
                +g(mink(D,5),mink(D,7))*Q(4,mink(D,6))
                -g(mink(D,6),mink(D,7))*Q(3,mink(D,5))
                -g(mink(D,6),mink(D,7))*Q(4,mink(D,5))
            )
            *g(mink(D,4),mink(D,7))
            *t(coad(Nc^2-1,6),cof(Nc,5),dind(cof(Nc,4)))
            *f(coad(Nc^2-1,7),coad(Nc^2-1,8),coad(Nc^2-1,9))
            *g(bis(D,0),bis(D,5))
            *g(bis(D,1),bis(D,4))
            *g(mink(D,2),mink(D,5))
            *g(mink(D,3),mink(D,6))
            *g(coad(Nc^2-1,2),coad(Nc^2-1,7))
            *g(coad(Nc^2-1,3),coad(Nc^2-1,8))
            *g(coad(Nc^2-1,6),coad(Nc^2-1,9))
            *g(cof(Nc,0),dind(cof(Nc,5)))
            *g(cof(Nc,4),dind(cof(Nc,1)))
            *gamma(bis(D,5),bis(D,4),mink(D,4))
            *vbar(1,bis(D,1))
            *u(0,bis(D,0))
            *ϵbar(2,mink(D,2))
            *ϵbar(3,mink(D,3))",
            default_namespace = "spenso"
        ),
        parse_lit!(
            -4 * TR
                ^ 2 * Nc * G
                ^ 4 * (Nc - 1) * (Nc + 1) * (D - 2)
                ^ -2 * (-2 * dot(Q(0), Q(1)) * dot(Q(2), Q(2)) + dot(Q(0), Q(1)) * dot(Q(2), Q(3))
                    - 3 * dot(Q(0), Q(1)) * dot(Q(2), Q(4))
                    - 2 * dot(Q(0), Q(1)) * dot(Q(3), Q(3))
                    - 3 * dot(Q(0), Q(1)) * dot(Q(3), Q(4))
                    - 3 * dot(Q(0), Q(1)) * dot(Q(4), Q(4))
                    + 2 * dot(Q(0), Q(2)) * dot(Q(1), Q(2))
                    - dot(Q(0), Q(2)) * dot(Q(1), Q(3))
                    + dot(Q(0), Q(2)) * dot(Q(1), Q(4))
                    - dot(Q(0), Q(3)) * dot(Q(1), Q(2))
                    + 2 * dot(Q(0), Q(3)) * dot(Q(1), Q(3))
                    + dot(Q(0), Q(3)) * dot(Q(1), Q(4))
                    + dot(Q(0), Q(4)) * dot(Q(1), Q(2))
                    + dot(Q(0), Q(4)) * dot(Q(1), Q(3))
                    + 2 * dot(Q(0), Q(4)) * dot(Q(1), Q(4))
                    + D * dot(Q(0), Q(1)) * dot(Q(2), Q(2))
                    - D * dot(Q(0), Q(1)) * dot(Q(2), Q(3))
                    + D * dot(Q(0), Q(1)) * dot(Q(2), Q(4))
                    + D * dot(Q(0), Q(1)) * dot(Q(3), Q(3))
                    + D * dot(Q(0), Q(1)) * dot(Q(3), Q(4))
                    + D * dot(Q(0), Q(1)) * dot(Q(4), Q(4))
                    - D * dot(Q(0), Q(2)) * dot(Q(1), Q(2))
                    + D * dot(Q(0), Q(2)) * dot(Q(1), Q(3))
                    + D * dot(Q(0), Q(3)) * dot(Q(1), Q(2))
                    - D * dot(Q(0), Q(3)) * dot(Q(1), Q(3))),
            default_namespace = "spenso"
        ),
    )
}

#[test]
fn compact_printing() {
    initialize();

    let (atom, _) = colored_matrix_element();
    println!(
        "{}",
        atom.printer(SpensoPrintSettings::compact().nice_symbolica())
    );

    println!(
        "{}",
        atom.printer(
            SpensoPrintSettings {
                parens: true,
                with_dim: false,
                commas: false,
                index_subscripts: false,
                symbol_scripts: false,
            }
            .nice_symbolica()
        )
    );

    println!(
        "{}",
        atom.printer(SpensoPrintSettings::typst().nice_symbolica())
    )
}

#[test]
fn t_structure() {
    println!("{}", CS.t_strct::<AbstractIndex>(3, 8));

    let _ = Atom::Zero.simplify_metrics();
}

#[test]
fn test_val() {
    initialize();
    let expr = parse_lit!(
        (G ^ 3
            * (g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                - g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9))))
            * g(dind(cof(Nc, 2)), cof(Nc, l(5)))
            * g(mink(4, l(0)), mink(4, l(6)))
            * g(mink(4, l(1)), mink(4, l(7)))
            * g(mink(4, l(4)), mink(4, l(8)))
            * g(mink(4, l(5)), mink(4, l(9)))
            * g(bis(4, l(2)), bis(4, l(5)))
            * g(bis(4, l(3)), bis(4, l(6)))
            * g(dind(cof(Nc, l(6))), cof(Nc, 3))
            * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, l(8)))
            * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, l(9)))
            * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, l(10)))
            * g(coad(Nc ^ 2 - 1, l(7)), coad(Nc ^ 2 - 1, l(11)))
            * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
            * t(coad(Nc ^ 2 - 1, l(7)), cof(Nc, l(6)), dind(cof(Nc, l(5))))
            * f(
                coad(Nc ^ 2 - 1, l(8)),
                coad(Nc ^ 2 - 1, l(11)),
                coad(Nc ^ 2 - 1, l(12))
            )
            * f(
                coad(Nc ^ 2 - 1, l(9)),
                coad(Nc ^ 2 - 1, l(10)),
                coad(Nc ^ 2 - 1, l(12))
            )
            * ubar(2, bis(4, l(2)))
            * v(3, bis(4, l(3)))
            * ϵ(0, mink(4, l(0)))
            * ϵ(1, mink(4, l(1)))
            * ϵbar(4, mink(4, l(4)))
            + G
            ^ 3 * (g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                - g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                * g(dind(cof(Nc, 2)), cof(Nc, l(5)))
                * g(mink(4, l(0)), mink(4, l(6)))
                * g(mink(4, l(1)), mink(4, l(7)))
                * g(mink(4, l(4)), mink(4, l(8)))
                * g(mink(4, l(5)), mink(4, l(9)))
                * g(bis(4, l(2)), bis(4, l(5)))
                * g(bis(4, l(3)), bis(4, l(6)))
                * g(dind(cof(Nc, l(6))), cof(Nc, 3))
                * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, l(8)))
                * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, l(9)))
                * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, l(10)))
                * g(coad(Nc ^ 2 - 1, l(7)), coad(Nc ^ 2 - 1, l(11)))
                * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                * t(coad(Nc ^ 2 - 1, l(7)), cof(Nc, l(6)), dind(cof(Nc, l(5))))
                * f(
                    coad(Nc ^ 2 - 1, l(8)),
                    coad(Nc ^ 2 - 1, l(10)),
                    coad(Nc ^ 2 - 1, l(12))
                )
                * f(
                    coad(Nc ^ 2 - 1, l(9)),
                    coad(Nc ^ 2 - 1, l(11)),
                    coad(Nc ^ 2 - 1, l(12))
                )
                * ubar(2, bis(4, l(2)))
                * v(3, bis(4, l(3)))
                * ϵ(0, mink(4, l(0)))
                * ϵ(1, mink(4, l(1)))
                * ϵbar(4, mink(4, l(4)))
                + G
            ^ 3 * (g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9)))
                - g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                * g(dind(cof(Nc, 2)), cof(Nc, l(5)))
                * g(mink(4, l(0)), mink(4, l(6)))
                * g(mink(4, l(1)), mink(4, l(7)))
                * g(mink(4, l(4)), mink(4, l(8)))
                * g(mink(4, l(5)), mink(4, l(9)))
                * g(bis(4, l(2)), bis(4, l(5)))
                * g(bis(4, l(3)), bis(4, l(6)))
                * g(dind(cof(Nc, l(6))), cof(Nc, 3))
                * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, l(8)))
                * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, l(9)))
                * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, l(10)))
                * g(coad(Nc ^ 2 - 1, l(7)), coad(Nc ^ 2 - 1, l(11)))
                * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                * t(coad(Nc ^ 2 - 1, l(7)), cof(Nc, l(6)), dind(cof(Nc, l(5))))
                * f(
                    coad(Nc ^ 2 - 1, l(8)),
                    coad(Nc ^ 2 - 1, l(9)),
                    coad(Nc ^ 2 - 1, l(12))
                )
                * f(
                    coad(Nc ^ 2 - 1, l(10)),
                    coad(Nc ^ 2 - 1, l(11)),
                    coad(Nc ^ 2 - 1, l(12))
                )
                * ubar(2, bis(4, l(2)))
                * v(3, bis(4, l(3)))
                * ϵ(0, mink(4, l(0)))
                * ϵ(1, mink(4, l(1)))
                * ϵbar(4, mink(4, l(4))))
            * (-G
                ^ 3 * (g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                    - g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9))))
                    * g(dind(cof(Nc, 3)), cof(Nc, r(6)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(dind(cof(Nc, r(5))), cof(Nc, 2))
                    * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, r(8)))
                    * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, r(9)))
                    * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, r(10)))
                    * g(coad(Nc ^ 2 - 1, r(7)), coad(Nc ^ 2 - 1, r(11)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    * t(coad(Nc ^ 2 - 1, r(7)), cof(Nc, r(5)), dind(cof(Nc, r(6))))
                    * f(
                        coad(Nc ^ 2 - 1, r(8)),
                        coad(Nc ^ 2 - 1, r(11)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * f(
                        coad(Nc ^ 2 - 1, r(9)),
                        coad(Nc ^ 2 - 1, r(10)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * u(2, bis(4, r(2)))
                    * vbar(3, bis(4, r(3)))
                    * ϵ(4, mink(4, r(4)))
                    * ϵbar(0, mink(4, r(0)))
                    * ϵbar(1, mink(4, r(1)))
                    - G
                ^ 3 * (g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                    - g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                    * g(dind(cof(Nc, 3)), cof(Nc, r(6)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(dind(cof(Nc, r(5))), cof(Nc, 2))
                    * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, r(8)))
                    * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, r(9)))
                    * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, r(10)))
                    * g(coad(Nc ^ 2 - 1, r(7)), coad(Nc ^ 2 - 1, r(11)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    * t(coad(Nc ^ 2 - 1, r(7)), cof(Nc, r(5)), dind(cof(Nc, r(6))))
                    * f(
                        coad(Nc ^ 2 - 1, r(8)),
                        coad(Nc ^ 2 - 1, r(10)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * f(
                        coad(Nc ^ 2 - 1, r(9)),
                        coad(Nc ^ 2 - 1, r(11)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * u(2, bis(4, r(2)))
                    * vbar(3, bis(4, r(3)))
                    * ϵ(4, mink(4, r(4)))
                    * ϵbar(0, mink(4, r(0)))
                    * ϵbar(1, mink(4, r(1)))
                    - G
                ^ 3 * (g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9)))
                    - g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                    * g(dind(cof(Nc, 3)), cof(Nc, r(6)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(dind(cof(Nc, r(5))), cof(Nc, 2))
                    * g(coad(Nc ^ 2 - 1, 0), coad(Nc ^ 2 - 1, r(8)))
                    * g(coad(Nc ^ 2 - 1, 1), coad(Nc ^ 2 - 1, r(9)))
                    * g(coad(Nc ^ 2 - 1, 4), coad(Nc ^ 2 - 1, r(10)))
                    * g(coad(Nc ^ 2 - 1, r(7)), coad(Nc ^ 2 - 1, r(11)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    * t(coad(Nc ^ 2 - 1, r(7)), cof(Nc, r(5)), dind(cof(Nc, r(6))))
                    * f(
                        coad(Nc ^ 2 - 1, r(8)),
                        coad(Nc ^ 2 - 1, r(9)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * f(
                        coad(Nc ^ 2 - 1, r(10)),
                        coad(Nc ^ 2 - 1, r(11)),
                        coad(Nc ^ 2 - 1, r(12))
                    )
                    * u(2, bis(4, r(2)))
                    * vbar(3, bis(4, r(3)))
                    * ϵ(4, mink(4, r(4)))
                    * ϵbar(0, mink(4, r(0)))
                    * ϵbar(1, mink(4, r(1)))),
        default_namespace = "spenso"
    );

    println!("{expr}");
    println!("Simplify_metrics");
    println!("{}", expr.clone().simplify_metrics());
    println!("Simplify_color");
    println!(
        "{:>}",
        expr.simplify_gamma()
            .simplify_color()
            .expand()
            .simplify_metrics()
            .to_dots()
    );
}

#[test]
fn ratio_simplify() {
    test_initialize();
    let expr = parse_lit!(
        G ^ 4 * ee
            ^ 2 * f(
                coad(ohoho, dummy(0)),
                coad(ohoho, dummy(1)),
                coad(ohoho, dummy(2))
            ) * t(
                coad(ohoho, dummy(0)),
                cof(ahaha, dummy(3)),
                dind(cof(ahaha, dummy(4)))
            ) * t(
                coad(ohoho, dummy(1)),
                cof(ahaha, dummy(5)),
                dind(cof(ahaha, dummy(3)))
            ) * t(
                coad(ohoho, dummy(2)),
                cof(ahaha, dummy(4)),
                dind(cof(ahaha, dummy(5)))
            ),
        default_namespace = "spenso"
    );

    insta::assert_snapshot!(
        expr.cook_indices()
            .simplify_color()
            .to_bare_ordered_string(),@"(-1𝑖*TR^2*ahaha^(-1)*g(cof(ahaha,dummy_3),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,dummy_4),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_5),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_5)))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_3)))+-1𝑖*TR^2*ahaha^(-1)*g(cof(ahaha,dummy_3),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_4),dind(cof(ahaha,dummy_5)))*g(cof(ahaha,dummy_5),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_3)))+-1𝑖*TR^2*ahaha^(-1)*g(cof(ahaha,dummy_3),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_4),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_5),dind(cof(ahaha,dummy_3)))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_5)))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))+-1𝑖*TR^2*ahaha^(-2)*g(cof(ahaha,dummy_3),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,dummy_4),dind(cof(ahaha,dummy_5)))*g(cof(ahaha,dummy_5),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_3)))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))+-1𝑖*TR^2*ahaha^(-2)*g(cof(ahaha,dummy_3),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,dummy_4),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_5),dind(cof(ahaha,dummy_3)))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_5)))+-1𝑖*TR^2*g(cof(ahaha,dummy_3),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_4),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_5),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_3)))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_5)))+1𝑖*TR^2*ahaha^(-1)*g(cof(ahaha,dummy_3),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,dummy_4),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_5),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_3)))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_5)))+1𝑖*TR^2*ahaha^(-1)*g(cof(ahaha,dummy_3),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_4),dind(cof(ahaha,dummy_5)))*g(cof(ahaha,dummy_5),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_3)))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))+1𝑖*TR^2*ahaha^(-1)*g(cof(ahaha,dummy_3),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_4),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_5),dind(cof(ahaha,dummy_3)))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_5)))+1𝑖*TR^2*ahaha^(-2)*g(cof(ahaha,dummy_3),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,dummy_4),dind(cof(ahaha,dummy_5)))*g(cof(ahaha,dummy_5),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_3)))+1𝑖*TR^2*ahaha^(-2)*g(cof(ahaha,dummy_3),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,dummy_4),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_5),dind(cof(ahaha,dummy_3)))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_5)))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))+1𝑖*TR^2*g(cof(ahaha,dummy_3),dind(cof(ahaha,j(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_4),dind(cof(ahaha,k(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,dummy_5),dind(cof(ahaha,i(dummy_0,dummy_1,dummy_2))))*g(cof(ahaha,i(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_4)))*g(cof(ahaha,j(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_5)))*g(cof(ahaha,k(dummy_0,dummy_1,dummy_2)),dind(cof(ahaha,dummy_3))))*G^4*ee^2"
    );
}

#[test]
fn simple() {
    test_initialize();
    // tc₄ г₁₀ г₈·tc₆ г₈ г₁₀·fc₄ c₀ c₂·fc₆ c₂ c₁
    let expr = parse_lit!(
        f(
            coad(Nc ^ 2 - 1, c4),
            coad(Nc ^ 2 - 1, c0),
            coad(Nc ^ 2 - 1, c2)
        ) * f(
            coad(Nc ^ 2 - 1, c6),
            coad(Nc ^ 2 - 1, c2),
            coad(Nc ^ 2 - 1, c0)
        ) * t(coad(Nc ^ 2 - 1, c4), cof(Nc, i0), dind(cof(Nc, i2)))
            * t(coad(Nc ^ 2 - 1, c6), cof(Nc, i2), dind(cof(Nc, i0))),
        default_namespace = "spenso"
    );

    println!("{}", expr.simplify_color());
    println!(
        "{}",
        expr.simplify_color()
            .replace(parse_lit!(spenso::Nc))
            .with(3)
    );
}

#[test]
fn minus_sign() {
    test_initialize();
    let expr1 = parse_lit!(
        f(coad(8, hedge(5)), coad(8, hedge(9)), coad(8, hedge(13)))
            * g(coad(8, hedge(12)), coad(8, hedge(13)))
            * g(coad(8, hedge(4)), coad(8, hedge(5)))
            * g(coad(8, hedge(8)), coad(8, hedge(9)))
            * g(cof(3, hedge(10)), dind(cof(3, hedge(11))))
            * g(cof(3, hedge(11)), dind(cof(3, hedge(15))))
            * g(cof(3, hedge(15)), dind(cof(3, hedge(14))))
            * g(cof(3, hedge(16)), dind(cof(3, hedge(17))))
            * g(cof(3, hedge(17)), dind(cof(3, hedge(7))))
            * g(cof(3, hedge(2)), dind(cof(3, hedge(3))))
            * g(cof(3, hedge(7)), dind(cof(3, hedge(6))))
            * t(
                coad(8, hedge(12)),
                cof(3, hedge(14)),
                dind(cof(3, hedge(16)))
            )
            * t(coad(8, hedge(4)), cof(3, hedge(6)), dind(cof(3, hedge(2))))
            * t(coad(8, hedge(8)), cof(3, hedge(3)), dind(cof(3, hedge(10)))),
        default_namespace = "spenso"
    );

    let expr2 = parse_lit!(
        f(coad(8, hedge(5)), coad(8, hedge(9)), coad(8, hedge(13)))
            * g(coad(8, hedge(12)), coad(8, hedge(13)))
            * g(coad(8, hedge(4)), coad(8, hedge(5)))
            * g(coad(8, hedge(8)), coad(8, hedge(9)))
            * g(cof(3, hedge(11)), dind(cof(3, hedge(10))))
            * g(cof(3, hedge(14)), dind(cof(3, hedge(15))))
            * g(cof(3, hedge(15)), dind(cof(3, hedge(11))))
            * g(cof(3, hedge(17)), dind(cof(3, hedge(16))))
            * g(cof(3, hedge(3)), dind(cof(3, hedge(2))))
            * g(cof(3, hedge(6)), dind(cof(3, hedge(7))))
            * g(cof(3, hedge(7)), dind(cof(3, hedge(17))))
            * t(
                coad(8, hedge(12)),
                cof(3, hedge(16)),
                dind(cof(3, hedge(14)))
            )
            * t(coad(8, hedge(4)), cof(3, hedge(2)), dind(cof(3, hedge(6))))
            * t(coad(8, hedge(8)), cof(3, hedge(10)), dind(cof(3, hedge(3)))),
        default_namespace = "spenso"
    );

    // println!(
    //     "{}",
    //     (expr1.cook_indices().canonize(AbstractIndex::Dummy)
    //         / (expr2.cook_indices().canonize(AbstractIndex::Dummy)))
    //     .cancel()
    // );
    println!(
        "{}\n",
        expr1
            .simplify_metrics()
            .cook_indices()
            .canonize(AbstractIndex::Dummy)
    );
    println!(
        "{}\n",
        expr2
            .simplify_metrics()
            .cook_indices()
            .canonize(AbstractIndex::Dummy)
    );
    assert!((expr1.simplify_color() + expr2.simplify_color()).is_zero());
}

mod failing {
    use super::*;

    #[test]
    #[should_panic]
    fn test_color_matrix_element() {
        initialize();
        let _q = symbol!("spenso::Q";Real);
        let spin_sum_rule = parse!(
            "
                g(coad(Nc^2-1, left(3)), coad(Nc^2-1, right(3)))
                    * g(coad(Nc^2-1, left(2)), coad(Nc^2-1, right(2)))
                    * g(cof(Nc, right(0)), dind(cof(Nc, left(0))))
                    * g(cof(Nc, left(1)), dind(cof(Nc, right(1))))",
            default_namespace = "spenso"
        );

        let amplitude_color = parse!(
            "
                t(coad(Nc^2-1, 6), cof(Nc, 5), dind(cof(Nc, 4)))
                    * f(coad(Nc^2-1, 7), coad(Nc^2-1, 8), coad(Nc^2-1, 9))
                    * g(coad(Nc^2-1, 2), coad(Nc^2-1, 7))
                    * g(coad(Nc^2-1, 3), coad(Nc^2-1, 8))
                    * g(coad(Nc^2-1, 6), coad(Nc^2-1, 9))
                    * g(cof(Nc, 0), dind(cof(Nc, 5)))
                    * g(cof(Nc, 4), dind(cof(Nc, 1)))",
            default_namespace = "spenso"
        );
        let amplitude_color_left = amplitude_color.wrap_indices(symbol!("spenso::left"));

        // return;
        let amplitude_color_right = amplitude_color
            .dirac_adjoint::<AbstractIndex>()
            .unwrap()
            .wrap_indices(symbol!("spenso::right"));
        println!("left{amplitude_color_left}");

        println!("right{amplitude_color_right}");
        let amp_squared_color = amplitude_color_left * spin_sum_rule * amplitude_color_right;
        let simplified_color = amp_squared_color.simplify_metrics().simplify_color();
        println!("simplified_color={}", simplified_color);

        let spin_sum_rule_src = parse_lit!(
            spenso::vbar(1, bis(D, left(1)))
                * spenso::v(1, bis(D, right(1)))
                * spenso::u(0, bis(D, left(0)))
                * spenso::ubar(0, bis(D, right(0)))
                * spenso::ϵbar(2, mink(D, left(2)))
                * spenso::ϵ(2, mink(D, right(2)))
                * spenso::ϵbar(3, mink(D, left(3)))
                * spenso::ϵ(3, mink(D, right(3))),
            default_namespace = "spenso"
        );

        let spin_sum_rule_trg = parse!(
            "
                1/4*1/(D-2)^2*
                (
                    (-1) * gamma(bis(D,left(1)),bis(D,right(1)),mink(D,1337))*Q(1,mink(D,1337))
                    * gamma(bis(D,right(0)),bis(D,left(0)),mink(D,1338))*Q(0,mink(D,1338))
                    * (-1) * g(mink(D,left(2)),mink(D,right(2)))
                    * (-1) * g(mink(D,left(3)),mink(D,right(3)))
                    )
                    * (
                        g(coad(Nc^2-1, left(3)), coad(Nc^2-1, right(3)))
                        * g(coad(Nc^2-1, left(2)), coad(Nc^2-1, right(2)))
                        * g(cof(Nc, right(0)), dind(cof(Nc, left(0))))
                        * g(cof(Nc, left(1)), dind(cof(Nc, right(1))))
                    )",
            default_namespace = "spenso"
        );

        let (amplitude, tgt) = colored_matrix_element();

        let amplitude_left = amplitude.wrap_indices(symbol!("spenso::left"));

        println!("Amplitude left:\n{}", amplitude_left.collect_factors());

        println!(
            "Amplitude left cooked:\n{}",
            amplitude_left.collect_factors().cook_indices()
        );
        let amplitude_right = amplitude.wrap_indices(symbol!("spenso::right"));

        println!("Amplitude right:\n{}", amplitude_right.factor());

        println!(
            "Amplitude right conj:\n{}",
            amplitude_right
                .dirac_adjoint::<AbstractIndex>()
                .unwrap()
                .factor()
        );

        let mut amp_squared =
            amplitude_left * amplitude_right.dirac_adjoint::<AbstractIndex>().unwrap();

        println!("Amplitude squared:\n{}", amp_squared.factor());

        // let dangling_atoms = amp_squared.list_dangling::<AbstractIndex>();

        // assert_eq!(dangling_atoms.len(), 8);

        amp_squared = amp_squared
            .expand()
            .replace(spin_sum_rule_src.to_pattern())
            .with(spin_sum_rule_trg.to_pattern());

        println!("Amplitude squared spin-summed:\n{}", amp_squared);

        let mut simplified_amp_squared = amp_squared.clone();

        simplified_amp_squared = simplified_amp_squared.simplify_color();

        println!(
            "Color-simplified amplitude squared:\n{}",
            simplified_amp_squared
        );

        simplified_amp_squared = simplified_amp_squared.simplify_gamma();

        println!(
            "Gamma+color-simplified amplitude squared:\n{}",
            simplified_amp_squared
        );

        simplified_amp_squared = simplified_amp_squared.to_dots();

        assert_eq!(
            tgt,
            simplified_amp_squared.factor(),
            "{}\nnot equal to\n{}",
            tgt,
            simplified_amp_squared.factor()
        );
    }

    #[test]
    #[should_panic]
    fn test_color_matrix_element_two() {
        initialize();

        let spin_sum_rule_src = parse_lit!(
            vbar(1, bis(D, left(1)))
                * v(1, bis(D, right(1)))
                * u(0, bis(D, left(0)))
                * ubar(0, bis(D, right(0)))
                * ϵbar(2, mink(D, left(2)))
                * ϵ(2, mink(D, right(2)))
                * ϵbar(3, mink(D, left(3)))
                * ϵ(3, mink(D, right(3))),
            default_namespace = "spenso"
        );

        let spin_sum_rule_trg = parse!(
            "
                1/4*1/(D-2)^2*
                (
                    (-1) * gamma(bis(D,left(1)),bis(D,right(1)),mink(D,1337))*Q(1,mink(D,1337))
                    * gamma(bis(D,right(0)),bis(D,left(0)),mink(D,1338))*Q(0,mink(D,1338))
                    * (-1) * g(mink(D,left(2)),mink(D,right(2)))
                    * (-1) * g(mink(D,left(3)),mink(D,right(3)))
                    )
                    ",
            default_namespace = "spenso"
        );

        let (amplitude, tgt) = colored_matrix_element();

        let amplitude_left = amplitude.wrap_dummies::<AbstractIndex>(symbol!("spenso::left"));

        println!("Amplitude left:\n{}", amplitude_left.collect_factors());

        let amplitude_right = amplitude.wrap_dummies::<AbstractIndex>(symbol!("spenso::right"));

        println!("Amplitude right:\n{}", amplitude_right.conj().factor());

        let mut amp_squared = amplitude_left * amplitude_right.conj();

        println!("Amplitude squared:\n{}", amp_squared.factor());

        let _spin_sum_pat = parse!(
            "gamma(bis(D,left(1)),bis(D,right(1)),mink(D,1337))",
            default_namespace = "spenso"
        )
        .to_pattern();

        let pols_pats: Vec<Pattern> = vec![
            function!(PS.vbar, RS.x__).into(),
            function!(PS.v, RS.x__).into(),
            function!(PS.u, RS.x__).into(),
            function!(PS.ubar, RS.x__).into(),
            function!(PS.eps, RS.x__).into(),
            function!(PS.ebar, RS.x__).into(),
        ];

        let mut pols_coefs = amp_squared.expand_in_patterns(&pols_pats);

        for (c, _) in &mut pols_coefs {
            println!("c:{c}");
            let r = c
                .replace(spin_sum_rule_src.to_pattern())
                .with(spin_sum_rule_trg.to_pattern());
            println!("r:{r}");
            *c = r;
        }

        amp_squared = pols_coefs.iter().fold(Atom::Zero, |a, (c, s)| a + c * s);

        println!("Amplitude squared spin-summed:\n{}", amp_squared.factor());

        let mut simplified_amp_squared = amp_squared.clone();

        simplified_amp_squared = simplified_amp_squared.simplify_color();

        println!(
            "Color-simplified amplitude squared:\n{}",
            simplified_amp_squared
        );

        simplified_amp_squared = simplified_amp_squared.simplify_gamma().simplify_metrics();

        println!(
            "Gamma+color-simplified amplitude squared:\n{}",
            simplified_amp_squared
        );

        simplified_amp_squared = simplified_amp_squared.to_dots();

        assert_eq!(
            tgt,
            simplified_amp_squared.factor(),
            "{:#}\nnot equal to\n{:#}",
            tgt,
            simplified_amp_squared.factor()
        );
    }
}
