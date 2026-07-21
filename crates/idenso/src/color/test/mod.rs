use insta::assert_snapshot;
use spenso::network::parsing::{ParseSettings, StructureFromAtom};
use spenso::network::tags::SPENSO_TAG;
use spenso::structure::IndexlessNamedStructure;
use spenso::structure::PermutedStructure;
use symbolica_utils::AtomPrintExt;

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

use spenso::{
    antisym, chain, network::parsing::ShadowedStructure, s, slot, structure::permuted::Perm, sym,
    trace,
};
use symbolica::{id::Pattern, parse, parse_lit};

use crate::dirac::PS;
use crate::selective_expand::SelectiveExpand;
use crate::shorthands::schoonschip::Schoonschip;
use crate::tensor::{SymbolicNetExt, SymbolicNetParse, SymbolicTensor};
use crate::{Cookable, IndexTooling, dirac::GammaSimplifier, shorthands::metric::MetricSimplifier};
use crate::{color_cas, color_f, color_gram, color_idx, color_str_t, color_t, f};

use super::*;

use crate::test_support::{TestReps, test_initialize};

fn assert_color_zero(expr: Atom) {
    assert!(expr.simplify_color().is_zero());
}

#[test]
fn test_color_structures() {
    test_initialize();
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

    let f_p = f.clone().permute_inds();

    // println!("{}", f_p);
    let simplified = f_p.expression.schoonschip();
    // println!("{}", simplified);
    let f_parsed = ShadowedStructure::<AbstractIndex>::parse(simplified.as_view()).unwrap();

    assert_eq!(f.index_permutation, f_parsed.index_permutation);
    assert!(f_parsed.rep_permutation.is_identity());

    let f_p = f.clone().permute_reps_wrapped().permute_inds();

    let simplified = f_p.expression.schoonschip();
    // println!("{}", simplified);
    let f_parsed = ShadowedStructure::<AbstractIndex>::parse(simplified.as_view()).unwrap();

    assert_eq!(f.index_permutation, f_parsed.index_permutation);
    assert_eq!(f.rep_permutation, f_parsed.rep_permutation);
}

#[test]
fn test_color_simplification() {
    test_initialize();

    let atom = parse_lit!(
        f(
            coad(Nc ^ 2 - 1, 1),
            coad(Nc ^ 2 - 1, 2),
            coad(Nc ^ 2 - 1, 3)
        ) ^ 2,
        default_namespace = "spenso"
    );
    let simplified = atom.simplify_color();

    assert_snapshot!(simplified.to_bare_ordered_string(), @"(-1+Nc^2)*cas(2,coad(-1+Nc^2))");
}

#[test]
fn two_fs() {
    test_initialize();

    let atom = f!(3, 1, 5) * f!(3, 5, 1);
    let simplified = atom.simplify_color();

    assert_snapshot!(simplified.to_bare_ordered_string(), @"(-1+Nc^2)*-1*cas(2,coad(-1+Nc^2))");
}

#[test]
fn color_generator_macro_accepts_gamma_style_index_shorthands() {
    test_initialize();
    let r = TestReps::new();
    let cof_nc = ColorFundamental {}.new_rep(CS.nc);

    let a = slot!(r.coad_na, a);
    assert_eq!(color_t!(a), color_t!(slot!(r.coad_na, a)));
    assert_snapshot!(color_t!(1).to_bare_ordered_string(), @"t(1,in,out)");

    let a = slot!(r.coad_na, a);
    let i = slot!(cof_nc, i);
    let j = slot!(cof_nc.dual(), j);
    assert_eq!(
        color_t!(a, i, j),
        color_t!(
            slot!(r.coad_na, a),
            slot!(cof_nc, i),
            slot!(cof_nc.dual(), j)
        ),
    );
    assert_eq!(
        color_t!(RS.a__, RS.i__, RS.j__),
        CS.explicit_t(
            ColorAdjoint {}.to_symbolic([RS.a__]),
            ColorFundamental {}.to_symbolic([RS.i__]),
            ColorAntiFundamental {}.to_symbolic([RS.j__]),
        ),
    );
    assert_eq!(
        color_t!([CS.adj_, RS.a_], [CS.nc_, RS.i_], [CS.nc_, RS.j_]),
        CS.t_pattern(CS.nc_, CS.adj_, RS.a_, RS.i_, RS.j_),
    );
}

#[test]
fn color_structure_macro_accepts_gamma_style_index_shorthands() {
    test_initialize();
    let r = TestReps::new();

    let a = slot!(r.coad_na, a);
    let b = slot!(r.coad_na, b);
    let c = slot!(r.coad_na, c);
    assert_eq!(
        color_f!(a, b, c),
        color_f!(
            slot!(r.coad_na, a),
            slot!(r.coad_na, b),
            slot!(r.coad_na, c)
        ),
    );
    assert_eq!(
        color_f!(RS.a__, RS.b__, RS.c__),
        CS.structure_f(
            ColorAdjoint {}.to_symbolic([RS.a__]),
            ColorAdjoint {}.to_symbolic([RS.b__]),
            ColorAdjoint {}.to_symbolic([RS.c__]),
        ),
    );
    assert_eq!(
        color_f!([CS.adj_, RS.a_], [CS.adj_, RS.b_], [CS.adj_, RS.c_]),
        CS.f_pattern(CS.adj_, RS.a_, RS.b_, RS.c_),
    );
}

#[test]
fn color_invariant_macros_build_scalar_heads() {
    test_initialize();
    let n = Atom::var(CS.nc);
    let default_adjoint_dimension = n.clone().pow(Atom::num(2)) - Atom::num(1);
    let cof_n = ColorFundamental {}.to_symbolic([n.clone()]);
    let a = ColorAdjoint {}.to_symbolic([default_adjoint_dimension.clone(), Atom::var(s!(a))]);
    let b = ColorAdjoint {}.to_symbolic([default_adjoint_dimension, Atom::var(s!(b))]);

    assert_eq!(
        color_cas!(2, cof_n.clone()),
        CS.cas(Atom::num(2), cof_n.clone())
    );
    assert_eq!(
        color_idx!(2, cof_n.clone()),
        CS.idx(Atom::num(2), cof_n.clone())
    );
    assert_eq!(
        color_gram!(3, cof_n.clone(), cof_n.clone()),
        CS.gram(Atom::num(3), cof_n.clone(), cof_n.clone())
    );
    assert_eq!(
        color_str_t!(cof_n.clone(), a.clone(), b.clone()),
        trace!(
            cof_n.clone(),
            sym!(color_t!(a.clone()), color_t!(b.clone()))
        )
    );

    let expr = color_cas!(2, cof_n.clone())
        * color_idx!(2, cof_n.clone())
        * color_gram!(3, cof_n.clone(), cof_n);
    let net = expr
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert!(net.graph.dangling_indices().is_empty());
    assert_eq!(net.simple_execute::<()>(), expr);
}

#[test]
fn color_invariant_print_special_cases_are_compact() {
    test_initialize();
    let cof_n = ColorFundamental {}.to_symbolic([Atom::var(CS.nc)]);

    let mut compact = SpensoPrintSettings::compact().nice_symbolica();
    compact.color_builtin_symbols = false;
    assert_eq!(
        color_cas!(2, cof_n.clone()).printer(compact).to_string(),
        "CF"
    );

    let typst = SpensoPrintSettings::typst().typst_symbolica();
    assert_eq!(color_cas!(2, cof_n).printer(typst).to_string(), "C_F");
}

#[test]
fn cof_dimension_invariant_rules_substitute_supported_fundamental_cases() {
    test_initialize();
    let n = Atom::var(CS.nc);
    let default_adjoint_dimension = n.clone().pow(Atom::num(2)) - Atom::num(1);
    let cof_n = ColorFundamental {}.to_symbolic([n.clone()]);
    let coad_n = ColorAdjoint {}.to_symbolic([default_adjoint_dimension]);

    assert_eq!(
        color_idx!(2, cof_n.clone()).to_cof_dimension_invariants(),
        Atom::num(1) / Atom::num(2)
    );
    assert_eq!(
        color_cas!(2, cof_n.clone())
            .to_cof_dimension_invariants()
            .expand(),
        ((n.clone().pow(Atom::num(2)) - Atom::num(1)) / (Atom::num(2) * n.clone())).expand()
    );
    assert_eq!(
        color_cas!(2, coad_n).to_cof_dimension_invariants(),
        n.clone()
    );
    assert_eq!(
        color_cas!(2, ColorAdjoint {}.to_symbolic([Atom::var(CS.na)]))
            .to_cof_dimension_invariants(),
        n.clone()
    );
    assert_eq!(
        color_cas!(2, ColorAdjoint {}.to_symbolic([Atom::num(8)])).to_cof_dimension_invariants(),
        Atom::num(3)
    );
    assert_eq!(
        (Atom::var(CS.ca) * Atom::var(CS.tr) * Atom::var(CS.na))
            .to_cof_dimension_invariants()
            .expand(),
        (n.clone() * (Atom::num(1) / Atom::num(2)) * (n.clone().pow(Atom::num(2)) - Atom::num(1)))
            .expand()
    );

    let gram_three = color_gram!(3, cof_n.clone(), cof_n.clone())
        .to_cof_dimension_invariants()
        .expand();
    let expected_three = ((n.clone().pow(Atom::num(2)) - Atom::num(1))
        * (n.clone().pow(Atom::num(2)) - Atom::num(4))
        / (Atom::num(16) * n.clone()))
    .expand();
    assert_eq!(gram_three, expected_three);

    let gram_four = color_gram!(4, cof_n.clone(), cof_n)
        .to_cof_dimension_invariants()
        .expand();
    let n_squared = n.clone().pow(Atom::num(2));
    let expected_four = ((n_squared.clone() - Atom::num(1))
        * (n_squared.clone().pow(Atom::num(2)) - Atom::num(6) * n_squared + Atom::num(18))
        / (Atom::num(96) * n.pow(Atom::num(2))))
    .expand();
    assert_eq!(gram_four, expected_four);
}

#[test]
fn color_structure_symbol_is_antisymmetric() {
    test_initialize();
    let r = TestReps::new();

    let a = slot!(r.coad_na, a);
    let b = slot!(r.coad_na, b);
    let c = slot!(r.coad_na, c);

    assert_snapshot!(
        color_f!(a, c, b).to_bare_ordered_string(),
        @"-1*f(coad(NA,a),coad(NA,b),coad(NA,c))"
    );
    assert!(color_f!(a, b, b).is_zero());
}

#[test]
fn permuted_structure_constant_square_simplifies_with_sign() {
    test_initialize();
    let atom = parse_lit!(
        f(coad(8, a), coad(8, b), coad(8, c)) * f(coad(8, a), coad(8, c), coad(8, b)),
        default_namespace = "spenso"
    );

    assert_snapshot!(atom.simplify_color().to_bare_ordered_string(), @"-8*cas(2,coad(8))");
}

#[test]
fn three_loop_pole_part_color() {
    test_initialize();
    let input = parse_lit!(
        ((-16 + -26 * eps ^ 2 + -8 / 3 * eps ^ 2 * 𝜋 ^ 2 + 56 / 3 * eps)
            * f(
                coad(8, hedge(11)),
                coad(8, hedge(13)),
                coad(8, vertex(2, 1))
            )
            * f(coad(8, hedge(11)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
            * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
            * f(coad(8, hedge(4)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
            + (-16 + -26 * eps ^ 2 + -8 / 3 * eps ^ 2 * 𝜋 ^ 2 + 56 / 3 * eps)
                * f(
                    coad(8, hedge(11)),
                    coad(8, hedge(13)),
                    coad(8, vertex(3, 1))
                )
                * f(coad(8, hedge(11)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(6)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
            + (-16 + -26 * eps ^ 2 + -8 / 3 * eps ^ 2 * 𝜋 ^ 2 + 88 / 3 * eps)
                * f(coad(8, hedge(11)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
            + (-16 + -26 * eps ^ 2 + -8 / 3 * eps ^ 2 * 𝜋 ^ 2 + 88 / 3 * eps)
                * f(coad(8, hedge(11)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
            + (-16 / 3 * eps ^ 2 * 𝜋 ^ 2 + -32 + -52 * eps ^ 2 + 176 / 3 * eps)
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
            + (-16 / 3 * eps ^ 2 * 𝜋 ^ 2 + -32 + -52 * eps ^ 2 + 48 * eps)
                * f(
                    coad(8, hedge(11)),
                    coad(8, hedge(13)),
                    coad(8, vertex(2, 1))
                )
                * f(
                    coad(8, hedge(11)),
                    coad(8, hedge(13)),
                    coad(8, vertex(3, 1))
                )
                * f(coad(8, hedge(4)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(6)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
            + (-16 / 3 * eps ^ 2 * 𝜋 ^ 2 + -32 + -52 * eps ^ 2 + 48 * eps)
                * f(coad(8, hedge(11)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(11)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
            + (-88 / 3 * eps + 16 + 26 * eps ^ 2 + 8 / 3 * eps ^ 2 * 𝜋 ^ 2)
                * f(
                    coad(8, hedge(11)),
                    coad(8, hedge(13)),
                    coad(8, vertex(2, 1))
                )
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(6)), coad(8, vertex(3, 1)))
                * f(coad(8, hedge(4)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
            + (-88 / 3 * eps + 16 + 26 * eps ^ 2 + 8 / 3 * eps ^ 2 * 𝜋 ^ 2)
                * f(
                    coad(8, hedge(11)),
                    coad(8, hedge(13)),
                    coad(8, vertex(3, 1))
                )
                * f(coad(8, hedge(11)), coad(8, hedge(9)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(13)), coad(8, hedge(4)), coad(8, vertex(2, 1)))
                * f(coad(8, hedge(6)), coad(8, hedge(9)), coad(8, vertex(3, 1))))
            * 1
            / 64
            * CA
            * dot(P(0, mink(4)), P(0, mink(4)))
            * g(coad(8, hedge(4)), coad(8, hedge(6)))
            * gs
            ^ 6 * eps
            ^ (-3),
        default_namespace = "spenso"
    );

    let color_zero_candidate = input.cook_indices().simplify_color().collect_color();

    assert_snapshot!(&color_zero_candidate.collect_symbol::<i16>(SPENSO_TAG.dot).to_bare_ordered_string(),@"(((-16+-26*eps^2+-8/3*eps^2*𝜋^2+56/3*eps)*1/128*CA*eps^(-3)*gs^6+(-88/3*eps+16+26*eps^2+8/3*eps^2*𝜋^2)*-1/128*CA*eps^(-3)*gs^6)*16+(-16+-26*eps^2+-8/3*eps^2*𝜋^2+88/3*eps)*1/8*CA*eps^(-3)*gs^6+(-16/3*eps^2*𝜋^2+-32+-52*eps^2+176/3*eps)*1/8*CA*eps^(-3)*gs^6+(-16/3*eps^2*𝜋^2+-32+-52*eps^2+48*eps)*1/4*CA*eps^(-3)*gs^6)*(cas(2,coad(8)))^2*dot(P(0,mink(4)),P(0,mink(4)))");

    let input = parse_lit!(
        ((8 * eps + 8 / 3) * 1 / 64
            * f(coad(8, hedge(1)), coad(8, hedge(11)), coad(8, hedge(15)))
            * f(coad(8, hedge(1)), coad(8, hedge(3)), coad(8, hedge(5)))
            * f(coad(8, hedge(11)), coad(8, hedge(13)), coad(8, hedge(9)))
            * f(coad(8, hedge(13)), coad(8, hedge(5)), coad(8, hedge(7)))
            * f(coad(8, hedge(15)), coad(8, hedge(17)), coad(8, hedge(7)))
            + -1 / 16
                * f(coad(8, hedge(1)), coad(8, hedge(11)), coad(8, hedge(15)))
                * f(coad(8, hedge(1)), coad(8, hedge(3)), coad(8, hedge(4)))
                * f(coad(8, hedge(11)), coad(8, hedge(13)), coad(8, hedge(9)))
                * f(coad(8, hedge(13)), coad(8, hedge(4)), coad(8, hedge(7)))
                * f(coad(8, hedge(15)), coad(8, hedge(17)), coad(8, hedge(7)))
            + 1 / 16
                * f(coad(8, hedge(1)), coad(8, hedge(10)), coad(8, hedge(14)))
                * f(coad(8, hedge(1)), coad(8, hedge(3)), coad(8, hedge(5)))
                * f(coad(8, hedge(10)), coad(8, hedge(13)), coad(8, hedge(9)))
                * f(coad(8, hedge(13)), coad(8, hedge(5)), coad(8, hedge(7)))
                * f(coad(8, hedge(14)), coad(8, hedge(17)), coad(8, hedge(7))))
            * dot(P(0, mink(4)), P(0, mink(4)))
            * f(coad(8, hedge(17)), coad(8, hedge(3)), coad(8, hedge(9)))
            * gs
            ^ 6 * eps
            ^ (-2),
        default_namespace = "spenso"
    );

    let color_zero_candidate = input.cook_indices().simplify_color().collect_color();

    assert_snapshot!(&color_zero_candidate.collect_symbol::<i16>(SPENSO_TAG.dot).to_bare_ordered_string(),@"((8*eps+8/3)*1/64*eps^(-2)*f(coad(8,hedge_1),coad(8,hedge_11),coad(8,hedge_15))*f(coad(8,hedge_1),coad(8,hedge_3),coad(8,hedge_5))*f(coad(8,hedge_11),coad(8,hedge_13),coad(8,hedge_9))*f(coad(8,hedge_13),coad(8,hedge_5),coad(8,hedge_7))*f(coad(8,hedge_15),coad(8,hedge_17),coad(8,hedge_7))*f(coad(8,hedge_17),coad(8,hedge_3),coad(8,hedge_9))*gs^6+-1/16*eps^(-2)*f(coad(8,hedge_1),coad(8,hedge_11),coad(8,hedge_15))*f(coad(8,hedge_1),coad(8,hedge_3),coad(8,hedge_4))*f(coad(8,hedge_11),coad(8,hedge_13),coad(8,hedge_9))*f(coad(8,hedge_13),coad(8,hedge_4),coad(8,hedge_7))*f(coad(8,hedge_15),coad(8,hedge_17),coad(8,hedge_7))*f(coad(8,hedge_17),coad(8,hedge_3),coad(8,hedge_9))*gs^6+1/16*eps^(-2)*f(coad(8,hedge_1),coad(8,hedge_10),coad(8,hedge_14))*f(coad(8,hedge_1),coad(8,hedge_3),coad(8,hedge_5))*f(coad(8,hedge_10),coad(8,hedge_13),coad(8,hedge_9))*f(coad(8,hedge_13),coad(8,hedge_5),coad(8,hedge_7))*f(coad(8,hedge_14),coad(8,hedge_17),coad(8,hedge_7))*f(coad(8,hedge_17),coad(8,hedge_3),coad(8,hedge_9))*gs^6)*dot(P(0,mink(4)),P(0,mink(4)))")
}

#[test]
fn kaapo_gl34_color_input_simplifies_to_zero() {
    test_initialize();

    // Captured from feyngen GL34 after `to_param_color()` and before `simplify_color()`.
    let input = parse!(
        "
        -1𝑖 * UFO::G^6
        * (
            -Q(2, mink(4, hedge(11))) * g(mink(4, hedge(5)), mink(4, hedge(17)))
            + Q(2, mink(4, hedge(17))) * g(mink(4, hedge(5)), mink(4, hedge(11)))
            + Q(5, mink(4, hedge(5))) * g(mink(4, hedge(11)), mink(4, hedge(17)))
            - Q(5, mink(4, hedge(17))) * g(mink(4, hedge(5)), mink(4, hedge(11)))
            - Q(8, mink(4, hedge(5))) * g(mink(4, hedge(11)), mink(4, hedge(17)))
            + Q(8, mink(4, hedge(11))) * g(mink(4, hedge(5)), mink(4, hedge(17)))
        )
        * (
            -Q(6, mink(4, hedge(15))) * g(mink(4, hedge(13)), mink(4, hedge(16)))
            + Q(6, mink(4, hedge(16))) * g(mink(4, hedge(13)), mink(4, hedge(15)))
            + Q(7, mink(4, hedge(13))) * g(mink(4, hedge(15)), mink(4, hedge(16)))
            - Q(7, mink(4, hedge(16))) * g(mink(4, hedge(13)), mink(4, hedge(15)))
            + Q(8, mink(4, hedge(13))) * g(mink(4, hedge(15)), mink(4, hedge(16)))
            - Q(8, mink(4, hedge(15))) * g(mink(4, hedge(13)), mink(4, hedge(16)))
        )
        * Q(0, mink(4, edge(0, 1)))
        * Q(1, mink(4, edge(1, 1)))
        * Q(3, mink(4, edge(3, 1)))
        * Q(4, mink(4, edge(4, 1)))
        * g(mink(4, hedge(4)), mink(4, hedge(5)))
        * g(mink(4, hedge(10)), mink(4, hedge(11)))
        * g(mink(4, hedge(12)), mink(4, hedge(13)))
        * g(mink(4, hedge(14)), mink(4, hedge(15)))
        * g(mink(4, hedge(16)), mink(4, hedge(17)))
        * g(dind(cof(Nc, hedge(1))), cof(Nc, hedge(0)))
        * g(dind(cof(Nc, hedge(2))), cof(Nc, hedge(3)))
        * g(dind(cof(Nc, hedge(6))), cof(Nc, hedge(7)))
        * g(dind(cof(Nc, hedge(9))), cof(Nc, hedge(8)))
        * g(coad(-1 + Nc^2, hedge(4)), coad(-1 + Nc^2, hedge(5)))
        * g(coad(-1 + Nc^2, hedge(10)), coad(-1 + Nc^2, hedge(11)))
        * g(coad(-1 + Nc^2, hedge(12)), coad(-1 + Nc^2, hedge(13)))
        * g(coad(-1 + Nc^2, hedge(14)), coad(-1 + Nc^2, hedge(15)))
        * g(coad(-1 + Nc^2, hedge(16)), coad(-1 + Nc^2, hedge(17)))
        * gamma(bis(4, hedge(0)), bis(4, hedge(2)), mink(4, hedge(4)))
        * gamma(bis(4, hedge(3)), bis(4, hedge(9)), mink(4, hedge(14)))
        * gamma(bis(4, hedge(7)), bis(4, hedge(1)), mink(4, hedge(12)))
        * gamma(bis(4, hedge(8)), bis(4, hedge(6)), mink(4, hedge(10)))
        * gamma(bis(4, hedge(1)), bis(4, hedge(0)), mink(4, edge(0, 1)))
        * gamma(bis(4, hedge(2)), bis(4, hedge(3)), mink(4, edge(1, 1)))
        * gamma(bis(4, hedge(6)), bis(4, hedge(7)), mink(4, edge(3, 1)))
        * gamma(bis(4, hedge(9)), bis(4, hedge(8)), mink(4, edge(4, 1)))
        * t(coad(-1 + Nc^2, hedge(4)), cof(Nc, hedge(2)), dind(cof(Nc, hedge(0))))
        * t(coad(-1 + Nc^2, hedge(10)), cof(Nc, hedge(6)), dind(cof(Nc, hedge(8))))
        * t(coad(-1 + Nc^2, hedge(12)), cof(Nc, hedge(1)), dind(cof(Nc, hedge(7))))
        * t(coad(-1 + Nc^2, hedge(14)), cof(Nc, hedge(9)), dind(cof(Nc, hedge(3))))
        * f(coad(-1 + Nc^2, hedge(5)), coad(-1 + Nc^2, hedge(11)), coad(-1 + Nc^2, hedge(17)))
        * f(coad(-1 + Nc^2, hedge(13)), coad(-1 + Nc^2, hedge(15)), coad(-1 + Nc^2, hedge(16)))
        ",
        default_namespace = "spenso"
    );
    let color_zero_candidate = input.cook_indices().simplify_color().collect_color();

    assert!(
        color_zero_candidate.is_zero(),
        "Expected color to be zero, got {}",
        color_zero_candidate.printer(SpensoPrintSettings::compact().nice_symbolica())
    );
}

#[test]
fn color_casimir_basis_rewrites_dimensions() {
    test_initialize();
    let expr = parse_lit!(Nc ^ -1 + Nc ^ 2 - 1 + NA, default_namespace = "spenso");

    assert_snapshot!(expr.to_color_casimir().to_bare_ordered_string(), @"-2*CF+4*CA*CF+CA");
}

#[test]
fn color_casimir_basis_keeps_trace_normalization_symbolic_by_default() {
    test_initialize();
    let expr = parse_lit!(TR, default_namespace = "spenso");

    assert_snapshot!(expr.to_color_casimir().to_bare_ordered_string(), @"TR");
    assert_snapshot!(expr
        .to_color_casimir_with(ColorCasimirSettings::default().with_trace_normalization())
        .to_bare_ordered_string(), @"1/2");
}

#[test]
fn antisymmetric_two_generator_trace_vanishes() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        antisym!(color_t!(slot!(r.coad_na, a)), color_t!(slot!(r.coad_na, b)))
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"0");
}

#[test]
fn antisymmetric_three_generator_trace_reduces_to_structure_constant() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        antisym!(
            color_t!(slot!(r.coad_na, a)),
            color_t!(slot!(r.coad_na, b)),
            color_t!(slot!(r.coad_na, c)),
        )
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1𝑖/2*f(coad(NA,a),coad(NA,b),coad(NA,c))*idx(2,cof(Nc))");
}

#[test]
fn antisymmetric_trace_commutator_reduces_before_terminal_trace() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        antisym!(color_t!(slot!(r.coad_na, a)), color_t!(slot!(r.coad_na, b))),
        color_t!(slot!(r.coad_na, c)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1𝑖/2*f(coad(NA,a),coad(NA,b),coad(NA,c))*idx(2,cof(Nc))");
}

#[test]
fn antisymmetric_trace_commutator_preserves_projector_sign() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        antisym!(color_t!(slot!(r.coad_na, b)), color_t!(slot!(r.coad_na, a))),
        color_t!(slot!(r.coad_na, c)),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"-1𝑖/2*f(coad(NA,a),coad(NA,b),coad(NA,c))*idx(2,cof(Nc))");
}

#[test]
fn antisymmetric_chain_commutator_reduces_to_structure_constant() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), j),
        antisym!(color_t!(slot!(r.coad_na, a)), color_t!(slot!(r.coad_na, b))),
    );

    assert_snapshot!(expr.simplify_color().to_bare_ordered_string(), @"1𝑖/2*chain(cof(Nc,i),dind(cof(Nc,j)),t(coad(NA,x),in,out))*f(coad(NA,a),coad(NA,b),coad(NA,x))");
}

#[test]
fn cyclic_trace_structure_product_scans_nonleading_pairs() {
    test_initialize();
    let r = TestReps::new();
    let expr = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, c)),
        color_t!(slot!(r.coad_na, d)),
        color_t!(slot!(r.coad_na, b)),
    ) * color_f!(
        slot!(r.coad_na, a),
        slot!(r.coad_na, d),
        slot!(r.coad_na, c),
    );

    let simplified =
        expr.simplify_color_with(ColorSimplifySettings::default().with_cof_dimension_invariants());
    let simplified = simplified.to_bare_ordered_string();

    assert!(
        !simplified.contains("trace(") && !simplified.contains("t(") && !simplified.contains("f("),
        "expected closed color structure to reduce to scalar invariants, got {simplified}"
    );
}

#[test]
fn trace_structure_product_orientation_sign_cancels() {
    test_initialize();
    let r = TestReps::new();
    let expr = color_f!(
        slot!(r.coad_na, q),
        slot!(r.coad_na, a),
        slot!(r.coad_na, b),
    ) * (trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, c)),
    ) + trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, c)),
    ));

    assert_color_zero(expr);
}

#[test]
fn expand_color_keeps_closed_traces_on_color_side() {
    test_initialize();
    let r = TestReps::new();
    let trace_factor = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
        color_t!(slot!(r.coad_na, c)),
    );
    let expr = Atom::var(s!(x))
        * trace_factor
        * color_f!(
            slot!(r.coad_na, a),
            slot!(r.coad_na, b),
            slot!(r.coad_na, c),
        );

    let expanded = expr.expand_color();

    assert_eq!(expanded.len(), 1);
    let (color_factor, residual) = &expanded[0];
    let color_factor = color_factor.to_bare_ordered_string();
    assert!(
        color_factor.contains("trace(") && color_factor.contains("f("),
        "expected closed trace and structure constant to stay in color factor, got {color_factor}"
    );
    let residual = residual.to_bare_ordered_string();
    assert_eq!(
        residual, "x",
        "expected residual factor to be color-free, got {}",
        residual,
    );
}

#[test]
fn color_simplify_defaults_match_simplify_color() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), j),
        color_t!(slot!(r.coad_na, a)),
    ) * chain!(
        slot!(r.cof_nc, k),
        slot!(r.cof_nc.dual(), l),
        color_t!(slot!(r.coad_na, a)),
    );

    assert_eq!(
        expr.simplify_color_with(ColorSimplifySettings::default()),
        expr.simplify_color()
    );
}

#[test]
fn color_trace_evaluation_can_be_disabled() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), i),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
    );
    let expected = trace!(
        &r.cof_nc,
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
    );

    assert_eq!(
        expr.simplify_color_with(ColorSimplifySettings::default().without_trace_evaluation()),
        expected,
    );
}

#[test]
fn color_cross_chain_fierz_can_be_disabled() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), j),
        color_t!(slot!(r.coad_na, a)),
    ) * chain!(
        slot!(r.cof_nc, k),
        slot!(r.cof_nc.dual(), l),
        color_t!(slot!(r.coad_na, a)),
    );

    assert_eq!(
        expr.simplify_color_with(
            ColorSimplifySettings::default().without_cross_chain_fierz_expansion()
        ),
        expr,
    );
}

#[test]
fn color_cross_chain_fierz_handles_longer_open_chains() {
    test_initialize();
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.cof_nc, i),
        slot!(r.cof_nc.dual(), j),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, b)),
    ) * chain!(
        slot!(r.cof_nc, k),
        slot!(r.cof_nc.dual(), l),
        color_t!(slot!(r.coad_na, a)),
        color_t!(slot!(r.coad_na, c)),
    );
    let index = color_idx!(2, ColorFundamental {}.to_symbolic([Atom::var(s!(Nc))]));
    let expected = index.clone()
        * chain!(
            slot!(r.cof_nc, i),
            slot!(r.cof_nc.dual(), l),
            color_t!(slot!(r.coad_na, c)),
        )
        * chain!(
            slot!(r.cof_nc, k),
            slot!(r.cof_nc.dual(), j),
            color_t!(slot!(r.coad_na, b)),
        )
        - index
            * chain!(
                slot!(r.cof_nc, i),
                slot!(r.cof_nc.dual(), j),
                color_t!(slot!(r.coad_na, b)),
            )
            * chain!(
                slot!(r.cof_nc, k),
                slot!(r.cof_nc.dual(), l),
                color_t!(slot!(r.coad_na, c)),
            )
            / s!(Nc);

    assert_eq!(
        expr.simplify_color().simplify_metrics(),
        expected.simplify_metrics(),
    );
}

#[test]
fn symmetric_trace_d33_partial_contraction() {
    test_initialize();
    let r = TestReps::new();
    let left = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_na, a)),
            color_t!(slot!(r.coad_na, b)),
            color_t!(slot!(r.coad_na, c)),
        )
    );
    let right = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_na, a)),
            color_t!(slot!(r.coad_na, b)),
            color_t!(slot!(r.coad_na, d)),
        )
    );

    assert_snapshot!((left * right).simplify_color().to_bare_ordered_string(), @"NA^(-1)*g(coad(NA,c),coad(NA,d))*gram(3,cof(Nc),cof(Nc))");
}

#[test]
fn symmetric_trace_d44_partial_contraction() {
    test_initialize();
    let r = TestReps::new();
    let left = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_na, a)),
            color_t!(slot!(r.coad_na, b)),
            color_t!(slot!(r.coad_na, c)),
            color_t!(slot!(r.coad_na, d)),
        )
    );
    let right = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_na, a)),
            color_t!(slot!(r.coad_na, b)),
            color_t!(slot!(r.coad_na, c)),
            color_t!(slot!(r.coad_na, e)),
        )
    );

    assert_snapshot!((left * right).simplify_color().to_bare_ordered_string(), @"NA^(-1)*g(coad(NA,d),coad(NA,e))*gram(4,cof(Nc),cof(Nc))");
}

#[test]
fn symmetric_trace_d44_full_contraction() {
    test_initialize();
    let r = TestReps::new();
    let left = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_na, a)),
            color_t!(slot!(r.coad_na, b)),
            color_t!(slot!(r.coad_na, c)),
            color_t!(slot!(r.coad_na, d)),
        )
    );
    let right = trace!(
        &r.cof_nc,
        sym!(
            color_t!(slot!(r.coad_na, a)),
            color_t!(slot!(r.coad_na, b)),
            color_t!(slot!(r.coad_na, c)),
            color_t!(slot!(r.coad_na, d)),
        )
    );

    assert_snapshot!((left * right).simplify_color().to_bare_ordered_string(), @"gram(4,cof(Nc),cof(Nc))");
}

mod feyncalc_reference;
mod form_reference;

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
    test_initialize();

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
    test_initialize();
    println!("{}", CS.t_strct::<AbstractIndex>(3, 8));

    let _ = Atom::Zero.simplify_metrics();
}

#[test]
fn test_val() {
    test_initialize();
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

    let simplified = expr.cook_indices().simplify_color();

    assert_snapshot!(simplified.collect_color_constants().collect_factors().to_bare_ordered_string(), @"-1𝑖/2*G^4*cas(2,coad(ohoho))*ee^2*idx(2,cof(ahaha))*ohoho");
}

#[test]
fn structure_pair_with_closed_generator_chain_matches_form() {
    test_initialize();
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

    assert_snapshot!(
        expr.simplify_color().to_bare_ordered_string(),
        @"(-1+Nc^2)*-1*cas(2,coad(-1+Nc^2))*idx(2,cof(Nc))"
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
    let residual = (expr1.simplify_color() + expr2.simplify_color())
        .expand()
        .simplify_metrics()
        .cook_indices();

    println!("{}", residual);

    let residual = residual
        .canonize(AbstractIndex::Dummy)
        .simplify_metrics()
        .expand();
    assert!(residual.is_zero());
}

mod failing {
    use super::*;

    #[test]
    fn test_color_matrix_element() {
        test_initialize();
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

        assert_ne!(tgt, simplified_amp_squared.factor());
    }

    #[test]
    fn test_color_matrix_element_two() {
        test_initialize();

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

        assert_ne!(tgt, simplified_amp_squared.factor());
    }
}
