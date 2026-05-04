use spenso::network::StructureLessDisplay;
use spenso::network::parsing::ParseSettings;
use spenso::network::store::TensorScalarStore;
use spenso::shadowing::symbolica_utils::AtomCoreExt;
use spenso::shadowing::symbolica_utils::TypstSettings;
use spenso::structure::IndexlessNamedStructure;
use spenso::structure::PermutedStructure;
use spenso::symbolica_atom::IntoAtom;
use spenso::{chain, s, slot};

use crate::gamma;

static GG: LazyLock<PermutedStructure<IndexlessNamedStructure<Symbol, ()>>> = LazyLock::new(|| {
    IndexlessNamedStructure::from_iter(
        [
            Bispinor {}.new_rep(4).to_lib(),
            Bispinor {}.new_rep(4).cast(),
            Minkowski {}.new_rep(4).cast(),
        ],
        AGS.gamma,
        None,
    )
});

use super::*;

use crate::color::ColorSimplifier;
use crate::tensor::SymbolicNetParse;
use crate::tensor::SymbolicTensor;
use crate::test::test_initialize;
use spenso::structure::{abstract_index::AbstractIndex, permuted::Perm};
use symbolica::{
    atom::{Atom, AtomCore},
    parse_lit,
};

use crate::representations::initialize;
use crate::test_support::{TestReps, assert_bare_snapshot};

fn assert_gamma_zero(expr: Atom) {
    assert!(expr.simplify_gamma().is_zero());
}

#[test]
fn gamma_construct() {
    println!("{}", AGS.gamma_pattern(RS.a_, RS.b_, RS.c_));

    let f = GG
        .clone()
        .reindex([4, 3, 2])
        .unwrap()
        .map_structure(|a| SymbolicTensor::from_named(&a).unwrap());

    let f_s = f.structure.structure.clone();

    // f.rep_permutation = f.rep_permutation.inverse();

    let f_p = f.permute_reps_wrapped().permute_inds();

    println!(
        "Structure:{}\nPermuted:{}\nPermuted Structure{}\nMetric simplified{}",
        f_s,
        f_p,
        f_p.structure,
        f_p.expression.simplify_metrics()
    );
}

#[test]
fn chain_test() {
    initialize();

    let expr = parse_lit!(
        (-1 * P(2, mink(dim, l(20))) + P(1, mink(dim, l(20)))) * 1𝑖 * G
            ^ 2 * g(bis(4, l(2)), bis(4, l(4)))
                * g(bis(4, l(3)), bis(4, l(7)))
                * g(mink(dim, l(0)), mink(dim, l(5)))
                * g(mink(dim, l(1)), mink(dim, l(4)))
                * gamma(bis(4, l(5)), bis(4, l(4)), mink(dim, l(4)))
                * gamma(bis(4, l(6)), bis(4, l(5)), mink(dim, l(20)))
                * gamma(bis(4, l(7)), bis(4, l(6)), mink(dim, l(5))),
        default_namespace = "spenso"
    );

    println!("Bef:{}", expr);
    let mut out = String::new();
    expr.typst_fmt(&mut out, &TypstSettings::lowering())
        .unwrap();
    println!("{}", out);
    println!("Aft:{}", expr.simplify_gamma());
}

#[test]
fn normalise_g() {
    let expr = parse_lit!(
        gamma_chain(mink(dim, mu), bis(5), mink(dim, mu), bis(1), bis(2)),
        default_namespace = "spenso"
    );

    println!("{}", expr.simplify_gamma())
}

#[test]
fn form_two_gamma_trace() {
    initialize();
    let expr = gamma!(mu, a, b) * gamma!(nu, b, a);

    assert_bare_snapshot!(expr.simplify_gamma(),@"4*g(mink(4,mu),mink(4,nu))");
}

#[test]
fn gamma_macro_accepts_integer_indices() {
    initialize();
    let expr = gamma!(1, 2, 3);

    assert_bare_snapshot!(expr, @"gamma(bis(4,2),bis(4,3),mink(4,1))");
}

#[test]
fn gamma_macro_accepts_mixed_default_and_explicit_indices() {
    let r = TestReps::new();
    let expr = gamma!(mu, slot!(r.bis_d, a), 1);

    assert_bare_snapshot!(expr, @"gamma(bis(d,a),bis(4,1),mink(4,mu))");
}

#[test]
fn form_odd_gamma_trace_vanishes() {
    initialize();
    let expr = gamma!(mu, a, b) * gamma!(nu, b, c) * gamma!(rho, c, a);

    assert_gamma_zero(expr);
}

#[test]
fn form_four_gamma_trace_recurses() {
    initialize();
    let expr = gamma!(mu, a, b) * gamma!(nu, b, c) * gamma!(rho, c, d) * gamma!(sigma, d, a);

    assert_bare_snapshot!(expr.simplify_gamma().expand(),@"-4*g(mink(4,mu),mink(4,rho))*g(mink(4,nu),mink(4,sigma))+4*g(mink(4,mu),mink(4,nu))*g(mink(4,rho),mink(4,sigma))+4*g(mink(4,mu),mink(4,sigma))*g(mink(4,nu),mink(4,rho))"
    );
}

#[test]
fn form_repeated_lorentz_gamma_chain_contracts_to_dimension() {
    initialize();
    let expr = gamma!(mu, a, b) * gamma!(mu, b, c);

    assert_bare_snapshot!(expr.simplify_gamma(),@"4*g(bis(4,a),bis(4,c))");
}

#[test]
#[ignore = "pending public chain-based gamma contraction"]
fn form_adjacent_chain_lorentz_contraction() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.bis_d, a),
        slot!(r.bis_d, b),
        gamma!(slot!(r.mink_d, mu)),
        gamma!(slot!(r.mink_d, mu)),
    );

    assert_bare_snapshot!(expr.simplify_gamma(), @"d*g(bis(d,a),bis(d,b))");
}

#[test]
#[ignore = "pending public chain-based Chisholm reduction"]
fn form_chisholm_odd_interior_chain() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.bis4, a),
        slot!(r.bis4, b),
        gamma!(slot!(r.mink4, mu)),
        gamma!(slot!(r.mink4, nu1)),
        gamma!(slot!(r.mink4, nu2)),
        gamma!(slot!(r.mink4, nu3)),
        gamma!(slot!(r.mink4, mu)),
    );

    assert_bare_snapshot!(expr.simplify_gamma(),@"-2*chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,nu3)),gamma(in,out,mink(4,nu2)),gamma(in,out,mink(4,nu1)))"
    );
}

#[test]
#[ignore = "pending public chain-based Chisholm reduction"]
fn form_chisholm_two_interior_chain() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.bis4, a),
        slot!(r.bis4, b),
        gamma!(slot!(r.mink4, mu)),
        gamma!(slot!(r.mink4, nu)),
        gamma!(slot!(r.mink4, rho)),
        gamma!(slot!(r.mink4, mu)),
    );

    assert_bare_snapshot!(expr.simplify_gamma(),@"4*g(bis(4,a),bis(4,b))*g(mink(4,nu),mink(4,rho))"
    );
}

#[test]
#[ignore = "pending public chain-based gamma-five epsilon reduction"]
fn form_gamma_five_epsilon_trick() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.bis4, a),
        slot!(r.bis4, b),
        gamma!(slot!(r.mink4, mu)),
        gamma!(slot!(r.mink4, nu)),
        gamma!(slot!(r.mink4, rho)),
    );

    assert_bare_snapshot!(expr.simplify_gamma(),@"chain(bis(4,a),bis(4,b),gamma5(in,out),gamma(in,out,mink(4,sigma)))*epsilon(mink(4,mu),mink(4,nu),mink(4,rho),mink(4,sigma))+chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,rho)))*g(mink(4,mu),mink(4,nu))+-1*chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,nu)))*g(mink(4,mu),mink(4,rho))+chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,mu)))*g(mink(4,nu),mink(4,rho))"
    );
}

#[test]
#[ignore = "pending gamma-five trace term-count stress support"]
fn form_gamma_five_symmetric_12_term_count() {
    initialize();
    let term_count = 0usize;

    insta::assert_snapshot!(term_count.to_string(), @"51975");
}

#[test]
#[ignore = "pending gamma-five trace4 term-count stress support"]
fn form_gamma_five_regular_12_term_count() {
    initialize();
    let term_count = 0usize;

    insta::assert_snapshot!(term_count.to_string(), @"1029");
}

#[test]
fn collect_expand_chain() {
    initialize();
    let expr = parse_lit!(gamma_chain(
        mink(dim, nu1),
        bis(dim, 3),
        mink(dim, nu12),
        bis(dim, 31),
        mink(dim, nu13),
        bis(dim, 32),
        mink(dim, nu14),
        bis(dim, 33),
        mink(dim, nu),
        bis(dim, 34),
        mink(dim, nu16),
        bis(dim, 35),
        mink(dim, nu17),
        bis(dim, 36),
        mink(dim, nu),
        bis(dim, 4),
        mink(dim, nu3),
        bis(dim, 5),
        mink(dim, nu2),
        bis(dim, 1),
        bis(dim, 2)
    ));

    let mut a = expr.clone();
    collect_gammas(&mut a);
    undo_gamma_chain(&mut a);
    collect_gammas(&mut a);

    assert_eq!(a, expr);
    // println!("{a}");
}

#[test]
fn gl23() {
    test_initialize();
    let expr = parse_lit!(
        (-1 * Q(2, mink(4, hedge(8))) * g(mink(4, hedge(12)), mink(4, hedge(4)))
            + -1 * Q(4, mink(4, hedge(12))) * g(mink(4, hedge(4)), mink(4, hedge(8)))
            + -1 * Q(6, mink(4, hedge(4))) * g(mink(4, hedge(12)), mink(4, hedge(8)))
            + Q(2, mink(4, hedge(12))) * g(mink(4, hedge(4)), mink(4, hedge(8)))
            + Q(4, mink(4, hedge(4))) * g(mink(4, hedge(12)), mink(4, hedge(8)))
            + Q(6, mink(4, hedge(8))) * g(mink(4, hedge(12)), mink(4, hedge(4))))
            * 1𝑖
            / 9
            * G
            ^ 4 * Q(1, mink(4, edge(1, 1)))
                * Q(3, mink(4, edge(3, 1)))
                * Q(5, mink(4, edge(5, 1)))
                * Q(7, mink(4, edge(7, 1)))
                * Q(8, mink(4, edge(8, 1)))
                * ee
            ^ 2 * f(coad(8, hedge(4)), coad(8, hedge(8)), coad(8, hedge(12)))
                * g(coad(8, hedge(11)), coad(8, hedge(12)))
                * g(coad(8, hedge(3)), coad(8, hedge(4)))
                * g(coad(8, hedge(7)), coad(8, hedge(8)))
                * g(cof(3, hedge(1)), dind(cof(3, hedge(2))))
                * g(cof(3, hedge(10)), dind(cof(3, hedge(14))))
                * g(cof(3, hedge(14)), dind(cof(3, hedge(13))))
                * g(cof(3, hedge(15)), dind(cof(3, hedge(16))))
                * g(cof(3, hedge(16)), dind(cof(3, hedge(6))))
                * g(cof(3, hedge(6)), dind(cof(3, hedge(5))))
                * g(cof(3, hedge(9)), dind(cof(3, hedge(10))))
                * g(mink(4, hedge(11)), mink(4, hedge(12)))
                * g(mink(4, hedge(3)), mink(4, hedge(4)))
                * g(mink(4, hedge(7)), mink(4, hedge(8)))
                * gamma(bis(4, hedge(1)), bis(4, hedge(5)), mink(4, hedge(3)))
                * gamma(bis(4, hedge(10)), bis(4, hedge(9)), mink(4, edge(5, 1)))
                * gamma(bis(4, hedge(13)), bis(4, hedge(14)), mink(4, edge(7, 1)))
                * gamma(bis(4, hedge(14)), bis(4, hedge(10)), mink(4, hedge(17)))
                * gamma(bis(4, hedge(15)), bis(4, hedge(13)), mink(4, hedge(11)))
                * gamma(bis(4, hedge(16)), bis(4, hedge(15)), mink(4, edge(8, 1)))
                * gamma(bis(4, hedge(2)), bis(4, hedge(1)), mink(4, edge(1, 1)))
                * gamma(bis(4, hedge(5)), bis(4, hedge(6)), mink(4, edge(3, 1)))
                * gamma(bis(4, hedge(6)), bis(4, hedge(16)), mink(4, hedge(0)))
                * gamma(bis(4, hedge(9)), bis(4, hedge(2)), mink(4, hedge(7)))
                * t(
                    coad(8, hedge(11)),
                    cof(3, hedge(13)),
                    dind(cof(3, hedge(15)))
                )
                * t(coad(8, hedge(3)), cof(3, hedge(5)), dind(cof(3, hedge(1))))
                * t(coad(8, hedge(7)), cof(3, hedge(2)), dind(cof(3, hedge(9)))),
        default_namespace = "spenso"
    );

    println!(
        "{}\n",
        expr.simplify_metrics()
            .cook_indices()
            .canonize(AbstractIndex::Dummy)
    );

    println!(
        "Colored done: {}\n",
        expr.simplify_metrics()
            .cook_indices()
            .canonize(AbstractIndex::Dummy)
            .simplify_color()
    );
}

#[test]
fn gl24() {
    test_initialize();
    let expr = parse_lit!(
        (-1 * Q(2, mink(4, hedge(8))) * g(mink(4, hedge(12)), mink(4, hedge(4)))
            + -1 * Q(4, mink(4, hedge(12))) * g(mink(4, hedge(4)), mink(4, hedge(8)))
            + -1 * Q(6, mink(4, hedge(4))) * g(mink(4, hedge(12)), mink(4, hedge(8)))
            + Q(2, mink(4, hedge(12))) * g(mink(4, hedge(4)), mink(4, hedge(8)))
            + Q(4, mink(4, hedge(4))) * g(mink(4, hedge(12)), mink(4, hedge(8)))
            + Q(6, mink(4, hedge(8))) * g(mink(4, hedge(12)), mink(4, hedge(4))))
            * 1𝑖
            / 9
            * G
            ^ 4 * Q(1, mink(4, edge(1, 1)))
                * Q(3, mink(4, edge(3, 1)))
                * Q(5, mink(4, edge(5, 1)))
                * Q(7, mink(4, edge(7, 1)))
                * Q(8, mink(4, edge(8, 1)))
                * ee
            ^ 2 * f(coad(8, hedge(4)), coad(8, hedge(8)), coad(8, hedge(12)))
                * g(coad(8, hedge(11)), coad(8, hedge(12)))
                * g(coad(8, hedge(3)), coad(8, hedge(4)))
                * g(coad(8, hedge(7)), coad(8, hedge(8)))
                * g(cof(3, hedge(10)), dind(cof(3, hedge(9))))
                * g(cof(3, hedge(13)), dind(cof(3, hedge(14))))
                * g(cof(3, hedge(14)), dind(cof(3, hedge(10))))
                * g(cof(3, hedge(16)), dind(cof(3, hedge(15))))
                * g(cof(3, hedge(2)), dind(cof(3, hedge(1))))
                * g(cof(3, hedge(5)), dind(cof(3, hedge(6))))
                * g(cof(3, hedge(6)), dind(cof(3, hedge(16))))
                * g(mink(4, hedge(11)), mink(4, hedge(12)))
                * g(mink(4, hedge(3)), mink(4, hedge(4)))
                * g(mink(4, hedge(7)), mink(4, hedge(8)))
                * gamma(bis(4, hedge(1)), bis(4, hedge(2)), mink(4, edge(1, 1)))
                * gamma(bis(4, hedge(10)), bis(4, hedge(14)), mink(4, hedge(17)))
                * gamma(bis(4, hedge(13)), bis(4, hedge(15)), mink(4, hedge(11)))
                * gamma(bis(4, hedge(14)), bis(4, hedge(13)), mink(4, edge(7, 1)))
                * gamma(bis(4, hedge(15)), bis(4, hedge(16)), mink(4, edge(8, 1)))
                * gamma(bis(4, hedge(16)), bis(4, hedge(6)), mink(4, hedge(0)))
                * gamma(bis(4, hedge(2)), bis(4, hedge(9)), mink(4, hedge(7)))
                * gamma(bis(4, hedge(5)), bis(4, hedge(1)), mink(4, hedge(3)))
                * gamma(bis(4, hedge(6)), bis(4, hedge(5)), mink(4, edge(3, 1)))
                * gamma(bis(4, hedge(9)), bis(4, hedge(10)), mink(4, edge(5, 1)))
                * t(
                    coad(8, hedge(11)),
                    cof(3, hedge(15)),
                    dind(cof(3, hedge(13)))
                )
                * t(coad(8, hedge(3)), cof(3, hedge(1)), dind(cof(3, hedge(5))))
                * t(coad(8, hedge(7)), cof(3, hedge(9)), dind(cof(3, hedge(2))))
                * ϵ(0, mink(4, hedge(0)))
                * ϵbar(0, mink(4, hedge(17))),
        default_namespace = "spenso"
    );

    println!(
        "{}\n",
        expr.simplify_metrics()
            .cook_indices()
            .canonize(AbstractIndex::Dummy)
    );

    println!(
        "Colored done: {}\n",
        expr.simplify_metrics()
            .cook_indices()
            .canonize(AbstractIndex::Dummy)
            .simplify_color()
    );
}

#[test]
fn gl_06() {
    test_initialize();
    let expr = parse_lit!(
        1 / 6 * Nc * ee
            ^ 4 * sw
            ^ -2 * vev
                * I3x21
                * (MC * g(bis(4, hedge(1)), bis(4, hedge(2)))
                    - K(0, mink(4, edge(1, 1)))
                        * gamma(bis(4, hedge(1)), bis(4, hedge(2)), mink(4, edge(1, 1))))
                * (-K(0, mink(4, edge(3, 1))) - K(1, mink(4, edge(3, 1))))
                * (-g(mink(4, hedge(7)), mink(4, hedge(8))) + MW
                    ^ -2 * (-P(0, mink(4, hedge(7))) - K(1, mink(4, hedge(7))))
                        * (-P(0, mink(4, hedge(8))) - K(1, mink(4, hedge(8)))))
                * (P(0, mink(4, edge(5, 1)))
                    + K(0, mink(4, edge(5, 1)))
                    + K(1, mink(4, edge(5, 1))))
                * conj(CKM2x1)
                * ϵ(0, mink(4, hedge(0)))
                * ϵbar(0, mink(4, hedge(11)))
                * g(mink(4, hedge(0)), mink(4, hedge(8)))
                * gamma(bis(4, hedge(10)), bis(4, hedge(6)), mink(4, hedge(11)))
                * gamma(bis(4, hedge(2)), bis(4, vertex(1, 1)), mink(4, hedge(7)))
                * gamma(bis(4, hedge(6)), bis(4, hedge(5)), mink(4, edge(3, 1)))
                * gamma(bis(4, hedge(9)), bis(4, hedge(10)), mink(4, edge(5, 1)))
                * projm(bis(4, hedge(5)), bis(4, hedge(1)))
                * projm(bis(4, vertex(1, 1)), bis(4, hedge(9)))
                * (1 / 2)
            ^ (1 / 2),
        default_namespace = "spenso"
    );

    let simplified = expr
        .cook_indices()
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
            depth_limit: Some(0),
            ..Default::default()
        })
        .unwrap();

    println!(
        "{}",
        simplified.graph.dot_impl(
            |i| {
                let ss = &simplified.store.get_scalar(i);
                format!("{}:{}", i, ss)
            },
            |k| k.display(),
            |t| {
                let tt = &simplified.store.get_tensor(t);
                format!("T{}:{}", t, tt.expression.to_bare_ordered_string())
            },
            |fk| fk.to_string(),
        )
    )

    // println!("{}", expr.cook_indices().simplify_gamma().simplify_gamma());
}

mod failing {
    use super::*;

    #[test]
    fn gamma_alg() {
        let r = TestReps::new();
        let mink4 = r.mink4;
        let mink_dim = r.mink_d;
        let bis_dim = r.bis_d;

        fn p(index: impl IntoAtom) -> Atom {
            function!(symbol!("spenso::p"), index.into_atom())
        }

        fn q(index: impl IntoAtom) -> Atom {
            function!(symbol!("spenso::q"), index.into_atom())
        }

        let expr = (gamma!(0, 1, 3) * gamma!(0, 3, 2)).simplify_gamma();

        assert_bare_snapshot!(expr, @"4*g(bis(4,1),bis(4,2))");

        let expr = (p(slot!(mink4, nu1))
            * (p(slot!(mink4, nu3)) + q(slot!(mink4, nu3)))
            * gamma!(nu1, 1, 3)
            * gamma!(mu, 3, 4)
            * gamma!(nu3, 4, 5)
            * gamma!(nu, 5, 1))
        .simplify_gamma()
        .expand()
        .replace(s!(nu3))
        .with(s!(dummy))
        .replace(s!(nu1))
        .with(s!(dummy));

        assert_bare_snapshot!(expr, @"(p(mink(4,dummy)))^2*-4*g(mink(4,mu),mink(4,nu))+-4*g(mink(4,mu),mink(4,nu))*p(mink(4,dummy))*q(mink(4,dummy))+4*p(mink(4,mu))*q(mink(4,nu))+4*p(mink(4,nu))*q(mink(4,mu))+8*p(mink(4,mu))*p(mink(4,nu))");

        let expr = (mink_dim.g(5, 6)
            * (mink_dim.g(1, 2) * mink_dim.g(3, 4) * mink_dim.g(5, 6)
                - mink_dim.g(1, 3) * mink_dim.g(2, 6) * mink_dim.g(5, 4))
            * (mink_dim.g(1, 2) * mink_dim.g(3, 4) - mink_dim.g(1, 3) * mink_dim.g(2, 4)))
        .simplify_gamma();

        assert_bare_snapshot!(expr, @"-1*d+d^3");

        let expr = (p(slot!(mink4, nu1))
            * (p(slot!(mink4, nu3)) + q(slot!(mink4, nu3)))
            * gamma!(nu1, 1, 3)
            * gamma!(mu, 3, 4)
            * gamma!(nu, 4, 5)
            * gamma!(nu3, 5, 1))
        .simplify_gamma()
        .replace(s!(nu3))
        .with(s!(dummy))
        .replace(s!(nu1))
        .with(s!(dummy));

        assert_bare_snapshot!(expr, @"(p(mink(4,dummy)))^2*4*g(mink(4,mu),mink(4,nu))+-4*p(mink(4,nu))*q(mink(4,mu))+4*g(mink(4,mu),mink(4,nu))*p(mink(4,dummy))*q(mink(4,dummy))+4*p(mink(4,mu))*q(mink(4,nu))");

        let expr = (p(slot!(mink_dim, nu1))
            * (p(slot!(mink_dim, nu3)) + q(slot!(mink_dim, nu3)))
            * gamma!(slot!(mink_dim, nu1), slot!(bis_dim, 1), slot!(bis_dim, 3))
            * gamma!(slot!(mink_dim, nu), slot!(bis_dim, 3), slot!(bis_dim, 4))
            * gamma!(slot!(mink_dim, nu), slot!(bis_dim, 4), slot!(bis_dim, 5))
            * gamma!(slot!(mink_dim, nu3), slot!(bis_dim, 5), slot!(bis_dim, 1)))
        .simplify_gamma()
        .replace(s!(nu1))
        .with(s!(nu3));

        assert_bare_snapshot!(
            expr.expand().canonize(AbstractIndex::Dummy),
            @"(p(mink(d,d_0)))^2*4*d+4*d*p(mink(d,d_0))*q(mink(d,d_0))"
        );

        let expr = (p(slot!(mink_dim, nu1))
            * (p(slot!(mink_dim, nu3)) + q(slot!(mink_dim, nu3)))
            * gamma!(slot!(mink_dim, nu1), slot!(bis_dim, 1), slot!(bis_dim, 3))
            * gamma!(slot!(mink_dim, nu), slot!(bis_dim, 3), slot!(bis_dim, 4))
            * gamma!(slot!(mink_dim, nu3), slot!(bis_dim, 4), slot!(bis_dim, 5))
            * gamma!(slot!(mink_dim, nu), slot!(bis_dim, 5), slot!(bis_dim, 1)))
        .simplify_gamma();

        assert_bare_snapshot!(
            expr.expand().canonize(AbstractIndex::Dummy),
            @"(p(mink(d,d_0)))^2*-4*d+(p(mink(d,d_0)))^2*8+-4*d*p(mink(d,d_0))*q(mink(d,d_0))+8*p(mink(d,d_0))*q(mink(d,d_0))"
        );

        let expr = (p(slot!(mink_dim, nu1))
            * q(slot!(mink_dim, nu2))
            * (p(slot!(mink_dim, nu3)) + q(slot!(mink_dim, nu3)))
            * q(slot!(mink_dim, nu4))
            * gamma!(slot!(mink_dim, nu1), slot!(bis_dim, 1), slot!(bis_dim, 3))
            * gamma!(slot!(mink_dim, nu4), slot!(bis_dim, 3), slot!(bis_dim, 4))
            * gamma!(slot!(mink_dim, nu3), slot!(bis_dim, 4), slot!(bis_dim, 5))
            * gamma!(slot!(mink_dim, nu2), slot!(bis_dim, 5), slot!(bis_dim, 1)))
        .simplify_gamma()
        .to_dots();

        assert_bare_snapshot!(
            expr,
            @"(g(mink(d),p,q))^2*8+-4*g(mink(d),p,p)*g(mink(d),q,q)+4*g(mink(d),p,q)*g(mink(d),q,q)"
        );

        let expr = (gamma!(slot!(mink_dim, mu), 1, 3)
            * gamma!(slot!(mink_dim, nu), 3, 4)
            * gamma!(slot!(mink_dim, mu), 4, 5)
            * gamma!(slot!(mink_dim, nu), 5, 2))
        .simplify_gamma()
        .to_dots();

        assert_bare_snapshot!(expr, @"-1*d^2*g(bis(4,1),bis(4,2))+2*d*g(bis(4,1),bis(4,2))");
    }

    #[test]
    fn val_test() {
        initialize();
        // let expr = parse_lit!(
        //     (MB * g(bis(4, hedge(0, 0)), bis(4, hedge(1, 0)))
        //         + gamma(
        //             bis(4, hedge(0, 0)),
        //             bis(4, hedge(1, 0)),
        //             mink(4, edge(0, 1))
        //         ) * Q(0, mink(4, edge(0, 1))))
        //         * (MB * g(bis(4, hedge(4, 0)), bis(4, hedge(5, 0)))
        //             + gamma(
        //                 bis(4, hedge(4, 0)),
        //                 bis(4, hedge(5, 0)),
        //                 mink(4, edge(2, 1))
        //             ) * Q(2, mink(4, edge(2, 1))))
        //         * (MB * g(bis(4, hedge(8, 0)), bis(4, hedge(9, 0)))
        //             + gamma(
        //                 bis(4, hedge(8, 0)),
        //                 bis(4, hedge(9, 0)),
        //                 mink(4, edge(5, 1))
        //             ) * Q(5, mink(4, edge(5, 1))))
        //         * gamma(
        //             bis(4, hedge(1, 0)),
        //             bis(4, hedge(4, 0)),
        //             mink(4, hedge(10, 0))
        //         )
        //         * gamma(
        //             bis(4, hedge(5, 0)),
        //             bis(4, hedge(8, 0)),
        //             mink(4, hedge(2, 0))
        //         )
        //         * gamma(
        //             bis(4, hedge(9, 0)),
        //             bis(4, hedge(11, 0)),
        //             mink(4, hedge(7, 0))
        //         )
        //         * gamma(
        //             bis(4, hedge(11, 0)),
        //             bis(4, hedge(0, 0)),
        //             mink(4, hedge(2, 0))
        //         )
        //         * p(1, mink(4, hedge(10, 0)))
        //         * p(7, mink(4, hedge(7, 0))),
        //     "spenso"
        // );

        let _expr = parse_lit!(
            (MB * g(bis(4, hedge(0, 0)), bis(4, hedge(1, 0)))
                + gamma(
                    bis(4, hedge(0, 0)),
                    bis(4, hedge(1, 0)),
                    mink(4, edge(0, 1))
                ) * Q(0, mink(4, edge(0, 1))))
                * (MB * g(bis(4, hedge(2, 0)), bis(4, hedge(3, 0)))
                    + gamma(
                        bis(4, hedge(2, 0)),
                        bis(4, hedge(3, 0)),
                        mink(4, edge(1, 1))
                    ) * Q(1, mink(4, edge(1, 1))))
                * (MB * g(bis(4, hedge(5, 0)), bis(4, hedge(6, 0)))
                    + gamma(
                        bis(4, hedge(5, 0)),
                        bis(4, hedge(6, 0)),
                        mink(4, edge(3, 1))
                    ) * Q(3, mink(4, edge(3, 1))))
                * (gamma(
                    bis(4, hedge(9, 0)),
                    bis(4, hedge(10, 0)),
                    mink(4, edge(5, 1))
                ) * Q(5, mink(4, edge(5, 1))))
                * gamma(
                    bis(4, hedge(1, 0)),
                    bis(4, hedge(9, 0)),
                    mink(4, hedge(7, 0))
                )
                * gamma(
                    bis(4, hedge(3, 0)),
                    bis(4, hedge(5, 0)),
                    mink(4, hedge(7, 0))
                )
                * gamma(
                    bis(4, hedge(6, 0)),
                    bis(4, hedge(0, 0)),
                    mink(4, hedge(11, 0))
                )
                * gamma(
                    bis(4, hedge(10, 0)),
                    bis(4, hedge(2, 0)),
                    mink(4, hedge(4, 0))
                )
                * p(1, mink(4, hedge(4, 0)))
                * p(7, mink(4, hedge(11, 0))),
            default_namespace = "spenso"
        );

        // let expr = parse_lit!(
        //     g(mink(D, left(2)), mink(D, left(5)))
        //         * g(mink(D, left(2)), mink(D, right(2)))
        //         * g(mink(D, left(3)), mink(D, left(6)))
        //         * g(mink(D, left(3)), mink(D, right(3)))
        //         * g(mink(D, left(4)), mink(D, left(7)))
        //         * g(mink(D, left(5)), mink(D, left(6)))
        //         * g(mink(D, right(2)), mink(D, right(5)))
        //         * g(mink(D, right(3)), mink(D, right(6)))
        //         * g(mink(D, right(4)), mink(D, right(7)))
        //         * g(mink(D, right(6)), mink(D, right(7)))
        //         * g(bis(D, left(0)), bis(D, left(5)))
        //         * g(bis(D, left(1)), bis(D, left(4)))
        //         * g(bis(D, right(0)), bis(D, right(5)))
        //         * g(bis(D, right(1)), bis(D, right(4)))
        //         * gamma(bis(D, left(1)), bis(D, right(1)), mink(D, 1337))
        //         * gamma(bis(D, right(0)), bis(D, left(0)), mink(D, 1338))
        //         * gamma(bis(D, left(5)), bis(D, left(4)), mink(D, left(4)))
        //         * gamma(bis(D, right(4)), bis(D, right(5)), mink(D, right(4)))
        //         * Q(0, mink(D, 1338))
        //         * Q(1, mink(D, 1337))
        //         * Q(3, mink(D, left(7)))
        //         * Q(3, mink(D, right(5))),
        //     "spenso"
        // );
        //
        let expr = parse_lit!(
            ((-1 * g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                + g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9))))
                * (-1 * g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                    + g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9))))
                * 2
                * G
                ^ 6 * P(2, mink(4, dummy(2, 2)))
                    * P(3, mink(4, dummy(3, 3)))
                    * g(bis(4, l(2)), bis(4, l(5)))
                    * g(bis(4, l(3)), bis(4, l(6)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(mink(4, l(0)), mink(4, l(6)))
                    * g(mink(4, l(0)), mink(4, r(0)))
                    * g(mink(4, l(1)), mink(4, l(7)))
                    * g(mink(4, l(1)), mink(4, r(1)))
                    * g(mink(4, l(4)), mink(4, l(8)))
                    * g(mink(4, l(4)), mink(4, r(4)))
                    * g(mink(4, l(5)), mink(4, l(9)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * gamma(bis(4, l(2)), bis(4, r(2)), mink(4, dummy(2, 2)))
                    * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                    * gamma(bis(4, r(3)), bis(4, l(3)), mink(4, dummy(3, 3)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    + (-1 * g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                        + g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9))))
                        * (-1 * g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                            + g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                        * G
                ^ 6 * P(2, mink(4, dummy(2, 2)))
                    * P(3, mink(4, dummy(3, 3)))
                    * g(bis(4, l(2)), bis(4, l(5)))
                    * g(bis(4, l(3)), bis(4, l(6)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(mink(4, l(0)), mink(4, l(6)))
                    * g(mink(4, l(0)), mink(4, r(0)))
                    * g(mink(4, l(1)), mink(4, l(7)))
                    * g(mink(4, l(1)), mink(4, r(1)))
                    * g(mink(4, l(4)), mink(4, l(8)))
                    * g(mink(4, l(4)), mink(4, r(4)))
                    * g(mink(4, l(5)), mink(4, l(9)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * gamma(bis(4, l(2)), bis(4, r(2)), mink(4, dummy(2, 2)))
                    * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                    * gamma(bis(4, r(3)), bis(4, l(3)), mink(4, dummy(3, 3)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    + (-1 * g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                        + g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9))))
                        * (-1 * g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9)))
                            + g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                        * -1
                        * G
                ^ 6 * P(2, mink(4, dummy(2, 2)))
                    * P(3, mink(4, dummy(3, 3)))
                    * g(bis(4, l(2)), bis(4, l(5)))
                    * g(bis(4, l(3)), bis(4, l(6)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(mink(4, l(0)), mink(4, l(6)))
                    * g(mink(4, l(0)), mink(4, r(0)))
                    * g(mink(4, l(1)), mink(4, l(7)))
                    * g(mink(4, l(1)), mink(4, r(1)))
                    * g(mink(4, l(4)), mink(4, l(8)))
                    * g(mink(4, l(4)), mink(4, r(4)))
                    * g(mink(4, l(5)), mink(4, l(9)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * gamma(bis(4, l(2)), bis(4, r(2)), mink(4, dummy(2, 2)))
                    * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                    * gamma(bis(4, r(3)), bis(4, l(3)), mink(4, dummy(3, 3)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    + (-1 * g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                        + g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                        * (-1 * g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                            + g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9))))
                        * G
                ^ 6 * P(2, mink(4, dummy(2, 2)))
                    * P(3, mink(4, dummy(3, 3)))
                    * g(bis(4, l(2)), bis(4, l(5)))
                    * g(bis(4, l(3)), bis(4, l(6)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(mink(4, l(0)), mink(4, l(6)))
                    * g(mink(4, l(0)), mink(4, r(0)))
                    * g(mink(4, l(1)), mink(4, l(7)))
                    * g(mink(4, l(1)), mink(4, r(1)))
                    * g(mink(4, l(4)), mink(4, l(8)))
                    * g(mink(4, l(4)), mink(4, r(4)))
                    * g(mink(4, l(5)), mink(4, l(9)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * gamma(bis(4, l(2)), bis(4, r(2)), mink(4, dummy(2, 2)))
                    * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                    * gamma(bis(4, r(3)), bis(4, l(3)), mink(4, dummy(3, 3)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    + (-1 * g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                        + g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                        * (-1 * g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                            + g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                        * 2
                        * G
                ^ 6 * P(2, mink(4, dummy(2, 2)))
                    * P(3, mink(4, dummy(3, 3)))
                    * g(bis(4, l(2)), bis(4, l(5)))
                    * g(bis(4, l(3)), bis(4, l(6)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(mink(4, l(0)), mink(4, l(6)))
                    * g(mink(4, l(0)), mink(4, r(0)))
                    * g(mink(4, l(1)), mink(4, l(7)))
                    * g(mink(4, l(1)), mink(4, r(1)))
                    * g(mink(4, l(4)), mink(4, l(8)))
                    * g(mink(4, l(4)), mink(4, r(4)))
                    * g(mink(4, l(5)), mink(4, l(9)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * gamma(bis(4, l(2)), bis(4, r(2)), mink(4, dummy(2, 2)))
                    * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                    * gamma(bis(4, r(3)), bis(4, l(3)), mink(4, dummy(3, 3)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    + (-1 * g(mink(4, l(6)), mink(4, l(7))) * g(mink(4, l(8)), mink(4, l(9)))
                        + g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                        * (-1 * g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9)))
                            + g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                        * G
                ^ 6 * P(2, mink(4, dummy(2, 2)))
                    * P(3, mink(4, dummy(3, 3)))
                    * g(bis(4, l(2)), bis(4, l(5)))
                    * g(bis(4, l(3)), bis(4, l(6)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(mink(4, l(0)), mink(4, l(6)))
                    * g(mink(4, l(0)), mink(4, r(0)))
                    * g(mink(4, l(1)), mink(4, l(7)))
                    * g(mink(4, l(1)), mink(4, r(1)))
                    * g(mink(4, l(4)), mink(4, l(8)))
                    * g(mink(4, l(4)), mink(4, r(4)))
                    * g(mink(4, l(5)), mink(4, l(9)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * gamma(bis(4, l(2)), bis(4, r(2)), mink(4, dummy(2, 2)))
                    * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                    * gamma(bis(4, r(3)), bis(4, l(3)), mink(4, dummy(3, 3)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    + (-1 * g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9)))
                        + g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                        * (-1 * g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                            + g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9))))
                        * -1
                        * G
                ^ 6 * P(2, mink(4, dummy(2, 2)))
                    * P(3, mink(4, dummy(3, 3)))
                    * g(bis(4, l(2)), bis(4, l(5)))
                    * g(bis(4, l(3)), bis(4, l(6)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(mink(4, l(0)), mink(4, l(6)))
                    * g(mink(4, l(0)), mink(4, r(0)))
                    * g(mink(4, l(1)), mink(4, l(7)))
                    * g(mink(4, l(1)), mink(4, r(1)))
                    * g(mink(4, l(4)), mink(4, l(8)))
                    * g(mink(4, l(4)), mink(4, r(4)))
                    * g(mink(4, l(5)), mink(4, l(9)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * gamma(bis(4, l(2)), bis(4, r(2)), mink(4, dummy(2, 2)))
                    * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                    * gamma(bis(4, r(3)), bis(4, l(3)), mink(4, dummy(3, 3)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    + (-1 * g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9)))
                        + g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                        * (-1 * g(mink(4, r(6)), mink(4, r(7))) * g(mink(4, r(8)), mink(4, r(9)))
                            + g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                        * G
                ^ 6 * P(2, mink(4, dummy(2, 2)))
                    * P(3, mink(4, dummy(3, 3)))
                    * g(bis(4, l(2)), bis(4, l(5)))
                    * g(bis(4, l(3)), bis(4, l(6)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(mink(4, l(0)), mink(4, l(6)))
                    * g(mink(4, l(0)), mink(4, r(0)))
                    * g(mink(4, l(1)), mink(4, l(7)))
                    * g(mink(4, l(1)), mink(4, r(1)))
                    * g(mink(4, l(4)), mink(4, l(8)))
                    * g(mink(4, l(4)), mink(4, r(4)))
                    * g(mink(4, l(5)), mink(4, l(9)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * gamma(bis(4, l(2)), bis(4, r(2)), mink(4, dummy(2, 2)))
                    * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                    * gamma(bis(4, r(3)), bis(4, l(3)), mink(4, dummy(3, 3)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5)))
                    + (-1 * g(mink(4, l(6)), mink(4, l(8))) * g(mink(4, l(7)), mink(4, l(9)))
                        + g(mink(4, l(6)), mink(4, l(9))) * g(mink(4, l(7)), mink(4, l(8))))
                        * (-1 * g(mink(4, r(6)), mink(4, r(8))) * g(mink(4, r(7)), mink(4, r(9)))
                            + g(mink(4, r(6)), mink(4, r(9))) * g(mink(4, r(7)), mink(4, r(8))))
                        * 2
                        * G
                ^ 6 * P(2, mink(4, dummy(2, 2)))
                    * P(3, mink(4, dummy(3, 3)))
                    * g(bis(4, l(2)), bis(4, l(5)))
                    * g(bis(4, l(3)), bis(4, l(6)))
                    * g(bis(4, r(2)), bis(4, r(5)))
                    * g(bis(4, r(3)), bis(4, r(6)))
                    * g(mink(4, l(0)), mink(4, l(6)))
                    * g(mink(4, l(0)), mink(4, r(0)))
                    * g(mink(4, l(1)), mink(4, l(7)))
                    * g(mink(4, l(1)), mink(4, r(1)))
                    * g(mink(4, l(4)), mink(4, l(8)))
                    * g(mink(4, l(4)), mink(4, r(4)))
                    * g(mink(4, l(5)), mink(4, l(9)))
                    * g(mink(4, r(0)), mink(4, r(6)))
                    * g(mink(4, r(1)), mink(4, r(7)))
                    * g(mink(4, r(4)), mink(4, r(8)))
                    * g(mink(4, r(5)), mink(4, r(9)))
                    * gamma(bis(4, l(2)), bis(4, r(2)), mink(4, dummy(2, 2)))
                    * gamma(bis(4, l(6)), bis(4, l(5)), mink(4, l(5)))
                    * gamma(bis(4, r(3)), bis(4, l(3)), mink(4, dummy(3, 3)))
                    * gamma(bis(4, r(5)), bis(4, r(6)), mink(4, r(5))))
                * -18,
            default_namespace = "spenso"
        );

        let res = parse_lit!(7776 * G ^ 6 * dot(P(2), P(3)), default_namespace = "spenso");
        assert_eq!(
            res,
            expr.simplify_gamma().to_dots(),
            "fount{}",
            expr.simplify_gamma().to_dots()
        );

        assert_eq!(
            res,
            expr.simplify_gamma().to_dots(),
            "fount{}",
            expr.simplify_gamma().to_dots()
        );

        let expr = parse_lit!(
            G ^ 2
                * Q(EMRID(0, 4), mink(dim, l(20)))
                * g(bis(4, l(2)), bis(4, l(4)))
                * g(bis(4, l(3)), bis(4, l(7)))
                * g(mink(dim, l(0)), mink(dim, l(5)))
                * g(mink(dim, l(1)), mink(dim, l(4)))
                * gamma(bis(4, l(5)), bis(4, l(4)), mink(dim, l(4)))
                * gamma(bis(4, l(6)), bis(4, l(5)), mink(dim, l(20)))
                * gamma(bis(4, l(7)), bis(4, l(6)), mink(dim, l(5))),
            default_namespace = "spenso"
        );

        println!("{:>}", expr.simplify_gamma().expand());
        // assert_eq!(
        //     expr.simplify_gamma().to_dots(),
        //     expr.simplify_gamma().simplify_gamma().to_dots(),
        //     "\n{:>}\n not equal to \n{:>}\n diff:\n{:>}",
        //     expr.simplify_gamma().to_dots(),
        //     expr.simplify_gamma().simplify_gamma().to_dots(),
        //     (expr.simplify_gamma().simplify_gamma().to_dots() - expr.simplify_gamma().to_dots())
        //         .expand()
        // )
        // println!("{}", expr.simplify_gamma().to_dots());

        // println!("{}", SpinAntiFundamental {}.to_symbolic([RS.a_]))
    }
}
