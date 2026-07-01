use super::*;
use insta::assert_snapshot;
use spenso::{g, p, q, trace};

#[test]
fn dirac_simplify_id1_repeated_d_dim_gamma_is_dimension() {
    let r = test_initialize();
    let expr = chain!(
        slot!(r.bis_d, i),
        slot!(r.bis_d, j),
        gamma!(slot!(r.mink_d, mu)),
        gamma!(slot!(r.mink_d, mu)),
    );

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"d*g(bis(d,i),bis(d,j))");
}

#[test]
fn dirac_simplify_id2_odd_interior_chain() {
    let r = test_initialize();
    let expr = chain!(
        slot!(r.bis4, a),
        slot!(r.bis4, b),
        gamma!(slot!(r.mink4, mu)),
        gamma!(slot!(r.mink4, nu)),
        gamma!(slot!(r.mink4, rho)),
        gamma!(slot!(r.mink4, sigma)),
        gamma!(slot!(r.mink4, mu)),
    );

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"-2*chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,sigma)),gamma(in,out,mink(4,rho)),gamma(in,out,mink(4,nu)))");
}

#[test]
fn dirac_simplify_id3_four_interior_chain() {
    let r = test_initialize();
    let expr = chain!(
        slot!(r.bis4, i),
        slot!(r.bis4, j),
        gamma!(slot!(r.mink4, mu)),
        gamma!(slot!(r.mink4, alpha)),
        gamma!(slot!(r.mink4, beta)),
        gamma!(slot!(r.mink4, rho)),
        gamma!(slot!(r.mink4, sigma)),
        gamma!(slot!(r.mink4, mu)),
    ) / 2;
    let simplified = expr.simplify_gamma();
    println!(
        "{}",
        simplified.spenso_print(&SpensoPrintSettings::compact())
    );
    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"(2*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,rho)),gamma(in,out,mink(4,beta)),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,sigma)))+2*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,sigma)),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,beta)),gamma(in,out,mink(4,rho))))*1/2");
    let simplified = expr.simplify_gamma_with(GammaSimplifySettings::canonical());
    println!(
        "{}",
        simplified.spenso_print(&SpensoPrintSettings::compact())
    );
    assert_snapshot!(&expr.simplify_gamma_with(GammaSimplifySettings::canonical()).to_bare_ordered_string(), @"((((-1*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,beta)),gamma(in,out,mink(4,rho)),gamma(in,out,mink(4,sigma)))+2*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,beta)))*g(mink(4,rho),mink(4,sigma)))*-1+2*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,rho)))*g(mink(4,beta),mink(4,sigma)))*-1+2*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,beta)),gamma(in,out,mink(4,rho)))*g(mink(4,alpha),mink(4,sigma)))*2+(((-1*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,beta)),gamma(in,out,mink(4,rho)),gamma(in,out,mink(4,sigma)))+2*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,rho)),gamma(in,out,mink(4,sigma)))*g(mink(4,alpha),mink(4,beta)))*-1+2*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,beta)),gamma(in,out,mink(4,sigma)))*g(mink(4,alpha),mink(4,rho)))*-1+2*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,sigma)))*g(mink(4,beta),mink(4,rho)))*2)*1/2");
}

#[test]
fn dirac_simplify_id4_slash_sandwich() {
    let r = test_initialize();
    let p = p!(r.mink4);
    let q = q!(r.mink4);
    let expr = Atom::var(s!(m))
        * gamma!(p.clone(), slot!(r.bis4, i), slot!(r.bis4, a))
        * gamma!(p.clone(), slot!(r.bis4, a), slot!(r.bis4, j))
        - gamma!(p.clone(), slot!(r.bis4, i), slot!(r.bis4, a))
            * gamma!(q.clone(), slot!(r.bis4, a), slot!(r.bis4, b))
            * gamma!(p.clone(), slot!(r.bis4, b), slot!(r.bis4, j));

    assert_snapshot!(expr.simplify_gamma().expand().to_bare_ordered_string(), @"-2*chain(bis(4,i),bis(4,j),gamma(in,out,p(mink(4))))*g(p(mink(4)),q(mink(4)))+chain(bis(4,i),bis(4,j),gamma(in,out,q(mink(4))))*g(p(mink(4)),p(mink(4)))+g(bis(4,i),bis(4,j))*g(p(mink(4)),p(mink(4)))*m");
}

#[test]
fn dirac_simplify_id5_gamma5_anticommutes_left() {
    let r = test_initialize();
    let expr = chain!(
        slot!(r.bis4, i),
        slot!(r.bis4, j),
        gamma5!(),
        gamma!(slot!(r.mink4, mu)),
    );

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"-1*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,mu)),gamma5(in,out))");
}

#[test]
fn dirac_simplify_id23_trace_evaluation_can_stay_disabled() {
    let r = test_initialize();
    let expr = trace!(
        r.bis4.to_symbolic([]),
        gamma!(slot!(r.mink4, a)),
        gamma!(slot!(r.mink4, b)),
        gamma!(slot!(r.mink4, a))
    ) + gamma!(a, 1, 2) * gamma!(a, 2, 3);

    assert_snapshot!(expr.simplify_gamma_with(GammaSimplifySettings::repeated_pairs().without_trace_evaluation()).to_bare_ordered_string(), @"4*g(bis(4,1),bis(4,3))+trace(bis(4),cyclic(gamma(in,out,mink(4,a)),gamma(in,out,mink(4,a)),gamma(in,out,mink(4,b))))");
}

#[test]
fn dirac_simplify_id24_odd_trace_vanishes() {
    let r = test_initialize();
    let expr = trace!(
        r.bis4.to_symbolic([]),
        gamma!(slot!(r.mink4, a)),
        gamma!(slot!(r.mink4, b)),
        gamma!(slot!(r.mink4, a))
    ) + chain!(
        slot!(r.bis4, i),
        slot!(r.bis4, j),
        gamma!(slot!(r.mink4, a)),
        gamma!(slot!(r.mink4, b)),
    );

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,a)),gamma(in,out,mink(4,b)))");
}

#[test]
fn dirac_simplify_id25_four_trace_recurses_beside_open_chain() {
    let r = test_initialize();
    let expr = trace!(
        r.bis4.to_symbolic([]),
        gamma!(slot!(r.mink4, a)),
        gamma!(slot!(r.mink4, b)),
        gamma!(slot!(r.mink4, c)),
        gamma!(slot!(r.mink4, d))
    ) + chain!(
        slot!(r.bis4, i),
        slot!(r.bis4, j),
        gamma!(slot!(r.mink4, a)),
        gamma!(slot!(r.mink4, b)),
        gamma!(slot!(r.mink4, c)),
        gamma!(slot!(r.mink4, d)),
    );

    assert_snapshot!(expr.simplify_gamma().expand().to_bare_ordered_string(), @"-4*g(mink(4,a),mink(4,c))*g(mink(4,b),mink(4,d))+4*g(mink(4,a),mink(4,b))*g(mink(4,c),mink(4,d))+4*g(mink(4,a),mink(4,d))*g(mink(4,b),mink(4,c))+chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,a)),gamma(in,out,mink(4,b)),gamma(in,out,mink(4,c)),gamma(in,out,mink(4,d)))");
}

#[test]
fn dirac_simplify_id30_dirac_order_anticommutator() {
    let r = test_initialize();
    let p = p!(&r.mink4);
    let expr = chain!(
        slot!(r.bis4, i),
        slot!(r.bis4, j),
        gamma!(slot!(r.mink4, nu)),
        gamma!(p.clone()),
    ) + chain!(
        slot!(r.bis4, i),
        slot!(r.bis4, j),
        gamma!(p.clone()),
        gamma!(slot!(r.mink4, nu)),
    ) - Atom::num(2)
        * g!(slot!(r.mink4, nu), p)
        * chain!(slot!(r.bis4, i), slot!(r.bis4, j); Vec::<Atom>::new());

    assert_snapshot!(expr.simplify_gamma_with(GammaSimplifySettings::canonical()).to_bare_ordered_string(), @"0");
}

#[test]
fn dirac_simplify_id36_empty_trace_is_spin_dimension() {
    let r = test_initialize();
    let expr_4 = trace!(r.bis4.to_symbolic([]); Vec::<Atom>::new());
    let expr_d = trace!(r.bis_d.to_symbolic([]); Vec::<Atom>::new());
    let expr_color = trace!(r.cof_nc.to_symbolic([]); Vec::<Atom>::new());
    let pair_d = trace!(
        r.bis_d.to_symbolic([]),
        gamma!(slot!(r.mink_d, a)),
        gamma!(slot!(r.mink_d, b))
    );

    assert_snapshot!(expr_4.simplify_gamma().to_bare_ordered_string(), @"4");
    assert_snapshot!(expr_d.simplify_gamma().to_bare_ordered_string(), @"d");
    assert_snapshot!(expr_color.simplify_gamma().to_bare_ordered_string(), @"Nc");
    assert_snapshot!(pair_d.simplify_gamma().to_bare_ordered_string(), @"d*g(mink(d,a),mink(d,b))");
}

#[test]
fn dirac_simplify_id40_repeated_gamma_and_two_trace() {
    let r = test_initialize();
    let p = p!(&r.mink4);
    let open = Atom::var(s!(c1))
        * gamma!(mu, slot!(r.bis4, i), slot!(r.bis4, a))
        * gamma!(p.clone(), slot!(r.bis4, a), slot!(r.bis4, b))
        * gamma!(mu, slot!(r.bis4, b), slot!(r.bis4, j))
        + Atom::var(s!(c1))
            * Atom::var(s!(m))
            * gamma!(mu, slot!(r.bis4, i), slot!(r.bis4, a))
            * gamma!(mu, slot!(r.bis4, a), slot!(r.bis4, j));
    let tr = Atom::var(s!(c2))
        * trace!(
            r.bis4.to_symbolic([]),
            gamma!(slot!(r.mink4, mu)),
            gamma!(slot!(r.mink4, nu))
        );

    assert_snapshot!((open + tr).simplify_gamma().expand().to_bare_ordered_string(), @"-2*c1*chain(bis(4,i),bis(4,j),gamma(in,out,p(mink(4))))+4*c1*g(bis(4,i),bis(4,j))*m+4*c2*g(mink(4,mu),mink(4,nu))");
}

#[test]
fn dirac_simplify_id45_slash_square_and_sandwich() {
    let r = test_initialize();
    let p = p!(&r.mink4);
    let slash_square = g!(p.clone(), p.clone())
        * chain!(
            slot!(r.bis4, i),
            slot!(r.bis4, j),
            gamma!(slot!(r.mink4, mu)),
        );
    let sandwich = chain!(
        slot!(r.bis4, i),
        slot!(r.bis4, j),
        gamma!(p.clone()),
        gamma!(slot!(r.mink4, mu)),
        gamma!(p),
    );

    assert_snapshot!(slash_square.simplify_gamma().to_bare_ordered_string(), @"chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,mu)))*g(p(mink(4)),p(mink(4)))");
    assert_snapshot!(sandwich.simplify_gamma().expand().to_bare_ordered_string(), @"-1*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,mu)))*g(p(mink(4)),p(mink(4)))+2*chain(bis(4,i),bis(4,j),gamma(in,out,p(mink(4))))*p(mink(4,mu))");
}

#[test]
fn dirac_simplify_id46_d_dim_repeated_gamma_slash_sum() {
    let r = test_initialize();
    let p = p!(&r.mink_d);
    let q = q!(&r.mink_d);
    let expr = gamma!(slot!(r.mink_d, mu), slot!(r.bis_d, i), slot!(r.bis_d, a))
        * gamma!(p.clone(), slot!(r.bis_d, a), slot!(r.bis_d, b))
        * gamma!(slot!(r.mink_d, mu), slot!(r.bis_d, b), slot!(r.bis_d, j))
        + gamma!(slot!(r.mink_d, mu), slot!(r.bis_d, i), slot!(r.bis_d, a))
            * gamma!(q.clone(), slot!(r.bis_d, a), slot!(r.bis_d, b))
            * gamma!(slot!(r.mink_d, mu), slot!(r.bis_d, b), slot!(r.bis_d, j))
        + Atom::var(s!(m))
            * gamma!(slot!(r.mink_d, mu), slot!(r.bis_d, i), slot!(r.bis_d, a))
            * gamma!(slot!(r.mink_d, mu), slot!(r.bis_d, a), slot!(r.bis_d, j));

    assert_snapshot!(expr.simplify_gamma().expand().to_bare_ordered_string(), @"-1*chain(bis(d,i),bis(d,j),gamma(in,out,p(mink(d))))*d+-1*chain(bis(d,i),bis(d,j),gamma(in,out,q(mink(d))))*d+2*chain(bis(d,i),bis(d,j),gamma(in,out,p(mink(d))))+2*chain(bis(d,i),bis(d,j),gamma(in,out,q(mink(d))))+d*g(bis(d,i),bis(d,j))*m");
}
