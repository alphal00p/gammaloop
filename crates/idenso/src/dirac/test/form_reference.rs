use super::*;
use insta::assert_snapshot;

#[test]
fn two_gamma_trace() {
    initialize();
    let expr = gamma!(mu, a, b) * gamma!(nu, b, a);

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"4*g(mink(4,mu),mink(4,nu))");
}

#[test]
fn odd_gamma_trace_vanishes() {
    initialize();
    let expr = gamma!(mu, a, b) * gamma!(nu, b, c) * gamma!(rho, c, a);
    assert!(expr.simplify_gamma().is_zero());
}

#[test]
fn four_gamma_trace_recurses() {
    initialize();
    let expr = gamma!(mu, a, b) * gamma!(nu, b, c) * gamma!(rho, c, d) * gamma!(sigma, d, a);

    assert_snapshot!(expr.simplify_gamma().expand().to_bare_ordered_string(), @"-4*g(mink(4,mu),mink(4,rho))*g(mink(4,nu),mink(4,sigma))+4*g(mink(4,mu),mink(4,nu))*g(mink(4,rho),mink(4,sigma))+4*g(mink(4,mu),mink(4,sigma))*g(mink(4,nu),mink(4,rho))");
}

#[test]
fn repeated_lorentz_gamma_chain_contracts_to_dimension() {
    initialize();
    let expr = gamma!(mu, a, b) * gamma!(mu, b, c);

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"4*g(bis(4,a),bis(4,c))");
}

#[test]
fn adjacent_chain_lorentz_contraction() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.bis_d, a),
        slot!(r.bis_d, b),
        gamma!(slot!(r.mink_d, mu)),
        gamma!(slot!(r.mink_d, mu)),
    );

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"d*g(bis(d,a),bis(d,b))");
}

#[test]
fn trace4gen_chisholm_odd_interior_chain() {
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

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"-2*chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,nu3)),gamma(in,out,mink(4,nu2)),gamma(in,out,mink(4,nu1)))");
}

#[test]
fn trace4gen_chisholm_two_interior_chain() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.bis4, a),
        slot!(r.bis4, b),
        gamma!(slot!(r.mink4, mu)),
        gamma!(slot!(r.mink4, nu)),
        gamma!(slot!(r.mink4, rho)),
        gamma!(slot!(r.mink4, mu)),
    );

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"4*g(bis(4,a),bis(4,b))*g(mink(4,nu),mink(4,rho))");
}

#[test]
#[ignore = "pending public chain-based gamma-five epsilon reduction"]
fn gamma_five_epsilon_trick() {
    let r = TestReps::new();
    let expr = chain!(
        slot!(r.bis4, a),
        slot!(r.bis4, b),
        gamma!(slot!(r.mink4, mu)),
        gamma!(slot!(r.mink4, nu)),
        gamma!(slot!(r.mink4, rho)),
    );

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"chain(bis(4,a),bis(4,b),gamma5(in,out),gamma(in,out,mink(4,sigma)))*epsilon(mink(4,mu),mink(4,nu),mink(4,rho),mink(4,sigma))+chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,rho)))*g(mink(4,mu),mink(4,nu))+-1*chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,nu)))*g(mink(4,mu),mink(4,rho))+chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,mu)))*g(mink(4,nu),mink(4,rho))");
}

#[test]
#[ignore = "pending gamma-five trace term-count stress support"]
fn gamma_five_symmetric_12_term_count() {
    initialize();
    let term_count = 0usize;

    assert_snapshot!(term_count.to_string(), @"51975");
}

#[test]
#[ignore = "pending gamma-five trace4 term-count stress support"]
fn gamma_five_regular_12_term_count() {
    initialize();
    let term_count = 0usize;

    assert_snapshot!(term_count.to_string(), @"1029");
}
