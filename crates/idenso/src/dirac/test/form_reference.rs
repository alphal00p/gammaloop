use super::*;
use insta::assert_snapshot;
use symbolica::atom::AtomView;

#[test]
fn two_gamma_trace() {
    test_initialize();
    let expr = gamma!(a, b, mu) * gamma!(b, a, nu);

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"4*g(mink(4,mu),mink(4,nu))");
}

#[test]
fn odd_gamma_trace_vanishes() {
    test_initialize();
    let expr = gamma!(a, b, mu) * gamma!(b, c, nu) * gamma!(c, a, rho);
    assert!(expr.simplify_gamma().is_zero());
}

#[test]
fn four_gamma_trace_recurses() {
    test_initialize();
    let expr = gamma!(a, b, mu) * gamma!(b, c, nu) * gamma!(c, d, rho) * gamma!(d, a, sigma);

    assert_snapshot!(expr.simplify_gamma().expand().to_bare_ordered_string(), @"-4*g(mink(4,mu),mink(4,rho))*g(mink(4,nu),mink(4,sigma))+4*g(mink(4,mu),mink(4,nu))*g(mink(4,rho),mink(4,sigma))+4*g(mink(4,mu),mink(4,sigma))*g(mink(4,nu),mink(4,rho))");
}

#[test]
fn repeated_lorentz_gamma_chain_contracts_to_dimension() {
    test_initialize();
    let expr = gamma!(a, b, mu) * gamma!(b, c, mu);

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"4*g(bis(4,a),bis(4,c))");
}

#[test]
fn adjacent_chain_lorentz_contraction() {
    let r = test_initialize();
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
    let r = test_initialize();
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
    let r = test_initialize();
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
fn gamma_five_epsilon_trick() {
    let r = test_initialize();
    let expr = chain!(
        slot!(r.bis4, a),
        slot!(r.bis4, b),
        gamma!(slot!(r.mink4, mu)),
        gamma!(slot!(r.mink4, nu)),
        gamma!(slot!(r.mink4, rho)),
    );

    assert_snapshot!(expr
        .simplify_gamma_with(GammaSimplifySettings::repeated_pairs().with_gamma5_epsilon_expansion())
        .to_bare_ordered_string(), @"-1*chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,nu)))*g(mink(4,mu),mink(4,rho))+-1*chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,sigma)),gamma5(in,out))*epsilon(mink(4,mu),mink(4,nu),mink(4,rho),mink(4,sigma))+chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,mu)))*g(mink(4,nu),mink(4,rho))+chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,rho)))*g(mink(4,mu),mink(4,nu))");
}

#[test]
fn gamma_five_four_gamma_trace_is_epsilon() {
    let r = test_initialize();
    let expr = trace!(
        r.bis4.to_symbolic([]),
        gamma5!(),
        gamma!(slot!(r.mink4, mu)),
        gamma!(slot!(r.mink4, nu)),
        gamma!(slot!(r.mink4, rho)),
        gamma!(slot!(r.mink4, sigma)),
    );

    assert_snapshot!(expr.simplify_gamma().to_bare_ordered_string(), @"4*epsilon(mink(4,mu),mink(4,nu),mink(4,rho),mink(4,sigma))");
}

#[test]
fn gamma_five_two_gamma_trace_vanishes() {
    let r = test_initialize();
    let expr = trace!(
        r.bis4.to_symbolic([]),
        gamma5!(),
        gamma!(slot!(r.mink4, mu)),
        gamma!(slot!(r.mink4, nu)),
    );

    assert!(expr.simplify_gamma().is_zero());
}

#[test]
fn short_trace_dispatch_covers_every_arity_and_gamma5() {
    let r = test_initialize();
    let ordinary_counts = [1, 3, 15, 105, 693, 4383, 26931];
    let axial_counts = [0, 1, 6, 33, 180, 1029, 6042];
    let indices = (0..14)
        .map(|i| {
            r.mink4
                .pattern(symbolica::symbol!(&format!("short_trace_mu{i}")))
        })
        .collect::<Vec<_>>();
    for n in 1..=14 {
        for axial in [false, true] {
            let mut factors = indices[..n].iter().map(|mu| gamma!(mu)).collect::<Vec<_>>();
            if axial {
                factors.insert(0, gamma5!());
            }
            let expr = trace!(r.bis4.to_symbolic([]); factors);
            let result = expr.simplify_gamma();
            let expected = if n % 2 == 1 {
                0
            } else if axial {
                axial_counts[n / 2 - 1]
            } else {
                ordinary_counts[n / 2 - 1]
            };
            let count = if result.is_zero() {
                0
            } else {
                result.expand().nterms()
            };
            assert_eq!(count, expected, "length {n}, axial {axial}");
            assert_eq!(
                result.simplify_gamma(),
                result,
                "trace result must be a fixed point"
            );
        }
    }
}

#[test]
fn chisholm_trace_reduction_is_cyclic_and_preserves_coefficients() {
    let r = test_initialize();
    let mu = gamma!(slot!(r.mink4, mu));
    let middle = [
        gamma!(slot!(r.mink4, nu)),
        gamma!(slot!(r.mink4, rho)),
        gamma!(slot!(r.mink4, alpha)),
    ];
    let tail = [
        gamma!(slot!(r.mink4, beta)),
        gamma!(slot!(r.mink4, sigma)),
        gamma!(slot!(r.mink4, tau)),
    ];
    let mut factors = [vec![mu.clone()], middle.to_vec(), vec![mu], tail.to_vec()].concat();
    let coefficient = parse_lit!((x + y) ^ 8);
    let expected = Atom::num(-2)
        * trace!(r.bis4.to_symbolic([]); middle.iter().rev().chain(tail.iter())).simplify_gamma();
    let expected = expected.expand();
    for _ in 0..factors.len() {
        let expr = &coefficient * trace!(r.bis4.to_symbolic([]); &factors);
        let result = expr.simplify_gamma();
        // The metric polynomial may be factored differently after rotation.
        // Compare trace bodies while requiring the spectator to stay intact.
        assert!(matches!(result.as_view(), AtomView::Mul(product)
            if product.iter().any(|factor| factor == coefficient.as_view())));
        assert_eq!((&result / &coefficient).expand(), expected);
        assert_eq!(result.simplify_gamma(), result);
        factors.rotate_left(1);
    }
}

#[test]
fn long_chisholm_interiors_use_shared_open_chain_identity() {
    let r = test_initialize();
    let mu = gamma!(slot!(r.mink4, mu));
    let interior = (0..6)
        .map(|i| {
            gamma!(
                r.mink4
                    .pattern(symbolica::symbol!(&format!("long_chisholm_mu{i}")))
            )
        })
        .collect::<Vec<_>>();
    for n in [5, 6] {
        let word = [vec![mu.clone()], interior[..n].to_vec(), vec![mu.clone()]].concat();
        let result = chain!(slot!(r.bis4, a), slot!(r.bis4, b); word).simplify_gamma();
        let expected = if n % 2 == 1 {
            -2 * chain!(slot!(r.bis4, a), slot!(r.bis4, b); interior[..n].iter().rev())
        } else {
            2 * chain!(slot!(r.bis4, a), slot!(r.bis4, b); interior[..n-1].iter().rev().chain([&interior[n-1]]))
                + 2 * chain!(slot!(r.bis4, a), slot!(r.bis4, b); [&interior[n-1]].into_iter().chain(interior[..n-1].iter()))
        };
        assert_eq!(result, expected);
    }
}

#[test]
fn symbolic_dimension_keeps_generic_trace_recursion() {
    let r = test_initialize();
    let factors = (0..10)
        .map(|i| {
            gamma!(
                r.mink_d
                    .pattern(symbolica::symbol!(&format!("generic_trace_mu{i}")))
            )
        })
        .collect::<Vec<_>>();
    let result = trace!(r.bis4.to_symbolic([]); factors).simplify_gamma();
    assert_eq!(result.expand().nterms(), 945);
}

#[test]
fn short_trace_terminal_shortcut_preserves_surrounding_contractions() {
    let r = test_initialize();
    let factors = (0..10)
        .map(|i| {
            gamma!(
                r.mink4
                    .pattern(symbolica::symbol!(&format!("terminal_trace_mu{i}")))
            )
        })
        .collect::<Vec<_>>();
    let expr = trace!(r.bis4.to_symbolic([]); factors);
    let coefficient = parse_lit!((x + y) ^ 8);
    let result = (&coefficient * &expr).simplify_gamma();
    assert!(matches!(result.as_view(), AtomView::Mul(product)
        if product.iter().any(|factor| factor == coefficient.as_view())));
    // Expand only the trace body; the scalar spectator remains factored.
    assert_eq!(
        (&result / &coefficient).expand(),
        expr.simplify_gamma().expand()
    );
    assert_eq!(result.simplify_gamma(), result);

    let mu = r.mink4.pattern(s!(mu));
    let nu = r.mink4.pattern(s!(nu));
    let alpha = r.mink4.pattern(s!(alpha));
    let beta = r.mink4.pattern(s!(beta));
    let contracted = g!(&mu, &nu)
        * trace!(
            r.bis4.to_symbolic([]),
            gamma!(&mu),
            gamma!(&alpha),
            gamma!(&nu),
            gamma!(&beta)
        );
    assert_eq!(
        contracted.simplify_gamma().expand().simplify_metrics(),
        -8 * g!(&alpha, &beta)
    );
}
