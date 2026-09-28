use super::*;
use crate::shorthands::{UndoShorthands, chain::Chain};
use insta::assert_snapshot;
use spenso::{g, network::tags::SPENSO_TAG as T, p, q, trace};

#[test]
fn dirac_simplify_id1_repeated_d_dim_gamma_is_dimension() {
    let r = test_initialize();
    let expr = chain!(
        slot!(r.bis_d, i),
        slot!(r.bis_d, j),
        gamma!(slot!(r.mink_d, mu)),
        gamma!(slot!(r.mink_d, mu)),
    );

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_gamma(crate::dirac::GammaSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"d*g(bis(d,i),bis(d,j))");
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

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_gamma(crate::dirac::GammaSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"-2*chain(bis(4,a),bis(4,b),gamma(in,out,mink(4,sigma)),gamma(in,out,mink(4,rho)),gamma(in,out,mink(4,nu)))");
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
    let simplified = crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
        .unwrap()
        .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
    println!(
        "{}",
        simplified.spenso_print(&SpensoPrintSettings::compact())
    );
    // Keep the original arithmetic oracle without distributing the numerator.
    crate::test_support::assert_factored_snapshot_eq(
        &simplified.to_bare_ordered_string(),
        "chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,rho)),gamma(in,out,mink(4,beta)),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,sigma)))+chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,sigma)),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,beta)),gamma(in,out,mink(4,rho)))",
    );
    let simplified = crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
        .unwrap()
        .simplify_gamma(GammaSimplifySettings::canonical())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression();
    println!(
        "{}",
        simplified.spenso_print(&SpensoPrintSettings::compact())
    );
    crate::test_support::assert_factored_snapshot_eq(
        &simplified.to_bare_ordered_string(),
        "(-4*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,beta)),gamma(in,out,mink(4,rho)),gamma(in,out,mink(4,sigma)))+-4*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,rho)))*g(mink(4,beta),mink(4,sigma))+-4*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,beta)),gamma(in,out,mink(4,sigma)))*g(mink(4,alpha),mink(4,rho))+4*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,beta)))*g(mink(4,rho),mink(4,sigma))+4*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,alpha)),gamma(in,out,mink(4,sigma)))*g(mink(4,beta),mink(4,rho))+4*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,beta)),gamma(in,out,mink(4,rho)))*g(mink(4,alpha),mink(4,sigma))+4*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,rho)),gamma(in,out,mink(4,sigma)))*g(mink(4,alpha),mink(4,beta)))*1/2",
    );
}

#[test]
fn dirac_simplify_id4_slash_sandwich() {
    let r = test_initialize();
    let p = p!(r.mink4);
    let q = q!(r.mink4);
    let expr = Atom::var(s!(m))
        * chain!(slot!(r.bis4, i), slot!(r.bis4, a), gamma!(p.clone()))
        * chain!(slot!(r.bis4, a), slot!(r.bis4, j), gamma!(p.clone()))
        - chain!(slot!(r.bis4, i), slot!(r.bis4, a), gamma!(p.clone()))
            * chain!(slot!(r.bis4, a), slot!(r.bis4, b), gamma!(q.clone()))
            * chain!(slot!(r.bis4, b), slot!(r.bis4, j), gamma!(p.clone()));

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_gamma(crate::dirac::GammaSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().expand().to_bare_ordered_string(), @"-2*chain(bis(4,i),bis(4,j),gamma(in,out,p(mink(4))))*g(p(mink(4)),q(mink(4)))+chain(bis(4,i),bis(4,j),gamma(in,out,q(mink(4))))*g(p(mink(4)),p(mink(4)))+g(bis(4,i),bis(4,j))*g(p(mink(4)),p(mink(4)))*m");
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

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_gamma(crate::dirac::GammaSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"-1*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,mu)),gamma5(in,out))");
}

#[test]
fn dirac_simplify_id23_trace_evaluation_can_stay_disabled() {
    let r = test_initialize();
    let expr = trace!(
        r.bis4.to_symbolic([]),
        gamma!(slot!(r.mink4, a)),
        gamma!(slot!(r.mink4, b)),
        gamma!(slot!(r.mink4, a))
    ) + gamma!(1, 2, a) * gamma!(2, 3, a);

    // These references combine independent tensors with different interfaces.
    // Simplify each summand through its checked boundary; such a sum itself
    // is deliberately not admitted as one tensor.
    let independent = expr;
    assert!(crate::tensor::SymbolicTensor::infer(independent.clone()).is_err());
    let symbolica::atom::AtomView::Add(terms) = independent.as_view() else {
        panic!("the reference must contain distinct tensor summands");
    };
    let simplified = Atom::add_many(terms.iter().map(|term| {
        crate::tensor::SymbolicTensor::infer(term.to_owned())
            .unwrap()
            .simplify_gamma(GammaSimplifySettings::repeated_pairs().without_trace_evaluation())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
    }));
    assert_snapshot!(simplified.to_bare_ordered_string(), @"4*g(bis(4,1),bis(4,3))+trace(bis(4),cyclic(gamma(in,out,mink(4,a)),gamma(in,out,mink(4,a)),gamma(in,out,mink(4,b))))");
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

    // These references combine independent tensors with different interfaces.
    // Simplify each summand through its checked boundary; such a sum itself
    // is deliberately not admitted as one tensor.
    let independent = expr;
    assert!(crate::tensor::SymbolicTensor::infer(independent.clone()).is_err());
    let symbolica::atom::AtomView::Add(terms) = independent.as_view() else {
        panic!("the reference must contain distinct tensor summands");
    };
    let simplified = Atom::add_many(terms.iter().map(|term| {
        crate::tensor::SymbolicTensor::infer(term.to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
    }));
    assert_snapshot!(simplified.to_bare_ordered_string(), @"chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,a)),gamma(in,out,mink(4,b)))");
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

    // These references combine independent tensors with different interfaces.
    // Simplify each summand through its checked boundary; such a sum itself
    // is deliberately not admitted as one tensor.
    let independent = expr;
    assert!(crate::tensor::SymbolicTensor::infer(independent.clone()).is_err());
    let symbolica::atom::AtomView::Add(terms) = independent.as_view() else {
        panic!("the reference must contain distinct tensor summands");
    };
    let simplified = Atom::add_many(terms.iter().map(|term| {
        crate::tensor::SymbolicTensor::infer(term.to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
    }));
    assert_snapshot!(simplified.to_bare_ordered_string(), @"-4*g(mink(4,a),mink(4,c))*g(mink(4,b),mink(4,d))+4*g(mink(4,a),mink(4,b))*g(mink(4,c),mink(4,d))+4*g(mink(4,a),mink(4,d))*g(mink(4,b),mink(4,c))+chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,a)),gamma(in,out,mink(4,b)),gamma(in,out,mink(4,c)),gamma(in,out,mink(4,d)))");
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

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned()).unwrap().simplify_gamma(GammaSimplifySettings::canonical()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"0");
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

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr_4).as_atom_view().to_owned()).unwrap().simplify_gamma(crate::dirac::GammaSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"4");
    assert_snapshot!(crate::tensor::SymbolicTensor::infer((expr_d).as_atom_view().to_owned()).unwrap().simplify_gamma(crate::dirac::GammaSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"d");
    let color = crate::tensor::SymbolicTensor::infer(expr_color.clone()).unwrap();
    assert_eq!(
        color
            .simplify_gamma(GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
        expr_color
    );
    assert_snapshot!(color.simplify_color(crate::color::ColorSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"Nc");
    assert_snapshot!(crate::tensor::SymbolicTensor::infer((pair_d).as_atom_view().to_owned()).unwrap().simplify_gamma(crate::dirac::GammaSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"d*g(mink(d,a),mink(d,b))");
}

#[test]
fn dirac_simplify_id40_repeated_gamma_and_two_trace() {
    let r = test_initialize();
    let p = p!(&r.mink4);
    let open = Atom::var(s!(c1))
        * gamma!(slot!(r.bis4, i), slot!(r.bis4, a), mu)
        * chain!(slot!(r.bis4, a), slot!(r.bis4, b), gamma!(p.clone()))
        * gamma!(slot!(r.bis4, b), slot!(r.bis4, j), mu)
        + Atom::var(s!(c1))
            * Atom::var(s!(m))
            * gamma!(slot!(r.bis4, i), slot!(r.bis4, a), mu)
            * gamma!(slot!(r.bis4, a), slot!(r.bis4, j), mu);
    let tr = Atom::var(s!(c2))
        * trace!(
            r.bis4.to_symbolic([]),
            gamma!(slot!(r.mink4, mu)),
            gamma!(slot!(r.mink4, nu))
        );

    // These references combine independent tensors with different interfaces.
    // Simplify each summand through its checked boundary; such a sum itself
    // is deliberately not admitted as one tensor.
    let independent = open + tr;
    assert!(crate::tensor::SymbolicTensor::infer(independent.clone()).is_err());
    let symbolica::atom::AtomView::Add(terms) = independent.as_view() else {
        panic!("the reference must contain distinct tensor summands");
    };
    let simplified = Atom::add_many(terms.iter().map(|term| {
        crate::tensor::SymbolicTensor::infer(term.to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression()
    }));
    assert_snapshot!(simplified.expand().to_bare_ordered_string(), @"-2*c1*chain(bis(4,i),bis(4,j),gamma(in,out,p(mink(4))))+4*c1*g(bis(4,i),bis(4,j))*m+4*c2*g(mink(4,mu),mink(4,nu))");
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

    assert_snapshot!(crate::tensor::SymbolicTensor::infer((slash_square).as_atom_view().to_owned()).unwrap().simplify_gamma(crate::dirac::GammaSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().to_bare_ordered_string(), @"chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,mu)))*g(p(mink(4)),p(mink(4)))");
    assert_snapshot!(crate::tensor::SymbolicTensor::infer((sandwich).as_atom_view().to_owned()).unwrap().simplify_gamma(crate::dirac::GammaSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().expand().to_bare_ordered_string(), @"-1*chain(bis(4,i),bis(4,j),gamma(in,out,mink(4,mu)))*g(p(mink(4)),p(mink(4)))+2*chain(bis(4,i),bis(4,j),gamma(in,out,p(mink(4))))*p(mink(4,mu))");
}

#[test]
fn dirac_simplify_id46_d_dim_repeated_gamma_slash_sum() {
    let r = test_initialize();
    let p = p!(&r.mink_d);
    let q = q!(&r.mink_d);
    let expr = gamma!(slot!(r.bis_d, i), slot!(r.bis_d, a), slot!(r.mink_d, mu))
        * chain!(slot!(r.bis_d, a), slot!(r.bis_d, b), gamma!(p.clone()))
        * gamma!(slot!(r.bis_d, b), slot!(r.bis_d, j), slot!(r.mink_d, mu))
        + gamma!(slot!(r.bis_d, i), slot!(r.bis_d, a), slot!(r.mink_d, mu))
            * chain!(slot!(r.bis_d, a), slot!(r.bis_d, b), gamma!(q.clone()))
            * gamma!(slot!(r.bis_d, b), slot!(r.bis_d, j), slot!(r.mink_d, mu))
        + Atom::var(s!(m))
            * gamma!(slot!(r.bis_d, i), slot!(r.bis_d, a), slot!(r.mink_d, mu))
            * gamma!(slot!(r.bis_d, a), slot!(r.bis_d, j), slot!(r.mink_d, mu));

    // Formal spin channels are admitted by the existing word boundary;
    // standalone gamma factories retain their four-dimensional spin ports.
    let chainified = expr.chainify(Bispinor {}.into());
    // The source already contains slash chains; compare both expressions at
    // the same explicit-matrix boundary without expanding their products.
    assert_eq!(
        chainified.undo_chain::<AbstractIndex>().unwrap(),
        expr.undo_chain::<AbstractIndex>().unwrap()
    );
    assert_snapshot!(crate::tensor::SymbolicTensor::infer(chainified).unwrap().simplify_gamma(crate::dirac::GammaSimplifySettings::default()).unwrap().resolved().unwrap().into_expression().expand().to_bare_ordered_string(), @"-1*chain(bis(d,i),bis(d,j),gamma(in,out,p(mink(d))))*d+-1*chain(bis(d,i),bis(d,j),gamma(in,out,q(mink(d))))*d+2*chain(bis(d,i),bis(d,j),gamma(in,out,p(mink(d))))+2*chain(bis(d,i),bis(d,j),gamma(in,out,q(mink(d))))+d*g(bis(d,i),bis(d,j))*m");
}

#[test]
fn chiral_projector_traces() {
    let r = test_initialize();
    let rep = r.bis4.to_symbolic([]);
    let mu = slot!(r.mink4, mu);
    let nu = slot!(r.mink4, nu);
    let rho = slot!(r.mink4, rho);
    let sigma = slot!(r.mink4, sigma);
    for (symbol, sign) in [(AGS.projp, 1), (AGS.projm, -1)] {
        let projector = function!(symbol, T.chain_in, T.chain_out);
        assert_eq!(
            crate::tensor::SymbolicTensor::infer(
                (trace!(rep.clone(), projector.clone()))
                    .as_atom_view()
                    .to_owned()
            )
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression(),
            Atom::num(2)
        );
        let expr = trace!(rep.clone(), projector.clone(), gamma!(mu), gamma!(nu));
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            2 * g!(mu, nu)
        );

        let ordinary = trace!(
            rep.clone(),
            gamma!(mu),
            gamma!(nu),
            gamma!(rho),
            gamma!(sigma)
        );
        let axial = trace!(
            rep.clone(),
            gamma5!(),
            gamma!(mu),
            gamma!(nu),
            gamma!(rho),
            gamma!(sigma)
        );
        let expr = trace!(
            rep.clone(),
            projector.clone(),
            gamma!(mu),
            gamma!(nu),
            gamma!(rho),
            gamma!(sigma)
        );
        let reduced = crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
            .unwrap()
            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
            .unwrap()
            .resolved()
            .unwrap()
            .into_expression();
        assert_eq!(
            (reduced.clone()
                - (crate::tensor::SymbolicTensor::infer((ordinary).as_atom_view().to_owned())
                    .unwrap()
                    .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                    .unwrap()
                    .resolved()
                    .unwrap()
                    .into_expression()
                    + Atom::num(sign)
                        * crate::tensor::SymbolicTensor::infer((axial).as_atom_view().to_owned())
                            .unwrap()
                            .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                            .unwrap()
                            .resolved()
                            .unwrap()
                            .into_expression())
                    / 2)
            .expand(),
            Atom::Zero
        );
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((reduced).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            reduced
        );
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(GammaSimplifySettings::default().without_trace_evaluation())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            expr
        );

        // A gamma flips chirality, so equal projectors separated by one gamma vanish.
        let expr = trace!(
            rep.clone(),
            projector.clone(),
            gamma!(mu),
            projector.clone(),
            gamma!(nu)
        );
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            Atom::Zero
        );
        let expr = trace!(rep.clone(), projector.clone(), projector.clone());
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((expr).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            Atom::num(2)
        );
        let mixed = trace!(rep.clone(), projector.clone(), gamma!(slot!(r.mink_d, mu)));
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((mixed).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            mixed
        );
        let dimensional = trace!(r.bis_d.to_symbolic([]), projector);
        assert_eq!(
            crate::tensor::SymbolicTensor::infer((dimensional).as_atom_view().to_owned())
                .unwrap()
                .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
                .unwrap()
                .resolved()
                .unwrap()
                .into_expression(),
            dimensional
        );
    }
    let plus = function!(AGS.projp, T.chain_in, T.chain_out);
    let minus = function!(AGS.projm, T.chain_in, T.chain_out);
    assert_eq!(
        crate::tensor::SymbolicTensor::infer(
            (trace!(rep.clone(), plus.clone(), minus.clone()))
                .as_atom_view()
                .to_owned()
        )
        .unwrap()
        .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression(),
        Atom::Zero
    );
    assert_eq!(
        crate::tensor::SymbolicTensor::infer(
            (trace!(rep, plus, gamma!(mu), minus, gamma!(nu)))
                .as_atom_view()
                .to_owned()
        )
        .unwrap()
        .simplify_gamma(crate::dirac::GammaSimplifySettings::default())
        .unwrap()
        .resolved()
        .unwrap()
        .into_expression(),
        2 * g!(mu, nu)
    );
}
