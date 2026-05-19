use crate::structure::representation::{LibraryRep, Minkowski};
#[allow(unused_imports)]
use crate::{
    aind, antisym, bracket, chain, chain_factor, cyclic, dind, dot, euc, g, lor, mink,
    network::tags::SPENSO_TAG,
    p, pure_scalar, q,
    shadowing::symbolica_utils::AtomCoreExt,
    shadowing::{ProjectorExpander, TensorCollectExt},
    sym, trace, trace_sym,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    symbol,
};

#[test]
fn tensor_macros_create_tagged_heads() {
    #[allow(unused_imports)]
    use crate::{tensor, vector};

    let tensor = tensor!(tensor_macro_test_head);
    let vector = vector!(vector_macro_test_head, 1);

    let AtomView::Fun(tensor) = tensor.as_view() else {
        panic!("tensor macro should produce a function");
    };
    assert!(tensor.get_symbol().has_tag(&SPENSO_TAG.tensor));
    assert!(!tensor.get_symbol().has_tag(&SPENSO_TAG.rank1));

    let AtomView::Fun(vector) = vector.as_view() else {
        panic!("vector macro should produce a function");
    };
    assert!(vector.get_symbol().has_tag(&SPENSO_TAG.tensor));
    assert!(vector.get_symbol().has_tag(&SPENSO_TAG.rank1));
}

#[test]
fn surface_macros_build_parser_syntax() {
    let mu = mink!(4, mu);
    let nu = mink!(4, nu);
    let p = p!(mink!(4));
    let q = q!(mink!(4));

    insta::assert_snapshot!(g!(mu.clone(), nu.clone()).to_bare_ordered_string(), @"g(mink(4,mu),mink(4,nu))");
    insta::assert_snapshot!(dot!(p.clone(), q.clone()).to_bare_ordered_string(), @"dot(p(mink(4)),q(mink(4)))");
    insta::assert_snapshot!(aind!(mu.clone(), nu.clone()).to_bare_ordered_string(), @"aind(mink(4,mu),mink(4,nu))");
    let dual = dind!(lor!(4, rho));
    let AtomView::Fun(dual) = dual.as_view() else {
        panic!("dind macro should produce a function");
    };
    assert_eq!(
        dual.get_symbol(),
        crate::structure::abstract_index::AIND_SYMBOLS.dind
    );
    insta::assert_snapshot!(pure_scalar!(Atom::num(1)).to_bare_ordered_string(), @"pure_scalar(1)");
    insta::assert_snapshot!(bracket!(p, q).to_bare_ordered_string(), @"bracket(p(mink(4)),q(mink(4)))");
    insta::assert_snapshot!(
        chain_factor!(factor, in, mu, out).to_bare_ordered_string(),
        @"factor(in,mink(4,mu),out)"
    );
}

#[test]
fn collect_tensors_keeps_scalar_products_factored() {
    let (a, b, c, d) = symbol!("a", "b", "c", "d");
    let p = p!(mink!(4));
    let expr = (Atom::var(a) + Atom::var(b)) * (Atom::var(c) + Atom::var(d)) * p.clone()
        + Atom::var(a) * p;

    insta::assert_snapshot!(
        expr.collect_tensors().to_bare_ordered_string(),
        @"((a+b)*(c+d)+a)*p(mink(4))"
    );
}

#[test]
fn collect_tensors_marks_chain_like_forms_as_maximal_factors() {
    let (a, b) = symbol!("a", "b");
    let mu = mink!(4, mu);
    let gamma_mu = chain_factor!(collect_test_gamma, in, mu, out);

    let chain_expr = chain!(mink!(4, i), mink!(4, j), gamma_mu.clone());
    insta::assert_snapshot!(
        (Atom::var(a) * chain_expr.clone() + Atom::var(b) * chain_expr)
            .collect_tensors()
            .to_bare_ordered_string(),
        @"(a+b)*chain(mink(4,i),mink(4,j),collect_test_gamma(in,mink(4,mu),out))"
    );

    let trace_expr = trace!(mink!(4), gamma_mu.clone());
    insta::assert_snapshot!(
        (Atom::var(a) * trace_expr.clone() + Atom::var(b) * trace_expr)
            .collect_tensors()
            .to_bare_ordered_string(),
        @"(a+b)*trace(mink(4),cyclic(collect_test_gamma(in,mink(4,mu),out)))"
    );

    let sym_expr = sym!(gamma_mu.clone());
    insta::assert_snapshot!(
        (Atom::var(a) * sym_expr.clone() + Atom::var(b) * sym_expr)
            .collect_tensors()
            .to_bare_ordered_string(),
        @"(a+b)*sym(collect_test_gamma(in,mink(4,mu),out))"
    );

    let cyclic_expr = cyclic!(gamma_mu);
    insta::assert_snapshot!(
        (Atom::var(a) * cyclic_expr.clone() + Atom::var(b) * cyclic_expr)
            .collect_tensors()
            .to_bare_ordered_string(),
        @"(a+b)*cyclic(collect_test_gamma(in,mink(4,mu),out))"
    );
}

#[test]
fn collect_rep_only_wraps_matching_representations() {
    let (a, b) = symbol!("a", "b");
    let mink_p = p!(mink!(4));
    let euc_q = q!(euc!(3));
    let expr = Atom::var(a) * mink_p.clone()
        + Atom::var(b) * mink_p
        + Atom::var(a) * euc_q.clone()
        + Atom::var(b) * euc_q;

    insta::assert_snapshot!(
        expr.collect_rep(LibraryRep::from(Minkowski {}))
            .to_bare_ordered_string(),
        @"(a+b)*p(mink(4))+a*q(euc(3))+b*q(euc(3))"
    );
}

#[test]
fn collect_metrics_only_wraps_metric_heads() {
    let (a, b) = symbol!("a", "b");
    let metric = g!(mink!(4, mu), mink!(4, nu));
    let vector = p!(mink!(4));
    let expr = Atom::var(a) * metric.clone()
        + Atom::var(b) * metric
        + Atom::var(a) * vector.clone()
        + Atom::var(b) * vector;

    insta::assert_snapshot!(
        expr.collect_metrics().to_bare_ordered_string(),
        @"(a+b)*g(mink(4,mu),mink(4,nu))+a*p(mink(4))+b*p(mink(4))"
    );
}

#[test]
fn chain_macro_accepts_iterable_factors() {
    let start = Atom::var(symbol!("start"));
    let end = Atom::var(symbol!("end"));
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));

    assert_eq!(
        chain!(start.clone(), end.clone(); vec![first.clone(), second.clone()]),
        chain!(start, end, first, second)
    );
}

#[test]
fn trace_macro_accepts_iterable_factors() {
    let rep = Atom::var(symbol!("rep"));
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));

    assert_eq!(
        trace!(rep.clone(); vec![first.clone(), second.clone()]),
        trace!(rep.clone(), first.clone(), second.clone())
    );
    assert_eq!(
        trace!(rep, first.clone(), second.clone()),
        trace!(Atom::var(symbol!("rep")), cyclic!(first, second))
    );
}

#[test]
fn sym_macro_canonicalizes_arguments() {
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));

    assert_eq!(sym!(second.clone(), first.clone()), sym!(first, second));
}

#[test]
fn cyclic_macro_canonicalizes_rotations() {
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));
    let third = Atom::var(symbol!("third"));

    assert_eq!(
        cyclic!(second.clone(), third.clone(), first.clone()),
        cyclic!(first.clone(), second.clone(), third.clone())
    );
    assert_ne!(
        cyclic!(first.clone(), third.clone(), second.clone()),
        cyclic!(second, third, first)
    );
}

#[test]
fn cyclic_unwraps_single_symmetric_projector() {
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));

    assert_eq!(
        cyclic!(sym!(first.clone(), second.clone())),
        sym!(first, second)
    );
}

#[test]
fn trace_sym_macro_builds_symmetric_trace() {
    let rep = Atom::var(symbol!("rep"));
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));

    assert_eq!(
        trace_sym!(rep.clone(), first.clone(), second.clone()),
        trace!(rep, cyclic!(sym!(first, second)))
    );
}

#[test]
fn expand_sym_in_chain_preserves_endpoints() {
    let start = Atom::var(symbol!("start"));
    let end = Atom::var(symbol!("end"));
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));

    let expanded = chain!(
        start.clone(),
        end.clone(),
        sym!(first.clone(), second.clone())
    )
    .expand_projectors();
    let expected = Atom::num(1) / Atom::num(2)
        * chain!(start.clone(), end.clone(), first.clone(), second.clone())
        + Atom::num(1) / Atom::num(2) * chain!(start, end, second, first);

    assert_eq!(expanded, expected);
}

#[test]
fn expand_antisym_in_trace_uses_cyclicity() {
    let rep = Atom::var(symbol!("rep"));
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));

    let expanded = trace!(rep.clone(), antisym!(first.clone(), second.clone())).expand_projectors();

    assert_eq!(expanded, Atom::Zero);
}

#[test]
fn expand_sym_in_trace_emits_cyclic_projectors() {
    let rep = Atom::var(symbol!("rep"));
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));
    let third = Atom::var(symbol!("third"));

    let expanded = trace!(
        rep.clone(),
        sym!(first.clone(), second.clone(), third.clone())
    )
    .expand_projectors();
    let expected = Atom::num(1) / Atom::num(2)
        * trace!(
            rep.clone(),
            cyclic!(first.clone(), second.clone(), third.clone())
        )
        + Atom::num(1) / Atom::num(2) * trace!(rep, cyclic!(first, third, second));

    assert_eq!(expanded, expected);
}

#[test]
fn expand_cyclic_in_chain_rotates_factors() {
    let start = Atom::var(symbol!("start"));
    let end = Atom::var(symbol!("end"));
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));
    let third = Atom::var(symbol!("third"));

    let expanded = chain!(
        start.clone(),
        end.clone(),
        cyclic!(first.clone(), second.clone(), third.clone())
    )
    .expand_projectors();
    let expected = Atom::num(1) / Atom::num(3)
        * chain!(
            start.clone(),
            end.clone(),
            first.clone(),
            second.clone(),
            third.clone()
        )
        + Atom::num(1) / Atom::num(3)
            * chain!(
                start.clone(),
                end.clone(),
                second.clone(),
                third.clone(),
                first.clone()
            )
        + Atom::num(1) / Atom::num(3) * chain!(start, end, third, first, second);

    assert_eq!(expanded, expected);
}

#[test]
fn expand_cyclic_in_trace_keeps_compact_cyclic_trace() {
    let rep = Atom::var(symbol!("rep"));
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));
    let third = Atom::var(symbol!("third"));

    let expanded = trace!(
        rep.clone(),
        cyclic!(first.clone(), second.clone(), third.clone())
    )
    .expand_projectors();
    let expected = trace!(rep, cyclic!(first, second, third));

    assert_eq!(expanded, expected);
}

#[test]
fn expand_antisym_in_chain_handles_canonicalization_sign() {
    let start = Atom::var(symbol!("start"));
    let end = Atom::var(symbol!("end"));
    let first = Atom::var(symbol!("first"));
    let second = Atom::var(symbol!("second"));

    let expanded = chain!(
        start.clone(),
        end.clone(),
        antisym!(second.clone(), first.clone())
    )
    .expand_projectors();
    let expected = Atom::num(-1) / Atom::num(2)
        * chain!(start.clone(), end.clone(), first.clone(), second.clone())
        + Atom::num(1) / Atom::num(2) * chain!(start, end, second, first);

    assert_eq!(expanded.expand(), expected.expand());
}
