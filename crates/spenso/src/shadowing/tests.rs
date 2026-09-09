use crate::structure::representation::{LibraryRep, Minkowski};
#[allow(unused_imports)]
use crate::{
    aind, antisym, bracket, chain, chain_factor, cyclic, dind, dot, euc, g, lor, mink,
    network::tags::SPENSO_TAG,
    p, pure_scalar, q,
    shadowing::{ProjectorExpander, TensorCollectExt},
    sym, trace, trace_sym,
};
use ::symbolica_utils::AtomPrintExt;
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
fn collect_tensors_keeps_shared_coefficient_outside_tensor_sum() {
    let (a, b, c) = symbol!(
        "shared_coefficient_a",
        "shared_coefficient_b",
        "shared_coefficient_c"
    );
    let coefficient = (Atom::var(a) + Atom::var(b)).pow(9) * (Atom::var(a) + Atom::var(c));
    let expression = coefficient * (p!(mink!(4)) + q!(mink!(4)));

    let collected = expression.collect_tensors();
    assert_eq!(collected, expression);
    assert_eq!(collected.collect_tensors(), collected);
}

#[test]
fn collect_rep_keeps_unrelated_traces_outside_tensor_sum() {
    let coefficient = trace_sym!(
        euc!(3),
        chain_factor!(collect_spectator_a, in, out),
        chain_factor!(collect_spectator_b, in, out),
    ) + trace!(euc!(3), chain_factor!(collect_spectator_c, in, out));
    let expression = coefficient * (p!(mink!(4)) + q!(mink!(4)));
    let rep = LibraryRep::from(Minkowski {});

    let collected = expression.collect_rep(rep);
    assert_eq!(collected, expression);
    assert_eq!(collected.collect_rep(rep), collected);
}

#[test]
fn collect_tensors_preserves_symbolic_powers_and_repeated_coefficients() {
    let (a, b, c) = symbol!("coefficient_a", "coefficient_b", "coefficient_c");
    let sum = Atom::var(a) + Atom::var(b);
    let coefficient = Atom::var(c).pow(sum.clone()) * sum.pow(3);
    let expression = coefficient.clone() * p!(mink!(4)) + coefficient * q!(mink!(4));

    let collected = expression.collect_tensors();
    assert_eq!(collected, expression);
    assert_eq!(collected.collect_tensors(), collected);
}

#[test]
fn collect_tensors_coefficient_aliases_do_not_capture_input_symbols() {
    let (reserved, a, b) = symbol!(
        "spenso::collect_coefficient_0",
        "coefficient_capture_a",
        "coefficient_capture_b"
    );
    let coefficient = (Atom::var(reserved) + Atom::var(a)).pow(3) * (Atom::var(a) + Atom::var(b));
    let expression = coefficient * p!(mink!(4));

    let collected = expression.collect_tensors();
    assert_eq!(collected, expression);
    assert!(collected.contains_symbol(reserved));
}

#[test]
fn collect_tensors_restores_linear_coefficient_cancellation() {
    let (a, b) = symbol!("coefficient_cancel_a", "coefficient_cancel_b");
    let tensor = p!(mink!(4));
    let expression =
        (Atom::var(a) + Atom::var(b)) * &tensor - Atom::var(a) * &tensor - Atom::var(b) * tensor;

    assert!(expression.collect_tensors().is_zero());
}

#[test]
fn collect_tensors_keeps_opaque_polynomial_coefficients_intact() {
    let (a, b) = symbol!("coefficient_zero_a", "coefficient_zero_b");
    let a = Atom::var(a);
    let b = Atom::var(b);
    let coefficient = (&a + &b).pow(2) - a.pow(2) - Atom::num(2) * &a * &b - b.pow(2);
    let expression = coefficient * p!(mink!(4));

    // Tensor collection groups tensor factors; it does not expand a protected
    // coefficient merely to prove a polynomial identity inside that factor.
    assert_eq!(expression.collect_tensors(), expression);
}

#[test]
fn collect_rep_callback_receives_complete_unaliased_tensor_payload() {
    use crate::shadowing::Collectable;
    use symbolica::atom::FunctionBuilder;

    let (a, b, c, target) = symbol!(
        "coefficient_callback_a",
        "coefficient_callback_b",
        "coefficient_callback_c",
        "coefficient_callback_tensor"
    );
    let sum = Atom::var(a) + Atom::var(b);
    let tensor = FunctionBuilder::new(target)
        .add_arg(sum.clone().pow(3))
        .add_arg(mink!(4, mu))
        .finish();
    let coefficient = Atom::var(c).pow(sum);
    let expression = coefficient.clone() * &tensor + Atom::var(a) * &tensor;
    let replacement = p!(mink!(4));
    let collected =
        expression.collect_rep_with_map(LibraryRep::from(Minkowski {}), |wrapped, _, out| {
            assert_eq!(wrapped, tensor.clone().wrap_in_collect().as_view());
            **out = replacement.clone();
        });

    assert_eq!(collected, (coefficient + Atom::var(a)) * replacement);
}

#[test]
fn collect_rep_callback_receives_complete_contractions_across_tensor_sum() {
    use crate::shadowing::Collectable;

    let (a, b) = symbol!("contracted_coefficient_a", "contracted_coefficient_b");
    let coefficient = (Atom::var(a) + Atom::var(b)).pow(3);
    let p = p!(mink!(4, mu));
    let q = q!(mink!(4, mu));
    let r = symbol!("contracted_callback_r").call(mink!(4, mu));
    let pq = (&p * &q).wrap_in_collect();
    let pr = (&p * &r).wrap_in_collect();
    let expression = &coefficient * p * (q + r);
    let mut seen = [0, 0];

    let result =
        expression.collect_rep_with_map(LibraryRep::from(Minkowski {}), |wrapped, _, out| {
            if wrapped == pq.as_view() {
                seen[0] += 1;
                **out = Atom::num(7);
            } else {
                assert_eq!(wrapped, pr.as_view());
                seen[1] += 1;
                **out = Atom::num(11);
            }
        });

    assert_eq!(seen, [1, 1]);
    assert_eq!(result, Atom::num(18) * coefficient);
}

#[test]
fn collect_rep_callback_receives_complete_compressed_tensor_powers() {
    use crate::shadowing::Collectable;

    let (a, b) = symbol!("powered_coefficient_a", "powered_coefficient_b");
    let coefficient = (Atom::var(a) + Atom::var(b)).pow(9);
    let tensor = p!(mink!(4, mu));
    let power = tensor.pow(3);
    let result = (coefficient.clone() * &power).collect_rep_with_map(
        LibraryRep::from(Minkowski {}),
        |wrapped, _, out| {
            assert_eq!(wrapped, power.clone().wrap_in_collect().as_view());
            **out = Atom::num(7);
        },
    );
    assert_eq!(result, Atom::num(7) * coefficient);
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
fn representation_collection_scan_matches_pattern_oracle() {
    use crate::{
        broadcast_symbol,
        shadowing::{collect::TensorCollectFilter, static_symbols::W_},
        structure::representation::Euclidean,
    };
    use symbolica::{
        atom::{FunctionBuilder, representation::FunView},
        function,
    };
    use symbolica_utils::ReplaceBuilderExt;

    // Retain the original matcher as an independent oracle for the direct scan,
    // including guards that apply only at the outer function boundary.
    fn pattern_matches(fun: FunView<'_>, reps: &[LibraryRep]) -> bool {
        let symbol = fun.get_symbol();
        if symbol == SPENSO_TAG.pure_scalar || symbol == SPENSO_TAG.bracket {
            return false;
        }
        if symbol.has_tag(&SPENSO_TAG.broadcast) {
            let args = fun.iter().collect::<Vec<_>>();
            return matches!(args.as_slice(), [AtomView::Fun(arg)] if pattern_matches(*arg, reps));
        }
        for arg in fun.iter() {
            for rep in reps {
                if arg.replace(function!(rep.symbol(), W_.a__)).matches() {
                    return true;
                }
            }
        }
        if symbol == SPENSO_TAG.chain {
            return fun
                .iter()
                .skip(2)
                .any(|arg| matches!(arg, AtomView::Fun(arg) if pattern_matches(arg, reps)));
        }
        if symbol == SPENSO_TAG.trace {
            return fun
                .iter()
                .skip(1)
                .any(|arg| matches!(arg, AtomView::Fun(arg) if pattern_matches(arg, reps)));
        }
        false
    }

    let mink = LibraryRep::from(Minkowski {});
    let euc = LibraryRep::from(Euclidean {});
    let (outer, nested, x) = symbol!("scan_outer", "scan_nested", "scan_x");
    let mink_arg = function!(mink.symbol(), 4);
    let payloads = [
        Atom::num(1),
        Atom::var(mink.symbol()),
        FunctionBuilder::new(mink.symbol()).finish(),
        mink_arg.clone(),
        function!(mink.symbol(), 4, Atom::var(x)),
        function!(euc.symbol(), 3),
        function!(nested, mink_arg.clone()),
        function!(SPENSO_TAG.pure_scalar, mink_arg.clone()),
        function!(SPENSO_TAG.bracket, mink_arg.clone()),
        mink_arg.clone().pow(2),
        mink_arg.clone().pow(-1),
        Atom::var(x).pow(mink_arg.clone()),
        mink_arg.clone() + Atom::var(x),
        mink_arg * Atom::var(x),
    ];
    let heads = [
        outer,
        SPENSO_TAG.pure_scalar,
        SPENSO_TAG.bracket,
        broadcast_symbol!(scan_broadcast),
        SPENSO_TAG.chain,
        SPENSO_TAG.trace,
    ];
    for head in heads {
        for payload in &payloads {
            for arity in 0..=3 {
                let mut builder = FunctionBuilder::new(head);
                for _ in 0..arity {
                    builder = builder.add_arg(payload);
                }
                let expression = builder.finish();
                // Linear heads may normalize to sums or scalar factors. Check
                // every resulting function, including the normalized branches.
                expression.visitor(&mut |atom| {
                    if let AtomView::Fun(fun) = atom {
                        for reps in [&[][..], &[mink][..], &[euc][..], &[mink, euc][..]] {
                            assert_eq!(
                                TensorCollectFilter::<0>::function_contains_rep(fun, reps),
                                pattern_matches(fun, reps),
                                "representation scan differs for {atom} in {expression} and {reps:?}",
                            );
                        }
                    }
                    true
                });
            }
        }
    }
}

#[test]
fn expand_rep_with_map_visits_expanded_collect_wrappers() {
    let (a, b, mapped_tensor) = symbol!("a", "b", "mapped_tensor");
    let expr = (Atom::var(a) + Atom::var(b)) * p!(mink!(4));

    let expanded = expr.expand_rep_with_map(LibraryRep::from(Minkowski {}), |_arg, _ctx, out| {
        **out = Atom::var(mapped_tensor);
    });

    insta::assert_snapshot!(
        expanded.to_bare_ordered_string(),
        @"(a+b)*mapped_tensor"
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
