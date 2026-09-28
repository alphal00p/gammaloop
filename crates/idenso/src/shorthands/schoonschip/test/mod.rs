use insta::assert_snapshot;
use spenso::{
    chain, g, mink,
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    p, q, slot,
    structure::{
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        representation::{Minkowski, RepName, Representation},
        slot::IsAbstractSlot,
    },
    trace, trace_sym,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    function, symbol,
};
use symbolica_utils::AtomPrintExt;

use crate::{
    representations::{Bispinor, ColorFundamental},
    test_support::{contracted_atom, test_initialize},
};

use super::{Schoonschip, SchoonschipSettings};

fn chain_in() -> Atom {
    Atom::var(T.chain_in)
}

fn chain_out() -> Atom {
    Atom::var(T.chain_out)
}

pub(super) fn is_slot(atom: &Atom) -> bool {
    is_slot_view(atom.as_view())
}

fn is_slot_view(atom: AtomView<'_>) -> bool {
    let AtomView::Fun(f) = atom else {
        return false;
    };

    if f.get_symbol().has_tag(&T.representation) {
        return f.get_nargs() == 2;
    }

    // Downstairs indices wrap the actual slot as `dind(slot)`.
    f.get_symbol() == AIND_SYMBOLS.dind
        && f.get_nargs() == 1
        && f.iter().next().is_some_and(is_slot_view)
}

#[test]
fn tensor_sum_square_uses_local_contractions() {
    test_initialize();
    let p = p!(1, mink!(4, 1));
    let q = q!(2, mink!(4, 1));
    let x = Atom::var(symbol!("square_spectator_x"));
    let y = Atom::var(symbol!("square_spectator_y"));
    let spectator = (x + y).pow(12);
    let square = &p * &p + Atom::num(2) * &p * &q + &q * &q;
    let expression = &spectator * ((p + q).pow(2) - square + Atom::num(1));
    let result = contracted_atom(expression.as_view()).unwrap().to_dots();
    assert!((result - spectator).expand().is_zero());
}

#[test]
fn network_sum_preserves_factored_spectators_logical_order_and_zero() {
    use crate::tensor::SymbolicTensor;
    use spenso::structure::partial::{PartialStructure, PartialStructureExt};

    test_initialize();
    let a = mink!(4, 91201);
    let b = mink!(4, 91202);
    let c = mink!(4, 91203);
    let tensor = spenso::tensor_symbol!("network_sum_ordered_tensor");
    let x = Atom::var(symbol!("network_sum_x"));
    let y = Atom::var(symbol!("network_sum_y"));
    let spectator = (x + y).pow(7);
    let contracted = &spectator * function!(tensor, &c, &b);
    let product = &spectator * g!(&a, &b) * function!(tensor, &c, &a);
    for (source, expected) in [
        (&product + &contracted, Atom::num(2) * &contracted),
        (&product - &contracted, Atom::Zero),
    ] {
        let mut typed = SymbolicTensor::<PartialStructure>::infer(source.clone()).unwrap();
        typed.structure =
            PartialStructure::from_logical_slots(typed.structure.logical_slots().into_iter().rev());
        assert_eq!(typed.structure.logical_slots().len(), 2);
        let result = contracted_atom(source.as_view()).unwrap();
        assert_eq!(result, expected);
        let result = typed.with_rewritten_expression(result).unwrap();
        assert_eq!(result.structure, typed.structure);
        assert_eq!(
            contracted_atom(result.expression.as_view()).unwrap(),
            expected
        );
    }
}

#[test]
fn network_sum_matches_ordered_single_term_callbacks() {
    use std::sync::{Arc, Mutex};

    test_initialize();
    let calls = Arc::new(Mutex::new(Vec::new()));
    let observed = Arc::clone(&calls);
    let tensor = spenso::tensor_symbol!(
        "network_sum_callback_tensor",
        norm = move |value, _| observed.lock().unwrap().push(value.to_owned())
    );
    let a = mink!(4, 91301);
    let b = mink!(4, 91302);
    let source = Atom::add_many((0..8).map(|i| g!(&a, &b) * function!(tensor, i, &a)));
    let AtomView::Add(sum) = source.as_view() else {
        panic!("the fixture must exercise whole-sum accumulation");
    };
    calls.lock().unwrap().clear();
    // This is the former driver's ordered prefix accumulation, with the same
    // single-pass reduction for each term.
    let mut expected = Atom::Zero;
    for term in sum.iter() {
        expected += contracted_atom(term).unwrap();
    }
    let transcript = calls.lock().unwrap().clone();
    assert!(!transcript.is_empty());
    calls.lock().unwrap().clear();
    let actual = contracted_atom(source.as_view()).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(*calls.lock().unwrap(), transcript);
}

#[test]
fn simple_dot() {
    test_initialize();
    let dim = symbol!("D");
    let mink: Representation<_> = Minkowski {}.new_rep(dim);
    let mink_d = mink!(dim);

    let p1 = p!(1, slot!(mink, 1));
    let p2 = p!(2, slot!(mink, 1));
    let p1_2 = p!(1, slot!(mink, 2));
    let p2_2 = p!(2, slot!(mink, 2));
    let p1_stripped = p!(1, mink_d.clone());
    let p2_stripped = p!(2, mink_d.clone());

    let q2 = q!(2, symbol!("bla"), slot!(mink, 1));
    let q2_2 = q!(2, symbol!("bla"), slot!(mink, 2));

    let q3 = q!(3, slot!(mink, 1));
    let q3_2 = q!(3, slot!(mink, 2));

    let result = contracted_atom((&p1 * &q2).as_view()).unwrap();
    assert_snapshot!(result.to_bare_ordered_string(),@"g(p(1,mink(D)),q(2,bla,mink(D)))");

    let result = (&p1 * &p1).schoonschip();
    assert_snapshot!(result.to_bare_ordered_string(), @"g(p(1,mink(D)),p(1,mink(D)))");

    let result = p!(1, &p2_stripped).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"g(p(1,mink(dim)),p(2,mink(dim)))");

    let result = g!(slot!(mink, 1), &p1_stripped).normalize_dots();
    // `mink!(dim)` names the literal dimension `dim`, whereas `mink` above uses
    // `D`. A mismatched compact slot must not acquire another explicit slot.
    assert_snapshot!(result.to_bare_ordered_string(), @"g(mink(D,1),p(1,mink(dim)))");

    let result = p1.clone().pow(Atom::num(4)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"(g(p(1,mink(D)),p(1,mink(D))))^2");

    let result = p1.clone().pow(Atom::num(3)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"g(p(1,mink(D)),p(1,mink(D)))*p(1,mink(D,1))");

    let metric = g!(slot!(mink, 1), slot!(mink, 2));
    let result = metric.clone().pow(Atom::num(4)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"D^2");

    let result = metric.pow(Atom::num(3)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"D*g(mink(D,1),mink(D,2))");

    let result = contracted_atom((&p1 * (&q2 + &p2 * &q3_2 * &q2_2)).as_view())
        .unwrap()
        .expand();
    assert_snapshot!(result.to_bare_ordered_string(),@"g(p(1,mink(D)),p(2,mink(D)))*g(q(2,bla,mink(D)),q(3,mink(D)))+g(p(1,mink(D)),q(2,bla,mink(D)))");

    let nested = &p1 * (&q2 + &p2 * (&q3_2 * &q2_2 + &p2_2 * &q2_2));
    let result = contracted_atom(nested.as_view()).unwrap();
    let rep = mink.to_symbolic([]);
    let p1_compact = p!(1, rep.clone());
    let p2_compact = p!(2, rep.clone());
    let q2_compact = q!(2, symbol!("bla"), rep.clone());
    let q3_compact = q!(3, rep);
    let expected = (g!(&p2_compact, &q2_compact) + g!(&q2_compact, &q3_compact))
        * g!(&p1_compact, &p2_compact)
        + g!(&p1_compact, &q2_compact);
    assert!((result - expected).expand().is_zero());

    let expr = (p1 + q3 * p1_2 * q2_2) * (q2 + p2);

    let result = contracted_atom(expr.as_view()).unwrap().expand();
    let result = result.to_bare_ordered_string();
    assert_snapshot!(result, @"g(p(1,mink(D)),p(2,mink(D)))+g(p(1,mink(D)),q(2,bla,mink(D)))+g(p(1,mink(D)),q(2,bla,mink(D)))*g(p(2,mink(D)),q(3,mink(D)))+g(p(1,mink(D)),q(2,bla,mink(D)))*g(q(2,bla,mink(D)),q(3,mink(D)))");
    assert!(!result.contains("mink(D,1)"));
}

#[test]
fn vakint_rank1_input_simplifies_to_dots() {
    test_initialize();
    let dim = symbol!("D");
    let mink: Representation<_> = Minkowski {}.new_rep(dim);

    let input = p!(0, slot!(mink, 2))
        * p!(0, slot!(mink, 5))
        * p!(4, slot!(mink, 7))
        * q!(0, slot!(mink, 2))
        * q!(0, slot!(mink, 5))
        * q!(0, slot!(mink, 7));

    let result = input.to_dots();

    assert_snapshot!(
        result.to_bare_ordered_string(),
        @"(dot(p(0,mink(D)),q(0,mink(D))))^2*dot(p(4,mink(D)),q(0,mink(D)))"
    );
}

#[test]
fn schoonschip_settings_substitutes_metric_slots_in_plain_functions() {
    test_initialize();
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let cof: Representation<_> = ColorFundamental {}.new_rep(symbol!("N"));
    let coaf = cof.dual();
    let p = T.rank_one_tensor_symbol("P");
    let f = symbol!("F");
    let p_stripped = function!(p, 1, mink.to_symbolic([]));

    let self_dual = ETS.metric(&p_stripped, slot!(mink, mu).to_atom())
        * function!(f, symbol!("x"), slot!(mink, mu).to_atom(), symbol!("y"));
    assert_snapshot!(self_dual.schoonschip().to_bare_ordered_string(), @"F(x,P(1,mink(D)),y)");

    let dualizable = ETS.metric(slot!(cof, i).to_atom(), slot!(coaf, j).to_atom())
        * function!(f, symbol!("x"), slot!(cof, j).to_atom(), symbol!("y"));
    assert_snapshot!(dualizable.schoonschip().to_bare_ordered_string(), @"F(x,cof(N,i),y)");
}

#[test]
fn slot_detection_uses_representation_tags() {
    test_initialize();
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let cof: Representation<_> = ColorFundamental {}.new_rep(symbol!("N"));
    let coaf = cof.dual();
    let p = T.rank_one_tensor_symbol("P");

    assert!(is_slot(&slot!(mink, mu).to_atom()));
    assert!(is_slot(&slot!(coaf, j).to_atom()));
    assert!(is_slot(&function!(
        AIND_SYMBOLS.dind,
        slot!(cof, j).to_atom()
    )));
    assert!(!is_slot(&mink.to_symbolic([])));
    assert!(!is_slot(&function!(p, 1, mink.to_symbolic([]))));
}

#[test]
fn chain_like_metric_simplification_is_opt_in() {
    test_initialize();
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let bis: Representation<_> = Bispinor {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let f = symbol!("F");
    let i = slot!(bis, 1).to_atom();
    let j = slot!(bis, 2).to_atom();
    let nu = slot!(mink, 2).to_atom();
    let p_stripped = function!(p, 1, mink.to_symbolic([]));
    let expr =
        ETS.metric(&nu, &p_stripped) * chain!(&i, &j, function!(f, chain_in(), chain_out(), &nu));

    assert_snapshot!(expr.schoonschip().to_bare_ordered_string(), @"P(1,mink(D,2))*chain(bis(D,1),bis(D,2),F(in,out,mink(D,2)))");

    let result =
        expr.schoonschip_with_settings(&SchoonschipSettings::default().with_chain_like_functions());
    assert_snapshot!(result.to_bare_ordered_string(), @"chain(bis(D,1),bis(D,2),F(in,out,P(1,mink(D))))");
}

#[test]
fn chain_like_metric_simplification_keeps_compact_scalar_products() {
    test_initialize();
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let bis: Representation<_> = Bispinor {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let q = T.rank_one_tensor_symbol("Q");
    let f = symbol!("F");
    let i = slot!(bis, i).to_atom();
    let j = slot!(bis, j).to_atom();
    let p_stripped = function!(p, 1, mink.to_symbolic([]));
    let q_stripped = function!(q, 1, mink.to_symbolic([]));

    let expr = ETS.metric(&p_stripped, &q_stripped)
        * chain!(&i, &j, function!(f, chain_in(), chain_out(), &p_stripped));
    let result =
        expr.schoonschip_with_settings(&SchoonschipSettings::default().with_chain_like_functions());
    assert_snapshot!(result.to_bare_ordered_string(), @"chain(bis(D,i),bis(D,j),F(in,out,P(1,mink(D))))*g(P(1,mink(D)),Q(1,mink(D)))");
}

#[test]
fn chain_like_rank1_schoonschip_is_settings_controlled() {
    test_initialize();
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let bis: Representation<_> = Bispinor {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let f = symbol!("F");
    let i = slot!(bis, 1).to_atom();
    let j = slot!(bis, 2).to_atom();
    let mu = slot!(mink, 1).to_atom();
    let expr = function!(p, 1, &mu) * chain!(&i, &j, function!(f, chain_in(), chain_out(), &mu));

    let without_chain_like = expr.schoonschip_with_settings(&SchoonschipSettings::default());
    assert_snapshot!(without_chain_like.to_bare_ordered_string(), @"P(1,mink(D,1))*chain(bis(D,1),bis(D,2),F(in,out,mink(D,1)))");

    let without_rank1 = expr.schoonschip_with_settings(
        &SchoonschipSettings::default()
            .with_chain_like_functions()
            .without_rank1_tensors(),
    );
    assert_snapshot!(without_rank1.to_bare_ordered_string(), @"P(1,mink(D,1))*chain(bis(D,1),bis(D,2),F(in,out,mink(D,1)))");

    let result =
        expr.schoonschip_with_settings(&SchoonschipSettings::default().with_chain_like_functions());
    assert_snapshot!(result.to_bare_ordered_string(), @"chain(bis(D,1),bis(D,2),F(in,out,P(1,mink(D))))");
}

#[test]
fn chain_like_metric_simplification_handles_traces() {
    test_initialize();
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let bis: Representation<_> = Bispinor {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let f = symbol!("F");
    let mu = slot!(mink, 1).to_atom();
    let p_stripped = function!(p, 1, mink.to_symbolic([]));

    let expr = ETS.metric(&mu, &p_stripped)
        * trace!(
            bis.to_symbolic([]),
            function!(f, chain_in(), chain_out(), &mu)
        );
    let result =
        expr.schoonschip_with_settings(&SchoonschipSettings::default().with_chain_like_functions());
    assert_snapshot!(result.to_bare_ordered_string(), @"trace(bis(D),cyclic(F(in,out,P(1,mink(D)))))");
}

#[test]
fn chain_like_metric_simplification_handles_symmetric_traces() {
    test_initialize();
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let bis: Representation<_> = Bispinor {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let f = symbol!("F");
    let mu = slot!(mink, 1).to_atom();
    let p_stripped = function!(p, 1, mink.to_symbolic([]));

    let expr = ETS.metric(&mu, &p_stripped)
        * trace_sym!(
            bis.to_symbolic([]),
            function!(f, chain_in(), chain_out(), &mu)
        );
    let result =
        expr.schoonschip_with_settings(&SchoonschipSettings::default().with_chain_like_functions());
    assert_snapshot!(result.to_bare_ordered_string(), @"trace(bis(D),sym(F(in,out,P(1,mink(D)))))");
}

#[test]
fn chain_like_rank1_schoonschip_handles_dual_traces() {
    test_initialize();
    let cof: Representation<_> = ColorFundamental {}.new_rep(symbol!("N"));
    let coaf = cof.dual();
    let bis: Representation<_> = Bispinor {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let f = symbol!("F");
    let i = slot!(cof, 1).to_atom();
    let dual_i = slot!(coaf, 1).to_atom();
    let settings = SchoonschipSettings::default().with_chain_like_functions();

    let expr = function!(p, 1, &i)
        * trace!(
            bis.to_symbolic([]),
            function!(f, chain_in(), chain_out(), &dual_i)
        );
    let result = expr.schoonschip_with_settings(&settings);
    assert_snapshot!(result.to_bare_ordered_string(), @"trace(bis(D),cyclic(F(in,out,P(1,cof(N)))))");

    let expr = function!(p, 1, &i)
        * trace_sym!(
            bis.to_symbolic([]),
            function!(f, chain_in(), chain_out(), &dual_i)
        );
    let result = expr.schoonschip_with_settings(&settings);
    assert_snapshot!(result.to_bare_ordered_string(), @"trace(bis(D),sym(F(in,out,P(1,cof(N)))))");
}

#[test]
fn chain_like_metric_simplification_handles_chain_endpoints() {
    test_initialize();
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let f = symbol!("F");
    let mu = slot!(mink, 1).to_atom();
    let nu = slot!(mink, 2).to_atom();
    let p_stripped = function!(p, 1, mink.to_symbolic([]));

    let expr =
        ETS.metric(&mu, &p_stripped) * chain!(&mu, &nu, function!(f, chain_in(), chain_out()));
    let result =
        expr.schoonschip_with_settings(&SchoonschipSettings::default().with_chain_like_functions());
    assert_snapshot!(result.to_bare_ordered_string(), @"chain(P(1,mink(D)),mink(D,2),F(in,out))");
}

#[test]
fn contraction_orders_preserve_the_same_factorized_value() {
    use crate::tensor::SymbolicTensor;
    test_initialize();
    let input = (spenso::p!(mink!(4, 94101)) + spenso::q!(mink!(4, 94101)))
        * (spenso::vector!(order_r, mink!(4, 94101)) + spenso::vector!(order_s, mink!(4, 94101)));
    let source = SymbolicTensor::infer(input.clone()).unwrap();
    let expected = input.expand().schoonschip().expand();
    for order in [[0, 1], [1, 0]] {
        let result = source
            .contract(crate::tensor::ContractionSettings::default().with_order(&order))
            .unwrap()
            .resolved()
            .unwrap();
        assert_eq!(result.expression.expand(), expected);
    }
}

#[test]
fn mixed_registered_and_ordinary_momentum_dot() {
    test_initialize();
    let ordinary = symbol!("mixed_dot_ordinary");
    let mink = Minkowski {}.new_rep(4);
    let rep = mink.to_symbolic([]);
    let vector = function!(ordinary, &rep);
    let compact = p!(0, &vector);
    let expected = spenso::dot!(p!(0, &rep), &vector);
    assert_eq!(compact.to_dots(), expected);
    assert_eq!(expected.to_dots(), expected);
    let expanded = p!(0, slot!(mink, mu)) * function!(ordinary, slot!(mink, mu).to_atom());
    assert_eq!(expanded.to_dots(), expected);
    // A scalar function of a vector is not itself a rank-one tensor.
    let scalar = function!(symbol!("mixed_dot_scalar"), vector);
    assert_eq!(scalar.normalize_dots(), scalar);
}

#[test]
fn explicit_metric_vectors_are_not_scalar_labels() {
    test_initialize();
    let mink = Minkowski {}.new_rep(4);
    let rep = mink.to_symbolic([]);
    let p = function!(symbol!("metric_dot_p"), &rep);
    let q = function!(symbol!("metric_dot_q"), &rep);
    let expr = g!(slot!(mink, mu), &p) * g!(slot!(mink, mu), &q);
    assert_eq!(expr.to_dots(), spenso::dot!(p, q));
}

#[test]
fn network_scalar_leaves_preserve_cleanup_scope_and_numeric_expansion() {
    use crate::tensor::{SymbolicNetExt, SymbolicNetParse};
    use spenso::network::parsing::{ParseSettings, ShorthandParsing};

    test_initialize();
    let parse = |source: &str| {
        Atom::parse(
            source,
            "scalar_leaf_test",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap()
    };
    let _ = symbol!("scalar_leaf_test::broadcast", tag = T.broadcast);
    for (source, factored, expanded) in [
        ("s*(-2*(x+y)+z)", "s*(-2*(x+y)+z)", "s*(-2*x-2*y+z)"),
        ("s*(spenso::bracket(x)+y)", "s*(x+y)", "s*(x+y)"),
        ("s*(spenso::bracket(x,y)+z)", "s*(x*y+z)", "s*(x*y+z)"),
        ("s*(spenso::pure_scalar(x)+y)", "s*(x+y)", "s*(x+y)"),
        // A bracket inside an ordinary function is opaque to network parsing.
        (
            "s*(f(spenso::bracket(x))+y)",
            "s*(f(spenso::bracket(x))+y)",
            "s*(f(spenso::bracket(x))+y)",
        ),
        (
            "s*(broadcast(spenso::bracket(x))+y)",
            "s*(broadcast(x)+y)",
            "s*(broadcast(x)+y)",
        ),
    ] {
        let input = parse(source);
        let network = input
            .as_view()
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::expand_all(),
                precontract_scalars: false,
                ..ParseSettings::default()
            })
            .unwrap();
        let result = network.simple_execute::<()>().unwrap();
        assert_eq!(result, parse(factored), "{source}");
        assert_eq!(
            result.expand(),
            parse(expanded).expand(),
            "{source}, explicit expansion"
        );
    }

    let opaque = parse("s*(spenso::pure_scalar(x)+y)");
    assert_eq!(contracted_atom(opaque.as_view()).unwrap(), opaque);
    assert!(
        parse("spenso::pure_scalar(x,y)")
            .as_view()
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings {
                shorthand_parsing: ShorthandParsing::Opaque,
                ..ParseSettings::default()
            })
            .is_err()
    );
}

#[test]
fn network_scalar_sums_keep_internal_contractions_and_factored_spectators() {
    test_initialize();
    let spectator = (Atom::var(symbol!("scalar_spectator_x"))
        + Atom::var(symbol!("scalar_spectator_y")))
    .pow(8);
    let input = &spectator
        * (p!(0, mink!(4, 1)) + p!(1, mink!(4, 1)))
        * (p!(2, mink!(4, 1)) + p!(3, mink!(4, 1)));
    let actual = contracted_atom(input.as_view()).unwrap();
    let expected = spectator
        * (g!(p!(0, mink!(4)), p!(2, mink!(4)))
            + g!(p!(0, mink!(4)), p!(3, mink!(4)))
            + g!(p!(1, mink!(4)), p!(2, mink!(4)))
            + g!(p!(1, mink!(4)), p!(3, mink!(4))));
    assert_eq!(actual.expand(), expected.expand());
    assert_ne!(
        actual,
        actual.expand(),
        "scalar spectator must stay factored"
    );
}

#[test]
fn network_scalar_leaves_retain_parser_errors_and_opaque_compact_forms() {
    test_initialize();
    let parse = |source: &str| {
        Atom::parse(
            source,
            "scalar_leaf_test",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap()
    };
    let _ = symbol!("scalar_leaf_test::broadcast", tag = T.broadcast);
    T.rank_one_tensor_symbol("scalar_leaf_test::v");
    use crate::tensor::SymbolicNetParse;
    use spenso::network::parsing::ParseSettings;
    {
        for (source, diagnostic) in [
            ("s*(spenso::bracket()+z)", "empty bracket expression"),
            ("s*(spenso::dot(x)+z)", "Invalid dot function"),
            ("s*(spenso::chain()+z)", "Too many arguments"),
            ("s*(spenso::trace()+z)", "Too many arguments"),
            ("s*(spenso::pure_scalar(x,y)+z)", "Too many arguments"),
            ("s*(broadcast(x,y)+z)", "Too many arguments"),
            (
                "s*(spenso::aind(spenso::mink(4,a),x)+z)",
                "Not a slot, is composite",
            ),
        ] {
            let input = parse(source);
            let error = input
                .as_view()
                .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
                .unwrap_err()
                .to_string();
            assert!(error.contains(diagnostic), "{source}: {error}");
        }
        let malformed = parse("s*(v(0,spenso::mink(4,a,b))+z)");
        assert!(contracted_atom(malformed.as_view()).is_err());
        // The raw parser retains these opaque scalar forms. This is distinct
        // from the typed constructor's stricter builtin signature admission.
        for source in [
            "s*(spenso::g(x)+z)",
            "s*(v(0,spenso::mink())+z)",
            "s*(v(0,spenso::dind(spenso::mink(4),x))+z)",
            "s*(v(0,spenso::mink(f(D)))+z)",
        ] {
            let input = parse(source);
            use crate::tensor::SymbolicNetExt;
            let parsed = input
                .as_view()
                .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
                .unwrap();
            assert_eq!(parsed.simple_execute::<()>().unwrap(), input, "{source}");
        }
    }
}

#[test]
fn network_scalar_roots_preserve_compact_forms_and_parser_scope() {
    test_initialize();
    let parse = |source: &str| {
        Atom::parse(
            source,
            "scalar_root_test",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap()
    };
    T.rank_one_tensor_symbol("scalar_root_test::v");
    {
        for source in [
            "v(0,spenso::mink(4))",
            "v(0,spenso::mink(4))^2",
            "v(0,spenso::mink(4))*v(1,spenso::mink(4))",
            "spenso::mink(4)",
            "spenso::dind(spenso::lor(4))",
            "(x+y)^8*(u+w)",
        ] {
            let input = parse(source);
            assert_eq!(contracted_atom(input.as_view()).unwrap(), input, "{source}");
        }
        for (source, factored, expanded) in [
            ("-2*(x+y)", "-2*(x+y)", "-2*x-2*y"),
            ("x+spenso::bracket(y)", "x+y", "x+y"),
            (
                "x+f(spenso::bracket(y))",
                "x+f(spenso::bracket(y))",
                "x+f(spenso::bracket(y))",
            ),
        ] {
            let input = parse(source);
            let result = contracted_atom(input.as_view()).unwrap();
            assert_eq!(result, parse(factored), "{source}");
            assert_eq!(
                result.expand(),
                parse(expanded).expand(),
                "{source}, explicit expansion"
            );
        }
        let input = parse("x+spenso::bracket()");
        use crate::tensor::SymbolicNetParse;
        use spenso::network::parsing::ParseSettings;
        let error = input
            .as_view()
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap_err()
            .to_string();
        assert!(error.contains("empty bracket expression"), "{error}");
    }
}
