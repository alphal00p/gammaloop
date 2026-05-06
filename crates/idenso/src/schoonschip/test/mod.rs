use insta::assert_snapshot;
use spenso::{
    chain,
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    shadowing::symbolica_utils::AtomCoreExt,
    slot,
    structure::{
        abstract_index::{AIND_SYMBOLS, AbstractIndex},
        representation::{Minkowski, RepName, Representation},
        slot::IsAbstractSlot,
    },
    trace,
};
use symbolica::{
    atom::{Atom, AtomCore},
    function, symbol,
};

use crate::representations::{Bispinor, ColorFundamental};

use super::{Schoonschip, SchoonschipSettings, chain_like::is_slot};

fn chain_in() -> Atom {
    Atom::var(T.chain_in)
}

fn chain_out() -> Atom {
    Atom::var(T.chain_out)
}

#[test]
fn simple_dot() {
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let q = T.rank_one_tensor_symbol("Q");

    let p1 = function!(p, 1, slot!(mink, 1).to_atom());
    let p2 = function!(p, 2, slot!(mink, 1).to_atom());
    let p1_2 = function!(p, 1, slot!(mink, 2).to_atom());
    let p2_2 = function!(p, 2, slot!(mink, 2).to_atom());
    let p1_stripped = function!(p, 1, mink.to_symbolic([]));
    let p2_stripped = function!(p, 2, mink.to_symbolic([]));

    let q2 = function!(q, 2, symbol!("bla"), slot!(mink, 1).to_atom());
    let q2_2 = function!(q, 2, symbol!("bla"), slot!(mink, 2).to_atom());

    let q3 = function!(q, 3, slot!(mink, 1).to_atom());
    let q3_2 = function!(q, 3, slot!(mink, 2).to_atom());

    let result = (&p1 * &q2)
        .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::full());
    assert_snapshot!(result.to_bare_ordered_string(),@"g(P(1,mink(D)),Q(2,bla,mink(D)))");

    let result = (&p1 * &p1).schoonschip();
    assert_snapshot!(result.to_bare_ordered_string(), @"g(P(1,mink(D)),P(1,mink(D)))");

    let result = function!(p, 1, &p2_stripped).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"g(P(1,mink(D)),P(2,mink(D)))");

    let result = ETS
        .metric(slot!(mink, 1).to_atom(), &p1_stripped)
        .normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"P(1,mink(D,1))");

    let result = p1.clone().pow(Atom::num(4)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"(g(P(1,mink(D)),P(1,mink(D))))^2");

    let result = p1.clone().pow(Atom::num(3)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"P(1,mink(D,1))*g(P(1,mink(D)),P(1,mink(D)))");

    let metric = ETS.metric(slot!(mink, 1).to_atom(), slot!(mink, 2).to_atom());
    let result = metric.clone().pow(Atom::num(4)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"D^2");

    let result = metric.pow(Atom::num(3)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"D*g(mink(D,1),mink(D,2))");

    let result = (&p1 * (&q2 + &p2 * &q3_2 * &q2_2))
        .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::partial());
    assert_snapshot!(result.to_bare_ordered_string(),@"g(P(1,mink(D)),P(2,mink(D)))*g(Q(2,bla,mink(D)),Q(3,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))");

    let result = (&p1 * (&q2 + &p2 * (&q3_2 * &q2_2 + &p2_2 * &q2_2)))
        .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::full());
    assert_snapshot!(result.to_bare_ordered_string(),@"(g(P(2,mink(D)),Q(2,bla,mink(D)))+g(Q(2,bla,mink(D)),Q(3,mink(D))))*g(P(1,mink(D)),P(2,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))");

    let expr = (p1 + q3 * p1_2 * q2_2) * (q2 + p2);

    let result =
        expr.schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::full());
    assert_snapshot!(result.to_bare_ordered_string(),@"(P(1,mink(D,1))+Q(3,mink(D,1))*g(P(1,mink(D)),Q(2,bla,mink(D))))*(P(2,mink(D,1))+Q(2,bla,mink(D,1)))");

    let result = expr.schoonschip_with_net::<false, true, AbstractIndex>(
        &SchoonschipSettings::full().with_expanded_contracted_sums(),
    );
    let result = result.to_bare_ordered_string();
    assert_snapshot!(result, @"g(P(1,mink(D)),P(2,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(P(2,mink(D)),Q(3,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(Q(2,bla,mink(D)),Q(3,mink(D)))");
    assert!(!result.contains("mink(D,1)"));

    let result = expr.schoonschip_with_net::<false, true, AbstractIndex>(
        &SchoonschipSettings::partial().with_expanded_contracted_sums(),
    );
    let result = result.to_bare_ordered_string();
    assert_snapshot!(result, @"g(P(1,mink(D)),P(2,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(P(2,mink(D)),Q(3,mink(D)))+g(P(1,mink(D)),Q(2,bla,mink(D)))*g(Q(2,bla,mink(D)),Q(3,mink(D)))");
    assert!(!result.contains("mink(D,1)"));
}

#[test]
fn schoonschip_settings_substitutes_metric_slots_in_plain_functions() {
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

    let result = expr.schoonschip_with_settings(
        &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
    );
    assert_snapshot!(result.to_bare_ordered_string(), @"chain(bis(D,1),bis(D,2),F(in,out,P(1,mink(D))))");
}

#[test]
fn chain_like_metric_simplification_handles_generic_nested_arguments() {
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let bis: Representation<_> = Bispinor {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let f = symbol!("F");
    let h = symbol!("H");
    let mu = slot!(mink, 1).to_atom();
    let i = slot!(bis, i).to_atom();
    let j = slot!(bis, j).to_atom();
    let p_stripped = function!(p, 1, mink.to_symbolic([]));

    let expr = ETS.metric(&mu, &p_stripped)
        * chain!(
            &i,
            &j,
            function!(f, chain_in(), chain_out(), symbol!("x"), function!(h, &mu))
        );
    let result = expr.schoonschip_with_settings(
        &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
    );
    assert_snapshot!(result.to_bare_ordered_string(), @"chain(bis(D,i),bis(D,j),F(in,out,x,H(P(1,mink(D)))))");
}

#[test]
fn chain_like_metric_simplification_keeps_compact_scalar_products() {
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
    let result = expr.schoonschip_with_settings(
        &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
    );
    assert_snapshot!(result.to_bare_ordered_string(), @"chain(bis(D,i),bis(D,j),F(in,out,P(1,mink(D))))*g(P(1,mink(D)),Q(1,mink(D)))");
}

#[test]
fn chain_like_metric_simplification_handles_traces() {
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
    let result = expr.schoonschip_with_settings(
        &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
    );
    assert_snapshot!(result.to_bare_ordered_string(), @"trace(bis(D),F(in,out,P(1,mink(D))))");
}

#[test]
fn chain_like_metric_simplification_handles_chain_endpoints() {
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let f = symbol!("F");
    let mu = slot!(mink, 1).to_atom();
    let nu = slot!(mink, 2).to_atom();
    let p_stripped = function!(p, 1, mink.to_symbolic([]));

    let expr =
        ETS.metric(&mu, &p_stripped) * chain!(&mu, &nu, function!(f, chain_in(), chain_out()));
    let result = expr.schoonschip_with_settings(
        &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
    );
    assert_snapshot!(result.to_bare_ordered_string(), @"chain(P(1,mink(D)),mink(D,2),F(in,out))");
}

#[test]
fn benchmark_modes_output() {
    let mink: Representation<_> = Minkowski {}.new_rep(symbol!("D"));
    let p = T.rank_one_tensor_symbol("P");
    let q = T.rank_one_tensor_symbol("Q");

    let p1 = function!(p, 1, slot!(mink, 1).to_atom());
    let p2 = function!(p, 2, slot!(mink, 1).to_atom());
    let p2_2 = function!(p, 2, slot!(mink, 2).to_atom());

    let q2 = function!(q, 2, symbol!("bla"), slot!(mink, 1).to_atom());
    let q2_2 = function!(q, 2, symbol!("bla"), slot!(mink, 2).to_atom());
    let q3_2 = function!(q, 3, slot!(mink, 2).to_atom());

    let expr = &p1 * (&q2 + &p2 * (&q3_2 * &q2_2 + &p2_2 * &q2_2));

    let single_pass_depth_one = expr
        .clone()
        .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::single_pass(
            Some(1),
        ))
        .to_bare_ordered_string();
    let depth_first_depth_one = expr
        .clone()
        .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::partial())
        .to_bare_ordered_string();
    let breadth_first_depth_one = expr
        .clone()
        .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::breadth_first(
            Some(1),
        ))
        .to_bare_ordered_string();
    let full_top = expr
        .schoonschip_with_net::<false, true, AbstractIndex>(&SchoonschipSettings::full())
        .to_bare_ordered_string();

    assert_ne!(single_pass_depth_one, full_top);
    assert_eq!(depth_first_depth_one, full_top);
    assert_eq!(breadth_first_depth_one, full_top);
}
