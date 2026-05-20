use insta::assert_snapshot;
use spenso::{
    chain, g, mink,
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    p, q,
    shadowing::symbolica_utils::AtomCoreExt,
    slot,
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

use crate::{
    representations::{Bispinor, ColorFundamental},
    test_support::test_initialize,
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

    let result =
        (&p1 * &q2).schoonschip_with_net::<false, AbstractIndex>(&SchoonschipSettings::full());
    assert_snapshot!(result.to_bare_ordered_string(),@"g(p(1,mink(D)),q(2,bla,mink(D)))");

    let result = (&p1 * &p1).schoonschip();
    assert_snapshot!(result.to_bare_ordered_string(), @"g(p(1,mink(D)),p(1,mink(D)))");

    let result = p!(1, &p2_stripped).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"g(p(1,mink(dim)),p(2,mink(dim)))");

    let result = g!(slot!(mink, 1), &p1_stripped).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"p(1,mink(dim),mink(D,1))");

    let result = p1.clone().pow(Atom::num(4)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"(g(p(1,mink(D)),p(1,mink(D))))^2");

    let result = p1.clone().pow(Atom::num(3)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"g(p(1,mink(D)),p(1,mink(D)))*p(1,mink(D,1))");

    let metric = g!(slot!(mink, 1), slot!(mink, 2));
    let result = metric.clone().pow(Atom::num(4)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"D^2");

    let result = metric.pow(Atom::num(3)).normalize_dots();
    assert_snapshot!(result.to_bare_ordered_string(), @"D*g(mink(D,1),mink(D,2))");

    let result = (&p1 * (&q2 + &p2 * &q3_2 * &q2_2))
        .schoonschip_with_net::<false, AbstractIndex>(&SchoonschipSettings::partial());
    assert_snapshot!(result.to_bare_ordered_string(),@"g(p(1,mink(D)),p(2,mink(D)))*g(q(2,bla,mink(D)),q(3,mink(D)))+g(p(1,mink(D)),q(2,bla,mink(D)))");

    let result = (&p1 * (&q2 + &p2 * (&q3_2 * &q2_2 + &p2_2 * &q2_2)))
        .schoonschip_with_net::<false, AbstractIndex>(&SchoonschipSettings::full());
    assert_snapshot!(result.to_bare_ordered_string(),@"(g(p(2,mink(D)),q(2,bla,mink(D)))+g(q(2,bla,mink(D)),q(3,mink(D))))*g(p(1,mink(D)),p(2,mink(D)))+g(p(1,mink(D)),q(2,bla,mink(D)))");

    let expr = (p1 + q3 * p1_2 * q2_2) * (q2 + p2);

    let result = expr.schoonschip_with_net::<false, AbstractIndex>(&SchoonschipSettings::full());
    assert_snapshot!(result.to_bare_ordered_string(),@"(g(p(1,mink(D)),q(2,bla,mink(D)))*q(3,mink(D,1))+p(1,mink(D,1)))*(p(2,mink(D,1))+q(2,bla,mink(D,1)))");

    let result = expr.schoonschip_with_net::<false, AbstractIndex>(
        &SchoonschipSettings::full().with_expanded_contracted_sums(),
    );
    let result = result.to_bare_ordered_string();
    assert_snapshot!(result, @"g(p(1,mink(D)),p(2,mink(D)))+g(p(1,mink(D)),q(2,bla,mink(D)))+g(p(1,mink(D)),q(2,bla,mink(D)))*g(p(2,mink(D)),q(3,mink(D)))+g(p(1,mink(D)),q(2,bla,mink(D)))*g(q(2,bla,mink(D)),q(3,mink(D)))");
    assert!(!result.contains("mink(D,1)"));

    let result = expr.schoonschip_with_net::<false, AbstractIndex>(
        &SchoonschipSettings::partial().with_expanded_contracted_sums(),
    );
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

    let result = expr.schoonschip_with_settings(
        &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
    );
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
    let result = expr.schoonschip_with_settings(
        &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
    );
    assert_snapshot!(result.to_bare_ordered_string(), @"chain(bis(D,i),bis(D,j),F(in,out,P(1,mink(D))))*g(P(1,mink(D)),Q(1,mink(D)))");
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
    let result = expr.schoonschip_with_settings(
        &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
    );
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
    let result = expr.schoonschip_with_settings(
        &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
    );
    assert_snapshot!(result.to_bare_ordered_string(), @"trace(bis(D),sym(F(in,out,P(1,mink(D)))))");
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
    let result = expr.schoonschip_with_settings(
        &SchoonschipSettings::single_pass(None).with_chain_like_functions(),
    );
    assert_snapshot!(result.to_bare_ordered_string(), @"chain(P(1,mink(D)),mink(D,2),F(in,out))");
}

#[test]
fn benchmark_modes_output() {
    test_initialize();
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
        .schoonschip_with_net::<false, AbstractIndex>(&SchoonschipSettings::single_pass(Some(1)))
        .to_bare_ordered_string();
    let depth_first_depth_one = expr
        .clone()
        .schoonschip_with_net::<false, AbstractIndex>(&SchoonschipSettings::partial())
        .to_bare_ordered_string();
    let breadth_first_depth_one = expr
        .clone()
        .schoonschip_with_net::<false, AbstractIndex>(&SchoonschipSettings::breadth_first(Some(1)))
        .to_bare_ordered_string();
    let full_top = expr
        .schoonschip_with_net::<false, AbstractIndex>(&SchoonschipSettings::full())
        .to_bare_ordered_string();

    assert_ne!(single_pass_depth_one, full_top);
    assert_eq!(depth_first_depth_one, full_top);
    assert_eq!(breadth_first_depth_one, full_top);
}
