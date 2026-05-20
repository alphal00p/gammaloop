use std::sync::LazyLock;

use spenso::{
    chain,
    network::{library::symbolic::ETS, tags::SPENSO_TAG as T},
    shadowing,
    tensors::parametric::atomcore::PatternReplacement,
    trace, trace_sym,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView},
    function,
    id::Replacement,
};

use crate::{W_, shorthands::metric::not_slot};

use super::{api::Schoonschip, settings::SchoonschipSettings};

static METRIC_FUNCTION_CONTRACTIONS: LazyLock<[Replacement; 3]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);
    let function_with_replacement = function!(W_.a_, W_.a___, W_.c_, W_.b___);

    [
        // g(i,j)*T(...,j,...)->T(...,i,...)
        Replacement::new(
            (function!(ETS.metric, &self_dual, W_.c_)
                * function!(W_.a_, W_.a___, &self_dual, W_.b___))
            .to_pattern(),
            function_with_replacement.clone(),
        ),
        // g(i,j)*T(...,d(j),...)->T(...,i,...)
        Replacement::new(
            (function!(ETS.metric, W_.c_, &dualizable)
                * function!(W_.a_, W_.a___, &dualizable_dual, W_.b___))
            .to_pattern(),
            function_with_replacement.clone(),
        ),
        // g(i,d(j))*T(...,j,...)->T(...,i,...)
        Replacement::new(
            (function!(ETS.metric, W_.c_, &dualizable_dual)
                * function!(W_.a_, W_.a___, &dualizable, W_.b___))
            .to_pattern(),
            function_with_replacement.clone(),
        ),
    ]
});

static METRIC_FUNCTION_CONTRACTIONS_ON_CHAIN: LazyLock<[Replacement; 3]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

    fn chain_<'a>(a: impl Into<AtomOrView<'a>>) -> Atom {
        chain!(W_.x_, W_.y_, W_.x___, a.into(), W_.y___)
    }
    let function_with_replacement = chain_(function!(W_.a_, W_.a___, W_.c_, W_.b___));

    [
        // g(i,j)*chain(x,y,...,T(...,j,...),...)->chain(x,y,...,T(...,i,...),...)
        Replacement::new(
            (function!(ETS.metric, &self_dual, W_.c_)
                * chain_(function!(W_.a_, W_.a___, &self_dual, W_.b___)))
            .to_pattern(),
            function_with_replacement.clone(),
        ),
        // g(i,j)*chain(x,y,...,T(...,d(j),...),...)->chain(x,y,...,T(...,i,...),...)
        Replacement::new(
            (function!(ETS.metric, W_.c_, &dualizable)
                * chain_(function!(W_.a_, W_.a___, &dualizable_dual, W_.b___)))
            .to_pattern(),
            function_with_replacement.clone(),
        ),
        // g(i,d(j))*chain(x,y,...,T(...,j,...),...)->chain(x,y,...,T(...,i,...),...)
        Replacement::new(
            (function!(ETS.metric, W_.c_, &dualizable_dual)
                * chain_(function!(W_.a_, W_.a___, &dualizable, W_.b___)))
            .to_pattern(),
            function_with_replacement.clone(),
        ),
    ]
});
static METRIC_FUNCTION_CONTRACTIONS_ON_TRACE: LazyLock<[Replacement; 6]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

    fn trace_cyclic_<'a>(a: impl Into<AtomOrView<'a>>) -> Atom {
        trace!(W_.y_, a.into(), W_.x___)
    }
    fn trace_sym_<'a>(a: impl Into<AtomOrView<'a>>) -> Atom {
        trace_sym!(W_.y_, a.into(), W_.x___)
    }
    let cyclic_function_with_replacement = trace_cyclic_(function!(W_.a_, W_.a___, W_.c_, W_.b___));
    let sym_function_with_replacement = trace_sym_(function!(W_.a_, W_.a___, W_.c_, W_.b___));

    [
        // g(i,j)*trace(y,...,T(...,j,...),...)->trace(y,...,T(...,i,...),...)
        Replacement::new(
            (function!(ETS.metric, &self_dual, W_.c_)
                * trace_cyclic_(function!(W_.a_, W_.a___, &self_dual, W_.b___)))
            .to_pattern(),
            cyclic_function_with_replacement.clone(),
        ),
        // g(i,j)*trace(y,...,T(...,d(j),...),...)->trace(y,...,T(...,i,...),...)
        Replacement::new(
            (function!(ETS.metric, W_.c_, &dualizable)
                * trace_cyclic_(function!(W_.a_, W_.a___, &dualizable_dual, W_.b___)))
            .to_pattern(),
            cyclic_function_with_replacement.clone(),
        ),
        // g(i,d(j))*trace(y,...,T(...,j,...),...)->trace(y,...,T(...,i,...),...)
        Replacement::new(
            (function!(ETS.metric, W_.c_, &dualizable_dual)
                * trace_cyclic_(function!(W_.a_, W_.a___, &dualizable, W_.b___)))
            .to_pattern(),
            cyclic_function_with_replacement.clone(),
        ),
        // g(i,j)*trace(y,sym(...,T(...,j,...),...))->trace(y,sym(...,T(...,i,...),...))
        Replacement::new(
            (function!(ETS.metric, &self_dual, W_.c_)
                * trace_sym_(function!(W_.a_, W_.a___, &self_dual, W_.b___)))
            .to_pattern(),
            sym_function_with_replacement.clone(),
        ),
        // g(i,j)*trace(y,sym(...,T(...,d(j),...),...))->trace(y,sym(...,T(...,i,...),...))
        Replacement::new(
            (function!(ETS.metric, W_.c_, &dualizable)
                * trace_sym_(function!(W_.a_, W_.a___, &dualizable_dual, W_.b___)))
            .to_pattern(),
            sym_function_with_replacement.clone(),
        ),
        // g(i,d(j))*trace(y,sym(...,T(...,j,...),...))->trace(y,sym(...,T(...,i,...),...))
        Replacement::new(
            (function!(ETS.metric, W_.c_, &dualizable_dual)
                * trace_sym_(function!(W_.a_, W_.a___, &dualizable, W_.b___)))
            .to_pattern(),
            sym_function_with_replacement,
        ),
    ]
});
/// f(...,i)*g(...,i)->g(f(...,rep)*g(...,rep)) but only when ... is not a slot.
/// This is more costly than just using the rank1 tag on f and g, so should be an optional setting.
static _DOT_PRODUCT: LazyLock<[Replacement; 1]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);

    [Replacement::new(
        (function!(W_.f_, W_.a___, &self_dual) * function!(W_.g_, W_.b___, &self_dual))
            .to_pattern(),
        ETS.metric(
            function!(W_.f_, W_.a___, &self_dual_stripped),
            function!(W_.g_, W_.b___, &self_dual_stripped),
        ),
    )
    .with_conditions(not_slot(W_.a___) & not_slot(W_.b___))]
});

static VECTOR_DOT_PRODUCTS: LazyLock<[Replacement; 2]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

    [
        // f(...,i)*g(...,i)->g(f(...,rep)*g(...,rep)) for tagged rank1 f,g
        Replacement::new(
            (T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual])
                * T.rank1_::<1, _>([&Atom::var(W_.a___), &self_dual]))
            .to_pattern(),
            ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                T.rank1_::<1, _>([&Atom::var(W_.a___), &self_dual_stripped]),
            ),
        ),
        // f(...,i)*g(...,d(i))->g(f(...,rep)*g(...,rep)) for tagged rank1 f,g
        Replacement::new(
            (T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable])
                * T.rank1_::<1, _>([&Atom::var(W_.a___), &dualizable_dual]))
            .to_pattern(),
            ETS.metric(
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
                T.rank1_::<1, _>([&Atom::var(W_.a___), &dualizable_stripped]),
            ),
        ),
    ]
});

static SCHOONSCHIP_VECTOR: LazyLock<[Replacement; 3]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);

    [
        // a(...,i,...)*p(...,i)->a(....,p(...,rep),...) for p a rank1 tagged function
        Replacement::new(
            (function!(W_.a_, W_.a___, &self_dual, W_.b___)
                * T.rank1_::<0, _>([Atom::var(W_.c___), self_dual.clone()]))
            .to_pattern(),
            function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                W_.b___
            ),
        ),
        // a(...,i,...)*p(...,d(i))->a(....,p(...,rep),...) for p a rank1 tagged function
        Replacement::new(
            (function!(W_.a_, W_.a___, &dualizable, W_.b___)
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_dual]))
            .to_pattern(),
            function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
                W_.b___
            ),
        ),
        // a(...,d(i),...)*p(...,i)->a(....,p(...,rep),...) for p a rank1 tagged function
        Replacement::new(
            (function!(W_.a_, W_.a___, &dualizable_dual, W_.b___)
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable]))
            .to_pattern(),
            function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([Atom::var(W_.c___), dualizable_stripped]),
                W_.b___
            ),
        ),
    ]
});

static SCHOONSCHIP_VECTOR_ON_CHAIN: LazyLock<[Replacement; 3]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);
    fn chain_<'a>(a: impl Into<AtomOrView<'a>>) -> Atom {
        chain!(W_.x_, W_.y_, W_.x___, a.into(), W_.y___)
    }

    [
        // chain(x,y,...,a(...,i,...),...)*p(...,i)->chain(x,y,...,a(....,p(...,rep),...),...) for p a rank1 tagged function
        Replacement::new(
            (chain_(function!(W_.a_, W_.a___, &self_dual, W_.b___))
                * T.rank1_::<0, _>([Atom::var(W_.c___), self_dual.clone()]))
            .to_pattern(),
            chain_(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                W_.b___
            )),
        ),
        // chain(x,y,...,a(...,i,...),...)*p(...,d(i))->chain(x,y,...,a(....,p(...,rep),...),...) for p a rank1 tagged function
        Replacement::new(
            (chain_(function!(W_.a_, W_.a___, &dualizable, W_.b___))
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_dual]))
            .to_pattern(),
            chain_(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
                W_.b___
            )),
        ),
        // chain(x,y,...,a(...,d(i),...),...)*p(...,i)->chain(x,y,...,a(....,p(...,rep),...),...) for p a rank1 tagged function
        Replacement::new(
            (chain_(function!(W_.a_, W_.a___, &dualizable_dual, W_.b___))
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable]))
            .to_pattern(),
            chain_(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([Atom::var(W_.c___), dualizable_stripped]),
                W_.b___
            )),
        ),
    ]
});

static SCHOONSCHIP_VECTOR_ON_TRACE: LazyLock<[Replacement; 6]> = LazyLock::new(|| {
    let self_dual = T.self_dual_::<0, _>([W_.d_, W_.i_]);
    let self_dual_stripped = T.self_dual_::<0, _>([W_.d_]);
    let dualizable = T.dualizable_::<0, _>([W_.d_, W_.i_]);
    let dualizable_stripped = T.dualizable_::<0, _>([W_.d_]);
    let dualizable_dual = T.dualizable_dual_::<0, _>([W_.d_, W_.i_]);
    fn trace_cyclic_<'a>(a: impl Into<AtomOrView<'a>>) -> Atom {
        trace!(W_.x_, a.into(), W_.z___)
    }
    fn trace_sym_<'a>(a: impl Into<AtomOrView<'a>>) -> Atom {
        trace_sym!(W_.x_, a.into(), W_.z___)
    }
    [
        // trace(rep1,...,a(...,i,...),...)*p(...,i)->trace(rep1,...,a(....,p(...,rep),...),...) for p a rank1 tagged function
        Replacement::new(
            (trace_cyclic_(function!(W_.a_, W_.a___, &self_dual, W_.b___))
                * T.rank1_::<0, _>([Atom::var(W_.c___), self_dual.clone()]))
            .to_pattern(),
            trace_cyclic_(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                W_.b___
            )),
        ),
        // trace(rep1,...,a(...,i,...),...)*p(...,d(i))->trace(rep1,...,a(....,p(...,rep),...),...) for p a rank1 tagged function
        Replacement::new(
            (trace_cyclic_(function!(W_.a_, W_.a___, &dualizable, W_.b___))
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_dual]))
            .to_pattern(),
            trace_cyclic_(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
                W_.b___
            )),
        ),
        // trace(rep1,...,a(...,d(i),...),...)*p(...,i)->trace(rep1,...,a(....,p(...,rep),...),...) for p a rank1 tagged function
        Replacement::new(
            (trace_cyclic_(function!(W_.a_, W_.a___, &dualizable, W_.b___))
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable]))
            .to_pattern(),
            trace_cyclic_(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([Atom::var(W_.c___), dualizable_stripped.clone()]),
                W_.b___
            )),
        ),
        // trace(rep1,sym(...,a(...,i,...),...))*p(...,i)->trace(rep1,sym(...,a(....,p(...,rep),...),...))
        Replacement::new(
            (trace_sym_(function!(W_.a_, W_.a___, &self_dual, W_.b___))
                * T.rank1_::<0, _>([Atom::var(W_.c___), self_dual]))
            .to_pattern(),
            trace_sym_(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &self_dual_stripped]),
                W_.b___
            )),
        ),
        // trace(rep1,sym(...,a(...,i,...),...))*p(...,d(i))->trace(rep1,sym(...,a(....,p(...,rep),...),...))
        Replacement::new(
            (trace_sym_(function!(W_.a_, W_.a___, &dualizable, W_.b___))
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_dual]))
            .to_pattern(),
            trace_sym_(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable_stripped]),
                W_.b___
            )),
        ),
        // trace(rep1,sym(...,a(...,d(i),...),...))*p(...,i)->trace(rep1,sym(...,a(....,p(...,rep),...),...))
        Replacement::new(
            (trace_sym_(function!(W_.a_, W_.a___, &dualizable, W_.b___))
                * T.rank1_::<0, _>([&Atom::var(W_.c___), &dualizable]))
            .to_pattern(),
            trace_sym_(function!(
                W_.a_,
                W_.a___,
                T.rank1_::<0, _>([Atom::var(W_.c___), dualizable_stripped]),
                W_.b___
            )),
        ),
    ]
});

struct SchoonschipWithSettings<'a> {
    settings: &'a SchoonschipSettings,
}

impl SchoonschipWithSettings<'_> {
    fn run(&self, view: AtomView<'_>) -> Atom {
        let mut current = view.to_owned();
        loop {
            let next = self.apply_once(current.as_view());
            if next == current {
                return next;
            }
            current = next;
        }
    }

    fn apply_once(&self, view: AtomView<'_>) -> Atom {
        let metric_simplified = view
            .normalize_dots()
            .to_owned()
            .replace_multiple_repeat(&*METRIC_FUNCTION_CONTRACTIONS);

        let simplified = if self.settings.schoonschip_rank1_tensors {
            metric_simplified
                .replace_multiple(&*VECTOR_DOT_PRODUCTS)
                .replace_multiple_repeat(&*SCHOONSCHIP_VECTOR)
        } else {
            metric_simplified
        };

        let simplified = if self.settings.simplify_chain_like_functions {
            self.apply_chain_like_rules(simplified)
        } else {
            simplified
        };

        simplified.normalize_dots()
    }

    fn apply_chain_like_rules(&self, expression: Atom) -> Atom {
        let (metric_trace_rules, metric_trace_sym_rules) =
            METRIC_FUNCTION_CONTRACTIONS_ON_TRACE.split_at(3);
        let (vector_trace_rules, vector_trace_sym_rules) = SCHOONSCHIP_VECTOR_ON_TRACE.split_at(3);

        let simplified = expression
            .replace_multiple_repeat(&*METRIC_FUNCTION_CONTRACTIONS_ON_CHAIN)
            .replace_multiple_repeat(metric_trace_rules)
            .replace_multiple_repeat(&*SCHOONSCHIP_VECTOR_ON_CHAIN)
            .replace_multiple_repeat(vector_trace_rules);

        if Self::contains_symmetric_projector(simplified.as_view()) {
            simplified
                .replace_multiple_repeat(metric_trace_sym_rules)
                .replace_multiple_repeat(vector_trace_sym_rules)
        } else {
            simplified
        }
    }

    fn contains_symmetric_projector(expr: AtomView<'_>) -> bool {
        match expr {
            AtomView::Fun(f) => {
                f.get_symbol() == *shadowing::SYM
                    || f.iter().any(Self::contains_symmetric_projector)
            }
            AtomView::Add(add) => add.iter().any(Self::contains_symmetric_projector),
            AtomView::Mul(mul) => mul.iter().any(Self::contains_symmetric_projector),
            AtomView::Pow(pow) => pow.iter().any(Self::contains_symmetric_projector),
            AtomView::Num(_) | AtomView::Var(_) => false,
        }
    }
}

pub(super) fn schoonschip_with_settings(
    view: AtomView<'_>,
    settings: &SchoonschipSettings,
) -> Atom {
    SchoonschipWithSettings { settings }.run(view)
}
