use std::sync::LazyLock;

use spenso::{
    chain,
    network::{library::symbolic::ETS, parsing::AtomStructureExt, tags::SPENSO_TAG as T},
    shadowing, trace, trace_sym,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomOrView, AtomView},
    function,
    id::Replacement,
};
use symbolica_utils::PatternReplacement;

use crate::{
    W_,
    shorthands::{bracket::BracketNormalizer, metric::not_slot},
};

use super::{
    api::Schoonschip, metric_contraction::MetricContraction, settings::SchoonschipSettings,
};

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
    .when(not_slot(W_.a___) & not_slot(W_.b___))]
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
            (trace_cyclic_(function!(W_.a_, W_.a___, &dualizable_dual, W_.b___))
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
            (trace_sym_(function!(W_.a_, W_.a___, &dualizable_dual, W_.b___))
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

pub(crate) struct SchoonschipWithSettings<'a> {
    pub(crate) settings: &'a SchoonschipSettings,
}

impl SchoonschipWithSettings<'_> {
    pub(crate) fn run(&self, view: AtomView<'_>) -> Atom {
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
        let normalized = BracketNormalizer::normalize(view).normalize_dots();
        let repeated_indices = normalized.has_repeated_explicit_indices();
        let metric_simplified = if repeated_indices {
            MetricContraction::run(
                normalized.as_view(),
                self.settings.simplify_chain_like_functions,
            )
        } else {
            normalized
        };

        let simplified = if self.settings.schoonschip_rank1_tensors {
            metric_simplified
                .replace_multiple(&*VECTOR_DOT_PRODUCTS)
                .replace_multiple_repeat(&*SCHOONSCHIP_VECTOR)
        } else {
            metric_simplified
        };

        // These rules also need two occurrences of an explicit slot. Reuse the
        // input scan: preceding contractions cannot create one when none existed.
        let simplified = if repeated_indices
            && self.settings.simplify_chain_like_functions
            && self.settings.schoonschip_rank1_tensors
        {
            self.apply_chain_like_vector_rules(simplified)
        } else {
            simplified
        };

        BracketNormalizer::normalize(simplified.normalize_dots().as_view())
    }

    fn apply_chain_like_vector_rules(&self, expression: Atom) -> Atom {
        let (trace_rules, symmetric_trace_rules) = SCHOONSCHIP_VECTOR_ON_TRACE.split_at(3);
        let simplified = expression
            .replace_multiple_repeat(&*SCHOONSCHIP_VECTOR_ON_CHAIN)
            .replace_multiple_repeat(trace_rules);

        if Self::contains_symmetric_projector(simplified.as_view()) {
            simplified.replace_multiple_repeat(symmetric_trace_rules)
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
