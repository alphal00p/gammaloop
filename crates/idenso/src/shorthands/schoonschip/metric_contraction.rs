//! Metric-first contraction: select one metric, then replace one compatible slot.

use std::sync::LazyLock;

use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    shadowing,
    structure::{
        representation::{LibraryRep, RepName},
        slot::{SlotMatch, SlotMatcher},
    },
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, representation::FunView},
    function,
    id::{MatchSettings, Pattern},
};

use crate::W_;

#[derive(Clone, Copy)]
struct Endpoint<'a> {
    representation: LibraryRep,
    dimension: AtomView<'a>,
    index: AtomView<'a>,
}

impl<'a> Endpoint<'a> {
    fn parse(value: AtomView<'a>, slots: &mut SlotMatcher) -> Option<Self> {
        let SlotMatch::Explicit(view) = slots.classify(value) else {
            return None;
        };
        Some(Self {
            representation: slots.representation(view).ok()?,
            dimension: view.dimension(),
            index: view.index(),
        })
    }

    fn matches(&self, other: &Self) -> bool {
        self.index == other.index && self.compatible(other)
    }

    fn compatible(&self, other: &Self) -> bool {
        self.dimension == other.dimension && self.representation.matches(&other.representation)
    }
}

static METRIC_AND_REMAINDER: LazyLock<Pattern> =
    LazyLock::new(|| (function!(ETS.metric, W_.a_, W_.b_) * Atom::var(W_.r___)).to_pattern());

pub(crate) struct MetricContraction;

impl MetricContraction {
    /// Contract explicit metrics without expanding sums or interpreting powers.
    /// Dot normalization, including metric traces and powers, remains a
    /// prerequisite stage owned by the caller.
    pub(crate) fn run(view: AtomView<'_>, chain_like: bool) -> Atom {
        let mut current = view.to_owned();
        let mut slots = SlotMatcher::default();
        loop {
            let next = current.replace_map(|atom, _, out| match atom {
                AtomView::Mul(_) => {
                    if let Some(contracted) = Self::contract_product(atom, chain_like, &mut slots) {
                        **out = contracted;
                    }
                }
                AtomView::Fun(_) if !matches!(slots.classify(atom), SlotMatch::Other) => {
                    out.set_from_view(&atom);
                }
                // Local contractions inside a function or power are valid;
                // replace_partner prevents contractions across those scopes.
                _ => {}
            });
            if next == current {
                return next;
            }
            current = next;
        }
    }

    fn contract_product(
        view: AtomView<'_>,
        chain_like: bool,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        let settings = MatchSettings::new()
            .min_level(0)
            .max_level(0)
            .level_is_tree_depth(true)
            .partial(false);
        let mut candidates = view.pattern_match(&METRIC_AND_REMAINDER, None, Some(&settings));
        while let Some(candidate) = candidates.next_detailed() {
            let captures = candidate.match_stack;
            let first = captures.get_atom(W_.a_)?;
            let second = captures.get_atom(W_.b_)?;
            // The symmetric head yields both orders; each candidate below
            // already tries both contraction directions.
            if first > second {
                continue;
            }
            let first_slot = Endpoint::parse(first, slots);
            let second_slot = Endpoint::parse(second, slots);
            if let (Some(first), Some(second)) = (first_slot, second_slot)
                && (!first.compatible(&second) || first.matches(&second))
            {
                continue;
            }
            if first_slot.is_none() && second_slot.is_none() {
                continue;
            }
            let remainder = captures.get(W_.r___)?.to_atom();
            for (source, replacement) in [(first_slot, second), (second_slot, first)] {
                if let Some(source) = source
                    && let Some(replaced) = Self::replace_partner(
                        remainder.as_view(),
                        source,
                        replacement,
                        chain_like,
                        slots,
                    )
                {
                    return Some(replaced);
                }
            }
        }
        None
    }

    fn replace_partner(
        remainder: AtomView<'_>,
        source: Endpoint<'_>,
        replacement: AtomView<'_>,
        chain_like: bool,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        let mut changed = false;
        let result = remainder.replace_map(|atom, _, out| match atom {
            AtomView::Fun(function) => {
                if !changed
                    && !function.get_symbol().is_scalar()
                    && let Some(replaced) =
                        Self::replace_function(function, source, replacement, chain_like, slots)
                {
                    changed = true;
                    **out = replaced;
                } else {
                    // Setting the original value prunes this subtree. Returning
                    // without setting would descend into dimensions and metadata.
                    out.set_from_view(&atom);
                }
            }
            AtomView::Add(_) | AtomView::Pow(_) => out.set_from_view(&atom),
            AtomView::Mul(_) | AtomView::Num(_) | AtomView::Var(_) => {}
        });
        changed.then_some(result)
    }

    fn replace_function(
        function: FunView<'_>,
        source: Endpoint<'_>,
        replacement: AtomView<'_>,
        chain_like: bool,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        let head = function.get_symbol();
        if head.is_scalar() || !matches!(slots.classify(function.as_view()), SlotMatch::Other) {
            return None;
        }
        for (position, argument) in function.iter().enumerate() {
            let replaced =
                if Endpoint::parse(argument, slots).is_some_and(|slot| source.matches(&slot)) {
                    Some(replacement.to_owned())
                } else if chain_like && head == SPENSO_TAG.chain && position >= 2 {
                    if let AtomView::Fun(factor) = argument {
                        Self::replace_function(factor, source, replacement, false, slots)
                    } else {
                        None
                    }
                } else if chain_like && head == SPENSO_TAG.trace && position == 1 {
                    if let AtomView::Fun(projector) = argument
                        && [*shadowing::CYCLIC, *shadowing::SYM].contains(&projector.get_symbol())
                    {
                        Self::replace_projector(projector, source, replacement, slots)
                    } else {
                        None
                    }
                } else {
                    None
                };
            if let Some(replaced) = replaced {
                return Some(Self::replace_argument(function, position, replaced));
            }
        }
        None
    }

    fn replace_projector(
        projector: FunView<'_>,
        source: Endpoint<'_>,
        replacement: AtomView<'_>,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        for (position, factor) in projector.iter().enumerate() {
            if let AtomView::Fun(function) = factor
                && let Some(replaced) =
                    Self::replace_function(function, source, replacement, false, slots)
            {
                return Some(Self::replace_argument(projector, position, replaced));
            }
        }
        None
    }

    fn replace_argument(function: FunView<'_>, position: usize, replacement: Atom) -> Atom {
        let mut result = FunctionBuilder::new(function.get_symbol());
        for (i, argument) in function.iter().enumerate() {
            result = result.add_arg(if i == position {
                replacement.as_view()
            } else {
                argument
            });
        }
        result.finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(source: &str) -> Atom {
        Atom::parse(
            source,
            "spenso",
            symbolica::parser::ParseSettings::symbolica(),
        )
        .unwrap()
    }

    #[test]
    fn metric_first_substitution_is_oriented_and_slot_aware() {
        crate::representations::initialize();
        for (input, expected) in [
            ("g(mink(4,a),mink(4,b))*T(mink(4,b))", "T(mink(4,a))"),
            (
                "g(mink(4,a),mink(4,b))*T(mink(4,a))*U(mink(4,b))",
                "T(mink(4,b))*U(mink(4,b))",
            ),
            ("g(lor(4,a),dind(lor(4,b)))*T(lor(4,b))", "T(lor(4,a))"),
            ("g(mink(4,a),P(mink(4)))*T(mink(4,a))", "T(P(mink(4)))"),
        ] {
            let expression = parse(input);
            let result = MetricContraction::run(expression.as_view(), false);
            assert_eq!(result, parse(expected), "{input}");
            assert_eq!(MetricContraction::run(result.as_view(), false), result);
        }
    }

    #[test]
    fn independent_contractions_inside_powers_and_functions_are_preserved() {
        crate::representations::initialize();
        for (input, expected) in [
            (
                "(1+g(mink(4,a),mink(4,b))*T(mink(4,a))*U(mink(4,b)))^2",
                "(1+T(mink(4,b))*U(mink(4,b)))^2",
            ),
            (
                "f(g(mink(4,a),mink(4,b))*T(mink(4,a))*U(mink(4,b)))",
                "f(T(mink(4,b))*U(mink(4,b)))",
            ),
        ] {
            let expression = parse(input);
            assert_eq!(
                MetricContraction::run(expression.as_view(), false),
                parse(expected),
                "{input}"
            );
        }
    }

    #[test]
    fn exact_large_index_payloads_do_not_alias_narrowed_integers() {
        crate::representations::initialize();
        let distinct =
            parse("g(mink(4,1),mink(4,2))*T(mink(4,4294967297))*U(mink(4,k))*V(mink(4,k))");
        assert_eq!(MetricContraction::run(distinct.as_view(), false), distinct);
        for (input, expected) in [
            (
                "g(mink(4,1),mink(4,4294967297))*T(mink(4,4294967297))",
                "T(mink(4,1))",
            ),
            (
                "g(mink(4,4294967297),mink(4,2))*T(mink(4,4294967297))",
                "T(mink(4,2))",
            ),
        ] {
            assert_eq!(
                MetricContraction::run(parse(input).as_view(), false),
                parse(expected)
            );
        }
    }

    #[test]
    fn symbolic_dimension_and_index_payloads_remain_exact() {
        crate::representations::initialize();
        for (input, expected) in [
            (
                "g(coad(Nc^2-1,a),coad(Nc^2-1,b))*T(coad(Nc^2-1,b))",
                "T(coad(Nc^2-1,a))",
            ),
            ("g(mink(4,f(a)),mink(4,b))*T(mink(4,f(a)))", "T(mink(4,b))"),
            (
                "g(mink(4,f(a)),mink(4,b))*T(mink(4,f(c)))",
                "g(mink(4,f(a)),mink(4,b))*T(mink(4,f(c)))",
            ),
        ] {
            assert_eq!(
                MetricContraction::run(parse(input).as_view(), false),
                parse(expected)
            );
        }
    }

    #[test]
    fn chain_endpoints_are_direct_slots_and_bodies_are_opt_in() {
        crate::representations::initialize();
        let endpoint = parse("g(bis(4,a),bis(4,b))*chain(bis(4,b),bis(4,c),F(in,out))");
        assert_eq!(
            MetricContraction::run(endpoint.as_view(), false),
            parse("chain(bis(4,a),bis(4,c),F(in,out))")
        );
        for (input, expected) in [
            (
                "g(mink(4,a),mink(4,b))*chain(bis(4,i),bis(4,j),F(in,out,mink(4,b)))",
                "chain(bis(4,i),bis(4,j),F(in,out,mink(4,a)))",
            ),
            (
                "g(mink(4,a),mink(4,b))*trace(bis(4),cyclic(F(in,out,mink(4,b))))",
                "trace(bis(4),cyclic(F(in,out,mink(4,a))))",
            ),
            (
                "g(mink(4,a),mink(4,b))*trace(bis(4),sym(F(in,out,mink(4,b))))",
                "trace(bis(4),sym(F(in,out,mink(4,a))))",
            ),
        ] {
            let expression = parse(input);
            assert_eq!(
                MetricContraction::run(expression.as_view(), false),
                expression,
                "{input}"
            );
            assert_eq!(
                MetricContraction::run(expression.as_view(), true),
                parse(expected),
                "{input}"
            );
        }
    }

    #[test]
    fn metric_candidates_continue_after_a_noncontracting_metric() {
        crate::representations::initialize();
        let expression = parse("g(mink(4,a),mink(4,b))*g(mink(4,c),mink(4,d))*T(mink(4,d))");
        assert_eq!(
            MetricContraction::run(expression.as_view(), false),
            parse("g(mink(4,a),mink(4,b))*T(mink(4,c))")
        );
    }

    #[test]
    fn metric_first_substitution_preserves_opaque_scopes() {
        crate::representations::initialize();
        for input in [
            "g(mink(4,a),mink(4,b))*T(mink(5,b))",
            "g(mink(4,a),lor(4,b))*T(lor(4,b))",
            "g(lor(4,a),lor(4,b))*T(lor(4,b))",
            "g(lor(4,a),lor(4,b))*T(dind(lor(4,b)))",
            "g(mink(4,a),mink(4,a))*T(mink(4,a))",
            "g(mink(4,a),mink(4,b))*f(T(mink(4,b)))",
            "g(mink(4,a),mink(4,b))*T(mink(4,f(mink(4,b))))",
            "g(mink(4,a),mink(4,b))*T(mink(4,b))^2",
            "g(mink(4,a),mink(4,b))*(T(mink(4,b))+U(mink(4,b)))",
            "g(P(mink(4)),Q(mink(4)))*T(P(mink(4)))",
        ] {
            let expression = parse(input);
            assert_eq!(
                MetricContraction::run(expression.as_view(), false),
                expression,
                "{input}"
            );
        }
    }
}
