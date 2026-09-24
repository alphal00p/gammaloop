//! Contract a metric or tagged vector by replacing one compatible tensor slot.

use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    shadowing,
    structure::{
        representation::{LibraryRep, RepName},
        slot::{SlotMatch, SlotMatcher},
    },
};
use symbolica::atom::{
    Atom, AtomCore, AtomView, FunctionBuilder,
    representation::{FunView, MulView},
};

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

pub(crate) struct SlotContraction;

impl SlotContraction {
    /// Contract explicit metrics and optional vectors without expanding sums or interpreting powers.
    /// Dot normalization, including metric traces and powers, remains a
    /// prerequisite stage owned by the caller.
    pub(crate) fn run(view: AtomView<'_>, chain_like: bool, rank_one: bool) -> Atom {
        let mut current = view.to_owned();
        let mut slots = SlotMatcher::default();
        loop {
            if let AtomView::Mul(product) = current.as_view()
                && let Some(next) =
                    Self::contract_product(product, chain_like, rank_one, &mut slots)
            {
                if next == current {
                    return next;
                }
                current = next;
            }
            // The root product was already checked. A visitor can prune slot
            // payloads without causing replace_map to rebuild their ancestors.
            let mut root = true;
            let mut found = false;
            current.visitor(&mut |atom| {
                let is_root = std::mem::replace(&mut root, false);
                if found {
                    return false;
                }
                match atom {
                    AtomView::Mul(product) if !is_root => {
                        found =
                            Self::contract_product_once(product, chain_like, rank_one, &mut slots)
                                .is_some();
                    }
                    AtomView::Fun(_) if !matches!(slots.classify(atom), SlotMatch::Other) => {
                        return false;
                    }
                    _ => {}
                }
                !found
            });
            if !found {
                return current;
            }
            let next = current.replace_map(|atom, _, out| match atom {
                AtomView::Mul(product) => {
                    if let Some(contracted) =
                        Self::contract_product(product, chain_like, rank_one, &mut slots)
                    {
                        **out = contracted;
                    }
                }
                AtomView::Fun(_) if !matches!(slots.classify(atom), SlotMatch::Other) => {
                    out.set_from_view(&atom);
                }
                _ => {}
            });
            if next == current {
                return next;
            }
            current = next;
        }
    }

    // Finish this product before rebuilding enclosing sums and functions. Each
    // step still normalizes its factors and restarts metric-first selection.
    fn contract_product(
        product: MulView<'_>,
        chain_like: bool,
        rank_one: bool,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        let mut current = Self::contract_product_once(product, chain_like, rank_one, slots)?;
        while let AtomView::Mul(product) = current.as_view() {
            let Some(next) = Self::contract_product_once(product, chain_like, rank_one, slots)
            else {
                break;
            };
            if next == current {
                break;
            }
            current = next;
        }
        Some(current)
    }

    fn contract_product_once(
        product: MulView<'_>,
        chain_like: bool,
        rank_one: bool,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        // Preserve metric-first contraction. Canonical products already expose
        // their factors, so finding a source needs no commutative pattern search.
        for (position, factor) in product.iter().enumerate() {
            let AtomView::Fun(function) = factor else {
                continue;
            };
            if function.get_symbol() != ETS.metric {
                continue;
            }
            let mut arguments = function.iter();
            if arguments.len() != 2 {
                continue;
            }
            let first = arguments.next().unwrap();
            let second = arguments.next().unwrap();
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
            for (source, replacement) in [(first_slot, second), (second_slot, first)] {
                if let Some(source) = source
                    && let Some(replaced) = Self::replace_partner(
                        product,
                        position,
                        source,
                        &|| replacement.to_owned(),
                        chain_like,
                        slots,
                    )
                {
                    return Some(replaced);
                }
            }
        }
        if rank_one {
            for (position, factor) in product.iter().enumerate() {
                let AtomView::Fun(function) = factor else {
                    continue;
                };
                let Some((slot_position, source)) = Self::rank_one_endpoint(function, slots) else {
                    continue;
                };
                // Construct the compact vector only after finding a partner.
                if let Some(replaced) = Self::replace_partner(
                    product,
                    position,
                    source,
                    &|| {
                        let stripped = source.representation.base().to_symbolic([source.dimension]);
                        Self::replace_argument(function, slot_position, stripped)
                    },
                    chain_like,
                    slots,
                ) {
                    return Some(replaced);
                }
            }
        }
        None
    }

    /// Tagged vectors have exactly one slot, in their final argument. Earlier
    /// arguments are scalar parameters; do not search inside their payloads.
    fn rank_one_endpoint<'a>(
        function: FunView<'a>,
        slots: &mut SlotMatcher,
    ) -> Option<(usize, Endpoint<'a>)> {
        if !function.get_symbol().has_tag(&SPENSO_TAG.rank1) {
            return None;
        }
        let argument = slots.vector_argument(function)?;
        Endpoint::parse(argument, slots).map(|endpoint| (function.get_nargs() - 1, endpoint))
    }

    fn replace_partner(
        product: MulView<'_>,
        source_position: usize,
        source: Endpoint<'_>,
        replacement: &impl Fn() -> Atom,
        chain_like: bool,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        for (position, factor) in product.iter().enumerate() {
            if position == source_position {
                continue;
            }
            // Sums, powers and function payloads cannot consume an outside
            // index. Only direct tensor slots and opted-in chain bodies can.
            let AtomView::Fun(function) = factor else {
                continue;
            };
            if let Some(replaced) =
                Self::replace_function(function, source, replacement, chain_like, slots)
            {
                return Some(
                    product
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| *i != source_position)
                        .map(|(i, factor)| {
                            if i == position {
                                replaced.as_view()
                            } else {
                                factor
                            }
                        })
                        .product(),
                );
            }
        }
        None
    }

    fn replace_function(
        function: FunView<'_>,
        source: Endpoint<'_>,
        replacement: &impl Fn() -> Atom,
        chain_like: bool,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        let head = function.get_symbol();
        if head.is_scalar() || !matches!(slots.classify(function.as_view()), SlotMatch::Other) {
            return None;
        }
        if head.has_tag(&SPENSO_TAG.rank1) {
            let (position, endpoint) = Self::rank_one_endpoint(function, slots)?;
            return source
                .matches(&endpoint)
                .then(|| Self::replace_argument(function, position, replacement()));
        }
        for (position, argument) in function.iter().enumerate() {
            let replaced =
                if Endpoint::parse(argument, slots).is_some_and(|slot| source.matches(&slot)) {
                    Some(replacement())
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
        replacement: &impl Fn() -> Atom,
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
    use crate::shorthands::schoonschip::{Schoonschip, SchoonschipSettings};

    fn initialize_vectors() {
        crate::representations::initialize();
        SPENSO_TAG.rank_one_tensor_symbol("spenso::slot_contraction_p");
        SPENSO_TAG.rank_one_tensor_symbol("spenso::slot_contraction_q");
        symbolica::symbol!("spenso::slot_contraction_scalar"; Scalar);
    }

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
            let result = SlotContraction::run(expression.as_view(), false, false);
            assert_eq!(result, parse(expected), "{input}");
            assert_eq!(SlotContraction::run(result.as_view(), false, false), result);
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
                SlotContraction::run(expression.as_view(), false, false),
                parse(expected),
                "{input}"
            );
        }
    }

    #[test]
    fn product_contractions_finish_inside_independent_scopes() {
        initialize_vectors();
        for (input, expected) in [
            (
                "g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))*T(mink(4,c))+g(mink(4,d),mink(4,e))*g(mink(4,e),mink(4,f))*U(mink(4,f))",
                "T(mink(4,a))+U(mink(4,d))",
            ),
            (
                "f(g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))*T(mink(4,c)))",
                "f(T(mink(4,a)))",
            ),
            (
                "(x+g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))*T(mink(4,c)))^2",
                "(x+T(mink(4,a)))^2",
            ),
            (
                "g(mink(4,a),mink(4,b))*slot_contraction_p(mink(4,b))*T(mink(4,a))",
                "T(slot_contraction_p(mink(4)))",
            ),
        ] {
            let result = SlotContraction::run(parse(input).as_view(), false, true);
            assert_eq!(result, parse(expected), "{input}");
            assert_eq!(SlotContraction::run(result.as_view(), false, true), result);
        }
    }

    #[test]
    fn contraction_probe_prunes_root_slot_payloads() {
        crate::representations::initialize();
        for input in [
            "mink(g(mink(4,a),mink(4,b))*T(mink(4,b)),index)",
            "mink(4,g(mink(4,a),mink(4,b))*T(mink(4,b)))",
            "dind(cof(N,index),g(mink(4,a),mink(4,b))*T(mink(4,b)))",
        ] {
            let expression = parse(input);
            assert_eq!(
                SlotContraction::run(expression.as_view(), true, true),
                expression
            );
        }
    }

    #[test]
    fn exact_large_index_payloads_do_not_alias_narrowed_integers() {
        crate::representations::initialize();
        let distinct =
            parse("g(mink(4,1),mink(4,2))*T(mink(4,4294967297))*U(mink(4,k))*V(mink(4,k))");
        assert_eq!(
            SlotContraction::run(distinct.as_view(), false, false),
            distinct
        );
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
                SlotContraction::run(parse(input).as_view(), false, false),
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
                SlotContraction::run(parse(input).as_view(), false, false),
                parse(expected)
            );
        }
    }

    #[test]
    fn chain_endpoints_are_direct_slots_and_bodies_are_opt_in() {
        crate::representations::initialize();
        let endpoint = parse("g(bis(4,a),bis(4,b))*chain(bis(4,b),bis(4,c),F(in,out))");
        assert_eq!(
            SlotContraction::run(endpoint.as_view(), false, false),
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
                SlotContraction::run(expression.as_view(), false, false),
                expression,
                "{input}"
            );
            assert_eq!(
                SlotContraction::run(expression.as_view(), true, false),
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
            SlotContraction::run(expression.as_view(), false, false),
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
                SlotContraction::run(expression.as_view(), false, false),
                expression,
                "{input}"
            );
        }
    }

    #[test]
    fn vectors_and_compact_metrics_share_slot_substitution() {
        initialize_vectors();
        for (input, expected) in [
            (
                "slot_contraction_p(label,mink(4,a))*T(mink(4,a))",
                "T(slot_contraction_p(label,mink(4)))",
            ),
            (
                "g(mink(4,a),slot_contraction_p(label,mink(4)))*T(mink(4,a))",
                "T(slot_contraction_p(label,mink(4)))",
            ),
            (
                "slot_contraction_p(label,mink(4,a))*slot_contraction_q(other,mink(4,a))",
                "g(slot_contraction_p(label,mink(4)),slot_contraction_q(other,mink(4)))",
            ),
            (
                "slot_contraction_p(label,dind(cof(N,a)))*T(cof(N,a))",
                "T(slot_contraction_p(label,cof(N)))",
            ),
            (
                "slot_contraction_p(label,cof(N,a))*T(dind(cof(N,a)))",
                "T(slot_contraction_p(label,cof(N)))",
            ),
            (
                "slot_contraction_p(coad(Nc^2-1,f(a)))*T(coad(Nc^2-1,f(a)))",
                "T(slot_contraction_p(coad(Nc^2-1)))",
            ),
        ] {
            let result = parse(input).schoonschip();
            assert_eq!(result, parse(expected), "{input}");
            assert_eq!(result.schoonschip(), result, "{input}");
        }
    }

    #[test]
    fn vector_parameters_and_contraction_scopes_remain_opaque() {
        initialize_vectors();
        for input in [
            "slot_contraction_p(mink(4,4294967297))*T(mink(4,1))",
            "slot_contraction_p(mink(4,a))*T(mink(5,a))",
            "slot_contraction_p(cof(N,a))*T(cof(N,a))",
            "slot_contraction_p(mink(4,a))*slot_contraction_scalar(mink(4,a))",
            "slot_contraction_p(mink(4,a))*T(f(mink(4,a)))",
            "slot_contraction_p(mink(4,a))*T(mink(4,f(mink(4,a))))",
            "slot_contraction_p(mink(4,a))*T(mink(4,a))^2",
            "slot_contraction_p(mink(4,a))*(T(mink(4,a))+U(mink(4,a)))",
            "slot_contraction_p(mink(4,a),mink(4,b))*T(mink(4,b))",
            "slot_contraction_p(mink(4,a,extra))*T(mink(4,a))",
            "slot_contraction_p(dind(cof(N,a),extra))*T(cof(N,a))",
        ] {
            let expression = parse(input);
            assert_eq!(expression.schoonschip(), expression, "{input}");
        }
        for (input, expected) in [
            (
                "slot_contraction_p(slot_contraction_scalar(mink(4,a)),mink(4,a))*T(mink(4,a))",
                "T(slot_contraction_p(slot_contraction_scalar(mink(4,a)),mink(4)))",
            ),
            (
                "slot_contraction_p(mink(4,a))*T(slot_contraction_scalar(mink(4,a)),mink(4,a))",
                "T(slot_contraction_scalar(mink(4,a)),slot_contraction_p(mink(4)))",
            ),
            (
                "(1+slot_contraction_p(mink(4,a))*T(mink(4,a)))^2",
                "(1+T(slot_contraction_p(mink(4))))^2",
            ),
        ] {
            assert_eq!(parse(input).schoonschip(), parse(expected), "{input}");
        }
    }

    #[test]
    fn vector_contraction_obeys_rank_one_and_chain_settings() {
        initialize_vectors();
        let expression = parse("slot_contraction_p(mink(4,a))*T(mink(4,a))");
        assert_eq!(
            expression.schoonschip_with_settings(
                &SchoonschipSettings::default().without_rank1_tensors(),
            ),
            expression,
        );
        for (input, expected) in [
            (
                "slot_contraction_p(mink(4,a))*chain(bis(4,i),bis(4,j),F(in,out,mink(4,a)))",
                "chain(bis(4,i),bis(4,j),F(in,out,slot_contraction_p(mink(4))))",
            ),
            (
                "slot_contraction_p(mink(4,a))*trace(bis(4),cyclic(F(in,out,mink(4,a))))",
                "trace(bis(4),cyclic(F(in,out,slot_contraction_p(mink(4)))))",
            ),
            (
                "slot_contraction_p(mink(4,a))*trace(bis(4),sym(F(in,out,mink(4,a))))",
                "trace(bis(4),sym(F(in,out,slot_contraction_p(mink(4)))))",
            ),
        ] {
            let expression = parse(input);
            assert_eq!(expression.schoonschip(), expression, "{input}");
            assert_eq!(
                expression.schoonschip_with_settings(
                    &SchoonschipSettings::default().with_chain_like_functions(),
                ),
                parse(expected),
                "{input}",
            );
            assert_eq!(
                expression.schoonschip_with_settings(
                    &SchoonschipSettings::default()
                        .with_chain_like_functions()
                        .without_rank1_tensors(),
                ),
                expression,
                "{input}",
            );
        }
        let input = parse("g(mink(4,a),mink(4,b))*slot_contraction_p(mink(4,b))");
        assert_eq!(
            input
                .schoonschip_with_settings(&SchoonschipSettings::default().without_rank1_tensors()),
            parse("slot_contraction_p(mink(4,a))"),
        );
    }
}
