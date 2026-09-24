//! Contract a metric or tagged vector by replacing one compatible tensor slot.

use ahash::AHashMap;

use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    shadowing,
    structure::{
        representation::{LibraryRep, RepName},
        slot::{SlotMatch, SlotMatcher},
    },
};
use symbolica::atom::{
    Atom, AtomCore, AtomOrView, AtomView, FunctionBuilder,
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

    fn same_index(&self, other: &Self) -> bool {
        self.index == other.index
            && self.dimension == other.dimension
            && self.representation.base() == other.representation.base()
    }
}

enum SlotReplacement {
    Absent,
    Replaced(Atom),
    Ambiguous,
}

pub(crate) struct SlotContraction;

impl SlotContraction {
    /// Contract explicit metrics and optional vectors while preserving scalar factors.
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
                // A tensor normalizer can expose another factor, for example
                // the sign of a reordered epsilon. Recheck the normalized
                // product before pruning the root from the nested walk.
                continue;
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
                        found = Self::find_contraction(
                            product.iter(),
                            chain_like,
                            rank_one,
                            &mut slots,
                        )
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

    // Normalize each changed tensor immediately, but rebuild the surrounding
    // product only after its contractions finish. Unchanged factors remain
    // borrowed, including potentially large scalar spectators.
    fn contract_product(
        product: MulView<'_>,
        chain_like: bool,
        rank_one: bool,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        if let Some(contracted) = Self::contract_metric_components(product, slots) {
            return Some(contracted);
        }
        let mut contraction = Self::find_contraction(product.iter(), chain_like, rank_one, slots)?;
        let mut factors: Vec<_> = product.iter().map(AtomOrView::View).collect();
        loop {
            let (source, partner, replacement) = contraction;
            if replacement.is_zero() {
                return Some(Atom::Zero);
            }
            factors[partner] = AtomOrView::Atom(replacement);
            factors.remove(source);
            let Some(next) = Self::find_contraction(
                factors.iter().map(AtomCore::as_atom_view),
                chain_like,
                rank_one,
                slots,
            ) else {
                return Some(Atom::mul_many(factors));
            };
            contraction = next;
        }
    }

    /// Internal indices in a metric path only transport its two boundary slots;
    /// a closed path contributes its dimension. Resolve those connections before
    /// constructing any new metric, without parsing unrelated tensor payloads.
    fn contract_metric_components(product: MulView<'_>, slots: &mut SlotMatcher) -> Option<Atom> {
        if product.iter().len() < 4 {
            return None;
        }
        let mut metrics = Vec::new();
        for (position, factor) in product.iter().enumerate() {
            let AtomView::Fun(function) = factor else {
                continue;
            };
            if function.get_symbol() != ETS.metric || function.get_nargs() != 2 {
                continue;
            }
            let mut arguments = function.iter();
            let arguments = [arguments.next().unwrap(), arguments.next().unwrap()];
            let [Some(first), Some(second)] = arguments.map(|a| Endpoint::parse(a, slots)) else {
                continue;
            };
            if first.compatible(&second) && !first.matches(&second) {
                metrics.push((position, arguments, [first, second]));
            }
        }
        if metrics.len() < 3 {
            return None;
        }
        let mut occurrences: AHashMap<_, (usize, usize, bool)> =
            AHashMap::with_capacity(2 * metrics.len());
        let mut neighbors = vec![[None::<usize>; 2]; metrics.len()];
        for (edge, (_, _, endpoints)) in metrics.iter().enumerate() {
            for (side, endpoint) in endpoints.iter().enumerate() {
                let key = (
                    endpoint.representation.base(),
                    endpoint.dimension,
                    endpoint.index,
                );
                if let Some((other_edge, other_side, repeated)) = occurrences.get_mut(&key) {
                    // More than two metric occurrences have no unambiguous
                    // path interpretation; retain ordered slot substitution.
                    if std::mem::replace(repeated, true) {
                        return None;
                    }
                    if endpoint.matches(&metrics[*other_edge].2[*other_side]) {
                        neighbors[edge][side] = Some(*other_edge);
                        neighbors[*other_edge][*other_side] = Some(edge);
                    }
                } else {
                    occurrences.insert(key, (edge, side, false));
                }
            }
        }
        // A tensor occurrence in addition to two metric endpoints makes the
        // incidence ambiguous. Do not turn that into a path contraction; the
        // ordered substitution below continues to own such expressions.
        let mut metric_positions = metrics.iter().map(|metric| metric.0).peekable();
        for (position, factor) in product.iter().enumerate() {
            if metric_positions.peek() == Some(&position) {
                metric_positions.next();
                continue;
            }
            let mut ambiguous = false;
            factor.visitor(&mut |atom| {
                if ambiguous {
                    return false;
                }
                if let Some(endpoint) = Endpoint::parse(atom, slots) {
                    let key = (
                        endpoint.representation.base(),
                        endpoint.dimension,
                        endpoint.index,
                    );
                    ambiguous = occurrences.get(&key).is_some_and(|entry| entry.2);
                    return false;
                }
                !matches!(atom, AtomView::Fun(f) if f.get_symbol().is_scalar())
                    && matches!(slots.classify(atom), SlotMatch::Other)
            });
            if ambiguous {
                return None;
            }
        }
        let mut visited = vec![false; metrics.len()];
        let mut removed = vec![false; product.iter().len()];
        let mut replacements = Vec::new();
        let mut pending = Vec::new();
        let mut component = Vec::new();
        let mut boundary = Vec::with_capacity(2);
        for start in 0..metrics.len() {
            if visited[start] {
                continue;
            }
            pending.push(start);
            component.clear();
            boundary.clear();
            while let Some(edge) = pending.pop() {
                if std::mem::replace(&mut visited[edge], true) {
                    continue;
                }
                component.push(edge);
                for (side, neighbor) in neighbors[edge].iter().enumerate() {
                    if let Some(neighbor) = neighbor {
                        pending.push(*neighbor);
                    } else {
                        boundary.push(metrics[edge].1[side]);
                    }
                }
            }
            if component.len() < 2 {
                continue;
            }
            let contracted = match boundary.as_slice() {
                [] => metrics[start].2[0].dimension.to_owned(),
                [first, second] => FunctionBuilder::new(ETS.metric)
                    .add_arg(*first)
                    .add_arg(*second)
                    .finish(),
                _ => unreachable!("a metric component is a path or a cycle"),
            };
            for &edge in &component {
                removed[metrics[edge].0] = true;
            }
            replacements.push(contracted);
        }
        if replacements.is_empty() {
            return None;
        }
        Some(Atom::mul_many(
            product
                .iter()
                .enumerate()
                .filter_map(|(i, factor)| (!removed[i]).then_some(AtomOrView::View(factor)))
                .chain(replacements.into_iter().map(AtomOrView::Atom)),
        ))
    }

    fn find_contraction<'a>(
        factors: impl Iterator<Item = AtomView<'a>> + Clone,
        chain_like: bool,
        rank_one: bool,
        slots: &mut SlotMatcher,
    ) -> Option<(usize, usize, Atom)> {
        // Preserve metric-first contraction. Canonical products already expose
        // their factors, so finding a source needs no commutative pattern search.
        for (position, factor) in factors.clone().enumerate() {
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
                        factors.clone(),
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
            for (position, factor) in factors.clone().enumerate() {
                let AtomView::Fun(function) = factor else {
                    continue;
                };
                let Some((slot_position, source)) = Self::rank_one_endpoint(function, slots) else {
                    continue;
                };
                // Construct the compact vector only after finding a partner.
                if let Some(replaced) = Self::replace_partner(
                    factors.clone(),
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
        // Substitution can leave a metric inside a factored sum while its
        // tensor partner sits outside, for example (g(a,b)*g(c,d) -
        // g(a,d)*g(c,b))*epsilon(a,e,f,h). Pass only that partner through the
        // sum; unrelated products retain their factorization. A normalized
        // metric with a compact vector becomes a rank-one source instead.
        for (source, factor) in factors.clone().enumerate() {
            if !matches!(factor, AtomView::Add(_)) {
                continue;
            }
            for (partner, tensor) in factors.clone().enumerate() {
                if let AtomView::Fun(function) = tensor
                    && function.get_symbol() != ETS.metric
                    && let Some(replacement) =
                        Self::contract_linear_source(factor, tensor, chain_like, rank_one, slots)
                {
                    return Some((partner, source, replacement));
                }
            }
        }
        None
    }

    fn contract_linear_source(
        expression: AtomView<'_>,
        partner: AtomView<'_>,
        chain_like: bool,
        rank_one: bool,
        slots: &mut SlotMatcher,
    ) -> Option<Atom> {
        match expression {
            AtomView::Fun(function)
                if function.get_symbol() == ETS.metric
                    || (rank_one && function.get_symbol().has_tag(&SPENSO_TAG.rank1)) =>
            {
                Self::find_contraction(
                    [expression, partner].into_iter(),
                    chain_like,
                    rank_one,
                    slots,
                )
                .map(|(_, _, result)| result)
            }
            AtomView::Add(sum) => {
                // Every branch must absorb the partner before it can be
                // removed from the surrounding product.
                let terms = sum
                    .iter()
                    .map(|term| {
                        Self::contract_linear_source(term, partner, chain_like, rank_one, slots)
                    })
                    .collect::<Option<Vec<_>>>()?;
                Some(Atom::add_many(terms))
            }
            AtomView::Mul(product) => {
                for (position, factor) in product.iter().enumerate() {
                    if let Some(replaced) =
                        Self::contract_linear_source(factor, partner, chain_like, rank_one, slots)
                    {
                        return Some(Atom::mul_many(product.iter().enumerate().map(
                            |(i, factor)| {
                                if i == position {
                                    replaced.as_view()
                                } else {
                                    factor
                                }
                            },
                        )));
                    }
                }
                None
            }
            _ => None,
        }
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

    fn replace_partner<'a>(
        factors: impl Iterator<Item = AtomView<'a>>,
        source_position: usize,
        source: Endpoint<'_>,
        replacement: &impl Fn() -> Atom,
        chain_like: bool,
        slots: &mut SlotMatcher,
    ) -> Option<(usize, usize, Atom)> {
        for (position, factor) in factors.enumerate() {
            if position == source_position {
                continue;
            }
            let result = match factor {
                AtomView::Fun(function) => {
                    Self::replace_function(function, source, replacement, chain_like, false, slots)
                }
                AtomView::Add(_) => {
                    Self::replace_linear(factor, source, replacement, chain_like, slots)
                }
                _ => SlotReplacement::Absent,
            };
            if let SlotReplacement::Replaced(replaced) = result {
                return Some((source_position, position, replaced));
            }
        }
        None
    }

    /// A sum exposes a slot only when every term contains it exactly once.
    /// Products pass it through one factor; all other factors remain borrowed.
    /// Powers and function metadata retain their independent scopes. Construct
    /// replacements inside this linear syntax without distributing products.
    fn replace_linear(
        expression: AtomView<'_>,
        source: Endpoint<'_>,
        replacement: &impl Fn() -> Atom,
        chain_like: bool,
        slots: &mut SlotMatcher,
    ) -> SlotReplacement {
        match expression {
            AtomView::Fun(function) => {
                Self::replace_function(function, source, replacement, chain_like, true, slots)
            }
            AtomView::Add(sum) => {
                let mut replaced = Vec::new();
                let mut absent = false;
                for term in sum.iter() {
                    match Self::replace_linear(term, source, replacement, chain_like, slots) {
                        SlotReplacement::Replaced(term) if !absent => replaced.push(term),
                        SlotReplacement::Absent if replaced.is_empty() => absent = true,
                        _ => return SlotReplacement::Ambiguous,
                    }
                }
                if absent {
                    SlotReplacement::Absent
                } else {
                    SlotReplacement::Replaced(Atom::add_many(replaced))
                }
            }
            AtomView::Mul(product) => {
                let mut changed = None;
                for (position, factor) in product.iter().enumerate() {
                    match Self::replace_linear(factor, source, replacement, chain_like, slots) {
                        SlotReplacement::Replaced(factor) if changed.is_none() => {
                            changed = Some((position, factor));
                        }
                        SlotReplacement::Absent => {}
                        _ => return SlotReplacement::Ambiguous,
                    }
                }
                let Some((position, replacement)) = changed else {
                    return SlotReplacement::Absent;
                };
                SlotReplacement::Replaced(Atom::mul_many(product.iter().enumerate().map(
                    |(i, factor)| {
                        if i == position {
                            replacement.as_view()
                        } else {
                            factor
                        }
                    },
                )))
            }
            AtomView::Pow(power) => {
                // A tensor power can expose a nonlinear occurrence. It cannot
                // be ignored when another factor offers the same index. Only
                // a source-free power is an opaque scalar spectator here.
                match Self::replace_linear(
                    power.get_base_exp().0,
                    source,
                    replacement,
                    chain_like,
                    slots,
                ) {
                    SlotReplacement::Absent => SlotReplacement::Absent,
                    _ => SlotReplacement::Ambiguous,
                }
            }
            _ => SlotReplacement::Absent,
        }
    }

    fn replace_function(
        function: FunView<'_>,
        source: Endpoint<'_>,
        replacement: &impl Fn() -> Atom,
        chain_like: bool,
        linear: bool,
        slots: &mut SlotMatcher,
    ) -> SlotReplacement {
        let head = function.get_symbol();
        if head.is_scalar() || !matches!(slots.classify(function.as_view()), SlotMatch::Other) {
            return SlotReplacement::Absent;
        }
        if head.has_tag(&SPENSO_TAG.rank1) {
            let Some((position, endpoint)) = Self::rank_one_endpoint(function, slots) else {
                return SlotReplacement::Absent;
            };
            return if source.matches(&endpoint) {
                SlotReplacement::Replaced(Self::replace_argument(function, position, replacement()))
            } else if linear && source.same_index(&endpoint) {
                SlotReplacement::Ambiguous
            } else {
                SlotReplacement::Absent
            };
        }
        let mut changed = None;
        for (position, argument) in function.iter().enumerate() {
            let replaced = if let Some(endpoint) = Endpoint::parse(argument, slots) {
                if source.matches(&endpoint) {
                    SlotReplacement::Replaced(replacement())
                } else if linear && source.same_index(&endpoint) {
                    // A same-variance occurrence cannot hide an internal pair
                    // while another occurrence consumes the external metric.
                    SlotReplacement::Ambiguous
                } else {
                    SlotReplacement::Absent
                }
            } else if chain_like && head == SPENSO_TAG.chain && position >= 2 {
                if let AtomView::Fun(factor) = argument {
                    Self::replace_function(factor, source, replacement, false, linear, slots)
                } else {
                    SlotReplacement::Absent
                }
            } else if chain_like && head == SPENSO_TAG.trace && position == 1 {
                if let AtomView::Fun(projector) = argument
                    && [*shadowing::CYCLIC, *shadowing::SYM].contains(&projector.get_symbol())
                {
                    Self::replace_projector(projector, source, replacement, linear, slots)
                } else {
                    SlotReplacement::Absent
                }
            } else {
                SlotReplacement::Absent
            };
            match replaced {
                SlotReplacement::Replaced(replaced) if !linear => {
                    return SlotReplacement::Replaced(Self::replace_argument(
                        function, position, replaced,
                    ));
                }
                SlotReplacement::Replaced(replaced) if changed.is_none() => {
                    changed = Some((position, replaced));
                }
                SlotReplacement::Absent => {}
                _ => return SlotReplacement::Ambiguous,
            }
        }
        match changed {
            Some((position, replacement)) => {
                SlotReplacement::Replaced(Self::replace_argument(function, position, replacement))
            }
            None => SlotReplacement::Absent,
        }
    }

    fn replace_projector(
        projector: FunView<'_>,
        source: Endpoint<'_>,
        replacement: &impl Fn() -> Atom,
        linear: bool,
        slots: &mut SlotMatcher,
    ) -> SlotReplacement {
        let mut changed = None;
        for (position, factor) in projector.iter().enumerate() {
            let AtomView::Fun(function) = factor else {
                continue;
            };
            match Self::replace_function(function, source, replacement, false, linear, slots) {
                SlotReplacement::Replaced(replaced) if !linear => {
                    return SlotReplacement::Replaced(Self::replace_argument(
                        projector, position, replaced,
                    ));
                }
                SlotReplacement::Replaced(replaced) if changed.is_none() => {
                    changed = Some((position, replaced));
                }
                SlotReplacement::Absent => {}
                _ => return SlotReplacement::Ambiguous,
            }
        }
        match changed {
            Some((position, replacement)) => {
                SlotReplacement::Replaced(Self::replace_argument(projector, position, replacement))
            }
            None => SlotReplacement::Absent,
        }
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
    fn linear_sum_partners_preserve_factorization_and_substitute_each_term() {
        initialize_vectors();
        for (input, expected) in [
            (
                "g(mink(D,a),mink(D,b))*(T(mink(D,b))+U(mink(D,b)))",
                "T(mink(D,a))+U(mink(D,a))",
            ),
            (
                "(x+y)^6*g(mink(4,a),mink(4,b))*(z*T(mink(4,b))+(x+y)*(U(mink(4,b))+V(mink(4,b))))",
                "(x+y)^6*(z*T(mink(4,a))+(x+y)*(U(mink(4,a))+V(mink(4,a))))",
            ),
            (
                "g(lor(D,a),dind(lor(D,b)))*(T(lor(D,b))+U(lor(D,b)))",
                "T(lor(D,a))+U(lor(D,a))",
            ),
            (
                "g(dind(lor(D,a)),lor(D,b))*(T(dind(lor(D,b)))+U(dind(lor(D,b))))",
                "T(dind(lor(D,a)))+U(dind(lor(D,a)))",
            ),
            (
                "slot_contraction_p(mink(4,a))*(T(mink(4,a))+U(mink(4,a)))",
                "T(slot_contraction_p(mink(4)))+U(slot_contraction_p(mink(4)))",
            ),
            (
                "g(mink(4,a),slot_contraction_p(label,mink(4)))*(T(mink(4,a))+U(mink(4,a)))",
                "T(slot_contraction_p(label,mink(4)))+U(slot_contraction_p(label,mink(4)))",
            ),
            (
                "g(mink(4,a),mink(4,b))*(chain(bis(4,i),bis(4,j),T(in,out,mink(4,b)))+chain(bis(4,i),bis(4,j),U(in,out,mink(4,b))))",
                "chain(bis(4,i),bis(4,j),T(in,out,mink(4,a)))+chain(bis(4,i),bis(4,j),U(in,out,mink(4,a)))",
            ),
            (
                "g(mink(4,a),mink(4,b))*(trace(bis(4),cyclic(T(in,out,mink(4,b))))+trace(bis(4),cyclic(U(in,out,mink(4,b)))))",
                "trace(bis(4),cyclic(T(in,out,mink(4,a))))+trace(bis(4),cyclic(U(in,out,mink(4,a))))",
            ),
            (
                "g(mink(4,4294967297),mink(4,f(x)))*(T(mink(4,f(x)))+U(mink(4,f(x))))",
                "T(mink(4,4294967297))+U(mink(4,4294967297))",
            ),
        ] {
            let expression = parse(input);
            let result = SlotContraction::run(expression.as_view(), true, true);
            assert_eq!(result, parse(expected), "{input}");
            assert_eq!(SlotContraction::run(result.as_view(), true, true), result);
        }
    }

    #[test]
    fn linear_metric_sums_contract_external_tensors_without_expanding_spectators() {
        initialize_vectors();
        let determinant = parse(
            "(g(mink(4,a),mink(4,b))*g(mink(4,c),mink(4,d))-g(mink(4,a),mink(4,d))*g(mink(4,c),mink(4,b)))*epsilon(mink(4,a),mink(4,e),mink(4,f),mink(4,h))",
        );
        let expected = parse(
            "g(mink(4,c),mink(4,d))*epsilon(mink(4,b),mink(4,e),mink(4,f),mink(4,h))-g(mink(4,c),mink(4,b))*epsilon(mink(4,d),mink(4,e),mink(4,f),mink(4,h))",
        );
        let result = SlotContraction::run(determinant.as_view(), false, false);
        // Epsilon normalization can move a common minus sign outside the
        // two-term tensor polynomial. There are no scalar spectators here.
        assert_eq!(result.expand(), expected.expand());
        assert_eq!(SlotContraction::run(result.as_view(), false, false), result);
        for (input, expected, rank_one) in [
            (
                "(x+y)^6*(z*g(mink(D,a),mink(D,b))+(u+v)*(g(mink(D,a),mink(D,c))+g(mink(D,a),mink(D,d))))*T(mink(D,a))",
                "(x+y)^6*(z*T(mink(D,b))+(u+v)*(T(mink(D,c))+T(mink(D,d))))",
                false,
            ),
            (
                "(g(lor(D,a),dind(lor(D,b)))+g(lor(D,c),dind(lor(D,b))))*T(lor(D,b))",
                "T(lor(D,a))+T(lor(D,c))",
                false,
            ),
            (
                "(g(mink(4,a),slot_contraction_p(mink(4)))+g(mink(4,a),slot_contraction_q(mink(4))))*T(mink(4,a))",
                "T(slot_contraction_p(mink(4)))+T(slot_contraction_q(mink(4)))",
                true,
            ),
        ] {
            let expression = parse(input);
            if rank_one {
                assert_eq!(
                    SlotContraction::run(expression.as_view(), false, false),
                    expression
                );
            }
            let result = SlotContraction::run(expression.as_view(), false, rank_one);
            assert_eq!(result, parse(expected), "{input}");
            assert_eq!(
                SlotContraction::run(result.as_view(), false, rank_one),
                result
            );
        }
        for input in [
            "(g(mink(4,a),mink(4,b))+g(mink(4,c),mink(4,d)))*T(mink(4,a))",
            "(g(mink(4,a),mink(4,b))+g(mink(5,a),mink(5,b)))*T(mink(4,a))",
            "(g(mink(4,a),mink(4,b))^3+g(mink(4,a),mink(4,c)))*T(mink(4,a))",
            "(g(mink(4,a),mink(4,b))+g(mink(4,a),mink(4,c)))*f(T(mink(4,a)))",
        ] {
            let expression = parse(input);
            assert_eq!(
                SlotContraction::run(expression.as_view(), false, false),
                expression
            );
        }
    }

    #[test]
    fn sum_substitution_rejects_missing_repeated_and_incompatible_slots() {
        initialize_vectors();
        for input in [
            "g(mink(4,a),mink(4,b))*(T(mink(4,b))+U(mink(4,c)))",
            "g(mink(4,a),mink(4,b))*(T(mink(4,b))+U(mink(5,b)))",
            "g(mink(4,a),mink(4,b))*(T(mink(4,b))+U(euc(4,b)))",
            "g(mink(4,a),mink(4,b))*(T(mink(4,b),mink(4,b))+U(mink(4,b)))",
            "g(mink(4,a),mink(4,b))*(T(mink(4,b))*V(mink(4,b))+U(mink(4,b)))",
            "g(lor(D,a),dind(lor(D,b)))*(T(lor(D,b))+U(dind(lor(D,b))))",
            "g(lor(D,a),dind(lor(D,b)))*(T(lor(D,b))*V(dind(lor(D,b)))+U(lor(D,b)))",
            "g(mink(4,a),mink(4,b))*(T(f(mink(4,b)))+U(f(mink(4,b))))",
            "g(mink(4,a),mink(4,b))*(T(mink(4,b))^3+U(mink(4,b)))",
            "g(mink(4,a),mink(4,b))*(T(mink(4,b))^3*U(mink(4,b))+V(mink(4,b)))",
            "g(mink(4,a),mink(4,b))*(slot_contraction_scalar(mink(4,b))+U(mink(4,b)))",
            "g(mink(4,a),mink(4,b))*(trace(bis(4),cyclic(T(in,out,mink(4,b)),U(in,out,mink(4,b))))+V(mink(4,b)))",
        ] {
            let expression = parse(input);
            assert_eq!(
                SlotContraction::run(expression.as_view(), true, true),
                expression,
                "{input}"
            );
        }
    }

    #[test]
    fn metric_components_transport_boundaries_and_close_loops() {
        crate::representations::initialize();
        for (input, expected) in [
            (
                "g(mink(D,a),mink(D,b))*g(mink(D,b),mink(D,c))*g(mink(D,c),mink(D,d))*g(mink(D,d),mink(D,e))*T(mink(D,e))",
                "T(mink(D,a))",
            ),
            (
                "g(mink(D,a),mink(D,b))*g(mink(D,b),mink(D,c))*g(mink(D,c),mink(D,d))*g(mink(D,d),mink(D,a))",
                "D",
            ),
            (
                "g(lor(D,a),dind(lor(D,b)))*g(lor(D,b),dind(lor(D,c)))*g(lor(D,c),dind(lor(D,d)))*g(lor(D,d),dind(lor(D,e)))*T(lor(D,e))",
                "T(lor(D,a))",
            ),
            (
                "g(dind(lor(D,a)),lor(D,b))*g(dind(lor(D,b)),lor(D,c))*g(dind(lor(D,c)),lor(D,d))*g(dind(lor(D,d)),lor(D,a))",
                "D",
            ),
            (
                "g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))*g(mink(5,a),mink(5,b))*g(mink(5,b),mink(5,c))",
                "g(mink(4,a),mink(4,c))*g(mink(5,a),mink(5,c))",
            ),
        ] {
            let result = SlotContraction::run(parse(input).as_view(), false, false);
            assert_eq!(result, parse(expected), "{input}");
            assert_eq!(SlotContraction::run(result.as_view(), false, false), result);
        }
    }

    #[test]
    fn ambiguous_metric_incidence_keeps_ordered_substitution() {
        crate::representations::initialize();
        for source in [
            "g(mink(4,a),mink(4,b))*g(mink(4,a),mink(4,c))*g(mink(4,a),mink(4,d))*g(mink(4,x),mink(4,y))",
            "g(mink(4,a),mink(4,b))*g(mink(4,b),mink(4,c))*g(mink(4,c),mink(4,d))*T(mink(4,b))",
        ] {
            let expression = parse(source);
            let AtomView::Mul(product) = expression.as_view() else {
                panic!("expected product");
            };
            assert!(
                SlotContraction::contract_metric_components(product, &mut SlotMatcher::default())
                    .is_none()
            );
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
